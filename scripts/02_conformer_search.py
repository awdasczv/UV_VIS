"""
02_conformer_search.py
----------------------
각 토토머의 컨포머 탐색 (요구사항 7, 8번).

두 가지 백엔드를 지원한다.
  * backend="crest"  : CREST(GFN2-xTB, metadynamics) 사용. Linux/WSL 에서만 가능.
  * backend="etkdg"  : RDKit ETKDGv3 로 다수 구조 생성 -> xTB(GFN2) 최적화
                       -> 에너지 + heavy-atom RMSD 로 중복 제거.
                       Windows 네이티브에서 CREST 없이 쓰는 대체 경로.

두 경우 모두 결과 형식은 같다:
  conformers/<tautomer>/conformers.xyz        (중복 제거된 앙상블)
  conformers/<tautomer>/ensemble.json         (에너지, 상대에너지, 볼츠만 가중치)
  conformers/<tautomer>/selected.json         (누적 가중치 기준 대표 컨포머)

실행 예:
  python scripts/02_conformer_search.py --backend etkdg --nconf 300
  python scripts/02_conformer_search.py --backend crest --crest-exe crest
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import (CONFORMERS, LOGS, MOL_CONFIG, ROOT, STRUCTURES, Geometry,
                       boltzmann_weights, read_multi_xyz, read_xyz,
                       rel_energies_kcal, save_checkpoint,
                       select_by_cumulative_weight)

TAUTOMERS = ["enolA", "enolB", "diketo"]


# --------------------------------------------------------------- xtb 호출
def find_xtb() -> str:
    exe = shutil.which("xtb")
    if exe:
        return exe
    # micromamba 환경 안
    cand = Path(sys.prefix) / "Library" / "bin" / "xtb.exe"
    if cand.exists():
        return str(cand)
    raise FileNotFoundError("xtb 실행파일을 찾을 수 없습니다.")


def xtb_optimize(geom: Geometry, workdir: Path, solvent: str | None = "ethanol",
                 gfn: str = "2", opt_level: str = "normal") -> tuple[Geometry, float]:
    """GFN2-xTB 구조 최적화. (최적화된 구조, 총에너지[hartree]) 반환."""
    workdir.mkdir(parents=True, exist_ok=True)
    inp = workdir / "in.xyz"
    geom.write(inp, "xtb input")

    cmd = [find_xtb(), "in.xyz", "--gfn", gfn, "--opt", opt_level, "--chrg", "0", "--uhf", "0"]
    if solvent:
        cmd += ["--alpb", solvent]

    env = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
               XTBPATH=str(workdir))
    # 주의: 한국어 Windows 는 기본 로케일이 cp949 라, text=True 만 주면 파이썬이
    # xtb 출력의 UTF-8 박스문자(e.g. U+2500)를 cp949 로 디코딩하려다
    # UnicodeDecodeError 로 리더 스레드가 죽고 stdout 이 None 이 된다.
    # 반드시 encoding/errors 를 명시해야 한다.
    res = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True,
                         encoding="utf-8", errors="replace", env=env)
    (workdir / "xtb.out").write_text((res.stdout or "") + "\n=== STDERR ===\n" +
                                     (res.stderr or ""),
                                     encoding="utf-8", errors="replace")

    optxyz = workdir / "xtbopt.xyz"
    if not optxyz.exists():
        raise RuntimeError(f"xtb 최적화 실패 (returncode={res.returncode}); "
                           f"로그: {workdir/'xtb.out'}")

    opt = read_xyz(optxyz)
    # xtbopt.xyz 의 comment 줄에 " energy: <E> ..." 형태로 총에너지가 들어있다
    energy = None
    for tok in opt.comment.replace(":", " ").split():
        try:
            v = float(tok)
            if -2000 < v < 0:      # 총에너지로 그럴듯한 값
                energy = v
                break
        except ValueError:
            continue
    if energy is None:             # 백업: stdout 파싱
        for ln in (res.stdout or "").splitlines():
            if "TOTAL ENERGY" in ln:
                energy = float(ln.split()[3])
    if energy is None:
        raise RuntimeError("xtb 총에너지를 파싱하지 못했습니다.")
    return opt, energy


# ------------------------------------------------------------- 중복 제거
def kabsch_rmsd(a: np.ndarray, b: np.ndarray) -> float:
    """무게중심 정렬 + Kabsch 회전 후 RMSD (원자 순서가 같다고 가정)."""
    a = a - a.mean(axis=0)
    b = b - b.mean(axis=0)
    h = a.T @ b
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    a_rot = (r @ a.T).T
    return float(np.sqrt(((a_rot - b) ** 2).sum() / len(a)))


def _symmetry_aware_rmsd_factory(template_mol):
    """
    RDKit 의 대칭 인식 RMSD (GetBestRMS) 를 쓰는 비교 함수를 만든다.

    왜 필요한가:
      단순 Kabsch RMSD 는 tert-butyl 의 세 메틸기 교환, 페닐 고리 180도 뒤집기,
      그리고 거울상 배치를 '서로 다른 구조'로 본다. 실제로는 같은 컨포머다.
      GetBestRMS 는 분자 그래프의 자기동형사상(automorphism)을 모두 시도해
      최소 RMSD 를 찾으므로 이런 가짜 중복을 제대로 걸러낸다.

    원자 순서는 RDKit -> xyz -> xtb -> xtbopt.xyz 전 과정에서 보존되므로
    좌표만 갈아끼우면 된다. 수소는 빼고 heavy atom 만 비교한다.
    """
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    from rdkit.Geometry import Point3D

    if template_mol is None:
        return None

    heavy_mol = Chem.RemoveHs(Chem.Mol(template_mol))
    heavy_idx = [a.GetIdx() for a in template_mol.GetAtoms()
                 if a.GetSymbol() != "H"]
    if heavy_mol.GetNumAtoms() != len(heavy_idx):
        return None                      # 매핑이 안 맞으면 포기하고 Kabsch 로

    def make(coords: np.ndarray):
        m = Chem.Mol(heavy_mol)
        m.RemoveAllConformers()
        conf = Chem.Conformer(m.GetNumAtoms())
        for k, ai in enumerate(heavy_idx):
            x, y, z = coords[ai]
            conf.SetAtomPosition(k, Point3D(float(x), float(y), float(z)))
        m.AddConformer(conf, assignId=True)
        return m

    def rmsd(a: np.ndarray, b: np.ndarray) -> float:
        # GetBestRMS 는 probe 를 정렬시키므로 사본을 넘긴다
        return float(rdMolAlign.GetBestRMS(make(a), make(b)))

    return rmsd


def deduplicate(geoms: list[Geometry], energies: list[float],
                e_thr_kcal: float = 0.10, rmsd_thr: float = 0.125,
                e_window_kcal: float = 6.0, template_mol=None):
    """
    CREST 기본값과 비슷한 기준으로 중복 제거.
      - 에너지 창 밖(기준 대비 e_window_kcal 이상) 구조는 버림
      - 에너지 차 < e_thr_kcal 이고 RMSD < rmsd_thr 이면 같은 컨포머
      - RMSD 는 가능하면 대칭 인식(GetBestRMS), 아니면 heavy-atom Kabsch
    """
    order = np.argsort(energies)
    geoms = [geoms[i] for i in order]
    energies = [energies[i] for i in order]

    rel = rel_energies_kcal(energies)
    keep_idx = [i for i in range(len(geoms)) if rel[i] <= e_window_kcal]

    sym_rmsd = _symmetry_aware_rmsd_factory(template_mol)
    heavy = [i for i, s in enumerate(geoms[0].symbols) if s != "H"]
    if sym_rmsd is None:
        print("    (경고) 대칭 인식 RMSD 사용 불가 -> heavy-atom Kabsch RMSD 로 대체")

    uniq_g, uniq_e = [], []
    for i in keep_idx:
        g, e = geoms[i], energies[i]
        dup = False
        for gu, eu in zip(uniq_g, uniq_e):
            de = abs(e - eu) * 627.5094740631
            if de >= e_thr_kcal:
                continue
            r = (sym_rmsd(g.coords, gu.coords) if sym_rmsd
                 else kabsch_rmsd(g.coords[heavy], gu.coords[heavy]))
            if r < rmsd_thr:
                dup = True
                break
        if not dup:
            uniq_g.append(g)
            uniq_e.append(e)
    return uniq_g, uniq_e


# ----------------------------------------------------------- 백엔드: ETKDG
def backend_etkdg(taut: str, nconf: int, solvent: str | None, outdir: Path):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem
    RDLogger.DisableLog("rdApp.warning")

    spec = json.loads(MOL_CONFIG.read_text(encoding="utf-8"))
    smi = next(t["smiles"] for t in spec["tautomers"] if t["id"] == taut)

    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    params = AllChem.ETKDGv3()
    params.randomSeed = 20260720
    # pruneRmsThresh 는 쓰지 않는다.
    # (0.5 를 줬더니 200개 요청이 6~7개로 잘려나갔다. RDKit 의 이 가지치기는
    #  수소를 포함한 전체 원자 RMSD 를 쓰고 대칭을 고려하지 않아 지나치게 공격적이다.
    #  중복 제거는 xTB 최적화 후에 에너지 + 대칭 인식 RMSD 로 제대로 한다.)
    params.pruneRmsThresh = -1.0
    params.numThreads = 0
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=nconf, params=params))
    AllChem.MMFFOptimizeMoleculeConfs(mol, mmffVariant="MMFF94s",
                                      maxIters=2000, numThreads=0)
    print(f"  ETKDG+MMFF: {len(cids)} 개 초기 구조")

    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    raw = [Geometry(syms, np.array(mol.GetConformer(c).GetPositions()), f"etkdg_{c}")
           for c in cids]

    # xTB 최적화 (체크포인트: 하나 끝날 때마다 저장)
    ck = outdir / "xtb_opt_checkpoint.json"
    done = json.loads(ck.read_text(encoding="utf-8")) if ck.exists() else {}
    geoms, energies = [], []
    t0 = time.time()
    for k, g in enumerate(raw):
        key = str(k)
        wd = outdir / "xtb" / f"c{k:04d}"
        if key in done and (wd / "xtbopt.xyz").exists():
            geoms.append(read_xyz(wd / "xtbopt.xyz"))
            energies.append(done[key])
            continue
        try:
            og, e = xtb_optimize(g, wd, solvent=solvent)
        except RuntimeError as err:
            print(f"    [경고] 구조 {k} xtb 실패: {err}")
            continue
        geoms.append(og)
        energies.append(e)
        done[key] = e
        ck.write_text(json.dumps(done), encoding="utf-8")
        if (k + 1) % 20 == 0:
            print(f"    xTB 최적화 {k+1}/{len(raw)}  ({time.time()-t0:.0f}s)")
    print(f"  xTB 최적화 완료: {len(geoms)} 개 성공, {time.time()-t0:.0f}s")
    return geoms, energies, mol


# ----------------------------------------------------------- 백엔드: CREST
def backend_crest(taut: str, crest_exe: str, solvent: str | None,
                  outdir: Path, nthreads: int):
    src = STRUCTURES / f"{taut}.xyz"
    wd = outdir / "crest_run"
    wd.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, wd / "in.xyz")

    cmd = [crest_exe, "in.xyz", "--gfn2", "-T", str(nthreads)]
    if solvent:
        cmd += ["--alpb", solvent]
    print("  실행:", " ".join(cmd))
    with (LOGS / f"crest_{taut}.log").open("w", encoding="utf-8") as log:
        subprocess.run(cmd, cwd=wd, stdout=log, stderr=subprocess.STDOUT, text=True)

    ens = wd / "crest_conformers.xyz"
    if not ens.exists():
        raise RuntimeError(f"CREST 결과가 없습니다: {ens}")
    geoms = read_multi_xyz(ens)
    energies = [float(g.comment.split()[0]) for g in geoms]
    print(f"  CREST 앙상블: {len(geoms)} 개")
    return geoms, energies, None


# ------------------------------------------------------------------ main
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--backend", choices=["etkdg", "crest"], default="etkdg")
    ap.add_argument("--nconf", type=int, default=200, help="ETKDG 초기 구조 개수")
    ap.add_argument("--solvent", default="ethanol", help="xtb ALPB 용매 ('none' 이면 기체상)")
    ap.add_argument("--crest-exe", default="crest")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--target-weight", type=float, default=0.90)
    ap.add_argument("--max-select", type=int, default=3)
    ap.add_argument("--tautomers", nargs="*", default=TAUTOMERS)
    args = ap.parse_args()

    solvent = None if args.solvent.lower() == "none" else args.solvent
    summary = {}

    for taut in args.tautomers:
        print(f"\n=== 컨포머 탐색: {taut} (backend={args.backend}, solvent={solvent}) ===")
        outdir = CONFORMERS / taut
        outdir.mkdir(parents=True, exist_ok=True)

        if args.backend == "crest":
            geoms, energies, tmpl = backend_crest(taut, args.crest_exe, solvent,
                                                  outdir, args.threads)
        else:
            geoms, energies, tmpl = backend_etkdg(taut, args.nconf, solvent, outdir)

        uniq_g, uniq_e = deduplicate(geoms, energies, template_mol=tmpl)
        w = boltzmann_weights(uniq_e)
        rel = rel_energies_kcal(uniq_e)
        print(f"  중복 제거 후: {len(uniq_g)} 개 고유 컨포머")

        # 앙상블 xyz 저장
        with (outdir / "conformers.xyz").open("w", encoding="utf-8") as fh:
            for i, (g, e) in enumerate(zip(uniq_g, uniq_e)):
                fh.write(f"{g.natoms}\n{taut}_c{i:03d} E={e:.10f} Eh "
                         f"dE={rel[i]:.4f} kcal/mol w={w[i]:.5f}\n")
                fh.write(g.to_xyz_block() + "\n")
        for i, g in enumerate(uniq_g):
            g.write(outdir / f"{taut}_c{i:03d}.xyz",
                    f"{taut} conformer {i} E={uniq_e[i]:.10f} Eh")

        ensemble = {
            "tautomer": taut, "backend": args.backend, "solvent_xtb_alpb": solvent,
            "level": "GFN2-xTB",
            "n_raw": len(geoms), "n_unique": len(uniq_g),
            "dedup_criteria": {"dE_kcal": 0.10, "rmsd_A": 0.125,
                               "rmsd_type": "RDKit GetBestRMS (대칭 인식, heavy atom)",
                               "energy_window_kcal": 6.0},
            "temperature_K": 298.15,
            "conformers": [
                {"conf_id": f"{taut}_c{i:03d}",
                 "energy_hartree": float(uniq_e[i]),
                 "rel_energy_kcalmol": float(rel[i]),
                 "boltzmann_weight": float(w[i])}
                for i in range(len(uniq_g))
            ],
        }
        save_checkpoint(outdir / "ensemble.json", ensemble)

        picked, cum = select_by_cumulative_weight(w, args.target_weight, args.max_select)
        sel = {
            "tautomer": taut,
            "target_cumulative_weight": args.target_weight,
            "max_conformers": args.max_select,
            "achieved_cumulative_weight": float(cum),
            "selected": [
                {"conf_id": f"{taut}_c{i:03d}",
                 "rel_energy_kcalmol": float(rel[i]),
                 "boltzmann_weight": float(w[i]),
                 "xyz": f"conformers/{taut}/{taut}_c{i:03d}.xyz"}
                for i in picked
            ],
        }
        save_checkpoint(outdir / "selected.json", sel)
        print(f"  대표 컨포머 {len(picked)} 개 선택, 누적 가중치 = {cum:.3f}")
        for s in sel["selected"]:
            print(f"    {s['conf_id']}  dE={s['rel_energy_kcalmol']:.2f} kcal/mol  "
                  f"w={s['boltzmann_weight']:.3f}")
        summary[taut] = {"n_unique": len(uniq_g), "cum_weight": cum,
                         "selected": [s["conf_id"] for s in sel["selected"]]}

    save_checkpoint(CONFORMERS / "search_summary.json", summary)
    print("\n완료. conformers/search_summary.json 저장")
    return 0


if __name__ == "__main__":
    sys.exit(main())
