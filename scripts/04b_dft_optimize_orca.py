"""
04b_dft_optimize_orca.py
------------------------
ORCA 로 DFT 구조 최적화 (요구사항 9번).

왜 이 단계가 필요한지가 계산으로 드러났다
--------------------------------------------
기저셋 비교 실험(enolA_c000, B3LYP, 에탄올 CPCM, GFN2-xTB 구조)에서:

    def2-SVP   (432 기저함수, 확산 없음)  334.0 nm
    6-31+G(d)  (481, 확산 있음)           340.6 nm
    def2-SVPD  (645, 확산 있음)           341.3 nm
    def2-TZVP  (845, 확산 없음)           340.0 nm

확산함수 유무와 무관하게 큰 기저셋은 모두 340~341 nm 로 수렴한다.
즉 기저셋 수렴값은 약 340.5 nm 다.

그런데 문헌(ACS Omega 2026)은 **같은 B3LYP/6-31+G(d)/PCM** 으로 360.8 nm 를
보고했다. 함수도 기저셋도 용매도 같은데 20 nm 가 차이난다.
남은 변수는 **구조** 하나뿐이다. 문헌은 DFT 최적화 구조를 썼고 우리는
GFN2-xTB 구조를 썼다.

물리적으로도 말이 된다. xTB 구조의 아릴 비틀림각은 27.1°/22.2° 인데,
평면에 가까울수록 pi 공액이 좋아져 흡수가 장파장으로 간다.
이 스크립트로 DFT 최적화를 해서 그 가설을 검증한다.

실행:
  .\scripts\run.ps1 scripts\04b_dft_optimize_orca.py
  .\scripts\run.ps1 scripts\04b_dft_optimize_orca.py --tautomers enolA --max-conformers 1
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orca_common import (build_input, classify_orca_failure, find_orca,
                         run_orca)
from qc_common import (CALCULATIONS, CONFORMERS, INPUTS, LOGS, Geometry,
                       load_checkpoint, read_xyz, save_checkpoint)

OUTROOT = CALCULATIONS / "01_dft_opt"
SOLVENT_NAMES = {"ethanol": "Ethanol", "methanol": "Methanol", "water": "Water"}


def read_orca_opt_xyz(workdir: Path, name: str = "opt") -> Geometry | None:
    """ORCA 가 남기는 최적화 결과 xyz (<name>.xyz) 를 읽는다."""
    p = workdir / f"{name}.xyz"
    return read_xyz(p) if p.exists() else None


def aryl_twists(geom: Geometry) -> dict | None:
    """
    최적화 전후로 발색단 평면성이 얼마나 변했는지 보기 위해
    아릴 고리 비틀림각을 다시 잰다. (02b 와 같은 정의)
    """
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import rdDetermineBonds, rdMolTransforms
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        return None

    tmp = Path(str(Path.cwd() / "_twist_tmp.xyz"))
    geom.write(tmp)
    try:
        raw = Chem.MolFromXYZFile(str(tmp))
        if raw is None:
            return None
        mol = Chem.Mol(raw)
        rdDetermineBonds.DetermineBonds(mol, charge=0)
    except Exception:                                    # noqa: BLE001
        return None
    finally:
        tmp.unlink(missing_ok=True)

    conf = mol.GetConformer()

    def tors(a, b, c, d):
        v = abs(rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)) % 180.0
        return round(min(v, 180.0 - v), 1)

    out = {}
    m = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[OX1]=[CX3]([c])[CX3]=[CX3]([c])[OX2H]"))
    if m:
        o_ket, c_ket, ar_ket, c_a, c_enol, ar_enol, o_enol = m[0]
        n1 = [x.GetIdx() for x in mol.GetAtomWithIdx(ar_ket).GetNeighbors()
              if x.GetIdx() != c_ket and x.GetSymbol() == "C"]
        n2 = [x.GetIdx() for x in mol.GetAtomWithIdx(ar_enol).GetNeighbors()
              if x.GetIdx() != c_enol and x.GetSymbol() == "C"]
        if n1:
            out["aryl_ketone_twist_deg"] = tors(o_ket, c_ket, ar_ket, n1[0])
        if n2:
            out["aryl_enol_twist_deg"] = tors(o_enol, c_enol, ar_enol, n2[0])
        p = lambda i: np.array(conf.GetAtomPosition(i))
        hs = [x.GetIdx() for x in mol.GetAtomWithIdx(o_enol).GetNeighbors()
              if x.GetSymbol() == "H"]
        if hs:
            out["chelate_H_to_O_A"] = round(
                float(np.linalg.norm(p(hs[0]) - p(o_ket))), 3)
    return out or None


def optimize_one(taut: str, conf_id: str, cfg: dict) -> dict:
    opt_cfg = cfg["orca"]["optimization"]
    orca_cfg = cfg["orca"]
    outdir = OUTROOT / taut / conf_id
    ck = outdir / "result.json"

    done = load_checkpoint(ck)
    if done and done.get("ok"):
        print(f"  [건너뜀] {taut}/{conf_id} 이미 완료 "
              f"(E = {done['energy_hartree']:.8f} Eh)")
        return done

    src = CONFORMERS / taut / f"{conf_id}.xyz"
    if not src.exists():
        print(f"  [건너뜀] {taut}/{conf_id}: 구조 없음 ({src})")
        return {"ok": False, "reason": "geometry missing"}

    g0 = read_xyz(src)
    solv = SOLVENT_NAMES.get(opt_cfg["solvent"]) if opt_cfg.get("pcm") else None
    print(f"  [최적화] {taut}/{conf_id}  "
          f"{opt_cfg['functional']}/{opt_cfg['basis']}  용매={solv or '없음'}")

    inp = build_input(
        g0.symbols, g0.coords,
        functional=opt_cfg["functional"], basis=opt_cfg["basis"],
        nstates=0, solvent=solv, optimize=True,
        nprocs=orca_cfg["nprocs"], maxcore_mb=orca_cfg["maxcore_mb"],
        rijcosx=orca_cfg["rijcosx"], aux_basis=orca_cfg["aux_basis"],
        comment=f"geometry optimization: {taut} / {conf_id}",
    )
    out_path, seconds = run_orca(inp, outdir, name="opt")
    text = out_path.read_text(encoding="utf-8", errors="replace")

    rec: dict = {
        "engine": "ORCA 6.1.1", "tautomer": taut, "conf_id": conf_id,
        "level": f"{opt_cfg['functional']}/{opt_cfg['basis']}",
        "pcm": bool(solv), "solvent": opt_cfg.get("solvent") if solv else None,
        "source_xyz": str(src), "wall_seconds": round(seconds, 1),
        "twists_before": aryl_twists(g0),
    }

    converged = "HURRAY" in text or "THE OPTIMIZATION HAS CONVERGED" in text
    geom = read_orca_opt_xyz(outdir, "opt")
    if not converged or geom is None:
        rec["ok"] = False
        cls = classify_orca_failure(out_path)
        rec["failure"] = cls
        data = load_checkpoint(LOGS / "failures.json") or {"failures": []}
        data["failures"].append({"tag": f"opt-orca:{taut}/{conf_id}",
                                 "engine": "ORCA", **cls})
        save_checkpoint(LOGS / "failures.json", data)
        print(f"    [실패] 수렴={converged}  분류={cls['code']}")
        save_checkpoint(ck, rec)
        return rec

    m = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", text)
    energy = float(m[-1]) if m else None
    nsteps = len(re.findall(r"GEOMETRY OPTIMIZATION CYCLE", text))

    geom.write(outdir / "optimized.xyz",
               f"{taut}/{conf_id} {rec['level']} E={energy:.10f} Eh")
    rec.update({
        "ok": True, "energy_hartree": energy, "n_opt_cycles": nsteps,
        "optimized_xyz": str(outdir / "optimized.xyz"),
        "twists_after": aryl_twists(geom),
    })
    print(f"    완료 {seconds:.0f}초, {nsteps} 사이클, E = {energy:.8f} Eh")
    b, a = rec["twists_before"], rec["twists_after"]
    if b and a:
        print("    발색단 기하 변화 (xTB -> DFT):")
        for k in sorted(set(b) | set(a)):
            if k in b and k in a:
                unit = "A" if k.endswith("_A") else "deg"
                print(f"      {k:26s} {b[k]:7.2f} -> {a[k]:7.2f} {unit} "
                      f"({a[k]-b[k]:+.2f})")
    save_checkpoint(ck, rec)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tautomers", nargs="*", default=["enolA", "enolB", "diketo"])
    ap.add_argument("--max-conformers", type=int, default=1)
    args = ap.parse_args()

    cfg = json.loads((INPUTS / "calc_config.json").read_text(encoding="utf-8"))
    if "optimization" not in cfg["orca"]:
        print("inputs/calc_config.json 의 orca 섹션에 optimization 설정이 없습니다.")
        return 1
    o = cfg["orca"]["optimization"]
    print(f"엔진: ORCA ({find_orca()})")
    print(f"수준: {o['functional']}/{o['basis']}  PCM={o.get('pcm')} "
          f"({o.get('solvent')})  병렬 {cfg['orca']['nprocs']}")

    summary = {}
    for taut in args.tautomers:
        sel = CONFORMERS / taut / "selected.json"
        if not sel.exists():
            print(f"[건너뜀] {taut}: selected.json 없음")
            continue
        confs = json.loads(sel.read_text(encoding="utf-8"))["selected"]
        confs = confs[:args.max_conformers]
        print(f"\n=== {taut}: {len(confs)} 개 구조 ===")
        summary[taut] = [optimize_one(taut, c["conf_id"], cfg) for c in confs]
        save_checkpoint(OUTROOT / "summary.json", summary)

    print("\n=== 요약 ===")
    for taut, recs in summary.items():
        for r in recs:
            if r.get("ok"):
                print(f"  {taut}/{r['conf_id']:14s} E={r['energy_hartree']:.8f} Eh "
                      f"({r['n_opt_cycles']} 사이클, {r['wall_seconds']:.0f}초)")
            else:
                print(f"  {taut}/{r.get('conf_id','?'):14s} 실패")
    return 0


if __name__ == "__main__":
    sys.exit(main())
