"""
04_dft_optimize.py
------------------
대표 컨포머의 DFT 구조 최적화 (요구사항 9번).

  - inputs/calc_config.json 의 active_profile 설정을 따른다.
  - conformers/<taut>/selected.json 에서 대표 컨포머를 읽는다.
  - scope="lowest_only" 면 토토머별 최저에너지 컨포머 1개만 DFT 재최적화하고
    나머지는 GFN2-xTB 구조를 그대로 쓴다 (비용 관리).
  - 구조 하나가 끝날 때마다 즉시 체크포인트를 저장하므로,
    중간에 죽어도 다시 실행하면 남은 것만 이어서 한다 (요구사항 17번).
  - 실패하면 원인을 분류해 logs/failures.json 에 기록한다 (요구사항 16번).

실행:  .\scripts\run.ps1 scripts\04_dft_optimize.py
       .\scripts\run.ps1 scripts\04_dft_optimize.py --tautomers enolA
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import psi4
from psi4_helpers import (classify_failure, init_psi4, make_molecule, pcm_block)
from qc_common import (CALCULATIONS, CONFORMERS, INPUTS, LOGS, Geometry,
                       load_checkpoint, read_xyz, save_checkpoint)

OUTROOT = CALCULATIONS / "01_dft_opt"


def load_profile() -> dict:
    cfg = json.loads((INPUTS / "calc_config.json").read_text(encoding="utf-8"))
    prof = cfg["profiles"][cfg["active_profile"]]
    prof["_profile_name"] = cfg["active_profile"]
    return prof


def record_failure(tag: str, err_text: str, extra: dict | None = None) -> dict:
    """실패를 분류해 logs/failures.json 에 누적 기록."""
    path = LOGS / "failures.json"
    data = load_checkpoint(path) or {"failures": []}
    cls = classify_failure(err_text)
    entry = {"tag": tag, "code": cls["code"], "remedy": cls["remedy"],
             "error_tail": err_text[-1200:]}
    if extra:
        entry.update(extra)
    data["failures"].append(entry)
    save_checkpoint(path, data)
    print(f"    [실패 분류] {cls['code']}")
    print(f"    [대응 방안] {cls['remedy']}")
    return entry


def optimize_one(taut: str, conf_id: str, xyz_path: Path, prof: dict) -> dict:
    """구조 하나를 DFT 로 최적화."""
    opt = prof["geometry_optimization"]
    outdir = OUTROOT / taut / conf_id
    outdir.mkdir(parents=True, exist_ok=True)
    ck = outdir / "result.json"

    done = load_checkpoint(ck)
    if done and done.get("ok"):
        print(f"  [건너뜀] {taut}/{conf_id} 이미 완료 "
              f"(E = {done['energy_hartree']:.8f} Eh)")
        return done

    print(f"  [최적화] {taut}/{conf_id}  "
          f"{opt['functional']}/{opt['basis']} PCM={opt['pcm']}")
    rec: dict = {"tautomer": taut, "conf_id": conf_id,
                 "level": f"{opt['functional']}/{opt['basis']}",
                 "pcm": opt["pcm"], "solvent": opt["solvent"] if opt["pcm"] else None,
                 "source_xyz": str(xyz_path)}
    t0 = time.time()
    try:
        psi4.core.clean()
        psi4.core.clean_options()
        init_psi4(LOGS / f"psi4_opt_{taut}_{conf_id}.out",
                  memory=prof["memory"], nthreads=prof["nthreads"])

        g = read_xyz(xyz_path)
        mol = make_molecule(g.symbols, g.coords)

        options = {
            "basis": opt["basis"],
            "scf_type": "df",
            "df_scf_guess": True,
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "maxiter": 250,
            "geom_maxiter": opt.get("geom_maxiter", 60),
            "g_convergence": "gau",
        }
        if opt["pcm"]:
            options.update({"pcm": True, "pcm_scf_type": "total"})
        psi4.set_options(options)
        if opt["pcm"]:
            # 구조 최적화(바닥상태)에는 평형 용매화를 쓴다.
            psi4.pcm_helper(pcm_block(opt["solvent"], nonequilibrium=False))

        energy = psi4.optimize(opt["functional"], molecule=mol)
        mol.update_geometry()

        # 최적화된 좌표 저장 (bohr -> angstrom)
        coords = mol.geometry().to_array() * 0.52917720859
        syms = [mol.symbol(i) for i in range(mol.natom())]
        geom = Geometry(syms, coords, f"{taut}/{conf_id} {rec['level']} E={energy:.10f}")
        geom.write(outdir / "optimized.xyz")

        rec.update({"ok": True, "energy_hartree": float(energy),
                    "seconds": round(time.time() - t0, 1),
                    "optimized_xyz": str((outdir / "optimized.xyz"))})
        print(f"    완료: E = {energy:.8f} Eh  ({rec['seconds']:.0f} s)")
    except Exception:                                    # noqa: BLE001
        err = traceback.format_exc()
        rec.update({"ok": False, "seconds": round(time.time() - t0, 1)})
        rec["failure"] = record_failure(f"opt:{taut}/{conf_id}", err,
                                        {"level": rec["level"]})
    save_checkpoint(ck, rec)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tautomers", nargs="*",
                    default=["enolA", "enolB", "diketo"])
    args = ap.parse_args()

    prof = load_profile()
    opt = prof["geometry_optimization"]
    print(f"프로파일: {prof['_profile_name']}  ({prof['description']})")
    if not opt.get("enabled", True):
        print("geometry_optimization.enabled = false -> DFT 최적화를 건너뜁니다.")
        return 0
    print(f"수준: {opt['functional']}/{opt['basis']}  PCM={opt['pcm']} "
          f"({opt['solvent']})  scope={opt['scope']}")

    summary = {}
    for taut in args.tautomers:
        sel_path = CONFORMERS / taut / "selected.json"
        if not sel_path.exists():
            print(f"[건너뜀] {taut}: {sel_path} 없음. 먼저 02_conformer_search.py 실행.")
            continue
        sel = json.loads(sel_path.read_text(encoding="utf-8"))
        targets = sel["selected"]
        if opt["scope"] == "lowest_only":
            targets = targets[:1]

        print(f"\n=== {taut}: {len(targets)} 개 구조 DFT 최적화 ===")
        recs = []
        for s in targets:
            xyz = Path(s["xyz"])
            if not xyz.is_absolute():
                xyz = CONFORMERS.parent / s["xyz"]
            recs.append(optimize_one(taut, s["conf_id"], xyz, prof))
        summary[taut] = recs
        save_checkpoint(OUTROOT / "summary.json", summary)

    print("\n=== DFT 최적화 요약 ===")
    for taut, recs in summary.items():
        for r in recs:
            status = f"E={r['energy_hartree']:.8f} Eh" if r.get("ok") else \
                     f"실패({r.get('failure', {}).get('code', '?')})"
            print(f"  {taut}/{r['conf_id']:16s} {status}  {r.get('seconds', 0):.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
