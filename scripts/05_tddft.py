"""
05_tddft.py
-----------
TD-DFT 들뜬상태 계산 (요구사항 10, 11번).

  - 계산 대상: 04_dft_optimize.py 의 최적화 구조가 있으면 그것, 없으면 GFN2-xTB 구조.
  - 각 (구조 x 이론수준 x 용매) 조합마다 계산하고 결과를 즉시 체크포인트로 저장한다.
    이미 끝난 조합은 자동으로 건너뛴다 -> 중간에 죽어도 처음부터 안 해도 된다.
  - 각 전이에 대해 에너지/파장/진동자세기/주요 오비탈 전이를 기록한다.
  - 실패 시 원인을 분류해 기록하고, 다음 조합으로 계속 진행한다.

실행:
  .\scripts\run.ps1 scripts\05_tddft.py
  .\scripts\run.ps1 scripts\05_tddft.py --tautomers enolA --levels b3lyp_631gd --solvents ethanol
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
from psi4_helpers import (classify_failure, init_psi4, make_molecule,
                          parse_tdscf_results, pcm_block)
from qc_common import (CALCULATIONS, CONFORMERS, INPUTS, LOGS, load_checkpoint,
                       read_xyz, save_checkpoint)

OPTROOT = CALCULATIONS / "01_dft_opt"
OUTROOT = CALCULATIONS / "02_tddft"


def load_profile() -> dict:
    cfg = json.loads((INPUTS / "calc_config.json").read_text(encoding="utf-8"))
    prof = cfg["profiles"][cfg["active_profile"]]
    prof["_profile_name"] = cfg["active_profile"]
    return prof


def record_failure(tag: str, err_text: str, extra: dict | None = None) -> dict:
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


def geometry_for(taut: str, conf_id: str) -> tuple[Path, str]:
    """DFT 최적화 구조가 있으면 그걸, 없으면 xTB 구조를 쓴다."""
    dft = OPTROOT / taut / conf_id / "optimized.xyz"
    if dft.exists():
        return dft, "DFT"
    xtb = CONFORMERS / taut / f"{conf_id}.xyz"
    return xtb, "GFN2-xTB"


def run_tddft(taut: str, conf_id: str, level: dict, solvent: str,
              prof: dict) -> dict:
    """한 조합(구조 x 수준 x 용매)의 TD-DFT."""
    td = prof["tddft"]
    tag = f"{taut}_{conf_id}_{level['id']}_{solvent}"
    outdir = OUTROOT / taut / conf_id / level["id"] / solvent
    outdir.mkdir(parents=True, exist_ok=True)
    ck = outdir / "result.json"

    done = load_checkpoint(ck)
    if done and done.get("ok"):
        b = done.get("brightest", {})
        print(f"  [건너뜀] {tag}  (lambda_max = {b.get('wavelength_nm')} nm)")
        return done

    xyz, geom_source = geometry_for(taut, conf_id)
    if not xyz.exists():
        print(f"  [건너뜀] {tag}: 구조 파일 없음 ({xyz})")
        return {"ok": False, "reason": "geometry missing"}

    use_pcm = solvent != "none"
    print(f"  [TD-DFT] {tag}  ({level['functional']}/{level['basis']}, "
          f"{td['n_states']} states, {geom_source} 구조)")

    rec: dict = {"tautomer": taut, "conf_id": conf_id,
                 "level_id": level["id"], "functional": level["functional"],
                 "basis": level["basis"], "solvent": solvent,
                 "pcm": use_pcm, "tda": td["tda"], "n_states": td["n_states"],
                 "geometry_source": geom_source, "geometry_file": str(xyz)}
    t0 = time.time()
    try:
        psi4.core.clean()
        psi4.core.clean_options()
        init_psi4(LOGS / f"psi4_td_{tag}.out", memory=prof["memory"],
                  nthreads=prof["nthreads"])

        g = read_xyz(xyz)
        mol = make_molecule(g.symbols, g.coords)

        options = {
            "basis": level["basis"],
            "scf_type": "df",
            "df_scf_guess": True,
            "e_convergence": 1e-8,
            "d_convergence": 1e-8,
            "maxiter": 250,
            "save_jk": True,
            "tdscf_maxiter": 120,
            "tdscf_r_convergence": 1e-4,
        }
        if use_pcm:
            options.update({"pcm": True, "pcm_scf_type": "total"})
        psi4.set_options(options)
        if use_pcm:
            # 수직 전이이므로 비평형 용매화
            psi4.pcm_helper(pcm_block(solvent, nonequilibrium=True,
                                      area=td.get("pcm_area", 1.0)))

        e_scf, wfn = psi4.energy(level["functional"], molecule=mol, return_wfn=True)
        t_scf = time.time() - t0
        nocc = wfn.nalpha()

        t1 = time.time()
        res = psi4.procrouting.response.scf_response.tdscf_excitations(
            wfn, states=td["n_states"], tda=td["tda"], triplets="NONE")
        t_td = time.time() - t1

        trans = parse_tdscf_results(res, nocc)
        bright = max(trans, key=lambda t: t["osc_strength"])
        rec.update({
            "ok": True,
            "scf_energy_hartree": float(e_scf),
            "n_basis": wfn.basisset().nbf(),
            "n_occupied": nocc,
            "scf_seconds": round(t_scf, 1),
            "tdscf_seconds": round(t_td, 1),
            "transitions": trans,
            "brightest": bright,
        })
        print(f"    완료 ({t_scf:.0f}s + {t_td:.0f}s)  "
              f"가장 센 전이 {bright['wavelength_nm']:.1f} nm "
              f"(f={bright['osc_strength']:.3f})  {bright['orbital_transitions_str']}")
    except Exception:                                     # noqa: BLE001
        err = traceback.format_exc()
        rec.update({"ok": False, "seconds": round(time.time() - t0, 1)})
        rec["failure"] = record_failure(f"tddft:{tag}", err,
                                        {"level": f"{level['functional']}/{level['basis']}",
                                         "solvent": solvent})
    save_checkpoint(ck, rec)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tautomers", nargs="*", default=["enolA", "enolB", "diketo"])
    ap.add_argument("--levels", nargs="*", default=None,
                    help="calc_config.json 의 level id. 생략하면 전부.")
    ap.add_argument("--solvents", nargs="*", default=None)
    ap.add_argument("--max-conformers", type=int, default=None)
    args = ap.parse_args()

    prof = load_profile()
    td = prof["tddft"]
    levels = td["levels"]
    if args.levels:
        levels = [l for l in levels if l["id"] in args.levels]
    solvents = args.solvents or td["solvents"]

    print(f"프로파일: {prof['_profile_name']}")
    print(f"수준: {[l['id'] for l in levels]}   용매: {solvents}   "
          f"상태 수: {td['n_states']}  TDA={td['tda']}")

    all_recs = []
    for taut in args.tautomers:
        sel_path = CONFORMERS / taut / "selected.json"
        if not sel_path.exists():
            print(f"[건너뜀] {taut}: selected.json 없음")
            continue
        sel = json.loads(sel_path.read_text(encoding="utf-8"))
        confs = sel["selected"]
        if args.max_conformers:
            confs = confs[:args.max_conformers]

        print(f"\n=== {taut}: 컨포머 {len(confs)} 개 x 수준 {len(levels)} "
              f"x 용매 {len(solvents)} = {len(confs)*len(levels)*len(solvents)} 계산 ===")
        for c in confs:
            for level in levels:
                for solv in solvents:
                    rec = run_tddft(taut, c["conf_id"], level, solv, prof)
                    rec["boltzmann_weight"] = c["boltzmann_weight"]
                    rec["rel_energy_kcalmol"] = c["rel_energy_kcalmol"]
                    all_recs.append(rec)
                    save_checkpoint(OUTROOT / "all_results.json",
                                    {"results": all_recs})

    ok = [r for r in all_recs if r.get("ok")]
    print(f"\n=== 완료: {len(ok)}/{len(all_recs)} 성공 ===")
    for r in ok:
        print(f"  {r['tautomer']:7s} {r['conf_id']:14s} {r['level_id']:18s} "
              f"{r['solvent']:8s} lambda_max={r['brightest']['wavelength_nm']:7.1f} nm "
              f"f={r['brightest']['osc_strength']:.3f}")
    fails = [r for r in all_recs if not r.get("ok")]
    if fails:
        print(f"\n실패 {len(fails)} 건 -> logs/failures.json 확인")
    return 0


if __name__ == "__main__":
    sys.exit(main())
