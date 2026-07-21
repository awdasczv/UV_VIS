"""
05b_tddft_orca.py
-----------------
ORCA 6.1 로 TD-DFT 본 계산 (요구사항 10, 11번).

Psi4 판(05_tddft.py)과 결과 형식이 같아서, 이후 단계(06 스펙트럼, 07 보고서)는
그대로 쓸 수 있다. 엔진만 갈아끼운 것이다.

  - 대상: conformers/<토토머>/selected.json 의 대표 컨포머
  - 조합: (구조 × 이론수준 × 용매) 마다 계산하고 즉시 체크포인트 저장.
    이미 끝난 조합은 자동으로 건너뛴다 -> 중간에 꺼도 처음부터 안 해도 된다.
  - 실패 시 원인을 분류해 logs/failures.json 에 기록하고 다음 조합으로 진행.

실행:
  .\scripts\run.ps1 scripts\05b_tddft_orca.py
  .\scripts\run.ps1 scripts\05b_tddft_orca.py --levels b3lyp_def2svp --solvents ethanol
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orca_common import (build_input, classify_orca_failure, find_orca,
                         mpi_available, parse_output, run_orca)
from qc_common import (CALCULATIONS, CONFORMERS, INPUTS, LOGS, load_checkpoint,
                       read_xyz, save_checkpoint)

OUTROOT = CALCULATIONS / "02_tddft_orca"

# ORCA 의 CPCM 용매 이름
SOLVENT_NAMES = {"ethanol": "Ethanol", "methanol": "Methanol", "water": "Water",
                 "acetonitrile": "Acetonitrile", "cyclohexane": "Cyclohexane",
                 "dmso": "DMSO"}


def load_cfg() -> dict:
    return json.loads((INPUTS / "calc_config.json").read_text(encoding="utf-8"))


def record_failure(tag: str, info: dict, extra: dict | None = None) -> dict:
    path = LOGS / "failures.json"
    data = load_checkpoint(path) or {"failures": []}
    entry = {"tag": tag, "engine": "ORCA", **info}
    if extra:
        entry.update(extra)
    data["failures"].append(entry)
    save_checkpoint(path, data)
    print(f"    [실패 분류] {info.get('code')}")
    print(f"    [대응 방안] {info.get('remedy')}")
    return entry


def geometry_for(taut: str, conf_id: str) -> tuple[Path, str]:
    """DFT 최적화 구조가 있으면 그것, 없으면 GFN2-xTB 구조."""
    dft = CALCULATIONS / "01_dft_opt" / taut / conf_id / "optimized.xyz"
    if dft.exists():
        return dft, "DFT"
    return CONFORMERS / taut / f"{conf_id}.xyz", "GFN2-xTB"


def run_one(taut: str, conf_id: str, level: dict, solvent: str,
            orca_cfg: dict, weight: float, rel_e: float) -> dict:
    td = orca_cfg["tddft"]
    tag = f"{taut}_{conf_id}_{level['id']}_{solvent}"
    outdir = OUTROOT / taut / conf_id / level["id"] / solvent
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

    g = read_xyz(xyz)
    solv_name = SOLVENT_NAMES.get(solvent) if solvent != "none" else None
    print(f"  [TD-DFT] {tag}")
    print(f"           {level['functional']}/{level['basis']}, "
          f"{td['n_states']}상태, 용매={solv_name or '없음'}, 구조={geom_source}")

    inp = build_input(
        g.symbols, g.coords,
        functional=level["functional"], basis=level["basis"],
        nstates=td["n_states"], tda=td["tda"], solvent=solv_name,
        nprocs=orca_cfg["nprocs"], maxcore_mb=orca_cfg["maxcore_mb"],
        rijcosx=orca_cfg["rijcosx"], aux_basis=orca_cfg["aux_basis"],
        comment=f"{taut} / {conf_id} / {level['id']} / solvent={solvent}",
    )

    t0 = time.time()
    out_path, seconds = run_orca(inp, outdir, name="td")
    rec: dict = {
        "engine": "ORCA 6.1.1",
        "tautomer": taut, "conf_id": conf_id,
        "level_id": level["id"], "functional": level["functional"],
        "basis": level["basis"], "solvent": solvent,
        "pcm": solvent != "none", "solvent_model": "CPCM" if solv_name else None,
        "tda": td["tda"], "n_states": td["n_states"],
        "rijcosx": orca_cfg["rijcosx"], "nprocs": orca_cfg["nprocs"],
        "geometry_source": geom_source, "geometry_file": str(xyz),
        "boltzmann_weight": weight, "rel_energy_kcalmol": rel_e,
        "wall_seconds": round(seconds, 1),
    }

    parsed = parse_output(out_path)
    if not parsed.get("terminated_normally") or not parsed.get("transitions"):
        rec["ok"] = False
        rec["failure"] = record_failure(
            f"tddft-orca:{tag}", classify_orca_failure(out_path),
            {"level": f"{level['functional']}/{level['basis']}",
             "solvent": solvent, "output": str(out_path)})
        save_checkpoint(ck, rec)
        return rec

    rec.update({
        "ok": True,
        "scf_energy_hartree": parsed["scf_energy_hartree"],
        "n_basis": parsed["n_basis"],
        "n_occupied": parsed["n_occupied"],
        "transitions": parsed["transitions"],
        "brightest": parsed["brightest"],
    })
    b = parsed["brightest"]
    lo = min(t["wavelength_nm"] for t in parsed["transitions"])
    rec["shortest_wavelength_nm"] = lo
    rec["covers_200nm"] = lo <= 200.0
    print(f"           완료 {seconds:.0f}초 | 최강 {b['wavelength_nm']:.1f} nm "
          f"(f={b['osc_strength']:.3f}) | {b['orbital_transitions_str'][:44]}")
    print(f"           최단파장 {lo:.1f} nm "
          f"({'200 nm 커버 OK' if lo <= 200 else '200 nm 미달'})")
    save_checkpoint(ck, rec)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tautomers", nargs="*", default=["enolA", "enolB", "diketo"])
    ap.add_argument("--levels", nargs="*", default=None)
    ap.add_argument("--solvents", nargs="*", default=None)
    ap.add_argument("--max-conformers", type=int, default=None)
    args = ap.parse_args()

    cfg = load_cfg()
    orca_cfg = cfg["orca"]
    td = orca_cfg["tddft"]
    levels = [l for l in td["levels"]
              if not args.levels or l["id"] in args.levels]
    solvents = args.solvents or td["solvents"]

    print(f"엔진   : ORCA  ({find_orca()})")
    print(f"병렬   : {orca_cfg['nprocs']} 프로세스, MPI 사용 가능 = {mpi_available()}")
    if orca_cfg["nprocs"] > 1 and not mpi_available():
        print("  [경고] MPI 를 못 찾았습니다. 직렬로 돌게 되어 약 2배 느려집니다.")
    print(f"수준   : {[l['id'] for l in levels]}")
    print(f"용매   : {solvents}")
    print(f"상태 수: {td['n_states']}  (TDA={td['tda']}, RIJCOSX={orca_cfg['rijcosx']})")

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

        n = len(confs) * len(levels) * len(solvents)
        print(f"\n=== {taut}: 컨포머 {len(confs)} × 수준 {len(levels)} "
              f"× 용매 {len(solvents)} = {n} 계산 ===")
        for c in confs:
            for level in levels:
                for solv in solvents:
                    rec = run_one(taut, c["conf_id"], level, solv, orca_cfg,
                                  c.get("weight_normalized",
                                        c.get("boltzmann_weight", 1.0)),
                                  c.get("rel_energy_kcalmol", 0.0))
                    all_recs.append(rec)
                    save_checkpoint(OUTROOT / "all_results.json",
                                    {"results": all_recs})

    ok = [r for r in all_recs if r.get("ok")]
    print(f"\n=== 완료: {len(ok)}/{len(all_recs)} 성공 ===")
    print(f"{'토토머':8s} {'컨포머':14s} {'수준':18s} {'용매':8s} "
          f"{'lambda_max':>10s} {'f':>7s} {'초':>7s}")
    for r in ok:
        print(f"{r['tautomer']:8s} {r['conf_id']:14s} {r['level_id']:18s} "
              f"{r['solvent']:8s} {r['brightest']['wavelength_nm']:10.1f} "
              f"{r['brightest']['osc_strength']:7.3f} {r['wall_seconds']:7.0f}")
    fails = [r for r in all_recs if not r.get("ok")]
    if fails:
        print(f"\n실패 {len(fails)} 건 -> logs/failures.json 확인")
    return 0


if __name__ == "__main__":
    sys.exit(main())
