"""
12_functional_benchmark.py
--------------------------
함수(functional) 벤치마크: 구조와 기저셋을 고정하고 함수만 바꿔 TD-DFT 를 돈다.

왜 필요한가 (DHHB 사례)
  DHHB 의 장파장 UVA 전이는 도너(NEt2)->억셉터(C=O) CT 성격을 강하게 시사하며,
  B3LYP 는 +37 nm 적색, CAM-B3LYP 는 -36 nm 청색으로 실험(354 nm)을 크게 벗어났다.
  두 함수가 실험을 감싸므로, 장거리 정확교환 비율이 그 사이인 함수들을 스크리닝해
  이 발색단 유형에 맞는 함수를 찾는다. 구조(c000 DFT-opt)와 기저셋(6-31+G(d)),
  용매(CPCM ethanol)를 고정해야 함수 차이만 분리해서 볼 수 있다.

무엇을 돌리는가
  함수 5종 x (TDA + full TD-DFT) = 10 계산. 전부 NTO(자연 전이 오비탈) 생성.
  NTO 로 'S1 번호'가 아니라 '같은 hole->particle 성격의 전이'를 함수 간 추적한다
  (함수가 바뀌면 상태 순서가 바뀌거나 다른 전이와 섞일 수 있기 때문).

실행:
  .\scripts\run.ps1 -Molecule dhhb scripts\12_functional_benchmark.py
  .\scripts\run.ps1 -Molecule dhhb scripts\12_functional_benchmark.py --functionals pbe0 m062x --modes tda
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orca_common import (build_input, classify_orca_failure, find_orca,
                         mpi_available, parse_nto, parse_output)
from qc_common import (CALCULATIONS, CONFORMERS, INPUTS, LOGS, MOL_CONFIG,
                       load_checkpoint, read_xyz, save_checkpoint)

OUTROOT = CALCULATIONS / "03_functional_benchmark"

# 벤치마크 함수 목록. hf_exchange 는 보고서용 메타데이터(함수의 정확교환 구성).
FUNCTIONALS = {
    "b3lyp":    {"orca": "B3LYP",     "hf_exchange": "전역 20%"},
    "pbe0":     {"orca": "PBE0",      "hf_exchange": "전역 25%"},
    "m062x":    {"orca": "M062X",     "hf_exchange": "전역 54%"},
    "camb3lyp": {"orca": "CAM-B3LYP", "hf_exchange": "범위분리 19% -> 장거리 65%"},
    "wb97xd4":  {"orca": "WB97X-D4",  "hf_exchange": "범위분리 ~17% -> 장거리 100%"},
}

BRIGHT_F = 0.10          # '밝은 전이' 판정 문턱 (11_component_viewer 와 동일)
NTO_STATES = "1,2,3,4,5,6"   # UVA CT 밴드는 모든 후보 함수에서 최저 6개 안에 있다


def pick_uva_band(transitions: list[dict]) -> dict | None:
    """밝은(f>=0.10) 전이 중 가장 장파장의 것 = UVA 밴드 후보. NTO 로 재확인한다."""
    bright = [t for t in transitions if t["osc_strength"] >= BRIGHT_F]
    return max(bright, key=lambda t: t["wavelength_nm"]) if bright else None


def run_one(species: str, conf_id: str, xyz: Path, func_id: str, mode: str,
            basis: str, solvent_name: str, nstates: int, orca_cfg: dict,
            tddft_extra: tuple = ()) -> dict:
    from orca_common import run_orca            # 지연 import (테스트 편의)
    func = FUNCTIONALS[func_id]
    tag = f"{conf_id}_{func_id}_{mode}"
    outdir = OUTROOT / conf_id / func_id / mode
    ck = outdir / "result.json"

    done = load_checkpoint(ck)
    if done and done.get("ok"):
        u = done.get("uva_band_transition") or {}
        print(f"  [건너뜀] {tag}  (UVA band = {u.get('wavelength_nm')} nm)")
        return done

    g = read_xyz(xyz)
    tda = mode == "tda"
    print(f"  [벤치마크] {tag}")
    print(f"             {func['orca']}/{basis}, {'TDA' if tda else 'full TD-DFT'}, "
          f"{nstates}상태, CPCM({solvent_name}), NTO {NTO_STATES}")

    inp = build_input(
        g.symbols, g.coords,
        functional=func["orca"], basis=basis,
        nstates=nstates, tda=tda, solvent=solvent_name,
        nprocs=orca_cfg["nprocs"], maxcore_mb=orca_cfg["maxcore_mb"],
        rijcosx=orca_cfg["rijcosx"], aux_basis=orca_cfg["aux_basis"],
        nto_states=NTO_STATES, tddft_extra=tddft_extra,
        comment=f"functional benchmark / {species} {conf_id} / "
                f"{func_id} / {mode} / CPCM({solvent_name})"
                + (f" / extra={list(tddft_extra)}" if tddft_extra else ""),
    )

    out_path, seconds = run_orca(inp, outdir, name="td")
    rec: dict = {
        "engine": "ORCA 6.1.1",
        "species": species, "conf_id": conf_id,
        "functional_id": func_id, "functional": func["orca"],
        "hf_exchange": func["hf_exchange"],
        "basis": basis, "mode": mode, "tda": tda,
        "solvent": solvent_name.lower(), "solvent_model": "CPCM",
        "n_states": nstates, "nto_states": NTO_STATES,
        "geometry_file": str(xyz), "geometry_source": "DFT",
        "wall_seconds": round(seconds, 1),
    }

    parsed = parse_output(out_path)
    if not parsed.get("terminated_normally") or not parsed.get("transitions"):
        rec["ok"] = False
        info = classify_orca_failure(out_path)
        rec["failure"] = info
        data = load_checkpoint(LOGS / "failures.json") or {"failures": []}
        data["failures"].append({"tag": f"funcbench:{tag}", "engine": "ORCA", **info,
                                 "output": str(out_path)})
        save_checkpoint(LOGS / "failures.json", data)
        print(f"    [실패] {info.get('code')} -> {info.get('remedy')}")
        save_checkpoint(ck, rec)
        return rec

    uva = pick_uva_band(parsed["transitions"])
    rec.update({
        "ok": True,
        "scf_energy_hartree": parsed["scf_energy_hartree"],
        "n_basis": parsed["n_basis"],
        "n_occupied": parsed["n_occupied"],
        "transitions": parsed["transitions"],
        "brightest": parsed["brightest"],
        "uva_band_transition": uva,
        "nto": parse_nto(out_path),
    })
    u = uva or {}
    print(f"             완료 {seconds:.0f}초 | UVA band {u.get('wavelength_nm', '?')} nm "
          f"(f={u.get('osc_strength', 0):.3f}, S{u.get('state', '?')}) | "
          f"{(u.get('orbital_transitions_str') or '')[:44]}")
    save_checkpoint(ck, rec)
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default=None,
                    help="기본값: config.json 의 principal_species")
    ap.add_argument("--conf", default=None,
                    help="기본값: selected.json 의 첫 번째(최저에너지) 컨포머")
    ap.add_argument("--functionals", nargs="*", default=list(FUNCTIONALS),
                    choices=list(FUNCTIONALS))
    ap.add_argument("--modes", nargs="*", default=["tda", "full"],
                    choices=["tda", "full"])
    ap.add_argument("--basis", default="6-31+G(d)")
    ap.add_argument("--solvent", default="Ethanol")
    ap.add_argument("--nstates", type=int, default=22)
    ap.add_argument("--tddft-extra", nargs="*", default=[],
                    help="%%tddft 블록에 추가할 줄들. 예) --tddft-extra "
                         "'maxiter 300' 'solver bhp22'  (M06-2X Davidson 미수렴 재시도용)")
    args = ap.parse_args()

    cfg = json.loads((INPUTS / "calc_config.json").read_text(encoding="utf-8"))
    orca_cfg = cfg["orca"]
    mol_cfg = json.loads(MOL_CONFIG.read_text(encoding="utf-8"))

    species = args.species or mol_cfg["principal_species"]
    if args.conf:
        conf_id = args.conf
    else:
        sel = json.loads((CONFORMERS / species / "selected.json")
                         .read_text(encoding="utf-8"))
        conf_id = sel["selected"][0]["conf_id"]

    xyz = CALCULATIONS / "01_dft_opt" / species / conf_id / "optimized.xyz"
    if not xyz.exists():
        print(f"[중단] DFT 최적화 구조가 없습니다: {xyz}")
        print("       04b_dft_optimize_orca.py 를 먼저 실행하세요.")
        return 1

    print(f"엔진   : ORCA  ({find_orca()})   MPI = {mpi_available()}")
    print(f"대상   : {species} / {conf_id}  (DFT-opt 구조 고정)")
    print(f"수준   : {args.functionals} x {args.modes} / {args.basis} / "
          f"CPCM({args.solvent}) / {args.nstates}상태")
    n = len(args.functionals) * len(args.modes)
    print(f"계산 수: {n}  (TDA 1건 실측 ~45-70분, full 은 그 ~2배. 체크포인트로 재개 가능)\n")

    recs = []
    # TDA 전체를 먼저 (빠른 쪽부터 전 함수의 경향을 확보), 그다음 full.
    for mode in args.modes:
        for func_id in args.functionals:
            rec = run_one(species, conf_id, xyz, func_id, mode,
                          args.basis, args.solvent, args.nstates, orca_cfg,
                          tddft_extra=tuple(args.tddft_extra))
            recs.append(rec)
            save_checkpoint(OUTROOT / "summary.json",
                            {"species": species, "conf_id": conf_id,
                             "basis": args.basis, "solvent": args.solvent.lower(),
                             "results": recs})

    ok = [r for r in recs if r.get("ok")]
    print(f"\n=== 완료 {len(ok)}/{len(recs)} ===")
    print(f"{'함수':12s} {'모드':5s} {'UVA band':>9s} {'f':>7s} {'전역최강':>9s} {'초':>7s}")
    for r in ok:
        u = r.get("uva_band_transition") or {}
        b = r.get("brightest") or {}
        print(f"{r['functional']:12s} {r['mode']:5s} "
              f"{u.get('wavelength_nm', 0):9.1f} {u.get('osc_strength', 0):7.3f} "
              f"{b.get('wavelength_nm', 0):9.1f} {r['wall_seconds']:7.0f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
