"""
03_test_single_point.py
-----------------------
가장 저렴한 "단일 구조" 테스트 계산 + 비용 측정.

목적
  1) 파이프라인 전체(xTB 최적화 -> Psi4 SCF -> TD-DFT -> PCM)가 실제로 도는지 확인
  2) 아보벤존(45원자) 크기에서 각 단계가 이 노트북에서 몇 초 걸리는지 측정해
     본 계산의 함수/기저셋을 근거 있게 고르기 위함 (요구사항 9번)

실행:  .\scripts\run.ps1 scripts\03_test_single_point.py
결과:  calculations/00_test/test_report.json , logs/psi4_test.out
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import (CALCULATIONS, EV_TO_NM, LOGS, STRUCTURES, read_xyz,
                       save_checkpoint)

import psi4

OUT = CALCULATIONS / "00_test"
OUT.mkdir(parents=True, exist_ok=True)

psi4.core.set_output_file(str(LOGS / "psi4_test.out"), False)
psi4.set_memory("6 GB")
psi4.set_num_threads(4)

PCM_ETHANOL = """
Units = Angstrom
Medium {
    SolverType = IEFPCM
    Solvent = Ethanol
    Nonequilibrium = True
}
Cavity {
    RadiiSet = Bondi
    Type = GePol
    Scaling = True
    Area = 0.3
    Mode = Implicit
}
"""


def make_mol(path: Path):
    g = read_xyz(path)
    return psi4.geometry(
        f"0 1\n{g.to_xyz_block()}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n"), g


def run_case(name: str, xyz: Path, functional: str, basis: str,
             pcm: bool, nstates: int, tda: bool = True) -> dict:
    """SCF + TD-DFT 한 세트를 돌리고 시간과 결과를 기록."""
    print(f"\n--- {name}: {functional}/{basis}  PCM={pcm}  states={nstates} ---")
    rec: dict = {"name": name, "functional": functional, "basis": basis,
                 "pcm": pcm, "nstates": nstates, "tda": tda}
    try:
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.set_output_file(str(LOGS / "psi4_test.out"), True)
        psi4.set_memory("6 GB")
        psi4.set_num_threads(4)

        mol, geom = make_mol(xyz)
        opts = {"basis": basis, "scf_type": "df", "df_scf_guess": True,
                "e_convergence": 1e-8, "d_convergence": 1e-8,
                "maxiter": 200, "save_jk": True, "print": 1}
        if pcm:
            opts.update({"pcm": True, "pcm_scf_type": "total"})
        psi4.set_options(opts)
        if pcm:
            psi4.pcm_helper(PCM_ETHANOL)

        nbf = psi4.core.BasisSet.build(mol, "ORBITAL", basis).nbf()
        rec["nbf"] = nbf
        print(f"    기저함수 개수 = {nbf}")

        t0 = time.time()
        e_scf, wfn = psi4.energy(functional, molecule=mol, return_wfn=True)
        t_scf = time.time() - t0
        rec.update({"scf_energy_hartree": e_scf, "scf_seconds": round(t_scf, 1)})
        print(f"    SCF  E = {e_scf:.8f} Eh   ({t_scf:.1f} s)")

        t0 = time.time()
        res = psi4.procrouting.response.scf_response.tdscf_excitations(
            wfn, states=nstates, tda=tda, triplets="NONE")
        t_td = time.time() - t0
        rec["tdscf_seconds"] = round(t_td, 1)

        trans = []
        for i, r in enumerate(res, start=1):
            ev = float(r["EXCITATION ENERGY"]) * 27.211386245988
            f_osc = float(r["OSCILLATOR STRENGTH (LEN)"])
            trans.append({"state": i, "energy_eV": round(ev, 4),
                          "wavelength_nm": round(EV_TO_NM / ev, 2),
                          "osc_strength": round(f_osc, 5)})
        rec["transitions"] = trans
        print(f"    TD-DFT ({t_td:.1f} s)")
        for t in trans:
            bar = "*" * int(min(40, t["osc_strength"] * 40))
            print(f"      S{t['state']:<2d} {t['energy_eV']:6.3f} eV  "
                  f"{t['wavelength_nm']:7.1f} nm  f={t['osc_strength']:.4f} {bar}")
        bright = max(trans, key=lambda t: t["osc_strength"])
        rec["brightest"] = bright
        print(f"    -> 가장 센 전이: {bright['wavelength_nm']:.1f} nm "
              f"(f={bright['osc_strength']:.4f})")
        rec["ok"] = True
    except Exception:                                  # noqa: BLE001
        rec["ok"] = False
        rec["error"] = traceback.format_exc()[-1500:]
        print("    [실패]", rec["error"].splitlines()[-1][:200])
    return rec


def main() -> int:
    # GFN2-xTB(ALPB ethanol) 로 최적화된 구조를 우선 사용한다.
    xyz = STRUCTURES / "enolA_xtbopt.xyz"
    if not xyz.exists():
        xyz = STRUCTURES / "enolA.xyz"
    if not xyz.exists():
        print("먼저 01_build_structures.py 를 실행하세요.")
        return 1
    print(f"사용 구조: {xyz.name}")

    report = {"structure": str(xyz.relative_to(xyz.parent.parent)),
              "note": "GFN2-xTB(ALPB ethanol) 최적화 구조에 대한 비용 측정용 단일점 계산",
              "cases": []}

    cases = [
        # (이름, 함수, 기저, PCM, 상태수)  - 싼 것부터
        ("T1_cheap_gas",   "B3LYP",     "6-31G",     False, 5),
        ("T2_cheap_pcm",   "B3LYP",     "6-31G",     True,  5),
        ("T3_pol_gas",     "B3LYP",     "6-31G(d)",  False, 8),
        ("T4_pol_pcm",     "B3LYP",     "6-31G(d)",  True,  8),
        ("T5_camb3lyp_pcm", "CAM-B3LYP", "6-31G(d)", True,  8),
    ]
    for name, func, basis, pcm, ns in cases:
        rec = run_case(name, xyz, func, basis, pcm, ns)
        report["cases"].append(rec)
        save_checkpoint(OUT / "test_report.json", report)   # 매 케이스마다 즉시 저장
        if not rec["ok"]:
            print("    (실패했지만 다음 케이스 계속 진행)")

    print("\n=== 비용 요약 ===")
    print(f"{'case':16s} {'nbf':>5s} {'SCF s':>8s} {'TDDFT s':>9s} {'lambda_max nm':>14s}")
    for c in report["cases"]:
        if c.get("ok"):
            print(f"{c['name']:16s} {c['nbf']:5d} {c['scf_seconds']:8.1f} "
                  f"{c['tdscf_seconds']:9.1f} {c['brightest']['wavelength_nm']:14.1f}")
        else:
            print(f"{c['name']:16s} {'-':>5s} {'FAIL':>8s}")
    save_checkpoint(OUT / "test_report.json", report)
    print(f"\n저장: {(OUT / 'test_report.json')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
