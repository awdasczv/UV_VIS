"""
08_lowcost_baseline.py
----------------------
"MINDO/3–TDA 대체" 저비용 기준선 계산 (요구사항 15번의 비교 대상).

왜 필요한가
  아보벤존에 대해 발표된 MINDO/3–TDA λmax 를 문헌에서 찾지 못했다.
  (inputs/experimental_reference.json 의 MINDO3_TDA.status = "NOT_FOUND")
  따라서 "저비용 계산 대비 얼마나 좋아졌는가" 를 말하려면 기준선을 직접 만들어야 한다.

무엇을 기준선으로 삼는가
  MINDO/3–TDA 는 (1) 반경험적 해밀토니안 + (2) TDA(= CIS) 형태의 들뜬상태 처리다.
  이 조합에 가장 가까우면서 우리 엔진으로 재현 가능한 것은
      HF / 최소기저(STO-3G) + CIS(=TDA)
  이다. 전자상관도 없고 기저셋도 최소라, 반경험적 수준의 오차 크기를 보여준다.
  비교를 위해 조금 더 나은 HF/6-31G 도 함께 계산한다.

  ※ 이것은 MINDO/3 "그 자체"가 아니다. 보고서에는 반드시
    "문헌에 MINDO/3 값이 없어 동급의 저비용 기준선을 자체 계산했다" 고 명시한다.

실행:  .\scripts\run.ps1 scripts\08_lowcost_baseline.py
"""
from __future__ import annotations

import json
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import psi4
from psi4_helpers import init_psi4, make_molecule, parse_tdscf_results
from qc_common import CALCULATIONS, INPUTS, LOGS, STRUCTURES, read_xyz, save_checkpoint

OUT = CALCULATIONS / "00_test"
OUT.mkdir(parents=True, exist_ok=True)

# (id, 방법, 기저셋, 설명)
BASELINES = [
    ("hf_sto3g_cis", "HF", "STO-3G",
     "HF/STO-3G + CIS(TDA). 반경험적 MINDO/3-TDA 와 가장 성격이 비슷한 저비용 조합."),
    ("hf_631g_cis", "HF", "6-31G",
     "HF/6-31G + CIS(TDA). 기저셋만 조금 키운 저비용 기준선."),
]
NSTATES = 8


def run_one(bid: str, method: str, basis: str, xyz: Path) -> dict:
    print(f"\n--- {bid}: {method}/{basis} + CIS(TDA), {NSTATES} 상태 ---")
    rec = {"id": bid, "method": method, "basis": basis, "nstates": NSTATES,
           "tda": True, "solvent": "none", "geometry": str(xyz)}
    t0 = time.time()
    try:
        psi4.core.clean()
        psi4.core.clean_options()
        init_psi4(LOGS / f"psi4_baseline_{bid}.out", memory="6 GB", nthreads=4)

        g = read_xyz(xyz)
        mol = make_molecule(g.symbols, g.coords)
        psi4.set_options({"basis": basis, "scf_type": "df", "df_scf_guess": True,
                          "e_convergence": 1e-8, "d_convergence": 1e-8,
                          "maxiter": 250, "save_jk": True,
                          "tdscf_maxiter": 120, "tdscf_r_convergence": 1e-4,
                          "reference": "rhf"})
        e, wfn = psi4.energy(method, molecule=mol, return_wfn=True)
        res = psi4.procrouting.response.scf_response.tdscf_excitations(
            wfn, states=NSTATES, tda=True, triplets="NONE")
        trans = parse_tdscf_results(res, wfn.nalpha())
        bright = max(trans, key=lambda t: t["osc_strength"])
        rec.update({"ok": True, "n_basis": wfn.basisset().nbf(),
                    "seconds": round(time.time() - t0, 1),
                    "transitions": trans,
                    "lambda_max_nm": bright["wavelength_nm"],
                    "osc_strength": bright["osc_strength"],
                    "level": f"{method}/{basis} + CIS(TDA), 기체상"})
        print(f"    기저함수 {rec['n_basis']}  ({rec['seconds']:.0f}s)")
        for t in trans:
            bar = "*" * int(min(30, t["osc_strength"] * 30))
            print(f"      S{t['state']:<2d} {t['energy_eV']:6.3f} eV "
                  f"{t['wavelength_nm']:7.1f} nm f={t['osc_strength']:.4f} {bar}")
        print(f"    -> lambda_max = {bright['wavelength_nm']:.1f} nm "
              f"(f={bright['osc_strength']:.3f})")
    except Exception:                                    # noqa: BLE001
        rec.update({"ok": False, "error": traceback.format_exc()[-1200:]})
        print("    실패:", rec["error"].splitlines()[-1][:160])
    return rec


def main() -> int:
    xyz = STRUCTURES / "enolA_xtbopt.xyz"
    if not xyz.exists():
        print("먼저 구조를 준비하세요 (01_build_structures.py + xtb 최적화).")
        return 1

    ref = json.loads((INPUTS / "experimental_reference.json").read_text(encoding="utf-8"))
    exp = ref["experimental"]["enol_band"]["primary_target"]["lambda_max_nm"]

    report = {
        "purpose": "MINDO/3-TDA 문헌값이 없어 동급의 저비용 기준선을 자체 계산한 것",
        "caveat": "이것은 MINDO/3 자체가 아니다. 반경험적 수준의 오차 크기를 보여주는 대용물이다.",
        "structure": str(xyz.name),
        "tautomer": "enolA",
        "experimental_nm": exp,
        "runs": [],
    }
    for bid, method, basis, desc in BASELINES:
        rec = run_one(bid, method, basis, xyz)
        rec["description"] = desc
        report["runs"].append(rec)
        save_checkpoint(OUT / "baseline_lowcost.json", report)

    ok = [r for r in report["runs"] if r.get("ok")]
    if ok:
        # 07_report.py 가 읽는 대표값: 가장 싼(첫 번째) 성공 결과
        primary = ok[0]
        report["lambda_max_nm"] = primary["lambda_max_nm"]
        report["level"] = primary["level"]
        save_checkpoint(OUT / "baseline_lowcost.json", report)

        print("\n=== 저비용 기준선 요약 ===")
        print(f"{'id':16s} {'nbf':>5s} {'초':>7s} {'lambda nm':>10s} "
              f"{'실험 대비':>10s}")
        for r in ok:
            print(f"{r['id']:16s} {r['n_basis']:5d} {r['seconds']:7.0f} "
                  f"{r['lambda_max_nm']:10.1f} {r['lambda_max_nm']-exp:+10.1f}")
        print(f"\n실험값(에탄올, 에놀): {exp} nm")
    return 0


if __name__ == "__main__":
    sys.exit(main())
