"""
03b_pcm_cost_probe.py
---------------------
PCM 이 TD-DFT 를 왜 그렇게 느리게 만드는지 원인을 짚고 최적화한다.

측정 결과 (아보벤존 45원자, B3LYP/6-31G = 251 기저함수, 5 상태, 4코어):
    용매 미적용  Davidson 1회 ≈  25 초
    PCM(Area=0.3) Davidson 1회 ≈ 440 초   -> 약 17배

원인 가설: PCM 은 Davidson 시행벡터마다 공동 표면의 겉보기 전하를 다시 풀어야 한다.
          비용은 tessera(표면 조각) 개수에 좌우되고, tessera 개수는 Cavity Area 로 정해진다.
          Area = 0.3 A^2 는 매우 촘촘한 격자다.

이 스크립트는 Area 를 바꿔가며
  - tessera 개수
  - SCF 시간 / TD-DFT 시간
  - 얻어지는 lambda_max
를 비교해서, 정확도를 얼마나 잃지 않고 얼마나 빨라지는지 정량화한다.

실행:  .\scripts\run.ps1 scripts\03b_pcm_cost_probe.py
"""
from __future__ import annotations

import re
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import psi4
from psi4_helpers import init_psi4, make_molecule, parse_tdscf_results, pcm_block
from qc_common import CALCULATIONS, LOGS, STRUCTURES, read_xyz, save_checkpoint

OUT = CALCULATIONS / "00_test"
OUT.mkdir(parents=True, exist_ok=True)

FUNCTIONAL = "B3LYP"
BASIS = "6-31G"
NSTATES = 5
# Area = 0.3 은 이미 03_test_single_point.py 에서 재봤다 (Davidson 1회 약 440초).
# 여기서는 더 성긴 격자만 시험한다.
AREAS = [1.0, 2.0]           # A^2 per tessera
INCLUDE_GAS = False          # 기체상도 T1 에서 이미 측정됨 (SCF 76s + TDDFT 149s)


def count_tesserae(logfile: Path) -> int | None:
    """psi4 출력에서 PCM 공동의 tessera(유한요소) 개수를 찾는다."""
    try:
        text = logfile.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    for pat in (r"Number of finite elements\s*=\s*(\d+)",
                r"finite elements\s*=\s*(\d+)",
                r"Cavity size:\s*(\d+)"):
        m = re.search(pat, text)
        if m:
            return int(m.group(1))
    return None


def run(area: float | None) -> dict:
    tag = "gas" if area is None else f"area{area}"
    log = LOGS / f"psi4_pcmcost_{tag}.out"
    print(f"\n--- {tag} ---")
    rec: dict = {"tag": tag, "area": area, "functional": FUNCTIONAL,
                 "basis": BASIS, "nstates": NSTATES}
    try:
        psi4.core.clean()
        psi4.core.clean_options()
        init_psi4(log, memory="6 GB", nthreads=4)

        g = read_xyz(STRUCTURES / "enolA_xtbopt.xyz")
        mol = make_molecule(g.symbols, g.coords)

        opts = {"basis": BASIS, "scf_type": "df", "df_scf_guess": True,
                "e_convergence": 1e-8, "d_convergence": 1e-8,
                "maxiter": 250, "save_jk": True,
                "tdscf_maxiter": 120, "tdscf_r_convergence": 1e-4}
        if area is not None:
            opts.update({"pcm": True, "pcm_scf_type": "total"})
        psi4.set_options(opts)
        if area is not None:
            psi4.pcm_helper(pcm_block("ethanol", nonequilibrium=True, area=area))

        t0 = time.time()
        e, wfn = psi4.energy(FUNCTIONAL, molecule=mol, return_wfn=True)
        t_scf = time.time() - t0

        t1 = time.time()
        res = psi4.procrouting.response.scf_response.tdscf_excitations(
            wfn, states=NSTATES, tda=True, triplets="NONE")
        t_td = time.time() - t1

        trans = parse_tdscf_results(res, wfn.nalpha())
        bright = max(trans, key=lambda t: t["osc_strength"])
        rec.update({"ok": True, "scf_seconds": round(t_scf, 1),
                    "tdscf_seconds": round(t_td, 1),
                    "n_tesserae": count_tesserae(log),
                    "scf_energy": float(e),
                    "brightest_nm": bright["wavelength_nm"],
                    "brightest_f": bright["osc_strength"],
                    "transitions": trans})
        print(f"    tessera={rec['n_tesserae']}  SCF={t_scf:.0f}s  "
              f"TDDFT={t_td:.0f}s  lambda_max={bright['wavelength_nm']:.1f} nm "
              f"(f={bright['osc_strength']:.3f})")
    except Exception:                                    # noqa: BLE001
        rec.update({"ok": False, "error": traceback.format_exc()[-1000:]})
        print("    실패:", rec["error"].splitlines()[-1][:160])
    return rec


def main() -> int:
    report = {"note": "PCM cavity Area 에 따른 TD-DFT 비용/정확도 절충 측정",
              "runs": []}
    for area in ([None] if INCLUDE_GAS else []) + AREAS:
        report["runs"].append(run(area))
        save_checkpoint(OUT / "pcm_cost_probe.json", report)

    print("\n=== 요약 ===")
    print(f"{'조건':10s} {'tessera':>8s} {'SCF s':>8s} {'TDDFT s':>9s} "
          f"{'lambda nm':>10s} {'gas 대비 이동':>12s}")
    gas = next((r for r in report["runs"] if r["tag"] == "gas" and r.get("ok")), None)
    # 기체상을 다시 안 돌렸으면 T1 에서 측정한 값을 기준으로 쓴다
    gas_nm = gas["brightest_nm"] if gas else 315.7
    for r in report["runs"]:
        if not r.get("ok"):
            print(f"{r['tag']:10s} 실패")
            continue
        shift = r["brightest_nm"] - gas_nm
        print(f"{r['tag']:10s} {str(r['n_tesserae']):>8s} {r['scf_seconds']:8.1f} "
              f"{r['tdscf_seconds']:9.1f} {r['brightest_nm']:10.1f} {shift:+12.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
