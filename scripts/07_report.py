"""
07_report.py
------------
최종 분석 보고서(results/report.md) 의 정량적 부분을 자동 생성한다 (요구사항 15번).

자동으로 계산하는 것
  - 토토머별/수준별/용매별 계산 lambda_max 와 실험값 대비 오차 (nm, eV, %)
  - 컨포머 평균(앙상블) vs 최저에너지 단일 구조의 차이
  - 용매 적용 vs 미적용의 차이
  - 함수/기저셋에 따른 차이
  - 저비용 기준선 대비 개선 정도
  - 실패 기록 요약

서술 부분(해석, 한계, 다음 개선점)은 이 파일이 만든 표를 근거로 사람이 덧붙인다.

실행:  .\scripts\run.ps1 scripts\07_report.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import (CALCULATIONS, CONFORMERS, EV_TO_NM, INPUTS, LOGS,
                       RESULTS, load_checkpoint)

TAUT_LABEL = {"enolA": "킬레이트 에놀 A", "enolB": "킬레이트 에놀 B",
              "diketo": "디케토"}


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def brightest_of(recs: list[dict], weighted: bool = True) -> dict | None:
    """
    한 그룹(같은 토토머·수준·용매)의 대표 lambda_max.
      weighted=True  : 볼츠만 가중치로 컨포머별 lambda_max 를 평균
      weighted=False : 가장 안정한 컨포머의 값
    """
    if not recs:
        return None
    if not weighted:
        r = min(recs, key=lambda x: x.get("rel_energy_kcalmol") or 0.0)
        return {"lambda_nm": r["brightest"]["wavelength_nm"],
                "f": r["brightest"]["osc_strength"],
                "n_conf": 1,
                "orbitals": r["brightest"].get("orbital_transitions_str", "")}
    w = np.array([r.get("boltzmann_weight") or 1.0 for r in recs], float)
    w = w / w.sum()
    lam = np.array([r["brightest"]["wavelength_nm"] for r in recs], float)
    f = np.array([r["brightest"]["osc_strength"] for r in recs], float)
    top = recs[int(np.argmax(w))]
    return {"lambda_nm": float((w * lam).sum()), "f": float((w * f).sum()),
            "n_conf": len(recs),
            "orbitals": top["brightest"].get("orbital_transitions_str", "")}


def err_row(calc_nm: float, exp_nm: float) -> tuple[float, float, float]:
    """(nm 오차, eV 오차, % 오차)"""
    d_nm = calc_nm - exp_nm
    d_ev = EV_TO_NM / calc_nm - EV_TO_NM / exp_nm
    return d_nm, d_ev, 100.0 * d_nm / exp_nm


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    ref = load(INPUTS / "experimental_reference.json")
    cfg = load(INPUTS / "calc_config.json")
    td = load_checkpoint(CALCULATIONS / "02_tddft" / "all_results.json")
    if not td:
        print("TD-DFT 결과가 없습니다. 먼저 05_tddft.py 를 실행하세요.")
        return 1
    results = [r for r in td["results"] if r.get("ok")]

    exp_enol = ref["experimental"]["enol_band"]["primary_target"]["lambda_max_nm"]
    exp_keto = ref["experimental"]["diketo_band"]["primary_target"]["lambda_max_nm"]
    exp_for = {"enolA": exp_enol, "enolB": exp_enol, "diketo": exp_keto}

    tauts = sorted({r["tautomer"] for r in results})
    levels = sorted({r["level_id"] for r in results})
    solvents = sorted({r["solvent"] for r in results})

    L = []
    A = L.append
    A("# 아보벤존 TD-DFT UV–Vis 계산 — 최종 분석 보고서")
    A("")
    A(f"생성 시각: {datetime.now():%Y-%m-%d %H:%M}  ")
    A(f"프로파일: `{cfg['active_profile']}`  ")
    A(f"성공한 TD-DFT 계산: {len(results)} 건")
    A("")
    A("## 1. 실험 기준값")
    A("")
    A("| 밴드 | 용매 | 실험 λmax (nm) | ε (M⁻¹cm⁻¹) | 출처 |")
    A("|---|---|---|---|---|")
    for v in ref["experimental"]["enol_band"]["values"]:
        A(f"| 에놀 | {v['solvent']} | {v['lambda_max_nm']} | "
          f"{v['epsilon_M-1cm-1'] or '-'} | {v['ref']} |")
    for v in ref["experimental"]["diketo_band"]["values"]:
        A(f"| 디케토 | {v['solvent']} | {v['lambda_max_nm']} | - | {v['ref']} |")
    A("")
    A(f"**주요 비교 대상**: 에놀 {exp_enol} nm (에탄올), 디케토 {exp_keto} nm")
    A("")

    # ---------------------------------------------- 2. 전체 결과 표
    A("## 2. 계산된 최대흡수파장과 실험값의 오차")
    A("")
    A("컨포머 앙상블(볼츠만 가중) 기준.")
    A("")
    A("| 토토머 | 수준 | 용매 | 계산 λmax (nm) | f | 실험 (nm) | 오차 (nm) | 오차 (eV) | 오차 (%) | 주요 전이 |")
    A("|---|---|---|---|---|---|---|---|---|---|")
    rows = []
    for taut in tauts:
        for lvl in levels:
            for solv in solvents:
                recs = [r for r in results if r["tautomer"] == taut
                        and r["level_id"] == lvl and r["solvent"] == solv]
                b = brightest_of(recs, weighted=True)
                if not b:
                    continue
                e = exp_for[taut]
                d_nm, d_ev, d_pc = err_row(b["lambda_nm"], e)
                solv_lbl = "기체상" if solv == "none" else solv
                A(f"| {TAUT_LABEL[taut]} | {lvl} | {solv_lbl} | {b['lambda_nm']:.1f} | "
                  f"{b['f']:.3f} | {e:.1f} | {d_nm:+.1f} | {d_ev:+.3f} | {d_pc:+.1f} | "
                  f"{b['orbitals']} |")
                rows.append({"tautomer": taut, "level": lvl, "solvent": solv,
                             "lambda_nm": b["lambda_nm"], "f": b["f"],
                             "exp_nm": e, "err_nm": d_nm, "err_eV": d_ev})
    A("")
    df = pd.DataFrame(rows)
    if not df.empty:
        df.to_csv(RESULTS / "lambda_max_summary.csv", index=False,
                  encoding="utf-8-sig")

    # ---------------------------------------------- 3. 용매 효과
    A("## 3. 용매 모델의 영향")
    A("")
    if "none" in solvents and len(solvents) > 1:
        solv2 = [s for s in solvents if s != "none"][0]
        A(f"| 토토머 | 수준 | 기체상 (nm) | {solv2} (nm) | 용매 이동 (nm) | 용매 이동 (eV) |")
        A("|---|---|---|---|---|---|")
        for taut in tauts:
            for lvl in levels:
                g = brightest_of([r for r in results if r["tautomer"] == taut
                                  and r["level_id"] == lvl and r["solvent"] == "none"])
                s = brightest_of([r for r in results if r["tautomer"] == taut
                                  and r["level_id"] == lvl and r["solvent"] == solv2])
                if not g or not s:
                    continue
                d = s["lambda_nm"] - g["lambda_nm"]
                dev = EV_TO_NM / s["lambda_nm"] - EV_TO_NM / g["lambda_nm"]
                A(f"| {TAUT_LABEL[taut]} | {lvl} | {g['lambda_nm']:.1f} | "
                  f"{s['lambda_nm']:.1f} | {d:+.1f} | {dev:+.3f} |")
    else:
        A("_기체상/용매 두 조건이 모두 계산되지 않아 비교 불가._")
    A("")

    # ---------------------------------------------- 4. 컨포머 평균 효과
    A("## 4. 컨포머 평균의 영향")
    A("")
    A("| 토토머 | 수준 | 용매 | 최저 단일 구조 (nm) | 앙상블 평균 (nm) | 차이 (nm) | 컨포머 수 |")
    A("|---|---|---|---|---|---|---|")
    for taut in tauts:
        for lvl in levels:
            for solv in solvents:
                recs = [r for r in results if r["tautomer"] == taut
                        and r["level_id"] == lvl and r["solvent"] == solv]
                if len(recs) < 1:
                    continue
                lo = brightest_of(recs, weighted=False)
                en = brightest_of(recs, weighted=True)
                A(f"| {TAUT_LABEL[taut]} | {lvl} | {solv} | {lo['lambda_nm']:.1f} | "
                  f"{en['lambda_nm']:.1f} | {en['lambda_nm']-lo['lambda_nm']:+.1f} | "
                  f"{en['n_conf']} |")
    A("")

    # ---------------------------------------------- 5. 토토머 선택 효과
    A("## 5. 토토머 선택의 영향")
    A("")
    ens = load(CONFORMERS / "search_summary.json")
    A("| 토토머 | 계산 λmax (nm) | 실험 (nm) | 오차 (nm) |")
    A("|---|---|---|---|")
    best_lvl = levels[0] if levels else None
    best_solv = "ethanol" if "ethanol" in solvents else (solvents[0] if solvents else None)
    for taut in tauts:
        b = brightest_of([r for r in results if r["tautomer"] == taut
                          and r["level_id"] == best_lvl and r["solvent"] == best_solv])
        if not b:
            continue
        d, _, _ = err_row(b["lambda_nm"], exp_for[taut])
        A(f"| {TAUT_LABEL[taut]} | {b['lambda_nm']:.1f} | {exp_for[taut]:.1f} | {d:+.1f} |")
    A("")

    # ---------------------------------------------- 6. 함수/기저셋 효과
    A("## 6. 함수와 기저셋 선택의 영향")
    A("")
    if len(levels) > 1:
        A("| 토토머 | 용매 | " + " | ".join(levels) + " | 최대-최소 차 (nm) |")
        A("|---|---|" + "---|" * (len(levels) + 1))
        for taut in tauts:
            for solv in solvents:
                vals = []
                for lvl in levels:
                    b = brightest_of([r for r in results if r["tautomer"] == taut
                                      and r["level_id"] == lvl and r["solvent"] == solv])
                    vals.append(b["lambda_nm"] if b else None)
                got = [v for v in vals if v is not None]
                if len(got) < 2:
                    continue
                cells = " | ".join(f"{v:.1f}" if v else "-" for v in vals)
                A(f"| {TAUT_LABEL[taut]} | {solv} | {cells} | {max(got)-min(got):.1f} |")
    else:
        A("_이론 수준이 하나뿐이라 비교 불가._")
    A("")

    # ---------------------------------------------- 7. 저비용 기준선 대비
    A("## 7. 저비용 기준선(MINDO/3–TDA 대체) 대비 개선 정도")
    A("")
    mindo = ref["previous_calculations"]["MINDO3_TDA"]
    if mindo.get("user_supplied_lambda_max_nm"):
        base = mindo["user_supplied_lambda_max_nm"]
        A(f"사용자 제공 MINDO/3–TDA 값: **{base} nm**")
    else:
        A("> **주의**: 아보벤존에 대해 발표된 MINDO/3–TDA λmax 를 문헌에서 찾지 못했다.")
        A("> 따라서 본 프로젝트에서 직접 계산한 저비용 기준선을 대신 사용한다.")
        base_rec = load_checkpoint(CALCULATIONS / "00_test" / "baseline_lowcost.json")
        base = base_rec.get("lambda_max_nm") if base_rec else None
        if base:
            A(f"자체 저비용 기준선: **{base:.1f} nm** ({base_rec.get('level','?')})")
    if base:
        b = brightest_of([r for r in results if r["tautomer"] == "enolA"
                          and r["level_id"] == best_lvl and r["solvent"] == best_solv])
        if b:
            e0 = abs(base - exp_enol)
            e1 = abs(b["lambda_nm"] - exp_enol)
            A("")
            A("| 방법 | λmax (nm) | 실험 대비 절대오차 (nm) |")
            A("|---|---|---|")
            A(f"| 저비용 기준선 | {base:.1f} | {e0:.1f} |")
            A(f"| 본 계산 ({best_lvl}, {best_solv}) | {b['lambda_nm']:.1f} | {e1:.1f} |")
            A("")
            if e0 > 0:
                A(f"→ 절대오차가 **{e0:.1f} nm → {e1:.1f} nm** 로 "
                  f"{100*(e0-e1)/e0:+.0f}% 변화했다.")
    A("")

    # ---------------------------------------------- 8. 선행 DFT 연구와 비교
    A("## 8. 선행 DFT/TD-DFT 연구와의 비교")
    A("")
    A("| 연구 | 수준 | 에놀 λmax (nm) | 디케토 λmax (nm) |")
    A("|---|---|---|---|")
    pc = ref["previous_calculations"]
    A(f"| {pc['B3LYP_6-31+Gd_PCM']['ref']} | {pc['B3LYP_6-31+Gd_PCM']['level']} | "
      f"{pc['B3LYP_6-31+Gd_PCM']['enol']['ethanol_nm']} (EtOH) | "
      f"{pc['B3LYP_6-31+Gd_PCM']['keto']['ethanol_nm']} (EtOH) |")
    A(f"| {pc['CAM-B3LYP_TZVP_gas']['ref']} | {pc['CAM-B3LYP_TZVP_gas']['level']} | "
      f"{pc['CAM-B3LYP_TZVP_gas']['enol_S1']['nm']} (기체상) | - |")
    A("| **본 계산** | 위 2장 표 참고 | | |")
    A("")

    # ---------------------------------------------- 9. 실패 기록
    A("## 9. 계산 실패와 대응 기록")
    A("")
    fails = load_checkpoint(LOGS / "failures.json")
    if fails and fails.get("failures"):
        A("| 대상 | 분류 | 대응 |")
        A("|---|---|---|")
        for f in fails["failures"]:
            A(f"| `{f['tag']}` | {f['code']} | {f['remedy'][:120]} |")
    else:
        A("기록된 실패 없음.")
    A("")

    # ---------------------------------------------- 10. 한계
    A("## 10. 남아 있는 계산적 한계")
    A("")
    for line in [
        "**선형응답 PCM 의 한계** — 연속체 모델은 용매를 균일한 유전체로 다룬다. "
        "에탄올·메탄올이 킬레이트 O–H 와 직접 만드는 분자간 수소결합, 그리고 그것이 "
        "토토머 평형과 전이 에너지에 주는 영향은 반영되지 않는다.",
        "**수직 전이 근사** — 진동 구조(Franck–Condon)를 계산하지 않고 가우시안 "
        "broadening 으로 대체했다. 실험 밴드의 비대칭성과 어깨는 재현되지 않는다.",
        "**TD-DFT 자체의 계통 오차** — 이 발색단에 대해 함수마다 결과가 크게 갈린다 "
        "(6장 표). 벤치마크 없이 한 함수만 쓰면 우연히 맞을 수는 있어도 근거가 약하다.",
        "**바닥상태 구조 기반** — 들뜬상태 구조 이완(형광, ESIPT)은 다루지 않았다.",
        "**토토머 존재비 미반영** — 계산은 각 토토머를 따로 다뤘고, 실제 용액의 "
        "에놀:디케토 비율(극성 용매에서 약 85:15~98:2)로 가중한 합성 스펙트럼은 별도 논의가 필요하다.",
    ]:
        A(f"- {line}")
    A("")

    A("## 11. 다음 계산에서 우선적으로 개선할 요소")
    A("")
    A("1. 더 빠른 TD-DFT 엔진(예: RIJCOSX 를 쓰는 ORCA, 또는 GPU 가속)으로 옮겨 "
      "def2-TZVP 급 기저셋과 범위분리 함수를 감당한다.")
    A("2. 명시적 용매 분자 1–2개를 넣은 미세용매화(microsolvation) 모형으로 "
      "수소결합 효과를 확인한다.")
    A("3. 진동 분해(Franck–Condon) 스펙트럼으로 밴드 모양을 재현한다.")
    A("4. 토토머 상대 자유에너지를 계산해 실제 용액 조성으로 가중한 합성 스펙트럼을 만든다.")
    A("")

    out = RESULTS / "report.md"
    out.write_text("\n".join(L), encoding="utf-8")
    print(f"보고서 -> {out}")
    print(f"요약 CSV -> {RESULTS / 'lambda_max_summary.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
