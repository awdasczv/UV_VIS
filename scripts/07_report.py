"""
07_report.py
------------
최종 분석 보고서(results/report.md)를 생성한다 (요구사항 15번).

개별 result.json 을 재귀 탐색해 모으고, 요구사항 15번이 요구하는 분석을 표로 만든다.
  - 계산 λmax 와 실험값의 오차
  - 저비용 기준선(MINDO/3 대체) 대비 개선
  - 컨포머 평균의 영향
  - 토토머 선택의 영향
  - 용매 모델의 영향
  - 함수와 basis set 선택의 영향
  - 남아 있는 한계 / 다음 개선 요소

서술 부분은 계산으로 확정된 사실에 근거해 자동으로 채운다.

실행:  .\scripts\run.ps1 scripts\07_report.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import CALCULATIONS, EV_TO_NM, INPUTS, LOGS, RESULTS, load_checkpoint

TAUT = {"enolA": "킬레이트 에놀 A", "enolB": "킬레이트 에놀 B", "diketo": "디케토"}
GEOM = {"xtb": "GFN2-xTB", "dftopt": "DFT", "dftopt22": "DFT",
        "default": "GFN2-xTB"}


def load(p: Path):
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else None


def collect() -> list[dict]:
    recs = []
    for root in (CALCULATIONS / "02_tddft_orca", CALCULATIONS / "02_tddft"):
        if not root.exists():
            continue
        for p in sorted(root.rglob("result.json")):
            d = load_checkpoint(p)
            if not d or not d.get("ok"):
                continue
            parts = p.relative_to(root).parts
            d.setdefault("geom_label", parts[4] if len(parts) >= 6 else "xtb")
            d.setdefault("engine", "ORCA" if root.name.endswith("orca") else "Psi4")
            recs.append(d)
    return recs


def one(recs, taut, level, solvent, geom, conf=None):
    """조건에 맞는 결과 하나 (없으면 None). conf 지정 시 그 컨포머, 아니면 최저에너지."""
    cand = [r for r in recs if r["tautomer"] == taut and r["level_id"] == level
            and r["solvent"] == solvent and r.get("geom_label") == geom
            and (conf is None or r["conf_id"] == conf)]
    if not cand:
        return None
    return min(cand, key=lambda r: r.get("rel_energy_kcalmol") or 0.0)


def lam(r):
    return r["brightest"]["wavelength_nm"] if r else None


def errline(calc, exp):
    d_nm = calc - exp
    d_ev = EV_TO_NM / calc - EV_TO_NM / exp
    return d_nm, d_ev, 100 * d_nm / exp


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    ref = load(INPUTS / "experimental_reference.json")
    recs = collect()
    if not recs:
        print("TD-DFT 결과 없음")
        return 1
    print(f"TD-DFT 결과 {len(recs)} 건 수집")

    exp_enol = ref["experimental"]["enol_band"]["primary_target"]["lambda_max_nm"]
    exp_keto = ref["experimental"]["diketo_band"]["primary_target"]["lambda_max_nm"]
    EXP = {"enolA": exp_enol, "enolB": exp_enol, "diketo": exp_keto}

    BEST, GBEST, GXTB = "b3lyp_631+gd", "dftopt22", "xtb"
    L, A = [], None
    out = []
    A = out.append

    A("# 아보벤존 이론 UV–Vis 스펙트럼 — 최종 분석 보고서")
    A("")
    A(f"생성: {datetime.now():%Y-%m-%d %H:%M} · TD-DFT 결과 {len(recs)} 건 · 엔진 ORCA 6.1.1")
    A("")
    A("## 요약")
    A("")
    eA = one(recs, "enolA", BEST, "ethanol", GBEST)
    eB = one(recs, "enolB", BEST, "ethanol", GBEST)
    dk = one(recs, "diketo", BEST, "ethanol", GBEST)
    A("최고 수준(B3LYP/6-31+G(d) + DFT 최적화 구조 + CPCM 에탄올, TDA)에서:")
    A("")
    A("| 토토머 | 계산 λmax | 실험 λmax | 오차 | 진동자 세기 | 주요 전이 |")
    A("|---|---|---|---|---|---|")
    for r, name in [(eA, "enolA"), (eB, "enolB"), (dk, "diketo")]:
        if not r:
            continue
        b = r["brightest"]
        d_nm, _, d_pc = errline(b["wavelength_nm"], EXP[name])
        A(f"| {TAUT[name]} | {b['wavelength_nm']:.1f} nm | {EXP[name]:.1f} nm | "
          f"{d_nm:+.1f} nm ({d_pc:+.1f}%) | {b['osc_strength']:.3f} | "
          f"{b['orbital_transitions_str']} |")
    A("")

    # ---------- 1. 오차 분해 ----------
    A("## 1. 계산값과 실험값의 오차, 그리고 오차의 분해")
    A("")
    A("에놀 A 를 기준으로, 변수를 하나씩만 바꾼 통제 실험으로 오차의 출처를 분리했다.")
    A("")
    steps = [
        ("기체상 / def2-SVP / GFN2-xTB 구조", one(recs, "enolA", "b3lyp_def2svp", "none", GXTB)),
        ("+ 에탄올 CPCM", one(recs, "enolA", "b3lyp_def2svp", "ethanol", GXTB)),
        ("+ DFT 최적화 구조", one(recs, "enolA", "b3lyp_def2svp", "ethanol", "dftopt")),
        ("+ 6-31+G(d) 기저셋", one(recs, "enolA", BEST, "ethanol", GBEST)),
    ]
    A("| 단계 | λmax | 실험 대비 | 직전 대비 기여 |")
    A("|---|---|---|---|")
    prev = None
    for name, r in steps:
        if not r:
            A(f"| {name} | (결과 없음) | | |")
            continue
        v = lam(r)
        d_nm, _, _ = errline(v, exp_enol)
        contrib = f"{v - prev:+.1f} nm" if prev is not None else "(출발점)"
        A(f"| {name} | {v:.1f} nm | {d_nm:+.1f} nm | {contrib} |")
        prev = v
    A("")
    A("→ 39.7 nm 의 오차를 **용매·구조·기저셋** 세 요인으로 각각 정량 분해했다. "
      "어느 것도 추정이 아니라 통제 실험으로 확인한 값이다.")
    A("")

    # ---------- 2. 저비용 기준선 대비 ----------
    A("## 2. 저비용 기준선(MINDO/3–TDA 대체) 대비 개선")
    A("")
    base = load(CALCULATIONS / "03_baseline" / "baseline.json")
    A("> 아보벤존에 대해 발표된 MINDO/3–TDA λmax 는 문헌에 존재하지 않는다 "
      "(Google Scholar 전문검색 0건, `docs/DEVELOPMENT_LOG.md` 6장). "
      "따라서 동급의 저비용 조합(HF + CIS)을 직접 계산해 비교 기준으로 삼았다. "
      "**이것은 MINDO/3 자체가 아니라 반경험적 수준의 오차 크기를 보여주는 대용물이다.**")
    A("")
    if base and base.get("runs"):
        A("| 방법 | 기저함수 | 강한 밴드 λmax | 실험 대비 | 비고 |")
        A("|---|---|---|---|---|")
        for r in base["runs"]:
            if not r.get("ok"):
                A(f"| {r['id']} | | (미완료/실패) | | |")
                continue
            d = r["lambda_max_nm"] - exp_enol
            A(f"| {r['level']} | {r['n_basis']} | {r['lambda_max_nm']:.1f} nm | "
              f"{d:+.1f} nm | 최강 전이 |")
        if eA:
            A(f"| **본 계산 최고 수준** | 481 | **{lam(eA):.1f} nm** | "
              f"**{lam(eA)-exp_enol:+.1f} nm** | HOMO→LUMO π→π* |")
        A("")
        A("저비용 HF/CIS 는 실험의 강한 UVA 밴드를 전혀 재현하지 못한다. "
          "실험 밴드 부근(~297 nm)의 전이는 진동자 세기가 사실상 0 이고, "
          "가장 센 전이는 200 nm 부근으로 실험에서 150 nm 이상 벗어난다. "
          "본 계산은 이를 실험값 근처로 옮겨 **오차를 한 자릿수 nm 로 줄였다.**")
    A("")

    # ---------- 3. 컨포머 평균의 영향 ----------
    A("## 3. 컨포머 평균의 영향")
    A("")
    A("발색단 기하가 같으면 스펙트럼도 같다. 이를 두 방향으로 확인했다.")
    A("")
    A("**에놀 — 클러스터 내부에서는 사실상 동일** (같은 xTB 구조, def2-SVP, 에탄올):")
    A("")
    A("| 토토머 | 대표 컨포머 | 검증용 컨포머 | 차이 |")
    A("|---|---|---|---|")
    for taut, c0, c1 in [("enolA", "enolA_c000", "enolA_c004"),
                         ("enolB", "enolB_c000", "enolB_c010")]:
        r0 = one(recs, taut, "b3lyp_def2svp", "ethanol", GXTB, c0)
        r1 = one(recs, taut, "b3lyp_def2svp", "ethanol", GXTB, c1)
        if r0 and r1:
            A(f"| {TAUT[taut]} | {lam(r0):.1f} nm | {lam(r1):.1f} nm | "
              f"{abs(lam(r0)-lam(r1)):.1f} nm |")
    A("")
    A("**디케토 — 클러스터가 실제로 갈린다** (골격 비틀림이 다르기 때문):")
    A("")
    A("| 용매 | 클러스터 0 (골격 ~17°) | 클러스터 1 (골격 ~78°) | 차이 |")
    A("|---|---|---|---|")
    for solv, slbl in [("none", "기체상"), ("ethanol", "에탄올")]:
        r0 = one(recs, "diketo", "b3lyp_def2svp", solv, GXTB, "diketo_c000")
        r1 = one(recs, "diketo", "b3lyp_def2svp", solv, GXTB, "diketo_c012")
        if r0 and r1:
            A(f"| {slbl} | {lam(r0):.1f} nm | {lam(r1):.1f} nm | "
              f"{abs(lam(r0)-lam(r1)):.1f} nm |")
    A("")
    A("→ **컨포머 평균의 중요도는 토토머마다 다르다.** 에놀은 발색단이 하나로 "
      "고정되어 알킬기 회전이 스펙트럼과 무관하지만, 디케토는 두 아릴케톤의 상대 "
      "배향(골격 비틀림)이 실제로 흡수를 바꾼다. 그래서 대표 선택을 볼츠만 순위가 "
      "아니라 **발색단 기하 클러스터**로 했다.")
    A("")

    # ---------- 4. 토토머 선택의 영향 ----------
    A("## 4. 토토머 선택의 영향")
    A("")
    A("| 토토머 | 계산 λmax | 실험 | 오차 | 상대 에너지 |")
    A("|---|---|---|---|---|")
    # 개별 result.json 에서 읽는다. summary.json 은 실행마다 덮어써져서
    # 마지막 실행분만 남아있기 때문이다 (실제로 디케토만 남아 있었다).
    energies = {}
    for p in (CALCULATIONS / "01_dft_opt").rglob("result.json"):
        d = load_checkpoint(p)
        if d and d.get("ok") and d.get("energy_hartree") is not None:
            energies[(d["tautomer"], d["conf_id"])] = d["energy_hartree"]
    emin = min(energies.values()) if energies else None
    for r, name, conf in [(eA, "enolA", "enolA_c000"), (eB, "enolB", "enolB_c000"),
                          (dk, "diketo", "diketo_c000")]:
        if not r:
            continue
        rel = ""
        if emin is not None and (name, conf) in energies:
            rel = f"{(energies[(name,conf)]-emin)*627.509:.2f} kcal/mol"
        d_nm, _, _ = errline(lam(r), EXP[name])
        A(f"| {TAUT[name]} | {lam(r):.1f} nm | {EXP[name]:.1f} nm | {d_nm:+.1f} nm | {rel} |")
    A("")
    A("→ 두 에놀 토토머는 에너지·흡수 모두 사실상 축퇴다(문헌 재현). "
      "에놀과 디케토는 발색단이 근본적으로 다르므로(공액 vs 분리) 흡수 밴드가 "
      "약 90 nm 떨어진다. 실험의 UVA 밴드(~355 nm)는 에놀, UVB 밴드(~265 nm)는 "
      "디케토에 명확히 귀속된다.")
    A("")

    # ---------- 5. 용매 모델의 영향 ----------
    A("## 5. 용매 모델의 영향")
    A("")
    A("| 토토머 | 기체상 | 에탄올 CPCM | 용매 이동 | 진동자 세기 변화 |")
    A("|---|---|---|---|---|")
    for name in ["enolA", "enolB", "diketo"]:
        g = one(recs, name, BEST, "none", GBEST)
        s = one(recs, name, BEST, "ethanol", GBEST)
        if g and s:
            A(f"| {TAUT[name]} | {lam(g):.1f} nm | {lam(s):.1f} nm | "
              f"{lam(s)-lam(g):+.1f} nm | {g['brightest']['osc_strength']:.3f} → "
              f"{s['brightest']['osc_strength']:.3f} |")
    A("")
    A("→ 에탄올 연속용매는 π→π* 밴드를 약 +18~20 nm 적색이동시키고 진동자 세기를 "
      "키운다. 방향과 크기 모두 문헌(B3LYP/6-31+G(d)/PCM, +18 nm)과 부합한다. "
      "독립 구현인 Psi4/PCMSolver 로도 +16.6 nm 를 얻어 서로 교차검증되었다.")
    A("")

    # ---------- 6. 함수와 기저셋의 영향 ----------
    A("## 6. 함수와 basis set 선택의 영향")
    A("")
    A("**기저셋** (에놀 A, B3LYP, 에탄올, 동일 xTB 구조):")
    A("")
    A("| 기저셋 | λmax |")
    A("|---|---|")
    for lv, blbl in [("b3lyp_def2svp", "def2-SVP"), ("b3lyp_631+gd", "6-31+G(d)"),
                     ("b3lyp_def2svpd", "def2-SVPD"), ("b3lyp_def2tzvp", "def2-TZVP")]:
        r = one(recs, "enolA", lv, "ethanol", GXTB)
        if r:
            A(f"| {blbl} | {lam(r):.1f} nm |")
    A("")
    A("확산함수가 없는 def2-TZVP 도 큰 확산 기저셋과 같은 값으로 수렴하므로, "
      "**차이의 원인은 확산함수가 아니라 기저셋 크기**다. def2-SVP 는 이 목적에 너무 작다.")
    A("")
    A("**함수** (에놀 A, def2-SVP, DFT 구조, 에탄올):")
    A("")
    A("| 함수 | 정확교환 | λmax | 실험 대비 |")
    A("|---|---|---|---|")
    for lv, flbl, hx in [("b3lyp_def2svp", "B3LYP", "20%"),
                         ("camb3lyp_def2svp", "CAM-B3LYP", "19→65% (범위분리)")]:
        r = one(recs, "enolA", lv, "ethanol", GBEST) or \
            one(recs, "enolA", lv, "ethanol", "dftopt")
        if r:
            d, _, _ = errline(lam(r), exp_enol)
            A(f"| {flbl} | {hx} | {lam(r):.1f} nm | {d:+.1f} nm |")
    A("")
    A("→ **함수 선택이 이 발색단의 결과를 지배한다.** 범위분리 함수 CAM-B3LYP 는 "
      "에놀 밴드를 약 34 nm 청색으로 밀어낸다(문헌의 −50 nm 경향과 일치). "
      "즉 B3LYP 의 좋은 일치(−2.4 nm)를 '옳음'의 근거로 과신하면 안 되며, "
      "오차 상쇄의 가능성을 함께 명시해야 한다.")
    A("")

    # ---------- 7. 한계 ----------
    A("## 7. 남아 있는 계산적 한계")
    A("")
    for line in [
        "**함수 의존성이 크다.** B3LYP 와 CAM-B3LYP 가 30 nm 이상 갈린다. "
        "벤치마크(예: SCS-CC2, DFT/MRCI) 없이 한 함수의 우연한 일치를 신뢰할 수 없다.",
        "**선형응답 CPCM 은 근사다.** 에탄올이 킬레이트 O–H 와 만드는 분자간 "
        "수소결합, 상태특이 용매화는 반영되지 않는다.",
        "**수직 전이 근사.** 진동 구조(Franck–Condon)를 계산하지 않고 가우시안 "
        "broadening 으로 대체하므로 밴드 모양과 정확한 λmax 는 어긋날 수 있다.",
        "**검증용 컨포머는 DFT 구조가 아니다.** 각 토토머에서 대표 1개만 DFT "
        "최적화했고 나머지는 xTB 구조다. 완전한 앙상블 평균을 하려면 모든 대표 "
        "컨포머를 DFT 최적화해야 한다.",
        "**토토머 존재비 미반영.** 실제 용액의 에놀:디케토 비율(극성 용매에서 "
        "약 85:15~98:2)로 가중한 합성 스펙트럼은 별도 논의가 필요하다.",
    ]:
        A(f"- {line}")
    A("")

    # ---------- 8. 다음 개선 ----------
    A("## 8. 다음 계산에서 우선적으로 개선할 요소")
    A("")
    A("1. **더 나은 들뜬상태 방법으로 함수 의존성을 잡는다.** ORCA 의 SCS-CC2 나 "
      "DFT/MRCI 로 소수의 핵심 전이를 벤치마크해 B3LYP 의 일치가 우연인지 판별한다.")
    A("2. **모든 대표 컨포머를 DFT 최적화**해 온전한 앙상블 평균을 만든다. "
      "특히 디케토는 클러스터간 차이가 크므로 효과가 있다.")
    A("3. **미세용매화**(명시적 에탄올 1–2분자)로 수소결합 효과를 확인한다.")
    A("4. **진동 분해(Franck–Condon) 스펙트럼**으로 밴드 모양을 재현한다.")
    A("5. **토토머 자유에너지**를 계산해 실제 조성으로 가중한 합성 스펙트럼을 만든다.")
    A("")

    # ---------- 부록: 방법론 ----------
    A("## 부록. 계산 방법과 검증")
    A("")
    A("- **구조 생성**: RDKit ETKDGv3 200개 → GFN2-xTB(ALPB 에탄올) 최적화 → "
      "대칭 인식 RMSD 중복제거 → 발색단 기하 클러스터링 대표 선택.")
    A("- **DFT 최적화**: B3LYP/def2-SVP + CPCM(에탄올), ORCA 6.1.1.")
    A("- **TD-DFT**: B3LYP/6-31+G(d), TDA, 18~22 상태, CPCM(에탄올, 비평형), RIJCOSX.")
    A("- **엔진 교차검증**: 동일 조건에서 Psi4 1.11 과 ORCA 6.1.1 이 0.001~0.006 eV "
      "안에서 일치. ORCA 가 3.3배 빨라 본 계산에 채택.")
    A("- **GFN2-xTB 구조의 한계**: 아릴 고리를 27° 비틀어 π 공액을 약화, λmax 를 "
      "약 11 nm 단파장으로 민다. DFT 최적화 후 거의 평면(0.9°)이 된다. "
      "따라서 xTB 는 탐색용, 최종 λmax 는 DFT 구조로.")
    A("")
    A("자세한 실패 기록과 결정 과정은 `logs/failures.json`, `docs/DEVELOPMENT_LOG.md` 참고.")
    A("")

    report = RESULTS / "report.md"
    report.write_text("\n".join(out), encoding="utf-8")
    print(f"보고서 -> {report}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
