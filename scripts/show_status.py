"""
show_status.py
--------------
프로젝트 전체 진행 상황을 한눈에 보여준다.
어떤 계산이 끝났고, 어떤 결과가 나왔고, 무엇이 남았는지.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CALC = ROOT / "calculations"
CONF = ROOT / "conformers"
EXP = {"enolA": 354.9, "enolB": 354.9, "diketo": 265.0}


def load(p: Path):
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:                                    # noqa: BLE001
        return None


print("=" * 78)
print("1. 컨포머 탐색")
print("=" * 78)
s = load(CONF / "selection_summary.json")
if s:
    for taut, v in s.items():
        print(f"  {taut:8s} 발색단 클러스터 {v['n_clusters']}종, "
              f"커버 가중치 {v['cum_weight']:.3f}, 대표 {v['selected']}")

print()
print("=" * 78)
print("2. DFT 구조 최적화")
print("=" * 78)
opts = sorted((CALC / "01_dft_opt").rglob("result.json"))
if not opts:
    print("  없음")
for p in opts:
    d = load(p)
    if not d or not d.get("ok"):
        print(f"  {p.parent.name}: 실패 또는 진행중")
        continue
    print(f"  {d['tautomer']:8s}/{d['conf_id']:14s} E={d['energy_hartree']:.8f} Eh  "
          f"{d['n_opt_cycles']:2d} 사이클  {d['wall_seconds']/60:.0f}분")
    b, a = d.get("twists_before"), d.get("twists_after")
    if b and a:
        for k in sorted(set(b) & set(a)):
            unit = "A" if k.endswith("_A") else "도"
            print(f"       {k:26s} {b[k]:7.2f} -> {a[k]:7.2f} {unit} "
                  f"({a[k]-b[k]:+.2f})")
    else:
        print("       (발색단 비틀림각 미측정 - 디케토는 에놀 패턴과 달라 별도 정의 필요)")

print()
print("=" * 78)
print("3. TD-DFT 결과")
print("=" * 78)
recs = []
root = CALC / "02_tddft_orca"
if root.exists():
    for p in sorted(root.rglob("result.json")):
        d = load(p)
        if d and d.get("ok"):
            parts = p.relative_to(root).parts
            d["geom_label"] = parts[4] if len(parts) >= 6 else "xtb"
            recs.append(d)

if not recs:
    print("  없음")
else:
    print(f"{'토토머':8s} {'컨포머':13s} {'수준':15s} {'용매':8s} {'구조':9s} "
          f"{'λmax':>7s} {'f':>6s} {'최단':>6s} {'오차':>7s}")
    print("-" * 92)
    for r in sorted(recs, key=lambda x: (x["tautomer"], x["geom_label"],
                                         x["level_id"], x["solvent"])):
        b = r["brightest"]
        exp = EXP.get(r["tautomer"], 0)
        gl = {"xtb": "xTB"}.get(r["geom_label"], r["geom_label"])
        print(f"{r['tautomer']:8s} {r['conf_id']:13s} {r['level_id']:15s} "
              f"{r['solvent']:8s} {gl:9s} {b['wavelength_nm']:7.1f} "
              f"{b['osc_strength']:6.3f} {r.get('shortest_wavelength_nm',0):6.1f} "
              f"{b['wavelength_nm']-exp:+7.1f}")

    # 최고 수준 결과만 추림
    best = [r for r in recs if r["geom_label"].startswith("dftopt")
            and r["solvent"] == "ethanol"]
    if best:
        print()
        print("  --- DFT 구조 + 에탄올 (최종 비교 대상) ---")
        for r in sorted(best, key=lambda x: x["tautomer"]):
            b = r["brightest"]
            exp = EXP.get(r["tautomer"], 0)
            print(f"    {r['tautomer']:8s} {r['level_id']:15s} "
                  f"{b['wavelength_nm']:7.1f} nm (실험 {exp:.1f}, "
                  f"오차 {b['wavelength_nm']-exp:+.1f} nm)  {b['orbital_transitions_str'][:40]}")

print()
print("=" * 78)
print("4. 남은 작업")
print("=" * 78)
todo = []
if not (CALC / "00_test" / "baseline_lowcost.json").exists():
    todo.append("저비용 기준선 (MINDO/3 대체)")
if not any(r["level_id"].startswith("camb3lyp") for r in recs):
    todo.append("CAM-B3LYP 함수 비교")
if not (ROOT / "results" / "spectra_all.csv").exists():
    todo.append("스펙트럼 생성 + CSV + 그래프")
if not (ROOT / "results" / "report.md").exists():
    todo.append("최종 분석 보고서")
for t in todo:
    print(f"  - {t}")
if not todo:
    print("  모두 완료")
