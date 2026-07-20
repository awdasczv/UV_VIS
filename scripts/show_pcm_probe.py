"""PCM cavity Area 비용/정확도 절충 결과를 표로 출력한다."""
import json
import sys
from pathlib import Path

p = Path(__file__).resolve().parent.parent / "calculations" / "00_test" / "pcm_cost_probe.json"
if not p.exists():
    print("결과 없음")
    sys.exit(0)

d = json.loads(p.read_text(encoding="utf-8"))
GAS_NM, GAS_SCF, GAS_TD = 315.7, 76.2, 149.4   # 03_test_single_point.py 의 T1

print(f"{'조건':16s} {'tessera':>8s} {'SCF s':>8s} {'TDDFT s':>9s} "
      f"{'lambda nm':>10s} {'f':>7s} {'기체상 대비':>12s}")
print("-" * 78)
print(f"{'기체상':16s} {'-':>8s} {GAS_SCF:8.1f} {GAS_TD:9.1f} "
      f"{GAS_NM:10.1f} {1.012:7.3f} {'기준':>12s}")
for r in d["runs"]:
    if not r.get("ok"):
        print(f"{r['tag']:16s} 실패")
        continue
    print(f"{'PCM Area '+str(r['area']):16s} {r['n_tesserae']:8d} "
          f"{r['scf_seconds']:8.1f} {r['tdscf_seconds']:9.1f} "
          f"{r['brightest_nm']:10.1f} {r['brightest_f']:7.3f} "
          f"{r['brightest_nm']-GAS_NM:+11.1f}nm")

ok = [r for r in d["runs"] if r.get("ok")]
if len(ok) >= 2:
    a, b = ok[0], ok[-1]
    print(f"\nArea {a['area']} vs {b['area']} 차이: "
          f"{abs(a['brightest_nm']-b['brightest_nm']):.2f} nm  "
          f"(TD-DFT 시간 {a['tdscf_seconds']:.0f}s vs {b['tdscf_seconds']:.0f}s)")

print("\n--- 전이 상세 (Area 1.0) ---")
for r in ok:
    if r["area"] != 1.0:
        continue
    for t in r["transitions"]:
        bar = "*" * int(min(30, t["osc_strength"] * 30))
        print(f"  S{t['state']:<2d} {t['energy_eV']:6.3f} eV {t['wavelength_nm']:7.1f} nm "
              f"f={t['osc_strength']:.4f} {t['orbital_transitions_str']:28s} {bar}")
