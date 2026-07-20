"""테스트 계산 비용 요약을 보기 좋게 출력한다."""
import json
import sys
from pathlib import Path

p = Path(__file__).resolve().parent.parent / "calculations" / "00_test" / "test_report.json"
if not p.exists():
    print("아직 결과 없음")
    sys.exit(0)

d = json.loads(p.read_text(encoding="utf-8"))
print(f"{'case':16s} {'nbf':>5s} {'SCF s':>8s} {'TDDFT s':>9s} {'합계 s':>8s} "
      f"{'lmax nm':>9s} {'f':>7s}")
print("-" * 70)
for c in d["cases"]:
    if c.get("ok"):
        tot = c["scf_seconds"] + c["tdscf_seconds"]
        b = c["brightest"]
        print(f"{c['name']:16s} {c['nbf']:5d} {c['scf_seconds']:8.1f} "
              f"{c['tdscf_seconds']:9.1f} {tot:8.1f} "
              f"{b['wavelength_nm']:9.1f} {b['osc_strength']:7.3f}")
    else:
        err = c.get("error", "").strip().splitlines()
        print(f"{c['name']:16s} {'FAILED':>25s}  {err[-1][:60] if err else ''}")

print()
for c in d["cases"]:
    if not c.get("ok"):
        continue
    print(f"[{c['name']}] {c['functional']}/{c['basis']} PCM={c['pcm']}")
    for t in c["transitions"]:
        bar = "*" * int(min(30, t["osc_strength"] * 30))
        print(f"   S{t['state']:<2d} {t['energy_eV']:6.3f} eV {t['wavelength_nm']:7.1f} nm "
              f"f={t['osc_strength']:.4f} {bar}")
