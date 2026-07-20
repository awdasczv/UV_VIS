"""완료된 TD-DFT 결과를 사람이 읽기 좋게 출력한다."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
base = ROOT / "calculations" / "02_tddft"

files = sorted(base.rglob("result.json"))
if not files:
    print("아직 완료된 TD-DFT 결과가 없습니다.")
    sys.exit(0)

for f in files:
    d = json.loads(f.read_text(encoding="utf-8"))
    if not d.get("ok"):
        print(f"\n[실패] {f.relative_to(base)}  "
              f"{d.get('failure', {}).get('code', '?')}")
        continue
    print(f"\n=== {d['tautomer']} / {d['conf_id']} / "
          f"{d['functional']}/{d['basis']} / 용매={d['solvent']} ===")
    print(f"    기저함수 {d['n_basis']}, 점유궤도 {d['n_occupied']}, "
          f"SCF {d['scf_seconds']:.0f}s + TD-DFT {d['tdscf_seconds']:.0f}s, "
          f"구조출처 {d['geometry_source']}")
    trans = d["transitions"]
    for t in trans:
        bar = "*" * int(min(34, t["osc_strength"] * 34))
        print(f"    S{t['state']:<2d} {t['energy_eV']:6.3f} eV {t['wavelength_nm']:7.1f} nm "
              f"f={t['osc_strength']:.4f}  {t['orbital_transitions_str'][:46]:46s} {bar}")
    b = d["brightest"]
    lo = min(t["wavelength_nm"] for t in trans)
    print(f"    -> 최강 전이 {b['wavelength_nm']:.1f} nm (f={b['osc_strength']:.3f})")
    print(f"    -> 계산된 최단파장 {lo:.1f} nm "
          f"({'200 nm 커버 OK' if lo <= 200 else '200 nm 까지 못 미침 - 상태 수 부족'})")
