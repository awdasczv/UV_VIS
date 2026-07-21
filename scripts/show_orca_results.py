"""완료된 ORCA TD-DFT 결과를 요약해서 보여준다."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
base = ROOT / "calculations" / "02_tddft_orca"
EXP = {"enolA": 354.9, "enolB": 354.9, "diketo": 265.0}

files = sorted(base.rglob("result.json"))
if not files:
    print("아직 완료된 결과가 없습니다.")
    sys.exit(0)

recs = []
for f in files:
    d = json.loads(f.read_text(encoding="utf-8"))
    if d.get("ok"):
        recs.append(d)
    else:
        print(f"[실패] {f.relative_to(base)}  "
              f"{d.get('failure', {}).get('code', '?')}")

print(f"완료 {len(recs)} 건\n")
print(f"{'토토머':8s} {'컨포머':13s} {'수준':17s} {'용매':8s} "
      f"{'λmax':>8s} {'f':>7s} {'최단nm':>7s} {'실험':>7s} {'오차':>7s} {'초':>6s}")
print("-" * 104)
for r in sorted(recs, key=lambda x: (x["tautomer"], x["solvent"], x["conf_id"])):
    b = r["brightest"]
    exp = EXP.get(r["tautomer"])
    err = b["wavelength_nm"] - exp if exp else 0.0
    print(f"{r['tautomer']:8s} {r['conf_id']:13s} {r['level_id']:17s} "
          f"{r['solvent']:8s} {b['wavelength_nm']:8.1f} {b['osc_strength']:7.3f} "
          f"{r.get('shortest_wavelength_nm', 0):7.1f} {exp:7.1f} {err:+7.1f} "
          f"{r['wall_seconds']:6.0f}")

# 용매 효과
print("\n=== 용매 효과 (기체상 -> 에탄올 CPCM) ===")
for taut in ["enolA", "enolB", "diketo"]:
    for lvl in sorted({r["level_id"] for r in recs}):
        g = [r for r in recs if r["tautomer"] == taut and r["solvent"] == "none"
             and r["level_id"] == lvl]
        s = [r for r in recs if r["tautomer"] == taut and r["solvent"] == "ethanol"
             and r["level_id"] == lvl]
        if g and s:
            gn = g[0]["brightest"]["wavelength_nm"]
            sn = s[0]["brightest"]["wavelength_nm"]
            gf = g[0]["brightest"]["osc_strength"]
            sf = s[0]["brightest"]["osc_strength"]
            print(f"  {taut:8s} {lvl:17s} {gn:7.1f} -> {sn:7.1f} nm "
                  f"({sn-gn:+6.1f} nm)   f {gf:.3f} -> {sf:.3f}")

print("\n=== 주요 전이 성격 ===")
for r in sorted(recs, key=lambda x: (x["tautomer"], x["solvent"])):
    b = r["brightest"]
    print(f"  {r['tautomer']:8s} {r['solvent']:8s} {b['wavelength_nm']:7.1f} nm  "
          f"{b['orbital_transitions_str']}")
