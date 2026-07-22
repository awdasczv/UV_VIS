"""저비용 기준선 결과 요약."""
import json
import sys
from pathlib import Path

p = Path(__file__).resolve().parent.parent / "calculations" / "03_baseline" / "baseline.json"
d = json.loads(p.read_text(encoding="utf-8"))
exp = d["experimental_nm"]
print(f"실험(에놀, 에탄올): {exp} nm\n")
for r in d["runs"]:
    if not r.get("ok"):
        print(f"{r['id']}: 실패/미완료")
        continue
    b = r["lambda_max_nm"]
    print(f"{r['level']:24s} 기저함수 {r['n_basis']:4d}  "
          f"최강밴드 {b:6.1f} nm (f={r['osc_strength']:.3f})  실험대비 {b-exp:+.1f} nm  "
          f"[{r['wall_seconds']:.0f}초]")
    # 실험 밴드(354.9 nm) 부근 전이가 어두운지 확인
    near = [t for t in r["transitions"] if 320 <= t["wavelength_nm"] <= 380]
    if near:
        for t in near:
            print(f"    실험밴드 부근: {t['wavelength_nm']:.1f} nm  f={t['osc_strength']:.4f}  "
                  f"{'(어두움)' if t['osc_strength'] < 0.05 else '(밝음)'}")
