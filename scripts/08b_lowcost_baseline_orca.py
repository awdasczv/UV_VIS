"""
08b_lowcost_baseline_orca.py
----------------------------
"MINDO/3-TDA 대체" 저비용 기준선을 ORCA 로 계산 (요구사항 15번의 비교 대상).

왜 필요한가
  아보벤존에 대해 발표된 MINDO/3-TDA lambda_max 를 문헌에서 못 찾았다
  (Google Scholar 전문검색 0건, docs/DEVELOPMENT_LOG.md 6장).
  "저비용 계산보다 얼마나 좋아졌는가"를 말하려면 기준선을 직접 만들어야 한다.

무엇을 기준선으로 삼는가
  MINDO/3-TDA = 반경험적 해밀토니안 + TDA(=CIS) 형태의 들뜬상태.
  이에 가장 가까우면서 우리 엔진으로 재현 가능한 것은 HF/최소기저 + CIS 다.
  전자상관도 없고 기저셋도 최소라 반경험적 수준의 오차 크기를 보여준다.
    baseline_min : HF/STO-3G + CIS   (가장 저렴, MINDO/3 성격에 가장 근접)
    baseline_hf  : HF/def2-SVP + CIS (기저셋만 키운 대조군)

  ※ 이것은 MINDO/3 자체가 아니다. 보고서에는 반드시
    "문헌에 MINDO/3 값이 없어 동급의 저비용 기준선을 자체 계산했다" 고 명시한다.

  비교 공정성을 위해 최고 수준과 같은 GFN2-xTB 구조(enolA_c000)를 쓴다.
  (최고 수준은 DFT 구조를 쓰지만, 저비용 방법의 '한계'를 보이는 것이 목적이므로
   구조까지 저비용 쪽 조건에 맞추지 않고, 오히려 같은 출발 구조에서 방법만
   바꿨을 때의 차이를 본다. xTB 구조는 두 경우 모두에 존재한다.)

실행:  .\scripts\run.ps1 scripts\08b_lowcost_baseline_orca.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orca_common import build_input, parse_output, run_orca
from qc_common import CALCULATIONS, CONFORMERS, INPUTS, read_xyz, save_checkpoint

OUT = CALCULATIONS / "03_baseline"
NSTATES = 12

# (id, 방법 키워드, 기저셋, 설명)
BASELINES = [
    ("hf_sto3g_cis", "HF", "STO-3G",
     "HF/STO-3G + CIS. 반경험적 MINDO/3-TDA 와 성격이 가장 가까운 저비용 조합."),
    ("hf_def2svp_cis", "HF", "def2-SVP",
     "HF/def2-SVP + CIS. 기저셋만 키운 저비용 대조군."),
]


def main() -> int:
    ref = json.loads((INPUTS / "experimental_reference.json").read_text(encoding="utf-8"))
    exp = ref["experimental"]["enol_band"]["primary_target"]["lambda_max_nm"]

    xyz = CONFORMERS / "enolA" / "enolA_c000.xyz"
    if not xyz.exists():
        print("먼저 컨포머를 준비하세요 (02d_extract_conformers.py).")
        return 1
    g = read_xyz(xyz)

    report = {
        "purpose": "MINDO/3-TDA 문헌값이 없어 동급 저비용 기준선을 자체 계산",
        "caveat": "이것은 MINDO/3 자체가 아니다. 반경험적 수준의 오차 크기를 보여주는 대용물이다.",
        "structure": "enolA_c000 (GFN2-xTB 구조, 최고 수준과 동일 출발 구조)",
        "experimental_nm": exp, "runs": [],
    }

    for bid, method, basis, desc in BASELINES:
        print(f"\n--- {bid}: {method}/{basis} + CIS, {NSTATES}상태, 기체상 ---")
        outdir = OUT / bid
        # ORCA 에서 HF + CIS: functional 자리에 HF, %tddft 대신 %cis 를 쓴다.
        # build_input 은 %tddft 를 쓰므로 HF 에도 통하도록 method=HF 로 넣고
        # TDA=true 로 두면 ORCA 가 CIS 로 처리한다.
        inp = build_input(
            g.symbols, g.coords, functional=method, basis=basis,
            nstates=NSTATES, tda=True, solvent=None,
            nprocs=4, maxcore_mb=2500, rijcosx=False,  # HF 소형이라 RIJCOSX 불필요
            comment=f"low-cost baseline: {bid}",
        )
        out_path, seconds = run_orca(inp, outdir, name="cis")
        parsed = parse_output(out_path)
        rec = {"id": bid, "method": method, "basis": basis, "description": desc,
               "wall_seconds": round(seconds, 1)}
        if parsed.get("terminated_normally") and parsed.get("transitions"):
            b = parsed["brightest"]
            rec.update({"ok": True, "n_basis": parsed["n_basis"],
                        "lambda_max_nm": b["wavelength_nm"],
                        "osc_strength": b["osc_strength"],
                        "level": f"{method}/{basis} + CIS, 기체상",
                        "transitions": parsed["transitions"]})
            print(f"    기저함수 {parsed['n_basis']}  ({seconds:.0f}초)")
            print(f"    lambda_max = {b['wavelength_nm']:.1f} nm (f={b['osc_strength']:.3f}) "
                  f"-> 실험 대비 {b['wavelength_nm']-exp:+.1f} nm")
        else:
            rec["ok"] = False
            print("    실패")
        report["runs"].append(rec)
        save_checkpoint(OUT / "baseline.json", report)

    ok = [r for r in report["runs"] if r.get("ok")]
    if ok:
        report["lambda_max_nm"] = ok[0]["lambda_max_nm"]
        report["level"] = ok[0]["level"]
        save_checkpoint(OUT / "baseline.json", report)
        print("\n=== 저비용 기준선 요약 (에놀 A, 실험 {:.1f} nm) ===".format(exp))
        for r in ok:
            print(f"  {r['id']:16s} {r['n_basis']:4d} 기저함수  "
                  f"{r['lambda_max_nm']:7.1f} nm  (실험 대비 {r['lambda_max_nm']-exp:+.1f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
