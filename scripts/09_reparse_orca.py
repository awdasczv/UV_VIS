"""
09_reparse_orca.py
------------------
이미 저장된 ORCA 출력(td.out)을 다시 파싱해 result.json 을 갱신한다.
재계산은 하지 않는다 (td.out 이 진실이므로).

왜 필요한가
  파서를 개선했다 (오비탈 인덱스 범위 검사 + 인위적 해 플래그).
  기존 result.json 은 옛 파서로 만든 것이라, 새 검증 필드가 없다.
  이 스크립트가 td.out 을 다시 읽어 transitions/flag/범위검사 결과를 최신화한다.

무엇을 보존하는가
  result.json 의 계산 메타데이터(wall_seconds, boltzmann_weight, geometry_source 등)는
  그대로 두고, 파싱으로 얻는 부분(transitions, brightest, n_basis, flag 등)만 덮어쓴다.

실행:  .\scripts\run.ps1 scripts\09_reparse_orca.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orca_common import EV_TO_NM, parse_output
from qc_common import CALCULATIONS, load_checkpoint, save_checkpoint

PARSED_KEYS = ["scf_energy_hartree", "n_basis", "n_occupied", "n_electrons",
               "transitions", "brightest", "n_out_of_range_orbital_indices"]


def main() -> int:
    roots = [CALCULATIONS / "02_tddft_orca", CALCULATIONS / "03_baseline"]
    n_updated, n_flagged, n_badidx = 0, 0, 0
    flagged_states = []

    for root in roots:
        if not root.exists():
            continue
        for ck in sorted(root.rglob("result.json")):
            td = ck.parent / "td.out"
            if not td.exists():
                td = ck.parent / "cis.out"      # 저비용 기준선
            if not td.exists():
                continue

            rec = load_checkpoint(ck)
            if rec is None:
                continue

            parsed = parse_output(td)
            if not parsed.get("transitions"):
                continue

            # 파싱 부분만 갱신
            for k in PARSED_KEYS:
                if k in parsed:
                    rec[k] = parsed[k]
            # 최단파장/커버 재계산
            lo = min(t["wavelength_nm"] for t in parsed["transitions"])
            rec["shortest_wavelength_nm"] = lo
            rec["covers_200nm"] = lo <= 200.0
            b = max(parsed["transitions"], key=lambda t: t["osc_strength"])
            rec["brightest"] = b

            # 플래그 집계
            flagged = [t for t in parsed["transitions"] if t.get("flag")]
            if flagged:
                n_flagged += len(flagged)
                tag = ck.parent.relative_to(root)
                for t in flagged:
                    flagged_states.append(f"{tag} state {t['state']} "
                                          f"({t['wavelength_nm']}nm f={t['osc_strength']}): "
                                          f"{t['flag']}")
            if parsed.get("n_out_of_range_orbital_indices"):
                n_badidx += parsed["n_out_of_range_orbital_indices"]

            save_checkpoint(ck, rec)
            n_updated += 1

    print(f"재파싱 완료: {n_updated} 개 result.json 갱신")
    print(f"범위 벗어난 오비탈 인덱스: {n_badidx} 개 (0 이면 파서 오독 없음)")
    print(f"인위적 해로 플래그된 전이: {n_flagged} 개")
    for s in flagged_states:
        print(f"  - {s}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
