"""
02d_extract_conformers.py
-------------------------
conformers/<토토머>/conformers.xyz 에서 개별 컨포머 xyz 파일을 복원한다.

왜 필요한가
  개별 컨포머 xyz 는 git 으로 추적하지 않는다. 같은 좌표가 앙상블 파일
  conformers.xyz 안에 전부 들어 있어 중복이기 때문이다 (파일 46개 절약).
  그런데 이후 단계(04 최적화, 05 TD-DFT)는 개별 파일 경로를 입력으로 받는다.
  저장소를 새로 clone 한 뒤에는 이 스크립트를 한 번 돌려서 복원하면 된다.

실행:  .\scripts\run.ps1 scripts\02d_extract_conformers.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import CONFORMERS, read_multi_xyz

TAUTOMERS = ["enolA", "enolB", "diketo"]


def main() -> int:
    total = 0
    for taut in TAUTOMERS:
        ens = CONFORMERS / taut / "conformers.xyz"
        if not ens.exists():
            print(f"[건너뜀] {taut}: {ens} 없음. 먼저 02_conformer_search.py 실행.")
            continue
        geoms = read_multi_xyz(ens)
        made, skipped = 0, 0
        for i, g in enumerate(geoms):
            out = CONFORMERS / taut / f"{taut}_c{i:03d}.xyz"
            if out.exists():
                skipped += 1
                continue
            g.write(out, g.comment or f"{taut} conformer {i}")
            made += 1
        total += made
        print(f"{taut:8s} 앙상블 {len(geoms)} 개 -> 새로 만듦 {made}, 이미 있음 {skipped}")
    print(f"\n총 {total} 개 복원")
    return 0


if __name__ == "__main__":
    sys.exit(main())
