"""
02c_select_representatives.py
-----------------------------
TD-DFT 를 돌릴 대표 컨포머를 고른다 (요구사항 8번).

왜 단순 볼츠만 순위로 안 뽑는가
--------------------------------
02b 의 분석 결과, 에놀형 컨포머 9~11개는 발색단(공액계) 기하가 사실상 동일했다.
차이는 tert-butyl 메틸기 회전 같은, 흡수 스펙트럼과 무관한 자유도뿐이었다.
그래서 볼츠만 가중치가 1/9 씩 고르게 퍼지고, 상위 3개를 뽑아도 누적 가중치가
38% 밖에 안 된다. 하지만 이는 "앙상블의 62%를 놓쳤다"는 뜻이 아니다.
빠진 컨포머들의 스펙트럼이 뽑힌 것과 같기 때문이다.

따라서 다음 두 단계로 고른다.
  1) 발색단 비틀림각으로 컨포머를 클러스터링한다 (스펙트럼상 구별되는 종류).
  2) 각 클러스터에서 최저에너지 구조를 대표로 뽑고, 클러스터 가중치가 큰 순으로
     목표 누적 가중치(기본 0.90)를 채운다. 최대 개수 제한도 그대로 지킨다.

이렇게 하면 계산 수는 줄면서 앙상블 커버리지는 오히려 올라간다.
비교를 위해 기존 볼츠만 순위 선택 결과는 selected_boltzmann.json 으로 보존한다.

실행:  .\scripts\run.ps1 scripts\02c_select_representatives.py
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import CONFORMERS, load_checkpoint, save_checkpoint


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-weight", type=float, default=0.90)
    ap.add_argument("--max-select", type=int, default=3)
    ap.add_argument("--intra-cluster-check", action="store_true",
                    help="클러스터가 1개뿐인 토토머에 대해, 같은 클러스터 안에서 "
                         "발색단 기하가 가장 다른 컨포머를 하나 더 넣는다. "
                         "'클러스터 안에서는 스펙트럼이 정말 같은가'를 실증하기 위함.")
    args = ap.parse_args()

    val = load_checkpoint(CONFORMERS / "validation_report.json")
    if not val:
        print("먼저 02b_validate_conformers.py 를 실행하세요.")
        return 1

    summary = {}
    for taut, info in val.items():
        clusters = info.get("chromophore_clusters") or []
        if not clusters:
            print(f"[건너뜀] {taut}: 클러스터 정보 없음")
            continue

        by_w = sorted(clusters, key=lambda c: -c["total_weight"])
        picked, cum = [], 0.0
        for c in by_w:
            picked.append(c)
            cum += c["total_weight"]
            if cum >= args.target_weight or len(picked) >= args.max_select:
                break

        # 대표 컨포머의 상세 정보 찾기
        detail = {e["conf_id"]: e for e in info["conformers"]}
        selected = []
        for c in picked:
            d = detail[c["representative"]]
            selected.append({
                "conf_id": c["representative"],
                "rel_energy_kcalmol": d["rel_energy_kcalmol"],
                # 이 대표가 대변하는 클러스터 전체의 가중치를 쓴다
                "boltzmann_weight": c["total_weight"],
                "own_boltzmann_weight": d["boltzmann_weight"],
                "cluster": c["cluster"],
                "cluster_size": len(c["members"]),
                "chromophore_dihedrals": c["dihedrals"],
                "xyz": f"conformers/{taut}/{c['representative']}.xyz",
            })
        # 클러스터가 하나뿐이면, 같은 클러스터 안에서 발색단 기하가 가장 다른
        # 컨포머를 하나 더 넣어 '정말 같은 스펙트럼인가'를 확인한다.
        if args.intra_cluster_check and len(clusters) == 1:
            rep_id = selected[0]["conf_id"]
            rep_d = selected[0]["chromophore_dihedrals"] or {}

            def dist(e):
                d = e.get("chromophore_dihedrals") or {}
                keys = [k for k in rep_d if rep_d.get(k) is not None
                        and d.get(k) is not None]
                return max((abs(rep_d[k] - d[k]) for k in keys), default=0.0)

            cands = [e for e in info["conformers"] if e["conf_id"] != rep_id]
            if cands:
                far = max(cands, key=dist)
                selected.append({
                    "conf_id": far["conf_id"],
                    "rel_energy_kcalmol": far["rel_energy_kcalmol"],
                    "boltzmann_weight": far["boltzmann_weight"],
                    "own_boltzmann_weight": far["boltzmann_weight"],
                    "cluster": 0,
                    "cluster_size": 1,
                    "role": "intra_cluster_check",
                    "max_dihedral_diff_from_rep_deg": round(dist(far), 1),
                    "chromophore_dihedrals": far.get("chromophore_dihedrals"),
                    "xyz": f"conformers/{taut}/{far['conf_id']}.xyz",
                })
                # 대표의 가중치에서 이 구조 몫을 떼어낸다
                selected[0]["boltzmann_weight"] -= far["boltzmann_weight"]

        # 가중치 재규격화 (뽑힌 클러스터들 안에서 합이 1이 되도록)
        tot = sum(s["boltzmann_weight"] for s in selected)
        for s in selected:
            s["weight_normalized"] = round(s["boltzmann_weight"] / tot, 5)

        old = CONFORMERS / taut / "selected.json"
        if old.exists() and not (CONFORMERS / taut / "selected_boltzmann.json").exists():
            shutil.copy(old, CONFORMERS / taut / "selected_boltzmann.json")

        out = {
            "tautomer": taut,
            "selection_method": "발색단 기하 클러스터별 대표 (에너지 최저)",
            "rationale": ("컨포머들의 볼츠만 가중치가 고르게 퍼지는 이유는 "
                          "스펙트럼과 무관한 알킬기 회전 때문이다. 발색단 기하로 묶어야 "
                          "적은 계산으로 앙상블을 제대로 대표할 수 있다."),
            "n_unique_conformers": info["n_conformers"],
            "n_chromophore_clusters": len(clusters),
            "target_cumulative_weight": args.target_weight,
            "max_conformers": args.max_select,
            "achieved_cumulative_weight": round(cum, 4),
            "selected": selected,
        }
        save_checkpoint(old, out)

        print(f"\n=== {taut} ===")
        print(f"  고유 컨포머 {info['n_conformers']} 개 -> 발색단 기하 {len(clusters)} 종류")
        print(f"  대표 {len(selected)} 개 선택, 누적 가중치 = {cum:.3f}")
        for s in selected:
            d = s["chromophore_dihedrals"] or {}
            ds = ", ".join(f"{k.replace('_deg','')}={v}" for k, v in d.items()
                           if v is not None)
            print(f"    {s['conf_id']:16s} 클러스터{s['cluster']} "
                  f"({s['cluster_size']}개 대표)  w={s['boltzmann_weight']:.3f}  {ds}")
        summary[taut] = {"n_clusters": len(clusters),
                         "cum_weight": round(cum, 4),
                         "selected": [s["conf_id"] for s in selected]}

    save_checkpoint(CONFORMERS / "selection_summary.json", summary)
    total = sum(len(v["selected"]) for v in summary.values())
    print(f"\n총 {total} 개 구조를 TD-DFT 대상으로 선택했다.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
