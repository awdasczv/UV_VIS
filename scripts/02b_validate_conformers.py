"""
02b_validate_conformers.py
--------------------------
컨포머 앙상블 검증.

xTB 최적화 과정에서 구조가 엉뚱하게 변하지 않았는지 확인한다. 특히
  - 결합 연결(connectivity)이 원래 SMILES 와 같은가?  (토토머가 바뀌지 않았는가)
  - 에놀형에서 킬레이트 O-H...O 수소결합이 유지되는가?
  - 디케토형에서 sp3 CH2 브리지가 유지되는가?
  - 원자가 겹치는 등 물리적으로 이상한 구조는 없는가?

이 검증을 안 하면, 예컨대 xTB 최적화 중 에놀 양성자가 옮겨가 다른 토토머가 되어도
모르고 TD-DFT 를 돌리게 된다.

실행:  .\scripts\run.ps1 scripts\02b_validate_conformers.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import CONFORMERS, MOL_CONFIG, read_xyz, save_checkpoint

RDLogger.DisableLog("rdApp.*")

# 공유결합 반지름 (Angstrom) - 최소 원자간 거리 점검용
COV_R = {"H": 0.31, "C": 0.76, "O": 0.66}


def mol_from_xyz(path: Path):
    """xyz 파일에서 결합을 추정해 RDKit 분자를 만든다."""
    raw = Chem.MolFromXYZFile(str(path))
    if raw is None:
        return None
    mol = Chem.Mol(raw)
    try:
        rdDetermineBonds.DetermineBonds(mol, charge=0)
    except Exception:                                   # noqa: BLE001
        return None
    return mol


def canonical(smiles_or_mol) -> str | None:
    try:
        m = (Chem.MolFromSmiles(smiles_or_mol)
             if isinstance(smiles_or_mol, str) else smiles_or_mol)
        if m is None:
            return None
        # 입체 정보는 빼고 골격/토토머 동일성만 본다
        m = Chem.RemoveHs(m)
        Chem.RemoveStereochemistry(m)
        return Chem.MolToSmiles(m)
    except Exception:                                   # noqa: BLE001
        return None


def check_close_contacts(geom, factor: float = 0.7) -> list[str]:
    """원자간 거리가 공유결합 반지름 합의 factor 배보다 짧으면 이상."""
    bad = []
    n = geom.natoms
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(geom.coords[i] - geom.coords[j]))
            lim = factor * (COV_R.get(geom.symbols[i], 0.7)
                            + COV_R.get(geom.symbols[j], 0.7))
            if d < lim:
                bad.append(f"{geom.symbols[i]}{i}-{geom.symbols[j]}{j} = {d:.2f} A")
    return bad


def chelate_distance(mol, geom) -> float | None:
    """에놀 O-H ... O=C 의 H...O 거리."""
    patt = Chem.MolFromSmarts("[OX1]=[CX3][CX3]=[CX3][OX2H]")
    m = mol.GetSubstructMatches(patt)
    if not m:
        return None
    o_ket, _, _, _, o_enol = m[0]
    # o_enol 에 붙은 H 찾기
    hs = [n.GetIdx() for n in mol.GetAtomWithIdx(o_enol).GetNeighbors()
          if n.GetSymbol() == "H"]
    if not hs:
        return None
    return float(np.linalg.norm(geom.coords[hs[0]] - geom.coords[o_ket]))


def chromophore_dihedrals(mol, geom) -> dict | None:
    """
    발색단(공액계)의 기하를 특징짓는 비틀림각을 잰다.

    왜 중요한가:
      컨포머끼리 에너지가 거의 같아도(예: tert-butyl 회전) 발색단 기하가 같으면
      흡수 스펙트럼도 사실상 같다. 반대로 아릴 고리가 공액면에서 얼마나 비틀렸는지는
      pi 공액 정도를 직접 바꾸므로 lambda_max 에 크게 영향을 준다.
      따라서 "볼츠만 가중치"보다 "발색단 기하가 몇 종류인가"가 스펙트럼에는 더 중요하다.

    반환: 각 아릴 고리의 비틀림각(도, 0 = 공액면과 같은 평면)
    """
    from rdkit.Chem import rdMolTransforms
    conf = mol.GetConformer() if mol.GetNumConformers() else None

    def tors(a, b, c, d):
        if conf is None:
            return None
        v = rdMolTransforms.GetDihedralDeg(conf, a, b, c, d)
        # 0/180 을 같은 평면으로 보고 0~90 으로 접는다
        v = abs(v) % 180.0
        return round(min(v, 180.0 - v), 1)

    out = {}
    # 에놀형: O=C-Ar  와  HO-C=C ... C-Ar
    enol = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[OX1]=[CX3]([c])[CX3]=[CX3]([c])[OX2H]"))
    if enol:
        o_ket, c_ket, ar_ket, c_a, c_enol, ar_enol, o_enol = enol[0]
        ao = [n.GetIdx() for n in mol.GetAtomWithIdx(ar_ket).GetNeighbors()
              if n.GetIdx() != c_ket and n.GetSymbol() == "C"]
        if ao:
            out["aryl_ketone_twist_deg"] = tors(o_ket, c_ket, ar_ket, ao[0])
        ae = [n.GetIdx() for n in mol.GetAtomWithIdx(ar_enol).GetNeighbors()
              if n.GetIdx() != c_enol and n.GetSymbol() == "C"]
        if ae:
            out["aryl_enol_twist_deg"] = tors(o_enol, c_enol, ar_enol, ae[0])
    else:
        # 디케토형: 두 개의 Ar-C(=O)-CH2
        dk = mol.GetSubstructMatches(
            Chem.MolFromSmarts("[c][CX3](=[OX1])[CX4H2][CX3](=[OX1])[c]"))
        if dk:
            ar1, c1, o1, ch2, c2, o2, ar2 = dk[0]
            n1 = [n.GetIdx() for n in mol.GetAtomWithIdx(ar1).GetNeighbors()
                  if n.GetIdx() != c1 and n.GetSymbol() == "C"]
            n2 = [n.GetIdx() for n in mol.GetAtomWithIdx(ar2).GetNeighbors()
                  if n.GetIdx() != c2 and n.GetSymbol() == "C"]
            if n1:
                out["aryl1_twist_deg"] = tors(o1, c1, ar1, n1[0])
            if n2:
                out["aryl2_twist_deg"] = tors(o2, c2, ar2, n2[0])
            out["backbone_C1_CH2_C3_deg"] = tors(c1, ch2, c2, o2)

    # 메톡시기 방향 (Ar-O-CH3 가 고리 평면과 이루는 각)
    ome = mol.GetSubstructMatches(Chem.MolFromSmarts("[c][OX2][CH3]"))
    if ome:
        ar, o, me = ome[0]
        nb = [n.GetIdx() for n in mol.GetAtomWithIdx(ar).GetNeighbors()
              if n.GetIdx() != o and n.GetSymbol() == "C"]
        if nb:
            out["methoxy_twist_deg"] = tors(nb[0], ar, o, me)
    return out or None


def cluster_by_chromophore(entries: list[dict], tol_deg: float = 12.0) -> list[list[int]]:
    """발색단 비틀림각이 tol_deg 이내면 스펙트럼상 같은 것으로 보고 묶는다."""
    keys = []
    for e in entries:
        d = e.get("chromophore_dihedrals") or {}
        keys.append(tuple(sorted((k, v) for k, v in d.items() if v is not None)))
    clusters: list[list[int]] = []
    reps: list[tuple] = []
    for i, k in enumerate(keys):
        placed = False
        for ci, rk in enumerate(reps):
            if len(rk) == len(k) and all(
                    abs(a[1] - b[1]) <= tol_deg for a, b in zip(rk, k)):
                clusters[ci].append(i)
                placed = True
                break
        if not placed:
            clusters.append([i])
            reps.append(k)
    return clusters


def main() -> int:
    spec = json.loads(MOL_CONFIG.read_text(encoding="utf-8"))
    expected = {t["id"]: canonical(t["smiles"]) for t in spec["tautomers"]}

    report, all_ok = {}, True
    for taut in [t["id"] for t in spec["tautomers"]]:
        ens_path = CONFORMERS / taut / "ensemble.json"
        if not ens_path.exists():
            print(f"[건너뜀] {taut}: ensemble.json 없음")
            continue
        ens = json.loads(ens_path.read_text(encoding="utf-8"))
        print(f"\n=== {taut} ({len(ens['conformers'])} 개 컨포머) ===")
        print(f"  기대 골격: {expected[taut]}")

        entries, n_bad = [], 0
        for c in ens["conformers"]:
            xyz = CONFORMERS / taut / f"{c['conf_id']}.xyz"
            if not xyz.exists():
                continue
            geom = read_xyz(xyz)
            mol = mol_from_xyz(xyz)
            got = canonical(mol) if mol is not None else None
            same = (got == expected[taut])
            contacts = check_close_contacts(geom)
            hb = chelate_distance(mol, geom) if mol is not None else None

            problems = []
            if mol is None:
                problems.append("결합 추정 실패")
            elif not same:
                problems.append(f"골격이 다름: {got}")
            if contacts:
                problems.append("원자 겹침: " + ", ".join(contacts[:3]))
            if taut.startswith("enol"):
                if hb is None:
                    problems.append("킬레이트 골격 미검출")
                elif not (1.3 <= hb <= 2.2):
                    problems.append(f"O-H...O = {hb:.2f} A (킬레이트 범위 밖)")

            if problems:
                n_bad += 1
                all_ok = False
                print(f"  [문제] {c['conf_id']}: {'; '.join(problems)}")
            dihe = chromophore_dihedrals(mol, geom) if mol is not None else None
            entries.append({"conf_id": c["conf_id"],
                            "rel_energy_kcalmol": c["rel_energy_kcalmol"],
                            "boltzmann_weight": c["boltzmann_weight"],
                            "skeleton_matches": same,
                            "canonical_smiles": got,
                            "chelate_H_to_O_A": round(hb, 3) if hb else None,
                            "chromophore_dihedrals": dihe,
                            "problems": problems})

        n = len(entries)
        print(f"  검증: {n - n_bad}/{n} 통과")
        if taut.startswith("enol") and n:
            hbs = [e["chelate_H_to_O_A"] for e in entries if e["chelate_H_to_O_A"]]
            if hbs:
                print(f"  O-H...O 거리: 최소 {min(hbs):.3f} / 평균 "
                      f"{sum(hbs)/len(hbs):.3f} / 최대 {max(hbs):.3f} A")

        # 발색단 기하 표
        print(f"  {'컨포머':16s} {'dE':>6s} {'w':>6s}  발색단 비틀림각(도)")
        for e in entries:
            d = e.get("chromophore_dihedrals") or {}
            s = "  ".join(f"{k.replace('_deg','').replace('_twist',''):14s}={v:5.1f}"
                          for k, v in d.items() if v is not None)
            print(f"  {e['conf_id']:16s} {e['rel_energy_kcalmol']:6.2f} "
                  f"{e['boltzmann_weight']:6.3f}  {s}")

        clusters = cluster_by_chromophore(entries)
        print(f"  -> 발색단 기하 기준으로 {len(clusters)} 종류로 묶임")
        cl_info = []
        for ci, idxs in enumerate(clusters):
            w = sum(entries[i]["boltzmann_weight"] for i in idxs)
            rep = min(idxs, key=lambda i: entries[i]["rel_energy_kcalmol"])
            cl_info.append({"cluster": ci, "total_weight": round(w, 4),
                            "representative": entries[rep]["conf_id"],
                            "members": [entries[i]["conf_id"] for i in idxs],
                            "dihedrals": entries[rep].get("chromophore_dihedrals")})
            print(f"     클러스터 {ci}: 가중치 합 {w:.3f}, 대표 {entries[rep]['conf_id']}, "
                  f"{len(idxs)} 개")

        report[taut] = {"expected_skeleton": expected[taut],
                        "n_conformers": n, "n_problems": n_bad,
                        "chromophore_clusters": cl_info,
                        "conformers": entries}

    save_checkpoint(CONFORMERS / "validation_report.json", report)
    print(f"\n검증 리포트 -> conformers/validation_report.json")
    print("전체:", "모두 통과" if all_ok else "문제 있음 (위 목록 확인)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
