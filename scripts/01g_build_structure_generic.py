"""
01g_build_structure_generic.py
------------------------------
분자-독립적 3D 구조 생성기.  molecules/<name>/config.json 의 SMILES 로부터
각 종(species)의 3D 구조를 만든다.  아보벤존 전용 01_build_structures.py 와 달리
킬레이트/토토머를 가정하지 않는다 (DHHB 같은 비-토토머 분자용).

하는 일
  1) SMILES -> 분자 (수소 명시적 추가)
  2) ETKDGv3 로 여러 초기 3D 구조 embed
  3) MMFF94s 로 사전 최적화, 최저에너지 선택
  4) 검증: 분자식 / 무거운원자수 / 수소수 / 원자가 / 형식전하
  5) (있으면) 오쏘-히드록시케톤 등 분자내 수소결합 거리 리포트 (강제하지 않음)
  6) molecules/<name>/structures/<id>.xyz + build_report.json 저장

분자는 환경변수 UV_MOLECULE 로 고른다.
실행:  .\scripts\run.ps1 -Molecule dhhb scripts\01g_build_structure_generic.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.warning")

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import MOL_CONFIG, MOLECULE, ROOT, STRUCTURES  # noqa: E402

N_EMBED = 40
SEED = 20260723


def build_3d(smiles: str):
    """SMILES -> MMFF94s 로 사전최적화된 최저에너지 3D conformer."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"SMILES 파싱 실패: {smiles}")
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    params.useSmallRingTorsions = True
    cids = list(AllChem.EmbedMultipleConfs(mol, numConfs=N_EMBED, params=params))
    if not cids:
        raise RuntimeError("3D embedding 실패")

    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    for cid in cids:
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", confId=cid, maxIters=2000)

    def energy(cid: int) -> float:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
        return ff.CalcEnergy()

    scored = sorted((energy(c), c) for c in cids)
    best_e, best_cid = scored[0]

    keep = Chem.Mol(mol)
    keep.RemoveAllConformers()
    keep.AddConformer(mol.GetConformer(best_cid), assignId=True)
    return keep, best_e, len(cids)


def hbond_report(mol, smarts: str, rng) -> dict | None:
    """donor O-H 와 억셉터 O=C 사이 최단 H...O 거리를 찾아 리포트 (강제하지 않음).

    smarts 는 도너 O(H) 로 시작하는 패턴을 가정한다 (예: 오쏘 히드록시아릴케톤).
    억셉터는 분자 내 모든 카보닐 산소 중 H 에 가장 가까운 것으로 잡는다.
    """
    patt = Chem.MolFromSmarts(smarts) if smarts else None
    if patt is None:
        return None
    matches = mol.GetSubstructMatches(patt)
    if not matches:
        return None
    conf = mol.GetConformer()
    p = lambda i: np.array(conf.GetAtomPosition(i))

    # 도너 히드록시 O 와 그 H
    donor_o = matches[0][0]
    h_idx = [n.GetIdx() for n in mol.GetAtomWithIdx(donor_o).GetNeighbors()
             if n.GetSymbol() == "H"]
    if not h_idx:
        return None
    h = h_idx[0]

    # 카보닐 산소들 (C=O)
    carbonyl_o = [m[0] for m in
                  mol.GetSubstructMatches(Chem.MolFromSmarts("[OX1]=[CX3]"))]
    if not carbonyl_o:
        return None
    dists = [(float(np.linalg.norm(p(h) - p(o))), o) for o in carbonyl_o]
    d, acc_o = min(dists)
    lo, hi = rng
    return {
        "donor_O_idx": int(donor_o), "donor_H_idx": int(h), "acceptor_O_idx": int(acc_o),
        "H_to_acceptorO_A": round(d, 3),
        "in_expected_range": bool(lo <= d <= hi),
        "expected_range_A": [lo, hi],
    }


def validate(mol, targets) -> tuple[dict, list[str]]:
    problems: list[str] = []
    formula = rdMolDescriptors.CalcMolFormula(mol)
    n_heavy = mol.GetNumHeavyAtoms()
    n_h = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H")

    if "formula" in targets and formula != targets["formula"]:
        problems.append(f"분자식 불일치: {formula} != {targets['formula']}")
    if "heavy_atoms" in targets and n_heavy != targets["heavy_atoms"]:
        problems.append(f"무거운 원자 수 불일치: {n_heavy} != {targets['heavy_atoms']}")
    if "hydrogens" in targets and n_h != targets["hydrogens"]:
        problems.append(f"수소 수 불일치: {n_h} != {targets['hydrogens']}")

    expected_valence = {"C": 4, "N": 3, "O": 2}
    for a in mol.GetAtoms():
        sym = a.GetSymbol()
        if sym in expected_valence and a.GetTotalValence() != expected_valence[sym]:
            problems.append(f"{sym}{a.GetIdx()} 원자가={a.GetTotalValence()} "
                            f"({expected_valence[sym]} 기대)")
        if a.GetFormalCharge() != 0:
            problems.append(f"원자 {a.GetIdx()} 형식전하 != 0")

    hb = None
    if targets.get("hbond_donor_smarts"):
        hb = hbond_report(mol, targets["hbond_donor_smarts"],
                          targets.get("hbond_OH_O_distance_range_A", [1.4, 2.2]))
        if hb and not hb["in_expected_range"]:
            problems.append(
                f"분자내 수소결합 O-H...O={hb['H_to_acceptorO_A']} A 가 기대범위 밖 "
                "(MMFF 는 과소평가 -> xTB/DFT 후 재확인)")

    report = {
        "formula": formula, "n_heavy_atoms": n_heavy, "n_hydrogens": n_h,
        "mol_weight": round(Descriptors.MolWt(mol), 4),
        "intramolecular_hbond": hb,
        "canonical_smiles_from_3d": Chem.MolToSmiles(Chem.RemoveHs(mol)),
    }
    return report, problems


def write_xyz(mol, path: Path, comment: str) -> None:
    conf = mol.GetConformer()
    lines = [str(mol.GetNumAtoms()), comment]
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():2s} {pos.x:14.8f} {pos.y:14.8f} {pos.z:14.8f}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    spec = json.loads(MOL_CONFIG.read_text(encoding="utf-8"))
    targets = spec.get("validation_targets", {})
    STRUCTURES.mkdir(parents=True, exist_ok=True)

    all_reports, any_problem = {}, False
    for t in spec["tautomers"]:
        print(f"\n=== {t['id']}: {t['label']} ===")
        print(f"    SMILES: {t['smiles']}")
        mol, e, n = build_3d(t["smiles"])
        print(f"    ETKDG conformer {n} 개 생성, MMFF94s 최저에너지 = {e:.3f} kcal/mol")

        rep, probs = validate(mol, targets)
        for k, v in rep.items():
            print(f"    {k}: {v}")
        if probs:
            any_problem = True
            print("    [검증 경고]")
            for p in probs:
                print(f"      - {p}")
        else:
            print("    [검증 통과]")

        xyz = STRUCTURES / f"{t['id']}.xyz"
        write_xyz(mol, xyz, f"{MOLECULE} {t['id']} | MMFF94s preopt | {t['smiles']}")
        print(f"    -> {xyz.relative_to(ROOT)}")

        rep.update({"smiles": t["smiles"], "label": t["label"],
                    "mmff94s_energy_kcalmol": round(e, 4),
                    "n_embedded": n, "problems": probs,
                    "xyz": str(xyz.relative_to(ROOT))})
        all_reports[t["id"]] = rep

    out = STRUCTURES / "build_report.json"
    out.write_text(json.dumps(all_reports, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n검증 리포트 저장: {out.relative_to(ROOT)}")
    print("전체 결과:", "경고 있음 (위 목록 확인)" if any_problem else "모든 구조 검증 통과")
    return 0


if __name__ == "__main__":
    sys.exit(main())
