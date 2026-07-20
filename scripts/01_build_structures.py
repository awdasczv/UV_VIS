"""
01_build_structures.py
----------------------
inputs/tautomers.json 의 SMILES 로부터 아보벤존 3개 토토머의 3D 구조를 만든다.

하는 일
  1) SMILES -> 분자 (수소 명시적으로 추가)
  2) ETKDGv3 로 여러 개의 초기 3D 구조 embed
  3) MMFF94s 로 사전 최적화, 가장 낮은 것 선택
  4) 검증: 분자식 / 원자 개수 / 결합 차수 / 수소 위치 / 에놀 O-H...O 분자내 수소결합 거리
  5) structures/<id>.xyz 저장 + structures/build_report.json 기록

실행:  <env>/python.exe scripts/01_build_structures.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdMolTransforms

RDLogger.DisableLog("rdApp.warning")

# 에놀 킬레이트 골격:  O=C - C(H) = C - O(H)
#                     Oket Cket  Ca    Cenol Oenol
CHELATE_SMARTS = Chem.MolFromSmarts("[OX1]=[CX3][CX3H0,CX3H1]=[CX3][OX2H]")

ROOT = Path(__file__).resolve().parent.parent
SPEC = ROOT / "inputs" / "tautomers.json"
OUTDIR = ROOT / "structures"
N_EMBED = 40
SEED = 20260720


def chelate_atoms(mol):
    """에놀 킬레이트 골격 원자 인덱스 (Oket, Cket, Ca, Cenol, Oenol, H) 반환."""
    m = mol.GetSubstructMatches(CHELATE_SMARTS)
    if not m:
        return None
    o_ket, c_ket, c_a, c_enol, o_enol = m[0]
    h = [n.GetIdx() for n in mol.GetAtomWithIdx(o_enol).GetNeighbors()
         if n.GetSymbol() == "H"]
    if not h:
        return None
    return o_ket, c_ket, c_a, c_enol, o_enol, h[0]


def enforce_chelate(mol, conf_id: int) -> bool:
    """
    6원 킬레이트 고리 (O-H...O=C) 가 닫히도록 회전 가능한 두 이면각을 명시적으로 세팅한다.

    ETKDG/MMFF 는 분자내 수소결합을 잘 못 잡아서 s-trans (열린) 형태를 내놓는 경우가 많다.
    (2026-07-20 실제로 O-H...O = 5.3 A 인 열린 구조가 나옴 -> 이 단계 추가함)

    고리가 평면 6원 고리이므로 고리를 따라가는 이면각은 모두 ~0 도가 되어야 한다:
      Cenol - Ca  - Cket - Oket   = 0   (Ca-Cket 결합 회전)
      Ca    - Cenol - Oenol - H   = 0   (Cenol-Oenol 결합 회전)
    Ca=Cenol 이중결합의 cis/trans 는 SMILES 에서 이미 지정되어 있으므로 건드리지 않는다.
    """
    idx = chelate_atoms(mol)
    if idx is None:
        return False
    o_ket, c_ket, c_a, c_enol, o_enol, h = idx
    conf = mol.GetConformer(conf_id)
    rdMolTransforms.SetDihedralDeg(conf, c_enol, c_a, c_ket, o_ket, 0.0)
    rdMolTransforms.SetDihedralDeg(conf, c_a, c_enol, o_enol, h, 0.0)
    return True


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
    if len(cids) == 0:
        raise RuntimeError("3D embedding 실패")

    idx = chelate_atoms(mol)
    props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")

    for cid in cids:
        if idx is None:
            AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s",
                                         confId=cid, maxIters=2000)
            continue
        o_ket, c_ket, c_a, c_enol, o_enol, h = idx
        enforce_chelate(mol, cid)
        # 1단계: 킬레이트 이면각을 묶어 두고 나머지를 이완
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
        ff.MMFFAddTorsionConstraint(c_enol, c_a, c_ket, o_ket, False, -5.0, 5.0, 200.0)
        ff.MMFFAddTorsionConstraint(c_a, c_enol, o_enol, h, False, -5.0, 5.0, 200.0)
        ff.Minimize(maxIts=2000)
        # 2단계: 구속 없이 자유 최적화 (킬레이트 골짜기 안에서 출발)
        AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s",
                                     confId=cid, maxIters=2000)

    def conf_energy(cid: int) -> float:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
        return ff.CalcEnergy()

    def hb_dist(cid: int) -> float:
        if idx is None:
            return 0.0
        conf = mol.GetConformer(cid)
        p = lambda i: np.array(conf.GetAtomPosition(i))
        return float(np.linalg.norm(p(idx[5]) - p(idx[0])))

    scored = [(conf_energy(c), c) for c in cids]
    if idx is not None:
        # 킬레이트가 실제로 닫힌 구조만 후보로 삼는다 (H...O < 2.2 A)
        closed = [(e, c) for e, c in scored if hb_dist(c) < 2.2]
        if closed:
            scored = closed
        else:
            print("    [경고] 킬레이트가 닫힌 구조를 하나도 못 찾음 -> 최저에너지 구조 사용")
    scored.sort()
    best_e, best_cid = scored[0]

    # 최저 conformer 만 남긴다
    keep = Chem.Mol(mol)
    keep.RemoveAllConformers()
    keep.AddConformer(mol.GetConformer(best_cid), assignId=True)
    return keep, best_e, len(cids)


def hbond_check(mol) -> dict | None:
    """에놀 킬레이트의 기하 검증: O-H...O 거리, O...O 거리, C=C 입체, 고리 평면성."""
    idx = chelate_atoms(mol)
    if idx is None:
        return None
    o_ket, c_ket, c_a, c_enol, o_enol, h = idx
    conf = mol.GetConformer()
    p = lambda i: np.array(conf.GetAtomPosition(i))
    return {
        "atoms": {"O_ketone": o_ket, "C_ketone": c_ket, "C_alpha": c_a,
                  "C_enol": c_enol, "O_enol": o_enol, "H_enol": h},
        "H_to_acceptorO_A": round(float(np.linalg.norm(p(h) - p(o_ket))), 3),
        "O_to_O_A": round(float(np.linalg.norm(p(o_enol) - p(o_ket))), 3),
        # C=C 이중결합 입체: OH 와 케톤이 cis 여야 킬레이트가 가능 (|각| ~ 0)
        "dihedral_Oenol_Cenol_Ca_Cket_deg": round(
            rdMolTransforms.GetDihedralDeg(conf, o_enol, c_enol, c_a, c_ket), 1),
        # 고리 골격 평면성
        "dihedral_Cenol_Ca_Cket_Oket_deg": round(
            rdMolTransforms.GetDihedralDeg(conf, c_enol, c_a, c_ket, o_ket), 1),
    }


def validate(mol, spec_tauto, targets) -> tuple[dict, list[str]]:
    """구조 검증. (리포트, 문제 목록) 반환"""
    problems: list[str] = []
    formula = rdMolDescriptors.CalcMolFormula(mol)
    n_heavy = mol.GetNumHeavyAtoms()
    n_h = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H")

    if formula != targets["formula"]:
        problems.append(f"분자식 불일치: {formula} != {targets['formula']}")
    if n_heavy != targets["heavy_atoms"]:
        problems.append(f"무거운 원자 수 불일치: {n_heavy} != {targets['heavy_atoms']}")
    if n_h != targets["hydrogens"]:
        problems.append(f"수소 수 불일치: {n_h} != {targets['hydrogens']}")

    # 원자별 원자가 검사 (RDKit sanitize 를 이미 통과했지만 명시적으로 한번 더)
    for a in mol.GetAtoms():
        if a.GetSymbol() == "C" and a.GetTotalValence() != 4:
            problems.append(f"C{a.GetIdx()} 원자가={a.GetTotalValence()} (4 이어야 함)")
        if a.GetSymbol() == "O" and a.GetTotalValence() != 2:
            problems.append(f"O{a.GetIdx()} 원자가={a.GetTotalValence()} (2 이어야 함)")
        if a.GetFormalCharge() != 0:
            problems.append(f"원자 {a.GetIdx()} 형식전하 != 0")

    # 작용기 개수
    n_ketone = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX1]=[CX3]")))
    n_enolOH = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX3]=[CX3]-[OX2H]"))) + \
               len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H]-[CX3]=[CX3]")))
    n_sp3_ch2 = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4H2](C=O)C=O")))

    is_enol = spec_tauto["id"].startswith("enol")
    if is_enol:
        if n_ketone != 1:
            problems.append(f"에놀형인데 C=O 개수가 {n_ketone} (1 이어야 함)")
        if n_enolOH < 1:
            problems.append("에놀형인데 C=C-OH 패턴을 찾지 못함")
        hb = hbond_check(mol)
        lo, hi = targets["enol_OH_O_distance_range_A"]
        if hb is None:
            problems.append("킬레이트 골격(O=C-C=C-OH)을 찾을 수 없음")
        else:
            if abs(hb["dihedral_Oenol_Cenol_Ca_Cket_deg"]) > 30.0:
                problems.append(
                    f"C=C 입체 오류: OH 와 케톤이 cis 가 아님 "
                    f"(이면각 {hb['dihedral_Oenol_Cenol_Ca_Cket_deg']} deg). "
                    "SMILES 의 / \\ 방향을 확인할 것.")
            if not (lo <= hb["H_to_acceptorO_A"] <= hi):
                problems.append(
                    f"O-H...O 거리 {hb['H_to_acceptorO_A']} A 가 킬레이트 범위 {lo}-{hi} A 밖 "
                    "(MMFF 는 분자내 수소결합을 과소평가함 -> xTB/DFT 후 재확인)")
    else:
        hb = None
        if n_ketone != 2:
            problems.append(f"디케토형인데 C=O 개수가 {n_ketone} (2 이어야 함)")
        if n_sp3_ch2 != 1:
            problems.append(f"디케토형인데 sp3 CH2 브리지가 {n_sp3_ch2} 개 (1 이어야 함)")

    report = {
        "formula": formula,
        "n_heavy_atoms": n_heavy,
        "n_hydrogens": n_h,
        "mol_weight": round(Descriptors.MolWt(mol), 4),
        "n_carbonyl": n_ketone,
        "n_enol_OH": n_enolOH,
        "n_sp3_CH2_bridge": n_sp3_ch2,
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
    spec = json.loads(SPEC.read_text(encoding="utf-8"))
    targets = spec["validation_targets"]
    OUTDIR.mkdir(parents=True, exist_ok=True)

    all_reports, any_problem = {}, False
    for t in spec["tautomers"]:
        print(f"\n=== {t['id']}: {t['label']} ===")
        print(f"    SMILES: {t['smiles']}")
        mol, e, n = build_3d(t["smiles"])
        print(f"    ETKDG conformer {n} 개 생성, MMFF94s 최저에너지 = {e:.3f} kcal/mol")

        rep, probs = validate(mol, t, targets)
        for k, v in rep.items():
            print(f"    {k}: {v}")
        if probs:
            any_problem = True
            print("    [검증 경고]")
            for p in probs:
                print(f"      - {p}")
        else:
            print("    [검증 통과]")

        xyz = OUTDIR / f"{t['id']}.xyz"
        write_xyz(mol, xyz, f"avobenzone {t['id']} | MMFF94s preopt | {t['smiles']}")
        print(f"    -> {xyz.relative_to(ROOT)}")

        rep.update({"smiles": t["smiles"], "label": t["label"],
                    "mmff94s_energy_kcalmol": round(e, 4),
                    "n_embedded": n, "problems": probs, "xyz": str(xyz.relative_to(ROOT))})
        all_reports[t["id"]] = rep

    out = OUTDIR / "build_report.json"
    out.write_text(json.dumps(all_reports, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n검증 리포트 저장: {out.relative_to(ROOT)}")
    print("전체 결과:", "경고 있음 (위 목록 확인)" if any_problem else "모든 구조 검증 통과")
    return 0


if __name__ == "__main__":
    sys.exit(main())
