"""
02e_geometry_analysis.py
------------------------
DFT 최적화된 모든 구조의 발색단 기하를 측정해 한 표로 정리한다.

04b 의 aryl_twists() 가 에놀 SMARTS 만 인식해서 디케토 비틀림각이 비어 있었다.
여기서 에놀/디케토 모두를 다루는 정의로 다시 측정하고, xTB 구조와 DFT 구조를
나란히 비교한다. 결과는 results/geometry_analysis.csv 와 콘솔.

측정하는 것
  에놀:   두 아릴 고리의 비틀림 (공액면 기준), 킬레이트 O-H...O 거리
  디케토: 두 아릴케톤의 비틀림, 그리고 O=C-CH2-C=O 골격 이면각(발색단 분리 정도)

실행:  .\scripts\run.ps1 scripts\02e_geometry_analysis.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdDetermineBonds, rdMolTransforms

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import CALCULATIONS, CONFORMERS, RESULTS, read_xyz

RDLogger.DisableLog("rdApp.*")


def mol_from_geom(geom):
    tmp = Path.cwd() / "_geomtmp.xyz"
    geom.write(tmp)
    try:
        raw = Chem.MolFromXYZFile(str(tmp))
        if raw is None:
            return None
        m = Chem.Mol(raw)
        rdDetermineBonds.DetermineBonds(m, charge=0)
        return m
    except Exception:                                    # noqa: BLE001
        return None
    finally:
        tmp.unlink(missing_ok=True)


def fold(angle: float) -> float:
    """0/180 을 같은 평면으로 보고 0~90 으로 접는다."""
    v = abs(angle) % 180.0
    return round(min(v, 180.0 - v), 1)


def measure(mol) -> dict:
    conf = mol.GetConformer()
    p = lambda i: np.array(conf.GetAtomPosition(i))
    tors = lambda a, b, c, d: fold(rdMolTransforms.GetDihedralDeg(conf, a, b, c, d))
    out: dict = {}

    # 에놀: O=C-C=C-OH
    enol = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[OX1]=[CX3]([c])[CX3]=[CX3]([c])[OX2H]"))
    if enol:
        o_ket, c_ket, ar_ket, c_a, c_enol, ar_enol, o_enol = enol[0]
        n1 = [x.GetIdx() for x in mol.GetAtomWithIdx(ar_ket).GetNeighbors()
              if x.GetIdx() != c_ket and x.GetSymbol() == "C"]
        n2 = [x.GetIdx() for x in mol.GetAtomWithIdx(ar_enol).GetNeighbors()
              if x.GetIdx() != c_enol and x.GetSymbol() == "C"]
        if n1:
            out["aryl_ketone_twist_deg"] = tors(o_ket, c_ket, ar_ket, n1[0])
        if n2:
            out["aryl_enol_twist_deg"] = tors(o_enol, c_enol, ar_enol, n2[0])
        hs = [x.GetIdx() for x in mol.GetAtomWithIdx(o_enol).GetNeighbors()
              if x.GetSymbol() == "H"]
        if hs:
            out["chelate_H_to_O_A"] = round(float(np.linalg.norm(
                p(hs[0]) - p(o_ket))), 3)
        out["form"] = "enol"
        return out

    # 디케토: Ar-C(=O)-CH2-C(=O)-Ar
    dk = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[c][CX3](=[OX1])[CX4H2][CX3](=[OX1])[c]"))
    if dk:
        ar1, c1, o1, ch2, c2, o2, ar2 = dk[0]
        n1 = [x.GetIdx() for x in mol.GetAtomWithIdx(ar1).GetNeighbors()
              if x.GetIdx() != c1 and x.GetSymbol() == "C"]
        n2 = [x.GetIdx() for x in mol.GetAtomWithIdx(ar2).GetNeighbors()
              if x.GetIdx() != c2 and x.GetSymbol() == "C"]
        if n1:
            out["aryl1_twist_deg"] = tors(o1, c1, ar1, n1[0])
        if n2:
            out["aryl2_twist_deg"] = tors(o2, c2, ar2, n2[0])
        # 두 카보닐 사이의 골격 비틀림 (발색단이 얼마나 서로 꺾여 있나)
        out["backbone_OCCO_deg"] = tors(o1, c1, c2, o2)
        out["form"] = "diketo"
        return out

    out["form"] = "unknown"
    return out


def main() -> int:
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows = []
    for opt_json in sorted((CALCULATIONS / "01_dft_opt").rglob("result.json")):
        import json
        d = json.loads(opt_json.read_text(encoding="utf-8"))
        if not d.get("ok"):
            continue
        taut, conf_id = d["tautomer"], d["conf_id"]
        xtb_xyz = CONFORMERS / taut / f"{conf_id}.xyz"
        dft_xyz = opt_json.parent / "optimized.xyz"

        m_xtb = mol_from_geom(read_xyz(xtb_xyz)) if xtb_xyz.exists() else None
        m_dft = mol_from_geom(read_xyz(dft_xyz)) if dft_xyz.exists() else None
        g_xtb = measure(m_xtb) if m_xtb else {}
        g_dft = measure(m_dft) if m_dft else {}

        print(f"\n=== {taut} / {conf_id}  ({g_dft.get('form','?')}) ===")
        keys = [k for k in g_dft if k != "form"]
        for k in keys:
            b = g_xtb.get(k)
            a = g_dft.get(k)
            if b is not None and a is not None:
                unit = "A" if k.endswith("_A") else "deg"
                print(f"  {k:24s} xTB {b:7.2f} -> DFT {a:7.2f} {unit} "
                      f"({a-b:+.2f})")
                rows.append({"tautomer": taut, "conformer": conf_id,
                             "form": g_dft["form"], "parameter": k,
                             "xtb": b, "dft": a, "change": round(a - b, 3)})
            elif a is not None:
                print(f"  {k:24s} DFT {a:7.2f}")
                rows.append({"tautomer": taut, "conformer": conf_id,
                             "form": g_dft["form"], "parameter": k,
                             "xtb": "", "dft": a, "change": ""})

    out = RESULTS / "geometry_analysis.csv"
    with out.open("w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=["tautomer", "conformer", "form",
                                           "parameter", "xtb", "dft", "change"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n저장: {out.relative_to(RESULTS.parent)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
