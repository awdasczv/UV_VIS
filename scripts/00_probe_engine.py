"""
00_probe_engine.py
------------------
계산 엔진(Psi4)이 "실제로" 무엇을 지원하는지 프로그램 출력으로 확인한다.
요구사항 9번: 함수/기저셋을 문서나 프로그램 출력으로 확인한 뒤 선택할 것.

확인 항목
  1) Psi4 버전, 애드온(pcmsolver/ddx/dftd3/dftd4 등) 활성 여부
  2) 후보 범함수(functional)들이 실제로 build 되어 있는지
  3) 후보 기저셋을 아보벤존에 대해 실제로 만들 수 있는지 + 기저함수 개수
  4) 작은 분자로 TD-DFT(TDA/RPA) 가 도는지
  5) 작은 분자로 PCM + SCF, 그리고 PCM + TD-DFT 조합이 되는지 (핵심 확인)
결과는 logs/engine_probe.json 과 콘솔에 남는다.
"""
from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

import psi4

ROOT = Path(__file__).resolve().parent.parent
LOGS = ROOT / "logs"
LOGS.mkdir(exist_ok=True)

# 주의: Windows 에서 psi4.core.be_quiet() 는 /dev/null 을 열려다 죽는다.
#       (RuntimeError: PsiOutStream: Failed to open file /dev/null)
#       대신 출력 파일을 명시적으로 지정한다.
psi4.core.set_output_file(str(LOGS / "psi4_probe.out"), False)
psi4.set_memory("2 GB")
psi4.set_num_threads(4)

report: dict = {"psi4_version": psi4.__version__}

# --------------------------------------------------- 1) 애드온 상태
print(f"Psi4 version: {psi4.__version__}")
addons = {}
for name in ["pcmsolver", "ddx", "dftd3", "dftd4", "s-dftd3", "gdma", "cppe",
             "libefp", "chemps2", "simint", "brianqc", "adcc", "cfour", "v2rdm_casscf"]:
    try:
        addons[name] = bool(psi4.addons(name))
    except Exception as exc:                       # noqa: BLE001
        addons[name] = f"query failed: {exc}"
report["addons"] = addons
print("애드온:", json.dumps(addons, ensure_ascii=False))

# --------------------------------------------------- 2) 범함수 지원 여부
CANDIDATE_FUNCTIONALS = [
    "B3LYP", "B3LYP-D3BJ", "PBE0", "PBE0-D3BJ",
    "CAM-B3LYP", "CAM-B3LYP-D3BJ",
    "WB97X-D", "WB97X-D3BJ", "WB97X-V", "WB97M-V",
    "M06-2X", "BHANDHLYP", "LRC-WPBEH", "HSE06",
]
h2o = psi4.geometry("""
0 1
O  0.000000  0.000000  0.117300
H  0.000000  0.757200 -0.469200
H  0.000000 -0.757200 -0.469200
units angstrom
symmetry c1
no_reorient
no_com
""")
psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "print": 0})

func_ok = {}
for f in CANDIDATE_FUNCTIONALS:
    try:
        psi4.energy(f"scf/sto-3g", dft_functional=f, molecule=h2o) if False else None
        e = psi4.energy(f, molecule=h2o)
        func_ok[f] = {"available": True, "test_energy_hartree": round(e, 8)}
    except Exception as exc:                       # noqa: BLE001
        func_ok[f] = {"available": False, "error": str(exc).splitlines()[0][:200]}
    psi4.core.clean()
report["functionals"] = func_ok
print("\n범함수:")
for k, v in func_ok.items():
    print(f"  {k:16s} {'OK' if v['available'] else 'X  ' + v.get('error','')[:90]}")

# --------------------------------------------------- 3) 기저셋 크기 (실제 분자)
CANDIDATE_BASES = ["6-31G(d)", "6-31+G(d)", "6-311G(d,p)",
                   "def2-SVP", "def2-SVPD", "def2-TZVP", "def2-TZVP(-f)", "def2-TZVPD"]
xyz = (ROOT / "structures" / "enolA.xyz")
basis_info = {}
if xyz.exists():
    lines = xyz.read_text(encoding="utf-8").splitlines()
    geom_block = "\n".join(lines[2:])
    mol = psi4.geometry(f"0 1\n{geom_block}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n")
    for b in CANDIDATE_BASES:
        try:
            bs = psi4.core.BasisSet.build(mol, "ORBITAL", b)
            basis_info[b] = {"available": True, "nbf": bs.nbf(),
                             "nshell": bs.nshell(), "puream": bs.has_puream()}
        except Exception as exc:                   # noqa: BLE001
            basis_info[b] = {"available": False, "error": str(exc).splitlines()[0][:200]}
    print("\n기저셋 (아보벤존 enolA, C20H22O3 = 45 원자):")
    for k, v in basis_info.items():
        print(f"  {k:16s} " + (f"nbf={v['nbf']:5d}  shells={v['nshell']}"
                               if v["available"] else "X " + v.get("error", "")[:80]))
else:
    print("\n[건너뜀] structures/enolA.xyz 가 아직 없어 기저셋 크기 측정 불가")
report["basis_sets_avobenzone"] = basis_info

# --------------------------------------------------- 4) TD-DFT 동작 확인
tdtest = {}
try:
    psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "print": 0,
                      "save_jk": True, "tdscf_states": 3, "tdscf_tda": True})
    e, wfn = psi4.energy("B3LYP", molecule=h2o, return_wfn=True)
    res = psi4.procrouting.response.scf_response.tdscf_excitations(
        wfn, states=3, tda=True)
    tdtest["gas_TDA"] = {"ok": True,
                         "excitation_energies_eV": [r["EXCITATION ENERGY"] * 27.211386
                                                    for r in res]}
except Exception as exc:                           # noqa: BLE001
    tdtest["gas_TDA"] = {"ok": False, "error": traceback.format_exc()[-600:]}
psi4.core.clean()

try:
    psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "print": 0, "save_jk": True})
    e, wfn = psi4.energy("B3LYP", molecule=h2o, return_wfn=True)
    res = psi4.procrouting.response.scf_response.tdscf_excitations(
        wfn, states=3, tda=False)
    tdtest["gas_RPA"] = {"ok": True,
                         "excitation_energies_eV": [r["EXCITATION ENERGY"] * 27.211386
                                                    for r in res]}
except Exception as exc:                           # noqa: BLE001
    tdtest["gas_RPA"] = {"ok": False, "error": traceback.format_exc()[-600:]}
psi4.core.clean()

# --------------------------------------------------- 5) PCM 확인 (핵심)
PCM_BLOCK = """
Units = Angstrom
Medium {
    SolverType = IEFPCM
    Solvent = Ethanol
}
Cavity {
    RadiiSet = Bondi
    Type = GePol
    Scaling = True
    Area = 0.3
    Mode = Implicit
}
"""
pcm = {}
try:
    psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "print": 0, "pcm": True,
                      "pcm_scf_type": "total"})
    psi4.pcm_helper(PCM_BLOCK)
    e_pcm = psi4.energy("B3LYP", molecule=h2o)
    pcm["scf_pcm"] = {"ok": True, "energy_hartree": round(e_pcm, 8)}
except Exception as exc:                           # noqa: BLE001
    pcm["scf_pcm"] = {"ok": False, "error": traceback.format_exc()[-800:]}
psi4.core.clean()

try:
    psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "print": 0, "pcm": True,
                      "pcm_scf_type": "total", "save_jk": True})
    psi4.pcm_helper(PCM_BLOCK)
    e, wfn = psi4.energy("B3LYP", molecule=h2o, return_wfn=True)
    res = psi4.procrouting.response.scf_response.tdscf_excitations(
        wfn, states=3, tda=True)
    pcm["tddft_pcm"] = {"ok": True,
                        "excitation_energies_eV": [r["EXCITATION ENERGY"] * 27.211386
                                                   for r in res]}
except Exception as exc:                           # noqa: BLE001
    pcm["tddft_pcm"] = {"ok": False, "error": traceback.format_exc()[-800:]}
psi4.core.clean()

# DDX (ddCOSMO) 도 확인
try:
    psi4.set_options({"basis": "sto-3g", "scf_type": "pk", "print": 0,
                      "ddx": True, "ddx_model": "cosmo",
                      "ddx_solvent": "ethanol", "ddx_radii_set": "bondi"})
    e_ddx = psi4.energy("B3LYP", molecule=h2o)
    pcm["scf_ddx"] = {"ok": True, "energy_hartree": round(e_ddx, 8)}
except Exception as exc:                           # noqa: BLE001
    pcm["scf_ddx"] = {"ok": False, "error": str(exc).splitlines()[-1][:300]}
psi4.core.clean()

report["tdscf"] = tdtest
report["solvation"] = pcm

print("\nTD-DFT / 용매 모델:")
print(f"  기체상 TDA : {'OK' if tdtest['gas_TDA']['ok'] else 'X'}")
print(f"  기체상 RPA : {'OK' if tdtest['gas_RPA']['ok'] else 'X'}")
print(f"  SCF + PCM  : {'OK' if pcm['scf_pcm']['ok'] else 'X'}")
print(f"  TDDFT + PCM: {'OK' if pcm['tddft_pcm']['ok'] else 'X'}")
print(f"  SCF + DDX  : {'OK' if pcm['scf_ddx']['ok'] else 'X'}")
if not pcm["tddft_pcm"]["ok"]:
    print("  -> TDDFT+PCM 오류 요약:", pcm["tddft_pcm"]["error"].splitlines()[-1][:200])

out = LOGS / "engine_probe.json"
out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"\n저장: {out}")
