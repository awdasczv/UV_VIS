"""
psi4_helpers.py
---------------
Psi4 를 다루면서 반복되는 부분을 모아둔 모듈.

  - Windows 안전한 초기화 (be_quiet() 가 /dev/null 때문에 죽는 문제 회피)
  - PCM 입력 블록 생성
  - 분자 만들기
  - TD-DFT 결과에서 주요 오비탈 전이 뽑아내기
  - 실패 원인 분류 (요구사항 16번: 무작정 낮추지 말고 원인부터 분류)
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import psi4

HARTREE_TO_EV = 27.211386245988
EV_TO_NM = 1239.841984

# PCMSolver 가 아는 용매 이름 (Psi4 문서/PCMSolver 내장 목록)
SOLVENT_NAMES = {
    "ethanol": "Ethanol",
    "methanol": "Methanol",
    "water": "Water",
    "acetonitrile": "Acetonitrile",
    "cyclohexane": "Cyclohexane",
    "dmso": "DimethylSulfoxide",
}


def init_psi4(output_file: Path, memory: str = "6 GB", nthreads: int = 4,
              append: bool = False) -> None:
    """Windows 에서 안전한 Psi4 초기화."""
    output_file.parent.mkdir(parents=True, exist_ok=True)
    psi4.core.set_output_file(str(output_file), append)
    psi4.set_memory(memory)
    psi4.set_num_threads(nthreads)


def pcm_block(solvent: str = "ethanol", nonequilibrium: bool = True,
              area: float = 1.0) -> str:
    """
    IEFPCM 입력 블록.

    Nonequilibrium = True 는 수직 전이용이다. 전자 들뜸은 순간적으로 일어나므로
    용매의 배향(느린 자유도)은 따라오지 못하고 전자 분극(빠른 자유도)만 반응한다.
    Psi4 문서가 수직 들뜸 계산에 이것을 쓰라고 명시하고 있다.

    Area(표면 조각 하나의 넓이) 기본값을 1.0 A^2 로 둔 근거:
      아보벤존에서 직접 수렴 시험을 했다 (calculations/00_test/pcm_cost_probe.json).
        Area 1.0 -> 표면조각 1024개, lambda_max 332.28 nm, TD-DFT 566초
        Area 2.0 -> 표면조각  750개, lambda_max 332.35 nm, TD-DFT 425초
      두 값의 차이가 0.07 nm 에 불과하므로 격자 밀도에 대해 수렴했다.
      Psi4 예제의 0.3 (표면조각 2149개) 은 이 목적에는 불필요하게 촘촘하다.
    """
    name = SOLVENT_NAMES.get(solvent.lower(), solvent)
    neq = "True" if nonequilibrium else "False"
    return f"""
Units = Angstrom
Medium {{
    SolverType = IEFPCM
    Solvent = {name}
    Nonequilibrium = {neq}
}}
Cavity {{
    RadiiSet = Bondi
    Type = GePol
    Scaling = True
    Area = {area}
    Mode = Implicit
}}
"""


def make_molecule(symbols, coords, charge: int = 0, mult: int = 1):
    body = "\n".join(f"{s} {x:.10f} {y:.10f} {z:.10f}"
                     for s, (x, y, z) in zip(symbols, coords))
    return psi4.geometry(
        f"{charge} {mult}\n{body}\nunits angstrom\nsymmetry c1\nno_reorient\nno_com\n")


# ------------------------------------------------- TD-DFT 결과 파싱
def dominant_orbital_transitions(res_entry: dict, nocc: int, n_top: int = 3,
                                 min_weight: float = 0.05) -> list[dict]:
    """
    TD-DFT 한 상태의 고유벡터에서 기여가 큰 오비탈 쌍을 뽑는다.

    Psi4 의 tdscf_excitations 결과에는 'RIGHT EIGENVECTOR ALPHA' 가
    (nocc x nvirt) 행렬로 들어있다. 원소 c_ia 의 제곱이 그 홀-입자 쌍의 기여도다.
    HOMO / LUMO 기준 상대 표기(HOMO-1 -> LUMO+2 등)로 바꿔 돌려준다.
    """
    key = None
    for k in ("RIGHT EIGENVECTOR ALPHA", "RIGHT EIGENVECTOR", "EIGENVECTOR ALPHA"):
        if k in res_entry:
            key = k
            break
    if key is None:
        return []

    vec = res_entry[key]
    try:
        arr = np.asarray(vec.to_array())
    except AttributeError:
        arr = np.asarray(vec)
    if arr.ndim != 2:
        return []

    weights = arr ** 2
    total = weights.sum()
    if total <= 0:
        return []
    weights = weights / total

    flat = np.dstack(np.unravel_index(np.argsort(-weights, axis=None), weights.shape))[0]
    out = []
    for i, a in flat[:n_top]:
        w = float(weights[i, a])
        if w < min_weight:
            break
        # i 는 occupied index (0 = 최하위), a 는 virtual index (0 = LUMO)
        occ_off = nocc - 1 - int(i)          # 0 이면 HOMO
        vir_off = int(a)                     # 0 이면 LUMO
        occ_lbl = "HOMO" if occ_off == 0 else f"HOMO-{occ_off}"
        vir_lbl = "LUMO" if vir_off == 0 else f"LUMO+{vir_off}"
        out.append({"from": occ_lbl, "to": vir_lbl, "weight": round(w, 4)})
    return out


def format_transitions(trans: list[dict]) -> str:
    """주요 전이를 사람이 읽을 문자열로."""
    if not trans:
        return ""
    return "; ".join(f"{t['from']}->{t['to']} ({t['weight']*100:.0f}%)" for t in trans)


def parse_tdscf_results(res, nocc: int) -> list[dict]:
    """tdscf_excitations 결과 리스트를 우리 형식의 dict 리스트로."""
    out = []
    for i, r in enumerate(res, start=1):
        ev = float(r["EXCITATION ENERGY"]) * HARTREE_TO_EV
        f_len = float(r.get("OSCILLATOR STRENGTH (LEN)", 0.0))
        f_vel = float(r.get("OSCILLATOR STRENGTH (VEL)", 0.0))
        orbs = dominant_orbital_transitions(r, nocc)
        out.append({
            "state": i,
            "energy_eV": round(ev, 5),
            "wavelength_nm": round(EV_TO_NM / ev, 3) if ev > 0 else None,
            "osc_strength": round(f_len, 6),
            "osc_strength_velocity": round(f_vel, 6),
            "orbital_transitions": orbs,
            "orbital_transitions_str": format_transitions(orbs),
        })
    return out


# ------------------------------------------------- 실패 원인 분류
FAILURE_PATTERNS = [
    (r"PsiException:\s*Could not converge SCF|SCF (?:iterations )?(?:did not|failed to) converge"
     r"|Could not converge SCF iterations",
     "SCF_NOT_CONVERGED",
     "SCF 가 수렴하지 않음. 대응: SOSCF/second-order SCF 켜기, maxiter 늘리기, "
     "damping/level shift 적용, guess 를 SAD->GWH/CORE 로 바꾸기, 기저셋을 낮춰 얻은 "
     "궤도를 초기 guess 로 재사용."),
    (r"not enough memory|insufficient memory|std::bad_alloc|MemoryError|"
     r"Cannot allocate|out of memory",
     "OUT_OF_MEMORY",
     "메모리 부족. 대응: psi4 memory 설정 낮추기, scf_type DF 사용(이미 사용중이면 "
     "df_ints_io=SAVE), 상태 수 줄이기, 기저셋 줄이기, 스레드 수 줄이기."),
    (r"Davidson.*(?:not converge|failed)|TDSCF.*did not converge|"
     r"maximum number of iterations.*reached",
     "EXCITED_STATE_NOT_CONVERGED",
     "들뜬상태(Davidson) 미수렴. 대응: tdscf_maxiter 늘리기, tdscf_r_convergence 완화, "
     "요청 상태 수 줄이기, TDA(tda=True) 로 바꾸기, guess 벡터 수 늘리기."),
    (r"Atoms.*too close|inter-atomic distance|Nuclear repulsion.*(?:inf|nan)|"
     r"Fatal Error.*geometry",
     "BAD_GEOMETRY",
     "구조 이상(원자 겹침 등). 대응: 입력 xyz 재검증, xTB 로 사전 최적화."),
    (r"BasisSetNotFound|Unable to find a basis set",
     "BASIS_NOT_FOUND",
     "해당 기저셋이 없음. 대응: 00_probe_engine.py 로 사용 가능한 이름 확인."),
    (r"OPTKING|optimization did not converge|Maximum number of steps exceeded",
     "GEOMETRY_OPT_NOT_CONVERGED",
     "구조 최적화 미수렴. 대응: geom_maxiter 늘리기, opt_coordinates 바꾸기, "
     "더 좋은 초기 구조(xTB 최적화 결과) 사용, 수렴 기준 완화."),
    (r"S matrix is not positive-definite|not positive.definite",
     "PCM_ILL_CONDITIONED",
     "PCM 공동 표면의 겹침행렬 S 가 양정치가 아니다. 원래 에러인데 PCMSolver 가 "
     "패치로 경고로 낮춘 것이므로 결과를 그대로 믿으면 안 된다. "
     "대응: Cavity Area 를 바꿔가며(1.0 / 0.6 / 1.5 / 0.3 / 2.0) 경고가 사라지는 값을 "
     "찾고 채택한 값을 기록할 것. 그래도 안 되면 RadiiSet 을 UFF 로 바꾸거나 "
     "Scaling=False 를 시도. 근본 해결은 ORCA 의 CPCM 사용."),
    (r"PCMSOLVER|pcmsolver.*error|cavity",
     "PCM_ERROR",
     "PCM 오류. 대응: cavity Area 키우기, RadiiSet 바꾸기, 용매 이름 철자 확인."),
]


def classify_failure(text: str) -> dict:
    """
    로그/예외 문자열에서 실패 원인을 분류한다.
    요구사항 16번: 실패하면 무조건 수준을 낮추지 말고 먼저 원인을 분류한다.
    """
    for pattern, code, remedy in FAILURE_PATTERNS:
        if re.search(pattern, text, re.IGNORECASE):
            return {"code": code, "remedy": remedy}
    return {"code": "UNKNOWN",
            "remedy": "자동 분류 실패. logs/ 의 원본 psi4 출력을 직접 확인할 것."}
