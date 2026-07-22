"""
orca_common.py
--------------
ORCA 6.1 을 다루는 공통 모듈. 입력 생성 / 실행 / 출력 파싱 / 실패 분류.

왜 ORCA 로 갈아탔는가
  같은 계산(아보벤존 45원자, B3LYP/6-31G(d), TDA 12상태, 기체상)을
  두 프로그램으로 돌려 비교한 결과:

    Psi4 4코어 : 812 초
    ORCA 1코어 : 494 초
    ORCA 4코어 : 249 초   -> Psi4 대비 3.3배

  그러면서 결과는 사실상 동일했다 (최강 전이 314.2 nm vs 314.3 nm, 0.001 eV 차이).
  독립적인 두 코드가 일치하므로 양쪽 다 신뢰할 수 있고, 이제 더 빠른 쪽을 쓴다.
  덤으로 Psi4/PCMSolver 에서 나던 'S matrix is not positive-definite' 경고를
  ORCA 의 CPCM 으로 우회할 수 있다.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import time
from pathlib import Path

HARTREE_TO_EV = 27.211386245988
EV_TO_NM = 1239.841984


# ------------------------------------------------------------------ 실행파일
def find_orca() -> str:
    """orca.exe 경로를 찾는다."""
    env = os.environ.get("ORCA_EXE")
    if env and Path(env).exists():
        return env
    which = shutil.which("orca")
    if which:
        return which
    for c in (r"C:\ORCA_6.1.1\orca.exe", r"C:\ORCA_6.1.0\orca.exe",
              r"C:\Program Files\ORCA\orca.exe"):
        if Path(c).exists():
            return c
    raise FileNotFoundError(
        "orca.exe 를 찾을 수 없습니다. scripts/run.ps1 이 ORCA 설치 경로를 "
        "PATH 에 넣는지 확인하세요.")


def mpi_available() -> bool:
    if shutil.which("mpiexec"):
        return True
    return Path(r"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe").exists()


# ------------------------------------------------------------------ 입력 생성
def build_input(symbols, coords, *, functional: str, basis: str,
                nstates: int = 0, tda: bool = True, solvent: str | None = None,
                nprocs: int = 4, maxcore_mb: int = 2500,
                optimize: bool = False, rijcosx: bool = True,
                aux_basis: str = "def2/J", extra_keywords: str = "",
                comment: str = "") -> str:
    """
    ORCA 입력 파일 내용을 만든다.

    rijcosx=True 가 ORCA 속도의 핵심이다. Coulomb 항은 RI(보조기저 전개),
    정확교환 항은 COSX(수치격자)로 근사한다. 하이브리드 범함수의 가장 비싼
    부분을 줄여 주며, 이 크기 분자에서 오차는 무시할 수준이다.

    solvent 를 주면 CPCM 을 켠다. ORCA 는 TD-DFT 수직 전이에 대해
    비평형 용매화를 자동으로 처리한다.
    """
    kw = ["!", functional, basis]
    if rijcosx:
        kw += ["RIJCOSX", aux_basis]
    kw.append("TightSCF")
    if optimize:
        kw.append("Opt")
    if solvent:
        kw.append(f"CPCM({solvent})")
    if extra_keywords:
        kw.append(extra_keywords)

    lines = []
    if comment:
        for ln in comment.splitlines():
            lines.append(f"# {ln}")
    lines.append(" ".join(kw))
    if nprocs > 1:
        lines.append(f"%pal nprocs {nprocs} end")
    lines.append(f"%maxcore {maxcore_mb}")
    lines.append("%scf")
    lines.append("  MaxIter 300")
    lines.append("end")
    if nstates > 0:
        lines.append("%tddft")
        lines.append(f"  nroots {nstates}")
        lines.append(f"  tda    {'true' if tda else 'false'}")
        lines.append("  maxdim 5")
        lines.append("end")
    lines.append("* xyz 0 1")
    for s, (x, y, z) in zip(symbols, coords):
        lines.append(f"{s:2s} {x:16.10f} {y:16.10f} {z:16.10f}")
    lines.append("*")
    return "\n".join(lines) + "\n"


# ------------------------------------------------------------------ 실행
def run_orca(inp_text: str, workdir: Path, name: str = "job") -> tuple[Path, float]:
    """
    ORCA 를 돌리고 (출력파일 경로, 소요초) 를 돌려준다.

    주의할 점 세 가지
      1) 병렬 실행 시 ORCA 는 반드시 **전체 경로**로 호출해야 한다.
         (하위 모듈을 자기 위치 기준으로 찾기 때문)
      2) 셸(cmd /c)을 거치지 않는다.
         subprocess 에 리스트를 넘기면 파이썬이 Windows 명령줄로 합치면서
         내부 따옴표를 백슬래시로 이스케이프하는데, cmd 가 그것을 해석하지 못해
         "The directory name is invalid" 로 죽는다. 실제로 겪은 문제다.
         exe 를 직접 실행하고 stdout 을 파이썬이 파일로 받는다.
      3) 출력을 **바이너리로** 받아 파일에 쓴다. 한국어 Windows 의 cp949
         로케일에서 파이썬이 출력을 디코딩하다 깨지는 문제를 원천 차단한다.
         (같은 이유로 xtb 에서 stdout 이 None 이 되는 버그를 이미 겪었다)
    """
    workdir.mkdir(parents=True, exist_ok=True)
    inp = workdir / f"{name}.inp"
    out = workdir / f"{name}.out"
    inp.write_text(inp_text, encoding="ascii", errors="replace")

    exe = find_orca()
    t0 = time.time()
    with out.open("wb") as fh:
        subprocess.run([exe, inp.name], cwd=str(workdir),
                       stdout=fh, stderr=subprocess.STDOUT, check=False)
    return out, time.time() - t0


# ------------------------------------------------------------------ 파싱
_RE_NEL = re.compile(r"Number of Electrons\s+NEL\s*\.+\s*(\d+)")
_RE_NBF = re.compile(r"Number of basis functions\s*\.+\s*(\d+)")
_RE_ENERGY = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
_RE_STATE = re.compile(r"^STATE\s+(\d+):\s+E=\s+(-?[\d.]+)\s+au\s+([\d.]+)\s+eV")
_RE_CONTRIB = re.compile(r"^\s*(\d+)([ab])\s*->\s*(\d+)([ab])\s*:\s*([\d.]+)")
_RE_RUNTIME = re.compile(r"TOTAL RUN TIME:\s*(.+)")


def _orbital_label(idx: int, homo: int) -> str:
    """ORCA 의 0-기반 MO 번호를 HOMO/LUMO 기준 상대표기로 바꾼다."""
    if idx <= homo:
        d = homo - idx
        return "HOMO" if d == 0 else f"HOMO-{d}"
    d = idx - homo - 1
    return "LUMO" if d == 0 else f"LUMO+{d}"


def _flag_state(orbs: list[dict], nmo: int | None, homo: int | None) -> str | None:
    """
    들뜬상태의 전이 기여로부터 '의심스러운/인위적 해'를 판정한다.

    ORCA TDA 는 요청한 상태 수(nroots)가 물리적으로 의미있는 저에너지
    다양체보다 많으면, 계수가 정확히 +-1.000 이고 유일한 기여가 가장 높은
    가상궤도로만 가는 '가짜 뿌리(spurious root)'를 만들어낸다. 진동자 세기는
    0 이라 스펙트럼에는 무해하지만, 이런 전이는 물리적 의미가 없으므로 표시해 둔다.

    반환: 문제 없으면 None, 있으면 사유 문자열.
    """
    if not orbs or nmo is None or homo is None:
        return None
    top_virtual = nmo - 1                       # 가장 높은 가상궤도 인덱스
    lumo = homo + 1
    n_virtual = nmo - lumo
    # 유일한(또는 거의 유일한) 기여가 가상궤도 공간의 최상단 2% 안쪽으로 가고
    # 그 가중치가 0.98 이상이면 인위적 해로 본다.
    dominant = max(orbs, key=lambda o: o["weight"])
    if dominant["weight"] >= 0.98:
        to = dominant["to"]
        if to.startswith("LUMO+"):
            off = int(to.split("+")[1])
            if lumo + off >= top_virtual - max(2, int(0.02 * n_virtual)):
                return (f"의심: 계수≈1 이 최상단 가상궤도({to})로만 감. "
                        "ORCA TDA 의 인위적 뿌리로 추정 (f≈0, 스펙트럼 무해).")
    return None


def parse_output(out_path: Path) -> dict:
    """ORCA 출력에서 필요한 것들을 뽑아낸다."""
    text = out_path.read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()

    res: dict = {"terminated_normally": "****ORCA TERMINATED NORMALLY****" in text}

    m = _RE_NEL.search(text)
    nel = int(m.group(1)) if m else None
    nocc = nel // 2 if nel else None
    homo = nocc - 1 if nocc else None          # ORCA MO 번호는 0부터 시작
    res["n_electrons"] = nel
    res["n_occupied"] = nocc

    m = _RE_NBF.search(text)
    res["n_basis"] = int(m.group(1)) if m else None
    nmo = res["n_basis"]                        # MO 개수 = 기저함수 개수
    n_bad_index = 0                             # 유효범위를 벗어난 오비탈 인덱스 개수

    ms = _RE_ENERGY.findall(text)
    res["scf_energy_hartree"] = float(ms[-1]) if ms else None

    m = _RE_RUNTIME.search(text)
    res["orca_reported_runtime"] = m.group(1).strip() if m else None

    # --- 들뜬상태별 주요 오비탈 전이 ---
    contribs: dict[int, list[dict]] = {}
    cur = None
    for ln in lines:
        sm = _RE_STATE.match(ln.strip())
        if sm:
            cur = int(sm.group(1))
            contribs[cur] = []
            continue
        if cur is not None:
            cm = _RE_CONTRIB.match(ln)
            if cm:
                i, _, a, _, w = cm.groups()
                occ_idx, vir_idx = int(i), int(a)
                # 인덱스 유효범위 검사: 0 <= occ <= HOMO,  LUMO <= vir <= 최고 MO.
                # 벗어나면 파서가 엉뚱한 숫자를 읽은 것이므로 버리고 센다.
                if homo is not None and nmo is not None:
                    if not (0 <= occ_idx <= homo and homo < vir_idx < nmo):
                        n_bad_index += 1
                        continue
                if homo is not None:
                    contribs[cur].append({
                        "from": _orbital_label(occ_idx, homo),
                        "to": _orbital_label(vir_idx, homo),
                        "weight": round(float(w), 4),
                    })
            elif ln.strip() == "" and contribs.get(cur):
                cur = None
    res["n_out_of_range_orbital_indices"] = n_bad_index

    # --- 흡수 스펙트럼 ---
    # ORCA 는 같은 형식의 표를 여러 번 찍는다 (전기 쌍극자 / 속도 쌍극자 /
    # CD 회전세기 ...). 우리가 쓸 것은 **전기 쌍극자 표 하나뿐**이므로
    # 그 블록만 읽고 표가 끝나면 즉시 멈춘다.
    trans = []
    idx = None
    for k, ln in enumerate(lines):
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in ln:
            idx = k
            break
    if idx is not None:
        started = False
        for ln in lines[idx + 1:idx + 400]:
            # 예: "  0-1A  ->  1-1A    3.759416   30321.7   329.8   0.129186760 ..."
            m = re.match(r"\s*\S+\s*->\s*(\d+)-\S+\s+([\d.]+)\s+([\d.]+)\s+"
                         r"([\d.]+)\s+([\d.eE+-]+)", ln)
            if m:
                started = True
                state = int(m.group(1))
                all_orbs = sorted(contribs.get(state, []),
                                  key=lambda d: -d["weight"])
                orbs = all_orbs[:3]
                flag = _flag_state(all_orbs, nmo, homo)
                entry = {
                    "state": state,
                    "energy_eV": round(float(m.group(2)), 5),
                    "wavelength_nm": round(float(m.group(4)), 3),
                    "osc_strength": round(float(m.group(5)), 6),
                    "orbital_transitions": orbs,
                    "orbital_transitions_str": "; ".join(
                        f"{o['from']}->{o['to']} ({o['weight']*100:.0f}%)" for o in orbs),
                }
                if flag:
                    entry["flag"] = flag
                trans.append(entry)
            elif started:
                # 표가 끝났다 (빈 줄 또는 구분선). 다음 블록은 읽지 않는다.
                break
    res["transitions"] = trans
    if trans:
        res["brightest"] = max(trans, key=lambda t: t["osc_strength"])
    return res


# ------------------------------------------------------------------ 실패 분류
ORCA_FAILURE_PATTERNS = [
    (r"SCF NOT CONVERGED AFTER|This wavefunction IS NOT FULLY CONVERGED",
     "SCF_NOT_CONVERGED",
     "SCF 미수렴. 대응: %scf MaxIter 늘리기, SOSCF 켜기, "
     "'! SlowConv' 또는 'VerySlowConv' 추가, Damp/Shift 조정, "
     "더 작은 기저셋의 궤도를 MORead 로 초기 guess 로 재사용."),
    (r"CIS/TDA.*not converged|Davidson.*not converged|TDDFT.*NOT CONVERGED",
     "EXCITED_STATE_NOT_CONVERGED",
     "들뜬상태 미수렴. 대응: %tddft 의 maxdim 키우기, MaxIter 늘리기, "
     "nroots 줄이기, TDA 사용."),
    (r"not enough memory|insufficient memory|MEMORY|out of memory|"
     r"Error.*allocat",
     "OUT_OF_MEMORY",
     "메모리 부족. 대응: %maxcore 낮추기(코어당 MB), nprocs 줄이기, "
     "기저셋/상태 수 줄이기."),
    (r"mpiexec.*not recognized|MPI.*fail|error.*mpi",
     "MPI_ERROR",
     "MPI 오류. 대응: Microsoft MPI 설치 확인, PATH 에 mpiexec 존재 확인, "
     "임시로 nprocs 1 로 직렬 실행."),
    (r"Error.*basis|BASIS.*not.*found|unknown basis",
     "BASIS_NOT_FOUND",
     "기저셋 이름 오류. ORCA 매뉴얼의 표기를 확인할 것."),
    (r"atoms.*too close|SEVERE.*geometry|coordinates",
     "BAD_GEOMETRY",
     "구조 이상. 입력 xyz 재검증, xTB 로 사전 최적화."),
    (r"ORCA finished by error termination in\s+(\w+)",
     "MODULE_ERROR",
     "특정 모듈에서 종료. 출력 파일에서 해당 모듈 직전 메시지를 확인할 것."),
]


def classify_orca_failure(out_path: Path | None, text: str | None = None) -> dict:
    if text is None:
        text = out_path.read_text(encoding="utf-8", errors="replace") if out_path \
            and out_path.exists() else ""
    for pattern, code, remedy in ORCA_FAILURE_PATTERNS:
        if re.search(pattern, text, re.IGNORECASE):
            return {"code": code, "remedy": remedy}
    return {"code": "UNKNOWN",
            "remedy": "자동 분류 실패. ORCA 출력 파일 끝부분을 직접 확인할 것."}
