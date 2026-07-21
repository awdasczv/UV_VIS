"""
qc_common.py
------------
프로젝트 전역에서 쓰는 공통 유틸.
  - 경로 상수
  - XYZ 읽기/쓰기
  - 볼츠만 가중치
  - 가우시안 broadening (stick -> 연속 스펙트럼)
  - 체크포인트 (JSON) 저장/로드  : 긴 계산이 중간에 죽어도 이어서 하기 위함
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
INPUTS = ROOT / "inputs"
STRUCTURES = ROOT / "structures"
CONFORMERS = ROOT / "conformers"
CALCULATIONS = ROOT / "calculations"
RESULTS = ROOT / "results"
LOGS = ROOT / "logs"

# 물리 상수
HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL = 627.5094740631
KB_KCAL = 1.987204259e-3      # kcal/(mol K)
EV_TO_NM = 1239.841984        # nm = EV_TO_NM / eV
T_DEFAULT = 298.15


# ------------------------------------------------------------------ XYZ I/O
@dataclass
class Geometry:
    symbols: list[str]
    coords: np.ndarray        # (N,3) angstrom
    comment: str = ""

    @property
    def natoms(self) -> int:
        return len(self.symbols)

    def to_xyz_block(self) -> str:
        """psi4 geometry 문자열 (단위 angstrom)"""
        return "\n".join(
            f"{s} {x:.10f} {y:.10f} {z:.10f}"
            for s, (x, y, z) in zip(self.symbols, self.coords)
        )

    def write(self, path: Path, comment: str | None = None) -> None:
        # 소수점 자리수는 to_xyz_block() 과 반드시 같아야 한다.
        # 예전에 여기만 8자리였는데, 앙상블 파일(10자리)에서 좌표를 읽어
        # 개별 파일로 다시 쓰면 반올림 경계에서 마지막 자리가 어긋났다.
        # (크기는 1e-8 A 로 물리적으로 무의미하지만, 개별 xyz 를 git 추적 대상에서
        #  빼고 앙상블에서 복원하는 방식을 쓰므로 왕복이 정확히 일치해야 한다.)
        c = self.comment if comment is None else comment
        lines = [str(self.natoms), c]
        lines += [f"{s:2s} {x:16.10f} {y:16.10f} {z:16.10f}"
                  for s, (x, y, z) in zip(self.symbols, self.coords)]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_xyz(path: Path) -> Geometry:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    n = int(lines[0].split()[0])
    comment = lines[1] if len(lines) > 1 else ""
    syms, xyz = [], []
    for ln in lines[2:2 + n]:
        p = ln.split()
        syms.append(p[0].capitalize())
        xyz.append([float(p[1]), float(p[2]), float(p[3])])
    return Geometry(syms, np.asarray(xyz, dtype=float), comment)


def read_multi_xyz(path: Path) -> list[Geometry]:
    """CREST/xtb 의 다중 구조 xyz (예: crest_conformers.xyz) 를 전부 읽는다."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    geoms, i = [], 0
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        n = int(lines[i].split()[0])
        comment = lines[i + 1].strip() if i + 1 < len(lines) else ""
        syms, xyz = [], []
        for ln in lines[i + 2:i + 2 + n]:
            p = ln.split()
            syms.append(p[0].capitalize())
            xyz.append([float(p[1]), float(p[2]), float(p[3])])
        geoms.append(Geometry(syms, np.asarray(xyz, float), comment))
        i += 2 + n
    return geoms


# -------------------------------------------------------- Boltzmann 가중치
def boltzmann_weights(energies_hartree, T: float = T_DEFAULT) -> np.ndarray:
    """절대에너지(hartree) 배열 -> 298.15 K 볼츠만 가중치 (합=1)"""
    e = np.asarray(energies_hartree, dtype=float)
    rel_kcal = (e - e.min()) * HARTREE_TO_KCAL
    w = np.exp(-rel_kcal / (KB_KCAL * T))
    return w / w.sum()


def rel_energies_kcal(energies_hartree) -> np.ndarray:
    e = np.asarray(energies_hartree, dtype=float)
    return (e - e.min()) * HARTREE_TO_KCAL


def select_by_cumulative_weight(weights, target: float = 0.90, max_n: int = 3):
    """
    가중치가 큰 순서로 누적합이 target 을 넘을 때까지 고른다.
    단, 최대 max_n 개까지만 (계산 비용 관리 - 요구사항 8번).
    반환: (선택된 인덱스 리스트, 실제 커버한 누적 가중치)
    """
    w = np.asarray(weights, dtype=float)
    order = np.argsort(-w)
    picked, cum = [], 0.0
    for idx in order:
        picked.append(int(idx))
        cum += float(w[idx])
        if cum >= target or len(picked) >= max_n:
            break
    return picked, cum


# ------------------------------------------------------------- 스펙트럼
def gaussian_spectrum(wavelengths_nm, osc_strengths, grid_nm,
                      fwhm_ev: float = 0.30, weight: float = 1.0):
    """
    stick spectrum -> 가우시안 broadening 된 몰흡광계수 epsilon(lambda).

    표준식 (Gaussian/ORCA 관례, 에너지축에서 대칭인 가우시안):
        eps(nu~) = 1.3062974e8 * sum_i  f_i / (sigma_cm) * exp(-((nu~ - nu~_i)/sigma_cm)^2)
    여기서 sigma 는 1/e 반폭. FWHM = 2*sqrt(ln2)*sigma.
    입력 선폭은 eV(FWHM)로 받고 내부에서 파수로 변환한다.
    단위: L mol^-1 cm^-1
    """
    grid_nm = np.asarray(grid_nm, dtype=float)
    lam = np.asarray(wavelengths_nm, dtype=float)
    f = np.asarray(osc_strengths, dtype=float)
    eps = np.zeros_like(grid_nm)
    if lam.size == 0:
        return eps

    nu_grid = 1.0e7 / grid_nm                  # cm^-1
    nu_i = 1.0e7 / lam                         # cm^-1
    fwhm_cm = fwhm_ev / HARTREE_TO_EV * 219474.6313632   # eV -> cm^-1
    sigma_cm = fwhm_cm / (2.0 * math.sqrt(math.log(2.0)))

    for nu0, fi in zip(nu_i, f):
        if fi <= 0:
            continue
        eps += 1.3062974e8 * (fi / sigma_cm) * np.exp(-(((nu_grid - nu0) / sigma_cm) ** 2))
    return eps * weight


# ---------------------------------------------------------- 체크포인트
def save_checkpoint(path: Path, data: dict) -> None:
    """계산 하나가 끝날 때마다 즉시 저장 (요구사항 17번)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=str),
                   encoding="utf-8")
    tmp.replace(path)


def load_checkpoint(path: Path) -> dict | None:
    p = Path(path)
    if not p.exists():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None
