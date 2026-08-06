"""
13_nto_analysis.py
------------------
함수 벤치마크(12번)의 NTO 를 정량 분석한다: 각 함수가 추적한 UVA 전이가 정말
'같은 성격'(도너 고리 hole -> 케톤/억셉터 쪽 particle)의 CT 전이인지 확인.

왜 필요한가
  함수가 바뀌면 상태 순서가 바뀌거나 전이가 섞일 수 있어, 'S1 번호'로는 같은
  전이를 추적할 수 없다. 오비탈 번호(HOMO->LUMO)도 함수마다 오비탈 모양이 달라
  불충분하다. NTO 는 전이를 단일 hole->particle 쌍으로 압축하므로, hole/particle
  밀도의 '공간 분포'를 조각(fragment)별로 적분해 비교하면 함수 간 추적이 된다.
  이것이 CT 성격의 정량 확인이기도 하다 (지금까지는 '시사' 단계).

무엇을 계산하는가 (전이마다)
  - orca_plot 로 주 NTO 쌍(hole/particle)의 cube 생성
  - 조각별 밀도 분율: donor(NEt2 고리 쪽) / bridge(케톤 C=O) / acceptor(에스터 고리 쪽)
    조각은 구조 연결성에서 자동 결정: 케톤 탄소(=O 1개 + 고리 C 2개)를 그래프에서
    제거했을 때 N 을 포함한 성분 = donor, 반대쪽 = acceptor.
  - hole-electron 중심 거리 d_he (Angstrom)
  - overlap 지표 Lambda = integral sqrt(rho_h * rho_e) dV  (0=완전 CT, 1=국소 전이)

실행 (벤치마크 완료 후):
  .\scripts\run.ps1 -Molecule dhhb scripts\13_nto_analysis.py
  .\scripts\run.ps1 -Molecule dhhb scripts\13_nto_analysis.py --functionals b3lyp camb3lyp --modes tda
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qc_common import CALCULATIONS, load_checkpoint, read_xyz, save_checkpoint

OUTROOT = CALCULATIONS / "03_functional_benchmark"
BOHR_TO_A = 0.529177210903
GRID_INTERVALS = 60          # cube 격자 간격 수 (중심/조각적분/overlap 에 충분)

COVALENT_R = {"H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "S": 1.05}


# ------------------------------------------------------------------ 조각 결정
def bonds_from_geometry(symbols, coords) -> list[set]:
    """공유결합 반지름 합 x 1.25 이내면 결합으로 본다. 인접 리스트 반환."""
    n = len(symbols)
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            rmax = 1.25 * (COVALENT_R.get(symbols[i], 0.8)
                           + COVALENT_R.get(symbols[j], 0.8))
            if np.linalg.norm(coords[i] - coords[j]) <= rmax:
                adj[i].add(j)
                adj[j].add(i)
    return adj


def find_fragments(symbols, coords) -> dict:
    """
    케톤 다리를 기준으로 donor/bridge/acceptor 조각을 자동 결정한다.

    케톤 탄소 = O 이웃이 정확히 1개(그 O 는 다른 이웃 없음)이고 C 이웃이 2개인 탄소.
    (에스터 탄소는 O 이웃이 2개라 걸리지 않는다.)
    그래프에서 케톤 탄소를 제거하면 두 성분으로 갈라진다:
    N(디에틸아미노)을 포함한 쪽 = donor, 반대쪽 = acceptor.
    """
    adj = bonds_from_geometry(symbols, coords)
    ketone_c = None
    for i, s in enumerate(symbols):
        if s != "C":
            continue
        o_nb = [j for j in adj[i] if symbols[j] == "O"]
        c_nb = [j for j in adj[i] if symbols[j] == "C"]
        if len(o_nb) == 1 and len(c_nb) == 2 and len(adj[o_nb[0]]) == 1:
            ketone_c = i
            break
    if ketone_c is None:
        raise RuntimeError("케톤 탄소를 찾지 못했습니다 (이 조각 정의는 다이아릴케톤용).")
    ketone_o = next(j for j in adj[ketone_c] if symbols[j] == "O")

    # 케톤 탄소를 뺀 그래프에서 연결 성분 탐색
    comp = [-1] * len(symbols)
    cid = 0
    for start in range(len(symbols)):
        if comp[start] != -1 or start == ketone_c:
            continue
        stack = [start]
        comp[start] = cid
        while stack:
            u = stack.pop()
            for v in adj[u]:
                if v != ketone_c and comp[v] == -1:
                    comp[v] = cid
                    stack.append(v)
        cid += 1

    n_atoms = [i for i, s in enumerate(symbols) if s == "N"]
    if len(n_atoms) != 1:
        raise RuntimeError(f"N 원자가 1개가 아닙니다 ({len(n_atoms)}개).")
    donor_comp = comp[n_atoms[0]]

    frags = {"donor": [], "bridge": [ketone_c, ketone_o], "acceptor": []}
    for i in range(len(symbols)):
        if i in frags["bridge"]:
            continue
        (frags["donor"] if comp[i] == donor_comp else frags["acceptor"]).append(i)
    return frags


# ------------------------------------------------------------------ cube
def gen_cube(nto_file: Path, orbital: int) -> Path:
    """orca_plot 대화형 메뉴를 stdin 으로 구동해 MO cube 를 만든다."""
    expect = nto_file.parent / f"{nto_file.stem}.mo{orbital}a.cube"
    if expect.exists():
        return expect
    import os
    orca_dir = Path(os.environ.get("ORCA_EXE", r"C:\ORCA_6.1.1\orca.exe")).parent
    cmds = (f"1\n1\n"                  # plot type = MO
            f"4\n{GRID_INTERVALS}\n"   # 격자 간격
            f"5\n7\n"                  # 출력 형식 = cube
            f"2\n{orbital}\n11\n"      # 오비탈 선택 + 생성
            f"12\n")                   # 종료
    r = subprocess.run([str(orca_dir / "orca_plot.exe"), nto_file.name, "-i"],
                       cwd=str(nto_file.parent), input=cmds.encode(),
                       capture_output=True, check=False)
    if not expect.exists():
        raise RuntimeError(f"cube 생성 실패: {expect.name}\n"
                           + r.stdout.decode(errors="replace")[-800:])
    return expect


def read_cube(path: Path):
    """Gaussian cube -> (원점, 축벡터(3,3) [Bohr], 값 3D 배열)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        fh.readline(); fh.readline()
        parts = fh.readline().split()
        natoms = int(parts[0])
        origin = np.array([float(x) for x in parts[1:4]])
        nvec, vecs = [], []
        for _ in range(3):
            p = fh.readline().split()
            nvec.append(int(p[0]))
            vecs.append([float(x) for x in p[1:4]])
        for _ in range(abs(natoms)):
            fh.readline()
        vals = np.array(fh.read().split(), dtype=float)
    # MO cube 는 원자 블록 뒤에 '오비탈 개수 + 인덱스' 줄이 하나 더 붙는다
    # (Gaussian cube 규약: NAtoms 가 음수일 때). 그 줄의 길이가 가변이므로
    # 앞에서 세지 말고 뒤에서 격자 크기만큼 잘라낸다.
    n = nvec[0] * nvec[1] * nvec[2]
    return origin, np.array(vecs), vals[-n:].reshape(nvec)


def density_metrics(cube_h: Path, cube_p: Path, symbols, coords_A, frags) -> dict:
    """hole/particle cube -> 조각 분율, 중심 거리, overlap Lambda."""
    o1, v1, psi_h = read_cube(cube_h)
    o2, v2, psi_p = read_cube(cube_p)
    if psi_h.shape != psi_p.shape or not np.allclose(v1, v2):
        raise RuntimeError("hole/particle cube 격자가 다릅니다.")
    rho_h, rho_p = psi_h ** 2, psi_p ** 2
    dV = abs(np.linalg.det(v1))                      # Bohr^3

    nx, ny, nz = rho_h.shape
    ii, jj, kk = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz),
                             indexing="ij")
    pts = (o1[None, :] + ii.reshape(-1, 1) * v1[0] + jj.reshape(-1, 1) * v1[1]
           + kk.reshape(-1, 1) * v1[2]) * BOHR_TO_A   # (Npts,3) Angstrom

    wh = rho_h.reshape(-1) / rho_h.sum()
    wp = rho_p.reshape(-1) / rho_p.sum()
    cen_h = (pts * wh[:, None]).sum(axis=0)
    cen_p = (pts * wp[:, None]).sum(axis=0)
    d_he = float(np.linalg.norm(cen_p - cen_h))

    # overlap Lambda = int sqrt(rho_h rho_p) dV (각각 1로 규격화한 밀도 기준)
    nh = rho_h / (rho_h.sum() * dV)
    npd = rho_p / (rho_p.sum() * dV)
    lam = float(np.sqrt(nh * npd).sum() * dV)

    # 조각 분율: 격자점을 가장 가까운 원자에 배정 (Voronoi/Becke 근사)
    # 메모리 절약을 위해 블록 단위로 처리
    atom_pos = np.asarray(coords_A)
    frag_of_atom = np.empty(len(symbols), dtype=int)
    keys = list(frags)                                # ['donor','bridge','acceptor']
    for fi, k in enumerate(keys):
        frag_of_atom[frags[k]] = fi
    frac_h = np.zeros(len(keys))
    frac_p = np.zeros(len(keys))
    B = 200_000
    for s in range(0, pts.shape[0], B):
        blk = pts[s:s + B]
        d2 = ((blk[:, None, :] - atom_pos[None, :, :]) ** 2).sum(axis=2)
        nearest = frag_of_atom[np.argmin(d2, axis=1)]
        for fi in range(len(keys)):
            m = nearest == fi
            frac_h[fi] += wh[s:s + B][m].sum()
            frac_p[fi] += wp[s:s + B][m].sum()

    return {
        "hole_fraction": {k: round(float(frac_h[i]), 4) for i, k in enumerate(keys)},
        "particle_fraction": {k: round(float(frac_p[i]), 4) for i, k in enumerate(keys)},
        "hole_centroid_A": [round(float(x), 3) for x in cen_h],
        "particle_centroid_A": [round(float(x), 3) for x in cen_p],
        "d_hole_electron_A": round(d_he, 3),
        "overlap_lambda": round(lam, 4),
    }


# ------------------------------------------------------------------ 메인
def analyze_run(run_dir: Path, frags, symbols, coords_A) -> dict | None:
    rec = load_checkpoint(run_dir / "result.json")
    if not rec or not rec.get("ok"):
        return None
    uva = rec.get("uva_band_transition")
    if not uva:
        return None
    state = uva["state"]
    nto_pairs = (rec.get("nto") or {}).get(str(state)) or []
    if not nto_pairs:
        print(f"  [경고] {run_dir}: 상태 {state} 의 NTO 없음 "
              f"(NTOStates 범위 밖일 수 있음)")
        return None
    dom = nto_pairs[0]
    nto_file = run_dir / f"td.s{state}.nto"
    if not nto_file.exists():
        print(f"  [경고] {nto_file} 없음")
        return None

    cube_h = gen_cube(nto_file, dom["hole"])
    cube_p = gen_cube(nto_file, dom["particle"])
    met = density_metrics(cube_h, cube_p, symbols, coords_A, frags)

    # CT 판정(정량): hole 이 donor 에 주로 있고, particle 의 donor 분율이
    # hole 보다 뚜렷이 줄며 bridge+acceptor 로 이동했는가.
    ct_shift = round((met["hole_fraction"]["donor"]
                      - met["particle_fraction"]["donor"]), 4)
    out = {
        "functional_id": rec["functional_id"], "functional": rec["functional"],
        "hf_exchange": rec["hf_exchange"], "mode": rec["mode"],
        "state": state,
        "wavelength_nm": uva["wavelength_nm"], "osc_strength": uva["osc_strength"],
        "orbital_transitions_str": uva.get("orbital_transitions_str"),
        "dominant_nto_occ": dom["occ"],
        "n_significant_nto_pairs": sum(1 for p in nto_pairs if p["occ"] >= 0.10),
        **met,
        "donor_depletion": ct_shift,
    }
    print(f"  {rec['functional']:12s} {rec['mode']:5s} S{state} "
          f"{uva['wavelength_nm']:6.1f} nm | NTO점유 {dom['occ']:.3f} | "
          f"hole(donor {met['hole_fraction']['donor']:.2f}) -> "
          f"particle(donor {met['particle_fraction']['donor']:.2f}, "
          f"bridge {met['particle_fraction']['bridge']:.2f}, "
          f"acceptor {met['particle_fraction']['acceptor']:.2f}) | "
          f"d={met['d_hole_electron_A']:.2f} A, Lambda={met['overlap_lambda']:.3f}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--conf", default=None, help="기본값: 03_ 아래 유일한 컨포머")
    ap.add_argument("--functionals", nargs="*", default=None)
    ap.add_argument("--modes", nargs="*", default=None)
    args = ap.parse_args()

    confs = sorted(p for p in OUTROOT.iterdir() if p.is_dir()) \
        if OUTROOT.exists() else []
    if args.conf:
        confs = [OUTROOT / args.conf]
    if not confs:
        print("벤치마크 결과가 없습니다. 12_functional_benchmark.py 를 먼저 실행하세요.")
        return 1

    all_out = []
    for conf_dir in confs:
        # 구조는 첫 result.json 의 geometry_file 에서
        first = next(conf_dir.rglob("result.json"), None)
        rec0 = load_checkpoint(first) if first else None
        if not rec0:
            continue
        g = read_xyz(Path(rec0["geometry_file"]))
        frags = find_fragments(g.symbols, g.coords)
        print(f"\n=== {conf_dir.name}: 조각 크기 donor {len(frags['donor'])}, "
              f"bridge {len(frags['bridge'])}, acceptor {len(frags['acceptor'])} 원자 ===")
        for run_dir in sorted(p.parent for p in conf_dir.rglob("result.json")):
            fid = run_dir.parent.name
            mode = run_dir.name
            if args.functionals and fid not in args.functionals:
                continue
            if args.modes and mode not in args.modes:
                continue
            out = analyze_run(run_dir, frags, g.symbols, g.coords)
            if out:
                out["conf_id"] = conf_dir.name
                all_out.append(out)

    if all_out:
        save_checkpoint(OUTROOT / "nto_analysis.json", {"results": all_out})
        print(f"\n저장 -> {OUTROOT / 'nto_analysis.json'}")
        print("\n해석 기준: donor_depletion(hole~particle 의 donor 분율 감소)이 크고")
        print("Lambda 가 작을수록 CT 성격이 강하다. 함수 간에 이 패턴이 같은 상태를")
        print("추적해야 '같은 전이'를 비교하는 것이다.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
