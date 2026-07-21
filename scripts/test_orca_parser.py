"""ORCA 출력 파서를 실제 벤치마크 출력으로 검증한다."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from orca_common import find_orca, mpi_available, parse_output

ROOT = Path(__file__).resolve().parent.parent
out = ROOT / "calculations" / "00_orca_test" / "bench_enolA_par" / "bench.out"

print("orca.exe :", find_orca())
print("MPI 사용 가능 :", mpi_available())
print()

r = parse_output(out)
print(f"정상 종료 : {r['terminated_normally']}")
print(f"전자 수   : {r['n_electrons']}  (점유궤도 {r['n_occupied']})")
print(f"기저함수  : {r['n_basis']}")
print(f"SCF 에너지: {r['scf_energy_hartree']}")
print(f"ORCA 보고 실행시간: {r['orca_reported_runtime']}")
print(f"\n전이 {len(r['transitions'])} 개")
print(f"{'상태':>4s} {'eV':>8s} {'nm':>8s} {'f':>8s}  주요 오비탈 전이")
for t in r["transitions"]:
    bar = "*" * int(min(28, t["osc_strength"] * 28))
    print(f"{t['state']:4d} {t['energy_eV']:8.3f} {t['wavelength_nm']:8.1f} "
          f"{t['osc_strength']:8.4f}  {t['orbital_transitions_str'][:44]:44s} {bar}")

b = r.get("brightest")
if b:
    print(f"\n최강 전이 : {b['wavelength_nm']:.1f} nm (f={b['osc_strength']:.3f})")
    print(f"           {b['orbital_transitions_str']}")

print("\n--- Psi4 결과와 대조 (같은 구조/수준/기체상) ---")
psi4 = {1: (3.765, 329.3, 0.1380), 2: (3.946, 314.2, 0.8494),
        3: (4.422, 280.4, 0.2566), 12: (5.683, 218.2, 0.0185)}
print(f"{'상태':>4s} {'Psi4 eV':>9s} {'ORCA eV':>9s} {'차이':>8s} "
      f"{'Psi4 nm':>9s} {'ORCA nm':>9s}")
for t in r["transitions"]:
    if t["state"] in psi4:
        pe, pn, pf = psi4[t["state"]]
        print(f"{t['state']:4d} {pe:9.3f} {t['energy_eV']:9.3f} "
              f"{t['energy_eV']-pe:+8.3f} {pn:9.1f} {t['wavelength_nm']:9.1f}")
