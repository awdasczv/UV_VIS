# 아보벤존(Avobenzone) 이론 UV–Vis 흡광 스펙트럼 계산

로컬 PC에서 아보벤존(butyl methoxydibenzoylmethane, CAS 70356-09-1)의
**킬레이트 에놀형 두 가지**와 **디케토형**의 UV–Vis 흡수 스펙트럼을
DFT / TD-DFT 수준으로 계산하고, 알려진 실험 최대흡수파장과 비교한다.

## 최종 결과 (요약)

최고 수준 **B3LYP/6-31+G(d) + DFT 최적화 구조 + CPCM(에탄올), TDA** (ORCA 6.1.1):

| 토토머 | 계산 λmax | 실험 λmax | 오차 | 주요 전이 |
|---|---|---|---|---|
| 킬레이트 에놀 A | 352.5 nm | 354.9 nm | **−2.4 nm** | HOMO→LUMO 98% (π→π*) |
| 킬레이트 에놀 B | 356.1 nm | 354.9 nm | **+1.2 nm** | HOMO→LUMO 96% |
| 디케토 | 268.7 nm | 265.0 nm | **+3.7 nm** | HOMO→LUMO+1 72% |

에놀 A 의 오차 39.7 nm 를 세 요인으로 분해했다:
**용매 +18.8 nm · 구조(xTB→DFT) +10.9 nm · 기저셋 +7.6 nm.**

저비용 기준선(MINDO/3 대체, HF/CIS)은 강한 밴드를 −107~−152 nm 로 완전히
빗나갔다. 본 계산이 이를 한 자릿수 nm 로 줄였다. 산출물은 [results/](results/) 참고
(최종 보고서: [results/report.md](results/report.md), 비교 그래프: [results/comparison.png](results/comparison.png)).

---

## 1. 계산 목적

기존의 저비용 반경험적 계산(MINDO/3–TDA 수준)보다 **한 단계 높은 이론 수준**에서

1. 아보벤존의 세 가지 토토머(에놀 A / 에놀 B / 디케토)를 각각 계산하고
2. 컨포머 앙상블을 볼츠만 평균하여
3. 에탄올 연속용매(PCM) 조건에서 TD-DFT 흡수 스펙트럼을 얻은 뒤
4. 실험값 **에놀 밴드 354.9 nm (에탄올)**, **디케토 밴드 265–269 nm** 과 비교한다.

특히 다음 네 가지가 λmax에 각각 얼마나 기여하는지 분리해서 본다.

| 요인 | 확인 방법 |
|---|---|
| 토토머 선택 | 에놀 A / 에놀 B / 디케토 를 따로 계산 |
| 컨포머 평균 | 최저에너지 단일 구조 vs 볼츠만 가중 앙상블 |
| 용매 모델 | 기체상 vs PCM(에탄올) |
| 함수·기저셋 | B3LYP / PBE0 / CAM-B3LYP / ωB97X-D × 여러 기저셋 |

---

## 2. 이 PC의 환경 (2026-07-20 점검)

| 항목 | 값 |
|---|---|
| OS | Windows 11 Pro 10.0.26200 |
| CPU | Intel Core i5-1135G7 — 물리 4코어 / 8스레드 |
| RAM | 15.7 GB |
| GPU | Intel Iris Xe (내장) — CUDA 없음, GPU 가속 불가 |
| Python | 3.11.9 (시스템), 3.10 (py 런처) |
| Conda | 없음 → 프로젝트 안에 micromamba 로 별도 구축 |
| WSL | WSL2 동작하나 리눅스 배포판은 없음 (docker-desktop 만) |
| 기존 양자화학 SW | **없음** (ORCA / Gaussian / xtb / CREST 모두 미설치) |

→ **GPU 가속 불가, 물리 4코어**. 이것이 이론 수준을 고르는 가장 큰 제약이다.

---

## 3. 계산 엔진 선정 근거

### 3.1 후보 조사 결과 (플랫폼 실제 확인)

| 도구 | Windows 네이티브 | 근거 |
|---|---|---|
| **Psi4** | ✅ win-64 | conda-forge 에 `win-64` 빌드 존재 (psi4 1.11) |
| **xTB** | ✅ win-64 | conda-forge `win-64`, GitHub 릴리스에 windows zip |
| **PCMSolver** | ✅ win-64 | psi4 win-64 패키지가 `pcmsolver 1.2.3` 을 의존성으로 포함 |
| PySCF | ❌ | 공식 문서: *"PySCF is not supported natively on Windows. You must use WSL."* |
| CREST | ❌ | conda-forge 에 linux/osx 만. GitHub 릴리스도 ubuntu 바이너리만 |
| DDX (ddCOSMO) | ❌ | conda-forge `pyddx` 에 win-64 빌드 없음 |
| ORCA | ❌ 미설치 | 학술용 무료지만 회원가입 후 **수동 다운로드** 필요 |

### 3.2 프로그램 출력으로 직접 확인한 것 (`scripts/00_probe_engine.py`)

`logs/engine_probe.json` 에 전체 결과가 있다.

* **사용 가능한 함수**: B3LYP, PBE0, CAM-B3LYP, ωB97X-D, ωB97X-V, ωB97M-V,
  M06-2X, BHandHLYP, LRC-wPBEh, HSE06 — 모두 실제 에너지 계산 성공.
  * `-D3BJ` 변종은 `s-dftd3.exe` 가 PATH에 없어 실패했었다 → `scripts/run.ps1` 이
    환경의 `Library\bin` 을 PATH에 넣어 해결.
  * ⚠️ ωB97X-V / ωB97M-V 는 VV10 성분 때문에 **Psi4 TDSCF 에서 사용 불가**(문서 명시 제한).
* **아보벤존(45원자) 기저함수 개수** — 비용의 직접 지표:

  | 기저셋 | 기저함수 수 |
  |---|---|
  | 6-31G(d) | 389 |
  | 6-31+G(d) | 481 |
  | def2-SVP | 432 |
  | def2-SVPD | 645 |
  | def2-TZVP | 845 |
  | def2-TZVPD | 1058 |

* **TD-DFT + PCM 조합이 실제로 작동함** ← 이 프로젝트의 핵심 확인 사항.
  Psi4 공식 문서는 PCM+TD-SCF를 *"RHF/UHF 에 대해 지원"* 이라고만 적고 있어
  TD-**DFT** + PCM 은 문서상 미검증 영역이다. 직접 시험한 결과:

  | 계산 | S1, S2, S3 (eV) |
  |---|---|
  | 기체상 TDA (B3LYP/STO-3G, H₂O) | 11.714, 13.976, 14.648 |
  | **PCM TDA (동일 조건)** | **11.982, 14.167, 14.820** |

  값이 실제로 달라지므로 PCM 이 TD-DFT 에 반영되고 있다.
  → **WSL 없이 Windows 네이티브로 요구사항을 모두 만족할 수 있다.**

### 3.3 최종 선택 — 계산 엔진을 ORCA 로 전환

프로젝트 도중 **ORCA 6.1.1** 을 설치해 Psi4 와 동일 조건으로 벤치마크했다.

| 엔진 | 코어 | 시간 (enolA, B3LYP/6-31G(d), TDA 12상태, 기체상) |
|---|---|---|
| Psi4 1.11 | 4 | 812 초 |
| ORCA 6.1.1 | 1 | 494 초 |
| **ORCA 6.1.1** | **4** | **249 초** (Psi4 대비 3.3배) |

그러면서 결과는 사실상 동일했다 (최강 전이 314.2 vs 314.3 nm, 0.001 eV 차이,
오비탈 귀속 HOMO→LUMO 74% vs 75%). **서로 독립적으로 개발된 두 코드가
일치하므로 양쪽 결과를 모두 신뢰할 수 있다.** 이후 본 계산은 더 빠른 ORCA 로 했다.
(ORCA 의 속도 이점은 RIJCOSX 근사 덕분. Psi4/PCMSolver 의 수치 경고도 ORCA CPCM 으로 우회.)

**최종 파이프라인:**

```
구조 생성      RDKit ETKDGv3 (200개) + MMFF94s  (킬레이트 이면각 명시적 세팅)
컨포머 탐색    GFN2-xTB (ALPB ethanol) — 구조당 약 1.8 초
              → 대칭 인식 RMSD 중복제거 → 발색단 기하 클러스터링 대표 선택
구조 최적화    ORCA B3LYP/def2-SVP + CPCM(ethanol)
들뜬상태       ORCA TD-DFT (TDA) + CPCM(ethanol, 비평형), RIJCOSX
              최고 수준: B3LYP/6-31+G(d), 18~22 상태
```

CREST 는 Windows 에서 못 쓰므로 **RDKit ETKDG(다수 구조) → GFN2-xTB 최적화 →
대칭 인식 RMSD 중복제거** 로 대체한다. xTB 가 구조당 2초 수준이라 수백 개를
돌려도 몇 분이면 끝난다.

> **주의**: GFN2-xTB 구조는 아릴 고리를 약 27° 비틀어 π 공액을 약화시켜 λmax 를
> 약 11 nm 단파장으로 민다. 따라서 xTB 는 **컨포머 탐색용**이고, 최종 λmax 는
> 반드시 **DFT 최적화 구조**로 얻어야 한다. (에놀에서 특히 크고, 디케토는 두
> 발색단이 이미 분리되어 영향이 작다 — `results/geometry_analysis.csv`)

---

## 4. 실측한 계산 비용

이 PC(4코어)에서 **아보벤존 45원자**를 실제로 돌려서 잰 값이다.
(`calculations/00_test/test_report.json`, `calculations/00_test/pcm_cost_probe.json`)

| 단계 | 조건 | 실측 시간 |
|---|---|---|
| GFN2-xTB 구조 최적화 | ALPB(에탄올), 41 스텝 수렴 | **1.8 초** |
| DFT SCF | B3LYP/6-31G (251 기저함수), 기체상 | 76 초 |
| TD-DFT (TDA, 5상태) | 위와 동일, 기체상 | 149 초 (Davidson 1회 ≈ 25초) |
| TD-DFT (TDA, 5상태) | 위와 동일, **PCM(에탄올), Area 0.3 Å²** | **Davidson 1회 ≈ 440초** |

### 여기서 나온 중요한 사실

1. **컨포머 탐색은 사실상 공짜다.** xTB 가 구조당 2초이므로 수백 개를 돌려도 몇 분이다.
   CREST 가 없어도 요구사항 7·8번을 충분히 만족할 수 있다.
2. **PCM 이 TD-DFT 를 약 17배 느리게 만든다.** 이것이 이 프로젝트의 진짜 병목이다.
   원인은 Davidson 시행벡터마다 공동(cavity) 표면의 겉보기 전하를 다시 풀어야 하기 때문이고,
   비용은 tessera(표면 조각) 개수에 좌우된다. `scripts/03b_pcm_cost_probe.py` 가
   `Cavity Area` 를 바꿔가며 비용과 정확도의 절충을 정량화한다.
3. 따라서 **무작정 큰 기저셋으로 가는 것보다, PCM 설정과 계산 엔진을 먼저 손보는 것이
   훨씬 효과가 크다.**

### 4c. PCM 비용의 실체

처음에는 PCM 이 TD-DFT 를 17배 느리게 만드는 것처럼 보였다. 그러나 그것은
**첫 Davidson 반복에 포함된 일회성 초기화 비용을 반복당 비용으로 잘못 외삽한 것**이었다.
총 소요시간으로 다시 재면:

| 조건 | 표면 조각 수 | SCF | TD-DFT (5상태) | λmax | f |
|---|---|---|---|---|---|
| 기체상 | — | 76 초 | 149 초 | 315.7 nm | 1.012 |
| PCM(에탄올), Area 1.0 Å² | 1024 | 325 초 | **566 초** | **332.3 nm** | **1.141** |

**PCM 은 약 3.8배 비쌀 뿐이며, 충분히 감당 가능하다.**
(첫 반복이 유난히 오래 걸리는 것은 공동 표면 관련 자료구조를 만드는 비용 때문이고,
이후 반복은 훨씬 빠르다.)

그리고 결과가 문헌과 잘 맞는다:

| | 본 계산 (B3LYP/6-31G) | 문헌 (B3LYP/6-31+G(d), ACS Omega 2026) |
|---|---|---|
| 용매에 의한 적색이동 | **+16.6 nm** | +18.0 nm |
| 진동자 세기 변화 | 1.012 → 1.141 | 0.991 → 1.138 |

전이 성격도 HOMO→LUMO 94% 의 깨끗한 π→π* 로, 킬레이트 에놀의 UVA 밴드와 일치한다.

**교훈**: "느리다"는 인상만으로 계산 수준을 낮추지 않고 총 비용을 실제로 재본 것이
계획을 크게 바꿨다. 이 프로젝트에서 실패·병목은 항상 먼저 측정하고 분류한다.

## 4b. 계산 수준 (확정)

| 단계 | 수준 |
|---|---|
| 컨포머 탐색 | RDKit ETKDGv3 (200 구조) → GFN2-xTB / ALPB(에탄올) → 대칭 인식 RMSD 중복제거 → 발색단 기하 클러스터링 |
| 구조 최적화 | ORCA **B3LYP/def2-SVP + CPCM(에탄올)** — 토토머별 최저 컨포머 |
| TD-DFT (본 계산) | ORCA **B3LYP/6-31+G(d)**, TDA, 18~22 상태, CPCM(에탄올, 비평형), RIJCOSX |
| 함수 비교 | **CAM-B3LYP/def2-SVP** (B3LYP 와 대조) |
| 기저셋 수렴 | def2-SVP / 6-31+G(d) / def2-SVPD / def2-TZVP |
| 저비용 기준선 | **HF/STO-3G + CIS**, **HF/def2-SVP + CIS** (MINDO/3 대체) |

**기저셋 근거**: 실측한 기저함수 개수(6-31G(d)=389, 6-31+G(d)=481, def2-TZVP=845)와
수렴 시험 결과, def2-SVP(432)는 너무 작아 λmax 를 6 nm 과소평가한다. 큰 기저셋
셋(6-31+G(d), def2-SVPD, def2-TZVP)이 모두 약 340.5 nm 로 수렴하므로
**6-31+G(d)** 를 본 계산 기저셋으로 채택했다(수렴값에 도달하면서 def2-TZVP 보다 저렴).

**함수 근거**: B3LYP 는 실험과 −2.4 nm, CAM-B3LYP 는 −34 nm(과도한 청색이동).
함수 선택이 이 발색단의 결과를 지배하므로 두 함수를 모두 계산해 명시했다.

### 더 높은 수준이 필요할 때 (데스크탑)

`inputs/calc_config.json` 의 `active_profile` 을 `desktop_5600x` 로 바꾸면
def2-SVP + 4개 함수(B3LYP / PBE0 / CAM-B3LYP / ωB97X-D)로 확장된다.

**하드웨어 조사 결론** (자세한 근거는 10장):

* Ryzen 5 5600X 는 이 노트북 대비 **약 3배** (PassMark 21,828 vs 9,317, 게다가
  노트북은 장시간 계산에서 12–15 W 로 열제한이 걸린다).
* RTX 5060 Ti 는 **기대만큼 도움이 되지 않는다.** 양자화학은 전부 배정밀도(FP64)로
  도는데, 이 카드의 FP64 성능은 약 370 GFLOP/s 로 FP32의 1/64 이며
  **5600X CPU(약 375 GFLOP/s)와 사실상 동급**이다.
  GPU4PySCF 가 A100 에서 보인 20배 같은 수치는 이 카드에 적용되지 않는다.
  RTX 50 시리즈에서의 벤치마크는 존재하지 않고, Blackwell 커널 최적화가
  미완성이라는 미해결 이슈도 있다. 현실적 기대치는 3–8배(미검증).
* **가장 효과적인 다음 수단은 ORCA 6.1 이다.** ORCA 는 GPU 를 전혀 쓰지 않지만
  (6.1 매뉴얼 전체에 "GPU"·"NVIDIA" 가 0회 등장) RIJCOSX 근사 덕분에
  같은 CPU 에서 Psi4 보다 훨씬 빠르고, CPCM/SMD 를 포함한 TD-DFT 구현이 성숙해 있다.
  Windows 네이티브 설치 파일이 있어 WSL 도 CUDA 도 필요 없다.
  **단, orcaforum 가입 후 수동 다운로드가 필요하다 (자동화 불가, 사용자 조치 사항).**

---

## 5. 프로젝트 구조

```
inputs/          토토머 정의(SMILES·근거), 실험 참조값, 계산 설정
structures/      토토머별 초기 3D 구조 + 검증 리포트
conformers/      컨포머 앙상블, 에너지·볼츠만 가중치, 대표 컨포머 선택
calculations/    DFT 최적화 / TD-DFT 원본 결과와 체크포인트
results/         CSV, PNG, 최종 분석 보고서
scripts/         모든 계산 스크립트 (재현용)
logs/            프로그램 원본 로그, 실패 분석
tools/           micromamba (git 추적 제외)
.mamba/          conda 환경 (git 추적 제외)
```

---

## 6. 설치

시스템 전역 설정은 건드리지 않는다. 모든 것이 프로젝트 폴더 안에 들어간다.

```powershell
.\scripts\setup_env.ps1
```

이 스크립트는

1. `tools/micromamba.exe` 를 내려받고 (단일 실행파일, 사용자 영역)
2. `.mamba/envs/qc` 에 conda-forge 로 환경을 만든다:
   `python=3.11 psi4 xtb rdkit numpy scipy pandas matplotlib`

**라이선스 / 수동 설치 구분**

| 항목 | 구분 |
|---|---|
| Psi4, xTB, RDKit, PCMSolver | 자동 설치, 오픈소스 (LGPL/GPL/BSD) — 별도 라이선스 절차 없음 |
| CREST | Windows 바이너리 없음 → WSL 필요 (**사용자 확인 후에만 진행**) |
| ORCA | 학술용 무료지만 orcaforum 가입 후 **수동 다운로드** — 자동화 불가 |
| Gaussian | 상용 라이선스 — 이 프로젝트에서는 사용하지 않음 |

---

## 7. 실행 방법

모든 스크립트는 `scripts/run.ps1` 래퍼로 실행한다. 이 래퍼가 micromamba 환경과
ORCA/MS-MPI 경로를 PATH 에 넣고, 스레드 수와 scratch 경로를 고정한다.

```powershell
# 1) 세 토토머의 3D 구조 생성 + 연결·수소·킬레이트 검증
.\scripts\run.ps1 scripts\01_build_structures.py

# 2) 컨포머 탐색 → 검증 → 발색단 클러스터링 대표 선택
.\scripts\run.ps1 scripts\02_conformer_search.py --backend etkdg --nconf 200
.\scripts\run.ps1 scripts\02b_validate_conformers.py
.\scripts\run.ps1 scripts\02c_select_representatives.py --intra-cluster-check

# 3) DFT 구조 최적화 (ORCA)
.\scripts\run.ps1 scripts\04b_dft_optimize_orca.py

# 4) TD-DFT 본 계산 (ORCA) — DFT 구조로 최고 수준
.\scripts\run.ps1 scripts\05b_tddft_orca.py --levels b3lyp_631+gd --geom-label dftopt22

# 5) 저비용 기준선 (MINDO/3 대체)
.\scripts\run.ps1 scripts\08b_lowcost_baseline_orca.py

# 6) 스펙트럼 + CSV + 그래프, 그리고 최종 보고서
.\scripts\run.ps1 scripts\06_build_spectra.py
.\scripts\run.ps1 scripts\07_report.py

# 진행 상황을 한눈에
.\scripts\run.ps1 scripts\show_status.py
```

모든 계산은 `(구조 × 이론수준 × 용매)` 조합 단위로 `result.json` 체크포인트를
남기고, 재실행 시 완료된 조합을 자동으로 건너뛴다. 중간에 꺼도 안전하다.

> Psi4 판 스크립트(`00_probe_engine.py`, `03_*`, `05_tddft.py`)도 저장소에 남아 있다.
> ORCA 로 전환하기 전의 엔진 능력 점검·비용 측정·교차검증에 쓰였다.

---

## 8. 선택 사항: 더 높은 수준으로 확장하는 세 가지 길

현재 파이프라인은 이것들 없이도 완결된다. 아래는 모두 **사용자 확인이 필요한**
시스템 수준 변경이거나 수동 설치다.

### 8.1 ORCA 6.1 (가장 권장)

TD-DFT + 용매를 이 크기 분자에서 실용적인 시간에 돌리는 가장 확실한 방법.

* 학술 목적 무료. `orcaforum.kofo.mpg.de` **가입 후 수동 다운로드** (자동화 불가).
* Windows 네이티브 설치 파일 존재 → WSL·CUDA 불필요.
* RIJCOSX 근사 + 성숙한 CPCM/SMD 구현.
* GPU 가속은 **전혀 없다** (6.1 매뉴얼에 "GPU"·"NVIDIA" 0회 등장).

### 8.2 WSL 리눅스 (CREST / PySCF / GPU4PySCF)

현재 이 PC 의 WSL 에는 `docker-desktop` 하나만 등록되어 있다. 이것은 Docker Desktop
전용 최소 배포판(90 MB)이라 일반 리눅스로 쓸 수 없다. Ubuntu 를 새로 설치해야 한다.

```powershell
wsl --install -d Ubuntu-24.04        # 약 2-3 GB, 재부팅이 필요할 수 있음
# 이후 WSL 안에서
#   micromamba create -n qc -c conda-forge crest xtb pyscf
```

**다만 Docker 가 이미 설치되어 동작 중이므로**(8 CPU / 8 GB 할당),
배포판을 새로 등록하지 않고 컨테이너로 리눅스 도구를 쓰는 방법도 있다.

### 8.3 GPU (GPU4PySCF) — 기대치를 낮출 것

RTX 5060 Ti 에서 실행하려면 WSL2 + CUDA 13 + `gpu4pyscf-cuda13x` 휠이 필요하다
(`cuda12x` 휠에는 sm_120 네이티브 코드가 없어 PTX JIT 로만 돌아간다).
TD-DFT + IEF-PCM 은 v1.4.1 부터 지원된다.

그러나 **소비자용 Blackwell 의 FP64 는 FP32 의 1/64 로 제한**되어 있어
이 카드의 배정밀도 성능(약 370 GFLOP/s)이 5600X CPU(약 375 GFLOP/s)와 사실상 같다.
RTX 50 시리즈에서 이 코드를 측정한 벤치마크는 존재하지 않으며,
Blackwell 커널 성능 이슈가 미해결로 남아 있다. 3–8배 정도가 현실적 기대치다(미검증).
또한 def2-TZVP 급에서는 DF 텐서가 약 14 GB 라 16 GB VRAM 에 다 들어가지 않아
시스템 RAM 을 쓰게 되므로 **32 GB 이상 RAM 을 권장**한다.

---

## 8b. 저장소에 무엇을 넣고 무엇을 뺐는가

계산은 용량이 큰 중간 산출물을 대량으로 만든다. 원칙은
**"계산을 재현하고 검증하는 데 필요한 것만 추적한다"** 이다.

### calculations/ 는 화이트리스트 방식

ORCA 는 계산 1건마다 `.gbw`(5.7 MB × 2), `.cis`(2.9 MB), `.densities`(2.9 MB),
`.diisao.tmp.0`(40 MB) 같은 파일을 쏟아낸다. 확장자를 하나씩 막는 방식은
계속 새는 것이 생기므로(실제로 `*.tmp` 패턴이 `.tmp.0` 을 못 걸러 40 MB 가
추적될 뻔했다) **전부 무시하고 필요한 것만 되살린다.**

| 확장자 | 추적 | 이유 |
|---|---|---|
| `.inp` | ✅ | 계산 입력. 이것만 있으면 그대로 재현 가능 |
| `.out` | ✅ | 사람이 읽는 전체 로그. 검증의 근거 |
| `.json` | ✅ | 파싱된 결과와 체크포인트 |
| `.xyz` | ✅ | 구조 |
| 나머지 전부 | ❌ | 파동함수·밀도·임시파일. 재계산으로 얻을 수 있다 |

### 개별 컨포머 xyz 는 추적하지 않는다

같은 좌표가 `conformers/<토토머>/conformers.xyz` 안에 전부 들어 있어 중복이다
(파일 46개). 저장소를 새로 clone 한 뒤에는 다음 한 줄로 복원한다.

```powershell
.\scripts\run.ps1 scripts\02d_extract_conformers.py
```

복원 결과가 원본과 **바이트 단위로 일치**함을 확인했다.
(처음에는 앙상블 파일이 소수점 10자리, 개별 파일이 8자리라 반올림 경계에서
마지막 자리가 어긋났다. 크기는 10⁻⁸ Å 로 물리적으로 무의미하지만,
복원 방식을 쓰는 이상 정확히 같아야 하므로 자리수를 맞췄다.)

### 결과

| | 파일 수 | 용량 |
|---|---|---|
| 정리 전 | 5,535 | 225 MB (.git) |
| 정리 후 | 133 | 3.1 MB |

## 9. 예상 산출물

| 파일 | 내용 |
|---|---|
| `results/spectra_all.csv` | 파장별 연속 스펙트럼 (모든 조건) |
| `results/transitions.csv` | 전이 상세 (토토머·컨포머·에너지·f·주요 오비탈 전이) |
| `results/conformer_energies.csv` | 컨포머 에너지 및 볼츠만 가중치 |
| `results/comparison.png` | 비교 그래프 |
| `results/calc_settings.json` | 모든 계산 설정 |
| `results/report.md` | 최종 분석 보고서 |

---

## 10. 알려진 한계

* **TD-DFT 는 이 발색단에 대해 계통 오차가 크다.**
  선행 연구에서 CAM-B3LYP/TZVP(진공)는 에놀 밴드를 약 +0.6 eV(약 −50 nm) 청색으로,
  B3LYP/6-31+G(d)/PCM 은 +3~+6 nm 적색으로 예측했다. 즉 **함수 선택이 결과를 지배**한다.
* **선형응답 PCM 은 근사**다. 상태특이(state-specific) 용매화나 명시적 용매 분자,
  특히 에탄올·메탄올이 킬레이트 O–H 와 만드는 분자간 수소결합은 반영되지 않는다.
* **수직 전이 근사**: 진동 구조(Franck–Condon)를 계산하지 않고 가우시안 broadening 으로
  대체하므로, 실험 스펙트럼의 밴드 모양과 λmax 는 정확히 일치하지 않을 수 있다.
* **이중 여기 / 전하이동 상태**는 TD-DFT 가 근본적으로 부정확하다.
* **비교 대상 MINDO/3–TDA 값이 문헌에 없다.** 아보벤존에 대해 발표된 반경험적
  계산 λmax 를 찾지 못했다 (`inputs/experimental_reference.json` 의 `MINDO3_TDA` 참고).
  비교하려면 사용자가 가진 값을 넣거나 저비용 기준선을 직접 계산해야 한다.
* GPU 가속이 없고 물리 4코어이므로 def2-TZVP 이상의 TD-DFT 는 현실적으로 무겁다.

---

## 11. 참고문헌

`inputs/experimental_reference.json` 의 `references` 항목에 DOI/URL 과 함께 정리되어 있다.
주요 출처:

* Vallejo et al., *Vitae* **18** (2011) 63 — 에탄올 354.9 nm, ε 34,000
* Mturi & Martincigh, *J. Photochem. Photobiol. A* **200** (2008) 410 — 용매별 λmax, 토토머 평형
* Kojić, Petković & Etinski, *J. Serb. Chem. Soc.* **81** (2016) 1393 — CAM-B3LYP/TZVP
* *ACS Omega* (2026), doi:10.1021/acsomega.5c09234 — B3LYP/6-31+G(d)/PCM
* Wong et al., *J. Phys. Chem. A* **124** (2020) — ωB97X-D/def2-SVP
