# 아보벤존(Avobenzone) 이론 UV–Vis 흡광 스펙트럼 계산

로컬 PC에서 아보벤존(butyl methoxydibenzoylmethane, CAS 70356-09-1)의
**킬레이트 에놀형 두 가지**와 **디케토형**의 UV–Vis 흡수 스펙트럼을
DFT / TD-DFT 수준으로 계산하고, 알려진 실험 최대흡수파장과 비교한다.

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

### 3.3 최종 선택

```
구조 생성      RDKit ETKDGv3 + MMFF94s  (킬레이트 이면각 명시적 세팅)
컨포머 탐색    GFN2-xTB (ALPB ethanol) — 45원자 1구조 최적화에 약 1.8 초
구조 최적화    DFT + PCM(ethanol)        ← 아래 4장
들뜬상태       TD-DFT (TDA/RPA) + PCM(ethanol, Nonequilibrium)
```

CREST 는 Windows 에서 못 쓰므로 **RDKit ETKDG(다수 구조) → GFN2-xTB 최적화 →
에너지 + heavy-atom RMSD 중복제거** 로 대체한다. xTB 가 구조당 2초 수준이라
수백 개 구조를 돌려도 몇 분이면 끝나므로 실용적인 대안이 된다.
(진짜 CREST 메타다이내믹스가 필요하면 WSL 옵션 — 8장 참고)

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

### 노트북(현재 PC)에서 실제로 수행하는 것

| 단계 | 수준 | 근거 |
|---|---|---|
| 컨포머 탐색 | RDKit ETKDGv3 (200 구조) → **GFN2-xTB / ALPB(에탄올)** 최적화 → 에너지 + heavy-atom RMSD 중복제거 → 298.15 K 볼츠만 가중치 | 구조당 1.8초로 사실상 공짜. CREST 부재를 실용적으로 대체 |
| 구조 | GFN2-xTB 최적화 구조 (+ 최저 컨포머는 DFT 재최적화) | DFT 최적화가 구조당 수 시간이라 전면 적용 불가 |
| TD-DFT (기체상) | **B3LYP/6-31G(d)** 및 **CAM-B3LYP/6-31G(d)**, TDA, 12 상태 | 두 함수를 비교해야 "함수 선택의 영향"을 말할 수 있다. 선행 연구에서 B3LYP 는 적색, CAM-B3LYP 는 청색으로 치우침이 알려져 있음 |
| TD-DFT (용매) | 위와 동일 + **IEF-PCM(에탄올), Nonequilibrium, Area 1.0 Å²** | 수직 전이이므로 비평형 용매화. Area 는 비용/정확도 절충으로 1.0 채택 |
| 저비용 기준선 | **HF/STO-3G + CIS(TDA)** | MINDO/3–TDA 문헌값이 없어 동급 저비용 조합을 자체 계산 |

**왜 6-31G(d) 인가**: 실측한 기저함수 개수(6-31G(d) = 389, def2-TZVP = 845)와
TD-DFT 비용이 기저함수 수의 3~4제곱에 가깝게 늘어나는 점을 고려하면,
def2-TZVP 는 이 노트북에서 PCM 과 함께 쓸 수 없다. 6-31G(d) 는 편극함수를 포함해
π→π* 전이를 정성적으로 옳게 기술하는 최소한의 수준이다.
선행 연구(ACS Omega 2026)가 6-31+G(d) 로 실험과 +3~+6 nm 오차를 얻었으므로
비교 기준으로도 적절하다.

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

모든 스크립트는 `scripts/run.ps1` 래퍼로 실행한다.
(환경의 `Library\bin` 을 PATH에 넣어 `s-dftd3` / `dftd4` / `xtb` 를 찾게 하고,
스레드 수와 scratch 경로를 고정한다.)

```powershell
# 0) 엔진 능력 점검 (함수/기저셋/TD-DFT/PCM 실제 동작 확인)
.\scripts\run.ps1 scripts\00_probe_engine.py

# 1) 세 토토머의 3D 구조 생성 + 연결·수소·킬레이트 검증
.\scripts\run.ps1 scripts\01_build_structures.py

# 2) 컨포머 탐색 + 볼츠만 가중치 + 대표 컨포머 선택
.\scripts\run.ps1 scripts\02_conformer_search.py --backend etkdg --nconf 200

# 3) 비용 측정용 최저비용 단일점 테스트
.\scripts\run.ps1 scripts\03_test_single_point.py
```

*(이후 단계는 추가되는 대로 여기에 기록한다.)*

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
