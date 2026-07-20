# 다음에 이어서 할 일 (인수인계)

작성: 2026-07-20. 계산 중단 시점 기준.

---

## 현재 진행 상태

| 단계 | 상태 |
|---|---|
| 환경 구축 (Psi4 1.11 + xTB 6.7.1 + RDKit) | ✅ 완료 |
| 세 토토머 구조 생성·검증 | ✅ 완료 |
| 엔진 능력 확인 (함수/기저셋/TD-DFT/PCM) | ✅ 완료 |
| 비용 측정 및 계산 수준 확정 | ✅ 완료 |
| 컨포머 탐색 (200 구조 × 3 토토머) | ✅ 완료 |
| 컨포머 검증 + 발색단 클러스터링 대표 선택 | ✅ 완료 (대표 6개) |
| **TD-DFT 본 계산** | ⏸ **1/12 완료** |
| 저비용 기준선 (HF/STO-3G + CIS) | ⬜ 미실행 |
| 스펙트럼·CSV·그래프 | ⬜ 미실행 |
| 최종 보고서 | ⬜ 미실행 |

완료된 TD-DFT: `calculations/02_tddft/enolA/enolA_c000/b3lyp_631gd/none/result.json`

---

## ⚠️ 미해결 문제 1: PCMSolver 의 S 행렬 특이성 (최우선)

### 증상

PCM 을 켠 TD-DFT 실행 중 stdout 에 다음이 2회 나온다.

```
S matrix is not positive-definite!
Consider changing the average area of the cavity finite elements.
Normally an error, this has been commuted to a warning via patch.
Please report this issue: https://github.com/PCMSolver/pcmsolver/issues/206
```

원래는 **에러**인데 패치로 경고로 낮춘 것이다. 계산은 진행되지만
**결과를 그대로 믿으면 안 된다.**

### 지금까지 파악한 것

경고는 `Area` 값만의 문제가 아니라 **구조에 따라 달라진다**.

| 로그 | 구조 | Area | 표면조각 | 경고 |
|---|---|---|---|---|
| `logs/03_test_single_point.log` | `enolA_xtbopt` | 0.3 | 2149 | **2회** |
| `logs/03b_pcm_cost_probe.log` | `enolA_xtbopt` | 1.0 | 1024 | 0회 |
| `logs/03b_pcm_cost_probe.log` | `enolA_xtbopt` | 2.0 | 750 | 0회 |
| `logs/05_tddft_stage1.log` | `enolA_c000` | 1.0 | 1047 | **2회** |

즉 같은 Area 1.0 이라도 구조가 `enolA_xtbopt` 면 깨끗하고 `enolA_c000` 이면
경고가 난다. GePol 공동 생성에서 표면 조각이 겹치거나 퇴화해 겹침행렬이
특이해지는 것으로 보인다.

**중요**: 앞서 "Area 1.0 과 2.0 의 λmax 차이가 0.07 nm 이므로 수렴했다"고
판단한 근거는 **경고가 없던 실행들**이었으므로 그 자체는 유효하다.
문제는 새 컨포머 구조에서 경고가 난다는 점이다.

### 해야 할 일 (순서대로)

1. **경고를 자동으로 감지**해서 실패로 분류한다.
   `scripts/psi4_helpers.py` 의 `FAILURE_PATTERNS` 에
   `S matrix is not positive-definite` 패턴과 `PCM_ILL_CONDITIONED` 코드를 추가.
   현재는 이 패턴이 없어서 조용히 지나간다.
   또 PCMSolver 는 이 메시지를 psi4 출력 파일이 아니라 **stdout 으로** 흘리므로,
   `05_tddft.py` 가 계산 중 stdout 을 캡처해 검사하도록 고쳐야 한다.

2. **Area 를 바꿔가며 경고가 사라지는 값을 자동 탐색**한다.
   예: `[1.0, 0.6, 1.5, 0.3, 2.0]` 순으로 시도하고 경고 없는 첫 값을 채택.
   채택한 Area 를 `result.json` 에 반드시 기록한다.

3. 그래도 안 되면 다음 대안을 시험한다.
   - `RadiiSet = UFF` (Bondi 대신)
   - `Scaling = False`
   - `Type = GePol` 의 `MinRadius` 조정 (추가 구 생성 억제)

4. **검증**: 경고 없는 조건에서 얻은 λmax 와, 경고가 났던 계산의 λmax 를 비교해
   실제로 결과가 달랐는지 확인하고 그 차이를 보고서에 적는다.
   (경고가 무해했을 수도 있다. 확인 없이 단정하지 말 것.)

---

## ⚠️ 미해결 문제 2: 들뜬상태 개수가 부족하다

완료된 enolA 기체상 계산에서 12개 상태의 최단파장이 **218.2 nm** 였다.
요구사항은 **200~450 nm 를 덮는 것**이므로 미달이다.

→ `inputs/calc_config.json` 의 `tddft.n_states` 를 12에서 **16~18** 로 올린다.
   다만 상태 수가 늘면 Davidson 비용이 늘어나므로, 먼저 기체상에서 몇 개면
   200 nm 에 닿는지 재고 결정할 것.

---

## 실측된 비용 (계획 수립에 사용할 것)

노트북 i5-1135G7 4코어 기준, 아보벤존 45원자.

| 계산 | 조건 | 시간 |
|---|---|---|
| GFN2-xTB 최적화 | ALPB(에탄올) | 1.8 초 |
| TD-DFT 기체상 | B3LYP/6-31G(d), 389 기저함수, 12상태 | SCF 158초 + TD 654초 = **약 14분** |
| TD-DFT PCM | 같은 조건 + IEF-PCM(에탄올) | 미완료. 기체상의 약 3.8배로 추정 → **약 50분** |

남은 계산량: 대표 구조 6개 × (기체상 + 에탄올) = 12건.
대략 **6~8 시간**. 여기에 CAM-B3LYP 비교와 저비용 기준선이 추가된다.

---

## 재개 방법

```powershell
cd c:\Users\ysm04\Documents\etc\uv_vis

# 1단계 TD-DFT 이어서 실행 (완료된 것은 자동으로 건너뜀)
.\scripts\run.ps1 scripts\05_tddft.py --levels b3lyp_631gd --max-conformers 1

# 저비용 기준선 (몇 분이면 끝남, 먼저 해도 좋다)
.\scripts\run.ps1 scripts\08_lowcost_baseline.py

# 전부 끝난 뒤
.\scripts\run.ps1 scripts\06_build_spectra.py
.\scripts\run.ps1 scripts\07_report.py
```

모든 계산은 `(구조 × 이론수준 × 용매)` 조합 단위로 `result.json` 체크포인트를
남기고, 재실행 시 완료된 조합을 자동으로 건너뛴다. 중간에 꺼도 안전하다.

---

## 결정 대기 중인 사항

### ORCA 6.1 설치 여부

가장 효과적인 다음 수단이다. 근거는 `docs/DEVELOPMENT_LOG.md` 5장.

* 얻는 것: RIJCOSX 로 TD-DFT 가속, def2-SVP/TZVP 실용화, 성숙한 CPCM/SMD,
  들뜬상태 해석적 gradient.
* 필요한 사용자 조치 두 가지 (자동화 불가):
  1. `orcaforum.kofo.mpg.de` 가입 후 ORCA 6.1 Windows 버전 수동 다운로드
  2. **Microsoft MPI** 설치 — 없으면 1코어로만 돈다
* 설치하면 **똑같은 조건(enolA, B3LYP/6-31G(d), PCM, 12상태)으로 벤치마크**해서
  Psi4 대비 실제 배수를 측정할 것. 현재 Psi4 기준선: 기체상 SCF 158초 + TD 654초.

### 데스크탑(Ryzen 5 5600X + RTX 5060 Ti) 이전

* CPU 는 약 3배 이득. **GPU 는 기대만큼 도움 안 됨**
  (소비자용 Blackwell 의 FP64 가 FP32 의 1/64 로 제한되어 5600X CPU 와 동급).
* 이전하려면 `inputs/calc_config.json` 의 `active_profile` 을
  `desktop_5600x` 로 바꾸고 `scripts/setup_env.ps1` 만 실행하면 된다.
