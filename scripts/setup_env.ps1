# =============================================================
#  setup_env.ps1 - 로컬(사용자 영역) 양자화학 계산 환경 구축
#  시스템 전역 설정을 바꾸지 않는다. 모든 것이 프로젝트 폴더 안에 설치된다.
#   - micromamba (단일 exe, 사용자 영역)  -> tools/
#   - conda-forge 환경 'qc'               -> .mamba/envs/qc
#     psi4 (DFT/TD-DFT), xtb (GFN2-xTB), rdkit (구조 생성), numpy/scipy/pandas/matplotlib
# =============================================================
$ErrorActionPreference = 'Stop'
$ProgressPreference    = 'SilentlyContinue'

$root   = Split-Path -Parent $PSScriptRoot
$tools  = Join-Path $root 'tools'
$mmroot = Join-Path $root '.mamba'
$mm     = Join-Path $tools 'micromamba.exe'

New-Item -ItemType Directory -Force -Path $tools, $mmroot | Out-Null

if (-not (Test-Path $mm)) {
    Write-Host "[1/3] micromamba 내려받는 중..."
    Invoke-WebRequest `
        -Uri 'https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-win-64.exe' `
        -OutFile $mm
}
& $mm --version

$env:MAMBA_ROOT_PREFIX = $mmroot

Write-Host "[2/3] conda-forge 환경 'qc' 생성 중 (수 분 소요)..."
& $mm create -y -n qc -c conda-forge --override-channels `
    "python=3.11" psi4 xtb rdkit numpy scipy pandas matplotlib

Write-Host "[3/3] 설치 확인"
$py = Join-Path $mmroot 'envs\qc\python.exe'
& $py -c "import psi4, rdkit, numpy, scipy, pandas, matplotlib; print('psi4', psi4.__version__); print('rdkit', rdkit.__version__)"
# xtb 는 배너를 stderr 로 찍는다. PowerShell 5.1 이 이를 오류로 오해하므로 감싸준다.
$xtb = Join-Path $mmroot 'envs\qc\Library\bin\xtb.exe'
$ver = (cmd /c "`"$xtb`" --version 2>&1") -join "`n"
($ver -split "`n" | Where-Object { $_ -match 'xtb version' })

Write-Host ""
Write-Host "환경 준비 완료:"
Write-Host "  python : $py"
Write-Host "  xtb    : $xtb"
