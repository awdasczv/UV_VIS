# =============================================================
#  run.ps1 - 프로젝트 환경에서 파이썬 스크립트를 실행하는 얇은 래퍼
#
#  하는 일
#    - .mamba/envs/qc 의 bin 경로들을 PATH 앞에 붙인다.
#      (이렇게 해야 psi4 가 s-dftd3 / dftd4 실행파일을 찾아 -D3BJ 보정을 쓸 수 있고,
#       xtb.exe 도 그냥 'xtb' 로 호출된다.)
#    - 스레드 수를 고정한다 (노트북 4코어).
#    - 네이티브 실행파일이 stderr 로 배너를 찍어도 PowerShell 이 오류로 오해하지 않게 한다.
#
#  사용법:
#    .\scripts\run.ps1 scripts\03_dft_optimize.py --tautomer enolA
# =============================================================
param(
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$Args
)

$root = Split-Path -Parent $PSScriptRoot
$env:MAMBA_ROOT_PREFIX = Join-Path $root '.mamba'
$envdir = Join-Path $root '.mamba\envs\qc'

$env:PATH = @(
    $envdir,
    (Join-Path $envdir 'Library\bin'),
    (Join-Path $envdir 'Library\usr\bin'),
    (Join-Path $envdir 'Scripts')
) -join ';' | ForEach-Object { "$_;$env:PATH" }

# 노트북 4코어 / 8스레드. 물리 코어 수에 맞춘다.
$env:OMP_NUM_THREADS   = '4'
$env:MKL_NUM_THREADS   = '4'
$env:PSI_SCRATCH       = Join-Path $root 'calculations\_scratch'
New-Item -ItemType Directory -Force -Path $env:PSI_SCRATCH | Out-Null

$py = Join-Path $envdir 'python.exe'
$quoted = ($Args | ForEach-Object { '"' + $_ + '"' }) -join ' '
cmd /c "`"$py`" $quoted 2>&1"
exit $LASTEXITCODE
