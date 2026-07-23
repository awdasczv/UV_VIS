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
#    .\scripts\run.ps1 scripts\05b_tddft_orca.py --tautomers enolA
#    .\scripts\run.ps1 -Molecule dhhb scripts\02_conformer_search.py --tautomers dhhb
#
#  -Molecule <이름> 은 환경변수 UV_MOLECULE 를 세팅한다. 파이썬 쪽 qc_common.py 가
#  이 값을 읽어 molecules/<이름>/ 아래로 모든 경로를 잡는다 (기본 avobenzone).
# =============================================================
param(
    [string]$Molecule = 'avobenzone',
    [Parameter(ValueFromRemainingArguments = $true)]
    [string[]]$Args
)

$root = Split-Path -Parent $PSScriptRoot

# -Molecule 은 위치 파라미터라, 이름 없이 부르면 첫 토큰(대개 스크립트 경로)을
# 삼켜버린다. $Molecule 이 실제 molecules/<name> 디렉터리가 아니면 그것은 분자명이
# 아니라 인자다 -> 인자 목록 앞으로 되돌리고 기본 분자를 쓴다.
# (명시적으로 '-Molecule dhhb' 를 준 경우엔 dhhb 가 실제 디렉터리이므로 그대로 유지)
# 실제 분자는 molecules/<name>/config.json 을 가진다. config.json 이 없으면
# 그 토큰은 분자명이 아니라 인자다 (대개 스크립트 경로).
if (-not (Test-Path (Join-Path $root "molecules\$Molecule\config.json"))) {
    if ($Molecule) { $Args = @($Molecule) + $Args }
    $Molecule = 'avobenzone'
}

$env:MAMBA_ROOT_PREFIX = Join-Path $root '.mamba'
$envdir = Join-Path $root '.mamba\envs\qc'
$env:UV_MOLECULE = $Molecule

$env:PATH = @(
    $envdir,
    (Join-Path $envdir 'Library\bin'),
    (Join-Path $envdir 'Library\usr\bin'),
    (Join-Path $envdir 'Scripts')
) -join ';' | ForEach-Object { "$_;$env:PATH" }

# ORCA 6.1 + Microsoft MPI
#   ORCA 는 병렬 실행 시 반드시 전체 경로로 호출해야 하고, 설치 폴더가 PATH 에
#   있어야 orca_* 하위 모듈을 찾는다. MS-MPI 는 설치 후 재로그인 전까지
#   PATH 에 안 올라오는 경우가 있어 여기서 직접 넣어준다.
foreach ($p in @('C:\ORCA_6.1.1', 'C:\Program Files\Microsoft MPI\Bin')) {
    if (Test-Path $p) { $env:PATH = "$p;$env:PATH" }
}
if (Test-Path 'C:\ORCA_6.1.1\orca.exe') { $env:ORCA_EXE = 'C:\ORCA_6.1.1\orca.exe' }

# 노트북 4코어 / 8스레드. 물리 코어 수에 맞춘다.
$env:OMP_NUM_THREADS   = '4'
$env:MKL_NUM_THREADS   = '4'
$env:PSI_SCRATCH       = Join-Path $root "molecules\$Molecule\calculations\_scratch"
New-Item -ItemType Directory -Force -Path $env:PSI_SCRATCH | Out-Null

$py = Join-Path $envdir 'python.exe'
$quoted = ($Args | ForEach-Object { '"' + $_ + '"' }) -join ' '
cmd /c "`"$py`" $quoted 2>&1"
exit $LASTEXITCODE
