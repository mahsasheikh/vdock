# vdock installer for native Windows 10/11 (x64).
# Run from the vina-docking folder:   powershell -ExecutionPolicy Bypass -File install-windows.ps1
param([string]$EnvName = "vina")
$ErrorActionPreference = "Stop"
$Here = Split-Path -Parent $MyInvocation.MyCommand.Path
Set-Location $Here
$VinaUrl = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.7/vina_1.2.7_win.exe"

function Find-Conda {
    $c = Get-Command conda -ErrorAction SilentlyContinue
    if ($c) { return $c.Source }
    foreach ($root in @("$env:USERPROFILE\miniforge3", "$env:LOCALAPPDATA\miniforge3",
                        "$env:ProgramData\miniforge3", "$env:USERPROFILE\miniconda3",
                        "$env:USERPROFILE\anaconda3")) {
        $exe = Join-Path $root "Scripts\conda.exe"
        if (Test-Path $exe) { return $exe }
    }
    return $null
}

Write-Host "== 1/4  Looking for conda (Miniforge) ==" -ForegroundColor Cyan
$conda = Find-Conda
if (-not $conda) {
    Write-Host "Miniforge not found. Trying to install it with winget ..."
    if (Get-Command winget -ErrorAction SilentlyContinue) {
        winget install -e --id CondaForge.Miniforge3 --accept-package-agreements --accept-source-agreements
        $conda = Find-Conda
    }
    if (-not $conda) {
        Write-Host "Please install Miniforge manually from https://conda-forge.org/download/ (Windows x86_64)," -ForegroundColor Yellow
        Write-Host "then re-run this script." -ForegroundColor Yellow
        exit 1
    }
}
Write-Host "Using $conda"

Write-Host "== 2/4  Creating/updating conda env '$EnvName' (several minutes) ==" -ForegroundColor Cyan
$envs = & $conda env list
if ($envs -match "^\s*$EnvName\s") {
    & $conda env update -n $EnvName -f env\environment-windows.yml
} else {
    & $conda env create -n $EnvName -f env\environment-windows.yml
}
if ($LASTEXITCODE -ne 0) { throw "conda environment creation failed" }

Write-Host "== 3/4  Downloading official AutoDock Vina 1.2.7 for Windows ==" -ForegroundColor Cyan
New-Item -ItemType Directory -Force -Path "$Here\bin" | Out-Null
$exe = "$Here\bin\vina.exe"
if (-not (Test-Path $exe)) {
    [Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12
    Invoke-WebRequest -Uri $VinaUrl -OutFile $exe -UseBasicParsing
}
if ((Get-Item $exe).Length -lt 500000) { throw "vina.exe download looks incomplete - delete bin\vina.exe and retry" }
Unblock-File $exe -ErrorAction SilentlyContinue

Write-Host "== 4/4  Verifying ==" -ForegroundColor Cyan
& $exe --version
& $conda run -n $EnvName python -c "import meeko, rdkit, molscrub, pdbfixer, gemmi, prolif; print('meeko', meeko.__version__, '| rdkit', rdkit.__version__)"
if ($LASTEXITCODE -ne 0) { throw "Python package check failed" }
& $conda run -n $EnvName python vdock.py -V

Write-Host ""
Write-Host "Installed successfully." -ForegroundColor Green
Write-Host "Next: open 'Miniforge Prompt', then:"
Write-Host "   conda activate $EnvName"
Write-Host "   cd `"$Here`""
Write-Host "   examples\run_1iep_redock.bat"
