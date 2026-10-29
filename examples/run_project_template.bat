@echo off
REM ==========================================================================
REM  vdock project template (Windows).  See GUIDE.md for what each step does.
REM  1) copy this file into a new project folder, e.g. C:\docking\projects\abl1\
REM  2) edit the SETTINGS block, 3) in Miniforge Prompt:  conda activate vina
REM  4) run it step by step:   run_project_template.bat sites   (then check!)
REM                            run_project_template.bat prep
REM                            run_project_template.bat dock
REM                            run_project_template.bat analyse
REM ==========================================================================
setlocal
REM ---------------- SETTINGS -------------------------------------------------
set VD=C:\docking\vina-docking\vdock.py
set PDB_ID=1IEP
set CHAIN=A
REM  box file chosen after the "sites" step (e.g. sites\lig_STI_A201_box.txt or sites\pocket01_box.txt)
set BOX=sites\lig_STI_A201_box.txt
REM  ligands: .smi (SMILES name per line), .sdf, .csv, .xlsx
set LIGANDS=ligands.smi
set EXH=32
set SCORING=vina
REM ---------------------------------------------------------------------------
set LOG=--log vdock.log

if "%1"=="sites" (
  python "%VD%" %LOG% fetch %PDB_ID% --outdir input || exit /b 1
  python "%VD%" %LOG% sites -i input\%PDB_ID%.pdb --chains %CHAIN% --outdir sites || exit /b 1
  echo.
  echo Now look at sites\sites.csv, choose a site, set BOX above, then run:  %~nx0 prep
  exit /b 0
)
if "%1"=="prep" (
  python "%VD%" %LOG% prep-receptor -i input\%PDB_ID%.pdb --chains %CHAIN% --ph 7.4 --box %BOX% --outdir receptor --name receptor || exit /b 1
  python "%VD%" %LOG% prep-ligands -i %LIGANDS% --ph 7.4 --outdir ligands || exit /b 1
  echo Check ligands\prep_failures.tsv (if present). Next:  %~nx0 dock
  exit /b 0
)
if "%1"=="dock" (
  python "%VD%" %LOG% dock --ligands ligands\pdbqt --receptor receptor\receptor.pdbqt --box %BOX% --scoring %SCORING% --exhaustiveness %EXH% --outdir docking || exit /b 1
  echo Next:  %~nx0 analyse
  exit /b 0
)
if "%1"=="analyse" (
  python "%VD%" summarize --dockdir docking --collapse_states || exit /b 1
  python "%VD%" interactions --receptor receptor\receptor_clean.pdb --dockdir docking || exit /b 1
  python "%VD%" report --receptor receptor\receptor_clean.pdb --dockdir docking --box %BOX% --title "%PDB_ID% docking" || exit /b 1
  start "" docking\report.html
  exit /b 0
)
echo Usage: %~nx0 sites ^| prep ^| dock ^| analyse
