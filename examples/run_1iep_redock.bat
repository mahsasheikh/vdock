@echo off
REM Validation: redock imatinib into c-Abl (PDB 1IEP). Run inside "conda activate vina".
REM Expected: score about -13 kcal/mol, "VALIDATION PASS" with RMSD below 1 A.
cd /d "%~dp0\1iep"
python ..\..\vdock.py validate --receptor 1iep_receptor_noH.pdb --ligand 1iep_ligand.sdf --skip_scrub --exhaustiveness 32 --seed 42 --outdir ..\..\runs\1iep_validation
