#!/usr/bin/env bash
# Validation: redock imatinib into c-Abl (PDB 1IEP), the official Vina tutorial system.
# Expected: Vina score ≈ -13 kcal/mol, top-pose RMSD < 1 Å  (validated: -13.28, 0.26 Å).
set -euo pipefail
cd "$(dirname "$0")/1iep"
python ../../vdock.py validate \
  --receptor 1iep_receptor_noH.pdb \
  --ligand   1iep_ligand.sdf \
  --skip_scrub \
  --exhaustiveness 32 --seed 42 \
  --outdir ../../runs/1iep_validation
