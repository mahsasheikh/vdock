#!/usr/bin/env bash
# vdock project template (Linux/macOS/WSL). Usage: bash run_project_template.sh sites|prep|dock|analyse
set -euo pipefail
# ---------------- SETTINGS ----------------
VD="$HOME/vina-docking/vdock.py"
PDB_ID=1IEP
CHAIN=A
BOX=sites/lig_STI_A201_box.txt     # chosen after the "sites" step
LIGANDS=ligands.smi
EXH=32
SCORING=vina                        # vina | vinardo | ad4
# ------------------------------------------
L="--log vdock.log"
case "${1:-}" in
  sites)   python "$VD" $L fetch $PDB_ID --outdir input
           python "$VD" $L sites -i input/$PDB_ID.pdb --chains $CHAIN --outdir sites
           echo "Choose a site in sites/sites.csv, set BOX, then: $0 prep" ;;
  prep)    AD4=""; [ "$SCORING" = ad4 ] && AD4="--ad4"
           python "$VD" $L prep-receptor -i input/$PDB_ID.pdb --chains $CHAIN --ph 7.4 --box $BOX --outdir receptor --name receptor $AD4
           python "$VD" $L prep-ligands -i $LIGANDS --ph 7.4 --outdir ligands ;;
  dock)    if [ "$SCORING" = ad4 ]; then
             python "$VD" $L dock --ligands ligands/pdbqt --scoring ad4 --maps receptor/receptor --exhaustiveness $EXH --outdir docking
           else
             python "$VD" $L dock --ligands ligands/pdbqt --receptor receptor/receptor.pdbqt --box $BOX --scoring $SCORING --exhaustiveness $EXH --outdir docking
           fi ;;
  analyse) python "$VD" summarize --dockdir docking --collapse_states
           python "$VD" interactions --receptor receptor/receptor_clean.pdb --dockdir docking
           python "$VD" report --receptor receptor/receptor_clean.pdb --dockdir docking --box $BOX --title "$PDB_ID docking" ;;
  *) echo "Usage: $0 sites|prep|dock|analyse"; exit 1 ;;
esac
