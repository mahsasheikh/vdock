# vdock: AutoDock Vina docking setup

A reproducible, validated docking pipeline built on the current upstream toolchain from the
CCSB/Forli lab (Scripps):

| Stage | Tool | Version |
|---|---|---|
| Docking engine | AutoDock Vina (CLI + Python) | 1.2.7 |
| Receptor & ligand PDBQT, SDF export | Meeko | 0.8.0 |
| Ligand protonation / tautomers / 3D | molscrub | 0.3.2 |
| Receptor clean-up & protonation | PDBFixer + OpenMM | 1.12 / 8.6 |
| AD4 grid maps (optional) | AutoGrid | 4.2.9 |
| Interaction fingerprints | ProLIF | 2.2.3 |
| 3D viewer in HTML report | 3Dmol.js (bundled, BSD) | 2.5.3 |
| Cheminformatics | RDKit, Python 3.11 | 2026.03 |

Everything comes from **conda-forge**, with molscrub (and Meeko on Apple Silicon) from PyPI. Legacy
MGLTools / `prepare_receptor4.py` and Open Babel PDBQT conversion are **not** used; upstream no
longer recommends them.

## Install

```bash
bash install.sh              # creates env "vina"; installs micromamba if no conda is present
micromamba activate vina     # or: conda activate vina
bash examples/run_1iep_redock.sh   # validation, ~3–10 min
```

Platforms: Linux x86-64/ARM, macOS Intel, and macOS Apple Silicon (separate env file, see
`env/`). **Windows:** run `install-windows.bat` (native; uses the official `vina_1.2.7_win.exe`) and
follow **SETUP-WINDOWS.md**. You can also use WSL2 with `install.sh`.

## Pipeline

```
fetch ─► sites/box ─► prep-receptor ─► prep-ligands ─► dock ─► summarize ─► interactions ─► report
                         (validate = redock a known complex; do it once per new target)
```

**New users: read [GUIDE.md](GUIDE.md).** It's a step-by-step walkthrough covering choosing a structure,
finding the site, preparation, docking, interpreting results, interactions, visualisation, and troubleshooting.

```bash
python vdock.py fetch 1IEP --outdir input
python vdock.py box --ligand ref_ligand.sdf --padding 5 --out box.txt            # or --structure + --residues A:315,318
python vdock.py prep-receptor -i input/1IEP.pdb --chains A --ph 7.4 --box box.txt --outdir receptor --name receptor
python vdock.py prep-ligands -i library.smi --ph 7.4 --outdir ligands
python vdock.py dock --ligands ligands/pdbqt --receptor receptor/receptor.pdbqt --box box.txt --exhaustiveness 32
python vdock.py summarize --dockdir docking --collapse_states                      # -> docking/results.csv
python vdock.py interactions --receptor receptor/receptor_clean.pdb --dockdir docking    # ProLIF contacts
python vdock.py report --receptor receptor/receptor_clean.pdb --dockdir docking --box box.txt   # HTML + PyMOL
```

`sites` lists co-crystallised ligands and detects pockets geometrically (grid buriedness, so it works
on apo structures and AlphaFold models). It writes a ready-to-use box file for each site.

Each step appends its parameters and exact tool versions to a `manifest.json` in its output folder,
for provenance. Docking can be **resumed**: ligands that already have output are skipped unless you
pass `--overwrite`. Per-ligand Vina logs are written to `docking/logs/`, and poses go to
`docking/poses/*.pdbqt` and `*.sdf`. The SDF files have correct bond orders.

`examples/run_project_template.{bat,sh}` is an editable 4-stage template for a new target (sites, prep, dock, analyse).

## Defaults and why (following current guidelines)

- **Receptor:** waters and heteroatoms are removed, non-standard residues are mapped to standard ones,
  missing heavy atoms are rebuilt, and hydrogens are added at pH 7.4 (PDBFixer). Meeko then matches
  each residue against chemical templates, so problems fail loudly rather than silently. Missing
  loops are *not* modelled unless you pass `--fill_gaps`. Check His/Cys states in your binding site,
  and override them with Meeko templates if needed (e.g. `--delete_residues`, or edit the clean PDB).
- **Ligands:** start from SMILES or SDF, never from PDB, because PDB has no bond orders.
  molscrub enumerates the relevant protonation states and tautomers at pH 7.4 and builds 3D coordinates
  (ETKDGv3). States that share a parent get `_s1`, `_s2` suffixes, and `--collapse_states` keeps
  the best-scoring state per parent.
- **Box:** the ligand or residue extent plus 5 Å per side, with a minimum edge of 18 Å. You'll get a
  warning above 30×30×30 Å³. `box.pdb` lets you view the box in PyMOL or ChimeraX.
- **Search:** `exhaustiveness 32` (the upstream tutorial value for reproducible poses), 9 modes,
  3 kcal/mol energy range, **fixed seed 42**. For large screens, 8–16 is acceptable for a first
  pass; re-dock your hits at 32.
- **Parallelism:** the default is `cores / min(cores, exhaustiveness)` concurrent runs. For big
  libraries, many single-thread jobs are faster (`--cpu_per_job 1`).
- **Scoring functions:** `vina` (default), `vinardo`, or `ad4` (`prep-receptor --ad4` builds the
  maps). Scores from different scoring functions **are not comparable**.
- **Validate first:** redock a co-crystallised ligand for every new target. A top-pose RMSD of
  2 Å or less means the protocol is acceptable. RMSD is computed in place (no re-alignment), is
  symmetry-aware, and is tolerant to protonation differences.

## Validation results (on this setup)

PDB 1IEP, c-Abl kinase domain with imatinib (the official Vina tutorial system):

| Run | Score (kcal/mol) | Upstream reference | Top-pose RMSD |
|---|---|---|---|
| Vina, crystal ligand, exh 32 | −13.28 | ≈ −13.2 | **0.26 Å** |
| AD4 maps, ligand from SMILES via molscrub, exh 8 | −14.66 | ≈ −14.7 | **0.47 Å** |
| Vina, ligand from SMILES via molscrub, exh 8 | −12.92 | – | **0.47 Å** |

In the 5-compound demo screen (`examples/1iep/expected_screen_results.csv`), imatinib ranks first,
ahead of ponatinib, dasatinib, ibuprofen and aspirin. That order is expected for this DFG-out structure.

## Caveats

Docking scores are for ranking; they are not binding free energies. Inspect top poses visually, and
use consensus or rescoring (or MD/FEP) before drawing conclusions. Metal ions, covalent ligands,
macrocycles, cofactors and explicit waters need extra setup: see the Meeko docs (`--flexres`,
reactive residues, `-w` hydrated docking, custom templates).

## Citations

- Eberhardt J, Santos-Martins D, Tillack AF, Forli S. *AutoDock Vina 1.2.0: new docking methods, expanded
  force field, and Python bindings.* J Chem Inf Model 2021, 61, 3891–3898.
- Trott O, Olson AJ. *AutoDock Vina: improving the speed and accuracy of docking…* J Comput Chem 2010, 31, 455–461.
- Forli S et al. *Computational protein–ligand docking and virtual drug screening with the AutoDock suite.*
  Nat Protoc 2016, 11, 905–919.
- Meeko and molscrub: https://github.com/forlilab/Meeko, https://github.com/forlilab/molscrub
- Bouysset C, Fiorucci S. *ProLIF: a library to encode molecular interactions as fingerprints.* J Cheminform 2021, 13, 72.
- Rego N, Koes D. *3Dmol.js: molecular visualization with WebGL.* Bioinformatics 2015, 31, 1322–1324.
- Eastman P et al. *OpenMM 8.* J Phys Chem B 2024 (PDBFixer).
