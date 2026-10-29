# vdock 1.1.1-safety — 3BM9 troubleshooting and patch notes

This is a **safety-focused source-code patch** to `vdock.py` from the uploaded `vina-docking-windows-v1.1.zip`. It is not a new docking scoring function and does not guarantee top-1 crystal-pose recovery.

## Your actual 3BM9 results

- The first run used **BXZ_ideal.sdf** to construct a box near `(0.675, 0.308, -0.185)` Å. This is not in the 3BM9 crystal ligand coordinate frame, so the initial ~−1.5 scores are unusable.
- The corrected box, center `(1.377, 35.302, 25.178)` Å and size `(18, 18, 18)` Å, was obtained from the crystal ligand.
- With the corrected box, the BXZ reference produces a near-native second-ranked pose: **pose 2 RMSD 0.447 Å, score −7.702**.
- Pose 1 scores **−7.773** but has RMSD **6.635 Å**. This is a **top-1 ranking failure, not a complete docking-search failure**.
- Native-start local optimization retains a crystal-like pose, RMSD **0.084 Å**, but this is not independent blind redocking evidence.

## Changes to the original `vdock.py`

1. `box --ligand ... --structure input/3BM9.cif`: checks whether the ligand is near protein atoms in the same coordinate system. A distant ideal ligand file is rejected instead of generating a misleading box. Without `--structure`, a warning is issued.
2. `dock`: checks whether the docking box is near the receptor PDBQT before running. This is a coarse safeguard, not a guarantee that you chose the biologically correct pocket.
3. `dock`: writes `run_signature.json` containing settings, file hashes, and Vina version. Existing poses without a matching signature are refused, preventing reuse under changed receptors, boxes, ligand inputs, or parameters. **Always use a NEW `--outdir` for a changed protocol.**
4. `dock`: errors when PDBQT output has no score, SDF export fails, or one or more dockings fail. Failed ligands are listed in `docking_failures.tsv` and the command exits nonzero.
5. `summarize --reference`: includes **top-1 RMSD**, **minimum RMSD across saved poses**, and **rank of the lowest-RMSD pose**, each labeled separately.
6. `validate`: uses the **globally best-scoring state** for top-1 validation, and independently reports the minimum RMSD among saved poses. It will still return a failure code when top-1 RMSD exceeds 2 Å, even if another pose recovers the crystal geometry.

The patch **does not** repair state population weighting, metal/cofactor modeling, crystallographic-water selection, or physical ranking errors. Those require explicit chemical modeling and independent benchmarking.

## Windows CMD: install patch without losing your existing files

Run in a CMD / Miniforge Prompt. Download the patched `vdock.py` and put it in a location you can access, for example `%USERPROFILE%\Downloads\vdock.py`.

```bat
cd /d C:\docking\vina-docking
copy vdock.py vdock_v1.1.0_backup.py
copy /Y "%USERPROFILE%\Downloads\vdock.py" vdock.py
conda activate vina
python vdock.py -V
```

The version should display `vdock 1.1.1-safety`.

If using the full patched ZIP instead, unzip to a **new folder**, then copy only its `vina-docking\vdock.py` into your existing installation (after backing up the original). Do not overwrite your working environment or the installed Vina executable.

## Validate 3BM9 with a fresh output directory

```bat
cd /d C:\docking\projects\target_3bm9
set VD=C:\docking\vina-docking\vdock.py

rem Verify the original crystal-instance ligand is spatially aligned:
python %VD% box --ligand input\reference_crystal.sdf --structure input\3BM9.cif --padding 5 --out box_checked.txt

rem Do not reuse a pre-patch validation output folder:
python %VD% validate --receptor input\3BM9.cif --ligand input\reference_crystal.sdf --chains A --exhaustiveness 32 --seed 42 --outdir validation_BXZ_111

type validation_BXZ_111\docking\results.csv
```

Expect the output to distinguish `rmsd_top_pose_A` (previously ~6.6 Å) from `rmsd_best_saved_pose_A` (previously ~0.45 Å). A different run can produce different scores and poses. **Top-1 may still FAIL**; this is scientifically honest. You can still examine the correct pose among the saved poses using `redocking_diagnostics.py`.

## What you should do before screening novel ligands

1. Keep the **correct crystal-coordinate box**; never derive a pocket box from CCD/ideal ligand coordinates.
2. Verify the reference SDF is the **3BM9 BXZ instance SDF** and matches the protein coordinate frame.
3. Inspect active-site water molecules, especially possible HSP90 water-mediated contacts. The standard receptor preparation removes them; `--keep_water` is **not** automatically a correct scientific fix.
4. Repeat BXZ redocking with independent Vina seeds, keeping the setup constant. Record top-1 and top-9 success separately.
5. Evaluate an independent set of known HSP90 actives and inactives/decoys and report enrichment before trusting any score ranking for unseen molecules.
6. Use a **fresh `--outdir`** for every change of receptor, box, ligand state, or scoring model.

A reference-guided pose choice can be appropriate for diagnosing BXZ itself, but **you cannot choose the lowest-RMSD pose for an unknown screening molecule** because its crystal pose is not known.

## Testing performed here

- Python syntax compilation of patched `vdock.py` passed.
- Crystal-aligned bundled 1IEP ligand passed the spatial box check (~1.8 Å closest protein atom).
- Artificially displaced 1IEP ligand (>100 Å) was rejected.
- Per-pose RMSD helper produced `0.0` Å for identical conformations and ~17.32 Å for a test pose translated by (10, 10, 10) Å.
- No complete new AutoDock Vina, molscrub, or Meeko run was executed in this environment. Windows installation and numerical results remain to be verified locally.
