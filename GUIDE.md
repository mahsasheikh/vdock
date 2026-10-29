# vdock user guide: from a new protein to interpreted results

Commands are written for **Windows (Miniforge Prompt)**. On Linux/macOS/WSL, swap `\` for `/`.
Every command has built-in help: `python vdock.py <command> -h`.

---

## 0. Every session

```bat
conda activate vina
cd C:\docking\projects\my_target
```

Throughout this guide, **`VD`** stands for the full path to `vdock.py`. Either type it in full
(e.g. `python C:\docking\vina-docking\vdock.py ...`) or set it once per session:

```bat
set VD=C:\docking\vina-docking\vdock.py
python %VD% -h
```

### Where to keep projects

Keep **projects outside OneDrive**, e.g. `C:\docking\projects\...`. A screen writes thousands of small
files, and OneDrive syncing them slows everything down and can lock files mid-run. Copy only the final
results into OneDrive.

Suggested layout (the project template creates it for you):

```
C:\docking\projects\abl1\
  ligands.smi            your input ligands
  input\1IEP.pdb         downloaded structure
  sites\                 binding-site analysis + box files
  receptor\              receptor_clean.pdb (with H), receptor.pdbqt
  ligands\               scrubbed.sdf, pdbqt\, prep_failures.tsv
  docking\               poses\, logs\, results.csv, interactions_*.csv, report.html
  vdock.log              log of everything you ran
```

The quickest start is to copy `examples\run_project_template.bat` into your project folder, edit the
SETTINGS block, and run it in four stages: `sites`, `prep`, `dock`, `analyse`. The sections
below explain each stage, so you know what you're running and what to check.

---

## 1. Choose the protein structure

Go to rcsb.org, search your protein (or its UniProt ID), and pick a structure using these criteria, in this order:

1. **Correct protein and species**, wild type (or the variant you care about).
2. **Holo structure**: crystallised *with a ligand in your site of interest*, ideally a ligand chemically similar to yours.
   The side chains are then already shaped for binding (apo structures often have collapsed pockets).
3. **Resolution ≤ 2.5 Å** (X-ray), and no missing residues *in the binding site*. On the RCSB page, check
   "Sequence" → unmodelled regions.
4. For kinases, GPCRs and similar targets, the **conformational state** matters (DFG-in/out, active/inactive). Dock into the state your ligands bind.

No experimental structure? Use an AlphaFold model (alphafold.ebi.ac.uk) as `input\model.pdb`. Only trust
regions with pLDDT > 70, and be aware that side chains in AlphaFold pockets are not ligand-optimised.

```bat
python %VD% fetch 1IEP --outdir input
```

You can also use your own PDB/mmCIF file. Just put it in `input\`.

---

## 2. Find and check the binding site

```bat
python %VD% sites -i input\1IEP.pdb --chains A --outdir sites
```

This prints a table and writes `sites\sites.csv`:

```
site                      source            #res  note
lig_STI_A201              co-crystallised     21  ligand, 37 heavy atoms
pocket01                  geometric           69  volume 913 Å^3, buriedness 0.82; overlaps STI A201
pocket02                  geometric           23  volume 190 Å^3, buriedness 0.78
```

- **`lig_...` sites** are ligands bound in the crystal structure. If your ligands target the same site,
  use this box (`sites\lig_STI_A201_box.txt`). It's the most reliable option.
- **`pocket..` sites** come from geometric detection: enclosed cavities in the protein surface, ranked by
  volume × buriedness. Use these for apo structures, AlphaFold models, or to look for other sites.
  A note like "overlaps STI A201" confirms that a pocket is the known site.
- Crystallisation additives (sulfate, glycerol, PEG, etc.) and ions are skipped. Add `--all_hetero` to list
  them, which matters if a metal ion is part of your active site.

**Check the site against the literature and UniProt.** On the UniProt entry, look under "Function → Sites"
and "Binding site" for catalytic and binding residues. The residues in `lining_residues` (in
`sites.csv`) should match.

**Look at it.** Open PyMOL (see §7) and run:

```
load input/1IEP.pdb
load sites/lig_STI_A201_box.pdb
load sites/pockets.pdb
show sticks, sites/*_box*
```

The box (8 corner atoms) should completely enclose the site, with a few Å of space around it.

**Box size rules of thumb:**
- The default is the ligand or pocket extent plus 5 Å on each side, with a minimum of 18 Å per edge.
- Ligands much larger than the reference ligand need a larger box: re-run `box` with `--padding 7`.
- Boxes larger than about 30 × 30 × 30 Å slow the search and make it less reliable. Raise `--exhaustiveness` (e.g. 64) or shrink the box.
- To set a box by hand or from residues:
  ```bat
  python %VD% box --structure input\1IEP.pdb --residues A:271,286,315,318,381 --out box.txt
  python %VD% box --center 15.2 53.9 16.9 --size 22 22 22 --out box.txt
  ```

---

## 3. Prepare the receptor

```bat
python %VD% prep-receptor -i input\1IEP.pdb --chains A --ph 7.4 --box sites\lig_STI_A201_box.txt --outdir receptor --name receptor
```

What it does:
1. **PDBFixer**: keeps the chains you list, removes waters, ligands and ions, converts modified residues
   (e.g. MSE → MET), rebuilds missing side-chain atoms, and adds hydrogens at the pH you give.
2. **Meeko**: matches every residue against chemical templates, assigns charges and atom types, and writes `receptor.pdbqt`.

Outputs: `receptor\receptor_clean.pdb` (the structure you analyse and visualise; it has hydrogens) and
`receptor\receptor.pdbqt` (the file you dock into).

Useful options:

| Need | Option |
|---|---|
| Keep crystallographic waters (rarely helps with Vina) | `--keep_water` |
| Model missing internal loops | `--fill_gaps` (check the result!) |
| Different pH (e.g. an acidic compartment) | `--ph 5.5` |
| Flexible side chains in the site | `--flexres A:315,A:381` (then dock with `--flex receptor\receptor_flex.pdbqt`) |
| A residue Meeko can't match | `--allow_bad_res` (drops it, so check it's not in the site) |
| AutoDock4 scoring (Linux/macOS/WSL only) | `--ad4` |

**Check:** open `receptor_meeko.log` if something fails. Histidine protonation states (HID/HIE/HIP) are
chosen by PDBFixer from the hydrogen-bonding network. If a His, Asp/Glu or Cys in your site has a known
special state (a catalytic His, a metal-binding Cys), check it.

---

## 4. Prepare the ligands

Input formats:

| Format | Example |
|---|---|
| `.smi` | one per line: `SMILES name` (separated by a space or tab, no header) |
| `.sdf` | from PubChem, ChemDraw, a vendor library, etc.; the name comes from the title line |
| `.csv` / `.xlsx` | SMILES in the first column (or `--column` in molscrub) |
| a single SMILES | `-i "CC(=O)Oc1ccccc1C(=O)O"` |

Example `ligands.smi`:

```
Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1 imatinib
CC(=O)Oc1ccccc1C(=O)O aspirin
```

For SMILES from PubChem, open the compound page and copy "Isomeric SMILES" under "Computed Descriptors".
**Define stereocentres.** Undefined ones are left to RDKit unless you add `--stereoisomers`, which enumerates them.

```bat
python %VD% prep-ligands -i ligands.smi --ph 7.4 --outdir ligands
```

What it does:
1. **molscrub** enumerates the relevant protonation states and tautomers at pH 7.4, generates 3D
   coordinates (ETKDGv3 + force-field minimisation), and fixes ring conformations.
   One input molecule can become several states: `name`, `name_s2`, `name_s3`.
2. **Meeko** assigns atom types, charges and rotatable bonds, then writes `ligands\pdbqt\*.pdbqt`.

**Check:**
- `ligands\prep_failures.tsv` and `ligands\scrub_failed.sdf` list molecules that failed, with the reason.
- The number of PDBQT files should be at least the number of input molecules.
- Already have protonated 3D structures (e.g. from your own QM work)? Add `--skip_scrub`.
- Fewer states (faster): `--skip_tautomers` and/or `--skip_acidbase`.

Practical limits: Vina handles up to about 32 rotatable bonds, but accuracy drops noticeably above about 10–12.
Macrocycles are treated as flexible by default; add `--rigid_macrocycles` to keep the input ring conformation.
Peptides above about 6 residues are beyond what rigid-receptor docking does well.

---

## 5. Validate the protocol for each new target (redocking)

Before trusting results for a new target, check that the setup can reproduce a known crystal pose:

1. Get the co-crystallised ligand as an SDF with correct bond orders **in the crystal coordinates**. On the
   rcsb.org structure page, find the ligand in the "Small Molecules" section and download its
   **instance** coordinates as SDF. Don't use the "ideal" coordinates: they aren't positioned in the
   binding site, so the RMSD would be meaningless. Don't use a PDB-format ligand either.
2. Run (the structure can still contain the ligand, since all hetero groups are removed from the receptor):

```bat
python %VD% validate --receptor input\1IEP.pdb --chains A --ligand STI_A_201.sdf --outdir validation
```

Ligand files without hydrogens (the usual case for RCSB files) are fine: molscrub adds them. This was
tested and passed with an RMSD of 0.45 Å.

**PASS (top-pose RMSD ≤ 2 Å)** means the structure, box and settings are fine.
**FAIL** has a few common causes:
- The ligand SDF uses "ideal" rather than "instance" coordinates (see step 1).
- The box is too small or off-centre.
- An important water, metal or cofactor was removed. A Zn²⁺ site, for example, needs the ion kept. That's an advanced case:
  ask for help, or see the Meeko docs.
- The pocket is flexible (an induced fit). Try `--flexres` on the gatekeeper or loop residues, or a different structure.

You can also check score reproducibility by docking the same ligand with seeds 1, 2 and 3 (`--seed`). Top scores
should agree within about 0.5 kcal/mol at exhaustiveness 32.

---

## 6. Dock

```bat
python %VD% dock --ligands ligands\pdbqt --receptor receptor\receptor.pdbqt --box sites\lig_STI_A201_box.txt --exhaustiveness 32 --outdir docking
```

| Option | Default | Advice |
|---|---|---|
| `--exhaustiveness` | 32 | 8–16 for a first-pass screen of thousands of molecules; re-dock the top hits at 32–64 |
| `--scoring` | vina | `vinardo` is an alternative function; don't mix scores from different functions |
| `--num_modes` | 9 | poses saved per ligand |
| `--seed` | 42 | fixed so the run is reproducible; record it in your methods |
| `--cpu_per_job` / `--jobs` | auto | for big libraries, `--cpu_per_job 1` runs one ligand per core (fastest) |
| `--overwrite` | off | re-runs ligands that are already docked |

**Resuming:** if a run is interrupted (power cut, closed window), run the same command again. Finished ligands are skipped.

**Timing on your 32-thread PC:** imatinib at exhaustiveness 32 took 38 s. A 1,000-compound screen at
exhaustiveness 16 with `--cpu_per_job 1` should take very roughly 1–3 hours. This depends a lot on ligand size and box size.

---

## 7. Results, interactions, visualisation

```bat
python %VD% summarize --dockdir docking --collapse_states
python %VD% interactions --receptor receptor\receptor_clean.pdb --dockdir docking
python %VD% report --receptor receptor\receptor_clean.pdb --dockdir docking --box sites\lig_STI_A201_box.txt --title "ABL1 screen"
start docking\report.html
```

### `docking\results.csv` (open it in Excel)

| Column | Meaning |
|---|---|
| `best_score_kcal_mol` | Vina score of the top pose. More negative is better. |
| `score_pose2` | score of the 2nd pose. If it's close to the best score but the pose is very different, the binding mode is uncertain. |
| `ligand_efficiency` | −score / heavy atoms. ≥ 0.3 is good; it corrects the bias that bigger molecules score better. |
| `rmsd_top_pose_A` | only with `--reference`: distance from a known pose |
| `parent` | with `--collapse_states`: the best protonation state/tautomer per input molecule |

### Interactions

`interactions` uses **ProLIF** to detect H-bonds, hydrophobic contacts, π-stacking, cation-π, salt bridges
(Cationic/Anionic), halogen bonds and metal contacts. It writes:

- `docking\interactions_detailed.csv`: one row per interaction (residue, type, distance, atoms).
- `docking\interactions_matrix.csv`: ligands × (residue:interaction) as 0/1. Sort or filter it in Excel.

**Interaction names are from the ligand's point of view:** `MET318.A:HBAcceptor` means the *ligand accepts* an
H-bond from Met318. `HBDonor` means the ligand donates.

Checking known key interactions across a screen (e.g. a kinase hinge H-bond):

```bat
python %VD% interactions --receptor receptor\receptor_clean.pdb --dockdir docking --key "MET318.A:HBAcceptor,THR315.A,ASP381.A"
```

Each ligand gets a count such as 2/3, plus the list of key contacts it makes. Poses missing an essential interaction
seen in the crystal structure deserve suspicion, even when they score well.

Analyse more than the top pose with `--n_poses 3`.

H-bond detection uses strict geometric cutoffs, so a borderline H-bond (e.g. 3.3 Å or a poor angle) may not be
listed, while PyMOL's looser criteria do show it. That's normal. Judge visually.

### `docking\report.html`: interactive 3D viewer

Double-click it. It opens in any browser, works offline, and can be emailed to colleagues. Click a ligand to see its top
pose, the binding-site residues, labelled interacting residues, H-bond/ionic contacts (dashed lines) and the
box. Checkboxes toggle the pocket surface, labels and so on.

### PyMOL (publication-quality images)

Install once:

```bat
conda install -n vina -c conda-forge pymol-open-source
```

Then open PyMOL with `pymol docking\report.pml`, or start PyMOL and use File → Run Script → `report.pml`.
It loads the receptor, box and all poses, shows the site residues and draws H-bonds (`hbonds`).
Useful PyMOL commands:

```
set state, 2                    # next pose of the selected ligand
distance hb, ligname, receptor, mode=2   # polar contacts for another ligand
ray 2400, 1800; png figure.png, dpi=300  # high-resolution image
```

UCSF ChimeraX (free for academics) also opens `receptor_clean.pdb` and the `*_out.sdf` files directly.

---

## 8. How to read the results (and how not to)

- **Docking scores are not binding affinities.** The typical error against experiment is about 2–3 kcal/mol.
  Use scores to **rank** and **prioritise**, not to quantify binding.
- **Always include controls:** a known binder (positive control), a few known non-binders or decoys, and,
  for a screen, the redocked co-crystal ligand. A hit is interesting if it scores close to or better than the known binder.
- **Look at every top hit.** Check that the pose sits in the pocket, makes the key interactions, has no buried
  unpaired polar groups, and has no strained geometry.
- **Watch the size bias.** Big molecules score better almost automatically. Compare ligand efficiency, too.
- **Watch pose stability.** If poses 1–3 are completely different binding modes with similar scores, the pose is uncertain.
- For decisions that matter, follow up with rescoring, consensus docking (e.g. vina + vinardo), MD simulation, or,
  ultimately, experiment.

**What to report in a paper:** vdock version and tool versions (in every `manifest.json`), PDB ID and chain,
protonation pH, box center and size, exhaustiveness, number of modes, seed, scoring function, redocking RMSD,
and the citations in `README.md`.

---

## 9. Troubleshooting

**Find the error first:**
1. The last lines printed in the terminal usually say what failed.
2. `vdock.log` (if you used `--log`) has the full history.
3. Per-step logs: `receptor\*_meeko.log`, `ligands\scrub.log`, `docking\logs\<ligand>.log`.
4. Add `-v` for verbose output, e.g. `python %VD% -v prep-receptor ...`.
5. `manifest.json` in each output folder records exactly which parameters and versions were used.

| Message / symptom | Cause → fix |
|---|---|
| `O sistema não pode encontrar o caminho especificado` / `No such file` | You're in the wrong folder. `cd` to your project, or use full paths. |
| `'conda' is not recognized` | Use **Miniforge Prompt**, not plain cmd/PowerShell. |
| `ModuleNotFoundError` / `... not found. Activate the environment` | Run `conda activate vina`. |
| `ERROR: receptor file not found` | An earlier step failed. Scroll up to its error. |
| prep-receptor: `template` / `No template matched` / missing atoms | A residue is broken or unusual. Use `--chains` to drop irrelevant chains, `--allow_bad_res` to drop it (if it's outside the site), or `--fill_gaps`. |
| prep-receptor: PDBFixer error on a residue or ligand line | The PDB file is malformed. Download the mmCIF instead (`fetch 1IEP --format cif`) and use `-i input\1IEP.cif`. |
| Many ligands in `prep_failures.tsv` | Usually an invalid SMILES, or a metal-containing compound / organometallic (not supported). The `reason` column says which; check the SMILES in a viewer. |
| Docking log says the ligand is outside the grid box | The box is too small or in the wrong place. Re-check §2. |
| Very positive scores (> 0) | The ligand doesn't fit, is clashing, or the box is wrong. View it in the report. |
| Docking very slow | The box is too large, exhaustiveness is high, or the project is in OneDrive. |
| `sdf export failed` in the dock log | The pose is fine (PDBQT exists); only SDF conversion failed. Usually an unusual ligand chemistry. |
| `interactions`: `Receptor has no hydrogens` | Pass `receptor\receptor_clean.pdb`, not the original PDB. |
| Windows blocks `vina.exe` | It's the official build. Allow it in SmartScreen or your antivirus, or re-run `install-windows.bat`. |
| Validation FAIL | See §5. |

**Starting over:** delete the step's output folder (e.g. `docking\`) and re-run that step. Steps don't
modify their input files.

**Reinstalling:** `conda env remove -n vina`, then run `install-windows.bat` again.

---

## 10. Cheat sheet

```bat
conda activate vina
set VD=C:\docking\vina-docking\vdock.py
python %VD% fetch 1IEP --outdir input
python %VD% sites -i input\1IEP.pdb --chains A --outdir sites
python %VD% prep-receptor -i input\1IEP.pdb --chains A --box sites\lig_STI_A201_box.txt --outdir receptor --name receptor
python %VD% prep-ligands -i ligands.smi --outdir ligands
python %VD% dock --ligands ligands\pdbqt --receptor receptor\receptor.pdbqt --box sites\lig_STI_A201_box.txt
python %VD% summarize --dockdir docking --collapse_states
python %VD% interactions --receptor receptor\receptor_clean.pdb --dockdir docking
python %VD% report --receptor receptor\receptor_clean.pdb --dockdir docking --box sites\lig_STI_A201_box.txt
start docking\report.html
```
