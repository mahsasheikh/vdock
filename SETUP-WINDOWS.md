# Setting up vdock on Windows 10/11

Total time: about 15–20 minutes, mostly downloads. You need an internet connection and about 3 GB of disk space.

## 1. Unpack

Right-click `vina-docking.tar.gz`, choose **Extract All**, and put the folder somewhere without
spaces or accents in the path, e.g. `C:\docking\vina-docking`.
(Older Windows versions without built-in .tar.gz support can use 7-Zip, or run
`tar -xzf vina-docking.tar.gz` in a terminal.)

## 2. Install

Double-click **`install-windows.bat`** in that folder. It:

1. finds Miniforge (conda-forge's Python distribution), or installs it with `winget` if it's missing.
   If `winget` isn't available, it tells you to download Miniforge from https://conda-forge.org/download/.
   Install it with the default options, then run the installer again;
2. creates a conda environment called `vina` with Meeko 0.8.0, molscrub, RDKit, PDBFixer and OpenMM;
3. downloads the official **AutoDock Vina 1.2.7 for Windows** from the Vina GitHub release into
   `bin\vina.exe`;
4. checks that everything works.

The run ends with **"Installed successfully."**

> If Windows SmartScreen or your antivirus asks about `vina.exe`, it's the official build from
> github.com/ccsb-scripps/AutoDock-Vina. Allow it.

## 3. Validate (do this once)

Open **Miniforge Prompt** from the Start menu and run:

```bat
conda activate vina
cd C:\docking\vina-docking
examples\run_1iep_redock.bat
```

After about 3–10 minutes the last line should read:

```
VALIDATION PASS: top-pose RMSD to crystal = 0.2x Å (criterion ≤ 2.0 Å)
```

## 4. Your own projects

Always work in **Miniforge Prompt** after `conda activate vina`. The commands are the same as on Linux;
see `README.md`:

```bat
python vdock.py fetch 1IEP --outdir myproj\input
python vdock.py box --structure myproj\input\1IEP.pdb --residues A:271,286,315,318,381 --out myproj\box.txt
python vdock.py prep-receptor -i myproj\input\1IEP.pdb --chains A --box myproj\box.txt --outdir myproj\receptor --name receptor
python vdock.py prep-ligands -i library.smi --outdir myproj\ligands
python vdock.py dock --ligands myproj\ligands\pdbqt --receptor myproj\receptor\receptor.pdbqt --box myproj\box.txt
python vdock.py summarize --dockdir docking --collapse_states
```

Or copy `examples\run_project_template.bat`, edit the SETTINGS at the top, and run it in stages. **GUIDE.md** walks through everything.

**Using all your CPU cores:** by default `dock` already splits the work across all logical cores.
For large libraries (hundreds of ligands or more), `--cpu_per_job 1` is usually fastest because each core docks its own ligand.
Results open in PyMOL or ChimeraX from `docking\poses\*.sdf`, together with `receptor\receptor_clean.pdb`
and `box.pdb`.

## Windows limitations

- **AD4 scoring (`--scoring ad4`) is unavailable** because AutoGrid has no Windows build. Vina and Vinardo
  scoring both work. If you need AD4, use WSL2 (below).
- Avoid spaces and non-ASCII characters in folder and file names. Some of the chemistry tools handle them badly.

## Alternative: WSL2 (identical to the Linux setup)

If you'd rather run the exact Linux environment (including AD4), open PowerShell as Administrator
and run `wsl --install`. Reboot, open **Ubuntu**, then:

```bash
cd /mnt/c/docking/vina-docking
bash install.sh
micromamba activate vina
bash examples/run_1iep_redock.sh
```

WSL2 uses all of your PC's cores at near-native speed.

## Troubleshooting

| Symptom | Fix |
|---|---|
| `conda` not recognised | Use **Miniforge Prompt**, not plain cmd/PowerShell. |
| `ERROR: AutoDock Vina executable not found` | `bin\vina.exe` is missing: re-run `install-windows.bat`. |
| `ModuleNotFoundError: meeko` | You forgot `conda activate vina`. |
| Script blocked by execution policy | Use `install-windows.bat` (it bypasses the policy just for this script). |
| Receptor prep fails on a residue | Check the PDB for broken residues/ligands; try `--allow_bad_res` or `--chains A`. |
