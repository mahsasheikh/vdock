#!/usr/bin/env python
"""
vdock - a reproducible AutoDock Vina 1.2.x docking pipeline.

Follows the current upstream (CCSB / Forli lab) protocol:
  receptor  : PDBFixer clean-up + protonation  ->  Meeko mk_prepare_receptor.py  -> PDBQT (+JSON)
  ligands   : molscrub (protonation states, tautomers, 3D)  ->  Meeko  -> PDBQT
  docking   : AutoDock Vina 1.2.x  (vina | vinardo | ad4 scoring)
  export    : Meeko mk_export.py  -> SDF with correct bond orders / formal charges

Subcommands
  fetch           download a structure from RCSB (PDB/mmCIF)
  box             define the search box (from a reference ligand, residues, or explicit center/size)
  prep-receptor   clean, protonate and convert a receptor to PDBQT
  prep-ligands    protonate / enumerate / embed ligands and convert them to PDBQT
  dock            dock all prepared ligands (parallel, resumable)
  summarize       collect results into a ranked CSV (+ optional RMSD to a reference pose)
  validate        redock a co-crystallised ligand and report pose RMSD (protocol sanity check)
  sites           find binding sites (co-crystal ligands + geometric pocket detection)
  interactions    H-bonds / hydrophobic / pi / ionic / halogen contacts of poses (ProLIF)
  report          interactive HTML viewer report + PyMOL script

Run `python vdock.py <subcommand> -h` for options.
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

__version__ = "1.1.1-safety"
log = logging.getLogger("vdock")


# ----------------------------------------------------------------------------- utils
def setup_logging(verbose: bool, logfile: Path | None = None) -> None:
    handlers: list[logging.Handler] = [logging.StreamHandler(sys.stderr)]
    if logfile:
        logfile.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(logfile, encoding="utf-8"))
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(message)s",
        datefmt="%H:%M:%S",
        handlers=handlers,
        force=True,
    )


def run(cmd: list[str], cwd: Path | None = None, logfile: Path | None = None, check: bool = True):
    log.debug("$ %s", " ".join(map(str, cmd)))
    p = subprocess.run([str(c) for c in cmd], cwd=cwd, capture_output=True, text=True)
    if logfile:
        logfile.write_text(f"$ {' '.join(map(str, cmd))}\n\n--- stdout ---\n{p.stdout}\n--- stderr ---\n{p.stderr}",
                           encoding="utf-8")
    if check and p.returncode != 0:
        raise RuntimeError(f"Command failed ({p.returncode}): {' '.join(map(str, cmd))}\n{p.stderr[-3000:]}")
    return p


HERE = Path(__file__).resolve().parent
ACTIVATE_HINT = "Activate the environment first (micromamba/conda activate vina)."


def require(tool: str) -> str:
    path = shutil.which(tool)
    if not path:
        sys.exit(f"ERROR: '{tool}' not found on PATH. {ACTIVATE_HINT}")
    return path


# Cross-platform tool resolution (Linux / macOS / native Windows).
def vina_exe() -> str:
    """Vina from PATH (conda-forge), else the official binary placed in ./bin by install.ps1."""
    for cand in (shutil.which("vina"), HERE / "bin" / "vina.exe", HERE / "bin" / "vina"):
        if cand and Path(cand).is_file():
            return str(cand)
    sys.exit(f"ERROR: AutoDock Vina executable not found. {ACTIVATE_HINT}")


def meeko_cli(name: str) -> list[str]:
    """Run Meeko CLIs as modules: works on every OS (Windows can't exec *.py scripts)."""
    try:
        __import__(f"meeko.cli.{name}")
    except ImportError:
        sys.exit(f"ERROR: Meeko is not installed. {ACTIVATE_HINT}")
    return [sys.executable, "-m", f"meeko.cli.{name}"]


def scrub_cli() -> list[str]:
    import sysconfig
    for cand in (shutil.which("scrub.py"), Path(sysconfig.get_path("scripts")) / "scrub.py",
                 Path(sys.prefix) / "Scripts" / "scrub.py", Path(sys.prefix) / "bin" / "scrub.py"):
        if cand and Path(cand).is_file():
            return [sys.executable, str(cand)]
    sys.exit(f"ERROR: molscrub (scrub.py) not found. {ACTIVATE_HINT}")


def tool_versions() -> dict:
    v = {"vdock": __version__, "python": sys.version.split()[0]}
    for mod in ("vina", "meeko", "rdkit", "molscrub", "pdbfixer", "openmm", "gemmi", "prody"):
        try:
            m = __import__(mod)
            v[mod] = getattr(m, "__version__", "installed")
        except Exception:
            v[mod] = None
    try:
        v["vina_cli"] = run([vina_exe(), "--version"], check=False).stdout.strip()
    except SystemExit:
        v["vina_cli"] = None
    return v


def write_manifest(path: Path, step: str, params: dict) -> None:
    """Append a provenance record (parameters + tool versions) for every pipeline step."""
    data = []
    if path.exists():
        try:
            data = json.loads(path.read_text())
        except Exception:
            data = []
    data.append({
        "step": step,
        "time": _dt.datetime.now().isoformat(timespec="seconds"),
        "params": {k: (str(v) if isinstance(v, Path) else v) for k, v in params.items() if not callable(v)},
        "versions": tool_versions(),
    })
    path.write_text(json.dumps(data, indent=2))


# ----------------------------------------------------------------------------- fetch
def cmd_fetch(a):
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    fmt = a.format.lower()
    url = f"https://files.rcsb.org/download/{a.pdb_id.upper()}.{fmt}"
    dest = out / f"{a.pdb_id.upper()}.{fmt}"
    log.info("Downloading %s", url)
    urllib.request.urlretrieve(url, dest)
    log.info("Saved %s", dest)


# ----------------------------------------------------------------------------- box
def _coords_from_file(path: Path) -> list[tuple[float, float, float]]:
    suf = path.suffix.lower()
    if suf in (".sdf", ".mol", ".mol2"):
        from rdkit import Chem
        if suf == ".mol2":
            mols = [Chem.MolFromMol2File(str(path), removeHs=False)]
        else:
            mols = [m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None]
        if not mols or mols[0] is None:
            raise ValueError(f"Could not read molecule(s) from {path}")
        xyz = []
        for m in mols:
            xyz += [tuple(p) for p in m.GetConformer().GetPositions()]
        return xyz
    # PDB / PDBQT: plain ATOM/HETATM parsing
    xyz = []
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if not xyz:
        raise ValueError(f"No atoms found in {path}")
    return xyz


def _coords_from_residues(structure: Path, selection: str) -> list[tuple[float, float, float]]:
    """selection like 'A:315,318,B:12' (chain:resnum list, Meeko-style syntax)."""
    wanted, chain = set(), None
    for tok in selection.split(","):
        tok = tok.strip()
        if ":" in tok:
            chain, tok = tok.split(":", 1)
        wanted.add((chain or "", tok))
    import gemmi
    st = gemmi.read_structure(str(structure))
    xyz = []
    for ch in st[0]:
        for res in ch:
            key = (ch.name, f"{res.seqid.num}{res.seqid.icode.strip()}")
            if key in wanted:
                xyz += [(at.pos.x, at.pos.y, at.pos.z) for at in res]
    if not xyz:
        raise ValueError(f"No atoms matched residue selection '{selection}' in {structure}")
    return xyz


def write_box(center, size, out_txt: Path, out_pdb: Path | None = None) -> None:
    out_txt.parent.mkdir(parents=True, exist_ok=True)
    out_txt.write_text(
        "center_x = {:.3f}\ncenter_y = {:.3f}\ncenter_z = {:.3f}\n"
        "size_x = {:.3f}\nsize_y = {:.3f}\nsize_z = {:.3f}\n".format(*center, *size)
    )
    if out_pdb:  # 8 corners, viewable in PyMOL/ChimeraX
        cx, cy, cz = center
        hx, hy, hz = (s / 2 for s in size)
        corners = [(cx + sx * hx, cy + sy * hy, cz + sz * hz) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
        lines = [f"HETATM{i+1:5d}  C   BOX X   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
                 for i, (x, y, z) in enumerate(corners)]
        edges = [(1, 2), (1, 3), (1, 5), (2, 4), (2, 6), (3, 4), (3, 7), (4, 8), (5, 6), (5, 7), (6, 8), (7, 8)]
        lines += [f"CONECT{a:5d}{b:5d}" for a, b in edges]
        out_pdb.write_text("\n".join(lines) + "\nEND\n")


def read_box(path: Path) -> tuple[list[float], list[float]]:
    kv = {}
    for line in Path(path).read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            kv[k.strip()] = float(v.split("#")[0])
    return ([kv["center_x"], kv["center_y"], kv["center_z"]], [kv["size_x"], kv["size_y"], kv["size_z"]])


def cmd_box(a):
    if a.center and a.size:
        center, size = a.center, a.size
    else:
        if a.ligand:
            xyz = _coords_from_file(Path(a.ligand))
            if a.structure:
                # A crystal-instance ligand must be spatially aligned with the receptor.
                # This catches CCD/ideal ligand SDFs, which are often near (0,0,0).
                import gemmi
                import numpy as np
                from scipy.spatial import cKDTree
                st = gemmi.read_structure(str(a.structure))
                protein_xyz = []
                for ch in st[0]:
                    for residue in ch:
                        if not gemmi.find_tabulated_residue(residue.name).is_amino_acid():
                            continue
                        for atom in residue:
                            if atom.element.name not in ("H", "D"):
                                protein_xyz.append((atom.pos.x, atom.pos.y, atom.pos.z))
                if not protein_xyz:
                    raise ValueError(f"No protein heavy atoms found in {a.structure}")
                shortest = float(cKDTree(np.asarray(protein_xyz)).query(np.asarray(xyz))[0].min())
                if shortest > 8.0:
                    raise ValueError(
                        f"COORDINATE-FRAME ERROR: reference ligand is {shortest:.1f} Å from "
                        f"the nearest protein atom in {a.structure}. Do not use an ideal/CCD SDF; "
                        "download the crystal INSTANCE coordinates or choose the co-crystal "
                        "ligand box created by `vdock sites`."
                    )
                log.info("Ligand/protein coordinate check OK: nearest distance %.2f Å", shortest)
            else:
                log.warning("No --structure was provided: cannot verify that the ligand SDF is "
                            "in the receptor crystal coordinate frame. Ideal/CCD SDFs are unsafe for boxes.")
        elif a.residues and a.structure:
            xyz = _coords_from_residues(Path(a.structure), a.residues)
        else:
            sys.exit("Give --center + --size, or --ligand, or --structure + --residues")
        lo = [min(c[i] for c in xyz) for i in range(3)]
        hi = [max(c[i] for c in xyz) for i in range(3)]
        center = [(l + h) / 2 for l, h in zip(lo, hi)]
        size = [max(h - l + 2 * a.padding, a.min_size) for l, h in zip(lo, hi)]
        if a.size:
            size = a.size
    out = Path(a.out)
    write_box(center, size, out, out.with_suffix(".pdb"))
    vol = size[0] * size[1] * size[2]
    log.info("Box center %s  size %s  (%.0f Å^3) -> %s", [round(float(c), 3) for c in center], [round(float(s), 2) for s in size], vol, out)
    if vol > 27000:
        log.warning("Box volume > 30x30x30 Å^3: consider raising --exhaustiveness (search space is large).")


# ----------------------------------------------------------------------------- receptor
def cmd_prep_receptor(a):
    meeko_cli("mk_prepare_receptor")
    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    name = a.name or Path(a.input).stem
    clean_pdb = outdir / f"{name}_clean.pdb"

    # mmCIF: OpenMM would use label_asym_id as chain ID, not the author chain IDs users see
    # ("chain A" on RCSB). Convert to PDB with gemmi, keeping author chain IDs + numbering.
    if Path(a.input).suffix.lower() in (".cif", ".mmcif"):
        import gemmi
        st = gemmi.read_structure(str(a.input))
        st.setup_entities()
        st.remove_alternative_conformations()
        if any(len(ch.name) > 1 for ch in st[0]):
            st.shorten_chain_names()
            log.warning("Multi-letter chain IDs were shortened for PDB format; check chain names in %s", name)
        conv = outdir / f"{name}_input.pdb"
        st.write_pdb(str(conv))
        log.info("Converted mmCIF -> %s (author chain IDs)", conv)
        a.input = str(conv)

    if a.skip_fix:
        shutil.copy(a.input, clean_pdb)
    else:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile

        log.info("PDBFixer: reading %s", a.input)
        fixer = PDBFixer(filename=str(a.input))
        if a.chains:
            keep = set(a.chains.split(","))
            drop = [i for i, c in enumerate(fixer.topology.chains()) if c.id not in keep]
            fixer.removeChains(drop)
        fixer.findMissingResidues()
        if not a.fill_gaps:
            fixer.missingResidues = {}  # do not model missing loops by default
        else:  # model only internal gaps, never termini
            chains = list(fixer.topology.chains())
            fixer.missingResidues = {k: v for k, v in fixer.missingResidues.items()
                                     if 0 < k[1] < len(list(chains[k[0]].residues()))}
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=a.keep_water)
        fixer.findMissingAtoms()
        n_missing = sum(len(v) for v in fixer.missingAtoms.values())
        if n_missing:
            log.info("PDBFixer: adding %d missing heavy atoms", n_missing)
        fixer.addMissingAtoms()
        log.info("PDBFixer: adding hydrogens at pH %.1f", a.ph)
        fixer.addMissingHydrogens(a.ph)
        with open(clean_pdb, "w") as fh:
            PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
        log.info("Clean receptor -> %s", clean_pdb)

    base = outdir / name
    cmd = meeko_cli("mk_prepare_receptor") + ["--read_pdb", clean_pdb, "-o", base, "-p", "-j"]
    if a.delete_residues:
        cmd += ["-d", a.delete_residues]
    if a.flexres:
        cmd += ["-f", a.flexres]
    if a.allow_bad_res:
        cmd += ["--delete_bad_res"]
    if a.box:
        center, size = read_box(Path(a.box))
        cmd += ["--box_center", *map(str, center), "--box_size", *map(str, size), "-v"]
        if a.ad4:
            cmd += ["-g"]
    elif a.ad4:
        sys.exit("--ad4 needs --box (grid maps are computed inside the box)")
    run(cmd, logfile=outdir / f"{name}_meeko.log")
    rigid = base.with_suffix(".pdbqt")
    if not rigid.exists():
        rigid = Path(f"{base}_rigid.pdbqt")
    log.info("Receptor PDBQT -> %s", rigid)

    if a.ad4:
        if not shutil.which("autogrid4"):
            sys.exit("ERROR: autogrid4 not found (no Windows build exists; use --scoring vina/vinardo there).")
        run(["autogrid4","-p", f"{name}.gpf", "-l", f"{name}.glg"], cwd=outdir)
        log.info("AD4 affinity maps written (use: vdock dock --scoring ad4 --maps %s)", base)

    write_manifest(outdir / "manifest.json", "prep-receptor", vars(a))


# ----------------------------------------------------------------------------- ligands
def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "lig"


def cmd_prep_ligands(a):
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    outdir = Path(a.outdir)
    pdbqt_dir = outdir / "pdbqt"
    pdbqt_dir.mkdir(parents=True, exist_ok=True)

    src = Path(a.input)
    if a.skip_scrub:
        sdf = src
        if src.suffix.lower() not in (".sdf", ".mol"):
            sys.exit("--skip_scrub requires an SDF/MOL input with explicit hydrogens and 3D coordinates")
    else:
        sdf = outdir / "scrubbed.sdf"
        cmd = scrub_cli() + [src,"-o", sdf, "--ph", str(a.ph), "--cpu", str(a.cpu),
               "--write_failed_mols", outdir / "scrub_failed.sdf"]
        if a.skip_tautomers:
            cmd.append("--skip_tautomers")
        if a.skip_acidbase:
            cmd.append("--skip_acidbase")
        if a.stereoisomers:
            cmd.append("--do_stereoisomers")
        log.info("molscrub: pH %.1f, enumerating states and 3D-embedding ...", a.ph)
        run(cmd, logfile=outdir / "scrub.log")

    prep = MoleculePreparation(rigid_macrocycles=a.rigid_macrocycles)
    supplier = Chem.SDMolSupplier(str(sdf), removeHs=False)
    seen: dict[str, int] = {}
    ok = fail = 0
    failures = []
    for i, mol in enumerate(supplier):
        if mol is None:
            fail += 1
            failures.append((f"record_{i}", "RDKit could not parse record"))
            continue
        base = _safe(mol.GetProp("_Name") if mol.HasProp("_Name") and mol.GetProp("_Name") else f"lig{i:06d}")
        # protonation states/tautomers from molscrub share a name -> suffix _s1, _s2 ...
        seen[base] = seen.get(base, 0) + 1
        name = base if seen[base] == 1 else f"{base}_s{seen[base]}"
        if seen[base] == 2:  # rename first one for consistency
            first = pdbqt_dir / f"{base}.pdbqt"
            if first.exists():
                first.rename(pdbqt_dir / f"{base}_s1.pdbqt")
        try:
            if mol.GetNumConformers() == 0 or not mol.GetConformer().Is3D():
                raise ValueError("no 3D coordinates")
            if any(at.GetAtomicNum() == 1 for at in mol.GetAtoms()) is False:
                raise ValueError("no explicit hydrogens (run without --skip_scrub)")
            setups = prep.prepare(mol)
            pdbqt, good, err = PDBQTWriterLegacy.write_string(setups[0])
            if not good:
                raise ValueError(err)
            (pdbqt_dir / f"{name}.pdbqt").write_text(pdbqt)
            ok += 1
        except Exception as e:  # keep going on bad molecules, report them
            fail += 1
            failures.append((name, str(e).strip().splitlines()[0][:200]))
    if failures:
        with open(outdir / "prep_failures.tsv", "w") as fh:
            fh.write("name\treason\n")
            for n, r in failures:
                fh.write(f"{n}\t{r}\n")
    log.info("Ligands prepared: %d OK, %d failed -> %s", ok, fail, pdbqt_dir)
    write_manifest(outdir / "manifest.json", "prep-ligands", vars(a))


# ----------------------------------------------------------------------------- docking
_RESULT_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)\s+(\S+)\s+(\S+)")


def parse_vina_out(pdbqt: Path) -> list[float]:
    return [float(m.group(1)) for m in _RESULT_RE.finditer(pdbqt.read_text())]


def _sha256(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for piece in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(piece)
    return h.hexdigest()


def _dock_signature(a, ligs: list[Path], vina_path: str) -> dict:
    # Reuse is only safe when molecular inputs and calculation settings match.
    return {
        "format": 1,
        "vina_executable": str(Path(vina_path).resolve()),
        "vina_version": run([vina_path, "--version"], check=False).stdout.strip(),
        "receptor_sha256": _sha256(a.receptor) if a.receptor else None,
        "box_sha256": _sha256(a.box) if a.box else None,
        "flex_sha256": _sha256(a.flex) if a.flex else None,
        "maps": a.maps,
        "scoring": a.scoring,
        "exhaustiveness": a.exhaustiveness,
        "num_modes": a.num_modes,
        "energy_range": a.energy_range,
        "seed": a.seed,
        "ligands": {p.name: _sha256(p) for p in ligs},
    }


def _dock_one(job: dict) -> dict:
    lig = Path(job["ligand"])
    name = lig.stem
    out_dir = Path(job["outdir"])
    out_pdbqt = out_dir / "poses" / f"{name}_out.pdbqt"
    out_sdf = out_dir / "poses" / f"{name}_out.sdf"
    log_file = out_dir / "logs" / f"{name}.log"
    if out_pdbqt.exists() and not job["overwrite"]:
        sc = parse_vina_out(out_pdbqt)
        if sc and out_sdf.exists() and out_sdf.stat().st_size > 0:
            return {"ligand": name, "status": "skipped (verified files)", "best": sc[0], "n_poses": len(sc)}
        log.warning("Incomplete previous output for %s: re-docking", name)

    cmd = [job["vina"], "--ligand", lig,"--out", out_pdbqt, "--scoring", job["scoring"],
           "--exhaustiveness", job["exhaustiveness"], "--num_modes", job["num_modes"],
           "--energy_range", job["energy_range"], "--seed", job["seed"], "--cpu", job["cpu"],
           "--verbosity", 1]
    if job["scoring"] == "ad4":
        cmd += ["--maps", job["maps"]]
    else:
        cmd += ["--receptor", job["receptor"], "--config", job["box"]]
        if job.get("flex"):
            cmd += ["--flex", job["flex"]]
    try:
        run(cmd, logfile=log_file)
        sc = parse_vina_out(out_pdbqt)
        # Export with correct bond orders (uses the SMILES Meeko stored in the PDBQT)
        exp = subprocess.run(job["mk_export"] + [str(out_pdbqt), "-s", str(out_sdf)], capture_output=True, text=True)
        if not sc:
            raise RuntimeError("Vina finished without parseable scores")
        if exp.returncode or not out_sdf.exists() or out_sdf.stat().st_size == 0:
            raise RuntimeError("SDF export failed: " + exp.stderr[-300:].strip())
        return {"ligand": name, "status": "ok", "best": sc[0], "n_poses": len(sc)}
    except Exception as e:
        return {"ligand": name, "status": f"FAILED: {str(e).splitlines()[0][:150]}", "best": None, "n_poses": 0}


def cmd_dock(a):
    vina_path = vina_exe()
    mk_export = meeko_cli("mk_export")
    lig_dir = Path(a.ligands)
    ligs = sorted(lig_dir.glob("*.pdbqt")) if lig_dir.is_dir() else [lig_dir]
    if not ligs:
        sys.exit(f"No .pdbqt ligands found in {lig_dir}")
    if a.scoring == "ad4" and not a.maps:
        sys.exit("--scoring ad4 needs --maps <receptor basename> (from prep-receptor --ad4)")
    if a.scoring != "ad4" and not (a.receptor and a.box):
        sys.exit("--receptor and --box are required for vina/vinardo scoring")

    for label, f in (("receptor", a.receptor), ("box", a.box), ("flex", a.flex)):
        if f and not Path(f).is_file():
            sys.exit(f"ERROR: {label} file not found: {f}  (did the previous step finish without errors?)")
    if a.maps and not Path(f"{a.maps}.maps.fld").is_file():
        sys.exit(f"ERROR: AD4 maps not found: {a.maps}.maps.fld  (run prep-receptor --ad4)")
    if a.scoring != "ad4":
        import numpy as np
        center, size = read_box(Path(a.box))
        receptor_xyz = np.asarray(_coords_from_file(Path(a.receptor)), dtype=float)
        if receptor_xyz.size == 0:
            sys.exit("ERROR: receptor PDBQT contains no recognizable atom coordinates")
        # Distance to nearest point on the bounding box, not just its center.
        gap = np.maximum(np.abs(receptor_xyz - np.asarray(center)) - np.asarray(size) / 2, 0)
        closest = float(np.sqrt((gap * gap).sum(axis=1)).min())
        if closest > 8.0:
            sys.exit(f"ERROR: docking box is {closest:.1f} Å away from the closest receptor atom. "
                     "Possible ideal-ligand/wrong-coordinate-frame box. "
                     "Compare box.txt with the 'sites' crystal-ligand box.")
        log.info("Receptor/box spatial check OK (closest atom %.2f Å from box)", closest)

    outdir = Path(a.outdir)
    signature_path = outdir / "run_signature.json"
    signature = _dock_signature(a, ligs, vina_path)
    previous_poses = (outdir / "poses").exists() and any((outdir / "poses").glob("*_out.pdbqt"))
    if previous_poses:
        if not signature_path.is_file():
            sys.exit("ERROR: existing poses lack run_signature.json. Use a NEW --outdir "
                     "to avoid mixing different receptors/boxes/settings.")
        try:
            prior = json.loads(signature_path.read_text(encoding="utf-8"))
        except Exception as exc:
            sys.exit(f"ERROR: cannot read existing docking signature: {exc}")
        if prior != signature:
            sys.exit("ERROR: receptor, box, ligand inputs, Vina version, or docking settings "
                     "changed since this output folder was created. Use a NEW --outdir. "
                     "--overwrite cannot safely remove outputs for ligands absent from the new input.")
    (outdir / "poses").mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    # Signature is recorded BEFORE dispatching jobs; interrupted runs can safely resume.
    signature_path.write_text(json.dumps(signature, indent=2), encoding="utf-8")

    ncores = os.cpu_count() or 1
    cpu = a.cpu_per_job or min(ncores, max(1, a.exhaustiveness))
    jobs = a.jobs or max(1, ncores // cpu)
    log.info("Docking %d ligand(s): %d parallel job(s) x %d CPU, exhaustiveness %d, scoring %s, seed %d",
             len(ligs), jobs, cpu, a.exhaustiveness, a.scoring, a.seed)

    base = dict(vina=vina_path, mk_export=mk_export, outdir=str(outdir), receptor=a.receptor, box=a.box, maps=a.maps, flex=a.flex,
                scoring=a.scoring, exhaustiveness=a.exhaustiveness, num_modes=a.num_modes,
                energy_range=a.energy_range, seed=a.seed, cpu=cpu, overwrite=a.overwrite)
    results = []
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        futs = [ex.submit(_dock_one, {**base, "ligand": str(l)}) for l in ligs]
        for n, f in enumerate(as_completed(futs), 1):
            r = f.result()
            results.append(r)
            log.info("[%d/%d] %-30s %s  best=%s", n, len(ligs), r["ligand"], r["status"], r["best"])
    n_fail = sum(r["status"].startswith("FAILED") for r in results)
    write_manifest(outdir / "manifest.json", "dock", {**vars(a), "cpu_per_job_used": cpu, "jobs_used": jobs})
    log.info("Done: %d docked, %d failed. Next: vdock summarize --dockdir %s", len(results) - n_fail, n_fail, outdir)
    if n_fail:
        failures_file = outdir / "docking_failures.tsv"
        with failures_file.open("w", encoding="utf-8") as f:
            f.write("ligand\tstatus\n")
            for r in sorted(results, key=lambda r: r["ligand"]):
                if r["status"].startswith("FAILED"):
                    f.write(f"{r['ligand']}\t{r['status']}\n")
        log.error("DOCKING INCOMPLETE: %d failed. See %s", n_fail, failures_file)
        return 1


# ----------------------------------------------------------------------------- summarize
def _rmsd_to_ref(pose_sdf: Path, ref: Path) -> float | None:
    """Symmetry-aware heavy-atom RMSD *without* re-alignment (in-place, as docking requires)."""
    # Graph matching ignores formal charges, H counts and bond orders, so a pose of a
    # differently protonated/tautomerised state still maps onto the crystal ligand.
    import numpy as np
    from rdkit import Chem

    def skeleton(m):
        rw = Chem.RWMol(Chem.RemoveHs(m, sanitize=False))
        for at in rw.GetAtoms():
            at.SetFormalCharge(0); at.SetIsAromatic(False); at.SetNoImplicit(True); at.SetNumExplicitHs(0)
        for b in rw.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE); b.SetIsAromatic(False)
        return rw.GetMol()

    try:
        if ref.suffix.lower() in (".sdf", ".mol"):
            refm = Chem.MolFromMolFile(str(ref), removeHs=False)
        else:
            refm = Chem.MolFromPDBFile(str(ref), removeHs=False)
        pose = next(m for m in Chem.SDMolSupplier(str(pose_sdf), removeHs=False) if m is not None)
        rs, ps = skeleton(refm), skeleton(pose)
        rx, px = rs.GetConformer().GetPositions(), ps.GetConformer().GetPositions()
        matches = ps.GetSubstructMatches(rs, uniquify=False, useChirality=False, maxMatches=10000)
        if not matches or rs.GetNumAtoms() != ps.GetNumAtoms():
            return None
        return float(min(np.sqrt(((px[list(mt)] - rx) ** 2).sum(1).mean()) for mt in matches))
    except Exception as e:
        log.debug("RMSD failed for %s: %s", pose_sdf, e)
        return None




def _all_pose_rmsds_to_ref(pose_sdf: Path, ref: Path) -> list[float]:
    """Per-pose in-frame heavy-atom RMSDs, symmetry-aware, no superposition."""
    import numpy as np
    from rdkit import Chem

    def skeleton(m):
        rw = Chem.RWMol(Chem.RemoveHs(m, sanitize=False))
        for at in rw.GetAtoms():
            at.SetFormalCharge(0)
            at.SetIsAromatic(False)
            at.SetNoImplicit(True)
            at.SetNumExplicitHs(0)
        for b in rw.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        return rw.GetMol()

    if ref.suffix.lower() in (".sdf", ".mol"):
        refmol = Chem.MolFromMolFile(str(ref), removeHs=False)
    else:
        refmol = Chem.MolFromPDBFile(str(ref), removeHs=False)
    if refmol is None:
        raise ValueError(f"Could not parse reference ligand: {ref}")
    rs = skeleton(refmol)
    rx = rs.GetConformer().GetPositions()
    values = []
    for pose in Chem.SDMolSupplier(str(pose_sdf), removeHs=False):
        if pose is None:
            raise ValueError(f"Unparseable docking pose in {pose_sdf}")
        ps = skeleton(pose)
        if ps.GetNumAtoms() != rs.GetNumAtoms():
            raise ValueError("Reference and docked pose have different heavy-atom counts")
        px = ps.GetConformer().GetPositions()
        matches = ps.GetSubstructMatches(rs, uniquify=False, useChirality=False, maxMatches=10000)
        if not matches:
            raise ValueError("No valid heavy-atom graph mapping for crystal pose")
        values.append(float(min(np.sqrt(((px[list(mt)] - rx) ** 2).sum(1).mean()) for mt in matches)))
    return values


def cmd_summarize(a):
    from rdkit import Chem
    d = Path(a.dockdir)
    rows = []
    for out in sorted((d / "poses").glob("*_out.pdbqt")):
        name = out.name[: -len("_out.pdbqt")]
        sc = parse_vina_out(out)
        if not sc:
            continue
        row = {"ligand": name, "best_score_kcal_mol": sc[0], "n_poses": len(sc),
               "score_pose2": sc[1] if len(sc) > 1 else None}
        sdf = out.with_suffix(".sdf")
        if sdf.exists():
            m = next((m for m in Chem.SDMolSupplier(str(sdf)) if m is not None), None)
            if m is not None:
                row["smiles"] = Chem.MolToSmiles(Chem.RemoveHs(m))
                row["heavy_atoms"] = m.GetNumHeavyAtoms()
                row["ligand_efficiency"] = round(-sc[0] / max(1, m.GetNumHeavyAtoms()), 3)
            if a.reference:
                try:
                    all_rmsds = _all_pose_rmsds_to_ref(sdf, Path(a.reference))
                    if not all_rmsds:
                        raise ValueError("No valid SDF poses")
                    row["rmsd_top_pose_A"] = round(all_rmsds[0], 3)
                    row["rmsd_best_saved_pose_A"] = round(min(all_rmsds), 3)
                    row["best_rmsd_pose_rank"] = 1 + min(range(len(all_rmsds)), key=all_rmsds.__getitem__)
                except (ValueError, OSError, StopIteration) as exc:
                    log.warning("RMSD unavailable for %s: %s", name, exc)
        rows.append(row)
    if not rows:
        sys.exit(f"No results in {d/'poses'}")
    # collapse protonation/tautomer states (_sN suffix) to the best-scoring one per parent
    if a.collapse_states:
        best = {}
        for r in rows:
            parent = re.sub(r"_s\d+$", "", r["ligand"])
            if parent not in best or r["best_score_kcal_mol"] < best[parent]["best_score_kcal_mol"]:
                best[parent] = {**r, "parent": parent}
        rows = list(best.values())
    rows.sort(key=lambda r: r["best_score_kcal_mol"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    cols = ["rank", "parent", "ligand", "best_score_kcal_mol", "score_pose2", "n_poses", "heavy_atoms",
            "ligand_efficiency", "rmsd_top_pose_A", "rmsd_best_saved_pose_A",
            "best_rmsd_pose_rank", "smiles"]
    cols = [c for c in cols if any(c in r for r in rows)]
    out_csv = Path(a.out or d / "results.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    log.info("Wrote %d rows -> %s", len(rows), out_csv)
    for r in rows[: a.top]:
        extra = f"  RMSD={r['rmsd_top_pose_A']} Å" if r.get("rmsd_top_pose_A") is not None else ""
        print(f"{r['rank']:>4}  {r['ligand']:<32} {r['best_score_kcal_mol']:>7.2f} kcal/mol{extra}")


# ----------------------------------------------------------------------------- validate
def cmd_validate(a):
    """Redock a reference ligand into its own receptor; report score + RMSD of the top pose."""
    work = Path(a.outdir)
    work.mkdir(parents=True, exist_ok=True)
    ns = argparse.Namespace
    box = work / "box.txt"
    cmd_box(ns(center=None, size=None, ligand=a.ligand, residues=None, structure=a.receptor,
               padding=a.padding, min_size=a.min_box, out=str(box)))
    cmd_prep_receptor(ns(input=a.receptor, outdir=str(work / "receptor"), name="receptor", chains=a.chains,
                         ph=a.ph, keep_water=False, fill_gaps=False, skip_fix=a.skip_fix, delete_residues=None,
                         flexres=None, allow_bad_res=a.allow_bad_res, box=str(box), ad4=False))
    cmd_prep_ligands(ns(input=a.ligand, outdir=str(work / "ligand"), ph=a.ph, cpu=1, skip_scrub=a.skip_scrub,
                        skip_tautomers=True, skip_acidbase=True, stereoisomers=False, rigid_macrocycles=False))
    rec = next((work / "receptor").glob("receptor*.pdbqt"))
    cmd_dock(ns(ligands=str(work / "ligand" / "pdbqt"), receptor=str(rec), box=str(box), maps=None, flex=None,
                scoring=a.scoring, exhaustiveness=a.exhaustiveness, num_modes=9, energy_range=3.0,
                seed=a.seed, cpu_per_job=None, jobs=1, overwrite=True, outdir=str(work / "docking")))
    cmd_summarize(ns(dockdir=str(work / "docking"), reference=a.ligand, collapse_states=False, out=None, top=5))
    res = list(csv.DictReader(open(work / "docking" / "results.csv")))
    valid = [r for r in res if r.get("rmsd_top_pose_A") and r.get("best_score_kcal_mol")]
    if not valid:
        print("\nVALIDATION INCOMPLETE: no calculable crystal RMSD for the docked ligand")
        return 1
    # Select the globally BEST-SCORED state/pose. Do not cherry-pick the smallest RMSD.
    top = min(valid, key=lambda r: float(r["best_score_kcal_mol"]))
    top_rmsd = float(top["rmsd_top_pose_A"])
    best_saved = min(float(r.get("rmsd_best_saved_pose_A") or float("inf")) for r in valid)
    verdict = "PASS" if top_rmsd <= 2.0 else "FAIL"
    print(f"\nVALIDATION {verdict}: globally top-scoring pose RMSD = {top_rmsd:.3f} Å (≤ 2.0 Å)")
    print(f"BEST-OF-SAVED-POSES RMSD = {best_saved:.3f} Å (different metric)")
    if top_rmsd > 2.0 and best_saved <= 2.0:
        print("POSE RECOVERY SUCCEEDED AMONG SAVED POSES, BUT TOP-1 RANKING FAILED.")
        print("Do NOT rerank screening compounds by crystal RMSD; unknown compounds have no crystal reference.")
    return 0 if verdict == "PASS" else 1


# ----------------------------------------------------------------------------- CLI
def main(argv=None):
    p = argparse.ArgumentParser(prog="vdock", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-V", "--version", action="version", version=f"vdock {__version__}")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--log", help="also write log to this file")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("fetch", help="download a PDB entry from RCSB")
    s.add_argument("pdb_id")
    s.add_argument("--format", default="pdb", choices=["pdb", "cif"])
    s.add_argument("--outdir", default="input")
    s.set_defaults(func=cmd_fetch)

    s = sub.add_parser("box", help="define the docking search box")
    s.add_argument("--ligand", help="reference ligand (sdf/mol/mol2/pdb/pdbqt): box envelops it")
    s.add_argument("--structure", help="PDB/mmCIF used with --residues; with --ligand, checks matching coordinate frame")
    s.add_argument("--residues", help="binding-site residues, e.g. 'A:315,318,381'")
    s.add_argument("--center", type=float, nargs=3, metavar=("X", "Y", "Z"))
    s.add_argument("--size", type=float, nargs=3, metavar=("SX", "SY", "SZ"))
    s.add_argument("--padding", type=float, default=5.0, help="Å added on each side of the ligand/residues [5]")
    s.add_argument("--min_size", type=float, default=18.0, help="minimum edge length in Å [18]")
    s.add_argument("--out", default="box.txt")
    s.set_defaults(func=cmd_box)

    s = sub.add_parser("prep-receptor", help="clean + protonate + PDBQT")
    s.add_argument("-i", "--input", required=True, help="receptor PDB/mmCIF")
    s.add_argument("--outdir", default="receptor")
    s.add_argument("--name", help="output basename [input stem]")
    s.add_argument("--chains", help="keep only these chains, e.g. 'A' or 'A,B'")
    s.add_argument("--ph", type=float, default=7.4)
    s.add_argument("--keep_water", action="store_true", help="keep crystallographic waters")
    s.add_argument("--fill_gaps", action="store_true", help="model missing internal loops with PDBFixer")
    s.add_argument("--skip_fix", action="store_true", help="input is already clean and fully protonated")
    s.add_argument("--delete_residues", help="Meeko -d selection, e.g. 'A:350'")
    s.add_argument("--flexres", help="flexible side chains, e.g. 'A:315,B:23'")
    s.add_argument("--allow_bad_res", action="store_true", help="drop residues that fail template matching")
    s.add_argument("--box", help="box.txt from `vdock box` (also writes Vina box config + box PDB)")
    s.add_argument("--ad4", action="store_true", help="also write GPF and run autogrid4 for AD4 scoring")
    s.set_defaults(func=cmd_prep_receptor)

    s = sub.add_parser("prep-ligands", help="molscrub + Meeko -> PDBQT")
    s.add_argument("-i", "--input", required=True, help=".smi/.sdf/.csv/.xlsx or a SMILES string")
    s.add_argument("--outdir", default="ligands")
    s.add_argument("--ph", type=float, default=7.4)
    s.add_argument("--cpu", type=int, default=os.cpu_count() or 1)
    s.add_argument("--skip_scrub", action="store_true", help="input SDF is already protonated + 3D")
    s.add_argument("--skip_tautomers", action="store_true")
    s.add_argument("--skip_acidbase", action="store_true")
    s.add_argument("--stereoisomers", action="store_true", help="enumerate undefined stereocentres")
    s.add_argument("--rigid_macrocycles", action="store_true")
    s.set_defaults(func=cmd_prep_ligands)

    s = sub.add_parser("dock", help="run AutoDock Vina on all prepared ligands")
    s.add_argument("--ligands", required=True, help="directory of .pdbqt (or single file)")
    s.add_argument("--receptor", help="receptor .pdbqt (rigid part if flexible)")
    s.add_argument("--flex", help="flexible side-chain .pdbqt")
    s.add_argument("--box", help="box.txt")
    s.add_argument("--maps", help="AD4 map basename (with --scoring ad4)")
    s.add_argument("--scoring", default="vina", choices=["vina", "vinardo", "ad4"])
    s.add_argument("--exhaustiveness", type=int, default=32, help="[32; upstream recommends ≥32 for reliable poses]")
    s.add_argument("--num_modes", type=int, default=9)
    s.add_argument("--energy_range", type=float, default=3.0)
    s.add_argument("--seed", type=int, default=42, help="fixed seed for reproducibility [42]")
    s.add_argument("--cpu_per_job", type=int, help="threads per Vina run [min(cores, exhaustiveness)]")
    s.add_argument("--jobs", type=int, help="parallel Vina runs [cores // cpu_per_job]")
    s.add_argument("--overwrite", action="store_true", help="redo ligands already docked")
    s.add_argument("--outdir", default="docking")
    s.set_defaults(func=cmd_dock)

    s = sub.add_parser("summarize", help="ranked CSV of results")
    s.add_argument("--dockdir", default="docking")
    s.add_argument("--reference", help="reference ligand SDF for RMSD (redocking)")
    s.add_argument("--collapse_states", action="store_true", help="keep best state per parent molecule")
    s.add_argument("--out", help="CSV path [<dockdir>/results.csv]")
    s.add_argument("--top", type=int, default=20)
    s.set_defaults(func=cmd_summarize)

    s = sub.add_parser("validate", help="redocking sanity check (RMSD ≤ 2 Å)")
    s.add_argument("--receptor", required=True, help="receptor PDB (ligand removed or not; hetero atoms are stripped)")
    s.add_argument("--ligand", required=True, help="crystal ligand SDF with correct bond orders")
    s.add_argument("--chains")
    s.add_argument("--ph", type=float, default=7.4)
    s.add_argument("--padding", type=float, default=5.0)
    s.add_argument("--min_box", type=float, default=18.0)
    s.add_argument("--scoring", default="vina", choices=["vina", "vinardo"])
    s.add_argument("--exhaustiveness", type=int, default=32)
    s.add_argument("--seed", type=int, default=42)
    s.add_argument("--skip_fix", action="store_true")
    s.add_argument("--skip_scrub", action="store_true", help="use crystal ligand protonation as given")
    s.add_argument("--allow_bad_res", action="store_true")
    s.add_argument("--outdir", default="validation")
    s.set_defaults(func=cmd_validate)

    # ---- analysis tools (vdock_analysis.py), imported lazily
    def _lazy(fn):
        def call(args):
            import vdock_analysis
            return getattr(vdock_analysis, fn)(args)
        return call

    s = sub.add_parser("sites", help="find binding sites: co-crystal ligands + geometric pockets")
    s.add_argument("-i", "--structure", required=True, help="PDB/mmCIF (raw from RCSB is fine)")
    s.add_argument("--chains", help="restrict to chains, e.g. 'A'")
    s.add_argument("--outdir", default="sites")
    s.add_argument("--padding", type=float, default=5.0)
    s.add_argument("--min_size", type=float, default=18.0)
    s.add_argument("--contact_cutoff", type=float, default=4.0, help="Å for lining residues [4]")
    s.add_argument("--max_pockets", type=int, default=5)
    s.add_argument("--min_buried", type=float, default=0.70, help="pocket enclosure threshold 0-1 [0.70]")
    s.add_argument("--all_hetero", action="store_true", help="also report ions and crystallisation additives")
    s.set_defaults(func=_lazy("cmd_sites"))

    s = sub.add_parser("interactions", help="protein-ligand interactions of docked poses (ProLIF)")
    s.add_argument("--receptor", required=True, help="receptor/<name>_clean.pdb from prep-receptor (has H)")
    s.add_argument("--dockdir", default="docking")
    s.add_argument("--poses", help="a single *_out.sdf instead of the whole dockdir")
    s.add_argument("--ligands", help="comma-separated ligand names to analyse (default: all)")
    s.add_argument("--n_poses", type=int, default=1, help="poses per ligand to analyse [1 = top pose]")
    s.add_argument("--key", help="key contacts to check, e.g. 'MET318.A:HBDonor,THR315.A'")
    s.add_argument("--top", type=int, default=15)
    s.add_argument("--outdir", help="[dockdir]")
    s.set_defaults(func=_lazy("cmd_interactions"))

    s = sub.add_parser("report", help="interactive HTML report (+ PyMOL script)")
    s.add_argument("--receptor", required=True, help="receptor/<name>_clean.pdb")
    s.add_argument("--dockdir", default="docking")
    s.add_argument("--box", help="box.txt (drawn in the viewer)")
    s.add_argument("--top", type=int, default=50, help="ligands to include [50]")
    s.add_argument("--title", default="Docking results")
    s.add_argument("--out", help="[<dockdir>/report.html]")
    s.set_defaults(func=_lazy("cmd_report"))

    for stream in (sys.stdout, sys.stderr):  # Windows consoles default to cp1252 (Å, ≤ ...)
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    a = p.parse_args(argv)
    setup_logging(a.verbose, Path(a.log) if a.log else None)
    rc = a.func(a)
    return rc or 0


if __name__ == "__main__":
    sys.exit(main())
