#!/usr/bin/env python3
r"""Automatically repair both intra- and inter-residue Meeko connectivity errors.

Run in your existing vdock vina environment. A separate repaired PDB is saved
only after a final whole-receptor connectivity scan succeeds (unless you
explicitly choose --allow-partial). Neither input is modified.

Windows CMD, from vdock's root:
    python relaxer_auto.py --original receptor\5I73\5I73.pdb --clean receptor\5I73\5I73_clean.pdb

For intra-residue errors, one residue moves at a time. For inter-residue
errors, the offending local residue cluster can move together. Newly added
heavy atoms and hydrogens are tried first, then restrained side chains.
Observed backbone heavy atoms stay fixed. Uses the installed Meeko function
for inter-residue bond perception. This is geometry repair, not MD or a
guarantee of all Meeko template matches.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import math
import re
import sys


@dataclass(frozen=True)
class PDBAtom:
    line_index: int
    serial: str
    name: str
    residue: str
    chain: str
    number: str
    insertion: str
    altloc: str
    element: str
    xyz: tuple[float, float, float]
    occupancy: float

    @property
    def key(self):
        return self.chain, self.number, self.insertion

    @property
    def label(self):
        return f"{self.name} (serial {self.serial})"

    @property
    def hydrogen(self):
        return self.element in {"H", "D"}


BACKBONE = {"N", "CA", "C", "O", "OXT", "OT1", "OT2"}


def say(message=""):
    print(message, flush=True)


def read_pdb(path: Path):
    """Read the first (and only, for a prepared input) PDB coordinate model."""
    lines = path.read_text(encoding="ascii", errors="replace").splitlines()
    atoms = []
    models = sum(line.startswith("MODEL ") for line in lines)
    if models > 1:
        raise ValueError(f"{path} has multiple models; supply a single-model PDB.")
    for i, line in enumerate(lines):
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        if len(line) < 54:
            raise ValueError(f"Malformed coordinate record at {path}:{i+1}")
        name = line[12:16].strip()
        element = line[76:78].strip().upper()
        # PDBFixer/OpenMM writes elements. Fallback is for ordinary protein atoms.
        if not element:
            element = re.sub(r"^[0-9]+", "", name)[0].upper()
        atoms.append(PDBAtom(
            i, line[6:11].strip(), name, line[17:20].strip(),
            line[21:22].strip(), line[22:26].strip(), line[26:27].strip(),
            line[16:17].strip(), element,
            (float(line[30:38]), float(line[38:46]), float(line[46:54])),
            float(line[54:60].strip() or 0.0),
        ))
    if not atoms:
        raise ValueError(f"No coordinate records in {path}")
    return lines, atoms


def parse_residue(selection):
    m = re.fullmatch(r"([^:]*):(-?\d+)([A-Za-z]?)", selection.strip())
    if not m:
        raise ValueError("Use --residue CHAIN:NUMBER, e.g. A:79 or A:79B.")
    return m.group(1), m.group(2), m.group(3)


def select_atoms(atoms, key):
    return [(i, a) for i, a in enumerate(atoms) if a.key == key]


def original_by_name(atoms):
    """Keep highest-occupancy original record per atom name."""
    groups = defaultdict(list)
    for atom in atoms:
        groups[atom.name].append(atom)
    return {name: max(group, key=lambda a: (a.occupancy, not a.altloc))
            for name, group in groups.items()}


def dist(p, q):
    return math.dist(p, q)


def pair(a, b):
    return tuple(sorted((a, b)))


def inferred_graph(atoms, positions, Chem, rdDetermineBonds):
    """Reproduce Meeko's distance-only within-residue RDKit bond inference."""
    rw = Chem.RWMol()
    for atom in atoms:
        symbol = atom.element.capitalize()
        if symbol == "D":
            symbol = "H"
        rd_atom = Chem.Atom(symbol)
        if atom.element == "D":
            rd_atom.SetIsotope(2)
        rw.AddAtom(rd_atom)
    mol = rw.GetMol()
    conf = Chem.Conformer(len(atoms))
    conf.Set3D(True)
    for i, xyz in enumerate(positions):
        conf.SetAtomPosition(i, tuple(float(x) for x in xyz))
    mol.AddConformer(conf, assignId=True)
    # This is the same connectivity function used by Meeko's
    # _aux_altloc_mol_build, not a hand-picked distance threshold.
    rdDetermineBonds.DetermineConnectivity(mol)
    edges = {pair(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
             for b in mol.GetBonds()}
    degrees = [len(a.GetNeighbors()) for a in mol.GetAtoms()]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and len(atom.GetNeighbors()) == 4:
            atom.SetFormalCharge(1)
    try:
        Chem.SanitizeMol(mol)
        error = None
    except Exception as exc:
        error = str(exc)
    return edges, degrees, error


def atom_origin(atom, originals):
    if atom.name in originals:
        return "present in original"
    if atom.hydrogen:
        return "absent in original (H added during preparation)"
    return "absent in original (likely rebuilt heavy atom)"


def report_graph(title, atoms, positions, expected, originals, Chem, rdDetermineBonds):
    say(f"\n--- {title} ---")
    inferred, degrees, error = inferred_graph(atoms, positions, Chem, rdDetermineBonds)
    extra = inferred - expected
    missing = expected - inferred
    say(f"RDKit connectivity: {len(inferred)} bonds; "
        f"standard residue topology: {len(expected)} bonds")
    say(f"Extra inferred bonds: {len(extra)}; missing expected bonds: {len(missing)}")
    say(f"Meeko-style RDKit sanitization: {'PASS' if error is None else 'FAIL: ' + error}")
    expected_degree = [sum(i in edge for edge in expected) for i in range(len(atoms))]
    suspects = {i for edge in extra | missing for i in edge}
    m = re.search(r"atom\s*#\s*(\d+)", error or "", flags=re.I)
    if m and int(m.group(1)) < len(atoms):
        suspects.add(int(m.group(1)))
    for i in sorted(suspects):
        atom = atoms[i]
        say(f"  atom #{i}: {atom.label}, {atom.element}, "
            f"inferred neighbors={degrees[i]}, expected={expected_degree[i]}; "
            f"{atom_origin(atom, originals)}")
    for label, edges in (("EXTRA", extra), ("MISSING", missing)):
        for i, j in sorted(edges):
            ai, aj = atoms[i], atoms[j]
            say(f"  {label}: #{i} {ai.label} -- #{j} {aj.label}: "
                f"{dist(positions[i], positions[j]):.3f} A")
            say(f"    {ai.name}: {atom_origin(ai, originals)}; "
                f"{aj.name}: {atom_origin(aj, originals)}")
            if ai.name in originals and aj.name in originals:
                org = dist(originals[ai.name].xyz, originals[aj.name].xyz)
                say(f"    Original heavy/atom distance: {org:.3f} A "
                    "(highest-occupancy records)")
    if error is None and not extra and not missing:
        say("  Geometry-based connectivity agrees with the standard residue topology.")
    return inferred, extra, missing, error


def report_close_contacts(all_atoms, all_xyz, selected, topo_bonds, np):
    """Report unusually close nonbonded neighbors; no decision based on cutoff alone."""
    changes = []
    selected_indices = {i for i, _ in selected}
    xyz = np.asarray(all_xyz, dtype=float)
    for i, a in selected:
        delta = xyz - xyz[i]
        dd = np.sqrt(np.sum(delta * delta, axis=1))
        for j in np.where(dd < 1.55)[0]:
            j = int(j)
            if j == i or pair(i, j) in topo_bonds:
                continue
            # Avoid duplicating intra-residue pairs.
            if j in selected_indices and j < i:
                continue
            b = all_atoms[j]
            # A short hydrogen-heavy contact can sometimes be legitimate;
            # report it for inspection rather than pruning atoms.
            threshold = 1.15 if a.hydrogen or b.hydrogen else 1.55
            if dd[j] < threshold:
                changes.append((float(dd[j]), a, b))
    changes.sort(key=lambda v: v[0])
    if changes:
        say(f"  Short nonbonded contacts involving the selected residue: {len(changes)}")
        for d, a, b in changes[:15]:
            say(f"    {a.name} -- {b.chain}:{b.number}{b.insertion} "
                f"{b.residue} {b.name}: {d:.3f} A")
        if len(changes) > 15:
            say("    (first 15 shown)")
    else:
        say("  No very short nonbonded contacts detected by the diagnostic cutoffs.")
    return changes


def update_coordinates(lines, atoms, new_xyz, mobile_indices):
    """Preserve all original lines, replacing only mobile atoms' XYZ fields."""
    out = list(lines)
    for i in mobile_indices:
        x, y, z = (float(v) for v in new_xyz[i])
        if not all(math.isfinite(v) for v in (x, y, z)):
            raise ValueError("Non-finite minimized coordinates.")
        coord = f"{x:8.3f}{y:8.3f}{z:8.3f}"
        if len(coord) != 24:
            raise ValueError("Coordinates outside PDB field width.")
        line_no = atoms[i].line_index
        line = out[line_no]
        out[line_no] = line[:30] + coord + line[54:]
    return out


def minimize(pdb, initial_xyz, atoms, selected, originals, mobile, args,
             mm, app, unit, np, cross_targets=(), original_maps=None):
    """Parameterize full protein, set other atom masses to zero, minimize."""
    # No constraints are used: zero-mass particles must not be coupled by
    # constraints to mobile particles. OpenMM uses chemical topology, not the
    # spurious distance-inferred bonds reported by RDKit.
    ff = app.ForceField(args.forcefield)
    system = ff.createSystem(
        pdb.topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.2 * unit.nanometer,
        constraints=None,
        rigidWater=False,
        removeCMMotion=False,
    )
    if system.getNumParticles() != len(atoms):
        raise ValueError("Force field generated extra particles; refusing ambiguous atom mapping.")
    mobile = set(mobile)
    for i in range(system.getNumParticles()):
        if i not in mobile:
            system.setParticleMass(i, 0 * unit.dalton)

    # Prefer the original heavy-atom geometry if a fallback must move
    # experimentally observed side-chain atoms.
    restraint = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", args.restraint)
    for p in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(p)
    def original_for(atom):
        return originals if original_maps is None else original_maps[atom.key]

    count = 0
    for i, atom in selected:
        if i in mobile and not atom.hydrogen and atom.name in original_for(atom):
            x0, y0, z0 = (v / 10.0 for v in initial_xyz[i])
            restraint.addParticle(i, [x0, y0, z0])
            count += 1
    if count:
        system.addForce(restraint)
    # A one-sided restraint on known false bonds helps avoid residual
    # geometry-based covalent contacts, including FF-excluded nearby pairs.
    if cross_targets:
        repel = mm.CustomBondForce("0.5*cross_k*step(r0-r)*(r0-r)^2")
        repel.addGlobalParameter("cross_k", args.cross_restraint)
        repel.addPerBondParameter("r0")
        for i, j, target_angstrom in cross_targets:
            if i in mobile or j in mobile:
                repel.addBond(i, j, [target_angstrom / 10.0])
        if repel.getNumBonds():
            system.addForce(repel)
            say(f"  Local false-bond separation restraints: {repel.getNumBonds()}")
    say(f"  Mobile atoms: {len(mobile)}; fixed atoms: {len(atoms)-len(mobile)}")
    say(f"  Restrained originally observed heavy atoms: {count}")
    for i, atom in selected:
        if i in mobile:
            say(f"    {atom.chain}:{atom.number}{atom.insertion} {atom.name}: {atom_origin(atom, original_for(atom))}")
    say(f"  Force field: {args.forcefield}; nonperiodic cutoff 12 A")
    say("  Note: this is a local geometry repair, not an equilibrated MD simulation.")

    platform = mm.Platform.getPlatformByName("CPU")
    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(system, integrator, platform, {"Threads": str(args.threads)})
    try:
        context.setPositions(np.asarray(initial_xyz) * unit.angstrom)
        before = context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            unit.kilojoule_per_mole)
        say(f"  Starting potential energy: {before:.3f} kJ/mol")
        mm.LocalEnergyMinimizer.minimize(
            context, args.tolerance * unit.kilojoule_per_mole / unit.nanometer,
            args.max_iterations,
        )
        state = context.getState(getPositions=True, getEnergy=True, getForces=True)
        after = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        result = np.array(state.getPositions(asNumpy=True).value_in_unit(unit.angstrom),
                          dtype=float)
        forces = np.array(state.getForces(asNumpy=True).value_in_unit(
            unit.kilojoule_per_mole / unit.nanometer), dtype=float)
        if not np.isfinite(result).all() or not math.isfinite(after) or not np.isfinite(forces).all():
            raise ValueError("Minimization produced non-finite values.")
        say(f"  Final potential energy:    {after:.3f} kJ/mol")
        say(f"  Mobile RMS force/component: "
            f"{np.sqrt(np.mean(forces[sorted(mobile)]**2)):.2f} kJ/mol/nm")
        fixed = sorted(set(range(len(atoms))) - mobile)
        if fixed:
            move = float(np.max(np.linalg.norm(result[fixed] - initial_xyz[fixed], axis=1)))
            say(f"  Largest fixed-atom shift: {move:.6f} A")
            if move > 1e-4:
                raise ValueError("A supposedly fixed atom moved; not writing output.")
        # Match what will actually be read back from a PDB with 3 decimal digits.
        rounded = np.array(initial_xyz, copy=True)
        for i in mobile:
            rounded[i] = np.array([float(f"{v:.3f}") for v in result[i]])
        for i, atom in selected:
            delta = dist(rounded[i], initial_xyz[i])
            if i in mobile:
                say(f"    displacement {atom.name}: {delta:.3f} A")
        return rounded
    finally:
        del context, integrator


@dataclass
class Residue:
    key: tuple
    selected: list
    atoms: list
    ids: list
    expected: set
    originals: dict
    unsupported: str | None = None

    @property
    def label(self):
        chain, number, insertion = self.key
        return f"{self.atoms[0].residue} {chain}:{number}{insertion}"


def make_residues(clean, raw, topo_bonds):
    """Index residues and their name-based OpenMM topology once."""
    clean_groups = defaultdict(list)
    raw_groups = defaultdict(list)
    for i, atom in enumerate(clean):
        clean_groups[atom.key].append((i, atom))
    for atom in raw:
        raw_groups[atom.key].append(atom)

    internal_bonds = defaultdict(set)
    for i, j in topo_bonds:
        if clean[i].key == clean[j].key:
            internal_bonds[clean[i].key].add((i, j))

    residues = []
    for key, selected in clean_groups.items():
        atoms = [a for _, a in selected]
        ids = [i for i, _ in selected]
        local = {i: k for k, i in enumerate(ids)}
        expected = {pair(local[i], local[j])
                    for i, j in internal_bonds[key]}
        originals = original_by_name(raw_groups[key])
        unsupported = None
        if len({a.residue for a in atoms}) != 1:
            unsupported = "ambiguous residue identity"
        elif len({a.name for a in atoms}) != len(atoms):
            unsupported = "duplicate atom names"
        elif any(a.altloc for a in atoms):
            unsupported = "unresolved alternate conformations"
        elif not expected:
            unsupported = "no standard intra-residue topology; may be ion or heterogen"
        elif any(a.hydrogen and not any(i in edge for edge in expected)
                 for i, a in enumerate(atoms)):
            unsupported = "OpenMM could not identify one or more hydrogen bonds"
        elif not any(a.hydrogen for a in atoms):
            unsupported = "residue has no explicit hydrogens"
        residues.append(Residue(key, selected, atoms, ids, expected,
                                originals, unsupported))
    return residues


def scan_residues(residues, xyz, Chem, rdDetermineBonds):
    """Scan every eligible residue; distinguish sanitization failures from graph differences."""
    bad = []
    unscannable = []
    for residue in residues:
        if residue.unsupported:
            unscannable.append((residue, residue.unsupported))
            continue
        try:
            edges, degrees, error = inferred_graph(
                residue.atoms, xyz[residue.ids], Chem, rdDetermineBonds)
            extra = edges - residue.expected
            missing = residue.expected - edges
            if extra or missing or error:
                bad.append((residue, extra, missing, error))
        except Exception as exc:
            unscannable.append((residue, f"connectivity inference failed: {exc}"))
    return bad, unscannable


def show_scan(title, count, bad, unscannable):
    say(f"\n=== {title} ===")
    say(f"Residues: {count}; consistent: {count-len(bad)-len(unscannable)}; "
        f"problematic: {len(bad)}; unscannable: {len(unscannable)}")
    for r, extra, missing, err in bad:
        say(f"  {r.label}: extra={len(extra)}, missing={len(missing)}; "
            f"sanitization={'FAIL: ' + err if err else 'PASS (graph mismatch)'}")
    for r, reason in unscannable:
        say(f"  NOT SCANNED: {r.label}: {reason}")


def contact_keys(contacts):
    return {tuple(sorted((a.line_index, b.line_index))) for _, a, b in contacts}


def repair_one(residue, xyz, clean, pdb, topo_bonds, args, mm, app, unit,
               np, Chem, rdDetermineBonds):
    """Return (verified updated coordinates, success) without modifying inputs."""
    r = residue
    say(f"\n{'=' * 65}")
    say(f"REPAIRING {r.label}")
    if not r.originals:
        say("  Residue is absent from the original PDB. Cannot distinguish "
            "added from observed atoms reliably; skipping.")
        return xyz, False

    if any(a.residue != r.atoms[0].residue for a in r.originals.values()):
        say("  NOTE: residue type differs from original (possible PDBFixer "
            "replacement); atom-name provenance may not be exact.")

    absent = [a for a in r.atoms if a.name not in r.originals]
    say(f"  Atoms absent from raw crystal: "
        f"{sum(not a.hydrogen for a in absent)} heavy, "
        f"{sum(a.hydrogen for a in absent)} hydrogen.")
    for a in absent:
        say(f"    {a.name}: {atom_origin(a, r.originals)}")

    report_graph("BEFORE RELAXATION", r.atoms, xyz[r.ids], r.expected,
                 r.originals, Chem, rdDetermineBonds)
    initial_contacts = report_close_contacts(clean, xyz, r.selected, topo_bonds, np)

    # Protect all originally observed heavy atoms in the first trial.
    mobile_added = {i for i, a in r.selected
                    if a.hydrogen or a.name not in r.originals}
    mobile_side = {i for i, a in r.selected
                   if a.hydrogen or a.name not in r.originals or a.name not in BACKBONE}
    stages = []
    if mobile_added:
        stages.append(("added heavy atoms + hydrogens only", mobile_added))
    if not args.no_sidechain_fallback and mobile_side and mobile_side != mobile_added:
        stages.append(("side chain + hydrogens; observed heavy atoms restrained", mobile_side))

    # Each stage starts from the CURRENT whole receptor, including earlier
    # accepted repairs. A rejected trial never modifies these coordinates.
    for title, mobile in stages:
        say(f"\n--- Trial: {title} ---")
        try:
            repaired = minimize(pdb, xyz, clean, r.selected, r.originals, mobile,
                                args, mm, app, unit, np)
            final = report_graph("AFTER RELAXATION", r.atoms, repaired[r.ids],
                                 r.expected, r.originals, Chem, rdDetermineBonds)
            final_contacts = report_close_contacts(clean, repaired, r.selected,
                                                    topo_bonds, np)
        except Exception as exc:
            say(f"  Trial failed: {exc}")
            continue
        new_contacts = contact_keys(final_contacts) - contact_keys(initial_contacts)
        if final[1] or final[2] or final[3] or new_contacts:
            say("  REJECTED: connectivity defect remains or a new very short "
                "nonbonded contact appeared.")
            continue
        # The old coordinates are left untouched by minimize; also explicitly
        # verify that this residue alone has moved before accepting the trial.
        outside = sorted(set(range(len(clean))) - set(r.ids))
        if outside and not np.array_equal(repaired[outside], xyz[outside]):
            say("  REJECTED: coordinates outside the selected residue changed.")
            continue
        say(f"  FIXED: {r.label}; verified connectivity.")
        return repaired, True

    say(f"  UNRESOLVED: {r.label}. Its input coordinates are unchanged.")
    return xyz, False


def repair_sequential(residues, xyz, clean, pdb, topo_bonds, args,
                      mm, app, unit, np, Chem, rdDetermineBonds):
    """Repair every detected residue serially and rescan after each pass."""
    fixed = set()
    current = np.array(xyz, copy=True)
    bad, unscannable = scan_residues(residues, current, Chem, rdDetermineBonds)
    show_scan("INITIAL WHOLE-RECEPTOR SCAN", len(residues), bad, unscannable)
    if args.diagnose_only:
        for r, _, _, _ in bad:
            say(f"\nDetails for {r.label}:")
            report_graph("DIAGNOSIS", r.atoms, current[r.ids], r.expected,
                         r.originals, Chem, rdDetermineBonds)
        return current, fixed, bad, unscannable

    for pass_no in range(1, args.max_passes + 1):
        if not bad:
            break
        say(f"\n=== PASS {pass_no}/{args.max_passes}: {len(bad)} residues ===")
        progress = 0
        for i, (r, _, _, _) in enumerate(bad, 1):
            say(f"\nResidue {i}/{len(bad)} in this pass: {r.label}")
            updated, ok = repair_one(r, current, clean, pdb, topo_bonds, args,
                                     mm, app, unit, np, Chem, rdDetermineBonds)
            if ok:
                current = updated
                fixed.add(r.key)
                progress += 1
        bad, unscannable = scan_residues(residues, current, Chem, rdDetermineBonds)
        show_scan(f"RESCAN AFTER PASS {pass_no}", len(residues), bad, unscannable)
        if not progress:
            say("No successful repair in this pass; stopping rather than retrying indefinitely.")
            break

    return current, fixed, bad, unscannable



class CrossConnectivity:
    """Use Meeko's installed inter-residue bond perception, not RDKit's
    intra-residue cutoff. Compare atom by atom against OpenMM topology."""

    def __init__(self, residues, clean, topo_bonds, Chem, np,
                 detector, covalent_radius):
        self.residues = {r.key: r for r in residues}
        self.clean = clean
        self.Chem = Chem
        self.np = np
        self.detector = detector
        self.radii = covalent_radius
        self.expected = {pair(i, j) for i, j in topo_bonds
                         if clean[i].key != clean[j].key}
        self.molecules = {}
        self.keymap = {}
        self.atom_numbers = {}
        for r in residues:
            key = f"{r.key[0]}:{r.key[1]}{r.key[2]}"
            rw = Chem.RWMol()
            for i, atom in r.selected:
                a = Chem.Atom("H" if atom.element == "D" else atom.element.capitalize())
                if a.GetAtomicNum() not in covalent_radius:
                    raise ValueError(f"Meeko has no covalent radius for {atom.element}.")
                self.atom_numbers[i] = a.GetAtomicNum()
                rw.AddAtom(a)
            mol = rw.GetMol()
            mol.AddConformer(Chem.Conformer(len(r.ids)))
            self.molecules[key] = (mol, r.atoms[0].residue)
            self.keymap[key] = r
        # Source uses 1.2 * sum of covalent radii. A generous search buffer
        # includes every residue potentially linked to a selected group.
        self.buffer = 2.0 * max(covalent_radius.values()) + 4.0

    def scan(self, xyz, active=None):
        """Return inferred extra and missing inter-residue atom pairs.
        For local scans include all nearby residues and all known partners.
        Final scans always call the actual Meeko detector on the whole protein.
        """
        np = self.np
        if active is None:
            include = set(self.residues)
            expected = self.expected
        else:
            active = set(active)
            active_ids = [i for key in active for i in self.residues[key].ids]
            lower = xyz[active_ids].min(axis=0) - self.buffer
            upper = xyz[active_ids].max(axis=0) + self.buffer
            nearby = np.flatnonzero(np.all((xyz >= lower) & (xyz <= upper), axis=1))
            include = active | {self.clean[int(i)].key for i in nearby}
            expected = {e for e in self.expected
                        if any(self.clean[i].key in active for i in e)}
            include |= {self.clean[i].key for e in expected for i in e}

        mols = {}
        for name, r in self.keymap.items():
            if r.key not in include:
                continue
            mol, resname = self.molecules[name]
            conf = mol.GetConformer()
            for local_i, i in enumerate(r.ids):
                conf.SetAtomPosition(local_i, tuple(float(v) for v in xyz[i]))
            mols[name] = (mol, resname)

        inferred = set()
        for (name1, name2), bonds in self.detector(mols).items():
            r1, r2 = self.keymap[name1], self.keymap[name2]
            if active is not None and r1.key not in active and r2.key not in active:
                continue
            for i, j in bonds:
                inferred.add(pair(r1.ids[i], r2.ids[j]))
        return inferred - expected, expected - inferred

    def threshold(self, i, j):
        return 1.2 * (self.radii[self.atom_numbers[i]] + self.radii[self.atom_numbers[j]])

    def show(self, title, issues, xyz):
        extra, missing = issues
        say(f"\n=== {title} ===")
        say(f"Extra inter-residue bonds: {len(extra)}; "
            f"missing topology inter-residue bonds: {len(missing)}")
        for label, edges in (("EXTRA", extra), ("MISSING", missing)):
            for i, j in sorted(edges):
                a, b = self.clean[i], self.clean[j]
                ra, rb = self.residues[a.key], self.residues[b.key]
                say(f"  {label}: {ra.label} {a.name} (serial {a.serial}) -- "
                    f"{rb.label} {b.name} (serial {b.serial}): "
                    f"{dist(xyz[i], xyz[j]):.3f} A "
                    f"(nominal Meeko cutoff {self.threshold(i,j):.3f} A)")
                say(f"    {a.name}: {atom_origin(a, ra.originals)}; "
                    f"{b.name}: {atom_origin(b, rb.originals)}")
                if a.name in ra.originals and b.name in rb.originals:
                    say(f"    Original distance: {dist(ra.originals[a.name].xyz, rb.originals[b.name].xyz):.3f} A")

    def groups(self, extra):
        """Connected components of residues linked by false bonds."""
        neighbors = defaultdict(set)
        for i, j in extra:
            ka, kb = self.clean[i].key, self.clean[j].key
            neighbors[ka].add(kb)
            neighbors[kb].add(ka)
        unseen = set(neighbors)
        while unseen:
            group = {min(unseen)}
            queue = list(group)
            while queue:
                key = queue.pop()
                for other in neighbors[key] - group:
                    group.add(other)
                    queue.append(other)
            unseen -= group
            yield group


def repair_cross_group(group, cross, xyz, clean, pdb, topo_bonds, args,
                       mm, app, unit, np, Chem, rdDetermineBonds):
    """Repair cross-residue clashes, trying each residue separately first.
    Accept only measurable improvement with no new cross/intra-residue
    connectivity defects; all atoms outside the chosen local region are fixed.
    """
    residues = [cross.residues[k] for k in sorted(group)]
    say("\n" + "=" * 65)
    say("REPAIRING CROSS-RESIDUE GROUP: " + ", ".join(r.label for r in residues))
    if any(r.unsupported or not r.originals for r in residues):
        say("  UNRESOLVED: unsupported residue or uncertain atom provenance.")
        return xyz, set(), False

    # A single-residue trial may be sufficient. A joint trial is the fallback.
    stages = []
    for selected_group in [{r.key} for r in residues] + [group]:
        rs = [cross.residues[k] for k in sorted(selected_group)]
        selected = [v for r in rs for v in r.selected]
        maps = {r.key: r.originals for r in rs}
        added = {i for i,a in selected if a.hydrogen or a.name not in maps[a.key]}
        side = {i for i,a in selected
                if a.hydrogen or a.name not in maps[a.key] or a.name not in BACKBONE}
        if added:
            stages.append((selected_group, selected, maps, added, "added atoms + H"))
        if side != added and side and not args.no_sidechain_fallback:
            stages.append((selected_group, selected, maps, side, "restrained side chains + H"))

    # Try every added-only option before permitting observed side-chain motion.
    stages.sort(key=lambda trial: trial[-1] != "added atoms + H")

    # Compare defects for the WHOLE implicated group, including outside
    # neighbors, irrespective of which subset moves in a trial.
    extra0, missing0 = cross.scan(xyz, group)
    if not extra0:
        return xyz, set(), False
    for active, selected, maps, mobile, title in stages:
        movable_pairs = [e for e in extra0 if any(i in mobile for i in e)]
        if not movable_pairs:
            continue
        say(f"\n--- Trial: {title}; " + ", ".join(cross.residues[k].label for k in sorted(active)) + " ---")
        targets = [(i,j,cross.threshold(i,j)+args.cross_margin)
                   for i,j in movable_pairs]
        try:
            trial = minimize(pdb, xyz, clean, selected, {}, mobile, args,
                             mm, app, unit, np, cross_targets=targets, original_maps=maps)
            extra1, missing1 = cross.scan(trial, group)
            # Preserve all previously correct inferred bonds and reject new ones.
            if extra1 - extra0 or missing1 - missing0:
                say("  REJECTED: new inter-residue connectivity defect.")
                continue
            if len(extra1) >= len(extra0):
                say("  REJECTED: no false inter-residue bond was removed.")
                continue
            local_checks = [inferred_graph(r.atoms, trial[r.ids], Chem, rdDetermineBonds)
                            for r in residues]
            if any(edges != r.expected or err for r, (edges, _, err)
                   in zip(residues, local_checks)):
                say("  REJECTED: intra-residue graph mismatch or sanitization failed.")
                continue
            # Check short nonbonded contacts; no new very short contact accepted.
            if contact_keys(report_close_contacts(clean, trial, selected, topo_bonds, np)) - contact_keys(
                    report_close_contacts(clean, xyz, selected, topo_bonds, np)):
                say("  REJECTED: new very short nonbonded contact.")
                continue
            changed = set(np.flatnonzero(np.any(trial != xyz, axis=1)).tolist())
            if not changed <= mobile:
                raise ValueError("Unexpected atom movement outside mobile atoms.")
            say(f"  ACCEPTED: false inter-residue bonds {len(extra0)} -> {len(extra1)}")
            return trial, {clean[i].key for i in changed}, True
        except Exception as exc:
            say(f"  Trial failed: {exc}")
    say("  UNRESOLVED: no safe local repair accepted; prior coordinates preserved.")
    return xyz, set(), False


def repair_cross_residues(cross, xyz, clean, pdb, topo_bonds, args,
                          mm, app, unit, np, Chem, rdDetermineBonds):
    current = np.array(xyz, copy=True)
    fixed = set()
    issues = cross.scan(current)
    cross.show("WHOLE-RECEPTOR INTER-RESIDUE SCAN", issues, current)
    if args.diagnose_only:
        return current, fixed, issues
    for n in range(1, args.max_passes+1):
        if not issues[0]:
            break
        say(f"\n=== INTER-RESIDUE REPAIR PASS {n}/{args.max_passes} ===")
        progress = False
        for group in list(cross.groups(issues[0])):
            current, updated, ok = repair_cross_group(
                group, cross, current, clean, pdb, topo_bonds, args,
                mm, app, unit, np, Chem, rdDetermineBonds)
            fixed |= updated
            progress |= ok
        issues = cross.scan(current)
        cross.show(f"INTER-RESIDUE RESCAN AFTER PASS {n}", issues, current)
        if not progress:
            say("No safe local repair in this pass; stopping.")
            break
    if issues[1]:
        say("NOTE: missing expected between-residue bonds are reported, not "
            "blindly pulled together. Verify chain gaps/disulfides/topology.")
    return current, fixed, issues


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--original", required=True, type=Path, help="Original crystal PDB.")
    p.add_argument("--clean", required=True, type=Path, help="Initial PDBFixer-cleaned, protonated PDB.")
    p.add_argument("--out", type=Path, help="Output PDB; default: <clean stem>_auto_relaxed.pdb.")
    p.add_argument("--diagnose-only", "--scan-only", action="store_true",
                   help="Scan every residue and list problems only.")
    p.add_argument("--max-passes", type=int, default=3,
                   help="Retry unresolved residues after other repairs (default 3).")
    p.add_argument("--allow-partial", action="store_true",
                   help="Save a clearly marked partial PDB if defects remain. Exit nonzero.")
    p.add_argument("--forcefield", default="amber14-all.xml")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--max-iterations", type=int, default=1000)
    p.add_argument("--tolerance", type=float, default=10.0,
                   help="Minimizer RMS force tolerance, kJ/mol/nm.")
    p.add_argument("--restraint", type=float, default=1000.0,
                   help="Restraint on observed side-chain atoms, kJ/mol/nm^2.")
    p.add_argument("--cross-restraint", type=float, default=20000.0,
                   help="One-sided separation restraint on false inter-residue bonds (kJ/mol/nm^2).")
    p.add_argument("--cross-margin", type=float, default=0.15,
                   help="Extra distance in angstrom above Meeko's nominal inter-residue cutoff.")
    p.add_argument("--no-sidechain-fallback", action="store_true",
                   help="Keep all originally observed heavy atoms fixed.")
    args = p.parse_args(argv)
    if min(args.threads, args.max_iterations, args.max_passes) < 1 or min(
            args.tolerance, args.restraint, args.cross_restraint, args.cross_margin) <= 0:
        p.error("Threads, iterations, passes, tolerance and restraint must be positive.")
    dest = args.out or args.clean.with_name(args.clean.stem + "_auto_relaxed.pdb")
    partial_dest = dest.with_name(dest.stem + "_PARTIAL" + dest.suffix)
    if not args.diagnose_only:
        prohibited = {args.original.resolve(), args.clean.resolve()}
        if dest.resolve() in prohibited or partial_dest.resolve() in prohibited:
            p.error("Output must be separate from both inputs.")
        if dest.exists() or (args.allow_partial and partial_dest.exists()):
            p.error("Output already exists. Choose another name with --out.")
    try:
        import numpy as np
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
        import openmm as mm
        from openmm import app, unit
        from meeko.polymer import find_inter_mols_bonds
        from meeko.utils.rdkitutils import covalent_radius
    except ImportError as exc:
        say(f"ERROR: {exc}. Activate your vdock vina environment.")
        return 2

    _, raw = read_pdb(args.original)
    lines, clean = read_pdb(args.clean)
    pdb = app.PDBFile(str(args.clean))
    oa = list(pdb.topology.atoms())
    xyz = np.array([a.xyz for a in clean], dtype=float)
    openmm_xyz = np.asarray(pdb.positions.value_in_unit(unit.angstrom), dtype=float)
    if len(oa) != len(clean) or openmm_xyz.shape != xyz.shape or not np.allclose(
            openmm_xyz, xyz, atol=1e-4, rtol=0):
        raise ValueError("OpenMM atom ordering/count differs from the cleaned PDB.")
    for i, (a, omatom) in enumerate(zip(clean, oa)):
        rid = str(omatom.residue.id).strip()
        ic = str(getattr(omatom.residue, "insertionCode", "")).strip()
        if (str(omatom.residue.chain.id).strip(), rid, ic) != a.key:
            raise ValueError(f"OpenMM residue mapping mismatch at atom index {i}.")

    topo_bonds = {pair(a.index, b.index) for a, b in pdb.topology.bonds()}
    residues = make_residues(clean, raw, topo_bonds)
    say("=== vdock automatic local residue repair v2 (intra + inter) ===")
    say(f"Original: {args.original}")
    say(f"Cleaned:  {args.clean}")
    say(f"Total atoms: {len(clean)}; residues: {len(residues)}")
    say("Atom provenance uses presence/absence of atom names in the raw PDB.")
    say("Single-residue repairs first, local residue groups if necessary. Accepted coordinates accumulate.")
    cross = CrossConnectivity(residues, clean, topo_bonds, Chem, np,
                              find_inter_mols_bonds, covalent_radius)

    current, fixed, bad, unscannable = repair_sequential(
        residues, xyz, clean, pdb, topo_bonds, args, mm, app, unit,
        np, Chem, rdDetermineBonds)

    # Cross-residue clashes may remain even if all local residue graphs passed.
    # Avoid proceeding with cross repairs if invalid intra-residue graphs remain.
    if args.diagnose_only or not (bad or unscannable):
        current, cross_fixed, inter_issues = repair_cross_residues(
            cross, current, clean, pdb, topo_bonds, args, mm, app, unit,
            np, Chem, rdDetermineBonds)
        fixed |= cross_fixed
    else:
        inter_issues = cross.scan(current)
        cross.show("INTER-RESIDUE SCAN (intra-residue repair incomplete)", inter_issues, current)
    if args.diagnose_only:
        say("\nScan complete. No files written.")
        return 1 if bad or unscannable or any(inter_issues) else 0

    bad, unscannable = scan_residues(residues, current, Chem, rdDetermineBonds)
    inter_issues = cross.scan(current)
    cross.show("FINAL MEEKO INTER-RESIDUE VERIFICATION", inter_issues, current)

    # All current coordinates already match the precision of the output PDB.
    # Identify only genuinely changed coordinates, and ensure that every
    # changed atom belongs to a verified, successfully repaired residue.
    changed = set(np.flatnonzero(np.any(current != xyz, axis=1)).tolist())
    if any(clean[i].key not in fixed for i in changed):
        raise ValueError("Unexpected movement outside successfully repaired residues.")
    show_scan("FINAL WHOLE-RECEPTOR VERIFICATION", len(residues), bad, unscannable)
    say(f"Successfully repaired residues: {len(fixed)}")
    for r in residues:
        if r.key in fixed:
            say(f"  {r.label}")
    say(f"Changed atom coordinates: {len(changed)}")

    unresolved = bool(bad or unscannable or any(inter_issues))
    if unresolved and not args.allow_partial:
        say("REPAIR INCOMPLETE: no output PDB written. Inspect the remaining residues.")
        say("Use --allow-partial only if you explicitly want a partial result.")
        return 1
    if unresolved:
        dest = partial_dest
        say("WARNING: writing an INCOMPLETE partial structure; it is not validated.")
    elif not changed:
        say("No connectivity defects found. No repair was needed; no output written.")
        return 0

    out_lines = update_coordinates(lines, clean, current, changed)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("x", encoding="ascii", newline="\n") as fh:
        fh.write("\n".join(out_lines) + "\n")
    say(f"\n{'PARTIAL' if unresolved else 'SUCCESS'}: {dest}")
    say("Neither input PDB was overwritten.")
    say("The scan verifies intra- and inter-residue connectivity and RDKit sanitization. "
        "It does not guarantee global force-field convergence or every Meeko template match.")
    if not unresolved:
        say("Next: run vdock prep-receptor --skip_fix on this repaired PDB.")
    return 1 if unresolved else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError) as exc:
        say(f"ERROR: {exc}")
        sys.exit(2)
