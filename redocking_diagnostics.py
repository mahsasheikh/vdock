#!/usr/bin/env python
"""Read-only vdock 1.1 redocking pose diagnostics.

Run from a vdock project directory:
    python redocking_diagnostics.py --reference input/reference_crystal.sdf \
        --dockdir validation_BXZ_corrected/docking --ligand BXZ \
        --prepared validation_BXZ_corrected/ligand/scrubbed.sdf

Measures in-place heavy-atom RMSD; it deliberately DOES NOT realign structures.
"""
from __future__ import annotations
import argparse
import csv
import re
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem


def load_mols(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Missing: {path.resolve()}")
    if path.suffix.lower() in ('.sdf', '.mol'):
        records = list(Chem.SDMolSupplier(str(path), removeHs=False))
    else:
        raise ValueError('SDF/MOL only')
    mols = [m for m in records if m is not None]
    if not mols:
        raise ValueError(f"No parseable molecules in {path}")
    if len(records) != len(mols):
        print(f"WARNING: {len(records) - len(mols)} unparsable record(s) in {path}")
    for m in mols:
        if not m.GetNumConformers():
            raise ValueError(f"A molecule has no coordinates: {path}")
    return mols


def heavy_coords(mol):
    conf = mol.GetConformer()
    idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    return np.array([list(conf.GetAtomPosition(i)) for i in idx], dtype=float)


def skeleton(m):
    rw = Chem.RWMol(Chem.RemoveHs(m, sanitize=False))
    for at in rw.GetAtoms():
        at.SetFormalCharge(0)
        at.SetIsAromatic(False)
        at.SetNoImplicit(True)
        at.SetNumExplicitHs(0)
    for bond in rw.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    return rw.GetMol()


def best_rmsd(pose, ref):
    # Follow vdock's symmetry-aware skeleton mapping; do NOT superpose.
    ps, rs = skeleton(pose), skeleton(ref)
    if ps.GetNumAtoms() != rs.GetNumAtoms():
        raise ValueError(f"Heavy atom counts disagree (pose {ps.GetNumAtoms()} vs reference {rs.GetNumAtoms()})")
    px, rx = ps.GetConformer().GetPositions(), rs.GetConformer().GetPositions()
    matches = ps.GetSubstructMatches(rs, uniquify=False, useChirality=False, maxMatches=10000)
    if not matches:
        raise ValueError("Ligand connectivity/atom mapping differs; RMSD undefined")
    val = min(float(np.sqrt(np.mean(np.sum((px[list(mapping)] - rx) ** 2, axis=1)))) for mapping in matches)
    return val, len(matches)


def descriptor(m):
    heavy = Chem.RemoveHs(m)
    charge = sum(a.GetFormalCharge() for a in m.GetAtoms())
    return f"{Chem.MolToSmiles(heavy)} (net charge {charge:+d}, heavy atoms {heavy.GetNumAtoms()})"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--reference', required=True, type=Path, help='Crystal instance SDF, not ideal CCD ligand')
    ap.add_argument('--dockdir', type=Path, help='vdock output directory containing poses/')
    ap.add_argument('--pose-sdf', type=Path, help='Alternative: analyze a directly exported pose SDF (e.g. local-only test)')
    ap.add_argument('--ligand', default='BXZ', help='Prepared ligand basename; default BXZ')
    ap.add_argument('--prepared', type=Path, help='Optional molscrub-generated ligand SDF to compare molecular states')
    ap.add_argument('--out', type=Path, help='Optional diagnostic CSV path')
    args = ap.parse_args()

    ref = load_mols(args.reference)[0]
    if args.pose_sdf:
        sdf = args.pose_sdf
        pdbqt = args.pose_sdf.with_suffix('.pdbqt')
    else:
        if args.dockdir is None:
            ap.error('Supply either --dockdir or --pose-sdf')
        sdf = args.dockdir / 'poses' / f'{args.ligand}_out.sdf'
        pdbqt = args.dockdir / 'poses' / f'{args.ligand}_out.pdbqt'
    docked = load_mols(sdf)
    if args.prepared:
        prepared = load_mols(args.prepared)
        print(f'Reference: {descriptor(ref)}')
        print(f'Prepared ligand states ({len(prepared)} records):')
        for i, m in enumerate(prepared, 1):
            print(f'  {i}: {descriptor(m)}')
    else:
        print(f'Reference: {descriptor(ref)}')

    scores = []
    if pdbqt.exists():
        scores = [float(x) for x in re.findall(r'REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)', pdbqt.read_text(errors='replace'))]
        if len(scores) != len(docked):
            print(f'WARNING: parsed {len(scores)} Vina scores, but {len(docked)} exported SDF poses; score-to-pose order may not match.')
    else:
        print('WARNING: corresponding PDBQT is missing; scores unavailable')

    ref_ctr = (heavy_coords(ref).min(axis=0) + heavy_coords(ref).max(axis=0))/2
    print('Reference bounding-box center:', ', '.join(f'{x:.3f}' for x in ref_ctr), 'Å')
    rows = []
    for i, m in enumerate(docked, 1):
        rmsd, nmatch = best_rmsd(m, ref)
        centroid_shift = float(np.linalg.norm(heavy_coords(m).mean(axis=0) - heavy_coords(ref).mean(axis=0)))
        rows.append({'pose': i, 'score_kcal_mol': scores[i-1] if i <= len(scores) and len(scores) == len(docked) else '',
                     'rmsd_A': round(rmsd, 3), 'centroid_distance_A': round(centroid_shift, 3), 'symmetry_mappings': nmatch})
    print('\nDocked poses (RMSD to crystal, no structural alignment):')
    print(f'{"Pose":>5} {"Score":>11} {"RMSD Å":>10} {"Centroid Å":>12} {"Maps":>6}')
    for r in rows:
        score = f'{r["score_kcal_mol"]:.3f}' if r['score_kcal_mol'] != '' else 'n/a'
        print(f'{r["pose"]:>5} {score:>11} {r["rmsd_A"]:>10.3f} {r["centroid_distance_A"]:>12.3f} {r["symmetry_mappings"]:>6}')
    best_geom = min(rows, key=lambda r: r['rmsd_A'])
    print(f'\nTop-ranked pose RMSD: {rows[0]["rmsd_A"]:.3f} Å')
    print(f'Minimum RMSD among saved poses: {best_geom["rmsd_A"]:.3f} Å (pose {best_geom["pose"]})')
    if rows[0]['rmsd_A'] <= 2.0:
        print('Interpretation: top-ranked pose recovers crystal geometry by the usual 2 Å heuristic.')
    elif best_geom['rmsd_A'] <= 2.0:
        print('Interpretation: a native-like pose exists among saved results, but it was not ranked first (scoring/ranking issue).')
    else:
        print('Interpretation: none of the saved poses is near-native; test sampling, receptor/ligand preparation, and score of native pose.')
    outfile = args.out or ((args.dockdir or sdf.parent) / f'{args.ligand}_pose_diagnostics.csv')
    outfile.parent.mkdir(parents=True, exist_ok=True)
    with outfile.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f'Saved: {outfile}')


if __name__ == '__main__':
    main()
