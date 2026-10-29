#!/usr/bin/env python3
"""Read-only structural and result checks for vdock v1.1 (Python 3.11, RDKit, Gemmi).

Examples (from a project directory; CMD):
  python C:\\docking\\vina-docking\\vdock_safety_checks.py ligand --sdf BXZ_ideal.sdf --structure input\\3BM9.cif --resname BXZ --chain A --resseq 1
  python C:\\docking\\vina-docking\\vdock_safety_checks.py hetero --structure input\\3BM9.cif --prepared receptor\\receptor_clean.pdb --resname BXZ --chain A --resseq 1
  python C:\\docking\\vina-docking\\vdock_safety_checks.py states --sdf ligands\\scrubbed.sdf --scores docking\\results_all_states.csv
  python C:\\docking\\vina-docking\\vdock_safety_checks.py rmsd --scores validation_seed42\\results.csv
"""
import argparse
import csv
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


def require_file(p):
    p = Path(p)
    if not p.is_file():
        raise SystemExit(f"FILE NOT FOUND: {p.resolve()}\nCurrent directory: {Path.cwd()}")
    return p


def parse_mols(path):
    from rdkit import Chem
    p = require_file(path)
    sup = Chem.SDMolSupplier(str(p), removeHs=False)
    records = []
    for i, mol in enumerate(sup, 1):
        if mol is None:
            print(f"WARNING: SDF record {i} could not be parsed")
        else:
            records.append(mol)
    if not records:
        raise SystemExit(f"No valid molecules in {p}")
    return records


def read_structure(path):
    import gemmi
    st = gemmi.read_structure(str(require_file(path)))
    if not st or len(st) == 0:
        raise SystemExit("Structure has no models")
    return st


def locate_ligand(st, resname, chain, resseq):
    matches = []
    for ch in st[0]:
        for res in ch:
            if res.name.upper() == resname.upper() and (not chain or chain == ch.name) and (resseq is None or res.seqid.num == resseq):
                matches.append((ch.name, res))
    if len(matches) != 1:
        msg = ", ".join(f"{c}:{r.name}{r.seqid.num}" for c,r in matches) or "none"
        raise SystemExit(f"Need exactly one reference residue {resname} (chain={chain}, resseq={resseq}); found {len(matches)}: {msg}")
    return matches[0]


def heavy_positions(res):
    import numpy as np
    xyz = [[at.pos.x, at.pos.y, at.pos.z] for at in res if at.element.name not in ('H', 'D')]
    return np.asarray(xyz, dtype=float)


def box_center(x):
    return (x.min(axis=0) + x.max(axis=0)) / 2.0


def cmd_ligand(args):
    import numpy as np
    mols = parse_mols(args.sdf)
    if len(mols) != 1:
        print(f"WARNING: {len(mols)} SDF records; 'vdock box' would combine all their coordinates.")
    m = mols[0]
    if not m.GetNumConformers():
        raise SystemExit("No 3D conformer in SDF")
    p = np.asarray(m.GetConformer().GetPositions(), dtype=float)
    heavy_indices = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() > 1]
    h = p[heavy_indices]
    center = box_center(h)
    print(f"SDF: {args.sdf} | {m.GetNumHeavyAtoms()} heavy / {m.GetNumAtoms()} total atoms")
    print(f"SDF heavy bounding-box center: {np.round(center,3)} A")
    print(f"SDF geometrical extent: {np.round(h.max(0)-h.min(0),3)} A")
    print(f"SDF has 3D conformer flag: {m.GetConformer().Is3D()}")
    if args.structure:
        st = read_structure(args.structure)
        chain, res = locate_ligand(st, args.resname, args.chain, args.resseq)
        rh = heavy_positions(res)
        rc = box_center(rh)
        diff = float(np.linalg.norm(center - rc))
        print(f"Crystal {res.name} {chain}:{res.seqid.num} heavy atoms: {len(rh)}")
        print(f"Crystal heavy bounding-box center: {np.round(rc,3)} A")
        print(f"Center separation: {diff:.2f} A")
        if len(rh) != len(h):
            print("WARNING: Heavy-atom count differs. Investigate SDF chemical identity.")
        if diff > 3.0:
            print("FAIL: SDF is not in the crystal reference frame. DO NOT use to build a crystal-ligand docking box or calculate in-place RMSD.")
        else:
            print("CENTER CHECK OK: centers overlap approximately; verify atom matching and pose, not only centroid.")


def cmd_hetero(args):
    import numpy as np
    st = read_structure(args.structure)
    lig_ch, lig_res = locate_ligand(st, args.resname, args.chain, args.resseq)
    lig_xyz = heavy_positions(lig_res)
    water_names = {'HOH', 'WAT', 'DOD', 'H2O'}
    metals = {'ZN', 'MG', 'CA', 'MN', 'FE', 'FE2', 'CO', 'NI', 'CU', 'CU1', 'CD', 'HG'}
    original = []
    for ch in st[0]:
        for res in ch:
            if res is lig_res:
                continue
            if not (res.het_flag == 'H' or res.name in water_names):
                continue
            x = heavy_positions(res)
            if not len(x):
                continue
            minimum = float(np.sqrt(((x[:,None,:]-lig_xyz[None,:,:])**2).sum(axis=2)).min())
            kind = 'water' if res.name in water_names else 'metal/ion' if res.name in metals else 'other heterogen'
            original.append((ch.name,res.name,res.seqid.num,res.seqid.icode.strip(),minimum,kind))
    counts = Counter(r[-1] for r in original)
    print(f"Original structure: {args.structure} | reference ligand {args.resname} {lig_ch}:{lig_res.seqid.num}")
    print(f"Heterogens: " + ", ".join(f"{k}={n}" for k,n in sorted(counts.items())))
    near = [r for r in original if r[4] <= args.cutoff]
    print(f"Within {args.cutoff:g} A of crystal ligand: {len(near)} heterogen groups")
    print(" KIND             CHAIN RESNAME RESSEQ     closest_A    cleanPDB  VinaPDBQT")
    prep_res = None
    pdbqt_res = None
    if args.pdbqt:
        pdbqt = require_file(args.pdbqt)
        pdbqt_res = set()
        for line in pdbqt.read_text(errors='replace').splitlines():
            if line.startswith(('ATOM  ', 'HETATM')):
                try:
                    pdbqt_res.add((line[21:22].strip(), line[17:20].strip(), int(line[22:26]), line[26:27].strip()))
                except ValueError:
                    continue
    if args.prepared:
        prep = read_structure(args.prepared)
        prep_res = {(ch.name, r.name, r.seqid.num, r.seqid.icode.strip()) for ch in prep[0] for r in ch}
    for ch,nam,num,ic,d,kind in sorted(near, key=lambda z:z[4]):
        kept = '-' if prep_res is None else ('YES' if (ch,nam,num,ic) in prep_res else 'NO')
        qkept = '-' if pdbqt_res is None else ('YES' if (ch,nam,num,ic) in pdbqt_res else 'NO')
        print(f" {kind:<16} {ch:<5} {nam:<7} {num}{ic:<3} {d:10.2f}    {kept:<8} {qkept}")
    for ch,nam,num,ic,d,kind in sorted((r for r in original if r[-1]=='metal/ion' and r[4] > args.cutoff), key=lambda z:z[4]):
        kept = '-' if prep_res is None else ('YES' if (ch,nam,num,ic) in prep_res else 'NO')
        qkept = '-' if pdbqt_res is None else ('YES' if (ch,nam,num,ic) in pdbqt_res else 'NO')
        print(f" {kind:<16} {ch:<5} {nam:<7} {num}{ic:<3} {d:10.2f}    {kept:<8} {qkept} (distant)")
    if args.prepared or args.pdbqt:
        print("Presence checks compare original residue labels; renumbering can cause false absences. Presence does NOT verify coordination, charges, or water orientation.")
    print("NOTE: The pipeline's default PDBFixer call removes all heterogens except waters if --keep_water is explicitly supplied;" )
    print("      retaining every water is NOT automatically scientifically correct.")


def make_state_names(mols):
    def safe(x): return re.sub(r'[^A-Za-z0-9_.-]+','_',x).strip('_') or 'lig'
    counts = Counter(safe(m.GetProp('_Name')) if m.HasProp('_Name') and m.GetProp('_Name') else f'lig{i:06d}' for i,m in enumerate(mols))
    seen = Counter()
    for i,m in enumerate(mols):
        base = safe(m.GetProp('_Name')) if m.HasProp('_Name') and m.GetProp('_Name') else f'lig{i:06d}'
        seen[base] += 1
        name = f'{base}_s{seen[base]}' if counts[base] > 1 else base
        yield base, name, m


def read_scores(path):
    p = require_file(path)
    with open(p, newline='', encoding='utf-8-sig') as fh:
        rows = list(csv.DictReader(fh))
    if not rows or 'ligand' not in rows[0] or 'best_score_kcal_mol' not in rows[0]:
        raise SystemExit("Scores CSV must contain ligand and best_score_kcal_mol. Re-run 'summarize' without --collapse_states.")
    return rows


def parent_name(row):
    return row.get('parent') or re.sub(r'_s\d+$','',row['ligand'])


def num(v):
    try: return float(v)
    except (ValueError, TypeError): return None


def cmd_states(args):
    from rdkit import Chem
    mols = parse_mols(args.sdf)
    groups = defaultdict(list)
    for base, name, m in make_state_names(mols):
        smi = Chem.MolToSmiles(Chem.RemoveHs(m), isomericSmiles=True)
        charge = sum(a.GetFormalCharge() for a in m.GetAtoms())
        groups[base].append((name,smi,charge,m.GetNumHeavyAtoms()))
    print(f"Total SDF records: {len(mols)} | parent names: {len(groups)}")
    print("NOTE: These are enumerated structures, not measured pH-state populations. A pKa tool is needed for fractions.")
    score_by_name = {r['ligand']:r for r in read_scores(args.scores)} if args.scores else {}
    for parent, records in sorted(groups.items()):
        print(f'\nParent: {parent} | {len(records)} records | {len(set(x[1] for x in records))} chemically distinct states')
        for name,smi,charge,heavy in records:
            r = score_by_name.get(name, {})
            sc = r.get('best_score_kcal_mol','-')
            print(f'  {name:<30} charge={charge:+d} heavy={heavy:<3d} vina={sc:>8}  {smi}')
        scored = [(n,num(score_by_name[n].get('best_score_kcal_mol'))) for n,*_ in records if n in score_by_name]
        scored = [(n,s) for n,s in scored if s is not None]
        if scored and len(records)>1:
            bestn,bests = min(scored,key=lambda x:x[1])
            print(f'  RAW SCORE BEST: {bestn} ({bests:+.2f}) — no microstate population weighting was applied.')
    missing = set(score_by_name) - {name for recs in groups.values() for name,*_ in recs}
    if missing: print(f"\nWARNING: {len(missing)} result IDs not matched to SDF-derived names (examples: {sorted(missing)[:5]}).")


def cmd_rmsd(args):
    rows = read_scores(args.scores)
    groups = defaultdict(list)
    for row in rows:
        score = num(row.get('best_score_kcal_mol'))
        rmsd = num(row.get('rmsd_top_pose_A'))
        if score is None: continue
        groups[parent_name(row)].append((row['ligand'],score,rmsd))
    any_rmsd=False
    for parent, recs in sorted(groups.items()):
        all_with_r = [r for r in recs if r[2] is not None]
        if not all_with_r: continue
        any_rmsd=True
        global_best = min(recs,key=lambda r:r[1])
        min_r = min(all_with_r,key=lambda r:r[2])
        print(f"Parent {parent} | {len(recs)} scored states")
        print(f"  BEST-SCORING state: {global_best[0]} score={global_best[1]:+.3f} RMSD={global_best[2]} A")
        print(f"  LOWEST-RMSD state:  {min_r[0]} score={min_r[1]:+.3f} RMSD={min_r[2]:.3f} A")
        if min_r[2] <= 2.0 and (global_best[2] is None or global_best[2] > 2.0):
            print("  WARNING: Minimum-RMSD PASS would be misleading: the global best-scoring state did NOT pass.")
        if global_best[2] is not None:
            print(f"  Top-score success <= 2 A: {global_best[2] <= 2.0}")
    if not any_rmsd: print("No RMSD values found. Use a crystal-instance reference SDF with 'vdock summarize --reference'.")
    print("NOTE: Scores/RMSDs concern pose 1 of each prepared state, not all nine reported poses.")


def main():
    p = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command',required=True)
    t=sub.add_parser('ligand',help='Check SDF coordinate alignment with crystal ligand')
    t.add_argument('--sdf',required=True);t.add_argument('--structure');t.add_argument('--resname',default='BXZ');t.add_argument('--chain',default='A');t.add_argument('--resseq',type=int,default=1);t.set_defaults(func=cmd_ligand)
    t=sub.add_parser('hetero',help='List nearby waters, ions, heterogens and loss upon preparation')
    t.add_argument('--structure',required=True);t.add_argument('--prepared');t.add_argument('--pdbqt');t.add_argument('--resname',default='BXZ');t.add_argument('--chain',default='A');t.add_argument('--resseq',type=int,default=1);t.add_argument('--cutoff',type=float,default=5.0);t.set_defaults(func=cmd_hetero)
    t=sub.add_parser('states',help='Enumerate ligand microstates and associated docking scores')
    t.add_argument('--sdf',required=True);t.add_argument('--scores');t.set_defaults(func=cmd_states)
    t=sub.add_parser('rmsd',help='Compare RMSD of best-scoring state with minimum across states')
    t.add_argument('--scores',required=True);t.set_defaults(func=cmd_rmsd)
    args=p.parse_args();args.func(args)

if __name__=='__main__': main()
