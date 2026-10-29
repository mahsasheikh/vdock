"""
vdock analysis tools: binding-site finding, protein-ligand interactions, visual reports.

  sites         list co-crystallised ligands + detect empty pockets; writes a box file per site
  interactions  H-bonds, hydrophobic, pi-stacking, salt bridges, halogen bonds ... (ProLIF)
  report        interactive HTML report (3Dmol.js, works offline) + PyMOL script
"""
from __future__ import annotations

import csv
import html
import json
import logging
import re
from pathlib import Path

import numpy as np

log = logging.getLogger("vdock")
HERE = Path(__file__).resolve().parent

# Hetero groups that are almost always crystallisation additives / buffer, not ligands.
ADDITIVES = set("""HOH WAT DOD H2O SO4 PO4 GOL EDO PEG PG4 PGE 1PE P6G MPD DMS DMSO ACT ACY ACE FMT
EOH IPA BME DTT TRS MES EPE HEPES CIT TAR MLI IMD BU3 NO3 SCN IOD BR CL NA K LI CS NH4 UNX UNL PE4
2PE 12P 15P NHE BTB CXS FLC PGO PGR TLA SIN AZI CO3 BCT""".split())
IONS = set("ZN MG CA MN FE FE2 CO NI CU CU1 CD HG".split())


# ============================================================================ sites
def _protein_atoms(st):
    """Heavy atoms of polymer residues (protein/nucleic) as (xyz array, list of residue labels)."""
    import gemmi
    xyz, lab = [], []
    for ch in st[0]:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if res.het_flag == "H" and not (info and info.is_amino_acid()):
                continue
            if res.name in ADDITIVES:
                continue
            for at in res:
                if at.element.name in ("H", "D"):
                    continue
                xyz.append((at.pos.x, at.pos.y, at.pos.z))
                lab.append(f"{res.name}{res.seqid.num}{res.seqid.icode.strip()}.{ch.name}")
    return np.array(xyz), lab


def _hetero_groups(st):
    import gemmi
    out = []
    for ch in st[0]:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if res.het_flag != "H" or res.name in ("HOH", "WAT", "DOD"):
                continue
            if info and info.is_amino_acid():  # modified residues (MSE etc.) are part of the chain
                continue
            heavy = [(a.pos.x, a.pos.y, a.pos.z) for a in res if a.element.name not in ("H", "D")]
            if not heavy:
                continue
            kind = "ion/metal" if res.name in IONS else "additive?" if res.name in ADDITIVES else "ligand"
            out.append({"resname": res.name, "chain": ch.name, "resnum": f"{res.seqid.num}{res.seqid.icode.strip()}",
                        "kind": kind, "xyz": np.array(heavy)})
    return out


def _fib_sphere(n=30):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


def detect_pockets(prot_xyz, spacing=1.0, probe=3.0, ray_len=12.0, n_rays=30, min_buried=0.70, min_volume=40):
    """Grid-based buriedness pocket detection (LIGSITE/PocketPicker-style).
    A grid point is a pocket point if it is empty (>= probe Å from protein heavy atoms) yet
    enclosed: >= min_buried of n_rays rays cast from it hit protein within ray_len Å.
    Connected pocket points are clustered; clusters are ranked by volume x mean buriedness."""
    from scipy.spatial import cKDTree
    from scipy import ndimage

    tree = cKDTree(prot_xyz)
    lo = prot_xyz.min(0) - 2
    hi = prot_xyz.max(0) + 2
    shape = np.ceil((hi - lo) / spacing).astype(int) + 1
    gx, gy, gz = (np.arange(s) * spacing + lo[i] for i, s in enumerate(shape))
    G = np.stack(np.meshgrid(gx, gy, gz, indexing="ij"), -1).reshape(-1, 3)
    dist, _ = tree.query(G, distance_upper_bound=probe + 6)
    occupied = (dist < 2.0).reshape(shape)                  # protein "solid" voxels for ray casting
    cand = np.where((dist >= probe) & (dist < probe + 4))[0]  # empty but near the surface
    if cand.size == 0:
        return []
    dirs = _fib_sphere(n_rays)
    steps = np.arange(spacing, ray_len + 1e-6, spacing)
    P = G[cand]
    hits = np.zeros(len(cand), dtype=np.int16)
    for d in dirs:
        hit = np.zeros(len(cand), dtype=bool)
        for s in steps:
            idx = np.round((P + d * s - lo) / spacing).astype(int)
            inb = np.all((idx >= 0) & (idx < shape), 1)
            ok = np.zeros(len(cand), dtype=bool)
            ok[inb] = occupied[idx[inb, 0], idx[inb, 1], idx[inb, 2]]
            hit |= ok
        hits += hit
    bur = hits / n_rays
    keep = cand[bur >= min_buried]
    if keep.size == 0:
        return []
    mask = np.zeros(G.shape[0], dtype=bool)
    mask[keep] = True
    bur_full = np.zeros(G.shape[0])
    bur_full[cand] = bur
    labels, n = ndimage.label(mask.reshape(shape), structure=np.ones((3, 3, 3)))
    labels = labels.reshape(-1)
    pockets = []
    for k in range(1, n + 1):
        ids = np.where(labels == k)[0]
        vol = len(ids) * spacing ** 3
        if vol < min_volume:
            continue
        pts = G[ids]
        pockets.append({"points": pts, "volume": vol, "buriedness": float(bur_full[ids].mean()),
                        "score": vol * float(bur_full[ids].mean())})
    pockets.sort(key=lambda p: -p["score"])
    return pockets


def _lining(prot_xyz, labels, pts, cutoff=4.0):
    from scipy.spatial import cKDTree
    t = cKDTree(pts)
    d, _ = t.query(prot_xyz, distance_upper_bound=cutoff)
    res = sorted({labels[i] for i in np.where(np.isfinite(d))[0]},
                 key=lambda s: (s.split(".")[-1], int(re.sub(r"\D", "", s.split(".")[0][3:]) or 0)))
    return res


def cmd_sites(a):
    import gemmi
    from vdock import write_box

    st = gemmi.read_structure(str(a.structure))
    st.setup_entities()
    if a.chains:
        keep = set(a.chains.split(","))
        for i in reversed(range(len(st[0]))):
            if st[0][i].name not in keep:
                del st[0][i]
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    prot_xyz, labels = _protein_atoms(st)
    rows = []

    def add_site(sid, source, xyz_site, note=""):
        lo, hi = xyz_site.min(0), xyz_site.max(0)
        center = (lo + hi) / 2
        size = np.maximum(hi - lo + 2 * a.padding, a.min_size)
        box = out / f"{sid}_box.txt"
        write_box([float(c) for c in center], [float(s) for s in size], box, box.with_suffix(".pdb"))
        lin = _lining(prot_xyz, labels, xyz_site, a.contact_cutoff)
        rows.append({"site": sid, "source": source, "note": note,
                     "center": " ".join(f"{c:.2f}" for c in center), "size": " ".join(f"{s:.1f}" for s in size),
                     "n_lining_residues": len(lin), "lining_residues": " ".join(lin), "box_file": box.name})

    # 1) co-crystallised hetero groups
    hets = _hetero_groups(st)
    for h in sorted(hets, key=lambda h: -len(h["xyz"])):
        if h["kind"] != "ligand" and not a.all_hetero:
            continue
        add_site(f"lig_{h['resname']}_{h['chain']}{h['resnum']}", "co-crystallised", h["xyz"],
                 f"{h['kind']}, {len(h['xyz'])} heavy atoms")

    # 2) geometric pockets (also useful for apo structures and to find allosteric sites)
    pockets = detect_pockets(prot_xyz, min_buried=a.min_buried)
    poc_pdb = []
    for i, p in enumerate(pockets[: a.max_pockets], 1):
        overl = ""
        for h in hets:
            if h["kind"] == "ligand":
                from scipy.spatial import cKDTree
                if cKDTree(p["points"]).query(h["xyz"], distance_upper_bound=2.0)[0].__lt__(2.0).mean() > 0.3:
                    overl = f"overlaps {h['resname']} {h['chain']}{h['resnum']}"
        add_site(f"pocket{i:02d}", "geometric", p["points"],
                 f"volume {p['volume']:.0f} Å^3, buriedness {p['buriedness']:.2f}" + (f"; {overl}" if overl else ""))
        for x, y, z in p["points"]:
            poc_pdb.append(f"HETATM{len(poc_pdb)+1:5d}  C   POC P{i:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{p['buriedness']:6.2f}           C")
    (out / "pockets.pdb").write_text("\n".join(poc_pdb) + "\nEND\n")

    with open(out / "sites.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()) if rows else ["site"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n{'site':<26}{'source':<17}{'#res':>5}  note")
    for r in rows:
        print(f"{r['site']:<26}{r['source']:<17}{r['n_lining_residues']:>5}  {r['note']}")
    print(f"\nDetails (lining residues, box) -> {out/'sites.csv'};  boxes -> {out}/*_box.txt;  pockets -> {out/'pockets.pdb'}")
    print("Pick a site, check it in the viewer, then use its box file with prep-receptor/dock.")


# ============================================================================ interactions
DEFAULT_INTERACTIONS = ["HBDonor", "HBAcceptor", "Hydrophobic", "PiStacking", "PiCation", "CationPi",
                        "Cationic", "Anionic", "XBDonor", "XBAcceptor", "MetalDonor", "MetalAcceptor"]


def _load_protein(receptor_pdb: Path):
    import prolif as plf
    from rdkit import Chem
    rdk = Chem.MolFromPDBFile(str(receptor_pdb), removeHs=False)
    if rdk is None:
        raise SystemExit(f"Could not read {receptor_pdb} with RDKit. Use the *_clean.pdb written by prep-receptor.")
    if not any(a.GetAtomicNum() == 1 for a in rdk.GetAtoms()):
        raise SystemExit("Receptor has no hydrogens; use receptor/<name>_clean.pdb from prep-receptor.")
    return plf.Molecule.from_rdkit(rdk)


def cmd_interactions(a):
    import prolif as plf
    from rdkit import Chem

    d = Path(a.dockdir)
    sdfs = [Path(a.poses)] if a.poses else sorted((d / "poses").glob("*_out.sdf"))
    if not sdfs:
        raise SystemExit(f"No pose SDF files found in {d/'poses'}")
    if a.ligands:
        wanted = set(a.ligands.split(","))
        sdfs = [s for s in sdfs if s.name[: -len("_out.sdf")] in wanted]
    log.info("Loading receptor %s", a.receptor)
    protein = _load_protein(Path(a.receptor))
    fp = plf.Fingerprint(DEFAULT_INTERACTIONS, count=False)

    detail_rows, matrix = [], {}
    for sdf in sdfs:
        name = sdf.name[: -len("_out.sdf")]
        try:
            poses = list(plf.sdf_supplier(str(sdf)))[: a.n_poses]
        except Exception as e:
            log.warning("%s: cannot read poses (%s)", name, e)
            continue
        for pi, lig in enumerate(poses, 1):
            ifp = fp.generate(lig, protein, metadata=True)
            key = f"{name}" if a.n_poses == 1 else f"{name}#{pi}"
            matrix[key] = {}
            for (lres, pres), inters in ifp.items():
                for itype, metas in inters.items():
                    matrix[key][f"{pres}:{itype}"] = 1
                    for m in metas:
                        li = m["indices"]["ligand"]
                        pj = m["indices"]["protein"]
                        lat = [lig.GetAtomWithIdx(i).GetSymbol() + str(i) for i in li]
                        patoms = []
                        resmol = protein[pres]
                        for j in pj:
                            at = resmol.GetAtomWithIdx(j)
                            info = at.GetPDBResidueInfo()
                            patoms.append(info.GetName().strip() if info else at.GetSymbol())
                        dist = m.get("distance")
                        detail_rows.append({"ligand": name, "pose": pi, "residue": str(pres), "interaction": itype,
                                            "distance_A": round(float(dist), 2) if dist is not None else "",
                                            "ligand_atoms": " ".join(lat), "protein_atoms": " ".join(patoms)})
    out_dir = Path(a.outdir or d)
    det = out_dir / "interactions_detailed.csv"
    with open(det, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["ligand", "pose", "residue", "interaction", "distance_A",
                                           "ligand_atoms", "protein_atoms"])
        w.writeheader()
        w.writerows(detail_rows)
    # ligand x (residue:interaction) 0/1 matrix, sorted by residue number
    cols = sorted({c for v in matrix.values() for c in v},
                  key=lambda c: (int(re.sub(r"\D", "", c.split(".")[0][3:]) or 0), c))
    mat = out_dir / "interactions_matrix.csv"
    with open(mat, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand"] + cols)
        for k, v in matrix.items():
            w.writerow([k] + [v.get(c, 0) for c in cols])
    # frequency of each contact across ligands
    freq = {c: sum(v.get(c, 0) for v in matrix.values()) for c in cols}
    print(f"\nInteractions for {len(matrix)} pose(s). Most frequent contacts:")
    for c, n in sorted(freq.items(), key=lambda x: -x[1])[: a.top]:
        print(f"  {c:<32} {n:>4} / {len(matrix)}")
    if a.key:
        keys = [k.strip() for k in a.key.split(",")]
        print("\nKey-interaction check:")
        for lig_key, v in matrix.items():
            have = [k for k in keys if any(c.startswith(k) for c in v)]
            print(f"  {lig_key:<30} {len(have)}/{len(keys)}  {' '.join(have)}")
    print(f"\n-> {det}\n-> {mat}")


# ============================================================================ report
def _read(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="replace")


def cmd_report(a):
    from vdock import read_box
    d = Path(a.dockdir)
    res_csv = d / "results.csv"
    if not res_csv.exists():
        raise SystemExit(f"{res_csv} not found: run `vdock summarize --dockdir {d}` first")
    rows = list(csv.DictReader(open(res_csv, encoding="utf-8")))[: a.top]
    inter = {}
    det = d / "interactions_detailed.csv"
    if det.exists():
        for r in csv.DictReader(open(det, encoding="utf-8")):
            if str(r["pose"]) == "1":
                inter.setdefault(r["ligand"], []).append(r)
    rec_pdb = _read(Path(a.receptor))
    rec_pdb = "\n".join(l for l in rec_pdb.splitlines() if l.startswith(("ATOM", "HETATM", "TER", "END")))
    box = None
    if a.box:
        c, s = read_box(Path(a.box))
        box = {"center": c, "size": s}
    ligs = []
    for r in rows:
        name = r["ligand"]
        sdf = d / "poses" / f"{name}_out.sdf"
        if not sdf.exists():
            continue
        first = _read(sdf).split("$$$$")[0] + "$$$$\n"
        ints = inter.get(name, [])
        ligs.append({"name": name, "rank": r.get("rank"), "score": r.get("best_score_kcal_mol"),
                     "le": r.get("ligand_efficiency", ""), "rmsd": r.get("rmsd_top_pose_A", ""),
                     "smiles": r.get("smiles", ""), "sdf": first,
                     "inter": [{"res": i["residue"], "type": i["interaction"], "d": i["distance_A"],
                                "lat": i["ligand_atoms"], "pat": i["protein_atoms"]} for i in ints]})
    if not ligs:
        raise SystemExit("No poses to show.")
    js3d = _read(HERE / "assets" / "3Dmol-min.js")
    data = json.dumps({"receptor": rec_pdb, "ligands": ligs, "box": box, "title": a.title})
    out = Path(a.out or d / "report.html")
    out.write_text(_REPORT_HTML.replace("/*3DMOL*/", js3d).replace("/*DATA*/", data.replace("</", "<\\/"))
                   .replace("{{TITLE}}", html.escape(a.title)), encoding="utf-8")
    log.info("HTML report (%d ligands) -> %s   (open it in any browser; works offline)", len(ligs), out)

    # PyMOL script
    pml = out.with_suffix(".pml")
    lines = [f"load {Path(a.receptor).resolve().as_posix()}, receptor",
             "hide everything", "show cartoon, receptor", "set cartoon_transparency, 0.2", "color grey80, receptor"]
    if a.box:
        lines.append(f"load {Path(a.box).with_suffix('.pdb').resolve().as_posix()}, box")
        lines += ["show sticks, box", "color yellow, box", "set stick_radius, 0.05, box"]
    for L in ligs:
        safe = re.sub(r"[^A-Za-z0-9_]", "_", L["name"])
        lines.append(f"load {(d/'poses'/(L['name']+'_out.sdf')).resolve().as_posix()}, {safe}, multiplex=0, state=1")
    first = re.sub(r"[^A-Za-z0-9_]", "_", ligs[0]["name"])
    lines += ["show sticks, not receptor and not box", "util.cbag not receptor and not box",
              f"select site, byres receptor within 5 of {first}", "show sticks, site and not elem H",
              "util.cbac site", f"distance hbonds, {first}, receptor, mode=2", "color magenta, hbonds",
              "label site and name CA, '%s%s' % (resn, resi)", "set label_size, 12",
              "hide sticks, elem H and not (elem H and neighbor elem N+O)",
              f"orient {first}", "disable all", "enable receptor", "enable box", f"enable {first}", "enable hbonds",
              "# Tip: tick other ligand objects in the right panel to compare poses; `set state, 2` shows the next pose."]
    pml.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log.info("PyMOL script -> %s   (File > Run Script..., or: pymol %s)", pml, pml.name)


_REPORT_HTML = r"""<!doctype html><html><head><meta charset="utf-8"><title>{{TITLE}}</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>
:root{--bg:#fafafa;--fg:#1d1d1f;--mut:#6b6b70;--line:#e3e3e8;--acc:#2f6fdb;--sel:#e8f0fe}
*{box-sizing:border-box}body{margin:0;font:14px/1.4 system-ui,Segoe UI,Roboto,sans-serif;background:var(--bg);color:var(--fg)}
header{padding:12px 16px;border-bottom:1px solid var(--line);display:flex;gap:16px;align-items:baseline;flex-wrap:wrap}
h1{font-size:17px;margin:0}header span{color:var(--mut)}
main{display:grid;grid-template-columns:minmax(300px,420px) 1fr;height:calc(100vh - 50px)}
#list{overflow:auto;border-right:1px solid var(--line)}table{border-collapse:collapse;width:100%}
th,td{padding:6px 8px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}
th{position:sticky;top:0;background:#f0f0f3;font-weight:600;cursor:pointer}td.n{text-align:right;font-variant-numeric:tabular-nums}
tr.row{cursor:pointer}tr.row:hover{background:#f3f3f6}tr.sel{background:var(--sel)}
#right{display:grid;grid-template-rows:1fr auto;min-width:0}#viewer{position:relative;width:100%;height:100%}
#panel{max-height:38vh;overflow:auto;border-top:1px solid var(--line);padding:8px 12px;background:#fff}
#panel h2{font-size:15px;margin:4px 0}.smi{font-family:ui-monospace,Consolas,monospace;font-size:12px;color:var(--mut);word-break:break-all}
.ctl{display:flex;gap:12px;flex-wrap:wrap;margin:6px 0}.ctl label{display:flex;gap:4px;align-items:center}
.tag{display:inline-block;padding:1px 6px;border-radius:8px;font-size:12px;background:#eee}
.HBDonor,.HBAcceptor{background:#dbe9ff}.Hydrophobic{background:#eaf5d9}.PiStacking,.PiCation,.CationPi{background:#f3e3ff}
.Cationic,.Anionic{background:#ffe1dc}.XBDonor,.XBAcceptor{background:#e0f7f4}
@media (max-width:800px){main{grid-template-columns:1fr;grid-template-rows:40vh 60vh;height:auto}#list{border-right:0}}
</style></head><body>
<header><h1>{{TITLE}}</h1><span id="sub"></span></header>
<main><div id="list"><table><thead><tr><th>#</th><th>Ligand</th><th>Score</th><th>LE</th><th>Contacts</th></tr></thead><tbody id="tb"></tbody></table></div>
<div id="right"><div id="viewer"></div><div id="panel"></div></div></main>
<script>/*3DMOL*/</script>
<script>
const D=/*DATA*/;
const viewer=$3Dmol.createViewer(document.getElementById('viewer'),{backgroundColor:'white'});
const rec=viewer.addModel(D.receptor,'pdb');
rec.setStyle({}, {cartoon:{color:'#c8c8d0',opacity:0.85}});
if(D.box){const c=D.box.center,s=D.box.size;viewer.addBox({center:{x:c[0],y:c[1],z:c[2]},dimensions:{w:s[0],h:s[1],d:s[2]},wireframe:true,color:'#d4a017'});}
let lig=null, shapes=[], labels=[], cur=0;
const opts={site:true,surf:false,hbond:true,labels:true,boxv:true};
document.getElementById('sub').textContent=D.ligands.length+' ligands · click a row to view its top pose';
const tb=document.getElementById('tb');
D.ligands.forEach((L,i)=>{const tr=document.createElement('tr');tr.className='row';tr.onclick=()=>show(i);
 tr.innerHTML=`<td class=n>${L.rank}</td><td>${L.name}</td><td class=n>${(+L.score).toFixed(2)}</td><td class=n>${L.le||''}</td><td class=n>${L.inter.length}</td>`;tb.appendChild(tr);});
function show(i){cur=i;const L=D.ligands[i];[...tb.children].forEach((r,j)=>r.classList.toggle('sel',j===i));
 if(lig)viewer.removeModel(lig);shapes.forEach(s=>viewer.removeShape(s));shapes=[];labels.forEach(l=>viewer.removeLabel(l));labels=[];viewer.removeAllSurfaces();
 lig=viewer.addModel(L.sdf,'sdf');
 const la=lig.selectedAtoms({});
 const polarH=new Set();la.forEach(a=>{if(a.elem!=='H')return;(a.bonds||[]).forEach(b=>{if(['N','O','S'].includes(la[b].elem))polarH.add(a.index);});});
 lig.setStyle({}, {stick:{colorscheme:'greenCarbon',radius:0.22}});
 lig.setStyle({predicate:a=>a.elem==='H'&&!polarH.has(a.index)},{});
 // residues with any atom within 5 A of the ligand
 const heavy=la.filter(a=>a.elem!=='H');const keys=new Set();
 rec.selectedAtoms({}).forEach(a=>{for(const b of heavy){const dx=a.x-b.x,dy=a.y-b.y,dz=a.z-b.z;if(dx*dx+dy*dy+dz*dz<25){keys.add(a.chain+':'+a.resi);break;}}});
 const nearSel={predicate:a=>keys.has(a.chain+':'+a.resi)};
 const cart={cartoon:{color:'#c8c8d0',opacity:0.55}};
 rec.setStyle({},cart);
 if(opts.site){rec.setStyle({predicate:a=>keys.has(a.chain+':'+a.resi)&&a.elem!=='H'},
   {cartoon:{color:'#c8c8d0',opacity:0.55},stick:{colorscheme:'whiteCarbon',radius:0.14}});}
 if(opts.surf)viewer.addSurface($3Dmol.SurfaceType.VDW,{opacity:0.5,color:'white'},{model:rec,predicate:nearSel.predicate});
 const ra=L.inter.map(x=>x.res);const seen=new Set();
 if(opts.labels){L.inter.forEach(x=>{if(seen.has(x.res))return;seen.add(x.res);const m=x.res.match(/^([A-Z]{3})(-?\d+)\.(\S+)$/);if(!m)return;
  const a=rec.selectedAtoms({resn:m[1],resi:+m[2],chain:m[3],atom:'CA'})[0];if(a)labels.push(viewer.addLabel(m[1]+m[2],{position:a,fontSize:11,backgroundOpacity:0.6}));});}
 if(opts.hbond){L.inter.filter(x=>/HB|XB|Cationic|Anionic|Metal/.test(x.type)).forEach(x=>{
  const m=x.res.match(/^([A-Z]{3})(-?\d+)\.(\S+)$/);if(!m)return;const li=+x.lat.split(' ')[0].replace(/\D/g,'');const pn=x.pat.split(' ')[0];
  const pa=rec.selectedAtoms({resn:m[1],resi:+m[2],chain:m[3],atom:pn})[0];const lA=la[li];
  if(pa&&lA)shapes.push(viewer.addCylinder({start:pa,end:lA,radius:0.06,dashed:true,fromCap:1,toCap:1,color:/HB/.test(x.type)?'#1f6feb':'#d63384'}));});}
 viewer.zoomTo({model:lig});viewer.zoom(0.9);viewer.render();panel(L);}
function panel(L){const groups={};L.inter.forEach(x=>(groups[x.res]=groups[x.res]||[]).push(x));
 let h=`<h2>#${L.rank} ${L.name} &nbsp; ${(+L.score).toFixed(2)} kcal/mol</h2><div class=smi>${L.smiles||''}</div>`;
 if(L.rmsd)h+=`<div>RMSD to reference: ${L.rmsd} Å</div>`;
 h+=`<div class=ctl>${['site:binding-site residues','surf:pocket surface','hbond:H-bond / ionic lines','labels:residue labels'].map(t=>{const[k,v]=t.split(':');return `<label><input type=checkbox data-k=${k} ${opts[k]?'checked':''}>${v}</label>`}).join('')}</div>`;
 if(!L.inter.length)h+='<p>No interaction data: run <code>vdock interactions</code> before <code>vdock report</code>.</p>';
 else h+='<table><tr><th>Residue</th><th>Interactions (distance Å)</th></tr>'+Object.entries(groups).map(([r,xs])=>`<tr><td>${r}</td><td>${xs.map(x=>`<span class="tag ${x.type}">${x.type}${x.d?' '+x.d:''}</span>`).join(' ')}</td></tr>`).join('')+'</table>';
 const p=document.getElementById('panel');p.innerHTML=h;p.querySelectorAll('input').forEach(c=>c.onchange=()=>{opts[c.dataset.k]=c.checked;show(cur);});}
show(0);window.addEventListener('resize',()=>viewer.resize());
</script></body></html>"""
