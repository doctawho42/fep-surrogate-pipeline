"""NIOCH cage screen — a calibrated, per-enantiomer target-hypothesis CAMPAIGN.

This is the operational debt: run the cage through docking and hand back ranked target
hypotheses + recommended assays. It is a REPORT, not a paper claim ("pending assays").

Honest by construction (the gate, Fig H, showed why): raw docking scores are not
comparable across pockets (they pick the greasiest pocket). So we DO NOT rank by raw
kcal/mol. We rank each target by the cage's **percentile among that target's known
binders** docked into the same pocket -- "the cage docks better than X% of known
T-inhibitors". That normalises the pocket-size artefact. Per enantiomer (invariant #6;
racemate ~ stronger binder). Caveats (promiscuous ureide motif, limited panel, docking
is approximate) are reported, not hidden.
"""
from __future__ import annotations

import csv
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from screen.dock import dock, prepare_target  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent.parent
DATA = ROOT / "data"
CAMP = DATA / "campaign"
N_REF = 15  # reference known-binders docked per target

# 11-target diverse panel (the gate panel + PDE10A, flagged by SwissTargetPrediction)
PANEL = [
    ("EGFR (kinase)", "CHEMBL203", "2ITY"),
    ("CDK2 (kinase)", "CHEMBL301", "1H00"),
    ("HIV-1 protease", "CHEMBL243", "1HSG"),
    ("Acetylcholinesterase", "CHEMBL220", "4EY7"),
    ("Estrogen receptor α", "CHEMBL206", "3ERT"),
    ("Carbonic anhydrase II", "CHEMBL205", "3HS4"),
    ("Factor Xa", "CHEMBL244", "2P16"),
    ("PDE5A", "CHEMBL1827", "1UDT"),
    ("Thrombin", "CHEMBL204", "1OYT"),
    ("Glucocorticoid receptor", "CHEMBL2034", "3BQD"),
    ("PDE10A", "CHEMBL4409", "3HR1"),
]
CAGE = {"R,R (given)": str(DATA / "cage" / "cage_given.sdf"),
        "S,S (mirror)": str(DATA / "cage" / "cage_enantiomer.sdf")}


def _load_cache():
    f = CAMP / "scores.csv"
    c = {}
    if f.exists():
        for lig, pdb, sc in csv.reader(open(f)):
            c[(lig, pdb)] = float(sc)
    return c


def _embed(smiles, path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3()
    p.randomSeed = 1
    if AllChem.EmbedMolecule(m, p) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(m)
    Chem.SDWriter(str(path)).write(m)
    return str(path)


def run(progress=print, n_workers=6):
    from concurrent.futures import ThreadPoolExecutor

    CAMP.mkdir(parents=True, exist_ok=True)
    (CAMP / "lig").mkdir(exist_ok=True)
    targets = {pdb: prepare_target(pdb, str(DATA / "pockets")) for _, _, pdb in PANEL}
    cache = _load_cache()

    # assemble jobs: cage (both enantiomers) x all pockets + N_REF known binders / own pocket
    jobs = []  # (lig_id, lig_sdf, pdb)
    cage_jobs = {}
    for ename, sdf in CAGE.items():
        for _, _, pdb in PANEL:
            cage_jobs[(ename, pdb)] = f"cage::{ename}"
            jobs.append((f"cage::{ename}", sdf, pdb))
    ref_ids = {pdb: [] for _, _, pdb in PANEL}
    for _tname, cid, pdb in PANEL:
        n = 0
        for mid, smi in list(csv.reader(open(DATA / "figH" / f"actives_{cid}.csv")))[:N_REF * 2]:
            if n >= N_REF:
                break
            sdfp = CAMP / "lig" / f"{mid}.sdf"
            s = str(sdfp) if sdfp.exists() else _embed(smi, sdfp)
            if s:
                jobs.append((mid, s, pdb))
                ref_ids[pdb].append(mid)
                n += 1

    todo = [j for j in jobs if (j[0], j[2]) not in cache]
    progress(f"docking {len(todo)} new pairs ({len(jobs)-len(todo)} cached) on {n_workers} workers")

    def run_one(j):
        lig_id, sdf, pdb = j
        return lig_id, pdb, dock(targets[pdb]["receptor"], sdf, targets[pdb]["ref_ligand"])

    cf = open(CAMP / "scores.csv", "a", newline="")
    w = csv.writer(cf)
    with ThreadPoolExecutor(max_workers=n_workers) as ex:
        for lig_id, pdb, sc in ex.map(run_one, todo):
            cache[(lig_id, pdb)] = sc
            w.writerow([lig_id, pdb, sc])
            cf.flush()
    cf.close()

    cage_scores = {k: cache[(cage_jobs[k], k[1])] for k in cage_jobs}
    ref = {pdb: [cache[(m, pdb)] for m in ids if cache.get((m, pdb), float("inf")) != float("inf")]
           for pdb, ids in ref_ids.items()}

    # calibrate: percentile of the cage among known binders (lower dock score = better)
    rows = []
    for tname, _, pdb in PANEL:
        dist = sorted(ref[pdb])
        out = {"target": tname, "pdb": pdb, "n_ref": len(dist)}
        for ename in CAGE:
            cs = cage_scores[(ename, pdb)]
            pct = 100.0 * sum(1 for r in dist if r > cs) / len(dist) if dist else float("nan")
            out[ename] = cs
            out[f"pct::{ename}"] = pct
        out["best_pct"] = max(out[f"pct::{e}"] for e in CAGE)
        rows.append(out)
    rows.sort(key=lambda r: -r["best_pct"])
    return rows, cage_scores, ref


if __name__ == "__main__":
    rows, _, _ = run()
    print("\nRanked target hypotheses (cage percentile among known binders):")
    for r in rows:
        es = "  ".join(f"{e} {r[e]:+.1f}({r['pct::'+e]:.0f}%)" for e in CAGE)
        print(f"  {r['target']:24s} best {r['best_pct']:5.0f}%   {es}")
