"""Broad CALIBRATED reverse docking — an honest hypothesis-class shortlist for the cage.

A ~33-target diverse-fold pocketome. The cage (both enantiomers, acetate + deacetylated
"active" form) is docked into every pocket and ranked NOT by raw kcal/mol (greasy-pocket
artefact, Fig H) but by its **percentile among that target's known binders** docked into
the same pocket (the calibration that killed the DHODH false-positive). Output: a ranked,
heavily-caveated shortlist of mechanistic hypotheses for experimental testing -- NOT a
target call.

Run:  python figs/make_reversedock.py   (docking cached to data/reversedock/, resumable;
assembles the report/figure from whatever is cached).
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from screen.cage import enantiomer_smiles  # noqa: E402
from screen.dock import dock, prepare_target  # noqa: E402

DATA = ROOT / "data" / "reversedock"
N_REF = 10
N_WORKERS = 6
AC = "CC(=O)O[C@]12C[C@H](C3=CC=CC=C3C1(F)F)C4=C(NC(=O)NC4=O)O2"
DE = "O[C@]12C[C@H](C3=CC=CC=C3C1(F)F)C4=C(NC(=O)NC4=O)O2"
FORMS = {"acetate R,R": AC, "acetate S,S": enantiomer_smiles(AC),
         "deacetyl R,R": DE, "deacetyl S,S": enantiomer_smiles(DE)}

# (name, chembl_target, pdb, explicit ligand resname or None)
PANEL = [
    ("EGFR", "CHEMBL203", "2ITY", None), ("CDK2", "CHEMBL301", "1H00", None),
    ("CDK6", "CHEMBL2508", "2EUF", None), ("ABL1", "CHEMBL1862", "2HYY", None),
    ("VEGFR2", "CHEMBL279", "4ASD", None), ("JAK2", "CHEMBL2971", "3KRR", None),
    ("GSK3b", "CHEMBL262", "1Q5K", None), ("p38a MAPK", "CHEMBL260", "3GCS", None),
    ("HIV protease", "CHEMBL243", "1HSG", None), ("Thrombin", "CHEMBL204", "1OYT", None),
    ("Factor Xa", "CHEMBL244", "2P16", None), ("BACE1", "CHEMBL4822", "2WJO", None),
    ("Cathepsin K", "CHEMBL268", "1MEM", None), ("DPP-4", "CHEMBL284", "2P8S", None),
    ("Renin", "CHEMBL286", "2V0Z", None), ("HIV integrase", "CHEMBL3471", "3NF6", None),
    ("Estrogen receptor a", "CHEMBL206", "3ERT", None), ("Androgen receptor", "CHEMBL1871", "2AM9", None),
    ("Glucocorticoid receptor", "CHEMBL2034", "3BQD", None), ("PPARg", "CHEMBL235", "2PRG", None),
    ("b2-adrenergic", "CHEMBL210", "2RH1", "CAU"), ("A2A adenosine", "CHEMBL251", "3EML", "ZMA"),
    ("D3 dopamine", "CHEMBL234", "3PBL", "ETQ"),
    ("Acetylcholinesterase", "CHEMBL220", "4EY7", None), ("MAO-B", "CHEMBL2039", "2V60", "C17"),
    ("Carbonic anhydrase II", "CHEMBL205", "3HS4", None), ("COX-2", "CHEMBL230", "3LN1", "CEL"),
    ("PDE5A", "CHEMBL1827", "1UDT", None), ("PDE10A", "CHEMBL4409", "3HR1", None),
    ("DHODH", "CHEMBL1966", "1D3G", "BRE"), ("DHFR", "CHEMBL202", "1U72", "MTX"),
    ("Thymidylate synthase", "CHEMBL1952", "1HVY", None), ("HSP90a", "CHEMBL3880", "1UYG", None),
    ("BRD4 bromodomain", "CHEMBL1163125", "3MXF", None),
]


def _style():
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                         "savefig.dpi": 300, "savefig.bbox": "tight"})


def fetch_actives(cid):
    f = DATA / f"actives_{cid}.csv"
    if f.exists():
        return [tuple(r) for r in csv.reader(open(f))]
    figh = ROOT / "data" / "figH" / f"actives_{cid}.csv"
    if figh.exists():
        return [tuple(r) for r in csv.reader(open(figh))]
    rows, off = [], 0
    while len(rows) < 200:
        u = (f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={cid}"
             f"&pchembl_value__gte=6.5&limit=1000&offset={off}")
        d = json.load(urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": "rd/1.0"}), timeout=60))
        for a in d["activities"]:
            if a.get("canonical_smiles") and a.get("molecule_chembl_id"):
                rows.append((a["molecule_chembl_id"], a["canonical_smiles"]))
        if not d["page_meta"].get("next"):
            break
        off += 1000
    DATA.mkdir(parents=True, exist_ok=True)
    csv.writer(open(f, "w", newline="")).writerows(rows)
    return rows


def embed(lig_id, smi):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    sdf = DATA / "lig" / f"{lig_id.replace('::', '_').replace(' ', '_').replace(',', '')}.sdf"
    if sdf.exists():
        return str(sdf)
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    m = Chem.AddHs(m)
    p = AllChem.ETKDGv3()
    p.randomSeed = 1
    if AllChem.EmbedMolecule(m, p) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(m)
    sdf.parent.mkdir(parents=True, exist_ok=True)
    Chem.SDWriter(str(sdf)).write(m)
    return str(sdf)


def run_docking(progress=print):
    DATA.mkdir(parents=True, exist_ok=True)
    targets, ok_panel = {}, []
    for name, cid, pdb, lig in PANEL:
        try:
            targets[pdb] = prepare_target(pdb, str(ROOT / "data" / "pockets"), lig_resname=lig)
            ok_panel.append((name, cid, pdb, lig))
        except Exception as e:
            progress(f"  drop {name} ({pdb}): {type(e).__name__}")
    cache = {}
    cf = DATA / "scores.csv"
    if cf.exists():
        for lig_id, pdb, sc in csv.reader(open(cf)):
            cache[(lig_id, pdb)] = float(sc)

    jobs, ref_ids = [], {pdb: [] for _, _, pdb, _ in ok_panel}
    for name, smi in FORMS.items():
        s = embed(f"cage::{name}", smi)
        for _, _, pdb, _ in ok_panel:
            if s and (f"cage::{name}", pdb) not in cache:
                jobs.append((f"cage::{name}", s, pdb))
    for name, cid, pdb, _ in ok_panel:
        n = 0
        for mid, smi in fetch_actives(cid)[:N_REF * 3]:
            if n >= N_REF:
                break
            s = embed(mid, smi)
            if not s:
                continue
            ref_ids[pdb].append(mid)
            n += 1
            if (mid, pdb) not in cache:
                jobs.append((mid, s, pdb))
    progress(f"{len(ok_panel)} targets; docking {len(jobs)} new pairs on {N_WORKERS} workers")

    def one(j):
        lig_id, sdf, pdb = j
        return lig_id, pdb, dock(targets[pdb]["receptor"], sdf, targets[pdb]["ref_ligand"])

    out = open(cf, "a", newline="")
    w = csv.writer(out)
    done = 0
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex:
        for lig_id, pdb, sc in ex.map(one, jobs):
            cache[(lig_id, pdb)] = sc
            w.writerow([lig_id, pdb, sc])
            out.flush()
            done += 1
            if done % 40 == 0:
                progress(f"  {done}/{len(jobs)}")
    out.close()
    return ok_panel, ref_ids, cache


def assemble():
    cache = {}
    for lig_id, pdb, sc in csv.reader(open(DATA / "scores.csv")):
        cache[(lig_id, pdb)] = float(sc)
    ref_ids = {pdb: [] for _, _, pdb, _ in PANEL}
    for _, cid, pdb, _ in PANEL:
        for mid, _smi in fetch_actives(cid)[:N_REF * 3]:
            if len(ref_ids[pdb]) >= N_REF:
                break
            ref_ids[pdb].append(mid)
    rows = []
    for name, _, pdb, _ in PANEL:
        dist = sorted(cache[(m, pdb)] for m in ref_ids[pdb]
                      if (m, pdb) in cache and np.isfinite(cache[(m, pdb)]))
        rec = {"t": name, "pdb": pdb, "n": len(dist)}
        for e in FORMS:
            cs = cache.get((f"cage::{e}", pdb), float("nan"))
            rec[e] = cs
            rec["p_" + e] = (100 * sum(1 for r in dist if r > cs) / len(dist)
                             if dist and np.isfinite(cs) else float("nan"))
        rec["best"] = max((rec["p_" + e] for e in FORMS if rec["p_" + e] == rec["p_" + e]), default=float("nan"))
        rows.append(rec)
    rows.sort(key=lambda r: -(r["best"] if r["best"] == r["best"] else -1))
    return rows


def write_outputs(rows):
    _style()
    have = [r for r in rows if r["n"] >= 5]
    fig, ax = plt.subplots(figsize=(6.8, max(3.5, 0.32 * len(have))))
    y = np.arange(len(have))[::-1]
    ax.barh(y, [r["best"] for r in have], color=["#D55E00" if r["best"] >= 90 else "#0072B2" for r in have])
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['t']} (n={r['n']})" for r in have], fontsize=7)
    ax.axvline(50, color="#999", ls=":", lw=1)
    ax.set_xlabel("best cage percentile among known binders (calibrated docking)")
    ax.set_title("Calibrated reverse docking — cage hypothesis shortlist (PRELIMINARY; orange=likely greasy-pocket artefact)",
                 fontsize=8, fontweight="bold")
    fig.tight_layout()
    for e in ("pdf", "png"):
        fig.savefig(str(ROOT / "figs" / ("reversedock_shortlist." + e)))
    lines = ["# Cage — broad calibrated reverse-docking shortlist (PRELIMINARY · hypotheses, not a target)",
             "",
             "Auto-generated by `make reversedock`. The cage is docked into a ~33-target",
             "diverse pocketome and ranked by its **percentile among each target's known",
             "binders** (calibrated; raw kcal/mol is meaningless across pockets — Fig H).",
             "**This is a hypothesis generator, not an answer** (the gate showed structure-based",
             "target-ID is unreliable; the DHODH test showed even a top hypothesis can be bottom-4%).",
             "",
             "| target | n_ref | best pctl | note |", "|---|---:|---:|---|"]
    for r in have:
        note = "**greasy-pocket artefact likely**" if r["best"] >= 90 else ("worth a focused assay" if r["best"] >= 70 else "")
        lines.append(f"| {r['t']} | {r['n']} | {r['best']:.0f}% | {note} |")
    lines += ["", "Caveats: docking approximate; large hydrophobic pockets score the rigid",
              "lipophilic cage high regardless (trust data-rich, moderate rows over thin 100%s);",
              "promiscuous ureide; panel finite. Per-enantiomer + deacetyl(active) forms tested.",
              "Confirm any hypothesis experimentally (panel / Cell Painting / chemoproteomics)."]
    (ROOT / "docs" / "reversedock_shortlist.md").write_text("\n".join(lines) + "\n")


def main():
    run_docking()
    rows = assemble()
    write_outputs(rows)
    print("wrote docs/reversedock_shortlist.md + figs/reversedock_shortlist.{pdf,png}")
    for r in rows[:12]:
        print(f"  {r['t']:24s} n={r['n']:>2} best={r['best'] if r['best']==r['best'] else 0:.0f}%")


if __name__ == "__main__":
    main()
