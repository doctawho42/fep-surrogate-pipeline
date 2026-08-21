"""Fig H — THE GATE: retrospective structure-based vs ligand-based target-ID.

Decides one paper vs two (docs/target_finding_plan.md §2). A diverse-fold panel of
targets (holo PDB pocket + ChEMBL actives). Each test ligand (scaffold-disjoint from its
target's training actives -- the orphan scenario) is docked against ALL panel pockets;
the true target should rank top. Compared to the ligand-shape baseline that failed in
Fig F (max Tanimoto to a target's training actives) and to random.

GATE: does structure (docking) recovery beat the ligand-shape baseline? Beats -> the
structure arm stands -> one paper (cage motivates, public data validates). Does not ->
cage to Paper 2, ship the methods spine as-is. Do not rescue the merge if it doesn't stand.

Run:  python figs/make_figH.py  (or `make figH`). Docking is parallelised and cached to
data/figH/ ; ChEMBL pulls cached too. Re-runs are cheap.
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
import urllib.request
import warnings
from concurrent.futures import ThreadPoolExecutor

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from screen.dock import dock, prepare_target  # noqa: E402
from screen.recovery import recovery_at_k, recovery_auroc  # noqa: E402

FIGDIR = ROOT / "figs"
DATA = ROOT / "data" / "figH"
SEED = 7
N_TEST = 5          # scaffold-disjoint test ligands per target
N_WORKERS = 6

# diverse-fold panel: (name, chembl_target_id, holo_pdb_id)
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
]
C_STRUCT, C_SHAPE, C_RAND, C_REF = "#0072B2", "#D55E00", "#999999", "#555555"


def _style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


def fetch_actives(chembl_id, cache, max_records=400):
    if cache.exists():
        return [tuple(r) for r in csv.reader(open(cache))]
    rows, offset = [], 0
    while len(rows) < max_records:
        url = (f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={chembl_id}"
               f"&pchembl_value__gte=6.5&limit=1000&offset={offset}")
        d = json.load(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "figH/1.0"}), timeout=60))
        for a in d["activities"]:
            if a.get("canonical_smiles") and a.get("molecule_chembl_id"):
                rows.append((a["molecule_chembl_id"], a["canonical_smiles"]))
        if not d["page_meta"].get("next"):
            break
        offset += 1000
    cache.parent.mkdir(parents=True, exist_ok=True)
    csv.writer(open(cache, "w", newline="")).writerows(rows)
    return rows


def build_queries():
    """Per target: dedup, drop cross-panel promiscuous, scaffold-split, take N_TEST
    scaffold-disjoint test ligands. Returns test queries + per-target train fingerprints."""
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator as fpg
    from rdkit.Chem.Scaffolds import MurckoScaffold

    raw = {}  # mol_id -> (smiles, set(target_idx))
    for ti, (_, cid, _) in enumerate(PANEL):
        for mid, smi in fetch_actives(cid, DATA / f"actives_{cid}.csv"):
            raw.setdefault(mid, [smi, set()])[1].add(ti)
    gen = fpg.GetMorganGenerator(radius=2, fpSize=2048)
    rng = np.random.default_rng(SEED)
    train_fps = {ti: [] for ti in range(len(PANEL))}
    by_target = {ti: [] for ti in range(len(PANEL))}
    for mid, (smi, tset) in raw.items():
        if len(tset) != 1:
            continue  # drop promiscuous
        ti = next(iter(tset))
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=m) or smi
        by_target[ti].append((mid, smi, scaf, gen.GetFingerprint(m)))

    queries = []  # (mol_id, smiles, true_target_idx)
    for ti, items in by_target.items():
        scaffs = sorted({s for _, _, s, _ in items})
        rng.shuffle(scaffs)
        n_test_scaf = max(1, int(0.3 * len(scaffs)))
        test_scaf = set(scaffs[:n_test_scaf])
        test_items, n_added = [], 0
        for mid, smi, scaf, fp in items:
            if scaf in test_scaf and n_added < N_TEST:
                test_items.append((mid, smi)); n_added += 1
            elif scaf not in test_scaf:
                train_fps[ti].append(fp)
        for mid, smi in test_items:
            queries.append((mid, smi, ti))
    return queries, train_fps


def embed_sdf(smiles, path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = 1
    if AllChem.EmbedMolecule(m, p) != 0:
        return None
    AllChem.MMFFOptimizeMolecule(m)
    Chem.SDWriter(str(path)).write(m)
    return str(path)


def main() -> None:
    _style()
    DATA.mkdir(parents=True, exist_ok=True)
    print("preparing pockets...")
    targets = [prepare_target(pdb, str(ROOT / "data" / "pockets")) for _, _, pdb in PANEL]
    print("building scaffold-disjoint queries...")
    queries, train_fps = build_queries()
    print(f"  {len(queries)} test ligands across {len(PANEL)} targets")

    # dock cache
    cache_f = DATA / "dock_scores.csv"
    cache = {}
    if cache_f.exists():
        for mid, pdb, sc in csv.reader(open(cache_f)):
            cache[(mid, pdb)] = float(sc)

    # embed test ligands once
    sdfs = {}
    for mid, smi, _ in queries:
        p = DATA / "lig" / f"{mid}.sdf"
        p.parent.mkdir(parents=True, exist_ok=True)
        sdfs[mid] = embed_sdf(smi, p) if not p.exists() else str(p)

    jobs = [(mid, ti2) for (mid, _, _) in queries for ti2 in range(len(PANEL))
            if (mid, PANEL[ti2][2]) not in cache and sdfs[mid]]
    print(f"docking {len(jobs)} (ligand x pocket) pairs on {N_WORKERS} workers...")

    def run_one(job):
        mid, ti2 = job
        tg = targets[ti2]
        return mid, PANEL[ti2][2], dock(tg["receptor"], sdfs[mid], tg["ref_ligand"])

    done = 0
    with ThreadPoolExecutor(max_workers=N_WORKERS) as ex, open(cache_f, "a", newline="") as cf:
        w = csv.writer(cf)
        for mid, pdb, sc in ex.map(run_one, jobs):
            cache[(mid, pdb)] = sc; w.writerow([mid, pdb, sc]); cf.flush()
            done += 1
            if done % 50 == 0:
                print(f"  {done}/{len(jobs)}")

    # score matrices
    from rdkit.Chem import DataStructs
    from rdkit.Chem import rdFingerprintGenerator as fpg
    from rdkit import Chem
    gen = fpg.GetMorganGenerator(radius=2, fpSize=2048)
    n_t = len(PANEL)
    struct = np.full((len(queries), n_t), np.inf)
    shape = np.zeros((len(queries), n_t))
    true_idx = np.array([ti for _, _, ti in queries])
    for q, (mid, smi, _) in enumerate(queries):
        qfp = gen.GetFingerprint(Chem.MolFromSmiles(smi))
        for ti2 in range(n_t):
            struct[q, ti2] = cache.get((mid, PANEL[ti2][2]), np.inf)
            tf = train_fps[ti2]
            shape[q, ti2] = max(DataStructs.BulkTanimotoSimilarity(qfp, tf), default=0.0) if tf else 0.0

    valid = np.isfinite(struct).all(1)
    struct, shape, true_idx = struct[valid], shape[valid], true_idx[valid]
    rec_s = recovery_at_k(struct, true_idx, lower_better=True)
    rec_l = recovery_at_k(shape, true_idx, lower_better=False)
    rec_r = np.array([k / n_t for k in range(1, n_t + 1)])
    auroc_s = recovery_auroc(struct, true_idx, lower_better=True)
    auroc_l = recovery_auroc(shape, true_idx, lower_better=False)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.6, 3.3))
    kk = np.arange(1, n_t + 1)
    axA.plot(kk, rec_s, "-o", color=C_STRUCT, ms=4, lw=1.7, label=f"structure/docking (AUROC {auroc_s:.2f})")
    axA.plot(kk, rec_l, "-s", color=C_SHAPE, ms=4, lw=1.7, label=f"ligand-shape (Fig F) (AUROC {auroc_l:.2f})")
    axA.plot(kk, rec_r, "--", color=C_RAND, lw=1.3, label="random")
    axA.set_xlabel("k (targets retrieved)"); axA.set_ylabel("top-k recovery")
    axA.set_title("A   hidden-target recovery (scaffold-disjoint)", loc="left", fontweight="bold")
    axA.set_xticks(kk); axA.legend(frameon=False, loc="lower right"); axA.set_ylim(0, 1.02)

    axB.bar([0, 1, 2], [rec_s[0], rec_l[0], rec_r[0]], color=[C_STRUCT, C_SHAPE, C_RAND], width=0.6)
    for x, v in zip([0, 1, 2], [rec_s[0], rec_l[0], rec_r[0]]):
        axB.text(x, v + 0.01, f"{v:.2f}", ha="center", fontsize=8)
    axB.set_xticks([0, 1, 2]); axB.set_xticklabels(["structure", "ligand-shape", "random"], rotation=12)
    axB.set_ylabel("top-1 recovery"); axB.set_ylim(0, 1.0)
    axB.set_title("B   THE GATE: top-1", loc="left", fontweight="bold")

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figH_structure_target_id.{ext}")
    print(f"\nwrote figH_structure_target_id.(pdf|png); n_test={valid.sum()}")
    print(f"  top-1 recovery: structure {rec_s[0]:.2f}  ligand-shape {rec_l[0]:.2f}  random {rec_r[0]:.2f}")
    print(f"  top-3 recovery: structure {rec_s[2]:.2f}  ligand-shape {rec_l[2]:.2f}")
    print(f"  recovery AUROC: structure {auroc_s:.3f}  ligand-shape {auroc_l:.3f}")
    gate = "PASS (structure > ligand-shape)" if rec_s[0] > rec_l[0] + 0.05 or auroc_s > auroc_l + 0.03 \
        else "FAIL (structure <= ligand-shape) -> cage to Paper 2"
    print(f"  GATE: {gate}")


if __name__ == "__main__":
    main()
