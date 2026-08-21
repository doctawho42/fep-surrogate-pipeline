"""Fig F — retrospective target identification (reverse screening). STRETCH / honest.

Hide the target, recover it from a library. Real ChEMBL ligand->target bioactivity;
Morgan fingerprints; SCAFFOLD-DISJOINT split (Bemis-Murcko, per plan §2). A multiclass
model predicts P(target | ligand); reverse screening ranks targets and asks whether the
true target is recovered top-k, plus the calibration of P(binds). Baseline = the SHAPE
similarity (max Tanimoto to each target's training actives). Falsifier: recovery <= shape.

Two regimes, because they tell different stories:
  EASY  8 DIVERSE single-protein targets (kinase / protease / GPCR / NR / ...): chemically
        very different -> similarity already solves it.
  HARD  8 within-family AMINERGIC GPCRs (5-HT/dopamine/histamine/muscarinic/adrenergic):
        shared basic-amine chemotypes -> similarity is confounded; target-specific SAR
        must do the work. This is the meaningful test.

Honest scope: ligand-based (not structure/pocket-based as the full plan envisions),
labels are bioactivity (not FEP), promiscuous ligands removed. A proof-of-concept of the
reverse-screening EVALUATION, not the structure-aware system.

Run:  python figs/make_figF.py   (or `make figF`). Caches ChEMBL pulls to data/ for
reproducibility (offline after first run).
"""
from __future__ import annotations

import csv
import json
import pathlib
import sys
import urllib.request
import warnings

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "src"))

FIGDIR = pathlib.Path(__file__).resolve().parent
DATA = FIGDIR.parent / "data"
SEED = 7

TARGETS_EASY = {
    "CHEMBL203": "EGFR (kinase)", "CHEMBL204": "Thrombin (protease)",
    "CHEMBL205": "Carbonic anhydrase II", "CHEMBL206": "Estrogen receptor α",
    "CHEMBL217": "Dopamine D2 (GPCR)", "CHEMBL220": "Acetylcholinesterase",
    "CHEMBL218": "Cannabinoid CB1 (GPCR)", "CHEMBL244": "Coagulation factor X",
}
TARGETS_HARD = {  # within-family aminergic GPCRs — confounded similarity
    "CHEMBL214": "5-HT1A", "CHEMBL224": "5-HT2A", "CHEMBL3155": "5-HT7",
    "CHEMBL217": "Dopamine D2", "CHEMBL234": "Dopamine D3", "CHEMBL231": "Histamine H1",
    "CHEMBL216": "Muscarinic M1", "CHEMBL1867": "α-2A adrenergic",
}
C_MODEL, C_SHAPE, C_RAND, C_REF = "#0072B2", "#D55E00", "#999999", "#555555"


def _style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 9, "axes.labelsize": 10, "axes.titlesize": 10,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 7.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


def fetch_target(tid, max_records=800):
    rows, offset = [], 0
    while len(rows) < max_records:
        url = (f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={tid}"
               f"&pchembl_value__gte=6.5&limit=1000&offset={offset}")
        req = urllib.request.Request(url, headers={"User-Agent": "figF/1.0"})
        d = json.load(urllib.request.urlopen(req, timeout=60))
        for a in d["activities"]:
            if a.get("canonical_smiles") and a.get("molecule_chembl_id"):
                rows.append((a["molecule_chembl_id"], a["canonical_smiles"], tid))
        if not d["page_meta"].get("next"):
            break
        offset += 1000
    return rows


def load_data(targets, cache):
    if cache.exists():
        with open(cache) as f:
            return [tuple(r) for r in csv.reader(f)][1:]
    DATA.mkdir(exist_ok=True)
    rows = []
    for tid in targets:
        r = fetch_target(tid)
        print(f"  fetched {len(r):4d} for {targets[tid]}")
        rows.extend(r)
    with open(cache, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["mol_id", "smiles", "target"]); w.writerows(rows)
    return rows


def build_dataset(targets, cache):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    rows = load_data(targets, cache)
    mol_targets: dict[str, set] = {}
    smiles: dict[str, str] = {}
    for mid, smi, tid in rows:
        mol_targets.setdefault(mid, set()).add(tid)
        smiles[mid] = smi
    keep = {m: next(iter(t)) for m, t in mol_targets.items() if len(t) == 1}  # drop promiscuous
    mols, y, scaffs = [], [], []
    for mid, tid in keep.items():
        mol = Chem.MolFromSmiles(smiles[mid])
        if mol is None:
            continue
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol) or smiles[mid]
        except Exception:
            scaf = smiles[mid]
        mols.append(mol); y.append(tid); scaffs.append(scaf)
    return mols, np.array(y), np.array(scaffs)


def scaffold_split(y, scaffs, test_frac=0.3, seed=SEED):
    rng = np.random.default_rng(seed)
    is_test = np.zeros(len(y), dtype=bool)
    for tid in np.unique(y):
        idx = np.where(y == tid)[0]
        groups: dict[str, list] = {}
        for i in idx:
            groups.setdefault(scaffs[i], []).append(i)
        gkeys = list(groups); rng.shuffle(gkeys)
        n_test, acc = int(test_frac * len(idx)), 0
        for k in gkeys:
            if acc >= n_test:
                break
            for i in groups[k]:
                is_test[i] = True
            acc += len(groups[k])
    return ~is_test, is_test


def run_benchmark(targets, cache, min_per=60):
    from rdkit.Chem import DataStructs
    from rdkit.Chem import rdFingerprintGenerator as fpg
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score

    mols, y, scaffs = build_dataset(targets, cache)
    gen = fpg.GetMorganGenerator(radius=2, fpSize=2048)
    X = np.array([gen.GetFingerprintAsNumPy(m) for m in mols])
    bits = [gen.GetFingerprint(m) for m in mols]

    labs, cnts = np.unique(y, return_counts=True)
    labs = [t for t, c in zip(labs, cnts) if c >= min_per]
    mask = np.isin(y, labs)
    X, y, scaffs = X[mask], y[mask], scaffs[mask]
    bits = [b for b, m in zip(bits, mask) if m]

    tr, te = scaffold_split(y, scaffs)
    clf = LogisticRegression(C=2.0, max_iter=2000).fit(X[tr], y[tr])
    classes = list(clf.classes_)
    proba = clf.predict_proba(X[te])
    y_te = y[te]

    def curve(scores):
        ranks = np.array([int(np.where(np.array(classes)[np.argsort(-scores[i])] == t)[0][0])
                          for i, t in enumerate(y_te)])
        return np.array([(ranks < k).mean() for k in range(1, len(classes) + 1)])

    def auroc(scores):
        a = [roc_auc_score((np.array(classes) == t).astype(int), scores[i])
             for i, t in enumerate(y_te) if 0 < (np.array(classes) == t).sum() < len(classes)]
        return float(np.mean(a))

    train_bits = {t: [bits[i] for i in np.where(tr)[0] if y[i] == t] for t in labs}
    te_idx = np.where(te)[0]
    shape = np.array([[max(DataStructs.BulkTanimotoSimilarity(bits[i], train_bits[t]), default=0.0)
                       for t in classes] for i in te_idx])

    pmax, chosen = proba.max(1), np.array(classes)[proba.argmax(1)]
    correct = (chosen == y_te).astype(float)
    return dict(
        classes=classes, n_test=int(te.sum()),
        model_curve=curve(proba), shape_curve=curve(shape),
        rand_curve=np.array([k / len(classes) for k in range(1, len(classes) + 1)]),
        model_auroc=auroc(proba), shape_auroc=auroc(shape),
        pmax=pmax, correct=correct,
    )


def ece_reliability(pmax, correct, nb=8):
    bins = np.linspace(0, 1, nb + 1)
    bi = np.clip(np.digitize(pmax, bins) - 1, 0, nb - 1)
    conf, acc, ece = [], [], 0.0
    for b in range(nb):
        m = bi == b
        if m.sum() >= 5:
            conf.append(pmax[m].mean()); acc.append(correct[m].mean())
            ece += m.mean() * abs(pmax[m].mean() - correct[m].mean())
    return np.array(conf), np.array(acc), ece


def main() -> None:
    _style()
    print("[EASY: diverse families]")
    easy = run_benchmark(TARGETS_EASY, DATA / "figF_chembl_easy.csv")
    print("[HARD: within-family aminergic GPCRs]")
    hard = run_benchmark(TARGETS_HARD, DATA / "figF_chembl_hard.csv")

    pmax = np.concatenate([easy["pmax"], hard["pmax"]])
    correct = np.concatenate([easy["correct"], hard["correct"]])
    conf, acc, ece = ece_reliability(pmax, correct)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.6, 3.3))

    # Panel A: top-1 recovery, easy vs hard, model vs shape
    groups = ["EASY\n(diverse)", "HARD\n(within-family)"]
    model_t1 = [easy["model_curve"][0], hard["model_curve"][0]]
    shape_t1 = [easy["shape_curve"][0], hard["shape_curve"][0]]
    rand_t1 = [easy["rand_curve"][0], hard["rand_curve"][0]]
    xg = np.arange(2); w = 0.36
    axA.bar(xg - w / 2, model_t1, w, color=C_MODEL, label="model")
    axA.bar(xg + w / 2, shape_t1, w, color=C_SHAPE, label="shape baseline")
    for i, r in enumerate(rand_t1):
        axA.plot([xg[i] - 0.5, xg[i] + 0.5], [r, r], color=C_RAND, ls="--", lw=1.2)
    for i in range(2):
        axA.text(xg[i] - w / 2, model_t1[i] + 0.01, f"{model_t1[i]:.2f}", ha="center", fontsize=7.5)
        axA.text(xg[i] + w / 2, shape_t1[i] + 0.01, f"{shape_t1[i]:.2f}", ha="center", fontsize=7.5)
    axA.set_xticks(xg); axA.set_xticklabels(groups)
    axA.set_ylabel("top-1 recovery"); axA.set_ylim(0, 1.08)
    axA.set_title("A   fingerprint model ties shape (both regimes)", loc="left", fontweight="bold")
    axA.legend(frameon=False, loc="lower left")
    axA.text(0.98, 0.04, "dashed = random", transform=axA.transAxes, fontsize=7, ha="right", color=C_RAND)

    # Panel B: P(binds) calibration (pooled)
    axB.plot([0, 1], [0, 1], "--", color=C_REF, lw=1.0)
    axB.plot(conf, acc, "-o", color=C_MODEL, ms=4, lw=1.7)
    axB.set_xlabel("predicted P(binds chosen target)"); axB.set_ylabel("empirical accuracy")
    axB.set_title("B   P(binds) calibration", loc="left", fontweight="bold")
    axB.text(0.05, 0.9, f"ECE = {ece:.3f}", transform=axB.transAxes, fontsize=8, color=C_MODEL)
    axB.set_xlim(0, 1); axB.set_ylim(0, 1)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figF_target_id.{ext}")
    print(f"\nwrote figF_target_id.(pdf|png) to {FIGDIR}")
    for name, b in [("EASY", easy), ("HARD", hard)]:
        print(f"  [{name}] n_test={b['n_test']}  top-1 model {b['model_curve'][0]:.2f} "
              f"shape {b['shape_curve'][0]:.2f} random {b['rand_curve'][0]:.2f}  | "
              f"AUROC model {b['model_auroc']:.3f} shape {b['shape_auroc']:.3f}")
    print(f"  pooled P(binds) ECE = {ece:.3f}")
    hv = "model > shape" if hard["model_curve"][0] > hard["shape_curve"][0] + 0.02 else "model ≈ shape"
    print(f"  HARD verdict: {hv}")


if __name__ == "__main__":
    main()
