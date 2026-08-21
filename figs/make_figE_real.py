"""Fig E (real-molecule illustration) -- chirality completeness on a concrete chiral drug.

The theorem is numbered in the manuscript, not here: the figure text deliberately carries no
theorem number, because an earlier renumbering (the conservation law was inserted ahead of it)
left a stale "Theorem 4" baked into this graphic while the LaTeX caption's \\ref stayed correct.

The synthetic Fig E proves the chirality-completeness statement on designed 4-point
tetrahedra. Here we instantiate it on a *real* molecule: thalidomide, the textbook case
where the two enantiomers have different biological activity. We embed one enantiomer in
3D (RDKit), reflect it to the exact mirror geometry, and read out both:

  * every O(3)-invariant ("even") readout is IDENTICAL on the pair -- all pairwise
    distances and all pairwise dot products (the Gram matrix) match to machine precision,
    so a norm-and-dot-product readout is provably enantiomer-blind;
  * the parity-odd 0o pseudoscalar chi = det[p1-p0, p2-p0, p3-p0] flips sign;
  * RDKit confirms the reflection swaps the CIP label R <-> S.

Illustrative (one molecule), reproducing the chirality-completeness proof on real chemistry.
Run:  PYTHONPATH=src python figs/make_figE_real.py    (or `make figEreal`). Deterministic.
"""
from __future__ import annotations

import io
import pathlib
import sys
import warnings

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")
ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from bar.chiral import even_features, signed_volume  # noqa: E402

FIGDIR = pathlib.Path(__file__).resolve().parent
EMBED_SEED = 0xC0FFEE
# thalidomide; one stereocentre in the glutarimide ring. R = [C@@H], S (mirror) = [C@H].
SMILES_R = "O=C1CC[C@@H](N2C(=O)c3ccccc3C2=O)C(=O)N1"
SMILES_S = "O=C1CC[C@H](N2C(=O)c3ccccc3C2=O)C(=O)N1"


def prep3d(smiles: str, seed: int = EMBED_SEED):
    """Embed a 3D conformer; return heavy-atom coordinates (N x 3)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    AllChem.MMFFOptimizeMolecule(m)
    return Chem.RemoveHs(m).GetConformer().GetPositions()


def cip_of(coords):
    """CIP label of the (single) stereocentre implied by a heavy-atom geometry."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(SMILES_R))
    AllChem.EmbedMolecule(m, randomSeed=EMBED_SEED)
    mh = Chem.RemoveHs(m)
    conf = mh.GetConformer()
    for i, p in enumerate(coords):
        conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))
    Chem.AssignStereochemistryFrom3D(mh)
    Chem.AssignStereochemistry(mh, cleanIt=True, force=True)
    return next((a.GetPropsAsDict().get("_CIPCode") for a in mh.GetAtoms()
                 if a.HasProp("_CIPCode")), "?")


def draw2d(smiles: str, px: int = 340):
    """Clean 2D structural depiction with wedge bonds + CIP annotation (RDKit cairo)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdCIPLabeler
    from rdkit.Chem.Draw import rdMolDraw2D
    m = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(m)
    rdCIPLabeler.AssignCIPLabels(m)
    d = rdMolDraw2D.MolDraw2DCairo(px, px)
    d.drawOptions().addStereoAnnotation = True
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    return plt.imread(io.BytesIO(d.GetDrawingText()), format="png")


def readouts(coords):
    V = coords - coords.mean(0)
    return dict(dist=even_features(coords), gram=V @ V.T, chi=signed_volume(coords))


def main() -> None:
    mpl.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                         "savefig.dpi": 300, "savefig.bbox": "tight"})
    coords = prep3d(SMILES_R)
    mirror = coords * np.array([-1.0, 1.0, 1.0])       # reflect x -> exact enantiomer geometry
    rR, rS = readouts(coords), readouts(mirror)
    dist_dmax = float(np.abs(rR["dist"] - rS["dist"]).max())
    gram_dmax = float(np.abs(rR["gram"] - rS["gram"]).max())
    cipR, cipS = cip_of(coords), cip_of(mirror)
    n_at, n_d = len(coords), len(rR["dist"])

    fig, (axR, axS, axT) = plt.subplots(1, 3, figsize=(8.6, 3.1),
                                        gridspec_kw={"width_ratios": [1, 1, 1.15]})
    for ax, smi, lab, chi in [(axR, SMILES_R, cipR, rR["chi"]), (axS, SMILES_S, cipS, rS["chi"])]:
        ax.imshow(draw2d(smi)); ax.set_axis_off()
        ax.set_title(f"({lab})-thalidomide   $\\chi={chi:+.2f}$", fontsize=10, fontweight="bold")
    axS.set_title(f"({cipS})-thalidomide (mirror)   $\\chi={rS['chi']:+.2f}$",
                  fontsize=10, fontweight="bold")

    axT.set_axis_off()
    txt = (f"chirality completeness on a real molecule\n({n_at} heavy atoms)\n\n"
           f"O(3)-invariant readouts, M vs mirror M$'$:\n"
           f"  all {n_d} pairwise distances\n     max|$\\Delta$| = {dist_dmax:.0e}\n"
           f"  all pairwise dot products (Gram)\n     max|$\\Delta$| = {gram_dmax:.0e}\n"
           f"  $\\Rightarrow$ an even readout is blind\n\n"
           f"parity-odd 0o pseudoscalar $\\chi$:\n"
           f"  $\\chi$(M) = {rR['chi']:+.3f},  $\\chi$(M$'$) = {rS['chi']:+.3f}\n"
           f"  sum = {rR['chi']+rS['chi']:+.0e}   (exact sign flip)\n\n"
           f"RDKit CIP: {cipR} $\\rightarrow$ {cipS}")
    axT.text(0.0, 1.0, txt, va="top", ha="left", fontsize=8.2, family="monospace",
             transform=axT.transAxes)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(FIGDIR / f"figE_chirality_real.{ext}")
    print(f"wrote figE_chirality_real.(pdf|png) to {FIGDIR}")
    print(f"\nthalidomide: {n_at} heavy atoms, {n_d} pairwise distances")
    print(f"  even (distance) readout    max|Delta(M,M')| = {dist_dmax:.2e}")
    print(f"  dot-product (Gram) readout max|Delta(M,M')| = {gram_dmax:.2e}")
    print(f"  0o pseudoscalar  chi(M)={rR['chi']:+.4f}  chi(M')={rS['chi']:+.4f}  "
          f"sum={rR['chi']+rS['chi']:+.2e}")
    print(f"  RDKit CIP: {cipR} -> {cipS}")


if __name__ == "__main__":
    main()
