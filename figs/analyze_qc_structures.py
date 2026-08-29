"""Structural context of the cycle-closure QC flags (Fig L, supplement).

Retrospective, data-internal structural check: the two systems the paper names -- BACE1
(``protonation-sensitive``) and BRD4 (``buried-water set``) -- are read straight from public
co-crystal structures to show the a-priori-known hard feature sits in the binding site:

  * BACE1 (PDB 4DJW): the catalytic aspartic-acid dyad hydrogen-bonds the inhibitor. Its
    protonation state is the textbook difficulty for aspartic-protease free energies.
  * BRD4 BD1 (PDB 3MXF, JQ1): the conserved acetyl-lysine-pocket water network and the
    KAc-recognition asparagine sit against the ligand -- the known BRD4 FEP challenge
    (Aldeghi 2016; Filippakopoulos 2010).

This is structural CONTEXT for why these systems are hard, not a per-edge causal claim (the
QC test is edge-level and blind to node-consistent bias). The flags coincide with the
independently-documented hard features. Public PDBs cached in data/pdb/.
Run: PYTHONPATH=src python figs/analyze_qc_structures.py   (or `make qcstruct`). Deterministic.
"""
from __future__ import annotations

import pathlib
import sys
import urllib.request
from collections import defaultdict

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import proj3d

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from paperstyle import (  # noqa: E402
    CATEGORICAL, INK, figsize, finish, legend, panel, tint, use_paper_style,
)

PDBDIR = ROOT / "data" / "pdb"
FIGDIR = pathlib.Path(__file__).resolve().parent
ION = {"HOH", "SO4", "PO4", "GOL", "EDO", "NA", "CL", "K", "MG", "ZN", "ACT", "DMS", "PEG",
       "BME", "TRS", "CO3", "IOD", "FMT", "MPD", "NO3", "CA"}
# ANNOTATION FIGURE -- see "The annotation-figure exemption" in paperstyle. The colours here
# annotate chemical OBJECTS, not methods: no estimator, no baseline and no method series appears
# in either panel, so the semantic palette does not apply and CATEGORICAL is drawn from instead.
# Three annotated objects, three obviously different hues -- the catalytic aspartate dyad (red),
# the crystallographic waters (blue, the structural-biology convention) and the recognition
# asparagine (green) -- with the ligand in neutral grey as context. The previous pass forced the
# semantic palette on: the dyad and the waters both went blue and the asparagine became a pale-
# blue square indistinguishable from a blue water circle at print size.
C_DYAD = CATEGORICAL[1]    # #EE6677 red   -- BACE1 catalytic Asp dyad
C_WATER = CATEGORICAL[0]   # #4477AA blue  -- crystallographic waters
C_ASN = CATEGORICAL[2]     # #228833 green -- BRD4 KAc-recognition Asn140
C_LIG = CATEGORICAL[5]     # #BBBBBB grey  -- the ligand, neutral context
VIEW = dict(elev=16, azim=40)   # one viewpoint for both panels


def fetch(pdb_id: str) -> pathlib.Path:
    p = PDBDIR / f"{pdb_id}.pdb"
    if not p.exists():
        PDBDIR.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pdb_id}.pdb", p)
    return p


def parse(path):
    at = []
    for line in open(path):
        if line[:6] in ("ATOM  ", "HETATM"):
            at.append(dict(rec=line[:6].strip(), resn=line[17:20].strip(), chain=line[21],
                           seq=int(line[22:26]), name=line[12:16].strip(),
                           xyz=np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])))
    return at


def ligand(atoms):
    grp = defaultdict(list)
    for a in atoms:
        if a["rec"] == "HETATM" and a["resn"] not in ION:
            grp[(a["resn"], a["chain"], a["seq"])].append(a)
    return max(grp.values(), key=len) if grp else []


def resacts(atoms, resn):
    r = defaultdict(list)
    for a in atoms:
        if a["rec"] == "ATOM" and a["resn"] == resn:
            r[(a["chain"], a["seq"])].append(a)
    return r


def min_dist(res_atoms, lig_pts):
    P = np.array([a["xyz"] for a in res_atoms])
    return float(np.min(np.linalg.norm(P[:, None] - lig_pts[None], axis=2)))


def _scatter(ax, pts, **kw):
    # depth shading is off on purpose: it fades a marker toward white by how far back it sits,
    # which under this article's colour rule would read as a tint, i.e. as a second meaning.
    ax.scatter(pts[:, 0], pts[:, 1], pts[:, 2], depthshade=False, **kw)


def _frame(ax, *clouds, pad=0.6, fill=0.86):
    """Fill the panel with one undistorted cube around the content, clipping nothing.

    A 3D axes left alone scales each axis to its own data range, which both distorts the
    coordinates and spends most of the panel on empty cube -- the reason this figure used to be
    a small blob of atoms in a field of white. One side length is used for all three axes, so an
    angstrom is the same length whichever way it points, and the side is the content's own
    extent plus ``pad`` angstrom of air.

    The magnification is then measured rather than guessed. At a given zoom the content projects
    to a box that the cube scales about the panel centre, so the zoom whose box occupies ``fill``
    of the panel follows from one measurement; two passes are taken because the projection is
    only approximately linear in the zoom. A hand-picked constant does not survive this: the
    value that filled the BACE1 panel cut a water in half in the BRD4 one.
    """
    P = np.vstack(clouds)
    lo, hi = P.min(0), P.max(0)
    mid = (lo + hi) / 2.0
    half = float(np.max(hi - lo)) / 2.0 + pad
    for i, setlim in enumerate((ax.set_xlim, ax.set_ylim, ax.set_zlim)):
        setlim(mid[i] - half, mid[i] + half)
    ax.set_axis_off()
    ax.view_init(**VIEW)
    zoom = 1.0
    for _ in range(2):
        ax.set_box_aspect((1, 1, 1), zoom=zoom)
        ax.figure.canvas.draw()
        d = ax.transData.transform(
            [proj3d.proj_transform(*q, ax.get_proj())[:2] for q in P])
        bb = ax.get_window_extent()
        rx = max(float(np.abs(d[:, 0] - (bb.x0 + bb.x1) / 2.0).max()), 1e-9)
        ry = max(float(np.abs(d[:, 1] - (bb.y0 + bb.y1) / 2.0).max()), 1e-9)
        zoom *= fill * min(bb.width / 2.0 / rx, bb.height / 2.0 / ry)
    ax.set_box_aspect((1, 1, 1), zoom=zoom)


def label3d(ax, xyz, text, offset, color):
    """Name a site with the text set clear of its own markers and tied back by a hairline.

    The name used to be drawn AT the residue's centroid, so it lay across the very markers it
    was naming. Here the anchor is projected to the panel's 2D space and the string is offset
    from it in points, which keeps the offset the same size whatever the cube's scale is.
    """
    x2, y2, _ = proj3d.proj_transform(*xyz, ax.get_proj())
    ax.annotate(
        text, xy=(x2, y2), xycoords="data", xytext=offset, textcoords="offset points",
        ha="left" if offset[0] > 0 else "right", va="center",
        fontsize=7.5, color=color, zorder=6, annotation_clip=False, clip_on=False,
        arrowprops=dict(arrowstyle="-", lw=0.6, color=tint(color, 0.45),
                        shrinkA=1.0, shrinkB=3.0),
    )


def main():
    use_paper_style()
    fig = plt.figure(figsize=figsize(2, 3.15))

    # ---- BACE1: catalytic aspartic dyad ----
    at = parse(fetch("4DJW"))
    lig = ligand(at); LP = np.array([a["xyz"] for a in lig])
    asp = sorted(((min_dist(v, LP), k, v) for k, v in resacts(at, "ASP").items()))[:2]
    axA = fig.add_subplot(1, 2, 1, projection="3d")
    dyad = [np.array([a["xyz"] for a in v if a["name"] in ("OD1", "OD2", "CG")]) for _, _, v in asp]
    _scatter(axA, LP, c=C_LIG, s=11, lw=0, alpha=0.85, label=f"ligand ({lig[0]['resn']})")
    _scatter(axA, np.vstack(dyad), c=C_DYAD, s=30, edgecolors=INK, linewidths=0.5,
             label="catalytic Asp dyad")
    _frame(axA, LP, np.vstack(dyad))
    panel(axA, "A", "BACE1 (PDB 4DJW)",
          f"catalytic Asp dyad, {asp[0][0]:.1f}/{asp[1][0]:.1f} Å to ligand")
    legend(axA, loc="upper left", bbox_to_anchor=(0.0, 0.0), scatterpoints=1,
           ncol=2, columnspacing=1.3, handletextpad=0.4)

    # ---- BRD4: conserved KAc-pocket waters + Asn140 ----
    at2 = parse(fetch("3MXF"))
    lig2 = ligand(at2); LP2 = np.array([a["xyz"] for a in lig2])
    wat = np.array([a["xyz"] for a in at2 if a["resn"] == "HOH"])
    dwat = np.min(np.linalg.norm(wat[:, None] - LP2[None], axis=2), axis=1)
    near = wat[dwat <= 4.5]
    asn = sorted(((min_dist(v, LP2), k, v) for k, v in resacts(at2, "ASN").items()))[0]
    an = np.array([a["xyz"] for a in asn[2] if a["name"] in ("ND2", "OD1", "CG")])
    axB = fig.add_subplot(1, 2, 2, projection="3d")
    _scatter(axB, LP2, c=C_LIG, s=11, lw=0, alpha=0.85, label=f"ligand ({lig2[0]['resn']})")
    _scatter(axB, near, c=C_WATER, s=26, edgecolors=INK, linewidths=0.5,
             label=f"waters ≤ 4.5 Å ({len(near)})")
    # a different object from the waters, so a different hue rather than a tint of theirs
    # the legend entry stays bare -- the "(KAc)" qualifier is already in the panel subtitle, and
    # spelling it out here pushed the third column past the right edge of the canvas.
    _scatter(axB, an, c=C_ASN, s=34, marker="s", edgecolors=INK, linewidths=0.5,
             label=f"Asn{asn[1][1]}")
    _frame(axB, LP2, near, an)
    panel(axB, "B", f"BRD4 BD1 (PDB 3MXF, {lig2[0]['resn']})",
          f"{len(near)} waters ≤ 4.5 Å; Asn{asn[1][1]} {asn[0]:.1f} Å (KAc)")
    legend(axB, loc="upper left", bbox_to_anchor=(0.0, 0.0), scatterpoints=1,
           ncol=3, columnspacing=1.3, handletextpad=0.4)

    # two reserved bands only -- the headings above, one legend line below, so that neither
    # ever lands on the point cloud; the cubes themselves run to the panel edges.
    fig.subplots_adjust(left=0.005, right=0.995, bottom=0.115, top=0.855, wspace=0.02)
    fig.canvas.draw()   # freeze the projections before anything is placed against them
    for (_, (_, seq), _), cc, off in zip(asp, dyad, ((-15.0, 9.0), (17.0, 6.0)), strict=True):
        label3d(axA, cc.mean(0), f"Asp{seq}", off, C_DYAD)
    label3d(axB, an.mean(0), f"Asn{asn[1][1]}", (17.0, -9.0), C_ASN)

    finish(fig, "figQC_structures", layout=None)
    print(f"\n[BACE1 4DJW] ligand {lig[0]['resn']}; catalytic aspartic dyad: "
          + ", ".join(f"{v[0]['resn']}{k[1]} @ {d:.2f} A" for d, k, v in asp))
    print(f"[BRD4 3MXF]  ligand {lig2[0]['resn']}; waters <=3.5A: {(dwat<=3.5).sum()}, "
          f"<=4.5A: {(dwat<=4.5).sum()}; conserved Asn{asn[1][1]} @ {asn[0]:.2f} A (KAc recognition)")


if __name__ == "__main__":
    main()
