"""Fig Design -- turning the observability map into a perturbation-network design rule.

Theorem 3 / D1 give a per-edge observability certificate: the curl-leverage
``h_e = 1 - w_e * Omega_e`` is the share of edge ``e``'s error that survives into the cycle-closure
residual, and ``delta*_e = sqrt(V_e / h_e)`` is the shift at unit noncentrality -- the resolution
at which closure can see that edge at all. Figure Hodge reports the map and stops: 48 benchmark
edges are bridges (``h_e = 0``) and carry no evidence at any magnitude. This script asks the
prospective question the map implies and never answers. For each benchmark system:

  (i)  how many edges must be ADDED, and between which ligands, so that no bridge remains
       (the network becomes 2-edge-connected), and
  (ii) how many further edges are needed so that the median ``delta*`` over the system's own
       edges falls below a target.

Nothing here runs molecular dynamics. It is graph arithmetic on the released replicate-0
benchmark topology and its reported per-edge standard errors.

PRE-REGISTRATION (fixed before the first run; nothing below was tuned afterwards)
--------------------------------------------------------------------------------
* Systems: the same 48 as Figure Hodge / Figure Lev -- replicate 0 of
  ``data/openfe_replicates/combined_pymbar4_edge_data.csv``, keeping systems with at least 3
  usable edges and at least one independent cycle (``dof >= 1``).
* ``delta*_e = sqrt(V_e / h_e)`` in kcal/mol, the same definition as Figure Hodge and
  Supplementary Table ``tab:auditability``. It is NOT a detection threshold: 80% power at
  alpha=0.05 needs 2.8 to 4.7 times it. The targets below are therefore targets on a resolution
  scale, not on a detectable error.
* Targets swept, stated before running: **median delta* over the system's original edges
  <= 1.0, 0.75 and 0.50 kcal/mol**, on top of step (i). 1.0 sits at the top of the band that
  relative-binding-free-energy decisions turn on; 0.5 is the sharpest resolution a third of the
  benchmark's edges already reach; 0.75 is the midpoint.
* Trajectory metric: the median of ``delta*_e`` over ALL of the system's original edges, with a
  bridge counted at ``delta* = inf``. That is the honest metric for a design rule, because the
  designer is choosing a network for all the perturbations they intend to run, not only for the
  ones that happen to be visible. It differs from the article's tabulated median, which is taken
  over auditable edges only and therefore silently drops the bridges; both are reported.
* Candidate additions: pairs of distinct ligands of that system NOT already directly connected.
  A repeat of an existing perturbation also puts that edge on a cycle and would remove the
  bridge, but it is excluded on purpose: Figure Lval measures per-edge standardized residuals
  correlating at r = +0.30 to +0.42 across independent replicates, so a systematic error that
  reproduces across repeats is exactly what a repeat cannot see. Only a distinct alchemical path
  closes a cycle against reproducible edge-level bias.
* Assumed variance of a new edge -- the one thing the arithmetic cannot know, since an edge's
  variance is not known before it is run: ``V_new = median_e V_e`` over that system's own
  replicate-0 edges. SENSITIVITY: the whole calculation is repeated at ``0.5 * V_new`` and
  ``2 * V_new`` and the change in the answer is reported.
* Minimum for (i): the bridge-block forest's augmentation number ``ceil((d + 2s)/2)`` (d = degree-1
  blocks, s = isolated blocks; 0 if the forest is a single block), the Eswaran--Tarjan lower bound.
  A construction achieving it is emitted and then VERIFIED to leave the graph connected with zero
  bridges, so the reported count is certified minimum rather than asserted. This step involves no
  variance at all and is therefore assumption-free.
* Greedy for (ii): at each step add the candidate pair that minimises the median ``delta*`` over
  the original edges, deterministic tie-break by ligand name. Budget cap 2E added edges -- a
  design that triples the campaign is not a design rule -- and, separately, the exact
  complete-graph value of the metric is computed, which decides whether a target is reachable at
  ANY topology within this design space.
* Reachability floor, exact: ``h_e <= 1`` always, so ``delta*_e >= se_e`` pointwise and the median
  can never fall below the system's median reported standard error, whatever is added.

Run:  PYTHONPATH=src python figs/make_figDesign.py    (or `make figDesign`).  Deterministic.
"""
from __future__ import annotations

import csv
import math
import pathlib
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from bar.leverage import bridges, curl_leverage  # noqa: E402
from bar.qc import gls_network  # noqa: E402
from paperstyle import (  # noqa: E402
    INK, MUTED, OURS, REF, check_min_type, figsize, finish, legend, panel, tint,
    use_paper_style,
)

DATA = ROOT / "data" / "openfe_replicates" / "combined_pymbar4_edge_data.csv"
DOC = ROOT / "docs" / "results_figDesign.md"
TABLE = ROOT / "docs" / "tab_design.tex"

TARGETS = (1.0, 0.75, 0.50)          # kcal/mol, median delta* over the original edges
VAR_SCALES = (0.5, 1.0, 2.0)         # sensitivity on the assumed variance of a new edge
BUDGET_MULT = 2                      # cap on added edges, as a multiple of E
FLAGGED = {"brd4", "bace", "faah", "cdk8", "hif2a", "p38"}   # the Fig L flag set


# ---------------------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------------------
def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def edge_val(r, k):
    cD = _f(r[f"complex_repeat_{k}_DG (kcal/mol)"])
    cd = _f(r[f"complex_repeat_{k}_dDG (kcal/mol)"])
    sD = _f(r[f"solvent_repeat_{k}_DG (kcal/mol)"])
    sd = _f(r[f"solvent_repeat_{k}_dDG (kcal/mol)"])
    if any(math.isnan(v) for v in (cD, cd, sD, sd)):
        return None
    return cD - sD, math.sqrt(cd * cd + sd * sd)


def load_systems():
    by = defaultdict(list)
    for r in csv.DictReader(DATA.open()):
        by[r["system name"]].append(r)
    out = {}
    for name, rows in sorted(by.items()):
        edges = [(r["ligand_A"], r["ligand_B"], *edge_val(r, 0)) for r in rows if edge_val(r, 0)]
        if len(edges) < 3:
            continue
        fit = gls_network(edges)
        if fit.dof < 1:
            continue
        out[name] = edges
    return out


# ---------------------------------------------------------------------------------------
# (i) bridge elimination: minimum augmentation to 2-edge-connectivity
# ---------------------------------------------------------------------------------------
def bridge_block_forest(nodes, pairs):
    """Return (block_of_node, blocks, forest_adj) for the bridge-block forest of ``pairs``."""
    br = bridges([(a, b, 0.0, 1.0) for a, b in pairs])
    adj = defaultdict(list)
    for i, (a, b) in enumerate(pairs):
        if i in br:
            continue
        adj[a].append(b)
        adj[b].append(a)
    block_of, blocks = {}, []
    for n in nodes:
        if n in block_of:
            continue
        comp, stack = [], [n]
        block_of[n] = len(blocks)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in block_of:
                    block_of[v] = len(blocks)
                    stack.append(v)
        blocks.append(sorted(comp, key=str))
    fadj = defaultdict(set)
    for i, (a, b) in enumerate(pairs):
        if i in br:
            fadj[block_of[a]].add(block_of[b])
            fadj[block_of[b]].add(block_of[a])
    for k in range(len(blocks)):
        fadj.setdefault(k, set())
    return block_of, blocks, fadj


def augmentation_lower_bound(nodes, pairs):
    """Eswaran--Tarjan bridge-connectivity augmentation number ceil((d + 2s)/2)."""
    _bo, blocks, fadj = bridge_block_forest(nodes, pairs)
    if len(blocks) == 1:
        return 0
    d = sum(1 for k in range(len(blocks)) if len(fadj[k]) == 1)
    s = sum(1 for k in range(len(blocks)) if len(fadj[k]) == 0)
    return -(-(d + 2 * s) // 2)


def _rep_pair(blocks, ka, kb, degree, existing):
    """Pick a representative ligand in each of two blocks, avoiding a duplicate of an
    existing perturbation. Highest degree first (the block's hub), name tie-break."""
    ca = sorted(blocks[ka], key=lambda n: (-degree[n], str(n)))
    cb = sorted(blocks[kb], key=lambda n: (-degree[n], str(n)))
    for u in ca:
        for v in cb:
            if u != v and frozenset((u, v)) not in existing:
                return u, v
    return ca[0], cb[0]          # only reachable if the two blocks are already saturated


def _tree_leaf_order(fadj, root):
    """Leaves of one forest component in DFS discovery order, rooted at ``root``."""
    order, seen, stack = [], {root}, [root]
    while stack:
        u = stack.pop()
        if len(fadj[u]) <= 1 and u != root:
            order.append(u)
        for v in sorted(fadj[u], reverse=True):
            if v not in seen:
                seen.add(v)
                stack.append(v)
    if len(fadj[root]) <= 1:
        order.append(root)
    return order


def augment_to_2ec(nodes, pairs):
    """Add edges until the graph is connected and bridgeless; return the added ligand pairs."""
    degree = defaultdict(int)
    for a, b in pairs:
        degree[a] += 1
        degree[b] += 1
    cur = list(pairs)
    existing = {frozenset(p) for p in cur}
    added = []
    while True:
        _bo, blocks, fadj = bridge_block_forest(nodes, cur)
        if len(blocks) == 1:
            break
        # forest components
        seen, comps = set(), []
        for k in range(len(blocks)):
            if k in seen:
                continue
            comp, stack = [], [k]
            seen.add(k)
            while stack:
                u = stack.pop()
                comp.append(u)
                for v in fadj[u]:
                    if v not in seen:
                        seen.add(v)
                        stack.append(v)
            comps.append(sorted(comp))
        if len(comps) > 1:
            # phase 1: chain the components, always spending a slot (a degree<=1 block)
            c0, c1 = comps[0], comps[1]
            ka = min((k for k in c0 if len(fadj[k]) <= 1), default=c0[0])
            kb = min((k for k in c1 if len(fadj[k]) <= 1), default=c1[0])
        else:
            # phase 2: one tree; pair leaves across the halves of a DFS leaf order
            comp = comps[0]
            root = max(comp, key=lambda k: (len(fadj[k]), -k))
            leaves = _tree_leaf_order(fadj, root)
            r = (len(leaves) + 1) // 2
            ka, kb = leaves[0], leaves[r % len(leaves)]
        u, v = _rep_pair(blocks, ka, kb, degree, existing)
        cur.append((u, v))
        existing.add(frozenset((u, v)))
        added.append((u, v))
        degree[u] += 1
        degree[v] += 1
    return added


# ---------------------------------------------------------------------------------------
# (ii) leverage arithmetic
# ---------------------------------------------------------------------------------------
def graph_components(nodes, pairs):
    """Connected components of the ligand graph, as lists of nodes."""
    adj = defaultdict(list)
    for a, b in pairs:
        adj[a].append(b)
        adj[b].append(a)
    seen, comps = set(), []
    for n in nodes:
        if n in seen:
            continue
        comp, stack = [], [n]
        seen.add(n)
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in adj[u]:
                if v not in seen:
                    seen.add(v)
                    stack.append(v)
        comps.append(sorted(comp, key=str))
    return comps


def bridgefree_plan(nodes, pairs):
    """Added pairs that leave NO BRIDGE, without joining separate components.

    An edge is auditable as soon as it lies on a cycle, and the closure fit already absorbs
    one offset per component, so bridge freedom -- not 2-edge-connectivity -- is the minimum
    that buys coverage. Each component is augmented on its own; the counts are summed.
    """
    out = []
    for comp in graph_components(nodes, pairs):
        members = set(comp)
        sub = [(a, b) for a, b in pairs if a in members]
        out.extend(augment_to_2ec(comp, sub))
    return out


def resistance_matrix(n, idx_edges):
    """Full effective-resistance matrix of the weighted graph (list of (i, j, w))."""
    L = np.zeros((n, n))
    for i, j, w in idx_edges:
        L[i, i] += w
        L[j, j] += w
        L[i, j] -= w
        L[j, i] -= w
    P = np.linalg.pinv(L)
    d = np.diag(P)
    return d[:, None] + d[None, :] - 2.0 * P


def delta_star(R, ia, ib, V):
    """delta*_e = sqrt(V_e / h_e) with h_e = 1 - Omega_e / V_e, on the given resistance matrix."""
    h = 1.0 - R[ia, ib] / V
    out = np.full(len(V), np.inf)
    good = h > 1e-12
    out[good] = np.sqrt(V[good] / h[good])
    return out


def greedy_trajectory(nodes, cur_pairs, cur_V, orig_slice, V_new, budget, stop_at):
    """Greedily add the pair that most reduces the median delta* over the original edges.

    Returns (trajectory, added_pairs). ``trajectory[k]`` is the median after k additions.
    """
    n = len(nodes)
    idx = {x: i for i, x in enumerate(nodes)}
    ia = np.array([idx[a] for a, _ in cur_pairs])
    ib = np.array([idx[b] for _, b in cur_pairs])
    V = np.array(cur_V, dtype=float)
    oa, ob = ia[orig_slice].copy(), ib[orig_slice].copy()
    oV = V[orig_slice].copy()
    have = {frozenset((idx[a], idx[b])) for a, b in cur_pairs}
    cand = np.array([(i, j) for i in range(n) for j in range(i + 1, n)
                     if frozenset((i, j)) not in have], dtype=int)
    w_new = 1.0 / V_new
    added, traj = [], []
    R = resistance_matrix(n, list(zip(ia, ib, 1.0 / V, strict=True)))
    traj.append(float(np.median(delta_star(R, oa, ob, oV))))
    while len(added) < budget and traj[-1] > stop_at and len(cand):
        A, B = cand[:, 0], cand[:, 1]
        term = 0.5 * (R[np.ix_(ob, A)] - R[np.ix_(oa, A)] - R[np.ix_(ob, B)] + R[np.ix_(oa, B)])
        Rp = R[oa, ob][:, None] - w_new * term ** 2 / (1.0 + w_new * R[A, B])[None, :]
        h = 1.0 - Rp / oV[:, None]
        with np.errstate(divide="ignore", invalid="ignore"):
            dd = np.where(h > 1e-12, np.sqrt(oV[:, None] / np.maximum(h, 1e-300)), np.inf)
        med = np.median(dd, axis=0)
        k = int(np.argmin(med))
        i, j = int(cand[k, 0]), int(cand[k, 1])
        added.append((nodes[i], nodes[j]))
        ia = np.append(ia, i)
        ib = np.append(ib, j)
        V = np.append(V, V_new)
        cand = np.delete(cand, k, axis=0)
        R = resistance_matrix(n, list(zip(ia, ib, 1.0 / V, strict=True)))
        traj.append(float(np.median(delta_star(R, oa, ob, oV))))
    return traj, added


def complete_graph_value(nodes, cur_pairs, cur_V, orig_slice, V_new):
    """Median delta* over the original edges when every missing ligand pair is added."""
    n = len(nodes)
    idx = {x: i for i, x in enumerate(nodes)}
    ia = [idx[a] for a, _ in cur_pairs]
    ib = [idx[b] for _, b in cur_pairs]
    V = list(cur_V)
    have = {frozenset((i, j)) for i, j in zip(ia, ib, strict=True)}
    for i in range(n):
        for j in range(i + 1, n):
            if frozenset((i, j)) not in have:
                ia.append(i)
                ib.append(j)
                V.append(V_new)
    R = resistance_matrix(n, list(zip(ia, ib, 1.0 / np.array(V), strict=True)))
    oa = np.array(ia)[orig_slice]
    ob = np.array(ib)[orig_slice]
    oV = np.array(V)[orig_slice]
    return float(np.median(delta_star(R, oa, ob, oV)))


# ---------------------------------------------------------------------------------------
# per-system driver
# ---------------------------------------------------------------------------------------
def design_system(name, edges, var_scale=1.0, want_traj=True):
    nodes = sorted({e[0] for e in edges} | {e[1] for e in edges}, key=str)
    pairs = [(a, b) for a, b, _, _ in edges]
    V = np.array([se ** 2 for _, _, _, se in edges])
    V_new = float(np.median(V)) * var_scale
    fit = gls_network(edges)

    h0 = curl_leverage(edges)
    n_bridge = int(np.sum(h0 <= 1e-9))
    d0 = np.full(len(edges), np.inf)
    aud = h0 > 1e-9
    d0[aud] = np.sqrt(V[aud] / h0[aud])
    med_auditable = float(np.median(d0[aud])) if aud.any() else math.nan
    med_all = float(np.median(d0))

    lb = augmentation_lower_bound(nodes, pairs)
    add2 = augment_to_2ec(nodes, pairs)
    comps = graph_components(nodes, pairs)
    addb = bridgefree_plan(nodes, pairs)
    lb_b = sum(augmentation_lower_bound(c, [(a, b) for a, b in pairs if a in set(c)])
               for c in comps)
    assert len(addb) == lb_b, (name, len(addb), lb_b)
    assert len(bridges([(a, b, 0.0, 1.0) for a, b in pairs + addb])) == 0, name
    # certify: the construction leaves the graph connected and bridgeless, and matches the bound
    aug_pairs = pairs + add2
    assert len(bridges([(a, b, 0.0, 1.0) for a, b in aug_pairs])) == 0, name
    _bo, blk, _fa = bridge_block_forest(nodes, aug_pairs)
    assert len(blk) == 1, name
    assert len(add2) == lb, (name, len(add2), lb)

    cur_V = list(V) + [V_new] * len(add2)
    orig_slice = np.arange(len(edges))
    # cross-check the incremental leverage against the released two-path curl_leverage
    R2 = resistance_matrix(len(nodes), [(nodes.index(a), nodes.index(b), 1.0 / v)
                                        for (a, b), v in zip(aug_pairs, cur_V, strict=True)])
    idx = {x: i for i, x in enumerate(nodes)}
    oa = np.array([idx[a] for a, _ in pairs])
    ob = np.array([idx[b] for _, b in pairs])
    h_inc = 1.0 - R2[oa, ob] / V
    h_ref = curl_leverage([(a, b, 0.0, math.sqrt(v))
                           for (a, b), v in zip(aug_pairs, cur_V, strict=True)])[:len(edges)]
    dev = float(np.max(np.abs(h_inc - h_ref)))

    med_after2 = float(np.median(delta_star(R2, oa, ob, V)))
    floor_pointwise = float(np.median(np.sqrt(V)))
    complete = complete_graph_value(nodes, aug_pairs, cur_V, orig_slice, V_new)

    budget = BUDGET_MULT * len(edges)
    traj, extra = ([med_after2], []) if not want_traj else greedy_trajectory(
        nodes, aug_pairs, cur_V, orig_slice, V_new, budget, min(TARGETS))
    need = {}
    for t in TARGETS:
        hit = next((k for k, m in enumerate(traj) if m <= t), None)
        if hit is not None:
            need[t] = hit
        elif complete > t:
            need[t] = "unreachable"
        else:
            need[t] = "over-budget"
    return {
        "system": name, "N": len(nodes), "E": len(edges), "dof": fit.dof,
        "n_bridge": n_bridge, "a2": len(add2), "lb": lb, "add2_pairs": add2,
        "c": len(comps), "ab": len(addb), "addb_pairs": addb,
        "med_auditable": med_auditable, "med_all": med_all, "med_after2": med_after2,
        "floor": floor_pointwise, "complete": complete, "V_new": V_new,
        "traj": traj, "need": need, "extra_pairs": extra, "leverage_dev": dev,
        "flagged": name in FLAGGED,
    }


# ---------------------------------------------------------------------------------------
# independent certification of the minima
# ---------------------------------------------------------------------------------------
def brute_force_minimality(systems, res, max_subsets=300_000):
    """Exhaustively confirm no set of ``a-1`` added edges achieves what ``a`` achieves.

    The Eswaran--Tarjan bound already proves the counts minimum, and the construction is
    checked against it; this is the belt-and-braces check that the bound was implemented
    right. Every system whose search space is small enough is checked, for both targets.
    """
    from itertools import combinations
    done, skipped = 0, []
    for r in res:
        edges = systems[r["system"]]
        nodes = sorted({e[0] for e in edges} | {e[1] for e in edges}, key=str)
        pairs = [(a, b) for a, b, _, _ in edges]
        have = {frozenset(p) for p in pairs}
        cand = [(u, v) for i, u in enumerate(nodes) for v in nodes[i + 1:]
                if frozenset((u, v)) not in have]
        for goal, want_connected in ((r["ab"], False), (r["a2"], True)):
            if goal == 0:
                continue
            k, n_sub = goal - 1, 1
            for i in range(goal - 1):
                n_sub = n_sub * (len(cand) - i) // (i + 1)
            if n_sub > max_subsets:
                skipped.append((r["system"], goal, n_sub))
                continue
            for combo in combinations(cand, k):
                trial = pairs + list(combo)
                if bridges([(a, b, 0.0, 1.0) for a, b in trial]):
                    continue
                if want_connected and len(graph_components(nodes, trial)) > 1:
                    continue
                raise AssertionError(f"{r['system']}: {k} added edges already suffice")
            done += 1
    return done, skipped


# ---------------------------------------------------------------------------------------
# figure
# ---------------------------------------------------------------------------------------
def _path_xy(r):
    """The design path: (added edges / E, median delta*) from as-built to the greedy end."""
    xs = [0.0]
    ys = [r["med_all"]]
    if r["a2"]:
        xs.append(r["a2"] / r["E"])
        ys.append(r["med_after2"])
    for k, m in enumerate(r["traj"]):
        if k == 0 and r["a2"]:
            continue
        xs.append((r["a2"] + k) / r["E"])
        ys.append(m)
    return xs, ys


def _log_ticks(axis, ticks):
    """Label exactly these decadal ticks on a log axis; silence the minor ones."""
    axis.set_major_locator(FixedLocator(ticks))
    axis.set_major_formatter(NullFormatter())
    axis.set_ticklabels([f"{t:g}" for t in ticks])
    axis.set_minor_locator(NullLocator())
    axis.set_minor_formatter(NullFormatter())


def _cum(values, grid, n):
    return np.array([sum(1 for v in values if v is not None and v <= x) / n for x in grid])


def make_figure(res, sens):
    use_paper_style()
    fig, axes = plt.subplots(1, 3, figsize=figsize(3, 2.95))
    n = len(res)
    grid = np.linspace(0.0, 2.0, 401)

    # ---- A: the design path, per system --------------------------------------------
    ax = axes[0]
    for r in res:
        if r["flagged"]:
            continue
        xs, ys = _path_xy(r)
        ax.plot(xs, ys, color=MUTED, lw=0.6, alpha=0.55, zorder=1)
    for r in res:
        if not r["flagged"]:
            continue
        xs, ys = _path_xy(r)
        ax.plot(xs, ys, color=OURS, lw=1.4, zorder=3)
        ax.plot([xs[0]], [ys[0]], marker="o", ms=3.0, color=OURS, zorder=4)
    for t in TARGETS:
        ax.axhline(t, color=REF, ls=":", lw=0.9, zorder=0.5)
    ax.set_yscale("log")
    ax.set_xlim(0, 2.0)
    ax.set_ylim(0.11, 5.2)
    ax.set_xlabel(r"edges added ($\times E$)")
    ax.set_ylabel(r"median $\delta^{*}$ (kcal/mol)")
    _log_ticks(ax.yaxis, [0.2, 0.5, 1.0, 2.0, 4.0])
    ax.annotate("targets", xy=(1.97, 1.06), ha="right", va="bottom", fontsize=7.5, color=REF)
    ax.plot([], [], color=OURS, lw=1.4, label="flagged")
    ax.plot([], [], color=MUTED, lw=0.6, label="other")
    legend(ax, loc="upper right")
    panel(ax, "A", "the cost of resolution", "each system's design path")

    # ---- B: how many systems reach each milestone, against budget --------------------
    ax = axes[1]
    ax.plot(grid, _cum([r["a2"] / r["E"] for r in res], grid, n), color=OURS, lw=1.6,
            label="2-edge-connected")
    styles = {1.0: ((4, 1.6), 0.30), 0.75: ((2, 1.4), 0.50), 0.50: ((1, 1.2), 0.68)}
    for t in TARGETS:
        dash, sh = styles[t]
        vals = [(r["a2"] + r["need"][t]) / r["E"] if isinstance(r["need"][t], int) else None
                for r in res]
        ax.plot(grid, _cum(vals, grid, n), color=tint(OURS, sh), lw=1.4, dashes=dash,
                label=rf"median $\delta^{{*}}\leq{t:g}$")
    lo = _cum([(r["a2"] + r["need"][1.0]) / r["E"] if isinstance(r["need"][1.0], int) else None
               for r in sens[0.5]], grid, n)
    hi = _cum([(r["a2"] + r["need"][1.0]) / r["E"] if isinstance(r["need"][1.0], int) else None
               for r in sens[2.0]], grid, n)
    ax.fill_between(grid, np.minimum(lo, hi), np.maximum(lo, hi), color=tint(OURS, 0.30),
                    alpha=0.30, lw=0, zorder=0.8)
    ax.set_xlim(0, 2.0)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel(r"budget: edges added ($\times E$)")
    ax.set_ylabel(f"fraction of the {n} systems")
    legend(ax, loc="lower right")
    panel(ax, "B", "what a budget buys", r"band: new-edge variance $\times$0.5 to $\times$2")

    # ---- C: the ceiling topology cannot pass ----------------------------------------
    ax = axes[2]
    for flag, colour, size, lab in ((False, MUTED, 12, "other systems"),
                                    (True, OURS, 20, "flagged systems")):
        pts = [r for r in res if r["flagged"] == flag]
        ax.scatter([r["med_auditable"] for r in pts], [r["complete"] for r in pts],
                   s=size, color=colour, edgecolors="none", label=lab,
                   zorder=3 if flag else 2)
    ax.plot([0.1, 4.2], [0.1, 4.2], color=REF, ls=":", lw=1.0, zorder=0.5)
    ax.axhline(1.0, color=REF, ls="--", lw=0.9, zorder=0.5)
    n_un = sum(1 for r in res if r["complete"] > 1.0)
    ax.annotate(f"{n_un} systems cannot reach 1.0\nat any topology", xy=(0.118, 1.07),
                ha="left", va="bottom", fontsize=7.5, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.11, 5.2)
    ax.set_ylim(0.11, 5.2)
    ticks = [0.2, 0.5, 1.0, 2.0, 4.0]
    _log_ticks(ax.xaxis, ticks)
    _log_ticks(ax.yaxis, ticks)
    ax.set_xlabel(r"median $\delta^{*}$ as built")
    ax.set_ylabel(r"median $\delta^{*}$, every pair run")
    legend(ax, loc="upper left")
    panel(ax, "C", "the floor under the rule", "topology buys the gap to $y=x$")

    offenders = check_min_type(fig)
    if offenders:
        raise SystemExit(f"type below the floor: {offenders}")
    finish(fig, "figDesign_network_design_rule")


# ---------------------------------------------------------------------------------------
# table + record
# ---------------------------------------------------------------------------------------
def _cell(v):
    if v == "unreachable":
        return "--"
    if v == "over-budget":
        return r"$>\!2E$"
    return f"${v}$"


def _sortkey(r):
    v = r["need"][1.0]
    return (-r["a2"], 0 if isinstance(v, int) else 1, v if isinstance(v, int) else 0, r["system"])


def write_table(res):
    rows = sorted(res, key=_sortkey)
    half = (len(rows) + 1) // 2
    left, right = rows[:half], rows[half:]
    body = []
    for i in range(half):
        cells = []
        for r in (left[i], right[i] if i < len(right) else None):
            if r is None:
                cells.append(" & " * 7 + " ")
                continue
            cells.append(
                r"\texttt{%s}%s & $%d$ & $%d$ & $%d$ & $%d$ & %s & %s & %s" % (
                    r["system"].replace("_", r"\_"),
                    r"$^\dagger$" if r["flagged"] else "",
                    r["E"], r["n_bridge"], r["ab"], r["a2"],
                    _cell(r["need"][1.0]), _cell(r["need"][0.75]), _cell(r["need"][0.50]),
                )
            )
        body.append(" & ".join(cells) + r" \\")
    head = (r"system & $E$ & br & $a_b$ & $a_2$ & $+1.0$ & $+0.75$ & $+0.5$ & "
            r"system & $E$ & br & $a_b$ & $a_2$ & $+1.0$ & $+0.75$ & $+0.5$ \\")
    tot_a2 = sum(r["a2"] for r in res)
    tot_ab = sum(r["ab"] for r in res)
    tot_br = sum(r["n_bridge"] for r in res)
    tot_e = sum(r["E"] for r in res)
    caption = (
        r"\textbf{The cost, in added edges, of a network that can be audited.} For the $48$ "
        r"benchmark systems carrying at least one independent cycle: $E$ edges and the number of "
        r"bridge edges (br), then two certified-minimum augmentation counts. $a_b$ is the "
        r"smallest number of new ligand pairs that leaves NO bridge, so that every edge lies on "
        r"a cycle and carries some evidence; it treats a system's connected components "
        r"separately, which is all auditability requires, since the closure fit already absorbs "
        r"one offset per component. $a_2$ additionally joins those components, leaving the whole "
        r"network 2-edge-connected. Both are pure topology and assume nothing about the new "
        r"edges. The last three columns are the FURTHER edges a greedy design adds, on top of "
        r"$a_2$, to bring the median $\delta^\ast$ over that system's own $E$ edges below $1.0$, "
        r"$0.75$ and $0.5$ kcal/mol, each new edge assumed to carry that system's median reported "
        r"variance; unlike $a_b$ and $a_2$ these are an achievable cost and not a proven minimum, "
        r"since the design is greedy. \texttt{--} marks a target no topology can reach, because "
        r"even measuring every remaining ligand pair leaves the median above it, and $>\!2E$ a "
        r"target not reached within the explored budget. Across the benchmark $%d$ added edges "
        r"remove all $%d$ bridges and $%d$ leave every system 2-edge-connected, $%.1f\%%$ and "
        r"$%.1f\%%$ of the $%d$ edges already run. $\delta^\ast$ is the shift at "
        r"unit noncentrality and not a detection threshold; $80\%%$ power at $\alpha=0.05$ needs "
        r"$2.8$ to $4.7$ times it, so these are targets on a resolution scale. Values are "
        r"replicate-0 and come from the same run of the released analysis code, "
        r"\texttt{make figDesign}. A dagger marks the six systems the cycle-closure test flags."
    ) % (tot_ab, tot_br, tot_a2, 100.0 * tot_ab / tot_e, 100.0 * tot_a2 / tot_e, tot_e)
    tex = "\n".join([
        r"\begin{table}[tp]\centering",
        r"\caption{" + caption + "}",
        r"\label{tab:design}",
        r"\footnotesize",
        r"\setlength{\tabcolsep}{2.4pt}",
        r"\begin{tabular}{@{}lrrrrrrr@{\hspace{6pt}}lrrrrrrr@{}}",
        r"\toprule",
        head,
        r"\midrule",
        *body,
        r"\bottomrule",
        r"\end{tabular}",
        r"\end{table}",
        "",
    ])
    TABLE.write_text(tex)
    print(f"wrote {TABLE}")


def _fmt_need(v):
    return "unreachable" if v == "unreachable" else (">2E" if v == "over-budget" else str(v))


def write_record(res, sens, certify):
    n = len(res)
    tot_e = sum(r["E"] for r in res)
    tot_br = sum(r["n_bridge"] for r in res)
    tot_a2 = sum(r["a2"] for r in res)
    tot_ab = sum(r["ab"] for r in res)
    n_with_br = sum(1 for r in res if r["n_bridge"])
    n_multi = sum(1 for r in res if r["c"] > 1)
    worse = [r for r in res if r["a2"] and r["med_after2"] > r["med_auditable"] + 1e-9]
    reach = {t: [r for r in res if isinstance(r["need"][t], int)] for t in TARGETS}
    unreach = {t: [r for r in res if r["need"][t] == "unreachable"] for t in TARGETS}
    over = {t: [r for r in res if r["need"][t] == "over-budget"] for t in TARGETS}
    dev = max(r["leverage_dev"] for r in res)

    def cost(t):
        """(n already there, n needing work, median and max cost among those needing work)."""
        tot = [(r, r["a2"] + r["need"][t]) for r in reach[t]]
        work = [(r, c) for r, c in tot if c > 0]
        if not work:
            return len(tot), 0, 0, 0, 0.0
        return (len(tot) - len(work), len(work), int(np.median([c for _, c in work])),
                int(max(c for _, c in work)),
                float(np.median([c / r["E"] for r, c in work])))

    def plural(k, word):
        return f"{k} {word}" + ("" if k == 1 else "s")

    lines = []
    A = lines.append
    A("# Results -- Fig Design: from the observability map to a network design rule")
    A("")
    A("**Figure:** `figs/figDesign_network_design_rule.{pdf,png}` | **Table:**")
    A("`docs/tab_design.tex` | **Reproduce:** `make figDesign` (or `PYTHONPATH=src python")
    A("figs/make_figDesign.py`). Deterministic; there is no randomness anywhere in this script.")
    A("Data: replicate 0 of `data/openfe_replicates/combined_pymbar4_edge_data.csv`, the released")
    A("OpenFE IndustryBenchmarks2024 set. Reuses `src/bar/leverage.py` (`curl_leverage`,")
    A("`bridges`) and `src/bar/qc.py` (`gls_network`). **No molecular dynamics is run and no")
    A("number already in the article changes.** This is graph arithmetic on the benchmark")
    A("topology and its reported per-edge standard errors, turning the observability map into the")
    A("prospective rule it implies but never states.")
    A("")
    A("![Fig Design](../figs/figDesign_network_design_rule.png)")
    A("")
    A("## The question")
    A("")
    A("Theorem 3 / D1 give the per-edge observability certificate `h_e = 1 - w_e*Omega_e` and the")
    A("resolution `delta*_e = sqrt(V_e/h_e)`, the shift at unit noncentrality. Fig Hodge reports")
    A("the map and stops: 48 benchmark edges are bridges, `h_e = 0`, carrying no evidence at any")
    A("magnitude. That map is the only prospective output the work has, and reporting a bound is")
    A("not the same as acting on it. So: **how many edges must be added, and between which")
    A("ligands, so that (i) no bridge remains and (ii) the median `delta*` falls below a target?**")
    A("")
    A("## Pre-registration (fixed before the first run)")
    A("")
    A("- **Targets swept, stated before running:** median `delta*` over the system's own edges")
    A("  `<= 1.0`, `0.75`, `0.50` kcal/mol, on top of the topology step.")
    A("- **Trajectory metric:** median of `delta*_e` over ALL of a system's original edges, a")
    A("  bridge counted at `delta* = inf`. The article's tabulated median is taken over auditable")
    A("  edges only and so drops the bridges silently; both are reported below, and the")
    A("  difference matters (see \"What the topology step does not buy\").")
    A("- **Candidates:** pairs of distinct ligands not already directly connected. A repeat of an")
    A("  existing perturbation would also put that edge on a cycle, and is excluded on purpose:")
    A("  Fig Lval measures per-edge standardized residuals correlating at r = +0.30 to +0.42")
    A("  across independent replicates, so an error that reproduces across repeats is exactly")
    A("  what a repeat cannot see. Only a distinct alchemical path closes a cycle against it.")
    A("- **Assumed variance of a new edge** (the one thing the arithmetic cannot know, since an")
    A("  edge's variance is not known before it is run): that system's own median reported")
    A("  variance, swept over 0.5x, 1x and 2x as the pre-registered sensitivity.")
    A("- **Budget cap** 2E added edges; separately, the exact complete-graph value decides")
    A("  whether a target is reachable at ANY topology.")
    A("- **Minimum for (i):** the Eswaran--Tarjan bound `ceil((d+2s)/2)` on the bridge-block")
    A("  forest, with a construction emitted and then verified bridgeless and connected, so the")
    A("  count is certified minimum rather than asserted.")
    A("")
    A("## (i) The topology step: assumption-free, and cheap")
    A("")
    A("Two different things get conflated by the phrase \"remove the bridges\", and they cost")
    A("different amounts, so both are reported.")
    A("")
    A("- `a_b`: the minimum added edges that leave **no bridge**. Each connected component is")
    A("  augmented on its own. This is all that auditability needs: an edge is observable as soon")
    A("  as it lies on a cycle, and the closure fit already absorbs one offset per component.")
    A("- `a_2`: the minimum added edges that additionally **join the components**, leaving the")
    A("  whole network 2-edge-connected. This is the quantity usually meant by")
    A("  \"2-edge-connected\", and it buys cross-component comparability as well as coverage.")
    A("")
    A("```")
    A(f"{n} systems, {tot_e} edges, {tot_br} bridges on {n_with_br} systems, "
      f"{n_multi} systems in more than one component")
    n_ab = sum(1 for r in res if r["ab"])
    n_a2 = sum(1 for r in res if r["a2"])
    A(f"a_b: {tot_ab} added edges leave no bridge  "
      f"({100.0*tot_ab/tot_e:.1f}% of the {tot_e} edges already run, on {n_ab} systems)")
    A(f"a_2: {tot_a2} added edges leave every system 2-edge-connected  "
      f"({100.0*tot_a2/tot_e:.1f}%, on {n_a2} systems)")
    A(f"both constructions == the Eswaran-Tarjan lower bound on {n}/{n} systems; every")
    A("  augmented network verified to have zero bridges, and every a_2 network verified connected")
    A(f"max |h incremental - curl_leverage| on the designed networks: {dev:.1e}")
    A(f"brute force: on {certify[0]} of the {certify[0] + len(certify[1])} "
      f"(system, target) pairs needing any edge, EVERY subset of one fewer edge was")
    A("  enumerated and none achieves the goal; the rest have too large a space to enumerate")
    A("```")
    A("")
    A(f"- **{tot_ab} added edges remove all {tot_br} bridges**, {100.0*tot_ab/tot_e:.1f}% more")
    A("  perturbations than the benchmark already ran; one edge can put several bridges on cycles")
    A(f"  at once, which is why {tot_ab} edges suffice for {tot_br} bridges. Making every system")
    A(f"  2-edge-connected as well costs {tot_a2}, i.e. {100.0*tot_a2/tot_e:.1f}%.")
    A("- Neither number assumes anything about the new edges. They are properties of the")
    A("  topology alone, and are certified minimum: the construction's count is checked against")
    A("  the Eswaran--Tarjan lower bound and the result is verified bridgeless. Independently")
    A("  of that bound, the script brute-forces the claim wherever the search space allows: for")
    A(f"  {certify[0]} of the {certify[0] + len(certify[1])} (system, target) pairs that need any")
    A("  edge at all, every subset of one fewer edge is enumerated and none achieves the goal.")
    A(f"- The most expensive single system needs "
      f"{max(r['ab'] for r in res)} for bridge freedom "
      f"({', '.join(r['system'] for r in res if r['ab'] == max(x['ab'] for x in res))}) and "
      f"{max(r['a2'] for r in res)} for 2-edge-connectivity "
      f"({', '.join(r['system'] for r in res if r['a2'] == max(x['a2'] for x in res))}).")
    A("- The per-system counts are in `docs/tab_design.tex`; the ligand pairs themselves are")
    A("  listed at the end of this record.")
    A("")
    A("## What the topology step does not buy")
    A("")
    A("Removing the bridges does not make a network sharper. It makes previously invisible edges")
    A("visible **at poor resolution**: an edge freshly placed on one long cycle has small `h_e`,")
    A("and `delta*_e = se_e / sqrt(h_e)` is correspondingly large. Measured over all of a")
    A("system's edges, the median `delta*` after the minimum augmentation is *worse* than the")
    A(f"article's auditable-edges-only median on {len(worse)} of the "
      f"{sum(1 for r in res if r['a2'])} systems that need any augmentation, because the bridges")
    A("the article's median excludes are now inside it. Step (i) buys coverage and step (ii)")
    A("buys resolution; they are separate purchases, and only (ii) costs a budget worth arguing")
    A("about.")
    A("")
    A("## (ii) Reaching a resolution target")
    A("")
    A("Greedy design on top of the `a_2` network: at each step add the ligand pair that most")
    A("reduces the median `delta*` over the system's original edges. Unlike (i), these counts are")
    A("an **achievable cost, not a proven minimum** -- the greedy is a heuristic and an optimal")
    A("design could be cheaper.")
    A("")
    A("```")
    for t in TARGETS:
        already, k_work, med_e, max_e, med_c = cost(t)
        A(f"median delta* <= {t:<4}: reachable on {len(reach[t])}/{n} systems  "
          f"| unreachable at any topology {len(unreach[t])} | over the 2E budget {len(over[t])}")
        A(f"                     {already} of those already there as built; the other {k_work} "
          f"need a median of {med_e} added edges ({med_c:.2f}E), at most {max_e}")
    A("```")
    A("")
    A("- The binding constraint is not topology, it is the edges' own standard errors. Because")
    A("  `h_e <= 1` always, `delta*_e >= se_e` pointwise: no design can push an edge's resolution")
    A("  below its own reported standard error. Measuring **every** remaining ligand pair still")
    A(f"  leaves {len(unreach[1.0])} of the {n} systems above 1.0 kcal/mol and "
      f"{len(unreach[0.5])} above 0.5.")
    A(f"- Where a target is reachable and not already met, the cost is real but not prohibitive: "
      f"a median of {plural(cost(1.0)[2], 'edge')} ({cost(1.0)[4]:.2f}E) for 1.0 kcal/mol and "
      f"{plural(cost(0.5)[2], 'edge')} ({cost(0.5)[4]:.2f}E) for 0.5.")
    A("- The systems that cannot be brought to 1.0 at any topology are small, or have large")
    A("  reported standard errors, or both. Two floors are at work: the pointwise one just")
    A("  stated, and a size floor, since a complete graph on `N` ligands with equal weights gives")
    A("  `h_e = 1 - 2/N`, so a small network cannot reach high leverage however densely it is")
    A("  wired.")
    A("")
    A("  | system | N | E | median se | median delta* as built | every pair run |")
    A("  |---|---|---|---|---|---|")
    for r in sorted(unreach[1.0], key=lambda r: -r["complete"]):
        A(f"  | `{r['system']}` | {r['N']} | {r['E']} | {r['floor']:.2f} | "
          f"{r['med_auditable']:.2f} | {r['complete']:.2f} |")
    A("")
    A("## Sensitivity to the assumed variance of a new edge")
    A("")
    A("An added edge's variance is not known before it is run, so everything in section (ii)")
    A("assumes one; section (i) assumes nothing and is untouched. Sweeping the assumption over a")
    A("factor of four, at the 1.0 kcal/mol target:")
    A("")
    A("| assumed new-edge variance | systems reaching 1.0 | median cost among those needing work |")
    A("|---|---|---|")
    for sc in VAR_SCALES:
        rr = sens[sc]
        ok = [(x, x["a2"] + x["need"][1.0]) for x in rr if isinstance(x["need"][1.0], int)]
        work = [(x, c) for x, c in ok if c > 0]
        med = int(np.median([c for _, c in work])) if work else 0
        medc = float(np.median([c / x["E"] for x, c in work])) if work else 0.0
        A(f"| {sc:g}x the system's median | {len(ok)}/{n} | {med} edges ({medc:.2f}E), "
          f"{len(work)} systems |")
    A("")
    base = {r["system"]: r["need"][1.0] for r in res}
    moved = []
    for sc in (0.5, 2.0):
        for x in sens[sc]:
            a, b = base[x["system"]], x["need"][1.0]
            if isinstance(a, int) and isinstance(b, int) and abs(a - b) > 1:
                moved.append((x["system"], sc, str(a), str(b)))
            elif isinstance(a, int) != isinstance(b, int):
                moved.append((x["system"], sc, _fmt_need(a), _fmt_need(b)))
    if moved:
        A("The verdict (reachable or not) and the count are stable to within a single edge on")
        A(f"{2*n - len(moved)} of the {2*n} system-by-scale comparisons. The exceptions:")
        A("")
        for name, sc, a, b in sorted(moved):
            A(f"- `{name}`: {a} at 1x, {b} at {sc:g}x")
        A("")
        A("These are the systems sitting just at a target, where a greedy path either finds a")
        A("cheap route or does not; they are also the clearest evidence that the greedy counts")
        A("are an upper bound on the cost and not the cost.")
    else:
        A("No system's answer moves by more than a single edge across that four-fold range.")
    A("")
    A("## What this can and cannot say")
    A("")
    A("- **Can:** the two topology counts and the ligand pairs achieving them are exact,")
    A("  certified minimum and assumption-free. The reachability verdict is exact too: it")
    A("  compares the target against the complete-graph value, which no topology can beat.")
    A("- **Cannot:** the edge counts in (ii) assume a variance for edges nobody has run, and the")
    A("  sensitivity above is the honest width of that assumption; they are also greedy, so they")
    A("  are an achievable cost rather than a minimum.")
    A("- **Cannot:** which ligand inside a block to attach is chosen here by graph degree. The")
    A("  graph is indifferent between the members of a block; chemistry is not, and a pair the")
    A("  graph likes may be a perturbation nobody can run. Treat the counts as the deliverable")
    A("  and the named pairs as one legal choice out of many.")
    A("- **Cannot:** `delta*` is the shift at unit noncentrality, not a detection threshold. 80%")
    A("  power at alpha=0.05 needs 2.8 to 4.7 times it, so a design reaching a median `delta*` of")
    A("  1.0 kcal/mol has a median edge that is *resolvable* at 1.0, not one that *detects* a 1.0")
    A("  kcal/mol error. The targets are on a resolution scale, and the same caveat the article")
    A("  attaches to its tabulated `delta*` attaches here.")
    A("- **Cannot:** none of this is validated prospectively. It is arithmetic on one benchmark's")
    A("  topology, and what it buys is observability, which Fig Lev already showed does not")
    A("  predict where reproducible error actually lands. The rule says where the instrument")
    A("  could see, not where the error will be.")
    A("")
    A("## The sentence a practitioner can act on")
    A("")
    A("> When planning a perturbation network, first spend the few edges that leave no bridge:")
    A(f"> across this benchmark that is {tot_ab} edges for {tot_br} bridges, "
      f"{100.0*tot_ab/tot_e:.1f}% more than was already")
    A("> run, it assumes nothing, and without it those edges carry no evidence at any magnitude.")
    A("> Then design for a target `delta*`, remembering that it can never fall below an edge's")
    A(f"> own standard error: a median of 1.0 kcal/mol costs a further "
      f"{plural(cost(1.0)[2], 'edge')} ({cost(1.0)[4]:.2f}E) on the")
    A(f"> {cost(1.0)[1]} systems that need any, and is out of reach at any topology on "
      f"{len(unreach[1.0])} of {n}, where the")
    A("> sampling budget and not the graph is what has to change.")
    A("")
    A("## The proposed edges")
    A("")
    A("Ligand pairs whose measurement leaves the network with no bridge (`a_b`), and, where they")
    A("differ, the pairs that additionally leave it 2-edge-connected (`a_2`). The counts are")
    A("minimum; the specific ligands are one achieving choice out of many, picked by graph degree")
    A("and not by chemistry.")
    A("")
    A("```")
    for r in sorted(res, key=lambda r: r["system"]):
        if not r["a2"] and not r["ab"]:
            continue
        A(f"{r['system']}  ({plural(r['n_bridge'], 'bridge')}, "
          f"{plural(r['c'], 'component')};  a_b = {r['ab']}, a_2 = {r['a2']})")
        if r["addb_pairs"]:
            A("  no bridge:")
            for u, v in r["addb_pairs"]:
                A(f"    {u}  --  {v}")
        if r["addb_pairs"] != r["add2_pairs"]:
            A("  2-edge-connected:")
            for u, v in r["add2_pairs"]:
                A(f"    {u}  --  {v}")
    A("```")
    DOC.write_text("\n".join(lines) + "\n")
    print(f"wrote {DOC}")


def main():
    systems = load_systems()
    res = [design_system(name, edges) for name, edges in systems.items()]
    sens = {s: ([design_system(name, edges, var_scale=s) for name, edges in systems.items()]
                if s != 1.0 else res) for s in VAR_SCALES}
    print(f"[design] {len(res)} systems, {sum(r['E'] for r in res)} edges, "
          f"{sum(r['n_bridge'] for r in res)} bridges, {sum(r['a2'] for r in res)} added for "
          f"2-edge-connectivity")
    for t in TARGETS:
        ok = [r for r in res if isinstance(r["need"][t], int)]
        print(f"[design] median delta* <= {t:g}: {len(ok)}/{len(res)} reachable, "
              f"{sum(1 for r in res if r['need'][t] == 'unreachable')} unreachable at any topology")
    certify = brute_force_minimality(systems, res)
    print(f"[design] brute-forced minimality on {certify[0]} (system, target) pairs; "
          f"{len(certify[1])} search spaces too large to enumerate")
    make_figure(res, sens)
    write_table(res)
    write_record(res, sens, certify)


if __name__ == "__main__":
    main()
