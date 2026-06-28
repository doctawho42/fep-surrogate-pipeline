"""
Phase-0 numerical verification of docs/bar_proofs.tex (v2, corrected).
Throwaway script: reproduce EVERY machine-checked number in the proof sheet,
and disambiguate pymbar 4's two BAR uncertainty methods against the sandwich
and against Monte-Carlo ground truth.

Convention (from the tex):
  Unified coords: x_i = W^f_i (forward),  x_j = -W^r_j (reverse).
  M = ln(n_f/n_r).   p(x; d) = sigma(x - d + M).
  Score S(d) = sum_{j in r} p(x_j) - sum_{i in f} (1 - p(x_i)).
  d_d S = -sum_all p(1-p) = -I(d).   Root S(dF_hat)=0.
"""
import numpy as np
from scipy.optimize import brentq

np.set_printoptions(precision=6, suppress=True)
sigmoid = lambda z: 1.0 / (1.0 + np.exp(-z))

def bar_score(d, xf, xr, M):
    pf = sigmoid(xf - d + M)
    pr = sigmoid(xr - d + M)
    return np.sum(pr) - np.sum(1.0 - pf)

def bar_solve(xf, xr, M):
    f = lambda d: bar_score(d, xf, xr, M)
    return brentq(f, -50, 50, xtol=1e-14, rtol=1e-14)

def info_terms(d, xf, xr, M):
    pf = sigmoid(xf - d + M); pr = sigmoid(xr - d + M)
    If = np.sum(pf * (1 - pf)); Ir = np.sum(pr * (1 - pr))
    return If, Ir, If + Ir

def ok(name, got, exp, tol, unit=""):
    rel = abs(got - exp) / (abs(exp) + 1e-30)
    flag = "PASS" if abs(got - exp) <= tol else "FAIL"
    print(f"  [{flag}] {name}: got={got:.6f} exp={exp:.6f} |Δ|={abs(got-exp):.2e} {unit}")
    return flag == "PASS"

results = []

# ============================================================
print("=" * 64)
print("CHECK 1 — O(1) backward / information-share gradient (Thm backward)")
print("=" * 64)
xf = np.array([0.0, 1.0]); xr = np.array([-0.5]); M = 0.0
d_hat = bar_solve(xf, xr, M)
If, Ir, I = info_terms(d_hat, xf, xr, M)
results.append(ok("dF_hat", d_hat, -0.59571, 1e-4))
results.append(ok("If/I (closed form)", If / I, 0.596825, 1e-5))
# IFT gradient via finite difference: shift all forward x by eps -> dF_hat/dmu_f
eps = 1e-6
d_p = bar_solve(xf + eps, xr, M); d_m = bar_solve(xf - eps, xr, M)
grad_f_fd = (d_p - d_m) / (2 * eps)
results.append(ok("dF/dmu_f (finite diff)", grad_f_fd, If / I, 1e-6))
# (a) derivative wrt the UNIFIED reverse coordinate mean: same IFT form => +Ir/I
grad_xr = (bar_solve(xf, xr + eps, M) - bar_solve(xf, xr - eps, M)) / (2 * eps)
results.append(ok("dF/d(x_r mean) = +Ir/I", grad_xr, Ir / I, 1e-6))
# (b) derivative wrt mu_r, where x_j = -(mu_r+eta) so dx_j/dmu_r = -1 => -Ir/I (theorem)
grad_mur = (bar_solve(xf, xr - eps, M) - bar_solve(xf, xr + eps, M)) / (2 * eps)
results.append(ok("dF/dmu_r (x_j=-(mu_r+eta)) = -Ir/I", grad_mur, -Ir / I, 1e-6))
results.append(ok("|share_f|+|share_r| = 1", abs(If/I) + abs(Ir/I), 1.0, 1e-12))
# curvature d_d S = -I via finite diff of the score
dS_fd = (bar_score(d_hat + eps, xf, xr, M) - bar_score(d_hat - eps, xf, xr, M)) / (2*eps)
results.append(ok("d_d S = -I", dS_fd, -I, 1e-5))

# ============================================================
print("=" * 64)
print("CHECK 2 — sandwich calibration vs naive 1/I (Thm calib)")
print("=" * 64)
# Gaussian work model: forward x~N(mf,s^2), reverse x~N(mr,s^2), mf-mr=s^2, dF=(mf+mr)/2.
# Separation in sigma units = (mf-mr)/s = s.  Choose symmetric means so dF=0.
def gaussian_mc(s, nf=20, nr=20, reps=2000, seed=0):
    rng = np.random.default_rng(seed)
    mf, mr = s**2 / 2, -s**2 / 2          # mf-mr=s^2, dF=(mf+mr)/2=0
    M = np.log(nf / nr)
    dF_true = (mf + mr) / 2
    d_hats, sand_se, naive_se = [], [], []
    pm_mbar_se, pm_bar_se = [], []
    from pymbar.other_estimators import bar as pmbar
    for r in range(reps):
        xf = rng.normal(mf, s, nf); xr = rng.normal(mr, s, nr)
        d = bar_solve(xf, xr, M)
        d_hats.append(d)
        pf = sigmoid(xf - d + M); pr = sigmoid(xr - d + M)
        If = np.sum(pf*(1-pf)); Ir = np.sum(pr*(1-pr)); I = If + Ir
        # B = n_f Var_f[p] + n_r Var_r[p]; use population variance (ddof=0) of p over each sample
        B = nf*np.var(pf, ddof=0) + nr*np.var(pr, ddof=0)
        sand_se.append(np.sqrt(B / I**2))
        naive_se.append(np.sqrt(1.0 / I))
        # pymbar: w_F = W^f = xf ; w_R = W^r = -xr  (since x_r = -W^r)
        rb = pmbar(xf, -xr, compute_uncertainty=True, uncertainty_method='BAR')
        rm = pmbar(xf, -xr, compute_uncertainty=True, uncertainty_method='MBAR')
        pm_bar_se.append(rb['dDelta_f']); pm_mbar_se.append(rm['dDelta_f'])
    d_hats = np.array(d_hats)
    return dict(
        emp_sd=d_hats.std(ddof=1), emp_mean=d_hats.mean(), dF_true=dF_true,
        sand=np.mean(sand_se), naive=np.mean(naive_se),
        pm_bar=np.mean(pm_bar_se), pm_mbar=np.mean(pm_mbar_se),
    )

print("\n  -- headline: separation = 1 sigma (s=1), nf=nr=20, 2000 reps --")
r1 = gaussian_mc(s=1.0, seed=12345)
print(f"    empirical SD(dF_hat) = {r1['emp_sd']:.4f}  (bias mean={r1['emp_mean']:+.4f}, true={r1['dF_true']:.1f})")
print(f"    sandwich se          = {r1['sand']:.4f}  ratio to emp = {r1['sand']/r1['emp_sd']:.3f}")
print(f"    naive 1/sqrt(I)      = {r1['naive']:.4f}  ratio to emp = {r1['naive']/r1['emp_sd']:.3f}")
print(f"    pymbar se ['MBAR']   = {r1['pm_mbar']:.4f}  ratio to emp = {r1['pm_mbar']/r1['emp_sd']:.3f}")
print(f"    pymbar se ['BAR']    = {r1['pm_bar']:.4f}  ratio to emp = {r1['pm_bar']/r1['emp_sd']:.3f}")
results.append(ok("emp SD ~ 0.160", r1['emp_sd'], 0.160, 0.012))
results.append(ok("sandwich/emp ~ 1.0", r1['sand']/r1['emp_sd'], 1.0, 0.06))
results.append(ok("naive/emp ~ 2.2", r1['naive']/r1['emp_sd'], 2.21, 0.20))

print("\n  -- overlap sweep: sep = 1.0, 1.7, 2.4 sigma --")
print("    sep    emp_sd   sand/emp  naive/emp  pmMBAR/emp  pmBAR/emp")
sweep = {}
for s in (1.0, 1.7, 2.4):
    rr = gaussian_mc(s=s, seed=777)
    sweep[s] = rr
    print(f"    {s:.1f}   {rr['emp_sd']:.4f}   {rr['sand']/rr['emp_sd']:.3f}     "
          f"{rr['naive']/rr['emp_sd']:.3f}      {rr['pm_mbar']/rr['emp_sd']:.3f}       {rr['pm_bar']/rr['emp_sd']:.3f}")
# expected naive/emp ~ 2.3, 1.5, 1.2 ; sandwich/emp ~ 1.0 across
results.append(ok("naive/emp @sep1.0 ~ 2.3", sweep[1.0]['naive']/sweep[1.0]['emp_sd'], 2.3, 0.25))
results.append(ok("naive/emp @sep2.4 ~ 1.2", sweep[2.4]['naive']/sweep[2.4]['emp_sd'], 1.2, 0.20))
results.append(ok("sand/emp  @sep2.4 ~ 1.0", sweep[2.4]['sand']/sweep[2.4]['emp_sd'], 1.0, 0.08))

# Which pymbar method equals the sandwich?
print("\n  -- pymbar method vs sandwich (per-rep mean se) --")
for s in (1.0, 1.7, 2.4):
    rr = sweep[s]
    print(f"    sep {s}: sandwich={rr['sand']:.4f}  pmMBAR={rr['pm_mbar']:.4f}  "
          f"pmBAR={rr['pm_bar']:.4f}  | MBAR matches sand: {abs(rr['pm_mbar']-rr['sand'])<0.01}")

# ============================================================
print("=" * 64)
print("CHECK 3 — Fisher–resistance correspondence (Thm fr)")
print("=" * 64)
import networkx as nx
# Triangle nodes 1,2,3; conductances edge(1,2)=1, edge(1,3)=0.5, edge(2,3)=2
nodes = [0, 1, 2]  # 0<->node1, 1<->node2, 2<->node3
edges = {(0,1): 1.0, (0,2): 0.5, (1,2): 2.0}
def laplacian(edges, n=3):
    L = np.zeros((n, n))
    for (i, j), w in edges.items():
        L[i,i]+=w; L[j,j]+=w; L[i,j]-=w; L[j,i]-=w
    return L
def eff_res(L, a, b):
    Lp = np.linalg.pinv(L)
    e = np.zeros(L.shape[0]); e[a]=1; e[b]=-1
    return e @ Lp @ e
L = laplacian(edges)
Omega12 = eff_res(L, 0, 1)
results.append(ok("Omega_12 (series-parallel)", Omega12, 0.714286, 1e-5))
# series-parallel closed form: g_direct=1, g_path = 1/(1/0.5+1/2)=0.4 -> Omega=1/1.4
results.append(ok("Omega_12 = 1/1.4", Omega12, 1/1.4, 1e-9))
# ker L = span{1}
evals, evecs = np.linalg.eigh(L)
v0 = evecs[:, 0]
results.append(ok("smallest eig(L) ~ 0", evals[0], 0.0, 1e-9))
results.append(ok("null vec ∝ 1 (|cos|)", abs(abs(v0 @ np.ones(3)/np.sqrt(3)/np.linalg.norm(v0))), 1.0, 1e-9))
# add edge along c=(1,2) i.e. (0,1) with g=0.7 -> Sherman-Morrison
g = 0.7
Omega_prime_sm = Omega12 / (1 + g * Omega12)
results.append(ok("Omega' = Omega/(1+gOmega)", Omega_prime_sm, 0.476190, 1e-5))
edges2 = dict(edges); edges2[(0,1)] += g
Omega12_prime_direct = eff_res(laplacian(edges2), 0, 1)
results.append(ok("Omega' direct-recompute", Omega12_prime_direct, Omega_prime_sm, 1e-9))
# single edge: Omega = 1/w = V_e
results.append(ok("single-edge Omega = V_e", eff_res(laplacian({(0,1):1.0}, n=2), 0, 1), 1.0, 1e-12))

# ============================================================
print("=" * 64)
print("CHECK 4 — chirality completeness (Thm chiral)")
print("=" * 64)
rng = np.random.default_rng(42)
M_pts = rng.normal(size=(4, 3))  # tetrahedron, 4 points in R^3
def chi(P):
    v1, v2, v3 = P[1]-P[0], P[2]-P[0], P[3]-P[0]
    return np.linalg.det(np.stack([v1, v2, v3], axis=1))
def chi_triple(P):  # triple product form v1.(v2 x v3)
    v1, v2, v3 = P[1]-P[0], P[2]-P[0], P[3]-P[0]
    return v1 @ np.cross(v2, v3)
chiM = chi(M_pts)
results.append(ok("chi = triple product", chi_triple(M_pts), chiM, 1e-12))
# random linear A: chi(A.M) = det(A) chi(M)
A = rng.normal(size=(3,3))
MA = M_pts @ A.T
results.append(ok("chi(A.M) = det(A) chi(M)", chi(MA), np.linalg.det(A)*chiM, 1e-9))
# reflection sigma=diag(-1,1,1): det=-1 -> chi flips sign
S = np.diag([-1.0, 1, 1])
results.append(ok("chi(reflect) = -chi", chi(M_pts @ S.T), -chiM, 1e-12))
# rotation (proper, det=+1) -> chi invariant
th = 0.7; R = np.array([[np.cos(th),-np.sin(th),0],[np.sin(th),np.cos(th),0],[0,0,1]])
results.append(ok("det(R)=+1", np.linalg.det(R), 1.0, 1e-12))
results.append(ok("chi(rotate) = chi", chi(M_pts @ R.T), chiM, 1e-9))
# O(3)-invariant distances identical for M and its mirror; chi differs
D_M = np.array([[np.linalg.norm(M_pts[i]-M_pts[j]) for j in range(4)] for i in range(4)])
Mp = M_pts @ S.T
D_Mp = np.array([[np.linalg.norm(Mp[i]-Mp[j]) for j in range(4)] for i in range(4)])
results.append(ok("pairwise dists M==mirror (max diff)", np.abs(D_M-D_Mp).max(), 0.0, 1e-12))
print(f"  [INFO] sgn chi(M)={np.sign(chiM):+.0f}  sgn chi(mirror)={np.sign(chi(Mp)):+.0f}  -> distances blind, 0o separates")

# ============================================================
print("=" * 64)
n_pass = sum(results); n_tot = len(results)
print(f"SUMMARY: {n_pass}/{n_tot} checks PASS")
print("=" * 64)
