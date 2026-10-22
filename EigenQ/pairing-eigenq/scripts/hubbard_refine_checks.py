#!/usr/bin/env python3
"""Final pre-submission pass of 2026-10-21 for paper III (the same checklist as
paper/morefeedback.md for paper I): convergence of the operating-point search
behind Tables 2 and 3 of paper III (notebook Section 14, figure 15).

The joint search of Section 13 used coarse control grids (ramp T in
{8,12,16,24,32,48}/t with n_s in {64,...,512} at eight sites; T in
{32,48,64,128}/t with n_s in {256,512,1024} at sixteen).  Here the ramp and
Path A (shared-flag coarse term, actual propagated coarse state) are
searched on FINE grids of T (steps of 1-2/t at eight sites, 2-4/t at sixteen)
with the step count on a geometric list (ratio 1.15) and then located to the integer by
bracketing + bisection at every target, and Path B on a finer T grid as well.
At eight sites every cycle number M = 0..6 is evaluated with the exact
quantized-time ensemble; at sixteen sites the refinement is restricted to
M = 0, which is complete for the ramp and Path A because one filter cycle
alone (1.8e5 CNOTs on average) costs more than their unfiltered minima at
every target (checked below), so no filtered point can undercut them.  Three
search levels are reported per route and target --- L0 (the Section-13 grid),
L1 (fine grid, geometric n_s), L2 (L1 + integer n_s) --- with the number of
distinct controls within 10 % of the minimum and the Path-A/ramp ratio.  The
coarse state of the 8->16 hierarchy is re-selected as the cheapest REFINED
eight-site unfiltered ramp below 3e-4 ('c') and below 1e-4 ('a').
Self-contained: rebuilds everything it needs from `pairinglib.hubbard`.
Usage: PYTHONPATH=src python3 scripts/hubbard_refine_checks.py [outdir]
"""
import sys as _sys, os, time, json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh
import pairinglib.hubbard as hb

OUT = _sys.argv[1] if len(_sys.argv) > 1 and not _sys.argv[1].startswith("-") else None   # None inside a notebook kernel
if OUT: matplotlib.use("Agg")
t0 = time.time()
c2f = lambda s: np.cos(np.pi * s / 2) ** 2
s2f = lambda s: np.sin(np.pi * s / 2) ** 2
U = 4.0
M_MAX = 6
RATIO = 1.15
TGT8 = (1e-3, 3e-4, 1e-4)
TGT16 = (1e-2, 3e-3, 1e-3, 3e-4)
# the minima of the Section-13 grids (paper III Tables 2 and 3 of 2026-10-20): level L0
L0 = {"8": {"U-ramp": (1.4e4, 2.8e4, 3.7e4), "Path A": (3.7e4, 6.7e4, 1.1e5), "Path B": (1.2e5, 1.6e5, 2.7e5)},
      "16": {"U-ramp": (7.8e4, 7.8e4, 1.6e5, 1.6e5), "Path A": (9.0e4, 1.5e5, 1.5e5, 2.8e5)}}
# the other rows of Table 3 (sixteen sites, 40 schedules; unchanged by this pass), for the figure
T3_OTHER = {"ADAPT-VQE": (4.8e4, 3.9e5, 7.9e5, 1.1e6), "pruned UCCSD + filter": (6.2e4, 4.1e5, 5.9e5, 9.8e5),
            "Path B": (1.1e5, 2.7e5, 5.3e5, 9.2e5)}


def geometric(lo, hi, ratio=RATIO):
    out = [lo]
    while out[-1] < hi:
        out.append(max(out[-1] + 1, int(np.ceil(out[-1] * ratio))))
    return out


class Route:
    """One adiabatic route: per control value a (Phi0, groups, coefs, extra cost) path; a cache of states."""
    def __init__(self, name, step_cnots, GS, path_fn):
        self.name, self.step, self.GS, self.path_fn = name, step_cnots, GS, path_fn
        self.cache = {}

    def state(self, ctrl, T, ns):
        key = (ctrl, float(T), int(ns))
        if key not in self.cache:
            self.cache[key] = self.path_fn(ctrl, float(T), int(ns))
        return self.cache[key]

    def inf0(self, ctrl, T, ns):
        return 1 - abs(np.vdot(self.GS, self.state(ctrl, T, ns))) ** 2

    def cprep(self, ctrl, ns):
        return self.extra(ctrl) + int(ns) * self.step

    def cost(self, ctrl, ns, M, acc, ccyc, reach):
        cp = self.cprep(ctrl, ns)
        return float(cp) if M == 0 else (cp + float(np.sum(reach[:M]) * ccyc)) / acc[M]


def min_steps(route, ctrl, T, target, M, ens, ns_lo, ns_hi):
    """Smallest n_s in [ns_lo, ns_hi] (bracketing by RATIO, then bisection) whose infidelity after M cycles
    is <= target; ens(psi) -> (inf array, acc array) or None for M=0.  Returns (n_s, inf, acc, reach) or None."""
    def ev(ns):
        if M == 0:
            return route.inf0(ctrl, T, ns), None, None
        e = ens(route.state(ctrl, T, ns))
        return e["inf"][M], e["acc"], e["acc"][:-1]
    lo, hi, ns = None, None, ns_lo
    while ns <= ns_hi:
        i, a, r_ = ev(ns)
        if i <= target:
            hi = (ns, i, a, r_); break
        lo = ns; ns = max(ns + 1, int(np.ceil(ns * RATIO)))
    if hi is None:
        return None
    if lo is None:
        return hi
    while hi[0] - lo > 1:
        mid = (lo + hi[0]) // 2
        i, a, r_ = ev(mid)
        if i <= target: hi = (mid, i, a, r_)
        else: lo = mid
    return hi


def levels(route, ctrls, Tgrid, ns_list, targets, ens, ccyc, Ms):
    """L1 cloud over (ctrl, T, n_s, M) and the L2 integer refinement at the targets."""
    cloud = []
    for ctrl in ctrls:
        print(f"    {route.name} control {ctrl}: T grid {Tgrid[0]:g}..{Tgrid[-1]:g}, {len(ns_list)} step counts   [{time.time()-t0:.0f} s]", flush=True)
        for T in Tgrid:
            for ns in ns_list:
                psi = route.state(ctrl, T, ns)
                if ens is None or Ms == (0,):
                    cloud.append((ctrl, T, ns, 0, route.inf0(ctrl, T, ns), float(route.cprep(ctrl, ns))))
                else:
                    e = ens(psi); reach = e["acc"][:-1]
                    for M in Ms:
                        cloud.append((ctrl, T, ns, M, e["inf"][M], route.cost(ctrl, ns, M, e["acc"], ccyc, reach)))
    res = {}
    for tgt in targets:
        ok = [p for p in cloud if p[4] <= tgt]
        if not ok:
            res[tgt] = None; continue
        best1 = min(ok, key=lambda p: p[5])
        per = {}
        for p in ok:
            k = (p[0], p[1], p[3])
            if k not in per or p[5] < per[k][5]: per[k] = p
        refined = []
        for p in sorted(per.values(), key=lambda p: p[5]):
            if p[5] > 1.4 * best1[5]: break
            ctrl, T, ns, M = p[0], p[1], p[2], p[3]
            hit = min_steps(route, ctrl, T, tgt, M, ens, max(ns_list[0], int(ns / RATIO)), ns)
            if hit is not None:
                refined.append((ctrl, T, hit[0], M, hit[1], route.cost(ctrl, hit[0], M, hit[2], ccyc, hit[3])))
        best2 = min(refined, key=lambda p: p[5])
        perc = {}
        for p in refined:
            k = (p[0], p[1])
            if k not in perc or p[5] < perc[k]: perc[k] = p[5]
        flat = [k for k, v in perc.items() if v <= 1.10 * best2[5]]
        res[tgt] = dict(L1=best1, L2=best2, flat=flat)
    return cloud, res


def fmt(p, route):
    ctrl, T, ns, M = p[0], p[1], p[2], p[3]
    if route == "U-ramp": return f"({T:g}, {ns}, {M})"
    if route == "Path A": return f"({T:g}, {ns}, {M}; {ctrl})"
    return f"({ctrl[0]:g}, {T:g}, {ns}, {M}; {ctrl[1]})"


def report(label, key, results, targets, ratio_pair=("Path A", "U-ramp")):
    print(f"\n{label}: best CNOTs per accepted preparation found at three search levels")
    print(f"  {'route':8s} {'target':>7s} | {'L0 Sec.13 grid':>14s} | {'L1 fine grid':>12s} {'point':>22s} | {'L2 +integer':>12s} {'point':>22s} | 1-F | controls within 10%")
    for route, res in results.items():
        for j, tgt in enumerate(targets):
            r_ = res[tgt]
            l0 = L0[key][route][j] if route in L0[key] else float("nan")
            if r_ is None:
                print(f"  {route:8s} {tgt:7.0e} | {l0:14.2e} | not reached"); continue
            fl = r_["flat"]; Ts = sorted(set(k[1] for k in fl))
            print(f"  {route:8s} {tgt:7.0e} | {l0:14.2e} | {r_['L1'][5]:12.4e} {fmt(r_['L1'], route):>22s} | {r_['L2'][5]:12.4e} {fmt(r_['L2'], route):>22s} | {r_['L2'][4]:.2e} | {len(fl)} (T = {Ts[0]:g}..{Ts[-1]:g})")
    a, b = ratio_pair
    print(f"  ratio {a} / {b} at the targets (L2): " + ", ".join(f"{tgt:.0e}: {results[a][tgt]['L2'][5]/results[b][tgt]['L2'][5]:.3f}" if results[a][tgt] and results[b][tgt] else f"{tgt:.0e}: -" for tgt in targets))
    print("  change L0 -> L2: " + ", ".join(f"{route} {tgt:.0e}: {100*(res[tgt]['L2'][5]/L0[key][route][j]-1):+.0f}%" for route, res in results.items() if route in L0[key] for j, tgt in enumerate(targets) if res[tgt]))


def frontier(cloud, tgts, Ms=None):
    return np.array([min((p[5] for p in cloud if p[4] <= tgt and (Ms is None or p[3] in Ms)), default=np.nan) for tgt in tgts])


# ============================================================ eight sites: chain 4 -> 8, (2,2), U = 4t, exact quantized-time ensemble
(Nc, bc), (Nf, bf), blocks, intra = hb.refine_chain(4)
Hf, _ = hb.H_hubbard(Nf, bf, 2, 2, U); Hd = Hf.toarray(); wf, Vf = np.linalg.eigh(Hd); Psih = Vf[:, 0]; E0 = wf[0]
P = hb.prolong(Nc, Nf, blocks, 2, 2); ps = hb.prolonged_state(Hf, P); El, Phi0ex = ps["El"], ps["Phi0"]
Deff, _ = hb.effective_gap(Phi0ex, wf, Vf); SIG = 5.0 / Deff; DT = 0.1
FF = np.linalg.eigh(hb.H_hubbard(Nf, bf, 2, 2, 0.0)[0].toarray())[1][:, 0]
terms = hb.hubbard_terms(Nf, bf, 2, 2, U); Kev, Kod, D = [t_.toarray() for t_ in terms]
tr8 = hb.trotter_evolution(terms, DT, order=2, target=Psih); E_eff = tr8["E_eff"]
st2 = hb.chain_step_cnots(Nf, bf, 2)["cnots"]; ct2 = hb.chain_step_cnots(Nf, bf, 2, controlled=True)["cnots"]
ca8 = hb.conditioned_term_cnots(Nc, bc, n_controls=8); CA_FLAG = ca8["flag"] + ca8["frame"]
CN_GIVENS = 2 * (Nf * (Nf - 1) // 2) * 2; CN_PROL = 2 * Nc * 2
cache8 = {}
ens8 = lambda psi: hb.exact_digital_ensemble(tr8, psi, Psih, SIG, M_MAX, cache=cache8)
e_ = ens8(FF.astype(complex)); ccyc8 = cache8["ns"] @ cache8["wn"] * ct2
print(f"chain 4->8, (2,2), U={U:.0f}: E0 = {E0:.6f}, D_eff = {Deff:.4f}, sigma = {SIG:.2f}/t, dt = {DT}, E_eff - E0 = {E_eff-E0:+.2e}, floor {tr8['floor']:.1e}; "
      f"step {st2} CNOTs (controlled {ct2}), cycle {ccyc8:.0f} CNOTs on average; Givens {CN_GIVENS}, prolongation {CN_PROL}, Path-A coarse term {CA_FLAG} (shared flag + frame)   [{time.time()-t0:.0f} s]")

# the actual four-site coarse states (as in Section 13: the cheapest four-site ramp below 1e-4 ('c') and 2e-5 ('a') on its grid)
N4, b4 = hb.chain(4); H4 = hb.H_hubbard(4, b4, 2, 2, U)[0].toarray(); GS4 = np.linalg.eigh(H4)[1][:, 0]
FF4 = np.linalg.eigh(hb.H_hubbard(4, b4, 2, 2, 0.0)[0].toarray())[1][:, 0]
g4 = [t_.toarray() for t_ in hb.hubbard_terms(4, b4, 2, 2, U)]; c4 = [lambda s: 1.0, lambda s: 1.0, s2f]
st4 = hb.chain_step_cnots(4, b4, 2)["cnots"]; CN_GIVENS4 = 2 * (4 * 3 // 2) * 2
grid4 = {}
for T4 in (8, 12, 16, 24):
    for ns4 in (32, 64, 128, 256, 512):
        psi4 = hb.trotter_refine(FF4, g4, c4, float(T4), ns4, order=2)
        grid4[(T4, ns4)] = dict(psi=psi4, inf=1 - abs(np.vdot(GS4, psi4)) ** 2, cost=CN_GIVENS4 + ns4 * st4)
coarse4 = {}
for tag, thr in (("c", 1e-4), ("a", 2e-5)):
    k4 = min((k for k, d in grid4.items() if d["inf"] < thr), key=lambda k: grid4[k]["cost"])
    v = P @ grid4[k4]["psi"]; coarse4[tag] = dict(Phi0=(v / np.linalg.norm(v)).astype(complex), cprep=grid4[k4]["cost"] + CN_PROL, T=k4[0], ns=k4[1], inf=grid4[k4]["inf"])
    print(f"  four-site coarse state '{tag}': ramp T={k4[0]}, n_s={k4[1]}, 1-F = {grid4[k4]['inf']:.1e}, {coarse4[tag]['cprep']} CNOTs incl. prolongation; prolonged p0 = {abs(np.vdot(Psih, coarse4[tag]['Phi0']))**2:.5f}")

mu = El + 2.0
tl, Ul, cst = hb.galerkin_params("chain", 1.0, U, 2 + 2)
Kc_e, Kc_o, Dc = [t_.toarray() for t_ in hb.hubbard_terms(4, b4, 2, 2, Ul, t=tl)]
# the same product formulas as Section 13 (groups [K_even, K_odd, D], the coarse term through P X P^T) in the
# Kronecker form of `hb.hubbard_groups` / `hb.ProjGroup` (validated against the dense `trotter_refine` to 1e-15; 70x faster)
groups8 = hb.hubbard_groups(Nf, bf, 2, 2, U)
gA = [hb.ProjGroup(P, Dc + (cst - mu) * np.eye(36)), hb.ProjGroup(P, Kc_e), hb.ProjGroup(P, Kc_o), groups8[0], groups8[1], hb.DiagGroup(groups8[2].d - mu)]; cA = [c2f] * 3 + [s2f] * 3
cR = [lambda s: 1.0, lambda s: 1.0, s2f]
_chk = hb.trotter_refine(FF, [Kev, Kod, D], cR, 16.0, 64, order=2)
assert abs(abs(np.vdot(_chk, hb.trotter_path_groups(FF, groups8, cR, 16.0, 64))) - 1) < 1e-12, "Kronecker groups must reproduce the dense product formula"

ramp8 = Route("U-ramp", st2, Psih, lambda ctrl, T, ns: hb.trotter_path_groups(FF, groups8, cR, T, ns)); ramp8.extra = lambda ctrl: CN_GIVENS
pathA8 = Route("Path A", st2 + 2 * CA_FLAG, Psih, lambda ctrl, T, ns: hb.trotter_path_groups(coarse4[ctrl]["Phi0"], gA, cA, T, ns)); pathA8.extra = lambda ctrl: coarse4[ctrl]["cprep"]
def _pathB(ctrl, T, ns):
    Lam, tag = ctrl; cB = [lambda s, L=Lam: 1 + L * np.cos(np.pi * s / 2) ** 2, lambda s: 1.0, lambda s: 1.0]
    return hb.trotter_path_groups(coarse4[tag]["Phi0"], groups8, cB, T, ns)
pathB8 = Route("Path B", st2, Psih, _pathB); pathB8.extra = lambda ctrl: coarse4[ctrl[1]]["cprep"]

ns8 = geometric(32, 768)
Ms8 = tuple(range(M_MAX + 1))
cloudR8, resR8 = levels(ramp8, (None,), np.arange(6, 48.01, 1.0), ns8, TGT8, ens8, ccyc8, Ms8)
print(f"  U-ramp: {len(cloudR8)} grid points evaluated   [{time.time()-t0:.0f} s]")
cloudA8, resA8 = levels(pathA8, ("c", "a"), np.arange(6, 24.01, 1.0), ns8, TGT8, ens8, ccyc8, Ms8)
print(f"  Path A: {len(cloudA8)} grid points evaluated   [{time.time()-t0:.0f} s]")
cloudB8, resB8 = levels(pathB8, [(L, tag) for L in (8.0, 16.0, 32.0) for tag in ("c", "a")], np.arange(12, 32.01, 2.0), geometric(64, 1024), TGT8, ens8, ccyc8, Ms8)
print(f"  Path B: {len(cloudB8)} grid points evaluated   [{time.time()-t0:.0f} s]")
res8 = {"U-ramp": resR8, "Path A": resA8, "Path B": resB8}
report("chain 4->8 (sixteen qubits), exact quantized-time ensemble", "8", res8, TGT8)
print(f"  [{time.time()-t0:.0f} s]")

# ============================================================ sixteen sites: chain 8 -> 16, (2,2), U = 4t, unfiltered (M = 0) refinement
(Nc16, bc16), (Nf16, bf16), blocks16, intra16 = hb.refine_chain(8)
Hf16, _ = hb.H_hubbard(Nf16, bf16, 2, 2, U)
w16, V16 = eigsh(Hf16, k=4, which="SA", tol=1e-10); o = np.argsort(w16); w16, V16 = w16[o], V16[:, o]; GS16 = V16[:, 0]
P16 = hb.prolong(Nc16, Nf16, blocks16, 2, 2, sparse=True); ps16 = hb.prolonged_state(Hf16, P16); Hl16, El16 = ps16["Hl"], ps16["El"]
pool16 = hb._uccsd_pool_kron(Nf16, bf16, 2, 2); FF16 = hb._to_site(pool16, pool16["ref"]).real
groups16 = hb.hubbard_groups(Nf16, bf16, 2, 2, U)
st16 = hb.chain_step_cnots(Nf16, bf16, 2)["cnots"]; ct16 = hb.chain_step_cnots(Nf16, bf16, 2, controlled=True)["cnots"]
ca16d = hb.conditioned_term_cnots(8, bc16, n_controls=16); ca16 = ca16d["flag"] + ca16d["frame"]; CN_GIVENS16 = 2 * (Nf16 * (Nf16 - 1) // 2) * 2; CN_PROL16 = 2 * 8 * 2
SIG16 = 5.0 / hb.effective_gap(ps16["Phi0"], w16, V16)[0]; DT16 = 0.1
nq, wq = hb.quantized_step_distribution(SIG16, DT16); ccyc16 = float(nq @ wq) * ct16
print(f"\nchain 8->16, (2,2), U={U:.0f}: E0 = {w16[0]:.6f}; step {st16} CNOTs (controlled {ct16}), Path-A coarse term {ca16} per application, Givens {CN_GIVENS16}, prolongation {CN_PROL16}; "
      f"sigma = {SIG16:.1f}/t, dt = {DT16}: a cycle costs {ccyc16:.3e} CNOTs on average   [{time.time()-t0:.0f} s]")
# the coarse state: the cheapest REFINED eight-site unfiltered ramp below 3e-4 ('c') and 1e-4 ('a')
coarse8 = {}
for tag, thr in (("c", 3e-4), ("a", 1e-4)):
    pts = [p for p in cloudR8 if p[3] == 0 and p[4] <= thr]
    best = min(pts, key=lambda p: p[5])
    hit = min_steps(ramp8, None, best[1], thr, 0, None, max(32, int(best[2] / RATIO)), best[2])
    psi8 = ramp8.state(None, best[1], hit[0]); v = P16 @ psi8
    coarse8[tag] = dict(Phi0=np.asarray(v, complex).ravel() / np.linalg.norm(v), cprep=int(CN_GIVENS + hit[0] * st2 + CN_PROL16), T=best[1], ns=hit[0], inf=hit[1])
    print(f"  eight-site coarse state '{tag}': ramp T={best[1]:g}, n_s={hit[0]}, 1-F = {hit[1]:.2e}, {coarse8[tag]['cprep']} CNOTs incl. prolongation (Section 13 used {27760+CN_PROL16 if tag=='c' else 36976+CN_PROL16}); prolonged p0 = {abs(np.vdot(GS16, coarse8[tag]['Phi0']))**2:.5f}")
mu16 = El16 + 2.0
Kc_e16, Kc_o16, Dc16 = [t_.toarray() for t_ in hb.hubbard_terms(8, bc16, 2, 2, Ul, t=tl)]
assert abs(Dc16 + (cst - mu16) * np.eye(784) + Kc_e16 + Kc_o16 - (Hl16 - mu16 * np.eye(784))).max() < 1e-10
gA16 = [hb.ProjGroup(P16, Dc16 + (cst - mu16) * np.eye(784)), hb.ProjGroup(P16, Kc_e16), hb.ProjGroup(P16, Kc_o16), groups16[0], groups16[1], hb.DiagGroup(groups16[2].d - mu16)]
cA16 = [c2f] * 3 + [s2f] * 3
ramp_coefs = [lambda s: 1.0, lambda s: 1.0, s2f]
ramp16 = Route("U-ramp", st16, GS16, lambda ctrl, T, ns: hb.trotter_path_groups(FF16, groups16, ramp_coefs, T, ns)); ramp16.extra = lambda ctrl: CN_GIVENS16
pathA16 = Route("Path A", st16 + 2 * ca16, GS16, lambda ctrl, T, ns: hb.trotter_path_groups(coarse8[ctrl]["Phi0"], gA16, cA16, T, ns)); pathA16.extra = lambda ctrl: coarse8[ctrl]["cprep"]
cloudR16, resR16 = levels(ramp16, (None,), np.arange(16, 80.01, 4.0), geometric(64, 768), TGT16, None, ccyc16, (0,))
print(f"  U-ramp (M=0): {len(cloudR16)} grid points evaluated   [{time.time()-t0:.0f} s]")
cloudA16, resA16 = levels(pathA16, ("c", "a"), np.arange(4, 40.01, 2.0), geometric(24, 512), TGT16, None, ccyc16, (0,))
print(f"  Path A (M=0): {len(cloudA16)} grid points evaluated   [{time.time()-t0:.0f} s]")
res16 = {"U-ramp": resR16, "Path A": resA16}
report("chain 8->16 (thirty-two qubits), unfiltered routes", "16", res16, TGT16)
print("  completeness of the M=0 refinement: the expected cost of any filtered point is at least its input plus one average cycle,"
      f" i.e. >= {CN_GIVENS16 + ccyc16:.2e} CNOTs even for the bare determinant, against the unfiltered minima "
      + ", ".join(f"{tgt:.0e}: {min(res16[r][tgt]['L2'][5] for r in res16 if res16[r][tgt]):.2e}" for tgt in TGT16))
print(f"  [{time.time()-t0:.0f} s]")

# ============================================================ figure 15: cost frontiers at eight and sixteen sites
tg8 = np.logspace(np.log10(5e-3), np.log10(5e-5), 25); tg16 = np.logspace(np.log10(2e-2), np.log10(2e-4), 25)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
a_ = ax[0]
for route, cloud, col in (("U-ramp", cloudR8, "k"), ("Path A", cloudA8, "C4"), ("Path B", cloudB8, "C6")):
    a_.loglog(frontier(cloud, tg8, Ms=(0,)), tg8, "--", color=col, lw=1, label=f"{route}, unfiltered ($M=0$)")
    a_.loglog(frontier(cloud, tg8), tg8, "-", color=col, lw=1.8, label=f"{route}, $M=0\\ldots6$ (fine grid)")
    for j, tgt in enumerate(TGT8):
        r_ = res8[route][tgt]
        if r_: a_.loglog([r_["L2"][5]], [tgt], "o", color=col, ms=6, mfc="none")
        a_.loglog([L0["8"][route][j]], [tgt], "s", color=col, ms=5, alpha=0.5)
for tgt in TGT8: a_.axhline(tgt, color="gray", ls=":", lw=0.8)
a_.set_title("chain $4\\to8$, $U=4t$ (exact quantized-time ensemble)"); a_.set_ylim(5e-5, 5e-3)
a_ = ax[1]
for route, cloud, col in (("U-ramp", cloudR16, "k"), ("Path A", cloudA16, "C4")):
    a_.loglog(frontier(cloud, tg16, Ms=(0,)), tg16, "-", color=col, lw=1.8, label=f"{route}, unfiltered (fine grid)")
    for j, tgt in enumerate(TGT16):
        r_ = res16[route][tgt]
        if r_: a_.loglog([r_["L2"][5]], [tgt], "o", color=col, ms=6, mfc="none")
        a_.loglog([L0["16"][route][j]], [tgt], "s", color=col, ms=5, alpha=0.5)
for (name, vals), col, mk in zip(T3_OTHER.items(), ("C3", "C1", "C6"), ("x", "s", "d")):
    a_.loglog(vals, TGT16, mk + ":", color=col, ms=5, label=f"{name} (Section 13 minima)")
for tgt in TGT16: a_.axhline(tgt, color="gray", ls=":", lw=0.8)
a_.set_title("chain $8\\to16$, $U=4t$"); a_.set_ylim(2e-4, 2e-2)
for a_ in ax:
    a_.set_xlabel("best CNOTs per accepted preparation found"); a_.set_ylabel("target infidelity $1-\\mathcal{F}_M$"); a_.grid(alpha=0.3, which="both"); a_.legend(fontsize=7, loc="upper right")
ax[0].text(0.02, 0.02, "open circles: level L2 at the targets; filled squares: the Section-13 grid (L0)", transform=ax[0].transAxes, fontsize=7)
plt.tight_layout()
if OUT:
    plt.savefig(os.path.join(OUT, "fig15.png"), dpi=150); print("wrote", OUT)
    with open(os.path.join(OUT, "hubbard_refine.json"), "w") as f:
        json.dump({key: {route: {f"{tgt:.0e}": dict(L1=[str(r_["L1"][0]), r_["L1"][1], r_["L1"][2], r_["L1"][3], r_["L1"][4], r_["L1"][5]],
                                                   L2=[str(r_["L2"][0]), r_["L2"][1], r_["L2"][2], r_["L2"][3], r_["L2"][4], r_["L2"][5]], flat=[str(c) for c in r_["flat"]])
                                 for tgt, r_ in res.items() if r_} for route, res in results.items()}
                   for key, results in (("8", res8), ("16", res16))}, f, indent=1)
else:
    plt.show()
print(f"[{time.time()-t0:.0f} s]")
