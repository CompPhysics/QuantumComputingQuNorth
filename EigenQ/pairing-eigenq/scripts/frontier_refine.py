#!/usr/bin/env python3
"""Final revision of 2026-10-21 (paper/morefeedback.md, Sec. 2): convergence of
the operating-point search behind Table tab:frontier of paper I (notebook
section EXT:REFINE, figure 25).

The joint search of `scripts/frontier_checks.py` used a coarse control grid
(T in steps of 2, three shifts, eight step counts).  Here both adiabatic
routes --- the HF ramp D -> H and the shift-controlled refinement k/2 -> k ---
are searched on a FINE grid of the evolution time T (steps of 0.25 or 0.5) and
of the shift mu_buf (ten / nine values), and the integer step count n_s is
located by bracketing + bisection at every target, for every cycle number
M = 0..6 with the exact averaged channel of the digital filter.  Three
levels of search are reported per route and target so that the convergence
of the "best cost found" is visible:

  L0  the coarse grid of Table tab:frontier (2026-10-20),
  L1  the fine (T, mu_buf) grid with n_s on a geometric list (ratio 1.15),
  L2  L1 plus the integer refinement of n_s at the best (T, mu_buf, M),

together with the number of (T, mu_buf) combinations whose best cost lies
within 10 % of the minimum (the flatness of the frontier) and the ratio of
the two routes at every target.  The ADAPT states are admitted at every
prefix (the four-operator prefix of the (4,4) run included), the single-layer
ansatz with the filter as before.  Figure 25 draws the cost frontier (best
cost found against the target infidelity) of every route at k = 8.
Usage: PYTHONPATH=src python3 scripts/frontier_refine.py [outdir]
"""
import sys as _sys, os, time, json
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pairinglib as pl
import pairinglib.seniority as sz

OUT = _sys.argv[1] if len(_sys.argv) > 1 and not _sys.argv[1].startswith("-") else None   # None inside a notebook kernel
if OUT: matplotlib.use("Agg")
t0 = time.time()
c2f = lambda s: np.cos(np.pi * s / 2) ** 2
s2f = lambda s: np.sin(np.pi * s / 2) ** 2
M_MAX = 6
RATIO = 1.15           # geometric spacing of the step-count list (level L1)

# The coarse-grid minima of Table tab:frontier (2026-10-20), per route and target: level L0.
L0 = {"k4g1": {"HF ramp": (9.47e3, 1.89e4, 2.55e4), "refinement": (9.92e3, 1.36e4, 1.98e4)},
      "k8g1": {"HF ramp": (1.09e5, 1.63e5, 2.18e5), "refinement": (1.29e5, 1.94e5, 2.58e5)},
      "k8g2": {"HF ramp": (2.18e5, 4.00e5, 5.83e5), "refinement": (2.58e5, 3.87e5, 5.70e5)}}


class Route:
    """One adiabatic route: a path (groups, coefs, start state) per control value."""
    def __init__(self, name, step_cnots, GS):
        self.name, self.step, self.GS = name, step_cnots, GS
        self.paths = {}                         # control -> (Phi0, groups, coefs)
        self.cache = {}                         # (control, n_s) -> psi

    def state(self, ctrl, T, ns):
        key = (ctrl, float(T), int(ns))
        if key not in self.cache:
            Phi0, groups, coefs = self.paths[ctrl]
            self.cache[key] = pl.trotter_path(Phi0, groups, coefs, float(T), int(ns), order=2)
        return self.cache[key]

    def inf0(self, ctrl, T, ns):
        return 1 - abs(np.vdot(self.GS, self.state(ctrl, T, ns))) ** 2

    def cost(self, ns, M, acc, ccyc_mean, reach):
        """CNOTs per accepted preparation, Eq. (cost) with early termination."""
        cprep = 2 * int(ns) * self.step
        if M == 0:
            return float(cprep)
        return (cprep + float(np.sum(reach[:M]) * ccyc_mean)) / acc[M]


def min_steps(route, ctrl, T, target, M, K, ns_lo, ns_hi):
    """Smallest integer n_s in [ns_lo, ns_hi] (bracketing by RATIO, then bisection)
    whose infidelity after M cycles is <= target.  Returns (n_s, inf, acc, reach)
    or None if ns_hi does not reach the target."""
    def ev(ns):
        if M == 0:
            return route.inf0(ctrl, T, ns), None, None
        tr = K.track(route.state(ctrl, T, ns), route.GS, M)
        return tr["inf"][M], tr["acc"], tr["acc"][:-1]
    lo, hi = None, None
    ns = ns_lo
    while ns <= ns_hi:
        i, a, r = ev(ns)
        if i <= target:
            hi = (ns, i, a, r); break
        lo = ns; ns = max(ns + 1, int(np.ceil(ns * RATIO)))
    if hi is None:
        return None
    if lo is None:
        return hi
    while hi[0] - lo > 1:
        mid = (lo + hi[0]) // 2
        i, a, r = ev(mid)
        if i <= target: hi = (mid, i, a, r)
        else: lo = mid
    return hi


def search_case(label, k, P, g, SIG, DT, targets, ramp_T, ref_mu, ref_T, ns_lo, ns_hi, adapt_min=2e-3):
    """Levels L1 and L2 for both adiabatic routes, plus the ADAPT prefixes and the
    single-layer ansatz with the exact channel; returns everything needed for the
    tables and the figure."""
    H = sz.H_pair(k, g, P); w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]; D, hops = sz.pair_terms(k, g, P)
    cn = pl.cnot_counts(k); STEP, CSTEP = cn["step_uncontrolled"], cn["step_controlled"]
    cterm = sz.coarse_term_cnots(k // 2, k, P)
    K = pl.AveragedFilterMap(sz.fixed_step_evolution([D] + hops, DT), E0, SIG, DT)
    ccyc = 2 * CSTEP * K.mean_steps
    print(f"\n{label}: dt={DT}, sigma={SIG}, step {STEP} CNOTs (+{cterm} coarse term), cycle {ccyc:.0f} CNOTs on average  [{time.time()-t0:.0f} s]")
    ramp = Route("HF ramp", STEP, GS)
    gR, cR = sz.ramp_terms(k, g, P); hf = sz.hf_state(k, P).astype(complex)
    for T in ramp_T:
        ramp.paths[(None, float(T))] = (hf, gR, cR)
    ref = Route("refinement", STEP + cterm, GS)
    for mb in ref_mu:
        r = sz.refine_pair(k // 2, k, P, g, mu_buf=mb, M=16)
        groups = [r["Hle"], D - r["mu"] * np.eye(len(D))] + hops; coefs = [c2f, s2f] + [s2f] * len(hops)
        for T in ref_T:
            ref.paths[(float(mb), float(T))] = (r["Phi0"].astype(complex), groups, coefs)

    # ---- level L1: geometric n_s list, every control, every M (exact channel)
    ns_list = [ns_lo]
    while ns_list[-1] < ns_hi:
        ns_list.append(max(ns_list[-1] + 1, int(np.ceil(ns_list[-1] * RATIO))))
    cloud = {}                                   # route -> list of (ctrl, T, ns, M, inf, cost)
    for route in (ramp, ref):
        pts = []
        for (mb, T) in route.paths:
            for ns in ns_list:
                tr = K.track(route.state((mb, T), T, ns), GS, M_MAX)
                reach = tr["acc"][:-1]
                for M in range(M_MAX + 1):
                    pts.append(((mb, T), T, ns, M, tr["inf"][M], route.cost(ns, M, tr["acc"], ccyc, reach)))
        cloud[route.name] = pts
        print(f"  {route.name}: {len(route.paths)} controls x {len(ns_list)} step counts x {M_MAX+1} cycle numbers evaluated  [{time.time()-t0:.0f} s]")

    # ---- levels L1 -> L2 at the targets
    result = {}
    for route in (ramp, ref):
        pts = cloud[route.name]
        result[route.name] = {}
        for tgt in targets:
            ok = [p for p in pts if p[4] <= tgt]
            if not ok:
                result[route.name][tgt] = None; continue
            best1 = min(ok, key=lambda p: p[5])
            # the best (ctrl, M) pairs within 10 % of L1, each refined to the integer n_s
            per_ctrl = {}
            for p in ok:
                key = (p[0], p[3])
                if key not in per_ctrl or p[5] < per_ctrl[key][5]: per_ctrl[key] = p
            cands = sorted(per_ctrl.values(), key=lambda p: p[5])
            refined = []
            for p in cands:
                if p[5] > 1.4 * best1[5]: break
                ctrl, T, ns, M = p[0], p[1], p[2], p[3]
                lo = max(ns_lo, int(ns / RATIO))
                hit = min_steps(route, ctrl, T, tgt, M, K, lo, ns)
                if hit is not None:
                    n2, i2, a2, r2 = hit
                    refined.append((ctrl, T, n2, M, i2, route.cost(n2, M, a2, ccyc, r2)))
            best2 = min(refined, key=lambda p: p[5])
            # flatness: distinct controls (T or (mu_buf, T)) whose refined cost lies within 10 % of the L2 minimum
            per_c = {}
            for p in refined:
                if p[0] not in per_c or p[5] < per_c[p[0]]: per_c[p[0]] = p[5]
            flat = [c for c, v in per_c.items() if v <= 1.10 * best2[5]]
            result[route.name][tgt] = dict(L1=best1, L2=best2, flat=flat, n_refined=len(refined))

    # ---- the other routes with the exact channel: ADAPT prefixes (all), single-layer ansatz + filter
    others = {}
    ad = sz.adapt_vqe(k, P, g, pool="generalized", grad_tol=1e-6, max_ops=40, H=H, target=GS)
    adapt_pts = []
    for h_ in ad["history"]:
        if h_["n_ops"] > 0 and h_["infidelity"] < adapt_min:
            tr = K.track(h_["psi"].astype(complex), GS, M_MAX); reach = tr["acc"][:-1]
            for M in range(M_MAX + 1):
                c = h_["cnots"] if M == 0 else (h_["cnots"] + float(np.sum(reach[:M]) * ccyc)) / tr["acc"][M]
                adapt_pts.append((h_["n_ops"], M, tr["inf"][M], c))
    others["ADAPT"] = adapt_pts
    a1 = sz.pair_double_ansatz(k, P, g, H=H)
    tr = K.track(a1["psi"].astype(complex), GS, M_MAX); reach = tr["acc"][:-1]
    others["single-layer ansatz"] = [(a1["cnots"], M, tr["inf"][M], a1["cnots"] if M == 0 else (a1["cnots"] + float(np.sum(reach[:M]) * ccyc)) / tr["acc"][M]) for M in range(M_MAX + 1)]
    print(f"  ADAPT prefixes with 1-F < {adapt_min:.0e}: {sorted(set(p[0] for p in adapt_pts))}; single-layer ansatz {a1['cnots']} CNOTs, 1-F = {tr['inf'][0]:.2e}  [{time.time()-t0:.0f} s]")
    return dict(label=label, k=k, g=g, targets=targets, cloud=cloud, result=result, others=others, ccyc=ccyc, ns_list=ns_list,
                ramp_T=list(ramp_T), ref_mu=list(ref_mu), ref_T=list(ref_T))


def fmt_pt(p, route):
    ctrl, T, ns, M = p[0], p[1], p[2], p[3]
    return (f"({T:g}, {ns}, {M})" if route == "HF ramp" else f"({ctrl[0]:g}, {T:g}, {ns}, {M})")


def report(case, key):
    label, targets, res = case["label"], case["targets"], case["result"]
    print(f"\n{label}: best CNOTs per accepted preparation found at three search levels")
    print(f"  {'route':12s} {'target':>7s} | {'L0 coarse grid':>14s} | {'L1 fine grid':>12s} {'(T[,mu_buf],n_s,M)':>24s} | {'L2 +integer n_s':>15s} {'(T[,mu_buf],n_s,M)':>24s} | 1-F | controls within 10%")
    for route in ("HF ramp", "refinement"):
        for j, tgt in enumerate(targets):
            r = res[route][tgt]
            if r is None:
                print(f"  {route:12s} {tgt:7.0e} | {L0[key][route][j]:14.3e} | not reached"); continue
            print(f"  {route:12s} {tgt:7.0e} | {L0[key][route][j]:14.3e} | {r['L1'][5]:12.4e} {fmt_pt(r['L1'], route):>24s} | {r['L2'][5]:15.4e} {fmt_pt(r['L2'], route):>24s} | {r['L2'][4]:.2e} | {len(r['flat'])}"
                  + (f" (T = {min(c[1] for c in r['flat']):g}..{max(c[1] for c in r['flat']):g}" + (f", mu_buf = {sorted(set(c[0] for c in r['flat']))})" if route == "refinement" else ")")))
    print("  ratio refinement / ramp at the targets (L2): " + ", ".join(
        f"{tgt:.0e}: {res['refinement'][tgt]['L2'][5]/res['HF ramp'][tgt]['L2'][5]:.3f}" if res["refinement"][tgt] and res["HF ramp"][tgt] else f"{tgt:.0e}: -" for tgt in targets))
    print("  change L0 -> L2: " + ", ".join(f"{route} {tgt:.0e}: {100*(res[route][tgt]['L2'][5]/L0[key][route][j]-1):+.0f}%" for route in ("HF ramp", "refinement") for j, tgt in enumerate(targets) if res[route][tgt]))
    for name, pts in case["others"].items():
        print(f"  {name}: " + ", ".join(
            (lambda ok: f"{tgt:.0e}: {min(ok, key=lambda p: p[3])[3]:.3e} ({'ops' if name=='ADAPT' else 'CNOTs'} {min(ok, key=lambda p: p[3])[0]}, M={min(ok, key=lambda p: p[3])[1]})" if ok else f"{tgt:.0e}: -")
            ([p for p in pts if p[2] <= tgt]) for tgt in targets))


def frontier_curve(pts, tgts, cost_idx, inf_idx, M_sel=None):
    out = []
    for tgt in tgts:
        ok = [p for p in pts if p[inf_idx] <= tgt and (M_sel is None or p[3 if cost_idx == 5 else 1] in M_sel)]
        out.append(min(p[cost_idx] for p in ok) if ok else np.nan)
    return np.array(out)


# ------------------------------------------------------------------ the three cases
cases = {}
cases["k4g1"] = search_case("k=4, N=4, g=1", 4, 2, 1.0, 4.0, 0.25, (2e-3, 5e-4, 2e-4),
                            ramp_T=np.arange(3.0, 12.01, 0.25), ref_mu=(0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 8.0),
                            ref_T=np.arange(3.0, 16.01, 0.5), ns_lo=6, ns_hi=160)
report(cases["k4g1"], "k4g1")
cases["k8g1"] = search_case("k=8, N=8, g=1", 8, 4, 1.0, 4.0, 0.1, (1e-3, 3e-4, 1e-4),
                            ramp_T=np.arange(3.5, 10.01, 0.25), ref_mu=(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 12.0, 16.0),
                            ref_T=np.arange(3.5, 12.01, 0.5), ns_lo=16, ns_hi=400)
report(cases["k8g1"], "k8g1")
cases["k8g2"] = search_case("k=8, N=8, g=2", 8, 4, 2.0, 1.5, 0.05, (1e-3, 3e-4, 1e-4),
                            ramp_T=np.arange(5.0, 16.01, 0.5), ref_mu=(4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 24.0, 32.0),
                            ref_T=np.arange(5.0, 16.01, 0.5), ns_lo=16, ns_hi=400)
report(cases["k8g2"], "k8g2")

# ------------------------------------------------------------------ a finer T grid (steps of 0.05) around the k=8, g=1 ramp minima
print("\nk=8, g=1, HF ramp on a T grid of step 0.05 around the L2 minima (integer n_s by bracketing + bisection, M=0):")
_H = sz.H_pair(8, 1.0, 4); _w, _V = np.linalg.eigh(_H); _GS = _V[:, 0]
_gR, _cR = sz.ramp_terms(8, 1.0, 4); _hf = sz.hf_state(8, 4).astype(complex)
_rampf = Route("HF ramp", pl.cnot_counts(8)["step_uncontrolled"], _GS)
for tgt, T0, T1 in ((1e-3, 4.5, 5.1), (3e-4, 5.0, 5.5)):
    best = None
    for T in np.arange(T0, T1 + 1e-9, 0.05):
        _rampf.paths[(None, float(T))] = (_hf, _gR, _cR)
        hit = min_steps(_rampf, (None, float(T)), T, tgt, 0, None, 16, 120)
        if hit is not None and (best is None or 2 * hit[0] * _rampf.step < best[0]):
            best = (2 * hit[0] * _rampf.step, float(T), hit[0], hit[1])
    print(f"  target {tgt:.0e}: best {best[0]:.4e} CNOTs at T={best[1]:.2f}, n_s={best[2]} (1-F = {best[3]:.2e}); L2 of the 0.25 grid: {cases['k8g1']['result']['HF ramp'][tgt]['L2'][5]:.4e}")

# ------------------------------------------------------------------ figure 25: the cost frontier at k=8
tg = np.logspace(np.log10(3e-3), np.log10(4e-5), 25)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for a_, key, ttl in ((ax[0], "k8g1", "$k=8$, $N=8$, $g=1$ ($\\delta t=0.1$, $\\sigma=4$)"), (ax[1], "k8g2", "$g=2$ ($\\delta t=0.05$, $\\sigma=1.5$)")):
    case = cases[key]
    for route, col in (("HF ramp", "k"), ("refinement", "C4")):
        pts = case["cloud"][route]
        a_.loglog(frontier_curve(pts, tg, 5, 4, M_sel={0}), tg, "--", color=col, lw=1, label=f"{route}, unfiltered ($M=0$)")
        a_.loglog(frontier_curve(pts, tg, 5, 4), tg, "-", color=col, lw=1.8, label=f"{route}, $M=0\\ldots6$ (fine grid)")
        for j, tgt in enumerate(case["targets"]):
            r = case["result"][route][tgt]
            if r: a_.loglog([r["L2"][5]], [tgt], "o", color=col, ms=6, mfc="none")
            a_.loglog([L0[key][route][j]], [tgt], "s", color=col, ms=5, alpha=0.5)
    a_.loglog(frontier_curve(case["others"]["ADAPT"], tg, 3, 2, M_sel={0}), tg, "x-", color="C3", ms=4, label="ADAPT prefixes, unfiltered")
    a_.loglog(frontier_curve(case["others"]["ADAPT"], tg, 3, 2), tg, ":", color="C3", label="ADAPT prefixes + filter")
    a_.loglog(frontier_curve(case["others"]["single-layer ansatz"], tg, 3, 2), tg, "s-", color="C1", ms=4, label="single-layer ansatz + filter")
    for tgt in case["targets"]: a_.axhline(tgt, color="gray", ls=":", lw=0.8)
    a_.set_xlabel("best CNOTs per accepted preparation found"); a_.set_ylabel("target infidelity $1-\\mathcal{F}_M$"); a_.set_title(ttl)
    a_.grid(alpha=0.3, which="both"); a_.legend(fontsize=7, loc="lower left" if key == "k8g1" else "upper right")
    a_.set_ylim(3e-5, 4e-3)
ax[0].text(0.02, 0.98, "open circles: level L2 at the targets\nfilled squares: coarse grid of 2026-10-20 (L0)", transform=ax[0].transAxes, va="top", fontsize=7)
plt.tight_layout()
if OUT:
    plt.savefig(os.path.join(OUT, "fig25.png"), dpi=150); print("wrote", OUT)
    with open(os.path.join(OUT, "frontier_refine.json"), "w") as f:
        json.dump({key: {route: {f"{tgt:.0e}": dict(L1=[str(r["L1"][0]), r["L1"][1], r["L1"][2], r["L1"][3], r["L1"][4], r["L1"][5]],
                                                   L2=[str(r["L2"][0]), r["L2"][1], r["L2"][2], r["L2"][3], r["L2"][4], r["L2"][5]], flat=[str(c) for c in r["flat"]])
                                 for tgt, r in res.items() if r} for route, res in case["result"].items()} for key, case in cases.items()}, f, indent=1)
else:
    plt.show()
print(f"[{time.time()-t0:.0f} s]")
