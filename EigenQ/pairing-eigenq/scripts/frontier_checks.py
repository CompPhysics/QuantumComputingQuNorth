#!/usr/bin/env python3
"""Pre-submission revision of 2026-10-20 (paper/pairing_pre_submission_critical_assessment.md):
the digital filter treated with its exact averaged post-selection map, and a
JOINT operating-point search for the adiabatic routes (notebook section
EXT:FRONTIER, figure 24).

  (a) The averaged map K(rho) = E_t[K_t rho K_t^+] of the digital filter is
      integrated with the boundary-split quadrature of `pl.half_gaussian_nodes`
      (the product formula changes at every |t| = n dt, so a global
      Gauss-Hermite rule is not converged; its 60/80-node values are listed for
      the record).  Dominant eigenvalue, floor, residual and the change under
      quadrature refinement, at k=4 (dt = 0.25, 0.5, 1) and k=8 (dt = 0.25,
      0.1, 0.05 at g=1; 0.1, 0.05 at g=2).
  (b) Because the cycle times are independent, the ensemble of accepted runs is
      rho_M = K^M(rho_in)/Tr, the acceptance Tr K^M(rho_in), and the expected
      cost per accepted preparation of Eq. (cost) is deterministic.  Every input
      of the matched comparisons (Tables II and V of paper I, the ADAPT states)
      is re-evaluated with the exact channel, and the adiabatic routes (HF ramp
      and shifted refinement) are searched JOINTLY over the time T, the shift
      mu_buf, the second-order step count n_s and the cycle number M = 0..6, the
      unfiltered (M = 0) points included, at both sizes and both couplings.
      Output: the minimum CNOTs per accepted preparation at the targets, per
      route (table tab:frontier of paper I), the early-termination size of
      Eq. (cost), and figure 24 (k=8, g=1 and g=2, every input incl. ADAPT).
Usage: PYTHONPATH=src python3 scripts/frontier_checks.py [outdir]
"""
import sys as _sys, os, time
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


def hermite_floor(Uf, E0, GS, sigma, n):
    """The global Gauss-Hermite rule used before this revision (for the record)."""
    from numpy.polynomial.hermite_e import hermegauss
    x, wq = hermegauss(n); wq = wq / np.sqrt(2 * np.pi)
    K = pl.AveragedFilterMap.__new__(pl.AveragedFilterMap)
    K.ts, K.ws = np.abs(sigma * x), wq
    K.Ks = [0.5 * (np.eye(len(GS)) + np.exp(1j * E0 * t) * Uf(t)) for t in K.ts]
    d = K.dominant(GS, iters=4000)
    return d["eigenvalue"], d["floor"]


# ------------------------------------------------------------------ (a) converged floors
print("(a) averaged post-selection map of the digital filter (second order, E_obj = E0): boundary-split quadrature (m = 6 and 12 Gauss-Legendre points per step interval, 8 sigma) versus the global Gauss-Hermite rule")
floors = {}
for (k, P, g, sigma, dts) in ((4, 2, 1.0, 4.0, (0.25, 0.5, 1.0)), (8, 4, 1.0, 4.0, (0.25, 0.1, 0.05)), (8, 4, 2.0, 1.5, (0.1, 0.05))):
    H = sz.H_pair(k, g, P); w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]; D, hops = sz.pair_terms(k, g, P)
    for dt in dts:
        Uf = sz.fixed_step_evolution([D] + hops, dt)
        out = pl.averaged_map_floor(Uf, E0, GS, sigma, dt)
        he = hermite_floor(Uf, E0, GS, sigma, 60 if k == 8 else 80) if dt >= 0.1 else (np.nan, np.nan)
        floors[(k, g, dt)] = out
        print(f"  k={k}, g={g:.0f}, sigma={sigma}, dt={dt:5.2f}: eigenvalue {out['eigenvalue']:.7f} (acceptance -{100*(1-out['eigenvalue']):.2f}%/cycle), floor {out['floor']:.4e}; "
              f"m=12: {out['floor_check']:.4e}, eigenvalue {out['eigenvalue_check']:.7f}; residual {out['residual']:.1e}, {out['iters']} iterations; mean steps/cycle {out['mean_steps']:.2f}"
              + (f"  [Gauss-Hermite {60 if k == 8 else 80} nodes: eigenvalue {he[0]:.4f}, floor {he[1]:.2e}]" if dt >= 0.1 else ""))
print(f"[{time.time()-t0:.0f} s]")


# ------------------------------------------------------------------ (b) exact-channel matched comparison and joint search
def evaluate(inputs, K, GS, cstep_ctrl):
    """inputs: name -> (psi, C_prep, route).  Exact channel: infidelity, acceptance,
    expected CNOTs per accepted preparation (early termination) and per full attempt."""
    rows = {}
    ccyc = 2 * cstep_ctrl * K.mean_steps
    for name, (psi, cprep, route) in inputs.items():
        tr = K.track(psi, GS, M_MAX)
        reach = tr["acc"][:-1]                                       # P_{m-1} (acc[0] = 1): probability that cycle m is executed
        num = cprep + np.cumsum(reach * ccyc)                        # early-terminated expected cost of an attempt
        cacc = np.concatenate([[float(cprep)], num / tr["acc"][1:]])
        catt = cprep + ccyc * np.arange(M_MAX + 1)                   # one attempt running all M cycles
        rows[name] = dict(p=1 - tr["inf"][0], cprep=cprep, inf=tr["inf"], acc=tr["acc"], cacc=cacc, catt=catt, route=route,
                          early=(catt[-1] - num[-1]) / catt[-1])
    return rows, ccyc


def minima(rows, targets):
    """Per route, the cheapest (input, M) reaching each target."""
    out = {}
    for tgt in targets:
        best = {}
        for name, d in rows.items():
            ok = np.where(d["inf"] <= tgt)[0]
            if ok.size:
                m = ok[np.argmin(d["cacc"][ok])]
                if d["route"] not in best or d["cacc"][m] < best[d["route"]][0]:
                    best[d["route"]] = (d["cacc"][m], name, int(m), d["inf"][m])
        out[tgt] = best
    return out


def report(label, rows, ccyc, targets, routes):
    print(f"\n{label}: exact averaged channel ({ccyc:.0f} CNOTs per cycle on average)")
    print(f"{'input':58s} {'p':>8s} {'C_prep':>8s} | 1-F_M, M=0..6 | P_6 | CNOTs per accepted preparation, M=0..6 | early-term.")
    for name, d in rows.items():
        print(f"{name:58s} {d['p']:8.5f} {d['cprep']:8d} | " + " ".join(f"{x:.2e}" for x in d["inf"]) + f" | {d['acc'][-1]:.3f} | " + " ".join(f"{c:.2e}" for c in d["cacc"]) + f" | {100*d['early']:.1f}%")
    mins = minima(rows, targets)
    print("  minimum CNOTs per accepted preparation (input, cycles) per route:")
    for route in routes:
        print(f"    {route:34s} " + " | ".join((f"{tgt:.0e}: {mins[tgt][route][0]:.2e} [{mins[tgt][route][1]}, M={mins[tgt][route][2]}, {mins[tgt][route][3]:.1e}]" if route in mins[tgt] else f"{tgt:.0e}: -") for tgt in targets))
    return mins


def adiabatic_grid(k, P, g, GS, D, hops, ramp_T, ramp_ns, ref_mu, ref_T, ref_ns, step, cterm):
    """Trotterized HF ramps and shifted refinements k/2 -> k over a grid of (T, n_s) and (mu_buf, T, n_s)."""
    inputs = {}
    gR, cR = sz.ramp_terms(k, g, P); hf = sz.hf_state(k, P).astype(complex)
    for T in ramp_T:
        for ns in ramp_ns:
            psi = pl.trotter_path(hf, gR, cR, float(T), ns, order=2)
            inputs[f"HF ramp, T={T}, n_s={ns}"] = (psi, 2 * ns * step, "HF ramp (grid)")
    kl = k // 2
    for mb in ref_mu:
        r = sz.refine_pair(kl, k, P, g, mu_buf=mb, M=160)
        groups = [r["Hle"], D - r["mu"] * np.eye(len(D))] + hops; coefs = [c2f, s2f] + [s2f] * len(hops)
        for T in ref_T:
            for ns in ref_ns:
                psi = pl.trotter_path(r["Phi0"].astype(complex), groups, coefs, float(T), ns, order=2)
                inputs[f"refined {kl}->{k}, mu_buf={mb}, T={T}, n_s={ns}"] = (psi, 2 * ns * (step + cterm), f"refinement {kl}->{k} (grid)")
    return inputs


all_rows = {}
# ---- k = 4, N = 4, g = 1 (the inputs of Table II plus the grid)
k, P, g, SIG, DT = 4, 2, 1.0, 4.0, 0.25
H = sz.H_pair(k, g, P); w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]; D, hops = sz.pair_terms(k, g, P)
cn = pl.cnot_counts(k); STEP, CSTEP = cn["step_uncontrolled"], cn["step_controlled"]
K4 = pl.AveragedFilterMap(sz.fixed_step_evolution([D] + hops, DT), E0, SIG, DT)
a4 = sz.pair_double_ansatz(4, P, g); ad4 = sz.adapt_vqe(4, P, g, pool="generalized", grad_tol=1e-6, max_ops=40, H=H, target=GS)
r3 = sz.refine_pair(3, 4, P, g, mu_buf=0.6, M=160)                           # prolonged exact k=3 state (two pair rotations, 96 CNOTs)
D3, hops3 = sz.pair_terms(3, g, P); Pm = r3["Pm"]
g34 = [Pm @ (D3 - r3["mu"] * np.eye(3)) @ Pm.T] + [Pm @ h @ Pm.T for h in hops3] + [D - r3["mu"] * np.eye(6)] + hops; c34 = [c2f] * 4 + [s2f] * 7
inputs = {"prolonged determinant (k=2)": (sz.hf_state(4, P).astype(complex), 0, "determinant + filter"),
          "four-pair ansatz (192)": (a4["psi"].astype(complex), 192, "single-layer ansatz"),
          "prolonged k=3 state (96)": (r3["Phi0"].astype(complex), 96, "prolonged coarse state + filter"),
          "ADAPT, 5 operators (240)": (ad4["history"][5]["psi"].astype(complex), 240, "ADAPT")}
for T, ns in ((8, 22), (8, 43), (8, 86)):
    psi = pl.trotter_path(r3["Phi0"].astype(complex), g34, c34, float(T), ns, order=2)
    inputs[f"refined 3->4, T={T}, n_s={ns} (+96)"] = (psi, 96 + 2 * ns * (STEP + sz.coarse_term_cnots(3, 4, P)), "refinement 3->4 (grid)")
inputs.update(adiabatic_grid(k, P, g, GS, D, hops, ramp_T=(4, 6, 8, 12), ramp_ns=(8, 11, 16, 22, 32, 43, 64, 86),
                             ref_mu=(0.6, 1.5, 3.0), ref_T=(6, 8, 12, 16, 30), ref_ns=(11, 16, 22, 32, 43, 64, 86, 128), step=STEP, cterm=sz.coarse_term_cnots(2, 4, P)))
rows4, ccyc4 = evaluate(inputs, K4, GS, CSTEP)
ROUTES = ["determinant + filter", "single-layer ansatz", "prolonged coarse state + filter", "ADAPT", "HF ramp (grid)", "refinement 2->4 (grid)", "refinement 3->4 (grid)"]
mins4 = report("k=4, N=4, g=1 (dt=0.25, sigma=4)", rows4, ccyc4, (2e-3, 5e-4, 2e-4), ROUTES)
all_rows["k4g1"] = (rows4, mins4, [])
print(f"[{time.time()-t0:.0f} s]")

# ---- k = 8, N = 8, g = 1 and 2
k, P = 8, 4
cn = pl.cnot_counts(k); STEP, CSTEP = cn["step_uncontrolled"], cn["step_controlled"]
A48, A68 = sz.coarse_term_cnots(4, 8, P), sz.coarse_term_cnots(6, 8, P)
ROUTES8 = ["determinant + filter", "single-layer ansatz", "prolonged coarse state + filter", "ADAPT", "HF ramp (grid)", "refinement 4->8 (grid)", "refinement 6->8 (coarse ansatz)"]
for g, SIG, DT, ref_mu, ref_T, ramp_T in ((1.0, 4.0, 0.1, (2.0, 4.0, 8.0), (6, 8, 12), (4, 6, 8)), (2.0, 1.5, 0.05, (4.0, 8.0, 16.0), (8, 12, 16), (8, 12, 16))):
    H = sz.H_pair(k, g, P); w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]; D, hops = sz.pair_terms(k, g, P)
    K8 = pl.AveragedFilterMap(sz.fixed_step_evolution([D] + hops, DT), E0, SIG, DT)
    a8 = sz.pair_double_ansatz(k, P, g, H=H); a6 = sz.pair_double_ansatz(6, P, g)
    ad = sz.adapt_vqe(k, P, g, pool="generalized", grad_tol=1e-6, max_ops=40, H=H, target=GS)
    r6 = sz.refine_pair(6, k, P, g, mu_buf=2.0, M=160, psi_low=a6["psi"]); Pm = r6["Pm"]; Dc, hopsc = sz.pair_terms(6, g, P)
    g68 = [Pm @ (Dc - r6["mu"] * np.eye(15)) @ Pm.T] + [Pm @ h @ Pm.T for h in hopsc] + [D - r6["mu"] * np.eye(70)] + hops; c68 = [c2f] * (1 + len(hopsc)) + [s2f] * (1 + len(hops))
    inputs = {"HF determinant": (sz.hf_state(k, P).astype(complex), 0, "determinant + filter"),
              "16-pair ansatz (768)": (a8["psi"].astype(complex), 768, "single-layer ansatz"),
              "prolonged k=6 ansatz (384)": (r6["Phi0"].astype(complex), 384, "prolonged coarse state + filter")}
    adapt_sel = []
    for tgt in (1e-3, 3e-4, 1e-4):
        hit = next((h_ for h_ in ad["history"] if h_["infidelity"] <= tgt), None)
        if hit and hit["n_ops"] not in adapt_sel: adapt_sel.append(hit["n_ops"])
    for h_ in ad["history"]:
        if h_["n_ops"] > 0 and h_["infidelity"] < 2e-3:
            inputs[f"ADAPT, {h_['n_ops']} operators ({h_['cnots']})"] = (h_["psi"].astype(complex), h_["cnots"], "ADAPT")
    if g == 1.0:
        for ns in (80, 160):
            psi = pl.trotter_path(r6["Phi0"].astype(complex), g68, c68, 8.0, ns, order=2)
            inputs[f"refined 6->8 (coarse ansatz), mu_buf=2, T=8, n_s={ns} (+384)"] = (psi, 384 + 2 * ns * (STEP + A68), "refinement 6->8 (coarse ansatz)")
    inputs.update(adiabatic_grid(k, P, g, GS, D, hops, ramp_T=ramp_T, ramp_ns=(20, 40, 60, 80, 120, 160, 240, 320),
                                 ref_mu=ref_mu, ref_T=ref_T, ref_ns=(40, 60, 80, 120, 160, 240, 320), step=STEP, cterm=A48))
    rows, ccyc = evaluate(inputs, K8, GS, CSTEP)
    mins = report(f"k=8, N=8, g={g:.0f} (dt={DT}, sigma={SIG})", rows, ccyc, (1e-3, 3e-4, 1e-4), ROUTES8)
    all_rows[f"k8g{g:.0f}"] = (rows, mins, [f"ADAPT, {n} operators ({48*n})" for n in adapt_sel])
    print(f"[{time.time()-t0:.0f} s]")

# ---- figure 24: k=8, g=1 and g=2, every input with the exact channel (grids: the frontier only)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for a_, key, ttl, lloc in ((ax[0], "k8g1", "$k=8$, $N=8$, $g=1$ ($\\delta t=0.1$, $\\sigma=4$)", "upper center"), (ax[1], "k8g2", "$g=2$ ($\\delta t=0.05$, $\\sigma=1.5$)", "lower left")):
    rows, mins, adapt_sel = all_rows[key]
    for route, mk, col in (("determinant + filter", "o-", "C0"), ("single-layer ansatz", "s-", "C1"), ("prolonged coarse state + filter", "d-", "C2"), ("refinement 6->8 (coarse ansatz)", "P-", "C6")):
        for j, (name, d) in enumerate((n, d) for n, d in rows.items() if d["route"] == route):
            c = d["cacc"][1:] if d["cprep"] == 0 else d["cacc"]; i = d["inf"][1:] if d["cprep"] == 0 else d["inf"]
            a_.loglog(c, np.maximum(i, 1e-7), mk, color=col, ms=4, alpha=1.0 if j == 0 else 0.5, label=(route + (" ($M=0\\ldots6$)" if route != "refinement 6->8 (coarse ansatz)" else "") if j == 0 else None))
    ad_rows = [(n, d) for n, d in rows.items() if d["route"] == "ADAPT"]
    a_.loglog([d["cprep"] for n, d in ad_rows], [max(d["inf"][0], 1e-7) for n, d in ad_rows], "x-", color="C3", ms=5, label="ADAPT prefixes, unfiltered ($M=0$)")
    for j, (n, d) in enumerate((n, d) for n, d in ad_rows if n in adapt_sel):
        a_.loglog(d["cacc"], np.maximum(d["inf"], 1e-7), "x:", color="C3", ms=4, alpha=0.7, label=("ADAPT + filter at the target crossings" if j == 0 else None))
    for route, mk, col, lab in (("HF ramp (grid)", "^", "k", "HF ramp: $(T,n_s,M)$ grid and its frontier"), ("refinement 4->8 (grid)", "v", "C4", "refinement $4\\to8$: $(T,\\mu_{\\rm buf},n_s,M)$ grid and its frontier")):
        pts = np.array([(d["cacc"][m], d["inf"][m], m) for d in rows.values() if d["route"] == route for m in range(M_MAX + 1)])
        a_.loglog(pts[:, 0], np.maximum(pts[:, 1], 1e-7), mk, color=col, ms=3, alpha=0.2)
        order = np.argsort(pts[:, 0]); best = np.inf; fx, fy = [], []
        for i_ in order:
            if pts[i_, 1] < best: best = pts[i_, 1]; fx.append(pts[i_, 0]); fy.append(pts[i_, 1])
        a_.loglog(fx, np.maximum(fy, 1e-7), mk + "-", color=col, ms=5, label=lab)
    for tgt in (1e-3, 3e-4, 1e-4): a_.axhline(tgt, color="gray", ls=":", lw=0.8)
    a_.set_xlabel("CNOTs per accepted preparation (input + filter)"); a_.set_ylabel("exact ensemble infidelity $1-\\mathcal{F}_M$"); a_.set_title(ttl); a_.legend(fontsize=6.5, loc=lloc); a_.grid(alpha=0.3, which="both")
    a_.set_ylim(1e-6, 1)
plt.tight_layout()
if OUT: plt.savefig(os.path.join(OUT, "fig24.png"), dpi=150); print("wrote", OUT)
else: plt.show()
print(f"[{time.time()-t0:.0f} s]")
