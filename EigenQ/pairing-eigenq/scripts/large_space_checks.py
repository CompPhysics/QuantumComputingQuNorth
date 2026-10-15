#!/usr/bin/env python3
"""A larger model space: N = 8 particles (four pairs) on k = 8 doubly degenerate
levels, sixteen qubits, seniority zero (notebook section EXT:LARGE, figs 21-22).

Everything is computed in the pair basis of `pairinglib.seniority` (the exact
seniority-zero block, 70 states for k = 8, P = 4), which is first validated
against the full-sector library at k = 4, N = 4 (block identity, every anchor
of the paper).  Then, for k = 8, P = 4:
  (a) spectrum, HF overlap and the restricted pair-double ansatz versus g;
      the closed form chi_u(0) = P(k-P) g^2 / (4 b^2) of the response norm;
  (b) exact-slice refinement chains k_low -> 8 for k_low = 4..7 and shifts
      mu_buf = 0.6, 1, 2 (T* table), the HF adiabatic ramp, and the coupling
      dependence of ramp, 4->8 and 6->8 transport (fig 21);
  (c) the circuit level: product-formula steps, the averaged post-selection
      map of the digital filter (dt_max = 0.25 is rejected at these norms;
      0.1 is used), the matched comparison of nine inputs to the same
      Trotterized filter at g = 1 with M = 0..6, CNOTs per attempt and per
      accepted preparation, bootstrap errors and minimum cost at targets, the
      same budget in the one-qubit-per-level pair encoding, and the g = 2
      comparison (fig 22).
Usage: PYTHONPATH=src python3 scripts/large_space_checks.py [outdir]
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

# ------------------------------------------------------------------ validation at k = 4, N = 4
k4, N4, G = 4, 4, 1.0
nq, st, ix = pl.build_sector(k4, N4); Hfull = pl.H_pairing_sparse(k4, G, N4, st, ix).toarray()
full, pair = sz.seniority_zero_states(k4, N4)
Hp4 = sz.H_pair(k4, G, 2); blk = np.zeros_like(Hp4); blk[np.ix_(pair, pair)] = Hfull[np.ix_(full, full)]
D4, hops4 = sz.pair_terms(k4, G, 2); tm4 = pl.model_terms(k4, G, N4)
dterm = 0.0
for i, T in enumerate([D4] + hops4):
    Bp = np.zeros_like(T); Bp[np.ix_(pair, pair)] = tm4[i][np.ix_(full, full)]; dterm = max(dterm, abs(Bp - T).max())
r4 = sz.refine_pair(2, 4, 2, G, mu_buf=0.6, M=160)
a4 = sz.pair_double_ansatz(4, 2, G)
psi25 = sz._evolve(r4["Phi0"], r4["eigs"], 25.0); c2 = np.abs(r4["V"].T @ psi25) ** 2; q = 0.5 * (1 + np.exp(-(r4["w"] - r4["w"][0]) ** 2 * 9 / 2))
c2f = lambda s: np.cos(np.pi * s / 2) ** 2; s2f = lambda s: np.sin(np.pi * s / 2) ** 2
p16 = abs(np.vdot(r4["Psih"], pl.trotter_path(r4["Phi0"], [r4["Hle"], D4 - r4["mu"] * np.eye(6)] + hops4, [c2f, s2f] + [s2f] * 6, 16.0, 43, order=2))) ** 2
print("validation of the pair basis at k=4, N=4, g=1 against the full-sector library:")
print(f"  |H_pair - seniority-zero block of H| = {abs(blk - Hp4).max():.1e}; |terms - blocks| = {dterm:.1e}; E0 = {np.linalg.eigvalsh(Hp4)[0]:.13f}")
print(f"  2->4: p0 = {r4['p0']:.13f}, p(25) = {sz.refine_pair_overlaps(r4, [25.0])[0]:.10f}, Trotterized T=16 n_s=43: {p16:.4f}; four-pair ansatz E = {a4['E']:.13f}, 1-F = {1 - abs(np.vdot(r4['Psih'], a4['psi']))**2:.4e}")
print("  exact-ensemble rodeo (T=25 input, sigma=3): " + ", ".join(f"M={M}: {np.sum(c2[1:]*q[1:]**M)/np.sum(c2*q**M):.3e}" for M in (1, 2, 6, 10)))

# ------------------------------------------------------------------ (a) k = 8, P = 4
k, P = 8, 4
print(f"\n(a) k=8, P=4 (N=8, 16 qubits), seniority-zero block of {len(sz.pair_basis(k, P)[0])} states; first-order step {pl.cnot_counts(k)['step_uncontrolled']} CNOTs "
      f"(controlled {pl.cnot_counts(k)['step_controlled']}); coarse term 4->8 {sz.coarse_term_cnots(4, 8, P)} CNOTs, 6->8 {sz.coarse_term_cnots(6, 8, P)} (Gray-code convention of paper I)")
print(f"{'g':>4s} {'E_HF':>7s} {'E_0':>9s} {'E_corr':>8s} {'gap':>6s} {'p_HF':>7s} | {'16-pair: E':>10s} {'1-F':>9s} | {'k=6 8-pair 1-F':>14s} {'p0(6->8)':>9s} | chi_u(0) b=0.6 (closed form)")
spec = {}
for g in (0.5, 1.0, 1.5, 2.0):
    H = sz.H_pair(k, g, P); w, V = np.linalg.eigh(H); hf = sz.hf_state(k, P)
    a = sz.pair_double_ansatz(k, P, g, H=H); a6 = sz.pair_double_ansatz(6, P, g); w6, V6 = np.linalg.eigh(sz.H_pair(6, g, P))
    r = sz.refine_pair(4, k, P, g, mu_buf=0.6, M=4); A, B = r["Hle"], r["Hhs"]; wa, Va = np.linalg.eigh(A); i0 = np.argmin(wa)
    chi = sum(abs(Va[:, n] @ (B - A) @ Va[:, i0]) ** 2 / (wa[n] - wa[i0]) ** 2 for n in range(len(wa)) if n != i0)
    r6 = sz.refine_pair(6, k, P, g, mu_buf=0.6, M=4)
    spec[g] = dict(E0=w[0], gap=w[1] - w[0], pHF=abs(V[:, 0] @ hf) ** 2, ans=a, inf_ans=1 - abs(np.vdot(V[:, 0], a["psi"])) ** 2, inf6=1 - abs(np.vdot(V6[:, 0], a6["psi"])) ** 2, p06=r6["p0"])
    print(f"{g:4.1f} {hf @ H @ hf:7.3f} {w[0]:9.5f} {w[0] - hf @ H @ hf:8.4f} {w[1] - w[0]:6.3f} {spec[g]['pHF']:7.4f} | {a['E']:10.5f} {spec[g]['inf_ans']:9.2e} | {spec[g]['inf6']:14.2e} {r6['p0']:9.4f} | {chi:.4f} ({P * (k - P) * g**2 / (4 * 0.36):.4f})")

# ------------------------------------------------------------------ (b) exact-slice refinement and the ramp
g = 1.0
H = sz.H_pair(k, g, P); w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]; hf = sz.hf_state(k, P)
D, hops = sz.pair_terms(k, g, P)
Ts = np.array([2, 4, 6, 8, 10, 12, 16, 20, 25, 30, 40, 50, 64])
print(f"\n(b) exact-slice refinement k_low -> 8 at g=1 (160 slices): T* = first grid value with 1-p < 1e-2 / 1e-3; grid {list(Ts)}")
chains = {}
for kl in (4, 5, 6, 7):
    for mb in (0.6, 1.0, 2.0):
        r = sz.refine_pair(kl, k, P, g, mu_buf=mb, M=160); ov = sz.refine_pair_overlaps(r, Ts)
        gmin = min(np.diff(np.linalg.eigvalsh(np.cos(np.pi * s / 2) ** 2 * r["Hle"] + np.sin(np.pi * s / 2) ** 2 * r["Hhs"])[:2])[0] for s in np.linspace(0, 1, 41))
        T2 = next((T for T, o in zip(Ts, ov) if 1 - o < 1e-2), None); T3 = next((T for T, o in zip(Ts, ov) if 1 - o < 1e-3), None)
        chains[(kl, mb)] = dict(ov=ov, p0=r["p0"], gmin=gmin, T2=T2, T3=T3)
        print(f"  {kl}->8, mu_buf={mb}: p0 = {r['p0']:.4f}, E_low = {r['E_low']:.4f}, min path gap {gmin:.3f}, T*(1e-2) = {T2}, T*(1e-3) = {T3};  p(8) = {ov[3]:.5f}, p(16) = {ov[6]:.5f}, p(30) = {ov[9]:.5f}")
print("  shift scan for the 4->8 chain (exact slices; T* for 1-p < 1e-3; second-order steps within 1e-3 of it; cost = 2 n_s (1360 + 254) CNOTs):")
shift_scan = {}
for gg in (1.0, 2.0):
    Hg = sz.H_pair(k, gg, P); wg_, Vg = np.linalg.eigh(Hg); GSg = Vg[:, 0]; Dg, hopsg = sz.pair_terms(k, gg, P)
    for mb in (0.6, 1.0, 2.0, 4.0, 8.0, 16.0):
        r = sz.refine_pair(4, k, P, gg, mu_buf=mb, M=160); ov = sz.refine_pair_overlaps(r, Ts)
        Tst = next((T for T, o in zip(Ts, ov) if 1 - o < 1e-3), None)
        groups = [r["Hle"], Dg - r["mu"] * np.eye(70)] + hopsg; coefs = [c2f, s2f] + [s2f] * len(hopsg)
        ns_, p_ = None, None
        if Tst is not None and Tst <= 25:
            pex = sz.refine_pair_overlaps(r, [float(Tst)])[0]
            for ns in (10, 20, 40, 80, 160, 320, 640):
                p_ = abs(np.vdot(GSg, pl.trotter_path(r["Phi0"], groups, coefs, float(Tst), ns, order=2))) ** 2
                if abs(p_ - pex) <= 1e-3: ns_ = ns; break
        shift_scan[(gg, mb)] = dict(ov=ov, Tst=Tst, ns=ns_, p=p_)
        print(f"    g={gg}, mu_buf={mb:4.1f}: chi_u(0) = {P*(k-P)*gg**2/(4*mb**2):7.3f}, T* = {Tst}, " + (f"n_s = {ns_} (p = {p_:.5f}), {2*ns_*(1360+254)} CNOTs" if ns_ else "(T* > 25, not Trotterized)"))
eigsR = sz._slice_eigs(D, H, 160); ovR = np.array([abs(np.vdot(GS, sz._evolve(hf, eigsR, T))) ** 2 for T in Ts])
gR = min(np.diff(np.linalg.eigvalsh((1 - u) * D + u * H)[:2])[0] for u in np.sin(np.pi * np.linspace(0, 1, 41) / 2) ** 2)
print(f"  HF ramp D -> H from the determinant (p0 = {abs(GS @ hf)**2:.4f}, min gap {gR:.3f}): " + " ".join(f"p({T})={o:.5f}" for T, o in zip(Ts, ovR) if T in (2, 4, 6, 8, 12, 16)))
# coupling dependence
Tg = np.array([2, 4, 6, 8, 10, 12, 16, 20, 25])
gdep = {}
for gg in (0.5, 1.0, 1.5, 2.0):
    Hg = sz.H_pair(k, gg, P); wg_, Vg = np.linalg.eigh(Hg); GSg = Vg[:, 0]; Dg, _ = sz.pair_terms(k, gg, P)
    eR = sz._slice_eigs(Dg, Hg, 160); ramp = [abs(np.vdot(GSg, sz._evolve(hf, eR, T))) ** 2 for T in Tg]
    rA = sz.refine_pair(4, k, P, gg, mu_buf=2.0, M=160); rB = sz.refine_pair(6, k, P, gg, mu_buf=2.0, M=160)
    gdep[gg] = dict(ramp=np.array(ramp), r48=sz.refine_pair_overlaps(rA, Tg), r68=sz.refine_pair_overlaps(rB, Tg), p06=rB["p0"])
    print(f"  g={gg}: 1-p at T=8/16/25: ramp {1-gdep[gg]['ramp'][3]:.1e}/{1-gdep[gg]['ramp'][6]:.1e}/{1-gdep[gg]['ramp'][8]:.1e}; 4->8 (mu_buf=2) {1-gdep[gg]['r48'][3]:.1e}/{1-gdep[gg]['r48'][6]:.1e}/{1-gdep[gg]['r48'][8]:.1e}; "
          f"6->8 (p0 {rB['p0']:.3f}) {1-gdep[gg]['r68'][3]:.1e}/{1-gdep[gg]['r68'][6]:.1e}/{1-gdep[gg]['r68'][8]:.1e}; 16-pair ansatz {spec[gg]['inf_ans']:.1e}")

fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
for (kl, mb), d in chains.items():
    if mb in (0.6, 2.0):
        ax[0].semilogy(Ts, np.maximum(1 - d["ov"], 1e-9), "o-" if mb == 2.0 else "s--", ms=3, color=f"C{kl-4}", label=f"${kl}\\to8$, $\\mu_{{\\rm buf}}={mb}$ ($p_0={d['p0']:.3f}$)")
ax[0].semilogy(Ts, np.maximum(1 - shift_scan[(1.0, 4.0)]["ov"], 1e-9), "-x", ms=5, color="C4", label="$4\\to8$, $\\mu_{\\rm buf}=4$")
ax[0].semilogy(Ts, np.maximum(1 - ovR, 1e-9), "k-^", ms=4, label="HF ramp $D\\to H$ (no shift)")
ax[0].set_xlabel("$T$"); ax[0].set_ylabel("$1-p(T)$"); ax[0].set_title("$k=8$, $N=8$, $g=1$: exact-slice refinement and the ramp"); ax[0].legend(fontsize=7, ncol=2); ax[0].grid(alpha=0.3, which="both")
for gg, d in gdep.items():
    ax[1].semilogy(Tg, np.maximum(1 - d["ramp"], 1e-9), "-^", ms=4, color=f"C{int(2*gg)-1}", label=f"ramp, $g={gg}$")
    ax[1].semilogy(Tg, np.maximum(1 - d["r48"], 1e-9), "--o", ms=3, color=f"C{int(2*gg)-1}", label=f"$4\\to8$ ($\\mu_{{\\rm buf}}=2$), $g={gg}$")
    ax[1].axhline(spec[gg]["inf_ans"], color=f"C{int(2*gg)-1}", ls=":", lw=0.8)
ax[1].set_xlabel("$T$"); ax[1].set_ylabel("$1-p(T)$"); ax[1].set_title("coupling dependence (dotted: 16-pair ansatz infidelity)"); ax[1].legend(fontsize=7, ncol=2); ax[1].grid(alpha=0.3, which="both")
plt.tight_layout()
if OUT: plt.savefig(os.path.join(OUT, "fig21.png"), dpi=150)
else: plt.show()

# ------------------------------------------------------------------ (c) circuit level
print(f"\n(c) circuit level at g=1  [{time.time()-t0:.0f} s]")
cn = pl.cnot_counts(k); CN_STEP, CN_CSTEP = cn["step_uncontrolled"], cn["step_controlled"]
pe = sz.pair_encoding_cnots(k); PE_STEP, PE_CSTEP = pe["step_uncontrolled"], pe["step_controlled"]
CN_A48, CN_A68 = sz.coarse_term_cnots(4, 8, P), sz.coarse_term_cnots(6, 8, P)
PE_A48 = 2 ** 4 - 2                                                   # pair encoding: phase conditioned on the 4 appended qubits
PE_A68 = sz.pair_encoding_cnots(6)["step_uncontrolled"] + sz.pair_encoding_cnots(6)["rotations"] * (2 ** 3 - 2)   # k=6 step, 2 controls per rotation
terms = [D] + hops
def n_needed(Phi0, groups, coefs, T, pex, tol=1e-3):
    out = []
    for ns in (10, 20, 40, 80, 160, 320):
        p = abs(np.vdot(GS, pl.trotter_path(Phi0, groups, coefs, T, ns, order=2))) ** 2; out.append((ns, p))
        if abs(p - pex) <= tol: return ns, p, out
    return ns, p, out
print("second-order steps reproducing the exact-slice overlap within 1e-3:")
rA2 = sz.refine_pair(4, k, P, g, mu_buf=2.0, M=160); gA2 = [rA2["Hle"], D - rA2["mu"] * np.eye(70)] + hops; cA = [c2f, s2f] + [s2f] * len(hops)
rA1 = sz.refine_pair(4, k, P, g, mu_buf=1.0, M=160); gA1 = [rA1["Hle"], D - rA1["mu"] * np.eye(70)] + hops
rA4 = sz.refine_pair(4, k, P, g, mu_buf=4.0, M=160); gA4 = [rA4["Hle"], D - rA4["mu"] * np.eye(70)] + hops
a6 = sz.pair_double_ansatz(6, P, g); Dc, hopsc = sz.pair_terms(6, g, P)
r6 = sz.refine_pair(6, k, P, g, mu_buf=2.0, M=160, psi_low=a6["psi"]); Pm = r6["Pm"]
g6 = [Pm @ (Dc - r6["mu"] * np.eye(15)) @ Pm.T] + [Pm @ h @ Pm.T for h in hopsc] + [D - r6["mu"] * np.eye(70)] + hops; c6 = [c2f] * (1 + len(hopsc)) + [s2f] * (1 + len(hops))
gR, cR = sz.ramp_terms(k, g, P)
steps = {}
for name, (Phi0, groups, coefs, T, pex, per_step, pe_step) in {
        "refined 4->8, mu_buf=2, T=8": (rA2["Phi0"], gA2, cA, 8.0, sz.refine_pair_overlaps(rA2, [8.0])[0], CN_STEP + CN_A48, PE_STEP + PE_A48),
        "refined 4->8, mu_buf=2, T=10": (rA2["Phi0"], gA2, cA, 10.0, sz.refine_pair_overlaps(rA2, [10.0])[0], CN_STEP + CN_A48, PE_STEP + PE_A48),
        "refined 4->8, mu_buf=1, T=25": (rA1["Phi0"], gA1, cA, 25.0, sz.refine_pair_overlaps(rA1, [25.0])[0], CN_STEP + CN_A48, PE_STEP + PE_A48),
        "refined 4->8, mu_buf=4, T=6": (rA4["Phi0"], gA4, cA, 6.0, sz.refine_pair_overlaps(rA4, [6.0])[0], CN_STEP + CN_A48, PE_STEP + PE_A48),
        "refined 6->8 (8-pair coarse ansatz, split coarse groups), mu_buf=2, T=8": (r6["Phi0"], g6, c6, 8.0, sz.refine_pair_overlaps(r6, [8.0])[0], CN_STEP + CN_A68, PE_STEP + PE_A68),
        "HF ramp, T=4": (hf, gR, cR, 4.0, ovR[1], CN_STEP, PE_STEP),
        "HF ramp, T=6": (hf, gR, cR, 6.0, ovR[2], CN_STEP, PE_STEP)}.items():
    ns, p, out = n_needed(Phi0.astype(complex), groups, coefs, T, pex)
    steps[name] = dict(ns=ns, p=p, psi=pl.trotter_path(Phi0.astype(complex), groups, coefs, T, ns, order=2), cost=2 * ns * per_step, pe_cost=2 * ns * pe_step)
    print(f"  {name}: exact {pex:.5f}; " + " ".join(f"n_s={n}: {pp:.5f}" for n, pp in out) + f"  -> n_s = {ns}, {steps[name]['cost']} CNOTs (pair encoding {steps[name]['pe_cost']})")

# averaged post-selection map of the digital filter (quadrature over half-Gaussian times, E_obj = E0)
from numpy.polynomial.hermite_e import hermegauss
xg, wq = hermegauss(60); wq = wq / np.sqrt(2 * np.pi)
def make_Ufun(terms_, DT_):
    eig_terms = [np.linalg.eigh(T) for T in terms_]; dim = terms_[0].shape[0]
    def Ustep(dt):
        U = np.eye(dim, dtype=complex)
        for lam, Vt in eig_terms: U = (Vt * np.exp(-1j * lam * dt / 2)) @ Vt.conj().T @ U
        for lam, Vt in reversed(eig_terms): U = (Vt * np.exp(-1j * lam * dt / 2)) @ Vt.conj().T @ U
        return U
    def Ufun(t):
        n = max(1, int(np.ceil(t / DT_))); return np.linalg.matrix_power(Ustep(t / n), n)
    return Ufun
def filter_floor(terms_, E0_, GS_, SIG_, DT_, iters=400):
    """Dominant eigenoperator of the averaged post-selection map K(rho) =
    E_t[K_t rho K_t^+] by power iteration (K_t = (1 + e^{iE0 t} U(t))/2 with the
    product-formula U), started from the ground-state projector."""
    Uf = make_Ufun(terms_, DT_); dim = len(GS_)
    Ks = [0.5 * (np.eye(dim) + np.exp(1j * E0_ * abs(SIG_ * xk)) * Uf(abs(SIG_ * xk))) for xk in xg]
    rho = np.outer(GS_, GS_).astype(complex); ev = 1.0
    for it in range(iters):
        new = sum(wk * (Kt @ rho @ Kt.conj().T) for Kt, wk in zip(Ks, wq))
        ev_new = np.real(np.trace(new)); rho = new / ev_new
        if it > 20 and abs(ev_new - ev) < 1e-10: break
        ev = ev_new
    return ev_new, 1 - np.real(GS_ @ rho @ GS_)
SIG = 4.0                                      # sigma * gap = 10.4 at g=1 (paper I: 4 x 2.3 = 9.2)
print(f"averaged post-selection map K of the digital filter (second order, sigma={SIG}, |t| half-Gaussian, E_obj = E0): dominant eigenvalue and asymptotic infidelity")
for DTm in (0.25, 0.1, 0.05):
    ev, fl = filter_floor(terms, E0, GS, SIG, DTm); print(f"  g=1, dt_max = {DTm}: eigenvalue {ev:.5f}, floor {fl:.2e}")
DT = 0.1; Ufun = make_Ufun(terms, DT)

# matched comparison, g = 1
NSAMP, M_MAX, SEED = 40, 6, 3
a8 = spec[1.0]["ans"]
inputs = {"HF determinant": (hf.astype(complex), 0, 0),
          "16-pair ansatz": (a8["psi"].astype(complex), a8["cnots"], 4 * 16),
          "prolonged k=6 ansatz": (r6["Phi0"].astype(complex), a6["cnots"], 4 * 8)}
for name in ("HF ramp, T=4", "HF ramp, T=6", "refined 4->8, mu_buf=2, T=8", "refined 4->8, mu_buf=2, T=10", "refined 4->8, mu_buf=1, T=25", "refined 4->8, mu_buf=4, T=6"):
    d = steps[name]; inputs[f"{name} (n_s={d['ns']})"] = (d["psi"], d["cost"], d["pe_cost"])
d = steps["refined 6->8 (8-pair coarse ansatz, split coarse groups), mu_buf=2, T=8"]
inputs[f"refined 6->8 (ansatz), mu_buf=2, T=8 (n_s={d['ns']})"] = (d["psi"], a6["cnots"] + d["cost"], 32 + d["pe_cost"])

def matched(inputs, Ufun, E0, GS, cstep_full, cstep_pe, label, DT, SIG):
    rng = np.random.default_rng(SEED); tsets = [np.abs(rng.normal(0.0, SIG, M_MAX)) for _ in range(NSAMP)]
    nsteps = np.array([[max(1, int(np.ceil(t / DT))) for t in ts] for ts in tsets])
    rows = {}; rngb = np.random.default_rng(5)
    for name, (v0, cprep, cprep_pe) in inputs.items():
        p = abs(np.vdot(GS, v0)) ** 2; F = np.zeros((NSAMP, M_MAX)); Acc = np.zeros((NSAMP, M_MAX))
        for s_, ts in enumerate(tsets): F[s_], Acc[s_] = pl.rodeo_track_U(v0, Ufun, ts, E0, GS)
        Fw = (F * Acc).sum(0) / Acc.sum(0); Am = Acc.mean(0)
        boots = np.array([1 - ((F * Acc)[idx].sum(0) / Acc[idx].sum(0)) for idx in (rngb.integers(0, NSAMP, NSAMP) for _ in range(400))]); rse = boots.std(0) / np.maximum(1 - Fw, 1e-30)
        reach = np.concatenate([np.ones((NSAMP, 1)), Acc[:, :-1]], axis=1)
        row = dict(p=p, cprep=cprep, inf=np.concatenate([[1 - p], 1 - Fw]), rse=np.concatenate([[0], rse]), acc=np.concatenate([[1], Am]))
        for key, cs, cp in (("full", cstep_full, cprep), ("pe", cstep_pe, cprep_pe)):
            cyc = nsteps * cs * 2
            row["cacc_" + key] = np.concatenate([[float(cp)], (cp + np.cumsum(reach * cyc, axis=1).mean(0)) / Am])
            row["catt_" + key] = np.concatenate([[float(cp)], cp + np.cumsum(cyc, axis=1).mean(0)])
        rows[name] = row
    print(f"\n{label}: Trotterized rodeo (second order, dt_max={DT}, sigma={SIG}, E_obj=E0, {NSAMP} time sets; {nsteps.mean():.1f} steps = {2*cstep_full*nsteps.mean():.0f} CNOTs per cycle; bootstrap rse in brackets)")
    print(f"{'input':52s} {'p':>8s} {'prep':>7s} | 1-F_M for M=0..6 | acc(6) | CNOTs per accepted preparation, M=0..6")
    for name, d in rows.items():
        print(f"{name:52s} {d['p']:8.5f} {d['cprep']:7d} | " + " ".join(f"{d['inf'][m]:.1e}({100*d['rse'][m]:2.0f}%)" for m in range(M_MAX + 1)) + f" | {d['acc'][-1]:.3f} | " + " ".join(f"{c:.2e}" for c in d["cacc_full"]))
    for key, lab in (("full", "2k-qubit encoding of paper I"), ("pe", "one-qubit-per-level pair encoding")):
        print(f"  minimum CNOTs per accepted preparation at a target infidelity ({lab}):")
        for tgt in (1e-3, 3e-4, 1e-4):
            best = []
            for name, d in rows.items():
                ok = np.where(d["inf"] <= tgt)[0]
                if ok.size: m = ok[np.argmin(d["cacc_" + key][ok])]; best.append((d["cacc_" + key][m], name, m, d["inf"][m]))
            best.sort(); print(f"    {tgt:.0e}: " + "; ".join(f"{nm} M={m} {c:.2e} ({i:.1e})" for c, nm, m, i in best))
    return rows

rows1 = matched(inputs, Ufun, E0, GS, CN_CSTEP, PE_CSTEP, "g = 1", DT, SIG)

# g = 2: the strong-coupling comparison (ramp, refinement 4->8 mu_buf=2 at T=25, the ansatz, the determinant); sigma from sigma*gap ~ 9
g2 = 2.0; H2 = sz.H_pair(k, g2, P); w2, V2 = np.linalg.eigh(H2); E02, GS2 = w2[0], V2[:, 0]; D2, hops2 = sz.pair_terms(k, g2, P); terms2 = [D2] + hops2
SIG2 = 1.5
print(f"\ng = 2 (E0 = {E02:.5f}, p_HF = {abs(GS2 @ hf)**2:.4f}, gap {w2[1]-w2[0]:.3f}, sigma = {SIG2}): digital-filter floors")
floors2 = {}
for DTm in (0.1, 0.05):
    ev, fl = filter_floor(terms2, E02, GS2, SIG2, DTm); floors2[DTm] = fl; print(f"  dt_max = {DTm}: eigenvalue {ev:.5f}, floor {fl:.2e}")
DT2 = 0.1 if floors2[0.1] < 1e-4 else 0.05; Ufun2 = make_Ufun(terms2, DT2)
print(f"  operating point at g=2: dt_max = {DT2}")
GS_save = GS; GS = GS2
rA = sz.refine_pair(4, k, P, g2, mu_buf=2.0, M=160); gA = [rA["Hle"], D2 - rA["mu"] * np.eye(70)] + hops2
rA8 = sz.refine_pair(4, k, P, g2, mu_buf=8.0, M=160); gA8 = [rA8["Hle"], D2 - rA8["mu"] * np.eye(70)] + hops2
gR2, cR2 = sz.ramp_terms(k, g2, P); eR2 = sz._slice_eigs(D2, H2, 160)
inputs2 = {"HF determinant": (hf.astype(complex), 0, 0), "16-pair ansatz": (spec[2.0]["ans"]["psi"].astype(complex), 768, 64)}
for name, (Phi0, groups, coefs, T, pex, per_step, pe_step) in {
        "HF ramp, T=8": (hf, gR2, cR2, 8.0, abs(np.vdot(GS2, sz._evolve(hf, eR2, 8.0))) ** 2, CN_STEP, PE_STEP),
        "HF ramp, T=12": (hf, gR2, cR2, 12.0, abs(np.vdot(GS2, sz._evolve(hf, eR2, 12.0))) ** 2, CN_STEP, PE_STEP),
        "refined 4->8, mu_buf=2, T=25": (rA["Phi0"], gA, cA, 25.0, sz.refine_pair_overlaps(rA, [25.0])[0], CN_STEP + CN_A48, PE_STEP + PE_A48),
        "refined 4->8, mu_buf=8, T=8": (rA8["Phi0"], gA8, cA, 8.0, sz.refine_pair_overlaps(rA8, [8.0])[0], CN_STEP + CN_A48, PE_STEP + PE_A48)}.items():
    ns, p, out = n_needed(Phi0.astype(complex), groups, coefs, T, pex)
    psi = pl.trotter_path(Phi0.astype(complex), groups, coefs, T, ns, order=2)
    inputs2[f"{name} (n_s={ns})"] = (psi, 2 * ns * per_step, 2 * ns * pe_step)
    print(f"  {name}: exact {pex:.5f}; n_s = {ns}: {p:.5f}, {2*ns*per_step} CNOTs")
rows2 = matched(inputs2, Ufun2, E02, GS2, CN_CSTEP, PE_CSTEP, "g = 2", DT2, SIG2)
GS = GS_save

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
mks = ["o-", "s-", "d-", "^-", "v-", "P-", "X-", "*-", "<-", ">-"]
for (name, d), mk in zip(rows1.items(), mks):
    ax[0].loglog(d["cacc_full"][1:] if d["cprep"] == 0 else d["cacc_full"], d["inf"][1:] if d["cprep"] == 0 else d["inf"], mk, ms=4, label=f"{name} ($p={d['p']:.4f}$)")
for tgt in (1e-3, 3e-4, 1e-4): ax[0].axhline(tgt, color="gray", ls=":", lw=0.8)
ax[0].set_xlabel("CNOTs per accepted preparation (input + filter)"); ax[0].set_ylabel("ensemble infidelity $1-\\mathcal{F}_M$"); ax[0].set_title(f"$k=8$, $N=8$, $g=1$: ten inputs to the same filter ($\\delta t={DT}$), $M=0\\ldots6$"); ax[0].legend(fontsize=6); ax[0].grid(alpha=0.3, which="both")
for (name, d), mk in zip(rows2.items(), mks):
    ax[1].loglog(d["cacc_full"][1:] if d["cprep"] == 0 else d["cacc_full"], d["inf"][1:] if d["cprep"] == 0 else d["inf"], mk, ms=4, label=f"{name} ($p={d['p']:.4f}$)")
for tgt in (1e-3, 3e-4, 1e-4): ax[1].axhline(tgt, color="gray", ls=":", lw=0.8)
ax[1].set_xlabel("CNOTs per accepted preparation (input + filter)"); ax[1].set_ylabel("$1-\\mathcal{F}_M$"); ax[1].set_title(f"$g=2$ ($\\delta t={DT2}$, $\\sigma={SIG2}$)"); ax[1].legend(fontsize=6); ax[1].grid(alpha=0.3, which="both")
plt.tight_layout()
if OUT: plt.savefig(os.path.join(OUT, "fig22.png"), dpi=150); print("wrote", OUT)
else: plt.show()
print(f"[{time.time()-t0:.0f} s]")
