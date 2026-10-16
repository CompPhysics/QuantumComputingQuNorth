#!/usr/bin/env python3
"""ADAPT-VQE comparator (notebook section EXT:ADAPT, fig 23).

ADAPT-VQE (Grimsley et al. 2019) in the seniority-zero pair basis of
`pairinglib.seniority`: the pool is either the P(k-P) occupied-pair ->
virtual-pair generators of the restricted UCCSD ('pairs') or all k(k-1)/2
pair-hopping generators ('generalized'); at every iteration the operator with
the largest energy gradient is appended (repeats allowed) and all amplitudes
are re-optimized.  Each generator is a 48-CNOT block in the 2k-qubit encoding
of paper I (4 CNOTs in the one-qubit-per-level pair encoding).
Cases: (k, N) = (4, 4) and (6, 4) at g = 1, and (8, 8) at g = 1 and g = 2.
For (8, 8) the ADAPT states at the first iteration below 1e-4 are also fed to
the Trotterized rodeo filter of Section 20 (M = 0..6).
Usage: PYTHONPATH=src python3 scripts/adapt_checks.py [outdir]
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
TARGETS = (1e-3, 3e-4, 1e-4, 1e-5)
cases = [(4, 2, 1.0), (6, 2, 1.0), (8, 4, 1.0), (8, 4, 2.0)]
runs = {}
print("ADAPT-VQE in the pair basis (gradient tolerance 1e-6, at most 40 operators; 48 CNOTs per operator in the 2k-qubit encoding, 4 in the pair encoding)")
print(f"{'(k,N,g)':>10s} {'pool':>12s} {'|pool|':>6s} {'ops':>4s} {'1-F final':>10s} | operators (CNOTs) to reach 1e-3 / 3e-4 / 1e-4 / 1e-5")
for (k, P, g) in cases:
    H = sz.H_pair(k, g, P); GS = np.linalg.eigh(H)[1][:, 0]
    for pool in ("pairs", "generalized"):
        r = sz.adapt_vqe(k, P, g, pool=pool, grad_tol=1e-6, max_ops=40, H=H, target=GS)
        h = r["history"]; runs[(k, P, g, pool)] = r
        reach = []
        for tgt in TARGETS:
            hit = next((d for d in h if d["infidelity"] <= tgt), None)
            reach.append(f"{hit['n_ops']} ({hit['cnots']})" if hit else "-")
        print(f"{'(%d,%d,%.0f)' % (k, 2*P, g):>10s} {pool:>12s} {len(r['pool']):6d} {len(r['ops']):4d} {h[-1]['infidelity']:10.1e} | " + " / ".join(reach))
        if k == 4:
            print(f"             selected: {r['ops']}  E = {r['E']:.10f}")
# the restricted single-layer ansatz of Secs. 19-20 for comparison
for (k, P, g) in cases:
    a = sz.pair_double_ansatz(k, P, g); GS = np.linalg.eigh(sz.H_pair(k, g, P))[1][:, 0]
    print(f"  single-layer restricted ansatz (k={k}, N={2*P}, g={g:.0f}): {a['n_params']} operators, {a['cnots']} CNOTs, 1-F = {1-abs(np.vdot(GS, a['psi']))**2:.2e}")
print(f"[{time.time()-t0:.0f} s]")

# ---- (8, 8): ADAPT states through the Trotterized filter of Section 20 (same settings), M = 0..6
k, P = 8, 4
NSAMP, M_MAX, SEED = 40, 6, 3
cn = pl.cnot_counts(k); CN_CSTEP = cn["step_controlled"]
print("\n(8, 8): ADAPT states (generalized pool, first iteration below 1e-4) through the Trotterized filter of Section 20")
adapt_rows = {}
for g, SIG, DT in ((1.0, 4.0, 0.1), (2.0, 1.5, 0.05)):
    H = sz.H_pair(k, g, P); w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]
    D, hops = sz.pair_terms(k, g, P); terms = [D] + hops
    eig_terms = [np.linalg.eigh(T) for T in terms]
    def Ustep(dt):
        U = np.eye(70, dtype=complex)
        for lam, Vt in eig_terms: U = (Vt * np.exp(-1j * lam * dt / 2)) @ Vt.conj().T @ U
        for lam, Vt in reversed(eig_terms): U = (Vt * np.exp(-1j * lam * dt / 2)) @ Vt.conj().T @ U
        return U
    def Ufun(t):
        n = max(1, int(np.ceil(t / DT))); return np.linalg.matrix_power(Ustep(t / n), n)
    r = runs[(k, P, g, "generalized")]; h = r["history"]
    hit = next(d for d in h if d["infidelity"] <= 1e-4); n_hit = hit["n_ops"]; psi = hit["psi"]; E = hit["E"]
    v0 = psi.astype(complex); cprep = 48 * n_hit
    rng = np.random.default_rng(SEED); tsets = [np.abs(rng.normal(0.0, SIG, M_MAX)) for _ in range(NSAMP)]
    nsteps = np.array([[max(1, int(np.ceil(t / DT))) for t in ts] for ts in tsets]); cyc = nsteps * CN_CSTEP * 2
    F = np.zeros((NSAMP, M_MAX)); Acc = np.zeros((NSAMP, M_MAX))
    for s_, ts in enumerate(tsets): F[s_], Acc[s_] = pl.rodeo_track_U(v0, Ufun, ts, E0, GS)
    Fw = (F * Acc).sum(0) / Acc.sum(0); Am = Acc.mean(0)
    reach = np.concatenate([np.ones((NSAMP, 1)), Acc[:, :-1]], axis=1)
    cacc = np.concatenate([[float(cprep)], (cprep + np.cumsum(reach * cyc, axis=1).mean(0)) / Am])
    inf = np.concatenate([[1 - abs(np.vdot(GS, v0)) ** 2], 1 - Fw])
    adapt_rows[g] = dict(n_ops=n_hit, cprep=cprep, inf=inf, acc=np.concatenate([[1.0], Am]), cacc=cacc, E=E)
    print(f"  g={g:.0f}: {n_hit} operators, {cprep} CNOTs, E = {E:.6f} (E0 {E0:.6f}), p = {1-inf[0]:.6f}; 1-F_M for M=0..6: " + " ".join(f"{x:.1e}" for x in inf)
          + f"; acceptance(6) {Am[-1]:.3f}; CNOTs per accepted preparation " + " ".join(f"{c:.1e}" for c in cacc))

fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
for (kk, PP, gg, pool), r in runs.items():
    h = r["history"]; n = [d["n_ops"] for d in h]; inf_ = [max(d["infidelity"], 1e-16) for d in h]
    ax[0].semilogy(n, inf_, "o-" if pool == "generalized" else "s--", ms=3, label=f"$k={kk}$, $N={2*PP}$, $g={gg:.0f}$, {pool} pool")
ax[0].set_xlabel("ADAPT iteration (number of pair-rotation operators)"); ax[0].set_ylabel("infidelity $1-F$"); ax[0].set_ylim(1e-9, 1)
ax[0].set_title("ADAPT-VQE in the seniority-zero pair basis"); ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3, which="both")
for (kk, PP, gg, pool), r in runs.items():
    if pool != "generalized": continue
    h = r["history"]; c = [max(d["cnots"], 1) for d in h]; inf_ = [max(d["infidelity"], 1e-16) for d in h]
    ax[1].loglog(c, inf_, "o-", ms=3, label=f"ADAPT, $k={kk}$, $N={2*PP}$, $g={gg:.0f}$")
for g, d in adapt_rows.items():
    ax[1].loglog(d["cacc"], np.maximum(d["inf"], 1e-16), "x:", ms=6, label=f"ADAPT state + Trotterized filter, $k=8$, $g={g:.0f}$ ($M=0\\ldots6$)")
for tgt in (1e-3, 3e-4, 1e-4): ax[1].axhline(tgt, color="gray", ls=":", lw=0.8)
ax[1].set_xlabel("CNOTs (2k-qubit encoding; filter cycles included for the dotted curves)"); ax[1].set_ylabel("$1-F$"); ax[1].set_ylim(1e-9, 1)
ax[1].set_title("cost versus accuracy of the adaptive variational route"); ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3, which="both")
plt.tight_layout()
if OUT: plt.savefig(os.path.join(OUT, "fig23.png"), dpi=150); print("wrote", OUT)
else: plt.show()
print(f"[{time.time()-t0:.0f} s]")
