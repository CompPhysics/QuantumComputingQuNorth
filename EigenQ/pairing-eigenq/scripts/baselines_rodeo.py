#!/usr/bin/env python3
"""Matched-observable baselines for the rodeo stage (k=4, N=4, g=1).

Feeds the SAME Trotterised rodeo filter (second order, fixed step dt=0.25,
sigma=4, |t_k| from a Gaussian) four inputs and tracks, cycle by cycle,
  * the ensemble infidelity of the accepted runs (acceptance-weighted mean of
    the post-selected fidelity; the uniform mean over time sets is also given),
  * the cumulative acceptance P_M,
  * CNOTs per accepted preparation = (preparation + rodeo)/P_M, with and
    without early termination after the first failed ancilla measurement.
Inputs: (i) the prolonged k=2 determinant, (ii) the full-space gate-level
UCCSD state, (iii) the Trotterised resolution-refined state at T=16 (T*),
(iv) the same at T=30 (the operating point of Table II).  The refinement cost
INCLUDES the embedded coarse term A = P(H_low-mu)P^T (14 CNOTs per
application for 2->4; see Sec. V.B of the paper).

Usage:  PYTHONPATH=src python3 scripts/baselines_rodeo.py [outdir]
Writes <outdir>/fig17.png and prints the table used in the paper.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
sys.path.insert(0, "src")
import pairinglib as pl

OUT = sys.argv[1] if len(sys.argv) > 1 else None   # directory for fig17.png; None -> plt.show()
N, G, k = 4, 1.0, 4
DT, SIG, NSAMP, M_MAX, SEED = 0.25, 4.0, 40, 6, 3
NS_DENS = 80/30.0                       # Trotter steps per unit adiabatic time (Table II)
CN_STEP, CN_CSTEP = 296, 416            # pairing step, uncontrolled / controlled (Eq. 9)
CN_COARSE = 14                          # 4-qubit conditional phase of the coarse term (2->4)
CN_UCCSD = 1312

nq, st, ix = pl.build_sector(k, N)
H = pl.H_pairing_sparse(k, G, N, st, ix).toarray()
w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]
terms = pl.model_terms(k, G, N)
Ufun = lambda t: pl.trotter_U_dt(k, G, N, t, DT, order=2, terms=terms)[0]

# ---- inputs ---------------------------------------------------------------
Pm = pl.prolong_matrix(2, k, N); pro = Pm @ np.array([1.0]); pro /= np.linalg.norm(pro)
E_u, psi_full = pl.gate_uccsd_state(k, N, G)
ucc = np.array([psi_full[I] for I in st]); ucc /= np.linalg.norm(ucc)
inputs = {}
inputs["prolonged determinant"] = (pro, 0)
inputs["UCCSD (full space)"] = (ucc, CN_UCCSD)
for T in (16, 30):
    ns = int(np.ceil(NS_DENS*T))
    psi, _ = pl.refine_state_trotter(N, k, float(T), ns, order=2)
    inputs[f"refined, T={T} (n_s={ns})"] = (psi, 2*ns*(CN_STEP + CN_COARSE))

# ---- common random times ---------------------------------------------------
rng = np.random.default_rng(SEED)
tsets = [np.abs(rng.normal(0.0, SIG, M_MAX)) for _ in range(NSAMP)]
nsteps = np.array([[max(1, int(np.ceil(t/DT))) for t in ts] for ts in tsets])   # (NSAMP, M)
cyc_cost = nsteps*CN_CSTEP*2                                                       # CNOTs per cycle

rows = {}
for name, (v0, cprep) in inputs.items():
    p = abs(np.vdot(GS, v0))**2
    F = np.zeros((NSAMP, M_MAX)); A = np.zeros((NSAMP, M_MAX))
    for s, ts in enumerate(tsets):
        F[s], A[s] = pl.rodeo_track_U(v0, Ufun, ts, E0, GS)
    Fmean = F.mean(0); Amean = A.mean(0); Fw = (F*A).sum(0)/A.sum(0)
    # cost per accepted preparation after M cycles
    cum = np.cumsum(cyc_cost, axis=1)                       # cost of M cycles, per time set
    no_et = (cprep + cum.mean(0))/Amean                     # run every cycle, then post-select
    # early termination: expected cost = prep + sum_m P(reach m) c_m, P(reach 1)=1
    reach = np.concatenate([np.ones((NSAMP, 1)), A[:, :-1]], axis=1)
    et = (cprep + np.cumsum(reach*cyc_cost, axis=1).mean(0))/Amean
    rows[name] = dict(p=p, cprep=cprep, F=Fmean, Fw=Fw, A=Amean, cost=no_et, cost_et=et)

# ---- table -----------------------------------------------------------------
print(f"k={k}, N={N}, g={G}: rodeo second order, dt={DT}, sigma={SIG}, {NSAMP} time sets (seed {SEED})")
print(f"{'input':28s} {'p':>9s} {'prep CNOT':>10s} | " + " | ".join(f"M={m+1}: 1-F_w  acc  CNOT/acc(ET)" for m in range(M_MAX) if m+1 in (1, 2, 4, 6)))
for name, r in rows.items():
    cells = [f"{1-r['Fw'][m]:.2e} {r['A'][m]:.4f} {r['cost_et'][m]:.0f}" for m in range(M_MAX) if m+1 in (1, 2, 4, 6)]
    print(f"{name:28s} {r['p']:9.6f} {r['cprep']:10d} | " + " | ".join(cells))
print("\nuniform mean of normalised fidelities vs acceptance-weighted (M=2, M=6):")
for name, r in rows.items():
    print(f"  {name:28s} 1-<F>: {1-r['F'][1]:.3e} / {1-r['F'][5]:.3e}   1-F_w: {1-r['Fw'][1]:.3e} / {1-r['Fw'][5]:.3e}"
          f"   cost no-ET vs ET at M=2: {r['cost'][1]:.0f} / {r['cost_et'][1]:.0f}")

# ---- figure: cost versus achieved infidelity ------------------------------
import matplotlib
if OUT: matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
mk = {"prolonged determinant": "o-", "UCCSD (full space)": "s-", "refined, T=16 (n_s=43)": "d-", "refined, T=30 (n_s=80)": "^-"}
Ms = np.arange(0, M_MAX+1)
for name, r in rows.items():
    inf = np.concatenate([[1-r['p']], 1-r['Fw']])
    cost = np.concatenate([[r['cprep']], r['cost_et']])
    m0 = 1 if r['cprep'] == 0 else 0            # a gate-free input has no M=0 point on a log axis
    ax[0].loglog(cost[m0:], inf[m0:], mk[name], label=f"{name} ($p={r['p']:.4f}$)")
    for m in (m0, 2, 6):
        ax[0].annotate(f"$M={m}$", (cost[m], inf[m]), textcoords="offset points", xytext=(5, -9 if m == 6 else 4), fontsize=7)
    ax[1].semilogy(Ms, inf, mk[name], label=name)
ax[0].set_xlabel("CNOTs per accepted preparation (input + filter, early termination)")
ax[0].set_ylabel(r"ensemble infidelity $1-\mathcal{F}_M$")
ax[0].set_title(r"cost versus accuracy, $k=4$, $N=4$, $g=1$ ($\delta t=0.25$, $\sigma=4$)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3, which="both")
ax[1].set_xlabel("rodeo cycles $M$ ($M=0$: unfiltered input)"); ax[1].set_ylabel(r"$1-\mathcal{F}_M$")
ax[1].set_title("same data versus the number of cycles"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, which="both")
plt.tight_layout()
if OUT:
    plt.savefig(os.path.join(OUT, "fig17.png"), dpi=150); print("wrote", os.path.join(OUT, "fig17.png"))
else:
    plt.show()

# ---- particle-hole model: UCCSD + rodeo at the Table V operating points ----
print("\nparticle-hole term: UCCSD input + Trotterised rodeo (M=2), k=4, N=4, g=1")
for f in (0.2, 0.5):
    Hf = pl.H_pairing_sparse(k, G, N, st, ix, f=f).toarray(); wf, Vf = np.linalg.eigh(Hf); E0f, GSf = wf[0], Vf[:, 0]
    E_uf, psi_f = pl.gate_uccsd_state(k, N, G, f=f); uf = np.array([psi_f[I] for I in st]); uf /= np.linalg.norm(uf)
    termsf = pl.model_terms(k, G, N, f=f); cnf = pl.cnot_counts(k, f=f)
    for dtm in (0.25, 0.125):
        Uf = lambda t, dtm=dtm: pl.trotter_U_dt(k, G, N, t, dtm, order=2, terms=termsf)[0]
        rng = np.random.default_rng(SEED); F2 = A2 = FA = 0.0; nst = 0
        for s in range(NSAMP):
            ts = np.abs(rng.normal(0.0, SIG, 2))
            Fc, Ac = pl.rodeo_track_U(uf, Uf, ts, E0f, GSf)
            FA += Fc[-1]*Ac[-1]; A2 += Ac[-1]; nst += sum(max(1, int(np.ceil(t/dtm))) for t in ts)
        inf = 1 - FA/A2; acc = A2/NSAMP; cost = (CN_UCCSD + nst/NSAMP*cnf['step_controlled']*2)/acc
        print(f"  f={f} dt={dtm}: 1-F(UCCSD)={1-abs(np.vdot(GSf, uf))**2:.2e}  UCCSD+rodeo M=2: 1-F={inf:.2e} acc={acc:.4f} "
              f"rodeo CNOT={nst/NSAMP*cnf['step_controlled']*2:.0f} per accepted={cost:.0f}")
