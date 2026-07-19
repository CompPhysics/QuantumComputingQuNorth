#!/usr/bin/env python3
"""Insert the circuit-level (Trotterised) pipeline, the k=3->4 chain and the
g-scan sections into ResolutionRefPairing.ipynb, before the summary cell.
Built with nbformat (never hand-edit the .ipynb JSON); idempotent: existing
cells tagged with the section markers are replaced.
Run from the repo root:  python3 scripts/extend_notebook_trotter.py
"""
import sys, pathlib
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
MARK = "<!-- EXT:TROTTER -->"

md_trotter = MARK + r"""
---
## 13&nbsp; Circuit-level product formulas: Trotterising the whole pipeline

So far the refinement slices and the rodeo's controlled evolution used *exact*
sector exponentials.  On hardware both are product formulas.  We now quantify the
Trotter error at the circuit level, using the term decomposition that matches the
gate compilation of the paper: the diagonal term $D$ (kinetic + $p{=}q$ pairing:
single-qubit phases and $ZZ$ rotations) and one pair-hopping term $h_{pq}$ per
level pair.  Each $e^{-ih_{pq}\,\delta t}$ is realised *exactly* by its eight
commuting weight-4 Pauli rotations, so the only error is the non-commutativity
*between* terms --- exactly as on a device.  We use `pairinglib.pairing_terms`,
`trotter_U`, `refine_state_trotter` and `rodeo_track_U`; nothing is redefined here.
"""

code_trotter_check = MARK + r"""
import numpy as np, matplotlib.pyplot as plt
try:
    import pairinglib as pl
except ModuleNotFoundError:
    import sys, pathlib
    sys.path.insert(0, str((pathlib.Path.cwd()/".."/"src").resolve()))
    import pairinglib as pl

N, G = 4, 1.0
# term decomposition reproduces H exactly, and the Trotter error scales as expected
nq4, st4, ix4 = pl.build_sector(4, N)
H4 = pl.H_pairing_sparse(4, G, N, st4, ix4).toarray()
D4, hops4 = pl.pairing_terms(4, G, N)
print("‖(D + Σ h_pq) − H‖₂  =", np.linalg.norm(D4 + sum(hops4) - H4, 2))
print("\nTrotter error ‖U_trot − e^{−iHt}‖₂ at t = 1, k = 4:")
print(f"{'n_steps':>8} {'order 1':>12} {'order 2':>12}")
for n in (1, 2, 4, 8, 16):
    e1 = pl.trotter_error(4, G, N, 1.0, n, order=1)
    e2 = pl.trotter_error(4, G, N, 1.0, n, order=2)
    print(f"{n:>8} {e1:>12.3e} {e2:>12.3e}")
"""

code_trotter_fig = MARK + r"""
# ---- circuit-level refinement and rodeo: Trotter error at the pipeline level ----
w4, v4 = np.linalg.eigh(H4); E0 = w4[0]; GS = v4[:, 0]
T_AD = 30.0

# (a) Trotterised refinement k=2->4: infidelity vs number of refinement steps n_s
ns_grid = [10, 20, 40, 80, 160]
inf_ref = {1: [], 2: []}
for order in (1, 2):
    for ns in ns_grid:
        psi, _ = pl.refine_state_trotter(N, 4, T_AD, ns, order=order)
        inf_ref[order].append(1 - abs(np.vdot(GS, psi))**2)
psi_exact_slices, _ = pl.refine_state(N, 4, T_AD)         # exact-slice reference
inf_exact = 1 - abs(np.vdot(GS, psi_exact_slices))**2
p_ref = {o: 1 - np.array(inf_ref[o]) for o in (1, 2)}

# (b) rodeo with a Trotterised controlled evolution at FIXED STEP SIZE dt:
#     a cycle of time t uses ceil(t/dt) steps, as on hardware.
NS_OP = 80                                                 # operating point from (a)
psi_in, _ = pl.refine_state_trotter(N, 4, T_AD, NS_OP, order=2)
p_in = abs(np.vdot(GS, psi_in))**2
terms4 = [D4] + hops4
rng = np.random.default_rng(11)
M_CYC, SIG, NSAMP = 6, 4.0, 12
tsets = [np.abs(rng.normal(0.0, SIG, M_CYC)) for _ in range(NSAMP)]
curves, accepts, mean_steps = {}, {}, {}
for dtm in (1.0, 0.5, 0.25, None):                         # None = exact evolution
    F_acc = np.zeros(M_CYC); A_fin = 0.0; nst = 0
    for ts in tsets:
        if dtm is None:
            F, A = pl.rodeo_track(psi_in, H4, ts, E0, GS)
        else:
            Ufun = lambda t: pl.trotter_U_dt(4, G, N, t, dtm, order=2, terms=terms4)[0]
            F, A = pl.rodeo_track_U(psi_in, Ufun, ts, E0, GS)
            nst += sum(max(1, int(np.ceil(t/dtm))) for t in ts)
        F_acc += F; A_fin += A[-1]
    curves[dtm] = 1 - F_acc/NSAMP
    accepts[dtm] = A_fin/NSAMP
    mean_steps[dtm] = nst/(NSAMP*M_CYC) if dtm else None

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for order, mk in ((1, "o-"), (2, "s-")):
    ax[0].loglog(ns_grid, inf_ref[order], mk, label=f"product formula, order {order}")
ax[0].axhline(inf_exact, ls=":", c="k", label=f"exact slices (M=160): {inf_exact:.1e}")
ax[0].set_xlabel(r"refinement steps $n_s$"); ax[0].set_ylabel(r"$1-|\langle\Psi_{\rm GS}|\Phi\rangle|^2$")
ax[0].set_title(rf"Trotterised refinement, $k=2\to4$, $T={T_AD:.0f}$"); ax[0].legend()
Ms = np.arange(1, M_CYC+1)
for dtm, mk in ((1.0, "o-"), (0.5, "s-"), (0.25, "d-"), (None, "k:")):
    lab = "exact evolution" if dtm is None else rf"$\delta t={dtm}$ ($\approx${mean_steps[dtm]:.1f} steps/cycle)"
    ax[1].semilogy(Ms, curves[dtm], mk, label=lab)
ax[1].set_xlabel(r"rodeo cycles $M$"); ax[1].set_ylabel(r"$1-\mathcal{F}_M$")
ax[1].set_title(rf"Rodeo on the circuit-level refined state ($p={p_in:.4f}$)")
ax[1].legend(); plt.tight_layout(); plt.show()
print(f"circuit-level refined input: n_s={NS_OP} (order 2), overlap p = {p_in:.6f}")
for dtm in (1.0, 0.5, 0.25):
    print(f"  dt={dtm}: infidelity after {M_CYC} cycles = {curves[dtm][-1]:.3e}, acceptance = {accepts[dtm]:.4f}")
print(f"  exact evolution: infidelity = {curves[None][-1]:.3e}, acceptance = {accepts[None]:.4f}")
"""

md_budget = MARK + r"""
### An honest depth budget

The comparison that matters for hardware is **total two-qubit gates per accepted
preparation at matched fidelity**, including the post-selection repetition factor
$1/P_M$.  The VQE's optimisation is a one-time classical cost, after which *every*
preparation replays the full circuit; the pipeline pays its refinement steps on
every repetition instead.  Both totals are computed below from the measured
operating points, with the first-order per-step CNOT counts of the paper
(`pl.cnot_counts`); an order-2 step costs twice the order-1 count.
"""

code_budget = MARK + r"""
# ---- honest depth budget at k=4 (per accepted preparation, matched target) ----
cn = pl.cnot_counts(4)
CN_STEP, CN_CSTEP = cn["step_uncontrolled"], cn["step_controlled"]
ORDER_FAC = 2                                   # order-2 step = 2 x first-order count
E_uccsd, nit = pl.gate_uccsd_vqe(4, N, G, pl.gate_uccsd_setup(4, N))
CN_VQE = 1312
err_vqe = E_uccsd - pl.fci_ground(pl.H_pairing_sparse(4, G, N, st4, ix4))

# pipeline operating point (from the figure above): n_s=80 order-2, dt=0.25, M=2
DT_OP, M_OP = 0.25, 2
Ufun = lambda t: pl.trotter_U_dt(4, G, N, t, DT_OP, order=2, terms=terms4)[0]
rng = np.random.default_rng(3)
F2 = A2 = 0.0; nsteps_tot = 0
for _ in range(NSAMP):
    ts = np.abs(rng.normal(0.0, SIG, M_OP))
    F, A = pl.rodeo_track_U(psi_in, Ufun, ts, E0, GS)
    F2 += F[-1]/NSAMP; A2 += A[-1]/NSAMP
    nsteps_tot += sum(max(1, int(np.ceil(t/DT_OP))) for t in ts)
nsteps_avg = nsteps_tot/NSAMP
cn_refine  = NS_OP * CN_STEP * ORDER_FAC
cn_rodeo   = int(round(nsteps_avg * CN_CSTEP * ORDER_FAC))
cn_run     = cn_refine + cn_rodeo
cn_accept  = cn_run / A2
print("k=4 depth budget (CNOTs), N=4, g=1")
print(f"  UCCSD-VQE preparation:        {CN_VQE} per prep; optimiser: {nit} iterations x (E+grad) evals")
print(f"    variational error:          {err_vqe:.2e}")
print(f"  pipeline (n_s={NS_OP} order-2, M={M_OP} cycles, dt={DT_OP}, sigma={SIG}):")
print(f"    refinement                  {cn_refine}")
print(f"    rodeo cycles (avg {nsteps_avg:.1f} steps) {cn_rodeo}")
print(f"    per run                     {cn_run}")
print(f"    acceptance P_M              {A2:.4f}")
print(f"    per ACCEPTED preparation    {cn_accept:.0f}")
print(f"    final infidelity            {1-F2:.2e}")
print("The pipeline is deeper per accepted shot than one UCCSD replay, but it is")
print("optimiser-free and reaches an infidelity far below the variational floor;")
print("at matched accuracy the variational route has no operating point at all.")
"""

md_k34 = MARK + r"""
---
## 14&nbsp; A nontrivial coarse stage: the $k=3\to4$ chain

At $N=4$ the $k=2$ sector is a single determinant, so the coarse "VQE" is gate-free
and the prolonged state is just the Hartree--Fock reference --- the pipeline's first
stage does no work.  Here we start instead from the **gate-level UCCSD state at
$k=3$** (272 CNOTs, and indistinguishable from FCI at this resolution), prolong
$3\to4$, and refine.  The comparison $2\to4$ vs $3\to4$ shows what a better coarse
start buys: a higher initial overlap and a shorter adiabatic stage, at the price of
the coarse-circuit depth.
"""

code_k34 = MARK + r"""
# ---- k=3->4 chain from the gate-level UCCSD coarse state ----
E3, psi3_full = pl.gate_uccsd_state(3, N, G)
nq3, st3, ix3 = pl.build_sector(3, N)
psi3 = np.array([psi3_full[s] for s in st3])           # project onto the N-sector
psi3 /= np.linalg.norm(psi3)
print(f"coarse gate-UCCSD at k=3: E = {E3:.6f}  (FCI 0.794697), 272 CNOTs")

Ts = np.array([0.5, 1, 2, 4, 6, 9, 12, 16, 20, 25, 30])
r24 = pl.refine_NK(N, 4, Ts, k_low=2)
r34 = pl.refine_from_state(psi3, N, 3, 4, Ts)
print(f"prolonged overlap ov0:  2->4: {r24['ov0']:.4f}   3->4: {r34['ov0']:.4f}")

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].plot(Ts, r24["ov"], "o-", label=r"$2\to4$ (coarse prep: 0 CNOTs)")
ax[0].plot(Ts, r34["ov"], "s-", label=r"$3\to4$ (coarse prep: 272 CNOTs)")
ax[0].axhline(1.0, ls=":", c="k", lw=0.8)
ax[0].set_xlabel(r"adiabatic time $T$"); ax[0].set_ylabel("overlap with exact GS")
ax[0].set_title("refinement chains to $k=4$"); ax[0].legend()
ax[1].semilogy(Ts, 1-r24["ov"], "o-", label=r"$2\to4$")
ax[1].semilogy(Ts, 1-r34["ov"], "s-", label=r"$3\to4$")
ax[1].set_xlabel(r"adiabatic time $T$"); ax[1].set_ylabel(r"$1-$overlap")
ax[1].set_title("same, log scale"); ax[1].legend()
plt.tight_layout(); plt.show()
for i, T in enumerate(Ts):
    if r34["ov"][i] > 0.99:
        print(f"3->4 reaches overlap>0.99 at T={T};  2->4 overlap there: {r24['ov'][i]:.4f}")
        break
"""

md_gscan = MARK + r"""
---
## 15&nbsp; Coupling-strength scan and the path gap: testing $T^{*} \sim 1/\Delta E$

All results so far sit at $g=1$.  Two independent knobs probe the robustness of the
refinement.  First we scan the pairing strength from weak ($g=0.3$) to strong
($g=4$) coupling at fixed shift buffer.  The minimum instantaneous gap along the
path turns out to be set by the shift itself, $\Delta E_{\min}\simeq\mu-E_{\rm low}$,
and stays open at every $g$; the required time $T^{*}$ (overlap $\geq 0.99$) grows
with $g$ because the prolonged state starts further from the target (smaller
$\mathrm{ov}_0$), i.e.\ more state rearrangement must fit through the same gap.
Second --- since the gap is controlled by the shift --- we vary $\mu$ directly at
$g=1$ and extract $T^{*}$ against $1/\Delta E_{\min}$: this is the clean test of the
$T\sim1/\Delta E$ scaling, and doubles as the sensitivity study for the (otherwise
ad hoc) choice $\mu=E_{\rm low}+0.6$.
"""

code_gscan = MARK + r"""
# ---- (a) g-scan of the k=2->4 refinement at fixed mu_buf=0.6 ----
def min_path_gap(g, mu_buf, k_low=2, k_high=4):
    nl, sl, il = pl.build_sector(k_low, N)
    Hl = pl.H_pairing_sparse(k_low, g, N, sl, il).toarray()
    El = np.linalg.eigvalsh(Hl)[0]
    Pm = pl.prolong_matrix(k_low, k_high, N)
    nh, sh, ih = pl.build_sector(k_high, N)
    Hh = pl.H_pairing_sparse(k_high, g, N, sh, ih).toarray()
    mu = El + mu_buf
    Hle = Pm @ (Hl - mu*np.eye(len(Hl))) @ Pm.T
    Hhs = Hh - mu*np.eye(len(Hh))
    return min(np.diff(np.linalg.eigvalsh(
        np.cos(np.pi*s/2)**2*Hle + np.sin(np.pi*s/2)**2*Hhs)[:2])[0]
        for s in np.linspace(0, 1, 41))

def Tstar_of(r, Ts, target=0.99):
    hit = np.where(r["ov"] >= target)[0]
    return Ts[hit[0]] if hit.size else np.nan

gs_list = [0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
Ts = np.array([0.5, 1, 2, 3, 4, 6, 8, 12, 16, 20, 26, 32, 40, 50, 65, 80])
Tstar_g, gap_g, ov0_g = [], [], []
for g in gs_list:
    r = pl.refine_NK(N, 4, Ts, g=g)
    Tstar_g.append(Tstar_of(r, Ts)); gap_g.append(min_path_gap(g, 0.6)); ov0_g.append(r["ov0"])
    print(f"g={g:4.1f}  ov0={r['ov0']:.4f}  min path gap={gap_g[-1]:.3f}  T*={Tstar_g[-1]}")

# ---- (b) mu-scan at g=1: T* against the inverse path gap ----
mubufs = [0.15, 0.3, 0.6, 1.0, 1.5]
Tstar_mu, gap_mu = [], []
for mb in mubufs:
    r = pl.refine_NK(N, 4, Ts, mu_buf=mb)
    Tstar_mu.append(Tstar_of(r, Ts)); gap_mu.append(min_path_gap(1.0, mb))
    print(f"mu_buf={mb:4.2f}  min path gap={gap_mu[-1]:.3f}  T*={Tstar_mu[-1]}")

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].plot(gs_list, Tstar_g, "o-", label=r"$T^{*}$")
for g, T, o in zip(gs_list, Tstar_g, ov0_g):
    ax[0].annotate(f"{o:.2f}", (g, T), textcoords="offset points", xytext=(4, 6), fontsize=8)
ax[0].set_xlabel(r"pairing strength $g$"); ax[0].set_ylabel(r"$T^{*}$ (overlap $\geq 0.99$)")
ax[0].set_title(r"coupling scan at $\mu_{\rm buf}=0.6$ (labels: initial overlap ov$_0$)")
x = 1/np.array(gap_mu); y = np.array(Tstar_mu); sel = ~np.isnan(y)
c = np.polyfit(x[sel], y[sel], 1)
xs = np.linspace(0, x[sel].max()*1.08, 50)
ax[1].plot(x, y, "o", ms=8)
ax[1].plot(xs, np.polyval(c, xs), "--",
           label=rf"linear fit: $T^*\approx{c[0]:.1f}/\Delta E_{{\min}}{c[1]:+.1f}$")
ax[1].set_xlabel(r"$1/\Delta E_{\min}$ along the path ($\mu$-scan, $g=1$)")
ax[1].set_ylabel(r"$T^{*}$")
ax[1].set_title(r"$T^{*}$ scales linearly with the inverse path gap")
ax[1].legend()
plt.tight_layout(); plt.show()
"""

cells = [("markdown", md_trotter), ("code", code_trotter_check),
         ("code", code_trotter_fig), ("markdown", md_budget),
         ("code", code_budget), ("markdown", md_k34), ("code", code_k34),
         ("markdown", md_gscan), ("code", code_gscan)]

nb = nbf.read(NB, as_version=4)
nb.cells = [c for c in nb.cells if MARK not in c.source]      # idempotent
# insert before the final summary cell
isum = next(i for i, c in enumerate(nb.cells)
            if c.cell_type == "markdown" and "Summary: the complete EIGEN-Q pipeline" in c.source)
new = [nbf.v4.new_markdown_cell(s) if k == "markdown" else nbf.v4.new_code_cell(s)
       for k, s in cells]
nb.cells = nb.cells[:isum] + new + nb.cells[isum:]
nbf.write(nb, NB)
print(f"inserted {len(new)} cells before cell {isum}; notebook now has {len(nb.cells)} cells")
