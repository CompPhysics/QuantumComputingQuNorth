#!/usr/bin/env python3
"""Insert the particle-hole (pair-breaking) section into
ResolutionRefPairing.ipynb, after the coupling-strength scan and before the
summary cell, and renumber the summary.  Built with nbformat (never hand-edit
the .ipynb JSON); idempotent: existing cells tagged with the section marker
are replaced.  Run from the repo root:  python3 scripts/extend_notebook_ph.py
Then execute the new cells:            python3 scripts/run_new_cells.py --mark "<!-- EXT:PH -->"
"""
import pathlib
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
MARK = "<!-- EXT:PH -->"

md_intro = MARK + r"""
---
## 16&nbsp; Beyond seniority: the pairing model with a particle-hole interaction

Everything so far rests on a Hamiltonian that is almost *too* well behaved.  The
pairing force moves complete pairs and does nothing else, so the seniority (the number
of unpaired particles) is conserved, the $N=4$ ground state lives in a
$\binom{k}{2}$-dimensional seniority-zero space, and an entire class of excitations is
absent: there is no matrix element between the reference determinant and any
one-particle-one-hole state.  Realistic nuclear and electronic Hamiltonians do have that
channel, and mean-field theories such as Tamm--Dancoff and RPA are built on it.

Following Chapter 4 of *Quantum mechanics for many-particle systems* we therefore add
one extra term, keeping the level scheme, the particle number $N=4$ and the resolutions
$k=2,3,4$ used above unchanged:
$$
\hat H=\hat H_0+\hat V_{\rm pair}+\hat V_{\rm ph},\qquad
\hat V_{\rm ph}=-\frac{f}{2}\sum_{p,q,r}\Big(\hat a^{\dagger}_{p\uparrow}\hat a^{\dagger}_{p\downarrow}
\hat a_{q\downarrow}\hat a_{r\uparrow}+{\rm h.c.}\Big).
$$
The triple sum contains three distinct processes: the $q=r$ pieces are the pairing term
again and merely renormalise $g\to g+2f$; the $q=p\neq r$ pieces carry a *single*
particle from one level to another (a $1p$--$1h$ excitation of seniority two); and the
pieces with $p,q,r$ all different remove two unpaired particles and deposit a pair
elsewhere (seniority-two $2p$--$2h$ states).  $\hat V_{\rm ph}$ still conserves
$N_\uparrow$ and $N_\downarrow$ separately, so $S_z$ remains good, but the seniority does
not: at $k=4$ the $S_z=0$ sector has $36$ determinants instead of the six pair
configurations, and the one-qubit-per-pair encoding is no longer available --- the
$2k$-qubit Jordan--Wigner register used throughout this notebook is now the *only*
option.  Two further consequences matter for the pipeline: the reference determinant is
no longer the Hartree--Fock solution (Brillouin's theorem fails, with
$\langle\Phi_i^a|\hat H|\Phi_0\rangle=\pm f/2$), so the UCCSD *singles*, identically zero
for the pure pairing force, switch on; and the reference energy becomes
$E_{\rm ref}=2-g-2f$.

We use the strengths $f=0$, $0.05$, $0.2$ and $0.5$ of the book, at $g=1$, so that every
exact number below can be cross-checked against its tables; `pairinglib` implements
the term as `H_pairing_sparse(..., f=f)` / `build_H_full(..., f=f)` (module `phmodel`)
and the refinement routines accept the same keyword.  Nothing is redefined here.
"""

code_bench = "# " + MARK + r"""
import numpy as np, matplotlib.pyplot as plt
try:
    import pairinglib as pl
except ModuleNotFoundError:
    import sys, pathlib
    sys.path.insert(0, str((pathlib.Path.cwd()/".."/"src").resolve()))
    import pairinglib as pl

N, G = 4, 1.0
FS = [0.0, 0.05, 0.2, 0.5]                       # pair-breaking strengths of the book

# ---- (a) exact and variational references at k=3 and k=4 ----
print("Pairing + particle-hole model, N=4, g=1  (E_ref = 2-g-2f; FCI in the full N-sector)")
print("="*96)
print(f"{'k':>2} {'f':>5} {'E_ref':>8} {'E_FCI':>12} {'E_UCCSD':>12} {'E_UCCSD-E_FCI':>14} "
      f"{'max|singles|':>13} {'sen-0 weight':>13} {'dim(Sz=0)':>10}")
bench = {}
for k in (3, 4):
    su = pl.setup_uccsd(k, N)
    sz0 = sum(1 for I in su['states']                                   # N_up = N_dn determinants
              if sum((I >> (2*k-1-j)) & 1 for j in range(0, 2*k, 2)) == N//2)
    for f in FS:
        H = pl.H_pairing_sparse(k, G, N, su['states'], su['index'], f=f)
        Hd = H.toarray(); w, V = np.linalg.eigh(Hd); fci = w[0]
        E, x, _, _ = pl.uccsd_vqe(su, H)
        sing = np.max(np.abs(x[:su['nS']])) if su['nS'] else 0.0
        wt = pl.seniority_weight(V[:, 0], k, su['states'])
        bench[(k, f)] = dict(fci=fci, uccsd=E, Eref=pl.E_ref(N, G, f), gs=V[:, 0], sen0=wt, gap=w[1]-w[0])
        print(f"{k:>2} {f:>5.2f} {pl.E_ref(N, G, f):>8.4f} {fci:>12.6f} {E:>12.6f} {E-fci:>14.2e} "
              f"{sing:>13.4f} {wt:>13.4f} {sz0:>10d}")
print("\nBook anchors (k=4): f=0.05 -> 0.45058234, f=0.2 -> -0.18348455, f=0.5 -> -1.69173670")

# ---- (b) what connects to the reference determinant (book Table 'pairingphcoupling') ----
print("\nDeterminants reached from the reference by one action of H (k=4, N=4):")
for f in (0.0, 0.05):
    tab = pl.ph_coupling_table(4, N, G, f)
    desc = ", ".join(f"{cnt} x {r}p-{r}h (seniority {s}) |H|={v}" for (r, s), (cnt, v) in sorted(tab.items()))
    print(f"  f={f:.2f}: {desc}")

# ---- (c) the price on a quantum computer: Pauli strings of the 2k-qubit Hamiltonian ----
print("\nPauli-string decomposition of H on 2k qubits (count, largest weight):")
print(f"{'k':>2} {'qubits':>7} {'pairing (f=0)':>16} {'pairing + p-h':>16}")
for k in (2, 3, 4):
    c0 = pl.pauli_count(pl.build_H_full(k, G)); c1 = pl.pauli_count(pl.build_H_full(k, G, f=0.2))
    print(f"{k:>2} {2*k:>7} {str(c0):>16} {str(c1):>16}")
"""

code_fig_scan = "# " + MARK + r"""
# ---- continuous f-scan at g=1: energies, pair breaking and the UCCSD error ----
f_grid = np.linspace(0.0, 0.6, 25)
curves = {k: dict(fci=[], uccsd=[], sen0=[]) for k in (3, 4)}
for k in (3, 4):
    su = pl.setup_uccsd(k, N); x0 = None
    for f in f_grid:
        H = pl.H_pairing_sparse(k, G, N, su['states'], su['index'], f=f)
        w, V = np.linalg.eigh(H.toarray())
        E, x0, _, _ = pl.uccsd_vqe(su, H, x0=x0)
        curves[k]['fci'].append(w[0]); curves[k]['uccsd'].append(E)
        curves[k]['sen0'].append(pl.seniority_weight(V[:, 0], k, su['states']))
Eref_grid = np.array([pl.E_ref(N, G, f) for f in f_grid])

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
axL.plot(f_grid, Eref_grid, 'k--', lw=1.4, label=r'reference determinant $E_{\rm ref}=1-2f$')
for k, c in [(3, 'C2'), (4, 'C3')]:
    axL.plot(f_grid, curves[k]['fci'], c, lw=1.8, label=f'FCI, $k={k}$')
    axL.plot(f_grid, curves[k]['uccsd'], c, ls=':', lw=1.8, label=f'UCCSD, $k={k}$')
for f in FS[1:]:
    axL.axvline(f, color='gray', lw=0.6, ls=':')
axL.set_xlabel('pair-breaking strength $f$'); axL.set_ylabel('ground-state energy')
axL.set_title('$N=4$, $g=1$: exact and UCCSD energies vs $f$'); axL.legend(fontsize=8); axL.grid(alpha=0.3)
for k, c in [(3, 'C2'), (4, 'C3')]:
    axR.semilogy(f_grid[1:], 1-np.array(curves[k]['sen0'][1:]), c, lw=1.8,
                 label=f'GS weight outside seniority zero, $k={k}$')
axR.semilogy(f_grid[1:], np.array(curves[4]['uccsd'][1:])-np.array(curves[4]['fci'][1:]),
             'C3', ls=':', lw=1.8, label='$E_{\\rm UCCSD}-E_{\\rm FCI}$, $k=4$  ($k=3$: exact)')
axR.set_ylim(1e-4, 0.3)
axR.set_xlabel('pair-breaking strength $f$'); axR.set_ylabel('weight  /  energy error')
axR.set_title('pair breaking and the variational error'); axR.legend(fontsize=8); axR.grid(alpha=0.3, which='both')
plt.tight_layout(); plt.show()
print("At k=3 UCCSD remains exact for every f (as it is at f=0, now with non-zero singles);")
print("at k=4 the UCCSD error grows with f, from %.1e at f=0 to %.1e at f=0.5, as the pair-breaking"
      % (bench[(4, 0.0)]['uccsd']-bench[(4, 0.0)]['fci'], bench[(4, 0.5)]['uccsd']-bench[(4, 0.5)]['fci']))
print("channels take up to %.1f%% of the ground state out of the seniority-zero space." % (100*(1-bench[(4, 0.5)]['sen0'])))
"""

md_refine = MARK + r"""
### Resolution refinement with pairs broken

The pipeline does not care where the correlations come from: we prolong the $k=2$ state
(at $N=4$ still a single determinant, now with energy $E_{\rm low}=E_{\rm ref}$) into the
$k=4$ space and evolve along the same path
$\cos^2(\pi s/2)\,P(\hat H_{\rm low}-\mu)P^\dagger+\sin^2(\pi s/2)(\hat H_{\rm high}-\mu)$
with the same shift rule $\mu=E_{\rm low}+0.6$, for each $f$.  Two things change with
$f$.  The target state acquires seniority-two components that the prolonged determinant
does not contain, so the initial overlap ${\rm ov}_0$ drops; and the effective pairing
strength $g+2f$ grows, which is the mechanism already seen in the $g$-scan of Section~15.
The minimum gap along the path, however, is still set by the shift, so the adiabatic
time $T^{*}$ needed for an overlap of $0.99$ grows only because the state has further to
travel, not because the path closes.  We also repeat the $k=3\to4$ chain from the
gate-level UCCSD state at $k=3$, whose singles now carry non-zero angles.
"""

code_refine = "# " + MARK + r"""
# ---- k=2->4 refinement for the four pair-breaking strengths ----
def min_path_gap_f(f, mu_buf=0.6, k_low=2, k_high=4):
    nl, sl, il = pl.build_sector(k_low, N)
    Hl = pl.H_pairing_sparse(k_low, G, N, sl, il, f=f).toarray(); El = np.linalg.eigvalsh(Hl)[0]
    Pm = pl.prolong_matrix(k_low, k_high, N)
    nh, sh, ih = pl.build_sector(k_high, N)
    Hh = pl.H_pairing_sparse(k_high, G, N, sh, ih, f=f).toarray()
    mu = El + mu_buf
    Hle = Pm @ (Hl - mu*np.eye(len(Hl))) @ Pm.T; Hhs = Hh - mu*np.eye(len(Hh))
    return min(np.diff(np.linalg.eigvalsh(np.cos(np.pi*s/2)**2*Hle + np.sin(np.pi*s/2)**2*Hhs)[:2])[0]
               for s in np.linspace(0, 1, 41))

def Tstar_of(r, Ts, target=0.99):
    hit = np.where(r["ov"] >= target)[0]
    return Ts[hit[0]] if hit.size else np.nan

Ts = np.array([0.5, 1, 2, 3, 4, 6, 8, 12, 16, 20, 26, 32, 40, 50, 65, 80])
ref_f, ref34_f = {}, {}
print("Refinement to k=4 (N=4, g=1, mu = E_low + 0.6), overlap target 0.99")
print("="*92)
print(f"{'f':>5} {'E_low':>7} {'E_FCI(k=4)':>11} {'gap(H_high)':>12} {'min path gap':>13} "
      f"{'ov0 2->4':>9} {'T* 2->4':>8} {'ov0 3->4':>9} {'T* 3->4':>8}")
for f in FS:
    r = pl.refine_NK(N, 4, Ts, f=f); ref_f[f] = r
    E3, psi3_full = pl.gate_uccsd_state(3, N, G, f=f)
    nq3, st3, ix3 = pl.build_sector(3, N)
    psi3 = np.array([psi3_full[s] for s in st3]); psi3 /= np.linalg.norm(psi3)
    r34 = pl.refine_from_state(psi3, N, 3, 4, Ts, f=f); ref34_f[f] = r34
    print(f"{f:>5.2f} {r['El']:>7.4f} {r['Eh']:>11.6f} {r['gap']:>12.4f} {min_path_gap_f(f):>13.4f} "
          f"{r['ov0']:>9.4f} {Tstar_of(r, Ts):>8.0f} {r34['ov0']:>9.4f} {Tstar_of(r34, Ts):>8.0f}")
print("\nCompare the g-scan of Section 15: the pure pairing model at g=2 (= g+2f for f=0.5) had ov0=0.6167, T*=32.")

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
cols = {0.0: 'C0', 0.05: 'C1', 0.2: 'C2', 0.5: 'C3'}
for f in FS:
    r = ref_f[f]
    axL.plot(Ts, r['ov'], 'o-', color=cols[f], lw=1.6, ms=4, label=f'$f={f}$  (guess {r["ov0"]:.3f})')
    axR.plot(Ts, r['en'] - r['Eh'], 'o-', color=cols[f], lw=1.6, ms=4, label=f'$f={f}$')
axL.axhline(0.99, color='gray', ls=':', lw=1); axL.axhline(1.0, color='k', lw=0.6)
axL.set_xlabel('total adiabatic time $T$'); axL.set_ylabel(r'overlap $|\langle\Psi_{\rm high}|\Phi(T)\rangle|^2$')
axL.set_title(r'$k=2\to4$ refinement with the particle-hole term'); axL.legend(fontsize=8); axL.grid(alpha=0.3)
axR.set_yscale('log'); axR.set_xlabel('total adiabatic time $T$')
axR.set_ylabel(r'$\langle\Phi(T)|H_{\rm high}|\Phi(T)\rangle - E_{\rm FCI}$')
axR.set_title('energy error (variational, falls toward FCI from above)'); axR.legend(fontsize=8); axR.grid(alpha=0.3, which='both')
plt.tight_layout(); plt.show()
"""

md_rodeo = MARK + r"""
### Rodeo filtering with a denser spectrum

The rodeo algorithm is the stage most exposed to the extra term: its energy scan sees
the *whole* spectrum of the $2k$-qubit Hamiltonian in the $N=4$ sector, which the
particle-hole term has made both denser and lower (the $36$-dimensional $S_z=0$ block
now mixes with the seniority-two states).  We repeat the two rodeo experiments of
Section~11 at the strongest coupling, $f=0.5$, where the gap to the first excited
state is in fact *larger* than in the pure pairing model ($2.86$ against $1.82$), and
then track fidelity and acceptance for all four strengths, always feeding the filter the
resolution-refined ($T=25$) state.
"""

code_rodeo_scan = "# " + MARK + r"""
# ---- rodeo energy scan at f=0.5: prolonged vs refined input ----
F_SCAN = 0.5
psi_ref, Hrod = pl.refine_state(N, 4, 25.0, f=F_SCAN)
Pm = pl.prolong_matrix(2, 4, N)
nl, sl, il = pl.build_sector(2, N); Hl = pl.H_pairing_sparse(2, G, N, sl, il, f=F_SCAN).toarray()
psi_pro = Pm @ np.linalg.eigh(Hl)[1][:, 0]; psi_pro = psi_pro/np.linalg.norm(psi_pro)
wR, VR = np.linalg.eigh(Hrod); E0 = wR[0]; GS = VR[:, 0]
p_ref = abs(np.vdot(GS, psi_ref))**2; p_pro = abs(np.vdot(GS, psi_pro))**2
print(f"f={F_SCAN}: E0 = {E0:.6f} (book: -1.69173670), gap = {wR[1]-wR[0]:.4f}, "
      f"input overlap p: prolonged {p_pro:.4f}, refined (T=25) {p_ref:.4f}")

def scan(v0, Es, sigma, H, M=4, n_sets=12, seed=0):
    rng = np.random.default_rng(seed); out = np.zeros(len(Es))
    for s in range(n_sets):
        for j, E in enumerate(Es):
            out[j] += pl.rodeo_allzero_prob(v0, H, rng.normal(0, sigma, M), E)
    return out/n_sets

Es_coarse = np.linspace(E0-0.6, 4.0, 90)
Es_fine = np.linspace(E0-0.6, E0+0.6, 161)
Pc_pro = scan(psi_pro, Es_coarse, 2.0, Hrod); Pc_ref = scan(psi_ref, Es_coarse, 2.0, Hrod)
Pf_ref = scan(psi_ref, Es_fine, 8.0, Hrod)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
axL.plot(Es_coarse, Pc_pro, 'o-', ms=3, color='C0', label='input: prolonged $\\Phi_0$ ($p$=%.3f)' % p_pro)
axL.plot(Es_coarse, Pc_ref, 's-', ms=3, color='C3', label='input: refined state ($p$=%.3f)' % p_ref)
for e in np.unique(np.round(wR[wR < 4.2], 4)):
    axL.axvline(e, color='gray', ls=':', lw=0.6)
axL.axhline(0.5**4, color='k', ls='--', lw=0.8, label='$2^{-M}$ background')
axL.set_xlabel('target energy $E$'); axL.set_ylabel(r'$P_{0^M}(E)$  (M=4, $\sigma=2$)')
axL.set_title(f'coarse scan, $f={F_SCAN}$ (dotted = exact eigenvalues)'); axL.legend(fontsize=8); axL.grid(alpha=0.3)
axR.plot(Es_fine, Pf_ref, 's-', ms=3, color='C3')
axR.axvline(E0, color='gray', ls=':', lw=1.0, label=f'$E_0={E0:.4f}$')
axR.set_xlabel('target energy $E$'); axR.set_ylabel(r'$P_{0^M}(E)$  (M=4, $\sigma=8$)')
axR.set_title('fine scan on the ground-state peak'); axR.legend(fontsize=8); axR.grid(alpha=0.3)
plt.tight_layout(); plt.show()
jpk = np.argmax(Pf_ref); y0, y1, y2 = Pf_ref[jpk-1:jpk+2]; dx = Es_fine[1]-Es_fine[0]
Epk = Es_fine[jpk] + 0.5*dx*(y0-y2)/(y0-2*y1+y2)
print(f"peak of fine scan at E = {Epk:.4f}  (exact E0 = {E0:.6f});  "
      f"{np.sum(wR < 4.0)} eigenvalues below E=4 against {np.sum(np.linalg.eigvalsh(pl.H_pairing_sparse(4, G, N, *pl.build_sector(4, N)[1:]).toarray()) < 4.0)} for f=0")
"""

code_rodeo_track = "# " + MARK + r"""
# ---- eigenstate preparation vs number of cycles, for the four strengths ----
M_max, n_sets, sigma = 10, 40, 3.0
Ms = np.arange(1, M_max+1)
track = {}
for f in FS:
    psi_in, Hh = pl.refine_state(N, 4, 25.0, f=f)
    w, V = np.linalg.eigh(Hh); E0f = w[0]; GSf = V[:, 0]
    p0 = abs(np.vdot(GSf, psi_in))**2
    rng = np.random.default_rng(7); F = np.zeros(M_max); A = np.zeros(M_max)
    for s in range(n_sets):
        ff, aa = pl.rodeo_track(psi_in, Hh, rng.normal(0, sigma, M_max), E0f, GSf); F += ff; A += aa
    track[f] = (F/n_sets, A/n_sets, p0)

fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))
for f in FS:
    F, A, p0 = track[f]
    axL.semilogy(Ms, 1-F, 'o-', color=cols[f], label=f'$f={f}$  ($p$={p0:.4f})')
    axL.semilogy(Ms, (1-p0)/p0*2.0**(-Ms), ':', color=cols[f], lw=1.1)
    axR.plot(Ms, A, 'o-', color=cols[f], label=f'$f={f}$'); axR.axhline(p0, color=cols[f], ls=':', lw=1.0)
axL.set_xlabel('number of rodeo cycles $M$'); axL.set_ylabel(r'infidelity $1-\mathcal{F}_M$')
axL.set_title(r'refined ($T=25$) input; dotted: $\frac{1-p}{p}2^{-M}$'); axL.legend(fontsize=8); axL.grid(alpha=0.3, which='both')
axR.set_xlabel('number of rodeo cycles $M$'); axR.set_ylabel('acceptance $P_{0^M}$')
axR.set_title(r'post-selection acceptance $\to p$ (dotted)'); axR.legend(fontsize=8); axR.grid(alpha=0.3)
plt.tight_layout(); plt.show()
print("Final fidelities (M = %d, averaged over %d time sets, sigma = %.1f):" % (M_max, n_sets, sigma))
for f in FS:
    F, A, p0 = track[f]
    print(f"  f={f:<5}: input p={p0:.5f}   F={F[-1]:.6f}   1-F={1-F[-1]:.1e}   acceptance={A[-1]:.4f}")
"""

md_trot = MARK + r"""
### The circuit-level budget with the particle-hole term

Section~13 Trotterised the pipeline for the pairing force with two kinds of exactly
exponentiable blocks: the diagonal term and one pair-hopping term per level pair.  The
particle-hole term adds blocks of its own.  Its $q=r$ pieces are pairing terms and only
rescale the existing rotation angles ($g\to g+2f$) at no cost; the remaining pieces,
one per $(p,q,r)$ with $q\neq r$, split into $2k(k-1)$ one-particle-one-hole blocks
($\hat a^\dagger_{p\uparrow}\hat n_{p\downarrow}\hat a_{r\uparrow}+{\rm h.c.}$ and its
spin partner; four Pauli strings each) and $k(k-1)(k-2)$ pair-breaking blocks
(eight strings each, like a double excitation).  Within every block the strings commute,
so each factor $e^{-i\hat T\delta t}$ is again realised exactly by its rotations, and the
product-formula error comes only from the non-commutativity *between* blocks --- $55$ of
them at $k=4$ instead of seven.  What changes on hardware is the depth of one step: the
$1p$--$1h$ strings carry Jordan--Wigner parity chains of weight up to $2k$, and the
step cost grows from $296$ to $2728$ CNOTs at $k=4$ (`pl.cnot_counts(4, f=f)`).  We use
`pl.model_terms`, `pl.trotter_error`, `pl.refine_state_trotter` and `pl.trotter_U_dt`
with the `f` keyword; nothing is redefined here.
"""

code_trot_terms = "# " + MARK + r"""
# ---- term decomposition with V_ph, Trotter error, and the per-step CNOT bill ----
nq4, st4, ix4 = pl.build_sector(4, N)
print("Term decomposition and product-formula error at k=4, N=4, g=1, t=1")
print("="*78)
print(f"{'f':>5} {'blocks':>7} {'‖Σ T − H‖':>11} | {'n=1 (o1)':>9} {'n=4 (o1)':>9} {'n=16 (o1)':>10} | "
      f"{'n=1 (o2)':>9} {'n=4 (o2)':>9} {'n=16 (o2)':>10}")
terms_f = {}
for f in (0.0, 0.2, 0.5):
    H = pl.H_pairing_sparse(4, G, N, st4, ix4, f=f).toarray()
    tm = pl.model_terms(4, G, N, f); terms_f[f] = tm
    e1 = [pl.trotter_error(4, G, N, 1.0, n, order=1, f=f, terms=tm) for n in (1, 4, 16)]
    e2 = [pl.trotter_error(4, G, N, 1.0, n, order=2, f=f, terms=tm) for n in (1, 4, 16)]
    print(f"{f:>5.2f} {len(tm):>7d} {np.abs(sum(tm)-H).max():>11.1e} | "
          + " ".join(f"{x:>9.2e}" for x in e1) + " | " + " ".join(f"{x:>9.2e}" for x in e2))

print("\nCNOTs per first-order Trotter step (ladders 2(w−1) per string; controlled: +2 per rotation)")
print(f"{'k':>2} {'qubits':>7} {'pairing: rot / CNOT / ctrl':>28} {'+ p-h: rot / CNOT / ctrl':>26} {'ratio':>6}")
for k in (3, 4):
    c0 = pl.cnot_counts(k); c1 = pl.cnot_counts(k, f=0.2)
    print(f"{k:>2} {2*k:>7} {c0['rotations']:>10} / {c0['step_uncontrolled']:>5} / {c0['step_controlled']:>5}"
          f" {c1['rotations']:>10} / {c1['step_uncontrolled']:>5} / {c1['step_controlled']:>5}"
          f" {c1['step_uncontrolled']/c0['step_uncontrolled']:>6.1f}")
pc = pl.ph_pauli_costs(4)
for fam, d in pc.items():
    print(f"  k=4 {fam} blocks: {d['terms']} terms, {d['rotations']} strings, {d['cnots']} CNOTs, largest weight {d['max_weight']}")
"""

code_trot_fig = "# " + MARK + r"""
# ---- Trotterised refinement and rodeo with the particle-hole term ----
T_AD = 30.0; ns_grid = [10, 20, 40, 80, 160]
inf_ref = {}; inf_exact = {}; GSf = {}; E0f = {}; Hf = {}
for f in (0.2, 0.5):
    Hh = pl.H_pairing_sparse(4, G, N, st4, ix4, f=f).toarray(); w, V = np.linalg.eigh(Hh)
    Hf[f] = Hh; E0f[f] = w[0]; GSf[f] = V[:, 0]
    for order in (1, 2):
        inf_ref[(f, order)] = [1 - abs(np.vdot(GSf[f], pl.refine_state_trotter(N, 4, T_AD, ns, order=order, f=f)[0]))**2
                               for ns in ns_grid]
    pe, _ = pl.refine_state(N, 4, T_AD, f=f); inf_exact[f] = 1 - abs(np.vdot(GSf[f], pe))**2

# rodeo with Trotterised controlled evolution at fixed step size dt, refined (n_s=80, order 2) input
NS_OP, M_CYC, SIG, NSAMP = 80, 6, 4.0, 12
rng = np.random.default_rng(11)
tsets = [np.abs(rng.normal(0.0, SIG, M_CYC)) for _ in range(NSAMP)]
psi_in = {}; p_in = {}; curves = {}; accepts = {}; mean_steps = {}
for f in (0.2, 0.5):
    psi_in[f], _ = pl.refine_state_trotter(N, 4, T_AD, NS_OP, order=2, f=f)
    p_in[f] = abs(np.vdot(GSf[f], psi_in[f]))**2
    w, V = np.linalg.eigh(Hf[f])
    for dtm in (1.0, 0.5, 0.25, 0.125, None):
        F_acc = np.zeros(M_CYC); A_fin = 0.0; nst = 0
        for ts in tsets:
            if dtm is None:
                Ufun = lambda t: (V*np.exp(-1j*w*t)) @ V.T
            else:
                Ufun = lambda t, dtm=dtm, f=f: pl.trotter_U_dt(4, G, N, t, dtm, order=2, terms=terms_f[f])[0]
                nst += sum(max(1, int(np.ceil(t/dtm))) for t in ts)
            Fc, Ac = pl.rodeo_track_U(psi_in[f], Ufun, ts, E0f[f], GSf[f])
            F_acc += Fc/NSAMP; A_fin += Ac[-1]/NSAMP
        curves[(f, dtm)] = F_acc; accepts[(f, dtm)] = A_fin; mean_steps[(f, dtm)] = nst/NSAMP

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for f, ls in ((0.2, '-'), (0.5, '--')):
    for order, c in ((1, 'C0'), (2, 'C3')):
        ax[0].loglog(ns_grid, inf_ref[(f, order)], 'o' + ls, color=c, label=f'$f={f}$, order {order}')
    ax[0].axhline(inf_exact[f], color='gray', ls=ls, lw=0.9, label=f'$f={f}$, exact slices')
ax[0].set_xlabel('refinement steps $n_s$'); ax[0].set_ylabel(r'$1-|\langle\Psi_0|\Phi\rangle|^2$')
ax[0].set_title(r'Trotterised refinement $k=2\to4$, $T=30$, with $V_{\rm ph}$'); ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3, which='both')
Ms = np.arange(1, M_CYC+1)
for dtm, c in ((1.0, 'C0'), (0.5, 'C1'), (0.25, 'C2'), (0.125, 'C4'), (None, 'k')):
    lab = 'exact evolution' if dtm is None else f'$\\delta t={dtm}$ ({mean_steps[(0.5, dtm)]/M_CYC:.1f} steps/cycle)'
    ax[1].semilogy(Ms, 1 - curves[(0.5, dtm)], 'o-' if dtm else ':', color=c, label=lab + f', acc {accepts[(0.5, dtm)]:.3f}')
ax[1].set_xlabel('rodeo cycles $M$'); ax[1].set_ylabel(r'$1-\mathcal{F}_M$')
ax[1].set_title(f'rodeo on the circuit-level refined input, $f=0.5$ ($p={p_in[0.5]:.4f}$)')
ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3, which='both')
plt.tight_layout(); plt.show()
for f in (0.2, 0.5):
    print(f"f={f}: circuit-level refined input (n_s=80, order 2): overlap p = {p_in[f]:.5f}  (exact slices: {1-inf_exact[f]:.5f})")
    for dtm in (1.0, 0.5, 0.25, 0.125, None):
        lab = 'exact evolution' if dtm is None else f'dt={dtm}'
        print(f"   {lab:>16}: infidelity after 2 cycles = {1-curves[(f, dtm)][1]:.3e}, after 6 = {1-curves[(f, dtm)][5]:.3e}, acceptance = {accepts[(f, dtm)]:.4f}")
"""

code_trot_budget = "# " + MARK + r"""
# ---- depth budget at k=4 with the particle-hole term (per accepted preparation) ----
ORDER_FAC, M_OP = 2, 2
print("k=4 depth budget with V_ph (CNOTs), N=4, g=1; refinement n_s=80 order 2, rodeo M=2")
print("="*96)
print(f"{'f':>5} {'dt':>6} {'step':>6} {'ctrl step':>9} {'refine':>8} {'rodeo':>8} {'per run':>8} {'P_M':>7} "
      f"{'per accepted':>12} {'1-F (M=2)':>10} {'UCCSD floor':>11}")
budget = {}
for f in (0.2, 0.5):
    cn = pl.cnot_counts(4, f=f); err_vqe = bench[(4, f)]['uccsd'] - bench[(4, f)]['fci']
    for dtm in (0.25, 0.125):
        rng = np.random.default_rng(3); F2 = A2 = 0.0; nsteps_tot = 0
        for _ in range(NSAMP):
            ts = np.abs(rng.normal(0.0, SIG, M_OP))
            Fc, Ac = pl.rodeo_track_U(psi_in[f], lambda t: pl.trotter_U_dt(4, G, N, t, dtm, order=2, terms=terms_f[f])[0],
                                      ts, E0f[f], GSf[f])
            F2 += Fc[-1]/NSAMP; A2 += Ac[-1]/NSAMP
            nsteps_tot += sum(max(1, int(np.ceil(t/dtm))) for t in ts)
        nsteps_avg = nsteps_tot/NSAMP
        cn_refine = NS_OP*cn['step_uncontrolled']*ORDER_FAC
        cn_rodeo = int(round(nsteps_avg*cn['step_controlled']*ORDER_FAC))
        cn_run = cn_refine + cn_rodeo; cn_acc = cn_run/A2
        budget[(f, dtm)] = dict(refine=cn_refine, rodeo=cn_rodeo, run=cn_run, acc=A2, per_acc=cn_acc, inf=1-F2, steps=nsteps_avg)
        print(f"{f:>5.2f} {dtm:>6} {cn['step_uncontrolled']:>6} {cn['step_controlled']:>9} {cn_refine:>8} {cn_rodeo:>8} "
              f"{cn_run:>8} {A2:>7.4f} {cn_acc:>12.0f} {1-F2:>10.2e} {err_vqe:>11.1e}")
print("Pairing-only reference (Section 13): 296 / 416 CNOTs per step, 47360 + ~23000 ≈ 7.1e4 per accepted prep, 1-F = 5.5e-4.")
print("The UCCSD-VQE circuit (1312 CNOTs, 26 parameters) is unchanged by V_ph; only its floor moves.")
"""

md_close = MARK + r"""
### What the particle-hole term teaches us

* **The pipeline is unchanged; only the Hamiltonian is.**  Prolongation, the adiabatic
  path, the shift rule and the rodeo filter were all written for the $2k$-qubit
  Jordan--Wigner register, and they run on the seniority-breaking Hamiltonian without
  modification.  For every $f$ the refined state reaches an overlap above $0.99$, the
  rodeo infidelity falls geometrically and the acceptance stays at the input overlap.
* **The cost moves, it does not explode.**  The minimum path gap is still the shift
  $\mu-E_{\rm low}=0.6$; $T^{*}$ grows from $16$ to $32$ between $f=0$ and $f=0.5$ because
  the prolonged determinant starts further from the target (${\rm ov}_0$ from $0.87$ to
  $0.60$), the same mechanism as raising $g$ to $g+2f$ in Section~15.  A better coarse
  stage buys this back: the gate-level $k=3$ UCCSD state, singles now active, starts at
  ${\rm ov}_0=0.85$ even at $f=0.5$.
* **The variational route degrades, the filtered route does not.**  At $k=4$ the UCCSD
  error grows by a factor $\sim20$ across the scan as the ground state leaks out of the
  seniority-zero space, whereas the refined-plus-rodeo state is limited only by the
  number of cycles.
* **The price on hardware is in the Hamiltonian, not in the state preparation.**  The
  number of Pauli strings grows from $61$ to $325$ at $k=4$ and their largest weight from
  four to eight (the $1p$--$1h$ strings drag a parity chain across the register).  Compiled
  block by block, one Trotter step costs $2728$ CNOTs instead of $296$ --- a factor
  $9.2$ --- while the product-formula error at fixed step count grows only by a factor
  two to four, because the new blocks carry the small coefficient $f/2$.  The same
  operating point ($n_s=80$, $M=2$, $\delta t=0.25$) therefore still works at $f=0.2$ and
  the budget per accepted preparation scales essentially with the step cost; at $f=0.5$
  the rodeo needs $\delta t=0.125$ to keep its acceptance, doubling the rodeo share.
"""

nb = nbf.read(NB, as_version=4)
# keep the stored outputs of tagged code cells whose source is unchanged
kept = {c.source: (c.get("outputs", []), c.get("execution_count"))
        for c in nb.cells if MARK in c.source and c.cell_type == "code"}
nb.cells = [c for c in nb.cells if MARK not in c.source]        # idempotent
summary_idx = next(i for i, c in enumerate(nb.cells)
                   if c.cell_type == "markdown" and "Summary: the complete EIGEN-Q pipeline" in c.source)
new_cells = [nbf.v4.new_markdown_cell(md_intro), nbf.v4.new_code_cell(code_bench),
             nbf.v4.new_code_cell(code_fig_scan), nbf.v4.new_markdown_cell(md_refine),
             nbf.v4.new_code_cell(code_refine), nbf.v4.new_markdown_cell(md_rodeo),
             nbf.v4.new_code_cell(code_rodeo_scan), nbf.v4.new_code_cell(code_rodeo_track),
             nbf.v4.new_markdown_cell(md_trot), nbf.v4.new_code_cell(code_trot_terms),
             nbf.v4.new_code_cell(code_trot_fig), nbf.v4.new_code_cell(code_trot_budget),
             nbf.v4.new_markdown_cell(md_close)]
for c in new_cells:
    if c.cell_type == "code" and c.source in kept:
        c["outputs"], c["execution_count"] = kept[c.source]
nb.cells[summary_idx:summary_idx] = new_cells

# renumber the summary and record the new section in it
summ = nb.cells[summary_idx + len(new_cells)]
summ.source = summ.source.replace("## 16&nbsp; Summary", "## 17&nbsp; Summary")
addition = r"""
5. **Beyond seniority (Section 16)** — adding the pair-breaking particle-hole interaction
   $\hat V_{\rm ph}$ of the many-body book (strengths $f=0.05,0.2,0.5$ at $g=1$) opens the
   $1p$--$1h$ channel, switches on the UCCSD singles and takes the ground state out of the
   seniority-zero space.  The same four stages run unchanged: the refinement time grows
   only through the lower initial overlap (the path gap is still set by the shift), the
   rodeo filter converges geometrically for every $f$, and the hardware cost appears in
   the Hamiltonian ($61\to325$ Pauli strings, $296\to2728$ CNOTs per Trotter step at
   $k=4$) rather than in the pipeline, whose circuit-level operating point survives.
"""
anchor = "   acceptance holding at $\\approx p$.\n"
if anchor in summ.source and "Beyond seniority (Section 16)" not in summ.source:
    summ.source = summ.source.replace(anchor, anchor + addition)
nbf.write(nb, NB)
print(f"inserted {len(new_cells)} cells before the summary (now cell {summary_idx + len(new_cells)}); "
      f"notebook has {len(nb.cells)} cells")
