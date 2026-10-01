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
  four to eight (the $1p$--$1h$ strings drag a parity chain across the register), so a
  Trotterised controlled evolution of $\hat H$ costs roughly five times more per step than
  for the pure pairing force; the pair-breaking strings do not fall into the two commuting
  groups of Section~13, and extending the circuit-level budget of Section~13 to $f\neq0$
  is the natural next step.
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
   the Hamiltonian ($61\to325$ Pauli strings at $k=4$) rather than in the pipeline.
"""
anchor = "   acceptance holding at $\\approx p$.\n"
if anchor in summ.source and "Beyond seniority (Section 16)" not in summ.source:
    summ.source = summ.source.replace(anchor, anchor + addition)
nbf.write(nb, NB)
print(f"inserted {len(new_cells)} cells before the summary (now cell {summary_idx + len(new_cells)}); "
      f"notebook has {len(nb.cells)} cells")
