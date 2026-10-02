#!/usr/bin/env python3
"""(a) The pipeline at other model-space sizes and particle numbers (exact-slice
refinement, exact controlled evolution) with the matched UCCSD comparator;
(b) a nuclear observable -- level occupation numbers n_p and the correlation
energy -- for the exact, UCCSD, prolonged, refined and refined+rodeo states at
k=4, N=4, g=1 (circuit-level pipeline of Table II).
Usage: PYTHONPATH=src python3 scripts/sizes_observables.py"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src")); sys.path.insert(0, "src")
import pairinglib as pl
from pairinglib._pairlib import apply_exc
from pairinglib.hamiltonian import _bit

G = 1.0

def uccsd_sector_state(k, N, H):
    su = pl.setup_uccsd(k, N)
    E, x, _, _ = pl.uccsd_vqe(su, H)
    psi = np.zeros(su['M']); psi[su['hf']] = 1.0
    for kk in range(su['P']): apply_exc(psi, x[kk], su['arrs'][kk])
    return E, psi, su['P']

def rodeo_ensemble(v0, w, V, E0, GS, M, sigma=3.0, n_sets=40, seed=7):
    """Acceptance-weighted ensemble infidelity and mean acceptance after M cycles
    (exact controlled evolution, applied in the eigenbasis of H)."""
    rng = np.random.default_rng(seed); FA = A = 0.0
    c0 = V.T @ v0                                    # amplitudes in the eigenbasis
    for s in range(n_sets):
        ts = rng.normal(0, sigma, M); c = c0.astype(complex); P = 1.0
        for t in ts:
            cw = 0.5*(1 + np.exp(1j*(E0 - w)*t))*c; p = float(np.vdot(cw, cw).real)
            c = cw/np.sqrt(p); P *= p
        FA += abs(c[0])**2 * P; A += P
    return 1 - FA/A, A/n_sets

def uccsd_ncnot(k, N):
    """CNOT count of the gate-level UCCSD circuit from the Jordan-Wigner weights:
    a single (i,a) is 2 strings of weight a-i+1, a double (i<j<a<b) is 8 strings
    of weight 4+(j-i-1)+(b-a-1); a weight-w string costs 2(w-1) CNOTs
    (reproduces 272 at k=3 and 1312 at k=4)."""
    s, d = pl.uccsd_pool(k, N)
    return sum(4*(a-i) for (i, a) in s) + sum(16*(4+(j-i-1)+(b-a-1)-1) for (i, j, a, b) in d)

print("(a) pipeline at other sizes, g=1: coarse space = filled determinant (k_low=N/2), "
      "refinement T=30 (exact slices), rodeo sigma=3, 40 time sets")
print(f"{'(k,N)':>7} {'qubits':>6} {'dim':>5} {'gap':>6} {'p0 prol.':>9} {'p refined':>10} "
      f"{'1-F UCCSD':>10} {'dE UCCSD':>9} {'#par':>5} {'#CNOT':>6} | {'refined+rodeo M=2 / M=6':>24} | {'UCCSD+rodeo M=2 / M=6':>22} | {'step CNOT':>9}")
for (k, N) in ((4, 4), (6, 4), (4, 6), (6, 6)):
    nq, st, ix = pl.build_sector(k, N)
    Hs = pl.H_pairing_sparse(k, G, N, st, ix); H = Hs.toarray()
    w, V = np.linalg.eigh(H); E0, GS, gap = w[0], V[:, 0], w[1]-w[0]
    r = pl.refine_NK(N, k, [30.0], k_low=N//2)
    psi_ref, _ = pl.refine_state(N, k, 30.0, k_low=N//2)
    E_u, psi_u, npar = uccsd_sector_state(k, N, Hs)
    ncnot = uccsd_ncnot(k, N)
    rr = [rodeo_ensemble(psi_ref, w, V, E0, GS, M)[0] for M in (2, 6)]
    ru = [rodeo_ensemble(psi_u, w, V, E0, GS, M)[0] for M in (2, 6)]
    step = pl.cnot_counts(k)['step_uncontrolled']
    print(f"({k},{N}) {2*k:>6} {len(st):>5} {gap:6.3f} {r['ov0']:9.4f} {r['ov'][0]:10.5f} "
          f"{1-abs(np.vdot(GS, psi_u))**2:10.2e} {E_u-E0:9.2e} {npar:5d} {ncnot:6d} | "
          f"{rr[0]:.1e} / {rr[1]:.1e}{'':>8} | {ru[0]:.1e} / {ru[1]:.1e}{'':>6} | {step:9d}")

# ---- (b) observables at k=4, N=4 ------------------------------------------
print("\n(b) level occupations n_p = sum_sigma <n_{p sigma}> and correlation energy at k=4, N=4, g=1")
k, N = 4, 4
nq, st, ix = pl.build_sector(k, N)
H = pl.H_pairing_sparse(k, G, N, st, ix).toarray(); w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]
occ = np.array([[sum(_bit(I, 2*p, nq) + _bit(I, 2*p+1, nq) for I in [I]) for p in range(k)] for I in st], float)
def nbar(v): return (abs(v)**2) @ occ
def energy(v): return float(np.real(np.vdot(v, H @ v)))
E_ref = 1.0
Pm = pl.prolong_matrix(2, k, N); pro = Pm @ np.array([1.0]); pro /= np.linalg.norm(pro)
E_u, psi_full = pl.gate_uccsd_state(k, N, G); ucc = np.array([psi_full[I] for I in st]); ucc /= np.linalg.norm(ucc)
ref30, _ = pl.refine_state_trotter(N, k, 30.0, 80, order=2)
terms = pl.model_terms(k, G, N)
Ufun = lambda t: pl.trotter_U_dt(k, G, N, t, 0.25, order=2, terms=terms)[0]
rng = np.random.default_rng(3); tsets = [np.abs(rng.normal(0, 4.0, 6)) for _ in range(40)]
def rodeo_avg_obs(v0, M):
    nb = np.zeros(k); en = 0.0; A = 0.0; FA = 0.0
    for ts in tsets:
        v = v0.astype(complex); P = 1.0
        for t in ts[:M]:
            wv = 0.5*(v + np.exp(1j*E0*t)*(Ufun(t) @ v)); p = float(np.vdot(wv, wv).real); v = wv/np.sqrt(p); P *= p
        nb += P*nbar(v); en += P*energy(v); A += P; FA += P*abs(np.vdot(GS, v))**2
    return nb/A, en/A, 1-FA/A
rows = [("exact ground state", nbar(GS), energy(GS), 0.0),
        ("reference determinant", nbar(pro), energy(pro), 1-abs(np.vdot(GS, pro))**2),
        ("UCCSD (26 parameters)", nbar(ucc), energy(ucc), 1-abs(np.vdot(GS, ucc))**2),
        ("refined, T=30, n_s=80", nbar(ref30), energy(ref30), 1-abs(np.vdot(GS, ref30))**2)]
for M in (2, 6):
    nb, en, inf = rodeo_avg_obs(ref30, M); rows.append((f"refined + rodeo M={M} (dt=0.25)", nb, en, inf))
for M in (2,):
    nb, en, inf = rodeo_avg_obs(ucc, M); rows.append((f"UCCSD + rodeo M={M} (dt=0.25)", nb, en, inf))
print(f"{'state':32s} {'n_0':>7} {'n_1':>7} {'n_2':>7} {'n_3':>7} {'E':>10} {'E-E_ref':>9} {'1-F':>9}")
for name, nb, en, inf in rows:
    print(f"{name:32s} " + " ".join(f"{x:7.4f}" for x in nb) + f" {en:10.6f} {en-E_ref:9.5f} {inf:9.2e}")
