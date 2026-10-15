#!/usr/bin/env python3
"""Independent checks of the quantitative claims in paper/article_correctness_assessment.md
(k=4, N=4, g=1 unless stated).  PYTHONPATH=src python3 scripts/verify_assessment.py"""
import numpy as np, sys, os
sys.path.insert(0, "src")
import pairinglib as pl
from pairinglib._exc import make_double, apply_exc_circuit, jw_create, jw_ann
from scipy.linalg import expm
from scipy.optimize import minimize

N, G, k = 4, 1.0, 4
nq, st, ix = pl.build_sector(k, N)
H = pl.H_pairing_sparse(k, G, N, st, ix).toarray()
w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]
Pm = pl.prolong_matrix(2, k, N); HF = Pm @ np.array([1.0])
print(f"[1] E0 = {E0:.13f}, p(HF) = {abs(GS @ HF)**2:.13f}  (assessment: 0.6355484735756, 0.8676375677137)")
print(f"    N=4 sector gap = {w[1]-w[0]:.4f}; seniority-zero weights of the lowest states:")
sen = np.array([pl.seniority_of(I, k) for I in st])
w0 = [np.sum(np.abs(V[sen == 0, i])**2) for i in range(len(w))]
first_s0 = next(i for i in range(1, len(w)) if w0[i] > 0.5)
print(f"    first excited state with seniority-zero weight: E_{first_s0}-E_0 = {w[first_s0]-w[0]:.5f} (assessment: 2.29983)")

# ---- [2] shift-scan fits
D = np.array([0.3, 0.6, 1.0, 1.499]); T = np.array([50, 16, 8, 6.0])
for name, X in (("a/D + b", 1/D), ("a/D^2 + b", 1/D**2)):
    A = np.vstack([X, np.ones_like(X)]).T; c, res, *_ = np.linalg.lstsq(A, T, rcond=None)
    print(f"[2] fit T* = {name}: a={c[0]:.3f} b={c[1]:.3f} SSR={np.sum((A@c-T)**2):.3f}")

# ---- [3] response norm chi_u(0) = g^2/b^2 and derivative couplings at s=0.25
El = pl.fci_ground(pl.H_pairing_sparse(2, G, N, *pl.build_sector(2, N)[1:]))
for b in (0.3, 0.6, 1.0, 1.5):
    mu = El + b; A = -b*np.outer(HF, HF); Bm = H - mu*np.eye(len(H))
    dH = Bm - A                                     # dH/du, u = sin^2(pi s/2)
    wa, Va = np.linalg.eigh(A); i0 = np.argmin(wa)
    chi = sum(abs(Va[:, n] @ dH @ Va[:, i0])**2/(wa[n]-wa[i0])**2 for n in range(len(wa)) if n != i0)
    print(f"[3] b={b}: chi_u(0) = {chi:.4f}  (g^2/b^2 = {G**2/b**2:.4f})")
b = 0.6; mu = El + b; A = -b*np.outer(HF, HF); Bm = H - mu*np.eye(len(H))
s = 0.25; Hs = np.cos(np.pi*s/2)**2*A + np.sin(np.pi*s/2)**2*Bm; ws, Vs = np.linalg.eigh(Hs)
dHs = (np.pi/2)*np.sin(np.pi*s)*(Bm - A)
cpl = np.abs(Vs[:, 1:].T @ dHs @ Vs[:, 0])
sen_w = np.array([np.sum(np.abs(Vs[sen == 0, i])**2) for i in range(len(ws))])
s0_idx = [i for i in range(len(ws)) if sen_w[i] > 0.5]
print(f"[3] s=0.25: couplings |<n|dH/ds|0>| to the seniority-zero states: " + ", ".join(f"{cpl[i-1]:.5f}" for i in s0_idx if i > 0))
print(f"    ||B-A|| full sector = {np.linalg.norm(Bm-A, 2):.7f}; in the seniority-zero block = "
      f"{np.linalg.norm((Bm-A)[np.ix_(sen==0, sen==0)], 2):.7f}  (assessment: 7.6205693)")

# ---- [4] slice convergence of the exact-slice refinement
for Tt in (25.0, 30.0):
    for Msl in (160, 1280):
        psi, _ = pl.refine_state(N, k, Tt, M=Msl)
        print(f"[4] T={Tt:.0f} slices={Msl}: p = {abs(np.vdot(GS, psi))**2:.10f}, <H> = {np.real(np.vdot(psi, H@psi)):.10f}")

# ---- [5] rodeo: sampled (seed 7, reset / advanced) vs exact ensemble, sigma=3
psi25, _ = pl.refine_state(N, k, 25.0); p25 = abs(np.vdot(GS, psi25))**2
c2 = np.abs(V.T @ psi25)**2; sig = 3.0
q = 0.5*(1 + np.exp(-(w-E0)**2*sig**2/2))
exact = {M: np.sum(c2[1:]*q[1:]**M)/np.sum(c2*q**M) for M in (1, 2, 6, 10)}
def sample(seed, advance_first):
    rng = np.random.default_rng(seed)
    if advance_first:
        for _ in range(40): rng.normal(0, sig, 10)
    F = np.zeros((40, 10)); Acc = np.zeros((40, 10))
    for s_ in range(40):
        F[s_], Acc[s_] = pl.rodeo_track(psi25, H, rng.normal(0, sig, 10), E0, GS)
    return 1-F.mean(0), 1-(F*Acc).sum(0)/Acc.sum(0)
unw, wgt = sample(7, False); unw2, wgt2 = sample(7, True)
print(f"[5] p(T=25) = {p25:.10f}")
for M in (1, 2, 6, 10):
    print(f"    M={M:2d}: sampled unweighted {unw[M-1]:.5e}  weighted {wgt[M-1]:.5e}  (generator advanced through Phi0 first: {wgt2[M-1]:.5e})   exact ensemble {exact[M]:.5e}")
# relative standard error of the mean residual weight for 40 schedules (Monte Carlo over 20000 schedules)
rng = np.random.default_rng(123); nS = 20000
res = np.zeros((nS, 10))
dE = (w - E0)[None, :]
for s_ in range(nS):
    ts = rng.normal(0, sig, 10)
    fac = np.cumprod(np.cos(dE*ts[:, None]/2)**2, axis=0)      # (M, states)
    res[s_] = (fac[:, 1:]*c2[None, 1:]).sum(1)
rse = res.std(0)/res.mean(0)/np.sqrt(40)
print("    relative standard error of the 40-schedule mean residual weight: " + ", ".join(f"M={M}: {100*rse[M-1]:.1f}%" for M in (2, 6, 10)) + "  (assessment: 12.7, 36.6, 85.5)")
Mneed = next(M for M in range(1, 40) if np.sum(c2[1:]*q[1:]**M)/np.sum(c2*q**M) <= 1e-6)
print(f"    cycles for ensemble infidelity <= 1e-6: {Mneed} (assessment: 12)")

# ---- [6] four-pair restricted ansatz (occupied pairs -> virtual pairs), gate-level, manuscript order
su = pl.gate_uccsd_setup(k, N)
singles, doubles = pl.uccsd_pool(k, N)
pair_idx = [m for m, (i, j, a, b) in enumerate(doubles) if j == i+1 and i % 2 == 0 and b == a+1 and a % 2 == 0]
terms = [su["terms"][len(singles)+m] for m in pair_idx]; gens = [su["gens"][len(singles)+m] for m in pair_idx]
Hfull = pl.build_H_full(k, G); hf = su["hf"]; dim = 2**nq
def state(th):
    psi = np.zeros(dim, complex); psi[hf] = 1.0
    for m in range(len(terms)): psi = apply_exc_circuit(psi, th[m], terms[m], nq)
    return psi
def energy(th):
    psi = state(th); return float(np.real(psi.conj() @ (Hfull @ psi)))
def grad(th):
    psi = np.zeros(dim, complex); psi[hf] = 1.0; stt = [psi.copy()]
    for m in range(len(terms)): psi = apply_exc_circuit(psi, th[m], terms[m], nq); stt.append(psi.copy())
    Hpsi = Hfull @ stt[-1]; gg = np.zeros(len(terms))
    for mu_ in range(len(terms)):
        v = gens[mu_] @ stt[mu_+1]
        for m in range(mu_+1, len(terms)): v = apply_exc_circuit(v, th[m], terms[m], nq)
        gg[mu_] = 2.0*float(np.real(v.conj() @ Hpsi))
    return gg
r4 = minimize(energy, np.zeros(len(terms)), jac=grad, method="L-BFGS-B", options={"ftol": 1e-14, "gtol": 1e-11})
psi4 = state(r4.x); sec4 = np.array([psi4[I] for I in st]); sec4 /= np.linalg.norm(sec4)
ncn = sum(2*(len(P)-1) for tm in terms for (P, _) in tm)
Eu, psiu = pl.gate_uccsd_state(k, N, G); secu = np.array([psiu[I] for I in st]); secu /= np.linalg.norm(secu)
print(f"[6] four-pair ansatz: {len(terms)} factors, {ncn} CNOTs, E = {r4.fun:.13f}, 1-F = {1-abs(np.vdot(GS, sec4))**2:.10e}  (assessment: 0.6369869810446, 1.7764427147e-4)")
print(f"    full 26-operator UCCSD: E = {Eu:.13f}, 1-F = {1-abs(np.vdot(GS, secu))**2:.10e}")

# ---- [7] 3->4 chain: ideal p at T=8 from the k=3 gate-UCCSD state
E3, psi3 = pl.gate_uccsd_state(3, N, G); n3, s3, i3 = pl.build_sector(3, N)
sec3 = np.array([psi3[I] for I in s3]); sec3 /= np.linalg.norm(sec3)
r = pl.refine_from_state(sec3, N, 3, 4, np.array([8.0, 9.0]))
print(f"[7] 3->4 chain from gate-UCCSD k=3: p0 = {r['ov0']:.6f}, p(T=8) = {r['ov'][0]:.6f}, p(T=9) = {r['ov'][1]:.6f}  (assessment: 0.963668, 0.990737)")

# ---- [8] digital rodeo map (second order, dt_max=0.25, sigma=4, |t| half-Gaussian): asymptotic infidelity
tm = pl.model_terms(k, G, N)
def Ut(t, dtm):
    return pl.trotter_U_dt(k, G, N, t, dtm, order=2, terms=tm)[0]
# restrict to the seniority-zero block (6 states), which the filter preserves
idx0 = np.where(sen == 0)[0]; d0 = len(idx0)
GS0 = GS[idx0]/np.linalg.norm(GS[idx0])
from numpy.polynomial.hermite_e import hermegauss
xg, wg = hermegauss(60)                      # Gauss-Hermite for N(0,1): E[f] = sum w f(x)/sqrt(2pi)
wg = wg/np.sqrt(2*np.pi)
for dtm in (0.25, 0.5):
    K = np.zeros((d0*d0, d0*d0), complex)
    for xk, wk in zip(xg, wg):
        t = abs(4.0*xk)
        Kt = 0.5*(np.eye(d0) + np.exp(1j*E0*t)*Ut(t, dtm)[np.ix_(idx0, idx0)])
        K += wk*np.kron(Kt, Kt.conj())
    lam, Vk = np.linalg.eig(K); j = np.argmax(np.abs(lam))
    rho = Vk[:, j].reshape(d0, d0); rho = rho/np.trace(rho); rho = 0.5*(rho + rho.conj().T)
    print(f"[8] dt_max={dtm}: leading eigenvalue {abs(lam[j]):.7f}, asymptotic infidelity 1-<GS|rho|GS> = {1-np.real(GS0 @ rho @ GS0):.5e}  (assessment: 0.9978208 / 1.1399e-4 at 0.25; 1.9345e-3 at 0.5)")
