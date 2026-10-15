#!/usr/bin/env python3
"""Revision checks prompted by paper/article_correctness_assessment.md (k=4, N=4, g=1).

(a) Shift scan: affine fits in 1/Delta and 1/Delta^2 to the T*(mu_buf) points, the
    first-order response norm chi_u(0) = g^2/b^2 of the k=2->4 path and chi(s) along
    the path for every shift; slice convergence of the exact-slice refinement.
(b) Rodeo statistics: the exact Gaussian-ensemble infidelity against the sampled
    40-schedule values (seed reset per input) with their standard errors from
    20 000 schedules; cycles needed for 1e-6.
(c) The four-pair restricted ansatz (occupied -> virtual pair rotations only):
    energy, infidelity and 192 CNOTs, and the vanishing of the other 22 UCCSD amplitudes.
(d) The k=3->4 chain on a fine time grid and at the counted circuit level, with the
    embedded coarse term split into its coarse Trotter groups; the HF adiabatic ramp
    D -> H_high as a baseline; the prolonged k=3 state as a filter input.
(e) The long-run state of the digital (product-formula) filter: leading eigenvalue
    and asymptotic infidelity of the averaged post-selection map at dt = 0.25, 0.5, 1.
(f) Matched-observable comparison of eight inputs to the same Trotterised filter,
    including M = 0, with per-attempt and per-accepted CNOTs and the minimum cost at
    target infidelities 2e-3, 5e-4, 2e-4.
Usage: PYTHONPATH=src python3 scripts/assessment_checks.py [outdir]
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src")); sys.path.insert(0, "src")
import pairinglib as pl
from pairinglib._exc import apply_exc_circuit
from scipy.optimize import minimize
from numpy.polynomial.hermite_e import hermegauss

import sys as _sys
OUT = _sys.argv[1] if len(_sys.argv) > 1 and not _sys.argv[1].startswith("-") else None   # None inside a notebook kernel
N, G, k = 4, 1.0, 4
nq, st, ix = pl.build_sector(k, N)
H = pl.H_pairing_sparse(k, G, N, st, ix).toarray()
w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]
sen = np.array([pl.seniority_of(I, k) for I in st])
Pm2 = pl.prolong_matrix(2, k, N); HF = Pm2 @ np.array([1.0])
El2 = pl.fci_ground(pl.H_pairing_sparse(2, G, N, *pl.build_sector(2, N)[1:]))
D, hops = pl.pairing_terms(k, G, N)
Ts_grid = np.array([0.5, 1, 2, 3, 4, 6, 8, 12, 16, 20, 26, 32, 40, 50, 65, 80])
trapz = getattr(np, "trapezoid", None) or np.trapz

# =============================================================== (a) shift scan, response norm
print("(a) shift scan at g=1 (T* = first grid time with overlap >= 0.99; grid of Fig. 11):")
mubufs = [0.3, 0.6, 1.0, 1.5]
Tst, gmin = [], []
for b in mubufs:
    r = pl.refine_NK(N, k, Ts_grid, mu_buf=b)
    hit = np.where(r["ov"] >= 0.99)[0]; Tst.append(Ts_grid[hit[0]])
    A = -b*np.outer(HF, HF); B = H - (El2 + b)*np.eye(len(H))
    gmin.append(min(np.diff(np.linalg.eigvalsh(np.cos(np.pi*s/2)**2*A + np.sin(np.pi*s/2)**2*B)[:2])[0] for s in np.linspace(0, 1, 41)))
    print(f"   mu_buf={b}: min gap {gmin[-1]:.3f}, T* = {Tst[-1]}")
Tst, gmin = np.array(Tst, float), np.array(gmin)
fits = {}
for name, X in (("1/Delta", 1/gmin), ("1/Delta^2", 1/gmin**2)):
    A_ = np.vstack([X, np.ones_like(X)]).T; c, *_ = np.linalg.lstsq(A_, Tst, rcond=None)
    fits[name] = (c, np.sum((A_@c - Tst)**2))
    print(f"   affine fit T* = a*{name} + b: a = {c[0]:.2f}, b = {c[1]:+.2f}, sum of squared residuals {fits[name][1]:.2f}")
print("   -> four threshold times on a coarse grid do not identify an exponent; neither law is excluded.")
print("\n    first-order response norm chi_u(0) = sum_n |<n|dH/du|0>|^2/(E_n-E_0)^2 at u = sin^2(pi s/2) = 0, versus g^2/b^2:")
chi_s = {}
for b in mubufs:
    mu = El2 + b; A = -b*np.outer(HF, HF); B = H - mu*np.eye(len(H)); dHu = B - A
    sgrid = np.linspace(0, 1, 81); chi = []
    for s in sgrid:
        u = np.sin(np.pi*s/2)**2
        ws, Vs = np.linalg.eigh((1-u)*A + u*B)
        chi.append(sum(abs(Vs[:, n] @ dHu @ Vs[:, 0])**2/(ws[n]-ws[0])**2 for n in range(1, len(ws))))
    chi_s[b] = (sgrid, np.array(chi))
    print(f"   b = {b}: chi_u(0) = {chi[0]:.4f} (g^2/b^2 = {G**2/b**2:.4f}); max_s chi_u(s) = {max(chi):.4f} at s = {sgrid[int(np.argmax(chi))]:.2f}")
b = 0.6; mu = El2 + b; A = -b*np.outer(HF, HF); B = H - mu*np.eye(len(H))
s_ = 0.25; ws, Vs = np.linalg.eigh(np.cos(np.pi*s_/2)**2*A + np.sin(np.pi*s_/2)**2*B)
dHs = (np.pi/2)*np.sin(np.pi*s_)*(B - A); cpl = np.abs(Vs.T @ dHs @ Vs[:, 0])
s0 = [n for n in range(1, len(ws)) if np.sum(np.abs(Vs[sen == 0, n])**2) > 0.5]
print(f"   at s = 0.25 (mu_buf = 0.6) the derivative couplings |<n|dH/ds|0>| to the five seniority-zero excited states are "
      + ", ".join(f"{cpl[n]:.4f}" for n in s0) + f"; ||B-A|| = {np.linalg.norm(B-A, 2):.4f} in the full sector and "
      f"{np.linalg.norm((B-A)[np.ix_(sen==0, sen==0)], 2):.4f} in the seniority-zero block")
w0w = [np.sum(np.abs(V[sen == 0, i])**2) for i in range(len(w))]
i_s0 = next(i for i in range(1, len(w)) if w0w[i] > 0.5)
print(f"   gaps of H_high: N=4 sector {w[1]-w[0]:.4f}; first seniority-zero excitation {w[i_s0]-w[0]:.5f}")
print("\n    slice convergence of the exact-slice refinement (mu_buf = 0.6):")
for Tt in (25.0, 30.0):
    for Msl in (160, 1280):
        psi, _ = pl.refine_state(N, k, Tt, M=Msl)
        print(f"   T = {Tt:.0f}, {Msl:4d} slices: p = {abs(np.vdot(GS, psi))**2:.10f}, <H> = {np.real(np.vdot(psi, H@psi)):.10f}")

# =============================================================== (b) rodeo statistics
print("\n(b) rodeo filter at E = E_0, sigma = 3, exact evolution, input = exact-slice refined state at T = 25:")
psi25, _ = pl.refine_state(N, k, 25.0); p25 = abs(np.vdot(GS, psi25))**2
c2 = np.abs(V.T @ psi25)**2; sig = 3.0
q = 0.5*(1 + np.exp(-(w-E0)**2*sig**2/2))
Ms = np.arange(1, 11)
exact = np.array([np.sum(c2[1:]*q[1:]**M)/np.sum(c2*q**M) for M in Ms])
rng = np.random.default_rng(7)                       # convention: the generator is reset to seed 7 for every input
F = np.zeros((40, 10)); Acc = np.zeros((40, 10))
for s_ in range(40):
    F[s_], Acc[s_] = pl.rodeo_track(psi25, H, rng.normal(0, sig, 10), E0, GS)
samp_w = 1 - (F*Acc).sum(0)/Acc.sum(0); samp_u = 1 - F.mean(0)
rng = np.random.default_rng(123); nS = 20000
res = np.zeros((nS, 10)); dE = (w - E0)[None, :]
for s_ in range(nS):
    ts = rng.normal(0, sig, 10)
    fac = np.cumprod(np.cos(dE*ts[:, None]/2)**2, axis=0)
    res[s_] = (fac[:, 1:]*c2[None, 1:]).sum(1)
rse40 = res.std(0)/res.mean(0)/np.sqrt(40)
print(f"   p = {p25:.7f}; {'M':>3s} {'exact ensemble':>15s} {'sampled (40, weighted)':>22s} {'sampled (unweighted)':>20s} {'rel. std. error of a 40-set mean':>32s}")
for M in Ms:
    print(f"   {M:3d} {exact[M-1]:15.4e} {samp_w[M-1]:22.4e} {samp_u[M-1]:20.4e} {100*rse40[M-1]:31.1f}%")
Mneed = next(M for M in range(1, 40) if np.sum(c2[1:]*q[1:]**M)/np.sum(c2*q**M) <= 1e-6)
print(f"   ensemble infidelity <= 1e-6 needs M = {Mneed} cycles for this p")
print("   detuned target E_obj = E_0 + d (exact ensemble, sigma = 3, six cycles, same refined input): infidelity and acceptance")
for d in (0.0, 0.0007, 0.05, 0.1, 0.2):
    qd = 0.5*(1 + np.exp(-(w-E0-d)**2*sig**2/2)); tot = c2*qd**6
    print(f"      d = {d:6.4f} (|d| sigma sqrt(M) = {abs(d)*sig*np.sqrt(6):.3f}): 1-F_6 = {1-tot[0]/tot.sum():.2e}, acceptance P_6 = {tot.sum():.4f}")
rodeo_stats = dict(p=p25, Ms=Ms, exact=exact, samp_w=samp_w, rse40=rse40)

# =============================================================== (c) four-pair ansatz
print("\n(c) the four-pair restricted ansatz (occupied pair -> virtual pair rotations, manuscript order):")
su = pl.gate_uccsd_setup(k, N); singles, doubles = pl.uccsd_pool(k, N)
pair_idx = [m for m, (i, j, a, b_) in enumerate(doubles) if j == i+1 and i % 2 == 0 and b_ == a+1 and a % 2 == 0]
terms4 = [su["terms"][len(singles)+m] for m in pair_idx]; gens4 = [su["gens"][len(singles)+m] for m in pair_idx]
Hfull = pl.build_H_full(k, G); hf = su["hf"]; dim = 2**nq
def _state(th, terms):
    psi = np.zeros(dim, complex); psi[hf] = 1.0
    for m in range(len(terms)): psi = apply_exc_circuit(psi, th[m], terms[m], nq)
    return psi
def _energy(th, terms):
    psi = _state(th, terms); return float(np.real(psi.conj() @ (Hfull @ psi)))
def _grad(th, terms, gens):
    psi = np.zeros(dim, complex); psi[hf] = 1.0; stt = [psi.copy()]
    for m in range(len(terms)): psi = apply_exc_circuit(psi, th[m], terms[m], nq); stt.append(psi.copy())
    Hpsi = Hfull @ stt[-1]; gg = np.zeros(len(terms))
    for mu_ in range(len(terms)):
        v = gens[mu_] @ stt[mu_+1]
        for m in range(mu_+1, len(terms)): v = apply_exc_circuit(v, th[m], terms[m], nq)
        gg[mu_] = 2.0*float(np.real(v.conj() @ Hpsi))
    return gg
r4 = minimize(lambda t: _energy(t, terms4), np.zeros(4), jac=lambda t: _grad(t, terms4, gens4), method="L-BFGS-B", options={"ftol": 1e-14, "gtol": 1e-11})
psi4 = _state(r4.x, terms4); sec4 = np.array([psi4[I] for I in st]); sec4 /= np.linalg.norm(sec4)
cn4 = sum(2*(len(P)-1) for tm in terms4 for (P, _) in tm)
rfull = minimize(lambda t: _energy(t, su["terms"]), np.zeros(su["P"]), jac=lambda t: _grad(t, su["terms"], su["gens"]), method="L-BFGS-B", options={"ftol": 1e-14, "gtol": 1e-11})
psiU = _state(rfull.x, su["terms"]); secU = np.array([psiU[I] for I in st]); secU /= np.linalg.norm(secU)
others = np.abs(np.delete(rfull.x, [len(singles)+m for m in pair_idx]))
print(f"   four factors, {cn4} CNOTs: E = {r4.fun:.10f}, 1-F = {1-abs(np.vdot(GS, sec4))**2:.4e}, amplitudes {np.round(r4.x, 5)}")
print(f"   full 26-operator UCCSD ({su['ncnot']} CNOTs): E = {rfull.fun:.10f}, 1-F = {1-abs(np.vdot(GS, secU))**2:.4e}; largest of the other 22 amplitudes {others.max():.1e}")

# the k=3 coarse state (three seniority-zero determinants) with two pair rotations instead of the 272-CNOT template
su3 = pl.gate_uccsd_setup(3, N); sg3, db3 = pl.uccsd_pool(3, N)
pair3 = [m for m, (i, j, a, b_) in enumerate(db3) if j == i+1 and i % 2 == 0 and b_ == a+1 and a % 2 == 0]
terms3 = [su3["terms"][len(sg3)+m] for m in pair3]; H3full = pl.build_H_full(3, G); hf3 = su3["hf"]
def _state3(th):
    psi = np.zeros(2**6, complex); psi[hf3] = 1.0
    for m in range(len(terms3)): psi = apply_exc_circuit(psi, th[m], terms3[m], 6)
    return psi
E3fn = lambda th: float(np.real(np.vdot(_state3(th), H3full @ _state3(th))))
r3 = min((minimize(E3fn, x0, method="BFGS", options={"gtol": 1e-12}) for x0 in (np.zeros(2), np.array([0.3, 0.2]), np.array([-0.3, 0.5]))), key=lambda r: r.fun)
n3, s3, i3 = pl.build_sector(3, N); H3 = pl.H_pairing_sparse(3, G, N, s3, i3).toarray(); w3, V3 = np.linalg.eigh(H3)
psi3b = _state3(r3.x); sec3b = np.array([psi3b[I] for I in s3]); sec3b /= np.linalg.norm(sec3b)
cn3 = sum(2*(len(P)-1) for tm in terms3 for (P, _) in tm)
print(f"   k=3 coarse state with two pair rotations ({cn3} CNOTs): E = {r3.fun:.10f} (FCI {w3[0]:.10f}), 1-F = {1-abs(np.vdot(V3[:,0], sec3b))**2:.1e}; the 272-CNOT template gives the same state")

# =============================================================== (d) 3->4 chain, HF ramp, prolonged k=3
print("\n(d) the k=3->4 chain from the gate-level UCCSD state at k=3:")
E3, psi3 = pl.gate_uccsd_state(3, N, G)
sec3 = np.array([psi3[I] for I in s3]); sec3 /= np.linalg.norm(sec3)
Tfine = np.array([4, 5, 6, 7, 8, 9, 10, 12, 16.0])
r34 = pl.refine_from_state(sec3, N, 3, 4, Tfine)
print("   exact slices: " + ", ".join(f"p({T:.0f}) = {o:.5f}" for T, o in zip(Tfine, r34["ov"])) + f"  (p0 = {r34['ov0']:.6f})")
T34 = Tfine[np.where(r34["ov"] >= 0.99)[0][0]]
print(f"   T* (overlap >= 0.99) on this grid: {T34:.0f}")
# split coarse groups: A = P (H_low - mu) P^T = P(D_low - mu)P^T + sum_pq P h_pq^low P^T
Pm3 = pl.prolong_matrix(3, k, N); Hl3 = pl.H_pairing_sparse(3, G, N, s3, i3).toarray()
El3 = float(np.real(sec3 @ Hl3 @ sec3)); mu3 = El3 + 0.6
D3, hops3 = pl.pairing_terms(3, G, N)
groups34 = [Pm3 @ (D3 - mu3*np.eye(len(D3))) @ Pm3.T] + [Pm3 @ h @ Pm3.T for h in hops3] + [D - mu3*np.eye(len(D))] + hops
c2f = lambda s: np.cos(np.pi*s/2)**2; s2f = lambda s: np.sin(np.pi*s/2)**2
coef34 = [c2f]*(1+len(hops3)) + [s2f]*(1+len(hops))
Phi34 = Pm3 @ sec3
print(f"   Trotterised at T = {T34:.0f} with the coarse term split into its {1+len(hops3)} coarse groups (second order):")
p34_tr = {}
for ns in (8, 11, 16, 22, 32, 48):
    psi = pl.trotter_path(Phi34, groups34, coef34, float(T34), ns, order=2)
    p34_tr[ns] = abs(np.vdot(GS, psi))**2
    print(f"      n_s = {ns:3d}: p = {p34_tr[ns]:.6f}")
CN_STEP, CN_CSTEP, CN_A2, CN_A3, CN_UCC3_TEMPLATE = 296, 416, 14, 350, 272
CN_UCC3 = cn3                                           # pruned two-rotation preparation of the same state
pex34 = r34["ov"][list(Tfine).index(T34)]
ns34 = next(ns for ns in (8, 11, 16, 22, 32, 48) if abs(p34_tr[ns] - pex34) <= 1 - pex34 and p34_tr[ns] >= 0.99)
print(f"   steps reproducing the exact-slice overlap within 1-p AND reaching the 0.99 target: n_s = {ns34}; cost 2*{ns34}*({CN_STEP}+{CN_A3}) + {CN_UCC3} = {2*ns34*(CN_STEP+CN_A3)+CN_UCC3} CNOTs (with the 272-CNOT template: {2*ns34*(CN_STEP+CN_A3)+CN_UCC3_TEMPLATE})")

print("\n    the HF adiabatic ramp H(s) = D + sin^2(pi s/2) sum_pq h_pq (no coarse term, no shift):")
groupsR = [D] + hops; coefR = [lambda s: 1.0] + [s2f]*len(hops)
def ramp_exact(T, Msl=160):
    psi = HF.astype(complex).copy(); dt = T/Msl
    eig = [np.linalg.eigh(D + np.sin(np.pi*((m+0.5)/Msl)/2)**2*sum(hops)) for m in range(Msl)]
    for lam, Vv in eig: psi = Vv @ (np.exp(-1j*lam*dt)*(Vv.conj().T @ psi))
    return psi
pR = {T: abs(np.vdot(GS, ramp_exact(T)))**2 for T in Tfine}
print("   exact slices: " + ", ".join(f"p({T:.0f}) = {pR[T]:.5f}" for T in Tfine))
TR = next(T for T in Tfine if pR[T] >= 0.99)
pR_tr = {}
for ns in (4, 6, 8, 11, 16, 22):
    pR_tr[ns] = abs(np.vdot(GS, pl.trotter_path(HF, groupsR, coefR, float(TR), ns, order=2)))**2
print(f"   T* = {TR:.0f}; Trotterised (second order): " + ", ".join(f"n_s={ns}: {pR_tr[ns]:.6f}" for ns in pR_tr))
nsR = next(ns for ns in pR_tr if abs(pR_tr[ns] - pR[TR]) <= 1 - pR[TR] and pR_tr[ns] >= 0.99)
print(f"   steps reproducing the exact-slice overlap within 1-p and reaching 0.99: n_s = {nsR}; cost 2*{nsR}*{CN_STEP} = {2*nsR*CN_STEP} CNOTs (no coarse term)")
D_gap_s0 = np.sort(np.diag(D)[sen == 0]); print(f"   gap of D within seniority zero at s = 0: {D_gap_s0[1]-D_gap_s0[0]:.3f} (ramp); shift-controlled gap of the refinement path: 0.600")

# =============================================================== (e) long-run state of the digital filter
print("\n(e) the averaged post-selection map K(rho) = E_t[K_t rho K_t^+] of the Trotterised filter (second order, |t| half-Gaussian, sigma = 4):")
tm = pl.model_terms(k, G, N)
idx0 = np.where(sen == 0)[0]; d0 = len(idx0); GS0 = GS[idx0]/np.linalg.norm(GS[idx0])
xg, wg = hermegauss(80); wg = wg/np.sqrt(2*np.pi)
floors = {}
for dtm in (0.25, 0.5, 1.0):
    K = np.zeros((d0*d0, d0*d0), complex)
    for xk, wk in zip(xg, wg):
        t = abs(4.0*xk)
        U = pl.trotter_U_dt(k, G, N, t, dtm, order=2, terms=tm)[0][np.ix_(idx0, idx0)]
        Kt = 0.5*(np.eye(d0) + np.exp(1j*E0*t)*U)
        K += wk*np.kron(Kt, Kt.conj())
    lam, Vk = np.linalg.eig(K); j = np.argmax(np.abs(lam))
    rho = Vk[:, j].reshape(d0, d0); rho = rho/np.trace(rho); rho = 0.5*(rho + rho.conj().T)
    floors[dtm] = (abs(lam[j]), 1 - np.real(GS0 @ rho @ GS0))
    print(f"   dt_max = {dtm:4.2f}: leading eigenvalue {floors[dtm][0]:.6f} (long-run acceptance decays by {100*(1-floors[dtm][0]):.2f}% per cycle), asymptotic infidelity {floors[dtm][1]:.3e}")

# =============================================================== (f) matched comparison, eight inputs, M = 0..6
print("\n(f) the same Trotterised filter (second order, dt = 0.25, sigma = 4, 40 time sets, seed 3) on eight inputs:")
DT, SIG, NSAMP, M_MAX, SEED = 0.25, 4.0, 40, 6, 3
Ufun = lambda t: pl.trotter_U_dt(k, G, N, t, DT, order=2, terms=tm)[0]
NS_DENS = 80/30.0
inputs = {}
inputs["prolonged determinant (k=2)"] = (HF.astype(complex), 0)
inputs["four-pair ansatz"] = (sec4.astype(complex), cn4)
inputs["UCCSD, 26 operators"] = (secU.astype(complex), su["ncnot"])
inputs["prolonged k=3 (two rotations)"] = (Phi34.astype(complex), CN_UCC3)
inputs[f"HF ramp, T={TR:.0f} (n_s={nsR})"] = (pl.trotter_path(HF, groupsR, coefR, float(TR), nsR, order=2), 2*nsR*CN_STEP)
for T in (16, 30):
    ns = int(np.ceil(NS_DENS*T)); psi, _ = pl.refine_state_trotter(N, k, float(T), ns, order=2)
    inputs[f"refined 2->4, T={T} (n_s={ns})"] = (psi, 2*ns*(CN_STEP + CN_A2))
inputs[f"refined 3->4, T={T34:.0f} (n_s={ns34})"] = (pl.trotter_path(Phi34, groups34, coef34, float(T34), ns34, order=2), CN_UCC3 + 2*ns34*(CN_STEP + CN_A3))
rng = np.random.default_rng(SEED)
tsets = [np.abs(rng.normal(0.0, SIG, M_MAX)) for _ in range(NSAMP)]
nsteps = np.array([[max(1, int(np.ceil(t/DT))) for t in ts] for ts in tsets]); cyc_cost = nsteps*CN_CSTEP*2
rows = {}
for name, (v0, cprep) in inputs.items():
    p = abs(np.vdot(GS, v0))**2
    Fm = np.zeros((NSAMP, M_MAX)); Am = np.zeros((NSAMP, M_MAX))
    for s_, ts in enumerate(tsets):
        Fm[s_], Am[s_] = pl.rodeo_track_U(v0, Ufun, ts, E0, GS)
    Fw = (Fm*Am).sum(0)/Am.sum(0); Amean = Am.mean(0)
    reach = np.concatenate([np.ones((NSAMP, 1)), Am[:, :-1]], axis=1)
    cost_acc = (cprep + np.cumsum(reach*cyc_cost, axis=1).mean(0))/Amean
    cost_att = cprep + np.cumsum(cyc_cost, axis=1).mean(0)                  # CNOTs of one attempt that runs all M cycles
    inf = np.concatenate([[1-p], 1-Fw]); acc = np.concatenate([[1.0], Amean])
    rows[name] = dict(p=p, cprep=cprep, inf=inf, acc=acc, cost_acc=np.concatenate([[float(cprep)], cost_acc]),
                      cost_att=np.concatenate([[float(cprep)], cost_att]))
print(f"   {'input':34s} {'p':>8s} {'prep':>6s} | " + " | ".join(f"M={m}: 1-F_w    acc   /attempt /accepted" for m in (0, 1, 2, 4, 6)))
for name, d in rows.items():
    print(f"   {name:34s} {d['p']:8.5f} {d['cprep']:6d} | " + " | ".join(f"{d['inf'][m]:.1e} {d['acc'][m]:.3f} {d['cost_att'][m]:8.0f} {d['cost_acc'][m]:8.0f}" for m in (0, 1, 2, 4, 6)))
print("\n    minimum CNOTs per accepted preparation reaching a target infidelity (over inputs and M = 0..6; floor of the filter ~1e-4):")
for tgt in (2e-3, 5e-4, 2e-4):
    best = []
    for name, d in rows.items():
        ok = np.where(d["inf"] <= tgt)[0]
        if ok.size:
            m = ok[np.argmin(d["cost_acc"][ok])]; best.append((d["cost_acc"][m], name, m, d["inf"][m]))
    best.sort()
    print(f"   target {tgt:.0e}: " + "; ".join(f"{nm} M={m} ({inf_:.1e}) {c:.0f}" for c, nm, m, inf_ in best[:4]))

# =============================================================== figures
import matplotlib
if OUT: matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
x = 1/gmin; xs = np.linspace(0.5, x.max()*1.08, 100)
ax[0].plot(x, Tst, "o", ms=8, label="$T^*$ (overlap $\\geq0.99$)")
c1, _ = fits["1/Delta"]; c2_, _ = fits["1/Delta^2"]
ax[0].plot(xs, c1[0]*xs + c1[1], "--", label=f"$a/\\Delta_{{\\min}}+b$: $a={c1[0]:.1f}$, $b={c1[1]:+.1f}$, SSR $={fits['1/Delta'][1]:.1f}$")
ax[0].plot(xs, c2_[0]*xs**2 + c2_[1], "-", label=f"$a/\\Delta_{{\\min}}^2+b$: $a={c2_[0]:.1f}$, $b={c2_[1]:+.1f}$, SSR $={fits['1/Delta^2'][1]:.2f}$")
ax[0].set_xlabel("$1/\\Delta_{\\min}$ (shift scan, $g=1$)"); ax[0].set_ylabel("$T^*$"); ax[0].set_title("two affine fits to the four shift-scan points"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
for b_, (sg, ch) in chi_s.items():
    ax[1].plot(sg, ch, label=f"$\\mu_{{\\rm buf}}={b_}$ ($\\chi_u(0)=g^2/b^2={G**2/b_**2:.2f}$)")
ax[1].set_yscale("log"); ax[1].set_xlabel("$s$"); ax[1].set_ylabel("$\\chi_u(s)=\\sum_{n>0}|\\langle n|\\partial_uH|0\\rangle|^2/(E_n-E_0)^2$")
ax[1].set_title("first-order response norm along the $k=2\\to4$ path"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3, which="both")
plt.tight_layout()
if OUT: plt.savefig(os.path.join(OUT, "fig19.png"), dpi=150)
else: plt.show()

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
mks = ["o-", "s-", "d-", "^-", "v-", "P-", "X-", "*-"]
for (name, d), mk_ in zip(rows.items(), mks):
    m0 = 1 if d["cprep"] == 0 else 0
    ax[0].loglog(d["cost_acc"][m0:], d["inf"][m0:], mk_, ms=5, label=f"{name} ($p={d['p']:.4f}$)")
    ax[1].semilogy(np.arange(0, M_MAX+1), d["inf"], mk_, ms=5, label=name)
for tgt in (2e-3, 5e-4, 2e-4):
    ax[0].axhline(tgt, color="gray", ls=":", lw=0.8)
ax[0].set_xlabel("CNOTs per accepted preparation (input + filter, early termination)"); ax[0].set_ylabel("ensemble infidelity $1-\\mathcal{F}_M$")
ax[0].set_title("eight inputs to the same Trotterised filter ($\\delta t=0.25$, $\\sigma=4$); $M=0\\ldots6$"); ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3, which="both")
ax[1].set_xlabel("rodeo cycles $M$ ($M=0$: unfiltered input)"); ax[1].set_ylabel("$1-\\mathcal{F}_M$"); ax[1].set_title("same data versus the number of cycles"); ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3, which="both")
plt.tight_layout()
if OUT: plt.savefig(os.path.join(OUT, "fig20.png"), dpi=150); print("wrote", OUT)
else: plt.show()
