#!/usr/bin/env python3
"""Sufficient conditions for the Trotterized refinement path (paper Eq. eq:Tbound),
evaluated for the k=2->4 chain at N=4, g=1 and compared with the empirical
refinement time T* and step count: (a) the adiabatic-time bound of Jansen,
Ruskai and Seiler in the form used by the nucleation paper (Sadia et al.,
arXiv:2609.37527, Theorem 1), T >= eps^{-1/2} int_0^1 [|H''|/Delta^2
+ 7|H'|^2/Delta^3] ds; (b) the first-order Trotter-step bound (Theorem 2),
n_s >= eta^{-1}(T Lambda_1/4 + T^2 Lambda_2/2), with the multi-group
generalisation Lambda_2 = max_s sum_{i<j} |[H_i(s),H_j(s)]| for the term
partition of the gate compilation ([A, D-mu, h_pq..., V_ph pieces]).
All norms are spectral norms in the N-particle sector, where the evolution lives.
Usage: PYTHONPATH=src python3 scripts/adiabatic_bound.py"""
import sys, os
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src")); sys.path.insert(0, "src")
import pairinglib as pl

N, G = 4, 1.0
K_LOW, K_HIGH, MU_BUF = 2, 4, 0.6
EPS = (1e-2, 1e-3)                                   # target infidelities
F_LIST = (0.0, 0.05, 0.2, 0.5)
Ts = np.array([0.5, 1, 2, 3, 4, 6, 8, 12, 16, 20, 26, 32, 40, 50, 65, 80])  # T* grid of Sec. 15/16

def opnorm(M):
    return float(np.linalg.norm(M, 2))

def path_pieces(f):
    """A = P(H_low-mu)P^T, B = H_high-mu and the Trotter groups of the compilation."""
    nl, sl, il = pl.build_sector(K_LOW, N)
    Hl = pl.H_pairing_sparse(K_LOW, G, N, sl, il, f=f).toarray(); El = np.linalg.eigvalsh(Hl)[0]
    Pm = pl.prolong_matrix(K_LOW, K_HIGH, N)
    nh, sh, ih = pl.build_sector(K_HIGH, N)
    Hh = pl.H_pairing_sparse(K_HIGH, G, N, sh, ih, f=f).toarray()
    mu = El + MU_BUF
    A = Pm @ (Hl - mu*np.eye(len(Hl))) @ Pm.T; B = Hh - mu*np.eye(len(Hh))
    tm = pl.model_terms(K_HIGH, G, N, f)
    groups = [A, tm[0] - mu*np.eye(len(B))] + tm[1:]      # [A, D-mu, h_pq..., V_ph...]
    return A, B, groups

def bound_table(f, s_grid=np.linspace(0, 1, 401), s_comm=np.linspace(0, 1, 41)):
    A, B, groups = path_pieces(f)
    nBA = opnorm(B - A)
    # H(s) = cos^2(pi s/2) A + sin^2(pi s/2) B  ->  H' = (pi/2) sin(pi s)(B-A), H'' = (pi^2/2) cos(pi s)(B-A)
    gap = np.array([np.diff(np.linalg.eigvalsh(np.cos(np.pi*s/2)**2*A + np.sin(np.pi*s/2)**2*B)[:2])[0]
                    for s in s_grid])
    dH1 = (np.pi/2)*np.abs(np.sin(np.pi*s_grid))*nBA
    dH2 = (np.pi**2/2)*np.abs(np.cos(np.pi*s_grid))*nBA
    integrand = dH2/gap**2 + 7*dH1**2/gap**3
    trapz = getattr(np, "trapezoid", None) or np.trapz
    I = trapz(integrand, s_grid)
    I_second = trapz(dH2/gap**2, s_grid)
    Lam1 = (np.pi/2)*nBA
    def coef(s):
        c2, s2 = np.cos(np.pi*s/2)**2, np.sin(np.pi*s/2)**2
        return [c2] + [s2]*(len(groups) - 1)
    Lam2 = 0.0
    for s in s_comm:
        c = coef(s); tot = 0.0
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                Ci, Cj = c[i]*groups[i], c[j]*groups[j]
                tot += opnorm(Ci @ Cj - Cj @ Ci)
        Lam2 = max(Lam2, tot)
    r = pl.refine_NK(N, K_HIGH, Ts, k_low=K_LOW, f=f)
    hit = np.where(r["ov"] >= 0.99)[0]; Tstar = Ts[hit[0]] if hit.size else np.nan
    # steps needed in practice at T*: the product formula reproduces the exact-slice
    # overlap to within the adiabatic error 1-p_exact of the path itself (the criterion of Sec. 16)
    p_exact = pl.refine_NK(N, K_HIGH, np.array([Tstar]), k_low=K_LOW, f=f)["ov"][0]
    nh, sh, ih = pl.build_sector(K_HIGH, N)
    Hh = pl.H_pairing_sparse(K_HIGH, G, N, sh, ih, f=f).toarray(); GS = np.linalg.eigh(Hh)[1][:, 0]
    n_prac = {}
    for order in (1, 2):
        n_prac[order] = None
        for n in (5, 10, 20, 40, 80, 160, 320, 640, 1280):
            psi, _ = pl.refine_state_trotter(N, K_HIGH, Tstar, n, k_low=K_LOW, order=order, f=f)
            if abs(abs(np.vdot(GS, psi))**2 - p_exact) <= 1 - p_exact:
                n_prac[order] = n; break
    out = dict(f=f, gap_min=gap.min(), nBA=nBA, I=I, I_second=I_second, Lam1=Lam1, Lam2=Lam2,
               Tstar=Tstar, p_exact=p_exact, n_prac=n_prac, s=s_grid, gap=gap, integrand=integrand)
    for eps in EPS:
        Tb = I/np.sqrt(eps); eta = np.sqrt(eps)/2
        out[f"Tb{eps:g}"] = Tb
        out[f"ns_at_Tstar{eps:g}"] = (Tstar*Lam1/4 + Tstar**2*Lam2/2)/eta
        out[f"ns_at_Tb{eps:g}"] = (Tb*Lam1/4 + Tb**2*Lam2/2)/eta
    return out

res = [bound_table(f) for f in F_LIST]
print(f"Sufficient conditions along the k={K_LOW}->{K_HIGH} refinement path (N={N}, g={G}, mu=E_low+{MU_BUF}); "
      "spectral norms in the N-particle sector")
print("="*118)
print(f"{'f':>5} {'Delta_min':>9} {'|B-A|':>7} {'int':>9} {'(H'' part)':>10} "
      f"{'T_b(1e-2)':>10} {'T_b(1e-3)':>10} {'T*':>4} {'T_b/T*':>8} "
      f"{'Lam1':>6} {'Lam2':>7} {'n_s(T*)':>9} {'n_s(T_b)':>9} {'p(T*)':>7} {'n1 prac':>8} {'n2 prac':>8}")
for r in res:
    print(f"{r['f']:5.2f} {r['gap_min']:9.3f} {r['nBA']:7.2f} {r['I']:9.1f} {r['I_second']:10.1f} "
          f"{r['Tb0.01']:10.3g} {r['Tb0.001']:10.3g} {r['Tstar']:4.0f} {r['Tb0.01']/r['Tstar']:8.0f} "
          f"{r['Lam1']:6.2f} {r['Lam2']:7.2f} {r['ns_at_Tstar0.01']:9.3g} {r['ns_at_Tb0.01']:9.3g} "
          f"{r['p_exact']:7.4f} {str(r['n_prac'][1]):>8} {str(r['n_prac'][2]):>8}")
print("\nT_b(eps): adiabatic-time bound for infidelity eps; T*: empirical time for overlap 0.99 (grid of Sec. 15);")
print("n_s(T): first-order Trotter-step bound at eta = sqrt(eps)/2, eps = 1e-2, evaluated at T* and at T_b;")
print("p(T*): exact-slice overlap at T*; n1/n2 prac: first-/second-order steps (from 5,10,20,...) that reproduce")
print("it to within the adiabatic error 1-p(T*) of the path itself (the criterion of Section 16).")
r0 = res[0]
print(f"\nf=0: the H'' term contributes {100*r0['I_second']/r0['I']:.1f}% of the integral; "
      f"7|H'|^2/Delta^3 dominates because Delta_min={r0['gap_min']:.3f} and |B-A|={r0['nBA']:.2f}.")

fig, (axL, axR) = plt.subplots(1, 2, figsize=(10, 3.6))
for r in res:
    axL.plot(r['s'], r['gap'], label=f"f={r['f']:g}")
    axR.semilogy(r['s'], r['integrand'], label=f"f={r['f']:g}")
axL.set_xlabel("s = t/T"); axL.set_ylabel(r"instantaneous gap $\Delta(s)$"); axL.legend(frameon=False)
axL.set_title(r"gap along the $2\to4$ path")
axR.set_xlabel("s = t/T"); axR.set_ylabel(r"$\|\ddot H\|/\Delta^2 + 7\|\dot H\|^2/\Delta^3$")
axR.set_title(r"integrand of the time bound"); axR.legend(frameon=False)
plt.tight_layout(); plt.show()
