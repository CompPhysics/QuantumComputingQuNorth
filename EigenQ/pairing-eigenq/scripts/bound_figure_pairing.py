#!/usr/bin/env python3
"""Pure-pairing (f=0) version of the adiabatic-bound figure (paper Fig. fig:bound):
instantaneous gap and the integrand of Eq. eq:Tbound along the k=2->4 path at
N=4, g=1, mu=E_low+0.6.  Same construction as scripts/adiabatic_bound.py (which
also covers the particle-hole strengths used in the companion paper); this script
only redraws the f=0 curves so that the pure-pairing article shows the pairing
model alone.  Writes paper/figs/fig_bound_pairing.png.
Usage: PYTHONPATH=src python3 scripts/bound_figure_pairing.py"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import pairinglib as pl

N, G, K_LOW, K_HIGH, MU_BUF, F = 4, 1.0, 2, 4, 0.6, 0.0

nl, sl, il = pl.build_sector(K_LOW, N)
Hl = pl.H_pairing_sparse(K_LOW, G, N, sl, il, f=F).toarray(); El = np.linalg.eigvalsh(Hl)[0]
Pm = pl.prolong_matrix(K_LOW, K_HIGH, N)
nh, sh, ih = pl.build_sector(K_HIGH, N)
Hh = pl.H_pairing_sparse(K_HIGH, G, N, sh, ih, f=F).toarray()
mu = El + MU_BUF
A = Pm @ (Hl - mu*np.eye(len(Hl))) @ Pm.T; B = Hh - mu*np.eye(len(Hh))
nBA = float(np.linalg.norm(B - A, 2))
s = np.linspace(0, 1, 401)
gap = np.array([np.diff(np.linalg.eigvalsh(np.cos(np.pi*x/2)**2*A + np.sin(np.pi*x/2)**2*B)[:2])[0] for x in s])
dH1 = (np.pi/2)*np.abs(np.sin(np.pi*s))*nBA
dH2 = (np.pi**2/2)*np.abs(np.cos(np.pi*s))*nBA
integrand = dH2/gap**2 + 7*dH1**2/gap**3
trapz = getattr(np, "trapezoid", None) or np.trapz
I = trapz(integrand, s)
print(f"f=0: Delta_min={gap.min():.3f}  |B-A|={nBA:.2f}  I={I:.1f}  T_b(1e-2)={I/np.sqrt(1e-2):.3g}  "
      f"peak at s={s[np.argmax(integrand)]:.2f}  gap(s=1)={gap[-1]:.3f}")

fig, (axL, axR) = plt.subplots(1, 2, figsize=(10, 3.6))
axL.plot(s, gap, color="C0")
axL.axhline(MU_BUF, ls=":", color="gray", lw=1)
axL.set_xlabel("s = t/T"); axL.set_ylabel(r"instantaneous gap $\Delta(s)$")
axL.set_title(r"gap along the $2\to4$ path (pairing model, $g=1$)")
axR.semilogy(s, integrand, color="C0")
axR.set_xlabel("s = t/T"); axR.set_ylabel(r"$\|\ddot H\|/\Delta^2 + 7\|\dot H\|^2/\Delta^3$")
axR.set_title(rf"integrand of the time bound ($\int_0^1 = {I:.0f}$)")
plt.tight_layout()
out = os.path.join(HERE, "..", "paper", "figs", "fig_bound_pairing.png")
plt.savefig(out, dpi=150); print("wrote", os.path.normpath(out))
