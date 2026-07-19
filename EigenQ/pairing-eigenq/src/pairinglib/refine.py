"""Prolongation operator and adiabatic resolution refinement.

The slice evolution exp(-i H(s_m) dt) is applied through the eigendecomposition
of H(s_m); since the s-grid is fixed, the decompositions are computed once and
reused for every total time T (identical results to scipy.linalg.expm, much
faster for repeated T scans)."""
import numpy as np
from scipy.linalg import expm
from .hamiltonian import build_sector, H_pairing_sparse, _bit


def _slice_eigs(Hle, Hhs, M):
    """Eigendecompositions of H(s_m) on the fixed midpoint s-grid."""
    out = []
    for m in range(M):
        s = (m + 0.5)/M
        H = np.cos(np.pi*s/2)**2*Hle + np.sin(np.pi*s/2)**2*Hhs
        out.append(np.linalg.eigh(H))
    return out


def _evolve(Phi0, eigs, T):
    psi = Phi0.astype(complex).copy()
    dt = T/len(eigs)
    for lam, V in eigs:
        psi = V @ (np.exp(-1j*lam*dt) * (V.conj().T @ psi))
    return psi

__all__ = ["prolong_matrix", "refine_state", "refine_NK", "refine_from_state"]

def prolong_matrix(k_low, k_high, N):
    """Basis-refinement embedding: new high-resolution levels start empty."""
    nl, sl, il = build_sector(k_low, N); nh, sh, ih = build_sector(k_high, N)
    Pm = np.zeros((len(sh), len(sl)))
    for Il in sl:
        Ih = 0
        for j in range(nl):
            if _bit(Il, j, nl): Ih |= 1 << (nh-1-j)
        Pm[ih[Ih], il[Il]] = 1.0
    return Pm

def refine_state(N, k_high, T, k_low=2, M=160, mu_buf=0.6, g=1.0):
    """Adiabatically refine the prolonged low-res ground state; return (psi, H_high)."""
    nl, sl, il = build_sector(k_low, N); Hl = H_pairing_sparse(k_low, g, N, sl, il).toarray()
    wl, vl = np.linalg.eigh(Hl); El = wl[0]; psi_low = vl[:, 0]
    Pm = prolong_matrix(k_low, k_high, N); Phi0 = Pm @ psi_low
    nh, sh, ih = build_sector(k_high, N); Hh = H_pairing_sparse(k_high, g, N, sh, ih).toarray()
    mu = El + mu_buf
    Hle = Pm @ (Hl - mu*np.eye(len(Hl))) @ Pm.T; Hhs = Hh - mu*np.eye(len(Hh))
    psi = _evolve(Phi0, _slice_eigs(Hle, Hhs, M), T)
    return psi/np.linalg.norm(psi), Hh

def refine_NK(N, k_high, Ts, k_low=2, M=160, mu_buf=0.6, g=1.0):
    """Overlap and physical energy of the refined state vs adiabatic time T."""
    nl, sl, il = build_sector(k_low, N); Hl = H_pairing_sparse(k_low, g, N, sl, il).toarray()
    wl, vl = np.linalg.eigh(Hl); El = wl[0]; psi_low = vl[:, 0]
    Pm = prolong_matrix(k_low, k_high, N); Phi0 = Pm @ psi_low
    nh, sh, ih = build_sector(k_high, N); Hh = H_pairing_sparse(k_high, g, N, sh, ih).toarray()
    wh, vh = np.linalg.eigh(Hh); Eh = wh[0]; Psih = vh[:, 0]
    mu = El + mu_buf
    Hle = Pm @ (Hl - mu*np.eye(len(Hl))) @ Pm.T; Hhs = Hh - mu*np.eye(len(Hh))
    eigs = _slice_eigs(Hle, Hhs, M)
    ov, en = [], []
    for T in Ts:
        psi = _evolve(Phi0, eigs, T)
        ov.append(abs(np.vdot(Psih, psi))**2)
        en.append(float(np.real(psi.conj() @ (Hh @ psi))))
    return dict(ov=np.array(ov), en=np.array(en), Eh=Eh, El=El, gap=wh[1]-wh[0],
                ov0=abs(np.vdot(Psih, Phi0))**2, E0=float(np.real(Phi0 @ Hh @ Phi0)))

def refine_from_state(psi_low, N, k_low, k_high, Ts, M=160, mu_buf=0.6, g=1.0):
    """Refinement starting from a SUPPLIED coarse state (e.g. a gate-level
    UCCSD-VQE state at k_low > 2), rather than the exact coarse ground state.
    Returns the same dictionary as refine_NK.  The shift mu uses the coarse
    variational energy <psi_low|H_low|psi_low>, which is all a device knows."""
    nl, sl, il = build_sector(k_low, N); Hl = H_pairing_sparse(k_low, g, N, sl, il).toarray()
    psi_low = np.asarray(psi_low, dtype=complex); psi_low = psi_low/np.linalg.norm(psi_low)
    El = float(np.real(psi_low.conj() @ (Hl @ psi_low)))
    Pm = prolong_matrix(k_low, k_high, N); Phi0 = Pm @ psi_low
    nh, sh, ih = build_sector(k_high, N); Hh = H_pairing_sparse(k_high, g, N, sh, ih).toarray()
    wh, vh = np.linalg.eigh(Hh); Eh = wh[0]; Psih = vh[:, 0]
    mu = El + mu_buf
    Hle = Pm @ (Hl - mu*np.eye(len(Hl))) @ Pm.T; Hhs = Hh - mu*np.eye(len(Hh))
    eigs = _slice_eigs(Hle, Hhs, M)
    ov, en = [], []
    for T in Ts:
        psi = _evolve(Phi0, eigs, T)
        ov.append(abs(np.vdot(Psih, psi))**2)
        en.append(float(np.real(psi.conj() @ (Hh @ psi))))
    return dict(ov=np.array(ov), en=np.array(en), Eh=Eh, El=El, gap=wh[1]-wh[0],
                ov0=abs(np.vdot(Psih, Phi0))**2,
                E0=float(np.real(Phi0.conj() @ (Hh @ Phi0))))
