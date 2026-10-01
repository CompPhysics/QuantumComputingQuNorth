"""Constant-pairing Hamiltonian, encoding, exact diagonalisation."""
import numpy as np
from scipy.sparse import csr_matrix
from ._pairlib import (_bit, _flip, _jw_sign, build_sector,
                       H_pairing_sparse as _H_pair_only, fci_ground, E_HF)
from .phmodel import H_ph_sparse, build_H_ph_full

__all__ = ["build_sector", "H_pairing_sparse", "fci_ground", "E_HF",
           "build_H_full", "_bit", "_flip", "_jw_sign"]

def H_pairing_sparse(k, g, N, states, index, delta=1.0, f=0.0):
    """Sparse Hamiltonian in the fixed-N sector.  f = 0 (default) is the
    constant-pairing model; f != 0 adds the pair-breaking particle-hole term
    V_ph of `phmodel` (g -> g + 2f in the seniority-zero block, plus 1p-1h and
    seniority-two 2p-2h couplings)."""
    H = _H_pair_only(k, g, N, states, index, delta)
    if f != 0.0: H = H + H_ph_sparse(k, f, N, states, index)
    return H

def build_H_full(k, g, delta=1.0, f=0.0):
    """H over the FULL 2^(2k) Fock space (used by HEA / gate-UCCSD).
    f != 0 adds the particle-hole term V_ph (see `phmodel`)."""
    if f != 0.0: return _build_H_full_pair(k, g, delta) + build_H_ph_full(k, f)
    return _build_H_full_pair(k, g, delta)

def _build_H_full_pair(k, g, delta=1.0):
    nq = 2*k; dim = 2**nq; gh = -0.5*g; rows, cols, vals = [], [], []
    for I in range(dim):
        rows.append(I); cols.append(I)
        vals.append(sum(delta*(j//2) for j in range(nq) if _bit(I, j, nq)))
        for p in range(k):
            for q in range(k):
                r0, r1 = 2*q, 2*q+1; c0, c1 = 2*p, 2*p+1
                if not _bit(I, r0, nq): continue
                s1 = _jw_sign(I, r0, nq); t = _flip(I, r0, nq)
                if not _bit(t, r1, nq): continue
                s2 = _jw_sign(t, r1, nq); t = _flip(t, r1, nq)
                if _bit(t, c1, nq): continue
                s3 = _jw_sign(t, c1, nq); t = _flip(t, c1, nq)
                if _bit(t, c0, nq): continue
                s4 = _jw_sign(t, c0, nq); J = _flip(t, c0, nq)
                rows.append(J); cols.append(I); vals.append(gh*s1*s2*s3*s4)
    H = csr_matrix((vals, (rows, cols)), shape=(dim, dim))
    return ((H + H.T)*0.5).toarray().real
