"""The pairing model with a pair-breaking particle-hole interaction.

    H = H_0 + V_pair + V_ph,
    V_ph = -(f/2) sum_{pqr} ( a+_{p up} a+_{p dn} a_{q dn} a_{r up} + h.c. ),

with H_0 and V_pair exactly the constant-pairing Hamiltonian of
`hamiltonian.py` (Hjorth-Jensen, *Quantum mechanics for many-particle
systems*, Chapter 4, Eq. 4-Vph).  The q = r pieces renormalise g -> g + 2f; the
q = p != r pieces move a single particle (1p-1h, seniority 2); the p, q, r all
different pieces deposit a pair while removing two unpaired particles (2p-2h,
seniority 2).  V_ph conserves N_up and N_dn separately but NOT the seniority,
so the ground state leaves the seniority-zero subspace as soon as f != 0.

Conventions are those of the rest of the package: qubit 2p carries spin-up of
level p, qubit 2p+1 spin-down, big-endian bit order, Jordan-Wigner signs.

Validated anchors (k = 4 levels, N = 4, delta = 1, g = 1; book Table 4.x):
    f = 0.00 -> E_0 =  0.63554847
    f = 0.05 -> E_0 =  0.45058234
    f = 0.20 -> E_0 = -0.18348455
    f = 0.50 -> E_0 = -1.69173670
and E_ref = 2 - g - 2f = 1 - 2f for the reference determinant.
"""
import numpy as np
from scipy.sparse import csr_matrix
from ._pairlib import _bit, _flip, _jw_sign, build_sector, E_HF

__all__ = ["H_ph_sparse", "build_H_ph_full", "E_ref", "seniority_of",
           "seniority_weight", "seniority_zero_mask", "ph_coupling_table"]


def _apply_ph_string(I, nq, p, q, r):
    """Apply a+_{p up} a+_{p dn} a_{q dn} a_{r up} to determinant I.
    Returns (J, sign) or (None, 0)."""
    o1, o2, o3, o4 = 2*r, 2*q + 1, 2*p + 1, 2*p      # annihilate, annihilate, create, create
    if not _bit(I, o1, nq): return None, 0
    s1 = _jw_sign(I, o1, nq); t = _flip(I, o1, nq)
    if not _bit(t, o2, nq): return None, 0
    s2 = _jw_sign(t, o2, nq); t = _flip(t, o2, nq)
    if _bit(t, o3, nq): return None, 0
    s3 = _jw_sign(t, o3, nq); t = _flip(t, o3, nq)
    if _bit(t, o4, nq): return None, 0
    s4 = _jw_sign(t, o4, nq); J = _flip(t, o4, nq)
    return J, s1*s2*s3*s4


def H_ph_sparse(k, f, N, states, index):
    """Sparse matrix of V_ph (Hermitian, including the h.c. term) in the
    fixed-N sector spanned by `states`."""
    nq = 2*k; M = len(states); rows, cols, vals = [], [], []
    for I in states:
        a = index[I]
        for p in range(k):
            for q in range(k):
                for r in range(k):
                    J, s = _apply_ph_string(I, nq, p, q, r)
                    if J is None: continue
                    rows.append(index[J]); cols.append(a); vals.append(-0.5*f*s)
    T = csr_matrix((vals, (rows, cols)), shape=(M, M))
    return T + T.T            # the Hermitian conjugate is a separate sum


def build_H_ph_full(k, f):
    """Dense V_ph over the FULL 2^(2k) Fock space (for gate-level circuits)."""
    nq = 2*k; dim = 2**nq; rows, cols, vals = [], [], []
    for I in range(dim):
        for p in range(k):
            for q in range(k):
                for r in range(k):
                    J, s = _apply_ph_string(I, nq, p, q, r)
                    if J is None: continue
                    rows.append(J); cols.append(I); vals.append(-0.5*f*s)
    T = csr_matrix((vals, (rows, cols)), shape=(dim, dim))
    return (T + T.T).toarray().real


def E_ref(N, g, f=0.0, delta=1.0):
    """Energy of the reference determinant (lowest N/2 levels doubly filled).
    V_ph contributes only through its q = r pieces, g -> g + 2f, so
    E_ref = E_HF(N, g) - f N/2  (= 1 - 2f for N = 4, g = delta = 1)."""
    return E_HF(N, g, delta) - f*(N//2)


def seniority_of(I, k):
    """Number of particles in determinant I that are not part of a complete pair."""
    nq = 2*k; s = 0
    for p in range(k):
        up, dn = _bit(I, 2*p, nq), _bit(I, 2*p + 1, nq)
        s += (up != dn)
    return s


def seniority_zero_mask(k, states):
    """Boolean mask over `states` selecting the seniority-zero determinants."""
    return np.array([seniority_of(I, k) == 0 for I in states])


def seniority_weight(psi, k, states):
    """Weight sum |c_I|^2 of a sector state vector on the seniority-zero
    determinants (1 for any pure-pairing eigenstate)."""
    m = seniority_zero_mask(k, states)
    psi = np.asarray(psi); w = np.abs(psi)**2
    return float(w[m].sum()/w.sum())


def ph_coupling_table(k, N, g, f, delta=1.0):
    """Determinants connected to the reference by ONE action of H, classified
    by excitation rank (number of particle-hole pairs) and seniority, with the
    matrix element.  Reproduces the book's Table 'pairingphcoupling'.
    Returns a dict {(rank, seniority): (count, sorted unique |matrix elements|)}."""
    from .hamiltonian import H_pairing_sparse
    nq, states, index = build_sector(k, N)
    H = H_pairing_sparse(k, g, N, states, index, delta=delta, f=f).tocsc()
    ref = sum(1 << (nq - 1 - j) for j in range(N)); a = index[ref]
    col = H[:, a].toarray().ravel()
    out = {}
    for b, hval in enumerate(col):
        if b == a or abs(hval) < 1e-14: continue
        I = states[b]
        holes = sum(1 for j in range(N) if not _bit(I, j, nq))
        key = (holes, seniority_of(I, k))
        cnt, vals = out.get(key, (0, []))
        out[key] = (cnt + 1, vals + [abs(hval)])
    return {key: (cnt, sorted({float(round(x, 10)) for x in v})) for key, (cnt, v) in out.items()}
