"""Trotterised evolution for the pairing Hamiltonian in the fixed-N sector.

The term decomposition mirrors the gate compilation used in the paper: a
diagonal part D (kinetic energy plus the p = q pairing term, compiled from
single-qubit phases and ZZ rotations) and one pair-hopping term h_pq per level
pair p < q (compiled as 8 mutually commuting weight-4 Pauli rotations, so each
factor exp(-i h_pq dt) is realised *exactly* by its circuit block).  The
Trotter error therefore comes only from the non-commutativity *between* terms,
exactly as it would on hardware.
"""
import numpy as np
from scipy.linalg import expm
from .hamiltonian import build_sector, H_pairing_sparse
from ._pairlib import _bit, _flip, _jw_sign
from .refine import prolong_matrix

__all__ = ["pairing_terms", "trotter_U", "trotter_U_dt", "trotter_error",
           "refine_state_trotter", "rodeo_track_U", "cnot_counts"]


def pairing_terms(k, g, N, delta=1.0):
    """Return (D, hops): diagonal term (dense matrix) and the list of
    pair-hopping terms h_pq (p < q), all in the fixed-N sector basis.
    Their sum equals H_pairing_sparse exactly."""
    nq, states, index = build_sector(k, N)
    M = len(states)
    D = np.zeros((M, M))
    for I in states:
        a = index[I]
        kin = sum(delta * (j // 2) for j in range(nq) if _bit(I, j, nq))
        ndouble = sum(1 for p in range(k)
                      if _bit(I, 2 * p, nq) and _bit(I, 2 * p + 1, nq))
        D[a, a] = kin - 0.5 * g * ndouble
    hops = []
    for p in range(k):
        for q in range(p + 1, k):
            h = np.zeros((M, M))
            for (pp, qq) in ((p, q), (q, p)):        # A+_pp A_qq and h.c.
                r0, r1 = 2 * qq, 2 * qq + 1
                c0, c1 = 2 * pp, 2 * pp + 1
                for I in states:
                    if not _bit(I, r0, nq):
                        continue
                    s1 = _jw_sign(I, r0, nq); t = _flip(I, r0, nq)
                    if not _bit(t, r1, nq):
                        continue
                    s2 = _jw_sign(t, r1, nq); t = _flip(t, r1, nq)
                    if _bit(t, c1, nq):
                        continue
                    s3 = _jw_sign(t, c1, nq); t = _flip(t, c1, nq)
                    if _bit(t, c0, nq):
                        continue
                    s4 = _jw_sign(t, c0, nq); J = _flip(t, c0, nq)
                    h[index[J], index[I]] += -0.5 * g * s1 * s2 * s3 * s4
            hops.append(0.5 * (h + h.T))
    return D, hops


def _eigterms(terms):
    """Eigendecompose each Hermitian term once, so every exponential factor
    exp(-i c T dt) is two matrix products instead of a fresh expm."""
    return [np.linalg.eigh(T) for T in terms]


def _factor(eig, c):
    lam, V = eig
    return (V * np.exp(-1j * lam * c)) @ V.conj().T


def _step_U(terms, dt, order=2, eigs=None):
    """One product-formula step exp(-i sum_i T_i dt) from term factors."""
    if eigs is None:
        eigs = _eigterms(terms)
    U = np.eye(terms[0].shape[0], dtype=complex)
    if order == 1:
        for e in eigs:
            U = _factor(e, dt) @ U
        return U
    # order 2 (symmetric): forward half-steps then backward half-steps
    for e in eigs:
        U = _factor(e, dt / 2) @ U
    for e in reversed(eigs):
        U = _factor(e, dt / 2) @ U
    return U


def trotter_U(k, g, N, t, n_steps, order=2, delta=1.0, terms=None):
    """Product-formula approximation to exp(-i H t) in the fixed-N sector."""
    if terms is None:
        D, hops = pairing_terms(k, g, N, delta)
        terms = [D] + hops
    Ustep = _step_U(terms, t / n_steps, order)
    U = np.linalg.matrix_power(Ustep, n_steps)
    return U


def trotter_U_dt(k, g, N, t, dt_max, order=2, delta=1.0, terms=None):
    """Product formula for exp(-i H t) with a FIXED step size: the number of
    steps is ceil(t/dt_max), as on hardware where the per-step depth is fixed
    and longer evolutions use more steps.  Returns (U, n_steps)."""
    n_steps = max(1, int(np.ceil(abs(t) / dt_max)))
    return trotter_U(k, g, N, t, n_steps, order, delta, terms), n_steps


def trotter_error(k, g, N, t, n_steps, order=2, delta=1.0):
    """Spectral-norm distance between the product formula and exp(-i H t)."""
    nq, states, index = build_sector(k, N)
    H = H_pairing_sparse(k, g, N, states, index, delta).toarray()
    Uex = expm(-1j * H * t)
    Utr = trotter_U(k, g, N, t, n_steps, order, delta)
    return float(np.linalg.norm(Utr - Uex, 2))


def refine_state_trotter(N, k_high, T, n_s, k_low=2, mu_buf=0.6, g=1.0,
                         order=2, delta=1.0, psi_low=None):
    """Adiabatic refinement with n_s product-formula steps (the circuit-level
    counterpart of refine_state, which uses exact slice exponentials).
    Each step splits H(s) into the embedded coarse term, the diagonal term and
    the pair-hopping terms, all realisable as the gate blocks of the paper.
    Returns (psi, H_high)."""
    nl, sl, il = build_sector(k_low, N)
    Hl = H_pairing_sparse(k_low, g, N, sl, il, delta).toarray()
    wl, vl = np.linalg.eigh(Hl)
    El = wl[0]
    if psi_low is None:
        psi_low = vl[:, 0]
    Pm = prolong_matrix(k_low, k_high, N)
    Phi0 = Pm @ psi_low
    nh, sh, ih = build_sector(k_high, N)
    Hh = H_pairing_sparse(k_high, g, N, sh, ih, delta).toarray()
    D, hops = pairing_terms(k_high, g, N, delta)
    mu = El + mu_buf
    Hle = Pm @ (Hl - mu * np.eye(len(Hl))) @ Pm.T
    Ds = D - mu * np.eye(len(D))
    base = [Hle, Ds] + hops                    # fixed matrices; only the s-dependent
    eigs = _eigterms(base)                     # prefactors change along the path
    psi = Phi0.astype(complex).copy()
    dt = T / n_s
    for m in range(n_s):
        s = (m + 0.5) / n_s
        c2, s2 = np.cos(np.pi * s / 2) ** 2, np.sin(np.pi * s / 2) ** 2
        coef = [c2, s2] + [s2] * len(hops)
        if order == 1:
            for e, c in zip(eigs, coef):
                psi = _factor(e, c * dt) @ psi
        else:
            for e, c in zip(eigs, coef):
                psi = _factor(e, c * dt / 2) @ psi
            for e, c in zip(reversed(eigs), reversed(coef)):
                psi = _factor(e, c * dt / 2) @ psi
    return psi / np.linalg.norm(psi), Hh


def rodeo_track_U(v0, Ufun, ts, E, target):
    """Rodeo fidelity/acceptance tracking with a user-supplied evolution
    Ufun(t) (e.g. a Trotterised controlled evolution)."""
    v = v0.astype(complex).copy()
    P = 1.0
    F, A = [], []
    for t in ts:
        w = 0.5 * (v + np.exp(1j * E * t) * (Ufun(t) @ v))
        p = float(np.vdot(w, w).real)
        v = w / np.sqrt(p)
        P *= p
        F.append(abs(np.vdot(target, v)) ** 2)
        A.append(P)
    return np.array(F), np.array(A)


def cnot_counts(k):
    """CNOT bookkeeping for one first-order Trotter step of the pairing H on
    2k qubits (paper Eq. for the step cost): k(k-1)/2 pair-hopping terms of 8
    weight-4 rotations (48 CNOTs each), k ZZ rotations (2 CNOTs), 2k phases.
    Controlling the step for the rodeo cycle adds 2 CNOTs per rotation."""
    n_hop_rot = 8 * (k * (k - 1) // 2)
    n_rot = n_hop_rot + k + 2 * k          # + ZZ rotations + phases
    uncontrolled = 48 * (k * (k - 1) // 2) + 2 * k
    controlled = uncontrolled + 2 * n_rot
    return dict(step_uncontrolled=uncontrolled, step_controlled=controlled,
                rotations=n_rot)
