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
from .phmodel import ph_term_sector, ph_term_full
from .pauli import pauli_decompose

__all__ = ["pairing_terms", "ph_terms", "model_terms", "trotter_U", "trotter_U_dt",
           "trotter_error", "refine_state_trotter", "rodeo_track_U", "cnot_counts",
           "ph_pauli_costs"]


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


def ph_terms(k, f, N):
    """The genuinely new pieces of the particle-hole term V_ph, as a list of
    dense sector matrices, one per (p, q, r) with q != r: 2k(k-1) one-particle-
    one-hole terms (q = p or r = p) followed by k(k-1)(k-2) pair-breaking
    2p-2h terms.  The q = r pieces are pairing terms (g -> g + 2f) and belong
    to `pairing_terms`.  Every piece is exactly exponentiable by commuting
    Pauli rotations, so a product formula over `model_terms` has its error in
    the non-commutativity BETWEEN terms only, as on hardware."""
    nq, states, index = build_sector(k, N)
    one, two = [], []
    for p in range(k):
        for q in range(k):
            for r in range(k):
                if q == r:
                    continue
                T = ph_term_sector(k, f, N, states, index, p, q, r)
                (one if (q == p or r == p) else two).append(T)
    return one + two


def model_terms(k, g, N, f=0.0, delta=1.0):
    """Term list [D, h_pq..., V_ph pieces...] whose sum is H_pairing_sparse(f=f):
    the diagonal term and pair-hopping terms at the renormalised strength
    g + 2f, then the 1p-1h and pair-breaking pieces of V_ph (f != 0 only)."""
    D, hops = pairing_terms(k, g + 2.0*f, N, delta)
    terms = [D] + hops
    if f != 0.0:
        terms += ph_terms(k, f, N)
    return terms


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


def trotter_U(k, g, N, t, n_steps, order=2, delta=1.0, terms=None, f=0.0):
    """Product-formula approximation to exp(-i H t) in the fixed-N sector
    (f != 0 adds the particle-hole pieces, see `model_terms`)."""
    if terms is None:
        terms = model_terms(k, g, N, f, delta)
    Ustep = _step_U(terms, t / n_steps, order)
    U = np.linalg.matrix_power(Ustep, n_steps)
    return U


def trotter_U_dt(k, g, N, t, dt_max, order=2, delta=1.0, terms=None, f=0.0):
    """Product formula for exp(-i H t) with a FIXED step size: the number of
    steps is ceil(t/dt_max), as on hardware where the per-step depth is fixed
    and longer evolutions use more steps.  Returns (U, n_steps)."""
    n_steps = max(1, int(np.ceil(abs(t) / dt_max)))
    return trotter_U(k, g, N, t, n_steps, order, delta, terms, f), n_steps


def trotter_error(k, g, N, t, n_steps, order=2, delta=1.0, f=0.0, terms=None):
    """Spectral-norm distance between the product formula and exp(-i H t)."""
    nq, states, index = build_sector(k, N)
    H = H_pairing_sparse(k, g, N, states, index, delta, f=f).toarray()
    Uex = expm(-1j * H * t)
    Utr = trotter_U(k, g, N, t, n_steps, order, delta, terms, f)
    return float(np.linalg.norm(Utr - Uex, 2))


def refine_state_trotter(N, k_high, T, n_s, k_low=2, mu_buf=0.6, g=1.0,
                         order=2, delta=1.0, psi_low=None, f=0.0):
    """Adiabatic refinement with n_s product-formula steps (the circuit-level
    counterpart of refine_state, which uses exact slice exponentials).
    Each step splits H(s) into the embedded coarse term, the diagonal term and
    the pair-hopping terms, all realisable as the gate blocks of the paper.
    Returns (psi, H_high)."""
    nl, sl, il = build_sector(k_low, N)
    Hl = H_pairing_sparse(k_low, g, N, sl, il, delta, f=f).toarray()
    wl, vl = np.linalg.eigh(Hl)
    El = wl[0]
    if psi_low is None:
        psi_low = vl[:, 0]
    Pm = prolong_matrix(k_low, k_high, N)
    Phi0 = Pm @ psi_low
    nh, sh, ih = build_sector(k_high, N)
    Hh = H_pairing_sparse(k_high, g, N, sh, ih, delta, f=f).toarray()
    tm = model_terms(k_high, g, N, f, delta)
    D, rest = tm[0], tm[1:]
    mu = El + mu_buf
    Hle = Pm @ (Hl - mu * np.eye(len(Hl))) @ Pm.T
    Ds = D - mu * np.eye(len(D))
    base = [Hle, Ds] + rest                    # fixed matrices; only the s-dependent
    eigs = _eigterms(base)                     # prefactors change along the path
    psi = Phi0.astype(complex).copy()
    dt = T / n_s
    for m in range(n_s):
        s = (m + 0.5) / n_s
        c2, s2 = np.cos(np.pi * s / 2) ** 2, np.sin(np.pi * s / 2) ** 2
        coef = [c2, s2] + [s2] * len(rest)
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


_PH_COST_CACHE = {}


def ph_pauli_costs(k):
    """Pauli content of the q != r pieces of V_ph on 2k qubits (independent of
    f): for the 1p-1h and the pair-breaking 2p-2h families, the number of
    terms, of Pauli strings (= rotations), the CNOT count sum_strings 2(w-1)
    of their ladders, and the largest weight.  Computed once per k from the
    full-space Pauli decomposition of every piece."""
    if k in _PH_COST_CACHE:
        return _PH_COST_CACHE[k]
    out = {}
    for p in range(k):
        for q in range(k):
            for r in range(k):
                if q == r:
                    continue
                fam = "1p1h" if (q == p or r == p) else "2p2h"
                d = out.setdefault(fam, dict(terms=0, rotations=0, cnots=0, max_weight=0))
                dec = pauli_decompose(ph_term_full(k, 1.0, p, q, r))
                ws = [sum(ch != "I" for ch in s) for s in dec]
                d["terms"] += 1; d["rotations"] += len(ws)
                d["cnots"] += sum(2 * (w - 1) for w in ws)
                d["max_weight"] = max(d["max_weight"], max(ws))
    _PH_COST_CACHE[k] = out
    return out


def cnot_counts(k, f=0.0):
    """CNOT bookkeeping for one first-order Trotter step on 2k qubits (paper
    Eq. for the step cost): k(k-1)/2 pair-hopping terms of 8 weight-4
    rotations (48 CNOTs each), k ZZ rotations (2 CNOTs), 2k phases.  With
    f != 0 the particle-hole pieces add their rotations and CNOT ladders
    (`ph_pauli_costs`; the q = r pieces only rescale the pairing angles).
    Controlling the step for the rodeo cycle adds 2 CNOTs per rotation."""
    n_hop_rot = 8 * (k * (k - 1) // 2)
    n_rot = n_hop_rot + k + 2 * k          # + ZZ rotations + phases
    uncontrolled = 48 * (k * (k - 1) // 2) + 2 * k
    if f != 0.0:
        for d in ph_pauli_costs(k).values():
            n_rot += d["rotations"]; uncontrolled += d["cnots"]
    controlled = uncontrolled + 2 * n_rot
    return dict(step_uncontrolled=uncontrolled, step_controlled=controlled,
                rotations=n_rot)
