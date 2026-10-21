"""Seniority-zero (pair-basis) implementation of the constant-pairing pipeline.

The pairing Hamiltonian, the basis-refinement prolongation, the refinement
path, the pair-hopping Trotter terms and the rodeo filter all conserve
seniority, so for an N-paired target the whole pipeline of paper I lives in
the seniority-zero block: P = N/2 pairs distributed over k doubly degenerate
levels, dimension C(k, P).  This module builds that block directly in the
pair-occupation basis (one bit per level), which makes model spaces such as
k = 8, P = 4 (70 states instead of the 12 870 determinants of the full
N = 8 sector on 16 qubits) cheap.  Every object is the exact restriction of
its full-sector counterpart in `pairinglib` (verified in tests/test_seniority.py):
the pair-hopping terms h_pq are the seniority-zero blocks of the 48-CNOT
gate blocks of paper I, so the product-formula dynamics is identical, and the
CNOT bookkeeping (`cnot_counts`, `coarse_term_cnots`) is that of the full
2k-qubit Jordan-Wigner compilation of paper I unless the pair encoding
(one qubit per level) is requested explicitly.

Conventions as in `pairinglib`: level p (0-based) has energy delta*p per
particle, H = sum_p 2 delta p n_p - (g/2) sum_{pq} P_p^+ P_q, the appended
levels of a refinement k_low -> k_high are the highest ones and start empty.
"""
import numpy as np
from itertools import combinations
from scipy.optimize import minimize
from .trotter import cnot_counts

__all__ = ["pair_basis", "H_pair", "pair_terms", "prolong_pair", "hf_state",
           "seniority_zero_states", "refine_pair", "refine_pair_overlaps",
           "pair_double_ansatz", "coarse_term_cnots", "pair_encoding_cnots",
           "ramp_terms"]


# ------------------------------------------------------------------ basis and Hamiltonian
def pair_basis(k, P):
    """Pair-occupation basis: bit p of a state is 1 if level p holds a pair.
    Returns (states, index); states are sorted integers."""
    states = sorted(sum(1 << p for p in c) for c in combinations(range(k), P))
    return states, {s: i for i, s in enumerate(states)}


def pair_terms(k, g, P, delta=1.0):
    """[D, h_pq (p < q) ...] in the pair basis: D = sum_p (2 delta p - g/2) n_p,
    h_pq = -(g/2)(P_p^+ P_q + h.c.).  Their sum is H_pair.  Pair operators
    carry no Jordan-Wigner sign (both fermions of a pair cross the same
    occupied orbitals), so every hopping element is -g/2."""
    states, index = pair_basis(k, P)
    dim = len(states)
    D = np.zeros((dim, dim))
    for s in states:
        D[index[s], index[s]] = sum(2.0 * delta * p - 0.5 * g for p in range(k) if (s >> p) & 1)
    hops = []
    for p in range(k):
        for q in range(p + 1, k):
            h = np.zeros((dim, dim))
            for s in states:
                if (s >> q) & 1 and not (s >> p) & 1:          # pair q -> p
                    t = s ^ (1 << q) ^ (1 << p)
                    h[index[t], index[s]] = h[index[s], index[t]] = -0.5 * g
            hops.append(h)
    return D, hops


def H_pair(k, g, P, delta=1.0):
    """Dense seniority-zero Hamiltonian, dimension C(k, P)."""
    D, hops = pair_terms(k, g, P, delta)
    return D + sum(hops) if hops else D


def hf_state(k, P):
    """The Hartree-Fock determinant (lowest P levels paired) in the pair basis."""
    states, index = pair_basis(k, P)
    v = np.zeros(len(states)); v[index[(1 << P) - 1]] = 1.0
    return v


def prolong_pair(k_low, k_high, P):
    """Basis-refinement embedding C(k_low,P) -> C(k_high,P): the appended levels
    k_low..k_high-1 start empty (same convention as `prolong_matrix`)."""
    sl, il = pair_basis(k_low, P); sh, ih = pair_basis(k_high, P)
    Pm = np.zeros((len(sh), len(sl)))
    for s in sl:
        Pm[ih[s], il[s]] = 1.0
    return Pm


def seniority_zero_states(k, N):
    """Map between the pair basis and the full N-sector of `build_sector`:
    returns (full_indices, pair_indices) such that full-sector vector
    components full_indices correspond to pair-basis components pair_indices."""
    from .hamiltonian import build_sector
    nq, st, ix = build_sector(k, N)
    ps, pix = pair_basis(k, N // 2)
    full, pair = [], []
    for I in st:
        occ = [(I >> (nq - 1 - j)) & 1 for j in range(nq)]
        if all(occ[2 * p] == occ[2 * p + 1] for p in range(k)):
            s = sum(1 << p for p in range(k) if occ[2 * p])
            full.append(ix[I]); pair.append(pix[s])
    return np.array(full), np.array(pair)


# ------------------------------------------------------------------ refinement (exact slices)
def _slice_eigs(Hle, Hhs, M):
    out = []
    for m in range(M):
        s = (m + 0.5) / M
        out.append(np.linalg.eigh(np.cos(np.pi * s / 2) ** 2 * Hle + np.sin(np.pi * s / 2) ** 2 * Hhs))
    return out


def _evolve(Phi0, eigs, T):
    psi = Phi0.astype(complex).copy(); dt = T / len(eigs)
    for lam, V in eigs:
        psi = V @ (np.exp(-1j * lam * dt) * (V.conj().T @ psi))
    return psi


def refine_pair(k_low, k_high, P, g, mu_buf=0.6, M=160, psi_low=None, delta=1.0):
    """Set up the shift-controlled projector path of paper I in the pair basis.
    Returns dict(Phi0, Hle, Hhs, eigs, H_high, w, V, E_low, mu, p0, Psih):
    `eigs` are the M midpoint-slice eigendecompositions, reusable for any T.
    With psi_low given, the coarse state is that vector and mu uses its
    variational coarse energy, as in `refine_from_state`."""
    Hl = H_pair(k_low, g, P, delta); Hh = H_pair(k_high, g, P, delta)
    if psi_low is None:
        wl, Vl = np.linalg.eigh(Hl); El = wl[0]; psi_low = Vl[:, 0]
    else:
        psi_low = np.asarray(psi_low, complex) / np.linalg.norm(psi_low)
        El = float(np.real(psi_low.conj() @ Hl @ psi_low))
    Pm = prolong_pair(k_low, k_high, P); Phi0 = Pm @ psi_low
    w, V = np.linalg.eigh(Hh); Psih = V[:, 0]
    mu = El + mu_buf
    Hle = Pm @ (Hl - mu * np.eye(len(Hl))) @ Pm.T; Hhs = Hh - mu * np.eye(len(Hh))
    return dict(Phi0=Phi0, Hle=Hle, Hhs=Hhs, eigs=_slice_eigs(Hle, Hhs, M), H_high=Hh, w=w, V=V,
                E_low=El, mu=mu, p0=abs(np.vdot(Psih, Phi0)) ** 2, Psih=Psih, Pm=Pm, Hl=Hl)


def refine_pair_overlaps(r, Ts):
    """Overlap p(T) with the fine ground state along the path of `refine_pair`."""
    return np.array([abs(np.vdot(r["Psih"], _evolve(r["Phi0"], r["eigs"], T))) ** 2 for T in Ts])


def fixed_step_evolution(terms, dt_max, order=2):
    """Ufun(t) = product formula for exp(-i sum terms t) with a FIXED maximum
    step dt_max: n = ceil(|t|/dt_max) steps of size t/n (as on hardware), each
    term exponentiated exactly from its eigendecomposition (done once)."""
    from .trotter import _eigterms, _step_U
    eigs = _eigterms(terms)
    def Ufun(t):
        n = max(1, int(np.ceil(abs(t) / dt_max)))
        return np.linalg.matrix_power(_step_U(terms, t / n, order, eigs), n)
    return Ufun


def ramp_terms(k, g, P, delta=1.0):
    """Groups and coefficients of the ordinary HF adiabatic ramp
    H(s) = D + sin^2(pi s/2) sum_pq h_pq (the baseline (v) of paper I)."""
    D, hops = pair_terms(k, g, P, delta)
    s2 = lambda s: np.sin(np.pi * s / 2) ** 2
    return [D] + hops, [lambda s: 1.0] + [s2] * len(hops)


# ------------------------------------------------------------------ restricted UCCSD (pair doubles)
def pair_double_ansatz(k, P, g, delta=1.0, maxiter=500, theta0=None, H=None):
    """The seniority-zero restricted UCCSD of paper I: the product of the
    P(k-P) occupied-pair -> virtual-pair rotations exp(theta (P_a^+ P_i - h.c.))
    applied to the HF determinant, in the order of the full ansatz (i outer,
    a inner), optimized with L-BFGS and the exact adjoint gradient.  Each
    factor is the exact seniority-zero block of a 48-CNOT double-excitation
    block.  Returns dict(E, theta, psi, n_params, cnots, pairs)."""
    states, index = pair_basis(k, P)
    dim = len(states)
    if H is None:
        H = H_pair(k, g, P, delta)
    pairs = [(i, a) for i in range(P) for a in range(P, k)]
    gens = []
    for (i, a) in pairs:
        Is, Js = [], []
        for s in states:
            if (s >> i) & 1 and not (s >> a) & 1:
                Is.append(index[s]); Js.append(index[s ^ (1 << i) ^ (1 << a)])
        gens.append((np.array(Is), np.array(Js)))

    def apply(th, psi, Is, Js):
        c, s = np.cos(th), np.sin(th)
        a, b = psi[Is].copy(), psi[Js].copy()
        out = psi.copy(); out[Is] = c * a - s * b; out[Js] = s * a + c * b
        return out

    def gen_apply(psi, Is, Js):                       # K = P_a^+ P_i - h.c.
        out = np.zeros_like(psi); out[Js] += psi[Is]; out[Is] -= psi[Js]
        return out

    hf = hf_state(k, P)

    def forward(th):
        st = [hf.copy()]
        for m, (Is, Js) in enumerate(gens):
            st.append(apply(th[m], st[-1], Is, Js))
        return st

    def energy(th):
        psi = forward(th)[-1]
        return float(psi @ H @ psi)

    def grad(th):
        st = forward(th); lam = H @ st[-1]; gg = np.zeros(len(th))
        for m in range(len(gens) - 1, -1, -1):
            Is, Js = gens[m]
            gg[m] = 2.0 * float(lam @ gen_apply(st[m + 1], Is, Js))
            lam = apply(-th[m], lam, Is, Js)
        return gg

    x0 = np.zeros(len(pairs)) if theta0 is None else np.asarray(theta0, float)
    res = minimize(energy, x0, jac=grad, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-15, "gtol": 1e-11})
    psi = forward(res.x)[-1]
    return dict(E=float(res.fun), theta=res.x, psi=psi / np.linalg.norm(psi), n_params=len(pairs),
                cnots=48 * len(pairs), pairs=pairs, nit=res.nit)


# ------------------------------------------------------------------ gate bookkeeping
def coarse_term_cnots(k_low, k_high, P, convention="gray"):
    """CNOTs per application of the embedded coarse term A = P(H_low-mu)P^T of a
    basis refinement k_low -> k_high in the 2k-qubit encoding of paper I,
    conditioned on the c = 2(k_high-k_low) appended qubits being empty.
    * Filled coarse model (P = k_low): A is a single conditioned phase.
      'gray': a c-qubit Gray-code phase polynomial, 2^c - 2 CNOTs (14 for the
      2->4 chain, 254 for 4->8); 'toffoli'/'flag': a clean-ancilla AND ladder of
      relative-phase Toffolis, 6(c-1) + 2 CNOTs (44 for 4->8).
    * Nontrivial coarse space: the k_low step with every rotation (single-qubit
      phases included, since they are no longer free) conditioned on the c
      qubits, PLUS the conditioned scalar phase of the constant part of
      H_low - mu, which is a physical relative phase after projection.
      'gray': 2^(c+1) - 2 extra CNOTs per rotation and 2^c - 2 for the scalar
      (350 for 3->4, 4 886 for 6->8); 'toffoli': a clean-ancilla ladder per
      rotation, 6(c-1) + 2 each; 'flag': the conjunction of the c controls is
      computed ONCE into a clean flag qubit, every rotation is singly
      controlled by it (2 CNOTs) and the flag is uncomputed,
      base + 2 n_rot + 6(c-1) + 2 (224 for 3->4, 1 028 for 6->8; needs c-1
      clean work qubits beyond the 2k+1 of the ancilla-free compilation)."""
    c = 2 * (k_high - k_low)
    if P == k_low:
        return 2 ** c - 2 if convention == "gray" else 6 * (c - 1) + 2
    base = cnot_counts(k_low); nrot = base["rotations"]; un = base["step_uncontrolled"]
    if convention == "gray":
        return un + nrot * (2 ** (c + 1) - 2) + (2 ** c - 2)
    if convention == "toffoli":
        return un + (nrot + 1) * (6 * (c - 1) + 2)
    if convention == "flag":
        return un + 2 * nrot + 6 * (c - 1) + 2
    raise ValueError(convention)


def pair_encoding_cnots(k):
    """CNOT bookkeeping of the one-qubit-per-level (pair) encoding, in which
    h_pq = -(g/4)(X_p X_q + Y_p Y_q): two weight-2 rotations (4 CNOTs) per
    level pair and single-qubit phases for D.  Same keys as `cnot_counts`."""
    nhop = k * (k - 1) // 2
    nrot = 2 * nhop + k
    un = 4 * nhop
    return dict(step_uncontrolled=un, step_controlled=un + 2 * nrot, rotations=nrot)


# ------------------------------------------------------------------ ADAPT-VQE in the pair basis
def _pair_generators(k, P, pairs):
    """Index arrays (Is, Js) of the pair-hopping generators K_{ia} = P_a^+ P_i - h.c.
    for the level pairs in `pairs` (states with a pair on i and none on a map to
    the states with the pair moved to a)."""
    states, index = pair_basis(k, P)
    gens = []
    for (i, a) in pairs:
        Is, Js = [], []
        for s_ in states:
            if (s_ >> i) & 1 and not (s_ >> a) & 1:
                Is.append(index[s_]); Js.append(index[s_ ^ (1 << i) ^ (1 << a)])
        gens.append((np.array(Is, int), np.array(Js, int)))
    return gens


def _apply_rot(th, psi, Is, Js):
    c, s_ = np.cos(th), np.sin(th)
    a, b = psi[Is].copy(), psi[Js].copy()
    out = psi.copy(); out[Is] = c * a - s_ * b; out[Js] = s_ * a + c * b
    return out


def _apply_gen(psi, Is, Js):
    out = np.zeros_like(psi); out[Js] += psi[Is]; out[Is] -= psi[Js]
    return out


def _optimize_product(H, hf, gens, theta0, maxiter=500):
    """Minimize <psi|H|psi> for psi = prod_m exp(theta_m K_m)|HF> (K_m applied in
    list order) with L-BFGS and the exact adjoint gradient."""
    def forward(th):
        st = [hf.copy()]
        for m, (Is, Js) in enumerate(gens):
            st.append(_apply_rot(th[m], st[-1], Is, Js))
        return st

    def energy(th):
        psi = forward(th)[-1]; return float(psi @ H @ psi)

    def grad(th):
        st = forward(th); lam = H @ st[-1]; gg = np.zeros(len(th))
        for m in range(len(gens) - 1, -1, -1):
            Is, Js = gens[m]
            gg[m] = 2.0 * float(lam @ _apply_gen(st[m + 1], Is, Js))
            lam = _apply_rot(-th[m], lam, Is, Js)
        return gg

    if len(gens) == 0:
        return np.zeros(0), float(hf @ H @ hf), hf.copy()
    res = minimize(energy, np.asarray(theta0, float), jac=grad, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-15, "gtol": 1e-11})
    psi = forward(res.x)[-1]
    return res.x, float(res.fun), psi / np.linalg.norm(psi)


def adapt_vqe(k, P, g, pool="generalized", grad_tol=1e-4, max_ops=None, delta=1.0, H=None,
              target=None, maxiter=500):
    """ADAPT-VQE (Grimsley et al. 2019) in the seniority-zero pair basis.
    Pool: 'pairs' = the P(k-P) occupied-pair -> virtual-pair generators of the
    restricted UCCSD; 'generalized' = all k(k-1)/2 pair-hopping generators
    P_q^+ P_p - h.c. (p < q), the seniority-zero analog of a generalized-doubles
    pool.  At every iteration the pool operator with the largest energy gradient
    |<psi|[H, K]|psi>| is appended (an operator may be selected again) and all
    amplitudes are re-optimized from the previous optimum.  Stops when the
    largest gradient falls below grad_tol or after max_ops operators.  Returns
    dict(history=[dict(n_ops, E, infidelity, max_grad, op, cnots, theta, psi)],
    theta, ops, psi, E, pool); cnots counts 48 CNOTs per pair rotation (2k-qubit
    encoding); history[n] holds the ADAPT state after n operators."""
    if H is None:
        H = H_pair(k, g, P, delta)
    if target is None:
        target = np.linalg.eigh(H)[1][:, 0]
    if pool == "pairs":
        pairs = [(i, a) for i in range(P) for a in range(P, k)]
    elif pool == "generalized":
        pairs = [(p, q) for p in range(k) for q in range(p + 1, k)]
    else:
        raise ValueError(pool)
    pool_gens = _pair_generators(k, P, pairs)
    hf = hf_state(k, P)
    if max_ops is None:
        max_ops = 4 * len(pairs)
    ops, theta, psi = [], np.zeros(0), hf.copy()
    E = float(hf @ H @ hf)
    hist = [dict(n_ops=0, E=E, infidelity=1 - abs(np.vdot(target, psi)) ** 2, max_grad=None, op=None, cnots=0)]
    while len(ops) < max_ops:
        Hpsi = H @ psi
        grads = np.array([2.0 * float(Hpsi @ _apply_gen(psi, Is, Js)) for (Is, Js) in pool_gens])
        m = int(np.argmax(np.abs(grads)))
        hist[-1]["max_grad"] = float(abs(grads[m]))
        if abs(grads[m]) < grad_tol:
            break
        ops.append(m)
        gens = [pool_gens[i] for i in ops]
        theta, E, psi = _optimize_product(H, hf, gens, np.concatenate([theta, [0.0]]), maxiter)
        hist.append(dict(n_ops=len(ops), E=E, infidelity=1 - abs(np.vdot(target, psi)) ** 2,
                         max_grad=None, op=pairs[m], cnots=48 * len(ops), theta=theta.copy(), psi=psi.copy()))
    return dict(history=hist, theta=theta, ops=[pairs[i] for i in ops], psi=psi, E=E, pool=pairs)
