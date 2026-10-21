"""Lattice refinement + rodeo for the Hubbard model (paper III).

Sector-exact implementation of the EIGEN-Q pipeline for the Hubbard model,
H = -t sum_<ij>s (c+_is c_js + h.c.) + U sum_i n_iu n_id, in the conserved
(N_up, N_dn) sector.  The ED conventions are those of paper3/Hubbard.ipynb:
spin-orbitals are Jordan-Wigner ordered in two blocks (all up orbitals, then
all down orbitals), one-species occupation masks are sorted uint64 integers,
and the sector Hamiltonian factorises as

    H = K_up (x) 1 + 1 (x) K_dn + U diag(double occupancy).

Lattice refinement (Bogner et al. 2026): every coarse site c owns a block
B(c) of fine sites and the prolongation maps the coarse orbital to the
normalized symmetric (bonding) combination of its block.  The Galerkin coarse
model H_low = P^T H_high P is again a Hubbard model (chain: (t,U)->(t/2,U/2)
minus t N; square: (t,U)->(t/2,U/4) minus 2 t N), see `galerkin_params`.

The refinement paths, the exact slice propagation, the ensemble rodeo
infidelity and the Trotterized versions mirror refine.py / rodeo.py /
trotter.py of the pairing model; the symmetry diagnostics (total spin,
reflection, C4v irreps) are new.
"""
import numpy as np
from math import comb
from itertools import combinations, product
from scipy.sparse import csr_matrix, coo_matrix, identity, kron, diags, issparse
from scipy.sparse.linalg import eigsh, expm_multiply
from scipy.optimize import minimize

_ONE = np.uint64(1)

__all__ = [
    "chain", "square", "refine_chain", "refine_square", "galerkin_params",
    "spin_basis", "hop_spin", "docc", "H_hubbard", "hubbard_terms", "ground",
    "prolong_spin", "prolong", "prolonged_state",
    "path_A", "path_B", "slice_eigs", "evolve", "refine_overlaps",
    "refine_sparse", "path_gap",
    "spectral_weights", "effective_gap", "rodeo_infidelity",
    "apply_Splus", "total_spin", "site_permutation", "reflection_chain",
    "c4v_elements", "c4v_weights", "c4v_irrep",
    "trotter_refine", "trotter_step_matrix", "trotter_evolution", "rodeo_track_U", "orbital_transform",
    "hubbard_pauli_counts", "chain_step_cnots", "conditioned_term_cnots",
    "uccsd_hubbard", "uccsd_pool_size", "uccsd_cnots", "hubbard_full_jw",
]

# ------------------------------------------------------------------ lattices
def chain(L, pbc=False):
    """Open (or periodic) chain: (number of sites, list of bonds (i, j))."""
    bonds = [(i, i + 1) for i in range(L - 1)]
    if pbc and L > 2:
        bonds.append((L - 1, 0))
    return L, bonds


def square(L, pbc=False):
    """L x L square lattice with site index i = x L + y (notebook convention)."""
    site = lambda x, y: (x % L) * L + (y % L)
    bonds = []
    for x in range(L):
        for y in range(L):
            if x + 1 < L or pbc:
                bonds.append((site(x, y), site(x + 1, y)))
            if y + 1 < L or pbc:
                bonds.append((site(x, y), site(x, y + 1)))
    return L * L, bonds


def refine_chain(L, b=2, pbc=False):
    """Coarse chain of L sites -> fine chain of bL sites; block(c) = {bc, ..., bc+b-1}.
    Returns ((Nc, bonds_c), (Nf, bonds_f), blocks, intra-block bonds)."""
    Nc, bc = chain(L, pbc)
    Nf, bf = chain(b * L, pbc)
    blocks = [list(range(b * c, b * c + b)) for c in range(L)]
    intra = [(f, f + 1) for blk in blocks for f in blk[:-1]]
    return (Nc, bc), (Nf, bf), blocks, intra


def refine_square(L, pbc=False):
    """Coarse L x L -> fine 2L x 2L; block(c=(x,y)) is the 2x2 plaquette at (2x,2y)."""
    Nc, bc = square(L, pbc)
    Nf, bf = square(2 * L, pbc)
    Lf = 2 * L
    site = lambda x, y: x * Lf + y
    blocks, intra = [], []
    for x in range(L):
        for y in range(L):
            blk = [site(2 * x, 2 * y), site(2 * x + 1, 2 * y),
                   site(2 * x, 2 * y + 1), site(2 * x + 1, 2 * y + 1)]
            blocks.append(blk)
            intra += [(blk[0], blk[1]), (blk[2], blk[3]), (blk[0], blk[2]), (blk[1], blk[3])]
    return (Nc, bc), (Nf, bf), blocks, intra


def galerkin_params(kind, t, U, N):
    """(t_low, U_low, constant) of the closed-form Galerkin coarse model
    H_low = H_Hubbard(t_low, U_low) + constant for a block refinement with the
    bonding prolongation: chain (b=2): (t/2, U/2, -tN); square (b=4): (t/2, U/4, -2tN)."""
    if kind == "chain":
        return t / 2, U / 2, -t * N
    if kind == "square":
        return t / 2, U / 4, -2 * t * N
    raise ValueError(kind)


# ------------------------------------------------------------------ ED (notebook conventions)
def spin_basis(N, n):
    """Sorted uint64 occupation masks of one spin species: n of N bits set."""
    S = np.fromiter((sum(1 << p for p in c) for c in combinations(range(N), n)),
                    dtype=np.uint64, count=comb(N, n))
    S.sort()
    return S


def hop_spin(N, bonds, S, amps=None):
    """One-species hopping matrix -sum_b t_b (c+_i c_j + h.c.) on the sorted masks S
    (t_b = 1 unless `amps` gives one amplitude per bond)."""
    dim = len(S)
    rows, cols, vals = [], [], []
    cols0 = np.arange(dim, dtype=np.int64)
    for k, (i, j) in enumerate(bonds):
        tb = 1.0 if amps is None else amps[k]
        if tb == 0.0:
            continue
        mn, mx = min(i, j), max(i, j)
        btw = ((_ONE << np.uint64(mx)) - _ONE) ^ ((_ONE << np.uint64(mn + 1)) - _ONE)
        for (c, d) in ((i, j), (j, i)):                 # hop d -> c
            bc = _ONE << np.uint64(c); bd = _ONE << np.uint64(d)
            m = ((S & bd) != 0) & ((S & bc) == 0)
            if not m.any():
                continue
            Sm = S[m]
            sgn = np.bitwise_count(Sm & btw).astype(np.int64) & 1
            rows.append(np.searchsorted(S, Sm ^ (bc | bd)))
            cols.append(cols0[m])
            vals.append(-tb * (1.0 - 2.0 * sgn))
    if not rows:
        return csr_matrix((dim, dim))
    return coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(dim, dim)).tocsr()


def docc(Sup, Sdn):
    """Double occupancy popcount(up & dn) in kron order n = i_up*len(Sdn) + i_dn."""
    return np.bitwise_count(Sup[:, None] & Sdn[None, :]).astype(np.float64).ravel()


def H_hubbard(N, bonds, Nup, Ndn, U, t=1.0, amps=None):
    """Hubbard Hamiltonian in the (Nup, Ndn) sector (csr) and the mask arrays (Sup, Sdn)."""
    Sup, Sdn = spin_basis(N, Nup), spin_basis(N, Ndn)
    Kup = t * hop_spin(N, bonds, Sup, amps)
    Kdn = t * hop_spin(N, bonds, Sdn, amps)
    H = kron(Kup, identity(len(Sdn), format="csr"), format="csr") \
        + kron(identity(len(Sup), format="csr"), Kdn, format="csr")
    if U != 0.0:
        H = H + diags(U * docc(Sup, Sdn))
    return H.tocsr(), (Sup, Sdn)


def _kron_both(Kup, Kdn):
    nu, nd = Kup.shape[0], Kdn.shape[0]
    return (kron(Kup, identity(nd, format="csr"), format="csr")
            + kron(identity(nu, format="csr"), Kdn, format="csr")).tocsr()


def hubbard_terms(N, bonds, Nup, Ndn, U, t=1.0, layers=None, amps=None):
    """Term decomposition of H for product formulas: one hopping term per bond
    *layer* (a set of disjoint bonds whose exponential is a product of commuting
    two-qubit blocks on a device) for both spins, then the diagonal interaction
    term.  Default layers for a chain: even bonds, odd bonds.  Returns the list
    of csr matrices [K_layer_1, K_layer_2, ..., U diag]; their sum is H."""
    Sup, Sdn = spin_basis(N, Nup), spin_basis(N, Ndn)
    if layers is None:
        layers = [[k for k in range(len(bonds)) if k % 2 == 0],
                  [k for k in range(len(bonds)) if k % 2 == 1]]
    terms = []
    for lay in layers:
        bl = [bonds[k] for k in lay]
        al = None if amps is None else [amps[k] for k in lay]
        terms.append(_kron_both(t * hop_spin(N, bl, Sup, al), t * hop_spin(N, bl, Sdn, al)))
    terms.append(diags(U * docc(Sup, Sdn)).tocsr())
    return terms


def ground(H, k=2, dense_cutoff=2500):
    """k lowest eigenpairs (dense eigh below `dense_cutoff`, ARPACK Lanczos above)."""
    dim = H.shape[0]
    if dim <= dense_cutoff:
        w, V = np.linalg.eigh(H.toarray() if issparse(H) else H)
        return w[:k], V[:, :k]
    w, V = eigsh(H, k=k, which="SA", tol=1e-12, maxiter=40000)
    o = np.argsort(w)
    return w[o], V[:, o]


# ------------------------------------------------------------------ prolongation
def _perm_parity(seq):
    seq = list(seq); p = 0
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                p ^= 1
    return p


def prolong_spin(Nc, Nf, blocks, n, sparse=False):
    """One-species prolongation (fine_dim x coarse_dim): coarse mode c ->
    |B(c)|^{-1/2} sum_{f in B(c)} c+_f.  For a coarse occupation pattern the fine
    amplitudes sit on every choice of one fine site per occupied block, with the
    sign of the permutation that sorts the chosen fine indices."""
    Sc, Sf = spin_basis(Nc, n), spin_basis(Nf, n)
    rows, cols, vals = [], [], []
    for ic, mc in enumerate(Sc):
        occ = [c for c in range(Nc) if (int(mc) >> c) & 1]
        amp = float(np.prod([1 / np.sqrt(len(blocks[c])) for c in occ])) if occ else 1.0
        for choice in product(*[blocks[c] for c in occ]):
            mf = sum(1 << f for f in choice)
            rows.append(int(np.searchsorted(Sf, np.uint64(mf)))); cols.append(ic)
            vals.append(-amp if _perm_parity(choice) else amp)
    P = coo_matrix((vals, (rows, cols)), shape=(len(Sf), len(Sc))).tocsr()
    return P if sparse else P.toarray()


def prolong(Nc, Nf, blocks, Nup, Ndn, sparse=False):
    """Sector prolongation P = P_up (x) P_dn (an isometry, P^T P = 1)."""
    Pu, Pd = prolong_spin(Nc, Nf, blocks, Nup, True), prolong_spin(Nc, Nf, blocks, Ndn, True)
    P = kron(Pu, Pd, format="csr")
    return P if sparse else P.toarray()


def prolonged_state(Hf, P, which=0):
    """Galerkin coarse model H_low = P^T H_high P, its eigenpairs, and the
    prolonged coarse eigenstate `which`.  Returns dict(Hl, wl, Vl, Phi0, El)."""
    Pc = csr_matrix(P) if not issparse(P) else P
    Hl = Pc.T @ (Hf @ Pc)
    Hl = Hl.toarray() if issparse(Hl) else np.asarray(Hl)
    wl, Vl = np.linalg.eigh(Hl)
    Phi0 = np.asarray(Pc @ Vl[:, which]).ravel()
    return dict(Hl=Hl, wl=wl, Vl=Vl, Phi0=Phi0, El=wl[which])


# ------------------------------------------------------------------ refinement paths
def path_A(Hf, Hl, P, mu):
    """Path A (projector path): H(s) = cos^2 P(H_low-mu)P^T + sin^2 (H_high-mu).
    Returns the dense pair (H_le, H_hs)."""
    Hf = Hf.toarray() if issparse(Hf) else Hf
    P = P.toarray() if issparse(P) else P
    Hle = P @ (Hl - mu * np.eye(len(Hl))) @ P.T
    return Hle, Hf - mu * np.eye(len(Hf))


def path_B(Hf, Kintra, Lam):
    """Path B (un-dimerization): H(s) = H_high + Lam cos^2 K_intra, with K_intra the
    intra-block hopping (both spins, t = 1 amplitude; the sign is that of hop_spin,
    i.e. an *extra* hopping of strength Lam).  Returns (H_le, H_hs)."""
    Hf = Hf.toarray() if issparse(Hf) else Hf
    Kintra = Kintra.toarray() if issparse(Kintra) else Kintra
    return Hf + Lam * Kintra, Hf


def slice_eigs(Hle, Hhs, M):
    """Eigendecompositions of H(s_m) on the midpoint grid s_m = (m+1/2)/M."""
    out = []
    for m in range(M):
        s = (m + 0.5) / M
        out.append(np.linalg.eigh(np.cos(np.pi * s / 2) ** 2 * Hle + np.sin(np.pi * s / 2) ** 2 * Hhs))
    return out


def evolve(Phi0, eigs, T):
    """Exact slice propagation of Phi0 for total time T (reuses `eigs` for every T)."""
    psi = Phi0.astype(complex).copy(); dt = T / len(eigs)
    for lam, V in eigs:
        psi = V @ (np.exp(-1j * lam * dt) * (V.conj().T @ psi))
    return psi


def refine_overlaps(Phi0, Hle, Hhs, Ts, target, M=96):
    """Overlap p(T) = |<target|psi(T)>|^2 along the path for every T in Ts."""
    eigs = slice_eigs(Hle, Hhs, M)
    return np.array([abs(np.vdot(target, evolve(Phi0, eigs, T))) ** 2 for T in Ts]), eigs


def refine_sparse(Phi0, Hle, Hhs, T, M=48):
    """Slice propagation with a sparse Krylov exponential (large sectors)."""
    Hle, Hhs = csr_matrix(Hle), csr_matrix(Hhs)
    psi = Phi0.astype(complex); dt = T / M
    for m in range(M):
        s = (m + 0.5) / M
        H = np.cos(np.pi * s / 2) ** 2 * Hle + np.sin(np.pi * s / 2) ** 2 * Hhs
        psi = expm_multiply(-1j * dt * H.tocsc(), psi)
    return psi


def path_gap(Hle, Hhs, ns=41, k=2):
    """Lowest k eigenvalues of H(s) on a uniform s grid -> (s, eigenvalue array)."""
    s_grid = np.linspace(0, 1, ns); W = []
    for s in s_grid:
        H = np.cos(np.pi * s / 2) ** 2 * Hle + np.sin(np.pi * s / 2) ** 2 * Hhs
        W.append(np.linalg.eigvalsh(H)[:k])
    return s_grid, np.array(W)


# ------------------------------------------------------------------ spectral diagnostics
def spectral_weights(psi, w, V):
    """|<E_k|psi>|^2 over the eigenpairs (w, V)."""
    return np.abs(V.conj().T @ psi) ** 2


def effective_gap(psi, w, V, tol=1e-10):
    """Delta_eff: lowest excitation energy (above w[0]) carrying weight > tol in psi."""
    c2 = spectral_weights(psi, w, V)
    idx = np.where(c2[1:] > tol)[0]
    if idx.size == 0:
        return np.inf, c2
    return float(w[1 + idx[0]] - w[0]), c2


def rodeo_infidelity(psi, w, V, Ms, sigma, E_obj=None):
    """Exact ensemble infidelity 1-F_M (paper-I Eq. fid) of the Gaussian-time rodeo
    filter with M cycles, width sigma, objective energy E_obj (default w[0]),
    for an input expanded in the eigenbasis (w, V).  Also returns the ensemble
    acceptance after M cycles (second array)."""
    c2 = spectral_weights(psi, w, V)
    E0 = w[0] if E_obj is None else E_obj
    q = 0.5 * (1 + np.exp(-(w - E0) ** 2 * sigma ** 2 / 2))
    inf, acc = [], []
    for M in Ms:
        tot = c2 * q ** M
        inf.append(float(np.sum(tot[1:]) / np.sum(tot)))
        acc.append(float(np.sum(tot)))
    return np.array(inf), np.array(acc)


# ------------------------------------------------------------------ symmetry diagnostics
def apply_Splus(psi, N, Nup, Ndn):
    """S^+ = sum_i c+_{i up} c_{i dn} applied to a sector vector: returns the
    vector in the (Nup+1, Ndn-1) sector (block Jordan-Wigner signs included)."""
    Sup, Sdn = spin_basis(N, Nup), spin_basis(N, Ndn)
    Sup2, Sdn2 = spin_basis(N, Nup + 1), spin_basis(N, Ndn - 1)
    Psi = np.asarray(psi).reshape(len(Sup), len(Sdn))
    out = np.zeros((len(Sup2), len(Sdn2)), dtype=Psi.dtype)
    nup_parity = (Nup & 1)                       # c_{i dn} passes all Nup up electrons
    for i in range(N):
        bi = _ONE << np.uint64(i); below = bi - _ONE
        mu = (Sup & bi) == 0; md = (Sdn & bi) != 0
        if not mu.any() or not md.any():
            continue
        su = 1.0 - 2.0 * (np.bitwise_count(Sup[mu] & below).astype(np.int64) & 1)
        sd = 1.0 - 2.0 * ((np.bitwise_count(Sdn[md] & below).astype(np.int64) + nup_parity) & 1)
        ru = np.searchsorted(Sup2, Sup[mu] | bi); rd = np.searchsorted(Sdn2, Sdn[md] ^ bi)
        out[np.ix_(ru, rd)] += (su[:, None] * sd[None, :]) * Psi[np.ix_(np.where(mu)[0], np.where(md)[0])]
    return out.ravel()


def total_spin(psi, N, Nup, Ndn):
    """<S^2> = Sz(Sz+1) + |S^+ psi|^2 for a normalized sector vector; returns (S(S+1), S)."""
    psi = np.asarray(psi); psi = psi / np.linalg.norm(psi)
    Sz = 0.5 * (Nup - Ndn)
    s2 = Sz * (Sz + 1) + (float(np.vdot(apply_Splus(psi, N, Nup, Ndn), apply_Splus(psi, N, Nup, Ndn)).real)
                          if Ndn > 0 else 0.0)
    S = 0.5 * (-1 + np.sqrt(1 + 4 * s2))
    return s2, S


def site_permutation(perm, N, Nup, Ndn, sparse=True):
    """Sector matrix of the fermionic site permutation c+_i -> c+_{perm[i]} (same map
    for both spins), i.e. a lattice symmetry operation; a signed permutation
    matrix P_up (x) P_dn."""
    def one(n):
        S = spin_basis(N, n); rows, cols, vals = [], [], []
        for a, m in enumerate(S):
            occ = [p for p in range(N) if (int(m) >> p) & 1]
            img = [perm[p] for p in occ]
            rows.append(int(np.searchsorted(S, np.uint64(sum(1 << q for q in img)))))
            cols.append(a); vals.append(-1.0 if _perm_parity(img) else 1.0)
        return coo_matrix((vals, (rows, cols)), shape=(len(S), len(S))).tocsr()
    G = kron(one(Nup), one(Ndn), format="csr")
    return G if sparse else G.toarray()


def reflection_chain(N):
    """Site map of the reflection i -> N-1-i of the open chain."""
    return [N - 1 - i for i in range(N)]


def c4v_elements(L):
    """The eight site maps of C4v acting on the L x L square lattice (i = xL + y)
    about its center, keyed by the standard class labels; and the character
    table of C4v for the five irreps (A1, A2, B1, B2, E)."""
    site = lambda x, y: x * L + y
    def mk(f):
        return [site(*f(x, y)) for x in range(L) for y in range(L)]
    c = L - 1
    ops = {
        "E":    mk(lambda x, y: (x, y)),
        "C4":   mk(lambda x, y: (c - y, x)),
        "C4^3": mk(lambda x, y: (y, c - x)),
        "C2":   mk(lambda x, y: (c - x, c - y)),
        "sx":   mk(lambda x, y: (c - x, y)),
        "sy":   mk(lambda x, y: (x, c - y)),
        "sd":   mk(lambda x, y: (y, x)),
        "sd'":  mk(lambda x, y: (c - y, c - x)),
    }
    classes = {"E": ["E"], "2C4": ["C4", "C4^3"], "C2": ["C2"], "2sv": ["sx", "sy"], "2sd": ["sd", "sd'"]}
    chars = {  # rows: irrep, columns: E, 2C4, C2, 2sv, 2sd
        "A1": [1, 1, 1, 1, 1], "A2": [1, 1, 1, -1, -1],
        "B1": [1, -1, 1, 1, -1], "B2": [1, -1, 1, -1, 1], "E": [2, 0, -2, 0, 0]}
    return ops, classes, chars


def c4v_weights(psi, L, Nup, Ndn):
    """Weight of a normalized sector vector in each C4v irrep,
    w_G = (d_G/|G|) sum_g chi_G(g)^* <psi|g|psi>; sums to one."""
    ops, classes, chars = c4v_elements(L)
    N = L * L
    psi = np.asarray(psi); psi = psi / np.linalg.norm(psi)
    ev = {}
    for name, perm in ops.items():
        ev[name] = float(np.vdot(psi, site_permutation(perm, N, Nup, Ndn) @ psi).real)
    cls_order = ["E", "2C4", "C2", "2sv", "2sd"]
    out = {}
    for irrep, ch in chars.items():
        d = ch[0]; tot = 0.0
        for ci, cl in enumerate(cls_order):
            tot += ch[ci] * sum(ev[g] for g in classes[cl])
        out[irrep] = d * tot / 8
    return out, ev


def c4v_irrep(psi, L, Nup, Ndn, tol=1e-6):
    """Label of the dominant irrep of psi (and its weight)."""
    w, _ = c4v_weights(psi, L, Nup, Ndn)
    k = max(w, key=w.get)
    return k, w[k]


# ------------------------------------------------------------------ Trotterized pipeline
def _eigterms(terms):
    return [np.linalg.eigh(T.toarray() if issparse(T) else T) for T in terms]


def _factor_apply(eig, c, psi):
    lam, V = eig
    return V @ (np.exp(-1j * lam * c) * (V.conj().T @ psi))


def trotter_refine(Phi0, groups, coefs, T, n_s, order=2, eigs=None):
    """Product-formula refinement: H(s) = sum_i coefs[i](s) groups[i], with n_s
    steps of the first- or second-order formula (the circuit-level counterpart of
    the exact slice propagation).  `coefs` is a list of callables s -> float."""
    if eigs is None:
        eigs = _eigterms(groups)
    psi = Phi0.astype(complex).copy(); dt = T / n_s
    for m in range(n_s):
        s = (m + 0.5) / n_s
        c = [f(s) for f in coefs]
        if order == 1:
            for e, ci in zip(eigs, c):
                psi = _factor_apply(e, ci * dt, psi)
        else:
            for e, ci in zip(eigs, c):
                psi = _factor_apply(e, ci * dt / 2, psi)
            for e, ci in zip(reversed(eigs), reversed(c)):
                psi = _factor_apply(e, ci * dt / 2, psi)
    return psi / np.linalg.norm(psi)


def trotter_step_matrix(terms, dt, order=2):
    """Dense matrix of one product-formula step exp(-i H dt), H = sum(terms)."""
    eigs = _eigterms(terms)
    dim = terms[0].shape[0]
    U = np.eye(dim, dtype=complex)
    seq = eigs if order == 1 else eigs + eigs[::-1]
    c = dt if order == 1 else dt / 2
    for lam, V in seq:
        U = (V * np.exp(-1j * lam * c)) @ (V.conj().T @ U)
    return U


def trotter_evolution(terms, dt, order=2, target=None):
    """Fixed-step product formula for exp(-i H t), diagonalized once: the step
    unitary U_step is Schur-decomposed (it is normal, so the Schur form is
    diagonal) and U_step^n = W e^{i n phi} W^+ is applied with two matrix-vector
    products whatever n.  Returns dict(Ufun, nsteps, tq, phases, W, E_eff, floor):
    Ufun(n) applies U_step^n to a vector; nsteps(t) = max(1, round(|t|/dt)) and
    tq(t) = nsteps(t) dt is the quantized evolution time; if `target` (the exact
    ground state) is given, E_eff = -phi_k/dt is the Trotter-effective ground
    energy (the eigenphase of the U_step eigenvector of largest overlap with the
    target, which is what an energy scan on the device locates) and `floor` =
    1 - |<target|w_k>|^2 is the infidelity floor of the Trotterized filter."""
    from scipy.linalg import schur
    Ust = trotter_step_matrix(terms, dt, order)
    T, W = schur(Ust, output="complex")
    phases = np.angle(np.diag(T))
    out = dict(phases=phases, W=W, dt=dt)
    out["nsteps"] = lambda t: max(1, int(round(abs(t) / dt)))
    out["tq"] = lambda t: out["nsteps"](t) * dt
    def Ufun(n):
        def apply(psi):
            return W @ (np.exp(1j * n * phases) * (W.conj().T @ psi))
        return apply
    out["Ufun"] = Ufun
    if target is not None:
        ov = np.abs(W.conj().T @ target) ** 2
        k = int(ov.argmax())
        out["E_eff"] = float(-phases[k] / dt); out["floor"] = float(1 - ov[k]); out["k"] = k
    return out


def rodeo_track_U(v0, Ufun, ts, E, target, dt=None):
    """Rodeo fidelity / cumulative acceptance per cycle with a Trotterized
    controlled evolution: Ufun(n) returns a callable applying n fixed steps of
    length dt (from `trotter_evolution`); the cycle time is quantized to n dt and
    the phase uses that time.  With dt=None, Ufun(t) is a callable for the exact
    time t (as in paper I)."""
    v = v0.astype(complex).copy(); P = 1.0; F, A = [], []
    for t in ts:
        if dt is None:
            w = 0.5 * (v + np.exp(1j * E * t) * Ufun(t)(v))
        else:
            n = max(1, int(round(abs(t) / dt))); tq = n * dt
            w = 0.5 * (v + np.exp(1j * E * tq) * Ufun(n)(v))
        p = float(np.vdot(w, w).real); v = w / np.sqrt(p); P *= p
        F.append(abs(np.vdot(target, v)) ** 2); A.append(P)
    return np.array(F), np.array(A)


# ------------------------------------------------------------------ gate counting
def hubbard_pauli_counts(N, bonds):
    """Pauli content of the Hubbard Hamiltonian on 2N qubits in the block
    Jordan-Wigner ordering (up orbitals 0..N-1, down orbitals N..2N-1):
    hopping strings (X X + Y Y, with Z strings between non-adjacent sites), ZZ
    interaction strings and single-Z terms.  Returns dict with the numbers of
    strings per class and the maximal weight."""
    nb = len(bonds)
    wmax = max(abs(i - j) + 1 for i, j in bonds) if bonds else 0
    return dict(hopping=4 * nb, interaction=N, single_Z=2 * N,
                total=4 * nb + 3 * N, max_weight=max(wmax, 2))


def chain_step_cnots(N, bonds, order=1, controlled=False, absorb_scalar=True):
    """CNOT count of one product-formula step of the Hubbard Hamiltonian with the
    Pauli-rotation primitive of paper I: a weight-w rotation costs 2(w-1) CNOTs
    (XX+YY on adjacent qubits: 2 rotations x 2 CNOTs), a ZZ rotation 2 CNOTs,
    phases are free; controlling a step adds 2 CNOTs per rotation.  Second order
    doubles the count.  At fixed particle number the sum of the 2N single-Z
    terms of the interaction is a scalar, so in a CONTROLLED step those 2N
    phases are a single phase on the control (absorbed into the objective-energy
    phase gate) and are not charged when absorb_scalar=True (the default since
    the revision of 2026-10-20; absorb_scalar=False gives the earlier counts,
    352 / 736 instead of 288 / 608 at 8 / 16 sites)."""
    hop_cnots = sum(2 * 2 * (abs(i - j)) * 2 for i, j in bonds)   # 2 strings/bond/spin, weight |i-j|+1, 2 spins
    nrot = 4 * len(bonds) + N + 2 * N
    nrot_ctrl = nrot - (2 * N if absorb_scalar else 0)
    cn = hop_cnots + 2 * N
    if controlled:
        cn += 2 * nrot_ctrl
    f = 2 if order == 2 else 1
    return dict(cnots=cn * f, rotations=nrot * f, rotations_controlled=nrot_ctrl * f)


def conditioned_term_cnots(Nc, bonds_c, n_controls):
    """CNOTs per application of the embedded coarse term P(H_low-mu)P^T of Path A:
    a coarse Hubbard step (bonding-mode register) in which every Pauli rotation
    (single-qubit phases included, since conditioned they are no longer free)
    is conditioned on the n_controls antibonding qubits being empty, plus the
    conditioned scalar phase of the constant part of H_low - mu (a physical
    relative phase after projection).  Conventions returned:
      'gray'    Gray-code phase polynomial over the controls per rotation
                (paper-I convention): 2^(c+1)-2 extra CNOTs per rotation,
                2^c-2 for the scalar;
      'toffoli' a clean-ancilla ladder of relative-phase Toffolis PER rotation:
                6(c-1)+2 extra CNOTs per rotation and for the scalar;
      'flag'    the conjunction of the c controls is computed ONCE into a clean
                flag qubit (6(c-1) CNOTs for the ladder and its uncomputation),
                every rotation is singly controlled by the flag (2 CNOTs) and
                the scalar is one controlled phase (2): base + 2 n_rot + 6(c-1)
                + 2, at the price of c-1 clean work qubits (122+2 = 124 CNOTs
                for 4->8, 266+2 = 268 for 8->16; revision of 2026-10-20);
      'frame'   the basis change between the site register and the
                bonding/antibonding frame around one application (the
                prolongation Givens network and its inverse, 2 x 4 N_c CNOTs),
                to be added to 'flag' unless merged with a neighboring block.
    The coarse evolution does not change the control qubits, so the flag is
    valid for superpositions and can be reused across the whole coarse block."""
    base = chain_step_cnots(Nc, bonds_c, 1)
    nrot = base["rotations"]
    c = n_controls
    return dict(rotations=nrot, base=base["cnots"],
                gray=base["cnots"] + nrot * (2 ** (c + 1) - 2) + (2 ** c - 2),
                toffoli=base["cnots"] + (nrot + 1) * (6 * (c - 1) + 2),
                flag=base["cnots"] + 2 * nrot + 6 * (c - 1) + 2,
                frame=2 * 4 * Nc)


# ------------------------------------------------------------------ sector-level UCCSD
def _one_body_spin(N, S, i, j):
    """c+_i c_j for one species on sorted masks S (i != j: hopping with JW sign;
    i == j: number operator)."""
    dim = len(S)
    if i == j:
        return diags(((S >> np.uint64(i)) & _ONE).astype(np.float64)).tocsr()
    bi, bj = _ONE << np.uint64(i), _ONE << np.uint64(j)
    mn, mx = min(i, j), max(i, j)
    btw = ((_ONE << np.uint64(mx)) - _ONE) ^ ((_ONE << np.uint64(mn + 1)) - _ONE)
    m = ((S & bj) != 0) & ((S & bi) == 0)
    if not m.any():
        return csr_matrix((dim, dim))
    Sm = S[m]
    sgn = 1.0 - 2.0 * (np.bitwise_count(Sm & btw).astype(np.int64) & 1)
    return coo_matrix((sgn, (np.searchsorted(S, Sm ^ (bi | bj)), np.where(m)[0])),
                      shape=(dim, dim)).tocsr()


def orbital_transform(N, n, phi):
    """One-species change of basis from site occupation masks to orbital
    occupation masks: Phi[I, J] = det(phi[sites(I), orbitals(J)]), so that
    |J>_orb = sum_I Phi[I, J] |I>_site; Phi is orthogonal."""
    S = spin_basis(N, n)
    occ = [[p for p in range(N) if (int(m) >> p) & 1] for m in S]
    Phi = np.zeros((len(S), len(S)))
    for J, oJ in enumerate(occ):
        for I, oI in enumerate(occ):
            Phi[I, J] = np.linalg.det(phi[np.ix_(oI, oJ)]) if n else 1.0
    return Phi


def uccsd_pool_size(N, Nup, Ndn):
    """Number of singles and doubles of a spin-resolved UCCSD pool (HF reference of
    Nup/Ndn lowest orbitals): singles n_occ n_virt per spin; doubles: same-spin
    C(n_occ,2) C(n_virt,2) per spin + opposite-spin products."""
    s = Nup * (N - Nup) + Ndn * (N - Ndn)
    d = comb(Nup, 2) * comb(N - Nup, 2) + comb(Ndn, 2) * comb(N - Ndn, 2) \
        + Nup * (N - Nup) * Ndn * (N - Ndn)
    return s, d


def uccsd_hubbard(N, bonds, Nup, Ndn, U, t=1.0, maxiter=500, H=None, verbose=False, keep=None, theta0=None):
    """Sector-level UCCSD-VQE for the Hubbard model with the free-fermion (U=0)
    orbitals as reference: |psi> = prod_mu exp(theta_mu (tau_mu - tau_mu^T)) |HF>,
    tau = spin-resolved singles and doubles in the orbital basis, applied as a
    product in a fixed order (singles first, then doubles).  The sector is
    transformed to the orbital occupation basis (one orthogonal determinant map
    per spin), where every generator connects pairs of determinants and its
    exponential is the exact Givens-type rotation
    exp(theta K) = 1 + sin(theta) K + (cos(theta) - 1) Pi, Pi = -K^2 — the same
    structure as the compiled excitation blocks of paper I.  Optimised with
    L-BFGS-B and the analytic adjoint gradient.  Returns dict(E, psi [site sector
    basis], theta, nit, E_hf, n_params, singles, doubles, eps, phi)."""
    if H is None:
        H, _ = H_hubbard(N, bonds, Nup, Ndn, U, t)
    H = H.toarray() if issparse(H) else np.asarray(H)
    T1 = np.zeros((N, N))
    for i, j in bonds:
        T1[i, j] -= t; T1[j, i] -= t
    eps, phi = np.linalg.eigh(T1)                       # orbitals phi[:, a], ascending energy
    Phi_u, Phi_d = orbital_transform(N, Nup, phi), orbital_transform(N, Ndn, phi)
    Phi = np.kron(Phi_u, Phi_d)
    Horb = Phi.T @ H @ Phi
    Su, Sd = spin_basis(N, Nup), spin_basis(N, Ndn)     # orbital occupation masks (same combinatorics)
    Iu, Id = identity(len(Su), format="csr"), identity(len(Sd), format="csr")
    Eu = lambda a, b: _one_body_spin(N, Su, a, b)
    Ed = lambda a, b: _one_body_spin(N, Sd, a, b)
    occ_u, vir_u = list(range(Nup)), list(range(Nup, N))
    occ_d, vir_d = list(range(Ndn)), list(range(Ndn, N))
    gens, singles, doubles = [], [], []
    for (i, a) in product(occ_u, vir_u):
        gens.append(kron(Eu(a, i), Id, format="csr")); singles.append(("u", i, a))
    for (i, a) in product(occ_d, vir_d):
        gens.append(kron(Iu, Ed(a, i), format="csr")); singles.append(("d", i, a))
    for (i, j) in combinations(occ_u, 2):
        for (a, b) in combinations(vir_u, 2):
            gens.append(kron(Eu(a, i) @ Eu(b, j), Id, format="csr")); doubles.append(("uu", i, j, a, b))
    for (i, j) in combinations(occ_d, 2):
        for (a, b) in combinations(vir_d, 2):
            gens.append(kron(Iu, Ed(a, i) @ Ed(b, j), format="csr")); doubles.append(("dd", i, j, a, b))
    for (i, a) in product(occ_u, vir_u):
        for (j, b) in product(occ_d, vir_d):
            gens.append(kron(Eu(a, i), Ed(b, j), format="csr")); doubles.append(("ud", i, j, a, b))
    if keep is not None:                                 # pruned pool: keep the listed operator indices (order kept)
        keep = sorted(keep); labels = singles + doubles
        gens = [gens[i] for i in keep]; labels = [labels[i] for i in keep]
        singles = [l for l in labels if len(l) == 3]; doubles = [l for l in labels if len(l) == 5]
    Ks = [(g - g.T).tocsr() for g in gens]
    Pis = []
    for K in Ks:
        Pi = (-(K @ K)).tocsr()
        d = Pi.diagonal()
        assert abs(Pi - diags(d)).max() < 1e-12 and np.all((abs(d) < 1e-12) | (abs(d - 1) < 1e-12)), "generator is not a pair rotation"
        Pis.append(d)
    def fac(m, th, psi):
        return psi + np.sin(th) * (Ks[m] @ psi) + (np.cos(th) - 1.0) * (Pis[m] * psi)
    # HF reference: the lowest Nup / Ndn orbitals occupied -> one orbital basis state
    ref = np.zeros(len(Su) * len(Sd), dtype=complex)
    ju = int(np.searchsorted(Su, np.uint64(sum(1 << a for a in range(Nup)))))
    jd = int(np.searchsorted(Sd, np.uint64(sum(1 << a for a in range(Ndn)))))
    ref[ju * len(Sd) + jd] = 1.0
    P = len(Ks)
    def forward(th):
        st = [ref]
        for m in range(P):
            st.append(fac(m, th[m], st[-1]))
        return st
    def energy(th):
        psi = forward(th)[-1]
        return float(np.real(np.vdot(psi, Horb @ psi)))
    def grad(th):
        st = forward(th); psi = st[-1]
        lam_ = Horb @ psi
        g = np.zeros(P)
        for m in range(P - 1, -1, -1):
            g[m] = 2.0 * float(np.real(np.vdot(lam_, Ks[m] @ st[m + 1])))
            lam_ = fac(m, -th[m], lam_)                  # U_m^+ = exp(-theta K)
        return g
    E_hf = energy(np.zeros(P))
    x0 = np.zeros(P) if theta0 is None else np.asarray(theta0, float)
    res = minimize(energy, x0, jac=grad, method="L-BFGS-B",
                   options={"ftol": 1e-14, "gtol": 1e-10, "maxiter": maxiter})
    psi = Phi @ forward(res.x)[-1]
    return dict(E=float(res.fun), psi=psi / np.linalg.norm(psi), theta=res.x, nit=int(res.nit),
                E_hf=float(E_hf), n_params=P, singles=singles, doubles=doubles, eps=eps, phi=phi)


def uccsd_cnots(N, singles, doubles):
    """CNOT estimate of the UCCSD circuit on 2N qubits (block JW ordering): a single
    excitation i->a of one spin compiles into 2 Pauli strings of weight a-i+1
    (2(w-1) CNOTs each); a double with sorted JW indices p<q<r<s into 8 strings
    of weight (q-p+1)+(s-r+1) --- the Z strings cancel in the gap between q and
    r, so the support is the two intervals, NOT the full span s-p+1 (which the
    counter used before the revision of 2026-10-20 and which over-charged
    opposite-spin doubles: 1720 -> 1464 CNOTs for the four-site template,
    30 864 -> 23 312 for eight sites, 335 360 -> 235 904 for sixteen); plus the
    Givens network that rotates the site register into the orbital basis,
    N(N-1)/2 rotations per spin of 2 CNOTs each (applied once, after the
    excitation blocks that act in the orbital basis)."""
    cn = 0
    for (_, i, a) in singles:
        cn += 2 * 2 * (abs(a - i))
    for (kind, i, j, a, b) in doubles:
        if kind == "ud":
            idx = sorted([i, a, j + N, b + N])
        else:
            idx = sorted([i, j, a, b])
        w = (idx[1] - idx[0] + 1) + (idx[3] - idx[2] + 1)
        cn += 8 * 2 * (w - 1)
    givens = 2 * (N * (N - 1) // 2) * 2
    return dict(excitations=cn, givens=givens, total=cn + givens)


def hubbard_full_jw(N, bonds, U, t=1.0):
    """Dense Hubbard Hamiltonian on all 2N qubits in the block Jordan-Wigner
    ordering (up orbital i -> qubit i, down orbital i -> qubit N+i), built from
    JW creation/annihilation operators; for Pauli-structure checks on small N."""
    from ._exc import jw_create, jw_ann
    nq = 2 * N; dim = 2 ** nq
    H = np.zeros((dim, dim), complex)
    for i, j in bonds:
        for s in (0, N):
            H -= t * (jw_create(i + s, nq) @ jw_ann(j + s, nq) + jw_create(j + s, nq) @ jw_ann(i + s, nq))
    for i in range(N):
        nu = jw_create(i, nq) @ jw_ann(i, nq); nd = jw_create(i + N, nq) @ jw_ann(i + N, nq)
        H += U * nu @ nd
    return H


# ------------------------------------------------------------------ large sectors: Kronecker-structured product formulas
class LayerGroup:
    """A hopping layer t*(K_up (x) 1 + 1 (x) K_dn) whose exponential factorizes,
    exp(-i c K) = exp(-i c K_up) (x) exp(-i c K_dn): applied to psi reshaped as a
    (dim_up x dim_dn) matrix at the cost of two small dense products.  Kup, Kdn
    are the one-species hopping matrices on the sorted masks (dense or sparse)."""
    def __init__(self, Kup, Kdn):
        Kup = Kup.toarray() if issparse(Kup) else np.asarray(Kup)
        Kdn = Kdn.toarray() if issparse(Kdn) else np.asarray(Kdn)
        self.wu, self.Vu = np.linalg.eigh(Kup); self.wd, self.Vd = np.linalg.eigh(Kdn)
        self.nu, self.nd = len(self.wu), len(self.wd)

    def apply(self, c, psi):
        X = psi.reshape(self.nu, self.nd)
        A = (self.Vu * np.exp(-1j * c * self.wu)) @ self.Vu.conj().T
        B = (self.Vd * np.exp(-1j * c * self.wd)) @ self.Vd.conj().T
        return (A @ X @ B.T).ravel()

    def matvec(self, psi):
        X = psi.reshape(self.nu, self.nd)
        Ku = (self.Vu * self.wu) @ self.Vu.conj().T; Kd = (self.Vd * self.wd) @ self.Vd.conj().T
        return (Ku @ X + X @ Kd.T).ravel()


class DiagGroup:
    """A diagonal term (the interaction U n_up n_dn): exp(-i c d) elementwise."""
    def __init__(self, d):
        self.d = np.asarray(d, float)

    def apply(self, c, psi):
        return np.exp(-1j * c * self.d) * psi

    def matvec(self, psi):
        return self.d * psi


class ProjGroup:
    """The embedded coarse term P X P^T (P a sparse isometry, X a dense coarse
    matrix): exp(-i c P X P^T) = 1 + P (exp(-i c X) - 1) P^T, exact and cheap."""
    def __init__(self, P, X):
        self.P = csr_matrix(P); self.PT = self.P.T.tocsr()
        self.w, self.V = np.linalg.eigh(np.asarray(X))

    def apply(self, c, psi):
        y = self.PT @ psi
        return psi + self.P @ (self.V @ ((np.exp(-1j * c * self.w) - 1.0) * (self.V.conj().T @ y)))

    def matvec(self, psi):
        return self.P @ (self.V @ (self.w * (self.V.conj().T @ (self.PT @ psi))))


def hubbard_groups(N, bonds, Nup, Ndn, U, t=1.0, layers=None):
    """[LayerGroup(even bonds), LayerGroup(odd bonds), DiagGroup(U docc)] — the
    product-formula groups of `hubbard_terms` in Kronecker form (any sector size
    whose one-species dimension is dense-diagonalizable, e.g. 120 for 16 sites
    at (2,2))."""
    Sup, Sdn = spin_basis(N, Nup), spin_basis(N, Ndn)
    if layers is None:
        layers = [[k for k in range(len(bonds)) if k % 2 == 0], [k for k in range(len(bonds)) if k % 2 == 1]]
    groups = []
    for lay in layers:
        bl = [bonds[k] for k in lay]
        groups.append(LayerGroup(t * hop_spin(N, bl, Sup), t * hop_spin(N, bl, Sdn)))
    groups.append(DiagGroup(U * docc(Sup, Sdn)))
    return groups


def trotter_path_groups(Phi0, groups, coefs, T, n_s, order=2):
    """Product-formula evolution along H(s) = sum_i coefs[i](s) groups[i] with
    n_s first- or second-order steps (midpoint rule in s), for group objects
    with an `apply(c, psi)` method.  Returns the normalized final state."""
    psi = np.asarray(Phi0, complex).copy(); dt = T / n_s
    for m in range(n_s):
        s = (m + 0.5) / n_s
        c = [f(s) for f in coefs]
        if order == 1:
            for g_, ci in zip(groups, c):
                psi = g_.apply(ci * dt, psi)
        else:
            for g_, ci in zip(groups, c):
                psi = g_.apply(ci * dt / 2, psi)
            for g_, ci in zip(reversed(groups), reversed(c)):
                psi = g_.apply(ci * dt / 2, psi)
    return psi / np.linalg.norm(psi)


def step_apply_groups(groups, dt, psi, order=2):
    """One product-formula step exp(-i sum(groups) dt) applied to psi."""
    if order == 1:
        for g_ in groups:
            psi = g_.apply(dt, psi)
        return psi
    for g_ in groups:
        psi = g_.apply(dt / 2, psi)
    for g_ in reversed(groups):
        psi = g_.apply(dt / 2, psi)
    return psi


def effective_step_state(groups, dt, psi0, E0, sigma, iters=60, seed=11, order=2, return_residual=False):
    """Trotter-effective ground state of the fixed-step filter for sectors too
    large to diagonalize the step unitary: iterate the quantized-time
    post-selection K_t = (1 + e^{i E t} U_step^n)/2 (n = t/dt, t half-Gaussian of
    width sigma) on psi0 while updating E to the eigenphase estimate
    -arg<psi|U_step|psi>/dt.  Returns (psi_eff, E_eff): the step eigenvector
    closest to the input and its energy; 1-|<target|psi_eff>|^2 is the floor.
    With return_residual=True also returns ||U_step psi - e^{-i E dt} psi||,
    the eigenvector residual that certifies the convergence of the iteration."""
    rng = np.random.default_rng(seed)
    psi = np.asarray(psi0, complex) / np.linalg.norm(psi0)
    E = E0
    for it in range(iters):
        Upsi = step_apply_groups(groups, dt, psi, order)
        E = -np.angle(np.vdot(psi, Upsi)) / dt
        n = max(1, int(round(abs(rng.normal(0.0, sigma)) / dt)))
        v = psi.copy()
        for _ in range(n):
            v = step_apply_groups(groups, dt, v, order)
        psi = 0.5 * (psi + np.exp(1j * E * n * dt) * v); psi /= np.linalg.norm(psi)
    Upsi = step_apply_groups(groups, dt, psi, order)
    E = -np.angle(np.vdot(psi, Upsi)) / dt
    if return_residual:
        return psi, E, float(np.linalg.norm(Upsi - np.exp(-1j * E * dt) * psi))
    return psi, E


def quantized_step_distribution(sigma, dt, xmax=8.0):
    """Probability w_n of n = max(1, round(|t|/dt)) steps per cycle for the
    half-Gaussian cycle time |t|, t ~ N(0, sigma^2); returns (n, w_n)."""
    from scipy.stats import norm
    nmax = int(xmax * sigma / dt) + 2
    ns = np.arange(1, nmax + 1)
    edges = (ns - 0.5) * dt; edges[0] = 0.0
    cdf = lambda x: 2 * norm.cdf(x / sigma) - 1
    w = cdf(np.append(edges[1:], np.inf)) - cdf(edges)
    return ns, w / w.sum()


def exact_digital_ensemble(tr, psi0, target, sigma, M, E_obj=None):
    """EXACT ensemble of accepted runs of the fixed-step digital filter
    (`trotter_evolution` output `tr`, cycle times quantized to n dt): because
    every cycle is a function of the same step unitary U = W e^{i phi} W^+, the
    accepted density matrix after M cycles is rho^(M)_ij = c_i c_j^* Q_ij^M in
    the step eigenbasis, with c = W^+ psi0, f_j(n) = [1 + e^{i n (phi_j + E dt)}]/2
    and Q_ij = E_n[f_i(n) f_j(n)^*] over the rounded Gaussian distribution of n.
    Returns dict(inf, acc, mean_steps): the infidelity 1 - <target|rho|target>/Tr rho
    (coherences kept: the exact and the digital eigenbases differ) and the
    acceptance Tr rho for M = 0..M, and the mean number of steps per cycle.
    Replaces the sampled 40-schedule averages wherever the step can be
    diagonalized (revision of 2026-10-20)."""
    W, ph, dt = tr["W"], tr["phases"], tr["dt"]
    if E_obj is None:
        E_obj = tr["E_eff"]
    ns, wn = quantized_step_distribution(sigma, dt)
    f = 0.5 * (1 + np.exp(1j * np.outer(ns, ph + E_obj * dt)))
    Q = (f.T * wn) @ f.conj()
    c = W.conj().T @ np.asarray(psi0, complex); g = W.conj().T @ np.asarray(target, complex)
    rho = np.outer(c, c.conj())
    inf = [1 - float((g.conj() @ rho @ g).real / np.trace(rho).real)]; acc = [1.0]
    for _ in range(M):
        rho = rho * Q
        tr_ = float(np.trace(rho).real)
        acc.append(tr_); inf.append(1 - float((g.conj() @ rho @ g).real) / tr_)
    return dict(inf=np.array(inf), acc=np.array(acc), mean_steps=float(np.sum(ns * wn)))


def rodeo_track_groups(v0, groups, ts, E, target, dt, order=2):
    """Rodeo fidelity/acceptance tracking with a fixed-step product-formula
    evolution built from group objects; cycle times are quantized to multiples
    of dt (the convention of `rodeo_track_U(..., dt=)`)."""
    v = np.asarray(v0, complex).copy(); P = 1.0; F, A = [], []
    for t in ts:
        n = max(1, int(round(abs(t) / dt))); u = v.copy()
        for _ in range(n):
            u = step_apply_groups(groups, dt, u, order)
        w = 0.5 * (v + np.exp(1j * E * n * dt) * u)
        p = float(np.vdot(w, w).real); v = w / np.sqrt(p); P *= p
        F.append(abs(np.vdot(target, v)) ** 2); A.append(P)
    return np.array(F), np.array(A)


# ------------------------------------------------------------------ large sectors: UCCSD and ADAPT-VQE with sparse operators
def _uccsd_pool_kron(N, bonds, Nup, Ndn, t=1.0):
    """Orbital-basis machinery shared by `uccsd_hubbard_kron` and `adapt_hubbard`:
    the free-fermion orbitals, the one-species orbital transforms Phi_u, Phi_d,
    the HF reference in the orbital occupation basis, and the spin-resolved
    singles/doubles pool as sparse generators K = tau - tau^T with their
    projectors Pi = -K^2 (diagonal, so every factor is an exact Givens-type
    rotation).  Returns a dict."""
    T1 = np.zeros((N, N))
    for i, j in bonds:
        T1[i, j] -= t; T1[j, i] -= t
    eps, phi = np.linalg.eigh(T1)
    Phi_u, Phi_d = orbital_transform(N, Nup, phi), orbital_transform(N, Ndn, phi)
    Su, Sd = spin_basis(N, Nup), spin_basis(N, Ndn)
    Iu, Id = identity(len(Su), format="csr"), identity(len(Sd), format="csr")
    Eu = lambda a, b: _one_body_spin(N, Su, a, b)
    Ed = lambda a, b: _one_body_spin(N, Sd, a, b)
    occ_u, vir_u = list(range(Nup)), list(range(Nup, N))
    occ_d, vir_d = list(range(Ndn)), list(range(Ndn, N))
    gens, labels = [], []
    for (i, a) in product(occ_u, vir_u):
        gens.append(kron(Eu(a, i), Id, format="csr")); labels.append(("u", i, a))
    for (i, a) in product(occ_d, vir_d):
        gens.append(kron(Iu, Ed(a, i), format="csr")); labels.append(("d", i, a))
    for (i, j) in combinations(occ_u, 2):
        for (a, b) in combinations(vir_u, 2):
            gens.append(kron(Eu(a, i) @ Eu(b, j), Id, format="csr")); labels.append(("uu", i, j, a, b))
    for (i, j) in combinations(occ_d, 2):
        for (a, b) in combinations(vir_d, 2):
            gens.append(kron(Iu, Ed(a, i) @ Ed(b, j), format="csr")); labels.append(("dd", i, j, a, b))
    for (i, a) in product(occ_u, vir_u):
        for (j, b) in product(occ_d, vir_d):
            gens.append(kron(Eu(a, i), Ed(b, j), format="csr")); labels.append(("ud", i, j, a, b))
    Ks = [(g - g.T).tocsr() for g in gens]
    Pis = [(-(K @ K)).diagonal() for K in Ks]
    ref = np.zeros(len(Su) * len(Sd), dtype=complex)
    ju = int(np.searchsorted(Su, np.uint64(sum(1 << a for a in range(Nup)))))
    jd = int(np.searchsorted(Sd, np.uint64(sum(1 << a for a in range(Ndn)))))
    ref[ju * len(Sd) + jd] = 1.0
    return dict(eps=eps, phi=phi, Phi_u=Phi_u, Phi_d=Phi_d, Ks=Ks, Pis=Pis, labels=labels, ref=ref,
                nu=len(Su), nd=len(Sd))


def _to_site(pool, psi_orb):
    X = psi_orb.reshape(pool["nu"], pool["nd"])
    return (pool["Phi_u"] @ X @ pool["Phi_d"].T).ravel()


def _to_orb(pool, psi_site):
    X = psi_site.reshape(pool["nu"], pool["nd"])
    return (pool["Phi_u"].T @ X @ pool["Phi_d"]).ravel()


def _product_state_opt(pool, H, idx, theta0, maxiter=500):
    """Optimize <psi|H|psi> for psi = prod_m exp(theta_m K_{idx[m]}) |HF> (orbital
    basis; H sparse in the site basis) with L-BFGS-B and the adjoint gradient.
    Returns (theta, E, psi_site)."""
    Ks, Pis = pool["Ks"], pool["Pis"]
    def fac(m, th, psi):
        return psi + np.sin(th) * (Ks[m] @ psi) + (np.cos(th) - 1.0) * (Pis[m] * psi)
    def forward(th):
        st = [pool["ref"]]
        for n_, m in enumerate(idx):
            st.append(fac(m, th[n_], st[-1]))
        return st
    def energy(th):
        ps = _to_site(pool, forward(th)[-1])
        return float(np.real(np.vdot(ps, H @ ps)))
    def grad(th):
        st = forward(th); ps = _to_site(pool, st[-1])
        lam_ = _to_orb(pool, H @ ps); g = np.zeros(len(idx))
        for n_ in range(len(idx) - 1, -1, -1):
            m = idx[n_]
            g[n_] = 2.0 * float(np.real(np.vdot(lam_, Ks[m] @ st[n_ + 1])))
            lam_ = fac(m, -th[n_], lam_)
        return g
    if len(idx) == 0:
        ps = _to_site(pool, pool["ref"]); return np.zeros(0), float(np.real(np.vdot(ps, H @ ps))), ps
    res = minimize(energy, np.asarray(theta0, float), jac=grad, method="L-BFGS-B",
                   options={"ftol": 1e-14, "gtol": 1e-10, "maxiter": maxiter})
    ps = _to_site(pool, forward(res.x)[-1])
    return res.x, float(res.fun), ps / np.linalg.norm(ps)


def uccsd_hubbard_kron(N, bonds, Nup, Ndn, U, t=1.0, maxiter=500, H=None, keep=None, theta0=None):
    """`uccsd_hubbard` for sectors too large for dense matrices: the same ansatz,
    reference, pool and ordering, with H sparse and the orbital transform
    applied as a Kronecker product.  Returns dict(E, psi, theta, nit, E_hf,
    n_params, singles, doubles, labels)."""
    if H is None:
        H, _ = H_hubbard(N, bonds, Nup, Ndn, U, t)
    H = csr_matrix(H)
    pool = _uccsd_pool_kron(N, bonds, Nup, Ndn, t)
    idx = list(range(len(pool["Ks"]))) if keep is None else sorted(keep)
    x0 = np.zeros(len(idx)) if theta0 is None else np.asarray(theta0, float)
    ps_ref = _to_site(pool, pool["ref"]); E_hf = float(np.real(np.vdot(ps_ref, H @ ps_ref)))
    theta, E, psi = _product_state_opt(pool, H, idx, x0, maxiter)
    labels = [pool["labels"][i] for i in idx]
    return dict(E=E, psi=psi, theta=theta, E_hf=E_hf, n_params=len(idx), labels=labels,
                singles=[l for l in labels if len(l) == 3], doubles=[l for l in labels if len(l) == 5])


def adapt_hubbard(N, bonds, Nup, Ndn, U, t=1.0, grad_tol=1e-4, max_ops=60, H=None, target=None, maxiter=400):
    """ADAPT-VQE for the Hubbard sector with the spin-resolved singles/doubles pool
    of `uccsd_hubbard` (operators may repeat): at every iteration the operator
    with the largest |<psi|[H,K]|psi>| is appended and all amplitudes are
    re-optimized from the previous optimum.  Returns dict(history=[dict(n_ops,
    E, infidelity, max_grad, op, psi)], ops, theta, psi, E, pool_size)."""
    if H is None:
        H, _ = H_hubbard(N, bonds, Nup, Ndn, U, t)
    H = csr_matrix(H)
    pool = _uccsd_pool_kron(N, bonds, Nup, Ndn, t); Ks = pool["Ks"]
    ps = _to_site(pool, pool["ref"]); E = float(np.real(np.vdot(ps, H @ ps)))
    ops, theta = [], np.zeros(0)
    inf = lambda v: (1 - abs(np.vdot(target, v)) ** 2) if target is not None else None
    hist = [dict(n_ops=0, E=E, infidelity=inf(ps), max_grad=None, op=None, psi=ps)]
    psi_orb = pool["ref"]
    while len(ops) < max_ops:
        lam_ = _to_orb(pool, H @ ps)
        grads = np.array([2.0 * float(np.real(np.vdot(lam_, K @ psi_orb))) for K in Ks])
        m = int(np.argmax(np.abs(grads))); hist[-1]["max_grad"] = float(abs(grads[m]))
        if abs(grads[m]) < grad_tol:
            break
        ops.append(m)
        theta, E, ps = _product_state_opt(pool, H, ops, np.concatenate([theta, [0.0]]), maxiter)
        # orbital-basis state for the next gradient
        st = pool["ref"]
        for n_, mm in enumerate(ops):
            st = st + np.sin(theta[n_]) * (Ks[mm] @ st) + (np.cos(theta[n_]) - 1.0) * (pool["Pis"][mm] * st)
        psi_orb = st
        hist.append(dict(n_ops=len(ops), E=E, infidelity=inf(ps), max_grad=None, op=pool["labels"][m], psi=ps))
    return dict(history=hist, ops=[pool["labels"][i] for i in ops], theta=theta, psi=ps, E=E, pool_size=len(Ks))
