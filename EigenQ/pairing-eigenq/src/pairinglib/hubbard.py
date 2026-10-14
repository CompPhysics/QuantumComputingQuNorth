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
normalised symmetric (bonding) combination of its block.  The Galerkin coarse
model H_low = P^T H_high P is again a Hubbard model (chain: (t,U)->(t/2,U/2)
minus t N; square: (t,U)->(t/2,U/4) minus 2 t N), see `galerkin_params`.

The refinement paths, the exact slice propagation, the ensemble rodeo
infidelity and the Trotterised versions mirror refine.py / rodeo.py /
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
    """Path B (un-dimerisation): H(s) = H_high + Lam cos^2 K_intra, with K_intra the
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
    """<S^2> = Sz(Sz+1) + |S^+ psi|^2 for a normalised sector vector; returns (S(S+1), S)."""
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
    about its centre, keyed by the standard class labels; and the character
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
    """Weight of a normalised sector vector in each C4v irrep,
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


# ------------------------------------------------------------------ Trotterised pipeline
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
    """Fixed-step product formula for exp(-i H t), diagonalised once: the step
    unitary U_step is Schur-decomposed (it is normal, so the Schur form is
    diagonal) and U_step^n = W e^{i n phi} W^+ is applied with two matrix-vector
    products whatever n.  Returns dict(Ufun, nsteps, tq, phases, W, E_eff, floor):
    Ufun(n) applies U_step^n to a vector; nsteps(t) = max(1, round(|t|/dt)) and
    tq(t) = nsteps(t) dt is the quantised evolution time; if `target` (the exact
    ground state) is given, E_eff = -phi_k/dt is the Trotter-effective ground
    energy (the eigenphase of the U_step eigenvector of largest overlap with the
    target, which is what an energy scan on the device locates) and `floor` =
    1 - |<target|w_k>|^2 is the infidelity floor of the Trotterised filter."""
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
    """Rodeo fidelity / cumulative acceptance per cycle with a Trotterised
    controlled evolution: Ufun(n) returns a callable applying n fixed steps of
    length dt (from `trotter_evolution`); the cycle time is quantised to n dt and
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


def chain_step_cnots(N, bonds, order=1, controlled=False):
    """CNOT count of one product-formula step of the Hubbard Hamiltonian with the
    Pauli-rotation primitive of paper I: a weight-w rotation costs 2(w-1) CNOTs
    (XX+YY on adjacent qubits: 2 rotations x 2 CNOTs), a ZZ rotation 2 CNOTs,
    phases are free; controlling a step adds 2 CNOTs per rotation.  Second order
    doubles the count."""
    hop_cnots = sum(2 * 2 * (abs(i - j)) * 2 for i, j in bonds)   # 2 strings/bond/spin, weight |i-j|+1, 2 spins
    nrot = 4 * len(bonds) + N + 2 * N
    cn = hop_cnots + 2 * N
    if controlled:
        cn += 2 * nrot
    return dict(cnots=cn * (2 if order == 2 else 1), rotations=nrot * (2 if order == 2 else 1))


def conditioned_term_cnots(Nc, bonds_c, n_controls):
    """CNOTs per application of the embedded coarse term P(H_low-mu)P^T of Path A:
    a coarse Hubbard step (bonding-mode register) in which every Pauli rotation is
    conditioned on the n_controls antibonding qubits being empty.  Two conventions
    are returned: 'gray' (Gray-code phase polynomial over the controls, the paper-I
    convention: 2^(c+1)-2 extra CNOTs per rotation, exact for few controls) and
    'toffoli' (a clean-ancilla Toffoli ladder with relative-phase Toffolis:
    6(c-1)+2 extra CNOTs per rotation)."""
    base = chain_step_cnots(Nc, bonds_c, 1)
    nrot = base["rotations"]
    c = n_controls
    return dict(rotations=nrot, base=base["cnots"],
                gray=base["cnots"] + nrot * (2 ** (c + 1) - 2),
                toffoli=base["cnots"] + nrot * (6 * (c - 1) + 2))


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


def uccsd_hubbard(N, bonds, Nup, Ndn, U, t=1.0, maxiter=500, H=None, verbose=False):
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
    res = minimize(energy, np.zeros(P), jac=grad, method="L-BFGS-B",
                   options={"ftol": 1e-14, "gtol": 1e-10, "maxiter": maxiter})
    psi = Phi @ forward(res.x)[-1]
    return dict(E=float(res.fun), psi=psi / np.linalg.norm(psi), theta=res.x, nit=int(res.nit),
                E_hf=float(E_hf), n_params=P, singles=singles, doubles=doubles, eps=eps, phi=phi)


def uccsd_cnots(N, singles, doubles):
    """CNOT estimate of the UCCSD circuit on 2N qubits (block JW ordering): a single
    excitation i->a of one spin compiles into 2 Pauli strings of weight a-i+1
    (2(w-1) CNOTs each), a double into 8 strings whose weight is the span of the
    four orbitals; plus the Givens network that rotates the site register into the
    orbital basis, N(N-1)/2 rotations per spin of 2 CNOTs each."""
    cn = 0
    for (_, i, a) in singles:
        cn += 2 * 2 * (abs(a - i))
    for (kind, i, j, a, b) in doubles:
        if kind == "ud":
            idx = [i, a, j + N, b + N]
        else:
            idx = [i, j, a, b]
        w = max(idx) - min(idx) + 1
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
