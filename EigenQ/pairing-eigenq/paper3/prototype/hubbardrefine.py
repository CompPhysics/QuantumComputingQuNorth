"""Resolution refinement + rodeo for the Hubbard model, sector-basis prototype.

Reuses the ED conventions of paper3/Hubbard.ipynb (up block then down block,
sorted uint64 masks, H = K_up x 1 + 1 x K_dn + U diag) and the refinement /
rodeo logic of pairinglib (refine.py, rodeo.py), generalized to arbitrary
coarse->fine lattice maps.

Lattice refinement: every coarse site c owns a block B(c) of fine sites; the
prolongation maps the coarse mode to the normalized symmetric (bonding)
combination of the block modes.  The coarse Hamiltonian is the Galerkin
projection H_low = P^dagger H_high P.
"""
import numpy as np
from math import comb
from itertools import combinations, product
from scipy.sparse import csr_matrix, coo_matrix, identity, kron, diags
from scipy.sparse.linalg import eigsh, expm_multiply

_ONE = np.uint64(1)

# ------------------------------------------------------------------ lattices
def chain(L, pbc=False):
    bonds = [(i, i+1) for i in range(L-1)]
    if pbc and L > 2: bonds.append((L-1, 0))
    return L, bonds

def square(L, pbc=False):
    site = lambda x, y: (x % L)*L + (y % L)
    bonds = []
    for x in range(L):
        for y in range(L):
            if x+1 < L or pbc: bonds.append((site(x, y), site(x+1, y)))
            if y+1 < L or pbc: bonds.append((site(x, y), site(x, y+1)))
    return L*L, bonds

def refine_chain(L, pbc=False):
    """coarse chain L -> fine chain 2L; block(c) = {2c, 2c+1}."""
    Nc, bc = chain(L, pbc); Nf, bf = chain(2*L, pbc)
    blocks = [[2*c, 2*c+1] for c in range(L)]
    intra = [(2*c, 2*c+1) for c in range(L)]
    return (Nc, bc), (Nf, bf), blocks, intra

def refine_square(L, pbc=False):
    """coarse LxL -> fine 2Lx2L; block(c=(x,y)) = 2x2 plaquette."""
    Nc, bc = square(L, pbc); Nf, bf = square(2*L, pbc)
    Lf = 2*L
    site = lambda x, y: x*Lf + y
    blocks, intra = [], []
    for x in range(L):
        for y in range(L):
            b = [site(2*x, 2*y), site(2*x+1, 2*y), site(2*x, 2*y+1), site(2*x+1, 2*y+1)]
            blocks.append(b)
            intra += [(b[0], b[1]), (b[2], b[3]), (b[0], b[2]), (b[1], b[3])]
    return (Nc, bc), (Nf, bf), blocks, intra

# ------------------------------------------------------------------ ED (notebook conventions)
def spin_basis(N, n):
    S = np.fromiter((sum(1 << p for p in c) for c in combinations(range(N), n)),
                    dtype=np.uint64, count=comb(N, n))
    S.sort(); return S

def hop_spin(N, bonds, S, amps=None):
    """-sum_b t_b (c^+_i c_j + h.c.) for one species on sorted masks S."""
    dim = len(S); rows, cols, vals = [], [], []
    cols0 = np.arange(dim, dtype=np.int64)
    for k, (i, j) in enumerate(bonds):
        tb = 1.0 if amps is None else amps[k]
        if tb == 0.0: continue
        mn, mx = min(i, j), max(i, j)
        btw = ((_ONE << np.uint64(mx)) - _ONE) ^ ((_ONE << np.uint64(mn + 1)) - _ONE)
        for (c, d) in ((i, j), (j, i)):
            bc = _ONE << np.uint64(c); bd = _ONE << np.uint64(d)
            m = ((S & bd) != 0) & ((S & bc) == 0)
            if not m.any(): continue
            Sm = S[m]
            sgn = np.bitwise_count(Sm & btw).astype(np.int64) & 1
            rows.append(np.searchsorted(S, Sm ^ (bc | bd))); cols.append(cols0[m])
            vals.append(-tb*(1.0 - 2.0*sgn))
    if not rows: return csr_matrix((dim, dim))
    return coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(dim, dim)).tocsr()

def docc(Sup, Sdn):
    return np.bitwise_count(Sup[:, None] & Sdn[None, :]).astype(np.float64).ravel()

def H_hubbard(N, bonds, Nup, Ndn, U, t=1.0, amps=None):
    Sup, Sdn = spin_basis(N, Nup), spin_basis(N, Ndn)
    a = None if amps is None else amps
    Kup = t*hop_spin(N, bonds, Sup, a); Kdn = t*hop_spin(N, bonds, Sdn, a)
    H = kron(Kup, identity(len(Sdn), format='csr'), format='csr') \
        + kron(identity(len(Sup), format='csr'), Kdn, format='csr')
    if U != 0.0: H = H + diags(U*docc(Sup, Sdn))
    return H.tocsr(), (Sup, Sdn)

def ground(H, k=2):
    dim = H.shape[0]
    if dim <= 2500:
        w, V = np.linalg.eigh(H.toarray()); return w[:k], V[:, :k]
    w, V = eigsh(H, k=k, which='SA', tol=1e-12, maxiter=20000)
    o = np.argsort(w); return w[o], V[:, o]

# ------------------------------------------------------------------ prolongation
def _perm_parity(seq):
    seq = list(seq); p = 0
    for i in range(len(seq)):
        for j in range(i+1, len(seq)):
            if seq[i] > seq[j]: p ^= 1
    return p

def prolong_spin(Nc, Nf, blocks, n):
    """One-species prolongation matrix (fine_dim x coarse_dim): coarse mode c ->
    (1/sqrt|B|) sum_{f in B(c)} c^+_f.  Sign = parity of sorting the chosen fine
    indices (product of creation operators in coarse order vs. ascending order)."""
    Sc, Sf = spin_basis(Nc, n), spin_basis(Nf, n)
    P = np.zeros((len(Sf), len(Sc)))
    for ic, mc in enumerate(Sc):
        occ = [c for c in range(Nc) if (int(mc) >> c) & 1]
        amp = np.prod([1/np.sqrt(len(blocks[c])) for c in occ]) if occ else 1.0
        for choice in product(*[blocks[c] for c in occ]):
            mf = sum(1 << f for f in choice)
            sgn = -1.0 if _perm_parity(choice) else 1.0
            P[np.searchsorted(Sf, np.uint64(mf)), ic] += sgn*amp
    return P

def prolong(Nc, Nf, blocks, Nup, Ndn):
    return np.kron(prolong_spin(Nc, Nf, blocks, Nup), prolong_spin(Nc, Nf, blocks, Ndn))

# ------------------------------------------------------------------ refinement
def slice_eigs(Hle, Hhs, M):
    out = []
    for m in range(M):
        s = (m + 0.5)/M
        out.append(np.linalg.eigh(np.cos(np.pi*s/2)**2*Hle + np.sin(np.pi*s/2)**2*Hhs))
    return out

def evolve(Phi0, eigs, T):
    psi = Phi0.astype(complex).copy(); dt = T/len(eigs)
    for lam, V in eigs:
        psi = V @ (np.exp(-1j*lam*dt)*(V.conj().T @ psi))
    return psi

def path_gap(Hle, Hhs, ns=41, sparse=False):
    gaps, s_grid = [], np.linspace(0, 1, ns)
    for s in s_grid:
        H = np.cos(np.pi*s/2)**2*Hle + np.sin(np.pi*s/2)**2*Hhs
        if sparse:
            w = np.sort(eigsh(csr_matrix(H), k=2, which='SA', tol=1e-10)[0])
        else:
            w = np.linalg.eigvalsh(H)[:2]
        gaps.append(w[1]-w[0])
    return s_grid, np.array(gaps)

def rodeo_infidelity(psi, w, V, Ms, sigma):
    """Exact ensemble infidelity (paper-1 Eq. fid) of the Gaussian-time rodeo
    filter with E_obj = w[0], for input psi expanded in eigenbasis (w, V)."""
    c2 = np.abs(V.conj().T @ psi)**2
    q = 0.5*(1 + np.exp(-(w - w[0])**2*sigma**2/2))
    out = []
    for M in Ms:
        num = np.sum(c2[1:]*q[1:]**M); out.append(num/(c2[0] + num))
    return np.array(out)
