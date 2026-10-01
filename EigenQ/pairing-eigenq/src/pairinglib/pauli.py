"""Pauli-string decomposition of a dense qubit Hamiltonian.

    H = sum_P c_P P,   c_P = Tr(P H) / 2^n,   P in {I,X,Y,Z}^n .

The transform is applied qubit by qubit (a 4x4 change of basis on each
(row-bit, column-bit) pair), so the cost is O(n 4^n) instead of 8^n: fine up
to the n = 8 qubits of the k = 4 problem.  Qubit 0 is the most significant
bit (big-endian), as everywhere in the package.
"""
import numpy as np

__all__ = ["pauli_decompose", "pauli_count"]

_I = np.eye(2, dtype=complex)
_X = np.array([[0, 1], [1, 0]], complex)
_Y = np.array([[0, -1j], [1j, 0]], complex)
_Z = np.array([[1, 0], [0, -1]], complex)
_LETTERS = "IXYZ"
# W[sigma, a*2+b] = P_sigma[b, a] / 2   so that  c = sum_{ab} W[sigma,(a,b)] H[a,b]
_W = np.array([[P[b, a]/2 for a in range(2) for b in range(2)]
               for P in (_I, _X, _Y, _Z)])


def pauli_decompose(H, tol=1e-10):
    """Return a dict {pauli_string: real coefficient} for a Hermitian dense H."""
    H = np.asarray(H, dtype=complex); dim = H.shape[0]
    n = int(round(np.log2(dim))); assert 2**n == dim
    T = H.reshape((2,)*(2*n))                        # axes a_0..a_{n-1}, b_0..b_{n-1}
    perm = [ax for q in range(n) for ax in (q, n + q)]  # interleave (a_q, b_q)
    T = np.transpose(T, perm).reshape((4,)*n)       # axis q holds (a_q, b_q)
    for q in range(n):                               # (a_q,b_q) -> sigma_q
        T = np.moveaxis(np.tensordot(_W, T, axes=([1], [q])), 0, q)
    out = {}
    for idx in zip(*np.nonzero(np.abs(T) > tol)):
        out["".join(_LETTERS[i] for i in idx)] = float(np.real(T[idx]))
    return out


def pauli_count(H, tol=1e-10, drop_identity=False):
    """(number of Pauli strings, largest Pauli weight) of H."""
    terms = pauli_decompose(H, tol)
    if drop_identity: terms.pop("I"*len(next(iter(terms))), None)
    weights = [sum(ch != "I" for ch in s) for s in terms]
    return len(terms), (max(weights) if weights else 0)
