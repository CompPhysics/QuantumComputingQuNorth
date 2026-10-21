"""The pair-basis (seniority-zero) module is the exact restriction of the
full-sector library, and the k=8, N=8 anchors of paper I's Sec. V.F."""
import numpy as np
import pairinglib as pl
import pairinglib.seniority as sz


def test_pair_basis_is_the_seniority_zero_block():
    k, N, g = 4, 4, 1.0
    nq, st, ix = pl.build_sector(k, N); H = pl.H_pairing_sparse(k, g, N, st, ix).toarray()
    full, pair = sz.seniority_zero_states(k, N)
    Hp = sz.H_pair(k, g, 2); blk = np.zeros_like(Hp); blk[np.ix_(pair, pair)] = H[np.ix_(full, full)]
    assert abs(blk - Hp).max() < 1e-12
    D, hops = sz.pair_terms(k, g, 2); tm = pl.model_terms(k, g, N)
    for i, T in enumerate([D] + hops):
        B = np.zeros_like(T); B[np.ix_(pair, pair)] = tm[i][np.ix_(full, full)]
        assert abs(B - T).max() < 1e-12
    assert abs(np.linalg.eigvalsh(Hp)[0] - 0.635548) < 1e-6
    # pipeline anchors of paper I reproduced in the pair basis
    r = sz.refine_pair(2, 4, 2, g, mu_buf=0.6, M=160)
    assert abs(r["p0"] - 0.8676375677) < 1e-9
    assert abs(sz.refine_pair_overlaps(r, [25.0])[0] - 0.9975712826) < 1e-9
    a = sz.pair_double_ansatz(4, 2, g)
    assert abs(a["E"] - 0.6369869810) < 1e-8 and a["cnots"] == 192
    assert abs(1 - abs(np.vdot(r["Psih"], a["psi"]))**2 - 1.7764e-4) < 1e-7


def test_k8_n8_anchors_and_bookkeeping():
    k, P, g = 8, 4, 1.0
    H = sz.H_pair(k, g, P); w, V = np.linalg.eigh(H); hf = sz.hf_state(k, P)
    assert H.shape == (70, 70) and abs(w[0] - 8.889170) < 1e-6 and abs(hf @ H @ hf - 10.0) < 1e-12
    assert abs(abs(V[:, 0] @ hf)**2 - 0.69024) < 1e-5 and abs(w[1] - w[0] - 2.5918) < 1e-4
    # response norm at s = 0 in closed form: chi_u(0) = P(k-P) g^2 / (4 b^2)
    b = 0.6; r = sz.refine_pair(4, k, P, g, mu_buf=b, M=4); A, B = r["Hle"], r["Hhs"]
    wa, Va = np.linalg.eigh(A); i0 = np.argmin(wa)
    chi = sum(abs(Va[:, n] @ (B - A) @ Va[:, i0])**2/(wa[n] - wa[i0])**2 for n in range(len(wa)) if n != i0)
    assert abs(chi - P*(k - P)*g**2/(4*b**2)) < 1e-9
    # gate bookkeeping of paper I reproduced and extended
    assert sz.coarse_term_cnots(2, 4, 2) == 14 and sz.coarse_term_cnots(3, 4, 2) == 350      # incl. the conditioned scalar phase
    assert sz.coarse_term_cnots(4, 8, 4) == 254 and sz.coarse_term_cnots(6, 8, 4) == 4886
    assert sz.coarse_term_cnots(4, 8, 4, "flag") == 44 and sz.coarse_term_cnots(3, 4, 2, "flag") == 224
    assert sz.coarse_term_cnots(6, 8, 4, "flag") == 1028
    assert pl.cnot_counts(8)["step_uncontrolled"] == 1360 and sz.pair_encoding_cnots(8)["step_uncontrolled"] == 112
    a = sz.pair_double_ansatz(k, P, g, H=H)
    assert a["n_params"] == 16 and a["cnots"] == 768 and abs(1 - abs(np.vdot(V[:, 0], a["psi"]))**2 - 3.12e-3) < 2e-5


def test_adapt_vqe_reproduces_restricted_ansatz_and_converges():
    k, P, g = 4, 2, 1.0
    H = sz.H_pair(k, g, P); GS = np.linalg.eigh(H)[1][:, 0]
    r = sz.adapt_vqe(k, P, g, pool="pairs", grad_tol=1e-6, H=H, target=GS)
    assert len(r["ops"]) == 4 and abs(r["E"] - 0.6369869810) < 1e-8          # the four-pair ansatz
    r2 = sz.adapt_vqe(k, P, g, pool="generalized", grad_tol=1e-6, H=H, target=GS)
    assert r2["history"][-1]["infidelity"] < 1e-12 and len(r2["ops"]) == 5     # exact with one more operator
    assert all(abs(np.linalg.norm(d["psi"]) - 1) < 1e-10 for d in r2["history"][1:])
