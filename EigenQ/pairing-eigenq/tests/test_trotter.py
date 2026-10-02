import numpy as np
import pairinglib as pl


def _H(k, g, N):
    nq, states, index = pl.build_sector(k, N)
    return pl.H_pairing_sparse(k, g, N, states, index).toarray()


def test_terms_sum_to_hamiltonian():
    for k in (3, 4):
        D, hops = pl.pairing_terms(k, 1.0, 4)
        H = _H(k, 1.0, 4)
        assert np.allclose(D + sum(hops), H, atol=1e-12)


def test_trotter_error_scaling():
    # order-1 ~ 1/n, order-2 ~ 1/n^2
    e1a = pl.trotter_error(4, 1.0, 4, 1.0, 4, order=1)
    e1b = pl.trotter_error(4, 1.0, 4, 1.0, 8, order=1)
    e2a = pl.trotter_error(4, 1.0, 4, 1.0, 4, order=2)
    e2b = pl.trotter_error(4, 1.0, 4, 1.0, 8, order=2)
    assert e1a / e1b > 1.7          # ~2
    assert e2a / e2b > 3.4          # ~4
    assert e2a < e1a


def test_trotterised_refinement_reaches_high_overlap():
    psi, Hh = pl.refine_state_trotter(4, 4, 30.0, 120, order=2)
    w, v = np.linalg.eigh(Hh)
    p = abs(np.vdot(v[:, 0], psi)) ** 2
    assert p > 0.99


def test_rodeo_with_trotterised_evolution_converges():
    # fixed Trotter STEP SIZE (steps per cycle ~ t/dt, as on hardware)
    psi, Hh = pl.refine_state_trotter(4, 4, 30.0, 120, order=2)
    w, v = np.linalg.eigh(Hh)
    gs = v[:, 0]
    E0 = w[0]
    rng = np.random.default_rng(7)
    ts = np.abs(rng.normal(0.0, 4.0, 6))
    D, hops = pl.pairing_terms(4, 1.0, 4)
    terms = [D] + hops

    def Ufun(t):
        return pl.trotter_U_dt(4, 1.0, 4, t, 0.25, order=2, terms=terms)[0]

    F, A = pl.rodeo_track_U(psi, Ufun, ts, E0, gs)
    assert F[-1] > 1 - 1e-3       # floor set by the O(dt^2) rotation of H_eff's GS
    assert A[-1] > 0.9


def test_refine_from_state_matches_refine_NK_for_exact_input():
    nl, sl, il = pl.build_sector(3, 4)
    Hl = pl.H_pairing_sparse(3, 1.0, 4, sl, il).toarray()
    wl, vl = np.linalg.eigh(Hl)
    r1 = pl.refine_from_state(vl[:, 0], 4, 3, 4, np.array([10.0]))
    r2 = pl.refine_NK(4, 4, np.array([10.0]), k_low=3)
    assert abs(r1["ov"][0] - r2["ov"][0]) < 1e-9
    assert abs(r1["ov0"] - r2["ov0"]) < 1e-9


def test_cnot_counts_match_paper():
    c = pl.cnot_counts(4)
    assert c["step_uncontrolled"] == 296     # 24 k(k-1) + 2k at k=4
    assert c["rotations"] == 60
    assert c["step_controlled"] == 416       # + 2 CNOTs per controlled rotation


# ---- particle-hole term at the circuit level ----
def test_model_terms_sum_to_H_with_ph():
    N = 4
    for k, f in ((3, 0.2), (4, 0.5)):
        nq, st, ix = pl.build_sector(k, N)
        H = pl.H_pairing_sparse(k, 1.0, N, st, ix, f=f).toarray()
        tm = pl.model_terms(k, 1.0, N, f)
        assert len(tm) == 1 + k*(k-1)//2 + 2*k*(k-1) + k*(k-1)*(k-2)
        assert np.max(np.abs(sum(tm) - H)) < 1e-12
    assert len(pl.model_terms(4, 1.0, N, 0.0)) == 7          # f=0: pairing blocks only


def test_ph_blocks_are_exactly_exponentiable():
    # every V_ph block is a sum of mutually commuting Pauli strings
    from pairinglib._exc import pauli_mat
    k = 3; n = 2*k
    for (p, q, r) in ((0, 0, 1), (1, 2, 1), (0, 1, 2), (2, 0, 1)):
        dec = pl.pauli_decompose(pl.ph_term_full(k, 0.3, p, q, r))
        mats = [pauli_mat({i: c for i, c in enumerate(s) if c != "I"}, n) for s in dec]
        assert len(dec) in (4, 8)
        assert max(np.abs(A @ B - B @ A).max() for A in mats for B in mats) < 1e-12


def test_trotter_error_with_ph_decreases():
    e = [pl.trotter_error(4, 1.0, 4, 1.0, n, order=2, f=0.2) for n in (1, 4, 16)]
    assert e[0] > e[1] > e[2] and e[2] < 5e-3


def test_cnot_counts_with_ph():
    c0, c1 = pl.cnot_counts(4), pl.cnot_counts(4, f=0.2)
    assert (c0['step_uncontrolled'], c0['step_controlled']) == (296, 416)
    assert (c1['step_uncontrolled'], c1['step_controlled'], c1['rotations']) == (2728, 3424, 348)
    pc = pl.ph_pauli_costs(4)
    assert pc['1p1h']['terms'] == 24 and pc['2p2h']['terms'] == 24
    assert pc['1p1h']['rotations'] == 96 and pc['2p2h']['rotations'] == 192


def test_refine_state_trotter_with_ph():
    N = 4; nq, st, ix = pl.build_sector(4, N)
    H = pl.H_pairing_sparse(4, 1.0, N, st, ix, f=0.2).toarray(); GS = np.linalg.eigh(H)[1][:, 0]
    psi, _ = pl.refine_state_trotter(N, 4, 30.0, 80, order=2, f=0.2)
    assert abs(np.vdot(GS, psi))**2 > 0.99
