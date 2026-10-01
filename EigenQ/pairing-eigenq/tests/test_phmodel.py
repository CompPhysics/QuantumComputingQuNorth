"""Particle-hole extension: anchors from Chapter 4 of the many-body book
(Table 'pairingph', Table 'pairingphcoupling', Table '4-paulicounts')."""
import numpy as np
import pairinglib as pl

BOOK = {0.00: 0.63554847, 0.05: 0.45058234, 0.20: -0.18348455, 0.50: -1.69173670}


def test_book_energies_k4():
    nq, st, ix = pl.build_sector(4, 4)
    for f, ref in BOOK.items():
        H = pl.H_pairing_sparse(4, 1.0, 4, st, ix, f=f)
        assert abs(pl.fci_ground(H) - ref) < 1e-7
        assert abs(pl.E_ref(4, 1.0, f) - (1.0 - 2*f)) < 1e-12


def test_f_zero_is_pure_pairing():
    nq, st, ix = pl.build_sector(3, 4)
    H0 = pl.H_pairing_sparse(3, 1.0, 4, st, ix)
    H1 = pl.H_pairing_sparse(3, 1.0, 4, st, ix, f=0.0)
    assert abs(H0 - H1).max() == 0.0
    assert abs(pl.build_H_full(3, 1.0) - pl.build_H_full(3, 1.0, f=0.0)).max() == 0.0


def test_full_space_matches_sector():
    nq, st, ix = pl.build_sector(3, 4)
    Hf = pl.build_H_full(3, 1.0, f=0.2)
    Hs = pl.H_pairing_sparse(3, 1.0, 4, st, ix, f=0.2).toarray()
    assert np.max(np.abs(Hf[np.ix_(st, st)] - Hs)) < 1e-12
    assert np.max(np.abs(Hf - Hf.T)) < 1e-12          # Hermitian


def test_coupling_table_and_seniority():
    # f = 0: only the four seniority-zero 2p-2h states, element g/2
    t0 = pl.ph_coupling_table(4, 4, 1.0, 0.0)
    assert t0 == {(2, 0): (4, [0.5])}
    # f = 0.05: 8 1p-1h and 8 seniority-two 2p-2h with f/2; pairing block (g+2f)/2
    t1 = pl.ph_coupling_table(4, 4, 1.0, 0.05)
    assert t1[(2, 0)] == (4, [0.55]) and t1[(1, 2)] == (8, [0.025]) and t1[(2, 2)] == (8, [0.025])
    nq, st, ix = pl.build_sector(4, 4)
    H = pl.H_pairing_sparse(4, 1.0, 4, st, ix, f=0.2).toarray()
    gs = np.linalg.eigh(H)[1][:, 0]
    w = pl.seniority_weight(gs, 4, st)
    assert 0.0 < w < 1.0                                # pairs are broken
    H0 = pl.H_pairing_sparse(4, 1.0, 4, st, ix).toarray()
    assert abs(pl.seniority_weight(np.linalg.eigh(H0)[1][:, 0], 4, st) - 1.0) < 1e-12


def test_pauli_counts_book_table():
    assert pl.pauli_count(pl.build_H_full(2, 1.0)) == (15, 4)
    assert pl.pauli_count(pl.build_H_full(2, 1.0, f=0.3)) == (27, 4)
    assert pl.pauli_count(pl.build_H_full(3, 1.0)) == (34, 4)
    assert pl.pauli_count(pl.build_H_full(3, 1.0, f=0.3)) == (118, 6)
    # decomposition reconstructs H
    H = pl.build_H_full(2, 1.0, f=0.3); n = 4
    from pairinglib._exc import pauli_mat
    R = sum(c*pauli_mat({q: ch for q, ch in enumerate(s) if ch != "I"}, n)
            for s, c in pl.pauli_decompose(H).items())
    assert np.max(np.abs(R - H)) < 1e-10


def test_refinement_with_ph_term():
    r = pl.refine_NK(4, 4, np.array([0.5, 25.0]), f=0.2)
    assert abs(r['E0'] - r['El']) < 1e-9             # prolongation keeps E_low
    assert abs(r['Eh'] - BOOK[0.20]) < 1e-7
    assert r['ov0'] < r['ov'][-1] and r['ov'][-1] > 0.99
    assert r['en'][-1] >= r['Eh'] - 1e-6


def test_uccsd_singles_switch_on():
    su = pl.setup_uccsd(4, 4)
    H = pl.H_pairing_sparse(4, 1.0, 4, su['states'], su['index'], f=0.2)
    E, x, _, _ = pl.uccsd_vqe(su, H)
    assert E >= BOOK[0.20] - 1e-9 and E - BOOK[0.20] < 0.02       # variational, close
    assert np.max(np.abs(x[:su['nS']])) > 1e-4                   # singles no longer vanish
