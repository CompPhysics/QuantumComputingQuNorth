import numpy as np
import pairinglib.hubbard as hb


def test_dimer_and_free_fermions():
    N, bonds = hb.chain(2)
    for U in (0.0, 4.0, 8.0):
        H, _ = hb.H_hubbard(N, bonds, 1, 1, U)
        E0 = np.linalg.eigvalsh(H.toarray())[0]
        assert abs(E0 - 0.5 * (U - np.sqrt(U**2 + 16))) < 1e-12
    N, bonds = hb.chain(6)
    T1 = np.zeros((N, N))
    for i, j in bonds:
        T1[i, j] = T1[j, i] = -1.0
    w1 = np.linalg.eigvalsh(T1)
    H, _ = hb.H_hubbard(N, bonds, 2, 1, 0.0)
    assert abs(np.linalg.eigvalsh(H.toarray())[0] - (w1[:2].sum() + w1[0])) < 1e-12


def test_terms_sum_to_H():
    N, bonds = hb.chain(5)
    H, _ = hb.H_hubbard(N, bonds, 2, 2, 3.0)
    terms = hb.hubbard_terms(N, bonds, 2, 2, 3.0)
    assert abs(sum(T.toarray() for T in terms) - H.toarray()).max() < 1e-12


def test_galerkin_closed_form_chain_and_square():
    # chain 4 -> 8, (2,2): H_low = Hubbard(4; t/2, U/2) - t N
    (Nc, bc), (Nf, bf), blocks, intra = hb.refine_chain(4)
    U = 4.0
    Hf, _ = hb.H_hubbard(Nf, bf, 2, 2, U)
    P = hb.prolong(Nc, Nf, blocks, 2, 2)
    assert abs(P.T @ P - np.eye(P.shape[1])).max() < 1e-12
    Hl = P.T @ Hf.toarray() @ P
    tl, Ul, c = hb.galerkin_params("chain", 1.0, U, 4)
    Hc, _ = hb.H_hubbard(Nc, bc, 2, 2, Ul, t=tl)
    assert abs(Hl - (Hc.toarray() + c * np.eye(len(Hl)))).max() < 1e-12
    # square 2x2 -> 4x4, (1,1): H_low = Hubbard(2x2; t/2, U/4) - 2 t N
    (Nc, bc), (Nf, bf), blocks, intra = hb.refine_square(2)
    Hf, _ = hb.H_hubbard(Nf, bf, 1, 1, U)
    P = hb.prolong(Nc, Nf, blocks, 1, 1)
    Hl = P.T @ Hf.toarray() @ P
    tl, Ul, c = hb.galerkin_params("square", 1.0, U, 2)
    Hc, _ = hb.H_hubbard(Nc, bc, 1, 1, Ul, t=tl)
    assert abs(Hl - (Hc.toarray() + c * np.eye(len(Hl)))).max() < 1e-12


def test_chain48_anchors():
    # prototype anchors (paper3/prototype/run_chain.log), U = 4, (2,2)
    (Nc, bc), (Nf, bf), blocks, intra = hb.refine_chain(4)
    Hf, _ = hb.H_hubbard(Nf, bf, 2, 2, 4.0)
    wf, Vf = np.linalg.eigh(Hf.toarray())
    P = hb.prolong(Nc, Nf, blocks, 2, 2)
    ps = hb.prolonged_state(Hf, P)
    assert abs(ps["El"] - (-4.976573)) < 2e-6
    assert abs(wf[0] - (-5.634427)) < 2e-6
    p0 = abs(Vf[:, 0] @ ps["Phi0"]) ** 2
    assert abs(p0 - 0.797) < 2e-3
    Deff, c2 = hb.effective_gap(ps["Phi0"], wf, Vf)
    assert abs(Deff - 0.588) < 2e-3 and abs(wf[1] - wf[0] - 0.223) < 2e-3
    # <Phi0|H_high|Phi0> = E_low
    assert abs(ps["Phi0"] @ Hf @ ps["Phi0"] - ps["El"]) < 1e-10
    # path A refinement, mu_buf = 2: p(16) > 0.9999
    Hle, Hhs = hb.path_A(Hf, ps["Hl"], P, ps["El"] + 2.0)
    ov, _ = hb.refine_overlaps(ps["Phi0"], Hle, Hhs, [16.0], Vf[:, 0], M=48)
    assert ov[0] > 0.9999


def test_symmetry_diagnostics():
    (Nc, bc), (Nf, bf), blocks, intra = hb.refine_chain(4)
    Hf, _ = hb.H_hubbard(Nf, bf, 2, 2, 4.0)
    wf, Vf = np.linalg.eigh(Hf.toarray())
    s2, S = hb.total_spin(Vf[:, 0], Nf, 2, 2)
    assert abs(s2) < 1e-8                       # singlet ground state
    R = hb.site_permutation(hb.reflection_chain(Nf), Nf, 2, 2)
    assert abs(abs(Vf[:, 0] @ (R @ Vf[:, 0])) - 1.0) < 1e-8   # reflection eigenstate
    assert abs((R @ R).toarray() - np.eye(R.shape[0])).max() < 1e-12   # involution
    # C4v: a (1,1) state on 2x2 at U=0 -> the ground state is A1
    N, bonds = hb.square(2)
    H, _ = hb.H_hubbard(N, bonds, 1, 1, 0.0)
    w, V = np.linalg.eigh(H.toarray())
    lab, wt = hb.c4v_irrep(V[:, 0], 2, 1, 1)
    assert lab == "A1" and abs(wt - 1.0) < 1e-10
    wts, _ = hb.c4v_weights(V[:, 0], 2, 1, 1)
    assert abs(sum(wts.values()) - 1.0) < 1e-10


def test_rodeo_formula_and_trotter():
    (Nc, bc), (Nf, bf), blocks, intra = hb.refine_chain(2)   # 2 -> 4, (1,1)
    Hf, _ = hb.H_hubbard(Nf, bf, 1, 1, 2.0)
    Hd = Hf.toarray(); wf, Vf = np.linalg.eigh(Hd)
    P = hb.prolong(Nc, Nf, blocks, 1, 1); ps = hb.prolonged_state(Hf, P)
    inf, acc = hb.rodeo_infidelity(ps["Phi0"], wf, Vf, [1, 4], 4.0)
    assert inf[1] < inf[0] < 1 - abs(Vf[:, 0] @ ps["Phi0"]) ** 2
    # Trotterised evolution converges to the exact one
    terms = hb.hubbard_terms(Nf, bf, 1, 1, 2.0)
    from scipy.linalg import expm
    v = ps["Phi0"].astype(complex)
    ex = expm(-1j * Hd * 1.5) @ v
    def err(dt, order):
        n = int(round(1.5 / dt))
        tr = hb.trotter_evolution(terms, dt, order=order, target=Vf[:, 0])
        assert abs(tr["tq"](1.5) - 1.5) < 1e-12
        return np.linalg.norm(tr["Ufun"](n)(v) - ex)
    e1 = [err(dt, 1) for dt in (0.1, 0.05)]; e2 = [err(dt, 2) for dt in (0.1, 0.05)]
    assert 1.8 < e1[0] / e1[1] < 2.2 and 3.6 < e2[0] / e2[1] < 4.4   # first / second order
    tr = hb.trotter_evolution(terms, 0.05, order=2, target=Vf[:, 0])
    assert tr["floor"] < 1e-4 and abs(tr["E_eff"] - wf[0]) < 1e-2        # Trotter-effective ground state
    # UCCSD on the 4-site chain, (2,2): below HF, above FCI, singlet
    N, bonds = hb.chain(4)
    H, _ = hb.H_hubbard(N, bonds, 2, 2, 4.0)
    w, V = np.linalg.eigh(H.toarray())
    r = hb.uccsd_hubbard(N, bonds, 2, 2, 4.0, maxiter=300)
    assert r["E_hf"] > r["E"] > w[0] - 1e-9
    assert r["E"] - w[0] < 0.05 and abs(hb.total_spin(r["psi"], N, 2, 2)[0]) < 1e-5


def test_uccsd_pruning_and_seed():
    # keep= restricts the pool (labels and CNOT count follow), theta0= seeds the optimizer
    N, bonds = hb.chain(4)
    H, _ = hb.H_hubbard(N, bonds, 2, 2, 4.0); Hd = H.toarray()
    full = hb.uccsd_hubbard(N, bonds, 2, 2, 4.0, maxiter=300, H=Hd)
    keep = [i for i, th in enumerate(full["theta"]) if abs(th) > 1e-3]
    pr = hb.uccsd_hubbard(N, bonds, 2, 2, 4.0, maxiter=300, H=Hd, keep=keep, theta0=full["theta"][keep])
    assert pr["n_params"] == len(keep) < full["n_params"]
    assert len(pr["singles"]) + len(pr["doubles"]) == len(keep)
    assert hb.uccsd_cnots(N, pr["singles"], pr["doubles"])["total"] < hb.uccsd_cnots(N, full["singles"], full["doubles"])["total"]
    assert abs(pr["E"] - full["E"]) < 1e-4                       # the dropped amplitudes were negligible
    assert abs(np.linalg.norm(pr["psi"]) - 1.0) < 1e-10
