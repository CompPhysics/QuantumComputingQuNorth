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


def test_kronecker_groups_match_dense_and_adapt():
    # the Kronecker-structured product formulas of the 16-site study reproduce the
    # dense machinery on the chain 2 -> 4 at (1,1): ramp, Path A (split coarse groups),
    # fixed-step filter track and the Trotter-effective ground state
    (Nc, bc), (Nf, bf), blocks, intra = hb.refine_chain(2)
    U = 2.0
    Hf, _ = hb.H_hubbard(Nf, bf, 1, 1, U); Hd = Hf.toarray(); wf, Vf = np.linalg.eigh(Hd); GS = Vf[:, 0]
    terms = hb.hubbard_terms(Nf, bf, 1, 1, U); groups = hb.hubbard_groups(Nf, bf, 1, 1, U)
    for g, t_ in zip(groups, terms):
        v = np.random.default_rng(0).normal(size=Hd.shape[0]) + 0j
        from scipy.linalg import expm
        assert np.linalg.norm(g.apply(0.3, v) - expm(-1j * 0.3 * t_.toarray()) @ v) < 1e-10
    pool = hb._uccsd_pool_kron(Nf, bf, 1, 1); FF = hb._to_site(pool, pool["ref"]).real
    s2 = lambda s: np.sin(np.pi * s / 2) ** 2
    dense = hb.trotter_refine(FF, [t_.toarray() for t_ in terms], [lambda s: 1.0, lambda s: 1.0, s2], 4.0, 16, order=2)
    kron = hb.trotter_path_groups(FF, groups, [lambda s: 1.0, lambda s: 1.0, s2], 4.0, 16, order=2)
    assert np.linalg.norm(dense - kron) < 1e-10
    P = hb.prolong(Nc, Nf, blocks, 1, 1, sparse=True); ps = hb.prolonged_state(Hf, P)
    tl, Ul, cst = hb.galerkin_params("chain", 1.0, U, Nc)
    Kc_e, Kc_o, Dc = [t_.toarray() for t_ in hb.hubbard_terms(Nc, bc, 1, 1, Ul, t=tl)]
    mu = ps["El"] + 2.0
    gA = [hb.ProjGroup(P, Dc + (cst - mu) * np.eye(Dc.shape[0])), hb.ProjGroup(P, Kc_e), hb.ProjGroup(P, Kc_o), groups[0], groups[1], hb.DiagGroup(groups[2].d - mu)]
    Pd = P.toarray()
    termsA = [Pd @ (Dc + (cst - mu) * np.eye(Dc.shape[0])) @ Pd.T, Pd @ Kc_e @ Pd.T, Pd @ Kc_o @ Pd.T] + [t_.toarray() for t_ in terms[:2]] + [terms[2].toarray() - mu * np.eye(Hd.shape[0])]
    c2 = lambda s: np.cos(np.pi * s / 2) ** 2
    cA = [c2] * 3 + [s2] * 3
    assert np.linalg.norm(hb.trotter_refine(ps["Phi0"], termsA, cA, 4.0, 16, order=2) - hb.trotter_path_groups(ps["Phi0"], gA, cA, 4.0, 16, order=2)) < 1e-10
    tr = hb.trotter_evolution(terms, 0.1, order=2, target=GS)
    ts = np.array([0.7, 1.3, 0.4])
    F1, A1 = hb.rodeo_track_U(ps["Phi0"].astype(complex), tr["Ufun"], ts, tr["E_eff"], GS, dt=0.1)
    F2, A2 = hb.rodeo_track_groups(ps["Phi0"], groups, ts, tr["E_eff"], GS, 0.1)
    assert np.allclose(F1, F2, atol=1e-10) and np.allclose(A1, A2, atol=1e-10)
    psi_eff, E_eff = hb.effective_step_state(groups, 0.1, GS, wf[0], 3.0, iters=40)
    assert abs(E_eff - tr["E_eff"]) < 1e-6 and abs((1 - abs(np.vdot(GS, psi_eff)) ** 2) - tr["floor"]) < 1e-6
    # sparse UCCSD agrees with the dense one; ADAPT improves monotonically and beats a random prefix
    N, bonds = hb.chain(4)
    Hs, _ = hb.H_hubbard(N, bonds, 2, 2, 4.0); w, V = np.linalg.eigh(Hs.toarray())
    dense_u = hb.uccsd_hubbard(N, bonds, 2, 2, 4.0, maxiter=300, H=Hs.toarray())
    kron_u = hb.uccsd_hubbard_kron(N, bonds, 2, 2, 4.0, maxiter=300, H=Hs)
    assert abs(dense_u["E"] - kron_u["E"]) < 1e-6 and kron_u["n_params"] == dense_u["n_params"]
    ad = hb.adapt_hubbard(N, bonds, 2, 2, 4.0, grad_tol=1e-5, max_ops=12, H=Hs, target=V[:, 0])
    inf = [d["infidelity"] for d in ad["history"]]
    assert all(b <= a + 1e-9 for a, b in zip(inf, inf[1:])) and inf[-1] < 1e-2 and ad["history"][-1]["E"] > w[0] - 1e-9


def test_galerkin_decomposition_with_unequal_site_and_particle_counts():
    """The constant of the Galerkin model is -t N_particles, not -t N_coarse: the
    embedded coarse groups of Path A must sum to P(H_low - mu)P^T also when the
    coarse site count (8) differs from the particle number (4) --- the case
    the 16-site study got wrong before the revision of 2026-10-20."""
    import numpy as np
    import pairinglib.hubbard as hb
    (Nc, bc), (Nf, bf), blocks, intra = hb.refine_chain(8)
    Hf, _ = hb.H_hubbard(Nf, bf, 2, 2, 4.0)
    P = hb.prolong(Nc, Nf, blocks, 2, 2, sparse=True); ps = hb.prolonged_state(Hf, P)
    mu = ps["El"] + 2.0
    tl, Ul, cst = hb.galerkin_params("chain", 1.0, 4.0, 2 + 2)
    Kc_e, Kc_o, Dc = [t_.toarray() for t_ in hb.hubbard_terms(8, bc, 2, 2, Ul, t=tl)]
    assert abs(Dc + (cst - mu) * np.eye(784) + Kc_e + Kc_o - (ps["Hl"] - mu * np.eye(784))).max() < 1e-12
    tl8, Ul8, cst8 = hb.galerkin_params("chain", 1.0, 4.0, 8)
    assert abs(cst8 - cst) == 4.0


def test_cnot_conventions_of_the_revision():
    import pairinglib.hubbard as hb
    # JW double-excitation weight = parity support of the two intervals, not the full span
    assert hb.uccsd_cnots(8, [], [("ud", 0, 0, 2, 2)])["excitations"] == 8 * 2 * (6 - 1)     # indices 0,2,8,10 -> weight 6
    pool = hb._uccsd_pool_kron(4, hb.chain(4)[1], 2, 2)
    lab = pool["labels"]
    assert hb.uccsd_cnots(4, [l for l in lab if len(l) == 3], [l for l in lab if len(l) == 5])["total"] == 1464
    # controlled step: the scalar single-Z sum is absorbed into the ancilla phase
    assert hb.chain_step_cnots(8, hb.chain(8)[1], 2, controlled=True)["cnots"] == 288
    assert hb.chain_step_cnots(8, hb.chain(8)[1], 2, controlled=True, absorb_scalar=False)["cnots"] == 352
    assert hb.chain_step_cnots(16, hb.chain(16)[1], 2, controlled=True)["cnots"] == 608
    # shared-flag coarse term (52 rotations at eight coarse sites, not 48)
    ca = hb.conditioned_term_cnots(8, hb.chain(8)[1], 16)
    assert ca["rotations"] == 52 and ca["flag"] == 72 + 2 * 52 + 6 * 15 + 2 and ca["frame"] == 64
    assert hb.conditioned_term_cnots(4, hb.chain(4)[1], 8)["flag"] == 124


def test_exact_digital_ensemble_matches_assessment_anchor():
    """Deterministic quantized-time ensemble at eight sites: the T=8 ramp with 64
    steps reaches 1.05045e-4 after six cycles (the sampled 40-schedule value
    was 8.9e-5, which wrongly passed the 1e-4 target)."""
    import numpy as np
    import pairinglib.hubbard as hb
    N, b = hb.chain(8); U = 4.0
    H, _ = hb.H_hubbard(N, b, 2, 2, U); w, V = np.linalg.eigh(H.toarray()); GS = V[:, 0]
    FF = np.linalg.eigh(hb.H_hubbard(N, b, 2, 2, 0.0)[0].toarray())[1][:, 0]
    Kev, Kod, D = [t_.toarray() for t_ in hb.hubbard_terms(N, b, 2, 2, U)]
    psi = hb.trotter_refine(FF, [Kev, Kod, D], [lambda s: 1.0, lambda s: 1.0, lambda s: np.sin(np.pi * s / 2) ** 2], 8.0, 64, order=2)
    tr = hb.trotter_evolution(hb.hubbard_terms(N, b, 2, 2, U), 0.1, order=2, target=GS)
    sigma = 5.0 / hb.effective_gap(FF, w, V)[0]
    out = hb.exact_digital_ensemble(tr, psi, GS, sigma, 6)
    assert abs(out["inf"][6] / 1.05045e-4 - 1) < 1e-4 and abs(out["acc"][6] - 0.99696) < 1e-4
    assert abs(out["inf"][0] - (1 - abs(np.vdot(GS, psi)) ** 2)) < 1e-12
