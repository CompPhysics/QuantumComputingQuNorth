import numpy as np
from scipy.linalg import expm
import pairinglib as pl

def test_ancilla_equals_compact_filter():
    rng = np.random.default_rng(0); H = np.diag([0.0, 1.0, 2.5]); n = 3
    v = rng.normal(size=n) + 1j*rng.normal(size=n); v /= np.linalg.norm(v)
    U = expm(-1j*H*0.7)
    a, pa = pl.rodeo_cycle_ancilla(v.copy(), U, 0.7, 0.5)
    b, pb = pl.rodeo_post0(v.copy(), U, 0.7, 0.5)
    assert np.max(np.abs(a - b)) < 1e-12 and abs(pa - pb) < 1e-12

def test_eq2_gaussian_average():
    sigma, dE = 3.0, 0.5
    ts = np.random.default_rng(1).normal(0, sigma, 400000)
    emp = np.mean(np.cos(dE*ts/2)**2)
    ana = 0.5*(1 + np.exp(-dE**2*sigma**2/2))
    assert abs(emp - ana) < 5e-3

def test_rodeo_polishes_refined_state():
    psi, Hh = pl.refine_state(4, 4, 25.0)
    w, V = np.linalg.eigh(Hh); E0, GS = w[0], V[:, 0]
    p0 = abs(np.vdot(GS, psi))**2
    ts = np.random.default_rng(2).normal(0, 3.0, 6)
    F, A = pl.rodeo_track(psi, Hh, ts, E0, GS)
    assert F[-1] > p0 and F[-1] > 0.999      # fidelity improves toward 1
    assert abs(A[-1] - p0) < 5e-3            # acceptance -> p


def test_averaged_map_boundary_split_quadrature_is_converged():
    """The averaged post-selection map of the digital filter (k=4, N=4, g=1,
    sigma=4, second order): the boundary-split quadrature reproduces the
    anchors to all quoted digits and is insensitive to its own refinement,
    whereas a global Gauss-Hermite rule is not (60 vs 80 nodes differ)."""
    import numpy as np
    import pairinglib as pl
    import pairinglib.seniority as sz
    H = sz.H_pair(4, 1.0, 2); w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]
    D, hops = sz.pair_terms(4, 1.0, 2)
    ts, ws = pl.half_gaussian_nodes(4.0, 0.25)
    assert abs(ws.sum() - 1) < 1e-12 and ts.min() > 0
    for dt, ev_ref, fl_ref in ((0.25, 0.99782077, 1.13988e-4), (0.5, 0.96512446, 1.93451e-3), (1.0, 0.62969672, 0.153936)):
        out = pl.averaged_map_floor(sz.fixed_step_evolution([D] + hops, dt), E0, GS, 4.0, dt)
        assert abs(out["eigenvalue"] - ev_ref) < 1e-7 and abs(out["floor"] / fl_ref - 1) < 2e-5
        assert abs(out["floor_check"] / out["floor"] - 1) < 1e-5 and out["residual"] < 1e-10
    # the unconverged global rule, for the record
    from numpy.polynomial.hermite_e import hermegauss
    Uf = sz.fixed_step_evolution([D] + hops, 0.5)
    fl = []
    for n in (60, 80):
        x, wq = hermegauss(n); wq = wq / np.sqrt(2 * np.pi)
        K = pl.AveragedFilterMap.__new__(pl.AveragedFilterMap)
        K.ts, K.ws = np.abs(4.0 * x), wq
        K.Ks = [0.5 * (np.eye(6) + np.exp(1j * E0 * t) * Uf(t)) for t in K.ts]
        fl.append(K.dominant(GS)["floor"])
    assert abs(fl[0] - fl[1]) > 1e-4


def test_averaged_map_superoperator_matches_node_loop():
    """The cached d^2 x d^2 superoperator of the averaged map gives the same
    exact-channel track (infidelity and acceptance) as the node loop, and the
    half-normal (positive-time) convention is the one implemented."""
    import numpy as np
    import pairinglib as pl
    import pairinglib.seniority as sz
    H = sz.H_pair(4, 1.0, 2); w, V = np.linalg.eigh(H); E0, GS = w[0], V[:, 0]
    D, hops = sz.pair_terms(4, 1.0, 2)
    K = pl.AveragedFilterMap(sz.fixed_step_evolution([D] + hops, 0.25), E0, 4.0, 0.25)
    assert K.ts.min() > 0                                   # t = |z|, z ~ N(0, sigma^2)
    psi = sz.hf_state(4, 2).astype(complex)
    a = K.track(psi, GS, 6, superop=False); b = K.track(psi, GS, 6, superop=True); c = K.track(psi, GS, 6)
    assert np.allclose(a["inf"], b["inf"], rtol=1e-10, atol=1e-14) and np.allclose(a["acc"], b["acc"], rtol=1e-12)
    assert np.allclose(a["inf"], c["inf"], rtol=1e-10, atol=1e-14)
    S = K.superoperator(); rho = np.outer(psi, psi.conj())
    assert np.allclose((S @ rho.reshape(-1)).reshape(6, 6), K.apply(rho), atol=1e-13)
