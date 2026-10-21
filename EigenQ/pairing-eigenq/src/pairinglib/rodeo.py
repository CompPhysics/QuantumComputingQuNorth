"""Rodeo algorithm: stochastic spectral filtering (Qian et al., EPJA 60, 151)."""
import numpy as np
from scipy.linalg import expm

__all__ = ["rodeo_post0", "rodeo_cycle_ancilla", "rodeo_sweep",
           "rodeo_allzero_prob", "rodeo_track", "half_gaussian_nodes",
           "AveragedFilterMap", "averaged_map_floor"]

def rodeo_post0(v, U, t, E):
    """Post-selected |0> branch of one cycle = (I + e^{iEt} U)/2 . Returns (state, prob)."""
    w = 0.5*(v + np.exp(1j*E*t)*(U @ v))
    p = float(np.vdot(w, w).real)
    return w/np.sqrt(p), p

def rodeo_cycle_ancilla(v, U, t, E):
    """Explicit ancilla circuit (H - ctrl-U - phase - H - measure|0>); same map as rodeo_post0."""
    D = len(v); J = np.zeros((D, 2), complex); J[:, 0] = v
    a0, a1 = J[:, 0], J[:, 1]; J = np.stack([(a0+a1), (a0-a1)], axis=1)/np.sqrt(2)
    J[:, 1] = U @ J[:, 1]; J[:, 1] = J[:, 1]*np.exp(1j*E*t)
    a0, a1 = J[:, 0], J[:, 1]; J = np.stack([(a0+a1), (a0-a1)], axis=1)/np.sqrt(2)
    w = J[:, 0]; p = float(np.vdot(w, w).real)
    return w/np.sqrt(p), p

def rodeo_sweep(v0, H, ts, E):
    """Apply len(ts) cycles; return (final state, cumulative acceptance array)."""
    v = v0.astype(complex).copy(); P = 1.0; Ps = []
    for t in ts:
        v, p = rodeo_post0(v, expm(-1j*H*t), t, E); P *= p; Ps.append(P)
    return v, np.array(Ps)

def rodeo_allzero_prob(v0, H, ts, E):
    """Joint all-zero acceptance probability for one time set at target E."""
    return rodeo_sweep(v0, H, ts, E)[1][-1]

def rodeo_track(v0, H, ts, E, target):
    """Return fidelity-with-`target` and cumulative acceptance after each cycle."""
    v = v0.astype(complex).copy(); P = 1.0; F, A = [], []
    for t in ts:
        v, p = rodeo_post0(v, expm(-1j*H*t), t, E); P *= p
        F.append(abs(np.vdot(target, v))**2); A.append(P)
    return np.array(F), np.array(A)


# ---------------------------------------------------------------- averaged post-selection map of a DIGITAL filter
def half_gaussian_nodes(sigma, dt_max, m=6, xmax=8.0):
    """Quadrature nodes and weights for E_t[f(|t|)] with t ~ N(0, sigma^2), for an
    integrand that is only piecewise smooth: a fixed-step product formula uses
    n(t) = ceil(|t|/dt_max) steps, so f changes discontinuously at every
    boundary |t| = n dt_max.  The half-Gaussian density 2 phi(t/sigma)/sigma is
    therefore integrated interval by interval between those boundaries, up to
    xmax*sigma, with an m-point Gauss-Legendre rule on each interval.  A global
    Gauss-Hermite rule is NOT converged for this integrand (a few tens of
    nodes change the dominant eigenoperator of the map qualitatively at
    sixteen qubits); refining m here changes nothing beyond 1e-10."""
    xg, wg = np.polynomial.legendre.leggauss(m)
    edges = np.arange(0.0, xmax * sigma + dt_max, dt_max)
    ts, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        t = 0.5 * (b - a) * xg + 0.5 * (b + a)
        ts.append(t)
        ws.append(0.5 * (b - a) * wg * 2.0 * np.exp(-0.5 * (t / sigma) ** 2) / (sigma * np.sqrt(2.0 * np.pi)))
    return np.concatenate(ts), np.concatenate(ws)


class AveragedFilterMap:
    """The averaged post-selection map K(rho) = E_t[K_t rho K_t^+] of the rodeo
    filter with a (digital) evolution Ufun(t), K_t = (1 + e^{i E t} Ufun(t))/2,
    on resonance E = E_obj, with |t| half-Gaussian of width sigma and the
    boundary-split quadrature of `half_gaussian_nodes`.  Because the cycle
    times are independent, the ensemble of accepted runs after M cycles is
    rho_M = K^M(rho_in)/Tr K^M(rho_in), its trace the acceptance, and the
    expected number of product-formula steps per cycle is E[ceil(|t|/dt)]."""

    def __init__(self, Ufun, E_obj, sigma, dt_max, m=6, xmax=8.0):
        self.ts, self.ws = half_gaussian_nodes(sigma, dt_max, m, xmax)
        self.Ks = []
        for t in self.ts:
            U = Ufun(t)
            self.Ks.append(0.5 * (np.eye(U.shape[0]) + np.exp(1j * E_obj * t) * U))
        self.mean_steps = float(np.sum(self.ws * np.ceil(self.ts / dt_max)))
        self.dim = self.Ks[0].shape[0]

    def apply(self, rho):
        out = np.zeros_like(rho, dtype=complex)
        for Kt, wk in zip(self.Ks, self.ws):
            out += wk * (Kt @ rho @ Kt.conj().T)
        return out

    def track(self, psi0, target, M):
        """Exact ensemble infidelity, acceptance and mean cycle count, M = 0..M.
        Returns dict(inf, acc) of length M+1 (acc[0] = 1)."""
        rho = np.outer(psi0, psi0.conj()).astype(complex)
        inf = [1 - float((target.conj() @ rho @ target).real)]
        acc = [1.0]
        for _ in range(M):
            rho = self.apply(rho)
            tr = float(np.trace(rho).real)
            acc.append(tr)
            inf.append(1 - float((target.conj() @ rho @ target).real) / tr)
        return dict(inf=np.array(inf), acc=np.array(acc))

    def dominant(self, target, iters=20000, tol=1e-13):
        """Dominant eigenoperator by power iteration from |target><target|.
        Returns dict(eigenvalue, floor = 1 - <target|rho|target>, iters,
        residual = max|K(rho) - lambda rho|, and the same floor/eigenvalue with
        a doubled Gauss-Legendre order as the convergence check)."""
        rho = np.outer(target, target.conj()).astype(complex)
        ev = 0.0
        for it in range(iters):
            new = self.apply(rho)
            ev_new = float(np.trace(new).real)
            new /= ev_new
            dstate = abs(new - rho).max()
            rho = new
            if it > 5 and abs(ev_new - ev) < tol and dstate < 1e-12:
                break
            ev = ev_new
        resid = abs(self.apply(rho) - ev_new * rho).max()
        return dict(eigenvalue=ev_new, floor=1 - float((target.conj() @ rho @ target).real),
                    iters=it + 1, residual=resid, rho=rho)


def averaged_map_floor(Ufun, E_obj, target, sigma, dt_max, m=6, xmax=8.0, check=True):
    """Dominant eigenvalue (long-run acceptance per cycle) and asymptotic
    infidelity (floor) of the averaged post-selection map of a digital filter,
    with a convergence check under refinement of the quadrature (m -> 2m).
    Returns dict(eigenvalue, floor, residual, iters, floor_check, eigenvalue_check)."""
    K = AveragedFilterMap(Ufun, E_obj, sigma, dt_max, m, xmax)
    d = K.dominant(target)
    out = dict(eigenvalue=d["eigenvalue"], floor=d["floor"], residual=d["residual"], iters=d["iters"],
               mean_steps=K.mean_steps)
    if check:
        d2 = AveragedFilterMap(Ufun, E_obj, sigma, dt_max, 2 * m, xmax).dominant(target)
        out.update(floor_check=d2["floor"], eigenvalue_check=d2["eigenvalue"])
    return out
