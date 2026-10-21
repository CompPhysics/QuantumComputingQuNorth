# Rodeo algorithm — formulas and conventions

Reference: K. Choi, D. Lee, J. Bonitati, Z. Qian, J. Watkins, PRL 127, 040505
(2021); Z. Qian et al., Eur. Phys. J. A 60, 151 (2024) (+ supplemental).

## One cycle (object register + 1 reused ancilla)
1. Hadamard on ancilla.
2. controlled-`e^{-iH t_k}` on the object (control = ancilla).
3. phase `P(E t_k)` on the ancilla (this shifts H -> H - E).
4. Hadamard on ancilla; measure; reuse via mid-circuit measurement.

Post-selecting ancilla |0> applies `(I + e^{iE t_k} e^{-iH t_k})/2` to the
object; an eigenstate |E_j> picks up `e^{i(E-E_j)t_k/2} cos[(E_j-E)t_k/2]`.

## Success probability
- Eigenstate input, given times:  `P = prod_k cos^2[ t_k (E_obj - E)/2 ]`   (Eq. 1)
- Gaussian times, RMS sigma:      `P = [ (1 + e^{-(E_obj-E)^2 sigma^2/2}) / 2 ]^M`  (Eq. 2)
- On resonance (E = E_obj) every factor is 1; far off resonance the per-cycle
  factor -> 1/2, so the background is `2^-M`. Energy resolution ~ 1/sigma.

## Convergence criterion (state preparation, E tuned to E_0)
- `P_M = p + sum_{k!=0} |c_k|^2 [ (1+e^{-(E_k-E_0)^2 sigma^2/2})/2 ]^M`, `p=|c_0|^2`.
- Fidelity `F_M = p / P_M ~ [1 + (1-p)/p 2^-M]^-1`; infidelity `1-F_M <~ (1-p)/p 2^-M`.
- Cycles for `1-delta`:  `M >~ log2[(1-p)/(p delta)]`.
- Acceptance `P_M -> p`, so repetitions ~ 1/p; total effort `O(|log delta|/(p eps))`.

## On hardware
`e^{-iHt}` is realised by Suzuki-Trotter over the 2k qubits (the rodeo then
filters eigenstates of the effective Trotterised Hamiltonian). The single-qubit
benchmark uses a Z-Y-Z Euler rotation `R_nhat(theta)=e^{-i theta nhat.sigma/2}`.
In simulation we apply the controlled evolution exactly inside the conserved
N-particle sector (see rodeo.py), isolating algorithmic behaviour from Trotter error.

## The averaged post-selection map of a DIGITAL filter (added 2026-10-20)

With a fixed-step product formula, `n(t) = ceil(|t|/dt)` steps of size `t/n`, the
post-selection operators `K_t = (1 + e^{iE t} U~(t))/2` do not share eigenvectors.
Because the cycle times are independent, the ensemble of accepted runs after `M`
cycles is `rho_M = K^M(rho_in) / Tr`, with `K(rho) = E_t[K_t rho K_t^+]`; the
long-run state is the dominant eigenoperator of `K`, its eigenvalue the long-run
acceptance per cycle, and the expected cost of Eq. (cost) is deterministic.
**Integrate the expectation piecewise between the boundaries |t| = n dt** (the
integrand is discontinuous there): `pl.half_gaussian_nodes(sigma, dt, m=6)`,
`pl.AveragedFilterMap(Ufun, E0, sigma, dt)` (`.track(psi, GS, M)` for the exact
ensemble infidelity/acceptance, `.dominant(GS)` for the floor) and
`pl.averaged_map_floor(...)` (floor + eigenvalue + residual + the change under
m -> 2m).  A global Gauss-Hermite rule with 60-80 nodes is NOT converged for
this integrand: at k=8, N=8, dt=0.25 it produced a spurious dominant state with
zero ground-state overlap (the converged floor is 2.3e-3 with eigenvalue 0.794).
Anchors (second order): k=4, g=1, sigma=4: dt=0.25 1.1399e-4 / 0.997821; dt=0.5
1.9345e-3 / 0.965124; dt=1 0.15394 / 0.629697.  k=8, g=1, sigma=4: dt=0.1
4.838e-5 / 0.993240; dt=0.05 3.10e-6.  k=8, g=2, sigma=1.5: dt=0.1 4.02e-4,
dt=0.05 2.73e-5 / 0.998873.  Build `Ufun` with `sz.fixed_step_evolution(terms, dt)`.
