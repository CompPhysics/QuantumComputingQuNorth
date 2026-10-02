# pairinglib API reference

## hamiltonian.py
- `build_sector(k, N) -> (nq, states, index)` — fixed-N Fock sector (bitstring list + lookup).
- `H_pairing_sparse(k, g, N, states, index, delta=1.0, f=0.0) -> csr_matrix` — sparse N-sector H (`f` adds V_ph).
- `fci_ground(H) -> float` — lowest eigenvalue (exact).
- `E_HF(N, g, delta=1.0) -> float` — Hartree-Fock energy.
- `build_H_full(k, g, delta=1.0, f=0.0) -> ndarray` — dense H over the full 2^(2k) space (HEA/gate-UCCSD).

## ccd.py
- `run_ccd(k, g, N, delta=1.0, ...) -> (E, iters)` — coupled-cluster doubles.

## uccsd.py  (operator-level / statevector)
- `uccsd_pool(k, N) -> (singles, doubles)`
- `setup_uccsd(k, N, delta=1.0) -> dict(nq, states, index, ...)`
- `uccsd_vqe(setup, H, n_trotter=1, ...) -> (E, x, energy_fn, grad_fn)`

## gates.py
- `apply_1q(psi, U, q, n)`, `apply_cnot(psi, c, t, n)`, `Ry(th)`, `Rz(th)`, `Rx(th)`
- `pauli_exp(psi, phi, pauli, n)` — exp(-i phi/2 P) via basis change + CNOT ladder + Rz.

## gate_uccsd.py  (gate-level circuits)
- `make_single(i, a, n)`, `make_double(i, j, a, b, n)` — Pauli terms of an excitation.
- `apply_exc_circuit(psi, theta, terms, n)`
- `gate_uccsd_setup(k, N) -> dict(nq, terms, gens, hf, P, ncnot, nrot)`
- `gate_uccsd_vqe(k, N, g, setup, delta=1.0, f=0.0) -> (E, nit)`
- `gate_uccsd_state(k, N, g, f=0.0, delta=1.0) -> (E, statevector)`

## refine.py
- `prolong_matrix(k_low, k_high, N) -> ndarray` — gate-free basis embedding.
- `refine_state(N, k_high, T, k_low=2, M=160, mu_buf=0.6, g=1.0, f=0.0) -> (psi, H_high)`
- `refine_NK(N, k_high, Ts, ..., f=0.0) -> dict(ov, en, Eh, El, gap, ov0, E0)`
- `refine_from_state(psi_low, N, k_low, k_high, Ts, M=160, mu_buf=0.6, g=1.0, f=0.0) -> same dict`

## rodeo.py
- `rodeo_post0(v, U, t, E) -> (state, prob)` — post-selected one-cycle filter.
- `rodeo_cycle_ancilla(v, U, t, E)` — explicit ancilla circuit (same map).
- `rodeo_sweep(v0, H, ts, E) -> (state, cumulative_acceptance)`
- `rodeo_allzero_prob(v0, H, ts, E) -> float` — for energy scans.
- `rodeo_track(v0, H, ts, E, target) -> (fidelity_array, acceptance_array)`

## phmodel.py  (pairing + particle-hole term, seniority-breaking)
- `H_ph_sparse(k, f, N, states, index) -> csr_matrix` — V_ph alone (Hermitian, h.c. included) in the N-sector.
- `build_H_ph_full(k, f) -> ndarray` — dense V_ph over the full 2^(2k) space.
- `E_ref(N, g, f=0.0, delta=1.0) -> float` — reference-determinant energy, E_HF - f N/2.
- `seniority_of(I, k) -> int`, `seniority_zero_mask(k, states) -> bool array`,
  `seniority_weight(psi, k, states) -> float` — weight of a sector vector on seniority-zero determinants.
- `ph_coupling_table(k, N, g, f, delta=1.0) -> {(rank, seniority): (count, [|H|])}` — determinants reached from the reference by one action of H.

## pauli.py
- `pauli_decompose(H_dense, tol=1e-10) -> {pauli_string: coeff}` — O(n 4^n) qubit-by-qubit transform, big-endian.
- `pauli_count(H_dense, tol=1e-10, drop_identity=False) -> (n_strings, max_weight)`.

## trotter.py  (circuit-level product formulas; `f` adds the particle-hole blocks)
- `pairing_terms(k, g, N, delta=1.0) -> (D, hops)`; `ph_terms(k, f, N) -> [1p-1h blocks..., 2p-2h blocks...]`
- `model_terms(k, g, N, f=0.0, delta=1.0) -> [D] + hops (at g+2f) + ph blocks` — sums to `H_pairing_sparse(f=f)`.
- `trotter_U(k, g, N, t, n_steps, order=2, delta=1.0, terms=None, f=0.0)`, `trotter_U_dt(..., dt_max, ..., f=0.0) -> (U, n_steps)`
- `trotter_error(k, g, N, t, n_steps, order=2, delta=1.0, f=0.0, terms=None) -> float`
- `refine_state_trotter(N, k_high, T, n_s, k_low=2, mu_buf=0.6, g=1.0, order=2, delta=1.0, psi_low=None, f=0.0) -> (psi, H_high)`
- `cnot_counts(k, f=0.0) -> dict(step_uncontrolled, step_controlled, rotations)`; `ph_pauli_costs(k) -> {family: dict(terms, rotations, cnots, max_weight)}`
