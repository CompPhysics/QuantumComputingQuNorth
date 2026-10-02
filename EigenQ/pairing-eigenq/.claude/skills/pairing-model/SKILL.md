---
name: pairing-model
description: >
  Physics conventions and the validated Python API for the constant-pairing
  Hamiltonian (nuclear/condensed-matter pairing model). Use this skill WHENEVER
  the task touches the pairing model, seniority-zero spectra, the
  delta/g Hamiltonian, Jordan-Wigner encoding of pairs, Hartree-Fock pairing
  energies, FCI/CCD/UCCSD benchmarks for pairing, or the `pairinglib` package.
  Consult it before writing any pairing physics so the conventions (level
  index, qubit ordering, sign of g, benchmark numbers) stay consistent across
  the notebook, the library, and the paper.
---

# Pairing model: conventions and API

This project studies the constant-pairing Hamiltonian as a controlled testbed
for the resolution-refinement + rodeo eigenstate-preparation pipeline. Always
use the conventions below; they are baked into `src/pairinglib` and the paper.

## Hamiltonian and conventions (do not deviate)
- `H = delta * sum_p p * sum_sigma n_{p,sigma} - (g/2) * sum_{p,q} a+_{p up} a+_{p dn} a_{q dn} a_{q up}`
- `delta = 1` throughout. Level `p` (0-indexed) has single-particle energy `delta*p`.
- Qubit index `= 2*p + sigma`, with `sigma = 0` (up), `1` (down); a `k`-level
  problem uses `2k` qubits. Big-endian: qubit `j` carries bit `(n-1-j)`.
- The force conserves seniority, so the `N`-paired ground state lies entirely in
  the seniority-zero subspace; the full `N`-sector FCI equals the seniority-0 result.
- Hartree-Fock fills the lowest `N/2` levels; `E_HF(N=4, g) = 2 - g`.
- "Basis refinement" means increasing `k` at FIXED `N` (the basis/qubit count
  grows; the particle number does not).

## Validated benchmark numbers (N=4, g=1) — use as regression anchors
| k | qubits | FCI       | CCD       | UCCSD     |
|---|--------|-----------|-----------|-----------|
| 2 | 4      | 1.000000  | 1.000000  | 1.000000  |
| 3 | 6      | 0.794697  | 0.794697  | 0.794697  |
| 4 | 8      | 0.635548  | 0.630443  | 0.636987  |
| 8 | 16     | 0.145559  | 0.097263  | 0.157856  |
UCCSD is variational (>= FCI); CCD over-binds (< FCI). Singles vanish (seniority).

## The particle-hole extension (seniority-breaking; `phmodel.py`)
- `H = H_0 + V_pair + V_ph`, `V_ph = -(f/2) sum_{pqr} (a+_{p up} a+_{p dn} a_{q dn} a_{r up} + h.c.)`
  (many-body book, Chapter 4).  `q=r` pieces renormalise `g -> g+2f`; `q=p!=r` is a
  1p-1h move (seniority 2); `p,q,r` distinct is a seniority-2 2p-2h process.
- Conserves N_up and N_dn (S_z good) but NOT seniority: the one-qubit-per-pair
  encoding is gone, use the 2k-qubit register and the full N-sector.  The ground
  state stays in the S_z=0 block (36 of 70 determinants at k=4, N=4).
- Switch on with `f=`: `H_pairing_sparse(..., f=f)`, `build_H_full(..., f=f)`,
  `refine_state/refine_NK/refine_from_state(..., f=f)`, `gate_uccsd_vqe/_state(..., f=f)`.
  `f=0` reproduces the pure pairing model bit for bit.
- Reference energy `E_ref(N, g, f) = E_HF - f N/2` (= 1 - 2f at N=4, g=1).  The
  reference determinant is no longer the HF solution (Brillouin fails,
  <Phi_i^a|H|Phi_0> = +-f/2) so UCCSD singles no longer vanish.
- Circuit level: `model_terms(k, g, N, f)` = `[D] + hops` at g+2f plus the q!=r
  pieces of V_ph (`ph_terms`: 2k(k-1) 1p-1h blocks of 4 strings, k(k-1)(k-2)
  pair-breaking blocks of 8 strings, each exactly exponentiable); `trotter_U`,
  `trotter_error`, `refine_state_trotter`, `trotter_U_dt`, `cnot_counts` take `f=`.
  `run_ccd` stays pairing-only (the reference is not the HF state when f != 0).

### Validated anchors with the particle-hole term (k=4, N=4, g=1) — book Table 'pairingph'
| f    | E_ref | FCI          | UCCSD (k=4) | GS seniority-0 weight |
|------|-------|--------------|-------------|-----------------------|
| 0.00 | 1.0   |  0.63554847  | 0.636987    | 1.0000 |
| 0.05 | 0.9   |  0.45058234  | 0.453024    | 0.9988 |
| 0.20 | 0.6   | -0.18348455  | -0.175303   | 0.9869 |
| 0.50 | 0.0   | -1.69173670  | -1.659594   | 0.9605 |
Coupling to the reference (k=4, f=0.05): 4 x 2p-2h seniority-0 with (g+2f)/2,
8 x 1p-1h and 8 x seniority-2 2p-2h with +-f/2 (`ph_coupling_table`).
Pauli strings on 2k qubits (`pauli_count(build_H_full(k, g, f=f))`):
k=2: 15 -> 27 (weight 4); k=3: 34 -> 118 (weight 6); k=4: 61 -> 325 (weight 8).
First-order Trotter step (`cnot_counts(k, f)`): k=3: 150 -> 790 CNOTs; k=4: 296 -> 2728
(controlled 416 -> 3424; rotations 60 -> 348; 1p-1h blocks 640 CNOTs, 2p-2h 1792).
Trotter error (k=4, t=1, order 2, 16 steps): 2.1e-3 / 3.8e-3 / 7.5e-3 at f=0/0.2/0.5.
`tests/test_phmodel.py` asserts all of these.

## Library API (import from `pairinglib`)
See `references/library-api.md` for full signatures. Do NOT redefine these in
notebooks or scripts — import them, so there is a single tested source of truth.

## Common pitfalls
- `np.eye(2, dtype=complex)` — the 2nd positional arg of `np.eye` is the column
  count, not the dtype.
- Never shadow stdlib modules with file names (e.g. a file called `re.py`).
- The dense JW pair-interaction operator `W = sum_{p,q} A+_p A_q` is already
  Hermitian; use `H = H_kin - (g/2) W` (do not add `W + W.dag`, that double counts).
