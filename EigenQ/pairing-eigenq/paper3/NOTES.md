# Paper 3 — refinement + rodeo on the Hubbard model (first draft, 2026-10-11)

Files: `eigenq_hubbard.tex` / `.pdf` (RevTeX, same author list and structure as
paper 1; red `[...]` = open items), `prototype/` (sector-exact code and the
JSON/log outputs behind every number in the draft), `Hubbard.ipynb` (ED
reference, unchanged).

## What the draft establishes (prototype, sector-exact, t = 1, OBC)

- Lattice refinement: coarse site → block of fine sites, prolongation = bonding
  combination (Givens rotations on a device). Galerkin coarse model verified to
  machine precision: 1D `Hubbard(L, t/2, U/2) − tN` (self-similar at fixed U/t),
  2D `Hubbard(L², t/2, U/4) − 2tN`.
- Chain 4→8, (2,2), U/t = 0,2,4,8: p0 = 0.78/0.81/0.80/0.75; Path A (μ_buf = 2t)
  p(T=8) > 0.997, p(T=16) > 0.9999 for all U, although the fine sector gap falls
  0.53 → 0.14. Reason: spin and reflection symmetry are conserved along the
  path; the lowest excitations carry zero weight; the effective gap Δ_eff =
  0.88/0.78/0.59/0.37.
- Path B (un-dimerization, extra intra-block hopping Λ cos²): refines as fast,
  but capped at 1 − p ≈ 4.8e-3 (Λ = 8t) / 1.4e-3 (Λ = 16t) ≈ (t/Λ)² — the
  prolonged state is not the exact g.s. of H_B(0). Rodeo removes it.
- Rodeo (σ = 5/Δ_eff, exact ensemble formula): prolonged input 4–5e-3 after 6
  cycles; refined A input (T = 16) ~1e-6 after 6 cycles.
- Filled dimer 2→4 (2,2): seed is ONE basis state; p0 = 0.80 … 0.48 (U = 0 … 8),
  p(16) > 0.9997. So 2→4→8 from a single basis state, no optimizer.
- 8→16 (2,2), U = 4: p0 = 0.869, fine gap 0.051, p(16) = 0.9993, p(32) = 0.99995.
- 2D 2×2→4×4: (1,1) fine (p0 ≈ 0.8–0.85, p(8) > 0.9998). (2,2), U = 4: coarse
  g.s. (d-wave, B1) has ZERO overlap with the fine g.s. (s-wave, A1); the
  coarse first excited state (0.07t higher, A1) has p0 = 0.647 and refines to
  0.99997 at T = 16 despite a fine gap of 0.09. Lesson: refine the coarse
  state in the target's symmetry sector.

## Still to do (marked in the tex)

1. Port `prototype/hubbardrefine.py` into the project notebook / pairinglib
   (generic refine/rodeo already exist there; add lattice prolongation and the
   Galerkin coarse builder); regenerate all numbers with figures.
2. UCCSD coarse solver on the L = 4 chain (8 qubits) → matched-observable
   comparison as in paper 1 Sec. V.C.
3. Circuit level: JW Pauli counts (52 strings, weight ≤ 2 for L_f = 8), Givens
   prolongation (16 CNOTs), cost of the conditioned coarse term vs Path B,
   Trotterized refinement and controlled evolution, CNOT budgets per accepted
   preparation, JRS bound with Δ_eff.
4. Symmetry diagnostics: S² and C4v characters in the sector basis; locate the
   coarse B1–A1 crossing vs U; Δ_eff for 8→16 and the 2D (2,2) case.
5. Composed chain 2→4→8→16 using refined (not exact) inputs; single-shot 2→8.
6. Decide on a t–J / strong-coupling section (ED reference exists in the
   notebook); Path-B schedule with Λ → ∞ at s = 0; PBC/APBC variants.
7. Figures (planned): pipeline schematic; p(T) for both paths and all U; spectral
   weights vs excitation energy (symmetry protection); dimer bootstrap; 2D
   overlap decomposition; rodeo M-scans; gate budgets.

## Prototype usage

```
python3 prototype/run_chain.py   # chain 4->8, both paths, rodeo  (~6 min)
python3 prototype/run2.py        # Δ_eff + σ-matched rodeo, dimer 2->4, 2D, 8->16 (~2 min)
```
`hubbardrefine.py` reuses the notebook's `spin_basis` / `hop_spin` / kron
structure; `prolong()` builds P in the sector basis for any block map.
