# Paper 3 — refinement + rodeo on the Hubbard model

## Status 2026-10-14: second draft, all numbers from the project notebook

Files: `eigenq_hubbard.tex` / `.pdf` (RevTeX, 12 pp., same author list as paper I;
red `[...]` = open items), `figs/fig1-8.png` (extracted from the notebook),
`eigenq_hubbard_v1_backup.tex` (the 2026-10-11 first draft), `Hubbard.ipynb` (ED
reference, unchanged), `prototype/` (superseded sector-exact prototype; the
notebook reproduces its numbers).

Pipeline: `src/pairinglib/hubbard.py` (library, `tests/test_hubbard.py`) →
`scripts/build_notebook_hubbard.py` → `notebooks/ResolutionRefHubbard.ipynb`
(28 cells, 8 figures, ~25 min clean run) → `scripts/check_notebook_hubbard.py`
(anchors) → `make figures-hubbard` → `make paper-hubbard` (needs revtex4-2, i.e.
`texlive-publishers`).  `make hubbard` runs the chain.

## What the notebook establishes (t = 1, OBC)

- Galerkin coarse model verified to machine precision: chain `Hubbard(L, t/2, U/2) − tN`,
  square `Hubbard(L², t/2, U/4) − 2tN`.
- Chain 4→8, (2,2), U/t = 0,2,4,8 (Table I = prototype numbers): p0 0.777/0.814/0.797/0.750;
  Path A (μ_buf = 2t) p(8) > 0.997, T* = 10 (1−p < 1e-3), p(16) > 0.9999 for all U;
  Δ 0.532/0.319/0.223/0.137, Δ_eff 0.879/0.780/0.588/0.366.  Path B ceilings
  1−p(64): Λ=8 4.8–7.3e-3, Λ=16 1.7–3.7e-3.
- Symmetry protection made explicit with S² (S⁺ map), reflection parity and C4v
  characters built in the sector basis: the prolonged state has weight 1.000000000000
  in the (S=0, even) sector; the states defining the sector gap are triplets/odd with
  weights ~1e-31.  Gap-within-the-sector along the path never approaches the sector
  gap (Fig. 3).  JRS bound with Δ vs Δ_eff: 5.7e4→1.2e4 (U=0) … 1.0e7→5.8e5 (U=8) vs T*=10.
- Dimer 2→4 (Table III = prototype).  Composed 2→4→8 with refined inputs (variational μ):
  p = 0.99988 (T=16 per step) vs 0.99994 exact input; first step T=8 caps the chain at
  its input overlap 0.9979.  Single-shot 2→8 (blocks of 4): p0 0.579, p(16) 0.9966,
  p(32) 0.99988.
- 8→16 (2,2), U=4: p0 0.8692, Δ 0.0510, Δ_eff 0.1379 (weights among the lowest 40),
  p(8)=0.9924, p(16)=0.99931, p(32)=0.99995; composed 2→4→8→16 (T=16 each): 0.99923.
- 2D 2×2→4×4: (1,1) all A1, p(16) ≥ 0.99998.  (2,2): CORRECTION of the first draft —
  the fine ground state is an (A2, S=1) TRIPLET for all U>0 (checked: same energy in
  the (3,1) sector), the coarse plaquette ground state a (B1, S=0) singlet, so the
  overlap is zero by spin AND point group.  The lowest coarse (A2,1) state (0.071t up
  at U=4) carries 0.647 and refines to 0.99997 at T=16 (Δ_eff = 1.77 in that sector
  vs sector gap 0.090); from the B1 state p = 3e-27.  No B1/A2 crossing vs U in the
  coarse model for U ∈ [1,8]: the crossing is a property of the filling.
- Rodeo (exact ensemble, σ = 5/Δ_eff): prolonged 4–5e-3 at M=6; Path A T=16 input
  0.5–1.4e-6 at M=6; UCCSD L=8 input 1−p = 4.7e-4/3.1e-3/1.2e-2 (U=2/4/8).
- UCCSD (sector-level, orbital basis, exact Givens-type factors, adjoint gradient):
  coarse L=4 26 params/1720 CNOTs, 1−F 1.8e-4/3.3e-3/1.7e-2; fine L=8 198 params/
  30864 CNOTs, 1−F 4.7e-4/3.1e-3/1.2e-2.  UCCSD coarse → refined saturates at ≈ the
  coarse fidelity (0.99669 at U=4, 0.98234 at U=8).
- Circuits: L=8 chain 52 strings weight ≤ 2 (verified by full-space Pauli
  decomposition at L=4); step CNOTs 72/144/352 (order 1 / order 2 / controlled 2);
  prolongation 16; Path-A coarse term 1088 (Toffoli ladder) or 12272 (Gray code)
  per application.  Trotterised refinement T=16, U=4: n_s = 128 (A), 256 (B Λ=8),
  512 (B Λ=16) second-order steps for |p−p_exact| ≤ 1e-3; first order not converged
  at 512.  Trotterised rodeo: fixed step diagonalised once, E_obj = Trotter-effective
  eigenphase, δt = 0.1 (floor 5.8e-5; δt = 0.25 gives 2.6e-3 and E_eff−E0 = 0.09).
  Matched comparison (Table VI/Fig. 8): Path B Λ=8 3.7e4 CNOTs (p 0.9921) ~ UCCSD
  3.1e4 (p 0.9969); Path A 3.0e5 (p 0.9998, coarse term dominates); prolonged needs
  6 cycles/1.5e5 for 5e-3.  Same verdict as paper I: not a gate saver, optimizer-free.

## Still open (red in the tex)

1. 2D: U=0 labels degenerate; add an (S,Γ)-resolved spectrum at another filling
   (e.g. (3,3)) and the periodic lattice (Lieb's theorem at half filling).
2. Compile the conditioned coarse term explicitly and verify against the sector
   exponential (as for the 2→4 phase in paper I); compare JW orderings on linear
   connectivity.
3. t–J / strong-coupling section: decide.
4. PBC/APBC variants; Path-B schedule with Λ→∞ at s=0; Path B at T<16.
5. Repository link; `HubbardNotebook` bibitem authors; acknowledgments copied from
   paper I (check).
6. Possible: UCCSD coarse solver with the target's quantum numbers imposed; optimized
   refinement (nucleation-style) to cut the 1e2–1e3 Trotter steps.
