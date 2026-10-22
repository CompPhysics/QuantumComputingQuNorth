# Paper 3 — refinement + rodeo on the Hubbard model

## Status 2026-10-22: sixth draft — final pre-submission pass (the paper-I checklist of `paper/morefeedback.md`)

Response: `response_to_final_pass.md`.  Backups of the fifth draft: `eigenq_hubbard_v5_backup.tex`,
`eigenq_hubbard_supp_v5_backup.tex`.  Main change: the operating-point search behind Tables 2–3
was not converged (coarse T grids).  New notebook Section 14 (EXT:HREFINE,
`scripts/hubbard_refine_checks.py`, self-contained; inserted by
`scripts/extend_notebook_hubbard_refine.py`, emitted by the builder; fig 15 = the cost frontiers,
now the main-text figure; figs 13–14 moved to SM Sec. S9 with Table S7): fine T grids, geometric
n_s (ratio 1.15) then integer by bisection, M=0..6 with the exact ensemble at eight sites, M=0 at
sixteen (complete: one cycle 1.76e5 > every unfiltered minimum).  Minima (CNOTs per accepted
preparation): eight sites 1e-3/3e-4/1e-4 — ramp 1.0e4/1.6e4/2.4e4, Path A 3.0e4/5.3e4/8.5e4
(×2.95–3.53), ADAPT 1.7e4/2.3e4/7.9e4, pruned UCCSD 5.1e4/9.0e4/–, Path B 9.8e4/1.4e5/2.1e5;
sixteen sites 1e-2/3e-3/1e-3/3e-4 — ramp 2.9e4/4.4e4/7.5e4/1.2e5, Path A 4.6e4/7.2e4/1.2e5/2.2e5
(×1.60–1.86; coarse states 16 416 / 24 192 CNOTs), ADAPT 4.8e4/3.9e5/7.9e5/1.1e6, pruned
6.2e4/4.1e5/5.9e5/9.8e5, Path B 1.1e5/2.7e5/5.3e5/9.2e5.  The "tie at 1e-3 at thirty-two
qubits" was a T-spacing artifact.  Also: positive quantized-time protocol stated, 1/p scoped
(Lin–Tong), availability section, SM built against `eigenq_hubbard_xr.aux` (hyperref+xr failed
on the imported \bibcite entries), `make check-hubbard` → 15 figures.  Library:
`hb.exact_digital_ensemble(cache=)`, complex eigenvectors in `hb._eigterms`.


## Status 2026-10-20: fifth draft — revision after `hubbard_pre_submission_critical_assessment.md`

Response: `response_to_pre_submission_assessment.md` (verdict per point + the rebuilt
numbers).  Backups of the fourth draft: `eigenq_hubbard_v4_backup.tex`,
`eigenq_hubbard_supp_v4_backup.tex`.  All ten points of the assessment were agreed with
and every numerical claim in it reproduced before the changes.

What changed in the code (`src/pairinglib/hubbard.py`, tests in `tests/test_hubbard.py`,
11 pass): the Galerkin constant of the embedded Path-A term at 8→16 is −t·N_particles
(`galerkin_params("chain", 1, U, 2+2)`; the old −8t spoiled every 32-qubit Path-A
number); `conditioned_term_cnots` returns the shared-flag cost (`flag` + `frame`:
156 / 332 per application for 4→8 / 8→16, 52 rotations at eight coarse sites) besides
the ladder (1132 / 4948); `chain_step_cnots(controlled=True)` absorbs the scalar single-Z
sum (288 / 608); `uccsd_cnots` uses parity-support weights (1464 / 23 312 / 235 904;
pruned 11 408 / 61 584); `exact_digital_ensemble` (ρ_ij = c_i c_j* Q_ij^M) gives the
eight-site filter exactly; `effective_step_state(return_residual=True)`.

Notebook: Section 13 (cells 32–35, figs 13–14; ~3.5 h on top of Sections 1–12) is the
joint operating-point search — (T, n_s, M) for the ramp, (Λ, T, n_s, M) for Path B,
(T, n_s, M) for Path A, M = 0 allowed, every ADAPT prefix — with the actually prepared
coarse state charged and propagated (two options per level: 4-site ramp 8 216 CNOTs /
6.7e-5 and 16 408 / 1.4e-5; 8-site ramp 27 760 / 1.8e-4 and 36 976 / 5.6e-5).  Sections
11–12 keep the old fixed-input comparisons (exact-coarse oracle, sampled filter) as
diagnostics.  36 cells, 14 figures; `check_notebook_hubbard.py ... 14`.  Long runs:
`scripts/run_notebook_checkpoint.py` (dill checkpoint after Section 12, `--resume-from
16 --source <rebuilt nb>`).

Results (Tables 2–3 of the paper, SM Tables S5–S6): see the response doc.  Headline:
the cheapest route from 3e-3 down at both sizes is an UNFILTERED ramp from the
determinant (8 sites 1.4e4 / 2.8e4 / 3.7e4 at 1e-3 / 3e-4 / 1e-4; 16 sites 7.8e4 /
7.8e4 / 1.56e5 / 1.56e5 at 1e-2 / 3e-3 / 1e-3 / 3e-4); Path A with the shared flag is
×2.4–3.1 at 8 sites (was ×16) and ties the ramp at 1e-3 at 16 sites (1.52e5 vs
1.56e5), ×1.15–1.9 elsewhere; ADAPT is cheapest only where it needs no cycle (1e-3,
3e-4 at 8 sites; 1e-2 at 16); Path B ×6–9 / ×1.4–5.9; filter cycles are not needed by
the adiabatic routes.  Bootstrap of the operating-point selection at 16 sites (400
resamples of 40 schedules): ramp/Path A 100 % stable.  Main text figures now 1, 2, 3,
5, 11, 13, 14 (10 and 12 superseded, kept in the notebook).

Still open: `\bibitem{Code}` URL/DOI, `Hubbard.ipynb` authors, `make paper-hubbard`
on the Mac (quantumarticle is not in the container; a shim compile passes without
errors or undefined references).

## Status 2026-10-16 (superseded by the 2026-10-20 revision; the 32-qubit Path-A numbers below carry the −8t bug): fourth draft — `quantumarticle` (Quantum), main text + Supplemental Material, 8→16 at the circuit level, ADAPT-VQE

Files: `eigenq_hubbard.tex` / `.pdf` (main text, `quantumarticle` two-column, 16 pp.
incl. refs; `booktabs` tables; unstructured abstract; cites paper I only — no
"companion"/PHcompanion references), `eigenq_hubbard_supp.tex` / `.pdf` (Supplemental
Material, Secs. S1–S8, Tables S1–S7, Figs. S1–S4; `xr` against `eigenq_hubbard.aux`, so
compile the main text first — `make paper-hubbard` does both), `figs/fig1-12.png`
(main text: 1, 2, 3, 5, 10, 11, 12; SM: 4, 6, 7, 9; fig8 superseded),
`eigenq_hubbard_v3_backup.tex` (the 2026-10-15 RevTeX/PRC third draft).

Structure (paper-I template): Intro with contributions (a) reference implementation
(b) matched benchmark at 8 AND 16 sites (c) transferable results — self-similar
Galerkin hierarchy, symmetry-protected Δ_eff and the wrong-sector diagnostics, filter
cost ∝ 1/Δ_eff with a fixed-step floor; Sec. 2 model/ED/diagnostics; Sec. 3 pipeline
(+ comparators: determinant, U-ramp, UCCSD, pruned, ADAPT-VQE); Sec. 4 sector-exact
results (4→8, symmetry protection, dimer hierarchy to 16 sites, plaquette, filter +
variational summary); Sec. 5 circuits + matched comparison at 8 sites (Table 2,
Fig. 10); Sec. 6 NEW: 8→16 at the circuit level (Fig. 11 transport + ADAPT, Table 3,
Fig. 12); Sec. 7 discussion.  SM: S1 JRS bound + response norm (Table S1, Fig. S1);
S2 dimer step + composition (Table S2, Fig. S2) + 8→16 exact; S3 plaquette (Table S3);
S4 UCCSD tables (Table S4) + ADAPT details + 8-site ADAPT through the filter; S5
Trotterized refinement at 8 sites (Fig. S3) + coarse-term cost conventions; S6 exact
ensemble filter (Fig. S4); S7 ingredients of the 32-qubit study (Kronecker product
formulas, effective step state, gate counts, step scans, ramp gap); S8 operating
points (Table S5 8 sites; Table S7 full 16-site data).

Notebook Section 12 (cells 30–31, figs 11–12; ~1.5 h on top of Sections 1–11; built
on `hb.hubbard_groups/LayerGroup/DiagGroup/ProjGroup/trotter_path_groups/
effective_step_state/rodeo_track_groups/uccsd_hubbard_kron/adapt_hubbard`), chain
8→16 at (2,2), U=4 (784→14 400 states):
- inputs: prolonged coarse GS p0 = 0.8692, determinant p0 = 0.7162, both Δ_eff =
  0.1379; sector gap 0.0510; gap along the U-ramp 0.165→0.051.
- exact slices: ramp p(4..64) = 0.82856/0.90373/0.94636/0.96970/0.98970/0.99591/
  0.99887/0.99919 (T=64 for 1e-3); Path A μ_buf=2: p(8/16/32) = 0.99237/0.99931/
  0.99995 (μ_buf=4: 0.99155/0.99879/0.99991; larger shift does not help); Path B
  Λ=8: p(16/32) = 0.99492/0.99531, Λ=16: 0.99500/0.99605 (ceiling, as at 8 sites).
- gates: 16-site step 304 (controlled 736), coarse term 4856 (Toffoli ladder, 16
  controls), Givens 480, prolongation 32.  Trotter steps: ramp T=16 n_s=128 (0.96927,
  39 392 CNOTs), T=32 n_s=256 (0.99556, 78 304; n_s=32 → 0.005!); Path A split T=16
  n_s=128 (0.99835, 1 282 080 + coarse prep), T=8 plateaus at 0.99079 (exact 0.99237);
  Path B Λ=8 T=16 n_s=256 (0.99476, 77 856 + coarse prep; n_s=32/64/128 →
  4e-5/0.27/0.987); Path B T=32 not Trotterized (same ceiling, would need 512).
- UCCSD 16 sites: 1022 par., 335 360 CNOTs, 1−F 6.98e-3; pruned 255 par., 87 776
  CNOTs, 7.00e-3.  ADAPT 16 sites (pool 1022): 187 ops (68 008 CNOTs) for 1e-2, 9.45e-3
  after 200.  ADAPT 8 sites (pool 198): 66/101/149/194 ops (11 712/15 528/22 816/
  30 000 CNOTs) for 1e-2/3e-3/1e-3/3e-4; the 149-op state through the Section-11
  filter: 1e-3/3e-4/1e-4 at 2.3e4 (M=0)/6.9e4 (M=2)/1.5e5 (M=5) — does NOT beat the
  ramp (1.9e4/4.0e4/9.1e4), unlike the pairing model.
- filter at 16 sites: σ = 5/Δ_eff = 36.2, δt = 0.1, E_eff − E0 = +4.92e-3, floor
  3.7e-5 (`effective_step_state`), mean 289 steps = 2.13e5 CNOTs per cycle (9× the
  8-site cycle).  Matched comparison (20 time sets, M = 0..6, coarse prep of the
  hierarchy charged as the 8-site ramp T=16, 18 544 CNOTs), min CNOTs per accepted
  preparation at 1e-2 / 3e-3 / 1e-3 / 3e-4: ramp T=32 7.8e4 (0) / 3.0e5 (1) / 7.9e5
  (3) / 1.3e6 (5); Path B Λ=8 T=16 9.6e4 / 3.2e5 / 8.1e5 / 1.3e6 (tie within the
  10–35 % bootstrap errors); ADAPT 187 ops 6.8e4 / 5.4e5 / 1.0e6 / 1.4e6; pruned UCCSD
  8.8e4 / 5.6e5 / 1.1e6 / 1.3e6; Path A T=16 1.3e6 / 1.3e6 / 1.5e6 (M=1) / 2.0e6
  (M=3); ramp T=16 5.2e5 / 1.0e6 / 1.4e6 / –; prolonged 1.3e6 / 1.4e6 / – / –;
  determinant 1.5e6 / – / – / –.
- VERDICT at 32 qubits: the ramp (now T=32) is still cheapest, but Path B ties with it
  at every target, the variational routes cost 1.3–1.8×, and Path A costs 1.6–1.9× at
  the two tightest targets instead of 16×; the cost is set by the filter cycles each
  input still needs (2.1e5 each), i.e. by how close the input starts — the quantity
  the refinement controls.  The trend 8→16 is the one the proposal predicted.

Pipeline: `make hubbard` (notebook per-cell timeout 7200 s; `make check-hubbard`
expects 12 figures; `make paper-hubbard` builds main + SM).  Standalone seed of
Section 12: `scripts/hubbard16_checks.py` (its filter costs count the controlled
step twice — the notebook is the reference).  `tests/test_hubbard.py` has a test of
the Kronecker machinery / effective step state / sparse UCCSD / ADAPT.

Open: repository URL / Zenodo DOI in `\bibitem{Code}` and `\bibitem{SM}`; authors of
`Hubbard.ipynb` (`\bibitem{HubbardNotebook}`); co-author sign-off; the abstract is
~330 words (Quantum has no limit, but it could be trimmed).

## Status 2026-10-15: third draft (revision after the paper-I assessment), US spelling

Files: `eigenq_hubbard.tex` / `.pdf` (RevTeX, 15 pp., structured PRC abstract, same
author list as paper I; red `[...]` = open items), `figs/fig1-10.png` (extracted from
the notebook; the paper uses 1–7, 9, 10 — fig8 is superseded by fig10),
`eigenq_hubbard_v1_backup.tex` (2026-10-11 first draft), `eigenq_hubbard_v2_backup.tex`
(2026-10-14 second draft), `Hubbard.ipynb` (ED reference, unchanged), `prototype/`
(superseded sector-exact prototype; the notebook reproduces its numbers).

### Revision 2026-10-15 (the points of `paper/article_correctness_assessment.md` that apply)
Notebook Section 11 (cells 28–29, figs 9–10), all at U = 4, chain 4→8:
- Path A counted with the coarse term SPLIT into its three conditioned coarse groups
  [P(Dc−μ)P†, P Kc_even P†, P Kc_odd P†]: n_s = 64/128/256/512 → 0.997299/0.999818/
  0.999941/0.999943 (one block: 0.997456/0.999831/0.999941/0.999943); n_s = 128 and the
  296 976-CNOT cost unchanged.  The split state is the one in the comparison.
- Pruned UCCSD: 72 of 198 amplitudes > 1e-2, 98 > 1e-3; keeping those 98 (12 singles,
  86 doubles) and re-optimizing (`uccsd_hubbard(keep=, theta0=)`): E = −5.618755,
  1−F = 3.145e-3 (full: −5.618758, 3.143e-3), 15 120 CNOTs vs 30 864.
- Ordinary adiabatic ramp H(s) = K + sin²(πs/2) D from the free-fermion determinant
  (Givens network 112 CNOTs): p0 = 0.7804, Δ_eff = 0.588; exact slices p(4) = 0.96878,
  p(8) = 0.99686, p(12) = 0.99960, p(16) = 0.99989, p(24) = 0.99999 (as fast as Path A,
  Fig. 9 right).  Trotterized: T=8 n_s=64 → 0.996891 (9 328 CNOTs); T=16 n_s=128 →
  0.999659 (18 544 CNOTs).
- Response norm χ_u(s) of Path A: χ_u(0) = 0.375/0.527/0.715/0.874 (U = 0/2/4/8),
  maximum at s = 0 for every U; at U=8 it rises again to 0.45 at s=1 (Fig. 9 left).
- Fixed-step filter (quantized times, E_obj = E_eff): long-run state = the step
  eigenvector closest to the ground state (floor 5.8e-5), acceptance tends to the
  input weight on it (does not decay).
- Matched comparison (Table VI/VII, Fig. 10): eight inputs, M = 0..6, CNOTs per attempt
  and per accepted, bootstrap rse (≤ 18 % unrefined at M=6, < 10 % refined; mean 68
  steps = 2.4e4 CNOTs per cycle).  Minimum cost per accepted preparation at targets
  1e-3 / 3e-4 / 1e-4: ramp T=16 1.9e4 (M=0) / 4.0e4 (M=1) / 9.1e4 (M=3); ramp T=8
  5.6e4 / 1.1e5 / 1.5e5; UCCSD pruned 6.2e4 / 1.2e5 / –; UCCSD full 7.8e4 / 1.3e5 / –;
  Path B Λ=8 1.4e5 (M=4) / 1.6e5 (M=5) / –; Path A 3.0e5 / 3.0e5 / 3.4e5 (M=2);
  prolonged and determinant reach none.  Only the ramp T=16 and Path A reach the floor
  within six cycles (others a factor 1.5–80 above it at M=6).
- VERDICT CHANGED vs the second draft: the ordinary adiabatic ramp from the
  free-fermion determinant is the cheapest route at every target (×3 vs pruned UCCSD,
  ×4–7 vs Path B, ×16 vs Path A) because at this size the determinant has the same p0
  (0.78 vs 0.80) and the same Δ_eff as the prolonged coarse state.  The case for
  refinement is the regime where the coarse eigenstate beats the determinant (larger
  L, stronger U) — to be demonstrated (8→16 at circuit level).
- Text: structured abstract; error budget (√ε_low + a_ad + η)²; "calibrated against
  the exact answer"; Eq. (fid) = exact ensemble on resonance, 2^−M a lower bound,
  sampled values carry errors; UCCSD ordering stated; isolated-branch condition
  verified by Fig. 3; JRS bound vs T* = different criteria; paper I's inverse-gap
  statement no longer "re-read" (it was withdrawn there); fig:budget caption
  corrected ("refined inputs reach the floor; unrefined do not"); setup costs /
  amortization / unoptimized parameters qualifications; bibliography in citation
  order (HubbardNotebook before Jansen2007).

Pipeline: `src/pairinglib/hubbard.py` (library, `tests/test_hubbard.py`) →
`scripts/build_notebook_hubbard.py` → `notebooks/ResolutionRefHubbard.ipynb`
(31 cells, 10 figures, ~45 min clean run) → `scripts/check_notebook_hubbard.py`
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
  per application.  Trotterized refinement T=16, U=4: n_s = 128 (A), 256 (B Λ=8),
  512 (B Λ=16) second-order steps for |p−p_exact| ≤ 1e-3; first order not converged
  at 512.  Trotterized rodeo: fixed step diagonalized once, E_obj = Trotter-effective
  eigenphase, δt = 0.1 (floor 5.8e-5; δt = 0.25 gives 2.6e-3 and E_eff−E0 = 0.09).
  Matched comparison (second draft, Fig. 8): Path B Λ=8 3.7e4 CNOTs (p 0.9921) ~ UCCSD
  3.1e4 (p 0.9969); Path A 3.0e5 (p 0.9998, coarse term dominates); prolonged needs
  6 cycles/1.5e5 for 5e-3.  Superseded by the eight-input comparison above.

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
7. From the assessment: device-level calibration protocol for (T, μ_buf, n_s, δt,
   E_obj); ADAPT-VQE proper and a matched shot budget; the 8→16 chain at the circuit
   level (does the coarse eigenstate beat the free-fermion determinant at larger L?);
   joint search over (T, μ_buf, Λ, n_s, δt) for the minimum-cost curves.
