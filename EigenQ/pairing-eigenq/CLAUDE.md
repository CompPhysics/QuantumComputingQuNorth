# CLAUDE.md — project memory

This repository produces a reproducible Jupyter notebook and a LaTeX article on
**resolution refinement + the rodeo algorithm for eigenstate preparation**,
demonstrated on the constant-pairing Hamiltonian (the EIGEN-Q pipeline).

## Source of truth
- Physics lives in `src/pairinglib/` (a tested package). **Import it**; never
  re-define library functions inside notebooks or scripts.
- The notebook (`notebooks/ResolutionRefPairing.ipynb`) is a build artifact.
- There are TWO papers, split on 2026-10-12 (the pre-split manuscript is
  `paper/eigenq_pairing_v3_backup.tex`):
  * Paper I, `paper/eigenq_pairing.tex` — the pure (seniority-conserving)
    pairing model only; cites paper II as `\cite{PHcompanion}`.
  * Paper II, `paper_ph/eigenq_pairing_ph.tex` — the pairing model plus the
    seniority-breaking particle-hole term V_ph (model, refinement, bounds,
    rodeo, circuit budget at f=0.05/0.2/0.5); self-contained, with a short
    recap of the pipeline and of the f=0 reference numbers; cites paper I as
    `\cite{PairingCompanion}`.  Its figures are copies of notebook figures
    12-16 and 18 (`make figures-ph`); paper I uses `figs/fig_bound_pairing.png`
    (f=0-only redraw of fig18, `scripts/bound_figure_pairing.py`).
  Both are generated from the notebook's results and figures. Every number in
  either paper must trace to a notebook output.  Keep the two consistent:
  the f=0 row of paper II's tables is paper I's result.
- Paper III, `paper3/eigenq_hubbard.tex` (RevTeX, PRC): LATTICE refinement of
  the Hubbard model (chain 2->4->8->16, square 2x2->4x4), same authors and
  pipeline.  It has its OWN notebook, `notebooks/ResolutionRefHubbard.ipynb`,
  built from scratch by `scripts/build_notebook_hubbard.py` (nbformat, 31 cells,
  10 figures, ~45 min clean run) on the module `src/pairinglib/hubbard.py`
  (`import pairinglib.hubbard as hb`; tests in `tests/test_hubbard.py`).  Figures
  are `paper3/figs/fig1-10.png` (`make figures-hubbard`; the paper uses 1-7, 9, 10);
  `make check-hubbard` (10 figures) asserts the anchors (chain 4->8 U=4: E_low -4.976573, E_high -5.634427;
  16 sites: -7.209813; 4x4 (2,2): -10.4197).  The ED reference and its
  conventions are `paper3/Hubbard.ipynb` (unchanged); `paper3/prototype/` is the
  superseded sector-exact prototype whose numbers the notebook reproduces;
  `paper3/NOTES.md` is the running status; drafts are kept as
  `paper3/eigenq_hubbard_v1_backup.tex` (prototype numbers) and `_v2_backup.tex`
  (2026-10-14, before the assessment revision).  Paper III was revised on
  2026-10-15 against the paper-I assessment (notebook Section 11, US spelling,
  structured PRC abstract).

## Skills (in `.claude/skills/`) — read the matching one before working
- `pairing-model`  — Hamiltonian conventions, encoding, benchmark anchors, API.
- `vqe-circuits`   — UCCSD/HEA, prolongation, refinement, and the rodeo algorithm.
- `notebook-builder` — assemble/execute/validate the notebook (scripts included).
- `physics-paper`  — extract figures, write, and compile the article.

## Pipeline (also wired into the Makefile)
1. `make install`   — editable install of `pairinglib` + deps.
2. `make test`      — pytest: library API + benchmark anchors.
3. `make notebook`  — execute the notebook in place (generous timeout).
4. `make check`     — validate the executed notebook (no errors, figures, anchors).
5. `make figures`   — extract notebook PNGs into `paper/figs/`.
6. `make paper`     — compile `paper/eigenq_pairing.tex` to PDF (paper I).
7. `make paper-ph`  — copy the V_ph figures and compile paper II; `make papers` builds both.
8. `make all`       — the whole chain end to end (papers I and II).
9. `make hubbard`   — paper III: build+execute the Hubbard notebook, check it, extract
   figures, compile `paper3/eigenq_hubbard.tex` (needs `texlive-publishers` for
   revtex4-2; the individual targets are notebook-hubbard / check-hubbard /
   figures-hubbard / paper-hubbard).

## Conventions
- Python: import `pairinglib as pl`. Keep new physics in the package + a test.
- Notebooks: built with `notebook-builder` scripts (nbformat), never hand-edited JSON.
- LaTeX: `article` class + `authblk`; authors ALPHABETICAL by surname.
- Benchmark anchors (N=4, g=1): FCI(k=3)=0.794697, FCI(k=4)=0.635548; rodeo
  acceptance -> p. With the particle-hole term (k=4): f=0.05 -> 0.45058234,
  f=0.2 -> -0.18348455, f=0.5 -> -1.69173670 (book Chapter 4). CI and
  `make check` assert these.
- The seniority-breaking particle-hole term `V_ph` (strength `f`) is switched on
  with the `f=` keyword everywhere (`phmodel.py`); `f=0` is the pure pairing model.
- Notebook extensions are tagged sections (`scripts/extend_notebook_*.py` +
  `scripts/run_new_cells.py --mark`); see the notebook-builder skill.  Sections:
  EXT:TROTTER, EXT:PH, EXT:BASELINES (matched-observable baselines, sizes,
  observables; cell code is read from `scripts/baselines_rodeo.py` and
  `scripts/sizes_observables.py`, which also run standalone), EXT:BOUND
  (sufficient adiabatic-time / Trotter-step conditions of the nucleation paper,
  arXiv:2609.37527, evaluated on the 2->4 path; code from
  `scripts/adiabatic_bound.py`; paper Sec. IV.E, Table tab:bound, fig18),
  EXT:ASSESS (revision of 2026-10-14 after `paper/article_correctness_assessment.md`;
  code from `scripts/assessment_checks.py`, inserted by
  `scripts/extend_notebook_assess.py`, which ALSO patches three original cells:
  the fine energy scan is centred on the coarse-scan maximum, the rodeo
  convergence cell resets seed 7 per input and plots exact-ensemble curves with
  standard errors, the shift-scan panel shows both affine fits; figs 19-20),
  EXT:LARGE (2026-10-15: four pairs on eight levels, N=8, k=8, 16 qubits,
  seniority zero; code from `scripts/large_space_checks.py`, inserted by
  `scripts/extend_notebook_large.py`; figs 21-22; paper I Sec. V.F, Tables
  tab:large/tab:shift/tab:largebudget; `make check` now expects 22 figures).
  `scripts/verify_assessment.py` recomputes every number claimed in the assessment.
- Seniority-zero pair basis: `pairinglib.seniority` (`import pairinglib.seniority as sz`;
  `pair_basis`, `H_pair`, `pair_terms` = [D, h_pq] (the seniority-zero blocks of
  `model_terms`, verified), `prolong_pair`, `refine_pair`/`refine_pair_overlaps`,
  `ramp_terms`, `pair_double_ansatz` (the restricted P(k-P)-parameter UCCSD, 48
  CNOTs per rotation), `coarse_term_cnots` (Gray-code convention: 2->4 14,
  3->4 348, 4->8 (P=4) 254, 6->8 4872), `pair_encoding_cnots` (one qubit per
  level: step 112/240 at k=8)); tests in `tests/test_seniority.py`.  Use it for
  any pure-pairing study with k > 4 (70 states at k=8, P=4 instead of 12 870).
- k=8, N=8, g=1 anchors (paper I Sec. V.F): E0 = 8.889170, E_HF = 10, gap 2.592,
  p_HF = 0.690; 16-pair ansatz 1-F 3.1e-3 (768 CNOTs); chi_u(0) = P(k-P) g^2/(4 b^2)
  (= 11.1 at b=0.6); T*(1e-3) = 50/25/8/6/8/10 for mu_buf = 0.6/1/2/4/8/16 (n_s =
  -/160/40/40/40/80); HF ramp exact p(6) = 0.99993, min gap 2 at all g; step
  1360/1856 CNOTs; digital filter: dt=0.25 REJECTED (map eigenvalue 0.96, floor 1),
  dt=0.1 floor 4.6e-5 (eigenvalue 0.9933), dt=0.05 3.1e-6; cycle 1.2e5 CNOTs
  (sigma=4); g=2: sigma=1.5, dt=0.05 (floor 2.5e-5), cycle 9.0e4.  Min cost per
  accepted prep at 1e-3/3e-4/1e-4 (g=1): ramp T=6 1.1e5/2.2e5/4.8e5; refined
  mu_buf=4 T=6 1.3e5/2.4e5/5.0e5; 16-pair ansatz 2.4e5/5.1e5/7.4e5; determinant
  none (1.0e-2 at M=6).  g=2: ramp T=8 2.2e5/4.0e5/6.0e5; refined mu_buf=8 T=8
  2.6e5/4.4e5/6.4e5; ansatz 5.5e5/-/-.  VERDICT: at 16 qubits refinement+filter
  beats the variational route (x1.5-2) and the filtered determinant, but the
  HF ramp stays cheaper by the coarse-term overhead (4-19%) at g=1 and g=2 once
  the shift is scaled with the size (mu_buf=0.6 of k=4 is 6x too slow at k=8).
  Pair encoding: all adiabatic costs /12, cycles /8, same ranking.
- Compare like with like: an energy error is NOT an infidelity
  (<H>-E0 >= gap*(1-F)).  UCCSD k=4, g=1 infidelities: 1.78e-4 (f=0),
  1.01e-3 (f=0.2), 3.72e-3 (f=0.5).  Every comparison with the pipeline uses the
  acceptance-weighted ensemble infidelity of the accepted runs.
- Gate budgets must include the embedded coarse term A = P(H_low-mu)P^T of the
  refinement path (2->4: 4-qubit conditional phase, 14 CNOTs/application;
  3->4: doubly-conditioned k=3 step, ~350 CNOTs).  `cnot_counts` counts the
  high-space step only.
- Bounds vs. practice (k=2->4, N=4, g=1): T_bound(eps=1e-2)=5.1e3 vs T*=16;
  first-order step bound 3.4e4 at T* vs 160 (order 1) / 40 (order 2) in practice
  (different criteria: operator norm vs one overlap — say so).
- Rodeo statistics: quote the EXACT Gaussian-ensemble infidelity (Eq. fid) as the
  reference; sampled 40-set values carry relative standard errors 13%/37%/82% at
  M=2/6/10; seed convention = `default_rng(7)` RESET per input.  Detuning study
  (exact ensemble, sigma=3, M=6): d=0.05/0.1/0.2 -> acceptance 0.965/0.873/0.596.
- Scaling claims: the mu-scan points (0.3,50),(0.6,16),(1.0,8),(1.5,6) fit
  17.05/D-8.42 (SSR 28) AND 4.13/D^2+4.16 (SSR 0.22) — NO exponent is identified;
  chi_u(0)=g^2/b^2 varies x25 across the scan. Never write "confirms 1/Delta" or
  "incompatible with 1/Delta^2".  Gaps: shift 0.6 / sector 1.823 / seniority-zero 2.300.
- Comparators that must appear in any gate-budget statement: the four-pair
  restricted ansatz (192 CNOTs, 1-F=1.78e-4, identical to the 26-operator UCCSD
  whose other 22 amplitudes are exactly 0; k=3 coarse state: two rotations, 96
  CNOTs, exact) and the HF adiabatic ramp D -> H (T=4, 11 second-order steps,
  6512 CNOTs, p=0.9978, no shift/coarse term).  Minimum cost at targets
  2e-3/5e-4/2e-4: four-pair 192 (M=0); HF ramp 1.7e4/4.0e4/6.4e4; refined T=16
  4.9e4/7.4e4/9.4e4.  Digital filter (dt=0.25): averaged-map floor 1.2e-4,
  eigenvalue 0.9979 (acceptance decays 0.2%/cycle); dt=0.5: 1.5e-3; dt=1: 0.21.
  3->4 chain: T*=8, 22 steps with the coarse term SPLIT into coarse groups
  (`pl.trotter_path`), 28 520 CNOTs.
- Review history: `paper/revisions.txt` + `paper/response_to_revisions.md`
  (2026-10-02) and `paper/article_correctness_assessment.md` +
  `paper/response_to_assessment.md` (2026-10-14; previous tex kept as
  `paper/eigenq_pairing_v4_backup.tex`; v5 = before the reference/abstract/US
  pass; v6 = before the N=8, k=8 section of 2026-10-15).  Paper II has NOT been re-audited
  against the 2026-10-14 points (fixed-cycle comparisons, "never couples").
- Spelling: US English in all three papers (Trotterized, diagonalized, centered,
  neighbor, un-dimerization, analog, ...).
- Hubbard conventions (`hubbard.py`): t=1, OBC, block JW ordering (up orbitals
  then down), sorted uint64 masks, H = K_up (x) 1 + 1 (x) K_dn + U diag; block
  prolongation = bonding combination; Galerkin coarse model chain (t/2, U/2) - tN,
  square (t/2, U/4) - 2tN; Path A = projector path with mu = E_low + 2t; Path B =
  un-dimerization with extra intra-block hopping Lambda cos^2 (ceiling
  1-p ~ (t/Lambda)^2); rodeo width sigma = 5/Delta_eff with Delta_eff the lowest
  excitation carrying weight (symmetry-protected gap, NOT the sector gap);
  symmetry labels from `total_spin` (S^+ map), `site_permutation` (reflection) and
  `c4v_weights` (character projection).  Key lesson: the (2,2) 4x4 fine ground
  state is an (A2, S=1) TRIPLET, the coarse plaquette ground state a (B1, S=0)
  singlet -> zero overlap; refine the lowest coarse state in the target's sector.
- Hubbard circuit level: fine chain L=8 has 52 weight-<=2 strings; step CNOTs
  72 (order 1) / 144 (order 2) / 352 (controlled order 2); Path-A coarse term
  1088 CNOTs/application (Toffoli-ladder convention; 12272 with the Gray-code
  convention of paper I); T=16 needs n_s = 128 (A) / 256 (B, Lambda=8) / 512
  (B, Lambda=16) second-order steps.  Trotterized rodeo: fixed step dt=0.1
  diagonalized once (`trotter_evolution`), E_obj = Trotter-effective ground
  energy, floor 5.8e-5 (dt=0.25 gives 2.6e-3 — too coarse for Hubbard norms);
  quantized times => long-run state is the step eigenvector, acceptance does
  not decay.  Path A is COUNTED with the coarse term split into its three
  conditioned coarse groups (n_s=128 unchanged, 296 976 CNOTs).  Baselines
  that must appear in every paper-III budget statement (U=4, chain 4->8):
  the ordinary U-ramp K -> K+D from the free-fermion determinant (Givens 112
  CNOTs; exact p(8)=0.99686, p(16)=0.99989; Trotterized T=8 n_s=64 9 328
  CNOTs p=0.99689, T=16 n_s=128 18 544 CNOTs p=0.99966) and the PRUNED UCCSD
  (98 of 198 amplitudes > 1e-3, `uccsd_hubbard(keep=, theta0=)`, 15 120 CNOTs,
  same 1-F 3.1e-3).  Minimum cost per accepted preparation at 1e-3/3e-4/1e-4:
  ramp T=16 1.9e4/4.0e4/9.1e4; pruned UCCSD 6.2e4/1.2e5/-; Path B (Lambda=8)
  1.4e5/1.6e5/-; Path A 3.0e5/3.0e5/3.4e5.  VERDICT: the U-ramp is the
  cheapest at every target (x3 pruned UCCSD, x4-7 Path B, x16 Path A) because
  the determinant has the same p0 and Delta_eff as the prolonged state at this
  size; refinement's case is larger L / stronger U, not demonstrated.  Response
  norm chi_u(0) = 0.375/0.527/0.715/0.874 (U=0/2/4/8), max at s=0.

## Definition of done
A change is done only when `make test` and `make check` pass and (if a paper
changed) `make paper` / `make paper-ph` / `make paper-hubbard` builds without
unresolved references; a change to `hubbard.py` or the Hubbard notebook also
needs `make check-hubbard` to pass.
