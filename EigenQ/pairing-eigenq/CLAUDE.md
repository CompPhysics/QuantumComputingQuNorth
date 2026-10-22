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
  * Paper I, `paper/eigenq_pairing.tex` (main text) + `paper/eigenq_pairing_supp.tex`
    (Supplemental Material, cross-referenced with `xr` against the main .aux;
    `make paper` builds both) — the pure (seniority-conserving) pairing model.
    Restructured 2026-10-16 (pre-split tex: `eigenq_pairing_v7_backup.tex`;
    `scripts/split_supplement.py` is the one-off restructuring script, kept for
    the record) around (a) the reference implementation, (b) the matched-accuracy
    benchmark at 8 and 16 qubits, (c) three transferable results (shift rule
    chi_u(0) = P(k-P) g^2/(4 b^2), averaged post-selection map of the digital
    filter, bounds vs practice).  Title: "Resolution refinement followed by the
    rodeo algorithm: a gate-level benchmark of eigenstate-preparation routes on
    the pairing Hamiltonian".  Paper I NO LONGER cites papers II or III (user
    decision 2026-10-16); the SM holds figs S1-S8 / tables S1-S5 (fig4, fig3,
    fig6, fig10, fig19, fig_bound, fig7, fig21; tab:bound, tab:ops, tab:sizes,
    tab:obs, tab:large).  Submission target: PRC.  Open before submission:
    repository URL + Zenodo DOI (bibitem `Code`, CITATION.cff), a LICENSE file,
    co-author sign-off.
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
- Paper III, `paper3/eigenq_hubbard.tex` (main text, `quantumarticle` class for
  Quantum, two-column, `booktabs` tables) + `paper3/eigenq_hubbard_supp.tex`
  (Supplemental Material, Secs. S1-S8, `xr` against the main .aux): LATTICE
  refinement of the Hubbard model (chain 2->4->8->16, square 2x2->4x4), written
  along the lines of paper I (reference implementation / matched benchmark at
  8 AND 16 sites / transferable results; ADAPT-VQE comparator; cites paper I,
  no "companion" references).  It has its OWN notebook,
  `notebooks/ResolutionRefHubbard.ipynb`, built from scratch by
  `scripts/build_notebook_hubbard.py` (nbformat, 36 cells, 14 figures; Sections
  1-11 ~70 min, Section 12 = the 8->16 chain at the circuit level + ADAPT,
  ~1 h, Section 13 = the pre-submission revision of 2026-10-20, ~3.5 h) on the
  module `src/pairinglib/hubbard.py` (`import pairinglib.hubbard as hb`; tests
  in `tests/test_hubbard.py`).  `scripts/run_notebook_checkpoint.py` executes
  it in-process with a dill checkpoint (`--checkpoint-after 15` = after
  Section 12) so that Section 13 can be fixed and resumed without the earlier
  hours; the builder takes an optional output path.  Figures are
  `paper3/figs/fig1-14.png` (`make figures-hubbard`; main text uses 1, 2, 3, 5,
  11, 13, 14; the SM uses 4, 6, 7, 9; figs 8, 10, 12 are superseded by 13-14);
  `make check-hubbard` (14 figures) asserts the anchors (chain 4->8 U=4: E_low
  -4.976573, E_high -5.634427; 16 sites: -7.209813; 4x4 (2,2): -10.4197).  The
  ED reference and its conventions are `paper3/Hubbard.ipynb` (unchanged);
  `paper3/prototype/` is the superseded sector-exact prototype; `paper3/NOTES.md`
  is the running status; drafts are kept as `paper3/eigenq_hubbard_v1_backup.tex`
  (prototype numbers), `_v2_backup.tex` (2026-10-14), `_v3_backup.tex`
  (2026-10-15 RevTeX/PRC third draft, structured abstract).  The standalone
  16-site study that seeded notebook Section 12 is `scripts/hubbard16_checks.py`
  (superseded; it counted 2x the controlled step per filter cycle and had the
  Galerkin-constant bug below; the notebook is the reference).
  PRE-SUBMISSION REVISION 2026-10-20 (`paper3/hubbard_pre_submission_critical_assessment.md`,
  response `paper3/response_to_pre_submission_assessment.md`, previous tex
  `eigenq_hubbard_v4_backup.tex` / `eigenq_hubbard_supp_v4_backup.tex`), all
  points agreed and confirmed numerically: (1) `galerkin_params(kind, t, U, N)`
  takes the PARTICLE number — Section 12 passed the coarse site count 8 instead
  of 4, i.e. an extra -4t cos^2 P P^T in the 16-site Path A (p(8,128) 0.99078 ->
  0.99237, p(16,128) 0.99835 -> 0.99930; the "plateau at T=8" was this bug, not
  the split coarse groups); (2) `uccsd_cnots` now uses the parity-support
  weight (q-p+1)+(s-r+1) of a double (templates 1720->1464, 30864->23312,
  335360->235904; pruned 8 sites 11408, 16 sites 61584; ADAPT 8 sites
  66/101/149/194 ops = 8832/11672/17168/22608, 16 sites 187 ops = 48008);
  (3) `conditioned_term_cnots` has 'flag' (shared flag, base + 2 n_rot + 6(c-1)
  + 2: 124 for 4->8, 268 for 8->16; 52 rotations, not 48) and 'frame' (32/64);
  Path A is costed flag+frame = 156 / 332 per application (was 1088 / 4856 with
  a ladder per rotation; `toffoli` now 1132/4948 incl. the scalar); (4)
  `chain_step_cnots(controlled=True)` absorbs the scalar single-Z sum
  (288 / 608 instead of 352 / 736; `absorb_scalar=False` for the old counts);
  (5) `exact_digital_ensemble` (rho_ij = c_i c_j^* Q_ij^M) replaces the sampled
  schedule averages at 8 sites (T=8 ramp n_s=64: 1.05045e-4 at M=6, NOT 8.9e-5);
  `effective_step_state(return_residual=True)` (16 sites: residual 1.9e-15);
  (6) the coarse state of the hierarchy is PREPARED and PROPAGATED (4-site ramp
  for 4->8, the cheapest 8-site unfiltered route for 8->16; two options each,
  'c' and 'a'), and the adiabatic routes are searched jointly over
  (T, n_s, M) / (Lambda, T, n_s, M) with M=0 allowed (notebook Section 13,
  figs 13-14, paper III Tables 2-3).  The old verdict ("ramp cheapest at every
  target, x16 Path A") is WITHDRAWN; the new numbers are in the "Hubbard
  8->16" entry below and in `paper3/response_to_pre_submission_assessment.md`.

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
6. `make paper`     — compile `paper/eigenq_pairing.tex` to PDF (paper I), then the SM
   against `paper/eigenq_pairing_xr.aux` (the main .aux without FirstPage/LastPage/
   LastBibItem/\bibcite, so that xr raises no multiply-defined warnings; the SM
   bibliography is self-contained).
7. `make paper-ph`  — copy the V_ph figures and compile paper II; `make papers` builds both.
8. `make all`       — the whole chain end to end (papers I and II).
9. `make hubbard`   — paper III: build+execute the Hubbard notebook (per-cell
   timeout 7200 s), check it (12 figures), extract figures, compile the main text
   and the SM (`quantumarticle` class, in texlive-publishers; the individual
   targets are notebook-hubbard / check-hubbard / figures-hubbard / paper-hubbard).

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
  tab:shift/tab:largebudget, SM tab:large/fig:large), EXT:ADAPT (2026-10-16:
  ADAPT-VQE comparator in the pair basis, `scripts/adapt_checks.py` via
  `scripts/extend_notebook_adapt.py`, fig 23; paper I Sec. V.D, Table tab:adapt;
  `make check` expected 23 figures), EXT:FRONTIER (2026-10-20, after
  `paper/pairing_pre_submission_critical_assessment.md`, response in
  `paper/response_to_pre_submission_assessment.md`, previous tex
  `eigenq_pairing_v8_backup.tex`: converged averaged map + exact-channel joint
  operating-point search, `scripts/frontier_checks.py` via
  `scripts/extend_notebook_frontier.py`, fig 24 = paper I fig:largebudget,
  Table tab:frontier), EXT:REFINE (2026-10-21, after `paper/morefeedback.md`,
  response in `paper/response_to_morefeedback.md`, previous tex
  `eigenq_pairing_v9_backup.tex` / `eigenq_pairing_supp_v9_backup.tex`:
  convergence of the operating-point search — fine (T, mu_buf) grids, integer
  n_s by bracketing + bisection at every target and M, three search levels
  L0/L1/L2 and the flatness count; `scripts/frontier_refine.py` via
  `scripts/extend_notebook_refine.py`, fig 25 = paper I fig:largebudget (the
  cost frontier), fig 24 moved to the SM (fig:frontiergrid), Table tab:frontier
  = level L2, SM Table tab:converge; `make check` now expects 25 figures).  The
  averaged map of the digital filter MUST be integrated piecewise between the
  step-count boundaries |t| = n dt (`pl.half_gaussian_nodes` /
  `pl.AveragedFilterMap` / `pl.averaged_map_floor`); never use a global
  Gauss-Hermite rule for it.  `AveragedFilterMap.track` applies the map through
  its cached d^2 x d^2 superoperator (`superoperator()`, 20x faster at d=70,
  identical numbers; node loop for d > 128), and `pl.trotter_path` applies each
  factor as two matrix-vector products (`_apply_factor`, 6x faster).  The cycle
  times of the DIGITAL filter are t = |z|, z ~ N(0, sigma^2) (half-normal,
  positive times) in the sampled schedules AND in the averaged map; the ideal
  filter uses signed draws (equivalent there).  A signed-time digital channel
  would differ (k=4 floors 1.1385e-4 vs 1.1399e-4 at dt=0.25, 0.129 vs 0.154
  at dt=1) — say which convention, never mix them.
  `scripts/verify_assessment.py` recomputes every number claimed in the assessment.
- Seniority-zero pair basis: `pairinglib.seniority` (`import pairinglib.seniority as sz`;
  `pair_basis`, `H_pair`, `pair_terms` = [D, h_pq] (the seniority-zero blocks of
  `model_terms`, verified), `prolong_pair`, `refine_pair`/`refine_pair_overlaps`,
  `ramp_terms`, `pair_double_ansatz` (the restricted P(k-P)-parameter UCCSD, 48
  CNOTs per rotation), `coarse_term_cnots` (Gray-code convention incl. the
  conditioned scalar phase: 2->4 14, 3->4 350, 4->8 (P=4) 254, 6->8 4886;
  convention="flag" = shared flag with c-1 clean work qubits: 224 / 44 / 1028), `pair_encoding_cnots` (one qubit per
  level: step 112/240 at k=8), `adapt_vqe(k, P, g, pool='pairs'|'generalized')`
  (gradient selection with repeats, warm-started L-BFGS re-optimization,
  history[n] holds the state after n operators)); tests in
  `tests/test_seniority.py`.  Use it for any pure-pairing study with k > 4
  (70 states at k=8, P=4 instead of 12 870).
- ADAPT-VQE anchors (generalized pool, 48 CNOTs/operator): (4,4) g=1: 4 ops =
  the four-pair ansatz, 5 ops (240 CNOTs) EXACT; (6,4): 9 ops (432) < 1e-5;
  (8,8) g=1: 19 (912) 1.5e-4, 23 (1104) 1.7e-5, 27 (1296) 4e-6; g=2: 25/29/31
  ops (1200/1392/1488) for 1e-3/3e-4/1e-4.  The 23-op state is BELOW the
  dt=0.1 filter floor (4.8e-5): the filter only degrades it; the g=2 31-op state
  (8.9e-5) is ABOVE the floor 2.7e-5 and improves x3 in six cycles.  VERDICT
  (final, 2026-10-20): the adaptive ansatz undercuts every other route by 2
  orders of magnitude at both sizes; with the adiabatic controls searched
  jointly (Table tab:frontier) ramp and refinement TIE within the coarse-term
  overhead at every size/coupling and need NO filter cycle at the targets.
- k=8, N=8, g=1 anchors (paper I Sec. V.F): E0 = 8.889170, E_HF = 10, gap 2.592,
  p_HF = 0.690; 16-pair ansatz 1-F 3.1e-3 (768 CNOTs); chi_u(0) = P(k-P) g^2/(4 b^2)
  (= 11.1 at b=0.6); T*(1e-3) = 50/25/8/6/8/10 for mu_buf = 0.6/1/2/4/8/16 (n_s =
  -/160/40/40/40/80); HF ramp exact p(6) = 0.99993, min gap 2 at all g; step
  1360/1856 CNOTs; digital filter (CONVERGED boundary-split quadrature,
  `pl.averaged_map_floor`, 2026-10-20): dt=0.25 REJECTED (eigenvalue 0.794,
  acceptance -21%/cycle, floor 2.3e-3 — NOT "zero overlap", that was a
  Gauss-Hermite artifact), dt=0.1 floor 4.8e-5 (eigenvalue 0.99324), dt=0.05
  3.1e-6; cycle 1.2e5 CNOTs (sigma=4); g=2: sigma=1.5, dt=0.1 4.0e-4, dt=0.05
  floor 2.7e-5 (0.99887), cycle 9.1e4.  Fixed-input (sampled, Table V) costs:
  ramp T=6 1.1e5/2.2e5/4.8e5; refined mu_buf=4 T=6 1.3e5/2.4e5/5.0e5; ansatz
  2.4e5/5.1e5/7.4e5.  JOINT SEARCH with the exact channel (Table tab:frontier,
  `scripts/frontier_checks.py`, fig24) at 1e-3/3e-4/1e-4, g=1: ramp T=6
  n_s=40/60/80 M=0 1.1e5/1.6e5/2.2e5; refined mu_buf=4 T=6 same n_s M=0
  1.3e5/1.9e5/2.6e5 (+19% = coarse term 254/1614); ansatz+filter
  2.4e5/4.9e5/7.4e5; ADAPT 912/912/1104.  g=2: ramp T=8 n_s=80 M=0/2/4
  2.2e5/4.0e5/5.8e5; refined mu_buf=8 (T=8,80,M=0) 2.6e5, mu_buf=16
  (T=8,120,M=0) 3.9e5, (M=2) 5.7e5; ADAPT 1200/1392/1488.  k=4 (2e-3/5e-4/2e-4):
  ramp (4,16,0) 9.5e3, (6,32,0) 1.9e4, (8,43,0) 2.5e4; refined mu_buf=3 T=6
  n_s=16/22/32 M=0 9.9e3/1.4e4/2.0e4.  VERDICT: the cheapest adiabatic operating
  points use NO filter cycle (a cycle = doubling n_s and buys less); ramp and
  refinement tie within the coarse-term overhead; mu_buf=0.6 of k=4 is 6x too
  slow at k=8 (matched thresholds: T*(1e-3) 40 -> 50, T*(1e-2) 16 -> 30).
  Pair encoding: all adiabatic costs /12, cycles /8, same ranking.
- Compare like with like: an energy error is NOT an infidelity
  (<H>-E0 >= gap*(1-F)).  UCCSD k=4, g=1 infidelities: 1.78e-4 (f=0),
  1.01e-3 (f=0.2), 3.72e-3 (f=0.5).  Every comparison with the pipeline uses the
  acceptance-weighted ensemble infidelity of the accepted runs.
- Gate budgets must include the embedded coarse term A = P(H_low-mu)P^T of the
  refinement path (2->4: 4-qubit conditional phase, 14 CNOTs/application;
  3->4: doubly-conditioned k=3 step incl. the conditioned scalar phase, 350
  CNOTs; ancilla-free reference compilation — quote the shared-flag alternative
  and its extra work qubits when discussing the penalty).  `cnot_counts` counts
  the high-space step only.
- Bounds vs. practice (k=2->4, N=4, g=1): T_bound(eps=1e-2)=5.1e3 vs T*=16;
  first-order step bound 3.4e4 at T* vs 160 (order 1) / 40 (order 2) in practice
  (different criteria: operator norm vs one overlap — say so).
- Rodeo statistics: quote the EXACT Gaussian-ensemble infidelity (Eq. fid) as the
  reference; sampled 40-set values carry relative standard errors 13%/37%/82% at
  M=2/6/10; seed convention = `default_rng(7)` RESET per input.  Detuning study
  (exact ensemble, sigma=3, M=6): d=0.05/0.1/0.2 -> acceptance 0.965/0.873/0.596.
- Scaling claims: the mu-scan points (0.3,50),(0.6,16),(1.0,8),(1.5,6) fit
  17.05/D-8.42 (SSR 28) AND 4.13/D^2+4.16 (SSR 0.22) — NO exponent is identified;
  chi_u(0)=g^2/b^2 (affine-path response) varies x25 across the scan, the
  scheduled chi_s = (pi/2 sin pi s)^2 chi_u only x3.4 (max 1.40/0.76/0.52/0.41).
  General identity: chi_u(0) = ||Q H_high P psi_low||^2 / b^2 (Eq. chi0);
  P(k-P)g^2/(4b^2) is its value for the filled pairing seed only. Never write "confirms 1/Delta" or
  "incompatible with 1/Delta^2".  Gaps: shift 0.6 / sector 1.823 / seniority-zero 2.300.
- Comparators that must appear in any gate-budget statement: the four-pair
  restricted ansatz (192 CNOTs, 1-F=1.78e-4, identical to the 26-operator UCCSD
  whose other 22 amplitudes are exactly 0; k=3 coarse state: two rotations, 96
  CNOTs, exact) and the HF adiabatic ramp D -> H (T=4, 11 second-order steps,
  6512 CNOTs, p=0.9978, no shift/coarse term).  Minimum cost at targets
  2e-3/5e-4/2e-4: four-pair 192 (M=0); HF ramp 1.7e4/4.0e4/6.4e4; refined T=16
  4.9e4/7.4e4/9.4e4 (ratios 2.9/1.9/1.5, NOT "2.5-3 at every target").  Digital
  filter (dt=0.25): averaged-map floor 1.14e-4, eigenvalue 0.99782 (acceptance
  decays 0.2%/cycle); dt=0.5: 1.9e-3 (0.9651); dt=1: 0.15 (0.6297) — boundary-split
  quadrature; the old Gauss-Hermite values 1.2e-4/1.5e-3/0.21 were unconverged.
  3->4 chain: T*=8, 22 steps with the coarse term SPLIT into coarse groups
  (`pl.trotter_path`), 28 520 CNOTs.
- Review history: `paper/revisions.txt` + `paper/response_to_revisions.md`
  (2026-10-02) and `paper/article_correctness_assessment.md` +
  `paper/response_to_assessment.md` (2026-10-14; previous tex kept as
  `paper/eigenq_pairing_v4_backup.tex`; v5 = before the reference/abstract/US
  pass; v6 = before the N=8, k=8 section of 2026-10-15; v7 = before the
  2026-10-16 restructuring into main text + Supplemental Material).  Paper II has NOT been re-audited
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
- Hubbard 8->16 at the circuit level (notebook Sections 12-13, paper III Sec. 6,
  Table 3, Figs. 11 and 14; U=4, (2,2), 14 400 states, 32 qubits): prolonged p0 =
  0.8692 vs determinant 0.7162, Delta_eff = 0.1379 for both, sector gap 0.051.
  Exact slices: ramp p(32/48/64) = 0.99591/0.99887/0.99919; Path A mu_buf=2
  p(8/16/32) = 0.99237/0.99931/0.99995 (larger shifts worse); Path B ceiling
  ~5e-3 (Lambda=8) / 4e-3 (Lambda=16).  Gates: step 304 (controlled 608), coarse
  term 332 (shared flag + frame; ladder 4948), Givens 480, prolongation 32;
  Trotter: ramp T=16/32 n_s=128/256, Path A T=8/16 n_s=64/128 (0.99228/0.99930,
  converged — no plateau once the Galerkin constant is -t N), Path B Lambda=8
  T=16 n_s=256.  UCCSD 1022 par. 235 904 CNOTs 6.98e-3 (pruned 255 par.,
  61 584); ADAPT 187 ops (48 008 CNOTs) for 1e-2, 9.45e-3 after 200.  Filter:
  sigma 36.2, dt 0.1, E_eff-E0 +4.92e-3, floor 3.7e-5, residual 1.9e-15,
  1.75e5 CNOTs/cycle (289 controlled steps).  Section-13 search (40 schedules,
  coarse state 'c' = 8-site ramp T=16 n_s=192 27 760 CNOTs 1.8e-4, 'a' = T=24
  n_s=256 36 976 5.6e-5), min CNOTs/accepted at 1e-2/3e-3/1e-3/3e-4: ramp
  7.8e4 (32,256,0)/7.8e4 (48,256,0)/1.56e5 (48,512,0)/1.56e5 (64,512,0);
  Path A 9.0e4/1.52e5/1.52e5/2.76e5 (all M=0, coarse c); ADAPT 4.8e4 (0)/
  3.9e5 (2)/7.9e5 (4)/1.1e6 (6); pruned 6.2e4/4.1e5/5.9e5/9.8e5; Path B
  1.06e5/2.7e5/5.3e5/9.2e5; prolonged 8.1e5/1.1e6/-/-; determinant 1.2e6/-.
  Bootstrap of the selection (400 resamples): ramp/Path A 100 % stable.
  Eight sites (exact ensemble, cycle 1.95e4) at 1e-3/3e-4/1e-4: ramp 1.4e4
  (12,96,0)/2.8e4 (16,192,0)/3.7e4 (24,256,0); ADAPT 1.7e4 (149)/2.3e4 (194)/
  7.9e4 (177, M=3); Path A 3.7e4/6.7e4/1.1e5; pruned 5.1e4/9.0e4/-; Path B
  1.2e5/1.6e5/2.7e5.  VERDICT: unfiltered ramp cheapest from 3e-3 down at both
  sizes; Path A (shared flag) x2.4-3.1 at 8 sites, tied at 1e-3 at 16 sites and
  x1.15-1.9 elsewhere; ADAPT cheapest only where it needs no cycle; Path B
  x6-9 / x1.4-5.9; no adiabatic route needs a filter cycle at these targets.

## Definition of done
A change is done only when `make test` and `make check` pass and (if a paper
changed) `make paper` / `make paper-ph` / `make paper-hubbard` builds without
unresolved references; a change to `hubbard.py` or the Hubbard notebook also
needs `make check-hubbard` to pass.
