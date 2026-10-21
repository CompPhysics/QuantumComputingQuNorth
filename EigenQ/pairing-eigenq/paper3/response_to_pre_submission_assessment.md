# Response to the pre-submission critical assessment of paper III (19 October 2026)

Assessment: `paper3/hubbard_pre_submission_critical_assessment.md`.  Every numerical claim was
reproduced with the repository code before any change (all matched to the digits quoted), the
code was corrected (`src/pairinglib/hubbard.py`, tests in `tests/test_hubbard.py`), the notebook
was rebuilt from scratch (`scripts/build_notebook_hubbard.py`, new Section 13) and the manuscript
and Supplemental Material were rewritten where the numbers changed.  Previous versions:
`paper3/eigenq_hubbard_v4_backup.tex`, `paper3/eigenq_hubbard_supp_v4_backup.tex`.

## Verdict per point

| # | Point | Verdict | What changed |
|---|-------|---------|--------------|
| 1 | Wrong constant (−8t instead of −4t) in the 32-qubit Path-A path | **Agree, confirmed** (p(8/t,128) 0.990785 → 0.992368; p(16/t,128) 0.998347 → 0.999298; exact slices 0.992368 at 192 slices). The S7 explanation of the "plateau" was wrong. | `galerkin_params(..., N_particles)` in notebook Sec. 12 + assert that the embedded groups sum to P(H_low−μ)P†; new test with N_coarse ≠ N_particles; `scripts/hubbard16_checks.py` fixed; all 32-qubit Path-A results, tables, figures and prose regenerated; the S7 "different limiting path" paragraph removed. |
| 2 | Conditioned coarse cost avoidable (shared flag); 48 vs 52 rotations; frame change | **Agree.** `conditioned_term_cnots` now returns `flag` (base + 2 n_rot + 6(c−1) + 2 for the conditioned scalar: 124 / 268) and `frame` (32 / 64); 52 rotations. | Path A costed with flag + frame throughout (156 / 332 CNOTs per application instead of 1 088 / 4 856); clean-work-qubit requirement (c−1 = 7 / 15) stated in Secs. 5.1, 7 and the hardware remarks; "never affordable" withdrawn. |
| 3 | JW double-excitation weight = parity support, not full span | **Agree, confirmed** (1 720 → 1 464, 30 864 → 23 312, 335 360 → 235 904). | `uccsd_cnots` corrected; every pruned/ADAPT operator list repriced; the orbital-basis order of the excitation blocks and the single final Givens network stated in Sec. 3.5. |
| 4 | Longer unfiltered ramps change the comparison; fixed inputs are not a frontier | **Agree, confirmed** (8 sites: T=24, n_s=192 → 1.53e-4 at 27 760; 16 sites: T=64, n_s=512 → 2.49e-4 at 156 128). | New Sec. 13 of the notebook: joint search over (T, n_s, M) for the ramp, (Λ, T, n_s, M) for Path B and (T, n_s, M) for Path A, M = 0 allowed, same grid density, with every ADAPT prefix admitted; Tables 2 and 3 replaced by the minima per route; Figs. 10 and 12 replaced by the frontier figures; abstract, Secs. 5.3, 6, 7 rewritten from the new tables. |
| 5 | Coarse preparation charged but its error not propagated; 18 544-CNOT ramp is 3.4e-4, not < 3e-4 | **Agree.** | The actual coarse state is prepared and propagated at both levels: a four-site ramp (< 1e-5) for 4→8, and the cheapest eight-site unfiltered route below 1e-4 for 8→16 (chosen from the eight-site frontier); the exact-coarse-oracle comparison of Sec. 12 is kept in the notebook, labelled as such. |
| 6 | Sound parts | Preserved. | — |
| 7 | χ_u is a metric, not leakage; "flat in size"/order-t gap overstated; transport vs input-support gap; continuum sentence; SM ‖B−A‖; projector distance; Path-B endpoints; plaquette Givens count | **Agree on all.** The UCCSD support below Δ_eff was checked: 3.1e-7 (full), 2.5e-7 (pruned). | Secs. 3.1–3.3, 4.2, 7 and SM S1 reworded accordingly (accessible gap 2t[cos π/(L+1) − cos 3π/(L+1)] closes with L; favorable observation, not a scaling law). |
| 8 | Statistics: use the exact digital ensemble at 8 sites (T=8 ramp: 1.05e-4, not 8.9e-5); residual for the 16-site eigenvector; bootstrap the operating-point selection | **Agree, confirmed** (1.05045e-4 reproduced). | `exact_digital_ensemble` (ρ_ij = c_i c_j* Q_ij^M, coherences kept) replaces the sampled averages at eight sites; at sixteen sites 40 shared schedules, bootstrap errors per input AND a bootstrap of the cost-minimizing (input, M) per route; `effective_step_state(return_residual=True)` reports ‖Uψ − e^{−iEδt}ψ‖; the superseded diagnostic script is marked. |
| 9 | Scope: amortized logical CNOT metric, ancilla counts, optimizer-free vs ED-calibrated, HVA comparator, scalar single-Z (352 → 288, 736 → 608), dilute family, noiseless benchmark | **Agree.** | Scalar absorbed in `chain_step_cnots(controlled=True)` (default); the metric, the qubit footprint (2L_f + 1 + c − 1 for Path A), the dilute fixed-N scope and the noiseless character are stated; the HVA is named as the untested natural comparator and the variational claims are restricted to the pools tested. |
| 10 | "Two orders of magnitude" (actually ×2.2); "ramp cheapest at every target" (ADAPT cheaper at 1e-2); 48 vs 52; "every U>0"; availability placeholders | **Agree.** | All fixed; the repository URL / DOI and the `Hubbard.ipynb` authorship remain to be supplied by the authors. |

## What the rebuilt comparison says

All numbers from notebook Section 13 (`notebooks/ResolutionRefHubbard.ipynb`, 36 cells, 14 figures,
`scripts/check_notebook_hubbard.py ... 14` passes).  Minimum CNOTs per accepted preparation per
route, with the operating point in parentheses:

**Eight sites (16 qubits), exact quantized-time ensemble, cycle 1.95e4 CNOTs** — targets 1e-3 / 3e-4 / 1e-4

| route | 1e-3 | 3e-4 | 1e-4 |
|---|---|---|---|
| U-ramp from the determinant | 1.4e4 (T=12, 96, M=0) | 2.8e4 (16, 192, 0) | 3.7e4 (24, 256, 0) |
| ADAPT-VQE | 1.7e4 (149 op., 0) | 2.3e4 (194, 0) | 7.9e4 (177, 3) |
| Path A (shared flag, incl. 8 216-CNOT coarse prep.) | 3.7e4 (12, 64, 0) | 6.7e4 (12, 128, 0) | 1.1e5 (16, 128, 'a', 2) |
| UCCSD pruned (98) + filter | 5.1e4 (2) | 9.0e4 (4) | — |
| UCCSD full (198) + filter | 6.3e4 (2) | 1.0e5 (4) | — |
| Path B | 1.2e5 (16, 24, 512, 2) | 1.6e5 (16, 24, 512, 4) | 2.7e5 (16, 32, 1024, 6) |
| prolonged / determinant + filter | — | — | — |

**Sixteen sites (32 qubits), 40 shared schedules, cycle 1.75e5 CNOTs** — targets 1e-2 / 3e-3 / 1e-3 / 3e-4
(coarse state 'c' = eight-site ramp T=16, n_s=192, 1.78e-4, 27 760 CNOTs; 'a' = T=24, n_s=256, 5.6e-5, 36 976)

| route | 1e-2 | 3e-3 | 1e-3 | 3e-4 |
|---|---|---|---|---|
| U-ramp from the determinant | 7.8e4 (32, 256, 0) | 7.8e4 (48, 256, 0) | 1.56e5 (48, 512, 0) | 1.56e5 (64, 512, 0) |
| Path A (shared flag, incl. coarse prep.) | 9.0e4 (8, 64, c, 0) | 1.52e5 (16, 128, c, 0) | 1.52e5 (16, 128, c, 0) | 2.76e5 (32, 256, c, 0) |
| ADAPT-VQE | 4.8e4 (187, 0) | 3.9e5 (187, 2) | 7.9e5 (187, 4) | 1.1e6 (187, 6) |
| UCCSD pruned (255) + filter | 6.2e4 (0) | 4.1e5 (2) | 5.9e5 (3) | 9.8e5 (5) |
| Path B | 1.06e5 (8, 16, 256, c, 0) | 2.7e5 (8, 16, 256, c, 1) | 5.3e5 (16, 32, 512, c, 2) | 9.2e5 (16, 32, 512, c, 4) |
| prolonged + filter | 8.1e5 (4) | 1.1e6 (6) | — | — |
| determinant + filter | 1.2e6 (6) | — | — | — |

Bootstrap of the operating-point selection (400 resamples of the 40 schedules): the ramp and
Path-A selections are stable in 100 % of resamples at every target, Path B 97–100 %, ADAPT 82–100 %,
pruned UCCSD 72–100 %, prolonged 61–100 %; the 5–95 % ranges of the minima are in SM Table S3.

**Conclusions that replaced the previous ones**

- The cheapest route from 3e-3 down at both sizes is an *unfiltered* adiabatic ramp from the
  free-fermion determinant (previously: "ramp cheapest at every target" — now false at 1e-3 and 3e-4
  at eight sites and at 1e-2 at sixteen, where ADAPT is cheaper).
- Path A with the shared-flag coarse term is the cheapest route that starts from the coarse
  eigenstate: ×2.4–3.1 the ramp at eight sites (previously ×16), tied at 1e-3 at sixteen sites
  (1.52e5 vs 1.56e5) and ×1.15–1.9 elsewhere (previously "×1.6–1.9 at the two tightest targets",
  reached only with the −8t bug and filter cycles).  "Never affordable" is withdrawn for the
  shared-flag compilation; it stands for the per-rotation ladder.
- Path B is no longer "tied with the ramp": ×6–9 at eight sites, ×1.4–5.9 at sixteen, because its
  ceiling forces cycles.
- Filter cycles are not needed by any adiabatic route at the targets tested; a cycle (2e4 / 1.8e5
  CNOTs) buys more as extra evolution.  The old framing "the cost of every route is set by the filter
  cycles it still needs" is replaced by "the cost is the transport, and the ramp's step is cheapest".
- The sixteen-site coarse eigenstate starts ahead (0.87 vs 0.72) and Path A transports to 1e-3 in a
  third of the ramp's time (16/t vs 48/t), but the ramp's 304-CNOT step against Path A's
  304 + 2×332 cancels the advantage.
- Single-layer UCCSD deteriorates by ×2.2 from eight to sixteen sites, not "two orders of magnitude".

**Manuscript length**: 1 085 → 1 171 lines of LaTeX (the frontier tables and the scope statements
of point 9 are the additions; the old fixed-input tables and the wrong plateau explanation are gone).
