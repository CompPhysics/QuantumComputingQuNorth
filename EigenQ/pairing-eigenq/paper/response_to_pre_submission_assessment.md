# Response to the pre-submission critical assessment (20 October 2026)

Assessment: `paper/pairing_pre_submission_critical_assessment.md`.  Every numerical claim
of the assessment was recomputed independently with the repository code before any change
was made; the recomputation lives in the new notebook section 22 (`scripts/frontier_checks.py`,
figure 24) and in the library (`pairinglib.averaged_map_floor`, `AveragedFilterMap`,
`seniority.fixed_step_evolution`, `seniority.coarse_term_cnots(convention="flag")`), with
tests in `tests/test_rodeo.py` and `tests/test_seniority.py`.  The previous main text is kept
as `paper/eigenq_pairing_v8_backup.tex`.

## Verdict per point

| # | Point | Verdict | What changed |
|---|-------|---------|--------------|
| 1 | Averaged digital-filter map unconverged (Gauss–Hermite); the sixteen-qubit "zero overlap" claim is a quadrature artifact | **Agree, confirmed.** Boundary-split Gauss–Legendre reproduces every number of the assessment (k=4: 1.1399e-4/0.99782, 1.9345e-3/0.96512, 0.15394/0.62970; k=8, dt=0.25: floor 2.2552e-3, eigenvalue 0.79450; dt=0.1: 4.838e-5/0.99324) and is stable under m=6→12 to 1e-10 with eigen-residuals < 1e-12.  The 60/80-node Hermite rule gives spurious dominant states at k=8 (floor 1.0 with 60 nodes, 0.9998 with 80). | Library routine with convergence check; Secs. V B, V E, discussion, Fig. 6 caption, SM Table S2 corrected; the zero-overlap explanation removed everywhere; dt=0.25 at k=8 is still rejected, now for the right reason (floor 2.3e-3 above every target, acceptance −21 %/cycle).  g=2 floors: 4.0e-4 (dt=0.1), 2.7e-5 (dt=0.05). |
| 2 | Size comparison mixes thresholds (16 vs 50) | **Agree.** Recomputed table matches (k=4/k=8 at 1e-3: T=40/50; at 1e-2: 16/30). | Sec. V E now states 40→50 at 1e-3 and 16→30 at 1e-2. |
| 3 | χ_s vs χ_u; the pair-counting formula is model specific, the Galerkin identity is not | **Agree.** χ_s max = 1.40/0.76/0.52/0.41 (factor 3.4) reproduced. | Sec. III C defines χ_u on the affine path and gives χ_s = (π/2 sin πs)² χ_u; Sec. IV D quotes both factors; Sec. V E now states the identity χ_u(0) = ‖Q H_high P ψ_low‖²/b² (new Eq. chi0) with P(k−P)g²/(4b²) as its pairing-seed value; abstract, introduction, discussion and SM Fig. S5 reworded; the ramp-gap sentence softened. |
| 4 | Fixed-preparation minima are not a cost frontier; "2.5–3 at every target" contradicts the tables (2.88/1.85/1.47) | **Agree.** The extra unfiltered ramps of the assessment are reproduced exactly (k=4: 1.3533e-3/9472, 2.0708e-4/18944, 8.9658e-5/25456; k=8: 1.2349e-4/163200, 8.6555e-5/217600).  The same holds for the refined route (k=8, μ_buf=4, T=6: n_s=60 → 1.29e-4 at 1.94e5, n_s=80 → 9.6e-5 at 2.58e5). | New joint search over (T, μ_buf, n_s, M), M=0 included, for both adiabatic routes at both sizes, with the EXACT averaged channel instead of 40 sampled schedules (new Table tab:frontier, Fig. 9 replaced by the channel-exact Fig. with ADAPT included).  Tables II and V are kept and labelled as fixed-preparation case studies; their minimum-cost columns are removed.  All ratios in abstract, introduction, Secs. V C/V E and conclusions rewritten from the new table. |
| 5 | Coarse-term penalty depends on compilation/workspace; scalar phase must be shown | **Agree.** `coarse_term_cnots` now includes the conditioned scalar phase (3→4: 348+2 = 350, as the paper already used; 6→8: 4872 → 4886) and offers the shared-flag convention (224 / 1028 / 44 CNOTs, needing c−1 clean work qubits). | Sec. V B states the ancilla-free reference compilation is not a lower bound and quotes the shared-flag counts; Sec. V E and the hardware paragraph carry the width caveat; the conclusion "a better coarse state costs more to embed than it saves" is now conditional on the compilation. |
| 6 | Operating-point crossings need uncertainty; use deterministic channel integration | **Agree.** The exact channel reproduces the assessment's checks (single-layer ansatz g=1, six cycles: 9.856e-5; ramp T=6 n_s=40, three cycles: 8.66e-5). | The new table and figure are channel-exact, so the crossings carry no sampling error; the sampled tables keep their bootstrap errors; the asymptote/accessibility qualifications are stated in Sec. V B. |
| 7 | ADAPT statements and Fig. 9 out of sync | **Agree.** g=2 ADAPT input (8.9e-5) is above the floor and improves ×3.2 in six cycles. | Sec. V D, Fig. 8 caption and the paragraph after Table V corrected; the sixteen-qubit figure regenerated with the ADAPT states; earlier ADAPT prefixes are admitted in the joint search. |
| 8 | Correct results to preserve | Preserved. | — |
| 9 | "Measurement-free"; early termination < 5 %; proxy wording; ramp gap vs coupling; SM Fig. S3 / Table S5 captions; Fig. 1 caption; repository URL; workspace vs 2k+1 | **Agree on all.** Early termination is 9.4 % for the determinant at M=6 (59.8k vs 66k), < 2 % for p > 0.99. | All wording fixed ("free of training measurements"; weights 3–8 in Fig. 1; SM captions; Table S5 locator).  The repository URL/DOI remain to be inserted by the authors (bibitem `Code`, `CITATION.cff`). |

## Points where we did not follow the assessment literally

- We did **not** replace Tables II and V by channel-exact numbers: they document the sampled
  protocol that a device would run (40 common schedules, bootstrap errors) and the new Table
  tab:frontier supersedes only their minimum-cost role.  Where both exist, the sampled and the
  exact values agree within the quoted bootstrap errors.
- We did not extend the joint search to the filter width σ and step δt: both were fixed by the
  spectral width (σΔ ≈ 10) and by the floor (δt such that the floor lies below the smallest
  target), and the text says so.
- The shared-flag compilation is quoted as a constructive estimate, as the assessment itself
  labels it; it was not validated by an explicit circuit simulation.

## What the revision changed in the conclusions

With a comparable step-count budget granted to the unfiltered adiabatic routes, the cheapest
adiabatic operating points at the targets 1e-3, 3e-4 and 1e-4 (k=8) and 2e-3, 5e-4, 2e-4 (k=4)
use **no filter cycle**: a filter cycle at sixteen qubits (1.2e5 CNOTs) costs as much as doubling
the step count of the adiabatic stage, and the latter buys more at these targets.  The ranking
is unchanged — ADAPT ≪ HF ramp < refinement (by the coarse-term overhead) < single-layer ansatz
+ filter < determinant + filter at sixteen qubits — but the statement "the filter cycles dominate
every budget" is replaced by the statement above, and the numerical factors are those of the new
table.  The filter's role at these sizes is the removal of the adiabatic/digitization tail below
~1e-4, between the unfiltered error and its own floor, and the diagnostic rejection of a
corrupted evolution.
