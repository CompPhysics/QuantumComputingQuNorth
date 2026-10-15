# Response to `article_correctness_assessment.md` (revision of 2026-10-14)

Every quantitative claim of the assessment was recomputed independently
(`scripts/verify_assessment.py`) and found correct.  The revision adds the
tagged notebook section `EXT:ASSESS` (Sec. 19, code in
`scripts/assessment_checks.py`, inserted by `scripts/extend_notebook_assess.py`,
figures 19–20), patches three original cells (energy scan, rodeo convergence,
shift scan), adds `pairinglib.trotter_path`, and rewrites the affected parts of
`paper/eigenq_pairing.tex` (previous version: `eigenq_pairing_v4_backup.tex`).
Numbers below are notebook outputs.

## Verified (assessment claim → our recomputation)

| Claim | Recomputed |
|---|---|
| E0 = 0.6355484735756, p(HF) = 0.8676375677137 | identical |
| seniority-zero excitation gap at s=1: 2.29983 (sector gap 1.823) | 2.29983 / 1.8231 |
| shift-scan fits: 17.05/Δ − 8.43 (SSR 28.2) vs 4.13/Δ² + 4.16 (SSR 0.22) | 17.05/Δ − 8.42 (28.18); 4.13/Δ² + 4.16 (0.22) |
| χ_u(0) = g²/b²: 11.11 … 0.444, maximum at the initial endpoint | identical for b = 0.3, 0.6, 1.0, 1.5; max at s = 0 |
| derivative couplings at s = 0.25 to the two upper pair-sector states 0.1115, 0.0815; ‖B−A‖ = 7.6206 in the seniority-zero block | 0.11152, 0.08148; 7.6206 (same as full sector) |
| slice convergence T=25: 0.9975712826 (160) → 0.9975698271 (1280); T=30 likewise | identical to all digits |
| sampled rodeo (seed 7 reset): 1.350e-3, 6.49e-4, 2.085e-5, 1.985e-6 vs exact 1.216e-3, 6.08e-4, 3.80e-5, 2.38e-6 | identical; notebook's advanced generator gave 5.45e-6 at M=10 |
| relative standard errors 12.7 / 36.6 / 85.5 % at M = 2/6/10 | 12.7 / 37.0 / 82.4 % (20 000 schedules) |
| 12 cycles for ensemble infidelity ≤ 1e-6 | 12 |
| four-pair ansatz: E = 0.6369869810446, 1−F = 1.7764e-4, 192 CNOTs | identical; the other 22 UCCSD amplitudes are exactly 0 |
| k=3 coarse state with two pair rotations: 96 CNOTs | 1−F = 0 to machine precision |
| 3→4 chain: p(T=8) = 0.990737; split-coarse symmetric formula 0.991172 at 22 steps; 28 696 CNOTs | 0.990738; 0.991172; 28 696 (272-CNOT template) / 28 520 (96) |
| digital-filter map: floor 1.14e-4, eigenvalue 0.99782 (δt 0.25); 1.93e-3 (δt 0.5) | 1.18e-4, 0.99786; 1.52e-3 (quadrature on the step function limits agreement to ~10–20 %); δt = 1: 0.21, eigenvalue 0.64 |
| cost table rows 23 000 / 48 700 / 69 700 / 94 000 / 117 000 | 22 979 / 48 700 / 69 707 / 94 035 / 117 028 |

## Changed in the paper

1. **Pipeline assumptions and error budget (Sec. III).** Factorized/disentangled UCCSD with its ordering stated; both harmonics in the parameter dependence; singles vanish *in this optimization*; coarse input written as an approximation with ε_low; prolongation identity conditioned on P†H_high P = H_low and separated from the high-overlap claim; isolated-branch condition and the shift as a relative penalty; the combined bound 1−F ≤ (√ε_low + a_ad + η)² (new Eq. budget); Λ₂ defined as the sum over pairs of groups; "optimizer-free in execution, calibrated against the exact answer" stated explicitly.
2. **Rodeo text (Sec. III.D).** 2^−M is a lower bound for the ensemble, not for one schedule; width 1/σ at fixed M and 1/(σ√M) otherwise; on-resonance/exact-evolution scope of "P_M → p", detuned retained weight p q₀^M; degenerate targets; schedules drawn per attempt; the O(log(1/δ)/(pε)) estimate is not the gate complexity of the pipeline.
3. **Scaling interpretation (Sec. IV.D, Fig. 11, new Fig. 19).** Both affine fits with intercepts and residuals; "incompatible with 1/Δ²" and "confirms 1/Δ" removed; the scans presented as sensitivity studies; χ_u(0) = g²/b² and χ along the path shown; the coupling-scan causal wording replaced; three gaps distinguished (shift-controlled 0.6, sector 1.823, seniority-zero 2.300).
4. **Bounds vs practice (Sec. IV.E).** "never couples" corrected (couplings 0.60, 0.37, 0.00, 0.11, 0.08 at s = 0.25; ‖B−A‖ identical in the accessible block); operator-norm vs overlap criteria of Table II labelled.
5. **Rodeo statistics (Sec. IV.F, Fig. 8).** Exact ensemble values quoted as the reference; sampled values with standard errors from 20 000 schedules; the 45 % discrepancy at M = 6 explained; seed convention fixed (reset per input) and the earlier 5.5e-6 acknowledged; "of order ten" cycles with 12 for 1e-6; fine scan centred on the coarse-scan peak with shot-noise uncertainty (±0.02 at 10³ shots, ±0.006 at 10⁴); detuning statement limited to the tested regime.
6. **3→4 chain (Sec. IV.C).** T* = 8 on a finer grid; validated with the coarse term split into its coarse groups (22 steps, p = 0.9912); pruned 96-CNOT coarse preparation; "within 10 %", not an optimized minimum.
7. **Digital floor (Sec. V.B).** The averaged map K(ρ) and its dominant eigenoperator replace the "mixture of effective ground states" picture; floors 1.2e-4 / 1.5e-3 / 0.21 and eigenvalues 0.9979 / 0.966 / 0.64; "never converges" reworded.
8. **Gate budget (Sec. V.C, Table III rewritten, Fig. 17 → Fig. 20).** Eight inputs incl. the four-pair ansatz, the prolonged k=3 state and the HF ramp (D → H, T = 4, 11 steps, 6 512 CNOTs, p = 0.9978); M = 0 included; CNOTs per attempt and per accepted preparation; Eq. (cost) with the expectation inside; minimum cost at targets 2e-3 / 5e-4 / 2e-4; the "same floor after six cycles" claim corrected; the fixed-cycle "factor three" withdrawn; setup costs and amortization discussed; new operating-point table (Table IV).
9. **ADAPT remarks (Sec. V.D).** The 192-CNOT result is a demonstrated restricted ansatz, distinguished from an unperformed ADAPT run.
10. **Abstract, introduction, discussion.** Conclusions restated: refinement is tied with filtering the determinant, costs 2.5–3× the HF ramp at every target and 200–500× the unfiltered four-pair preparation; the inverse-gap statement withdrawn; "none of the stages uses seniority" qualified; outlook lists pair-encoded and symmetry-projected baselines and equal compilation effort.

## Not done (open)

- ADAPT-VQE proper; pair-encoded (one qubit per level) variational/filter baseline; symmetry-projected BCS input.
- A global search over refinement time, shift and step size for the minimum-cost curve; shot-budget model for training and energy search.
- Hardware routing on a concrete connectivity; rotation-synthesis / fault-tolerant resource model.
- A device-implementable calibration protocol for (T, μ, n_s).
- Companion paper II has not been re-audited against these points (its Sec. on bounds and its fixed-cycle comparisons inherit items 4 and 8).
