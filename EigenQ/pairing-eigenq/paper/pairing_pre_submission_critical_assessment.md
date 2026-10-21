# Critical assessment before submission: pairing resolution refinement and rodeo

Reviewed 20 October 2026. Sources: the current 22-page main article, *Resolution refinement followed by the rodeo algorithm: a gate-level benchmark of eigenstate-preparation routes on the pairing Hamiltonian*, and its eight-page supplemental material, both dated 16 October 2026. Page references below refer to the PDF/printed page numbers. The corresponding local implementation and executed notebook were also inspected. The manuscript, supplement and source repository were not modified.

This is a fresh assessment of these versions. Earlier concerns were checked against the current text rather than assumed to remain open. Several have now been explicitly addressed. Findings from the separate Hubbard manuscript were not presumed to apply to this paper.

## Overall assessment

**I recommend revision before submission, but the underlying construction is scientifically sound and the manuscript can become a useful controlled benchmark.** The model Hamiltonian, prolongation, shifted path, coarse-error budget, ideal rodeo ensemble formula, and accepted-preparation cost formula withstand the checks performed. The paper now distinguishes ideal filtering from digital filtering, includes meaningful competing inputs, and is appropriately candid that the small adaptive variational circuits win in the studied sectors.

The remaining issues are consequential rather than merely stylistic:

1. Some central digital-filter asymptotes are obtained from unconverged numerical averaging. At sixteen qubits this produces a qualitatively wrong claim that the dominant filtered state has zero ground-state overlap.
2. The size-dependent refinement-time comparison mixes two different fidelity thresholds.
3. The discussion switches between two differently parameterized response metrics without carrying the schedule factor through the interpretation.
4. The reported cost minima omit useful unfiltered digital ramp operating points, and several stated cost ratios contradict the tables.
5. The penalty attributed to a nontrivial coarse Hamiltonian depends strongly on the selected control compilation and clean-workspace budget.
6. Some figures, captions and abstract-level statements are inconsistent with the updated ADAPT results.

These issues do not invalidate the qualitative conclusion that small symmetry-adapted ADAPT circuits are economical on this pairing testbed. They do prevent the present numerical factors and some explanations of digital filtering from serving as definitive conclusions.

## 1. Highest priority: converge the averaged digital-filter map

**Locations:** main Section V B, p. 12; Section V E, pp. 15–16; discussion, p. 17; Figure 6; supplement Table S2. Implementation: the notebook's digital-map checks, corresponding to `assessment_checks.py` and the large-space analysis.

The mathematical formulation is correct. For the implemented variable-step evolution,

\[
\widetilde K_t=\frac12\left(I+e^{iE_{\rm obj}t}\widetilde U(t)\right),
\qquad
\mathcal K(\rho)=\mathbb E_t[\widetilde K_t\rho\widetilde K_t^\dagger],
\]

and independent times imply

\[
\rho_M=\frac{\mathcal K^M(\rho_{\rm in})}
{\operatorname{Tr}\mathcal K^M(\rho_{\rm in})}.
\]

A dominant positive eigenoperator describes the limiting accepted state, subject to the usual uniqueness/accessibility conditions. Its eigenvalue controls asymptotic survival per cycle. This is the appropriate description here because each cycle uses

\[
n(t)=\lceil |t|/\delta t\rceil,
\qquad \text{step size}=t/n(t),
\]

so different cycle times generally do not share one effective Hamiltonian's eigenvectors.

**The problem is the numerical evaluation of the expectation, not this map formulation.** The code uses a small fixed Gauss-Hermite rule, with 80 nodes in the small case and 60 in the large case. The ceiling rule changes the product formula at every boundary |t|=nδt. A global Gaussian quadrature rule is not automatically converged for that piecewise-defined, oscillatory integrand. In the large case, the sampled times can create artificial near-resonances with excited states.

### Eight system qubits

For k=4, N=4, g=1, σ=4, Eobj=E0, independent integration split at every ceiling boundary gives:

| Maximum step δt | Manuscript limiting infidelity | Boundary-split result | Manuscript eigenvalue | Boundary-split result |
|---:|---:|---:|---:|---:|
| 0.25 | about 1.2×10^-4 | 1.13988×10^-4 | 0.9979 | 0.99782077 |
| 0.5 | 1.5×10^-3 | 1.93451×10^-3 | 0.966 | 0.96512446 |
| 1.0 | 0.21 | 0.153936 | 0.64 | 0.62969672 |

The published values are reproduced by the original quadrature. Doubling the within-interval Gauss-Legendre order from 12 to 24 leaves the replacement values stable to approximately machine precision. The differences at δt=0.5 and 1.0 are too large to describe as rounding.

### Sixteen system qubits: a qualitative change

For k=8, N=8, g=1, σ=4, δt=0.25, Section V E states that the dominant eigenoperator has zero ground-state overlap and eigenvalue about 0.960. The 60-node calculation reproduces that behavior but is not converged with respect to either integration or, in some runs, state iteration. Increasing the Gaussian rule changes the apparent dominant branch drastically.

Splitting the integral at the ceiling boundaries instead gives

\[
\lambda_{\max}\simeq0.794498775,
\qquad
1-\langle\Psi_0|\rho_\infty|\Psi_0\rangle
\simeq2.25523\times10^{-3}.
\]

Six and twelve Gauss-Legendre points on each interval give the same values to the precision relevant here; the final normalized-state change is below 10^-10. The integration extends to eight Gaussian widths. The propagator used in this check was separately compared with the manuscript's product-formula routine at representative times, agreeing to approximately 10^-14.

Thus the limiting state's ground-state weight is about **0.997745**, rather than zero. The step is still too coarse for the paper's tight accuracy targets, and acceptance still decays substantially. What fails is the specific explanation that the continuous Gaussian ensemble locks onto a completely different state with zero target overlap. That is a quadrature artifact in the present calculation.

At the working δt=0.1, an interval-split check gives an asymptotic infidelity about 4.838×10^-5 and eigenvalue 0.9932403, close to but distinguishable from the quoted 4.6×10^-5 and 0.9933.

**Required revision:** recompute all reported map eigenvalues and floors with demonstrated integration convergence. Split at the ceiling boundaries or use another method shown to converge for this integrand. Report both a state/eigenoperator residual and convergence under integration refinement. A stable trace or eigenvalue estimate alone is not enough. Replace the false zero-overlap explanation in the abstract/discussion where it informs general claims about rejecting corrupted evolution.

The finite-schedule Table II/V tracks were generated separately; a wrong map asymptote does not by itself prove those tracks are wrong. Their statistical precision and limiting interpretation must be checked separately.

## 2. The size comparison uses unmatched accuracy thresholds

**Locations:** Section V E, p. 15, “What changes with the size”; compare Section IV D, pp. 8–9, and Table IV, p. 17.

The text says that with buffer b=0.6, the refinement time grows to T*=50 for infidelity below 10^-3 “instead of 16” at the smaller size. But T=16 at k=4 is the threshold for overlap 0.99, meaning infidelity 10^-2. These are different targets.

Independent pair-basis propagation of the stated path with 160 midpoint slices gives:

| T | k=4, N=4 infidelity | k=8, N=8 infidelity |
|---:|---:|---:|
| 16 | 9.0430×10^-3 | 3.4226×10^-2 |
| 25 | 2.4287×10^-3 | 1.00065×10^-2 |
| 30 | 1.2428×10^-3 | 5.1234×10^-3 |
| 40 | 3.7353×10^-4 | 1.4862×10^-3 |
| 50 | 1.3567×10^-4 | 5.3952×10^-4 |

On the stated grid, the matched 10^-3 comparison is **40→50**, not 16→50. At 10^-2 it is **16→30**. The endpoint response does grow fourfold at fixed buffer; that algebraic fact is independent of this timing comparison. The within-large-system improvement from tuning the shift also remains valid.

**Required revision:** use one target for every size comparison and distinguish these comparisons from the separate improvement obtained by retuning the shift. Do not use the unmatched thresholds to infer a size-dependent runtime factor.

## 3. Distinguish the affine-path and scheduled response metrics

**Locations:** Section III C, p. 4; Section IV D, pp. 8–9; supplement S3 and Figure S5.

Section III C defines a response using ∂sH, but the later calculation uses ∂uH with u=sin²(πs/2). These quantities are related by

\[
\chi_s(s)=\left[\frac\pi2\sin(\pi s)\right]^2\chi_u(s).
\]

The exact result χu(0)=g²/b² for the filled k=2 seed is correct. It is not the initial value of the quantity defined with ∂sH: χs vanishes at both scheduled endpoints.

| Buffer b | Maximum χu | Maximum χs |
|---:|---:|---:|
| 0.3 | 11.1111 | 1.3993 |
| 0.6 | 2.7778 | 0.7640 |
| 1.0 | 1.0000 | 0.5205 |
| 1.5 | 0.4444 | 0.4131 |

The factor 25 quoted across the scan belongs to the affine-path response. The scheduled maxima change by about 3.39. This does not rescue a claimed inverse-gap scaling law: the response still varies, and the text correctly says four threshold points do not identify an exponent. It does require consistent definitions when connecting the calculation to the cited argument, which is derived for a linear interpolation in the [original resolution-refinement supplement](https://arxiv.org/html/2511.14732v2).

The response quantity is a squared eigenvector response, or quantum metric. Its value is not by itself a dynamical leakage probability. Runtime also depends on the schedule and oscillatory gap-dependent phases. The present text is more careful about this than earlier versions; preserve that distinction.

### What generalizes beyond this model

The formula P(k−P)g²/(4b²) counts pair transfers from the filled coarse seed of this particular constant-pairing model. The general Galerkin identity for an exact coarse eigenstate is

\[
\chi_u(0)=\frac{\|QH_{\rm high}P|\psi_{\rm low}\rangle\|^2}{b^2},
\qquad Q=I-PP^\dagger.
\]

The response-based diagnostic generalizes. The specific pair-counting coefficient does not apply to arbitrary interactions or arbitrary correlated coarse states. Replace statements that the exact pairing formula itself is model independent with this distinction. Keeping an endpoint response of order one is a useful tested shift rule, rather than a sufficient general theorem that fixes the full transport time.

Also soften the statement that a minimum ramp gap of 2δ means the ramp does not slow with g. A gap alone cannot establish that; derivative matrix elements and the schedule still matter. The manuscript correctly makes this point for the refinement path and should apply it consistently to the ramp.

## 4. The resource comparison is improved, but is not a full cost frontier

**Locations:** Section V C, pp. 12–13; Section V E, pp. 15–18; Tables II and V; abstract and discussion.

Several elements are now fair and should be retained: the common observable is ground-state infidelity; M=0 is included; every input receives the same specified filter; repeated preparations and early failures are counted; and setup costs are explicitly separated from execution costs. The ordinary ramp and small pair/ADAPT ansätze are substantive baselines.

However, minimizing over filter cycle count for a few fixed digital preparations does not establish the least cost of the underlying method. The step counts are often chosen by reproducing an exact-slice overlap within a tolerance around 10^-3, while final targets can be 10^-4. Paying for more accurate unfiltered evolution can be cheaper than paying for subsequent controlled filtering.

Additional ordinary ramps, evaluated using the current pairing code and its own CNOT convention, demonstrate this directly:

| System | Unfiltered ramp | Infidelity | CNOTs | Relevant listed filtered-ramp cost |
|---|---|---:|---:|---:|
| k=4, N=4, g=1 | T=6, ns=16 | 1.35330×10^-3 | 9,472 | 17,000 for target 2×10^-3 |
| k=4, N=4, g=1 | T=8, ns=32 | 2.07082×10^-4 | 18,944 | 40,000 for target 5×10^-4 |
| k=4, N=4, g=1 | T=8, ns=43 | 8.96578×10^-5 | 25,456 | 64,000 for target 2×10^-4 |
| k=8, N=8, g=1 | T=6, ns=60 | 1.23487×10^-4 | 163,200 | 220,000 for target 3×10^-4 |
| k=8, N=8, g=1 | T=6, ns=80 | 8.65545×10^-5 | 217,600 | 480,000 for target 10^-4 |

These are actual product-formula fidelities against the exact ground state, not adiabatic estimates. The original anchors are reproduced: k=4, T=4, ns=11 gives infidelity 2.21843×10^-3 and 6,512 CNOTs; k=8, T=6, ns=40 gives 3.53046×10^-4 and 108,800 CNOTs. In the last two new examples the time is unchanged: finer digital evolution alone supplies the improvement.

The extra points are constructive counterexamples to the stated minima, not claims that this small additional grid finds global optima. At g=2 the checks did not establish an equally strong unfiltered advantage over the best listed filtered route at every target; do not extrapolate the g=1 result into an unsupported universal statement.

**Required revision:** search T, ns and M jointly, with M=0 allowed, and grant comparable calibration freedom to both adiabatic routes. Include shift and filter step/width in the disclosed search where feasible. Alternatively, label the tables explicitly as fixed-preparation case studies and restrict the abstract conclusions accordingly.

The strongest surviving execution-cost conclusion is that short symmetry-adapted ADAPT circuits are much cheaper in these particular sectors. What needs rebuilding is the quantitative comparison of filtering versus better digital preparation, and the inference that expensive filter cycles necessarily dominate an optimized adiabatic budget.

### A direct arithmetic inconsistency

The current Section V C reports ramp costs 17k, 40k and 64k at the three eight-qubit targets, and refinement costs 49k, 74k and 94k. Their ratios are

\[
2.88,\qquad1.85,\qquad1.47.
\]

The repeated claim “2.5–3 times more at every target” is therefore false even before broadening the parameter search. Replace it with target-specific numbers or the actual range, and update the abstract, introduction and conclusion together.

## 5. The embedded-coarse penalty depends strongly on compilation and workspace

**Locations:** Sections IV C, V B and V E; Table V; the coarse-term gate counter.

The quoted generic Pauli-rotation costs are internally meaningful reference-compiler costs. They are not lower bounds on what the embedded coarse Hamiltonian must cost.

For a nontrivial coarse space, every coarse rotation is conditioned on the same appended qubits being empty. The coarse evolution never changes those control qubits. A circuit can therefore compute their conjunction once into clean workspace, reuse a single flag across the entire coarse block, and uncompute it once. The relative phases of the compute/uncompute ladder cancel under the same conditions as the original ladder construction.

Using the relative-phase-Toffoli convention already available in the code gives the constructive shared-flag estimate

\[
C_{\rm block}=C_{\rm base}+2n_{\rm rot}+6(c-1).
\]

| Coarse application | Base CNOTs | Rotations | Appended controls | Quoted/implemented count | Shared-flag estimate |
|---|---:|---:|---:|---:|---:|
| 3→4 | 150 | 33 | 2 | about 350 | 222 |
| 6→8 | 732 | 138 | 4 | 4,872 | 1,026 |

At the Table V 6→8 operating point, 80 second-order steps and the 384-CNOT coarse preparation would then cost about **382,144** CNOTs instead of 997,504, under those workspace assumptions and before further optimizations. This is an estimate to validate with an explicit circuit, not a same-width replacement result. A c-control ladder needs c−1 clean work qubits; one can reuse the rodeo ancilla during refinement, but the four-control block needs additional workspace.

The filled 4→8 coarse model is already a single conditioned phase, so there is no repeated block condition to share. Its 254-CNOT Gray-code count is consistent with the specified ancilla-free convention. A clean-ancilla alternative in the code costs 44 CNOTs instead. For the 40-step preparation, the overhead over the ramp would then be about 3.2%, rather than 18.7%, at a larger workspace budget.

These alternatives do not show that refinement wins. They do show that the conclusion “a useful nontrivial coarse state costs more to embed than it saves” is conditional on the chosen compilation and allowed width. Report that tradeoff explicitly rather than attributing the whole overhead to the refinement principle.

### Complete conditioned phases need to be shown

The coarse gate counter prices hopping, Z and ZZ rotations. The scalar coefficient of Hlow−μ becomes a physical conditional phase after projection and cannot simply be discarded as a global phase. The 3→4 counter gives 348 while the paper uses 350; an additional phase on the two appended controls could account for that difference. For other blocks, demonstrate how the corresponding scalar phase is included or merged into the diagonal compilation. This is a request for a complete circuit check, not a confirmed missing-gate correction of a specific size. A shared-flag implementation handles this scalar as a phase on the flag.

The present pairing UCCSD counter uses actual Pauli dictionaries, so the incorrect full-span double-excitation weight found in the separate Hubbard implementation does **not** apply here. Do not import that issue into this assessment.

## 6. Statistics, target crossings, and what a fidelity floor means

The ideal Gaussian ensemble formula, Eq. (10), is correct. The paper also now explicitly reports the large schedule-sampling uncertainty at late cycles, uses acceptance weighting, distinguishes sampled from ideal values, and discusses energy detuning. These are substantial strengths.

The digital operating-point decisions still need uncertainty attached to the decision, rather than only to individual points. For example, the single-layer g=1 ansatz in Table V is assigned a six-cycle 10^-4 operating point at a sampled infidelity about 9.8×10^-5, very close to the threshold. A modest change can move the selected minimum by one expensive cycle.

An exact averaged-channel calculation of this current operating point gives 9.85644×10^-5 after six cycles and approximately 740,901 CNOTs per accepted state, consistent with the table. Thus this check does **not** establish that the listed crossing is false. It establishes how narrow the margin is, and that deterministic averaging is feasible in the 70-dimensional pair sector. The corresponding ramp T=6, ns=40 gives 8.65968×10^-5 after three cycles and approximately 477,015 accepted-preparation CNOTs, also consistent with its table entry.

Where practical, replace the finite-schedule digital averages by deterministic channel integration. Otherwise bootstrap the complete paired comparison, including target crossings and minimization over M. State which percent-level rankings survive that uncertainty. Shared random times are a good variance-reduction choice but do not by themselves certify an operating-point ordering.

Two mathematical qualifications should remain explicit:

- A limiting fidelity floor is an asymptote, not a theorem that every finite cycle must worsen an input already below that floor. Read such claims from the actual finite-cycle track.
- A dominant eigenoperator needs to be both numerically converged and accessible from the input under the channel. State uniqueness/accessibility assumptions when presenting the map result generally.

## 7. ADAPT statements and figures are not fully synchronized

**Locations:** Section V D, p. 14; Figure 8, p. 16; Section V E/discussion, pp. 17–18; Table V and Figure 9.

At g=1, the selected ADAPT input has infidelity 1.7×10^-5, below the digital-filter floor near 4.8×10^-5. The displayed filtering track degrades that already good state.

At g=2, however, the selected input has infidelity **8.9×10^-5**, above the quoted floor 2.5×10^-5. Table V shows an improvement to 2.8×10^-5 after six cycles, a factor about 3.18. Thus the statements that both inputs are below the floor, and that the filter cannot improve both, are false. The p. 14 “factor two in six cycles” is also inconsistent with the table.

The economically relevant conclusion can be retained: both unfiltered adaptive inputs already meet the displayed 10^-4 target at vastly lower cost than filtering. That statement does not require claiming that filtering cannot improve their fidelity.

Figure 9 is stale relative to Table V. Its left title/legend show ten inputs without ADAPT, and the right legend shows six without ADAPT; the text and table now contain eleven and seven. Regenerate the plot so it includes the adaptive winner that drives the headline conclusion. Figure 8 does include adaptive states, but does not repair the mismatch in the matched-comparison figure.

For a true method-level frontier, also allow earlier adaptive prefixes at looser targets. Table III already supplies some cheaper qualifying prefixes than the single selected final ADAPT input used for every target in Table V. This would strengthen rather than weaken the adaptive baseline.

## 8. Results that are correct and should be preserved

The following checks support the manuscript's scientific foundation:

1. **Pairing convention and spectra.** An independent full 70-state fermionic construction at k=4, N=4, g=1 gives E0=0.6355484735756, N-sector gap 1.8230702612753, seniority-zero gap 2.2998329531153, and HF ground-state overlap 0.8676375677137. These match the article. The convention EHF=2−g is correct for Eq. (1).
2. **Seniority and prolongation.** For the stated nested bare Hamiltonians with the same spacing and coupling, P†HhighP=Hlow. The current text correctly separates this identity from the independent question of large ground-state overlap. Pairing dynamics and the term-level product formulas preserve the relevant seniority sector.
3. **Shifted path.** Equations (3)–(4) correctly retain the appended-register projector. The initial-gap formula and saturation at a coarse excitation gap are correct. The warning that a shift alone does not establish an isolated branch along the whole path is appropriate.
4. **Coarse-state error budget.** Equation (5), (sqrt(εlow)+aad+η)^2, follows from the triangle inequality for the final component orthogonal to the target. The text correctly says unitary refinement does not purify coarse-state error and that different error sources must share the total budget.
5. **Adiabatic bound.** Equation (6) has the correct ε^-1/2 dependence for an infidelity target; boundary derivative terms vanish for the sin² schedule. Independent evaluation gives ||B−A||=7.6205693, maximum derivative norm 11.9703623 and integral 510.8425, consistent with Table S1. This is the [Jansen–Ruskai–Seiler bound](https://arxiv.org/pdf/quant-ph/0603175) applied to the stated path, with the hypotheses and ideal/digital distinction now clearly separated.
6. **Ideal filtering.** The Gaussian success factors, acceptance-weighted ensemble fidelity and large-gap factor-of-two reduction per cycle are correct. A fixed schedule need not obey the ensemble's 2^-M lower bound; the current text makes this distinction.
7. **Execution-cost accounting.** Equation (12) correctly includes preparation retries and the reach probability inside the expectation over cycle costs, so early failures avoid paying for later cycles. It is an appropriate expected-work metric for repeated accepted preparations.
8. **Gate accounting versus depth.** The paper correctly distinguishes gate-by-gate validation of primitives from sector-level application of compiled blocks, and explicitly says its numerical outputs are CNOT counts rather than scheduled device depths. Preserve that distinction.
9. **Limits of the testbed.** The current text acknowledges that fixed-N basis growth is not a many-particle scaling study, that the main seed can be a single determinant, and that a small seniority-adapted variational pool is unusually favorable here. Those qualifications make the comparison more credible.

## 9. Smaller corrections and scope qualifications

- **“Measurement-free” pipeline:** the refinement transport requires no training measurements, but rodeo preparation explicitly uses ancilla measurement and post-selection. Replace unqualified “measurement-free route” in the abstract and conclusion with “no variational training measurements during refinement,” or similarly precise wording. Keep the distinction from parameter calibration and energy search.
- **Early-termination size:** the claim in Section V C that early termination changes costs by less than 5% for every input is inconsistent with the determinant's six-cycle Table II values. Catt≈66k, P6≈0.854 and accepted cost≈70k imply an actual early-terminated numerator around 59.8k, roughly 9% below 66k. Equation (12) remains correct; the statement of negligible effect needs revision.
- **Accepted coherent cost:** the unweighted cost of executing every scheduled cycle is a useful full-length proxy, but is not exactly the mean coherent gate count conditioned on acceptance. That conditional average weights schedules by their acceptance probabilities. Retain the word proxy or report the conditional quantity when discussing accepted-run coherence requirements.
- **Gap versus coupling:** the ramp's nearly constant minimum gap cannot by itself prove coupling-independent runtime. Avoid using an argument for the ramp that the paper correctly rejects for refinement.
- **Supplement Figure S3:** its caption still says the convergence is on a roughly inverse-gap timescale. Qualify this to match the main text's conclusion that the scan does not identify a scaling exponent.
- **Supplement Table S5 caption:** it refers to a bottom shift-scan section that is not present in that table; the shift scan is in main Table IV. Correct the locator.
- **Figure 1 caption:** “160 of them” suggests 160 identical weight-four blocks, although the full 1,312-CNOT template contains strings with differing weights. Say 160 Pauli rotations with their actual weights.
- **Repository availability:** fill in the repository URL and archival DOI in Ref. [60], and archive the exact executed source used for the revised tables and figures. A publisher-supplied supplemental URL placeholder is a separate journal-production convention.
- **Hardware scope:** the current noise and logical-qubit discussion is appropriately cautious. Additional clean-workspace alternatives must be accompanied by updated qubit counts; they cannot silently be inserted into the claimed 2k+1-qubit specification.

## 10. What to do before submission

### Necessary corrections

1. Replace the unconverged digital-map integrals, verify state/eigenoperator residuals, and correct the false sixteen-qubit zero-overlap explanation.
2. Recalculate all size/runtime comparisons at the same infidelity target and define χs and χu consistently.
3. Expand the adiabatic operating-point search to include finer unfiltered ramps; then update resource tables and all abstract-level ratios together.
4. Correct the ADAPT below-floor statements, early-termination claim, stale Figure 9 and supplemental locators.
5. Provide a complete conditioned-coarse circuit accounting, including scalar phases and declared workspace, or clearly restrict the claimed penalties to the chosen reference compiler.
6. Add convergence/statistical evidence for operating-point crossings and complete the code archive references.

### Useful extensions, rather than mandatory barriers to a scoped paper

- Exercise a genuinely approximate nontrivial coarse state and report the full coarse-preparation/transport/filter cost with its actual error propagated.
- Test one seniority-breaking interaction if making claims about regimes where structured variational inputs become difficult. A benchmark paper on the seniority-conserving model need not solve that extension if its scope is explicit.
- Apply comparable compilation effort to all methods, including merged half steps, optimized excitation blocks, pair encoding and concrete connectivity.
- For hardware-efficiency claims, supply a common shot-cost model and an amortization study for training, energy search and calibration. The present execution-only comparison can stand as a distinct metric if described accordingly.

The most defensible framing is a reproducible small-model benchmark of refinement, filtering and competing preparations, together with a response diagnostic and a correct averaged-channel treatment of digital filtering. The paper does not need to demonstrate a refinement speedup to be scientifically useful. It does need numerical convergence and conclusions that follow quantitatively from the comparisons actually performed.

## Scope and evidence

The assessment used the two current PDFs, including visual inspection of the relevant complete pages, and read-only inspection of their supporting code and notebook. It included an independently written small-system Hamiltonian/theory check, reproductions using the provided propagation routines, independent boundary-split channel integration, and additional digital-ramp operating points. It did not rerun the full variational optimization campaign for every size, transpile all routes to a physical device, or simulate hardware noise. Compiler alternatives above are labeled constructive estimates pending complete circuit validation.

Numerical evidence and check scripts are retained in the chat workspace under `work/pairing-review-current/`, including `theory_check.json`, `filter_quadrature24.json`, `large_filter_interval.json`, `large_filter_interval01.json`, `finite_digital_ensemble.json`, `ramp_frontier.json`, and the associated scripts. These scratch files support the reported checks; the source manuscript and repository remain unchanged.

Reviewed PDFs:

- `/Users/mhjensen/Software/QuantumComputingQuNorth/EigenQ/pairing-eigenq/paper/eigenq_pairing.pdf`
- `/Users/mhjensen/Software/QuantumComputingQuNorth/EigenQ/pairing-eigenq/paper/eigenq_pairing_supp.pdf`
