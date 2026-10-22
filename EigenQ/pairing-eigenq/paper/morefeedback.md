# Critical pre-submission assessment of the revised pairing manuscript

**Reviewed 21 October 2026.** Main article: *Resolution refinement followed by the rodeo algorithm: a gate-level benchmark of eigenstate-preparation routes on the pairing Hamiltonian*, dated 20 October 2026, 24 pages. Supplemental Material: same date, eight pages. Page numbers below refer to these PDFs. The October 14 assessment was used as a checklist; the intervening October 20 assessment and response were also consulted. Their conclusions were checked rather than assumed.

## Recommendation

**The article is substantially closer to submission, and the central mathematical construction is sound. I recommend one focused numerical and editorial revision before submitting the present version to Physical Review C.** Most of the important October 14 objections have been addressed. The remaining substantive problem is the strength of the resource-optimization conclusions: the newly added search is too coarse to support several quoted minimum costs, small percentage margins, and claims about when filtering is cost-optimal. There is also a definite ADAPT prefix-search error in Table V.

These findings do not require abandoning the method or rebuilding the study. The paper can make a defensible contribution as a carefully controlled benchmark of state-preparation routes, including an informative negative result. Its strongest supported message is that resolution refinement improves the input overlap, but this pairing testbed offers inexpensive symmetry-adapted variational preparations and no established preparation-cost advantage for refinement over a well-tuned ordinary ramp. The precise ramp/refinement rankings and filter-cycle choices need a more converged comparison.

I would distinguish **scientific correctness of the construction**, which is now in good shape, from **readiness of the headline quantitative comparison**, which still needs work. I found no remaining mathematical defect comparable to the earlier missing error budget or unsupported inference of an inverse-gap exponent.

## Scope and evidence

The review covered both PDFs, their associated LaTeX source, the earlier assessment, relevant calculation routines, and selected primary literature. PDF pages were rendered and inspected, with full-page checks of the principal equations, the new frontier table, discussion, and supplemental bound table. The manuscript and repository were not modified.

Independent checks reconstructed the seniority-zero Hamiltonians, exact spectra, midpoint adiabatic evolutions, and second-order pair-hopping product formulas. Selected new operating points were then cross-checked with the repository's own routines. The small-system digital-filter channel and sufficient bounds were reconstructed independently. ADAPT runs were repeated to check the prefix and sixteen-qubit comparisons.

This was not a clean rerun of every notebook cell or figure, a full audit of all bibliography entries, or a simulation of hardware noise, routing, fault-tolerant rotation synthesis, or training shots. In particular, the sixteen-qubit digital-map eigenoperator was not independently recomputed here. Statements below distinguish reproduced results from conclusions that depend on the reported calculations.

## 1. What has been resolved since October 14

The changes are substantial and should be retained.

| Earlier concern | Current assessment |
|---|---|
| Approximate coarse states treated as exact | Addressed in Sections III A/C and Eq. (5): coarse infidelity, adiabatic leakage, and digitization share an amplitude budget. Refinement is correctly described as unitary rather than purifying. |
| Prolongation assumed to preserve energy generally | Addressed in III B: the additional Galerkin identity is stated, with expectation-value preservation separated from fine-space fidelity. |
| Shift assumed to guarantee the entire path remains gapped | Addressed in III C: branch isolation is an explicit condition, and the shift is correctly interpreted as a relative penalty. |
| Multigroup Trotter bound ambiguously used a largest commutator | Addressed: the sum of pairwise commutator norms is explicit. |
| Optimizer-free execution confused with practical calibration | Addressed: exact-answer calibration of time, shift, and step count is acknowledged; an operational stopping protocol is not claimed as demonstrated. |
| Four-point scan claimed to establish inverse-gap scaling | Addressed in IV D/S3: the inverse-square fit and changing response are acknowledged, and the scan is no longer treated as an exponent measurement. |
| Ideal and digital filtering conflated | Addressed in V B and Table S2: the averaged digital channel, biased asymptote, and decaying acceptance are distinguished from exact on-resonance filtering. |
| Forty schedules treated as precise at very small infidelity | Substantially addressed: accepted-ensemble weighting, statistical uncertainty, and deterministic channel evaluation are now used in the appropriate places. |
| Variational preparation charged for an unnecessarily large template | Addressed: the 192-CNOT four-pair comparator is explicit, and ADAPT is actually run. |
| Ordinary ramp and direct filtering of a useful coarse state omitted | Addressed: both are included. |
| Nontrivial coarse evolution not simulated at its counted splitting | Addressed for the reported nontrivial chain, including the pruned coarse preparation and conditioned scalar phase. |
| Fixed-cycle comparisons mistaken for matched-accuracy comparisons | Major improvement: zero filter cycles and a joint control search are included. Search convergence remains the principal unresolved issue. |
| Gate counts, coherent depth, retries, and setup mixed together | Substantially addressed: CNOT counts and accepted-preparation work are distinguished from depth and one-time costs. |
| Orbital growth at fixed particle number overinterpreted | Improved by the four-pair/eight-level calculation. This is a useful fixed-filling enlargement, while still insufficient to establish asymptotic scaling. |

The ideal-filter qualifications are also now correct: the ensemble background differs from a fixed schedule; energy degeneracy selects an eigenspace; detuning reduces target survival; and the peak width depends on both time width and cycle count. The coarse-to-fine energy identity, combined error budget, and digital-channel formulation withstand direct examination.

## 2. Main remaining issue: the cost frontier is not converged

**Locations:** Table V, p. 18; Section V C, pp. 13–14; Section V E, pp. 17–18; abstract, introduction, and discussion/conclusions.

Table V does disclose finite control grids. Consequently, most entries can be legitimate minima *within those grids*. The concern is that the prose goes further: the paths are said to travel at the same speed and differ only by coarse-term overhead, and particular filtered points are used to describe the cheapest route at a target accuracy. The grid is too sparse to justify these interpretations at the quoted precision.

Simple omitted evolution times already produce materially cheaper preparations. The following use the paper's own midpoint schedule, term order, unmerged second-order CNOT convention, and encoding. All have **zero filter cycles**, so neither random-schedule error nor digital-channel quadrature enters their infidelity or cost.

| System and route | Target infidelity | Added operating point | Independently reproduced infidelity | CNOTs | Table V cost at that target |
|---|---:|---|---:|---:|---:|
| Eight qubits, HF ramp, g=1 | 2×10⁻³ | T=4.5, n_s=11, M=0 | 1.63292×10⁻³ | 6,512 | 9,472 |
| Same | 5×10⁻⁴ | T=4.75, n_s=16, M=0 | 4.34611×10⁻⁴ | 9,472 | 18,944 |
| Same | 2×10⁻⁴ | T=4.75, n_s=22, M=0 | 1.57748×10⁻⁴ | 13,024 | 25,456 |
| Sixteen qubits, HF ramp, g=1 | 3×10⁻⁴ | T=5, n_s=40, M=0 | 2.23472×10⁻⁴ | 108,800 | 163,200 |
| Same | 10⁻⁴ | T=5.2, n_s=60, M=0 | 6.15426×10⁻⁵ | 163,200 | 217,600 |
| Sixteen qubits, refinement, g=2 | 10⁻⁴ | b=16, T=10, n_s=160, M=0 | 8.94159×10⁻⁵ | 516,480 | About 570,000, with M=2 |

Every row was calculated with an independent pair-basis implementation and confirmed with the authors' routines. The cost formulas are 2 n_s×296 for the eight-qubit ramp, 2 n_s×1,360 for the sixteen-qubit ramp, and 2 n_s×(1,360+254) for the stated sixteen-qubit refinement. These examples are constructive counterexamples to treating the published envelope as converged; they are not claims to have found global optima.

There are three consequences.

1. **The small ramp/refinement margins are unstable.** At eight qubits, the revised ramp points are cheaper than the published refinement winners at the two tighter thresholds. The published refinement advantage therefore cannot be interpreted without giving both routes comparable local optimization. Even on the original grid, calling 14k versus 19k and 20k versus 25k a “tie” needs an explicit meaning: these are deterministic differences of roughly 20–26%, not statistical ties.
2. **The exact 19% coarse-term explanation is conditional on the selected operating points.** For the published sixteen-qubit g=1 winners, the step counts coincide, so their cost ratio is indeed set by the chosen coarse-term compilation. Nearby ramp times change the required step count. The numerical observation does not establish dynamical equivalence of the two paths or a general coarse-term-only separation.
3. **The g=2 tightest target does not require filtering at the advertised best cost.** Adding T=10 gives an unfiltered refinement with infidelity below 10⁻⁴ and cost below the listed two-cycle preparation. This does not prove that filtering is never useful. It shows that its claimed optimal role at this point depends on the coarse search grid.

**Recommended repair:** perform a local refinement of T and b around every competitive target crossing for both adiabatic routes, and search integer step counts near those crossings. Include every M from zero to the stated cap. Show how the best-found cost changes under a finer search. The small six- and seventy-dimensional blocks make this a bounded numerical task.

A certified global optimum is unnecessary. “Best cost found under the stated search” is a suitable description, provided the search is sufficiently converged for the conclusions drawn. If retaining the current grid, restrict the article to a finite-grid benchmark and remove narrow quantitative claims that imply general optimized rankings. The existing restriction to fixed filter width and step can remain explicit; this review does not demand a universal optimization of every possible filter.

## 3. A definite Table V error: the eight-qubit ADAPT prefix

**Locations:** Table III, p. 16; Table V, p. 18; the control-search implementation.

Table V lists generalized-pool ADAPT as five operators/240 CNOTs at all three eight-qubit targets, while its caption says prefixes are considered. The four-operator prefix already has

\[
1-F=1.7764427147\times10^{-4},\qquad C=4\times48=192.
\]

It meets every listed target: 2×10⁻³, 5×10⁻⁴, and 2×10⁻⁴. Table III and Section V D are consistent with that result. The frontier script admits only the five-operator state for this small-system ADAPT input, unlike its treatment of larger-system prefixes.

**Correct the eight-qubit Table V ADAPT row to 192 CNOTs/four operators, with M=0, and admit its prefixes in the actual search.** Five operators/240 CNOTs remains appropriate for the tighter 10⁻⁴ target, which Table III considers but the eight-qubit panel of Table V does not.

This correction strengthens ADAPT and makes it tied with the four-pair preparation at these thresholds. It does not threaten the broad variational comparison.

## 4. Specify the random-time protocol consistently

**Locations:** III D/Figure 2, p. 6; V B, p. 12; V C, pp. 12–14; the digital-map and sampled-filter implementations.

The prose sometimes says cycle times are Gaussian without specifying a sign convention. The implemented digital tracks and deterministic map consistently use positive times,

\[
t=|z|,\qquad z\sim\mathcal N(0,\sigma^2).
\]

This is a half-normal distribution. The existing calculations are internally consistent with that protocol. The issue is specification, not a contradiction between the reported sampled and deterministic implementations.

With ideal evolution, the relevant survival factors are even in t, so the ideal population and fidelity formulas remain unchanged by taking the absolute value. With product-formula evolution, the entire density-operator channel need not be the same. For signed Gaussian times, reducing integration to positive t requires

\[
\mathcal K_{\rm signed}(\rho)=\int_0^\infty f_{\rm half}(t)
\frac{\widetilde K_t\rho\widetilde K_t^\dagger+
\widetilde K_{-t}\rho\widetilde K_{-t}^\dagger}{2}\,dt.
\]

The implemented positive-time map has only the first channel term. Independent calculations illustrate the distinction: at eight qubits, σ=4, and step ceiling 1, the positive-time asymptotic infidelity is about 0.153936, whereas signed Gaussian times give about 0.129160. At the useful step ceiling 0.25 the difference is much smaller: approximately 1.13988×10⁻⁴ versus 1.13846×10⁻⁴.

**The straightforward repair is to define the implemented positive-time convention explicitly and use it consistently in captions and methods.** No rerun is needed simply to describe the existing protocol correctly. If the authors intend a signed-time algorithm instead, both signs must be included and the corresponding digital results recomputed. Do not claim the two channels are interchangeable merely because the ideal cosine-squared survival factor is even.

## 5. Results independently reproduced and worth preserving

The remaining criticism should not obscure the amount that now checks out.

**Basic model and digital evolution.** For k=4, N=4, g=1, independent diagonalization gives E₀=0.635548473575598 and HF squared overlap 0.867637567713702. The T=30, 80-step digital refinement gives p=0.997720914627350. The published ramp and refinement anchors at both sizes were also reproduced, including the sixteen-qubit g=1 ramp at T=6, n_s=40 and the g=2 ramp at T=8, n_s=80.

**Sufficient bounds.** Supplemental Table S1 agrees with an independent calculation: ||B−A||=7.6205693, the Eq. (6) integral is 510.84253, Λ₁=11.97036, Λ₂=13, and the first-order bound at T=16 and η=0.05 requires about 34,238 steps. The reported separation between the sufficient bound and empirical performance is real. Be more cautious about attributing the entire separation to one cause: the supplement correctly says its magnitude alone does not isolate the mechanism, while the main discussion sounds more definitive.

**Digital-channel floors.** Piecewise integration and direct diagonalization in the small-system pair sector reproduce the revised positive-time asymptotes:

| Maximum step | Dominant eigenvalue | Asymptotic infidelity |
|---:|---:|---:|
| 0.25 | 0.9978207657 | 1.1398755×10⁻⁴ |
| 0.5 | 0.9651244641 | 1.9345114×10⁻³ |
| 1 | 0.6296967173 | 0.1539357 |

The dominant-eigenoperator description is appropriate under the stated uniqueness and accessibility conditions. The revised paper properly separates an asymptotic floor from a finite-cycle track, including inputs that start below that floor. The earlier unconverged large-system quadrature explanation has been withdrawn; this review independently checks the small-system replacement, not every large-system map value.

**Ideal filtering and schedule statistics.** The exact accepted-ensemble infidelities from the T=25 input reproduce the reported values, including about 3.804×10⁻⁵ after six ideal cycles, 2.378×10⁻⁶ after ten, and 5.944×10⁻⁷ after twelve. An independent second-moment calculation gives relative standard errors of approximately 12.7%, 36.6%, and 85.6% for a forty-schedule residual-weight mean at M=2, 6, and 10. This supports the manuscript's warning about large relative uncertainty at late cycles. Those percentages concern schedule averaging; they are not a complete hardware shot-noise model.

**Larger-space ingredients.** All four coupling rows of Table S5 were independently reproduced for the ground energy, seniority-zero gap, HF overlap, and selected exact-slice ramp/refinement overlaps. Thus the sixteen-qubit extension has substantive numerical support. It should remain a benchmark at these sizes and parameters, rather than evidence of asymptotic advantage.

**ADAPT.** The sixteen-qubit generalized-pool results reproduce: 19/23/27 operators at g=1 give infidelities about 1.516×10⁻⁴, 1.702×10⁻⁵, and 4.145×10⁻⁶; at g=2, 25/29/31 operators give about 9.049×10⁻⁴, 1.926×10⁻⁴, and 8.896×10⁻⁵. The resulting preparation circuits are much cheaper than the tested adiabatic/filter routes. Their exact-gradient training and symmetry-adapted pool remain important conditions, already acknowledged in the paper.

## 6. Interpretation and presentation corrections

These are smaller than the frontier issue, but should be fixed in the same pass.

1. **Rewrite the abstract's ADAPT cost statement.** At target 10⁻⁴ the reported costs are 240 CNOTs at eight qubits and approximately 1,100–1,500 at sixteen. “About a thousand at both sizes” is inaccurate. “Two orders below every other route” also needs the target, system, and compilation specified; at the looser eight-qubit Table V targets the four-pair preparation has the same corrected 192-CNOT cost. Prefer actual costs or scoped ratios over a universal phrase.
2. **Correct the averaging-bound arithmetic on p. 6.** The bound 1/p−1 is about 15% at p≈0.87, rather than “a few percent.” A smaller observed discrepancy can be reported separately, but it does not follow from that bound. Around p=0.99 the bound is approximately 1%, not strictly below 1% for every p>0.99.
3. **Keep the general response identity distinct from a global runtime theorem.** Equation (13) is sound under its exact-coarse-state/Galerkin assumptions. It is an initial-response identity. Keeping it of order one is useful tested shift guidance; it does not alone bound the maximum response or runtime along the whole path. The paper now makes this distinction and should preserve it consistently in its headline language.
4. **Replace the “only locally” analogy to nucleation on p. 21.** The pairing force connects all pairs of levels, and the coarse term contains a many-qubit projector. The algorithms share incremental adiabatic construction; this example does not share the same spatial locality as lattice nucleation. “Incrementally” would be more accurate.
5. **Qualify the introductory universal 1/p cost statement.** Repetition without amplitude amplification has the stated inverse-success-probability cost. It is not a universal lower bound covering all QPE/QSP/state-preparation algorithms under every access model. Coherent amplitude amplification changes the overlap dependence under additional assumptions; see [Lin and Tong, *Near-optimal ground state preparation*](https://quantum-journal.org/papers/q-2020-12-14-372/). This does not require adding every algorithm as a numerical baseline.
6. **Keep compilation-dependent margins conditional.** The reference CNOT accounting is transparent, and the ancilla/workspace caveat is now explicit. Shared flags, merged half-steps, optimized excitation circuits, and routing can change small margins. A 19% reference-compiler difference should not become a hardware-independent conclusion. No new routing study is needed for the present scoped gate-count benchmark.
7. **Temper the hardware-specification language.** The paper correctly labels counts rather than scheduled depths and excludes several implementation costs. Calling the workload a logical gate-count benchmark is supported. A complete machine specification would additionally require connectivity, scheduling, rotation-synthesis accuracy, and a noise/error-correction model. Hardware claims should also be dated to their cited experiments rather than presented as a comprehensive statement about present capabilities.

## 7. Reproducibility and submission presentation

The code-availability reference still contains repository URL and Zenodo DOI placeholders. Complete the actual access information before claiming that readers can obtain and reproduce the implementation. A frozen release or commit should identify the exact version used for the final manuscript. An archive is useful; the essential immediate point is truthful, usable availability information.

Include an explicit Data Availability statement covering the numerical data and code, as requested by the [current PRC author instructions](https://journals.aps.org/prc/authors). The current Code Availability paragraph can be incorporated into that statement. Do not invent an archive identifier before the archive exists. The publisher-supplied Supplemental Material URL placeholder is a different matter and is normal before publication.

The PDFs are readable and the principal tables and figures are present. The existing compilation logs contain duplicate citation/label warnings in the supplement and float-placement warnings in the main article. These did not establish a scientific error or an obviously missing item in the rendered PDFs. Resolve them during the final build and check cross-document links. The supplement should be intelligible and correctly numbered when opened independently.

The discussion currently repeats several conclusions and caveats. A modest shortening would help readers identify the new contribution. Retain the ideal/digital operating-point table, the response distinction, and the matched-accuracy definition; they prevent real misunderstandings. Preserve the nuclear-observable checks in Table S4 and consider drawing one concrete occupation/correlation-energy result into the main discussion to make the nuclear-structure relevance more visible.

## 8. Suitability for Physical Review C and a practical route to submission

**PRC is a plausible venue for a carefully framed nuclear many-body benchmark.** Its [published criteria](https://journals.aps.org/prc/about) emphasize a significant, substantive addition to nuclear physics. This is an editorial assessment, not a prediction of acceptance.

The paper should lead with what the present study adds beyond the already published [resolution-refinement proposal](https://doi.org/10.1016/j.physletb.2026.140363) and the related [nucleation work](https://arxiv.org/abs/2609.37527): an explicit counted implementation, an informative comparison against strong inputs, a response-based explanation of shift sensitivity, and a quantified digital-filter bias. Simply combining known refinement and filtering procedures would be a weaker novelty case. The controlled comparison and the lessons from its failure to outperform simpler alternatives make the stronger case.

The integrable seniority-conserving model and very small classical blocks limit the generality, but are legitimate choices for validation. A seniority-breaking or realistic shell-model example would strengthen significance; it is not a necessary correction to make the present benchmark mathematically publishable. Nor should every alternative encoding, variational algorithm, or filter be demanded before submitting a clearly limited study.

A focused final revision should:

1. Refine the control search around the competitive points, quantify its stability, and regenerate Table V and the associated envelope figure.
2. Fix the eight-qubit ADAPT prefix row and include its prefixes in the implementation.
3. Specify the positive-time protocol consistently.
4. Update every dependent abstract, percentage ratio, filter-cycle claim, and concluding statement from the final comparison.
5. Complete usable data/code access information and check the final independent supplement build.

After those changes, I would regard the article as a credible candidate for submission as a controlled benchmark with conditional resource conclusions. In its current form I would revise first, because the quantitative comparison is central enough that the readily reproducible counterexamples are likely to attract justified referee criticism. The evidence supports a focused repair, rather than a fundamental objection to the physics or the hybrid construction.

---

**Source versions reviewed:** `paper/eigenq_pairing.pdf` and `paper/eigenq_pairing_supp.pdf`, both dated 20 October 2026, in `/Users/mhjensen/Software/QuantumComputingQuNorth/EigenQ/pairing-eigenq/`; `paper/article_correctness_assessment.md` dated 14 October; associated source and selected calculation routines. External primary-source links above were checked during this assessment. Independent calculations are review checks, not additions to the authors' reported dataset.
