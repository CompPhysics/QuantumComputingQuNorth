# Critical assessment before submission: Hubbard resolution refinement and rodeo

Reviewed 19 October 2026. This assessment covers the 16-page main manuscript, *Lattice refinement followed by the rodeo algorithm: a gate-level benchmark of eigenstate-preparation routes on the Hubbard chain*, and its seven-page supplemental material, dated 16 October 2026. Page references below refer to PDF pages. The corresponding local notebook, generator script, Hubbard module and relevant tests were inspected read-only. No manuscript or repository source was changed.

## Recommendation

**Revise before submission.** The Galerkin construction, symmetry obstruction, small-system spectra, and fixed-step digital rodeo formulation provide a sound scientific foundation. Several numerical anchors were reproduced. However, the resource conclusions cannot presently be treated as a reliable comparison of the methods' achievable costs.

Four issues require correction and new results: a wrong constant in the 32-qubit Path-A implementation; unnecessarily repeated multi-control computation that greatly inflates Path A's quoted cost; an error in the Jordan-Wigner weights used to price double excitations; and omission of longer, substantially cheaper unfiltered ramp operating points. A fifth issue is the use of exact coarse input states while charging for an approximate preparation whose error is consequential at the reported targets.

These issues act in different directions. Correcting Path A's compilation makes refinement more competitive; extending the ordinary ramp baseline makes unfiltered preparation more competitive. The appropriate response is to rebuild the cost-versus-accuracy frontiers, rather than predict the final winner from the present tables.

## 1. Confirmed implementation error in the sixteen-site refinement

**Locations:** main Section 6, p. 12; supplement S7, pp. 5–6; Table 3, Figure 7 and Table S6 for the Path-A input. Source: `scripts/build_notebook_hubbard.py`, around line 934, the corresponding notebook cell, and the diagnostic script `scripts/hubbard16_checks.py`.

The circuit-level calculation calls

```python
tl, Ul, cst = hb.galerkin_params("chain", 1.0, U, 8)
```

The fourth argument of `galerkin_params` is **particle number**, not coarse site count. Here there are four particles, so the constant in Eq. (4) is −4t. Passing 8 gives −8t. The direct Galerkin matrix and the split coarse generator therefore differ by −4t times the coarse identity. After embedding, the implemented path differs from Eq. (6) by

\[
\delta H_A(s)=-4t\cos^2(\pi s/2)PP^\dagger.
\]

This is not an irrelevant global phase: it changes the relative energy of the prolonged subspace and its complement. A global constant in a coarse Hamiltonian becomes physically consequential when embedded with a projector unless the full shift convention is adjusted consistently.

The published discrepancy was reproduced, followed by a calculation changing only this argument from 8 to 4:

| Sixteen-site Path A, U=4t | Published implementation | Corrected constant |
|---|---:|---:|
| T=8/t, 128 steps: ground-state overlap | 0.990784757 | 0.992367850 |
| T=8/t, 512 steps | 0.990788581 | 0.992368342 |
| T=16/t, 128 steps | 0.998346666 | 0.999298441 |
| T=16/t, 256 steps | 0.998357777 | 0.999332878 |

For comparison, direct propagation of Eq. (6) at T=8/t gives 0.992373704 with 48 midpoint slices and 0.992368302 with 192 slices. The corrected 512-step result agrees with the converged direct result to about 4×10^-8 in overlap. The coarse matrix identity is restored to about 2×10^-15.

Thus S7's explanation that splitting noncommuting coarse groups follows a different limiting path is **incorrect**. For a consistent decomposition, the product formula converges to the same time-dependent evolution as its step shrinks. Here the persistent discrepancy comes from different Hamiltonians. The existing small-system checks miss the mistake because their coarse site count happens to equal their particle count.

**Required revision:** fix the particle-number argument, add a test with unequal site and particle counts that checks the sum of the embedded generators against Eq. (6), and regenerate all affected 32-qubit Path-A states, filter results, operating points and prose. The direct Galerkin spectra and the exact-slice results are not invalidated by this particular bug. The corrected overlap is evidence for the repair, not a replacement for rerunning the complete benchmark.

## 2. The large conditioned-coarse cost is substantially avoidable

**Locations:** Sections 5.1–5.3 and 7; supplement S5 and S7; `conditioned_term_cnots` in `src/pairinglib/hubbard.py`.

The cost model computes and uncomputes the same conjunction of antibonding controls separately for every coarse Pauli rotation. Every rotation in a coarse block has the same controls, and the coarse evolution does not change those controls. The conjunction can therefore be computed once into a clean flag, reused throughout the coarse block, and uncomputed once.

Within the manuscript's own relative-phase-Toffoli ladder convention, the quoted cost

\[
C_{\rm base}+n_{\rm rot}[6(c-1)+2]
\]

can be replaced for a coarse block by

\[
C_{\rm base}+2n_{\rm rot}+6(c-1).
\]

| Coarse block | Base CNOTs | Rotations | Controls | Quoted block cost | Shared-flag block cost |
|---|---:|---:|---:|---:|---:|
| Four coarse sites, eight fine sites | 32 | 24 | 8 | 1,088 | 122 |
| Eight coarse sites, sixteen fine sites | 72 | 52 | 16 | 4,856 | 266 |

The 4,856 value actually uses **52 rotations**, consistent with the module, whereas the main text and S7 say 48. That textual count should also be corrected.

The flag construction is valid for superpositions: its control register remains unchanged throughout the controlled coarse block. Relative phases in the compute/uncompute ladder cancel under the same assumptions as the original ladder construction. Fine evolution can change antibonding occupations, so the flag must be uncomputed before intervening fine evolution and recomputed for the next coarse block.

There is also a basis-conversion cost to make explicit. An embedded coarse block in the site register has the form

\[
G\left[|0_a\rangle\langle0_a|\otimes e^{-iX\tau}
+(I-|0_a\rangle\langle0_a|)\otimes I\right]G^\dagger,
\]

where G rotates between bonding/antibonding and site modes. Unless absorbed elsewhere in the complete circuit, G and its inverse cost 32 CNOTs per coarse application at eight fine sites and 64 at sixteen. The present scalar cost formulas do not visibly include these conversions around every application.

For orientation, keeping the manuscript's other counts and its 128-step construction, shared-flag counting gives 49,680 CNOTs for the eight-site Path-A preparation before those frame conversions, or about 57,872 if charged separately around both coarse blocks in each step. At sixteen sites, including the stated 18,544 coarse-preparation charge, the corresponding illustrative figures are 125,584 and 141,968, rather than 1,300,624. These are transparent constructive estimates; an explicit complete circuit must settle phase, ancilla and cancellation details before they become final table entries.

**Implication:** the claims that a conditioned coarse term is “never affordable” and that its factor-sixteen disadvantage cannot be moved by compilation are not defensible with the present implementation. This is a major opportunity for the proposed method, as well as a fairness issue. Report the clean work-qubit requirement: a system size of 16 or 32 qubits is not the complete device footprint of the chosen ladder circuit.

## 3. Double-excitation CNOT counts use the wrong Pauli weight

**Locations:** Section 3.5; Section 5; variational columns of Tables 2, 3, S4–S6; `uccsd_cnots`, around lines 699–716.

For four distinct sorted Jordan-Wigner indices p<q<r<s, a double-excitation Pauli word has support on the two intervals [p,q] and [r,s]. The Z factors cancel in the central gap. Its weight is

\[
w=(q-p+1)+(s-r+1),
\]

not the full span s−p+1 used in the cost function. For example, an excitation involving indices 0,2,8,10 has eight surviving words of weight 6; the current function charges weight 11. This disproportionately penalizes opposite-spin doubles in the block spin ordering.

Under the stated eight-strings-per-double and 2(w−1)-CNOT primitive, the full-template totals change as follows:

| Ansatz | Current CNOT count | Corrected parity-support count |
|---|---:|---:|
| Four-site UCCSD | 1,720 | 1,464 |
| Eight-site UCCSD | 30,864 | 23,312 |
| Sixteen-site UCCSD | 335,360 | 235,904 |

The variational state itself is unaffected by this counting error. Reprice every selected pruned and ADAPT operator list, not just the full templates, and regenerate the accepted-preparation costs.

A separate presentation correction concerns the orbital basis. The simulated state has the form Φ U_orb|reference occupation⟩, where Φ maps orbital coefficients to the site basis. A straightforward circuit applies the orbital excitation blocks before the final Φ Givens network. The prose saying site-register excitation blocks act after the determinant-preparing Givens network needs an explicit conjugation or a corrected order. A single final basis rotation can retain the intended count; it should be shown rather than assumed.

## 4. Longer ordinary ramps change the cost comparison materially

**Locations:** Sections 5.3, 6 and 7; Tables 2–3 and S5–S6.

The listed ramp inputs use T=8/t and 16/t at eight sites, and T=16/t and 32/t at sixteen sites. They are not the full cost-versus-accuracy frontier of ordinary adiabatic preparation. In particular, Section 6 already shows longer exact-slice ramp times but does not include their digital versions as competitors in Table 3.

Additional unfiltered ramps were evaluated with the same Hubbard group propagator, schedule, reference determinant and CNOT convention as the manuscript:

| Fine sites | Ramp time | Second-order steps | Unfiltered infidelity | CNOTs |
|---|---:|---:|---:|---:|
| 8 | 24/t | 192 | 1.53257×10^-4 | 27,760 |
| 8 | 24/t | 384 | 2.17583×10^-5 | 55,408 |
| 16 | 48/t | 384 | 6.87792×10^-4 | 117,216 |
| 16 | 64/t | 512 | 2.48805×10^-4 | 156,128 |
| 16 | 64/t | 1,024 | 1.24744×10^-4 | 311,776 |
| 16 | 128/t | 1,024 | 9.45124×10^-5 | 311,776 |

These are actual noiseless product-formula state fidelities against ED, with no filter or post-selection. They are not inferred from an adiabatic error estimate. They are feasible circuit constructions under the paper's assumptions, although the particular time/step choices are not claimed optimal. Some benefit from finite-step dynamics; at fixed achieved fidelity, a physical digital baseline should be allowed that benefit.

The consequences are substantial:

- Eight-site target 3×10^-4: 27,760 replaces a listed best ramp-plus-filter cost of about 40,000.
- Eight-site target 10^-4: 55,408 replaces about 91,000, and reaches appreciably higher fidelity.
- Sixteen-site target 10^-3: 117,216 replaces about 790,000.
- Sixteen-site target 3×10^-4: 156,128 replaces about 1,300,000.

Even restricting step counts to powers of two does not remove the last counterexample: T=64/t and 512 steps lie within the larger time/step scans discussed in the paper. All resource figures here still use the manuscript's unoptimized fine-step cost.

Thus “filter cycles dominate every route at thirty-two qubits” is a consequence of restricting the listed initial preparations, not an established property of the preparation problem. Similarly, a claimed tie between Path B and the ramp refers to those selected inputs. It does not show that Path B ties the best ramp cost at the target.

**Required revision:** minimize accepted-preparation cost jointly over T, step count, shift or Λ, filter step, width, and cycle count, with M=0 allowed. Plot the lower envelope for each method. The corrected shared-flag Path A and correctly priced variational states must enter the same comparison. A rigorous global optimum is unnecessary; a disclosed, adequately broad and comparably allocated parameter search is necessary.

## 5. Coarse preparation is charged but its actual error is not propagated

**Locations:** Table 3 caption; Section 6; supplement S7, p. 5; contrast the coarse-error discussion in Sections 3.2 and 4.5.

The 18,544-CNOT eight-site ramp is described as reaching 3×10^-4, yet its reported overlap is 0.999659468, giving infidelity 3.40532×10^-4. It fails that threshold. The statement that this residual is below all sixteen-site targets is also false: the tightest target is 3×10^-4.

More fundamentally, the simulated sixteen-site prolonged and refined inputs use the exact coarse eigenstate while paying for this approximate ramp state. Actual coarse-state propagation gives:

| Sixteen-site input/result | Exact coarse state | Charged approximate coarse preparation |
|---|---:|---:|
| Prolonged overlap | 0.869188994 | 0.868470417 |
| Published wrong-constant Path A, T=16/t, 128 steps | 0.998346666 | 0.997784182 |
| Corrected Path A, T=16/t, 128 steps | 0.999298441 | 0.998892827 |
| Path B, Λ=8t, T=16/t, 256 steps | about 0.99476 | 0.994138045 |

For corrected Path A, coarse-state error changes whether the unfiltered input meets 10^-3. This is precisely the mechanism the manuscript itself explains elsewhere. It cannot be neglected just because the coarse error looks small in isolation. Subsequent filtering can remove part of it, but that must be calculated for the actual state.

At eight sites, Table 2 charges the prolonged state only 16 CNOTs and omits preparation of the four-site coarse eigenstate. State explicitly whether that comparison grants an exact coarse-state oracle; preferably include a realizable coarse preparation consistently at both sizes.

## 6. What is sound in the physics and mathematics

The following central claims survived algebraic and numerical checks:

1. **Galerkin coarse Hamiltonian.** Equations (4) and (5) are correct for the stated bonding block map: chain parameters t/2,U/2 with constant −tN; square parameters t/2,U/4 with constant −2tN. Projection identities were checked to machine precision.
2. **Symmetry preservation and the plaquette obstruction.** Reflection/C4v intertwining and the total-spin diagnostic are sound. At U=4t the coarse square ground state is a B1 singlet, E=−9.414213562, while the fine square ground state is an A2 triplet, E=−10.419659046. The prolonged wrong-symmetry overlap was reproduced at approximately 7×10^-28. A symmetry-preserving filter cannot repair exact zero support on the target.
3. **Chain spectral interpretation.** At eight sites and U=4t, the first gaps 0.223489t and 0.459157t are triplets, while the accessible even-singlet gap is 0.587888t. The quoted p0=0.797409 is reproduced.
4. **Ideal rodeo formula.** Equation (8) correctly describes the fidelity of the acceptance-weighted ensemble for independent Gaussian times, exact controlled evolution and on-resonance objective energy. It does not describe every individual random schedule.
5. **Fixed-step digital filtering.** Quantizing every cycle to integer powers of one step unitary gives common eigenvectors across cycles. On resonance with an isolated eigenphase, the accepted state tends to that digital eigenvector and acceptance tends to its input weight. At eight sites, the reconstructed effective energy shift is +0.013519712t and the limiting infidelity is 5.76435×10^-5, agreeing with the paper. This is a well-founded part of the study.
6. **Accepted-cost normalization.** The published accepted-preparation estimator includes preparation retries, early termination and schedule-specific survival weights. It is a much better comparison metric than multiplying a total cycle cost by a bare 1/p without accounting for early failures.
7. **JRS bound.** The endpoint terms vanish for the implemented sin² schedule, and ε^-1/2 is appropriate for an infidelity guarantee. Independent symmetry-restricted calculations reproduce Table S1's arithmetic. The table remains a numerical quadrature evaluation of a rigorous bound, rather than a certified enclosure of all discretization errors. See [Jansen, Ruskai and Seiler, Theorem 3](https://arxiv.org/pdf/quant-ph/0603175).

## 7. Mathematical interpretations that should be tightened

### Response metric versus dynamical leakage

Main Section 4.2 and S1 define

\[
\chi_u=\sum_{n>0}\frac{|\langle n|\partial_uH|0\rangle|^2}{(E_n-E_0)^2}.
\]

This is the ground-state quantum metric, or squared eigenvector response. Its integral alone does not set leading adiabatic leakage. First-order transition amplitudes also contain u'(s) and oscillatory phases involving T times the integrated gaps. Rescaling H by a constant leaves χu unchanged but changes the dynamics at fixed T. Moreover u'(0)=0, so a maximum of χu at the endpoint is not a maximum of the actual eigenvector speed under this schedule.

The reported initial values 0.375, 0.527083, 0.714586 and 0.874011 are correct. Keep them as useful sensitivity diagnostics; remove the claim that their integral directly determines leakage or explains runtime by itself.

### Size-independent time and an order-t protected gap

The statements “flat in size” and a persistent order-t gap exceed the evidence. The reported threshold time changes from 10/t for 4→8 to a first tested successful value 16/t for 8→16, while the accessible gap shrinks from 0.588t to 0.138t.

At U=0 and fixed two particles per spin, an accessible even-singlet excitation has energy

\[
2t\left[\cos\frac{\pi}{L+1}-\cos\frac{3\pi}{L+1}\right]
\sim \frac{8\pi^2t}{(L+1)^2}.
\]

It reproduces the eight-site value 0.879385t and closes with size. Symmetry excludes forbidden states; it does not supply an asymptotically nonzero accessible gap. Present the timing as a favorable observation over the tested cases, not a scaling law.

### Transport gap versus input-support gap

Use distinct definitions for the minimum gap along the path within its invariant symmetry sector and the smallest final excitation carrying weight in a particular filter input. They can coincide in the examples but are not identical concepts. Spin-resolved UCCSD is not automatically SU(2)- or parity-adapted. Check each input's actual support before asserting the same protected gap for all inputs. An eight-site check found only approximately 10^-7 weight below the chosen gap in the full/pruned UCCSD states, so there is no demonstrated significant bias at those displayed targets; the qualification remains necessary at other sizes and under noise.

### Continuum interpretation

Section 3.2 says a continuum-scaling coarse model would “change only p0.” Equation (6) shows otherwise: changing the coarse Hamiltonian changes the full path, shift convention, gaps and finite-time dynamics. If only the input is changed while holding the Galerkin path fixed, it need not be an initial eigenstate.

The fixed-U/t hierarchy is a valid algebraic lattice hierarchy. It is not, without additional parameter/spacing analysis, a continuum convergence result for one fixed physical model. For a conventional one-dimensional contact discretization, t scales as a^-2 and U as a^-1 at fixed mass and coupling. These are different requirements from preserving U/t in Eq. (4). No new continuum study is essential if the scope is stated accurately.

### Smaller mathematical corrections

- In S1 the displayed definition of ||B−A|| omits −μ(I−PP†). The notebook uses the correct shifted expression, so this is a textual correction rather than evidence of wrong Table S1 values.
- For a rank-one tracked branch, projector-distance squared equals the state's infidelity. “Every state in the branch” versus “one initial state” does not explain the loose bound here. Conservative norm estimates do; the table also compares theoretical ε=10^-2 against empirical 10^-3.
- Path B approaches Path A's prolonged initial eigenstate as Λ grows; at finite Λ their initial Hamiltonians and exact initial eigenstates differ. Avoid introducing them as literally identical endpoints.
- One Givens rotation per block/spin applies to two-site chain blocks. A uniform four-site plaquette mode generally needs three rotations per block/spin.

## 8. Statistical precision and reproducibility

Forty schedules at eight sites and twenty at sixteen give useful exploratory estimates, and reporting bootstrap errors is commendable. However, hard pass/fail targets and minimum cycle counts can change within those uncertainties.

There is a concrete example. Table S5 says the eight-site T=8/t ramp reaches 10^-4 after six cycles. The forty-schedule estimate is 8.92595×10^-5. Deterministic evaluation of the same quantized-time digital ensemble gives **1.05045×10^-4**, failing that target. This does not reverse the ranking among the listed fixed inputs, but demonstrates why the operating-point decisions need uncertainty treatment.

At eight sites the exact digital ensemble is inexpensive once the step is diagonalized. If Uδ|uj⟩=e^{iφj}|uj⟩, let fj(n)=[1+exp(in(φj+Eobjδt))]/2 and Qij=E_n[fi(n)fj(n)*], using the actual rounded Gaussian distribution. Then

\[
\rho_{ij}^{(M)}=c_i c_j^*Q_{ij}^M.
\]

Compute fidelity with the exact ground state from this full density matrix, retaining coherences because the exact and digital eigenbases differ. This removes schedule sampling at eight sites. At sixteen sites use more schedules or deterministic channel methods where practical, and bootstrap the paired comparison and operating-point selection itself. Shared random times are helpful; separate inputwise relative errors do not quantify uncertainty in the cost-minimizing cycle count.

For the iterated sixteen-site digital eigenvector, report a final eigenvector residual ||Uδψ−exp(−iEeffδt)ψ||. Agreement of the method with the eight-site Schur result is useful validation, but a residual establishes convergence in the larger calculation directly.

Make the notebook and script versions agree. Some surrounding explanatory text and an older diagnostic script use superseded parameter/cost conventions. One notebook field labeled cost per attempt omits early-termination weights, although the accepted-cost field used by the paper is correct. Archive the exact executed source, parameter record and outputs used to generate the revised figures.

## 9. Fairness and scope of the comparator set

The paper makes a serious effort at fair comparison: common targets, common filter, explicit post-selection cost, declared exclusion of setup/training measurements, and an ordinary ramp baseline. The comparison is appropriately candid about unfavorable refinement results. Those are strengths worth retaining.

Nevertheless, the current methods are priced with different avoidable inefficiencies, and the selected operating points do not establish a full frontier. In addition to Sections 2–4 above:

- State that the central metric is amortized logical CNOT count per accepted state, not total experimental work or circuit depth. Give ancilla counts and connectivity assumptions.
- Preserve the disclosed distinction between optimizer-free execution and ED-assisted calibration. The benchmark does not yet establish an end-to-end optimizer-free, calibration-free hardware protocol.
- Include a Hamiltonian variational ansatz as a natural Hubbard comparator, or restrict claims to the particular UCCSD/ADAPT pools tested. Such ansätze are established primary-literature baselines: [Wecker, Hastings and Troyer](https://arxiv.org/abs/1507.08969), [Cade et al.](https://arxiv.org/abs/1912.06007). Their existence does not prove they beat refinement here; it makes a conclusion about variational methods generally too broad without testing them.
- The known fixed particle number makes the sum of the interaction's single-Z terms a scalar. Absorbing that scalar into the objective phase reduces the stated controlled second-order cost from 352 to 288 at eight sites and 736 to 608 at sixteen, under the same primitive. This is another useful common compilation improvement.
- The demonstration doubles site count while keeping four particles. It probes a progressively dilute family, not scalability at fixed density or the half-filled strongly correlated Hubbard regime. Either retain this explicit scope or add at least one fixed-filling case.
- A 32-qubit register with a 14,400-dimensional sector is a useful circuit-development test, but is far from a claim of classical intractability. State clearly that these are noiseless algorithmic benchmarks. Hardware accessibility cannot be inferred from register size alone at 10^4–10^6 logical CNOTs.

## 10. Direct contradictions and submission cleanup

1. **“Single-layer ansatz has lost two orders of magnitude.”** The stated infidelities are 3.143×10^-3 at eight sites and about 6.98×10^-3 at sixteen, a deterioration of about 2.2, not 100. Correct the abstract, introduction and discussion, or specify a different comparison supported by data.
2. **“Ramp is cheapest at every target.”** Table 3 itself gives ADAPT at 68,008 CNOTs for 10^-2, below the ramp's 78,304. Restrict even the current-table claim to the tighter targets. After repricing and extending the search, rewrite the ranking entirely.
3. **Forty-eight versus fifty-two coarse rotations.** Correct the sixteen-site description to match the actual block count.
4. **Scope of fine-plaquette spin claim.** “Every U>0” is stronger than a study over the stated U/t=1–8 values. Restrict to the tested range unless proved separately.
5. **Availability references.** Replace the repository URL and Zenodo DOI placeholders and complete the authorship/reference information for the Hubbard ED notebook. A publisher-inserted supplement URL placeholder can follow journal convention; missing repository identity should not survive submission.
6. **Presentation.** Add uncertainty bars or bands to the resource comparison plots where schedule uncertainty matters. Fix the step-convergence explanation rather than merely deleting its warning. Keep main/supplement definitions and rounded tolerances consistent.

## 11. Recommended revision order

1. Repair the sixteen-site embedded Hamiltonian and add decomposition checks outside the coincidental Nparticles=Ncoarse cases.
2. Write explicit complete circuits for the shared-control coarse block and orbital-basis variational preparation. Correct JW double weights, frame changes, ancilla counts and phase handling.
3. Propagate the actual charged coarse-preparation state at each hierarchy level.
4. Recompute cost frontiers, allowing longer ramps and M=0. Compare methods at the same achieved fidelity with a comparable search budget.
5. Remove schedule noise where exact digital averaging is available; attach uncertainty to remaining operating-point decisions.
6. Correct the response-metric, gap and scaling interpretations, and replace the abstract/discussion with conclusions supported by the rebuilt tables.
7. Add a Hubbard-specific variational comparator and a fixed-filling test if broad competitiveness/scalability claims are retained. These are extensions rather than prerequisites for a carefully scoped small-system benchmark.

The defensible scientific contribution is a transparent study of how a concrete Galerkin hierarchy, symmetry selection, transport construction and digital filtering interact. The present evidence supports that contribution. It does not yet settle the economic advantage or disadvantage of refinement, because both the preparation frontiers and compilation costs require the corrections above.

## Audit scope and numerical evidence

The review included visual inspection of the relevant equations, figures and tables in both PDFs, algebraic checks, code inspection and targeted numerical reproductions. It did not rerun the entire variational optimization campaign, synthesize every complete circuit in a hardware transpiler, or perform noisy hardware simulations. Reproductions used the supplied model builder with independently constructed checks; they are not a second independently written Hubbard ED implementation. The numerical counterexamples are therefore strong internal consistency and benchmark tests, with that limit on independence stated explicitly.

The main and supplemental PDFs are located at:

- `/Users/mhjensen/Software/QuantumComputingQuNorth/EigenQ/pairing-eigenq/paper3/eigenq_hubbard.pdf`
- `/Users/mhjensen/Software/QuantumComputingQuNorth/EigenQ/pairing-eigenq/paper3/eigenq_hubbard_supp.pdf`

Scratch evidence is retained in the current chat workspace under `work/hubbard-review/`, including `path_a_reproduction.json`, `long_ramps.json`, `prepared_coarse.json`, `exactdigital_result.json`, and the corresponding check scripts. The original article, supplement, notebook and code remain unchanged.
