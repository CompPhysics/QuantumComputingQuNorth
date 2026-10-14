# Integrated assessment of the pairing manuscript

Reviewed: October 14, 2026. Manuscript dated October 12, 2026, “Resolution refinement followed by the rodeo algorithm for eigenstate preparation: implementation and assessment on the pairing Hamiltonian.”

This document combines the assessments of mathematical and physical correctness, the hybrid pipeline, the interpretation of the results, and the fairness of the circuit-resource comparisons.

## Overall assessment

The central construction is physically sensible and the main numerical comparison in Section V is reproducible. The manuscript's negative conclusion is supported: at the tested pairing-model sizes, refinement is not a cheaper source of a high-overlap input than a small variational circuit. The work is useful as an implementation and controlled benchmark.

Section V is substantially fair as a comparison of four specified inputs under one specified, relatively simple compilation and filter. It is not a comparison of the best available implementations of the competing methods, and it is not a hardware circuit-depth study. Several statements should be corrected before publication, particularly the inference of inverse-gap scaling, the claimed common accuracy after six cycles, the treatment of the UCCSD preparation cost, and the distinction between practical calibration and an optimizer-free circuit.

Scope: I inspected the PDF and the accompanying source/notebook and selected implementation files, checked primary literature, and independently reconstructed the six-dimensional seniority-zero problem for k=4, N=4, g=1. I reproduced the Section V filter and refinement numbers without using the package's state-preparation routines. The audit also checks time-slice convergence, coarse-input error propagation, the shift-dependent response norm, the Section IV filtering statistics, and the long-run digital filtering map. I did not rerun every figure, audit the companion manuscript, or perform a hardware-noise simulation.

## 1. Physical model and common accuracy measures

**Pairing Hamiltonian and reference energy, Section II.** With the stated convention -g/2 multiplying the pair-scattering sum, the filled two-level determinant has E_HF=2-g. At g=1, independent diagonalization gives E0=0.6355484735756 and squared overlap of that determinant with the ground state p=0.8676375677137. These agree with the manuscript.

**Energy versus fidelity.** The inequality <H>-E0 >= Delta(1-F) is correct. Comparing energy error with infidelity would be invalid; the present comparison uses infidelity and explicitly retracts that earlier argument. The use of occupation numbers as an additional diagnostic is valuable: a general observable can have an error linear in the state-vector error, whereas the energy error about an eigenstate is quadratic.

**Seniority and scaling.** “None of the stages uses seniority” should mean that the register and circuit implementations do not exploit seniority for compression. The favorable dynamics and available small input circuits do rely on the pairing Hamiltonian's seniority conservation. Fixed-N=4 growth in orbital count is a basis-size study and cannot establish many-particle complexity at fixed filling, as the manuscript elsewhere correctly states.

## 2. Section III: correctness of the hybrid pipeline

### 2.1 Coarse-state preparation, excitation circuits, and input error

The ordered product of unitary excitation factors is a valid normalized, particle-number-preserving ansatz. The Pauli factors within each of the stated single/double excitation generators commute, so their circuit decomposition does not introduce an approximation within one excitation. The finite ordered product differs, in general, from the exponential of the summed cluster operator. It is therefore best called a factorized or disentangled UCCSD ansatz; specify the ordering as part of its definition.

The inspected analytic-gradient implementation differentiates that ordered product correctly. A generic fermionic excitation generator has eigenvalues 0 and +/-i, so the parameter dependence of an expectation value can contain both first and second harmonics. A simple single-frequency two-shift rule is not generally sufficient. This is a sound reason to use the analytic derivative here; “beyond the first factor” describes a common way those extra harmonics appear, rather than a universal theorem about every later factor.

The zero-single-amplitude observation is consistent with the symmetry and the zero-initialized optimization used in the calculations. Phrase it as a result of this optimization on the seniority-preserving parameter manifold. Symmetry alone does not prove that every unconstrained optimum of an ordered ansatz has zero seniority-changing amplitudes. Also, unitary doubles and classical CCD are distinct approximations; their common absence of singles does not make their equations or variational properties equivalent. Section IV's different UCCSD/CCD energies are therefore physically sensible.

The substantive missing assumption is that VQE ordinarily supplies an approximation to the coarse ground state, whereas Section III subsequently reasons as if its input were an exact eigenstate. Write

\[
|\psi_{\rm low}\rangle=c_0|E_0^{\rm low}\rangle+\sum_{j>0}c_j|E_j^{\rm low}\rangle,
\qquad \epsilon_{\rm low}=1-|c_0|^2.
\]

If the ground branch remains isolated, ideal adiabatic continuation transports its weight. In the perfectly adiabatic limit it gives fine-space ground-state fidelity |c0|², rather than automatically giving unity. Refinement is unitary and does not purify a genuinely imperfect coarse input. The subsequent rodeo stage can remove the residual excited-state component if the desired overlap remains nonzero.

This distinction does not undermine the reported main benchmarks: the k=2 start is exactly one determinant, and the k=3 coarse preparation is essentially exact. It does mean these calculations do not yet test the pipeline's robustness to a non-negligible coarse-solver error.

A useful combined guarantee is

\[
1-F_{\rm refine,digital}
\leq\left(\sqrt{\epsilon_{\rm low}}+a_{\rm ad}+\eta\right)^2,
\]

where a_ad bounds leakage amplitude for an exact coarse-ground-state start and eta bounds the digitization error in operator norm. This follows by decomposing the input into its ground and orthogonal components and applying the triangle inequality to the final component orthogonal to the target. Eq. (5) supplies a_ad <= I/T for the stated ideal path, with I its integral. The bound is conservative but makes every source of preparation error explicit. For an exact coarse input and a target total infidelity epsilon, one sufficient allocation is a_ad <= sqrt(epsilon)/2 and eta <= sqrt(epsilon)/2. Assigning the full epsilon to ideal adiabatic infidelity before adding eta=sqrt(epsilon)/2 does not certify the same total epsilon. Once the actual output overlap p is known, Eq. (9) can use that p directly for an exact rodeo filter.

### 2.2 Prolongation and the energy identity

The correct general identity is

\[
\langle P\psi_{\rm low}|H_{\rm high}|P\psi_{\rm low}\rangle
=\langle\psi_{\rm low}|H_{\rm low}|\psi_{\rm low}\rangle
\quad\text{if } P^\dagger H_{\rm high}P=H_{\rm low}.
\]

This operator identity holds for the nested bare pairing Hamiltonians used here, with unchanged spacing and coupling. It does not follow from the definition of prolongation alone. It can fail for a renormalized effective coarse Hamiltonian, changed interaction parameters, or a different discretization. Equality with an exact E_low additionally requires an exact coarse eigenstate; with a variational state it equals the variational coarse energy.

Nor does preservation of the coarse expectation value establish a high fine-space overlap. That is an empirical property of the chosen coarse/fine Hamiltonians, or a consequence of additional energy-gap information. The reported overlaps support it in this example. State the identity and the observed quality of the guess as separate claims.

### 2.3 The shifted path, adiabatic guarantee, and digitization budget

The path and the conditioned coarse term in Eqs. (3)-(4) are correct. The initial ground-state gap is correctly given by the minimum of the coarse excitation gap and mu-E0_low, with the latter alone relevant for a one-dimensional coarse space.

Taking mu>E0_low establishes the correct initial ordering below the complementary zero-energy sector. It does not by itself rule out an interior crossing or prove the desired continuation. The manuscript should explicitly condition the final-ground-state statement on the ground branch remaining isolated throughout the path. Section IV's spectral calculations are the numerical check of that condition for the tested examples.

If the shift is set as mu=E_var+b, the initial gap is

\[
\Delta(0)=\min\{\Delta_{\rm low},\ b+E_{\rm var}-E_0^{\rm low}\}.
\]

It is simply b in the exact one-dimensional example, but not automatically for every variational coarse state. The shift is not merely a global energy change along the interpolation: rewriting its contribution gives a global -mu I plus a relative penalty mu cos²(pi s/2)(I-PP†). That relative penalty is why changing mu changes the path dynamics while leaving the final physical spectrum unchanged.

Eq. (5) is a valid ideal-continuous-evolution bound when the initial state is the exact coarse eigenstate and the gap assumptions hold. The sine-squared schedule removes its endpoint first-derivative terms. Combining that guarantee with coarse-input and digitization errors requires the explicit budget above, rather than assigning the complete target error separately to multiple stages.

There is also a definition mismatch in the first-order Trotter condition. The cited nucleation theorem is stated for two groups. For the many groups used here, a sufficient extension uses

\[
\Lambda_2=\max_s\sum_{i<j}\|[H_i(s),H_j(s)]\|.
\]

Section IV/Table II uses this sum correctly. Section III's “largest commutator norm between the groups” is ambiguous and is insufficient if read as the largest individual pair. Give the sum explicitly and describe this as an extension of the two-group bound. The midpoint-freezing contribution T Lambda1/(4 n_s) is consistent with the usual Duhamel bound.

“Optimizer-free” accurately describes execution of the refinement circuit after its controls are chosen. The reported T* is determined from exact fine-space fidelity, and the small successful product-formula step is chosen by simulation. Therefore, the claim that the schedule is set using coarse-space information alone is not demonstrated. Report this as an exact-benchmark calibration and discuss an implementable stopping/calibration protocol. For a nontrivial coarse space, coarse-state training is also a setup cost.

### 2.4 The rodeo filter, random ensembles, and resolution

The Kraus operator K0(t)=[I+exp(iEt)exp(-iHt)]/2, Gaussian attenuation q_j=[1+exp(-sigma²(E_j-E)²/2)]/2, and acceptance-weighted accepted-ensemble fidelity are correct under exact evolution and independent randomized cycle times. For a nondegenerate target on resonance, the target weight remains p and the total acceptance approaches p. The bound using q_max has the correct inequality direction.

For an exactly degenerate energy, filtering selects the corresponding eigenspace, not a unique vector. A component with zero initial target overlap cannot be created by an exact function of H.

The distinction between an unweighted mean of normalized post-selected fidelities and the fidelity of the accepted ensemble is important and correctly made. Specify that schedules are sampled per attempted preparation if this accepted-ensemble operational interpretation is intended. Drawing one schedule and repeating until it succeeds is a different sampling procedure.

The sentence about the 2^-M background needs correction. Since q_j>=1/2, the Gaussian-averaged retained fraction q_j^M is at least 2^-M. A particular fixed random schedule has no such lower bound, because a cosine factor can vanish. Distinguish the ensemble statement from an individual schedule.

At fixed M it is reasonable to describe the energy width as proportional to 1/sigma. If M also varies, near a peak

\[
q_j(E)^M\simeq\exp[-M\sigma^2(E-E_j)^2/4],
\]

so the width scales as 1/(sigma sqrt(M)). Likewise, the statement that the target weight p is preserved and acceptance tends to p is an on-resonance, exact-evolution statement. At detuning d=E-E0, its averaged weight is p q0(E)^M and decays with M. Conditional fidelity can still improve if the desired eigenvalue is closest to E, but detuning toward another eigenvalue can select a different state. Maintaining substantial target acceptance requires approximately |d| sigma sqrt(M) << 1.

Finally, the quoted O(log(1/delta)/(p epsilon)) estimate concerns ideal controlled-evolution effort at a resolution sufficient to isolate the target. It is not the complete CNOT complexity of the implemented pipeline. Input preparation, target-energy search, Hamiltonian-simulation accuracy, and any reduction of the product-formula step size with requested precision remain to be counted. Section V's accuracy floor is direct evidence of why that distinction matters.

## 3. Section IV: interpretation and consistency of the results

### 3.1 Coarse-solver comparison

The unitary ansatz is variational, so its energies remaining above FCI are expected. Classical CCD can overbind because its energy is not a variational upper bound. Agreement of the gate and operator implementations supports the correctness of the chosen ansatz compilation. It does not establish that the optimizer found the global optimum over every available parameter or ordering.

The reported HEA difficulty is evidence about the tested ansatz, starts, optimizer, and budget. It is not by itself a diagnosis of barren plateaus or a general failure of hardware-efficient methods. Section V's conclusion should rest on the measured structured comparator and resource counts. The manuscript's explicit distinction between basis-size growth at fixed N and many-particle scaling is correct and should be retained.

### 3.2 Refinement fidelity and numerical convergence

The “exact-slice” calculation still approximates a continuously changing Hamiltonian by midpoint slices. Independent convergence checks show this is a small numerical error at the reported operating points:

| Total time | Slices | Ground-state population | Physical energy |
|---|---:|---:|---:|
| T=25 | 160 | 0.9975712826 | 0.6432272773 |
| T=25 | 1,280 | 0.9975698271 | 0.6432320131 |
| T=30 | 160 | 0.9987571886 | 0.6396536865 |
| T=30 | 1,280 | 0.9987561057 | 0.6396573241 |

The changes in population are about 1-1.5e-6, much smaller than the approximately 1e-3 adiabatic infidelities. Thus the principal ideal-refinement numbers are well resolved. They are not the T=30, 80-step product-formula state in Section V, whose population is 0.9977209. The difference is expected and is explicitly quantified there.

Energy approaching FCI from above is consistent with the variational principle for any normalized prepared state. Observed improvement with time is plausible; monotonicity at every time is not guaranteed by the adiabatic theorem and should not be implied beyond the plotted data. An open path gap establishes applicability of adiabatic reasoning, not the particular inverse-gap cost law.

### 3.3 The nontrivial coarse stage

The k=3 coarse state exercises a real preparation stage and gives the reported prolonged overlap of about 0.963668. Its essential exactness means it does not test an imperfect coarse solver. The shorter sweep and the additional conditioned-coarse-Hamiltonian cost are mutually consistent with Sections III and V: better initial overlap need not mean lower total gate cost.

The stated threshold T*=9 is grid-dependent. An independent ideal calculation already reaches population 0.990737 at T=8. A symmetric product formula splitting both coarse and fine Hamiltonians into their actual commuting blocks gives population 0.991172 at T=8 with 22 steps. Under the manuscript's unoptimized counting, this is approximately 28,696 CNOTs including the full-template coarse preparation, versus 26,660 for the reported 2-to-4, T=16 operating point. These support comparable costs; they are not globally optimized minima or precise universal ratios.

There is an implementation distinction to make explicit. The k=3 embedded coarse Hamiltonian contains noncommuting terms. Exponentiating the entire embedded A as one exact block, as the generic numerical routine does, is not identical to controlling the 33 rotations of one coarse first-order product-formula step. A complete circuit validation must include that coarse splitting. My check of the flattened symmetric formula still reaches the stated target, so this is a missing validation step rather than evidence that the chain fails.

The 272-CNOT coarse-state template can also be pruned: the three-dimensional seniority-zero coarse space admits exact preparation of any real normalized vector with two occupied-to-virtual pair rotations, costing 96 CNOTs in the stated construction. This reinforces the relevance of a prolonged-coarse-state-plus-rodeo baseline.

### 3.4 Coupling and shift scans: what they establish about scaling

The coupling scan does establish that the nearly fixed minimum gap does not determine the required time by itself. Its longer times correlate with loss of prolonged overlap. But changing g also changes transition matrix elements and eigenvector-path geometry. It is not an experiment isolating overlap as the unique cause. In fact, the shift scan changes the required time while keeping the initial overlap unchanged.

Replace “grows because the starting overlap degrades” with “grows alongside the loss of prolonged overlap, consistent with a larger change of the eigenstate along the path.” Quantify that change using response norms or path length if a causal interpretation is desired. Similarly, a higher coarse start at strong coupling is a reasonable candidate, but its superiority requires a comparison that includes the cost of preparing it and evolving its conditioned Hamiltonian. Section IV C itself explains why shorter sweep time does not establish that superiority.

**Shift-scan fit and its interpretation.**

Section IV D, Fig. 6, and the conclusion claim that the observed timing is incompatible with an inverse-square gap law. The evidence does not support this claim.

The four finite timing points are approximately (Delta_min,T*)=(0.3,50), (0.6,16), (1.0,8), (1.499,6). The plotted linear fit includes an intercept:

T* approximately 17.05/Delta_min - 8.43.

The prose drops the intercept and writes approximately 17/Delta_min, which does not reproduce its own operating point T*=16 at Delta_min=0.6. Moreover, fitting the same rounded points to a/Delta_min²+b gives approximately

T* approximately 4.13/Delta_min² + 4.16.

The sums of squared time residuals are about 28.2 for the affine inverse-gap fit and 0.22 for the affine inverse-square fit, with two fit parameters in each case. This does not establish inverse-square scaling either: four threshold times from a coarse time grid, with the path changing as mu changes, cannot identify an asymptotic exponent. It does directly contradict the assertion that inverse-square behavior is excluded.

Replace “confirms” and “incompatible” with an empirical sensitivity statement. Compute the transition susceptibility chi(s)=sum_{n>0}|<n|partial_s H|0>|²/(E_n-E0)² and its behavior across the path family if a scaling argument is desired. A finite value for every finite system is not sufficient; its dependence on the gap and system size must be controlled. The cited resolution-refinement derivation uses a linear schedule, whereas this paper uses sin²; transfer of an asymptotic scaling claim should state the relevant regularity and uniformity assumptions.

The physical fixed-N gap of about 1.823 is a valid conservative gap. For this pure pairing calculation, the dynamically accessible seniority-zero excitation gap at the endpoint is about 2.29983. Distinguish an all-sector bound from the transition-relevant gap. The artificial shift-controlled initial gap is a third quantity and should remain separately identified.

**The response factor also changes with the shift.**

An analytic check of the hypothesis invoked in Section III reinforces this interpretation. Let b=mu-E_low for the k=2 to 4, N=4 path and use its linear geometric coordinate u=sin²(pi s/2). Then

\[
H(u)=(1-u)A+uB,\quad
A=-b|\mathrm{HF}\rangle\langle\mathrm{HF}|,\quad
B=H_{\rm high}-\mu I.
\]

At u=0 there are four orthogonal pair-transfer configurations coupled to HF, each with amplitude -g/2. Every orthogonal state of A has energy zero and hence excitation gap b. Consequently the first-order eigenvector-response norm used in the cited linear-path argument is exactly

\[
\chi_u(0)=\sum_{n>0}
\frac{|\langle n|\partial_uH|\mathrm{HF}\rangle|^2}{(E_n-E_0)^2}
=\frac{g^2}{b^2}.
\]

For g=1, changing b from 0.3 to 1.5 changes this quantity from 11.11 to 0.4444, a factor of 25. An independent calculation along the six-dimensional eigenvector path gives those same values and finds the maximum at the initial endpoint for the tested shifts.

Thus the response factor is not held uniform in the gap parameter varied in Section IV. The sine-squared schedule makes partial_s H vanish at s=0, but it traverses the same eigenvector path; that endpoint cancellation is not a proof of the missing uniform scaling assumption. This does not establish a different asymptotic exponent. It shows why the cited conditional inverse-gap argument cannot be tested by attributing these timing changes to the gap alone.

### 3.5 Sufficient bounds versus empirical requirements

The numerical bound is a valid sufficient guarantee under its assumptions, and a large separation from empirical requirements is credible. A loose sufficient bound is not a lower bound on the resources required by another method.

Table II compares an operator-norm Trotter guarantee with an empirical criterion involving one final overlap. These are different error criteria. Reproducing the scalar ground-state overlap does not certify that the digitized evolution approximates the ideal unitary in operator norm, or even that the two output vectors are close. Retain the qualitative observation that rigorous bounds are conservative, but label the reported overshoot factors with their different criteria or add an empirical check at the same norm tolerance.

The explanation that the highest states “never couple” to the tracked state is literally incorrect for this example. Restricting entirely to the accessible six-dimensional seniority-zero space already gives ||B-A||=7.6205693, matching the quoted full-sector norm. At s=0.25, derivative couplings from the ground state to the upper two accessible pair-sector eigenstates are about 0.11152 and 0.08148, both nonzero.

A supported explanation is that global spectral norms ignore the small, state-specific transition matrix elements and their individual gaps and may therefore be conservative. The magnitude of the bound alone does not prove a unique cause of its looseness or verify the inverse-gap mechanism.

### 3.6 Rodeo convergence, statistics, and energy search

The ideal filtering numbers can be reproduced with the stated refined input and an appropriate specified random seed. They do not agree to better than 1% with the exact Gaussian-ensemble infidelity:

| Cycles M | Sampled unweighted infidelity | Sampled acceptance-weighted infidelity | Exact Gaussian-ensemble infidelity |
|---:|---:|---:|---:|
| 1 | 1.34977e-3 | 1.35017e-3 | 1.21584e-3 |
| 2 | 6.48861e-4 | 6.49099e-4 | 6.08287e-4 |
| 6 | 2.08501e-5 | 2.08521e-5 | 3.80397e-5 |
| 10 | 1.98510e-6 | 1.98513e-6 | 2.37756e-6 |

Here p=0.9975712826, sigma=3, and the sample uses 40 schedules with seed 7 reset for the refined input, reproducing the article's numbers. The six-cycle weighted infidelity is 45.2% below the analytic expectation. Comparing fidelities close to unity can conceal that discrepancy; the relevant comparison is the infidelity actually plotted and discussed.

Changing weighting does not resolve it. At this p, the relative effect of weighting on the sampled infidelity is bounded by 1/p-1, about 0.244%. The discrepancy is ordinary finite sampling of products of random attenuation factors. Independent Gaussian-moment calculations give approximate relative standard errors of the mean residual weight of 12.7%, 36.6%, and 85.5% at M=2, 6, and 10 for 40 schedules. Therefore “±30% statistical precision” is not justified at the smallest errors. Use analytic ensemble curves plus uncertainty intervals for the finite-schedule points.

There is also seed-policy drift: the main notebook advances its generator through the determinant input before the refined input and reports about 5.45e-6 at ten cycles; resetting the seed for the refined input gives the article's 1.99e-6. Both are plausible samples. Select and document one reproducible convention for text, figures, and notebook. For the actual p, Eq. (9) requires 12 cycles to reach ensemble infidelity at most 1e-6. “Of order ten” is reasonable; exactly ten does not reach that target in the ensemble.

The fine energy scan is a valid local filter-response illustration, but its grid is centered on the known exact E0 in the notebook. Consequently the 0.6352 interpolation does not independently demonstrate a device-ready coarse-to-fine search or establish four-decimal precision. To demonstrate the search, choose the fine interval from the coarse result alone and include probability-estimation uncertainty. The small numerical target-energy error is fully consistent with the negligible effect on the short exact-evolution filter in the reported detuning test. The broader statement that acceptance diagnoses detuning should remain limited to the tested regime; it is not a certificate that the selected state is the intended one.

## 4. Section V: correctness and fairness of the resource comparison

### 4.1 Circuit primitives and the scope of block simulation

A weight-w Pauli rotation has a standard 2(w-1)-CNOT parity-ladder construction. The cancellation of internal Jordan-Wigner strings in pair transfers is correct, as is the stated 48-CNOT construction from eight commuting weight-four rotations. Controlling only the central rotation is valid when the surrounding circuit is uncomputed. The uncontrolled first-order cost 24k(k-1)+2k gives 296 at k=4; adding controls to 60 rotations gives 416. The 14-CNOT diagonal synthesis of the four-appended-qubit conditional phase is consistent with the stated construction.

The block simulation is mathematically legitimate: each complete block preserves the relevant subspace and equals the corresponding circuit unitary there. It should continue to be described as block simulation, rather than a full gate-by-gate execution of all long circuits.

### 4.2 The pruned variational comparator and ADAPT claims

The 1,312-CNOT count charges for the full 26-operator template even though this model admits a much smaller pair-only preparation. Section V D recognizes this but treats 192 CNOTs as an untested ADAPT prediction.

An independent optimization of the four occupied-to-virtual pair-excitation factors, in the manuscript's order, gives:

- E=0.6369869810446;
- 1-F=1.7764427147e-4;
- four excitation factors, or 4 x 48=192 CNOTs in the manuscript's ladder compilation.

These reproduce the reported full-template UCCSD energy and fidelity to the displayed precision. Thus, at this operating point, a 192-CNOT pair-only comparator can be directly demonstrated. This is an explicit restricted variational ansatz calculation, not an ADAPT run.

Distinguish this numerical result from the general ADAPT assertion. A zero gradient for a seniority-changing operator at a seniority-zero state follows from the symmetry selection rule. It does not prove that ADAPT selects exactly four factors, selects each once, or achieves a specified error. ADAPT can reuse operators, and a generalized pool has more operators than the occupied-to-virtual subset. Nor does a zero gradient on the symmetry-preserving parameter manifold prove that every unconstrained variational optimum has zero symmetry-changing amplitudes.

### 4.3 Resources at matched fidelity, including zero filter cycles

A fixed-cycle comparison is informative, but a cost ratio at equal cycle number is not a ratio at matched fidelity. The variational input already has infidelity below 2e-4 before filtering. Under that target, it costs 1,312 CNOTs as reported, or 192 after the demonstrated pruning. Adding two expensive rodeo cycles is unnecessary.

Some reproduced operating points, excluding training, energy search, routing, and hardware noise:

| Preparation | Rodeo cycles | Final infidelity | CNOTs per accepted preparation |
|---|---:|---:|---:|
| Full-template UCCSD | 0 | 1.78e-4 | 1,312 |
| Four-pair restricted ansatz | 0 | 1.78e-4 | 192 |
| Full-template UCCSD plus filter | 2 | 1.30e-4 | about 23,000 |
| Refined T=16 plus filter | 2 | 1.68e-3 | about 48,700 |
| Unrefined determinant plus filter | 6 | 1.78e-3 | about 69,700 |
| Refined T=16 plus filter | 6 | 1.74e-4 | about 94,000 |
| Refined T=30 plus filter | 6 | 1.48e-4 | about 117,000 |

The four-pair result uses an independent optimization. The other rows reproduce manuscript values, with slightly more digits. These are operating points, not proven global optima. The paper's factor of roughly three compares states at different achieved errors; it should not be presented as a universal matched-accuracy ratio.

Plot the minimum cost that reaches each requested infidelity, minimizing over cycle count, refinement time, shift, and product-formula step size within a clearly specified search. Keep the error targets above the demonstrated numerical uncertainty. At least include targets around 2e-3, 5e-4, and 2e-4. If a target is below a chosen simulation-error floor, reduce the step size and count the additional gates.

### 4.4 The accuracy floor and the approximate digital filter

Section V C and the Fig. 13 caption say that all inputs reach the approximately 1.1-1.5e-4 floor after six cycles. Table III instead gives 1.8e-3 for the determinant and 1.7e-4 for T=16. The determinant still carries substantial unfiltered excited-state weight at M=6.

Use wording such as: “The high-overlap inputs approach a common product-formula-limited accuracy; the determinant requires additional cycles to approach that regime.” Avoid treating six-cycle points as measured asymptotes. A longer independent calculation retaining the paper's first six random times and appending additional times brings the determinant near 1.8e-4 at ten cycles. That extension is illustrative and is not a statistically converged new benchmark.

With exact evolution on resonance, every K_t leaves the ground-state amplitude unchanged and all filters share H's eigenvectors. Thus P_M tends to p and ideal infidelity can vanish.

For the product-formula evolution U_tilde(t), define

\[
\widetilde K_t=\frac{I+e^{iEt}\widetilde U(t)}2,
\qquad
\mathcal K(\rho)=\mathbb E_t[\widetilde K_t\rho\widetilde K_t^\dagger].
\]

For independent schedules per attempt, the accepted ensemble is

\[
\rho_M=\frac{\mathcal K^M(\rho_{\rm in})}
{\operatorname{Tr}\mathcal K^M(\rho_{\rm in})}.
\]

The true ground state is generally no longer preserved by every approximate filter. A dominant eigenoperator of this trace-decreasing map can yield a biased long-run state while acceptance continues to decay. This gives a precise interpretation of the Section V accuracy floor; it is not generally projection onto one effective Hamiltonian or a literal mixture of the separate filters' ground states.

As an additional independent analysis, deterministic quadrature of this map in the six-dimensional pair sector, using Section V's positive half-Gaussian times, sigma=4, and step ceiling 0.25, gives asymptotic infidelity about 1.1399e-4 and leading map eigenvalue 0.9978208. Doubling quadrature order leaves these figures unchanged to the precision needed here. This supports the approximate floor claimed in V, but also implies that long-run acceptance decays rather than approaching the ideal p. At step ceiling 0.5 the analogous floor is about 1.9345e-3, also consistent with the reported behavior.

These are additional calculations, not manuscript data or hardware results. They show how to reconcile the ideal theorem with the digital results. The description “the fidelity never converges” at step ceiling 1 should likewise mean that it fails to approach the desired accuracy in the displayed range; convergence to a biased state is not excluded.

### 4.5 Coherent depth, accepted-state cost, and setup costs

The expected early-termination cost is

C_M = [C_prep + sum_m E(P_{m-1} c_m)] / E(P_M).

This is the right repeated-preparation metric when the time schedule is randomized per attempt. The implementation uses the product inside the expectation. The manuscript's factorization is valid in expectation when the current cycle time is independent of earlier times; writing the expression above avoids ambiguity.

Add separate resource columns for two-qubit count per attempt, two-qubit depth of a surviving attempt, expected total count per accepted preparation, ancillas, and measurement/reset rounds. Rejected attempts increase total work but do not concatenate into the coherent depth of the successful state. VQE training repetitions likewise do not increase the depth of its final preparation circuit.

The natural amortized comparison for Q requested accepted states is C_total(Q)=C_setup+Q*C_accepted. Estimate training measurements, energy-search measurements, and calibration separately. The eight analytic-gradient optimizer iterations in a noiseless simulation do not determine a hardware shot budget without an estimator, grouping scheme, shot precision, gradient implementation, and stopping criterion.

The energy scan is a common setup cost only if the same scan input and search protocol are explicitly shared. Running the scan with each candidate input charges different preparation costs and gives different peak statistics.

### 4.6 Compilation fairness, alternative encodings, and missing baselines

The chosen ladder compilation makes the comparison transparent. It does not give minimal implementation costs. Merge adjacent exponentials of the same term in symmetric product formulas, remove identity/zero-angle operations, and apply the same compilation effort to every route. Optimized fermionic excitation circuits can improve on independent Pauli ladders; pair hopping in the Hamiltonian can benefit as well, so the net ranking needs measurement rather than assumption.

The seniority-zero encoding uses one qubit per level, with fixed pair number. Retaining the 2k spin-orbital register is reasonable to demonstrate a method intended for seniority-breaking extensions. It does not establish the lowest resource cost for the pure pairing Hamiltonian. A pair-encoded variational and filtering baseline would make that distinction explicit. Symmetry-projected BCS-like inputs are another relevant pairing-specific comparison.

Add the prolonged k=3 coarse state followed directly by the same filter. The paper has already prepared this state, and its p=0.9637 overlap tests whether refinement is needed after a useful coarse solution. Also add an ordinary adiabatic ramp from the diagonal/HF Hamiltonian, optionally followed by rodeo. At the principal k=2 start the coarse state is a determinant, so this comparator is necessary to isolate the benefit of the specific refinement path from generic adiabatic preconditioning.

It is unnecessary to implement every QPE/QSP/QETU method for a carefully scoped benchmark. However, “optimizer-free geometric convergence” is not unique to rodeo, and an unamplified post-selection cost proportional to 1/p is not a universal lower bound for all eigenstate-preparation algorithms. Amplitude amplification changes the overlap dependence under additional coherent-access assumptions. Restrict general claims accordingly.

Hardware routing needs separate treatment. The two qubits within each time-reversed pair are adjacent in the chosen ordering; the two different levels in a general pair transfer need not be adjacent on a line, square lattice, or heavy-hex device. Cancellation of Jordan-Wigner strings does not remove connectivity costs. A logical workload also needs a rotation-synthesis tolerance and fault-tolerant resource model before being a hardware specification.

## 5. Coherence of Sections III, IV, and V

### 5.1 Assumptions and operating points

| Claim or result | Relevant assumptions and operating point | Supported interpretation |
|---|---|---|
| III coarse-state preparation | Factorized unitary ansatz, optimized numerically | Valid circuit family; coarse fidelity must be supplied or bounded |
| III prolongation | Nested bare Hamiltonians, empty appended modes | Gate-free embedding and the stated expectation-value identity |
| III adiabatic guarantee | Exact coarse eigenstate, isolated branch, continuous ideal evolution | Conditional convergence bound; not the complete digital error budget |
| IV B refinement | T=25 or 30, exact midpoint-slice exponentials, typically 160 slices | Well-resolved ideal-refinement demonstration |
| IV F filtering | T=25 refined input, sigma=3, exact controlled evolution, 40 schedules | Illustration of geometric ideal filtering, with finite-sampling uncertainty |
| V B-C circuit comparison | Product-formula refinement at T=16 or 30; sigma=4; filter step ceiling 0.25 | Cost and error of a specific digital implementation |

The different fidelities are not inherently contradictory because these are different operating points. The paper should not use the ideal T=25, sigma=3 residual errors as if they were achieved by the counted T=30 or T=16 digital circuits. An explicit table like this would make the logical progression clear.

### 5.2 Why successful preconditioning need not save gates

The ideal large-gap expression is

\[
1-F_M\simeq\frac{1-p}{p}\,2^{-M}.
\]

Therefore increasing the input overlap reduces the number of required cycles logarithmically in its odds ratio. Moving from the determinant p=0.86764 to the ideal T=25 refined input p=0.99757 saves roughly six ideal cycles at a sufficiently demanding common target. The failure/retry advantage is modest because the determinant already has high overlap: asymptotic success is about 0.868 instead of approximately 0.998. Whether those saved cycles repay the many thousands of refinement gates is a separate resource question.

This is exactly why Section III's useful preconditioning mechanism and Section V's lack of a demonstrated advantage can both be true. There is no physical contradiction. The pipeline's value depends on the cost of obtaining the improved p relative to simply filtering or cheaply improving the original input. The already highly accurate variational input makes that trade unfavorable here.

## 6. Consolidated revision priorities

1. **Complete the pipeline assumptions and error budget.** Include coarse-input error, the nested-Hamiltonian condition, the isolated-branch condition, and the multigroup commutator sum. Distinguish an optimizer-free circuit from the calibration used to choose its controls.
2. **Revise the scaling interpretation.** Retain the intercept in the reported fit, remove the claim that inverse-square behavior is excluded, and treat the coupling and shift scans as sensitivity studies unless the additional response assumptions are controlled.
3. **Correct statistical reporting and reproducibility.** Compare infidelities with the exact ensemble prediction, add uncertainty intervals, and standardize the random schedules across the text, figures, and notebook.
4. **Use the demonstrated pair-only comparator.** Include the 192-CNOT four-pair preparation and distinguish it from an unperformed ADAPT calculation. Prune the coarse preparation as well when comparing optimized resource costs.
5. **Compare minimum cost at specified accuracy.** Allow zero rodeo cycles and optimize the preparation/filter choices within a stated search. Avoid fixed-cycle ratios presented as matched-accuracy results.
6. **Separate resource currencies.** Report coherent depth, gates per attempt, gates per accepted state, setup/training measurements, energy search, and calibration separately; include amortization when many states are requested.
7. **Make the ideal/digital distinction explicit.** Include an operating-point table like the one in Section 5.1 of this assessment, correct the six-cycle common-floor claim, and explain the digital filter's loss of exact target preservation.
8. **Validate the nontrivial coarse evolution at the counted circuit level.** Split its noncommuting coarse terms rather than treating the complete embedded coarse exponential as exact at the cost of one product-formula step.
9. **Add the most informative missing baselines.** Test the prolonged k=3 state plus the same filter and an ordinary adiabatic ramp plus filter; otherwise limit the efficiency claims to the implementations actually compared.
10. **Specify compilation and hardware assumptions.** Apply comparable circuit simplification to every route, distinguish a general spin-orbital implementation from pair-encoded costs, and either schedule a hardware circuit or retitle Section V as a gate-cost comparison.
11. **Qualify the energy-search demonstration.** Select the fine interval using only the coarse scan if claiming operational discovery, and report the uncertainty in the estimated peak.

The retained conclusion is that a valid conditional pipeline is demonstrated accurately in ideal simulation, its digitization has a measurable cost and accuracy floor, and this pairing example offers no demonstrated preparation-cost advantage over inexpensive symmetry-adapted inputs. The corrections concern assumptions, statistical reporting, and efficiency/scaling interpretations rather than a failure of the underlying pairing physics or the exact rodeo filter.

## Primary literature checked

- Choi et al., [Rodeo Algorithm for Quantum Computing](https://arxiv.org/abs/2009.04092), PRL 127, 040505 (2021).
- Bogner et al., [Quantum State Preparation with Resolution Refinement](https://arxiv.org/html/2511.14732v2), including the assumptions in its supplemental scaling derivation.
- Jansen, Ruskai, and Seiler, [Bounds for the adiabatic approximation with applications to quantum computation](https://arxiv.org/pdf/quant-ph/0603175), especially Theorem 3.
- Sadia et al., [Digital quantum state preparation by nucleation](https://arxiv.org/html/2609.37527v1), especially its continuous adiabatic and two-group digitization bounds.
- Grimsley et al., [An adaptive variational algorithm for exact molecular simulations on a quantum computer](https://www.nature.com/articles/s41467-019-10988-2), including operator reuse in ADAPT.
- Yordanov et al., [Efficient quantum circuits for quantum computational chemistry](https://doi.org/10.1103/PhysRevA.102.062612).
- Lacroix, [Symmetry-Assisted Preparation of Entangled Many-Body States on a Quantum Computer](https://doi.org/10.1103/PhysRevLett.125.230502).
- Lin and Tong, [Near-optimal ground state preparation](https://quantum-journal.org/papers/q-2020-12-14-372/).
- Dong, Lin, and Tong, [Ground-State Preparation and Energy Estimation on Early Fault-Tolerant Quantum Computers via Quantum Eigenvalue Transformation of Unitary Matrices](https://doi.org/10.1103/PRXQuantum.3.040305).
