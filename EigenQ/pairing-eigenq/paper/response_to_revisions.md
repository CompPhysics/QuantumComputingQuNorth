# Response to `revisions.txt` — what was checked, what was changed, what is still open

Date: 2026-10-02.  Manuscript: `paper/eigenq_pairing.tex` (previous version kept as
`paper/eigenq_pairing_v1_backup.tex`).  All numbers below were recomputed with
`pairinglib` (new tagged notebook section `EXT:BASELINES`, cells 83–86, and the
standalone scripts `scripts/baselines_rodeo.py`, `scripts/sizes_observables.py`).

## Verdict on the review

Every quantitative claim in the review was reproduced; none was found to be wrong.
The three substantive points (energy error vs. infidelity, missing baselines,
uncounted coarse term) are correct and change the conclusions of the paper.  The
manuscript has been rewritten around a matched-observable comparison, and the
central claim is now stated as a negative-plus-qualification result: at the sizes
we can validate exactly, refinement followed by the rodeo filter does *not* beat
filtering a full-space UCCSD input, and is roughly tied with filtering the
unrefined prolonged state; its merit is that it is optimizer- and measurement-free.

## Point by point

### 1. Energy error vs. infidelity — **confirmed, fixed**
Reproduced with the gate-level circuit (`gate_uccsd_state`), k=4, N=4, g=1:

| f | ΔE (UCCSD) | actual 1−F (UCCSD) | reviewer |
|---|---|---|---|
| 0 | 1.4385e-3 | 1.7764e-4 | 1.7764e-4 ✓ |
| 0.05 | 2.441e-3 | 3.03e-4 | — |
| 0.2 | 8.181e-3 | 1.0144e-3 | 1.0144e-3 ✓ |
| 0.5 | 3.2143e-2 | 3.7202e-3 | 3.7202e-3 ✓ |

Changes: abstract rewritten; Table I discussion adds the spectral inequality
⟨H⟩−E₀ ≥ Δ(1−F) and the actual infidelity; Table III (particle-hole) gains a
1−F_UCCSD column; Table II and Table V now list infidelities only; the
"one-to-three orders of magnitude", "eightfold/sixfold" statements are removed;
Sec. VI.A states that the UCCSD errors are optimizer results, not an expressivity
proof.

### 2. Baselines — **confirmed, new figure and table added (Fig. 17, Table II)**
Same Trotterized filter (order 2, δt=0.25, σ=4, 40 common time sets, seed 3),
acceptance-weighted ensemble infidelity, cost per accepted preparation with early
termination (Eq. 10 of the paper), coarse term included in the refinement cost:

| input | p | C_prep | M=2: 1−F / P / CNOT | M=6: 1−F / P / CNOT |
|---|---|---|---|---|
| prolonged determinant | 0.868 | 0 | 2.9e-2 / 0.887 / 2.3e4 | 1.8e-3 / 0.854 / 7.0e4 |
| UCCSD full space | 0.99982 | 1312 | 1.3e-4 / 0.996 / 2.3e4 | 1.1e-4 / 0.986 / 6.8e4 |
| refined T*=16 (n_s=43) | 0.9922 | 26 660 | 1.7e-3 / 0.989 / 4.9e4 | 1.7e-4 / 0.979 / 9.4e4 |
| refined T=30 (n_s=80) | 0.9977 | 49 600 | 5.3e-4 / 0.994 / 7.2e4 | 1.5e-4 / 0.985 / 1.2e5 |

The reviewer's diagnostic (UCCSD+rodeo 1.35e-4 / ~22 500 CNOTs; refined 6.4e-4 /
~68 800) is reproduced within sampling noise.  All inputs converge to the same
~1.1–1.5e-4 floor set by the product formula at δt=0.25.  Also computed and
reported: UCCSD+rodeo with the particle-hole term (Table V, last column), a
size/particle-number table for (k,N)=(4,4),(6,4),(4,6),(6,6) (Table III, Sec. V.E)
— same ranking at every size — and level occupation numbers + correlation energy
(Table IV).  The T*=16 refinement is cheaper than T=30 at equal fidelity beyond
M=4, which is stated (the old operating point was not economical).
Not done: ADAPT-VQE / pair-adapted comparator (stated as the natural next step);
a joint optimization of T and M beyond the two T values shown; uncertainty bands
(40 time sets; scatter quoted as ±30 % at the smallest infidelities).

### 3. The coarse term A = P(H_low−μ)P† — **confirmed, circuit given and counted**
`refine_state_trotter` evolves A as an exact matrix factor; `cnot_counts` did not
count it.  Verified numerically (`A` has rank 1 for 2→4, rank 15 for 3→4, supported
only on states with the appended orbitals empty; exp(−icAδt) equals the
conditional phase on |HF⟩).  Circuit: for 2→4, X⊗4·C³P(φ)·X⊗4 = 15-term Z phase
polynomial → 14 CNOTs + 15 R_z along a Gray code (34 CNOTs with naive ladders);
for 3→4, every one of the 33 rotations of a k=3 step gets a double control (6 CNOTs
each) → ≈350 CNOTs per application.  Consequences: Table II refinement 47 360 →
49 600 (+5 %); Table V refinement 436 480 → 438 720; **Sec. IV.C reversed**: the
3→4 chain's shorter T* (9 vs 16) is offset by the conditioned coarse term, 2.7e4
vs 3.1e4 CNOTs — no gate saving (the old "∼2×10⁴ CNOTs saved" claim is removed).
Sec. V now states explicitly that product-formula evolutions are simulated block by
block (exact block exponentials in the N sector), that this is not a gate-by-gate
execution of every rotation/CNOT, and that CNOT counts are quoted, not scheduled
depths.

### 4. Convergence headline, bound, ensemble, detuning — **confirmed, fixed**
Fig. 8 reproduced (T=25 input, p=0.9976, σ=3, 40 sets): 1−F = 1.4e-3, 6.5e-4,
2.1e-5, 2.0e-6 at M = 1, 2, 6, 10.  Abstract/Sec. IV.E now say two cycles →
5–6e-4, ten cycles → 2e-6.  Eq. (7) replaced by the exact ensemble expression and
the q_max bound; 2⁻ᴹ identified as the Δσ≫1 limit (here Δ_min σ ≥ 5, so the two
coincide numerically).  The dotted lines are now called "analytic ensemble
expectation", and the ±30 % scatter of the 40-set average is stated.  Ensemble
specified: the notebook averages normalized fidelities uniformly over time sets;
the acceptance-weighted average differs by <1 % (refined) and a few % (prolonged)
— both are now quoted.  Scan cost stated (≈2×10³ four-cycle runs for the fine
scan).  Detuning study added: using the scan estimate E=0.6352 instead of E₀
changes nothing; for E−E₀ = 0.05/0.1/0.2 the six-cycle infidelity moves by
10 %/35 %/×2.7 while acceptance falls to 0.971/0.895/0.656 — acceptance no longer
tends to p, as the reviewer notes.

### 5. Inverse-gap and scalability claims — **confirmed, qualified**
Sec. III.C now states the assumption of Ref. [Bogner2026] (bounded first-order
correction / transition matrix elements).  Sec. IV.D: the μ-scan changes the
operator as well as the gap → "empirical trend for this path family", not an
isolation of gap dependence; shift-controlled initial gap vs. physical gap
distinguished.  Shift-buffer "ceiling" statement corrected (initial gap =
min(E₁ˡᵒʷ−E₀ˡᵒʷ, μ−E₀ˡᵒʷ), saturating at the coarse excitation gap).  Fig. 2/3
described as a UCCSD-only basis-size study (sector ∝k⁴, seniority-zero ∝k²), not
many-particle scaling.  Universal "overlap vanishes at scale" softened to "for a
fixed unoptimized reference".  The pipeline was extended to (6,4), (4,6), (6,6)
with exact slices — presented as a controlled small-system demonstration.

### 6. Positioning — **confirmed, rewritten**
Intro now cites: Choi et al. 2021 (adiabatic preconditioning in the original rodeo
paper), Patkowski et al., PRA 113, 052442 (2026) / arXiv:2510.19039 (fusion
method), and the explicit "refinement then rodeo" recommendation of
Bogner et al. 2026.  "First combination … closing the overlap loophole" removed;
contribution list corrected to five items and reframed (implementation for a
fermionic basis hierarchy, quantitative benchmark, seniority-breaking test,
documented resource analysis, reproducible notebook).  Title changed to
"…: implementation and assessment on the pairing Hamiltonian".

### Technical and literature corrections — all applied
- Brillouin: reference is not the stationary HF determinant; theorem does not
  fail; reference relaxation vs. correlation distinguished (optimized-reference
  comparison not done — flagged).
- RPA: statement removed; Richardson-model pp-RPA literature cited
  (Hirsch et al., Ann. Phys. 296, 187 (2002); N. Dinh Dang, PRC 74, 034326 (2006)).
- 2k-qubit encoding: "only option" → convenient general choice; compression via
  N and S_z noted.
- ADAPT wrong-symmetry inputs: energy filter cannot fix zero overlap, wrong N
  sector, or exact degeneracies — qualified.
- Ref. Aychet-Claisse et al.: attribution corrected (pure pairing Hamiltonian, g
  drives the transition; V_ph not part of that model).
- Ref. [42] → PhD thesis of E. A. Ruiz Guzmán, Univ. Paris-Saclay (2023), arXiv:2310.17996.
- Ref. [44] → Eur. Phys. J. A 61, 263 (2025).
- Trotter step convention h = t/⌈|t|/δt_max⌉ stated; "one common effective
  Hamiltonian" qualified; self-correction not claimed as a guarantee.
- Seven qubits → 2→3 chain (trivial coarse state) on 7 qubits or 3→4 chain on 9.
- "four contributions … (v)" fixed.
- Bogner2026 bibitem: full author list added.

### Presentation / hardware
Vendor-by-vendor paragraph replaced by a requirements paragraph (qubits, ancilla
reuse, number of arbitrary-angle rotations, gate budgets 2e4–1e6, fault counts at
demonstrated 2Q error rates).  Kim 2023 now cited only for the Eagle-class error
rate; Helios cited from arXiv:2511.05465; Rigetti press release and the Google
QEC paper dropped (not load-bearing); neutral-atom demonstrations described as
error detection/partial correction with mid-circuit readout.

## Still open (needs author decisions or new work)
1. ADAPT-VQE / pair-adapted comparator (the cheapest route by our count: 192 CNOTs + 2 cycles).
2. Joint optimization of (T, M) and statistical uncertainty bands (more time sets).
3. Optimized-reference (true HF) comparison for the particle-hole model.
4. Whether to keep the particle-hole section at its present length given the
   corrected conclusion (it now confirms the negative result rather than a win).
5. Code-availability GitHub link; `Chapter 4 of Ref. [MHJbook]` is still a
   lecture-notes citation.
6. A regime where refinement wins (poor variational input) is asserted as the
   open question, not demonstrated — a reviewer may still ask for one.

## Repository changes
- `paper/eigenq_pairing.tex` (revised), `paper/eigenq_pairing.pdf` (rebuilt, 23 pp.),
  `paper/eigenq_pairing_v1_backup.tex` (previous version), `paper/figs/fig17.png` (new).
- `notebooks/ResolutionRefPairing.ipynb`: 4 new cells tagged `<!-- EXT:BASELINES -->`
  (executed; `check_notebook.py … 17` passes).
- `scripts/baselines_rodeo.py`, `scripts/sizes_observables.py`,
  `scripts/extend_notebook_baselines.py` (new).
- `Makefile`: `check` target now expects 17 figures.
- `CLAUDE.md`: new anchors and the "compare infidelities, not energy errors" rule.
