# Critical assessment of `eigenq_pairing.tex` (draft of 2026-07-13)

Target journals considered: Physical Review C, Physics Letters B.

## Overall verdict

The manuscript is technically sound and unusually well validated (every number traces
to exact diagonalisation; gate counts are explicit and verified). It is at the level
of a strong, complete *draft*. It is not yet submittable, for three reasons: (i) the
novelty relative to Bogner et al., PLB 875, 140363 (2026) is not sharply enough
delineated; (ii) the demonstrated instance is too easy — the coarse VQE stage is
trivial (k=2 is a single determinant), only one coupling (g=1) and one particle
number are shown, and all evolutions are exact rather than Trotterised; (iii) the
bibliography ignores essentially all prior quantum-computing work on the pairing
model and the post-2021 rodeo literature, which a PRC referee will notice
immediately. All three are fixable with the existing codebase.

## Journal fit

- The documentclass currently says `prresearch` — inconsistent with the stated PRC/PLB target; decide and switch.
- **PRC** is the better fit: methodology-plus-benchmark papers on schematic nuclear models (Lipkin, pairing, seniority) have a long track record there, length is unconstrained, and the referees will be the right community. Requires strengthening the nuclear-structure framing and the prior-art coverage (see references below).
- **PLB** is harder to argue: the companion resolution-refinement paper is *already* in PLB, so a second PLB letter demonstrating the same idea on another model invites an "incremental" verdict. If PLB is preferred, the letter must be cut to ~6 pages and led by the one genuinely new quantitative result: the gate-level depth budget showing the rodeo-terminated pipeline beating full-space UCCSD–VQE by orders of magnitude at higher fidelity.
- A third option worth considering given the current class file: PR Research or PRX Quantum, where the algorithm-centric framing needs no change.

## What is actually new here (state it explicitly, early)

Referees will ask "what is beyond Ref. [Bogner2026] and Refs. [Choi2021, Qian2024]?"
The introduction should answer in one paragraph: (a) first *combination* of
resolution refinement with rodeo filtering, closing the overlap loophole of
projective filters; (b) first application of refinement to a *basis* (single-particle
level) hierarchy in a fermionic pairing problem, with a gate-free prolongation;
(c) complete gate-level compilation and a quantitative depth comparison against the
variational route; (d) reproducible reference implementation. Currently these points
are scattered and implicit.

## Required scientific additions (ordered by referee impact)

1. **Make the coarse stage nontrivial.** At N=4, k=2 the "cheap VQE" prepares a
   single Slater determinant with zero gates — the pipeline's first stage does no
   work in the demonstrated instance, and the prolonged state is just the HF state.
   Add at least one chain where the coarse solver matters, e.g. k=3 → k=4 (coarse
   UCCSD at 272 CNOTs, prolongation, refinement, rodeo) and/or N=6 with k=3 coarse.
   Compare refined overlap and total gate cost of k=2→4 vs k=3→4 — this is also an
   interesting result in itself (where in the hierarchy to start).

2. **Scan the coupling g.** Only g=1 is shown. The method's claims live or die on
   the spectral gap along the refinement path; show overlap/energy vs T for weak
   (g≈0.2), intermediate and strong (g≈2–4) pairing, where the low-lying gap
   shrinks. A plot of required T* vs 1/ΔE(g) would directly substantiate the
   T ~ 1/ΔE claim, which is currently asserted from one data point.

3. **Deliver the promised Trotter-error quantification.** Sec. III.D says the
   Trotter error is "quantified in Sec. VI", but the Discussion only makes
   qualitative remarks — an internal inconsistency a referee will flag. Run the
   refinement and the controlled evolution with an actual first/second-order Trotter
   decomposition at several n_s, n_T and show the induced infidelity and energy
   shift, confirming the "self-correcting" claim numerically. This also turns the
   symbolic entries n_s, n_T of Table III into concrete numbers, without which the
   depth budget of Sec. V.B cannot be checked.

4. **Make the depth comparison fair and explicit.** The current comparison
   (1312 CNOTs × optimisation loop vs one pipeline pass) can be attacked: the VQE
   optimisation is a one-time classical cost, after which *each* preparation is 1312
   CNOTs, while each pipeline preparation costs n_s×296 + M×420·n_T CNOTs — for
   n_s of order tens this is *not* obviously smaller. Report total CNOTs per
   accepted preparation for both routes, at matched final fidelity, including the
   1/P_M repetition factor. The honest statement (pipeline wins on *accuracy per
   coherent depth* and removes the optimiser, not necessarily on single-shot depth)
   is still a strong one.

5. **Shot noise and sampling.** All results are statevector-level. Add one figure or
   paragraph with finite-shot estimates of the acceptance and of the energy peak
   position in the rodeo scan (the acceptance is a Bernoulli estimate — this is
   cheap to simulate) so the hardware claims are grounded.

6. **Justify μ = E_low + 0.6.** The energy shift placing tracked states below the
   null space of P† is ad hoc as written; give the criterion (e.g. μ between E_low
   plus the coarse gap and the bottom of the null-space band) and a one-line
   sensitivity statement.

7. **Exact solvability.** The constant-pairing Hamiltonian is Richardson-integrable;
   this must be acknowledged (Richardson 1963; Dukelsky–Pittel–Sierra RMP 2004),
   both as the reason it is an ideal benchmark and to preempt "why not solve it
   classically" — the answer (transparent testbed for a scalable pipeline) should be
   stated.

8. **Magne — factual correction.** Magne is QuNorth's machine in **Copenhagen,
   Denmark** (Atom Computing hardware, Microsoft software stack; ~50 logical /
   1200+ physical qubits, delivered Jan 2026). The text/abstract implies it in a way
   that a Nordic referee will read as Norwegian; also verify and cite its
   mid-circuit-measurement capability before claiming the pipeline targets it, and
   add a citable reference or URL footnote for the system.

## Presentation and mechanics

9. Abstract is ~350 words and reads as a summary; both PRC and PLB will want it
   tightened (aim ≤ 200 words, drop the gate-count sentence or compress it).
10. Figures are numbered fig3–fig8 in `figs/`; regenerate with publication settings
    (fontsize matching column width, vector PDF rather than PNG if possible —
    `make figures` currently extracts PNGs).
11. Switch from inline `thebibliography` to BibTeX (`refs.bib` already exists in
    `paper/`), add DOIs and arXiv IDs throughout.
12. **Uncited bibitems**: `LeeReview2023`, `LeeUCC2019`, `NielsenChuang` are in the
    bibliography but never cited in the text — cite them (LeeReview2023 belongs in
    the introduction next to the rodeo references; LeeUCC2019 in Sec. III.A) or
    remove.
13. Code availability: "provided as a single Jupyter notebook" is not enough for
    APS/Elsevier in 2026 — deposit the repository (GitHub + Zenodo DOI) and cite it.
14. Add funding/grant numbers to the acknowledgments (QuNorth flagship, EIGEN-Q,
    RCN/FRIB grants as applicable); PRC requires them.
15. Author list: CLAUDE.md convention (alphabetical by surname) is satisfied — keep
    it when adding authors.

## Reference update (concrete additions)

Prior quantum computing on pairing/schematic nuclear models — currently absent, and
several referees will come from this community:

- E. Ovrum and M. Hjorth-Jensen, *Quantum computation algorithm for many-body studies*, arXiv:0705.1928 (2007) — the first pairing-Hamiltonian quantum algorithm, and by one of the present authors.
- Z.-H. Jiang and J. Pei (author list to be checked against the journal page), *Quantum computing of the pairing Hamiltonian at finite temperature*, Phys. Rev. C **107**, 044308 (2023).
- P. Zhang, D. Lacroix et al., *Neutron-proton pairing correlations described on quantum computers*, Phys. Rev. C **110**, 064320 (2024).
- D. Lacroix et al., *Symmetry breaking and restoration for many-body problems treated on quantum computers*, arXiv:2310.17996 (review, 2023).
- A. Pérez-Obiol et al., *Nuclear shell-model simulation in digital quantum computers*, Sci. Rep. **13** (2023), arXiv:2302.03641.
- M. Q. Hlatshwayo et al., *Simulating excited states of the Lipkin model on a quantum computer*, Phys. Rev. C **106**, 024319 (2022), arXiv:2203.01478.

Rodeo-algorithm literature beyond 2021/2024 — needed because the paper uses Gaussian
time sampling, which is no longer state of the art:

- T. D. Cohen and H. Oh, *Optimizing the rodeo projection algorithm*, Phys. Rev. A **108**, 032422 (2023) — deterministic time choices beat Gaussian sampling.
- M. Patkowski et al., *Improved rodeo algorithm performance for spectral functions and state preparation*, arXiv:2602.05978 (2026) — geometric time series; note Patkowski is a co-author of the resolution-refinement paper, so omitting this would look odd.
- (verify) M. Bee-Lindgren et al., *Rodeo algorithm with controlled reversal gates*, arXiv:2208.13557 — hardware-oriented variant; check citation data before adding.
- (optional) Qudit rodeo formulation, arXiv:2603.16049 (2026).
- (verify) I. Stetcu, A. Baroni, J. Carlson, projection/state-preparation algorithm on quantum computers, Phys. Rev. C (2022) — closely related projective filtering from the nuclear community; confirm exact reference.

Competing state-preparation filters — the introduction claims QPE/rodeo are the
projective options; referees may ask about polynomial filters:

- L. Lin and Y. Tong, *Near-optimal ground state preparation*, Quantum **4**, 372 (2020).
- Y. Dong, L. Lin, Y. Tong, *Ground-state preparation ... via quantum eigenvalue transformation of unitaries (QETU)*, PRX Quantum **3**, 040305 (2022).
  One sentence positioning the rodeo algorithm against QSP/QETU filters (rodeo: no block-encoding, mid-circuit-measurement friendly) would strengthen Sec. III.D.

Context/reviews:

- M. Larocca et al., *Barren plateaus in variational quantum computing*, Nat. Rev. Phys. **7**, 174 (2025) — cite alongside McClean2018; also soften the barren-plateau motivation slightly (8–16 qubit UCCSD is not in the barren-plateau regime; the honest argument is optimiser cost + ansatz error floor, which the paper itself demonstrates).
- C. W. Bauer, Z. Davoudi, N. Klco, M. J. Savage, *Quantum simulation for high-energy physics* / Nat. Rev. Phys. **5**, 420 (2023), and D. Beck et al., *Quantum information science and technology for nuclear physics*, arXiv:2303.00113 — for the PRC framing.

Exact solvability:

- R. W. Richardson, Phys. Lett. **3**, 277 (1963); R. W. Richardson and N. Sherman, Nucl. Phys. **52**, 221 (1964).
- J. Dukelsky, S. Pittel, G. Sierra, Rev. Mod. Phys. **76**, 643 (2004).

Existing references verified: Bogner2026 = PLB **875**, 140363 (2026) / arXiv:2511.14732 ✓; Choi2021, Qian2024, Larocca-adjacent entries all check out. All "(verify)" items above should be checked against INSPIRE/ADS before submission — do not trust reference details from memory or secondary sources.

## Suggested revision order

1. Items 3 + 4 (Trotterised pipeline, honest depth budget) — new notebook cells, feeds Table III and one new figure.
2. Items 1 + 2 (nontrivial coarse stage, g-scan) — two new figures.
3. Reference overhaul + Richardson framing + Magne correction (items 7, 8, 12, and the list above).
4. Compression pass (abstract, figures, BibTeX) and journal decision.

After each change: `make test && make check`; after paper edits: `make paper`.
