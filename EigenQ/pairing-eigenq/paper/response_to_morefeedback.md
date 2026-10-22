# Response to the final pre-submission assessment (21 October 2026)

Assessment: `paper/morefeedback.md`.  Every numerical claim of the assessment was recomputed
with the repository code before any change was made.  The recomputation and the new search
live in notebook section 23 (`scripts/frontier_refine.py`, figure 25) and in the library
(`AveragedFilterMap.superoperator`, a faster `trotter_path`), with a new test in
`tests/test_rodeo.py`; the previous main text and Supplemental Material are kept as
`paper/eigenq_pairing_v9_backup.tex` and `paper/eigenq_pairing_supp_v9_backup.tex`.

## Verdict per point

| § | Point | Verdict | What changed |
|---|-------|---------|--------------|
| 1 | Resolved since 14 October | Retained as is. | — |
| 2 | The cost frontier is not converged: the coarse control grid cannot support percent-level rankings, the "exact 19 %" coarse-term explanation, or the filter-cycle claim at g=2 | **Agree, confirmed.** All six counterexample rows reproduce to every quoted digit (k=4 ramp: 1.63292e-3 at (4.5, 11), 4.34611e-4 at (4.75, 16), 1.57748e-4 at (4.75, 22); k=8 ramp 2.23472e-4 at (5, 40), 6.15426e-5 at (5.2, 60); k=8 g=2 refinement 8.94159e-5 at (16, 10, 160)). | New search on fine grids of T (steps 0.25/0.5) and μ_buf (nine/ten values), integer n_s by bracketing + bisection at every target and every M = 0..6, exact channel; three search levels L0 (old grid) / L1 (fine grid) / L2 (integer refinement) and the number of controls within 10 % of the minimum are reported (SM Table S6, Fig. 25).  Table V rebuilt from L2 and relabelled "best cost found under the stated search".  Results: every adiabatic minimum of the coarse grid falls (k=4: −31…−51 %; k=8, g=1: −25…−40 %; g=2: −8…−21 %), the integer refinement at the fine grid moves the minima by ≤ 12 %, a 0.05-step T grid around the g=1 ramp minima returns the same values, and the route ratio refinement/ramp becomes 0.85–0.95 at k=4, 1.22–1.42 at k=8 g=1 and 0.86–1.01 at g=2 — i.e. the "19 % coarse-term" separation is not the whole story at g=1 (the shifted path also needs a somewhat longer time at the tighter targets) and at g=2 the refinement reaches every target unfiltered while the ramp needs one and three cycles, as the assessment anticipated.  Abstract, introduction, Secs. V C / V E and the discussion rewritten accordingly; "exactly the conditioned phase of the coarse term", "travel at the same speed" and the g=2 filter-cycle statement withdrawn or made conditional. |
| 3 | Eight-qubit ADAPT row of Table V: the four-operator prefix (192 CNOTs) meets every listed target | **Agree.** 1−F = 1.7764427e-4 at 4 operators, reproduced. | Table V row corrected to 192 (4) at the three k=4 targets; `frontier_checks.py` and `frontier_refine.py` admit every ADAPT prefix at both sizes (the k=4 search had admitted only the exact five-operator state); the text notes that the four-pair preparation and the four-operator ADAPT prefix are the same cost. |
| 4 | Random-time protocol: the implemented digital filter uses t = \|z\| (half-normal); signed and positive-time channels differ with a product formula | **Agree.** Reproduced: positive-time floors 1.139876e-4 / 0.1539359 versus signed 1.138456e-4 / 0.1291596 at δt = 0.25 / 1 (k=4, σ=4). | Sec. III D now states the convention: signed Gaussian draws for the ideal filter (where only even functions of t enter), positive times t = \|z\| for every digital result, sampled and averaged alike; Fig. 2 caption, Secs. V B and V C updated; the signed-time channel and the size of the difference are quoted in Sec. V B.  No rerun was needed. |
| 5 | Results independently reproduced | Preserved. | The supplement's "magnitude alone does not isolate the mechanism" wording is now also used in the discussion (bounds vs practice). |
| 6.1 | Abstract ADAPT cost statement inaccurate | **Agree.** | Replaced by the actual costs with target and size: 192 CNOTs at eight qubits (the same as the four-pair preparation) and 1 100–1 500 at sixteen for 10⁻⁴, with the scoped ratio to the adiabatic routes. |
| 6.2 | Averaging-bound arithmetic on p. 6 | **Agree.** 1/p − 1 = 15 % at p = 0.87, ≈ 1 % at p = 0.99. | Corrected; the observed differences (notebook Sec. 17: 2.935e-2 vs 2.888e-2 at M=2 for the determinant, i.e. 1.6 %) are stated as staying well below the bound. |
| 6.3 | Keep Eq. (13) an initial-response identity, not a runtime theorem | **Agree.** | Abstract now says the identity "gives a tested rule for how the shift must scale"; the discussion already had the qualification and keeps it. |
| 6.4 | "Only locally" nucleation analogy | **Agree.** | "Incrementally"; a sentence states that the increments are not spatially local (pairing force connects all levels; the coarse term is a many-qubit projector).  Also in the introduction. |
| 6.5 | Universal 1/p cost statement | **Agree.** | Scoped to repetition-and-post-selection filters; the 1/√p of coherent amplitude amplification (Lin and Tong 2020) noted under its access-model assumptions. |
| 6.6 | Compilation-dependent margins | **Agree.** | Every percentage between the two adiabatic routes is now tied to "the reference compilation" and to the stated search; the 19 % figure is quoted as a property of the chosen operating points, not as a general separation. |
| 6.7 | Hardware-specification language | **Agree.** | "Logical gate counts, not scheduled depths"; the ingredients a machine specification would add are listed; error rates dated to the cited experiments (2023, 2023/2025); "offered as its specification" → "give its logical gate counts". |
| 7 | Reproducibility and presentation | **Agree.** | "Data and code availability" section (notebook with executed outputs = the data set; package, notebook, tests and LaTeX sources in the repository and the archived release corresponding to the submitted version).  The URL, the release tag and the DOI remain placeholders **to be filled by the authors**: we did not invent identifiers.  SM build: the duplicate-label/citation warnings came from `xr` importing the main `.aux` (FirstPage/LastPage/LastBibItem labels and \bibcite entries that both documents define); `make paper` now builds the SM against a filtered copy (`eigenq_pairing_xr.aux`) and the SM bibliography is self-contained (Kandala 2017 added) — zero multiply-defined warnings.  Discussion shortened (the concluding paragraph no longer restates the ranking); one occupation/correlation-energy result of Table S4 drawn into the discussion. |
| 8 | Route to submission | Items 1–5 done, within the limits stated for item 7. | — |

## Points where we did not follow the assessment literally

- The search is still a finite one and we say so: "best cost found under the stated search".
  The convergence evidence is the L0 → L1 → L2 sequence and the flatness count, not a
  certificate.  The filter width σ and step δt remain fixed by the spectral width and by the
  floor, as before.
- We kept the sampled fixed-input Tables II and V (now VI) as the record of the protocol a
  device would run; Table V (tab:frontier) is the only table that carries minimum costs.
- The signed-time channel is quoted for comparison only; the paper's protocol is the
  positive-time one that was implemented throughout, so no digital number changed.

## What the refined search changed in the conclusions

All adiabatic minima of Table V are lower than before, by 8–51 %, and every quoted ratio was recomputed from the refined search (SM Table S6, notebook Section 23, figure 25 = the cost frontier of the main text; the old grid figure moved to the SM as Fig. S9).  The statements that changed:

- k=4 (targets 2e-3/5e-4/2e-4): ramp 6.5e3/9.5e3/1.2e4 (T=4.25–4.75, M=0), refinement 6.2e3/8.1e3/1.1e4 (μ_buf=2–2.5, T=4.5–5, M=0): a tie within 15 %, now slightly in the refinement's favor (was: refinement 5 % dearer at 2e-3, 20–26 % cheaper below — "a tie in which the shift-set gap outweighs the coarse-term overhead" was overinterpreting the coarse grid).
- k=8, g=1 (1e-3/3e-4/1e-4): ramp 7.1e4/9.8e4/1.4e5 (T=4.75–5.25, n_s=26/36/50, M=0), refinement 8.7e4/1.2e5/1.9e5 (μ_buf=3–4, T=5–6.5, n_s=27/37/60, M=0): the ramp is cheaper by 22 %, 22 % and 42 %, of which the coarse term accounts for 19 %; "travel at the same speed and differ only by the embedded coarse term" withdrawn.
- k=8, g=2: ramp 2.0e5 (M=0)/3.6e5 (M=1)/5.3e5 (M=3) at T=8; refinement 2.0e5/3.1e5/4.6e5 unfiltered at μ_buf=16–24, T=7–8.5: a tie at 1e-3 and 12–14 % in the refinement's favor below it; the earlier "refinement needs two cycles at 1e-4" is gone (the assessment's T=10 point, 5.16e5, is itself beaten by T=8.5, n_s=141, 4.55e5).
- The only adiabatic operating points that still use filter cycles are the g=2 ramp ones; "the cheapest adiabatic operating points use no filter cycle" is now stated with that exception.
- ADAPT: 192 CNOTs (4 operators) at the three k=4 targets, 912/912/1104 and 1200/1392/1488 at k=8: 80–120× (g=1) and 170–310× (g=2) below the cheapest adiabatic route; the abstract quotes these costs instead of "about a thousand at both sizes, two orders of magnitude below every other route".
- Margins between the two adiabatic routes are stated with the resolution of the search (≈10 %) and the compilation; the surviving rankings are the orders of magnitude (ADAPT ≪ adiabatic routes < single-layer ansatz + filter < determinant + filter at sixteen qubits), not the percentages.

Library changes behind the search: `AveragedFilterMap.superoperator()` (cached d²×d² matrix; `track` is 20× faster at d=70, numbers identical to the node loop — tested) and `trotter_path` applying each factor as two matrix-vector products (6× faster, identical numbers).  Build: `make test` (42 tests) and `make check` (25 figures) pass; `make paper` builds the main text and the SM without unresolved references or multiply-defined labels.
