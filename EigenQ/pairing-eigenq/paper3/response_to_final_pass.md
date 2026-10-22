# Paper III: final pre-submission pass (22 October 2026)

The checklist of the final assessment of paper I (`paper/morefeedback.md`, 21 October 2026)
applied to paper III, point by point.  The new computation lives in notebook Section 14
(`scripts/hubbard_refine_checks.py`, inserted by `scripts/extend_notebook_hubbard_refine.py`
and emitted by `scripts/build_notebook_hubbard.py` on a rebuild; figure 15) and in the library
(`hb.exact_digital_ensemble(..., cache=)`, `hb._eigterms` storing complex eigenvectors); the
previous main text and Supplemental Material are kept as `paper3/eigenq_hubbard_v5_backup.tex`
and `paper3/eigenq_hubbard_supp_v5_backup.tex`.

## Verdict per point

| § (paper-I assessment) | Point, transferred to paper III | Verdict | What changed |
|---|---|---|---|
| 2 | The cost frontier is not converged: coarse control grids cannot carry the quoted ratios ("2.4–3.1 times the ramp", "tied at 10⁻³ at thirty-two qubits", "cheaper by 1.9 and 1.8") | **Applies, confirmed.** Fine T grids (steps of 1–2/t at eight sites, 2–4/t at sixteen), nine to eighteen step counts on a geometric list, integer n_s by bracketing + bisection at every target; every M = 0..6 with the exact quantized-time ensemble at eight sites; M = 0 at sixteen sites, which is complete for the ramp and Path A because one average cycle (1.76e5 CNOTs) exceeds their unfiltered minima at every target.  Three search levels L0/L1/L2 and the number of controls within 10 % reported (SM Sec. S9, Table S7, Figs. S5–S6; main-text Fig. 6 = the cost frontiers). | Every adiabatic minimum of Tables 2–3 falls: at eight sites by 28–44 % (ramp), 20–22 % (Path A), 15–23 % (Path B); at sixteen sites by 25–63 % (ramp) and 18–52 % (Path A).  The integer refinement at the fine grid moves the minima by ≤ 11 %, one to six controls lie within 10 % of each unfiltered minimum, and the two routes separate uniformly: Path A / ramp = 2.95–3.53 at eight sites (was "2.4–3.1") and 1.60–1.86 at sixteen (was "tied at 10⁻³, within 1.9 elsewhere" — the tie was a T-spacing artifact: the ramp's minimum moved from (48, 512) to (44, 244) and halved).  Tables 2–3 rebuilt from L2 and relabelled "best cost found under the stated search"; abstract, introduction, Secs. 5.3/6 and the discussion rewritten; the two grid figures moved to the SM (Figs. S5–S6) and replaced by the cost-frontier figure. |
| 3 | ADAPT prefixes admitted | Already the case in Section 13 (every prefix below 2e-2 at eight sites; every tenth prefix, the 1e-2 crossing and the last state at sixteen). | No change. |
| 4 | Random-time protocol specified consistently | **Applies.** The implemented protocol is t = \|z\|, z ~ N(0, σ²), quantized to n = max(1, round(t/δt)) steps. | Sec. 3.4 and Sec. 5.3 now state it; the sign is immaterial for the exact-ensemble formula and the quantized protocol is the one every digital number refers to. |
| 6.1 | Headline ADAPT/ratio statements | **Applies.** | Abstract, Sec. 5.3, Sec. 6 and the discussion rewritten from the L2 minima with the resolution of the search stated; ratios given as ranges. |
| 6.5 | Universal 1/p statement | **Applies.** | Scoped to repetition-and-post-selection filters; Lin–Tong 1/√p noted (new reference). |
| 6.6 | Compilation-dependent margins | **Applies.** | Every ratio between routes is tied to the shared-flag compilation and to the stated search. |
| 6.7 | Hardware language | Paper III already labels its metric "amortized logical CNOTs on an ideal all-to-all device, not a scheduled depth"; the outlook's "within reach of current devices" is now "in qubit number". | Minor wording. |
| 7 | Reproducibility; SM build | **Applies.** The SM did not build in the container at all: `xr` imported the main `.aux`'s `\bibcite` entries and hyperref failed on them (the Mac build had silently produced duplicate-citation warnings).  `make paper-hubbard` now builds the SM against `eigenq_hubbard_xr.aux` (main `.aux` without `\bibcite`); the SM bibliography is self-contained.  "Code and data availability" extended (executed notebook = data set; repository + archived release); the URL/tag/DOI remain placeholders **to be filled by the authors**. | — |

## What the refined search changed in the conclusions

- Eight sites (1e-3/3e-4/1e-4): ramp 1.0e4/1.6e4/2.4e4 (T=13–19/t, unfiltered; was 1.4e4/2.8e4/3.7e4); Path A 3.0e4/5.3e4/8.5e4 (3.0–3.5× the ramp; was 2.4–3.1×); ADAPT 1.7e4/2.3e4/7.9e4 is now 1.5–3.3× the ramp and never cheaper ("0.8 times the ramp at the second target" withdrawn); pruned UCCSD 5–6×, Path B 9–10× (filtered, 9.8e4/1.4e5/2.1e5).
- Sixteen sites (1e-2/3e-3/1e-3/3e-4): ramp 2.9e4/4.4e4/7.5e4/1.2e5 (T=28–60/t, n_s=94–391, unfiltered; was 7.8e4/7.8e4/1.6e5/1.6e5), Path A 4.6e4/7.2e4/1.2e5/2.2e5 (T=8–26/t; was 9.0e4/1.5e5/1.5e5/2.8e5), ratio 1.60/1.62/1.65/1.86; the eight-site coarse state of the hierarchy is now the refined ramp (16 416 / 24 192 CNOTs with the prolongation instead of 27 792 / 37 008).  ADAPT at 1e-2 (4.8e4) is no longer the cheapest input (1.65× the ramp); below it the variational routes cost 8–11× the ramp (was 4–7×); Path B 4–8× (was 1.4–5.9×).
- "An unfiltered ramp is the cheapest route from 3e-3 down" → "at every target from 1e-2 to 1e-4"; "tied with it at 10⁻³ at thirty-two qubits" and "the cheaper step of the ramp cancels the advantage" withdrawn; the surviving statement is that the projector path narrows the gap from ×3–3.5 to ×1.6–1.9 with the size but does not close it, and the margins are stated with the ≈10 % resolution of the search and the shared-flag compilation.
- The M=0 restriction of the sixteen-site refinement is justified in the text and checked in the notebook: the expected cost of any filtered point is ≥ 1.76e5 CNOTs (one average cycle, even for the bare determinant), above every unfiltered minimum of Table 3.

Build: 11 Hubbard tests pass; notebook Section 14 executed in-process (≈ 1.3 h; `make check-hubbard` → 15 figures); `make paper-hubbard` builds main text (18 pp.) and SM (10 pp.) without unresolved references or multiply-defined labels.
