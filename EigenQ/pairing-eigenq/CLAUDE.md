# CLAUDE.md — project memory

This repository produces a reproducible Jupyter notebook and a LaTeX article on
**resolution refinement + the rodeo algorithm for eigenstate preparation**,
demonstrated on the constant-pairing Hamiltonian (the EIGEN-Q pipeline).

## Source of truth
- Physics lives in `src/pairinglib/` (a tested package). **Import it**; never
  re-define library functions inside notebooks or scripts.
- The notebook (`notebooks/ResolutionRefPairing.ipynb`) is a build artifact.
- There are TWO papers, split on 2026-10-12 (the pre-split manuscript is
  `paper/eigenq_pairing_v3_backup.tex`):
  * Paper I, `paper/eigenq_pairing.tex` — the pure (seniority-conserving)
    pairing model only; cites paper II as `\cite{PHcompanion}`.
  * Paper II, `paper_ph/eigenq_pairing_ph.tex` — the pairing model plus the
    seniority-breaking particle-hole term V_ph (model, refinement, bounds,
    rodeo, circuit budget at f=0.05/0.2/0.5); self-contained, with a short
    recap of the pipeline and of the f=0 reference numbers; cites paper I as
    `\cite{PairingCompanion}`.  Its figures are copies of notebook figures
    12-16 and 18 (`make figures-ph`); paper I uses `figs/fig_bound_pairing.png`
    (f=0-only redraw of fig18, `scripts/bound_figure_pairing.py`).
  Both are generated from the notebook's results and figures. Every number in
  either paper must trace to a notebook output.  Keep the two consistent:
  the f=0 row of paper II's tables is paper I's result.

## Skills (in `.claude/skills/`) — read the matching one before working
- `pairing-model`  — Hamiltonian conventions, encoding, benchmark anchors, API.
- `vqe-circuits`   — UCCSD/HEA, prolongation, refinement, and the rodeo algorithm.
- `notebook-builder` — assemble/execute/validate the notebook (scripts included).
- `physics-paper`  — extract figures, write, and compile the article.

## Pipeline (also wired into the Makefile)
1. `make install`   — editable install of `pairinglib` + deps.
2. `make test`      — pytest: library API + benchmark anchors.
3. `make notebook`  — execute the notebook in place (generous timeout).
4. `make check`     — validate the executed notebook (no errors, figures, anchors).
5. `make figures`   — extract notebook PNGs into `paper/figs/`.
6. `make paper`     — compile `paper/eigenq_pairing.tex` to PDF (paper I).
7. `make paper-ph`  — copy the V_ph figures and compile paper II; `make papers` builds both.
8. `make all`       — the whole chain end to end.

## Conventions
- Python: import `pairinglib as pl`. Keep new physics in the package + a test.
- Notebooks: built with `notebook-builder` scripts (nbformat), never hand-edited JSON.
- LaTeX: `article` class + `authblk`; authors ALPHABETICAL by surname.
- Benchmark anchors (N=4, g=1): FCI(k=3)=0.794697, FCI(k=4)=0.635548; rodeo
  acceptance -> p. With the particle-hole term (k=4): f=0.05 -> 0.45058234,
  f=0.2 -> -0.18348455, f=0.5 -> -1.69173670 (book Chapter 4). CI and
  `make check` assert these.
- The seniority-breaking particle-hole term `V_ph` (strength `f`) is switched on
  with the `f=` keyword everywhere (`phmodel.py`); `f=0` is the pure pairing model.
- Notebook extensions are tagged sections (`scripts/extend_notebook_*.py` +
  `scripts/run_new_cells.py --mark`); see the notebook-builder skill.  Sections:
  EXT:TROTTER, EXT:PH, EXT:BASELINES (matched-observable baselines, sizes,
  observables; cell code is read from `scripts/baselines_rodeo.py` and
  `scripts/sizes_observables.py`, which also run standalone), EXT:BOUND
  (sufficient adiabatic-time / Trotter-step conditions of the nucleation paper,
  arXiv:2609.37527, evaluated on the 2->4 path; code from
  `scripts/adiabatic_bound.py`; paper Sec. IV.E, Table tab:bound, fig18).
- Compare like with like: an energy error is NOT an infidelity
  (<H>-E0 >= gap*(1-F)).  UCCSD k=4, g=1 infidelities: 1.78e-4 (f=0),
  1.01e-3 (f=0.2), 3.72e-3 (f=0.5).  Every comparison with the pipeline uses the
  acceptance-weighted ensemble infidelity of the accepted runs.
- Gate budgets must include the embedded coarse term A = P(H_low-mu)P^T of the
  refinement path (2->4: 4-qubit conditional phase, 14 CNOTs/application;
  3->4: doubly-conditioned k=3 step, ~350 CNOTs).  `cnot_counts` counts the
  high-space step only.
- Bounds vs. practice (k=2->4, N=4, g=1): T_bound(eps=1e-2)=5.1e3 vs T*=16;
  first-order step bound 3.4e4 at T* vs 160 (order 1) / 40 (order 2) in practice.
- Review history: `paper/revisions.txt` (external review) and
  `paper/response_to_revisions.md` (what was verified/changed/open).

## Definition of done
A change is done only when `make test` and `make check` pass and (if a paper
changed) `make paper` / `make paper-ph` builds without unresolved references.
