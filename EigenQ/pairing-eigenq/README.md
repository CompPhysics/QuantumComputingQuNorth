# pairing-eigenq

Reproducible pipeline for **resolution refinement + the rodeo algorithm** applied
to the constant-pairing Hamiltonian — a controlled demonstration of the EIGEN-Q
eigenstate-preparation workflow, producing Jupyter notebooks and three LaTeX
articles: paper I on the pure pairing model, paper II on the pairing model
with a seniority-breaking particle-hole interaction, and paper III on the
lattice refinement of the Hubbard model.

```
pairing-eigenq/
├── CLAUDE.md                 # project memory for Claude Code
├── .claude/skills/           # 4 Agent Skills (pairing-model, vqe-circuits,
│                             #   notebook-builder, physics-paper)
├── src/pairinglib/           # tested Python package (single source of truth)
├── notebooks/                # ResolutionRefPairing.ipynb, ResolutionRefHubbard.ipynb (build artifacts)
├── paper/                    # paper I: eigenq_pairing.tex (main text) + eigenq_pairing_supp.tex (Supplemental Material) + figs/
├── paper_ph/                 # paper II: eigenq_pairing_ph.tex + figs/ (copied from paper/figs)
├── paper3/                   # paper III: eigenq_hubbard.tex + eigenq_hubbard_supp.tex (quantumarticle) + figs/, Hubbard.ipynb (ED reference)
├── tests/                    # pytest: API + benchmark anchors
├── scripts/pipeline.py       # one-command orchestration
├── Makefile                  # install / test / notebook / check / figures / paper / paper-ph / all
└── .github/workflows/ci.yml  # runs tests + builds notebook + paper
```

## Quickstart
```bash
make install        # pip install -e . + jupyter/latex deps
make test           # library + benchmark anchors
make all            # notebook -> check -> figures -> paper I (main text + Supplemental Material) and paper II PDFs
make paper-ph       # paper II only (copies the V_ph figures first)
make hubbard        # paper III: Hubbard notebook (~2.5 h) -> check -> figures -> main text + SM PDFs (quantumarticle)
```

## Working with Claude Code
Open the repo with Claude Code; the skills in `.claude/skills/` load
automatically and `CLAUDE.md` provides project memory. Example requests:
- "Add a Trotterized controlled-evolution variant to the rodeo chapter and
  quantify the Trotter error." (uses `vqe-circuits` + `notebook-builder`)
- "Regenerate the paper figures and rebuild the PDF." (uses `physics-paper`)
- "Extend the four-particle analysis to N=6." (uses `pairing-model` + library)
- "Run the pipeline with the particle-hole term at f=0.3." (pass `f=` to the
  library calls; see `scripts/extend_notebook_ph.py` for the full section)

Authors (alphabetical): Bogner, Glittum, Hergert, Hjorth-Jensen, Lange, LaRose,
Lee, Massel, Sadia. See `paper/eigenq_pairing.tex` for affiliations.

## Reproducibility and data availability
Every number and figure in paper I and its Supplemental Material is an output
of `notebooks/ResolutionRefPairing.ipynb` (22 figures; `make check` asserts the
benchmark anchors), built on the tested package `src/pairinglib`; paper III's
numbers come from `notebooks/ResolutionRefHubbard.ipynb`.  A clean run of the
pairing notebook takes about 35 minutes on a laptop.  The figure files used by
the papers are the extracted notebook outputs in `paper/figs/` and
`paper3/figs/`.  On publication the repository is archived on Zenodo; cite the
archived version (see `CITATION.cff`) together with the article.
