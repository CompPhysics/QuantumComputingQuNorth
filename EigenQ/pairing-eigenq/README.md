# pairing-eigenq

Reproducible pipeline for **resolution refinement + the rodeo algorithm** applied
to the constant-pairing Hamiltonian — a controlled demonstration of the EIGEN-Q
eigenstate-preparation workflow, producing a Jupyter notebook and two LaTeX
articles: paper I on the pure pairing model and paper II on the pairing model
with a seniority-breaking particle-hole interaction.

```
pairing-eigenq/
├── CLAUDE.md                 # project memory for Claude Code
├── .claude/skills/           # 4 Agent Skills (pairing-model, vqe-circuits,
│                             #   notebook-builder, physics-paper)
├── src/pairinglib/           # tested Python package (single source of truth)
├── notebooks/                # ResolutionRefPairing.ipynb (build artifact)
├── paper/                    # paper I: eigenq_pairing.tex + figs/ (generated)
├── paper_ph/                 # paper II: eigenq_pairing_ph.tex + figs/ (copied from paper/figs)
├── tests/                    # pytest: API + benchmark anchors
├── scripts/pipeline.py       # one-command orchestration
├── Makefile                  # install / test / notebook / check / figures / paper / paper-ph / all
└── .github/workflows/ci.yml  # runs tests + builds notebook + paper
```

## Quickstart
```bash
make install        # pip install -e . + jupyter/latex deps
make test           # library + benchmark anchors
make all            # notebook -> check -> figures -> both paper PDFs
make paper-ph       # paper II only (copies the V_ph figures first)
```

## Working with Claude Code
Open the repo with Claude Code; the skills in `.claude/skills/` load
automatically and `CLAUDE.md` provides project memory. Example requests:
- "Add a Trotterised controlled-evolution variant to the rodeo chapter and
  quantify the Trotter error." (uses `vqe-circuits` + `notebook-builder`)
- "Regenerate the paper figures and rebuild the PDF." (uses `physics-paper`)
- "Extend the four-particle analysis to N=6." (uses `pairing-model` + library)
- "Run the pipeline with the particle-hole term at f=0.3." (pass `f=` to the
  library calls; see `scripts/extend_notebook_ph.py` for the full section)

Authors (alphabetical): Bogner, Glittum, Hergert, Hjorth-Jensen, Lange, LaRose,
Lee, Massel, Sadia. See `paper/eigenq_pairing.tex` for affiliations.
