NB     := notebooks/ResolutionRefPairing.ipynb
TEX    := paper/eigenq_pairing.tex
TEXSM  := paper/eigenq_pairing_supp.tex
FIGS   := paper/figs
TEXPH  := paper_ph/eigenq_pairing_ph.tex
FIGSPH := paper_ph/figs
PHFIGS := fig12 fig13 fig14 fig15 fig16 fig18
NBH    := notebooks/ResolutionRefHubbard.ipynb
TEXH   := paper3/eigenq_hubbard.tex
TEXHSM := paper3/eigenq_hubbard_supp.tex
FIGSH  := paper3/figs
SKILLS := .claude/skills
PY     := python3

.PHONY: install test notebook check figures figures-ph paper paper-ph papers all clean \
        notebook-hubbard check-hubbard figures-hubbard paper-hubbard hubbard

install:
	$(PY) -m pip install -e ".[notebook,dev]"

test:
	$(PY) -m pytest -q

notebook:
	$(PY) $(SKILLS)/notebook-builder/scripts/execute_notebook.py $(NB) 1800

check:
	$(PY) $(SKILLS)/notebook-builder/scripts/check_notebook.py $(NB) 25

figures:
	$(PY) $(SKILLS)/physics-paper/scripts/extract_figures.py $(NB) $(FIGS)
	PYTHONPATH=src $(PY) scripts/bound_figure_pairing.py

# the particle-hole paper reuses the notebook figures of the EXT:PH/EXT:BOUND sections
figures-ph: figures
	mkdir -p $(FIGSPH)
	for f in $(PHFIGS); do cp $(FIGS)/$$f.png $(FIGSPH)/; done

paper: figures
	bash $(SKILLS)/physics-paper/scripts/build_paper.sh $(TEX)
	# the SM resolves main-text references with xr; the filtered copy of the main .aux drops the
	# labels/citations that both documents define (FirstPage, LastPage, LastBibItem, \bibcite)
	grep -v -e FirstPage -e LastPage -e LastBibItem -e bibcite paper/eigenq_pairing.aux > paper/eigenq_pairing_xr.aux
	bash $(SKILLS)/physics-paper/scripts/build_paper.sh $(TEXSM)   # Supplemental Material (xr against eigenq_pairing_xr.aux)

paper-ph: figures-ph
	bash $(SKILLS)/physics-paper/scripts/build_paper.sh $(TEXPH)

papers: paper paper-ph

# ---- paper III: lattice refinement of the Hubbard model (own notebook, built by scripts/build_notebook_hubbard.py)
notebook-hubbard:
	$(PY) scripts/build_notebook_hubbard.py
	$(PY) $(SKILLS)/notebook-builder/scripts/execute_notebook.py $(NBH) 7200

check-hubbard:
	$(PY) scripts/check_notebook_hubbard.py $(NBH) 15

figures-hubbard:
	$(PY) $(SKILLS)/physics-paper/scripts/extract_figures.py $(NBH) $(FIGSH)

paper-hubbard: figures-hubbard
	bash $(SKILLS)/physics-paper/scripts/build_paper.sh $(TEXH)
	# the SM resolves main-text references with xr against a copy of the main .aux without the \bibcite entries
	# (hyperref + xr otherwise fail on citations the SM defines itself)
	grep -v -e bibcite paper3/eigenq_hubbard.aux > paper3/eigenq_hubbard_xr.aux
	bash $(SKILLS)/physics-paper/scripts/build_paper.sh $(TEXHSM)   # Supplemental Material (xr against eigenq_hubbard_xr.aux)

hubbard: notebook-hubbard check-hubbard paper-hubbard

all: test notebook check papers

clean:
	rm -rf paper/*.aux paper/*.log paper/*.out paper_ph/*.aux paper_ph/*.log paper_ph/*.out paper3/*.aux paper3/*.log paper3/*.out paper3/*.bbl paper3/*Notes.bib __pycache__ src/*.egg-info
