NB     := notebooks/ResolutionRefPairing.ipynb
TEX    := paper/eigenq_pairing.tex
FIGS   := paper/figs
TEXPH  := paper_ph/eigenq_pairing_ph.tex
FIGSPH := paper_ph/figs
PHFIGS := fig12 fig13 fig14 fig15 fig16 fig18
SKILLS := .claude/skills
PY     := python3

.PHONY: install test notebook check figures figures-ph paper paper-ph papers all clean

install:
	$(PY) -m pip install -e ".[notebook,dev]"

test:
	$(PY) -m pytest -q

notebook:
	$(PY) $(SKILLS)/notebook-builder/scripts/execute_notebook.py $(NB) 1800

check:
	$(PY) $(SKILLS)/notebook-builder/scripts/check_notebook.py $(NB) 18

figures:
	$(PY) $(SKILLS)/physics-paper/scripts/extract_figures.py $(NB) $(FIGS)
	PYTHONPATH=src $(PY) scripts/bound_figure_pairing.py

# the particle-hole paper reuses the notebook figures of the EXT:PH/EXT:BOUND sections
figures-ph: figures
	mkdir -p $(FIGSPH)
	for f in $(PHFIGS); do cp $(FIGS)/$$f.png $(FIGSPH)/; done

paper: figures
	bash $(SKILLS)/physics-paper/scripts/build_paper.sh $(TEX)

paper-ph: figures-ph
	bash $(SKILLS)/physics-paper/scripts/build_paper.sh $(TEXPH)

papers: paper paper-ph

all: test notebook check papers

clean:
	rm -rf paper/*.aux paper/*.log paper/*.out paper_ph/*.aux paper_ph/*.log paper_ph/*.out __pycache__ src/*.egg-info
