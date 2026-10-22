#!/usr/bin/env python3
"""Insert the tagged section EXT:REFINE (final revision of 2026-10-21 after
paper/morefeedback.md: convergence of the operating-point search behind
Table tab:frontier; code from scripts/frontier_refine.py, figure 25) before
the Summary cell of notebooks/ResolutionRefPairing.ipynb.  Built with
nbformat; idempotent.
  python3 scripts/extend_notebook_refine.py
  PYTHONPATH=src python3 scripts/run_new_cells.py --mark "<!-- EXT:REFINE -->"
  make check    # figure count 24 -> 25
"""
import pathlib, re
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
SCRIPTS = ROOT / "scripts"
MARK = "<!-- EXT:REFINE -->"

PREAMBLE = """import numpy as np, matplotlib.pyplot as plt
try:
    import pairinglib as pl
except ModuleNotFoundError:
    import sys, pathlib
    sys.path.insert(0, str((pathlib.Path.cwd()/".."/"src").resolve()))
    import pairinglib as pl
"""


def cell_from_script(name):
    src = (SCRIPTS / name).read_text()
    src = re.sub(r'^#!.*\n', '', src)
    src = re.sub(r'^"""[\s\S]*?"""\n', '', src)
    src = re.sub(r'^(import numpy as np|import matplotlib\.pyplot as plt|import pairinglib as pl|sys\.path\.insert.*)\n',
                 '', src, flags=re.M)
    return "# " + MARK + "\n" + PREAMBLE + src


md = MARK + r"""
---
## 23&nbsp; Convergence of the operating-point search

The joint search of Section 22 used a coarse control grid ($T$ in steps of $2$, three shifts,
eight step counts), and the final pre-submission assessment of 21 October 2026
(`paper/morefeedback.md`) showed by constructive examples --- e.g. the unfiltered HF ramp at
$T=4.5$, $n_s=11$ ($1.63\times10^{-3}$ for $6\,512$ CNOTs) against the $9\,472$ CNOTs of the
coarse grid at $k=4$ --- that the quoted minima were minima *of that grid*, too coarse to carry
percent-level rankings of the two adiabatic routes.  Here both routes are searched on a fine
grid of the time $T$ (steps of $0.25$ or $0.5$) and of the shift $\mu_{\rm buf}$ (nine or ten
values), with the step count $n_s$ on a geometric list (ratio $1.15$) and then refined to the
integer by bracketing and bisection at every target and every cycle number $M=0,\ldots,6$, with
the exact averaged channel of the digital filter (now applied through its cached
$d^2\times d^2$ superoperator, `AveragedFilterMap.superoperator`, twenty times faster than the
node loop at $d=70$).  Three search levels are reported per route and target --- the coarse grid
of Section 22 (L0), the fine grid (L1) and the integer refinement (L2) --- together with the
number of distinct controls whose refined cost lies within $10\%$ of the minimum, i.e. how flat
the frontier is, and the ratio of the two routes.  The ADAPT states are admitted at every
prefix (the four-operator prefix of the $(4,4)$ run, $192$ CNOTs, was missing from the Section-22
search).  Figure 25 draws the cost frontier at $k=8$: the best cost found against the target
infidelity, unfiltered and with up to six cycles, for every route.  These are the numbers of
Table `tab:frontier` and of Table S6 of paper I.
"""

cells = [("markdown", md), ("code", cell_from_script("frontier_refine.py"))]
nb = nbf.read(NB, as_version=4)
nb.cells = [c for c in nb.cells if MARK not in c.source]      # idempotent
isum = next(i for i, c in enumerate(nb.cells)
            if c.cell_type == "markdown" and "Summary: the complete EIGEN-Q pipeline" in c.source)
new = [nbf.v4.new_markdown_cell(s) if k == "markdown" else nbf.v4.new_code_cell(s) for k, s in cells]
nb.cells = nb.cells[:isum] + new + nb.cells[isum:]
nbf.write(nb, NB)
print(f"inserted {len(new)} cells before cell {isum}; notebook now has {len(nb.cells)} cells")
