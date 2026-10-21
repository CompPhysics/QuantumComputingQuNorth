#!/usr/bin/env python3
"""Insert the tagged section EXT:FRONTIER (pre-submission revision of 2026-10-20:
the converged averaged post-selection map and the joint operating-point search
with the exact channel; code from scripts/frontier_checks.py, figure 24) before
the Summary cell of notebooks/ResolutionRefPairing.ipynb.  Built with nbformat;
idempotent.
  python3 scripts/extend_notebook_frontier.py
  PYTHONPATH=src python3 scripts/run_new_cells.py --mark "<!-- EXT:FRONTIER -->"
  make check    # figure count 23 -> 24
"""
import pathlib, re
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
SCRIPTS = ROOT / "scripts"
MARK = "<!-- EXT:FRONTIER -->"

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
## 22&nbsp; The digital filter as an exact averaged channel, and a joint operating-point search

Two revisions follow from the pre-submission assessment of 20 October 2026
(`paper/pairing_pre_submission_critical_assessment.md`).

**(a) A converged averaged map.**  The long-run state of the digital filter is the dominant
eigenoperator of $\mathcal K(\rho)=\mathbb E_t[\tilde K_t\rho\tilde K_t^\dagger]$.  The
fixed-step product formula uses $n(t)=\lceil|t|/\delta t\rceil$ steps, so the integrand
changes discontinuously at every $|t|=n\,\delta t$ and a *global* Gauss–Hermite rule (60–80
nodes, used in Sections 19–20) is not converged for it: at sixteen qubits it even produced a
spurious dominant state with zero ground-state overlap.  `pl.averaged_map_floor` now integrates
the half-Gaussian piecewise between the boundaries (Gauss–Legendre, refinement $m=6\to12$ as the
check) and reports the eigen-equation residual.  The corrected floors are $1.1\times10^{-4}$,
$1.9\times10^{-3}$ and $0.15$ at $k=4$ for $\delta t=0.25,0.5,1$ and, at $k=8$, $g=1$,
$2.3\times10^{-3}$ with eigenvalue $0.79$ at $\delta t=0.25$ (rejected because the acceptance
decays by 21% per cycle and the floor lies above every target, *not* because the ground state
is lost), $4.8\times10^{-5}$ at $\delta t=0.1$ and $3.1\times10^{-6}$ at $0.05$.

**(b) The exact channel and a joint search.**  Since the cycle times are independent, the
ensemble of accepted runs after $M$ cycles is $\mathcal K^M(\rho_{\rm in})/\mathrm{Tr}$, the
acceptance is the trace, and the expected cost per accepted preparation of Eq. (cost) of the
paper is deterministic --- no schedule sampling is needed.  We re-evaluate every input of the
matched comparisons with it and search the two adiabatic routes *jointly* over the time $T$,
the shift $\mu_{\rm buf}$, the step count $n_s$ and the cycle number $M=0,\ldots,6$, the
unfiltered points included, at $k=4$ ($g=1$) and $k=8$ ($g=1,2$).  The minima per route and
target are the table `tab:frontier` of the paper; figure 24 shows the sixteen-qubit frontier
with the ADAPT states included.
"""

cells = [("markdown", md), ("code", cell_from_script("frontier_checks.py"))]
nb = nbf.read(NB, as_version=4)
nb.cells = [c for c in nb.cells if MARK not in c.source]      # idempotent
isum = next(i for i, c in enumerate(nb.cells)
            if c.cell_type == "markdown" and "Summary: the complete EIGEN-Q pipeline" in c.source)
new = [nbf.v4.new_markdown_cell(s) if k == "markdown" else nbf.v4.new_code_cell(s) for k, s in cells]
nb.cells = nb.cells[:isum] + new + nb.cells[isum:]
nbf.write(nb, NB)
print(f"inserted {len(new)} cells before cell {isum}; notebook now has {len(nb.cells)} cells")
