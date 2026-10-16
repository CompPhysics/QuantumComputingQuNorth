#!/usr/bin/env python3
"""Insert the tagged section EXT:ADAPT (the ADAPT-VQE comparator; code from
scripts/adapt_checks.py, figure 23) before the Summary cell of
notebooks/ResolutionRefPairing.ipynb.  Built with nbformat; idempotent.
  python3 scripts/extend_notebook_adapt.py
  PYTHONPATH=src python3 scripts/run_new_cells.py --mark "<!-- EXT:ADAPT -->"
  make check    # figure count 22 -> 23
"""
import pathlib, re
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
SCRIPTS = ROOT / "scripts"
MARK = "<!-- EXT:ADAPT -->"

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
## 21&nbsp; The ADAPT-VQE comparator

The matched comparisons of Sections 19 and 20 used a *restricted, single-layer* variational
ansatz (the occupied-pair$\to$virtual-pair rotations, each once).  ADAPT-VQE (Grimsley *et al.*,
2019) builds the ansatz adaptively: at every iteration the pool operator with the largest energy
gradient $|\langle\psi|[H,K]|\psi\rangle|$ is appended --- an operator may be selected again --- and
all amplitudes are re-optimized.  In the seniority-zero pair basis of `pairinglib.seniority` the
natural pools are the $P(k-P)$ occupied$\to$virtual pair rotations of the restricted ansatz and
the *generalized* pool of all $k(k-1)/2$ pair-hopping generators $P_q^\dagger P_p-{\rm h.c.}$;
each generator is a 48-CNOT block in the $2k$-qubit encoding (4 CNOTs in the pair encoding).  We
run both pools for $(k,N)=(4,4)$ and $(6,4)$ at $g=1$ and for $(8,8)$ at $g=1$ and $2$, record the
operators and CNOTs needed to reach the targets of Section 20, and feed the $(8,8)$ ADAPT states
at the first iteration below $10^{-4}$ to the same Trotterized filter (figure 23).
"""

cells = [("markdown", md), ("code", cell_from_script("adapt_checks.py"))]
nb = nbf.read(NB, as_version=4)
nb.cells = [c for c in nb.cells if MARK not in c.source]      # idempotent
isum = next(i for i, c in enumerate(nb.cells)
            if c.cell_type == "markdown" and "Summary: the complete EIGEN-Q pipeline" in c.source)
new = [nbf.v4.new_markdown_cell(s) if k == "markdown" else nbf.v4.new_code_cell(s) for k, s in cells]
nb.cells = nb.cells[:isum] + new + nb.cells[isum:]
nbf.write(nb, NB)
print(f"inserted {len(new)} cells before cell {isum}; notebook now has {len(nb.cells)} cells")
