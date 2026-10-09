#!/usr/bin/env python3
"""Insert the 'sufficient conditions versus practice' study (paper Sec. IV.E,
Eq. eq:Tbound) into ResolutionRefPairing.ipynb, before the summary cell.  The
cell code is read from scripts/adiabatic_bound.py, so the notebook and the
standalone script cannot drift apart.  Built with nbformat (never hand-edit the
.ipynb JSON); idempotent: cells tagged EXT:BOUND are replaced.
Run from the repo root:
  python3 scripts/extend_notebook_bound.py
  python3 scripts/run_new_cells.py --mark "<!-- EXT:BOUND -->"
  make figures   # -> paper/figs/fig18.png (the 18th figure in cell order)
"""
import pathlib, re
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
SCRIPTS = ROOT / "scripts"
MARK = "<!-- EXT:BOUND -->"

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
    src = re.sub(r'^"""[\s\S]*?"""\n', '', src)                   # module docstring
    src = re.sub(r'^(import numpy as np|import matplotlib\.pyplot as plt|import pairinglib as pl|sys\.path\.insert.*)\n',
                 '', src, flags=re.M)
    src = src.replace('import sys, os\n', '')
    return "# " + MARK + "\n" + PREAMBLE + src


md = MARK + r"""
---
## 18&nbsp; Sufficient conditions versus practice along the refinement path

The interpolation $H(s)=\cos^2(\pi s/2)\,A+\sin^2(\pi s/2)\,B$ with
$A=P(H_{\rm low}-\mu)P^\dagger$ and $B=H_{\rm high}-\mu$ is the $\sin^2$ schedule
for which the adiabatic bound of Jansen, Ruskai and Seiler has no boundary terms.
The nucleation paper (Sadia *et al.*, arXiv:2609.37527) turns it into explicit
sufficient conditions for a Trotterised adiabatic stage: a total time
$T\ge\epsilon^{-1/2}\int_0^1[\|\ddot H\|/\Delta^2+7\|\dot H\|^2/\Delta^3]\,ds$
guarantees a final infidelity $\le\epsilon$, and
$n_s\ge\eta^{-1}(T\Lambda_1/4+T^2\Lambda_2/2)$ first-order Trotter steps keep the
operator error of the digitised sweep below $\eta$, with $\Lambda_1=\max_s\|\dot H\|$
and $\Lambda_2$ the largest total commutator norm between the Trotter groups
(here $\sum_{i<j}\|[H_i,H_j]\|$ over the groups of the gate compilation,
$[A,\,D-\mu,\,h_{pq},\ldots,V_{\rm ph}\ {\rm pieces}]$).  Along our path
$\dot H=\tfrac{\pi}{2}\sin(\pi s)(B-A)$ and $\ddot H=\tfrac{\pi^2}{2}\cos(\pi s)(B-A)$,
so everything reduces to $\|B-A\|$ and the instantaneous gap $\Delta(s)$; all norms
are spectral norms in the $N$-particle sector, where the evolution lives.  We
evaluate both conditions for the $k=2\to4$ chain at $N=4$, $g=1$, for the pure
pairing model and the three particle-hole strengths of Section 16, and compare
them with the empirical $T^*$ (overlap $0.99$) and with the number of product-formula
steps that reproduce the exact-slice overlap at $T^*$ to within the adiabatic error
of the path itself.
"""

cells = [("markdown", md), ("code", cell_from_script("adiabatic_bound.py"))]

nb = nbf.read(NB, as_version=4)
nb.cells = [c for c in nb.cells if MARK not in c.source]      # idempotent
isum = next(i for i, c in enumerate(nb.cells)
            if c.cell_type == "markdown" and "Summary: the complete EIGEN-Q pipeline" in c.source)
new = [nbf.v4.new_markdown_cell(s) if k == "markdown" else nbf.v4.new_code_cell(s)
       for k, s in cells]
nb.cells = nb.cells[:isum] + new + nb.cells[isum:]
nbf.write(nb, NB)
print(f"inserted {len(new)} cells before cell {isum}; notebook now has {len(nb.cells)} cells")
