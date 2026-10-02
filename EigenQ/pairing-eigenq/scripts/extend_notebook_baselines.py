#!/usr/bin/env python3
"""Insert the matched-observable baseline study (paper Sec. V.C), the
size/particle-number extension and the observables table (paper Sec. V.E) into
ResolutionRefPairing.ipynb, before the summary cell.  The cell code is read
from scripts/baselines_rodeo.py and scripts/sizes_observables.py, so the
notebook and the standalone scripts cannot drift apart.  Built with nbformat
(never hand-edit the .ipynb JSON); idempotent: cells tagged EXT:BASELINES are
replaced.
Run from the repo root:
  python3 scripts/extend_notebook_baselines.py
  python3 scripts/run_new_cells.py --mark "<!-- EXT:BASELINES -->"
  make figures   # -> paper/figs/fig17.png (the 17th figure in cell order)
"""
import pathlib, re
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
SCRIPTS = ROOT / "scripts"
MARK = "<!-- EXT:BASELINES -->"

PREAMBLE = """import numpy as np, matplotlib.pyplot as plt
try:
    import pairinglib as pl
except ModuleNotFoundError:
    import sys, pathlib
    sys.path.insert(0, str((pathlib.Path.cwd()/".."/"src").resolve()))
    import pairinglib as pl
"""


def cell_from_script(name):
    """Strip the module docstring, the sys.path/import boilerplate and the
    command-line handling of a script; what remains is the notebook cell."""
    src = (SCRIPTS / name).read_text()
    src = re.sub(r'^#!.*\n', '', src)
    src = re.sub(r'^"""[\s\S]*?"""\n', '', src)                   # module docstring
    src = re.sub(r'^(import numpy as np|import pairinglib as pl|sys\.path\.insert.*)\n',
                 '', src, flags=re.M)
    src = src.replace('import sys, os\n', 'import os\n')
    src = src.replace('OUT = sys.argv[1] if len(sys.argv) > 1 else None', 'OUT = None')
    return "# " + MARK + "\n" + PREAMBLE + src


md_base = MARK + r"""
---
## 17&nbsp; Which input should the filter get?  Matched-observable baselines

The value of the refinement stage is decided by comparing *inputs to the same
filter* on the *same observable*.  An energy error is not an infidelity: for a
normalised state $\langle H\rangle-E_0=\sum_{j>0}|c_j|^2(E_j-E_0)\ge\Delta(1-F)$,
so the UCCSD energy error $1.4\times10^{-3}$ at $k=4$ caps its infidelity at
$7.9\times10^{-4}$; the actual value is $1.8\times10^{-4}$.  Below, one and the
same Trotterised rodeo filter (order 2, $\delta t=0.25$, $\sigma=4$, 40 common
time sets) is fed (i) the prolonged determinant, (ii) the full-space gate-level
UCCSD state, (iii) the refined state at $T^*=16$ and (iv) at $T=30$.  We report
the acceptance-weighted ensemble infidelity, the acceptance, and the CNOTs per
accepted preparation including early termination and the embedded coarse term
$A=P(H_{\rm low}-\mu)P^\dagger$ (a 4-qubit conditional phase, 14 CNOTs per
application for the $2\to4$ chain).  The last block repeats the UCCSD+rodeo
comparison with the particle-hole term (paper Table V).
"""

md_sizes = MARK + r"""
### 17.1&nbsp; Other sizes and particle numbers; level occupations

The same comparison with exact-slice refinement and exact controlled evolution
for $(k,N)=(4,4),(6,4),(4,6),(6,6)$ (paper Table III), followed by the level
occupation numbers $n_p$ and the correlation energy of the prepared states at
$k=4$, $N=4$ (paper Table IV).
"""

cells = [("markdown", md_base), ("code", cell_from_script("baselines_rodeo.py")),
         ("markdown", md_sizes), ("code", cell_from_script("sizes_observables.py"))]

nb = nbf.read(NB, as_version=4)
nb.cells = [c for c in nb.cells if MARK not in c.source]      # idempotent
isum = next(i for i, c in enumerate(nb.cells)
            if c.cell_type == "markdown" and "Summary: the complete EIGEN-Q pipeline" in c.source)
new = [nbf.v4.new_markdown_cell(s) if k == "markdown" else nbf.v4.new_code_cell(s)
       for k, s in cells]
nb.cells = nb.cells[:isum] + new + nb.cells[isum:]
nbf.write(nb, NB)
print(f"inserted {len(new)} cells before cell {isum}; notebook now has {len(nb.cells)} cells")
