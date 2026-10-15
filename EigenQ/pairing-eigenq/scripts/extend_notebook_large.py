#!/usr/bin/env python3
"""Insert the tagged section EXT:LARGE (a larger model space: N = 8 particles on
k = 8 levels, sixteen qubits, seniority zero; code from
scripts/large_space_checks.py, figures 21-22) before the Summary cell of
notebooks/ResolutionRefPairing.ipynb.  Built with nbformat; idempotent.
  python3 scripts/extend_notebook_large.py
  PYTHONPATH=src python3 scripts/run_new_cells.py --mark "<!-- EXT:LARGE -->"
  make check    # figure count 20 -> 22
"""
import pathlib, re
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
SCRIPTS = ROOT / "scripts"
MARK = "<!-- EXT:LARGE -->"

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
## 20&nbsp; A larger model space: four pairs on eight levels ($N=8$, $k=8$, sixteen qubits)

Does a larger, more correlated model space favor refinement + rodeo?  We repeat the matched
comparison of Section 19 for $N=8$ particles on $k=8$ doubly degenerate levels, staying in the
seniority-zero sector of the target.  Every stage of the pipeline conserves seniority (the
Hamiltonian, the embedding, the path, the exactly exponentiated pair-hopping blocks and the
filter), so the whole computation lives in the pair basis of `pairinglib.seniority` --- the exact
seniority-zero block, $\binom{8}{4}=70$ states instead of the $12\,870$ determinants of the
$N=8$ sector --- which is first validated against the full-sector library at $k=4$, $N=4$ (block
identity for $H$ and for every Trotter term, and every anchor of the paper).  The gate counts
remain those of the $2k$-qubit Jordan--Wigner compilation of paper I (a first-order step costs
$1\,360$ CNOTs at $k=8$, the coarse term of the $4\to8$ chain is a phase conditioned on the
eight appended qubits, $254$ CNOTs with the Gray-code convention of Section 17); the
one-qubit-per-level pair encoding is given alongside.  The study covers (a) the spectrum, the
HF overlap and the restricted pair-double ansatz for $g=0.5$--$2$ and the closed form
$\chi_u(0)=P(k-P)\,g^2/(4b^2)$ of the response norm, (b) exact-slice refinement chains
$k_{\rm low}\to8$ for $k_{\rm low}=4$--$7$ and shifts $\mu_{\rm buf}=0.6,1,2$, the HF ramp and
the coupling dependence of the transports (figure 21), and (c) the circuit level: product-formula
steps, the averaged post-selection map of the digital filter ($\delta t=0.25$ is rejected at
these norms), nine inputs to the same Trotterized filter at $g=1$ and five at $g=2$, with $M=0$,
CNOTs per attempt and per accepted preparation, bootstrap errors and the minimum cost at three
targets (figure 22).
"""

cells = [("markdown", md), ("code", cell_from_script("large_space_checks.py"))]
nb = nbf.read(NB, as_version=4)
nb.cells = [c for c in nb.cells if MARK not in c.source]      # idempotent
isum = next(i for i, c in enumerate(nb.cells)
            if c.cell_type == "markdown" and "Summary: the complete EIGEN-Q pipeline" in c.source)
new = [nbf.v4.new_markdown_cell(s) if k == "markdown" else nbf.v4.new_code_cell(s) for k, s in cells]
nb.cells = nb.cells[:isum] + new + nb.cells[isum:]
nbf.write(nb, NB)
print(f"inserted {len(new)} cells before cell {isum}; notebook now has {len(nb.cells)} cells")
