#!/usr/bin/env python3
"""Refresh, IN PLACE, the code cell of a tagged notebook section whose code is a
copy of a script (EXT:ASSESS <- assessment_checks.py, EXT:LARGE <-
large_space_checks.py, EXT:ADAPT <- adapt_checks.py, EXT:FRONTIER <-
frontier_checks.py) after the script changed, keeping the cell position and
the section numbering; then re-execute it with run_new_cells.py.
  python3 scripts/refresh_notebook_cells.py EXT:ASSESS EXT:LARGE
  PYTHONPATH=src python3 scripts/run_new_cells.py --mark "<!-- EXT:ASSESS -->"
"""
import sys, pathlib, re
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
SCRIPTS = ROOT / "scripts"
SOURCES = {"EXT:ASSESS": "assessment_checks.py", "EXT:LARGE": "large_space_checks.py",
           "EXT:ADAPT": "adapt_checks.py", "EXT:FRONTIER": "frontier_checks.py"}
PREAMBLE = """import numpy as np, matplotlib.pyplot as plt
try:
    import pairinglib as pl
except ModuleNotFoundError:
    import sys, pathlib
    sys.path.insert(0, str((pathlib.Path.cwd()/".."/"src").resolve()))
    import pairinglib as pl
"""


def cell_from_script(name, mark):
    src = (SCRIPTS / name).read_text()
    src = re.sub(r'^#!.*\n', '', src)
    src = re.sub(r'^"""[\s\S]*?"""\n', '', src)
    src = re.sub(r'^(import numpy as np|import matplotlib\.pyplot as plt|import pairinglib as pl|sys\.path\.insert.*)\n',
                 '', src, flags=re.M)
    return "# <!-- " + mark + " -->\n" + PREAMBLE + src


nb = nbf.read(NB, as_version=4)
for mark in sys.argv[1:]:
    tag = f"<!-- {mark} -->"
    idx = [i for i, c in enumerate(nb.cells) if c.cell_type == "code" and tag in c.source]
    assert len(idx) == 1, (mark, idx)
    nb.cells[idx[0]].source = cell_from_script(SOURCES[mark], mark)
    nb.cells[idx[0]].outputs = []
    print(f"refreshed {mark} (cell {idx[0]}) from scripts/{SOURCES[mark]}")
nbf.write(nb, NB)
