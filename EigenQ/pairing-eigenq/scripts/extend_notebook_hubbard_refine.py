#!/usr/bin/env python3
"""Insert the tagged Section 14 (EXT:HREFINE: convergence of the operating-point
search, code from scripts/hubbard_refine_checks.py, figure 15) into the EXECUTED
notebooks/ResolutionRefHubbard.ipynb before its Summary cell, and renumber the
Summary to 15 — the same cells that scripts/build_notebook_hubbard.py now
emits on a full rebuild, so the executed notebook can be extended without the
six-hour rerun.  Built with nbformat; idempotent.
  python3 scripts/extend_notebook_hubbard_refine.py
  PYTHONPATH=src python3 scripts/run_new_cells.py --nb notebooks/ResolutionRefHubbard.ipynb --mark "<!-- EXT:HREFINE -->"
  make check-hubbard    # figure count 14 -> 15
"""
import pathlib, sys
import nbformat as nbf
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import build_notebook_hubbard as bh          # for MD_HREFINE / _cell_from_script / MARK (the module only defines cells on import)

NB = ROOT / "notebooks/ResolutionRefHubbard.ipynb"
nb = nbf.read(NB, as_version=4)
nb.cells = [c for c in nb.cells if bh.MARK_HREFINE not in c.source]          # idempotent
isum = next(i for i, c in enumerate(nb.cells) if c.cell_type == "markdown" and "Summary" in c.source and c.source.lstrip().startswith("---\n## 1"))
nb.cells[isum].source = nb.cells[isum].source.replace("## 14&nbsp; Summary", "## 15&nbsp; Summary")
new = [nbf.v4.new_markdown_cell(bh.MD_HREFINE), nbf.v4.new_code_cell(bh._cell_from_script("hubbard_refine_checks.py"))]
nb.cells = nb.cells[:isum] + new + nb.cells[isum:]
nbf.write(nb, NB)
print(f"inserted {len(new)} cells before cell {isum}; notebook now has {len(nb.cells)} cells")
