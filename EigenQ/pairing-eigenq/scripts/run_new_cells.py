#!/usr/bin/env python3
"""Execute the EXT:TROTTER cells of ResolutionRefPairing.ipynb in-process and
attach their real outputs (stdout + PNG figures) to the notebook via nbformat.
Usage: run_new_cells.py [first_idx last_idx]   (indices into the new-cell list;
earlier new cells are always re-executed silently to rebuild kernel state).
The remaining (pre-existing) cells keep their stored outputs; a full clean
re-execution is `make notebook`."""
import sys, io, base64, pathlib, contextlib
import nbformat as nbf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
MARK = "<!-- EXT:TROTTER -->"

nb = nbf.read(NB, as_version=4)
new_code = [i for i, c in enumerate(nb.cells)
            if c.cell_type == "code" and MARK in c.source]
lo, hi = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, len(new_code) - 1)

maxec = max([c.get("execution_count") or 0 for c in nb.cells if c.cell_type == "code"])
ns = {}
figs_b64 = []


def grab_show(*a, **k):
    for num in plt.get_fignums():
        buf = io.BytesIO()
        plt.figure(num).savefig(buf, format="png", dpi=110, bbox_inches="tight")
        figs_b64.append(base64.b64encode(buf.getvalue()).decode())
    plt.close("all")


plt.show = grab_show

for j, idx in enumerate(new_code):
    src = nb.cells[idx].source.replace(MARK, "")
    figs_b64 = []
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        exec(compile(src, f"<newcell{j}>", "exec"), ns)
    grab_show()
    if lo <= j <= hi:
        outs = []
        if buf.getvalue():
            outs.append(nbf.v4.new_output("stream", name="stdout", text=buf.getvalue()))
        for b64 in figs_b64:
            outs.append(nbf.v4.new_output("display_data", data={"image/png": b64}))
        nb.cells[idx]["outputs"] = outs
        nb.cells[idx]["execution_count"] = maxec + j + 1
        print(f"cell {j} (nb index {idx}): {len(outs)} outputs")

nbf.write(nb, NB)
print("notebook updated")
