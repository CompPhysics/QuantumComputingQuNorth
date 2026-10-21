#!/usr/bin/env python3
"""Execute a notebook's code cells in-process (stdout + figures captured into the
notebook, as run_new_cells.py does), with a dill session checkpoint so that a
failing late cell can be fixed and re-run WITHOUT repeating the hours of the
earlier cells.  The result is the same as a clean kernel run (same code, same
seeds); `make notebook-hubbard` remains the clean reference run.
  python3 scripts/run_notebook_checkpoint.py NB --checkpoint-after K --session S.pkl
      executes cells 0..end, dumping the session after code cell index K (the index
      counts CODE cells); on an exception the notebook is written with the outputs
      so far and the session file (if K was reached) is kept.
  python3 scripts/run_notebook_checkpoint.py NB --resume-from K+1 --session S.pkl [--source SRC_NB]
      loads the session and executes code cells K+1..end, taking the code of those
      cells from SRC_NB (default: NB itself, e.g. after a rebuild by the builder
      script) while the outputs of cells 0..K are kept from NB.
"""
import sys, io, base64, pathlib, contextlib, argparse, traceback
import nbformat as nbf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import dill

ap = argparse.ArgumentParser()
ap.add_argument("nb"); ap.add_argument("--checkpoint-after", type=int, default=None)
ap.add_argument("--resume-from", type=int, default=None); ap.add_argument("--session", default=None)
ap.add_argument("--source", default=None); ap.add_argument("--cwd", default=None)
args = ap.parse_args()
NB = pathlib.Path(args.nb)
nb = nbf.read(NB, as_version=4)
src_nb = nbf.read(args.source, as_version=4) if args.source else nb
codes = [i for i, c in enumerate(nb.cells) if c.cell_type == "code"]
src_codes = [i for i, c in enumerate(src_nb.cells) if c.cell_type == "code"]
assert len(codes) == len(src_codes), "source and target notebooks must have the same number of code cells"
if args.cwd:
    import os; os.chdir(args.cwd)
figs_b64 = []
def grab_show(*a, **k):
    for num in plt.get_fignums():
        buf = io.BytesIO(); plt.figure(num).savefig(buf, format="png", dpi=110, bbox_inches="tight")
        figs_b64.append(base64.b64encode(buf.getvalue()).decode())
    plt.close("all")
plt.show = grab_show
start = 0
dill.settings["recurse"] = True
if args.resume_from is not None:
    with open(args.session, "rb") as f:
        ns = dill.load(f)
    start = args.resume_from
else:
    ns = {}
for j, (idx, sidx) in enumerate(zip(codes, src_codes)):
    if j < start:
        continue
    src = src_nb.cells[sidx].source
    nb.cells[idx].source = src
    figs_b64 = []; buf = io.StringIO(); err = None
    print(f"--- code cell {j} (nb index {idx})", flush=True)
    try:
        with contextlib.redirect_stdout(buf):
            exec(compile(src, f"<cell{j}>", "exec"), ns)
        grab_show()
    except Exception:
        err = traceback.format_exc()
    outs = []
    if buf.getvalue():
        outs.append(nbf.v4.new_output("stream", name="stdout", text=buf.getvalue()))
    for b64 in figs_b64:
        outs.append(nbf.v4.new_output("display_data", data={"image/png": b64}))
    if err:
        outs.append(nbf.v4.new_output("error", ename="Exception", evalue=err.strip().split("\n")[-1], traceback=err.split("\n")))
    nb.cells[idx]["outputs"] = outs; nb.cells[idx]["execution_count"] = j + 1
    nbf.write(nb, NB)
    print(buf.getvalue()[-2000:], flush=True)
    if err:
        print(err, flush=True); print(f"FAILED at code cell {j}; notebook written with outputs so far", flush=True); sys.exit(1)
    if args.checkpoint_after is not None and j == args.checkpoint_after and args.session:
        with open(args.session, "wb") as f:
            dill.dump(ns, f)
        print(f"session checkpoint written after code cell {j}: {args.session}", flush=True)
print("notebook executed and written", flush=True)
