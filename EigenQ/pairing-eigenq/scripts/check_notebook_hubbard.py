#!/usr/bin/env python3
"""Validate the executed Hubbard notebook (paper III): no errors, >= min figures,
anchor numbers present.  Usage: check_notebook_hubbard.py NB [min_figs]"""
import sys, json
nb = json.load(open(sys.argv[1])); min_figs = int(sys.argv[2]) if len(sys.argv) > 2 else 8
errs, figs, text = [], 0, []
for i, c in enumerate(nb["cells"]):
    if c["cell_type"] != "code": continue
    for o in c.get("outputs", []):
        if o.get("output_type") == "error": errs.append((i, o.get("ename")))
        if o.get("output_type") == "display_data" and "image/png" in o.get("data", {}): figs += 1
        if "text" in o: text.append("".join(o["text"]))
blob = "\n".join(text)
anchors = ["-4.976573", "-5.634427",      # chain 4->8, U=4: E_low (Galerkin), E_high (ED)
           "-7.209813",                   # chain 16 sites (2,2), U=4: E_high
           "-10.4197",                    # 4x4 (2,2), U=4: E_high
           "(U - sqrt(U^2+16t^2))/2",     # dimer validation ran
           "sums to H"]                   # term decomposition validated
missing = [a for a in anchors if a not in blob]
ok = not errs and figs >= min_figs and not missing
print(f"errors={errs or 'none'}  figures={figs} (min {min_figs})  missing_anchors={missing or 'none'}")
print("CHECK PASSED" if ok else "CHECK FAILED")
sys.exit(0 if ok else 1)
