#!/usr/bin/env python3
"""Revision of 2026-10-14 (paper/article_correctness_assessment.md): insert the
tagged section EXT:ASSESS (code from scripts/assessment_checks.py) before the
Summary cell, and patch three original cells so that the notebook's own numbers
match the article's conventions:
  * the fine rodeo energy scan is centered on the COARSE-scan peak, not on the
    exact E0, and reports a shot-noise uncertainty of the located peak;
  * the rodeo-convergence cell resets the random generator (seed 7) for every
    input, plots the exact Gaussian-ensemble curve of Eq. (fid) and error bars
    from the standard error of a 40-schedule mean;
  * the shift-scan panel shows both affine fits (1/Delta and 1/Delta^2).
Built with nbformat; idempotent.  Run from the repo root:
  python3 scripts/extend_notebook_assess.py
  make notebook && make check    # full clean run; figure count 18 -> 20
"""
import pathlib, re
import nbformat as nbf

ROOT = pathlib.Path(__file__).resolve().parents[1]
NB = ROOT / "notebooks/ResolutionRefPairing.ipynb"
SCRIPTS = ROOT / "scripts"
MARK = "<!-- EXT:ASSESS -->"

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
    src = src.replace('import sys, os\n', '')
    return "# " + MARK + "\n" + PREAMBLE + src


md = MARK + r"""
---
## 19&nbsp; Revision checks: scaling, statistics, the pruned comparator and the missing baselines

An independent assessment of the article (October 2026) questioned several
interpretations rather than the physics; this section settles each point numerically.
(a) The four shift-scan points of Section 15 are fitted with an intercept to both
$a/\Delta_{\min}+b$ and $a/\Delta_{\min}^2+b$; the response norm
$\chi_u(s)=\sum_{n>0}|\langle n|\partial_uH|0\rangle|^2/(E_n-E_0)^2$ of the path, which the
inverse-gap argument of Bogner *et al.* requires to be uniformly bounded, is evaluated along
the path for every shift (at $s=0$ it is exactly $g^2/b^2$ with $b=\mu-E_{\rm low}$).  The
derivative couplings of the tracked state at $s=0.25$ and the convergence of the exact-slice
propagation in the number of slices are also recorded.  (b) The sampled rodeo infidelities of
Section 11 are compared with the exact Gaussian-ensemble expression and given their standard
errors.  (c) The seniority-zero reference admits an exact four-parameter preparation: the four
occupied-pair$\to$virtual-pair rotations of the UCCSD pool, 192 CNOTs; we optimize it and show
that it reproduces the full 26-operator UCCSD state, whose other 22 amplitudes vanish.
(d) The $k=3\to4$ chain is re-timed on a finer grid and Trotterized with the embedded coarse
term split into its own coarse groups (the counted circuit), and two baselines that the
comparison of Section 17 lacked are added: the prolonged $k=3$ UCCSD state fed directly to
the filter, and an ordinary adiabatic ramp from the diagonal (HF) Hamiltonian,
$H(s)=D+\sin^2(\pi s/2)\sum_{pq}h_{pq}$, which has no shift and no coarse term.  (e) The
long-run state of the digital filter is obtained from the leading eigenvector of the averaged
post-selection map $\mathcal K(\rho)=\mathbb E_t[\tilde K_t\rho\tilde K_t^\dagger]$ in the
seniority-zero block.  (f) Eight inputs are sent through the same Trotterized filter, including
$M=0$, with CNOTs per attempt and per accepted preparation, and the minimum cost reaching three
target infidelities is tabulated.
"""

cells = [("markdown", md), ("code", cell_from_script("assessment_checks.py"))]

nb = nbf.read(NB, as_version=4)
nb.cells = [c for c in nb.cells if MARK not in c.source]      # idempotent
isum = next(i for i, c in enumerate(nb.cells)
            if c.cell_type == "markdown" and "Summary: the complete EIGEN-Q pipeline" in c.source)
new = [nbf.v4.new_markdown_cell(s) if k == "markdown" else nbf.v4.new_code_cell(s)
       for k, s in cells]
nb.cells = nb.cells[:isum] + new + nb.cells[isum:]

# ---------------------------------------------------------------- patch the energy-scan cell
for c in nb.cells:
    if c.cell_type == "code" and "Es_fine   = np.linspace(E0-0.6, E0+0.6, 161)" in c.source:
        src = c.source
        src = src.replace(
            "Es_coarse = np.linspace(-0.5, 4.0, 60)\nEs_fine   = np.linspace(E0-0.6, E0+0.6, 161)\nPc_pro = scan(psi_prolonged, Es_coarse, 2.0, Hrod)\nPc_ref = scan(psi_refined,   Es_coarse, 2.0, Hrod)\nPf_ref = scan(psi_refined,   Es_fine,   8.0, Hrod)",
            "Es_coarse = np.linspace(-0.5, 4.0, 60)\nPc_pro = scan(psi_prolonged, Es_coarse, 2.0, Hrod)\nPc_ref = scan(psi_refined,   Es_coarse, 2.0, Hrod)\n"
            "# the fine window is chosen from the COARSE scan alone (its highest point), not from the exact E0\nE_coarse_pk = Es_coarse[np.argmax(Pc_ref)]\n"
            "Es_fine   = np.linspace(E_coarse_pk-0.6, E_coarse_pk+0.6, 161)\nPf_ref = scan(psi_refined,   Es_fine,   8.0, Hrod)")
        src = src.replace(
            'print(f"peak of fine scan at E = {Epk:.4f}  (exact E0 = {E0:.6f})")',
            'print(f"coarse-scan maximum at E = {E_coarse_pk:.4f}; fine window [{Es_fine[0]:.3f}, {Es_fine[-1]:.3f}] (grid spacing {Es_fine[1]-Es_fine[0]:.4f})")\n'
            'print(f"peak of fine scan at E = {Epk:.4f}  (exact E0 = {E0:.6f}; error {Epk-E0:+.4f})")\n'
            '# shot noise: every acceptance is estimated from n_shots runs; locate the peak by a quadratic fit\n'
            '# over the points above 0.8 x max, with binomial noise added, 200 times\n'
            'def peak_quadratic(y):\n'
            '    m = y >= 0.8*y.max(); c = np.polyfit(Es_fine[m], y[m], 2); return -c[1]/(2*c[0])\n'
            'rng_sh = np.random.default_rng(11)\n'
            'for n_shots in (10**3, 10**4):\n'
            '    pk = [peak_quadratic(rng_sh.binomial(n_shots, np.clip(Pf_ref, 0, 1))/n_shots) for _ in range(200)]\n'
            '    print(f"  with {n_shots:5d} shots per energy: peak {np.mean(pk):.4f} +/- {np.std(pk):.4f} (quadratic fit over the top 20% of the scan; noise-free fit {peak_quadratic(Pf_ref):.4f})")')
        c.source = src
        print("patched: energy-scan cell")

# ---------------------------------------------------------------- patch the rodeo-convergence cell
for c in nb.cells:
    if c.cell_type == "code" and "M_max, n_sets, sigma = 10, 40, 3.0" in c.source and "rng = np.random.default_rng(7)" in c.source:
        src = c.source
        src = src.replace(
            "M_max, n_sets, sigma = 10, 40, 3.0\nrng = np.random.default_rng(7)\n",
            "M_max, n_sets, sigma = 10, 40, 3.0\n# convention (article, Sec. IV.F): the generator is RESET to seed 7 for every input,\n# so each input sees the same 40 Gaussian time sets\n")
        src = src.replace(
            "results = {}\nfor name, v0, p0, c in inputs:\n    F = np.zeros(M_max); A = np.zeros(M_max)\n    for s in range(n_sets):\n        f, a = rodeo_track(v0, rng.normal(0, sigma, M_max), E0, Hrod, GS)\n        F += f; A += a\n    results[name] = (F/n_sets, A/n_sets, p0, c)",
            "results = {}\nqk = 0.5*(1 + np.exp(-(wR - E0)**2*sigma**2/2))          # Gaussian-averaged single-cycle factors\n"
            "for name, v0, p0, c in inputs:\n    rng = np.random.default_rng(7)\n    Fs = np.zeros((n_sets, M_max)); As = np.zeros((n_sets, M_max))\n"
            "    for s in range(n_sets):\n        Fs[s], As[s] = rodeo_track(v0, rng.normal(0, sigma, M_max), E0, Hrod, GS)\n"
            "    c2 = np.abs(VR.T @ v0)**2; Ms_ = np.arange(1, M_max+1)\n"
            "    exact = np.array([np.sum(c2[1:]*qk[1:]**M)/np.sum(c2*qk**M) for M in Ms_])   # Eq. (fid), exact ensemble\n"
            "    # standard error of a 40-schedule mean of the residual weight, from 20000 schedules\n"
            "    rng2 = np.random.default_rng(123); nS = 20000; res = np.zeros((nS, M_max)); dE = (wR - E0)[None, :]\n"
            "    for s in range(nS):\n        ts = rng2.normal(0, sigma, M_max)\n        res[s] = (np.cumprod(np.cos(dE*ts[:, None]/2)**2, axis=0)[:, 1:]*c2[None, 1:]).sum(1)\n"
            "    rse = res.std(0)/res.mean(0)/np.sqrt(n_sets)\n"
            "    results[name] = (Fs.mean(0), As.mean(0), p0, c, 1-(Fs*As).sum(0)/As.sum(0), exact, rse)")
        src = src.replace(
            "for name, (F, A, p0, c) in results.items():\n    axL.semilogy(Ms, 1-F, 'o-', color=c, label=name)\n    axL.semilogy(Ms, (1-p0)/p0 * 2.0**(-Ms), ':', color=c, lw=1.2)",
            "for name, (F, A, p0, c, infw, exact, rse) in results.items():\n    axL.errorbar(Ms, infw, yerr=infw*rse, fmt='o-', color=c, capsize=3, label=name + ' (sampled, 40 sets)')\n    axL.semilogy(Ms, exact, ':', color=c, lw=1.5, label='exact ensemble, Eq. (fid)')")
        src = src.replace("axL.set_title('Exponential convergence  (dotted: $\\\\frac{1-p}{p}2^{-M}$)')",
                          "axL.set_title('Exponential convergence (acceptance-weighted; bars: std. error of a 40-set mean)')")
        src = src.replace("for name, (F, A, p0, c) in results.items():\n    axR.plot(Ms, A, 'o-', color=c, label=name)",
                          "for name, (F, A, p0, c, infw, exact, rse) in results.items():\n    axR.plot(Ms, A, 'o-', color=c, label=name)")
        src = src.replace(
            'print("Final fidelities (M = %d, averaged over %d time sets):" % (M_max, n_sets))\nfor name, (F, A, p0, c) in results.items():\n    print(f"  {name:32s}: F={F[-1]:.6f}  1-F={1-F[-1]:.1e}  acceptance={A[-1]:.4f}")',
            'print("Acceptance-weighted ensemble infidelity: sampled (40 time sets, seed 7 reset per input) vs the exact Gaussian ensemble, Eq. (fid)")\n'
            'for name, (F, A, p0, c, infw, exact, rse) in results.items():\n'
            '    print(f"  {name}:")\n'
            '    for M in (1, 2, 4, 6, 10):\n'
            '        print(f"     M={M:2d}: sampled {infw[M-1]:.2e} (rel. std. error {100*rse[M-1]:.0f}%)   exact {exact[M-1]:.2e}   acceptance {A[M-1]:.4f}")')
        c.source = src
        print("patched: rodeo-convergence cell")

# ---------------------------------------------------------------- patch the shift-scan panel (both fits)
for c in nb.cells:
    if c.cell_type == "code" and "mubufs = [0.15, 0.3, 0.6, 1.0, 1.5]" in c.source:
        src = c.source
        src = src.replace(
            'c = np.polyfit(x[sel], y[sel], 1)\nxs = np.linspace(0, x[sel].max()*1.08, 50)\nax[1].plot(x, y, "o", ms=8)\nax[1].plot(xs, np.polyval(c, xs), "--",\n           label=rf"linear fit: $T^*\\approx{c[0]:.1f}/\\Delta E_{{\\min}}{c[1]:+.1f}$")',
            'c = np.polyfit(x[sel], y[sel], 1)\nA2 = np.vstack([x[sel]**2, np.ones(sel.sum())]).T; c2_ = np.linalg.lstsq(A2, y[sel], rcond=None)[0]\n'
            'ssr1 = np.sum((np.polyval(c, x[sel]) - y[sel])**2); ssr2 = np.sum((A2 @ c2_ - y[sel])**2)\n'
            'xs = np.linspace(0.5, x[sel].max()*1.08, 50)\nax[1].plot(x, y, "o", ms=8)\n'
            'ax[1].plot(xs, np.polyval(c, xs), "--",\n           label=rf"affine in $1/\\Delta E_{{\\min}}$: $T^*\\approx{c[0]:.1f}/\\Delta E_{{\\min}}{c[1]:+.1f}$ (SSR {ssr1:.1f})")\n'
            'ax[1].plot(xs, c2_[0]*xs**2 + c2_[1], "-",\n           label=rf"affine in $1/\\Delta E_{{\\min}}^2$: $T^*\\approx{c2_[0]:.1f}/\\Delta E_{{\\min}}^2{c2_[1]:+.1f}$ (SSR {ssr2:.2f})")\n'
            'print(f"fits with intercept: a/D + b: a={c[0]:.2f} b={c[1]:+.2f} SSR={ssr1:.2f};  a/D^2 + b: a={c2_[0]:.2f} b={c2_[1]:+.2f} SSR={ssr2:.2f}  (four points on a coarse grid: no exponent is identified)")')
        src = src.replace('ax[1].set_title(r"$T^{*}$ scales linearly with the inverse path gap")',
                          'ax[1].set_title(r"$T^{*}$ versus the inverse path gap: both affine fits")')
        c.source = src
        print("patched: shift-scan cell")

nbf.write(nb, NB)
print(f"inserted {len(new)} cells before cell {isum}; notebook now has {len(nb.cells)} cells")
