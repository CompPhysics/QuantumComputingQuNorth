#!/usr/bin/env python3
"""The 8 -> 16 chain at the circuit level (notebook Section 12, figs 11-12).

Chain of 16 sites, (2, 2), U = 4 (14 400 states; one-species dimension 120, so
every product-formula factor is a Kronecker product of 120 x 120 exponentials,
`hb.hubbard_groups`, and the embedded coarse term is applied exactly through the
784-dimensional coarse space, `hb.ProjGroup`).  The question left open by
Section 11: does the coarse eigenstate beat the free-fermion determinant once
the lattice is larger?  Steps: (a) inputs and spectra; (b) exact-slice ramp and
refinement with a shift scan; (c) product-formula steps and gate counts;
(d) UCCSD (1 022 parameters), its pruned version and ADAPT-VQE with the same
pool; (e) the Trotterized filter at these norms and the matched comparison.
Usage: PYTHONPATH=src python3 scripts/hubbard16_checks.py [outdir]
"""
import sys as _sys, os, time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix, identity
from scipy.sparse.linalg import eigsh
import pairinglib as pl
import pairinglib.hubbard as hb

OUT = _sys.argv[1] if len(_sys.argv) > 1 and not _sys.argv[1].startswith("-") else None   # None inside a notebook kernel
if OUT: matplotlib.use("Agg")
t0 = time.time()
U = 4.0
(Nc, bc), (Nf, bf), blocks, intra = hb.refine_chain(8)
Hf, (Sup, Sdn) = hb.H_hubbard(Nf, bf, 2, 2, U)
w16, V16 = eigsh(Hf, k=40, which="SA", tol=1e-10); o = np.argsort(w16); w16, V16 = w16[o], V16[:, o]; GS = V16[:, 0]; E0 = w16[0]
P16 = hb.prolong(Nc, Nf, blocks, 2, 2, sparse=True); ps = hb.prolonged_state(Hf, P16)
Hl = ps["Hl"]; El = ps["El"]; Phi0 = ps["Phi0"]
pool = hb._uccsd_pool_kron(Nf, bf, 2, 2); FF = hb._to_site(pool, pool["ref"]).real
Deff_p, c2p = hb.effective_gap(Phi0, w16, V16); Deff_f, c2f_ = hb.effective_gap(FF, w16, V16)
print(f"(a) chain 8->16, (2,2), U={U:.0f}: E_low = {El:.6f}, E_high = {E0:.6f}, sector gap {w16[1]-w16[0]:.4f}")
print(f"    prolonged coarse ground state: p0 = {c2p[0]:.4f}, D_eff = {Deff_p:.4f}, weight on the lowest 40 states {c2p.sum():.4f}")
print(f"    free-fermion determinant:      p0 = {c2f_[0]:.4f}, D_eff = {Deff_f:.4f}, weight on the lowest 40 states {c2f_.sum():.4f}   [{time.time()-t0:.0f} s]")

# ---- (b) exact slices (Krylov, 48 slices): ramp from the determinant, Path A with a shift scan
groups = hb.hubbard_groups(Nf, bf, 2, 2, U)
Kfree = hb.H_hubbard(Nf, bf, 2, 2, 0.0)[0]                      # the hopping part alone
Ts = (4, 8, 12, 16, 24, 32)
ramp_ex = {}
for T in Ts:
    psi = hb.refine_sparse(FF, Kfree, Hf, float(T), M=48); ramp_ex[T] = abs(np.vdot(GS, psi)) ** 2
print("(b) exact slices: U-ramp from the determinant: " + " ".join(f"p({T})={ramp_ex[T]:.5f}" for T in Ts) + f"   [{time.time()-t0:.0f} s]")
pathA_ex = {}
for mb in (2.0, 4.0, 8.0):
    mu = El + mb
    Hle = P16 @ csr_matrix(Hl - mu * np.eye(len(Hl))) @ P16.T; Hhs = Hf - mu * identity(Hf.shape[0])
    for T in (8, 16):
        psi = hb.refine_sparse(Phi0, Hle, Hhs, float(T), M=48); pathA_ex[(mb, T)] = abs(np.vdot(GS, psi)) ** 2
    print(f"    Path A, mu_buf={mb:.0f}: p(8) = {pathA_ex[(mb, 8)]:.5f}, p(16) = {pathA_ex[(mb, 16)]:.5f}   [{time.time()-t0:.0f} s]")

# ---- (c) product formulas and gate counts
c2 = lambda s: np.cos(np.pi * s / 2) ** 2; s2 = lambda s: np.sin(np.pi * s / 2) ** 2
st2 = hb.chain_step_cnots(Nf, bf, 2)["cnots"]; ct2 = hb.chain_step_cnots(Nf, bf, 2, controlled=True)["cnots"]
ca16 = hb.conditioned_term_cnots(8, bc, n_controls=16)["toffoli"]
CN_GIVENS16 = 2 * (Nf * (Nf - 1) // 2) * 2; CN_PROL16 = 2 * 8 * 2
print(f"(c) 16-site step: {st2} CNOTs (second order), controlled {ct2}; coarse term (8-site step conditioned on 16 qubits, Toffoli ladder) {ca16}; Givens network {CN_GIVENS16}; prolongation 8->16 {CN_PROL16}")
def n_needed(Phi, grp, coefs, T, pex, ns_list=(32, 64, 128, 256, 512), tol=1e-3):
    out = []
    for ns in ns_list:
        p = abs(np.vdot(GS, hb.trotter_path_groups(Phi, grp, coefs, T, ns, order=2))) ** 2; out.append((ns, p))
        if abs(p - pex) <= tol: return ns, p, out
    return ns, p, out
ramp_groups = groups; ramp_coefs = [lambda s: 1.0, lambda s: 1.0, s2]
trot = {}
for T in (8, 16):
    ns, p, out = n_needed(FF, ramp_groups, ramp_coefs, float(T), ramp_ex[T])
    trot[("ramp", T)] = dict(ns=ns, p=p, psi=hb.trotter_path_groups(FF, ramp_groups, ramp_coefs, float(T), ns, order=2), cost=CN_GIVENS16 + ns * st2)
    print(f"    ramp T={T}: exact {ramp_ex[T]:.5f}; " + " ".join(f"{n}:{pp:.5f}" for n, pp in out) + f" -> n_s={ns}, {trot[('ramp', T)]['cost']} CNOTs   [{time.time()-t0:.0f} s]")
tl, Ul, cst = hb.galerkin_params("chain", 1.0, U, 8)
Kc_e, Kc_o, Dc = [t_.toarray() for t_ in hb.hubbard_terms(8, bc, 2, 2, Ul, t=tl)]
MB = 4.0 if pathA_ex[(4.0, 8)] >= pathA_ex[(2.0, 8)] else 2.0          # the better shift at T = 8
mu = El + MB
gA = [hb.ProjGroup(P16, Dc + (cst - mu) * np.eye(784)), hb.ProjGroup(P16, Kc_e), hb.ProjGroup(P16, Kc_o), groups[0], groups[1], hb.DiagGroup(groups[2].d - mu)]
cA = [c2] * 3 + [s2] * 3
for T in (8, 16):
    ns, p, out = n_needed(Phi0, gA, cA, float(T), pathA_ex[(MB, T)])
    trot[("A", T)] = dict(ns=ns, p=p, psi=hb.trotter_path_groups(Phi0, gA, cA, float(T), ns, order=2), cost=CN_PROL16 + ns * (st2 + 2 * ca16))
    print(f"    Path A (split coarse groups, mu_buf={MB:.0f}) T={T}: exact {pathA_ex[(MB, T)]:.5f}; " + " ".join(f"{n}:{pp:.5f}" for n, pp in out) + f" -> n_s={ns}, {trot[('A', T)]['cost']} CNOTs (+ coarse preparation)   [{time.time()-t0:.0f} s]")

# ---- (d) variational inputs at 16 sites
uf = hb.uccsd_hubbard_kron(Nf, bf, 2, 2, U, maxiter=300, H=Hf)
inf_u = 1 - abs(np.vdot(GS, uf["psi"])) ** 2; th = np.abs(uf["theta"])
keep = [i for i in range(len(th)) if th[i] > 1e-3]
up = hb.uccsd_hubbard_kron(Nf, bf, 2, 2, U, maxiter=300, H=Hf, keep=keep, theta0=uf["theta"][keep])
cn_full = hb.uccsd_cnots(Nf, uf["singles"], uf["doubles"])["total"]; cn_pr = hb.uccsd_cnots(Nf, up["singles"], up["doubles"])["total"]
print(f"(d) UCCSD (16 sites): {uf['n_params']} parameters, {cn_full} CNOTs, E = {uf['E']:.6f} (E0 {E0:.6f}), 1-F = {inf_u:.3e}; "
      f"pruned (|theta| > 1e-3): {up['n_params']} parameters, {cn_pr} CNOTs, 1-F = {1-abs(np.vdot(GS, up['psi']))**2:.3e}   [{time.time()-t0:.0f} s]")
ad = hb.adapt_hubbard(Nf, bf, 2, 2, U, grad_tol=1e-6, max_ops=200, H=Hf, target=GS)
hist = ad["history"]
def cn_of(ops):
    return hb.uccsd_cnots(Nf, [l for l in ops if len(l) == 3], [l for l in ops if len(l) == 5])["total"]
ad_reach = {}
for tgt in (1e-2, 3e-3, 1e-3, 3e-4):
    hit = next((d for d in hist if d["infidelity"] <= tgt), None)
    if hit: ad_reach[tgt] = (hit["n_ops"], cn_of(ad["ops"][:hit["n_ops"]]), hit["psi"])
print(f"    ADAPT-VQE (pool {ad['pool_size']}, repeats allowed): {len(ad['ops'])} operators, 1-F = {hist[-1]['infidelity']:.2e}; operators (CNOTs incl. Givens) to reach "
      + ", ".join(f"{tgt:.0e}: {ad_reach[tgt][0]} ({ad_reach[tgt][1]})" if tgt in ad_reach else f"{tgt:.0e}: -" for tgt in (1e-2, 3e-3, 1e-3, 3e-4)) + f"   [{time.time()-t0:.0f} s]")
# the same at 8 sites for the comparison of sizes
H8, _ = hb.H_hubbard(8, bc, 2, 2, U); w8, V8 = np.linalg.eigh(H8.toarray()); GS8 = V8[:, 0]
ad8 = hb.adapt_hubbard(8, bc, 2, 2, U, grad_tol=1e-6, max_ops=200, H=H8, target=GS8)
h8 = ad8["history"]
print("    ADAPT-VQE at 8 sites (pool 198): " + ", ".join(f"{tgt:.0e}: {next((d['n_ops'] for d in h8 if d['infidelity'] <= tgt), '-')}" for tgt in (1e-2, 3e-3, 1e-3, 3e-4)) + f" operators; final 1-F {h8[-1]['infidelity']:.2e} after {len(ad8['ops'])}")

# ---- (e) the Trotterized filter at 16 sites and the matched comparison
SIG = 5.0 / Deff_p; DT = 0.1
psi_eff, E_eff = hb.effective_step_state(groups, DT, GS, E0, SIG, iters=80)
floor = 1 - abs(np.vdot(GS, psi_eff)) ** 2
print(f"(e) filter: sigma = 5/D_eff = {SIG:.1f}, dt = {DT}: E_eff - E0 = {E_eff-E0:+.2e}, floor {floor:.1e} (iterated quantized-time filter); "
      f"mean steps per cycle {SIG*np.sqrt(2/np.pi)/DT:.0f}, {2*ct2*SIG*np.sqrt(2/np.pi)/DT:.2e} CNOTs per cycle   [{time.time()-t0:.0f} s]")
NSAMP, M_MAX, SEED = 20, 6, 3
rng = np.random.default_rng(SEED); tsets = [np.abs(rng.normal(0.0, SIG, M_MAX)) for _ in range(NSAMP)]
nsteps = np.array([[max(1, int(round(t / DT))) for t in ts] for ts in tsets]); cyc = nsteps * ct2 * 2
# coarse preparation for the hierarchy: the cheapest 8-site route of Section 11 reaching 3e-4 (the U-ramp at T = 16, 18 544 CNOTs)
CN_COARSE8 = 18544
inputs = {"free-fermion determinant": (FF, CN_GIVENS16),
          "prolonged 8-site ground state (coarse prep. 18 544)": (Phi0, CN_COARSE8 + CN_PROL16),
          f"UCCSD pruned ({up['n_params']} par.)": (up["psi"], cn_pr),
          f"U-ramp T=8 (n_s={trot[('ramp', 8)]['ns']})": (trot[("ramp", 8)]["psi"], trot[("ramp", 8)]["cost"]),
          f"U-ramp T=16 (n_s={trot[('ramp', 16)]['ns']})": (trot[("ramp", 16)]["psi"], trot[("ramp", 16)]["cost"]),
          f"Path A (split) T=8 (n_s={trot[('A', 8)]['ns']})": (trot[("A", 8)]["psi"], CN_COARSE8 + trot[("A", 8)]["cost"]),
          f"Path A (split) T=16 (n_s={trot[('A', 16)]['ns']})": (trot[("A", 16)]["psi"], CN_COARSE8 + trot[("A", 16)]["cost"])}
if 1e-3 in ad_reach:
    inputs[f"ADAPT ({ad_reach[1e-3][0]} op.)"] = (ad_reach[1e-3][2], ad_reach[1e-3][1])
rows = {}; rngb = np.random.default_rng(5)
for name, (v0, cprep) in inputs.items():
    v0 = np.asarray(v0, complex); p = abs(np.vdot(GS, v0)) ** 2
    F = np.zeros((NSAMP, M_MAX)); Acc = np.zeros((NSAMP, M_MAX))
    for s_, ts in enumerate(tsets): F[s_], Acc[s_] = hb.rodeo_track_groups(v0, groups, ts, E_eff, GS, DT)
    Fw = (F * Acc).sum(0) / Acc.sum(0); Am = Acc.mean(0)
    boots = np.array([1 - ((F * Acc)[idx].sum(0) / Acc[idx].sum(0)) for idx in (rngb.integers(0, NSAMP, NSAMP) for _ in range(400))]); rse = boots.std(0) / np.maximum(1 - Fw, 1e-30)
    reach = np.concatenate([np.ones((NSAMP, 1)), Acc[:, :-1]], axis=1)
    rows[name] = dict(p=p, cprep=cprep, inf=np.concatenate([[1 - p], 1 - Fw]), rse=np.concatenate([[0], rse]), acc=np.concatenate([[1], Am]),
                      cacc=np.concatenate([[float(cprep)], (cprep + np.cumsum(reach * cyc, axis=1).mean(0)) / Am]))
    print(f"    {name:52s} p={p:.5f} prep={cprep:7d} | " + " ".join(f"{rows[name]['inf'][m]:.1e}({100*rows[name]['rse'][m]:2.0f}%)" for m in range(M_MAX + 1)) + f" | acc(6) {Am[-1]:.3f} | cost " + " ".join(f"{c:.1e}" for c in rows[name]["cacc"]) + f"   [{time.time()-t0:.0f} s]")
print("    minimum CNOTs per accepted preparation at a target infidelity:")
for tgt in (1e-2, 3e-3, 1e-3, 3e-4):
    best = []
    for name, d in rows.items():
        ok = np.where(d["inf"] <= tgt)[0]
        if ok.size: m = ok[np.argmin(d["cacc"][ok])]; best.append((d["cacc"][m], name, m, d["inf"][m]))
    best.sort(); print(f"      {tgt:.0e}: " + "; ".join(f"{nm} M={m} {c:.2e} ({i:.1e})" for c, nm, m, i in best))

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].semilogy([d["n_ops"] for d in h8], [max(d["infidelity"], 1e-12) for d in h8], "s--", ms=3, label="ADAPT, 8 sites (pool 198)")
ax[0].semilogy([d["n_ops"] for d in hist], [max(d["infidelity"], 1e-12) for d in hist], "o-", ms=3, label=f"ADAPT, 16 sites (pool {ad['pool_size']})")
ax[0].axhline(inf_u, color="C1", ls=":", label=f"single-layer UCCSD, 16 sites ({uf['n_params']} par.)")
ax[0].set_xlabel("ADAPT iteration (operators)"); ax[0].set_ylabel("$1-F$"); ax[0].set_title("ADAPT-VQE on the Hubbard chain, $U=4t$, (2,2)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3, which="both")
mks = ["o-", "s-", "d-", "^-", "v-", "P-", "X-", "*-"]
for (name, d), mk in zip(rows.items(), mks):
    ax[1].loglog(d["cacc"], d["inf"], mk, ms=4, label=f"{name} ($p={d['p']:.4f}$)")
for tgt in (3e-3, 1e-3, 3e-4): ax[1].axhline(tgt, color="gray", ls=":", lw=0.8)
ax[1].set_xlabel("CNOTs per accepted preparation (input + filter)"); ax[1].set_ylabel("$1-\\mathcal{F}_M$"); ax[1].set_title(f"chain $8\\to16$, $U=4t$: inputs to the same filter ($\\delta t={DT}$, $\\sigma={SIG:.0f}$), $M=0\\ldots6$"); ax[1].legend(fontsize=6); ax[1].grid(alpha=0.3, which="both")
plt.tight_layout()
if OUT: plt.savefig(os.path.join(OUT, "fig11.png"), dpi=150); print("wrote", OUT)
else: plt.show()
print(f"[{time.time()-t0:.0f} s]")
