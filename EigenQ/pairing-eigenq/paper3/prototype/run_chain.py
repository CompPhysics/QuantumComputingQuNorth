import numpy as np, json, time
from hubbardrefine import *

def study(L, Nup, Ndn, U, mu_bufs=(1.0,2.0), Lams=(8.0,16.0), Ts=None, M=96, sigma=None, label=""):
    (Nc,bc),(Nf,bf),blocks,intra = refine_chain(L)
    Hf,_ = H_hubbard(Nf, bf, Nup, Ndn, U); Hf = Hf.toarray()
    wf, Vf = np.linalg.eigh(Hf); Eh = wf[0]; Psih = Vf[:,0]
    P = prolong(Nc, Nf, blocks, Nup, Ndn)
    Hl = P.T @ Hf @ P; wl, Vl = np.linalg.eigh(Hl); El = wl[0]; psil = Vl[:,0]
    Phi0 = P @ psil
    ov0 = abs(Psih @ Phi0)**2
    res = dict(label=label, L=L, Nup=Nup, Ndn=Ndn, U=U, dim_low=Hl.shape[0], dim_high=Hf.shape[0],
               E_low=El, E_high=Eh, gap_low=wl[1]-wl[0], gap_high=wf[1]-wf[0],
               E_prolonged=float(Phi0 @ Hf @ Phi0), ov0=ov0, paths={})
    if Ts is None: Ts = np.array([0,1,2,3,4,6,8,10,12,16,20,24,32,40,48,64,80])
    # Path A
    for mb in mu_bufs:
        mu = El + mb
        Hle = P @ (Hl - mu*np.eye(len(Hl))) @ P.T; Hhs = Hf - mu*np.eye(len(Hf))
        s, g = path_gap(Hle, Hhs, 21)
        eigs = slice_eigs(Hle, Hhs, M)
        ov = [abs(np.vdot(Psih, evolve(Phi0, eigs, T)))**2 for T in Ts]
        res['paths'][f"A_mu{mb}"] = dict(min_gap=float(g.min()), s_min=float(s[g.argmin()]), gap_s0=float(g[0]), gap_s1=float(g[-1]), Ts=Ts.tolist(), ov=ov)
    # Path B (orbital/dimer lowering): extra intra-block hopping Lambda*cos^2
    Nf_, _ = chain(2*L)
    Sup, Sdn = spin_basis(Nf, Nup), spin_basis(Nf, Ndn)
    Kin = (kron(hop_spin(Nf, intra, Sup), identity(len(Sdn)), format='csr') + kron(identity(len(Sup)), hop_spin(Nf, intra, Sdn), format='csr')).toarray()
    for Lam in Lams:
        Hle = Hf + Lam*Kin; Hhs = Hf
        s, g = path_gap(Hle, Hhs, 21)
        w0 = np.linalg.eigh(Hle); ov_start = abs(np.vdot(w0[1][:,0], Phi0))**2
        eigs = slice_eigs(Hle, Hhs, M)
        ov = [abs(np.vdot(Psih, evolve(Phi0, eigs, T)))**2 for T in Ts]
        res['paths'][f"B_Lam{Lam}"] = dict(min_gap=float(g.min()), s_min=float(s[g.argmin()]), gap_s0=float(g[0]), gap_s1=float(g[-1]), ov_gs0_vs_prolonged=float(ov_start), Ts=Ts.tolist(), ov=ov)
    # rodeo: exact ensemble infidelity for prolonged input and for refined inputs
    if sigma is None: sigma = 5.0/(wf[1]-wf[0])
    Ms = list(range(1,13))
    res['rodeo'] = dict(sigma=float(sigma), Ms=Ms, prolonged=rodeo_infidelity(Phi0.astype(complex), wf, Vf, Ms, sigma).tolist())
    # refined inputs at a few T for the best path A (mu 1.0) and B (Lam 8)
    for key, (Hle, Hhs) in {"A_mu1.0": (P @ (Hl - (El+1.0)*np.eye(len(Hl))) @ P.T, Hf-(El+1.0)*np.eye(len(Hf))), "B_Lam8.0": (Hf+8.0*Kin, Hf)}.items():
        eigs = slice_eigs(Hle, Hhs, M)
        for T in (8, 16, 32):
            psi = evolve(Phi0, eigs, T)
            res['rodeo'][f"{key}_T{T}"] = dict(p=float(abs(np.vdot(Psih,psi))**2), infid=rodeo_infidelity(psi, wf, Vf, Ms, sigma).tolist())
    return res

out = []
t0=time.time()
for U in [0.0, 2.0, 4.0, 8.0]:
    r = study(4, 2, 2, U, label="chain4to8_N22")
    out.append(r)
    print(f"\n=== chain 4->8, (2,2), U={U} ===  dims {r['dim_low']}->{r['dim_high']}  [{time.time()-t0:.0f}s]")
    print(f" E_low={r['E_low']:.6f} E_high={r['E_high']:.6f} E_prol={r['E_prolonged']:.6f} gap_low={r['gap_low']:.4f} gap_high={r['gap_high']:.4f} ov0={r['ov0']:.4f}")
    for k,v in r['paths'].items():
        ovs = " ".join(f"{T}:{o:.4f}" for T,o in zip(v['Ts'],v['ov']) if T in (4,8,16,32,64))
        extra = f" ov(gs(s=0),prol)={v['ov_gs0_vs_prolonged']:.4f}" if 'ov_gs0_vs_prolonged' in v else ""
        print(f"  {k:10s} min_gap={v['min_gap']:.4f}@s={v['s_min']:.3f} gap0={v['gap_s0']:.3f} gap1={v['gap_s1']:.3f}{extra} | ov {ovs}")
    rd = r['rodeo']
    print(f"  rodeo sigma={rd['sigma']:.2f}: prolonged 1-F M=1,2,4,6,10: " + " ".join(f"{rd['prolonged'][m-1]:.1e}" for m in (1,2,4,6,10)))
    for k,v in rd.items():
        if isinstance(v, dict): print(f"    {k}: p={v['p']:.4f} 1-F M=1,2,4,6: " + " ".join(f"{v['infid'][m-1]:.1e}" for m in (1,2,4,6)))
json.dump(out, open("chain4to8.json","w"), indent=1)
