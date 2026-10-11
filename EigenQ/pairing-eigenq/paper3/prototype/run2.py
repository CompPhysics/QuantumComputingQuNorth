import numpy as np, json, time
from scipy.sparse import csr_matrix, kron, identity
from scipy.sparse.linalg import expm_multiply, eigsh
from hubbardrefine import *
t0=time.time()
R = {}
# ---------- (1) effective gaps + rodeo with sigma from Delta_eff, chain 4->8 (2,2)
R['eff'] = {}
for U in [0.0,2.0,4.0,8.0]:
    (Nc,bc),(Nf,bf),blocks,intra = refine_chain(4)
    Hf,_ = H_hubbard(Nf,bf,2,2,U); Hf=Hf.toarray(); wf,Vf=np.linalg.eigh(Hf)
    P = prolong(Nc,Nf,blocks,2,2); Hl=P.T@Hf@P; wl,Vl=np.linalg.eigh(Hl); Phi0=P@Vl[:,0]
    c2 = np.abs(Vf.T@Phi0)**2
    idx = np.where(c2>1e-10)[0]; Deff = wf[idx[1]]-wf[0]
    low = [(float(wf[k]-wf[0]), float(c2[k])) for k in range(1,8)]
    # refined states T=16 (A mu2, B Lam16) with sigma=5/Deff
    sig = 5.0/Deff
    Sup,Sdn = spin_basis(Nf,2),spin_basis(Nf,2)
    Kin = (kron(hop_spin(Nf,intra,Sup),identity(len(Sdn)),format='csr')+kron(identity(len(Sup)),hop_spin(Nf,intra,Sdn),format='csr')).toarray()
    out = dict(Deff=float(Deff), gap=float(wf[1]-wf[0]), low=low, sigma=float(sig), ov0=float(c2[0]),
               rodeo_prol=rodeo_infidelity(Phi0.astype(complex),wf,Vf,list(range(1,13)),sig).tolist())
    mu = wl[0]+2.0
    for key,(Hle,Hhs) in {"A_mu2":(P@(Hl-mu*np.eye(len(Hl)))@P.T, Hf-mu*np.eye(len(Hf))), "B_Lam16":(Hf+16.0*Kin,Hf)}.items():
        eigs = slice_eigs(Hle,Hhs,96)
        for T in (16,32):
            psi = evolve(Phi0,eigs,T)
            out[f"{key}_T{T}"] = dict(p=float(abs(np.vdot(Vf[:,0],psi))**2), infid=rodeo_infidelity(psi,wf,Vf,list(range(1,13)),sig).tolist())
    R['eff'][U]=out
    print(f"U={U}: gap={out['gap']:.4f} Deff={Deff:.4f} sigma={sig:.2f} low(dE,w)={[(round(a,3),f'{b:.1e}') for a,b in low]}", flush=True)
    for k,v in out.items():
        if isinstance(v,dict): print(f"   {k}: p={v['p']:.5f} 1-F(M=1,2,4,6)="+" ".join(f"{v['infid'][m-1]:.1e}" for m in (1,2,4,6)), flush=True)
    print(f"   prolonged 1-F(M=1,2,4,6,10)="+" ".join(f"{out['rodeo_prol'][m-1]:.1e}" for m in (1,2,4,6,10)), f"[{time.time()-t0:.0f}s]", flush=True)

# ---------- (2) filled dimer 2->4 (2,2): coarse dim 1
R['dimer'] = {}
for U in [0.0,2.0,4.0,8.0]:
    (Nc,bc),(Nf,bf),blocks,intra = refine_chain(2)
    Hf,_ = H_hubbard(Nf,bf,2,2,U); Hf=Hf.toarray(); wf,Vf=np.linalg.eigh(Hf)
    P = prolong(Nc,Nf,blocks,2,2); Hl=P.T@Hf@P; El=Hl[0,0]; Phi0=P[:,0]
    c2=np.abs(Vf.T@Phi0)**2; idx=np.where(c2>1e-10)[0]; Deff=wf[idx[1]]-wf[0]
    mu=El+2.0; Hle=P@(Hl-mu)@P.T; Hhs=Hf-mu*np.eye(len(Hf))
    s,g = path_gap(Hle,Hhs,21); eigs=slice_eigs(Hle,Hhs,96)
    Ts=[2,4,8,16,32]; ov=[float(abs(np.vdot(Vf[:,0],evolve(Phi0,eigs,T)))**2) for T in Ts]
    R['dimer'][U]=dict(E_low=float(El),E_high=float(wf[0]),ov0=float(c2[0]),gap=float(wf[1]-wf[0]),Deff=float(Deff),min_gap=float(g.min()),Ts=Ts,ov=ov)
    print(f"dimer 2->4 U={U}: E_low={El:.4f} E_high={wf[0]:.4f} ov0={c2[0]:.4f} gap={wf[1]-wf[0]:.3f} Deff={Deff:.3f} mingap={g.min():.3f} ov@T={list(zip(Ts,[round(o,4) for o in ov]))}", flush=True)

# ---------- (3) 2D 2x2 -> 4x4 (1,1) and (2,2 prolonged overlap only)
R['sq'] = {}
for U in [0.0,4.0,8.0]:
    (Nc,bc),(Nf,bf),blocks,intra = refine_square(2)
    Hf,_ = H_hubbard(Nf,bf,1,1,U); Hf=Hf.toarray(); wf,Vf=np.linalg.eigh(Hf)
    P = prolong(Nc,Nf,blocks,1,1); Hl=P.T@Hf@P; wl,Vl=np.linalg.eigh(Hl); Phi0=P@Vl[:,0]
    c2=np.abs(Vf.T@Phi0)**2; idx=np.where(c2>1e-10)[0]; Deff=wf[idx[1]]-wf[0]
    mu=wl[0]+2.0; Hle=P@(Hl-mu*np.eye(len(Hl)))@P.T; Hhs=Hf-mu*np.eye(len(Hf))
    s,g=path_gap(Hle,Hhs,21); eigs=slice_eigs(Hle,Hhs,96)
    Ts=[4,8,16,32,64]; ov=[float(abs(np.vdot(Vf[:,0],evolve(Phi0,eigs,T)))**2) for T in Ts]
    R['sq'][U]=dict(E_low=float(wl[0]),E_high=float(wf[0]),ov0=float(c2[0]),gap=float(wf[1]-wf[0]),Deff=float(Deff),min_gap=float(g.min()),Ts=Ts,ov=ov,dim_low=Hl.shape[0],dim_high=Hf.shape[0])
    print(f"2x2->4x4 (1,1) U={U}: dims {Hl.shape[0]}->{Hf.shape[0]} E_low={wl[0]:.4f} E_high={wf[0]:.4f} ov0={c2[0]:.4f} gap={wf[1]-wf[0]:.3f} Deff={Deff:.3f} mingap={g.min():.3f} ov={list(zip(Ts,[round(o,4) for o in ov]))} [{time.time()-t0:.0f}s]", flush=True)
# (2,2): ED + prolonged overlap only (dim 14400, sparse)
for U in [4.0]:
    (Nc,bc),(Nf,bf),blocks,intra = refine_square(2)
    Hf,_ = H_hubbard(Nf,bf,2,2,U); w,V = eigsh(Hf,k=3,which='SA',tol=1e-10); o=np.argsort(w); w=w[o]; V=V[:,o]
    P = prolong(Nc,Nf,blocks,2,2); Pc=csr_matrix(P); Hl=(Pc.T@(Hf@Pc)).toarray(); wl,Vl=np.linalg.eigh(Hl); Phi0=P@Vl[:,0]
    ov0=abs(V[:,0]@Phi0)**2
    R['sq'][f"22_U{U}"]=dict(E_low=float(wl[0]),E_high=float(w[0]),ov0=float(ov0),gap=float(w[1]-w[0]),dim_low=Hl.shape[0],dim_high=Hf.shape[0])
    print(f"2x2->4x4 (2,2) U={U}: dims {Hl.shape[0]}->{Hf.shape[0]} E_low={wl[0]:.4f} E_high={w[0]:.4f} ov0={ov0:.4f} gap={w[1]-w[0]:.3f} [{time.time()-t0:.0f}s]", flush=True)

# ---------- (4) chain 8->16 (2,2), U=4: ED + prolonged overlap + refinement at T=16 (sparse Krylov, 48 slices)
(Nc,bc),(Nf,bf),blocks,intra = refine_chain(8)
U=4.0
Hf,_ = H_hubbard(Nf,bf,2,2,U); w,V = eigsh(Hf,k=3,which='SA',tol=1e-10); o=np.argsort(w); w=w[o]; V=V[:,o]
P = prolong(Nc,Nf,blocks,2,2); Pc=csr_matrix(P); Hl=(Pc.T@(Hf@Pc)).toarray(); wl,Vl=np.linalg.eigh(Hl); Phi0=P@Vl[:,0]
ov0=abs(V[:,0]@Phi0)**2
print(f"chain 8->16 (2,2) U={U}: dims {Hl.shape[0]}->{Hf.shape[0]} E_low={wl[0]:.4f} E_high={w[0]:.4f} ov0={ov0:.4f} gap={w[1]-w[0]:.3f} [{time.time()-t0:.0f}s]", flush=True)
mu=wl[0]+2.0
Hle = Pc@csr_matrix(Hl-mu*np.eye(len(Hl)))@Pc.T; Hhs = Hf - mu*identity(Hf.shape[0])
res816={}
for T in (16,32):
    Msl=48; psi=Phi0.astype(complex); dt=T/Msl
    for m in range(Msl):
        s=(m+0.5)/Msl; H=np.cos(np.pi*s/2)**2*Hle+np.sin(np.pi*s/2)**2*Hhs
        psi=expm_multiply(-1j*dt*H.tocsc(),psi)
    ovT=abs(np.vdot(V[:,0],psi))**2; res816[T]=float(ovT)
    print(f"   T={T}: ov={ovT:.5f} [{time.time()-t0:.0f}s]", flush=True)
R['c816']=dict(E_low=float(wl[0]),E_high=float(w[0]),ov0=float(ov0),gap=float(w[1]-w[0]),ov=res816,dim_low=Hl.shape[0],dim_high=Hf.shape[0])
json.dump(R,open('run2.json','w'),indent=1,default=str)
print("DONE")
