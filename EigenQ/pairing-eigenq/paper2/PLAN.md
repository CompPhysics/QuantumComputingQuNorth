# Paper 2 — refinement + rodeo for the two-fragment pair-transfer model

Status 2026-10-10: exploration only, no calculations. Companion inputs:
`paper2/twontransfer.pdf` (two-fragment model, "the draft") and
`paper/eigenq_pairing.tex` (refinement + rodeo on the pairing model, "paper 1").
Skeleton: `paper2/refinement_rodeo_transfer.tex` (RevTeX, compiles; red `[...]`
marks everything still to be computed or decided).

## 1. What the draft gives us

- Two seniority-zero fragments, each a Richardson-type pairing Hamiltonian with
  realistic ε_p and degeneracies Ω_p ≤ 6; G_X = 20.5/A MeV (0.177 MeV for Sn,
  0.342 MeV for Ni). Entrance: |GS_T(2)⟩ ⊗ |GS_P(5)⟩. Collision: Gaussian pulse
  of a uniform pair-hopping term between fragments; pair number conserved.
- Pair-slot encoding → **strictly two-local** qubit Hamiltonian
  H = E_c + Σ ζ_i Z_i + Σ J_ij (XX+YY). 116Sn+60Ni: 30 qubits (11 T + 19 P),
  30 Z terms, 435 XY terms (55 within T, 171 within P, 209 cross, the cross
  ones only during the pulse). 96Zr+40Ca: 26 qubits, dim 145.
- Dimensions (checked): d_T(2)=9, d_P(5)=154, full M=7 space 5554; exit sector
  (3,4): 15×93; fixed-M slot sector for a gate-level simulator: C(30,7)=2.0e6.
- Observables: P_κ(D), P_{+1}/P_0 vs data; R_2; S_{T|P} = S_num + S_cfg with
  S_cfg tiny (0.3% at D=13.47 fm, ≤2% at the smallest D).
- Reference numbers at D=13.4669 fm: P_0=0.97209, P_{+1}=0.02206,
  P_{-1}=0.00566, P_{+2}=1.96e-4, S=0.20622 bits (S_cfg=4.31e-4).
- τ = 2.251e-22 s = 0.342 ħ/MeV; the ±5τ pulse is 3.42 ħ/MeV long.

## 2. Where paper-1 machinery plugs in (the ideas, ranked)

1. **Entrance-state preparation by orbital-by-orbital refinement from a filled
   configuration.** "Resolution" = active-orbital set at fixed M_X. Prolongation
   = append empty slots (gate-free). Coarse seeds are single basis states:
   Ni: 2 pairs in the 2 slots of 2p3/2; Sn: 5 pairs in the 4+1 slots of
   1g7/2+3s1/2. Chains (slots / sym. dim): T: 2/1 → 6/5 → 11/9;
   P: 5/1 → 7/5 → 13/26 → 19/154. Zero variational parameters anywhere.
   This is exactly the merit paper 1 found for refinement (optimizer- and
   measurement-free), in a setting where that merit is decisive.
2. **Two refinement paths.** Path A = paper-1 projector path; its coarse term
   P(H_low−μ)P† = Π_0 H_high Π_0 is a multi-controlled operator over the *new*
   slots (up to 6 for 1h11/2) — the cost driver identified in paper 1, worse
   here. Path B ("orbital lowering"): keep the full two-local H and sweep the
   new orbitals' ε_p down from +Λ: H(s) = H_high + Λ cos²(πs/2) Σ_new n_i.
   Same end points up to O((GΩ/Λ)²), every step costs a collision step, and it
   reads as "releasing the frozen-orbital approximation adiabatically".
   A vs B (T* × cost per step) is a genuinely new result.
3. **Rodeo filtering of each fragment** (11 and 19 qubits independently, one
   shared ancilla), E_obj = E_X^(0)(M_X); a rodeo *scan* doubles as an on-device
   check that refinement landed on the right state.
4. **The collision is the same primitive.** Second-order Trotter of the
   two-local H with midpoint rule in G_PT(t). Key trick: interaction picture
   w.r.t. the commuting Z part (ε ~ 10 MeV ≫ G ≤ 0.35 MeV) keeps the
   Hamiltonian two-local (XX+YY rotates into XY−YX) and makes the Trotter error
   O(G²) instead of O(Gε) — probably very few steps for a 3.4 ħ/MeV pulse.
   Also: Trotter leaks out of the orbital-symmetric subspace (unphysical
   configurations with different pairing energy) — a specific, measurable
   diagnostic, and the rodeo filter removes it.
5. **Rodeo spectroscopy of the exit channel (the new use).** After the pulse
   [H_P, H_T] = 0, so filter both fragments jointly: a 2D scan in (E_T, E_P)
   whose peak heights converge to the exclusive probabilities |c_{κab}|²;
   reading M_T by Z in the same shot assigns κ and breaks accidental
   cross-sector degeneracies. Gives: (a) g.s.→g.s. channel separately from the
   inclusive P_{+1} compared with Sn+Ni data; (b) the 0⁺ strength in 42Ca that
   the inclusive Zr+Ca data lump together (the draft had to put the +1 centroid
   at Q=0 by hand); (c) rare channels measured through the plateau of
   P_M vs M — background 2^−M(1−p) below p at M≈6 for p≈0.02, below 1% of p at
   M≈13; M≈19 for P_{+2}≈2e-4 (shot overhead 1/p, unavoidable).
6. **A two-sided, device-measurable bracket on the entanglement entropy:**
   S_num ≤ S_{T|P} ≤ S_diag, with S_diag the Shannon entropy of the diagonal of
   ρ_T in the H_T eigenbasis — which is exactly what the H_T-only rodeo scan
   measures (P_{κa} = ⟨E_T^a|ρ_T|E_T^a⟩). Lower bound from Z readout, upper
   bound from the same rodeo circuit. Since S_cfg is tiny the bracket should
   be tight; where it isn't, its width measures intra-sector coherence.
   (Check literature on diagonal/observational entropy bounds before claiming
   novelty.)
7. **Adiabatic bound on the non-elastic probability.** The pulse starts and
   ends at H_P+H_T, so the JRS bound paper 1 evaluated for the refinement path
   bounds 1−P_000 ≈ 0.03 for the collision. Expected to be loose (paper 1:
   2–3 orders), but a clean link between the two papers.
8. Optional: static di-nuclear ground state at frozen peak coupling via
   "coupling refinement" H(s) = H_P + H_T + s·G_pk V_PT (P = 1, two-local).

## 3. What is genuinely new relative to paper 1 (so it is not a rerun)

realistic spectra and Ω_p > 1 · two coupled registers · dynamics ·
rodeo as exit-channel spectroscopy · entropy bracket · two-local encoding
(vs 325 JW strings of weight ≤ 8 at 8 qubits) · experiment-anchored observable
P_2n(D) · Path B.

## 4. Risks / open questions to settle numerically first

- Overlap of the filled Sn seed (1g7/2+3s1/2) with GS: 3s1/2 is only 48%
  occupied and 2d3/2 29%, so the first projectile step may carry a low overlap
  → long T* or an extra chain step. (Paper 1's k=2→4 had 0.868.)
- Minimum gaps along the paths: set by orbital spacings, not pairing — could
  make refinement *fast* in ħ/MeV, but near-degeneracies among 154 states in
  Sn at step 3 need checking.
- Path B: dependence on Λ and on the schedule; initial infidelity O((GΩ/Λ)²).
- Exit-channel scan: spectral density of H_P+H_T over the M=7 space (5554
  product levels); σ must resolve the smallest relevant spacing (Q=1.30 MeV
  between sector ground states; pairing excitations ~ MeV). M_T readout is
  the safety net.
- Shot overhead for the +2 channel (≈5e3 shots per accepted event).
- Whether to keep a variational comparator (pUCCD in XY form). Suggest yes,
  briefly, in the matched-observable style of paper 1 — but the pitch is
  "no optimizer anywhere between the filled configuration and the readout".
- Author list / ordering across the two author sets; journal (PRC again?).

## 5. Suggested order of calculations (when we start)

1. Generalize `pairinglib` to Ω_p > 1, arbitrary ε_p, two registers,
   time-dependent J_ij; symmetric-subspace exact solver (dims 9/154/5554) +
   fixed-M slot-sector statevector simulator (C(19,5)=11628 for Sn,
   C(30,7)=2.0e6 for the collision). Reproduce the draft's reference numbers.
2. 60Ni chain (11 qubits): spectra, gaps, overlaps, T* for Paths A and B,
   rodeo 1−F vs M. Cheapest complete refinement study.
3. 116Sn chain (19 qubits): same; decide chain granularity from the overlaps.
4. 96Zr+40Ca end-to-end (26 qubits, target starts in |0…0⟩, only Zr needs
   preparation): Trotterized pulse vs DOP853; (E_T,E_P) scan; 42Ca 0⁺ strength;
   entropy bracket. This is the natural first full demonstration.
5. 116Sn+60Ni end-to-end (30 qubits): P_κ(D) overlay on the draft's Fig. 2,
   exclusive vs inclusive, bracket vs D, budgets.
6. JRS bounds (paths and pulse), pUCCD comparator, gate-budget tables.

Figures planned (8) and tables (3) are listed in the skeleton's Results section.
