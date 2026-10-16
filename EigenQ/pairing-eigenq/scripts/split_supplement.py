#!/usr/bin/env python3
"""One-off restructuring of paper I (2026-10-16): reframe the manuscript around
(a) the validated reference implementation, (b) the matched-accuracy benchmark
and (c) the three transferable results; move verification material to a
Supplemental Material file; remove the references to the companion papers.
Input: paper/eigenq_pairing_v7_backup.tex (the pre-split manuscript).
Output: paper/eigenq_pairing.tex (main text) and paper/eigenq_pairing_supp.tex.
Run from the repo root: python3 scripts/split_supplement.py
"""
import re, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
src = (ROOT / "paper/eigenq_pairing_v7_backup.tex").read_text()
s = src


def cut_env(text, label):
    """Remove the float environment carrying `label`; return (text, block)."""
    for m in re.finditer(r"\\begin\{(figure\*?|table\*?)\}.*?\\end\{\1\}\n?", text, re.S):
        if f"\\label{{{label}}}" in m.group(0):
            return text[:m.start()] + text[m.end():], m.group(0)
    raise KeyError(label)


def cut_between(text, start, end):
    i = text.index(start); j = text.index(end, i)
    return text[:i] + text[j:], text[i:j]


def rep(text, a, b, n=1):
    assert text.count(a) == n, (a[:70], text.count(a))
    return text.replace(a, b)


moved = {}
for lab in ("fig:methods", "fig:scaling", "fig:fourparticle", "fig:chains", "fig:assess",
            "fig:scan", "tab:ops", "fig:large", "tab:large"):
    s, moved[lab] = cut_env(s, lab)
s, moved["sec:res-bound"] = cut_between(s, "\\subsection{Sufficient conditions versus practice}",
                                        "\\subsection{Rodeo filtering of the refined state}")
s, moved["sec:sizes"] = cut_between(s, "\\subsection{Other model-space sizes, particle numbers, and observables}",
                                    "\\subsection{A larger model space: four pairs on eight levels}")

# ------------------------------------------------------------------ header, title, macros
s = rep(s, """%  Paper I of two: the pure (seniority-conserving) pairing model.  The
%  seniority-breaking particle-hole extension is the companion paper
%  paper_ph/eigenq_pairing_ph.tex (cited as \\cite{PHcompanion}).""",
"""%  Main text; the Supplemental Material is eigenq_pairing_supp.tex (built with
%  xr against this file's .aux).  Restructured 2026-10-16 around the reference
%  implementation, the matched-accuracy benchmark and the transferable results;
%  the pre-split manuscript is eigenq_pairing_v7_backup.tex.""")
s = rep(s, """\\title{Resolution refinement followed by the rodeo algorithm for
eigenstate preparation: implementation and assessment on the pairing Hamiltonian}""",
"""\\title{Resolution refinement followed by the rodeo algorithm: a gate-level
benchmark of eigenstate-preparation routes on the pairing Hamiltonian}""")
s = rep(s, "\\newcommand{\\Eobj}{E_{\\mathrm{obj}}}\n",
        "\\newcommand{\\Eobj}{E_{\\mathrm{obj}}}\n\\newcommand{\\SM}{Supplemental Material~\\cite{SM}}\n")

# ------------------------------------------------------------------ abstract
i = s.index("\\begin{abstract}"); j = s.index("\\end{abstract}") + len("\\end{abstract}")
s = s[:i] + r"""\begin{abstract}
\textbf{Background:} Projective spectral filters such as the rodeo algorithm
prepare eigenstates with an accuracy that improves exponentially in the
number of cycles, at a cost set by the squared overlap $p$ of the input
with the target.  Several routes can supply that input --- variational
circuits, an ordinary adiabatic ramp, or resolution refinement, the
adiabatic transport of an eigenstate from a coarse to a fine
single-particle basis along a shift-controlled path --- and they have not
been compared on equal terms.
\textbf{Purpose:} To provide a validated, gate-counted reference
implementation of resolution refinement followed by the rodeo algorithm on
the exactly solvable constant-pairing Hamiltonian, and to benchmark it at
matched accuracy against the simplest alternative inputs to the same
filter: the Hartree--Fock determinant, restricted and full unitary
coupled-cluster singles-and-doubles (UCCSD) states, and an adiabatic ramp
from the Hartree--Fock Hamiltonian.
\textbf{Methods:} All stages are compiled to single-qubit rotations and
controlled-NOT (CNOT) gates, with second-order product formulas for the
refinement and for the controlled evolution of the filter, and are simulated
exactly within the particle-number sector against full configuration
interaction, coupled-cluster doubles, UCCSD and the analytic
success-probability formulas of the filter.  The figure of merit is the
ground-state infidelity of the accepted ensemble against the CNOT count per
accepted preparation, minimized over the number of filter cycles, for four
particles on four levels (eight qubits) and four pairs on eight levels
(sixteen qubits).
\textbf{Results:} Every stage behaves as predicted: refinement raises the
overlap of the prolonged determinant from $0.87$ to $0.998$ along a gapped
path with no trained parameter, and the filter removes the remainder
geometrically, in agreement with the exact Gaussian-ensemble expression
derived here.  Three results transfer beyond the model.  The response norm
of the refinement path at its start is $P(k-P)g^2/(4b^2)$ for $P$ pairs on
$k$ levels and shift $b$, which fixes how the shift must scale with size and
coupling.  The long-run state and acceptance of the digital filter are the
dominant eigenoperator of an averaged post-selection map, which sets the
accuracy floor ($1.2\times10^{-4}$ at eight qubits) and rejects a too-coarse
product formula outright.  The available sufficient adiabatic-time and
Trotter-step conditions overestimate the empirical cost by two to three
orders of magnitude.  At eight qubits a four-parameter pair-rotation ansatz
of $192$ CNOTs reaches an infidelity of $1.8\times10^{-4}$ without any
filter and the ramp beats refinement by a factor $2.5$--$3$ at every target;
at sixteen qubits the ansatz degrades to $3\times10^{-3}$, filter cycles
become eleven times dearer, and refinement plus filter is $1.5$--$2$ times
cheaper than the variational route while staying within $4$--$19\%$ of the
ramp.
\textbf{Conclusions:} Resolution refinement followed by the rodeo algorithm
is a validated, optimizer-free and measurement-free route to a
high-fidelity correlated eigenstate.  On the pairing model it is never the
cheapest route --- the adiabatic ramp from the same determinant is, by the
cost of the embedded coarse term --- but it overtakes the variational
preparation as the space grows.  The benchmark, its baselines and its shift
rule fix the terms on which the method's advantages for larger spaces and
poorer variational inputs must be demonstrated.
\end{abstract}""" + s[j:]

# ------------------------------------------------------------------ introduction
i = s.index("Relative to the resolution-refinement proposal of Ref.~\\cite{Bogner2026}, the")
j = s.index("The pairing force conserves seniority, and that symmetry is what keeps the")
s = s[:i] + r"""Relative to the resolution-refinement proposal of Ref.~\cite{Bogner2026}, the
rodeo literature~\cite{Choi2021,Qian2024} and the fusion
method~\cite{Patkowski2025fusion}, which established adiabatic preconditioning
of the filter on lattice and spin systems, the present work contributes
three things.  (a) A reference implementation: refinement for a hierarchy of
\emph{single-particle bases} in a fermionic pairing problem, where the
prolongation is gate-free and the embedded coarse Hamiltonian is a conditioned
operator whose circuit we give and count, with every stage --- Trotterized
refinement and Trotterized controlled evolution included --- compiled to
rotations and CNOTs, validated against exact results at eight and at sixteen
qubits, and reproducible from a single executable notebook.  (b) A benchmark
on equal terms: the refinement-plus-filter chain is compared with the same
filter applied to the unrefined input, to an ordinary adiabatic ramp, to a
prolonged coarse variational state and to full-space variational states ---
the UCCSD template and the restricted pair-rotation ansatz that reproduces
it --- with the ground-state infidelity as the common observable, matched
evolution accuracy, a two-qubit-gate budget per accepted preparation that
separates one-time training, repeated preparation and the filter, and the
minimum cost at a given target accuracy, filter cycles included or not.
(c) Three results that transfer beyond the model: the closed-form response
norm of the refinement path, $P(k-P)g^2/(4b^2)$, which fixes how the shift
of the path must scale with the size and the coupling; the characterization
of the digital filter's long-run state and acceptance by the dominant
eigenoperator of an averaged post-selection map, which sets its accuracy
floor and explains why a too-coarse product formula is rejected; and the
evaluation of the rigorous sufficient conditions on the adiabatic time and
on the Trotter step count of Ref.~\cite{Sadia2026} against the resources
actually needed.  The outcome of (b) is stated plainly.  At eight qubits
refinement followed by filtering does not outperform filtering a variational
input, is roughly tied with filtering the unrefined prolonged state at
matched fidelity, and costs $2.5$--$3$ times more than an ordinary adiabatic
ramp from the Hartree--Fock Hamiltonian.  At sixteen qubits
(Sec.~\ref{sec:large}) the variational input degrades and the filter cycles
become expensive, so refinement plus filter is then $1.5$--$2$ times cheaper
than the variational route --- but the plain ramp from the same determinant
stays cheaper by the cost of the embedded coarse term.  The case for the
refinement stage rests on its freedom from classical optimization and
training measurements, and on coarse eigenstates that are better than the
determinant, which this hierarchy does not exercise.

""" + s[j:]
s = rep(s, """throughout the present paper, as does most of the quantum-algorithm
literature on pairing cited above.  In a companion paper~\\cite{PHcompanion}
we remove the simplification by adding a seniority-breaking particle-hole
interaction and repeat the complete assessment --- exact and variational
references, refinement, filtering and the circuit-level budget --- on the
resulting Hamiltonian, which carries the one-particle--one-hole channel of
realistic two-body interactions.""", """throughout the present paper, as does most of the quantum-algorithm
literature on pairing cited above; a seniority-breaking interaction, which
opens the one-particle--one-hole channel of realistic two-body forces and
degrades the single-reference variational input, is the natural next test
of the benchmark and is discussed in Sec.~\\ref{sec:discussion}.""")

# ------------------------------------------------------------------ pipeline section: pointers to the SM
s = rep(s, """commutator norms over all pairs of groups --- the natural extension of the
two-group statement of Ref.~\\cite{Sadia2026}, which is what
Sec.~\\ref{sec:res-bound} evaluates.""", """commutator norms over all pairs of groups --- the natural extension of the
two-group statement of Ref.~\\cite{Sadia2026}.""")
s = rep(s, """coarse-space estimate of $\\chi$, and we have not demonstrated that.
Sec.~\\ref{sec:res-bound} evaluates both conditions for our path and
quantifies the overshoot.""", """coarse-space estimate of $\\chi$, and we have not demonstrated that.
Both conditions are evaluated for our path in Sec.~S4 of the \\SM{} and
summarized in Sec.~\\ref{sec:res-gscan}.""")

# ------------------------------------------------------------------ results: coarse solver
s = rep(s, """Figure~\\ref{fig:methods} compares all solvers (exact, CCD, statevector and gate-level
UCCSD, and a hardware-efficient ansatz~\\cite{Kandala2017}) as a function of $k$.  The hardware-efficient
ansatz is exact at $k=2$ and reaches FCI at $k=3$ only after a multi-start
optimization, and struggles at $k=4$ --- the symmetry-blind landscape that motivates
the physically structured UCCSD circuit.  Figure~\\ref{fig:scaling} shows that the same
UCCSD solver remains a tight variational bound all the way to $k=8$ (sixteen qubits),
at the cost of a CNOT count that grows steeply ($1312$ CNOTs already at $k=4$); the
compiled circuits and their depths are analyzed in Sec.~\\ref{sec:circuits}.""",
"""Figures~S1 and S2 of the \\SM{} compare all solvers (exact, CCD, statevector
and gate-level UCCSD, and a hardware-efficient ansatz~\\cite{Kandala2017}) as
a function of $k$: the hardware-efficient ansatz is exact at $k=2$, reaches
FCI at $k=3$ only after a multi-start optimization and fails at $k=4$ --- the
symmetry-blind landscape that motivates the physically structured UCCSD
circuit --- while the UCCSD solver remains a tight variational bound all the
way to $k=8$ (sixteen qubits), at the cost of a CNOT count that grows steeply
($1312$ CNOTs already at $k=4$; Sec.~\\ref{sec:circuits}).""")
s = rep(s, "$0.635548$) at $T=30$ (Fig.~\\ref{fig:fourparticle}).", "$0.635548$) at $T=30$ (Fig.~S3 of the \\SM{}).")
s = rep(s, "quantity a device can know.  Figure~\\ref{fig:chains} compares the two chains.",
        "quantity a device can know.  Figure~S4 of the \\SM{} compares the two chains.")
s = rep(s, "the limitation explicit (Fig.~\\ref{fig:assess}, right).", "the limitation explicit (Fig.~S5 of the \\SM{}).")

# the sufficient-conditions subsection is replaced by a summary paragraph at the end of the coupling-scan subsection
s = rep(s, """is $2.0$, refines the same determinant in a quarter of the time.

""", """is $2.0$, refines the same determinant in a quarter of the time.

\\paragraph{Sufficient conditions versus practice.}
The conditions of Sec.~\\ref{sec:refine} can be evaluated exactly for the
path of Eq.~\\eqref{eq:path}, since $\\dot H$ and $\\ddot H$ are both
proportional to $B-A$ with $A=P(\\Hlow-\\mu)P^\\dagger$ and $B=\\Hhigh-\\mu$;
Sec.~S4 of the \\SM{} gives the integrand, the gap along the path and the
numbers.  With $\\norm{B-A}=7.62$ and an integral of $511$ in
Eq.~\\eqref{eq:Tbound}, the time condition demands $T\\ge5.1\\times10^{3}$
for an infidelity of $10^{-2}$ where the computed overlap reaches $0.99$ at
$T^{*}=16$, an overshoot of $3\\times10^{2}$; the first-order step condition
demands $3.4\\times10^{4}$ steps along the $T^{*}$ path where $160$
first-order or $40$ second-order steps reproduce the exact-slice overlap, a
factor $2\\times10^{2}$.  The two step counts answer different questions ---
the bound certifies an operator-norm error of the digitized unitary, one
reproduced overlap does not --- and the overshoot has an identifiable
source: the bound pays for the full norm $\\norm{B-A}$ and the global gap,
whereas the adiabatic error is governed by the individual derivative
couplings of the tracked state, which at $s=0.25$ are $0.60$, $0.37$,
$0.00$, $0.11$ and $0.08$ to the five accessible excited states.  These are
the fermionic counterpart of the ratios Ref.~\\cite{Sadia2026} reports for
nucleation of an Ising chain ($6\\times10^{1}$ and $10^{7}$): the conditions
certify convergence but do not predict cost, and the refinement times and
step counts used below are therefore read off the computed overlap with
the exact ground state.

""")

# ------------------------------------------------------------------ results: rodeo (energy scan shortened)
s = rep(s, """product-formula evolution, and Table~\\ref{tab:ops} keeps the two apart.
Scanning the target energy $E$ (Fig.~\\ref{fig:scan}; $M=4$ cycles, $12$ time sets
per energy) produces a peak at the ground-state energy; a coarse scan
($\\sigma=2$, $60$ energies on $[-0.5,4]$) resolves it from the excited band,
and a fine scan ($\\sigma=8$, $161$ energies with spacing $0.0075$ on a
window of width $1.2$ centered on the highest point of the coarse scan,
$E=0.644$, i.e.\\ chosen without knowledge of $E_0$) locates it at
$E=0.6362$ against the exact $E_0=0.635548$ by three-point parabolic
interpolation at the grid maximum.  On a device every acceptance is a shot
estimate: adding binomial noise for $10^{3}$ and $10^{4}$ shots per energy
and locating the peak by a quadratic fit through the points above $80\\%$ of
the maximum gives $0.637\\pm0.001$ and $0.6373\\pm0.0004$, so with this scan
the systematic error of the peak estimator ($\\sim2\\times10^{-3}$, from the
non-parabolic shape of the peak over the fitted window) exceeds the shot
noise, and the agreement to $10^{-3}$ is a property of the scan design
rather than a demonstrated four-decimal precision.  Starting from the
refined state, essentially all spectral weight sits in the single
ground-state peak; starting from $\\Phi_0$, the small residual weight on
excited states is visible as secondary structure above the $2^{-M}$
background.  The scan is part of the cost of the method: the fine scan alone
evaluates $161$ energies with $12$ time sets each, i.e.\\ about $2\\times10^{3}$
four-cycle runs, each repeated over enough shots to estimate an acceptance
probability.""", """product-formula evolution, and Table~S2 of the \\SM{} keeps the two apart.
Scanning the target energy $E$ ($M=4$ cycles, $12$ time sets per energy;
Fig.~S7 and Sec.~S5 of the \\SM{}) produces a peak at the ground-state
energy: a coarse scan ($\\sigma=2$) resolves it from the excited band and a
fine scan ($\\sigma=8$, centered on the coarse-scan maximum, i.e.\\ chosen
without knowledge of $E_0$) locates it at $E=0.6362$ against the exact
$0.635548$; with $10^{3}$--$10^{4}$ shots per energy the peak estimator
returns $0.637\\pm0.001$, its systematic error exceeding the shot noise.  The
scan is part of the cost of the method --- about $2\\times10^{3}$ four-cycle
runs for the fine scan alone, each repeated over enough shots to estimate an
acceptance probability.""")

# ------------------------------------------------------------------ circuits: HEA remark, operating points, sizes
s = rep(s, "ground state (Fig.~\\ref{fig:methods}).  Depth is not the only currency:",
        "ground state (Fig.~S1 of the \\SM{}).  Depth is not the only currency:")
s = rep(s, """Table~\\ref{tab:ops} records the operating points behind the different
fidelities quoted in this paper, which are not contradictory because they
belong to different stages of idealization:""", """Table~S2 of the \\SM{} records the operating points behind the different
fidelities quoted in this paper, which are not contradictory because they
belong to different stages of idealization:""")
s = rep(s, """For the pairing model with $N=4$ that
regime is not reached at any size we can validate exactly
(Sec.~\\ref{sec:sizes}); with four pairs on eight levels it is
(Sec.~\\ref{sec:large}),""", """For the pairing model with $N=4$ that
regime is not reached at $k=4$ or $k=6$, nor at $N=6$ (Sec.~S7 of the
\\SM{}); with four pairs on eight levels it is (Sec.~\\ref{sec:large}),""")
s = rep(s, """The sizes of Table~\\ref{tab:sizes} keep the variational input good and the
filter cycles cheap.""", """The sizes of Table~\\ref{tab:depth}, and the $N=4$ and $N=6$ systems on
$k=4$ and $6$ levels of Sec.~S7 of the \\SM{} (where the ranking is the
same), keep the variational input good and the filter cycles cheap.""")
s = rep(s, """computation can be carried out in the pair basis, one bit per level, which we
validated against the full-sector implementation at $k=4$, $N=4$: the
Hamiltonian and every Trotter term coincide with the seniority-zero blocks to
machine precision, and every anchor of this paper (the energy
$0.6355484736$, $p_0=0.8676375677$, $p(25)=0.9975712826$, the Trotterized
$p=0.9922$ at $T^*=16$, the four-pair ansatz $E=0.6369869810$ and
$1-F=1.776\\times10^{-4}$, the exact-ensemble infidelities $1.216\\times10^{-3}$,
$6.08\\times10^{-4}$, $3.80\\times10^{-5}$ and $2.38\\times10^{-6}$) is
reproduced.""", """computation can be carried out in the pair basis, one bit per level, which we
validated against the full-sector implementation at $k=4$, $N=4$ (Sec.~S8
of the \\SM{}): the Hamiltonian and every Trotter term coincide with the
seniority-zero blocks to machine precision and every anchor of this paper is
reproduced.""")
s = rep(s, "Table~\\ref{tab:large} lists the ingredients.  The HF overlap falls from",
        "Table~S5 of the \\SM{} lists the ingredients.  The HF overlap falls from")
s = rep(s, "one curve (Fig.~\\ref{fig:large}, left): the transport is controlled by the",
        "one curve (Fig.~S8 of the \\SM{}, left): the transport is controlled by the")
s = rep(s, "coupling, so it does not slow down with $g$ (Fig.~\\ref{fig:large}, right),",
        "coupling, so it does not slow down with $g$ (Fig.~S8 of the \\SM{}, right),")

# ------------------------------------------------------------------ discussion
i = s.index("The pairing model lets us exhibit the complete refinement-plus-filter chain")
j = s.index("\\paragraph{Relation to symmetry-restoration and subspace approaches.}")
s = s[:i] + r"""What this work delivers is, first, a reference implementation: the
complete refinement-plus-filter chain with every stage benchmarked against
exact results --- a gate-level UCCSD coarse solver, a gate-free prolongation
with a counted conditioned coarse term, a gapped adiabatic refinement, and a
rodeo filter whose ensemble infidelity follows Eq.~\eqref{eq:fid} --- all
validated a second time at the circuit level with product formulas, at
eight and at sixteen qubits.  The two algorithms are complementary in the
intended sense: refinement raises the overlap from $0.87$ to $0.998$, and
the filter removes the adiabatic tail geometrically.  The register and the
circuits of all stages --- prolongation, the adiabatic path, the shift rule
and the rodeo cycle --- are written for the generic $2k$-qubit register and
do not exploit seniority for compression; the favorable dynamics and the
small variational inputs available here do rely on the seniority
conservation of the pairing force.

Second, three results that do not depend on the model.  The response norm
of the refinement path at its start, $\chi_u(0)=P(k-P)g^2/(4b^2)$, says how
the shift must grow with the size and the coupling for the transport to
keep its speed --- the shift that is right at $k=4$ is six times too slow at
$k=8$ --- and shows why a scan of the refinement time against the
shift-controlled gap cannot identify a scaling law, since the scan varies
the quantity that the inverse-gap argument assumes bounded.  The digital
filter's long-run state and acceptance are the dominant eigenoperator and
eigenvalue of the averaged post-selection map, which sets the accuracy
floor ($1.2\times10^{-4}$ at $\delta t=0.25$ for eight qubits,
$4.6\times10^{-5}$ at $\delta t=0.1$ for sixteen) and the slow geometric
decay of the acceptance, and whose eigenoperator loses the ground state
altogether when the step is too coarse for the Hamiltonian norm, so that the
filter rejects the corrupted evolution --- a useful diagnostic, not a
guarantee against smaller coherent errors.  And the rigorous sufficient
conditions on the adiabatic time and the Trotter step count overestimate the
resources actually needed by two to three orders of magnitude, for an
identifiable reason: they pay for the global norm and gap where the dynamics
pays only for the couplings of the tracked state.

Third, and this is the finding the benchmark was built for, the ranking of
the routes.  Compared on the same observable --- the ground-state
infidelity --- and with the same Trotterized filter, refinement followed by
the rodeo algorithm is, at $N=4$ and $N=6$ on four and six levels, tied with
filtering the unrefined prolonged determinant, costs $2.5$--$3$ times more
than an ordinary adiabatic ramp from the Hartree--Fock Hamiltonian at every
target accuracy, and is beaten by two orders of magnitude by the
four-parameter pair-rotation ansatz, which needs no filter at all
(Table~\ref{tab:depth} and Sec.~S7 of the \SM{}).  At four pairs on eight
levels (Sec.~\ref{sec:large}) the variational ansatz is no longer good, the
cycles are eleven times dearer, and refinement plus filter overtakes both
the variational route and the filtered determinant by factors $1.5$--$2$
and more; but the ordinary ramp from the same determinant stays ahead by
the cost of the conditioned coarse term, at $g=1$ and at $g=2$, once the
shift is scaled with the size.  The reason is simple, and
Sec.~\ref{sec:depth} states it as the ideal large-gap law
$1-\mathcal F_M\simeq\tfrac{1-p}{p}2^{-M}$: a better input saves cycles only
logarithmically in its odds ratio, so whether a preparation stage repays its
gates depends on the price of a cycle relative to the price of the stage,
and the size dependence of the verdict follows from the size dependence of
that ratio.  What the refinement offers in exchange is that it is
optimizer-free and measurement-free during preparation, and that its cost is
set by the Hamiltonian evolution rather than by the trainability of a
landscape.  Whether that trade pays off against a variational preparation
is, as Sec.~\ref{sec:large} shows, a question of size and coupling, and the
answer turns in favor of the adiabatic routes already at sixteen qubits;
whether the \emph{shift-controlled} path pays off against the plain ramp is
a question about coarse eigenstates that are better than the determinant,
which a hierarchy seeded by the determinant cannot answer, and about
hardware on which the training measurements of a variational loop
dominate.  We regard the present work as the controlled baseline against
which such claims should be measured, not as evidence for them.

\paragraph{Limitations.}
Four limitations bound the claims.  The comparison is at the level of
coherent two-qubit gate counts for an ideal device; no noise model enters,
and the acceptance statistics of the filter, which diagnose a corrupted
evolution, are the only place where an imperfect device is visible.  The
costs are per accepted preparation and exclude the one-time costs --- the
training measurements of the variational inputs, the energy scan of the
filter, and the calibration of the refinement controls --- which have no
common shot model and are amortized over the number of states requested.
The refinement controls $(T,\mu,n_s)$, the filter step $\delta t$ and its
target energy are calibrated here against the exact answer; an
implementable protocol would set $\mu$ from the coarse spectrum and the
closed-form response norm, read $T$ and $n_s$ from the convergence of the
measured energy, and locate the target by the scan, with the acceptance as
the check, and we have not demonstrated it.  And the variational comparator
is a restricted ansatz chosen by the symmetry selection rule, not an
adaptive construction; an ADAPT-VQE run proper is the next comparator to
add.

""" + s[j:]
s = rep(s, """the pairing force conserves seniority; the particle-hole term of the
companion paper~\\cite{PHcompanion}, which is not part of the model of
Ref.~\\cite{AychetClaisse2026}, is the ingredient that forces one out of it
and onto the full $2k$-qubit register used throughout this work.""", """the pairing force conserves seniority; a seniority-breaking particle-hole
term, which is not part of the model of Ref.~\\cite{AychetClaisse2026}, is
the ingredient that would force one out of it and onto the full $2k$-qubit
register used throughout this work.""")
s = rep(s, """against bounds of $T\\ge2.6\\times10^{3}$ and $N\\ge2\\times10^{10}$.  Section~\\ref{sec:res-bound} makes the same comparison for the refinement
path: the time bound of Eq.~\\eqref{eq:Tbound} overshoots the empirical
$T^{*}$ by a factor $3\\times10^{2}$ and the step bound overshoots
the steps actually needed by a factor $2\\times10^{2}$
(Table~\\ref{tab:bound}), because""", """against bounds of $T\\ge2.6\\times10^{3}$ and $N\\ge2\\times10^{10}$.  Section~\\ref{sec:res-gscan} and Sec.~S4 of the \\SM{} make the same comparison for the refinement
path: the time bound of Eq.~\\eqref{eq:Tbound} overshoots the empirical
$T^{*}$ by a factor $3\\times10^{2}$ and the step bound overshoots
the steps actually needed by a factor $2\\times10^{2}$, because""")
s = rep(s, """with one cycle to $10^{5}$ for the refined input with six cycles (the companion paper~\\cite{PHcompanion} gives the
corresponding budget with the particle-hole term), on top of the energy scan""",
"""with one cycle to $10^{5}$ for the refined input with six cycles at eight
qubits, and ten times these at sixteen, on top of the energy scan""")
s = rep(s, """treated with a nontrivial coarse space --- which is where the trainability
argument for refinement would be tested rather than asserted; the companion
paper~\\cite{PHcompanion} takes the first step in this direction by adding a
seniority-breaking particle-hole interaction, which degrades the
single-reference ansatz with the strength of the pair breaking.  Third,""",
"""treated with a nontrivial coarse space --- which is where the trainability
argument for refinement would be tested rather than asserted; the first
step in this direction is a seniority-breaking particle-hole interaction,
which degrades the single-reference ansatz with the strength of the pair
breaking and forces the full $2k$-qubit register.  Third,""")
s = rep(s, """cost that, for the pairing model at these sizes, exceeds that of filtering a
cheap variational input.""", """cost that, for the pairing model, exceeds that of an adiabatic ramp from
the same determinant at every size studied and that of a cheap variational
input at eight but not at sixteen qubits.""")

# ------------------------------------------------------------------ code availability, bibliography
s = rep(s, """is provided as a tested Python package and a single Jupyter notebook that
reproduces every figure in this article and in the companion
paper~\\cite{PHcompanion}.
Give GitHub link when everything is ready.""", """is provided as a tested Python package and a single Jupyter notebook that
reproduces every figure and number in this article and in its
\\SM{}~\\cite{Code}.""")
i = s.index("\\bibitem{PHcompanion}"); j = s.index("\\bibitem", i + 1)
s = s[:i] + s[j:]
s = rep(s, "\\end{thebibliography}", """\\bibitem{SM} See Supplemental Material at [URL will be inserted by
publisher] for the coarse-solver comparison at larger basis size, the
refinement energies and the $3\\to4$ chain, the shift-scan fits and the
response norm along the path, the evaluation of the sufficient
adiabatic-time and Trotter-step conditions, the rodeo energy scan and its
statistics, the operating points behind the quoted fidelities, the other
model-space sizes, particle numbers and observables, and the validation of
the pair-basis implementation used at sixteen qubits.
\\bibitem{Code} Repository \\texttt{pairing-eigenq}, [URL to be inserted];
archived version: Zenodo, [DOI to be inserted].
\\end{thebibliography}""")
assert "PHcompanion" not in s and "companion" not in s.lower(), [l for l in s.splitlines() if "companion" in l.lower()]
(ROOT / "paper/eigenq_pairing.tex").write_text(s)

# ================================================================== Supplemental Material
head = r"""\documentclass[aps,prc,reprint,superscriptaddress,floatfix,nofootinbib]{revtex4-2}
% Supplemental Material for "Resolution refinement followed by the rodeo
% algorithm: a gate-level benchmark of eigenstate-preparation routes on the
% pairing Hamiltonian".  Built by scripts/split_supplement.py; cross-references
% to the main text resolve through xr against eigenq_pairing.aux (compile the
% main text first).
\usepackage{amsmath,amssymb}
\usepackage{physics}
\usepackage{graphicx}
\usepackage{hyperref}
\usepackage{xr}
\externaldocument{eigenq_pairing}
\newcommand{\Hlow}{H_{\mathrm{low}}}
\newcommand{\Hhigh}{H_{\mathrm{high}}}
\newcommand{\Eobj}{E_{\mathrm{obj}}}
\renewcommand{\thesection}{S\arabic{section}}
\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\theequation}{S\arabic{equation}}

\begin{document}
\title{Supplemental Material for ``Resolution refinement followed by the
rodeo algorithm: a gate-level benchmark of eigenstate-preparation routes on
the pairing Hamiltonian''}
"""
i = src.index("\\author{Scott Bogner}"); j = src.index("\\date{\\today}")
authors = src[i:j]
head += authors + "\\date{\\today}\n\\maketitle\n\n"

body = r"""
This Supplemental Material collects the verification material behind the
main text: the comparison of coarse solvers at larger basis size
(Sec.~\ref{sm:solvers}), the refinement energies and the $k=3\to4$ chain
(Sec.~\ref{sm:refine}), the shift-scan fits and the response norm along the
path (Sec.~\ref{sm:shift}), the evaluation of the sufficient adiabatic-time
and Trotter-step conditions (Sec.~\ref{sm:bound}), the rodeo energy scan
and its statistics (Sec.~\ref{sm:scan}), the operating points behind the
fidelities quoted in the main text (Sec.~\ref{sm:ops}), the other
model-space sizes, particle numbers and observables (Sec.~\ref{sm:sizes}),
and the validation of the pair-basis implementation with the ingredients
of the sixteen-qubit study (Sec.~\ref{sm:large}).  Equation, section,
figure and table numbers without the prefix S refer to the main text.
All numbers are produced by the notebook of the code repository.

\section{Coarse solvers at larger basis size}
\label{sm:solvers}
Figure~\ref{fig:methods} compares all coarse solvers of Sec.~\ref{sec:res-uccsd}
of the main text --- exact diagonalization, coupled-cluster doubles,
statevector and gate-level UCCSD, and a hardware-efficient ansatz --- as a
function of the number of levels $k$ at $N=4$, $g=1$, and
Fig.~\ref{fig:scaling} follows the gate-level UCCSD solver to $k=8$ (sixteen
qubits).  The hardware-efficient ansatz is exact at $k=2$, reaches FCI at
$k=3$ only after a multi-start optimization with a particle-number penalty,
and fails at $k=4$; the UCCSD solver stays variational and tight throughout,
at a CNOT count that grows steeply ($1\,312$ at $k=4$).

""" + moved["fig:methods"] + "\n" + moved["fig:scaling"] + r"""
\section{Refinement energies and the $k=3\to4$ chain}
\label{sm:refine}
Figure~\ref{fig:fourparticle} shows the overlap and the physical energy
$\expval{\Hhigh}{\Phi(T)}$ of the refined state versus the adiabatic time
$T$ for the $k=2\to3$ and $2\to4$ chains at $N=4$, $g=1$
(Sec.~\ref{sec:res-refine} of the main text), and Fig.~\ref{fig:chains}
compares the $2\to4$ chain with the $3\to4$ chain that starts from the
gate-level UCCSD state at $k=3$ (Sec.~\ref{sec:res-k34}).

""" + moved["fig:fourparticle"] + "\n" + moved["fig:chains"] + r"""
\section{Shift scan: fits and response norm}
\label{sm:shift}
Figure~\ref{fig:assess} shows the two affine fits to the four shift-scan
points of Fig.~\ref{fig:gapscan} of the main text, in $1/\Delta E_{\min}$
and in $1/\Delta E_{\min}^2$, with their sums of squared residuals ($28$
and $0.2$), and the first-order response norm $\chi_u(s)$ along the
$k=2\to4$ path for the four shifts; at $s=0$ it equals $g^2/b^2$ with
$b=\mu-E_{\rm low}$, the $k=4$, $N=4$ case of the closed form
$P(k-P)g^2/(4b^2)$ of Sec.~\ref{sec:large} of the main text.

""" + moved["fig:assess"] + "\n"

bound = moved["sec:res-bound"]
bound = bound.replace("\\subsection{Sufficient conditions versus practice}\n\\label{sec:res-bound}",
                      "\\section{Sufficient conditions versus practice}\n\\label{sm:bound}")
bound = rep(bound, """two unitaries, are close.  (The companion paper~\\cite{PHcompanion} repeats this evaluation
with the particle-hole term switched on, where the time bound is insensitive
to the pair breaking because the minimum gap is set by the shift, while the
step bound grows with the number of Trotter groups.)""", "two unitaries, are close.")
bound = bound.replace("The conditions quoted in Sec.~\\ref{sec:refine} can be evaluated exactly",
                      "The conditions quoted in Sec.~\\ref{sec:refine} of the main text can be evaluated exactly")
bound = bound.replace("gate compilation of Sec.~\\ref{sec:trotterized} (the", "gate compilation of Sec.~\\ref{sec:trotterized} of the main text (the")
body += bound

body += r"""
\section{The rodeo energy scan}
\label{sm:scan}
Figure~\ref{fig:scan} shows the energy scans of Sec.~\ref{sec:res-rodeo} of
the main text ($M=4$ cycles, $12$ time sets per energy).  The coarse scan
($\sigma=2$, $60$ energies on $[-0.5,4]$) resolves the ground-state peak
from the excited band; the fine scan ($\sigma=8$, $161$ energies with
spacing $0.0075$ on a window of width $1.2$ centered on the highest point
of the coarse scan, $E=0.644$, i.e.\ chosen without knowledge of $E_0$)
locates the peak at $E=0.6362$ against the exact $E_0=0.635548$ by
three-point parabolic interpolation at the grid maximum.  On a device every
acceptance is a shot estimate: adding binomial noise for $10^{3}$ and
$10^{4}$ shots per energy and locating the peak by a quadratic fit through
the points above $80\%$ of the maximum gives $0.637\pm0.001$ and
$0.6373\pm0.0004$, so with this scan the systematic error of the peak
estimator ($\sim2\times10^{-3}$, from the non-parabolic shape of the peak
over the fitted window) exceeds the shot noise, and the agreement to
$10^{-3}$ is a property of the scan design rather than a demonstrated
four-decimal precision.  Starting from the refined state, essentially all
spectral weight sits in the single ground-state peak; starting from the
prolonged determinant, the small residual weight on excited states is
visible as secondary structure above the $2^{-M}$ background.  The fine
scan alone evaluates $161$ energies with $12$ time sets each, about
$2\times10^{3}$ four-cycle runs, each repeated over enough shots to estimate
an acceptance probability.

""" + moved["fig:scan"] + r"""
\section{Operating points}
\label{sm:ops}
Table~\ref{tab:ops} records the operating points behind the different
fidelities quoted in the main text, which belong to different stages of
idealization: the ideal filter of Sec.~\ref{sec:res-rodeo} and the digital
pipeline of Secs.~\ref{sec:trotterized}--\ref{sec:depth} differ in
refinement time, slice treatment, filter width and evolution, and the
residual errors of the former are not achieved by the counted circuits of
the latter.

""" + moved["tab:ops"] + "\n"

sizes = moved["sec:sizes"]
sizes = sizes.replace("\\subsection{Other model-space sizes, particle numbers, and observables}\n\\label{sec:sizes}",
                      "\\section{Other model-space sizes, particle numbers, and observables}\n\\label{sm:sizes}")
sizes = sizes.replace("To see whether the ranking of Table~\\ref{tab:depth} is an accident of the",
                      "To see whether the ranking of Table~\\ref{tab:depth} of the main text is an accident of the")
sizes = rep(sizes, """poor, and Sec.~\\ref{sec:large} shows that four pairs on eight levels is
already such a space.  The first step
beyond the pairing force --- opening the one-particle--one-hole channel, which
switches on the UCCSD singles and degrades the single-reference ansatz with
the strength of the pair breaking --- is taken in the companion
paper~\\cite{PHcompanion}.""", """poor, and Sec.~\\ref{sec:large} of the main text shows that four pairs on
eight levels is already such a space.""")
sizes = sizes.replace("correlation energy $E-E_{\\rm ref}$ of the states of Table~\\ref{tab:depth} at",
                      "correlation energy $E-E_{\\rm ref}$ of the states of Table~\\ref{tab:depth} of the main text at")
body += sizes

body += r"""
\section{The pair basis and the sixteen-qubit ingredients}
\label{sm:large}
Every stage of the pipeline conserves seniority --- the Hamiltonian, the
embedding, the path, the exactly exponentiated pair-hopping blocks and the
filter --- so for an $N$-paired target the whole computation can be carried
out in the pair basis (one bit per level, $\binom{k}{N/2}$ states).  The
pair-basis implementation used for Sec.~\ref{sec:large} of the main text
was validated against the full-sector code at $k=4$, $N=4$: the Hamiltonian
and every Trotter term coincide with the seniority-zero blocks of their
full-sector counterparts to machine precision, and every anchor of the main
text is reproduced --- the ground-state energy $0.6355484736$, the prolonged
overlap $p_0=0.8676375677$, the refined overlap $p(25)=0.9975712826$, the
Trotterized $p=0.9922$ at $T^*=16$ with $43$ second-order steps, the
four-pair ansatz $E=0.6369869810$ with $1-F=1.776\times10^{-4}$, and the
exact-ensemble infidelities $1.216\times10^{-3}$, $6.08\times10^{-4}$,
$3.80\times10^{-5}$ and $2.38\times10^{-6}$ after $1$, $2$, $6$ and $10$
cycles.  The gate counts of the sixteen-qubit study are those of the
$2k$-qubit compilation of Sec.~\ref{sec:circuits} of the main text.

Table~\ref{tab:large} lists the ingredients of the sixteen-qubit study for
$g=0.5$ to $2$, and Fig.~\ref{fig:large} the exact-slice transport along
the refinement chains $k_{\rm low}\to8$ and along the HF ramp, with its
coupling dependence.

""" + moved["tab:large"] + "\n" + moved["fig:large"] + "\n"

# bibliography of the SM: the items actually cited in the SM, copied from the main text
cited = sorted(set(re.findall(r"\\cite\{([^}]*)\}", body)))
keys = sorted({k.strip() for c in cited for k in c.split(",")})
bibitems = {}
for m in re.finditer(r"\\bibitem\{([^}]*)\}(.*?)(?=\\bibitem\{|\\end\{thebibliography\})", src, re.S):
    bibitems[m.group(1)] = m.group(0).rstrip() + "\n"
bib = "\\begin{thebibliography}{99}\n" + "".join(bibitems[k] for k in keys if k in bibitems) + "\\end{thebibliography}\n"
missing = [k for k in keys if k not in bibitems]
assert not missing, missing
(ROOT / "paper/eigenq_pairing_supp.tex").write_text(head + body + "\n" + bib + "\\end{document}\n")
print("main:", len(s.splitlines()), "lines; SM:", len((head + body + bib).splitlines()), "lines; SM cites", keys)


# ================================================================== second pass: trims (revision artifacts, duplicated prose)
s = (ROOT / "paper/eigenq_pairing.tex").read_text()
sm = (ROOT / "paper/eigenq_pairing_supp.tex").read_text()
s = rep(s, """ensemble values unless stated otherwise.  (The acceptance weighting itself
changes the sampled infidelity by at most $1/p-1=0.24\\%$ at this $p$; an
earlier version of the notebook advanced one generator through both inputs
and reported $5.5\\times10^{-6}$ at ten cycles, an equally valid sample.)  The""",
"""ensemble values unless stated otherwise (the acceptance weighting itself
changes the sampled infidelity by at most $1/p-1=0.24\\%$ at this $p$).  The""")
s = rep(s, """$0.998$ to $0.965$, $0.873$ and $0.596$ (exact ensemble values, Sec.~19 of
the notebook; an earlier version quoted sampled values) --- in""",
"""$0.998$ to $0.965$, $0.873$ and $0.596$ (exact ensemble values) --- in""")
s = rep(s, """than refining it first.  The cost ratios quoted in an earlier version of
this work at a fixed cycle number compared states at different achieved
errors and are superseded by these matched-accuracy figures; none of the
operating points is a globally optimized minimum, since refinement time,
shift and step size were not searched.  The claim that the hybrid undercuts
the variational preparation by orders of magnitude, which compared the UCCSD
\\emph{energy error} with the pipeline \\emph{infidelity}, does not survive the
comparison of like with like.""",
"""than refining it first.  Comparisons at a fixed cycle number would compare
states at different achieved errors, and a comparison of the UCCSD
\\emph{energy error} with the pipeline \\emph{infidelity} would favor the
pipeline by orders of magnitude; neither survives the comparison of like
with like, and none of the operating points is a globally optimized
minimum, since refinement time, shift and step size were not searched.""")
# HEA paragraph -> SM
i = s.index("\\paragraph{Hardware-efficient ansatz.}"); j = s.index("\\paragraph{Rodeo cycle.}")
hea = s[i:j]; s = s[:i] + s[j:]
hea = hea.replace("(Fig.~S1 of the \\SM{})", "(Fig.~\\ref{fig:methods})")
sm = rep(sm, "at a CNOT count that grows steeply ($1\\,312$ at $k=4$).\n", "at a CNOT count that grows steeply ($1\\,312$ at $k=4$).\n\n" + hea)
s = rep(s, """Figure~\\ref{fig:circ-vqe} shows both levels of this hierarchy.  Panel (a) is
the VQE at the ansatz level: the product of excitation unitaries acts on the HF
register, the energy is estimated from repeated measurements, and a classical
optimizer updates the $P$ amplitudes $\\boldsymbol\\theta$ before the circuit is
run again --- it is this feedback loop, wrapped around the full-depth state
preparation, that multiplies the $1312$-CNOT cost into the total gate budget of
the traditional route.  Panel (b) is the compiled interior of a single
weight-$4$ Pauli rotation, the sub-block that occurs eight times in every
pair-hopping double and $160$ times in the $k=4$ preparation; its only
parametrized element is the central $R_z(\\phi)$, with $\\phi$ fixed by the
amplitude of the parent excitation operator.""",
"""Figure~\\ref{fig:circ-vqe} shows both levels of this hierarchy: the VQE at
the ansatz level, whose classical feedback loop wrapped around the
full-depth preparation multiplies the $1312$-CNOT cost into the total gate
budget of the variational route, and the compiled interior of one
weight-$4$ Pauli rotation, the sub-block that occurs $160$ times in the
$k=4$ preparation.""")
s = rep(s, """The
contrast with Fig.~\\ref{fig:circ-vqe}(a) is the structural point of this
section: the rodeo circuit contains \\emph{no} variational parameters and closes
no optimization loop.  The times $t_k$ are drawn once from a Gaussian, the
target energy $E$ is a classical setting, and the algorithm is driven entirely
by the measurement record --- a run is accepted when all $M$ ancilla outcomes
read $0$, which applies the filter of Sec.~\\ref{sec:rodeo} cycle by cycle.
Where the VQE re-executes its full-depth circuit to steer $\\boldsymbol\\theta$,
the rodeo circuit is executed as printed, and its only system-size-dependent
element is the controlled evolution.  The controlled evolution is""",
"""The
rodeo circuit contains no variational parameters and closes no optimization
loop: the times $t_k$ are drawn once from a Gaussian, the target energy $E$
is a classical setting, and a run is accepted when all $M$ ancilla outcomes
read $0$, which applies the filter of Sec.~\\ref{sec:rodeo} cycle by cycle.
The controlled evolution is""")
# hardware paragraph shortened
i = s.index("\\paragraph{Hardware requirements.}"); j = s.index("\\paragraph{Concluding remarks and perspectives.}")
s = s[:i] + r"""\paragraph{Hardware requirements.}
The Trotterized pipeline quantifies the dominant coherent error of an ideal
device, and the requirements it implies are explicit: $2k+1$ qubits ($9$ at
$k=4$, $17$ at $k=8$); a mid-circuit measurement with reset and reuse of one
ancilla; arbitrary-angle $R_z$ rotations; and a coherent budget that ranges
from $192$ two-qubit gates for the cheapest route of Table~\ref{tab:depth}
through $1.7\times10^{4}$ per accepted preparation for the HF ramp with one
cycle to $10^{5}$ for the refined input with six cycles at eight qubits,
and ten times these at sixteen, on top of the energy scan and the
repetitions needed to estimate acceptance probabilities.  At the two-qubit
error rates of large-scale experiments today --- of order $10^{-2}$ on the
superconducting processor of Ref.~\cite{Kim2023} and $(0.8$--$2)\times10^{-3}$
on trapped-ion processors~\cite{Moses2023,Ransford2025} --- a $10^{4}$-gate
run suffers of order ten faults even on the best of these, so none of
these budgets is executable with physical qubits without substantial error
mitigation.  Early logical-qubit demonstrations on neutral-atom arrays have
executed circuits on tens of logical qubits with error detection and
mid-circuit readout~\cite{Bluvstein2024,Reichardt2024}, and the rodeo
algorithm has been run in two-qubit form on trapped-ion and superconducting
hardware~\cite{BeeLindgren2022,Qian2024}.  A proportionate near-term target
is a mitigated $k=3$ demonstration --- the $2\to3$ chain on seven qubits or
the $3\to4$ chain on nine, with a short refinement schedule and one or two
rodeo cycles at coarse $\delta t$, i.e.\ $10^{3}$--$10^{4}$ two-qubit gates
--- where the collapse of the acceptance under a corrupted evolution doubles
as a hardware diagnostic; the full $k=4$ budget is a workload for the first
generation of logical-qubit machines such as the neutral-atom system planned
by the QuNorth initiative~\cite{QuNorth}, and Table~\ref{tab:depth} is
offered as its specification.

""" + s[j:]
(ROOT / "paper/eigenq_pairing.tex").write_text(s)
(ROOT / "paper/eigenq_pairing_supp.tex").write_text(sm)
print("trim pass done")


# ================================================================== third pass: circuit figures move to where they are first cited; abstract trim
s = (ROOT / "paper/eigenq_pairing.tex").read_text()
s, f1 = cut_env(s, "fig:circ-vqe")
s, f2 = cut_env(s, "fig:circ-rodeo")
anchor1 = "what an imperfect input costs.\n\n\\subsection{The prolongation operator}"
assert s.count(anchor1) == 1
s = s.replace(anchor1, "what an imperfect input costs.\n\n" + f1 + "\n\\subsection{The prolongation operator}")
anchor2 = "\\paragraph{Convergence criterion and the random-time ensemble.}"
assert s.count(anchor2) == 1
s = s.replace(anchor2, f2 + "\n" + anchor2)
s = rep(s, """and are simulated
exactly within the particle-number sector against full configuration
interaction, coupled-cluster doubles, UCCSD and the analytic
success-probability formulas of the filter.""", """and are simulated
exactly within the particle-number sector against full configuration
interaction, coupled-cluster doubles and the analytic success-probability
formulas of the filter.""")
s = rep(s, """Every stage behaves as predicted: refinement raises the
overlap of the prolonged determinant from $0.87$ to $0.998$ along a gapped
path with no trained parameter, and the filter removes the remainder
geometrically,""", """Refinement raises the overlap of the prolonged determinant from $0.87$ to
$0.998$ along a gapped path with no trained parameter, and the filter
removes the remainder geometrically,""")
(ROOT / "paper/eigenq_pairing.tex").write_text(s)
a = s[s.index("\\begin{abstract}"):s.index("\\end{abstract}")]
print("abstract words:", len(re.sub(r"\\textbf\{[A-Za-z]*:\}", "", a).split()))
