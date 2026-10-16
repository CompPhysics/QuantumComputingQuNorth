#!/usr/bin/env python3
"""One-off (2026-10-16): write the ADAPT-VQE comparator (notebook Sec. 21,
fig 23) into paper I: Sec. V.D rewritten with Table tab:adapt and Fig.
fig:adapt, rows in Table tab:largebudget, and the verdicts in Secs. V.C, V.F,
the abstract, the introduction and the discussion."""
import pathlib, re
ROOT = pathlib.Path(__file__).resolve().parents[1]
p = ROOT / "paper/eigenq_pairing.tex"; s = p.read_text()

def rep(a, b, n=1):
    global s
    assert s.count(a) == n, (a[:70], s.count(a)); s = s.replace(a, b)

# ---------------- Sec. V.D: the adaptive comparator
i = s.index("\\subsection{Remarks on adaptive and improved VQE}"); j = s.index("\\subsection{A larger model space: four pairs on eight levels}")
s = s[:i] + r"""\subsection{The adaptive variational comparator}
\label{sec:adapt}
The steep UCCSD gate count is not irreducible, and much recent work aims
precisely at shallower variational circuits (see Ref.~\cite{Tilly2022} for a
review).  ADAPT-VQE~\cite{Grimsley2019} grows the ansatz one operator at a
time, at each step appending the pool operator with the largest energy
gradient $\abs{\expval{[H,K]}{\psi}}$ --- an operator may be selected again
--- and re-optimizing all amplitudes.  The pairing model makes the pool
transparent: because the Hamiltonian conserves seniority and the reference
is seniority-zero, the gradients of all singles and of all pair-breaking
doubles vanish identically on the seniority-zero manifold, so the pool is
either the $\frac N2(k-\frac N2)$ occupied-pair$\to$virtual-pair rotations
of the restricted ansatz of Sec.~\ref{sec:depth} or the \emph{generalized}
pool of all $k(k-1)/2$ pair-hopping generators $P_q^\dagger P_p-{\rm h.c.}$,
each a $48$-CNOT block in the $2k$-qubit encoding.  We ran ADAPT-VQE with
both pools, in the seniority-zero pair basis of Sec.~\ref{sec:large} and
with the exact gradient, for $(k,N)=(4,4)$ and $(6,4)$ at $g=1$ and for
$(8,8)$ at $g=1$ and $2$ (Table~\ref{tab:adapt}, Fig.~\ref{fig:adapt}).
At $k=4$ the restricted pool reproduces the four-pair ansatz exactly ---
four operators, $192$ CNOTs, $1.78\times10^{-4}$ --- and the generalized
pool adds the occupied--occupied pair rotation $(0,1)$ as a fifth operator
and reaches the exact ground state ($1-F<10^{-15}$, $240$ CNOTs): the
six-dimensional seniority-zero space is spanned.  At sixteen qubits the
adaptive ansatz does what the single-layer one cannot: with the
generalized pool $19$ operators ($912$ CNOTs) reach $1.5\times10^{-4}$,
$23$ ($1\,104$ CNOTs) reach $1.7\times10^{-5}$ and $27$ ($1\,296$ CNOTs)
reach $4\times10^{-6}$ at $g=1$, and $25$, $29$ and $31$ operators
($1\,200$--$1\,488$ CNOTs) reach $10^{-3}$, $3\times10^{-4}$ and $10^{-4}$
at $g=2$, where the single-layer ansatz stops at $3\times10^{-3}$ and
$4\times10^{-2}$; the restricted pool with repeated operators is slower
($29$ operators for $10^{-4}$ at $g=1$, no $3\times10^{-4}$ within $40$ at
$g=2$).  Fed to the Trotterized filter of Sec.~\ref{sec:large}, the
$1\,104$-CNOT state at $g=1$ starts \emph{below} the floor of the digital
filter ($1.7\times10^{-5}$ against $4.6\times10^{-5}$), so the filter can
only degrade it, and the $1\,488$-CNOT state at $g=2$ gains a factor two in
six cycles at four hundred times the cost.  Qubit-ADAPT-VQE~\cite{Tang2021}
would push further by drawing the pool from individual Pauli strings.  The
usual qualifications apply: the amplitudes are trained noise-free with
exact gradients, each iteration measures $\abs{\rm pool}$ gradients and
re-optimizes the ansatz, and these one-time costs are excluded from the
gate budget as they are for UCCSD; and the pool is the symmetry-adapted one
that the seniority conservation of the pairing force makes available.  A
different route to a shallower circuit --- keeping the Trotterized
refinement path as the ansatz structure and optimizing a few of its layers,
as in the optimized nucleation of Ref.~\cite{Sadia2026} --- is taken up in
Sec.~\ref{sec:discussion}.

\begin{table}[t]
\caption{\label{tab:adapt}ADAPT-VQE in the seniority-zero pair basis
(gradient tolerance $10^{-6}$, at most $40$ operators; $48$ CNOTs per
operator in the $2k$-qubit encoding, $4$ in the pair encoding): the pool
size and the number of operators (CNOTs) at which the infidelity first
falls below $10^{-3}$, $3\times10^{-4}$, $10^{-4}$ and $10^{-5}$; a dash
means not within $40$ operators.  The last column gives the infidelity of
the single-layer restricted ansatz used in Tables~\ref{tab:depth} and
\ref{tab:largebudget}.}
\footnotesize
\begin{ruledtabular}
\begin{tabular}{ccc cccc c}
$(k,N,g)$ & pool & size & $10^{-3}$ & $3\times10^{-4}$ & $10^{-4}$ & $10^{-5}$ & single layer\\
\colrule
$(4,4,1)$ & pairs & 4 & 4 (192) & 4 (192) & --- & --- & $1.8\times10^{-4}$\\
$(4,4,1)$ & generalized & 6 & 4 (192) & 4 (192) & 5 (240) & 5 (240) & \\
$(6,4,1)$ & pairs & 8 & 8 (384) & 9 (432) & 11 (528) & 13 (624) & $7.6\times10^{-4}$\\
$(6,4,1)$ & generalized & 15 & 8 (384) & 9 (432) & 9 (432) & 9 (432) & \\
$(8,8,1)$ & pairs & 16 & 17 (816) & 28 (1\,344) & 29 (1\,392) & --- & $3.1\times10^{-3}$\\
$(8,8,1)$ & generalized & 28 & 19 (912) & 19 (912) & 23 (1\,104) & 27 (1\,296) & \\
$(8,8,2)$ & pairs & 16 & 37 (1\,776) & --- & --- & --- & $4.2\times10^{-2}$\\
$(8,8,2)$ & generalized & 28 & 25 (1\,200) & 29 (1\,392) & 31 (1\,488) & --- & \\
\end{tabular}
\end{ruledtabular}
\end{table}

\begin{figure*}[t]
\centering
\includegraphics[width=\linewidth]{figs/fig23.png}
\caption{ADAPT-VQE in the seniority-zero pair basis.  Left: infidelity
versus the number of operators for the restricted (dashed) and the
generalized (solid) pool at $(k,N)=(4,4)$, $(6,4)$ and $(8,8)$.  Right:
infidelity against CNOTs for the generalized pool; the dotted curves feed
the $(8,8)$ ADAPT states at the first iteration below $10^{-4}$ to the
Trotterized filter of Sec.~\ref{sec:large} ($M=0$ to $6$), which cannot
improve on an input that already lies below its floor.}
\label{fig:adapt}
\end{figure*}

""" + s[j:]

# ---------------- Sec. V.C: the k=4 verdict mentions the exact five-operator state
rep(r"""unfiltered four-pair ansatz, $192$ CNOTs at $1.8\times10^{-4}$, followed by
the same state in its $1\,312$-CNOT template; the cheapest route through the""",
r"""unfiltered four-pair ansatz, $192$ CNOTs at $1.8\times10^{-4}$, followed by
the same state in its $1\,312$-CNOT template (and the adaptive ansatz of
Sec.~\ref{sec:adapt} reaches the exact state with one more operator, $240$
CNOTs); the cheapest route through the""")

# ---------------- Sec. V.F: table rows and verdict
rep(r"""\multicolumn{9}{l}{$g=2$ ($E_0=2.48659$, cycle $9.0\times10^{4}$ CNOTs)}\\""",
r"""ADAPT, generalized pool ($23$ op.) & $0.99998$ & $1\,104$ & $1.7\times10^{-5}$ & $4.1\times10^{-5}$ & $4.8\times10^{-5}$ & $1.1\times10^{3}$ (0) & $1.1\times10^{3}$ (0) & $1.1\times10^{3}$ (0)\\
\colrule
\multicolumn{9}{l}{$g=2$ ($E_0=2.48659$, cycle $9.0\times10^{4}$ CNOTs)}\\""")
rep(r"""ref.\ $4{\to}8$, $\mu_{\rm buf}{=}8$, $T{=}8$ ($n_s{=}80$) & $0.9992$ & $258\,240$ & $7.6\times10^{-4}$ & $1.8\times10^{-4}$ & $3.5\times10^{-5}$ & $2.6\times10^{5}$ (0) & $4.4\times10^{5}$ (2) & $6.4\times10^{5}$ (4)\\
\end{tabular}""",
r"""ref.\ $4{\to}8$, $\mu_{\rm buf}{=}8$, $T{=}8$ ($n_s{=}80$) & $0.9992$ & $258\,240$ & $7.6\times10^{-4}$ & $1.8\times10^{-4}$ & $3.5\times10^{-5}$ & $2.6\times10^{5}$ (0) & $4.4\times10^{5}$ (2) & $6.4\times10^{5}$ (4)\\
ADAPT, generalized pool ($31$ op.) & $0.99991$ & $1\,488$ & $8.9\times10^{-5}$ & $4.2\times10^{-5}$ & $2.8\times10^{-5}$ & $1.5\times10^{3}$ (0) & $1.5\times10^{3}$ (0) & $1.5\times10^{3}$ (0)\\
\end{tabular}""")
rep(r"""comparison of Sec.~\ref{sec:depth} with ten inputs at $g=1$ and six at $g=2$
(forty common time sets, $M=0$ to $6$, bootstrap standard errors below
$10\%$ for every input with $p>0.99$).  Three things have changed relative
to $k=4$.  First, the unfiltered variational preparation no longer meets
any target: the $16$-pair ansatz needs two, four and six cycles for""",
r"""comparison of Sec.~\ref{sec:depth} with eleven inputs at $g=1$ and seven at
$g=2$ (forty common time sets, $M=0$ to $6$, bootstrap standard errors
below $10\%$ for every input with $p>0.99$).  Three things have changed
relative to $k=4$ among the fixed-ansatz and adiabatic routes; the adaptive
ansatz, discussed last, stands apart.  First, the unfiltered single-layer
variational preparation no longer meets any target: the $16$-pair ansatz
needs two, four and six cycles for""")
rep(r"""Refinement followed by the filter is now \emph{cheaper}
than the variational route at every target --- $1.3\times10^{5}$,""",
r"""Refinement followed by the filter is now \emph{cheaper}
than this variational route at every target --- $1.3\times10^{5}$,""")
rep(r"""The answer to the question that motivates this subsection is therefore
split.  A larger and more correlated space does favor refinement followed
by the filter over the variational preparation, and over filtering the bare
determinant, because the single-reference ansatz degrades and the cycles
become expensive while the adiabatic transport, properly shifted, does
not slow down.  It does not favor the shift-controlled projector path over
the ordinary adiabatic ramp from the same determinant, which is cheaper by
the coarse-term overhead at every size and coupling studied; the specific""",
r"""The adaptive variational ansatz of Sec.~\ref{sec:adapt} changes the scale
of the whole table: $23$ operators ($1\,104$ CNOTs) reach $1.7\times10^{-5}$
at $g=1$ and $31$ ($1\,488$ CNOTs) reach $8.9\times10^{-5}$ at $g=2$,
without any filter cycle, two to three orders of magnitude below every
other route at every target, and below the floor of the digital filter, so
that the filter has nothing to add to them.  Among the routes that need no
classical optimization the ramp is the cheapest; among all routes the
adaptive ansatz is, at both sizes and both couplings.

The answer to the question that motivates this subsection is therefore
split three ways.  A larger and more correlated space does favor
refinement followed by the filter over a \emph{fixed} single-layer
variational preparation, and over filtering the bare determinant, because
the single-reference ansatz degrades and the cycles become expensive while
the adiabatic transport, properly shifted, does not slow down.  It does not
favor it over an \emph{adaptively} grown ansatz, which reaches the targets
with of order a thousand CNOTs because the seniority-zero pool is small and
exactly the right one.  And it does not favor the shift-controlled
projector path over the ordinary adiabatic ramp from the same determinant,
which is cheaper by the coarse-term overhead at every size and coupling
studied; the specific""")

# ---------------- abstract
rep(r"""at sixteen qubits the ansatz degrades to $3\times10^{-3}$, filter cycles
become eleven times dearer, and refinement plus filter is $1.5$--$2$ times
cheaper than the variational route while staying within $4$--$19\%$ of the
ramp.""", r"""at sixteen qubits the single-layer ansatz degrades to $3\times10^{-3}$,
filter cycles become eleven times dearer, and refinement plus filter is
$1.5$--$2$ times cheaper than that route while staying within $4$--$19\%$
of the ramp --- but an ADAPT-VQE ansatz grown from the seniority-zero pool
reaches $10^{-4}$ with about a thousand CNOTs at both sizes, below the
floor of the digital filter.""")
rep(r"""On the pairing model it is never the
cheapest route --- the adiabatic ramp from the same determinant is, by the
cost of the embedded coarse term --- but it overtakes the variational
preparation as the space grows.""", r"""On the pairing model it is never the
cheapest route: among the optimizer-free routes the adiabatic ramp from the
same determinant is, by the cost of the embedded coarse term, and the
adaptive variational ansatz undercuts every filtered route by two to three
orders of magnitude at both sizes.""")

# ---------------- introduction outcome
rep(r"""At sixteen qubits
(Sec.~\ref{sec:large}) the variational input degrades and the filter cycles
become expensive, so refinement plus filter is then $1.5$--$2$ times cheaper
than the variational route --- but the plain ramp from the same determinant
stays cheaper by the cost of the embedded coarse term.""", r"""At sixteen qubits
(Sec.~\ref{sec:large}) the single-layer variational input degrades and the
filter cycles become expensive, so refinement plus filter is then
$1.5$--$2$ times cheaper than that route --- but the plain ramp from the
same determinant stays cheaper by the cost of the embedded coarse term, and
an adaptively grown ansatz (Sec.~\ref{sec:adapt}) reaches the same accuracy
with a thousand CNOTs and no filter.""")
rep(r"""the UCCSD template and the restricted pair-rotation ansatz that reproduces
it --- with the ground-state infidelity as the common observable,""",
r"""the UCCSD template, the restricted pair-rotation ansatz that reproduces it,
and an ADAPT-VQE ansatz grown from the seniority-zero pool --- with the
ground-state infidelity as the common observable,""")

# ---------------- discussion
rep(r"""At four pairs on eight
levels (Sec.~\ref{sec:large}) the variational ansatz is no longer good, the
cycles are eleven times dearer, and refinement plus filter overtakes both
the variational route and the filtered determinant by factors $1.5$--$2$
and more; but the ordinary ramp from the same determinant stays ahead by
the cost of the conditioned coarse term, at $g=1$ and at $g=2$, once the
shift is scaled with the size.""", r"""At four pairs on eight
levels (Sec.~\ref{sec:large}) the single-layer ansatz is no longer good,
the cycles are eleven times dearer, and refinement plus filter overtakes
both that route and the filtered determinant by factors $1.5$--$2$ and
more; but the ordinary ramp from the same determinant stays ahead by the
cost of the conditioned coarse term, at $g=1$ and at $g=2$, once the shift
is scaled with the size, and an ADAPT-VQE ansatz grown from the
seniority-zero pool (Sec.~\ref{sec:adapt}) reaches $10^{-4}$ with
$1\,100$--$1\,500$ CNOTs at both couplings, two to three orders of
magnitude below every filtered route and below the floor of the digital
filter.""")
rep(r"""Whether that trade pays off against a variational preparation
is, as Sec.~\ref{sec:large} shows, a question of size and coupling, and the
answer turns in favor of the adiabatic routes already at sixteen qubits;
whether""", r"""Whether that trade pays off against a fixed variational preparation
is, as Sec.~\ref{sec:large} shows, a question of size and coupling, and the
answer turns in favor of the adiabatic routes already at sixteen qubits;
against an adaptive one it does not pay off at any size studied, as long as
a small, symmetry-adapted pool exists and its gradients can be measured;
whether""")
rep(r"""And the variational comparator
is a restricted ansatz chosen by the symmetry selection rule, not an
adaptive construction; an ADAPT-VQE run proper is the next comparator to
add.""", r"""And the adaptive comparator is trained with exact gradients on a
symmetry-adapted pool that the pairing force makes available; a pool that
must include seniority-breaking operators, or gradients estimated from
shots, would raise its one-time cost and may change its convergence.""")
rep(r"""First, comparators and compilation: an ADAPT-VQE run proper, a
pair-encoded (one qubit per level) variational and filtering baseline and a
symmetry-projected BCS-like input,""", r"""First, comparators and compilation: a qubit-ADAPT run with a Pauli-string
pool, a pair-encoded (one qubit per level) variational and filtering
baseline and a symmetry-projected BCS-like input,""")
rep(r"""cost that, for the pairing model, exceeds that of an adiabatic ramp from
the same determinant at every size studied and that of a cheap variational
input at eight but not at sixteen qubits.""", r"""cost that, for the pairing model, exceeds that of an adiabatic ramp from
the same determinant at every size studied, that of a fixed variational
input at eight but not at sixteen qubits, and that of an adaptively grown
variational input by two to three orders of magnitude at both.""")
p.write_text(s)
a = s[s.index("\\begin{abstract}"):s.index("\\end{abstract}")]
print("abstract words:", len(re.sub(r"\\textbf\{[A-Za-z]*:\}", "", a).split()))
