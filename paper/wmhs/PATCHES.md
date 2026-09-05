# Patches for the four `\input` files

> ## STATUS — 2026-08-31, after the P0/P1 pass
>
> **The four `\input` files are now in this repository** (`sections/`), so the
> "find → replace against the compiled PDF" framing below is historical. Edit
> the files directly. `refs.bib` is here too, and the paper builds with zero
> undefined citations or references.
>
> **Landed.** A1 (1.37% → 1.3516%) · A2 (the bound collapses at both ends) ·
> A3 (cut to the claim the bound alone supports — the per-arm ρ table still does
> not exist) · A4 (dropped `66 of 190` and `3.9×`; kept the 64→16 flip) ·
> B1 (control ranges → −1.000…−0.951 and −1.000…−0.952) · B2 (63% → 60%) ·
> C1 (Appendix B written) · C3/G1 (benchmark appendix) · D1 · D2 ·
> E1/G4 (`withdrawn.tex` carries no `\section`) · F (nine bib entries) ·
> H1b (Figure 1 regenerated; the `ok = 120` annotation is gone) ·
> H2 (bench re-run from a committed tree).
>
> **Bibliography verified 2026-08-31.** All eight added entries checked against
> the publishers' own records — author order, title, venue, volume, issue, pages.
> All correct as written. cacoa was the flagged one: still a preprint, so the
> bioRxiv entry is right rather than stale.
>
> **Also landed, beyond the patch list:** the sweep behind Figure 1 is now
> committed as `src/harness/calibration_gap.py` and reproduces the previously
> published tables exactly; the provenance check counts untracked files and
> names them in the sidecar, with a regression test that fails against the old
> version; Appendix Fig. 2 has a table and a `make_fig3.py` behind it.
>
> **G2 — checked, and the premise had already changed.** G2 says cut three of
> five Appendix A items. The list is already down to three: the two "ordinary
> implementation errors" were cut before this pass. Of the three left, G2
> explicitly says *keep* items 1 and 2, and item 3 (the provenance stamp) is now
> the one defect in the paper that is closed — fixed, regression-tested, and the
> reason §3's numbers could be re-derived. Cutting three would delete everything
> G2 wanted kept. Did what the cut was *for* instead: replaced the defensive
> framing ("we record them because they circulated internally") with what the
> reader gets, and corrected the count in `withdrawn.tex` — three items are
> described, two are named and not described.
>
> **Still open.** G3 (the abstention figure is gone; nothing references it) · the
> two §4 defects remain live *by choice* and the appendix says so · fill in
> `\artifacturl` and run `check_anonymity.sh` · full paper vs 4-page abstract.
>
> ## 2026-09-05 — the third withdrawn guard
>
> **Landed.** A third guard joins §3: the premise check behind the silencing
> reading passed on detection rate while its own controls moved 1.19–1.64-fold
> in expression, because at a reference detection rate of 0.967–1.000 the
> statistic has nowhere to go. **The tolerance sat above the statistic's
> attainable range — the confound diagnostic's defect, in a guard we wrote to
> this paper's own rule, while writing this paper.** The premise then acquired a
> third verdict for the same reason the estimator has one, and returns
> *undecided* on all three cohorts. §5 and the conclusion carry it.
>
> Two counts were wrong and are fixed: the conclusion claimed "two of this
> paper's three statistics" while listing three, and `withdrawn.tex` said two
> guards. Both are now checked rather than proofread — `tests/test_paper_numbers.py`
> re-derives every figure in the new paragraph from
> `results/*/coexpression_silencing*.parquet` and asserts it against the literal
> string in the `.tex`. All eleven assertions were mutation-tested: each one
> fails when its number is perturbed. Editing prose or table alone now breaks CI.
>
> **Page budget — NOT verified.** `neurips_2026.sty` is not vendored, so the
> real page count could not be measured. Against a geometry-matched stub the
> full build grows one page (10 → 11 stub pages) and the short build's main text
> is unchanged (5 → 5). The stub runs ~1.4× longer than the real style, so the
> estimate is roughly 7 → 7.7 of 9 real pages. **Drop the official style in and
> run `./build.sh` before submitting.**


`main.tex` is edited and complete. These are the changes I could not make because
the files are not in the repo. Each is find → replace against the text as it
appears in the compiled PDF.

---

## A · `sections/withdrawn.tex` (§4)

### A1 — BLOCKING. The crossing is 1.3516%, not 1.37%

**Find:** `Below $p = 1.37\%$ the 0.20 tolerance sits above the bound`
**Replace:** `Below $p = 1.3516\%$ the 0.20 tolerance sits above the bound`

Solving `sqrt(3p(1−p)) = 0.20` gives p = 0.013516. At p = 1.37% the bound is
**0.2013**, which is *above* 0.20 — so at the figure as printed the check *can*
still fire, and the sentence asserts the opposite. Also fix
`docs/gate_memo_w2.md:1310`, which carries the same value.

### A2 — the bound collapses at *both* ends, and you tell only half

**Find:** `and got filed as cleanest of four.`
**Replace:** append —

```latex
The same collapse happens at the other end. On the coarsest resolution every
epithelial cell is called mature, prevalence is exactly $1.000$, and the bound is
exactly zero: the correlation is undefined on all 64 patient-by-axis rows and the
check reports nothing at all. \textbf{Degenerate labels are invisible to this
check at either extreme}, and the two extremes are the coarsest and the finest
resolution we ran.
```

Verified on `results/2026-09-04_3f6c07e/depth_confound_reference.parquet`:
64 epithelial rows, `worst_rho` NaN on all 64. This is a strict upgrade — it
covers a resolution the paper currently leaves out, and it makes §4's "the blind
spot sits exactly where the science is" into the sharper "the check is blind at
both ends of its own range, and the science lives at one of them."

### A3 — the normalised-ranking paragraph has no committed table behind it

`depth_confound_reference.parquet` carries only `worst_rho`, the max over arms.
There is **no per-arm ρ column anywhere in `results/`**, so `0.179`, `0.174`,
`0.175`, `0.185`, `0.092`, `0.094` cannot be checked, and the reference-arm
values that *are* committed are 0.188 / 0.174 / 0.171 — so `0.179` and `0.175`
do not match them either.

Two options, and I'd take the second unless the corrected diagnostic can be run
and committed today:

**Keep it** — run the corrected diagnostic, commit the parquet with both arms'
ρ, bound and prevalence, and quote from that.

**Cut to the qualitative claim**, which follows from the bound alone:

```latex
Normalising by the attainable bound changes the ranking. Two resolutions whose
ceilings differ by a factor of five are not comparable on the raw statistic at
all: the rarest resolution's ceiling of 0.066 and the coarser resolutions'
ceilings near 0.37 mean the same $|\rho|$ carries five times the evidence at one
resolution as at the other. Cleanest-of-four was a ranking of ceilings.
```

### A4 — `66 of 190` and `up to 3.9×` are unbacked; `64 of 64 → 16` is fine

**Find:** `Correcting the pairing moves 66 of 190 patient-by-rung rows, by up to 3.9$\times$. The reachability flag took the maximum too, and`
**Replace:** `Correcting the pairing takes the reachability flag with it, and`

The reachability flip (64 of 64 → 16) reproduces exactly and should stay. The
row-movement figures need the same missing per-arm column as A3.

---

## B · `sections/responsible.tex` (§5)

### B1 — BLOCKING. −1.11 is not a possible value

**Find:** `moves by $-0.57$ to $-0.82$ relative to its own baseline, against $-0.62$ to $-1.11$ for the compositional controls, and the ranges overlap.`
**Replace:**

```latex
moves by $-0.951$ to $-1.000$ relative to its own baseline, against $-0.952$ to
$-1.000$ for the compositional controls: the retained-gene control sits
\emph{inside} the compositional controls' range almost exactly.
```

The quantity is `m_T/m_N − 1`, bounded below by −1 for non-negative means, so
−1.11 cannot occur. Recomputed on
`results/2026-08-28_6f81018/decomposition_summary.parquet`, tier A/B/D. The
corrected numbers make your claim **stronger**, not weaker — "overlap" becomes
"inside".

While you are there: tier B's range is $-1.000$ to $+0.146$, which is the
sign change the section already asserts, now with a number.

### B2 — 63% should be 60%, and the input needs five minutes

**Find:** `still tracks depth at 63\% of the correlation its own prevalence permits`
**Replace:** `still tracks depth at 60\% of the correlation its own prevalence permits`

The committed ρ is **0.31** (`decomposition_lee_smc.meta.json`,
`docs/open_decisions.md:3196`); 0.33 appears in no table. 0.31 / 0.52 = 0.596.
One of the two sources is stale — worth finding which before you ship either
number.

---

## C · `sections/appendix.tex`

### C1 — BLOCKING. Appendix B is an empty heading. Fill it with this.

`\section{The statistic that could not see the estimator}` has no body; Figure 3
floats above it, orphaned. `main.tex` §4 now points at it as
`Appendix~\ref{app:generator}`, so it needs a label and a body. This is also the
paper's most venue-native material — a transferable warning to anyone validating
an in-silico trial simulator — so it is worth writing rather than deleting:

```latex
\section{The statistic that could not see the estimator}
\label{app:generator}

The failure in \S\ref{sec:blind}\emph{(i)} has a general form worth stating apart
from our instance, because the setup that produces it is the default one.

Let $\theta$ be the parameter a simulator is asked for, $\hat{\theta}$ the
estimator's output, and $\theta_{\text{realised}}$ the value the draw actually
took --- which differs from $\theta$ by sampling alone. A recovery curve plots
$\hat{\theta}/\theta$. Now suppose the estimator is a deterministic function of
the same sufficient statistics the generator used to construct the draw. Then
$\hat{\theta} = \theta_{\text{realised}}$ identically, not approximately, and
\[
\frac{\hat{\theta}}{\theta}
\;=\; \frac{\theta_{\text{realised}}}{\theta},
\]
in which the estimator does not appear. The curve is a diagnostic of the
generator's sampling behaviour. It will sit near 1, tighten with $n$, and look
exactly like a validated estimator, because sampling error does shrink with $n$
--- and it would do all of that if the estimator were replaced by any other
estimator with the same property.

\paragraph{Why this is easy to walk into.} Oracle arms are built precisely so the
estimator sees clean inputs, and the cleanest available input is the summary the
generator just computed. The condition is a property of the plumbing rather than
of the model, so it survives code review, does not depend on the estimator being
wrong, and leaves no trace in the curve. In our case the curve took the closed
form $1/(1-s)$ scaled by the realised-over-requested mature fraction, with
\emph{zero variance} across all 50 replicates at each shift --- a constant, read
at the time as a two-fold over-estimate in the dangerous direction.

\paragraph{The one-line check.} Compare the estimate against the \emph{realised}
truth rather than the parametric truth, and test whether the difference is
identically zero:
\begin{center}
\texttt{assert not np.allclose(est - truth\_realised, 0.0)}
\end{center}
If it is zero the recovery curve cannot see the estimator, and the arm needs a
generator whose realised summary statistics the estimator does not receive.
The statistic is not worthless --- it still catches an estimator that fails to
reproduce its input, and the design ships a broken-estimator control returning
$0.0$ for which the ratio comes out zero. It simply cannot certify the arm it was
being used to certify, and the design document never said so.
```

Figure 3 then belongs inside this section, and its existing caption needs no
change.

### C2 — BLOCKING. Appendix A claims regression tests that are not in the repo

**Find:** `The guards in §4 and the two defects in §2 now carry regression tests failing against the code as it stood`

Neither defect is fixed on this branch and the calibration test does not exist:

- `src/harness/depth_confound.py:189,199` still takes `max(rhos)` and
  `max(ceilings)` as independent maxima — the exact defect §4 describes.
- `src/harness/calibration.py:123-124` still lets `NaN` comparisons coerce to
  `False`, so abstentions still score as failures.
- There is no `tests/test_calibration.py`. The assertion message the appendix
  quotes verbatim — *"coverage moved by 0.50× when only undefined intervals were
  added"* — appears nowhere in `tests/`.

The bound itself *is* tested (`tests/test_depth_confound.py:269`, four
prevalences). Either land the two fixes with their tests, or change the sentence
to the future tense. Quoting an assertion message that does not exist is the
single most checkable claim in the paper.

### C3 — new appendix for the benchmark referenced from §3

`main.tex` now cites `Appendix~\ref{app:bench}`. Add:

```latex
\section{What the alternatives do where the estimand does not exist}
\label{app:bench}

Six worlds, 200 replicates, 2{,}000 cells per arm, seed 20260829. Mature cells
express $\mathrm{Poisson}(\mu)$ in reference tissue and $\mathrm{Poisson}(\mu s)$
in diseased tissue; immature cells express nothing. The arm mean is
$f \times$ (mean among mature), which is exactly what
Equation~\ref{eq:kitagawa} splits, so the truth is available in closed form
rather than by simulation. The world that matters sets the diseased mature
fraction to zero, so no diseased mature cell survives and the intrinsic estimand
has no referent.

\begin{center}
\small
\begin{tabular}{@{}lccc@{}}
\toprule
\textbf{Method} & \textbf{returned a number} & \textbf{false confidence} &
\textbf{median $|$intrinsic$|$ invented} \\
\midrule
gated (ours)        & $0/200$   & $0.00$ & --- \\
composition-only    & $0/200$   & $0.00$ & --- (no intrinsic arm) \\
ungated ablation    & $200/200$ & $1.00$ & 7.99 \\
pseudobulk DE       & $200/200$ & $1.00$ & 7.99 \\
naive $\Delta$ mean & $200/200$ & $1.00$ & 20.00 \\
\bottomrule
\end{tabular}
\end{center}

\texttt{composition-only} scores $0.00$ for a different reason: it has no
intrinsic arm, which is inapplicability rather than caution, and crediting it as
a refusal would reward a method for not competing. The ungated ablation and
pseudobulk DE coincide here because at a diseased mature fraction of zero both
reduce algebraically to $-(f_N m_N)$; they diverge as soon as any diseased mature
cell exists.

Scored where the estimand \emph{does} exist and a real intrinsic effect is
present, the gate declines 88 of 600 available answers, all in the
$\sim$20-cell regime at the cutpoint, and median absolute error is
\emph{lower} with the gate than without ($0.111$ against $0.122$) because the
dropped cases are the noisy low-count ones.

\paragraph{Limitations.} Everything here is synthetic and the generative model is
exactly the one Equation~\ref{eq:kitagawa} assumes, which is a favourable
setting. It tests the gate, not robustness to a misspecified model, and says
nothing about real tissue. The competitors are faithful reimplementations rather
than the published software. The cutpoint is our own pre-committed rule, so the
detection rate would move if it moved. Only the intrinsic term's estimability is
tested.
```

**Before this appendix can be quoted**, two fixes in `submission/`:

1. `submission/results/bench.meta.json` records `git_branch: w1/threshold-sweep`,
   `git_sha: 730bec0`, `git_dirty: true`. HEAD is `0ec0891` on
   `submission/competitor-bench`. The recorded sha does not identify the code
   that ran. One clean re-run on the branch fixes it; the tables come back
   bit-identical.
2. `submission/FINDINGS.md` reports signed errors `+0.017` and `+0.007`;
   `sensitivity_where_estimable.parquet` says `−0.005054` and `−0.002227` —
   opposite signs. Absolute errors match to three decimals, so it is the signed
   column alone, carried over from an earlier run.

---

## D · `sections/llm.tex`

### D1 — BLOCKING. Remove the draft marker

**Find:** `\textbf{DRAFT. Requires author sign-off before submission.}` → delete.

### D2 — BLOCKING. The provenance sentence is not true of Figures 1 and 3

**Find:** `Every number here gets read from a versioned result table by the figure code rather than transcribed, and each got independently re-derived before reaching the text.`

I searched all 33 directories under `results/`: **there is no calibration or
cutpoint-sweep parquet anywhere.** Figure 1 is the paper's central figure and
Figure 3 is its key appendix figure, and neither has a versioned table. The
`(0.36, 0.31)` row is quoted from an internal report, not a table.

This matters more here than it would in another paper, because Appendix A item 5
*tells the reviewer* that this exact calibration was one of three results
produced by an uncommitted script, stamped clean against a commit that did not
contain it. You have handed a reproducibility-minded reviewer the map and then
claimed the property it lacks.

**Either** run the sweep under the project's provenance rule and commit the
parquet — then the sentence stands as written — **or** replace with:

```latex
Except where a number is attributed to an internal report, the figure code reads
every number in this paper from a result table carrying a commit hash and a fixed
seed. The calibration sweep behind Figures~\ref{fig:calibration}
and~\ref{fig:generator} is re-derived from committed code at a fixed seed but its
output table is not yet versioned under that rule, for the reason
Appendix~\ref{app:ours} gives. Every number was independently re-derived before
reaching the text.
```

Shipping the original sentence *and* Appendix A item 5 is the one combination
that does not work.

---

## E · Structural changes in `main.tex` that touch your section files

### E1 — `withdrawn.tex` lost its own `\section{}` header

It used to sit under `\section{Two guards that could not do their job}`. It now
sits inside **§4, "Why our validation statistics could not see any of this"**,
after the three-statistics list and the general-form paragraph. Those two guards
*are* checks that could not fire, so they belong under that heading rather than
in a section of their own, and the merge buys back the space the new §2 material
costs.

**Check `withdrawn.tex` does not open with its own `\section` or
`\subsection`.** If it opens with `\paragraph{}`s it drops straight in. If it
opens with a section header, either delete that line or demote it to
`\paragraph{Two guards that could not do their job.}`.

### E2 — three labels must resolve

`main.tex` now references:

| label | lives in | status |
|---|---|---|
| `app:bench` | new appendix, patch **C3** | you must paste it |
| `app:generator` | Appendix B, patch **C1** | you must paste it |
| `app:limitations` | your existing Appendix C | should already resolve |

The `(Appendix ??)` you saw was `app:bench` with no target. Both new appendices
are written out above; paste them and the references resolve.

### E3 — the page budget

§2 grew by roughly half a page (the named methods, the quantified error, the
non-circularity paragraph). Paid for by merging the old §4 header into §4, by
compressing the calibration setup, and by moving the modal-outcome result into
§2 as a paragraph rather than its own section. If it still overruns, the honest
cut is **Figure~\ref{fig:abstention}** — its two headline numbers (100\%/86\%,
and the six *not identifiable* verdicts) survive in the prose, and it is the only
figure whose deletion costs no argument. Cut the figure before cutting the
non-circularity paragraph.

---
---

# ROUND 2 — supersedes parts of the above

`main.tex` was restructured after review. A1, A2, A3, A4, B1, B2, C1, D1, D2 are
**unchanged and still needed**. C3 is **superseded by G1** below.

## F · `refs.bib` — nine new entries

`main.tex` now cites related work and conformal prediction. **Verify every field
before submitting** — I am confident in the authors and the claims attributed to
each, less so in some volume/page numbers.

```bibtex
@article{crowell2020, title={muscat detects subpopulation-specific state
  transitions from multi-sample multi-condition single-cell transcriptomics
  data}, author={Crowell, Helena L and Soneson, Charlotte and Germain,
  Pierre-Luc and Calini, Daniela and Collin, Ludovic and Raposo, Catarina and
  Malhotra, Dheeraj and Robinson, Mark D}, journal={Nature Communications},
  volume={11}, number={1}, pages={6077}, year={2020}}

@article{buettner2021, title={scCODA is a Bayesian model for compositional
  single-cell data analysis}, author={B{\"u}ttner, Maren and Ostner, Johannes
  and M{\"u}ller, Christian L and Theis, Fabian J and Schubert, Benjamin},
  journal={Nature Communications}, volume={12}, number={1}, pages={6876},
  year={2021}}

@article{dann2022, title={Differential abundance testing on single-cell data
  using k-nearest neighbor graphs}, author={Dann, Emma and Henderson, Neil C and
  Teichmann, Sarah A and Morgan, Michael D and Marioni, John C}, journal={Nature
  Biotechnology}, volume={40}, number={2}, pages={245--253}, year={2022}}

@article{phipson2022, title={propeller: testing for differences in cell type
  proportions in single cell data}, author={Phipson, Belinda and Sim, Choon Boon
  and Porrello, Enzo R and Hewitt, Alex W and Powell, Joseph and Oshlack,
  Alicia}, journal={Bioinformatics}, volume={38}, number={20},
  pages={4720--4726}, year={2022}}

@article{zhao2021, title={Detection of differentially abundant cell
  subpopulations in scRNA-seq data}, author={Zhao, Jun and Jaffe, Ariel and
  Li, Henry and Lindenbaum, Ofir and Sefik, Esen and Jackson, Ruaidhr{\'\i} and
  Cheng, Xiuyuan and Flavell, Richard A and Kluger, Yuval}, journal={PNAS},
  volume={118}, number={22}, year={2021}}

% VERIFY THE VENUE. cacoa was on bioRxiv for a long time before journal
% publication; the repo pins it from GitHub (kharchenkolab/cacoa).
@article{petukhov2022, title={Case-control analysis of single-cell RNA-seq
  studies}, author={Petukhov, Viktor and van den Brand, Teun and Biederstedt,
  Evan and others}, journal={bioRxiv}, year={2022}}

@book{vovk2005, title={Algorithmic Learning in a Random World}, author={Vovk,
  Vladimir and Gammerman, Alexander and Shafer, Glenn}, publisher={Springer},
  year={2005}}

@article{angelopoulos2023, title={A gentle introduction to conformal prediction
  and distribution-free uncertainty quantification}, author={Angelopoulos,
  Anastasios N and Bates, Stephen}, journal={Foundations and Trends in Machine
  Learning}, year={2023}}
```

## G · appendix and section changes

### G1 — `app:bench` appendix, SUPERSEDING patch C3

The benchmark now carries a sixth method and a third result. Use this instead:

```latex
\section{The competitor benchmark}
\label{app:bench}

Six worlds, 200 replicates, 2{,}000 cells per arm, seed 20260829. Mature cells
express $\mathrm{Poisson}(\mu)$ in reference tissue and $\mathrm{Poisson}(\mu s)$
in diseased tissue; immature cells express nothing. The arm mean is
$f \times$ (mean among mature), exactly what Equation~\ref{eq:kitagawa} splits,
so the truth is available in closed form rather than by simulation. The world
that matters sets the diseased mature fraction to zero: no diseased mature cell
survives and the intrinsic estimand has no referent.

Six methods, of which five run. \texttt{cacoa} ships as an adapter that reports
its own unavailability rather than being silently absent --- it installs from
GitHub and is not in the pinned environment --- so it appears as a named skip
with a reason.

\paragraph{The width gate.} Abstain when the 95\% half-width on the intrinsic
term exceeds the detectable effect, both taken from the same design document as
the count cutpoint, so neither gate is tuned against the other here. The
per-cell moments are recovered from the full arm without needing the mature-cell
mask, because immature cells express exactly zero and contribute nothing to
either the sum or the sum of squares. \textbf{Finiteness is checked before the
comparison.} Written the obvious way --- comparing first --- the gate returns
\texttt{False} on an undefined standard error and declines to abstain in exactly
the world that motivates it.

\paragraph{Limitations.} Everything here is synthetic, and the generative model
is exactly the one Equation~\ref{eq:kitagawa} assumes, which is a favourable
setting. It tests the gate, not robustness to a misspecified model. The
competitors are faithful reimplementations, not the published software. The
naive difference of arm means and pseudobulk DE answer a slightly different
question, so their bias is partly unfair to them; the false-confidence result
does not depend on them and the ablation does not depend on them at all.
Result~3 is a statement about these regimes: the width criterion fails to bind
partly because a halving of per-cell output is large relative to sampling error
at 2{,}000 cells per arm, and a smaller detectable effect or a heavier-tailed
expression model would move the binding point. What survives regardless is that
the finiteness guard, not the threshold, catches the undefined estimand.
```

### G2 — Appendix A: cut three of the five

The reviewer is right that forty lines of your own implementation errors reads as
less careful, not more. **Keep items 3 and 4** --- a population definition
selected by file order, and a join that validated its keys and never its counts.
Both have conceptual content and both are failures a careful group would repeat.
**Cut items 1, 2 and 5**; they are ordinary bugs.

One interaction to note before you cut item 5: it is the entry that explains
*why* the calibration table is not committed under the project's provenance rule,
which patch **D2** leans on. If item 5 goes, D2's replacement wording should
simply state the table is not yet versioned, without the back-reference to
Appendix~\ref{app:ours}.

### G3 — the abstention-rate figure is no longer referenced

`main.tex` folds that result into prose at the end of \S\ref{sec:calibration}.
Both headline numbers survive (28/28 and 24/28; the six *not identifiable*
verdicts). **Delete `figures/fig2_abstention.pdf` and its float** --- that space
pays for the new Table~\ref{tab:bench} and the related-work and conformal
paragraphs. If you would rather keep the figure, the paragraph is the thing to
cut, not the table.

### G4 — `withdrawn.tex` now sits under \S5 with no header of its own

Unchanged from patch E1: check it does not open with its own `\section`.

## H · what I could NOT do, and it matters

### H1 — RESOLVED. The gap is closed, and it changed the answer

The four Lee/SMC files are downloaded to `data/raw/lee/` and **sha256-verified
against `data/manifest.csv`** (all four exact). The sweep runs in seconds, not
hours. Results are committed to
`results/2026-08-31_0ec0891/calibration_gap_{cutpoints,bins}.parquet`.

13 seeds x 2 draw pools x 2 grids, everything else held fixed:

| pool | grid | returned a cutpoint | `ok` range |
|---|---|---|---|
| pooled | committed | 0 / 13 | — (max discrimination 0.78 vs 0.80 target) |
| pooled | **extended** | **0 / 13** | — |
| reference | committed | 13 / 13 | **100 – 200** |
| reference | **extended** | 13 / 13 | **45 – 90** |

The committed grid varies a mature *fraction* against 2,000 fixed cells, so the
counts it can reach are {0, 20, 40, 100, 200, 400, 800} — nothing between 40 and
100. Adding nine fractions inside that interval (reaching 50, 60, 70, 80, 90,
110, 120, 130, 150) and changing nothing else moves the answer to 45–90. **The
two ranges do not overlap, and the new one contains the committed 50.**

So the paper can no longer say the procedure returns 120. \S4 now says the
narrower true thing: the pooled draw fails outright, and the reference draw does
not identify a cutpoint precisely enough at 50 replicates to confirm or reject
50. The crossing wanders with the seed because coverage and discrimination are
nearly flat across the interval. That is the same lesson as the recovery curve,
in a second statistic — which is why it strengthens the paper rather than
denting it.

### H1b — BLOCKING: Figure 1 now contradicts its own caption

`figures/fig1_calibration.pdf` has `calibrated ok = 120` annotated on the right
panel with a dotted vertical line at 120. **That annotation is now wrong** and
the text says so three paragraphs later. I do not have the figure code, so this
one is yours. Either regenerate the right panel from
`results/2026-08-31_0ec0891/calibration_gap_bins.parquet` — filtering
`pool == "reference"`, and showing the 45–90 band across seeds rather than a
single dotted line — or drop the annotation and the line entirely and let the
caption carry it. Do not ship the figure as it stands.

### H2 — the benchmark's provenance stamp is still dirty

`bench.meta.json` now records the right branch and sha (`0ec0891`,
`submission/competitor-bench`) but `git_dirty: true`, because the width-gate
change is uncommitted. Commit `submission/competitors.py`,
`submission/FINDINGS.md` and `tests/test_submission_bench.py`, then re-run
`python -m submission.run_bench`. The tables are deterministic, so only the
stamp changes.
