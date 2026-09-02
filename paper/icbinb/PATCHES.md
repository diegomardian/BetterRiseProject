# ICBINB submission — defect list and fixes

> Written 2026-09-01, one day before the deadline, from an independent review
> pass: the paper read twice, built against a real NeurIPS style file, and every
> headline number re-derived against the committed parquet tables.
>
> **What holds.** The build is clean (0 errors, 0 undefined, 0 overfull). The
> main text ends at the bottom of page 8 and References opens page 9, so the
> limit is met with no margin. The √(3p(1−p)) bound is correct analytically and
> numerically, including the finite-*n* case. Fifteen of the paper's headline
> numbers reproduce **exactly** off the committed tables. That is a better
> verification rate than most submitted papers would survive.
>
> **What follows is what does not hold**, in priority order, each with its exact
> location and a paste-ready fix. Word costs are stated because the paper has
> under 20 words of slack; **P1-1 frees ~18 words and pays for the rest.**

---

## STATUS — 2026-09-01, all of the above applied

**The paper builds clean and the main text still occupies 8 pages**, References
at the top of page 9. `check_anonymity.sh` passes. `pytest` is at **1200 passed**
(up from 1187) with the same 22 pre-existing failures — `diptest`, `lifelines`
and `anndata` missing from the local env, unrelated and present beforehand.
`ruff check` clean.

**Landed.** P0-1 (typo) · P0-2 (rung names, fixed in `make_figs.py` at zero word
cost, plus `finest`→`rarest` and `rung`→`resolution`) · P1-1 (unbacked §3
numbers dropped, in main text and Appendix D) · P1-2 (tier-B abundances) ·
P1-3 (Appendix C's third escape) · P1-4 (conformal cited) · P2-1 (both modules
repaired, with failing-input tests) · P2-4 (`_tables.py` sorts by timestamp) ·
P3-1 (empty figure row) · P3-3 (§5 tightened) · **P4-1 (the rare-label
experiment, run on both cohorts and written to `results/`)** · S4 (four uncited
bib entries removed) · the framing change (§1 now names the assistant) · the
abstract now carries the demonstration rather than only the proof.

**Two sentences were cut because the new evidence made them false**, not to save
space: *"what we cannot claim is that we caught it biting"* and *"that is the
theorem talking, not this cohort"*. Both were true when written.

### The page check caught me the way it catches the paper

Adding P4-1 pushed ~81 words onto page 9. A checker that located the page
containing "References" still reported **8 pages** — because the heading had
moved *down* page 9, not onto page 10. `build.sh` measures how far into the page
the heading sits and would have caught it; my scratch checker did not. That is
this document's own subject, committed while writing it. Fixed, and the README
now says to measure the spill rather than the heading.

After the three cuts the budget is not zero but roughly **80–99 words** — the
cuts overshot the overflow. That slack paid for the framing change and the
abstract sentence, and about 40 words remain.

### Two recommendations in this document were wrong

**A5 was wrong.** I wrote that Table 2 "resolves to a markdown file" because
`results/` holds no `*bench*` table. It does not: `submission/run_bench.py`
writes `refusal_table.parquet` and `sensitivity_where_estimable.parquet` with a
`bench.meta.json` carrying sha `54686bc`, `git_dirty: false`, seed 20260829 and
200 replicates, and every Table 2 number reproduces off them. The separation
from `results/` is deliberate and documented in the runner's docstring — that
directory is governed by the frozen schema and holds results about colorectal
cancer, and these are simulated. **No change made.**

**S1 is wrong and is withdrawn.** I recommended dropping the two degenerate
annotation resolutions, arguing it would remove "the weak 58-degenerate-rows
argument". Those 58 rows are load-bearing for §3's *other* finding — that a check
with no abstention state reports absence of evidence as evidence of absence.
Dropping the resolutions would delete the evidence for one of the two failures
the paper reports hitting twice. **Not done, and it should not be.** The empty
figure row it was partly aimed at is fixed directly instead (P3-1).

### Still open, and why

- **GSE178341 is not on disk** (1.2 GB `.h5`, in the manifest). This is why P1-1
  was a deletion rather than a re-derivation. Fetching it is the prerequisite for
  ever restoring those four numbers.
- **The new table was written from a dirty tree**, its own job untracked — the
  sidecar names it, which is the provenance check working as designed. Re-run
  after committing: `python -m src.reference.jobs.depth_confound_rare_labels`.
  The numbers will not change; the stamp will.
- **No artifact link.** Unchanged, and now the weakest remaining position:
  `tests/test_checks_can_fail.py` is exactly the thing that would be worth
  linking to.

---

## Priority 0 — a reviewer will see this

### P0-1 · Duplicated word on page 8

**Where:** [main.tex:512-513](main.tex#L512-L513)

Renders on page 8, line 273 of the PDF, in the second-to-last paragraph of main
text:

> "Comparisons against `NaN` return false, so **undefined an undefined interval**
> counted as a coverage *and* a discrimination failure"

**Fix.** Delete the stray `undefined` at the end of line 512.

```latex
holding one. Comparisons against \texttt{NaN} return false, so
an undefined interval counted as a coverage \emph{and} a discrimination failure,
```

**Cost:** −1 word.

---

### P0-2 · The figure legend uses four names the paper never defines

**Where:** [figures/fig_ceiling.pdf](figures/fig_ceiling.pdf) right panel, drawn
by [make_figs.py:80-95](make_figs.py#L80-L95)

Figure 1's right panel is labelled `best4`, `lineage`, `crypt_position`,
`epithelial`. **None of those four strings appears anywhere in `main.tex` or
`sections/`.** The prose says "rarest / middle / coarsest / finest resolution."
A reader cannot map the figure to the argument that depends on it.

Compounding it: §3 calls the rarest rung *"the rarest resolution"*
([main.tex:229](main.tex#L229), [:250](main.tex#L250), [:271](main.tex#L271),
[:282](main.tex#L282)) and §7 calls the same rung *"the finest resolution"*
([main.tex:437](main.tex#L437)). And "rung" — pure internal vocabulary — leaks
into the body once at [main.tex:286](main.tex#L286).

**Fix, three edits.**

1. [main.tex:229-230](main.tex#L229-L230), name them once at first use:

```latex
$|\rho| \geq 0.20$. It never fired on the rarest annotation resolution
(\texttt{best4} in Figure~\ref{fig:ceiling}), and a results sidecar recorded
that resolution as the cleanest of four.
```

2. [main.tex:437](main.tex#L437), `finest` → `rarest`:

```latex
compositional question gets sharpest. At the rarest resolution 28 of 28 patients
```

3. [main.tex:286](main.tex#L286), `rung` → `resolution`:

```latex
40 of theirs. One resolution's surviving rows against another's complete ones is the
```

**Cost:** +4 words.

**Cheaper alternative if the budget bites:** relabel the axis inside
`make_figs.py` instead — `{"best4": "rarest", "lineage": "middle", ...}` — which
costs zero words and fixes the figure at source. This is the better fix if you
have the five minutes to regenerate.

---

## Priority 1 — the claim is not supported by anything in `results/`

### P1-1 · §3's primary-cohort repair numbers have no committed table behind them

**Where:** [main.tex:269-274](main.tex#L269-L274)

Four numbers:
- "moves **66 of the 190** rows carrying a correlation on the primary cohort"
- "by up to **3.9×**"
- "flips the reachability flag … from **64 of 64** … down to 16"
- "put **21 of 166** rows *above their own bound*"

All four need per-arm ρ **on GSE178341**. It does not exist:

| Table | Cohort | Has per-arm ρ? | Has prevalence? |
|---|---|---|---|
| `depth_confound_reference.parquet` (256 rows) | GSE178341 | ✗ `worst_rho` only | ✗ no column |
| `depth_confound_per_arm.parquet` (160 rows) | **GSE132465** | ✓ | ✓ |

And [src/reference/jobs/depth_confound_per_arm.py:112-121](../../src/reference/jobs/depth_confound_per_arm.py#L112-L121)
calls `load_lee_cohort()` — it runs on the replication cohort only. GSE178341's
1.2 GB `.h5` is in `data/manifest.csv` but **is not on disk**, so this cannot be
regenerated before the deadline.

The denominators do check out (190 finite rows, 64 at `best4`), which is why this
reads as solid. The *deltas* do not.

`docs/paper_number_audit.md:126-128` already flagged all four as **unbacked**.
`paper/wmhs/PATCHES.md:11` records the decision — *"A4 (dropped `66 of 190` and
`3.9×`; kept the 64→16 flip)"* — **the sibling paper already removed them.** The
fix simply never propagated to this draft.

**Fix.** Replace [main.tex:268-275](main.tex#L268-L275):

```latex
over the two arms, so the statistic and its bound sometimes came from different
arms. Correcting the pairing flips the reachability flag at the rarest
resolution from ``capable of firing'' on \textbf{64 of 64} combinations down to
16 (Appendix~\ref{app:collapse} sets out the grid). We then made the same error
a third time inside the analysis of the fix, which put rows \emph{above their
own bound} --- a violation visible the moment the points are
```

Then in [sections/appendix.tex:140-143](sections/appendix.tex#L140-L143), drop
190 and 166 since nothing cites them any more:

```latex
differ because it is filtered differently for different questions. 64 is one
resolution across both annotation axes on the primary cohort (GSE178341,
$n = 32$). The per-arm counts out
```

**Cost:** **−18 words.** This is the edit that funds P0-2 and P1-3.

> Keep 64→16: it is a reachability flag computable from prevalence alone, and the
> `best4` rung's prevalences are recoverable. The other three are not.

---

### P1-2 · §4's tier-B abundances mix two reads, and one number is in neither

**Where:** [main.tex:318-327](main.tex#L318-L327)

> "MLH1 at 0.040, SFRP1 and SFRP2 at zero, against **24.9 and 19.4** for GUCA2A
> and GUCA2B … The ratio between the tiers runs from roughly **490-fold to
> 900-fold** across the eight strata"

Re-derived, per-stratum patient medians:

| Quoted | Where it actually comes from |
|---|---|
| MLH1 0.040 | matched read, `lineage` · `opposite_lineage` |
| GUCA2A 24.9 | **unmatched** read, `lineage` · `opposite_lineage` (24.897) |
| GUCA2B 19.4 | **neither read, any stratum.** Nearest: 17.78 unmatched, 17.94 matched |

So two of three come from different tables and the third is unreproducible.

The **490–900-fold** range reproduces on no aggregation I tried:

| Read | Aggregation | Range over 8 strata × 2 A-genes |
|---|---|---|
| unmatched | median | 238 … 6177 |
| unmatched | mean | 457 … 3224 |
| matched | median | 292 … ∞ (MLH1 hits 0) |
| matched | mean | 449 … 4563 |

**The qualitative claim is safe on every read** — tier B is two to three orders
of magnitude below tier A, and the minimum over all four reads is 238-fold.
Only the quoted interval is wrong.

**Fix.** Quote one stratum from one read, and state the floor rather than a
range:

```latex
arm (MLH1 at 0.040, SFRP1 and SFRP2 at zero), against 25.0 and 17.9 for GUCA2A
and GUCA2B, the two abundant compositional controls, on the depth-matched middle
resolution. The ratio between the tiers is never below 200-fold in any of the
eight strata, on the matched read or the unmatched one, so neither the floor nor
the population definition produces it. The contrast asks the estimator to separate
```

**Cost:** −15 words. Every number in it verified against
`results/2026-08-29_4c2a3a9/decomposition_summary_matched.parquet`.

---

### P1-3 · The "forced move" premise has an unaddressed escape

**Where:** [main.tex:120-128](main.tex#L120-L128) and
[sections/appendix.tex:100-125](sections/appendix.tex#L100-L125)

The chain is: leakage invariant bars target genes → nothing remains but
absence-of-opposite-markers → absence reads dropout. Appendix C rebuts count
splitting and reference-based annotators. It never rebuts the obvious third
option: **label the mature state with positive markers that are not targets.**

`config/panel.yaml` tier E holds AQP8, CA1, CA2, CA4, CEACAM7, SLC26A3, FABP1,
PIGR, KRT20, LGALS4, VIL1, SATB2, MUC2, TFF3 — a complete positive
mature-colonocyte marker set. The move is forced *because the panel swallowed
them all*, not because nature forces it.

This is the first thing a single-cell reviewer will ask, and the paper has no
prepared answer. The answer exists — it just is not written down.

**Fix.** One sentence into Appendix C, which is free (appendices do not count):

```latex
\paragraph{Positive markers that are not targets.} The obvious third escape is
labelling the mature state by positive markers held outside the measured set ---
AQP8, CEACAM7, SLC26A3, KRT20 and the rest of the panel's exploratory tier. Our
panel forecloses it by construction, since those genes are in the panel and the
leakage invariant withholds the whole panel rather than the outcome gene alone.
That is a design choice and a defensible one --- a label built from genes the
same disease process regulates re-imports the confound through a second door ---
but it is ours, not a property of the assay, and a panel scoped to a single
outcome gene would not face it.
```

---

### P1-4 · Conformal prediction is in the bibliography and not cited

**Where:** [refs.bib](refs.bib) — `vovk2005` and `angelopoulos2023` are present
and **uncited**. So are `phipson2022`, `burton2006`, `morris2019`.

§7 grounds abstention in Chow's error-reject tradeoff and selective
classification, and stops there. Conformal prediction is the framework a 2026
reviewer will expect to see named in any "return an interval or refuse"
discussion. The reference is already sitting in the file.

**Fix.** [main.tex:427-430](main.tex#L427-L430):

```latex
the error-reject tradeoff \citep{chow1970} and selective classification
\citep{elyaniv2010,geifman2017,angelopoulos2023}, not a bespoke
```

**Cost:** 0 words. Pre-empts a likely review question for one citation key.

---

## Priority 2 — invisible to reviewers, real for the artifact

These matter because the paper's thesis is *a check you cannot falsify is not a
check*. Nobody will catch them in review. They are wrong anyway.

### P2-1 · Two repairs described in the past tense are still live in `src/`

**§3 claims:** *"The diagnostic now persists one row per arm, correlation and
bound and prevalence and depth from the same cells."*

[src/harness/depth_confound.py:189-199](../../src/harness/depth_confound.py#L189-L199)
still computes them as independent maxima:

```python
worst_rho = max(rhos) if rhos else float("nan")
...
ceiling = max(ceilings) if ceilings else float("nan")
```

That is the mispairing bug, unrepaired. The *new* job
(`depth_confound_per_arm.py`) does it correctly, so the paper's sentence is true
of one module and false of the one still shipping.

**§7 claims:** the NaN-scored-as-failure defect, in the past tense.

[src/harness/calibration.py:122-123](../../src/harness/calibration.py#L122-L123)
still has no NaN guard:

```python
covered = (rows["ci_low"] <= truth) & (truth <= rows["ci_high"])
excludes_zero = (rows["ci_low"] > 0) | (rows["ci_high"] < 0)
```

A NaN interval → `False` → counts as a coverage **and** a discrimination failure.
Exactly the described bug.

**Fix.** Either repair both modules, or change the paper's tense. Repairing is
better and is maybe twenty lines:

```python
has_ci = np.isfinite(rows["ci_low"]) & np.isfinite(rows["ci_high"])
covered = has_ci & (rows["ci_low"] <= truth) & (truth <= rows["ci_high"])
excludes_zero = has_ci & ((rows["ci_low"] > 0) | (rows["ci_high"] < 0))
```

…with the abstained rows counted and reported in their own column rather than
folded into a rate. **Add a regression test that fails against the old version**
— the paper's own rule, applied to its own repair.

---

### P2-2 · Figure 1's table was produced by uncommitted code

`results/2026-08-31_0042d33/depth_confound_per_arm.meta.json`:

```json
"git_dirty": true,
"git_untracked": ["src/reference/jobs/depth_confound_per_arm.py"],
```

The script was committed later (`d94e2a1`), but the table was never regenerated
against the clean tree. So the paper's one honest figure rests on a table stamped
dirty — which is Appendix A item 5, recurring inside the paper that reports it.

**Fix.** Re-run the job on the clean tree and re-bundle. Two minutes, and the
provenance check you built will now say what you want it to say.

---

### P2-3 · Table 2 has no table in `results/`

Table 2's numbers live in `submission/FINDINGS.md` and `submission/bench.py`, not
under `results/`. Every other number in the paper resolves to a versioned parquet
with a commit hash and a seed. This one resolves to a markdown file.

**Fix.** Have `submission/bench.py` emit `competitor_bench.parquet` through
`src/common/provenance.py` like everything else.

---

### P2-4 · `_tables.py` picks "newest" by string sort

[_tables.py](_tables.py) resolves the newest results directory
lexicographically, so within one date it picks by commit-SHA string rather than
by time. Live for the `pilot_*` tables (oldest of four same-day runs wins). No
ICBINB figure reads them, so this is not a paper bug — but it is Appendix A's
item 3 still shipping in the tool that draws the paper's figures.

**Fix.** Sort on `utc_timestamp` from the sidecar, falling back to mtime.

---

## Priority 3 — presentation

| # | Issue | Fix |
|---|---|---|
| P3-1 | Figure 1's right panel renders an **empty `epithelial` row** (bound is 0 there, so no ratio exists). Reads as a plotting bug. | Drop the rung from `present` in [make_figs.py:88](make_figs.py#L88), or annotate it "undefined — bound is zero". |
| P3-2 | Figure 2's grey shading is **unexplained**; Figure 3's identical shading is explained. | One clause in the Fig. 2 caption. Captions are inside the page budget — take it from P1-1's savings. |
| P3-3 | Figures 2 and 3 waste the left third of the axis on *n* = 1…20 with no data. | `set_xlim(left=15)` in `../wmhs/make_fig1.py` / `make_fig3.py`. Regenerate there, copy across — never edit in place. |
| P3-4 | §5 is ten lines and its result is "both arms were withdrawn before the analysis ran." It has a section number and no result. | Fold into §8 as a takeaway, or promote Appendix E's simulation into it. Folding frees ~half a page for P4-1 below. |
| P3-5 | Title says "four checks", Table 1 lists five guards. Reconcilable but costs the reader a beat. | Table 1's caption already says "Every one passed. One measured something." Make it explicit: "…the leakage invariant holds; the other four could not fail." |

---

## Priority 4 — the substantive weakness, and the fix is cheap

### P4-1 · §3's blind spot is unpopulated in the data, and it does not have to be

This is the paper's largest exposure. §3 proves a check goes blind below
p = 1.3516%, then reports that **no row in the cohort lies there**. The 58 blind
rows are degenerate — 18 at p exactly 0, 40 at p exactly 1 — where "a threshold
cannot fire" is arithmetic, not a finding. Figure 1's left panel shows an empty
shaded region. The paper concedes it: *"what we cannot claim is that we caught it
biting."*

**The regime is populated. Your own rungs just never reach it.**

`data/raw/lee/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz` — already
on disk, sha256-verified in `data/manifest.csv` — carries **36 published cell
subtypes**, a dozen of them under 1.35% prevalence, including the mature
epithelial populations this study is about.

**I ran it, on two cohorts.** Full-matrix depth pass plus per-(subtype, patient,
arm) Spearman against depth. Two minutes of compute each, no download, no new
data:

**GSE132465 (SMC), 10 paired patients**

| | rows | median \|ρ\| | median \|ρ\|/bound | flagged at 0.20 |
|---|---|---|---|---|
| **p < 1.3516% (cannot fire)** | **309** | 0.035 | **0.422** | **0 (0%)** |
| p ≥ 1.3516% (can fire) | 241 | 0.131 | 0.385 | 72 (30%) |

**GSE144735 (KUL3), 6 paired patients — different centre, already on disk**

| | rows | median \|ρ\| | median \|ρ\|/bound | flagged at 0.20 |
|---|---|---|---|---|
| **p < 1.3516% (cannot fire)** | **261** | 0.052 | **0.504** | **0 (0%)** |
| p ≥ 1.3516% (can fire) | 226 | 0.104 | 0.316 | 49 (22%) |

- **570 rows in the blind band across two cohorts, every one genuinely rare
  rather than degenerate** — against the current draft's 58, all degenerate.
- **250 of the 570 sit above half their attainable bound.** On SMC 31 sit above
  80%; the maxima are **0.994** and **0.963** — rows at 99% and 96% of the
  strongest depth-label association their prevalence permits, reported clean.
- **Zero of 570 flagged**, by construction.
- On SMC the two groups are **statistically indistinguishable in how confounded
  they are relative to what is attainable** (Mann–Whitney p = 0.20). On KUL3 the
  blind rows are **more** confounded than the visible ones (0.504 vs 0.316).
  Same confounding or worse. One group can never be flagged.
- Spread over 36 and 40 subtypes and every patient in both cohorts, so it is not
  one population.
- The mature epithelial populations the decomposition is about are inside the
  band and reported clean: Mature Enterocytes type 1 (median ratio 0.63), type 2
  (0.52), Stem-like/TA (0.61 SMC, 0.90 KUL3), Goblet cells (0.28, 0.39).

**And it replicates.** §4 already argues that a structural failure replicating is
the null, not the finding — here that argument works *for* you: two cohorts,
different centres, agreeing that the check is blind exactly where the rare
populations are.

That converts §3 from *"the theorem says this happens; our data has no rows
there"* into *"here are 309 rows there, 117 with strong confounding, and the
diagnostic reports clean on every one."* It is the difference between a proof
with a null and a proof with a demonstration.

**Honest caveat to state in the paper:** these are the published subtype
annotations, not your maturity label. That is a feature — the claim is about the
*statistic*, not your construct, and demonstrating it on labels you did not build
is stronger, not weaker. Say so in one sentence.

**LANDED** as `src/reference/jobs/depth_confound_rare_labels.py`, emitting
`depth_confound_rare_labels.parquet` through `src/common/provenance.py`:

    python -m src.reference.jobs.depth_confound_rare_labels

It runs both cohorts through `depth_confound_report` — the same function, cells
and depth definition as `depth_confound_per_arm`, so only the label differs and
the two tables are directly comparable.

**Word cost:** roughly 60 for a proper paragraph. P1-1 gives back 18 and P1-2
gives back 15; folding §5 (P3-4) gives back the rest with room to spare.

---

## Suggested order under deadline

| Time | Do |
|---|---|
| 5 min | P0-1 (typo), P1-4 (conformal cite) |
| 20 min | P1-1 (drop unbacked numbers — **do this first, it funds the rest**), P1-2 (tier-B numbers) |
| 20 min | P0-2 via `make_figs.py` relabel + regenerate; P3-1, P3-2 while you are in there |
| 20 min | P1-3 (Appendix C paragraph — free, appendices are unlimited) |
| 2 hr | **P4-1** — the one change that materially raises the paper's ceiling |
| after | P2-1 … P2-4, which nobody will catch and which are wrong anyway |

If you only have an hour: P0-1, P1-1, P1-2, P0-2. Those four remove everything a
reviewer could catch you on.

---
---

# Part 2 — what to change about the project

The paper is a methods post-mortem. The repo is a full research programme.
Most of the repo does not serve the paper, and the parts that do are the newest
and thinnest. That mismatch is the root of several defects above.

## Add

### A1 · The rare-label diagnostic — `src/reference/jobs/depth_confound_rare_labels.py`

P4-1 above. Validated on two cohorts, both already on disk. **This is the single
highest-value change available and it needs no new data.** It converts the
paper's central theorem from an unpopulated claim into a two-cohort
demonstration.

### A2 · Promote GSE144735 (KUL3) to a first-class cohort

`data/raw/lee/GSE144735_*` is on disk, sha256-verified in the manifest, and
**already supported** by [src/estimator/lee_io.py](../../src/estimator/lee_io.py)
— which handles its `Border` class and its 6 excluded patients. 27,414 cells,
6 patients, **all six with matched normal and tumour**, 40 subtypes.

The paper claims two cohorts and calls replication "the standard evidence that a
finding is not an artefact of one dataset." A third paired cohort, at zero
ingestion cost, makes that argument stronger everywhere it appears — and A1
already replicates on it.

Nothing in `results/` reads KUL3 today. That is the largest unused asset in the
repository.

### A3 · `tests/test_checks_can_fail.py` — the paper's thesis, as code

This is the change that would most improve the paper's *standing*, as opposed to
its content.

The paper's remedy is: **test each check against an input that forces it to
fail, and treat a guard with no such test as untested.** The repo has 44 test
files and 231 `pytest.raises` assertions — good discipline, but those assert that
*errors* raise. The guards do not raise; they return **verdicts**. A verdict test
that only exercises the pass path is precisely the paper's subject.

One file, one parametrised test per guard:

| Guard | The input that must trip it |
|---|---|
| leakage invariant | a marker list containing a target gene, in the namespace the guard actually reads |
| depth correlation check | a label perfectly separated by depth at a prevalence where the bound permits firing |
| depth correlation check | the same at p = 0.005, asserting the check **cannot** fire — the bound, as a test |
| tier separation | a synthetic panel where tier D is genuinely retained, asserting the tiers separate |
| recovery curve | an estimator that is *not* a function of the realised statistics, asserting the curve moves |
| abstention gate | an all-NaN interval, asserting it lands in the abstain column and not in a failure rate |
| provenance stamp | a repo with an uncommitted producing script |

Appendix A already tells you five of these inputs — they are the five defects.
Writing them as tests turns the appendix from a confession into an artifact, and
lets §8 say "and here is the suite" rather than "and you should do this."

### A4 · Bring GSE178341 back to disk

The "primary cohort" is in `data/manifest.csv` and **not on the filesystem**
(1.2 GB `.h5`). Every primary-cohort claim is therefore unreproducible today,
which is why P1-1 exists. This is the structural cause of the paper's worst
sourcing gap.

### A5 · Emit `competitor_bench.parquet` into `results/`

Table 2 is the paper's best evidence and the only table that does not resolve to
a versioned parquet with a commit hash and a seed — it resolves to
`submission/FINDINGS.md`. Route `submission/bench.py` through
`src/common/provenance.py` like everything else.

### A6 · An anonymised artifact repo

For a paper whose thesis is *a check you cannot falsify is not a check*, shipping
with no code link is a self-inflicted wound. `git log` shows the URL was
deliberately dropped (`bf51365 paper: drop the artifact URL, and stop promising a
release`) — understandable, but it is now the first thing a reviewer will ask
about. With A3 in place you would have something worth linking to.

---

## Subtract

### S1 · Two of the four annotation resolutions

Appendix G already says it, and no reader will notice:

- `epithelial` is **degenerate by construction** — every scored cell belongs to
  the type in question, so the compositional term is exactly zero.
- `crypt_position` and `lineage` are **the same partition on most patients** —
  "one measurement reported twice rather than two agreeing measurements."

So two of four rungs carry no information, and they are inflating every
denominator in §3. Drop them and:

- the 160 per-arm rows become 80, of which none are degenerate
- the 58-that-cannot-fire disappears entirely, which is *good* — it was the weak
  version of the argument, and A1 replaces it with 570 rows that mean something
- Figure 1's empty `epithelial` row (P3-1) disappears
- "the three medians are not over comparable populations" (P0-2's caption
  problem) becomes two medians over comparable ones

Reporting four resolutions when two are known-degenerate is itself a small
instance of the paper's own thesis.

### S2 · The bulk arm's footprint in the paper

`src/bulk/` is **21 modules and 2,881 lines of test** — TCGA-CDR, survival, purity,
batch confounding, clinical coverage, GSE39582 replication. It produces §5:
ten lines, whose result is *"both arms were withdrawn before the analysis ran."*

That is the worst evidence-to-effort ratio in the paper, and §5 has a section
number and no result. Fold it into §8 as a takeaway (P3-4). Keep the code — it
is fine work and it belongs to the project — but stop giving it a section that
advertises how little came back.

Folding §5 also frees roughly half a page, which is where A1's paragraph lives.

### S3 · The unbacked §3 numbers

P1-1. `paper/wmhs/PATCHES.md` already records the decision to drop them; it
simply never reached this draft.

### S4 · Five uncited bibliography entries

`phipson2022`, `burton2006`, `morris2019`, `vovk2005`, `angelopoulos2023`. Cite
`angelopoulos2023` (P1-4 — it pre-empts a likely review question) and delete the
other four.

### S5 · The pre-registration rhetoric that MS4A12 undercuts

The paper leans hard on pre-registration — "frozen in week 0," "before any
analysis code existed," "an exact-value assertion in the test suite." Then §4
reveals the pre-registered retained control retains 4.1%, because nobody checked
a well-documented fact about MS4A12 in colorectal cancer.

A CRC-literate reviewer will read those two things together and conclude the
pre-registration was ceremonial. The paper's honesty about it is genuinely
admirable and it is still the weakest link in the chain.

**Do not subtract the disclosure.** Subtract the *asymmetry*: every time
pre-registration is invoked as a credential, add the clause that it did not
protect the one thing it was pointed at. §4 already lands this once
("we chose it from the prior expectation and never checked the depletion"). §1
and the ethics statement still invoke pre-registration unqualified.

---

## One framing change worth considering

The scope risk is real: there is no learned model anywhere in this paper. It is
a Kitagawa decomposition, a Spearman threshold, a marker-score labeller and a
bootstrap, at a workshop called **Failure Modes of AI in Biology**. §1's last
paragraph answers with *assay* properties, which is an answer about biology, not
about AI.

The strongest available answer is already in the repo and is currently confined
to a required disclosure. [sections/llm.tex](sections/llm.tex) says:

- an LLM assistant wrote the estimator jobs, the sweep, the diagnostics and the
  figures
- **two of the five defects in Appendix A came from assistant-written code**
- the §3 mispairing was *introduced, repaired incorrectly, then reproduced a
  third time inside the analysis of its own repair*
- two of the reported failures surfaced through a standing adversarial audit by
  an independently prompted agent

That is a 2026-relevant story about AI in biology that nobody else at the
workshop will have: **what happens to scientific checks when the code appears
faster than the scrutiny does.** The paper's own sentence for it is already
written — *"Assistance changes how fast code and analysis appear. It leaves
untouched which checks go structurally blind, and multiplies the cost of writing
one."*

Moving that sentence into §1 and letting Table 1 carry a column for *"was this
guard's code assistant-written?"* would cost very few words and would convert the
paper's weakest reviewer objection into its most distinctive contribution.

I would not rewrite the paper the night before the deadline. I would move one
sentence.
