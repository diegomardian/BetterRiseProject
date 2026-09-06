# Pre-registration — the decomposition on adenoma, at the rungs it supports

**Written:** 2026-09-06 · **Author:** W1 (Bode) · **Status:** proposed, awaiting
team ratification · **Closes** `docs/NEXT_AVENUES.md` §1a / proposal A.

> **This document is only worth what its timestamp is worth.** It fixes the
> estimand, the population each term is computed on, the interval, the
> comparator set and the falsifiers **before any decomposition has been run on
> this cohort.** Ratify it, amend it, or reject it; do not edit the prediction
> after results exist.

> ### Amendment 1 — 2026-09-06, before ratification and before any run
>
> Three corrections from review. **No result exists yet**, so these are
> amendments to a proposal, not edits to a prediction; the original text is in
> `099984f` and what changed is stated here rather than silently applied.
>
> 1. **§4 said the Student-t interval "needs nobody's approval". That was
>    false.** `coerce_results` raises `SchemaViolation` on any column outside
>    `REQUIRED_COLUMNS`, citing invariant 3 by name, and the decomposition
>    writes through exactly that path (`run_decomposition.py:198`). So
>    `t_ci_low`/`t_ci_high` cannot be columns of the schema frame. §4 now
>    pre-commits a **companion table** instead, which really does need nobody's
>    approval.
> 2. **§1's agreement claim was too tight and its geometric mean needed a
>    qualifier.** The unfiltered geometric mean is *degenerate* for both target
>    genes. Restated with the filter named and the real bound.
> 3. **§3.1 presented a settled decision as a new pre-commitment, and §3.3
>    omitted the compositional arm's own cutpoint.** The denominator rule is
>    **open decision #14**, already implemented in
>    `labels.mature_cell_counts`, which also already emits the per-arm shares
>    §3.1 asked for. What is genuinely new is only the *gate*. And the
>    compositional term has its own estimability rule — `COMPOSITIONAL_CUTPOINTS`,
>    decision #22 — which §3.3 never named.

---

## 1 · Why this is the project's deliverable and not another reading

The README's bar is *"decompose differentiation-marker loss into compositional,
cell-intrinsic, and not estimable — per patient, at several annotation
resolutions."* **That has never been produced on a substrate where it is
identifiable.** On carcinoma the estimator collapses algebraically:

    intrinsic / compositional = (f_N / Δf) × (m_T/m_N − 1)

As a gene's surviving per-cell mean → 0 the bracket → −1 and the ratio becomes
`−(f_N/Δf)`, a property of the *cell fractions*, identical for every gene on the
same labels. On GSE178341 that constant is −5.85 and five panel genes are
indistinguishable by it. The week-0 falsification rule fired and its
consequence — no gene-specific mechanism claim from the decomposition — stands.

**On adenoma the collapse does not fire.** Re-checked 2026-09-06 against
`results/2026-09-05_d869bdd/icbi_adenoma.parquet`:

| gene | m_T/m_N (`lineage`) | bracket |
|---|---|---|
| KRT8 | 0.946 | −0.054 |
| ACTB | 0.834 | −0.166 |
| CDX2 | 0.791 | −0.209 |
| EPCAM | 0.737 | −0.263 |
| MS4A12 | 0.383 | −0.617 |
| GUCA2A | 0.374 | −0.626 |

An order of magnitude of spread and nothing near −1.

**The cohort ratios above are pooled, and the per-patient summaries agree with
them — with the filter named.** Over per-patient ratios that are **positive and
finite**, the median and geometric mean agree with the pooled ratio to within
**0.10** (worst: MS4A12, median 0.480 against pooled 0.383), and the *ordering* —
the two targets far below the other four — is identical on all three.

The filter is not cosmetic. **Three of 44 patients have an exactly-zero arm mean
for GUCA2A and for MS4A12**, so the unfiltered geometric mean of those two genes
is **degenerate**: undefined where a denominator is zero, zero where a numerator
is. That is the same estimability problem §3.3 is about, showing up in a summary
statistic, and it is why the filter has to be stated rather than assumed.

*The arithmetic mean of per-patient ratios agrees with nothing* — it reads
MS4A12 at 0.841 against a pooled 0.383, because a ratio with a small denominator
explodes and MS4A12's normal-arm mean tops out at 8.0 CP10K. It is the obvious
first thing to compute, it is an estimator artefact, and it is not a finding.

**And this needs nothing that failed.** No MLH1, no premise resolution, no
instrument sensitivity. `docs/prereg_g2_mlh1_within_stratum.md` closed the
instrument-validation question negatively; this avenue does not depend on it.

---

## 2 · The estimand, and the rungs — a curve, not a point

**Per-patient Kitagawa split** from `src.estimator.kitagawa.decompose_cohort`:
`compositional`, `intrinsic`, `interaction`, under all three weightings
(`normal`, `tumour`, `doubly_robust`), never folded together (invariant 7).

**Reported as a CURVE across granularity rungs.** `config/labeling_axes.yaml` is
frozen and says so:

> *The split is annotation-relative by construction: a compositional change
> becomes an expression change purely by redrawing cluster boundaries. Report
> the split as a CURVE across all four rungs. A single point estimate would
> present a modelling choice as a measurement.*

**This is the one place §1a's own framing had to be corrected.** §1a is titled
"the decomposition at `best4`". A `best4`-only reading is exactly the single
point estimate the frozen file forbids, and it would present the choice of
annotation resolution as a measurement of biology. **This runs at every rung the
cohort supports.**

| rung | committed patients | status in this design |
|---|---|---|
| `epithelial` | not yet scored | **scored, and expected degenerate** |
| `lineage` | 44 | scored |
| `crypt_position` | not yet scored | **scored** |
| `best4` | 20 | scored |

`epithelial` is included *because* it is degenerate — every epithelial cell is
mature there, so Δ(mature fraction) measures only epithelium-vs-rest and
everything inside the epithelium lands in the intrinsic term. `RUNG_SPECS`
calls that "the lower bound of the granularity curve and it is supposed to look
degenerate — that is what it demonstrates." **A curve whose lower bound is
absent cannot show it is a curve.**

---

## 3 · The three populations, fixed here, because two of them are not obvious

`decompose_cohort` needs ten columns. Six are identifiers and two —
`mean_normal` / `mean_tumour` (as `cp10k_*`) and `n_cells_mature` (as
`n_tumour`, correct: the tumour arm is what `classify_estimability` reads) — are
already on the committed rows. `frac_mature_normal` / `frac_mature_tumour` are
not.

**But they are already implemented, and Amendment 1 corrects this section's
claim that they were not.** `src.reference.labels.mature_cell_counts` returns,
per (patient, **tissue**, axis, rung): `n_cells_mature`, `n_cells_epithelial`,
`n_cells_unresolved`, `n_cells_resolved`, **`mature_fraction`** and
**`unresolved_fraction`**. Tissue *is* the arm, so both fractions this design
needs are a pivot of an existing function, not a new computation. **The build is
to CALL it in the adenoma path and emit its output beside the deltas** — not to
derive anything.

### 3.1 The denominator, and the quarter of the epithelium that leaves it

**The denominator is the RESOLVED epithelium of that arm** — cells whose label
is one of the rung's own bins, excluding `non_epithelial` and excluding
`unresolved_depth`.

**This is not a new pre-commitment and the first version of this section wrongly
implied it was.** It is **open decision #14**, already settled and already
implemented: `mature_cell_counts` computes `n_cells_resolved =
n_cells_epithelial − n_cells_unresolved` and divides by it, on the stated
ground that *"a cell that could not be labelled is not a cell measured to be
immature."* This design **adopts** that decision; it does not make it.

**And the excluded share is not data.** Every one of the thirty patients in the
MLH1 run reported *exactly* 25.0% of epithelial cells as `unresolved_depth` —
351/1405, 404/1614, 1103/4411, thirty for thirty. That is `DEPTH_QUANTILE = 0.25`:
the depth target is the 25th percentile, so a quarter falls below it by
construction.

**The risk this creates, stated before the numbers.** The cut is applied in
`assign_labels`, per patient, over **both arms pooled**, and **before depth
matching**. If the arms differ in depth pre-matching, the exclusions concentrate
in the shallower arm, the two denominators differ for a purely technical reason,
and **the compositional term measures sequencing depth.** That is the failure
mode this arm is most exposed to and no committed table can rule it out —
`depth_ratio` runs 0.97–1.05 on all 44 patients but is measured *after*
matching.

**The gate below is what IS new here.** The per-arm shares already exist as
`unresolved_fraction`; nothing has ever compared them between arms or acted on
the difference.

- `unresolved_fraction` is pivoted to `unresolved_share_normal` /
  `unresolved_share_tumour` and **emitted on every row**, per patient per rung.
- **Tolerance, fixed now: 0.05.** If `|unresolved_share_normal −
  unresolved_share_tumour| > 0.05` for a patient, that patient's
  **compositional term is reported with `depth_confounded = True`** and is
  excluded from the cohort summary, counted in the report. Not dropped
  silently — the arm that runs out first would otherwise bias the result in the
  direction of the hypothesis.
- If **more than a third of patients** at a rung are flagged, **that rung's
  compositional reading is NOT REPORTED**, and the reason is named.

### 3.2 Which population each term is computed on

`rows_for_patient` selects mature cells, depth-matches **those**, and computes
`cp10k_*` on the survivors. The epithelial denominator cannot be in that
matched set, because matching is defined on the mature cells alone. So:

**Pre-committed: `mean_*` is post-matching; `frac_mature_*` is pre-matching, on
the QC-surviving naive epithelium of that arm.** Both are recorded on the row
(`matched_basis="mean"`, `unmatched_basis="fraction"`).

**This is a known asymmetry and it is declared rather than discovered.**
`build_decomposition_summary` records learning the same lesson the hard way —
matching a wider population and intersecting afterwards leaves the analysed
subset unmatched. The alternative, matching the whole epithelium, would change
what `cp10k_*` means and make this reading non-comparable with every committed
table. **The asymmetry is the lesser cost and it is on the row.**

### 3.3 Estimability — TWO rules, reported separately, never folded

**Amendment 1 adds the compositional half; the first version named only the
intrinsic one.** `src/harness/positivity.py` is explicit that these are two
questions:

| term | gates on | rule | standing |
|---|---|---|---|
| **intrinsic** | `n_cells_mature` | `CUTPOINTS`, `classify_estimability` | **PROVISIONAL** — `execution_plan.md` §4, awaiting W2 calibration |
| **compositional** | `n_cells_resolved` | `COMPOSITIONAL_CUTPOINTS`, `classify_compositional_estimability` | **decision #22**, pre-committed 2026-08-27, *not* provisional |

Both are `ok=50, wide=20`; the symmetry is *"a default, not a finding"* and the
two must be recalibrated separately. `estimability_verdicts` reports them
separately and **this design does not fold them into one column.** A row may
have an estimable compositional term and an unestimable intrinsic one — that is
the ordinary case at a starved rung, and invariant 1 governs both: undefined is
`None`, never `0.0`.

**The compositional arm's estimability cannot be pre-computed here**, and that
is worth saying plainly: `n_cells_resolved` is not on any committed table, so
unlike the intrinsic mix below it will be known only at run time. It is
pre-committed to be **reported**, per rung, in the same table as the intrinsic
mix.

#### The intrinsic mix, which at `best4` is most of the cohort

`CUTPOINTS` is `ok=50, wide=20`, and its own `source` field says
**`provisional (execution_plan.md §4)`** — never calibrated. Applied to the
committed mature-cell counts:

| rung | ok | wide_interval | not_estimable |
|---|---|---|---|
| `lineage` (n=44) | 38 | 5 | 1 |
| `best4` (n=20) | **4** | **14** | 2 |

**At `best4`, 70% of patients are `wide_interval` and only 4 are `ok`.** That is
not a reason to skip the rung — the intrinsic term is still written there — but
it is a reason not to quote `best4` as though it were `lineage`.

Pre-committed: **the estimability mix is reported beside every rung's summary**,
and a rung whose `ok` count is below 5 is labelled **`weak`** in the results and
in any prose quoting it. `best4` will be `weak` at n=4 unless the re-run moves
it. **`not_estimable` rows carry `intrinsic = None`, never `0.0`** (invariant 1,
enforced by `src/schema.py`'s writer assertion, not by review); their
compositional term is still estimable and the row is **not dropped**.

---

## 4 · The interval

**Two intervals, both reported, neither replacing the other.**

`ci_low` / `ci_high` **stay as the schema decided.** `docs/open_decisions.md`
#10 puts `bootstrap_over_patients`'s percentile band in that slot — W2 proposed,
W4 confirmed 2026-08-22 — and `src/estimator/` is W4's under CONTRIBUTING §2.
**This design does not change a settled decision or another workstream's file.**

**A Student-t interval is carried BESIDE it, in a COMPANION TABLE** —
`adenoma_decomposition_t_intervals.parquet`, keyed by the schema's
`KEY_COLUMNS` — because the schema band is miscalibrated at these n and by a
known amount
(`docs/HANDOFF.md` §3a — the percentile bootstrap of a mean is
`z·sqrt((n−1)/n)/t(n−1)` times the correct width, a function of n alone):

| rung | n | schema band's real false-positive rate |
|---|---|---|
| `lineage` | 44 | **5.9%** |
| `best4` | 20 | **7.1%** |

**Amendment 1 corrects how this is done, and the claim made for it.** The first
version put `t_ci_low` / `t_ci_high` on the schema frame and called that
"additive, needs nobody's approval." **It is neither.** `coerce_results` raises
`SchemaViolation` on any column outside `REQUIRED_COLUMNS` — *"The schema is
frozen — adding a field needs a PR with two approvals (CLAUDE.md invariant 3)"* —
and the decomposition writes through exactly that call
(`run_decomposition.py:198`). A Student-t column on the schema frame does not
need approval so much as it **fails to write at all**.

So: **a separate versioned parquet, keyed by `KEY_COLUMNS`**, in the same shape
the ICBI jobs already use for non-schema artifacts. The schema frame stays
percentile-compliant, no frozen file is touched, and the join is exact. *That*
needs nobody's approval.

**The alternative, if the team wants the columns real:** a schema amendment PR
adding `t_ci_low` / `t_ci_high` — two approvals and a written reason under
invariant 3. Heavier, and it makes every downstream reader aware of the second
interval rather than leaving it in a companion. **Decide this at ratification,
not after the split exists**, because the choice determines which table the
paper's numbers are quoted from.

**Every cross-rung or cross-gene claim in the write-up is made on the Student-t
interval**; the percentile band is carried for schema compliance and
comparability with committed tables.

---

## 5 · The comparator set — the whole design, and it is fixed here

**Score the full panel in the same run, and report every pair.** Not
`GUCA2A − X`.

`GUCA2A`, `MS4A12` (targets) · `CDX2` (intestinal identity) · `EPCAM`
(epithelial) · `ACTB`, `KRT8` (housekeeping).

**Why this is load-bearing.** Run against housekeeping alone, an intrinsic term
for GUCA2A will look like a clean gene-specific result. §1a says it and the
corrected specificity table showed it: the panel separates into **two blocks**,
and GUCA2A is not separable from MS4A12 within its own. **Identifiable is not
the same as gene-specific**, and a decomposition that only ever compared the
target to a housekeeping gene could not tell the difference.

The first specificity table made exactly this mistake — it reported only
`GUCA2A − X`, so the claim that intestinal identity is retained, which is a
statement about where **CDX2** sits relative to the **controls**, had no row
behind it. See `docs/HANDOFF.md` §6d.

### What is explicitly OPEN and assumed in neither direction

**Whether GUCA2A separates from MS4A12.** The two rungs disagree on the ratio
scale — `lineage` 0.374 / 0.383 (indistinguishable), `best4` 0.517 / 0.385
(further apart) — and neither has an interval. §1a's "again indistinguishable"
is a `lineage` statement and **must not be carried to `best4`**; nor may its
negation. **This reading is what would settle it**, and it is pre-registered as
open.

---

## 6 · What would falsify it, and what each branch commits to

| branch | consequence |
|---|---|
| The bracket collapses toward −1 after all — the per-patient `i/c` ratios cluster on a single constant across genes | **The estimand is not identifiable here either.** Same verdict as carcinoma, on a second substrate, and the decomposition arm of this project closes for good. |
| Identifiable, and the intrinsic term separates GUCA2A from **housekeeping and CDX2 and MS4A12** | The strongest available result: gene-specific intrinsic loss. |
| Identifiable, GUCA2A separates from housekeeping and CDX2 but **not from MS4A12** | **The expected outcome on current evidence.** A tier-level intrinsic result — terminal differentiation down, identity retained — and explicitly **not** a GUCA2A-specific one. |
| The compositional term is depth-confounded in more than a third of patients at a rung (§3.1) | That rung's compositional reading is **not reported.** |
| The direction reverses across rungs | **A labelling artefact, not biology** — which is what the rung curve exists to detect, and is a reportable negative. |

**No result from this is a mechanism claim about survivorship.** GUCA2A-high
cells having been preferentially destroyed is not transcript-detectable, and a
compositional/intrinsic split does not bear on it. That caveat survives every
branch above.

---

## 7 · Standing

**This is the project's original deliverable, and the first substrate on which
it can be produced.** It does not revive G2, it does not depend on the
instrument's sensitivity, and it does not rest on the premise resolving — it
rests on the algebra not collapsing, which is measured and committed.

It is a **new primary result**, so it needs the team, not just W1.

---

## RESULT

*Not run.* Requires the ICBI atlas (cluster-only) and a re-run of the adenoma
scoring path that **calls `labels.mature_cell_counts`** and emits its output
(`mature_fraction`, `unresolved_fraction`, `n_cells_resolved`) pivoted by arm,
plus `epithelial` and `crypt_position` alongside `lineage` and `best4`.

**One decision is owed at ratification, not after:** §4's companion table
against a schema amendment. It determines which artifact the paper's intervals
are quoted from.

*Nothing above may be edited when it is.*
