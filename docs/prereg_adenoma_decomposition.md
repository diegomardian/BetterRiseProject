# Pre-registration — the decomposition on adenoma, at the rungs it supports

**Written:** 2026-09-06 · **Author:** W1 (Bode) · **Status:** proposed, awaiting
team ratification · **Closes** `docs/NEXT_AVENUES.md` §1a / proposal A.

> **This document is only worth what its timestamp is worth.** It fixes the
> estimand, the population each term is computed on, the interval, the
> comparator set and the falsifiers **before any decomposition has been run on
> this cohort.** Ratify it, amend it, or reject it; do not edit the prediction
> after results exist.

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

An order of magnitude of spread and nothing near −1. The median and geometric
mean of the per-patient ratios agree with these cohort ratios to within ~0.08.
(*The arithmetic mean of per-patient ratios does not* — it reads MS4A12 at 0.841,
because a ratio with a small denominator explodes. That is an estimator
artefact, it is the obvious first thing to compute, and it is not a finding.)

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
already on the committed rows. **`frac_mature_normal` / `frac_mature_tumour` are
not, and the whole compositional arm is those two numbers.**

### 3.1 The denominator, and the quarter of the epithelium that leaves it

**Pre-committed: the denominator is the RESOLVED epithelium of that arm** —
cells whose label is one of the rung's own bins, excluding `non_epithelial` and
excluding `unresolved_depth`.

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

**So it is measured and gated:**

- `unresolved_share_normal` and `unresolved_share_tumour` are **emitted on every
  row**, per patient per rung.
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

### 3.3 Estimability, which at `best4` is most of the cohort

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

**A Student-t interval is carried BESIDE it**, in `t_ci_low` / `t_ci_high`,
because the schema band is miscalibrated at these n and by a known amount
(`docs/HANDOFF.md` §3a — the percentile bootstrap of a mean is
`z·sqrt((n−1)/n)/t(n−1)` times the correct width, a function of n alone):

| rung | n | schema band's real false-positive rate |
|---|---|---|
| `lineage` | 44 | **5.9%** |
| `best4` | 20 | **7.1%** |

Additive, needs nobody's approval, and lets a reader see both. **Every
cross-rung or cross-gene claim in the write-up is made on the Student-t
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
scoring path emitting `frac_mature_*`, `unresolved_share_*`, and the two
additional rungs.

*Nothing above may be edited when it is.*
