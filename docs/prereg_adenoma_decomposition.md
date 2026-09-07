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
> ### Amendment 3 — 2026-09-06, from calibrating the cutpoints, still before any run
>
> **§3.3 said the intrinsic cutpoints were provisional. They are worse than
> provisional: on the project's own pre-registered criteria they are not
> identifiable on either available carcinoma cohort.** The calibration had
> already been run four times (`results/2026-08-31_350e546`, `_593e11e`,
> `2026-09-03_e4a00b3` ×2) and no document read those tables. Re-run here on a
> grid dense enough to resolve the answer (`results/2026-09-06_77fee05`,
> W2_HANDOFF §5 item 6, previously not started):
>
> | cohort | pool | returned | `ok` | `wide` |
> |---|---|---|---|---|
> | SMC | `pooled` | **0/8** | — | — |
> | SMC | `reference` | 8/8 | 70 | **42** |
> | KUL3 | `pooled` | **0/8** | — | — |
> | KUL3 | `reference` | 8/8 | 800 (400–800) | **70** |
>
> The `pooled` draw pool returns **no cutpoint on any seed of either cohort on
> any grid** — max discrimination 0.750/0.795 against a 0.80 target. That is G4,
> reproducible, and not a grid-resolution artefact now that the grid reaches 5.
>
> **What this does to §3.3's rung table, and it is severe for one rung:**
>
> | rung | provisional 50/20 | SMC 70/42 | KUL3 400–800/70 |
> |---|---|---|---|
> | `lineage` (n=44) | 43/44 keep an intrinsic term | 39/44 | 36/44 |
> | `best4` (n=20) | 18/20 | **6/20** | **0/20** |
>
> **`lineage` is robust and `best4`'s intrinsic arm is not.** Whether that rung
> yields any intrinsic estimate at all depends on which carcinoma cohort the
> cutpoint was calibrated against. §6 gains a falsifier for it below.
>
> **The compositional arm is untouched.** It gates on `n_cells_resolved` under
> `COMPOSITIONAL_CUTPOINTS` — decision #22, *not* provisional — and the adenoma
> resolved-epithelium counts are large. So `best4` still contributes a
> compositional point to the curve. That is exactly why §3.3 reports the two
> rules separately, and it is the first time the separation has paid.
>
> **`positivity.CUTPOINTS` is NOT changed by this design.** It is W2's file, and
> on this evidence the honest value is not 70 and not 400 — it is unresolved.
> The reading is reported under **all three** candidate cutpoint sets instead.
>
> ### Amendment 2 — 2026-09-06, from building it, still before any run
>
> **§3.1's gate was mis-specified and is replaced by a sensitivity analysis.**
> Found by running the scoring path on a synthetic cohort whose two arms have
> identical composition by construction, so any Δ(mature fraction) is pure
> artefact. Two things came out, and they point opposite ways:
>
> * **The unresolved share IS endogenous** — driven by the effect under study,
>   not independent of it. With the targets collapsing in the tumour arm its
>   cells carry fewer counts, more fall below the depth target, and the arm
>   shares split **0.229 normal / 0.270 tumour**. With the collapse removed
>   both sit at ~0.25 (0.252 / 0.248). That is the confound §3.1 feared and it
>   is real.
> * **It does not propagate to the estimand.** The mature fraction came out
>   **0.500 / 0.502** against a true 0.500 — the exclusion removes cells from
>   numerator and denominator in a way that roughly preserves the ratio.
>
> So a threshold on `unresolved_arm_gap` gates on a **diagnostic rather than on
> the harm**, and at 0.05 it would fire constantly: the median gap is 0.046
> with the effect and **0.028 under the pure null**, so the tolerance sat inside
> the noise band. A rule that discards a rung's compositional reading on that
> basis would throw away good data — the opposite failure from the one it was
> written to prevent, and one this project has no excuse for after §3a.
>
> **Replaced by: compute the fraction BOTH ways and compare.** `mature_cell_counts`
> already carries `n_cells_epithelial`, so the denominator open decision #14
> *rejected* is free. Both decompositions are run and the comparison is the
> gate. See §3.1.
>
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

**AMENDMENT 2 REPLACES THE THRESHOLD GATE THIS SECTION ORIGINALLY CARRIED.**
It was a tolerance of 0.05 on `|unresolved_share_normal − unresolved_share_tumour|`,
excluding flagged patients and killing a rung above a third. Building the path
showed that rule to be wrong in both halves — see the amendment note above: the
gap runs 0.028 **under a pure null**, so 0.05 sits inside the noise band, and
the mature fraction is not measurably biased by the asymmetry anyway
(0.500/0.502 against a true 0.500). It gated on a diagnostic instead of on the
harm, and it would have discarded good rungs.

**What replaces it: the fraction is computed BOTH ways and the decomposition is
run on both.**

- `frac_mature_*` — denominator `n_cells_resolved`, **open decision #14**, the
  primary and the one the schema table carries.
- `frac_mature_*_all_epithelial` — denominator `n_cells_epithelial`, i.e. the
  unresolved cells counted as immature. **This is the denominator decision #14
  rejected**, computed here as the sensitivity arm. It is free: no extra
  counting, only the other denominator.

**Pre-committed reading rule.** Both decompositions are reported. If the
compositional term's sign and its interval's exclusion of zero **agree** across
the two denominators, the reading does not depend on decision #14 and is
reported without qualification. **Where they disagree, that contrast is
reported as denominator-dependent and no unqualified claim is made from it** —
the same discipline `adenoma_specificity_disagreements.parquet` already applies
to the choice of statistic.

- `unresolved_share_normal` / `unresolved_share_tumour` / `unresolved_arm_gap`
  are still **emitted on every row**, as diagnostics. They are no longer a gate.
  They are worth reading because the asymmetry is endogenous — it is caused by
  the expression change under study — and a future substrate with a larger
  panel share of the library could make it propagate where it does not here.

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

Pre-committed: **the estimability mix is reported beside every rung's summary,
under EVERY candidate cutpoint set — provisional 50/20, SMC-dense 70/42 and
KUL3-dense 400/70 — not under one.** Amendment 3 is why: there is no single
calibrated value to report against, and picking one would present a choice
between two disagreeing cohorts as a measurement.

A rung whose `ok` count is below 5 is labelled **`weak`** in the results and in
any prose quoting it. **A rung that keeps an intrinsic term for fewer than half
its patients under any candidate is labelled `intrinsic_not_supported` at that
candidate**, and no unqualified intrinsic claim is made from it. On present
counts `best4` will carry that label under both calibrated candidates and
`lineage` under none. `best4` will be `weak` at n=4 unless the re-run moves
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
| The compositional term's sign or zero-exclusion **flips between the two denominators** (§3.1) | That contrast is **denominator-dependent**: reported as such, no unqualified claim. |
| The direction reverses across rungs | **A labelling artefact, not biology** — which is what the rung curve exists to detect, and is a reportable negative. |
| `best4` keeps an intrinsic term for 0–6 of 20 patients under the calibrated cutpoints (§3.3, Amendment 3) | **Expected, and pre-committed as not a failure of this design.** `best4`'s contribution to the curve is its **compositional** point; its intrinsic arm is reported with the calibration sensitivity attached and carries no unqualified claim. |

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

## RESULT — run 2026-09-06. Nothing above was edited.

`results/2026-09-06_5791c01/`, over `results/2026-09-06_765eb29/icbi_adenoma.parquet`
(Chen_2021_Cell, 44 patients, four rungs, mature fractions emitted).

### The falsifier did not fire: the estimand is IDENTIFIABLE here

§6's first row — *"the per-patient `i/c` ratios cluster on a single constant
across genes"* — is what killed the decomposition on carcinoma, where the ratio
collapsed onto **−5.85** for every gene on the same labels. At `lineage` the
median ratio runs:

| CDX2 | KRT8 | ACTB | EPCAM | MS4A12 | GUCA2A |
|---|---|---|---|---|---|
| 0.330 | 0.481 | 0.671 | 0.724 | 1.583 | 2.029 |

Nothing near −5.85, and a six-fold spread across genes. **The project's original
estimand is computable on this substrate.** That is the claim §1 was written to
test and it holds.

### And it is a TIER, not a gene — §6's third row, pre-registered as expected

On the intrinsic share, `lineage`, 43 patients, Student-t:

| gene | role | share [95%] |
|---|---|---|
| CDX2 | identity | 0.502 [0.411, 0.593] |
| EPCAM | epithelial | 0.518 [0.435, 0.602] |
| ACTB | control | 0.546 [0.457, 0.636] |
| KRT8 | control | 0.550 [0.459, 0.640] |
| **MS4A12** | **target** | **0.709 [0.639, 0.779]** |
| **GUCA2A** | **target** | **0.715 [0.653, 0.776]** |

**On this statistic** all eight cross-block contrasts exclude zero and neither
within-block contrast does — the four-gene block returns 0 of 6, and
`GUCA2A − MS4A12` is **−0.011**, containing zero. `best4` (n=18) reproduces it.

> ### But the statistic was not pre-specified, and it matters — see "How much of this depends on the scale" below.
> Run on all three defensible scale-free constructions
> (`results/2026-09-06_5f70bb3/`), the cross-block count is **8/8 on
> `share_abs`, 7/8 on `ratio` and 6/8 on `share_signed` at `lineage`** — and at
> **`best4` it is 8/8, 1/8 and 3/8.** The two-block reading is well supported at
> `lineage` and **is not robust at `best4`.**

This is §6's third row exactly: **identifiable, GUCA2A separating from
housekeeping and from CDX2 but NOT from MS4A12 — a tier-level intrinsic result
and explicitly not a gene-specific one.** It was named as the expected outcome
before the run, and §5 exists precisely so that this could not be reported as a
GUCA2A result.

**CDX2 sits with the controls on this statistic too**, so *terminal
differentiation down, intestinal identity retained* now rests on the
decomposition as well as on the corrected specificity reading — **two different
estimands, agreeing.**

### The arithmetic check passed

`epithelial` returned a compositional term of **exactly 0.000** for all six
genes, and an interaction of exactly 0.000. Every epithelial cell is mature at
that rung, so Δf is identically zero and the whole change must land in the
intrinsic term. That is arithmetic rather than data, and it is the cheapest test
that the new mature-fraction code is right. It is also why §2 refused a
`best4`-only reading: a curve whose lower bound is absent cannot show it is a
curve.

### The curve has THREE distinct points, not four

**`crypt_position` collapsed onto `lineage`.** Its mature fraction and mature
cell count are *identical* for **41 of 44 patients**. The run log says why, once
per patient: *"supports only 2 of 3 bins — scores are tied across a quantile
boundary. Using ('crypt_bottom', 'crypt_top')."* The tertile split could not be
formed on this data, so it degenerated to a two-bin split — which is
`lineage`'s median split.

So the granularity curve §2 required is **epithelial → lineage(≈crypt_position)
→ best4**, three points, and this document should not be quoted as reporting
four. It does not weaken the curve's purpose — the split still moves with
annotation resolution, which is what the frozen axes file demands be shown —
but a reader counting rungs would otherwise over-count the evidence.

### §3.3's separation paid, on its first use

**Compositional estimability is `ok` for 20 of 20 patients at `best4`**, where
the intrinsic arm is 4 `ok` / 14 `wide_interval` / 2 `not_estimable` on the
provisional cutpoints and **0 `ok` on either calibrated candidate** (Amendment
3). The two rules genuinely answer different questions, and folding them would
have reported `best4`'s compositional term as unavailable when it is fully
estimable.

**So `best4`'s contribution to the curve is its compositional point.** Its
intrinsic share is reported and reproduces `lineage`'s block structure, but it
carries the Amendment 3 caveat and no unqualified claim.

### How much of this depends on the scale — measured, and it is not all robust

`results/2026-09-06_5f70bb3/`, by `adenoma_decomposition_scales.py`. Every
ordered pair on three scale-free statistics; **74 contrasts disagree between
them, 44 of those cross-block.** Cross-block contrasts excluding zero, out of 8:

| rung | `share_abs` | `ratio` | `share_signed` |
|---|---|---|---|
| `lineage` | **8/8** | **7/8** | 6/8 |
| `best4` | **8/8** | **1/8** | 3/8 |

**At `lineage` the two-block separation survives the choice** — 8/8 and 7/8 on
the two statistics that keep the blocks internally homogeneous (`share_abs` and
`ratio` both return 0 of 7 within-block). `share_signed` returns 4 of 7
within-block, so on that construction the blocks are not clean and its 6/8 is
not evidence for the same structure.

**At `best4` it does not survive.** `ratio` — the form
`docs/NEXT_AVENUES.md` §1a states the identifiability claim in — returns **1 of
8**. So `best4`'s block structure rests on `share_abs` alone, and `share_abs` is
the statistic chosen after seeing the output. **No two-block claim should be
made at `best4`.** That rung's contribution stands as its compositional point
and its estimability result, both of which are scale-free by construction.

**`GUCA2A − MS4A12` contains zero on 2 of 3 at `lineage` and 3 of 3 at
`best4`.** The "not gene-specific" conclusion is the most robust thing here —
which is the right way round, since it is the conclusion that withholds a claim.
The exception is `share_signed` at `lineage` (−0.184, excluding zero), and it is
recorded rather than dropped.

### What is NOT claimed, and one gap in this document

**Seven contrasts are denominator-dependent** and carry no unqualified claim
(`adenoma_decomposition_denominator_disagreements.parquet`). None is a
cross-block contrast on the share statistic, so the reading above does not
depend on open decision #14.

**§5 fixed the comparator SET and did not fix the cross-gene STATISTIC**, and
the section above measures what that cost. The raw terms are in each gene's own
CP10K units, so comparing their magnitudes across genes is the very error the
detection-scale correction was written to stop; a scale-free statistic is
required and **none was pre-specified**. The one the headline quotes was chosen
after seeing the output.

**The consequence is not uniform.** At `lineage` the reading survives every
construction that keeps the blocks internally homogeneous. At `best4` it does
not, and that is a real retraction of the sentence this document first carried
("`best4` reproduces the same structure") — it reproduces it on one statistic
out of three. **A successor design must fix the statistic before the run**, and
the agreement table is what stands in for a pre-specification that does not
exist, not a substitute for one.

**Nothing here bears on survivorship.** §7 stands unchanged: GUCA2A-high cells
having been preferentially destroyed is not transcript-detectable, and a
compositional/intrinsic split does not address it.

### Two bugs the real data found, both in code written for this

Neither appeared on the synthetic fixture. Both are committed with failing
inputs.

1. **The schema write failed outright.** `bootstrap_over_patients` is long-form
   by term, so merging it onto `decompose_cohort`'s wide output fanned every
   patient row into three and carried `term`/`n_boot` into a frozen schema.
   `attach_intrinsic_ci` exists for exactly this and encodes open_decisions
   #10's real choice.
2. **The denominator-disagreement detector fired on noise** — the same defect as
   the threshold gate Amendment 2 replaced it with, reappearing in the
   replacement. It flagged 46 contrasts, nearly all at `epithelial` where the
   compositional term is exactly 0.0 one way and −0.001 the other:
   `np.sign(0.0)` is 0 and `np.sign(-0.001)` is −1. A sign flip now counts only
   where both intervals exclude zero. Seven survive.
