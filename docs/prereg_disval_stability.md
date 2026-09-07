# Pre-registration — is the adenoma result driven by one specimen collection?

**Written:** 2026-09-06 · **Author:** W1 (Bode) · **Status:** proposed ·
**Depends on** [prereg_decomposition_statistic.md](prereg_decomposition_statistic.md)
(`ac7eca1`) and [prereg_adenoma_decomposition.md](prereg_adenoma_decomposition.md).

> **Committed before the split is computed.** The statistic was fixed in
> `ac7eca1`, the decomposition is committed in
> `results/2026-09-06_5791c01/`, and the patient→collection map is derivable
> from the cached obs. Nothing below has been run.

---

## 1 · What this is, and the two things it is not

`Chen_2021_Cell` is **four specimen collections, not one.** From the atlas's
`dataset` column, verified 2026-09-06:

| set | cells | scored patients at `lineage` |
|---|---|---|
| `VUMC_HTAN_discovery` | 40,058 | **15** |
| `VUMC_HTAN_validation` | 55,384 | **15** |
| `VUMC_HTAN_cohort3` | 73,724 | **13** |
| `VUMC_HTAN_CRC` | 8,113 | 0 (carcinoma, excluded upstream) |
| *in both discovery and validation* | | **1** (`HTA11_866`) |

Discovery and validation are described by the source publication as two
independent specimen collections about a year apart, and the patient sets are
**disjoint but for one**. So the 44-patient `lineage` reading can be re-cut into
halves that share no patients and were collected at different times, with the
paired design and the estimand untouched.

**It is NOT an independent replication.** Same lab, same platform, same
COLON MAP population, and — decisively — **the same cells, already analysed.**
A split of data whose pooled answer is known cannot confirm that answer. The
"single-cohort" qualifier on the adenoma result **stays**, and this document
does not license removing it.

**It is NOT a substitute for Becker.** B1 remains the only route to a second
substrate, whatever its cohort problems.

**What it IS: a test of one specific failure mode** — that the finding is an
artefact of a single collection batch. That is a real and common way for a
single-cohort result to be wrong, it is the failure mode this data can address,
and no other available analysis addresses it. If the result holds in two
collections a year apart, batch-drivenness is excluded. Nothing else is.

## 2 · Why this exists now: B1 cannot close the statistic gap

`prereg_becker_replication.md` was the designated closing path for avenue A's
post-hoc-statistic problem. Its cohort does not support it: **only four Becker
donors carry both a polyp and unaffected tissue** (Amendment 2 there). At n=4
this project's own interval fires 18.8% under a true null, and that is the arm
the MLH1 pre-registration refused to take any verdict from.

So the statistic gap needs a closing path that does not depend on Becker. This
is the strongest one the committed data supports, and it is laptop-runnable.

## 3 · The analysis, entirely inherited

Nothing new. That is again the point.

- **Statistic:** `log_ratio = log(|intrinsic| / |compositional|)`, **fixed in
  `ac7eca1`** before any of this, with the agreement rule against `share_abs`,
  `share_signed` and `ratio`.
- **Interval:** Student-t over patients. **Not** the percentile bootstrap —
  and note that splitting costs *power, not calibration*: Student-t is 5.0% at
  every n here, where the percentile bootstrap would be 7.7–8.4%.
- **Rung:** `lineage` only. `best4` carries no claim.
- **Weighting:** doubly robust primary, all three reported.
- **Source:** `results/2026-09-06_5791c01/adenoma_decomposition.parquet`. No
  re-scoring, no cluster, no new cells.

**The shared patient is excluded from both halves.** `HTA11_866` appears in
discovery and validation; assigning it to either would put the same patient on
both sides of a comparison meant to be disjoint. It is dropped, named in the
report, and the halves are therefore **15 and 15**.

**`cohort3` is reported as a third set**, not merged into either. It is 13
patients from a later collection and it costs nothing to carry.

## 4 · What is being tested, exactly

The four contrasts that survived the agreement rule at `lineage` and are the
adenoma result's primary claim:

**`GUCA2A` against `ACTB`, `CDX2`, `EPCAM`, `KRT8`** — unanimous on all four
statistics at n=43.

## 5 · What would falsify it

| branch | consequence |
|---|---|
| All four exclude zero in **both** halves | **Batch-drivenness excluded.** The result is not an artefact of one specimen collection. It remains single-cohort, single-lab, single-platform. |
| All four hold in one half, fail in the other | **Ambiguous by design at this n**, and it must not be read as a failure. See §6 — the halves are 1.8× wider and a real effect can miss. Report both halves and the pooled result; claim nothing new. |
| Any contrast **reverses sign** in either half | **A real problem.** Width does not flip signs. This would say the pooled result is a mixture of opposing collections, and the adenoma reading would need withdrawing to the `m_T/m_N` ratio table. |
| Fewer than 4 contrasts hold in either half **and** the effect sizes shrink toward zero in one | Consistent with the pooled result being driven by the other. Report as such. |

## 6 · The power cost, stated before the numbers

At n=15 the Student-t half-width is **1.80×** what it is at n=43
(t=2.145/√15 against t=2.018/√43); at n=13, **1.96×**.

**So a contrast whose n=43 interval sits within a factor of ~1.8 of zero can
fail in a half for width alone.** Of the four, `GUCA2A − CDX2` was the largest
at n=43 on `share_abs` (+0.213) and is the most likely to survive; the smallest
is the most likely to miss. **A miss is therefore weak evidence and a
sign-reversal is strong evidence**, which is why §5 treats them differently.

This asymmetry is stated here so that a half-failure cannot later be read as a
refutation, and a half-success cannot be read as independent confirmation.

## 7 · Standing

**A stability check, pre-registered, on a statistic fixed beforehand and a split
by a variable unrelated to the outcome.** It cannot make the adenoma result
confirmatory — the same cells produced it. It can exclude one specific way of
being wrong, and that is worth having and worth not overselling.

---

## RESULT — computed 2026-09-06, after this document was committed

`results/2026-09-06_705dd5b/`. The split came out **15 / 15 / 13** exactly as §1
predicted, with `HTA11_866` excluded from both halves per §3.

### Verdict: AMBIGUOUS AT THIS N — §5's third branch

`log_ratio` (load-bearing, fixed in `ac7eca1`), GUCA2A against each control:

| contrast | discovery (15) | validation (15) | cohort3 (12) | pooled (42) |
|---|---|---|---|---|
| − ACTB | **+0.881 \*** | **+0.952 \*** | **+0.899 \*** | **+0.912 \*** |
| − CDX2 | +0.648 | **+1.442 \*** | +1.149 | **+1.075 \*** |
| − EPCAM | +0.833 | **+1.520 \*** | +0.834 | **+1.079 \*** |
| − KRT8 | +1.004 | **+1.027 \*** | +0.545 | **+0.881 \*** |

`*` excludes zero. Validation holds 4 of 4; discovery and cohort3 hold 1 of 4.

### The strong falsifier did not fire, and that is the substantive result

**Zero sign reversals, on all four statistics, in all three collections.** §5
reserved that branch as the one width cannot explain, and nothing tripped it.

**Every point estimate is positive and they are mutually consistent** —
discovery 0.648 to 1.004, validation 0.952 to 1.520, cohort3 0.545 to 1.149,
against a pooled 0.881 to 1.079. The halves agree on direction and on rough
magnitude. What they differ on is whether the interval clears zero, which is
exactly what §6 said to expect at 1.80× the width.

**`GUCA2A − ACTB` holds independently in all three collections.** That is the
only contrast to do so, and it is the one with the least between-half spread.

### What this does and does not establish

**Does:** the four contrasts are not carried by one collection batch. Two
independent collections a year apart, plus a third, all give the same sign and
comparable magnitude, and no half contradicts the pooled reading.

**Does not:** exclude batch-drivenness *formally*. §5's first branch required 4
of 4 in both halves and that did not happen. On the pre-registered reading this
is **ambiguous, not confirmed** — and §5 fixed in advance that it must not be
read as a failure either.

**And it remains what §1 said it was.** Not an independent replication: same
lab, same platform, same population, same cells. **The single-cohort qualifier
on the adenoma result stays.**

### One inconsistency this run exposed, fixed before the numbers above

The first version applied a Student-t interval to every statistic, including
`ratio` — which `adenoma_decomposition_scales.py` summarises by the **median
with a rank-based interval**, because it is heavy-tailed and one patient with a
near-zero compositional term sends its mean anywhere. That mismatch produced
**the only two sign flips in the first run**, both on `ratio`, both artefacts of
the summary rather than properties of the data. With the pre-registered
treatment applied consistently, **all four statistics show zero flips.**

Worth recording because a sign flip is §5's *strongest* branch: an inconsistent
summary in the reading job would have tripped the one verdict that width cannot
explain.
