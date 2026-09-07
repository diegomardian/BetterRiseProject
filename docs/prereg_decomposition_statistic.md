# Pre-registration — the cross-gene statistic for the decomposition

**Written:** 2026-09-06 · **Author:** W1 (Bode) · **Status:** proposed ·
**Closes** the gap recorded in
[prereg_adenoma_decomposition.md](prereg_adenoma_decomposition.md) §5 / RESULT.

> **What this document can and cannot buy, stated first.** It fixes the
> statistic **before** it has been computed on anything — the git history is the
> evidence, and this commit precedes the run. But I have already seen three
> other statistics' results on adenoma, so **for adenoma this is not
> confirmatory and does not claim to be.** It is confirmatory for the
> **replication**, which has not been run and whose substrate has not been
> fetched. On adenoma it buys one thing only: a statistic chosen on stated
> grounds rather than on its answer, applied to a question already answered
> three other ways, so that a reader can see whether the reading moves.

---

## 1 · The gap this closes

`prereg_adenoma_decomposition.md` §5 fixed the comparator **set** — score every
pair of the six-gene panel — and did not fix the **scale**. That mattered:

- The raw `compositional` and `intrinsic` terms are in each gene's own CP10K
  units. ACTB sits near 17.6 and MS4A12 near 2.2, so comparing their magnitudes
  across genes is the error `src/reference/detection_scale.py` exists to stop,
  one estimand over.
- The statistic the RESULT quoted (`share_abs`) was chosen **after** seeing the
  output, and `results/2026-09-06_5f70bb3/` measured what that cost: at `best4`
  the cross-block count is 8/8 on `share_abs`, 3/8 on `share_signed` and **1/8
  on `ratio`**.

## 2 · The choice, and the grounds — which are not "which one won"

**Load-bearing: `log_ratio = log(|intrinsic| / |compositional|)`.**

Chosen on three grounds stated here, none of which is its answer, because its
answer is not yet known to anyone:

1. **The project already made this decision one estimand over, for the same
   reason.** `detection_scale.py` exists because a *difference* of two
   proportions is not comparable across genes at different baselines, and its
   remedy was to move to a **log fold change**. The decomposition presents the
   identical problem — a difference of two terms, each scaling with the gene's
   own abundance — and the identical remedy applies. Picking a *different*
   solution here than the one already adopted there would need a reason, and
   there isn't one.
2. **It is the only candidate that is symmetric and unbounded in both
   directions.** `share_abs` and `share_signed` are bounded in [0,1] and [−1,1],
   so a gene whose intrinsic term dominates by 3× and one where it dominates by
   30× are compressed toward the same value; the log ratio separates them
   linearly. Compression at the top of the range is precisely where the two
   target genes live.
3. **It is well behaved where `ratio` is not.** `ratio` is what
   `docs/NEXT_AVENUES.md` §1a states the collapse in, and it is heavy-tailed:
   a patient whose compositional term is near zero sends it anywhere, which is
   why `adenoma_decomposition_scales.py` had to summarise it by a median with a
   rank-based interval. Taking the log of the magnitude removes exactly that
   pathology while keeping the quantity `ratio` is about.

**Reported beside it, never as the claim:** `share_abs`, `share_signed` and
`ratio`, on the existing agreement table. **A contrast on which `log_ratio` and
at least one other statistic disagree carries no unqualified claim**, the same
rule `adenoma_specificity_disagreements.parquet` already applies.

### What it costs, admitted here

`log_ratio` **discards the sign of both terms.** A gene whose intrinsic and
compositional terms have opposite signs — the two mechanisms opposing rather
than compounding — reads the same as one where they agree. That is a real loss
and it is why `share_signed` stays on the table beside it. It is accepted
because the question this statistic answers is *"which of the two mechanisms is
larger"*, and that is a question about magnitudes.

`log_ratio` is **undefined when either term is exactly zero.** Those rows are
dropped and **counted in the report**, never floored to a small constant. At the
`epithelial` rung the compositional term is exactly zero for every gene by
construction, so **`log_ratio` is undefined at `epithelial` for the whole
panel** — that rung contributes its degeneracy demonstration and no cross-gene
contrast, which is what it was included for.

## 3 · The analysis, fixed

- **Estimand.** Per-patient `log(|intrinsic| / |compositional|)`, from the
  committed decomposition, at every rung and weighting.
- **Contrast.** Every ordered pair of the six-gene panel, paired within patient.
- **Interval.** Student-t over patients. Not the percentile bootstrap:
  `docs/HANDOFF.md` §3a, and n is 43 and 18 here.
- **Blocks, as already claimed.** {KRT8, ACTB, EPCAM, CDX2} and
  {MS4A12, GUCA2A}.

## 4 · What would falsify the existing reading

| branch | consequence |
|---|---|
| `lineage` cross-block 8/8, within-block 0/7 | The `lineage` two-block reading holds on a statistic chosen on stated grounds. It remains **post-hoc for adenoma** and becomes testable on the replication. |
| `lineage` cross-block below 8/8, **or** any within-block contrast excluding zero | **The two-block reading is scale-dependent at `lineage` too**, not only at `best4`, and the RESULT's headline must be withdrawn to the `m_T/m_N` ratio table — which needs no decomposition and no scale choice. |
| `best4` still fails | Expected. Already recorded; changes nothing. |
| `GUCA2A − MS4A12` excludes zero | The "not gene-specific" conclusion becomes scale-dependent. It currently contains zero on 2 of 3 at `lineage` and 3 of 3 at `best4`. |

**The second row is a live possibility and it is why this is worth running.**
`log_ratio` is not `share_abs` with better manners — it is unbounded where
`share_abs` compresses, and the two targets sit at the compressed end.

## 5 · Standing

**For adenoma: a fourth reading of an already-answered question, on a statistic
chosen before it was computed.** It cannot make that result confirmatory and
does not claim to.

**For the replication: confirmatory.** The statistic is fixed here, before the
substrate has been fetched, and the replication's own pre-registration will
reference this document rather than restate it.

---

## RESULT

*Not computed. This document is committed before the run.*
