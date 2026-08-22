# Pre-registration — refined tier-B test (MLH1), gate criterion G2

**Written:** 2026-08-22 · **Author:** W1 (Bode) · **Status:** proposed, awaiting
team ratification · **Closes:** [open decision #10](open_decisions.md)

> **This document is only worth what its timestamp is worth.** It states a
> directional prediction, the analysis that tests it, and what would falsify it,
> *before any target-gene expression has been examined*. Ratify it, amend it, or
> reject it — but do not edit the prediction after results exist. Git history is
> the evidence; a change to the prediction after week 3 should be read as a new
> document, not a correction of this one.

---

## 1 · The claim

Tier B genes are predicted to lose expression **cell-intrinsically** — the mature
cells are still there and have gone quiet — rather than compositionally. MLH1 is
the tier B gene with a known silencing mechanism, so it admits a sharper test
than "does it come out intrinsic on average".

**Prediction.** Among patients whose tumours are MMR-deficient:

| Stratum | Mechanism | Predicted intrinsic MLH1 loss |
|---|---|---|
| `mlh1_methylated` | promoter hypermethylation silences MLH1 transcription | **High** |
| `mlh1_intact_mmrd` | MMRd via MSH2 / MSH6 / PMS2; MLH1 transcription untouched | **Near zero** |
| `mmr_proficient` | no MMR defect | Near zero |
| `mlh1_deficient_unmethylated` | MLH1 protein lost without methylation | **Reported separately, no prediction** |

**Direction is the whole test.** The prediction is not "the two groups differ" —
it is that `mlh1_methylated` shows the larger intrinsic loss. A difference in the
other direction falsifies it.

---

## 2 · Why this is stronger than the MMR contrast already pre-registered

1. **It holds MMR status fixed.** `MLH1Meth` is strictly nested inside MMRd (22 of
   22, zero MMRp), so the methylated-vs-intact comparison sits entirely within
   MMR-deficient patients. It is not a restatement of the MMRd-vs-MMRp contrast.
2. **The negative control is mechanistic, not statistical.** `mlh1_intact_mmrd`
   patients reach the same MSI-H phenotype through a different gene. Same disease
   biology, same selective pressure, MLH1 specifically spared. If intrinsic loss
   is really about MLH1 silencing, it should be absent there — and if our
   intrinsic term fires anyway, the term is measuring something else.
3. **It is a within-tier test.** Tier separation (A compositional, B intrinsic, D
   neither) can be satisfied by an estimator that merely tracks expression level.
   This cannot.

---

## 3 · Cohort composition, fixed as of this document

From `assign_mlh1_strata()` over all 62 patients, emitted as a versioned artifact
by `src/reference/jobs/emit_cohort_table.py` **before any expression was
examined**. That parquet, not this table, is the record; this is a copy for
reading.

| Stratum | n | matched | expected after positivity |
|---|---|---|---|
| `mlh1_methylated` | 22 | **12** | fewer |
| `mlh1_intact_mmrd` | 10 | **7** | fewer |
| `mlh1_deficient_unmethylated` | 2 | 2 | — |
| `mmr_proficient` | 28 | 15 | — |

**C115 and C132 are excluded from the contrast on purpose.** Their IHC reads
"MLH1 and PMS2 deficient" — MLH1 protein lost *without* methylation, so plausibly
a germline MLH1 variant whose transcript may or may not survive nonsense-mediated
decay. On `MMRStatus` and `MLH1Status` alone they are indistinguishable from the
negative-control group, and including them would dilute the arm where near-zero
loss is predicted. They are reported as their own stratum.

---

## 4 · The analysis, fixed in advance

- **Estimand.** Per-patient intrinsic term for MLH1 from
  `src.estimator.kitagawa.decompose_cohort`, at every granularity rung and
  labelling axis the run produces.
- **Contrast.** `mlh1_methylated` against `mlh1_intact_mmrd`. Matched patients
  only — an unmatched patient has no compositional or intrinsic term
  (decision #9).
- **Reported as** a difference with an interval, **not as a hypothesis test.**
  12 against 7 is the ceiling and positivity will reduce both arms, possibly to
  6 against 3. §8.4 puts interaction contrasts at roughly 4x the primary; there
  is no reading of this cohort on which a p-value from it means anything.
- **Weighting.** All three (normal, tumour, doubly robust), never folded together.
- **Patients with `estimability="not_estimable"` are excluded from the contrast
  and counted in the report.** Their intrinsic term is `None`, not `0.0`
  (invariant 1), and dropping them silently would bias the arm that runs out of
  mature cells first — which is the arm the prediction is about.

## 5 · What would falsify it

- `mlh1_intact_mmrd` shows intrinsic MLH1 loss comparable to or greater than
  `mlh1_methylated`.
- The direction reverses across granularity rungs, which would mean the result is
  a labelling artifact rather than biology.
- Fewer than 3 patients per arm survive positivity, in which case **this test
  returns no evidence** — neither for nor against. Say that; do not report the
  point estimate as though the interval were narrow.

## 6 · Standing

**Supporting evidence for G2, not its primary basis.** G2's primary test remains
tier separation across all matched patients. This is a second, mechanistically
sharper line that costs nothing to commit to now and is worth nothing committed
to later.

G2 is a pre-committed gate criterion, so **this needs the team, not just W1.**
