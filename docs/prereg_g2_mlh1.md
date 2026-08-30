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

---

## RESULT — run 2026-08-29, three weeks after this document was written

Executed by `src/reference/jobs/run_g2_mlh1_contrast.py` on the **depth-matched**
decomposition (#24.1). Nothing above was edited.

### The pre-registered direction held

**18 of 18 evaluable combinations as predicted**, every interval excluding zero,
direction consistent across all three usable rungs. Neither §5 falsifier fired.

| `lineage`, matched | methylated | intact-MMRd | difference |
|---|---|---|---|
| `stem_pole`, doubly robust | −0.0172 | +0.0469 | **−0.0641** |
| `opposite_lineage`, doubly robust | −0.0094 | +0.0372 | **−0.0466** |

By §1's criterion — *"direction is the whole test"* — this passes.

### But the predicted magnitude is absent, and the reason is the panel

§1 predicted **high** intrinsic loss in `mlh1_methylated` and **near zero** in
`mlh1_intact_mmrd`. What is there is a methylated arm at **−0.017** and a control
arm at **+0.047** — both close to zero, with the separation coming as much from
the control being positive as from the treatment being negative.

The expression scale explains it:

| gene | mean_normal | mean_tumour | loss within mature cells |
|---|---|---|---|
| GUCA2A | **24.499** | 0.508 | **97.9%** |
| MLH1 | **0.039** | 0.022 | 43.6% |

**MLH1 sits ~600× below GUCA2A — roughly one count per 250,000 UMIs.** The whole
contrast rests on distinguishing 0.039 from 0.022 in per-cell means, in the
regime where detection is largely a function of sequencing depth. At that level
±0.05 is what noise looks like, which is why the control arm reads +0.047 rather
than 0.

### What this says about G2's tier-B failure

**Tier B could not have validated the estimator whatever the biology did.**

Tier B is the *intrinsic* control — the arm whose job is to demonstrate that the
estimator can see silencing when silencing is present. It was populated with a
gene too lowly expressed to measure per-cell in 10x data. G2's tier-B arm
"showing nothing" is therefore a statement about MLH1's abundance in this assay,
not about MLH1 silencing in these tumours.

That is not a failure of this analysis and not a mistake anyone could have caught
before measuring. It is a design constraint discovered by running the design.

### Standing, unchanged from §6

**Supporting evidence for G2, not its primary basis. G2 failed as pre-registered
and this does not change that.** What it adds is the reason one of the three arms
failed, which is more useful than the failure alone.

### Two caveats that belong with the number

1. **`n_intact_mmrd` = 4**, one above §5's floor of 3. Depth matching cost the
   control arm two patients.
2. **The exclusion was not neutral.** Usable fractions: methylated 11 of 20
   (55%), intact-MMRd 4 of 10 (40%) — the control arm was thinned harder, which
   is the asymmetry §4 said to watch for. On the *unmatched* read the rates were
   near-identical (39% against 40%); matching introduced the gap.

### What is NOT claimed

That MLH1 silencing was observed. A directional difference between two arms that
both sit inside the noise band of a barely-detected transcript is consistent with
the mechanism and with several other things. **The GUCA2A result stands on its
own expression level; the control that was meant to corroborate it never had the
dynamic range to do so.**
