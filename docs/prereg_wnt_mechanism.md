# Pre-registration — does the differentiation fall track Wnt, inside surviving mature cells?

**Written:** 2026-09-06 · **Author:** W1 (Bode) · **Status:** proposed, awaiting
team ratification · **Follows** [prereg_adenoma_decomposition.md](prereg_adenoma_decomposition.md).

> **What this is for.** Avenue A established *that* GUCA2A's loss is
> cell-intrinsic in adenoma — the cells are still there and have turned their
> output down. It says nothing about *why*. This is the cheapest available test
> of a mechanism, on cells already scored, and it is Stage 3's per-cell test
> which has never been run on human tissue.

---

## 1 · The claim

Adenoma is a Wnt-driven lesion — APC loss is the initiating event in FAP and in
most sporadic polyps. If the terminal-differentiation programme is being held
down by Wnt signalling rather than merely lost, then **among cells that are
still labelled mature, the ones with more Wnt target activity should have less
differentiation output.**

**Prediction.** Within the mature cells of the polyp arm, the per-cell Wnt-target
score is **negatively** associated with GUCA2A and MS4A12 expression, after
conditioning (§3), and **not** with ACTB or KRT8.

**The controls carry the test.** A negative association with the targets and
none with housekeeping is the only pattern that means anything. An association
with everything is a technical artefact, and §5 says so in advance.

## 2 · The Wnt score, and invariant 8

**Invariant 8 is explicit and this design follows it exactly:** *"CTNNB1 /
TCF7L2 transcript level is not Wnt activity. Use a target signature (AXIN2,
NKD1, RNF43, NOTUM, TCF7); drop ASCL2/LGR5 when the stem axis is in play."*

**Signature: AXIN2, NKD1, RNF43, NOTUM, TCF7.** Five genes, no CTNNB1, no
TCF7L2, and **no ASCL2 or LGR5** — the stem axis *is* in play (`AXIS =
"stem_pole"` uses LGR5 and ASCL2), so including them would put the labeller's
own markers into the predictor and make the test circular by construction.

**Leakage checked, not assumed.** None of the five appears on the frozen panel
(`config/panel.yaml`) or in any labelling axis (`config/labeling_axes.yaml`).
Verified 2026-09-06. A guard asserts it at run time in the shape of
`build_signature()`'s, and the failing input is committed.

*One thing worth knowing about the panel:* **SFRP1 and SFRP2 are tier B panel
genes and they are Wnt antagonists.** They are not in the six-gene scoring set
and not in this signature, so nothing is circular — but a future design that
scored them alongside a Wnt score would be, and this is where that is written
down.

**Invariant 8 also requires reporting the Wnt/stem correlation when the stem
axis is in play.** It is reported, per patient, as a first-class output and not
a diagnostic. If the Wnt score is simply the stem score by another name, this
test is measuring maturity and §5's fourth row applies.

## 3 · The estimand, and the two things it must be conditioned on

**Per patient**, within the mature cells of the polyp arm: the **partial
Spearman correlation** between the Wnt-target score and each panel gene's
expression, **conditioning on**

1. **the maturity score** (the `stem_pole` axis score the labeller already
   computes), and
2. **log library depth.**

Then aggregate across patients — **the patient is the unit of inference**
(invariant 5) — with a **Student-t** interval, not the percentile bootstrap
(`docs/HANDOFF.md` §3a; n is 44 here).

### Why each conditioner is not optional

**Maturity.** Both the Wnt score and GUCA2A vary with differentiation state, and
the mature bin is a *bin*, not a point — there is residual maturity variation
inside it. An unconditioned correlation would recover maturity and report it as
Wnt. This is the single most likely way to get a false positive here.

**Depth.** Two per-cell expression scores in the same cell are correlated
through library size even after CP10K normalisation, because detection is
depth-dependent. **This is the reason the housekeeping controls exist**: ACTB
and KRT8 are the direct measurement of how much correlation survives
conditioning for purely technical reasons. Whatever they show is the floor, and
a target-gene correlation is only interpretable as the amount by which it
exceeds that floor.

**Spearman rather than Pearson**, because per-cell counts are zero-inflated and
heavy-tailed and a Pearson correlation over them is dominated by the few cells
that fired — the same argument that put this project on detection rather than
means in the first place.

## 4 · Cohort and scope

`Chen_2021_Cell`, the same 44 patients, the same cells, the same labels as
`results/2026-09-06_765eb29/`. **`lineage` only.** Not `best4` — its intrinsic
arm carries no claim (`prereg_adenoma_decomposition.md` Amendment 3 and RESULT),
and running a mechanism test against a rung whose result is retracted would be
asking why something is true that has not been established.

**The polyp arm is where the prediction lives**, and the normal arm is reported
beside it as the within-patient reference: Wnt activity should be low and
unstructured there, so a correlation of the same size in normal tissue is
evidence the whole thing is technical.

**This needs a cluster run.** No committed table carries per-cell values — they
are all per (patient, gene) aggregates — so this is a new pass over the atlas
already on disk. It needs **no new data** and no download, which is what makes
it the cheap option, but it is not laptop-runnable.

## 5 · What would falsify it, and what each branch commits to

| branch | verdict |
|---|---|
| Targets negative, housekeeping not, after conditioning | **Consistent with Wnt holding the programme down.** Associational, in one cohort, and not causal — see §6. |
| Targets and housekeeping both negative | **Technical.** A residual depth or quality gradient the conditioners did not remove. No mechanism claim; report the size as a measurement of the floor. |
| Nothing anywhere | **No detectable Wnt association at this depth.** A real negative: the mature cells' output does not track Wnt target activity within patients. |
| Targets **positive** | Falsifies the prediction as stated. Report it; do not re-describe it as "Wnt-independent". |
| Wnt score correlates with the maturity score above 0.7 | **The test is measuring maturity.** Report the correlation and withhold the mechanism reading — invariant 8's reporting requirement is what catches this. |

## 6 · What this cannot show, at any outcome

**Direction.** A within-cell association between Wnt targets and low
differentiation output is equally consistent with Wnt suppressing the programme
and with less-differentiated cells having more Wnt tone for other reasons. This
is a correlation inside a cross-section; nothing here is causal.

**Survivorship.** Unchanged and unaddressed, as everywhere else in this project.

**Anything about carcinoma.** Adenoma only. The premise does not hold in
carcinoma (`docs/HANDOFF.md` §6g) and this test is not gated on the premise, but
that is not licence to run it where the arms are not comparable.

## 7 · Standing

**Exploratory-but-pre-registered**, and it is the first mechanism test in the
project. It follows a positive result rather than rescuing a negative one, which
is the right order. It needs the team because it introduces a new signature to
the analysis — even though invariant 8 fixed that signature in week 0, which is
most of why this is cheap to propose now.

---

## RESULT — run 2026-09-06. Nothing above was edited.

`results/2026-09-06_0d73b33/`. 43 patients, both arms, `lineage`.

### The gate passed: the Wnt score is not the maturity score

Mean partial r(Wnt, maturity) over 86 patient-arms = **−0.060**, against §5's
0.70 ceiling. **SEPARABLE.** Invariant 8's required report is satisfied and the
mechanism reading is not withheld on those grounds.

### The answer is §5's second branch: TECHNICAL. No mechanism claim follows.

Partial ρ, conditioned on maturity and log depth:

| gene | role | polyp | normal |
|---|---|---|---|
| ACTB | control | −0.032 | −0.062 |
| KRT8 | control | −0.049 | −0.038 |
| EPCAM | epithelial | −0.041 | −0.027 |
| CDX2 | identity | +0.007 | +0.005 |
| MS4A12 | identity | −0.073 | −0.050 |
| **GUCA2A** | **target** | **−0.038** | **−0.045** |

**GUCA2A sits inside the technical floor.** The controls span −0.049 to −0.032
in the polyp arm and GUCA2A is −0.038 — not beyond them, *between* them. And it
is the same value in the **normal** arm (−0.045), where Wnt activity should be
low and unstructured. §5 named that pattern in advance: everything weakly
negative including the controls, in both arms, is a residual gradient the
conditioners did not remove.

**`excludes_zero` is true for almost every row and means almost nothing here.**
Those intervals exclude zero because 43 patients with hundreds of cells each
estimate a tiny correlation precisely, not because the correlation is large.
|ρ| ≈ 0.04 is about 0.16% of variance. Precision without magnitude.

### The over-conditioning objection, tested and dead

**A real design risk, and it is not in §3 because I did not think of it before
the run:** conditioning on maturity is correct if maturity is a *confounder*,
but if the causal chain is Wnt → less maturity → less GUCA2A then maturity is a
**mediator**, and conditioning on it blocks exactly the path being looked for. A
null would then be an artefact of the design rather than a fact about the data.

The job computed `unconditioned_rho` alongside, which settles it:

| gene | unconditioned | conditioned | shift |
|---|---|---|---|
| **GUCA2A** | **−0.037** | **−0.038** | **−0.001** |
| ACTB | −0.016 | −0.032 | −0.016 |
| KRT8 | −0.043 | −0.049 | −0.007 |
| MS4A12 | −0.065 | −0.073 | −0.008 |
| CDX2 | +0.015 | +0.007 | −0.008 |

**Every shift is at most 0.016 and GUCA2A's is 0.001. There was no association
to block.** The unconditioned correlations are already at the floor, and
unconditioned GUCA2A (−0.037) sits inside the unconditioned control range
(−0.016 to −0.043) exactly as the conditioned one does. **The null is not an
over-conditioning artefact**, and it survives an objection the design did not
anticipate.

### What this does and does not say

**Does:** within the surviving mature cells of adenoma, per-cell Wnt-target tone
does **not** track per-cell differentiation output beyond a technical floor.
That is a real negative and it is §5's second row, taken.

**Does not:** that Wnt is uninvolved. This is a *within-cell, within-patient*
correlation among cells already labelled mature. A mechanism operating at the
clone, crypt or tissue level — where whole lesions differ in Wnt tone and the
mature cells inside them differ correspondingly — produces exactly this null,
because the between-patient and between-lesion variation is what a within-patient
correlation removes by construction. **This test was never able to see that**,
and saying so is not a retreat: §6 already said direction was out of reach, and
this is the same limit one level up.

**The controls earned their place.** Without ACTB and KRT8 scored through the
identical path, GUCA2A at −0.038 with an interval excluding zero would have read
as a finding. The floor is the only reason it reads as nothing.

### Standing

**A clean negative on the first mechanism test this project has run**, on a
statistic whose confound was measured rather than assumed, robust to an
objection raised after the design was fixed. It removes one candidate mechanism
at the per-cell level and leaves the tissue-level version untested.
