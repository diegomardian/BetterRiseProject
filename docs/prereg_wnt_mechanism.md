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

## RESULT

*Not run. Requires a cluster pass over the ICBI atlas — no new data.*
