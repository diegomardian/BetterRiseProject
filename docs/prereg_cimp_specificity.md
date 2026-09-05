# Is the CIMP association with GUCA2A locus-specific, or axis-level?

**Written:** 2026-09-04, before the test is run.
**Arm:** W3 (bulk). `src/bulk/` is W3's under CONTRIBUTING §2 — this needs sign-off
or a W3 owner, not a silent edit.
**Status:** screening. Nothing here can establish a mechanism; see §6.

## 1. What is already known, and why this is not a blind pre-registration

**The stratum medians are already committed and already read.** They are in
`gse39582_premise_bimodality.parquet`, they were printed to a log on 2026-09-02,
and the direction is not in doubt:

| gene | CIMP+ (n=91) | CIMP− (n=405) | gap | fold |
|---|---|---|---|---|
| CDX2   | 6.549 | 7.517 | +0.968 | 1.96× |
| GUCA2A | 5.869 | 6.523 | +0.654 | 1.57× |
| MS4A12 | 4.001 | 4.532 | +0.531 | 1.44× |

So this document pre-registers a **test and a decision rule**, not ignorance of
the direction. Saying otherwise would be theatre. What is genuinely unknown is
whether the two contrasts below have intervals that exclude zero — which is the
only thing the run adds, and the only reason to do it.

## 2. The question

GUCA2A is lower in the methylator phenotype. Two readings:

- **Locus-specific silencing.** CIMP methylation reaches the GUCA2A promoter.
  Predicts GUCA2A falls further than markers of the same population.
- **Axis-level suppression.** CIMP+ tumours carry less differentiated
  epithelium, or suppress the CDX2→GUCA2A axis upstream. Predicts GUCA2A falls
  *with* its population's other markers, not beyond them.

Only the first would justify the 450k methylation leg. That is what this
screening is for: closing a route cheaply, or not.

## 3. Two references, on purpose, and they are expected to disagree

**CDX2** — upstream transcription factor, drives GUCA2A. The mechanistically
relevant reference. From the medians, GUCA2A falls **0.314 less** than CDX2.

**MS4A12** — colonocyte-restricted, no regulatory relationship to GUCA2A. The
population reference. From the medians, GUCA2A falls **0.123 more**.

**The references already point opposite ways, and that is the expected finding.**
This is the panel's tier structure doing the job it was frozen in week 0 to do.
A single pooled interaction would hide it.

## 4. Estimands

Per sample, the within-sample difference — pairing the two genes inside one
array cancels per-sample loading and platform effects, and is lower variance
than subtracting two independently estimated stratum medians:

    D_ref(s) = log2 GUCA2A(s) − log2 ref(s)

    interaction(ref) = mean[D_ref | CIMP+] − mean[D_ref | CIMP−]

for `ref ∈ {CDX2, MS4A12}`. Stratified bootstrap over samples — resample the 91
and the 405 separately — 10,000 draws, percentile interval, seed from
`src.common.provenance.DEFAULT_SEED`.

The unit of inference is the sample. GSE39582 is one sample per patient, so
invariant 5 is satisfied rather than bypassed.

**Second read, same pass:** the same contrasts adjusted for MMR status and
tumour location, both already in the parsed metadata (dMMR n=75, proximal
n=224). Categorical adjustment, no new compute.

## 5. Decision rule, committed before the run

**GUCA2A-specific silencing is supported only if BOTH interactions favour
GUCA2A with intervals excluding zero.** One reference agreeing is not support —
that is the configuration the medians already show, and calling it support would
be choosing the reference after seeing the answer.

Anything else reads as: *no robust GUCA2A-locus-specific CIMP signal; the
methylator association is axis-level or compositional.*

**Expected outcome, stated in advance:** not specific. CDX2 interaction near
−0.31, MS4A12 near +0.12, and the second of those may well contain zero — in
which case even the weak arm evaporates.

**Escalation, committed in advance:** if both interactions favour GUCA2A with
intervals excluding zero, that is the surprise, and it triggers ESTIMATE-based
purity adjustment on GPL570 (legitimate — ESTIMATE was built for Affymetrix) and
a full re-read before anything is claimed. Purity does **not** enter otherwise;
paying for it up front to test a hypothesis the medians already argue against
is not a good use of it.

## 5a. Correction to §5, 2026-09-04, before any run

**§5's two expected values are sign-flipped.** It says *"CDX2 interaction near
−0.31, MS4A12 near +0.12"*. Under §4's own estimand they are the other way
round:

| reference | D given CIMP+ | D given CIMP− | interaction |
|---|---|---|---|
| CDX2   | 5.869 − 6.549 = −0.680 | 6.523 − 7.517 = −0.994 | **+0.314** |
| MS4A12 | 5.869 − 4.001 = +1.868 | 6.523 − 4.532 = +1.991 | **−0.123** |

§3 and §4 are correct, and so is the implementation — `run_cimp_specificity`
computes exactly §4, and its tests use +0.31 and −0.12. The error is confined to
§5's prose, and §5 contradicts itself on it: *"the second of those may well
contain zero — even the weak arm evaporates"* only parses if MS4A12 is the
weakly **favourable** arm, which requires it to be negative.

**Corrected expectation, unchanged in substance:** CDX2 near **+0.314** (GUCA2A
falls *less* than CDX2 — against specificity), MS4A12 near **−0.123** (GUCA2A
falls *more* — weakly for it), and the MS4A12 interval may well contain zero, in
which case even that arm evaporates. One reference favourable out of two, so the
§5 decision rule returns **NOT SPECIFIC** either way. The verdict is untouched.

Recorded as a correction rather than an edit to §5. The value of writing this
before the run is entirely in not being able to revise it afterwards, and a
silent fix — even to a typo — spends that. Flagged now because the run will
return +0.31/−0.12 and would otherwise look like the pre-committed expectation
had flipped, which is the confusion this document exists to prevent.

## 6. What this cannot establish, in any outcome

Bulk expression is fraction × per-cell mean. "Lower in CIMP+" is equally
consistent with methylation silencing, with CIMP+ carrying fewer colonocytes,
and with destruction having preferentially removed GUCA2A-high cells. Adjusting
for other markers controls the general effect; **it cannot control
colonocyte-specific compositional loss.** This is the same non-identifiability
the single-cell arm reported on 2026-09-04 across three cohorts, one level up.

So a null here closes a route. A positive here would not open one — it would
promote the question to the methylation leg, where bulk 450k has the identical
two-factor structure and only cell-type-resolved methylation escapes it.

**And note the asymmetry in the likely result:** CDX2 drives GUCA2A, so "GUCA2A
falls less than CDX2" is not a clean falsification of silencing. It is the
parsimonious reading that CIMP+ suppresses the axis upstream. What the screening
closes is *locus-specific promoter silencing as the bulk-level story* — which is
precisely the claim that would have justified the 450k work.

## 7. Result

Run 2026-09-05 on GSE39582, 3 genes x 585 samples, sha `9203809`, clean tree.
Table: `results/2026-09-05_9203809/gse39582_cimp_specificity.parquet`.

**Verdict: NOT SPECIFIC.** Zero of two references favour GUCA2A.

| reference | contrast | 95% CI | clears zero |
|---|---|---|---|
| CDX2   | **+0.544** | [+0.219, +0.878] | yes |
| MS4A12 | −0.147 | [−0.514, +0.205] | **no** |

The §5 decision rule needed both references favourable with intervals excluding
zero. One argues against locus-specificity and the other is silent.

**The outcome is cleaner than the one predicted.** §5a expected the two
references to disagree — CDX2 against, MS4A12 weakly for. MS4A12's interval
contains zero, so the weak arm does not merely fail to reach significance, it
carries no information at all. The reading is not "the references disagree" but
"one reference argues against specificity and the other says nothing."

**A discrepancy worth recording rather than smoothing.** §5a pre-committed CDX2
"near +0.314". The run returns **+0.544** — same sign, about 1.7x the magnitude.
This is not an error: +0.314 was computed from stratum *medians*, and the
estimand is the *mean of within-sample paired differences*. Different statistics
on a skewed distribution, and the pairing changes the variance structure as
well. The direction and the verdict are unaffected. It is recorded because a
pre-registration that only reports the predictions it got right is not doing
anything.

**Robustness across the adjusted reads.** The CDX2 contrast is positive in all
four MMR and location strata and clears zero in three:

| reference | stratum | contrast | 95% CI | clears zero | n (CIMP+/−) |
|---|---|---|---|---|---|
| CDX2 | dMMR | +0.976 | [+0.218, +1.716] | yes | 42/27 |
| CDX2 | pMMR | +0.487 | [+0.045, +0.947] | yes | 41/339 |
| CDX2 | proximal | +0.812 | [+0.392, +1.224] | yes | 74/130 |
| CDX2 | distal | +0.327 | [−0.258, +0.911] | no | 17/275 |
| MS4A12 | *(all four)* | −0.445 … +0.207 | all straddle 0 | **0 of 4** | |

The single non-clearing CDX2 stratum is distal, where only 17 samples are CIMP+.
Same direction, too few to resolve — underpowered rather than contradictory.
MS4A12 clears zero in none of the four, so its silence is not an artefact of the
pooled analysis.

**Two caveats on the strata.** They are not independent: CIMP+ tumours are
predominantly proximal (74 of 91) and predominantly dMMR, so these are
overlapping views of largely the same samples rather than four tests. And the
adjustment is stratification, not regression, so it controls the marginal
association and not the joint one.

**What this closes.** Locus-specific promoter silencing of GUCA2A is not the
bulk-level story in this cohort. GUCA2A falls in the methylator phenotype, but
it falls *less* than the transcription factor upstream of it and no differently
from the other colonocyte-restricted marker. That is what upstream suppression
of the CDX2 axis, or a shift in differentiated content, predicts.

So the 450k methylation leg has no hypothesis left to test at bulk level, which
was the only reason to consider it. Per §6, that leg would in any case have
carried the same fraction-times-mean ambiguity one level down.

**What this does not establish.** Everything in §6 stands. No outcome here
separates silencing from colonocyte-specific compositional loss, and this is a
null that closes a route rather than a finding about mechanism. It is consistent
with the single-cell arm's 2026-09-04 result across three cohorts — that
silencing versus destruction is not identifiable on this data — and adds a
population-scale instance of the same limit rather than resolving it.
