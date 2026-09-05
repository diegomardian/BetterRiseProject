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

*(appended after the run, as `prereg_g2_mlh1.md` does.)*
