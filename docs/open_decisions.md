# Open decisions

Things the scaffold could not decide on its own, recorded so they get decided
deliberately rather than discovered in week 9. Close each with a line in the
weekly meeting and a PR that updates this file.

Open questions of a scientific kind live in
[execution_plan.md §10](../execution_plan.md#10-open-questions). This file is
for decisions that block or shape code.

---

## 1 · MUC2 and TFF3 are both targets and labels — OPEN

**Raised:** scaffold, 2026-08-15 · **Owner:** W1 + W2 · **Needed by:** week 3
(W1 builds labels in weeks 3–4)

MUC2 and TFF3 sit in **tier E** (targets) and in **labelling axis 2**
(opposite lineage: MUC2, TFF3, SPDEF, ITLN1). Invariant 2 says target genes
never appear in labels or the reference matrix. Read across the whole panel,
axis 2 loses half its markers and the "agreement across structurally different
axes" argument weakens considerably.

The scaffold takes the **narrow** reading, without deciding the question:
`build_signature()` excludes the target set *for the run in question*, not the
whole panel. Consequence as implemented — **a run testing MUC2 or TFF3 may not
use axis 2.** `tests/test_freeze.py` pins the collision set at exactly
`{MUC2, TFF3}` so a future panel or axis edit that adds a new one fails loudly.

Options:

| | Approach | Cost |
|---|---|---|
| a | Keep the narrow reading. MUC2/TFF3 results carry axes 1, 3 only. | Two tier-E genes get thinner axis coverage. Cheapest. |
| b | Drop MUC2 and TFF3 from tier E. | Panel edit — 2 approvals. Tier E is exploratory, so the loss is small. |
| c | Rebuild axis 2 on SPDEF, ITLN1 and non-panel goblet markers. | Axis edit — 2 approvals, and axes were frozen first on purpose. |

**Recommendation: (a) now, revisit at the gate.** It costs nothing today and
the gate will show whether axis 2 is carrying weight.

---

## 2 · Who owns the shared gene index — OPEN

**Raised:** scaffold, 2026-08-15 · **Owner:** W1 + W3 · **Needed by:** week 1

W1 emits S matrices on a fixed gene index; W3 emits bulk on the same index;
integration is a join. The plan does not say who *produces* the index, and both
need it in week 1. See [config/gene_index/README.md](../config/gene_index/README.md).

**Recommendation: W1 emits it in week 1** from the GSE178341 feature table,
committed as `config/gene_index/gene_index_1.0.0.txt`. W3 reindexes onto it
rather than the reverse — single-cell features are the more constrained set. If
W1's ingest slips, W3 emits a provisional index from the GDC gene model and W1
conforms. Decide in the week-1 meeting; do not let both build one.

---

## 3 · Symbol vs. Ensembl ID on the shared index — OPEN

**Raised:** scaffold, 2026-08-15 · **Owner:** W1 + W3 · **Needed by:** week 1

Follows from #2 and is the usual way a join silently loses 8% of genes. TCGA
STAR counts arrive as Ensembl IDs with versions; the panel and both labelling
axes are written as symbols; single-cell references are usually symbols.

**Recommendation: Ensembl IDs (unversioned) as the index, symbols as a mapped
column.** Panel and axis lookups resolve through the map, which is committed
alongside the index. Unmapped panel genes are a week-1 finding, not a week-9
surprise.

---

## 4 · COAD and READ: pool or stratify — OPEN

**Raised:** execution_plan.md §4 (W3 gotchas) · **Owner:** W3 · **Needed by:** week 3

Different treatment patterns. The plan says decide now and write down the
reason. It is here so the reason gets written somewhere findable.

---

## 5 · CODEOWNERS handles and branch protection — OPEN

**Raised:** scaffold, 2026-08-15 · **Owner:** W2 · **Needed by:** week 1

[`.github/CODEOWNERS`](../.github/CODEOWNERS) has placeholder handles, so
GitHub currently ignores most of it and the "two approvals for frozen code" rule
is convention rather than mechanism. Fill in the handles and turn on branch
protection for `main`.

---

## 6 · Which cohort W2 ingests — OVERTAKEN BY EVENTS, see §8

**Raised:** W2 planning, 2026-08-15 · **Resolved:** 2026-08-15 by W4 shipping
`src/estimator/lee_io.py`, which loads **both** cohorts (SMC and KUL3) behind one
`load_lee_cohort(which=...)`. The split proposed below did not happen and no
longer needs to — nobody duplicated anything, because W4 got there first.

What replaces it is a narrower question about the *shape* of the artifact rather
than who downloads it. See **§8** below.

<details>
<summary>Original proposal, kept for the record</summary>

W2 needs real cells in week 1 so the harness is never queued behind W1's
five-patient pilot. W4 owns both Lee cohorts (GSE132465 SMC, GSE144735 KUL3) and
does its own ingest and ambient correction in weeks 1–2. If W2 ingests Lee too,
two people run the same QC in the same fortnight.

**Proposal: split the pair rather than duplicate it.**

| | Cohort | Does |
|---|---|---|
| W4 | GSE132465 (SMC) | ingest, QC, ambient correction, labels — as planned |
| W2 | GSE144735 (KUL3) | same pipeline *shape*, then hands the cleaned object to W4 |

W4 needs both cohorts eventually and gets one for free. W2 gets real cells in
week 1 without a redundant fortnight. Pipeline shape is coordinated; code is not
shared yet — the same rule W4 already has with W1.

**Fallback if W4 objects:** W2 waits for W4's SMC object and spends week 1 on the
design spec and the data-free plumbing (`common/io.py`, `harness/results.py`,
the deconvolution adapter protocol, the NNLS baseline) — all of which are
already done, so the cost of the fallback is low.

**Note on independence:** sharing a cohort between W2 and W4 does not compromise
the harness. The harness needs cells with labels, not a dataset nobody else has
touched. What it does require is that the patients held out for a pseudobulk
sample are held out — which is enforced in `generate_pseudobulk`, not assumed.

</details>

---

## 7 · Detectable effect for cutpoint calibration — OPEN

**Raised:** W2 design spec, 2026-08-15 · **Owner:** W2 + W4 · **Needed by:** week 4,
before the attenuation sweep runs

The cutpoints are derived from the smallest mature-cell count at which the
estimator can still discriminate a pre-registered effect. The spec pre-registers
`s = 0.5` — a halving of per-cell output — because it sits mid-band in the
published attenuation range (×0.6–0.8).

If W4 expects MLH1 silencing to be closer to complete (`s` near 0), a smaller
effect is easier to detect, the cutpoints get looser and G4 gets easier to pass.
If the realistic effect is milder, the reverse. **This number must be fixed
before the sweep is looked at**, which is what makes the cutpoints derived rather
than chosen. See [harness_design_spec.md §4](harness_design_spec.md).

---

## 8 · The harness needs raw counts; `lee_io` emits CP10K — OPEN

**Raised:** W2, 2026-08-15, reviewing `w4/estimator-core` · **Owner:** W2 + W4 ·
**Needed by:** week 2, before the harness runs on real cells

`LeeCohort.expression` is **CP10K-normalised** and restricted to
genes-of-interest. That is exactly right for `decompose_cohort`, which consumes
per-patient summary means and needs a linear depth-normalised scale for the
Kitagawa identity to be additive. It is the wrong artifact for the harness,
which needs:

| | `lee_io` gives | harness needs |
|---|---|---|
| scale | CP10K floats | **raw integer counts** — binomial thinning is defined on counts |
| gene width | genes-of-interest | **the full index** — deconvolution wants 500–2000 markers |

Handing the CP10K frame to `generate_pseudobulk` used to truncate silently:
`astype(int64)` turns every value below 1.0 into 0, destroying exactly the
low-expressing cells the near-zero mature-cell edge cases are made of. **The
generator now refuses a non-integer matrix** with a message naming CP10K, so this
can no longer happen quietly — but refusing is not the same as being unblocked.

Options:

| | Approach | Cost |
|---|---|---|
| a | `lee_io` grows a `return_raw_counts=True` / second accessor on `LeeCohort` | Small, and W4 already streams the raw matrix — the counts exist before the CP10K step at `lee_io.py:322` |
| b | W2 writes its own reader against the same manifest files | Duplicates the parsing and the QC that `ingest.qc_flags` already does |
| c | `lee_io` caches a raw-count parquet next to the CP10K frame | Both arms served, one parse, costs disk |

**Recommendation: (a).** It is a keyword argument on code W4 has already
written, and the raw counts are in hand at the point of normalisation. Needs
W4's agreement since it is their module.

**Until it lands**, the harness runs on synthetic fixtures. That is not a block
on the week-5 gate — G3 is a statement about the estimator's algebra and
synthetic cells test it fine — but the attenuation curve wants real
expression distributions, so this should not drift past week 3.

---

## 9 · `doubly_robust` folds the interaction into both arms — OPEN

**Raised:** W2, 2026-08-15, reviewing `d9c08c0` · **Owner:** W4 + W2 ·
**Needed by:** before any decomposition result is written with
`weighting="doubly_robust"`

CLAUDE.md invariant 7: *the interaction term is reported separately, never
folded into either arm.* The pooled-reference split in
`kitagawa._doubly_robust_split` reports `interaction = 0.0` — but the cross term
has not vanished, it has been distributed **50/50 into the other two arms**:

```
comp_pooled = Δf·(m_n+m_t)/2 = Δf·m_n + Δf·Δm/2 = comp_normal + interaction/2
intr_pooled = (f_n+f_t)/2·Δm = f_n·Δm + Δf·Δm/2 = intr_normal + interaction/2
```

Verified numerically on `f_n=0.40, f_t=0.10, m_n=10.0, m_t=4.0`:

| weighting | compositional | intrinsic | interaction |
|---|---|---|---|
| normal | −3.0000 | −2.4000 | **+1.8000** |
| tumour | −1.2000 | −0.6000 | −1.8000 |
| doubly_robust | −2.1000 | −1.5000 | **0.0000** |

`−2.1 = −3.0 + 0.9` and `−1.5 = −2.4 + 0.9`, and `0.9` is exactly half the
interaction term. Here that shrinks the intrinsic term by **37.5%, toward
zero** — the direction of our prior hypothesis, which README calls out as the
worst way for a result to move.

This is not a claim that W4 is wrong. Kline (2011) is a real citation, the
pooled reference is a real estimator, and W4's docstring is explicit that it is
a first cut pending cell-level AIPW. The problem is narrower:

1. It is what invariant 7 was written to forbid, and invariants change by PR
   with two approvals and a written reason — not by implementation.
2. Writing `0.0` into the schema's `interaction` column reads as "no interaction
   here" when the honest statement is "interaction distributed into the arms."
   That is the invariant-1 failure mode wearing a different hat.

Options:

| | Approach | Cost |
|---|---|---|
| a | Amend invariant 7 to "never folded silently" and require the pooled split to record what it absorbed | Frozen-doc PR, 2 approvals. Keeps a citable estimator. |
| b | Keep `doubly_robust` but report the normal-weighting interaction alongside, not `0.0` | Schema already has the column; the value stops being a lie |
| c | Rename it — it is a pooled-reference split, not AIPW — and defer true doubly-robust to cell-level data | Cheapest, most honest about what it is; W4's docstring already says this |

**Recommendation: (c) plus (b).** Call it what it is, and put the real cross
term in the interaction column rather than zero. Neither needs new maths.

---

## 10 · Which interval goes in the schema's `ci_low`/`ci_high` — W2 PROPOSES, W4 TO CONFIRM

**Raised:** W4 in `bootstrap_over_patients`' docstring · **Answered by:** W2,
2026-08-16 · **Owner:** W2 + W4 · **Needed by:** week 5

W4 wrote, and was right to:

> *"which quantity belongs there (this cohort-level bootstrap band on
> `intrinsic`, a per-patient interval from within-patient cell-count
> uncertainty, or a value read off the hierarchical model instead) is an open
> call for whoever owns the week 4–5 hierarchical-model deliverable to make
> deliberately."*

**The two are different estimands and the project needs both.**

| Question | Sampling unit | Where |
|---|---|---|
| What is the cohort's intrinsic loss, and how sure are we? | **patients** (invariant 5) | W4's `bootstrap_over_patients` |
| Does *this* patient have enough mature cells for their own estimate to mean anything? | **cells**, by construction | W2's `harness/interval.py` |

The cutpoints have to be calibrated on the second. A cohort-level band
broadcast onto every patient row is *identical* for a patient with 800 mature
cells and one with 21, so coverage against `n_cells_mature` is flat and no
cutpoint exists. Measured on the new interval, CI width runs 2.583 at 5 mature
cells to 0.596 at 800 — a 4.3× range, which is the dependence a cutpoint needs.

**This is not a violation of invariant 5.** Invariant 5 exists because
resampling cells to make a *population* claim inflates n by the cells-per-patient
count. That failure mode requires the claim to be about patients. A
within-patient statement has cells as its sample by construction. The argument
is in `src/harness/interval.py`'s module docstring in full, so a reader who finds
the resampling loop finds the reasoning next to it.

**Proposal for the schema slot:** keep W4's current choice — the cohort-level
intrinsic band — because the per-patient row is read as part of a population
result, and say which one it is wherever it is presented, as
`attach_intrinsic_ci` already instructs. The within-patient interval is a
harness-side quantity for calibration and does not need a schema slot.

**W4: flag it if you disagree.** The only thing that would change is which
number lands in `ci_low`/`ci_high`; both intervals exist and are tested either
way.

---

## Closed

*(none yet — move entries here with the date and the decision, do not delete them)*
