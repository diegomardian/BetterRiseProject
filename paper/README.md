# ICBINB submission — revision

A revision of the workshop draft, made against reviewer criteria for
ICBINB 2026 (*Failure Modes of AI in Biology*).

## Read this first: these files do not build on their own

The canonical draft **lives outside this repository** — there is no `.tex`,
`.bib` or `.sty` anywhere else in the tree. Missing here and required to
compile: `neurips_2026.sty`, `refs.bib`, `sections/llm.tex`,
`sections/appendix.tex`, and `figures/*.pdf`.

So treat this directory as **a patch against wherever the draft actually
lives**, not as the paper.

> **If the canonical copy is in Overleaf, apply these and close the PR rather
> than merging it.** Two sources of truth for one artifact is the failure this
> project has now hit three times — two gene indexes, two labellers, two lock
> functions — and a paper is a worse place for it than code, because nothing
> here fails a test when the copies diverge.
>
> If instead the team wants the paper under version control alongside the result
> tables it cites, then merge this **and delete the external copy in the same
> sitting.** Either is fine. Having both is not.

## How to apply

| File | Action |
|---|---|
| `main.tex` | replaces the existing `main.tex` |
| `appendix_withdrawn.tex` | append to `sections/appendix.tex` |
| `refs_additions.bib` | append to `refs.bib` — **verify every entry first** |

## What changed, and why

**Anonymised.** The style option is `dblblindworkshop` and the draft printed
four names, an affiliation and placeholder contact fields. That is a desk reject
before any reviewer reads the content, and it was the single highest-priority
fix. `\author{}` now, restore for camera-ready.

**New §3, "Where the failure stops".** The largest change. The draft reported
four failure modes and no case where the machinery works, which leaves a
reviewer unable to separate *"these checks are structurally blind"* from
*"this pipeline never had signal to find."* §3 tests §2's mechanism on the bulk
arm, where the precondition holds and the mechanism does not transfer:

- arms unbalanced in depth as in §2 — 36.5M vs 48.4M median library,
  *p* = 5e-10 — and **the sign is reversed** from the single-cell cohort
- non-detection flat at ≤1.6% per library-size decile, against the
  80.8% → 8.5% monotone gradient the single-cell cohort shows
- library size explains *r²* = 0.046, and runs the wrong way for a dropout
  artefact
- the distributional test the arm's conclusion rests on moves 0.982 → 0.992
  conditioned on library size

This converts "everything we built failed" into "here is where this failure
mode applies and here is where it stops", which is the difference between a
lab notebook and a finding.

**§2/§3 tension resolved.** §2's headline was ρ = +0.33 against a bound of 0.52,
while §4 is a section about how pairing a correlation with the wrong arm's bound
is an error the project made three times. §2 now states that both numbers come
from the depleted arm, and says why the sentence names the arm.

**Withdrawn numbers moved to the appendix.** They were mid-§2, inviting *"if
your most-quoted number was irreproducible, why trust the rest?"* and never
answering it. The appendix version answers it directly, and §8 gains the lesson
as a named takeaway.

**§8 gains a prospective case.** Every failure in the draft was found *after* a
check reported clean. The Stage 4 pre-registration hit the same class of defect
— a statistic confounded with a nuisance variable — and was caught **before the
analysis ran**, because the specification was written down in a form that could
be objected to. That is the only forward-looking example in the paper and §8 is
where it belongs.

**Related work.** Half a paragraph conceding that depth confounding is known,
then stating the delta precisely: not the confound, but that the standard guards
against it are structurally unable to detect it.

**Contributions list** in §1, leading with the bound, which is the most portable
result and previously sat third.

**Abstract** rebuilt around the bound and the boundary case. 231 words; cut the
last two sentences if the page budget bites.

**One fix not on the review list:** `fig:tiers` was defined but never referenced
in the body. Now cited in §5.

## Before submitting — open items

1. **Verify the five new bib entries.** They were written from memory of the
   literature, not from resolved DOIs. `squair2021confronting` especially: it is
   closest to our framing and the one a reviewer is most likely to know. If an
   entry does not say what the paragraph claims, **cut the citation** rather
   than rewording the paper around it.

2. **Confirm §2's ρ = 0.33 and bound = 0.52 come from the same arm.** The
   sentence now asserts they do, because §4 requires it. If they do not, that
   number needs recomputing before the claim stands. This is exactly the pairing
   error §4 documents, so it would be an unfortunate one to ship.

3. **Decide where the paper lives** — see the warning above.

## Provenance of the new numbers in §3

All from the bulk arm's committed tables and reproducible from the repo:

| Claim | Source |
|---|---|
| 624 tumours, 51 normals; median library 36.5M vs 48.4M; *p* = 5e-10 | `data/processed/bulk/tcga_counts_0.9.0.parquet` + `sample_manifest.tsv` |
| non-detection ≤1.6% per decile; *r²* = 0.046 | same, with `tcga_log2cpm_0.9.0.parquet` |
| dip 0.982 → 0.992 | `results/2026-08-17_cc06981/tcga_premise_bimodality.parquet`, residualised |
| 80.8% → 8.5% single-cell gradient | issue #44 |
| tier-D gene lost ~200-fold in bulk | `results/notes/w3.1_gdc_ingest.md`, recorded 2026-08-17 |
| §8 prospective case: R² 0.891 → 0.000; slope 0.742 vs 0.989 | issue #54 and `config/stage4_prespecification.yaml` |
