# Amendment 1 to the pre-registered MMR contrast — the tumour-arm definition

**Written:** 2026-08-24 · **Author:** W1 (Bode) · **Status:** proposed, needs the
team · **Amends:** the MMRd-vs-MMRp contrast (execution_plan.md §8.4) and
[`docs/prereg_g2_mlh1.md`](prereg_g2_mlh1.md)

> **No expression has been examined.** This amendment responds to a measured
> property of the *method*, discovered while validating malignancy calling, and
> is written before any decomposition has run. The same change made after seeing
> an MMR result would be indefensible and — from the outside — indistinguishable
> from this one. **Its credibility is entirely in the timing**, which is why it
> is a dated document rather than an edit to the original.

---

## 1 · What was measured

62 patients through inferCNV. Restricted to the 36 with a matched normal:

| stratum | separable | median enrichment | malignant fraction of tumour epithelium* |
|---|---|---|---|
| MMR-proficient | **15 / 15** | 7.38 | **24.3%** |
| MMRd (all strata) | 15 / 20 | ~2.5 | **7.2%** |

\* among separable patients only — i.e. where the caller worked in both arms.

**Two distinct biases, and the second was not anticipated.**

1. **Patient level.** Every MMRp patient yields a separable aneuploid
   population; a quarter of MMRd patients do not. Probability of all five
   failures falling outside MMRp by chance: **0.048** one-sided.
2. **Cell level.** *Even among patients where the caller works*, MMRd tumours
   have 3.4x fewer epithelial cells called malignant, tracking their lower
   aneuploidy exactly.

Both follow from one fact: **MMR-deficient tumours are near-diploid.** Copy-number
inference has less to find in them. This is a property of MSI-H biology, not a
defect in inferCNV, and no alternative caller escapes it — CopyKAT infers copy
number from expression too.

## 2 · Why this breaks the contrast as originally specified

The decomposition depends on which cells constitute the tumour arm. Under
malignancy filtering that arm is defined **differently in each stratum**: all of
MMRp survives, three-quarters of MMRd does, and within survivors MMRd
contributes a third as many cells. The comparison is then not MMRd against
MMRp. It is *all MMRp tumours* against *the aneuploid subset of MMRd tumours,
thinned further at the cell level*.

**A consistent bias cancels in a difference. A differential one does not.**

## 3 · What was considered and rejected

| | approach | why not |
|---|---|---|
| a | Filter both arms, as originally specified | Confounded at two levels, in a known direction, along the axis being tested. |
| b | Restrict to separable patients in both arms | Does not help: selecting separable MMRd patients selects the *aneuploid* ones, while selecting separable MMRp patients selects everyone. Same bias, smaller n. |
| c | **Use the unfiltered (sample-of-origin) arm for the MMR contrast only** | Consistent definition across strata, but unbiased *only if* true tumour purity does not differ by MMR — and **the only measurement of purity available is the biased caller**. Picking this asserts something unverifiable. |
| d | Match strata on aneuploidy | Distributions barely overlap (medians 7.4 against 2.5); matching discards most of the cohort. |
| e | Drop the contrast | Discards a pre-registered analysis over a problem that can be measured instead. |

## 4 · The amendment

**Report the MMR contrast under BOTH tumour-arm definitions, and pre-commit to
what agreement and disagreement mean.**

- **Definition A — filtered.** Tumour arm = epithelial cells called malignant.
  Patients with no separable population are `not_estimable`.
- **Definition B — unfiltered.** Tumour arm = all epithelial cells from tumour
  samples. Every matched patient contributes.

**The pre-committed reading:**

- **They agree in direction, and in magnitude to within the interval** → the MMR
  conclusion is robust to the definition, and is reported with both numbers.
- **They disagree** → the MMR contrast is **not identifiable from this data**,
  and is reported as such. Not as a caveat on the filtered number, and not by
  choosing whichever is more interesting.

This is the project's own three-way framing — compositional / intrinsic / not
estimable — applied one level up, to a modelling choice rather than to cell
counts. It is the same move already made for the depth target in decision #14,
where an estimate that changes sign across the nuisance parameter is reported as
unidentified rather than quoted at one arbitrary setting.

**Neither definition is the "real" one.** A is precise and differentially
biased; B is consistent and imprecise. Their agreement is evidence *because*
they fail in unrelated ways.

## 5 · What this does not fix

**Tumour purity may genuinely differ by MMR status**, and definition B would
inherit that. It cannot be checked here: the only purity measurement available
in this cohort is the caller this amendment exists because of. If W3's TCGA
purity estimates can speak to whether MSI-H tumours differ systematically in
epithelial purity, that is worth knowing — but it is a different cohort and
would be supporting evidence, not a fix.

**Stated plainly: definition B is unbiased with respect to the artifact we
measured, and untested with respect to a biological difference we cannot
measure.**

## 6 · Scope

This amendment applies **only to comparisons between MMR strata**. The primary
analysis — tier separation across all matched patients — is not a between-strata
comparison, so the differential loss does not confound it, and it uses
definition A with `not_called` patients reported as `not_estimable`.

Open decision #15 carries the underlying measurement.
