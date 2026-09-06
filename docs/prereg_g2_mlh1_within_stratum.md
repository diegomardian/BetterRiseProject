# Pre-registration — MLH1 as a positive control for the instrument

**Written:** 2026-09-06 · **Author:** W1 (Bode) · **Status:** proposed, awaiting
team ratification · **Supersedes the design of** [prereg_g2_mlh1.md](prereg_g2_mlh1.md),
whose difference-in-differences is not available on this data (§2).

> **This document is only worth what its timestamp is worth.** It states a
> directional prediction, the analysis that tests it, the interval that analysis
> will report, and what would falsify it — **before any MLH1 expression has been
> read on this pipeline**. The calibration and power tables it cites were
> committed first, in `cf3c76f`, for the same reason. Ratify it, amend it, or
> reject it; do not edit the prediction after results exist.

---

## 1 · The question, which is about the instrument and not about cancer

Every null this project has produced rests on a detection statistic whose
sensitivity to real transcriptional silencing has never been shown:

| result | verdict |
|---|---|
| coexpression premise, 3 original cohorts | UNRESOLVED |
| coexpression premise, 13 ICBI studies, 122 patients | UNRESOLVED |
| bulk CIMP specificity screen | NOT SPECIFIC |
| adenoma panel, gene-specificity of the fall | not gene-specific |

**A negative result from an instrument of unknown sensitivity is not evidence of
absence. It is not evidence of anything.** That is the gap this closes.

MLH1 is the one gene in this cohort whose silencing is known **from an assay
rather than from expression**. Pelka's atlas rows carry
`MLH1_promoter_methylation_status` for all 62 patients, patient-level, from
methylation calling; promoter hypermethylation silencing MLH1 transcription is
established biology and not a hypothesis of this project.

**Prediction.** Among patients whose MLH1 promoters are methylated, MLH1
detection **falls** in the mature cells of the tumour arm relative to the mature
cells of their own normal arm — a negative mean `delta_cloglog`.

**Direction is the whole test.** A rise falsifies it (§5).

---

## 2 · What changed since [prereg_g2_mlh1.md](prereg_g2_mlh1.md), and why

That document pre-registered a **difference-in-differences**: intrinsic MLH1
loss in `mlh1_methylated` against `mlh1_intact_mmrd`, patients who reach the
same MSI-H phenotype through MSH2/MSH6/PMS2 with MLH1 transcription untouched.
It was the stronger design and it is **not available on this data.**

Only **29 of Pelka's 62 patients** survive the ICBI pipeline's own filters
(naive fraction, ≥100 epithelial cells per arm, QC, labellability, ≥40 mature
cells, depth matching). Joined to the pre-registered strata at `lineage`:

| stratum | patients scored | median mature cells per arm |
|---|---|---|
| `mlh1_methylated` | **10** | 262 |
| `mlh1_intact_mmrd` | **4** | 127 |
| `mmr_proficient` | 14 | 182 |
| `mlh1_deficient_unmethylated` | 1 | 230 |

**Four is the same number the original prereg reached after depth matching**, by
a completely different route — the atlas reprocessing rather than the GSE178341
path. Two independent pipelines agreeing on four says it is a property of the
cohort, not of anyone's filters. The mechanistic negative control does not
exist at usable size on either route.

**And stratifying cannot fix the dilution, because the stratum you would
stratify into has four patients in it.**

*A check that had to be run before this could be written.* The atlas's
per-cell annotation and the week-0 clinical strata are independent derivations
of the same fact, and they **agree exactly** on all 62 patients — 22 `meth`
against 22 `mlh1_methylated`, no crossings in either direction. A disagreement
would have meant the arm this reading is about is not the arm the
pre-registration named. `strata_for()` re-runs that comparison in the job and
records crossings in the sidecar.

---

## 3 · The three arms, fixed, with different standing

They are never merged and the table names the standing of each on every row.

| arm | n | standing |
|---|---|---|
| `mlh1_methylated` | 10 | **primary.** Powered — see §6. |
| `mlh1_unmethylated` | 19 | **secondary, CONFOUNDED.** |
| `mlh1_intact_mmrd` | 4 | **pre-registered control, UNDERPOWERED.** |

**The secondary arm's confound, stated before it is read.** It is 14
MMR-proficient patients, 4 intact-MMRd and 1 MLH1-deficient-unmethylated, so it
mixes methylation status with **MSI status**. It is not the mechanistic control
the original prereg wanted and it will not be described as one. It is reported
because at n=19 it can say something about strong silencing (§6) where the n=4
arm cannot say anything at all, and because a confounded comparison named as
confounded is worth more than an uninterpretable one named as a control.

**The n=4 arm carries no verdict, in either direction.** At n=4 this project's
own percentile interval excludes zero **18.8%** of the time under a true null
(§4). Its numbers are still emitted — suppressing them would make the table's
shape depend on the result, and a reader could not then check the claim that the
arm is uninformative.

---

## 4 · The interval, chosen by measurement and fixed here

**This reading does not use the interval the rest of the repository uses.**
Justified before the fact by `results/2026-09-06_7685baa/interval_calibration.parquet`.

The percentile bootstrap over patients — used by `premise_holds`, `summarise`,
`specificity` and `control_log2_interval` — is **narrower than it claims at
small n**, and the mechanism is arithmetic with no data in it. The bootstrap
distribution of a mean has standard deviation `s·sqrt((n−1)/n)`, so the
percentile interval is approximately `mean ± z·s·sqrt((n−1)/n)/√n` where the
correct interval is `mean ± t(n−1)·s/√n`. The ratio is `z·sqrt((n−1)/n)/t(n−1)`:
**a function of n alone.**

| n | width vs correct | false-positive rate, closed form | measured |
|---|---|---|---|
| 4 | 0.53× | **18.8%** | 19.1% |
| 10 | 0.82× | **9.6%** | 10.7% |
| 19 | 0.91× | 7.3% | 8.3% |
| 20 | 0.91× | 7.1% | 7.3% |
| 44 | 0.96× | 5.9% | 6.2% |

Nominal is 5%. The closed form is a **floor** — it assumes the per-patient
values are normal enough for the bootstrap mean to be; the measured excess above
it, positive at every n, is the skew of a rare transcript's delta. Two
independent routes to one number.

**Neither sophisticated repair works.** BCa is *worse* — 17 of 20 cells
miscalibrated against the plain percentile's 14 — because its bias correction
and jackknife acceleration are themselves estimated from the same ten numbers.
More bootstrap replicates do not help: `N_BOOTSTRAP` is already 10,000 and the
error is not Monte Carlo error.

**So this reading reports the Student-t interval**, which was calibrated in
**0 of 20** cells miscalibrated, worst 5.9%, at both a rare gene's abundance and
a common one's.

**This generalises past MLH1 and is not a claim about this gene.** Every
percentile-bootstrap number in this repository carries this at whatever n its
table ran at. What it does to existing results is in
[HANDOFF.md §6d](HANDOFF.md); the committed numbers are not being restated.

---

## 5 · What would falsify it, and what each branch commits to

Fixed here, taken by `instrument_verdict()` against a table rather than by a
reader against a paragraph.

| branch | verdict | pre-committed consequence |
|---|---|---|
| MLH1 falls, interval excludes zero | **INSTRUMENT SEES KNOWN SILENCING** | The statistic has demonstrated sensitivity to a silencing event known from an assay. The project's nulls become evidence of absence rather than absence of evidence — **without** a mechanistic control arm behind them. |
| MLH1 **rises**, interval excludes zero | **WRONG DIRECTION — falsified** | An instrument that fires the wrong way on a known event is not calibrated. The nulls are not rehabilitated. |
| interval contains zero | **INSTRUMENT DOES NOT SEE IT** | Informative against **strong** silencing, **not** informative against moderate — §6. The other nulls stay uninformative: an instrument that cannot see a known event cannot be cited for not seeing an unknown one. |
| premise fails in the methylated arm | **UNINTERPRETABLE** | Not a negative result. The cells compared are not established to be the same population. |
| fewer than 5 of 10 patients carry any MLH1+ cell in the normal arm | **NOT ESTIMABLE** | The delta would be mostly the boundary correction. Under **invariant 1** this is not zero silencing and must never be written as `0.0`. |

**The premise is checked in the methylated arm only.** A premise holding over
all 29 patients does not license a claim about these 10.

*And that check is anti-conservative here, which is recorded now rather than
discovered later.* `premise_holds` is an **equivalence** test — it asks whether
the interval fits *inside* a tolerance. A narrower interval fits more easily, so
at n=10 the percentile bootstrap pushes that check **toward HOLDS**. A HOLDS at
this n is weaker than the same word at n=44. It is left on the percentile
bootstrap so this reading's premise verdict stays comparable to every other
premise verdict in the project; the caveat travels in the sidecar.

---

## 6 · Power, on the interval that will actually be reported

From `results/2026-09-06_7685baa/mlh1_power.parquet`. Student-t, n=10, the real
per-patient cell counts and depths, MLH1 at 0.039 CP10K → ~3.2% detection →
about **8 positive cells per patient per arm**.

| silencing | τ = 0 | τ = 0.2 | τ = 0.4 |
|---|---|---|---|
| none (H₀) | 4.6% | 5.4% | 4.7% |
| **50%** | 75.9% | **73.7%** | 65.4% |
| **75%** | 99.3% | **99.3%** | 98.6% |
| 90% | 100% | 99.9% | 99.9% |

**τ is between-patient heterogeneity and it is measured, not assumed.**
`results/2026-09-06_7685baa/interval_heterogeneity.parquet`, on Pelka, net of
binomial sampling: ACTB 0.203, KRT8 0.197, EPCAM 0.289, CDX2 0.496 — and τ
**rises as baseline detection falls**. MLH1 at ~0.03 detection sits off the
bottom of that range, so 0.4 is the honest pessimistic end and **τ = 0 is not a
row anyone may quote**.

**Two things this table corrects in the analysis that proposed this reading.**

1. **The earlier figure of 86% at 50% silencing came from the percentile
   bootstrap** — the interval this same design rejects for firing at 9.6% under
   the null. Power and coverage are properties of one method and may not be
   taken from two. On the interval actually being reported it is **73.7%**.
2. **That simulation had τ = 0 in it.** Only binomial noise, so it could not
   come out underpowered from patient-to-patient variation, because there was
   none. That is a check that cannot fail, in a power calculation.

`power_curve()` now has no code path that emits power without the false-positive
rate of its own method, measured on the same generator, on the same row.

**So what this reading is powered for, stated in advance:** it will detect
silencing of **75% or more** (≥98.6% at every τ) and is **marginal at 50%**
(65–76%). A null here is evidence against strong silencing and **is not**
evidence against moderate silencing, and §5's third branch must be read with
that sentence attached.

### What the two-sample contrasts could do, for the record

Neither is the primary reading. Measured on the same generator, Welch t-test:

| contrast | 50% | 75% | 90% |
|---|---|---|---|
| vs `mlh1_intact_mmrd` (n=4) — **the pre-registered DiD** | 22% | 55% | 82% |
| vs `mlh1_unmethylated` (n=19) — confounded | 59% | 97% | ~100% |

**This is the quantitative form of "the DiD is not available."** At 75%
silencing the pre-registered control arm is a coin flip. The confounded arm is
not — which is why §3 reports it as secondary rather than discarding it.

---

## 7 · What this cannot say, at any outcome

**That the response is specific to promoter methylation.** That needs the
mechanistic negative control, and §2 is the finding that this cohort does not
have one. A fall in the methylated arm shows the statistic responds to a gene
whose promoter is methylated; it does not show the response is *because of* the
methylation. `mlh1_unmethylated` bears on it and is confounded.

**That MLH1 silencing was observed in a form the decomposition could use.** This
is the coexpression statistic inside a fixed label, not the Kitagawa split. The
decomposition's algebraic collapse (HANDOFF §2) is untouched by anything here.

**Anything about survivorship.** GUCA2A-high or MLH1-high cells having been
preferentially destroyed is not transcript-detectable and no result here bears
on it.

---

## 8 · Standing

**This is a test of the instrument, and it is the first one this project has
had.** It is not a G2 gate criterion and does not revive G2 — G2 failed as
pre-registered and [prereg_g2_mlh1.md](prereg_g2_mlh1.md)'s closing sections
stand unchanged.

Its value is asymmetric and both directions are useful:

- **If MLH1 falls:** every null in this project becomes a null *from an
  instrument with demonstrated sensitivity*, which is a materially stronger
  claim than any of them currently carries.
- **If it does not:** the nulls stay uninformative, and that reframes the
  project — the write-up would then be about a method whose sensitivity could
  not be established, which is a different and more honest paper than one
  reporting four negative results.

**Because it can reframe the project, it needs the team, not just W1.**

---

## RESULT

*Not run. Requires the 30 GB ICBI atlas, which is cluster-only:*

    qsub src/reference/jobs/mlh1_positive_control.sh

*Nothing above may be edited when it is.*
