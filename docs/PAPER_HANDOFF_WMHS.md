# Paper handoff — WMHS @ NeurIPS 2026, deadline 1 Sept AoE

**Written 2026-08-29 by the W2 agent.** You are taking over to produce a
submission in **~3 days**. Read this, then [CLAUDE.md](../CLAUDE.md), then
[gate_memo_w2.md](gate_memo_w2.md) — the memo is the results inventory and it
contains retractions you must not quote past.

---

## 1 · The honest verdict up front

**Do we have enough results? Yes — for a measurement-validity / abstention paper.
No — for any paper making a biological claim, and that is settled rather than a
matter of more work.**

G2 (the control-tier criterion) **fails on both cohorts, on every tier**. The
negative control MS4A12 is indistinguishable from the positive control tier A
(−0.662 vs −0.653 to −0.90 relative intrinsic loss). Per `execution_plan.md` §5 a
G2 failure means *methods and validation paper, no biological claim*. **Do not
write a paper that says differentiation-marker loss in CRC is intrinsic or
compositional.** The controls do not license it and a reviewer who reads the
appendix will find that.

What we do have is a strong, unusual, honest methods contribution. See §3.

## 2 · The venue, and how well we actually fit

**WMHS @ NeurIPS 2026** — *"World Models for High-Stakes Health: Reliable Clinical
Trial Simulation and Intervention-Aware Reasoning."*
<https://wmhs-neurips.github.io/WMHS/>

| | |
|---|---|
| Deadline | **1 September 2026, AoE** |
| Format | Extended abstract **≤4 pages** main text, or full paper ≤9 |
| References / appendices | **Do not count** toward the limit |
| Archival | **Non-archival.** Under review elsewhere is fine; can be published elsewhere after |
| Review | **Double-blind.** Anonymise names, affiliations, acknowledgments; cite own prior work in third person |
| Mandatory | **A responsible-use statement covering limitations and impacts. Missing it warrants desk rejection.** |
| Submit | OpenReview: `https://openreview.net/group?id=NeurIPS.cc/2026/Workshop/WMHS` |

### Be clear-eyed: we fit ONE topic bullet, not the workshop's core

The workshop's centre of gravity is clinical trial simulation, patient world
models, EHR foundation models, counterfactual treatment effects, target trial
emulation. **We are none of those.** We are single-cell and bulk transcriptomics
in colorectal cancer.

We fit here, verbatim from their topic list:

> *"Uncertainty quantification, calibration, **abstention**, ambiguity, and
> **selective prediction**"*

and secondarily:

> *"Benchmarking, robustness, and validation against real-world evidence …"*
> *"Causal representation learning, causal inference … for intervention-aware
> simulation"*

**Recommended framing, and this is the whole pitch.** Lead with abstention as a
first-class output in high-stakes health measurement; make the single-cell
decomposition the *instance*, not the subject.

The bridge to their audience that actually works: **the compositional-vs-intrinsic
confound is case-mix confounding.** Our estimator is the Kitagawa /
demographic-standardisation decomposition, the same tool epidemiology uses to
separate "the population changed" from "the rate changed" — which is exactly what
synthetic control arms and target trial emulation must handle. A world model of
tissue that reports cell-type-resolved expression change inherits this
identifiability problem whether or not it names it. Say that in the intro and the
paper becomes legible to that room.

**Do not oversell the bridge.** One honest paragraph is better than pretending to
be a trial-simulation paper.

### Scope recommendation: 4-page extended abstract

Three days, double-blind, mandatory responsible-use statement, and a domain
off the workshop's centre. A tight 4 pages that nails abstention will do better
than a sprawling 9 that reads like an oncology methods paper. It is non-archival,
so the full version can go to a methods venue later.

## 3 · The contribution, in the order it should be argued

1. **A decomposition that returns `not estimable` as a first-class third
   outcome.** Every comparable method returns a number regardless. The rule is
   pre-registered (`n_cells_mature ≥ 50` ok, 20–49 wide, <20 abstain) and lives
   in one place both the estimator and the harness import.
2. **The abstention rule is calibrated on simulation and validated on real
   cells.** Chosen on synthetic data in week 5; on real held-out patients,
   recovery of the intrinsic term is **0.85–1.08 across the whole `ok` band** and
   degrades to 1.18 in `wide_interval`. Below the band the estimator does not
   fail gracefully — it **over**-estimates, which is the dangerous direction.
   This is a calibration curve for an abstention rule, measured, not asserted.
3. **Abstention is the modal outcome where the question is sharpest.** At the
   finest granularity rung, **28/28 patients abstain** on the primary cohort, on
   both axes and under both tumour-arm definitions. The finest rung is where the
   compositional question is best posed and it is exactly where nothing is
   estimable.
4. **Identifiability depends on n, and we measured how.** At n=10 the
   compositional band contains zero; at n=27 it excludes zero
   ([−8.05, −2.44]). Same pipeline, same rungs. The abstention machinery told us
   what n was needed *before* the data arrived, and at that n the term appears.
5. **Diagnostics for when a measurement is not separable from a technical
   confound** — depth confounding, with a two-condition test (does the call track
   depth *within* an arm; are the arms matched), an arm-matching correction, and
   a **prevalence ceiling** (`max |ρ| = sqrt(3p(1−p))`) that makes a naive
   correlation test unable to fire below p = 1.37%.
6. **A pre-registered criterion that provably fails when the hypothesis is
   TRUE** (G1), withdrawn rather than repaired, with the reasoning recorded. This
   is a genuinely interesting point about pre-registration in high-stakes
   measurement and it is a *strength*, not an embarrassment.

## 4 · Numbers you may quote — and the ones you may NOT

### ⛔ RETRACTED. Do not quote these. Gate memo §17 is the retraction record.

- **§12's original G4 table.** It silently used one of two tumour-arm
  definitions (`.first()` took `filtered` in 224/224 combinations by row order).
  Superseded by `results/2026-08-28_0be41e6/g4_verdict_gse178341_by_arm.parquet`.
- **"Tier A's loss is predominantly intrinsic"** as a finding about tier A. It is
  an arithmetic identity in the `m_t → 0` limit, numerically identical for tier
  D, and it reverses at `best4` — the rung tier A is defined on.
- **"best4 is the cleanest rung."** Its |ρ| was at the prevalence ceiling; it is
  marginally the *most* affected rung once normalised.
- **"a defensible PASS at n=28 allows up to 12."** `largest_clean_pass(28) = 8`.
- **§14's "+25.1% gap"** as stated. It is cell-pooled and violates invariant 5.
  Per patient it is +0.170, 95% CI **[−0.011, +0.342]**, containing zero.

### ✅ Current and safe

| claim | number | source |
|---|---|---|
| Abstention at the finest rung, primary cohort | 28/28 both arms, CI [0.879, 1.000] | `2026-08-28_0be41e6/g4_verdict_gse178341_by_arm.parquet` |
| Coarse rungs, pre-registered verdict | `not_identifiable` (arms disagree) | `..._amendment1.parquet` |
| Recovery inside the `ok` band | 0.85–1.08; 1.18 in `wide_interval` | memo §15 |
| Compositional band, primary, n=27 | GUCA2A [−8.05, −2.44] excl. 0 | `2026-08-29_08a6aa0/..._bands.parquet` |
| Intrinsic band, primary | GUCA2A [−27.61, −14.45] excl. 0 | same |
| G2 failure, replicated | MS4A12 −0.662 (Pelka) vs −0.662 (Lee) | memo §18.3 |
| Ambient artefact is one-directional | manufactures intrinsic ~4.6% of compositional; **never** compositional | `2026-08-26_09f0bc3/ambient_sensitivity_sweep.parquet` |
| Depth confound, before/after one labeller | ρ −0.92 → −0.31; gap 46.1% → 25.1% | memo §13.1 |
| Gap survives arm matching | +20.4% (mine), +0.182 (independent matcher) | memo §14.2, §17.3 |
| Prevalence ceiling | `sqrt(3p(1−p))`; 0.20 unreachable below p=1.37% | memo §17.4 |

## 5 · MUST DO before writing: reconcile two decompositions

**There are two independent GSE178341 decompositions and they do not use the same
population.**

| | patients | arms | GUCA2A comp / intr |
|---|---|---|---|
| W2 (`2026-08-29_08a6aa0`) | 28 (27 unfiltered) | **both, carried separately** | −2.53 / −13.15 |
| W1 (`2026-08-28_8965a6f`) | 32 | **no `tumour_arm` column** | −2.80 / −14.79 |

They agree on sign and magnitude, which is real independent verification of the
core result and worth a sentence in the paper. **But W1's table has the defect
§17.1 fixed** — no arm column at all, so either one arm was used and not
recorded, or they were collapsed. Prereg amendment 1 forbids choosing between
them silently.

**Resolve this first.** Either establish that W1's used one arm and say which, or
regenerate from `w1/decomposition-summary` with the arms carried. Do not publish
a number until the two agree on a population.

## 6 · Remaining runs, with costs

| task | cost | needed for the paper? |
|---|---|---|
| Reconcile the two decompositions (§5) | ~30 min | **YES — blocking** |
| Depth-matched decomposition on the primary cohort | ~20 min | **Strongly recommended.** `lineage` carries the depth caveat unmatched; on Lee matching moved the compositional term 30% and left the intrinsic one unchanged |
| Figures (2) | 2–3 h | yes |
| Responsible-use statement | 30 min | **mandatory — desk reject without it** |
| Anonymisation pass | 30 min | **mandatory** |

Everything else is done. There is **no** further data to collect and no further
method to build.

Suggested figures, both from existing parquet:
1. **The calibration curve** — recovery vs mature-cell count, with the
   pre-registered cutpoint marked and the `wide`/`abstain` bands shaded. This is
   the paper's core claim in one panel.
2. **Abstention rate by granularity rung and cohort size** — shows abstention
   rising to 100% at the finest rung, and the n=10 vs n=27 resolution difference.

## 7 · Hard constraints you must not break

- **`None` is not `0.0`** (CLAUDE.md invariant 1). An abstained estimate is
  `None`. The whole paper rests on this; do not let a figure impute zeros.
- **Report both tumour arms.** Prereg amendment 1: disagreement means *not
  identifiable*, not a choice.
- **Estimate per study, meta-analyse, never pool** (invariant 4). Lee and Pelka
  stay separate; the replication is the point.
- **Bootstrap over patients, not cells** (invariant 5). §14's error was exactly
  this; don't repeat it.
- **The interaction term is reported separately** (invariant 7). It is non-zero
  and its band excludes zero — GUCA2A [+2.03, +7.55].
- **No biological claim.** §1.
- **Anonymise.** The repo is public under the author's name. Do not link it, or
  use an anonymised mirror.

## 8 · Things a reviewer will find, so put them in yourself

Workshops reward candour and this project has an unusual amount of it. Consider a
short "what we got wrong" paragraph — it is genuinely part of the contribution:

- Five guards that could not fail were shipped and withdrawn (namespace-blind
  leakage checks, an acceptance test comparing two disjoint identifier spaces).
- A published result table was computed over the wrong population and caught by
  an adversarial audit, not by review.
- A pre-registered gate criterion was withdrawn because it fails when the
  hypothesis is true.
- The maturity labels on one cohort manufactured a 46-point effect out of
  sequencing depth, in the hypothesised direction, before it was caught.

That list is the strongest available evidence that the abstention machinery is
load-bearing rather than decorative.

## 9 · Where things are

- Results inventory + all retractions: [gate_memo_w2.md](gate_memo_w2.md),
  read §11–§18 and **§17 especially**
- Decisions: [open_decisions.md](open_decisions.md), #19–#24
- G1's withdrawal: [g1_withdrawal_case.md](g1_withdrawal_case.md)
- Abstention rule: `src/harness/positivity.py`
- Confound diagnostics: `src/harness/depth_confound.py`
- Gate operating characteristics: `src/harness/gate_cost.py`
- Tests are the specification: 1,139 passing; `tests/test_positivity.py`,
  `test_depth_confound.py`, `test_gate_cost.py`
- No paper infrastructure exists yet — no `.tex`, no `.bib`, no figure scripts.
