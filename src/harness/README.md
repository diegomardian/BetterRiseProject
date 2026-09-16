# W2 — Method

> **Taking over W2?** Start with
> [docs/W2_HANDOFF_TO_NEXT_AGENT.md](../../docs/W2_HANDOFF_TO_NEXT_AGENT.md) —
> current state, what is left, and the mistakes worth not repeating.

**Owner:** strongest ML/stats person · **Env:** `env/w2_harness.yml` → `conda activate brp-w2`
**Branch prefix:** `w2/…` · **Blocked by:** W1's five-patient pilot, end of week 2

The most underrated stream. **The only place true ground truth exists**, and the
stream that adjudicates the week-5 gate.

## What you deliver

| Wk | Task | Done when |
|----|------|-----------|
| 1 | Harness **design spec**, reviewed by all four | Written before code. Defines what "known ground truth" means here. |
| 2–3 | Pseudobulk generator: held-out patients, known fractions, known per-cell shifts | Can generate arbitrary (composition, intrinsic) pairs — **including near-zero mature-cell edge cases** |
| 3–4 | Deconvolution bake-off: ν-SVR, CIBERSORTx, MuSiC, BayesPrism, NNLS baseline | Ranked on fraction recovery (r, RMSE) against harness ground truth |
| 4 | 11-gene vs. full-signature comparison, quantified | Settles the Executive Brief claim empirically rather than by argument |
| 4–5 | **Attenuation curve**: intrinsic recovery from pseudobulk vs. mature-cell fraction | The §2.2 calibration curve — a publishable object on its own |
| 5 | Calibrated positivity cutpoints replacing the provisional ones | Derived from where CI width crosses a stated threshold, **not chosen** |

You also own the negative controls: housekeeping genes (should show neither
term) and within-patient label permutation (should destroy both).

## Week 1 is a document, not code

The design spec comes first and all four people review it. It defines what
"known ground truth" means here — which is the load-bearing assumption of the
entire gate. It is at
[docs/harness_design_spec.md](../../docs/harness_design_spec.md) and is **draft
for review**: §7 lists what still needs W1/W3/W4 sign-off, and §4 pre-registers
the three numbers the week-5 cutpoints are derived from.

## What is built

| File | State |
|---|---|
| [truth.py](truth.py) | **done** — analytic Kitagawa terms, parametric and realised truth, identity assertion |
| [pseudobulk.py](pseudobulk.py) | **done** — patient holdout, composition draw, multiplicative shift on integer counts |
| [controls.py](controls.py) | **done** — within-patient permutation, housekeeping negatives |
| [positivity.py](positivity.py) | **done** — two cutpoints (intrinsic provisional, compositional pre-committed as #22), plus G4 with a Wilson interval and a `resolvable` flag |
| [results.py](results.py) | **done** — seven harness table shapes, written via `src.common.io` |
| [deconvolve/](deconvolve/) | protocol, NNLS baseline and ν-SVR **done**; MuSiC wk 3–4; CIBERSORTx, BayesPrism staged |
| [attenuation.py](attenuation.py) | **done** — the §2.2 sweep, oracle + bulk arms |
| [bulk_recovery.py](bulk_recovery.py) | **done** — the thing invariant 6 forbids, measured rather than used. `reference_profiles` now refuses a call that does not say which of the two matrices it is building |
| [calibration.py](calibration.py) | **done** — cutpoints derived from pre-registered criteria; needs CIs attached |
| [bakeoff.py](bakeoff.py) | **done** — ranking, plus the signature-width comparison |
| [interval.py](interval.py) | **done** — per-patient CI; read its docstring on invariant 5 |
| [g1_amendment.py](g1_amendment.py) | **done** — W2's ratification of prereg amendment 2 (#37): five worlds, can the proposed G1 fail? |
| [gate_cost.py](gate_cost.py) | **done** — the gate re-costed at 10 and 36 patients rather than the plan's 60 |
| [ambient_sensitivity.py](ambient_sensitivity.py) | **done** — what residual soup does to a known decomposition; the artefact is one-directional |
| [depth_confound.py](depth_confound.py) | **done** — is the maturity call measuring maturity or sequencing depth? On Lee it is depth (#44) |
| `ingest.py` | unblocked once W4 merges `w2/lee-raw-counts` (open_decisions #8) |

## The two arms, and why both

```python
from src.harness.attenuation import SweepConfig, SweepGrid, run_sweep, summarise_sweep
sweep = run_sweep(SweepConfig(counts, cell_type, patient_id, genes, "GUCA2A"),
                  SweepGrid(), seed=20260815)
summarise_sweep(sweep)
```

**oracle** runs W4's `decompose()` on cell-level summary statistics — the
reliable half, and G3. **bulk** deconvolves, backs the mature mean out of bulk,
then decomposes — the half invariant 6 forbids using for results, run only to
measure how far it can be pushed. That measurement *is* §2.2.

Reporting bulk without oracle beside it confounds estimator error with bulk
attenuation, and those have opposite consequences: fix the estimator, or don't
use bulk for this.

Read `null_arm_recovers_zero()` before any other number in a sweep. It checks
the **parametric** null is exactly zero. The **realised** null is not zero and
must not be asserted to be — normal and tumour are different draws of cells, so
their empirical means differ even when nothing was silenced. Use
`null_arm_noise_ratio()` for that, which is the meaningful version of the
question.

`bulk_overconfidence()` counts rows where bulk returned a number and the truth
was "not estimable". On the first synthetic run that was 41 of 80. See
[the design spec §7](../../docs/harness_design_spec.md) for the validation run
and — importantly — what it does not show.

```python
from src.harness import generate_pseudobulk, patient_holdout

train, held = patient_holdout(patient_ids, n_held_out=5, seed=1)
sample = generate_pseudobulk(
    counts, cell_type, patient_id, genes,
    composition_normal={"mature_colonocyte": 0.4, ...},
    composition_tumour={"mature_colonocyte": 0.02, ...},   # sweep this to zero
    shift={"GUCA2A": 0.5},                                  # 1.0 is the exact null
    held_out_patients=held, n_cells=5000, seed=1,
)
sample.truth.parametric["GUCA2A"]["normal"]["intrinsic"]
sample.truth.realised["GUCA2A"]["normal"]["intrinsic"]
```

**Both truths are recorded on every sample.** Parametric is what we asked for;
realised is what the drawn cells actually have. Recovery against realised
isolates estimator bias, recovery against parametric also carries sampling
noise. Report one and a sampling artefact reads as estimator bias.

## The cutpoints are yours, and W4 imports them

[positivity.py](positivity.py) holds the provisional values (n≥50 / 20≤n<50 /
n<20). W4 calls `classify_estimability()` rather than reimplementing the rule,
so when you recalibrate at week 5 you edit `CUTPOINTS` and W4 inherits it.
Keep `PROVISIONAL` around so the gate memo can show both.

```python
from src.harness import classify_estimability, gate_g4_verdict
classify_estimability(3)          # 'not_estimable'  -> intrinsic MUST be None
gate_g4_verdict([120, 80, 4, 2])  # the G4 numbers, with the pre-committed consequence
```

### There are two cutpoints, on two different counts

Decision #22 gives the **compositional** arm its own gate, on `n_cells_resolved`
— the denominator `mature_fraction` is computed on — at the same 50/20. Without
it a fraction computed on 9% of the epithelium reads identically to one computed
on 90% (issue #36).

```python
from src.harness import estimability_verdicts, classify_counts_frame
estimability_verdicts(n_cells_mature=3, n_cells_resolved=800)
# intrinsic 'not_estimable', compositional 'ok' — two claims, not one
classify_counts_frame(w1_mature_cell_counts)   # both columns, per row
```

**Do not fold them into one verdict by taking the worse.** "The fraction is
imprecise" and "there are too few mature cells to ask about expression" are
different findings and the second is the contribution. `src/schema.py` is frozen
and its `estimability` enum can only carry the intrinsic one, so the
compositional verdict lives in harness tables and the gate memo until the gate
decides otherwise.

Mature cells are a subset of resolved cells, so the compositional gate only ever
binds where the intrinsic gate has **already** flagged the row — 0 rows against
108 on W1's full counts table. The caveat published with the pre-commitment said
the opposite; decision #22 carries the correction and the numbers.

## `reference_profiles` will not guess which matrix you want

The counts frames in this project **carry the panel on purpose** —
`LeeCohort.raw_counts` has to, because the generator applies a multiplicative
shift to the target gene. That made an invariant-2 violation one forgotten
keyword away, silently, back when `exclude_genes` defaulted to `()`.

There is no safe default, so there is no default:

```python
reference_profiles(counts, cell_type, genes, exclude_genes=[target])  # fractions
reference_profiles(counts, cell_type, genes, include_targets=True)    # target profile
reference_profiles(counts, cell_type, genes)                          # ReferenceSeamError
```

An empty `exclude_genes` is refused too, for the reason `build_signature`
refuses an empty target set: it excludes nothing while reading as enforcement.

## What you adjudicate at the gate

| Gate | Your role |
|------|-----------|
| **G1** ambient correction did not eliminate the intrinsic signal | Cross-check W1's retention-vs-abundance statistic |
| **G2** control tiers separate on pilot data | If MLH1 specifically fails: harness passes → detection floor in the data; harness fails → broken estimator. **Both reportable.** |
| **G3** estimator recovers known ground truth on pseudobulk | **Entirely yours.** |
| **G4** <50% of patients below the positivity threshold | Your cutpoints decide this |

## Design notes

- Hold out **patients**, not cells. The patient is the unit of inference
  (CLAUDE.md invariant 5) and it applies to the harness too.
- The near-zero mature-cell cases are the point, not an afterthought. A
  generator that only produces comfortable mixtures validates only the easy half
  — and the third segment is the contribution.
- Include a passthrough case: zero shift must return zero intrinsic, not noise.
- The attenuation bias is **directional** — it shrinks the intrinsic term and
  leaves the compositional term intact, pushing toward "compositional," which is
  also our prior hypothesis. A result that confirms your expectation for
  methodological reasons is the worst kind of result. Measure it; do not argue
  about it.
