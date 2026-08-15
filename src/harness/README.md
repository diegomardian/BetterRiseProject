# W2 — Method

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
entire gate. Put it in [`docs/`](../../docs/) and link it from the PR.

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
