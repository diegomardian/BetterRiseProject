# The competitor benchmark — what it measures and what it found

**Run 2026-08-29 · seed 20260829 · 6 worlds x 200 replicates · 2,000 cells per arm**
Code: [`competitors.py`](competitors.py), [`bench.py`](bench.py) · Tables: [`results/`](results/)
Reproduce: `python -m submission.run_bench` — bit-identical across processes
(world seeds are CRC32, not `hash()`; see `bench.world_seed`)

> **Everything here is synthetic.** Every number comes from simulated cells whose
> truth is known in closed form. Nothing in this directory is a result about
> colorectal cancer.

---

## The question

`README.md:83-86` says, of every existing method:

> "Every existing method returns a number. None flags that the intrinsic estimate
> is meaningless in a tumour with no mature cells left."

That is the project's whole novelty claim, and it is an **empirical claim about
other methods** that had never been measured. This benchmark measures it.

The headline is therefore **not accuracy**. It is what a method does where the
estimand *does not exist* — a tumour with no mature cells, where "how much does
each surviving mature cell make" has no referent.

## The setup

A mature cell expresses `Poisson(mu)` in normal tissue and `Poisson(mu*s)` in
tumour; an immature cell expresses nothing. The arm mean is then
`f x (mean among mature)`, which is exactly what Kitagawa splits — so
`harness.truth.analytic_terms` gives the truth in closed form rather than by
simulation, and `s = 1` yields an intrinsic term of *exactly* 0.0.

Six worlds. The one that matters is **`annihilated`**: `frac_mature_tumour = 0`,
so no mature tumour cell survives and the intrinsic estimand is undefined.

Five methods run; a sixth is reported as skipped.

## Result 1 — where the estimand does not exist

| method | can refuse | returned a number | **false confidence** | median \|intrinsic\| invented |
|---|---|---|---|---|
| `kitagawa+positivity` (ours) | yes | 0 / 200 | **0.00** | — (refused all 200) |
| `composition-only` | no | 0 / 200 | **0.00** | — (*no intrinsic arm at all*) |
| `kitagawa-no-gate` | no | 200 / 200 | **1.00** | 7.99 |
| `pseudobulk-de` | no | 200 / 200 | **1.00** | 7.99 |
| `naive-delta-mean` | no | 200 / 200 | **1.00** | 20.00 |

**The claim holds, and the invented numbers are large** — on a scale where the
normal-tissue per-cell mean is 20.0, the naive method reports the entire 20-point
drop as silencing in cells that do not exist.

Two readings that the table deliberately keeps apart:

- **`composition-only` scores 0.00 for a different reason.** It has no intrinsic
  arm (`estimates_intrinsic = False`, `n_refused = 0`). That is *inapplicability*,
  not caution, and crediting it as a refusal would reward a method for not
  competing. Milo/scCODA-shaped methods sit here.
- **`pseudobulk-de` and `kitagawa-no-gate` return identical values.** Not a bug:
  at `f_t = 0` both reduce algebraically to `-(f_n x mean_n)`. They diverge
  normally (−3.99 vs −7.58) as soon as any mature tumour cell exists.

## Result 2 — the counterweight

Refusing is trivially achievable by refusing always. So the same methods are
scored where the estimand **does** exist and a real intrinsic effect is present:

| method | detection rate | median abs error | median signed error |
|---|---|---|---|
| `kitagawa-no-gate` | 1.000 | 0.122 | +0.017 |
| `kitagawa+positivity` (ours) | **0.853** | **0.111** | +0.007 |
| `pseudobulk-de` | 1.000 | 3.435 | −3.435 |
| `naive-delta-mean` | 1.000 | 6.329 | −6.329 |
| `composition-only` | 0.000 | — | — |

**The trade is 14.7% of detections for 100% of the false confidence.** The gate
declines to answer in 88 of 600 cases where an answer was available — all of them
in the `depleted_wide` regime, ~20 mature cells, right at the cutpoint.

And the gate does not cost accuracy where it does answer: median absolute error
is *slightly lower* with the gate (0.111 vs 0.122), because the cases it drops
are the noisy low-n ones.

`naive-delta-mean` and `pseudobulk-de` are biased by −6.33 and −3.43 because they
never standardise: they report compositional loss as silencing. That is a
different failure from false confidence and is included so the two are not
conflated.

## Why the ablation is the load-bearing comparison

`kitagawa-no-gate` is identical arithmetic with the estimability gate removed —
it differs from ours in exactly one way. Every other method differs in several at
once, so none of them can isolate what the gate is responsible for. The pair
`kitagawa+positivity` vs `kitagawa-no-gate` is the actual argument; the rest is
context showing the ungated behaviour is what the field does.

## Limitations — stated plainly

1. **Fully synthetic**, and the generative model is *exactly* the model Kitagawa
   assumes. This is a favourable setting for the estimator. It tests the **gate**,
   not robustness to a misspecified model, and it cannot speak to real data.
2. **The competitors are faithful reimplementations, not the published software.**
   `cacoa` ships as an adapter that reports its own unavailability
   (`r-devtools` is pinned; cacoa installs from GitHub and is not installed) rather
   than being silently absent. A reviewer is entitled to discount the comparison
   accordingly.
3. **`naive-delta-mean` and `pseudobulk-de` answer a slightly different question**
   by design, so their bias in Result 2 is partly unfair to them. The
   false-confidence result in Result 1 does not depend on that, and the ablation
   does not depend on them at all.
4. **The cutpoint is not this benchmark's.** `n < 20` comes from
   `harness.positivity.CUTPOINTS`, the project's own pre-committed rule, so our
   method is scored against a threshold chosen before this benchmark existed —
   but it is still *our* threshold, and 0.835 would move if it moved.
5. **One estimand.** Only the intrinsic term's estimability is tested. The
   compositional arm has its own degeneracy problem, which this does not touch.
