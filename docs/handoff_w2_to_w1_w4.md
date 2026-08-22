# Handoff — what W2 needs from W1 and W4

**From:** W2 (Method / harness) · **Date:** 2026-08-22
**State of W2:** weeks 3–5 complete on synthetic cells. The harness, the §2.2
attenuation sweep, the per-patient interval, cutpoint calibration, the bake-off
and the negative controls all run and are tested. What is missing is real cells
and two decisions.

Each item says what is needed, how W2 will check it, and what W2 does with it.
The acceptance criteria are guesses about your side, not requirements on how you
get there — if one is wrong or costlier than it looks, say so rather than
building the wrong thing.

**Correction up front.** An earlier draft of this said the pilot had not been
run, on the evidence that nothing in `results/` was W1's. That was wrong:
`src/reference/jobs/run_pilot.py` exists, ten pilot arms have been run, and the
kappa is measured. Apologies — the ask below is narrower as a result.

---

# W1 — one artifact, one team decision

## W1-A · The S matrices — **DELIVERED AND ACCEPTED** (PR #17, 2026-08-22)

All four passed every published check. Nothing to redo on the matrices themselves.

| rung | genes x types | 500–2000 | no target leak | non-epithelial cols | on shared index |
|---|---|---|---|---|---|
| epithelial | 800 x 5 | PASS | PASS | PASS | 800/800 |
| lineage | 800 x 6 | PASS | PASS | PASS | 800/800 |
| crypt_position | 800 x 6 | PASS | PASS | PASS | 800/800 |
| best4 | 800 x 6 | PASS | PASS | PASS | 800/800 |

Two deviations from the spec, both fine and worth recording:

- Versioned `0.1.0-pilot`, not `1.0.0`, and written to `results/{date}_{sha}/`
  rather than `s_matrix_path()`'s `data/processed/reference/`. **Better than what
  was asked** — it carries provenance. W2 reads them by glob; if you want
  `s_matrix_path()` to find them later, that helper needs a pointer, which is
  W2's to change.
- Index is `gene_index_0.9.0`, not a `1.0.0`. Fine for a pilot; open decisions
  #2/#3 still need settling before a full-cohort version.

### But the rung *structure* has two problems, and one is a real bug

Found while checking the matrices. Neither is a defect in the S matrices; both
are in how the rungs are defined, and both hit the four-resolution curve.

**1 · The epithelial rung's compositional term is structurally zero.**

`mature_fraction` is 1.000000 for every patient, both arms, both axes, so
Δ(mature fraction) = 0.000000 exactly:

```
opposite_lineage/epithelial   max |delta mature fraction| = 0.000000
       stem_pole/epithelial   max |delta mature fraction| = 0.000000
```

At this rung "mature" is *all resolved epithelium*, and the denominator is also
resolved epithelium — so the fraction is 1 by construction and the compositional
arm cannot move. The rung contributes a guaranteed zero to the curve.

This looks like a denominator choice rather than a deep problem. If the epithelial
rung's denominator were **all cells** rather than resolved epithelial cells, it
would measure epithelial content of the sample, which genuinely differs between
normal and tumour and is a meaningful compositional quantity. **W1: worth a look
— it is a small change with a large effect on what the coarsest rung means.**

**2 · The degeneracy is axis-specific, which is better news than reported.**

`lineage` and `crypt_position` are identical on `stem_pole` — confirmed
independently: their S-matrix columns are bit-identical,
`corr(differentiated, crypt_top) = 1.0000`, `np.allclose` True. But on
`opposite_lineage` they separate:

| axis | lineage | crypt_position | separate? |
|---|---|---|---|
| stem_pole | 0.5996 | 0.5996 | **no — identical** |
| opposite_lineage | 0.3711 | 0.5326 | **yes** |

So the curve is not lost, it is **axis-dependent**: two usable interior points on
axis 2, one on axis 1. Combined with the harness result below — identical
partitions give bit-identical estimates, genuinely different ones separate at
33% — the picture is complete, and axis 2 is where the granularity contribution
lives.

**3 · 37% of epithelial cells are unresolved** (median `unresolved_fraction`
0.367, up to 0.65 in one arm). Not a blocker, but it is a third of the compartment
excluded before any biology, and it should be stated wherever the mature fraction
is reported.

---

## W1-A (original ask, kept for the record)

`build_signature()` works, `_select_markers` is implemented, the pilot has run.
What does not exist yet is the emitted artifact W2 and W3 both join against.

```python
from src.common.paths import s_matrix_path
from src.common.panel import granularity_rungs, panel_genes
from src.reference.signature import build_signature

for rung in granularity_rungs():        # epithelial, lineage, crypt_position, best4
    S = build_signature(
        expression,                      # cells x genes, log-normalised
        cell_type_at(rung),
        target_genes=panel_genes(),      # the whole panel — invariant 2
        gene_index=shared_gene_index,    # config/gene_index/
        n_genes=1000,                    # must be 500..2000 (§2.1 error #4)
        require_non_epithelial=True,
    )
    S.to_parquet(s_matrix_path(rung, "1.0.0"))
```

**W2's acceptance check** — exactly this, nothing hidden:

| Check | Why |
|---|---|
| File at `s_matrix_path(rung, "1.0.0")` | the frozen path contract |
| Index ⊆ the shared gene index | integration is a join, not a negotiation |
| **No panel gene in the index** | invariant 2 — `build_signature` asserts it; belt and braces |
| Columns include stromal, immune, endothelial | bulk CRC is 30–60% non-epithelial; without them stromal signal is absorbed arbitrarily — the CMS4 failure mode |
| 500 ≤ rows ≤ 2000 | ν-SVR robustness is dimensional |

**Emit what you have, even if it is one rung.** Given the degeneracy below, a
single `lineage` matrix unblocks the bake-off and the real attenuation curve.
Four rungs is the eventual target, not the gate for starting.

**Note for the gene index:** you flagged the deposit is hg19
(`GRCh37_liftover_v28`) while TCGA is GRCh38. That lands on open decisions #2/#3
and it is W3's problem as much as yours — but the S matrix inherits whichever
index you build on, so it needs settling before `1.0.0` means anything.

## W1-B · Rung degeneracy — a team decision W2 has a stake in

You reported that `lineage` and `crypt_position` are **the same partition on
axis 1** in all ten pilot arms, and `best4` is unusable (sens 0.04).

That is not just a W1 problem. **It hits a headline deliverable.** README design
decision 3 and §6.2 promise the split reported "as a curve across four
resolutions", with the divergence between them as a contribution. Two identical
rungs and one unusable one leaves a curve with two real points.

W2's read, offered as input rather than a verdict:

- This is a **finding**, not a failure, and it is publishable as one. "The
  granularity knob does less than the field assumes, and here is the measurement"
  is a real result — it is the same shape as the non-identifiability finding the
  project already treats as a headline-if-true.
- It is also **exactly what the harness is for**. Before deciding the rungs are
  degenerate in the biology, W2 can check whether the *estimator* separates them
  when the truth says it should: generate pseudobulk where lineage and
  crypt_position genuinely differ, and see whether the decomposition recovers the
  difference. If it does, the degeneracy is real biology or real labelling; if it
  does not, it is us. **W2 will run this** — it needs no new data.
- Your "choose the depth target by kappa, not by the flag" reads right to me. If
  the sweep gives a target where the two rungs *do* separate, that is worth
  knowing before anyone concludes they cannot.

### The harness check is done — the estimator is cleared

Run on a synthetic cohort with a genuinely nested structure: a coarse rung that
pools `crypt_top` and `crypt_mid` into "mature", and a fine rung that calls only
`crypt_top` mature. `src/harness/rungs.py`, tested in `tests/test_rungs.py`.

**A · when the rungs genuinely differ:**

| rung | n mature | mean normal | mean tumour | compositional | intrinsic |
|---|---|---|---|---|---|
| lineage (coarse) | 750 | 43.88 | 14.40 | −6.58 | **−11.79** |
| crypt_position (fine) | 150 | 79.77 | 40.44 | −11.97 | **−7.87** |

Relative gap: **33% on the intrinsic term, 45% on the compositional term.**

**B · when the two rungs name the same cell types** — W1's observed case:

```
intrinsic: -11.794400 vs -11.794400    absolute gap 0.00e+00
```

**Conclusion: the estimator is not the explanation.** It resolves a real rung
difference at a third of the effect size, and returns bit-identical answers only
when the partitions are genuinely identical. So the degeneracy W1 measured is a
statement about the labelling, and it can be reported as one.

**A second thing fell out, and it is the more interesting one.** Between the two
rungs the compositional and intrinsic terms move in *opposite* directions —
coarse (−6.58, −11.79), fine (−11.97, −7.87) — while the total is roughly
conserved (−18.4 vs −19.8). The granularity choice does not add or remove loss;
it **reallocates it between the two mechanisms**. That is §6.2's "if it swings,
that divergence is the contribution", demonstrated on a case where the ground
truth is known. It also means a single-rung result would present a modelling
choice as a measurement, exactly as README design decision 3 warns.

**What W2 needs from you:** nothing blocking. Worth knowing whether your
depth-target rerun produces a target where the rungs stop being the same
partition — if it does, this reallocation is measurable on real data and it is a
result rather than a caveat.

---

# W4 — two decisions, no code

Both are questions W4 raised and deliberately left open. Neither needs
implementation; W2 needs a position so the gate memo can state one.

## W4-A · Decision #9 — `doubly_robust` folds the interaction into both arms

CLAUDE.md invariant 7: *the interaction term is reported separately, never folded
into either arm.* The pooled-reference split reports `interaction = 0.0`, but the
cross term has been distributed 50/50 into the other two arms. Verified on
`f_n=0.40, f_t=0.10, m_n=10.0, m_t=4.0`:

| weighting | compositional | intrinsic | interaction |
|---|---|---|---|
| normal | −3.0000 | −2.4000 | **+1.8000** |
| tumour | −1.2000 | −0.6000 | −1.8000 |
| doubly_robust | −2.1000 | −1.5000 | **0.0000** |

`−2.1 = −3.0 + 0.9` and `−1.5 = −2.4 + 0.9`; `0.9` is exactly half the
interaction. Here it shrinks the intrinsic term **37.5% toward zero** — the
direction README names as the worst way for a result to move.

This is not a claim the implementation is wrong. Kline (2011) is real, the pooled
reference is a real estimator, and the docstring is explicit that it is a first
cut. The problem is narrower: invariant 7 was written to forbid exactly this, and
invariants change by PR with a written reason, not by implementation.

**W2's recommendation:** rename it to what it is — a pooled-reference split, not
AIPW — and put the real cross term in the `interaction` column instead of `0.0`.
Neither needs new maths.

**"No, keep it" with a written reason is an equally good outcome.** What W2
cannot do is write a gate memo that leaves it unstated.

## W4-B · Decision #10 — which interval fills the schema's CI slot

`attach_intrinsic_ci` broadcasts a cohort-level band onto every patient row, and
says so. W2 built a *within-patient* interval for cutpoint calibration, because a
cohort band is identical at 800 mature cells and at 21 — the coverage curve is
flat and no cutpoint exists.

**W2's proposal: keep your current choice**, the cohort-level intrinsic band,
because the patient row is read as part of a population result. The within-patient
interval stays harness-side and needs no schema slot.

`src/harness/interval.py`'s docstring argues why resampling cells there is not an
invariant-5 violation: invariant 5 governs *population* inference; this is a
within-patient statement where cells are the sample by construction. Two
estimands, both real, neither replacing the other.

**Confirm or object** — only the number in `ci_low`/`ci_high` changes.

## W4-C · Review a change merged without you

`w2/lee-raw-counts` added `keep_raw_counts` and `extra_genes` to
`load_lee_cohort`, in your module, merged on the repo owner's instruction because
it blocked every harness run on real cells. No existing default changed. It
deserves your eyes retrospectively.

---

# Cross-cutting

## X-1 · `open_decisions.md` has duplicate numbers

There are now **two #9s and two #10s** — W2's (`doubly_robust`, interval) and
W1's (matched normal, MLH1 pre-registration). An earlier collision was resolved
by renumbering W1's #8 to #11; this one has not been. References like "open
decision #9" are ambiguous, and both pairs are cited from code comments.

**Proposal:** renumber the later-dated pair, note it at each head as was done for
#11, and chase the references. W2 will do it if nobody objects — a docs change,
but touching another workstream's entries wants a nod.

## X-2 · The dev-dependency list keeps breaking main

Twice in three days `[project.optional-dependencies]` lost a line to conflict
resolution — `h5py`, then anndata's workaround interacting with a `provenance`
bug. Both fixed; main is green at 822 passed / 4 skipped, verified in a fresh
CI-equivalent venv. [CODEOWNERS](../.github/CODEOWNERS) still has placeholder
handles, so nothing requires a second reviewer on the one file every workstream
edits. Open decision #5, and cheap.

---

# What W2 owes you back

Commitments, not requests.

## To W1, on the MLH1 pre-registration

W2 is a named owner and had not answered. **Position: yes, pre-register it, and
do it before any expression is looked at.** The stratified version is strictly
stronger than "does MLH1 come out intrinsic on average", because it makes a
directional prediction with MMR status held fixed — and a prediction that could
fail is worth more than an average that cannot. Two conditions:

1. The predictions land in the repo with a git sha **before** the strata are
   crossed with expression. That is the entire value.
2. `mlh1_methylated` has 12 matched patients. Against W2's provisional cutpoints
   that supports a per-patient intrinsic estimate but is thin for a between-strata
   contrast — report it as an estimate with an interval, not a test. Same
   treatment §6.2 gives the MMR contrast.

## To W1, on n=36 — and it is worse than 36

Taken, and it moves W2's side of the gate. G4 asks whether fewer than 50% of
patients fall below the positivity threshold; **the denominator is 36, not 62** —
and your "~30 unsorted in both arms" is the number that actually constrains the
compositional arm. The gate memo will use those rather than the larger figure.

This also raises a question W2 should answer and has not: the provisional
cutpoints were calibrated per-patient, but with n≈30 the *patient-level* bootstrap
W4 runs has far less to resample than the plan assumed. W2 will quantify what
n=30 does to the meta-analytic interval and report it at the gate.

## To W4, on the harness

The bake-off is ready to rank your decomposition against known truth the moment
an S matrix exists. Preliminary G3 on synthetic cells passes: the oracle arm
recovers the known split within 1% wherever the mature compartment is non-empty,
and the broken-estimator control fails as it should.

## Timeline

W2 is not idle while waiting. The harness runs on the real Lee cohorts through
`load_lee_cohort(..., keep_raw_counts=True)`, so the attenuation curve will be on
real expression distributions before the S matrices land. What they change is
that G2 and G3 become statements about the **primary** cohort rather than the
replication one.
