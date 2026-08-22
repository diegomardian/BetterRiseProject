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

## W1-A · The S matrices — *the blocking item*

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

**What W2 needs from you:** nothing immediately. Say whether you want the harness
check above before or after your depth-target rerun.

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
