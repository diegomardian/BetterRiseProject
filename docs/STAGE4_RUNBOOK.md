# Stage 4 — running the variance arm

Everything below is committed and tested. What is **not** here is the TCGA bulk
matrix, which is gitignored and lives with the data. That is the only reason
this is a cluster run rather than a laptop one.

Read `docs/HANDOFF.md` first if you are arriving cold.

---

## Before you start: what is already decided

| | |
|---|---|
| Pre-specification | `config/stage4_prespecification.yaml`, **locked 2026-09-05** |
| Prediction | GUCA2A percentile ≤ 0.95 **and** CDX2 > 0.95, in both cohorts |
| Statistic | percentile within an abundance-matched null. **Never a raw R²** |
| Instrument gate | non-epithelial fraction vs (1 − ABSOLUTE purity), r ≥ 0.5, **STOP on failure** |
| Pooling | never. Per cohort, per rung (invariant 4) |

The lock is enforced in code: `require_locked_prespec` refuses a proposed spec
and `run_stage4_variance` calls it before anything else.

---

## The one thing that will bite you

**The committed S matrices are on the wrong scale for this.** They are the mean
of `log1p(CP10K)` — correct for marker selection, wrong for a linear mixture.
TCGA TPM is a linear mixture.

Run them as-is against linear bulk and NNLS returns
`mature_colonocyte_fraction` as **exactly 0.0 on every sample**. That column is
the entire predictor. A constant predictor gives every gene an R² of 0, which is
what the pre-registered arm expects for GUCA2A.

**The instrument gate does not catch this.** It reads the non-epithelial
aggregate, which comes back at r = 0.881 in the same run where the predictor is
a constant. Epithelium-vs-everything survives the mismatch; the
epithelial-internal split does not, and that is the only part Stage 4 uses.

So the driver refuses the mismatch up front, and `--linearise-reference` is
required rather than default. `expm1` of a mean-of-log is a **geometric** mean,
biased low by Jensen, worst for the most dispersed genes. It is a repair, not
the fix. The fix is a linearly-built reference from W1 — see the last section.

---

## The short version

Two jobs, in order. From the repo on the cluster:

```
export BRP_DATA_DIR=/projectnb/rise-batteries/bode/guanylin/data
qsub src/bulk/ingest_cluster.sh     # only if tcga_tpm_1.0.0.parquet is absent
qsub src/bulk/stage4_cluster.sh     # after the first one finishes
```

**Check the version before assuming you need the first job.** The matrices on
disk as of 2026-09-05 were `*_0.9.0.parquet`, and `ingest build` defaults to
`PROVISIONAL_VERSION`, which is still `0.9.0`. Everything downstream asks for
1.0.0 by name — the S matrices are on it, Stage 4 reads it,
`run_purity_conditioned` reads it. A 0.9.0 matrix does not fail to join. It
joins on a different gene set.

```
ls $BRP_DATA_DIR/processed/bulk/tcga_tpm_*.parquet
```

It runs fractions → the predictor check (enforced, not eyeballed) → the purity
producer → gate-and-arm on all four rungs, refuses before any compute if the
tree is dirty or an input is missing, and treats "the gate failed" as a result
rather than an error. Read `logs/brp_stage4.o<jobid>` when it lands. It does
**not** commit anything.

The step-by-step below is what that script does, for when you want to run a
piece of it by hand.

## Setup

```bash
ssh <you>@scc1.bu.edu
cd $BRP_PROJECT_ROOT          # if unset: cd /projectnb/rise-batteries/bode/guanylin
git fetch origin && git checkout submission/competitor-bench && git pull --ff-only
```

**Do not `git rebase` on this branch.** It carries stamped results; a rebase
orphans the commits their sidecars point at. That has already cost a recovery
once. `git pull --ff-only` or nothing.

```bash
module load miniconda
conda activate brp-w3
export BRP_DATA_DIR=/projectnb/rise-batteries/bode/guanylin/data
python -c "import src.bulk.deconvolution; print('ok')"
```

Confirm the bulk matrices exist before queueing anything:

```bash
ls -la $BRP_DATA_DIR/processed/bulk/tcga_{tpm,log2cpm}_1.0.0.parquet
```

If they are missing, build them first (this is the long step):

```bash
python -m src.bulk.ingest build
```

---

## Step 1 — fractions

```bash
python -m src.bulk.run_deconvolution \
    --rung all \
    --linearise-reference \
    --bulk $BRP_DATA_DIR/processed/bulk/tcga_tpm_1.0.0.parquet
```

Writes `results/<date>_<sha>/stage4_fractions.parquet` and
`stage4_predictor_checks.parquet`.

**Read the predictor checks before going on.** Every row must say `usable`. A
row saying `refused` means that method's fraction is a constant and nothing may
be regressed on it.

Exit codes: `0` fine · `2` scale mismatch refused (you forgot
`--linearise-reference`) · `3` no method produced a usable predictor — that is a
Stage 4 result and both tables are still written.

---

## Step 2 — gate, then arm

```bash
python -m src.bulk.run_stage4_variance --rung lineage
```

It resolves the newest committed fractions and purity tables on its own. To pin
them, pass `--fractions` and `--purity`.

Add covariates when you have them — without the flag you get the **marginal** R²
only, and the locked spec asks for the adjusted one too:

```bash
python -m src.bulk.run_stage4_variance --rung lineage \
    --covariates $BRP_DATA_DIR/processed/bulk/stage4_covariates.parquet
```

The frame needs `sample_id` plus `stage, age, sex, msi_status, purity, site`
(the `expression_models` context) and `plate`. Plate enters as fixed dummies:
the covariate lock's "too many for fixed effects" is an events-per-df argument
about Cox, and OLS on 624 samples has ~590 residual df left.

Exit codes: `0` a verdict was reached · `3` no usable predictor · `4` the
instrument gate failed, **which is itself the Stage 4 result** · `5` no verdict
could be formed.

Then run every rung. The spec's quotability rule applies: lineage and
crypt_position remain depth-confounded per PR #49 and must not be quoted. **If
no rung is both quotable and estimable, that is the Stage 4 result.**

```bash
for rung in epithelial lineage crypt_position best4; do
    python -m src.bulk.run_stage4_variance --rung $rung
done
```

---

## How to read what comes back

`stage4_variance_verdicts.parquet`, one row per (rung, method, R² kind, gene).

- `primary_verdict` — `confirmed` / `disconfirmed` / `indeterminate`, applied
  verbatim from the lock.
- `r_squared_kind` — **both `partial` and `marginal` are reported.** The lock
  asks for both and its arms say "R-squared" without saying which. That gap is
  real; picking one after seeing the numbers is what pre-specification prevents.
  If they disagree the driver warns and both stand.
- `negative_controls` — `breached` means every R² is an upper bound (ACTB/GAPDH)
  or the primary is indeterminate rather than confirmed (matched-null median).

**`indeterminate` is not a weak result.** A dead instrument produces exactly it,
so read `stage4_predictor_checks` and `stage4_instrument_gate` beside it. That
pairing is the difference between "the biology did not separate" and "the
measurement did not work".

And the lock's own words, which survive any result: bulk has no independent read
on cell-type identity (invariant 6). A confirmed prediction is **consistency**
with the compositional account, never mechanism.

---

## While you are on the cluster: the purity producer

Unrelated to Stage 4, same data, one command. It closes the repo's last
uncommitted-producer case and lets two dirty tables finally be deleted.

```bash
python -m src.bulk.run_purity_conditioned
```

Then `pytest tests/test_bulk_purity_conditioned.py` — the test that pins the
2026-08-18 tables as the only copies is written to **skip** once a clean twin
exists, and that skip is the signal it is safe to
`git rm results/2026-08-18_7c49e99/tcga_{premise_purity_conditioned,purity_expression_association}.*`.

## GSE39582

Affymetrix, not RNA-seq. The spec allows reporting that the S matrix cannot be
applied rather than forcing it. Do that — do not cross-platform-force a
signature built from 10x data onto an array cohort to fill a replication cell.

---

## What would actually fix the scale problem

A linearly-built reference from W1: the arithmetic mean of CP10K per cell type,
rather than `expm1` of the mean of `log1p`. It needs the cells, so it needs the
cluster and the `brp-w1` env. `src/reference/signature.py:build_signature`
selects markers on the log scale, which is right and should not change — what is
needed is a second emission of the *profile* on the linear scale, on the same
gene index and the same marker set.

Until that exists, every table from this chain carries
`reference_scale: linear_cp10k` and a `reference_note` saying it was derived by
`expm1`. Do not drop that note when writing anything up.
