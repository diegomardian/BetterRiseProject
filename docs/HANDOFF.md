# Where this is, for someone arriving cold

**Written 2026-09-05, branch `submission/competitor-bench`.** Read
[CLAUDE.md](../CLAUDE.md) first for the invariants; this says what state they
are in and what to do next.

---

## 1. The one-paragraph state

The project asks whether a differentiation marker's loss in colorectal cancer is
*compositional* (the cells left) or *intrinsic* (the cells stayed and went
quiet). **As of 2026-09-05 the answer is that this is not identifiable on this
data**, established twice by independent routes at two scales. That is a real,
pre-committed result, not a stalled analysis. Two papers exist; the WMHS one is
a methods paper about validation statistics that cannot fail, and it is the near
deliverable.

---

## 2. What is established, and how confident to be

**The decomposition cannot separate the mechanisms on this panel.** Algebraic,
not statistical:

    intrinsic / compositional = (f_N / Δf) x (m_T/m_N − 1)

As a gene's surviving per-cell mean → 0 the bracket → −1 and the ratio collapses
onto `−(f_N/Δf)`, a property of the *cell fractions*, identical for every gene on
the same labels. On GSE178341 that constant is **−5.85**, and GUCA2A (5.67),
GUCA2B (5.80), OTOP2 (5.85), CA7 (5.83) and MS4A12 (5.58) are not distinguishable
by it. Tier A was pre-registered as compositional and tier D as neither; both
return ~99% intrinsic. **The week-0 falsification rule fired and its
pre-committed consequence stands: no gene-specific mechanism claim from the
decomposition.**

**The coexpression reading cannot rescue it.** Built to sidestep that algebra by
measuring per-cell detection inside a fixed label rather than a variance split.
Its premise — that the diseased cells are still the same kind of cell — is
**UNRESOLVED on all three cohorts**, the control intervals straddling the
tolerance:

| cohort | n | control | interval | tolerance |
|---|---|---|---|---|
| GSE132465 | 7 | KRT8 | [−1.250, −0.289] | 0.5 |
| GSE144735 | 5 | ACTB | [+0.012, +1.112] | 0.5 |
| GSE178341 | 30 | ACTB | [+0.216, +0.640] | 0.5 |

Four times the patients halved the interval and it still straddles. Closing the
gap needs roughly **250–300 patients**. Table: `results/2026-09-04_975cf5c/`.

**The bulk arm agrees, by a different route.** Pre-registered CIMP screening
(`docs/prereg_cimp_specificity.md`, `results/2026-09-05_9203809/`):
**NOT SPECIFIC**, 0 of 2 references. GUCA2A falls *less* than CDX2
(+0.544 [+0.219, +0.878]) and no differently from MS4A12 (−0.147, contains
zero). So locus-specific promoter silencing is not the bulk-level story, and
**the 450k methylation leg has no hypothesis left to test.**

**The bulk arm reproduces across platforms.** 15 tables re-derived on Linux from
a freshly downloaded 3.1 GB cohort against the Windows originals: 11 bit-identical,
4 differing only at floating point (max *relative* 2.0e-16 to 9.3e-13) plus one
`int32→int64`. Zero gene-model drift.

---

## 3. The recurring defect, which is also the paper's thesis

**A check that cannot fail reports success.** It has now been found nine times,
including three times inside guards written to prevent it. Assume the next one
exists.

| where | the check | why it could not fail |
|---|---|---|
| recovery curve | estimator vs known truth | estimator cancels; curve measures the generator |
| calibration grid | where criteria first hold | grid cannot reach the crossing |
| depth confound | ρ vs tolerance 0.20 | √(3p(1−p)) ceiling below the tolerance |
| calibration rates | coverage/discrimination | `nan > x` is False, abstention scored as failure |
| width gate | interval vs threshold | same coercion |
| **premise control** | ACTB/KRT8 detection | **saturated at ~1.00, nowhere to fall** |
| **premise verdict** | point estimate vs tolerance | **no interval; flipped with the seed** |
| **invariant 1** | `None` is not `0.0` | validating writer used by 1 of 26 call sites |
| **invariant 2** | targets absent from signature | classifier with no reject option: a positional index read as "symbol" |

The last four were found this week. The invariant-2 guard has now been fixed
**three times**; each fix covered the case just found, and the next was always
the input nobody had loaded yet.

**The rule the repo works to:** a guard needs a committed input that forces it to
fail. `tests/test_checks_can_fail.py` holds those — 21 of them, 23 with
parametrisation. If you add a
guard, add its failing input. If you cannot construct one, that is the finding.

---

## 4. Traps that have each cost real time

**`BRP_DATA_DIR` unset does not fail — it reads the wrong disk.** `paths.py`
falls back to `REPO_ROOT/data`, which exists and is nearly empty. Two cluster
jobs died on this. The coexpression wrapper now refuses; other jobs do not.

**Export both variables in every cluster shell**, or `#$ -V` carries an
environment without them:

    export BRP_PROJECT_ROOT=/projectnb/rise-batteries/bode/guanylin
    export BRP_DATA_DIR=$BRP_PROJECT_ROOT/data

**Never pool cohorts.** Invariant 4, and it bites: a pooled coexpression mean
read −0.416 with an interval excluding zero, which neither cohort supported
alone.

**Do not rebase a pushed branch carrying stamped results.** The stamps point at
commits the rebase destroys. That is how `e5ebdc3` died (now closed).

**`newest()` in `paper/*/_tables.py` resolves by mtime**, so a stale table from
today beats a good one from today. Not fixed.

**Ownership** (CONTRIBUTING §2): `src/bulk/` is W3, `src/estimator/` is W4,
`src/harness/` is W2. Two files here were changed across that line on explicit
instruction, flagged in their commit messages.

---

## 5. Data: what is where

| | laptop | cluster |
|---|---|---|
| Lee GSE132465/GSE144735 | yes | yes (fetched 2026-09-04, sha-verified) |
| GSE178341 (371k cells) | **no** | yes |
| TCGA + GSE39582 bulk | **no** | yes |
| `results/` (~120 tables) | yes | yes |

`data/manifest.csv` carries every file's url and sha256 and is the only record
that travels — verify downloads against it. GSE178341 needs the cluster: 32–64 GB.

---

## 6. What to do next

**The WMHS paper is the outstanding deliverable. Deadline 15 September 2026,
AoE.** As of 2026-09-05 it carries the week's findings: `allow_dirty` and
`extra_meta` were already in Appendix A, and the saturated control plus the
three-cohort UNRESOLVED result went into §3 as the *third* withdrawn guard,
with §5 and the conclusion carrying the non-resolution (commit `e9f0ef9`).

Two things about the paper are still open, and one is blocking:

- **BLOCKING — the page limits are unverified.** `neurips_2026.sty` is
  deliberately not vendored, so `./build.sh` cannot run. Against a
  geometry-matched stub the full build grew one page and the short build's main
  text did not move; the stub runs ~1.4× long, so the estimate is 7 → 7.7 of 9
  real pages. **Download the official style, run `./build.sh`, then
  `./check_anonymity.sh`.** Both must pass before submission.
- The cross-platform reproduction (11/15 bulk tables bit-identical
  Windows→Linux, 4 differing at ≤9.3e-13 relative) is still unclaimed anywhere
  in the paper. It is a real result and it has no home yet.

`tests/test_paper_numbers.py` now ties the prose to its tables: every figure in
the new paragraph is re-derived from `results/*/coexpression_silencing*.parquet`
and asserted against the literal string in the `.tex`, so a re-run table and a
stale sentence can no longer diverge silently. All eleven assertions are
mutation-tested. **If you re-run that job, this test tells you what to edit.**

## 6a. Stage 4 — built 2026-09-05, ready to run

The variance arm is implemented and the whole chain runs end to end in the test
suite. See **`docs/STAGE4_RUNBOOK.md`** for the exact cluster commands.

- `config/stage4_prespecification.yaml` is **locked** (2026-09-05). Only the six
  lock fields changed; a test parses both sides of that commit and fails if the
  prediction moved with it.
- `src/bulk/deconvolution.py` → fractions, `instrument.py` → the pre-committed
  gate, `variance_arm.py` → the percentile-within-matched-null statistic,
  `run_deconvolution.py` and `run_stage4_variance.py` → the two drivers.

**Read this before running it.** The committed S matrices are mean-of-log1p and
TCGA TPM is a linear mixture. Deconvolve one against the other and NNLS returns
`mature_colonocyte_fraction` as **exactly 0.0 on every sample** — and the
pre-committed instrument gate **passes anyway** at r = 0.881, because it reads
the non-epithelial aggregate rather than the epithelial-internal split the
analysis uses. A constant predictor gives every gene R² = 0, which is what the
pre-registered arm expects for GUCA2A. `--linearise-reference` is required and
is an approximation; the real fix is a linearly-built reference from W1.

Blocked on the cluster only: the TCGA matrices are gitignored and live there.

## 6b. ICBI — sized 2026-09-05, and two committed bugs found doing it

`paired_sample_summary` reported **0 paired patients in all 49 studies** and
`platform_summary` reported **0 plate-based cells**. Both were vocabulary
misses — the atlas writes "primary tumor"/"adjacent normal" and "Smart-seq2";
the code tested for "tumor"/"normal" and "smartseq2". The real numbers are
**229 paired patients across 24 studies** and **24,136 plate-based cells**. The
tests could not catch it: the fixture was written from the code's vocabulary,
not the atlas's. Fixed, and a new guard caught the same defect a third time.

Sizing (`src/reference/jobs/icbi_premise_feasibility.py`): at ≥100 epithelial
cells per arm the atlas carries **14 studies / 136 paired patients, of which 12
studies and 85 patients are new**. 96.9% raw counts, so invariant 4 is
satisfiable; 15 reference genomes, which the shared index has to absorb.

**Its table is not stamped yet** — `results/` writes refuse a dirty tree, and
`paper/icbinb/` has uncommitted edits. Commit or stash those and re-run.

## 6c. Then, in cost order:

1. ~~Finish the Phase-5 plumbing check~~ — **done.** S matrix → synthetic bulk →
   NNLS + nu-SVR → stamped parquet runs end to end in
   `tests/test_bulk_deconvolution.py`.
2. **Housekeeping** — 18 dirty tables across 7 directories, every one now
   superseded by a clean twin; one `git rm`. And
   `tcga_premise_purity_conditioned` / `tcga_purity_expression_association` have
   **no committed producer** (`purity_conditioned_check()` exists and is tested,
   nothing calls it) — writing that driver is the last instance of Appendix A
   item 3.
3. **Stage-4 deconvolution → survival**, cluster, W3/W2 scope. Fractions only
   (invariant 6).

**To actually settle the mechanism question you need different data, not more of
it.** Every design here conditions cell identity on transcription, which is
circular. Ranked: an **adenoma/early-lesion single-cell cohort** (runs the
existing pipeline at `best4`, the resolution where the question is posed and
where the median is 3 mature cells — substrate fix, no new method); **segmented
spatial** (Xenium/CosMx, not Visium — a spot is still a mixture); then
cell-type-resolved methylation.

---

## 7. Running things

    pytest -q                      # expect 1239 passed, 22 env-only failures
    ruff check src tests submission

The 22 are `anndata`, `diptest`, `lifelines` absent locally. That count is the
baseline; anything else is yours.

Cluster jobs are SGE — `qsub`, not `./`:

    qsub src/reference/jobs/coexpression_silencing.sh

A job writing results needs a **clean tree**; `write_versioned_table` refuses
otherwise, and the wrapper refuses earlier so you find out in seconds rather
than after the compute.
