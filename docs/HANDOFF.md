# Where this is, for someone arriving cold

**Written 2026-09-05, branch `submission/competitor-bench`.** Read
[CLAUDE.md](../CLAUDE.md) first for the invariants; this says what state they
are in and what to do next.

---

## 1. The one-paragraph state

The project asks whether a differentiation marker's loss in colorectal cancer is
*compositional* (the cells left) or *intrinsic* (the cells stayed and went
quiet). **The answer is that this is not identifiable on any data currently
available**, established by FOUR independent routes — the decomposition's
algebraic collapse, the coexpression premise at 13 studies, the Stage 4
deconvolution gate, and the bulk CIMP screen. That is a real, pre-committed
result, not a stalled analysis.

**Three analysis paths ran to completion on 2026-09-05** — Stage 4 (§6a), the
13-study coexpression meta (§6b), and the adenoma reading (§6d). Nothing is
half-finished.

Path C is the one positive: the premise HOLDS, `best4` is estimable for the
first time, and GUCA2A falls inside the surviving mature population — but
indistinguishably from a generic maturity marker, so no gene-specific claim
follows. **`docs/NEXT_AVENUES.md` reviews what is left**, and two items there
outrank a data hunt: the decomposition may be identifiable on adenoma, and the
instrument has never had a positive control though one is sitting in the atlas.

Two papers exist; the WMHS one is a methods paper about validation statistics
that cannot fail, and it is the near deliverable — **deadline 15 September
2026**.

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

**The coexpression reading cannot rescue it, and scale does not fix that.**
Built to sidestep the algebra by measuring per-cell detection inside a fixed
label rather than a variance split. Its premise — that the diseased cells are
still the same kind of cell — was UNRESOLVED on all three original cohorts, so
it was taken to the ICBI atlas: **13 studies, 122 patients**
(`results/2026-09-05_3380d15/`, meta at `results/2026-09-05_61ba221/`).

| control | pooled | I² | verdict |
|---|---|---|---|
| ACTB | +0.152 [−0.013, +0.317] | 62.8% | HOLDS |
| KRT8 | −0.453 | **87.6%** | UNRESOLVED |

Every control must hold, so the premise does not. **KRT8's per-study estimates
run from −1.177 to +0.088 — the studies disagree about whether their own two
arms are comparable.** That is not a precision problem that more patients fix.
Per-study verdicts agree: 3 HOLDS, 1 REFUSED, 7 UNRESOLVED, 2 UNDEFINED.

**This closes the "more data" question.** It was the one blocker that looked
like a power problem, and at four times the studies and four times the patients
it is a disagreement problem instead.

*A correction worth knowing, because it is the kind of inference to avoid.* On
the three original cohorts ACTB pooled to +0.487 with I² = 0.0%, and that was
read as "the cohorts agree, so the failure is precision" with a prediction that
ACTB could only resolve "by a hair". At k=11 ACTB pools to **+0.152** and holds
comfortably — the point estimate moved 0.335, because three similar cohorts
(range +0.431 to +0.586) were not a sample of the fourteen (range −0.589 to
+0.532). **I² near zero on a small, similar set is not evidence of homogeneity
in the population.**

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

**A check that cannot fail reports success.** It has now been found **thirteen**
times, including five times inside guards written to prevent it, and twice
inside guards written *during* this work. Assume the next one exists.

The four found on 2026-09-05, all in code written that day:

- The Stage 4 predictor check read `sd` as a proxy for "does this column carry
  information". On real data it **ranked backwards**: `best4`/nnls had 4
  informative samples of 675 and a *higher* sd than a column with 32, because
  its survivors were extreme. Both passed. Fixed by reading the non-zero share.
- The A1 verdict text still said a gap justified a rebuild after the premise
  behind that had been retracted — a claim that could no longer fire, left able
  to fire.
- The ICBI validation bar **documented** an ACTB log2 check and **implemented**
  a detection check. For a saturated control those are not near-substitutes:
  the premise can flip UNRESOLVED→REFUSED while the detection intervals still
  overlap and the bar returns PASS. Demonstrated, not argued.
- A test fixture seeded from `hash(gene)` — Python randomises string hashing per
  process, so the test exercised different data on every run.

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

**The ICBI atlas's `/X` is log1p-normalised. Raw counts are in
`layers["counts"]`.** `adata.X` is what any obvious code reaches for, and
detection at ≥1 UMI against log values is wrong while nothing raises. Measured:
`/X` runs 0.2795–5.404 over 3525 distinct values; `layers/counts` runs 1–289
over 155, starting `[1,1,1,1,2,1,2,5]`.

**Identifier spaces have bitten this repo four times.** S matrices are Ensembl,
Lee's GEO matrices are symbols, and the ICBI atlas's `/var/_index` is Ensembl
while its symbols live in a separate `GeneSymbol` column. Each time the symptom
was an *empty intersection reported as a finding*, never an error. When a lookup
returns nothing, suspect the identifier space before the data.

**A conda env file is not a conda env.** `env/w2_harness.yml` declares
scikit-learn; `brp-w2` has never been created on the cluster. A job hardcoding
it died. Wrappers now try several and report which one they got.

**Two filesystems, two quotas, and `/scratch` is node-local.**
`/projectnb/rise-batteries` and `/project/rise-batteries` are 50 GB each and
both are near full (§5). `/scratch` is `/dev/sda8` — per-node and purged, so a
file written by one job is invisible to the next. Same class as the login-node
`/tmp` trap.

---

## 5. Data: what is where

| | laptop | cluster | where |
|---|---|---|---|
| Lee GSE132465/GSE144735 | yes | yes | `$BRP_DATA_DIR/raw/lee/` |
| GSE178341 (371k cells) | **no** | yes | `$BRP_DATA_DIR/raw/GSE178341/` |
| TCGA bulk **1.0.0** | **no** | yes | `$BRP_DATA_DIR/processed/bulk/` |
| GSE39582 | **no** | yes | `$BRP_DATA_DIR/processed/bulk/` |
| **ICBI atlas, 30.44 GiB** | **no** | yes | **`/project/rise-batteries/bode/icbi/`** |
| ICBI obs cache | yes | yes | `$BRP_DATA_DIR/interim/icbi_obs.parquet` |
| `results/` (138 tables) | yes | yes | in git |

`data/manifest.csv` carries every file's url and sha256 and is the only record
that travels — verify downloads against it.

**The ICBI atlas is the one artifact NOT under `BRP_DATA_DIR`.**
`/projectnb` had 15 GB free against a 50 GB quota when it was fetched, so it
went to `/project` instead. Jobs reading it take `BRP_ICBI_DIR`:

    export BRP_ICBI_DIR=/project/rise-batteries/bode/icbi

Its manifest row records that. **Both quotas are now tight** — `/project` was at
45/50 GB after the fetch and `/projectnb` at 40/50 — so check `pquota` before
downloading anything else.

**TCGA is built at index 1.0.0.** `ingest build` defaults to
`PROVISIONAL_VERSION`, which is still `"0.9.0"`; a build without
`--version 1.0.0` silently produces matrices on a different gene set that do
not fail to join, they join wrongly. `src/bulk/ingest_cluster.sh` pins it.

---

## 6. What to do next

**Read this first: there is nothing half-finished to resume.** Both analysis
paths ran to completion and terminated in pre-committed negative results (§6a,
§6b). Every table is committed and clean-stamped; the repo has zero dirty
tables and no uncommitted producers. The two open items are a **write-up**
(§6e) and a **data hunt** (§6d).

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

## 6a. Stage 4 — RUN, and it terminated. Nothing to resume.

**Result: the instrument gate failed on every estimable rung. No R² is
reported.** That is the locked prespec's pre-committed consequence, taken.
Tables: `results/2026-09-05_d358109/`.

| rung | usable predictor? | gate (threshold 0.5) |
|---|---|---|
| epithelial | none — no maturity call | never reached |
| **best4** | **none — both methods degenerate** | never reached |
| lineage | nusvr only | **0.462** — fails |
| crypt_position | nusvr only | **0.479** — fails |

Only 2 of 8 (rung, method) pairs produced a usable predictor, both nu-SVR; NNLS
was degenerate at every rung. `best4` — the rung matching GUCA2A's biology —
returns an exactly-zero mature fraction on 95–99% of tumours and never reaches
the gate at all.

**The reference was correct when this ran, so the failure is not a scale
artifact.** This is worth reading carefully because an earlier version of this
file said the opposite. `run_full_reference` accumulates ONE pseudo-cell per
cell type carrying that type's *summed* counts (`run_full_reference.py:314`),
so the committed profile is `log1p(CP10K(summed))` — no within-type averaging,
no Jensen gap, and `expm1` inverts it **exactly** (4.7e-06, float32 noise;
pinned by `test_expm1_is_an_exact_inverse_for_the_committed_construction`). The
run passed `--linearise-reference`, which for these matrices *is* the exact
linear scale.

**Therefore a W1 linear rebuild is NOT justified** — it would emit `expm1` of
what already exists. `build_signature_sparse(profile_scale="linear")` exists and
is tested, and is correct for any future build that averages over *cells*; it is
redundant for this one.

`docs/STAGE4_RUNBOOK.md` has the commands if it ever needs re-running.

## 6b. Path B (ICBI) — RUN, and it terminated. Nothing to resume.

**Full detail: `docs/ICBI_RUNBOOK.md`.** The result is in §2. Three pieces were built and all three are committed:

| piece | module |
|---|---|
| extraction | `src/reference/icbi_slice.py` — CSR row reads from `layers/counts` |
| adaptation | `src/reference/jobs/icbi_coexpression.py` + `.sh` |
| meta-analysis | `src/harness/meta.py`, `src/reference/jobs/coexpression_meta.py` |

**`src/harness/meta.py` is new and general.** Invariant 4 has demanded
random-effects meta-analysis in three documents since week 0 and nothing
implemented it; this is that estimator (DerSimonian–Laird, heterogeneity gating
the verdict, Higgins–Thompson prediction interval). It is reusable by anything.

**The one thing that would silently ruin a re-run: `/X` is log1p-normalised.**
Raw counts are in `layers["counts"]`. `adata.X` is what any obvious code reaches
for, and detection at ≥1 UMI against log values is wrong while nothing raises.
`assert_raw_counts` checks the values, not the layer name.

The adaptation was validated before the other thirteen ran: `Pelka_2021_Cell`
IS GSE178341, and the ICBI path reproduced its committed result — same premise
verdict, overlapping ACTB interval, GUCA2A drift **0.0096** against a 0.15 bar.
That check runs in the job and exits 4 on failure.

### ICBI sizing, and two committed bugs found doing it

`paired_sample_summary` reported **0 paired patients in all 49 studies** and
`platform_summary` reported **0 plate-based cells**. Both were vocabulary
misses — the atlas writes "primary tumor"/"adjacent normal" and "Smart-seq2";
the code tested for "tumor"/"normal" and "smartseq2". Real numbers: **229 paired
patients across 24 studies**, **24,136 plate-based cells**. The tests could not
catch it because the fixture was written from the code's vocabulary rather than
the atlas's. Sizing tables: `results/2026-09-05_d241b35/`.

## 6c. Housekeeping — done

- Phase-5 plumbing check: **done**, runs end to end in
  `tests/test_bulk_deconvolution.py`.
- **Zero dirty tables across 138.** 16 of the 18 were verified frame-identical
  to their clean twins before deletion; the last two were the same two that had
  no producer, so `src/bulk/run_purity_conditioned.py` was written first and run
  on the cluster — it reproduced them (all 4 association rows within 0.0003 r²,
  all 32 conditioned rows at 100% verdict agreement) — and only then were the
  originals deleted. **There is no uncommitted-producer case left.**
- The **survival arm** was deliberately never built. The deconvolution leg
  closed, so there are no fractions to drive it, and the locked prespec
  excludes it by design (`not_prespecified`). It would need its own
  pre-specification.

## 6d. Path C — RUN, and it is the one positive result. See `docs/ICBI_RUNBOOK.md`.

`Chen_2021_Cell` **is** the Vanderbilt/HTAN polyp atlas (`dataset` reads
`VUMC_HTAN_*`, sample ids `HTA11_*`), and it was already inside the 30 GB object
fetched for path B. No data hunt was needed. 44 patients with a matched polyp
and their own normal.

**Three firsts, all on `results/2026-09-05_3a1af9f/`:**

- **The premise HOLDS** — 44 patients at `lineage`, 20 at `best4`. It had never
  held anywhere in this project.
- **`best4` is estimable** — 32 of 44 patients cleared the mature-cell floor,
  median ~85 cells. Carcinoma had a median of 3, and the Stage 4 run returned an
  exactly-zero mature fraction on 95–99% of tumours there.
- **GUCA2A falls inside the surviving mature population**: −0.174
  [−0.243, −0.109] at `lineage`, and the within-patient paired difference
  against ACTB is −0.139 [−0.197, −0.081], excluding zero.

**And then the discriminating test says it is not gene-specific.**

| contrast | lineage |
|---|---|
| GUCA2A − ACTB (control) | −0.139 [−0.197, −0.081] excludes zero |
| **GUCA2A − MS4A12 (identity)** | **−0.009 [−0.064, +0.045] contains zero** |

MS4A12 is a mature-colonocyte marker with no silencing story attached, and it
falls −0.165 against GUCA2A's −0.174. What is down is the **mature programme**,
not GUCA2A. That is the week-0 falsification rule firing a second time in a
second place: when the target and the identity markers return the same answer,
no gene-specific claim follows.

**The limit this exposes, which is the paper's thesis again.** The premise check
compares ACTB and KRT8 between arms — it asks whether these are still *cells*.
It cannot ask whether they are still equally *mature* cells, because housekeeping
does not vary with maturity, and a control that cannot move with the thing in
question cannot certify it. Two readings survive the premise holding —
coordinate down-regulation of the mature programme, or a mature label admitting
less-mature cells in the polyp arm — and this statistic cannot separate them.

## 6d-bis. What is left to try: `docs/NEXT_AVENUES.md`

Reviewed against the data, with the dead options marked. Two things from it
belong here because they change priorities:

**The decomposition may be identifiable on adenoma.** The collapse that killed
it in carcinoma (`i/c → −(f_N/Δf)`) does NOT fire here: m_T/m_N runs 0.37 to 0.95
across the panel and the bracket runs −0.05 to −0.63, nothing near −1. The
project's original method may work on this substrate. Same run, cells already
loaded.

**The instrument has never had a positive control, and one is available.** The
atlas carries `MLH1_promoter_methylation_status` on 240,630 of Pelka's cells,
patient-level, from an assay rather than from expression. MLH1 silencing in
methylated patients is a known event. Asking whether the detection statistic can
see it tests **the instrument, not the biology** — and every null this project
has produced (UNRESOLVED ×3, UNRESOLVED at 13 studies, not-specific on adenoma)
rests on an instrument whose sensitivity to real silencing has never been shown.
If it cannot see MLH1, those nulls are uninformative rather than evidence.

## 6e. Different data, not more of it — still true for survivorship

Every design here conditions cell identity on transcription, which is circular,
and §2 now shows four independent routes terminating. **More of the same data
has been tried and does not work** — that is what the 13-study result settles.

Ranked, all data-hunt-gated:

1. **Adenoma / early-lesion single-cell cohort** (HTAN, Vanderbilt). Highest
   leverage and needs no new method: it makes `best4` estimable — the
   resolution where the question is actually posed, and where the current
   median is 3 mature cells — and adds a normal→adenoma→carcinoma gradient
   constraining *when* loss happens. **Nobody has sized it yet.** The obvious
   first step is a read-only feasibility check in the shape of
   `src/reference/jobs/icbi_premise_feasibility.py`: does it carry both arms,
   the panel genes, and mature cells at `best4` resolution?
2. **Segmented spatial** (Xenium/CosMx, *not* Visium — a spot is still a
   mixture). Breaks the transcript circularity; survivorship still limits it.
3. **Cell-type-resolved methylation** (EPISCORE-style). The one assay that
   separates "surviving colonocytes are methylated" from "unmethylated
   survivors were death-selected". Weeks of new machinery. Note the CIMP screen
   weakened the methylation prior but did **not** close this — GUCA2A silencing
   need not be CIMP-tied.

**The caveat that survives every one of these.** Survivorship — GUCA2A-high
cells having been preferentially destroyed — is **not transcript-detectable**.
A resolved premise makes the detection reading interpretable; it never rules
that out. Longitudinal tracing would, and does not exist in public human data.

## 6f. The write-up

Two terminal results landed on 2026-09-05 and neither is in either paper yet:

- **Stage 4's gate failure on a verified-correct reference** (§6a). Note the
  arc: it was suspected of being a scale artifact, the suspicion was tested,
  and the reference turned out already correct — so the result strengthened.
- **The premise unresolvable at 13 studies because studies disagree** (§2),
  which is a materially stronger statement than "undecided at 3".

For the WMHS paper specifically, `tests/test_paper_numbers.py` ties its prose to
its tables: every figure in the third-guard paragraph is re-derived from
`results/*/coexpression_silencing*.parquet` and asserted against the literal
string in the `.tex`. All eleven assertions are mutation-tested. **If you re-run
that job, this test tells you exactly what to edit.**

---

## 7. Running things

    pytest -q                      # expect 1420 passed, 22 env-only failures
    ruff check src tests submission

The 22 are `anndata`, `diptest`, `lifelines` absent locally. That count is the
baseline; anything else is yours.

Cluster jobs are SGE — `qsub`, not `./`. Every wrapper refuses a dirty tree, an
unset `BRP_DATA_DIR`, and missing inputs *before* the compute rather than after:

    export BRP_PROJECT_ROOT=/projectnb/rise-batteries/bode/guanylin
    export BRP_DATA_DIR=$BRP_PROJECT_ROOT/data
    export BRP_ICBI_DIR=/project/rise-batteries/bode/icbi

    qsub src/bulk/ingest_cluster.sh                      # TCGA at index 1.0.0
    qsub src/bulk/stage4_cluster.sh                      # the whole Stage 4 chain
    qsub src/reference/jobs/fetch_icbi_atlas.sh          # 30.44 GiB, resumable
    qsub -v BRP_ICBI_STUDY=Pelka_2021_Cell \
         src/reference/jobs/icbi_coexpression.sh         # validate, then:
    qsub -v BRP_ICBI_STUDY=all src/reference/jobs/icbi_coexpression.sh

**All of these have already been run.** They are here for re-derivation, not
because anything is pending.

Local-only, no cluster needed — it reads committed tables:

    python -m src.reference.jobs.coexpression_meta        # newest ICBI run
    python -m src.reference.jobs.coexpression_meta \
        --deltas results/2026-09-04_975cf5c/coexpression_silencing.parquet

That second form is the three-cohort dry run, and it is worth understanding
before trusting the meta layer: it reproduces the answer this project already
knew (UNRESOLVED at k=3), which is how the combiner was checked.

Jobs do **not** commit their own tables. Read them, then
`git add results/<dir> && git commit && git push`.

**Merging, never rebasing.** `git pull --no-rebase` on this branch — the
sidecars point at specific shas and a rebase orphans them.
