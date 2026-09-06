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
first time, and inside the surviving mature cells of adenomas the
**terminal-differentiation tier is down while intestinal identity is retained**.
On the log fold-change scale the panel separates into **two blocks** —
{KRT8, ACTB, EPCAM, CDX2} and {MS4A12, GUCA2A} — with all eight cross-block
contrasts excluding zero and CDX2 indistinguishable from housekeeping. GUCA2A is
not separable within its block, so no gene-specific claim follows. **The first
version of that reading compared raw detection deltas across genes, which is not
a comparison the statistic supports** (§3, §6d); the conclusion survived the
correction, its stated evidence did not. **`docs/NEXT_AVENUES.md` reviews what is left**, and two items there
outrank a data hunt: the decomposition may be identifiable on adenoma, and the
instrument has never had a positive control though one is sitting in the atlas.

**The positive control RAN on 2026-09-06 and returned UNINTERPRETABLE — the
premise does not hold in the arm it is about (§6g).** ACTB moved 1.36× in the
methylated arm's tumour cells, and `premise_holds` was biased *toward* passing
by the very interval defect below and failed anyway. **The instrument's
sensitivity to a known silencing event remains untested, and a FIFTH
independent route says it cannot be tested on available data:** the premise has
held in exactly one cohort (adenoma, `Chen_2021_Cell`), that cohort carries no
MLH1 methylation annotation, and `Pelka_2021_Cell` is the only study in the
whole 49-study atlas that does. The cohort with the ground truth cannot support
the comparison; the cohort that supports it has no ground truth. **So every null
this project has produced stays uninformative rather than becoming evidence of
absence.** And the statistic every
reading ends in has been calibrated for the first time: **the percentile
bootstrap over patients is 0.82× the width it claims at n=10 and 0.53× at n=4,
by a closed form with no data in it** (§3a). Nothing is retracted; `best4`
contrasts are ~7% tests rather than 5% ones.

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

**A check that cannot fail reports success.** It has now been found **seventeen**
times, including six times inside guards written to prevent it, and twice
inside guards written *during* this work. Assume the next one exists.

**The seventeenth is not a check that could not fail — it is a quantity with no
check at all**, which is the same failure one step earlier: the MLH1 secondary
arm was defined twice, as 15 and as 19, and nothing compared them until a
cluster log printed one against a document holding the other. §6g.

**The fifteenth and sixteenth are not guards at all, and that is the news.**
They are an *interval* and a *power calculation* — the two places a project
states how sure it is. Same shape, same silence: an interval narrower than it
claims turns noise into a finding exactly as a guard that cannot fail turns
absence into a green light, and neither raises. See §3a; the measurement is
`results/2026-09-06_a0483ae/`.

The fourteenth is the most consequential yet, because it sat under the project's
**only positive result**. `specificity()` compared six genes' detection deltas to
each other; those genes span baseline detection 0.36 to 0.98, and the
sensitivity of a proportion depends on where it sits. Its guard could not have
caught it — the test fixture named `delta_detect` directly, so no baseline rate
existed anywhere in it and a gene at 0.98 and a gene at 0.44 were the same
object. The violating input is now committed: one uniform thinning of 0.75
across six baselines, no gene-specificity at all, and **detection reports 26 of
30 pairwise contrasts as real**. See §6d.

The five found on 2026-09-05, all in code written that day:

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
| **gene specificity** | one gene's detection delta vs another's | **a proportion's sensitivity scales with its baseline; the fixture had no baselines in it** |
| **the interval itself** | percentile bootstrap over patients | **0.82x the correct width at n=10, 0.53x at n=4 — a normal quantile where a t quantile is needed** |
| **the power calculation** | simulated power at n=10 | **only binomial noise in it, so it could not come out underpowered from patient heterogeneity** |
| **the secondary arm's size** | *nothing compared them* | **one arm, two definitions in two files — 15 and 19 — surfaced only when the cluster printed one against a document holding the other** |

The last four were found this week. The invariant-2 guard has now been fixed
**three times**; each fix covered the case just found, and the next was always
the input nobody had loaded yet.

**The rule the repo works to:** a guard needs a committed input that forces it to
fail. `tests/test_checks_can_fail.py` holds those — 38 of them, 40 with
parametrisation. If you add a
guard, add its failing input. If you cannot construct one, that is the finding.

---

## 3a. The interval is 0.82x the width it claims, and that is arithmetic

**Every reading in this project ends in a percentile bootstrap over patients** —
`premise_holds`, `summarise`, `specificity`, `control_log2_interval`, all at
`N_BOOTSTRAP = 10_000`. At n=44 that is fine. Below about n=20 it is not, and
**the mechanism is closed form with no data in it.**

The bootstrap distribution of a mean has standard deviation `s·sqrt((n−1)/n)` —
the plug-in, divide-by-n one. So the percentile interval is approximately

    mean ± z·s·sqrt((n−1)/n)/√n     against the correct     mean ± t(n−1)·s/√n

and the ratio is **`z·sqrt((n−1)/n)/t(n−1)`: a function of n alone.** Two errors
running the same way — a normal quantile where a t quantile is needed, and a
biased standard deviation where an unbiased one is needed.

| n | width vs correct | false positives, closed form | measured |
|---|---|---|---|
| 4 | 0.53× | 18.8% | 19.1% |
| 10 | 0.82× | 9.6% | 10.7% |
| 19 | 0.91× | 7.3% | 8.3% |
| 20 | 0.91× | 7.1% | 7.3% |
| 44 | 0.96× | 5.9% | 6.2% |

Nominal is 5%. **Simulation and arithmetic agree to within 1.1pp at every n, and
the excess is positive everywhere** — that gap is the skew a normal
approximation drops, so the closed form is a *floor* and the two routes describe
one phenomenon. It is not about MLH1 being rare, not about `cloglog`, and not
about this project: it is a property of the estimator at that n.

**Neither sophisticated repair works.** BCa is *worse* (17 of 20 cells
miscalibrated against the percentile's 14) because its bias correction and
jackknife acceleration are estimated from the same ten numbers. Raising
`N_BOOTSTRAP` does nothing — the error is not Monte Carlo error, so more
replicates estimate the wrong thing more precisely. **The plain Student-t
interval is calibrated at every n measured**: 0 of 20 cells, worst 5.9%.

`src/reference/interval_calibration.py` is the module;
`results/2026-09-06_a0483ae/` is the measurement.

### What it does to the committed tables — and what is NOT being restated

**Nothing is retracted and 138 tables are not being re-run.**

- **The adenoma `lineage` reading is n=44 → 5.9%.** Fine, and its margins are
  large: all eight cross-block contrasts exclude zero with room to spare.
- **The `best4` reading is n=20 → 7.1%.** A recorded caveat, not a retraction —
  and note it is the rung where the thinning null already fires
  (GRADIENT IS ABUNDANCE, 68.7% variance explained), so the block structure
  there was already resting on the load-bearing scale alone. **Quote `best4`
  contrasts as ~7% tests, not 5% ones.**
- **The 13-study meta-analysis is unaffected.** `src/harness/meta.py` is
  DerSimonian–Laird with a Higgins–Thompson prediction interval; it is not a
  percentile bootstrap over patients.
- **`premise_holds` is the one to think about, and it goes the wrong way.** It is
  an *equivalence* test — it asks whether the interval fits *inside* a tolerance
  — so a **narrower** interval fits more easily and the check is
  **anti-conservative** at small n. A HOLDS at n=10 is weaker than the same word
  at n=44. The premise results that matter (13 studies, k=11) are UNRESOLVED
  anyway, and a check biased toward HOLDS returning UNRESOLVED is a stronger
  negative, not a weaker one.

**The rule going forward:** any new reading below about n=20 reports the
Student-t interval and says so. Existing tables keep their numbers and gain this
paragraph.

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

**Read this first.** Both analysis paths ran to completion and terminated in
pre-committed negative results (§6a, §6b). Every table is committed and
clean-stamped; the repo has zero dirty tables and no uncommitted producers.

**One thing is now waiting rather than finished: the MLH1 positive control
(§6g).** Its code, pre-registration and calibration are committed and it needs a
single `qsub` against the cluster's atlas. It is the highest-leverage item here
because it is the first test of the instrument this project has ever had, and
both of its outcomes change how the write-up reads.

The other open items are a **write-up** (§6f) and a **data hunt** (§6e).

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

**Three firsts.** Per-patient deltas `results/2026-09-05_d869bdd/`
(`3a1af9f` holds the same deltas byte-for-byte but its *summary* is the
superseded pooled-rungs one reading `n_patients = 64`; do not cite it):

- **The premise HOLDS** — 44 patients at `lineage`, 20 at `best4`. It had never
  held anywhere in this project.
- **`best4` is estimable** — 32 of 44 patients cleared the mature-cell floor,
  median ~85 cells. Carcinoma had a median of 3, and the Stage 4 run returned an
  exactly-zero mature fraction on 95–99% of tumours there.
- **GUCA2A falls inside the surviving mature population**: −0.174
  [−0.243, −0.109] at `lineage`, over 44 patients, with the premise holding.
  That is a within-gene, between-arm comparison and the detection scale is
  fine for it. It is comparing that number to ANOTHER GENE's that needs the
  log fold-change scale — see below.

**What is down is a TIER, not the gene — and the first reading of that was
right for the wrong reason.** The gradient table this section used to carry
ranked the panel by raw detection delta. Detection deltas are **not comparable
between genes**: the sensitivity of a proportion depends on its baseline, and
this panel spans 0.36 to 0.98. The re-read is
`results/2026-09-05_9c43f4f/`, on the log fold-change scale
(`src/reference/detection_scale.py`), and it is a **two-block separation rather
than a gradient**:

| lineage, n=44, cloglog | KRT8 | ACTB | EPCAM | CDX2 |
|---|---|---|---|---|
| **MS4A12** | +0.761 \* | +0.581 \* | +0.688 \* | +0.659 \* |
| **GUCA2A** | +0.691 \* | +0.510 \* | +0.617 \* | +0.589 \* |

`*` excludes zero. **All eight cross-block contrasts exclude zero**;
`MS4A12 − GUCA2A` is −0.071 and contains zero. And the contrast the conclusion
actually turns on, which no committed row previously reported: **`CDX2 − KRT8`
is −0.102 and contains zero**, as are `CDX2 − ACTB` (+0.079) and `CDX2 − EPCAM`
(−0.028). CDX2 is not distinguishable from housekeeping.

So: **the cells retain intestinal identity and have lost terminal
differentiation output.** Same conclusion as before, now resting on a
comparison that supports it — and a sharper shape, because a *gradient* is
what a uniform thinning also produces and two blocks is not.

*Three things to carry, all of which change how this is quoted.*

**On the committed detection statistic the identity claim reverses.**
`CDX2 − KRT8` = −0.057 [−0.092, −0.021], excluding zero: read that way, CDX2 is
down relative to housekeeping and "identity retained" is false. Six contrasts
change verdict with the statistic and
`adenoma_specificity_disagreements.parquet` names them. `GUCA2A − MS4A12`
contains zero on **all three** statistics at **both** rungs, so the
"no gene-specific claim" conclusion does not depend on the choice.

**The `best4` numbers are ~7% tests, not 5% ones.** Everything in this section
rests on the percentile bootstrap over patients, and that interval is 0.91× the
width it claims at n=20 — a false-positive rate of 7.1%, by the closed form in
§3a. At `lineage` (n=44) it is 5.9% and the margins are large, so the two-block
reading is unaffected. At `best4` (n=20) it is a **recorded caveat, not a
retraction**: the committed numbers stand, and a contrast quoted from that rung
should be quoted as a 7% test. Note this compounds with the fact that `best4` is
also where the thinning null fires, below.

**The thinning null fires at one rung and not the other**, which is why it is in
the run rather than done once by hand. Fitting ONE common fold change across all
six genes — the most gene-unspecific model there is:

| rung | verdict | c | variance explained |
|---|---|---|---|
| lineage | STRUCTURE SURVIVES | 0.677 | 27.0% |
| **best4** | **GRADIENT IS ABUNDANCE** | 0.701 | **68.7%** |

At `best4` the raw detection ordering is **not readable as biology** and the
block structure there rests on the load-bearing scale alone.

**Neither block is perfectly flat, and the exceptions are informative.** At
`lineage` ACTB separates from KRT8 (−0.181) and from EPCAM (−0.107) — a control
moving relative to another control, inside the premise tolerance but not
nothing. At `best4` the four-gene block IS flat (0 of 12 contrasts exclude
zero), but CDX2 no longer separates from MS4A12 (+0.332) or GUCA2A (+0.253)
either: at n=20 it is unresolved between the blocks rather than assigned to one.

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

**The instrument has never had a positive control. One is now built,
pre-registered, and waiting on a cluster run.** The atlas carries
`MLH1_promoter_methylation_status` on 240,630 of Pelka's cells, patient-level,
from an assay rather than from expression. MLH1 silencing in methylated patients
is a known event. Asking whether the detection statistic can see it tests **the
instrument, not the biology** — and every null this project has produced
(UNRESOLVED ×3, UNRESOLVED at 13 studies, not-specific on adenoma) rests on an
instrument whose sensitivity to real silencing has never been shown. If it
cannot see MLH1, those nulls are uninformative rather than evidence.

**See §6g.** The design is fixed in `docs/prereg_g2_mlh1_within_stratum.md` and
the calibration it depends on is committed. What is left is one `qsub`.

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

## 6g. The MLH1 positive control — RAN 2026-09-06. UNINTERPRETABLE.

**Full account: `docs/prereg_g2_mlh1_within_stratum.md`, RESULT section.**
Tables `results/2026-09-06_4b1afca/` — **still on the cluster, not yet
committed.** 29 of 30 eligible patients scored; arms 10 methylated / 4
intact-MMRd / 15 unmethylated.

**The verdict is §5's first branch and it was pre-committed.**

    UNRESOLVED: control ACTB +0.443 [+0.116, +0.737] on log2 expression
    (tolerance 0.5)

ACTB is detected in 98.4%/99.3% of cells, so the premise assesses it on log2
expression rather than detection; it moved 1.36×. **The interval correction in
§3a does not rescue it** — the honest interval is *wider*, roughly
[+0.065, +0.821], straddling the tolerance by more. And the bias runs the useful
way: `premise_holds` is an equivalence test, so a too-narrow interval makes
HOLDS *easier*. The check was biased toward passing and failed anyway.

**MLH1 did fall, only in the methylated arm, and it does not count.**
−1.666 [−2.453, −0.880] methylated; +0.191 unmethylated; +0.517 intact-MMRd;
`patients_with_signal` 10/10 exactly as sized. That is the pre-registered
pattern. **The gate is ordered before the reading precisely so a result this
suggestive cannot be promoted by whoever wanted it**, and §5 binds. No claim
about MLH1 silencing follows.

### The result that outlives the verdict

**The positive control is not available on any currently available data.** The
loop is closed:

| | |
|---|---|
| the instrument can be validated only where | **the premise holds** |
| the premise has held in exactly one place | **adenoma, `Chen_2021_Cell`, n=44** |
| `Chen_2021_Cell` MLH1 annotation | **none — every cell null** |
| studies in the 49-study atlas carrying it | **`Pelka_2021_Cell` only** (22 meth / 40 no_meth) |
| and Pelka is | **where the premise just failed** |

**A fifth independent route to §2**, and the one that closes the
instrument-validation question rather than a biological one. Its consequence is
§5's third branch: an instrument whose sensitivity cannot be established cannot
be cited for what it fails to see, so **the project's nulls stay uninformative
rather than becoming evidence of absence.** That is weaker than a passing
positive control would have bought, and it is what the data supports.

### Two things to do

1. **Commit the cluster tables.** `results/2026-09-06_4b1afca/` is the only
   uncommitted result in the project. Until it lands, the RESULT section quotes
   numbers with no artifact behind them.
2. **Avenue 1a is next** — the decomposition at `best4` on `Chen_2021`, the
   project's original estimand, on the one cohort where the premise holds.
   Re-checked 2026-09-06 and it survives: the algebraic collapse does not fire
   on adenoma, and estimability at `best4` is 20/20 for GUCA2A and 19/20 for
   MS4A12. It needs the adenoma job to emit `frac_mature_*` (no committed table
   carries them) plus a cluster run. `docs/NEXT_AVENUES.md` §1a.

### A defect found by running it

**The secondary arm had two definitions, 15 and 19, and no check.** `arm_of()`
broke out `mlh1_intact_mmrd`, so the reported arm is 15; `interval_calibration`
sized the same arm as `mlh1_stratum != "mlh1_methylated"`, which is 19. The
prereg's §3 inherited the 19. Both were right about different questions given
one name in two files. **Not a check that could not fail — a quantity with two
definitions and no check at all**, which is the same failure one step earlier.
`src/reference/mlh1_arms.py` is now the single definition, with three failing
inputs. Re-sized at n=15 (`results/2026-09-06_3bd168f/`): immaterial — 96.3% vs
96.9% power at 75% silencing — and it does not touch the verdict, which rests on
the n=10 arm.

---

## 7. Running things

    pytest -q                      # expect 1485 passed, 22 env-only failures
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

    qsub src/reference/jobs/mlh1_positive_control.sh      # §6g -- RAN 2026-09-06

**All of these have already been run.** The MLH1 run's tables are the one
uncommitted result in the project: `results/2026-09-06_4b1afca/`, on the
cluster. They are here for re-derivation, not
because anything is pending.

Local-only, no cluster needed — it reads committed tables:

    python -m src.reference.jobs.adenoma_specificity      # the corrected path C read
    python -m src.reference.jobs.interval_calibration     # §3a, ~6 min
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
