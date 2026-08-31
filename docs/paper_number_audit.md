# The ICBINB draft: every number, checked

**W1 → whoever rewrites the paper · 2026-08-30**
Checked against `origin/main` at `66548a6`. Every figure below was recomputed from
`results/` or by executing the committed harness. Nothing was transcribed from a
document, including from this project's own memos.

Verdicts use the project's own three-state scale, because the same distinction
applies: **held** = reproduced; **unbacked** = no committed table carries the
columns needed to check it; **wrong** = contradicted by what I can compute. A claim
that could not be checked is not a claim that failed.

**Tally: 21 held · 11 unbacked · 3 wrong · 5 things the draft misses.**

---

## 1. Read this paragraph if you read nothing else

Three numbers are wrong and need changing: **1.37% → 1.3516%**, **4,400 → 2,800**,
and *"agree to three decimals"* → *"agree to within 3%, exactly as the algebra
predicts."* Eleven more have no committed table behind them, which collides with the
ethics section's promise that every number comes from one — so either run them under
invariant 10 or drop the promise. Two additions would materially strengthen the
paper: **Week 0's falsification rule**, which shows the §4 failure was predicted
rather than discovered, and **the competitor benchmark**, which is the only place the
paper's positive proposal is tested rather than asserted. And check **Figure 2's
y-axis**, because −1.11 is not a possible value of the quantity its caption names.

---

## 2. The three wrong numbers

### 2.1 · §3 — the crossover is 1.3516%, not 1.37%

Solving `sqrt(3p(1-p)) = 0.20` gives **p = 1.3516%**. At p = 1.37% the bound is
**0.2013**, which is *above* the tolerance — so at the stated figure the check
can still fire, and the sentence claims the opposite.

In a paper arguing that thresholds must be checked against their own reachability,
stating an unreachable threshold slightly wrong is the worst available error. It is
also in `docs/gate_memo_w2.md:1310`, so fix both.

The bound itself is correct and worth keeping — see §3.2 below.

### 2.2 · §5 — 2,800 runs, not 4,400

`src/harness/attenuation.py` commits to:

```
DEFAULT_MATURE_FRACTIONS = (0.40, 0.20, 0.10, 0.05, 0.02, 0.01, 0.0)   # 7
DEFAULT_SHIFTS           = (1.0, 0.8, 0.5, 0.25)                       # 4
n_replicates             = 50
ARMS                     = ("oracle", "bulk")                          # 2
```

7 × 4 × 50 × 2 = **2,800**. 4,400 implies 11 mature fractions, which is not the
committed grid. Either the run used a non-default grid — in which case say so and
commit it — or the count is wrong.

### 2.3 · §4 — "agree to three decimals" is 0 of 31 patients

Computed on `results/2026-08-28_6f81018/decomposition_summary.parquet`, tier A/B/D,
`lineage` · `stem_pole`:

| | |
|---|---|
| Patients where all 8 genes agree to 3 decimals | **0 of 31** |
| Patients where all 8 agree to 2 decimals | 1 of 31 |
| Median spread across genes | 2.69 — **93% of the value itself** |

Restricted to strongly depleted genes (`m_T/m_N < 0.05`), which is the fair reading:
0 of 29 to three decimals, 2 of 29 to two, **median spread 0.078 ≈ 2.7%**.

**The section's own algebra predicts exactly this.** With
`intr/comp = (f_N/Δf)(m_T/m_N − 1)` and `m_T/m_N` ranging over [0, 0.05], the bracket
ranges over [−1, −0.95] — a 5% spread, so ~2.7% observed is the prediction landing.
Rewriting the sentence as *"agree to within 3%, exactly the spread the algebra
predicts from the residual m_T/m_N"* makes the section **self-consistent** instead of
merely correct-in-spirit. As written, the paper claims a precision its own theory
forbids.

### 2.4 · Figure 2's y-axis — −1.11 is out of range

The caption says *"relative intrinsic change per gene."* That quantity is
`m_T/m_N − 1`, which is bounded below by **−1** whenever means are non-negative. The
stated tier-A minimum of **−1.11** cannot be a value of it.

Either the axis is a different quantity than the caption names, or the number is
wrong. Both are fixable; neither survives a careful reviewer.

For reference, what I compute on the committed tables (patient-median of
`m_T/m_N − 1`, 8 strata = 4 rungs × 2 axes, unmatched and matched):

| Tier | Range |
|---|---|
| A · compositional control | −0.952 … −1.000 |
| B · intrinsic control | −1.000 … **+0.146** |
| D · retained control | −0.951 … −1.000 |

The **qualitative claim is confirmed and stronger than the draft states** — tier D
sits inside tier A's range almost exactly, and tier B does change sign. Only the
numeric ranges are off.

---

## 3. The ledger

### 3.1 · §2 — Depth contaminates the label

| Claim | Verdict | What I get |
|---|---|---|
| Arms 2.36× → ~1.0× after matching | **held** | 2.36× → 1.02× |
| Matching costs three quarters of the cells | **held** | 1,616 of 6,372 kept — 75% discarded |
| Rows go 36/0/4 → 16/16/8 | **held** | exact, `lineage` · `stem_pole` · `normal` |
| Withdrawn: 46 pp, ρ = −0.92, arms 4.3× | **held** | all three documented as withdrawn |
| ρ = +0.33 against a ceiling of 0.52 | **unbacked** | committed value is **0.31**; 0.33 appears nowhere |

On the last row: `decomposition_lee_smc.meta.json` and `docs/open_decisions.md:3196`
both say `|rho| 0.31`. One of the two is stale — worth five minutes to find which.

The ceiling of 0.52 implies a prevalence of 10.02%, which is plausible and internally
consistent; 0.33/0.52 = 0.635, so the "63%" is right *given* 0.33.

### 3.2 · §3 — The confound check couldn't fire

| Claim | Verdict | What I get |
|---|---|---|
| Bound is `sqrt(3p(1−p))` | **held** | matches perfect-separation Spearman to 4 dp at every prevalence tested |
| 8 cells in 5,564 → bound 0.066 | **held** | 0.06563 |
| Reachability flips 64 of 64 → 16 | **held** | exact — 16 of 56 rows reachable, × 2 axes |
| The independent-maxima pairing bug | **held** | still live on `main`, see §4 |
| Reference arm at 0.17–0.18 | **held** | 0.188 / 0.174 / 0.171 reproduced exactly |
| Below p = 1.37% the tolerance is unreachable | **wrong** | see §2.1 |
| Depleted arm 0.6% prevalence, bound 0.137 | **unbacked** | 0.137 needs p = 0.63%. Tables give 0.42% (unfiltered) or **0.00%** (filtered) |
| Under tolerance on 32 of 48 | **unbacked** | no denominator of 48 exists; nearest are 40/56, 46/56, 86/112 |
| 0.185 vs 0.092 / 0.094 on the depleted arm | **unbacked** | per-arm ρ is not persisted in any table |
| 66 of 190 rows move, up to 3.9× | **unbacked** | needs per-arm ρ |
| 21 of 166 rows above their own bound | **unbacked** | needs per-arm ρ |

**The bound is the paper's best result and it is solid.** I verified it numerically
against perfect separation at p = 0.005, 0.0137, 0.05, 0.10, 0.30, 0.50 — agreement
to four decimals every time, including the ceiling of 0.866 at p = 0.5 that surprises
people. Keep this; it is the section's spine.

The unbacked rows are all the *same* gap: `depth_confound_reference.parquet` carries
only `worst_rho`, the max over arms. Everything requiring per-arm ρ is uncheckable
until the corrected diagnostic is run and committed.

### 3.3 · §4 — The control panel can't separate

| Claim | Verdict | What I get |
|---|---|---|
| Tier D lands inside tier A's range | **held** | D −0.951…−1.000 · A −0.952…−1.000 |
| Tier B changes sign between strata | **held** | −1.000 … +0.146 |
| Ratio = `(f_N/Δf)(m_T/m_N − 1)` | **held** | algebra confirmed |
| Replicates on Lee to two decimals | **held** | MS4A12 = −1.0000 across all eight Lee strata |
| Tier D −0.57…−0.82, tier A −0.62…−1.11 | **wrong** | see §2.4 |
| Ratios agree to three decimals per patient | **wrong** | see §2.3 |
| Eight strata cross two QC floors | **unbacked** | no QC-floor column in any committed table |

On the last row: the draft's strata are *QC floor × population definition × matched* =
8. The committed tables stratify by *rung × axis* = 8. Same count, different axes. If
the QC-floor runs exist they are not in `results/`.

**One structural note for the rewrite.** *"This estimator was never going to separate
the tiers, at any sample size"* holds **for genes with `m_T/m_N ≈ 0`**. Tier D was
chosen precisely because it was expected *not* to have that profile — so the primary
finding is empirical (tier D turned out depleted too), and the algebra explains why
that makes the ratio uninformative *afterwards*. State it in that order and the seam
closes. As written, a reviewer will find it.

### 3.4 · §5 — The validation curve measures the simulator

| Claim | Verdict | What I get |
|---|---|---|
| Oracle estimate reproduces realised truth exactly | **held** | exactly 0.0 across 1,400 oracle runs I executed |
| Zero variance across replicates | **held** | std = 1e-16 at every shift |
| Values are `1/(1−s)` times a constant | **held** | the draft's three numbers are exactly **481/480 × 1/(1−s)** |
| Shifts 0.8 / 0.5 / 0.25 | **held** | `DEFAULT_SHIFTS = (1.0, 0.8, 0.5, 0.25)` |
| 4,400 runs | **wrong** | see §2.2 |

The closed-form claim is exact to six decimals:

```
5.010417 = (481/480) × 1/(1−0.80)
2.004167 = (481/480) × 1/(1−0.50)
1.336111 = (481/480) × 1/(1−0.25)
```

**This is the strongest verified claim in the paper.** I ran the committed sweep and
the oracle arm reproduces the realised truth to *exactly* zero across every run. The
circularity is real, demonstrable in one line, and the section is correct as written
apart from the run count.

### 3.5 · §6 — Abstention, and its own calibration

| Claim | Verdict | What I get |
|---|---|---|
| 28 of 28 abstain under one definition, 24 of 28 under the other | **held** | exact, both axes, filtered vs unfiltered |
| Criteria: 0.5 shift, 90% coverage, 80% discrimination | **held** | `PREREGISTERED`, verbatim |
| NaN comparisons score abstentions as failures | **held** | reproduced — 0.470/0.500 shipped vs 0.940/1.000 defined |
| "Exactly half" | **held** | ratio exactly 2.000 on both rates |
| The committed cutpoint of 50 is not recovered | **held** | my run returns `ok=40, wide=40` |
| Returns 120 under one draw pool, nothing under the other | **unbacked** | needs the real cohort; no committed table |
| Published row (0.36, 0.31) vs (0.72, 0.62) | **unbacked** | mechanism verified; exact values cohort-specific |

The 28/28 and 24/28 come straight out of
`results/2026-08-28_0be41e6/g4_verdict_gse178341_by_arm.parquet`: `best4` filtered
tumour arm, `n_below_threshold = 28` of 28 on both axes; unfiltered, 24 of 28. Clean.

---

## 4. Two bugs the paper describes are still live on `main`

The draft writes about both in the past tense — *"we fixed the diagnostic"*, *"the
diagnostic now persists both arms' correlation"*. **Neither fix is in the repository.**

### 4.1 · The independent-maxima pairing — `src/harness/depth_confound.py:189, 199`

```python
worst_rho = max(rhos)          # max over arms
ceiling   = max(ceilings)      # max over arms — INDEPENDENTLY
rho_vs_ceiling = worst_rho / ceiling
reachable = ceiling >= rho_tolerance
```

Exactly the defect §3 describes. `tolerance_is_reachable` takes the *largest* ceiling
across arms, so a normal arm at 5% prevalence (ceiling 0.381) certifies the test as
capable of firing on a tumour arm at 0.4% (ceiling 0.112). That is the 64/64 → 16
flip, and it still ships.

### 4.2 · The NaN coercion — `src/harness/calibration.py:123–124`

```python
covered       = (rows["ci_low"] <= truth) & (truth <= rows["ci_high"])
excludes_zero = (rows["ci_low"] > 0) | (rows["ci_high"] < 0)
```

Both return `False` on a null interval, and `.mean()` then counts abstentions as
failures. The guard on line 102 only raises when *every* CI is missing, so a
partially-abstaining bin passes straight through.

I reproduced it: 10 valid replicates at coverage 1.00 plus 10 abstentions returns
**0.50**, and `n_replicates` still reads **20**, so the row looks fully populated.
That is the paper's own point about a coercion arriving as a boolean in a rate rather
than a number in a column — and it is still true of the shipped code.

### 4.3 · Why this matters for submission, not just hygiene

The ethics section claims:

> The figure code reads every number in this paper from a versioned result table
> rather than a transcription, each table carries a commit hash and a fixed seed.

For most of §3, all of §5 and most of §6, **no such table exists in `results/`**.
That sentence is the one a reproducibility-minded reviewer will test first. Either
run those analyses under invariant 10 and commit the parquet, or soften the sentence.
Submitting with both is the one option that doesn't work.

---

## 5. Five things the paper has and doesn't know it has

### 5.1 · The bound collapses at *both* ends, and §3 tells only the rare half

`sqrt(3p(1−p)) → 0` as p → 1 just as it does as p → 0. On the `epithelial` rung,
prevalence is exactly **1.000 on 222 of 222 rows** — every cell is called mature — so
the ceiling is exactly zero and the check is structurally incapable of firing there
too. `depth_confound_reference.parquet` confirms it: `worst_rho` is **NaN on all 64
epithelial rows**.

The paper's claim is currently *"rare labels are invisible to this check."* The true
claim is *"**degenerate** labels are invisible, at either extreme"* — which is
stronger, covers the coarsest rung as well as the finest, and explains a second
published result the team already withdrew. This is the single biggest free upgrade
available to the paper.

### 5.2 · The depleted arm is worse than §3 says

On the **filtered** population definition, `best4` tumour-arm prevalence is a median
of **exactly 0.000** — 25 of 28 patients have no mature tumour cells at all — so the
ceiling is exactly 0 and 46 of 56 rows cannot fire. The draft's "runs 0.6%"
understates its own finding.

### 5.3 · §5's statistic is undefined on the null arm

At `shift = 1.0` the parametric intrinsic term is exactly 0, so the recovery ratio
divides by zero and returns **−inf**. The design's own null control is the one place
the validation statistic cannot be evaluated. Same shape as everything else in the
paper, and free to add.

### 5.4 · The calibration collapses the three-state rule to two

`calibrate_cutpoints` returns **`ok = 40, wide = 40`** on my synthetic run — the two
cutpoints land on the same value, so the middle state disappears. If that reproduces
on the real cohort it is a sharper version of §6's point than *"the number came back
120 instead of 50"*, because it means the rule stops being three-state at all.

### 5.5 · §4's algebra predicts the precision the draft claims to observe

Covered in §2.3. Worth restating because it changes the section's status: the fix
doesn't weaken §4, it makes the theory and the measurement agree.

---

## 6. What to add — Week 0

There is no "Phase 0" in the repo. The concept is **Week 0 — shared setup**,
`execution_plan.md §3`: half a day, everyone present, three things frozen before any
analysis code existed — the panel, the labelling axes, the output schema. And this,
written into `config/panel.yaml` itself:

> **FALSIFICATION RULE:** if tiers A, B and D all return the same answer, the
> estimator is broken and no biological claim may be made.

**That is the rule §4 reports firing**, and the ethics section's "we make no
biological claim" is the consequence it pre-committed. The draft never mentions that
any of it was written down in advance.

It is enforced, not merely documented:

- Frozen in `config/panel.yaml` with a `frozen_on: week-0` stamp
- Mirrored by an exact-value assertion in `tests/test_freeze.py` — CI goes red if
  anyone edits the panel without editing the test
- Named in the logic of three separate modules (`src/harness/g1_amendment.py`,
  `src/bulk/gene_index.py`, `tests/test_reference_gse178341.py`)

And `README.md`'s "Every branch ends in a result" table pre-committed the outcome of
**both** branches that actually fired:

| Pre-registered branch | Pre-committed outcome | Fired |
|---|---|---|
| Tiers A/B/D don't separate | A methods and validation paper; no biological claim | yes — §4 |
| Most patients fall below the positivity threshold | A non-identifiability finding with diagnostics — the headline, not a caveat | yes — §6 |

**Why this changes how the paper reads.** Without Week 0, the draft is a team
reporting that its controls failed — which a reviewer may read as post-hoc
rationalisation of a null result. With Week 0, it is a team that wrote down in
advance what would make them abandon the biological claim, watched it happen, and
abandoned the claim. That is the difference between a failure report and a
demonstration that pre-registration works in this setting.

It also earns §7 a sixth takeaway that no other section does: **write the
falsification rule into the config file, not the protocol document, and assert it in
CI.**

**Suggested placement:** a short §1.5, or a paragraph at the head of §4, quoting the
rule verbatim with its `frozen_on: week-0` stamp. ~120 words, and the best
credibility-per-word in the document.

---

## 7. What to add — the competitor benchmark

Branch `submission/competitor-bench`, 21 tests, seed 20260829, six worlds × 200
replicates. I re-ran it: **6,000 rows, bit-identical** to the committed tables.

It measures the project's actual novelty claim — `README.md:83-86`, *"every existing
method returns a number; none flags that the intrinsic estimate is meaningless in a
tumour with no mature cells left"* — which had never been tested.

**Where the estimand does not exist** (`annihilated` world, no mature tumour cell
survives):

| Method | Returned a number | False confidence | Median \|intrinsic\| invented |
|---|---|---|---|
| `kitagawa+positivity` (ours) | 0 / 200 | **0.00** | — (refused all 200) |
| `composition-only` | 0 / 200 | 0.00 | — (*no intrinsic arm at all*) |
| `kitagawa-no-gate` | 200 / 200 | **1.00** | 7.99 |
| `pseudobulk-de` | 200 / 200 | **1.00** | 7.99 |
| `naive-delta-mean` | 200 / 200 | **1.00** | 20.00 |

On a scale where the normal-tissue per-cell mean is 20.0, the naive method reports
the *entire* 20-point drop as silencing in cells that do not exist.

**The counterweight** — refusing is trivially achievable by refusing always, so the
same methods are scored where the estimand *does* exist:

| Method | Detection rate | Median abs error |
|---|---|---|
| `kitagawa-no-gate` | 1.000 | 0.122 |
| `kitagawa+positivity` (ours) | **0.853** | **0.111** |
| `pseudobulk-de` | 1.000 | 3.435 |
| `naive-delta-mean` | 1.000 | 6.329 |

**The trade is 14.7% of detections for 100% of the false confidence** — and the gate
*improves* accuracy where it does answer, because the cases it drops are the noisy
low-n ones.

The load-bearing comparison is the **ablation**: `kitagawa-no-gate` is identical
arithmetic with the gate removed, differing in exactly one way. Every other method
differs in several at once.

**Why the paper needs it.** §6 currently *asserts* that abstention is the right
response to the four failures. The benchmark *demonstrates* it, against named
alternatives, with a stated price. It converts §6 from a design preference into a
measured trade. It is synthetic and the generative model is the one Kitagawa assumes
— `submission/FINDINGS.md` says so in its own limitations — but that limits what it
proves, not whether it belongs.

**Before quoting it, fix two things** (details in §9): the provenance record points
at a dirty tree, and `FINDINGS.md`'s signed errors disagree with its own parquet.

---

## 8. Results outside the failure catalogue

The paper is **not** just a failure catalogue, and framing it that way undersells it
twice. It is the pre-registered outcome of two named branches, with an abstention
mechanism no competitor has. Two more results exist that the draft excludes:

**The GUCA2A decomposition.** Loss is predominantly intrinsic on both cohorts, and
the direction survives depth matching — intrinsic −26.32, band [−79.7, −16.4],
essentially unmoved from −25.68 unmatched. The ethics section affirmatively disclaims
it. Given §4 shows the control panel cannot validate the decomposition, that is a
defensible call — but it should be a stated decision, not a silence.

**The bulk arm.** GUCA2A at 1.40% / 1.72% of normal against CDX2 at 94.7% / 84.8% —
a ~55-fold divergence between two markers of the same population, replicated across
TCGA and GSE39582. Fractions only, per invariant 6. Not mentioned in the draft.

---

## 9. Code state

Two branches, **neither pushed**. They are independent — the benchmark imports only
`src.estimator.kitagawa`, `src.harness.truth`, `src.harness.positivity` and
`src.common.provenance`, and touches nothing the other branch changes — so either can
merge alone.

| Branch | Contents | Tests | `ruff check` |
|---|---|---|---|
| `w4/lee-depth-quantile` | Depth target was 0.10 on Lee and 0.25 elsewhere — fix plus `n_counts`; a committed driver for the Lee decomposition | 1127 pass | clean |
| `submission/competitor-bench` | The benchmark, findings, three result tables | 1131 pass | clean |

Both carry the same 22 pre-existing failures — `diptest`, `lifelines` and `anndata`
absent from the local env, unrelated to either change and present beforehand. `main`
is byte-identical to `origin/main` at `66548a6` and the tree is clean.

### Three things to settle before either PR opens

**The ownership boundary.** `CONTRIBUTING.md` §2 gives W1 `src/reference/` and lists
`src/estimator/` under "do not edit". The `w4/` branch edits `src/estimator/lee_io.py`
and adds `run_lee_decomposition.py`. Naming a branch `w4/` does not make it W4's
commit — this needs W4 to own the change or say in writing that they accept it.

**The benchmark's provenance record points at a dirty tree.**
`submission/results/bench.meta.json` carries `git_dirty: true`,
`git_branch: w1/threshold-sweep`, `git_sha: 730bec0` — none of which is where this
code now lives, so the recorded sha does not identify the code that ran. Invariant 10
asks for more. The fix is one re-run on the branch; I verified the tables come back
bit-identical, so nothing about the result changes.

**`FINDINGS.md`'s signed errors do not match its own parquet.** The document reports
`+0.017` and `+0.007`; `sensitivity_where_estimable.parquet` says `−0.005054` and
`−0.002227` — **opposite signs**. The absolute errors match to three decimals, so it
is the signed column alone, probably carried over from an earlier run. Fix before
anyone quotes it in the paper.

Two smaller notes: CI runs `ruff check src tests`, so `submission/` is never linted by
CI, though `pytest -q` does collect its tests. And `submission/` is a new top-level
directory outside the `data/ src/ results/ env/` repo contract — nobody owns it, so
there is no violation, but it deserves a sentence in the PR body rather than leaving
reviewers to notice.

---

## 10. How this was checked

- **Result tables.** Recomputed from parquet under `results/`, joined where the
  paper's quantity spans two tables (e.g. prevalence from `mature_cell_counts_full`
  against ρ from `depth_confound_reference`).
- **The bound.** Verified numerically against perfect separation at six prevalences,
  200,000 cells each — agreement to four decimals throughout.
- **§5's oracle identity and §6's NaN coercion.** Verified by *executing* the
  committed harness on a synthetic cohort built from the test fixtures, not by
  reading the code. `data/` is empty locally, so anything requiring the real cohort
  is marked **unbacked** rather than guessed at.
- **The benchmark.** Re-run end to end and compared row-for-row against the committed
  parquet.

Where a number could not be checked, it is marked **unbacked** and the missing column
is named. None of the eleven is asserted to be wrong.
