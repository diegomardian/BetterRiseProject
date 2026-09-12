# Decision record — submission scope, and the next project

**2026-09-11. This is a decision, not an analysis.** It records an external
audit of the checkout at `1b3e9b2`, the corrections that changed the decision,
and what is and is not funded for the next month. Where it conflicts with
`docs/HANDOFF.md` §6f or `docs/NEXT_AVENUES.md`, **this document is the later
decision.**

Verification status of every numerical claim is in §7. The audit's file
references were all checked and all exist; `1b3e9b2` is a real commit on this
branch's history.

---

## 1. The decision, in four lines

- **Yes** to a modest journal submission in one month, as a
  **methods-and-limitations paper targeting BMC Bioinformatics**.
- **No** to spending that month turning the adenoma result into a standalone
  mechanistic CRC paper. The observation deserves publication; the available
  evidence does not establish that surviving mature cells were silenced.
- **Next project: auditing published single-cell validation pipelines, with
  learned selection of informative stress tests.** Several months, not another
  four-week rescue.
- **Stop the active search for another CRC mechanism result.**

The stated rationale for the pivot: the comparative advantage here is the
experience of constructing counterexamples and correcting scientific claims,
not privileged CRC data.

## 2. Corrections that materially affect the decision

Each one weakens a claim this repository has been making. They are listed with
what the committed artifact actually says.

| repository claim | what the artifact says |
|---|---|
| "~110 versioned tables; ten invariants" | **283 Parquet tables at `1b3e9b2`** (287 at `ec93216`), **13.5 MB**, and **eleven** invariants. |
| "Carcinoma collapses onto **−5.85**" | The paper-number tests use depth-matched inputs: the cohort-summary limit is **+6.94**, with six near-floor genes at **6.19–6.63**, a **1.07-fold** spread. **The negative constant in HANDOFF is stale and sign-inconsistent with the tests.** |
| "Adenoma spreads 11.5-fold (0.33–2.03 / 0.27–3.13)" | Correct **for the formula applied to cohort-averaged inputs**, giving 0.272–3.126. **It is not the spread of the patient-level decomposition summaries.** Different aggregations; label them separately. |
| "`best4` has 0/20 usable" | 0/20 are `ok` under either calibrated candidate, but **6/20 retain an intrinsic estimate under SMC's 70/42 cutpoints**, against 0/20 under KUL3's 400–800/70. **"Not ok" and "no estimate" are different.** |
| "PFI HR 0.9697, n=624" | HR and interval reproduce; the fitted ABSOLUTE-adjusted PFI model is **n=489, 124 events**. **624 was the feasibility population**, not the model's. |
| "Neither paper contains the newer findings" (HANDOFF §6f) | **Stale.** The current WMHS source already carries adenoma non-collapse, the three-grid analysis, the interval material and newer simulator experiments. **They cannot be budgeted as untouched novelty.** |
| "GUCA2A is on no stock panel" | **Xenium Colon v1 carries GUCA2A and MS4A12.** The verified obstacle is that **none of the four audited panels supplies the required target/control combination** — not universal absence. (`HANDOFF` §1 states this correctly; shorter paraphrases do not.) |
| "Everything listed is on disk or manifest-retrievable" | **13 of 703 manifest entries are present locally**, 304.8 MB, all 13 matching recorded size and SHA-256. The atlas, Pelka matrix and TCGA inputs are absent locally. **Crowell files are not in this checkout's manifest** despite HANDOFF saying otherwise. |

Local disk is ~87 GiB free; **9.8 GB is retained as the planning ceiling**
because it describes the recorded cluster constraint, which the audit could not
re-verify.

`icbinb-bio-submission.pdf` was **not found** in the checkout, so the exact
submitted PDF and the submission status of the ICBI-NB paper are **unverified**.

## 3. The methods paper — scope, and what it may claim

**Central claim, as fixed:**

> In one colorectal single-cell analysis workflow, executable counterexamples
> expose validation checks that miss their stated failure conditions, while
> calibration and sensitivity analyses delimit which descriptive estimates
> remain reportable **without establishing cell-intrinsic silencing**.

The contribution is **the executable audit and its measured consequences**.
Neither "we found many bugs" nor "bootstrap intervals can under-cover" is
sufficient novelty alone.

**What the adenoma evidence supports — stronger than a discarded negative,
weaker than "IS identifiable":**

- Across 44 paired patients the ratios of cohort means reproduce: GUCA2A 0.374,
  MS4A12 0.383, CDX2 0.791, EPCAM 0.737, ACTB 0.834, KRT8 0.946. A substantial
  descriptive pattern **in the selected, depth-matched populations**.
- Mean absolute intrinsic shares are 0.715 (GUCA2A, n=43) and 0.709 (MS4A12,
  n=42) under symmetric weighting. These are summaries of `|i|/(|i|+|c|)` —
  **not** estimates that 71% of biological loss was caused by silencing. The
  statistic discards signs and excludes the interaction from its denominator.
- Six of eight cross-block contrasts survive four-statistic agreement at
  `lineage`; **one of eight** at `best4`. **Failure to distinguish the two
  target genes does not establish their equivalence.**
- Crowell reproduces (GUCA2A −1.307 [−1.872, −0.741], MS4A12 −1.241 [−1.874,
  −0.609], negative in 7/7). **Supporting evidence for a tissue-level pattern**,
  not independent replication of the decomposition or its mechanism.
- Passing the algebraic falsifier shows **this particular degeneracy is
  absent**. It does not establish causal identifiability: which cells enter a
  transcript-defined population, composition within the coarse label, capture
  and survival can all still move its mean.
- Guanylin loss in early lesions **is not itself new** — the published
  APC/guanylin work already provides biological evidence here. The new
  contribution must be the defensible decomposition or measurement framework,
  which is exactly where the remaining limitations bind.

**Is the single-cohort adenoma result publishable?** Yes, as a clearly
exploratory descriptive analysis — single-cohort status alone does not preclude
publication. **For this group, this claim and this month: only inside the
methods-and-limitations paper.** Not as a standalone biological submission to a
journal expecting a mechanistic advance.

### The three reviewer objections that should drive the work

| objection | closable in one month? |
|---|---|
| "You have shown a change in selected cells, not silencing in retained mature cells." | **Structurally unclosable** with available deposits. Crowell supplies no control; Becker fails the design; MLH1 does not validate the CRC instrument. **Remove the causal claim.** |
| "Your positive result depends on post-hoc statistics, labels and estimability rules." | **Partially.** Publish all fixed weightings, denominators, exclusions and statistic disagreements; recompute multiplicity-sensitive summaries. No retrospective specification or DIS/VAL split restores independent confirmation. |
| "This is familiar testing advice, an internal case study, and overlapping workshop material." | **Partially.** Needs a reusable tool, blinded challenges beyond the fixtures used to build it, and a precise workshop-overlap table. Field-wide prevalence needs the external audit. |

### Two conceptual corrections required before submission

1. **The bootstrap width expression is a normal-approximation comparison for the
   sample mean**, not an exact distribution-free formula for percentile-bootstrap
   coverage. The committed experiments do reproduce **14/20 miscalibrated
   percentile cells, 17/20 BCa, 0/20 Student-t**, with Student-t rejection
   3.00–5.93%. That supports a result **about the tested generators**. It does
   not establish that BCa is generally worse, that the approximation is a
   universal lower bound, or that Student-t is universally calibrated.
2. **The WMHS argument overstates what equality to "realised truth" proves.** A
   sample mean equals the realised sample mean exactly; simulation can still
   validly measure bias, variance and coverage against a *population* mean, and
   a deterministic function of the same sufficient statistics need not be the
   same function. The claim to audit is **which assumption or performance
   property this validation cannot challenge**, not "zero residual means invalid
   validation." Likewise **a ratio of median absolute residuals is not a variance
   decomposition** and does not support "88% of the curve is generator noise."

## 4. The four-week plan

Filenames are **proposed outputs, not results already obtained**. Written under
the normal date/SHA directory with fixed seeds.

| week | work | deliverable and purpose |
|---|---|---|
| **1** | Two people classify the 21 ledger entries into logical impossibility, low power, implementation error, provenance/reporting error. A third constructs previously unseen clean and violating cases; the fourth audits workshop overlap and the mathematical claims. Include a correctly calibrated sample-mean simulation as a clean control. | `claim_check_inventory.parquet`, `blinded_guard_challenge.parquet`, overlap matrix. Show the audit separates a missed scientific failure from a legitimate algebraic identity. **The 55 existing tests are development material, not held-out evaluation.** |
| **2** | Re-run interval and cutpoint experiments on the existing Lee files and committed patient summaries. **Keep SMC and KUL3 separate.** Add skewed/heavy-tailed and boundary regimes; compare original checks, basic input checks and the proposed audit. Quantify Monte Carlo uncertainty. Report calibration crossings as **brackets, not precise thresholds**. | `interval_stress_calibration.parquet`, `cutpoint_crossing_brackets.parquet`, `guard_challenge_results.parquet`. |
| **3** | Recompute Chen summaries from committed patient inputs: all three weightings, all required rungs, both denominators, all four statistics, patient counts and missingness. Reproduce DIS/VAL and Crowell **without opening new endpoints**. | `adenoma_claim_sensitivity.parquet`, `estimability_attrition.parquet`, `crowell_summary_reproduction.parquet`. **Label as reproduction from derived inputs, not raw-data replication.** |
| **4** | Package the audit for another researcher; a teammate reproduces it in a fresh environment. Write the paper with the adenoma example **subordinate** to the validation question. Generate every reported number from tables; disclose workshop overlap. **Resolve the recorded review debt for shared-code changes.** | Submission-ready manuscript, executable examples, reproduction instructions. Submit one paper. |

**Resources.** No GPU, no new large dataset. Lee matrices ~185.5 MB compressed;
results 13.5 MB. Budget **under 2 GB** additional working storage excluding a
fresh environment; process synthetic matrices sequentially. **Measure the first
calibration batch before scheduling the full run** — the expanded experiments
are untimed; the existing targeted tests took seconds.

**Week-one stop rule.** If the blinded challenge finds **zero additional
decision-relevant failures** beyond ordinary finite-value, identifier and
input-validation checks — once legitimate estimator identities and estimand
mismatches are excluded — **withdraw the methods submission**, finish the
workshop record and redirect. A longer taxonomy without demonstrated additional
utility is not enough.

## 5. The next project — four candidates, ranked

### 1 · Learn which stress tests expose misleading validation in published pipelines — **COMMIT TO THIS**

> Can a learned testing policy find reproducible, scientifically consequential
> validation failures in previously unseen pipelines more efficiently than fixed
> checklists and random testing?

It converts existing experience into a testable capability. **A manual audit
establishes individual failures more cheaply; ML becomes necessary only when
selecting among expensive, interacting stress conditions is a measurable search
problem.**

Learn a surrogate from pipeline characteristics and stress conditions to a
*verified failure* — the published validation still looks acceptable while a
separately specified target property fails. Conditions vary donor heterogeneity,
label contamination, dependence, rare-population prevalence, generator
misspecification. **The model selects the next expensive experiment; it does not
decide whether the biological claim is true.**

Audit subjects: Splatter, scDesign3, scMultiSim, each **against its own stated
capabilities** — a simulator is not wrong for omitting what it never promised.
Systema as a later, deliberately different transfer task. **These are candidate
subjects, not methods found defective.**

Data: the two checksum-verified Lee cohorts, the 13.5 MB result corpus,
generated experiments. A 2,000 × 1,000 pilot is 8 MB per dense float32 matrix;
retain seeds and summaries, not every matrix. **1–3 GB.** Add the 79.8 MB
Replogle pseudobulk only for questions its aggregation supports.

**The novelty bar is substantial** — simulator realism, single-cell false
discoveries, perturbation-metric confounding and adaptive falsification are all
already published. "We applied Bayesian optimization to find bad examples" does
not clear it. What would: an explicit **published claim → promised validation →
counterexample → changed conclusion** chain; published implementations and
reproducible repairs; **held-out pipeline and failure families**, not new seeds;
equal compute against manual checklists, random search and ordinary Bayesian
optimization; **independent confirmation after adaptive search**, so searching
until a failure appears does not manufacture significance.

**The 55 tests are seed examples, not 55 independent training observations. The
21 ledger entries are a starting taxonomy, not a measured field-wide rate.**

- Target: **Genome Biology**. Nature Methods only with broader impact and
  particularly strong external examples. **Neither bar is cleared today.**
- Killer risk: the failures are trivial input errors or known simulator
  limitations and learned search adds nothing.
- Week-one test: pin two external packages; build a small prospectively
  specified challenge space; compare surrogate-guided search against random
  search and a fixed checklist at equal budgets; repeat counterexamples on
  independent seeds. **Proposed** investment threshold: ≥2× independently
  confirmed failure regimes at equal cost without more false accusations — a
  business decision, not a result.
- Timeline: **12–16 weeks minimum; 4–6 calendar months** alongside other work.
  If ML adds nothing but the audit changes important conclusions, **publish the
  audit without an ML novelty claim.**

### 2 · Perturbation prediction on Replogle — real ML, poor current bet

The worthwhile question is generalization to **unseen perturbations or
contexts**; randomly splitting cells from the same perturbations does not
establish it. Must beat no-change, mean-response and linear baselines on
perturbation-specific endpoints. The generic foundation-model comparison is
crowded — strong simple-baseline results and a 27-method / 29-dataset benchmark
already exist.

Disk, from the repo's own gate: K562 essential pseudobulk 79.8 MB (present);
K562 essential single-cell **10.66 GB** (already over the ceiling); RPE1
single-cell 8.70 GB (leaves ~1.1 GB); K562 genome-wide **65.83 GB**
(incompatible). **These are files within one resource, not interchangeable
datasets.**

Week-one test: establish the strongest simple baselines on a **locked
perturbation-held-out** task before model development; if performance vanishes
when shared response and target leakage are controlled, **stop**. Timeline 4–6
months minimum. **Rule out as the main pivot.**

The existing Replogle sensitivity protocol is a **different project and stays
distinct**: it can test a narrow measurement limitation; it cannot validate CRC
silencing, remove survivorship, or become a prediction paper by adding a neural
network.

### 3 · Clinical selective prediction — feasible data, insufficient differentiation

Worthwhile question: can a predictor hold a specified error level while
retaining useful coverage across groups and shifts? Showing that confidence
thresholding improves accuracy among retained cases is routine. Compare against
calibrated confidence and established risk-control procedures; adaptive
risk-controlled selection has substantial prior art.

Pilot: UCI Diabetes 130-US Hospitals, 101,766 encounters, patient identifier,
~3 MB download / 18.3 MB CSV, CC BY 4.0. **Split by patient; define prediction
time before selecting features. The collection's name does not by itself provide
a hospital-held-out evaluation.**

Target **Machine Learning for Healthcare**. Week-one test: learned selection
against calibrated logistic/boosted-tree confidence and random deferral at
matched coverage; check gains persist **within predefined groups** rather than
arising from refusing an entire difficult group. 8–12 weeks bounded, 4–6 months
for clinical validation. **Backup, not preferred.**

### 4 · H&E foundation models — stop the biological branch

Retain only optional engineering work. An image model is necessary to map
morphology to an independently measured molecular endpoint from the same
specimen; it is **not** necessary to rediscover a diagnostic label assigned from
the same image.

- **UniToPatho**: 8,669 patches / 292 patients, locked 161/43/88 split; archive
  **48.37 GiB**, outside the repo in Google Drive. Streaming reduces local
  storage, not transfer. **Morphology labels, not the molecular endpoint.**
- **MHIST**: 3,152 tiles, 333 MB; annotations lack the patient/specimen
  molecular join and permission for derived modelling is unresolved. **Already
  closed for the molecular question.**
- **HTAN/CRDC**: 18 specimen-exact H&E–VCF matches are an **inventory**, not a
  sized callable cohort. Slide-storage cost and endpoint prevalence unverified.

Target Medical Image Analysis for a substantive advance with external
validation; frozen embeddings plus a linear head and ordinary abstention do not
clear it. Week-one test: inspect **only** the specimen join, callable endpoint
counts and feasible patient/site splits; if those fail, **stop before encoder
selection**. The present record already fails to establish the substrate.

## 6. What is funded, and what is not

**Funded this month:** the bounded BMC Bioinformatics methods submission, and
the first external-audit pilot.

**Not funded:** the active search for another CRC mechanism result; the new
predictor; the H&E ambitions — unless their week-one gates produce evidence that
changes this ranking.

**Kept as the longer programme:** simulation-based validation, with its unit of
evidence made **a verified failure of a stated scientific claim in external
software**.

### Consequences for open items in this repository

- `docs/prereg_chen_lesion_subtype.md` (item 3, unlocked) is **not funded this
  month**. It stays unlocked. Its §4 conclusion — that no compositional test is
  available on any obtainable label — is consistent with §3 here and is part of
  why the mechanistic paper is not pursued.
- `docs/DATA_HUNT_2026-09-10.md` §7's ordering is superseded: step 1 (recording
  the Blomain disagreement) is already done; steps 2–6 are CRC mechanism work
  and are **out of scope this month**.
- `docs/wes_subtype_plan.md` Synapse certification is **not chased**.

## 7. Verification status — what was checked here, on `ec93216`

`VERIFIED` means recomputed in this checkout on 2026-09-11.

| claim | status |
|---|---|
| `1b3e9b2` is a real commit on this history | **VERIFIED** |
| every file path the audit cites exists | **VERIFIED** (17/17) |
| 703 manifest entries | **VERIFIED** |
| 13.5 MB Parquet corpus; 287 tables at `ec93216` | **VERIFIED** (audit's 283 was at `1b3e9b2`) |
| PFI model is n=489 / 124 events, HR 0.9697 | **VERIFIED** |
| carcinoma limit +6.94, six genes 6.19–6.63, 1.07-fold | **VERIFIED** in `tests/test_paper_numbers.py:351`; the sign is positive |
| adenoma 0.27–3.13, 11.5-fold, from **cohort-averaged** inputs | **VERIFIED** in the same file at :369 — `_ratio_terms` takes block means |
| Xenium Colon v1 carries GUCA2A **and** MS4A12, no control | **VERIFIED** in `panel_coverage.parquet` |
| 14/20 percentile, 17/20 BCa, 0/20 Student-t miscalibrated; Student-t 3.00–5.93% | **VERIFIED** exactly |
| 6/8 cross-block at `lineage`, 1/8 at `best4` | **VERIFIED with a caveat**: 6/8 holds under `doubly_robust` and `normal` weighting, **5/8 under `tumour`**. The headline is weighting-dependent and should be quoted with its weighting. |
| `best4` 6/20 retain an intrinsic estimate under SMC vs 0/20 under KUL3 | **NOT VERIFIED HERE** — needs the estimability rule in code, not just the inputs table. Treat as the audit's finding. |
| 13 of 703 manifest entries present locally, 304.8 MB, hashes match | **NOT RE-VERIFIED** — the audit ran it; it is cheap to repeat. |
| intrinsic shares 0.715 (n=43) / 0.709 (n=42) | **NOT RE-VERIFIED** |
| 9.8 GB cluster ceiling | **NOT VERIFIABLE** from this machine (~87 GiB free locally) |
| ICBI-NB submitted PDF / submission status | **UNVERIFIED — the file is not in the checkout** |

**Note on the working tree.** This was written on branch `wmhs/build-fixes` at
`ec93216`, with uncommitted edits to `paper/wmhs/sections/{appendix,bench,
calibration}.tex` and untracked LaTeX build artifacts. Only this document and
the pointer in `HANDOFF.md` were committed.
