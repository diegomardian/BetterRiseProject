# Screening queue — prevalence audit

Built 2026-09-13 by the rule in `prereg_prevalence_audit.md` §2.4, **before any
repository's code was opened**. Committed separately from the screening
verdicts so that the order cannot be adjusted after the fact.

Queries as amended by A-1. Within each query, results are in the source's
default relevance order, recorded verbatim. The queue interleaves round-robin
(result #1 of every query, then #2 of every query, …); duplicates are dropped
at their second appearance and the slot passes to the next result of that
query.

## Raw query yields

### Q1b (S1, `synthetic control arm`) — 2 results

1. `OdyOSG/SyTrial` (1★) — "Synthetic Control Arm Construction Using OMOP CDM and Advanced Causal Inference"
2. `mondek5854/Aevolia` (0★) — "AI-native CRO for rare disease clinical trials — synthetic control arm generator…"

### Q2b (S1, `external control arm`) — 2 results

1. `mariannawicks-sys/nsclc-external-control-arm` (0★) — no description
2. `natsousa/ecasim` (0★) — "External Control Arm with Propensity Score Methods"

### Q3b (S1, `target trial emulation`) — 8 results

1. `bldestavola/TTE-Short-Course` (20★)
2. `kathoffman/steroids-trial-emulation` (73★)
3. `Alvin-RB/antipsychotics_tte_cprd` (8★)
4. `cjyarnell/TargetTrial_ThresholdsForIMV` (12★)
5. `surajraj99/Corticosteroids-in-Patients-with-Critical-Illness` (8★)
6. `nakashimatakaya/vasopressin-timing-target-trial-emulation` (3★)
7. `minji79/TargetTrialEmulation` (0★)
8. `UQualityCausal/TargetTrialEmulation` (1★)

### Q4 (S1, `plasmode simulation`, **unamended**) — 8 results

1. `cran/Plasmode` (2★) — read-only CRAN mirror of the `Plasmode` R package
2. `inab-certh/PlasmodeSimulation` (0★)
3. `PamelaShaw/PlasmodeSimulation` (1★) — "Paper exploring proper methods for plasmode simulation to study a treatment effect"
4. `smlemons/plasmode` (0★)
5. `carolynlou/hcct` (0★)
6. `detal9/LongitudinalPlasmode` (1★)
7. `ehsanx/hdps-tmle-sim` (3★)
8. `MariekeStolte/CompPlasmodeParamHDClassif` (0★)

### Q5b (S1, `causal estimator simulation`) — **0 results**

Recorded as a zero-yield query. Not replaced again (A-1).

### Q6b (S2, `"true_ate"`) — 15 file hits / 10 distinct repositories

1. `brycewang-stanford/Auto-Empirical-Research-Skills` (`benchmark/check_benchmark.py`, `benchmark/tasks/qte-recovery.toml`, `benchmark/tasks/cate-recovery.toml`)
2. `Netflix-Skunkworks/oci-agent` (`evals/smoketest/eval.py`)
3. `lfiaschi/bayesian-autoresearcher` (`problems/ihdp/prepare.py`, `problems/lalonde/prepare.py`, `problems/twins/prepare.py`, `runner.py`)
4. `brycewang-stanford/StatsPAI` (`src/statspai/utils/dgp.py`)
5. `KDL-umass/papers` (`causal-eval-rct-icml-2021/merge-nn-results.R`)
6. `pymc-labs/pathmc` (`docs/examples/01-foundations/causal_do_operator.qmd`)
7. `Hemanthyadv/causal-inference-platform` (`API.md`)
8. `sandroama/pricing-lab` (`docs/API.md`)
9. `CausalAILab/NeuralCausalModels` (`src/run/pipeline.py`)
10. `facebookresearch/data_decomposition` (`experiments/utils.py`)

### Q7b (S2, `"true_effect"`) — 15 file hits / 8 distinct repositories

1. `rohitg00/ai-engineering-from-scratch` (`phases/01-math-foundations/15-statistics-for-ml/code/statistics.py`)
2. `easystats/bayestestR` (`vignettes/web_only/indicesEstimationComparison.Rmd`)
3. `brycewang-stanford/Auto-Empirical-Research-Skills` — duplicate of Q6b #1
4. `easystats/easystats` (`publications/ludecke_----_priors/data/#01 - Simulations.R`)
5. `Netflix-Skunkworks/oci-agent` — duplicate of Q6b #2
6. `brycewang-stanford/StatsPAI` — duplicate of Q6b #4
7. `igerber/diff-diff` (`diff_diff/prep_dgp.py`, `diff_diff/power.py`)
8. `py-why/causaltune` (`causaltune/datasets.py`)

### Q8 (S3, web, **unamended**) — papers, then their announced code

Top results and their code routes:

1. arXiv 2504.11740, *A cautionary note for plasmode simulation studies in the
   setting of causal inference* (Shaw, Gruber, Williamson, Desai, Shortreed,
   Krakauer, Nelson, van der Laan). **Abstract page announces no code URL.**
   The repository `PamelaShaw/PlasmodeSimulation` is reached independently by
   Q4 #3 and enters the queue there.
2. `causl` R package Plasmode vignette (mrcieu r-universe / Oxford) →
   source repository `rje42/causl`.
3. BMC Med Res Methodol 10.1186/s12874-023-02062-9, *Longitudinal plasmode
   algorithms…* → `detal9/LongitudinalPlasmode`, already Q4 #6.

### Q9 (S3, web, **unamended**) — **0 examinable candidates**

1. arXiv 2507.16048, *Evaluating virtual-control-augmented trials for
   reproducing treatment effect from original RCTs* (Fernandes, Porcher, Tran,
   Petit). Full text checked for a code/data availability statement:
   **NO CODE URL**. Fails I1 before any screening of content.
2.–8. Commercial and editorial pages (Medidata, IntuitionLabs, Nova In Silico,
   Quibim, TalkingHealthTech), an unrelated IEEE robotics paper on a "virtual
   arm" rehabilitation model, and a medRxiv case study with no code URL. None
   is a code artifact.

This is itself a finding and is carried into the coverage limits: the query
aimed squarely at the paper's own domain — virtual control arms — returned no
public code at all.

### Q10 (S3, web, **unamended**) — GitHub links in result order

1. `Zhengxian-Fan/target-trial-emulation`
2. `kathoffman/steroids-trial-emulation` — duplicate of Q3b #2
3. `didierbrassard/NuAge_protocol`
4. `Causal-LDA/TrialEmulation`
5. `andrew-yiu/Predictive-causal-inference`
6. `dsgelab/trial_emulations_genetics`

---

## The queue

| # | repository / artifact | query |
|---|---|---|
| 1 | `OdyOSG/SyTrial` | Q1b |
| 2 | `mariannawicks-sys/nsclc-external-control-arm` | Q2b |
| 3 | `bldestavola/TTE-Short-Course` | Q3b |
| 4 | `cran/Plasmode` | Q4 |
| 5 | `brycewang-stanford/Auto-Empirical-Research-Skills` | Q6b |
| 6 | `rohitg00/ai-engineering-from-scratch` | Q7b |
| 7 | `rje42/causl` | Q8 |
| 8 | `Zhengxian-Fan/target-trial-emulation` | Q10 |
| 9 | `mondek5854/Aevolia` | Q1b |
| 10 | `natsousa/ecasim` | Q2b |
| 11 | `kathoffman/steroids-trial-emulation` | Q3b |
| 12 | `inab-certh/PlasmodeSimulation` | Q4 |
| 13 | `Netflix-Skunkworks/oci-agent` | Q6b |
| 14 | `easystats/bayestestR` | Q7b |
| 15 | `detal9/LongitudinalPlasmode` | Q8 |
| 16 | `didierbrassard/NuAge_protocol` | Q10 |
| 17 | `Alvin-RB/antipsychotics_tte_cprd` | Q3b |
| 18 | `PamelaShaw/PlasmodeSimulation` | Q4 |
| 19 | `lfiaschi/bayesian-autoresearcher` | Q6b |
| 20 | `easystats/easystats` | Q7b |
| 21 | `Causal-LDA/TrialEmulation` | Q10 |
| 22 | `cjyarnell/TargetTrial_ThresholdsForIMV` | Q3b |
| 23 | `smlemons/plasmode` | Q4 |
| 24 | `brycewang-stanford/StatsPAI` | Q6b |
| 25 | `igerber/diff-diff` | Q7b |
| 26 | `andrew-yiu/Predictive-causal-inference` | Q10 |
| 27 | `surajraj99/Corticosteroids-in-Patients-with-Critical-Illness` | Q3b |
| 28 | `carolynlou/hcct` | Q4 |
| 29 | `KDL-umass/papers` | Q6b |
| 30 | `py-why/causaltune` | Q7b |
| 31 | `dsgelab/trial_emulations_genetics` | Q10 |
| 32 | `nakashimatakaya/vasopressin-timing-target-trial-emulation` | Q3b |
| 33 | `ehsanx/hdps-tmle-sim` | Q4 |
| 34 | `pymc-labs/pathmc` | Q6b |
| 35 | `minji79/TargetTrialEmulation` | Q3b |
| 36 | `MariekeStolte/CompPlasmodeParamHDClassif` | Q4 |
| 37 | `Hemanthyadv/causal-inference-platform` | Q6b |
| 38 | `UQualityCausal/TargetTrialEmulation` | Q3b |
| 39 | `CausalAILab/NeuralCausalModels` | Q6b |
| 40 | `sandroama/pricing-lab` | Q6b |
| 41 | `facebookresearch/data_decomposition` | Q6b |

The queue is walked in this order. Screening verdicts for every entry walked
are recorded in `prevalence_audit_result.md`, eligible or not.
