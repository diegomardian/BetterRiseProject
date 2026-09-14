"""Emit the prevalence-audit tables as versioned parquet.

The audit itself is human work: the classification rule is pre-registered in
``docs/prereg_prevalence_audit.md``, the screening queue in
``docs/prevalence_audit_queue.md``, and the verdicts were reached by reading
each repository's code. This module is only the recorder — it puts those
verdicts on disk in repo shape (``results/{date}_{sha7}/`` plus a provenance
sidecar) so that a reader can join them to a commit.

Two tables:

``prevalence_audit_verdicts``
    One row per repository **examined and classified** under the
    pre-registered rule. Carries the estimator, the reference, where each is
    computed, the verdict and its basis.

``prevalence_audit_screening``
    One row per queue entry **walked but not examined**, with the inclusion
    criterion it failed. The prereg (§2.4 step 3) requires a recorded verdict
    for every entry walked, eligible or not — a screening trail that only
    lists the hits is the shape a cherry-picked audit takes.

The seed is ``0`` and means it: no sampling, no RNG. Invariant 10 asks for a
seed to be stated, not invented, and ``0`` is the honest answer for a table
whose only randomness lives in GitHub's relevance ranking, which is recorded
by query and position instead.
"""

from __future__ import annotations

import pandas as pd

from src.common.io import write_versioned_table

#: Deterministic table. See the module docstring.
SEED = 0

SEARCH_DATE = "2026-09-13"

# ---------------------------------------------------------------------------
# Examined and classified
# ---------------------------------------------------------------------------

VERDICTS: list[dict] = [
    {
        "queue_position": 5,
        "query_id": "Q6b",
        "name": "brycewang-stanford/Auto-Empirical-Research-Skills",
        "url": "https://github.com/brycewang-stanford/Auto-Empirical-Research-Skills",
        "commit_examined": "d31f17f0e26a2086a3efa6d982d7fc089fed33ce",
        "language": "Python",
        "domain": "benchmark suite grading agent-written causal-estimator pipelines",
        "estimator_description": (
            "cross-fitted partialling-out DML: residualise y and d on the four "
            "controls with out-of-fold linear nuisances over 4 folds, then "
            "sum(rd*ry)/sum(rd*rd)"
        ),
        "estimator_location": "benchmark/lib/dml.py:103-119 (wired at benchmark/reference_pipeline.py:186)",
        "reference_description": (
            "full-sample OLS of y on [1, d, x1..x4], coefficient on d, "
            "recomputed by the checker from the same CSV"
        ),
        "reference_location": "benchmark/lib/dml.py:89-93 (consumed at benchmark/check_benchmark.py:442)",
        "reference_provenance": "3c",
        "step4_first_difference": "4.4 nuisance estimation (out-of-fold vs full-sample), then 4.2 weighting",
        "verdict": "NO",
        "verdict_basis": "4.4",
        "code_executed": "yes",
        "max_residual": 3.975e-14,
        "residual_note": (
            "numerical equality on the shipped draw; resolved algebraically, "
            "not functional reuse — see justification"
        ),
        "resolves_unclear": "",
        "justification": (
            "This is the only repository in the audit whose reference is "
            "recomputed from the realised simulated draw rather than from a "
            "generating parameter, so it is the only one to reach Step 4. "
            "Running the shipped code gives true_theta = 1.4999999999999607 and "
            "dml_theta = 1.5000000000000004, a residual of 4.0e-14: numerical "
            "equality, which the paper's own screen says must be investigated "
            "algebraically rather than read as reuse. The investigation shows "
            "it is not reuse. generate() builds y = THETA*d + _g(x) with _g "
            "exactly linear in the controls and no noise (benchmark/lib/dml.py:"
            "48-58), so for any OLS fit of y on the controls the coefficient "
            "vector is exactly THETA*beta_d + beta_g; the residuals therefore "
            "satisfy ry = THETA*rd pointwise, in-fold and out-of-fold alike, "
            "and every ratio sum(rd*ry)/sum(rd*rd) returns THETA regardless of "
            "the fold partition. The equality is a property of this noiseless "
            "exactly-linear data file, not of the two functionals. Adding "
            "Gaussian outcome noise to the same rows and changing nothing else "
            "separates them immediately: true_theta = 1.599471, dml_theta = "
            "1.474678, residual 0.1248. Cross-fitting versus a full-sample fit "
            "is precisely the Step 4.4 difference. Verdict NO. Two further "
            "observations, recorded because they bear on how the repository is "
            "read rather than on its verdict. First, the benchmark separates "
            "its checks by name: the golds titled 'recovers-true' compare the "
            "candidate's number to true_theta, while those titled 'honest-*' "
            "compare it to the checker's own run of the same estimator and are "
            "described in-repo as ensuring 'no fabricated numbers' "
            "(benchmark/tasks/dml-recovery.toml:54-62). The second family is a "
            "reproduction check with a stated anti-fabrication purpose, and the "
            "paper's §blind draws exactly this distinction — 'equality "
            "identifies redundancy' — so naming them apart is the opposite of "
            "conflating them. A reader who nonetheless took an 'honest-*' gold "
            "as the validation would reach YES on that gold; the prereg "
            "classifies the primary experiment, which is the recovery gold. "
            "Second, the estimator under test is normally the candidate "
            "pipeline's own code, not in this repository; the shipped "
            "reference-dml candidate exists so the benchmark runs out of the "
            "box (benchmark/reference_pipeline.py:2-7)."
        ),
    },
    {
        "queue_position": 7,
        "query_id": "Q8",
        "name": "rje42/causl",
        "url": "https://github.com/rje42/causl",
        "commit_examined": "0659ed8fff846d0df323353bdbc0c524326585ca",
        "language": "R",
        "domain": "R package for specifying, simulating from and fitting causal models (frugal parametrisation)",
        "estimator_description": (
            "maximum-likelihood fit of the frugal/marginal causal model, "
            "fit_causl(dat, family=c(1,1,1)), taking out$pars$y$beta as the "
            "estimate of the causal parameter"
        ),
        "estimator_location": "tests/testthat/test-fit.R:9 and :14 (also vignettes/Introduction.Rmd:133)",
        "reference_description": (
            "pars$y$beta = c(0, 0.5), the generating coefficient vector written "
            "down before any data existed and handed to causalSamp() as the DGP"
        ),
        "reference_location": "tests/testthat/test-fit.R:4 (compared at :18-19)",
        "reference_provenance": "3a",
        "step4_first_difference": "not reached — Step 3 decides",
        "verdict": "NO",
        "verdict_basis": "3a",
        "code_executed": "no",
        "max_residual": float("nan"),
        "residual_note": "",
        "resolves_unclear": "",
        "justification": (
            "The generating parameters are fixed at tests/testthat/test-fit.R:2-5 "
            "before any draw: pars <- list(z=..., x=..., y=list(beta=c(0,0.5), "
            "phi=0.5), cop=...). Data are drawn with causalSamp(1e2, par=pars, "
            "family=fam) and the model is fitted with fit_causl(dat, "
            "family=c(1,1,1)). The check is a Wald-type statistic against the "
            "generating vector: expect_lt(sum((out$pars$y$beta - pars$y$beta)^2 / "
            "out$pars$y$beta_sandwich^2), qchisq(0.99, df=2)). The reference is "
            "the requested parameter and depends on no realised quantity, so it "
            "cannot be a functional of the draw. The vignette states the same "
            "comparison in prose: 'the parameter estimates are all well within "
            "two standard errors of the values we initially put in' "
            "(vignettes/Introduction.Rmd:136-137). This is the correct practice "
            "the paper's analysis recommends, and the test can fail."
        ),
    },
    {
        "queue_position": 8,
        "query_id": "Q10",
        "name": "Zhengxian-Fan/target-trial-emulation",
        "url": "https://github.com/Zhengxian-Fan/target-trial-emulation",
        "commit_examined": "35983ac2d7e3b855a68152332fb2e6ff8667949e",
        "language": "Python",
        "domain": "target trial emulation on CPRD, with a plasmode-style benchmarking arm",
        "estimator_description": (
            "three estimators of the treatment effect on survival: 1:1 "
            "propensity-score matching with a caliper then a cluster-robust Cox "
            "fit; stabilised IPTW Cox over ten trimming thresholds; TMLE for the "
            "risk ratio with an influence-curve variance — pooled across "
            "imputations by Rubin's rules"
        ),
        "estimator_location": "estimate_treatment_effects.py:73-152 (pooled at :178-195)",
        "reference_description": (
            "CONFIG['theta_treat'] = -0.5, the log hazard ratio multiplied into "
            "the hazard when the simulated outcome is generated"
        ),
        "reference_location": "simulate.py:14 (entering the hazard at :68)",
        "reference_provenance": "3a",
        "step4_first_difference": "not reached — Step 3 decides",
        "verdict": "NO",
        "verdict_basis": "3a",
        "code_executed": "no",
        "max_residual": float("nan"),
        "residual_note": "",
        "resolves_unclear": "",
        "justification": (
            "simulate.py keeps the real cohort's covariates and regenerates "
            "treatment and outcome: treatment is drawn from a confounded "
            "propensity (simulate.py:57) and the hazard carries "
            "exp(CONFIG['theta_treat'] * treatment) with theta_treat = -0.5 "
            "(simulate.py:14, :68). The README describes the intended use "
            "exactly: 'simulates Treatment and Outcome based on a known Hazard "
            "Ratio ... to check if it recovers the true effect size'. The "
            "reference is therefore a constant in the generator's configuration "
            "and no functional of the draw can equal it. Recorded as a limit on "
            "the evidence rather than on the verdict: the comparison itself is "
            "not automated. estimate_treatment_effects.py prints pooled hazard "
            "and risk ratios (:181, :188, :195) and the reader compares them to "
            "exp(-0.5) by eye; no line of code computes a bias, a residual or a "
            "tolerance. That leaves nothing for an equality screen to reuse, and "
            "it also means the recovery check is not run by the test suite."
        ),
    },
    {
        "queue_position": 10,
        "query_id": "Q2b",
        "name": "natsousa/ecasim",
        "url": "https://github.com/natsousa/ecasim",
        "commit_examined": "53a201ebac1ebca39ef63fec1230c737d96951de",
        "language": "R",
        "domain": "external control arm for a single-arm oncology trial, propensity-score methods",
        "estimator_description": (
            "Cox proportional-hazards fit for the treatment term, run "
            "unweighted, after IPTW/stabilised weights, and after 2:1 "
            "propensity-score matching"
        ),
        "estimator_location": "R/analyze.R (fit_cox), exercised at tests/testthat/test-core.R:35-42 and analysis/report.qmd",
        "reference_description": (
            "log_hr_true <- log(0.7), the treatment coefficient in the "
            "exponential survival DGP"
        ),
        "reference_location": "R/simulate.R:61 (entering the linear predictor at :62-63)",
        "reference_provenance": "3a",
        "step4_first_difference": "not reached — Step 3 decides",
        "verdict": "NO",
        "verdict_basis": "3a",
        "code_executed": "no",
        "max_residual": float("nan"),
        "residual_note": "",
        "resolves_unclear": "",
        "justification": (
            "simulate_rwd() draws covariates, assigns the trial arm by ranking a "
            "latent eligibility score so that the arms are deliberately "
            "confounded (R/simulate.R:54-58), and generates exponential event "
            "times whose rate carries log_hr_true * treatment with log_hr_true "
            "<- log(0.7) (R/simulate.R:60-66). The comparison appears in the "
            "analysis report as prose: 'the estimate is closer to the true "
            "generative HR of 0.70' (analysis/report.qmd:195). The reference is "
            "a constant fixed by the generator, so Step 3(a) decides and Step 4 "
            "is not reached. Recorded as a limit on the evidence: the package's "
            "own tests check structural properties — that the cohort is "
            "confounded in the intended direction (test-core.R:10-15), that "
            "matching reduces the age imbalance (:26-33), that fit_cox returns a "
            "positive hazard ratio for the treatment term (:35-42) — and none "
            "asserts recovery of 0.7 within a tolerance. The recovery claim is "
            "therefore made in the report rather than by the test suite."
        ),
    },
    {
        "queue_position": 18,
        "query_id": "Q4",
        "name": "PamelaShaw/PlasmodeSimulation",
        "url": "https://github.com/PamelaShaw/PlasmodeSimulation",
        "commit_examined": "7ae26ee73509eaa0da02c033da76bd333fea6daf",
        "language": "R",
        "domain": (
            "plasmode simulation for causal inference; code for arXiv:2504.11740, "
            "'A cautionary note for plasmode simulation studies in the setting of "
            "causal inference'"
        ),
        "estimator_description": (
            "eleven estimators of the ATE (and of RR, EY0, EY1, CATE) applied to "
            "each plasmode replicate: unadjusted, matching, IPTW bounded and "
            "unbounded, TMLE bounded and unbounded, outcome-regression, "
            "PS-covariate, and PS stratification at 5, 10 and 20 strata"
        ),
        "estimator_location": "plasmode_simulations/00_plasmode_utils.R:257-263 (one_sim), driven at :307-320",
        "reference_description": (
            "psi0.sim = ey1 - ey0, where ey1 and ey0 standardise the fitted "
            "outcome model over the ORIGINAL cohort's covariates with treatment "
            "forced to 1 and to 0; likewise RR0 = ey1/ey0 and CATE0 = the "
            "generating coefficient itself"
        ),
        "reference_location": "plasmode_simulations/00_plasmode_utils.R:278-288 (compared at :340-349 via calcResults, defined :213-221)",
        "reference_provenance": "3a and 3b",
        "step4_first_difference": "not reached — Step 3 decides",
        "verdict": "NO",
        "verdict_basis": "3a, 3b",
        "code_executed": "no",
        "max_residual": float("nan"),
        "residual_note": "",
        "resolves_unclear": "",
        "justification": (
            "The reference is computed once, before the replication loop begins "
            "and therefore before any draw exists: ey1 <- mean(pred.Y(A = "
            "rep(1, nrow(original_data)), W = original_data[,-(1:2)], "
            "outcome_model = outcome_model)) and the matching ey0 with A = 0, "
            "then psi0.sim <- ey1 - ey0 (00_plasmode_utils.R:278-286). Three "
            "things put this outside Step 4 at once. It is evaluated on the "
            "original cohort's covariates, all nrow(original_data) of them, not "
            "on the resampled draw of size n. It forces treatment to 1 and to 0 "
            "for every record, so it reads counterfactual predictions the "
            "estimator never sees. And it is computed from the fitted "
            "outcome_model that defines the DGP, not from the simulated "
            "outcomes. CATE0 <- coefficients(outcome_model)[2] (:288) is the "
            "generating coefficient outright. The comparison is an explicit "
            "bias calculation over B replicates — calcResults returns bias = "
            "colMeans(est - psi0) together with median bias, percent bias and "
            "bias-to-SE (:213-221) — and the replicate estimates enter it as a "
            "B-by-11 matrix (:340-349), so the residual is neither zero nor "
            "constrained to be. This is the practice the paper's analysis "
            "recommends, in a repository whose paper is about how to define the "
            "plasmode truth correctly. Coverage note: the underlying cohort is "
            "not shareable, so the code was read and not run."
        ),
    },
    {
        "queue_position": 24,
        "query_id": "Q6b",
        "name": "brycewang-stanford/StatsPAI",
        "url": "https://github.com/brycewang-stanford/StatsPAI",
        "commit_examined": "dd4640e26b5fe32a791aa5c95af2589010ca549d",
        "language": "Python",
        "domain": "benchmark of causal estimators (RCT, DiD, RDD, IV, …) inside a statistics toolkit",
        "estimator_description": (
            "the toolkit's estimator for each design, run once on the full "
            "simulated dataset to give a point estimate"
        ),
        "estimator_location": "src/statspai/smart/benchmark.py:238-251",
        "reference_description": (
            "spec['true_effect'] = 0.5, the same scalar passed to the generator "
            "as its effect argument"
        ),
        "reference_location": "src/statspai/smart/benchmark.py:48-104 (spec table); generator at src/statspai/utils/dgp.py:95, :101, :119",
        "reference_provenance": "3a",
        "step4_first_difference": "not reached — Step 3 decides",
        "verdict": "NO",
        "verdict_basis": "3a",
        "code_executed": "no",
        "max_residual": float("nan"),
        "residual_note": "",
        "resolves_unclear": "",
        "justification": (
            "Each benchmark spec carries dgp_kwargs with an effect of 0.5 and a "
            "true_effect field of 0.5 (benchmark.py:48-62), and the generator "
            "writes the same scalar onto the frame: unit_effects <- "
            "np.full(n_units, effect), y = unit_fe + time_fe + unit_effects*d + "
            "noise, df.attrs['true_effect'] = effect (dgp.py:95, :101, :119). "
            "The comparison is bias = point - spec['true_effect'] "
            "(benchmark.py:251-252), aggregated as a mean absolute bias across "
            "replicates (:293). The reference is the requested parameter, so "
            "Step 3(a) decides. One detail worth stating because it is the "
            "distinction the paper's analysis turns on: under heterogeneous=True "
            "the unit effects are effect + N(0, 0.3) (dgp.py:93), so the "
            "realised sample's average effect differs from 0.5 by sampling "
            "noise, and the reported bias mixes estimator error with that "
            "realised-versus-requested gap. That is a property of using the "
            "requested parameter as the reference, and it is the opposite "
            "failure from the one this audit looks for."
        ),
    },
    {
        "queue_position": 25,
        "query_id": "Q7b",
        "name": "igerber/diff-diff",
        "url": "https://github.com/igerber/diff-diff",
        "commit_examined": "ff03d08c6762306b7b69ad388c3ea7bb4cb11c4a",
        "language": "Python",
        "domain": "difference-in-differences library (Callaway–Sant'Anna, synthetic DiD, honest DiD, event studies)",
        "estimator_description": (
            "Callaway–Sant'Anna group-time ATT(g,t), its event-study and overall "
            "aggregations, under doubly-robust, outcome-regression and IPW "
            "estimands"
        ),
        "estimator_location": "tests/test_csdid_ported.py:96, :155, :168, :193-200 and throughout",
        "reference_description": (
            "hard-coded DGP constants — ATT = 1 in the Callaway–Sant'Anna "
            "generator, 3 (with 1.5x dose response) in the two-group generator, "
            "and event-study ATT(e=k) = k by construction"
        ),
        "reference_location": "tests/helpers/csdid_dgp.py:283, :310-323; asserted at tests/test_csdid_ported.py:96, :342, :376",
        "reference_provenance": "3a",
        "step4_first_difference": "not reached — Step 3 decides",
        "verdict": "NO",
        "verdict_basis": "3a",
        "code_executed": "no",
        "max_residual": float("nan"),
        "residual_note": "",
        "resolves_unclear": "",
        "justification": (
            "The generator builds both potential-outcome paths and selects the "
            "observed one — Y_treated[:, 2] += att_rand with att_rand ~ N(3, 1), "
            "Y_treated[:, 3] += 1.5 * att_rand, then Y_obs = where(is_treated_now, "
            "Y_treated, Y_untreated) (tests/helpers/csdid_dgp.py:311-330). The "
            "assertions compare the estimate to the DGP constant with an "
            "explicit tolerance, e.g. assert abs(gt_effects[first_key]['effect'] "
            "- 1.0) < 0.5 (tests/test_csdid_ported.py:96) and assert "
            "abs(att_e2 - 3.0) < 1.0 (:342). The reference is the requested "
            "parameter; Step 3(a) decides. Noted as a secondary observation: the "
            "suite labels a second tier that asserts 'strict tolerance against "
            "R's exact output' (:4-8). That is an agreement check against an "
            "independent implementation rather than against a simulated truth, "
            "and the file separates the two tiers by name in its own docstring, "
            "so the recovery check and the port check are not conflated."
        ),
    },
]

# ---------------------------------------------------------------------------
# Walked but not examined
# ---------------------------------------------------------------------------

SCREENING: list[dict] = [
    (1, "Q1b", "OdyOSG/SyTrial", "https://github.com/OdyOSG/SyTrial", "I4",
     "R package for synthetic control arms (TMLE, G-computation, IPTW). Its "
     "examples simulate a cohort and run the estimators, but a grep for truth, "
     "true_*, bias, rmse and coverage across R/, README.md and the .Rd files "
     "returns nothing: no quantity is compared to a target, so there is no "
     "validation to audit."),
    (2, "Q2b", "mariannawicks-sys/nsclc-external-control-arm",
     "https://github.com/mariannawicks-sys/nsclc-external-control-arm", "E3",
     "The four analysis scripts are placeholders: 01_data_sim.R contains the 16 "
     "bytes '01_data_sim.R', and the other three likewise contain only their own "
     "path. The simulated CSV, the figures and the Cox results are committed but "
     "the producing code is absent, so neither estimator nor reference can be "
     "located."),
    (3, "Q3b", "bldestavola/TTE-Short-Course",
     "https://github.com/bldestavola/TTE-Short-Course", "I2",
     "Short-course teaching materials: lecture PDFs, Stata do-files and an R "
     "practical operating on fixed .dta files. No data-generating process, and "
     "no target value to recover."),
    (4, "Q4", "cran/Plasmode", "https://github.com/cran/Plasmode", "I3, I4",
     "The Plasmode R package is a generator, not an evaluation: it returns "
     "simulated datasets plus RR and RD but applies no estimator and makes no "
     "comparison. Recorded because the reference it hands downstream users is "
     "directly relevant — RR[sim] <- mean(p_1)/mean(p_0) where p_1 and p_0 are "
     "computed from the true coefficient vector bnew on the sampled design "
     "matrix with the exposure column set to 1 and to 0 (R/PlasmodeBin.R:174-180). "
     "That never touches the realised outcomes ynew, so a downstream estimator "
     "cannot be the same functional of it."),
    (6, "Q7b", "rohitg00/ai-engineering-from-scratch",
     "https://github.com/rohitg00/ai-engineering-from-scratch", "I5, E2",
     "A statistics-for-ML teaching module inside a general AI-engineering "
     "curriculum. Its ab_test_simulator generates two groups and runs a t-test "
     "but is not a trial, control-arm, target-trial or plasmode simulation. Had "
     "it been included, the verdict would have been NO: true_effect is the "
     "generator's own argument (statistics.py:396, :402, :412)."),
    (9, "Q1b", "mondek5854/Aevolia", "https://github.com/mondek5854/Aevolia", "I1, I2",
     "Two HTML files and a README describing a synthetic-control-arm service. No "
     "source code, no simulation."),
    (11, "Q3b", "kathoffman/steroids-trial-emulation",
     "https://github.com/kathoffman/steroids-trial-emulation", "I2, I4",
     "Target-trial-emulation tutorial with sequentially doubly robust estimation. "
     "The README states the demo data are simulated (n=500), but the generating "
     "code is not in the repository — code/analysis.R loads data/dat_demo.rds — "
     "and no target value is compared to."),
    (12, "Q4", "inab-certh/PlasmodeSimulation",
     "https://github.com/inab-certh/PlasmodeSimulation", "I3, I4",
     "An OMOP-CDM plasmode sample generator: it fits Poisson models on real "
     "cohorts, modifies coefficients to inject known rate ratios, and emits "
     "samples. Its seven exported functions are all cohort extraction, model "
     "fitting and sample generation; no estimator is applied and nothing is "
     "compared."),
    (13, "Q6b", "Netflix-Skunkworks/oci-agent",
     "https://github.com/Netflix-Skunkworks/oci-agent", "I2",
     "Evaluates causal estimators against the ACIC benchmark, but the "
     "data-generating process is the ACIC distribution and is not in this "
     "repository — evals/smoketest/eval.py reads zymu_{r}.csv. Recorded because "
     "its reference is the correct practice in the purest form: true_ate = "
     "float(ite.mean()) with ite = zymu['mu1'] - zymu['mu0'] (eval.py:34-38), "
     "i.e. the realised sample mean of both potential outcomes, which the "
     "estimator never sees. Excluding it drops a NO, which shrinks the audit "
     "rather than flattering it."),
    (14, "Q7b", "easystats/bayestestR", "https://github.com/easystats/bayestestR", "I3, I5",
     "The matching vignette compares Bayesian point-estimate indices (median, "
     "mean, MAP) against the simulated true_effect of a regression coefficient. "
     "It is an index-estimation comparison, not a causal contrast, and not in "
     "any of the audit's domains."),
    (15, "Q8", "detal9/LongitudinalPlasmode",
     "https://github.com/detal9/LongitudinalPlasmode", "I3, I4",
     "Four generator scripts for longitudinal plasmode simulation, including "
     "counterfactual variants that force a treatment regime via the A.fixed "
     "matrix (plasmode_counterfactual_parametric.R:188-192, :227-230). No "
     "estimator is applied and no comparison is made inside the repository; the "
     "forced-regime datasets are the raw material a user would compute a truth "
     "from."),
    (16, "Q10", "didierbrassard/NuAge_protocol",
     "https://github.com/didierbrassard/NuAge_protocol", "I2, I4",
     "A published protocol for a target trial emulation. The 'simulation' data "
     "are Health Canada's simulated composite diets used as exemplar inputs, not "
     "a data-generating process with a known effect, and no recovery is checked."),
    (17, "Q3b", "Alvin-RB/antipsychotics_tte_cprd",
     "https://github.com/Alvin-RB/antipsychotics_tte_cprd", "I2",
     "Analysis code for a published CPRD target trial emulation. Real data "
     "throughout; no data-generating process."),
    (19, "Q6b", "lfiaschi/bayesian-autoresearcher",
     "https://github.com/lfiaschi/bayesian-autoresearcher", "I2",
     "All four problems (ihdp, lalonde, nhefs, twins) read fixed public "
     "datasets; no data-generating process in the repository. Same disposition "
     "as oci-agent, and the same correct reference: true_ate = float((df['mu1'] "
     "- df['mu0']).mean()) (problems/ihdp/prepare.py:26)."),
    (20, "Q7b", "easystats/easystats", "https://github.com/easystats/easystats", "I3, I5",
     "The matching file is a Bayesian-priors simulation in a publications "
     "directory. No causal contrast; outside the audit's domains."),
    (21, "Q10", "Causal-LDA/TrialEmulation",
     "https://github.com/Causal-LDA/TrialEmulation", "I4",
     "The TrialEmulation R package does contain a data-generating process — "
     "data_gen_censored(), implementing Young and Tchetgen Tchetgen (2014), with "
     "all_treat and all_control switches that force a regime 'as if in an RCT' "
     "(R/data_simulation.R:9-11). But the simulated data are used only for "
     "software tests and vignette demonstrations; no test or vignette compares "
     "an estimate to a target value, so there is no simulation-based validation "
     "in the repository to classify."),
    (22, "Q3b", "cjyarnell/TargetTrial_ThresholdsForIMV",
     "https://github.com/cjyarnell/TargetTrial_ThresholdsForIMV", "I2",
     "Cohort construction and analysis for a target trial emulation on MIMIC and "
     "AmsterdamUMCdb. The files named *_validate* validate the confounder and "
     "outcome models, not the recovery of a simulated causal effect; there is no "
     "data-generating process with a known effect."),
    (23, "Q4", "smlemons/plasmode", "https://github.com/smlemons/plasmode", "I3, I5",
     "Ships as a source tarball whose contents are plasmodeItem, plasmodePerson, "
     "plasmodeData and fnr — 'plasmode' in the psychometric item-response sense, "
     "not the causal-inference one. No treatment effect, no estimator under "
     "test."),
    (26, "Q10", "andrew-yiu/Predictive-causal-inference",
     "https://github.com/andrew-yiu/Predictive-causal-inference", "I2",
     "R code for causal predictive inference applied to the Cattaneo birthweight "
     "dataset, downloaded from Stata Press. No data-generating process in the "
     "repository."),
]


def build_verdicts() -> pd.DataFrame:
    df = pd.DataFrame(VERDICTS)
    df.insert(0, "search_date", SEARCH_DATE)
    return df.sort_values("queue_position").reset_index(drop=True)


def build_screening() -> pd.DataFrame:
    df = pd.DataFrame(
        SCREENING,
        columns=["queue_position", "query_id", "name", "url", "criterion_failed", "note"],
    )
    df.insert(0, "search_date", SEARCH_DATE)
    df["screening_verdict"] = "screened, not examined"
    return df.sort_values("queue_position").reset_index(drop=True)


def main(allow_dirty: bool = False) -> None:
    verdicts = build_verdicts()
    screening = build_screening()

    counts = verdicts["verdict"].value_counts().to_dict()
    shared = {
        "audit_search_date": SEARCH_DATE,
        "prereg": "docs/prereg_prevalence_audit.md",
        "prereg_amendments": ["A-1 (2026-09-13): relaxed six zero-yield query strings"],
        "queue": "docs/prevalence_audit_queue.md",
        "queue_length": 41,
        "queue_walked_through_position": 26,
        "n_examined": int(len(verdicts)),
        "n_screened_not_examined": int(len(screening)),
        "verdict_counts": counts,
        "n_reaching_step4": int((verdicts["reference_provenance"] == "3c").sum()),
        "single_rater": True,
        "prevalence_rate_computed": False,
    }

    write_versioned_table(
        verdicts,
        name="prevalence_audit_verdicts",
        seed=SEED,
        notes=(
            "Repositories examined and classified under the pre-registered rule "
            "in docs/prereg_prevalence_audit.md §3. One row per repository."
        ),
        extra_meta=shared | {"table_role": "verdicts"},
        allow_dirty=allow_dirty,
    )
    write_versioned_table(
        screening,
        name="prevalence_audit_screening",
        seed=SEED,
        notes=(
            "Queue entries walked but not examined, with the inclusion criterion "
            "each failed. Required by prereg §2.4 step 3."
        ),
        extra_meta=shared | {"table_role": "screening"},
        allow_dirty=allow_dirty,
    )


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--allow-dirty", action="store_true")
    main(**vars(p.parse_args()))
