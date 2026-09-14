# Reproducing the WMHS audit examples

These commands exercise the new examples without any biological data or GPU.
They require a clean checkout for provenance-stamped output. For exploratory
runs from a modified tree, add `--allow-dirty`; the resulting sidecar will
record that state and should not be treated as submission evidence.

## Fresh environment

From the repository root, using Python 3.11:

```sh
python -m venv .venv-wmhs
.venv-wmhs/bin/python -m pip install --upgrade pip
.venv-wmhs/bin/python -m pip install -e '.[dev]'
```

The `dev` extra intentionally requires `lifelines`. Its comparison test is not
silently skipped when the package is absent.

## Independent Cox comparison

Run the fixed four-horizon comparison between the repository's Cox fits and
`lifelines.CoxPHFitter`:

```sh
.venv-wmhs/bin/python -m src.harness.lifelines_check \
  --results-dir /tmp/wmhs-reproduction
```

The command exits nonzero if either implementation returns an invalid value or
if any paired coefficient differs by more than the prospectively fixed
tolerance of `1e-5`. This one synthetic comparison checks implementation
agreement; it does not estimate failure prevalence across clinical pipelines.

## Residual/performance clean control

Generate the prospectively specified null and non-null comparison:

```sh
.venv-wmhs/bin/python -m src.harness.clean_control \
  --seed 20260913 --replicates 2000 --patients 200 \
  --results-dir /tmp/wmhs-reproduction
```

The output reports attempted and valid draws, numerical equality or departure
from the predeclared empirical reference, bias, RMSE, interval coverage and its
Monte Carlo standard error. The equality class contains both calibrated and
deliberately narrow intervals; the departure class contains both a valid
estimator and a known biased estimator.

## Calibration uncertainty

The fixed-pair sensitivity run reuses identical outer samples at three inner
bootstrap budgets and then performs patient omission as a separate analysis:

```sh
.venv-wmhs/bin/python -m src.harness.calibration_sensitivity \
  --cohort both --budgets 200 1000 5000 \
  --results-dir /tmp/wmhs-reproduction
```

The SMC run enumerates all 45 two-patient holdout pairs and KUL3 enumerates all
15. Its output includes per-count coverage and discrimination, crossing
brackets, and the number and identities of source patients whose omission
changes each candidate. These are cohort-conditional sensitivity results, not
population confidence intervals.

## Controlled grid and binning comparison

The final design runs the union of all count settings once for each cohort,
draw pool and seed, then subsets those same draws into committed, extended and
dense grids under adaptive, fixed and per-count summaries:

```sh
.venv-wmhs/bin/python -m src.harness.controlled_grid \
  --cohort both --replicates 200 --n-boot 200 --seeds 8 \
  --results-dir /tmp/wmhs-reproduction
```

This is a long CPU run (the 2026-09-13 pilot projected roughly 200 CPU-hours).
Use scheduled compute and a clean committed SHA. For a smoke test only, set
`--replicates 1 --n-boot 10 --seeds 1` and label the output as a pilot.

## Focused verification

```sh
.venv-wmhs/bin/python -m pytest \
  tests/test_lifelines_check.py \
  tests/test_clean_control.py \
  tests/test_trial_survival.py -q
```

For all implementation and paper gates, run:

```sh
.venv-wmhs/bin/python -m ruff check .
.venv-wmhs/bin/python -m pytest -q
cd paper/wmhs && ./build.sh && ./check_anonymity.sh
```
