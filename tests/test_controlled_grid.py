from __future__ import annotations

import numpy as np
import pandas as pd

from src.harness.controlled_grid import compare_grid_subsets, crossing_brackets


def _sweep() -> pd.DataFrame:
    rows = []
    for count in (10, 20, 40, 80):
        for replicate in range(10):
            valid = not (count == 10 and replicate < 5)
            half = 2.0 / np.sqrt(count)
            rows.append(
                {
                    "arm": "oracle",
                    "shift": 0.5,
                    "replicate": replicate,
                    "frac_mature_tumour": count / 2000,
                    "n_cells_mature": count,
                    "intrinsic_true_parametric": -1.0,
                    "ci_low": -1.0 - half if valid else np.nan,
                    "ci_high": -1.0 + half if valid else np.nan,
                }
            )
    return pd.DataFrame(rows)


def test_comparison_reports_all_binning_readings_and_counts_abstentions():
    rates = compare_grid_subsets(
        _sweep(),
        subsets={"small": (0.005, 0.01, 0.02, 0.04)},
        fixed_bin_edges=(0, 20, 50, 100),
    )
    assert set(rates["binning"]) == {"adaptive_bins", "fixed_bins", "per_count"}
    per_count = rates[rates["binning"] == "per_count"]
    assert per_count["n_replicates"].sum() == 40
    assert (per_count["n_attempted"] == per_count["n_replicates"]).all()
    assert per_count["n_abstained"].sum() == 5
    assert {"coverage_mc_se", "discrimination_mc_se"} <= set(rates)


def test_crossings_are_brackets_or_explicitly_missing_never_endorsements():
    rates = compare_grid_subsets(
        _sweep(),
        subsets={"small": (0.005, 0.01, 0.02, 0.04)},
        fixed_bin_edges=(0, 20, 50, 100),
    ).assign(cohort="synthetic", pool="reference", seed=1)
    brackets = crossing_brackets(rates)
    assert set(brackets["status"]) <= {
        "bracketed",
        "lower_bound_unobserved",
        "no_crossing",
        "coincident_cutpoints",
    }
    assert not brackets["status"].str.contains("pass|valid|endors", case=False).any()
