"""Stage 4's instrument gate, against the inputs that must defeat it.

The gate's force comes entirely from purity being measured independently of the
expression the fractions were deconvolved from. Most of these tests are about
the ways that independence can be lost while the gate keeps returning a number.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.bulk.instrument import (
    INDEPENDENT_PURITY_METHOD,
    InstrumentError,
    gate_verdict,
    independent_purity,
    run_instrument_check,
)
from src.common.paths import RESULTS_DIR


def _purity(n: int = 40, *, method: str = INDEPENDENT_PURITY_METHOD,
            expression_derived: bool = False, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "barcode": [f"TCGA-XX-{i:04d}-01A-11R-AAAA-07" for i in range(n)],
        "method": method,
        "purity": rng.uniform(0.3, 0.9, n),
        "expression_derived": expression_derived,
    })


def _fractions(purity: pd.DataFrame, *, r_target: float, method: str = "nnls",
               seed: int = 2) -> pd.DataFrame:
    """Fractions whose non-epithelial column tracks (1 - purity) at ~r_target."""
    rng = np.random.default_rng(seed)
    signal = 1.0 - purity["purity"].to_numpy()
    noise = rng.normal(0, 1, len(signal))
    z = lambda v: (v - v.mean()) / v.std()  # noqa: E731
    mixed = r_target * z(signal) + np.sqrt(max(0.0, 1 - r_target**2)) * z(noise)
    return pd.DataFrame({
        "sample_id": purity["barcode"],
        "method": method,
        "non_epithelial_fraction": 0.4 + 0.1 * z(mixed),
        "mature_colonocyte_fraction": rng.uniform(0.1, 0.5, len(signal)),
    })


# ---------------------------------------------------------------------------
# It measures what it says it measures


def test_a_tracking_fraction_passes_and_a_random_one_does_not():
    purity = _purity()
    good = run_instrument_check(_fractions(purity, r_target=0.9), purity,
                                method="nnls", rung="lineage")
    assert good.passed and good.r >= 0.5

    bad = run_instrument_check(_fractions(purity, r_target=0.05, seed=9), purity,
                               method="nnls", rung="lineage")
    assert not bad.passed
    assert "STOP" in bad.detail


def test_the_threshold_is_the_prespec_s_and_binds_at_it():
    """r just under 0.5 fails, just over passes. A gate whose threshold is not
    load-bearing is a report, not a gate."""
    purity = _purity()
    over = run_instrument_check(_fractions(purity, r_target=0.75), purity,
                                method="nnls", rung="lineage")
    under = run_instrument_check(_fractions(purity, r_target=0.2, seed=13), purity,
                                 method="nnls", rung="lineage")
    assert over.threshold == under.threshold == 0.5
    assert over.passed and not under.passed


# ---------------------------------------------------------------------------
# The independence the gate's logic rests on


def test_expression_derived_purity_is_refused_even_when_the_method_name_matches():
    """The trap: 675 expression-derived calls against 556 independent ones.

    Correlating a fraction deconvolved FROM expression against a purity score
    derived FROM expression passes because both read the same signal. The
    filter is on the property, not only on the name, so a method that is
    expression-derived is excluded whatever it is called.
    """
    purity = _purity(method=INDEPENDENT_PURITY_METHOD, expression_derived=True)
    with pytest.raises(InstrumentError, match="expression-derived"):
        independent_purity(purity)


def test_a_differently_named_independent_method_is_not_silently_accepted():
    with pytest.raises(InstrumentError, match="no 'absolute' purity calls"):
        independent_purity(_purity(method="aran_cpe"))


def test_a_purity_table_without_the_independence_column_is_refused():
    purity = _purity().drop(columns=["expression_derived"])
    with pytest.raises(InstrumentError, match="independent of expression"):
        independent_purity(purity)


def test_independent_purity_keeps_only_the_absolute_rows():
    mixed = pd.concat([
        _purity(20, method="absolute", expression_derived=False),
        _purity(30, method="estimate_affy_extrapolated", expression_derived=True, seed=5),
        _purity(25, method="aran_cpe", expression_derived=False, seed=6),
    ])
    kept = independent_purity(mixed)
    assert set(kept["method"]) == {"absolute"}
    assert len(kept) == 20


# ---------------------------------------------------------------------------
# Undefined is not a pass


def test_a_constant_fraction_is_a_failure_not_a_pass():
    """`nan >= 0.5` is False, so a constant column must not reach the comparison.

    This is the coercion the paper documents: an undefined statistic compared
    against a threshold returns False in numpy either way, and a reader cannot
    tell "did not correlate" from "could not be evaluated".
    """
    purity = _purity()
    fractions = _fractions(purity, r_target=0.9)
    fractions["non_epithelial_fraction"] = 0.4
    result = run_instrument_check(fractions, purity, method="nnls", rung="lineage")
    assert not result.passed
    assert np.isnan(result.r)
    assert "constant" in result.detail and "undefined" in result.detail


def test_too_few_matched_samples_is_unevaluable_not_failed():
    purity = _purity(2)
    with pytest.raises(InstrumentError, match="unevaluable"):
        run_instrument_check(_fractions(purity, r_target=0.9), purity,
                             method="nnls", rung="lineage")


def test_no_overlap_between_fractions_and_purity_is_refused():
    purity = _purity()
    fractions = _fractions(purity, r_target=0.9)
    fractions["sample_id"] = [f"NOPE-{i}" for i in range(len(fractions))]
    with pytest.raises(InstrumentError, match="cannot be evaluated"):
        run_instrument_check(fractions, purity, method="nnls", rung="lineage")


def test_duplicate_aliquots_do_not_double_weight_a_patient():
    """Two RNA aliquots of one patient are one patient's worth of evidence."""
    purity = _purity(30)
    doubled = pd.concat([purity, purity.assign(
        barcode=purity["barcode"].str.replace("-01A-11R-AAAA-07", "-01A-22R-BBBB-07",
                                              regex=False))])
    single = run_instrument_check(_fractions(purity, r_target=0.8), purity,
                                  method="nnls", rung="lineage")
    with_dupes = run_instrument_check(_fractions(purity, r_target=0.8), doubled,
                                      method="nnls", rung="lineage")
    assert with_dupes.n_samples == single.n_samples


# ---------------------------------------------------------------------------
# The verdict across methods


def test_one_passing_method_lets_stage4_proceed_and_none_does_not():
    purity = _purity()
    good = run_instrument_check(
        _fractions(purity, r_target=0.9, method="nusvr"), purity,
        method="nusvr", rung="lineage",
    )
    bad = run_instrument_check(_fractions(purity, r_target=0.05, seed=9), purity,
                               method="nnls", rung="lineage")
    proceed, message = gate_verdict([good, bad])
    assert proceed and "1 of 2" in message

    stopped, message = gate_verdict([bad])
    assert not stopped and message.startswith("STOP")

    never_run, message = gate_verdict([])
    assert not never_run and "never run" in message


# ---------------------------------------------------------------------------
# Against the committed purity table


def test_the_committed_table_has_absolute_calls_that_are_not_expression_derived():
    """556 ABSOLUTE against 675 expression-derived, on the real artifact."""
    matches = sorted(RESULTS_DIR.glob("*/tcga_purity.parquet"))
    if not matches:
        pytest.skip("no committed purity table")
    purity = pd.read_parquet(matches[-1])
    kept = independent_purity(purity)
    assert len(kept) > 100
    assert not kept["expression_derived"].any()
    assert set(kept["method"]) == {INDEPENDENT_PURITY_METHOD}
    # The one that must not be used is the larger of the two.
    derived = purity[purity["expression_derived"].astype(bool)]
    assert len(derived) > len(kept), (
        "the expression-derived calls are no longer the larger set, which was "
        "the reason this filter is easy to get wrong -- re-check the docstring"
    )
