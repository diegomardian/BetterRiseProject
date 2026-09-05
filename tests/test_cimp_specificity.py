"""The CIMP screening, against the inputs that must break it.

Pre-registered in docs/prereg_cimp_specificity.md. The decision rule is the part
worth testing hardest: it exists to remove a degree of freedom, and a rule that
can be satisfied by whichever reference happens to agree removes nothing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.bulk.run_cimp_specificity import (
    REFERENCES,
    TARGET,
    paired_differences,
    stratum_contrast,
    verdict,
)


def _rows(**contrasts) -> pd.DataFrame:
    """One unadjusted row per reference, with the contrast and whether it clears."""
    return pd.DataFrame([
        {"reference": ref, "adjustment": "none", "contrast": c,
         "ci_low": c - 0.05, "ci_high": c + 0.05, "excludes_zero": abs(c) > 0.05}
        for ref, c in contrasts.items()
    ])


# ---------------------------------------------------------------------------
# The decision rule
# ---------------------------------------------------------------------------


def test_one_agreeable_reference_is_not_support():
    """The violating input is the configuration the committed medians show.

    GUCA2A falls less than CDX2 and more than MS4A12. A rule satisfied by
    either reference alone would license "specific" from exactly this, by
    quoting MS4A12 and not CDX2 -- choosing the reference after seeing the
    answer, which is what the two-reference design exists to prevent.
    """
    supported, reading = verdict(_rows(CDX2=+0.31, MS4A12=-0.12))
    assert not supported
    assert "NOT SPECIFIC" in reading
    assert "CDX2" in reading and "MS4A12" in reading, "both references get named"


def test_both_references_agreeing_is_support_and_escalates():
    supported, reading = verdict(_rows(CDX2=-0.40, MS4A12=-0.35))
    assert supported
    assert "SPECIFIC" in reading and "purity" in reading


def test_a_large_contrast_that_does_not_clear_zero_is_not_support():
    """Direction is not evidence. A wide interval straddling zero is undecided,
    and the rule must not read the point estimate through it."""
    rows = _rows(CDX2=-0.40, MS4A12=-0.35)
    rows.loc[:, ["ci_low", "ci_high", "excludes_zero"]] = [-1.2, 0.4, False]
    assert not verdict(rows)[0]


def test_a_missing_reference_is_undefined_not_support():
    supported, reading = verdict(_rows(CDX2=-0.40))
    assert not supported and "UNDEFINED" in reading


# ---------------------------------------------------------------------------
# Why the estimand is a within-sample difference
# ---------------------------------------------------------------------------


def test_pairing_inside_a_sample_cancels_a_per_sample_scale_factor():
    """The reason for the design, stated as a test.

    Give every sample its own loading factor, hitting both genes equally. The
    paired difference must not move; a difference of separately-taken stratum
    means does, because the factor does not cancel across strata unless it
    happens to balance.
    """
    rng = np.random.default_rng(0)
    n = 60
    clean = pd.DataFrame(
        {f"s{i}": {TARGET: 6.0, "CDX2": 7.0} for i in range(n)}
    )
    loading = pd.Series(rng.normal(0, 1.5, n), index=clean.columns)
    noisy = clean + loading

    assert paired_differences(clean, "CDX2").std() == pytest.approx(0.0, abs=1e-12)
    assert paired_differences(noisy, "CDX2").std() == pytest.approx(0.0, abs=1e-12), (
        "a factor hitting both genes equally must divide out of the pair"
    )
    # The unpaired quantity does carry it, which is what the pairing avoids.
    assert noisy.loc[TARGET].std() > 1.0


def test_the_strata_are_resampled_separately():
    """A pooled resample of 91+405 can return splits the design never had.
    Stratified draws keep each arm's size fixed, so the interval describes this
    design rather than a hypothetical one."""
    rng = np.random.default_rng(1)
    d = pd.Series(np.concatenate([rng.normal(-0.5, 0.2, 91), rng.normal(0.0, 0.2, 405)]))
    pos = np.zeros(len(d), bool)
    pos[:91] = True
    neg = ~pos
    out = stratum_contrast(d, pos, neg, seed=7)
    assert out["n_positive"] == 91 and out["n_negative"] == 405
    assert out["ci_low"] < out["contrast"] < out["ci_high"]
    assert out["excludes_zero"], "a real half-log2 separation must clear zero"


def test_a_stratum_too_small_to_resample_returns_no_interval():
    d = pd.Series([0.1, 0.2, 0.3, 0.4, 0.5])
    pos = np.array([True, True, False, False, False])
    out = stratum_contrast(d, pos, ~pos, seed=1)
    assert out["n_positive"] == 2
    assert not np.isfinite(out["ci_low"]) and not out["excludes_zero"]


def test_the_sign_convention_is_the_one_the_prereg_states():
    """Negative means the target fell FURTHER in CIMP+ than the reference --
    the direction locus-specific silencing predicts."""
    d = pd.Series([-1.0] * 20 + [0.0] * 20)      # target lower in the first arm
    pos = np.array([True] * 20 + [False] * 20)
    assert stratum_contrast(d, pos, ~pos, seed=1)["contrast"] < 0


def test_a_reference_absent_from_the_panel_raises():
    frame = pd.DataFrame({f"s{i}": {TARGET: 6.0} for i in range(5)})
    with pytest.raises(KeyError, match="CDX2"):
        paired_differences(frame, "CDX2")


def test_every_declared_reference_has_a_stated_role():
    """A reference without a reason is a reference chosen for its answer."""
    assert len(REFERENCES) == 2
    assert all(role.strip() for role in REFERENCES.values())


# ---------------------------------------------------------------------------
# run(), and the alignment it depends on
# ---------------------------------------------------------------------------


def _cohort(n_pos: int = 91, n_neg: int = 405, *, shuffle_metadata: bool = False):
    """Synthetic GSE39582 carrying the committed gaps.

    Expression is genes x samples; metadata is samples x fields. The two are
    joined by sample id and nothing checks that they are in the same order,
    which is the point of the shuffled variant.

    The arm sizes are the real ones, 91 and 405. At 40/80 the expected MS4A12
    contrast of -0.123 sits inside one standard error and the sign is not
    reliable -- which is a fair preview of why that arm's interval on real data
    may well contain zero.
    """
    rng = np.random.default_rng(4)
    samples = [f"s{i}" for i in range(n_pos + n_neg)]
    gaps = {"GUCA2A": 0.654, "CDX2": 0.968, "MS4A12": 0.531}
    base = {"GUCA2A": 6.523, "CDX2": 7.517, "MS4A12": 4.532}
    loading = rng.normal(0, 1.2, len(samples))
    values = {}
    for gene, gap in gaps.items():
        v = np.full(len(samples), base[gene]) + loading + rng.normal(0, 0.15, len(samples))
        v[:n_pos] -= gap
        values[gene] = v
    expression = pd.DataFrame(values, index=samples).T
    metadata = pd.DataFrame({
        "dataset": ["discovery"] * len(samples),
        "cimp.status": ["+"] * n_pos + ["-"] * n_neg,
        "mmr.status": rng.choice(["dMMR", "pMMR"], len(samples), p=[0.15, 0.85]),
        "tumor.location": rng.choice(["proximal", "distal"], len(samples)),
    }, index=samples)
    if shuffle_metadata:
        metadata = metadata.sample(frac=1.0, random_state=7)
    return expression, metadata


def test_run_produces_both_references_and_the_adjusted_reads():
    from src.bulk.run_cimp_specificity import run

    table = run(*_cohort())
    unadjusted = table[table["adjustment"] == "none"]
    assert set(unadjusted["reference"]) == set(REFERENCES)
    assert set(table["adjustment"]) == {"none", "mmr.status", "tumor.location"}
    # Direction reproduces the committed structure: CDX2 positive, MS4A12 negative.
    signs = unadjusted.set_index("reference")["contrast"]
    assert signs["CDX2"] > 0, "GUCA2A falls less than CDX2 in CIMP+"
    assert signs["MS4A12"] < 0, "GUCA2A falls more than MS4A12 in CIMP+"
    assert not verdict(table)[0]


def test_run_is_wrong_if_metadata_row_order_is_trusted_over_sample_ids():
    """The violating input: metadata whose rows are shuffled relative to
    expression.

    `strata()` builds boolean masks positionally from the metadata frame, and
    those masks index the expression columns. Nothing asserts the two share an
    order, so a shuffled metadata silently assigns CIMP status to the wrong
    samples -- and the run still completes, reporting a contrast near zero
    because the strata are now random. Exactly the class of silent misalignment
    this repository keeps finding.
    """
    from src.bulk.run_cimp_specificity import run

    aligned = run(*_cohort())
    shuffled = run(*_cohort(shuffle_metadata=True))

    a = aligned[aligned["adjustment"] == "none"].set_index("reference")["contrast"]
    b = shuffled[shuffled["adjustment"] == "none"].set_index("reference")["contrast"]
    assert abs(b["CDX2"]) < abs(a["CDX2"]) / 2, (
        "shuffling metadata must destroy the signal; if it does not, the strata "
        "were not being read from the metadata at all"
    )
