"""Combining per-study coexpression readings.

The dry run is the important test: against the committed three-cohort table the
combiner must return the answer this project already knows -- UNRESOLVED at
k = 3 -- because each cohort came back UNRESOLVED alone. A combiner that cannot
reproduce a known answer has not been checked.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.common.paths import RESULTS_DIR
from src.harness.meta import MIN_STUDIES
from src.reference.jobs.coexpression_meta import (
    CONTROLS,
    meta_control,
    meta_detection,
    meta_premise,
    per_study_stats,
)
from src.reference.jobs.coexpression_silencing import CONTROL_LOG2_TOLERANCE

COMMITTED = RESULTS_DIR / "2026-09-04_975cf5c" / "coexpression_silencing.parquet"


def _study(study: str, gene: str, *, log2: float, n: int, spread: float = 0.3,
           seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame({
        "study_id": study, "gene": gene,
        "patient_id": [f"{study}-p{i}" for i in range(n)],
        "log2_cp10k_ratio": rng.normal(log2, spread, n),
        "delta_detect": rng.normal(-0.4, 0.05, n),
    })


# ---------------------------------------------------------------------------
# The dry run: reproduce what we already know


def test_the_three_committed_cohorts_return_undecided():
    """The baseline the fourteen studies exist to beat."""
    if not COMMITTED.exists():
        pytest.skip("the three-cohort table is not committed")
    deltas = pd.read_parquet(COMMITTED)
    verdict, detail, controls, per_study = meta_premise(deltas)

    assert verdict == "UNRESOLVED", detail
    assert set(controls["gene"]) == set(CONTROLS)
    assert set(controls["verdict"]) == {"UNRESOLVED"}
    assert set(per_study["study_id"]) == {"GSE132465", "GSE144735", "GSE178341"}


def test_the_committed_cohorts_are_homogeneous_so_the_failure_is_precision():
    """I^2 near zero means the studies AGREE and simply lack precision.

    That is the whole argument for going to fourteen: if they disagreed, more
    of them would not help and the I^2 gate would refuse the pooled reading
    however many there were.
    """
    if not COMMITTED.exists():
        pytest.skip("the three-cohort table is not committed")
    deltas = pd.read_parquet(COMMITTED)
    for gene in CONTROLS:
        row, _ = meta_control(deltas, gene)
        assert row["i_squared"] < 0.5, (
            f"{gene} heterogeneity is {row['i_squared']:.1%}; the three cohorts "
            f"no longer agree, and more studies would not be the fix"
        )
        assert row["homogeneous"]


def test_actb_sits_just_inside_the_tolerance_and_that_is_the_hard_case():
    """Recorded before the fourteen land, because predicting it after is worthless.

    ACTB pools to about +0.487 against a 0.5 tolerance -- 2.6% below the line.
    More studies narrow an interval; they do not move a point estimate. So ACTB
    can only resolve to HOLDS by a hair, while KRT8 at about -0.37 has real
    room. If the fourteen come back UNRESOLVED, this is the reason to check
    first.
    """
    if not COMMITTED.exists():
        pytest.skip("the three-cohort table is not committed")
    deltas = pd.read_parquet(COMMITTED)
    actb, _ = meta_control(deltas, "ACTB")
    assert abs(actb["pooled"]) < CONTROL_LOG2_TOLERANCE
    assert CONTROL_LOG2_TOLERANCE - abs(actb["pooled"]) < 0.05, (
        "ACTB is no longer on the tolerance boundary; the note above is stale"
    )


def test_the_standard_error_agrees_with_premise_holds_bootstrap():
    """The meta weights by sem; the per-study check reports a bootstrap interval.

    Two answers to one question if they disagree, so they are pinned to agree.
    """
    if not COMMITTED.exists():
        pytest.skip("the three-cohort table is not committed")
    from src.reference.jobs.icbi_coexpression import control_log2_interval

    deltas = pd.read_parquet(COMMITTED)
    block = deltas[deltas["study_id"] == "GSE178341"]
    stats = per_study_stats(block, "ACTB").iloc[0]
    _, lo, hi = control_log2_interval(block, "ACTB", seed=7)
    boot_se = (hi - lo) / (2 * 1.96)
    assert stats["se"] == pytest.approx(boot_se, rel=0.25), (
        f"sem {stats['se']:.4f} against bootstrap {boot_se:.4f} -- the meta "
        f"would weight by a different quantity than the per-study check reports"
    )


# ---------------------------------------------------------------------------
# Every control must hold


def test_one_refused_control_refuses_the_premise():
    """ACTB fine, KRT8 far beyond tolerance. Pooling 'the worst per study'
    would have hidden which one."""
    frames = []
    for i in range(4):
        frames.append(_study(f"S{i}", "ACTB", log2=0.05, n=8, seed=i))
        frames.append(_study(f"S{i}", "KRT8", log2=2.0, n=8, spread=0.1, seed=i + 10))
    verdict, detail, controls, _ = meta_premise(pd.concat(frames, ignore_index=True))
    assert verdict == "REFUSED"
    assert "KRT8" in detail
    assert controls.set_index("gene").loc["ACTB", "verdict"] == "HOLDS"


def test_both_controls_inside_tolerance_holds():
    # Seed from the gene's POSITION, never hash(gene): Python randomises string
    # hashing per process, so that fixture generated different data on every
    # run and the test passed or failed by luck. Invariant 10 is about results,
    # and a test whose input is not reproducible is the same defect one level
    # down.
    frames = []
    for i in range(5):
        for j, gene in enumerate(CONTROLS):
            frames.append(_study(f"S{i}", gene, log2=0.02, n=10, spread=0.08,
                                 seed=100 * j + i))
    verdict, _, controls, _ = meta_premise(pd.concat(frames, ignore_index=True))
    assert verdict == "HOLDS"
    assert set(controls["verdict"]) == {"HOLDS"}


def test_too_few_studies_is_undecided_rather_than_an_error():
    frames = [_study(f"S{i}", gene, log2=0.02, n=8, seed=i)
              for i in range(MIN_STUDIES - 1) for gene in CONTROLS]
    verdict, detail, controls, _ = meta_premise(pd.concat(frames, ignore_index=True))
    assert verdict == "UNRESOLVED"
    assert str(MIN_STUDIES) in " ".join(controls["detail"])


def test_a_study_below_the_patient_floor_is_excluded_from_the_pool():
    """One patient's mean IS that patient; premise_holds refuses it and so does
    this."""
    frames = [_study(f"S{i}", "ACTB", log2=0.02, n=8, seed=i) for i in range(3)]
    frames.append(_study("TINY", "ACTB", log2=5.0, n=2, seed=99))
    row, per_study = meta_control(pd.concat(frames, ignore_index=True), "ACTB")
    assert bool(per_study.set_index("study_id").loc["TINY", "below_patient_floor"])
    assert row["k_studies"] == 3
    assert abs(row["pooled"]) < 0.5, "the two-patient study leaked into the pool"


# ---------------------------------------------------------------------------
# The gate


def test_detection_is_only_read_when_the_premise_holds():
    """The order is the whole point: a marker falling inside a population that
    has itself changed is not silencing."""
    import inspect

    from src.reference.jobs import coexpression_meta

    source = inspect.getsource(coexpression_meta.main)
    assert 'if verdict == "HOLDS":' in source
    assert source.index('verdict, detail, controls') < source.index("meta_detection")


def test_meta_detection_pools_the_target_genes():
    frames = [_study(f"S{i}", "GUCA2A", log2=0.0, n=8, seed=i) for i in range(5)]
    got = meta_detection(pd.concat(frames, ignore_index=True))
    row = got[got["gene"] == "GUCA2A"]
    assert len(row) == 1
    assert row.iloc[0]["k_studies"] == 5
    assert row.iloc[0]["pooled"] == pytest.approx(-0.4, abs=0.05)


def test_the_fixtures_are_reproducible():
    """A test whose input changes between runs passes or fails by luck.

    `hash(gene)` was seeding a fixture here, and Python randomises string
    hashing per process unless PYTHONHASHSEED is pinned -- so the same test
    exercised different data every time it ran.
    """
    import inspect

    from tests import test_coexpression_meta as module

    # Only SEED ARGUMENTS, and only outside comments. A scan for "hash(" over
    # the whole module matched this function's own docstring and then its own
    # implementation -- a check firing on itself, which is a fair warning that
    # source inspection was the wrong instrument for the job.
    offenders = [
        line.strip()
        for line in inspect.getsource(module).splitlines()
        if "seed=" in line.split("#", 1)[0] and "hash(" in line.split("#", 1)[0]
    ]
    assert not offenders, (
        f"a fixture is seeded from hash(): {offenders}. String hashing is "
        f"randomised per process, so the fixture is not reproducible."
    )
    first = _study("S0", "ACTB", log2=0.2, n=8, seed=4)
    second = _study("S0", "ACTB", log2=0.2, n=8, seed=4)
    pd.testing.assert_frame_equal(first, second)
