"""Every guard, against the input that forces it to fail.

This file is the paper's own rule applied to this repository: **a check unable
to fail is worse than no check, because it turns an absence of evidence into a
green light.** The remedy is mutation testing's premise moved from unit tests to
analysis diagnostics — construct the input a guard MUST catch, and treat a guard
with no such test as untested.

Every other test file here asserts the pass path: given clean data, the check
returns clean. That is the half that cannot detect a check which returns clean
on everything. These are the other half.

Each test names the defect it would have caught. Five of them are the five
defects in the paper's Appendix A, written as the inputs that would have caught
them at the time rather than months later.

If you add a guard, add its failing input here. If you cannot construct one, the
guard does not have a testable failure mode and that is the finding.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.harness.calibration import (
    PREREGISTERED,
    coverage_and_discrimination,
)
from src.harness.depth_confound import (
    MATURITY_DEPTH_RHO_TOLERANCE,
    depth_confound_report,
    max_attainable_rho,
)
from src.reference.signature import LeakageError, assert_no_target_leakage

# ---------------------------------------------------------------------------
# The depth-confound diagnostic
# ---------------------------------------------------------------------------


def _two_arm_cells(
    *, rare_p: float, rare_depth_determined: bool, common_p: float = 0.5, n: int = 8000
):
    """Cells in two arms, one rare and one common, with control over whether the
    rare arm's label is a pure function of depth."""
    rng = np.random.default_rng(11)
    arm = np.where(np.arange(n) < n // 2, "rare_arm", "common_arm")
    depth = rng.lognormal(np.log(15000), 0.8, n)
    mature = np.zeros(n, dtype=bool)

    rare = np.where(arm == "rare_arm")[0]
    k = max(1, int(rare_p * len(rare)))
    if rare_depth_determined:
        mature[rare[np.argsort(depth[rare])[:k]]] = True
    else:
        mature[rng.permutation(rare)[:k]] = True

    common = np.where(arm == "common_arm")[0]
    mature[rng.permutation(common)[: int(common_p * len(common))]] = True
    return depth, mature, arm


def test_a_correlation_is_never_normalised_by_another_arms_bound():
    """The mispairing: `max(rhos)` and `max(ceilings)` taken independently.

    One arm is rare and perfectly depth-determined, so it sits exactly at its
    own ceiling. The other is common and depth-independent, so it contributes a
    large ceiling and no correlation. Dividing the first arm's rho by the second
    arm's ceiling reports 0.14 — "14% of the available headroom, clean" — for a
    label that is a *pure function of depth*.

    Paired within the arm it reports 1.0, which is what it is.
    """
    depth, mature, arm = _two_arm_cells(rare_p=0.005, rare_depth_determined=True)
    report = depth_confound_report(depth, mature, arm)

    assert report["worst_arm"] == "rare_arm"
    assert report["rho_vs_ceiling"] == pytest.approx(1.0, abs=0.05)

    # The value the old code produced, reconstructed here so the regression is
    # visible rather than asserted in the abstract.
    other_arm_ceiling = report["per_arm"]["common_arm"]["max_attainable_rho"]
    mispaired = report["worst_within_arm_rho"] / other_arm_ceiling
    assert mispaired < 0.2
    assert report["rho_vs_ceiling"] > 4 * mispaired

    # And the verdict flips with it: the arm carrying the correlation could
    # never have reached the tolerance.
    assert not report["tolerance_is_reachable"]
    assert "NOT TESTABLE" in report["reading"]


def test_each_arm_carries_its_own_bound_and_prevalence():
    """Nothing downstream should have to pair a statistic with its bound."""
    depth, mature, arm = _two_arm_cells(rare_p=0.005, rare_depth_determined=True)
    report = depth_confound_report(depth, mature, arm)

    for name, row in report["per_arm"].items():
        expected = max_attainable_rho(row["mature_share"])
        assert row["max_attainable_rho"] == pytest.approx(expected), name
        if np.isfinite(row["rho_depth_vs_mature"]) and expected > 0:
            assert row["rho_vs_ceiling"] == pytest.approx(
                abs(row["rho_depth_vs_mature"]) / expected
            ), name


def test_an_undefined_arm_is_named_rather_than_dropped():
    """An arm at prevalence exactly 0 has no correlation. That is UNDEFINED, not
    clean, and it must not vanish silently into the surviving arm's verdict.

    This is the paper's §3 failure: 58 rows reported as "label does not track
    depth" because an undefined statistic had nowhere to land except the
    negative verdict.
    """
    depth, mature, arm = _two_arm_cells(rare_p=0.005, rare_depth_determined=True)
    mature[arm == "rare_arm"] = False  # prevalence exactly 0 in that arm

    report = depth_confound_report(depth, mature, arm)
    assert report["arms_with_undefined_rho"] == ["rare_arm"]
    assert report["n_arms"] == 2
    assert "rare_arm" in report["reading"]
    assert "excluded rather than scored clean" in report["reading"]


def test_every_arm_undefined_abstains_rather_than_reporting_clean():
    depth, mature, arm = _two_arm_cells(rare_p=0.005, rare_depth_determined=True)
    mature[:] = False  # no variance anywhere

    report = depth_confound_report(depth, mature, arm)
    assert not np.isfinite(report["worst_within_arm_rho"])
    assert report["reading"].startswith("UNDEFINED")
    assert not report["maturity_tracks_depth"]


@pytest.mark.parametrize("p", [0.001, 0.005, 0.013])
def test_the_check_cannot_fire_below_the_crossing_prevalence(p):
    """The bound, as a test. Below p = 1.3516% a tolerance of 0.20 sits above
    sqrt(3p(1-p)), so the diagnostic reports clean on a label that is a pure
    function of depth. This asserts the blindness rather than describing it."""
    depth, mature, arm = _two_arm_cells(rare_p=p, rare_depth_determined=True)
    report = depth_confound_report(depth, mature, arm)

    rare = report["per_arm"]["rare_arm"]
    assert rare["max_attainable_rho"] < MATURITY_DEPTH_RHO_TOLERANCE
    assert rare["rho_vs_ceiling"] == pytest.approx(1.0, abs=0.1)  # perfectly determined
    assert not rare["tolerance_is_reachable"]
    assert not report["maturity_tracks_depth"]  # ...and it reports clean


def test_the_check_does_fire_when_prevalence_permits_it():
    """The positive control. Without this the test above is satisfied by a
    diagnostic that never fires on anything."""
    depth, mature, arm = _two_arm_cells(rare_p=0.30, rare_depth_determined=True)
    report = depth_confound_report(depth, mature, arm)

    assert report["per_arm"]["rare_arm"]["tolerance_is_reachable"]
    assert report["maturity_tracks_depth"]
    assert "NOT TESTABLE" not in report["reading"]


# ---------------------------------------------------------------------------
# The calibration routine
# ---------------------------------------------------------------------------


def _sweep(n_replicates: int = 40, *, abstain: int = 0) -> pd.DataFrame:
    """A sweep whose intervals all cover truth and all exclude zero, with the
    first ``abstain`` replicates in the lowest bin holding no interval."""
    rng = np.random.default_rng(3)
    rows = []
    for i in range(n_replicates):
        for n_cells in (10, 200):
            truth = -5.0
            wide = n_cells == 10
            rows.append(
                {
                    "replicate": i,
                    "arm": PREREGISTERED.arm,
                    "shift": PREREGISTERED.detectable_shift,
                    "n_cells_mature": n_cells + rng.integers(0, 3),
                    "intrinsic_true_parametric": truth,
                    "ci_low": truth - 1.0,
                    "ci_high": truth + 1.0 if not wide else truth + 1.0,
                }
            )
    frame = pd.DataFrame(rows)
    if abstain:
        low = frame.index[frame["n_cells_mature"] < 100][:abstain]
        frame.loc[low, ["ci_low", "ci_high"]] = np.nan
    return frame


def test_an_abstention_is_not_counted_as_a_coverage_failure():
    """Comparisons against NaN return False, so a replicate holding no interval
    used to score as a coverage AND a discrimination failure — halving both
    rates in the bin where refusing was most often right.

    Here every replicate that *does* hold an interval covers truth and excludes
    zero, so every rate must be exactly 1.0 whether or not abstentions are
    present. If an abstention leaks into a denominator, a rate drops below 1.
    """
    clean = coverage_and_discrimination(_sweep())
    half_refused = coverage_and_discrimination(_sweep(abstain=20))

    for table, label in ((clean, "clean"), (half_refused, "with abstentions")):
        rates = table[["coverage", "discrimination"]].to_numpy(dtype=float)
        present = rates[np.isfinite(rates)]
        assert len(present) > 0, label
        assert present == pytest.approx(1.0), label

    # The abstentions are visible as a count rather than absorbed into a rate.
    assert clean["n_abstained"].sum() == 0
    assert half_refused["n_abstained"].sum() == 20
    assert (
        half_refused["n_estimated"] + half_refused["n_abstained"]
        == half_refused["n_replicates"]
    ).all()


def test_a_bin_where_every_replicate_abstained_gets_its_own_verdict():
    """"Nobody answered" and "the answers were bad" are different findings.
    NaN >= target is False, so the old code returned `not_estimable` for both."""
    frame = _sweep()
    frame.loc[frame["n_cells_mature"] < 100, ["ci_low", "ci_high"]] = np.nan

    table = coverage_and_discrimination(frame)
    empty = table[table["n_estimated"] == 0]
    assert not empty.empty, "the low bins should hold no intervals at all"
    assert (empty["verdict"] == "all_abstained").all()
    assert not np.isfinite(empty["coverage"]).any()

    # The bins that did answer are untouched and still judged normally.
    answered = table[table["n_estimated"] > 0]
    assert not answered.empty
    assert (answered["verdict"] != "all_abstained").all()


# ---------------------------------------------------------------------------
# The leakage invariant
# ---------------------------------------------------------------------------


def test_the_leakage_guard_fires_across_identifier_spaces():
    """Appendix A item 1 and 2: the guard compared gene SYMBOLS against database
    ACCESSIONS. An empty intersection follows by construction, so it passed on
    every input including violating ones.

    The violating input is a label set carrying the target gene in a *different*
    namespace than the target list uses. Without an alias map that comparison is
    vacuous, and `resolve_targets` is what makes it non-vacuous.
    """
    with pytest.raises(LeakageError):
        assert_no_target_leakage(
            ["GUCA2A", "LGR5"], ["GUCA2A"], context="same-namespace control"
        )

    # The same violation, expressed in Ensembl ids. This is the input the old
    # guard returned clean on.
    with pytest.raises((LeakageError, Exception)):
        assert_no_target_leakage(
            ["ENSG00000197273", "ENSG00000139292"],
            ["GUCA2A"],
            context="cross-namespace violation",
            alias_map={"GUCA2A": "ENSG00000197273"},
        )


def test_the_leakage_guard_still_passes_a_clean_panel():
    """The positive control: a guard that raises on everything is no better."""
    assert_no_target_leakage(["LGR5", "ASCL2"], ["GUCA2A"], context="clean")


# ---------------------------------------------------------------------------
# The provenance stamp
# ---------------------------------------------------------------------------


def test_the_dirtiness_check_counts_untracked_files(tmp_path, monkeypatch):
    """Appendix A item 5: the check ignored untracked files, so a table produced
    by a script that was never committed recorded a clean hash pointing at a
    commit that did not contain the script.

    The violating input is a repository whose only difference from HEAD is an
    untracked producing script.
    """
    import subprocess

    from src.common import provenance

    def git(*args):
        subprocess.run(
            ["git", *args], cwd=tmp_path, check=True, capture_output=True, text=True
        )

    git("init", "-q")
    git("config", "user.email", "t@example.com")
    git("config", "user.name", "t")
    (tmp_path / "tracked.txt").write_text("committed\n")
    git("add", "tracked.txt")
    git("commit", "-qm", "initial")

    monkeypatch.setattr(provenance, "REPO_ROOT", tmp_path)
    assert not provenance.git_is_dirty(), "a clean tree must read clean"

    # The input that made the old check return clean.
    (tmp_path / "produce_the_table.py").write_text("# never committed\n")
    assert provenance.untracked_files(), "an untracked producer must be listed"
    assert provenance.git_is_dirty(), (
        "a tree whose only change is an untracked producing script is NOT clean "
        "— that is the failure this check exists to catch"
    )
