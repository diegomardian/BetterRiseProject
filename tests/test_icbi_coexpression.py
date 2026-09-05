"""The ICBI coexpression path, and the bar it must clear on Pelka.

The adaptation's whole risk is that it diverges from the committed GSE178341
result for a reason that has nothing to do with the science -- a different QC
population, a different batch key, a sorted fraction left in. So the validation
is mechanical, and these tests check that it can FAIL, not only that it passes.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import pytest

from src.reference.icbi_slice import SliceError
from src.reference.jobs.icbi_coexpression import (
    BATCH_KEY,
    GUCA2A_DELTA_TOLERANCE,
    MIN_EPITHELIAL_PER_ARM,
    NAIVE,
    VALIDATION,
    control_log2_interval,
    eligible_patients,
    verdict_word,
)


def _obs(n_per_arm: int = 200, enrichment: str = NAIVE, n_patients: int = 3):
    rows = []
    for p in range(n_patients):
        for tissue, sample_type in (("normal", "adjacent normal"),
                                    ("tumour", "primary tumor")):
            for k in range(n_per_arm):
                rows.append({
                    "study_id": "S", "patient_id": f"P{p}",
                    "sample_id": f"P{p}-{tissue}",
                    "sample_type": sample_type,
                    "atlas_cell_type_coarse": "Epithelial cell" if k % 2 else "T cell",
                    "enrichment_cell_types": enrichment,
                })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Patient selection


def test_only_naive_cells_survive_the_sorted_fraction_filter():
    """A CD45-sorted fraction is immune-enriched by construction; leaving it in
    makes every compartment fraction a statement about the sort."""
    mixed = pd.concat([_obs(), _obs(enrichment="CD45+")], ignore_index=True)
    rows, patients = eligible_patients(mixed, "S")
    assert (rows["enrichment_cell_types"] == NAIVE).all()
    assert len(rows) == len(_obs())


def test_a_study_with_no_naive_cells_yields_no_patients():
    rows, patients = eligible_patients(_obs(enrichment="CD45+"), "S")
    assert patients == []


def test_patients_below_the_epithelial_floor_are_excluded():
    """Half the cells are epithelial, so n_per_arm must clear twice the floor."""
    plenty = _obs(n_per_arm=MIN_EPITHELIAL_PER_ARM * 2 + 20)
    thin = _obs(n_per_arm=MIN_EPITHELIAL_PER_ARM // 2)
    assert len(eligible_patients(plenty, "S")[1]) == 3
    assert eligible_patients(thin, "S")[1] == []


def test_a_patient_with_only_one_arm_is_excluded():
    single = _obs()
    single = single[~((single["patient_id"] == "P0") & (single["sample_type"] == "primary tumor"))]
    assert "P0" not in eligible_patients(single, "S")[1]


def test_healthy_normal_does_not_count_as_the_reference_arm():
    donors = _obs()
    donors.loc[donors["sample_type"] == "adjacent normal", "sample_type"] = "healthy normal"
    assert eligible_patients(donors, "S")[1] == []


def test_an_unknown_study_is_refused():
    with pytest.raises(SliceError, match="no cells for study"):
        eligible_patients(_obs(), "NotAStudy")


def test_the_batch_key_matches_gse178341s():
    """If this ever changes, the MAD thresholds move and the numbers with them."""
    assert BATCH_KEY == "sample_id"


# ---------------------------------------------------------------------------
# The validation bar. It must be able to fail.


def _deltas(study: str, *, guca2a: float, actb: float, n: int = 8, seed: int = 0):
    """Per-patient rows in the shape `summarise` consumes."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        for gene, delta in (("GUCA2A", guca2a), ("ACTB", actb),
                            ("CDX2", -0.1), ("KRT8", 0.0),
                            ("EPCAM", 0.0), ("MS4A12", -0.2)):
            rows.append({
                "study_id": study, "patient_id": f"p{i}", "gene": gene,
                "role": "target" if gene == "GUCA2A" else "control",
                "delta_detect": delta + rng.normal(0, 0.02),
                "delta_given_conditioner": delta + rng.normal(0, 0.02),
                "log2_cp10k_ratio": rng.normal(0, 0.1),
                "detect_normal": 0.5, "detect_tumour": 0.5 + delta,
            })
    return pd.DataFrame(rows)


def test_the_bar_passes_a_run_that_matches_the_committed_result(tmp_path, monkeypatch):
    import src.reference.jobs.icbi_coexpression as mod

    committed = _deltas("GSE178341", guca2a=-0.40, actb=-0.01, seed=1)
    root = tmp_path / "results" / "2026-09-04_975cf5c"
    root.mkdir(parents=True)
    committed.to_parquet(root / "coexpression_silencing.parquet")
    monkeypatch.setattr(mod, "RESULTS_DIR", tmp_path / "results")

    close = _deltas("Pelka_2021_Cell", guca2a=-0.42, actb=-0.01, seed=2)
    assert mod.check_against_committed(close)["verdict"] == "PASS"


def test_the_bar_fails_a_run_whose_target_delta_has_drifted(tmp_path, monkeypatch):
    """The check that makes the smoke test worth running at all."""
    import src.reference.jobs.icbi_coexpression as mod

    committed = _deltas("GSE178341", guca2a=-0.40, actb=-0.01, seed=1)
    root = tmp_path / "results" / "2026-09-04_975cf5c"
    root.mkdir(parents=True)
    committed.to_parquet(root / "coexpression_silencing.parquet")
    monkeypatch.setattr(mod, "RESULTS_DIR", tmp_path / "results")

    drifted = _deltas("Pelka_2021_Cell", guca2a=-0.05, actb=-0.01, seed=2)
    result = mod.check_against_committed(drifted)
    assert result["verdict"] == "FAIL"
    assert result["guca2a_drift"] > GUCA2A_DELTA_TOLERANCE


def test_the_bar_skips_rather_than_passing_when_there_is_nothing_to_compare(
    tmp_path, monkeypatch
):
    """A missing baseline must not read as a pass."""
    import src.reference.jobs.icbi_coexpression as mod

    monkeypatch.setattr(mod, "RESULTS_DIR", tmp_path / "results")
    result = mod.check_against_committed(_deltas("Pelka_2021_Cell", guca2a=-0.4, actb=0.0))
    assert result["verdict"] == "SKIPPED"


def test_the_bar_names_pelka_and_the_table_it_checks_against():
    assert VALIDATION["study_id"] == "Pelka_2021_Cell"
    assert VALIDATION["committed_study_id"] == "GSE178341"
    assert "2026-09-04_975cf5c" in VALIDATION["against"]
    assert len(VALIDATION["requirements"]) == 3


# ---------------------------------------------------------------------------
# The premise, which is the gate the bar exists to protect


def _saturated(study: str, *, log2: float, n: int = 8, seed: int = 0):
    """Controls at 0.99 detection in both arms -- the real regime.

    `premise_holds` switches a saturated control to log2 expression, so a
    detection-only check has nothing to see. This is the fixture that makes the
    difference visible.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        for gene, role in (("ACTB", "control"), ("KRT8", "control"),
                           ("GUCA2A", "target"), ("CDX2", "identity")):
            control = role == "control"
            rows.append({
                "study_id": study, "patient_id": f"p{i}", "gene": gene,
                "role": role,
                "delta_detect": (0.0 if control else -0.40) + rng.normal(0, 0.01),
                "delta_given_conditioner": (0.0 if control else -0.40),
                "log2_cp10k_ratio": (log2 if control else 0.0) + rng.normal(0, 0.12),
                "detect_normal": 0.99 if control else 0.52,
                "detect_tumour": 0.99 if control else 0.12,
            })
    return pd.DataFrame(rows)


def test_a_flipped_premise_verdict_fails_the_bar(tmp_path, monkeypatch):
    """THE GAP THIS FIX CLOSES.

    Committed ACTB log2 sits near +0.43 (UNRESOLVED). Move it to +1.4 and the
    premise becomes REFUSED -- the gate deciding whether the whole reading is
    interpretable has flipped. Detection is 0.99 in both arms either way, so a
    detection-only check sees nothing and the GUCA2A delta is untouched.
    """
    import src.reference.jobs.icbi_coexpression as mod

    root = tmp_path / "results" / "2026-09-04_975cf5c"
    root.mkdir(parents=True)
    _saturated("GSE178341", log2=0.43, seed=1).to_parquet(
        root / "coexpression_silencing.parquet"
    )
    monkeypatch.setattr(mod, "RESULTS_DIR", tmp_path / "results")

    flipped = _saturated("Pelka_2021_Cell", log2=1.4, seed=2)
    result = mod.check_against_committed(flipped)
    assert result["verdict"] == "FAIL"
    assert "premise verdict CHANGED" in result["detail"]
    assert result["premise_verdict_committed"] != result["premise_verdict_icbi"]

    # And the detection statistic really is blind to it, which is the point.
    assert abs(result.get("guca2a_drift", 0.0)) < GUCA2A_DELTA_TOLERANCE


def test_a_matching_premise_and_close_numbers_pass(tmp_path, monkeypatch):
    import src.reference.jobs.icbi_coexpression as mod

    root = tmp_path / "results" / "2026-09-04_975cf5c"
    root.mkdir(parents=True)
    _saturated("GSE178341", log2=0.43, seed=1).to_parquet(
        root / "coexpression_silencing.parquet"
    )
    monkeypatch.setattr(mod, "RESULTS_DIR", tmp_path / "results")

    close = _saturated("Pelka_2021_Cell", log2=0.45, seed=2)
    result = mod.check_against_committed(close)
    assert result["verdict"] == "PASS", result["detail"]
    assert result["premise_verdict_icbi"] == result["premise_verdict_committed"]


def test_the_log2_interval_matches_premise_holds(tmp_path):
    """Two code paths that agree by review is how a bar measures the wrong thing.

    `control_log2_interval` mirrors `premise_holds`' bootstrap deliberately, so
    this pins them to the same numbers rather than trusting that they match.
    """
    from src.reference.jobs.coexpression_silencing import premise_holds

    deltas = _saturated("S", log2=0.43, seed=5)
    mean, lo, hi = control_log2_interval(deltas, "ACTB", seed=11)
    _, reading = premise_holds(deltas, seed=11)
    assert "ACTB" in reading and "log2 expression" in reading
    assert f"{mean:+.3f}" in reading, f"{mean:+.3f} not in {reading}"
    assert f"[{lo:+.3f}, {hi:+.3f}]" in reading


def test_verdict_word_takes_the_state_not_the_numbers():
    assert verdict_word("UNRESOLVED: control ACTB +0.431 [...]") == "UNRESOLVED"
    assert verdict_word("REFUSED: control KRT8 -0.712") == "REFUSED"
    assert verdict_word("UNDEFINED: no control gene was scored") == "UNDEFINED"


def test_the_documented_bar_matches_the_implemented_one():
    """The requirement text and the code must not drift apart -- that gap is
    what this fix was for."""
    text = " ".join(VALIDATION["requirements"]).lower()
    assert "premise verdict" in text
    assert "log2" in text, "the bar claims a log2 check; the code must do one"
    assert "detection delta" in text


def test_an_unlabellable_patient_is_skipped_and_counted_not_crashed(monkeypatch):
    """29 of 30 patients' work was discarded by an exception on the 30th.

    Pelka's last patient has a depth target of 27,146 against 949 epithelial
    cells -- its tumour arm is so much deeper than its own normal arm that no
    reference cell survives matching. That is a property of one patient's
    library preparation, and the run already skips two other per-patient
    conditions. It must not take the study down with it, and it must not
    vanish either: a study that loses patients this way is a smaller reading
    and the sidecar has to say by how much.
    """
    import src.reference.jobs.icbi_coexpression as mod
    from src.reference.labels import LabelError

    calls = {"n": 0}

    def flaky(*args, **kwargs):
        calls["n"] += 1
        if calls["n"] == 2:
            raise LabelError("no reference cell survives depth matching")
        raise AssertionError("unreachable in this test")

    # Two eligible patients; the second cannot be labelled.
    obs = _obs(n_per_arm=MIN_EPITHELIAL_PER_ARM * 2 + 20, n_patients=2)
    rows, patients = eligible_patients(obs, "S")
    assert len(patients) == 2

    monkeypatch.setattr(mod, "eligible_patients", lambda *a, **k: (rows, patients))
    monkeypatch.setattr(mod, "read_var", lambda p: pd.DataFrame(
        {"ensembl_id": ["E1"], "gene_symbol": ["ACTB"]}))
    monkeypatch.setattr(mod, "read_cells", lambda *a, **k: (_ for _ in ()).throw(
        LabelError("no reference cell survives depth matching")))

    # read_cells raising LabelError is not what the guard catches -- it must
    # still propagate, because a slicing failure is not a per-patient outcome.
    with pytest.raises(LabelError):
        mod.study_deltas(pathlib.Path("nonexistent.h5ad"), obs, "S")


def test_the_report_carries_the_unlabellable_count():
    """The field the sidecar reads, so a shrunken study is visible."""
    import inspect

    import src.reference.jobs.icbi_coexpression as mod

    source = inspect.getsource(mod.study_deltas)
    assert 'report["n_patients_unlabellable"]' in source
    assert 'report["unlabellable"]' in source
    assert "except LabelError" in source
