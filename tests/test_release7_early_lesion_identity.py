from __future__ import annotations

import pytest

from src.reference.jobs.release7_early_lesion_identity import (
    EXPECTED_R7_PARTICIPANTS,
    Release7IdentityError,
    identity_table,
    validate_specification,
)


def _release7_ids() -> tuple[str, ...]:
    return tuple(f"HTA11_R7_{number:03d}" for number in range(EXPECTED_R7_PARTICIPANTS))


def test_identity_audit_declares_gene_free_metadata_measurements():
    assert validate_specification() == ()


def test_disjoint_fixed_universe_is_independent_of_cached_chen():
    table, verdict = identity_table(_release7_ids(), ("HTA11_CHEN_001", "HTA11_CHEN_002"))
    assert verdict == "INDEPENDENT OF CACHED CHEN"
    assert len(table) == EXPECTED_R7_PARTICIPANTS
    assert not table["in_cached_chen"].any()


def test_one_overlap_changes_the_verdict_and_leaves_the_participant_visible():
    release7 = _release7_ids()
    table, verdict = identity_table(release7, ("HTA11_CHEN_001", release7[17]))
    assert verdict == "OVERLAPS CACHED CHEN"
    overlap = table.loc[table["in_cached_chen"]]
    assert overlap["HTAN_Participant_ID"].tolist() == [release7[17]]
    assert overlap["identity_status"].tolist() == ["overlaps_cached_chen"]


def test_identity_audit_refuses_a_silently_changed_source_universe():
    with pytest.raises(Release7IdentityError, match="verified unique 55"):
        identity_table(_release7_ids()[:-1], ())
