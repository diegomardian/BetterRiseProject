from __future__ import annotations

import json

import pandas as pd
import pytest

from src.reference.he_crdc import (
    N_EXACT_HE_VCF_CANDIDATES,
    CrcdCrosswalkError,
    CrcdManifestError,
    read_exact_he_crosswalk,
    signed_url_from_response,
    single_object_id,
)
from src.reference.jobs import he_molecular_crdc_header


def test_single_object_id_reads_portal_manifest(tmp_path) -> None:
    manifest = tmp_path / "gen3_manifest.json"
    manifest.write_text(json.dumps([{"object_id": "dg.4DFC/abc"}]), encoding="utf-8")
    assert single_object_id(manifest) == "dg.4DFC/abc"


@pytest.mark.parametrize("payload", [[], [{"object_id": "one"}, {"object_id": "two"}], [{}]])
def test_single_object_id_refuses_ambiguous_or_missing_manifest(tmp_path, payload) -> None:
    manifest = tmp_path / "gen3_manifest.json"
    manifest.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(CrcdManifestError):
        single_object_id(manifest)


def test_signed_url_from_response_accepts_documented_url_field() -> None:
    assert signed_url_from_response({"url": "https://signed.example/object"}) == "https://signed.example/object"


def test_signed_url_from_response_does_not_accept_ambiguous_response() -> None:
    with pytest.raises(CrcdManifestError, match="exactly one HTTPS URL"):
        signed_url_from_response({"first": "https://one.example", "second": "https://two.example"})


def _exact_crosswalk() -> pd.DataFrame:
    biospecimens = [f"HTA11_{index}_2000001011" for index in range(N_EXACT_HE_VCF_CANDIDATES)]
    return pd.DataFrame(
        {
            "biospecimen_id": biospecimens,
            "HTAN_Data_File_ID": [f"file-{index}" for index in range(N_EXACT_HE_VCF_CANDIDATES)],
            "HTAN_Assayed_Biospecimen_ID": biospecimens,
            "HTAN_Originating_Biospecimen_ID": biospecimens,
            "drs_uri": [f"dg.4DFC/{index}" for index in range(N_EXACT_HE_VCF_CANDIDATES)],
        }
    )


def test_read_exact_he_crosswalk_requires_fixed_unique_exact_rows(tmp_path) -> None:
    path = tmp_path / "crosswalk.csv"
    expected = _exact_crosswalk()
    expected.to_csv(path, index=False)
    actual = read_exact_he_crosswalk(path)
    pd.testing.assert_frame_equal(
        actual,
        expected.astype("string").sort_values("biospecimen_id").reset_index(drop=True),
    )


@pytest.mark.parametrize("mutation", ["too_short", "duplicate", "non_exact", "missing_drs"])
def test_read_exact_he_crosswalk_refuses_nonfixed_or_inexact_input(tmp_path, mutation: str) -> None:
    table = _exact_crosswalk()
    if mutation == "too_short":
        table = table.iloc[:-1]
    elif mutation == "duplicate":
        table.loc[1, "biospecimen_id"] = table.loc[0, "biospecimen_id"]
    elif mutation == "non_exact":
        table.loc[0, "HTAN_Assayed_Biospecimen_ID"] = "different"
    else:
        table.loc[0, "drs_uri"] = ""
    path = tmp_path / "crosswalk.csv"
    table.to_csv(path, index=False)
    with pytest.raises(CrcdCrosswalkError):
        read_exact_he_crosswalk(path)


def test_probe_crosswalk_carries_the_exact_h_and_e_file_id(tmp_path, monkeypatch) -> None:
    path = tmp_path / "crosswalk.csv"
    _exact_crosswalk().to_csv(path, index=False)

    def fake_probe_object(*, object_id: str, biospecimen_id: str, credentials):
        return pd.DataFrame(
            [
                {
                    "biospecimen_id": biospecimen_id,
                    "crdc_object_id": object_id,
                    "header_status": "parsed",
                }
            ]
        )

    monkeypatch.setattr(he_molecular_crdc_header, "probe_object", fake_probe_object)
    actual = he_molecular_crdc_header.probe_crosswalk(mapping_csv=path, credentials=tmp_path)
    assert len(actual) == N_EXACT_HE_VCF_CANDIDATES
    assert actual["HTAN_Data_File_ID"].is_unique
    assert actual["crdc_object_id"].str.startswith("dg.").all()
