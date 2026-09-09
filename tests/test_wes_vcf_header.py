"""The WES header probe must stop before a variant body can be interpreted."""

from __future__ import annotations

import pandas as pd
import pytest

from src.reference.jobs.wes_subtype_vcf_header import (
    _download_eligibility,
    _header_verdict,
    _presigned_url,
    validate_specification,
)
from src.reference.wes_vcf_header import WESHeaderError, exact_vcf_candidates, parse_vcf_header


def _crosswalk() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "patient_id": ["P1"],
            "scRNA_biospecimen_id": ["HTA11_1_2000001011"],
            "specimen_exact": [True],
            "crosswalk_status": ["specimen_exact_unique"],
            "matched_vcf_entity_ids": ["syn1"],
            "matched_vcf_data_file_ids": ["HTA11_0_1"],
        }
    )


def test_header_probe_declares_the_independent_genotype_channel_before_access():
    assert validate_specification() == ()


def test_signed_url_lookup_uses_a_metadata_bundle_before_the_file_handle_service(monkeypatch):
    class FakeSynapse:
        def restPOST(self, uri, body):
            assert uri == "/entity/syn1/bundle2"
            assert '"includeEntity": true' in body
            return {"entity": {"dataFileHandleId": "99"}}

    def fake_handle(handle_id, entity_id, **kwargs):
        assert (handle_id, entity_id) == ("99", "syn1")
        return {"preSignedURL": "https://example.test/ranged"}

    monkeypatch.setattr(
        "synapseclient.api.file_services.get_file_handle_for_download", fake_handle
    )
    assert _presigned_url(FakeSynapse(), "syn1") == "https://example.test/ranged"


def test_download_eligibility_separates_account_certification_from_file_requirements():
    class FakeSynapse:
        def restGET(self, uri):
            assert uri == "/entity/syn1/permissions"
            return {
                "canView": True,
                "canDownload": False,
                "isCertifiedUser": False,
                "isCertificationRequired": True,
            }

        def restPOST(self, uri, body):
            return {"restrictionInformation": {"hasUnmetAccessRequirement": False}}

    assert _download_eligibility(FakeSynapse(), "syn1") == {
        "can_view": True,
        "can_download": False,
        "is_certified_user": False,
        "is_certification_required": True,
        "has_unmet_access_requirement": False,
    }


def test_all_download_ineligible_candidates_receive_a_specific_nonvariant_verdict():
    assert _header_verdict(
        parsed=pd.Series([False, False]), ineligible=pd.Series([True, True])
    ).startswith("ACCOUNT NOT DOWNLOAD-ELIGIBLE")


def test_exact_candidates_refuse_joined_or_ambiguous_entity_identifiers():
    crosswalk = _crosswalk()
    crosswalk.loc[0, "matched_vcf_entity_ids"] = "syn1 | syn2"
    with pytest.raises(WESHeaderError, match="one unjoined identifier"):
        exact_vcf_candidates(crosswalk)


def test_exact_candidates_refuse_nonexact_row_even_if_it_has_an_entity():
    crosswalk = _crosswalk()
    crosswalk.loc[0, "crosswalk_status"] = "participant_vcf_not_specimen_exact"
    with pytest.raises(WESHeaderError, match="specimen_exact_unique"):
        exact_vcf_candidates(crosswalk)


def test_header_parser_records_schema_but_never_pretends_two_columns_are_matched_normal():
    result = parse_vcf_header(
        [
            "##fileformat=VCFv4.2\n",
            "##reference=GRCh38\n",
            "##source=GATK\n",
            "##FILTER=<ID=LowQual,Description=low quality>\n",
            "##INFO=<ID=AF,Number=A,Type=Float,Description=allele frequency>\n",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=depth>\n",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=allelic depths>\n",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
        ]
    )
    assert result["vcf_fileformat"] == "VCFv4.2"
    assert result["reference_declarations"] == "GRCh38"
    assert result["has_dp_format"] and result["has_ad_format"] and result["has_af_info"]
    assert result["n_sample_columns"] == 2
    assert result["sample_layout"] == "two_sample_header_roles_unresolved"
    assert result["matched_normal_status"] == "not_determined_from_header"


def test_header_parser_refuses_body_row_before_or_in_place_of_header_terminator():
    with pytest.raises(WESHeaderError, match="variant record"):
        parse_vcf_header(["##fileformat=VCFv4.2\n", "1\t1\t.\tA\tT\t.\tPASS\t.\n"])


def test_header_parser_refuses_no_header_terminator():
    with pytest.raises(WESHeaderError, match="#CHROM"):
        parse_vcf_header(["##fileformat=VCFv4.2\n"])
