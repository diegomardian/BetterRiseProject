"""Bounded Synapse VCF-header capability check for exact Chen WES candidates.

This is Step 2 of ``docs/wes_subtype_plan.md``.  It authenticates only for the
35 specimen-exact VCF entities recorded by the committed Step-1 crosswalk,
requests a byte range capped at 262,144 bytes, and stops at the ``#CHROM``
line.  It does not download a VCF, inspect a variant record, select BRAF/APC
rules, or calculate a transcript subgroup result.

    python -m pip install -e '.[a2]'
    synapse config
    python -m src.reference.jobs.wes_subtype_vcf_header --no-write

Run without ``--no-write`` only from a clean, committed tree to write the
technical-capability artifact.  The resulting schema informs the later
pre-specification; it is not a pass to read variant records.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

from src.common.io import read_versioned_table, write_versioned_table
from src.common.label_provenance import Measurement, check_no_circular_claim, provenance_meta
from src.common.paths import RESULTS_DIR
from src.common.provenance import DEFAULT_SEED
from src.reference.a2_synapse import login_synapse
from src.reference.wes_vcf_header import (
    HEADER_BYTE_BUDGET,
    WESHeaderError,
    exact_vcf_candidates,
    parse_vcf_header,
)

DEFAULT_CROSSWALK = (
    RESULTS_DIR / "2026-09-09_0bf9734" / "wes_subtype_provenance_crosswalk.parquet"
)
EXPECTED_EXACT_CANDIDATES = 35

LABEL_PROVENANCE = Measurement(
    modality="sample_annotation",
    assay="HTAN exact scRNA-to-WES biospecimen provenance crosswalk",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="genotype",
    assay="bounded header schema of exact-match HTAN Release-7 WES VCFs",
    genes=(),
)


def validate_specification() -> tuple[str, ...]:
    return check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)


def _read_crosswalk(path: Path) -> pd.DataFrame:
    table, meta = read_versioned_table(path)
    if meta.get("vcf_content_read") is not False:
        raise WESHeaderError("Step-1 crosswalk is not the expected metadata-only artifact")
    candidates = exact_vcf_candidates(table)
    if len(candidates) != EXPECTED_EXACT_CANDIDATES:
        raise WESHeaderError(
            f"expected {EXPECTED_EXACT_CANDIDATES} exact WES candidates, found {len(candidates)}"
        )
    return candidates


def _presigned_url(syn: object, entity_id: str) -> str:
    """Get a signed download URL without invoking Synapse's full-download path."""
    try:
        from synapseclient.api.file_services import get_file_handle_for_download
    except ImportError as exc:  # pragma: no cover - optional operational dependency
        raise WESHeaderError(
            "synapseclient is not installed; run `pip install -e '.[a2]'`"
        ) from exc
    entity = syn.get(entity_id, downloadFile=False)
    properties = getattr(entity, "properties", {}) or {}
    handle_id = properties.get("dataFileHandleId") or getattr(entity, "dataFileHandleId", None)
    if not handle_id:
        raise WESHeaderError(f"Synapse entity {entity_id} has no dataFileHandleId")
    result = get_file_handle_for_download(
        str(handle_id), entity_id, synapse_client=syn  # type: ignore[arg-type]
    )
    url = result.get("preSignedURL")
    if not isinstance(url, str) or not url:
        raise WESHeaderError(f"Synapse did not return a signed URL for {entity_id}")
    return url


def bounded_vcf_header(url: str, *, byte_budget: int = HEADER_BYTE_BUDGET) -> tuple[list[str], int]:
    """Stream only through the terminating ``#CHROM`` line of a ranged response."""
    request = Request(
        url,
        headers={
            "Range": f"bytes=0-{byte_budget - 1}",
            "Accept-Encoding": "identity",
        },
    )
    with urlopen(request, timeout=60) as response:  # noqa: S310 - Synapse signed URL
        content_range = response.headers.get("Content-Range", "")
        if response.status != 206 or not content_range.startswith("bytes 0-"):
            raise WESHeaderError("Synapse storage did not honor the bounded byte-range request")
        lines: list[str] = []
        current = bytearray()
        bytes_read = 0
        while bytes_read < byte_budget:
            chunk = response.read(1)
            if not chunk:
                break
            bytes_read += 1
            current.extend(chunk)
            if chunk != b"\n":
                continue
            try:
                line = current.decode("utf-8")
            except UnicodeDecodeError as exc:
                raise WESHeaderError("VCF header is not UTF-8 text") from exc
            lines.append(line)
            current.clear()
            if line.startswith("#CHROM"):
                return lines, bytes_read
        raise WESHeaderError("bounded read ended before a newline-terminated #CHROM header")


def probe_header(syn: object, row: pd.Series) -> dict[str, object]:
    """Produce a result row; access or schema failure remains explicit evidence."""
    base = {
        "patient_id": row["patient_id"],
        "scRNA_biospecimen_id": row["scRNA_biospecimen_id"],
        "vcf_data_file_id": row["matched_vcf_data_file_ids"],
        "vcf_entity_id": row["matched_vcf_entity_ids"],
        "header_byte_budget": HEADER_BYTE_BUDGET,
    }
    try:
        lines, bytes_read = bounded_vcf_header(_presigned_url(syn, str(base["vcf_entity_id"])))
        return {
            **base,
            "header_status": "parsed",
            "header_bytes_read": bytes_read,
            "error": "",
            **parse_vcf_header(lines),
        }
    except Exception as exc:  # pragma: no cover - depends on authenticated remote state
        return {
            **base,
            "header_status": "unresolved",
            "header_bytes_read": pd.NA,
            "error": str(exc),
            "vcf_fileformat": "",
            "reference_declarations": "",
            "caller_declarations": "",
            "filter_declaration_ids": "",
            "info_declaration_ids": "",
            "format_declaration_ids": "",
            "n_sample_columns": pd.NA,
            "sample_layout": "unresolved",
            "matched_normal_status": "unresolved",
            "has_dp_info": pd.NA,
            "has_dp_format": pd.NA,
            "has_ad_format": pd.NA,
            "has_af_info": pd.NA,
            "has_af_format": pd.NA,
        }


def probe_headers(syn: object, candidates: pd.DataFrame) -> pd.DataFrame:
    """Probe every exact candidate independently; one failure cannot disappear."""
    return pd.DataFrame.from_records([probe_header(syn, row) for _, row in candidates.iterrows()])


def header_summary(result: pd.DataFrame) -> pd.DataFrame:
    parsed = result["header_status"].eq("parsed")
    layout_counts = result["sample_layout"].value_counts(sort=False).sort_index()
    return pd.DataFrame(
        [
            {
                "n_specimen_exact_vcf_candidates": len(result),
                "n_unique_patients": result["patient_id"].nunique(),
                "n_headers_parsed": int(parsed.sum()),
                "n_headers_unresolved": int((~parsed).sum()),
                "n_with_dp_field": int(
                    (result.loc[parsed, "has_dp_info"].eq(True)
                    | result.loc[parsed, "has_dp_format"].eq(True)).sum()
                ),
                "n_with_ad_or_af_field": int(
                    (result.loc[parsed, "has_ad_format"].eq(True)
                    | result.loc[parsed, "has_af_info"].eq(True)
                    | result.loc[parsed, "has_af_format"].eq(True)).sum()
                ),
                "sample_layout_rollup": " | ".join(
                    f"{layout}:{count}" for layout, count in layout_counts.items()
                ),
                "verdict": (
                    "HEADER CAPABILITY OBSERVED — lock Step 3 before variant records"
                    if bool(parsed.all())
                    else "CONTENT ACCESS OR HEADER CAPABILITY UNRESOLVED — do not read variants"
                ),
            }
        ]
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosswalk", type=Path, default=DEFAULT_CROSSWALK)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    validate_specification()
    result = probe_headers(login_synapse(), _read_crosswalk(args.crosswalk))
    summary = header_summary(result)
    print(summary.to_string(index=False))
    print(
        "SYNAPSE VCF HEADER ONLY — ranged bytes through #CHROM; "
        "no variant record, molecular label, or transcript subgroup was read."
    )
    if args.no_write:
        return 0
    meta = {
        "plan": "docs/wes_subtype_plan.md",
        "candidate_rule": "Step-1 specimen_exact_unique only",
        "header_byte_budget_per_entity": HEADER_BYTE_BUDGET,
        "vcf_variant_records_read": False,
        "molecular_endpoint_selected": False,
        "subgroup_outcome_calculated": False,
        **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
    }
    for table, name in (
        (result, "wes_subtype_vcf_header"),
        (summary, "wes_subtype_vcf_header_summary"),
    ):
        path = write_versioned_table(
            table,
            name,
            seed=args.seed,
            results_dir=args.results_dir,
            allow_dirty=args.allow_dirty,
            extra_meta=meta,
        )
        print(f"wrote {path}")
    return 0 if result["header_status"].eq("parsed").all() else 5


if __name__ == "__main__":
    sys.exit(main())
