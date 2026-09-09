"""Pure helpers for the bounded, header-only Chen WES capability check.

This module deliberately knows nothing about BRAF, APC, genomic positions, or
variant records.  Step 2 is allowed to establish only whether an exact-match
VCF can be opened and which technical fields its header declares.  The label
rules are a later, separately locked step.
"""

from __future__ import annotations

import re
from collections.abc import Iterable

import pandas as pd

HEADER_BYTE_BUDGET = 262_144


class WESHeaderError(RuntimeError):
    """The bounded capability check cannot make a trustworthy header claim."""


def exact_vcf_candidates(crosswalk: pd.DataFrame) -> pd.DataFrame:
    """Return exactly one VCF entity per Step-1 specimen-exact lesion.

    A joined ID, a duplicate entity, or a non-exact row would make an
    authenticated header result look more lesion-specific than it is.  Refuse
    instead of choosing one of several possible WES objects.
    """
    required = {
        "patient_id",
        "scRNA_biospecimen_id",
        "specimen_exact",
        "crosswalk_status",
        "matched_vcf_entity_ids",
        "matched_vcf_data_file_ids",
    }
    missing = sorted(required - set(crosswalk.columns))
    if missing:
        raise WESHeaderError(f"Step-1 crosswalk is missing required columns {missing}")

    selected = crosswalk.loc[crosswalk["specimen_exact"].eq(True)].copy()
    if selected.empty:
        raise WESHeaderError("Step-1 crosswalk contains no specimen-exact WES candidates")
    if not selected["crosswalk_status"].eq("specimen_exact_unique").all():
        raise WESHeaderError("only specimen_exact_unique rows may reach the header probe")
    for column in ("matched_vcf_entity_ids", "matched_vcf_data_file_ids"):
        values = selected[column].fillna("").astype(str).str.strip()
        if values.eq("").any() or values.str.contains(r"\s\|\s", regex=True).any():
            raise WESHeaderError(f"{column} must contain one unjoined identifier per candidate")
        selected[column] = values
    if selected["scRNA_biospecimen_id"].duplicated().any():
        raise WESHeaderError("a scRNA biospecimen appears more than once in the exact-match arm")
    if selected["matched_vcf_entity_ids"].duplicated().any():
        raise WESHeaderError("one VCF entity is linked to multiple exact-match biospecimens")
    return selected.sort_values(
        ["patient_id", "scRNA_biospecimen_id", "matched_vcf_entity_ids"], ignore_index=True
    )


_DECLARATION_ID = re.compile(r"(?:<|,)ID=([^,>]+)")


def _declaration_id(line: str) -> str:
    match = _DECLARATION_ID.search(line)
    return match.group(1) if match else ""


def _joined(values: Iterable[str]) -> str:
    return " | ".join(sorted({value for value in values if value}))


def _sample_layout(n_sample_columns: int) -> str:
    if n_sample_columns == 0:
        return "no_sample_columns"
    if n_sample_columns == 1:
        return "single_sample_header"
    if n_sample_columns == 2:
        return "two_sample_header_roles_unresolved"
    return "multiple_sample_header_roles_unresolved"


def parse_vcf_header(lines: Iterable[str]) -> dict[str, object]:
    """Parse VCF declarations through ``#CHROM`` and reject every body line.

    ``lines`` must end immediately after the header terminator.  A body row is
    not harmless here: allowing it would erase the promised boundary between
    technical reconnaissance and outcome-bearing variant inspection.
    """
    fileformat = ""
    references: list[str] = []
    sources: list[str] = []
    filters: list[str] = []
    infos: list[str] = []
    formats: list[str] = []
    header_columns: list[str] | None = None

    for line in lines:
        if line.startswith("##fileformat="):
            fileformat = line.removeprefix("##fileformat=").strip()
        elif line.startswith("##reference="):
            references.append(line.removeprefix("##reference=").strip())
        elif line.startswith("##source="):
            sources.append(line.removeprefix("##source=").strip())
        elif line.startswith("##FILTER=<"):
            filters.append(_declaration_id(line))
        elif line.startswith("##INFO=<"):
            infos.append(_declaration_id(line))
        elif line.startswith("##FORMAT=<"):
            formats.append(_declaration_id(line))
        elif line.startswith("#CHROM"):
            header_columns = line.rstrip("\r\n").split("\t")
            break
        elif line.startswith("#"):
            continue
        else:
            raise WESHeaderError("variant record encountered before #CHROM header terminator")

    if header_columns is None:
        raise WESHeaderError("bounded read ended before a #CHROM header terminator")
    fixed = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    if header_columns[:8] != fixed:
        raise WESHeaderError("#CHROM header does not contain the required VCF fixed columns")
    if len(header_columns) > 8 and header_columns[8] != "FORMAT":
        raise WESHeaderError("VCF sample columns appear without a FORMAT column")

    n_samples = max(len(header_columns) - 9, 0)
    info_ids = _joined(infos)
    format_ids = _joined(formats)
    return {
        "vcf_fileformat": fileformat,
        "reference_declarations": _joined(references),
        "caller_declarations": _joined(sources),
        "filter_declaration_ids": _joined(filters),
        "info_declaration_ids": info_ids,
        "format_declaration_ids": format_ids,
        "n_sample_columns": n_samples,
        "sample_layout": _sample_layout(n_samples),
        # A header's sample names and count cannot safely establish which, if
        # any, sample is matched normal.  Preserve that uncertainty explicitly.
        "matched_normal_status": "not_determined_from_header",
        "has_dp_info": "DP" in infos,
        "has_dp_format": "DP" in formats,
        "has_ad_format": "AD" in formats,
        "has_af_info": "AF" in infos,
        "has_af_format": "AF" in formats,
    }
