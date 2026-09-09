"""Read bounded TIFF headers from exact HTA11 H&E--VCF--CRDC objects.

This job accepts either the one-object ``gen3_manifest.json`` produced by the
HTAN portal or the fixed 18-row BigQuery crosswalk exported as CSV. It uses the
official Gen3 SDK only to request a short-lived signed URL, then requests no
more than 65,536 bytes per entity across TIFF header regions the file declares.
It never downloads an image, opens a VCF, selects an endpoint, or fits a model.

The result is either a one-object pilot or an all-candidate header audit. It
does not itself establish that a deposited TIFF is the scanner's original
full-resolution object.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.reference.he_crdc import (
    read_exact_he_crosswalk,
    signed_url_from_response,
    single_object_id,
)
from src.reference.he_tiff_header import parse_tiff_header_ranges

CRDC_ENDPOINT = "https://nci-crdc.datacommons.io"
HEADER_BYTES = 65_536


def crdc_signed_url(*, object_id: str, credentials: Path) -> str:
    """Request one short-lived URL without emitting the credential or URL."""
    if not credentials.is_file():
        raise FileNotFoundError(f"CRDC credentials not found: {credentials}")
    try:
        from gen3.auth import Gen3Auth
        from gen3.file import Gen3File
    except ImportError as exc:  # pragma: no cover - optional operational dependency
        message = "Gen3 SDK is not installed; run `python -m pip install -e '.[he]'`"
        raise RuntimeError(message) from exc
    auth = Gen3Auth(CRDC_ENDPOINT, refresh_file=str(credentials))
    response = Gen3File(CRDC_ENDPOINT, auth).get_presigned_url(object_id)
    return signed_url_from_response(response)


def header_range(url: str, offset: int, length: int) -> bytes:
    """Request one declared TIFF region, without entering the full-download path."""
    request = Request(
        url,
        headers={"Range": f"bytes={offset}-{offset + length - 1}", "Accept-Encoding": "identity"},
    )
    with urlopen(request, timeout=60) as response:  # noqa: S310 - short-lived CRDC signed URL
        content_range = response.headers.get("Content-Range", "")
        expected = f"bytes {offset}-"
        if response.status != 206 or not content_range.startswith(expected):
            raise RuntimeError("CRDC storage did not honor the requested TIFF byte range")
        return response.read(length)


def probe_object(*, object_id: str, biospecimen_id: str, credentials: Path) -> pd.DataFrame:
    """Return one explicit row; remote failure is unresolved, not absent, image."""
    base: dict[str, object] = {"biospecimen_id": biospecimen_id, "crdc_object_id": object_id}
    try:
        url = crdc_signed_url(object_id=object_id, credentials=credentials)
        header, bytes_requested = parse_tiff_header_ranges(
            lambda offset, length: header_range(url, offset, length),
            byte_budget=HEADER_BYTES,
        )
        return pd.DataFrame(
            [
                {
                    **base,
                    "header_status": "parsed",
                    "header_bytes_requested": bytes_requested,
                    "header_byte_budget": HEADER_BYTES,
                    "error": "",
                    **header,
                }
            ]
        )
    except Exception as exc:  # pragma: no cover - depends on authenticated remote state
        return pd.DataFrame(
            [
                {
                    **base,
                    "header_status": "unresolved",
                    "header_bytes_requested": pd.NA,
                    "header_byte_budget": HEADER_BYTES,
                    "error": str(exc),
                    "tiff_kind": "",
                    "pixel_width": pd.NA,
                    "pixel_height": pd.NA,
                    "resolution_unit": "",
                    "x_microns_per_pixel": pd.NA,
                    "y_microns_per_pixel": pd.NA,
                }
            ]
        )


def probe(*, manifest: Path, biospecimen_id: str, credentials: Path) -> pd.DataFrame:
    """Probe the one object selected through the HTAN portal."""
    return probe_object(
        object_id=single_object_id(manifest),
        biospecimen_id=biospecimen_id,
        credentials=credentials,
    )


def probe_crosswalk(*, mapping_csv: Path, credentials: Path) -> pd.DataFrame:
    """Probe every member of the fixed exact-match crosswalk, independently."""
    crosswalk = read_exact_he_crosswalk(mapping_csv)
    rows = []
    for row in crosswalk.itertuples(index=False):
        result = probe_object(
            object_id=row.drs_uri,
            biospecimen_id=row.biospecimen_id,
            credentials=credentials,
        )
        result.insert(1, "HTAN_Data_File_ID", row.HTAN_Data_File_ID)
        rows.append(result)
    records = [record for row in rows for record in row.to_dict(orient="records")]
    return pd.DataFrame.from_records(records)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--manifest", type=Path)
    source.add_argument("--crosswalk-csv", type=Path)
    parser.add_argument("--biospecimen-id")
    parser.add_argument("--credentials", required=True, type=Path)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    if args.manifest:
        if not args.biospecimen_id:
            parser.error("--biospecimen-id is required with --manifest")
        result = probe(
            manifest=args.manifest,
            biospecimen_id=args.biospecimen_id,
            credentials=args.credentials,
        )
    else:
        if args.biospecimen_id:
            parser.error("--biospecimen-id is only valid with --manifest")
        result = probe_crosswalk(mapping_csv=args.crosswalk_csv, credentials=args.credentials)
    parsed = result["header_status"].eq("parsed")
    print(result.to_string(index=False))
    print(
        "CRDC TIFF HEADER ONLY — no full image or VCF was downloaded. "
        f"headers parsed: {int(parsed.sum())}/{len(parsed)}"
    )
    if args.no_write:
        return 0
    path = write_versioned_table(
        result,
        "he_molecular_crdc_header",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta={
            "prereg": "docs/he_molecular_data_gate.md",
            "header_byte_budget_per_entity": HEADER_BYTES,
            "full_image_downloaded": False,
            "vcf_content_downloaded": False,
            "endpoint_selected": False,
            "all_candidate_headers_parsed": bool(parsed.all()),
        },
    )
    print(f"wrote {path}")
    return 0 if parsed.all() else 5


if __name__ == "__main__":
    sys.exit(main())
