"""Bounded, authenticated TIFF-header check for candidate HTA11 H&E images.

The job resolves each exact premalignant H&E candidate's open Synapse entity,
then requests no more than 65,536 bytes from the signed file URL. It extracts
first-page TIFF dimensions and physical pixel scale when those tags occur in
that prefix. It never downloads an image, opens a VCF, selects an endpoint, or
fits a model.

    pip install -e '.[a2]'
    synapse config
    python -m src.reference.jobs.he_molecular_image_headers \
      --files /path/to/files.tsv \
      --biospecimens /path/to/biospecimens.tsv \
      --cases /path/to/cases.tsv \
      --no-write
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any
from urllib.request import Request, urlopen

import pandas as pd

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.reference.a2_synapse import login_synapse
from src.reference.he_molecular_gate import build_inventory
from src.reference.he_tiff_header import TiffHeaderError, parse_tiff_header

HEADER_BYTES = 65_536


def _read(path: Path, *, name: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{name} metadata not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def candidate_images(attrition: pd.DataFrame) -> pd.DataFrame:
    """Return only exact premalignant image entities, never all atlas H&E files."""
    rows = attrition[attrition["has_exact_level3_vcf"] & attrition["candidate_premalignant"]]
    output = rows.loc[:, ["biospecimen_id", "participant_id", "image_synapse_id"]].rename(
        columns={"image_synapse_id": "entity_id"}
    )
    if output.empty or output["entity_id"].eq("").any():
        raise TiffHeaderError("no exact premalignant H&E candidate image entities")
    if output.duplicated(["biospecimen_id", "entity_id"]).any():
        raise TiffHeaderError("candidate image extraction emitted a duplicate row")
    return output.sort_values(["biospecimen_id", "entity_id"], ignore_index=True)


def _header_prefix(url: str) -> bytes:
    request = Request(
        url, headers={"Range": f"bytes=0-{HEADER_BYTES - 1}", "Accept-Encoding": "identity"}
    )
    with urlopen(request, timeout=60) as response:  # noqa: S310 - Synapse-signed URL
        return response.read(HEADER_BYTES)


def probe_headers(syn: Any, candidates: pd.DataFrame) -> pd.DataFrame:
    """Return per-image header facts; one unparseable entity remains explicit."""
    try:
        from synapseclient.api.file_services import get_file_handle_for_download
        from synapseclient.operations import FileOptions, get
    except ImportError as exc:  # pragma: no cover - depends on the optional extra
        raise RuntimeError("synapseclient is not installed; run `pip install -e '.[a2]'`") from exc

    rows: list[dict[str, object]] = []
    for candidate in candidates.itertuples(index=False):
        base = candidate._asdict()
        try:
            entity = get(
                synapse_id=candidate.entity_id,
                file_options=FileOptions(download_file=False),
                synapse_client=syn,
            )
            file_handle_id = getattr(entity, "data_file_handle_id", None)
            if not file_handle_id:
                raise TiffHeaderError("Synapse file entity has no dataFileHandleId")
            download = get_file_handle_for_download(
                str(file_handle_id), candidate.entity_id, synapse_client=syn
            )
            handle = download["fileHandle"]
            header = parse_tiff_header(_header_prefix(download["preSignedURL"]))
            rows.append(
                {
                    **base,
                    "header_status": "parsed",
                    "file_size_bytes": int(handle.get("contentSize", 0)),
                    "header_bytes_requested": HEADER_BYTES,
                    "error": "",
                    **header,
                }
            )
        except Exception as exc:  # pragma: no cover - depends on authenticated remote state
            rows.append(
                {
                    **base,
                    "header_status": "unresolved",
                    "file_size_bytes": pd.NA,
                    "header_bytes_requested": HEADER_BYTES,
                    "error": str(exc),
                    "tiff_kind": "",
                    "pixel_width": pd.NA,
                    "pixel_height": pd.NA,
                    "resolution_unit": "",
                    "x_microns_per_pixel": pd.NA,
                    "y_microns_per_pixel": pd.NA,
                }
            )
    return pd.DataFrame(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--files", required=True, type=Path)
    parser.add_argument("--biospecimens", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    attrition, _ = build_inventory(
        _read(args.files, name="Files"),
        _read(args.biospecimens, name="Biospecimen"),
        _read(args.cases, name="Case"),
    )
    result = probe_headers(login_synapse(), candidate_images(attrition))
    parsed = result["header_status"].eq("parsed")
    print(result.to_string(index=False))
    print(
        "TIFF HEADER ONLY — no full image or VCF was downloaded. "
        f"candidate headers parsed: {int(parsed.sum())}/{len(result)}"
    )
    if args.no_write:
        return 0
    path = write_versioned_table(
        result,
        "he_molecular_image_headers",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta={
            "prereg": "docs/he_molecular_data_gate.md",
            "header_bytes_requested_per_entity": HEADER_BYTES,
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
