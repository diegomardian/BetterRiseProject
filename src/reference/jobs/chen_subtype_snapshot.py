"""Pin the open cBioPortal HTA11 responses that `prereg_chen_lesion_subtype` reads.

    python -m src.reference.jobs.chen_subtype_snapshot --no-write

Gate condition 3 of `docs/prereg_chen_lesion_subtype.md` §9. A public API is a
moving target: cBioPortal restates a study whenever the deposit is revised, and
an analysis that queries it live cannot say afterwards which version it read.
This job writes the three responses to ``data/raw/cbioportal_hta11/`` and emits
the checksum table that goes into ``data/manifest.csv``.

It reads no decomposition value and joins nothing to a label. Pinning the
labels is the step that has to happen *before* anyone can be tempted to look.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from src.common.io import write_versioned_table
from src.common.label_provenance import (
    Measurement,
    check_no_circular_claim,
    provenance_meta,
)
from src.common.paths import RAW_DIR
from src.common.provenance import DEFAULT_SEED

STUDY_ID = "crc_hta11_htan_2021"
API = "https://www.cbioportal.org/api"
SNAPSHOT_DIR = RAW_DIR / "cbioportal_hta11"

#: Invariant 11. The labels are a pathologist's reading of morphology; the
#: claim this snapshot will later serve is a transcript decomposition. Declared
#: here because the population is chosen the moment these rows are frozen.
LABEL_PROVENANCE = Measurement(
    modality="morphology",
    assay="Chen 2021 pathologist polyp diagnosis via cBioPortal HTA11 clinical attributes",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="transcript",
    assay="Chen 2021 snRNA-seq adenoma decomposition (lineage rung)",
    genes=(),
)


class SnapshotError(RuntimeError):
    """The pinned snapshot cannot be established as written."""


def _get(path: str, *, timeout: int) -> bytes:
    with urllib.request.urlopen(API + path, timeout=timeout) as response:  # noqa: S310
        return response.read()


def _post(path: str, payload: dict[str, Any], *, timeout: int) -> bytes:
    request = urllib.request.Request(
        API + path,
        data=json.dumps(payload).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:  # noqa: S310
        return response.read()


def fetch_snapshot(*, timeout: int = 240) -> dict[str, bytes]:
    """The three responses the pre-registration names, and nothing else."""
    return {
        "study.json": _get(f"/studies/{STUDY_ID}", timeout=timeout),
        "clinical_sample.json": _get(
            f"/studies/{STUDY_ID}/clinical-data"
            "?clinicalDataType=SAMPLE&projection=DETAILED&pageSize=5000",
            timeout=timeout,
        ),
        "mutations.json": _post(
            f"/molecular-profiles/{STUDY_ID}_mutations/mutations/fetch?projection=DETAILED",
            {"sampleListId": f"{STUDY_ID}_all"},
            timeout=timeout,
        ),
    }


def _counts(name: str, raw: bytes) -> dict[str, Any]:
    """Row counts a later reader can check the snapshot against."""
    payload = json.loads(raw)
    if name == "study.json":
        return {
            "n_records": 1,
            "n_samples": payload.get("allSampleCount"),
            "n_patients": None,
        }
    if name == "clinical_sample.json":
        return {
            "n_records": len(payload),
            "n_samples": len({row["sampleId"] for row in payload}),
            "n_patients": len({row["patientId"] for row in payload}),
        }
    return {
        "n_records": len(payload),
        "n_samples": len({row["sampleId"] for row in payload}),
        "n_patients": None,
    }


def snapshot_table(responses: dict[str, bytes], *, written: bool) -> pd.DataFrame:
    rows = []
    for name in sorted(responses):
        raw = responses[name]
        # Refuse before parsing: an empty body is a failed fetch, not a document.
        if not raw:
            raise SnapshotError(f"{name} came back empty; that is not a snapshot")
        rows.append(
            {
                "file": name,
                "path": f"data/raw/cbioportal_hta11/{name}",
                "sha256": hashlib.sha256(raw).hexdigest(),
                "bytes": len(raw),
                "written": written,
                **_counts(name, raw),
            }
        )
    return pd.DataFrame(rows)


def manifest_rows(frame: pd.DataFrame, *, downloaded_on: str) -> pd.DataFrame:
    """The `data/manifest.csv` rows, emitted rather than appended blindly."""
    source = {
        "study.json": f"{API}/studies/{STUDY_ID}",
        "clinical_sample.json": (
            f"{API}/studies/{STUDY_ID}/clinical-data"
            "?clinicalDataType=SAMPLE&projection=DETAILED&pageSize=5000"
        ),
        "mutations.json": (
            f"{API}/molecular-profiles/{STUDY_ID}_mutations/mutations/fetch"
            "?projection=DETAILED [POST sampleListId=" + STUDY_ID + "_all]"
        ),
    }
    return pd.DataFrame(
        {
            "path": frame["path"],
            "sha256": frame["sha256"],
            "bytes": frame["bytes"],
            "source_url": frame["file"].map(source),
            "accession": STUDY_ID,
            "downloaded_on": downloaded_on,
            "downloaded_by": "bode",
            "workstream": "W1",
            "notes": "prereg_chen_lesion_subtype gate 3; open cBioPortal snapshot; "
            "labels only, no decomposition value read",
        }
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--timeout", type=int, default=240)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    # Invariant 11: refuse before the first read, not at the writer.
    check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)

    responses = fetch_snapshot(timeout=args.timeout)
    if not args.no_write:
        SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
        for name, raw in responses.items():
            (SNAPSHOT_DIR / name).write_bytes(raw)

    frame = snapshot_table(responses, written=not args.no_write)
    print(frame.to_string(index=False))
    print("\nmanifest rows:")
    rows = manifest_rows(frame, downloaded_on="2026-09-10")
    print(rows.to_csv(index=False, header=False), end="")
    print("LABELS ONLY — no decomposition value was read and nothing was joined.")

    if args.no_write:
        return 0
    path = write_versioned_table(
        frame,
        "chen_subtype_cbioportal_snapshot",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta={
            **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
            "prereg": "docs/prereg_chen_lesion_subtype.md",
            "gate_condition": 3,
            "study_id": STUDY_ID,
            "decomposition_value_read": False,
            "label_joined_to_value": False,
        },
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
