"""Audit whether the fixed HTAN Release-7 source overlaps cached Chen IDs.

This reads only a BigQuery participant-ID CSV and the cached ICBI ``obs``
metadata. It never opens a Release-7 H5AD, requests Synapse metadata, or reads
an expression layer.

Run after committing this job, from the repository root:

    python -m src.reference.jobs.release7_early_lesion_identity \
      --participants-csv "$HOME/Downloads/job_<BigQuery job ID>.csv"
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.label_provenance import Measurement, check_no_circular_claim, provenance_meta
from src.common.paths import INTERIM_DIR
from src.common.provenance import DEFAULT_SEED

log = logging.getLogger(__name__)

CHEN_STUDY_ID = "Chen_2021_Cell"
CHEN_ID_PREFIX = f"{CHEN_STUDY_ID}."
EXPECTED_R7_PARTICIPANTS = 55
EXPECTED_CHEN_PARTICIPANTS = 106
PARTICIPANT_COLUMN = "HTAN_Participant_ID"

R7_H5AD_IDS = ("HTA11_0_14002", "HTA11_0_14004")
R7_SYNAPSE_IDS = ("syn53710088", "syn53710094")

LABEL_PROVENANCE = Measurement(
    modality="sample_annotation",
    assay="HTAN Release-7 participant identifiers from BigQuery provenance metadata",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="sample_annotation",
    assay="cached Chen_2021_Cell participant identifiers",
    genes=(),
)


class Release7IdentityError(ValueError):
    """The fixed source-identity audit cannot be evaluated honestly."""


def validate_specification() -> tuple[str, ...]:
    """Declare that this identity audit makes no gene-level claim."""
    return check_no_circular_claim(labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE)


def read_release7_participants(path: Path) -> tuple[str, ...]:
    """Read the one-column, fixed BigQuery export without guessing its schema."""
    if not path.exists():
        raise Release7IdentityError(f"Release-7 participant export not found: {path}")
    frame = pd.read_csv(path)
    if list(frame.columns) != [PARTICIPANT_COLUMN]:
        raise Release7IdentityError(
            f"expected exactly one {PARTICIPANT_COLUMN!r} column, got {list(frame.columns)!r}"
        )
    values = frame[PARTICIPANT_COLUMN]
    if values.isna().any() or (values.astype(str).str.strip() == "").any():
        raise Release7IdentityError("Release-7 participant export contains a missing or blank ID")
    ids = tuple(values.astype(str).str.strip())
    if len(set(ids)) != len(ids):
        raise Release7IdentityError("Release-7 participant export contains duplicated IDs")
    if len(ids) != EXPECTED_R7_PARTICIPANTS:
        raise Release7IdentityError(
            f"Release-7 export has {len(ids)} participants, expected {EXPECTED_R7_PARTICIPANTS}. "
            "Do not silently change the source cohort after inspection."
        )
    return tuple(sorted(ids))


def cached_chen_participants(path: Path) -> tuple[str, ...]:
    """Return all cached Chen participants, including the carcinoma subset."""
    if not path.exists():
        raise Release7IdentityError(f"cached ICBI metadata not found: {path}")
    obs = pd.read_parquet(path, columns=["study_id", "patient_id"])
    patient_ids = obs.loc[obs["study_id"].astype(str) == CHEN_STUDY_ID, "patient_id"]
    if patient_ids.empty:
        raise Release7IdentityError(f"cached metadata contains no rows for {CHEN_STUDY_ID}")
    values = patient_ids.astype(str)
    malformed = sorted(set(values[~values.str.startswith(CHEN_ID_PREFIX)]))
    if malformed:
        raise Release7IdentityError(
            f"cached {CHEN_STUDY_ID} participant IDs lack prefix {CHEN_ID_PREFIX!r}: {malformed}"
        )
    ids = tuple(sorted(set(values.str.removeprefix(CHEN_ID_PREFIX))))
    if len(ids) != EXPECTED_CHEN_PARTICIPANTS:
        raise Release7IdentityError(
            f"cached Chen participant universe has {len(ids)} IDs, expected "
            f"{EXPECTED_CHEN_PARTICIPANTS}. Refuse a source-identity comparison "
            "against a changed cache."
        )
    return ids


def identity_table(
    release7_ids: Sequence[str], chen_ids: Sequence[str]
) -> tuple[pd.DataFrame, str]:
    """Return one auditable row per source participant and the fixed verdict."""
    source = tuple(release7_ids)
    if len(source) != EXPECTED_R7_PARTICIPANTS or len(set(source)) != len(source):
        raise Release7IdentityError(
            "source IDs must be the verified unique 55-participant universe"
        )
    chen = set(chen_ids)
    table = pd.DataFrame({PARTICIPANT_COLUMN: source})
    table["in_cached_chen"] = table[PARTICIPANT_COLUMN].isin(chen)
    table["identity_status"] = table["in_cached_chen"].map({
        False: "not_in_cached_chen", True: "overlaps_cached_chen",
    })
    verdict = (
        "INDEPENDENT OF CACHED CHEN"
        if not table["in_cached_chen"].any()
        else "OVERLAPS CACHED CHEN"
    )
    return table, verdict


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--participants-csv", required=True, type=Path)
    parser.add_argument("--cache", type=Path, default=INTERIM_DIR / "icbi_obs.parquet")
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    validate_specification()
    release7_ids = read_release7_participants(args.participants_csv)
    chen_ids = cached_chen_participants(args.cache)
    table, verdict = identity_table(release7_ids, chen_ids)
    overlapping = table.loc[table["in_cached_chen"], PARTICIPANT_COLUMN].tolist()
    path = write_versioned_table(
        table,
        "release7_early_lesion_identity",
        seed=DEFAULT_SEED,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        notes="metadata-only source-identity audit; no Release-7 H5AD or expression data read",
        extra_meta={
            "prereg": "docs/prereg_release7_early_lesion_feasibility.md",
            "verdict": verdict,
            "n_release7_participants": len(release7_ids),
            "n_cached_chen_participants": len(chen_ids),
            "n_overlapping_participants": len(overlapping),
            "overlapping_participants": overlapping,
            "release7_h5ad_ids": list(R7_H5AD_IDS),
            "release7_synapse_ids": list(R7_SYNAPSE_IDS),
            "file_content_read": False,
            "what_this_licenses": (
                "Only the next metadata-access gate for the fixed Release-7 source."
            ),
            "what_this_does_not_license": (
                "H5AD access, expression analysis, a normal-to-lesion contrast, a "
                "three-stage claim, or an AD-versus-SSL claim."
            ),
            **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
        },
    )
    log.info("wrote %s", path)
    log.info("%s — %d Release-7 IDs; %d cached Chen IDs; %d overlap", verdict,
             len(release7_ids), len(chen_ids), len(overlapping))
    if overlapping:
        log.info("overlapping IDs: %s", ", ".join(overlapping))
        return 5
    return 0


if __name__ == "__main__":
    sys.exit(main())
