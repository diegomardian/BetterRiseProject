"""Outcome-free inventory helpers for the Chen/HTAN MxIF A2 extension.

The public MILWRM deposit is a processed pixel product, not a single-cell
table.  These helpers keep three distinctions explicit before an A2 statistic
is chosen:

* an archive entry is not necessarily described in the paper's sample table;
* an exact participant match is not an exact biospecimen match; and
* multiple regions from one participant do not increase the inferential n.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from src.common.paths import CONFIG_DIR

REGION_MANIFEST = CONFIG_DIR / "a2_mxif_regions.csv"
AVENUE_A_DATASET = "VUMC_HTAN_validation"
AVENUE_A_DISCOVERY_DATASET = "VUMC_HTAN_discovery"
AVENUE_A_PATIENT_PREFIX = "Chen_2021_Cell."

REGION_COLUMNS = (
    "participant_id",
    "batch_name",
    "slide_region",
    "tissue_category",
    "broad_precancer_type",
    "metadata_status",
    "archive_entry",
)

MXIF_CHANNELS = (
    "BCATENIN", "CD20", "CD3D", "CD4_", "CD68", "CD8", "CGA",
    "COLLAGEN", "ERBB2", "FOXP3", "HLAA", "LYSOZYME", "MUC2",
    "NAKATPASE", "OLFM4", "PANCK", "PCNA", "PEGFR", "PSTAT3", "SMA",
    "SOX9", "VIMENTIN", "GACTIN", "CDX2", "MUC5AC", "DAPI",
)

_ENTRY = re.compile(
    r"^Upload_data/colon_polyp_images/"
    r"(?P<batch>HTA11_[^/]+?)_ADJ_(?P<region>region_\d{3})_downsampled\.npz$"
)


class A2MxIFInventoryError(ValueError):
    """The public-product inventory is internally inconsistent."""


def read_region_manifest(path: str | Path = REGION_MANIFEST) -> pd.DataFrame:
    """Read and validate the fixed public-product inventory."""
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = sorted(set(REGION_COLUMNS) - set(frame.columns))
    if missing:
        raise A2MxIFInventoryError(f"region manifest is missing columns {missing}")
    frame = frame.loc[:, REGION_COLUMNS].copy()
    validate_region_manifest(frame)
    return frame


def validate_region_manifest(frame: pd.DataFrame) -> None:
    """Refuse mismatched identifiers and labels that look more certain than they are."""
    if frame.empty:
        raise A2MxIFInventoryError("region manifest is empty")
    for column in ("archive_entry", "slide_region"):
        duplicated = frame.loc[frame[column].duplicated(), column].tolist()
        if duplicated:
            raise A2MxIFInventoryError(f"duplicate {column}: {duplicated[:3]}")

    allowed_status = {"supplementary_table_1", "archive_only"}
    bad_status = sorted(set(frame["metadata_status"]) - allowed_status)
    if bad_status:
        raise A2MxIFInventoryError(f"unknown metadata_status values {bad_status}")

    labelled = frame["metadata_status"].eq("supplementary_table_1")
    if frame.loc[labelled, ["tissue_category", "broad_precancer_type"]].eq("").any().any():
        raise A2MxIFInventoryError("Supplementary Table 1 rows must carry both labels")
    if frame.loc[~labelled, ["tissue_category", "broad_precancer_type"]].ne("").any().any():
        raise A2MxIFInventoryError(
            "archive-only rows may not inherit labels absent from Supplementary Table 1"
        )
    if not set(frame.loc[labelled, "tissue_category"]) <= {"normal", "tumor"}:
        raise A2MxIFInventoryError("tissue_category must be normal or tumor")
    if not set(frame.loc[labelled, "broad_precancer_type"]) <= {"AD", "SSL"}:
        raise A2MxIFInventoryError("broad_precancer_type must be AD or SSL")

    for row in frame.itertuples(index=False):
        match = _ENTRY.fullmatch(row.archive_entry)
        if match is None:
            raise A2MxIFInventoryError(f"unrecognised archive entry {row.archive_entry!r}")
        expected_region = f"{match.group('batch')}_{match.group('region')}"
        if match.group("batch") != row.batch_name or expected_region != row.slide_region:
            raise A2MxIFInventoryError(
                f"archive/metadata identifier mismatch for {row.archive_entry!r}"
            )
        if not row.batch_name.startswith(f"{row.participant_id}_"):
            raise A2MxIFInventoryError(
                f"participant {row.participant_id!r} does not own batch {row.batch_name!r}"
            )


def participant_inventory(regions: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    """One metadata-only row per patient; never one inferential row per region."""
    validate_region_manifest(regions)
    required_obs = {"dataset", "patient_id", "sample_type"}
    missing = sorted(required_obs - set(obs.columns))
    if missing:
        raise A2MxIFInventoryError(f"ICBI obs is missing columns {missing}")

    labelled = regions[regions["metadata_status"].eq("supplementary_table_1")]
    types = (
        labelled.groupby("participant_id")["broad_precancer_type"]
        .agg(lambda x: sorted(set(x)))
    )
    mixed = types[types.map(len).ne(1)]
    if len(mixed):
        raise A2MxIFInventoryError(
            f"participants map to multiple precancer types: {mixed.index.tolist()}"
        )

    per_patient = regions.groupby("participant_id", sort=True).agg(
        n_archive_regions=("archive_entry", "size"),
        n_labeled_regions=("metadata_status", lambda x: int((x == "supplementary_table_1").sum())),
        n_unlabeled_archive_regions=("metadata_status", lambda x: int((x == "archive_only").sum())),
    )
    per_patient["broad_precancer_type"] = types.map(lambda x: x[0])
    batch_counts = labelled.pivot_table(
        index=["participant_id", "batch_name"],
        columns="tissue_category",
        values="slide_region",
        aggfunc="size",
        fill_value=0,
    )
    batch_counts["normal"] = batch_counts.get("normal", 0)
    batch_counts["tumor"] = batch_counts.get("tumor", 0)
    batch_counts["is_pair"] = batch_counts["normal"].gt(0) & batch_counts["tumor"].gt(0)
    patient_counts = batch_counts.groupby(level="participant_id").agg(
        n_normal_regions=("normal", "sum"),
        n_tumor_regions=("tumor", "sum"),
        n_mxif_paired_batches=("is_pair", "sum"),
    )
    per_patient = per_patient.join(patient_counts)
    per_patient["mxif_within_batch_pair"] = per_patient["n_mxif_paired_batches"].gt(0)

    avenue = obs.loc[obs["dataset"].eq(AVENUE_A_DATASET), ["patient_id", "sample_type"]].copy()
    avenue["participant_id"] = avenue["patient_id"].astype(str).str.removeprefix(
        AVENUE_A_PATIENT_PREFIX
    )
    arms = avenue.groupby("participant_id")["sample_type"].agg(set)
    per_patient["avenue_a_patient_match"] = per_patient.index.isin(arms.index)
    per_patient["avenue_a_has_polyp"] = per_patient.index.map(
        lambda patient: "polyp" in arms.get(patient, set())
    )
    per_patient["avenue_a_has_healthy_normal"] = per_patient.index.map(
        lambda patient: "healthy normal" in arms.get(patient, set())
    )
    per_patient["avenue_a_has_both_arms"] = (
        per_patient["avenue_a_has_polyp"] & per_patient["avenue_a_has_healthy_normal"]
    )
    discovery_patients = set(
        obs.loc[obs["dataset"].eq(AVENUE_A_DISCOVERY_DATASET), "patient_id"]
        .astype(str)
        .str.removeprefix(AVENUE_A_PATIENT_PREFIX)
    )
    per_patient["avenue_a_in_discovery"] = per_patient.index.isin(discovery_patients)
    # HTA11_866 is the Chen DIS/VAL shared patient. This metadata flag is not a
    # split assignment and must never be used to promote A2's participant
    # overlap into an independent replication.
    per_patient["avenue_a_disval_shared"] = (
        per_patient["avenue_a_patient_match"] & per_patient["avenue_a_in_discovery"]
    )
    per_patient["crosswalk_level"] = "participant_only"
    per_patient["specimen_exact"] = False
    per_patient["conventional_ad_mxif_pair"] = (
        per_patient["broad_precancer_type"].eq("AD")
        & per_patient["mxif_within_batch_pair"]
    )
    return per_patient.reset_index()


def inventory_counts(regions: pd.DataFrame, participants: pd.DataFrame) -> dict[str, int]:
    """The Step-1 cardinalities, kept separate from any biological result."""
    return {
        "archive_regions": int(len(regions)),
        "labeled_regions": int(regions["metadata_status"].eq("supplementary_table_1").sum()),
        "unlabeled_archive_regions": int(regions["metadata_status"].eq("archive_only").sum()),
        "participants": int(len(participants)),
        "mxif_within_batch_pairs": int(participants["mxif_within_batch_pair"].sum()),
        "conventional_ad_mxif_pairs": int(participants["conventional_ad_mxif_pair"].sum()),
        "avenue_a_patient_matches": int(participants["avenue_a_patient_match"].sum()),
        "avenue_a_both_arm_patients": int(participants["avenue_a_has_both_arms"].sum()),
        "avenue_a_disval_shared_patients": int(participants["avenue_a_disval_shared"].sum()),
    }
