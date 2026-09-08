"""Outcome-blind Wnt detection gate for Becker's paired polyp donors.

Run on the cluster where the checksum-pinned GSE201348 files live:

    python -m src.reference.jobs.becker_lesion_wnt_gate \
      --tar "$BRP_DATA_DIR/raw/becker/GSE201348_RAW.tar" \
      --series-matrix "$BRP_DATA_DIR/raw/becker/GSE201348_series_matrix.txt.gz"

This is deliberately not an analysis. It reads only the five invariant-8 Wnt
targets in paired donors' polyp nuclei and stops after deciding whether a later
lesion-within-donor specification has a measurable predictor.

Exit status 0 means the gate passed. Exit status 5 means a durable no-substrate
artifact was written. Other failures are infrastructure or unrecognised-input
errors and do not claim a biological or assay result.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.label_provenance import Measurement, check_no_circular_claim, provenance_meta
from src.common.provenance import DEFAULT_SEED
from src.reference.jobs.coexpression_silencing import DETECTION_MIN_UMI
from src.reference.wnt_score import MIN_CELLS_FOR_CORRELATION, SIGNATURE

log = logging.getLogger(__name__)

EXPECTED_PAIRED_DONORS = ("A001", "A002", "A014", "A015")
MIN_DETECTION_RATE = 0.01
MIN_DETECTED_NUCLEI = MIN_CELLS_FOR_CORRELATION

LABEL_PROVENANCE = Measurement(
    modality="sample_annotation",
    assay="Becker GEO disease-stage metadata",
    genes=(),
)
CLAIM_PROVENANCE = Measurement(
    modality="transcript",
    assay="Becker GSE201348 snRNA-seq",
    genes=SIGNATURE,
)


class WntDetectionGateError(ValueError):
    """The fixed input cannot support the detection gate."""

    def __init__(self, message: str, *, failure_kind: str = "unrecognised_input"):
        super().__init__(message)
        self.failure_kind = failure_kind


def validate_specification() -> tuple[str, ...]:
    """Run invariant 11 before checking a path or reading a matrix."""
    return check_no_circular_claim(
        labels=LABEL_PROVENANCE, claim=CLAIM_PROVENANCE
    )


def signature_index(symbols: Sequence[str]) -> dict[str, int]:
    """Locate every frozen Wnt target exactly once in a feature-symbol vector."""
    values = np.asarray([str(value) for value in symbols])
    index: dict[str, int] = {}
    missing: list[str] = []
    duplicated: list[str] = []
    for gene in SIGNATURE:
        hits = np.flatnonzero(values == gene)
        if hits.size == 0:
            missing.append(gene)
        elif hits.size > 1:
            duplicated.append(gene)
        else:
            index[gene] = int(hits[0])
    if missing or duplicated:
        raise WntDetectionGateError(
            "Wnt signature lookup refused: "
            + (f"missing symbols {missing}. " if missing else "")
            + (f"duplicated symbols {duplicated}." if duplicated else ""),
            failure_kind="identifier_assay",
        )
    return index


def _column(values, column: int) -> np.ndarray:
    selected = values[:, column]
    return (np.asarray(selected.toarray()).ravel()
            if hasattr(selected, "toarray") else np.asarray(selected).ravel())


def summarise_blocks(blocks: Sequence[dict[str, object]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarise raw polyp-nucleus detection by donor and frozen signature gene.

    ``blocks`` is deliberately a simple boundary for synthetic forcing tests and
    the tar reader. Every block must name its donor, physical lesion, sequencing
    row, counts matrix, and feature symbols.
    """
    if not blocks:
        raise WntDetectionGateError("no paired-donor polyp blocks were supplied")

    records: list[dict[str, object]] = []
    reference_symbols: np.ndarray | None = None
    for block in blocks:
        donor = str(block["donor"])
        lesion = str(block["sample_id"])
        gsm = str(block["gsm"])
        counts = block["counts"]
        symbols = np.asarray([str(value) for value in block["symbols"]])
        if reference_symbols is None:
            reference_symbols = symbols
        elif not np.array_equal(symbols, reference_symbols):
            raise WntDetectionGateError(
                f"{gsm} has a different feature order from the first included "
                "sample; columns cannot be combined safely.",
                failure_kind="identifier_assay",
            )
        index = signature_index(symbols)
        for gene, column in index.items():
            values = _column(counts, column)
            records.append({
                "donor": donor,
                "sample_id": lesion,
                "gsm": gsm,
                "gene": gene,
                "n_nuclei": int(values.size),
                "n_detected": int((values >= DETECTION_MIN_UMI).sum()),
            })

    raw = pd.DataFrame(records)
    grouped = raw.groupby(["donor", "gene"], sort=True, observed=True)
    rows = []
    for (donor, gene), group in grouped:
        n_nuclei = int(group["n_nuclei"].sum())
        n_detected = int(group["n_detected"].sum())
        detection = n_detected / n_nuclei if n_nuclei else float("nan")
        rows.append({
            "donor": donor,
            "gene": gene,
            "n_lesions": int(group["sample_id"].nunique()),
            "n_sequence_rows": int(group["gsm"].nunique()),
            "n_nuclei": n_nuclei,
            "n_detected": n_detected,
            "detection": detection,
            "passes_detected_nuclei": n_detected >= MIN_DETECTED_NUCLEI,
            "passes_detection_rate": detection >= MIN_DETECTION_RATE,
        })
    by_donor = pd.DataFrame(rows).sort_values(["donor", "gene"], ignore_index=True)
    by_donor["passes"] = (
        by_donor["passes_detected_nuclei"] & by_donor["passes_detection_rate"]
    )

    observed = tuple(sorted(by_donor["donor"].unique()))
    if observed != EXPECTED_PAIRED_DONORS:
        raise WntDetectionGateError(
            f"paired donors {observed}, expected {EXPECTED_PAIRED_DONORS}; "
            "do not silently run a smaller cohort."
        )
    expected_rows = len(EXPECTED_PAIRED_DONORS) * len(SIGNATURE)
    if len(by_donor) != expected_rows:
        raise WntDetectionGateError(
            f"{len(by_donor)} donor-by-gene rows, expected {expected_rows}; "
            "a missing row is not a failed measurement."
        )

    summary = by_donor.groupby("gene", sort=True, observed=True).agg(
        n_paired_donors=("donor", "nunique"),
        nuclei_total=("n_nuclei", "sum"),
        detected_total=("n_detected", "sum"),
        donors_passing=("passes", "sum"),
    ).reset_index()
    summary["detection_pooled"] = summary["detected_total"] / summary["nuclei_total"]
    summary["passes_all_paired_donors"] = (
        summary["donors_passing"] == len(EXPECTED_PAIRED_DONORS)
    )
    return by_donor, summary


def gate_verdict(by_donor: pd.DataFrame) -> dict[str, object]:
    """The fixed consequence; no Wnt or differentiation outcome is consulted."""
    failing = by_donor.loc[~by_donor["passes"], ["donor", "gene"]]
    if failing.empty:
        return {
            "verdict": "PASSES DETECTION GATE",
            "detail": (
                "every frozen Wnt target clears both detection floors in every "
                "paired donor. This licenses a separate lesion-level Wnt "
                "pre-specification, not a fitted analysis."
            ),
        }
    names = ", ".join(f"{row.donor}/{row.gene}" for row in failing.itertuples())
    return {
        "verdict": "NO SUBSTRATE AT THIS snRNA RESOLUTION",
        "detail": (
            f"{names} fail one or both fixed detection floors. This is a "
            "measurement limit, not a Wnt null; do not fit a reduced signature."
        ),
    }


def read_paired_polyp_blocks(tar: Path, series_matrix: Path) -> list[dict[str, object]]:
    """Read only paired-donor polyp triplets, retaining replicate rows visibly."""
    from src.reference.becker_io import (
        gene_symbols,
        paired_donors,
        read_series_matrix,
        read_triplet,
        sample_files,
    )

    metadata = read_series_matrix(series_matrix)
    paired = tuple(sorted(paired_donors(metadata)))
    if paired != EXPECTED_PAIRED_DONORS:
        raise WntDetectionGateError(
            f"metadata paired donors {paired}, expected {EXPECTED_PAIRED_DONORS}"
        )
    files = sample_files(tar).merge(
        metadata, on="gsm", how="inner", suffixes=("_tar", "")
    )
    selected = files[(files["arm"] == "tumour") & files["donor"].isin(paired)]
    if selected.empty:
        raise WntDetectionGateError("no paired-donor polyp triplets after metadata join")
    mismatched = selected.loc[
        selected["sample_id_tar"] != selected["sample_id"]
    ]
    if not mismatched.empty:
        names = ", ".join(
            f"{row.gsm}: tar={row.sample_id_tar}, metadata={row.sample_id}"
            for row in mismatched.itertuples()
        )
        raise WntDetectionGateError(
            f"tar and metadata sample identifiers disagree: {names}",
        )
    incomplete = selected.loc[~selected["complete"]]
    if not incomplete.empty:
        affected = ", ".join(
            f"{row.donor}/{row.gsm}/{row.sample_id}"
            for row in incomplete.sort_values(["donor", "sample_id", "gsm"]).itertuples()
        )
        raise WntDetectionGateError(
            "required paired-donor polyp triplet is incomplete: "
            f"{affected}. No sample was silently dropped.",
            failure_kind="incomplete_triplet",
        )

    blocks: list[dict[str, object]] = []
    for _, row in selected.sort_values(["donor", "sample_id", "gsm"]).iterrows():
        counts, _, features = read_triplet(tar, row)
        blocks.append({
            "donor": row["donor"],
            "sample_id": row["sample_id"],
            "gsm": row["gsm"],
            "counts": counts,
            "symbols": gene_symbols(features),
        })
    return blocks


def failure_outcome(error: WntDetectionGateError) -> dict[str, str]:
    """Classify the fixed no-substrate failures that must leave an artifact."""
    if error.failure_kind == "identifier_assay":
        verdict = "NO SUBSTRATE — identifier/assay failure"
    elif error.failure_kind == "incomplete_triplet":
        verdict = "NO SUBSTRATE — incomplete assay triplet"
    else:
        raise error
    return {"verdict": verdict, "detail": str(error), "failure_kind": error.failure_kind}


def result_meta(outcome: dict[str, object]) -> dict[str, object]:
    """Metadata common to passing, detection-failure, and assay-failure records."""
    return {
        "prereg": "docs/prereg_becker_lesion_wnt_detection_gate.md",
        "input": "GSE201348 paired-donor polyp nuclei only",
        "paired_donors": list(EXPECTED_PAIRED_DONORS),
        "signature": list(SIGNATURE),
        "detection_min_umi": DETECTION_MIN_UMI,
        "min_detected_nuclei_per_donor": MIN_DETECTED_NUCLEI,
        "min_detection_rate_per_donor": MIN_DETECTION_RATE,
        "unit": "donor; lesion and sequencing-row counts are descriptive only",
        "verdict": outcome,
        "exploratory": False,
        "pre_registered": True,
        **provenance_meta(LABEL_PROVENANCE, CLAIM_PROVENANCE),
    }


def write_failure_artifact(
    outcome: dict[str, str], *, seed: int, results_dir: Path | None, allow_dirty: bool
) -> Path:
    """Persist a fixed assay/identifier no-substrate result before exit status 5."""
    frame = pd.DataFrame([outcome])
    return write_versioned_table(
        frame, "becker_lesion_wnt_detection_failure", seed=seed,
        results_dir=results_dir, allow_dirty=allow_dirty, extra_meta=result_meta(outcome),
    )


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tar", required=True, type=Path, help="GSE201348_RAW.tar")
    parser.add_argument("--series-matrix", required=True, type=Path,
                        help="GSE201348_series_matrix.txt.gz")
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    validate_specification()
    for path in (args.tar, args.series_matrix):
        if not path.exists():
            raise SystemExit(f"{path} not found")

    try:
        by_donor, summary = summarise_blocks(
            read_paired_polyp_blocks(args.tar, args.series_matrix)
        )
    except WntDetectionGateError as error:
        outcome = failure_outcome(error)
        log.info("wrote %s", write_failure_artifact(
            outcome, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty,
        ))
        log.info("%s\n%s", outcome["verdict"], outcome["detail"])
        return 5

    outcome = gate_verdict(by_donor)
    meta = result_meta(outcome)
    for frame, name in (
        (by_donor, "becker_lesion_wnt_detection_by_donor"),
        (summary, "becker_lesion_wnt_detection_summary"),
    ):
        log.info("wrote %s", write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    log.info("%s\n%s", outcome["verdict"], outcome["detail"])
    return 0 if outcome["verdict"] == "PASSES DETECTION GATE" else 5


if __name__ == "__main__":
    sys.exit(main())
