"""GSE178341 ingest and the week-1 verification guards. W1.

Four things have to be true before a single cell is analysed, and each of them
is cheaper to check in week 1 than to discover in week 4:

1. **The counts are raw.** Normalised or log-transformed values look fine in a
   histogram and quietly invalidate every downstream count-based method.
   :func:`assert_raw_counts`.
2. **The droplets are unfiltered.** SoupX and CellBender both learn the ambient
   profile from *empty* droplets. A GEO deposit that carries only the filtered,
   cell-containing matrix cannot support either — which takes week 2 and gate
   criterion G1 with it. :func:`assert_unfiltered_droplets`.
3. **The matched normals are there.** If normals are sparse, the compositional
   term loses its reference and the project changes shape
   (execution_plan.md §8.2). :func:`matched_normal_report`.
4. **Every file is recorded.** A download nobody put in the manifest is a result
   nobody can reproduce (data/README.md). :func:`append_manifest_row`.

The guards are deliberately plain functions over a count matrix rather than over
AnnData, so they are testable without the single-cell stack and reusable on the
Lee cohorts by W4. AnnData is imported lazily, inside the loaders.

Matrix orientation is **cells x genes** throughout — the AnnData convention.
10x `.mtx` files ship genes x barcodes, so transpose on read. Orientation cannot
be detected reliably (a raw droplet matrix has far more barcodes than genes, a
filtered one far fewer), so it is required rather than guessed.
"""

from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np

from src.common.paths import MANIFEST_PATH

#: Barcodes below this many UMIs are treated as empty droplets. 100 is the
#: conventional soup-profile floor: high enough to exclude real cells, low
#: enough that the ambient profile is estimated from a large population.
EMPTY_DROPLET_UMI_THRESHOLD = 100

#: A raw 10x matrix is overwhelmingly empty droplets — typically >95% of the
#: ~737k-barcode whitelist. Anything below half is not a raw droplet matrix.
MIN_EMPTY_DROPLET_FRACTION = 0.5

#: Row sums this uniform mean the matrix was normalised to a fixed library size
#: (CPM, CP10K) and then possibly rounded back to integers.
MAX_LIBRARY_SIZE_CV = 0.01

MANIFEST_FIELDS = (
    "path", "sha256", "bytes", "source_url", "accession",
    "downloaded_on", "downloaded_by", "workstream", "notes",
)


class IngestError(AssertionError):
    """A dataset failed a week-1 precondition. Stop; do not work around it."""


# ---------------------------------------------------------------------------
# Matrix helpers — work on numpy arrays and scipy sparse alike
# ---------------------------------------------------------------------------


def _row_sums(matrix: Any) -> np.ndarray:
    """Per-barcode total counts, for dense or sparse input."""
    return np.asarray(matrix.sum(axis=1)).ravel()


def _values(matrix: Any) -> np.ndarray:
    """The stored values. For sparse input the implicit zeros are irrelevant to
    every check here — they are non-negative, finite and integral by definition."""
    data = getattr(matrix, "data", None)
    return np.asarray(matrix if data is None else data).ravel()


# ---------------------------------------------------------------------------
# 1. Raw counts
# ---------------------------------------------------------------------------


def assert_raw_counts(matrix: Any, *, context: str = "the count matrix") -> None:
    """Raise unless `matrix` holds raw integer counts.

    Catches the three things that arrive in a file labelled "counts":
    log1p-transformed values (non-integral), CPM/CP10K values (non-integral),
    and normalised-then-rounded values (integral, but with near-constant library
    sizes). GSE178341's supplementary structure is awkward enough that this is a
    real risk, not a formality.
    """
    values = _values(matrix)
    if values.size == 0:
        raise IngestError(f"{context} is empty.")

    if not np.all(np.isfinite(values)):
        raise IngestError(f"{context} contains NaN or infinite values.")

    if values.min() < 0:
        raise IngestError(
            f"{context} contains negative values (min {values.min():.4g}). "
            f"That is scaled or centred data, not counts."
        )

    if not np.allclose(values, np.rint(values)):
        offenders = values[~np.isclose(values, np.rint(values))][:3]
        raise IngestError(
            f"{context} contains non-integer values (e.g. {offenders}). "
            f"This is normalised or log-transformed data, not raw counts. Every "
            f"count-based step downstream — ambient correction, doublet "
            f"detection, inferCNV — assumes raw."
        )

    totals = _row_sums(matrix)
    if totals.size > 1 and totals.mean() > 0:
        cv = float(totals.std() / totals.mean())
        if cv < MAX_LIBRARY_SIZE_CV:
            raise IngestError(
                f"{context} has integer values but near-constant library sizes "
                f"(CV {cv:.2e}, mean {totals.mean():.0f}). That is normalised "
                f"data rounded back to integers, which is not the same thing as "
                f"raw counts."
            )


# ---------------------------------------------------------------------------
# 2. Unfiltered droplets — the check nothing else in the repo makes
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DropletProfile:
    """What the barcode-count distribution says about how a matrix was filtered."""

    n_barcodes: int
    n_empty: int
    empty_fraction: float
    median_umi: float
    max_umi: float

    @property
    def looks_unfiltered(self) -> bool:
        return self.empty_fraction >= MIN_EMPTY_DROPLET_FRACTION

    def summary(self) -> str:
        return (
            f"{self.n_barcodes:,} barcodes, {self.n_empty:,} below "
            f"{EMPTY_DROPLET_UMI_THRESHOLD} UMI ({self.empty_fraction:.1%}), "
            f"median {self.median_umi:,.0f} UMI, max {self.max_umi:,.0f}"
        )


def droplet_profile(
    matrix: Any, *, empty_threshold: int = EMPTY_DROPLET_UMI_THRESHOLD
) -> DropletProfile:
    """Describe the barcode-count distribution. Cheap, and worth logging always."""
    totals = _row_sums(matrix)
    n_empty = int((totals < empty_threshold).sum())
    return DropletProfile(
        n_barcodes=int(totals.size),
        n_empty=n_empty,
        empty_fraction=float(n_empty / totals.size) if totals.size else 0.0,
        median_umi=float(np.median(totals)) if totals.size else 0.0,
        max_umi=float(totals.max()) if totals.size else 0.0,
    )


def assert_unfiltered_droplets(
    matrix: Any,
    *,
    context: str = "the count matrix",
    empty_threshold: int = EMPTY_DROPLET_UMI_THRESHOLD,
    min_empty_fraction: float = MIN_EMPTY_DROPLET_FRACTION,
) -> DropletProfile:
    """Raise unless `matrix` still contains empty droplets. Returns the profile.

    This is the week-1 blocking check for the ambient-correction arm. Both SoupX
    and CellBender estimate the soup from barcodes that contain no cell; a
    filtered matrix has had exactly those removed. "Raw counts" and "raw
    droplets" are different properties and a GEO deposit can satisfy the first
    while failing the second.

    If this raises, the fallbacks in order are: find the original submission's
    unfiltered HDF5s, ask the authors, or move G1's evidence onto the
    plate-based ambient-free subset (execution_plan.md §8.2). Report the answer
    at the week-1 meeting either way — a G1 that cannot be run is a week-1
    finding, not a week-5 surprise.
    """
    profile = droplet_profile(matrix, empty_threshold=empty_threshold)
    if profile.empty_fraction < min_empty_fraction:
        raise IngestError(
            f"{context} looks cell-filtered, not raw: {profile.summary()}. "
            f"SoupX and CellBender both learn the ambient profile from empty "
            f"droplets, so neither can run on this, and gate criterion G1 loses "
            f"its foundation. Find the unfiltered matrix before week 2 — see "
            f"src/reference/ingest.py:assert_unfiltered_droplets for the fallbacks."
        )
    return profile


# ---------------------------------------------------------------------------
# 3. Cohort structure — the week-1 deliverables
# ---------------------------------------------------------------------------


def cells_by_patient_and_tissue(patient_id: Any, tissue: Any) -> Any:
    """Cell counts by patient x tissue. The week-1 "done when" (§4, W1 row 1)."""
    import pandas as pd

    df = pd.DataFrame({"patient_id": list(patient_id), "tissue": list(tissue)})
    table = (
        df.groupby(["patient_id", "tissue"], observed=True)
        .size()
        .unstack("tissue", fill_value=0)
        .sort_index()
    )
    table.columns.name = None
    return table


def matched_normal_report(
    patient_id: Any,
    tissue: Any,
    *,
    normal_label: str = "normal",
    tumour_label: str = "tumour",
    min_cells: int = 1,
) -> Any:
    """Per-patient tumour/normal availability, with a `matched` flag.

    execution_plan.md §8.2 lists matched-normal completeness as a week-1
    check that can invalidate the plan: the compositional term is
    Delta(mature fraction) against the patient's own normal, so a patient
    without one contributes to no compositional estimate. Run this in week 1 and
    report the count, not the assumption.
    """
    table = cells_by_patient_and_tissue(patient_id, tissue)
    for column in (normal_label, tumour_label):
        if column not in table.columns:
            table[column] = 0
    report = table.loc[:, [tumour_label, normal_label]].copy()
    report["matched"] = (report[tumour_label] >= min_cells) & (
        report[normal_label] >= min_cells
    )
    return report


# ---------------------------------------------------------------------------
# 4. Manifest
# ---------------------------------------------------------------------------


def sha256_file(path: str | Path, *, chunk_bytes: int = 1 << 20) -> str:
    """Streaming sha256, so this works on a 12 GB matrix."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(chunk_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def append_manifest_row(
    path: str | Path,
    *,
    source_url: str,
    accession: str,
    downloaded_by: str,
    workstream: str = "W1",
    notes: str = "",
    manifest_path: Path | None = None,
) -> dict[str, str]:
    """Record a downloaded file in data/manifest.csv. Idempotent on path.

    Call this in the same PR as the code that reads the file (CONTRIBUTING §4).
    """
    manifest = Path(manifest_path or MANIFEST_PATH)
    target = Path(path)
    row = {
        "path": str(target),
        "sha256": sha256_file(target),
        "bytes": str(target.stat().st_size),
        "source_url": source_url,
        "accession": accession,
        "downloaded_on": date.today().isoformat(),
        "downloaded_by": downloaded_by,
        "workstream": workstream,
        "notes": notes,
    }

    existing: list[dict[str, str]] = []
    if manifest.exists() and manifest.stat().st_size:
        with open(manifest, newline="", encoding="utf-8") as handle:
            existing = [r for r in csv.DictReader(handle) if r["path"] != row["path"]]

    manifest.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(MANIFEST_FIELDS))
        writer.writeheader()
        writer.writerows([*existing, row])
    return row


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def verify_ingest(
    matrix: Any,
    *,
    context: str = "the count matrix",
    require_unfiltered: bool = True,
) -> DropletProfile:
    """Run every matrix-level week-1 guard. Call this on anything freshly read.

    `require_unfiltered=False` is for the cell-filtered object you work from
    *after* ambient correction — never for the input to SoupX or CellBender.
    """
    assert_raw_counts(matrix, context=context)
    profile = droplet_profile(matrix)
    if require_unfiltered:
        assert_unfiltered_droplets(matrix, context=context)
    return profile


def read_10x_mtx(directory: str | Path, *, prefix: str = "", verify: bool = True) -> Any:
    """Read a 10x mtx triplet into AnnData (cells x genes) and verify it.

    Kept thin on purpose: GSE178341's supplementary layout is awkward and the
    exact filenames are not known until the download lands, so the parsing that
    is genuinely dataset-specific belongs in the ingest script, not here. What
    this guarantees is that nothing skips the guards.
    """
    import scanpy as sc

    adata = sc.read_10x_mtx(Path(directory), prefix=prefix or None, cache=False)
    if verify:
        verify_ingest(adata.X, context=str(directory))
    return adata
