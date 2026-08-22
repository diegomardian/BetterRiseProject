"""Read the ICBI CRC atlas metadata WITHOUT downloading the atlas. W1, week 1.

execution_plan.md §8.2, boxed: *"Pull the ICBI atlas metadata table only — not
the 4.27M-cell object. That single table gives you the real sample size before
you commit any compute."*

The problem is that no metadata-only table is published. https://crc.icbi.at/
offers ``final_crc_atlas-adata.h5ad`` and that is **32.7 GB**. The Zenodo record
holds only multi-gigabyte archives, and the repo's ``tables/`` directory has a
schema but no sample sheet.

An ``.h5ad`` is HDF5, though, and HDF5 is random-access: ``/obs`` lives in its own
group and can be read on its own. The server sends ``Accept-Ranges: bytes``, so
this fetches the byte ranges HDF5 asks for and nothing else. Reading a dozen obs
columns costs a few hundred MB against 32.7 GB — the spirit of "metadata table
only", implemented.

What the metadata answers
-------------------------
Four live questions at once, from ``tables/reference_meta.yaml``'s vocabulary:

``platform``
    ``10x 5' v2``, ``10x 3' v2``, ``10x 3' v3``, **``smartseq2``**, ``bd rhapsody``.
    The plate-based subset is the one that matters twice over: plate protocols
    have essentially no ambient soup, so an intrinsic signal surviving there is
    strong evidence it is not contamination (§8.2) — that is G1's fallback after
    open decision #8. And they sequence far deeper per cell, which is the only
    route to a maturity call on markers as sparsely detected as axis 1's, where
    GSE178341 leaves 49% of epithelium tied (open decision #14).
``sample_type``
    ``tumor`` / ``normal`` / ``polyp`` / ``metastasis`` — paired availability, the
    real sample size, against the 36 of 62 this project has on GSE178341 (#9).
``enrichment_cell_types``
    ``naive`` / ``CD45+`` / ``CD45-`` — the same sorting problem as decision #11,
    at atlas scale.
``matrix_type``
    ``raw counts`` / ``processed counts`` / ``log norm`` — which studies could in
    principle support ambient correction at all.
"""

from __future__ import annotations

import io
from typing import Any, Final

import numpy as np
import pandas as pd

ATLAS_URL: Final[str] = "https://crc.icbi.at/h5ad/final_crc_atlas-adata.h5ad"

#: Obs columns worth pulling. Every one maps to a live decision; adding more is
#: cheap but each is ~17 MB of range reads at 4.27M cells, so it is not free.
DEFAULT_COLUMNS: Final[tuple[str, ...]] = (
    # identity
    "study_id", "dataset", "sample_id", "patient_id",
    # design
    "sample_type", "platform", "platform_fine", "enrichment_cell_types",
    "matrix_type", "reference_genome", "suspension_type", "tissue_cell_state",
    # cell annotation
    "atlas_cell_type_coarse",
    # PER-CELL DEPTH. The point of the exercise: n_genes and total_counts let us
    # measure whether a plate-based subset actually resolves the sparse axis-1
    # markers, rather than assuming it from the protocol's reputation.
    "n_genes", "total_counts", "pct_counts_mito",
    # already-computed doublet calls, and the tier-B stratifier
    "SOLO_doublet_status", "MLH1_promoter_methylation_status",
    "microsatellite_status",
)

#: Platforms with essentially no ambient soup and much deeper per-cell coverage.
PLATE_PLATFORMS: Final[tuple[str, ...]] = ("smartseq2",)

#: 4 MB blocks. HDF5 issues many small reads; without block caching this would be
#: thousands of HTTP requests instead of a few hundred.
BLOCK_SIZE: Final[int] = 4 * 1024 * 1024


class ICBIError(RuntimeError):
    """The atlas metadata could not be read."""


class HTTPRangeFile(io.RawIOBase):
    """A seekable read-only file over HTTP, backed by Range requests.

    h5py accepts any object with ``read``/``seek``/``tell`` via
    ``driver="fileobj"``. Blocks are cached because HDF5 walks B-trees with many
    small reads, and one request per read would be unusable.
    """

    def __init__(self, url: str, *, block_size: int = BLOCK_SIZE, session: Any = None):
        import requests

        self.url = url
        self.block_size = block_size
        self._session = session or requests.Session()
        self._pos = 0
        self._cache: dict[int, bytes] = {}
        self.bytes_fetched = 0

        head = self._session.head(url, allow_redirects=True, timeout=60)
        head.raise_for_status()
        if head.headers.get("Accept-Ranges") != "bytes":
            raise ICBIError(
                f"{url} does not advertise byte ranges; the whole object would "
                f"have to be downloaded."
            )
        self.size = int(head.headers["Content-Length"])

    # -- io plumbing --------------------------------------------------------
    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self._pos

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        base = {io.SEEK_SET: 0, io.SEEK_CUR: self._pos, io.SEEK_END: self.size}[whence]
        self._pos = max(0, min(self.size, base + offset))
        return self._pos

    def _block(self, index: int) -> bytes:
        if index not in self._cache:
            start = index * self.block_size
            stop = min(start + self.block_size, self.size) - 1
            response = self._session.get(
                self.url, headers={"Range": f"bytes={start}-{stop}"}, timeout=120
            )
            response.raise_for_status()
            self._cache[index] = response.content
            self.bytes_fetched += len(response.content)
        return self._cache[index]

    def readinto(self, buffer) -> int:  # noqa: ANN001 - io protocol
        data = self.read(len(buffer))
        buffer[: len(data)] = data
        return len(data)

    def read(self, size: int = -1) -> bytes:
        if size < 0:
            size = self.size - self._pos
        size = min(size, self.size - self._pos)
        out = bytearray()
        while len(out) < size:
            index = (self._pos + len(out)) // self.block_size
            block = self._block(index)
            offset = (self._pos + len(out)) - index * self.block_size
            out += block[offset : offset + (size - len(out))]
        self._pos += len(out)
        return bytes(out)


def _decode(group: Any, name: str) -> np.ndarray:
    """Read one obs column, resolving anndata's categorical encoding."""
    import h5py

    node = group[name]
    if isinstance(node, h5py.Group):          # categorical: codes + categories
        codes = node["codes"][:]
        categories = node["categories"][:]
        categories = np.asarray(
            [c.decode() if isinstance(c, bytes) else c for c in categories], dtype=object
        )
        out = np.full(codes.shape, None, dtype=object)
        valid = codes >= 0
        out[valid] = categories[codes[valid]]
        return out
    values = node[:]
    if values.dtype.kind in "SO":
        return np.asarray(
            [v.decode() if isinstance(v, bytes) else v for v in values], dtype=object
        )
    return values


def read_atlas_obs(
    url: str = ATLAS_URL,
    columns: Any = DEFAULT_COLUMNS,
    *,
    block_size: int = BLOCK_SIZE,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Read the atlas's per-cell metadata over HTTP. Returns ``(obs, report)``.

    `report` carries the bytes actually fetched, so the claim that this reads the
    metadata rather than the object is checkable rather than asserted.
    """
    import h5py

    handle = HTTPRangeFile(url, block_size=block_size)
    with h5py.File(handle, mode="r", driver="fileobj") as h5:
        if "obs" not in h5:
            raise ICBIError("no /obs group — is this an .h5ad?")
        obs_group = h5["obs"]
        available = list(obs_group.keys())
        wanted = [c for c in columns if c in available] if columns else available
        missing = [c for c in (columns or []) if c not in available]
        data = {name: _decode(obs_group, name) for name in wanted}

    obs = pd.DataFrame(data)
    report = {
        "n_cells": int(len(obs)),
        "columns_read": wanted,
        "columns_missing": missing,
        "columns_available": available,
        "bytes_fetched": handle.bytes_fetched,
        "object_size": handle.size,
        "fraction_fetched": handle.bytes_fetched / handle.size if handle.size else 0.0,
    }
    return obs, report


# ---------------------------------------------------------------------------
# The four questions §8.2 wants answered
# ---------------------------------------------------------------------------


def platform_summary(obs: pd.DataFrame) -> pd.DataFrame:
    """Cells, samples, patients and studies per sequencing platform.

    The row to read is ``smartseq2``: plate protocols have essentially no soup,
    so they are G1's fallback after decision #8, and they sequence deep enough to
    resolve the sparse axis-1 markers that leave 49% of GSE178341's epithelium
    tied (decision #14).
    """
    return _group_summary(obs, "platform").assign(
        plate_based=lambda d: d["platform"].isin(PLATE_PLATFORMS)
    )


def paired_sample_summary(obs: pd.DataFrame) -> pd.DataFrame:
    """Patients with both tumour and normal, per study. The real sample size.

    Compare against the 36 of 62 this project has on GSE178341 (decision #9).
    """
    if not {"patient_id", "sample_type"} <= set(obs.columns):
        raise ICBIError("obs needs patient_id and sample_type")
    frame = obs.loc[:, [c for c in ("study_id", "patient_id", "sample_type") if c in obs]]
    has = (
        frame.assign(one=1)
        .pivot_table(
            index=[c for c in ("study_id", "patient_id") if c in frame],
            columns="sample_type", values="one", aggfunc="size", fill_value=0,
        )
        .reset_index()
    )
    for column in ("tumor", "normal"):
        if column not in has.columns:
            has[column] = 0
    has["paired"] = (has["tumor"] > 0) & (has["normal"] > 0)
    keys = [c for c in ("study_id",) if c in has]
    if not keys:
        return has
    return (
        has.groupby(keys, observed=True)
        .agg(n_patients=("paired", "size"), n_paired=("paired", "sum"))
        .reset_index()
        .sort_values("n_paired", ascending=False, ignore_index=True)
    )


def epithelial_fraction(obs: pd.DataFrame) -> pd.DataFrame:
    """Epithelial share per study. Thin epithelium means a thin mature arm."""
    column = next(
        (c for c in ("atlas_cell_type_coarse", "cell_type_coarse", "cell_type_coarse_study")
         if c in obs.columns),
        None,
    )
    if column is None:
        raise ICBIError("obs has no coarse cell-type column")
    frame = obs.assign(
        is_epi=obs[column].astype(str).str.lower().str.contains("epitheli")
    )
    keys = [c for c in ("study_id", "platform") if c in frame.columns]
    return (
        frame.groupby(keys, observed=True)
        .agg(n_cells=("is_epi", "size"), n_epithelial=("is_epi", "sum"))
        .assign(epithelial_fraction=lambda d: d["n_epithelial"] / d["n_cells"])
        .reset_index()
        .sort_values("n_epithelial", ascending=False, ignore_index=True)
    )


def enrichment_summary(obs: pd.DataFrame) -> pd.DataFrame:
    """Cells by enrichment. ``naive`` is unsorted — decision #11 at atlas scale."""
    return _group_summary(obs, "enrichment_cell_types")


def _group_summary(obs: pd.DataFrame, column: str) -> pd.DataFrame:
    if column not in obs.columns:
        raise ICBIError(f"obs has no {column} column")
    aggregations: dict[str, Any] = {"n_cells": (column, "size")}
    for name, source in (
        ("n_samples", "sample_id"),
        ("n_patients", "patient_id"),
        ("n_studies", "study_id"),
    ):
        if source in obs.columns:
            aggregations[name] = (source, "nunique")
    return (
        obs.groupby(column, observed=True)
        .agg(**aggregations)
        .reset_index()
        .sort_values("n_cells", ascending=False, ignore_index=True)
    )


def depth_by_platform(obs: pd.DataFrame) -> pd.DataFrame:
    """Per-cell sequencing depth by platform. **The number that decides #14.**

    Axis 1 rests on five sparsely detected markers, and on GSE178341 that leaves
    49% of epithelium sharing one score — no gradient, so no graded maturity
    call. Depth is the lever: a cell with more detected genes is less likely to
    drop all five to zero.

    Plate protocols are *reputed* to sequence far deeper than droplet ones. This
    measures it on the actual atlas rather than taking the reputation on trust,
    which matters because the answer decides whether pulling axis 3 forward from
    week 13+ is necessary or merely tidy.

    Reports the median and quartiles of ``n_genes`` and ``total_counts``, plus
    the share of cells clearing 5,000 detected genes — roughly where five sparse
    markers stop dropping out together.
    """
    needed = {"platform", "n_genes"}
    missing = needed - set(obs.columns)
    if missing:
        raise ICBIError(f"obs is missing {sorted(missing)}")

    frame = obs.copy()
    frame["n_genes"] = pd.to_numeric(frame["n_genes"], errors="coerce")
    if "total_counts" in frame.columns:
        frame["total_counts"] = pd.to_numeric(frame["total_counts"], errors="coerce")

    aggregations: dict[str, Any] = {
        "n_cells": ("n_genes", "size"),
        "median_genes": ("n_genes", "median"),
        "q25_genes": ("n_genes", lambda s: s.quantile(0.25)),
        "q75_genes": ("n_genes", lambda s: s.quantile(0.75)),
        "share_over_5k_genes": ("n_genes", lambda s: float((s >= 5000).mean())),
    }
    if "total_counts" in frame.columns:
        aggregations["median_counts"] = ("total_counts", "median")

    out = (
        frame.groupby("platform", observed=True)
        .agg(**aggregations)
        .reset_index()
        .sort_values("median_genes", ascending=False, ignore_index=True)
    )
    out["plate_based"] = out["platform"].isin(PLATE_PLATFORMS)
    return out
