"""Pull one patient's cells out of the 30 GB ICBI atlas, without loading it.

The atlas is a single h5ad of 4,264,929 cells x 28,476 genes. Both ``/X`` and
``/layers/counts`` are CSR, which is cell-major -- so a set of cells is a set of
contiguous row slices, and that is the one access pattern this layout is good
at. Nothing here materialises more than the requested rows.

**READ ``layers['counts']``, NEVER ``X``.** Measured on the published object:

    /X                min 0.2795, max 5.404, 3525 distinct values  -> log1p
    /layers/counts    min 1, max 289, 155 distinct, [1,1,1,1,2,1,2,5]  -> raw

``adata.X`` is what any obvious code reaches for, and the coexpression reading
is built on raw integer counts -- detection at >= 1 UMI, CP10K, depth matching.
Run it against log1p values and every number is wrong while nothing raises. So
:func:`read_cells` defaults to the counts layer and
:func:`assert_raw_counts` checks the values rather than trusting the name.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


class SliceError(RuntimeError):
    """The atlas cannot be sliced as asked."""


#: The layer holding raw integer counts. See the module docstring.
COUNTS_LAYER = "layers/counts"

#: `atlas_cell_type_coarse` -> this project's compartment vocabulary.
#:
#: `Cancer cell` maps to epithelial deliberately: in the diseased arm the mature
#: colonocytes the reading is about carry that label, and dropping it takes most
#: of the tumour-arm epithelium with it (7 usable patients against 136).
COMPARTMENT_MAP: dict[str, str] = {
    "Epithelial cell": "epithelial",
    "Cancer cell": "epithelial",
    "T cell": "immune",
    "B cell": "immune",
    "Plasma cell": "immune",
    "Myeloid cell": "immune",
    "Granulocyte": "immune",
    "NK cell": "immune",
    "ILC": "immune",
    "Mast cell": "immune",
    "Stromal cell": "stromal",
    "Schwann cell": "stromal",
}

#: `sample_type` -> the two arms. Everything else is not this contrast:
#: `healthy normal` is a different donor, and metastasis / polyp / blood /
#: lymph node are different questions.
TISSUE_MAP: dict[str, str] = {
    "primary tumor": "tumour",
    "adjacent normal": "normal",
}


def assert_raw_counts(matrix, *, context: str) -> None:
    """Refuse anything that is not raw integer counts.

    Checks the VALUES, not the layer name. Trust the layer, then verify it --
    a mislabelled or reprocessed atlas would otherwise hand log1p values to a
    detection statistic that reports them happily.
    """
    data = matrix.data if hasattr(matrix, "data") else np.asarray(matrix).ravel()
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        raise SliceError(f"{context} holds no finite values")
    if not np.allclose(finite, np.round(finite)):
        raise SliceError(
            f"{context} is not raw integer counts: min {finite.min():.4f}, "
            f"max {finite.max():.4f}, {len(np.unique(finite))} distinct values. "
            f"The atlas's /X is log1p-normalised; read {COUNTS_LAYER!r} instead. "
            f"Detection at >= 1 UMI against log values is wrong and silent."
        )
    if finite.min() < 0:
        raise SliceError(f"{context} holds negative counts")


def read_var(path: str | Path) -> pd.DataFrame:
    """Gene table. The symbols are in ``GeneSymbol``; ``_index`` is Ensembl.

    Looking only at the index reports every gene symbol absent from an atlas
    that contains them all -- the identifier-space error this repo records
    three times. Both columns are returned so a caller cannot pick wrong
    silently.
    """
    import h5py

    from src.reference.icbi import _decode

    with h5py.File(str(path), "r") as h5:
        if "var" not in h5:
            raise SliceError(f"{path} has no /var group")
        var = h5["var"]
        index_key = var.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        frame = pd.DataFrame({
            "ensembl_id": [str(v) for v in _decode(var, index_key)],
        })
        for key in ("GeneSymbol", "var_names"):
            if key in var:
                frame["gene_symbol"] = [str(v) for v in _decode(var, key)]
                break
        else:
            raise SliceError(
                f"{path}/var has no GeneSymbol or var_names column, so genes "
                f"cannot be resolved to symbols. Keys: {list(var.keys())}"
            )
    return frame


def read_cells(
    path: str | Path,
    row_index: np.ndarray,
    *,
    layer: str = COUNTS_LAYER,
    n_genes: int | None = None,
):
    """Rows ``row_index`` of a CSR h5ad group, as a scipy CSR matrix.

    Reads each requested row's ``[indptr[i], indptr[i+1])`` slab from ``data``
    and ``indices``. Rows are read in sorted order -- HDF5 fancy-indexing a
    dataset requires increasing indices, and out-of-order reads on a compressed
    dataset thrash the chunk cache -- then permuted back to the caller's order,
    so ``result[j]`` is always ``row_index[j]``.
    """
    import h5py
    import scipy.sparse as sp

    row_index = np.asarray(row_index, dtype=np.int64)
    if row_index.size == 0:
        raise SliceError("no rows requested")

    with h5py.File(str(path), "r") as h5:
        if layer not in h5:
            raise SliceError(
                f"{path} has no {layer!r}. Available: {list(h5.keys())}. "
                f"Do NOT fall back to /X -- it is log1p-normalised."
            )
        group = h5[layer]
        for key in ("data", "indices", "indptr"):
            if key not in group:
                raise SliceError(f"{layer} is not a CSR group (missing {key!r})")
        indptr = group["indptr"][:]
        width = n_genes
        if width is None:
            shape = group.attrs.get("shape")
            if shape is None:
                raise SliceError(f"{layer} has no shape attribute; pass n_genes")
            width = int(np.asarray(shape).ravel()[1])

        order = np.argsort(row_index, kind="stable")
        ordered = row_index[order]
        if ordered[-1] >= len(indptr) - 1:
            raise SliceError(
                f"row {ordered[-1]} is past the end of a {len(indptr) - 1}-row matrix"
            )

        starts, stops = indptr[ordered], indptr[ordered + 1]
        lengths = (stops - starts).astype(np.int64)
        total = int(lengths.sum())
        data = np.empty(total, dtype=np.float32)
        cols = np.empty(total, dtype=np.int32)
        at = 0
        for start, stop, length in zip(starts, stops, lengths, strict=True):
            if length == 0:
                continue
            data[at:at + length] = group["data"][start:stop]
            cols[at:at + length] = group["indices"][start:stop]
            at += length

    new_indptr = np.concatenate([[0], np.cumsum(lengths)]).astype(np.int64)
    ordered_matrix = sp.csr_matrix(
        (data, cols, new_indptr), shape=(len(ordered), width)
    )
    # Undo the sort so the caller's row order is preserved.
    inverse = np.empty_like(order)
    inverse[order] = np.arange(order.size)
    return ordered_matrix[inverse]


def compartments(cell_type: pd.Series) -> pd.Series:
    """Atlas coarse labels -> epithelial / immune / stromal.

    An unmapped label becomes ``pd.NA`` and is REPORTED rather than silently
    binned. A new cell type quietly counted as immune would shift every
    compartment fraction and no check downstream reads the raw label again.
    """
    mapped = cell_type.astype(str).map(COMPARTMENT_MAP)
    unmapped = sorted(set(cell_type.astype(str)[mapped.isna()]))
    if unmapped:
        log.warning("unmapped cell types dropped: %s", unmapped)
    return mapped


def arms(sample_type: pd.Series) -> pd.Series:
    """`sample_type` -> normal / tumour, with everything else left NA."""
    return sample_type.astype(str).map(TISSUE_MAP)
