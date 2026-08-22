#!/usr/bin/env python
"""Dump the structure of an HDF5 matrix file. W1, week 1.

GSE178341's supplementary structure is awkward and the README budgets a day just
for parsing. Before writing a loader, find out what is actually in the file:
group layout, dataset shapes, dtypes, and whether the matrix is CSR, CSC, dense,
10x-style, or a MATLAB export.

    python src/reference/jobs/inspect_h5.py data/raw/GSE178341/*.h5

Read-only and cheap — it reads metadata and a small corner of the data, never
the whole matrix. Safe on a login node.
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np

MAX_DEPTH = 4
PREVIEW = 5


def describe(name: str, obj, depth: int = 0) -> None:
    pad = "  " * depth
    if isinstance(obj, h5py.Dataset):
        info = f"{pad}{name}  shape={obj.shape}  dtype={obj.dtype}"
        if obj.compression:
            info += f"  [{obj.compression}]"
        print(info)
        # A tiny sample tells you integer-vs-float and the identifier flavour
        # (symbols vs Ensembl) without reading anything large.
        try:
            if obj.size and obj.ndim == 1:
                head = obj[:PREVIEW]
                if head.dtype.kind in "SO":
                    head = [x.decode() if isinstance(x, bytes) else x for x in head]
                print(f"{pad}    head: {list(head)}")
        except Exception as exc:  # a corrupt or exotic dataset must not stop the dump
            print(f"{pad}    (unreadable: {exc})")
    else:
        print(f"{pad}{name}/")
        if depth < MAX_DEPTH:
            for key in obj:
                describe(key, obj[key], depth + 1)


def main(paths: list[str]) -> int:
    if not paths:
        print(__doc__)
        return 1

    for raw in paths:
        path = Path(raw)
        if not path.exists():
            print(f"!! missing: {path}")
            continue

        print("=" * 72)
        print(f"{path}   ({path.stat().st_size:,} bytes)")
        print("=" * 72)

        with h5py.File(path, "r") as handle:
            if handle.attrs:
                print("root attrs:")
                for key, value in handle.attrs.items():
                    print(f"  {key} = {value!r}")
            for key in handle:
                describe(key, handle[key])

            # Sparse matrices are the common case and the one worth calling out,
            # because the orientation decides whether cells are rows or columns.
            for key in handle:
                node = handle[key]
                if isinstance(node, h5py.Group) and {"data", "indices", "indptr"} <= set(node):
                    shape = node.attrs.get("shape") or node.get("shape")
                    shape = np.asarray(shape[:]) if hasattr(shape, "__len__") else shape
                    fmt = node.attrs.get("encoding-type") or node.attrs.get("h5sparse_format")
                    print(f"\n>> sparse matrix at /{key}: shape={shape} format={fmt!r}")
                    sample = node["data"][: min(1000, node["data"].shape[0])]
                    integral = bool(np.allclose(sample, np.rint(sample)))
                    print(f"   first 1000 values integral: {integral}")
                    print(f"   nnz: {node['data'].shape[0]:,}")
        print()

    print("Paste this output back and it becomes src/reference/ingest.py's loader.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
