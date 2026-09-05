"""Can six genes be pulled from the atlas without downloading 32.7 GB?

    python -m src.reference.jobs.probe_icbi_layout

The coexpression reading needs ACTB, KRT8, EPCAM, CDX2, MS4A12 and GUCA2A over
the epithelial cells of ~136 paired patients. That is a tiny slice of a 32.7 GB
object, and `HTTPRangeFile` already reads `/obs` over range requests touching
about 0.1% of it.

Whether the same works for expression depends entirely on how `X` is stored, and
that is a fact about the file rather than something to reason about:

    CSR (cells x genes)  the anndata default. Cell-major, so each cell is
                         contiguous -- which is the right way round, since we
                         want a subset of CELLS. Still not worth a partial read:
                         HDF5 compresses in chunks, so a range read pulls whole
                         chunks anyway.
    CSC (genes x cells)  six columns are six contiguous slices. A few hundred
                         megabytes, and the disk problem disappears.
    dense                shape and dtype decide it; a dense float32 of 4.26M x
                         ~30k is ~500 GB, so it will not be dense.

This probe reads the group metadata ONLY -- shapes, dtypes, encoding attributes,
and the var index. It does not read `X` itself. It exists because
`/projectnb/rise-batteries` has 15 GB free against a 50 GB quota, and the
alternative to knowing this is a 33 GB download onto a filesystem with 2 GB of
margin.
"""

from __future__ import annotations

import argparse
import logging
import sys
from collections.abc import Sequence

from src.reference.icbi import ATLAS_URL, HTTPRangeFile, ICBIError, _decode

log = logging.getLogger(__name__)

#: The genes the coexpression reading scores. Six of roughly thirty thousand.
WANTED: tuple[str, ...] = ("ACTB", "KRT8", "EPCAM", "CDX2", "MS4A12", "GUCA2A")


def describe(node) -> str:
    """Shape, dtype and encoding for one h5 node, without reading its data."""
    import h5py

    if isinstance(node, h5py.Dataset):
        return f"dataset shape={node.shape} dtype={node.dtype}"
    attrs = {k: node.attrs[k] for k in node.attrs}
    encoding = attrs.get("encoding-type", "?")
    shape = attrs.get("shape", "?")
    keys = list(node.keys())
    return f"group encoding={encoding!r} shape={shape} keys={keys}"


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", default=ATLAS_URL)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    import h5py

    handle = HTTPRangeFile(args.url)
    log.info("opening %s over range requests (no download)", args.url)

    with h5py.File(handle, mode="r", driver="fileobj") as h5:
        log.info("\ntop level: %s", list(h5.keys()))

        for name in ("X", "raw", "layers"):
            if name in h5:
                log.info("\n/%s  ->  %s", name, describe(h5[name]))
                node = h5[name]
                if hasattr(node, "keys"):
                    for key in list(node.keys())[:8]:
                        log.info("    /%s/%s  ->  %s", name, key, describe(node[key]))

        if "var" not in h5:
            raise ICBIError("no /var group -- gene names cannot be located")
        var = h5["var"]
        log.info("\n/var keys: %s", list(var.keys()))

        index_key = var.attrs.get("_index", "_index")
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        if index_key not in var:
            log.error("var index %r not found; keys are %s", index_key, list(var.keys()))
            return 2

        # LOOK IN EVERY COLUMN THAT COULD HOLD A SYMBOL, not just the index.
        # `_index` here is Ensembl ids, and `/var` carries `GeneSymbol`
        # separately -- checking only the index reported all six genes absent
        # from an atlas that contains all six. Same identifier-space class as
        # the two the WMHS appendix records, made inside a probe written to
        # stop exactly this kind of guessing.
        candidates = [index_key] + [
            k for k in ("GeneSymbol", "gene_symbol", "var_names", "ensembl", "Geneid")
            if k in var and k != index_key
        ]
        log.info("\nthe six the coexpression reading needs, per candidate column:")
        for key in candidates:
            try:
                values = _decode(var, key)
            except Exception as exc:            # noqa: BLE001 - report, do not stop
                log.info("  %-14s unreadable: %s", key, exc)
                continue
            lookup = {str(s): i for i, s in enumerate(values)}
            hits = {g: lookup[g] for g in WANTED if g in lookup}
            log.info("  %-14s %d/%d found  e.g. %s", key, len(hits), len(WANTED),
                     list(values[:3]))
            if hits:
                for gene in WANTED:
                    log.info("      %-8s %s", gene,
                             hits.get(gene, "NOT FOUND"))

        # The verdict.
        encoding = ""
        if "X" in h5 and hasattr(h5["X"], "attrs"):
            raw_encoding = h5["X"].attrs.get("encoding-type", "")
            encoding = (raw_encoding.decode() if isinstance(raw_encoding, bytes)
                        else str(raw_encoding))
        log.info("\n%s", "=" * 68)
        if "csc" in encoding.lower():
            log.info("CSC: gene-major. Six genes are six contiguous column slices.")
            log.info("A subset pull is feasible -- no 32.7 GB download needed.")
        elif "csr" in encoding.lower():
            log.info("CSR: cell-major, so each CELL is contiguous. That is the")
            log.info("right way round for us -- the coexpression reading wants")
            log.info("the epithelial cells of the paired patients, about 16%% of")
            log.info("the rows, not all of them.")
            log.info("")
            log.info("A partial read is therefore possible in principle, but not")
            log.info("worth engineering: HDF5 compresses in chunks, so a range")
            log.info("read pulls whole chunks rather than the rows you asked for,")
            log.info("and the download is a one-off that serves all 14 studies.")
            log.info("Fetch it. That needs ~35 GB -- /project/rise-batteries,")
            log.info("not /projectnb.")
        else:
            log.info("X encoding is %r. Read the shapes above before deciding;", encoding)
            log.info("a `raw` or `layers` entry may be stored differently from X.")
        log.info("%s", "=" * 68)

        report = handle.report() if hasattr(handle, "report") else {}
        if report:
            log.info("bytes fetched by this probe: %s", report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
