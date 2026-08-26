"""The shared gene index. W1 emits it; W3 conforms to it. Week 1.

execution_plan.md §3.4: "W1 emits `S_matrix_{rung}_{version}.parquet` with a
fixed gene index. W3 emits bulk on the same index. Integration is then a join,
not a negotiation." This module is that index.

Two decisions it settles
------------------------
Open decision #2 asks who produces it — W1, from the GSE178341 feature table,
because single-cell features are the more constrained set. Open decision #3 asks
what it is keyed on, and the answer is **unversioned Ensembl IDs**, with symbols
as a mapped column:

- TCGA STAR counts arrive as versioned Ensembl (``ENSG00000141510.16``)
- this deposit arrives versioned *and* CellRanger-suffixed
  (``ENSG00000243485.5_4``)
- the frozen panel and both labelling axes are written as symbols

The unversioned ID is the only form all three share. Joining on symbols instead
is the usual way a join silently loses several percent of genes, because symbols
are renamed between annotation releases while Ensembl IDs are stable.

The build W3 needs to know about
--------------------------------
GSE178341 is aligned to **GRCh37_liftover_v28** — an hg19 liftover of GENCODE
v28. TCGA STAR counts are GRCh38. Ensembl IDs are mostly stable across builds, so
the join works, but it is not a clean match: genes retired or reassigned between
builds will not join, and that loss is real rather than a bug. The genome is
recorded in the map file so nobody has to remember.

Duplicates
----------
CellRanger's ``_N`` suffix exists because the same gene appears more than once in
the reference. Stripping it collapses those features onto one Ensembl ID, so the
index has one row per gene and the collisions are counted rather than hidden —
whoever builds a matrix on this index must decide whether to sum the duplicate
features or drop all but one, and that decision belongs to them, not here.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from src.common.panel import axis_genes, panel_genes
from src.common.paths import CONFIG_DIR

#: Written into the map file so the build travels with the index.
GENOME_NOTE = (
    "GSE178341 is aligned to GRCh37_liftover_v28 (hg19 liftover of GENCODE v28). "
    "TCGA STAR counts are GRCh38. Ensembl IDs are mostly stable across builds so "
    "the join works, but genes retired or reassigned between builds will not "
    "join. See docs/open_decisions.md #3."
)

MAP_COLUMNS = ("ensembl_id", "gene_symbol", "feature_id", "genome", "n_features")


class GeneIndexError(ValueError):
    """The gene index could not be built from the input given."""


def build_gene_index(var: pd.DataFrame) -> tuple[list[str], pd.DataFrame, dict[str, Any]]:
    """Collapse a feature table onto one row per unversioned Ensembl ID.

    `var` is :func:`src.reference.ingest.read_gse178341_index`'s second return
    value, carrying ``ensembl_id``, ``gene_symbol``, ``feature_id`` and
    ``genome``.

    Returns ``(index, mapping, report)``. `index` is the ordered identifier list
    written one per line; `mapping` carries the symbol and provenance; `report`
    holds the counts worth reading before anyone builds on it.
    """
    required = {"ensembl_id", "gene_symbol", "feature_id"}
    missing = required - set(var.columns)
    if missing:
        raise GeneIndexError(f"var is missing column(s): {sorted(missing)}")
    if var.empty:
        raise GeneIndexError("var is empty")

    frame = var.reset_index(drop=True).copy()
    frame["ensembl_id"] = frame["ensembl_id"].astype(str)
    frame["gene_symbol"] = frame["gene_symbol"].astype(str)

    counts = frame.groupby("ensembl_id", observed=True).size().rename("n_features")
    # First occurrence wins, so the order follows the deposit's own feature order
    # and is reproducible.
    mapping = (
        frame.drop_duplicates("ensembl_id", keep="first")
        .merge(counts, on="ensembl_id")
        .loc[:, [c for c in MAP_COLUMNS if c in frame.columns or c == "n_features"]]
        .reset_index(drop=True)
    )

    report = {
        "n_features_in": int(len(frame)),
        "n_genes_out": int(len(mapping)),
        "n_collapsed": int(len(frame) - len(mapping)),
        "n_duplicated_ids": int((counts > 1).sum()),
        "genomes": sorted(frame["genome"].astype(str).unique())
        if "genome" in frame.columns
        else [],
    }
    return list(mapping["ensembl_id"]), mapping, report


def check_panel_coverage(mapping: pd.DataFrame) -> pd.DataFrame:
    """Which frozen panel and axis genes resolve to an Ensembl ID.

    **A week-1 finding, not a week-9 surprise.** A panel gene absent from the
    index cannot be tested at all, and an axis marker absent from it silently
    weakens the labels. Returns one row per gene with its ID or a blank.
    """
    symbol_to_id: dict[str, str] = {}
    for symbol, ensembl in zip(
        mapping["gene_symbol"].astype(str), mapping["ensembl_id"].astype(str), strict=True
    ):
        symbol_to_id.setdefault(symbol, ensembl)

    rows = []
    for source, genes in (
        ("panel", panel_genes()),
        ("axis:stem_pole", axis_genes("stem_pole")),
        ("axis:opposite_lineage", axis_genes("opposite_lineage")),
    ):
        for gene in genes:
            rows.append(
                {
                    "source": source,
                    "gene_symbol": gene,
                    "ensembl_id": symbol_to_id.get(gene, ""),
                    "found": gene in symbol_to_id,
                }
            )
    return pd.DataFrame(rows)


def read_shared_index(path: Path | str) -> list[str]:
    """Read another arm's emitted index — one unversioned Ensembl ID per line."""
    lines = [
        line.strip()
        for line in Path(path).read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]
    if not lines:
        raise GeneIndexError(f"{path} holds no identifiers")
    versioned = [g for g in lines if "." in g]
    if versioned:
        raise GeneIndexError(
            f"{path} carries versioned identifiers ({versioned[0]!r} and "
            f"{len(versioned) - 1} more). Decision #3 keys the shared index on "
            f"UNVERSIONED Ensembl IDs; strip the suffix before intersecting."
        )
    return lines


def intersect_gene_index(
    mapping: pd.DataFrame, shared: list[str]
) -> tuple[list[str], pd.DataFrame, dict[str, Any]]:
    """Intersect this deposit's features with another arm's index.

    **The operative index for deconvolution is the intersection**, not either
    side's full feature set — open decision #2, corrected. A gene present in
    only one matrix cannot be deconvolved from the other.

    The decisive argument is CLAUDE.md invariant 1, in its general form. A gene
    in the bulk index with no single-cell counts has an **unknown** signature,
    not a zero one. Carrying it on the shared index and filling the S-matrix row
    with zeros asserts "this cell type does not express this gene", which is a
    measurement nobody made. `None` is not `0.0`, and an absent row is not a
    silent row.

    **Sorted ascending**, matching the convention `gene_index_0.9.0.txt` already
    uses, so the output is byte-identical whichever arm runs it. That is what
    makes the ownership question in #2 stop mattering: the file is reproducible
    from two committed inputs rather than owned by whoever ran it first.

    Raises if any frozen panel gene fails to survive the join — #2 says a panel
    gene missing from the intersection is "a reason to revisit the index, not to
    drop the gene", so this refuses rather than emitting a quietly weaker panel.
    """
    ours = mapping["ensembl_id"].astype(str)
    shared_set = set(shared)
    keep = sorted(set(ours) & shared_set)
    if not keep:
        raise GeneIndexError(
            "the two indices share no identifier — check both are unversioned "
            "Ensembl IDs on the same key"
        )

    out_mapping = (
        mapping[ours.isin(shared_set)]
        .sort_values("ensembl_id", kind="stable")
        .reset_index(drop=True)
    )

    coverage = check_panel_coverage(out_mapping)
    lost = coverage[~coverage["found"] & (coverage["source"] == "panel")]
    if len(lost):
        raise GeneIndexError(
            f"{len(lost)} frozen panel gene(s) do not survive the intersection: "
            f"{sorted(lost['gene_symbol'])}. Decision #2: a panel gene missing "
            f"from the intersection is a reason to revisit the index, not to "
            f"drop the gene. Refusing to emit."
        )

    report: dict[str, Any] = {
        "n_ours": int(ours.nunique()),
        "n_shared": len(shared_set),
        "n_intersection": len(keep),
        "n_ours_only": int(ours.nunique() - len(keep)),
        "n_shared_only": len(shared_set) - len(keep),
        "panel_found": int(
            (coverage["found"] & (coverage["source"] == "panel")).sum()
        ),
        "panel_total": int((coverage["source"] == "panel").sum()),
    }
    return keep, out_mapping, report


def write_gene_index(
    index: list[str],
    mapping: pd.DataFrame,
    *,
    version: str,
    config_dir: Path | None = None,
    genome_note: str = GENOME_NOTE,
) -> tuple[Path, Path]:
    """Write ``gene_index_{version}.txt`` and its ``.map.tsv``. Never overwrites.

    Versioning is the point: config/gene_index/README.md says to bump the version
    and commit a new file rather than editing one in place, because results
    reference the index through the sha they were written under and an in-place
    edit makes every earlier result unreproducible without saying so.
    """
    base = Path(config_dir) if config_dir is not None else CONFIG_DIR / "gene_index"
    base.mkdir(parents=True, exist_ok=True)
    index_path = base / f"gene_index_{version}.txt"
    map_path = base / f"gene_index_{version}.map.tsv"

    for path in (index_path, map_path):
        if path.exists():
            raise GeneIndexError(
                f"{path} already exists. Bump the version rather than editing in "
                f"place — earlier results reference this file by version."
            )

    index_path.write_text("\n".join(index) + "\n", encoding="utf-8")
    with open(map_path, "w", encoding="utf-8") as handle:
        handle.write(f"# gene index {version}\n")
        for line in genome_note.splitlines():
            handle.write(f"# {line}\n")
        mapping.to_csv(handle, sep="\t", index=False)
    return index_path, map_path


def read_gene_index(version: str, config_dir: Path | None = None) -> list[str]:
    """Read an emitted index back. The one W3 and build_signature() both call."""
    base = Path(config_dir) if config_dir is not None else CONFIG_DIR / "gene_index"
    path = base / f"gene_index_{version}.txt"
    if not path.exists():
        raise GeneIndexError(f"{path} does not exist; emit it first")
    return [line for line in path.read_text(encoding="utf-8").split("\n") if line]
