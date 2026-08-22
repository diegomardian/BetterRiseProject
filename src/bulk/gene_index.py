"""The shared gene index, built from the GDC gene model. W3.

WHY THIS IS IN src/bulk/ AND NOT src/reference/
-----------------------------------------------
[docs/open_decisions.md #2] says W1 emits ``gene_index_1.0.0.txt`` from the
GSE178341 feature table in week 1 and W3 conforms, with a written fallback:

    "If W1's ingest slips, W3 emits a provisional index from the GDC gene model
     and W1 conforms."

W1's week-1 work landed — ingest guards, per-batch QC, ambient estimation,
labels — without the index, and nothing under ``src/`` handles Ensembl IDs at
all. The fallback has fired. This module is it.

The artifact is versioned **0.9.0**, not 1.0.0, to say exactly that: provisional,
W3-built, superseded the moment W1 emits theirs. Bumping to 1.0.0 is a decision
for the week-1 meeting, not a side effect of this file.

THE KEY IS THE UNVERSIONED ENSEMBL ID
-------------------------------------
Open decision #3, recommendation adopted here. TCGA STAR counts arrive as
versioned Ensembl (``ENSG00000000003.15``); a single-cell reference built from a
different GENCODE release carries different version suffixes for the same gene.
Joining on the versioned string silently drops those genes — the usual way a
join loses ~8% of the index. So the version suffix is *stripped into its own
column* and the join key is the bare ``ENSG``.

Symbols are a mapped column, never the key. The panel and both labelling axes
are written as symbols (``config/panel.yaml``), so every panel lookup resolves
through the map, and an unresolvable panel gene is a week-1 finding rather than
a week-9 surprise.

THE INDEX CONTAINS PANEL GENES. THAT IS DELIBERATE
--------------------------------------------------
``build_signature()`` refuses to run if a target gene appears in the index it is
handed (``src/reference/signature.py:96``, pinned by
``tests/test_leakage.py:45``). But W3's bulk matrix *must* carry GUCA2A and CDX2
— they are the outcome variables for the week-2 premise check and for the Stage
4 variance question. Both requirements are satisfiable: the index carries every
gene, and W1 passes a filtered view. :func:`target_free_index` is that view.

Removing panel genes from the index itself would make W3's own deliverable
unrepresentable on the shared index. Do not do it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from src.common.panel import panel_genes
from src.common.paths import gene_index_path

#: The version W3 emits under the open-decision-#2 fallback. W1's index, if and
#: when it lands, is 1.0.0 and supersedes this one.
PROVISIONAL_VERSION = "0.9.0"

#: GENCODE puts pseudoautosomal-region genes in twice: once on X, once on Y with
#: this suffix and an identical symbol. Dropping the _PAR_Y copy is the standard
#: deterministic resolution and it is what keeps ENSG unique.
PAR_Y_SUFFIX = "_PAR_Y"

#: STAR emits four alignment-summary rows above the gene rows in the same table.
#: They are not genes and must never reach a matrix.
STAR_SUMMARY_ROWS = ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")

#: Column order of the committed map file.
MAP_COLUMNS = (
    "ensembl_id",
    "ensembl_version",
    "gene_symbol",
    "gene_type",
    "symbol_ambiguous",
    "on_panel",
)


class GeneIndexError(RuntimeError):
    """The gene model could not be turned into a usable index."""


def gene_index_map_path(version: str) -> Path:
    """``config/gene_index/gene_index_{version}.map.tsv``.

    ``src/common/paths.py`` resolves the ``.txt`` index but has no accessor for
    the map documented alongside it in ``config/gene_index/README.md``. Defined
    here rather than added to shared code unilaterally (CONTRIBUTING §2); it
    should move to ``common/paths.py`` once someone else needs it.
    """
    return gene_index_path(version).with_name(f"gene_index_{version}.map.tsv")


# ---------------------------------------------------------------------------
# Identifier handling
# ---------------------------------------------------------------------------


def strip_version(gene_id: str) -> tuple[str, str]:
    """``"ENSG00000000003.15"`` -> ``("ENSG00000000003", "15")``.

    The ``_PAR_Y`` suffix rides on the *identifier*, not the version, so it is
    preserved in the stem: ``ENSG00000182378.14_PAR_Y`` ->
    ``("ENSG00000182378_PAR_Y", "14")``. Those rows are dropped by
    :func:`build_gene_index`, but stripping has to be lossless first or the drop
    cannot be audited.
    """
    par_y = gene_id.endswith(PAR_Y_SUFFIX)
    core = gene_id[: -len(PAR_Y_SUFFIX)] if par_y else gene_id
    stem, _, version = core.partition(".")
    if par_y:
        stem += PAR_Y_SUFFIX
    return stem, version


# ---------------------------------------------------------------------------
# Reading the GDC gene model
# ---------------------------------------------------------------------------


def read_star_counts(path: str | Path) -> pd.DataFrame:
    """Read one GDC STAR-counts TSV, dropping the alignment-summary rows.

    The file carries a ``# gene-model: GENCODE vNN`` comment, then a header, then
    four ``N_*`` summary rows, then one row per gene. Nothing here assumes a
    fixed header position — the header is located by its first field, because a
    GDC release that adds a comment line should not silently shift every column.
    """
    path = Path(path)
    header_row = None
    with open(path, encoding="utf-8") as handle:
        for i, line in enumerate(handle):
            if line.startswith("gene_id\t"):
                header_row = i
                break
    if header_row is None:
        raise GeneIndexError(
            f"{path} has no 'gene_id' header line. That is not a GDC STAR-counts "
            f"file — check the download rather than adjusting the parser."
        )

    df = pd.read_csv(path, sep="\t", skiprows=header_row, dtype={"gene_id": "string"})
    df = df.loc[~df["gene_id"].isin(STAR_SUMMARY_ROWS)].reset_index(drop=True)
    if df.empty:
        raise GeneIndexError(f"{path} has a header but no gene rows.")
    return df


def gene_model_version(path: str | Path) -> str | None:
    """The ``# gene-model: GENCODE vNN`` line, verbatim, or None if absent.

    Recorded in the manifest. Which GENCODE release the index came from is the
    single fact that determines whether a future join works.
    """
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("#"):
                return None
            if "gene-model" in line:
                return line.lstrip("#").strip()
    return None


# ---------------------------------------------------------------------------
# Building the index
# ---------------------------------------------------------------------------


@dataclass
class FilterReport:
    """How many genes each step removed. W3.1's "report how many genes are lost".

    Kept as a running list rather than a summary so the write-up can quote the
    steps in order. A filtering rule nobody can reconstruct is not documented.
    """

    steps: list[tuple[str, int, int]] = field(default_factory=list)

    def record(self, label: str, before: int, after: int) -> None:
        self.steps.append((label, before, after))

    def to_frame(self) -> pd.DataFrame:
        df = pd.DataFrame(self.steps, columns=["step", "n_before", "n_after"])
        df["n_lost"] = df["n_before"] - df["n_after"]
        return df

    def summary(self) -> str:
        return "\n".join(
            f"  {label:<34} {before:>7,} -> {after:>7,}  (-{before - after:,})"
            for label, before, after in self.steps
        )


def build_gene_index(
    star: pd.DataFrame,
    *,
    keep_biotypes: tuple[str, ...] | None = None,
) -> tuple[pd.DataFrame, FilterReport]:
    """Turn a STAR-counts frame into the index + symbol map.

    Parameters
    ----------
    star:
        Output of :func:`read_star_counts`. Only the annotation columns are
        read; expression columns are ignored.
    keep_biotypes:
        Biotypes to retain. ``None`` keeps everything, which is the default and
        the recommendation — **no expression or biotype filtering at ingest**.
        A threshold baked into the matrix is invisible downstream; biotype
        restriction belongs at signature-construction time, where W1 and W2 can
        see it. The parameter exists so that choice is explicit rather than
        unavailable.

    Returns
    -------
    (index, report)
        ``index`` has one row per unique unversioned ENSG, columns per
        :data:`MAP_COLUMNS`, sorted by ``ensembl_id`` so the committed file is
        byte-stable across runs.
    """
    required = {"gene_id", "gene_name", "gene_type"}
    missing = sorted(required - set(star.columns))
    if missing:
        raise GeneIndexError(f"STAR frame is missing annotation column(s): {missing}")

    report = FilterReport()
    n0 = len(star)

    stems_versions = [strip_version(g) for g in star["gene_id"].astype(str)]
    index = pd.DataFrame(
        {
            "ensembl_id": [s for s, _ in stems_versions],
            "ensembl_version": [v for _, v in stems_versions],
            "gene_symbol": star["gene_name"].astype("string"),
            "gene_type": star["gene_type"].astype("string"),
        }
    )
    report.record("gene rows in the gene model", n0, len(index))

    # --- PAR_Y: deterministic, and the reason ENSG can be a key ------------
    before = len(index)
    index = index.loc[~index["ensembl_id"].str.endswith(PAR_Y_SUFFIX)].copy()
    report.record("drop GENCODE _PAR_Y duplicates", before, len(index))

    # --- biotype, only if explicitly asked for -----------------------------
    if keep_biotypes is not None:
        before = len(index)
        index = index.loc[index["gene_type"].isin(keep_biotypes)].copy()
        report.record(f"biotype in {list(keep_biotypes)}", before, len(index))

    # --- a residual duplicate ENSG is a broken gene model, not a filter ----
    dupes = index["ensembl_id"].duplicated(keep=False)
    if dupes.any():
        example = sorted(index.loc[dupes, "ensembl_id"].unique())[:5]
        raise GeneIndexError(
            f"{int(dupes.sum())} rows share an unversioned Ensembl ID after "
            f"dropping _PAR_Y, e.g. {example}. The key is not unique, so the "
            f"join contract cannot hold. Investigate the gene model — do not "
            f"deduplicate blindly."
        )

    # --- symbol ambiguity: flagged, never collapsed ------------------------
    # Several ENSGs legitimately share a symbol. Collapsing them would fabricate
    # a gene; dropping them would lose real signal. Flag and keep both.
    index["symbol_ambiguous"] = index["gene_symbol"].duplicated(keep=False)
    index["on_panel"] = index["gene_symbol"].isin(set(panel_genes()))

    index = (
        index.loc[:, list(MAP_COLUMNS)]
        .sort_values("ensembl_id", kind="stable")
        .reset_index(drop=True)
    )
    return index, report


# ---------------------------------------------------------------------------
# Persistence
# ---------------------------------------------------------------------------


def write_gene_index(
    index: pd.DataFrame,
    *,
    version: str = PROVISIONAL_VERSION,
    config_dir: Path | None = None,
) -> tuple[Path, Path]:
    """Write ``gene_index_{version}.txt`` and its ``.map.tsv``. Returns both paths.

    Refuses to overwrite. ``config/gene_index/README.md``: bump the version and
    commit a new file, never edit one in place — an in-place edit makes every
    earlier result unreproducible without saying so.
    """
    idx_path = gene_index_path(version)
    map_path = gene_index_map_path(version)
    if config_dir is not None:
        idx_path = Path(config_dir) / idx_path.name
        map_path = Path(config_dir) / map_path.name

    for path in (idx_path, map_path):
        if path.exists():
            raise GeneIndexError(
                f"{path} already exists. Bump the version rather than editing an "
                f"index in place (config/gene_index/README.md)."
            )

    idx_path.parent.mkdir(parents=True, exist_ok=True)
    idx_path.write_text("\n".join(index["ensembl_id"]) + "\n", encoding="utf-8")
    index.to_csv(map_path, sep="\t", index=False)
    return idx_path, map_path


def load_gene_index(version: str = PROVISIONAL_VERSION) -> list[str]:
    """The ordered identifier list. This is what a matrix is reindexed onto."""
    path = gene_index_path(version)
    if not path.exists():
        raise GeneIndexError(
            f"{path} does not exist. Build it with "
            f"`python -m src.bulk.ingest gene-index`, or take W1's 1.0.0 if it "
            f"has landed (docs/open_decisions.md #2)."
        )
    return path.read_text(encoding="utf-8").split()


def load_gene_index_map(version: str = PROVISIONAL_VERSION) -> pd.DataFrame:
    """The identifier -> symbol map committed alongside the index."""
    path = gene_index_map_path(version)
    if not path.exists():
        raise GeneIndexError(f"{path} does not exist; build the index first.")
    return pd.read_csv(path, sep="\t", dtype={"ensembl_id": "string"})


# ---------------------------------------------------------------------------
# Symbol resolution — the part the panel needs
# ---------------------------------------------------------------------------


def resolve_symbols(
    index: pd.DataFrame, symbols: list[str]
) -> tuple[dict[str, str], list[str], dict[str, list[str]]]:
    """Resolve gene symbols to unversioned Ensembl IDs against the map.

    Returns ``(resolved, unmapped, ambiguous)``:

    - ``resolved``   symbol -> a single ensembl_id
    - ``unmapped``   symbols with no row in the index at all
    - ``ambiguous``  symbol -> every candidate id, for those with more than one

    Ambiguous symbols are **not** resolved automatically. A panel gene that maps
    to two Ensembl IDs is a decision, not a tie-break: pin it explicitly in the
    caller and write down why. Tier A/B/C/D are eleven genes; that is cheap.
    """
    lookup = index.loc[index["gene_symbol"].isin(symbols)]
    resolved: dict[str, str] = {}
    ambiguous: dict[str, list[str]] = {}
    for symbol, rows in lookup.groupby("gene_symbol", observed=True):
        ids = sorted(rows["ensembl_id"].tolist())
        if len(ids) == 1:
            resolved[str(symbol)] = ids[0]
        else:
            ambiguous[str(symbol)] = ids
    unmapped = sorted(set(symbols) - set(resolved) - set(ambiguous))
    return resolved, unmapped, ambiguous


def panel_resolution_report(index: pd.DataFrame) -> pd.DataFrame:
    """One row per panel gene: did it resolve, and to what.

    Run this in week 1. An unmapped tier-A or tier-B gene is a blocker for the
    week-2 premise check and for the falsification rule, and it is far cheaper
    to find now than after the matrix is built.
    """
    from src.common.panel import tier_of

    symbols = panel_genes()
    resolved, unmapped, ambiguous = resolve_symbols(index, symbols)
    rows = []
    for symbol in symbols:
        if symbol in resolved:
            status, ids = "resolved", [resolved[symbol]]
        elif symbol in ambiguous:
            status, ids = "ambiguous", ambiguous[symbol]
        else:
            status, ids = "unmapped", []
        rows.append(
            {
                "gene_symbol": symbol,
                "tier": tier_of(symbol),
                "status": status,
                "n_candidates": len(ids),
                "ensembl_id": ";".join(ids),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# The view W1 needs
# ---------------------------------------------------------------------------


def target_free_index(index_ids: list[str], index: pd.DataFrame, targets: list[str]) -> list[str]:
    """The index with the run's target genes removed, for ``build_signature()``.

    ``build_signature()`` asserts that no target gene appears in the index it is
    given (``src/reference/signature.py:96``). The committed index deliberately
    *does* contain them, because W3's bulk matrix needs GUCA2A and CDX2 as
    outcome variables. This bridges the two: W1 calls

        build_signature(..., gene_index=target_free_index(ids, idx, targets))

    and invariant 2 holds at the point where it matters — the reference matrix —
    without amputating the shared index for everyone else.
    """
    if not targets:
        raise ValueError(
            "targets is empty. An empty target set silently disables invariant 2 "
            "— pass the panel genes under test."
        )
    resolved, _, ambiguous = resolve_symbols(index, targets)
    drop = set(resolved.values())
    for ids in ambiguous.values():
        drop.update(ids)
    return [g for g in index_ids if g not in drop]
