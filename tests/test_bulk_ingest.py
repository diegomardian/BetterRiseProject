"""W3.1 — the gene index, the scale guards, and the barcode.

Three things are worth testing here and the rest is plumbing:

1. **The join key survives a version bump.** The whole point of stripping the
   Ensembl version is that a GENCODE release change does not silently drop genes.
2. **The scale guards actually fire.** An assertion nobody has seen fail is a
   comment.
3. **The barcode parses into the right fields.** W3.4's confounding tests rest
   entirely on TSS and plate being right, and a wrong-but-plausible parse would
   read as "plate is not confounded with stage".
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.bulk.gdc import (
    GDCError,
    barcode_frame,
    build_sample_manifest,
    deduplicate_aliquots,
    parse_barcode,
    reconcile_counts,
)
from src.bulk.gene_index import (
    GeneIndexError,
    build_gene_index,
    gene_model_version,
    panel_resolution_report,
    read_star_counts,
    resolve_symbols,
    strip_version,
    target_free_index,
    write_gene_index,
)
from src.bulk.normalise import (
    ScaleError,
    assert_counts,
    assert_linear_scale,
    assert_log_scale,
    assert_tpm,
    counts_to_log2_cpm,
    renormalise_tpm,
)
from src.common.panel import panel_genes

STAR_HEADER = (
    "gene_id\tgene_name\tgene_type\tunstranded\tstranded_first\tstranded_second\t"
    "tpm_unstranded\tfpkm_unstranded\tfpkm_uq_unstranded"
)


def _star_file(tmp_path, rows, *, gene_model="GENCODE v36", name="sample.tsv"):
    """Write a STAR-counts TSV shaped like the GDC's, summary rows included."""
    lines = [f"# gene-model: {gene_model}", STAR_HEADER]
    for code in ("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous"):
        lines.append(f"{code}\t\t\t100\t100\t100\t\t\t")
    for gene_id, symbol, biotype, count, tpm in rows:
        lines.append(
            f"{gene_id}\t{symbol}\t{biotype}\t{count}\t{count}\t{count}\t{tpm}\t0.0\t0.0"
        )
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


DEFAULT_ROWS = [
    ("ENSG00000000003.15", "TSPAN6", "protein_coding", 100, 25.0),
    ("ENSG00000160182.3", "GUCA2A", "protein_coding", 4000, 900.0),
    ("ENSG00000165556.9", "CDX2", "protein_coding", 800, 50.0),
    ("ENSG00000076242.13", "MLH1", "protein_coding", 600, 20.0),
    ("ENSG00000182378.14_PAR_Y", "PLCXD1", "protein_coding", 7, 1.0),
    ("ENSG00000182378.14", "PLCXD1", "protein_coding", 9, 4.0),
]


# ---------------------------------------------------------------------------
# Identifiers
# ---------------------------------------------------------------------------


def test_strip_version_separates_the_join_key_from_the_release():
    assert strip_version("ENSG00000000003.15") == ("ENSG00000000003", "15")


def test_strip_version_keeps_par_y_on_the_stem():
    """_PAR_Y rides on the identifier, not the version. Stripping must be
    lossless before the drop, or the drop cannot be audited."""
    assert strip_version("ENSG00000182378.14_PAR_Y") == ("ENSG00000182378_PAR_Y", "14")


def test_the_key_is_stable_across_a_gencode_version_bump():
    """The reason the index is unversioned. Same gene, two releases, one key."""
    v36, _ = strip_version("ENSG00000160182.3")
    v44, _ = strip_version("ENSG00000160182.4")
    assert v36 == v44 == "ENSG00000160182"


# ---------------------------------------------------------------------------
# Reading the gene model
# ---------------------------------------------------------------------------


def test_read_star_counts_drops_the_alignment_summary_rows(tmp_path):
    star = read_star_counts(_star_file(tmp_path, DEFAULT_ROWS))
    assert len(star) == len(DEFAULT_ROWS)
    assert not star["gene_id"].str.startswith("N_").any()


def test_read_star_counts_finds_the_header_wherever_it_is(tmp_path):
    """A GDC release that adds a comment line must not shift every column."""
    path = _star_file(tmp_path, DEFAULT_ROWS)
    path.write_text("# extra comment\n" + path.read_text(encoding="utf-8"), encoding="utf-8")
    assert len(read_star_counts(path)) == len(DEFAULT_ROWS)


def test_read_star_counts_rejects_a_file_that_is_not_star_output(tmp_path):
    path = tmp_path / "wrong.tsv"
    path.write_text("a\tb\n1\t2\n", encoding="utf-8")
    with pytest.raises(GeneIndexError, match="no 'gene_id' header"):
        read_star_counts(path)


def test_gene_model_version_is_recorded(tmp_path):
    assert "v36" in gene_model_version(_star_file(tmp_path, DEFAULT_ROWS))


# ---------------------------------------------------------------------------
# The index
# ---------------------------------------------------------------------------


def test_build_gene_index_drops_par_y_and_keeps_the_key_unique(tmp_path):
    index, report = build_gene_index(read_star_counts(_star_file(tmp_path, DEFAULT_ROWS)))
    assert not index["ensembl_id"].duplicated().any()
    assert not index["ensembl_id"].str.endswith("_PAR_Y").any()
    assert ("drop GENCODE _PAR_Y duplicates", 6, 5) in report.steps


def test_the_index_contains_panel_genes_on_purpose(tmp_path):
    """W3's bulk matrix needs GUCA2A and CDX2 as OUTCOME variables — the week-2
    premise check is about their distribution. Amputating them from the shared
    index to satisfy build_signature() would make W3's own deliverable
    unrepresentable. W1 filters at the call site instead."""
    index, _ = build_gene_index(read_star_counts(_star_file(tmp_path, DEFAULT_ROWS)))
    assert index.loc[index["gene_symbol"] == "GUCA2A", "on_panel"].all()
    assert {"GUCA2A", "CDX2"} <= set(index["gene_symbol"])


def test_target_free_index_is_what_build_signature_gets(tmp_path):
    """The bridge. src/reference/signature.py:96 refuses an index containing a
    target; tests/test_leakage.py:45 pins that. Both hold at once only because
    W1 passes this view rather than the raw index."""
    index, _ = build_gene_index(read_star_counts(_star_file(tmp_path, DEFAULT_ROWS)))
    ids = index["ensembl_id"].tolist()
    filtered = target_free_index(ids, index, ["GUCA2A", "CDX2"])
    assert "ENSG00000160182" not in filtered  # GUCA2A
    assert "ENSG00000165556" not in filtered  # CDX2
    assert "ENSG00000076242" in filtered  # MLH1 — not a target for this run
    assert len(filtered) == len(ids) - 2


def test_target_free_index_refuses_an_empty_target_set(tmp_path):
    """An empty target set silently disables invariant 2 — the same guard
    build_signature() already makes."""
    index, _ = build_gene_index(read_star_counts(_star_file(tmp_path, DEFAULT_ROWS)))
    with pytest.raises(ValueError, match="silently disables invariant 2"):
        target_free_index(index["ensembl_id"].tolist(), index, [])


def test_duplicate_symbols_are_flagged_never_collapsed(tmp_path):
    """Two Ensembl IDs, one symbol. Collapsing fabricates a gene; dropping loses
    signal. Keep both, flag both, and make the caller decide."""
    rows = [*DEFAULT_ROWS, ("ENSG00000999999.1", "TSPAN6", "protein_coding", 5, 1.0)]
    index, _ = build_gene_index(read_star_counts(_star_file(tmp_path, rows)))
    tspan = index.loc[index["gene_symbol"] == "TSPAN6"]
    assert len(tspan) == 2
    assert tspan["symbol_ambiguous"].all()


def test_ambiguous_symbols_are_not_auto_resolved(tmp_path):
    rows = [*DEFAULT_ROWS, ("ENSG00000999999.1", "GUCA2A", "protein_coding", 5, 1.0)]
    index, _ = build_gene_index(read_star_counts(_star_file(tmp_path, rows)))
    resolved, unmapped, ambiguous = resolve_symbols(index, ["GUCA2A", "MLH1"])
    assert "GUCA2A" not in resolved
    assert set(ambiguous["GUCA2A"]) == {"ENSG00000160182", "ENSG00000999999"}
    assert resolved["MLH1"] == "ENSG00000076242"
    assert unmapped == []


def test_unmapped_panel_genes_are_a_week_one_finding(tmp_path):
    """The report exists so an unresolvable tier-A gene surfaces now, not in
    week 9. Most of the panel is absent from this fixture, which is the point:
    the report says so rather than the pipeline shrugging."""
    index, _ = build_gene_index(read_star_counts(_star_file(tmp_path, DEFAULT_ROWS)))
    report = panel_resolution_report(index)
    assert len(report) == len(panel_genes())
    assert set(report.loc[report["status"] == "resolved", "gene_symbol"]) == {
        "GUCA2A", "CDX2", "MLH1",
    }
    assert report.loc[report["gene_symbol"] == "GUCA2A", "tier"].item() == "A"


def test_write_gene_index_refuses_to_overwrite(tmp_path):
    """Bump the version, never edit in place — an in-place edit makes every
    earlier result unreproducible without saying so."""
    index, _ = build_gene_index(read_star_counts(_star_file(tmp_path, DEFAULT_ROWS)))
    out = tmp_path / "config"
    out.mkdir()
    idx_path, map_path = write_gene_index(index, version="0.9.0", config_dir=out)
    assert idx_path.read_text(encoding="utf-8").split() == index["ensembl_id"].tolist()
    assert len(pd.read_csv(map_path, sep="\t")) == len(index)
    with pytest.raises(GeneIndexError, match="already exists"):
        write_gene_index(index, version="0.9.0", config_dir=out)


# ---------------------------------------------------------------------------
# Scale guards
# ---------------------------------------------------------------------------


def _counts(n_samples=3, n_genes=40, seed=0):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        rng.integers(0, 5000, size=(n_samples, n_genes)),
        index=[f"S{i}" for i in range(n_samples)],
        columns=[f"ENSG{i:011d}" for i in range(n_genes)],
    )


def test_assert_counts_rejects_normalised_input():
    with pytest.raises(ScaleError, match="non-integer"):
        assert_counts(_counts() / 7.0)


def test_assert_linear_scale_rejects_a_log_matrix():
    """The failure this prevents is silent: deconvolution on log input returns
    fractions that are wrong and plausible."""
    with pytest.raises(ScaleError, match="looks log-transformed"):
        assert_linear_scale(counts_to_log2_cpm(_counts()))


def test_assert_log_scale_rejects_a_linear_matrix():
    with pytest.raises(ScaleError, match="looks linear"):
        assert_log_scale(_counts().astype(float))


def test_the_two_guards_are_not_both_satisfiable():
    """A matrix cannot pass both. If it could, neither guard means anything."""
    counts = _counts()
    log2_cpm = counts_to_log2_cpm(counts)
    assert_log_scale(log2_cpm)
    with pytest.raises(ScaleError):
        assert_linear_scale(log2_cpm)


def test_log2_cpm_keeps_zeros_at_zero():
    counts = _counts()
    counts.iloc[0, 0] = 0
    assert counts_to_log2_cpm(counts).iloc[0, 0] == 0.0


def test_subsetting_tpm_breaks_the_sum_and_renormalising_fixes_it():
    """The specific reason GDC's TPM cannot simply be sliced onto the shared
    index: TPM is computed over the full gene set and a subset no longer sums
    to 1e6, so samples stop being comparable."""
    rng = np.random.default_rng(1)
    full = pd.DataFrame(rng.random((4, 50)), index=[f"S{i}" for i in range(4)])
    full = full.div(full.sum(axis=1), axis=0) * 1e6
    assert_tpm(full)

    subset = full.iloc[:, :30]
    with pytest.raises(ScaleError, match="must be renormalised"):
        assert_tpm(subset)
    assert_tpm(renormalise_tpm(subset))


# ---------------------------------------------------------------------------
# Barcodes — W3.4 depends on these being right
# ---------------------------------------------------------------------------


def test_parse_barcode_splits_every_technical_field():
    bc = parse_barcode("TCGA-A6-2670-01A-01R-1410-07")
    assert bc.patient_id == "TCGA-A6-2670"
    assert (bc.tss, bc.plate, bc.centre, bc.analyte) == ("A6", "1410", "07", "R")
    assert bc.sample_type_name == "primary_tumour"
    assert bc.is_tumour and not bc.is_normal_adjacent


def test_parse_barcode_reads_normal_adjacent():
    bc = parse_barcode("TCGA-AA-3697-11A-01R-1723-07")
    assert bc.sample_type_name == "normal_adjacent"
    assert bc.is_normal_adjacent and not bc.is_tumour


def test_parse_barcode_rejects_a_short_identifier():
    """A sample-level barcode has no plate or centre, and W3.4 needs both."""
    with pytest.raises(GDCError, match="not a full TCGA aliquot barcode"):
        parse_barcode("TCGA-A6-2670-01A")


def test_patient_id_is_the_join_key_to_tcga_cdr():
    """One patient, tumour and normal. Both rows carry the same patient_id —
    that is what the CDR join is on, and what makes the pairing question
    (open decision: normals are not matched to all tumours) answerable."""
    frame = barcode_frame(["TCGA-A6-2670-01A-01R-1410-07", "TCGA-A6-2670-11A-01R-1410-07"])
    assert frame["patient_id"].nunique() == 1
    assert set(frame["sample_type_name"]) == {"primary_tumour", "normal_adjacent"}


def test_deduplicate_keeps_one_aliquot_per_patient_and_sample_type():
    """The patient is the unit of inference (invariant 5). Two aliquots of one
    tumour would double-weight that patient in every downstream model."""
    files = pd.DataFrame(
        {
            "file_id": ["f1", "f2", "f3"],
            "barcode": [
                "TCGA-A6-2670-01B-01R-1410-07",
                "TCGA-A6-2670-01A-01R-1410-07",
                "TCGA-A6-2670-11A-01R-1410-07",
            ],
            "project": ["TCGA-COAD"] * 3,
        }
    )
    kept, dropped = deduplicate_aliquots(build_sample_manifest(files))
    assert len(kept) == 2  # one tumour, one normal
    assert len(dropped) == 1
    # Vial A wins over vial B, deterministically and without looking at depth.
    assert kept.loc[kept["sample_type"] == "01", "vial"].item() == "A"
    assert dropped["vial"].item() == "B"


def test_reconciliation_counts_samples_and_patients_separately():
    """They differ whenever a patient contributed both a tumour and a normal,
    which is exactly the number the portal comparison has to reconcile."""
    files = pd.DataFrame(
        {
            "file_id": ["f1", "f2", "f3"],
            "barcode": [
                "TCGA-A6-2670-01A-01R-1410-07",
                "TCGA-A6-2670-11A-01R-1410-07",
                "TCGA-AG-3574-01A-01R-1723-07",
            ],
            "project": ["TCGA-COAD", "TCGA-COAD", "TCGA-READ"],
        }
    )
    table = reconcile_counts(build_sample_manifest(files))
    coad_tumour = table.query("project == 'TCGA-COAD' and sample_type_name == 'primary_tumour'")
    assert coad_tumour["n_samples"].item() == 1
    assert table["n_samples"].sum() == 3


def test_build_sample_manifest_refuses_a_file_with_no_barcode():
    files = pd.DataFrame({"file_id": ["f1"], "barcode": [None], "project": ["TCGA-COAD"]})
    with pytest.raises(GDCError, match="no aliquot barcode"):
        build_sample_manifest(files)
