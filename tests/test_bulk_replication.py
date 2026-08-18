"""W3.8 — the replication cohort loader.

The correctness property that matters here is **positional alignment**. GEO
stores sample characteristics as one line per field, with values in column
order, and nothing in the file ties a value back to its sample. Get the
alignment wrong and every patient silently receives someone else's MMR status —
the analysis still runs, the numbers still look plausible, and the conclusion is
garbage. :func:`test_characteristics_stay_attached_to_their_own_sample` is the
guard.
"""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pytest

from src.bulk.replication import (
    ReplicationError,
    collapse_probes_to_symbols,
    fold_change_vs_normal,
    parse_platform_table,
    parse_series_matrix,
    strata,
)

TAB = "\t"


def _characteristics(field: str, values: list[str]) -> str:
    """One GEO characteristics line: values in sample column order."""
    return "!Sample_characteristics_ch1" + TAB + TAB.join(f'"{field}: {v}"' for v in values)


def _row(cells: list[str]) -> str:
    return TAB.join(cells)


SERIES = "\n".join(
    [
        _row(["!Series_title", '"A small cohort"']),
        _row(["!Sample_title", '"S1"', '"S2"', '"S3"', '"S4"']),
        _row(["!Sample_geo_accession", '"GSM1"', '"GSM2"', '"GSM3"', '"GSM4"']),
        _characteristics(
            "dataset", ["discovery", "validation", "Non Tumoral", "discovery"]
        ),
        _characteristics("mmr.status", ["dMMR", "pMMR", "N/A", "pMMR"]),
        _characteristics("cimp.status", ["+", "-", "N/A", "-"]),
        _characteristics("tumor.location", ["proximal", "distal", "N/A", "distal"]),
        "!series_matrix_table_begin",
        _row(['"ID_REF"', '"GSM1"', '"GSM2"', '"GSM3"', '"GSM4"']),
        _row(['"1007_s_at"', "5.0", "6.0", "11.0", "6.5"]),
        _row(['"1053_at"', "2.0", "2.5", "9.0", "3.0"]),
        _row(['"117_at"', "7.0", "7.5", "7.2", "7.4"]),
        _row(['"121_at"', "1.0", "1.2", "1.1", "1.3"]),
        "!series_matrix_table_end",
        "",
    ]
)

PLATFORM = "\n".join(
    [
        "^PLATFORM = GPL570",
        "#ID = Affymetrix Probe Set ID",
        _row(["ID", "GB_ACC", "Gene Symbol", "Gene Title"]),
        _row(["1007_s_at", "U48705", "GUCA2A", "guanylate cyclase activator 2A"]),
        _row(["1053_at", "M87338", "GUCA2A", "guanylate cyclase activator 2A"]),
        _row(["117_at", "X51757", "CDX2", "caudal type homeobox 2"]),
        _row(["121_at", "X69699", "PAX8 /// PAX8-AS1", "multi-mapping probe"]),
        "",
    ]
)


@pytest.fixture
def series_path(tmp_path):
    path = tmp_path / "series_matrix.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(SERIES)
    return path


@pytest.fixture
def platform_path(tmp_path):
    path = tmp_path / "GPL570_table.txt"
    path.write_text(PLATFORM, encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# Series matrix
# ---------------------------------------------------------------------------


def test_characteristics_stay_attached_to_their_own_sample(series_path):
    """THE test. GEO gives no per-sample key for characteristics — only column
    order. A misalignment gives every patient someone else's MMR status and
    nothing raises."""
    _, metadata = parse_series_matrix(series_path)
    assert metadata.loc["GSM1", "mmr.status"] == "dMMR"
    assert metadata.loc["GSM2", "mmr.status"] == "pMMR"
    assert metadata.loc["GSM1", "cimp.status"] == "+"
    assert metadata.loc["GSM2", "cimp.status"] == "-"
    assert metadata.loc["GSM1", "dataset"] == "discovery"
    assert metadata.loc["GSM3", "dataset"] == "Non Tumoral"


def test_geo_missing_becomes_nan_not_a_category(series_path):
    """An "N/A" left in place becomes a third MMR level that a model will fit."""
    _, metadata = parse_series_matrix(series_path)
    assert pd.isna(metadata.loc["GSM3", "mmr.status"])
    assert pd.isna(metadata.loc["GSM3", "cimp.status"])
    assert "N/A" not in set(metadata["mmr.status"].dropna())


def test_expression_table_is_probes_by_samples(series_path):
    expression, metadata = parse_series_matrix(series_path)
    assert list(expression.columns) == list(metadata.index)
    assert expression.shape == (4, 4)
    assert expression.loc["1007_s_at", "GSM3"] == 11.0


def test_a_file_with_no_table_is_refused(tmp_path):
    path = tmp_path / "broken.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write('!Series_title\t"nothing here"\n')
    with pytest.raises(ReplicationError, match="no expression table"):
        parse_series_matrix(path)


# ---------------------------------------------------------------------------
# Platform annotation
# ---------------------------------------------------------------------------


def test_multi_mapping_probes_are_dropped_not_assigned_to_the_first_gene(platform_path):
    """Affymetrix writes these as "A /// B". Taking the first would attribute a
    shared probe's signal to one arbitrary gene."""
    mapping = parse_platform_table(platform_path)
    assert "121_at" not in mapping.index
    assert set(mapping) == {"GUCA2A", "CDX2"}


def test_platform_table_header_is_located_not_assumed(platform_path):
    """The file has comment lines above the header; hard-coding a skip count
    would break on any release that adds one."""
    assert parse_platform_table(platform_path).loc["117_at"] == "CDX2"


# ---------------------------------------------------------------------------
# Probe collapse
# ---------------------------------------------------------------------------


def test_duplicate_probes_collapse_to_the_highest_mean(series_path, platform_path):
    """Two probes map to GUCA2A. Averaging would dilute a real signal with a
    dead probe; summing would inflate it. Keep the dominant one."""
    expression, _ = parse_series_matrix(series_path)
    genes, counts = collapse_probes_to_symbols(
        expression, parse_platform_table(platform_path)
    )
    assert list(genes.index) == ["CDX2", "GUCA2A"]
    # 1007_s_at has the higher mean, so it wins over 1053_at.
    assert genes.loc["GUCA2A", "GSM1"] == 5.0
    assert counts["probes"] == 4
    assert counts["probes_with_a_unique_symbol"] == 3
    assert counts["unique_genes"] == 2


# ---------------------------------------------------------------------------
# Strata
# ---------------------------------------------------------------------------


def test_strata_separate_tumour_from_normal_mucosa(series_path):
    _, metadata = parse_series_matrix(series_path)
    masks = strata(metadata)
    assert masks["tumour"].sum() == 3  # discovery + validation
    assert masks["normal_mucosa"].sum() == 1


def test_cimp_and_mmr_strata_are_tumour_only(series_path):
    """The normal-mucosa sample has no MMR call; it must not leak into either
    subgroup."""
    _, metadata = parse_series_matrix(series_path)
    masks = strata(metadata)
    assert masks["tumour|dMMR"].sum() == 1
    assert masks["tumour|pMMR"].sum() == 2
    assert masks["tumour|CIMP+"].sum() == 1
    assert (masks["tumour|CIMP+"] & masks["normal_mucosa"]).sum() == 0


# ---------------------------------------------------------------------------
# Fold change
# ---------------------------------------------------------------------------


def test_fold_change_is_computed_on_the_log_scale(series_path, platform_path):
    """Both cohorts are log2, so a difference of medians is a log2 fold change.
    That is what makes an array cohort comparable with an RNA-seq one even
    though the units are not."""
    expression, metadata = parse_series_matrix(series_path)
    genes, _ = collapse_probes_to_symbols(expression, parse_platform_table(platform_path))
    result = fold_change_vs_normal(genes, strata(metadata), "GUCA2A")
    assert result["normal_median_log2"] == 11.0
    assert result["log2_fold_change"] == pytest.approx(6.0 - 11.0, abs=1e-6)
    assert result["fold_change"] == pytest.approx(2**-5.0, abs=1e-4)


def test_a_missing_gene_returns_nan_rather_than_raising(series_path, platform_path):
    expression, metadata = parse_series_matrix(series_path)
    genes, _ = collapse_probes_to_symbols(expression, parse_platform_table(platform_path))
    result = fold_change_vs_normal(genes, strata(metadata), "NOT_A_GENE")
    assert np.isnan(result["fold_change"])
