"""Reading GSE201348, against a deposit built to the shape that was verified.

The fixtures here mirror what the GEO listing and series matrix actually showed
on 2026-09-06 — 10x triplets in a tar, arm labels only in the series matrix,
`A001-C-007` sample naming, and a replicate pair. The tests that matter are the
ones where getting it wrong is silent:

* an unknown `disease stage` must STOP the run, because a dropped sample is a
  smaller cohort reported as the same cohort;
* the matrix must come out cells x genes, because a wrong transpose gives every
  cell another cell's counts and nothing raises;
* the symbol column must be column 1, because reading column 0 reports every
  panel gene as absent — the identifier-space error this repo has made four
  times.
"""

from __future__ import annotations

import gzip
import io
import pathlib
import re
import tarfile
import tempfile

import numpy as np
import pandas as pd
import pytest

from src.reference.becker_io import (
    DISEASE_STAGE_MAP,
    BeckerError,
    gene_symbols,
    pooling_key,
    read_series_matrix,
    read_triplet,
    sample_files,
)

SAMPLES = [
    ("GSM0000001", "A001-C-007", "CRC", None),
    ("GSM0000002", "A001-C-014", "Polyp", None),
    ("GSM0000003", "A001-C-023", "Unaffected", None),
    ("GSM0000004", "A002-C-010", "Polyp", "Replicate1"),
    ("GSM0000005", "A002-C-010", "Polyp", "Replicate2"),
    ("GSM0000006", "A002-C-016", "Unaffected", None),
]


def _series_matrix(tmp_path, stages=None):
    stages = stages or [s[2] for s in SAMPLES]
    titles = [f'"{s[1]}{", " + s[3] if s[3] else ""}, snRNAseq"' for s in SAMPLES]
    lines = [
        "!Series_title\t\"Becker\"",
        "!Sample_title\t" + "\t".join(titles),
        "!Sample_geo_accession\t" + "\t".join(f'"{s[0]}"' for s in SAMPLES),
        "!Sample_characteristics_ch1\t" + "\t".join('"tissue: Colon"' for _ in SAMPLES),
        "!Sample_characteristics_ch1\t" + "\t".join(
            f'"disease stage: {v}"' for v in stages),
        "!Sample_characteristics_ch1\t" + "\t".join(
            '"familial adenomatous_polyposis: Y"' for _ in SAMPLES),
        "!Sample_characteristics_ch1\t" + "\t".join('"Sex: M"' for _ in SAMPLES),
        "!series_matrix_table_begin",
        "!series_matrix_table_end",
    ]
    path = tmp_path / "GSE201348_series_matrix.txt.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")
    return path


def _tar(tmp_path, *, n_genes=6, n_cells=5, transpose_wrong=False):
    """A tar of 10x triplets. Matrix Market from CellRanger is GENES x CELLS."""
    path = tmp_path / "GSE201348_RAW.tar"
    genes = ["ACTB", "KRT8", "EPCAM", "CDX2", "MS4A12", "GUCA2A"][:n_genes]
    with tarfile.open(path, "w") as tar:
        for gsm, sample, *_ in SAMPLES:
            rng = np.random.default_rng(abs(hash(gsm)) % 2**31)
            data = rng.integers(0, 5, (n_genes, n_cells))
            if transpose_wrong:
                data = data.T
            rows, cols = data.shape
            mm = ["%%MatrixMarket matrix coordinate integer general",
                  f"{rows} {cols} {int((data > 0).sum())}"]
            for i in range(rows):
                for j in range(cols):
                    if data[i, j]:
                        mm.append(f"{i+1} {j+1} {data[i, j]}")
            for kind, payload in (
                ("matrix.mtx", "\n".join(mm) + "\n"),
                ("barcodes.tsv", "\n".join(f"CELL{i}" for i in range(n_cells)) + "\n"),
                ("features.tsv", "\n".join(
                    f"ENSG{i:08d}\t{g}\tGene Expression" for i, g in enumerate(genes)) + "\n"),
            ):
                blob = gzip.compress(payload.encode())
                info = tarfile.TarInfo(f"{gsm}_{sample}_{kind}.gz")
                info.size = len(blob)
                tar.addfile(info, io.BytesIO(blob))
    return path


# ---------------------------------------------------------------------------
# The metadata, which is the only place the arms live
# ---------------------------------------------------------------------------


def test_the_arms_come_out_as_the_series_matrix_says():
    with tempfile.TemporaryDirectory() as d:
        frame = read_series_matrix(_series_matrix(pathlib.Path(d)))
    assert len(frame) == len(SAMPLES)
    by_sample = frame.set_index("gsm")
    assert by_sample.loc["GSM0000002", "arm"] == "tumour"      # Polyp
    assert by_sample.loc["GSM0000003", "arm"] == "normal"      # Unaffected
    assert by_sample.loc["GSM0000001", "arm"] is None          # CRC, excluded
    assert set(frame["donor"]) == {"A001", "A002"}


def test_replicates_are_identified_and_share_a_sample_id():
    with tempfile.TemporaryDirectory() as d:
        frame = read_series_matrix(_series_matrix(pathlib.Path(d)))
    reps = frame[frame["replicate"].notna()]
    assert len(reps) == 2
    assert set(reps["replicate"]) == {"Replicate1", "Replicate2"}
    assert reps["sample_id"].nunique() == 1, (
        "the two replicates are one physical sample and must collapse under "
        "sample_id — Amendment 1 pools them"
    )


def test_an_unknown_disease_stage_stops_the_run():
    """THE ONE THAT MATTERS. A dropped sample is a smaller cohort reported as
    the same cohort, and nothing about that is visible downstream."""
    stages = ["CRC", "Polyp", "Unaffected", "Polyp", "Polyp", "Adenocarcinoma"]
    with tempfile.TemporaryDirectory() as d:
        path = _series_matrix(pathlib.Path(d), stages=stages)
        with pytest.raises(BeckerError, match="unknown disease stage"):
            read_series_matrix(path)


def test_crc_is_excluded_by_a_recorded_decision_not_by_omission():
    assert "CRC" in DISEASE_STAGE_MAP
    assert DISEASE_STAGE_MAP["CRC"] is None


def test_a_sample_whose_donor_cannot_be_read_stops_the_run():
    from src.reference import becker_io
    with tempfile.TemporaryDirectory() as d:
        path = _series_matrix(pathlib.Path(d))
        original = becker_io._SAMPLE
        try:
            becker_io._SAMPLE = re.compile(r"^(?P<donor>ZZZ)-(?P<sample>.+)$")
            with pytest.raises(BeckerError, match="donor-sample pattern"):
                read_series_matrix(path)
        finally:
            becker_io._SAMPLE = original


# ---------------------------------------------------------------------------
# The matrices
# ---------------------------------------------------------------------------


def test_the_matrix_comes_out_cells_by_genes():
    """A wrong transpose gives every cell another cell's counts, silently."""
    with tempfile.TemporaryDirectory() as d:
        tar = _tar(pathlib.Path(d), n_genes=6, n_cells=5)
        files = sample_files(tar)
        counts, barcodes, features = read_triplet(tar, files.iloc[0])
    assert counts.shape == (5, 6), "cells x genes"
    assert len(barcodes) == 5
    assert len(features) == 6


def test_a_matrix_of_the_wrong_orientation_is_caught_not_transposed_blindly():
    with tempfile.TemporaryDirectory() as d:
        tar = _tar(pathlib.Path(d), n_genes=6, n_cells=5, transpose_wrong=True)
        files = sample_files(tar)
        with pytest.raises(BeckerError, match="another cell's counts"):
            read_triplet(tar, files.iloc[0])


def test_every_sample_in_the_tar_has_a_complete_triplet():
    with tempfile.TemporaryDirectory() as d:
        files = sample_files(_tar(pathlib.Path(d)))
    assert len(files) == len(SAMPLES)
    assert files["complete"].all()


def test_symbols_come_from_column_one_not_column_zero():
    """Column 0 is Ensembl. Reading it reports every panel gene as absent."""
    features = pd.DataFrame({0: ["ENSG00000075624"], 1: ["ACTB"], 2: ["Gene Expression"]})
    assert gene_symbols(features)[0] == "ACTB"

    with pytest.raises(BeckerError, match="column 0 is Ensembl"):
        gene_symbols(pd.DataFrame({0: ["ENSG00000075624"]}))


# ---------------------------------------------------------------------------
# The estimand decision
# ---------------------------------------------------------------------------


def test_pooling_has_no_default_because_the_two_are_different_estimands():
    frame = pd.DataFrame({"donor": ["A001", "A001"],
                          "sample_id": ["A001-C-014", "A001-C-023"]})
    assert list(pooling_key(frame, pool_by="donor")) == ["A001", "A001"]
    assert pooling_key(frame, pool_by="lesion").nunique() == 2
    with pytest.raises(BeckerError, match="different estimands"):
        pooling_key(frame, pool_by="patient")


# ---------------------------------------------------------------------------
# The gate, driven through the real deposit shape
# ---------------------------------------------------------------------------


def test_inspect_reports_the_cohort_without_applying_the_mapping_silently():
    """A human must see the arm counts before any analysis does."""
    from src.reference.jobs.becker_feasibility import inspect_deposit

    with tempfile.TemporaryDirectory() as d:
        root = pathlib.Path(d)
        report = inspect_deposit(_tar(root), _series_matrix(root))

    assert report["n_samples_in_tar"] == len(SAMPLES)
    assert report["arm_counts"]["tumour"] == 3      # 3 Polyp
    assert report["arm_counts"]["normal"] == 2      # 2 Unaffected
    assert report["disease_stage_counts"]["CRC"] == 1
    assert report["n_donors"] == 2
    assert report["n_donors_with_both_arms"] == 2
    assert report["replicate_samples"] == ["A002-C-010"]
    assert report["first_sample_shape_cells_by_genes"] == [5, 6]
    assert set(report["panel_genes_found"]) == {
        "ACTB", "KRT8", "EPCAM", "CDX2", "MS4A12", "GUCA2A"}
    assert report["in_tar_not_metadata"] == []
    assert report["in_metadata_not_tar"] == []


def test_the_deposit_reader_drops_crc_and_pools_as_amendment_1_fixed():
    from src.reference.jobs.becker_feasibility import _read_deposit

    with tempfile.TemporaryDirectory() as d:
        root = pathlib.Path(d)
        tar, series = _tar(root), _series_matrix(root)
        counts, index, keys, absent = _read_deposit(tar, series, pool_by="donor")
        # 5 scored samples (1 CRC dropped) x 5 cells
        assert counts.shape == (25, 6)
        assert absent == []
        assert set(keys) == {"A001", "A002"}, "pooled per donor"

        _, _, lesion_keys, _ = _read_deposit(tar, series, pool_by="lesion")
        assert len(set(lesion_keys)) == 4, (
            "per lesion, with the two replicates of A002-C-010 collapsing to "
            "one sample_id"
        )


def test_a_sample_with_a_different_gene_index_stops_the_stack():
    """Concatenating across a changed reference misaligns every gene, silently."""
    from src.reference.jobs.becker_feasibility import FeasibilityError, _read_deposit

    with tempfile.TemporaryDirectory() as d:
        root = pathlib.Path(d)
        series = _series_matrix(root)
        # rebuild the tar with one sample carrying a shuffled feature list
        path = root / "GSE201348_RAW.tar"
        genes = ["ACTB", "KRT8", "EPCAM", "CDX2", "MS4A12", "GUCA2A"]
        with tarfile.open(path, "w") as tar_out:
            for n, (gsm, sample, *_rest) in enumerate(SAMPLES):
                order = list(reversed(genes)) if n == 2 else genes
                mm = ["%%MatrixMarket matrix coordinate integer general",
                      "6 5 6"] + [f"{i+1} 1 3" for i in range(6)]
                for kind, payload in (
                    ("matrix.mtx", "\n".join(mm) + "\n"),
                    ("barcodes.tsv", "\n".join(f"C{i}" for i in range(5)) + "\n"),
                    ("features.tsv", "\n".join(
                        f"ENSG{i:08d}\t{g}\tGene Expression"
                        for i, g in enumerate(order)) + "\n"),
                ):
                    blob = gzip.compress(payload.encode())
                    info = tarfile.TarInfo(f"{gsm}_{sample}_{kind}.gz")
                    info.size = len(blob)
                    tar_out.addfile(info, io.BytesIO(blob))

        with pytest.raises(FeasibilityError, match="different gene index"):
            _read_deposit(path, series, pool_by="donor")
