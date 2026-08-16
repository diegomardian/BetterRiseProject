"""Week-1 ingest guards. W1.

The one that matters most is `TestUnfilteredDroplets`: SoupX and CellBender both
need empty droplets, a GEO deposit can be raw-but-filtered, and discovering that
in week 2 costs the ambient-correction arm and gate criterion G1 with it.
"""

from __future__ import annotations

import csv

import numpy as np
import pytest
from scipy import sparse

from src.reference.ingest import (
    EMPTY_DROPLET_UMI_THRESHOLD,
    MANIFEST_FIELDS,
    DropletProfile,
    IngestError,
    append_manifest_row,
    assert_raw_counts,
    assert_unfiltered_droplets,
    cells_by_patient_and_tissue,
    droplet_profile,
    matched_normal_report,
    sha256_file,
    verify_ingest,
)

RNG = np.random.default_rng(20260815)


def raw_counts(n_cells=200, n_genes=50, lam=3.0):
    """Poisson counts — integral, variable library size."""
    return RNG.poisson(lam, size=(n_cells, n_genes)).astype(np.int64)


def droplet_matrix(n_cells=50, n_empty=950, n_genes=40):
    """A raw droplet matrix: a few real cells in a sea of empty barcodes."""
    cells = RNG.poisson(60.0, size=(n_cells, n_genes))
    empties = RNG.poisson(0.4, size=(n_empty, n_genes))
    return np.vstack([cells, empties]).astype(np.int64)


class TestRawCounts:
    def test_poisson_counts_pass(self):
        assert_raw_counts(raw_counts())

    def test_sparse_counts_pass(self):
        assert_raw_counts(sparse.csr_matrix(raw_counts()))

    def test_log1p_is_rejected(self):
        with pytest.raises(IngestError, match="non-integer"):
            assert_raw_counts(np.log1p(raw_counts()))

    def test_cpm_is_rejected(self):
        counts = raw_counts().astype(float)
        cpm = counts / counts.sum(axis=1, keepdims=True) * 1e6
        with pytest.raises(IngestError, match="non-integer"):
            assert_raw_counts(cpm)

    def test_normalised_then_rounded_is_rejected(self):
        """Integral, but every library size identical. Not raw counts."""
        counts = raw_counts().astype(float)
        scaled = np.rint(counts / counts.sum(axis=1, keepdims=True) * 10_000)
        with pytest.raises(IngestError, match="near-constant library sizes"):
            assert_raw_counts(scaled)

    def test_negative_values_are_rejected(self):
        centred = raw_counts().astype(float) - 3.0
        with pytest.raises(IngestError, match="negative"):
            assert_raw_counts(centred)

    def test_nan_is_rejected(self):
        counts = raw_counts().astype(float)
        counts[0, 0] = np.nan
        with pytest.raises(IngestError, match="NaN or infinite"):
            assert_raw_counts(counts)

    def test_empty_matrix_is_rejected(self):
        with pytest.raises(IngestError, match="empty"):
            assert_raw_counts(np.zeros((0, 0)))

    def test_context_appears_in_the_message(self):
        with pytest.raises(IngestError, match="GSE178341/tumour"):
            assert_raw_counts(np.log1p(raw_counts()), context="GSE178341/tumour")


class TestUnfilteredDroplets:
    """The week-1 blocking check for the ambient-correction arm."""

    def test_raw_droplet_matrix_passes(self):
        profile = assert_unfiltered_droplets(droplet_matrix())
        assert profile.looks_unfiltered
        assert profile.n_barcodes == 1000

    def test_cell_filtered_matrix_is_rejected(self):
        """Raw counts, but the empty droplets have been removed."""
        cells_only = RNG.poisson(60.0, size=(500, 40)).astype(np.int64)
        assert_raw_counts(cells_only)  # passes the other guard — that is the point
        with pytest.raises(IngestError, match="looks cell-filtered"):
            assert_unfiltered_droplets(cells_only)

    def test_rejection_names_the_consequence(self):
        cells_only = RNG.poisson(60.0, size=(500, 40)).astype(np.int64)
        with pytest.raises(IngestError) as exc:
            assert_unfiltered_droplets(cells_only)
        message = str(exc.value)
        assert "CellBender" in message and "SoupX" in message and "G1" in message

    def test_profile_counts_empty_droplets(self):
        profile = droplet_profile(droplet_matrix(n_cells=50, n_empty=950))
        assert profile.n_barcodes == 1000
        assert profile.n_empty == pytest.approx(950, abs=20)
        assert 0.9 < profile.empty_fraction < 1.0

    def test_threshold_is_configurable(self):
        matrix = droplet_matrix()
        loose = droplet_profile(matrix, empty_threshold=1)
        strict = droplet_profile(matrix, empty_threshold=10_000)
        assert loose.n_empty < strict.n_empty
        assert strict.n_empty == 1000

    def test_summary_is_human_readable(self):
        summary = droplet_profile(droplet_matrix()).summary()
        assert "barcodes" in summary and str(EMPTY_DROPLET_UMI_THRESHOLD) in summary

    def test_sparse_input_works(self):
        assert_unfiltered_droplets(sparse.csr_matrix(droplet_matrix()))


class TestVerifyIngest:
    def test_both_guards_run(self):
        profile = verify_ingest(droplet_matrix())
        assert isinstance(profile, DropletProfile)

    def test_filtered_input_allowed_only_when_asked(self):
        cells_only = RNG.poisson(60.0, size=(300, 40)).astype(np.int64)
        with pytest.raises(IngestError):
            verify_ingest(cells_only)
        assert verify_ingest(cells_only, require_unfiltered=False).n_barcodes == 300

    def test_raw_count_guard_still_applies(self):
        with pytest.raises(IngestError, match="non-integer"):
            verify_ingest(np.log1p(droplet_matrix()), require_unfiltered=False)


class TestCohortTables:
    PATIENTS = ["P1"] * 5 + ["P2"] * 3 + ["P3"] * 4
    TISSUE = ["tumour"] * 3 + ["normal"] * 2 + ["tumour"] * 3 + ["tumour"] * 4

    def test_cell_counts_by_patient_and_tissue(self):
        table = cells_by_patient_and_tissue(self.PATIENTS, self.TISSUE)
        assert table.loc["P1", "tumour"] == 3
        assert table.loc["P1", "normal"] == 2
        assert table.loc["P3", "tumour"] == 4

    def test_matched_normal_report_flags_unmatched_patients(self):
        report = matched_normal_report(self.PATIENTS, self.TISSUE)
        assert bool(report.loc["P1", "matched"]) is True
        assert bool(report.loc["P2", "matched"]) is False
        assert int(report["matched"].sum()) == 1

    def test_missing_tissue_column_does_not_crash(self):
        """A cohort with no normals at all still produces a report."""
        report = matched_normal_report(["P1", "P1"], ["tumour", "tumour"])
        assert report.loc["P1", "normal"] == 0
        assert not report.loc["P1", "matched"]


class TestManifest:
    def _write(self, tmp_path, name="counts.mtx", content=b"abc"):
        path = tmp_path / name
        path.write_bytes(content)
        return path

    def test_sha256_matches_hashlib(self, tmp_path):
        import hashlib

        path = self._write(tmp_path, content=b"GSE178341")
        assert sha256_file(path) == hashlib.sha256(b"GSE178341").hexdigest()

    def test_row_is_appended_with_the_frozen_header(self, tmp_path):
        manifest = tmp_path / "manifest.csv"
        data = self._write(tmp_path)
        append_manifest_row(
            data,
            source_url="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/",
            accession="GSE178341",
            downloaded_by="bode",
            manifest_path=manifest,
        )
        with open(manifest, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        assert list(rows[0]) == list(MANIFEST_FIELDS)
        assert rows[0]["accession"] == "GSE178341"
        assert rows[0]["bytes"] == "3"
        assert rows[0]["workstream"] == "W1"

    def test_reregistering_a_file_replaces_its_row(self, tmp_path):
        """Re-downloading must not leave two rows with different checksums."""
        manifest = tmp_path / "manifest.csv"
        data = self._write(tmp_path, content=b"v1")
        common = dict(
            source_url="u", accession="GSE178341", downloaded_by="bode",
            manifest_path=manifest,
        )
        append_manifest_row(data, **common)
        data.write_bytes(b"version-two")
        row = append_manifest_row(data, **common)

        with open(manifest, newline="", encoding="utf-8") as handle:
            rows = list(csv.DictReader(handle))
        assert len(rows) == 1
        assert rows[0]["sha256"] == row["sha256"]

    def test_distinct_files_accumulate(self, tmp_path):
        manifest = tmp_path / "manifest.csv"
        common = dict(
            source_url="u", accession="GSE178341", downloaded_by="bode",
            manifest_path=manifest,
        )
        append_manifest_row(self._write(tmp_path, "a.mtx"), **common)
        append_manifest_row(self._write(tmp_path, "b.mtx"), **common)
        with open(manifest, newline="", encoding="utf-8") as handle:
            assert len(list(csv.DictReader(handle))) == 2
