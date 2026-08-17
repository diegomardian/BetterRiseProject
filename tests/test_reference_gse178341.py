"""The GSE178341 loader. W1.

Built against the real deposited file's structure, recorded in
src/reference/ingest.py: 10x CellRanger HDF5 v2, CSC, genes x barcodes, float64
integral counts, feature ids like "ENSG00000243485.5_4", and barcodes like
"C103_T_1_1_0_c1_v2_id-AAACCTGCATGCTAGT".

The fixture reproduces that layout in miniature so the loader is tested without
1.1 GB of download.
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy import sparse

from src.reference.ingest import (
    IngestError,
    check_chemistry_agreement,
    normalise_feature_id,
    parse_barcode,
    patient_cohort_table,
    read_gse178341,
    read_gse178341_index,
    read_gse178341_metadata,
)

h5py = pytest.importorskip("h5py")

#: (patient, tissue, n_cells) — mirrors the real file's grouping by sample.
SAMPLES = [("C1", "T", 5), ("C1", "N", 4), ("C2", "T", 6), ("C3", "T", 3), ("C3", "N", 3)]
N_GENES = 12
RNG = np.random.default_rng(20260817)


def _barcodes() -> list[str]:
    out = []
    for patient, tissue, n in SAMPLES:
        for _ in range(n):
            tag = "".join(RNG.choice(list("ACGT"), size=16))
            out.append(f"{patient}_{tissue}_1_1_0_c1_v2_id-{tag}")
    return out


@pytest.fixture
def deposit(tmp_path):
    """A miniature GSE178341: returns (path, dense_genes_by_cells)."""
    barcodes = _barcodes()
    n_cells = len(barcodes)
    # Lambda chosen so each cell clears EMPTY_DROPLET_UMI_THRESHOLD (100 UMI)
    # comfortably. A shallower fixture makes every synthetic cell register as an
    # empty droplet, which silently inverts the cell_filtered assertion below.
    dense = (RNG.poisson(25.0, size=(N_GENES, n_cells))).astype(np.float64)
    csc = sparse.csc_matrix(dense)

    path = tmp_path / "GSE178341_crc10x_full_c295v4_submit.h5"
    with h5py.File(path, "w") as handle:
        handle.attrs["filetype"] = np.array([b"matrix"])
        handle.attrs["version"] = np.array([2], dtype=np.int32)
        handle.attrs["chemistry_description"] = np.array([b"Single Cell 3' v3"])
        group = handle.create_group("matrix")
        group.create_dataset("barcodes", data=np.array(barcodes, dtype="S48"))
        group.create_dataset("data", data=csc.data.astype(np.float64))
        group.create_dataset("indices", data=csc.indices.astype(np.int32))
        group.create_dataset("indptr", data=csc.indptr.astype(np.int32))
        group.create_dataset("shape", data=np.array([N_GENES, n_cells], dtype=np.int32))
        features = group.create_group("features")
        features.create_dataset("_all_tag_keys", data=np.array([b"genome"]))
        features.create_dataset(
            "id",
            data=np.array([f"ENSG{i:011d}.{i % 9 + 1}_{i % 4 + 1}" for i in range(N_GENES)],
                          dtype="S32"),
        )
        features.create_dataset(
            "name", data=np.array([f"GENE{i}" for i in range(N_GENES)], dtype="S16")
        )
        features.create_dataset(
            "genome", data=np.array([b"GRCh37_liftover_v28"] * N_GENES)
        )
        features.create_dataset(
            "feature_type", data=np.array([b"Gene Expression"] * N_GENES)
        )
    return path, dense


class TestFeatureIds:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("ENSG00000243485.5_4", "ENSG00000243485"),
            ("ENSG00000237613.2_2", "ENSG00000237613"),
            ("ENSG00000141510.16", "ENSG00000141510"),
            ("ENSG00000141510", "ENSG00000141510"),
            (b"ENSG00000186092.6_4", "ENSG00000186092"),
        ],
    )
    def test_version_and_dedup_suffix_are_stripped(self, raw, expected):
        assert normalise_feature_id(raw) == expected

    def test_a_symbol_is_left_alone(self):
        """Gene symbols must survive unharmed — GUCA2A has no version suffix."""
        assert normalise_feature_id("GUCA2A") == "GUCA2A"
        assert normalise_feature_id("MS4A12") == "MS4A12"


class TestBarcodeParsing:
    def test_real_barcode_from_the_deposit(self):
        fields = parse_barcode("C103_T_1_1_0_c1_v2_id-AAACCTGCATGCTAGT")
        assert fields["patient_id"] == "C103"
        assert fields["tissue"] == "tumour"
        assert fields["sample_id"] == "C103_T_1_1_0_c1_v2"
        assert fields["cell_barcode"] == "AAACCTGCATGCTAGT"
        assert fields["chemistry"] == "v2"

    def test_normal_tissue_code(self):
        assert parse_barcode("C103_N_1_1_0_c1_v2_id-AAACCTGCATGCTAGT")["tissue"] == "normal"

    @pytest.mark.parametrize("code,region", [("TA", "A"), ("TB", "B")])
    def test_two_region_tumours_count_as_tumour(self, code, region):
        """C130_TA and C130_TB are two tumour regions from one patient. If these
        do not map to 'tumour', the patient's tumour arm reads as empty and the
        matched-normal count silently undercounts."""
        fields = parse_barcode(f"C130_{code}_1_1_0_c1_v2_id-AAACCTGAGAGGACGG")
        assert fields["tissue"] == "tumour"
        assert fields["tumour_region"] == region
        assert fields["sample_id"] == f"C130_{code}_1_1_0_c1_v2"

    def test_plain_tumour_has_no_region(self):
        assert parse_barcode("C103_T_1_1_0_c1_v2_id-ACGT")["tumour_region"] == ""

    def test_regions_stay_distinct_samples(self):
        """Both are tumour, but TA and TB must not merge into one sample."""
        a = parse_barcode("C130_TA_1_1_0_c1_v2_id-ACGT")
        b = parse_barcode("C130_TB_1_1_0_c1_v2_id-ACGT")
        assert a["tissue"] == b["tissue"] == "tumour"
        assert a["sample_id"] != b["sample_id"]

    def test_unknown_tissue_code_passes_through(self):
        """Better an unmapped code than a silently wrong label."""
        assert parse_barcode("C1_X_1_1_0_c1_v2_id-ACGT")["tissue"] == "X"

    def test_trailing_lane_suffix_is_tolerated(self):
        assert parse_barcode("C1_T_1_1_0_c1_v2_id-ACGTACGT-1")["cell_barcode"] == "ACGTACGT"

    def test_bytes_input_works(self):
        assert parse_barcode(b"C1_T_1_1_0_c1_v2_id-ACGT")["patient_id"] == "C1"

    def test_unparseable_barcode_raises(self):
        """A cell with no patient must never enter the cohort silently."""
        with pytest.raises(ValueError, match="does not match the GSE178341 scheme"):
            parse_barcode("AAACCTGCATGCTAGT-1")


class TestIndex:
    def test_obs_and_var_shapes(self, deposit):
        path, dense = deposit
        obs, var = read_gse178341_index(path)
        assert len(obs) == dense.shape[1]
        assert len(var) == N_GENES

    def test_obs_carries_the_sample_metadata(self, deposit):
        path, _ = deposit
        obs, _ = read_gse178341_index(path)
        assert set(obs["patient_id"]) == {"C1", "C2", "C3"}
        assert set(obs["tissue"]) == {"tumour", "normal"}
        assert (obs["chemistry"] == "v2").all()

    def test_var_is_keyed_on_the_normalised_ensembl_id(self, deposit):
        path, _ = deposit
        _, var = read_gse178341_index(path)
        assert all(i.startswith("ENSG") and "." not in i and "_" not in i for i in var.index)
        assert "gene_symbol" in var.columns
        assert (var["genome"] == "GRCh37_liftover_v28").all()

    def test_index_does_not_read_the_matrix(self, deposit):
        """The cheap path: week-1 tables without loading 9 GB."""
        path, _ = deposit
        obs, _ = read_gse178341_index(path)
        counts = obs.groupby(["patient_id", "tissue"], observed=True).size()
        assert counts[("C1", "tumour")] == 5
        assert counts[("C1", "normal")] == 4
        assert counts[("C2", "tumour")] == 6

    def test_feeds_the_week_one_tables(self, deposit):
        from src.reference.ingest import cells_by_patient_and_tissue, matched_normal_report

        path, _ = deposit
        obs, _ = read_gse178341_index(path)
        table = cells_by_patient_and_tissue(obs["patient_id"], obs["tissue"])
        assert table.loc["C1", "normal"] == 4

        report = matched_normal_report(obs["patient_id"], obs["tissue"])
        assert bool(report.loc["C1", "matched"])
        assert not bool(report.loc["C2", "matched"])  # tumour only


#: Barcodes covering the cases the real cohort contains: a matched patient, a
#: tumour-only patient, and a two-region patient with a normal.
COHORT_BARCODES = (
    ["C1_T_1_1_0_c1_v2_id-" + "A" * 16] * 5
    + ["C1_N_1_1_0_c1_v2_id-" + "C" * 16] * 4
    + ["C2_T_1_1_0_c1_v3_id-" + "G" * 16] * 6
    + ["C130_TA_1_1_0_c1_v2_id-" + "T" * 16] * 3
    + ["C130_TB_1_1_0_c1_v2_id-" + "AC" * 8] * 2
    + ["C130_N_1_1_0_c1_v2_id-" + "AG" * 8] * 4
)


def _cohort_obs():
    import pandas as pd

    return pd.DataFrame([parse_barcode(b) for b in COHORT_BARCODES])


def _cohort_metadata(obs, *, break_chemistry: bool = False):
    import pandas as pd

    rows = []
    for _, cell in obs.iterrows():
        rows.append(
            {
                "PID": cell["patient_id"],
                "PatientTypeID": f"{cell['patient_id']}_{cell['tissue_code']}",
                "SPECIMEN_TYPE": cell["tissue_code"],
                "SINGLECELL_TYPE": "SC3P" + cell["chemistry"],
                "MMRStatus": "MMRp" if cell["patient_id"] == "C1" else "MMRd",
                "MLH1Status": "MLH1NoMeth" if cell["patient_id"] == "C1" else "MLH1Meth",
                "Sex": "M",
                "SOURCE_HOSPITAL": "MGH",
            }
        )
    frame = pd.DataFrame(rows, index=obs.index)
    if break_chemistry:
        frame.loc[frame.index[0], "SINGLECELL_TYPE"] = "SC3Pv3"
    return frame


class TestCohortTable:
    """The week-1 cohort deliverable, including the TA/TB correction."""

    def test_two_region_patient_counts_as_matched(self):
        """C130 has TA + TB + normal. Before the TISSUE_CODES fix its tumour arm
        read as empty and it was scored unmatched."""
        table = patient_cohort_table(_cohort_obs())
        assert table.loc["C130", "n_tumour"] == 5      # 3 TA + 2 TB
        assert table.loc["C130", "n_normal"] == 4
        assert bool(table.loc["C130", "matched"])
        assert table.loc["C130", "tumour_regions"] == "A,B"

    def test_tumour_only_patient_is_unmatched(self):
        table = patient_cohort_table(_cohort_obs())
        assert table.loc["C2", "n_normal"] == 0
        assert not bool(table.loc["C2", "matched"])

    def test_matched_count_is_the_real_power(self):
        table = patient_cohort_table(_cohort_obs())
        assert int(table["matched"].sum()) == 2       # C1 and C130, not C2

    def test_regions_are_separate_samples(self):
        table = patient_cohort_table(_cohort_obs())
        assert table.loc["C130", "n_samples"] == 3    # TA, TB, N

    def test_clinical_variables_join_per_patient(self):
        obs = _cohort_obs()
        table = patient_cohort_table(obs, _cohort_metadata(obs))
        assert table.loc["C1", "MMRStatus"] == "MMRp"
        assert table.loc["C1", "MLH1Status"] == "MLH1NoMeth"
        assert table.loc["C130", "MLH1Status"] == "MLH1Meth"


class TestChemistryCrossCheck:
    def test_agreement_returns_nothing(self):
        obs = _cohort_obs()
        assert len(check_chemistry_agreement(obs, _cohort_metadata(obs))) == 0

    def test_disagreement_is_reported(self):
        """If the barcode parse and the metatables disagree, every per-batch
        threshold built on the parse is suspect."""
        obs = _cohort_obs()
        bad = check_chemistry_agreement(obs, _cohort_metadata(obs, break_chemistry=True))
        assert len(bad) == 1
        assert bad.iloc[0]["from_barcode"] == "v2"
        assert bad.iloc[0]["from_metadata"] == "v3"

    def test_missing_column_raises(self):
        obs = _cohort_obs()
        meta = _cohort_metadata(obs).drop(columns=["SINGLECELL_TYPE"])
        with pytest.raises(IngestError, match="SINGLECELL_TYPE"):
            check_chemistry_agreement(obs, meta)


class TestMetadataLoader:
    def test_reads_and_categorises(self, tmp_path):
        obs = _cohort_obs()
        path = tmp_path / "meta.csv.gz"
        meta = _cohort_metadata(obs)
        meta.insert(0, "cellID", [f"cell{i}" for i in range(len(meta))])
        meta.to_csv(path, index=False)

        loaded = read_gse178341_metadata(path)
        assert loaded.index.name == "barcode"
        assert "MLH1Status" in loaded.columns
        assert str(loaded["MMRStatus"].dtype) == "category"

    def test_all_columns_keeps_everything(self, tmp_path):
        obs = _cohort_obs()
        path = tmp_path / "meta.csv.gz"
        meta = _cohort_metadata(obs)
        meta.insert(0, "cellID", [f"cell{i}" for i in range(len(meta))])
        meta["UnlistedExtra"] = 1
        meta.to_csv(path, index=False)
        assert "UnlistedExtra" in read_gse178341_metadata(path, all_columns=True).columns
        assert "UnlistedExtra" not in read_gse178341_metadata(path).columns


def test_fixture_cells_look_like_cells_not_empty_droplets(deposit):
    """Guards the fixture itself, and runs without anndata.

    The droplet-profile assertion in TestFullLoad depends on synthetic cells
    clearing the 100-UMI floor. Those tests skip wherever anndata is absent, so
    without this the fixture can drift shallow and the failure only appears on
    the cluster — which is exactly how it was first found.
    """
    from src.reference.ingest import EMPTY_DROPLET_UMI_THRESHOLD, droplet_profile

    _, dense = deposit
    per_cell = dense.sum(axis=0)
    assert per_cell.min() > EMPTY_DROPLET_UMI_THRESHOLD, (
        f"fixture cells carry {per_cell.min():.0f} UMI, below the "
        f"{EMPTY_DROPLET_UMI_THRESHOLD} floor — raise the Poisson lambda"
    )
    profile = droplet_profile(dense.T)
    assert profile.empty_fraction == 0.0
    assert not profile.looks_unfiltered


class TestColumnReader:
    """The run-coalescing CSC reader, tested without anndata.

    This is where a subtle bug would live: the loader turns a scattered column
    selection into a handful of contiguous HDF5 reads, and getting the output
    offsets wrong would silently misplace counts rather than crash.
    """

    def _read(self, path, columns):
        from src.reference.ingest import _read_csc_columns

        with h5py.File(path, "r") as handle:
            group = handle["matrix"]
            n_genes = int(group["shape"][0])
            indptr = group["indptr"][:]
            data, indices, new_indptr = _read_csc_columns(
                group, indptr, np.asarray(columns), dtype="float64"
            )
        return sparse.csc_matrix(
            (data, indices, new_indptr), shape=(n_genes, len(columns))
        ).toarray()

    @pytest.mark.parametrize(
        "columns",
        [
            [0],                       # single column
            [0, 1, 2, 3, 4],           # one contiguous run
            [0, 5, 10],                # three isolated columns
            list(range(9)) + list(range(15, 21)),  # two runs, gap between
            [18, 19, 20],              # the final run, at the array's edge
            list(range(21)),           # everything
        ],
    )
    def test_matches_a_direct_dense_slice(self, deposit, columns):
        path, dense = deposit
        columns = sorted(columns)
        np.testing.assert_allclose(self._read(path, columns), dense[:, columns])

    def test_empty_column_is_handled(self, deposit):
        """A barcode with zero counts must not corrupt neighbouring offsets."""
        path, dense = deposit
        with h5py.File(path, "r+") as handle:
            # Blank column 3 by rewriting the matrix with it zeroed.
            csc = sparse.csc_matrix(
                (handle["matrix/data"][:], handle["matrix/indices"][:],
                 handle["matrix/indptr"][:]),
                shape=tuple(int(x) for x in handle["matrix/shape"][:]),
            ).toarray()
            csc[:, 3] = 0
            new = sparse.csc_matrix(csc)
            del handle["matrix/data"], handle["matrix/indices"], handle["matrix/indptr"]
            handle["matrix"].create_dataset("data", data=new.data.astype(np.float64))
            handle["matrix"].create_dataset("indices", data=new.indices.astype(np.int32))
            handle["matrix"].create_dataset("indptr", data=new.indptr.astype(np.int32))
            dense = csc

        got = self._read(path, [2, 3, 4])
        np.testing.assert_allclose(got, dense[:, [2, 3, 4]])
        assert got[:, 1].sum() == 0


class TestFullLoad:
    def test_matrix_is_cells_by_genes_and_values_match(self, deposit):
        pytest.importorskip("anndata")
        path, dense = deposit
        adata = read_gse178341(path)
        assert adata.shape == (dense.shape[1], N_GENES)
        np.testing.assert_allclose(adata.X.toarray(), dense.T)

    def test_patient_subset_reads_only_those_cells(self, deposit):
        pytest.importorskip("anndata")
        path, dense = deposit
        adata = read_gse178341(path, patients=["C1"])
        assert adata.shape == (9, N_GENES)          # 5 tumour + 4 normal
        assert set(adata.obs["patient_id"]) == {"C1"}
        np.testing.assert_allclose(adata.X.toarray(), dense.T[:9])

    def test_noncontiguous_subset_is_correct(self, deposit):
        """C1 and C3 are separated by C2 — the run-coalescing path."""
        pytest.importorskip("anndata")
        path, dense = deposit
        adata = read_gse178341(path, patients=["C1", "C3"])
        assert adata.shape == (15, N_GENES)         # 9 + 6
        expected = np.vstack([dense.T[:9], dense.T[15:]])
        np.testing.assert_allclose(adata.X.toarray(), expected)

    def test_dtype_is_float32_by_default(self, deposit):
        pytest.importorskip("anndata")
        path, _ = deposit
        assert read_gse178341(path).X.dtype == np.float32

    def test_droplet_profile_records_that_it_is_cell_filtered(self, deposit):
        """Decision #8 travels with the object rather than living in someone's head."""
        pytest.importorskip("anndata")
        path, _ = deposit
        adata = read_gse178341(path)
        assert adata.uns["droplet_profile"]["cell_filtered"] is True
        assert "open_decisions #8" in adata.uns["droplet_profile"]["note"]

    def test_raw_count_guard_still_runs(self, deposit):
        """Filtered droplets are tolerated; non-integer values are not."""
        pytest.importorskip("anndata")
        path, _ = deposit
        with h5py.File(path, "r+") as handle:
            handle["matrix/data"][:] = np.log1p(handle["matrix/data"][:])
        with pytest.raises(IngestError, match="non-integer"):
            read_gse178341(path)

    def test_unknown_patient_raises(self, deposit):
        pytest.importorskip("anndata")
        path, _ = deposit
        with pytest.raises(IngestError, match="not in the file"):
            read_gse178341(path, patients=["C1", "NOPE"])
