"""Slicing the atlas: the row extractor, and the layer it must not read.

The extractor is worth testing hard because its failure mode is quiet -- a
wrong row order gives every cell someone else's expression, and reading /X
instead of the counts layer gives log1p values to a detection statistic that
reports them without complaint.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.icbi_slice import (
    COUNTS_LAYER,
    SliceError,
    arms,
    assert_raw_counts,
    compartments,
    read_cells,
    read_var,
)

h5py = pytest.importorskip("h5py")
sp = pytest.importorskip("scipy.sparse")


@pytest.fixture
def atlas(tmp_path):
    """A miniature h5ad shaped like the real one: CSR X (log1p) and counts."""
    rng = np.random.default_rng(3)
    n_cells, n_genes = 60, 40
    dense = rng.poisson(2.0, size=(n_cells, n_genes)).astype(np.float32)
    counts = sp.csr_matrix(dense)
    logged = counts.copy()
    logged.data = np.log1p(logged.data)

    path = tmp_path / "atlas.h5ad"
    with h5py.File(path, "w") as h5:
        for name, matrix in (("X", logged), ("layers/counts", counts)):
            group = h5.create_group(name)
            group.attrs["encoding-type"] = "csr_matrix"
            group.attrs["shape"] = np.array([n_cells, n_genes])
            group.create_dataset("data", data=matrix.data)
            group.create_dataset("indices", data=matrix.indices.astype(np.int64))
            group.create_dataset("indptr", data=matrix.indptr.astype(np.int64))
        var = h5.create_group("var")
        var.attrs["_index"] = "_index"
        for key, values in (
            ("_index", [f"ENSG{i:08d}" for i in range(n_genes)]),
            ("GeneSymbol", [f"SYM{i}" for i in range(n_genes)]),
        ):
            var.create_dataset(key, data=np.array(values, dtype="S32"))
    return path, dense


def test_rows_come_back_in_the_order_asked_for(atlas):
    """The quiet failure: HDF5 needs increasing indices, so the extractor sorts.
    Forgetting to permute back gives every cell someone else's expression."""
    path, dense = atlas
    wanted = np.array([17, 3, 41, 3, 58, 0])
    got = read_cells(path, wanted).toarray()
    assert got.shape == (len(wanted), dense.shape[1])
    for j, row in enumerate(wanted):
        np.testing.assert_array_equal(got[j], dense[row])


def test_a_single_row_and_a_full_pass_both_work(atlas):
    path, dense = atlas
    np.testing.assert_array_equal(read_cells(path, [7]).toarray()[0], dense[7])
    every = read_cells(path, np.arange(len(dense))).toarray()
    np.testing.assert_array_equal(every, dense)


def test_the_counts_layer_is_raw_and_X_is_not(atlas):
    """The distinction the whole module exists to preserve."""
    path, _ = atlas
    assert_raw_counts(read_cells(path, [1, 2, 3]), context="counts layer")
    with pytest.raises(SliceError, match="not raw integer counts"):
        assert_raw_counts(read_cells(path, [1, 2, 3], layer="X"), context="/X")


def test_reading_a_missing_layer_refuses_rather_than_falling_back(atlas):
    path, _ = atlas
    with pytest.raises(SliceError, match="Do NOT fall back"):
        read_cells(path, [0], layer="layers/nope")


def test_out_of_range_rows_are_refused(atlas):
    path, dense = atlas
    with pytest.raises(SliceError, match="past the end"):
        read_cells(path, [len(dense)])


def test_no_rows_is_refused(atlas):
    path, _ = atlas
    with pytest.raises(SliceError, match="no rows"):
        read_cells(path, [])


def test_var_returns_both_identifier_spaces(atlas):
    path, _ = atlas
    var = read_var(path)
    assert list(var.columns) == ["ensembl_id", "gene_symbol"]
    assert var["ensembl_id"].iloc[0].startswith("ENSG")
    assert var["gene_symbol"].iloc[0] == "SYM0"


# ---------------------------------------------------------------------------
# The vocabulary maps


def test_cancer_cells_count_as_epithelium():
    """Excluding them takes most of the tumour arm with it: 7 patients vs 136."""
    got = compartments(pd.Series(["Cancer cell", "Epithelial cell", "T cell"]))
    assert list(got) == ["epithelial", "epithelial", "immune"]


def test_an_unmapped_cell_type_is_dropped_not_silently_binned(caplog):
    got = compartments(pd.Series(["Epithelial cell", "Tuft cell of the future"]))
    assert got.isna().sum() == 1
    assert "unmapped cell types" in caplog.text.lower()


def test_only_the_two_arms_map_and_healthy_normal_does_not():
    """A different donor is not the same patient's adjacent tissue."""
    got = arms(pd.Series(
        ["primary tumor", "adjacent normal", "healthy normal", "metastasis"]
    ))
    assert list(got[:2]) == ["tumour", "normal"]
    assert got[2:].isna().all()


def test_the_counts_layer_constant_is_not_X():
    assert COUNTS_LAYER == "layers/counts"


# ---------------------------------------------------------------------------
# The adenoma arms, and the label that cost a cohort


def test_healthy_normal_is_a_reference_arm_for_the_adenoma_reading():
    """The fifth vocabulary error in this repo, and the first inside a verdict.

    `healthy normal` was excluded because it can mean a different donor. But
    the reading only ever compares two arms OF ONE PATIENT and every caller
    groups by patient_id first, so a `healthy normal` reachable from a patient
    who also has a polyp is that patient's own mucosa. Reading the LABEL put
    the VUMC/HTAN cohort's usable pairs at 0; reading it per patient puts them
    at 44.
    """
    from src.reference.icbi_slice import ADENOMA_TISSUE_MAP

    got = arms(pd.Series(
        ["polyp", "healthy normal", "adjacent normal", "primary tumor"]
    ), ADENOMA_TISSUE_MAP)
    assert list(got[:3]) == ["tumour", "normal", "normal"]
    assert pd.isna(got.iloc[3]), "carcinoma is not the adenoma reading's arm"


def test_the_carcinoma_mapping_is_unchanged():
    """The adenoma fix must not move the published 13-study result."""
    got = arms(pd.Series(["primary tumor", "adjacent normal", "healthy normal", "polyp"]))
    assert list(got[:2]) == ["tumour", "normal"]
    assert got[2:].isna().all(), (
        "healthy normal or polyp entered the CARCINOMA contrast; the committed "
        "13-study result was computed without them"
    )


def test_allowing_healthy_normal_adds_no_carcinoma_patients():
    """Asserted against the real obs, because it is what makes the fix safe.

    If any of the 14 candidate studies gained a patient from the looser
    reference set, the committed carcinoma result would be incomplete rather
    than merely conservative.
    """

    from src.common.paths import INTERIM_DIR, RESULTS_DIR

    cache = INTERIM_DIR / "icbi_obs.parquet"
    candidates = sorted(RESULTS_DIR.glob("*/icbi_premise_candidate_studies.parquet"))
    if not cache.exists() or not candidates:
        pytest.skip("the cached obs or the candidate table is absent")

    obs = pd.read_parquet(cache)
    studies = set(pd.read_parquet(candidates[-1])["study_id"])
    epithelial = {"Epithelial cell", "Cancer cell"}
    gained = {}
    for study in studies:
        block = obs[(obs["study_id"] == study)
                    & (obs["enrichment_cell_types"].astype(str) == "naive")]
        cells = block[block["atlas_cell_type_coarse"].astype(str).isin(epithelial)]
        counts = cells.groupby(["patient_id", "sample_type"]).size().unstack(fill_value=0)
        for col in ("primary tumor", "adjacent normal", "healthy normal"):
            if col not in counts:
                counts[col] = 0
        strict = ((counts["primary tumor"] >= 100) & (counts["adjacent normal"] >= 100)).sum()
        loose = (
            (counts["primary tumor"] >= 100)
            & (counts[["adjacent normal", "healthy normal"]].max(axis=1) >= 100)
        ).sum()
        if loose != strict:
            gained[study] = int(loose - strict)
    assert not gained, (
        f"the looser reference set adds carcinoma patients {gained}; the "
        f"committed 13-study result would be incomplete, not conservative"
    )
