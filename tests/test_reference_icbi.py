"""Reading the ICBI atlas metadata without downloading the atlas. W1.

§8.2 says to take the metadata table before committing compute. No such table is
published — the smallest published artifact is a 32.7 GB h5ad — so this reads
/obs over HTTP range requests. The tests build a small local h5ad in the same
shape and drive the same code path against it, so the decoding and the summaries
are covered without touching the network.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd
import pytest

from src.reference.icbi import (
    NORMAL_SAMPLE_TYPES,
    PLATE_PLATFORMS,
    TUMOUR_SAMPLE_TYPES,
    ICBIError,
    _normalise_label,
    depth_by_platform,
    enrichment_summary,
    epithelial_fraction,
    paired_sample_summary,
    platform_summary,
)

h5py = pytest.importorskip("h5py")

# THE ATLAS'S OWN VOCABULARY, not ours. This fixture used to say "smartseq2"
# and "tumor"/"normal" -- the strings the implementation was testing for -- so
# fixture and code agreed with each other and both disagreed with the data. The
# summaries returned 0 plate-based cells and 0 paired patients on the real
# atlas, against 24,136 and 229, and every test here passed.
# `platform` is the coarse column ("10x 3p"); "10x 3' v2" lives in
# `platform_fine`. The fixture mixed the two, which the vocabulary guard below
# caught on its first run -- the same defect a third time.
PLATFORMS = ["10x 3p", "10x 5p", "Smart-seq2"]
DEPTHS = {"10x 3p": 2_500, "10x 5p": 1_800, "Smart-seq2": 9_000}
TUMOUR, NORMAL = "primary tumor", "adjacent normal"


def _atlas(tmp_path, n_per_platform=60):
    """An h5ad shaped like the atlas: categorical obs columns, numeric depth."""
    rows = []
    for platform in PLATFORMS:
        for index in range(n_per_platform):
            patient = f"P{index % 6}"
            # Arm must vary WITHIN a patient, so it keys on index // 6 rather
            # than index % 2 — the latter is locked to patient identity, which
            # made every patient single-armed and nothing paired.
            # P0-P3 get both arms; P4-P5 are tumour-only.
            sample_type = (
                TUMOUR if index % 6 >= 4
                else (NORMAL if (index // 6) % 2 == 0 else TUMOUR)
            )
            rows.append(
                {
                    "study_id": f"study_{platform.replace(' ', '_')}",
                    "sample_id": f"{patient}_{sample_type}_{platform}",
                    "patient_id": f"{platform}_{patient}",
                    "sample_type": sample_type,
                    "platform": platform,
                    "enrichment_cell_types": "naive" if index % 3 else "CD45+",
                    "atlas_cell_type_coarse": "Epithelial cell" if index % 2 else "T cell",
                    "n_genes": float(DEPTHS[platform] + index),
                    "total_counts": float(DEPTHS[platform] * 4 + index),
                }
            )
    frame = pd.DataFrame(rows)

    path = tmp_path / "atlas.h5ad"
    with h5py.File(path, "w") as h5:
        obs = h5.create_group("obs")
        for column in frame.columns:
            values = frame[column]
            if values.dtype == object:
                categories = sorted(values.unique())
                codes = values.map({c: i for i, c in enumerate(categories)}).to_numpy()
                group = obs.create_group(column)
                group.create_dataset("codes", data=codes.astype(np.int32))
                group.create_dataset(
                    "categories", data=np.array(categories, dtype="S64")
                )
            else:
                obs.create_dataset(column, data=values.to_numpy())
        h5.create_group("var")
    return path, frame


@pytest.fixture
def atlas(tmp_path):
    return _atlas(tmp_path)


def _read_local(path, columns):
    """Drive the same decoding as read_atlas_obs, from a local file."""
    from src.reference.icbi import _decode

    with h5py.File(path, "r") as h5:
        available = list(h5["obs"].keys())
        wanted = [c for c in columns if c in available]
        return pd.DataFrame({name: _decode(h5["obs"], name) for name in wanted})


class TestDecoding:
    def test_categorical_columns_are_resolved_to_labels(self, atlas):
        """anndata stores categoricals as codes plus categories. Reading the
        codes and forgetting to map them would silently give integers."""
        path, frame = atlas
        obs = _read_local(path, ["platform", "sample_type"])
        assert set(obs["platform"]) == set(PLATFORMS)
        assert set(obs["sample_type"]) <= {TUMOUR, NORMAL}

    def test_numeric_columns_survive(self, atlas):
        path, frame = atlas
        obs = _read_local(path, ["n_genes"])
        assert obs["n_genes"].dtype.kind == "f"
        assert obs["n_genes"].min() >= min(DEPTHS.values())

    def test_missing_columns_are_skipped_not_fatal(self, atlas):
        path, _ = atlas
        obs = _read_local(path, ["platform", "not_a_real_column"])
        assert list(obs.columns) == ["platform"]


class TestSummaries:
    def _obs(self, atlas):
        path, _ = atlas
        return _read_local(
            path,
            ["study_id", "sample_id", "patient_id", "sample_type", "platform",
             "enrichment_cell_types", "atlas_cell_type_coarse", "n_genes",
             "total_counts"],
        )

    def test_platform_summary_flags_the_plate_subset(self, atlas):
        """The row that decides whether axis 1 is measurable anywhere."""
        out = platform_summary(self._obs(atlas))
        plate = out[out["plate_based"]]
        assert len(plate) == 1
        assert _normalise_label(plate.iloc[0]["platform"]) in PLATE_PLATFORMS
        assert int(plate.iloc[0]["n_cells"]) == 60

    def test_depth_by_platform_ranks_plate_highest(self, atlas):
        """Measured, not assumed from the protocol's reputation."""
        out = depth_by_platform(self._obs(atlas))
        assert out.iloc[0]["platform"] == "Smart-seq2"
        assert bool(out.iloc[0]["plate_based"])
        assert out.iloc[0]["median_genes"] > out.iloc[-1]["median_genes"]

    def test_depth_reports_the_share_clearing_five_thousand_genes(self, atlas):
        """Roughly where five sparse markers stop dropping out together."""
        out = depth_by_platform(self._obs(atlas)).set_index("platform")
        assert out.loc["Smart-seq2", "share_over_5k_genes"] == 1.0
        assert out.loc["10x 5p", "share_over_5k_genes"] == 0.0

    def test_paired_summary_counts_patients_with_both_arms(self, atlas):
        out = paired_sample_summary(self._obs(atlas))
        assert {"n_patients", "n_paired"} <= set(out.columns)
        assert (out["n_paired"] <= out["n_patients"]).all()
        # P0-P3 have both arms in every platform's study; P4-P5 are tumour-only.
        assert set(out["n_paired"]) == {4}
        assert set(out["n_patients"]) == {6}

    def test_enrichment_summary_separates_sorted_from_naive(self, atlas):
        out = enrichment_summary(self._obs(atlas)).set_index("enrichment_cell_types")
        assert "naive" in out.index and "CD45+" in out.index

    def test_epithelial_fraction_is_a_proportion(self, atlas):
        out = epithelial_fraction(self._obs(atlas))
        assert out["epithelial_fraction"].between(0, 1).all()
        assert out["epithelial_fraction"].max() == pytest.approx(0.5, abs=0.05)

    def test_missing_column_raises_rather_than_returning_nonsense(self, atlas):
        obs = self._obs(atlas).drop(columns=["platform"])
        with pytest.raises(ICBIError, match="platform"):
            platform_summary(obs)


def test_read_atlas_obs_rejects_a_server_without_range_support(monkeypatch):
    """Without byte ranges the whole 32.7 GB object would come down. Refuse."""
    import src.reference.icbi as icbi

    class Response:
        headers: dict = {}

        def raise_for_status(self):
            return None

    class Session:
        def head(self, *args, **kwargs):
            return Response()

    monkeypatch.setattr(
        icbi, "HTTPRangeFile", icbi.HTTPRangeFile
    )
    with pytest.raises(ICBIError, match="byte ranges"):
        icbi.HTTPRangeFile("https://example.invalid/x.h5ad", session=Session())


# ---------------------------------------------------------------------------
# The guard that would have caught both defects


class TestVocabulary:
    """Both bugs were one bug: a label the atlas does not use.

    `PLATE_PLATFORMS` held "smartseq2" against an atlas that writes
    "Smart-seq2", and `paired_sample_summary` tested for "tumor"/"normal"
    against "primary tumor"/"adjacent normal". Neither matched anything, both
    aggregated to zero, and zero was printed as a finding: "plate-based cells
    available: 0", and 0 paired patients in all 49 studies.

    The tests could not catch it because the FIXTURE was written from the code's
    vocabulary rather than the atlas's, so it agreed with the bug.
    """

    def test_spelling_variants_all_resolve_to_the_same_platform(self):
        for spelling in ("Smart-seq2", "smartseq2", "SMART-Seq2", "smart seq2"):
            assert _normalise_label(spelling) in PLATE_PLATFORMS, spelling

    def test_a_vocabulary_that_matches_nothing_raises_instead_of_returning_zero(self):
        """The input that made the old code return 'no paired patients'."""
        from src.reference.icbi import paired_sample_summary

        obs = pd.DataFrame({
            "study_id": ["s"] * 4,
            "patient_id": ["p1", "p1", "p2", "p2"],
            "sample_type": ["carcinoma", "healthy normal", "carcinoma", "healthy normal"],
        })
        with pytest.raises(ICBIError, match="Refusing rather than returning zeros"):
            paired_sample_summary(obs)

    def test_healthy_normal_is_not_a_matched_normal(self):
        """A different donor is not the same patient's adjacent tissue."""
        assert "healthy normal" not in NORMAL_SAMPLE_TYPES
        assert "metastasis" not in TUMOUR_SAMPLE_TYPES

    def test_the_fixture_speaks_the_real_atlas_vocabulary(self, atlas):
        """Pins the fixture to the cached obs, if one is present.

        This is the check that was missing. A fixture built from the
        implementation can only confirm the implementation; it has to be built
        from the data.
        """
        cache = pathlib.Path("data/interim/icbi_obs.parquet")
        if not cache.exists():
            pytest.skip("no cached atlas obs to compare against")
        real = pd.read_parquet(cache, columns=["sample_type", "platform"])
        _path, frame = atlas
        assert set(frame["sample_type"]) <= set(real["sample_type"].unique()), (
            "the fixture uses sample_type labels the atlas does not"
        )
        assert set(frame["platform"]) <= set(real["platform"].unique()), (
            "the fixture uses platform labels the atlas does not"
        )

    def test_the_real_atlas_has_paired_patients_and_plate_cells(self):
        """Against the cached obs: 229 and 24,136, not 0 and 0."""
        from src.reference.icbi import paired_sample_summary, platform_summary

        cache = pathlib.Path("data/interim/icbi_obs.parquet")
        if not cache.exists():
            pytest.skip("no cached atlas obs")
        obs = pd.read_parquet(cache)
        paired = paired_sample_summary(obs)
        assert int(paired["n_paired"].sum()) == 229
        assert int((paired["n_paired"] > 0).sum()) == 24

        platforms = platform_summary(obs)
        plate = platforms[platforms["plate_based"]]
        assert int(plate["n_cells"].sum()) == 24136
        assert set(plate["platform"]) == {"Smart-seq2", "scTrio-seq2", "SMARTer (C1)"}
