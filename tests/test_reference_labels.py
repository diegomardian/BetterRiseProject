"""Differentiation labels. W1.

Two families of test matter most:

- **Invariant 2.** No marker set that touches a label may contain a panel gene,
  including W1's own rung markers. A leaked label breaks the study exactly as a
  leaked reference does, and `build_signature()` cannot see labels.
- **Non-overwriting columns.** Eight label columns must coexist. The whole
  granularity-and-axis analysis (§6.2) is impossible if one clobbers another.

The synthetic cohort has a real gradient in it: stem-like cells high in the axis-1
markers, mature cells low, goblet cells high in the axis-2 markers, plus a BEST4+
population and non-epithelial cells. So the labels have something true to recover
rather than only shapes to check.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.common.panel import axis_genes, panel_genes
from src.reference.labels import (
    BEST4_MARKERS,
    NON_EPITHELIAL,
    RUNG_SPECS,
    TRANSCRIPT_AXES,
    LabelError,
    assign_labels,
    cell_type_vector,
    describe_labels,
    label_column,
    label_columns,
    mature_cell_counts,
    mature_mask,
    maturity_score,
    maturity_summary,
    score_markers,
)
from src.reference.signature import LeakageError

RNG = np.random.default_rng(20260817)

#: Target set used by most tests. Tier A (GUCA2A, GUCA2B, OTOP2, CA7) does not
#: collide with either axis, so both axes are usable — which is the ordinary case.
#: Runs testing MUC2 or TFF3 are the exception; see TestNoPanelLeakage.
TARGETS = ["GUCA2A", "GUCA2B", "OTOP2", "CA7"]

STEM = list(axis_genes("stem_pole"))              # LGR5, ASCL2, MKI67, OLFM4, SMOC2
GOBLET = list(axis_genes("opposite_lineage"))     # MUC2, TFF3, SPDEF, ITLN1
BEST4 = sorted(BEST4_MARKERS)
FILLER = [f"FILL{i:03d}" for i in range(40)]
GENES = STEM + GOBLET + BEST4 + FILLER

#: cells per population, per sample
POPULATIONS = {"stem": 60, "mature": 60, "goblet": 30, "best4": 10, "immune": 40}
SAMPLES = ("s1", "s2")


def _cohort():
    """Two samples, each with a real differentiation gradient."""
    blocks, compartment, sample_id, truth = [], [], [], []
    for sample_index, sample in enumerate(SAMPLES):
        depth = 5_000 * (1 + sample_index)          # samples differ in depth
        for population, n in POPULATIONS.items():
            profile = np.full(len(GENES), 1.0)
            if population == "stem":
                for gene in STEM:
                    profile[GENES.index(gene)] = 60.0
            elif population == "goblet":
                for gene in GOBLET:
                    profile[GENES.index(gene)] = 60.0
            elif population == "best4":
                for gene in BEST4:
                    profile[GENES.index(gene)] = 80.0
            # "mature" and "immune" stay at baseline for the marker sets
            profile = profile / profile.sum()
            counts = RNG.poisson(np.outer(np.full(n, depth), profile))
            blocks.append(counts.astype(np.int64))
            compartment += ["immune" if population == "immune" else "epithelial"] * n
            sample_id += [sample] * n
            truth += [population] * n
    return (
        np.vstack(blocks),
        np.array(compartment, dtype=object),
        np.array(sample_id, dtype=object),
        np.array(truth, dtype=object),
    )


@pytest.fixture(scope="module")
def cohort():
    return _cohort()


# ---------------------------------------------------------------------------
# Invariant 2 — the one that matters
# ---------------------------------------------------------------------------


class TestNoPanelLeakage:
    def test_a_panel_gene_in_a_marker_set_raises(self, cohort):
        matrix, *_ = cohort
        with pytest.raises(LeakageError, match="GUCA2A"):
            score_markers(matrix, GENES, ["LGR5", "GUCA2A"], context="test labels",
                          target_genes=TARGETS)

    def test_the_best4_rung_markers_are_panel_clean(self):
        """BEST4 and SPIB are absent from the panel on purpose (§3.2's sequencing
        constraint) so they remain available as labels. Verify it stayed true."""
        assert not (set(panel_genes()) & set(BEST4_MARKERS))

    def test_every_rung_spec_marker_set_is_panel_clean(self):
        targets = set(panel_genes())
        for rung, spec in RUNG_SPECS.items():
            if spec.markers:
                assert not (targets & set(spec.markers)), f"{rung} leaks"

    def test_axis_one_is_panel_clean(self):
        assert not (set(panel_genes()) & set(axis_genes("stem_pole")))

    def test_axis_two_collides_with_the_panel_as_documented(self):
        """MUC2 and TFF3 are tier E AND axis 2 — open decision #1.

        The narrow reading is implemented: build_signature() excludes the target
        set for the run in question, so axis 2 is usable EXCEPT for runs testing
        MUC2 or TFF3. This test pins that the collision is still exactly those
        two, so a future edit that adds a third fails loudly.
        """
        collision = set(panel_genes()) & set(axis_genes("opposite_lineage"))
        assert collision == {"MUC2", "TFF3"}

    def test_axis_two_scores_when_the_target_is_not_muc2_or_tff2(self, cohort):
        """A run testing GUCA2A may use axis 2; only MUC2/TFF3 runs may not."""
        matrix, *_ = cohort
        scores = score_markers(
            matrix, GENES, axis_genes("opposite_lineage"),
            context="axis 2", target_genes=["GUCA2A"],
        )
        assert scores.shape[0] == matrix.shape[0]

    def test_axis_two_refuses_when_muc2_is_the_target(self, cohort):
        matrix, *_ = cohort
        with pytest.raises(LeakageError, match="MUC2"):
            score_markers(
                matrix, GENES, axis_genes("opposite_lineage"),
                context="axis 2", target_genes=["MUC2"],
            )


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


class TestScoring:
    def test_stem_cells_score_least_mature_on_axis_one(self, cohort):
        matrix, _, _, truth = cohort
        scores = maturity_score(matrix, GENES, "stem_pole", target_genes=TARGETS)
        assert scores[truth == "stem"].mean() < scores[truth == "mature"].mean()

    def test_goblet_cells_score_least_mature_on_axis_two(self, cohort):
        matrix, _, _, truth = cohort
        scores = maturity_score(matrix, GENES, "opposite_lineage", target_genes=TARGETS)
        assert scores[truth == "goblet"].mean() < scores[truth == "mature"].mean()

    def test_the_two_axes_are_not_the_same_ordering(self, cohort):
        """Structurally different axes. If they agreed perfectly, the
        agreement-across-axes argument would be vacuous."""
        matrix, _, _, _ = cohort
        a = maturity_score(matrix, GENES, "stem_pole", target_genes=TARGETS)
        b = maturity_score(matrix, GENES, "opposite_lineage", target_genes=TARGETS)
        assert abs(np.corrcoef(a, b)[0, 1]) < 0.95

    def test_depth_normalisation_is_applied(self, cohort):
        """Chemistry differs across samples, so a raw-count score would rank a
        deep cell above a shallow one for technical reasons alone."""
        matrix, _, sample_id, truth = cohort
        scores = maturity_score(matrix, GENES, "stem_pole", target_genes=TARGETS, normalise=True)
        mature = truth == "mature"
        by_sample = [scores[mature & (sample_id == s)].mean() for s in SAMPLES]
        assert abs(by_sample[0] - by_sample[1]) < 0.5

    def test_absent_markers_raise_with_a_naming_hint(self, cohort):
        matrix, *_ = cohort
        with pytest.raises(LabelError, match="open decision #3"):
            score_markers(matrix, GENES, ["ENSG00000139292"], context="test",
                          target_genes=TARGETS)

    def test_partially_present_marker_set_still_scores(self, cohort):
        matrix, *_ = cohort
        scores = score_markers(matrix, GENES, ["LGR5", "NOT_A_REAL_GENE"], context="test",
                               target_genes=TARGETS)
        assert np.isfinite(scores).all()

    def test_non_transcript_axis_is_refused(self, cohort):
        matrix, *_ = cohort
        for axis in ("chromatin", "spatial"):
            with pytest.raises(LabelError, match="not a transcript-based axis"):
                maturity_score(matrix, GENES, axis, target_genes=TARGETS)


# ---------------------------------------------------------------------------
# The label grid
# ---------------------------------------------------------------------------


class TestLabelGrid:
    def test_eight_columns_all_coexist(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=sample_id, target_genes=TARGETS
        )
        assert list(labels.columns) == label_columns()
        assert len(labels.columns) == 8
        assert not labels.isna().any().any()

    def test_no_column_overwrites_another(self, cohort):
        """§4's 'labels stored as separate columns, never overwriting each
        other'. Different axes must give different label vectors."""
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        a = labels[label_column("stem_pole", "crypt_position")].astype(str)
        b = labels[label_column("opposite_lineage", "crypt_position")].astype(str)
        assert not a.equals(b)

    def test_row_order_and_length_are_preserved(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        assert len(labels) == matrix.shape[0]

    def test_index_is_carried_through_for_anndata_alignment(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        barcodes = [f"cell{i}" for i in range(matrix.shape[0])]
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=sample_id,
            target_genes=TARGETS, index=barcodes,
        )
        assert list(labels.index) == barcodes

    def test_non_epithelial_cells_are_never_given_a_maturity_call(self, cohort):
        matrix, compartment, sample_id, truth = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        immune = truth == "immune"
        for column in labels.columns:
            assert (labels.loc[immune, column].astype(str) == NON_EPITHELIAL).all()
            assert not (labels.loc[~immune, column].astype(str) == NON_EPITHELIAL).any()

    def test_rung_bin_counts_follow_the_specs(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        for rung, spec in RUNG_SPECS.items():
            values = set(labels[label_column("stem_pole", rung)].astype(str)) - {NON_EPITHELIAL}
            assert values <= set(spec.bins), rung

    def test_granularity_curve_is_monotone_in_mature_fraction(self, cohort):
        """The point of the rungs. The coarsest calls all epithelium mature; the
        finest calls ~5% mature. If this ordering broke, the granularity curve
        would be meaningless."""
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        epithelial = compartment == "epithelial"
        fractions = {
            rung: mature_mask(labels, "stem_pole", rung)[epithelial].mean()
            for rung in RUNG_SPECS
        }
        assert fractions["epithelial"] == pytest.approx(1.0)
        assert fractions["epithelial"] > fractions["lineage"] > fractions["best4"]
        assert fractions["best4"] < 0.15

    def test_bins_are_computed_within_sample_not_pooled(self, cohort):
        """Samples differ in depth here. Pooled quantiles would let the deeper
        sample's library size decide the shallower one's labels."""
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        mature = mature_mask(labels, "stem_pole", "crypt_position")
        epithelial = compartment == "epithelial"
        shares = [
            mature[epithelial & (sample_id == s)].mean() for s in SAMPLES
        ]
        assert abs(shares[0] - shares[1]) < 0.15

    def test_unknown_rung_raises(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        with pytest.raises(LabelError, match="no RungSpec"):
            assign_labels(
                matrix, GENES, compartment=compartment, sample_id=sample_id,
                target_genes=TARGETS, rungs=["epithelial", "villus_tip"],
            )

    def test_length_mismatch_raises(self, cohort):
        matrix, compartment, _, _ = cohort
        with pytest.raises(LabelError, match="entries for"):
            assign_labels(matrix, GENES, compartment=compartment, sample_id=["s1"],
                          target_genes=TARGETS)

    def test_no_epithelium_raises_with_the_expected_value(self, cohort):
        matrix, _, sample_id, _ = cohort
        with pytest.raises(LabelError, match="assign_compartments"):
            assign_labels(
                matrix, GENES,
                compartment=["immune"] * matrix.shape[0], sample_id=sample_id,
                target_genes=TARGETS,
            )


# ---------------------------------------------------------------------------
# What W2 and W4 consume
# ---------------------------------------------------------------------------


class TestDownstreamInterface:
    def test_mature_mask_selects_the_most_mature_bin(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        mask = mature_mask(labels, "stem_pole", "lineage")
        values = labels.loc[mask, label_column("stem_pole", "lineage")].astype(str)
        assert set(values) == {RUNG_SPECS["lineage"].mature}

    def test_mature_mask_rejects_a_missing_column(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=sample_id,
            target_genes=TARGETS, axes=["stem_pole"], rungs=["lineage"],
        )
        with pytest.raises(LabelError, match="not in labels"):
            mature_mask(labels, "opposite_lineage", "lineage")

    def test_counts_have_one_row_per_patient_tissue_axis_rung(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        patient = np.where(sample_id == "s1", "P1", "P2")
        tissue = np.where(sample_id == "s1", "tumour", "normal")
        counts = mature_cell_counts(labels, patient_id=patient, tissue=tissue)

        assert set(counts["labeling_axis"]) == set(TRANSCRIPT_AXES)
        assert set(counts["granularity_rung"]) == set(RUNG_SPECS)
        assert len(counts) == 2 * len(TRANSCRIPT_AXES) * len(RUNG_SPECS)

    def test_counts_expose_n_cells_mature_for_the_schema(self, cohort):
        """n_cells_mature is a frozen schema column and drives positivity."""
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        patient = np.full(matrix.shape[0], "P1")
        counts = mature_cell_counts(labels, patient_id=patient, tissue=sample_id)
        assert {"n_cells_mature", "n_cells_epithelial", "mature_fraction"} <= set(counts.columns)
        assert (counts["n_cells_mature"] <= counts["n_cells_epithelial"]).all()
        assert counts["mature_fraction"].between(0, 1).all()

    def test_counts_exclude_non_epithelial_cells_from_the_denominator(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        counts = mature_cell_counts(
            labels, patient_id=np.full(matrix.shape[0], "P1"), tissue=sample_id
        )
        n_epithelial = int((compartment == "epithelial").sum())
        for rung in RUNG_SPECS:
            subset = counts[counts["granularity_rung"] == rung]
            assert int(subset["n_cells_epithelial"].sum()) / len(TRANSCRIPT_AXES) == n_epithelial

    def test_describe_labels_is_a_readable_summary(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        summary = describe_labels(labels)
        assert {"column", "label", "n_cells"} == set(summary.columns)
        assert summary["n_cells"].sum() == len(labels) * len(labels.columns)


def test_every_rung_spec_records_its_rationale():
    """The rung partitions are modelling choices, not measurements. Each must say
    why it is drawn where it is, because changing it changes the split."""
    for rung, spec in RUNG_SPECS.items():
        assert len(spec.rationale) > 80, rung
        assert spec.mature == spec.bins[-1]


def test_rung_names_match_the_frozen_config():
    from src.common.panel import granularity_rungs

    assert set(RUNG_SPECS) == set(granularity_rungs())


def test_sparse_and_dense_agree(cohort):
    from scipy import sparse

    matrix, compartment, sample_id, _ = cohort
    dense = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
    spars = assign_labels(
        sparse.csr_matrix(matrix), GENES, compartment=compartment, sample_id=sample_id,
        target_genes=TARGETS,
    )
    pd.testing.assert_frame_equal(dense, spars)


class TestTargetGenesIsRequired:
    """Invariant 2 must be impossible to skip by omission.

    `build_signature()` refuses an empty target set for this reason; label
    construction now does the same. A default would bury open decision #1: the
    whole panel makes axis 2 unusable, a permissive default silently disables the
    invariant.
    """

    def test_omitting_target_genes_is_a_type_error(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        with pytest.raises(TypeError):
            assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id)

    @pytest.mark.parametrize("empty", [[], (), None, set()])
    def test_empty_target_genes_is_refused(self, cohort, empty):
        matrix, compartment, sample_id, _ = cohort
        with pytest.raises(LabelError, match="silently disables"):
            assign_labels(
                matrix, GENES, compartment=compartment, sample_id=sample_id,
                target_genes=empty,
            )

    def test_whole_panel_as_target_makes_axis_two_unusable(self, cohort):
        """Documents decision #1's consequence rather than hiding it. If the team
        adopts option (b) and drops MUC2/TFF3 from tier E, this stops raising —
        update the test then, do not delete it."""
        matrix, compartment, sample_id, _ = cohort
        with pytest.raises(LeakageError, match="MUC2"):
            assign_labels(
                matrix, GENES, compartment=compartment, sample_id=sample_id,
                target_genes=panel_genes(),
            )

    def test_a_muc2_run_uses_axis_one_only(self, cohort):
        """The documented workaround: pass axes=['stem_pole'] for MUC2/TFF3 runs."""
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=sample_id,
            target_genes=["MUC2"], axes=["stem_pole"],
        )
        assert list(labels.columns) == label_columns(axes=["stem_pole"])
        assert len(labels.columns) == 4

    def test_whole_panel_works_for_axis_one(self, cohort):
        """Axis 1 is panel-clean, so even the strictest target set is fine there."""
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=sample_id,
            target_genes=panel_genes(), axes=["stem_pole"],
        )
        assert not labels.isna().any().any()


class TestW2Interface:
    """The pseudobulk generator identifies mature cells by string equality."""

    def test_cell_type_vector_renames_the_mature_bin(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment,
                               sample_id=sample_id, target_genes=TARGETS)
        vector = cell_type_vector(labels, "stem_pole", "best4")
        assert "mature_colonocyte" in set(vector)
        assert RUNG_SPECS["best4"].mature not in set(vector)
        assert NON_EPITHELIAL in set(vector)

    def test_it_matches_the_mature_mask(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment,
                               sample_id=sample_id, target_genes=TARGETS)
        for rung in RUNG_SPECS:
            vector = cell_type_vector(labels, "stem_pole", rung)
            np.testing.assert_array_equal(
                vector == "mature_colonocyte", mature_mask(labels, "stem_pole", rung)
            )

    def test_generate_pseudobulk_accepts_it(self, cohort):
        """End-to-end against W2's actual function, not a mock of it."""
        from src.harness.pseudobulk import generate_pseudobulk

        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment,
                               sample_id=sample_id, target_genes=TARGETS)
        cell_type = cell_type_vector(labels, "stem_pole", "lineage")
        other = RUNG_SPECS["lineage"].bins[0]

        sample = generate_pseudobulk(
            matrix, cell_type, sample_id, GENES,
            composition_normal={"mature_colonocyte": 0.6, other: 0.4},
            composition_tumour={"mature_colonocyte": 0.2, other: 0.8},
            shift={"FILL000": 0.5},
            held_out_patients=["s1"],
            n_cells=200,
            seed=1,
        )
        assert sample is not None

    def test_custom_mature_label_is_honoured(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment,
                               sample_id=sample_id, target_genes=TARGETS)
        vector = cell_type_vector(labels, "stem_pole", "best4", mature_label="MATURE")
        assert "MATURE" in set(vector)

    def test_unknown_column_raises(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment,
                               sample_id=sample_id, target_genes=TARGETS,
                               axes=["stem_pole"])
        with pytest.raises(LabelError, match="not in labels"):
            cell_type_vector(labels, "opposite_lineage", "best4")


class TestW4Interface:
    """maturity_summary must drop into decompose_cohort's required columns."""

    def _summary(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment,
                               sample_id=sample_id, target_genes=TARGETS)
        patient = np.full(matrix.shape[0], "P1")
        tissue = np.where(sample_id == "s1", "tumour", "normal")
        return maturity_summary(
            labels, patient_id=patient, tissue=tissue, study_id="GSE178341"
        )

    def test_supplies_every_label_derived_column(self, cohort):
        summary = self._summary(cohort)
        assert {
            "patient_id", "study_id", "granularity_rung", "labeling_axis",
            "frac_mature_normal", "frac_mature_tumour", "n_cells_mature",
        } <= set(summary.columns)

    def test_only_gene_level_means_are_left_to_the_caller(self, cohort):
        """Everything else decompose_cohort needs is present."""
        from src.estimator.kitagawa import _SUMMARY_COLUMNS

        summary = self._summary(cohort)
        outstanding = set(_SUMMARY_COLUMNS) - set(summary.columns)
        assert outstanding == {"gene", "mean_normal", "mean_tumour"}

    def test_decompose_cohort_runs_once_means_are_attached(self, cohort):
        """End-to-end against W4's actual estimator.

        Note the fan-out: each input row becomes three output rows, one per
        weighting (normal, tumour, doubly_robust). §4's estimator definition
        requires all three be reported, never folded together.
        """
        from src.estimator.kitagawa import decompose_cohort

        summary = self._summary(cohort).assign(
            gene="GUCA2A", mean_normal=2.0, mean_tumour=0.4
        )
        out = decompose_cohort(summary)
        assert len(out) == len(summary) * 3
        assert set(out["weighting"]) == {"normal", "tumour", "doubly_robust"}
        assert {"compositional", "intrinsic", "estimability"} <= set(out.columns)

    def test_the_invariant_one_path_survives_the_round_trip(self, cohort):
        """Rows the estimator cannot estimate must carry None, not 0.0."""
        from src.estimator.kitagawa import decompose_cohort

        summary = self._summary(cohort).assign(
            gene="GUCA2A", mean_normal=2.0, mean_tumour=0.4
        )
        out = decompose_cohort(summary)
        unestimable = out[out["estimability"] == "not_estimable"]
        if len(unestimable):
            assert unestimable["intrinsic"].isna().all()
            assert not (unestimable["intrinsic"] == 0.0).any()

    def test_fractions_are_proportions(self, cohort):
        summary = self._summary(cohort)
        for column in ("frac_mature_normal", "frac_mature_tumour"):
            assert summary[column].between(0, 1).all()

    def test_n_cells_mature_is_the_tumour_arm(self, cohort):
        """Positivity is about mature cells surviving in the TUMOUR — the arm
        that can run out. The intrinsic term is tumour-mature-fraction weighted."""
        summary = self._summary(cohort)
        assert (summary["n_cells_mature"] == summary["n_cells_mature_tumour"]).all()
