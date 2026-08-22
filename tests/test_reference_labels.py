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
    UNRESOLVED,
    LabelError,
    annotation_concordance,
    assign_labels,
    axis_tie_fraction,
    cell_type_vector,
    describe_labels,
    label_column,
    label_columns,
    label_depth_confounding,
    mature_cell_counts,
    mature_mask,
    maturity_score,
    maturity_summary,
    maturity_within_depth_strata,
    rung_degeneracy,
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
            values = set(labels[label_column("stem_pole", rung)].astype(str))
            values -= {NON_EPITHELIAL, UNRESOLVED}
            assert values <= set(spec.bins), rung

    def test_granularity_curve_is_monotone_in_mature_fraction(self, cohort):
        """The point of the rungs. The coarsest calls all epithelium mature; the
        finest calls ~5% mature. If this ordering broke, the granularity curve
        would be meaningless."""
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment, sample_id=sample_id,
                               target_genes=TARGETS)
        # Denominator is the RESOLVED epithelium — cells dropped by depth
        # matching are not immature, they are unmeasured.
        fractions = {}
        for rung in RUNG_SPECS:
            column = labels[label_column("stem_pole", rung)].astype(str).to_numpy()
            resolved = ~np.isin(column, [NON_EPITHELIAL, UNRESOLVED])
            fractions[rung] = mature_mask(labels, "stem_pole", rung)[resolved].mean()
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


class TestTheMeasurementTravelsWithTheNumber:
    """A mature fraction is not interpretable without what defined "mature".

    On axis 1 the mature bin is the block of cells with no marker detected at
    the depth target — a detection gate, not the median split the rung name
    implies. Two studies gating at different depths produce numbers that look
    comparable and are not, and invariant 4 has us meta-analyse across studies.
    So the depth and the definition ride on the row, like the git sha does.
    """

    def _summary(self, cohort, **kwargs):
        matrix, compartment, sample_id, _ = cohort
        patient = np.full(matrix.shape[0], "P1")
        tissue = np.where(sample_id == "s1", "tumour", "normal")
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=sample_id,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            **kwargs,
        )
        return maturity_summary(
            labels, patient_id=patient, tissue=tissue, study_id="GSE178341",
            depth_target=kwargs.get("depth_target"),
        )

    def test_the_depth_target_is_recorded(self, cohort):
        summary = self._summary(cohort, depth_target=2_000.0)
        assert (summary["depth_target"] == 2_000.0).all()
        assert (
            summary["mature_definition"]
            == "no axis marker detected at matched depth"
        ).all()

    def test_an_unmatched_run_says_so_rather_than_inventing_a_depth(self, cohort):
        summary = self._summary(cohort, depth_target=0)
        assert summary["depth_target"].isna().all() or (
            summary["depth_target"] == 0
        ).all()
        assert (summary["mature_definition"] != "").all()

    def test_extra_columns_do_not_break_W4(self, cohort):
        """decompose_cohort checks for missing columns, not for extra ones —
        but that is a promise worth pinning, since this adds three."""
        from src.estimator.kitagawa import decompose_cohort

        summary = self._summary(cohort, depth_target=2_000.0).assign(
            gene="GUCA2A", mean_normal=2.0, mean_tumour=0.4
        )
        out = decompose_cohort(summary)
        assert len(out) == len(summary) * 3


class TestDegeneracyIsCarriedOnTheRow:
    """The granularity curve must be buildable without recomputing degeneracy.

    A rung that collapsed onto a coarser rung's boundary is not a second point
    on the curve. The consumer building that curve is W4, which never sees the
    labels — only this summary — so the flag has to be on the row.
    """

    def _summary(self, **kwargs):
        matrix, compartment, patient, tissue = TestTieCollapse()._tied(**kwargs)
        # _tied is single-arm; give it a tumour arm so both fractions exist.
        tissue = np.array(
            ["normal"] * (len(tissue) // 2)
            + ["tumour"] * (len(tissue) - len(tissue) // 2),
            dtype=object,
        )
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["lineage", "crypt_position"],
        )
        return maturity_summary(
            labels, patient_id=patient, tissue=tissue, study_id="S",
            axes=["stem_pole"], rungs=["lineage", "crypt_position"],
        )

    def test_the_finer_rung_names_the_coarser_one_it_duplicates(self):
        summary = self._summary()
        finer = summary[summary["granularity_rung"] == "crypt_position"]
        assert len(finer)
        assert (finer["degenerate_with"] == "lineage").all()

    def test_the_coarser_rung_is_not_itself_flagged(self):
        """crypt_position duplicates lineage, not the other way round — the
        finer rung is the one claiming a resolution the data does not support."""
        summary = self._summary()
        coarser = summary[summary["granularity_rung"] == "lineage"]
        assert coarser["degenerate_with"].isna().all()

    def test_a_well_spread_axis_flags_nothing(self):
        summary = self._summary(n_tied=0, n_spread=900)
        assert summary["degenerate_with"].isna().all()


class TestZeroEpitheliumGroup:
    """A (patient, tissue) group with no epithelium must yield NaN, not a crash.

    Found on the real pilot: np.where evaluates both branches, so guarding the
    division did not prevent it and pandas raised ZeroDivisionError on the
    integer denominator.
    """

    def _labels_with_an_empty_group(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment,
                               sample_id=sample_id, target_genes=TARGETS)
        # P2/normal is made entirely non-epithelial.
        patient = np.where(compartment == "immune", "P2", "P1")
        tissue = np.where(compartment == "immune", "normal", "tumour")
        return labels, patient, tissue

    def test_counts_do_not_crash(self, cohort):
        labels, patient, tissue = self._labels_with_an_empty_group(cohort)
        counts = mature_cell_counts(labels, patient_id=patient, tissue=tissue)
        assert len(counts) > 0

    def test_the_empty_group_reports_nan_not_zero(self, cohort):
        """NaN, not 0.0 — the same distinction invariant 1 makes. A fraction of
        zero would claim 'no mature cells'; the truth is 'no cells to ask'."""
        labels, patient, tissue = self._labels_with_an_empty_group(cohort)
        counts = mature_cell_counts(labels, patient_id=patient, tissue=tissue)
        empty = counts[counts["n_cells_epithelial"] == 0]
        assert len(empty) > 0
        assert empty["mature_fraction"].isna().all()
        assert not (empty["mature_fraction"] == 0.0).any()

    def test_populated_groups_are_unaffected(self, cohort):
        labels, patient, tissue = self._labels_with_an_empty_group(cohort)
        counts = mature_cell_counts(labels, patient_id=patient, tissue=tissue)
        populated = counts[counts["n_cells_epithelial"] > 0]
        assert populated["mature_fraction"].notna().all()
        assert populated["mature_fraction"].between(0, 1).all()


class TestNonRangeIndex:
    """Labels carry barcodes as their index in real use. The counting code must
    not align on it.

    Missed on the first pass: the fixture never passed `index=`, so `labels` and
    the internal keys frame both had RangeIndexes and aligned by accident. On the
    real pilot, labels were barcode-indexed, the epithelial column aligned to
    nothing, and every n_cells_epithelial came back 0 with every mature_fraction
    NaN — silently, because 0 is a legal count.
    """

    def _barcoded(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        barcodes = [f"C1_T_1_1_0_c1_v2_id-{i:016d}" for i in range(matrix.shape[0])]
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=sample_id,
            target_genes=TARGETS, index=barcodes,
        )
        patient = np.full(matrix.shape[0], "P1")
        tissue = np.where(sample_id == "s1", "tumour", "normal")
        return labels, patient, tissue, compartment

    def test_epithelial_counts_are_not_zero(self, cohort):
        labels, patient, tissue, compartment = self._barcoded(cohort)
        counts = mature_cell_counts(labels, patient_id=patient, tissue=tissue)
        assert (counts["n_cells_epithelial"] > 0).all()
        assert int(counts[counts["granularity_rung"] == "epithelial"]
                   ["n_cells_epithelial"].sum()) == int(
            (compartment == "epithelial").sum()) * len(TRANSCRIPT_AXES)

    def test_mature_fraction_is_computed(self, cohort):
        labels, patient, tissue, _ = self._barcoded(cohort)
        counts = mature_cell_counts(labels, patient_id=patient, tissue=tissue)
        assert counts["mature_fraction"].notna().all()
        assert counts["mature_fraction"].between(0, 1).all()

    def test_the_epithelial_rung_is_a_fraction_of_one(self, cohort):
        """Coarsest rung calls all epithelium mature, so the fraction is exactly 1."""
        labels, patient, tissue, _ = self._barcoded(cohort)
        counts = mature_cell_counts(labels, patient_id=patient, tissue=tissue)
        coarsest = counts[counts["granularity_rung"] == "epithelial"]
        assert (coarsest["mature_fraction"] == 1.0).all()

    def test_matches_the_rangeindex_result(self, cohort):
        """Index must not change any number."""
        matrix, compartment, sample_id, _ = cohort
        patient = np.full(matrix.shape[0], "P1")
        tissue = np.where(sample_id == "s1", "tumour", "normal")
        plain = assign_labels(matrix, GENES, compartment=compartment,
                              sample_id=sample_id, target_genes=TARGETS)
        barcoded, _, _, _ = self._barcoded(cohort)
        pd.testing.assert_frame_equal(
            mature_cell_counts(plain, patient_id=patient, tissue=tissue),
            mature_cell_counts(barcoded, patient_id=patient, tissue=tissue),
        )

    def test_maturity_summary_also_survives_it(self, cohort):
        labels, patient, tissue, _ = self._barcoded(cohort)
        summary = maturity_summary(
            labels, patient_id=patient, tissue=tissue, study_id="GSE178341"
        )
        assert len(summary) > 0
        assert summary["frac_mature_normal"].notna().all()
        assert summary["frac_mature_tumour"].notna().all()


class TestCompositionalSignalIsRecoverable:
    """THE test this module was missing.

    Every earlier test checked shapes, monotonicity across rungs, and that the
    axes disagreed. None checked that a known compositional difference could be
    RECOVERED — which is the only thing the labels exist to enable. Its absence
    is why per-sample quantile binning, which pins the mature fraction to the
    quantile and makes Delta(mature fraction) identically zero, survived to real
    data.
    """

    def _paired(self, frac_mature_normal: float, frac_mature_tumour: float,
                n: int = 1200, depth: int = 5_000):
        """One patient, two arms, with a KNOWN difference in mature fraction.

        "Mature" cells carry no stem-marker expression; "stem" cells carry a lot.
        Only the mixing proportion differs between arms.
        """
        stem_profile = np.full(len(GENES), 1.0)
        for gene in STEM:
            stem_profile[GENES.index(gene)] = 60.0
        mature_profile = np.full(len(GENES), 1.0)

        blocks, tissue = [], []
        for arm, frac in (("normal", frac_mature_normal), ("tumour", frac_mature_tumour)):
            n_mature = int(round(n * frac))
            for profile, count in ((mature_profile, n_mature), (stem_profile, n - n_mature)):
                if count == 0:
                    continue
                p = profile / profile.sum()
                blocks.append(RNG.poisson(np.outer(np.full(count, depth), p)).astype(np.int64))
                tissue += [arm] * count
        matrix = np.vstack(blocks)
        tissue = np.array(tissue, dtype=object)
        return (
            matrix,
            np.full(matrix.shape[0], "epithelial", dtype=object),
            np.full(matrix.shape[0], "P1", dtype=object),
            tissue,
        )

    def _fraction(self, frac_normal, frac_tumour, rung="lineage"):
        matrix, compartment, patient, tissue = self._paired(frac_normal, frac_tumour)
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=[rung],
        )
        counts = mature_cell_counts(
            labels, patient_id=patient, tissue=tissue, axes=["stem_pole"], rungs=[rung]
        ).set_index("tissue")
        return counts.loc["normal", "mature_fraction"], counts.loc["tumour", "mature_fraction"]

    def test_a_real_compositional_loss_is_detected(self):
        """Normal 50% mature, tumour 10%. The tumour fraction must come back
        materially lower — this is the compositional term."""
        normal, tumour = self._fraction(0.50, 0.10)
        assert tumour < normal - 0.2, f"normal={normal:.3f} tumour={tumour:.3f}"

    def test_no_difference_is_reported_when_there_is_none(self):
        """Equal composition must NOT manufacture a difference."""
        normal, tumour = self._fraction(0.50, 0.50)
        assert abs(tumour - normal) < 0.1

    def test_a_compositional_gain_is_detected_too(self):
        """The estimator must not be one-directional."""
        normal, tumour = self._fraction(0.30, 0.70)
        assert tumour > normal + 0.2

    @pytest.mark.parametrize("rung", ["lineage", "crypt_position"])
    def test_the_signal_survives_at_every_binned_rung(self, rung):
        normal, tumour = self._fraction(0.60, 0.15, rung=rung)
        assert tumour < normal - 0.15, f"{rung}: normal={normal:.3f} tumour={tumour:.3f}"

    def test_the_delta_tracks_the_injected_magnitude(self):
        """A bigger true loss must produce a bigger measured loss."""
        small_n, small_t = self._fraction(0.50, 0.40)
        large_n, large_t = self._fraction(0.50, 0.05)
        assert (large_n - large_t) > (small_n - small_t)

    def test_per_sample_quantiles_would_have_failed_this(self):
        """Pins the regression. The reference arm sits near the intended
        quantile; the tumour arm must be free to leave it."""
        normal, tumour = self._fraction(0.50, 0.10)
        assert abs(normal - 0.5) < 0.15      # reference defines the cut
        assert tumour < 0.35                  # tumour is NOT pinned to it


class TestTieCollapse:
    """Coincident quantile cuts must not silently delete a bin.

    On the pilot, stem_pole/crypt_position came back with only crypt_bottom and
    crypt_top — crypt_middle never appeared — and the result was byte-identical
    to the lineage rung. Cells with zero counts across all five axis-1 markers
    share a score, so both tertile boundaries landed inside that tie.
    """

    def _tied(self, n_tied: int = 900, n_spread: int = 100):
        """A cohort where most epithelial cells have no axis-1 marker counts."""
        spread = np.full(len(GENES), 1.0)
        for i, gene in enumerate(STEM):
            spread[GENES.index(gene)] = 20.0 * (i + 1)
        flat = np.full(len(GENES), 1.0)
        for gene in STEM:
            flat[GENES.index(gene)] = 0.0

        blocks = []
        for profile, n in ((flat, n_tied), (spread, n_spread)):
            p = profile / profile.sum()
            blocks.append(RNG.poisson(np.outer(np.full(n, 5000), p)).astype(np.int64))
        matrix = np.vstack(blocks)
        total = matrix.shape[0]
        return (
            matrix,
            np.full(total, "epithelial", dtype=object),
            np.full(total, "P1", dtype=object),
            np.array(["normal"] * total, dtype=object),
        )

    def test_the_diagnostic_reports_the_tie(self):
        matrix, compartment, _, _ = self._tied()
        stats = axis_tie_fraction(
            matrix, GENES, "stem_pole", target_genes=TARGETS,
            epithelial=compartment == "epithelial",
        )
        assert stats["tied_fraction"] > 0.5
        assert stats["largest_tied_block"] >= 900

    def test_no_bin_is_silently_empty(self):
        """Whatever bins appear, none of the declared ones may be missing without
        the fallback having been taken."""
        matrix, compartment, patient, tissue = self._tied()
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["crypt_position"],
        )
        values = set(labels[label_column("stem_pole", "crypt_position")].astype(str))
        values -= {NON_EPITHELIAL, UNRESOLVED}
        assert values <= set(RUNG_SPECS["crypt_position"].bins)
        # The extremes must both survive the fallback.
        assert RUNG_SPECS["crypt_position"].bins[0] in values
        assert RUNG_SPECS["crypt_position"].mature in values

    def test_the_extremes_are_kept_not_the_middle(self):
        """A binary fallback must span bottom-to-top, not collapse onto adjacent
        bins — the mature end has to stay the mature end."""
        matrix, compartment, patient, tissue = self._tied()
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["crypt_position"],
        )
        column = labels[label_column("stem_pole", "crypt_position")].astype(str)
        if "crypt_middle" not in set(column):
            assert {"crypt_bottom", "crypt_top"} <= set(column)

    def test_an_unresolvable_axis_still_labels_every_cell(self):
        """All markers zero everywhere: no gradient exists at all."""
        matrix = RNG.poisson(5, size=(200, len(GENES))).astype(np.int64)
        for gene in STEM:
            matrix[:, GENES.index(gene)] = 0
        compartment = np.full(200, "epithelial", dtype=object)
        tissue = np.full(200, "normal", dtype=object)
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue,
            patient_id=np.full(200, "P1", dtype=object),
            axes=["stem_pole"], rungs=["crypt_position"],
        )
        assert labels.notna().all().all()

    def test_a_well_spread_axis_keeps_all_three_bins(self):
        """The fallback must not fire when the data does support three bins."""
        matrix, compartment, patient, tissue = self._tied(n_tied=0, n_spread=900)
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["crypt_position"],
        )
        values = set(labels[label_column("stem_pole", "crypt_position")].astype(str))
        assert set(RUNG_SPECS["crypt_position"].bins) <= values


class TestRungDegeneracy:
    """Two rungs that draw the same boundary are one point on the curve, not two.

    This is not hypothetical: on the pilot at `depth_quantile=0.10`, stem_pole's
    `lineage` and `crypt_position` returned identical mature sets in all ten
    arms, because both of crypt_position's tertile cuts landed inside the block
    of cells with no stem-marker counts. Nothing raised — each rung is
    individually well-formed — and the granularity curve silently reported three
    resolutions where the data supported two.
    """

    def _labels(self, rungs=("lineage", "crypt_position"), **kwargs):
        matrix, compartment, patient, tissue = TestTieCollapse()._tied(**kwargs)
        return assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=list(rungs),
        )

    def test_it_catches_the_collapse_it_exists_to_catch(self):
        """The pilot's actual failure, reproduced."""
        labels = self._labels()
        report = rung_degeneracy(
            labels, axes=["stem_pole"], rungs=["lineage", "crypt_position"]
        )
        collapsed = report[report["identical"]]
        assert len(collapsed) == 1
        assert set(collapsed.iloc[0][["rung_a", "rung_b"]]) == {
            "lineage", "crypt_position"
        }
        assert collapsed.iloc[0]["jaccard"] == pytest.approx(1.0)

    def test_a_well_spread_axis_is_not_degenerate(self):
        """No false alarm when the score really does support three bins."""
        labels = self._labels(n_tied=0, n_spread=900)
        report = rung_degeneracy(
            labels, axes=["stem_pole"], rungs=["lineage", "crypt_position"]
        )
        assert not report["identical"].any()
        assert (report["jaccard"] < 1.0).all()

    def test_it_covers_every_pair_once(self):
        rungs = ["epithelial", "lineage", "crypt_position"]
        labels = self._labels(rungs=rungs)
        report = rung_degeneracy(labels, axes=["stem_pole"], rungs=rungs)
        assert len(report) == 3
        pairs = {frozenset(r) for r in report[["rung_a", "rung_b"]].to_numpy()}
        assert len(pairs) == 3

    def test_missing_columns_are_skipped_not_raised(self):
        """A run that labelled only some rungs still gets a report."""
        labels = self._labels()
        report = rung_degeneracy(
            labels, axes=["stem_pole"], rungs=["lineage", "best4", "crypt_position"]
        )
        assert set(report["rung_a"]) | set(report["rung_b"]) == {
            "lineage", "crypt_position"
        }


class TestDepthConfounding:
    """Zero counts stay zero after depth normalisation, so a shallow cell is more
    likely to be called mature on a sparsely detected axis. On the pilot, axis 1's
    mature bin was exactly the tied block of cells with no stem-marker counts.
    """

    def _cohort_with_depth_split(self, deep_depth=20_000, shallow_depth=1_500):
        """Stem-marker-positive cells sequenced deep, marker-free cells shallow —
        the confound in its purest form."""
        stem = np.full(len(GENES), 1.0)
        for gene in STEM:
            stem[GENES.index(gene)] = 40.0
        flat = np.full(len(GENES), 1.0)
        for gene in STEM:
            flat[GENES.index(gene)] = 0.0

        blocks = []
        for profile, depth, n in ((stem, deep_depth, 400), (flat, shallow_depth, 400)):
            p = profile / profile.sum()
            blocks.append(RNG.poisson(np.outer(np.full(n, depth), p)).astype(np.int64))
        matrix = np.vstack(blocks)
        total = matrix.shape[0]
        return (
            matrix,
            np.full(total, "epithelial", dtype=object),
            np.full(total, "P1", dtype=object),
            np.full(total, "normal", dtype=object),
        )

    def _run(self, matrix, compartment, patient, tissue):
        from src.reference.qc import cell_qc_metrics

        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["lineage"],
        )
        metrics = cell_qc_metrics(matrix, GENES, batch=tissue)
        return label_depth_confounding(labels, metrics, axes=["stem_pole"],
                                       rungs=["lineage"])

    def test_a_depth_confound_is_detected(self):
        report = self._run(*self._cohort_with_depth_split())
        assert len(report) == 1
        assert bool(report.iloc[0]["flagged"])
        assert report.iloc[0]["counts_ratio"] < 1.0   # mature cells are shallower

    def test_balanced_depth_is_not_flagged(self):
        report = self._run(*self._cohort_with_depth_split(20_000, 20_000))
        assert not bool(report.iloc[0]["flagged"])
        assert report.iloc[0]["counts_ratio"] == pytest.approx(1.0, abs=0.25)

    def test_reports_both_counts_and_genes(self):
        report = self._run(*self._cohort_with_depth_split())
        assert {"median_counts_mature", "median_counts_other",
                "median_genes_mature", "median_genes_other"} <= set(report.columns)

    def test_metrics_length_mismatch_raises(self, cohort):
        import pandas as pd

        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment,
                               sample_id=sample_id, target_genes=TARGETS)
        with pytest.raises(LabelError, match="rows for"):
            label_depth_confounding(labels, pd.DataFrame({"n_counts": [1], "n_genes": [1]}))


class TestDepthMatchingRemovesTheConfound:
    """The pilot's headline problem, and the proof the fix works.

    Unmatched, axis 1's mature cells came back four times shallower than its
    non-mature cells (median 4,791 counts against 18,829) — the maturity call was
    substantially a depth measurement. Binomial thinning to a common depth should
    drive that ratio back toward 1.
    """

    def _confounded_cohort(self, n_per_group: int = 400):
        """Depth varies INDEPENDENTLY of the true population.

        Half the cells are truly stem-like, half truly mature, and each cell is
        sequenced deep or shallow at random. The markers are deliberately sparse
        — about 0.5 expected counts per marker in a shallow cell — so a shallow
        stem cell often drops out to zero and gets miscalled mature. Any
        association between the maturity call and depth is therefore pure
        artifact, which is what depth matching has to remove.

        A fixture where depth and biology move together (which was my first
        attempt) cannot test this: there is nothing to separate.
        """
        stem = np.full(len(GENES), 1.0)
        for gene in STEM:
            stem[GENES.index(gene)] = 0.01      # sparse, like the real markers
        mature = np.full(len(GENES), 1.0)
        for gene in STEM:
            mature[GENES.index(gene)] = 0.0

        rng = np.random.default_rng(4242)
        blocks = []
        for profile in (stem, mature):
            p = profile / profile.sum()
            depths = rng.choice([3_000, 20_000], size=n_per_group)   # independent
            blocks.append(rng.poisson(np.outer(depths, p)).astype(np.int64))
        matrix = np.vstack(blocks)
        total = matrix.shape[0]
        return (
            matrix,
            np.full(total, "epithelial", dtype=object),
            np.full(total, "P1", dtype=object),
            np.full(total, "normal", dtype=object),
        )

    def _miscalls_by_depth(self, depth_target):
        """Ground truth is available here, so measure the property that matters:
        are the errors depth-dependent?

        Depth matching does not reduce the error rate — it makes the error
        UNBIASED with respect to depth. Unmatched, every miscalled stem cell is a
        shallow one; matched, the miscalls split evenly.
        """
        matrix, compartment, patient, tissue = self._confounded_cohort()
        totals = matrix.sum(axis=1)
        truth = np.array(["stem"] * 400 + ["mature"] * 400, dtype=object)

        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["lineage"], depth_target=depth_target,
        )
        called_mature = mature_mask(labels, "stem_pole", "lineage")
        wrong = (truth == "stem") & called_mature      # stem cells called mature
        shallow = totals < 10_000
        return int((wrong & shallow).sum()), int((wrong & ~shallow).sum())

    def test_unmatched_errors_are_all_in_shallow_cells(self):
        """The confound in its purest form: dropout, not biology."""
        shallow, deep = self._miscalls_by_depth(depth_target=0)
        assert shallow > 0
        assert deep == 0, f"expected all miscalls shallow, got {shallow}/{deep}"

    def test_matching_makes_the_errors_depth_balanced(self):
        """Not fewer errors — unbiased ones."""
        shallow, deep = self._miscalls_by_depth(depth_target=2_500)
        assert shallow > 0 and deep > 0
        assert abs(shallow - deep) / (shallow + deep) < 0.35, f"{shallow}/{deep}"

    def test_a_lower_target_balances_them_further(self):
        shallow, deep = self._miscalls_by_depth(depth_target=1_500)
        assert abs(shallow - deep) / (shallow + deep) < 0.2, f"{shallow}/{deep}"

    def test_the_auc_statistic_survives_bimodal_depth(self):
        """counts_ratio flips its median when a bin splits near 50/50 on a
        bimodal depth distribution. depth_auc must not."""
        from src.reference.qc import cell_qc_metrics

        matrix, compartment, patient, tissue = self._confounded_cohort()
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["lineage"], depth_target=1_500,
        )
        report = label_depth_confounding(
            labels, cell_qc_metrics(matrix, GENES, batch=tissue),
            axes=["stem_pole"], rungs=["lineage"],
        )
        assert abs(float(report.iloc[0]["depth_auc"]) - 0.5) < 0.1

    def test_a_strong_monotone_confound_is_still_caught(self):
        """The pilot's case: a 4x median gap, not a median flip. Both statistics
        should fire."""
        from src.reference.qc import cell_qc_metrics

        stem = np.full(len(GENES), 1.0)
        for gene in STEM:
            stem[GENES.index(gene)] = 40.0
        flat = np.full(len(GENES), 1.0)
        for gene in STEM:
            flat[GENES.index(gene)] = 0.0
        blocks = []
        for profile, depth, n in ((stem, 20_000, 400), (flat, 3_000, 400)):
            p = profile / profile.sum()
            blocks.append(RNG.poisson(np.outer(np.full(n, depth), p)).astype(np.int64))
        matrix = np.vstack(blocks)
        total = matrix.shape[0]
        compartment = np.full(total, "epithelial", dtype=object)
        tissue = np.full(total, "normal", dtype=object)
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue,
            patient_id=np.full(total, "P1", dtype=object),
            axes=["stem_pole"], rungs=["lineage"], depth_target=0,
        )
        report = label_depth_confounding(
            labels, cell_qc_metrics(matrix, GENES, batch=tissue),
            axes=["stem_pole"], rungs=["lineage"],
        )
        assert bool(report.iloc[0]["flagged"])
        assert float(report.iloc[0]["counts_ratio"]) < 0.75

    def test_shallow_cells_become_unresolved_not_immature(self):
        """A cell that could not be measured is not a cell measured to be
        immature — open decision #14's distinction, in code."""
        matrix, compartment, patient, tissue = self._confounded_cohort()
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["lineage"], depth_target=10_000,
        )
        values = labels[label_column("stem_pole", "lineage")].astype(str)
        assert (values == UNRESOLVED).sum() > 0
        # Unresolved cells are neither mature nor counted as the immature bin.
        assert not mature_mask(labels, "stem_pole", "lineage")[values == UNRESOLVED].any()

    def test_counts_report_the_unresolved_share(self):
        matrix, compartment, patient, tissue = self._confounded_cohort()
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["lineage"], depth_target=10_000,
        )
        counts = mature_cell_counts(
            labels, patient_id=patient, tissue=tissue,
            axes=["stem_pole"], rungs=["lineage"],
        )
        row = counts.iloc[0]
        assert row["n_cells_unresolved"] > 0
        assert row["unresolved_fraction"] > 0
        assert row["n_cells_resolved"] + row["n_cells_unresolved"] == row["n_cells_epithelial"]
        # The fraction is over resolved cells, so it is not diluted by them.
        assert row["mature_fraction"] <= 1.0

    def test_matching_is_deterministic(self):
        """Thinning is random; the seed must make it reproducible."""
        matrix, compartment, patient, tissue = self._confounded_cohort()
        kwargs = dict(
            compartment=compartment, sample_id=tissue, target_genes=TARGETS,
            tissue=tissue, patient_id=patient, axes=["stem_pole"],
            rungs=["lineage"], depth_target=2_500,
        )
        a = assign_labels(matrix, GENES, seed=7, **kwargs)
        b = assign_labels(matrix, GENES, seed=7, **kwargs)
        pd.testing.assert_frame_equal(a, b)

    def test_the_compositional_signal_survives_matching(self):
        """Matching must not destroy real signal while removing the artifact."""
        from tests.test_reference_labels import TestCompositionalSignalIsRecoverable as T

        maker = T()
        matrix, compartment, patient, tissue = maker._paired(0.50, 0.10)
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["lineage"], depth_target=2_000,
        )
        counts = mature_cell_counts(
            labels, patient_id=patient, tissue=tissue,
            axes=["stem_pole"], rungs=["lineage"],
        ).set_index("tissue")
        normal = counts.loc["normal", "mature_fraction"]
        tumour = counts.loc["tumour", "mature_fraction"]
        assert tumour < normal - 0.2, f"normal={normal:.3f} tumour={tumour:.3f}"


def test_depth_diagnostic_excludes_unresolved_cells():
    """UNRESOLVED cells are defined by low depth, so including them in the
    comparison guarantees a perfect association. The pilot reported
    counts_ratio 11.8 and AUC 1.000 for the epithelial rung, which was the
    diagnostic measuring its own definition rather than the labels.
    """
    from src.reference.qc import cell_qc_metrics

    rng = np.random.default_rng(11)
    profile = np.full(len(GENES), 1.0)
    for gene in STEM:
        profile[GENES.index(gene)] = 0.01
    p = profile / profile.sum()
    depths = np.concatenate([np.full(300, 20_000), np.full(300, 900)])
    matrix = rng.poisson(np.outer(depths, p)).astype(np.int64)
    total = matrix.shape[0]
    compartment = np.full(total, "epithelial", dtype=object)
    tissue = np.full(total, "normal", dtype=object)

    labels = assign_labels(
        matrix, GENES, compartment=compartment, sample_id=tissue,
        target_genes=TARGETS, tissue=tissue,
        patient_id=np.full(total, "P1", dtype=object),
        axes=["stem_pole"], rungs=["epithelial"], depth_target=5_000,
    )
    assert (labels[label_column("stem_pole", "epithelial")].astype(str)
            == UNRESOLVED).sum() > 0

    report = label_depth_confounding(
        labels, cell_qc_metrics(matrix, GENES, batch=tissue),
        axes=["stem_pole"], rungs=["epithelial"],
    )
    # Every resolved cell is mature at this rung, so there is no comparison to
    # make and the row is omitted rather than reported as a perfect confound.
    assert len(report) == 0 or float(report.iloc[0]["depth_auc"]) < 0.99


class TestTieDiagnosticMatchesAssignLabels:
    """axis_tie_fraction's sweep must report what a real run would produce.

    _thin_to_depth clips the thinning probability at 1, so a cell below the
    target passes through unthinned. Leaving those in the tie computation mixes
    thinned and unthinned cells and reports a number no actual run would give —
    assign_labels drops them as unresolved_depth.
    """

    def _cohort(self):
        rng = np.random.default_rng(31)
        profile = np.full(len(GENES), 1.0)
        for gene in STEM:
            profile[GENES.index(gene)] = 0.02
        p = profile / profile.sum()
        depths = np.concatenate([np.full(300, 20_000), np.full(300, 800)])
        matrix = rng.poisson(np.outer(depths, p)).astype(np.int64)
        total = matrix.shape[0]
        return (
            matrix,
            np.full(total, "epithelial", dtype=object),
            np.full(total, "P1", dtype=object),
            np.full(total, "normal", dtype=object),
        )

    def test_the_diagnostic_excludes_sub_target_cells(self):
        matrix, compartment, _, _ = self._cohort()
        stats = axis_tie_fraction(
            matrix, GENES, "stem_pole", target_genes=TARGETS,
            epithelial=compartment == "epithelial", depth_target=5_000,
        )
        # Only the 300 deep cells clear 5,000, so that is what was scored.
        assert stats["n_cells"] == 300

    def test_its_cell_count_matches_what_assign_labels_resolves(self):
        matrix, compartment, patient, tissue = self._cohort()
        target = 5_000.0
        stats = axis_tie_fraction(
            matrix, GENES, "stem_pole", target_genes=TARGETS,
            epithelial=compartment == "epithelial", depth_target=target,
        )
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["lineage"], depth_target=target,
        )
        column = labels[label_column("stem_pole", "lineage")].astype(str)
        resolved = int((~column.isin([NON_EPITHELIAL, UNRESOLVED])).sum())
        assert stats["n_cells"] == resolved

    def test_no_target_means_no_exclusion(self):
        matrix, compartment, _, _ = self._cohort()
        stats = axis_tie_fraction(
            matrix, GENES, "stem_pole", target_genes=TARGETS,
            epithelial=compartment == "epithelial",
        )
        assert stats["n_cells"] == 600

    def test_an_unreachable_target_raises(self):
        matrix, compartment, _, _ = self._cohort()
        with pytest.raises(LabelError, match="no cell reaches"):
            axis_tie_fraction(
                matrix, GENES, "stem_pole", target_genes=TARGETS,
                epithelial=compartment == "epithelial", depth_target=1e9,
            )


class TestWithinDepthStrata:
    """Separating the two explanations label_depth_confounding cannot.

    A maturity call driven purely by dropout is all-or-nothing inside a narrow
    depth band. One driven by biology keeps a mix in every band.
    """

    def _labels_and_metrics(self, matrix, compartment, patient, tissue):
        from src.reference.qc import cell_qc_metrics

        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue, patient_id=patient,
            axes=["stem_pole"], rungs=["lineage"], depth_target=0,
        )
        return labels, cell_qc_metrics(matrix, GENES, batch=tissue)

    def test_a_purely_technical_call_shows_a_steep_depth_gradient(self):
        """Marker expression identical in every cell; only depth varies. The
        mature fraction must then fall steeply from shallow to deep.

        Note what this does NOT show: dropout is stochastic, so even here every
        stratum keeps a mix. That is why this function cannot separate technical
        from biological confounding — annotation_concordance does that."""
        rng = np.random.default_rng(19)
        profile = np.full(len(GENES), 1.0)
        for gene in STEM:
            profile[GENES.index(gene)] = 0.01
        p = profile / profile.sum()
        depths = rng.integers(600, 40_000, size=800)
        matrix = rng.poisson(np.outer(depths, p)).astype(np.int64)
        total = matrix.shape[0]
        compartment = np.full(total, "epithelial", dtype=object)
        tissue = np.full(total, "normal", dtype=object)
        labels, metrics = self._labels_and_metrics(
            matrix, compartment, np.full(total, "P1", dtype=object), tissue
        )
        out = maturity_within_depth_strata(labels, metrics)
        shallow = float(out.iloc[0]["mature_fraction"])
        deep = float(out.iloc[-1]["mature_fraction"])
        assert shallow > deep + 0.3, f"shallow={shallow:.2f} deep={deep:.2f}"

    def test_a_biological_call_keeps_a_mix_in_every_stratum(self):
        """Two true populations present at every depth. The call must then
        still separate cells inside each band."""
        rng = np.random.default_rng(23)
        stem = np.full(len(GENES), 1.0)
        for gene in STEM:
            stem[GENES.index(gene)] = 3.0          # well detected, not sparse
        mature = np.full(len(GENES), 1.0)
        for gene in STEM:
            mature[GENES.index(gene)] = 0.0

        blocks = []
        for profile in (stem, mature):
            p = profile / profile.sum()
            depths = rng.integers(600, 40_000, size=400)   # same depth range
            blocks.append(rng.poisson(np.outer(depths, p)).astype(np.int64))
        matrix = np.vstack(blocks)
        total = matrix.shape[0]
        compartment = np.full(total, "epithelial", dtype=object)
        tissue = np.full(total, "normal", dtype=object)
        labels, metrics = self._labels_and_metrics(
            matrix, compartment, np.full(total, "P1", dtype=object), tissue
        )
        out = maturity_within_depth_strata(labels, metrics)
        # Both populations present at every depth, so the gradient stays shallow.
        spread = out["mature_fraction"].max() - out["mature_fraction"].min()
        assert spread < 0.4, out["mature_fraction"].tolist()

    def test_it_reports_one_row_per_stratum_with_ranges(self):
        rng = np.random.default_rng(5)
        matrix = rng.poisson(3, size=(400, len(GENES))).astype(np.int64)
        total = matrix.shape[0]
        compartment = np.full(total, "epithelial", dtype=object)
        tissue = np.full(total, "normal", dtype=object)
        labels, metrics = self._labels_and_metrics(
            matrix, compartment, np.full(total, "P1", dtype=object), tissue
        )
        out = maturity_within_depth_strata(labels, metrics, n_strata=5)
        assert {"stratum", "n_cells", "counts_low", "counts_high",
                "mature_fraction"} == set(out.columns)
        assert (out["counts_high"] >= out["counts_low"]).all()

    def test_metrics_mismatch_raises(self, cohort):
        matrix, compartment, sample_id, _ = cohort
        labels = assign_labels(matrix, GENES, compartment=compartment,
                               sample_id=sample_id, target_genes=TARGETS)
        with pytest.raises(LabelError, match="rows for"):
            maturity_within_depth_strata(
                labels, pd.DataFrame({"n_counts": [1.0]})
            )


class TestAnnotationConcordance:
    """Signal versus dropout noise — the question the depth diagnostics cannot
    answer, because stochastic dropout produces a mix at every depth just as
    real variation would.

    A call made of noise has nothing to agree with. One measuring real maturity
    tracks an independently-derived annotation.
    """

    def _labels(self, separable: bool):
        """separable=True: two real populations. False: pure dropout noise."""
        rng = np.random.default_rng(41)
        stem = np.full(len(GENES), 1.0)
        for gene in STEM:
            stem[GENES.index(gene)] = 3.0 if separable else 0.01
        mature = np.full(len(GENES), 1.0)
        for gene in STEM:
            mature[GENES.index(gene)] = 0.0 if separable else 0.01

        blocks, truth = [], []
        for name, profile in (("cE01 (Stem/TA-like)", stem), ("cE07 (Mature)", mature)):
            p = profile / profile.sum()
            depths = rng.integers(2_000, 30_000, size=400)
            blocks.append(rng.poisson(np.outer(depths, p)).astype(np.int64))
            truth += [name] * 400
        matrix = np.vstack(blocks)
        total = matrix.shape[0]
        compartment = np.full(total, "epithelial", dtype=object)
        tissue = np.full(total, "normal", dtype=object)
        labels = assign_labels(
            matrix, GENES, compartment=compartment, sample_id=tissue,
            target_genes=TARGETS, tissue=tissue,
            patient_id=np.full(total, "P1", dtype=object),
            axes=["stem_pole"], rungs=["lineage"], depth_target=0,
        )
        return labels, np.array(truth, dtype=object)

    def test_a_real_call_agrees_with_the_annotation(self):
        labels, annotation = self._labels(separable=True)
        out = annotation_concordance(labels, annotation)
        assert out["kappa"] > 0.5, out
        assert out["informative"]

    def test_a_dropout_call_does_not(self):
        """Identical marker rates in both populations, so any 'maturity' the
        call finds is Poisson noise. Kappa must be near zero."""
        labels, annotation = self._labels(separable=False)
        out = annotation_concordance(labels, annotation)
        assert abs(out["kappa"]) < 0.2, out
        assert not out["informative"]

    def test_it_reports_the_full_two_by_two(self):
        labels, annotation = self._labels(separable=True)
        out = annotation_concordance(labels, annotation)
        assert {"n_cells", "agreement", "sensitivity", "specificity",
                "kappa", "informative"} <= set(out)
        assert 0.0 <= out["agreement"] <= 1.0

    def test_kappa_not_raw_agreement_is_the_verdict(self):
        """Raw agreement is inflated when one class dominates, which is exactly
        the situation at the best4 rung."""
        labels, annotation = self._labels(separable=False)
        out = annotation_concordance(labels, annotation)
        assert out["agreement"] > out["kappa"]

    def test_length_mismatch_raises(self):
        labels, _ = self._labels(separable=True)
        with pytest.raises(LabelError, match="entries for"):
            annotation_concordance(labels, ["cE01 (Stem/TA-like)"])
