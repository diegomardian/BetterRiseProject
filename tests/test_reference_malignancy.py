"""Malignant vs. normal epithelium. W1.

The test that matters is `TestNormalEpitheliumValidation`. execution_plan.md §4
makes "normal epithelium not misread as tumour" the done-when for this stage, and
it is only meaningful because the CNV baseline is non-epithelial — a population
used as the reference is non-malignant by construction, so validating on it would
prove nothing.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.reference.malignancy import (
    DIPLOID_COMPARTMENTS,
    MIN_REFERENCE_CELLS,
    MalignancyError,
    assign_cnv_roles,
    call_malignancy,
    cleanup_infercnv_run,
    infercnv_reference_groups,
    read_infercnv_score_table,
    read_infercnv_scores,
    run_infercnv,
    select_cnv_reference,
    validate_normal_epithelium,
    write_infercnv_inputs,
)

RNG = np.random.default_rng(20260818)


def cohort(
    n_reference: int = 200,
    n_normal_epi: int = 300,
    n_tumour_epi: int = 300,
    malignant_shift: float = 3.0,
    patient: str = "P1",
):
    """Diploid reference cells, normal epithelium that matches them, and tumour
    epithelium shifted upward in CNV score."""
    reference = RNG.normal(1.0, 0.15, n_reference)
    normal_epi = RNG.normal(1.0, 0.15, n_normal_epi)
    tumour_epi = RNG.normal(1.0 + malignant_shift * 0.15, 0.15, n_tumour_epi)

    return pd.DataFrame(
        {
            "cnv_score": np.concatenate([reference, normal_epi, tumour_epi]),
            "compartment": (
                ["immune"] * n_reference
                + ["epithelial"] * (n_normal_epi + n_tumour_epi)
            ),
            "tissue": (
                ["normal"] * (n_reference + n_normal_epi) + ["tumour"] * n_tumour_epi
            ),
            "patient_id": [patient] * (n_reference + n_normal_epi + n_tumour_epi),
        }
    )


class TestReferenceSelection:
    def test_diploid_compartments_are_additional_categories_not_the_baseline(self):
        """They supplement the matched-normal baseline rather than replacing it.
        Kept separate so inferCNV's per-category bounding can work."""
        assert "epithelial" not in DIPLOID_COMPARTMENTS
        assert set(DIPLOID_COMPARTMENTS) == {"immune", "stromal", "endothelial"}

    def test_counts_reference_cells_per_patient(self):
        df = cohort()
        out = select_cnv_reference(df["compartment"], patient_id=df["patient_id"])
        assert int(out.loc[0, "n_reference"]) == 200
        assert int(out.loc[0, "n_epithelial"]) == 600
        assert bool(out.loc[0, "usable"])

    def test_a_patient_without_enough_reference_is_marked_unusable(self):
        df = cohort(n_reference=10)
        out = select_cnv_reference(df["compartment"], patient_id=df["patient_id"])
        assert not bool(out.loc[0, "usable"])

    def test_all_reference_compartments_count(self):
        df = pd.DataFrame(
            {
                "compartment": ["immune", "stromal", "endothelial", "epithelial"],
                "patient_id": ["P1"] * 4,
            }
        )
        out = select_cnv_reference(df["compartment"], patient_id=df["patient_id"],
                                   min_cells=3)
        assert int(out.loc[0, "n_reference"]) == 3


class TestCalls:
    def test_tumour_epithelium_is_called_malignant(self):
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        tumour = calls[df["tissue"].to_numpy() == "tumour"]
        assert (tumour["call"].astype(str) == "malignant").mean() > 0.65

    def test_the_strict_threshold_trades_sensitivity_for_specificity(self):
        """MALIGNANT_QUANTILE = 0.99 is deliberate. At a 3-sigma separation it
        recovers about three quarters of malignant cells, and that is the right
        trade: a false MALIGNANT call moves a normal cell into the tumour arm,
        which is the direction that inflates the compositional term and the
        direction of the prior hypothesis.
        """
        df = cohort()
        strict = call_malignancy(
            df["cnv_score"], compartment=df["compartment"],
            patient_id=df["patient_id"], quantile=0.99,
        )
        loose = call_malignancy(
            df["cnv_score"], compartment=df["compartment"],
            patient_id=df["patient_id"], quantile=0.80,
        )
        is_tumour = df["tissue"].to_numpy() == "tumour"
        is_normal_epi = (df["tissue"].to_numpy() == "normal") & (
            df["compartment"].to_numpy() == "epithelial"
        )
        sensitivity = [
            float((c[is_tumour]["call"].astype(str) == "malignant").mean())
            for c in (strict, loose)
        ]
        false_positive = [
            float((c[is_normal_epi]["call"].astype(str) == "malignant").mean())
            for c in (strict, loose)
        ]
        assert sensitivity[1] > sensitivity[0]        # loose catches more
        assert false_positive[1] > false_positive[0]  # and costs specificity
        assert false_positive[0] < 0.05               # strict keeps it low

    def test_reference_cells_are_labelled_reference_not_non_malignant(self):
        """They defined the baseline; calling them non-malignant is circular."""
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        reference = calls[calls["compartment"] == "immune"]
        assert (reference["call"].astype(str) == "reference").all()

    def test_confidence_scales_with_distance_past_the_threshold(self):
        weak = cohort(malignant_shift=1.5)
        strong = cohort(malignant_shift=6.0)
        out = []
        for df in (weak, strong):
            calls = call_malignancy(
                df["cnv_score"], compartment=df["compartment"],
                patient_id=df["patient_id"],
            )
            tumour = calls[df["tissue"].to_numpy() == "tumour"]
            out.append(float(tumour["confidence"].median()))
        assert out[1] > out[0]

    def test_threshold_is_per_patient(self):
        """A baseline pooled across patients would fold germline copy-number
        variation and per-patient capture differences into the call."""
        noisy = cohort(patient="P2")
        noisy["cnv_score"] = noisy["cnv_score"] * 3.0
        df = pd.concat([cohort(patient="P1"), noisy], ignore_index=True)
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        thresholds = calls.groupby("patient_id", observed=True)["threshold"].first()
        assert thresholds["P2"] > thresholds["P1"]

    def test_a_patient_without_a_reference_is_not_called(self):
        """Never guessed. A cell with no baseline has no malignancy status."""
        df = cohort(n_reference=5)
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        assert (calls["call"].astype(str) == "not_called").all()
        assert calls["confidence"].isna().all()

    def test_length_mismatch_raises(self):
        df = cohort()
        with pytest.raises(MalignancyError, match="lengths differ"):
            call_malignancy(df["cnv_score"][:10], compartment=df["compartment"],
                            patient_id=df["patient_id"])

    def test_empty_input_raises(self):
        with pytest.raises(MalignancyError, match="no cells"):
            call_malignancy([], compartment=[], patient_id=[])


class TestNormalEpitheliumValidation:
    """§4's done-when: normal epithelium must not be misread as tumour."""

    def test_clean_calls_pass(self):
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        report = validate_normal_epithelium(calls, tissue=df["tissue"])
        assert bool(report.loc[0, "passed"])
        assert report.loc[0, "specificity"] > 0.95

    def test_normal_epithelium_miscalled_as_tumour_fails(self):
        """If this passed silently, every downstream number would be computed
        over a tumour arm contaminated with normal cells."""
        df = cohort()
        normal_epi = (df["tissue"] == "normal") & (df["compartment"] == "epithelial")
        df.loc[normal_epi, "cnv_score"] += 3.0        # make them look aneuploid
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        report = validate_normal_epithelium(calls, tissue=df["tissue"])
        assert not bool(report.loc[0, "passed"])

    def test_reference_cells_are_excluded_from_the_check(self):
        """They are in normal tissue but are the baseline, so counting them
        would inflate specificity toward 1 automatically."""
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        report = validate_normal_epithelium(calls, tissue=df["tissue"])
        assert int(report.loc[0, "n_normal_epithelial"]) == 300   # not 500

    def test_per_patient_rows(self):
        df = pd.concat([cohort(patient="P1"), cohort(patient="P2")], ignore_index=True)
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        report = validate_normal_epithelium(calls, tissue=df["tissue"])
        assert set(report["patient_id"]) == {"P1", "P2"}

    def test_no_normal_epithelium_raises_rather_than_passing_vacuously(self):
        df = cohort(n_normal_epi=0)
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        with pytest.raises(MalignancyError, match="nothing to validate"):
            validate_normal_epithelium(calls, tissue=df["tissue"])

    def test_tissue_length_mismatch_raises(self):
        df = cohort()
        calls = call_malignancy(
            df["cnv_score"], compartment=df["compartment"], patient_id=df["patient_id"]
        )
        with pytest.raises(MalignancyError, match="entries for"):
            validate_normal_epithelium(calls, tissue=["normal"])


class TestInferCNVWiring:
    """The two ways this can be quietly circular, plus the hg19 trap.

    inferCNV itself is not exercised here — it is R, it is slow, and it is not
    a pip dependency. What IS testable is everything that decides whether its
    answer means anything: which cells become the baseline, which are held back
    to score it, and whether the coordinates match the deposit's build.
    """

    def _roles(self, n_ref=120, n_hold=40, n_query=100):
        role = (
            ["reference_normal_epi"] * n_ref
            + ["holdout_normal_epi"] * n_hold
            + ["query"] * n_query
            + ["reference_diploid"] * 60
        )
        compartment = ["epithelial"] * (n_ref + n_hold + n_query) + (
            ["immune"] * 30 + ["stromal"] * 30
        )
        return pd.DataFrame({"role": role, "compartment": compartment})

    def _counts(self, roles, n_genes=25):
        return RNG.poisson(3, size=(len(roles), n_genes)).astype(np.int64)

    def test_the_holdout_is_never_a_reference_group(self):
        """It exists to be scored out-of-sample. As a reference it would be
        validating the baseline against itself."""
        groups = infercnv_reference_groups(self._roles())
        assert "holdout_normal_epi" not in groups
        assert "query" not in groups

    def test_matched_normal_epithelium_is_the_ONLY_reference(self):
        """One group, not four. inferCNV's STEP 08 runs with use_bounds=TRUE,
        which zeroes observation deviation falling inside the range of the
        reference-group means — and immune/stromal/endothelial/epithelial means
        differ for ordinary cell-type reasons. On the pilot that made 25-30% of
        values exactly 1 and inverted the ordering in four of five patients."""
        groups = infercnv_reference_groups(self._roles())
        assert groups == ["reference_normal_epi"]

    def test_diploid_groups_are_the_fallback_when_there_is_no_matched_normal(self):
        """A mismatched reference is the honest cost of having no better one.
        assign_cnv_roles calls this the diploid_only strategy and flags it."""
        roles = pd.DataFrame({
            "role": ["query"] * 50 + ["reference_diploid"] * 60,
            "compartment": ["epithelial"] * 50 + ["immune"] * 30 + ["stromal"] * 30,
        })
        groups = infercnv_reference_groups(roles)
        assert groups == ["ref_immune", "ref_stromal"]

    def test_the_diploid_fallback_still_keeps_compartments_separate(self):
        roles = pd.DataFrame({
            "role": ["reference_diploid"] * 3,
            "compartment": ["immune", "stromal", "endothelial"],
        })
        assert infercnv_reference_groups(roles) == [
            "ref_endothelial", "ref_immune", "ref_stromal"
        ]

    def test_no_reference_refuses_rather_than_running(self):
        roles = pd.DataFrame({
            "role": ["query"] * 10, "compartment": ["epithelial"] * 10
        })
        with pytest.raises(MalignancyError, match="no reference cells"):
            infercnv_reference_groups(roles)

    def test_a_grch38_gene_file_is_refused_by_absence(self, tmp_path):
        """Required with no default, and the error says why: this deposit is
        hg19, and a GRCh38 file produces arm-level artifacts that look exactly
        like real CNVs."""
        roles = self._roles()
        with pytest.raises(MalignancyError, match="hg19"):
            run_infercnv(
                self._counts(roles), [f"G{i}" for i in range(25)], roles,
                gene_position_file=tmp_path / "missing.txt",
                out_dir=tmp_path / "out", dry_run=True,
            )

    def test_dry_run_writes_inputs_and_returns_the_command(self, tmp_path):
        roles = self._roles()
        positions = tmp_path / "hg19_gene_pos.txt"
        positions.write_text("G0\tchr1\t1\t100\n")
        out = run_infercnv(
            self._counts(roles), [f"G{i}" for i in range(25)], roles,
            gene_position_file=positions, out_dir=tmp_path / "out", dry_run=True,
        )
        assert out["ran"] is False
        assert out["command"][0] == "Rscript"
        assert out["counts"].exists() and out["annotations"].exists()
        assert out["script"].exists()

    def test_the_generated_R_names_only_reference_groups(self, tmp_path):
        roles = self._roles()
        positions = tmp_path / "hg19_gene_pos.txt"
        positions.write_text("G0\tchr1\t1\t100\n")
        out = run_infercnv(
            self._counts(roles), [f"G{i}" for i in range(25)], roles,
            gene_position_file=positions, out_dir=tmp_path / "out", dry_run=True,
        )
        script = out["script"].read_text()
        assert "reference_normal_epi" in script
        assert "holdout_normal_epi" not in script
        # One reference group when a matched normal exists — passing the
        # diploid compartments alongside it is what bounded the signal away.
        assert "ref_immune" not in script

    def test_the_matrix_is_genes_by_cells(self, tmp_path):
        """inferCNV's orientation, not AnnData's. Transposed silently, this
        would run and return nonsense."""
        from scipy import io as sio

        roles = self._roles()
        paths = write_infercnv_inputs(
            self._counts(roles), [f"G{i}" for i in range(25)], roles,
            out_dir=tmp_path,
        )
        matrix = sio.mmread(paths["matrix"])
        assert matrix.shape == (25, len(roles))

    def test_the_matrix_is_written_sparse(self, tmp_path):
        """Dense would be 7.8 GB for the largest pilot patient before pandas
        takes its copy, and hours to write for a >90%-zero matrix."""
        roles = self._roles()
        paths = write_infercnv_inputs(
            self._counts(roles), [f"G{i}" for i in range(25)], roles,
            out_dir=tmp_path,
        )
        assert paths["matrix"].name.endswith(".mtx")
        assert (tmp_path / "genes.tsv").exists()
        assert (tmp_path / "barcodes.tsv").exists()
        assert paths["counts"] == tmp_path

    def test_the_R_script_reads_the_mtx_rather_than_a_path(self, tmp_path):
        """CreateInfercnvObject does read.table() on whatever path it gets and
        has no 10x-directory reader — a directory fails with 'not a regular
        file'. The matrix is loaded in R and passed as an object."""
        roles = self._roles()
        positions = tmp_path / "hg19.txt"
        positions.write_text("G0\tchr1\t1\t100\n")
        out = run_infercnv(
            self._counts(roles), [f"G{i}" for i in range(25)], roles,
            gene_position_file=positions, out_dir=tmp_path / "out", dry_run=True,
        )
        script = out["script"].read_text()
        assert "readMM" in script
        assert "raw_counts_matrix = counts" in script

    def test_duplicate_gene_symbols_are_collapsed(self, tmp_path):
        """The deposit maps several Ensembl IDs to one symbol, and the
        gene-order file is keyed on symbol — a duplicated row name gives
        inferCNV two positions for one gene."""
        from scipy import io as sio

        roles = self._roles()
        names = ["A", "B", "A", "C", "B"] + [f"G{i}" for i in range(20)]
        paths = write_infercnv_inputs(
            self._counts(roles), names, roles, out_dir=tmp_path,
        )
        written = (tmp_path / "genes.tsv").read_text().splitlines()
        symbols = [line.split("\t")[0] for line in written]
        assert len(symbols) == len(set(symbols))
        assert sio.mmread(paths["matrix"]).shape[0] == len(symbols)

    def test_genes_tsv_carries_the_symbol_in_both_columns(self, tmp_path):
        """The gene-order file is keyed on SYMBOL. A mismatch there does not
        error — it silently drops every gene from the inference."""
        roles = self._roles()
        write_infercnv_inputs(
            self._counts(roles), [f"G{i}" for i in range(25)], roles,
            out_dir=tmp_path,
        )
        first = (tmp_path / "genes.tsv").read_text().splitlines()[0]
        assert first.split("\t") == ["G0", "G0"]

    def test_unusable_cells_are_dropped_not_written(self, tmp_path):
        roles = self._roles()
        roles.loc[:9, "role"] = "unusable"
        paths = write_infercnv_inputs(
            self._counts(roles), [f"G{i}" for i in range(25)], roles,
            out_dir=tmp_path,
        )
        annotations = pd.read_csv(paths["annotations"], sep="\t", header=None)
        assert "unusable" not in set(annotations[1])

    def test_the_R_writes_the_score_itself(self, tmp_path):
        """no_plot=TRUE skips the step that writes infercnv.observations.txt,
        so the score comes off the final object instead of a 400 MB text
        matrix produced only to be reduced to one number per cell."""
        roles = self._roles()
        positions = tmp_path / "hg19.txt"
        positions.write_text("G0\tchr1\t1\t100\n")
        out = run_infercnv(
            self._counts(roles), [f"G{i}" for i in range(25)], roles,
            gene_position_file=positions, out_dir=tmp_path / "out", dry_run=True,
        )
        script = out["script"].read_text()
        assert "colMeans((expr - 1)^2)" in script
        assert "cnv_scores.csv" in script

    def test_the_score_table_joins_the_role_back_on(self, tmp_path):
        """call_malignancy needs the REFERENCE cells too — they set the
        threshold — so the group has to come back with the scores."""
        (tmp_path / "cnv_scores.csv").write_text(
            "cell,cnv_score\nc0,0.01\nc1,0.40\n"
        )
        (tmp_path / "annotations.tsv").write_text(
            "c0\tref_immune\nc1\tquery\n"
        )
        table = read_infercnv_score_table(tmp_path)
        assert set(table["group"]) == {"ref_immune", "query"}
        assert table.loc[table["cell"] == "c1", "cnv_score"].iloc[0] == 0.40

    def test_a_missing_score_file_says_where_to_look(self, tmp_path):
        with pytest.raises(MalignancyError, match="infercnv_R.log"):
            read_infercnv_score_table(tmp_path)

    def test_cleanup_refuses_when_the_result_is_missing(self, tmp_path):
        """A failed run keeps everything — the intermediates are how it gets
        diagnosed, and cnv_scores.csv is the only thing that cannot be
        recomputed without re-running the inference."""
        (tmp_path / "01_incoming_data.infercnv_obj").write_bytes(b"x" * 1000)
        freed = cleanup_infercnv_run(tmp_path)
        assert freed == 0
        assert (tmp_path / "01_incoming_data.infercnv_obj").exists()

    def test_cleanup_keeps_the_result_and_the_provenance(self, tmp_path):
        for name in ("cnv_scores.csv", "annotations.tsv", "run_infercnv.R",
                     "infercnv_R.log", "genes.tsv", "barcodes.tsv"):
            (tmp_path / name).write_text("x")
        for name in ("01_incoming_data.infercnv_obj", "22_denoise.infercnv_obj",
                     "preliminary.infercnv_obj", "matrix.mtx"):
            (tmp_path / name).write_bytes(b"x" * 1000)

        freed = cleanup_infercnv_run(tmp_path)
        assert freed == 4000
        assert (tmp_path / "cnv_scores.csv").exists()
        assert (tmp_path / "annotations.tsv").exists()
        assert (tmp_path / "run_infercnv.R").exists()
        assert not (tmp_path / "01_incoming_data.infercnv_obj").exists()
        assert not (tmp_path / "matrix.mtx").exists()

    def test_keep_final_is_opt_in(self, tmp_path):
        """Several hundred MB for the largest patient, times 62."""
        (tmp_path / "cnv_scores.csv").write_text("x")
        (tmp_path / "run.final.infercnv_obj").write_bytes(b"x" * 500)

        cleanup_infercnv_run(tmp_path, keep_final=True)
        assert (tmp_path / "run.final.infercnv_obj").exists()
        cleanup_infercnv_run(tmp_path)
        assert not (tmp_path / "run.final.infercnv_obj").exists()

    def test_scores_feed_call_malignancy(self, tmp_path):
        """Mean squared deviation from 1, and it must land on the scale
        call_malignancy thresholds against."""
        observations = tmp_path / "obs.txt"
        genes, cells = 10, 6
        frame = pd.DataFrame(
            np.full((genes, cells), 1.3),
            index=[f"G{i}" for i in range(genes)],
            columns=[f"c{i}" for i in range(cells)],
        )
        frame.to_csv(observations, sep=" ")
        scores = read_infercnv_scores(observations)
        assert len(scores) == cells
        assert scores.iloc[0] == pytest.approx(0.09)


def test_min_reference_cells_is_documented_not_arbitrary():
    assert MIN_REFERENCE_CELLS >= 20


class TestHeldOutReference:
    """The corrected design: matched normal epithelium as baseline, part held out.

    The first version swapped the reference to immune/stromal so the validation
    would be non-circular. That fixed validation by breaking inference — the
    inferCNV docs are explicit that using a different cell type as reference
    makes the method read cell-type differences as copy number, because
    co-regulated genes cluster on chromosomes. Holding out part of the matched
    normal keeps the right baseline AND gives an out-of-sample check.
    """

    def _mixed(self, n_normal_epi=300, n_diploid=200, n_tumour_epi=400, patient="P1"):
        return pd.DataFrame(
            {
                "compartment": (
                    ["epithelial"] * n_normal_epi
                    + ["immune"] * n_diploid
                    + ["epithelial"] * n_tumour_epi
                ),
                "tissue": (
                    ["normal"] * (n_normal_epi + n_diploid) + ["tumour"] * n_tumour_epi
                ),
                "patient_id": [patient] * (n_normal_epi + n_diploid + n_tumour_epi),
            }
        )

    def test_matched_normal_is_the_baseline_not_the_diploid_cells(self):
        df = self._mixed()
        roles, report = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"]
        )
        assert report.loc[0, "strategy"] == "matched_normal"
        assert (roles["role"] == "reference_normal_epi").sum() > 0

    def test_part_of_the_normal_epithelium_is_held_out(self):
        df = self._mixed(n_normal_epi=300)
        roles, _ = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"],
            holdout_fraction=0.30,
        )
        held = int((roles["role"] == "holdout_normal_epi").sum())
        used = int((roles["role"] == "reference_normal_epi").sum())
        assert held + used == 300
        assert held == pytest.approx(90, abs=2)

    def test_held_out_cells_are_disjoint_from_the_baseline(self):
        """The whole point. A cell in both would make the check circular again."""
        df = self._mixed()
        roles, _ = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"]
        )
        assert not (
            (roles["role"] == "holdout_normal_epi")
            & (roles["role"] == "reference_normal_epi")
        ).any()

    def test_diploid_cells_are_kept_as_a_separate_category(self):
        """Not merged into the baseline: inferCNV bounds the log fold change by
        per-category means, which is what suppresses cell-type false positives."""
        df = self._mixed()
        roles, _ = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"]
        )
        assert (roles["role"] == "reference_diploid").sum() == 200

    def test_tumour_epithelium_is_the_query(self):
        df = self._mixed()
        roles, _ = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"]
        )
        assert (roles["role"] == "query").sum() == 400

    def test_a_patient_without_matched_normal_falls_back_and_is_flagged(self):
        """26 of 62 patients on this cohort (decision #9). Their calls come from
        a weaker method and must not be pooled silently."""
        df = self._mixed(n_normal_epi=0)
        roles, report = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"]
        )
        assert report.loc[0, "strategy"] == "diploid_only"
        assert (roles["role"] == "query").sum() == 400
        assert (roles["role"] == "holdout_normal_epi").sum() == 0

    def test_a_patient_with_neither_is_unusable(self):
        df = self._mixed(n_normal_epi=0, n_diploid=5)
        roles, report = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"]
        )
        assert report.loc[0, "strategy"] == "none"
        assert (roles["role"] == "unusable").all()

    def test_assignment_is_deterministic(self):
        df = self._mixed()
        a, _ = assign_cnv_roles(df["compartment"], tissue=df["tissue"],
                                patient_id=df["patient_id"], seed=5)
        b, _ = assign_cnv_roles(df["compartment"], tissue=df["tissue"],
                                patient_id=df["patient_id"], seed=5)
        pd.testing.assert_frame_equal(a, b)

    def test_strategies_are_decided_per_patient(self):
        df = pd.concat(
            [self._mixed(patient="P1"), self._mixed(n_normal_epi=0, patient="P2")],
            ignore_index=True,
        )
        _, report = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"]
        )
        strategies = report.set_index("patient_id")["strategy"]
        assert strategies["P1"] == "matched_normal"
        assert strategies["P2"] == "diploid_only"

    def test_validation_uses_only_the_held_out_cells(self):
        df = self._mixed(n_normal_epi=300, n_diploid=200, n_tumour_epi=400)
        roles, _ = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"],
            holdout_fraction=0.30,
        )
        scores = np.where(df["tissue"] == "tumour",
                          RNG.normal(1.5, 0.15, len(df)),
                          RNG.normal(1.0, 0.15, len(df)))
        calls = call_malignancy(
            scores, compartment=df["compartment"], patient_id=df["patient_id"]
        )
        report = validate_normal_epithelium(
            calls, tissue=df["tissue"], role=roles["role"]
        )
        # ~90 held out, not all 300 normal epithelial cells.
        assert int(report.loc[0, "n_normal_epithelial"]) == pytest.approx(90, abs=5)

    def test_validation_refuses_when_roles_show_nothing_was_held_out(self):
        """Silently validating on baseline cells is the failure this exists to
        prevent."""
        df = self._mixed(n_normal_epi=0)
        roles, _ = assign_cnv_roles(
            df["compartment"], tissue=df["tissue"], patient_id=df["patient_id"]
        )
        scores = RNG.normal(1.0, 0.15, len(df))
        calls = call_malignancy(
            scores, compartment=df["compartment"], patient_id=df["patient_id"]
        )
        with pytest.raises(MalignancyError, match="no held-out normal epithelium"):
            validate_normal_epithelium(calls, tissue=df["tissue"], role=roles["role"])


class TestGenePositions:
    """The gene-order file decides where inferCNV thinks every gene is. Getting
    it wrong does not fail — it produces chromosome-arm artifacts that look
    exactly like real copy-number events."""

    def _parse(self, text):
        import io
        import sys
        from pathlib import Path

        sys.path.insert(0, str(Path("src/reference/jobs").resolve()))
        from fetch_gene_positions import parse_gtf, sort_rows

        return sort_rows(parse_gtf(io.StringIO(text)))

    def _gene(self, chrom, start, end, name):
        return (
            f'{chrom}\tHAVANA\tgene\t{start}\t{end}\t.\t+\t.\t'
            f'gene_name "{name}"; gene_type "protein_coding";\n'
        )

    def test_only_gene_rows_are_kept(self):
        text = self._gene("chr1", 100, 200, "A") + (
            'chr1\tHAVANA\texon\t100\t200\t.\t+\t.\tgene_name "A";\n'
        )
        assert len(self._parse(text)) == 1

    def test_scaffolds_and_chrM_are_excluded(self):
        """A contig with a handful of genes contributes noise with no positional
        meaning, and chrM has no copy number in the relevant sense."""
        text = (
            self._gene("chr1", 100, 200, "A")
            + self._gene("chrM", 1, 10, "MT-CO1")
            + self._gene("GL000191.1", 1, 10, "SCAF")
        )
        assert [r[0] for r in self._parse(text)] == ["A"]

    def test_pseudoautosomal_duplicates_are_dropped(self):
        """PAR genes appear on both X and Y. A gene at two positions makes the
        smoothing window ambiguous."""
        text = self._gene("chr1", 100, 200, "A") + self._gene(
            "chrY", 5, 9, "XG_PAR_Y"
        )
        assert [r[0] for r in self._parse(text)] == ["A"]

    def test_a_repeated_symbol_keeps_the_first_occurrence(self):
        text = self._gene("chr1", 100, 200, "A") + self._gene("chr1", 900, 950, "A")
        rows = self._parse(text)
        assert len(rows) == 1
        assert rows[0][2] == 100

    def test_output_is_in_genomic_order(self):
        """inferCNV walks the file top to bottom as the genome."""
        text = (
            self._gene("chr2", 50, 80, "B")
            + self._gene("chr1", 900, 950, "C")
            + self._gene("chr1", 100, 200, "A")
            + self._gene("chrX", 10, 20, "D")
        )
        assert [r[0] for r in self._parse(text)] == ["A", "C", "B", "D"]
