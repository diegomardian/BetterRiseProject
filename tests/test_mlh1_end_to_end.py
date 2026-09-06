"""The whole job, against a synthetic atlas whose MLH1 truth we set.

WHY THIS EXISTS. Every other test here checks a decision the job makes once the
cells are scored. None of them touches the path that does the scoring --
`read_var` -> `read_cells` -> QC -> `assign_labels` -> depth matching ->
`rows_for_patient` -> the interval -> the verdict. That path runs on a 30 GB
atlas on a cluster, and until this file existed, the first execution of it would
have been the real one.

Two runs, and the point is that they must come out DIFFERENT:

* an atlas where MLH1 is thinned to 0.15x in the tumour arm of methylated
  patients and nowhere else -- the job must find it, in the right arm, in the
  right direction;
* an atlas with no MLH1 effect anywhere -- the job must not.

A pipeline that returned "silencing detected" on both would pass any test that
only ran the first, and that is the failure mode this repository is named after.
The synthetic atlas is built in the real one's shape: symbols in
`var/GeneSymbol` with Ensembl in `_index`, raw counts in `layers/counts` as CSR,
per-cell lognormal size factors so the MAD thresholds have real spread to work
with.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

h5py = pytest.importorskip("h5py")
sparse = pytest.importorskip("scipy.sparse")

from src.reference.jobs.coexpression_silencing import premise_holds  # noqa: E402
from src.reference.jobs.icbi_coexpression import load_obs, study_deltas  # noqa: E402
from src.reference.jobs.mlh1_positive_control import (  # noqa: E402
    PRIMARY_STRATUM,
    SECONDARY_STRATUM,
    TARGET,
    arm_reading,
    instrument_verdict,
    strata_for,
)

AXIS = ["LGR5", "ASCL2", "MKI67", "OLFM4", "SMOC2"]
TIER_A = ["GUCA2A", "GUCA2B", "OTOP2", "CA7"]
PANEL = ["ACTB", "KRT8", "EPCAM", "CDX2", "MS4A12"]
BEST4 = ["BEST4", "CFTR", "HES4", "SPIB"]
GOBLET = ["MUC2", "TFF3", "SPDEF", "ITLN1"]
OTHER = [TARGET, "MT-CO1", "MT-ND1", "RPL13A", "RPS6"]
FILLER = [f"FILL{i}" for i in range(320)]
GENES = AXIS + TIER_A + PANEL + BEST4 + GOBLET + OTHER + FILLER

# Real Pelka ids, so the join to the committed cohort table resolves.
METHYLATED = ["C106", "C109", "C110", "C111", "C114"]
UNMETHYLATED = ["C103", "C104", "C105", "C107", "C112"]

#: MLH1's per-cell Poisson rate in the normal arm.
#:
#: DELIBERATELY HIGHER THAN THE REAL COHORT'S ~0.032 (~3% detection), and the
#: reason matters. At the real rate a patient carries ~8 MLH1+ cells per arm,
#: and reproducing that here would mean either a very large synthetic cohort or
#: a test that fails a fifth of the time -- because the real design is powered
#: at 73.7% for a 50% effect, and a unit test may not be.
#:
#: So this test is NOT a power measurement and must not be read as one. Power
#: is answered by `results/*/mlh1_power.parquet`, which simulates the real
#: cell counts and depths thousands of times and is the right instrument for
#: it. What this file checks is the PLUMBING: that the scoring path runs, that
#: an effect put into one arm comes out of that arm with the right sign, and
#: that no effect comes out as no effect. Conflating the two would give a
#: flaky test and a power figure that rests on one draw.
MLH1_RATE = 0.15
N_CELLS = 420


def _rates(mature: bool, arm: str, patient: str, is_best4: bool,
           silencing: float) -> dict[str, float]:
    r = dict.fromkeys(FILLER, 12.0)
    r |= {"ACTB": 28.0, "KRT8": 22.0, "RPL13A": 18.0, "RPS6": 15.0,
          "MT-CO1": 9.0, "MT-ND1": 7.0, "EPCAM": 12.0, "CDX2": 1.2}
    for g in AXIS:
        r[g] = 0.4 if mature else 6.0
    for g in TIER_A + ["MS4A12"]:
        r[g] = (14.0 if mature else 0.5) * (0.05 if arm == "tumour" else 1.0)
    for g in BEST4:
        r[g] = 8.0 if (mature and is_best4) else 0.3
    for g in GOBLET:
        r[g] = 5.0 if mature else 1.0
    rate = MLH1_RATE
    if arm == "tumour" and patient in METHYLATED:
        rate *= silencing
    r[TARGET] = rate
    return r


def _synthetic_atlas(tmp_path, silencing: float, seed: int = 7):
    """Write an atlas + obs cache in the ICBI atlas's shape. Returns both paths."""
    rng = np.random.default_rng(seed)
    rows, blocks = [], []
    for patient in METHYLATED + UNMETHYLATED:
        for arm, sample_type in (("normal", "adjacent normal"),
                                 ("tumour", "primary tumor")):
            for i in range(N_CELLS // 2):
                r = _rates(i % 2 == 0, arm, patient, i % 8 == 0, silencing)
                size = float(rng.lognormal(0.0, 0.35))
                blocks.append(
                    rng.poisson(np.array([r[g] for g in GENES]) * size)
                    .astype(np.float32)
                )
                rows.append({
                    "study_id": "Pelka_2021_Cell", "dataset": "synthetic",
                    "sample_id": f"{patient}_{arm}",
                    "patient_id": f"Pelka_2021_Cell.{patient}",
                    "sample_type": sample_type, "platform": "10x",
                    "platform_fine": "10x_3p_v2",
                    "enrichment_cell_types": "naive",
                    "matrix_type": "raw counts", "reference_genome": "gencode32",
                    "suspension_type": "cell", "tissue_cell_state": "fresh",
                    "atlas_cell_type_coarse": (
                        "Epithelial cell" if i % 10 != 9 else "T cell"),
                    "n_genes": 0, "total_counts": 0.0, "pct_counts_mito": 0.05,
                    "SOLO_doublet_status": "singlet",
                    "MLH1_promoter_methylation_status": (
                        "meth" if patient in METHYLATED else "no_meth"),
                    "microsatellite_status": (
                        "MSI-H" if patient in METHYLATED else "MSS"),
                })

    matrix = sparse.csr_matrix(np.vstack(blocks))
    obs = pd.DataFrame(rows)
    obs["n_genes"] = np.diff(matrix.indptr)
    obs["total_counts"] = np.asarray(matrix.sum(axis=1)).ravel()

    atlas = tmp_path / f"atlas_{silencing}.h5ad"
    obs_path = tmp_path / f"obs_{silencing}.parquet"
    obs.to_parquet(obs_path, index=False)
    with h5py.File(str(atlas), "w") as h5:
        var = h5.create_group("var")
        var.attrs["_index"] = "_index"
        var.create_dataset("_index", data=np.array(
            [f"ENSG{i:08d}" for i in range(len(GENES))], dtype="S20"))
        var.create_dataset("GeneSymbol", data=np.array(GENES, dtype="S20"))
        layer = h5.create_group("layers/counts")
        layer.create_dataset("data", data=matrix.data.astype(np.float32))
        layer.create_dataset("indices", data=matrix.indices.astype(np.int32))
        layer.create_dataset("indptr", data=matrix.indptr.astype(np.int64))
        layer.attrs["shape"] = np.array(matrix.shape, dtype=np.int64)
    return atlas, obs_path


def _run(tmp_path, silencing: float) -> tuple[dict, tuple[bool, str], pd.DataFrame]:
    atlas, obs_path = _synthetic_atlas(tmp_path, silencing)
    obs = load_obs(obs_path, atlas)
    strata = strata_for(obs)
    deltas, _ = study_deltas(atlas, obs, "Pelka_2021_Cell", seed=1,
                            reading="carcinoma", extra_genes=(TARGET,))
    assert not deltas.empty, "the scoring path produced nothing"
    deltas["short_id"] = deltas["patient_id"].astype(str).str.split(".").str[-1]
    deltas = deltas.merge(strata.drop(columns=["patient_id"], errors="ignore"),
                          on="short_id", how="left")
    assert "patient_id" in deltas.columns

    primary = deltas[deltas["mlh1_stratum"] == "mlh1_methylated"]
    row = arm_reading(primary, gene=TARGET, arm=PRIMARY_STRATUM, seed=1)
    premise = premise_holds(primary, seed=1)
    return row, premise, deltas


@pytest.mark.slow
def test_the_job_recovers_a_silencing_effect_it_was_given(tmp_path):
    """0.15x in the methylated tumour arm. True log fold change ln(0.15) = -1.90."""
    row, premise, deltas = _run(tmp_path, silencing=0.15)
    assert premise[0], f"premise must hold on clean synthetic data: {premise[1]}"
    assert row["mean_delta_cloglog"] < 0
    assert row["excludes_zero"]
    # Attenuated toward zero by the boundary rule at ~8 positive cells; the
    # point is the direction and rough magnitude, not the third decimal.
    assert row["mean_delta_cloglog"] == pytest.approx(np.log(0.15), abs=0.6)
    assert instrument_verdict(row, premise)["verdict"] == \
        "INSTRUMENT SEES KNOWN SILENCING"

    # ... and NOT in the arm where nothing was done to it.
    other = deltas[deltas["mlh1_stratum"] != "mlh1_methylated"]
    control = arm_reading(other, gene=TARGET, arm=SECONDARY_STRATUM, seed=1)
    assert not control["excludes_zero"], (
        "the effect was injected into the methylated arm only; firing in the "
        "comparison arm would mean the strata are not what the job thinks"
    )


@pytest.mark.slow
def test_the_job_does_not_fire_when_there_is_nothing_to_find(tmp_path):
    """THE HALF THAT MAKES THE OTHER ONE MEAN SOMETHING.

    Same atlas, same genes, same everything -- MLH1 simply unthinned. A
    pipeline that reported silencing here would report it on the real data too,
    and the first test would not have noticed.
    """
    row, premise, _ = _run(tmp_path, silencing=1.0)
    assert premise[0], f"premise must hold on clean synthetic data: {premise[1]}"
    # Assert the reading is ESTIMABLE first. Otherwise "did not fire" could be
    # the not-estimable branch, and the test would pass by never having asked
    # the question -- which is the failure mode this file is guarding against.
    assert row["patients_with_signal"] == len(METHYLATED), (
        "every synthetic patient should carry MLH1+ cells at this rate; if not, "
        "the null test is passing through NOT ESTIMABLE without ever forming "
        "the comparison"
    )
    verdict = instrument_verdict(row, premise)["verdict"]
    assert verdict == "INSTRUMENT DOES NOT SEE IT", (
        f"expected the estimable no-effect branch, got {verdict}: "
        f"{row['mean_delta_cloglog']:+.3f} "
        f"[{row['ci_low']:+.3f}, {row['ci_high']:+.3f}]"
    )
