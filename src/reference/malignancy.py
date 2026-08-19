"""Malignant vs. normal epithelium. W1, weeks 2-3.

Why this matters more than it looks
-----------------------------------
Until it exists, the "tumour" arm of the decomposition is *tumour-sample
epithelium*, which contains a great deal of non-malignant epithelium. The
contrast is then sample-of-origin, not malignant-versus-normal — a different
claim from the one the project makes, and one that would quietly attribute
normal cells sitting in a tumour to the tumour.

Choice of reference, and why it departs from the plan's wording
---------------------------------------------------------------
execution_plan.md §4 says to use *matched normal epithelium* as the inferCNV
reference. The same row's "done when" is that **normal epithelium is not misread
as tumour** — and those two cannot both hold, because a population used as the
CNV baseline is non-malignant by construction. The check would validate nothing.

So the reference here is the **non-epithelial compartments** — immune, stromal
and endothelial — which are reliably diploid and are not the population under
test. Malignancy is then called on *all* epithelium including the normal
samples, and normal-sample epithelium coming back non-malignant becomes a real,
non-circular validation. :func:`validate_normal_epithelium` is that check.

Per patient, never pooled: CNV inference compares a cell against a baseline, and
a baseline built across patients would fold germline copy-number variation and
per-patient capture differences into the malignancy call.

What is implemented, and what is not
------------------------------------
The reference selection, the call from CNV scores, the confidence, and the
validation are here and tested. Running inferCNV itself is an R call against real
matrices whose settings need tuning on real data — stubbed, the same pattern as
``_select_markers`` and ``flag_doublets``. ``src/reference/jobs/infercnv.sh``
drives it per patient as an SGE array job.
"""

from __future__ import annotations

from typing import Any, Final

import numpy as np
import pandas as pd

#: Compartments used as the diploid CNV baseline. Not epithelium, so the
#: "normal epithelium is not misread as tumour" check stays non-circular.
REFERENCE_COMPARTMENTS: Final[tuple[str, ...]] = ("immune", "stromal", "endothelial")

#: A baseline built from fewer cells than this is too noisy to call against.
MIN_REFERENCE_CELLS: Final[int] = 50

#: Quantile of the reference CNV score above which a cell is called malignant.
#: Deliberately strict: a false malignant call moves a normal cell into the
#: tumour arm, which is the direction that inflates the compositional term.
MALIGNANT_QUANTILE: Final[float] = 0.99

#: Below this share of normal-sample epithelium called non-malignant, the calls
#: are not trustworthy. Normal samples do contain some aneuploid cells, so the
#: bar is not 100%.
MIN_NORMAL_SPECIFICITY: Final[float] = 0.90

CALLS: Final[tuple[str, ...]] = ("malignant", "non_malignant", "reference", "not_called")


class MalignancyError(ValueError):
    """Malignancy could not be called from the input given."""


def select_cnv_reference(
    compartment: Any,
    *,
    patient_id: Any,
    min_cells: int = MIN_REFERENCE_CELLS,
) -> pd.DataFrame:
    """Per-patient reference availability. Run before inferCNV, not after.

    Returns one row per patient with the reference cell count and whether it
    clears `min_cells`. A patient without enough diploid cells cannot have
    malignancy called at all, and that is a cohort fact worth knowing before
    hours of CNV inference rather than after.
    """
    frame = pd.DataFrame(
        {
            "patient_id": [str(p) for p in patient_id],
            "compartment": [str(c) for c in compartment],
        }
    )
    frame["is_reference"] = frame["compartment"].isin(REFERENCE_COMPARTMENTS)
    out = (
        frame.groupby("patient_id", observed=True)
        .agg(
            n_cells=("compartment", "size"),
            n_reference=("is_reference", "sum"),
            n_epithelial=("compartment", lambda s: int((s == "epithelial").sum())),
        )
        .reset_index()
    )
    out["usable"] = out["n_reference"] >= min_cells
    return out


def call_malignancy(
    cnv_score: Any,
    *,
    compartment: Any,
    patient_id: Any,
    quantile: float = MALIGNANT_QUANTILE,
    min_cells: int = MIN_REFERENCE_CELLS,
) -> pd.DataFrame:
    """Per-cell malignancy call **with confidence**, thresholded per patient.

    `cnv_score` is whatever inferCNV or CopyKAT emits per cell — conventionally
    the mean squared deviation of the smoothed CNV profile from the reference.
    Higher means more aneuploid.

    The threshold is the `quantile` of the *reference* cells' own scores within
    that patient, so it adapts to how noisy that patient's data is rather than
    being a number carried between datasets. `confidence` is the score's distance
    above the threshold in reference standard deviations — a cell just over the
    line is not the same as one far past it, and the decomposition should be able
    to tell.

    Cells in a patient without a usable reference are ``not_called``, never
    guessed. Reference cells themselves are labelled ``reference`` rather than
    ``non_malignant``: they defined the baseline, so calling them non-malignant
    would be circular.
    """
    scores = np.asarray(cnv_score, dtype=float)
    comp = np.asarray([str(c) for c in compartment], dtype=object)
    patients = np.asarray([str(p) for p in patient_id], dtype=object)
    if not (scores.shape[0] == comp.shape[0] == patients.shape[0]):
        raise MalignancyError(
            f"lengths differ: cnv_score {scores.shape[0]}, compartment "
            f"{comp.shape[0]}, patient_id {patients.shape[0]}"
        )
    if scores.size == 0:
        raise MalignancyError("no cells to call")

    is_reference = np.isin(comp, REFERENCE_COMPARTMENTS)
    call = np.full(scores.shape, "not_called", dtype=object)
    confidence = np.full(scores.shape, np.nan, dtype=float)
    threshold = np.full(scores.shape, np.nan, dtype=float)

    for patient in pd.unique(patients):
        here = patients == patient
        reference = here & is_reference
        if int(reference.sum()) < min_cells:
            continue
        cut = float(np.quantile(scores[reference], quantile))
        spread = float(np.std(scores[reference])) or 1.0
        epithelial = here & (comp == "epithelial")

        threshold[here] = cut
        confidence[epithelial] = (scores[epithelial] - cut) / spread
        call[epithelial] = np.where(scores[epithelial] > cut, "malignant", "non_malignant")
        call[reference] = "reference"

    return pd.DataFrame(
        {
            "patient_id": patients,
            "compartment": comp,
            "cnv_score": scores,
            "threshold": threshold,
            "confidence": confidence,
            "call": pd.Categorical(call, categories=CALLS),
        }
    )


def validate_normal_epithelium(
    calls: pd.DataFrame, *, tissue: Any, min_specificity: float = MIN_NORMAL_SPECIFICITY
) -> pd.DataFrame:
    """**The check that makes the calls believable.** Per patient.

    Epithelium from a patient's *normal* sample should come back overwhelmingly
    non-malignant. Because the CNV baseline is non-epithelial, this is a genuine
    out-of-sample test rather than a restatement of the reference.

    execution_plan.md §4 lists it as the "done when" for this stage. If it fails,
    stop — every downstream compositional and intrinsic number would be computed
    over a tumour arm contaminated with normal cells, or a normal arm stripped of
    them.
    """
    tissue_arr = np.asarray([str(t) for t in tissue], dtype=object)
    if tissue_arr.shape[0] != len(calls):
        raise MalignancyError(f"tissue has {tissue_arr.shape[0]} entries for {len(calls)} cells")

    frame = calls.assign(tissue=tissue_arr)
    normal = frame[
        (frame["tissue"] == "normal")
        & frame["call"].astype(str).isin(["malignant", "non_malignant"])
    ]
    if normal.empty:
        raise MalignancyError(
            "no called epithelium in normal samples — nothing to validate against"
        )

    out = (
        normal.assign(ok=normal["call"].astype(str) == "non_malignant")
        .groupby("patient_id", observed=True)
        .agg(n_normal_epithelial=("ok", "size"), n_non_malignant=("ok", "sum"))
        .reset_index()
    )
    out["specificity"] = out["n_non_malignant"] / out["n_normal_epithelial"]
    out["passed"] = out["specificity"] >= min_specificity
    return out


def run_infercnv(*args: Any, **kwargs: Any) -> Any:
    """inferCNV against real matrices. W1, weeks 2-3. Not scaffolded further.

    Needs the actual per-patient matrices, and its window size, cutoff and HMM
    settings are judgement calls over real data rather than a formula — the same
    reason ``_select_markers`` and ``flag_doublets`` are stubs.

    Drive it with ``src/reference/jobs/infercnv.sh``, which submits one array
    task per patient. Pass the non-epithelial compartments as the reference
    group, not the matched normal epithelium: see this module's docstring for
    why the plan's wording cannot be taken literally without making the
    validation circular.

    Cross-check the result with CopyKAT and report the concordance — §4 asks for
    the cross-check, not a winner.
    """
    raise NotImplementedError(
        "W1 — inferCNV needs the real per-patient matrices. "
        "See src/reference/jobs/infercnv.sh."
    )
