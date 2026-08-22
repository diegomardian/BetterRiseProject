"""Malignant vs. normal epithelium. W1, weeks 2-3.

Why this matters more than it looks
-----------------------------------
Until it exists, the "tumour" arm of the decomposition is *tumour-sample
epithelium*, which contains a great deal of non-malignant epithelium. The
contrast is then sample-of-origin, not malignant-versus-normal — a different
claim from the one the project makes, and one that would quietly attribute
normal cells sitting in a tumour to the tumour.

Choice of reference
-------------------
execution_plan.md §4 says to use *matched normal epithelium*, and the same row's
"done when" is that **normal epithelium is not misread as tumour**. Taken
literally those conflict: a population used as the CNV baseline is non-malignant
by construction, so validating on it proves nothing.

The resolution is not to change the reference but to **hold part of it out**.
Matched normal epithelium is the right baseline — inferCNV compares smoothed
expression along the genome, and because co-regulated genes cluster on
chromosomes, a reference of a *different cell type* produces spurious CNV. The
inferCNV documentation names this directly: there are acknowledged caveats when
using epithelium against immune cells as reference, because the method finds
cell-type differences and reads them as copy number. Immune/stromal references
are common in practice, but they are the fallback, not the better option.

So, per patient:

- a random **70%** of normal-sample epithelium becomes the CNV baseline;
- the held-out **30%** is scored like any query cell, and must come back
  non-malignant — genuinely out-of-sample, because it never defined the baseline;
- the non-epithelial compartments are supplied as **additional reference
  categories**, not merged into one. inferCNV bounds the log fold change by the
  per-category means, which is what suppresses cell-type-specific false
  positives; pooling them into a single reference throws that away.

Patients with no matched normal — 26 of 62 on this cohort (open decision #9) —
cannot do this. They fall back to a diploid-only reference and are flagged, since
their calls come from a different and weaker method and should not be pooled with
the rest without saying so.

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

#: Diploid compartments, supplied to inferCNV as SEPARATE reference categories.
#: Kept separate on purpose: inferCNV bounds the log fold change by the per-
#: category means, which is what suppresses cell-type-specific false positives.
#: Merging them into one reference discards that protection.
DIPLOID_COMPARTMENTS: Final[tuple[str, ...]] = ("immune", "stromal", "endothelial")

#: Backwards-compatible alias. Prefer DIPLOID_COMPARTMENTS.
REFERENCE_COMPARTMENTS: Final[tuple[str, ...]] = DIPLOID_COMPARTMENTS

#: Share of normal-sample epithelium held out of the baseline so the
#: "normal epithelium is not misread as tumour" check is out-of-sample.
HOLDOUT_FRACTION: Final[float] = 0.30

#: Per-patient reference strategy.
STRATEGIES: Final[tuple[str, ...]] = ("matched_normal", "diploid_only", "none")

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


def assign_cnv_roles(
    compartment: Any,
    *,
    tissue: Any,
    patient_id: Any,
    holdout_fraction: float = HOLDOUT_FRACTION,
    min_cells: int = MIN_REFERENCE_CELLS,
    seed: int = 20260101,
) -> pd.DataFrame:
    """Assign every cell a role for CNV inference. Run this before inferCNV.

    Roles:

    ``reference_normal_epi``
        Normal-sample epithelium forming the baseline. Cell-type matched to the
        query, which is what keeps cell-type differences out of the CNV signal.
    ``holdout_normal_epi``
        Normal-sample epithelium deliberately excluded from the baseline, so the
        validation is out-of-sample. This is the population
        :func:`validate_normal_epithelium` scores.
    ``reference_diploid``
        Immune/stromal/endothelial, passed as separate additional categories.
    ``query``
        Tumour-sample epithelium — the cells being called.
    ``unusable``
        Cells in a patient with no viable reference at all.

    Also returns, per patient, which `strategy` was used. A patient with enough
    matched normal epithelium gets ``matched_normal``; one without falls back to
    ``diploid_only`` and is flagged, because those calls come from a weaker
    method and must not be silently pooled with the rest.
    """
    comp = np.asarray([str(c) for c in compartment], dtype=object)
    tis = np.asarray([str(t) for t in tissue], dtype=object)
    pat = np.asarray([str(p) for p in patient_id], dtype=object)
    if not (comp.shape[0] == tis.shape[0] == pat.shape[0]):
        raise MalignancyError(
            f"lengths differ: compartment {comp.shape[0]}, tissue {tis.shape[0]}, "
            f"patient_id {pat.shape[0]}"
        )
    if not 0.0 < holdout_fraction < 1.0:
        raise MalignancyError(f"holdout_fraction must be in (0, 1), got {holdout_fraction}")

    rng = np.random.default_rng(seed)
    role = np.full(comp.shape, "unusable", dtype=object)
    epithelial = comp == "epithelial"
    diploid = np.isin(comp, DIPLOID_COMPARTMENTS)

    rows = []
    for patient in pd.unique(pat):
        here = pat == patient
        normal_epi = np.flatnonzero(here & epithelial & (tis == "normal"))
        n_diploid = int((here & diploid).sum())

        # Enough matched normal epithelium to both build a baseline and hold
        # some back? min_cells applies to the baseline AFTER the holdout.
        needed = int(np.ceil(min_cells / (1.0 - holdout_fraction)))
        if normal_epi.size >= needed:
            shuffled = rng.permutation(normal_epi)
            n_hold = int(round(normal_epi.size * holdout_fraction))
            role[shuffled[:n_hold]] = "holdout_normal_epi"
            role[shuffled[n_hold:]] = "reference_normal_epi"
            role[here & diploid] = "reference_diploid"
            role[here & epithelial & (tis != "normal")] = "query"
            strategy = "matched_normal"
        elif n_diploid >= min_cells:
            # No usable matched normal. Diploid-only reference, and every
            # epithelial cell becomes a query — including normal-sample ones,
            # which then serve as the (weaker) validation.
            role[here & diploid] = "reference_diploid"
            role[here & epithelial] = "query"
            strategy = "diploid_only"
        else:
            strategy = "none"

        rows.append(
            {
                "patient_id": patient,
                "strategy": strategy,
                "n_normal_epithelial": int(normal_epi.size),
                "n_diploid": n_diploid,
                "n_epithelial": int((here & epithelial).sum()),
                "usable": strategy != "none",
            }
        )

    return pd.DataFrame(
        {"patient_id": pat, "compartment": comp, "tissue": tis, "role": role}
    ), pd.DataFrame(rows)


def select_cnv_reference(
    compartment: Any,
    *,
    patient_id: Any,
    min_cells: int = MIN_REFERENCE_CELLS,
) -> pd.DataFrame:
    """Diploid reference availability per patient. Cheap pre-flight check.

    Kept for the case where only compartment labels are to hand. When tissue is
    available prefer :func:`assign_cnv_roles`, which uses matched normal
    epithelium as the baseline and holds part of it out.
    """
    frame = pd.DataFrame(
        {
            "patient_id": [str(p) for p in patient_id],
            "compartment": [str(c) for c in compartment],
        }
    )
    frame["is_reference"] = frame["compartment"].isin(DIPLOID_COMPARTMENTS)
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
    calls: pd.DataFrame,
    *,
    tissue: Any,
    role: Any = None,
    min_specificity: float = MIN_NORMAL_SPECIFICITY,
) -> pd.DataFrame:
    """**The check that makes the calls believable.** Per patient.

    Epithelium from a patient's *normal* sample should come back overwhelmingly
    non-malignant.

    Pass `role` from :func:`assign_cnv_roles` and the check runs on the
    **held-out** normal epithelium only — cells that never entered the baseline,
    so the test is genuinely out-of-sample. Without `role` it falls back to all
    normal epithelium, which is only valid when the baseline was diploid-only;
    under a matched-normal baseline that fallback is partly circular and the
    specificity it reports is optimistic.

    execution_plan.md §4 lists it as the "done when" for this stage. If it fails,
    stop — every downstream compositional and intrinsic number would be computed
    over a tumour arm contaminated with normal cells, or a normal arm stripped of
    them.
    """
    tissue_arr = np.asarray([str(t) for t in tissue], dtype=object)
    if tissue_arr.shape[0] != len(calls):
        raise MalignancyError(f"tissue has {tissue_arr.shape[0]} entries for {len(calls)} cells")

    frame = calls.assign(tissue=tissue_arr)
    if role is not None:
        role_arr = np.asarray([str(r) for r in role], dtype=object)
        if role_arr.shape[0] != len(calls):
            raise MalignancyError(
                f"role has {role_arr.shape[0]} entries for {len(calls)} cells"
            )
        frame = frame.assign(role=role_arr)
        eligible = frame["role"] == "holdout_normal_epi"
        if not eligible.any():
            raise MalignancyError(
                "no held-out normal epithelium. Either assign_cnv_roles fell back "
                "to diploid_only for every patient, or roles were not passed "
                "through — validating on baseline cells would be circular."
            )
    else:
        eligible = frame["tissue"] == "normal"

    normal = frame[eligible & frame["call"].astype(str).isin(["malignant", "non_malignant"])]
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
    task per patient. Take the reference groups from :func:`assign_cnv_roles`:
    ``reference_normal_epi`` as the primary baseline, and each diploid
    compartment as its own additional reference category so inferCNV's
    per-category bounding can suppress cell-type false positives. Never pass the
    held-out cells — they exist precisely so the validation is out-of-sample.

    Cross-check the result with CopyKAT and report the concordance — §4 asks for
    the cross-check, not a winner.
    """
    raise NotImplementedError(
        "W1 — inferCNV needs the real per-patient matrices. "
        "See src/reference/jobs/infercnv.sh."
    )
