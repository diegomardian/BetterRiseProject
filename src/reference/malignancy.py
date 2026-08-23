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

from pathlib import Path
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
    role: Any = None,
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

    # WHICH CELLS SET THE THRESHOLD. Pass `role` and it is the patient's own
    # copy-neutral EPITHELIUM; without it, the diploid compartments.
    #
    # This matters more than it looks. Once the CNV baseline is matched normal
    # epithelium, immune and stromal cells are no longer references — they are
    # other cell types scored against an epithelial baseline, so they deviate
    # for ordinary cell-type reasons and score ABOVE the tumour. Thresholding on
    # them puts the cut above almost the entire tumour distribution: on the
    # pilot that called 21 of 2,259 tumour epithelial cells malignant, against
    # ~20% of them sitting above the copy-neutral 90th percentile.
    #
    # `holdout_normal_epi` is the right population — same cell type, same
    # patient, copy-neutral, and never in the baseline. Note that specificity
    # computed on those same cells is then partly circular; validate on
    # `reference_normal_epi` instead, which is disjoint from it.
    if role is not None:
        role_arr = np.asarray([str(r) for r in role], dtype=object)
        if role_arr.shape[0] != scores.shape[0]:
            raise MalignancyError(
                f"role has {role_arr.shape[0]} entries for {scores.shape[0]} cells"
            )
        is_reference = role_arr == "holdout_normal_epi"
    else:
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
        # Cells that DEFINED the threshold are labelled reference, not
        # non_malignant — calling them non-malignant would be circular.
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


#: Fraction of a copy-neutral population expected above its own 90th percentile.
#: The null for :func:`cnv_separation`, and it is 0.10 rather than 0 — a patient
#: sitting at 0.10 has no excess aneuploid population at all.
SEPARATION_NULL: Final[float] = 0.10

#: Enrichment over that null below which no malignancy call is made. 1.5x is a
#: judgement, not a derivation: it says the aneuploid population must be half
#: again what chance produces before a threshold is drawn through it.
MIN_SEPARATION_ENRICHMENT: Final[float] = 1.5


def cnv_separation(
    cnv_score: Any,
    *,
    group: Any,
    patient_id: Any,
    comparator: str = "holdout_normal_epi",
    query: str = "query",
    null: float = SEPARATION_NULL,
    min_enrichment: float = MIN_SEPARATION_ENRICHMENT,
) -> pd.DataFrame:
    """Is there an aneuploid population to call at all? Per patient.

    **The precondition for malignancy calling, and it is not always met.**
    :func:`call_malignancy` will happily threshold any distribution; this asks
    whether there is anything on the other side of the threshold.

    Compares the query against the patient's **held-out normal epithelium** —
    same cell type, same patient, never in the CNV baseline, so the only
    difference left is copy number. Not against the diploid compartments: once
    the baseline is matched normal epithelium, immune and stromal cells deviate
    from it for ordinary cell-type reasons, and comparing to them measures the
    wrong thing.

    The statistic is the fraction of query cells above the comparator's 90th
    percentile, read as enrichment over a null of **0.10** — by construction a
    tenth of copy-neutral cells sit above their own 90th percentile.

    **Expect this to fail for some patients, and expect it to fail
    non-randomly.** MMR-deficient tumours are characteristically near-diploid,
    so there is genuinely no aneuploid population to find. That is a fact about
    the tumour, not a defect in the run, and it is why `separable` is reported
    rather than a threshold being lowered until something appears. See open
    decision #15 — the failures are expected to concentrate in one arm of a
    pre-registered contrast.
    """
    scores = np.asarray(cnv_score, dtype=float)
    groups = np.asarray([str(g) for g in group], dtype=object)
    patients = np.asarray([str(p) for p in patient_id], dtype=object)
    if not (scores.shape[0] == groups.shape[0] == patients.shape[0]):
        raise MalignancyError(
            f"lengths differ: cnv_score {scores.shape[0]}, group "
            f"{groups.shape[0]}, patient_id {patients.shape[0]}"
        )

    rows = []
    for patient in pd.unique(patients):
        here = patients == patient
        q = scores[here & (groups == query)]
        c = scores[here & (groups == comparator)]
        if q.size == 0 or c.size == 0:
            rows.append({
                "patient_id": patient, "n_query": int(q.size),
                "n_comparator": int(c.size), "fraction_above": float("nan"),
                "enrichment": float("nan"), "separable": False,
                "reason": "no query or no held-out comparator",
            })
            continue
        cut = float(np.quantile(c, 1.0 - null))
        fraction = float((q > cut).mean())
        enrichment = fraction / null if null else float("nan")
        separable = bool(enrichment >= min_enrichment)
        rows.append({
            "patient_id": patient, "n_query": int(q.size),
            "n_comparator": int(c.size), "fraction_above": fraction,
            "enrichment": enrichment, "separable": separable,
            "reason": "" if separable else
                      f"enrichment {enrichment:.2f}x below {min_enrichment}x",
        })
    return pd.DataFrame(rows)


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


#: inferCNV's expression cutoff. **0.1 is the documented 10x value**; the
#: package's own guidance is 1 for Smart-seq2 and 0.1 for droplet data, because
#: droplet counts are far sparser per cell. Not a free parameter, and not tuned.
INFERCNV_CUTOFF_10X: Final[float] = 0.1

#: Genes per smoothing window. inferCNV's default. Larger windows suppress noise
#: and blur focal events; 101 is the setting the published analyses use, so it is
#: kept unless there is a reason on this data to move it.
INFERCNV_WINDOW: Final[int] = 101


def write_infercnv_inputs(
    counts: Any,
    gene_names: Any,
    roles: pd.DataFrame,
    *,
    out_dir: Any,
    barcodes: Any = None,
) -> dict[str, Path]:
    """Write the three files inferCNV reads. Returns their paths.

    inferCNV wants a genes x cells count matrix, a two-column annotation file
    mapping each cell to a group, and a gene-position file. The first two are
    written here; the third is a **reference file the caller supplies**, because
    it must match the deposit's genome build and inventing coordinates would be
    worse than failing.

    The annotation groups come straight from :func:`assign_cnv_roles`, with one
    deliberate detail: ``reference_diploid`` is split into one group per
    compartment (``ref_immune``, ``ref_stromal``, ``ref_endothelial``) rather
    than pooled. inferCNV bounds the log fold change by each reference
    category's own mean, so keeping them separate is what suppresses cell-type
    differences being read as copy number — the failure that made this project
    reverse its reference design once already.

    ``holdout_normal_epi`` cells are written as **observations, not
    references**. They exist to be scored out-of-sample by
    :func:`validate_normal_epithelium`; passing them as reference would make
    that check circular.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    required = {"role", "compartment"}
    missing = required - set(roles.columns)
    if missing:
        raise MalignancyError(f"roles is missing column(s): {sorted(missing)}")
    if len(roles) != counts.shape[0]:
        raise MalignancyError(
            f"roles has {len(roles)} rows for {counts.shape[0]} cells"
        )

    names = np.asarray([str(g) for g in gene_names], dtype=object)
    if names.shape[0] != counts.shape[1]:
        raise MalignancyError(
            f"gene_names has {names.shape[0]} entries for {counts.shape[1]} genes"
        )

    # One row per symbol. The deposit carries several Ensembl IDs mapping to the
    # same symbol, and the gene-order file is keyed on symbol — a duplicated
    # row name gives inferCNV two genomic positions for one gene, which makes
    # the smoothing window ambiguous. First occurrence wins, matching how
    # fetch_gene_positions.py resolves the same collision on the other side.
    _, first = np.unique(names, return_index=True)
    gene_keep = np.zeros(names.shape, dtype=bool)
    gene_keep[np.sort(first)] = True
    dropped = int((~gene_keep).sum())
    if dropped:
        print(f"note: {dropped:,} duplicate gene symbols dropped, keeping the "
              f"first occurrence of each ({int(gene_keep.sum()):,} remain)")
    names = names[gene_keep]

    if barcodes is None:
        barcodes = [f"cell{i}" for i in range(counts.shape[0])]
    cells = np.asarray([str(b) for b in barcodes], dtype=object)

    group = np.asarray(
        [
            f"ref_{c}" if r == "reference_diploid" else r
            for r, c in zip(
                roles["role"].astype(str),
                roles["compartment"].astype(str),
                strict=True,
            )
        ],
        dtype=object,
    )
    keep = group != "unusable"
    if not keep.any():
        raise MalignancyError(
            "every cell is 'unusable' — this patient has no viable CNV "
            "reference. Leave malignancy not_called rather than running."
        )

    subset = counts[keep]

    # SPARSE, in 10x layout. The dense route is not merely slower: the largest
    # pilot patient is ~22,000 cells x 43,113 genes, which is 7.8 GB as float64
    # before pandas takes its copy, and the TSV would be ~950 million numbers.
    # Writing that costs hours before inferCNV starts, and a matrix this shape
    # is >90% zeros. CreateInfercnvObject reads a directory of
    # matrix.mtx / barcodes.tsv / genes.tsv directly.
    from scipy import io as sio
    from scipy import sparse

    subset = subset[:, gene_keep]
    matrix = subset.T if sparse.issparse(subset) else sparse.csr_matrix(
        np.asarray(subset).T
    )
    matrix_path = out / "matrix.mtx"
    sio.mmwrite(str(matrix_path), matrix.tocoo())

    # Symbol in BOTH columns on purpose. 10x convention is id then symbol, and
    # inferCNV's reader takes the first — but the gene-order file is keyed on
    # SYMBOL, and a mismatch there does not error, it silently drops every gene
    # from the inference.
    (out / "genes.tsv").write_text(
        "".join(f"{g}\t{g}\n" for g in names)
    )
    (out / "barcodes.tsv").write_text(
        "".join(f"{b}\n" for b in cells[keep])
    )

    annotation_path = out / "annotations.tsv"
    pd.DataFrame({"cell": cells[keep], "group": group[keep]}).to_csv(
        annotation_path, sep="\t", header=False, index=False
    )
    return {"counts": out, "annotations": annotation_path,
            "matrix": matrix_path}


def infercnv_reference_groups(roles: pd.DataFrame) -> list[str]:
    """The group names to pass as ``ref_group_names``. **One group when a matched
    normal exists — not four.**

    ``query`` is what is being called and ``holdout_normal_epi`` is scored
    out-of-sample, so neither may ever appear here.

    **Why a single group, reversing an earlier decision.** The first version
    passed matched normal epithelium *plus* each diploid compartment as its own
    category, reasoning that inferCNV's per-category bounding would suppress
    cell-type differences being read as copy number. The bounding is real, and
    it over-suppresses: at STEP 08 inferCNV runs
    ``subtract_ref_expr_from_obs`` with ``use_bounds=TRUE``, which takes the
    range of the reference-group means per gene and zeroes any observation
    deviation falling *inside* it. Immune, stromal, endothelial and epithelial
    means differ for ordinary cell-type reasons, so that range is wide for most
    genes.

    Measured on the pilot with four groups: **25-30% of all values became
    exactly 1**, the interquartile range collapsed to 0.989-1.009, and tumour
    epithelium scored *below* the diploid reference in four of five patients —
    the reference cells keeping their cell-type deviation from the pooled mean
    while the observations had theirs bounded away.

    The defence against cell-type artifacts is a **cell-type-matched
    reference**, which matched normal epithelium already is. Bounding across
    dissimilar types is a weaker substitute that costs the signal.

    Diploid compartments remain the reference only when a patient has no usable
    matched normal — ``assign_cnv_roles``' ``diploid_only`` strategy — where a
    mismatched reference is the honest cost of having no better one, and those
    calls are flagged rather than pooled.
    """
    roles_seen = set(roles["role"].astype(str))
    if "reference_normal_epi" in roles_seen:
        return ["reference_normal_epi"]

    groups = {
        f"ref_{compartment}"
        for role, compartment in zip(
            roles["role"].astype(str), roles["compartment"].astype(str), strict=True
        )
        if role == "reference_diploid"
    }
    if not groups:
        raise MalignancyError(
            "no reference cells in roles — assign_cnv_roles() reported no "
            "viable strategy, so inferCNV must not be run for this patient."
        )
    return sorted(groups)


def run_infercnv(
    counts: Any,
    gene_names: Any,
    roles: pd.DataFrame,
    *,
    gene_position_file: Any,
    out_dir: Any,
    barcodes: Any = None,
    cutoff: float = INFERCNV_CUTOFF_10X,
    window_length: int = INFERCNV_WINDOW,
    denoise: bool = False,
    hmm: bool = False,
    analysis_mode: str = "subclusters",
    threads: int = 8,
    seed: int = 20260101,
    cleanup: bool = True,
    keep_final: bool = False,
    dry_run: bool = False,
) -> dict[str, Any]:
    """Run inferCNV for ONE patient. W1, weeks 2-3.

    Per patient, never pooled: a baseline built across patients folds germline
    copy-number variation and per-patient capture differences into the
    malignancy call. Drive it with ``src/reference/jobs/infercnv.sh``, which
    submits one array task per patient.

    Parameters worth knowing about rather than accepting:

    `cutoff`
        **0.1, the package's documented value for droplet data** (1 is for
        Smart-seq2). Not tuned here, and it should not be tuned to make the
        calls look better.
    `analysis_mode`
        ``"subclusters"`` rather than ``"samples"``. A tumour is not one clone,
        and a per-sample mode averages subclones together — which pulls
        low-burden malignant cells toward the reference and costs sensitivity
        exactly where this project needs it, since a mature-looking tumour cell
        with few CNVs is the case that decides compositional vs intrinsic.
        Slower. That is the trade being made deliberately.
    `denoise`
        **Off by default, and that is a correction.** Denoising sets every value
        within 1.5 SD of the reference mean *to* the mean. It is a plotting aid
        — it makes heatmaps legible — and it destroys the quantity this pipeline
        actually uses, because the per-cell score is the mean squared deviation
        from 1 and denoising sets most of those deviations to zero.

        The pilot ran with it on and produced scores around 3e-4, implying
        per-gene deviations of ~0.017 where a real copy-number change moves the
        smoothed residual by 0.1-0.3. Four of five patients also came back with
        tumour epithelium scoring *below* the diploid reference, which is
        consistent with denoising being applied to the observations while the
        reference keeps its noise. Turn it on only if someone wants the
        heatmaps, and do not score off a denoised object.
    `hmm`
        **Off by default.** The HMM produces discrete CNV state calls and
        roughly triples the runtime, and this pipeline does not use them:
        :func:`call_malignancy` thresholds a continuous score against each
        patient's own reference cells. Turn it on if someone wants CNV states as
        a result in their own right.
    `gene_position_file`
        Required, no default. inferCNV needs genomic coordinates and **this
        deposit is GRCh37_liftover_v28 (hg19)** — a GRCh38 gene-order file would
        misplace genes and produce chromosome-arm artifacts that look exactly
        like real CNVs. Supply an hg19/GENCODE-v19-era file and add its
        ``data/manifest.csv`` row in the same PR (CONTRIBUTING §4).

    With `dry_run`, writes the inputs and returns the command without executing
    it, so the wiring is testable without R installed.
    """
    if analysis_mode not in {"samples", "subclusters"}:
        raise MalignancyError(
            f"analysis_mode must be 'samples' or 'subclusters', got {analysis_mode!r}"
        )
    positions = Path(gene_position_file)
    if not positions.exists():
        raise MalignancyError(
            f"gene_position_file {positions} does not exist. inferCNV needs "
            f"genomic coordinates, and this deposit is hg19 "
            f"(GRCh37_liftover_v28) — a GRCh38 file would misplace genes and "
            f"produce arm-level artifacts indistinguishable from real CNVs. "
            f"Fetch an hg19 gene-order file and record it in data/manifest.csv."
        )

    out = Path(out_dir)
    paths = write_infercnv_inputs(
        counts, gene_names, roles, out_dir=out, barcodes=barcodes
    )
    references = infercnv_reference_groups(roles)

    script = out / "run_infercnv.R"
    script.write_text(_INFERCNV_R_TEMPLATE.format(
        matrix=paths["matrix"],
        genes=out / "genes.tsv",
        barcodes=out / "barcodes.tsv",
        annotations=paths["annotations"],
        positions=positions,
        ref_groups=", ".join(f'"{g}"' for g in references),
        cutoff=cutoff,
        window=window_length,
        denoise="TRUE" if denoise else "FALSE",
        hmm="TRUE" if hmm else "FALSE",
        analysis_mode=analysis_mode,
        threads=int(threads),
        seed=int(seed),
        out_dir=out,
    ))
    command = ["Rscript", str(script)]

    result: dict[str, Any] = {
        "command": command,
        "reference_groups": references,
        "out_dir": out,
        **paths,
        "script": script,
    }
    if dry_run:
        result["ran"] = False
        return result

    import subprocess

    # STREAMED, not captured. subprocess.run(capture_output=True) buffers until
    # the process exits, which for a run measured in hours means no way to tell
    # a working job from a stuck one until it is over. Each line is echoed as it
    # arrives — so `tail -f` on the SGE job log follows inferCNV live — and kept
    # for the saved log and the error message.
    log_path = out / "infercnv_R.log"
    lines: list[str] = []
    with log_path.open("w") as log:
        log.write(f"$ {' '.join(command)}\n\n")
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            log.flush()
            lines.append(line)
        returncode = process.wait()

    result["ran"] = True
    result["returncode"] = returncode
    result["stdout"] = "".join(lines)
    result["stderr"] = ""

    if returncode != 0:
        tail = "".join(lines[-40:])
        raise MalignancyError(
            f"inferCNV failed (exit {returncode}). Full log: {log_path}\n{tail}"
        )

    # Echo the lines that say whether the run means anything, rather than
    # capturing them and discarding them on success. The gene count AFTER the
    # gene-order intersection is the one that matters: a large drop means the
    # symbol join failed, and inferCNV does not error on that — it infers from
    # whatever survived.
    for line in lines:
        if line.startswith(("matrix:", "after gene-order intersection:")):
            print(f"  {line.rstrip()}")

    if cleanup:
        result["freed_bytes"] = cleanup_infercnv_run(out, keep_final=keep_final)
    return result


#: Files a finished run must keep. Everything else inferCNV writes is a staged
#: intermediate it uses to resume a crashed run, and nothing downstream reads
#: them.
INFERCNV_KEEP: Final[tuple[str, ...]] = (
    "cnv_scores.csv",      # the result
    "annotations.tsv",     # the roles, needed to join reference vs query
    "genes.tsv",
    "barcodes.tsv",
    "run_infercnv.R",      # provenance: the exact settings used
    "infercnv_R.log",      # provenance: what inferCNV actually did
    "infercnv_subclusters.observation_groupings.txt",
)


def cleanup_infercnv_run(out_dir: Any, *, keep_final: bool = False) -> int:
    """Delete a finished run's intermediates. Returns bytes freed.

    **Not tidiness — feasibility.** inferCNV writes every pipeline stage to disk
    so a crashed run can resume, and on this data each stage is the size of the
    expression matrix. One pilot patient produced ~16 GB across five patients,
    and the project filesystem has 55 GB in total. A 62-patient run would fill
    it somewhere around patient fifteen.

    Refuses to delete anything if ``cnv_scores.csv`` is missing, because that is
    the only artifact that cannot be recomputed without re-running the
    inference. A run that failed keeps everything, so it can be diagnosed.

    ``run.final.infercnv_obj`` is deleted by default. It is the only route to
    recomputing the score without re-running, which sounds worth keeping until
    you multiply it by 62 patients — for the largest it is several hundred MB.
    Pass `keep_final` when disk allows and you expect the score definition to
    change.
    """
    out = Path(out_dir)
    if not (out / "cnv_scores.csv").exists():
        print(f"note: {out.name} has no cnv_scores.csv — keeping everything so "
              f"the failure can be diagnosed.")
        return 0

    keep = set(INFERCNV_KEEP)
    if keep_final:
        keep.add("run.final.infercnv_obj")

    freed = 0
    for path in out.iterdir():
        if path.name in keep or not path.is_file():
            continue
        freed += path.stat().st_size
        path.unlink()
    if freed:
        print(f"  cleaned {out.name}: freed {freed / 1e9:.1f} GB")
    return freed


#: The R side. Kept as a template rather than a checked-in .R file so the
#: settings and their reasons live next to each other — the parameters are the
#: judgement calls, and separating them from the docstring is how they drift.
_INFERCNV_R_TEMPLATE = """\
# Generated by src/reference/malignancy.py:run_infercnv — do not edit by hand.
suppressPackageStartupMessages(library(infercnv))
set.seed({seed})

# Load the sparse matrix here rather than handing CreateInfercnvObject a path.
# It does read.table() on whatever path it is given and has no 10x-directory
# reader in this version, so a directory fails with "not a regular file". A
# dense TSV would be ~950 million numbers for the largest patient, so the matrix
# is read from Matrix Market and passed as an object instead.
suppressPackageStartupMessages(library(Matrix))
counts <- Matrix::readMM("{matrix}")
counts <- as(counts, "CsparseMatrix")
rownames(counts) <- read.delim("{genes}", header = FALSE)$V1
colnames(counts) <- read.delim("{barcodes}", header = FALSE)$V1
cat("matrix:", nrow(counts), "genes x", ncol(counts), "cells\n")

obj <- CreateInfercnvObject(
  raw_counts_matrix = counts,
  annotations_file  = "{annotations}",
  gene_order_file   = "{positions}",
  ref_group_names   = c({ref_groups}),
  delim             = "\t"
)
cat("after gene-order intersection:", nrow(obj@expr.data), "genes\n")

result <- infercnv::run(
  obj,
  cutoff              = {cutoff},
  out_dir             = "{out_dir}",
  window_length       = {window},
  cluster_by_groups   = TRUE,
  analysis_mode       = "{analysis_mode}",
  denoise             = {denoise},
  HMM                 = {hmm},
  num_threads         = {threads},
  no_plot             = TRUE
)

# no_plot = TRUE skips the final plot, and in inferCNV that is also the step
# that writes infercnv.observations.txt / infercnv.references.txt. Rather than
# turn plotting back on — a 400 MB text matrix and a heatmap nobody reads, per
# patient — take the score straight off the final object.
#
# Mean squared deviation from 1 across genes, which is the conventional per-cell
# CNV score and exactly what call_malignancy() thresholds. Squared, not
# absolute: a cell with a few large deviations is more plausibly aneuploid than
# one with many tiny ones, and the absolute deviation treats those the same.
expr <- result@expr.data

# Diagnostics, so the score's scale can be checked rather than assumed. A
# copy-neutral cell sits near 1 and a real copy-number change moves the smoothed
# residual by 0.1-0.3. If the spread here is an order of magnitude smaller, or
# most values are exactly 1, the matrix has been flattened — which is what
# denoise = TRUE does — and the score is measuring residual noise.
cat("expr range:", signif(range(expr), 4), "\n")
cat("expr quantiles:", signif(quantile(expr, c(0.01, 0.25, 0.5, 0.75, 0.99)), 4), "\n")
cat("fraction exactly 1:", signif(mean(expr == 1), 4), "\n")

scores <- colMeans((expr - 1)^2)
write.csv(
  data.frame(cell = names(scores), cnv_score = as.numeric(scores)),
  file = file.path("{out_dir}", "cnv_scores.csv"),
  row.names = FALSE
)
cat("wrote cnv_scores.csv for", length(scores), "cells\n")
"""


def read_infercnv_score_table(out_dir: Any) -> pd.DataFrame:
    """Per-cell CNV score and role from a finished inferCNV run.

    Reads ``cnv_scores.csv`` — written by the R step, which computes the mean
    squared deviation from 1 directly off the final object — and joins the role
    each cell was given, so :func:`call_malignancy` has both the scores and the
    reference cells that set the threshold.

    Preferred over :func:`read_infercnv_scores`. That one parses
    ``infercnv.observations.txt``, which inferCNV only writes when plotting is
    on: a ~400 MB text matrix per patient, produced solely to be read back and
    reduced to one number per cell.
    """
    directory = Path(out_dir)
    scores_path = directory / "cnv_scores.csv"
    if not scores_path.exists():
        raise MalignancyError(
            f"no {scores_path}. Either the run did not finish, or it predates "
            f"the R step that writes it — check {directory / 'infercnv_R.log'}."
        )
    scores = pd.read_csv(scores_path)

    annotations_path = directory / "annotations.tsv"
    if annotations_path.exists():
        annotations = pd.read_csv(
            annotations_path, sep="\t", header=None, names=["cell", "group"]
        )
        scores = scores.merge(annotations, on="cell", how="left")
    return scores


def read_infercnv_scores(
    observations_file: Any, *, reference_file: Any = None
) -> pd.Series:
    """Per-cell CNV score from inferCNV's output, shaped for
    :func:`call_malignancy`.

    inferCNV writes a genes x cells matrix of smoothed residual expression,
    centred so a copy-neutral cell sits near 1. The conventional per-cell score
    is the **mean squared deviation from 1** across genes — higher means more
    aneuploid — which is what `call_malignancy` expects.

    Squared, not absolute: a cell with a few large deviations is more plausibly
    aneuploid than one with many tiny ones, and the mean absolute deviation
    treats those the same.

    Pass `reference_file` (``infercnv.references.txt``) to score the reference
    cells on the same scale in the same call. They are needed — the threshold is
    a quantile of the reference cells' own scores — and reading them separately
    is how the two end up on different scalings.
    """
    def _score(path: Any) -> pd.Series:
        frame = pd.read_csv(path, sep=r"\s+", index_col=0)
        return ((frame - 1.0) ** 2).mean(axis=0)

    scores = _score(observations_file)
    if reference_file is not None:
        scores = pd.concat([scores, _score(reference_file)])
    if scores.empty:
        raise MalignancyError(f"no cells in {observations_file}")
    return scores.rename("cnv_score")
