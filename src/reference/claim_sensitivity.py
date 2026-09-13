"""Every fixed choice the adenoma reading could depend on, in one table.

    from src.reference.claim_sensitivity import (
        contrasts_by_denominator, estimability_attrition,
    )

WHY THIS IS THE SUBMISSION'S ROBUSTNESS TABLE. The external audit's second
reviewer objection is *"your positive result depends on post-hoc statistics,
labels and estimability rules"* (``DECISION_2026-09-11_scope_and_pivot.md`` §3).
The response it names is partial and specific: **publish all fixed weightings,
denominators, exclusions and statistic disagreements**. This module assembles
them so the prose can point at one artifact instead of asserting robustness.

WHAT IS VARIED. Each is a design choice the project made, and a reader should be
able to see whether the reading survives each one:

``denominator``  ``resolved`` (decision #14's primary, mature cells over cells
                 resolved to a maturity bin) against ``all_epithelial`` (the
                 denominator #14 rejected). The two decompositions genuinely
                 differ -- not a rounding difference.
``weighting``    ``normal``, ``tumour`` and ``doubly_robust``; the terms are
                 written so each one closes its own three-term identity.
``statistic``    the four scale-free constructions in
                 ``adenoma_decomposition_scales``: ``log_ratio``
                 (pre-registered, load-bearing), ``share_abs``,
                 ``share_signed``, ``ratio``.
``rung``         the four granularities; the curve is the point, not one rung.

WHAT IS **NOT** HERE, AND WHY. Detection-scale statistics (``cloglog``,
``log2_cp10k``, ``detection``) are a different estimand measured on a different
scale, and ``adenoma_specificity_disagreements.parquet`` already reports where
they disagree. Folding them in would compare two estimands and call it a
sensitivity analysis. The two-sample power and interval work are likewise
separate artifacts.

A ``None`` intrinsic term is **not** a zero here. ``scale_free`` drops those rows
and the attrition table counts them, so a gene that could not be estimated never
enters a contrast as a zero (invariant 1).
"""

from __future__ import annotations

import pandas as pd

from src.reference.jobs.adenoma_decomposition_scales import (
    LOAD_BEARING,
    STATISTICS,
    contrasts,
    scale_free,
)


class SensitivityError(ValueError):
    """A decomposition frame that cannot carry the sensitivity table."""


def _require(frame: pd.DataFrame, columns: set[str], what: str) -> None:
    missing = sorted(columns - set(frame.columns))
    if missing:
        raise SensitivityError(f"{what} is missing {missing}")


def contrasts_by_denominator(
    primary: pd.DataFrame, all_epithelial: pd.DataFrame, *, seed: int
) -> pd.DataFrame:
    """Every contrast on every statistic, under both denominators.

    A thin loop over the machinery that already exists: the same
    ``scale_free`` and ``contrasts`` the primary reading uses, run a second time
    on the rejected denominator. Nothing about the statistic computation is
    re-implemented, so a change there cannot silently desynchronise this table.
    """
    frames = []
    for denominator, split in (
        ("resolved", primary), ("all_epithelial", all_epithelial)
    ):
        _require(
            split,
            {"patient_id", "gene", "granularity_rung", "weighting",
             "intrinsic", "compositional"},
            f"{denominator} decomposition",
        )
        values = scale_free(split)
        for statistic in STATISTICS:
            block = contrasts(values, statistic=statistic, seed=seed)
            # The cohort the contrast was drawn from, so a reader can see the
            # denominator's own attrition in the same row as the estimate.
            cohort = (
                split.groupby("granularity_rung")["patient_id"].nunique()
            )
            block["n_patients_cohort"] = block["granularity_rung"].map(cohort)
            block["n_missing"] = (
                block["n_patients_cohort"] - block["n_patients"]
            )
            frames.append(block.assign(denominator=denominator))
    out = pd.concat(frames, ignore_index=True)
    out["load_bearing"] = out["statistic"].eq(LOAD_BEARING)
    ordered = [
        "denominator", "granularity_rung", "weighting", "statistic",
        "contrast", "gene", "other", "cross_block", "summary",
        "n_patients", "n_patients_cohort", "n_missing",
        "centre", "ci_low", "ci_high", "excludes_zero", "load_bearing",
    ]
    return out[ordered].sort_values(
        ["denominator", "granularity_rung", "weighting", "statistic", "contrast"],
        ignore_index=True,
    )


def estimability_attrition(
    primary: pd.DataFrame,
    all_epithelial: pd.DataFrame,
    compositional: pd.DataFrame,
) -> pd.DataFrame:
    """Where patients leave the analysis, per denominator, rung, weighting, gene.

    Two different gates drop rows and they are counted separately:

    - the **intrinsic** gate on ``n_cells_mature`` (``estimability``), which is a
      property of the patient and the rung and so is identical under both
      denominators and all three weightings;
    - the **compositional** gate on ``n_cells_resolved`` (decision #22), carried
      on its own committed table and joined here by (patient, gene, rung).

    ``n_log_ratio_undefined`` counts the rows the load-bearing statistic cannot
    read: ``log(|i| / |c|)`` is undefined where either term is exactly zero, and
    at ``epithelial`` the compositional term is exactly zero for the whole panel
    by construction. That is the pre-registered behaviour, not a defect.
    """
    _require(
        compositional,
        {"patient_id", "gene", "granularity_rung", "compositional_estimability"},
        "compositional estimability",
    )
    rows = []
    for denominator, split in (
        ("resolved", primary), ("all_epithelial", all_epithelial)
    ):
        _require(
            split,
            {"patient_id", "gene", "granularity_rung", "weighting",
             "intrinsic", "compositional", "estimability"},
            f"{denominator} decomposition",
        )
        values = scale_free(split)
        values["log_ratio_defined"] = values["log_ratio"].notna()
        frame = split.merge(
            values[["patient_id", "gene", "granularity_rung", "weighting",
                    "log_ratio_defined"]],
            on=["patient_id", "gene", "granularity_rung", "weighting"],
            how="left",
        ).merge(
            compositional[["patient_id", "gene", "granularity_rung",
                            "compositional_estimability"]],
            on=["patient_id", "gene", "granularity_rung"],
            how="left",
        )
        for (rung, weighting, gene), block in frame.groupby(
            ["granularity_rung", "weighting", "gene"], observed=True
        ):
            # A left merge can turn the indicator object-dtyped where no row
            # matched. Go through pandas' nullable boolean rather than
            # ``fillna(False)`` on object, which downcasts (and is deprecated);
            # a missing match is "not defined", not a crash.
            defined = (
                block["log_ratio_defined"].astype("boolean").fillna(False)
                .astype(bool)
            )
            rows.append({
                "denominator": denominator,
                "granularity_rung": rung,
                "weighting": weighting,
                "gene": gene,
                "n_patient_gene_rows": int(len(block)),
                "n_patients": int(block["patient_id"].nunique()),
                "n_intrinsic_present": int(block["intrinsic"].notna().sum()),
                "n_intrinsic_missing": int(block["intrinsic"].isna().sum()),
                "n_estimability_ok": int((block["estimability"] == "ok").sum()),
                "n_estimability_wide": int(
                    (block["estimability"] == "wide_interval").sum()
                ),
                "n_estimability_not_estimable": int(
                    (block["estimability"] == "not_estimable").sum()
                ),
                "n_compositional_ok": int(
                    (block["compositional_estimability"] == "ok").sum()
                ),
                "n_compositional_wide": int(
                    (block["compositional_estimability"] == "wide_interval").sum()
                ),
                "n_compositional_missing": int(
                    block["compositional_estimability"].isna().sum()
                ),
                "n_log_ratio_defined": int(defined.sum()),
                "n_log_ratio_undefined": int((~defined).sum()),
            })
    return pd.DataFrame(rows).sort_values(
        ["denominator", "granularity_rung", "weighting", "gene"],
        ignore_index=True,
    )


def cross_block_survival(sensitivity: pd.DataFrame) -> pd.DataFrame:
    """For each fixed-choice cell: how many cross-block contrasts exclude zero.

    The claim under audit is that the targets separate from the block of
    controls/identity, so the headline is a cross-block count. A cell here is one
    (denominator, rung, weighting, statistic); ``n_contrasts`` halves the ordered
    pairs, matching how the prose counts unordered ones.
    """
    frame = sensitivity[sensitivity["cross_block"]].copy()
    grouped = (
        frame.groupby(
            ["denominator", "granularity_rung", "weighting", "statistic"],
            observed=True,
        )
        .agg(
            n_contrasts=("excludes_zero", lambda s: int(len(s)) // 2),
            n_excluding_zero=("excludes_zero",
                              lambda s: int(s.sum()) // 2),
            n_patients_min=("n_patients", "min"),
            n_patients_max=("n_patients", "max"),
        )
        .reset_index()
    )
    grouped["share_excluding_zero"] = (
        grouped["n_excluding_zero"] / grouped["n_contrasts"]
    )
    grouped["load_bearing"] = grouped["statistic"].eq(LOAD_BEARING)
    return grouped.sort_values(
        ["denominator", "granularity_rung", "weighting", "statistic"],
        ignore_index=True,
    )
