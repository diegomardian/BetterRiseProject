#!/usr/bin/env python
"""The Kitagawa decomposition on a Lee cohort — SMC or KUL3. W4.

WHY THIS FILE EXISTS
--------------------
The tables in ``results/2026-08-28_033cd51/`` and ``results/2026-08-28_85615b4/``
were produced by a driver that was never committed — their sidecars record
``platform: Windows-11`` and a branch that has since merged. Every number the
project quotes from the Lee arm therefore came from a script nobody can run
again. This is that script, and running it is how the k=2 arm stops being a
reproducibility gap.

WHAT IS DELIBERATELY NOT HERE
-----------------------------
**No estimator.** ``decompose_cohort``, ``bootstrap_over_patients`` and
``attach_intrinsic_ci`` are called unmodified. Nothing here re-derives a term,
a cut point or an interval.

**No labeller.** ``lee_io.load_lee_cohort`` calls W1's ``assign_labels``
(decision #13, one labeller). This driver does not touch maturity.

DEPTH MATCHING IS THE PRIMARY READ — decision #24.1
---------------------------------------------------
Matched **within patient**, **after labelling**, and **on the scored
epithelium** — the same three choices W1's ``build_decomposition_summary``
makes, for the same reasons. Subsampling before labelling would move each
patient's own cut points; matching the pooled population would let a patient
with one deep arm stand in for a patient with two shallow ones.

Lee's arms are ~2.4x apart in median depth after the floor, so the unmatched
read cannot separate a compositional signal from a depth difference. ``--no-
match-depth`` reproduces the older unmatched tables; the table NAME says which
read it came from, rather than the distinction living only in a sidecar.

DEPTH QUANTILE 0.25, NOT 0.10
-----------------------------
``load_lee_cohort`` now defaults to 0.25, matching every GSE178341 job. This
used to fall through to ``assign_labels``' own 0.10, which put the two cohorts'
compositional terms on different scales — see that function's docstring.
Invariant 4 pools these estimates, so the gate has to be the same number.

    python -m src.estimator.run_lee_decomposition --cohort smc --match-depth
    python -m src.estimator.run_lee_decomposition --cohort kul3 --match-depth
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.panel import granularity_rungs, tier_genes
from src.common.provenance import set_global_seeds
from src.estimator.kitagawa import (
    attach_intrinsic_ci,
    bootstrap_over_patients,
    decompose_cohort,
)
from src.estimator.lee_io import (
    LeeCohort,
    build_gene_rung_axis_summary,
    load_lee_cohort,
)
from src.harness.depth_confound import match_arm_depth
from src.reference.labels import NON_EPITHELIAL, TRANSCRIPT_AXES, UNRESOLVED, label_column
from src.schema import coerce_results, write_results

#: Its own seed, distinct from DEFAULT_SEED, so a Lee run and a GSE178341 run
#: are not accidentally correlated through a shared bootstrap stream.
SEED = 20260829

#: Tiers A-D. Tier E is exploratory and is not decomposed.
TARGET_TIERS = ("A", "B", "C", "D")

#: Both transcript axes. `chromatin` and `spatial` are schema-legal but have no
#: data on Lee, so asking for them would emit empty rows rather than a refusal.
AXES = TRANSCRIPT_AXES


def _scorable(labels: pd.DataFrame) -> np.ndarray:
    """Cells that carry a real maturity call — epithelial and depth-resolved.

    `assign_labels` writes UNRESOLVED from `epithelial & ~resolvable` and
    NON_EPITHELIAL from the compartment, and neither mask depends on the axis
    or the rung. So this set is identical across all eight combinations, and
    matching on it keeps the rungs comparable to each other rather than giving
    each rung its own subsample.
    """
    first = labels[label_column(AXES[0], granularity_rungs()[0])].astype(str).to_numpy()
    return (first != UNRESOLVED) & (first != NON_EPITHELIAL)


def _subset(cohort: LeeCohort, keep: pd.Index) -> LeeCohort:
    """The same cohort restricted to `keep`, every frame in step."""
    return LeeCohort(
        study_id=cohort.study_id,
        cells=cohort.cells.loc[keep],
        expression=cohort.expression.loc[keep],
        labels=cohort.labels.loc[keep],
        axis_gene_coverage=cohort.axis_gene_coverage,
        excluded_patients=cohort.excluded_patients,
        n_border_cells=cohort.n_border_cells,
        label_compartment=cohort.label_compartment,
        raw_counts=(
            cohort.raw_counts.loc[keep] if len(cohort.raw_counts) else cohort.raw_counts
        ),
    )


def depth_matched_index(cohort: LeeCohort, *, seed: int) -> tuple[pd.Index, dict[str, int]]:
    """Per-patient depth-matched subsample of the scored epithelium.

    Returns the surviving index and a before/after cell count. A patient
    without both arms among its scorable cells contributes NOTHING rather than
    contributing unmatched — an unmatchable patient is not a matched patient,
    and letting it through is how a depth confound survives a matched read.
    """
    labels, cells = cohort.labels, cohort.cells
    scorable = _scorable(labels)
    keep_positions: list[np.ndarray] = []
    n_before = int(scorable.sum())

    for patient in sorted(cells["patient_id"].unique()):
        in_patient = (cells["patient_id"] == patient).to_numpy() & scorable
        idx = np.flatnonzero(in_patient)
        if idx.size == 0:
            continue
        arm = cells["tissue"].to_numpy()[idx]
        if len(set(arm.tolist())) != 2:
            continue  # one-armed after scoring: cannot be matched, so excluded
        depth = cells["n_counts"].to_numpy(dtype=float)[idx]
        keep_positions.append(idx[match_arm_depth(depth, arm, seed=seed)])

    if not keep_positions:
        raise SystemExit(
            "no patient had two scorable arms — nothing to depth-match. "
            "Check the depth floor and the compartment restriction before rerunning."
        )
    kept = np.concatenate(keep_positions)
    if kept.size == 0:
        # match_arm_depth keeps min(n_normal, n_tumour) per pooled-quantile
        # bin, so an empty result means every bin emptied one arm. There are
        # exactly two ways that happens and they need different responses, so
        # the message names both rather than guessing.
        #
        # Discovered by running this: a fixture with Gaussian depths 4x apart
        # is disjoint in practice, and the disjoint case is NOT a bug — two
        # arms with no common support have nothing to match onto, and no
        # amount of subsampling repairs that comparison.
        raise SystemExit(
            f"depth matching retained 0 of {n_before} scored cells — every depth bin "
            "emptied one arm. Either (a) the patients have fewer scored cells than "
            "the 20 bins, so the cohort is too small to match, or (b) the arms do "
            "not overlap in depth at all, in which case they have no common "
            "distribution to match onto and the contrast is confounded beyond "
            "repair. Check the per-arm depth quantiles before deciding which. "
            "--no-match-depth reports the unmatched read; it does not fix (b)."
        )
    return labels.index[np.sort(kept)], {
        "n_scorable_before_matching": n_before,
        "n_scorable_after_matching": int(kept.size),
    }


def main() -> int:
    set_global_seeds(SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", choices=("smc", "kul3"), default="smc")
    parser.add_argument("--raw-dir", default=None)
    parser.add_argument("--n-boot", type=int, default=2000)
    parser.add_argument("--allow-dirty", action="store_true")
    # A flag rather than a default, so the older unmatched tables stay
    # reproducible and the two reads stay comparable.
    parser.add_argument("--no-match-depth", action="store_true")
    args = parser.parse_args()
    match_depth = not args.no_match_depth

    genes = sorted({g for tier in TARGET_TIERS for g in tier_genes(tier)})
    raw_dir = Path(args.raw_dir) if args.raw_dir else None
    print(f"cohort: {args.cohort} · {len(genes)} genes · depth-matched: {match_depth}")

    cohort = load_lee_cohort(
        args.cohort, target_genes=genes, axes=AXES, raw_dir=raw_dir, depth_quantile=0.25
    )
    print(f"  study_id: {cohort.study_id}")
    print(f"  {len(cohort.cells):,} QC-passed cells · {cohort.cells.patient_id.nunique()} patients")
    if cohort.excluded_patients:
        print(f"  excluded (no matched arms): {', '.join(cohort.excluded_patients)}")
    # MUC2 is absent from both Lee matrices, so `opposite_lineage` scores on 3
    # of its 4 markers there and is NOT the same axis it is on GSE178341.
    # Printed rather than assumed: a meta-analysis that pools it across the two
    # cohorts is pooling two different measurements.
    if cohort.axis_gene_coverage.get("missing"):
        print(f"  MISSING markers: {sorted(cohort.axis_gene_coverage['missing'])}")

    matching: dict[str, int] = {}
    if match_depth:
        before = cohort.cells.groupby("tissue")["n_counts"].median()
        keep, matching = depth_matched_index(cohort, seed=SEED)
        cohort = _subset(cohort, keep)
        after = cohort.cells.groupby("tissue")["n_counts"].median()
        print(
            f"  depth matched: {matching['n_scorable_before_matching']:,} -> "
            f"{matching['n_scorable_after_matching']:,} scored cells"
        )
        print(f"    median depth ratio {_ratio(before):.2f}x -> {_ratio(after):.2f}x")

    summary = build_gene_rung_axis_summary(cohort, genes=genes, axes=AXES)
    if summary.empty:
        raise SystemExit("no patient produced a summary row — nothing to decompose")
    print(f"\nsummary: {len(summary):,} rows · {summary.patient_id.nunique()} patients")

    result = decompose_cohort(summary)
    print(f"decompose_cohort: {len(result):,} rows (one per weighting)")
    print(f"running bootstrap_over_patients, n_boot={args.n_boot} — patients, not cells")
    boot = bootstrap_over_patients(summary, n_boot=args.n_boot, seed=SEED)
    merged = attach_intrinsic_ci(result, boot)

    print("\n" + "=" * 68)
    print("ESTIMABILITY, BY RUNG")
    print("=" * 68)
    print(pd.crosstab(merged.granularity_rung, merged.estimability).to_string())

    suffix = "_depth_matched" if match_depth else ""
    stem = f"decomposition_lee_{args.cohort}{suffix}"
    written = write_results(
        coerce_results(merged),
        stem,
        seed=SEED,
        allow_dirty=args.allow_dirty,
        notes=(
            f"Kitagawa decomposition, {cohort.study_id} ({args.cohort.upper()}), tiers A-D, "
            f"all rungs, axes {list(AXES)}. Estimator is decompose_cohort, called not "
            "reimplemented; labeller is W1's assign_labels (decision #13).\n\n"
            f"Depth-matched: {match_depth} (decision #24.1) — within patient, after "
            "labelling, on the scored epithelium. A patient without two scorable arms "
            "contributes nothing rather than contributing unmatched.\n\n"
            "depth_quantile=0.25, matching every GSE178341 job, so the compositional "
            "terms are on a comparable scale for meta-analysis (invariant 4)."
        ),
    )
    print(f"\nwrote {written}")

    bands = write_versioned_table(
        boot,
        f"{stem}_bands",
        seed=SEED,
        allow_dirty=args.allow_dirty,
        notes=(
            "Cohort-level patient bootstrap bands — all THREE terms, not just the "
            "intrinsic one the schema slot carries. This is the meta-analysis input "
            "(invariant 4): `theta` is the point estimate on the observed patients, "
            "`se` the SD of the bootstrap draws.\n\n"
            "READ n_patients_contributing BEFORE POOLING `intrinsic`: it is a mean "
            "over the patients who RETAINED mature cells, since decompose_cohort "
            "nulls the intrinsic term for every not_estimable patient. n_patients is "
            "the unselected denominator."
        ),
        extra_meta={
            "cohort": args.cohort,
            "study_id": cohort.study_id,
            "n_boot": args.n_boot,
            "depth_matched": match_depth,
            "depth_quantile": 0.25,
            "schema_table": str(written.name),
            "axis_gene_coverage": {k: list(v) for k, v in cohort.axis_gene_coverage.items()},
            **matching,
        },
    )
    print(f"wrote {bands}")
    return 0


def _ratio(medians: pd.Series) -> float:
    """Tumour:normal median depth, or nan if an arm is missing."""
    if "normal" not in medians or "tumour" not in medians or not medians["normal"]:
        return float("nan")
    return float(medians["tumour"] / medians["normal"])


if __name__ == "__main__":
    raise SystemExit(main())
