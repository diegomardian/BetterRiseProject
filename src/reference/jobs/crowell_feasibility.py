"""Is the intrinsic reading askable on WTx spatial data? The feasibility read.

    python -m src.reference.jobs.crowell_feasibility --object <232.h5ad> --inspect
    python -m src.reference.jobs.crowell_feasibility --object <232.h5ad> \\
        --domain-column <from --inspect>

Pre-registered in ``docs/prereg_crowell_feasibility.md``, committed in
``cd7e4a5`` **before anything was downloaded**. Read §6 and §7 before reading any
number this produces.

WHAT MAKES THIS DEPOSIT WORTH THE TIME, in one line: **the label is histology.**
REF / TVA / CRC are pathologist-drawn shapes over H&E, not a function of any
transcript, so a domain assignment cannot be circular with respect to GUCA2A.
`prereg_becker_replication.md`'s RESULT is what the alternative costs — a mature
label built from the mature-colonocyte program made the one arm where GUCA2A
outran its controls unreadable, because the label and the target are the same
program.

THE FLOOR IS PER-PROBE, AND THAT IS NOT A DETAIL. The obvious floor — "fraction
of cells carrying at least one count on any of the 50 negative probes" — is a
union over 50 features and a gene's detection is one feature. Comparing them
overstates the floor by roughly the probe count and would fail every gene. The
comparable quantity is the MEAN PER-PROBE detection, and that is what
``negative_floor`` returns. The union rate is reported too, as a per-cell
contamination indicator, under a name that says which it is.

THE SEPARATION IS ON THE DETECTION SCALE. ``cloglog(p) = log(mu)``, so a
difference is a log ratio of expected UMIs and IS comparable across genes with
different baselines. A ratio of probabilities is not, and this repository has
now reintroduced that error twice (`docs/HANDOFF.md` §6d, and again in
``becker_feasibility``'s first ``fold_vs_chen``).

WHAT THIS CANNOT DECIDE, restated because the specificity control is good enough
to be quoted as something it is not: the negative probes bound the per-cell
false-POSITIVE floor. They say nothing about the false-negative rate, so §6g
stands and **no per-cell "GUCA2A-low = silenced" claim is licensed by anything
here**. Prereg §7.
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.reference.crowell_io import (
    OBS_NEGATIVE_FEATURES,
    OBS_QC_PASS,
    SECTIONS,
    CrowellError,
    check_counts_are_integers,
    control_features,
    domain_vocabulary,
    floor_from_obs,
    open_section,
    qc_pass_mask,
    require_controls,
    require_counts,
)
from src.reference.jobs.coexpression_silencing import DETECTION_MIN_UMI, GENE_ROLES

log = logging.getLogger(__name__)

#: A gene is usable in a domain when its expected-UMI rate is at least this
#: many times the per-probe false-positive floor, on the detection scale.
#: log(3): below 3x the floor, more than a third of the detected signal could
#: be the false-positive process. Fixed here before any Crowell byte was read.
MIN_LOG_SEPARATION = float(np.log(3.0))

#: The gene whose failure ends the spatial reading, as in B1.
CRITICAL_GENE = "GUCA2A"

#: Below this many cells a domain is not scored. A detection rate over a
#: handful of cells is not a detection rate.
MIN_CELLS_PER_DOMAIN = 200

def _mu(p):
    return -np.log1p(-np.clip(np.asarray(p, dtype=float), 0.0, 1 - 1e-12))


def find_panel(var_names) -> tuple[dict[str, int], list[str]]:
    """Panel symbol -> column. Whole-transcriptome, so no identifier gymnastics
    is expected — but absence is still reported rather than assumed away."""
    names = np.asarray([str(v) for v in var_names])
    index = {}
    for gene in GENE_ROLES:
        hit = np.flatnonzero(names == gene)
        if hit.size:
            index[gene] = int(hit[0])
    return index, sorted(set(GENE_ROLES) - set(index))


def _detection(matrix, column: int, min_umi: int) -> float:
    values = matrix[:, column]
    values = (np.asarray(values.todense()).ravel()
              if hasattr(values, "todense") else np.asarray(values))
    return float((np.asarray(values, dtype=float) >= min_umi).mean())


def negative_floor(matrix, control_indices: np.ndarray, *,
                   min_umi: int = DETECTION_MIN_UMI) -> dict[str, float]:
    """The per-cell false-positive floor, two ways, named so they cannot be swapped.

    ``per_probe_mean`` is the comparable one: the mean over control probes of
    each probe's own detection rate. One probe against one gene.

    ``any_probe`` is a union over all of them. It is a useful per-cell
    contamination indicator and it is NOT the floor a single gene should be
    compared against — it is larger by roughly the probe count and would fail
    everything.
    """
    if control_indices.size == 0:
        raise CrowellError("no control probes; the floor would be zero")
    rates = np.array([_detection(matrix, int(c), min_umi)
                      for c in control_indices], dtype=float)
    block = matrix[:, np.sort(control_indices)]
    block = (np.asarray(block.todense()) if hasattr(block, "todense")
             else np.asarray(block))
    any_hit = float((np.asarray(block, dtype=float) >= min_umi).any(axis=1).mean())
    return {
        "floor_per_probe_mean": float(rates.mean()),
        "floor_per_probe_max": float(rates.max()),
        "n_control_probes": int(rates.size),
        "any_probe_union_rate": any_hit,
    }


def per_domain_table(matrix, panel_index: dict[str, int], domains,
                     control_indices: np.ndarray | None = None, *,
                     obs: pd.DataFrame | None = None,
                     n_negative_probes: int | None = None,
                     min_umi: int = DETECTION_MIN_UMI) -> pd.DataFrame:
    """One row per (domain, gene): detection, the floor, separation, depth.

    The floor comes from ``control_indices`` when the control probes are
    features of the matrix, and from ``obs`` when they were summarised per cell
    before the object was written — which is what the Crowell deposit does. One
    of the two is required; there is no path to a floor of zero.
    """
    if control_indices is None and obs is None:
        raise CrowellError(
            "give either control_indices (probes in var) or obs (probes "
            "summarised per cell). A floor of zero is a gate every gene clears."
        )
    # Missing histopathology is not a fourth domain. The first QC-corrected
    # Crowell run stringified missing `typ` values into the literal label
    # "nan" and emitted 20,393 cells under it. That row did not affect the TVA
    # verdict, but it is not a biological group and must not enter the table.
    domain_values = pd.Series(domains, dtype="object")
    has_domain = domain_values.notna().to_numpy()
    domain_labels = domain_values.astype("string")
    counts_per_cell = np.asarray(matrix.sum(axis=1)).ravel()
    genes_per_cell = np.asarray((matrix > 0).sum(axis=1)).ravel()

    rows = []
    for domain in sorted(domain_labels.loc[has_domain].unique()):
        mask = domain_labels.eq(domain).fillna(False).to_numpy(dtype=bool)
        n_cells = int(mask.sum())
        if n_cells < MIN_CELLS_PER_DOMAIN:
            log.info("  %s: %d cells, below %d — not scored",
                     domain, n_cells, MIN_CELLS_PER_DOMAIN)
            continue
        block = matrix[mask]
        floor = (floor_from_obs(obs.loc[mask], n_negative_probes=n_negative_probes)
                 if control_indices is None
                 else negative_floor(block, control_indices, min_umi=min_umi))
        for gene, column in panel_index.items():
            detection = _detection(block, column, min_umi)
            rows.append({
                "domain": domain,
                "gene": gene,
                "role": GENE_ROLES[gene],
                "n_cells": n_cells,
                "detection": detection,
                # THE DEPTH ROW. Becker's mature label carried 2.04x the arm's
                # median UMIs and lifted every gene; the pass was library size.
                # A cross-domain comparison inherits any depth difference, and
                # this is where a reader sees it before reading anything else.
                "median_counts_per_cell": float(np.median(counts_per_cell[mask])),
                "median_genes_per_cell": float(np.median(genes_per_cell[mask])),
                **floor,
                # cloglog(p) = log(mu): a log ratio of expected UMIs, comparable
                # across genes with different baselines. A ratio of
                # probabilities is not, and this repo has made that error twice.
                # A gene detected in NO cell is not "separated by zero" — it
                # is not separated at all, and -inf says so. Computed
                # explicitly rather than let through np.log, which would emit a
                # divide-by-zero warning and invite someone to clamp it into a
                # finite number that reads like a measurement.
                "log_separation": (
                    float("-inf") if detection <= 0.0 else float(
                        np.log(_mu(detection)
                               / max(_mu(floor["floor_per_probe_mean"]), 1e-12)))),
            })
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["min_log_separation"] = MIN_LOG_SEPARATION
    frame["usable"] = frame["log_separation"] >= MIN_LOG_SEPARATION
    return frame.sort_values(["domain", "log_separation"], ascending=[True, False],
                             ignore_index=True)


def verdict(table: pd.DataFrame, *, adenoma_domain: str | None = None) -> dict:
    """Prereg §5: is the question askable here, and where?

    ``adenoma_domain`` is passed by a human who has read ``--inspect``. It is
    not guessed: the deposit's words are a measurement, and this project has
    put a cohort's usable pairs at zero once by reading a label instead.
    """
    if table.empty:
        return {"verdict": "NOT ESTIMABLE",
                "detail": "no domain carried enough cells to score."}

    controls = table[table["role"] == "control"]
    if len(controls) and not controls["usable"].any():
        return {
            "verdict": "READ REFUSED — THE INSTRUMENT IS NOT MEASURING THIS TISSUE",
            "detail": (
                "no control gene separates from the false-positive floor in any "
                "domain. Prereg §6's fourth branch: this is not a statement "
                "about the target."
            ),
        }

    critical = table[table["gene"] == CRITICAL_GENE]
    if critical.empty:
        return {"verdict": "NOT ESTIMABLE",
                "detail": f"{CRITICAL_GENE} is absent from the object."}

    usable_domains = sorted(critical.loc[critical["usable"], "domain"])
    if not usable_domains:
        return {
            "verdict": "NOT ASKABLE HERE",
            "detail": (
                f"{CRITICAL_GENE} does not separate from the false-positive "
                f"floor in any domain (best "
                f"{critical['log_separation'].max():+.3f} against a bar of "
                f"{MIN_LOG_SEPARATION:+.3f}). Prereg §6: this is a statement "
                f"about in-situ sensitivity, NOT about the biology — the same "
                f"distinction §3 of the Becker prereg drew, where it held."
            ),
        }

    if adenoma_domain is None:
        return {
            "verdict": "ASKABLE — DOMAIN NAMING NOT SUPPLIED",
            "detail": (
                f"{CRITICAL_GENE} separates in {usable_domains}. Which of these "
                f"is the adenoma is a reading of the deposit's own vocabulary "
                f"and is not guessed here; pass --adenoma-domain after "
                f"--inspect."
            ),
        }
    if adenoma_domain not in usable_domains:
        return {
            "verdict": "NOT ASKABLE IN THE ADENOMA",
            "detail": (
                f"{CRITICAL_GENE} separates in {usable_domains} but not in "
                f"{adenoma_domain!r}, which is the domain the estimand is "
                f"defined on. The spatial reading of avenue A cannot proceed "
                f"on this section."
            ),
        }
    return {
        "verdict": "ASKABLE IN THE ADENOMA",
        "detail": (
            f"{CRITICAL_GENE} separates from the false-positive floor in "
            f"{adenoma_domain!r} and in {usable_domains}. **This licenses the "
            f"QUESTION, not an answer.** Prereg §6: a fall from reference to "
            f"adenoma would be understated by a field-affected reference and is "
            f"therefore safe; NO fall is uninterpretable and may not be quoted "
            f"as a negative. Prereg §7: nothing here licenses a per-cell "
            f"'GUCA2A-low = silenced' claim."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--object", type=Path, required=True,
                        help="one section .h5ad; prereg §2 names 232.h5ad")
    parser.add_argument("--inspect", action="store_true",
                        help="report the vocabulary and map nothing. Run first.")
    parser.add_argument("--domain-column", default=None,
                        help="obs column holding the domain, from --inspect. "
                             "No default, deliberately.")
    parser.add_argument("--adenoma-domain", default=None,
                        help="the value in that column naming the adenoma, "
                             "from --inspect. Not guessed.")
    parser.add_argument("--layer", default=None, help="counts layer, if not X")
    parser.add_argument("--n-negative-probes", type=int, default=None,
                        help="negative-probe count for the per-probe floor. "
                             "Inferred from max(nFeature_negprobes) if absent, "
                             "which is CONSERVATIVE: too few probes divides by "
                             "too little and raises the floor. The paper "
                             "reports 50.")
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    adata = open_section(args.object)
    n_cells_before_qc = int(adata.n_obs)
    pass_qc = qc_pass_mask(adata.obs)
    n_cells_after_qc = int(pass_qc.sum())
    var_names = [str(v) for v in adata.var_names]
    panel_index, absent = find_panel(var_names)
    controls = control_features(var_names)

    if args.inspect:
        log.info("%s\nINSPECTION — assume nothing, map nothing\n%s",
                 "=" * 72, "=" * 72)
        log.info("  file                      %s", args.object)
        log.info("  shape                     %d cells x %d features",
                 adata.n_obs, adata.n_vars)
        log.info("  QC pass (%s=True)        %d of %d cells",
                 OBS_QC_PASS, n_cells_after_qc, n_cells_before_qc)
        log.info("  layers                    %s", list(adata.layers.keys()))
        log.info("  panel genes located       %d of %d%s",
                 len(panel_index), len(GENE_ROLES),
                 f" (absent: {absent})" if absent else "")
        log.info("  negative probes           %d via %s",
                 controls["n_negative"], controls["negative_patterns_matched"])
        log.info("  false codes               %d via %s",
                 controls["n_false_code"], controls["false_code_patterns_matched"])
        vocab = domain_vocabulary(adata.obs)
        log.info("  obs columns               %s", vocab["all_obs_columns"])
        for column, values in vocab["candidate_columns"].items():
            log.info("  vocabulary::%-14s %s", column, values)
        if vocab["expected_words_seen"]:
            log.info("  prereg words seen in      %s", vocab["expected_words_seen"])
        else:
            log.warning(
                "  NONE of %s appears in any candidate column. The domain may "
                "live in\n  a shape file rather than obs — Zenodo "
                "10.5281/zenodo.15552301 carries the\n  Napari inputs. Do NOT "
                "map a column to REF/TVA/CRC on resemblance.",
                list(("ref", "tva", "crc")))
        log.info(
            "\n  NEXT: choose --domain-column and --adenoma-domain from the "
            "vocabulary above,\n  BY A HUMAN WHO HAS SEEN IT. Reading a label "
            "rather than the grouping once put\n  Chen_2021's usable pairs at "
            "zero when the true number was 44.")
        return 0

    # The deposit defines `fil` as whether a cell passed quality control. The
    # first feasibility run accidentally used every row (298,151 rather than
    # the 278,691 retained cells) and is invalidated by this line. Filtering is
    # mandatory rather than an option because the alternative population is
    # not the one Crowell released for analysis.
    adata = adata[pass_qc]
    log.info("QC: retained %d of %d cells with obs[%r] == True",
             n_cells_after_qc, n_cells_before_qc, OBS_QC_PASS)

    if not panel_index:
        raise CrowellError(
            "no panel gene matched var_names. This is a whole-transcriptome "
            "deposit, so suspect the identifier space before concluding the "
            "genes are absent — that error has been made four times here."
        )
    if not args.domain_column:
        raise SystemExit(
            "--domain-column is required and has no default. Run --inspect and "
            "read it off the vocabulary. Guessing it is how a label gets read "
            "instead of a grouping."
        )
    if args.domain_column not in adata.obs.columns:
        raise SystemExit(
            f"{args.domain_column!r} is not an obs column; got "
            f"{list(adata.obs.columns)}")
    n_cells_missing_domain = int(adata.obs[args.domain_column].isna().sum())
    if n_cells_missing_domain:
        log.info("domain: excluding %d QC-passing cells with missing obs[%r]",
                 n_cells_missing_domain, args.domain_column)
    # WHERE THE FLOOR COMES FROM. Section 232 carries no control probes in var
    # -- they were summarised into obs before the object was written. Both
    # paths are real; neither may be skipped, because a floor of zero is a gate
    # every gene clears.
    floor_source = "var"
    if controls["n_negative"] == 0:
        if OBS_NEGATIVE_FEATURES not in adata.obs.columns:
            require_controls(controls)          # raises, with the reason
        floor_source = "obs"
        log.info("control probes are not features of X; taking the floor from "
                 "obs['%s'] instead", OBS_NEGATIVE_FEATURES)
    else:
        require_controls(controls)

    matrix = adata.layers[args.layer] if args.layer else adata.X
    if adata.isbacked:
        matrix = adata.to_memory().layers[args.layer] if args.layer \
            else adata.to_memory().X
    require_counts(check_counts_are_integers(matrix))

    table = per_domain_table(
        matrix, panel_index, adata.obs[args.domain_column].to_numpy(),
        control_indices=(None if floor_source == "obs"
                         else np.asarray(controls["negative_indices"], dtype=int)),
        obs=adata.obs if floor_source == "obs" else None,
        n_negative_probes=args.n_negative_probes)
    log.info("\n%s\nPER-DOMAIN DETECTION AGAINST THE FALSE-POSITIVE FLOOR\n%s",
             "=" * 72, "=" * 72)
    if not table.empty:
        log.info("%s", table[["domain", "gene", "role", "n_cells", "detection",
                              "floor_per_probe_mean", "log_separation",
                              "usable"]].to_string(index=False))
        log.info("\n  DEPTH PER DOMAIN — the confound that already fooled this "
                 "project once")
        log.info("%s", table.drop_duplicates("domain")[
            ["domain", "n_cells", "median_counts_per_cell",
             "median_genes_per_cell", "any_probe_union_rate"]].to_string(index=False))

    outcome = verdict(table, adenoma_domain=args.adenoma_domain)
    log.info("\n%s\nVERDICT\n%s", "=" * 72, "=" * 72)
    log.info("  %s\n  %s", outcome["verdict"], outcome["detail"])
    log.info(
        "\n  READ PREREG §6 AND §7 BEFORE QUOTING THIS. §6: every reference "
        "region is\n  adjacent normal from a cancer patient and there is no "
        "healthy arm, so a fall\n  is SAFE and no fall is UNINTERPRETABLE. §7: "
        "the negative probes are a\n  SPECIFICITY control — they bound false "
        "positives, not false negatives, so\n  §6g stands and no per-cell "
        "silencing claim is licensed.")

    meta = {
        "prereg": "docs/prereg_crowell_feasibility.md",
        "prereg_committed_in": "cd7e4a5, before anything was downloaded",
        "section": str(args.object.name),
        "sections_named_in_prereg": SECTIONS,
        "domain_column": args.domain_column,
        "adenoma_domain": args.adenoma_domain,
        "label_is_histological": (
            "REF/TVA/CRC are pathologist-drawn shapes over H&E, not a function "
            "of any transcript, so the label cannot be circular with respect to "
            "GUCA2A. That is why this deposit is worth the time — see the "
            "Becker prereg RESULT for what a transcript label costs."
        ),
        "floor_is_per_probe": (
            "mean over control probes of each probe's own detection rate. The "
            "union over all probes is reported separately as a contamination "
            "indicator and is NOT the comparable floor."
        ),
        "min_log_separation": MIN_LOG_SEPARATION,
        "floor_source": floor_source,
        "n_negative_probes": args.n_negative_probes,
        "qc_filter_column": OBS_QC_PASS,
        "n_cells_before_qc": n_cells_before_qc,
        "n_cells_after_qc": n_cells_after_qc,
        "n_cells_qc_failed_excluded": n_cells_before_qc - n_cells_after_qc,
        "n_cells_missing_domain_excluded": n_cells_missing_domain,
        "negative_probes_are_specificity_only": True,
        "does_not_license": (
            "any per-cell 'GUCA2A-low = silenced' claim. The false-negative "
            "rate is unmeasured here as everywhere else (§6g)."
        ),
        "n_patients_in_deposit": 7,
        "width_penalty_vs_n43": 3.01,
        "verdict": outcome,
        "exploratory": False,
        "pre_registered": True,
    }
    if not table.empty:
        log.info("\nwrote %s", write_versioned_table(
            table, "crowell_feasibility", seed=args.seed,
            results_dir=args.results_dir, allow_dirty=args.allow_dirty,
            extra_meta=meta,
        ))
    return 0 if outcome["verdict"].startswith("ASKABLE") else 5


if __name__ == "__main__":
    sys.exit(main())
