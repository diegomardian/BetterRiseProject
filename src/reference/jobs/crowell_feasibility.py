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
    check_depth_is_consistent,
    depth_from_obs,
    read_gene_columns,
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
                     min_umi: int = DETECTION_MIN_UMI,
                     detection: dict[str, np.ndarray] | None = None,
                     counts_per_cell: np.ndarray | None = None,
                     genes_per_cell: np.ndarray | None = None) -> pd.DataFrame:
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
    # PRECOMPUTED PATH. `to_memory()` on a backed CSC materialises every
    # non-zero as float64 -- 549M values, 4.39 GB, and it killed section 110.
    # The job needs five columns of 18,878 and a per-cell depth that obs
    # already carries, so both are read without the matrix and passed in.
    if matrix is None and detection is None:
        # The caller skipped both paths. `main` chose between them on
        # `adata.isbacked`, which is unreliable on a subset view, and this
        # function dereferenced the None instead of saying so.
        raise CrowellError(
            "no source of counts: matrix is None and no precomputed detection "
            "was given. A table built from that would be a table about nothing."
        )
    if detection is not None and matrix is None:
        if counts_per_cell is None or genes_per_cell is None:
            raise CrowellError(
                "precomputed detection needs counts_per_cell and "
                "genes_per_cell too; a depth row of zeros reads as 'these "
                "domains have identical capture', which is a claim."
            )
    # Missing histopathology is not a fourth domain. The first QC-corrected
    # Crowell run stringified missing `typ` values into the literal label
    # "nan" and emitted 20,393 cells under it. That row did not affect the TVA
    # verdict, but it is not a biological group and must not enter the table.
    domain_values = pd.Series(domains, dtype="object")
    has_domain = domain_values.notna().to_numpy()
    domain_labels = domain_values.astype("string")
    # Only computed from the matrix when the caller did not supply them. This
    # pair was unconditional for three attempts on section 110: an edit meant
    # to guard it silently matched nothing, so `main` passed the precomputed
    # arrays in and they were overwritten by `None.sum()` on the next line.
    if counts_per_cell is None:
        counts_per_cell = np.asarray(matrix.sum(axis=1)).ravel()
    if genes_per_cell is None:
        genes_per_cell = np.asarray((matrix > 0).sum(axis=1)).ravel()

    rows = []
    for domain in sorted(domain_labels.loc[has_domain].unique()):
        mask = domain_labels.eq(domain).fillna(False).to_numpy(dtype=bool)
        n_cells = int(mask.sum())
        if n_cells < MIN_CELLS_PER_DOMAIN:
            log.info("  %s: %d cells, below %d — not scored",
                     domain, n_cells, MIN_CELLS_PER_DOMAIN)
            continue
        block = None if matrix is None else matrix[mask]
        floor = (floor_from_obs(obs.loc[mask], n_negative_probes=n_negative_probes)
                 if control_indices is None
                 else negative_floor(block, control_indices, min_umi=min_umi))
        for gene, column in panel_index.items():
            gene_detection = (float(detection[gene][mask].mean())
                              if detection is not None
                              else _detection(block, column, min_umi))
            rows.append({
                "domain": domain,
                "gene": gene,
                "role": GENE_ROLES[gene],
                "n_cells": n_cells,
                "detection": gene_detection,
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
                    float("-inf") if gene_detection <= 0.0 else float(
                        np.log(_mu(gene_detection)
                               / max(_mu(floor["floor_per_probe_mean"]), 1e-12)))),
            })
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["min_log_separation"] = MIN_LOG_SEPARATION
    frame["usable"] = frame["log_separation"] >= MIN_LOG_SEPARATION
    return frame.sort_values(["domain", "log_separation"], ascending=[True, False],
                             ignore_index=True)


def control_referenced_change(table: pd.DataFrame, reference_domain: str,
                              target_domain: str) -> pd.DataFrame:
    """EXPLORATORY. Each gene's separation change, against the control's.

    **Not specified in prereg §5, and not a missing line of pre-registered
    code.** It is a new analysis on a new estimand: Becker's
    ``enrichment_audit`` compares transcript-defined subsets WITHIN an arm,
    this compares ACROSS histological domains. Amendment 2.

    **It carries no interval and it must not be given one.** The band is a
    min/max over the control-role genes present — one, since ACTB is absent
    from this deposit — and a range over one or two genes contains no
    patient-level uncertainty. The biological unit here is one patient.
    """
    wide = table.pivot(index="gene", columns="domain", values="log_separation")
    if reference_domain not in wide.columns or target_domain not in wide.columns:
        return pd.DataFrame()
    roles = table.drop_duplicates("gene").set_index("gene")["role"]
    change = (wide[target_domain] - wide[reference_domain]).rename("change")
    out = change.to_frame().join(roles)
    band = out.loc[out["role"] == "control", "change"]
    out["control_band_low"] = float(band.min()) if len(band) else float("nan")
    out["control_band_high"] = float(band.max()) if len(band) else float("nan")
    out["n_control_genes"] = int(len(band))
    out["opposite_to_controls"] = (
        out["change"] < 0) & (out["control_band_low"] > 0)
    out["reference_domain"] = reference_domain
    out["target_domain"] = target_domain
    out["exploratory"] = True
    out["carries_no_interval"] = True
    return out.reset_index().sort_values("change", ignore_index=True)


def _directional_read(critical: pd.DataFrame, controls: pd.DataFrame,
                      reference_domain: str | None, adenoma_domain: str | None
                      ) -> tuple[dict[str, object], list[str]]:
    """§6's pre-specified direction and the global-capture control, once.

    Extracted because the first version computed them only where the target sat
    BELOW the bar, and in that branch ``no_fall`` cannot occur — a target above
    the bar in the reference and below it in the adenoma has fallen by
    construction. An unreachable branch is a check that cannot fail. §6's
    "no fall is uninterpretable" case lives where BOTH domains clear the bar,
    and it is now evaluated there too.
    """
    fields: dict[str, object] = {}
    detail: list[str] = []
    if not (reference_domain and adenoma_domain):
        return fields, detail
    sep = critical.set_index("domain")["log_separation"]
    if reference_domain not in sep.index or adenoma_domain not in sep.index:
        return fields, detail

    fell = bool(sep[adenoma_domain] < sep[reference_domain])
    fields["prespecified_directional_read"] = "fall_observed" if fell else "no_fall"
    detail.append(
        f"§6's pre-specified direction: {CRITICAL_GENE} "
        f"{sep[reference_domain]:+.3f} in {reference_domain!r} -> "
        f"{sep[adenoma_domain]:+.3f} in {adenoma_domain!r} — "
        + ("a FALL, which §6 fixed in advance as the safe direction because a "
           "field-affected reference understates it."
           if fell else
           "NO fall, which §6 fixed in advance as UNINTERPRETABLE: a reference "
           "already depleted cannot show a further fall, and this may not be "
           "quoted as a negative."))

    by_domain = controls.set_index(["gene", "domain"])["log_separation"]
    worse = []
    for gene in controls["gene"].unique():
        try:
            worse.append(bool(by_domain[(gene, adenoma_domain)]
                              < by_domain[(gene, reference_domain)]))
        except KeyError:
            continue
    if worse:
        n = int(controls["gene"].nunique())
        fields["global_sensitivity_control"] = (
            "worse_in_adenoma" if all(worse) else "not_worse")
        detail.append(
            f"Global capture in {adenoma_domain!r} is "
            f"{'WORSE' if all(worse) else 'not worse'} than in "
            f"{reference_domain!r}, over {n} control-role gene(s). **That is "
            f"global capture only. It does not establish {CRITICAL_GENE}-"
            f"specific sensitivity or any false-negative rate** (prereg §7), "
            f"and with ACTB absent it rests on one gene (Amendment 1 §3).")
    return fields, detail


def verdict(table: pd.DataFrame, *, adenoma_domain: str | None = None,
            reference_domain: str | None = None, patient_n: int | None = None
            ) -> dict:
    """Prereg §5 and Amendment 2: the components, never one word.

    §6's branch table anticipates "separates in REF and falls in TVA" and "fails
    to separate in ANY domain". The outcome on section 231 is neither —
    separates in REF, below the bar in TVA — and the first implementation
    collapsed that to ``NOT ASKABLE IN THE ADENOMA``, which is true of the
    per-cell reading and silent about a direction §6 had pre-specified.

    So the components are returned as fields. **This never returns ASKABLE on
    an outcome where the target sits below the bar in the adenoma**: the target
    did not become independently measurable per cell and no relabelling makes
    it so.
    """
    base: dict[str, object] = {
        "per_cell_feasibility": "not_assessed",
        "prespecified_directional_read": "not_assessed",
        "global_sensitivity_control": "not_assessed",
        "control_referenced_pattern": "not_computed",
        "patient_n": patient_n,
    }
    if table.empty:
        return {**base, "verdict": "NOT ESTIMABLE",
                "detail": "no domain carried enough cells to score."}

    controls = table[table["role"] == "control"]
    n_controls = int(controls["gene"].nunique())
    if n_controls and not controls["usable"].any():
        return {**base,
                "verdict": "READ REFUSED — THE INSTRUMENT IS NOT MEASURING THIS TISSUE",
                "per_cell_feasibility": "failed",
                "global_sensitivity_control": "absent",
                "detail": (
                    "no control-role gene separates from the false-positive "
                    "floor in any domain. Prereg §6's fourth branch — which "
                    "Amendment 1 notes now rests on KRT8 alone, ACTB being "
                    "absent. This is not a statement about the target."
                )}

    critical = table[table["gene"] == CRITICAL_GENE]
    if critical.empty:
        return {**base, "verdict": "NOT ESTIMABLE",
                "detail": f"{CRITICAL_GENE} is absent from the object."}

    usable_domains = sorted(critical.loc[critical["usable"], "domain"])
    if not usable_domains:
        return {**base, "verdict": "NOT ASKABLE HERE",
                "per_cell_feasibility": "failed",
                "detail": (
                    f"{CRITICAL_GENE} clears the bar in no domain (best "
                    f"{critical['log_separation'].max():+.3f} against "
                    f"{MIN_LOG_SEPARATION:+.3f}). Prereg §6: a statement about "
                    f"in-situ sensitivity, NOT about the biology."
                )}

    if adenoma_domain is None:
        return {**base, "verdict": "DOMAIN NAMING NOT SUPPLIED",
                "detail": (
                    f"{CRITICAL_GENE} clears the bar in {usable_domains}. Which "
                    f"is the adenoma is a reading of the deposit's vocabulary "
                    f"and is not guessed; pass --adenoma-domain after --inspect."
                )}

    sep = critical.set_index("domain")["log_separation"]
    controls_by_domain = (controls.set_index(["gene", "domain"])["log_separation"]
                          if n_controls else None)

    if adenoma_domain not in usable_domains:
        fields = dict(base)
        fields["per_cell_feasibility"] = "failed"
        detail = [
            f"{CRITICAL_GENE} clears the bar in {usable_domains} and sits BELOW "
            f"it in {adenoma_domain!r} ({sep.get(adenoma_domain, float('nan')):+.3f} "
            f"against {MIN_LOG_SEPARATION:+.3f}). **Below the limit of "
            f"quantification, NOT censored** — the bar is analyst-chosen "
            f"(Amendment 1 §6), and nothing here builds a background-adjusted "
            f"bound.",
            "**The per-cell spatial reading is NOT licensed.**",
        ]
        extra_fields, extra_detail = _directional_read(
            critical, controls, reference_domain, adenoma_domain)
        fields.update(extra_fields)
        detail += extra_detail
        detail.append(
            f"n = {patient_n if patient_n is not None else 'unrecorded'} "
            f"patient(s). Supportive, NOT confirmatory.")
        return {
            **fields,
            "verdict": ("TARGET BELOW THE ADENOMA USABILITY BAR — "
                        "PER-CELL READING NOT LICENSED; "
                        "DIRECTIONAL READ SUPPORTIVE, NOT CONFIRMATORY"),
            "detail": " ".join(detail),
        }

    fields = dict(base)
    fields["per_cell_feasibility"] = "passed"
    detail = [
        f"{CRITICAL_GENE} clears the bar in {adenoma_domain!r} and in "
        f"{usable_domains}. **This licenses the QUESTION, not an answer.**",
    ]
    extra_fields, extra_detail = _directional_read(
        critical, controls, reference_domain, adenoma_domain)
    fields.update(extra_fields)
    detail += extra_detail
    detail.append(
        f"§7: nothing here licenses a per-cell '{CRITICAL_GENE}-low = silenced' "
        f"claim. n = {patient_n if patient_n is not None else 'unrecorded'} "
        f"patient(s).")
    return {**fields, "verdict": "ASKABLE IN THE ADENOMA",
            "detail": " ".join(detail)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--object", type=Path, required=True,
                        help="one section .h5ad; prereg §2 names 232.h5ad")
    parser.add_argument("--inspect", action="store_true",
                        help="report the vocabulary and map nothing. Run first.")
    parser.add_argument("--domain-column", default=None,
                        help="obs column holding the domain, from --inspect. "
                             "No default, deliberately.")
    parser.add_argument("--adenoma-domain", default=None, nargs="+",
                        help="the value(s) naming the adenoma, from --inspect. "
                             "Several are POOLED AT THE CELL LEVEL into one "
                             "observation — multisection Amendment 1: one "
                             "patient is one observation (invariant 5), and "
                             "pooling weights by cell count. Not guessed.")
    parser.add_argument("--reference-domain", default=None, nargs="+",
                        help="the value(s) naming the reference mucosa, from "
                             "--inspect. Enables §6's pre-specified "
                             "directional read and the exploratory "
                             "control-referenced pattern. Not guessed.")
    parser.add_argument("--patient-n", type=int, default=None,
                        help="biological units behind this run. Invariant 5. "
                             "One section of one block is n=1 and the verdict "
                             "says so.")
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

    # THE MATRIX IS NOT LOADED. `to_memory()` on a backed CSC materialises every
    # non-zero as float64 -- 549M values and 4.39 GB on section 110, which is
    # what killed the first attempt. Five columns of 18,878 are read directly
    # from the sparse arrays, and the per-cell depth comes from obs.
    # THE MATRIX IS NOT LOADED. `to_memory()` on a backed CSC materialises every
    # non-zero as float64 -- 549M values and 4.39 GB on section 110, which is
    # what killed the first attempt. Five columns of 18,878 are read directly
    # from the sparse arrays, and the per-cell depth comes from obs.
    #
    # THE PATH IS NOT CHOSEN BY `adata.isbacked`. That was the first version and
    # it is unreliable on a view: after `adata = adata[pass_qc]` it reported
    # False while `adata.X` also came back None, so BOTH branches were skipped
    # and per_domain_table was handed two Nones. The column reader either works
    # on this file or raises, and that is the condition.
    matrix = None
    detection = counts_per_cell = genes_per_cell = None
    depth_check: dict[str, object] = {}
    if not args.layer:
        try:
            detection = read_gene_columns(args.object, panel_index,
                                          min_umi=DETECTION_MIN_UMI)
        except CrowellError as exc:
            log.info("column reader unavailable (%s); falling back to the "
                     "matrix", exc)
            detection = None

    if detection is not None:
        # read_gene_columns reads the FILE, which is pre-QC; `adata` has
        # already been subset to obs['fil'] == True. The vectors must be put
        # through the same mask or they describe a different population --
        # 694,553 cells against 668,061 on section 110, which is how this was
        # found.
        for gene, vector in detection.items():
            if vector.size != n_cells_before_qc:
                raise CrowellError(
                    f"{gene}'s column is {vector.size} cells and the file holds "
                    f"{n_cells_before_qc}. The column reader and the object "
                    f"disagree about the population; do not align them by "
                    f"truncating."
                )
        detection = {gene: vector[pass_qc] for gene, vector in detection.items()}
        counts_per_cell, genes_per_cell = depth_from_obs(adata.obs)
        if any(v.size != counts_per_cell.size for v in detection.values()):
            raise CrowellError(
                "detection vectors and the depth arrays describe different "
                "numbers of cells after QC."
            )
        depth_check = check_depth_is_consistent(detection, genes_per_cell)
        if not depth_check["consistent"]:
            raise CrowellError(
                f"{depth_check['violations']} cells detect more panel genes "
                f"than obs['nFeature_RNA'] says they detect in total. That "
                f"column describes a different feature space from X, so the "
                f"depth row would not describe this matrix."
            )
        log.info("read %d gene columns directly; depth from obs "
                 "(%d cells checked, %d inconsistent)", len(detection),
                 depth_check["cells_checked"], depth_check["violations"])
    else:
        matrix = (adata.to_memory().layers[args.layer] if args.layer
                  else adata.to_memory().X)
        if matrix is None:
            raise CrowellError(
                "the matrix fallback produced None. Neither the column reader "
                "nor the object yielded counts, and a table built from that "
                "would be a table about nothing."
            )
        require_counts(check_counts_are_integers(matrix))

    # Exactly one path must be live. The first version could skip both silently
    # and hand per_domain_table a pair of Nones.
    if (matrix is None) == (detection is None):
        raise CrowellError(
            f"expected exactly one source of counts; matrix is "
            f"{'set' if matrix is not None else 'None'} and detection is "
            f"{'set' if detection is not None else 'None'}."
        )

    # POOLING, multisection Amendment 1. Several sub-domains of one class in one
    # block are one patient's observation, and §5 already fixed "pool the cells
    # before log_separation, never average after". The parts are scored too and
    # marked exploratory; they contribute nothing to n.
    raw_domains = adata.obs[args.domain_column].to_numpy()
    pooled_domains = raw_domains.astype(object).copy()
    pooled_names: dict[str, str] = {}
    for group in (args.reference_domain or [], args.adenoma_domain or []):
        if len(group) > 1:
            label = "+".join(sorted(group))
            for member in group:
                pooled_domains[raw_domains == member] = label
                pooled_names[member] = label
            log.info("pooling %s -> %r (Amendment 1: one patient, one "
                     "observation)", sorted(group), label)
    adenoma_label = ("+".join(sorted(args.adenoma_domain))
                     if args.adenoma_domain else None)
    reference_label = ("+".join(sorted(args.reference_domain))
                       if args.reference_domain else None)

    _precomputed = dict(detection=detection, counts_per_cell=counts_per_cell,
                        genes_per_cell=genes_per_cell)
    table = per_domain_table(
        matrix, panel_index, pooled_domains,
        control_indices=(None if floor_source == "obs"
                         else np.asarray(controls["negative_indices"], dtype=int)),
        obs=adata.obs if floor_source == "obs" else None,
        n_negative_probes=args.n_negative_probes, **_precomputed)
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

    outcome = verdict(table, adenoma_domain=adenoma_label,
                      reference_domain=reference_label,
                      patient_n=args.patient_n)

    # The parts, separately, so a reader can see whether the lesions agree.
    parts = pd.DataFrame()
    if pooled_names:
        parts = per_domain_table(
            matrix, panel_index, raw_domains,
            control_indices=(None if floor_source == "obs"
                             else np.asarray(controls["negative_indices"], dtype=int)),
            obs=adata.obs if floor_source == "obs" else None,
            n_negative_probes=args.n_negative_probes, **_precomputed)
        parts = parts[parts["domain"].isin(pooled_names)].copy()
        parts["pooled_into"] = parts["domain"].map(pooled_names)
        parts["exploratory"] = True
        parts["contributes_to_n"] = False
        log.info("\n%s\nSUB-DOMAINS SEPARATELY — EXPLORATORY, NOT IN n\n%s",
                 "=" * 72, "=" * 72)
        log.info("%s", parts[["domain", "pooled_into", "gene", "n_cells",
                              "detection", "log_separation"]].to_string(index=False))

    # EXPLORATORY, Amendment 2. Emitted separately, labelled in its own column,
    # and carrying no interval — the band is a min/max over the control-role
    # genes present (one, ACTB being absent) and holds no patient-level
    # uncertainty.
    pattern = pd.DataFrame()
    if reference_label and adenoma_label:
        pattern = control_referenced_change(table, reference_label, adenoma_label)
        if not pattern.empty:
            outcome["control_referenced_pattern"] = "exploratory_two_block" if (
                pattern.loc[pattern.gene == CRITICAL_GENE,
                            "opposite_to_controls"].any()) else "exploratory"
            log.info("\n%s\nCONTROL-REFERENCED CHANGE — EXPLORATORY, POST-HOC, "
                     "NO INTERVAL\n%s", "=" * 72, "=" * 72)
            log.info("%s", pattern[["gene", "role", "change", "control_band_low",
                                    "control_band_high", "n_control_genes",
                                    "opposite_to_controls"]].to_string(index=False))
            log.info("  Amendment 2: a new estimand, not a missing line of "
                     "pre-registered code. It compares ACROSS histological\n"
                     "  domains, where Becker's audit compared transcript-defined "
                     "subsets WITHIN an arm.")
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
        "matrix_loaded": matrix is not None,
        "depth_source": "obs" if matrix is None else "matrix",
        "depth_consistency_check": depth_check or None,
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
        "patient_n_this_run": args.patient_n,
        "reference_domain": args.reference_domain,
        # The RESOLVED labels, which is what the table carries. `*_domain`
        # above is the CLI request and the two are not the same object once
        # sub-domains are pooled.
        "adenoma_domain_label": adenoma_label,
        "reference_domain_label": reference_label,
        "pooled_subdomains": pooled_names or None,
        "pooling_rule": (
            "multisection Amendment 1 — same-class sub-domains are pooled at "
            "the CELL level into one observation per block (invariant 5). "
            "Pooling weights by cell count; the parts are scored separately "
            "and marked exploratory."
        ),
        "amendments": ["Amendment 1 — floor, obs controls, ACTB absent, the "
                       "§2 fallback DEVIATION, below-quantification not "
                       "censored", "Amendment 2 — verdict decomposed into "
                       "fields; control-referenced pattern is exploratory"],
        "width_penalty_vs_n43": 3.01,
        "verdict": outcome,
        "exploratory": False,
        "pre_registered": True,
    }
    if not table.empty:
        tables = [(table, "crowell_feasibility")]
        if not pattern.empty:
            tables.append((pattern, "crowell_control_referenced_exploratory"))
        if not parts.empty:
            tables.append((parts, "crowell_subdomains_exploratory"))
        for frame, name in tables:
            log.info("wrote %s", write_versioned_table(
                frame, name, seed=args.seed, results_dir=args.results_dir,
                allow_dirty=args.allow_dirty, extra_meta=meta,
            ))
    return 0 if outcome["verdict"].startswith("ASKABLE") else 5


if __name__ == "__main__":
    sys.exit(main())
