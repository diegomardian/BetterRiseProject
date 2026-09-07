"""§4's statistic, across blocks. The DiD that was being computed in a terminal.

    python -m src.reference.jobs.crowell_multisection

Pre-registered in ``docs/prereg_crowell_multisection.md``. **A thin local read
over committed per-section tables** — it re-reads every
``crowell_feasibility.parquet`` under ``results/``, forms §4's difference in
differences, and aggregates over blocks. No h5ad, no cluster.

WHY IT EXISTS. The per-section job emits ``log_separation`` per (domain, gene).
§4's statistic is a step further:

    Delta(gene) =  d(gene, adenoma) - d(gene, reference)
    DiD(gene)   =  Delta(gene) - Delta(CONTROL_GENE)

and after two sections that was being assembled by hand at a prompt. **A
pre-registered statistic computed ad hoc is a post-hoc statistic**, so it lives
here, reads the sidecars for which domain was which, and cannot quietly change
between sections.

THE FLOOR CANCELS, WHICH IS THE POINT. ``log_separation`` is
``log(mu_gene) - log(mu_floor)``, so the DiD is
``dlog(mu_gene) - dlog(mu_control)`` and the false-positive floor drops out
entirely. That matters here: section 110's floor moves 0.0102 -> 0.0173 between
domains, which would otherwise sit inside every number.

BLOCKS, NOT FILES. ``231`` and ``232`` are FOV ranges of one tissue block and
contribute **one** observation (invariant 5, prereg §5). ``SECTION_TO_BLOCK``
carries the mapping, read from the deposit's own ``metadata.txt``.

WHAT IT REFUSES. Below ``MIN_STUDIES`` blocks it reports the per-block values
side by side and **no interval** — prereg §6, and ``meta.MIN_STUDIES`` is the
same 3 this project has repeatedly declined to go under. At 3 or 4 the interval
is reported and pre-committed to be uninformative; the width table in §6 is
emitted beside it so the reader sees what they are looking at.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.common.io import write_versioned_table
from src.common.paths import RESULTS_DIR
from src.harness.meta import MIN_STUDIES
from src.reference.crowell_io import SECTION_TO_BLOCK, CrowellError
from src.reference.jobs.coexpression_silencing import GENE_ROLES
from src.reference.jobs.crowell_feasibility import MIN_LOG_SEPARATION

log = logging.getLogger(__name__)


def _mu(p):
    return -np.log1p(-np.clip(np.asarray(p, dtype=float), 0.0, 1 - 1e-12))

#: The control the DiD is taken against. ACTB is absent from this deposit, so
#: this is the only `control`-role gene present — feasibility Amendment 1 §3.
CONTROL_GENE = "KRT8"

#: The target whose fall is the prediction, and the identity marker that
#: discriminates "a tier moved" from "this gene moved" (prereg §7).
TARGET_GENE = "GUCA2A"
DISCRIMINATOR_GENE = "CDX2"


def resolve_label(recorded, present: set[str]) -> str | None:
    """The domain label in the TABLE, from what the sidecar recorded.

    The sidecar stores the CLI request and the table stores the pooled label,
    and they are not the same object: section 231 recorded the string
    ``'231_TVA'`` while section 110 recorded
    ``['110_TVA1', '110_TVA2', '110_TVA3']`` for a table whose domain is
    ``'110_TVA1+110_TVA2+110_TVA3'``. Both forms resolve here, and a resolution
    that is not actually in the table **raises** rather than being guessed at —
    a near-miss would silently score the wrong domain.
    """
    if recorded is None:
        return None
    if isinstance(recorded, str):
        candidates = [recorded]
    else:
        members = sorted(str(v) for v in recorded)
        candidates = ["+".join(members)] + members if len(members) > 1 else members
    for candidate in candidates:
        if candidate in present:
            return candidate
    raise CrowellError(
        f"sidecar records {recorded!r}; the table carries {sorted(present)}. "
        f"Tried {candidates}. Do not guess which domain was meant."
    )


def _sections(results_dir: Path) -> list[tuple[Path, dict]]:
    """Every committed per-section feasibility table, with its sidecar."""
    out = []
    for table in sorted(results_dir.glob("*/crowell_feasibility.parquet")):
        sidecar = table.with_suffix("").with_suffix(".meta.json")
        if not sidecar.exists():
            sidecar = table.parent / "crowell_feasibility.meta.json"
        if not sidecar.exists():
            log.warning("%s has no sidecar; skipped — the sidecar is where the "
                        "reference and adenoma domains are recorded", table)
            continue
        out.append((table, json.loads(sidecar.read_text())))
    return out


def per_block_did(results_dir: Path) -> pd.DataFrame:
    """§4's DiD, one row per (block, gene). Supersedes runs of the same block."""
    rows = []
    for table, meta in _sections(results_dir):
        section = str(meta.get("section", "")).replace(".h5ad", "")
        block = SECTION_TO_BLOCK.get(section)
        adenoma, reference = meta.get("adenoma_domain"), meta.get("reference_domain")
        if not (adenoma and reference):
            log.info("  %s: no reference/adenoma pair recorded — skipped",
                     table.parent.name)
            continue
        if block is None:
            raise CrowellError(
                f"section {section!r} is not in SECTION_TO_BLOCK. Blocks are "
                f"the unit of inference and an unmapped section cannot be "
                f"counted — add it from the deposit's metadata.txt."
            )
        frame = pd.read_parquet(table)
        present = set(frame["domain"].astype(str))
        # Amendment 2's diagnosis, RECORDED rather than re-derived. Every value
        # is already implied by `detection` and `median_counts_per_cell` in the
        # per-section table, but a number that has to be recomputed by hand to
        # be checked is a number that gets recomputed wrong: this document's
        # first version of it used Becker's KRT8 in place of Crowell's.
        detection = frame.set_index(["gene", "domain"])["detection"]
        depth = frame.drop_duplicates("domain").set_index("domain")[
            "median_counts_per_cell"]
        adenoma = resolve_label(adenoma, present)
        reference = resolve_label(reference, present)
        if not (adenoma and reference):
            continue
        sep = frame.set_index(["gene", "domain"])["log_separation"]
        if (CONTROL_GENE, adenoma) not in sep or (CONTROL_GENE, reference) not in sep:
            log.warning("  %s: %s missing in one domain — skipped",
                        table.parent.name, CONTROL_GENE)
            continue
        # §5 RULE 3, WHICH THIS CODE DID NOT IMPLEMENT UNTIL NOW. The rule is
        # "a KRT8 separation clearing MIN_LOG_SEPARATION in BOTH" — the
        # sensitivity gate, and it is about the control, never the target. The
        # first version checked only that KRT8 EXISTS in both domains. It never
        # bit (KRT8's worst separation across six blocks is 2.804 against a bar
        # of 1.099), but a pre-registered rule the code does not enforce is
        # this repository's own defect class, and an inclusion rule that cannot
        # exclude is a check that cannot fail.
        control_seps = (sep[(CONTROL_GENE, reference)], sep[(CONTROL_GENE, adenoma)])
        if min(control_seps) < MIN_LOG_SEPARATION:
            log.warning(
                "  %s: EXCLUDED by §5 rule 3 — %s separates at %+.3f / %+.3f "
                "against a bar of %+.3f. The sensitivity gate is about the "
                "control, and this block fails it.",
                table.parent.name, CONTROL_GENE, control_seps[0],
                control_seps[1], MIN_LOG_SEPARATION)
            continue
        control_delta = sep[(CONTROL_GENE, adenoma)] - sep[(CONTROL_GENE, reference)]
        for gene in frame["gene"].unique():
            if (gene, adenoma) not in sep or (gene, reference) not in sep:
                continue
            delta = sep[(gene, adenoma)] - sep[(gene, reference)]
            dlog_mu = float(np.log(_mu(detection[(gene, adenoma)])
                                   / _mu(detection[(gene, reference)])))
            rows.append({
                "block": block, "section": section, "run": table.parent.name,
                "gene": gene, "role": GENE_ROLES.get(gene),
                "delta": float(delta),
                "control_delta": float(control_delta),
                # Amendment 2: the epithelial group rises at or above the depth
                # term and the targets do not. Compare dlog_mu with log_depth.
                "detection_reference": float(detection[(gene, reference)]),
                "detection_adenoma": float(detection[(gene, adenoma)]),
                "dlog_mu": dlog_mu,
                "median_depth_reference": float(depth[reference]),
                "median_depth_adenoma": float(depth[adenoma]),
                "log_depth_ratio": float(np.log(depth[adenoma] / depth[reference])),
                "rises_faster_than_depth": bool(
                    dlog_mu > np.log(depth[adenoma] / depth[reference])),
                # Amendment 4: below the negative-probe floor the observed
                # signal is smaller than what noise alone supplies, so the true
                # rate is consistent with ZERO and this DiD is a bound, not a
                # point. Flagged per block; the consequence is in `aggregate`.
                "below_floor_adenoma": bool(sep[(gene, adenoma)] < 0),
                "below_floor_reference": bool(sep[(gene, reference)] < 0),
                "target_below_floor": bool(
                    sep[(gene, adenoma)] < 0 or sep[(gene, reference)] < 0),
                # Amendment 5: how close to the floor the WORSE of the two
                # inputs sat. The flag above is a hard `< 0` and no threshold
                # moves after seeing data — this is reported so a reader can
                # see a DiD built on a value that is technically above the
                # floor and practically at it. Block 221's CDX2 reference is
                # +0.003, a signal 0.3% above noise, and CDX2 is §7's decisive
                # row.
                "min_separation": float(min(sep[(gene, adenoma)],
                                            sep[(gene, reference)])),
                # The floor cancels here: log_sep is log(mu_gene) - log(mu_floor),
                # so this is dlog(mu_gene) - dlog(mu_control).
                "did": float(delta - control_delta),
                "adenoma_domain": adenoma, "reference_domain": reference,
            })
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    # One row per block: a block re-run supersedes its earlier run.
    frame = frame.sort_values("run").drop_duplicates(["block", "gene"], keep="last")
    return frame.sort_values(["gene", "block"], ignore_index=True)


def _interval(values: np.ndarray) -> tuple[float, float] | None:
    """Student-t, or None below ``MIN_STUDIES``."""
    n = int(values.size)
    if n < MIN_STUDIES:
        return None
    se = float(values.std(ddof=1) / np.sqrt(n))
    crit = float(stats.t.ppf(0.975, n - 1))
    return float(values.mean() - crit * se), float(values.mean() + crit * se)


def aggregate(per_block: pd.DataFrame) -> pd.DataFrame:
    """Student-t over blocks — and no interval below ``MIN_STUDIES``.

    Not the percentile bootstrap: at n=7 it is 0.742x the correct width and
    fires 12.0% under a true null (`docs/HANDOFF.md` §3a). Prereg §4.
    """
    rows = []
    for gene, group in per_block.groupby("gene", sort=True):
        values = group["did"].to_numpy(float)
        n = int(values.size)
        row: dict[str, object] = {
            "gene": gene, "role": GENE_ROLES.get(gene), "n_blocks": n,
            "mean_did": float(values.mean()) if n else float("nan"),
            "blocks": "; ".join(sorted(group["block"])),
            "all_same_sign": bool(n and (np.sign(values) == np.sign(values[0])).all()),
        }
        # Amendment 4's leave-out sensitivity, computed for every gene so the
        # consequence cannot be applied selectively.
        clean = group.loc[~group["target_below_floor"].fillna(False), "did"].to_numpy(float)
        row["n_blocks_below_floor"] = int(n - clean.size)
        row["min_separation_any_block"] = float(group["min_separation"].min())
        row["n_blocks_within_0p1_of_floor"] = int(
            (group["min_separation"] < 0.1).sum())
        clean_ci = _interval(clean)
        row["ci_low_excluding_floored"] = clean_ci[0] if clean_ci else None
        row["ci_high_excluding_floored"] = clean_ci[1] if clean_ci else None

        if n < MIN_STUDIES:
            # Invariant 1 one layer out: an interval that was refused is not an
            # interval of zero width.
            row.update({
                "ci_low": None, "ci_high": None, "interval": "refused",
                "why": (f"{n} block(s), below MIN_STUDIES={MIN_STUDIES}. Prereg "
                        f"§6: report the per-block values side by side and no "
                        f"interval."),
            })
        else:
            se = float(values.std(ddof=1) / np.sqrt(n))
            crit = float(stats.t.ppf(0.975, n - 1))
            row.update({
                "ci_low": float(values.mean() - crit * se),
                "ci_high": float(values.mean() + crit * se),
                "interval": "student_t",
                "se": se, "t_crit": crit,
                "width_vs_avenue_a_n43": float(
                    (crit / np.sqrt(n)) / (stats.t.ppf(0.975, 42) / np.sqrt(43))),
                "excludes_zero": bool(
                    (values.mean() - crit * se) * (values.mean() + crit * se) > 0),
                "excludes_zero_excluding_floored": (
                    bool(clean_ci[0] * clean_ci[1] > 0) if clean_ci else None),
                "why": ("pre-committed to be uninformative at n<5; the width "
                        "column is there so a reader sees it" if n < 5 else ""),
            })
        rows.append(row)
    return pd.DataFrame(rows)


def verdict(summary: pd.DataFrame) -> dict:
    """Prereg §7, whose decisive row is the DISCRIMINATOR, not the target."""
    if summary.empty:
        return {"verdict": "NOTHING TO AGGREGATE", "detail": "no block yielded a DiD."}
    by_gene = summary.set_index("gene")
    n = int(by_gene["n_blocks"].max())
    if TARGET_GENE not in by_gene.index:
        return {"verdict": "NOT ESTIMABLE",
                "detail": f"{TARGET_GENE} has no DiD in any block."}

    target = by_gene.loc[TARGET_GENE]
    disc = by_gene.loc[DISCRIMINATOR_GENE] if DISCRIMINATOR_GENE in by_gene.index else None

    if n < MIN_STUDIES:
        return {
            "verdict": f"BELOW THE INTERVAL FLOOR — {n} BLOCK(S), REPORTED SIDE BY SIDE",
            "detail": (
                f"{TARGET_GENE} mean DiD {target['mean_did']:+.3f} over {n} "
                f"block(s), all same sign: {bool(target['all_same_sign'])}. "
                f"**No interval, by §6.** "
                + (f"{DISCRIMINATOR_GENE} mean DiD {disc['mean_did']:+.3f} — §7's "
                   f"discriminator, and at this n it discriminates nothing yet. "
                   if disc is not None else "")
                + "Sign agreement across a handful of blocks is not evidence of "
                  "an effect size; it is what it says."
            ),
        }

    # Amendment 4: an interval whose sign turns on blocks where the target is
    # below the noise floor is INDETERMINATE, not a positive. Binding, not a
    # caveat — the Becker prereg records a check that warned and did not bind,
    # and it returned the most optimistic verdict in its table.
    floored = int(target.get("n_blocks_below_floor", 0) or 0)
    with_floored = bool(target.get("excludes_zero", False))
    without_floored = target.get("excludes_zero_excluding_floored")
    if floored and without_floored is not None and with_floored != without_floored:
        return {"verdict": "INDETERMINATE — THE SIGN TURNS ON BELOW-FLOOR BLOCKS",
                "detail": (
                    f"{TARGET_GENE} DiD {target['mean_did']:+.3f} "
                    f"[{target['ci_low']:+.3f}, {target['ci_high']:+.3f}] over "
                    f"{n} blocks "
                    f"{'excludes' if with_floored else 'includes'} zero; "
                    f"dropping the {floored} block(s) where the target sits "
                    f"below the negative-probe floor it "
                    f"{'excludes' if without_floored else 'includes'} it "
                    f"[{target['ci_low_excluding_floored']:+.3f}, "
                    f"{target['ci_high_excluding_floored']:+.3f}]. Below the "
                    f"floor the observed signal is smaller than what noise "
                    f"alone supplies, so that DiD is a BOUND and not a point, "
                    f"and a mean mixing bounds with points estimates neither. "
                    f"Amendment 4, fixed before this aggregate was computed. "
                    f"**Neither interval is the answer.**"
                )}

    # §7's THIRD ROW, which the code did not implement: "DiD(GUCA2A) > 0 in
    # most patients -> the single-block observation was an artefact. Report it
    # and stop." Without this branch a positive target would fall through to
    # "A TIER MOVED" or "TARGET FALLS...", both of which would be false
    # statements about the direction. Latent -- GUCA2A is negative in all six
    # blocks -- and a falsifier that cannot fire is not a falsifier.
    if bool(target.get("excludes_zero", False)) and target["mean_did"] > 0:
        return {"verdict": "TARGET RISES — §7's THIRD ROW; THE EARLIER BLOCKS WERE AN ARTEFACT",
                "detail": (
                    f"{TARGET_GENE} DiD {target['mean_did']:+.3f} "
                    f"[{target['ci_low']:+.3f}, {target['ci_high']:+.3f}] over "
                    f"{n} blocks excludes zero **POSITIVE** — the target rises "
                    f"against its control from reference to adenoma, which is "
                    f"the opposite of §7's prediction. §7: report it and stop. "
                    f"Any earlier block-level fall was not the effect."
                )}

    # Amendment 6: the same rule for the discriminator. §7 rests on TWO
    # conclusions and Amendment 4 protected only one of them.
    if disc is not None:
        d_floored = int(disc.get("n_blocks_below_floor", 0) or 0)
        d_with = bool(disc.get("excludes_zero", False))
        d_without = disc.get("excludes_zero_excluding_floored")
        if d_floored and d_without is not None and d_with != d_without:
            return {"verdict": "INDETERMINATE — THE DISCRIMINATOR TURNS ON "
                               "BELOW-FLOOR BLOCKS",
                    "detail": (
                        f"{DISCRIMINATOR_GENE} {disc['mean_did']:+.3f} "
                        f"[{disc['ci_low']:+.3f}, {disc['ci_high']:+.3f}] "
                        f"{'excludes' if d_with else 'includes'} zero; dropping "
                        f"the {d_floored} block(s) where it sits below the "
                        f"negative-probe floor it "
                        f"{'excludes' if d_without else 'includes'} it. §7's "
                        f"verdict rests on the discriminator as much as on the "
                        f"target, so a discriminator resting on a bound is the "
                        f"same defect. Amendment 6."
                    )}

    if not bool(target.get("excludes_zero", False)):
        return {"verdict": "NO CLAIM",
                "detail": (
                    f"{TARGET_GENE} DiD {target['mean_did']:+.3f} "
                    f"[{target['ci_low']:+.3f}, {target['ci_high']:+.3f}] over "
                    f"{n} blocks, includes zero. At {target.get('width_vs_avenue_a_n43', float('nan')):.2f}x "
                    f"avenue A's width this is WEAK evidence and **must not be "
                    f"quoted as a negative** (§7)."
                )}

    if disc is not None and bool(disc.get("excludes_zero", False)) and \
            np.sign(disc["mean_did"]) == np.sign(target["mean_did"]):
        return {"verdict": "A TIER MOVED, NOT A GENE",
                "detail": (
                    f"{TARGET_GENE} {target['mean_did']:+.3f} and "
                    f"{DISCRIMINATOR_GENE} {disc['mean_did']:+.3f} both exclude "
                    f"zero in the same direction. §7: the effect is not "
                    f"gene-specific. This is the same discriminator avenue A "
                    f"used, and it has fired."
                )}

    if disc is None:
        # The branch below reads the discriminator, and §7 makes it the decisive
        # row. Without it there is no verdict to give, and dereferencing None
        # here is what the first version did.
        return {"verdict": "TARGET FALLS, DISCRIMINATOR ABSENT — NO VERDICT",
                "detail": (
                    f"{TARGET_GENE} DiD {target['mean_did']:+.3f} "
                    f"[{target['ci_low']:+.3f}, {target['ci_high']:+.3f}] over "
                    f"{n} blocks excludes zero, but {DISCRIMINATOR_GENE} has no "
                    f"DiD in this run. §7 makes the discriminator the decisive "
                    f"row — a target falling alone is consistent with a whole "
                    f"tier moving — so there is nothing to conclude."
                )}
    # Amendment 7. The branch above requires the discriminator to clear zero in
    # the SAME direction as the target, so this one is reached in two different
    # states: the discriminator not clearing zero, and the discriminator
    # clearing it in the OPPOSITE direction. The first version asserted "does
    # not exclude zero" in both, which is false in the second — and CDX2's
    # lower bound has run -0.251, -0.156, -0.020 across four, five and six
    # blocks, so the second state is one block away.
    if bool(disc.get("excludes_zero", False)):
        return {"verdict": "TARGET FALLS AND THE DISCRIMINATOR MOVES THE OTHER WAY",
                "detail": (
                    f"{TARGET_GENE} DiD {target['mean_did']:+.3f} "
                    f"[{target['ci_low']:+.3f}, {target['ci_high']:+.3f}] over "
                    f"{n} blocks; {DISCRIMINATOR_GENE} "
                    f"{disc['mean_did']:+.3f} "
                    f"[{disc['ci_low']:+.3f}, {disc['ci_high']:+.3f}] excludes "
                    f"zero in the OPPOSITE direction. §7 asks whether a tier "
                    f"moved; a discriminator moving against the target is a "
                    f"stronger answer to that than one merely failing to clear "
                    f"zero. **It is not a stronger claim about mechanism** — "
                    f"Amendment 2 still applies, and a rising "
                    f"{DISCRIMINATOR_GENE} is equally consistent with more "
                    f"{DISCRIMINATOR_GENE}-positive cells in the lesion. "
                    f"Amendment 7."
                )}

    return {"verdict": "TARGET FALLS AND THE DISCRIMINATOR DOES NOT",
            "detail": (
                f"{TARGET_GENE} DiD {target['mean_did']:+.3f} "
                f"[{target['ci_low']:+.3f}, {target['ci_high']:+.3f}] over {n} "
                f"blocks; {DISCRIMINATOR_GENE} "
                f"{disc['mean_did']:+.3f} "
                f"[{disc['ci_low']:+.3f}, {disc['ci_high']:+.3f}] does not "
                f"exclude zero. "
                f"**Still not silencing and still not per-cell** — feasibility "
                f"§7 is untouched, and §8 of this document lists what remains "
                f"undecided."
            )}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    per_block = per_block_did(args.results_dir)
    if per_block.empty:
        log.error("no committed per-section table carried both a reference and "
                  "an adenoma domain.")
        return 4

    log.info("\n%s\n§4's DiD, PER BLOCK — Delta(gene) - Delta(%s)\n%s",
             "=" * 72, CONTROL_GENE, "=" * 72)
    wide = per_block.pivot(index="gene", columns="block", values="did")
    log.info("%s", wide.to_string(float_format=lambda v: f"{v:+8.3f}"))

    summary = aggregate(per_block)
    log.info("\n%s\nOVER BLOCKS\n%s", "=" * 72, "=" * 72)
    cols = [c for c in ("gene", "role", "n_blocks", "mean_did", "ci_low",
                        "ci_high", "interval", "all_same_sign") if c in summary]
    log.info("%s", summary[cols].to_string(index=False))

    outcome = verdict(summary)
    log.info("\n%s\nVERDICT\n%s", "=" * 72, "=" * 72)
    log.info("  %s\n  %s", outcome["verdict"], outcome["detail"])
    log.info(
        "\n  §4's control is %s, the ONLY control-role gene in this deposit, "
        "and it is an\n  epithelial keratin — it rises with epithelial FRACTION "
        "as well as with capture.\n  So this statistic does not separate "
        "'mature cells are gone' from 'mature cells\n  are present and "
        "silenced'. See multisection Amendment 2.", CONTROL_GENE)

    meta = {
        "prereg": "docs/prereg_crowell_multisection.md",
        "statistic": f"DiD(gene) = Delta(gene) - Delta({CONTROL_GENE}), "
                     f"on log_separation; the false-positive floor cancels",
        "control_gene": CONTROL_GENE,
        "control_is_also_epithelial": True,
        "discriminator": DISCRIMINATOR_GENE,
        "unit": "tissue block, not section (invariant 5; 231 and 232 are one)",
        "interval": f"Student-t over blocks; refused below {MIN_STUDIES}",
        "does_not_separate": (
            "compositional from intrinsic — the control rises with epithelial "
            "fraction. Amendment 2."
        ),
        "verdict": outcome,
        "exploratory": False,
        "pre_registered": True,
    }
    for frame, name in ((per_block, "crowell_multisection_per_block"),
                        (summary, "crowell_multisection_summary")):
        log.info("wrote %s", write_versioned_table(
            frame, name, seed=args.seed, results_dir=args.results_dir,
            allow_dirty=args.allow_dirty, extra_meta=meta,
        ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
