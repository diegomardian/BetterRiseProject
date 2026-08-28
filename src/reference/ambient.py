"""Ambient RNA: measuring it, and correcting it. W1, week 2.

GSE178341 as deposited has had its empty droplets removed (docs/open_decisions.md
#8), so the usual route to a soup profile — average the cell-free droplets — is
unavailable. CellBender cannot run at all. SoupX runs only if handed a profile
estimated some other way, and DecontX infers contamination from cluster
structure.

That makes the *measurement* in this module more important than the correction.
:func:`contamination_fraction` estimates ambient contamination **without any
empty droplets**, using genes that cannot be expressed in the cells you are
looking at. If a colonocyte shows haemoglobin counts, no colonocyte made them —
they came from the soup, and their size relative to the soup's own composition
is the contamination rate.

Why that matters here specifically. Both surviving correction methods infer the
ambient profile from the cells themselves, which risks absorbing the very signal
this project measures: a rare mature cell's genuine low-level GUCA2A is exactly
what a cluster-based method may reassign to "contamination from elsewhere". An
estimate built from impossible genes is independent of that machinery, so it can
audit the correction rather than inherit its assumptions. It also feeds gate
criterion G1 directly.

Matrix orientation is **cells x genes**, as in ``ingest.py`` and ``qc.py``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Final

import numpy as np
import pandas as pd

from src.common.panel import panel_genes
from src.reference.signature import assert_no_target_leakage

#: Genes that a cell of the given compartment cannot transcribe. Residual counts
#: in these are ambient by construction.
#:
#: Chosen to be lineage-exclusive and highly expressed in *some* population, so
#: they are well represented in the soup and give a stable ratio. Haemoglobin is
#: the workhorse — erythrocytes are abundant, lyse readily, and no epithelial
#: cell makes HBB.
#:
#: None of these may be panel genes; asserted at call time (CLAUDE.md invariant
#: 2). Using a target gene to measure contamination would be circular.
#: The epithelial set is deliberately broad. On the pilot, the original six
#: (haemoglobin, PTPRC, IGKC, IGHG1) held under 0.2% of the soup in 10 of 23
#: samples, because their share is dominated by how many erythrocytes and plasma
#: cells that particular dissociation happened to contain — so those samples came
#: back unestimable. The pan-leukocyte transcripts below are present in every
#: sample and abundant, which stabilises the denominator.
#:
#: Chosen for being **constitutively** haematopoietic. Deliberately excluded:
#: CD74 and HLA-DRA (interferon-inducible in epithelium), and the fibroblast
#: collagens (COL1A1/COL1A2/COL3A1), which are abundant and tempting but can be
#: switched on in partial-EMT tumour cells — exactly the state this project
#: studies, so using them would put the measurement inside the phenomenon.
IMPOSSIBLE_GENES: Final[dict[str, frozenset[str]]] = {
    "epithelial": frozenset(
        {
            # erythroid
            "HBB", "HBA1", "HBA2",
            # plasma / B
            "IGKC", "IGHG1", "MS4A1",
            # pan-leukocyte, constitutive and highly expressed
            "PTPRC", "SRGN", "LAPTM5", "CORO1A", "CD52", "ARHGDIB",
            # T
            "CD3D", "CD3E", "CD2",
            # myeloid
            "TYROBP", "FCER1G", "AIF1",
        }
    ),
    "immune": frozenset({"HBB", "HBA1", "HBA2", "EPCAM", "KRT8", "KRT18", "KRT19"}),
    "stromal": frozenset({"HBB", "HBA1", "HBA2", "EPCAM", "PTPRC", "KRT8"}),
    "endothelial": frozenset({"HBB", "HBA1", "HBA2", "EPCAM", "PTPRC", "KRT8"}),
}


class AmbientError(ValueError):
    """Contamination could not be estimated from the input given."""


def _row_sums(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel()


def _gene_positions(gene_names: Any, wanted: Any) -> np.ndarray:
    names = list(gene_names)
    lookup = {str(n): i for i, n in enumerate(names)}
    return np.array([lookup[g] for g in wanted if g in lookup], dtype=int)


def soup_profile_from_cells(matrix: Any, gene_names: Any) -> pd.Series:
    """Per-gene share of total counts, pooled over all cells.

    A stand-in for the empty-droplet profile. The soup is largely mRNA from
    cells that lysed during dissociation, so the tissue's overall composition
    approximates it — imperfectly, because it also contains the real biological
    signal we are trying to measure.

    This is what you pass to SoupX's ``setSoupProfile()`` when there is no raw
    droplet table to call ``estimateSoup()`` on. Compute it **per sample**: the
    soup is a property of one dissociation, not of the study.

    Pass the **whole sample, every compartment** — never a single population.
    The impossible-gene estimator divides the observed rate by this profile's
    share of those genes, and if the profile is built only from the cells being
    tested, both terms are the same number and the ratio degenerates to 1.0.
    The erythrocytes and leukocytes that actually transcribe HBB and PTPRC are
    what put those genes in the profile at all.

    Because the pooled average includes the tested cells, the estimate is biased
    slightly **upward** in proportion to how much of the sample they represent —
    they dilute the profile at exactly the genes the ratio depends on. Negligible
    when the tested population is a minority; worth stating when it is not.
    """
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    grand = totals.sum()
    if grand <= 0:
        raise AmbientError("matrix has no counts; cannot build a soup profile")
    return pd.Series(totals / grand, index=list(gene_names), name="soup_fraction")


def contamination_fraction(
    matrix: Any,
    gene_names: Any,
    *,
    cell_mask: Any,
    impossible: Any = None,
    compartment: str = "epithelial",
    soup: pd.Series | None = None,
    min_soup_share: float = 0.002,
    max_plausible: float = 0.95,
) -> float:
    """Ambient contamination in the masked cells, from impossible genes alone.

    The estimator, which is the one SoupX uses for its marker-based route::

        observed = counts of impossible genes in these cells / their total counts
        expected = share of those genes in the soup, if the cell were pure soup
        rho      = observed / expected

    A cell of this compartment cannot make these transcripts, so everything
    observed is ambient. Dividing by the soup's own share of them converts that
    into "what fraction of this cell's counts came from soup".

    Parameters
    ----------
    cell_mask:
        Boolean, selecting cells known to belong to `compartment`. Get it from a
        malignancy or lineage call — not from the panel.
    soup:
        Soup profile over the same genes. Defaults to
        :func:`soup_profile_from_cells` on the whole matrix, which is the
        approximation available without empty droplets. The default is
        deliberately the **whole matrix**, not the masked cells — see that
        function's note on why a single-compartment profile degenerates.
    min_soup_share:
        Refuse to estimate when the impossible genes hold less than this share of
        the soup — they cannot discriminate if they are barely in it.
    max_plausible:
        Refuse when the raw ratio reaches this. A ratio at or above 1 says the
        masked cells carry as much impossible-gene signal as pure soup would
        give, which is not contamination but a **violated assumption**: doublets,
        or cells misassigned to this compartment. Observed on C107_T_1_1_0, which
        returned exactly 1.000 with 543 epithelial cells in a 2,737-cell sample.
        Clipping that to 1.0 and calling it "100% ambient" would be a fabricated
        number; refusing is the honest answer.

    Returns a fraction in [0, 1]. Report it per sample and log it next to the
    correction — it is the number that says whether the correction was
    necessary and whether it went far enough.
    """
    genes = frozenset(impossible) if impossible is not None else IMPOSSIBLE_GENES.get(compartment)
    if not genes:
        raise AmbientError(
            f"no impossible-gene set for compartment {compartment!r}; "
            f"known: {sorted(IMPOSSIBLE_GENES)}"
        )

    mask = np.asarray(cell_mask, dtype=bool)
    if mask.shape[0] != matrix.shape[0]:
        raise AmbientError(f"cell_mask has {mask.shape[0]} entries for {matrix.shape[0]} cells")
    if not mask.any():
        raise AmbientError("cell_mask selects no cells")

    positions = _gene_positions(gene_names, sorted(genes))
    if positions.size == 0:
        raise AmbientError(
            f"none of the impossible genes {sorted(genes)} are in the matrix. "
            f"Check gene naming — symbols vs Ensembl IDs (open decision #3)."
        )

    # AFTER the naming check, not before. Using a target gene to measure
    # contamination would be circular — but a gene set that is not in the matrix
    # at all is a naming bug, and that diagnosis is the more useful one. Checking
    # an invariant on identifiers the matrix has never seen tests nothing.
    assert_no_target_leakage(genes, panel_genes(), context="the impossible-gene set")

    if soup is None:
        soup = soup_profile_from_cells(matrix, gene_names)
    expected = float(np.asarray(soup)[positions].sum())
    if expected <= 0:
        raise AmbientError(
            "the impossible genes have zero share of the soup, so contamination "
            "is not identifiable from them. Pick genes that are abundant "
            "somewhere in the tissue."
        )
    if expected < min_soup_share:
        raise AmbientError(
            f"the impossible genes hold only {expected:.4%} of the soup, below the "
            f"{min_soup_share:.4%} floor. This happens when the masked population "
            f"dominates the sample: the pooled profile is then mostly those same "
            f"cells, so observed and expected converge and the ratio runs to 1.0 "
            f"regardless of the truth. Estimate the soup where the source "
            f"populations are actually present, or report this sample as not "
            f"estimable — do not report the 1.0."
        )

    selected = matrix[mask, :]
    observed_counts = float(np.asarray(selected[:, positions].sum()).ravel().sum())
    total_counts = float(_row_sums(selected).sum())
    if total_counts <= 0:
        raise AmbientError("masked cells have no counts")

    ratio = (observed_counts / total_counts) / expected
    if ratio >= max_plausible:
        raise AmbientError(
            f"raw ratio {ratio:.3f} means the masked cells hold as much "
            f"impossible-gene signal as pure soup would give. That is not a "
            f"contamination level — it says the impossible-gene assumption is "
            f"violated for these cells (doublets, or a wrong compartment call). "
            f"Report this sample as not estimable rather than as ~100% ambient."
        )
    return float(max(ratio, 0.0))


def contamination_by_sample(
    matrix: Any,
    gene_names: Any,
    *,
    sample_id: Any,
    cell_mask: Any,
    compartment: str = "epithelial",
    min_soup_share: float = 0.002,
) -> pd.DataFrame:
    """Per-sample contamination estimate. One row per sample.

    Per sample, not pooled: the soup is a property of a single dissociation, and
    a study-wide number would average away exactly the variation that decides
    which samples to distrust.
    """
    samples = pd.Series(list(sample_id), name="sample_id")
    if len(samples) != matrix.shape[0]:
        raise AmbientError(f"sample_id has {len(samples)} entries for {matrix.shape[0]} cells")
    mask = np.asarray(cell_mask, dtype=bool)

    rows: list[dict[str, Any]] = []
    for sample, index in samples.groupby(samples).groups.items():
        positions = np.asarray(index, dtype=int)
        in_sample = np.zeros(matrix.shape[0], dtype=bool)
        in_sample[positions] = True
        selected = in_sample & mask
        if not selected.any():
            rows.append(
                {"sample_id": sample, "n_cells": 0, "mask_share": 0.0,
                 "soup_share": np.nan, "contamination": np.nan}
            )
            continue
        # Soup estimated within the sample — see the docstring above.
        soup = soup_profile_from_cells(matrix[in_sample, :], gene_names)
        genes = frozenset(IMPOSSIBLE_GENES.get(compartment, ()))
        share = float(np.asarray(soup)[_gene_positions(gene_names, sorted(genes))].sum())
        row = {
            "sample_id": sample,
            "n_cells": int(selected.sum()),
            "mask_share": float(selected.sum() / max(int(in_sample.sum()), 1)),
            "soup_share": share,
        }
        try:
            row["contamination"] = contamination_fraction(
                matrix,
                gene_names,
                cell_mask=selected,
                compartment=compartment,
                soup=soup,
                min_soup_share=min_soup_share,
            )
        except AmbientError:
            # Not estimable in this sample. None, never a number — the same
            # distinction CLAUDE.md invariant 1 makes for the intrinsic term.
            row["contamination"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values("sample_id", ignore_index=True)


#: Contamination above which a sample leaves the compositional arm.
#: **Open decision #16, committed 2026-08-23 before counting what it removes.**
#: Above roughly a tenth of counts being ambient, the per-cell marker detection
#: that axis 1's maturity call depends on is materially perturbed — that call is
#: a detection gate at a matched depth of ~3,300 UMIs, so ~330 ambient counts is
#: not a rounding error in it. Having chosen not to correct, a sample whose
#: ambient share approaches the effects being estimated cannot be rescued by a
#: caveat.
MAX_CONTAMINATION: Final[float] = 0.10


def ambient_exclusions(
    contamination: pd.DataFrame, *, threshold: float = MAX_CONTAMINATION
) -> pd.DataFrame:
    """Which samples leave the compositional arm, and what that costs.

    **The threshold was fixed before this was written** (#16). A threshold
    chosen after seeing how many patients it removes is not a threshold, so this
    reports the cost rather than offering to revise the number.

    Samples with no estimate are **kept**, not excluded. An unmeasurable
    contamination is not a high one, and treating "we could not tell" as "too
    dirty" would silently remove the sparsest samples — which are already the
    ones least able to spare cells.

    Returns one row per sample with `excluded` and the reason.
    """

    missing = {"sample_id", "contamination"} - set(contamination.columns)
    if missing:
        raise AmbientError(f"contamination is missing column(s): {sorted(missing)}")

    out = contamination.loc[:, ["sample_id", "contamination"]].copy()
    measured = out["contamination"].notna()
    out["excluded"] = measured & (out["contamination"] > threshold)
    out["reason"] = np.where(
        out["excluded"], f"contamination above {threshold:.0%}",
        np.where(measured, "", "no estimate — kept, not excluded"),
    )
    return out


def differential_contamination(
    contamination: pd.DataFrame, *, warn_at: float = 0.05
) -> pd.DataFrame:
    """Is the soup dirtier in tumour than in matched normal? **Per patient.**

    The third instance of one question. :func:`src.reference.qc.
    differential_retention` asks it of the mitochondrial cap and
    ``labels.differential_resolution`` asks it of the depth floor: **anything
    that affects the two arms unequally moves Delta(mature fraction), because
    that difference IS the compositional term.**

    Ambient RNA is the same shape of problem and worse in one respect. Soup is
    enriched for whatever is abundant in the dissociation, and a tumour sample's
    abundant populations differ from its matched normal's — so an asymmetry here
    is expected rather than surprising, and it does not cancel.

    Which way it biases depends on what is in the soup. Ambient *mature* marker
    transcripts give immature cells false mature counts; ambient stem markers do
    the reverse. Either way a patient whose tumour arm is materially dirtier
    than their normal arm has a Delta that is partly a dissociation artifact.

    Returns one row per patient with both arms' contamination, the gap, and
    `flagged`. **If this flags widely, open decision #16's exclusion rule should
    be per-patient on the gap rather than per-sample on the level** — a patient
    whose two arms are equally dirty is far less compromised than one whose arms
    differ by ten points, even at a higher absolute level.
    """

    required = {"sample_id", "contamination"}
    missing = required - set(contamination.columns)
    if missing:
        raise AmbientError(f"contamination is missing column(s): {sorted(missing)}")

    frame = contamination.copy()
    if "tissue" not in frame.columns or "patient_id" not in frame.columns:
        from src.reference.ingest import parse_barcode

        parsed = [parse_barcode(f"{s}_id-AAAA") for s in frame["sample_id"]]
        frame["patient_id"] = [p["patient_id"] for p in parsed]
        frame["tissue"] = [p["tissue"] for p in parsed]

    wide = frame.pivot_table(
        index="patient_id", columns="tissue", values="contamination",
        aggfunc="median", observed=True,
    )
    for arm in ("tumour", "normal"):
        if arm not in wide.columns:
            wide[arm] = float("nan")

    out = wide.loc[:, ["tumour", "normal"]].copy()
    out.columns = ["contamination_tumour", "contamination_normal"]
    out["difference"] = out["contamination_tumour"] - out["contamination_normal"]
    out["flagged"] = out["difference"].abs() > warn_at
    # Both arms needed. A patient with one arm has no compositional term at all
    # (open decision #9), so their asymmetry is undefined rather than zero.
    out["both_arms"] = (
        out["contamination_tumour"].notna() & out["contamination_normal"].notna()
    )
    return out.reset_index()


# ---------------------------------------------------------------------------
# Correction
# ---------------------------------------------------------------------------


def compare_retention(
    soupx: pd.DataFrame, decontx: pd.DataFrame, *, sample_id: str = ""
) -> pd.DataFrame:
    """Join two methods' per-gene retention. **One row per gene per sample.**

    The week-2 deliverable is "per-gene retention table; correlation between
    methods reported" — a comparison, not a winner.

    **Per sample, not averaged across them.** A gene that loses 40% of its
    counts in one dissociation and 2% in another is a different finding from one
    that consistently loses 20%, and a cohort mean erases the difference. The
    soup is a property of a single dissociation, so its effect on a given gene
    is too.

    Retention is `counts_after / counts_before`, so **1.0 means untouched and
    lower means more removed**. A gene both methods strip hard is soup by
    agreement; a gene only one strips is a disagreement worth looking at before
    either number is trusted.
    """

    for name, frame in (("soupx", soupx), ("decontx", decontx)):
        missing = {"gene", "retention"} - set(frame.columns)
        if missing:
            raise AmbientError(f"{name} is missing column(s): {sorted(missing)}")

    out = soupx[["gene", "retention"]].rename(
        columns={"retention": "retention_soupx"}
    ).merge(
        decontx[["gene", "retention"]].rename(
            columns={"retention": "retention_decontx"}
        ),
        on="gene", how="inner",
    )
    if out.empty:
        raise AmbientError(
            "no genes in common between the two retention tables — check both "
            "ran on the same matrix"
        )
    out["difference"] = out["retention_soupx"] - out["retention_decontx"]
    out["sample_id"] = sample_id
    return out


def retention_agreement(comparison: pd.DataFrame) -> dict[str, Any]:
    """How much the two methods agree, per sample. Spearman, not Pearson.

    Retention is bounded at 1 and piles up there — most genes are barely
    touched — so the distribution is heavily skewed and a Pearson correlation
    would be dominated by a handful of hard-stripped genes. Rank correlation
    asks the question that matters: **do the two methods strip the same genes
    hardest**, whatever the absolute amounts.

    `median_difference` carries the sign separately, because two methods can
    rank genes identically while one removes twice as much.
    """

    frame = comparison.dropna(subset=["retention_soupx", "retention_decontx"])
    if len(frame) < 3:
        raise AmbientError(
            f"only {len(frame)} gene(s) with both estimates — too few to "
            f"correlate"
        )
    spearman = float(
        frame["retention_soupx"].corr(frame["retention_decontx"], method="spearman")
    )
    return {
        "sample_id": frame["sample_id"].iloc[0] if "sample_id" in frame else "",
        "n_genes": int(len(frame)),
        "spearman": spearman,
        "median_retention_soupx": float(frame["retention_soupx"].median()),
        "median_retention_decontx": float(frame["retention_decontx"].median()),
        "median_difference": float(frame["difference"].median()),
        # Agreement on WHICH genes, not on how much. Two methods can rank
        # identically while one removes twice as much, and that is a different
        # kind of disagreement from ranking them differently.
        "agree": bool(spearman > 0.5),
    }


def _write_sparse(matrix, gene_names, barcodes, out_dir):
    """matrix.mtx / genes.tsv / barcodes.tsv, genes x cells. Sparse throughout.

    Dense would be the size of the expression matrix per sample, written so it
    can be read straight back — the same waste the inferCNV path already
    avoided.
    """
    from scipy import io as sio
    from scipy import sparse

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    names = np.asarray([str(g) for g in gene_names], dtype=object)
    _, first = np.unique(names, return_index=True)
    keep = np.zeros(names.shape, dtype=bool)
    keep[np.sort(first)] = True   # one row per symbol; SoupX keys on rownames

    subset = matrix[:, keep]
    m = subset.T if sparse.issparse(subset) else sparse.csr_matrix(
        np.asarray(subset).T
    )
    sio.mmwrite(str(out / "matrix.mtx"), m.tocoo())
    (out / "genes.tsv").write_text("".join(f"{g}\n" for g in names[keep]))
    (out / "barcodes.tsv").write_text(
        "".join(f"{b}\n" for b in map(str, barcodes))
    )
    return out, names[keep]


def run_soupx(
    matrix: Any,
    gene_names: Any,
    *,
    barcodes: Any,
    clusters: Any,
    soup_profile: Any,
    out_dir: Any,
    contamination: float | None = None,
    dry_run: bool = False,
) -> dict[str, Any]:
    """SoupX in degraded mode, for ONE sample. W1, week 2.

    **Degraded because there is no choice.** SoupX normally estimates the soup
    from empty droplets via ``estimateSoup()``, and open decision #8 established
    that GSE178341 ships none — GEO ran dropletUtils upstream. The channel is
    therefore built with ``calcSoupProfile = FALSE`` and handed a profile from
    :func:`soup_profile_from_cells` through ``setSoupProfile()``. That is a
    supported path, not a workaround, but it does mean the profile comes from
    cells rather than from ambient droplets and is only as good as the
    assumption that highly-expressed genes dominate both.

    Per sample, never pooled: the soup is a property of one dissociation.

    **The output is a per-gene retention table, not a corrected matrix.** Open
    decision #16 committed to measuring rather than correcting — at a 2.2%
    cohort median, DecontX-style correction risks absorbing genuine low-level
    marker expression in a rare population, which is this project's signal. So
    ``adjustCounts()`` runs, retention is computed from it, and the corrected
    matrix is discarded. Writing 62 of them would also not fit on the project
    filesystem.

    `contamination`
        Pass a value to fix rho instead of calling ``autoEstCont()``. Leave it
        None and SoupX estimates rho itself from clusters and marker genes.
        **Compare whatever comes back against
        :func:`contamination_fraction`** — they estimate the same quantity by
        unrelated routes, and disagreement is a result rather than something to
        reconcile.
    """
    out, kept = _write_sparse(matrix, gene_names, barcodes, out_dir)

    labels = np.asarray([str(c) for c in clusters], dtype=object)
    if labels.shape[0] != matrix.shape[0]:
        raise AmbientError(
            f"clusters has {labels.shape[0]} entries for {matrix.shape[0]} cells"
        )
    if len(set(labels)) < 2 and contamination is None:
        raise AmbientError(
            "autoEstCont needs more than one cluster to find marker genes. "
            "Pass contamination= to fix rho, or supply real cluster labels."
        )
    (out / "clusters.tsv").write_text("".join(f"{c}\n" for c in labels))

    # Collapse duplicate symbols BEFORE reindexing. soup_profile_from_cells
    # returns a Series indexed by gene symbol, and this deposit maps several
    # Ensembl IDs onto one symbol — so the index has duplicates and reindex()
    # refuses outright. Summing is the right collapse: the profile is a share
    # of the soup, and two rows for one symbol are two parts of the same share.
    profile = pd.Series(soup_profile)
    if profile.index.has_duplicates:
        profile = profile.groupby(level=0).sum()
    profile = profile.reindex(kept).fillna(0.0)
    if float(profile.sum()) <= 0:
        raise AmbientError(
            "the soup profile is empty over the genes in this matrix — "
            "setSoupProfile would have nothing to work with"
        )
    profile.to_csv(out / "soup_profile.csv", header=["est"], index_label="gene")

    script = out / "run_soupx.R"
    script.write_text(_SOUPX_R_TEMPLATE.format(
        out_dir=out,
        rho=("NULL" if contamination is None else f"{float(contamination)}"),
    ))
    command = ["Rscript", str(script)]
    result: dict[str, Any] = {"command": command, "out_dir": out, "script": script}
    if dry_run:
        result["ran"] = False
        return result

    import subprocess

    lines: list[str] = []
    with (out / "soupx_R.log").open("w") as log:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            lines.append(line)
        code = process.wait()
    result["ran"] = True
    result["returncode"] = code
    if code != 0:
        raise AmbientError(
            f"SoupX failed (exit {code}). Log: {out / 'soupx_R.log'}\n"
            + "".join(lines[-30:])
        )
    return result


#: The R side. Kept beside the reasoning rather than as a separate .R file, for
#: the same reason as the inferCNV template: the settings ARE the judgement
#: calls, and separating them from the docstring is how the two drift.
_SOUPX_R_TEMPLATE = """\
# Generated by src/reference/ambient.py:run_soupx — do not edit by hand.
suppressPackageStartupMessages({{library(SoupX); library(Matrix)}})

counts <- Matrix::readMM(file.path("{out_dir}", "matrix.mtx"))
counts <- as(counts, "CsparseMatrix")
rownames(counts) <- readLines(file.path("{out_dir}", "genes.tsv"))
colnames(counts) <- readLines(file.path("{out_dir}", "barcodes.tsv"))
clusters <- readLines(file.path("{out_dir}", "clusters.tsv"))
names(clusters) <- colnames(counts)

# calcSoupProfile = FALSE: there are no empty droplets to estimate from
# (open decision #8), so the profile is supplied from cells instead.
sc <- SoupChannel(counts, counts, calcSoupProfile = FALSE)
prof <- read.csv(file.path("{out_dir}", "soup_profile.csv"), row.names = 1)
prof$counts <- prof$est * sum(counts)
sc <- setSoupProfile(sc, prof)
sc <- setClusters(sc, clusters)

rho <- {rho}
if (is.null(rho)) {{
  sc <- autoEstCont(sc, doPlot = FALSE)
  rho <- mean(sc$metaData$rho)
}} else {{
  sc <- setContaminationFraction(sc, rho)
}}
cat("rho:", signif(rho, 4), "\n")

adjusted <- adjustCounts(sc, roundToInt = FALSE)

# Per-gene retention, then throw the corrected matrix away. Decision #16 is to
# MEASURE, not correct, and 62 corrected matrices would not fit on this
# filesystem anyway.
before <- Matrix::rowSums(counts)
after  <- Matrix::rowSums(adjusted)
write.csv(
  data.frame(
    gene = rownames(counts),
    counts_before = as.numeric(before),
    counts_after  = as.numeric(after),
    retention     = as.numeric(ifelse(before > 0, after / before, NA))
  ),
  file = file.path("{out_dir}", "soupx_retention.csv"), row.names = FALSE
)
cat("wrote soupx_retention.csv for", nrow(counts), "genes\n")
"""


def run_decontx(
    matrix: Any,
    gene_names: Any,
    *,
    barcodes: Any,
    clusters: Any,
    out_dir: Any,
    seed: int = 20260101,
    dry_run: bool = False,
) -> dict[str, Any]:
    """DecontX (``bioconductor-celda``), for ONE sample. W1, week 2.

    **The second method, and it exists because CellBender cannot run here.**
    execution_plan.md originally asked for "SoupX and CellBender, both,
    compared". CellBender learns the ambient profile from *empty droplets* and
    open decision #8 established that GSE178341 ships none — the deposit is
    already cell-called. DecontX needs no empties at all: it models each cell as
    a mixture of its own cluster's expression distribution and contamination
    from the other clusters.

    That independence from SoupX is the point. SoupX here runs in degraded mode
    off a profile computed from cells, so the two methods share an assumption
    only in the loosest sense — which makes their **agreement informative and
    their disagreement more so**. §4 asks for the comparison, not a winner.

    **The caveat that has to travel with every DecontX number.** Because it
    defines contamination as counts resembling *other clusters*, it can absorb
    genuine low-level expression of a marker in a rare population — which is
    exactly this project's signal. A mature marker expressed weakly in a small
    surviving population is indistinguishable, to this model, from soup. That is
    why open decision #16 reports rather than corrects, and why gate criterion
    G1 exists as an independent check rather than a confirmation.

    Emits a per-gene retention table and DecontX's own per-cell contamination
    estimate; the corrected matrix is computed and discarded, as in
    :func:`run_soupx`.
    """
    out, kept = _write_sparse(matrix, gene_names, barcodes, out_dir)

    labels = np.asarray([str(c) for c in clusters], dtype=object)
    if labels.shape[0] != matrix.shape[0]:
        raise AmbientError(
            f"clusters has {labels.shape[0]} entries for {matrix.shape[0]} cells"
        )
    if len(set(labels)) < 2:
        raise AmbientError(
            "DecontX defines contamination as counts resembling OTHER "
            "clusters, so one cluster leaves it nothing to compare against. "
            "Supply real cluster labels."
        )
    (out / "clusters.tsv").write_text("".join(f"{c}\n" for c in labels))

    script = out / "run_decontx.R"
    script.write_text(_DECONTX_R_TEMPLATE.format(out_dir=out, seed=int(seed)))
    command = ["Rscript", str(script)]
    result: dict[str, Any] = {"command": command, "out_dir": out, "script": script}
    if dry_run:
        result["ran"] = False
        return result

    import subprocess

    lines: list[str] = []
    with (out / "decontx_R.log").open("w") as log:
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            print(line, end="", flush=True)
            log.write(line)
            lines.append(line)
        code = process.wait()
    result["ran"] = True
    result["returncode"] = code
    if code != 0:
        raise AmbientError(
            f"DecontX failed (exit {code}). Log: {out / 'decontx_R.log'}\n"
            + "".join(lines[-30:])
        )
    return result


_DECONTX_R_TEMPLATE = """\
# Generated by src/reference/ambient.py:run_decontx — do not edit by hand.
suppressPackageStartupMessages({{library(celda); library(Matrix)}})
set.seed({seed})

counts <- Matrix::readMM(file.path("{out_dir}", "matrix.mtx"))
counts <- as(counts, "CsparseMatrix")
rownames(counts) <- readLines(file.path("{out_dir}", "genes.tsv"))
colnames(counts) <- readLines(file.path("{out_dir}", "barcodes.tsv"))
clusters <- readLines(file.path("{out_dir}", "clusters.tsv"))

# No background argument: DecontX needs no empty droplets, which is the whole
# reason it replaces CellBender here (open decision #8).
res <- decontX(counts, z = clusters, verbose = FALSE)
cat("mean per-cell contamination:", signif(mean(res$contamination), 4), "\n")

write.csv(
  data.frame(cell = colnames(counts), contamination = res$contamination),
  file = file.path("{out_dir}", "decontx_contamination.csv"), row.names = FALSE
)

# Per-gene retention, then discard the corrected matrix — decision #16 is to
# measure, not correct, and 62 of these would not fit on the filesystem.
before <- Matrix::rowSums(counts)
after  <- Matrix::rowSums(res$decontXcounts)
write.csv(
  data.frame(
    gene = rownames(counts),
    counts_before = as.numeric(before),
    counts_after  = as.numeric(after),
    retention     = as.numeric(ifelse(before > 0, after / before, NA))
  ),
  file = file.path("{out_dir}", "decontx_retention.csv"), row.names = FALSE
)
cat("wrote decontx_retention.csv for", nrow(counts), "genes\n")
"""
