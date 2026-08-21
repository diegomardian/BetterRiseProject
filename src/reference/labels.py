"""Differentiation labels: axes 1 and 2, four granularity rungs. W1, weeks 3-4.

What this produces, and why it looks the way it does
----------------------------------------------------
The decomposition asks "what fraction of cells are mature, and how much does each
make". Answering it requires deciding which cells are mature **without using the
gene under test** (CLAUDE.md invariant 2), and the answer is not unique: it
depends on which gene set defines maturity (the *axis*) and how finely the
epithelium is partitioned (the *rung*). README design decision 3 is explicit that
this is the point — "a compositional change becomes an expression change purely by
re-drawing cluster boundaries" — so the split is reported as a curve across
resolutions rather than as one number.

So labels are a **grid**: two transcript-based axes x four rungs = eight label
columns, all coexisting, none overwriting another. `label_{axis}_{rung}`.

    label_stem_pole_epithelial          label_opposite_lineage_epithelial
    label_stem_pole_lineage             label_opposite_lineage_lineage
    label_stem_pole_crypt_position      label_opposite_lineage_crypt_position
    label_stem_pole_best4               label_opposite_lineage_best4

Axes 3 (chromatin, spatial) are not transcript-based and are week 13+.

How maturity is scored
----------------------
Each axis has a frozen marker set (`config/labeling_axes.yaml`). Both axes score
*immaturity*, for different reasons, so maturity is the negated score:

- **stem_pole** (LGR5, ASCL2, MKI67, OLFM4, SMOC2): high = close to the stem pole,
  therefore less differentiated.
- **opposite_lineage** (MUC2, TFF3, SPDEF, ITLN1): high = goblet/secretory,
  therefore not on the absorptive maturation path this project follows.

The second is a weaker maturity proxy than the first, and deliberately so. The
claim in README design decision 2 is *agreement across structurally different
axes*, not that either axis is correct. Divergence between them is a finding.

What is frozen and what is not
------------------------------
The **axes are frozen** and come from `src.common.panel.axis_genes()` — never
retyped here. The **rung definitions below are W1's proposal** and are not frozen:
how many bins, and which bin counts as mature, are modelling choices. They are
parameterised and documented for that reason. Changing them changes the split by
design, which is the analysis, not a bug.

`best4` is marker-gated rather than score-binned, because BEST4+ is a discrete
population rather than a point on a gradient. BEST4 and SPIB are deliberately
absent from the frozen panel (execution_plan.md §3.2's sequencing constraint) so
they remain available as labels.

Invariant 2 is enforced on every marker set at call time, including W1's own rung
markers — not only on the frozen axes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Final

import numpy as np
import pandas as pd

from src.common.panel import axis_genes, granularity_rungs
from src.reference.signature import assert_no_target_leakage

#: Transcript-based axes. Chromatin and spatial are axis 3 and carry no markers.
TRANSCRIPT_AXES: Final[tuple[str, ...]] = ("stem_pole", "opposite_lineage")

#: Epithelial cells too shallow to score comparably once depth is matched.
#: Kept as their own label rather than folded into the least-mature bin: a cell
#: that could not be measured is not a cell measured to be immature, and open
#: decision #14 turns on exactly that distinction.
UNRESOLVED: Final[str] = "unresolved_depth"

#: Label given to cells outside the epithelium. They are not scored — a maturity
#: value for a fibroblast is meaningless — but they keep a row so the label frame
#: aligns with the AnnData it came from.
NON_EPITHELIAL: Final[str] = "non_epithelial"

#: BEST4+ markers. W1's proposal, not frozen. BEST4 and SPIB are absent from the
#: panel on purpose so they can be used here; CFTR and HES4 are co-expressed in
#: the same population and are also panel-clean.
BEST4_MARKERS: Final[frozenset[str]] = frozenset({"BEST4", "SPIB", "CFTR", "HES4"})

#: Fraction of epithelium called BEST4+ at that rung. Tier A's role note says the
#: BEST4+ program is "under 5% of epithelium even in healthy colon", so the gate
#: is a quantile rather than an absolute score, applied per sample.
BEST4_QUANTILE: Final[float] = 0.95


@dataclass(frozen=True)
class RungSpec:
    """How one granularity rung partitions the epithelium.

    `bins` are ordered least- to most-mature, so `bins[-1]` is the mature end.
    `markers` is set only for marker-gated rungs.
    """

    name: str
    bins: tuple[str, ...]
    markers: frozenset[str] | None
    rationale: str

    @property
    def mature(self) -> str:
        return self.bins[-1]


#: W1's proposal for the four frozen rung *names*. The names are frozen
#: (config/labeling_axes.yaml); these partitions are not.
RUNG_SPECS: Final[dict[str, RungSpec]] = {
    "epithelial": RungSpec(
        name="epithelial",
        bins=("epithelial",),
        markers=None,
        rationale=(
            "Coarsest rung: the whole epithelium is one population, so every "
            "epithelial cell counts as mature. Delta(mature fraction) then "
            "measures only epithelial-vs-non-epithelial shifts and everything "
            "within the epithelium lands in the intrinsic term. This is the "
            "lower bound of the granularity curve and it is supposed to look "
            "degenerate — that is what it demonstrates."
        ),
    ),
    "lineage": RungSpec(
        name="lineage",
        bins=("stem_like", "differentiated"),
        markers=None,
        rationale=(
            "Median split of the axis maturity score within each sample. Two "
            "bins, so roughly the absorptive/secretory/stem distinction the "
            "frozen rung description names, but derived from one axis's markers "
            "rather than from a joint clustering — using the other axis's genes "
            "here would collapse the two axes into one and destroy the "
            "agreement-across-axes argument."
        ),
    ),
    "crypt_position": RungSpec(
        name="crypt_position",
        bins=("crypt_bottom", "crypt_middle", "crypt_top"),
        markers=None,
        rationale=(
            "Tertiles of the axis maturity score within each sample, as a "
            "transcriptional proxy for crypt-bottom-to-top ordering. A proxy, "
            "not a measurement: real crypt position needs the spatial axis "
            "(axis 3, week 13+), and agreement between the two is what would "
            "license calling this positional."
        ),
    ),
    "best4": RungSpec(
        name="best4",
        bins=("other_epithelial", "best4"),
        markers=BEST4_MARKERS,
        rationale=(
            "Marker-gated, not score-binned: BEST4+ is a discrete population, "
            "not a point on a gradient. Top 5% of the BEST4-marker score within "
            "each sample, matching tier A's note that the program is under 5% of "
            "epithelium even in healthy colon. Finest rung, so the mature "
            "population is smallest and most of the change should appear as "
            "compositional — the upper bound of the granularity curve."
        ),
    ),
}


class LabelError(ValueError):
    """Labels could not be constructed from the input given."""


def label_column(axis: str, rung: str) -> str:
    """The column name for one (axis, rung) pair. `label_{axis}_{rung}`."""
    return f"label_{axis}_{rung}"


def label_columns(
    axes: Any = TRANSCRIPT_AXES, rungs: Any = None
) -> list[str]:
    """Every label column name, in a stable order."""
    rungs = list(rungs) if rungs is not None else granularity_rungs()
    return [label_column(a, r) for a in axes for r in rungs]


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def _positions(gene_names: Any, wanted: Any) -> tuple[np.ndarray, list[str]]:
    names = [str(n) for n in gene_names]
    lookup: dict[str, int] = {}
    for index, name in enumerate(names):
        lookup.setdefault(name, index)
    found = [g for g in sorted(wanted) if g in lookup]
    return np.array([lookup[g] for g in found], dtype=int), found


def _thin_to_depth(
    subset: np.ndarray, totals: np.ndarray, target: float, seed: int
) -> np.ndarray:
    """Binomial thinning of marker counts to a common sequencing depth.

    Each count is kept with probability ``target / total``, so every cell's
    marker counts are drawn as if it had been sequenced to `target` UMIs. The
    denominator is then the same for all cells and depth stops being a variable.

    This is the fix for the confound the pilot exposed: axis 1's mature cells
    came back four times shallower than its non-mature cells (median 4,791 counts
    against 18,829), because zero counts stay zero after CP10K normalisation, and
    because one stochastic count in a 1,000-UMI cell outranks ten in a
    20,000-UMI cell. Normalising a ratio cannot undo that; matching the depth can.
    """
    probability = np.clip(target / np.maximum(totals, 1.0), 0.0, 1.0)
    rng = np.random.default_rng(seed)
    return rng.binomial(subset.astype(np.int64), probability[:, None]).astype(float)


def score_markers(
    expression: Any,
    gene_names: Any,
    markers: Any,
    *,
    context: str,
    target_genes: Any,
    normalise: bool = True,
    depth_target: float | None = None,
    seed: int = 20260101,
) -> np.ndarray:
    """Per-cell score: mean z-scored, depth-normalised log expression of `markers`.

    Depth normalisation matters here. Chemistry is mixed across GSE178341's
    samples (v2 and v3 have different capture efficiency), so a raw-count score
    would rank a v3 cell above a v2 cell for technical reasons alone.

    `target_genes` is **required and has no default**, mirroring
    `build_signature()`'s refusal of an empty target set. Which genes count as
    targets for a given run is open decision #1, and a default would bury that
    decision: the whole panel makes axis 2 unusable, while a permissive default
    would silently disable invariant 2. The caller states it.
    """
    if target_genes is None or not list(target_genes):
        raise LabelError(
            "target_genes is empty. Pass the genes under test for this run — an "
            "empty target set silently disables CLAUDE.md invariant 2. Which "
            "genes to pass is open decision #1: the narrow reading is the target "
            "set for THIS run, not the whole panel."
        )
    assert_no_target_leakage(markers, target_genes, context=context)

    positions, found = _positions(gene_names, markers)
    if positions.size == 0:
        raise LabelError(
            f"{context}: none of {sorted(markers)} are in the matrix. Check gene "
            f"naming — this deposit carries symbols in var['gene_symbol'] and "
            f"Ensembl IDs in the index (open decision #3)."
        )
    missing = sorted(set(map(str, markers)) - set(found))

    subset = expression[:, positions]
    subset = np.asarray(subset.todense() if hasattr(subset, "todense") else subset, dtype=float)
    totals = np.asarray(expression.sum(axis=1), dtype=float).ravel()

    if depth_target is not None:
        # Depth-matched: thin to a common depth, then use that fixed denominator.
        subset = _thin_to_depth(subset, totals, depth_target, seed)
        subset = np.log1p(subset / depth_target * 1e4)
    elif normalise:
        safe = np.where(totals == 0, 1.0, totals).reshape(-1, 1)
        subset = np.log1p(subset / safe * 1e4)

    centre = subset.mean(axis=0)
    spread = subset.std(axis=0)
    spread[spread == 0] = 1.0
    scores = ((subset - centre) / spread).mean(axis=1)

    if missing:
        # Loud but not fatal: a partially present marker set still scores, and
        # which genes were dropped belongs in the run log.
        print(f"note: {context} missing {missing}; scored on {found}")
    return scores


def maturity_score(
    expression: Any, gene_names: Any, axis: str, *, target_genes: Any,
    normalise: bool = True, depth_target: float | None = None, seed: int = 20260101,
) -> np.ndarray:
    """Per-cell maturity along `axis`. Higher is more mature.

    Both transcript axes score *immaturity* — proximity to the stem pole, or
    commitment to the secretory lineage — so the score is negated. Markers come
    from the frozen `config/labeling_axes.yaml`, never from a local list.
    """
    if axis not in TRANSCRIPT_AXES:
        raise LabelError(
            f"{axis!r} is not a transcript-based axis; known: {TRANSCRIPT_AXES}. "
            f"chromatin and spatial are axis 3 and carry no markers (week 13+)."
        )
    markers = axis_genes(axis)
    if not markers:
        raise LabelError(f"axis {axis!r} has no markers in config/labeling_axes.yaml")
    return -score_markers(
        expression, gene_names, markers, context=f"axis {axis!r} labels",
        target_genes=target_genes, normalise=normalise,
        depth_target=depth_target, seed=seed,
    )


# ---------------------------------------------------------------------------
# Rung assignment
# ---------------------------------------------------------------------------


def _cut_points(
    values: np.ndarray, reference: np.ndarray, quantiles: np.ndarray
) -> np.ndarray | None:
    """Cut points taken from the reference cells only. None if unusable."""
    subset = values[reference]
    if subset.size == 0 or np.allclose(subset, subset[0]):
        return None
    return np.quantile(subset, quantiles)


def _bin_against_reference(
    values: np.ndarray,
    groups: np.ndarray,
    bins: tuple[str, ...],
    epithelial: np.ndarray,
    reference: np.ndarray,
) -> np.ndarray:
    """Bin every epithelial cell against cut points derived from `reference`.

    **This is the correction that makes the compositional term measurable.**

    An earlier version computed quantiles within each sample and binned that
    sample against its own cuts. That forces the mature fraction to equal the
    quantile in every sample by construction — a within-sample quantile cannot
    express a between-sample difference — so Delta(mature fraction) was pinned at
    zero and the compositional term, which is the whole project, could not move.
    Observed on the pilot: every `opposite_lineage` arm returned exactly 0.500 at
    the lineage rung and 0.333 at crypt_position.

    Cuts now come from the reference population — the patient's own **normal**
    tissue — and the same absolute threshold is applied to their tumour. Normal
    defines what mature looks like; tumour is then free to differ, which is
    precisely the paired contrast the decomposition is built on.

    `groups` is the unit the reference is taken within, normally the patient. A
    group with no usable reference gets the least-mature bin rather than an
    invented gradient — visible as a degenerate label, not as a plausible number.
    """
    out = np.full(values.shape, NON_EPITHELIAL, dtype=object)
    if len(bins) == 1:
        out[epithelial] = bins[0]
        return out

    quantiles = np.linspace(0, 1, len(bins) + 1)[1:-1]
    for group in pd.unique(groups):
        here = epithelial & (groups == group)
        if not here.any():
            continue
        cuts = _cut_points(values, here & reference, quantiles)
        if cuts is None:
            out[here] = bins[0]
            continue

        # Coincident cut points collapse a bin. It happens when a large block of
        # cells share a score — on axis 1, every cell with zero counts across all
        # five stem markers — so a quantile boundary lands inside the tie and the
        # next one lands on the same value. searchsorted can then never return
        # the middle index and that bin silently never appears: the pilot's
        # stem_pole crypt_position came back as two bins, with crypt_middle
        # absent and the result identical to the lineage rung.
        #
        # Use the distinct boundaries only and keep the extremes. A binary split
        # honestly reported beats a three-way split with an empty middle.
        distinct = np.unique(cuts)
        index = np.searchsorted(distinct, values[here], side="right")
        if distinct.size == len(bins) - 1:
            usable = np.asarray(bins, dtype=object)
        else:
            usable = np.asarray((bins[0], *bins[-distinct.size:]), dtype=object)
            print(
                f"note: group {group!r} supports only {distinct.size + 1} of "
                f"{len(bins)} bins — scores are tied across a quantile boundary. "
                f"Using {tuple(usable)}."
            )
        out[here] = usable[index]
    return out


def _gate_against_reference(
    values: np.ndarray,
    groups: np.ndarray,
    bins: tuple[str, ...],
    epithelial: np.ndarray,
    reference: np.ndarray,
    quantile: float,
) -> np.ndarray:
    """Marker gate, thresholded on the reference cells only.

    Same correction as :func:`_bin_against_reference`: a per-sample top-5% gate
    returns 5% in every sample whatever the biology, which is why the pilot's
    `best4` fraction was 0.050 in all ten arms.
    """
    out = np.full(values.shape, NON_EPITHELIAL, dtype=object)
    for group in pd.unique(groups):
        here = epithelial & (groups == group)
        if not here.any():
            continue
        subset = values[here & reference]
        if subset.size == 0:
            out[here] = bins[0]
            continue
        threshold = np.quantile(subset, quantile)
        out[here] = np.where(values[here] >= threshold, bins[-1], bins[0])
    return out


def assign_labels(
    expression: Any,
    gene_names: Any,
    *,
    compartment: Any,
    sample_id: Any,
    target_genes: Any,
    tissue: Any = None,
    patient_id: Any = None,
    reference_tissue: str = "normal",
    depth_target: float | None = None,
    depth_quantile: float = 0.10,
    seed: int = 20260101,
    axes: Any = TRANSCRIPT_AXES,
    rungs: Any = None,
    normalise: bool = True,
    index: Any = None,
) -> pd.DataFrame:
    """Build every `label_{axis}_{rung}` column. One row per cell, in input order.

    Parameters
    ----------
    compartment:
        Per-cell compartment from `ingest.assign_compartments()`. Only
        ``"epithelial"`` cells are scored; everything else is labelled
        ``non_epithelial`` in every column. Scoring a fibroblast for
        differentiation would put nonsense into the denominator.
    sample_id:
        Retained for provenance and for grouping when `patient_id` is absent.
    tissue, patient_id, reference_tissue:
        **How the cut points are set, and the reason the compositional term is
        measurable at all.** Thresholds are taken from the cells where
        ``tissue == reference_tissue`` — the patient's own normal — and the same
        absolute cut is applied to their tumour. Grouping is by `patient_id` when
        given, else by `sample_id`.

        Do not go back to per-sample quantiles. Binning each sample against its
        own quantiles pins the mature fraction to that quantile in every sample,
        so Delta(mature fraction) is identically zero by construction. The pilot
        showed it plainly: every `opposite_lineage` arm returned 0.500 at the
        lineage rung and 0.333 at crypt_position, and every `best4` arm 0.050.

        With `tissue=None` the reference falls back to all epithelial cells,
        which is only correct when there is no tumour/normal contrast to make.
    depth_target, depth_quantile:
        **Depth matching.** Pass a target depth, or leave it None to use the
        `depth_quantile` of epithelial totals. Marker counts are binomially
        thinned to that depth so every cell is scored as if sequenced equally.

        Without it the maturity call partly measures sequencing depth: on the
        pilot, axis 1's mature cells were four times shallower than its
        non-mature cells. Cells whose own depth is below the target cannot be
        thinned up and are labelled `unresolved_depth` — measurable in principle,
        not measurable here.

        Set `depth_target=0` to disable matching entirely, which reproduces the
        confounded behaviour and should only be used to demonstrate it.
    target_genes:
        The genes under test for this run. Required, no default (invariant 2).
        **Labels are therefore specific to a target set**, which is the price of
        open decision #1's narrow reading: MUC2 and TFF3 are simultaneously tier-E
        targets and axis-2 markers, so a run testing either cannot use axis 2 —
        pass ``axes=["stem_pole"]`` for those, and expect a LeakageError if you
        forget.

    Every column coexists; nothing is overwritten. That is what makes the
    granularity-and-axis curve computable later (execution_plan.md §6.2).
    """
    rungs = list(rungs) if rungs is not None else granularity_rungs()
    unknown = [r for r in rungs if r not in RUNG_SPECS]
    if unknown:
        raise LabelError(f"no RungSpec for {unknown}; known: {sorted(RUNG_SPECS)}")

    n_cells = expression.shape[0]
    compartment = np.asarray([str(c) for c in compartment], dtype=object)
    groups = np.asarray([str(s) for s in sample_id], dtype=object)
    for name, arr in (("compartment", compartment), ("sample_id", groups)):
        if arr.shape[0] != n_cells:
            raise LabelError(f"{name} has {arr.shape[0]} entries for {n_cells} cells")

    if patient_id is not None:
        groups = np.asarray([str(p) for p in patient_id], dtype=object)
        if groups.shape[0] != n_cells:
            raise LabelError(f"patient_id has {groups.shape[0]} entries for {n_cells} cells")

    epithelial = compartment == "epithelial"
    if tissue is None:
        reference = epithelial.copy()
    else:
        tissue_arr = np.asarray([str(t) for t in tissue], dtype=object)
        if tissue_arr.shape[0] != n_cells:
            raise LabelError(f"tissue has {tissue_arr.shape[0]} entries for {n_cells} cells")
        reference = epithelial & (tissue_arr == reference_tissue)
        if not reference.any():
            raise LabelError(
                f"no epithelial cells with tissue == {reference_tissue!r}. Cut "
                f"points come from the reference arm; without it the mature "
                f"fraction cannot be compared between arms."
            )

    if not epithelial.any():
        raise LabelError(
            "no cells are labelled 'epithelial'. Pass compartments from "
            "ingest.assign_compartments(), whose epithelial value is 'epithelial'."
        )

    totals = np.asarray(expression.sum(axis=1), dtype=float).ravel()
    if depth_target is None:
        depth_target = float(np.quantile(totals[epithelial], depth_quantile))
    matched = depth_target > 0
    if matched:
        resolvable = epithelial & (totals >= depth_target)
        if not resolvable.any():
            raise LabelError(
                f"no epithelial cell reaches the depth target {depth_target:.0f}. "
                f"Lower depth_quantile, or pass depth_target explicitly."
            )
        dropped = int(epithelial.sum() - resolvable.sum())
        if dropped:
            print(
                f"note: {dropped:,} of {int(epithelial.sum()):,} epithelial cells "
                f"are below the depth target {depth_target:,.0f} and are labelled "
                f"{UNRESOLVED!r} — not scored, not counted as immature."
            )
        # The reference must also be depth-matched, or the cut points come from a
        # differently-sequenced population than the cells being cut.
        reference = reference & resolvable
        if not reference.any():
            raise LabelError(
                "no reference cell survives depth matching; the cut points would "
                "come from nothing. Lower depth_quantile."
            )
    else:
        resolvable = epithelial

    frame = pd.DataFrame(index=pd.RangeIndex(n_cells) if index is None else pd.Index(index))
    best4_score: np.ndarray | None = None

    for axis in axes:
        scores = maturity_score(
            expression, gene_names, axis, target_genes=target_genes,
            normalise=normalise,
            depth_target=depth_target if matched else None, seed=seed,
        )
        for rung in rungs:
            spec = RUNG_SPECS[rung]
            if spec.markers is None:
                labels = _bin_against_reference(
                    scores, groups, spec.bins, resolvable, reference
                )
            else:
                if best4_score is None:
                    best4_score = score_markers(
                        expression, gene_names, spec.markers,
                        context=f"rung {rung!r} labels",
                        target_genes=target_genes, normalise=normalise,
                        depth_target=depth_target if matched else None, seed=seed,
                    )
                labels = _gate_against_reference(
                    best4_score, groups, spec.bins, resolvable, reference, BEST4_QUANTILE
                )
            labels[epithelial & ~resolvable] = UNRESOLVED
            frame[label_column(axis, rung)] = pd.Categorical(
                labels, categories=[NON_EPITHELIAL, UNRESOLVED, *spec.bins]
            )
    return frame


# ---------------------------------------------------------------------------
# What downstream consumers need
# ---------------------------------------------------------------------------


def mature_mask(labels: pd.DataFrame, axis: str, rung: str) -> np.ndarray:
    """Boolean: which cells count as mature at this (axis, rung).

    The mature category is `RUNG_SPECS[rung].mature`, which is the most-mature
    bin. This is the mask the compositional and intrinsic terms are both built
    on, so it is a single function rather than a convention each caller repeats.
    """
    column = label_column(axis, rung)
    if column not in labels.columns:
        raise LabelError(f"{column} not in labels; have {list(labels.columns)}")
    if rung not in RUNG_SPECS:
        raise LabelError(f"no RungSpec for {rung!r}")
    return (labels[column].astype(str) == RUNG_SPECS[rung].mature).to_numpy()


def mature_cell_counts(
    labels: pd.DataFrame,
    *,
    patient_id: Any,
    tissue: Any,
    axes: Any = TRANSCRIPT_AXES,
    rungs: Any = None,
) -> pd.DataFrame:
    """Mature-cell counts per (patient, tissue, axis, rung). Long form.

    This is `n_cells_mature` in the frozen output schema, and it is what decides
    positivity — whether a patient's intrinsic term is estimable at all. The
    thresholds themselves belong to W2 (`src/harness/positivity.py`); this
    supplies the counts, deliberately without applying them.
    """
    rungs = list(rungs) if rungs is not None else granularity_rungs()
    keys = pd.DataFrame(
        {
            "patient_id": [str(p) for p in patient_id],
            "tissue": [str(t) for t in tissue],
        }
    )
    if len(keys) != len(labels):
        raise LabelError(f"patient_id has {len(keys)} entries for {len(labels)} cells")

    rows: list[pd.DataFrame] = []
    for axis in axes:
        for rung in rungs:
            mature = mature_mask(labels, axis, rung)
            # .to_numpy() on BOTH: `keys` has a RangeIndex while `labels` is
            # indexed by barcode, so assigning a Series here aligns on index,
            # matches nothing, and silently yields an all-NaN column.
            column = labels[label_column(axis, rung)].astype(str).to_numpy()
            epithelial = column != NON_EPITHELIAL
            unresolved = column == UNRESOLVED
            grouped = (
                keys.assign(mature=mature, epithelial=epithelial, unresolved=unresolved)
                .groupby(["patient_id", "tissue"], observed=True)
                .agg(
                    n_cells_mature=("mature", "sum"),
                    n_cells_epithelial=("epithelial", "sum"),
                    n_cells_unresolved=("unresolved", "sum"),
                )
                .reset_index()
            )
            grouped["labeling_axis"] = axis
            grouped["granularity_rung"] = rung
            # Denominator is the RESOLVED epithelium. Cells dropped by depth
            # matching are reported separately rather than counted as immature —
            # a cell that could not be measured is not a cell measured to be
            # immature (open decision #14).
            grouped["n_cells_resolved"] = (
                grouped["n_cells_epithelial"] - grouped["n_cells_unresolved"]
            )
            # np.where evaluates BOTH branches, so a guard around the division
            # does not prevent it. Blank the denominator instead.
            denominator = grouped["n_cells_resolved"].astype(float)
            grouped["mature_fraction"] = (
                grouped["n_cells_mature"].astype(float)
                / denominator.where(denominator > 0)
            )
            # How much of the epithelium the fraction could not speak for. A
            # large value means the fraction is bounded, not measured.
            epithelium = grouped["n_cells_epithelial"].astype(float)
            grouped["unresolved_fraction"] = (
                grouped["n_cells_unresolved"].astype(float)
                / epithelium.where(epithelium > 0)
            )
            rows.append(grouped)
    return pd.concat(rows, ignore_index=True)


def describe_labels(labels: pd.DataFrame) -> pd.DataFrame:
    """Cell counts per label value, one row per (column, value). For run logs."""
    rows = []
    for column in labels.columns:
        counts = labels[column].astype(str).value_counts()
        for value, count in counts.items():
            rows.append({"column": column, "label": value, "n_cells": int(count)})
    return pd.DataFrame(rows)


def cell_type_vector(
    labels: pd.DataFrame,
    axis: str,
    rung: str,
    *,
    mature_label: str = "mature_colonocyte",
) -> np.ndarray:
    """A ``cell_type`` array in the form W2's pseudobulk generator consumes.

    ``src.harness.pseudobulk.generate_pseudobulk`` takes ``cell_type: Sequence[str]``
    and identifies the mature population by string equality against its
    ``mature_label`` argument, which defaults to ``"mature_colonocyte"``. The bin
    names in :data:`RUNG_SPECS` are rung-specific (``differentiated``,
    ``crypt_top``, ``best4``) and deliberately so, so this renames whichever bin
    is mature at this rung to `mature_label` and leaves every other value alone.

    Use it rather than hand-mapping at the call site: which bin counts as mature
    is a modelling choice recorded in `RUNG_SPECS`, and duplicating that mapping
    is how the two drift apart.
    """
    column = label_column(axis, rung)
    if column not in labels.columns:
        raise LabelError(f"{column} not in labels; have {list(labels.columns)}")
    values = labels[column].astype(str).to_numpy().copy()
    values[values == RUNG_SPECS[rung].mature] = mature_label
    return values


def maturity_summary(
    labels: pd.DataFrame,
    *,
    patient_id: Any,
    tissue: Any,
    study_id: str,
    axes: Any = TRANSCRIPT_AXES,
    rungs: Any = None,
) -> pd.DataFrame:
    """Per-patient mature fractions, shaped for W4's estimator.

    ``src.estimator.kitagawa.decompose_cohort`` requires ``patient_id``,
    ``study_id``, ``gene``, ``granularity_rung``, ``labeling_axis``,
    ``frac_mature_normal``, ``frac_mature_tumour``, ``mean_normal``,
    ``mean_tumour`` and ``n_cells_mature``. This supplies everything that depends
    only on labels; the caller adds ``gene``, ``mean_normal`` and ``mean_tumour``,
    which need the target gene's expression.

    ``n_cells_mature`` is the **tumour** mature count, because the intrinsic term
    is `tumour mature fraction x Delta(per-cell mean)` — positivity is about how
    many mature cells survive in the tumour, which is the arm that can run out.
    Patients missing either arm are dropped: without both, neither term is
    identifiable (open decision #9).
    """
    counts = mature_cell_counts(
        labels, patient_id=patient_id, tissue=tissue, axes=axes, rungs=rungs
    )
    wide = counts.pivot_table(
        index=["patient_id", "labeling_axis", "granularity_rung"],
        columns="tissue",
        values=["n_cells_mature", "n_cells_epithelial", "mature_fraction"],
        observed=True,
    )
    wide.columns = [f"{value}_{tissue}" for value, tissue in wide.columns]
    wide = wide.reset_index()

    needed = ["mature_fraction_normal", "mature_fraction_tumour",
              "n_cells_mature_tumour"]
    missing = [c for c in needed if c not in wide.columns]
    if missing:
        raise LabelError(
            f"cannot build a summary: {missing} absent. Both tumour and normal "
            f"arms are required — a patient with only one contributes to neither "
            f"term (open decision #9)."
        )

    out = wide.rename(
        columns={
            "mature_fraction_normal": "frac_mature_normal",
            "mature_fraction_tumour": "frac_mature_tumour",
        }
    )
    out["n_cells_mature"] = out["n_cells_mature_tumour"].astype("Int64")
    out["study_id"] = study_id
    out = out.dropna(subset=["frac_mature_normal", "frac_mature_tumour"])
    return out.reset_index(drop=True)


def axis_tie_fraction(
    expression: Any,
    gene_names: Any,
    axis: str,
    *,
    target_genes: Any,
    epithelial: Any = None,
    normalise: bool = True,
    depth_target: float | None = None,
    seed: int = 20260101,
) -> dict[str, float]:
    """How much of the maturity score is a single tied block. **Run this first.**

    Axis 1's markers (LGR5, ASCL2, MKI67, OLFM4, SMOC2) are sparsely detected, so
    a large share of epithelial cells carry zero counts for all five and share an
    identical score. Quantile boundaries cannot split a tie, which is what
    collapsed `crypt_position` into two bins on the pilot.

    The deeper issue is interpretive, not mechanical: a cell with no detected
    stem markers might be genuinely differentiated, or might simply be shallow.
    Scoring it as maximally mature is inference from absence of evidence, and how
    to treat those cells is a decision for the team rather than a default — see
    the note in the module docstring.

    Returns the tied fraction, the size of the largest tied block, and the number
    of distinct score values, so the problem is a number rather than an
    impression.
    """
    scores = maturity_score(
        expression, gene_names, axis, target_genes=target_genes,
        normalise=normalise, depth_target=depth_target, seed=seed,
    )
    keep = (
        np.ones(scores.shape, dtype=bool)
        if epithelial is None
        else np.asarray(epithelial, dtype=bool)
    )
    if depth_target:
        # Cells below the target cannot be thinned up — _thin_to_depth clips the
        # probability at 1 and leaves them alone — so including them here mixes
        # thinned and unthinned cells and reports a tie fraction that no actual
        # run would produce. assign_labels drops them as UNRESOLVED; match it.
        totals = np.asarray(expression.sum(axis=1), dtype=float).ravel()
        keep = keep & (totals >= depth_target)
        if not keep.any():
            raise LabelError(
                f"no cell reaches the depth target {depth_target:,.0f}"
            )
    scores = scores[keep]
    if scores.size == 0:
        raise LabelError("no cells to score")

    values, counts = np.unique(np.round(scores, 10), return_counts=True)
    largest = int(counts.max())
    return {
        "n_cells": int(scores.size),
        "n_distinct_scores": int(values.size),
        "largest_tied_block": largest,
        "tied_fraction": float(largest / scores.size),
        "resolvable_fraction": float(1.0 - largest / scores.size),
    }


def label_depth_confounding(
    labels: pd.DataFrame,
    metrics: pd.DataFrame,
    *,
    axes: Any = TRANSCRIPT_AXES,
    rungs: Any = None,
    warn_ratio: float = 1.25,
    warn_auc: float = 0.10,
) -> pd.DataFrame:
    """Is the maturity call tracking sequencing depth rather than biology?

    **Run this before quoting any compositional number.**

    Both axes score on sparsely detected markers, and zero counts stay zero after
    depth normalisation — a shallow cell is more likely to have none of the five
    stem markers and so be called mature. On the pilot, axis 1's mature bin was
    *exactly* the tied block of cells with no stem-marker counts at all (7,593 of
    16,955), which makes the concern concrete rather than theoretical.

    If mature cells are systematically shallower than non-mature ones, the mature
    fraction is partly a depth measurement, and Delta(mature fraction) between
    two arms with different depth is partly an artifact. `ratio` below 1 means
    mature cells are shallower; `flagged` marks a gap beyond `warn_ratio` in
    either direction.

    Two statistics, because neither is sufficient alone. `counts_ratio` compares
    medians and is what caught the pilot's 4x gap, but it is fragile when depth is
    bimodal — a bin split near 50/50 flips its median to whichever side holds one
    extra cell. `depth_auc` is the rank probability that a mature cell is deeper
    than a non-mature one: 0.5 means no association, and it is unaffected by
    bimodality. Either exceeding its tolerance sets `flagged`.

    Returns one row per (axis, rung).
    """
    rungs = list(rungs) if rungs is not None else granularity_rungs()
    for column in ("n_counts", "n_genes"):
        if column not in metrics.columns:
            raise LabelError(f"metrics needs a {column} column (from cell_qc_metrics)")
    if len(metrics) != len(labels):
        raise LabelError(f"metrics has {len(metrics)} rows for {len(labels)} cells")

    counts = np.asarray(metrics["n_counts"], dtype=float)
    genes = np.asarray(metrics["n_genes"], dtype=float)

    def _auc(values: np.ndarray, positive: np.ndarray) -> float:
        """P(a mature cell is deeper than a non-mature one). 0.5 = no association."""
        order = np.argsort(values, kind="mergesort")
        ranks = np.empty(values.size, dtype=float)
        ranks[order] = np.arange(1, values.size + 1)
        # Average ranks within ties so a tied block cannot fake an association.
        unique, inverse, counts_ = np.unique(values, return_inverse=True, return_counts=True)
        sums = np.zeros(unique.size)
        np.add.at(sums, inverse, ranks)
        ranks = (sums / counts_)[inverse]
        n1, n0 = int(positive.sum()), int((~positive).sum())
        if n1 == 0 or n0 == 0:
            return float("nan")
        return float((ranks[positive].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

    rows = []
    for axis in axes:
        for rung in rungs:
            column = labels[label_column(axis, rung)].astype(str).to_numpy()
            # UNRESOLVED cells are excluded from BOTH sides. They are defined by
            # having low depth, so leaving them in the comparison guarantees a
            # perfect association and reports AUC 1.0 for the epithelial rung —
            # a property of the diagnostic, not of the labels.
            scored = ~np.isin(column, [NON_EPITHELIAL, UNRESOLVED])
            mature = column == RUNG_SPECS[rung].mature
            other = scored & ~mature
            if not mature.any() or not other.any():
                continue
            median_mature = float(np.median(counts[mature]))
            median_other = float(np.median(counts[other]))
            ratio = median_mature / median_other if median_other else np.nan
            area = _auc(counts[scored], mature[scored])
            rows.append(
                {
                    "labeling_axis": axis,
                    "granularity_rung": rung,
                    "n_mature": int(mature.sum()),
                    "median_counts_mature": median_mature,
                    "median_counts_other": median_other,
                    "counts_ratio": ratio,
                    "median_genes_mature": float(np.median(genes[mature])),
                    "median_genes_other": float(np.median(genes[other])),
                    "depth_auc": area,
                    "flagged": bool(
                        (np.isfinite(ratio)
                         and (ratio < 1 / warn_ratio or ratio > warn_ratio))
                        or (np.isfinite(area) and abs(area - 0.5) > warn_auc)
                    ),
                }
            )
    return pd.DataFrame(rows)


def maturity_within_depth_strata(
    labels: pd.DataFrame,
    metrics: pd.DataFrame,
    *,
    axis: str = "stem_pole",
    rung: str = "lineage",
    n_strata: int = 10,
) -> pd.DataFrame:
    """The maturity call's depth gradient, stratum by stratum. Descriptive.

    Reports the mature fraction inside each depth decile. A steep monotone slope
    — near 1 in the shallowest stratum, near 0 in the deepest — means the call
    is largely tracking depth; a flat profile means it is not.

    **It does not separate technical from biological confounding**, and an
    earlier version of this docstring claimed it did. Dropout is stochastic: at
    a fixed depth some cells lose all five sparse markers by chance and some do
    not, so even a purely technical call keeps a mix inside every stratum. That
    mix is noise, and this function cannot tell it from signal.

    For that question use :func:`annotation_concordance` against an independent
    annotation. A call made of dropout noise will not agree with one; a call
    measuring real maturity will.
    """
    if "n_counts" not in metrics.columns:
        raise LabelError("metrics needs an n_counts column (from cell_qc_metrics)")
    if len(metrics) != len(labels):
        raise LabelError(f"metrics has {len(metrics)} rows for {len(labels)} cells")

    column = labels[label_column(axis, rung)].astype(str).to_numpy()
    scored = ~np.isin(column, [NON_EPITHELIAL, UNRESOLVED])
    if not scored.any():
        raise LabelError("no scored cells")

    counts = np.asarray(metrics["n_counts"], dtype=float)[scored]
    mature = (column == RUNG_SPECS[rung].mature)[scored]

    edges = np.quantile(counts, np.linspace(0, 1, n_strata + 1))
    edges = np.unique(edges)
    stratum = np.clip(np.searchsorted(edges, counts, side="right") - 1, 0, len(edges) - 2)

    rows = []
    for index in range(len(edges) - 1):
        here = stratum == index
        if not here.any():
            continue
        fraction = float(mature[here].mean())
        rows.append(
            {
                "stratum": index,
                "n_cells": int(here.sum()),
                "counts_low": float(edges[index]),
                "counts_high": float(edges[index + 1]),
                "mature_fraction": fraction,
            }
        )
    return pd.DataFrame(rows)


def annotation_concordance(
    labels: pd.DataFrame,
    annotation: Any,
    *,
    axis: str = "stem_pole",
    rung: str = "lineage",
    immature_pattern: str = "stem|TA-like|prolif",
) -> dict[str, Any]:
    """Does the maturity call agree with an INDEPENDENT annotation?

    **The test that separates signal from dropout noise.**

    :func:`label_depth_confounding` shows the call is associated with depth, and
    :func:`maturity_within_depth_strata` shows how steeply — but neither can say
    whether what remains is biology or noise, because dropout is stochastic and
    produces a mix at every depth just as real variation would.

    Agreement with an independent annotation can. A call made of dropout noise
    has nothing to agree with; one measuring real maturity will track a
    well-informed annotation built from many genes and a clustering.

    For GSE178341 that annotation is the authors' ``cl295v11SubFull``, whose
    epithelial subsets are named — ``cE01 (Stem/TA-like)``,
    ``cE03 (Stem/TA-like prolif)``. Cells whose annotation matches
    `immature_pattern` are treated as immature.

    **This is validation, not labelling.** Their clustering is transcriptional
    and may have used panel genes, so it must never become a label
    (CLAUDE.md invariant 2). Using it to ask whether our own independent call
    recovers the same structure is a different act, and the direction matters:
    agreement is evidence our call is real; disagreement is evidence it is not.

    Returns the 2x2 counts, agreement, sensitivity, specificity, and Cohen's
    kappa — kappa because raw agreement is inflated when one class dominates.
    """
    column = label_column(axis, rung)
    if column not in labels.columns:
        raise LabelError(f"{column} not in labels")
    values = labels[column].astype(str).to_numpy()
    scored = ~np.isin(values, [NON_EPITHELIAL, UNRESOLVED])
    if not scored.any():
        raise LabelError("no scored cells to compare")

    reference = pd.Series(annotation).astype(str).to_numpy()
    if reference.shape[0] != len(labels):
        raise LabelError(
            f"annotation has {reference.shape[0]} entries for {len(labels)} cells"
        )

    ours_mature = (values == RUNG_SPECS[rung].mature)[scored]
    theirs_immature = pd.Series(reference[scored]).str.contains(
        immature_pattern, case=False, regex=True, na=False
    ).to_numpy()
    theirs_mature = ~theirs_immature

    tp = int((ours_mature & theirs_mature).sum())
    tn = int((~ours_mature & theirs_immature).sum())
    fp = int((ours_mature & theirs_immature).sum())
    fn = int((~ours_mature & theirs_mature).sum())
    total = tp + tn + fp + fn
    if total == 0:
        raise LabelError("no overlapping cells")

    agreement = (tp + tn) / total
    expected = (
        ((tp + fp) * (tp + fn) + (tn + fn) * (tn + fp)) / (total * total)
    )
    kappa = (agreement - expected) / (1 - expected) if expected < 1 else float("nan")
    return {
        "n_cells": total,
        "n_mature_ours": tp + fp,
        "n_mature_theirs": tp + fn,
        "agreement": float(agreement),
        "sensitivity": float(tp / (tp + fn)) if (tp + fn) else float("nan"),
        "specificity": float(tn / (tn + fp)) if (tn + fp) else float("nan"),
        "kappa": float(kappa),
        # Kappa near 0 means the call carries no more information than chance
        # about an independently-derived maturity annotation.
        "informative": bool(np.isfinite(kappa) and kappa > 0.2),
    }
