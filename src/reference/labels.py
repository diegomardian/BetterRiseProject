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


def score_markers(
    expression: Any,
    gene_names: Any,
    markers: Any,
    *,
    context: str,
    target_genes: Any,
    normalise: bool = True,
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

    if normalise:
        totals = np.asarray(expression.sum(axis=1), dtype=float).reshape(-1, 1)
        totals[totals == 0] = 1.0
        subset = np.log1p(subset / totals * 1e4)

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
    normalise: bool = True,
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
    )


# ---------------------------------------------------------------------------
# Rung assignment
# ---------------------------------------------------------------------------


def _bin_within_groups(
    values: np.ndarray, groups: np.ndarray, bins: tuple[str, ...], epithelial: np.ndarray
) -> np.ndarray:
    """Quantile-bin `values` into `bins`, computed separately within each group.

    Per sample, not pooled: depth and composition differ between samples, so a
    pooled quantile would let one sample's library size decide another's labels.
    """
    out = np.full(values.shape, NON_EPITHELIAL, dtype=object)
    if len(bins) == 1:
        out[epithelial] = bins[0]
        return out

    edges = np.linspace(0, 1, len(bins) + 1)[1:-1]
    for group in pd.unique(groups):
        here = epithelial & (groups == group)
        if not here.any():
            continue
        subset = values[here]
        if np.allclose(subset, subset[0]):
            # No variation to bin on. Assign the least-mature bin rather than
            # inventing a gradient.
            out[here] = bins[0]
            continue
        cuts = np.quantile(subset, edges)
        out[here] = np.asarray(bins, dtype=object)[np.searchsorted(cuts, subset, side="right")]
    return out


def _gate_within_groups(
    values: np.ndarray, groups: np.ndarray, bins: tuple[str, ...],
    epithelial: np.ndarray, quantile: float,
) -> np.ndarray:
    """Call the top `quantile` of `values` the mature bin, per group."""
    out = np.full(values.shape, NON_EPITHELIAL, dtype=object)
    for group in pd.unique(groups):
        here = epithelial & (groups == group)
        if not here.any():
            continue
        threshold = np.quantile(values[here], quantile)
        out[here] = np.where(values[here] >= threshold, bins[-1], bins[0])
    return out


def assign_labels(
    expression: Any,
    gene_names: Any,
    *,
    compartment: Any,
    sample_id: Any,
    target_genes: Any,
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
        Quantile bins and the BEST4 gate are computed **within** each sample.
        Pooling them would let one sample's depth decide another's labels, and
        chemistry differs across samples here.
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

    epithelial = compartment == "epithelial"
    if not epithelial.any():
        raise LabelError(
            "no cells are labelled 'epithelial'. Pass compartments from "
            "ingest.assign_compartments(), whose epithelial value is 'epithelial'."
        )

    frame = pd.DataFrame(index=pd.RangeIndex(n_cells) if index is None else pd.Index(index))
    best4_score: np.ndarray | None = None

    for axis in axes:
        scores = maturity_score(
            expression, gene_names, axis, target_genes=target_genes, normalise=normalise
        )
        for rung in rungs:
            spec = RUNG_SPECS[rung]
            if spec.markers is None:
                labels = _bin_within_groups(scores, groups, spec.bins, epithelial)
            else:
                if best4_score is None:
                    best4_score = score_markers(
                        expression, gene_names, spec.markers,
                        context=f"rung {rung!r} labels",
                        target_genes=target_genes, normalise=normalise,
                    )
                labels = _gate_within_groups(
                    best4_score, groups, spec.bins, epithelial, BEST4_QUANTILE
                )
            frame[label_column(axis, rung)] = pd.Categorical(
                labels, categories=[NON_EPITHELIAL, *spec.bins]
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
            grouped = (
                keys.assign(mature=mature, epithelial=labels[label_column(axis, rung)]
                            .astype(str)
                            .ne(NON_EPITHELIAL))
                .groupby(["patient_id", "tissue"], observed=True)
                .agg(n_cells_mature=("mature", "sum"), n_cells_epithelial=("epithelial", "sum"))
                .reset_index()
            )
            grouped["labeling_axis"] = axis
            grouped["granularity_rung"] = rung
            # np.where evaluates BOTH branches, so a guard around the division
            # does not prevent it. Blank the denominator instead: a group with no
            # epithelium has no mature fraction, and NaN says so.
            denominator = grouped["n_cells_epithelial"].astype(float)
            grouped["mature_fraction"] = (
                grouped["n_cells_mature"].astype(float)
                / denominator.where(denominator > 0)
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
