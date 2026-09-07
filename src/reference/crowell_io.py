"""Reading the Crowell WTx CosMx deposit, assuming as little as possible.

Zenodo ``10.5281/zenodo.15574384``, CC-BY-4.0. Whole-transcriptome CosMx SMI:
19,867 genes over >99.5% of annotated protein-coding, on FFPE sections carrying
reference mucosa, tubulovillous adenoma and carcinoma **within the same section
from the same patient**. Seven patients, eight colon sections plus a lymph node.

Pre-registered in ``docs/prereg_crowell_feasibility.md``, committed in
``cd7e4a5`` before anything was downloaded.

WHAT THIS MODULE REFUSES TO ASSUME, and why each one has a name.

*The domain vocabulary.* The pre-registration names REF / TVA / CRC because the
paper does, and the paper's words are not necessarily the object's words.
:func:`domain_vocabulary` reports every candidate column and its values and
**maps nothing**. `ADENOMA_TISSUE_MAP` was written against one atlas's words and
reading a label rather than the grouping once put Chen_2021's usable pairs at
**zero** when the true number was 44.

*Which section carries which domain.* ``232`` is read first because it is the
smallest file, not because anyone has seen a TVA on it. If it carries none, the
pre-registration names ``231`` as the fallback, so the choice is not made after
seeing a number.

*The negative-probe naming.* CosMx writes negative probes and false codes as
features alongside real genes, and the prefix convention varies by release.
:func:`control_features` looks for several and reports which matched, because
this is the per-cell false-positive floor the whole gate is measured against —
a floor built from an empty match would be zero, and a zero floor makes every
gene look infinitely separated.

*That the counts are counts.* :func:`check_counts_are_integers` refuses a matrix
that has already been normalised. A detection rate off normalised data is not a
detection rate, and nothing about that raises on its own.

WHAT IT DOES ASSUME, stated so it can be checked: that the file is an ``.h5ad``
readable by ``anndata``, and that one row is one cell. Both are checked on read
and raise if false.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)


class CrowellError(RuntimeError):
    """The deposit does not have the shape this module verified."""


#: The section named in the pre-registration, and its fallback. Sizes and md5s
#: are from the Zenodo API listing, recorded so a re-download is checkable.
SECTIONS: dict[str, dict[str, object]] = {
    "232": {"bytes": 226_800_000, "md5": "f2271485bae235439267c1bd13f1a7c9",
            "role": "primary — smallest of eight"},
    "231": {"bytes": 534_800_000, "md5": "fdd226ec7632a3d45a4f9dd46c24a981",
            "role": "fallback, named in prereg §2 before any read"},
}

#: obs columns that might carry the pathologist's domain. Candidates, not a
#: mapping: the job prints what it finds and a human chooses.
DOMAIN_COLUMN_CANDIDATES: tuple[str, ...] = (
    "domain", "domains", "region", "regions", "annotation", "anno",
    "histo", "histology", "pathology", "path_anno", "tissue", "tissue_type",
    "niche", "compartment", "label", "sample_type", "condition",
)

#: Words the pre-registration expects to find among those values. Used ONLY to
#: report whether they appear — never to filter, rename or map.
EXPECTED_DOMAIN_WORDS: tuple[str, ...] = (
    "ref", "tva", "crc", "adenoma", "carcinoma", "normal", "mucosa",
    "tumor", "tumour", "ln", "lymph",
)

#: Prefixes CosMx has used for control features across releases. Reported with
#: which one matched, because an empty match makes the floor zero.
NEGATIVE_PROBE_PATTERNS: tuple[str, ...] = (
    r"^Negative", r"^NegPrb", r"^NegControlProbe", r"^SystemControl",
)
FALSE_CODE_PATTERNS: tuple[str, ...] = (
    r"^FalseCode", r"^NegControlCodeword", r"^Falsecode", r"^UnassignedCodeword",
)


def open_section(path: str | Path, *, backed: bool = True):
    """Open one section's ``.h5ad``. Backed by default — these are 227 MB up."""
    # The path is checked FIRST. Importing anndata first made a missing file
    # report "anndata is required", which sends the reader to fix an
    # environment when what is actually wrong is that nothing was downloaded.
    path = Path(path)
    if not path.exists():
        raise CrowellError(
            f"{path} not found. The prereg names 232.h5ad (226.8 MB, md5 "
            f"f2271485bae235439267c1bd13f1a7c9) from Zenodo "
            f"10.5281/zenodo.15574384; record it in data/manifest.csv on "
            f"download."
        )
    try:
        import anndata
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise CrowellError(
            "anndata is required to read this deposit; env/w1_reference.yml "
            "pins it."
        ) from exc
    adata = anndata.read_h5ad(str(path), backed="r" if backed else None)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise CrowellError(f"{path} is {adata.shape}; nothing to read.")
    return adata


def _match_any(names: np.ndarray | list[str], patterns: tuple[str, ...]) -> dict[str, np.ndarray]:
    """Which pattern matched, and the indices it matched. Reported, not merged."""
    names = np.asarray([str(n) for n in names])
    out: dict[str, np.ndarray] = {}
    for pattern in patterns:
        hit = np.flatnonzero(np.array([bool(re.match(pattern, n)) for n in names]))
        if hit.size:
            out[pattern] = hit
    return out


def control_features(var_names: np.ndarray | list[str]) -> dict[str, object]:
    """Locate the negative probes and false codes, and say which convention hit.

    THE FLOOR IS BUILT FROM THESE. An empty match returns zero controls, a zero
    floor, and a gate every gene passes — a check that cannot fail, which is
    this repository's signature defect. So the result carries the matched
    patterns and the caller is expected to refuse an empty one.
    """
    negatives = _match_any(var_names, NEGATIVE_PROBE_PATTERNS)
    false_codes = _match_any(var_names, FALSE_CODE_PATTERNS)
    neg_idx = (np.unique(np.concatenate(list(negatives.values())))
               if negatives else np.array([], dtype=int))
    false_idx = (np.unique(np.concatenate(list(false_codes.values())))
                 if false_codes else np.array([], dtype=int))
    return {
        "negative_patterns_matched": sorted(negatives),
        "false_code_patterns_matched": sorted(false_codes),
        "negative_indices": neg_idx,
        "false_code_indices": false_idx,
        "n_negative": int(neg_idx.size),
        "n_false_code": int(false_idx.size),
    }


def require_controls(controls: dict[str, object], *, minimum: int = 5) -> None:
    """Refuse a run whose false-positive floor would be built from nothing."""
    if int(controls["n_negative"]) < minimum:
        raise CrowellError(
            f"only {controls['n_negative']} negative probes matched "
            f"{list(NEGATIVE_PROBE_PATTERNS)}; the paper reports 50. The "
            f"per-cell false-positive floor is measured from these, so an "
            f"empty or near-empty match gives a floor of zero and a gate every "
            f"gene clears. Inspect var_names and add the convention this "
            f"release uses before running the gate."
        )


#: Where this deposit actually keeps its false-positive floor. The control
#: probes are NOT features of X — they were summarised per cell into obs before
#: the object was written, so `control_features` finds nothing in var_names and
#: `require_controls` correctly refuses. The floor is still measurable, just
#: from here.
OBS_NEGATIVE_COUNT = "nCount_negprobes"
OBS_NEGATIVE_FEATURES = "nFeature_negprobes"
OBS_FALSECODE_COUNT = "nCount_falsecode"
OBS_FALSECODE_FEATURES = "nFeature_falsecode"
#: Deposit-defined QC inclusion flag. The Zenodo record describes ``fil`` as
#: "logical flag indicating whether or not a cell passed quality control".
#: A production read over every row silently mixes retained cells with cells
#: the data producer rejected, so there is deliberately no missing-column
#: fallback.
OBS_QC_PASS = "fil"


def qc_pass_mask(obs: pd.DataFrame, *, column: str = OBS_QC_PASS) -> np.ndarray:
    """Return the deposit's QC-pass mask, refusing anything non-logical.

    AnnData may round-trip logical metadata as booleans or as the strings
    ``"True"``/``"False"``. Both encodings are accepted; missing values and
    any third value are refused rather than treated as failures or passes.
    """
    if column not in obs.columns:
        raise CrowellError(
            f"obs has no {column!r}. The Crowell deposit defines it as the "
            "cell-level QC-pass flag; do not run the feasibility gate on an "
            "unfiltered population."
        )
    values = obs[column]
    if values.isna().any():
        raise CrowellError(
            f"obs[{column!r}] contains missing values; QC inclusion is undefined "
            "for those cells."
        )
    logical = values.astype(str).str.strip().str.lower()
    unexpected = sorted(set(logical) - {"true", "false"})
    if unexpected:
        raise CrowellError(
            f"obs[{column!r}] is not logical; unexpected values {unexpected}. "
            "Do not guess which cells passed QC."
        )
    mask = logical.eq("true").to_numpy(dtype=bool)
    if not mask.any():
        raise CrowellError(f"obs[{column!r}] marks zero cells as passing QC.")
    return mask


def floor_from_obs(obs: pd.DataFrame, *, n_negative_probes: int | None = None
                   ) -> dict[str, object]:
    """The per-probe false-positive floor, from obs summaries rather than var.

    ``nFeature_negprobes`` is the number of DISTINCT negative probes detected in
    a cell, so the mean of it divided by the probe count is the mean per-probe
    detection rate — the same quantity :func:`negative_floor` computes from a
    feature matrix, and the one comparable to a single gene's detection.

    ``n_negative_probes`` is inferred as the maximum observed if not supplied.
    **That is deliberately conservative**: an underestimated probe count divides
    by too little, which makes the floor too HIGH and the separation too SMALL,
    so a gene that clears it clears it for real. The paper reports 50; passing
    it explicitly is better and the sidecar records which was used.
    """
    missing = [c for c in (OBS_NEGATIVE_COUNT, OBS_NEGATIVE_FEATURES)
               if c not in obs.columns]
    if missing:
        raise CrowellError(
            f"obs carries no {missing}; the false-positive floor cannot be "
            f"built from it either. Do not proceed with a floor of zero — that "
            f"is a gate every gene clears."
        )
    features = pd.to_numeric(obs[OBS_NEGATIVE_FEATURES], errors="coerce")
    observed_max = int(np.nanmax(features.to_numpy())) if len(features) else 0
    n_probes = int(n_negative_probes or observed_max)
    if n_probes <= 0:
        raise CrowellError(
            f"cannot establish a negative-probe count ({OBS_NEGATIVE_FEATURES} "
            f"max is {observed_max}); the floor would divide by zero."
        )
    per_probe = float(np.nanmean(features.to_numpy()) / n_probes)
    out: dict[str, object] = {
        "floor_per_probe_mean": per_probe,
        "n_control_probes": n_probes,
        "probe_count_source": ("supplied" if n_negative_probes
                               else f"inferred from max({OBS_NEGATIVE_FEATURES})"),
        "probe_count_is_conservative": n_negative_probes is None,
        "any_probe_union_rate": float(np.nanmean((features > 0).to_numpy())),
        "mean_negative_counts_per_cell": float(
            np.nanmean(pd.to_numeric(obs[OBS_NEGATIVE_COUNT],
                                     errors="coerce").to_numpy())),
    }
    if OBS_FALSECODE_FEATURES in obs.columns:
        fc = pd.to_numeric(obs[OBS_FALSECODE_FEATURES], errors="coerce")
        out["false_code_union_rate"] = float(np.nanmean((fc > 0).to_numpy()))
    return out


#: Per-cell totals the deposit already carries. Used for the DEPTH ROW ONLY --
#: never for detection and never for the floor, both of which are computed from
#: the matrix itself.
OBS_TOTAL_COUNTS = "nCount_RNA"
OBS_TOTAL_FEATURES = "nFeature_RNA"


def read_gene_columns(path: str | Path, col_indices: dict[str, int], *,
                      min_umi: int = 1) -> dict[str, np.ndarray]:
    """Per-cell detection for a handful of genes, read straight from the h5ad.

    **Why this exists.** ``adata.to_memory()`` on a backed CSC materialises
    every non-zero as float64: for Crowell section 110 that is 549M values,
    4.39 GB of data plus 2.19 GB of indices, and it killed the run. The job
    needs **five columns of 18,878** — 0.027% of the matrix — and CSC stores a
    column contiguously, so reading five of them is cheap.

    Returns one boolean vector per gene, length ``n_obs``. Nothing dense of
    matrix size is ever allocated.

    CSR is refused rather than served slowly: extracting columns from CSR means
    touching every row, and the fallback that "works" would quietly reintroduce
    the allocation this function exists to avoid.
    """
    import h5py

    out: dict[str, np.ndarray] = {}
    with h5py.File(str(path), "r") as handle:
        node = handle["X"]
        encoding = node.attrs.get("encoding-type", "")
        if isinstance(encoding, bytes):
            encoding = encoding.decode()
        shape = tuple(node.attrs["shape"]) if "shape" in node.attrs else None
        if encoding != "csc_matrix":
            raise CrowellError(
                f"X is {encoding!r}; this reader handles 'csc_matrix', where a "
                f"column is contiguous and five of them are cheap. Extracting "
                f"columns from CSR means touching every row, and a fallback "
                f"that 'works' would reintroduce the multi-GB allocation this "
                f"function exists to avoid."
            )
        n_obs = int(shape[0]) if shape else int(node["indptr"].shape[0])
        indptr = node["indptr"][:]
        for gene, column in col_indices.items():
            start, stop = int(indptr[column]), int(indptr[column + 1])
            hit = np.zeros(n_obs, dtype=bool)
            if stop > start:
                rows = node["indices"][start:stop]
                vals = node["data"][start:stop]
                hit[rows[np.asarray(vals) >= min_umi]] = True
            out[gene] = hit
    return out


def depth_from_obs(obs: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Per-cell total counts and detected features, from the deposit's own obs.

    **Descriptive only.** The depth row exists so a reader can see whether two
    domains differ in capture before reading anything else; it never enters a
    detection rate, a floor or a verdict. Taken from obs rather than by summing
    the matrix because summing means loading it, and loading it is what
    exhausted memory on section 110.

    Raises if either column is absent — a depth row silently filled with zeros
    would read as "these domains have identical capture", which is a claim.
    """
    missing = [c for c in (OBS_TOTAL_COUNTS, OBS_TOTAL_FEATURES)
               if c not in obs.columns]
    if missing:
        raise CrowellError(
            f"obs carries no {missing}, so the per-domain depth row cannot be "
            f"built without loading the whole matrix. Do not substitute zeros: "
            f"a flat depth row reads as 'these domains have identical capture'."
        )
    counts = pd.to_numeric(obs[OBS_TOTAL_COUNTS], errors="coerce").to_numpy(float)
    features = pd.to_numeric(obs[OBS_TOTAL_FEATURES], errors="coerce").to_numpy(float)
    return counts, features


def check_depth_is_consistent(detection: dict[str, np.ndarray],
                              features_per_cell: np.ndarray) -> dict[str, object]:
    """A cell cannot detect more panel genes than it detects genes in total.

    Weak, cheap and non-vacuous: it uses only the columns already read, and it
    fires if ``nFeature_RNA`` describes a different feature space from ``X``.
    A stronger check means summing the matrix, which is the thing being avoided.
    """
    stacked = np.vstack([v.astype(int) for v in detection.values()])
    panel_hits = stacked.sum(axis=0)
    violations = int((panel_hits > features_per_cell).sum())
    return {
        "cells_checked": int(panel_hits.size),
        "violations": violations,
        "consistent": violations == 0,
        "max_panel_hits": int(panel_hits.max()) if panel_hits.size else 0,
    }


def domain_vocabulary(obs: pd.DataFrame, *, max_distinct: int = 50) -> dict[str, object]:
    """Every LOW-CARDINALITY obs column and its values. MAPS NOTHING.

    The pre-registration expects REF / TVA / CRC. Whether those words are in
    this object is a measurement, and this function reports it rather than
    assuming it — `becker_feasibility --inspect` exists for the same reason and
    found a fifth arm nobody had planned for.

    **THE SELECTION IS BY CARDINALITY, NOT BY NAME, AND THAT IS A CORRECTION.**
    The first version filtered against :data:`DOMAIN_COLUMN_CANDIDATES`, a
    hardcoded list of readable English names. The Crowell deposit names its
    annotations ``typ``, ``roi``, ``ctx``, ``ist``, ``lv1``, ``lv2``, ``trj``,
    ``jst`` — none of which is in that list, so the first inspection reported
    exactly one candidate column and silently hid eight others. A module whose
    docstring says "assume nothing" cannot select by guessing what a column
    will be called. Any column with few enough distinct values to be a
    vocabulary is reported; the name list is kept only to mark the ones that
    were *expected*.
    """
    found: dict[str, object] = {}
    for column in obs.columns:
        try:
            values = obs[column].astype(str)
        except Exception:  # pragma: no cover - exotic dtypes
            continue
        n_distinct = int(values.nunique())
        if (
            n_distinct > max_distinct
            or n_distinct <= 1 and column.lower() not in DOMAIN_COLUMN_CANDIDATES
        ):
            # A constant column is not a vocabulary — unless it is one of the
            # named candidates, where "all cells are TVA" is the finding.
            if not (n_distinct <= 1 and column.lower() in DOMAIN_COLUMN_CANDIDATES):
                continue
        counts = values.value_counts()
        found[column] = counts.to_dict()

    matches = {
        column: sorted({w for w in EXPECTED_DOMAIN_WORDS
                        for v in (values if isinstance(values, dict) else {})
                        if w in str(v).lower()})
        for column, values in found.items()
    }
    return {
        "candidate_columns": found,
        "expected_words_seen": {c: m for c, m in matches.items() if m},
        "expected_by_name": [c for c in found if c.lower() in DOMAIN_COLUMN_CANDIDATES],
        "all_obs_columns": list(obs.columns),
    }


def check_counts_are_integers(matrix, *, n_sample: int = 20_000,
                              seed: int = 0) -> dict[str, object]:
    """Refuse a matrix that has already been normalised.

    A detection rate — "fraction of cells with >= 1 count" — is meaningless on
    log-normalised or scaled values, and it does not raise on its own: it
    returns a plausible number. Sampled rather than exhaustive because these
    matrices are millions of cells wide.
    """
    rng = np.random.default_rng(seed)
    n_rows = matrix.shape[0]
    rows = rng.choice(n_rows, size=min(n_sample, n_rows), replace=False)
    block = matrix[np.sort(rows)]
    values = (block.data if hasattr(block, "data")
              else np.asarray(block).ravel())
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return {"integral": True, "n_checked": 0, "note": "no non-zero values sampled"}
    integral = bool(np.allclose(values, np.rint(values)))
    return {
        "integral": integral,
        "n_checked": int(values.size),
        "max": float(values.max()),
        "min_nonzero": float(values[values > 0].min()) if (values > 0).any() else 0.0,
    }


def require_counts(report: dict[str, object]) -> None:
    if not report.get("integral", False):
        raise CrowellError(
            f"the matrix is not integer counts (max {report.get('max')}, "
            f"min non-zero {report.get('min_nonzero')}). A detection rate off "
            f"normalised values is not a detection rate and it will not raise "
            f"on its own — find the raw layer, or say in the sidecar which "
            f"layer was used and why."
        )
