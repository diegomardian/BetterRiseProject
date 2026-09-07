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


def _match_any(names: "np.ndarray | list[str]", patterns: tuple[str, ...]) -> dict[str, np.ndarray]:
    """Which pattern matched, and the indices it matched. Reported, not merged."""
    names = np.asarray([str(n) for n in names])
    out: dict[str, np.ndarray] = {}
    for pattern in patterns:
        hit = np.flatnonzero(np.array([bool(re.match(pattern, n)) for n in names]))
        if hit.size:
            out[pattern] = hit
    return out


def control_features(var_names: "np.ndarray | list[str]") -> dict[str, object]:
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


def domain_vocabulary(obs: pd.DataFrame) -> dict[str, object]:
    """Every candidate domain column and its values. MAPS NOTHING.

    The pre-registration expects REF / TVA / CRC. Whether those words are in
    this object is a measurement, and this function reports it rather than
    assuming it — `becker_feasibility --inspect` exists for the same reason and
    found a fifth arm nobody had planned for.
    """
    found: dict[str, object] = {}
    for column in obs.columns:
        if column.lower() not in DOMAIN_COLUMN_CANDIDATES:
            continue
        values = obs[column].astype(str)
        counts = values.value_counts()
        if len(counts) > 50:
            found[column] = {"n_distinct": int(len(counts)), "too_many_to_list": True}
            continue
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
