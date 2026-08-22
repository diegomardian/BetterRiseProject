"""W3.3 — tumour purity. ESTIMATE, ABSOLUTE, and the gap between them.

WHY THIS IS PYTHON AND NOT R
----------------------------
ESTIMATE ships as an R package (GPL-2, Yoshihara/Kim/Verhaak, R-Forge v1.0.13).
R does not run in this environment, so :func:`estimate_scores` reimplements
``estimateScore.R`` directly. That is a real risk — a reimplementation that
quietly differs from the canonical tool is worse than no tool — so it is
mitigated the only way that counts: **the package ships a worked example with
its own reference output, and `tests/test_bulk_purity.py` asserts we reproduce
it to six decimal places.** If that test passes, this is ESTIMATE.

The signatures and the common-gene list are *inputs*, downloaded to
``data/raw/estimate/`` and recorded in the manifest like any other external
data. They are deliberately **not** committed: they are GPL-2 material from
someone else's package, and vendoring them into this repo raises a licensing
question nobody needs.

THE PURITY FORMULA IS AFFYMETRIX-ONLY. THIS MATTERS
---------------------------------------------------
``estimateScore.R`` computes ``TumorPurity`` **only** when
``platform == "affymetrix"``. For agilent and illumina it emits the three
scores and stops. The conversion

    purity = cos(0.6049872018 + 0.0001467884 * ESTIMATEScore)

was fit on Affymetrix arrays against ABSOLUTE, and the authors declined to
extend it. TCGA STAR counts are Illumina RNA-seq. Applying the formula anyway is
common in the literature and is an extrapolation off the platform it was
calibrated on.

So :func:`estimate_scores` returns the three scores for any platform, and
:func:`affymetrix_purity` is a separate call you have to make deliberately. What
it returns is labelled ``estimate_affy_extrapolated``, never plain "purity".

AND ESTIMATE IS PARTLY CIRCULAR WITH OUR OUTCOME
------------------------------------------------
ESTIMATE infers purity *from expression* — stromal and immune signature
enrichment. Our outcome is also expression. Conditioning a GUCA2A analysis on an
expression-derived covariate removes some of the signal along with the
confounder, and the direction is not knowable in advance.

**ABSOLUTE is the more independent estimate** — it is called from copy-number
and SNP data, not expression. Where both exist, prefer ABSOLUTE. The ``method``
column exists so the two are never silently coalesced, which is what the brief
asks for and what makes the sensitivity analysis possible at all.
"""

from __future__ import annotations

import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

#: R-Forge source package. Pinned by version so the signatures cannot move
#: under us; recorded in data/manifest.csv with a sha256 when downloaded.
ESTIMATE_VERSION = "1.0.13"
ESTIMATE_URL = f"http://r-forge.r-project.org/src/contrib/estimate_{ESTIMATE_VERSION}.tar.gz"

#: Coefficients from estimateScore.R. Fit on Affymetrix. See module docstring.
AFFY_PURITY_INTERCEPT = 0.6049872018
AFFY_PURITY_SLOPE = 0.0001467884

SIGNATURE_NAMES = ("StromalSignature", "ImmuneSignature")


class PurityError(RuntimeError):
    """Purity could not be computed or joined."""


# ---------------------------------------------------------------------------
# The ESTIMATE inputs, as downloaded data
# ---------------------------------------------------------------------------


def fetch_estimate_package(dest_dir: str | Path) -> Path:
    """Download the ESTIMATE source tarball. Idempotent. Returns the path."""
    import requests

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / f"estimate_{ESTIMATE_VERSION}.tar.gz"
    if target.exists() and target.stat().st_size > 0:
        return target
    response = requests.get(ESTIMATE_URL, timeout=300)
    response.raise_for_status()
    target.write_bytes(response.content)
    return target


def _extract(tarball: str | Path, member_suffix: str) -> str:
    with tarfile.open(tarball, "r:gz") as tf:
        for member in tf.getmembers():
            if member.name.endswith(member_suffix):
                fh = tf.extractfile(member)
                if fh is None:
                    break
                return fh.read().decode("utf-8")
    raise PurityError(f"{member_suffix} not found in {tarball}")


def load_signatures(tarball: str | Path) -> dict[str, list[str]]:
    """The two 141-gene signatures from ``SI_geneset.gmt``."""
    text = _extract(tarball, "inst/extdata/SI_geneset.gmt")
    sets: dict[str, list[str]] = {}
    for line in text.strip().splitlines():
        fields = [f for f in line.rstrip("\n").split("\t") if f]
        sets[fields[0]] = fields[2:]
    missing = [n for n in SIGNATURE_NAMES if n not in sets]
    if missing:
        raise PurityError(f"signature(s) {missing} missing from SI_geneset.gmt")
    return sets


def load_common_genes(tarball: str | Path) -> list[str]:
    """The 10,412 genes ESTIMATE restricts to before scoring."""
    text = _extract(tarball, "inst/extdata/common_genes.txt")
    rows = [line.split("\t") for line in text.strip().splitlines()]
    header, body = rows[0], rows[1:]
    col = header.index("GeneSymbol")
    return [r[col] for r in body]


# ---------------------------------------------------------------------------
# The algorithm — a direct port of estimateScore.R
# ---------------------------------------------------------------------------


def filter_common_genes(expression: pd.DataFrame, common_genes: list[str]) -> pd.DataFrame:
    """Restrict a genes x samples matrix to ESTIMATE's common-gene set.

    A port of ``filterCommonGenes.R``, which does ``merge(common_genes, input,
    by="GeneSymbol")``. Two details are load-bearing and neither is obvious:

    - **The scores are defined on the filtered matrix, not the raw one.** The
      first step of ``estimateScore`` divides ranks by the gene count, so
      running on 17,256 genes instead of 10,412 shifts every score. The
      package's own reference output is post-filter.
    - **R's ``merge`` sorts by the join column.** Row order changes the
      tie-breaking inside ``order()`` and therefore the enrichment score, so the
      result is sorted alphabetically to match.
    """
    keep = expression.index.intersection(pd.Index(sorted(set(common_genes))))
    return expression.loc[keep].sort_index()


def _enrichment_score(values: np.ndarray, in_set: np.ndarray) -> float:
    """One sample, one gene set. The ssGSEA integral, exactly as ESTIMATE does it.

    Note ``ES = sum(RES)`` rather than the maximum deviation used by classical
    GSEA. That is ESTIMATE's choice and it is why these scores run to thousands
    rather than to one.
    """
    order = np.argsort(-values, kind="stable")  # R: order(x, decreasing=TRUE)
    correl = np.abs(values[order]) ** 0.25
    tag = in_set[order].astype(float)
    n_hit = float(tag.sum())
    n_miss = float(values.size - n_hit)
    if n_hit == 0 or n_miss == 0:
        return float("nan")
    sum_correl = float(correl[tag == 1].sum())
    f_miss = np.cumsum((1.0 - tag) / n_miss)
    f_hit = np.cumsum(tag * correl / sum_correl)
    return float(np.sum(f_hit - f_miss))


def estimate_scores(
    expression: pd.DataFrame,
    signatures: dict[str, list[str]],
) -> pd.DataFrame:
    """Stromal, immune and ESTIMATE scores. One row per sample.

    Parameters
    ----------
    expression:
        **genes (rows) x samples (columns)**, indexed by HGNC symbol, already
        restricted to the common-gene set. Any monotone transform of expression
        gives identical scores — the first step is a per-sample rank — so TPM and
        log2-CPM agree here by construction. That is a property worth knowing
        rather than a licence to be careless: see ``test_scale_invariance``.

    Returns
    -------
    DataFrame indexed by sample with StromalScore, ImmuneScore, ESTIMATEScore.
    No purity column. Purity is :func:`affymetrix_purity`, called deliberately.
    """
    if expression.empty:
        raise PurityError("expression matrix is empty")

    # Per-sample rank normalisation, then scale by gene count. estimateScore.R.
    n_genes = expression.shape[0]
    ranked = np.apply_along_axis(
        lambda col: rankdata(col, method="average"), 0, expression.to_numpy(dtype=float)
    )
    ranked = 10000.0 * ranked / n_genes

    gene_index = pd.Index(expression.index)
    out: dict[str, list[float]] = {}
    for name in SIGNATURE_NAMES:
        overlap = gene_index.isin(set(signatures[name]))
        if not overlap.any():
            raise PurityError(f"no genes from {name} are present in the matrix")
        out[name] = [
            _enrichment_score(ranked[:, j], overlap) for j in range(ranked.shape[1])
        ]

    scores = pd.DataFrame(
        {
            "StromalScore": out["StromalSignature"],
            "ImmuneScore": out["ImmuneSignature"],
        },
        index=pd.Index(expression.columns, name="barcode"),
    )
    scores["ESTIMATEScore"] = scores["StromalScore"] + scores["ImmuneScore"]
    return scores


def affymetrix_purity(estimate_score: pd.Series) -> pd.Series:
    """The Affymetrix-calibrated purity conversion, applied off-platform.

    Returns NaN where the conversion goes negative, as estimateScore.R does.
    The name is long on purpose. This is not "purity"; it is one platform's
    calibration extrapolated onto another, and every caller should be made to
    look at that.
    """
    purity = np.cos(AFFY_PURITY_INTERCEPT + AFFY_PURITY_SLOPE * estimate_score)
    return purity.where(purity >= 0)


# ---------------------------------------------------------------------------
# Preparing our matrix for ESTIMATE
# ---------------------------------------------------------------------------


def to_symbol_matrix(
    expression: pd.DataFrame,
    gene_map: pd.DataFrame,
    common_genes: list[str],
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Samples x Ensembl -> genes x symbol, restricted to ESTIMATE's common set.

    Duplicate symbols are collapsed to the **highest mean-expression** Ensembl
    ID rather than summed or averaged. Summing would inflate a symbol that
    happens to carry a pseudogene; averaging would dilute it. Picking the
    dominant transcript is deterministic, does not invent a value, and matters
    little because everything downstream is rank-based.

    Returns the matrix and a step-by-step gene count, because "how many of the
    10,412 did we actually match" is the number that says whether the scores
    are comparable to published ones.
    """
    counts = {"ensembl_in_matrix": expression.shape[1]}

    symbols = gene_map.set_index("ensembl_id")["gene_symbol"]
    usable = expression.columns.intersection(symbols.index)
    sub = expression.loc[:, usable]
    sub.columns = symbols.loc[usable].to_numpy()
    counts["mapped_to_a_symbol"] = sub.shape[1]

    keep = sub.columns.isin(set(common_genes))
    sub = sub.loc[:, keep]
    counts["in_estimate_common_genes"] = sub.shape[1]

    # Collapse duplicate symbols: keep the highest-mean column for each.
    means = sub.mean(axis=0)
    order = np.argsort(-means.to_numpy(), kind="stable")
    sub = sub.iloc[:, order]
    sub = sub.loc[:, ~sub.columns.duplicated(keep="first")]
    counts["unique_symbols"] = sub.shape[1]
    counts["common_genes_total"] = len(set(common_genes))

    # filter_common_genes does the alphabetical sort that R's merge implies.
    return filter_common_genes(sub.T, common_genes), counts


# ---------------------------------------------------------------------------
# ABSOLUTE, and the comparison
# ---------------------------------------------------------------------------


#: Aran's Supplementary Data 1 carries a two-line title before the header.
ARAN_HEADER_ROW = 3
ARAN_URL = (
    "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9971"
    "/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx"
)

#: ABSOLUTE is called from DNA, our expression is RNA, so the two never share an
#: aliquot barcode. They do share the sample: TCGA-A6-2670-01A. Sixteen
#: characters — patient, sample type and vial — is the finest key on which a
#: DNA and an RNA aliquot of the same tumour agree.
SAMPLE_KEY_LEN = 16


def sample_key(barcode: str) -> str:
    """``TCGA-A6-2670-01A-01R-1410-07`` -> ``TCGA-A6-2670-01A``."""
    return str(barcode)[:SAMPLE_KEY_LEN]


def fetch_aran_table(dest_dir: str | Path) -> Path:
    """Download Aran et al. 2015 Supplementary Data 1. Idempotent."""
    import requests

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    target = dest_dir / "aran2015_ncomms9971_supplementary_data_1.xlsx"
    if target.exists() and target.stat().st_size > 0:
        return target
    response = requests.get(ARAN_URL, timeout=300)
    response.raise_for_status()
    target.write_bytes(response.content)
    return target


def load_aran_purity(path: str | Path) -> pd.DataFrame:
    """Aran et al. 2015 pan-cancer purity table, keyed at sample level.

    Columns are ``Sample ID, Cancer type, ESTIMATE, ABSOLUTE, LUMP, IHC, CPE``.
    We take ABSOLUTE (the brief's request, and the one estimate not derived from
    expression) and keep Aran's own ESTIMATE column too — it is an independent
    check on our port against real TCGA data rather than against the package's
    ten-sample toy.

    Note ABSOLUTE is sparsely populated in this table; that sparsity *is* the
    coverage number the brief asks for, so it is preserved rather than filled.
    """
    path = Path(path)
    raw = pd.read_excel(path, header=ARAN_HEADER_ROW)
    expected = {"Sample ID", "Cancer type", "ESTIMATE", "ABSOLUTE", "CPE"}
    missing = sorted(expected - set(raw.columns))
    if missing:
        raise PurityError(
            f"{path.name} is missing column(s) {missing}; found {list(raw.columns)}. "
            f"The header row may have moved — it is row {ARAN_HEADER_ROW} in the "
            f"published file."
        )
    out = pd.DataFrame(
        {
            "sample_key": raw["Sample ID"].astype(str).str.strip().map(sample_key),
            "cancer_type": raw["Cancer type"].astype(str).str.strip(),
            "aran_absolute": pd.to_numeric(raw["ABSOLUTE"], errors="coerce"),
            "aran_estimate": pd.to_numeric(raw["ESTIMATE"], errors="coerce"),
            "aran_cpe": pd.to_numeric(raw["CPE"], errors="coerce"),
        }
    )
    return out.drop_duplicates("sample_key").reset_index(drop=True)


def assemble_purity_table(
    scores: pd.DataFrame,
    absolute: pd.DataFrame | None,
    *,
    manifest: pd.DataFrame,
) -> pd.DataFrame:
    """Long-format purity: one row per (sample, method). Never coalesced.

    The brief is explicit — *"emit one purity value per sample with a method
    column recording its source. Do not silently coalesce the two."* A single
    "purity" column filled from ABSOLUTE where available and ESTIMATE elsewhere
    would make the coverage difference invisible and would mix an
    expression-derived covariate with an independent one inside the same model
    term.
    """
    rows = []
    affy = affymetrix_purity(scores["ESTIMATEScore"])
    for barcode in scores.index:
        rows.append(
            {
                "barcode": barcode,
                "method": "estimate_affy_extrapolated",
                "purity": float(affy.loc[barcode]) if pd.notna(affy.loc[barcode]) else None,
                "stromal_score": float(scores.loc[barcode, "StromalScore"]),
                "immune_score": float(scores.loc[barcode, "ImmuneScore"]),
                "estimate_score": float(scores.loc[barcode, "ESTIMATEScore"]),
                "expression_derived": True,
            }
        )

    if absolute is not None:
        # Join on the 16-char sample key: ABSOLUTE is called from DNA and our
        # expression is RNA, so the aliquot barcodes never match.
        keyed = absolute.set_index("sample_key")
        for barcode in scores.index:
            key = sample_key(barcode)
            if key not in keyed.index:
                continue
            row = keyed.loc[key]
            for method, column, derived in (
                ("absolute", "aran_absolute", False),
                ("aran_cpe", "aran_cpe", False),
            ):
                value = row.get(column)
                if pd.notna(value):
                    rows.append(
                        {
                            "barcode": barcode,
                            "method": method,
                            "purity": float(value),
                            "stromal_score": None,
                            "immune_score": None,
                            "estimate_score": None,
                            "expression_derived": derived,
                        }
                    )

    table = pd.DataFrame(rows)
    meta = manifest.reindex(table["barcode"])
    for col in ("patient_id", "project", "sample_type", "sample_type_name"):
        if col in meta.columns:
            table[col] = meta[col].to_numpy()
    return table.sort_values(["barcode", "method"]).reset_index(drop=True)


def agreement(table: pd.DataFrame) -> dict[str, float | int]:
    """Pearson and Spearman between the two methods, plus ABSOLUTE's coverage."""
    wide = table.pivot_table(index="barcode", columns="method", values="purity")
    n_estimate = int(wide.get("estimate_affy_extrapolated", pd.Series(dtype=float)).notna().sum())
    n_absolute = int(wide.get("absolute", pd.Series(dtype=float)).notna().sum())
    out: dict[str, float | int] = {
        "n_samples": int(len(wide)),
        "n_estimate": n_estimate,
        "n_absolute": n_absolute,
        "absolute_coverage": round(n_absolute / len(wide), 4) if len(wide) else float("nan"),
    }
    if "absolute" in wide.columns and "estimate_affy_extrapolated" in wide.columns:
        both = wide.dropna(subset=["absolute", "estimate_affy_extrapolated"])
        out["n_both"] = int(len(both))
        if len(both) >= 3:
            out["pearson_r"] = round(
                float(both["absolute"].corr(both["estimate_affy_extrapolated"])), 4
            )
            out["spearman_rho"] = round(
                float(
                    both["absolute"].corr(
                        both["estimate_affy_extrapolated"], method="spearman"
                    )
                ),
                4,
            )
            out["median_difference"] = round(
                float((both["absolute"] - both["estimate_affy_extrapolated"]).median()), 4
            )
    return out
