"""W3.8 — independent replication of the premise check. GSE39582.

WHY A SECOND COHORT, AND WHY THIS ONE
-------------------------------------
W3.2 found bulk GUCA2A loss to be continuous rather than bimodal, and the team
is about to reshape the bulk arm around that. A finding that redirects a project
should be replicated *before* the pivot, not after.

GSE39582 (Marisa et al., French CIT consortium) is the strongest available test
of the specific things that could be wrong:

- **Different patients**, different country, different collection era.
- **Different measurement technology.** Affymetrix HG-U133 Plus 2.0 microarray
  against TCGA's Illumina RNA-seq. Different noise, different dynamic range,
  and — critically — a different detection floor. An artifact of one platform
  will not survive the other.
- **Richly annotated**: MMR status, and *CIMP status*, which TCGA's GDC clinical
  does not carry.

CIMP IS THE REASON THIS COHORT IS WORTH THE EFFORT
---------------------------------------------------
If a discrete "off" state exists anywhere, it is in the CpG island methylator
phenotype — promoter hypermethylation is a switch, not a dial. TCGA let us
stratify by MSI; this cohort lets us stratify by the methylation phenotype
itself. That is a mechanistically motivated place to look for the bimodality
W3.2 did not find, rather than a fishing expedition.

INVARIANT 4
-----------
*Estimate per study, then meta-analyse. Never pool.* The two cohorts are never
concatenated. Each is analysed independently with **the same test code** — this
module reuses ``premise.assess`` rather than reimplementing it, so a difference
between cohorts cannot come from a difference in the analysis.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd

GEO_SERIES_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE39nnn/GSE39582/matrix/"
    "GSE39582_series_matrix.txt.gz"
)
GEO_PLATFORM_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc=GPL570&targ=self&form=text&view=data"
)

#: Sample characteristics worth keeping. GEO stores them as "key: value".
KEEP_CHARACTERISTICS = (
    "dataset",
    "Sex",
    "age.at.diagnosis (year)",
    "tnm.stage",
    "tumor.location",
    "mmr.status",
    "cimp.status",
    "os.event",
    "os.delay (months)",
    "rfs.event",
    "rfs.delay",
)

#: GEO writes missing as this. Left alone it becomes a category.
GEO_MISSING = "N/A"


class ReplicationError(RuntimeError):
    """The replication cohort could not be loaded."""


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------


def fetch(url: str, target: str | Path) -> Path:
    """Download once. Idempotent."""
    import requests

    target = Path(target)
    if target.exists() and target.stat().st_size > 0:
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=900) as response:
        response.raise_for_status()
        tmp = target.with_suffix(target.suffix + ".partial")
        with open(tmp, "wb") as handle:
            for chunk in response.iter_content(chunk_size=1 << 20):
                handle.write(chunk)
        tmp.replace(target)
    return target


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_series_matrix(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read a GEO series matrix. Returns (expression probes x samples, metadata).

    The characteristics block is positional: every ``!Sample_characteristics_ch1``
    line carries one field for all samples in column order, and the field name is
    the part before the colon. Parsing by position rather than by matching names
    per sample is what keeps a sample's MMR status attached to that sample.
    """
    path = Path(path)
    characteristics: dict[str, list[str]] = {}
    accessions: list[str] = []
    titles: list[str] = []

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Sample_geo_accession"):
                accessions = _values(line)
            elif line.startswith("!Sample_title"):
                titles = _values(line)
            elif line.startswith("!Sample_characteristics_ch1"):
                values = _values(line)
                if not values:
                    continue
                key = values[0].split(":", 1)[0].strip()
                cleaned = [
                    v.split(": ", 1)[1].strip() if ": " in v else v.strip() for v in values
                ]
                characteristics.setdefault(key, cleaned)
        else:
            raise ReplicationError(f"{path} has no expression table")

        expression = pd.read_csv(
            handle, sep="\t", index_col=0, comment="!", low_memory=False
        )

    if not accessions:
        raise ReplicationError("no sample accessions found")
    metadata = pd.DataFrame({"geo_accession": accessions})
    if titles:
        metadata["title"] = titles
    for key in KEEP_CHARACTERISTICS:
        if key in characteristics:
            metadata[key] = characteristics[key]
    metadata = metadata.replace(GEO_MISSING, np.nan).set_index("geo_accession")

    expression.index.name = "probe"
    expression.columns = [c.strip().strip('"') for c in expression.columns]
    shared = [c for c in expression.columns if c in set(metadata.index)]
    if len(shared) < len(metadata) * 0.9:
        raise ReplicationError(
            f"only {len(shared)} of {len(metadata)} samples matched between the "
            f"metadata and the expression table"
        )
    return expression.loc[:, shared], metadata.loc[shared]


def _values(line: str) -> list[str]:
    return [v.strip().strip('"') for v in line.rstrip("\n").split("\t")[1:]]


def parse_platform_table(path: str | Path) -> pd.Series:
    """Probe id -> gene symbol from a GEO platform table."""
    path = Path(path)
    header_line = None
    with open(path, encoding="utf-8", errors="replace") as handle:
        for i, line in enumerate(handle):
            if line.startswith("ID\t"):
                header_line = i
                break
    if header_line is None:
        raise ReplicationError(f"{path} has no platform table header")

    table = pd.read_csv(
        path, sep="\t", skiprows=header_line, usecols=["ID", "Gene Symbol"],
        dtype=str, low_memory=False,
    )
    table = table.dropna(subset=["Gene Symbol"])
    # Affymetrix writes multi-mapping probes as "A /// B". A probe that cannot be
    # attributed to one gene is dropped rather than assigned to the first.
    single = table.loc[~table["Gene Symbol"].str.contains("///", regex=False)]
    return single.set_index("ID")["Gene Symbol"].str.strip()


def collapse_probes_to_symbols(
    expression: pd.DataFrame, probe_to_symbol: pd.Series
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Probes x samples -> genes x samples, keeping the highest-mean probe.

    Same rule as ``purity.to_symbol_matrix`` uses for duplicate Ensembl IDs:
    averaging would dilute a real signal with a dead probe, summing would
    inflate it, and picking the dominant probe invents nothing.
    """
    counts = {"probes": int(len(expression))}
    mapped = expression.loc[expression.index.intersection(probe_to_symbol.index)]
    counts["probes_with_a_unique_symbol"] = int(len(mapped))

    symbols = probe_to_symbol.loc[mapped.index]
    order = np.argsort(-mapped.mean(axis=1).to_numpy(), kind="stable")
    mapped = mapped.iloc[order]
    symbols = symbols.iloc[order]
    keep = ~symbols.duplicated(keep="first")
    out = mapped.loc[keep.to_numpy()]
    out.index = pd.Index(symbols.loc[keep.to_numpy()].to_numpy(), name="gene_symbol")
    counts["unique_genes"] = int(len(out))
    return out.sort_index(), counts


# ---------------------------------------------------------------------------
# Strata — the mechanistically motivated ones
# ---------------------------------------------------------------------------


def strata(metadata: pd.DataFrame) -> dict[str, np.ndarray]:
    """Boolean masks over samples, in the same spirit as the TCGA strata.

    ``dataset`` separates the 566 tumours from the 19 non-tumoural mucosa
    samples. The CIMP and MMR splits are the ones this cohort exists to add.
    """
    tumour = metadata["dataset"].isin(["discovery", "validation"]).to_numpy()
    out: dict[str, np.ndarray] = {
        "tumour": tumour,
        "normal_mucosa": (metadata["dataset"] == "Non Tumoral").to_numpy(),
    }
    if "mmr.status" in metadata:
        out["tumour|dMMR"] = tumour & (metadata["mmr.status"] == "dMMR").to_numpy()
        out["tumour|pMMR"] = tumour & (metadata["mmr.status"] == "pMMR").to_numpy()
    if "cimp.status" in metadata:
        out["tumour|CIMP+"] = tumour & (metadata["cimp.status"] == "+").to_numpy()
        out["tumour|CIMP-"] = tumour & (metadata["cimp.status"] == "-").to_numpy()
    if "tumor.location" in metadata:
        out["tumour|proximal"] = tumour & (metadata["tumor.location"] == "proximal").to_numpy()
        out["tumour|distal"] = tumour & (metadata["tumor.location"] == "distal").to_numpy()
    return out


def fold_change_vs_normal(
    expression: pd.DataFrame, masks: dict[str, np.ndarray], gene: str
) -> dict[str, float]:
    """Tumour/normal fold change on a log2 array, for comparison with TCGA.

    Both cohorts are log2, so a difference of means is a log2 fold change and the
    ratios are comparable even though the underlying units are not.
    """
    if gene not in expression.index:
        return {"gene": gene, "fold_change": float("nan")}
    values = expression.loc[gene]
    tumour = float(values[masks["tumour"]].median())
    normal = float(values[masks["normal_mucosa"]].median())
    return {
        "gene": gene,
        "tumour_median_log2": round(tumour, 3),
        "normal_median_log2": round(normal, 3),
        "log2_fold_change": round(tumour - normal, 3),
        "fold_change": round(float(2 ** (tumour - normal)), 4),
    }
