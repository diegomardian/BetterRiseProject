"""Reading GSE201348 — 72 10x triplets in a tar, and metadata that is not in it.

THE SHAPE OF THIS DEPOSIT, verified 2026-09-06 against the GEO FTP listing and
the series matrix rather than assumed:

    GSE201348_RAW.tar              1.2 GB, 72 samples x 3 files
      GSM6061645_A001-C-007_barcodes.tsv.gz
      GSM6061645_A001-C-007_features.tsv.gz
      GSM6061645_A001-C-007_matrix.mtx.gz
    GSE201348_series_matrix.txt.gz  the ONLY place the arm labels live

**No metadata is in the tar.** A filename gives a GSM, a donor and a sample --
it does not say whether the sample is a polyp or unaffected mucosa. That lives
in the series matrix's ``disease stage`` characteristic, which is why this
module needs both files and refuses to guess from either alone.

WHY THE VOCABULARY IS A DICT AND UNKNOWN VALUES RAISE. ``ADENOMA_TISSUE_MAP``
was written against the ICBI atlas's words, and reading a label rather than the
patient grouping once put Chen_2021's usable pairs at **zero** when the true
number was 44. The lesson taken was not "be careful with labels", it was: the
mapping is a decision recorded in one place, and a value nobody has seen stops
the run instead of being silently dropped. A dropped sample is a smaller cohort
reported as the same cohort.

THE TWO ESTIMAND DECISIONS ARE NOT MADE HERE. Pooling per donor and pooling
technical replicates are fixed in ``docs/prereg_becker_replication.md``
Amendment 1, before this module existed. This implements them; it does not
choose them, and ``pool_by`` has no default for exactly that reason.
"""

from __future__ import annotations

import gzip
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

#: `disease stage` -> the arm, read off the series matrix on 2026-09-06.
#: `CRC` is present in the deposit and deliberately excluded: this is the
#: adenoma reading, and carcinoma is where four routes already terminated.
#: Excluding it is a decision, so it is spelled out rather than achieved by
#: omission.
DISEASE_STAGE_MAP: dict[str, str | None] = {
    "Polyp": "tumour",
    "Unaffected": "normal",
    "CRC": None,
    # A FIFTH VALUE, found in the deposit on 2026-09-06 and NOT merged with
    # `Unaffected`. `Normal` appears only in B001/B004 -- donors with FAP=N and
    # NO POLYPS. It is a separate healthy person's tissue, where `Unaffected` is
    # a FAP patient's own uninvolved mucosa from the same colon as their polyps.
    #
    # Merging them would silently convert a paired within-patient design into a
    # cross-donor one, which is the estimand drift Amendment 1 exists to
    # prevent. It maps to its own arm so that anything wanting it has to ask.
    "Normal": "healthy_donor",
}

#: The arms a PAIRED reading may use. `healthy_donor` is deliberately absent:
#: those donors have no polyps, so pairing against them is cross-donor by
#: construction. Becker Amendment 2 refuses that rescue explicitly.
PAIRED_ARMS: frozenset[str] = frozenset({"tumour", "normal"})

#: Characteristics the series matrix carries, and the column each becomes.
CHARACTERISTIC_FIELDS: dict[str, str] = {
    "disease stage": "disease_stage",
    "familial adenomatous_polyposis": "fap",
    "Sex": "sex",
    "tissue": "tissue",
}

#: FOUR naming schemes, verified against the sample table on 2026-09-06. The
#: donor is the unit of inference under invariant 5 and Amendment 1, so reading
#: it wrong misassigns every cell in a sample.
#:
#:   A001-C-007    FAP series          -> donor A001
#:   B001-A-301    healthy donors      -> donor B001
#:   F007          single polyps       -> donor F007 (the sample IS the donor)
#:   CRC1_8810     sporadic carcinoma  -> donor CRC1  (excluded anyway)
#:
#: Tried in order, first match wins. The first version of this carried only the
#: first pattern and REFUSED the other seven samples rather than dropping them,
#: which is why they were found at all.
_SAMPLE_PATTERNS: tuple[re.Pattern, ...] = (
    re.compile(r"^(?P<donor>[AB]\d+)-(?P<sample>[A-Z]-\d+)$"),
    re.compile(r"^(?P<donor>CRC\d+)_(?P<sample>\d+)$"),
    re.compile(r"^(?P<donor>F\d+[A-Z]?)$"),
)


def parse_sample_id(sample_id: str) -> str | None:
    """The donor a sample belongs to, or ``None`` if no scheme matches.

    ``None`` is not a fallback -- the caller raises on it. A sample whose donor
    cannot be read cannot be placed in a paired design, and guessing is how
    every cell in it gets attributed to the wrong person.
    """
    for pattern in _SAMPLE_PATTERNS:
        match = pattern.match(sample_id)
        if match:
            return match.group("donor")
    return None


class BeckerError(ValueError):
    """The deposit does not have the shape this module verified."""


def _unquote(value: str) -> str:
    return value.strip().strip('"').strip()


def read_series_matrix(path: str | Path) -> pd.DataFrame:
    """One row per GSM: donor, sample, replicate, arm, and the covariates.

    Parses the ``!Sample_*`` header lines. The expression block below them is
    empty for this series and is not read.

    **Unknown ``disease stage`` values raise.** Every sample must land in a
    known arm or be a recognised exclusion; a value nobody has seen means the
    deposit is not the one this was verified against.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    titles: list[str] = []
    geo_ids: list[str] = []
    characteristics: list[list[str]] = []

    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.partition("\t")
            values = [_unquote(v) for v in rest.rstrip("\n").split("\t")]
            if key == "!Sample_title":
                titles = values
            elif key == "!Sample_geo_accession":
                geo_ids = values
            elif key == "!Sample_characteristics_ch1":
                characteristics.append(values)

    if not titles:
        raise BeckerError(
            f"{path} has no !Sample_title line. This is meant to be a GEO "
            f"series matrix; if it is not, the arm labels are somewhere else "
            f"and nothing downstream can proceed."
        )
    if geo_ids and len(geo_ids) != len(titles):
        raise BeckerError(
            f"{len(geo_ids)} accessions against {len(titles)} titles. The "
            f"series matrix is column-aligned by sample, so a mismatch means "
            f"every row would carry another sample's metadata."
        )

    frame = pd.DataFrame({"title": titles})
    if geo_ids:
        frame["gsm"] = geo_ids

    for row in characteristics:
        if len(row) != len(titles):
            continue
        parsed = [v.split(": ", 1) for v in row]
        names = {p[0] for p in parsed if len(p) == 2}
        if len(names) != 1:
            continue
        name = names.pop()
        column = CHARACTERISTIC_FIELDS.get(name)
        if column:
            frame[column] = [p[1] if len(p) == 2 else None for p in parsed]

    missing = sorted(set(CHARACTERISTIC_FIELDS.values()) - set(frame.columns))
    if "disease_stage" in missing:
        raise BeckerError(
            f"no 'disease stage' characteristic in {path}. That field is the "
            f"ONLY place the polyp/unaffected distinction lives in this "
            f"deposit — the tar carries none of it."
        )

    # "A001-C-007, Replicate2, snRNAseq" -> sample_id, replicate
    parts = frame["title"].str.split(",").apply(lambda p: [x.strip() for x in p])
    frame["sample_id"] = parts.str[0]
    frame["replicate"] = parts.apply(
        lambda p: next((x for x in p[1:] if x.lower().startswith("replicate")), None)
    )

    frame["donor"] = frame["sample_id"].map(parse_sample_id)
    unparsed = frame.loc[frame["donor"].isna(), "sample_id"].tolist()
    if unparsed:
        raise BeckerError(
            f"{unparsed} match none of the four known naming schemes "
            f"(A001-C-007, B001-A-301, F007, CRC1_8810). The donor is the unit "
            f"of inference (invariant 5, Amendment 1), so a sample whose donor "
            f"cannot be read cannot be placed — and guessing attributes every "
            f"cell in it to the wrong person."
        )

    unknown = sorted(set(frame["disease_stage"]) - set(DISEASE_STAGE_MAP))
    if unknown:
        raise BeckerError(
            f"unknown disease stage(s) {unknown}. Known: "
            f"{sorted(DISEASE_STAGE_MAP)}. A value nobody has seen must stop "
            f"the run rather than be dropped — a dropped sample is a smaller "
            f"cohort reported as the same cohort."
        )
    frame["arm"] = frame["disease_stage"].map(DISEASE_STAGE_MAP)
    return frame


def sample_files(tar_path: str | Path) -> pd.DataFrame:
    """The tar's members, grouped into one row per sample.

    A sample needs all three of barcodes / features / matrix. An incomplete
    triplet is named and excluded rather than half-read.
    """
    kinds = {"barcodes": "barcodes", "features": "features", "matrix": "matrix"}
    rows: dict[str, dict] = {}
    with tarfile.open(str(tar_path), "r") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            name = Path(member.name).name
            match = re.match(
                r"^(?P<gsm>GSM\d+)_(?P<sample>.+?)_(?P<kind>barcodes|features|matrix)\.",
                name,
            )
            if not match:
                continue
            key = match.group("gsm")
            row = rows.setdefault(key, {"gsm": key, "sample_id": match.group("sample")})
            row[kinds[match.group("kind")]] = member.name

    frame = pd.DataFrame(list(rows.values()))
    if frame.empty:
        raise BeckerError(
            f"no GSM-prefixed 10x triplet in {tar_path}. Expected members like "
            f"'GSM6061645_A001-C-007_matrix.mtx.gz'."
        )
    complete = frame[["barcodes", "features", "matrix"]].notna().all(axis=1)
    frame["complete"] = complete
    return frame


def read_triplet(tar_path: str | Path, row: pd.Series):
    """One sample's matrix, as ``(counts, barcodes, features)``.

    ``counts`` is cells x genes CSR, which is the orientation everything in
    this repository expects. Matrix Market from CellRanger is genes x cells, so
    it is transposed here **once**, at the boundary, rather than by each caller.
    """
    from scipy.io import mmread
    from scipy.sparse import csr_matrix

    with tarfile.open(str(tar_path), "r") as tar:
        with gzip.open(tar.extractfile(row["matrix"]), "rb") as handle:
            matrix = mmread(handle)
        with gzip.open(tar.extractfile(row["barcodes"]), "rt") as handle:
            barcodes = [line.strip() for line in handle if line.strip()]
        with gzip.open(tar.extractfile(row["features"]), "rt") as handle:
            features = pd.read_csv(handle, sep="\t", header=None)

    counts = csr_matrix(matrix).T.tocsr()
    if counts.shape[0] != len(barcodes):
        raise BeckerError(
            f"{row['gsm']}: {counts.shape[0]} rows after transpose against "
            f"{len(barcodes)} barcodes. Matrix Market from CellRanger is genes "
            f"x cells; if this deposit is already cells x genes the transpose "
            f"here is wrong and every cell carries another cell's counts."
        )
    if counts.shape[1] != len(features):
        raise BeckerError(
            f"{row['gsm']}: {counts.shape[1]} columns against {len(features)} "
            f"features."
        )
    return counts, barcodes, features


def gene_symbols(features: pd.DataFrame) -> np.ndarray:
    """Symbols from a CellRanger features table.

    Column 0 is the Ensembl id and column 1 the symbol. **Reading column 0 and
    calling it a symbol is the identifier-space error this repository has made
    four times**, so the choice is explicit and the fallback is to raise rather
    than to use whatever is there.
    """
    if features.shape[1] < 2:
        raise BeckerError(
            f"features table has {features.shape[1]} column(s); CellRanger "
            f"writes at least (ensembl_id, symbol, type). Without the symbol "
            f"column the panel cannot be located, and column 0 is Ensembl — "
            f"using it would report every panel gene as absent."
        )
    return features.iloc[:, 1].astype(str).to_numpy()


def paired_donors(metadata: pd.DataFrame) -> pd.Index:
    """Donors carrying BOTH a tumour and a normal arm. The real cohort size.

    THE NUMBER BECKER AMENDMENT 2 IS ABOUT. 72 samples across the deposit, and
    four donors with both arms — the rest are polyps with no same-donor
    reference, or healthy donors with no polyps. ``healthy_donor`` does not
    count: pairing a polyp against a different person's normal tissue is
    cross-donor by construction, and Amendment 2 refuses that rescue.
    """
    paired = metadata[metadata["arm"].isin(PAIRED_ARMS)]
    counts = paired.groupby("donor")["arm"].nunique()
    return counts[counts == len(PAIRED_ARMS)].index


def pooling_key(metadata: pd.DataFrame, *, pool_by: str) -> pd.Series:
    """The identifier samples are pooled on. NO DEFAULT, on purpose.

    ``docs/prereg_becker_replication.md`` Amendment 1 fixed both options before
    this module existed:

    ``donor``   pooled per donor. **Primary and confirmatory** — it reproduces
                Chen_2021's one-arm-per-patient shape, and a replication that
                changes the estimand is not a replication.
    ``lesion``  one unit per sample. **Secondary and exploratory** — it keeps
                between-lesion variation, which is arguably the better design
                and is the wrong choice for a replication.

    Technical replicates pool under both: they are one physical sample
    sequenced twice, and ``sample_id`` already collapses them.
    """
    if pool_by == "donor":
        return metadata["donor"].astype(str)
    if pool_by == "lesion":
        return metadata["sample_id"].astype(str)
    raise BeckerError(
        f"pool_by={pool_by!r}; Amendment 1 fixed 'donor' (primary) and "
        f"'lesion' (secondary). There is no default because the two are "
        f"different estimands."
    )
