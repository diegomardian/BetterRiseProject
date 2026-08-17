"""GDC query, download and the TCGA barcode. W3.

Two jobs, kept separate because only one of them needs the network:

1. **Query and download** STAR counts for TCGA-COAD and TCGA-READ.
2. **Parse the barcode** into the technical factors W3.4 tests for confounding —
   TSS, plate, sequencing centre, analyte. These are not in a metadata file
   anywhere; they are encoded in the barcode itself, and W3.4 cannot run without
   them.

Everything below the download boundary is a pure function over strings and
frames, so the barcode logic and the sample manifest are testable offline. That
matters more than usual here: the confounding tests in week 3 rest entirely on
these fields being parsed right, and a silent off-by-one in the plate field
would look exactly like "plate is not confounded with stage".
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import pandas as pd

GDC_API = "https://api.gdc.cancer.gov"

#: The two projects. Whether they are pooled or stratified in analysis is
#: docs/open_decisions.md #4 and is NOT decided here — both are ingested into
#: one matrix carrying a `project` column, which keeps that decision open.
PROJECTS = ("TCGA-COAD", "TCGA-READ")

#: TCGA sample-type codes. 01 and 11 are the brief's scope ("primary tumour and
#: normal-adjacent"); the rest are counted and excluded, never silently dropped.
SAMPLE_TYPES = {
    "01": "primary_tumour",
    "02": "recurrent_tumour",
    "05": "additional_new_primary",
    "06": "metastatic",
    "11": "normal_adjacent",
}
DEFAULT_SAMPLE_TYPES = ("01", "11")

#: TCGA-A6-2670-01A-01R-1410-07
BARCODE_RE = re.compile(
    r"^TCGA-(?P<tss>[A-Z0-9]{2})-(?P<participant>[A-Z0-9]{4})"
    r"-(?P<sample_type>\d{2})(?P<vial>[A-Z])"
    r"-(?P<portion>\d{2})(?P<analyte>[A-Z])"
    r"-(?P<plate>[A-Z0-9]{4})-(?P<centre>\d{2})$"
)


class GDCError(RuntimeError):
    """A GDC query or download did not return what was asked for."""


# ---------------------------------------------------------------------------
# The barcode
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Barcode:
    """A parsed TCGA aliquot barcode.

    ``patient_id`` is the TCGA-XX-XXXX prefix and is the join key to TCGA-CDR
    and to every clinical table. It is *not* the barcode — one patient can have
    several aliquots, which is what :func:`deduplicate_aliquots` is for.
    """

    barcode: str
    patient_id: str
    tss: str
    participant: str
    sample_type: str
    sample_type_name: str
    vial: str
    portion: str
    analyte: str
    plate: str
    centre: str

    @property
    def is_tumour(self) -> bool:
        return self.sample_type == "01"

    @property
    def is_normal_adjacent(self) -> bool:
        return self.sample_type == "11"


def parse_barcode(barcode: str) -> Barcode:
    """Parse a full aliquot barcode. Raises on anything that is not one.

    Strict by design. A barcode that does not match is either a different
    identifier type (sample-level, case-level) or a corrupted download, and
    both should stop the ingest rather than propagate a null TSS into the
    week-3 confounding tests.
    """
    match = BARCODE_RE.match(barcode.strip())
    if match is None:
        raise GDCError(
            f"{barcode!r} is not a full TCGA aliquot barcode. Expected the "
            f"28-character form TCGA-A6-2670-01A-01R-1410-07 — a shorter "
            f"identifier has no plate or centre, and W3.4 needs both."
        )
    parts = match.groupdict()
    return Barcode(
        barcode=barcode.strip(),
        patient_id=f"TCGA-{parts['tss']}-{parts['participant']}",
        sample_type_name=SAMPLE_TYPES.get(parts["sample_type"], "other"),
        **parts,
    )


def barcode_frame(barcodes: list[str]) -> pd.DataFrame:
    """Parse many barcodes into the sample manifest's technical columns."""
    return pd.DataFrame([asdict(parse_barcode(b)) for b in barcodes])


# ---------------------------------------------------------------------------
# Deduplication — one aliquot per (patient, sample type)
# ---------------------------------------------------------------------------


def deduplicate_aliquots(manifest: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Keep one aliquot per (patient, sample_type). Returns (kept, dropped).

    Some patients have several sequenced aliquots of the same sample. Leaving
    them all in double-counts those patients in every downstream model, and the
    patient is the unit of inference (CLAUDE.md invariant 5).

    **The rule, stated so it can be checked rather than trusted:** prefer the
    earliest vial (A before B), then the lowest plate, then the lexicographically
    smallest barcode. It is arbitrary but it is deterministic and it does not
    look at expression — a rule that picked, say, the highest-depth aliquot
    would let the data choose the sample, which is a subtle way to bias a cohort.
    """
    ordered = manifest.sort_values(
        ["patient_id", "sample_type", "vial", "plate", "barcode"], kind="stable"
    )
    keep_mask = ~ordered.duplicated(subset=["patient_id", "sample_type"], keep="first")
    return (
        ordered.loc[keep_mask].reset_index(drop=True),
        ordered.loc[~keep_mask].reset_index(drop=True),
    )


# ---------------------------------------------------------------------------
# Query
# ---------------------------------------------------------------------------


def _filters(projects: tuple[str, ...], sample_types: tuple[str, ...]) -> dict[str, Any]:
    """GDC filter payload for STAR gene-counts files on the given projects."""
    type_names = [SAMPLE_TYPES[t] for t in sample_types if t in SAMPLE_TYPES]
    gdc_type_names = {
        "primary_tumour": "Primary Tumor",
        "normal_adjacent": "Solid Tissue Normal",
        "recurrent_tumour": "Recurrent Tumor",
        "metastatic": "Metastatic",
        "additional_new_primary": "Additional - New Primary",
    }
    def term(field: str, value: list[str]) -> dict[str, Any]:
        return {"op": "in", "content": {"field": field, "value": value}}

    return {
        "op": "and",
        "content": [
            term("cases.project.project_id", list(projects)),
            term("data_category", ["Transcriptome Profiling"]),
            term("data_type", ["Gene Expression Quantification"]),
            term("analysis.workflow_type", ["STAR - Counts"]),
            term("cases.samples.sample_type", [gdc_type_names[t] for t in type_names]),
        ],
    }


def query_star_files(
    projects: tuple[str, ...] = PROJECTS,
    sample_types: tuple[str, ...] = DEFAULT_SAMPLE_TYPES,
    *,
    page_size: int = 1000,
) -> pd.DataFrame:
    """Ask the GDC which STAR-counts files exist. Returns file_id, name, barcode.

    Network call. Kept as the only one in the module so everything else can be
    tested offline.
    """
    import requests

    fields = ",".join(
        [
            "file_id",
            "file_name",
            "file_size",
            "md5sum",
            "cases.project.project_id",
            "cases.samples.portions.analytes.aliquots.submitter_id",
        ]
    )
    rows: list[dict[str, Any]] = []
    start = 0
    while True:
        response = requests.post(
            f"{GDC_API}/files",
            json={
                "filters": _filters(projects, sample_types),
                "fields": fields,
                "format": "JSON",
                "size": page_size,
                "from": start,
            },
            timeout=120,
        )
        response.raise_for_status()
        payload = response.json()["data"]
        hits = payload["hits"]
        rows.extend(_flatten_hit(h) for h in hits)
        pagination = payload.get("pagination", {})
        start += len(hits)
        if not hits or start >= pagination.get("total", 0):
            break

    if not rows:
        raise GDCError(
            f"GDC returned no STAR-counts files for {list(projects)}. Either the "
            f"filters drifted with an API change or the query is wrong — do not "
            f"proceed with a partial cohort."
        )
    return pd.DataFrame(rows)


def _flatten_hit(hit: dict[str, Any]) -> dict[str, Any]:
    case = (hit.get("cases") or [{}])[0]
    samples = case.get("samples") or [{}]
    portions = (samples[0].get("portions") or [{}])[0]
    analytes = (portions.get("analytes") or [{}])[0]
    aliquots = (analytes.get("aliquots") or [{}])[0]
    return {
        "file_id": hit.get("file_id"),
        "file_name": hit.get("file_name"),
        "file_size": hit.get("file_size"),
        "md5sum": hit.get("md5sum"),
        "project": (case.get("project") or {}).get("project_id"),
        "barcode": aliquots.get("submitter_id"),
    }


def download_file(file_id: str, destination: str | Path, *, chunk_bytes: int = 1 << 20) -> Path:
    """Fetch one file from the GDC data endpoint. Idempotent — skips if present."""
    import requests

    destination = Path(destination)
    if destination.exists() and destination.stat().st_size > 0:
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(f"{GDC_API}/data/{file_id}", stream=True, timeout=300) as response:
        response.raise_for_status()
        tmp = destination.with_suffix(destination.suffix + ".partial")
        with open(tmp, "wb") as handle:
            for chunk in response.iter_content(chunk_size=chunk_bytes):
                handle.write(chunk)
        tmp.replace(destination)
    return destination


# ---------------------------------------------------------------------------
# Sample manifest
# ---------------------------------------------------------------------------


def build_sample_manifest(files: pd.DataFrame) -> pd.DataFrame:
    """The W3.1 deliverable: barcode, project, sample type, TSS, plate, centre, analyte.

    Takes the query result and joins the parsed barcode onto it. Everything
    W3.4 needs to test technical confounding comes from here.
    """
    missing = files["barcode"].isna()
    if missing.any():
        raise GDCError(
            f"{int(missing.sum())} file(s) came back with no aliquot barcode. "
            f"Without it there is no patient id and no plate — investigate the "
            f"query rather than dropping the rows."
        )
    parsed = barcode_frame(files["barcode"].tolist())
    manifest = files.reset_index(drop=True).join(parsed.drop(columns=["barcode"]))
    return manifest.sort_values(["patient_id", "sample_type", "barcode"]).reset_index(drop=True)


def reconcile_counts(manifest: pd.DataFrame) -> pd.DataFrame:
    """Sample counts by project and sample type, for the portal reconciliation.

    W3.1's "done when" requires this to match the GDC portal. Print it, compare
    it by eye against the portal's own facet counts, and record the comparison
    in the note — a silent shortfall of forty samples is invisible in a matrix
    and fatal in a survival model.
    """
    table = (
        manifest.groupby(["project", "sample_type_name"], observed=True)
        .agg(n_samples=("barcode", "size"), n_patients=("patient_id", "nunique"))
        .reset_index()
    )
    return table.sort_values(["project", "sample_type_name"]).reset_index(drop=True)
