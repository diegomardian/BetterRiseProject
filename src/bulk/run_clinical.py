"""W3.5 — build the curated clinical table. CLI, idempotent.

    python -m src.bulk.run_clinical fetch
    python -m src.bulk.run_clinical build

Nothing is dropped except redacted cases. Endpoint exclusions are recorded as
per-endpoint flags so a patient missing DSS still contributes to PFI.
"""

from __future__ import annotations

import argparse
import sys

import pandas as pd

from src.bulk.cdr import (
    CDR_FILENAME,
    CDR_URL,
    ENDPOINTS,
    add_usability_flags,
    build_curated_table,
    cohort_reconciliation,
    event_summary,
    fetch_cdr,
    load_cdr,
    reconciliation,
    stage_disagreement,
)
from src.bulk.clinical import build_clinical_table, fetch_cases, fetch_msi
from src.bulk.gdc import read_manifest
from src.common.io import write_versioned_table
from src.common.paths import PROCESSED_DIR, RAW_DIR
from src.reference.ingest import append_manifest_row

SEED = 20260818

#: Set from --allow-dirty. Default False. These jobs used to pass
#: allow_dirty=True unconditionally, so the bulk arm could not write a
#: clean provenance stamp even from a spotless tree -- which is why every
#: committed bulk table records git_dirty: true.
ALLOW_DIRTY = False
CDR_DIR = RAW_DIR / "tcga_cdr"
GDC_CACHE = RAW_DIR / "gdc_clinical"
BULK = PROCESSED_DIR / "bulk"


def step_fetch(downloaded_by: str = "jeremy749") -> None:
    path = fetch_cdr(CDR_DIR)
    append_manifest_row(
        path,
        source_url=CDR_URL,
        accession="TCGA-CDR_Liu2018",
        downloaded_by=downloaded_by,
        workstream="W3",
        notes="TCGA Clinical Data Resource, Liu et al. 2018 Cell, Supplemental Table S1",
    )
    print(f"wrote {path} and a manifest row")


def step_build() -> None:
    path = CDR_DIR / CDR_FILENAME
    if not path.exists():
        raise SystemExit(f"{path} missing — run `fetch` first.")

    cdr = add_usability_flags(load_cdr(path))
    print(f"TCGA-CDR COAD/READ: {len(cdr)} patients")

    gdc_clinical = build_clinical_table(
        fetch_cases(cache=GDC_CACHE / "cases.json"),
        fetch_msi(cache=GDC_CACHE / "msi.json"),
    )
    curated = build_curated_table(cdr, gdc_clinical)

    manifest = read_manifest(BULK / "sample_manifest.tsv")
    expression_patients = sorted(
        manifest.loc[manifest["sample_type"] == "01", "patient_id"].unique()
    )

    print("\n=== cohort reconciliation ===")
    cohort = cohort_reconciliation(curated, expression_patients)
    print(cohort.to_string(index=False))

    print("\n=== per-endpoint exclusions (nothing dropped globally) ===")
    recon = reconciliation(cdr)
    with pd.option_context("display.width", 160):
        print(recon.to_string(index=False))

    print("\n=== events on usable patients ===")
    events = event_summary(cdr)
    with pd.option_context("display.width", 200, "display.max_colwidth", 60):
        print(events.to_string(index=False))

    print("\n=== covariate completeness, patients with expression ===")
    with_expr = curated.loc[curated["patient_id"].isin(expression_patients)]
    rows = []
    for column in ("stage", "age", "sex", "msi_status", "site", "treatment_outcome_first_course"):
        if column not in with_expr:
            continue
        present = int(with_expr[column].notna().sum())
        rows.append(
            {
                "covariate": column,
                "n_annotated": present,
                "n_patients": len(with_expr),
                "coverage": round(present / len(with_expr), 4) if len(with_expr) else None,
                "n_levels": int(with_expr[column].nunique(dropna=True)),
            }
        )
    completeness = pd.DataFrame(rows)
    print(completeness.to_string(index=False))

    disagree = stage_disagreement(curated)
    print(f"\nCDR/GDC stage disagreements: {len(disagree)}")
    if len(disagree):
        print(disagree.head(10).to_string(index=False))

    print("\n=== endpoint roles (invariant 9) ===")
    for endpoint, (role, note) in ENDPOINTS.items():
        print(f"  {endpoint:<4} {role:<10} {note}")

    BULK.mkdir(parents=True, exist_ok=True)
    curated.to_csv(BULK / "clinical_curated.tsv", sep="\t", index=False)
    for frame, name, note in (
        (recon, "tcga_cdr_exclusions", "W3.5 per-endpoint exclusions and reasons"),
        (events, "tcga_cdr_event_summary", "W3.5 events on usable patients"),
        (cohort, "tcga_cdr_cohort_reconciliation", "W3.5 CDR vs expression cohort"),
        (completeness, "tcga_covariate_completeness", "W3.5 covariate completeness"),
    ):
        write_versioned_table(frame, name=name, seed=SEED, notes=note, allow_dirty=ALLOW_DIRTY)
    print(f"\nwrote {BULK / 'clinical_curated.tsv'} and four results tables")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m src.bulk.run_clinical")
    sub = parser.add_subparsers(dest="command", required=True)
    f = sub.add_parser("fetch", help="download the TCGA-CDR workbook")
    f.add_argument("--downloaded-by", default="jeremy749")
    sub.add_parser("build", help="curate the clinical table")

    parser.add_argument(
        "--allow-dirty", action="store_true",
        help="write from a dirty tree; the recorded sha will not reproduce it",
    )
    args = parser.parse_args(argv)

    global ALLOW_DIRTY

    ALLOW_DIRTY = args.allow_dirty
    if args.command == "fetch":
        step_fetch(args.downloaded_by)
    else:
        step_build()
    return 0


if __name__ == "__main__":
    sys.exit(main())
