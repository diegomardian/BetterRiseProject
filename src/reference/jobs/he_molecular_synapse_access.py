"""Probe candidate HTA11 WES VCF metadata after an authenticated Synapse login.

This calls ``Synapse.get(..., downloadFile=False)`` for the exact-biospecimen,
premalignant candidates emitted by ``he_molecular_gate``. It never downloads or
opens a VCF, selects an endpoint, or fits a model.

    pip install -e '.[a2]'
    synapse config
    python -m src.reference.jobs.he_molecular_synapse_access \
      --files /path/to/files.tsv \
      --biospecimens /path/to/biospecimens.tsv \
      --cases /path/to/cases.tsv \
      --no-write

Use ``--no-write`` first. A successful metadata read establishes only that the
authenticated account can resolve the entity; it does not establish VCF content
download, endpoint callability, image resolution, or model eligibility.
Run ``synapse config`` to save a default profile before this job; ``synapse
login`` only verifies supplied credentials.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.reference.a2_synapse import login_synapse
from src.reference.he_molecular_gate import build_inventory, candidate_synapse_entities


def _read(path: Path, *, name: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{name} metadata not found: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def candidate_entities(attrition: pd.DataFrame) -> pd.DataFrame:
    """Explode only exact premalignant candidate IDs; blank IDs are a refusal."""
    entities = candidate_synapse_entities(attrition)
    participants = attrition.loc[:, ["biospecimen_id", "participant_id"]].drop_duplicates()
    return entities.merge(
        participants, on="biospecimen_id", how="left", validate="many_to_one"
    ).loc[:, ["biospecimen_id", "participant_id", "entity_id"]]


def probe_metadata(syn: Any, candidates: pd.DataFrame) -> pd.DataFrame:
    """Resolve entity metadata only; an error stays evidence, not an absence."""
    rows: list[dict[str, str]] = []
    for row in candidates.itertuples(index=False):
        try:
            entity = syn.get(row.entity_id, downloadFile=False)
            properties = getattr(entity, "properties", {}) or {}
            rows.append(
                {
                    **row._asdict(),
                    "metadata_status": "readable",
                    "entity_name": str(properties.get("name", row.entity_id)),
                    "entity_type": str(properties.get("concreteType", "unknown")),
                    "error": "",
                }
            )
        except Exception as exc:  # pragma: no cover - authenticated remote state
            rows.append(
                {
                    **row._asdict(),
                    "metadata_status": "unresolved",
                    "entity_name": "",
                    "entity_type": "",
                    "error": str(exc),
                }
            )
    return pd.DataFrame(rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--files", required=True, type=Path)
    parser.add_argument("--biospecimens", required=True, type=Path)
    parser.add_argument("--cases", required=True, type=Path)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    attrition, _ = build_inventory(
        _read(args.files, name="Files"),
        _read(args.biospecimens, name="Biospecimen"),
        _read(args.cases, name="Case"),
    )
    candidates = candidate_entities(attrition)
    result = probe_metadata(login_synapse(), candidates)
    readable = result["metadata_status"].eq("readable").all()
    print(result.to_string(index=False))
    print(
        "SYNAPSE METADATA ONLY — no VCF content or image was downloaded. "
        f"candidate entities readable: {int(result['metadata_status'].eq('readable').sum())}/"
        f"{len(result)}"
    )
    if args.no_write:
        return 0

    path = write_versioned_table(
        result,
        "he_molecular_synapse_access",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        extra_meta={
            "prereg": "docs/he_molecular_data_gate.md",
            "probe": "Synapse.get(downloadFile=False)",
            "vcf_content_downloaded": False,
            "endpoint_selected": False,
            "all_candidate_metadata_readable": bool(readable),
        },
    )
    print(f"wrote {path}")
    return 0 if readable else 5


if __name__ == "__main__":
    sys.exit(main())
