"""Emit the classified claim-check ledger as a versioned table.

Week 1, deliverable 1 of ``docs/DECISION_2026-09-11_scope_and_pivot.md`` §4.

    python -m src.reference.jobs.claim_check_inventory --no-write

Reads ``config/claim_check_ledger.yaml`` and nothing else. It touches no cohort,
no atlas and no result table, so it is laptop-runnable and outcome-blind by
construction. It classifies entries; it does not decide whether the underlying
scientific claims were true.

The reconciliation between the enumerated count and the count ``HANDOFF`` §3
asserts is reported in the sidecar and **is currently failing** (+2). That is
the intended state: the job refuses to write only if the mismatch is
undocumented.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.reference.claim_check import (
    LEDGER_PATH,
    agreement,
    assert_reconciles_or_says_why,
    build_inventory,
    load_ledger,
)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ledger", type=Path, default=LEDGER_PATH)
    parser.add_argument(
        "--rater-b",
        type=Path,
        default=None,
        help="a second rater's judgements, keyed by id. Absent means agreement "
        "is not established and the sidecar says so.",
    )
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    doc = load_ledger(args.ledger)
    rec = assert_reconciles_or_says_why(doc)
    inventory = build_inventory(doc)

    second = None
    if args.rater_b is not None:
        second = build_inventory(load_ledger(args.rater_b))
    pairs = agreement(inventory, second)

    by_class = inventory["classification"].value_counts().to_dict()
    by_verdict = inventory["verdict"].value_counts().to_dict()
    n_forcing = int(inventory["has_forcing_input"].sum())

    print(f"enumerated {rec.enumerated} entries, asserted {rec.asserted} ({rec.delta:+d})")
    print(f"  classes:  {by_class}")
    print(f"  verdicts: {by_verdict}")
    print(f"  with a committed forcing input: {n_forcing} of {len(inventory)}")
    if pairs is None:
        print("  inter-rater agreement: NOT ESTABLISHED (one rater)")
    else:
        print(
            "  classification agreement: "
            f"{int(pairs['classification_agrees'].sum())}/{len(pairs)}; "
            f"verdict agreement: {int(pairs['verdict_agrees'].sum())}/{len(pairs)}"
        )

    if args.no_write:
        print("--no-write: nothing written")
        return 0

    path = write_versioned_table(
        inventory,
        "claim_check_inventory",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        notes=(
            "The claim-check ledger enumerated and classified. Week 1 deliverable 1. "
            "Classification is a judgement by a single rater against criteria fixed "
            "in the ledger file before any entry was classified."
        ),
        extra_meta={
            "ledger_path": str(args.ledger.relative_to(args.ledger.parents[1])),
            "enumerated_count": rec.enumerated,
            "asserted_count": rec.asserted,
            "count_delta": rec.delta,
            "counts_reconcile": rec.reconciles,
            "reconciliation_note": rec.note,
            "by_classification": by_class,
            "by_verdict": by_verdict,
            "n_with_forcing_input": n_forcing,
            "n_raters": 1 if second is None else 2,
            "inter_rater_agreement_established": second is not None,
            "what_this_licenses": (
                "a described and classified inventory of the ledger's own entries"
            ),
            "what_this_does_not_license": (
                "a field-wide failure rate, a claim that the classification is "
                "reproducible between raters, or a count of independent training "
                "observations"
            ),
        },
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
