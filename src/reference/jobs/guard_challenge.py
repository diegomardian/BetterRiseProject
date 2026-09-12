"""Run the constructed challenge cases and emit the result table.

Week 1, deliverable 2 of ``docs/DECISION_2026-09-11_scope_and_pivot.md`` §4.

    python -m src.reference.jobs.guard_challenge --no-write

READ THIS BEFORE READING THE NUMBER
-----------------------------------
The table reports ``n_counting_toward_stop_rule``. **That figure is not the
week-1 stop rule and must not be read as it.** Every case in ``CASES`` was
written by someone who had already read the ledger, so a zero here means "the
harness reproduced known defects and correctly declined to count an identity",
not "a blinded search found nothing".

The sidecar carries ``is_held_out_evaluation: false`` for exactly this reason. A
number that resembles the decision statistic but is not it is the shape of the
defect this whole programme is about, and it is not going to be introduced here
by accident.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence
from pathlib import Path

from src.common.io import write_versioned_table
from src.common.provenance import DEFAULT_SEED
from src.reference.guard_challenge import CASES, run_challenges, stop_rule_summary


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, default=None)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--no-write", action="store_true")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)

    frame = run_challenges(CASES)
    summary = stop_rule_summary(frame)

    print(f"{summary['n_cases']} cases, {summary['n_clean_controls']} of them clean controls")
    print(f"  behaved as sealed:        {summary['n_as_expected']}/{summary['n_cases']}")
    print(f"  algebraic identities:     {summary['n_algebraic_identity']} (excluded by construction)")
    print(f"  input-validation only:    {summary['n_input_validation']} (excluded by construction)")
    print(f"  revised after falsification: {summary['n_revised_after_falsification']}")
    print(f"  counting toward stop rule: {summary['n_counting_toward_stop_rule']}")
    if not summary["any_case_blind_to_the_ledger"]:
        print(
            "  NOT A HELD-OUT EVALUATION: every case was constructed with the "
            "ledger in view. The stop rule cannot be read from this table."
        )

    if args.no_write:
        print("--no-write: nothing written")
        return 0

    path = write_versioned_table(
        frame,
        "blinded_guard_challenge",
        seed=args.seed,
        results_dir=args.results_dir,
        allow_dirty=args.allow_dirty,
        notes=(
            "Constructed inputs run against checks with the expected outcome sealed "
            "in the case definition. Week 1 deliverable 2. Development material: "
            "every case was written with the ledger in view."
        ),
        extra_meta={
            **summary,
            "is_held_out_evaluation": False,
            "stop_rule_readable_from_this_table": False,
            "why_not": (
                "every case has saw_ledger=true, so the set measures whether the "
                "harness discriminates, not whether a blinded search finds anything"
            ),
            "protocol_for_a_blinded_set": "docs/week1_challenge_protocol.md",
            "what_this_licenses": (
                "the claim that the harness separates a missed scientific failure "
                "from an algebraic identity and from ordinary input validation"
            ),
            "what_this_does_not_license": (
                "a failure rate, a field-wide prevalence, or the week-1 stop-rule "
                "decision"
            ),
        },
    )
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
