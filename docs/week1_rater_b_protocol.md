# Rater B protocol — classifying the claim-check ledger independently

**You are the second of the two classifiers in week 1 of
[DECISION_2026-09-11_scope_and_pivot.md](DECISION_2026-09-11_scope_and_pivot.md)
§4.** Rater A has already classified all 24 entries. Your job is to classify
them again **without seeing rater A's answers**, so that the agreement between
you means something.

If you read rater A's classifications first, stop — the measurement is gone and
no amount of care recovers it. There is no partial credit here.

## What you are given, and what you must not open

| open this | do not open |
|---|---|
| `config/claim_check_ledger.yaml`, fields `id`, `short_name`, `the_check`, `why_it_could_not_fail`, `source_anchor`, `found_by` | the same file's `classification`, `ambiguous_with`, `rationale`, `verdict` fields |
| `docs/HANDOFF.md` §3 and the sections each entry anchors to | `results/*/claim_check_inventory.parquet` |
| the code and tests each entry names | any commit message in this branch describing the classification |

The transcription fields are shared on purpose. **You are not re-deriving what
happened; you are judging what kind of failure it was.** Disagreeing about the
facts would be a different problem and would make the agreement statistic
uninterpretable.

Generate your view with:

```bash
python - <<'EOF'
import yaml
d = yaml.safe_load(open("config/claim_check_ledger.yaml"))
for e in d["entries"]:
    print(f"{e['id']}  {e['short_name']}")
    print(f"    check: {e['the_check']}")
    print(f"    why:   {e['why_it_could_not_fail']}")
    print(f"    found: {e['found_by']}  source: {e['source_anchor']}")
EOF
```

## The four classes, verbatim from the file

Use these and only these. They were fixed before any entry was classified.

- **`logical_impossibility`** — no dataset could have made the check fail. The
  compared quantity is bounded away from the threshold, the comparison is
  vacuous, or the criterion is unreachable by construction.
- **`low_power`** — the check could fail on some data, but not at the sample
  size, noise level or effect size it was actually run at.
- **`implementation_error`** — the code computes something other than what was
  specified. A correct specification, wrongly executed.
- **`provenance_reporting_error`** — two statements of one quantity were never
  compared, documentation and implementation diverged, or prose overrode an
  artifact that was correct.

**The hard cases are real and you should expect several.** An interval using a
normal quantile where a t quantile belongs is an implementation error whose
*consequence* is under-coverage — is that `implementation_error` or `low_power`?
There is no right answer supplied here. Record your primary class, put the
runner-up in `ambiguous_with`, and say why in `rationale`. Entries where the two
raters split on exactly these are the most informative rows in the table.

## The three verdicts

- **`scientific_failure`** — a check meant to be able to reject a scientific
  claim, that could not. Counts toward the stop rule.
- **`algebraic_identity`** — the "failure" is a mathematical identity behaving
  as identities do. Does **not** count.
- **`disputed`** — contested. Requires a `dispute_note`; the loader refuses a
  disputed verdict without one.

**L01 is the entry to think hardest about.** Whether a recovery curve whose
estimator cancels is a missed failure or an identity is the disagreement the
whole WMHS paper turns on. Form your own view before reading anyone else's.

## Producing your file

Write `config/claim_check_ledger_rater_b.yaml` with the **same structure** as
rater A's — the loader validates both identically — carrying every field.
Transcription fields must match rater A's exactly; only the four judgement
fields are yours. Then:

```bash
python -m src.reference.jobs.claim_check_inventory --rater-b config/claim_check_ledger_rater_b.yaml
```

The job prints classification and verdict agreement and writes both into the
sidecar. `inter_rater_agreement_established` flips to `true`, and the inventory
stops being one opinion with criteria attached.

## What to do about disagreement

**Do not reconcile silently.** Two raters who quietly converge have produced one
opinion with extra steps. Keep both files, report the agreement rate, and
discuss only the rows that split — in writing, in the paper, with the reasons on
both sides. A disagreement rate of zero on 24 judgement calls is evidence the
protocol leaked, not evidence the classes are clean.
