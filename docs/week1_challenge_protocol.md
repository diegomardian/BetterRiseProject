# Challenger protocol — constructing cases that count

**You are the third person in week 1 of
[DECISION_2026-09-11_scope_and_pivot.md](DECISION_2026-09-11_scope_and_pivot.md)
§4**, and your cases are the only ones the stop rule can be read from.

The twelve cases already in `src/reference/guard_challenge.py` do **not** count.
Every one carries `saw_ledger=true`: they were written by someone who had read
`docs/HANDOFF.md` §3 and the existing guard suite, so they demonstrate that the
harness discriminates and nothing more. That is the same objection the decision
record makes against counting the 53 existing guard tests as evidence, and it
applies to those twelve exactly as well.

## The one rule that makes your cases count

**Do not read `docs/HANDOFF.md` §3, `tests/test_checks_can_fail.py`,
`config/claim_check_ledger.yaml`, or `src/reference/guard_challenge.py` before
you finish constructing.**

Read instead:

- the checks themselves — `src/harness/`, `src/reference/`, `src/schema.py`;
- what each check's docstring says it detects;
- the pre-registrations in `docs/prereg_*.md`, which state what each analysis
  claims to be able to reject.

Your question for each check is: **what input would this check have to catch to
deserve its place, and does it?**

## Constructing a case

Each case is a `ChallengeCase` and every field is sealed before you run it.

```python
ChallengeCase(
    case_id="B01",
    target_check="src.harness.<module>.<function>",
    description="one paragraph: what the input is and why the check must catch it",
    expected="fires",              # or "passes" — SEAL THIS FIRST
    consequence="decision_relevant",
    constructed_by="challenger_b",
    saw_ledger=False,              # only true if you broke the rule above
    run=lambda: ...,               # returns True when the check fired
)
```

**`expected` is sealed before `run` is executed.** This is not ceremony. Rater
A's case C05 had its expectation falsified on first run because a prevalence was
chosen without computing the ceiling — and that error is now visible in the
table as `revised_after_falsification` with a required note. If you change an
input after seeing an outcome, you set those two fields. The harness refuses a
revision with no note.

**Roughly a third of your cases should be clean controls** — inputs where the
check *should* pass, and passing is correct. A challenge set made only of
violations cannot show the audit discriminates; it can only show it fires.

## Classifying the consequence, which is where the stop rule lives

- **`decision_relevant`** — a miss would change a scientific conclusion. Only
  these count.
- **`input_validation`** — a miss produces a crash or a malformed table, not a
  wrong conclusion. NaNs, dtypes, identifiers, empty frames. Ordinary software
  hygiene, and the stop rule excludes it by name.
- **`algebraic_identity`** — the check reports something mathematically
  necessary. Excluded by construction, because counting identities as findings
  is how an audit inflates its hit rate.

**Be strict with yourself on the first category.** The temptation is to classify
a crash as decision-relevant because finding it felt productive. Ask: if this
had shipped, would a number in a paper have been wrong, or would a job have
failed loudly? Only the first counts.

## Running, and reading the answer

```bash
python -m src.reference.jobs.guard_challenge --no-write   # dry run
python -m src.reference.jobs.guard_challenge              # writes the table
```

The stop-rule statistic is `n_counting_toward_stop_rule` computed over cases
with `saw_ledger=false`. The committed table currently reports it as `0` with
`is_held_out_evaluation: false` beside it — **that zero is not the stop rule and
must not be read as one.**

## The stop rule, stated so it can actually fire

> If the blinded challenge finds **zero additional decision-relevant failures**
> beyond ordinary finite-value, identifier and input-validation checks — once
> legitimate estimator identities and estimand mismatches are excluded —
> **withdraw the methods submission**, finish the workshop record and redirect.

Two things follow that are easy to get wrong:

1. **A case that behaves as sealed is not a finding.** Reproducing a defect you
   predicted is the harness working, not the audit discovering.
2. **Finding nothing is a real outcome and it is allowed to happen.** The rule
   exists because a longer taxonomy without demonstrated additional utility is
   not a paper. If your honest answer is zero, write zero.
