# Week 1 — what was produced, what was found, and what is blocked

**2026-09-11. Branch `w1/claim-check-audit`, worked in a separate worktree at
`BetterRiseProject-audit` because the main checkout carries concurrent WMHS
edits and the result writer refuses a dirty tree.**

Week 1 of [DECISION_2026-09-11_scope_and_pivot.md](DECISION_2026-09-11_scope_and_pivot.md)
§4. Three deliverables were named; all three exist. Two of the four roles the
plan assigns to people are **not** discharged and cannot be by one worker — §5
says which and what unblocks them.

---

## 1. The three tables

| deliverable | table | what it says |
|---|---|---|
| 1 · classified ledger | `results/2026-09-11_*/claim_check_inventory.parquet` | 24 entries against an asserted 21 |
| 2 · challenge set | `results/2026-09-11_*/blinded_guard_challenge.parquet` | 12 cases, 6 clean controls, all behaving as sealed |
| 3 · overlap matrix | `results/2026-09-11_*/workshop_overlap.parquet` | 14 rows; novelty share of the overlap rows **0.20** |

## 2. The ledger's own count has no check — and that is the first finding

`docs/HANDOFF.md` §3 asserts in prose that a check unable to fail has been found
**twenty-one** times. Nothing enumerated them. The ordinals are scattered across
§3, §6g, §6k, §6m and §6n; §7 still says seventeen; the §3 table holds fifteen
rows while the prose around it describes more.

Transcribing every entry the repository describes gives **23**. It is recorded,
not resolved — resolving it means editing §3, which belongs to the write-up. Two
candidate sources are visible in the text and both are in the file's
`reconciliation_note`: the prose says *"the five found on 2026-09-05"* and then
lists **four** bullets, so one of that five is already a table row; and L23
carries no ordinal because §6n calls it "a sibling worth reporting beside" L22
rather than numbering it. Under both readings the count is 21 exactly — a
plausible reconstruction that no document states, so it is not adopted.

**The reconciliation guard refuses a silent mismatch and permits a documented
one.** That is the whole design: a count asserted in one place and enumerated
differently in another, with nothing comparing them, is the ledger's own subject
one level out.

### The classification

| class | n |
|---|---|
| `implementation_error` | 10 |
| `logical_impossibility` | 6 |
| `provenance_reporting_error` | 6 |
| `low_power` | 2 |

23 `scientific_failure`, 1 `disputed`, and **15 of 24 name a committed forcing
input** — nine do not, which is a measurement the workshop paper does not have.

**L01 is the disputed entry and it is the consequential one.** Whether a
recovery curve whose estimator cancels is a missed scientific failure or an
algebraic identity is what the WMHS paper's thesis rests on, and the 2026-09-11
audit reads it the other way. Recorded as unresolved. §4 below narrows it.

## 3. L24 — the first entry this audit found rather than transcribed

Found by running new code, not by review.

`write_versioned_table` validated `extra_meta` for reserved provenance keys
**after** `df.to_parquet`. A collision therefore raised with the table already
on disk and **no sidecar beside it**. Every provenance guard in the suite
iterates `*.meta.json` — including
`test_every_result_sidecar_names_a_commit_that_exists` — so an orphan parquet is
not checked and failed, it is **not checked**. Invariant 10 satisfied by
absence.

It is not hypothetical. The repo-wide scan added with the fix found **five
committed parquets with no sidecar**:

- four retracted `0.1.0-pilot` S matrices, kept on purpose beside
  `RETRACTED_s_matrices.md` — the retraction is the record;
- **`2026-08-17_29e8a04/tcga_sample_reconciliation.parquet`**, a genuine orphan
  in a directory where both its siblings carry sidecars, with nothing explaining
  why it does not.

**Fixed**: `extra_meta` is validated before anything reaches disk, and the
sidecar is now written *before* the parquet, so a half-written pair fails loudly
rather than silently. **This touches `src/common/io.py`, which is shared code —
the change is small and covered by two forcing inputs, and it is owed a review
this branch has not had.** `KNOWN_ORPHAN_PARQUETS` ratchets the five, and a
second test fails if an entry gains a sidecar and is left in the list.

## 4. The mathematical audit found four claims the paper owes edits on

From `config/workshop_overlap.yaml`, rows `M01`–`M04`. Three are contradicted.

**M01 — "An exact zero proves blindness"** (`blind.tex` L163). Overstated as
written. An exact zero establishes that *this comparison* is uninformative about
the estimator, not that the validation is: simulation can still measure bias,
variance and coverage against a population parameter the sample does not
reproduce exactly. **The abstract's conditional form is already correct** — "if
the estimator is a deterministic function of the sufficient statistics the
generator drew" — and `src/harness/calibration.py` already draws the distinction
in its coverage comment, explaining why it measures against parametric rather
than realised truth. So the defect is in one sentence of prose, not in the
method or the code.

**That materially narrows L01's dispute**, and it is the most useful thing week
1 produced for the write-up: the disagreement is not "is the thesis wrong" but
"does one sentence claim more than the algebra". Case pair **C03/C04** in the
challenge set demonstrates it — the same estimator, zero residual against
realised truth and real residuals against parametric truth.

**M02** — the bootstrap width ratio is a normal-approximation comparison for the
sample mean, not a closed form for coverage. `HANDOFF` §3a already calls it "a
floor", which is the right word and the one the paper should use.

**M03 — "88% of its recovery curve is generator noise."** ρ is a ratio of median
absolute deviations; `1 − ρ` is not a variance share. The sentence immediately
before it already states the correct reading, so the fix is to **delete the
percentage**, not to add machinery.

**M04** — the six-of-eight cross-block headline reproduces under `doubly_robust`
and `normal` weighting and gives **five of eight under `tumour`**. It must be
quoted with its weighting.

## 5. What one worker cannot discharge, and what unblocks it

The plan assigns week 1 to four people. Two of those roles are structurally
closed to a single worker, and pretending otherwise would reproduce the defect
this whole programme is about.

**Two independent classifiers.** `agreement()` returns `None` for one rater
rather than a table of perfect self-agreement, and the sidecar carries
`inter_rater_agreement_established: false`. A rater compared with themselves is
a check that cannot fail, and that would have been entry twenty-five.
Unblocked by [week1_rater_b_protocol.md](week1_rater_b_protocol.md).

**A blinded challenger.** All twelve committed cases carry `saw_ledger=true`.
The job prints, and the sidecar records, `is_held_out_evaluation: false` and
`stop_rule_readable_from_this_table: false`. **The `0` in
`n_counting_toward_stop_rule` is not the stop rule** — it means the harness
reproduced known defects and correctly declined to count an identity. A number
that resembles the decision statistic but is not it is exactly the shape of
defect being catalogued, so it is labelled everywhere it appears.
Unblocked by [week1_challenge_protocol.md](week1_challenge_protocol.md).

## 6. What the overlap matrix says about the submission

Of ten overlap rows: **5 already published, 3 extended, 2 new.** Novelty share
**0.20**. The two new rows are the classified inventory and the harness's
identity-versus-failure discrimination. The recovery-curve thesis, the
realised-draw diagnostic, the interval calibration grid, the
carcinoma-versus-adenoma non-collapse and the forcing-input rule are all the
workshop paper's, cited rather than claimed.

**This is the honest answer to reviewer objection 3 and it is not a comfortable
one.** A reviewer who finds undeclared overlap will not believe anything else in
the paper, so the matrix declares it rather than minimising it.

## 7. Verification

- 2,060 tests pass. 21 fail, all `ModuleNotFoundError` for `diptest`,
  `lifelines` and `anndata`, absent from the interpreter used
  (`/opt/anaconda3`), none touching this work.
- 45 new tests across three files, each validation carrying a forcing input.
- All three tables written from a clean tree: `git_dirty: false`, seeds and shas
  recorded.

## 8. Reading order for whoever picks this up

1. This file.
2. [week1_rater_b_protocol.md](week1_rater_b_protocol.md) and
   [week1_challenge_protocol.md](week1_challenge_protocol.md) — the two blocked
   roles. **Assign these before doing anything else in week 2**; both are
   spoiled by reading the artifacts first.
3. `config/claim_check_ledger.yaml` and `config/workshop_overlap.yaml` — the
   judgements, as reviewable data.
4. §4 above, for the four paper edits week 2 does not need to wait for.
