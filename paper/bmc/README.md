# BMC Bioinformatics — methods-and-limitations submission

**Target:** BMC Bioinformatics, *Methodology* article.
**Governing decision:** [docs/DECISION_2026-09-11_scope_and_pivot.md](../../docs/DECISION_2026-09-11_scope_and_pivot.md).
Where this README and that document disagree, **the decision record wins** and
this file is the stale one.

**Branch:** `w1/bmc-manuscript`, cut from `w1/claim-check-audit` because every
number in this paper resolves against that branch's result tables.

## The claim, as fixed

Quoted from the decision record §3. **Do not paraphrase it into something
stronger while drafting:**

> In one colorectal single-cell analysis workflow, executable counterexamples
> expose validation checks that miss their stated failure conditions, while
> calibration and sensitivity analyses delimit which descriptive estimates
> remain reportable **without establishing cell-intrinsic silencing**.

The contribution is **the executable audit and its measured consequences**.
Neither "we found many bugs" nor "bootstrap intervals can under-cover" is
sufficient novelty alone.

**The adenoma example is subordinate to the validation question.** It is the
worked instance that makes the audit concrete, not the result being sold.

## What this paper may not say

Each of these is a claim the audit specifically retired. They are listed because
they are the sentences that will try to come back during drafting.

| forbidden | what is supportable instead |
|---|---|
| mature cells were silenced | a descriptive pattern in the selected, depth-matched populations |
| "71% of loss was caused by silencing" | mean absolute intrinsic **share** `\|i\|/(\|i\|+\|c\|)` is 0.715 / 0.709 — signs discarded, interaction excluded from the denominator |
| BCa is generally worse; Student-t is calibrated | a result **about the tested generators** |
| the bootstrap width expression is an exact coverage formula | a normal-approximation comparison for the sample mean |
| a zero residual means the validation is invalid | which assumption or performance property that validation cannot challenge |
| ρ is a variance decomposition | a ratio of median absolute departures; medians do not add |
| Crowell replicates the decomposition | supporting evidence for a tissue-level pattern, reproduced from derived inputs |
| the algebraic falsifier establishes identifiability | that *this particular degeneracy* is absent |
| GUCA2A is on no stock panel | none of the four audited panels supplies the target/control combination |

Statistics that are weighting- or denominator-dependent must be **quoted with
their weighting and denominator**. 6/8 cross-block at `lineage` holds under
`doubly_robust` and `normal`; it is 5/8 under `tumour`.

## Every number comes from a table

No figure is transcribed and no number is typed from memory. `_tables.py`
resolves a result table by name using **`newest_by_time`**, not `newest` —
the audit found 23 table names whose same-date runs disagree between the two
rules, 9 of them with different content, and it produced a wrong Crowell
comparison before it was caught. `paper/wmhs/_tables.py` still uses the
name-order rule; do not copy it.

The backbone is `results/*/reproduction_manifest.parquet` (newest by sidecar
time) — one row per deliverable with its producing command, sha, seed, row
count and a three-state `status`.

## Submission is gated on something not yet dischargeable

**The week-one stop rule is undischarged.** The decision record §4:

> If the blinded challenge finds **zero additional decision-relevant failures**
> beyond ordinary finite-value, identifier and input-validation checks — once
> legitimate estimator identities and estimand mismatches are excluded —
> **withdraw the methods submission.**

`guard_challenge_results` has `status=blocked` in the manifest, and
`blinded_guard_challenge` is committed but also `status=blocked`: every case has
`saw_ledger=true`, so `is_held_out_evaluation=false` and **the stop rule is not
readable from it.** The blinded challenger and rater B are unassigned.

The manifest's `status` is deliberately not a boolean. A table that exists but
cannot be used for its purpose is `blocked`, not `present`; `table_committed`
keeps the raw file fact separately. A `present=True` for the challenge would read
as success for the one thing the submission is gated on.

Drafting the manuscript is not the same as deciding to submit it. This paper
can be written to the point where only the challenge is missing; it **cannot be
submitted** until the rule is discharged or the decision record is amended in
writing.

## Two other obligations that are easy to forget

- **Workshop overlap must be disclosed.** `workshop_overlap.parquet` is the
  matrix. The WMHS source already carries adenoma non-collapse, the three-grid
  analysis, the interval material and the simulator experiments, so they
  **cannot be presented as untouched novelty**.
- **Reproduction paths are repo-relative.** `build_manifest` now records
  `results/<date>_<sha>/<table>.parquet`, never an absolute worktree path, so a
  fresh checkout of this branch resolves the index. (The `w1/claim-check-audit`
  branch still carries the older absolute-path producer; this fix should flow
  back there on the next merge.)

## Layout

```
main.tex          the manuscript; \input list only
sections/         one file per BMC section, each opening with its own
                  sourcing contract as a LaTeX comment
_tables.py        resolve a result table by name, by sidecar time
refs.bib          bibliography
make_fig1.py      Figure 1 — the interval stress surface
make_fig2.py      Figure 2 — cutpoint crossings as brackets
figures/          the two figure PDFs, committed and regenerated by the scripts
build.sh          builds and FAILS on any LaTeX error or undefined reference
                  (WMHS measured a page count on a 17-error build; never again)
reproduce.sh      clone into a temp dir with no raw data, rerun the
                  committed-input jobs, and check every output byte-for-byte
```

BMC ships a `bmcart` document class. It is **not vendored here** — download it
from the journal and drop it beside `main.tex`, the same arrangement WMHS uses
for the NeurIPS style. Until then `main.tex` builds under `article` so the prose
can be worked on, and **the built PDF is not submission geometry.**
