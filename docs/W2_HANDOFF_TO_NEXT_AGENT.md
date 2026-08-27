# W2 handoff — for whoever picks this up next

**Written 2026-08-26 by the outgoing W2 agent, running low on context.**
You are taking over **W2 (Method / harness)**. Read this, then
[CLAUDE.md](../CLAUDE.md), then [src/harness/README.md](../src/harness/README.md).

Everything below is either verified or explicitly flagged as unverified. Where I
was wrong, it says so — those entries are the most useful ones.

---

## 1 · Start here — the five-minute orientation

```bash
git pull
export PATH="$PATH:/c/Program Files/GitHub CLI"   # gh needs this on this machine
pytest -q                                          # ~3 min, expect ~909 passed
ruff check src tests
gh issue list --state open
gh pr list --state open
```

**Verify in a clean venv before claiming CI is green.** This has caught four
separate breaks and your local conda env will lie to you:

```bash
python -m venv /tmp/ci && /tmp/ci/bin/pip install -e ".[dev]" && /tmp/ci/bin/pytest -q
```

The project decomposes differentiation-marker loss in colorectal cancer into
**compositional** (mature cells gone), **cell-intrinsic** (cells present but
silenced), and **not estimable** (too few mature cells to ask). The third
segment is the contribution — every other method returns a number regardless.

W2 owns the **simulation harness**: the only place true ground truth exists, and
the stream that adjudicates the week-5 gate (G3 entirely; G2 and G4 read against
W2's cutpoints).

---

## 2 · What W2 has built, and what it is worth

All in `src/harness/`, all tested, all on `main`.

| Module | What it does | Trust level |
|---|---|---|
| `truth.py` | Analytic Kitagawa terms; parametric vs realised truth | solid |
| `pseudobulk.py` | Generator: patient holdout, composition draw, multiplicative shift on integer counts | solid |
| `interval.py` | **Within-patient** bootstrap CI — the only kind a positivity cutpoint can calibrate on | solid, read its docstring on invariant 5 |
| `calibration.py` | Derives cutpoints from pre-registered coverage/discrimination | works; numbers not yet meaningful (see §4) |
| `attenuation.py` | The §2.2 sweep, oracle + bulk arms | solid |
| `bulk_recovery.py` | The thing invariant 6 forbids *using*, measured | solid |
| `bakeoff.py` | Deconvolution ranking, signature-width comparison | solid |
| `controls.py` | Negative controls incl. the label-blind reference | solid |
| `rungs.py` | Does the estimator separate granularity rungs? | solid |
| `positivity.py` | The cutpoints. **W4 imports `classify_estimability` — do not fork it.** | solid |

### Results that are real and worth quoting

- **Bulk over-reports where the truth is undefined.** 41 of 80 `not_estimable`
  rows got a confident number, median |intrinsic| 103.5. Structural: bulk cannot
  count cells, so it cannot apply a positivity rule *or* form a per-patient
  interval. This is the project's thesis, measured, and it survives every caveat.
- **The estimator separates rungs when the truth says it should** — 33% relative
  gap on intrinsic, 45% on compositional; bit-identical (`0.00e+00`) when the
  partitions are the same. So W1's observed degeneracy is about the *labels*, not
  the estimator.
- **Granularity reallocates rather than changes the total.** Coarse
  (−6.58, −11.79), fine (−11.97, −7.87), total roughly conserved. That is §6.2's
  "divergence is the contribution", demonstrated where truth is known.
- **G4 flips on population choice**: matched-only (36) = 16.7% PASS; mixed (62) =
  51.6% FAIL. Decision #19 settled it as matched-only.

---

## 3 · Things I got wrong — read these first

**I published an acceptance check that could not fail.** I verified W1's four S
matrices for panel-gene leakage by intersecting `panel_genes()` (symbols) with an
S-matrix index that is **Ensembl IDs**. Two disjoint namespaces, so it returned
`none` every time and I reported PASS. GUCA2A, CDX2, SFRP1 and SFRP2 are in all
four; `best4` also has GUCA2B, CA7, OTOP2. **See issue #35 — it is open and it is
the highest-priority correctness item on W2's plate.**

The general lesson, which cost this project four separate bugs in a week: **a
guard that cannot fail is worse than no guard.** `assert_no_target_leakage` in
`src/reference/signature.py` is namespace-blind and `build_signature` calls it
four times, all vacuous.

**The permutation control was wrong twice** before it was right. v1 shuffled
labels before generating (silencing a random subset — a real effect, not a null).
v2 shuffled at estimation time and read 63%, which is *also* correct because
silencing 40% of cells moves the mean of any random subset. Only v3 — testing
against a **label-blind reference arm** — is a real control. If you are tempted
to assert a permutation control goes to zero, re-read `controls.py`'s docstring.

**Coverage was scored against the wrong truth** and came back a flat 1.0,
because the oracle estimate reproduces realised truth exactly, so I was asking
whether a percentile interval contains its own centre. It must be scored against
**parametric** truth.

---

## 4 · The bug family that keeps recurring

Four instances in one week, three workstreams, **all leaning toward the project's
hypothesis**:

1. **W4:** pooled MAD retention cut tumour epithelium 62% vs 88.5% normal —
   epithelium is 3.9× deeper than immune, immune sets the median, so the filter
   discarded epithelial cells for being epithelial. Inflates apparent
   compositional loss.
2. **W4:** label thresholds computed over all compartments; non-epithelial cells
   score as maximally mature on the inverted axis and drag the cut.
3. **W2:** `gate_g4_verdict` over a mixed matched/unmatched population — flips
   the gate.
4. **W1 (caught pre-emptively):** testing `unresolved_fraction` on 168 rows that
   are 28 patients counted six times.

**Shape: a cutoff or statistic computed over a mixed population, then applied to
a subgroup.** Never crashes; produces plausible wrong numbers.

**W2's code has been audited** — every other aggregation is computed *within* a
defined group. **W3 and W1 have not been audited by anyone.** That is an open
task and I think it is the highest-value unglamorous work left.

---

## 5 · What is left for W2, in priority order

| # | Task | State | Blocked by |
|---|---|---|---|
| 1 | **Fix the namespace-blind leakage guard** (#35) — resolve both sides through `config/gene_index/*.map.tsv`, raise on unresolvable IDs; add a regression test that puts an Ensembl ID in a symbol-indexed matrix | **not started, highest priority** | nothing |
| 2 | **Compositional-arm cutpoint** (#36) — I publicly pre-committed to `n_cells_resolved`: ok ≥50, wide 20–49, not_estimable <20, as decision **#21**. **Implement exactly that; do not re-choose the number.** | committed, not implemented | nothing |
| 3 | **Re-cost the gate at real n** — bootstrap CI width at n=10 (SMC) / 36 (GSE178341 matched) vs the assumed 60. I promised this to W4 twice. | not started | nothing |
| 4 | **Harden the `raw_counts` seam** — `LeeCohort.raw_counts` deliberately contains the panel, so it is one `reference_profiles()` call from an invariant-2 violation | not started | nothing |
| 5 | Real-data attenuation curve | not started | cell-level raw counts |
| 6 | Recalibrate cutpoints on a denser 5–50 grid | not started | 5 |
| 7 | Ambient-sensitivity sweep for G1 (CellBender cannot run — no empty droplets) | not started | nothing |
| 8 | Final gate memo | drafted, synthetic only | 2, 3, 5, 6 |

`docs/gate_memo_w2.md` exists and is **marked DRAFT ON SYNTHETIC DATA**. Keep
that marking until real cells are in it.

---

## 6 · Cross-workstream state

**W1 (Bode)** — very productive. Pilot run, four S matrices, inferCNV malignancy,
ambient across 62 patients. Open: **#14 blocks the headline** (neither labelling
method is a clean maturity measure); the epithelial rung's compositional term is
**structurally zero** (`mature_fraction` = 1.0 by construction — a denominator
choice); `lineage` ≡ `crypt_position` on `stem_pole` but **not** on
`opposite_lineage` (0.3711 vs 0.5326), so the granularity curve lives on axis 2.

**W3 (Jeremy)** — furthest along. Headline: **bulk GUCA2A loss is continuous, not
bimodal** (dip p 0.851/0.919/0.982), replicated in GSE39582 across a different
platform. This does *not* threaten the decomposition — ours is continuous-valued
already — it kills classifying patients from bulk.

**W4** — decisions #9/#10 resolved and merged. `doubly_robust` renamed to a
pooled-reference split; the interaction column now carries the real cross term,
so those three columns no longer sum to total (use `ADDITIVE_WEIGHTINGS` /
`identity_residual()`, do not hard-code). Deliberately produced **no decomposition
results** pending #14 — a defensible call the gate should know about.

**Cohort sizes are much smaller than the plan assumes.** GSE178341: 36 matched of
62 (~30 unsorted in both arms). SMC: 10 paired, not 23. Everything downstream was
designed for ~60. **Nobody has re-costed the gate.**

---

## 7 · Process facts that will save you time

- **`docs/open_decisions.md` has eight duplicate section numbers** (9–16 each
  appear twice) because three workstreams took the next free number
  independently. References by number are ambiguous. W2 offered a renumbering
  pass; W4 has no objection, W1 has not answered. Use ≥20 for new entries.
- **CODEOWNERS still has placeholder handles.** The one file every workstream
  edits (`[project.optional-dependencies]`) has broken `main` twice. Raised four
  times, still open as decision #5. It is the repo owner's to fix.
- **Any new third-party import must go in the `dev` extra.**
  `tests/test_dependencies.py` walks the AST and will catch you, including lazy
  imports inside functions.
- Never push to `main` without a clean-venv run. Branch names are `w2/…`.

---

## 8 · Open questions I did not get answers to

1. Should W4 emit decomposition numbers now against current definitions, with
   `quotable=False`, so the gate has something? (I said yes to the cheap version.)
2. Who audits W1's and W3's code for the §4 bug family?
3. Does the gate need formal re-costing at 10 / 36 patients — and who decides?
4. Is the four-resolution curve becoming a two-point finding an acceptable
   headline change? W2 and W4 both think it is publishable as a finding; nobody
   has decided it on purpose.
5. NeurIPS: **WMHS** (Atlanta, deadline **Sept 1**, 4pp) is the best fit — its CFP
   names *abstention* and *selective prediction*. Backup: **AI & Science**
   (Atlanta, Sept 7, 4/8pp), which lists "measuring properties that admit no
   ground truth". Nobody has decided whether to submit.
