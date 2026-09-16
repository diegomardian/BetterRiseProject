# Pre-registration: prevalence audit of functional reuse in simulation-based validation of causal effect estimators

**Written:** 2026-09-13
**Author:** Diego Mardian (with Claude Opus 5)
**Branch:** `wmhs/prevalence-audit`
**Status at time of writing:** no repository, package, or supplement has been
opened, searched for, or read. The only material consulted is our own
`paper/wmhs/sections/full/blind.tex` and `.../refdesign.tex`. This document is
committed before the first search is run. Every subsequent change to it is a
dated amendment (§7).

---

## 1. The question

Our paper (§`sec:blind`) shows that a simulation-based validation can be
structurally vacuous. Let `D` be a draw from a simulator requested at parameter
`θ`. Let `T(D) = h(S(D))` be the validation reference computed on that draw, and
`θ̂(D)` the estimator under test. If

    θ̂(D) = h(S(D)) = T(D)

*identically*, then the residual `r(D) = θ̂(D) − T(D) ≡ 0` by algebra, the
"recovery curve" `θ̂(D)/θ` reduces to `T(D)/θ` — a property of the generator
alone — and the validation cannot fail.

The paper states plainly: *"Our oracle wiring is an observed feature of this
pipeline, not evidence about its prevalence in virtual-control-arm practice."*

**This audit asks whether anyone else's public code does the same thing.**

It is a prevalence question, not a criticism. The paper's own analysis
establishes that computing the reference from the generating parameters, or
from data the estimator cannot see, is the *correct* practice; finding that
everyone does it correctly is a publishable answer and we commit to publishing
it (§6).

### 1.1 What does NOT count as the failure mode

Fixed in advance, because this is where a sloppy audit manufactures a result:

- **Sharing a sample is not sharing a functional.** The estimator and the
  reference being computed from the same simulated dataset is a *necessary*
  condition, not a sufficient one.
- **Sharing summary statistics is not sharing a functional.** Different
  functions of the same summaries need not agree. `blind.tex` is explicit:
  *"Sharing summaries alone is insufficient."*
- **The canonical trap is our own OLS-vs-standardisation example.** OLS with
  stratum dummies and G-computation over the same strata both consistently
  estimate the same homogeneous effect from the same records, and both are
  weighted averages of the same stratum-level contrasts. They are *not* the
  same functional: standardisation weights strata by prevalence `n_g/n`, OLS by
  normalised `n_g p̂_g(1−p̂_g)`. Maximum residual `0.0663` at n=5,000, not zero.
  Any candidate YES that rests on "they both take a mean" is to be checked
  against this example before it is recorded.
- **"Both are unbiased for the same estimand" is irrelevant.** The failure mode
  is an algebraic identity on every draw, not agreement in expectation.

---

## 2. Search protocol

**Search date:** 2026-09-13. All searches run on this date; result ordering and
star counts are as of this date and are recorded per query.

### 2.1 Sources, in the order they will be searched

**S1 — GitHub repository search** (via `gh search repos`, falling back to
`https://github.com/search?type=repositories&q=...` through web fetch if the CLI
is unavailable or unauthenticated).

**S2 — GitHub code search** (via `gh search code`, or web fetch fallback), used
to reach validation scripts inside repositories whose description does not
advertise simulation.

**S3 — Web search for method papers with public code**, restricted to papers
that announce a code URL: queries run against a general web search engine, then
the announced repository is opened. This is the route to CRAN/Bioconductor
packages, Zenodo deposits, and journal supplements, which GitHub search does not
index uniformly.

We will not use private code, code behind a registration wall, or code we
cannot read in full.

### 2.2 Exact queries

Run verbatim. If a query returns nothing, that is recorded, not silently
replaced.

| id | source | query string |
|----|--------|--------------|
| Q1 | S1 | `synthetic control arm clinical trial simulation` |
| Q2 | S1 | `external control arm simulation` |
| Q3 | S1 | `target trial emulation simulation` |
| Q4 | S1 | `plasmode simulation` |
| Q5 | S1 | `causal inference estimator simulation benchmark` |
| Q6 | S2 | `"true_ate" simulate` |
| Q7 | S2 | `"true_effect" simulation estimator` |
| Q8 | S3 | `plasmode simulation causal effect estimator validation code github` |
| Q9 | S3 | `synthetic control arm virtual control arm simulation validation reproducible code` |
| Q10 | S3 | `target trial emulation simulation study R package github true effect` |

Additional queries may be added **only** as a dated amendment stating why, and
any repository reached by an added query is flagged as such in the results
table so a reader can discount it.

### 2.3 Inclusion criteria

A candidate is **eligible for examination** if ALL of the following hold, judged
from the repository itself (not from its README's claims alone):

- **I1 — Public code.** The full source is readable without registration, at a
  URL we can cite, with a resolvable commit or release version.
- **I2 — Simulation present.** The repository contains code that *generates*
  synthetic outcome data (a data-generating process it controls), rather than
  only resampling a fixed public dataset without generating outcomes.
- **I3 — Causal effect estimator under test.** The repository applies at least
  one estimator of a causal contrast (ATE, ATT, hazard ratio, risk difference,
  RMST difference, survival difference, or a decomposition term interpreted
  causally) to that simulated data.
- **I4 — Validation comparison present.** The simulated estimate is compared
  against some quantity described as truth, target, reference, or requested
  effect — as bias, error, recovery, coverage, RMSE, or a plotted
  estimate-vs-truth relationship. *Without a comparison there is nothing to
  audit.*
- **I5 — Domain.** The simulation is in trial simulation, synthetic/external
  control arms, target trial emulation, plasmode simulation, or a directly
  adjacent causal-estimator evaluation (e.g. a benchmark suite for ATE
  estimators). General-purpose ML benchmarks with no causal contrast are out.

**Exclusions, stated in advance:**

- **E1** — Our own repository (`BetterRiseProject`) and any fork of it. The
  paper already reports its own wiring; including it would be circular.
- **E2** — Repositories that are only a thin wrapper, teaching notebook, or
  homework exercise with no validation step (fails I4).
- **E3** — Repositories whose validation code is not in the repository (e.g.
  results are checked in but the producing script is absent). These are
  recorded as *screened, not examinable* with the reason, and they count
  toward coverage limits, not toward N.

### 2.4 Selection rule — fixed before looking

Cherry-picking is the obvious failure of an audit like this, so the order is
fixed now:

1. For each query Q1…Q10 **in numerical order**, retrieve results in the
   source's **default relevance order** and additionally record the star count.
2. Build a single **screening queue** by interleaving: take result #1 of Q1,
   #1 of Q2, … #1 of Q10, then #2 of Q1, and so on (round-robin), so no single
   query dominates. Duplicates are dropped at their first appearance.
3. Walk the queue **in order**. For each entry, apply I1–I5 and E1–E3 and record
   the screening verdict *for every entry walked*, eligible or not.
4. Examine eligible repositories in queue order until the target N (§5) is
   reached or the queue is exhausted.
5. **No entry may be skipped because it looks unpromising.** If an eligible
   repository is not examined, the reason is recorded in the results file.
6. If the queue is exhausted before N is reached, that is the finding: the
   population of publicly auditable examples is smaller than the audit's
   target, and we report N as achieved.

Deviation from this order is permitted only for a mechanical reason (the
repository is unreachable, exceeds what we can read, or is a duplicate), and
every deviation is logged with its reason in `docs/prevalence_audit_result.md`.

---

## 3. Classification rule

Written so a second person can apply it to the same repository without
consulting us. Applied **per repository**; if a repository contains several
independent validation experiments, the primary one (the one its README or
paper foregrounds) is classified, and any others are noted in the
justification.

### 3.0 Step 0 — Locate the comparison

Find the code that (a) produces a simulated dataset `D`, (b) applies the
estimator to `D` producing `θ̂`, and (c) compares `θ̂` to a reference `T`.
Record file and line for (b) and (c). If (c) does not exist, the repository
fails I4 and is *not classified* — it leaves the denominator.

### 3.1 Step 1 — Write down the estimator functional

State, in one line of mathematics or pseudocode, the exact map from the realised
records to `θ̂`: which rows, which columns, which weights, which conditioning
set, which nuisance models, which exclusions.

### 3.2 Step 2 — Write down the reference

Same, for `T`. Record where the values in `T` come from: a variable set by the
generator, a column of `D`, a separate function call, a stored file.

### 3.3 Step 3 — Provenance triage (this decides most cases)

| condition | verdict |
|---|---|
| **(a)** `T` is a function of the generator's *parameters* only (the requested effect `θ`, a coefficient in the DGP, a fixed constant) and does not depend on the realised draw | **NO** |
| **(b)** `T` is computed from records or variables the estimator does **not** consume — both potential outcomes `Y(0)`,`Y(1)` for every unit, latent uncensored event times, a held-out sample, an independent replicate, the full population from which `D` was drawn | **NO** |
| **(c)** `T` is computed from the same observed records the estimator consumes | go to **Step 4** |
| **(d)** the provenance of `T` cannot be determined from the repository | **UNCLEAR** |

Note on (b): a reference averaged over the *realised* sample's individual
causal effects (`mean(Y1_i − Y0_i)` over the drawn units) is a realised
quantity, but the estimator never sees `Y0` for treated units. It is therefore a
**different functional of a different input** and is **NO**. This is the
standard correct practice and we will say so.

### 3.4 Step 4 — Functional comparison

`T` and `θ̂` are **the same functional** iff `T(D) = θ̂(D)` holds *identically*
— for every draw in the simulation's support, by algebra, not in expectation
and not asymptotically. Check, in this order, and stop at the first difference:

- **4.1 Input records.** Same rows and columns? Same filtering, same handling of
  empty cells, same missing-data rule?
- **4.2 Weighting.** Same weights? Prevalence weights, variance weights,
  propensity weights and unweighted means are different functionals.
  (*Our OLS-vs-standardisation example fails here.*)
- **4.3 Conditioning set.** Same covariates/strata? Saturated versus parametric
  adjustment on the same covariates is a different functional.
- **4.4 Nuisance estimation.** Cross-fitting, sample splitting, regularisation,
  smoothing, or any tuning breaks the identity even when the target matches.
- **4.5 Observation process.** Same censoring/truncation/horizon handling?

**Verdict YES** requires *affirmative* evidence, one of:

- **Y-i** *Shared code path*: the reference and the estimate are produced by the
  same function, or the reference is a copy/alias of the estimate variable.
- **Y-ii** *Textual identity*: the reference expression is the estimator's
  expression applied to the realised draw, differing only in variable names.
- **Y-iii** *Derived identity*: a written algebraic argument that the two
  expressions coincide on every draw (the saturated-IPW-equals-standardisation
  case). The derivation must appear in the justification, with its conditions
  (e.g. "requires nonempty arms in every stratum").

**Verdict NO** if Step 3 gives (a) or (b), or if any of 4.1–4.5 identifies a
difference that makes `T` and `θ̂` differ on at least one realisable draw. A
single constructed counterexample draw is sufficient evidence.

**Verdict UNCLEAR** if the reference is computed from the realised observed
draw (3c) but the repository does not determine whether it is the same
functional — for example the nuisance model is fitted by an external library
call whose behaviour we cannot pin down, the reference is read from a data file
whose producing code is absent, or the relevant branch depends on configuration
we cannot resolve. For every UNCLEAR we state **exactly what additional
information would resolve it** (a named file, a config value, a run of the
code).

### 3.5 Optional confirmatory numerics

Where a repository runs cheaply and deterministically, we may execute its
validation and record `max |θ̂ − T|` over paired finite outputs, with the
tolerance stated. This is **confirmatory only**: a numerical zero triggers an
algebraic investigation, it does not by itself produce a YES, and a nonzero
residual does not by itself produce a NO (it may be floating-point-visible only
at some sample sizes). Execution is **not required** for any verdict, and
whether it was attempted is recorded per repository.

### 3.6 Reporting discipline

We describe only what the code does. We do not characterise the intent,
competence, or care of any author. Where a repository does the correct thing we
say so explicitly. Quotations of code are verbatim with file and line.

---

## 4. Recorded fields

For every **examined** repository:

`name`, `url`, `commit_or_version`, `domain`, `estimator_description`,
`estimator_location` (file:line), `reference_description`, `reference_location`
(file:line), `reference_provenance` (3a/3b/3c/3d), `step4_first_difference`,
`verdict` (YES/NO/UNCLEAR), `verdict_basis` (Y-i/Y-ii/Y-iii/3a/3b/4.1–4.5/…),
`code_executed` (yes/no), `max_residual` (if executed), `resolves_unclear`
(what would settle it), `queue_position`, `query_id`, `justification`.

For every **screened but not examined** entry: `name`, `url`, `query_id`,
`queue_position`, `screening_verdict`, `criterion_failed`.

---

## 5. Target N

- **Target: 5 repositories examined and classified.**
- **Floor: 3.** Below 3 the audit is reported as inconclusive on coverage
  grounds, with the queue walked recorded in full so a reader can see why.
- **More if cheap:** if additional eligible repositories in queue order can be
  classified without additional search, they are examined and included. The
  stopping point is recorded. We do not stop early because the counts look
  interesting in one direction.

---

## 6. What would make this a null result — and the commitment to report it

The audit is **null** if, among the repositories examined, **zero** receive a
YES verdict.

We commit, in advance, to reporting that outcome in full and in the paper's
own terms: *"We examined N public repositories meeting pre-registered criteria;
k were classifiable; 0 exhibited functional reuse."* We further commit to
reporting the natural strengthening of a null — that the examined repositories
computed their reference from generating parameters or from oracle/held-out
data, i.e. that they did the correct thing — as a positive statement about that
practice, not as a hedge.

A null result **does not** invalidate the paper's analysis, which is an
algebraic claim about what a recovery curve measures and is established
independently of how often the wiring occurs. It does bear directly on the
scope of the paper's implicit claim that the failure mode matters beyond one
lab, and we will say which way it cuts. Conversely, one or more YES verdicts
would support that scope claim, and we would report the count without
extrapolating a rate to the field: a queue of ten queries is not a random
sample of practice, and no prevalence percentage will be quoted.

Outcomes that are **not** reportable as findings, fixed now so they cannot be
adopted later: a YES obtained by relaxing §1.1; a count taken after dropping an
examined repository from the denominator; a prevalence rate.

---

## 7. Amendments

The classification rule (§3) and the selection rule (§2.4) may not be changed
after a repository has been seen, except as a **dated amendment in this file**
that states (i) the date, (ii) what changed, (iii) which repository prompted it,
(iv) why the original rule was inadequate. Any amendment is **re-applied to
every repository already examined**, and the results table records both the
pre-amendment and post-amendment verdict for those repositories.

Amendments are appended below. If this section is empty at the end of the
audit, no amendment was made.

### Amendment log

#### A-1 — 2026-09-13 — relaxed query strings after six of ten returned zero

**What changed.** Six of the ten pre-registered query strings (Q1, Q2, Q3, Q5,
Q6, Q7) returned **zero** results. This is not a statement about the field; it
is GitHub's search semantics. `gh search repos` matches multi-word queries as a
conjunction over repository *name, description and topics* only — it does not
index READMEs — so `synthetic control arm clinical trial simulation` requires
all six words in one short description. `gh search code` likewise ANDs terms, so
`"true_ate" simulate` requires both tokens in the same file.

Verified at the time of amendment: `gh search code '"true_ate"'` returns
results, `gh search code '"true_ate" simulate'` returns none.

**Prompted by.** No repository. The amendment is prompted by a mechanical
property of the search tool, observed before any repository's code was opened.
At the time of writing this amendment, the only repository-level information
seen was the name, description and star count returned by Q4 and Q10 in their
own result lists. No file inside any repository had been read.

**Why the original was inadequate.** The original strings were written as
natural-language topic descriptions rather than as index queries. Left
unamended, the audit would have screened from Q4 and the three web queries
alone, losing the synthetic-control-arm, external-control-arm and target-trial
domains entirely — the domains the paper's claim is actually about.

**The amended strings.** Each is the original with the non-indexing terms
removed; no new topic is introduced and no repository was selected into a
query.

| id | replaces | source | amended query string | n results |
|----|----------|--------|----------------------|-----------|
| Q1b | Q1 | S1 | `synthetic control arm` | 2 |
| Q2b | Q2 | S1 | `external control arm` | 2 |
| Q3b | Q3 | S1 | `target trial emulation` | 8 |
| Q5b | Q5 | S1 | `causal estimator simulation` | 0 |
| Q6b | Q6 | S2 | `"true_ate"` | 15 hits / 10 repos |
| Q7b | Q7 | S2 | `"true_effect"` | 15 hits / 8 repos |

Q4, Q8, Q9 and Q10 ran **as pre-registered** and are unamended. Q5b still
returns zero and is recorded as a zero-yield query, not replaced again.

**Re-application to repositories already examined.** None had been examined.
The amendment therefore requires no re-classification. Every queue entry
records the query that produced it, so a reader can discount anything reached
by an amended query (§2.2).

**What did not change.** The classification rule (§3), the inclusion criteria
(§2.3), the round-robin selection rule (§2.4), the target N (§5) and the
null-result commitment (§6) are untouched.

---

## 8. Known limitations, acknowledged in advance

- **Single rater.** One reader applies the rule; there is no second independent
  rater and therefore no inter-rater agreement statistic. The rule is written
  for reproduction by a second person, which is the mitigation, not a
  substitute.
- **Not a random sample.** GitHub relevance ordering is opaque and favours
  popular and well-described repositories. The queue is a convenience sample
  with a fixed order, which prevents cherry-picking but does not deliver
  representativeness. **No prevalence rate will be computed.**
- **Reading, not running.** Verdicts rest on code reading. Reading can miss a
  dynamic code path; §3.5 is optional and will not be available for most
  repositories.
- **Language coverage.** Whatever the queue delivers. If it is all R, or all
  Python, that is recorded as a coverage limit.
- **Publication and platform bias.** Code that is public, on GitHub, and
  English-described. Industry virtual-control-arm pipelines are largely not
  public, and that is the population the paper's claim is really about. This
  limit is stated in the result and is not recoverable by any amount of
  searching.
