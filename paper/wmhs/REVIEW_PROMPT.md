# Review-and-revise prompt — WMHS @ NeurIPS 2026

Hand the model this file plus the repository (or `main.pdf`,
`sections/full/*.tex` and `sections/*.tex` if it cannot run a build). Everything
below the line is written to be pasted as-is.

---

You are reviewing, scoring and then revising a submission to the NeurIPS 2026
workshop **"World Models for High-Stakes Health: Reliable Clinical Trial
Simulation and Intervention-Aware Reasoning."** The source is in `paper/wmhs/`.

Do the three jobs in order, and **do not let the later ones contaminate the
earlier ones**: review it cold before you score it, and score it before you
decide what to change. Do not soften the review because you will have to fix
what you find.

## Venue facts that should shape your judgment

- **Non-archival, double-blind.** Deadline 15 September 2026; assume four days.
- **9 pages of main text** for the full paper, **4** for the extended abstract.
  References and appendices are excluded, but the main text must be
  self-contained.
- A **responsible-use statement covering limitations and impact is mandatory**;
  a submission without one is desk-rejected.
- Workshop bars are lower than the main conference's. The question is not "is
  this a NeurIPS main-track paper" but "is it correct, and will it generate a
  useful conversation in the room." Weigh novelty accordingly, and **say
  explicitly if you find yourself applying a main-track bar.**

## What the paper claims

That the standard way of validating a virtual control arm or in-silico trial
pipeline — simulate from a known truth, run the estimator, plot recovered
against true — can be structurally incapable of seeing the estimator. When the
estimator is a deterministic function of the same sufficient statistics the
generator drew from, it cancels *identically*, and the recovery ratio becomes
realised-truth over requested-truth: a property of the generator alone. The
claimed contribution is not the phenomenon but a one-line check (compare against
the **realised** draw; assert the difference is not identically zero), evidence
that estimators from different literatures fail it, a quantitative refinement
for the non-degenerate case, and the finding that the conclusion the blind curve
had certified does not survive a correct calibration.

---

# Job 1 — Review

**Correctness first.** The central claim is an algebraic identity, so check it
rather than admiring it. Does the identity follow under the stated condition? Is
the claimed equivalence between G-computation and saturated-propensity IPW
right? The Frisch-Waugh argument about OLS's variance-weighted estimand? The
censoring result — that the check holds against observed data and degrades
against latent event times? Flag anything overstated relative to the evidence
actually in the paper. **Being wrong is the only thing a paper on this subject
cannot survive.**

**Then novelty.** The paper preempts the obvious objection — that this is the
*inverse crime*, or Meng's *congeniality*, or the plasmode-benchmark warnings —
by arguing the degenerate case is a kind rather than a degree. Decide for
yourself whether that distinction carries a paper, and say so plainly.

**Then fit.** The paper's origin is single-cell genomics; the trial-simulator
experiment generalises it. Does the trial-simulation material carry enough
weight for this audience, or does the paper read as a different paper wearing
this workshop's hat? Name the page where a reviewer's topical judgment gets
made.

**Then what is missing.** Reach your own conclusions; these are only places
worth looking, not findings:
- What the paper claims about how *common* this failure is, and what that claim
  actually rests on.
- Whether existing instruments (e.g. simulation-based calibration) would catch
  this, and whether the paper says.
- Whether a reader could apply the check to their own pipeline on Monday, or
  whether the paper stops at the principle.
- What the absence of an artifact link costs at a double-blind workshop.

**Do not spend review on these.** They are settled, and flagging them tells the
authors nothing:
- Page limits, LaTeX build, and the mandatory responsible-use statement: both
  builds compile with zero errors and fit (full 8 of 9, short 4 of 4).
- "Only one instance / one gene / two cohorts" — stated repeatedly, in the
  limitations.
- "You make no biological claim" — deliberate and pre-committed; a falsification
  rule fired and the authors took it.
- Numbers disagreeing between the two builds — impossible, they share sources.

# Job 2 — Score

Give **a single integer 1-10 for how likely this is to be accepted at this
workshop as it stands.** Not how good it is. How likely it is to get in.

```
1-3   below the bar / desk-reject risk
4-5   plausible reject, depends on the reviewer draw
6-7   probably in
8-9   comfortably in, talk or spotlight candidate
10    best-paper contention
```

**State the calibration you used in one sentence** — what fraction of
submissions to a workshop like this you are assuming get accepted — because a
score without that is unreadable. Then give the three things most likely to sink
it, ranked, each with the objection in a reviewer's own voice, one sentence on
whether it is *fair*, and the smallest edit that would defuse it.

# Job 3 — Implement

Now make the changes. Rank candidate edits by reviewer-score gained per page of
space and hour of work, and implement in that order until you run out of either.

**Hard constraints. Violating any of these is worse than making no edit.**

1. **Never invent or adjust a number.** Every figure in the prose traces to a
   versioned table under `results/` with a commit sha and a fixed seed.
   `pytest tests/test_paper_numbers.py` checks prose against those tables — run
   it. If a number looks wrong, report it; do not correct it by guessing.
2. **No new experiments** unless you state plainly that the result is worth
   missing the deadline for.
3. **The builds now have distinct sources.** `main.tex` reads
   `sections/full/*.tex`; `main_short.tex` reads `sections/*.tex`. Review both,
   and require every shared numerical claim to agree. Do not assume an edit to
   one reaches the other.
4. **Space is not free.** The full build has about one page of slack. **The
   short build has none — it is at exactly 4 of 4.** Anything you add to a
   shared file must either be `\iffull`-gated or paid for with a cut.
5. **Do not modify `neurips_2024.sty` or `neurips_2026.sty`,** and add no
   geometry, spacing or float overrides. The style file says tweaking it risks
   desk rejection; meet the page budget by writing less.
6. **Do not cut the responsible-use statement** (`sections/responsible.tex`,
   rendered as section 5). Its absence is a desk rejection.
7. **Double-blind.** Introduce no identifying string, and no artifact URL unless
   you are also anonymising what it points at.
8. Figures are regenerated by `make_fig*.py`, never transcribed. If a caption
   and a figure would disagree after your edit, regenerate or leave it alone.

**Verify before you report.** Run `./build.sh` — it builds both and fails on
page overrun, on undefined references, and on LaTeX errors. **Read the `errors`
line before the page line**: a build with errors silently drops content and
under-reports its own length, which has already produced one bad page
measurement in this repo's history. Then run `./check_anonymity.sh`; it must
exit 0.

`build.sh` needs `neurips_2026.sty`, which is deliberately not committed. If it
is missing, the four-line recipe that stands in for it (the official 2024 style
behind a shim) is in `paper/wmhs/README.md` under "Two builds, one source."

---

# What to hand back

**A. The review.** Correctness, novelty, fit, what is missing — in that order.
Quote the paper when you criticise it. Distinguish what you *verified* from what
you are guessing. If you catch yourself writing a criticism the paper already
answers, go read that answer instead.

**B. The score,** with its calibration sentence and the three ranked risks.

**C. What you changed** — a table of file, what it said, what it says now, and
which reviewer objection it defuses. Note anything you decided *not* to do and
why, especially edits you judged too risky this close to the deadline.

**D. Anything factually or algebraically wrong,** quoted verbatim, with the
correction — and say whether you fixed it or only flagged it. **Leave this
section empty if you found nothing.** Do not manufacture an entry, and do not
pad any other section either; a short honest review is worth more than a long
one.

**E. The build result** — the literal `errors`, `undefined` and page lines for
both builds, and the exit status of `check_anonymity.sh`. If you could not
build, say so instead of assuming your edits fit.

**F. Your score again, after the edits,** with one sentence on what moved it. If
nothing moved it, say that.
