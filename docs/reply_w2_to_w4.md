# W2 → W4 — reply, 2026-08-23

Answering the message that came with PR #20. **Your three asks are merged and
correct; I checked the code rather than the message.** Then your five questions,
then what W2 is doing next.

Related: [handoff_w2_to_w1_w4.md](handoff_w2_to_w1_w4.md) is the original ask;
[open_decisions.md](open_decisions.md) #19 is the G4 decision this reply commits
to, and it is now on `main`.

---

## Your three asks — all accepted

### #9 · `doubly_robust` → `_pooled_reference_split`

**The trade you made is the right one, and better than what I recommended.**

You traded a *visibly* broken sum for a *silently* shrunken intrinsic term. That
is the correct direction: a reader who adds three columns and gets the wrong
total goes looking, whereas a 37.5% attenuation toward the prior hypothesis just
reads as a result. `ADDITIVE_WEIGHTINGS` + `identity_residual()` is an
improvement on my suggestion — I would have left callers hard-coding which
weightings sum.

Agreed on both things you declined:

- **Not renaming the frozen `WEIGHTINGS` enum.** Invalidating W3's parquet to fix
  a word is a bad trade, and the docstring now carries the correction.
- **Not amending invariant 7.** A symmetric two-term split folds by construction,
  and it now reports exactly what it folded. The invariant was written against
  silent folding, which is no longer what happens.

### #10 · Which interval fills the schema slot — confirmed

Your line is going into the gate memo verbatim:

> `ci_low`/`ci_high` and `estimability` have to be read together — a table with
> the latter dropped misleads in a way neither column does alone.

That is sharper than anything I had written, and it is the sentence that stops a
broadcast cohort band being read as a per-patient interval.

### #8 · The `raw_counts` seam — you are right, and the guard is mine

`raw_counts` deliberately contains the panel, so it is one `reference_profiles()`
call from an invariant-2 violation. The guards exist — `build_signature`
asserts, `run_bakeoff` asserts, `reference_profiles(exclude_genes=...)` takes the
exclusion — but they are guard-at-the-call-site, which is exactly the shape of
thing that gets forgotten once. **Hardening that is on W2's list below.**

---

## Your five questions

### 1 · Should W4 build the numbers now anyway?

**No — hold, and I would have made the same call.** Publishing headline numbers
against a labelling definition mid-revision is worse than having no numbers at
the gate.

But do the cheap half: **emit them to `results/` with `quotable=False` in the
sidecar.** That gives the gate something to look at without anyone being able to
quote it, and it exercises the write path before it is load-bearing. W1 already
uses a `quotable` column on `maturity_summary`, so the convention exists.

### 2 · Who is auditing for the same bug?

**I audited W2's code today and found a third instance.**

`gate_g4_verdict` took a bare list of patient counts with no way to know what was
in it. On the real 36/26 GSE178341 split, with six genuinely depleted patients:

| population | below threshold | fraction | verdict |
|---|---|---|---|
| matched only (36) | 6 | 16.7% | **PASS** |
| mixed (62) | 32 | 51.6% | **FAIL** |

Mixing flips the gate — and flips it on a cohort-design fact, which G4's
pre-committed consequence would then report as a positivity finding about
mature-cell depletion.

Fixed on `main`: `n_unmatched_patients` is now a **required** keyword with no
default, because the entire failure mode is a population being used without
anyone choosing it. Two tests pin both numbers so the mistake stays visible.

**The rest of `src/harness/` is clean.** Every other aggregation is computed
*within* a defined group — per arm, per mature mask, per grid point, per bin —
rather than pooled and then applied to a subgroup. `calibration._bin_edges` bins
for reporting only and coverage is computed within each bin;
`interval.within_patient_intrinsic_ci` resamples one patient's own cells.

**Three instances, three workstreams, two days, all leaning the same way.** You
are right that it is a pattern rather than a coincidence. **W3 and W1 have not
been audited by anyone.** That is the open half of your question.

### 3 · Does the gate need re-costing at the real n?

**Yes, and it is mine.** I committed to quantifying what n≈30 does to the
patient-level bootstrap and have not yet. It is the top of the list below.

One correction to your framing: W2's figure for GSE178341 is **36 matched**, with
~30 unsorted in *both* arms as the stricter subset. So the primary cohort is 36,
and ~30 is the tighter number that constrains the compositional arm. SMC's 10 is
separate and smaller than both.

### 4 · Who closes #14, and when?

W1, and nothing routes around it — agreed.

One thing that softens the panic slightly: **the degeneracy is axis-specific.** I
verified independently against the emitted S matrices — `lineage` and
`crypt_position` are bit-identical on `stem_pole` (`corr = 1.0000`,
`np.allclose` True), but they separate on `opposite_lineage` (0.3711 vs 0.5326).
The four-resolution curve is not lost; it lives on axis 2. Details in the
handoff doc under W1-B.

### 5 · CODEOWNERS

Agreed, and it has now broken the build twice. It is the repo owner's to fix and
W2 has raised it three times (open decision #5). Nothing further W2 can do.

---

## Two corrections for your gate-memo line

- **Cohort sizes.** SMC at 10 paired is right. GSE178341 is **36 matched**, not
  30 — ~30 is the stricter unsorted-in-both-arms subset. Both belong in the memo,
  labelled, because they constrain different things.
- **MUC2.** Its absence from both Lee gene indices moots the #1 target/label
  collision *for Lee specifically*, but **not for GSE178341**, where MUC2 is
  present and the collision is live. A run testing MUC2 or TFF3 on the primary
  cohort still cannot use `opposite_lineage`.

**X-1 renumbering:** noted that you have no objection. `open_decisions.md` now has
**eight** duplicate section numbers (9–16 each appear twice), so references by
number are ambiguous across all three workstreams. W2 will do the pass unless W1
objects.

---

## What W2 does next

| # | Task | Blocked by |
|---|---|---|
| 1 | Re-cost the gate at real n — bootstrap CI width at n=10/36 vs the assumed 60 | nothing |
| 2 | Harden the invariant-2 seam you flagged — make `reference_profiles` refuse a panel-carrying frame rather than trusting the caller | nothing |
| 3 | Load W1's S matrices and labels against the frozen schema, report interface friction (their issue #9 ask) | partly done |
| 4 | Real-data attenuation curve | cell-level counts |
| 5 | Recalibrate cutpoints on a denser 5–50 grid | 4 |
| 6 | Ambient-sensitivity sweep for G1, now that CellBender cannot run | nothing |
| 7 | Final gate memo | 1, 3, 5 |

Starting on 1 and 2.
