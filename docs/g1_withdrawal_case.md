# The case for withdrawing G1 as a gate criterion

**Written:** 2026-08-28 · **Author:** W2 (Method / harness) · **Status:** proposed,
needs the team · **Concerns:** `execution_plan.md` §5 G1, decision #17,
[prereg amendment 2](prereg_amendment_2_g1_tier_d.md),
[#37](https://github.com/diegomardian/BetterRiseProject/issues/37),
[#46](https://github.com/diegomardian/BetterRiseProject/issues/46)

> **No G1 number exists on real expression.** `src/reference/checks.py` still
> returns `not_estimable` and tier A has deliberately never been measured. Every
> claim below comes from simulation or from arithmetic. That ordering is the only
> thing that makes this document worth reading, and it is the same ordering that
> made amendment 2 worth reading.

---

## 1 · The proposal

**G1 ceases to be a gate criterion.** It is not replaced by a repaired version,
and its pre-committed consequence does not fire. The ambient question is answered
instead by the evidence in §4, which already exists and is stronger than G1 would
have been.

`execution_plan.md` §5 and the G1 row of the gate table are amended to say so, in
writing, with this document as the reason.

## 2 · Why G1 cannot be used as written

**G1 fails when the project's hypothesis is true.** W1 established this on #46 and
then, to their considerable credit, demolished their own proposed fix for it.

Threshold 2 asks MS4A12's within-abundance-bin percentile of M to be ≥ 0.50.
MS4A12 is **colonocyte-restricted** — `config/panel.yaml` says so in the same
line that calls it "frequently maintained". Tier A's genes are mature-restricted
too. So when mature colonocytes deplete, which *is* the central claim, every
mature-restricted gene's mean falls together, and both tiers sink relative to a
comparison set that is **not** cell-type-restricted.

Measured across three populations for M, on three worlds:

| construction | composition only | **composition + tier A silenced** | both silenced |
|---|---|---|---|
| *owed* | FAIL | **PASS** | FAIL |
| all-epithelium M | FAIL ✅ | **FAIL** ❌ | FAIL ✅ |
| within-mature M | FAIL ✅ | **FAIL** ❌ | FAIL ✅ |
| maturity-stratified M | FAIL ✅ | **FAIL** ❌ | FAIL ✅ |

The middle column is the world where the project is right and nothing else is
wrong. G1 fails it under every population tried.

**The defect is structural, not parametric.** MS4A12 is scored against
abundance-matched genes that do not share its cell-type restriction, so any
residual composition moves MS4A12 far more than it moves its comparison set. That
follows from the panel and the transcriptome, not from anyone's simulation
parameters.

## 3 · Why it cannot simply be repaired

**Threshold 2 is frozen by having been seen to fail.** A pre-registered threshold
that is adjusted after someone has watched it fail is no longer pre-registered,
whatever the adjustment is called. W1 stopped for exactly this reason after two
of their own constructions failed, and was right to.

**The same constraint binds W2, and binds it harder.** W2 is the workstream that
*ratified* amendment 2 on [PR #40](https://github.com/diegomardian/BetterRiseProject/pull/40).
A third construction proposed by the ratifier, after seeing two fail, would carry
the least credibility of anyone's.

**W2's ratification missed this, and the reason is instructive.** The five worlds
in `src/harness/g1_amendment.py` vary per-gene fold change against a flat
background. **None of them has a mature/immature compartment at all**, so in none
of them does a gene's mean over all epithelium differ from its mean within mature
cells. The ratification tested whether the statistic detects *abundance
dependence* — which it does, |ρ| = 0.94 on pure soup — and never tested whether it
distinguishes *depletion* from *silencing*, which is the question G1 exists to
answer. A harness that cannot represent the confound cannot rule it out, and
saying "ratified" over that gap was W2's error.

That is now fixed in the sense that matters: it is written down, and
`tests/test_g1_amendment.py` will be extended with a compartment-structured world
so the gap cannot reopen silently. It does not un-ratify anything, and it does not
license W2 to propose the replacement.

## 4 · What answers the ambient question instead

G1 asks: *is the intrinsic signal ambient contamination?* Three independent lines
of evidence already bear on that, and together they are stronger than a
percentile rule on one gene.

**4.1 · The harness has measured what residual ambient actually does**
(gate memo §10, `ambient_sensitivity_sweep.parquet`). At decision #16's 10%
exclusion cap — the worst contamination the cohort admits by design:

- real terms retain **94%** of their value;
- a compositional-only world acquires an apparent intrinsic term worth **4.6%**
  of its compositional one;
- and an intrinsic-only world acquires a compositional term of **exactly zero**,
  at every level tested.

The artefact is **one-directional and structurally so**: the compositional term
is a function of the mature-fraction *difference*, contamination moves means
rather than fractions, so **ambient can invent silencing and cannot invent
depletion**. That is a sharper statement than G1 was ever going to produce, it is
measured against known truth, and it bounds the damage rather than testing a
proxy for it.

**4.2 · The depth floor works, measured cohort-wide**
(decision #14, `depth_confound_reference.parquet`). W1 ran W2's
`depth_confound_report` unmodified over their own labels, 32 patients × 2 axes ×
4 rungs: median |ρ| between the maturity call and sequencing depth is **0.13**,
against **−0.92** for an unthinned labeller on Lee. Same instrument, and the
difference is the depth matching.

**4.3 · The plate-based route is untouched by any of this.** The ICBI atlas's
plate-based subset has essentially no soup, so an intrinsic signal surviving
there is evidence that is not contamination — and it does not depend on a
correction, a threshold, or a percentile.

**What none of this gives is a binary.** That is the honest cost of the proposal
and §6 states it.

## 5 · Why not the alternatives

| | option | why not |
|---|---|---|
| a | **Fire G1's consequence as it stands** | The consequence is *"paper becomes a caution about a widely-run analysis; pivot to snRNA-seq and spatial"* — abandoning the approach. Firing it on a criterion proven unable to pass when the hypothesis is true would be the single worst outcome available, and it would be *caused by the pre-registration* rather than prevented by it. |
| b | **Repair threshold 2** | §3. It is frozen by having failed. |
| c | **A fresh pre-registration** — new criterion and its test worlds specified before any measurement | Defensible, and the honest version of "keep an ambient gate". It costs time the project does not obviously have, and it needs someone who has *not* watched two constructions fail to write it. Offered as the alternative if the team wants a binary. |
| d | **Leave `checks.py` returning `not_estimable` and say nothing** | The gate then has no ambient criterion and nobody has said so out loud. That is the state today and it is the one thing worse than withdrawing deliberately. |

## 6 · What the project loses, stated plainly

- **G1 was the check on decision #15's "measure and report rather than correct".**
  Without it, the ambient limitation is characterised by §4 rather than tested by
  a rule with a pre-committed consequence. That is weaker in kind, and the paper
  should say so in its limitations rather than in a footnote.
- **A reviewer may reasonably ask why the ambient criterion was dropped.** The
  answer has to be this document, and it has to include that it was dropped
  *before* tier A was measured. If that ordering ever stops being true, the
  proposal should be refused.
- **Three of four gate criteria then rest on things other than a threshold.** G4
  is answered and resolvable (memo §12); G3 is preliminary on synthetic; G2 has
  never been run. The gate becomes less mechanical than intended, which is a real
  loss and not one to talk around.

## 7 · What the team is being asked to decide

1. **Does G1 stop being gate-bearing?** (§1) — or does the project want option
   (c), a fresh pre-registration written by someone unexposed to the failed
   constructions?
2. **If it stops, is §4 the accepted answer to the ambient question**, and does
   its weaker form get written into the paper's limitations now rather than at
   submission?
3. **Does `checks.py` keep returning `not_estimable`** — W2 says yes, and that it
   should keep saying so permanently rather than being deleted, because the file
   is the record that the criterion existed and was retired for cause.

Whatever is decided, it should be decided **before tier A is measured on real
expression**. After that, every option on this page becomes unfalsifiable.
