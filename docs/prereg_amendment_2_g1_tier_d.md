# Amendment 2 to G1 — tier D holds one gene, and a correlation needs more than one

**Written:** 2026-08-25 · **Author:** W1 (Bode) · **Status:** proposed, needs the
team · **Amends:** decision #17 in [`docs/open_decisions.md`](open_decisions.md)
and G1 in `execution_plan.md` §4

> **No G1 number exists.** `src/reference/checks.py` was written today and has
> never been run against expression — not on the pilot, not on the cohort. The
> defect below was found by reading the frozen panel, not by reading a result.
> **The credibility of this amendment is entirely in that timing**, which is why
> it is a dated document rather than an edit to #17.
>
> An amendment that changes a pre-registered statistic is the most suspect move
> available to an analyst. It is defensible here only because the original is
> *arithmetically impossible*, not because it gave an unwelcome answer.

---

## 1 · What was measured

The frozen panel, `config/panel.yaml`, unchanged since week 0:

| tier | role | genes | n |
|---|---|---|---|
| A | compositional control | GUCA2A, GUCA2B, OTOP2, CA7 | 4 |
| B | intrinsic control | MLH1, SFRP1, SFRP2 | 3 |
| **D** | **retained control (the negative control)** | **MS4A12** | **1** |

Decision #17 justifies its 0.2 separation threshold as *"the smallest separation
that survives n≈8 genes per tier."* **No tier has eight genes. Tier D has one.**

## 2 · Why this breaks G1 as specified

#17 specifies "Spearman correlation between gene abundance and apparent loss,
computed within each panel tier separately."

- **A Spearman correlation over one gene does not exist.** Not "is noisy" —
  is undefined. #17's threshold 1, *"|ρ| > 0.5 in tier D"*, cannot be evaluated.
- Threshold 2, *"the three tier correlations fall within 0.2 of one another"*,
  needs three correlations. Two of the three are unestimable.
- Tier A at n=4 is computable but barely: a rank correlation over four points
  takes only 11 distinct values (−1.0 to 1.0 in steps of 0.2), and with 24
  permutations its smallest attainable two-sided p is **0.083** — reached only
  by a *perfect* monotonic relationship. Tier A cannot reach significance at any
  conventional level, whatever the data says.

**Tier D is the half of the gate that carries the falsification logic.** Tier D
genes are chosen to have no differentiation story, so a strong abundance-loss
relationship *there* is the soup with no biological reading left. Losing tier D
does not weaken G1; it removes the part that could fail.

`checks.py` returns `not_estimable` rather than a number, and
`tests/test_reference_checks.py::test_frozen_panel_cannot_currently_evaluate_g1`
pins that behaviour, so the gate cannot be cleared by a tier nobody could
measure.

## 3 · What was considered and rejected

| | approach | why not |
|---|---|---|
| a | Add genes to tier D | The panel is frozen (CLAUDE.md invariant 3), and more to the point: **choosing a control set after discovering the old one is too small to fail is choosing a control set to pass.** Whatever genes were added, the test would no longer be pre-registered. |
| b | Use tier E as the negative control | Tier E is explicitly *"Excluded from the falsification rule"* with *"no pre-registered expectation"*, and its genes (AQP8, CA1, SLC26A3, …) are maturity markers **expected to be lost**. It is not a negative control; it is a second tier A. |
| c | Bootstrap or permute within tier D | Resampling one gene returns that gene. No amount of resampling creates a second data point. |
| d | Drop G1 | It is the check on decision #15's "measure and report rather than correct". Without it, the ambient limitation is asserted rather than tested. |
| e | **Move the correlation genome-wide and locate the tier genes on it** | Recommended. See §4. |

## 4 · Proposed replacement

**The insight that makes this work: with one gene you cannot compute a
correlation, but you can compute a percentile.** MS4A12's position among ~39,236
genes is precisely defined; its correlation with itself is not.

**G1a — the abundance–loss trend, genome-wide.** Spearman between gene abundance
(mean expression across the cohort) and apparent loss (Δ per-cell mean, tumour
minus normal) over **all genes on the shared index**, per study, never pooled
(CLAUDE.md invariant 4). This is the trend the soup produces. It is expected to
be strongly positive, and that is not a finding.

**G1b — where the panel sits on it.** Fit the trend, take each gene's residual,
and express each panel gene's residual as a **percentile of the genome-wide
residual distribution**. A gene at the 50th percentile behaves exactly like a
generic gene of its abundance. A gene at the 95th is losing far more than
abundance predicts.

The tier structure then makes a directional prediction that can fail:

- **Tier A** (compositional targets, genuinely lost) should sit **high** — loss
  in excess of what abundance explains.
- **Tier D** (MS4A12, retained control) should sit **at or below the middle** —
  it is the gene chosen to be *kept*.

### Pre-committed thresholds

**Committed 2026-08-25, before any residual has been computed.**

**G1 FAILS if any of:**

1. **Tier A's median residual percentile < 0.80.** If the compositional targets
   do not lose more than four in five genes of comparable abundance, the loss
   this project measures is abundance.
2. **MS4A12's residual percentile > 0.50.** The gene chosen to be retained
   showing above-median excess loss means the measurement is not specific to
   the biology.
3. **Tier A median − MS4A12 < 0.30 in percentile units.** The falsification rule
   of `config/panel.yaml` in its own terms: if the tiers do not separate, the
   estimator is broken and no biological claim may be made.

**G1 PASSES if all three clear.**

Tier B is **reported but not gate-bearing**. MLH1 is broadly expressed across
colonic epithelium, so its compositional term is structurally near zero (panel
tier B, `role`) and its residual carries a different meaning from tier A's.
n=3 could not support a threshold in any case.

### Why 0.80, 0.50 and 0.30

- **0.50 for tier D** is not a tuned number — it is the definition of "behaves
  like a typical gene". A retained control above the median is the failure.
- **0.80 for tier A** is a deliberately modest bar. These are the four genes the
  project's entire premise says are lost; requiring them to beat four in five
  comparable genes is the weakest form of that claim that still has content.
- **0.30 separation** is 0.80 − 0.50, so it carries no independent information
  when the first two hold; it exists to fail the case where both tiers drift up
  together, which is what a soup artefact would look like.

**None of these was chosen by looking at a G1 result, because none exists.**

## 5 · What is unchanged

- **#17's correction still stands.** Both statistics run: retention vs abundance
  (`execution_plan.md` §4, the named criterion) and loss vs abundance. This
  amendment changes the *unit of analysis* from within-tier to genome-wide plus
  tier position; it does not re-open which quantity is correlated.
- The `indeterminate` outcome in `g1_verdict()` is retained. #17's original rule
  leaves a gap and so may this one; surfacing it beats guessing.
- **A failure still means what #17 said it means.** Not that the project is
  wrong — that *this cohort cannot separate the signal from the soup*, and the
  honest report is the non-identifiability with its diagnostics.

## 6 · What the team must ratify

1. That the statistic moves genome-wide, and the reason is impossibility rather
   than inconvenience.
2. The three thresholds in §4, **before `checks.py` runs on real expression**.
   After it runs, the choice is unfalsifiable.
3. Whether tier B stays reported-but-not-gate-bearing.

Until 1–3 are ratified, `checks.py` returns `not_estimable` and G1 is undecided.
That is the correct state, not a blocker to route around.
