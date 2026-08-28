# Amendment 2 to G1 — the statistic cannot fail informatively, and tier D holds one gene

**Written:** 2026-08-25 · **Author:** W1 (Bode) · **Status:** **RATIFIED by W2
2026-08-27** ([#37](https://github.com/diegomardian/BetterRiseProject/issues/37)),
conditional on §4a below; open pending a second reader on W2's #40 ·
**Amends:** decision #17 in [`docs/open_decisions.md`](open_decisions.md)
and G1 in `execution_plan.md` §4

> **No G1 number exists.** `src/reference/checks.py` was written today and has
> never been run against expression — not on the pilot, not on the cohort. Both
> defects below were found by reading the frozen panel and by simulating the
> statistic on synthetic data, not by reading a result.
> **The credibility of this amendment is entirely in that timing**, which is why
> it is a dated document rather than an edit to #17.
>
> An amendment that changes a pre-registered statistic is the most suspect move
> available to an analyst. It is defensible here only because the original is
> *arithmetically broken* — in one case impossible to compute, in the other
> guaranteed to fail whatever the biology — not because it gave an unwelcome
> answer.

---

## 1 · Defect one: the statistic fails by construction

#17 specifies Spearman between **gene abundance** (mean expression across the
cohort) and **apparent loss** (Δ per-cell mean, tumour minus normal).

**An absolute difference scales with the thing it is differencing.** A gene
averaging 100 counts can lose 30; a gene averaging 0.1 cannot. So "loss"
measured as a raw Δ carries abundance inside it, and correlating the two
measures the units.

Simulated, 20,000 genes over four orders of magnitude, `seed=1`. **Every gene
loses exactly 30% in tumour** — uniform, proportional, no abundance dependence
anywhere in the truth:

| statistic | ρ | |
|---|---|---|
| #17 as written — cohort abundance vs absolute Δ | **−0.997** | artifact |
| normal-arm abundance vs absolute Δ | **−0.998** | artifact |
| **A vs log₂ fold change** (proposed) | **−0.038** | correct |

G1 fails at `|ρ| > 0.5` in tier D. **Under #17's statistic tier D returns
ρ ≈ −1 on data with no abundance dependence at all**, so G1 as pre-registered
cannot pass. A gate that always fails carries exactly as much information as one
that always passes — and had tier D held eight genes, this would have gone
unnoticed and been reported as a real result.

There is a second, subtler coupling. Measuring abundance "across the cohort"
folds the tumour arm into the predictor, and the tumour arm is where the loss
being tested lives. Measuring it on the normal arm instead is worse, not better:
loss is then bounded by its own predictor, since a gene with a normal-arm mean of
zero cannot show loss.

**The fix for both is the standard MA construction**, symmetric between arms and
scale-free:

- **M** = log₂(tumour mean / normal mean) — the apparent loss
- **A** = ½(log₂ tumour mean + log₂ normal mean) — the abundance

## 2 · Defect two: tier D holds one gene

The frozen panel, `config/panel.yaml`, unchanged since week 0:

| tier | role | genes | n |
|---|---|---|---|
| A | compositional control | GUCA2A, GUCA2B, OTOP2, CA7 | 4 |
| B | intrinsic control | MLH1, SFRP1, SFRP2 | 3 |
| **D** | **retained control (the negative control)** | **MS4A12** | **1** |

#17 justifies its 0.2 separation threshold as *"the smallest separation that
survives n≈8 genes per tier."* **No tier has eight genes. Tier D has one.**

- **A Spearman correlation over one gene does not exist.** Not "is noisy" — is
  undefined.
- **And it fails in the passing direction.** W2 sharpened this during
  ratification, and it is worse than "cannot be evaluated": `scipy.stats.spearmanr`
  on one observation returns **`nan` and does not raise**, and `abs(nan) > 0.5`
  is `False`. So #17's rule — *fail if |ρ| > 0.5 in tier D* — would not have
  errored. **It would have returned PASS**, for the one tier whose entire job is
  to be able to fail. Same defect as the leakage guard in issue #35, this time
  inside the gate. `checks.py` never reaches `spearmanr` for a tier below
  `MIN_GENES_PER_TIER`, and a test pins that.
- Threshold 2 needs three correlations; two of the three are unestimable.
- Tier A at n=4 is computable but barely: the statistic takes only 11 distinct
  values (−1.0 to 1.0 in steps of 0.2) and across 24 permutations its smallest
  attainable two-sided p is **0.083**, reached only by a *perfect* monotonic
  relationship. Tier A cannot reach significance at any conventional level,
  whatever the data says.

**Tier D is the half of the gate that carries the falsification logic.** Its
genes are chosen to have no differentiation story, so an abundance–loss
relationship *there* is the soup with no biological reading left. Losing tier D
does not weaken G1; it removes the part that could fail.

`checks.py` returns `not_estimable` rather than a number, and
`tests/test_reference_checks.py::test_frozen_panel_cannot_currently_evaluate_g1`
pins that, so the gate cannot be cleared by a tier nobody could measure.

## 3 · What was considered and rejected

| | approach | why not |
|---|---|---|
| a | Add genes to tier D | The panel is frozen (CLAUDE.md invariant 3), and more to the point: **choosing a control set after discovering the old one is too small to fail is choosing a control set to pass.** Whatever genes were added, the test would no longer be pre-registered. |
| b | Use tier E as the negative control | Tier E is explicitly *"Excluded from the falsification rule"* with *"no pre-registered expectation"*, and its genes (AQP8, CA1, SLC26A3, …) are maturity markers **expected to be lost**. It is not a negative control; it is a second tier A. |
| c | Bootstrap or permute within tier D | Resampling one gene returns that gene. No amount of resampling creates a second data point. |
| d | Keep the Δ statistic, regress abundance out of it | Regressing the predictor out of the response and then correlating the two is circular — and it still leaves tier D at n=1. |
| e | Drop G1 | It is the check on decision #15's "measure and report rather than correct". Without it, the ambient limitation is asserted rather than tested. |
| f | **Move to M-vs-A genome-wide, and locate the panel genes on it** | Recommended. See §4. |

## 4 · Proposed replacement

**The insight that makes tier D workable: with one gene you cannot compute a
correlation, but you can compute a percentile.** MS4A12's position among ~39,236
genes is precisely defined; its correlation with itself is not.

**G1a — the abundance–loss trend, genome-wide.** Spearman between **A** and **M**
as defined in §1, over all genes on the shared index, **per study, never pooled**
(CLAUDE.md invariant 4). This is the trend the soup produces. Restricted to genes
with a non-zero mean in **both** arms — a log ratio against zero is not a number,
and a pseudocount there invents the very quantity being measured.

**G1b — where the panel sits on it.** Genes are grouped into **20 bins of equal
size by A** (≈1,960 genes each), and each panel gene is reported as the
**percentile of its M within its own abundance bin**.

Binned rather than residuals-from-a-fit, deliberately: Spearman is a rank
correlation and produces **no fitted values**, so "residual from the trend" has
no definition until someone picks a functional form — and picking one is a free
parameter nobody pre-registered. A within-bin percentile assumes nothing about
the shape of the trend, and it is exactly what the thresholds below are already
phrased in terms of: *genes of comparable abundance*.

Verified on the §1 simulation: under uniform proportional loss the median
within-bin percentile is **0.500**, i.e. flat, as it must be.

The tier structure then makes a directional prediction that can fail. Percentiles
are of M, so **low = more lost**:

- **Tier A** (compositional targets, genuinely lost) should sit **low**.
- **Tier D** (MS4A12, retained control) should sit **near the middle** — it is
  the gene chosen to be *kept*.

### Pre-committed thresholds

**Committed 2026-08-25, before any M or A has been computed on real data.**

**G1 FAILS if any of:**

1. **Tier A's median within-bin percentile > 0.20.** If the compositional targets
   are not lost more than four in five genes of comparable abundance, the loss
   this project measures is abundance.
2. **MS4A12's within-bin percentile < 0.50.** The gene chosen to be retained
   sitting below the middle — more lost than a typical gene of its abundance —
   means the measurement is not specific to the biology.
3. **MS4A12 − tier A median < 0.30 in percentile units.** The falsification rule
   of `config/panel.yaml` in its own terms: if the tiers do not separate, the
   estimator is broken and no biological claim may be made.

**G1 PASSES if all three clear.**

Tier B is **reported but not gate-bearing**. MLH1 is broadly expressed across
colonic epithelium, so its compositional term is structurally near zero (panel
tier B, `role`) and its M carries a different meaning from tier A's. n=3 could
not support a threshold in any case.

### 4a · The premise G1's power rests on — W2's ratification condition

**Threshold 2 asks MS4A12's within-bin percentile to be ≥ 0.50, and a gene that
is unchanged against an unchanged background sits at 0.50 by definition.**

W2 measured the consequence across five worlds with known truth, 60 replicates
each (`src/harness/g1_amendment.py`):

| world | truth | P(PASS) | owed |
|---|---|---|---|
| broad loss, tier A hardest, MS4A12 kept | claim is true | **1.000** | PASS |
| **tier A gone, nothing else moves** | **also true** | **0.517** | PASS |
| every gene loses 30% | no biology | 0.017 | FAIL |
| loss is abundance and nothing else | pure soup | 0.017 | FAIL |
| MS4A12 lost as hard as tier A | tiers don't separate | 0.000 | FAIL |

Null false-pass rate 2.6%. Per threshold: tier A 5.1%, tier D **50.1%**,
separation 21.2%.

**So G1's power is not a property of G1.** It depends on there being a broad loss
background for the retained control to stand out against. Against one, MS4A12
sits at 0.890 ± 0.015 and P(PASS) = 1.000. Without one — a world where the
project's claim is true but *only* tier A moves — the gate is a coin flip.

That premise is the project's own and probably holds. **It was never stated**, and
it is stated here because if it fails, G1 rejects a true signal half the time for
reasons that have nothing to do with ambient RNA.

**A G1 FAIL must be read with this attached.** Amendment 2 made tier D
*computable*; it did not make n = 1 *powerful*. That is not a reason to reject it
— what it replaces returned ρ ≈ −1 on a null and silently passed tier D — but it
is a limit that belongs beside the thresholds rather than in a footnote.

**No numbers were changed during ratification.** Adjusting a pre-committed
threshold while ratifying it is the move this amendment exists to prevent; W2's
test pins 0.20 / 0.50 / 0.30 so neither side can drift them later.

### One property on the credit side

Within-bin percentiles are **invariant to library-size normalisation** — maximum
movement 0.016 under CP10K — because a global rescale of one arm shifts every M
by the same constant and leaves within-bin ranks untouched. A mean-based rule
would not have been. W2 found this; it was not claimed here, and it is a real
advantage of the rank construction.

### Why 0.20, 0.50 and 0.30

- **0.50 for tier D** is not a tuned number — it is the definition of "behaves
  like a typical gene of its abundance". A retained control below it is failure.
- **0.20 for tier A** is a deliberately modest bar. These are the four genes the
  project's entire premise says are lost; requiring them to beat four in five
  abundance-matched genes is the weakest form of that claim with any content.
- **0.30 separation** is 0.50 − 0.20, so it adds nothing when the first two hold.
  It exists to fail the case where both tiers drift together, which is what an
  ambient artefact would look like.

**None was chosen by looking at a G1 result, because none exists.**

## 5 · What is unchanged

- **#17's own correction still stands.** Both quantities are reported: retention
  vs abundance (`execution_plan.md` §4, the named criterion) and loss vs
  abundance. This amendment changes how loss and abundance are *measured* and the
  unit of analysis; it does not re-open which pair is compared.
- The `indeterminate` outcome in `g1_verdict()` is retained. #17's rule left a
  gap and so may this one; surfacing it beats guessing.
- **A failure still means what #17 said it means.** Not that the project is
  wrong — that *this cohort cannot separate the signal from the soup*, and the
  honest report is the non-identifiability with its diagnostics.

## 6 · What the team must ratify

1. That loss becomes a **log₂ fold change** and abundance the **MA average**.
   This is the substantive change, and §1 is the evidence for it.
2. That the unit of analysis moves genome-wide with within-bin percentiles, and
   the reason is impossibility rather than inconvenience.
3. The three thresholds in §4, **before `checks.py` runs on real expression**.
   After it runs, the choice is unfalsifiable.
4. Whether tier B stays reported-but-not-gate-bearing.

Until these are ratified, `checks.py` returns `not_estimable` and G1 is
undecided. That is the correct state, not a blocker to route around.

## Reproducing §1

```
python - <<'PY'
import numpy as np
from scipy.stats import spearmanr
rng = np.random.default_rng(1)
n_genes, n_cells = 20000, 3000
abundance = 10 ** rng.uniform(-2, 2, n_genes)
normal = rng.poisson(abundance[:, None] * np.ones(n_cells)).mean(axis=1)
tumour = rng.poisson(0.7 * abundance[:, None] * np.ones(n_cells)).mean(axis=1)
eps = 1e-3
print(spearmanr(0.5*(tumour+normal), tumour-normal).statistic)          # -0.997
print(spearmanr(0.5*(np.log2(tumour+eps)+np.log2(normal+eps)),
                np.log2((tumour+eps)/(normal+eps))).statistic)          # -0.038
PY
```
