# Pre-registration — does the REF→lesion fall hold across patients?

**Written:** 2026-09-07 · **Author:** W1 (Bode) · **Status:** proposed ·
**Depends on** [prereg_crowell_feasibility.md](prereg_crowell_feasibility.md)
and its Amendments 1–2.

> **This is not blind and saying otherwise would be the defect it guards
> against.** Section `231` has been read. §2 states exactly what is already
> known from it. Everything below is fixed before any *further* section is
> downloaded, and §3 names the sections and the order in advance so the choice
> cannot be made after seeing which ones help.

---

## 1 · What this is, and why it needs its own document

`prereg_crowell_feasibility.md` §1 says it covers *"a feasibility read and
nothing else."* It answered that: on section `231`, GUCA2A clears the
usability bar in the reference domain and sits below it in the adenoma, the
per-cell reading is not licensed, and §6's pre-specified direction is observed
at **n = 1**.

Extending that to more sections is **not** more feasibility. It is an analysis
with a patient-level estimand, and it needs its own falsifiers, its own n, and
its own statement of what a null would mean. Running it under the feasibility
document would be scope drift of exactly the kind this project has recorded
twice this week.

**It is still NOT a replication of avenue A.** Different platform, different
estimand, and a maximum of seven patients against 43.

**It is still NOT a per-cell reading.** Prereg §7 stands: the negative probes
bound false positives, not false negatives. Nothing here licenses a
"GUCA2A-low = silenced" claim about any cell.

## 2 · What is already known, declared

From `results/2026-09-07_6e7e93e/`, section `231`, block `24H12439_A4`, n = 1
patient. All of it post-hoc relative to this document.

| gene | role | REF→TVA | REF→CRC |
|---|---|---|---|
| KRT8 | **control** | +0.642 | +0.181 |
| CDX2 | identity | +0.629 | +0.191 |
| EPCAM | epithelial | +0.609 | +0.513 |
| **MS4A12** | identity | **−0.729** | **−0.847** |
| **GUCA2A** | **target** | **−0.871** | **−0.854** |

Changes in `log_separation` on the detection scale. **CDX2 tracked the control
in both lesion domains while the control itself moved 3.5× between them; the
two targets went negative in both.**

**One patient, one control gene, no interval.** That is the observation this
document is designed to test, and it is the reason the document must fix its
rule before more data arrives.

## 3 · The sections, the order, and the disk protocol — fixed here

Eight `.h5ad` files exist. From `metadata.txt` (Zenodo `10.5281/zenodo.15550908`)
they cover **seven distinct blocks**, `231` and `232` being FOV ranges of one:

| sid | block (`pid`) | MB |
|---|---|---|
| 110 | B-2080151-01-18 | 1510.9 |
| 120 | B-1970164-01-04 | 2284.5 |
| 210 | B1891803_1-5 | 1129.1 |
| 221 | B1914093_B1-2 | 701.6 |
| 222 | 23H29642_A4 | 240.3 |
| 231 + 232 | 24H12439_A4 | 534.8 + 226.8 · **already read** |
| 242 | 23H42952_A19 | 1101.1 |

`241` (block `24H06944_A16`) is in `metadata.txt` and **has no `.h5ad` in the
deposit**. The paper describes eight colon sections plus one lymph node; which
sid is the node is not established, and §5's inclusion rule handles it without
needing to know.

**They do not all fit.** `/projectnb` had 8.61 GB free before Crowell and 7.85
GB after `231`+`232`; the remaining six total **6.97 GB**, leaving 0.88 GB. The
largest single file is 2.28 GB and fits comfortably.

**So the protocol is one at a time, in the order above, and the `.h5ad` is
deleted after its table is written.** Fixed here so that no section is skipped
because the ones already run looked good, and no section is added because they
did not. The table and sidecar are the durable artifact; the h5ad is
re-downloadable from a DOI with a recorded md5.

## 4 · The statistic — a difference in differences, on the detection scale

Per section, per gene, with `d` the `log_separation` from
`crowell_feasibility`:

    Delta(gene)  =  d(gene, lesion) - d(gene, reference)
    DiD(gene)    =  Delta(gene) - Delta(CONTROL)

**The control is `KRT8` and it is the only one.** ACTB is absent from this
deposit (feasibility Amendment 1 §3). EPCAM is `epithelial` in the frozen
`GENE_ROLES` and is **reported but never used as a control** — invariant 3
freezes the panel and the axes, and widening a control set after the fact to
make a band look sturdier is precisely the move this document exists to
prevent.

**Why a difference and not a band.** With one control gene a "band" is a point,
and a point has no width to sit inside. A difference is honest about that: it
subtracts the control's movement and leaves one number per gene per patient.

**Within-section noise is negligible and between-patient noise is the whole
question.** Detection is estimated over 60,000–120,000 cells per domain, so the
binomial standard error at p ≈ 0.36 is ≈ 0.002. Invariant 5 makes the patient
the unit; each qualifying section contributes **one** `DiD` per gene, and the
interval is a **Student-t over patients**. Not the percentile bootstrap — at
n = 7 it is 0.742× the correct width and fires 12.0% under a true null.

## 5 · Inclusion, fixed before any section is inspected

A section enters the analysis when it carries, in the domain column identified
by `--inspect`:

1. a **reference** domain with ≥ `MIN_CELLS_PER_DOMAIN` cells, and
2. at least one **lesion** domain (adenoma or carcinoma) with ≥ the same, and
3. a `KRT8` separation clearing `MIN_LOG_SEPARATION` in **both**.

**Rule 3 is the sensitivity gate and it is about the control, never the
target.** A section where GUCA2A is below the bar still enters — that is the
expected outcome in lesion domains and excluding it would select on the
outcome.

**Where a section carries both an adenoma and a carcinoma domain, the adenoma
is primary and the carcinoma is reported beside it.** Avenue A's estimand is
identifiable on adenoma; §2's carcinoma agreement is a bonus and is labelled
secondary in its own column.

**Blocks, not files, are the unit.** `231` and `232` are one block and
contribute one observation; where both carry the same domain the cells are
pooled before `log_separation` is computed, never averaged after.

## 6 · What n buys, before any of it is run

| qualifying patients | half-width `t/√n` | × avenue A's n=43 |
|---|---|---|
| 7 | 0.9248 | **3.01×** |
| 6 | 1.0494 | 3.41× |
| 5 | 1.2417 | 4.03× |
| 4 | 1.5912 | 5.17× |
| 3 | 2.4841 | **8.07×** |

**Below 3 no interval is reported at all** — `meta.MIN_STUDIES` is 3 and this
project has an n=2 arm it has repeatedly refused to take a verdict from.
**At 3 or 4 the interval is reported and pre-committed to be uninformative**;
it is there so a reader sees the width, not so a claim can rest on it.

## 7 · What would falsify it

The prediction is `DiD(GUCA2A) < 0`: the target falls relative to the control
from reference to lesion.

| outcome | reading |
|---|---|
| `DiD(GUCA2A)` interval excludes zero, negative | The fall is not one patient. **Still not silencing and still not per-cell** — §7 of the feasibility prereg is untouched. |
| `DiD(GUCA2A)` includes zero | No claim. At these widths (§6) this is weak evidence and **must not be quoted as a negative**. |
| `DiD(GUCA2A) > 0` in most patients | The `231` observation was a single-patient artefact. Report it and stop. |
| **`DiD(CDX2)` also negative and comparable** | **The effect is not gene-specific — it is a tier.** This is the discriminator, and it is the same one avenue A used: CDX2 sitting with the control is what separates "terminal differentiation falls" from "everything epithelial falls". |
| `DiD(KRT8)` non-zero | Impossible by construction — it is identically zero and is reported as an arithmetic check, not a result. |

**`DiD(CDX2)` is the row that matters most**, because a negative `DiD(GUCA2A)`
alone is consistent with the whole differentiation programme moving, which is
not the claim.

## 8 · What this cannot decide, unchanged

- **Per-cell silencing.** The false-negative rate is unmeasured (§7 there).
- **Whether the reference is field-affected.** No healthy-donor arm exists in
  this deposit. §6 there fixed the asymmetry: a fall is understated and
  therefore safe; no fall is uninterpretable. **That carries over verbatim.**
- **Avenue A's replication.** Different platform, different estimand, n ≤ 7.
- **Anything about the lymph node**, which is excluded by §5's rules without
  needing to be identified.

## 9 · Standing

**A patient-level directional test on public CC-BY data, pre-registered before
seven of the eight sections are downloaded, with the one already read declared
in §2.** It can show that a single-patient observation is or is not shared
across a cohort. It cannot make it a silencing claim, and §7's CDX2 row is what
stops it being quoted as one.

---

## Amendment 1 — 2026-09-07, after inspecting `110` and before its gate is run

§5 does not say what to do with **more than one lesion sub-domain of the same
class in one block**, and section `110` has three: `110_TVA1` (328,901 cells),
`110_TVA2` (81,948) and `110_TVA3` (57,207), beside `110_REF` (114,910) and
`110_CRC` (36,893). Fixed here before any `110` gate output exists.

### 1 · Same-class sub-domains are pooled at the cell level, one per block

§5 already states the rule for `231`/`232`: *"where both carry the same domain
the cells are pooled before `log_separation` is computed, never averaged
after."* Three adenoma regions in one block are the same case — **one patient,
one adenoma observation** — and invariant 5 makes the patient the unit
regardless of how many lesions were sampled.

**The consequence is stated rather than hidden: pooling weights by cell count,
so `110_TVA1` supplies 70.3% of that patient's adenoma estimate**, against
17.5% and 12.2%. Averaging the three would give each 33.3%. **Pooling is chosen
because §5 already chose it**, not because it is obviously better — an average
over three regions of one patient would be an unweighted mean of three
quantities with very different precisions, and neither choice is a patient-level
estimator of anything.

**The sub-domains are also scored separately and reported `exploratory`**, so a
reader can see whether the three lesions agree. They contribute nothing to `n`.

### 2 · Transcript-defined labels are NOT domain candidates. Named, because one is tempting.

`110` carries `lv1`, `lv2`, `jst` and `ist` — a full cell-type annotation
including **`epi.entero` (11,943 cells)**, `epi.entero-like_ACSF2`,
`epi.stem-like_LGR5` and `epi.fetal-like_MMP7`.

**`epi.entero` is a mature-enterocyte label built from transcripts, and GUCA2A
is part of that programme.** Using it to define "mature cells" and then
measuring GUCA2A inside it is the exact circularity that made the Becker mature
label's one positive arm unreadable (`prereg_becker_replication.md` RESULT), and
it is what the proposed amendment to invariant 2 is about.

**The domain column is `typ` and only `typ`.** It is histology. Feasibility
prereg §4 is the whole reason this deposit was chosen and it is not negotiable
per section.

### 3 · The "expected words" line in `--inspect` is a substring matcher

It reported `jst: ['ln']` and `lv2: ['lymph']` for section `110`. **There is no
lymph-node tissue in `110`** — those matched immune cell-type names such as
`DC.lymphoid`. The line reports where the pre-registration's vocabulary
*appears as a substring*, nothing more, and it must not be read as a domain
inventory. Recorded because a reader could take it for one.

### 4 · What Amendment 1 does not change

`n` is unchanged: `110` is one block and contributes one observation. The
statistic, the control, the inclusion rules, the width table and §7's falsifiers
are untouched. **No threshold moved.**

---

## Amendment 2 — 2026-09-07, after two blocks: the control is also epithelial

§4 takes the DiD against `KRT8` because ACTB is absent and it is the only
`control`-role gene in the deposit. **`KRT8` is an epithelial keratin, and it
does not behave like a capture control here.**

> **Corrected 2026-09-07.** The first version of this table gave `KRT8` as
> **+1.746** in block 231. That was **0.154186** — Becker's normal-arm KRT8
> detection — used in place of Crowell 231_REF's **0.360144**: a number from
> another dataset, in a scratch calculation, propagated into a
> pre-registration. Found by re-deriving the table from the committed parquet,
> which is where every value below now comes from. The correction changes the
> magnitude and not the direction.

`dlog(mu)` from reference to adenoma, per block. Every value is re-derivable
from the committed `crowell_feasibility.parquet` — `detection` and
`median_counts_per_cell` are columns in it.

| block | depth | `log(depth)` | KRT8 | CDX2 | EPCAM | MS4A12 | GUCA2A |
|---|---|---|---|---|---|---|---|
| 24H12439_A4 (231) | 1.43× | **+0.356** | +0.766 | +0.754 | +0.733 | −0.605 | −0.747 |
| B-2080151-01-18 (110) | 1.97× | **+0.676** | +0.941 | +1.156 | +1.322 | +0.182 | +0.102 |

**The whole epithelial group rises at or above the depth term and the two
targets do not.** KRT8, CDX2 and EPCAM come in at 2.1–2.2× the log-depth in
231 and 1.4–2.0× in 110; MS4A12 and GUCA2A **fall** in 231 and rise a fifth of
the depth term in 110. So `KRT8` is tracking epithelial *fraction*, which is
higher in a lesion than in mucosa, as well as capture — and so are `CDX2` and
`EPCAM`, the latter `epithelial` by role.

### What that costs, stated plainly

**The DiD does not separate compositional from intrinsic.** A negative
`DiD(GUCA2A)` says the target did not rise with the epithelial signal, and that
is equally consistent with:

- **compositional** — the mature cells that carry GUCA2A are a smaller share of
  the epithelium in the lesion, or
- **intrinsic** — mature cells are present and the gene is off in them.

Separating them needs a maturity label. Feasibility prereg §4 requires that
label be **histological**, and this deposit has no histological maturity
annotation — only lesion domains. The transcript-defined `lv2` labels
(`epi.entero` and the rest) are barred by multisection Amendment 1 §2 for the
circularity reason.

**So §7's falsifiers stand and §8 gains a line**: whichever way the DiD comes
out, this design cannot say which of the two mechanisms produced it. It can say
whether the target moves against its control across patients, and that is all.

**No threshold, statistic or falsifier changes.** This is a limitation being
written down where the claim is made, before an interval exists to over-read.

---

## Amendment 3 — 2026-09-07, at three blocks: §7 has no branch for what happened

At n = 3 the interval stops being refused, and the gene that excludes zero is
**MS4A12, not GUCA2A**. §7's branch table is written entirely around
`DiD(GUCA2A)` and `DiD(CDX2)`. It does not contain this outcome, and the
verdict function returns `NO CLAIM` — correct, because §7's second branch keys
on the target — while a different gene's interval sits clear of zero
underneath it. Recorded before anyone quotes that.

| gene | mean DiD | 95% CI | excludes zero |
|---|---|---|---|
| **MS4A12** | −0.951 | **[−1.855, −0.047]** | **yes** |
| GUCA2A | −0.966 | [−2.195, +0.262] | no |
| CDX2 | +0.265 | [−0.491, +1.021] | no |
| EPCAM | +0.112 | [−0.468, +0.692] | no |

### It is not a bigger effect. It is a tighter one.

The two means are the same to within 0.015. `MS4A12` clears because its
per-block spread is smaller — sd 0.364 against 0.495, so |mean|/sd is 2.61
against 1.95. At n = 3 and `t = 4.303` that is the whole difference.

### Why it must not be read as a positive result

**Four informative genes, no multiplicity control, and none was pre-specified.**
Under a global null at α = 0.05 with four genes, `P(at least one excludes
zero)` is **0.19**. §7 fixes falsifiers for two named genes and says nothing
about the family, so a single unanticipated gene clearing is inside what chance
supplies.

**n = 3 is the smallest this project reports anything at**, and §6 pre-committed
that at 3 or 4 the interval "is reported and pre-committed to be uninformative;
it is there so a reader sees the width, not so a claim can rest on it." The
width here is **8.07×** avenue A's.

**And MS4A12 is the secondary claim, not the primary.** In the adenoma work
GUCA2A's contrasts are the claim and MS4A12 failing costs the secondary
reading. A secondary marker clearing while the primary does not is not the
primary claim arriving by another route.

### What this changes

**Nothing in the design.** No threshold, statistic, gene or falsifier moves, and
`verdict()` continues to key on the target. §7 gains a row saying that a
non-target gene excluding zero is **not** a branch of it, and that reporting it
requires the three qualifiers above in the same breath.

**The honest summary at three blocks is `NO CLAIM`**, and MS4A12 is a line in
the table rather than a finding.

---

## Amendment 4 — 2026-09-07, written BEFORE the four-block aggregate is computed

Block `B1891803_1-5` (section 210) has `GUCA2A` **below the negative-probe
floor** in both lesion domains — detection 0.0118 against a floor of 0.0167 in
the adenoma, 0.0176 against 0.0236 in the carcinoma. `log_separation` is
negative there. That is new: in every previous block the target was below the
usability *bar* but still above the *floor*.

**This is being written down before the aggregate is run, because the aggregate
turns on it.** Three blocks give `DiD(GUCA2A)` = −0.966 [−2.196, +0.264],
including zero. Adding 210's −2.373 gives −1.318 [−2.609, −0.027], **excluding
zero** — and the block that flips it is the one whose target value is not a
measurement.

### 1 · The block is included. Excluding it would select on the outcome.

§5's rule 3 gates on the **control**, and `KRT8` clears the bar in both of
210's domains (3.121 and 4.077). §5 says in terms that a block where the target
is low still enters. Dropping 210 because its target is inconveniently low is
the move §5 exists to forbid, and the fact that dropping it changes the answer
is exactly why the rule was written in advance.

### 2 · But a below-floor observation is a BOUND, not a point

At detection below the false-positive floor, the observed signal is smaller than
what the noise process alone supplies. The true rate is consistent with **zero**,
and the DiD computed from it is limited by where the floor sits, not by the
biology. It **understates** the fall — the same conservative direction as §6 —
but its magnitude is not an effect size, and a mean over a mixture of bounds and
points estimates neither.

### 3 · The rule, fixed here

**Each block is flagged `target_below_floor`, and the summary carries
`n_blocks_below_floor`.** Then:

> **If the target's interval excludes zero, and does NOT exclude zero once
> blocks with `target_below_floor` are dropped, the result is
> `INDETERMINATE`** — not a positive. The converse holds too: an interval that
> excludes zero only *without* the floored blocks is equally indeterminate.

This is a leave-out sensitivity with a **pre-committed consequence**, applied
mechanically to whatever comes out. It is not an outcome-based exclusion: both
intervals are computed and reported, and neither is "the" answer when they
disagree.

**Why a consequence rather than a caveat.** `docs/prereg_becker_replication.md`
RESULT records a check that warned and did not bind, and it returned the most
optimistic verdict in its table. A caveat beside a number that excludes zero
gets dropped the first time the number is quoted.

### 4 · What this does not change

No threshold, statistic, gene or inclusion rule moves. `MIN_LOG_SEPARATION`,
the control, §7's falsifiers and §6's width table are untouched. Amendment 4
adds a **flag, a count, and one branch** to the verdict — and it is written
before the number it governs exists.

---

## Amendment 5 — 2026-09-07, before the five-block aggregate: at the floor is not below it

Block `B1914093_B1-2` (section 221) has **no** below-floor value for the target
— `GUCA2A` separates at +0.951 in the reference and +0.439 in the adenoma, both
positive, so Amendment 4's flag does not fire. But **`CDX2`'s reference
separation is +0.0033**: a signal **0.3% above the noise floor**. Its change is
+1.805, the largest of any gene in any block, and almost all of that is the
reference sitting on the floor rather than the adenoma rising.

**`CDX2` is §7's decisive row**, and a floor-inflated `CDX2` biases toward the
*favourable* verdict — "the discriminator does not fall" is easier to reach when
the discriminator's reference had nowhere to go but up.

### No threshold moves. The flag stays `< 0`.

Amendment 4's rule is a hard `log_separation < 0`, and **changing it to a margin
now, having seen that 221's CDX2 sits just above it, is the tuning this whole
document exists to prevent.** A rule that gets a tolerance added the first time a
value lands near its edge is not a rule.

### What is added is reporting, not a rule

Each block/gene carries `min_separation` — the worse of its two inputs — and the
summary carries `min_separation_any_block` and `n_blocks_within_0p1_of_floor`.
A reader can then see which DiDs rest on a value that is technically above the
floor and practically at it, and apply their own judgement, which is what a
number is for.

**Read `CDX2` at five blocks with this in front of it.** If `CDX2`'s interval
sits clear of zero on the strength of block 221, that is a statement about where
221's reference happened to sit.

---

## Amendment 6 — 2026-09-07, before the six-block aggregate: the rule was asymmetric

Block `23H29642_A4` (section 222) has **`CDX2` below the negative-probe floor in
its reference** — separation −0.113. Amendment 4's consequence is written only
for the **target**: *"If the target's interval excludes zero, and does NOT
exclude zero once blocks with `target_below_floor` are dropped…"*. Nothing
covers the discriminator.

**That is an asymmetry and it is not defensible.** §7's verdict rests on two
conclusions, not one — the target excluding zero *and* the discriminator not
excluding it. If either turns on a value that is a bound rather than a point,
the verdict turns on a bound.

### The rule, extended symmetrically

> **If `DiD(CDX2)`'s conclusion — excludes zero or does not — changes when the
> blocks where `CDX2` is below the floor are dropped, the result is
> `INDETERMINATE`**, on the same terms and with the same wording as Amendment 4.

### Declared: I can see which way this one points, and it is being written anyway

A below-floor **reference** inflates that block's DiD **upward**. `CDX2`'s DiD
in 222 is +0.995, its largest but one, and dropping it would pull `CDX2`'s mean
**down** — away from excluding zero, which is the direction §7's pass already
needs. So this rule is unlikely to fire at six blocks, and I am not pretending
otherwise.

It is written regardless, for two reasons. **A rule added only when it might
bite is the tuning this document exists to prevent** — the test of a symmetric
rule is that it goes in when it is inconvenient *and* when it is idle. And
block `242` is not yet read: `CDX2`'s floored blocks could look different at
seven.

### What does not change

No threshold, statistic, gene, inclusion rule or width table moves. Amendment 5
stands: the flag is a hard `< 0` and near-floor values are reported through
`min_separation`, not folded into the flag. This adds **one branch**, symmetric
to one that already exists, before the aggregate it governs is computed.

---

## Amendment 7 — 2026-09-07, before the seventh block: §7 has no third state

`DiD(CDX2)`'s lower bound has run **−0.251, −0.156, −0.020** across four, five
and six blocks. At six it is **0.02 from excluding zero**, in the direction
opposite to the target.

§7's table has two states for the discriminator: it clears zero *with* the
target (a tier moved), or it does not clear zero (the target's fall is
gene-specific). **It has no state for the discriminator clearing zero in the
opposite direction**, and one more block may produce it.

### The code asserted something that would have been false

`verdict()`'s final branch is reached whenever the discriminator does not clear
zero **in the same direction as the target** — which includes clearing it in the
other direction. Its message said *"CDX2 … does not exclude zero"* in both
cases. At seven blocks that could have printed a false statement into a
committed sidecar. Both branches now print the discriminator's interval.

### The third state, and what it does and does not mean

> **`TARGET FALLS AND THE DISCRIMINATOR MOVES THE OTHER WAY`** — the target's
> DiD excludes zero negative and the discriminator's excludes zero positive.

**What it strengthens.** §7 asks whether a tier moved. A discriminator moving
*against* the target answers that more firmly than one merely failing to clear
zero, and it is consistent with published work reporting CDX2 **gain** in
colorectal adenomas.

**What it does not strengthen — and this is the part that will be dropped if it
is not written here.** It is **not** a stronger claim about mechanism.
Amendment 2 is unchanged: `KRT8` is an epithelial keratin, and a rising `CDX2`
is equally consistent with **more CDX2-positive cells in the lesion** as with
anything intrinsic. A cleaner separation between two genes is not a separation
between composition and silencing.

**And it is not a new falsifier.** §7's prediction was and remains
`DiD(GUCA2A) < 0` with `DiD(CDX2)` not tracking it. This names a state that was
already inside "the discriminator does not track the target" and gives it its
own words, because a verdict that cannot say what it saw will be paraphrased by
someone who was not there.

Written before block `242` is read.

---

## RESULT — five of seven blocks, 2026-09-07

`results/2026-09-07_daa3cbd/`. Supersedes the four-block table
(`results/2026-09-07_a9b3d6d/`), which is kept and must not be quoted.

### Verdict: TARGET FALLS AND THE DISCRIMINATOR DOES NOT — §7's pass branch

| gene | per-block DiD | mean | 95% CI | |
|---|---|---|---|---|
| **GUCA2A** | −1.513, −0.838, −0.548, −2.373, −1.417 | **−1.338** | **[−2.212, −0.464]** | **excludes zero** |
| **MS4A12** | −1.371, −0.759, −0.724, −2.687, −1.256 | **−1.359** | **[−2.348, −0.370]** | **excludes zero** |
| CDX2 | −0.012, +0.215, +0.591, −0.002, +0.899 | +0.338 | [−0.156, +0.832] | includes |
| EPCAM | −0.033, +0.381, −0.012, +0.052, +0.387 | +0.155 | [−0.108, +0.417] | includes |

**GUCA2A and MS4A12 are negative in all five blocks; CDX2 and EPCAM are not
consistently signed and neither clears zero.** That is §7's pass: the target
falls and the discriminator does not, so the effect is not simply a tier moving.

### Every guard was checked, and none of them was bypassed

**Amendment 4's leave-out holds.** GUCA2A excludes zero with the below-floor
block (`B1891803_1-5`) at [−2.212, −0.464] **and without it** at
[−1.815, −0.343]. The sign does not turn on a bound, so the INDETERMINATE branch
correctly did not fire. At four blocks it did.

**Amendment 5's near-floor block does not manufacture it.** `B1914093_B1-2`'s
CDX2 reference sat **+0.003** above the noise floor, which inflates its CDX2 DiD
upward and therefore *helps* "the discriminator does not fall". Dropping that
block entirely: CDX2 **+0.198 [−0.251, +0.647]**, still including zero; GUCA2A
**−1.318 [−2.609, −0.027]**, still excluding it. **The verdict survives removing
the block that biased toward it.**

**`n_blocks_within_0p1_of_floor` is 1 for CDX2 and 1 for GUCA2A**, and both are
named above rather than left in a column.

### What it is not — and this is most of what there is to say

**Not silencing, and not per-cell.** Feasibility §7 is untouched: the negative
probes bound false positives, not false negatives, so no block licenses a
"GUCA2A-low = silenced" claim about any cell.

**Not a separation of compositional from intrinsic.** Amendment 2 caps this at
any n: `KRT8` is the only control-role gene and it is an epithelial keratin, so
`DiD(GUCA2A) < 0` says *the target did not rise with the epithelial signal*.
That is equally consistent with the mature cells being a smaller share of the
lesion's epithelium and with them being present and silenced.

**Not endorsed by §6 either.** §6 pre-commits n=3 and n=4 to be uninformative.
It says nothing about n=5, which is **not** the same as calling it informative —
the width is still **4.03×** avenue A's.

**No multiplicity control, and none was pre-specified.** Two of four informative
genes clear zero in the same direction, which is a stronger pattern than
Amendment 3's single gene at n=3 — and the two that clear are both mature
markers and therefore correlated, so the naive family-wise arithmetic does not
apply and no adjusted number is quoted here.

**And it is five of seven.** §3 fixes the remaining blocks (`222`, `242`) and
its order. Stopping here would be stopping at the first pass, which is what §3
exists to prevent.

---

## RESULT — four of seven blocks, 2026-09-07 (SUPERSEDED by the five-block table above)

`results/2026-09-07_a9b3d6d/`. Blocks `24H12439_A4` (231), `B-2080151-01-18`
(110), `B-1970164-01-04` (120), `B1891803_1-5` (210).

### Verdict: INDETERMINATE — the sign turns on below-floor blocks

| gene | per-block DiD | mean | 95% CI | |
|---|---|---|---|---|
| **GUCA2A** | −1.513, −0.838, −0.548, **−2.373** | −1.318 | **[−2.609, −0.027]** | excludes zero |
| GUCA2A, floored block dropped | −1.513, −0.838, −0.548 | −0.966 | **[−2.195, +0.262]** | includes zero |
| MS4A12 | −1.371, −0.759, −0.724, −2.687 | −1.385 | [−2.845, +0.075] | includes |
| CDX2 | −0.012, +0.215, +0.591, +0.404 | +0.198 | [−0.251, +0.647] | includes |
| EPCAM | −0.033, +0.381, −0.012, +0.052 | +0.097 | [−0.210, +0.404] | includes |

**Amendment 4 fired.** The four-block interval excludes zero; drop the one block
where the target sits below the negative-probe floor and it does not. That
block's DiD is a **bound, not a point** — its observed detection (0.0118) is
below what the noise process alone supplies (0.0167) — and a mean mixing bounds
with points estimates neither. **Neither interval is the answer.**

The rule was fixed and pushed in `db1664b` **before this aggregate was
computed**, precisely because the three-block result already showed which way
210 would push it.

### Amendment 3's caution was right, and one more block proved it

At three blocks **MS4A12 excluded zero** and GUCA2A did not. Amendment 3 said
that was a line in a table and not a finding: four informative genes, no
multiplicity control, `P(≥1 excludes | global null) = 0.19`, and a secondary
marker clearing while the primary does not is not the primary claim arriving by
another route.

**At four blocks MS4A12 includes zero** — [−2.845, +0.075]. The finding that was
not quoted has evaporated, which is what the caution was for.

### What holds across all four blocks

`GUCA2A` and `MS4A12` are negative in every block; `CDX2` and `EPCAM` are not
consistently signed. That is the two-block pattern surviving to four, and it is
**sign agreement, not an effect size** — §6's width at n=4 is 5.17× avenue A's.

### What this does not say

Everything in §8, unchanged, plus Amendment 2: the control is an epithelial
keratin, so **this statistic does not separate compositional from intrinsic**.
And §7's per-cell limit is untouched — no block licenses a "GUCA2A-low =
silenced" claim about any cell.

**Three blocks remain**: `221`, `222`, `242`. §3's order and its
one-at-a-time disk protocol are unchanged.
