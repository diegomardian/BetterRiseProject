# Pre-registration — can the intrinsic reading be asked on WTx spatial data?

**Written:** 2026-09-07 · **Author:** W1 (Bode) · **Status:** proposed ·
**Depends on** [prereg_adenoma_decomposition.md](prereg_adenoma_decomposition.md)
and [prereg_becker_replication.md](prereg_becker_replication.md) (its RESULT is
why several sections below exist).

> **Written before a single byte is downloaded.** Nothing in the Crowell deposit
> has been read. §6's asymmetry, §5's gate and §3's width penalty are fixed
> here, and the section to read is named in §2 from the Zenodo file listing
> alone.

---

## 1 · What this is, and the three things it is not

`Crowell 2025` (bioRxiv `10.1101/2025.06.23.660674`) is **whole-transcriptome
CosMx SMI** — 19,867 genes, >99.5% of annotated protein-coding — on FFPE
sections carrying **reference mucosa, tubulovillous adenoma and carcinoma within
the same section, from the same patient.** Seven treatment-naïve patients, eight
colon sections plus one lymph node, ~3.4M retained cells. Public on Zenodo under
CC-BY-4.0.

**This document covers a feasibility read and nothing else**: does the panel
detect at usable rates inside each anatomical domain?

### Not a replication of avenue A

Avenue A is a per-patient Kitagawa decomposition at n=43. This is n=7 on a
different platform. §3 fixes what n=7 can carry, and it is not that.

### Not a rescue of the instrument question

`docs/HANDOFF.md` §6g closed the per-cell sensitivity question five ways. **In
situ imaging makes per-cell sensitivity worse than droplet, not better.** §7
states exactly what this deposit's controls do and do not license, because the
temptation to read them as the missing positive control is strong and wrong.

### Not a spatial niche analysis

Avenue B needs a per-cell "GUCA2A-low" call. That needs a per-cell
false-**negative** rate, which nothing here supplies (§7). B's objection is
untouched by this document.

## 2 · What is read, named before the download

**`232.h5ad`, 226.8 MB, md5 `f2271485bae235439267c1bd13f1a7c9`**, from Zenodo
`10.5281/zenodo.15574384`. The smallest of eight by a wide margin — the next is
`222.h5ad` at 240.3 MB and the largest, `120.h5ad`, is 2,284.5 MB. It fits in
the ~8.6 GB free on `/projectnb` without touching the ICBI atlas.

**Whether `232` carries TVA cells is NOT assumed.** The naming suggests slide 2,
region 3, and the H&E record suggests the TVA shapes sit there; neither is a
measurement. **The job reports the domain vocabulary it finds and maps
nothing**, exactly as `becker_feasibility --inspect` does, and for the same
reason: reading a label rather than the grouping once put Chen_2021's usable
pairs at zero when the true number was 44. If `232` carries no TVA, the named
fallback is **`231.h5ad`, 534.8 MB, md5 `fdd226ec7632a3d45a4f9dd46c24a981`** —
named here so the choice is not made after seeing a number.

## 3 · What n=7 costs, stated before any estimate

Invariant 5: the patient is the unit. **3.4M cells buy within-patient precision
and nothing across patients.**

| n | half-width `t(n-1)/√n` | × the n=43 width |
|---|---|---|
| 43 (avenue A) | 0.3078 | 1.00× |
| 15 (DIS/VAL halves) | 0.5538 | 1.80× |
| **7 (here)** | **0.9248** | **3.01×** |

**A contrast that cleared zero at n=43 must be about three times larger to clear
it here** — and 3.01× is well beyond the 1.80× that already made half the
DIS/VAL contrasts miss.

**And n=7 is one patient above the floor where the inverse-variance weight has
two moments** (`src/reference/meta_calibration.py`: `E[w]` finite for n≥4,
`Var[w]` for n≥6). The arithmetic says *barely*.

**The interval is Student-t. Not the percentile bootstrap.** At n=7 the
percentile bootstrap is **0.742× the correct width and excludes zero 12.0% of
the time under a true null** against a nominal 5% (`interval_calibration`,
`docs/HANDOFF.md` §3a). That is the largest miscalibration this project has
measured at any n it has used.

**So the primary object is the per-patient `m_T/m_N` ratio table** — scale-free,
comparable across genes, needing no decomposition and no denominator. Any
cross-patient contrast is pre-committed to be wide and is reported as such.

## 4 · The label is anatomical, and that is the point

Every maturity label this project has built has been made of transcripts, and
`prereg_becker_replication.md`'s RESULT is what that costs: the mature label
there was built from the mature-colonocyte program, GUCA2A belongs to that
program, and the one arm where GUCA2A outran its controls became unreadable
because selecting cells that kept the markers selects cells that kept the
target.

**Here the domain label is histology.** REF / TVA / CRC are pathologist-drawn
shapes over an H&E image. They are not a function of any transcript, so a
maturity or lesion assignment made from them **cannot be circular with respect
to GUCA2A by construction.** That is the rarest property in this project and it
is the main reason this deposit is worth the time.

Invariant 2 bars a target from its own label. Becker's RESULT proposes extending
it to genes **co-regulated with** the target — a `src/schema.py`-adjacent change
needing a PR and two approvals, **not taken here**. This document simply does
not rely on a transcript label at all, which satisfies both the current
invariant and the proposed one.

## 5 · The gate, fixed here

For each of the six panel genes, in each domain:

1. **detection** — fraction of cells with ≥1 count
2. **the negative-probe floor** — fraction of cells with ≥1 count on any of the
   50 negative probes, and the false-code rate over the 2,745 false codes
3. **separation** — the gene's detection against that floor, on the detection
   scale (`cloglog(p) = log(mu)`, `docs/HANDOFF.md` §6d), because a ratio of
   probabilities is not comparable across genes with different baselines and
   this repository has now reintroduced that error twice
4. **depth** — median counts and median genes per cell, per domain

**A gene is usable in a domain when its detection separates from that domain's
negative-probe floor by `MIN_LOG_SEPARATION` on the detection scale.**

**The decision that matters is GUCA2A and MS4A12 in the TVA.** GUCA2A failing
there ends the spatial reading of avenue A's estimand, exactly as it ended B1.

**Depth is reported per domain because it is the confound that already fooled
this project once.** Becker's mature label carried 2.04× the arm's median UMIs
and lifted every gene; the pass was library size. If the TVA and REF domains
differ in depth, every cross-domain comparison inherits it, and the per-domain
depth row is what lets a reader see that before reading anything else.

## 6 · The asymmetry, pre-committed: a positive is safe, a null is not

**Every REF region in this deposit is adjacent normal from a patient who has a
carcinoma. There is no healthy-donor arm.**

`prereg_becker_replication.md`'s RESULT is why this is written down before the
numbers. In Becker, FAP "unaffected" mucosa was GUCA2A-depleted **2.38× beyond
its own controls** relative to true healthy-donor colon — and that was only
visible *because* Becker kept a healthy-donor arm as a fifth stage. Crowell has
no such arm, so the same depletion could be present here and be invisible.

| outcome | reading |
|---|---|
| GUCA2A separates in REF and falls in TVA | **Safe.** A field-affected reference understates the fall, so the true effect is at least this large. |
| GUCA2A does not fall from REF to TVA | **Uninterpretable, and must not be reported as a negative.** A reference already depleted cannot show a further fall. |
| GUCA2A fails to separate from the floor in any domain | A statement about **in-situ sensitivity**, not about biology — the same distinction B1's §3 drew, and it held there. |
| The controls (ACTB, KRT8) fail to separate | The read is refused outright: the instrument is not measuring this tissue. |

**This asymmetry is fixed here so that a null cannot later be quoted as
evidence of absence.**

## 7 · What the negative probes are, and what they are not

The 50 negative probes and 2,745 false codes bound the **per-cell
false-POSITIVE floor**: how often the instrument reports a transcript that is
not there. **That is a specificity control, and it is the best instrument
control any dataset in this project has carried.**

**It is not a sensitivity control and it is not §6g's missing positive
control.** A negative probe cannot tell you what a real GUCA2A transcript you
failed to capture would have looked like. The per-cell false-**negative** rate
remains unmeasured here, as it was on every other substrate.

**So this deposit does not license any per-cell "GUCA2A-low = silenced" claim**,
which is avenue B's whole premise. Pre-committed, because the specificity
control is genuinely good and that is exactly what makes it temptingly quotable
as the thing it is not.

## 8 · What this read cannot decide

- **Per-cell sensitivity** (§7), and therefore B's niche claim.
- **Whether avenue A replicates.** n=7, different platform, §3.
- **Whether the REF domain is field-affected.** No healthy arm exists here; §6
  handles it by asymmetry rather than by measurement.
- **Anything about the other seven sections.** One section is read. A second is
  a separate decision with its own cost.

## 9 · Standing

**A feasibility read, pre-registered, on public CC-BY data, deciding one thing:
whether the six panel genes are measurable inside pathologist-drawn domains.**
It cannot make any claim about silencing. It can tell us whether the question is
askable here at all, which is what B1 spent a week failing to establish on a
different deposit.

---

## Amendment 1 — 2026-09-07, after the inspection and before the verdict is re-read

Five corrections. Three are things the deposit turned out not to be; two are
errors in this document. **One is a deviation that already happened and is
recorded as one rather than absorbed.**

### 1 · The floor is the MEAN PER-PROBE rate, not the union over 50

§5 item 2 specified *"fraction of cells with ≥1 count on any of the 50 negative
probes."* That is a union over 50 features and a gene's detection is one
feature; it overstates the floor by roughly the probe count and would fail every
gene including the controls. **The gate uses the mean per-probe detection
rate**, which is the comparable quantity. The union rate is retained and
reported as a **contamination QC indicator** under its own name.

The implementation has done this since `07d5300` and §5 was never amended to
match, so for one commit-range the pre-registration and its own code specified
different statistics. **The code was right and the document was stale**; both
are now the per-probe rate.

### 2 · The control probes are in `obs`, not in `var`

Zero features matched the negative-probe naming, correctly — this deposit
summarised the controls per cell (`nCount_negprobes`, `nFeature_negprobes`,
`nCount_falsecode`) before writing the object. `require_controls` refused the
run rather than proceeding on a floor of zero, which is the behaviour §5
depends on. The floor is now built from those obs columns, with the probe count
supplied as 50 from the paper and recorded in the sidecar.

### 3 · ACTB is absent, so the control branch rests on ONE gene

The whole-transcriptome space carries 18,936 features and ACTB is not one of
them — as it was not on the CosMx 1K or 6K panels. **§6's fourth branch was
written as "the controls (ACTB, KRT8) fail to separate" and now concerns KRT8
alone.** EPCAM is `epithelial` in `GENE_ROLES`, not `control`, and must not be
counted as a second control to make the branch look wider than it is.

**This materially weakens every sensitivity statement this document can make.**
One control gene establishes that global capture happened; it cannot establish
anything gene-specific, and §7 already said so.

### 4 · DEVIATION: the fallback fired, for the opposite reason than written

§2 said: *"If `232` carries no TVA, the named fallback is `231`."* **`232`
carries ONLY TVA** — `tissue` is `TVA` for all 130,814 cells. The file actually
read was `231`, chosen after the inspection because `232` lacked a **reference**
domain, which is not the condition §2 states.

**The file was pre-named and the reason was not.** That is a post-inspection
deviation, and it is recorded here as one. It is a mild deviation — `231` was
named in this document before any download and no third file was considered —
but a fallback whose trigger is rewritten after looking is exactly the mechanism
this project distrusts, and calling it "the prereg named 231" without the rest
would be the overstatement.

### 5 · "Within the same section" describes the tissue, not the files

§1 quotes the paper correctly: REF, TVA and CRC occur within one physical
section. **The deposit's objects are split by field-of-view range**, so one
`.h5ad` need not carry all three. `metadata.txt` (Zenodo `10.5281/zenodo.15550908`)
shows `231` and `232` as FOV ranges 1–146 and 147–217 of the same block
`24H12439_A4`, same run and replicate. `231` carries all three domains; `232`
carries one.

### 6 · Below-quantification, NOT censored

A target under `MIN_LOG_SEPARATION` has not been censored by the assay. It fell
below **an analyst-chosen usability bar** set in this document. The correct
description is **near-floor / below the limit of quantification**, and nothing
here constructs a background-adjusted bound that would license the stronger
word. Recorded because "censored" was used in review and is a technical claim
this design does not support.

### What Amendment 1 does NOT change

The gate, the bar, the primary object, §6's asymmetry and §7's
specificity-only limit are untouched. **No threshold moved and no branch was
added to make an outcome pass.**

---

## Amendment 2 — the verdict is decomposed, because §6 has no branch for what happened

§6's branch table anticipates *"GUCA2A separates in REF and falls in TVA"* and
*"fails to separate from the floor in ANY domain."* The observed outcome is
**neither**: GUCA2A separates in REF and sits **below the bar** in TVA. A single
label cannot carry that without discarding one half of it, and the first
implementation discarded the half that mattered — it returned `NOT ASKABLE IN
THE ADENOMA`, which is true of the per-cell reading and silent about a
directional signal §6 had pre-specified.

**So the verdict is emitted as fields, not as one word:**

| field | meaning |
|---|---|
| `per_cell_feasibility` | did the target clear the usability bar in the adenoma domain |
| `prespecified_directional_read` | §6's REF→adenoma direction, at section level |
| `global_sensitivity_control` | whether the control-role gene's separation is worse in the adenoma domain |
| `control_referenced_pattern` | **exploratory**, post-hoc, and labelled in its own column |
| `patient_n` | the biological unit, which here is **1** |

**`verdict()` does not return `ASKABLE IN THE ADENOMA` on this outcome.** The
target did not become independently measurable per cell and no relabelling makes
it so.

### The control-referenced comparison is an addition, not a repair

Referencing each gene's REF→domain change against the control's is **not
specified in §5** and is not a missing line of pre-registered code. It is a new
exploratory analysis on a new estimand — Becker's `enrichment_audit` compares
transcript-defined subsets *within* an arm, and this compares *across
histological domains*, which is a different question. It is emitted with
`exploratory = True` and it carries **no interval**: a min/max band over one or
two control genes has no patient-level uncertainty in it, and the biological
unit here is one patient.

### CDX2 is used continuously, never as a label

CDX2 missed the bar in the TVA by **0.0031** on a bar of 1.0986. A binary flag
at the fourth decimal is not a measurement. CDX2's separation is reported as a
number and enters the exploratory pattern continuously; the `usable` flag is not
read for it.

---

## RESULT

**Inspection and first gate run 2026-09-07; gate result pending a required QC
rerun.** Sections 232 and the pre-registered fallback 231 were downloaded and
their checksums matched. Section 231 carries the explicit histopathological
domains `231_REF`, `231_TVA`, and `231_CRC` in `typ`.

The first gate run is **invalidated, not interpreted**: it used all 298,151
rows, including 19,460 cells with `fil=False`. The deposit's own Zenodo record
defines `fil` as whether a cell passed quality control, so the analysis
population is the 278,691 cells with `fil=True`. The runner now requires this
field, refuses non-logical or missing values, filters before computing either
detection or the negative-probe floor, and records the before/after counts in
the sidecar. This correction was made before reading a filtered result.

The first QC-corrected rerun retained 278,691 cells and confirmed that all
19,460 QC failures lacked a `typ` assignment: every annotated REF/TVA/CRC row
was bit-for-bit unchanged. It also exposed a second output defect before the
table was committed: 20,393 QC-passing cells with missing `typ` were
stringified into a literal `"nan"` domain and scored. Missing histopathology is
not a fourth domain. The runner now excludes those cells from the per-domain
table and records their count in the sidecar. This does not change any named
domain or the gate verdict, but the table is rerun once more so the published
artifact contains only real histopathological domains.

Neither superseded table nor its verdict should be quoted. **A third
supersession followed**: Amendment 2 replaced the single-word verdict, because
`NOT ASKABLE IN THE ADENOMA` is true of the per-cell reading and silent about a
direction §6 had pre-specified. Rerun with the reference domain and the
biological n supplied, so both halves are emitted:

```bash
python -m src.reference.jobs.crowell_feasibility \
  --object "$BRP_DATA_DIR/raw/crowell/231.h5ad" \
  --domain-column typ \
  --adenoma-domain 231_TVA \
  --reference-domain 231_REF \
  --n-negative-probes 50 \
  --patient-n 1
```

### What the superseded numbers say, replayed through the corrected verdict

Not a result. **Three tables exist for section `231` and two are superseded**: `results/2026-09-07_d3d1369/` (pre-QC, pre-`nan`-exclusion, old verdict) and `results/2026-09-07_194b669/` (post-QC, old verdict). **`results/2026-09-07_6e7e93e/` is the one to quote.** Its named-domain numbers are bit-identical to both, which confirms that neither the QC filter nor the `nan` exclusion moved a real domain — only the verdict changed. Recorded so the rerun is not the first
time anyone sees the shape of it.

| field | value |
|---|---|
| `per_cell_feasibility` | **failed** — GUCA2A +0.700 in `231_TVA` against a bar of +1.099 |
| `prespecified_directional_read` | **fall_observed** — +1.571 in `231_REF` → +0.700 in `231_TVA` |
| `global_sensitivity_control` | **not_worse** — over **one** control-role gene, ACTB being absent |
| `patient_n` | **1** |

**Below the limit of quantification, not censored.** The per-cell spatial
reading is not licensed and no relabelling makes GUCA2A independently
measurable there.

**The exploratory control-referenced pattern, `exploratory=True`, no interval.**
REF→TVA change in separation: KRT8 **+0.642** (the entire control band, one
gene), CDX2 +0.629, EPCAM +0.609, MS4A12 **−0.729**, GUCA2A **−0.871**. The two
targets move opposite in sign to the control while CDX2 sits beside it — the
same two-block shape avenue A found on Chen, on a different platform. **At n=1,
post-hoc, on a band of one gene, that is a hypothesis with a number attached and
nothing more.** It may not be quoted as replication or as silencing.
