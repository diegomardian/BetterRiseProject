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

## RESULT

*Not run. Only section `231` has been read, and that under the feasibility
document.*
