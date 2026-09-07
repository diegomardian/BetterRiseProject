# Pre-registration — replicating the adenoma decomposition on Becker FAP

**Written:** 2026-09-06 · **Author:** W1 (Bode) · **Status:** proposed, awaiting
team ratification · **Depends on**
[prereg_decomposition_statistic.md](prereg_decomposition_statistic.md) (`ac7eca1`)
and [prereg_adenoma_decomposition.md](prereg_adenoma_decomposition.md).

> **This is the one open item on avenue A that has a closing path.** The others —
> `best4`'s intrinsic arm, the three-point curve, survivorship — are properties
> of the data, not of the code, and no amount of re-analysis touches them. The
> statistic being post-hoc *is* closable, and it closes by running the same
> design on a second substrate with the statistic already fixed. It was fixed in
> `ac7eca1`, before this document and before the substrate has been fetched.

> ### Amendment 1 — 2026-09-06, from the GEO metadata, before any download completes
>
> `GSE201348`'s series matrix resolves §6's unverified items and raises two
> design questions the original text did not cover. Both are decided here,
> **before the loader exists**, because both change the estimand.
>
> **The format risk is dead.** 72 samples, standard 10x triplets
> (`barcodes.tsv.gz` / `features.tsv.gz` / `matrix.mtx.gz`), 1.2 GB in one
> `_RAW.tar`. Not Seurat objects. The loader is small.
>
> **The arm vocabulary, read rather than assumed** — from
> `Sample_characteristics_ch1`, `disease stage`:
>
> | value | maps to |
> |---|---|
> | `Polyp` | **tumour** (diseased arm) |
> | `Unaffected` | **normal** (reference arm) |
> | `CRC` | **excluded** — see below |
>
> Donor is the `Sample_title` prefix (`A001-C-007` → `A001`). FAP status is a
> per-donor characteristic. **No metadata lives in the tar**, so the series
> matrix is a required second download and the mapping above is the whole of
> what §3's "decision for a human who has seen the list" amounts to.
>
> **1. Multiple polyps per donor — POOL PER DONOR, and the reason is that this
> is a replication.** Donor A001 alone carries 6 polyps, 3 unaffected and 1 CRC.
> Chen_2021 gave one polyp and one normal per patient, so invariant 5's "the
> patient is the unit" was unambiguous there; here it is not.
>
> Two estimands are available and **they are not interchangeable**: pooling all
> of a donor's polyps into one arm reproduces Chen_2021's shape exactly, while
> computing per lesion and aggregating within donor preserves between-lesion
> variation. The second is arguably the better *design* — between-lesion
> variation is exactly where a field-effect or Wnt-tone mechanism would live.
> **It is the wrong choice here anyway**, because a replication that changes the
> estimand is not a replication, and B1's entire purpose is to test whether
> Chen_2021's result holds with the statistic pre-fixed.
>
> So: **pooled-per-donor is primary and confirmatory. Per-lesion is secondary
> and exploratory**, reported beside it and never substituted for it. If the two
> disagree, that disagreement is a finding about lesion heterogeneity and not a
> reason to prefer whichever replicates.
>
> **2. Technical replicates — POOLED.** `A002-C-010` appears as `Replicate1` and
> `Replicate2`: one physical sample sequenced twice. Their cells are
> concatenated, because they are the same biological unit and the unit of
> inference is the donor. Dropping one discards data; averaging is wrong at the
> cell level. **Which samples had replicates is recorded in the report**, since
> a pooled replicate has more cells than its neighbours and that shows up in
> every per-sample count.
>
> **3. `CRC` samples are excluded, and what that forgoes is worth naming.** This
> is the adenoma reading; carcinoma is where four routes already terminated.
> But Becker carries **CRC, polyp and unaffected tissue from the same donors**,
> which is a normal→polyp→carcinoma gradient at better than the n=3 that
> `docs/NEXT_AVENUES.md` §1c is limited to. That is a real bonus avenue and it
> is **explicitly out of scope here** — it needs its own design, and folding it
> in would be the estimand drift this amendment exists to prevent.

---

## 1 · What is being replicated, stated exactly

Not "the adenoma result". **The six cross-block contrasts that survived
`prereg_decomposition_statistic.md`'s own agreement rule at `lineage`**, and
within those, the four that matter most:

| claim | status on Chen_2021 (n=43) | replication target |
|---|---|---|
| **GUCA2A separates from ACTB, CDX2, EPCAM, KRT8** | **4 of 4, unanimous on all four statistics** | **primary** |
| MS4A12 separates from ACTB, KRT8 | 2 of 4 (fails vs CDX2, EPCAM) | secondary |
| `GUCA2A − MS4A12` contains zero | 3 of 4 statistics | primary, and it is a *null* |
| `i/c` does not collapse onto a constant | bracket −0.05 to −0.63, nothing near −1 | primary |
| Two blocks at `best4` | **1 of 8 survives — retracted** | **not replicated; no claim to test** |

**The primary claim is about GUCA2A specifically, not about a tier**, because
the agreement rule left GUCA2A's four contrasts unanimous and MS4A12's split.
That is narrower than the RESULT's first headline and it is the version that
survived its own scrutiny.

## 2 · Why Becker, and one thing about it that is not a coincidence

`GSE201348` (snRNA-seq) and `GSE201349` (scATAC), Becker et al. 2022 — FAP and
sporadic polyps with matched normal, many polyps per patient.

**It is the substrate this project froze its third labelling axis against, in
week 0.** `config/labeling_axes.yaml` — frozen, two approvals to change — names
axis 3 as *"A different measurement — chromatin accessibility"*, source
**"Becker/Chang multiome"**, with the caveat *"Not transcript-based, and
therefore the strongest defence against label leakage."*

So this is not a new avenue chosen because the old one worked. It is the
substrate the design named before any result existed, and it has never been
fetched. The scATAC arm additionally answers the circularity objection every
reading in this project carries — that cell identity is conditioned on
transcription — which no amount of scRNA replication can.

**This pre-registration covers the snRNA arm only.** The chromatin axis needs
its own design and is not smuggled in here.

## 3 · The feasibility gate, which comes first and can stop everything

**Nuclear RNA is the risk, and it is quantified rather than feared.** GUCA2A and
MS4A12 are cytoplasmic transcripts; snRNA-seq samples nuclei. The panel's
baseline detection in Chen_2021's normal-arm mature cells:

| ACTB | KRT8 | EPCAM | CDX2 | GUCA2A | **MS4A12** |
|---|---|---|---|---|---|
| 0.984 | 0.958 | 0.900 | 0.822 | 0.437 | **0.363** |

MS4A12 is the floor. A protocol that halves cytoplasmic detection puts it near
0.18; one that quarters it puts it near 0.09.

### Pre-committed gate, per gene and not globally

A gene enters the replication only if, in the **normal** arm's mature cells:

1. detection ≥ **0.10**, and
2. a **non-zero** per-cell mean in ≥ **75%** of patients.

**Genes failing the gate are dropped and named. Contrasts involving them are not
attempted** — not attempted and reported as absent, rather than attempted and
reported as null, because a gene that is not measured has not been found equal
to anything.

**Per gene, not globally, because the claims are not equally exposed.** MS4A12
failing costs the secondary claim and leaves the primary intact. GUCA2A failing
ends the replication:

| gate outcome | consequence |
|---|---|
| GUCA2A fails | **The replication cannot be run.** Report the detection table as a measurement about snRNA-seq and stop. Not a negative result about the biology. |
| MS4A12 fails, GUCA2A passes | Primary claim testable, secondary not. `GUCA2A − MS4A12` untestable, so the "not gene-specific" null is **not** replicated and stays single-cohort. |
| ≥1 control fails | The comparator set is smaller. State which, and that cross-block counts are out of fewer than 8. |
| all pass | Full design. |

**This gate is read before any decomposition is computed**, and its table is
committed before the reading runs.

## 4 · The analysis, entirely inherited

Nothing new is specified here. That is the point.

- **Scoring:** `icbi_coexpression`-shaped — same QC, same labeller, same depth
  matching, same `--collect-fractions`. A different loader.
- **Decomposition:** `decompose_cohort`, all three weightings, two denominators
  (`prereg_adenoma_decomposition.md` §3.1 as amended).
- **Statistic:** `log_ratio`, **already fixed** in `ac7eca1`, with the agreement
  rule against `share_abs`, `share_signed` and `ratio`.
- **Interval:** Student-t over patients. Not the percentile bootstrap
  (`docs/HANDOFF.md` §3a).
- **Estimability:** both rules, reported separately, never folded — and under
  **all three** candidate cutpoint sets, because §3a-bis showed there is no
  single calibrated value.
- **Rungs:** whatever the cohort supports. `epithelial` included for its
  degeneracy check; `crypt_position` reported only if it forms three bins, since
  on Chen_2021 it collapsed onto `lineage` for 41 of 44 patients.

## 5 · What would falsify the adenoma reading

| branch | consequence |
|---|---|
| GUCA2A's four contrasts replicate, unanimous | **The primary claim holds on two substrates with the statistic pre-fixed on the second.** This is what the avenue is for. |
| GUCA2A's contrasts fail to replicate | **The Chen_2021 result is single-cohort and does not generalise.** It is not thereby wrong, but it stops being quotable without the qualifier, and the honest fallback is the `m_T/m_N` ratio table, which needs no decomposition. |
| `GUCA2A − MS4A12` **excludes** zero here | The "not gene-specific" null does not replicate. Report both; do not average them. |
| `i/c` collapses toward −1 | The estimand is not identifiable on this substrate, and identifiability is a property of the cohort rather than of adenoma-versus-carcinoma. That would materially weaken §1 of the adenoma prereg. |
| The premise does not hold | Not a negative result. The decomposition does not gate on the premise, but a substrate where the two arms are not comparable limits what the split means, and it is reported. |

## 6 · What cannot be checked from here, and must be before anyone runs this

**Stated as unverified rather than asserted.** No accession, file layout or size
below has been confirmed against GEO from this session:

1. ~~**The accessions and their contents.**~~ **`GSE201348` VERIFIED 2026-09-06**
   against the GEO FTP listing: 72 samples, 1.2 GB. The cell counts in the
   proposal (201,884) are still unchecked and will be known when the tar is
   read. `GSE201349` (scATAC) has NOT been verified and is out of scope.
2. ~~**The file format.**~~ **RESOLVED 2026-09-06 — Amendment 1.** Standard 10x
   triplets, 72 samples, 1.2 GB. Not Seurat. The loader is small.
3. **Disk — MEASURED 2026-09-06, and it scopes this design.** `pquota`:

       /project/rise-batteries      45.01 / 50 GB   ->  4.99 GB free
       /projectnb/rise-batteries    40.17 / 50 GB   ->  9.83 GB free

   They are separate filesystems, so the usable figure is **9.83 GB**, not the
   sum. That is the constraint this document already lives within by covering
   **the snRNA arm only**: a 200k-nucleus matrix plausibly fits, and the scATAC
   arm (447k cells, fragment files) does not — not by a little.

   **So axis 3 is not merely unscheduled, it is unaffordable at current
   occupancy.** The chromatin arm needs space freed or another filesystem before
   it can be designed, and that is a storage decision rather than a science one.
   The 30 GB ICBI atlas on `/project` is the obvious candidate to move or drop,
   and it should not be dropped while §6g and avenue A still read from it.

**No download job is written in this commit**, because writing one against an
unverified layout is how a cluster job dies at 3am.

## 7 · Standing

**Confirmatory for the statistic, exploratory for the substrate.** The statistic
was fixed before the data was fetched and this document was written after that
commit; the substrate is new and its feasibility is unknown.

It needs the team: it is a new cohort, a new download against tight quotas, and
it is the test that decides whether avenue A's result is a finding or a
single-cohort observation.

---

## RESULT

*Not run. Blocked on §6 — accession verification, file format, and disk.*
