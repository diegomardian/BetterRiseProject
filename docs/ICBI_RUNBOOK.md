# Path B — the coexpression premise at ICBI scale

> ## THIS HAS ALREADY BEEN RUN, AND IT TERMINATED.
>
> **2026-09-05: the premise is UNRESOLVED at 13 studies and 122 patients**, and
> not for lack of power — the studies disagree with each other. Tables:
> `results/2026-09-05_3380d15/` (per-study), `results/2026-09-05_61ba221/`
> (meta). Summary in `docs/HANDOFF.md` §2 and §6b.
>
> Everything below is retained for re-derivation and for anyone extending the
> reading to another atlas.

---

## What the question was

The decomposition cannot separate compositional from intrinsic loss on this
panel (the algebra collapses — `docs/HANDOFF.md` §2). The coexpression reading
sidesteps that by measuring per-cell detection *inside* a fixed label. But it
only means anything if the two arms' surviving mature cells are the same kind of
cell, and that premise came back **UNRESOLVED on all three original cohorts**.

That was the one blocker that looked like a *power* problem — the intervals
straddled the tolerance because there were 7, 5 and 30 patients. So it was worth
taking to a bigger sample. The ICBI CRC atlas carries 14 candidate studies at
≥100 epithelial cells per arm, against the 3 we had.

## The answer

| control | pooled | I² | verdict |
|---|---|---|---|
| ACTB | +0.152 [−0.013, +0.317] | 62.8% | HOLDS |
| KRT8 | −0.453 | **87.6%** | UNRESOLVED |

Every control must hold, so the premise does not. **KRT8's per-study estimates
run −1.177 to +0.088.** Per-study verdicts: 3 HOLDS, 1 REFUSED, 7 UNRESOLVED,
2 UNDEFINED.

Because the premise does not hold, **the per-gene detection deltas are not read
at the meta level.** That gate is the point: a marker falling inside a
population that has itself changed is not silencing, and that is as true of 13
studies as of one.

**This settles "can more data make it identifiable". It cannot.** Four times the
studies and four times the patients converted a precision problem into a
disagreement problem.

### One inference to avoid, because it was made here

On the three original cohorts ACTB pooled to +0.487 with **I² = 0.0%**, which
was read as "the cohorts agree, so the failure is precision rather than
disagreement", with a prediction that ACTB could only resolve "by a hair". At
k=11 ACTB pools to **+0.152** and holds comfortably — the point estimate moved
0.335. Three similar cohorts (range +0.431 to +0.586) were not a sample of the
fourteen (range −0.589 to +0.532).

**I² near zero on a small, similar set is not evidence of homogeneity in the
population.** It can just mean the set is small and similar.

---

## Path C: the adenoma reading on the same atlas

`Chen_2021_Cell` IS the Vanderbilt/HTAN polyp atlas. 44 patients with a matched
polyp and their own normal; the reference arm is labelled `healthy normal` and
is patient-matched, established by `patient_id` rather than by the label.

    qsub -v BRP_ICBI_STUDY=Chen_2021_Cell,BRP_ICBI_ARMS=adenoma \
         src/reference/jobs/icbi_coexpression.sh

Scores `lineage` and `best4`. The read bar (`ADENOMA_READ_BAR`) is printed at
the top of the log and stored in the sidecar, composed before the run.

**Result** (`results/2026-09-05_3a1af9f/`): the premise HOLDS at both rungs —
the first time anywhere in this project — and `best4` is estimable, 32 of 44
patients, median ~85 mature cells against carcinoma's median of 3. GUCA2A falls
−0.174 [−0.243, −0.109] inside the mature population.

**But it is not gene-specific.** GUCA2A − ACTB is −0.139 [−0.197, −0.081] and
excludes zero; GUCA2A − MS4A12 is −0.009 [−0.064, +0.045] and contains zero. The
mature programme is down, not GUCA2A. See `docs/HANDOFF.md` §6d.

**Two things the adenoma run also establishes**, both in
`docs/NEXT_AVENUES.md`: the decomposition's algebraic collapse does NOT fire
here (m_T/m_N runs 0.37–0.95), so the original method may be identifiable; and
`summarise` must be called PER RUNG — calling it once across both produced rows
reading `n_patients = 64`, which is 44 plus 20 and neither the same patients nor
the same cells.

---

## The three pieces

| piece | module | notes |
|---|---|---|
| extraction | `src/reference/icbi_slice.py` | CSR row reads, vocabulary maps |
| adaptation | `src/reference/jobs/icbi_coexpression.py` + `.sh` | same statistic, new loader |
| meta-analysis | `src/harness/meta.py` + `src/reference/jobs/coexpression_meta.py` | DerSimonian–Laird |

`src/harness/meta.py` is **general and reusable**. Invariant 4 has required
random-effects meta-analysis in three documents since week 0 and nothing
implemented it until now.

---

## Re-running it

```bash
export BRP_DATA_DIR=/projectnb/rise-batteries/bode/guanylin/data
export BRP_ICBI_DIR=/project/rise-batteries/bode/icbi

# 1. the atlas, if absent. 30.44 GiB, resumable, ~25 min.
qsub src/reference/jobs/fetch_icbi_atlas.sh

# 2. ONE study first. Pelka_2021_Cell IS GSE178341, the only one of the
#    fourteen with a committed result to check against.
qsub -v BRP_ICBI_STUDY=Pelka_2021_Cell src/reference/jobs/icbi_coexpression.sh

# 3. only if that PASSES — the job exits 4 if it does not.
qsub -v BRP_ICBI_STUDY=all src/reference/jobs/icbi_coexpression.sh

# 4. local, over committed tables. No cluster.
python -m src.reference.jobs.coexpression_meta
```

The validation on Pelka is mechanical, not a matter of judgement: it reads
`results/2026-09-04_975cf5c/`, requires the same premise verdict word, an
overlapping ACTB **log2** interval, and a GUCA2A detection delta within 0.15.
Observed drift was **0.0096**.

---

## Things that will bite you

**`/X` is log1p-normalised. Raw counts are in `layers["counts"]`.** `adata.X` is
what any obvious code reaches for, and detection at ≥1 UMI against log values is
wrong while nothing raises. Measured: `/X` 0.2795–5.404 over 3525 distinct
values; `layers/counts` 1–289 over 155, starting `[1,1,1,1,2,1,2,5]`.
`assert_raw_counts` checks values, not the layer name.

**`/var/_index` is Ensembl; the symbols are in `GeneSymbol`.** Reading only the
index reports every gene absent from an atlas containing them all. That
identifier-space error has now happened four times in this repo, once inside a
probe written to prevent it.

**The atlas lives on `/project`, not under `BRP_DATA_DIR`.** `/projectnb` had
15 GB free of a 50 GB quota. `BRP_ICBI_DIR` moves only this download. **Never
`/scratch`** — it is `/dev/sda8`, node-local and purged, so a file one job
writes is invisible to the next.

**Slicing is positional.** `load_obs` compares the cached obs length against the
atlas's own row count, because a stale cache would hand every cell another
cell's metadata, silently.

**`read_cells` sorts rows then permutes back.** HDF5 fancy-indexing needs
increasing indices; forgetting the second half gives every cell someone else's
expression.

**Sorted fractions must be excluded.** GSE178341 filters
`PROCESSING_TYPE == UNSORTED`; the atlas's analogue is
`enrichment_cell_types == "naive"`. Pelka is 210,667 naive against 130,019
CD45+, and a CD45-sorted fraction is immune-enriched by construction. Note
`Liu_2024_Cancer_Res` contributed **nothing** — 0 of its 176,531 cells are
naive.

**QC runs on the FULL patient block, all compartments.** The epithelial subset
is only the scoring population. Extracting epithelium alone changes the QC basis
and would have made the Pelka comparison diverge for a reason unrelated to
whether the adaptation was right.

**Per-study, one at a time, slices discarded.** The 136 qualifying patients hold
1,398,508 cells across all compartments — ~9.85 GB gzipped against ~10 GB free.
Nothing intermediate is kept.

---

## If you extend this

The premise check switches statistic per study: **detection** (tolerance 0.10)
where the control has headroom, **log2 expression** (tolerance 0.5) where it is
saturated. The meta layer pools log2 uniformly across all studies. That is a
deliberate choice — one statistic, comparable across studies — but it means
**the per-study verdict tally does not predict the meta verdict**, and the two
should not be read as answering the same question.

Two studies (`Ji_2024`, `Zheng_2022`) fall below `MIN_PREMISE_PATIENTS = 3` and
are excluded from the pool rather than contributing a two-patient estimate.
`Liu_2024` contributes nothing at all. So k=11, not 13.

**The caveat scale cannot touch.** Survivorship — GUCA2A-high cells having been
preferentially destroyed — is not transcript-detectable. A resolved premise
would make the detection reading interpretable; it would never rule that out.
