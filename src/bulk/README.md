# W3 — Bulk & clinical

**Owner:** can be a strong analyst without scRNA-seq background · **Env:** `env/w3_bulk.yml` → `conda activate brp-w3`
**Branch prefix:** `w3/…` · **Blocked by:** nothing — fully independent for five weeks

Laptop-fine. The single-cell arm needs 32–64 GB and a GPU; TCGA-COAD/READ does
not.

## What you deliver

| Wk | Task | Done when |
|----|------|-----------|
| 1 | GDC ingest, COAD + READ, STAR counts; normalise (**keep both** TPM and log-CPM) | Matrix on the shared gene index |
| 2 | **Premise check: distribution of GUCA2A and CDX2 in TCGA-COAD** | Histogram + bimodality test. Report to the whole team in week 2. |
| 2–3 | Tumour purity: ESTIMATE; pull precomputed ABSOLUTE calls where available | Purity per sample, with method noted |
| 3 | Batch and technical structure: plate, TSS, sequencing batch | Documented; confounding with stage and MMR tested |
| 3–4 | Clinical table from TCGA-CDR — DSS and PFI primary, OS secondary | Curated, censoring rules explicit |
| 4–5 | Covariate set pre-specified and **locked**: stage, age, sex, MMR/MSI, purity, site | Written down before any survival model runs |
| 5 | Baseline survival models on clinical covariates alone | Sanity check — if stage is not prognostic, something is wrong upstream |

## The week-2 premise check can redirect the project

"Assuming bulk GUCA2A is negligible" is an untested premise (Executive Brief
error #6). You test it. If loss is **continuous rather than bimodal**, the
two-type classification dissolves into a regression and the project changes
shape. That is a week-2 finding, not a week-12 one — bring it to the standing
meeting the moment you have the histogram.

## Survival endpoints

**DSS and PFI from TCGA-CDR (Liu et al. 2018).** OS is secondary and is reported
as such — COAD OS is heavily contaminated by non-cancer death (CLAUDE.md
invariant 9). Do not let OS become the headline because it is the column that
was easiest to pull.

## The shared gene index

Your bulk matrix and W1's S matrices sit on the same index
(`config/gene_index/`), which is what makes integration a join rather than a
negotiation. If a gene symbol mapping decision comes up (deprecated symbols,
multi-mapping, version drift), it affects W1 too — raise it, do not decide it
alone.

## Gotchas

- COAD and READ have **different treatment patterns**. Decide now whether to
  pool or stratify, and write the reason into the PR.
- MSI status in TCGA is incompletely annotated. Check coverage **before**
  committing to it as the pre-registered subgroup variable — if it is sparse,
  the subgroup variable changes.
- Normal-adjacent samples in TCGA-COAD are few and **not matched to all
  tumours**. Do not assume pairing.
- Purity is a pre-specified covariate in the survival models, not a robustness
  check bolted on later.

## Stage 4, later

Deconvolve for **mature-colonocyte fraction only**. No cell-type-specific
expression imputation from bulk, anywhere, ever (CLAUDE.md invariant 6) — bulk
gives fractions at r ≈ 0.92, and intrinsic estimates come back attenuated
×0.6–0.8 in the direction of our prior hypothesis. The Stage 4 question is how
much of the variance in bulk GUCA2A and CDX2 the fraction alone explains.
