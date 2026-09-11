# Chen lesion-subtype labels — terms and provenance gate

**Status:** gate conditions 1–6 of `docs/prereg_chen_lesion_subtype.md` §9 are
discharged. Terms reviewed 2026-09-10; label provenance confirmed at source;
the snapshot is pinned and checksummed; the identifier join reproduces every
frozen arm; both label families are declared and pass; the attrition table is
produced and inspected. Condition 7 is waived and recorded.

This gate licenses **joining a pre-registered label to a decomposition value**.
It does not license reading an outcome, and it does not turn a morphology label
into a biological endpoint.

## 1 · Data-use terms — condition 1

A public API is not a licence, so the terms were read rather than assumed.

**cBioPortal.** The portal's user documentation states: *"Unless otherwise
noted, data in cBioPortal are available under the ODC Open Database License
with no restrictions on the use of the data, as long as you properly give
attribution to the original studies."* It further states that *"There are some
studies that restrict the commercial use of the data, but that will be
explicitly mentioned in the study information."*

**This study carries no such note.** The pinned `study.json` records
`publicStudy: true`, `readPermission: true`, an empty `groups` field, and a
description containing no licence, consent, embargo, commercial or dbGaP
language. The default therefore applies.

**What that obliges.** Attribution to the original study — Chen et al., *Cell*
2021, PMID 34910928 — in any output, plus the cBioPortal platform citation
(Cerami et al.). ODbL governs **redistribution**; this project redistributes
nothing, and the snapshot under `data/raw/` is gitignored and checksummed in
`data/manifest.csv` rather than committed.

**What it does not clear.** The Synapse-hosted HTA11 entities are a separate
question with separate terms, unaffected by this review: `syn23520239` refuses
an authenticated read, and the Level 3 WES VCFs report `canDownload: false` for
a certified account with no unmet access requirement. Nothing here changes
either. The open cBioPortal deposit is a different distribution of overlapping
content, and only it is licensed by this section.

## 2 · Label provenance at source — condition 2

The condition is explicit that an attribute name is not evidence. Two things
were checked.

**The paper's methods.** Chen et al. state that *"Polyps were histologically
categorized by two pathologists into two subtypes: ADs consisting of tubular
ADs (TAs) and tubulovillous ADs (TVAs), or serrated polyps (SERs)"*, and that
hyperplastic polyps were further subdivided into goblet cell-rich and
microvesicular forms. `POLYP_TYPE` and `POLYP_SUBTYPE` are that reading of
H&E-stained sections. The label is not computed from the transcripts.

**One hazard, found and closed.** The same methods reclassify specimens whose
histology was *unconfirmed* using the single-cell data — a transcriptomically
assigned polyp type, which for those specimens would be exactly the circularity
invariant 11 forbids.

It does not reach the arms, and the evidence is positive rather than an absence:

- No `UNC` or unconfirmed value appears anywhere in the pinned snapshot.
  `POLYP_TYPE` takes only `AD` (31), `SER` (19), `Unknown` (11).
- **All 27 arm lesions carry a specific histologic subtype, and none is
  `Unknown`** — 13 Tubular Adenoma, 3 Tubulovillous Adenoma, 4 Sessile Serrated
  Lesion, 3 Microvesicular Hyperplastic Polyp, 2 Hyperplastic Polyp, 2 Goblet
  cell-rich Hyperplastic Polyp. That granularity is a morphological reading;
  a transcriptomic reclassification does not return "goblet cell-rich".
- The 11 `Unknown` specimens carry `Unknown` at subtype as well, consistent with
  the deposit never writing the transcriptomic reclassification back into the
  field. Rule 3 excludes them regardless.

**The guard could not have caught this, and that is the point of the
condition.** `check_no_circular_claim` short-circuits on an empty gene overlap,
so a transcript-derived label declaring no genes passes it silently. The
protection is the source check above and rule 3's exclusion, not the assertion.
`tests/test_chen_subtype_labels.py` pins that limitation so it is not later
mistaken for coverage.

## 3 · Snapshot and identifiers — conditions 3 and 4

`data/raw/cbioportal_hta11/`, checksummed in `data/manifest.csv`, written by
`src/reference/jobs/chen_subtype_snapshot.py`:
`study.json` 1,841 B; `clinical_sample.json` 561,950 B (1,067 records, 61
samples, **52 patients**); `mutations.json` 31,212,971 B (39,235 calls, 30
biospecimens). Checksums were identical across a dry run and the written run.

`src/reference/jobs/chen_subtype_id_provenance.py` joins on full-string
equality only, validated one-to-one. **Four prefix families hold more than one
specimen** — `HTA11_866` carries both `_2000001011` and `_3004761011` — so a
prefix match would merge distinct lesions and inflate an arm.

Everything §3 froze reproduces (`results/2026-09-10_461d20b/`): AD 16 labelled
/ 15 analysable / 13 patients, SER 11 / 10 / 9, the conflicted patient
`HTA11_6801`, and the concordant multi-lesion patients `HTA11_6818`,
`HTA11_8622`, `HTA11_866`.

## 4 · Label families declared — condition 5

`src/reference/chen_subtype_labels.py` declares both against one fixed claim,
in one file, because the question a reviewer must answer is whether the two are
independent of the claim in *different* ways — which is invisible when the
declarations sit in separate jobs. Family P is morphology; family G is genotype
over `APC`/`BRAF`, neither of which is a panel gene (invariant 2). Both pass
with an empty overlap.

**What neither declaration establishes.** Invariant 11 tests provenance, not
whether a label is independent of the *endpoint*. Family P is a reading of
tissue architecture and the compositional term measures tissue architecture in
another modality; the guard cannot see that and will not complain. §4 of the
pre-registration is what handles it, by giving the two terms different standing.

## 5 · Attrition — condition 6

`results/2026-09-10_dea2021/`, produced and inspected **before** any arm
contrast, which is what the condition's ordering is for. D2 §6a is the standing
reason to look: there the dropped participants' GUCA2A ran *lower* than the
retained set's, in the direction of the outcome.

**Coverage.** 23 of 44 patients covered, 21 not — 8 carrying `Unknown` and 13
with no cBioPortal record at all. The conflicted patient counts as covered:
they are reached by the label and removed by rule 1, which is not attrition.

**The comparison does not show a systematic difference.** No standardised
difference anywhere in the table reaches 0.4. The largest is
`KRT8`/compositional/normal at −0.391; the rest of the top five run
`CDX2`/intrinsic +0.385, `KRT8`/compositional/doubly_robust −0.377,
`GUCA2A`/compositional/tumour −0.360, `CDX2`/intrinsic +0.348.

**On the primary term — intrinsic, `lineage`, doubly robust — GUCA2A is the
flattest gene in the panel**: covered −3.344 against uncovered −3.384, a
standardised difference of **0.007**. The gene the design turns on is, as far
as this can show, unrelated to whether a patient carries a label. That is the
opposite of the D2 pattern and it is the reason the contrast does not need to
be reported as conditional on coverage.

**One asymmetry, and it is worth stating plainly.** All 18 not-estimable
intrinsic values — one patient across six genes and three weightings — fall in
the **uncovered** group. Every covered patient has an estimable intrinsic term
at `lineage`. Two consequences: the uncovered comparison group is 20 estimable
patients rather than 21, and both arms will carry a not-estimable count of zero
when §5 reports them (invariant 1 still requires the count be shown).

**No threshold was applied**, because none was pre-committed. Inventing one
after seeing the table is the choice this project refuses elsewhere, so the
magnitudes are reported and the judgement is recorded here as a judgement.

## 6 · What remains
- **Condition 7 — the interval departure.** W2 review was **waived by the
  project owner on 2026-09-10**, recorded here rather than dropped. The
  requirement was written because `HANDOFF` §6k records three harness changes
  landing without one; this is a fourth exception by explicit decision, not an
  oversight, and CONTRIBUTING §2–3 still says what it says.

## Sources

- cBioPortal user documentation, data usage and citation — <https://docs.cbioportal.org/user-guide/faq/>
- Chen et al., *Cell* 2021, PMID 34910928 — <https://pmc.ncbi.nlm.nih.gov/articles/PMC8941949/>
- Pinned study record, `data/raw/cbioportal_hta11/study.json`
