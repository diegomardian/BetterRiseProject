# ML abstention roadmap

**Written:** 2026-09-08 · **Status:** proposed engineering branch. This is
separate from the adenoma mechanism result. It does not reopen the closed H&E
molecular-prediction gate or turn a morphology-derived diagnosis into a
biological endpoint.

## Decision

Prioritize one small, public pathology benchmark: use a frozen pathology encoder
to test whether a pre-committed selective-prediction rule abstains more often on
images where pathologists disagreed. The empirical benchmark is the priority.
A methods survey is an accompanying component, not a standalone result, because
the current repository benchmark measures simulated estimability behaviour and
does not establish what published pathology systems report.

The core question is narrower than “does uncertainty help?”:

> On a held-out image set, does a rule fixed without test labels or annotator
> agreement defer the cases that independent pathologists found most ambiguous,
> while reporting its coverage and selective error honestly?

That is a real ML result, but it is an engineering result. The clinical label is
read from the image. It therefore cannot validate intrinsic expression loss,
the Chen decomposition, or a molecular endpoint. `CONTRIBUTING.md` §6 governs
that boundary.

## Ranked sources

| rank | source | verified strengths | constraint and role |
|---:|---|---|---|
| 1 | [MHIST](https://bmirds.github.io/MHIST/) | 3,152 H&E polyp images, HP-versus-SSA majority labels, and seven-pathologist vote information for every image | The delivered research-use agreement permits non-commercial research use but prohibits modification and derivative works. Written permission is required before embeddings, model fitting, or derived-model publication. The metadata also lacks patient/WSI identifiers, so any later permitted result is image-level only. |
| 2 | [UniToPatho](https://github.com/EIDOSLAB/UNITOPATHO) | 9,536 H&E patches from 292 slides, six polyp/dysplasia classes, 0.4415 microns/pixel, and one slide per patient | Download access and supplied split/group fields need an inventory, but the source explicitly states one slide per patient. This is the first patient-contained engineering benchmark and external check after MHIST. It has no multi-pathologist agreement outcome. |
| 3 | IMP-CRS2024 | Published work describes a large colorectal biopsy/polypectomy WSI collection with a held-out test set | Access terms, downloadable files, label definitions, and split provenance are not yet verified. Treat it as a data-access investigation, not an external test set in hand. |
| excluded for this branch | NCT-CRC-HE-100K, TCGA COAD/READ, CRAG, DigestPath, LC25000 | useful colorectal pathology benchmarks | They do not supply the polyp ambiguity construct. They may later test engineering transport, but cannot substitute for MHIST or a polyp patient-level test. |

MHIST’s official page confirms the seven-reader voting field and explicitly
requires a research-use request before download. UniToPatho’s repository
confirms both the six-class polyp task and the one-slide-per-patient structure.
Neither fact licenses a biological claim.

## WP-A — MHIST agreement-aware abstention benchmark

### Feasibility gate before a model

1. Obtain and record the MHIST research-use agreement. Keep the downloaded
   manifest, checksums, annotation schema, any supplied split, and license outside
   Git. After the files are placed in the managed raw-data location, add their
   checksums and access date to `data/manifest.csv`. The delivered agreement permits
   non-commercial research use but prohibits modification and
   derivative works. Therefore checksum/schema inspection is allowed, while image
   embeddings, model fitting, and publishing derived outputs remain **LICENSE
   PENDING** until the data publisher gives written permission for those exact uses.
2. Audit whether the delivered annotations contain patient, slide, or source
   image grouping. If not, the project may report only an image-level split. It
   must not call the evaluation patient-held-out. If no split is supplied, the
   later pre-registration must fix its random split and seed before training.
3. Audit label and vote fields before choosing an analysis: each image needs a
   majority label and a valid seven-reader vote count. Missing or contradictory
   votes are retained in an attrition table, not silently imputed.
4. Audit one exact frozen encoder release for research-use availability,
   checkpoint hash, preprocessing, and hardware feasibility. Do not select the
   encoder by test performance or replace an unavailable model after seeing
   results.

### Locked design after feasibility passes

Use one frozen encoder and one linear classification head. Fit the head only on
the training split fixed after the feasibility audit. Select temperature
calibration and a single abstention threshold only on the validation split. The
test split remains untouched until those choices are written and locked.

The primary prediction target is the supplied majority HP/SSA label. Annotator
agreement is an external difficulty annotation, not a feature, training weight,
threshold input, or selection criterion. On the untouched test split report:

- ordinary accuracy and calibration without abstention;
- coverage and selective error at the fixed threshold;
- the same two selective metrics within pre-declared agreement strata; and
- a matched-coverage random-deferral reference, so a lower selective error is
  not credited merely because fewer images were answered.

The primary construct check is whether the fixed rule defers low-consensus
images more often than high-consensus images. It is not proof that deferral is
clinically correct, nor does it make the majority vote a gold-standard biology.
An outcome in the opposite direction is a useful failure of the proposed
construct, not a reason to retune the threshold.

## WP-B — UniToPatho patient-contained extension

Start only after WP-A’s data gate is documented. First inventory the actual
download: license, manifest, class counts, group key, and supplied split. Lock
one task before feature extraction; do not choose a binary grouping after
inspecting results. Keep every patch from a slide/patient in one split. Use the
same frozen-encoder, validation-only threshold, and test-only reporting rules
as WP-A.

This tests transport of abstention metrics to a physically scaled,
patient-contained polyp dataset. It cannot test the human-agreement construct
because UniToPatho has one diagnostic label per patch. It remains an engineering
benchmark.

## WP-C — survey, bounded and falsifiable

The survey should not claim that “the field does not abstain” until it measures
that proposition. Register a search protocol before screening: databases,
search terms, time window, CRC/polyp pathology scope, inclusion criteria, and
the exact coding categories below. Every included paper receives one row:

```text
paper_id, task, data source, external test, uncertainty reported,
calibration reported, abstention_or_deferral_rule, threshold_source,
coverage_reported, selective_error_reported, reviewer_notes
```

Distinguish uncertainty visualisation from a decision rule, a decision rule
from a validation-fixed threshold, and abstention from unreported low
confidence. The result is a count with an auditable denominator. It becomes a
methods companion to WP-A only if the protocol is followed; otherwise it stays
as background research.

## Existing project branches and housekeeping

The HTA11 H&E molecular-prediction branch remains **NOT LICENSED**. It needs a
specimen-exact independent molecular endpoint, credible full-resolution
evidence, endpoint prevalence, and patient/site split support before any model.
Public morphology labels cannot fill any of those gaps.

Release-7 HTAN epithelial H5ADs are not an independent early-lesion cohort:
all 55 supplied Release-7 participant IDs overlap IDs in the cached Chen
universe. The
source-identity audit in `src/reference/jobs/release7_early_lesion_identity.py`
must be committed and run from a clean tree to persist that result. It is not an
ML substrate for independent replication.

The CRDC H&E header audit is already committed. The remaining cleanup items are
to mark the superseded `he_molecular_image_headers.py` path as retired and to
resolve the repository-wide pre-existing Ruff debt in a separate hygiene change;
neither should delay WP-A’s data-access feasibility gate.

## Immediate sequence

1. Commit the Release-7 identity audit currently in the worktree and run it
   cleanly to preserve the 55/55 overlap result.
2. Obtain written MHIST permission specifically covering image embeddings,
   frozen-encoder feature extraction, fitting an academic non-commercial model,
   and publishing aggregate performance results without releasing images or
   embeddings. Do not download a model or begin training until then.
3. Once that permission is recorded, write and lock the WP-A feasibility
   inventory and model specification before opening the image pixels.
4. Execute WP-A. Only after its test artifact is complete, inventory UniToPatho
   for WP-B.
5. In parallel but not as a substitute, write the WP-C survey protocol.

## Source notes

- [MHIST official dataset page](https://bmirds.github.io/MHIST/) documents the
  image count, seven-pathologist annotation field, and research-use request.
  Cite the dataset paper as: Jerry Wei, Arief Suriawinata, Bing Ren, Xiaoying
  Liu, Mikhail Lisovsky, Louis Vaickus, Charles Brown, Michael Baker, Naofumi
  Tomita, Lorenzo Torresani, Jason Wei, and Saeed Hassanpour, “A Petri Dish for
  Histopathology Image Analysis”, *International Conference on Artificial
  Intelligence in Medicine (AIME)*, 12721:11–24, 2021.
- [UniToPatho official repository](https://github.com/EIDOSLAB/UNITOPATHO)
  documents 9,536 patches from 292 WSIs, one slide per patient, its six classes,
  physical sampling, and its IEEE DataPort route.
- [PathBench’s external evaluation](https://www.nature.com/articles/s41551-025-01516-3)
  illustrates why a frozen-encoder comparison must control data overlap and
  report external evaluation carefully.
- [UAD-FM](https://doi.org/10.1038/s41746-025-02149-1) confirms that uncertainty
  and calibration are active CRC-pathology research areas. That makes a generic
  uncertainty claim insufficient; WP-A’s agreement-based construct is the
  differentiator to test rather than assume.
