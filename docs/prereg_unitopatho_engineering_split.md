# UniToPatho engineering benchmark — task and split contract

**Written and locked:** 2026-09-10, before image pixels, embeddings, model
selection, calibration, or test labels were read. **Scope:** engineering-only;
this is not a biological, molecular, clinical, or diagnostic claim.

## Dataset boundary

The only permitted source is the reconciled `deephealth-uc2-800` UniToPatho
configuration. It has 8,669 image paths from 292 source-derived slide groups.
`7000` and `7000_224` are prohibited: they are separate source configurations
and their overlapping slide groups have incompatible supplied split assignments.

IEEE DataPort's access language was reviewed on 2026-09-10. It states that
datasets are CC BY and may be copied, analyzed, or otherwise used with
attribution, clearing data licensing for this engineering benchmark. This is
not a biological-data license because the dataset carries no molecular endpoint.
It also does not select an encoder or permit the protocol to be tuned after
features are read.

## Locked prediction task

The task is the source's exact six-class morphology label:

```text
HP, NORM, TA.HG, TA.LG, TVA.HG, TVA.LG
```

No binary or three-class regrouping is permitted. The label is a morphology
annotation read from the image and is not an independent endpoint.

## Locked unit and split

The inferential unit is one source-derived `slide_group`, which the source
documents as one whole-slide image per patient. Every patch from a group stays
in one model split.

1. The 88 groups supplied as `test` remain the untouched test set.
2. Every source `validation` group remains validation.
3. For each class separately, validation contains `ceil(0.20 × non-test
   groups)`. Additional groups come only from source `training` and are selected
   by ascending SHA-256 of `unitopatho-development-v1:{slide_group}`.
4. All other source-training groups remain training.

The resulting fixed group counts are 161 training, 43 validation, and 88 test.
Validation group allocation is HP 7, NORM 3, TA.HG 4, TA.LG 20, TVA.HG 3, and
TVA.LG 6. The source test set is not used to select the task, split, encoder,
calibration, abstention threshold, or any hyperparameter.

## Execution boundary

`src/reference/jobs/unitopatho_development_split.py` produces the durable
metadata-only assignment after the manifest is SHA-256 checked. It is not a
modeling authorization. Before any feature extraction, a subsequent amendment
must name one frozen encoder and its immutable revision, the linear-head fitting
procedure, calibration method, abstention score and validation-only threshold
rule, and the final held-out test metrics.
