# UniToPatho engineering benchmark — data gate

**Status:** manifest structure audited; the 48.37-GiB selected archive is held
outside the repository in the researcher's Google Drive. Dataset terms were
reviewed on 2026-09-10: the IEEE DataPort access language states that datasets
are available under CC BY and may be copied, analyzed, or used for other
purposes with attribution. This clears data licensing for the engineering
benchmark. It does not turn its morphology labels into a biological endpoint or
replace the remaining pre-feature-extraction protocol amendment.

## Scope

This is WP-B in [`ml_abstention_roadmap.md`](ml_abstention_roadmap.md): an
engineering-only benchmark for a frozen H&E encoder and validation-fixed
abstention rule. It cannot be promoted to a biological project result. Its
labels are diagnoses read from the classified tissue and it carries no
specimen-exact independent molecular endpoint.

The selected configuration is the source's `deephealth-uc2-800` manifest,
`unitopath-public-800.yml`, rather than a merger of all published image
configurations. `800` is the source configuration name, not the delivered PNG
width: the selected patch paths specify 1,812x1,812-pixel crops. The source
README describes 9,536 H&E patches from 292 patient-distinct whole-slide
images at 0.4415 microns/pixel. The selected configuration contains 8,669
patches and all 292 source-derived slide keys.
The separate `deephealth-uc2-7000` configuration overlaps 243 slide keys and
assigns 17 of them differently between its supplied splits. Combining them
would create unreviewable leakage, so it is prohibited.

The source's repository is the provenance record for the public manifest and
the dataset's suggested citation is Barbano et al., *ICIP* 2021,
doi:10.1109/ICIP42928.2021.9506198.

## Gate

`src/reference/jobs/unitopatho_feasibility.py` reads only the downloaded
800-pixel YAML manifest. It rejects a changed SHA-256, an unknown class,
unassigned or multiply assigned patch indices, path/label disagreement, or a
derived slide key that crosses training/validation/test. It emits counts by
split and class and an explicit summary that the work remains engineering-only.

Run the gate only after downloading the public YAML manifest from the official
repository; it does not need the dataset images:

```bash
curl -L --fail https://raw.githubusercontent.com/EIDOSLAB/UNITOPATHO/main/unitopath-public-800.yml \
  -o "$HOME/Downloads/unitopath-public-800.yml"
python -m src.reference.jobs.unitopatho_feasibility \
  --manifest "$HOME/Downloads/unitopath-public-800.yml"
```

## Conditions still required before modeling

1. **Complete — retain the evidence.** Save the IEEE DataPort access-language
   record alongside local access notes. It states that datasets are CC BY and
   may be copied, analyzed, or otherwise used with attribution. Cite the
   dataset and the Barbano et al. source paper in any output. The GitHub
   repository's MIT code license is not the data license.
2. Record the dataset version, archive checksum, and image-to-manifest
   reconciliation without changing the selected 800-pixel configuration. The
   Google Drive listing itself is enough for this structural check; it does not
   open an image. On the cluster, after configuring a read-only `rclone`
   remote, run:

   ```bash
   mkdir -p "$BRP_DATA_DIR/interim/unitopatho"
   curl -L --fail \
     https://raw.githubusercontent.com/EIDOSLAB/UNITOPATHO/main/unitopath-public-800.yml \
     -o "$BRP_DATA_DIR/interim/unitopatho/unitopath-public-800.yml"
   rclone lsf -R "Meowsers:800" \
     > "$BRP_DATA_DIR/interim/unitopatho/unitopatho_800_listing.txt"
   python -m src.reference.jobs.unitopatho_archive_reconciliation \
     --manifest "$BRP_DATA_DIR/interim/unitopatho/unitopath-public-800.yml" \
     --listing "$BRP_DATA_DIR/interim/unitopatho/unitopatho_800_listing.txt"
   ```

   The expected listing has 8,671 objects: 8,669 manifest image paths plus
   `train.csv` and `test.csv`. Never mix the distinct `7000` configuration into
   this listing. Capture the shown `rclone lsf -R "Meowsers:800"` command
   verbatim: it supplies bare root-relative file paths, while any directory
   markers are ignored by the gate.
3. Lock one task and a patient-contained train/validation/test protocol before
   feature extraction. The task and split are locked below; a subsequent,
   pre-feature-extraction amendment must name the frozen encoder and evaluation
   procedure. The encoder remains frozen; all threshold selection is
   validation-only; the test set is opened once for coverage, selective error,
   calibration, and the ordinary non-abstaining baseline.

   The metadata-only task/split contract is now
   [`prereg_unitopatho_engineering_split.md`](prereg_unitopatho_engineering_split.md).
   After this code is committed, run:

   ```bash
   python -m src.reference.jobs.unitopatho_development_split \
     --manifest "$BRP_DATA_DIR/interim/unitopatho/unitopath-public-800.yml"
   ```

   This writes the deterministic 161/43/88 patient-contained assignment. It
   does not authorize or perform image access.

Failure of condition 1 stops before image download or modeling. It is not a
reason to substitute MHIST or to claim an independent molecular result.
