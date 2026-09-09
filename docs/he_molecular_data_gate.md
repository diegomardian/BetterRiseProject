# H&E molecular-prediction data gate

**Status: INVENTORY IN PROGRESS — file-level linkage observed; no model, image
download, endpoint choice, or biological result is licensed.**

This is the implementation record for `NEXT_AVENUES.md` item 2b.  It is a
separate morphology-to-molecular-label project, not a validation of avenue A,
not a repair for Crowell's control limitation, and not a classifier for a
diagnosis read from the same slide.

## Fixed gate

A candidate passes only when all four conditions below are evidenced in a
versioned metadata inventory, before an encoder or model is selected.

1. Licensed, full-resolution H&E/WSI files are downloadable.
2. Every retained image has a **specimen-exact** join to an independently
   assayed molecular call.  A patient or case identifier alone does not pass.
3. The proposed molecular endpoint is named before fitting, and its prevalence,
   missingness, and complete-case count are emitted from that join.
4. The retained data support patient-held-out evaluation and, when multiple
   acquisition sites exist, site-held-out evaluation.

Failure of any condition is `NO SUBSTRATE`; it does not permit replacement of
the molecular endpoint by polyp subtype, dysplasia, or another label read from
the H&E image.

## Candidate inventory, 2026-09-08

| Candidate | What is verified | Gate state | Reason / next metadata-only action |
|---|---|---|---|
| **HTA11 / COLON MAP** | The exported Files table contains **50** H&E Level-2 records (30 TIFF, 20 OME-TIFF), all labelled `CRDC-GC/SB-CGC (open access)`. It contains **174** unique Bulk-DNA Level-3 VCF biospecimens (253 VCF records), and **38/50 H&E biospecimen IDs exactly match** a VCF biospecimen ID. The Biospecimen export contains all 38 matches: **18** are `Premalignant` (14) or `Atypia - hyperplasia` (4), while **20** are `Primary` and excluded from the polyp estimand. Of **25** premalignant H&E biospecimens, 18 match a VCF and **7 do not**. | **NOT LICENSED — only in-scope lead** | Exact biospecimen linkage is evidenced, but full-resolution is not established by this metadata (Level 2 is not a resolution), candidate VCF access is `Synapse` and has not been authenticated, and the Case export records both diagnosis and site as `Not Reported` for all 18 candidates. No endpoint has been pre-specified or counted. |
| **TCGA COAD/READ** | GDC makes diagnostic slide images available and TCGA barcodes retain sample/portion lineage. | **EXCLUDED** | Carcinoma, not a polyp substrate. It may be an engineering benchmark, never evidence for this polyp estimand. A matching barcode prefix also would need an explicit portion-level crosswalk before any separate benchmark. |
| **SurGen** | Public 40× colorectal WSIs with KRAS, NRAS, BRAF, MMR/MSI, and survival annotations. | **EXCLUDED** | Carcinoma, not polyp; the public description does not establish the needed image-to-assay specimen crosswalk. It may be a future engineering benchmark only. |
| Public polyp image-only sets | H&E images and morphology labels exist. | **EXCLUDED** | No independently assayed molecular endpoint. Diagnostic-label imitation is out of scope. |

The first row is deliberately not called a pass merely because the same study
generated H&E and WES. The exact 18-specimen premalignant H&E–VCF overlap is
the right kind of linkage evidence, but the paper describes different physical
processing paths (fresh material for scRNA-seq and formalin-fixed material for
diagnosis/WES). The exported Case records do not establish the retained
diagnosis or site structure, and an endpoint must be fixed before any VCF is
read.

## Required inventory fields

The eventual input table must contain, at minimum:

```text
candidate_id, patient_id, biospecimen_id, image_file_id, image_synapse_id,
image_data_file_id, image_access, image_format, image_resolution_um_per_px,
molecular_file_id, molecular_access, molecular_assay,
molecular_biospecimen_id, endpoint_name, endpoint_value, endpoint_callable,
acquisition_site
```

The gate derives, rather than assumes: unique patients, unique sites,
image count, image–molecular intersection on `biospecimen_id`, endpoint
prevalence among complete cases, endpoint missingness, and the proposed split
eligibility.  It must retain an attrition row for every exclusion.

## Next action

The first Files-table export resolved the enforced assay vocabulary: `H&E`
at `Level 2`, `Bulk DNA` at `Level 3` in `vcf` format, image access labelled
`CRDC-GC/SB-CGC (open access)`, and molecular access labelled `Synapse`.
It does **not** report microns-per-pixel, so it cannot establish that the H&E
files satisfy the full-resolution clause. `Synapse` is not treated as dbGaP,
but it is also not assumed to be downloadable until an authenticated metadata
or file request succeeds.

The Biospecimen export resolved the first disease-scope split: 18
premalignant/atypia-hyperplasia exact matches are candidates, whereas 20
primary matches are excluded; 7 further premalignant H&E biospecimens have no
exact VCF. The Case export is present for all 18 candidates but supplies neither
diagnosis nor acquisition site. No image file is downloaded. The result is
therefore a durable `NOT LICENSED` metadata finding, not a source pass.

A later, independently documented source can reopen this gate only by providing
the missing resolution, molecular-access, clinical, and endpoint evidence. It
may not replace the molecular target with a morphology-derived diagnosis.

## Reproducible metadata gate

`src/reference/jobs/he_molecular_gate.py` implements the exact-ID join and its
attrition table. It reads only the three portal TSV exports:

```bash
python -m src.reference.jobs.he_molecular_gate \
  --files /path/to/files.tsv \
  --biospecimens /path/to/biospecimens.tsv \
  --cases /path/to/cases.tsv
```

The pre-write run on the 2026-09-08 exports reproduced the inventory above:

| H&E files | exact H&E–VCF biospecimens | premalignant H&E | premalignant exact matches | premalignant no VCF | candidate VCF access | reported candidate sites |
|---:|---:|---:|---:|---:|---|---:|
| 50 | 38 | 25 | 18 | 7 | `Synapse` | 0 |

All 18 candidate Case rows have `Primary Diagnosis` and `Site of Resection or
Biopsy` equal to `Not Reported`. The job consequently refuses to promote
missing site information into evidence of a single-site cohort. It also refuses
to turn the file table's access label into a successful molecular download or a
Level-2 designation into a resolution value. It has no endpoint argument by
design: an endpoint must be selected in a separate, locked pre-specification
before a VCF is opened. This pre-write run produced no result artifact because
the code and documentation were still uncommitted. A committed run that writes
`NOT LICENSED` persists both tables and exits **5**; `--no-write` exits 0
because it is only local inspection.

### Authenticated molecular-access determination

`src/reference/jobs/he_molecular_synapse_access.py` takes the candidate Synapse
IDs from the exact premalignant H&E–VCF join and calls
`Synapse.get(downloadFile=False)` on each one. It never downloads or opens a
VCF. Run it only after a CLI login, first with `--no-write`:

```bash
pip install -e '.[a2]'
synapse config
python -m src.reference.jobs.he_molecular_synapse_access \
  --files /path/to/files.tsv \
  --biospecimens /path/to/biospecimens.tsv \
  --cases /path/to/cases.tsv \
  --no-write
```

Use the default profile when prompted and supply a View-scope personal access
token. `synapse login` merely verifies a token; `synapse config` persists the
profile that the job's silent client login reads.

`readable` means the authenticated account can resolve metadata for that entity;
`unresolved` preserves the exact remote error and is not silently interpreted
as absence. Either outcome still leaves the resolution, clinical-metadata, and
endpoint gates closed.

After the inspection output is reviewed, rerun the access command without
`--no-write` to persist its table. Supply that exact parquet table to the main
gate to clear only the pending Synapse-*metadata* reason:

```bash
python -m src.reference.jobs.he_molecular_gate \
  --files /path/to/files.tsv \
  --biospecimens /path/to/biospecimens.tsv \
  --cases /path/to/cases.tsv \
  --synapse-access /path/to/he_molecular_synapse_access.parquet \
  --no-write
```

The main gate refuses an artifact whose candidate entities differ from its own
exact H&E–VCF crosswalk. Metadata readability is not treated as permission to
download or open VCF content.

### Native H&E header check

`src/reference/jobs/he_molecular_image_headers.py` is the next bounded check.
For each of the 18 exact premalignant H&E candidates, it obtains a signed URL
for the already-open image entity and requests only bytes `0–65535`. It reports
the first TIFF directory's pixel dimensions and microns-per-pixel if the needed
tags occur in that prefix. It does not download a slide or assert that a native
TIFF is the scanner's original full-resolution object.

```bash
python -m src.reference.jobs.he_molecular_image_headers \
  --files /path/to/files.tsv \
  --biospecimens /path/to/biospecimens.tsv \
  --cases /path/to/cases.tsv \
  --no-write
```

Only an explicit native-scale result plus provenance that the released Level-2
object is the required full-resolution image can resolve condition 1. A missing
or unparseable tag is retained as `unresolved`, never filled with a default.

## Sources consulted

* HTAN's portal guide documents that Vanderbilt exposes H&E and Bulk DNA as
  assay categories, and that file, case, and biospecimen metadata can be
  explored/exported: <https://docs.humantumoratlas.org/data_access/portal/>.
* Chen et al. describe the COLON MAP polyp collection, FFPE tumour WES, and
  HTAN/Synapse release: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8941949/>.
* GDC's barcode documentation distinguishes sample, portion, aliquot, and
  slide identifiers: <https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/>.
* The SurGen data note describes its public colorectal WSI and biomarker
  collection: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12569769/>.
