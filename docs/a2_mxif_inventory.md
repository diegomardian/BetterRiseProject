# A2 Step 1: Chen/HTAN MxIF product and crosswalk inventory

**Inventory frozen 2026-09-07. No marker outcome was read.** This is the data-
availability step for A2, the same-cohort protein extension of avenue A. It is
not B-style replication and it is not C's GUCA2A-specific spatial test.

The point of this step is to decide the analytic object before choosing a
statistic. In particular, a pixel array, a cell-intensity table, and a raw OME-
TIFF do not license the same spatial unit.

## 1. Sources checked

| source | what it says or contains | access on 2026-09-07 |
|---|---|---|
| [Chen et al., Cell 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8941949/) | COLON MAP paper; MANDO segmentation produced per-cell marker intensities, morphology, size and location; data availability names `syn23564801`, `syn23630431`, `syn23520239` and HTAN `HTA11` | article open; the three Synapse entities do not allow anonymous READ |
| [HTAN HTA11](https://data.humantumoratlas.org/center/hta11) | official Vanderbilt atlas landing page | open landing page; file export still routes through Synapse |
| [MILWRM paper](https://www.nature.com/articles/s42003-024-06281-8) | later pixel-level analysis of the same Vanderbilt colon MxIF cohort | open |
| [MILWRM Zenodo record](https://zenodo.org/records/10557593) | one `upload_data.zip`, 4,915,599,049 bytes, MD5 `9f89acceda8232b5044460e23c5671cc` | open; ZIP directory and one 1.5 MB processed member inspected by HTTP range, not a full download |
| [MILWRM Supplementary Information](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs42003-024-06281-8/MediaObjects/42003_2024_6281_MOESM1_ESM.pdf) | Supplementary Table 1 gives region-level normal/tumor and AD/SSL labels | open |
| [Synapse access thread](https://www.synapse.org/Synapse%3Asyn5637528/discussion/threadId%3D11828) | an earlier request for these exact three entities was directed to the Sage service desk | open; confirms that requesting access is the supported next route |

No Synapse token or remembered login is present on this machine. No attempt was
made to bypass its access control.

## 2. The public processed product

The Zenodo archive contains **42 colon-region NPZ files from 15 HTA11
participants**. One member was fetched by an exact byte range and inspected:

```
img   float64  [height, width, 26]  scaled marker intensities
ch    unicode  [26]                 channel names
mask  float64  [height, width]      binary tissue mask
```

The 26 channels are:

```
BCATENIN CD20 CD3D CD4_ CD68 CD8 CGA COLLAGEN ERBB2 FOXP3 HLAA
LYSOZYME MUC2 NAKATPASE OLFM4 PANCK PCNA PEGFR PSTAT3 SMA SOX9
VIMENTIN GACTIN CDX2 MUC5AC DAPI
```

This covers every A2 protein marker under discussion: CDX2, OLFM4, SOX9,
MUC5AC, and DAPI. The MILWRM methods say the images were downsampled 16-fold to
**5.6 micrometres per pixel**. They are therefore a reusable, normalized pixel
product.

They are **not MANDO cell tables**. `mask` is a tissue mask, not a labelled cell
segmentation; there are no cell IDs, per-cell intensities, nuclei, crypt
positions, or cell boundaries in an NPZ. A crypt-position/per-cell CDX2
equivalence analysis cannot be described as a join over this product.

## 3. Region metadata, including the mismatch

Supplementary Table 1 describes 38 of the 42 archive regions: 18 normal and 20
tumour. Four archive entries do not appear in that table and remain unlabelled:

- `HTA11_10167_0000_01_01_region_004`
- `HTA11_7862_0000_02_02_region_003`
- `HTA11_8920_0000_02_02_region_002`
- `HTA11_8920_0000_02_02_region_003`

They may not inherit a label from another region in the same batch. The fixed
row-level transcription is [the committed region manifest](../config/a2_mxif_regions.csv),
and its validator enforces that rule.

The 38 labelled regions represent 7 AD and 8 SSL participants. Eight
participants have at least one labelled normal and one labelled tumour region
within the same MxIF batch:

| type | participant IDs | paired n |
|---|---|---:|
| conventional adenoma (AD) | `HTA11_10623`, `HTA11_6298`, `HTA11_7862` | **3** |
| sessile serrated lesion (SSL) | `HTA11_10167`, `HTA11_4255`, `HTA11_7956`, `HTA11_8099`, `HTA11_8622` | **5** |

The remaining seven participants have labelled tumour regions only in this
public processed subset. Multiple regions do not increase the inferential n.

## 4. Avenue-A crosswalk

All 15 MxIF participant IDs match `Chen_2021_Cell` participants in the local
`VUMC_HTAN_validation` metadata exactly. `HTA11_866` is additionally present in
`VUMC_HTAN_discovery`: it is the original Chen DIS/VAL shared participant, not
an extra A2 patient or an independent replication. Eleven of the 15 validation
matches have both a polyp and healthy-normal scRNA sample; four are polyp-only
there (`HTA11_10623`, `HTA11_10711`, `HTA11_7179`, `HTA11_9341`).

This is an **exact participant-level crosswalk only**. The MxIF batches use
identifiers such as `HTA11_10623_0000_01_01`, while the ICBI scRNA sample uses
an identifier such as `HTA11_10623_2000001011`. No source inspected here proves
that these are the same biospecimen or adjacent sections. Every emitted row
therefore carries `crosswalk_level = participant_only` and
`specimen_exact = False`.

## 5. Step-1 decision

There is a real, immediately usable processed MxIF product, but it licenses a
**pixel-level A2**, not the proposed per-cell crypt-position equivalence test.
The conventional-adenoma within-batch paired n in this public subset is only
three. That is useful for pipeline validation and pictures; it is not a
confirmatory Student-t equivalence sample.

Therefore:

1. Do not compute CDX2 retention from these arrays yet. Choosing a pixel-domain
   statistic after learning that cell tables are absent would change the
   estimand.
2. Request Synapse access and inventory filenames/annotations first, without
   downloading raw images. The object sought is a MANDO-derived table or cell
   segmentation carrying cell coordinates and per-cell marker intensities.
3. If no such object exists, decide explicitly whether A2 becomes a new image-
   segmentation analysis. That is where ML enters as implementation.
4. Only after that choice, pre-commit the spatial unit, positive control,
   biological equivalence margin, patient aggregation, and treatment of AD
   versus SSL.

The authenticated next step is specified in
[A2 Synapse access and metadata inventory](a2_synapse_access_request.md). Its
job requests entity metadata only and refuses to treat inaccessible descendants
as absent products.

Run the deterministic inventory locally with:

```bash
python -m src.reference.jobs.a2_mxif_inventory --no-write
```

After the code and manifest are committed, remove `--no-write` to emit the two
versioned provenance-stamped metadata tables.
