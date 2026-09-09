# Chen WES molecular-subtype gate

**Status: STEP 1 COMPLETE; STEP 2 HEADER-CAPABILITY CHECK NEXT.** The committed
crosswalk found 35 specimen-exact WES candidates across 29 of the fixed 44
lineage patients. This plan can produce a source/callability verdict. It does
not license assigning a molecular subtype or recomputing the avenue-A
transcript result until its later steps have been separately locked.

## Question and boundary

The proposed secondary question is whether the already-computed Chen avenue-A
decomposition differs between **molecular-signature-defined** polyp arms. It
is a same-cohort, modality-orthogonal heterogeneity analysis, not an
independent replication and not a claim about every conventional adenoma or
every SSL.

The transcript-side universe is the committed 44-patient, `lineage`-rung Chen
adenoma result. The molecular label, if it becomes available, is produced by
FFPE WES rather than by the transcripts used to estimate the decomposition.
The two measurements are declared separately under invariant 11. Invariant 4
forbids pooling this within-study subgroup result with another cohort;
invariant 7 requires a subtype-by-term interaction to remain separate from its
compositional and intrinsic terms.

Chen et al. reported APC alterations in conventional adenomas/TVAs and BRAF
V600E enrichment in SSLs. That is the biological premise for this gate, **not**
an assignment rule or a prediction of the retained arm sizes. [Chen et al.,
*Cell* 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8941949/).

## Step 1 — specimen-exact provenance crosswalk

`src/reference/jobs/wes_subtype_lineage_manifest.py` first emits one row for
each deposited **polyp scRNA biospecimen** belonging to the fixed 44-patient
lineage universe. It reads only existing decomposition identifiers and cached
ICBI sample metadata; no HTAN table or VCF is accessed.

```bash
python -m src.reference.jobs.wes_subtype_lineage_manifest --no-write
```

The job retains multiple polyp samples from one participant as multiple
candidate lesions. It extracts the full deposited HTA11 biospecimen token but
does not parse or compare its numeric suffix. A participant overlap or a suffix
overlap is **not** a crosswalk.

The next, still-required input is a normalized HTAN `id_provenance_r7` export
with, for every candidate VCF file, its file ID, parent file ID, entity ID,
participant ID, assayed biospecimen ID, and originating biospecimen ID. The
crosswalk passes a lesion only if the full scRNA polyp biospecimen ID equals a
documented assayed or originating WES biospecimen ID through that chain. Each
non-match remains in the output with one of: `no_participant_vcf`,
`participant_vcf_not_specimen_exact`, `ambiguous_provenance`, or
`missing_provenance_field`.

No suffix rule, filename heuristic, or manual reconciliation can upgrade a
row. If no specimen-exact WES arm remains, the route ends **NO SUBSTRATE AT
THIS RESOLUTION**.

### Required BigQuery export

The observed `id_provenance_r7` schema supplies the needed file, participant,
assayed-biospecimen, and originating-biospecimen fields. Run this read-only
query in the user's BigQuery project, then download the results as
`htan_vanderbilt_vcf_provenance_r7.csv`:

```sql
SELECT
  HTAN_Data_File_ID AS vcf_data_file_id,
  Filename AS vcf_filename,
  entityId AS vcf_entity_id,
  HTAN_Parent_Data_File_ID AS vcf_parent_data_file_id,
  HTAN_Participant_ID AS htan_participant_id,
  HTAN_Assayed_Biospecimen_ID AS vcf_assayed_biospecimen_id,
  HTAN_Originating_Biospecimen_ID AS vcf_originating_biospecimen_id,
  Biospecimen_Path AS biospecimen_path
FROM `isb-cgc-bq.HTAN_versioned.id_provenance_r7`
WHERE HTAN_Center = 'HTAN Vanderbilt'
  AND REGEXP_CONTAINS(LOWER(Filename), r'\.vcf(?:\.gz)?$')
ORDER BY htan_participant_id, vcf_data_file_id;
```

This queries provenance metadata only. It does not authenticate to Synapse or
read any VCF content. The local reader
`src/reference/jobs/wes_subtype_provenance_crosswalk.py` rejects altered column
names, duplicated VCF file IDs, non-VCF filenames, and all participant-only or
suffix-only matches.

## Step 2 — access and technical capability, not variant outcomes

Only exact Step-1 rows may be access-probed. Authentication must establish
whether the VCF **content** can be read; the earlier `Synapse.get(...,
downloadFile=False)` result established entity-metadata readability only.

A bounded header/schema inspection may read `##` header lines and `#CHROM` to
record the reference genome, caller, FILTER/INFO/FORMAT declarations, and the
number of sample columns. It cannot infer matched-normal status from that count
or from an unrecorded sample name; that status remains explicitly unresolved
until a later locked rule can establish it. It must not read variant records.
Its sole decision is whether the source has the fields needed for a
pre-specified label. Missing fields end the route; they do not authorize
changing the label definition afterward.

`src/reference/jobs/wes_subtype_vcf_header.py` implements this boundary. It
requests at most 262,144 bytes from each exact-match entity via a Synapse
byte-range request and closes the stream as soon as `#CHROM` is encountered.
It refuses a storage response that ignores the range. Its versioned output is a
technical schema artifact, not a variant call or an endpoint screen:

```bash
python -m pip install -e '.[a2]'
synapse config
python -m src.reference.jobs.wes_subtype_vcf_header --no-write
```

A variant-only VCF cannot prove a negative call at a site it does not list.
Therefore the first analysis is positive-only unless a gVCF, coverage table, or
aligned-read source provides pre-specified locus-level negative-call evidence.

## Step 3 — lock label rules before a variant record is read

The only candidate retained arms are:

| arm | assignment mechanism |
|---|---|
| `serrated_positive` | a passing BRAF p.V600E call |
| `conventional_positive` | a passing truncating APC call |

`APC` missense-only, KRAS-only, dual-positive, neither-positive, and uncallable
lesions are **unclassified**. They are not reassigned from pathology, silently
dropped, or treated as molecularly negative. The final artifact reports every
unclassified category and its count (invariant 1).

Before any VCF body is read, the locked specification must state:

1. the exact reference build and transcript/functional-annotation source;
2. `FILTER` rule, VAF floor, depth/evidence rule, and any strand/FFPE-artifact
   rule;
3. the APC consequence classes eligible as truncating (frameshift, nonsense,
   and/or canonical splice), including treatment of homopolymer-associated
   calls;
4. matched-normal handling and the wording permitted for tumor-only calls; and
5. a numeric minimum-arm and precision rule in estimator units.

Those numbers are intentionally not filled in here: the header capability
inspection must first establish which fields exist. This document is a plan,
not a retroactive preregistration.

## Step 4 — outcome-blind callability gate

After the Step-3 specification is locked, read only enough WES content to
assign `serrated_positive`, `conventional_positive`, or a recorded
unclassified reason. Emit specimen-level attrition, patient-level arm counts,
content-access status, matched-normal status, and all callability reasons.

If either molecular-positive arm fails the precommitted size/precision rule,
write **NO SUBSTRATE AT THIS RESOLUTION**. This is a measurement limitation;
it is not evidence that the biological subtypes are equivalent.

## Step 5 — separately locked subtype-stratification analysis

Only a Step-4 pass licenses a new pre-registration for a patient-level,
same-cohort subtype-by-term comparison of the existing decomposition. It must
state the exact estimand, controls, patient aggregation, Student-t interval,
and the rule for multiple terms before an arm-specific transcript value is
calculated. It will be reported as a modality-orthogonal, same-cohort
secondary result—not as independent replication or universal AD-versus-SSL
biology.
