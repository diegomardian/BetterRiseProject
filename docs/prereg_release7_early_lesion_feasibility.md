# Release-7 HTAN source-identity audit

**Written and run:** 2026-09-08 · **Status:** SOURCE IDENTITY GATE FAILED —
the supplied Release-7 participant export overlaps cached Chen completely.
**Scope:** determine whether two Release-7 HTAN Vanderbilt epithelial H5ADs are
an independent normal-to-early-lesion cohort. This audit reads no Release-7
H5AD, Synapse metadata, or expression layer.

## 1 · Fixed source and source-identity rule

The only candidate count objects are HTAN Release 7 Level-4 epithelial H5ADs:

| arm | HTAN data-file ID | Synapse entity | expected participants |
|---|---|---|---:|
| Discovery | `HTA11_0_14002` | `syn53710088` | 30 |
| Validation | `HTA11_0_14004` | `syn53710094` | 26 |

The candidate source universe is the 55 distinct `HTAN_Participant_ID` values
returned by the fixed BigQuery query in §4. It would be independent of the
cached `Chen_2021_Cell` object only if all 55 values were absent from that
object's full participant set, including its carcinoma arm. The cached set is
defined by stripping the exact `Chen_2021_Cell.` prefix from `patient_id` for
all rows whose `study_id` is `Chen_2021_Cell`; the earlier 44-patient NewCo
proxy is not an independence test.

The identity job requires exactly 55 nonmissing, unique exported IDs and the
verified 106 cached Chen IDs. An overlap writes an `OVERLAPS CACHED CHEN`
artifact and exits with status 5. It is a source-identity result, not evidence
about expression.

**Result:** the supplied BigQuery export contains 55 unique IDs. All 55 overlap
the cached 106-participant Chen universe. The source therefore fails the
independence rule. This is a cohort re-release or reprocessing candidate, not
an independent replication.

## 2 · Consequence

Release 7 cannot strengthen Avenue A by independent replication. The reported
normal/early-lesion pairing counts may be useful only for a separately
pre-registered **same-cohort reprocessing** question that identifies what is
new about the Release-7 objects and avoids treating correlated participants as
new evidence. That is not designed or authorized here.

The audit does not license H5AD access, a normal-to-lesion contrast, a
three-stage claim, an AD-versus-SSL claim, a target-gene direction, or pooled
inference. It does not establish that the H5ADs are unusable for every purpose;
it establishes that they are not an independent cohort.

## 3 · Fixed BigQuery query for the source-ID export

```sql
SELECT DISTINCT HTAN_Participant_ID
FROM `isb-cgc-bq.HTAN_versioned.id_provenance_r7`
WHERE HTAN_Center = 'HTAN Vanderbilt'
  AND HTAN_Data_File_ID IN ('HTA11_0_14002', 'HTA11_0_14004')
ORDER BY HTAN_Participant_ID;
```
