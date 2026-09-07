# Stage 4 carcinoma resolvability audit

## Scope

This is an audit of the committed `best4` rows in
`results/2026-09-05_d358109/stage4_fractions.parquet`. It does not estimate a
biological mature-colonocyte fraction and it does not reopen the locked Stage 4
analysis.

## Recomputed denominator

The `best4` subset contains 1,350 rows from 675 unique TCGA sample IDs, with
two deconvolution methods per sample.

| quantity | count | denominator | share |
| --- | ---: | ---: | ---: |
| Exactly-zero method-by-sample estimates | 1,314 | 1,350 rows | 97.333% |
| Tumours zero under both methods | 641 | 675 tumours | 94.963% |
| Tumours with any nonzero method estimate | 34 | 675 tumours | 5.037% |

The often-shorthand phrase “97.3% of tumours” is therefore false. The 97.3%
quantity is a method-by-sample-row rate.

## Estimability inconsistency

All 1,350 `best4` rows have `estimability = 'estimated'` and an empty
`estimability_reason`, including the 1,314 zero estimates. This conflicts with
the metadata's stated invariant that a rung without a maturity call should be
`None` and `not_estimable`, rather than `0.0`.

The zero concentration is an artifact-level resolvability diagnostic. It is not
a measurement of tumour biology and cannot support a claim about the prevalence
of mature colonocytes in carcinoma. Any future quotation must name the source
artifact, distinguish row from tumour denominators, and explain the mismatch
between the numeric estimate and its estimability fields.
