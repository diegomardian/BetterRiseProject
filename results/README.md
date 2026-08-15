# results/

Versioned parquet. **These are committed** — they are small and they are the
record of what we found.

```
results/2026-08-15_a1b2c3d/
    lee_decomposition.parquet
    lee_decomposition.meta.json      git sha, branch, dirty flag, seed, versions, row counts
```

## Write them through the writer, always

```python
from src.schema import empty_results_frame, write_results

df = empty_results_frame()
# ... fill it ...
write_results(df, name="lee_decomposition", seed=20260815)
```

Never `df.to_parquet()`. The writer is what:

- validates against the frozen schema in [`src/schema.py`](../src/schema.py)
- **asserts that no `not_estimable` row carries a numeric intrinsic term** —
  `None` is not `0.0`, CLAUDE.md invariant 1
- stamps the git sha and the seed onto the sidecar — invariant 10
- refuses to write from a dirty working tree, because a sha that does not
  reproduce the table is worse than no sha. Pass `allow_dirty=True` only for
  scratch runs that nobody else will see.

## Reading someone else's result

```python
from src.schema import read_results
df = read_results("results/2026-08-15_a1b2c3d/lee_decomposition.parquet")
```

Validates on the way in, so a table that drifted from the contract fails loudly
at your desk rather than quietly in the figure.

## Naming

`{cohort}_{what}.parquet` — `lee_decomposition`, `tcga_fractions`,
`harness_attenuation_curve`, `pelka_decomposition_lineage`. The date and sha are
in the directory name; do not repeat them in the filename.

## One result per gate criterion

By week 5 there should be a file here answering each of G1–G4 with numbers, not
a slide. The gate is decided against pre-committed criteria — the criteria are
in [execution_plan.md §5](../execution_plan.md#5-week-5-gate) and they were
written before the data.
