# env/

One conda env per workstream, pinned. Four people, four sets of heavy
dependencies — CellBender wants a GPU and a specific torch, `lifelines` does
not, and neither should be able to break the other.

```bash
conda env create -f env/w2_harness.yml
conda activate brp-w2
pip install -e ".[dev]"      # from the repo root — installs src/ so imports work
```

| File | Env name | Workstream |
|------|----------|------------|
| `w1_reference.yml` | `brp-w1` | scanpy, scDblFinder/scrublet, SoupX, CellBender, inferCNV |
| `w2_harness.yml` | `brp-w2` | scanpy, scikit-learn, the deconvolution bake-off |
| `w3_bulk.yml` | `brp-w3` | pandas, lifelines, GDC tooling. Laptop-fine. |
| `w4_estimator.yml` | `brp-w4` | scanpy, statsmodels, bootstrap/hierarchical models |

## Rules

- **Pin what you add.** `package=1.2.3`, not `package`. An unpinned env is a
  result that reproduces on your machine and nowhere else.
- Add the dependency to *your* env file only. If two workstreams need the same
  new package, they both pin it — duplication is cheaper than coupling.
- Anything that must be identical across all four (numpy, pandas, pyarrow) is
  pinned in [`pyproject.toml`](../pyproject.toml) instead, because the shared
  code in `src/schema.py` and `src/common/` runs in all four envs.
- Export after you change an env, so the file matches reality:
  ```bash
  conda env export --no-builds --from-history > env/w2_harness.yml
  ```

## R

W1 and W4 may need R for inferCNV, CopyKAT and SoupX. Those are in the conda
envs via `r-base` and bioconda. If you end up managing R separately, say so in
the weekly meeting — a second, undeclared environment is exactly how a result
stops reproducing.

## Compute

32–64 GB RAM and GPU access for the single-cell arm (CellBender, inferCNV on
371k cells). The bulk/TCGA arm is laptop-fine. "Standard consumer laptop" was
wrong for the single-cell arm and is corrected in both documents. **Confirm the
machine before week 1** — it is the only hard infrastructure dependency.
