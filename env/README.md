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

## MKL 2026 breaks numpy 1.26 — read this before creating an env

**Symptom:** any BLAS or LAPACK call dies with a native Windows error
`0xC06D007F` and no Python traceback. `numpy.linalg.lstsq`, plain `matmul`,
`scipy.optimize.nnls` and `sklearn`'s KMeans all crash the interpreter;
`scipy.stats` and pandas are fine, so it looks like a broken package rather than
a broken environment. `threadpoolctl.threadpool_info()` returning `[]` is the
giveaway — numpy cannot load its BLAS at all.

**Cause:** conda-forge now resolves `mkl` to 2026.x, which is ABI-incompatible
with the `numpy=1.26.4` pinned in [`pyproject.toml`](../pyproject.toml).
`0xC06D007F` is the delay-load "module not found" code. Nothing in this repo is
at fault, and it was not reproducible until the envs were built fresh in
August 2026.

**Fix — per environment, until someone pins it properly:**

```bash
conda install -n brp-wN "libblas=*=*openblas"
```

Verify:

```bash
python -c "import numpy as np; print(np.linalg.lstsq(np.eye(3), np.ones(3), rcond=None)[0])"
```

**This affects all four env files, not just W3.** Anyone creating an env today
gets MKL 2026 and hits it. The durable fix is a `libblas=*=*openblas` line in
each `wN_*.yml`, or bumping numpy — but that is a `pyproject.toml` change
touching all four workstreams, so it wants a decision rather than a drive-by.
Raised by W3, 2026-08-18, after it cost most of a day.

**It is not a code bug.** W2's NNLS baseline was suspected and is innocent; with
OpenBLAS the full suite passes, 408 tests.
