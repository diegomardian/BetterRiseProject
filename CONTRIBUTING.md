# Contributing

How four people work in this repo in parallel without stepping on each other.
Read [CLAUDE.md](CLAUDE.md) for the invariants and
[execution_plan.md](execution_plan.md) for who does what when.

---

## 1. Setup, once

```bash
git clone https://github.com/diegomardian/BetterRiseProject.git
cd BetterRiseProject

# Your workstream's env. W1=reference, W2=harness, W3=bulk, W4=estimator.
conda env create -f env/w2_harness.yml
conda activate brp-w2

# Installs src/ as an editable package so `from src.schema import ...` works
# from anywhere, including notebooks.
pip install -e ".[dev]"

pytest          # must be green before you write anything
```

If `pytest` is red on a fresh clone, that is a repo bug — say so in the weekly
meeting, do not work around it.

---

## 2. Where your code goes

| You are | You own | Do not edit |
|---------|---------|-------------|
| **W1** Reference | `src/reference/` | `src/harness/` `src/bulk/` `src/estimator/` |
| **W2** Method | `src/harness/` | `src/reference/` `src/bulk/` `src/estimator/` |
| **W3** Bulk & clinical | `src/bulk/` | `src/reference/` `src/harness/` `src/estimator/` |
| **W4** Estimator | `src/estimator/` | `src/reference/` `src/harness/` `src/bulk/` |

Everyone reads `src/schema.py` and `src/common/`. **Nobody edits them alone** —
see §5.

Need something from another workstream? Don't reach into their module and don't
copy their code. Ask for it to be exposed in their `__init__.py`, or take the
artifact (parquet) instead of the function. Interfaces are files on disk, not
Python imports.

---

## 3. Branches, commits, PRs

**Branch name carries the workstream.** This is how the team sees what everyone
is working on at a glance in the GitHub branch list:

```
w1/ingest-gse178341
w2/pseudobulk-generator
w3/tcga-gdc-ingest
w4/kitagawa-unit-tests
shared/schema-add-field      # touches frozen code — needs 2 approvals
```

```bash
git checkout main && git pull
git checkout -b w2/pseudobulk-generator
# ... work ...
git push -u origin w2/pseudobulk-generator
gh pr create --fill
```

Rules:

- **Never push to `main`.** All changes land through a PR.
- Push your branch **daily**, even mid-work. An unpushed branch is invisible to
  the other three, and week-5 depends on nobody being surprised.
- Small PRs. A week-long branch is a merge conflict with a delay fuse.
- Rebase on `main` before asking for review: `git fetch && git rebase origin/main`.
- One approval to merge into `main` — **two** if the PR touches frozen code (§5).

---

## 4. Data never goes in git

`data/` is gitignored in full. What travels between machines is
[`data/manifest.csv`](data/manifest.csv): one row per file, with a sha256 and
where it came from. If you download something, add the row in the same PR as the
code that reads it. Read [data/README.md](data/README.md) before your first
download.

Results *do* go in git. They are parquet, they are small, and they are the
record of what we found. Write them with `src.schema.write_results()` — never
`df.to_parquet()` directly — so every file carries a git sha and a seed. See
[results/README.md](results/README.md).

---

## 5. Frozen code — PR + two approvals + a written reason

| Frozen | Why |
|--------|-----|
| `src/schema.py` | The output contract. Invariant 1 lives in its writer. |
| `config/panel.yaml` | The panel is a control set. Changing it retro-fits the falsification rule. |
| `config/labeling_axes.yaml` | Labels were frozen *before* the panel expanded, on purpose. |

Each is mirrored by an exact-value assertion in `tests/test_freeze.py`. Editing
the config alone turns CI red — that is the friction working, not a broken test.
Change both files, in one PR, titled `shared/...`, with the reason in the PR body.

---

## 6. Before you open a PR

```bash
pytest
ruff check src tests
```

**Your local env is not CI's env.** CI installs only `pip install -e ".[dev]"` —
the four thin shared deps plus the test tooling — not your workstream's conda
env. So a module that imports `scipy`, `statsmodels` or `scikit-learn` will
import fine on your machine and kill the CI run at *collection*, taking the whole
suite with it.

If you add such an import and a test touches that module, add the package to
`[project.optional-dependencies] dev` in `pyproject.toml` when it is light and
pip-installable. Do not add the heavy ones (CellBender, inferCNV, anything from
`r-base`) and do not reach for `pytest.importorskip` — that converts a missing
dependency into silently untested code, which is worse than a red build.

To check the way CI will see it:

```bash
python -m venv /tmp/ci && /tmp/ci/bin/pip install -e ".[dev]" && /tmp/ci/bin/pytest -q
```

- [ ] Tests pass, including `tests/test_freeze.py`
- [ ] No `intrinsic = 0.0` where the honest answer is `None` (invariant 1)
- [ ] Any new result written through `write_results()`, so it carries a sha
- [ ] Random seed fixed and logged, not left to the default
- [ ] New data files have a `data/manifest.csv` row
- [ ] Branch is `wN/...` or `shared/...`
- [ ] **Adversarially audited, if this PR publishes a result table or quotes a
      gate verdict** (decision #24.4). Brief an independent agent to *refute*
      the claims, not to confirm them, and act on what it returns.

**Why that last one is not optional.** One audit found a result table already
published from a statistic computed over the wrong population — the dedup line
that caused it cited invariant 5 and issue #36 by name and looked correct to the
person who wrote it. Ad-hoc review catches the errors you suspect. This project's
recurring failure is the errors that look right, and those need someone whose job
is to disagree.

---

## 7. Cadence

Per [execution_plan.md §11](execution_plan.md#11-standing-meetings): weekly
30 min, blockers only. Week 2 is the W3 premise check and the W1→W2 pilot
handoff. Week 5 is the gate. If your handoff is going to slip, say it in week 1,
not week 4 — §8.3 has a pre-agreed answer for a slipped pilot (move the gate to
week 7; do not compress W2).
