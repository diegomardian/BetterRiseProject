## What and why

<!-- One paragraph. What changed, and what question it answers. -->

**Workstream:** W1 / W2 / W3 / W4 / shared
**Week / task from execution_plan.md:**

---

## Checks

- [ ] `pytest` and `ruff check src tests` pass locally
- [ ] I only touched my own workstream's directory (or this is a `shared/` PR — see below)
- [ ] Any new data file has a row in `data/manifest.csv` with a sha256
- [ ] Any new result was written with `src.schema.write_results()` and carries a sha + seed
- [ ] Random seed fixed and logged, not left to a library default

## Invariants (tick the ones this PR touches)

- [ ] **1 · `None` is not `0.0`** — no unestimable intrinsic term written as zero
- [ ] **2 · No target-gene leakage** into labels or the reference matrix
- [ ] **4 · Per-study estimates, meta-analysed** — nothing pooled across datasets
- [ ] **5 · Bootstrap over patients**, not cells
- [ ] **6 · No cell-type-specific expression imputation from bulk** — fractions only
- [ ] **7 · Interaction term reported separately**, never folded into either arm
- [ ] **9 · DSS/PFI primary**, OS secondary
- [ ] N/A — this PR touches none of them

---

## Frozen code — delete this section unless it applies

Touching `src/schema.py`, `config/panel.yaml` or `config/labeling_axes.yaml`
needs **two approvals and a written reason** (CLAUDE.md invariant 3).

**Reason:**

**What breaks downstream if this merges:**

- [ ] `tests/test_freeze.py` updated in this same PR
- [ ] Two reviewers requested
