# notebooks/

One folder per workstream: `w1/ w2/ w3/ w4/`. Name them so the order is obvious
and the owner is obvious: `w3_02_premise_check_guca2a.ipynb`.

**Notebooks are for looking, not for producing.** Anything that generates a
result another workstream depends on — an S matrix, the bulk matrix, a
decomposition — belongs in `src/`, imported by the notebook, not defined in it.
A function that only exists in a notebook cannot be tested, cannot be reviewed
line-by-line in a PR, and will not survive the week-5 gate.

```python
from src.common import panel_genes, set_global_seeds
from src.schema import write_results

set_global_seeds(20260815)
```

Results go through `write_results()` even from a notebook, so they carry a git
sha and a seed like everything else.

## Before committing one

- **Clear outputs.** `jupyter nbconvert --clear-output --inplace notebooks/w2/foo.ipynb`
  Committed outputs make the diff unreadable and can leak patient-level data
  into git, where it is very hard to remove.
- Notebooks are marked `-diff` in `.gitattributes` — review the code you moved
  into `src/`, not the JSON.
