# config/gene_index/

The fixed gene index. **This is what makes integration a join and not a
negotiation** (CLAUDE.md, repo layout).

```
gene_index_1.0.0.txt      one identifier per line, ordered, no header
gene_index_1.0.0.map.tsv  identifier -> symbol, plus whatever else is needed
```

W1's S matrices and W3's bulk matrix are both reindexed onto the same file.
Load it through [`src/common/paths.py`](../../src/common/paths.py):

```python
from src.common.paths import gene_index_path
index = gene_index_path("1.0.0").read_text().split()
```

## Not yet produced

Two things must be settled in **week 1**, before either arm builds a matrix —
they are tracked as [open decisions #2 and #3](../../docs/open_decisions.md):

1. **Who emits it.** Recommendation: W1, from the GSE178341 feature table. W3
   conforms. Whoever it is, only one of them builds it.
2. **Symbols or Ensembl IDs.** Recommendation: unversioned Ensembl IDs as the
   index with symbols as a mapped column — TCGA STAR counts are versioned
   Ensembl, the panel and both axes are symbols, and an unmanaged mapping is the
   usual way a join silently loses ~8% of genes.

## Versioning

Bump the version and commit a new file; never edit one in place. Results
reference the index version through the sha they were written under, and an
in-place edit makes every earlier result unreproducible without saying so.

## One rule that is not negotiable

**Target genes must not be in the index used to build the reference matrix.**
`build_signature()` asserts this and will refuse to run. The index itself may
contain them for other purposes, but the call that builds an S matrix passes the
target set and gets checked. See CLAUDE.md invariant 2.
