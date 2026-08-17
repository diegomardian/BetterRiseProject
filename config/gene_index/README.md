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

## What exists: `0.9.0`, provisional, built by W3

W1's week-1 work landed without the index, so the fallback written into
[open decision #2](../../docs/open_decisions.md) fired and W3 built one from the
GDC gene model:

```bash
python -m src.bulk.ingest gene-index
```

[`src/bulk/gene_index.py`](../../src/bulk/gene_index.py) emits both files. The
version is **0.9.0 and not 1.0.0 on purpose** — provisional, superseded the
moment W1 emits theirs. Whether W1 conforms to it or replaces it is a weekly
decision; what must not happen is two indexes existing.

**Key:** unversioned Ensembl ID, per [decision #3](../../docs/open_decisions.md).
The version suffix lives in its own column, because the same gene carries a
different suffix in a different GENCODE release and a versioned key drops those
genes silently. Symbols are mapped, never a key.

## Versioning

Bump the version and commit a new file; never edit one in place. Results
reference the index version through the sha they were written under, and an
in-place edit makes every earlier result unreproducible without saying so.

## One rule that is not negotiable

**Target genes must not be in the index used to build the reference matrix.**
`build_signature()` asserts this and will refuse to run. See CLAUDE.md
invariant 2.

**The committed index does contain them, and must.** W3's bulk matrix carries
GUCA2A and CDX2 as *outcome* variables — the week-2 premise check is about their
distribution. So the rule binds on the argument passed to `build_signature()`,
not on the file:

```python
from src.bulk.gene_index import load_gene_index, load_gene_index_map, target_free_index

ids = load_gene_index("0.9.0")
build_signature(..., gene_index=target_free_index(ids, load_gene_index_map("0.9.0"), targets))
```

Note that `signature.py:96` currently asserts against the *whole* index, so
passing `ids` directly raises. That mismatch is
[open decision #12](../../docs/open_decisions.md) — the wrapper above is the
holding position.
