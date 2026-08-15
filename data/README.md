# data/

**Nothing in here is committed.** `data/` is gitignored in full; only this
README and [`manifest.csv`](manifest.csv) travel with the repo.

```
data/raw/         exactly as downloaded. Never edited, never written to. Read-only by habit.
data/interim/     intermediate objects. Deletable — anything here must be reproducible from raw/ by code in src/.
data/processed/   the artifacts other workstreams consume. S matrices, harmonised bulk, clinical tables.
```

If you want the data somewhere else — a scratch disk, a shared drive, not inside
your repo folder — set `BRP_DATA_DIR` and every path in
[`src/common/paths.py`](../src/common/paths.py) follows:

```bash
export BRP_DATA_DIR=/scratch/brp/data     # bash
$env:BRP_DATA_DIR = "D:\brp\data"          # PowerShell
```

## The manifest is the contract

Four people, four machines. The only way anyone knows they have the same bytes
you have is [`manifest.csv`](manifest.csv). **Add the row in the same PR as the
code that reads the file.** A download nobody recorded is a result nobody can
reproduce.

```
path,sha256,bytes,source_url,accession,downloaded_on,downloaded_by,workstream,notes
```

Get a checksum:

```bash
sha256sum data/raw/GSE178341/counts.mtx.gz          # linux/mac
Get-FileHash data\raw\GSE178341\counts.mtx.gz -Algorithm SHA256   # PowerShell
```

## Where things go

| Accession | What | Owner | Lands in |
|---|---|---|---|
| GSE178341 | Pelka 2021 — 371k cells, 62 patients, matched normal. **Primary.** | W1 | `raw/GSE178341/` |
| ICBI CRC atlas | 4.27M cells, 650 patients, 48 studies. **Metadata table first, object later.** | W1 | `raw/icbi/` |
| HTAN / Vanderbilt polyp atlas | Conventional vs. serrated, crypt-top colonocytes pre-annotated | W1 | `raw/htan/` |
| GSE132465 / GSE144735 | Lee SMC + KUL3, matched normal. Replication. | W4 | `raw/lee/` |
| Joanito 2022 | iCMS subtypes with Wnt/MYC annotation | W4 | `raw/joanito/` |
| TCGA-COAD/READ | STAR counts, GDC | W3 | `raw/tcga/` |
| TCGA-CDR | Survival — DSS and PFI primary | W3 | `raw/tcga_cdr/` |
| Becker/Chang multiome | Chromatin axis, week 13+ | W4 | `raw/multiome/` |
| Visium HD / Xenium | Spatial validation, week 13+ | W1 | `raw/spatial/` |

## Two things worth checking before you download 371k cells

- **Verify you have raw counts, not normalised values.** GSE178341's
  supplementary structure is awkward; budget a day just for parsing.
- **Pull the ICBI metadata table only, first.** That one table gives the real
  sample size — paired tumour/normal counts, epithelial fraction by study,
  platform mix — before you commit any compute.
