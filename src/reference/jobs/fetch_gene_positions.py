"""Build inferCNV's gene-order file, matched to THIS deposit's build. W1, week 2.

inferCNV needs genomic coordinates for every gene, and getting the build wrong
does not fail — it produces chromosome-arm artifacts that look exactly like real
copy-number events. So the build is not a detail to settle later.

**Which file, and why not the usual one.** GSE178341's feature IDs carry
``GRCh37_liftover_v28``: GENCODE **v28** gene models lifted onto **GRCh37
(hg19)** coordinates. inferCNV's own documentation points people at a
``gencode_v19_gene_pos.txt``, and v19 is the wrong answer here — it is nine
releases earlier, so every gene added between v19 and v28 would be dropped from
the CNV inference silently. GENCODE publishes the exact match as a "GRCh37
mapping" release, and that is what this fetches.

    python src/reference/jobs/fetch_gene_positions.py
    python src/reference/jobs/fetch_gene_positions.py --gtf /path/to/local.gtf.gz

Emits ``$BRP_DATA_DIR/raw/gene_order_hg19.txt`` and prints the
``data/manifest.csv`` row, which goes in the same PR as the code that reads it
(CONTRIBUTING §4).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import os
import re
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

#: GENCODE v28 mapped onto GRCh37 — the deposit's own build string, not v19.
GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_28/GRCh37_mapping/gencode.v28lift37.annotation.gtf.gz"
)

#: Autosomes plus the sex chromosomes, in the order inferCNV expects them.
#: Scaffolds and chrM are excluded deliberately: inferCNV infers CNV by
#: smoothing along a genomic axis, and a contig with a handful of genes on it
#: contributes noise with no positional meaning. chrM has no copy number in the
#: relevant sense at all.
CHROMOSOMES = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])

_GENE_NAME = re.compile(r'gene_name "([^"]+)"')


def parse_gtf(handle) -> list[tuple[str, str, int, int]]:
    """Gene rows only: ``(symbol, chromosome, start, end)``.

    Takes the **first** occurrence of each symbol. GENCODE carries a handful of
    symbols on more than one contig — pseudoautosomal genes appear on both X and
    Y as ``_PAR_Y`` entries, and a few names are reused — and a gene at two
    positions makes the smoothing window ambiguous. Keeping the first is
    reproducible; averaging the positions would be inventing a locus.
    """
    seen: set[str] = set()
    rows: list[tuple[str, str, int, int]] = []
    for line in handle:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        chromosome = fields[0]
        if chromosome not in CHROMOSOMES:
            continue
        match = _GENE_NAME.search(fields[8])
        if not match:
            continue
        symbol = match.group(1)
        if symbol.endswith("_PAR_Y") or symbol in seen:
            continue
        seen.add(symbol)
        rows.append((symbol, chromosome, int(fields[3]), int(fields[4])))
    return rows


def sort_rows(rows):
    """Genomic order. inferCNV walks the file top to bottom as the genome."""
    order = {name: i for i, name in enumerate(CHROMOSOMES)}
    return sorted(rows, key=lambda r: (order[r[1]], r[2], r[3]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", help="local GTF instead of downloading")
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw"
    data.mkdir(parents=True, exist_ok=True)
    out = Path(args.out) if args.out else data / "gene_order_hg19.txt"

    if args.gtf:
        source = Path(args.gtf)
        print(f"reading {source}")
        opener = gzip.open if source.suffix == ".gz" else open
        with opener(source, "rt") as handle:
            rows = parse_gtf(handle)
        source_url = str(source)
    else:
        import urllib.request

        print(f"downloading {GTF_URL}")
        print("  (~40 MB; this is the GRCh37 mapping of v28, matching the "
              "deposit's\n   GRCh37_liftover_v28 tag — NOT the v19 file "
              "inferCNV's docs suggest)")
        with urllib.request.urlopen(GTF_URL) as response:  # noqa: S310
            with gzip.open(response, "rt") as handle:
                rows = parse_gtf(handle)
        source_url = GTF_URL

    rows = sort_rows(rows)
    if not rows:
        raise SystemExit("no gene rows parsed — is this a GENCODE GTF?")

    with out.open("w") as handle:
        for symbol, chromosome, start, end in rows:
            handle.write(f"{symbol}\t{chromosome}\t{start}\t{end}\n")

    digest = hashlib.sha256(out.read_bytes()).hexdigest()
    size = out.stat().st_size
    print(f"\nwrote {out}")
    print(f"  {len(rows):,} genes across {len({r[1] for r in rows})} chromosomes")
    print("\nAdd this row to data/manifest.csv, in the same PR as the code that "
          "reads it:\n")
    # Repo-relative, not the absolute cluster path. The manifest travels with
    # the repo and every machine points BRP_DATA_DIR somewhere different, so an
    # absolute path here names a file nobody else has.
    try:
        relative = out.relative_to(Path(os.environ.get("BRP_DATA_DIR", "data")))
        manifest_path = Path("data") / relative
    except ValueError:
        manifest_path = out
    print(
        f"{manifest_path},{digest},{size},{source_url},GENCODE-v28lift37,"
        f"{date.today().isoformat()},W1,"
        f"inferCNV gene order. GENCODE v28 mapped to GRCh37 — matches "
        f"GSE178341's GRCh37_liftover_v28 tag. NOT v19: nine releases earlier "
        f"would drop genes silently.,"
    )
    print(
        "\nThen point the array job at it:\n"
        f"  export BRP_GENE_POSITIONS={out}\n"
        "  qsub -t 1-5 src/reference/jobs/infercnv.sh patients_pilot.txt"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
