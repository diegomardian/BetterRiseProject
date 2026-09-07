"""Are this project's genes on an in-situ probe set at all? — avenue C1's gate.

    python -m src.reference.jobs.panel_coverage

`docs/NEXT_AVENUES.md` ranked item 7: *"Whether the six panel genes are on a
Xenium/CosMx probe set is a lookup, not a pipeline. If they are not, C1 as
conceived is dead and the long-term plan changes today."* This is that lookup,
written down so it is checkable rather than remembered.

WHY IT IS A JOB AND NOT A NOTE. The answer decides whether a whole avenue
exists, and a vendor changes a panel between revisions. So the gene lists are
recorded in ``data/manifest.csv`` with sha256 and source URL, the match runs
from those files, and the result is a versioned table. A future reader can
re-run it against a newer revision and see what moved.

THE MATCH IS EXACT, ON THE SYMBOL COLUMN ONLY, AND THAT IS DELIBERATE. The
first pass at this lookup searched the whole vendor sheet against a HAND-WRITTEN
alias list and reported **MS4A12 PRESENT** in the 6K panel. It is not. The alias
list wrongly gave ``CD20L4`` to MS4A12; ``CD20L4`` belongs to **MS4A7**, which
is on the panel, so a wrong alias plus a whole-sheet search returned a gene that
is not there. **Presence is the answer that keeps an avenue alive**, so it is
the direction a sloppy match must not be allowed to fail in.

The fix is not "search fewer columns" — :func:`alias_audit` shows a whole-sheet
search for the exact symbol disagrees with the strict match on **nothing** here.
The fix is that no alias mapping is used at all: the match is exact against
``Gene Symbol(s)``, split on separators, and `alias_audit` is carried so that a
future reader can see the strict answer is not an artefact of looking in one
column. A gene genuinely listed only under an alias would show up there.

THE 10X HALF IS NOW CHECKED. Both 10x panel metadata CSVs were downloaded on
2026-09-07, pinned in ``data/manifest.csv``, and read by exact symbol column.
They are deliberately treated like the CosMx files rather than inferred from a
marketing-page gene count.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import re
import sys
from pathlib import Path

import pandas as pd

from src.common.io import write_versioned_table
from src.common.paths import DATA_DIR, RESULTS_DIR
from src.reference.jobs.coexpression_silencing import GENE_ROLES

log = logging.getLogger(__name__)

#: Invariant 8's Wnt target signature, carried because D1 uses it and a spatial
#: follow-up would need it on the same slide.
WNT_SIGNATURE = ("AXIN2", "NKD1", "RNF43", "NOTUM", "TCF7")

#: Panels whose gene list is in hand. ``sha256`` pins the revision.
PANELS = {
    "CosMx_Human_Universal_1K": {
        "path": "raw/panels/CosMx_Human_Universal_1K_gene_target_list.xlsx",
        "sheet": "Gene and Probe Details",
        "symbol_column": "Gene Symbol(s)",
        "plex": 1000,
        "vendor_doc": "LBL-11178-04",
        "sha256": "941832fddbc58c52821a7202b780e19b14cd418d5554c4cc1ae780d54ace9cc7",
    },
    "CosMx_Human_6K_Discovery": {
        "path": "raw/panels/CosMx_Human_6K_Discovery_gene_list.xlsx",
        "sheet": "Gene and Probe Details",
        "symbol_column": "Gene Symbol(s)",
        "plex": 6000,
        "vendor_doc": "LBL-11190-04",
        "sha256": "14395559815719c62eee8593149e52e75a0b038ce9728fa37f5818c747d3bcf9",
    },
    "Xenium_Human_Colon_v1": {
        "path": "raw/panels/Xenium_Human_Colon_v1_panel_metadata.csv",
        "format": "csv",
        "symbol_column": "Gene",
        "plex": 322,
        "vendor_doc": "Xenium Human Colon v1 panel metadata",
        "sha256": "68950fb07b11205fa3602762c3011b8b1a14914b044328dec6939b8c1be0da39",
    },
    "Xenium_Prime_5K_Human_Pan_Tissue": {
        "path": "raw/panels/Xenium_Prime_5K_Human_Pan_Tissue_panel_metadata.csv",
        "format": "csv",
        "symbol_column": "gene_name",
        "plex": 5001,
        "vendor_doc": "Xenium Prime 5K Human Pan Tissue panel metadata",
        "sha256": "833ddb3008f2eb5e1b39053a988869c74db82cbaa8124d932965459f3169c58a",
    },
}

_SPLIT = re.compile(r"[,;|/]")


class PanelError(RuntimeError):
    """A panel file is missing or is not the revision that was checked."""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_panel_frame(spec: dict, data_dir: Path) -> pd.DataFrame:
    """Read one verified vendor panel, preserving its declared symbol column."""
    path = data_dir / spec["path"]
    if not path.exists():
        raise PanelError(
            f"{path} absent. It is gitignored and recorded in "
            f"data/manifest.csv — re-fetch from the source_url there."
        )
    digest = sha256(path)
    if digest != spec["sha256"]:
        raise PanelError(
            f"{path} is sha256 {digest}, not the {spec['sha256']} this result "
            f"was computed against. The vendor revised the panel; re-run and "
            f"record the new answer rather than trusting the old one."
        )
    if spec.get("format", "xlsx") == "csv":
        frame = pd.read_csv(path)
    else:
        frame = pd.read_excel(path, sheet_name=spec["sheet"], header=1)
    column = spec["symbol_column"]
    if column not in frame.columns:
        raise PanelError(f"{path} has no {column!r}; got {list(frame.columns)}")
    return frame


def read_symbols(spec: dict, data_dir: Path) -> set[str]:
    """Exact gene symbols from one vendor sheet, with the revision verified."""
    frame = read_panel_frame(spec, data_dir)
    column = spec["symbol_column"]
    symbols: set[str] = set()
    for value in frame[column].dropna():
        symbols |= {t.strip().upper() for t in _SPLIT.split(str(value)) if t.strip()}
    return symbols


def alias_audit(spec: dict, data_dir: Path, genes: tuple[str, ...]) -> set[str]:
    """What a search of EVERY column for the exact symbol would have claimed.

    Reported, not used. On the panels in hand it agrees with the strict match
    everywhere, which is the point of carrying it: the strict answer is not an
    artefact of reading one column. It does NOT reproduce this lookup's original
    error, which came from a hand-written alias table rather than from the
    sheet — no alias mapping is used anywhere in this module, deliberately.
    """
    frame = read_panel_frame(spec, data_dir)
    blob = {str(v).strip().upper() for c in frame.columns for v in frame[c].dropna()}
    exploded: set[str] = set()
    for value in blob:
        exploded |= {t.strip() for t in _SPLIT.split(value) if t.strip()}
    return {g for g in genes if g.upper() in exploded}


def coverage(data_dir: Path) -> pd.DataFrame:
    genes = tuple(GENE_ROLES) + WNT_SIGNATURE
    roles = {**GENE_ROLES, **{g: "wnt_target" for g in WNT_SIGNATURE}}
    rows = []
    for panel, spec in PANELS.items():
        symbols = read_symbols(spec, data_dir)
        permissive = alias_audit(spec, data_dir, genes)
        for gene in genes:
            present = gene.upper() in symbols
            rows.append({
                "panel": panel,
                "plex": spec["plex"],
                "vendor_doc": spec["vendor_doc"],
                "gene": gene,
                "role": roles[gene],
                "status": "present" if present else "absent",
                "present": present,
                "permissive_match_would_say": gene in permissive,
                "match_disagrees": bool(gene in permissive) != present,
                "n_symbols_in_panel": len(symbols),
            })
    return pd.DataFrame(rows)


def verdict(table: pd.DataFrame) -> dict[str, str]:
    """C1's gate: can the design be run on a stock panel?

    The design needs a target, at least one identity marker to make the claim
    gene-specific rather than tier-level, and a control. Missing any one of
    those roles ends the stock-panel design; a target alone is not a runnable
    discriminator.
    """
    if table.empty:
        return {"verdict": "NOT CHECKED",
                "detail": "no panel list in hand; nothing was determined."}
    lines = []
    runnable = []
    for panel, group in table.groupby("panel", sort=True):
        got = group[group["status"] == "present"]["gene"].tolist()
        missing_roles = [
            role for role in ("target", "identity", "control")
            if not (group[(group["role"] == role) & (group["status"] == "present")].shape[0])
        ]
        if not missing_roles:
            runnable.append(panel)
        lines.append(
            f"{panel}: {len(got)}/{len(group)} present"
            + (f", MISSING REQUIRED ROLE(S) ({', '.join(missing_roles)})"
               if missing_roles else ", RUNNABLE")
        )
    if not runnable:
        return {
            "verdict": "C1 NOT RUNNABLE ON A STOCK PANEL",
            "detail": (
                "no checked stock panel contains a target, an identity marker, "
                "and a control. C1 needs a custom panel or a different spatial "
                "design; a target-only panel does not license the claim. "
                + "; ".join(lines)
            ),
        }
    return {"verdict": "C1 RUNNABLE ON STOCK PANEL(S)",
            "detail": f"runnable: {', '.join(runnable)}; " + "; ".join(lines)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    parser.add_argument("--seed", type=int, default=20260907)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    table = coverage(args.data_dir)
    log.info("\n%s\nPANEL COVERAGE — exact match on the vendor symbol column\n%s",
             "=" * 72, "=" * 72)
    wide = table.pivot(index=["gene", "role"], columns="panel", values="status")
    log.info("%s", wide.to_string())

    disagree = table[table["match_disagrees"]]
    log.info("\nAlias-permissive match would have disagreed on %d gene/panel "
             "pair(s): %s", len(disagree),
             ", ".join(f"{r.gene}@{r.panel}" for r in disagree.itertuples()) or "none")

    outcome = verdict(table)
    log.info("\n%s\nVERDICT\n%s", "=" * 72, "=" * 72)
    log.info("  %s\n  %s", outcome["verdict"], outcome["detail"])

    meta = {
        "what_this_is": "avenue C1's feasibility gate — a lookup, recorded",
        "ranked_item": "docs/NEXT_AVENUES.md 'Still open, ranked' item 7",
        "match_rule": (
            "exact, on the vendor's symbol column only; no alias mapping is "
            "used. A hand-written alias table in the first pass wrongly gave "
            "CD20L4 to MS4A12 (it is MS4A7's) and manufactured a PRESENT. "
            "permissive_match_would_say records what a whole-sheet search for "
            "the exact symbol claims — here, the same thing."
        ),
        "panels_checked": {k: v["vendor_doc"] for k, v in PANELS.items()},
        "panel_sha256": {k: v["sha256"] for k, v in PANELS.items()},
        "verdict": outcome,
        "exploratory": False,
        "pre_registered": False,
    }
    log.info("\nwrote %s", write_versioned_table(
        table, "panel_coverage", seed=args.seed, results_dir=args.results_dir,
        allow_dirty=args.allow_dirty, extra_meta=meta,
    ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
