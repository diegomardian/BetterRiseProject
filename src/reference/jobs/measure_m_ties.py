#!/usr/bin/env python
"""How tied is M? A question raised in review, answered before G1 runs. W1.

**Why this exists.** Reviewing W2's ratification of prereg amendment 2
([PR #40](https://github.com/diegomardian/BetterRiseProject/pull/40)) turned up
a defect in `within_bin_percentile`: its docstring promises average ranks and it
computes ordinal ones, so tied genes get percentiles decided by their **position
in the array** rather than by their value.

That would be cosmetic almost anywhere else. It is not here, because it lands on
**threshold 2 — MS4A12 >= 0.50** — and MS4A12 is the *retained* control, the gene
amendment 2 expects to sit at M ~ 0. It is therefore the gene most likely to fall
inside a tie mass at exactly zero. Simulated at a 20% tied fraction, MS4A12's
percentile moved 0.417 -> 0.564 purely by array position, which straddles the
threshold.

**What was NOT shown in that review** is whether real data is materially tied.
Only the mechanism and its magnitude were demonstrated. This job settles it, and
it is W1's data to settle it with.

The honest outcomes, stated before running:

- **Ties under ~1%** — the finding is theoretical. Say so on #40, keep the rank
  fix anyway because it costs nothing, and stop worrying about threshold 2 on
  this ground.
- **Ties above ~10%** — the fix is load-bearing, and threshold 2 carries two
  independent sources of arbitrariness: this one and the coin-flip power W2
  already measured on `isolated_tier_a_loss`.

**This is not a G1 run.** It computes M and A because they are what ties live in,
and deliberately does not evaluate any G1 threshold. `checks.py` still returns
``not_estimable``, and it stays that way until the team closes #37.

    python src/reference/jobs/measure_m_ties.py
    python src/reference/jobs/measure_m_ties.py --patients C122 C165
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from src.common.io import write_versioned_table  # noqa: E402
from src.common.panel import tier_genes  # noqa: E402
from src.common.provenance import DEFAULT_SEED, set_global_seeds  # noqa: E402
from src.reference.gene_index import read_gene_index  # noqa: E402
from src.reference.ingest import (  # noqa: E402
    assign_compartments,
    read_gse178341,
    read_gse178341_clusters,
    read_gse178341_index,
)

#: Amendment 2 §4: 20 equal-size bins by A.
N_ABUNDANCE_BINS = 20

#: Below this many epithelial cells in an arm, a per-gene mean is too noisy to
#: read and the patient contributes nothing. Same floor the full run uses.
MIN_CELLS_PER_ARM = 50

GENE_INDEX_VERSION = "1.0.0"


def _arm_mean(matrix, mask: np.ndarray) -> np.ndarray:
    """Per-gene mean over the selected cells, without densifying the matrix."""
    return np.asarray(matrix[mask].mean(axis=0)).ravel()


def within_bin_percentile(a: np.ndarray, m: np.ndarray, *, average: bool) -> np.ndarray:
    """Percentile of M within equal-size abundance bins.

    ``average=True`` is what amendment 2's rule means and what W2's docstring
    says; ``average=False`` reproduces the ordinal behaviour actually shipped in
    `src/harness/g1_amendment.py`. Both are computed so the *difference* is
    measured on real data rather than argued about.
    """
    order = np.argsort(a, kind="stable")
    bin_of = np.empty(len(a), dtype=int)
    bin_of[order] = np.arange(len(a)) * N_ABUNDANCE_BINS // len(a)

    out = np.empty(len(a), dtype=float)
    for b in range(N_ABUNDANCE_BINS):
        selected = bin_of == b
        values = m[selected]
        if average:
            from scipy.stats import rankdata

            ranks = rankdata(values, method="average") - 1.0
        else:
            ranks = values.argsort().argsort().astype(float)
        out[selected] = (ranks + 0.5) / len(values)
    return out


def measure_patient(
    matrix, ensembl, symbols, tissue, compartment, index_set
) -> dict | None:
    """M, A and the tie structure for one patient. None if either arm is thin."""
    epithelial = compartment == "epithelial"
    on_index = np.array([g in index_set for g in ensembl], dtype=bool)
    tumour = epithelial & (tissue == "tumour")
    normal = epithelial & (tissue == "normal")
    if tumour.sum() < MIN_CELLS_PER_ARM or normal.sum() < MIN_CELLS_PER_ARM:
        return None

    mean_t = _arm_mean(matrix, tumour)
    mean_n = _arm_mean(matrix, normal)

    # Amendment 2 §4: restricted to genes with a non-zero mean in BOTH arms. A
    # log ratio against zero is not a number, and a pseudocount there invents
    # the quantity being measured.
    usable = on_index & (mean_t > 0) & (mean_n > 0)
    if usable.sum() < N_ABUNDANCE_BINS:
        return None

    t, n = mean_t[usable], mean_n[usable]
    m = np.log2(t / n)
    a = 0.5 * (np.log2(t) + np.log2(n))
    # SYMBOLS, not Ensembl ids. The panel is written as symbols and this
    # matrix is keyed on Ensembl, so matching tier_genes() against the
    # Ensembl array would find nothing and report MS4A12 as absent on every
    # patient — the identifier-space defect from issue #35, one more time.
    gene_symbols = np.asarray(symbols)[usable]

    # The tie structure. Exact equality is the thing that matters: ordinal ranks
    # only diverge from average ranks on exact ties.
    _, counts = np.unique(m, return_counts=True)
    tied = int(counts[counts > 1].sum())
    exact_zero = int((m == 0.0).sum())

    ordinal = within_bin_percentile(a, m, average=False)
    averaged = within_bin_percentile(a, m, average=True)
    gap = np.abs(ordinal - averaged)

    row = {
        "n_genes": int(usable.sum()),
        "n_cells_tumour": int(tumour.sum()),
        "n_cells_normal": int(normal.sum()),
        "frac_tied": tied / int(usable.sum()),
        "frac_exactly_zero": exact_zero / int(usable.sum()),
        "largest_tie_cluster": int(counts.max()),
        "max_percentile_gap": float(gap.max()),
        "median_percentile_gap": float(np.median(gap)),
    }

    # Where the retained control actually sits, under both rules. This is the
    # number threshold 2 reads, and the whole reason the review flagged it.
    for symbol, ids in (("MS4A12", tier_genes("D")),):
        hit = np.isin(gene_symbols, ids)
        if hit.any():
            j = int(np.flatnonzero(hit)[0])
            row[f"{symbol}_pct_ordinal"] = float(ordinal[j])
            row[f"{symbol}_pct_average"] = float(averaged[j])
            row[f"{symbol}_M"] = float(m[j])
        else:
            row[f"{symbol}_pct_ordinal"] = None
            row[f"{symbol}_pct_average"] = None
            row[f"{symbol}_M"] = None
    return row


def _report(out: pd.DataFrame) -> None:
    print("\n" + "=" * 64)
    print("HOW TIED IS M?")
    print("=" * 64)
    print(out[["frac_tied", "frac_exactly_zero", "largest_tie_cluster",
               "max_percentile_gap"]].describe(
                   percentiles=[0.5, 0.9]).to_string())

    median_tied = float(out["frac_tied"].median())
    print(f"\nmedian tied fraction across patients: {median_tied:.2%}")
    if median_tied < 0.01:
        print(
            "\n  UNDER 1% — the review finding is theoretical on this cohort.\n"
            "  Say so on #40. Keep the rankdata fix anyway: it costs nothing and\n"
            "  the docstring already claims it."
        )
    elif median_tied > 0.10:
        print(
            "\n  ABOVE 10% — the rank fix is LOAD-BEARING. Threshold 2 then\n"
            "  carries two independent sources of arbitrariness: tie order, and\n"
            "  the coin-flip power W2 measured on isolated_tier_a_loss. Both\n"
            "  belong in the amendment beside the threshold."
        )
    else:
        print(
            "\n  BETWEEN 1% AND 10% — real but not dominant. Report the number\n"
            "  rather than a verdict, and let #40 decide with it in hand."
        )

    ms = out.dropna(subset=["MS4A12_pct_ordinal"])
    print("\n" + "=" * 64)
    print("WHERE MS4A12 SITS — the number threshold 2 reads")
    print("=" * 64)
    if len(ms):
        moved = (
            (ms["MS4A12_pct_ordinal"] >= 0.50) != (ms["MS4A12_pct_average"] >= 0.50)
        )
        print(f"  patients with MS4A12 on the index: {len(ms)}")
        print(f"  median percentile, ordinal: {ms['MS4A12_pct_ordinal'].median():.3f}")
        print(f"  median percentile, average: {ms['MS4A12_pct_average'].median():.3f}")
        print(f"  patients where the two rules land on OPPOSITE sides of 0.50: "
              f"{int(moved.sum())} of {len(ms)}")
        if moved.any():
            print(
                "\n  !! For those patients the rank convention alone decides\n"
                "     threshold 2. That is the review's claim, measured."
            )
    else:
        print("  MS4A12 is not on the index in any patient — check the join")

    print("\nNOTE: this is not a G1 result. No threshold was evaluated; "
          "checks.py\n      still returns not_estimable until #37 closes.")


def main() -> int:
    set_global_seeds(DEFAULT_SEED)
    parser = argparse.ArgumentParser()
    parser.add_argument("--patients", nargs="*", default=None)
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()

    data = Path(os.environ.get("BRP_DATA_DIR", "data")) / "raw" / "GSE178341"
    h5 = data / "GSE178341_crc10x_full_c295v4_submit.h5"
    clusters = read_gse178341_clusters(
        data / "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"
    )
    index_set = set(read_gene_index(GENE_INDEX_VERSION))
    print(f"gene index {GENE_INDEX_VERSION}: {len(index_set):,} genes")

    obs, _var = read_gse178341_index(h5)
    patients = args.patients or sorted(obs["patient_id"].unique())
    print(f"{len(patients)} patients")

    rows = []
    for i, patient in enumerate(patients, 1):
        adata = read_gse178341(h5, patients=[patient])
        compartment = assign_compartments(clusters).reindex(adata.obs.index).to_numpy()
        row = measure_patient(
            adata.X,
            adata.var["ensembl_id"].astype(str).to_numpy(),
            adata.var["gene_symbol"].astype(str).to_numpy(),
            adata.obs["tissue"].astype(str).to_numpy(),
            compartment,
            index_set,
        )
        if row is None:
            print(f"[{i}/{len(patients)}] {patient} — not paired or too thin, skipped")
            continue
        row["patient_id"] = patient
        rows.append(row)
        ms4a12 = row["MS4A12_pct_average"]
        where = "n/a" if ms4a12 is None else f"{ms4a12:.3f}"
        print(f"[{i}/{len(patients)}] {patient} — {row['n_genes']:,} genes, "
              f"{row['frac_tied']:.2%} tied, {row['frac_exactly_zero']:.2%} at "
              f"M=0, MS4A12 pct {where}")

    if not rows:
        raise SystemExit("no patient produced a paired measurement")
    out = pd.DataFrame(rows)
    _report(out)

    path = write_versioned_table(
        out, "m_tie_structure", seed=DEFAULT_SEED, allow_dirty=args.allow_dirty,
        notes=(
            "Tie structure of M = log2(tumour/normal) per patient, on gene index "
            "1.0.0, epithelial cells only, genes non-zero in both arms. Answers a "
            "question raised reviewing PR #40: W2's within_bin_percentile computes "
            "ORDINAL ranks while its docstring promises average ranks, which makes "
            "a tied gene's percentile depend on array position. MS4A12 is the gene "
            "most exposed, because threshold 2 reads it and it is expected at M~0. "
            "NOT a G1 result — no threshold is evaluated here."
        ),
        extra_meta={
            "n_patients": int(len(out)),
            "median_frac_tied": float(out["frac_tied"].median()),
            "gene_index_version": GENE_INDEX_VERSION,
        },
    )
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
