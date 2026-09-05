"""Stage 4 end to end on synthetic data: S matrix -> fractions -> gate -> verdict.

The whole chain, on the laptop, with a truth we set. Its job is to prove the
pipeline runs and that each gate stops it where it should, BEFORE the cluster
run where none of that can be checked cheaply.

The synthetic world is built so the pre-registered prediction has a known right
answer: CDX2 is made compositional (its expression is a function of the mature
fraction) and GUCA2A is not (it is noise at a low abundance). If the arm cannot
recover that on data we constructed, it cannot be read on data we did not.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src.bulk.deconvolution import LOG1P_CP10K, linearise, load_reference
from src.common.paths import RESULTS_DIR

COMMITTED = RESULTS_DIR / "2026-08-26_63ead2e"
N_SAMPLES = 220
SEED = 20260905


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    """A cohort with known fractions, purity, and a known compositional gene."""
    s_path = COMMITTED / "S_matrix_lineage_1.0.0.parquet"
    if not s_path.exists():
        pytest.skip("S matrix not committed")
    from src.bulk.gene_index import load_gene_index_map, resolve_symbols

    rng = np.random.default_rng(SEED)
    tmp = tmp_path_factory.mktemp("stage4")
    reference = linearise(load_reference(s_path, rung="lineage", scale=LOG1P_CP10K))
    types = list(reference.matrix.columns)

    samples = [f"TCGA-XX-{i:04d}-01A-11R-AAAA-07" for i in range(N_SAMPLES)]
    truth = pd.DataFrame(
        rng.dirichlet([3.0 if t == "differentiated" else 1.5 for t in types], N_SAMPLES),
        columns=types, index=samples,
    )

    # Bulk TPM: a linear mixture, renormalised.
    values = truth.to_numpy() @ reference.matrix.to_numpy().T
    tpm = pd.DataFrame(values, index=samples, columns=reference.matrix.index)
    tpm = tpm.div(tpm.sum(axis=1), axis=0) * 1e6
    tpm_path = tmp / "tcga_tpm_1.0.0.parquet"
    tpm.to_parquet(tpm_path)

    # ABSOLUTE purity: (1 - purity) tracks the non-epithelial fraction, so the
    # pre-committed gate has something real to detect.
    non_epi = truth[["immune", "stromal", "endothelial"]].sum(axis=1)
    purity = pd.DataFrame({
        "barcode": samples,
        "method": "absolute",
        "purity": np.clip(1.0 - (non_epi + rng.normal(0, 0.04, N_SAMPLES)), 0.05, 0.95),
        "expression_derived": False,
    })
    purity_path = tmp / "tcga_purity.parquet"
    purity.to_parquet(purity_path)

    # log2-CPM outcome matrix on the real gene index, so matched_null_genes has
    # a real `on_panel` column and a real abundance distribution to draw from.
    index_map = load_gene_index_map("1.0.0")
    ids, _, _ = resolve_symbols(index_map, ["GUCA2A", "CDX2", "ACTB", "GAPDH"])
    pool = [g for g in index_map["ensembl_id"] if g not in ids.values()][:900]
    columns = [*ids.values(), *pool]

    mature = truth["differentiated"].to_numpy()
    z = (mature - mature.mean()) / mature.std()
    expression = pd.DataFrame(index=samples, columns=columns, dtype=float)
    # Every gene sits near log2CPM 3, so the matched null has candidates within
    # the spec's +/-0.25 window; the target genes are placed in that band too.
    expression[pool] = 3.0 + rng.normal(0, 1.0, (N_SAMPLES, len(pool)))
    expression[ids["CDX2"]] = 3.0 + 1.4 * z + rng.normal(0, 0.6, N_SAMPLES)   # compositional
    expression[ids["GUCA2A"]] = 3.0 + rng.normal(0, 1.0, N_SAMPLES)           # not
    expression[ids["ACTB"]] = 3.0 + rng.normal(0, 1.0, N_SAMPLES)             # control
    expression[ids["GAPDH"]] = 3.0 + rng.normal(0, 1.0, N_SAMPLES)
    expression_path = tmp / "tcga_log2cpm_1.0.0.parquet"
    expression.to_parquet(expression_path)

    return {
        "tmp": tmp, "tpm": tpm_path, "purity": purity_path,
        "expression": expression_path, "truth": truth, "ids": ids,
    }


@pytest.fixture(scope="module")
def fractions(world):
    """Step 1: deconvolution. Returns the results directory it wrote into."""
    from src.bulk.run_deconvolution import main

    results = world["tmp"] / "results"
    code = main([
        "--rung", "lineage", "--bulk", str(world["tpm"]),
        "--s-matrix-dir", str(COMMITTED), "--linearise-reference",
        "--results-dir", str(results), "--allow-dirty",
    ])
    assert code == 0
    return results


def test_step1_recovers_the_mature_fraction_it_was_given(world, fractions):
    table = pd.read_parquet(sorted(fractions.glob("*/stage4_fractions.parquet"))[0])
    truth = world["truth"]["differentiated"]
    for method in table["method"].unique():
        got = table[table["method"] == method].set_index("sample_id")
        r = np.corrcoef(got.loc[truth.index, "mature_colonocyte_fraction"], truth)[0, 1]
        assert r > 0.5, f"{method} recovers the mature fraction at only r={r:.3f}"


def test_the_whole_chain_runs_and_reaches_the_pre_registered_verdict(world, fractions):
    """Step 2: gate, then arm. The verdict must be the one we built in."""
    from src.bulk.run_stage4_variance import main

    out = world["tmp"] / "verdicts"
    code = main([
        "--fractions", str(sorted(fractions.glob("*/stage4_fractions.parquet"))[0]),
        "--predictor-checks", str(sorted(fractions.glob("*/stage4_predictor_checks.parquet"))[0]),
        "--expression", str(world["expression"]),
        "--purity", str(world["purity"]),
        "--rung", "lineage", "--results-dir", str(out), "--allow-dirty",
    ])
    assert code == 0, f"the chain stopped with code {code}"

    gate = pd.read_parquet(sorted(out.glob("*/stage4_instrument_gate.parquet"))[0])
    assert gate["passed"].any()
    assert (gate["purity_method"] == "absolute").all()

    verdicts = pd.read_parquet(sorted(out.glob("*/stage4_variance_verdicts.parquet"))[0])
    assert set(verdicts["gene"]) == {"GUCA2A", "CDX2"}
    assert set(verdicts["r_squared_kind"]) == {"marginal", "partial"}

    marginal = verdicts[verdicts["r_squared_kind"] == "marginal"].set_index(["method", "gene"])
    for method in verdicts["method"].unique():
        assert marginal.loc[(method, "CDX2"), "exceeds_null"], (
            "CDX2 was built compositional and must exceed its matched null"
        )
        assert not marginal.loc[(method, "GUCA2A"), "exceeds_null"], (
            "GUCA2A was built as noise and must not exceed its matched null"
        )
        assert marginal.loc[(method, "CDX2"), "primary_verdict"] == "confirmed"

    nulls = pd.read_parquet(sorted(out.glob("*/stage4_matched_nulls.parquet"))[0])
    assert nulls["null_gene"].nunique() >= 20
    assert not nulls["target_gene"].isin(["ACTB", "GAPDH"]).any()


def test_the_gate_stops_the_chain_before_any_r_squared_exists(world, fractions):
    """Purity that carries no signal must STOP the run, not merely warn."""
    from src.bulk.run_stage4_variance import main

    rng = np.random.default_rng(3)
    broken = pd.read_parquet(world["purity"])
    broken["purity"] = rng.uniform(0.3, 0.9, len(broken))   # unrelated to composition
    broken_path = world["tmp"] / "purity_broken.parquet"
    broken.to_parquet(broken_path)

    out = world["tmp"] / "stopped"
    code = main([
        "--fractions", str(sorted(fractions.glob("*/stage4_fractions.parquet"))[0]),
        "--expression", str(world["expression"]), "--purity", str(broken_path),
        "--rung", "lineage", "--results-dir", str(out), "--allow-dirty",
    ])
    assert code == 4, "a failed instrument check must stop the analysis"
    assert not list(out.glob("*/stage4_variance_verdicts.parquet")), (
        "an R-squared table was written despite the gate failing. The locked "
        "spec says no R-squared is reported -- not computed then suppressed."
    )
    gate = pd.read_parquet(sorted(out.glob("*/stage4_instrument_gate.parquet"))[0])
    assert not gate["passed"].any()


def test_an_expression_derived_purity_call_cannot_satisfy_the_gate(world, fractions):
    """The 675-barcode trap, at the driver level rather than the unit level."""
    from src.bulk.run_stage4_variance import main

    derived = pd.read_parquet(world["purity"])
    derived["method"] = "estimate_affy_extrapolated"
    derived["expression_derived"] = True
    path = world["tmp"] / "purity_derived.parquet"
    derived.to_parquet(path)

    out = world["tmp"] / "derived"
    code = main([
        "--fractions", str(sorted(fractions.glob("*/stage4_fractions.parquet"))[0]),
        "--expression", str(world["expression"]), "--purity", str(path),
        "--rung", "lineage", "--results-dir", str(out), "--allow-dirty",
    ])
    assert code == 4, (
        "the gate accepted an expression-derived purity call. It would then be "
        "correlating a fraction derived from expression against a purity score "
        "derived from expression, and pass because both read the same signal."
    )


def test_the_run_refuses_a_proposed_prespecification(world, fractions, monkeypatch):
    """The lock is load-bearing: unlock it and nothing runs."""
    import src.bulk.run_stage4_variance as driver
    from src.bulk.prespec import load_prespec

    proposed = dict(load_prespec(), status="proposed")
    monkeypatch.setattr(driver, "load_prespec", lambda *a, **k: proposed)
    with pytest.raises(Exception, match="not 'locked'"):
        driver.main([
            "--fractions", str(sorted(fractions.glob("*/stage4_fractions.parquet"))[0]),
            "--expression", str(world["expression"]), "--purity", str(world["purity"]),
            "--rung", "lineage", "--results-dir", str(world["tmp"] / "x"), "--allow-dirty",
        ])


def _two_rung_fractions(world, fractions) -> Path:
    """The lineage fractions, relabelled under a second rung.

    Each rung's S matrix selects its own 800 genes, so synthetic bulk built for
    one rung cannot be deconvolved against another. This test is about the WRITE
    path -- whether every rung that reaches the gate survives into the committed
    table -- so the second rung is manufactured rather than deconvolved.
    """
    source = sorted(fractions.glob("*/stage4_fractions.parquet"))[0]
    frame = pd.read_parquet(source)
    second = frame.assign(granularity_rung="crypt_position")
    path = world["tmp"] / "two_rung_fractions.parquet"
    pd.concat([frame, second], ignore_index=True).to_parquet(path)
    return path


def test_every_rung_survives_into_one_gate_table(world, fractions):
    """The defect the 2026-09-05 run committed: only the last rung was kept.

    The cluster script looped in the shell, calling the driver once per rung.
    Each call wrote `stage4_instrument_gate` under the same name into the same
    {date}_{sha} directory, so lineage and crypt_position were overwritten by
    best4 and their numbers survived only in a log file. The committed evidence
    for a four-rung result covered one rung.

    `--rung all` loops inside the driver and writes once.
    """
    from src.bulk.run_stage4_variance import main

    out = world["tmp"] / "allrungs"
    main([
        "--fractions", str(_two_rung_fractions(world, fractions)),
        "--expression", str(world["expression"]), "--purity", str(world["purity"]),
        "--rung", "all", "--results-dir", str(out), "--allow-dirty",
    ])
    gate = pd.read_parquet(sorted(out.glob("*/stage4_instrument_gate.parquet"))[0])
    assert gate["granularity_rung"].nunique() > 1, (
        f"the gate table carries only {gate['granularity_rung'].unique()}. Every "
        f"rung that reached the gate must survive into the committed table, not "
        f"just the last one."
    )
    assert {"granularity_rung", "method", "pearson_r", "passed"} <= set(gate.columns)

    verdicts = sorted(out.glob("*/stage4_variance_verdicts.parquet"))
    if verdicts:
        assert pd.read_parquet(verdicts[0])["granularity_rung"].nunique() > 1


def test_a_shell_loop_over_rungs_would_have_lost_all_but_the_last(world, fractions):
    """The demonstration, so the fix is justified rather than merely asserted."""
    from src.bulk.run_stage4_variance import main

    path = _two_rung_fractions(world, fractions)
    out = world["tmp"] / "looped"
    for rung in ("lineage", "crypt_position"):
        main([
            "--fractions", str(path),
            "--expression", str(world["expression"]), "--purity", str(world["purity"]),
            "--rung", rung, "--results-dir", str(out), "--allow-dirty",
        ])
    gate = pd.read_parquet(sorted(out.glob("*/stage4_instrument_gate.parquet"))[0])
    assert set(gate["granularity_rung"]) == {"crypt_position"}, (
        "a shell loop no longer overwrites, so this test's premise is stale"
    )
