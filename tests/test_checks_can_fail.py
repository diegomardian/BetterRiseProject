"""Every guard, against the input that forces it to fail.

This file is the paper's own rule applied to this repository: **a check unable
to fail is worse than no check, because it turns an absence of evidence into a
green light.** The remedy is mutation testing's premise moved from unit tests to
analysis diagnostics — construct the input a guard MUST catch, and treat a guard
with no such test as untested.

Every other test file here asserts the pass path: given clean data, the check
returns clean. That is the half that cannot detect a check which returns clean
on everything. These are the other half.

Each test names the defect it would have caught. Five of them are the five
defects in the paper's Appendix A, written as the inputs that would have caught
them at the time rather than months later.

If you add a guard, add its failing input here. If you cannot construct one, the
guard does not have a testable failure mode and that is the finding.
"""

from __future__ import annotations

import ast
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest

from src.common.io import ReservedMetaKeyError, write_versioned_table
from src.common.paths import REPO_ROOT, RESULTS_DIR
from src.harness.calibration import (
    PREREGISTERED,
    coverage_and_discrimination,
)
from src.harness.depth_confound import (
    MATURITY_DEPTH_RHO_TOLERANCE,
    depth_confound_report,
    max_attainable_rho,
)
from src.reference.signature import LeakageError, assert_no_target_leakage
from src.schema import (
    REQUIRED_COLUMNS,
    SchemaViolation,
    coerce_results,
    validate_results,
)

# ---------------------------------------------------------------------------
# The depth-confound diagnostic
# ---------------------------------------------------------------------------


def _two_arm_cells(
    *, rare_p: float, rare_depth_determined: bool, common_p: float = 0.5, n: int = 8000
):
    """Cells in two arms, one rare and one common, with control over whether the
    rare arm's label is a pure function of depth."""
    rng = np.random.default_rng(11)
    arm = np.where(np.arange(n) < n // 2, "rare_arm", "common_arm")
    depth = rng.lognormal(np.log(15000), 0.8, n)
    mature = np.zeros(n, dtype=bool)

    rare = np.where(arm == "rare_arm")[0]
    k = max(1, int(rare_p * len(rare)))
    if rare_depth_determined:
        mature[rare[np.argsort(depth[rare])[:k]]] = True
    else:
        mature[rng.permutation(rare)[:k]] = True

    common = np.where(arm == "common_arm")[0]
    mature[rng.permutation(common)[: int(common_p * len(common))]] = True
    return depth, mature, arm


def test_a_correlation_is_never_normalised_by_another_arms_bound():
    """The mispairing: `max(rhos)` and `max(ceilings)` taken independently.

    One arm is rare and perfectly depth-determined, so it sits exactly at its
    own ceiling. The other is common and depth-independent, so it contributes a
    large ceiling and no correlation. Dividing the first arm's rho by the second
    arm's ceiling reports 0.14 — "14% of the available headroom, clean" — for a
    label that is a *pure function of depth*.

    Paired within the arm it reports 1.0, which is what it is.
    """
    depth, mature, arm = _two_arm_cells(rare_p=0.005, rare_depth_determined=True)
    report = depth_confound_report(depth, mature, arm)

    assert report["worst_arm"] == "rare_arm"
    assert report["rho_vs_ceiling"] == pytest.approx(1.0, abs=0.05)

    # The value the old code produced, reconstructed here so the regression is
    # visible rather than asserted in the abstract.
    other_arm_ceiling = report["per_arm"]["common_arm"]["max_attainable_rho"]
    mispaired = report["worst_within_arm_rho"] / other_arm_ceiling
    assert mispaired < 0.2
    assert report["rho_vs_ceiling"] > 4 * mispaired

    # And the verdict flips with it: the arm carrying the correlation could
    # never have reached the tolerance.
    assert not report["tolerance_is_reachable"]
    assert "NOT TESTABLE" in report["reading"]


def test_each_arm_carries_its_own_bound_and_prevalence():
    """Nothing downstream should have to pair a statistic with its bound."""
    depth, mature, arm = _two_arm_cells(rare_p=0.005, rare_depth_determined=True)
    report = depth_confound_report(depth, mature, arm)

    for name, row in report["per_arm"].items():
        expected = max_attainable_rho(row["mature_share"])
        assert row["max_attainable_rho"] == pytest.approx(expected), name
        if np.isfinite(row["rho_depth_vs_mature"]) and expected > 0:
            assert row["rho_vs_ceiling"] == pytest.approx(
                abs(row["rho_depth_vs_mature"]) / expected
            ), name


def test_an_undefined_arm_is_named_rather_than_dropped():
    """An arm at prevalence exactly 0 has no correlation. That is UNDEFINED, not
    clean, and it must not vanish silently into the surviving arm's verdict.

    This is the paper's §3 failure: 58 rows reported as "label does not track
    depth" because an undefined statistic had nowhere to land except the
    negative verdict.
    """
    depth, mature, arm = _two_arm_cells(rare_p=0.005, rare_depth_determined=True)
    mature[arm == "rare_arm"] = False  # prevalence exactly 0 in that arm

    report = depth_confound_report(depth, mature, arm)
    assert report["arms_with_undefined_rho"] == ["rare_arm"]
    assert report["n_arms"] == 2
    assert "rare_arm" in report["reading"]
    assert "excluded rather than scored clean" in report["reading"]


def test_every_arm_undefined_abstains_rather_than_reporting_clean():
    depth, mature, arm = _two_arm_cells(rare_p=0.005, rare_depth_determined=True)
    mature[:] = False  # no variance anywhere

    report = depth_confound_report(depth, mature, arm)
    assert not np.isfinite(report["worst_within_arm_rho"])
    assert report["reading"].startswith("UNDEFINED")
    assert not report["maturity_tracks_depth"]


@pytest.mark.parametrize("p", [0.001, 0.005, 0.013])
def test_the_check_cannot_fire_below_the_crossing_prevalence(p):
    """The bound, as a test. Below p = 1.3516% a tolerance of 0.20 sits above
    sqrt(3p(1-p)), so the diagnostic reports clean on a label that is a pure
    function of depth. This asserts the blindness rather than describing it."""
    depth, mature, arm = _two_arm_cells(rare_p=p, rare_depth_determined=True)
    report = depth_confound_report(depth, mature, arm)

    rare = report["per_arm"]["rare_arm"]
    assert rare["max_attainable_rho"] < MATURITY_DEPTH_RHO_TOLERANCE
    assert rare["rho_vs_ceiling"] == pytest.approx(1.0, abs=0.1)  # perfectly determined
    assert not rare["tolerance_is_reachable"]
    assert not report["maturity_tracks_depth"]  # ...and it reports clean


def test_the_check_does_fire_when_prevalence_permits_it():
    """The positive control. Without this the test above is satisfied by a
    diagnostic that never fires on anything."""
    depth, mature, arm = _two_arm_cells(rare_p=0.30, rare_depth_determined=True)
    report = depth_confound_report(depth, mature, arm)

    assert report["per_arm"]["rare_arm"]["tolerance_is_reachable"]
    assert report["maturity_tracks_depth"]
    assert "NOT TESTABLE" not in report["reading"]


# ---------------------------------------------------------------------------
# The calibration routine
# ---------------------------------------------------------------------------


def _sweep(n_replicates: int = 40, *, abstain: int = 0) -> pd.DataFrame:
    """A sweep whose intervals all cover truth and all exclude zero, with the
    first ``abstain`` replicates in the lowest bin holding no interval."""
    rng = np.random.default_rng(3)
    rows = []
    for i in range(n_replicates):
        for n_cells in (10, 200):
            truth = -5.0
            wide = n_cells == 10
            rows.append(
                {
                    "replicate": i,
                    "arm": PREREGISTERED.arm,
                    "shift": PREREGISTERED.detectable_shift,
                    "n_cells_mature": n_cells + rng.integers(0, 3),
                    "intrinsic_true_parametric": truth,
                    "ci_low": truth - 1.0,
                    "ci_high": truth + 1.0 if not wide else truth + 1.0,
                }
            )
    frame = pd.DataFrame(rows)
    if abstain:
        low = frame.index[frame["n_cells_mature"] < 100][:abstain]
        frame.loc[low, ["ci_low", "ci_high"]] = np.nan
    return frame


def test_an_abstention_is_not_counted_as_a_coverage_failure():
    """Comparisons against NaN return False, so a replicate holding no interval
    used to score as a coverage AND a discrimination failure — halving both
    rates in the bin where refusing was most often right.

    Here every replicate that *does* hold an interval covers truth and excludes
    zero, so every rate must be exactly 1.0 whether or not abstentions are
    present. If an abstention leaks into a denominator, a rate drops below 1.
    """
    clean = coverage_and_discrimination(_sweep())
    half_refused = coverage_and_discrimination(_sweep(abstain=20))

    for table, label in ((clean, "clean"), (half_refused, "with abstentions")):
        rates = table[["coverage", "discrimination"]].to_numpy(dtype=float)
        present = rates[np.isfinite(rates)]
        assert len(present) > 0, label
        assert present == pytest.approx(1.0), label

    # The abstentions are visible as a count rather than absorbed into a rate.
    assert clean["n_abstained"].sum() == 0
    assert half_refused["n_abstained"].sum() == 20
    assert (
        half_refused["n_estimated"] + half_refused["n_abstained"]
        == half_refused["n_replicates"]
    ).all()


def test_a_bin_where_every_replicate_abstained_gets_its_own_verdict():
    """"Nobody answered" and "the answers were bad" are different findings.
    NaN >= target is False, so the old code returned `not_estimable` for both."""
    frame = _sweep()
    frame.loc[frame["n_cells_mature"] < 100, ["ci_low", "ci_high"]] = np.nan

    table = coverage_and_discrimination(frame)
    empty = table[table["n_estimated"] == 0]
    assert not empty.empty, "the low bins should hold no intervals at all"
    assert (empty["verdict"] == "all_abstained").all()
    assert not np.isfinite(empty["coverage"]).any()

    # The bins that did answer are untouched and still judged normally.
    answered = table[table["n_estimated"] > 0]
    assert not answered.empty
    assert (answered["verdict"] != "all_abstained").all()


# ---------------------------------------------------------------------------
# The leakage invariant
# ---------------------------------------------------------------------------


def test_the_leakage_guard_fires_across_identifier_spaces():
    """Appendix A item 1 and 2: the guard compared gene SYMBOLS against database
    ACCESSIONS. An empty intersection follows by construction, so it passed on
    every input including violating ones.

    The violating input is a label set carrying the target gene in a *different*
    namespace than the target list uses. Without an alias map that comparison is
    vacuous, and `resolve_targets` is what makes it non-vacuous.
    """
    with pytest.raises(LeakageError):
        assert_no_target_leakage(
            ["GUCA2A", "LGR5"], ["GUCA2A"], context="same-namespace control"
        )

    # The same violation, expressed in Ensembl ids. This is the input the old
    # guard returned clean on.
    with pytest.raises((LeakageError, Exception)):
        assert_no_target_leakage(
            ["ENSG00000197273", "ENSG00000139292"],
            ["GUCA2A"],
            context="cross-namespace violation",
            alias_map={"GUCA2A": "ENSG00000197273"},
        )


def test_the_leakage_guard_still_passes_a_clean_panel():
    """The positive control: a guard that raises on everything is no better."""
    assert_no_target_leakage(["LGR5", "ASCL2"], ["GUCA2A"], context="clean")


# ---------------------------------------------------------------------------
# The provenance stamp
# ---------------------------------------------------------------------------


def test_the_dirtiness_check_counts_untracked_files(tmp_path, monkeypatch):
    """Appendix A item 5: the check ignored untracked files, so a table produced
    by a script that was never committed recorded a clean hash pointing at a
    commit that did not contain the script.

    The violating input is a repository whose only difference from HEAD is an
    untracked producing script.
    """
    import subprocess

    from src.common import provenance

    def git(*args):
        subprocess.run(
            ["git", *args], cwd=tmp_path, check=True, capture_output=True, text=True
        )

    git("init", "-q")
    git("config", "user.email", "t@example.com")
    git("config", "user.name", "t")
    (tmp_path / "tracked.txt").write_text("committed\n")
    git("add", "tracked.txt")
    git("commit", "-qm", "initial")

    monkeypatch.setattr(provenance, "REPO_ROOT", tmp_path)
    assert not provenance.git_is_dirty(), "a clean tree must read clean"

    # The input that made the old check return clean.
    (tmp_path / "produce_the_table.py").write_text("# never committed\n")
    assert provenance.untracked_files(), "an untracked producer must be listed"
    assert provenance.git_is_dirty(), (
        "a tree whose only change is an untracked producing script is NOT clean "
        "— that is the failure this check exists to catch"
    )


# ---------------------------------------------------------------------------
# The committed results themselves
#
# Everything above tests a guard against a synthetic input. These two walk the
# committed artefacts instead, because both invariants they cover are enforced
# by whichever writer a job happened to call rather than by anything that
# inspects the result afterwards. `schema.write_results` validates invariant 1
# and has one production call site; `io.write_versioned_table` checks the seed
# and the tree but never looks at the estimability column, and has twenty-six.
# A guard reachable only by opting into it is the shape this file exists to
# reject.
# ---------------------------------------------------------------------------


def _committed_tables() -> list[Path]:
    """Every committed parquet carrying the columns invariant 1 constrains."""
    out = []
    for path in sorted(RESULTS_DIR.glob("*/*.parquet")):
        try:
            names = set(pq.read_schema(path).names)
        except Exception:  # not a readable parquet; a different test's problem
            continue
        if "estimability" in names and "intrinsic" in names:
            out.append(path)
    return out


def test_no_committed_result_writes_a_number_where_it_means_not_estimable():
    """CLAUDE.md invariant 1, checked against the results rather than the writer.

    `None` is not `0.0`: a downstream average cannot tell "we measured zero"
    from "we did not ask". The assertion enforcing that lives in
    `schema.write_results`, which one production job calls. Every other job
    writes through `io.write_versioned_table`, which never looks at the column.
    So the invariant is currently a property of the code path, and this makes it
    a property of the artefacts.

    The invariant is asserted on every table carrying the two columns, not only
    on frames matching the full results schema. `threshold_sweep` carries them
    without `granularity_rung` or the CI pair, and it is exactly as capable of
    writing a zero where it means "did not ask"; an earlier version of this test
    ran `validate_results` and so skipped it with a schema complaint instead of
    checking the thing that matters.
    """
    tables = _committed_tables()
    assert tables, (
        "no committed table carries an estimability column, so this check "
        "cannot fail and is therefore not a check. If the results moved, point "
        "it at where they went."
    )
    for path in tables:
        frame = pd.read_parquet(path)
        where = path.relative_to(RESULTS_DIR.parent)

        refused = frame["estimability"] == "not_estimable"
        wrote_a_number = refused & frame["intrinsic"].notna()
        assert not wrote_a_number.any(), (
            f"{where}: {int(wrote_a_number.sum())} not_estimable row(s) carry a "
            f"number. None is not 0.0 (invariant 1)."
        )

        answered = frame["estimability"] == "ok"
        lost = answered & frame["intrinsic"].isna()
        assert not lost.any(), (
            f"{where}: {int(lost.sum())} row(s) marked 'ok' have no estimate. "
            f"Either the estimate exists or the row is not 'ok'."
        )

        # Tables that are full results frames get the whole frozen contract too.
        if set(REQUIRED_COLUMNS).issubset(frame.columns):
            try:
                validate_results(coerce_results(frame))
            except SchemaViolation as exc:  # pragma: no cover - the failure is the point
                pytest.fail(f"{where}: {exc}")


def test_the_invariant_1_check_fires_on_a_zero_standing_in_for_none():
    """The violating input: one not-estimable row carrying 0.0 instead of null.

    This is the single most likely route to a wrong conclusion in this project,
    and it is one keystroke from correct.
    """
    tables = _committed_tables()
    frame = coerce_results(pd.read_parquet(tables[0]))
    not_est = frame.index[frame["estimability"] == "not_estimable"]
    assert len(not_est), "need a not_estimable row to corrupt"

    frame.loc[not_est[0], "intrinsic"] = 0.0
    with pytest.raises(SchemaViolation, match="None is not 0.0"):
        validate_results(frame)


#: Result sidecars whose commit no longer resolves, with the reason. This is a
#: ratchet, not an amnesty: a new dead sha fails the test, and an entry that
#: starts resolving again must be deleted (the test checks that too). The
#: mechanism is structural rather than careless — CONTRIBUTING §3 asks for a
#: rebase before review, provenance stamps the pre-rebase hash, and the hash
#: dies at merge. Every result written on a branch is exposed to it.
KNOWN_UNRESOLVABLE_SHAS: dict[str, str] = {
    # Empty, and that is the point. The one entry this held --
    # e5ebdc330a66, depth_confound_reference, written on a branch that was
    # rebased before review -- was closed by re-deriving the table on a sha
    # that resolves. The re-run came back bit-identical, so the numbers were
    # never in doubt; only the record of what produced them was. A new entry
    # here needs the same reason written down, and the ratchet below deletes
    # any entry that starts resolving again.
}


def _git(*args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        ["git", *args], cwd=REPO_ROOT, capture_output=True, text=True
    )


def _sha_resolves(sha: str) -> bool:
    return _git("cat-file", "-t", sha).returncode == 0


def _sidecar_shas() -> list[tuple[Path, str]]:
    out = []
    for path in sorted(RESULTS_DIR.glob("*/*.meta.json")):
        meta = json.loads(path.read_text())
        sha = meta.get("git_sha") or meta.get("git_sha_short")
        if sha:
            out.append((path, sha))
    return out


def test_every_result_sidecar_names_a_commit_that_exists():
    """CLAUDE.md invariant 10 says every result carries the git sha that
    produced it. Nothing checked that the sha still resolves, so the invariant
    could be satisfied by a string.

    This is the sibling of Appendix A item 5. That defect was a stamp blind to
    uncommitted code and its repair counts untracked files; this is the same
    invariant defeated by a second route the repair does not cover, and one
    number in the paper travels through it.
    """
    if _git("rev-parse", "--is-shallow-repository").stdout.strip() == "true":
        pytest.skip(
            "shallow clone: history is absent, so every historical sha would "
            "read as dead. CI checks out with fetch-depth: 0 so this runs there."
        )

    dead = [(p, s) for p, s in _sidecar_shas() if not _sha_resolves(s)]
    unexpected = [(p, s) for p, s in dead if s not in KNOWN_UNRESOLVABLE_SHAS]
    assert not unexpected, "result sidecars naming commits that do not exist:\n" + "\n".join(
        f"  {p.relative_to(REPO_ROOT)} -> {s}" for p, s in unexpected
    )

    # The ratchet only tightens: an entry that resolves again is stale bookkeeping.
    resurrected = [s for s in KNOWN_UNRESOLVABLE_SHAS if _sha_resolves(s)]
    assert not resurrected, (
        f"{resurrected} resolve now — delete them from KNOWN_UNRESOLVABLE_SHAS "
        f"rather than carrying an exemption nothing needs"
    )


def test_the_sha_check_fires_on_a_commit_that_was_never_written():
    """The violating input: a well-formed hash naming nothing.

    A sidecar recording this would satisfy invariant 10 as written — the field
    is populated — and identify no code at all.
    """
    fabricated = "0" * 40
    assert not _sha_resolves(fabricated), "a fabricated sha must not resolve"
    assert _sha_resolves("HEAD"), "the check must still pass a commit that exists"


# ---------------------------------------------------------------------------
# The provenance record itself
# ---------------------------------------------------------------------------


def test_extra_meta_cannot_overwrite_the_provenance_record(tmp_path):
    """The violating input is the one that actually shipped.

    `write_versioned_table` merged `extra_meta` straight over the provenance
    record, so a caller could rewrite the fields invariant 10 exists to record.
    Three committed GSE39582 tables carry `platform: "GPL570"` for exactly this
    reason -- the GEO platform overwrote the OS string, and those sidecars no
    longer say what machine produced them.

    The same hole accepted `git_sha`, which would have let a table claim a
    commit it was not built from -- and the sha guard above would then have
    validated the claim.
    """
    frame = pd.DataFrame({"a": [1, 2]})

    with pytest.raises(ReservedMetaKeyError, match="platform"):
        write_versioned_table(
            frame, "reserved_key_probe", seed=1, results_dir=tmp_path,
            allow_dirty=True, extra_meta={"platform": "GPL570"},
        )

    with pytest.raises(ReservedMetaKeyError, match="git_sha"):
        write_versioned_table(
            frame, "reserved_key_probe", seed=1, results_dir=tmp_path,
            allow_dirty=True, extra_meta={"git_sha": "0" * 40},
        )

    # A caller annotating the result rather than rewriting its record is fine.
    path = write_versioned_table(
        frame, "reserved_key_probe", seed=1, results_dir=tmp_path,
        allow_dirty=True, extra_meta={"geo_platform": "GPL570", "cohort": "GSE39582"},
    )
    meta = json.loads(path.with_suffix(".meta.json").read_text())
    assert meta["geo_platform"] == "GPL570"
    assert meta["platform"] != "GPL570", "the OS record must survive annotation"


def test_no_job_hardcodes_the_dirty_tree_bypass():
    """`allow_dirty` is a scratch-run escape hatch, not a default.

    Every job in `src/bulk/` used to pass `allow_dirty=True` as a literal, so
    the bulk arm could not write a clean provenance stamp even from a spotless
    tree. That is why all fifteen committed bulk tables record
    `git_dirty: true`: not carelessness at the keyboard, but a guard switched
    off at the call site. Re-running them on a clean tree would have produced
    dirty stamps again.

    The flag must come from the caller. A literal `True` here is the defect.
    """
    offenders = []
    for path in sorted(Path("src").rglob("*.py")):
        tree = ast.parse(path.read_text(), filename=str(path))
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            for kw in node.keywords:
                if kw.arg == "allow_dirty" and isinstance(kw.value, ast.Constant):
                    if kw.value.value is True:
                        offenders.append(f"{path}:{kw.value.lineno}")
    assert not offenders, (
        "allow_dirty=True is hardcoded at:\n  " + "\n  ".join(offenders)
        + "\nThread it from the CLI instead; a job that cannot write a clean "
          "stamp makes invariant 10 unsatisfiable for everything it produces."
    )


def test_the_manifest_records_portable_paths_not_one_machine_s(tmp_path, monkeypatch):
    """The violating input: an absolute path under a data directory that lives
    outside the repository.

    `data/manifest.csv` is, in its own words, the only record of the raw data
    that travels. `BRP_DATA_DIR` routinely points elsewhere -- on the cluster it
    is an absolute path on the project filesystem -- and `append_manifest_row`
    recorded whatever it was handed. Five rows naming one machine's scratch disk
    reached a commit that way.
    """
    from src.reference import ingest

    monkeypatch.setattr(ingest, "DATA_DIR", tmp_path / "data")
    absolute = tmp_path / "data" / "raw" / "gse39582" / "series.txt.gz"
    absolute.parent.mkdir(parents=True)
    absolute.write_text("x")

    assert ingest.manifest_key(absolute) == "data/raw/gse39582/series.txt.gz", (
        "an absolute path under BRP_DATA_DIR must be recorded as data/..."
    )
    assert ingest.manifest_key("data/raw/lee/already_relative.txt") == \
        "data/raw/lee/already_relative.txt"

    # A scratch path outside every known root passes through unchanged: tests
    # write there legitimately. The portability guarantee is enforced on the
    # committed artefact by test_manifest.py, which is what caught the real bug.
    outside = tmp_path / "somewhere_else" / "stray.txt"
    outside.parent.mkdir(parents=True)
    outside.write_text("x")
    assert ingest.manifest_key(outside) == str(outside)


def test_the_lee_loader_reads_the_configured_data_directory(tmp_path, monkeypatch):
    """The violating input: BRP_DATA_DIR pointing outside the checkout.

    `load_lee_cohort` built its default path from REPO_ROOT/"data" rather than
    RAW_DIR, so it ignored BRP_DATA_DIR entirely. On a laptop with the variable
    unset the two are the same directory and nothing shows; on the cluster it
    went looking for Lee under the repository, found an almost-empty data/, and
    died -- twice, forty seconds into a job whose inputs were all present on the
    project filesystem the wrapper had just verified.

    A default that is correct only where two paths coincide is not a default.
    """
    from src.common import paths
    from src.estimator import lee_io

    elsewhere = tmp_path / "project" / "data" / "raw"
    monkeypatch.setattr(paths, "RAW_DIR", elsewhere)

    with pytest.raises((FileNotFoundError, OSError)) as caught:
        lee_io.load_lee_cohort("smc", target_genes=["GUCA2A"])

    named = str(caught.value)
    assert str(elsewhere) in named, (
        f"the loader must look under the configured data directory; it looked "
        f"at {named!r}"
    )
    assert "BetterRiseProject/data/raw" not in named, (
        "it fell back to the repository checkout despite BRP_DATA_DIR"
    )


def test_the_leakage_guard_refuses_a_positional_index():
    """The violating input is the committed S matrix loaded the obvious way.

    `_identifier_space` was a binary test: matches the Ensembl pattern, or else
    it is a symbol. A DataFrame read from parquet carries a RangeIndex, so
    `S.index` is [0, 1, 2, ...] -- not Ensembl, therefore "symbol", compared
    against symbol targets, intersecting nothing, passing.

    That is issue #35's vacuous pass surviving the fix for issue #35, in the
    identifier space the committed `S_matrix_lineage_1.0.0.parquet` actually
    presents: its gene ids are in a COLUMN, and nothing made the caller notice.
    A classifier with no reject option always answers, and its answer is
    worthless exactly where it matters.
    """
    from src.reference.signature import LeakageGuardError, _identifier_space

    assert _identifier_space([0, 1, 2]) == "unrecognised"
    assert _identifier_space(["0", "1"]) == "unrecognised"
    assert _identifier_space([]) == "unrecognised"
    # Real identifiers must keep working, hyphens and digits included.
    assert _identifier_space(["ENSG00000000971"]) == "ensembl"
    assert _identifier_space(["GUCA2A", "MT-CO1", "C1orf43", "HLA-A"]) == "symbol"

    with pytest.raises(LeakageGuardError, match="not gene identifiers"):
        assert_no_target_leakage(pd.RangeIndex(800), ["GUCA2A"], context="probe")


def test_the_committed_s_matrix_cannot_be_checked_by_accident():
    """End to end on the artefact itself: load it as anyone would, and the
    guard must refuse rather than certify."""
    import yaml

    from src.reference.signature import LeakageGuardError

    path = RESULTS_DIR / "2026-08-26_63ead2e" / "S_matrix_lineage_1.0.0.parquet"
    if not path.exists():                      # pragma: no cover - artefact moved
        pytest.skip(f"{path} absent")
    signature = pd.read_parquet(path)
    targets = sorted(yaml.safe_load(
        (REPO_ROOT / "config" / "panel.yaml").read_text())["tiers"]["A"]["genes"])

    with pytest.raises(LeakageGuardError):
        assert_no_target_leakage(signature.index, targets, context="the bake-off signature")

    # Promoted to the identifier column it actually carries, the check runs.
    assert_no_target_leakage(
        signature.set_index("gene").index, targets, context="probe",
        alias_map={t: f"ENSG_ABSENT_{t}" for t in targets},
    )


# ---------------------------------------------------------------------------
# Cross-gene comparison of detection deltas
#
# The defect: `specificity()` compared six genes' detection deltas to each
# other. A proportion's sensitivity depends on its baseline -- the panel spans
# 0.36 to 0.98 -- so the ranking was substantially a ranking by abundance, and
# the guard could not have noticed because its fixture named `delta_detect`
# directly and had no baseline rate anywhere in it.
#
# These are the inputs that force it. The ground truth in each is a UNIFORM
# thinning: one fold change applied to every gene, no gene-specificity at all.
# ---------------------------------------------------------------------------


def _thinned_panel(spec, *, n_patients=30, n_cells=400, noise=0.02, seed=11):
    """`{gene: (baseline_detection, fold_change)}` -> per-patient rows.

    Detection is DERIVED from a baseline and a thinning, never named. That is
    the whole point: a test can hold the fold change fixed across genes and vary
    only the baseline, which is what "the gradient might be abundance" means and
    what the previous fixture could not express.
    """
    from src.reference.detection_scale import delta_cloglog

    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_patients):
        for gene, (baseline, fold) in spec.items():
            p_n = float(np.clip(baseline + rng.normal(0, noise), 0.01, 0.99))
            p_t = float(np.clip(1.0 - (1.0 - p_n) ** fold + rng.normal(0, noise),
                                0.001, 0.999))
            rows.append({
                "granularity_rung": "lineage", "patient_id": f"p{i}", "gene": gene,
                "n_normal": n_cells, "n_tumour": n_cells,
                "detect_normal": p_n, "detect_tumour": p_t,
                "delta_detect": p_t - p_n,
                "log2_cp10k_ratio": float(np.log2(fold) + rng.normal(0, noise)),
            })
    frame = pd.DataFrame(rows)
    frame["delta_cloglog"] = delta_cloglog(frame)
    return frame


def test_a_uniform_thinning_is_not_reported_as_gene_specific():
    """One fold change, six baselines, no biology -- and detection sees a tier.

    This is the adenoma reading's defect with the answer known in advance. Every
    gene here is thinned by exactly 0.75, so every true pairwise difference is
    zero. The detection deltas nonetheless spread from -0.03 to -0.10 purely
    because the genes sit at different baselines, and read across genes that is
    a gradient with a story attached to it.
    """
    from src.reference.detection_scale import (
        UNIFORM_THINNING_R2_CEILING,
        uniform_thinning_null,
    )
    from src.reference.jobs.icbi_coexpression import specificity

    panel = _thinned_panel({
        "ACTB": (0.985, 0.75), "KRT8": (0.955, 0.75), "EPCAM": (0.90, 0.75),
        "CDX2": (0.82, 0.75), "GUCA2A": (0.44, 0.75), "MS4A12": (0.36, 0.75),
    })
    contrasts = specificity(panel, seed=1)

    detection = contrasts[contrasts["statistic"] == "detection"]
    assert detection["excludes_zero"].sum() >= 20, (
        "the input no longer forces the defect: detection is supposed to report "
        "most of these 30 contrasts as real when none of them are"
    )

    load_bearing = contrasts[contrasts["load_bearing"]]
    assert load_bearing["excludes_zero"].sum() <= 6, (
        "the load-bearing scale is inheriting the abundance confound"
    )
    # The contrast the reading actually turns on, between two genes far apart in
    # abundance, must not survive on the load-bearing scale.
    target = load_bearing.set_index("contrast").loc["GUCA2A - CDX2"]
    assert not target["excludes_zero"]

    _, verdict = uniform_thinning_null(panel)
    assert verdict["verdict"] == "GRADIENT IS ABUNDANCE"
    assert verdict["variance_explained"] >= UNIFORM_THINNING_R2_CEILING


def test_a_saturated_control_is_not_compared_naively_on_detection():
    """ACTB at 0.99 cannot fall, so everything falls further than it.

    `SATURATION_CEILING` already exists because of this, and `premise_holds`
    switches a saturated control to log2 expression. The specificity table never
    switched, so `GUCA2A - ACTB` excluding zero was guaranteed by ACTB's
    abundance rather than earned by GUCA2A's behaviour.
    """
    from src.reference.jobs.icbi_coexpression import specificity

    panel = _thinned_panel(
        {"ACTB": (0.99, 0.70), "GUCA2A": (0.45, 0.70), "CDX2": (0.82, 0.70)},
        noise=0.015, seed=7,
    )
    got = specificity(panel, seed=2)
    row = got[(got["statistic"] == "detection")
              & (got["contrast"] == "GUCA2A - ACTB")].iloc[0]
    assert row["excludes_zero"], (
        "the input no longer forces the defect -- detection is supposed to call "
        "this real when both genes were thinned identically"
    )
    for statistic in ("cloglog", "log2_cp10k"):
        clean = got[(got["statistic"] == statistic)
                    & (got["contrast"] == "GUCA2A - ACTB")].iloc[0]
        assert not clean["excludes_zero"], (
            f"{statistic} inherited the saturation artefact"
        )


def test_a_table_without_the_abundance_free_columns_claims_nothing_load_bearing():
    """The failure mode of a fix like this is quietly reverting to the old scale.

    A caller handing in a frame that carries only `delta_detect` -- every table
    committed before this change -- must get diagnostic rows and no load-bearing
    verdict, rather than detection silently being promoted because it is the
    only column present.
    """
    from src.reference.jobs.icbi_coexpression import specificity

    panel = _thinned_panel({"ACTB": (0.98, 0.7), "GUCA2A": (0.45, 0.7),
                            "CDX2": (0.8, 0.7)})
    old_shape = panel.drop(columns=["delta_cloglog", "log2_cp10k_ratio"])
    got = specificity(old_shape, seed=3)
    assert not got.empty
    assert set(got["statistic"]) == {"detection"}
    assert not got["load_bearing"].any()
    assert set(got["standing"]) == {"diagnostic only -- NOT read across genes"}


def test_the_boundary_rule_is_a_function_of_cell_count():
    """A 15-cell zero and a 500-cell zero must not claim the same evidence.

    cloglog is undefined at 0 and 1, and the obvious repair is a fixed epsilon.
    That would let an arm with 15 cells and an arm with 500 report the same
    number from the same observed rate, which is the "a check that cannot tell
    them apart" shape in a transform rather than a guard.
    """
    from src.reference.detection_scale import cloglog_rate

    sparse = float(cloglog_rate(0.0, 15))
    dense = float(cloglog_rate(0.0, 500))
    assert np.isfinite(sparse) and np.isfinite(dense)
    assert dense < sparse, (
        "a zero over 500 cells must sit lower than a zero over 15 -- it is the "
        "stronger evidence of absence"
    )
    assert float(cloglog_rate(1.0, 500)) > float(cloglog_rate(1.0, 15))


def test_the_thinning_null_abstains_when_baselines_do_not_spread():
    """One free parameter over six near-identical points explains everything.

    If every gene sits at the same abundance the null has nothing to
    discriminate on, and an R² from that fit would be an artefact reported as a
    verdict. The diagnostic has to say UNDEFINED rather than return clean.
    """
    from src.reference.detection_scale import (
        MIN_BASELINE_SPREAD,
        uniform_thinning_null,
    )

    flat = _thinned_panel({
        "ACTB": (0.50, 0.75), "KRT8": (0.51, 0.75), "EPCAM": (0.49, 0.75),
        "CDX2": (0.50, 0.75), "GUCA2A": (0.51, 0.60), "MS4A12": (0.49, 0.60),
    }, noise=0.005)
    table, verdict = uniform_thinning_null(flat)
    assert verdict["verdict"] == "UNDEFINED"
    assert table.empty
    assert verdict["baseline_spread"] < MIN_BASELINE_SPREAD


# ---------------------------------------------------------------------------
# An interval that overstates its own confidence
#
# The defect: every reading in this project ends in a percentile bootstrap over
# patients, and at the n a small stratum runs at that interval is narrower than
# it claims -- 0.82x the correct width at n=10, 0.53x at n=4. Nothing raises. A
# verdict reads "excludes zero" at a real false-positive rate of 9.6% or 18.8%
# against a nominal 5%.
#
# It is the repository's own defect one layer up: a check that cannot fail
# turns absence of evidence into a green light, and an interval narrower than
# it claims turns noise into a finding. Both are silent.
#
# The second defect, which is subtler and which HAPPENED: quoting power from
# one interval method while planning to report another's. The percentile
# bootstrap's power at 50% silencing on the MLH1 cohort is 86%; the Student-t
# interval that would actually be reported gives 74% on the same generator.
# Taking the first number and the second interval overstates the design, and
# there is no error message for it because both numbers are correct about
# something.
# ---------------------------------------------------------------------------


def _mlh1_cohort(n_patients: int = 10):
    return np.full(n_patients, 262), np.full(n_patients, 8000.0)


def test_the_calibration_check_reports_miscalibrated_on_the_interval_in_use():
    """The percentile bootstrap at n=10, under a null with no effect in it.

    THE INPUT THAT FORCES THE VERDICT. A calibration function that returned
    CALIBRATED on everything would be indistinguishable from this one on clean
    data, and clean data is all the repository had until this was measured. The
    generator here has a fold change of exactly 1.0 -- the two arms are the same
    distribution -- so every exclusion of zero is a false positive by
    construction.
    """
    from src.reference.interval_calibration import (
        calibration_verdict,
        rejection_rate,
        simulate_deltas,
    )

    n_cells, depth = _mlh1_cohort()

    def no_effect_at_all(rng):
        return simulate_deltas(
            n_cells=n_cells, depth=depth, cp10k=3.0,
            fold_change=1.0, tau=0.0, rng=rng,
        )

    rate = rejection_rate(
        no_effect_at_all, method="percentile",
        rng=np.random.default_rng(4), n_trials=600,
    )
    assert rate > 0.05, (
        "the percentile bootstrap over-rejects at n=10 and this fixture must "
        "reproduce that, not assert it"
    )
    assert calibration_verdict(rate) == "MISCALIBRATED"


def test_the_calibration_check_still_passes_a_calibrated_interval():
    """The other half: it must not call everything miscalibrated.

    A verdict function that returned MISCALIBRATED on all input would 'catch'
    the defect above and be worth nothing.
    """
    from src.reference.interval_calibration import (
        CALIBRATED_METHOD,
        calibration_verdict,
        rejection_rate,
        simulate_deltas,
    )

    n_cells, depth = _mlh1_cohort()
    rate = rejection_rate(
        lambda rng: simulate_deltas(
            n_cells=n_cells, depth=depth, cp10k=3.0,
            fold_change=1.0, tau=0.0, rng=rng,
        ),
        method=CALIBRATED_METHOD, rng=np.random.default_rng(4), n_trials=600,
    )
    assert calibration_verdict(rate) == "CALIBRATED"


def test_power_quoted_from_one_interval_beside_another_is_refused():
    """The defect that happened, as the table that carries it.

    Power 0.86 came from the percentile bootstrap; the design reports the
    Student-t interval, whose false-positive rate is on the same row. Every
    individual number here is real. The row is still not a claim about a
    design, because the two halves are about different intervals.
    """
    from src.reference.interval_calibration import (
        CalibrationError,
        check_power_carries_its_own_calibration,
    )

    mixed = pd.DataFrame([
        {"cohort": "mlh1_methylated", "tau": 0.2, "cp10k_normal": 0.039,
         "method": "percentile", "power": 0.86, "false_positive_rate": 0.096},
        {"cohort": "mlh1_methylated", "tau": 0.2, "cp10k_normal": 0.039,
         "method": "student_t", "power": 0.74, "false_positive_rate": 0.050},
    ])
    with pytest.raises(CalibrationError, match="mixes interval methods"):
        check_power_carries_its_own_calibration(mixed)


def test_one_method_carrying_two_false_positive_rates_is_refused():
    """Same method, two rates: they were measured on different generators.

    This is what assembling a table from a simulation someone ran last week and
    one run today looks like, and the power figures beside them are then not
    comparable to each other.
    """
    from src.reference.interval_calibration import (
        CalibrationError,
        check_power_carries_its_own_calibration,
    )

    stitched = pd.DataFrame([
        {"cohort": "c", "tau": 0.0, "cp10k_normal": 0.039, "method": "student_t",
         "power": 0.79, "false_positive_rate": 0.045},
        {"cohort": "c", "tau": 0.0, "cp10k_normal": 0.039, "method": "student_t",
         "power": 0.99, "false_positive_rate": 0.061},
    ])
    with pytest.raises(CalibrationError, match="different false-positive rates"):
        check_power_carries_its_own_calibration(stitched)


def test_a_power_table_without_its_calibration_is_refused():
    from src.reference.interval_calibration import (
        CalibrationError,
        check_power_carries_its_own_calibration,
    )

    bare = pd.DataFrame([{"cohort": "c", "tau": 0.0, "cp10k_normal": 0.039,
                          "method": "student_t", "power": 0.86}])
    with pytest.raises(CalibrationError, match="false_positive_rate"):
        check_power_carries_its_own_calibration(bare)


def test_heterogeneity_is_not_the_spread_of_the_deltas():
    """The naive implementation, forced to be wrong.

    ``tau`` is the between-patient variation NET of binomial sampling. The
    obvious implementation -- the standard deviation of the per-patient deltas
    -- reports a large number on data generated with no patient-to-patient
    variation whatsoever, because sampling noise alone produces spread. A power
    calculation fed that number would be pessimistic; fed zero where there IS
    heterogeneity it would be optimistic, which is the direction that matters
    and the direction the first MLH1 power statement went.

    Here the truth is tau = 0 exactly: every patient has the same underlying
    rate in both arms, and all spread is binomial.
    """
    from src.reference.interval_calibration import heterogeneity_tau

    rng = np.random.default_rng(5)
    n_patients, cells, p = 40, 250, 0.40
    frame = pd.DataFrame({
        "gene": "SYNTH", "n_normal": cells, "n_tumour": cells,
        "detect_normal": rng.binomial(cells, p, n_patients) / cells,
        "detect_tumour": rng.binomial(cells, p, n_patients) / cells,
    })
    out = heterogeneity_tau(frame, seed=5).iloc[0]
    assert out["observed_sd"] > 0.08, "sampling noise alone must produce spread"
    assert out["tau"] == pytest.approx(0.0, abs=0.05), (
        "reporting the raw spread as tau would attribute pure binomial noise "
        "to biology"
    )


def test_heterogeneity_does_not_report_zero_when_patients_really_differ():
    """The complement: a floor at zero must not become a floor on everything."""
    from src.reference.interval_calibration import heterogeneity_tau

    rng = np.random.default_rng(6)
    n_patients, cells, p = 40, 250, 0.40
    mu = -np.log1p(-p)
    log_fc = rng.normal(0.0, 0.5, n_patients)
    frame = pd.DataFrame({
        "gene": "SYNTH", "n_normal": cells, "n_tumour": cells,
        "detect_normal": rng.binomial(cells, p, n_patients) / cells,
        "detect_tumour": rng.binomial(
            cells, 1 - np.exp(-mu * np.exp(log_fc)), n_patients) / cells,
    })
    assert heterogeneity_tau(frame, seed=6).iloc[0]["tau"] > 0.3


# ---------------------------------------------------------------------------
# Scoring a gene outside the standing panel
#
# The MLH1 positive control adds a seventh gene to a six-gene scoring path.
# Two ways that goes wrong silently: a gene with no declared role lands in a
# results table as an unlabelled row, and a gene that is already on the panel
# gets scored twice with the second copy winning every groupby downstream.
# ---------------------------------------------------------------------------


def test_a_gene_scored_without_a_declared_role_is_refused():
    from src.reference.jobs.coexpression_silencing import rows_for_patient

    n = 120
    tissue = np.array(["normal"] * (n // 2) + ["tumour"] * (n // 2))
    rng = np.random.default_rng(9)
    counts = {g: rng.poisson(2.0, n).astype(float)
              for g in ("ACTB", "KRT8", "EPCAM", "CDX2", "MS4A12", "GUCA2A")}
    counts["MLH1"] = rng.poisson(0.05, n).astype(float)
    with pytest.raises(KeyError, match="no role for"):
        rows_for_patient(
            study_id="S", patient="P", counts=counts,
            depth=np.full(n, 6000.0), tissue=tissue, seed=1,
        )


def test_an_extra_gene_already_on_the_panel_is_refused():
    """Scoring GUCA2A 'extra' would put two rows for it on every patient."""
    import inspect

    from src.reference.jobs.icbi_coexpression import study_deltas

    source = inspect.getsource(study_deltas)
    assert "already on the standing panel" in source, (
        "the overlap guard is what keeps a duplicate gene from silently "
        "winning every groupby; if it moved, this test must move with it"
    )
    with pytest.raises(Exception, match="already on the standing panel"):
        study_deltas(
            Path("/nonexistent/atlas.h5ad"), pd.DataFrame(), "Pelka_2021_Cell",
            extra_genes=("GUCA2A",),
        )


# ---------------------------------------------------------------------------
# One arm, two definitions, no check
#
# The defect, found when the cluster printed `mlh1_unmethylated  15` against a
# pre-registration that said 19. `mlh1_positive_control.arm_of` sent
# `mlh1_intact_mmrd` to its own arm, so the secondary arm was 15;
# `interval_calibration` sized the same arm as `mlh1_stratum != methylated`,
# which is 19. Both were correct answers to different questions that had been
# given one name in two files.
#
# This is not a check that could not fail. It is a quantity with two
# definitions and NO check -- which is the same failure one step earlier, and
# the reason the panel and the labelling axes are loaded from one place.
# ---------------------------------------------------------------------------


def test_the_mlh1_arms_have_exactly_one_definition():
    """No consumer may re-derive an arm from `mlh1_stratum` by hand.

    The failing input is the source itself: a module that partitions the cohort
    with its own comparison instead of importing `arm_of` is how 15 and 19 got
    into two files. Both consumers are read, so adding a third that hand-rolls
    the split fails here rather than in a document six weeks later.
    """
    consumers = [
        REPO_ROOT / "src" / "reference" / "jobs" / "mlh1_positive_control.py",
        REPO_ROOT / "src" / "reference" / "jobs" / "interval_calibration.py",
    ]
    for path in consumers:
        source = path.read_text(encoding="utf-8")
        assert "from src.reference.mlh1_arms import" in source, (
            f"{path.name} must take the arms from src.reference.mlh1_arms"
        )
        # AST, NOT TEXT. The first version of this guard grepped for the
        # string and fired on the COMMENT that explains the defect -- a check
        # that cannot tell code from prose about code, which is a false
        # positive of exactly the kind this file exists to eliminate in the
        # other direction. Comparisons are read as comparisons.
        for node in ast.walk(ast.parse(source)):
            if not isinstance(node, ast.Compare):
                continue
            operands = [node.left, *node.comparators]
            names = {n.value for n in operands
                     if isinstance(n, ast.Constant) and isinstance(n.value, str)}
            if "mlh1_methylated" in names and any(
                isinstance(op, (ast.NotEq, ast.Eq)) for op in node.ops
            ):
                raise AssertionError(
                    f"{path.name} line {node.lineno} partitions the cohort by "
                    f"comparing against 'mlh1_methylated' directly. That is how "
                    f"n=19 reached the pre-registration for an arm the reading "
                    f"job emits at n=15: the comparison counts the four "
                    f"intact-MMRd patients that are reported separately. Use "
                    f"arm_of()."
                )


def test_the_arms_are_disjoint_and_exhaust_the_cohort():
    """Three arms, every stratum in exactly one, nothing counted twice.

    A stratum landing in two arms would report the same patients under two
    standings -- one 'powered', one 'UNDERPOWERED' -- and the reader could not
    tell that the second was a subset of the first.
    """
    from src.reference.mlh1_arms import ARMS, arm_of

    strata = ("mlh1_methylated", "mlh1_intact_mmrd", "mmr_proficient",
              "mlh1_deficient_unmethylated")
    assignments = [arm_of(s) for s in strata]
    assert set(assignments) <= set(ARMS)
    assert len(assignments) == len(strata), "arm_of must be a function, not a fan-out"
    assert arm_of("mlh1_intact_mmrd") != arm_of("mmr_proficient"), (
        "the pre-registered control arm must not be folded into the secondary "
        "arm; it is reported separately because it carries no verdict"
    )


def test_both_consumers_size_the_secondary_arm_the_same_way():
    """The two files, run over the same cohort table, must agree on every arm.

    THE INPUT THAT WOULD HAVE CAUGHT IT. This is the comparison nobody made:
    the reading job's arm sizes against the sizing job's, on the committed
    data, before either number reached a document.
    """
    import glob

    from src.reference.mlh1_arms import arm_of

    deltas = sorted(RESULTS_DIR.glob("*/icbi_coexpression.parquet"))
    cohorts = sorted(glob.glob(str(RESULTS_DIR / "*" / "cohort_table.parquet")))
    if not deltas or not cohorts:
        pytest.skip("needs the committed Pelka deltas and cohort table")

    scored = pd.read_parquet(deltas[-1])
    scored = scored[scored["study_id"] == "Pelka_2021_Cell"].copy()
    scored["short_id"] = scored["patient_id"].astype(str).str.split(".").str[-1]
    table = pd.read_parquet(cohorts[-1])
    joined = (scored.drop_duplicates("short_id")
              .merge(table[["patient_id", "mlh1_stratum"]],
                     left_on="short_id", right_on="patient_id",
                     suffixes=("", "_c")))

    by_arm = joined["mlh1_stratum"].map(arm_of).value_counts()
    # The superseded hand-rolled split, kept as the thing that must NOT be used.
    hand_rolled = int((joined["mlh1_stratum"] != "mlh1_methylated").sum())
    assert hand_rolled != by_arm.get("mlh1_unmethylated", 0), (
        "if these ever coincide the fixture has stopped exercising the defect"
    )
    assert sum(by_arm) == len(joined), "the arms must exhaust the scored cohort"
