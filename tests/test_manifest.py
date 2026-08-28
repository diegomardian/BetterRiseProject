"""`data/manifest.csv` is the only thing that travels — `data/` is gitignored in full.

Issue #43: five W3 rows carried an absolute Windows path, which is unresolvable
anywhere but the machine that wrote them. One of them was TCGA-CDR, which
invariant 9's DSS and PFI endpoints come from. The cause looks like
``join(RAW_DIR, already_absolute_path)``, which is a mistake nobody would make
twice on purpose and exactly the kind a test catches for free.

WHAT THIS DOES NOT CHECK, AND WHY
---------------------------------
Issue #43 lists two other problems that are **not** asserted here:

1. **GSE178341 and the ICBI atlas have no rows at all** (#43 item 1). Asserting
   "the primary dataset is recorded" would fail on `main` today, and the fix is
   W1's — they are the only ones who can compute the sha256 of a file they have.
2. **Some rows are shifted by one column** (#43 item 3) — the four Lee rows and
   `gene_order_hg19.txt` carry note text in `workstream` and a workstream code in
   `downloaded_by`. Those are W4's and W1's rows respectively.

Asserting either would turn `main` red over someone else's data, which is not
what a test added in a W3 PR should do. They are named here so this file is not
mistaken for a full validation of the manifest.
"""

from __future__ import annotations

import re

import pandas as pd
import pytest

from src.common.paths import MANIFEST_PATH

COLUMNS = (
    "path",
    "sha256",
    "bytes",
    "source_url",
    "accession",
    "downloaded_on",
    "downloaded_by",
    "workstream",
    "notes",
)


@pytest.fixture(scope="module")
def manifest() -> pd.DataFrame:
    return pd.read_csv(MANIFEST_PATH, dtype=str, keep_default_na=False)


def test_the_manifest_exists_and_is_not_empty(manifest):
    assert len(manifest) > 0


def test_the_header_is_the_documented_one(manifest):
    assert tuple(manifest.columns) == COLUMNS


def test_every_path_is_repo_relative(manifest):
    """The failure in #43. A drive letter resolves on exactly one machine."""
    absolute = manifest.loc[manifest["path"].str.contains(":", regex=False), "path"]
    assert absolute.empty, (
        f"{len(absolute)} row(s) carry an absolute path: {list(absolute)[:3]}. "
        f"Record the repo-relative path — data/ is gitignored, so the manifest "
        f"is the only thing another machine can act on."
    )


def test_no_path_uses_backslashes(manifest):
    """Three of the four workstreams are not on Windows."""
    windows = manifest.loc[manifest["path"].str.contains("\\\\", regex=True), "path"]
    assert windows.empty, f"{len(windows)} row(s) use backslashes: {list(windows)[:3]}"


def test_every_path_is_under_data(manifest):
    stray = manifest.loc[~manifest["path"].str.startswith("data/"), "path"]
    assert stray.empty, f"{len(stray)} row(s) are not under data/: {list(stray)[:3]}"


def test_no_file_is_recorded_twice(manifest):
    """Two rows for one file means two checksums that can disagree later."""
    duplicated = manifest.loc[manifest["path"].duplicated(), "path"]
    assert duplicated.empty, f"duplicate rows for: {list(duplicated)[:3]}"


def test_every_sha256_is_a_sha256(manifest):
    pattern = re.compile(r"[0-9a-f]{64}")
    bad = manifest.loc[
        ~manifest["sha256"].map(lambda s: bool(pattern.fullmatch(s))), "path"
    ]
    assert bad.empty, f"{len(bad)} row(s) have a malformed sha256: {list(bad)[:3]}"


def test_every_size_is_a_positive_integer(manifest):
    bad = manifest.loc[
        ~manifest["bytes"].map(lambda s: s.isdigit() and int(s) > 0), "path"
    ]
    assert bad.empty, f"{len(bad)} row(s) have a bad byte count: {list(bad)[:3]}"


def test_every_row_records_where_it_came_from(manifest):
    """`data/README.md`: a download nobody recorded is a result nobody can
    reproduce. A row with no source is a row that cannot be re-fetched."""
    bare = manifest.loc[manifest["source_url"].str.strip() == "", "path"]
    assert bare.empty, f"{len(bare)} row(s) have no source_url: {list(bare)[:3]}"
