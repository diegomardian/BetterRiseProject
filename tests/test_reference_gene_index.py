"""The shared gene index. W1 emits it, W3 conforms to it.

The join between W1's S matrices and W3's bulk is the project's integration
point, and §3.4 wants it to be "a join, not a negotiation". These tests pin the
things that would make it a negotiation: the key, the collision policy, and
whether the frozen panel actually resolves.
"""

from __future__ import annotations

import pandas as pd
import pytest

from src.common.panel import axis_genes, panel_genes
from src.reference.gene_index import (
    GENOME_NOTE,
    GeneIndexError,
    build_gene_index,
    check_panel_coverage,
    read_gene_index,
    write_gene_index,
)


def var_table(extra: list[tuple[str, str]] | None = None) -> pd.DataFrame:
    """A feature table shaped like GSE178341's, including its dedup suffixes."""
    rows = [
        ("ENSG00000000001.5_1", "GUCA2A"),
        ("ENSG00000000002.2_2", "MLH1"),
        ("ENSG00000000003.6", "LGR5"),
        ("ENSG00000000004.1_3", "MUC2"),
        ("ENSG00000000005.4", "EPCAM"),
    ] + (extra or [])
    return pd.DataFrame(
        {
            "feature_id": [r[0] for r in rows],
            "gene_symbol": [r[1] for r in rows],
            "genome": ["GRCh37_liftover_v28"] * len(rows),
            "ensembl_id": [r[0].split(".")[0] for r in rows],
        }
    )


class TestBuild:
    def test_one_row_per_unversioned_ensembl_id(self):
        index, mapping, report = build_gene_index(var_table())
        assert len(index) == 5
        assert all(i.startswith("ENSG") and "." not in i and "_" not in i for i in index)
        assert report["n_genes_out"] == 5

    def test_duplicate_features_collapse_and_are_counted(self):
        """CellRanger's _N suffix means the same gene appears more than once.
        Stripping it collapses them; the collision must be counted, not hidden —
        whoever builds a matrix has to decide whether to sum or drop."""
        var = var_table(extra=[("ENSG00000000001.5_2", "GUCA2A")])
        index, mapping, report = build_gene_index(var)
        assert len(index) == 5                      # not 6
        assert report["n_collapsed"] == 1
        assert report["n_duplicated_ids"] == 1
        assert int(mapping.set_index("ensembl_id").loc["ENSG00000000001", "n_features"]) == 2

    def test_order_follows_the_deposit_and_is_reproducible(self):
        var = var_table()
        first, _, _ = build_gene_index(var)
        second, _, _ = build_gene_index(var)
        assert first == second
        assert first[0] == "ENSG00000000001"

    def test_genome_build_is_reported(self):
        _, _, report = build_gene_index(var_table())
        assert report["genomes"] == ["GRCh37_liftover_v28"]

    def test_missing_column_raises(self):
        with pytest.raises(GeneIndexError, match="missing column"):
            build_gene_index(pd.DataFrame({"ensembl_id": ["ENSG1"]}))

    def test_empty_input_raises(self):
        with pytest.raises(GeneIndexError, match="empty"):
            build_gene_index(var_table().iloc[:0])


class TestPanelCoverage:
    def test_reports_which_frozen_genes_resolve(self):
        _, mapping, _ = build_gene_index(var_table())
        coverage = check_panel_coverage(mapping)
        found = coverage.set_index("gene_symbol")["found"]
        assert bool(found["GUCA2A"])
        assert bool(found["LGR5"])

    def test_absent_panel_genes_are_flagged_not_dropped(self):
        """A panel gene absent from the index cannot be tested at all. It must
        appear in the report as unfound rather than simply be missing from it."""
        _, mapping, _ = build_gene_index(var_table())
        coverage = check_panel_coverage(mapping)
        assert set(coverage[coverage["source"] == "panel"]["gene_symbol"]) == set(
            panel_genes()
        )
        assert (~coverage["found"]).any()          # most of the panel is absent here

    def test_both_axes_are_checked(self):
        _, mapping, _ = build_gene_index(var_table())
        coverage = check_panel_coverage(mapping)
        assert set(coverage["source"]) == {
            "panel", "axis:stem_pole", "axis:opposite_lineage",
        }
        stem = coverage[coverage["source"] == "axis:stem_pole"]
        assert set(stem["gene_symbol"]) == set(axis_genes("stem_pole"))


class TestWriteAndRead:
    def test_round_trip(self, tmp_path):
        index, mapping, _ = build_gene_index(var_table())
        write_gene_index(index, mapping, version="1.0.0", config_dir=tmp_path)
        assert read_gene_index("1.0.0", config_dir=tmp_path) == index

    def test_the_genome_note_travels_with_the_map(self, tmp_path):
        """W3 must not have to remember that this is hg19."""
        index, mapping, _ = build_gene_index(var_table())
        _, map_path = write_gene_index(index, mapping, version="1.0.0", config_dir=tmp_path)
        text = map_path.read_text()
        assert "GRCh37" in text and "GRCh38" in text
        assert GENOME_NOTE.split(".")[0] in text

    def test_overwriting_a_version_is_refused(self, tmp_path):
        """Results reference the index by version. Editing one in place makes
        every earlier result unreproducible without saying so."""
        index, mapping, _ = build_gene_index(var_table())
        write_gene_index(index, mapping, version="1.0.0", config_dir=tmp_path)
        with pytest.raises(GeneIndexError, match="Bump the version"):
            write_gene_index(index, mapping, version="1.0.0", config_dir=tmp_path)

    def test_reading_a_missing_version_raises(self, tmp_path):
        with pytest.raises(GeneIndexError, match="does not exist"):
            read_gene_index("9.9.9", config_dir=tmp_path)

    def test_the_map_is_tab_separated_with_a_stable_header(self, tmp_path):
        index, mapping, _ = build_gene_index(var_table())
        _, map_path = write_gene_index(index, mapping, version="1.0.0", config_dir=tmp_path)
        body = [ln for ln in map_path.read_text().splitlines() if not ln.startswith("#")]
        assert body[0].split("\t")[:2] == ["ensembl_id", "gene_symbol"]


def test_the_index_is_what_build_signature_consumes(tmp_path):
    """build_signature() asserts target genes are absent from the index it is
    handed. Confirm the emitted form is directly usable."""
    from src.reference.signature import (
        LeakageError,
        LeakageGuardError,
        assert_no_target_leakage,
    )

    index, mapping, _ = build_gene_index(var_table())
    write_gene_index(index, mapping, version="1.0.0", config_dir=tmp_path)
    loaded = read_gene_index("1.0.0", config_dir=tmp_path)

    # This test previously asserted the guard PASSED here, reasoning that
    # "Ensembl-keyed, so symbol-named panel genes cannot collide with it."
    # That is the bug (issue #35): the guard could not collide, so it could not
    # fire, and GUCA2A shipped in all four 0.1.0-pilot S matrices while every
    # call site reported success. A check that cannot fail must refuse.
    with pytest.raises(LeakageGuardError, match="cannot check invariant 2"):
        assert_no_target_leakage(loaded, panel_genes(), context="the shared gene index")

    # Given the translation it can do its job — and on an index built from a
    # table that carries a panel gene, it must fire.
    alias = dict(
        zip(
            mapping["gene_symbol"].astype(str),
            mapping["ensembl_id"].astype(str),
            strict=True,
        )
    )
    present = [g for g in panel_genes() if g in alias and alias[g] in set(loaded)]
    if present:
        with pytest.raises(LeakageError, match="invariant 2"):
            assert_no_target_leakage(
                loaded, present, context="the shared gene index", alias_map=alias
            )


class TestIntersection:
    """Open decision #2, corrected: 1.0.0 IS the intersection."""

    @staticmethod
    def _panel_ids():
        """Every frozen panel gene, so the coverage guard is satisfied."""
        from src.common.panel import panel_genes
        panel = panel_genes()
        return [f"ENSGP{i:04d}" for i in range(len(panel))], panel

    @classmethod
    def _mapping(cls, ids, symbols=None):
        """A mapping carrying the whole panel PLUS the ids under test."""
        import pandas as pd
        pids, panel = cls._panel_ids()
        return pd.DataFrame({
            "ensembl_id": list(ids) + pids,
            "gene_symbol": (symbols or [f"G{i}" for i in range(len(ids))]) + panel,
        })

    def test_keeps_only_shared_ids_and_sorts(self):
        from src.reference.gene_index import intersect_gene_index
        pids, _ = self._panel_ids()
        mapping = self._mapping(["ENSG3", "ENSG1", "ENSG2"])
        index, out, report = intersect_gene_index(
            mapping, ["ENSG2", "ENSG3", "ENSG9"] + pids
        )
        assert index == sorted(index), "must be sorted ascending"
        assert [g for g in index if g.startswith("ENSG")] == sorted(
            ["ENSG2", "ENSG3"] + pids
        )
        assert list(out["ensembl_id"]) == index
        assert report["n_intersection"] == 2 + len(pids)
        assert report["n_ours_only"] == 1
        assert report["n_shared_only"] == 1

    def test_byte_identical_regardless_of_input_order(self):
        """What makes 'who emits it' stop mattering."""
        from src.reference.gene_index import intersect_gene_index
        pids, _ = self._panel_ids()
        a, _, _ = intersect_gene_index(
            self._mapping(["ENSG3", "ENSG1", "ENSG2"]),
            ["ENSG1", "ENSG2", "ENSG3"] + pids,
        )
        b, _, _ = intersect_gene_index(
            self._mapping(["ENSG1", "ENSG2", "ENSG3"]),
            list(reversed(["ENSG3", "ENSG2", "ENSG1"] + pids)),
        )
        assert a == b

    def test_refuses_when_a_panel_gene_is_lost(self):
        from src.reference.gene_index import GeneIndexError, intersect_gene_index
        pids, _ = self._panel_ids()
        mapping = self._mapping([])
        with pytest.raises(GeneIndexError, match="do not survive"):
            intersect_gene_index(mapping, pids[1:])   # drop one panel gene

    def test_disjoint_indices_raise(self):
        from src.reference.gene_index import GeneIndexError, intersect_gene_index
        with pytest.raises(GeneIndexError, match="share no identifier"):
            intersect_gene_index(self._mapping(["ENSG1"]), ["ENSG2"])

    def test_versioned_shared_index_is_rejected(self, tmp_path):
        from src.reference.gene_index import GeneIndexError, read_shared_index
        path = tmp_path / "v.txt"
        path.write_text("ENSG00000141510.16\nENSG00000000003\n")
        with pytest.raises(GeneIndexError, match="versioned"):
            read_shared_index(path)

    def test_reads_a_clean_index(self, tmp_path):
        from src.reference.gene_index import read_shared_index
        path = tmp_path / "i.txt"
        path.write_text("# comment\nENSG1\n\nENSG2\n")
        assert read_shared_index(path) == ["ENSG1", "ENSG2"]

    def test_w3_index_carries_every_frozen_panel_gene(self):
        """Independent check of W3's 23/23 claim, bulk side.

        The single-cell side needs the deposit, which lives on the cluster.
        This half is checkable from the committed file.
        """
        from pathlib import Path

        import pandas as pd

        from src.common.panel import panel_genes
        f = Path("config/gene_index/gene_index_0.9.0.map.tsv")
        if not f.exists():
            pytest.skip("0.9.0 not present")
        symbols = set(pd.read_csv(f, sep="\t", comment="#")["gene_symbol"].astype(str))
        assert not [g for g in panel_genes() if g not in symbols]
