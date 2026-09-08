from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix

import src.reference.jobs.becker_lesion_wnt_gate as gate
from src.reference.jobs.becker_lesion_wnt_gate import (
    EXPECTED_PAIRED_DONORS,
    MIN_DETECTED_NUCLEI,
    WntDetectionGateError,
    gate_verdict,
    read_paired_polyp_blocks,
    signature_index,
    summarise_blocks,
    validate_specification,
)
from src.reference.wnt_score import SIGNATURE


def _blocks(*, detected_per_gene: int = 40) -> list[dict[str, object]]:
    symbols = list(SIGNATURE) + ["OTHER"]
    blocks = []
    for donor in EXPECTED_PAIRED_DONORS:
        counts = np.zeros((100, len(symbols)), dtype=int)
        counts[:detected_per_gene, :len(SIGNATURE)] = 1
        blocks.append({
            "donor": donor,
            "sample_id": f"{donor}-C-001",
            "gsm": f"GSM-{donor}",
            "counts": csr_matrix(counts),
            "symbols": symbols,
        })
    return blocks


def test_all_signature_genes_must_be_exact_feature_matches():
    with pytest.raises(WntDetectionGateError, match="missing symbols"):
        signature_index(SIGNATURE[:-1])
    with pytest.raises(WntDetectionGateError, match="duplicated symbols"):
        signature_index(list(SIGNATURE) + [SIGNATURE[0]])


def test_the_pre_read_provenance_guard_accepts_the_declared_non_circular_spec():
    """Exercise the exact call the CLI makes before it checks paths or reads data."""
    assert validate_specification() == ()


def test_gate_passes_only_when_every_donor_gene_row_clears_both_floors():
    by_donor, summary = summarise_blocks(_blocks())
    assert len(by_donor) == len(EXPECTED_PAIRED_DONORS) * len(SIGNATURE)
    assert (by_donor["n_detected"] >= MIN_DETECTED_NUCLEI).all()
    assert by_donor["passes"].all()
    assert summary["passes_all_paired_donors"].all()
    assert gate_verdict(by_donor)["verdict"] == "PASSES DETECTION GATE"


def test_one_sparse_donor_gene_is_no_substrate_not_a_reduced_signature():
    by_donor, _ = summarise_blocks(_blocks(detected_per_gene=1))
    outcome = gate_verdict(by_donor)
    assert not by_donor["passes"].any()
    assert outcome["verdict"] == "NO SUBSTRATE AT THIS snRNA RESOLUTION"
    assert "reduced signature" in outcome["detail"]


def test_changed_feature_order_refuses_to_combine_nuclei():
    blocks = _blocks()
    blocks[-1]["symbols"] = list(reversed(SIGNATURE)) + ["OTHER"]
    with pytest.raises(WntDetectionGateError, match="different feature order"):
        summarise_blocks(blocks)


def test_incomplete_selected_triplet_is_named_not_silently_dropped(monkeypatch, tmp_path):
    metadata = pd.DataFrame({
        "gsm": ["GSM-A001", "GSM-CRC"],
        "donor": ["A001", "CRC1"],
        "arm": ["tumour", None],
        "sample_id": ["A001-C-001", "CRC1_8810"],
    })
    files = pd.DataFrame({
        "gsm": ["GSM-A001", "GSM-CRC"], "complete": [False, True],
        "barcodes": ["barcodes", "barcodes"],
        "features": ["features", "features"],
        "matrix": [np.nan, "matrix"],
        # The CRC name is intentionally inconsistent: it is outside the gate's
        # selected cohort and must not obscure the selected triplet failure.
        "sample_id": ["A001-C-001", "wrong-ignored-sample"],
    })
    monkeypatch.setattr("src.reference.becker_io.read_series_matrix", lambda _: metadata)
    monkeypatch.setattr("src.reference.becker_io.paired_donors", lambda _: EXPECTED_PAIRED_DONORS)
    monkeypatch.setattr("src.reference.becker_io.sample_files", lambda _: files)
    with pytest.raises(WntDetectionGateError, match="No sample was silently dropped") as raised:
        read_paired_polyp_blocks(tmp_path / "raw.tar", tmp_path / "series.txt.gz")
    assert raised.value.failure_kind == "incomplete_triplet"
    outcome = gate.failure_outcome(raised.value)
    assert outcome["verdict"] == "NO SUBSTRATE — incomplete assay triplet"


def test_identifier_failure_writes_durable_no_substrate_artifact(monkeypatch, tmp_path):
    tar = tmp_path / "raw.tar"
    series = tmp_path / "series.txt.gz"
    tar.touch()
    series.touch()
    written: dict[str, object] = {}

    def refuse(*_args, **_kwargs):
        raise WntDetectionGateError("missing symbols ['TCF7']", failure_kind="identifier_assay")

    def record(frame, name, **kwargs):
        written.update({"frame": frame, "name": name, "meta": kwargs["extra_meta"]})
        return tmp_path / f"{name}.parquet"

    monkeypatch.setattr(gate, "read_paired_polyp_blocks", refuse)
    monkeypatch.setattr(gate, "write_versioned_table", record)
    assert gate.main(["--tar", str(tar), "--series-matrix", str(series)]) == 5
    assert written["name"] == "becker_lesion_wnt_detection_failure"
    assert written["frame"].iloc[0]["verdict"] == "NO SUBSTRATE — identifier/assay failure"
    assert written["meta"]["verdict"]["failure_kind"] == "identifier_assay"
