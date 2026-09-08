from __future__ import annotations

import numpy as np
import pytest
from scipy.sparse import csr_matrix

from src.reference.jobs.becker_lesion_wnt_gate import (
    EXPECTED_PAIRED_DONORS,
    MIN_DETECTED_NUCLEI,
    WntDetectionGateError,
    gate_verdict,
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
