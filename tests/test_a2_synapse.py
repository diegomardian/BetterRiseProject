"""The A2 access-gated inventory must remain metadata-only and complete."""

from __future__ import annotations

import pytest

from src.reference.a2_synapse import A2SynapseError, candidate_product, walk_metadata


class _Entity:
    def __init__(self, entity_id, name, kind="folder", parent=None):
        self.properties = {
            "id": entity_id,
            "name": name,
            "concreteType": kind,
            "parentId": parent,
        }


class _Synapse:
    def __init__(self):
        self.entities = {
            "syn_root": _Entity("syn_root", "Chen MxIF"),
            "syn_cells": _Entity("syn_cells", "MANDO_cell_intensities.csv", "file", "syn_root"),
            "syn_raw": _Entity("syn_raw", "raw_image.ome.tif", "file", "syn_root"),
        }
        self.annotations = {
            "syn_cells": {"dataType": ["single-cell quantification"]},
            "syn_raw": {"dataType": ["multiplex image"]},
        }
        self.children = {
            "syn_root": [{"id": "syn_cells"}, {"id": "syn_raw"}],
        }
        self.download_flags = []

    def get(self, entity_id, downloadFile=True):
        self.download_flags.append(downloadFile)
        return self.entities[entity_id]

    def get_annotations(self, entity_id):
        return self.annotations.get(entity_id, {})

    def getChildren(self, parent):
        if parent not in self.children:
            raise RuntimeError("entity is not a container")
        return self.children[parent]


def test_walk_is_recursive_and_never_requests_file_content():
    syn = _Synapse()
    got = walk_metadata(syn, roots=("syn_root",)).set_index("entity_id")
    assert set(got.index) == {"syn_root", "syn_cells", "syn_raw"}
    assert syn.download_flags and not any(syn.download_flags)
    assert bool(got.loc["syn_cells", "candidate_cell_product"])
    assert not bool(got.loc["syn_raw", "candidate_cell_product"])
    assert got.loc["syn_cells", "root_entity_id"] == "syn_root"


def test_candidate_screen_uses_annotations_as_well_as_the_filename():
    assert candidate_product("opaque.bin", {"assay": ["cell segmentation"]})
    assert not candidate_product("opaque.bin", {"assay": ["raw imaging"]})


def test_inaccessible_child_refuses_a_partial_inventory():
    class Broken(_Synapse):
        def get(self, entity_id, downloadFile=True):
            if entity_id == "syn_cells":
                raise RuntimeError("403 forbidden")
            return super().get(entity_id, downloadFile=downloadFile)

    with pytest.raises(A2SynapseError, match="do not treat it as absent"):
        walk_metadata(Broken(), roots=("syn_root",))
