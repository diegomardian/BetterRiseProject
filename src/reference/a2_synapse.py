"""Read-only metadata inventory for the access-gated Chen MxIF products.

This module intentionally has no module-level ``synapseclient`` import.  The
public test suite should not need an account, and the only code path that needs
the client is the authenticated command-line job.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from typing import Any

import pandas as pd

# Named in Chen et al.'s data-availability statement.  They are roots to
# inventory, not permissions to download their descendants.
CHEN_SYNAPSE_ROOTS = ("syn23564801", "syn23630431", "syn23520239")

_CANDIDATE = re.compile(
    r"(?:mand[ou0]|seg(?:ment|mentation)?|mask|cell|intens|feature|"
    r"quantif|table|csv|tsv|parquet|h5ad)",
    flags=re.IGNORECASE,
)


class A2SynapseError(RuntimeError):
    """The authenticated, metadata-only inventory cannot be completed."""


def _as_text(value: Any) -> str:
    """A stable searchable rendering of scalar/list annotation values."""
    if value is None:
        return ""
    if isinstance(value, (list, tuple, set)):
        return " | ".join(map(str, value))
    return str(value)


def candidate_product(name: str, annotations: dict[str, Any]) -> bool:
    """Flag a possible cell-derived product for review, never as a conclusion."""
    haystack = " ".join([name, *map(_as_text, annotations.values())])
    return bool(_CANDIDATE.search(haystack))


def walk_metadata(syn: Any, roots: Iterable[str] = CHEN_SYNAPSE_ROOTS) -> pd.DataFrame:
    """Recursively list entity metadata, with no file-content request.

    ``Synapse.get(..., downloadFile=False)`` retrieves entity properties only;
    ``getChildren`` returns child metadata.  An inaccessible child is an error,
    not an absent product: a partial directory is exactly the failure mode this
    inventory is designed to expose.
    """
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    pending: list[tuple[str, str | None]] = [(str(root), None) for root in roots]
    while pending:
        entity_id, root_id = pending.pop()
        if entity_id in seen:
            continue
        seen.add(entity_id)
        try:
            entity = syn.get(entity_id, downloadFile=False)
            annotations = dict(syn.get_annotations(entity_id) or {})
        except Exception as exc:  # pragma: no cover - exercised through the CLI
            raise A2SynapseError(
                f"cannot read metadata for {entity_id}; do not treat it as absent: {exc}"
            ) from exc

        properties = getattr(entity, "properties", {}) or {}
        name = str(properties.get("name", getattr(entity, "name", entity_id)))
        entity_type = str(properties.get("concreteType", properties.get("type", "unknown")))
        current_root = root_id or entity_id
        rows.append({
            "root_entity_id": current_root,
            "entity_id": entity_id,
            "parent_id": properties.get("parentId"),
            "entity_name": name,
            "entity_type": entity_type,
            "annotations_json": {
                str(key): _as_text(value) for key, value in sorted(annotations.items())
            },
            "candidate_cell_product": candidate_product(name, annotations),
        })
        try:
            children = list(syn.getChildren(parent=entity_id))
        except Exception as exc:  # File entities are leaves; other failures are not benign.
            if "not a container" in str(exc).lower() or "cannot have children" in str(exc).lower():
                children = []
            else:  # pragma: no cover - depends on the remote service
                raise A2SynapseError(
                    f"cannot list children for {entity_id}; inventory is partial: {exc}"
                ) from exc
        pending.extend((str(child["id"]), current_root) for child in children)

    frame = pd.DataFrame(rows)
    if frame.empty:
        raise A2SynapseError("Synapse returned no metadata rows")
    if frame["entity_id"].duplicated().any():
        raise A2SynapseError("an entity was emitted more than once")
    return frame.sort_values(["root_entity_id", "entity_name", "entity_id"], ignore_index=True)


def login_synapse() -> Any:
    """Return an authenticated official client without accepting credentials in code."""
    try:
        import synapseclient
    except ImportError as exc:  # pragma: no cover - depends on the optional extra
        raise A2SynapseError(
            "synapseclient is not installed; install the A2 extra with "
            "`pip install -e '.[a2]'` in the W1 environment"
        ) from exc
    syn = synapseclient.Synapse()
    syn.login(silent=True)
    return syn
