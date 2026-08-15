"""W3 — Bulk & clinical. Owns TCGA-COAD/READ. Fully independent for weeks 1-5.

Produces: harmonised bulk matrix on the SHARED gene index (the same one W1's
S matrices use — integration is then a join, not a negotiation), tumour purity,
and a curated clinical table.

Survival endpoints are DSS and PFI from TCGA-CDR (Liu et al. 2018). OS is
secondary — COAD OS is heavily contaminated by non-cancer death
(CLAUDE.md invariant 9).

The week-2 premise check is the highest-leverage thing in this workstream: if
GUCA2A loss in TCGA-COAD is continuous rather than bimodal, the two-type
classification dissolves into a regression and the project changes shape.
Report it to the whole team in week 2, not week 5.
"""

__all__: list[str] = []
