"""The MLH1 reading's arms, defined ONCE.

WHY THIS FILE EXISTS. It did not, and the arms were therefore defined twice:
``mlh1_positive_control.arm_of`` sent ``mlh1_intact_mmrd`` to its own arm, so
the secondary arm was **15** patients, while ``interval_calibration`` sized that
same arm as ``mlh1_stratum != "mlh1_methylated"``, which is **19**. The
pre-registration's §3 table and its two-sample power figures inherited the 19.

Neither number was wrong on its own terms. They were answers to two different
questions that had been given one name, in two files, with nothing comparing
them -- and the disagreement surfaced only when the cluster run printed
``mlh1_unmethylated  15`` against a document that said 19.

That is this repository's signature defect in its plainest form: not a check
that could not fail, but a quantity with two definitions and no check at all.
The remedy is the same one the panel and the labelling axes already use -- one
definition, imported by everyone, and a test that the importers agree.

The constants live here rather than in either consumer because
``mlh1_positive_control`` imports ``interval_calibration`` for the interval, so
the dependency cannot run the other way without a cycle.
"""

from __future__ import annotations

#: The arm the reading is about: promoters known methylated from an assay.
PRIMARY_STRATUM = "mlh1_methylated"

#: Everything not methylated and not the pre-registered control arm. CONFOUNDED
#: with MMR status -- it is mostly MMR-proficient patients -- and reported as
#: such, never as the mechanistic control.
SECONDARY_STRATUM = "mlh1_unmethylated"

#: The control the original difference-in-differences was built on. Reported,
#: never read: it is four patients, where the project's usual interval excludes
#: zero 18.8% of the time under a true null.
UNDERPOWERED_STRATUM = "mlh1_intact_mmrd"

ARMS: tuple[str, ...] = (PRIMARY_STRATUM, SECONDARY_STRATUM, UNDERPOWERED_STRATUM)


def arm_of(stratum: str) -> str:
    """Which reported arm a pre-registered stratum belongs to.

    The three arms are DISJOINT and exhaust the cohort. In particular
    ``mlh1_intact_mmrd`` is **not** part of the secondary arm: it is broken out
    because the original prereg named it a mechanistic control, and counting it
    in both places would report four patients twice under two standings.
    """
    if stratum == PRIMARY_STRATUM:
        return PRIMARY_STRATUM
    if stratum == UNDERPOWERED_STRATUM:
        return UNDERPOWERED_STRATUM
    return SECONDARY_STRATUM
