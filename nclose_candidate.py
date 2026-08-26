"""Small compatibility layer for stage-01 NClose candidates.

The legacy graph representation stores pairs in query/path order under a
contig-name key.  That order is meaningful and must not be replaced by the
sorted node-index key used by tracking and reporting code.

This module deliberately does not own node rows or persisted artifact schemas.
It provides an import-safe candidate view, lossless legacy adapters, and a
minimal order-preserving filter helper.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Hashable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field


CTG_RPTCHR = 16
CTG_RPTCASE = 17
CTG_CENSAT = 18

CandidateIdentity = tuple[str, tuple[int, int]]
RejectReason = Callable[["NCloseCandidate"], str | None]


@dataclass(frozen=True)
class NCloseCandidate:
    """One graph NClose with path order kept separate from canonical keys."""

    contig_name: str
    path_pair: tuple[int, int]
    origin: str = "legacy_unknown"
    synthetic: bool = False
    signatures: Mapping[str, Hashable] = field(default_factory=dict)
    provenance: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        pair = tuple(int(node_idx) for node_idx in self.path_pair)
        if len(pair) != 2:
            raise ValueError(f"Invalid NClose path pair: {self.path_pair!r}")
        object.__setattr__(self, "contig_name", str(self.contig_name))
        object.__setattr__(self, "path_pair", pair)
        object.__setattr__(self, "origin", str(self.origin))
        object.__setattr__(self, "synthetic", bool(self.synthetic))
        object.__setattr__(self, "signatures", dict(self.signatures))
        object.__setattr__(self, "provenance", dict(self.provenance))

    @property
    def identity(self) -> CandidateIdentity:
        """Return the exact legacy owner and path-ordered pair."""

        return self.contig_name, self.path_pair

    @property
    def event_key(self) -> tuple[int, int]:
        """Return the sorted node-index key used by tracking artifacts."""

        return tuple(sorted(self.path_pair))


@dataclass(frozen=True)
class NCloseRejection:
    """The first filter stage that removed a candidate."""

    candidate: NCloseCandidate
    stage: str
    reason: str

    def __post_init__(self) -> None:
        if not str(self.stage):
            raise ValueError("NClose rejection stage cannot be empty")
        if not str(self.reason):
            raise ValueError("NClose rejection reason cannot be empty")
        object.__setattr__(self, "stage", str(self.stage))
        object.__setattr__(self, "reason", str(self.reason))


def iter_legacy_candidates(
    nclose_nodes: Mapping,
    origin_by_identity: Mapping[CandidateIdentity, str] | None = None,
):
    """Yield candidates without changing legacy mapping, list, or pair order."""

    origins = origin_by_identity or {}
    for contig_name_value, pair_list in nclose_nodes.items():
        contig_name = str(contig_name_value)
        for pair_value in pair_list:
            pair = tuple(int(node_idx) for node_idx in pair_value)
            identity = (contig_name, pair)
            yield NCloseCandidate(
                contig_name=contig_name,
                path_pair=pair,
                origin=origins.get(identity, "legacy_unknown"),
            )


def candidates_to_legacy(
    candidates: Iterable[NCloseCandidate],
) -> defaultdict[str, list[tuple[int, int]]]:
    """Return the current ``defaultdict(list)`` graph representation."""

    nclose_nodes = defaultdict(list)
    for candidate in candidates:
        nclose_nodes[candidate.contig_name].append(candidate.path_pair)
    return nclose_nodes


def endpoint_metadata(
    candidate: NCloseCandidate,
    contig_data: Sequence[Sequence],
) -> tuple[dict, dict]:
    """Read exact repeat/CENSAT labels in candidate path order.

    Labels are read live from ``contig_data`` so synthetic-node annotation
    cannot make a cached candidate snapshot stale.
    """

    endpoints = []
    for node_idx in candidate.path_pair:
        node = contig_data[node_idx]
        if len(node) <= CTG_CENSAT:
            raise ValueError(
                f"NClose node {node_idx} has {len(node)} columns; "
                f"expected at least {CTG_CENSAT + 1}"
            )
        endpoints.append(
            {
                "node_idx": node_idx,
                "repeat_chrom": node[CTG_RPTCHR],
                "repeat_case": node[CTG_RPTCASE],
                "censat": node[CTG_CENSAT],
            }
        )
    return tuple(endpoints)


def apply_nclose_filter(
    candidates: Iterable[NCloseCandidate],
    stage: str,
    reject_reason: RejectReason,
) -> tuple[list[NCloseCandidate], list[NCloseRejection]]:
    """Apply one stable filter; ``None`` means keep, a string means reject."""

    kept = []
    rejected = []
    for candidate in candidates:
        reason = reject_reason(candidate)
        if reason is None:
            kept.append(candidate)
            continue
        if not isinstance(reason, str) or not reason:
            raise ValueError(
                "NClose reject_reason must return None or a non-empty string"
            )
        rejected.append(
            NCloseRejection(candidate=candidate, stage=stage, reason=reason)
        )
    return kept, rejected
