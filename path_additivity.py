"""Helpers for resolving compositional NClose path conflicts.

For three matrix columns A, B, and AB, the relationship is additive only when
the complete NClose usage multiset of AB is the element-wise sum of A and B.
Using a multiset (``Counter``), rather than a plain set, preserves repeated
traversals of the same breakend.  The chained A/BC/AB pattern is likewise
exact: AB = A + B and BC = B + C for non-empty, disjoint A and BC.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from math import isfinite
from typing import Hashable, Mapping, Sequence


EventKey = Hashable
UsageSignature = frozenset[tuple[EventKey, int]]


@dataclass(frozen=True)
class AdditivePathTriple:
    """One path-column triple satisfying usage(AB) = usage(A) + usage(B)."""

    path_a: int
    path_b: int
    path_ab: int
    weight_a: float
    weight_b: float
    weight_ab: float

    @property
    def total_weight(self) -> float:
        return self.weight_a + self.weight_b + self.weight_ab


@dataclass(frozen=True)
class ChainedPathTriple:
    """One strict A/BC/AB chain satisfying AB=A+B and BC=B+C."""

    path_a: int
    path_bc: int
    path_ab: int
    weight_a: float
    weight_bc: float
    weight_ab: float

    @property
    def total_weight(self) -> float:
        return self.weight_a + self.weight_bc + self.weight_ab


@dataclass(frozen=True)
class PathConflict:
    """Two removable side paths and the composite path they conflict with.

    ``path_b`` is the B path for ``A_B_AB`` and the BC path for
    ``A_BC_AB``.  The role fields make that distinction explicit in logs and
    audit output while allowing both conflict types to share one removal
    policy and one globally weight-sorted queue.
    """

    conflict_type: str
    path_a: int
    path_b: int
    path_ab: int
    role_a: str
    role_b: str
    role_ab: str
    weight_a: float
    weight_b: float
    weight_ab: float

    @property
    def total_weight(self) -> float:
        return self.weight_a + self.weight_b + self.weight_ab


@dataclass(frozen=True)
class RemovalChoice:
    """The path column(s) selected for permanent removal from one triple."""

    removed_columns: tuple[int, ...]
    reason: str


def usage_signature(usage: Mapping[EventKey, int]) -> UsageSignature:
    """Return a hashable positive-count representation of one path usage."""

    items = []
    for event_key, raw_count in usage.items():
        count = int(raw_count)
        if count != raw_count:
            raise ValueError(
                f"NClose usage count must be integral: {event_key!r}={raw_count!r}"
            )
        if count < 0:
            raise ValueError(
                f"NClose usage count must be non-negative: {event_key!r}={count}"
            )
        if count > 0:
            items.append((event_key, count))
    return frozenset(items)


def add_usage_signatures(
    usage_a: UsageSignature,
    usage_b: UsageSignature,
) -> UsageSignature:
    """Add two NClose usage multisets element by element."""

    combined = Counter(dict(usage_a))
    combined.update(dict(usage_b))
    return usage_signature(combined)


def subtract_usage_signatures(
    whole: UsageSignature,
    part: UsageSignature,
) -> UsageSignature | None:
    """Return ``whole - part``, or ``None`` when part is not a submultiset."""

    whole_counter = Counter(dict(whole))
    part_counter = Counter(dict(part))
    if any(count > whole_counter[event_key]
           for event_key, count in part_counter.items()):
        return None
    whole_counter.subtract(part_counter)
    return usage_signature(whole_counter)


def intersect_usage_signatures(
    usage_a: UsageSignature,
    usage_b: UsageSignature,
) -> UsageSignature:
    """Return the element-wise intersection of two usage multisets."""

    return usage_signature(Counter(dict(usage_a)) & Counter(dict(usage_b)))


def find_additive_path_triples(
    candidate_columns: Sequence[int] | set[int],
    path_usage: Sequence[Mapping[EventKey, int]],
    weights: Sequence[float],
) -> list[AdditivePathTriple]:
    """Find all distinct-column additive triples, ordered by weight sum.

    A and B are treated symmetrically and must have different, non-empty
    usage multisets. They may share an event: if they do, AB must contain the
    sum of both occurrence counts. Exact-usage duplicate columns remain
    distinct paths, so each matching AB column is returned separately.
    """

    columns = sorted(set(map(int, candidate_columns)))
    if any(column < 0 or column >= len(path_usage) for column in columns):
        raise IndexError("candidate path column is outside path_usage")
    if any(column >= len(weights) for column in columns):
        raise IndexError("candidate path column is outside weights")

    signature_by_column = {
        column: usage_signature(path_usage[column])
        for column in columns
    }
    columns_by_signature: dict[UsageSignature, list[int]] = defaultdict(list)
    for column, signature in signature_by_column.items():
        if signature:
            columns_by_signature[signature].append(column)

    triples = []
    for offset, path_a in enumerate(columns):
        usage_a = signature_by_column[path_a]
        if not usage_a:
            continue
        for path_b in columns[offset + 1:]:
            usage_b = signature_by_column[path_b]
            if not usage_b or usage_a == usage_b:
                continue
            usage_ab = add_usage_signatures(usage_a, usage_b)
            for path_ab in columns_by_signature.get(usage_ab, ()):
                if path_ab in {path_a, path_b}:
                    continue
                triples.append(AdditivePathTriple(
                    path_a=path_a,
                    path_b=path_b,
                    path_ab=path_ab,
                    weight_a=float(weights[path_a]),
                    weight_b=float(weights[path_b]),
                    weight_ab=float(weights[path_ab]),
                ))

    triples.sort(key=lambda triple: (
        triple.total_weight,
        triple.path_a,
        triple.path_b,
        triple.path_ab,
    ))
    return triples


def find_chained_path_triples(
    candidate_columns: Sequence[int] | set[int],
    path_usage: Sequence[Mapping[EventKey, int]],
    weights: Sequence[float],
) -> list[ChainedPathTriple]:
    """Find strict A/BC/AB triples, ordered by their three-path weight sum.

    A candidate must admit the exact multiset decomposition::

        AB = A + B
        BC = B + C

    A, B, and C must all be non-empty.  A and BC must be disjoint, which also
    guarantees that the complete AB/BC overlap is exactly B.  This strictness
    prevents an incidental shared breakend from being interpreted as a path
    chain.  Repeated breakends are compared by their full occurrence counts.
    """

    columns = sorted(set(map(int, candidate_columns)))
    if any(column < 0 or column >= len(path_usage) for column in columns):
        raise IndexError("candidate path column is outside path_usage")
    if any(column >= len(weights) for column in columns):
        raise IndexError("candidate path column is outside weights")

    signature_by_column = {
        column: usage_signature(path_usage[column])
        for column in columns
    }

    triples = []
    for path_a in columns:
        usage_a = signature_by_column[path_a]
        if not usage_a:
            continue

        for path_ab in columns:
            if path_ab == path_a:
                continue
            usage_ab = signature_by_column[path_ab]
            usage_b = subtract_usage_signatures(usage_ab, usage_a)
            if not usage_b:
                continue

            for path_bc in columns:
                if path_bc in {path_a, path_ab}:
                    continue
                usage_bc = signature_by_column[path_bc]
                if not usage_bc:
                    continue
                if intersect_usage_signatures(usage_a, usage_bc):
                    continue

                usage_c = subtract_usage_signatures(usage_bc, usage_b)
                if not usage_c:
                    continue
                if intersect_usage_signatures(usage_ab, usage_bc) != usage_b:
                    continue

                triples.append(ChainedPathTriple(
                    path_a=path_a,
                    path_bc=path_bc,
                    path_ab=path_ab,
                    weight_a=float(weights[path_a]),
                    weight_bc=float(weights[path_bc]),
                    weight_ab=float(weights[path_ab]),
                ))

    triples.sort(key=lambda triple: (
        triple.total_weight,
        triple.path_a,
        triple.path_bc,
        triple.path_ab,
    ))
    return triples


def find_path_conflicts(
    candidate_columns: Sequence[int] | set[int],
    path_usage: Sequence[Mapping[EventKey, int]],
    weights: Sequence[float],
) -> list[PathConflict]:
    """Return A/B/AB and A/BC/AB conflicts in one weight-sorted queue."""

    conflicts = [
        PathConflict(
            conflict_type='A_B_AB',
            path_a=triple.path_a,
            path_b=triple.path_b,
            path_ab=triple.path_ab,
            role_a='A',
            role_b='B',
            role_ab='AB',
            weight_a=triple.weight_a,
            weight_b=triple.weight_b,
            weight_ab=triple.weight_ab,
        )
        for triple in find_additive_path_triples(
            candidate_columns, path_usage, weights
        )
    ]
    conflicts.extend(
        PathConflict(
            conflict_type='A_BC_AB',
            path_a=triple.path_a,
            path_b=triple.path_bc,
            path_ab=triple.path_ab,
            role_a='A',
            role_b='BC',
            role_ab='AB',
            weight_a=triple.weight_a,
            weight_b=triple.weight_bc,
            weight_ab=triple.weight_ab,
        )
        for triple in find_chained_path_triples(
            candidate_columns, path_usage, weights
        )
    )
    conflicts.sort(key=lambda conflict: (
        conflict.total_weight,
        conflict.path_a,
        conflict.path_b,
        conflict.path_ab,
        conflict.conflict_type,
    ))
    return conflicts


def choose_removal(
    path_a: int,
    path_b: int,
    *,
    a_feasible: bool,
    b_feasible: bool,
    a_error: float,
    b_error: float,
    both_feasible: bool | None = None,
) -> RemovalChoice:
    """Apply the requested A/B removal policy for one additive triple.

    ``a_error`` is the global error after removing A, and likewise for B.
    When a forced single choice is required, the lower-error removal wins.
    Exact ties are resolved by the smaller removed column index.
    """

    path_a = int(path_a)
    path_b = int(path_b)
    if path_a == path_b:
        raise ValueError("A and B must be different path columns")
    if not isfinite(a_error) and not isfinite(b_error):
        raise ValueError("at least one single-path removal must have a finite error")

    def lower_error_single(reason: str) -> RemovalChoice:
        candidates = [
            (float(a_error), path_a),
            (float(b_error), path_b),
        ]
        _, selected = min(candidates, key=lambda item: (item[0], item[1]))
        return RemovalChoice((selected,), reason)

    if a_feasible and b_feasible:
        if both_feasible is None:
            raise ValueError(
                "both_feasible is required when both single removals are feasible"
            )
        if both_feasible:
            return RemovalChoice(
                tuple(sorted((path_a, path_b))),
                "joint_removal_feasible",
            )
        return lower_error_single("lower_error_after_joint_failure")

    if a_feasible:
        return RemovalChoice((path_a,), "only_a_feasible")
    if b_feasible:
        return RemovalChoice((path_b,), "only_b_feasible")
    return lower_error_single("lower_error_forced")
