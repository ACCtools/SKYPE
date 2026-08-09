"""Cluster reconstructed paths by shared breakend-graph reference spans.

The compact graph path records chromosome/direction transitions but stage 30
normally reduces each continuous reference run to a length-only karyotype
piece.  This module retains the genomic coordinates of those runs so paths
that reuse at least ``min_overlap`` bases of the same reference interval can
be shown as one ambiguity component.

Intervals use the same zero-based, half-open convention as PAF target spans.
Strand is retained for plotting but deliberately ignored for overlap tests: a
forward and reverse traversal of the same reference interval still share that
reference sequence.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from typing import Hashable, Iterable, Mapping, Sequence


NODE_NAME = 1
CHR_CHANGE_IDX = 2
DIR_CHANGE_IDX = 3

CTG_DIR = 4
CHR_NAM = 5
CHR_STR = 7
CHR_END = 8


@dataclass(frozen=True, order=True)
class ReferenceSpan:
    """One ordered, continuous reference run traversed by a graph path."""

    chrom: str
    start: int
    end: int
    strand: str
    source: str = "graph"

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True, order=True)
class ReferenceOverlap:
    """One same-chromosome intersection between two path coverages."""

    chrom: str
    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start


@dataclass(frozen=True)
class PathOverlapEdge:
    """A qualifying edge in the path reference-overlap graph."""

    path_a: Hashable
    path_b: Hashable
    overlaps: tuple[ReferenceOverlap, ...]

    @property
    def total_overlap(self) -> int:
        return sum(overlap.length for overlap in self.overlaps)

    @property
    def maximum_overlap(self) -> int:
        return max((overlap.length for overlap in self.overlaps), default=0)


@dataclass(frozen=True)
class ReferencePathCluster:
    """One non-singleton connected component of overlapping paths."""

    cluster_id: str
    members: tuple[Hashable, ...]
    edges: tuple[PathOverlapEdge, ...]


@dataclass(frozen=True)
class SharedReferenceBlock:
    """Atomic reference interval carried by two or more clustered paths."""

    block_id: str
    chrom: str
    start: int
    end: int
    carriers: tuple[Hashable, ...]

    @property
    def length(self) -> int:
        return self.end - self.start


def _padded_compact_path(index_path: Sequence[Sequence]) -> list[tuple]:
    """Copy a compact path and pad its endpoint records without mutating it."""

    path = [tuple(step) for step in index_path]
    if len(path) < 3:
        raise ValueError("A compact graph path must contain two endpoints and a real node")
    if len(path[0]) < 4:
        path[0] = (0,) + path[0]
    if len(path[-1]) < 4:
        path[-1] = (0,) + path[-1]
    if not isinstance(path[0][NODE_NAME], str) or not path[0][NODE_NAME]:
        raise ValueError(f"Invalid compact-path start endpoint: {path[0]!r}")
    if not isinstance(path[1][NODE_NAME], int):
        raise ValueError(f"Compact path has no first real node: {path[1]!r}")
    if not isinstance(path[-2][NODE_NAME], int):
        raise ValueError(f"Compact path has no final real node: {path[-2]!r}")
    return path


def extract_ordered_reference_spans(
    index_path: Sequence[Sequence],
    ppc_data: Sequence[Sequence],
    *,
    forced_break_edges: Iterable[tuple] = (),
    interrupt_spans_by_edge: Mapping[tuple, Iterable[ReferenceSpan]] | None = None,
) -> tuple[ReferenceSpan, ...]:
    """Recover continuous reference runs between compact-path breakends.

    ``forced_break_edges`` can be used for graph-only deletion edges whose
    chromosome/direction counters do not increment.  Edge keys have the stage
    30 form ``((graph_state, previous_node), (graph_state, current_node))``.
    ``interrupt_spans_by_edge`` inserts coordinate-bearing CTG_IN fragments in
    the same position where Virtual SKY inserts their chromosome pieces.
    """

    path = _padded_compact_path(index_path)
    forced_break_edges = set(forced_break_edges)
    interrupt_spans_by_edge = {
        edge: tuple(edge_spans)
        for edge, edge_spans in (interrupt_spans_by_edge or {}).items()
    }

    direction = "+" if path[0][NODE_NAME].endswith("f") else "-"
    first_node = ppc_data[path[1][NODE_NAME]]
    chrom = str(first_node[CHR_NAM])
    anchor = int(first_node[CHR_STR] if direction == "+" else first_node[CHR_END])
    spans: list[ReferenceSpan] = []

    def append_span(boundary: int) -> None:
        start, end = sorted((anchor, int(boundary)))
        if end > start:
            spans.append(ReferenceSpan(chrom, start, end, direction))

    for index in range(1, len(path) - 1):
        previous_id = path[index - 1][NODE_NAME]
        current_id = path[index][NODE_NAME]
        if not isinstance(previous_id, int) or not isinstance(current_id, int):
            continue

        edge_key = (
            (path[index - 1][0], previous_id),
            (path[index][0], current_id),
        )
        transition = (
            path[index][CHR_CHANGE_IDX] > path[index - 1][CHR_CHANGE_IDX]
            or path[index][DIR_CHANGE_IDX] > path[index - 1][DIR_CHANGE_IDX]
            or edge_key in forced_break_edges
            or edge_key in interrupt_spans_by_edge
        )
        if not transition:
            continue

        previous_node = ppc_data[previous_id]
        append_span(
            previous_node[CHR_END] if direction == "+" else previous_node[CHR_STR]
        )
        spans.extend(interrupt_spans_by_edge.get(edge_key, ()))

        current_node = ppc_data[current_id]
        if current_id > previous_id:
            direction = str(current_node[CTG_DIR])
        else:
            direction = "-" if current_node[CTG_DIR] == "+" else "+"
        chrom = str(current_node[CHR_NAM])
        anchor = int(
            current_node[CHR_STR] if direction == "+" else current_node[CHR_END]
        )

    last_node = ppc_data[path[-2][NODE_NAME]]
    append_span(last_node[CHR_END] if direction == "+" else last_node[CHR_STR])
    return tuple(spans)


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    """Return the union of zero-based half-open intervals."""

    merged: list[list[int]] = []
    for raw_start, raw_end in sorted((int(start), int(end)) for start, end in intervals):
        if raw_end <= raw_start:
            continue
        if merged and raw_start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], raw_end)
        else:
            merged.append([raw_start, raw_end])
    return tuple((start, end) for start, end in merged)


def merged_span_coverage(
    spans: Iterable[ReferenceSpan],
) -> dict[str, tuple[tuple[int, int], ...]]:
    """Union path spans by chromosome; repeated traversals count once."""

    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for span in spans:
        by_chrom[span.chrom].append((span.start, span.end))
    return {
        chrom: merge_intervals(intervals)
        for chrom, intervals in by_chrom.items()
        if intervals
    }


def intersect_merged_intervals(
    intervals_a: Sequence[tuple[int, int]],
    intervals_b: Sequence[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    """Intersect two sorted, internally disjoint interval sequences."""

    intersections: list[tuple[int, int]] = []
    index_a = 0
    index_b = 0
    while index_a < len(intervals_a) and index_b < len(intervals_b):
        start = max(intervals_a[index_a][0], intervals_b[index_b][0])
        end = min(intervals_a[index_a][1], intervals_b[index_b][1])
        if end > start:
            intersections.append((start, end))
        if intervals_a[index_a][1] < intervals_b[index_b][1]:
            index_a += 1
        else:
            index_b += 1
    return tuple(intersections)


def path_reference_overlaps(
    spans_a: Iterable[ReferenceSpan],
    spans_b: Iterable[ReferenceSpan],
    *,
    minimum_overlap: int = 1,
) -> tuple[ReferenceOverlap, ...]:
    """Return qualifying intersections of individual continuous graph runs.

    The length cutoff is applied before overlapping/repeated intersections are
    unioned.  Consequently two separate 60-kb traversals cannot combine into
    one 120-kb cluster edge, while repeated traversal of the same qualifying
    interval is still reported only once.
    """

    if minimum_overlap <= 0:
        raise ValueError("minimum_overlap must be positive")
    spans_a = tuple(spans_a)
    spans_b = tuple(spans_b)
    qualifying_by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for span_a in spans_a:
        for span_b in spans_b:
            if span_a.chrom != span_b.chrom:
                continue
            start = max(span_a.start, span_b.start)
            end = min(span_a.end, span_b.end)
            if end - start >= minimum_overlap:
                qualifying_by_chrom[span_a.chrom].append((start, end))
    return tuple(
        ReferenceOverlap(chrom, start, end)
        for chrom in sorted(qualifying_by_chrom)
        for start, end in merge_intervals(qualifying_by_chrom[chrom])
    )


def is_full_reference_path(
    spans: Iterable[ReferenceSpan],
    *,
    has_nclose: bool,
    chromosome_lengths: Mapping[str, int],
    minimum_ratio: float,
) -> tuple[bool, str | None, float]:
    """Classify the existing no-NClose, >=reference-ratio normal paths."""

    coverage = merged_span_coverage(spans)
    ratios = []
    for chrom, intervals in coverage.items():
        chromosome_length = int(chromosome_lengths.get(chrom, 0))
        if chromosome_length <= 0:
            continue
        covered_bases = sum(end - start for start, end in intervals)
        ratios.append((covered_bases / chromosome_length, chrom))
    if not ratios:
        return False, None, 0.0
    ratio, chrom = max(ratios)
    is_single_chrom_reference = len(coverage) == 1
    return (
        not has_nclose
        and is_single_chrom_reference
        and ratio >= float(minimum_ratio)
    ), chrom, ratio


def build_path_overlap_edges(
    path_spans: Mapping[Hashable, Iterable[ReferenceSpan]],
    *,
    minimum_overlap: int,
) -> tuple[PathOverlapEdge, ...]:
    """Build graph edges when any one contiguous overlap meets the cutoff."""

    if minimum_overlap <= 0:
        raise ValueError("minimum_overlap must be positive")
    path_keys = sorted(path_spans, key=_stable_key)
    edges = []
    for left_index, path_a in enumerate(path_keys):
        for path_b in path_keys[left_index + 1:]:
            qualifying = path_reference_overlaps(
                path_spans[path_a], path_spans[path_b],
                minimum_overlap=minimum_overlap,
            )
            if qualifying:
                edges.append(PathOverlapEdge(path_a, path_b, qualifying))
    return tuple(edges)


def _stable_key(value: Hashable) -> tuple[str, str]:
    if isinstance(value, int) and not isinstance(value, bool):
        return "00-int", f"{value:020d}"
    return type(value).__name__, repr(value)


def connected_path_components(
    path_keys: Iterable[Hashable],
    edges: Iterable[PathOverlapEdge],
) -> tuple[tuple[Hashable, ...], ...]:
    """Return deterministic connected components, including singletons."""

    ordered_keys = sorted(set(path_keys), key=_stable_key)
    parent = {path: path for path in ordered_keys}

    def find(path):
        root = path
        while parent[root] != root:
            root = parent[root]
        while parent[path] != path:
            next_path = parent[path]
            parent[path] = root
            path = next_path
        return root

    def union(path_a, path_b) -> None:
        root_a = find(path_a)
        root_b = find(path_b)
        if root_a == root_b:
            return
        if _stable_key(root_b) < _stable_key(root_a):
            root_a, root_b = root_b, root_a
        parent[root_b] = root_a

    edge_list = tuple(edges)
    for edge in edge_list:
        if edge.path_a not in parent or edge.path_b not in parent:
            raise ValueError("Overlap edge contains a path absent from path_keys")
        union(edge.path_a, edge.path_b)

    components: dict[Hashable, list[Hashable]] = defaultdict(list)
    for path in ordered_keys:
        components[find(path)].append(path)
    return tuple(sorted(
        (
            tuple(sorted(component, key=_stable_key))
            for component in components.values()
        ),
        key=lambda component: (-len(component), tuple(_stable_key(x) for x in component)),
    ))


def cluster_reference_paths(
    path_spans: Mapping[Hashable, Iterable[ReferenceSpan]],
    *,
    minimum_overlap: int,
) -> tuple[tuple[ReferencePathCluster, ...], tuple[Hashable, ...]]:
    """Return labeled non-singleton clusters and unclustered singleton paths."""

    edges = build_path_overlap_edges(path_spans, minimum_overlap=minimum_overlap)
    components = connected_path_components(path_spans, edges)
    clusters = []
    singletons = []
    for component in components:
        if len(component) < 2:
            singletons.extend(component)
            continue
        member_set = set(component)
        component_edges = tuple(
            edge
            for edge in edges
            if edge.path_a in member_set and edge.path_b in member_set
        )
        clusters.append(ReferencePathCluster("", component, component_edges))

    clusters.sort(
        key=lambda cluster: (
            -len(cluster.members),
            tuple(_stable_key(member) for member in cluster.members),
        )
    )
    clusters = [
        ReferencePathCluster(f"C{index:02d}", cluster.members, cluster.edges)
        for index, cluster in enumerate(clusters, start=1)
    ]
    return tuple(clusters), tuple(sorted(singletons, key=_stable_key))


def build_shared_reference_blocks(
    cluster: ReferencePathCluster,
    path_spans: Mapping[Hashable, Iterable[ReferenceSpan]],
) -> tuple[SharedReferenceBlock, ...]:
    """Split qualifying overlap regions into carrier-aware reference blocks."""

    missing_paths = set(cluster.members) - set(path_spans)
    if missing_paths:
        raise ValueError(f"Missing path spans for cluster members: {missing_paths}")
    chromosomes = sorted({
        overlap.chrom
        for edge in cluster.edges
        for overlap in edge.overlaps
    })
    raw_blocks: list[tuple[str, int, int, tuple[Hashable, ...]]] = []
    for chrom in chromosomes:
        boundaries = sorted({
            coordinate
            for edge in cluster.edges
            for overlap in edge.overlaps
            if overlap.chrom == chrom
            for coordinate in (overlap.start, overlap.end)
        })
        for start, end in zip(boundaries, boundaries[1:]):
            if end <= start:
                continue
            carrier_set = set()
            for edge in cluster.edges:
                if any(
                    overlap.chrom == chrom
                    and overlap.start <= start
                    and end <= overlap.end
                    for overlap in edge.overlaps
                ):
                    carrier_set.update((edge.path_a, edge.path_b))
            carriers = tuple(
                path for path in cluster.members if path in carrier_set
            )
            if len(carriers) < 2:
                continue
            if (
                raw_blocks
                and raw_blocks[-1][0] == chrom
                and raw_blocks[-1][2] == start
                and raw_blocks[-1][3] == carriers
            ):
                previous = raw_blocks[-1]
                raw_blocks[-1] = (chrom, previous[1], end, carriers)
            else:
                raw_blocks.append((chrom, start, end, carriers))

    return tuple(
        SharedReferenceBlock(
            f"R{index:02d}", chrom, start, end, carriers
        )
        for index, (chrom, start, end, carriers) in enumerate(raw_blocks, start=1)
    )
