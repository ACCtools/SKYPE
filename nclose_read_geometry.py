"""Pure coordinate helpers for strand-aware NClose raw-read counting."""

from __future__ import annotations

from collections.abc import Iterable


def _strand_step(strand: str) -> int:
    if strand == "+":
        return 1
    if strand == "-":
        return -1
    raise ValueError(f"Unsupported strand: {strand!r}")


def _side_index(side_index: int) -> int:
    side_index = int(side_index)
    if side_index not in (0, 1):
        raise ValueError(f"NClose side index must be 0 or 1: {side_index}")
    return side_index


def split_nclose_flanks(query_gap: int, window: int) -> tuple[int, int] | None:
    """Split the non-gap part of a fixed query window across both endpoints."""

    query_gap = max(0, int(query_gap))
    window = int(window)
    if window < 0:
        raise ValueError(f"NClose count window must be non-negative: {window}")
    if query_gap > window:
        return None
    remaining = window - query_gap
    first = remaining // 2
    return first, remaining - first


def reference_support_points(
    shared_anchor: int,
    strand: str,
    side_index: int,
    span: int,
    chrom_length: int,
) -> tuple[int, int]:
    """Return query-ordered normal-reference points sharing a junction anchor.

    The first NClose endpoint starts at its d2 anchor and continues along its
    alignment strand.  The second endpoint approaches its d2 anchor along its
    alignment strand.  Only the non-shared point is clipped at a chromosome
    boundary, so the d1/d3 route always retains the exact d2 anchor.
    """

    side_index = _side_index(side_index)
    step = _strand_step(strand)
    shared_anchor = int(shared_anchor)
    span = int(span)
    chrom_length = int(chrom_length)
    if span < 0:
        raise ValueError(f"Reference support span must be non-negative: {span}")
    if chrom_length < 1:
        raise ValueError(f"Chromosome length must be positive: {chrom_length}")
    if not 1 <= shared_anchor <= chrom_length:
        raise ValueError(
            f"Shared d2 anchor {shared_anchor} is outside chromosome length "
            f"{chrom_length}"
        )

    if side_index == 0:
        query_start = shared_anchor
        query_end = shared_anchor + step * span
        query_end = max(1, min(query_end, chrom_length))
    else:
        query_start = shared_anchor - step * span
        query_start = max(1, min(query_start, chrom_length))
        query_end = shared_anchor
    return query_start, query_end


def reference_support_span_rows(
    pair_id: int,
    count_idx: int,
    shared_anchor_row: list,
    strand: str,
    side_index: int,
    span: int,
    chrom_length: int,
) -> list[list]:
    """Build two BAM-anchor rows for a strand-aware d1 or d3 route."""

    if len(shared_anchor_row) < 8:
        raise ValueError("Malformed shared d2 anchor row")
    chrom = str(shared_anchor_row[3])
    shared_anchor = int(shared_anchor_row[6])
    query_start, query_end = reference_support_points(
        shared_anchor,
        strand,
        side_index,
        span,
        chrom_length,
    )
    fetch_start, fetch_end = sorted((query_start, query_end))
    expected_positive = strand == "+"
    return [
        [
            int(pair_id), int(count_idx), 1, chrom,
            fetch_start, fetch_end, query_start, expected_positive,
        ],
        [
            int(pair_id), int(count_idx), 2, chrom,
            fetch_start, fetch_end, query_end, expected_positive,
        ],
    ]


def aligned_depth_interval(
    breakend_coord: int,
    strand: str,
    side_index: int,
    flank: int,
    chrom_length: int,
) -> tuple[int, int] | None:
    """Return the 1-based inclusive depth interval on the aligned NClose side."""

    side_index = _side_index(side_index)
    _strand_step(strand)  # validate strand before selecting its genomic side
    breakend_coord = int(breakend_coord)
    flank = int(flank)
    chrom_length = int(chrom_length)
    if flank <= 0 or chrom_length < 1:
        return None

    # PAF target coordinates are 0-based half-open.  Convert exactly the
    # aligned-side half-open interval to a 1-based inclusive PanDepth interval.
    lower_reference_side = (
        (side_index == 0 and strand == "+")
        or (side_index == 1 and strand == "-")
    )
    if lower_reference_side:
        start0, end0 = breakend_coord - flank, breakend_coord
    else:
        start0, end0 = breakend_coord, breakend_coord + flank
    start0 = max(0, min(start0, chrom_length))
    end0 = max(0, min(end0, chrom_length))
    if end0 <= start0:
        return None
    return start0 + 1, end0


def interval_weighted_mean(
    depth_rows: Iterable[tuple[int, int, float]],
    interval: tuple[int, int] | None,
) -> float | None:
    """Length-weighted mean for inclusive PanDepth rows overlapping interval."""

    if interval is None:
        return None
    start, end = map(int, interval)
    weighted_sum = 0.0
    weight_sum = 0
    for win_start, win_end, mean_depth in depth_rows:
        overlap_start = max(start, int(win_start))
        overlap_end = min(end, int(win_end))
        if overlap_end < overlap_start:
            continue
        weight = overlap_end - overlap_start + 1
        weighted_sum += float(mean_depth) * weight
        weight_sum += weight
    if weight_sum == 0:
        return None
    return weighted_sum / weight_sum
