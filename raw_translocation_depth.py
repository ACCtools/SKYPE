"""Depth helpers for raw-read reciprocal-NClose filtering."""

from __future__ import annotations

import math

from denoised_relative_error import estimate_noise_sigma_by_chromosome


RAW_TRANSLOCATION_DEPTH_FLANK = 500_000
DEPTH_BALANCED_NOISE_SIGMA_MULTIPLIER = 1.5


def weighted_mean_depth(depth_by_chrom, chrom, start, end):
    """Return the length-weighted mean over one inclusive reference interval."""

    if chrom not in depth_by_chrom:
        return None
    start = int(start)
    end = int(end)
    if end < start:
        return None

    weighted_sum = 0.0
    weight_sum = 0
    for window_start, window_end, mean_depth in depth_by_chrom[chrom]:
        overlap_start = max(start, int(window_start))
        overlap_end = min(end, int(window_end))
        if overlap_end < overlap_start:
            continue
        weight = overlap_end - overlap_start + 1
        weighted_sum += float(mean_depth) * weight
        weight_sum += weight

    if weight_sum == 0:
        return None
    return weighted_sum / weight_sum


def breakpoint_flank_depths(
    depth_by_chrom,
    chrom,
    inner_start,
    inner_end,
    chrom_length,
    flank=RAW_TRANSLOCATION_DEPTH_FLANK,
):
    """Measure disjoint one-sided depth windows around a breakpoint interval."""

    inner_start, inner_end = sorted((int(inner_start), int(inner_end)))
    chrom_length = int(chrom_length)
    flank = int(flank)
    if flank < 1 or chrom_length < 1:
        return {"front": None, "back": None}

    front_start = max(1, inner_start - flank)
    front_end = min(chrom_length, inner_start - 1)
    back_start = max(1, inner_end + 1)
    back_end = min(chrom_length, inner_end + flank)
    return {
        "front": weighted_mean_depth(
            depth_by_chrom,
            chrom,
            front_start,
            front_end,
        ),
        "back": weighted_mean_depth(
            depth_by_chrom,
            chrom,
            back_start,
            back_end,
        ),
    }


def build_raw_point_depths(
    candidates,
    depth_by_chrom,
    chromosome_lengths,
    flank=RAW_TRANSLOCATION_DEPTH_FLANK,
):
    """Build point-a/point-b flank depths from reciprocal candidate sides."""

    point_depths = {}
    for candidate in candidates:
        pair_id = int(candidate["pair_id"])
        side_records = list(candidate.get("side_records", ()))
        for side_index, point_name in enumerate(("point_a", "point_b")):
            if side_index >= len(side_records):
                point_depths.setdefault(pair_id, {})[point_name] = {
                    "front": None,
                    "back": None,
                }
                continue
            side = side_records[side_index]
            chrom = side["chrom"]
            point_depths.setdefault(pair_id, {})[point_name] = \
                breakpoint_flank_depths(
                    depth_by_chrom,
                    chrom,
                    side["inner_st"],
                    side["inner_nd"],
                    chromosome_lengths.get(chrom, 0),
                    flank=flank,
                )
    return point_depths


def depth_pair_is_balanced(depth_pair, absolute_diff_threshold):
    """Judge one breakpoint's flank pair against the stage-23 noise threshold."""

    if not depth_pair:
        return False
    front = depth_pair.get("front")
    back = depth_pair.get("back")
    if front is None or back is None:
        return False
    front = float(front)
    back = float(back)
    if not (math.isfinite(front) and math.isfinite(back)):
        return False
    low, high = sorted((front, back))
    if low < 0:
        return False
    if low == 0:
        return high == 0
    return high - low < float(absolute_diff_threshold)


def breakpoint_is_depth_balanced(
    depth_by_chrom,
    chrom,
    inner_start,
    inner_end,
    chrom_length,
    absolute_diff_threshold,
    flank=RAW_TRANSLOCATION_DEPTH_FLANK,
):
    """Measure one breakpoint's flanks and judge them in a single call.

    Pass ``inner_start == inner_end`` for a point breakpoint.
    """

    depth_pair = breakpoint_flank_depths(
        depth_by_chrom,
        chrom,
        inner_start,
        inner_end,
        chrom_length,
        flank=flank,
    )
    return depth_pair_is_balanced(depth_pair, absolute_diff_threshold)


def global_depth_noise_threshold(
    depth_by_chrom,
    multiplier=DEPTH_BALANCED_NOISE_SIGMA_MULTIPLIER,
):
    """Fallback threshold when stage 01 does not pass its censat-aware one."""

    coordinates = []
    depths = []
    for chrom, windows in depth_by_chrom.items():
        for window_start, _, mean_depth in windows:
            depth = float(mean_depth)
            if not math.isfinite(depth):
                continue
            coordinates.append((str(chrom), int(window_start)))
            depths.append(depth)
    sigma = estimate_noise_sigma_by_chromosome(coordinates, depths)["__global__"]
    return multiplier * sigma


def depth_pair_mean(depth_pair):
    if not depth_pair:
        return None
    front = depth_pair.get("front")
    back = depth_pair.get("back")
    if front is None or back is None:
        return None
    if not (math.isfinite(front) and math.isfinite(back)):
        return None
    return (float(front) + float(back)) / 2
