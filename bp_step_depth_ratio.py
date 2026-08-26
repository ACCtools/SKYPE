"""Breakpoint step/depth-ratio annotations for native SKYPE calls."""

from __future__ import annotations

from collections import defaultdict

import numpy as np


BP_STEP_DEPTH_RATIO_B = "BP_STEP_DEPTH_RATIO_B"
BP_STEP_DEPTH_RATIO_PREDICT_B = "BP_STEP_DEPTH_RATIO_PREDICT_B"
DEFAULT_STEP_WINDOW = 500_000
MIN_EVENT_DEPTH = 1e-6


def _point_overlaps(intervals_by_chrom, chrom, coord):
    return any(
        max(int(start), int(coord)) < min(int(end), int(coord) + 1)
        for start, end in intervals_by_chrom.get(chrom, ())
    )


class BreakpointStepDepthRatio:
    """Calculate the former stage-23 signed ratio without filtering calls."""

    def __init__(
        self,
        coordinates,
        observed_depth,
        predicted_depth,
        censat_intervals,
        window=DEFAULT_STEP_WINDOW,
    ):
        self.coordinates = list(coordinates)
        self.observed_depth = np.asarray(observed_depth)
        self.predicted_depth = np.asarray(predicted_depth)
        self.censat_intervals = censat_intervals
        self.window = int(window)

        expected_shape = (len(self.coordinates),)
        if self.observed_depth.shape != expected_shape:
            raise ValueError(
                "Observed clean-depth shape does not match its coordinates: "
                f"{self.observed_depth.shape} != {expected_shape}"
            )
        if self.predicted_depth.shape != expected_shape:
            raise ValueError(
                "Predicted clean-depth shape does not match its coordinates: "
                f"{self.predicted_depth.shape} != {expected_shape}"
            )

        self.bins_by_chrom = defaultdict(list)
        for index, (chrom, start) in enumerate(self.coordinates):
            self.bins_by_chrom[chrom].append((index, int(start)))
        for bins in self.bins_by_chrom.values():
            bins.sort(key=lambda item: item[1])

    def _ratio(self, signal, chrom, coord, expected_high_side, event_depth):
        event_depth = float(event_depth)
        coord = int(coord)
        if event_depth < MIN_EVENT_DEPTH:
            return None
        if expected_high_side not in {"left", "right"}:
            return None
        if _point_overlaps(self.censat_intervals, chrom, coord):
            return None

        bins = self.bins_by_chrom.get(chrom, ())
        left = [index for index, start in bins if coord - self.window <= start < coord]
        right = [index for index, start in bins if coord <= start < coord + self.window]
        if not left or not right:
            return None

        left_mean = float(signal[left].mean())
        right_mean = float(signal[right].mean())
        if expected_high_side == "left":
            step = left_mean - right_mean
        else:
            step = right_mean - left_mean
        return step / event_depth

    def endpoint(self, chrom, coord, expected_high_side, event_depth):
        """Return observed and predicted ratios for one breakpoint."""

        return (
            self._ratio(
                self.observed_depth, chrom, coord, expected_high_side, event_depth
            ),
            self._ratio(
                self.predicted_depth, chrom, coord, expected_high_side, event_depth
            ),
        )

    def pair(self, endpoints, event_depth):
        """Return INFO-ready observed/predicted arrays in endpoint order."""

        values = [
            self.endpoint(chrom, coord, expected_high_side, event_depth)
            for chrom, coord, expected_high_side in endpoints
        ]
        if len(values) != 2:
            raise ValueError(f"Expected two breakpoints, found {len(values)}")
        return {
            BP_STEP_DEPTH_RATIO_B: [values[0][0], values[1][0]],
            BP_STEP_DEPTH_RATIO_PREDICT_B: [values[0][1], values[1][1]],
        }


def bnd_expected_high_side(direction, endpoint_role):
    """Return the reference side covered by an emitted BND anchor."""

    if direction not in {"+", "-"}:
        return None
    if endpoint_role == "exit":
        return "left" if direction == "+" else "right"
    if endpoint_role == "entry":
        return "right" if direction == "+" else "left"
    raise ValueError(f"Unknown BND endpoint role: {endpoint_role}")


def add_ratio_info(info, ratio_info, reverse=False):
    """Add the two ratio arrays, optionally making the second endpoint local."""

    if ratio_info is None:
        return
    for info_id in (
        BP_STEP_DEPTH_RATIO_B,
        BP_STEP_DEPTH_RATIO_PREDICT_B,
    ):
        values = list(ratio_info[info_id])
        info[info_id] = list(reversed(values)) if reverse else values
