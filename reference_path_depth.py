"""Observed, fitted, and per-path depth tracks for reference-sharing clusters.

Stage 22 builds the NNLS design matrix ``A`` from one simulated depth file per
separated path piece (``21_pat_depth/<key>.win.stat.gz``) and the observed
target ``B`` from the cell line's own window statistics.  The reference-span
cluster figures only knew each path's fitted scalar depth, so this module
recovers the same per-window quantities the solver actually saw:

``observed``   the cell line's real coverage per window, in copy-number units
``model``      ``A @ w`` for every fitted column (``predict_B*.npy``)
``path``       one member path's own ``A`` column scaled by its fitted weight

All three share the window grid of the main stat file, so a cluster member's
contribution can be read against the reference depth at the same coordinate.

Window statistics are 1-based inclusive; reference spans are zero-based
half-open like PAF targets.  Conversions happen at the boundaries below.
"""

from __future__ import annotations

import logging
import os
import pickle as pkl
from dataclasses import dataclass
from typing import Hashable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd


WIN_STAT_COLUMNS = (
    "chr", "st", "nd", "length", "covsite", "totaldepth", "cov", "meandepth",
)


@dataclass(frozen=True)
class DepthBinIndex:
    """Depth windows of the main stat file, kept in file order."""

    chroms: tuple[str, ...]
    starts: np.ndarray
    ends: np.ndarray
    row_by_key: Mapping[tuple[str, int], int]
    rows_by_chrom: Mapping[str, np.ndarray]

    @property
    def size(self) -> int:
        return len(self.starts)

    def rows_for_span(self, chrom: str, start: int, end: int) -> np.ndarray:
        """Rows of ``chrom`` intersecting the zero-based half-open span."""

        rows = self.rows_by_chrom.get(chrom)
        if rows is None or end <= start:
            return np.empty(0, dtype=np.int64)
        keep = (self.ends[rows] > int(start)) & (self.starts[rows] < int(end))
        return rows[keep]

    def overlap_lengths(self, rows: np.ndarray, start: int, end: int) -> np.ndarray:
        """Base overlap of each row with the zero-based half-open span."""

        if rows.size == 0:
            return np.empty(0, dtype=np.float64)
        left = np.maximum(self.starts[rows], int(start))
        right = np.minimum(self.ends[rows], int(end))
        return np.maximum(right - left, 0).astype(np.float64)


@dataclass(frozen=True)
class RegionDepthSummary:
    """Length-weighted depth attribution over one reference interval."""

    label: str
    chrom: str
    start: int
    end: int
    observed_cn: float
    observed_sd: float
    model_cn: float | None
    member_cn: Mapping[Hashable, float]
    carriers: tuple[Hashable, ...]

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def cluster_cn(self) -> float:
        return float(sum(self.member_cn.values()))

    @property
    def other_cn(self) -> float | None:
        if self.model_cn is None:
            return None
        return max(0.0, self.model_cn - self.cluster_cn)


@dataclass(frozen=True)
class ClusterDepthTracks:
    """Per-window observed/model/member depth for one cluster, in CN units."""

    bin_index: DepthBinIndex
    observed_cn: np.ndarray
    model_cn: np.ndarray | None
    member_cn: Mapping[Hashable, np.ndarray]
    cluster_rows: np.ndarray
    chrom_rows: Mapping[str, np.ndarray]

    @property
    def cluster_cn(self) -> np.ndarray:
        total = np.zeros(self.bin_index.size, dtype=np.float64)
        for values in self.member_cn.values():
            total += values
        return total

    @property
    def other_cn(self) -> np.ndarray | None:
        if self.model_cn is None:
            return None
        return np.maximum(self.model_cn - self.cluster_cn, 0.0)


def build_depth_bin_index(depth_frame: pd.DataFrame) -> DepthBinIndex:
    """Index the main stat windows once for every later coordinate lookup."""

    chroms = tuple(str(value) for value in depth_frame["chr"].to_numpy())
    # win.stat rows are 1-based inclusive; store them as zero-based half-open.
    starts = depth_frame["st"].to_numpy(dtype=np.int64) - 1
    ends = depth_frame["nd"].to_numpy(dtype=np.int64)
    row_by_key = {
        (chrom, int(start) + 1): row
        for row, (chrom, start) in enumerate(zip(chroms, starts))
    }
    rows_by_chrom: dict[str, list[int]] = {}
    for row, chrom in enumerate(chroms):
        rows_by_chrom.setdefault(chrom, []).append(row)
    return DepthBinIndex(
        chroms=chroms,
        starts=starts,
        ends=ends,
        row_by_key=row_by_key,
        rows_by_chrom={
            chrom: np.asarray(rows, dtype=np.int64)
            for chrom, rows in rows_by_chrom.items()
        },
    )


def observed_copy_number(
    depth_frame: pd.DataFrame, mean_depth: float
) -> np.ndarray:
    """Convert the cell line's own window coverage into copy-number units."""

    if mean_depth <= 0:
        raise ValueError("mean_depth must be positive")
    return depth_frame["meandepth"].to_numpy(dtype=np.float64) / mean_depth * 2.0


def load_model_copy_number(
    prefix: str,
    fig_prefix: str,
    bin_index: DepthBinIndex,
    mean_depth: float,
) -> np.ndarray | None:
    """Restore ``predict_B`` to main-stat window order, in copy-number units.

    Stage 22 stores ``B`` as the filtered windows followed by the excluded
    (CENSAT/over-coverage) windows, so ``predict_B`` uses that same permutation.
    ``23_input.pkl`` keeps the filtered half; the excluded half is whatever the
    main stat file has left, in file order.
    """

    predict_path = os.path.join(prefix, f"predict_B{fig_prefix}.npy")
    matrix_input_path = os.path.join(prefix, "23_input.pkl")
    if not (os.path.isfile(predict_path) and os.path.isfile(matrix_input_path)):
        logging.warning(
            "Reference cluster depth: no fitted depth track (%s missing)",
            predict_path if not os.path.isfile(predict_path) else matrix_input_path,
        )
        return None

    with open(matrix_input_path, "rb") as handle:
        filtered_keys = pkl.load(handle)[0]
    filtered_keys = [(str(chrom), int(start)) for chrom, start in filtered_keys]
    filtered_set = set(filtered_keys)
    excluded_keys = [
        key for key in bin_index.row_by_key if key not in filtered_set
    ]
    ordered_keys = filtered_keys + excluded_keys

    predicted = np.maximum(np.load(predict_path), 0.0)
    if len(ordered_keys) != bin_index.size or len(predicted) != bin_index.size:
        logging.warning(
            "Reference cluster depth: fitted depth length mismatch "
            "(%d windows, %d ordered keys, %d predictions)",
            bin_index.size, len(ordered_keys), len(predicted),
        )
        return None

    model = np.zeros(bin_index.size, dtype=np.float64)
    rows = np.asarray([bin_index.row_by_key[key] for key in ordered_keys])
    model[rows] = predicted
    return model / mean_depth * 2.0


class PathDepthProfiles:
    """Rebuild a path's ``A`` matrix row from its separated depth pieces."""

    def __init__(
        self,
        prefix: str,
        path_key_ints: Mapping[str, Sequence[int]],
        bin_index: DepthBinIndex,
    ) -> None:
        self._depth_folder = os.path.join(prefix, "21_pat_depth")
        self._path_key_ints = path_key_ints
        self._bin_index = bin_index
        self._piece_cache: dict[int, tuple[np.ndarray, np.ndarray]] = {}
        self._missing_pieces: set[int] = set()

    def _piece(self, key_int: int) -> tuple[np.ndarray, np.ndarray]:
        """Return (row indices, depth) of one simulated path piece."""

        cached = self._piece_cache.get(key_int)
        if cached is not None:
            return cached
        stat_path = os.path.join(self._depth_folder, f"{key_int}.win.stat.gz")
        if not os.path.isfile(stat_path):
            self._missing_pieces.add(key_int)
            empty = (np.empty(0, dtype=np.int64), np.empty(0, dtype=np.float64))
            self._piece_cache[key_int] = empty
            return empty
        frame = pd.read_csv(
            stat_path, compression="gzip", comment="#", sep="\t",
            names=list(WIN_STAT_COLUMNS),
        )
        rows = []
        values = []
        for chrom, start, depth in zip(frame["chr"], frame["st"], frame["meandepth"]):
            row = self._bin_index.row_by_key.get((str(chrom), int(start)))
            if row is None or depth == 0:
                continue
            rows.append(row)
            values.append(float(depth))
        piece = (
            np.asarray(rows, dtype=np.int64),
            np.asarray(values, dtype=np.float64),
        )
        self._piece_cache[key_int] = piece
        return piece

    def unit_profile(self, path: str) -> np.ndarray:
        """Per-window depth this path carries at one fitted copy.

        Values are traversal multiplicities: 1.0 where the path crosses the
        window once, 2.0 where it doubles back over it.
        """

        key_ints = self._path_key_ints.get(path)
        if key_ints is None:
            raise KeyError(f"Path has no separated depth pieces: {path}")
        profile = np.zeros(self._bin_index.size, dtype=np.float64)
        for key_int in key_ints:
            rows, values = self._piece(int(key_int))
            if rows.size:
                np.add.at(profile, rows, values)
        return profile

    @property
    def missing_pieces(self) -> tuple[int, ...]:
        return tuple(sorted(self._missing_pieces))


def build_cluster_depth_tracks(
    *,
    members: Sequence[Hashable],
    path_metadata: Mapping,
    path_spans: Mapping,
    bin_index: DepthBinIndex,
    observed_cn: np.ndarray,
    model_cn: np.ndarray | None,
    profiles: PathDepthProfiles,
) -> ClusterDepthTracks:
    """Scale each member's ``A`` column by its fitted copy number."""

    member_cn: dict[Hashable, np.ndarray] = {}
    for path in members:
        metadata = path_metadata.get(path, {})
        location = str(metadata.get("location", path))
        depth_n = float(metadata.get("depth_n", 0.0))
        member_cn[path] = profiles.unit_profile(location) * depth_n

    chrom_row_sets: dict[str, set[int]] = {}
    cluster_row_set: set[int] = set()
    for path in members:
        for span in path_spans[path]:
            rows = bin_index.rows_for_span(span.chrom, span.start, span.end)
            cluster_row_set.update(rows.tolist())
            chrom_row_sets.setdefault(span.chrom, set()).update(rows.tolist())
    return ClusterDepthTracks(
        bin_index=bin_index,
        observed_cn=observed_cn,
        model_cn=model_cn,
        member_cn=member_cn,
        cluster_rows=np.asarray(sorted(cluster_row_set), dtype=np.int64),
        chrom_rows={
            chrom: np.asarray(sorted(rows), dtype=np.int64)
            for chrom, rows in chrom_row_sets.items()
        },
    )


def summarize_region(
    *,
    label: str,
    chrom: str,
    start: int,
    end: int,
    tracks: ClusterDepthTracks,
    carriers: Sequence[Hashable] = (),
) -> RegionDepthSummary:
    """Length-weighted depth means over one reference interval."""

    bin_index = tracks.bin_index
    rows = bin_index.rows_for_span(chrom, start, end)
    weights = bin_index.overlap_lengths(rows, start, end)
    total_weight = float(weights.sum())

    def weighted_mean(values: np.ndarray) -> float:
        if total_weight <= 0:
            return 0.0
        return float(np.dot(values[rows], weights) / total_weight)

    observed_mean = weighted_mean(tracks.observed_cn)
    if total_weight > 0:
        deviation = tracks.observed_cn[rows] - observed_mean
        observed_sd = float(
            np.sqrt(max(np.dot(deviation ** 2, weights) / total_weight, 0.0))
        )
    else:
        observed_sd = 0.0
    return RegionDepthSummary(
        label=label,
        chrom=chrom,
        start=int(start),
        end=int(end),
        observed_cn=observed_mean,
        observed_sd=observed_sd,
        model_cn=None if tracks.model_cn is None else weighted_mean(tracks.model_cn),
        member_cn={
            path: weighted_mean(values) for path, values in tracks.member_cn.items()
        },
        carriers=tuple(carriers),
    )


def summarize_path_span(
    *,
    path: Hashable,
    spans: Iterable,
    tracks: ClusterDepthTracks,
) -> RegionDepthSummary:
    """Depth attribution over the merged reference footprint of one path."""

    bin_index = tracks.bin_index
    row_set: set[int] = set()
    chroms: set[str] = set()
    for span in spans:
        row_set.update(bin_index.rows_for_span(span.chrom, span.start, span.end).tolist())
        chroms.add(span.chrom)
    rows = np.asarray(sorted(row_set), dtype=np.int64)
    weights = (bin_index.ends[rows] - bin_index.starts[rows]).astype(np.float64)
    total_weight = float(weights.sum())

    def weighted_mean(values: np.ndarray) -> float:
        if total_weight <= 0:
            return 0.0
        return float(np.dot(values[rows], weights) / total_weight)

    observed_mean = weighted_mean(tracks.observed_cn)
    return RegionDepthSummary(
        label=str(path),
        chrom="+".join(sorted(chroms)),
        start=0,
        end=int(total_weight),
        observed_cn=observed_mean,
        observed_sd=0.0,
        model_cn=None if tracks.model_cn is None else weighted_mean(tracks.model_cn),
        member_cn={
            member: weighted_mean(values) for member, values in tracks.member_cn.items()
        },
        carriers=(path,),
    )
