"""Per-cluster renderer for shared-reference SKYPE path components.

Four layouts are available through ``modes``:

``plain``        the original branch graph and contribution legend
``lane``         the branch graph with a reference-depth ribbon inside each
                 chromosome lane, so a route and the coverage it must explain
                 share one x position
``profile``      the branch graph above full copy-number axes, one per
                 chromosome, with observed / fitted / stacked member depth
``attribution``  region-level bars answering how much of the observed depth
                 each member contributes over every shared block
"""

from __future__ import annotations

import csv
import glob
import math
import os
from pathlib import Path
from typing import Hashable, Mapping, Sequence

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np

from reference_path_clustering import (
    ReferencePathCluster,
    ReferenceSpan,
    SharedReferenceBlock,
    build_shared_reference_blocks,
)
from reference_path_depth import (
    ClusterDepthTracks,
    RegionDepthSummary,
    summarize_path_span,
    summarize_region,
)
from virtual_sky_plotting import CHR_COLORS


M = 1_000_000
BLOCK_COLORS = (
    "#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00",
    "#56B4E9", "#7A5195", "#EF5675", "#2F4B7C", "#59A14F",
)
DEPTH_MODES = ("plain", "lane", "profile", "attribution")
MODE_STEMS = {
    "plain": "virtual_sky_cluster_refspan",
    "lane": "virtual_sky_cluster_refspan_lane",
    "profile": "virtual_sky_cluster_refspan_profile",
    "attribution": "virtual_sky_cluster_refspan_attrib",
}
OBSERVED_COLOR = "#9AA3AE"
MODEL_COLOR = "#1F2933"
OTHER_COLOR = "#CBD2D9"
LANE_LEFT = 0.14
LANE_RIGHT = 0.97


def _short_path(path: str) -> str:
    parts = Path(path).parts
    return "/".join(parts[-2:]) if len(parts) >= 2 else path


def _format_coordinate(value: int) -> str:
    return f"{value / M:.3f}"


def _path_color(index: int):
    return plt.get_cmap("tab10")(index % 10)


def _path_label(path: Hashable) -> str:
    return f"P{int(path):04d}"


def _metadata_value(path_metadata: Mapping, path, key, default=None):
    return path_metadata.get(path, {}).get(key, default)


def _chromosome_order(cluster: ReferencePathCluster, path_spans: Mapping) -> list[str]:
    """Keep chromosomes close to their first traversal in the component."""

    order = []
    seen = set()
    for path in cluster.members:
        for span in path_spans[path]:
            if span.chrom not in seen:
                order.append(span.chrom)
                seen.add(span.chrom)
    return order


def _chromosome_bounds(cluster: ReferencePathCluster, path_spans: Mapping) -> dict:
    bounds = {}
    for path in cluster.members:
        for span in path_spans[path]:
            if span.chrom not in bounds:
                bounds[span.chrom] = [span.start, span.end]
            else:
                bounds[span.chrom][0] = min(bounds[span.chrom][0], span.start)
                bounds[span.chrom][1] = max(bounds[span.chrom][1], span.end)
    return {chrom: tuple(values) for chrom, values in bounds.items()}


def _coordinate_mapper(bounds: Mapping, left=LANE_LEFT, right=LANE_RIGHT):
    def map_coordinate(chrom: str, coordinate) -> float:
        start, end = bounds[chrom]
        if end <= start:
            return (left + right) / 2
        return left + (np.asarray(coordinate) - start) / (end - start) * (right - left)
    return map_coordinate


def _nice_ceiling(value: float) -> float:
    """Round a copy-number axis limit up to a readable tick value."""

    if not math.isfinite(value) or value <= 0:
        return 2.0
    if value <= 10:
        return float(max(2, math.ceil(value)))
    if value <= 40:
        return float(2 * math.ceil(value / 2))
    return float(5 * math.ceil(value / 5))


def _depth_ceiling(
    tracks: ClusterDepthTracks,
    chromosome_order: Sequence[str],
    bounds: Mapping,
) -> float:
    """One copy-number scale for every lane, so bands stay comparable."""

    observed_peaks = []
    hard_peaks = [2.0]
    cluster_cn = tracks.cluster_cn
    for chrom in chromosome_order:
        rows = tracks.bin_index.rows_for_span(chrom, *bounds[chrom])
        if rows.size == 0:
            continue
        observed_peaks.append(np.percentile(tracks.observed_cn[rows], 99.0))
        hard_peaks.append(float(cluster_cn[rows].max()))
        if tracks.model_cn is not None:
            hard_peaks.append(float(np.percentile(tracks.model_cn[rows], 99.0)))
    peak = max(hard_peaks + observed_peaks) if observed_peaks else max(hard_peaks)
    return _nice_ceiling(peak * 1.08)


def _step_series(
    tracks: ClusterDepthTracks, rows: np.ndarray, values: Sequence[np.ndarray]
) -> tuple[np.ndarray, list[np.ndarray]]:
    """Window edges plus ``step='post'`` aligned copies of each track."""

    starts = tracks.bin_index.starts[rows].astype(float)
    edges = np.append(starts, float(tracks.bin_index.ends[rows][-1]))
    return edges, [np.append(value[rows], value[rows][-1]) for value in values]


def _member_stack(
    tracks: ClusterDepthTracks, members: Sequence[Hashable], rows: np.ndarray
) -> list[np.ndarray]:
    return [tracks.member_cn[path] for path in members]


def _axes_heading(ax, title: str, subtitle: str) -> None:
    """Title plus caption at fixed point offsets, immune to axes stretching."""

    ax.set_title(title, loc="left", fontsize=14, fontweight="bold", pad=24)
    ax.annotate(
        subtitle, xy=(0.0, 1.0), xycoords="axes fraction",
        xytext=(0, 7), textcoords="offset points",
        ha="left", va="bottom", fontsize=9.2, color="#555555",
        annotation_clip=False,
    )


# --------------------------------------------------------------------------- #
# Branch graph
# --------------------------------------------------------------------------- #


def _draw_chromosome_route_graph(
    ax,
    cluster: ReferencePathCluster,
    path_spans: Mapping,
    shared_blocks: Sequence[SharedReferenceBlock],
    path_metadata: Mapping,
    *,
    tracks: ClusterDepthTracks | None = None,
    lane_spacing: float = 1.0,
    compact: bool = False,
) -> None:
    """Draw depth-weighted paths branching across chromosome-coordinate lanes."""

    members = list(cluster.members)
    chromosome_order = _chromosome_order(cluster, path_spans)
    bounds = _chromosome_bounds(cluster, path_spans)
    map_coordinate = _coordinate_mapper(bounds)
    chromosome_y = {
        chrom: (len(chromosome_order) - 1 - index) * lane_spacing
        for index, chrom in enumerate(chromosome_order)
    }
    max_depth_n = max(
        (float(_metadata_value(path_metadata, path, "depth_n", 0.0)) for path in members),
        default=1.0,
    ) or 1.0
    if len(members) == 1:
        path_offset = {members[0]: 0.0}
    else:
        path_offset = {
            path: (index - (len(members) - 1) / 2) * min(0.055, 0.22 / len(members))
            for index, path in enumerate(members)
        }

    ribbon_top = -0.18
    ribbon_bottom = ribbon_top - (0.58 * lane_spacing if tracks is not None else 0.0)
    axis_bottom = ribbon_bottom - 0.34 if tracks is not None else -0.72
    depth_ceiling = None
    if tracks is not None:
        depth_ceiling = _depth_ceiling(tracks, chromosome_order, bounds)

    ax.set_xlim(0, 1)
    ax.set_ylim(
        axis_bottom - 0.28,
        (len(chromosome_order) - 1) * lane_spacing + 0.92,
    )
    ax.axis("off")
    if not compact:
        _axes_heading(
            ax, "Chromosome-coordinate branch graph",
            "Each chromosome has its own x scale; curves are ordered path junctions.",
        )

    for chrom in chromosome_order:
        y = chromosome_y[chrom]
        chrom_color = CHR_COLORS.get(chrom, "#A0A0A0")
        ax.add_patch(patches.FancyBboxPatch(
            (0.012, y - 0.115), 0.096, 0.23,
            boxstyle="round,pad=0.005,rounding_size=0.018",
            facecolor="white", edgecolor=chrom_color, linewidth=2.2,
            zorder=8,
        ))
        ax.text(
            0.060, y, chrom, ha="center", va="center",
            fontsize=9.5, fontweight="bold", zorder=9,
        )
        ax.plot(
            [LANE_LEFT, LANE_RIGHT], [y, y], color=chrom_color,
            linewidth=7.0, alpha=0.16, solid_capstyle="round", zorder=0,
        )
        if tracks is not None:
            _draw_lane_depth_ribbon(
                ax, tracks, members, chrom, bounds, map_coordinate,
                baseline=y + ribbon_bottom, height=ribbon_top - ribbon_bottom,
                depth_ceiling=depth_ceiling,
            )
        start, end = bounds[chrom]
        midpoint = (start + end) // 2
        tick_top = y + (ribbon_bottom if tracks is not None else -0.09)
        for coordinate, x, alignment in (
            (start, LANE_LEFT, "left"),
            (midpoint, (LANE_LEFT + LANE_RIGHT) / 2, "center"),
            (end, LANE_RIGHT, "right"),
        ):
            if tracks is None:
                ax.plot([x, x], [y - 0.09, y + 0.09], color="#777777", linewidth=0.6)
            ax.text(
                x, tick_top - 0.06, f"{coordinate / M:.1f} Mb",
                ha=alignment, va="top", fontsize=7.0, color="#666666",
            )

    # Shared reference is the common trunk; individual path lanes sit on top.
    label_levels: dict[str, list[float]] = {}
    for block_index, block in enumerate(shared_blocks):
        y = chromosome_y[block.chrom]
        x_start = float(map_coordinate(block.chrom, block.start))
        x_end = float(map_coordinate(block.chrom, block.end))
        x_center = (x_start + x_end) / 2
        color = BLOCK_COLORS[block_index % len(BLOCK_COLORS)]
        chrom_levels = label_levels.setdefault(block.chrom, [])
        label_level = next(
            (
                level
                for level, previous_x in enumerate(chrom_levels)
                if x_center - previous_x >= 0.055
            ),
            len(chrom_levels),
        )
        if label_level == len(chrom_levels):
            chrom_levels.append(x_center)
        else:
            chrom_levels[label_level] = x_center
        ax.plot(
            [x_start, x_end], [y, y], color="white",
            linewidth=15, solid_capstyle="butt", zorder=1,
        )
        ax.plot(
            [x_start, x_end], [y, y], color=color,
            linewidth=10, alpha=0.34, solid_capstyle="butt", zorder=2,
        )
        ax.text(
            x_center, y + 0.17 + 0.085 * label_level,
            block.block_id, ha="center", va="bottom",
            fontsize=8.0, color=color, fontweight="bold", zorder=9,
        )

    # Junction curves are drawn first so the chromosome traversals read as nodes.
    for path_index, path in enumerate(members):
        spans = tuple(path_spans[path])
        depth_n = float(_metadata_value(path_metadata, path, "depth_n", 0.0))
        linewidth = 1.4 + 4.2 * depth_n / max_depth_n
        color = _path_color(path_index)
        offset = path_offset[path]
        for junction_index, (left_span, right_span) in enumerate(zip(spans, spans[1:])):
            if (
                left_span.source == "interrupt"
                and right_span.source == "interrupt"
                and left_span.chrom == right_span.chrom
            ):
                continue
            left_exit_coordinate = left_span.end if left_span.strand == "+" else left_span.start
            right_entry_coordinate = right_span.start if right_span.strand == "+" else right_span.end
            left_point = (
                float(map_coordinate(left_span.chrom, left_exit_coordinate)),
                chromosome_y[left_span.chrom] + offset,
            )
            right_point = (
                float(map_coordinate(right_span.chrom, right_entry_coordinate)),
                chromosome_y[right_span.chrom] + offset,
            )
            if left_span.chrom == right_span.chrom:
                curvature = 0.28 + 0.05 * (path_index % 3)
            else:
                direction = -1 if (path_index + junction_index) % 2 else 1
                curvature = direction * (0.08 + 0.025 * path_index)
            ax.add_patch(patches.FancyArrowPatch(
                left_point, right_point,
                arrowstyle="-|>", mutation_scale=7 + linewidth,
                connectionstyle=f"arc3,rad={curvature}",
                color=color, linewidth=linewidth, alpha=0.72,
                shrinkA=1.5, shrinkB=1.5, zorder=3,
            ))
            ax.scatter(
                [left_point[0], right_point[0]], [left_point[1], right_point[1]],
                s=11 + linewidth * 2, facecolor=color, edgecolor="white",
                linewidth=0.5, zorder=7,
            )

        for span_index, span in enumerate(spans):
            y = chromosome_y[span.chrom] + offset
            entry_coordinate = span.start if span.strand == "+" else span.end
            exit_coordinate = span.end if span.strand == "+" else span.start
            entry_x = float(map_coordinate(span.chrom, entry_coordinate))
            exit_x = float(map_coordinate(span.chrom, exit_coordinate))
            next_span = spans[span_index + 1] if span_index + 1 < len(spans) else None
            suppress_interrupt_arrow = (
                span.source == "interrupt"
                and next_span is not None
                and next_span.source == "interrupt"
                and next_span.chrom == span.chrom
            )
            ax.add_patch(patches.FancyArrowPatch(
                (entry_x, y), (exit_x, y),
                arrowstyle="-" if suppress_interrupt_arrow else "-|>",
                mutation_scale=6 + linewidth,
                color=color, linewidth=linewidth, alpha=0.78,
                shrinkA=0, shrinkB=0, zorder=5,
            ))
        if spans:
            first = spans[0]
            last = spans[-1]
            first_coordinate = first.start if first.strand == "+" else first.end
            last_coordinate = last.end if last.strand == "+" else last.start
            ax.scatter(
                [float(map_coordinate(first.chrom, first_coordinate))],
                [chromosome_y[first.chrom] + offset],
                marker="o", s=32, facecolor="white", edgecolor=color,
                linewidth=1.6, zorder=10,
            )
            ax.scatter(
                [float(map_coordinate(last.chrom, last_coordinate))],
                [chromosome_y[last.chrom] + offset],
                marker="s", s=29, facecolor=color, edgecolor="white",
                linewidth=0.7, zorder=10,
            )

    if not compact:
        caption = "Chromosome-specific coordinates (independent Mb scales)"
        if tracks is not None:
            caption += (
                f"  ·  depth ribbon 0–{depth_ceiling:.0f}N: grey = observed reference "
                "depth, colours = stacked member contribution, dashed = all fitted paths"
            )
        ax.text(
            0.555, axis_bottom - 0.14, caption,
            ha="center", va="center", fontsize=9, color="#444444",
        )


def _draw_lane_depth_ribbon(
    ax,
    tracks: ClusterDepthTracks,
    members: Sequence[Hashable],
    chrom: str,
    bounds: Mapping,
    map_coordinate,
    *,
    baseline: float,
    height: float,
    depth_ceiling: float,
) -> None:
    """Fill one lane's depth strip with observed, fitted, and member depth."""

    start, end = bounds[chrom]
    rows = tracks.bin_index.rows_for_span(chrom, start, end)
    if rows.size == 0:
        return

    def to_y(values):
        return baseline + np.clip(np.asarray(values) / depth_ceiling, 0, 1) * height

    stack_tracks = _member_stack(tracks, members, rows)
    model = tracks.model_cn
    series = [tracks.observed_cn] + stack_tracks + ([model] if model is not None else [])
    edges, stepped = _step_series(tracks, rows, series)
    x = np.clip(map_coordinate(chrom, edges), LANE_LEFT, LANE_RIGHT)
    observed = stepped[0]
    member_values = stepped[1:1 + len(stack_tracks)]
    model_values = stepped[-1] if model is not None else None

    ax.plot(
        [LANE_LEFT, LANE_RIGHT], [baseline, baseline],
        color="#B7BDC6", linewidth=0.7, zorder=1,
    )
    ax.fill_between(
        x, baseline, to_y(observed), step="post",
        facecolor=OBSERVED_COLOR, alpha=0.55, linewidth=0, zorder=2,
    )
    ax.step(
        x, to_y(observed), where="post",
        color="#6B7280", linewidth=0.7, alpha=0.9, zorder=4,
    )
    cumulative = np.zeros_like(observed)
    for member_index, values in enumerate(member_values):
        upper = cumulative + values
        ax.fill_between(
            x, to_y(cumulative), to_y(upper), step="post",
            facecolor=_path_color(member_index), alpha=0.72,
            linewidth=0, zorder=3,
        )
        cumulative = upper
    if model_values is not None:
        ax.step(
            x, to_y(model_values), where="post",
            color=MODEL_COLOR, linewidth=0.9, linestyle=(0, (3, 1.6)),
            alpha=0.85, zorder=5,
        )
    for level in (2.0, depth_ceiling):
        if level > depth_ceiling:
            continue
        ax.plot(
            [LANE_LEFT, LANE_RIGHT], [to_y(level)] * 2,
            color="#9CA3AF", linewidth=0.5, linestyle=(0, (1, 3)), zorder=1,
        )
        ax.text(
            LANE_LEFT - 0.006, float(to_y(level)), f"{level:.0f}N",
            ha="right", va="center", fontsize=6.2, color="#8A8F98",
        )


# --------------------------------------------------------------------------- #
# Copy-number profile panels
# --------------------------------------------------------------------------- #


def _draw_depth_profile_panel(
    ax,
    tracks: ClusterDepthTracks,
    members: Sequence[Hashable],
    chrom: str,
    bounds: Mapping,
    shared_blocks: Sequence[SharedReferenceBlock],
    *,
    depth_ceiling: float,
    show_legend: bool,
) -> None:
    """One chromosome's observed / fitted / stacked-member copy number."""

    start, end = bounds[chrom]
    rows = tracks.bin_index.rows_for_span(chrom, start, end)
    chrom_color = CHR_COLORS.get(chrom, "#A0A0A0")
    ax.set_xlim(start / M, end / M)
    ax.set_ylim(0, depth_ceiling)
    ax.set_ylabel("CN (N)", fontsize=8)
    ax.tick_params(axis="both", labelsize=7.4)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#B7BDC6")
    ax.spines["bottom"].set_color("#B7BDC6")
    ax.set_yticks(np.linspace(0, depth_ceiling, min(5, int(depth_ceiling) + 1)))
    ax.axhline(2.0, color="#9CA3AF", linewidth=0.6, linestyle=(0, (1, 3)), zorder=1)
    ax.text(
        0.004, 0.94, chrom, transform=ax.transAxes, ha="left", va="top",
        fontsize=9.5, fontweight="bold", color="#1F2933",
        bbox={
            "boxstyle": "round,pad=0.24", "facecolor": "white",
            "edgecolor": chrom_color, "linewidth": 1.6,
        },
        zorder=12,
    )

    for block_index, block in enumerate(shared_blocks):
        if block.chrom != chrom:
            continue
        color = BLOCK_COLORS[block_index % len(BLOCK_COLORS)]
        ax.axvspan(
            block.start / M, block.end / M,
            facecolor=color, alpha=0.10, linewidth=0, zorder=0,
        )
        ax.text(
            (block.start + block.end) / 2 / M, depth_ceiling * 0.985,
            block.block_id, ha="center", va="top", fontsize=7.2,
            color=color, fontweight="bold", zorder=11,
        )

    if rows.size == 0:
        return

    stack_tracks = _member_stack(tracks, members, rows)
    model = tracks.model_cn
    series = [tracks.observed_cn] + stack_tracks + ([model] if model is not None else [])
    edges, stepped = _step_series(tracks, rows, series)
    x = edges / M
    observed = stepped[0]
    member_values = stepped[1:1 + len(stack_tracks)]

    ax.fill_between(
        x, 0, observed, step="post", facecolor=OBSERVED_COLOR, alpha=0.5,
        linewidth=0, zorder=2, label="observed reference depth",
    )
    ax.step(
        x, observed, where="post", color="#6B7280", linewidth=0.8,
        alpha=0.9, zorder=4,
    )
    cumulative = np.zeros_like(observed)
    for member_index, (path, values) in enumerate(zip(members, member_values)):
        upper = cumulative + values
        ax.fill_between(
            x, cumulative, upper, step="post",
            facecolor=_path_color(member_index), alpha=0.75, linewidth=0,
            zorder=3, label=f"{_path_label(path)} contribution",
        )
        cumulative = upper
    ax.step(
        x, cumulative, where="post", color="#3F4A56", linewidth=0.8,
        alpha=0.8, zorder=5,
    )
    if model is not None:
        ax.step(
            x, stepped[-1], where="post", color=MODEL_COLOR, linewidth=1.0,
            linestyle=(0, (3, 1.6)), alpha=0.9, zorder=6,
            label="all fitted paths (A·w)",
        )
    if show_legend:
        ax.legend(
            loc="upper right", fontsize=6.8, framealpha=0.92, ncol=2,
            borderpad=0.35, handlelength=1.5, columnspacing=0.9,
        )


# --------------------------------------------------------------------------- #
# Region attribution
# --------------------------------------------------------------------------- #


def _cluster_region_summaries(
    cluster: ReferencePathCluster,
    path_spans: Mapping,
    shared_blocks: Sequence[SharedReferenceBlock],
    tracks: ClusterDepthTracks,
) -> tuple[list[RegionDepthSummary], list[RegionDepthSummary]]:
    """Depth attribution per shared block and per whole member footprint."""

    block_summaries = [
        summarize_region(
            label=block.block_id, chrom=block.chrom,
            start=block.start, end=block.end,
            tracks=tracks, carriers=block.carriers,
        )
        for block in shared_blocks
    ]
    path_summaries = [
        summarize_path_span(path=path, spans=path_spans[path], tracks=tracks)
        for path in cluster.members
    ]
    return block_summaries, path_summaries


def _draw_region_attribution(
    ax,
    cluster: ReferencePathCluster,
    shared_blocks: Sequence[SharedReferenceBlock],
    block_summaries: Sequence[RegionDepthSummary],
    path_summaries: Sequence[RegionDepthSummary],
    tracks: ClusterDepthTracks,
) -> None:
    """Stacked member depth against the observed depth of each region."""

    members = list(cluster.members)
    rows: list[tuple[str, str, RegionDepthSummary]] = []
    for block, summary in zip(shared_blocks, block_summaries):
        rows.append((
            "block",
            f"{block.block_id}  {block.chrom}:"
            f"{block.start / M:.1f}–{block.end / M:.1f} Mb",
            summary,
        ))
    for path, summary in zip(members, path_summaries):
        rows.append((
            "path",
            f"{_path_label(path)} whole footprint  ({summary.end / M:.1f} Mb)",
            summary,
        ))

    _axes_heading(
        ax, "Depth attribution per region",
        "Bars stack each member's fitted contribution; hatched grey is every other "
        "fitted path; ◆ marks the observed reference depth (±1 s.d. across windows).",
    )
    if not rows:
        ax.axis("off")
        ax.text(0.5, 0.5, "No shared reference block in this cluster",
                ha="center", va="center", fontsize=10, color="#888888")
        return

    y_positions = np.arange(len(rows))[::-1].astype(float)
    max_value = 0.0
    for _kind, _label, summary in rows:
        candidates = [summary.cluster_cn, summary.observed_cn + summary.observed_sd]
        if summary.model_cn is not None:
            candidates.append(summary.model_cn)
        max_value = max(max_value, *candidates)
    axis_max = _nice_ceiling(max_value * 1.18)

    for y, (kind, label, summary) in zip(y_positions, rows):
        height = 0.52 if kind == "block" else 0.40
        left = 0.0
        for member_index, path in enumerate(members):
            value = summary.member_cn.get(path, 0.0)
            if value <= 0:
                continue
            ax.barh(
                y, value, left=left, height=height,
                color=_path_color(member_index), alpha=0.85,
                edgecolor="white", linewidth=0.6, zorder=3,
            )
            if value >= axis_max * 0.055:
                ax.text(
                    left + value / 2, y, f"{value:.2f}",
                    ha="center", va="center", fontsize=7.0,
                    color="white", fontweight="bold", zorder=5,
                )
            left += value
        other = summary.other_cn
        if other is not None and other > 0.01:
            ax.barh(
                y, other, left=left, height=height,
                color=OTHER_COLOR, alpha=0.9, edgecolor="white",
                linewidth=0.6, hatch="///", zorder=3,
            )
            left += other
        if summary.observed_sd > 0:
            ax.plot(
                [
                    max(summary.observed_cn - summary.observed_sd, 0.0),
                    summary.observed_cn + summary.observed_sd,
                ],
                [y, y], color=MODEL_COLOR, linewidth=1.0, alpha=0.55, zorder=6,
            )
        ax.scatter(
            [summary.observed_cn], [y], marker="D", s=44,
            facecolor="white", edgecolor=MODEL_COLOR, linewidth=1.5, zorder=7,
        )
        explained = (
            summary.cluster_cn / summary.observed_cn * 100
            if summary.observed_cn > 0 else 0.0
        )
        ax.text(
            axis_max * 0.995, y,
            f"obs {summary.observed_cn:.2f}N · Σcluster {summary.cluster_cn:.2f}N"
            f" ({explained:.0f}%)",
            ha="right", va="center", fontsize=7.2, color="#3F4A56", zorder=8,
        )
        if kind == "path":
            ax.axhline(y + 0.5, color="#E4E7EB", linewidth=0.7, zorder=0)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [label for _kind, label, _summary in rows], fontsize=8.0,
    )
    for tick, (kind, _label, _summary) in zip(ax.get_yticklabels(), rows):
        tick.set_fontweight("bold" if kind == "block" else "normal")
        tick.set_color("#1F2933" if kind == "block" else "#5A6472")
    ax.set_xlim(0, axis_max)
    ax.set_xlabel("copy number (N)", fontsize=9)
    ax.tick_params(axis="x", labelsize=7.6)
    ax.set_ylim(-0.8, len(rows) - 0.2)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#B7BDC6")
    ax.axvline(2.0, color="#9CA3AF", linewidth=0.6, linestyle=(0, (1, 3)), zorder=1)
    ax.grid(axis="x", color="#EDEFF2", linewidth=0.6, zorder=0)
    ax.set_axisbelow(True)


def _draw_observed_model_scatter(
    ax, tracks: ClusterDepthTracks, shared_blocks: Sequence[SharedReferenceBlock]
) -> None:
    """Per-window agreement over the reference footprint of this cluster."""

    ax.set_title(
        "Observed vs fitted depth (cluster windows)",
        loc="left", fontsize=10.5, fontweight="bold", pad=8,
    )
    rows = tracks.cluster_rows
    if tracks.model_cn is None or rows.size == 0:
        ax.axis("off")
        ax.text(
            0.5, 0.5, "predict_B unavailable", ha="center", va="center",
            fontsize=9, color="#888888",
        )
        return

    shared_rows: set[int] = set()
    for block in shared_blocks:
        shared_rows.update(
            tracks.bin_index.rows_for_span(block.chrom, block.start, block.end).tolist()
        )
    in_block = np.asarray([int(row) in shared_rows for row in rows])
    observed = tracks.observed_cn[rows]
    model = tracks.model_cn[rows]
    # A handful of CENSAT windows carry extreme coverage; clip so the bulk stays legible.
    limit = _nice_ceiling(
        max(np.percentile(observed, 99.5), model.max()) * 1.08
    )

    ax.scatter(
        observed[~in_block], model[~in_block], s=5, alpha=0.30,
        facecolor="#8A94A6", edgecolor="none", label="outside shared blocks",
    )
    ax.scatter(
        observed[in_block], model[in_block], s=7, alpha=0.55,
        facecolor="#D55E00", edgecolor="none", label="inside shared blocks",
    )
    ax.plot([0, limit], [0, limit], color=MODEL_COLOR, linewidth=0.9,
            linestyle=(0, (3, 1.6)), alpha=0.8)
    residual = float(np.sqrt(np.mean((model - observed) ** 2)))
    clipped = int((observed > limit).sum())
    summary = f"RMSE {residual:.2f}N · {rows.size} windows"
    if clipped:
        summary += f"\n{clipped} window(s) above {limit:.0f}N not shown"
    ax.text(
        0.03, 0.96, summary,
        transform=ax.transAxes, ha="left", va="top", fontsize=7.8,
        color="#3F4A56",
    )
    ax.set_xlim(0, limit)
    ax.set_ylim(0, limit)
    ax.set_xlabel("observed CN (N)", fontsize=8.4)
    ax.set_ylabel("fitted CN (N)", fontsize=8.4)
    ax.tick_params(labelsize=7.2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.legend(loc="lower right", fontsize=7.0, framealpha=0.9, borderpad=0.3)


# --------------------------------------------------------------------------- #
# Legend
# --------------------------------------------------------------------------- #


def _draw_cluster_legend(
    ax,
    cluster: ReferencePathCluster,
    shared_blocks: Sequence[SharedReferenceBlock],
    path_metadata: Mapping,
    *,
    block_summaries: Sequence[RegionDepthSummary] = (),
    path_summaries: Sequence[RegionDepthSummary] = (),
    show_depth_lines: bool = True,
) -> None:
    """List path depths and exact shared spans beside the route graph."""

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    members = list(cluster.members)
    max_depth_n = max(
        (float(_metadata_value(path_metadata, path, "depth_n", 0.0)) for path in members),
        default=1.0,
    ) or 1.0
    has_depth = bool(block_summaries or path_summaries) and show_depth_lines
    line_count = 3 if has_depth else 2
    item_count = len(members) + len(shared_blocks)
    step = min(0.094 * line_count / 2, 0.76 / max(item_count, 1))
    path_summary_by_path = {
        path: summary for path, summary in zip(members, path_summaries)
    }
    block_summary_by_id = {
        summary.label: summary for summary in block_summaries
    }

    ax.set_title("Path contributions", loc="left", fontsize=13, fontweight="bold", pad=14)
    y = 0.94
    for path_index, path in enumerate(members):
        depth_n = float(_metadata_value(path_metadata, path, "depth_n", 0.0))
        raw_depth = float(_metadata_value(path_metadata, path, "raw_depth", 0.0))
        location = str(_metadata_value(path_metadata, path, "location", path))
        color = _path_color(path_index)
        linewidth = 1.4 + 4.2 * depth_n / max_depth_n
        ax.plot(
            [0.02, 0.18], [y, y], color=color, linewidth=linewidth,
            solid_capstyle="round",
        )
        ax.scatter([0.02], [y], marker="o", s=28, facecolor="white", edgecolor=color)
        ax.scatter([0.18], [y], marker="s", s=25, facecolor=color, edgecolor="white")
        ax.text(
            0.22, y + step * (0.24 if has_depth else 0.16),
            f"{_path_label(path)}  {_short_path(location)}",
            ha="left", va="center", fontsize=8.6, fontweight="bold",
        )
        ax.text(
            0.22, y - step * (0.02 if has_depth else 0.22),
            f"{depth_n:.2f}N contribution · {raw_depth:.2f}x fitted",
            ha="left", va="center", fontsize=7.7, color="#555555",
        )
        summary = path_summary_by_path.get(path) if has_depth else None
        if summary is not None:
            share = (
                summary.member_cn.get(path, 0.0) / summary.observed_cn * 100
                if summary.observed_cn > 0 else 0.0
            )
            ax.text(
                0.22, y - step * 0.28,
                f"over own span: observed {summary.observed_cn:.2f}N · "
                f"this path {share:.0f}%",
                ha="left", va="center", fontsize=7.2, color="#7B8794",
            )
        y -= step

    y -= step * 0.20
    ax.text(
        0.02, y, "Shared reference spans", ha="left", va="center",
        fontsize=11, fontweight="bold",
    )
    y -= step * 0.65
    for block_index, block in enumerate(shared_blocks):
        color = BLOCK_COLORS[block_index % len(BLOCK_COLORS)]
        total_depth = sum(
            float(_metadata_value(path_metadata, path, "depth_n", 0.0))
            for path in block.carriers
        )
        carrier_text = ", ".join(_path_label(path) for path in block.carriers)
        ax.add_patch(patches.Rectangle(
            (0.02, y - 0.012), 0.055, 0.024,
            facecolor=color, edgecolor="none",
        ))
        ax.text(
            0.10, y + step * (0.24 if has_depth else 0.16),
            (
                f"{block.block_id}  {block.chrom}:"
                f"{_format_coordinate(block.start)}–{_format_coordinate(block.end)} Mb"
            ),
            ha="left", va="center", fontsize=8.2, fontweight="bold", color=color,
        )
        ax.text(
            0.10, y - step * (0.02 if has_depth else 0.22),
            f"{block.length / M:.2f} Mb · {carrier_text} · Σ {total_depth:.2f}N",
            ha="left", va="center", fontsize=7.3, color="#555555",
        )
        summary = block_summary_by_id.get(block.block_id) if has_depth else None
        if summary is not None:
            model_text = (
                "" if summary.model_cn is None
                else f" · all fitted {summary.model_cn:.2f}N"
            )
            ax.text(
                0.10, y - step * 0.28,
                f"observed {summary.observed_cn:.2f}N{model_text}",
                ha="left", va="center", fontsize=7.2, color="#7B8794",
            )
        y -= step

    footer = "○ path start    ■ path end    arrow = traversal direction"
    if has_depth:
        footer += "\ngrey fill = observed reference depth    dashed = all fitted paths"
    ax.text(
        0.02, max(0.025, y - step * 0.10), footer,
        ha="left", va="bottom", fontsize=7.6, color="#555555",
    )


# --------------------------------------------------------------------------- #
# Figure assembly
# --------------------------------------------------------------------------- #


def _finish_figure(
    fig, cell_line: str, cluster: ReferencePathCluster,
    minimum_overlap: int, full_reference_ratio: float, subtitle_extra: str,
) -> None:
    fig.suptitle(
        f"{cell_line} reference-sharing path cluster {cluster.cluster_id}",
        fontsize=18, fontweight="bold", y=0.985,
    )
    fig.text(
        0.5, 0.940,
        (
            f"same-reference overlap ≥ {minimum_overlap / 1000:.0f} kb · "
            f"no-NClose full reference ≥ {full_reference_ratio:.2f} excluded"
            f"{subtitle_extra}"
        ),
        ha="center", va="top", fontsize=10.5, color="#555555",
    )
    fig.text(
        0.5, 0.018,
        (
            "Interpretation: shared chromosome lanes split into alternative reconstructed routes; "
            "this is an overlap/ambiguity component, not haplotype phasing."
        ),
        ha="center", va="bottom", fontsize=10, color="#6F42C1", fontweight="bold",
    )


def _save_figure(fig, output_prefix: str, stem: str) -> tuple[str, str]:
    png_path = os.path.join(output_prefix, f"{stem}.png")
    pdf_path = os.path.join(output_prefix, f"{stem}.pdf")
    fig.savefig(png_path, dpi=180, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def render_reference_path_cluster(
    *,
    output_prefix: str,
    cell_line: str,
    cluster: ReferencePathCluster,
    path_spans: Mapping,
    path_metadata: Mapping,
    minimum_overlap: int,
    full_reference_ratio: float,
    mode: str = "plain",
    tracks: ClusterDepthTracks | None = None,
) -> dict:
    """Render one reference-sharing ambiguity cluster as PNG and PDF."""

    if mode not in DEPTH_MODES:
        raise ValueError(f"Unknown reference cluster figure mode: {mode}")
    if mode != "plain" and tracks is None:
        raise ValueError(f"Figure mode {mode!r} needs depth tracks")

    shared_blocks = build_shared_reference_blocks(cluster, path_spans)
    chromosome_order = _chromosome_order(cluster, path_spans)
    chromosome_count = len(chromosome_order)
    member_count = len(cluster.members)
    block_summaries: Sequence[RegionDepthSummary] = ()
    path_summaries: Sequence[RegionDepthSummary] = ()
    if tracks is not None:
        block_summaries, path_summaries = _cluster_region_summaries(
            cluster, path_spans, shared_blocks, tracks
        )

    legend_height = 3.4 + 0.48 * (member_count + len(shared_blocks))
    if mode == "plain":
        figure_height = max(7.5, 3.4 + 1.15 * chromosome_count, legend_height)
        fig, axes = plt.subplots(
            1, 2, figsize=(19, figure_height),
            gridspec_kw={"width_ratios": [0.72, 0.28], "wspace": 0.10},
        )
        _draw_chromosome_route_graph(
            axes[0], cluster, path_spans, shared_blocks, path_metadata
        )
        _draw_cluster_legend(axes[1], cluster, shared_blocks, path_metadata)
        fig.subplots_adjust(top=0.84, bottom=0.10, left=0.04, right=0.98)
        subtitle_extra = ""

    elif mode == "lane":
        figure_height = max(9.0, 3.4 + 2.15 * chromosome_count, legend_height + 1.2)
        fig, axes = plt.subplots(
            1, 2, figsize=(19, figure_height),
            gridspec_kw={"width_ratios": [0.72, 0.28], "wspace": 0.10},
        )
        _draw_chromosome_route_graph(
            axes[0], cluster, path_spans, shared_blocks, path_metadata,
            tracks=tracks, lane_spacing=1.9,
        )
        _draw_cluster_legend(
            axes[1], cluster, shared_blocks, path_metadata,
            block_summaries=block_summaries, path_summaries=path_summaries,
        )
        fig.subplots_adjust(top=0.84, bottom=0.10, left=0.04, right=0.98)
        subtitle_extra = " · lane ribbons show reference depth under each route"

    elif mode == "profile":
        branch_share = max(2.4, 0.95 * chromosome_count)
        depth_share = 1.45 * chromosome_count
        figure_height = max(9.5, 3.0 + branch_share + depth_share, legend_height)
        fig = plt.figure(figsize=(19, figure_height))
        outer = fig.add_gridspec(
            2, 2, width_ratios=[0.72, 0.28],
            height_ratios=[branch_share, depth_share],
            wspace=0.10, hspace=0.24,
        )
        _draw_chromosome_route_graph(
            fig.add_subplot(outer[0, 0]), cluster, path_spans,
            shared_blocks, path_metadata,
        )
        _draw_cluster_legend(
            fig.add_subplot(outer[:, 1]), cluster, shared_blocks, path_metadata,
            block_summaries=block_summaries, path_summaries=path_summaries,
        )
        bounds = _chromosome_bounds(cluster, path_spans)
        depth_ceiling = _depth_ceiling(tracks, chromosome_order, bounds)
        depth_grid = outer[1, 0].subgridspec(chromosome_count, 1, hspace=0.40)
        for chrom_index, chrom in enumerate(chromosome_order):
            panel = fig.add_subplot(depth_grid[chrom_index])
            _draw_depth_profile_panel(
                panel, tracks, list(cluster.members), chrom, bounds, shared_blocks,
                depth_ceiling=depth_ceiling, show_legend=chrom_index == 0,
            )
            if chrom_index == 0:
                _axes_heading(
                    panel, "Reference depth explained by this cluster",
                    "Same Mb window as the lane above; grey is the cell line's own "
                    "coverage, colours stack each member's fitted contribution.",
                )
            if chrom_index == chromosome_count - 1:
                panel.set_xlabel("chromosome coordinate (Mb)", fontsize=8.6)
        fig.subplots_adjust(top=0.86, bottom=0.07, left=0.045, right=0.98)
        subtitle_extra = " · copy-number panels share each lane's Mb window"

    else:  # attribution
        row_count = member_count + len(shared_blocks)
        figure_height = max(9.0, 4.0 + 1.0 * chromosome_count + 0.46 * row_count)
        fig = plt.figure(figsize=(19, figure_height))
        outer = fig.add_gridspec(
            2, 2, width_ratios=[0.66, 0.34],
            height_ratios=[max(2.4, 1.0 * chromosome_count), max(3.0, 0.46 * row_count)],
            wspace=0.16, hspace=0.20,
        )
        _draw_chromosome_route_graph(
            fig.add_subplot(outer[0, 0]), cluster, path_spans,
            shared_blocks, path_metadata,
        )
        _draw_region_attribution(
            fig.add_subplot(outer[1, 0]), cluster, shared_blocks,
            block_summaries, path_summaries, tracks,
        )
        # Region bars already carry the observed/fitted numbers, so the legend
        # here stays two-line and fits the half-height right column.
        _draw_cluster_legend(
            fig.add_subplot(outer[0, 1]), cluster, shared_blocks, path_metadata,
            block_summaries=block_summaries, path_summaries=path_summaries,
            show_depth_lines=False,
        )
        _draw_observed_model_scatter(
            fig.add_subplot(outer[1, 1]), tracks, shared_blocks
        )
        fig.subplots_adjust(top=0.86, bottom=0.08, left=0.10, right=0.97)
        subtitle_extra = " · region bars quantify each member's share of observed depth"

    _finish_figure(
        fig, cell_line, cluster, minimum_overlap, full_reference_ratio, subtitle_extra
    )
    png_path, pdf_path = _save_figure(
        fig, output_prefix, f"{MODE_STEMS[mode]}_{cluster.cluster_id}"
    )
    return {
        "cluster_id": cluster.cluster_id,
        "mode": mode,
        "png_path": png_path,
        "pdf_path": pdf_path,
        "shared_block_count": len(shared_blocks),
    }


def _remove_stale_cluster_figures(output_prefix: str) -> None:
    """Remove only older files emitted under this renderer's exact prefixes."""

    for stem in MODE_STEMS.values():
        for extension in ("png", "pdf"):
            for path in glob.glob(os.path.join(output_prefix, f"{stem}_C*.{extension}")):
                os.unlink(path)


def write_reference_path_cluster_tables(
    *,
    output_prefix: str,
    clusters: Sequence[ReferencePathCluster],
    singletons: Sequence,
    path_metadata: Mapping,
    excluded_full_reference: Sequence[Mapping],
) -> dict:
    """Write auditable membership, exclusion, and exact-overlap tables."""

    membership_path = os.path.join(output_prefix, "reference_path_clusters.tsv")
    overlap_path = os.path.join(output_prefix, "reference_path_cluster_overlaps.tsv")
    cluster_by_path = {
        path: cluster.cluster_id
        for cluster in clusters
        for path in cluster.members
    }
    singleton_set = set(singletons)
    membership_fields = [
        "status", "cluster_id", "path_column", "path", "raw_depth_x",
        "depth_N", "nclose_count", "full_reference_chrom",
        "full_reference_ratio",
    ]
    rows = []
    for path, metadata in sorted(path_metadata.items(), key=lambda item: int(item[0])):
        if path in cluster_by_path:
            status = "clustered"
            cluster_id = cluster_by_path[path]
        elif path in singleton_set:
            status = "singleton"
            cluster_id = ""
        else:
            continue
        rows.append({
            "status": status,
            "cluster_id": cluster_id,
            "path_column": int(path),
            "path": metadata.get("location", ""),
            "raw_depth_x": f"{float(metadata.get('raw_depth', 0.0)):.10g}",
            "depth_N": f"{float(metadata.get('depth_n', 0.0)):.10g}",
            "nclose_count": int(metadata.get("nclose_count", 0)),
            "full_reference_chrom": "",
            "full_reference_ratio": "",
        })
    for excluded in excluded_full_reference:
        rows.append({
            "status": "excluded_full_reference",
            "cluster_id": "",
            "path_column": int(excluded["path_column"]),
            "path": excluded.get("location", ""),
            "raw_depth_x": f"{float(excluded.get('raw_depth', 0.0)):.10g}",
            "depth_N": f"{float(excluded.get('depth_n', 0.0)):.10g}",
            "nclose_count": int(excluded.get("nclose_count", 0)),
            "full_reference_chrom": excluded.get("full_reference_chrom", ""),
            "full_reference_ratio": f"{float(excluded.get('full_reference_ratio', 0.0)):.10g}",
        })
    rows.sort(key=lambda row: int(row["path_column"]))
    with open(membership_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=membership_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    overlap_fields = [
        "cluster_id", "path_a_column", "path_b_column", "chrom",
        "start", "end", "overlap_bp", "overlap_Mb",
    ]
    with open(overlap_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=overlap_fields, delimiter="\t")
        writer.writeheader()
        for cluster in clusters:
            for edge in cluster.edges:
                for overlap in edge.overlaps:
                    writer.writerow({
                        "cluster_id": cluster.cluster_id,
                        "path_a_column": int(edge.path_a),
                        "path_b_column": int(edge.path_b),
                        "chrom": overlap.chrom,
                        "start": overlap.start,
                        "end": overlap.end,
                        "overlap_bp": overlap.length,
                        "overlap_Mb": f"{overlap.length / M:.6f}",
                    })
    return {
        "membership_path": membership_path,
        "overlap_path": overlap_path,
    }


def write_reference_path_cluster_depth_table(
    *, output_prefix: str, depth_rows: Sequence[Mapping]
) -> str:
    """Write the numbers behind every depth figure so they stay auditable."""

    depth_path = os.path.join(output_prefix, "reference_path_cluster_depth.tsv")
    fields = [
        "cluster_id", "region_type", "region", "chrom", "start", "end",
        "length_bp", "observed_CN", "observed_CN_sd", "fitted_all_paths_CN",
        "cluster_CN", "other_paths_CN", "cluster_share_of_observed",
        "path_column", "path_CN",
    ]
    with open(depth_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(depth_rows)
    return depth_path


def _depth_table_rows(
    cluster: ReferencePathCluster,
    block_summaries: Sequence[RegionDepthSummary],
    path_summaries: Sequence[RegionDepthSummary],
) -> list[dict]:
    rows = []
    grouped = (
        ("shared_block", block_summaries),
        ("path_footprint", path_summaries),
    )
    for region_type, summaries in grouped:
        for summary in summaries:
            share = (
                summary.cluster_cn / summary.observed_cn
                if summary.observed_cn > 0 else 0.0
            )
            base = {
                "cluster_id": cluster.cluster_id,
                "region_type": region_type,
                "region": (
                    summary.label if region_type == "shared_block"
                    else _path_label(summary.carriers[0])
                ),
                "chrom": summary.chrom,
                "start": summary.start,
                "end": summary.end,
                "length_bp": summary.length,
                "observed_CN": f"{summary.observed_cn:.6g}",
                "observed_CN_sd": f"{summary.observed_sd:.6g}",
                "fitted_all_paths_CN": (
                    "" if summary.model_cn is None else f"{summary.model_cn:.6g}"
                ),
                "cluster_CN": f"{summary.cluster_cn:.6g}",
                "other_paths_CN": (
                    "" if summary.other_cn is None else f"{summary.other_cn:.6g}"
                ),
                "cluster_share_of_observed": f"{share:.6g}",
            }
            for path in cluster.members:
                rows.append({
                    **base,
                    "path_column": int(path),
                    "path_CN": f"{summary.member_cn.get(path, 0.0):.6g}",
                })
    return rows


def render_reference_path_clusters(
    *,
    output_prefix: str,
    cell_line: str,
    clusters: Sequence[ReferencePathCluster],
    singletons: Sequence,
    path_spans: Mapping,
    path_metadata: Mapping,
    excluded_full_reference: Sequence[Mapping],
    minimum_overlap: int,
    full_reference_ratio: float,
    modes: Sequence[str] = ("plain",),
    cluster_tracks: Mapping[str, ClusterDepthTracks] | None = None,
) -> dict:
    """Write tables and one figure per non-singleton reference component."""

    _remove_stale_cluster_figures(output_prefix)
    cluster_tracks = cluster_tracks or {}
    figures = []
    depth_rows: list[dict] = []
    for cluster in clusters:
        tracks = cluster_tracks.get(cluster.cluster_id)
        for mode in modes:
            if mode != "plain" and tracks is None:
                continue
            figures.append(render_reference_path_cluster(
                output_prefix=output_prefix,
                cell_line=cell_line,
                cluster=cluster,
                path_spans=path_spans,
                path_metadata=path_metadata,
                minimum_overlap=minimum_overlap,
                full_reference_ratio=full_reference_ratio,
                mode=mode,
                tracks=tracks,
            ))
        if tracks is not None:
            shared_blocks = build_shared_reference_blocks(cluster, path_spans)
            block_summaries, path_summaries = _cluster_region_summaries(
                cluster, path_spans, shared_blocks, tracks
            )
            depth_rows.extend(
                _depth_table_rows(cluster, block_summaries, path_summaries)
            )

    tables = write_reference_path_cluster_tables(
        output_prefix=output_prefix,
        clusters=clusters,
        singletons=singletons,
        path_metadata=path_metadata,
        excluded_full_reference=excluded_full_reference,
    )
    if depth_rows:
        tables["depth_path"] = write_reference_path_cluster_depth_table(
            output_prefix=output_prefix, depth_rows=depth_rows
        )
    return {
        "cluster_count": len(clusters),
        "singleton_count": len(singletons),
        "excluded_full_reference_count": len(excluded_full_reference),
        "figures": figures,
        **tables,
    }
