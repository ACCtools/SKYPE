"""Reusable Circos coverage renderer shared by SKYPE pipeline modes.

This module intentionally accepts already prepared, mode-specific data.  The
normal graph pipeline and the full-assembly pipeline therefore share the same
rendering code without coupling their event inference logic.
"""

from __future__ import annotations

import os
from collections import defaultdict
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from pycirclize import Circos
from matplotlib.lines import Line2D
from matplotlib.projections.polar import PolarAxes
from pycirclize.track import Track


M = 1_000_000
ABS_MAX_COVERAGE_RATIO = 3


def line_track_plot(
    line_track: Track,
    x: Sequence[float] | np.ndarray,
    y: Sequence[float] | np.ndarray,
    *,
    vmin: float = 0,
    vmax: float | None = None,
    arc: bool = True,
    **kwargs,
) -> None:
    """Plot a clipped line on a pycirclize track."""

    if len(x) != len(y):
        raise ValueError(f"List length is not match ({len(x)=}, {len(y)=})")

    rad = list(map(line_track.x_to_rad, x))
    vmax = max(y) if vmax is None else vmax
    r = [line_track._y_to_r(value, vmin, vmax) for value in y]

    if arc:
        plot_rad, plot_r = line_track._to_arc_radr(rad, r)
    else:
        plot_rad, plot_r = rad, r

    def plot_line(ax: PolarAxes) -> None:
        line, = ax.plot(plot_rad, plot_r, **kwargs)
        clip_patch = plt.Circle(
            (0, 0),
            max(line_track.r_plot_lim),
            transform=ax.transProjectionAffine + ax.transAxes,
        )
        clip_patch.set_visible(False)
        ax.add_patch(clip_patch)
        line.set_clip_path(clip_patch)

    line_track._plot_funcs.append(plot_line)


def line_track_circos_color(
    line_track: Track,
    x: Sequence[float] | np.ndarray,
    y1: Sequence[float] | np.ndarray,
    *,
    vmin: float = 0,
    vmax: float | None = None,
    target_value: Sequence[float],
    target_color: Sequence[str],
    **kwargs,
) -> None:
    """Fill a circular track with threshold-dependent colors."""

    y2 = 0
    if any(np.diff(target_value) < 0):
        raise ValueError("target_value must be sorted in ascending order")
    if len(target_color) != len(target_value) + 1:
        raise ValueError(
            "The number of target_color must be equal to len(target_value) + 1"
        )

    rad = list(map(line_track.x_to_rad, x))
    y_all = list(y1) + [y2]
    y2_values = [float(y2)] * len(x)
    vmin = min(y_all) if vmin is None else vmin
    vmax = max(y_all) if vmax is None else vmax
    plot_r2 = np.asarray([line_track._y_to_r(value, vmin, vmax) for value in y2_values])
    plot_r = np.asarray([line_track._y_to_r(value, vmin, vmax) for value in y1])
    plot_rad = rad

    def plot_fill_between(ax: PolarAxes) -> None:
        y1_arr = np.asarray(y1)
        segments = [(None, target_value[0])]
        segments.extend(
            (target_value[index - 1], target_value[index])
            for index in range(1, len(target_value))
        )
        segments.append((target_value[-1], None))

        for (lower, upper), color in zip(segments, target_color, strict=True):
            if lower is None:
                mask = y1_arr < upper
            elif upper is None:
                mask = y1_arr >= lower
            else:
                mask = (y1_arr >= lower) & (y1_arr < upper)
            mask[1:] = mask[1:] | mask[:-1]
            line = ax.fill_between(
                plot_rad,
                plot_r,
                plot_r2,
                where=mask,
                color=color,
                **kwargs,
            )
            clip_patch = plt.Circle(
                (0, 0),
                max(line_track.r_plot_lim),
                transform=ax.transProjectionAffine + ax.transAxes,
            )
            clip_patch.set_visible(False)
            ax.add_patch(clip_patch)
            line.set_clip_path(clip_patch)

    line_track._plot_funcs.append(plot_fill_between)


def rebin_depth_dataframe(frame: pd.DataFrame, n: int) -> pd.DataFrame:
    """Rebin depth rows exactly as the stage-31 coverage renderer expects."""

    rows = []
    for chrom, chrom_frame in frame.groupby("chr"):
        chrom_frame = chrom_frame.sort_values("st").reset_index(drop=True)
        for index in range(0, len(chrom_frame), n):
            chunk = chrom_frame.iloc[index:index + n]
            total_length = chunk["length"].sum()
            rows.append(
                {
                    "chr": chrom,
                    "st": chunk["st"].iloc[0],
                    "nd": chunk["nd"].iloc[-1],
                    "meandepth": np.sum(chunk["totaldepth"]) / total_length,
                }
            )
    return pd.DataFrame(rows)


def smooth_depth_for_coords(
    frame: pd.DataFrame,
    coords: Sequence[tuple[str, int]],
    n: int = 10,
) -> np.ndarray:
    """Return n-window-smoothed depth in the caller's coordinate order."""

    smoothed_by_coord: dict[tuple[str, int], float] = {}
    for chrom, chrom_frame in frame.groupby("chr"):
        chrom_frame = chrom_frame.sort_values("st").reset_index(drop=True)
        for index in range(0, len(chrom_frame), n):
            chunk = chrom_frame.iloc[index:index + n]
            total_length = chunk["length"].sum()
            value = float(np.sum(chunk["totaldepth"]) / total_length)
            for start in chunk["st"]:
                smoothed_by_coord[(str(chrom), int(start))] = value
    return np.asarray(
        [smoothed_by_coord[(str(chrom), int(start))] for chrom, start in coords],
        dtype=float,
    )


def render_total_coverage_circos(
    *,
    output_prefix: str,
    fig_prefix: str,
    chromosome_lengths: Mapping[str, int],
    cytobands_path: str,
    observed_depth: pd.DataFrame,
    prediction_coords: Sequence[tuple[str, int]],
    predicted_depth: np.ndarray,
    prediction_labels: np.ndarray,
    mean_depth: float,
    fragment_depth_per_chrom: Mapping[str, tuple[str, float]] | None = None,
    centromere_fragment_meta: Mapping[str, Mapping[str, int]] | None = None,
    breakend_links: Sequence[
        tuple[tuple[str, int], tuple[str, int], float, str]
    ] = (),
    telomere_arrows: Sequence[tuple[str, int, int, float, int]] = (),
) -> dict[str, str]:
    """Render the canonical stage-31 circular coverage figure.

    Event inference is deliberately outside this function.  A coverage-only
    mode passes empty ``breakend_links`` and ``telomere_arrows`` while keeping
    the exact same tracks, colors, legends, sizing, and output format.
    """

    fragment_depth_per_chrom = fragment_depth_per_chrom or {}
    centromere_fragment_meta = centromere_fragment_meta or {}
    prediction_coords = [(str(chrom), int(start)) for chrom, start in prediction_coords]
    predicted_depth = np.asarray(predicted_depth)
    prediction_labels = np.asarray(prediction_labels)
    if len(prediction_coords) != len(predicted_depth):
        raise ValueError("Prediction coordinate/depth length mismatch")
    if len(prediction_labels) != len(predicted_depth):
        raise ValueError("Prediction label/depth length mismatch")

    circos = Circos(dict(chromosome_lengths), space=3)
    circos.add_cytoband_tracks((95, 100), cytobands_path)
    target_color = ["blue", "red", "gray"]

    for sector in circos.sectors:
        sector.text(sector.name, r=107, size=8)
        sector.get_track("cytoband").xticks_by_interval(
            25 * M,
            label_size=6,
            tick_length=1,
            label_orientation="vertical",
            label_formatter=lambda value: f"{int(value / M)}",
        )
        sector.get_track("cytoband").xticks_by_interval(
            5 * M,
            tick_length=0.5,
            label_formatter=lambda value: "",
        )

        sector_depth = observed_depth.query(f'chr == "{sector.name}"')
        sector_starts = []
        sector_indices = []
        for index, (chrom, start) in enumerate(prediction_coords):
            if chrom == sector.name:
                sector_starts.append(start)
                sector_indices.append(index)
        sector_starts = np.asarray(sector_starts)
        sector_indices = np.asarray(sector_indices, dtype=int)
        order = np.argsort(sector_starts)
        sector_starts = sector_starts[order]
        sector_indices = sector_indices[order]

        cn_track = sector.add_track((60, 93))
        cn_track.axis()
        cn_track.grid(y_grid_num=7, zorder=-2)
        max_start = max(sector_depth["st"])
        cn_track.line(
            [0, max_start],
            [0.5 * mean_depth] * 2,
            color="blue",
            vmax=ABS_MAX_COVERAGE_RATIO * mean_depth,
            zorder=-1,
            alpha=0.5,
        )
        cn_track.line(
            [0, max_start],
            [2 * mean_depth] * 2,
            color="red",
            vmax=ABS_MAX_COVERAGE_RATIO * mean_depth,
            zorder=-1,
            alpha=0.5,
        )
        line_track_plot(
            cn_track,
            x=sector_starts,
            y=predicted_depth[sector_indices],
            vmax=ABS_MAX_COVERAGE_RATIO * mean_depth,
            color="green",
            linewidth=0.5,
            zorder=10,
        )
        line_track_plot(
            cn_track,
            x=sector_depth["st"].to_list(),
            y=sector_depth["meandepth"].to_list(),
            vmax=ABS_MAX_COVERAGE_RATIO * mean_depth,
            color="black",
            linewidth=0.2,
        )
        line_track_circos_color(
            cn_track,
            x=sector_depth["st"].to_list(),
            y1=sector_depth["meandepth"].to_list(),
            vmax=ABS_MAX_COVERAGE_RATIO * mean_depth,
            target_color=["#6baed6", "#969696", "#fb6a4a"],
            target_value=[0.5 * mean_depth, 2 * mean_depth],
        )

        if sector.name in fragment_depth_per_chrom:
            side, fragment_weight = fragment_depth_per_chrom[sector.name]
            info = centromere_fragment_meta[sector.name]
            if side == "right":
                start_low, start_high = info["mid"], info["chr_len"]
            else:
                start_low, start_high = 0, info["mid"]
            cn_track.line(
                [start_low, start_low, start_high, start_high],
                [0, fragment_weight, fragment_weight, 0],
                color="purple",
                linestyle="--",
                linewidth=0.7,
                vmax=ABS_MAX_COVERAGE_RATIO * mean_depth,
                zorder=11,
            )

        color_track = sector.add_track((53, 58))
        color_track.axis()
        line_track_circos_color(
            color_track,
            x=sector_starts,
            y1=prediction_labels[sector_indices],
            vmax=0.5,
            target_color=target_color,
            target_value=[2, 4],
        )

    offset = 0.03
    lower = 0
    upper = mean_depth * 2
    cmap = sns.color_palette("rocket_r", as_cmap=True)
    line_style_by_event = {
        "breakend": "-",
        "inversion": "-.",
        "indel": "--",
        "amplicon": (0, (5, 2)),
        "virtual_inv": (0, (2, 2)),
    }
    for index, (location1, location2, copy_depth, event_type) in enumerate(
        sorted(breakend_links, key=lambda value: value[2])
    ):
        normalized_copy_depth = np.clip(
            (copy_depth - lower) / (upper - lower), 0, 1
        )
        if copy_depth >= 0:
            circos.link_line(
                location1,
                location2,
                color=cmap(normalized_copy_depth),
                linestyle=line_style_by_event[event_type],
                lw=1,
                zorder=index,
            )

    circos.colorbar(
        vmin=0,
        vmax=2,
        bounds=(1.01 + offset, 0.825, 0.02, 0.1),
        cmap=cmap,
        label="Breakend CN",
        orientation="vertical",
        colorbar_kws=dict(
            format=mpl.ticker.FixedFormatter(["0", "2N", "4N"]),
            ticks=mpl.ticker.FixedLocator([0, 1, 2]),
        ),
        label_kws=dict(fontsize=9),
        tick_kws=dict(labelsize=8),
    )

    arrows_by_chrom = defaultdict(list)
    for arrow in telomere_arrows:
        arrows_by_chrom[arrow[0]].append(arrow)
    for sector in circos.sectors:
        track = sector.add_track((58.2, 59.8))
        for _chrom, start, end, copy_depth, zorder in arrows_by_chrom[sector.name]:
            normalized_copy_depth = np.clip(
                (copy_depth - lower) / (upper - lower), 0, 1
            )
            try:
                track.arrow(
                    start,
                    end,
                    fc=cmap(normalized_copy_depth),
                    zorder=zorder,
                )
            except Exception:
                # Preserve the original renderer's out-of-range tolerance.
                pass

    fig = circos.plotfig(figsize=(10, 10), dpi=500)
    legend_location = "upper left"
    cn_legend = circos.ax.legend(
        handles=[
            Line2D([], [], color="black", label="Real CN", linestyle="-"),
            Line2D([], [], color="green", label="SKYPE CN", linestyle="-"),
        ],
        loc=legend_location,
        bbox_to_anchor=(0.85 + offset, 1.03),
        fontsize=8,
        title="CN lines",
        handlelength=2,
    )
    circos.ax.add_artist(cn_legend)

    prediction_legend = circos.ax.legend(
        loc=legend_location,
        handles=[
            mpl.patches.Patch(color=color, label=label)
            for color, label in zip(
                target_color, ["Sucess", "Failed", "Ignored"], strict=True
            )
        ],
        bbox_to_anchor=(0.97 + offset, 1.03),
        fontsize=8,
        title="Predict label",
        handlelength=2,
    )
    circos.ax.add_artist(prediction_legend)

    breakend_legend = circos.ax.legend(
        loc=legend_location,
        handles=[
            Line2D([], [], color="black", label="Breakend", linestyle="-"),
            Line2D([], [], color="black", label="Inversion", linestyle="-."),
            Line2D([], [], color="black", label="Indel", linestyle="--"),
            Line2D(
                [], [], color="black", label="Amplicon", linestyle=(0, (5, 2))
            ),
            Line2D(
                [],
                [],
                color="black",
                label="Virtual inversion",
                linestyle=(0, (2, 2)),
            ),
        ],
        bbox_to_anchor=(0.85 + offset, 0.94),
        fontsize=8,
        title="Breakend lines",
        handlelength=2,
    )
    circos.ax.add_artist(breakend_legend)

    pdf_path = os.path.join(output_prefix, f"total_cov{fig_prefix}.pdf")
    png_path = os.path.join(output_prefix, f"total_cov{fig_prefix}.png")
    fig.savefig(pdf_path)
    fig.savefig(png_path)
    plt.close(fig)
    return {"pdf_path": pdf_path, "png_path": png_path}
