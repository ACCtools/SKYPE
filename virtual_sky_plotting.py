"""Reusable canonical Virtual SKY renderer for normal and assembly pipelines."""

from __future__ import annotations

import math
import os
from typing import Mapping

import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator


M = 1_000_000
KARYOTYPE_MIN_SEGMENT_LENGTH = 1_000_000

CHR_COLORS = {
    "chr1": "#FFC000", "chr2": "#FF0000", "chr3": "#00B050",
    "chr4": "#0070C0", "chr5": "#7030A0", "chr6": "#92D050",
    "chr7": "#C00000", "chr8": "#00B0F0", "chr9": "#FFFF00",
    "chr10": "#F08080", "chr11": "#48D1CC", "chr12": "#D2B48C",
    "chr13": "#B0C4DE", "chr14": "#FFA07A", "chr15": "#CD5C5C",
    "chr16": "#E9967A", "chr17": "#BDB76B", "chr18": "#D8BFD8",
    "chr19": "#BC8F8F", "chr20": "#FFD700", "chr21": "#7FFFD4",
    "chr22": "#ADFF2F", "chrX": "#6495ED", "chrY": "#EE82EE",
}

_NO_JUNCTION_READ_COUNT = object()
_NO_JUNCTION_READ_DEPTH = object()


def weighted_vaf_read_depth(record) -> float | None:
    """Return local-coverage/VAF-derived junction depth in raw coverage units."""

    counts = record.get("read_counts", {})
    try:
        junction_reads = int(counts.get("d2", 0))
        normal_reads_a = int(counts.get("d1", 0))
        normal_reads_b = int(counts.get("d3", 0))
    except (TypeError, ValueError):
        return None

    estimate = record.get("read_count_depth_estimate", {})
    vaf = record.get("vaf", {})
    side_data = []
    for side, normal_reads in (("a", normal_reads_a), ("b", normal_reads_b)):
        try:
            local_depth = float(estimate.get(f"chr_{side}_local_depth"))
        except (TypeError, ValueError):
            continue
        if not math.isfinite(local_depth):
            continue

        try:
            side_vaf = float(vaf.get(f"chr_{side}"))
        except (TypeError, ValueError):
            total_reads = junction_reads + normal_reads
            if total_reads <= 0:
                continue
            side_vaf = junction_reads / total_reads
        if not math.isfinite(side_vaf):
            continue

        support = junction_reads + normal_reads
        if support > 0:
            side_data.append((local_depth * side_vaf, support))

    total_support = sum(support for _depth, support in side_data)
    if total_support <= 0:
        return None
    return sum(depth * support for depth, support in side_data) / total_support


def chromosome_sort_key(chrom: str) -> tuple[int, str]:
    if chrom.startswith("chr"):
        suffix = chrom[3:]
        if suffix.isdigit():
            return int(suffix), chrom
        if suffix == "X":
            return 24, chrom
        if suffix == "Y":
            return 25, chrom
    return 1_000_000_000, chrom


def chrom_to_iscn(chrom: str) -> str:
    return chrom[3:] if chrom.startswith("chr") else chrom


def virtual_event_label(event_type: str) -> str:
    return {"d": "del", "i": "ins", "v": "inv"}[event_type]


def mark_overlapping_texts_with_arrows(ax, texts, min_gap=5) -> None:
    texts_sorted = sorted(texts, key=lambda text: text.get_position()[1])
    for index in range(1, len(texts_sorted)):
        previous = texts_sorted[index - 1]
        current = texts_sorted[index]
        _x, previous_y = previous.get_position()
        x, current_y = current.get_position()
        if current_y - previous_y < min_gap:
            new_y = previous_y + min_gap
            current.set_position((x, new_y))
            ax.annotate(
                "",
                xy=(x, new_y),
                xycoords="data",
                xytext=(5, current_y),
                textcoords="data",
                arrowprops=dict(arrowstyle="-|>", color="black", lw=0.5),
            )


def plot_virtual_chromosome(
    ax,
    segments_data,
    maxh,
    cell_col,
    default_cell_col,
    label=None,
    karyotype_str=None,
    event_labels=None,
    junction_read_counts=None,
    junction_read_depths=None,
) -> None:
    _path, segments = segments_data
    x_center = 0
    width = 10
    total_length = sum(segment_length for _chrom, segment_length in segments)
    radius = width / 2.0
    round_radius = min(radius, total_length / 4)
    clip_patch = patches.FancyBboxPatch(
        (x_center - radius, 0),
        width,
        total_length,
        boxstyle=f"round,pad=0,rounding_size={round_radius}",
        facecolor="none",
        edgecolor="none",
    )
    ax.add_patch(clip_patch)

    current_y = 0
    text_objects = []
    text_x = x_center + radius * 2
    text_args = dict(ha="left", va="center", fontsize=10, color="black")
    event_labels = event_labels or []
    junction_read_counts = junction_read_counts or []
    junction_read_depths = junction_read_depths or []
    event_positions = [event_y for event_y, _label in event_labels]
    for index, (real_chrom, segment_length) in enumerate(segments):
        rect = patches.Rectangle(
            (x_center - radius, current_y),
            width,
            segment_length,
            facecolor=CHR_COLORS.get(real_chrom[0], "gray"),
            edgecolor="none",
        )
        rect.set_clip_path(clip_patch)
        ax.add_patch(rect)
        annotated = any(abs(current_y - event_y) < 1e-6 for event_y in event_positions)
        if 0 < index < len(segments) and not annotated:
            previous_chrom = segments[index - 1][0]
            if real_chrom[0] != previous_chrom[0]:
                text_label = (
                    f"t({chrom_to_iscn(previous_chrom[0])};"
                    f"{chrom_to_iscn(real_chrom[0])})"
                )
            else:
                text_label = (
                    f"t({chrom_to_iscn(previous_chrom[0])}{previous_chrom[1]};"
                    f"{chrom_to_iscn(real_chrom[0])}{real_chrom[1]})"
                )
            matching_read_count = next(
                (
                    read_count
                    for boundary_index, read_count in junction_read_counts
                    if index == boundary_index
                ),
                _NO_JUNCTION_READ_COUNT,
            )
            matching_read_depth = next(
                (
                    read_depth
                    for boundary_index, read_depth in junction_read_depths
                    if index == boundary_index
                ),
                _NO_JUNCTION_READ_DEPTH,
            )
            if matching_read_depth is not _NO_JUNCTION_READ_DEPTH:
                read_depth_text = (
                    "NA" if matching_read_depth is None
                    else f"{float(matching_read_depth):.1f}x"
                )
                text_label = f"{text_label} : {read_depth_text}"
            elif matching_read_count is not _NO_JUNCTION_READ_COUNT:
                read_count_text = (
                    "NA" if matching_read_count is None
                    else str(int(matching_read_count))
                )
                text_label = f"{text_label} : {read_count_text}"
            text_objects.append(ax.text(text_x, current_y, text_label, **text_args))
        current_y += segment_length

    for event_y, event_label in event_labels:
        text_objects.append(
            ax.text(text_x, max(0, min(100, event_y)), event_label, **text_args)
        )
    mark_overlapping_texts_with_arrows(ax, text_objects, min_gap=5)
    ax.add_patch(
        patches.FancyBboxPatch(
            (x_center - radius, 0),
            width,
            total_length,
            boxstyle=f"round,pad=0,rounding_size={round_radius}",
            facecolor="none",
            edgecolor="black",
            linewidth=1.2,
            clip_on=False,
        )
    )
    if label:
        ax.text(x_center, -5, label, ha="center", va="top", fontsize=10)
    if karyotype_str:
        ax.text(x_center, -10, karyotype_str, ha="center", va="top", fontsize=10)
    ax.set_xlim(0, 60 * cell_col / default_cell_col)
    ax.set_ylim(0, 100)
    ax.axis("off")


def plot_indel(
    ax,
    indel,
    maxh,
    cell_col,
    default_cell_col,
    chromosome_length,
    label=None,
) -> None:
    x_center = 0
    width = 10
    real_chrom = indel[4]
    radius = width / 2.0
    round_radius = min(radius, chromosome_length / 4)
    chrom_color = CHR_COLORS.get(real_chrom, "gray")
    background = "white" if indel[0] in {"i", "v"} else chrom_color
    foreground = chrom_color if indel[0] in {"i", "v"} else "white"
    clip_patch = patches.FancyBboxPatch(
        (x_center - radius, 0), width, chromosome_length,
        boxstyle=f"round,pad=0,rounding_size={round_radius}",
        facecolor="none", edgecolor="none",
    )
    ax.add_patch(clip_patch)
    ax.add_patch(
        patches.FancyBboxPatch(
            (x_center - radius, 0), width, chromosome_length,
            boxstyle=f"round,pad=0,rounding_size={round_radius}",
            facecolor=background, edgecolor="black", linestyle="dashed",
            linewidth=1.2, clip_on=False,
        )
    )
    y1 = indel[1] / maxh * 100
    y2 = indel[2] / maxh * 100
    middle_y = (y1 + y2) / 2
    rect = patches.Rectangle(
        (x_center - radius, y1), width, y2 - y1,
        facecolor=foreground, edgecolor="none",
    )
    rect.set_clip_path(clip_patch)
    ax.add_patch(rect)
    x_start = x_center + radius
    x_end = x_center + 2 * radius
    ax.plot([x_start, x_end], [y1, middle_y], linewidth=1, color="black")
    ax.plot([x_start, x_end], [y2, middle_y], linewidth=1, color="black")
    event_label = virtual_event_label(indel[0])
    text = ax.text(
        x_end + radius / 5,
        middle_y,
        f"{event_label}({chrom_to_iscn(indel[-2])})",
        ha="left", va="center", fontsize=10, color="black",
    )
    mark_overlapping_texts_with_arrows(ax, [text], min_gap=5)
    if label:
        ax.text(x_center, -5, label, ha="center", va="top", fontsize=10)
    ax.text(
        x_center, -10, f"{event_label}({chrom_to_iscn(real_chrom)})",
        ha="center", va="top", fontsize=10,
    )
    ax.set_xlim(0, 60 * cell_col / default_cell_col)
    ax.set_ylim(0, 100)
    ax.axis("off")


def plot_centromere_fragment(
    ax, frag, maxh, cell_col, default_cell_col, chromosome_length_norm, label=None
) -> None:
    chrom, side, midpoint, chromosome_length = frag
    x_center = 0
    width = 10
    radius = width / 2.0
    round_radius = (
        min(radius, chromosome_length_norm / 4)
        if chromosome_length_norm > 0 else radius
    )
    clip_patch = patches.FancyBboxPatch(
        (x_center - radius, 0), width, chromosome_length_norm,
        boxstyle=f"round,pad=0,rounding_size={round_radius}",
        facecolor="none", edgecolor="none",
    )
    ax.add_patch(clip_patch)
    ax.add_patch(
        patches.FancyBboxPatch(
            (x_center - radius, 0), width, chromosome_length_norm,
            boxstyle=f"round,pad=0,rounding_size={round_radius}",
            facecolor="white", edgecolor="black", linestyle="dashed",
            linewidth=1.2, clip_on=False,
        )
    )
    midpoint_norm = midpoint / chromosome_length * chromosome_length_norm
    y1, y2 = (
        (midpoint_norm, chromosome_length_norm)
        if side == "right" else (0, midpoint_norm)
    )
    rect = patches.Rectangle(
        (x_center - radius, y1), width, y2 - y1,
        facecolor=CHR_COLORS.get(chrom, "gray"), edgecolor="none",
    )
    rect.set_clip_path(clip_patch)
    ax.add_patch(rect)
    middle_y = (y1 + y2) / 2
    x_start = x_center + radius
    x_end = x_center + 2 * radius
    ax.plot([x_start, x_end], [y1, middle_y], linewidth=1, color="black")
    ax.plot([x_start, x_end], [y2, middle_y], linewidth=1, color="black")
    arm = "q" if side == "right" else "p"
    text = ax.text(
        x_end + radius / 5, middle_y, f"{chrom}{arm}",
        ha="left", va="center", fontsize=10, color="black",
    )
    mark_overlapping_texts_with_arrows(ax, [text], min_gap=5)
    if label:
        ax.text(x_center, -5, label, ha="center", va="top", fontsize=10)
    ax.text(
        x_center, -10, f"i({chrom_to_iscn(chrom)})({arm}10)",
        ha="center", va="top", fontsize=10,
    )
    ax.set_xlim(0, 60 * cell_col / default_cell_col)
    ax.set_ylim(0, 100)
    ax.axis("off")


def plot_scale_bar(ax, chrom, maxh, chromosome_lengths: Mapping[str, int]) -> None:
    ax.set_xlim(-5, 5)
    ax.axis(True)
    ax.set_ylim(0, maxh / M)
    ax.set_ylabel("Scale bar (Mb)", rotation=90)
    ax.yaxis.set_major_locator(MultipleLocator(25))
    ax.yaxis.set_minor_locator(MultipleLocator(5))
    ax.tick_params(axis="y", which="major", length=5)
    ax.tick_params(axis="y", which="minor", length=3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines["left"].set_position(("data", 0))
    ax.xaxis.set_visible(False)
    if chrom in chromosome_lengths:
        ax.vlines(
            x=0.5, ymin=0, ymax=chromosome_lengths[chrom] / M,
            colors=CHR_COLORS[chrom], linewidth=4, alpha=0.6,
        )


def plot_chromosome_name(ax, chrom: str) -> None:
    ax.text(
        0, 0.5, chrom, transform=ax.transAxes,
        ha="center", va="center", fontsize=15, rotation=90,
    )


def render_karyotype_diagram(
    *,
    output_prefix: str,
    cell_line: str,
    chromosome_lengths: Mapping[str, int],
    grouped_norm_data,
    display_indel,
    virtual_inv_display,
    fragment_display,
    maxh,
    path_depth_n,
    path_karyotype,
    path_event_labels=None,
    path_junction_read_counts=None,
    path_junction_read_depths=None,
    path_depth_raw=None,
    fig_prefix="",
) -> dict:
    """Render the exact stage-30 Virtual SKY layout from prepared path data."""

    path_event_labels = path_event_labels or {}
    path_junction_read_counts = path_junction_read_counts or {}
    path_junction_read_depths = path_junction_read_depths or {}
    show_raw_path_depth = path_depth_raw is not None
    path_depth_raw = path_depth_raw or {}
    columns = 10
    display_chroms = sorted(
        set(grouped_norm_data) | set(display_indel), key=chromosome_sort_key
    )
    rows = 0
    for chrom in display_chroms:
        count = len(grouped_norm_data.get(chrom, [])) + len(display_indel.get(chrom, []))
        rows += (count - 1) // columns + 1 if count else 0
    if fragment_display:
        rows += (len(fragment_display) - 1) // columns + 1
    if virtual_inv_display:
        rows += (len(virtual_inv_display) - 1) // columns + 1
    rows = max(rows, 1)

    cell_col = 3
    cell_row = 4.5
    default_cell_col = 1.8
    prefix_ratios = [0.5, 1]
    width_ratios = prefix_ratios + [cell_col] * columns
    height_ratios = [1.5] + [cell_row] * rows
    fig, axes = plt.subplots(
        len(height_ratios), len(width_ratios),
        figsize=(sum(width_ratios), sum(height_ratios)),
        gridspec_kw=dict(
            width_ratios=width_ratios, height_ratios=height_ratios,
            wspace=0, hspace=0.2,
        ),
    )
    for axis_row in axes:
        for axis in axis_row:
            axis.axis("off")

    grouped_items = [
        (chrom, grouped_norm_data.get(chrom, [])) for chrom in display_chroms
    ]
    current_row = 1
    for chrom, data_list in grouped_items:
        indels = display_indel.get(chrom, [])
        first_row = current_row
        plot_chromosome_name(axes[current_row][0], chrom)
        for index, data in enumerate(data_list):
            path = data[0]
            depth_n = float(path_depth_n[path])
            depth_label = f"{round(depth_n, 2)}N"
            if show_raw_path_depth and path in path_depth_raw:
                depth_label = f"{float(path_depth_raw[path]):.2f}x"
            plot_virtual_chromosome(
                axes[index // columns + current_row][index % columns + len(prefix_ratios)],
                data, maxh, cell_col, default_cell_col,
                label=depth_label,
                karyotype_str=path_karyotype.get(path) or "",
                event_labels=path_event_labels.get(path, []),
                junction_read_counts=path_junction_read_counts.get(path, []),
                junction_read_depths=path_junction_read_depths.get(path, []),
            )
        for offset, indel in enumerate(indels):
            index = offset + len(data_list)
            plot_indel(
                axes[index // columns + current_row][index % columns + len(prefix_ratios)],
                indel, maxh, cell_col, default_cell_col,
                chromosome_lengths[chrom] / maxh * 100,
                label=f"{round(indel[3], 2)}N",
            )
        count = len(data_list) + len(indels)
        current_row += (count - 1) // columns + 1 if count else 0
        for row_index in range(first_row, current_row):
            plot_scale_bar(
                axes[row_index][len(prefix_ratios) - 1],
                chrom, maxh, chromosome_lengths,
            )

    if virtual_inv_display:
        first_row = current_row
        plot_chromosome_name(axes[current_row][0], "virtual inv")
        for index, indel in enumerate(virtual_inv_display):
            chrom = indel[4]
            plot_indel(
                axes[index // columns + current_row][index % columns + len(prefix_ratios)],
                indel, maxh, cell_col, default_cell_col,
                chromosome_lengths[chrom] / maxh * 100,
                label=f"{round(indel[3], 2)}N",
            )
        current_row += (len(virtual_inv_display) - 1) // columns + 1
        for row_index in range(first_row, current_row):
            plot_scale_bar(
                axes[row_index][len(prefix_ratios) - 1],
                "virtual inv", maxh, chromosome_lengths,
            )

    if fragment_display:
        first_row = current_row
        plot_chromosome_name(axes[current_row][0], "cen frag")
        for index, fragment in enumerate(fragment_display):
            chrom, side, midpoint, chromosome_length, depth_n = fragment
            plot_centromere_fragment(
                axes[index // columns + current_row][index % columns + len(prefix_ratios)],
                (chrom, side, midpoint, chromosome_length),
                maxh, cell_col, default_cell_col,
                chromosome_length / maxh * 100,
                label=f"{round(depth_n, 2)}N",
            )
        current_row += (len(fragment_display) - 1) // columns + 1
        for row_index in range(first_row, current_row):
            plot_scale_bar(
                axes[row_index][len(prefix_ratios) - 1],
                "cen", maxh, chromosome_lengths,
            )

    legend_handles = [
        patches.Patch(
            facecolor=color, edgecolor="black", linewidth=0.5, label=chrom
        )
        for chrom, color in sorted(CHR_COLORS.items(), key=lambda item: chromosome_sort_key(item[0]))
    ]
    grid = axes[0, 0].get_gridspec()
    for axis in axes[0]:
        axis.remove()
    legend_columns = 8
    reorder = lambda values: sum(
        (values[index::legend_columns] for index in range(legend_columns)), []
    )
    legend_axis = fig.add_subplot(grid[0, :])
    legend_axis.legend(
        handles=reorder(legend_handles),
        labels=reorder([handle.get_label() for handle in legend_handles]),
        ncol=legend_columns, loc="center", fontsize=9, frameon=True,
        handlelength=1.0, handleheight=1.0, labelspacing=0.8,
    )
    title = f"{cell_line} Virtual SKY result"
    if show_raw_path_depth:
        title += "\npath = fitted depth (x); junction = VAF-derived depth (x)"
    legend_axis.set_title(title, fontsize=15)
    legend_axis.axis("off")

    pdf_path = os.path.join(output_prefix, f"virtual_sky{fig_prefix}.pdf")
    png_path = os.path.join(output_prefix, f"virtual_sky{fig_prefix}.png")
    fig.savefig(pdf_path)
    fig.savefig(png_path)
    plt.close(fig)

    karyotype_rows = []
    for chrom, data_list in grouped_items:
        for path, _normalized in data_list:
            iscn = path_karyotype.get(path)
            if iscn is not None:
                karyotype_rows.append((iscn, path_depth_n[path]))
        for indel in display_indel.get(chrom, []):
            event_type, start, end, depth_n, event_chrom, _path = indel
            if abs(end - start) >= KARYOTYPE_MIN_SEGMENT_LENGTH:
                karyotype_rows.append(
                    (f"{virtual_event_label(event_type)}({chrom_to_iscn(event_chrom)})", depth_n)
                )
    for event_type, start, end, depth_n, chrom, _name in virtual_inv_display:
        if abs(end - start) >= KARYOTYPE_MIN_SEGMENT_LENGTH:
            karyotype_rows.append(
                (f"{virtual_event_label(event_type)}({chrom_to_iscn(chrom)})", depth_n)
            )
    for chrom, side, _midpoint, _length, depth_n in fragment_display:
        arm = "q" if side == "right" else "p"
        karyotype_rows.append((f"i({chrom_to_iscn(chrom)})({arm}10)", depth_n))

    karyotype_path = os.path.join(output_prefix, f"karyotype{fig_prefix}.txt")
    with open(karyotype_path, "wt", encoding="utf-8") as handle:
        handle.write("karyotype\tdepth\n")
        for iscn, depth_n in karyotype_rows:
            handle.write(f"{iscn}\t{round(float(depth_n), 2)}N\n")
    return {
        "displayed_paths": sum(len(data) for data in grouped_norm_data.values()),
        "karyotype_rows": karyotype_rows,
        "pdf_path": pdf_path,
        "png_path": png_path,
        "karyotype_path": karyotype_path,
    }
