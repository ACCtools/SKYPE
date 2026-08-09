"""Compare fitted NClose depth with depth inferred from raw-read counts.

NCloses are grouped by co-occurrence in a displayed cluster path.  Each
connected component receives its own auxiliary figure, plus a multi-page PDF
and machine-readable TSV/JSON summaries.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import pickle
from collections import defaultdict

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.stats import spearmanr

from nclose_tracking import calculate_event_weights, load_event_catalog, load_path_usage
from skype_utils import load_pipeline_mode, pipeline_mode_is_karyotype


TARGET_DEPTH_N = 0.2
ROUGH_RATIO_LOW = 0.5
ROUGH_RATIO_HIGH = 2.0
COUNT_RESULT_PKL = "cluster_nclose_count_result.pkl"
SUMMARY_TSV = "nclose_read_depth_comparison.tsv"
METRICS_JSON = "nclose_read_depth_metrics.json"
GROUP_DIR = "nclose_read_depth_groups"
GROUP_PDF = "nclose_read_depth_groups.pdf"
SCATTER_PNG = "nclose_read_depth_scatter.png"
SCATTER_PDF = "nclose_read_depth_scatter.pdf"

SHORT_STATUS = {
    "AGREE_WITHIN_2X": "AGREE",
    "DISAGREE_OVER_2X": "DISAGREE",
    "UNASSESSABLE_REPEAT_OR_LOWMAPQ": "UNASSESSABLE",
    "UNCOUNTABLE_QUERY_GAP": "UNCOUNTABLE",
    "HIGH_QUALITY_ZERO": "ZERO",
    "FILTERED_CURRENT_DEPTH": "FILTERED",
    "NO_READ_DEPTH_ESTIMATE": "NO_ESTIMATE",
}

SCATTER_STYLES = {
    "ASSESSABLE": {
        "marker": "o", "color": "#2166ac", "label": "Assessable",
    },
    "ASSESSABLE_DISCORDANT": {
        "marker": "D", "color": "#b2182b", "label": "Assessable, >25% from y=x",
    },
    "REPEAT_OR_LOWMAPQ": {
        "marker": "^", "color": "#f4a261", "label": "Low MAPQ / unassessable",
    },
    "ZERO_SUPPORT": {
        "marker": "s", "color": "#d73027", "label": "High-quality zero support",
    },
    "OTHER_PROBLEM": {
        "marker": "P", "color": "#7b3294", "label": "Other problem",
    },
}


def median_depth(path: str) -> float:
    opener = gzip.open if str(path).endswith(".gz") else open
    values = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 8 and fields[0] != "chrM":
                values.append(float(fields[7]))
    if not values:
        raise ValueError(f"No usable depth rows in {path}")
    return float(np.median(values))


def chrom_label(chrom: str) -> str:
    return chrom[3:] if str(chrom).startswith("chr") else str(chrom)


def event_label(event: dict) -> str:
    start_chrom = chrom_label(event.get("start_chr", "?"))
    end_chrom = chrom_label(event.get("end_chr", "?"))
    if start_chrom != end_chrom:
        return f"t({start_chrom};{end_chrom})"
    return (
        f"t({start_chrom}{event.get('start_dir', '?')};"
        f"{end_chrom}{event.get('end_dir', '?')})"
    )


def optional_float(value):
    if value is None:
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if np.isfinite(value) else None


def normalized_read_depth(depth_estimate: dict, haploid_depth: float) -> float | None:
    """Normalize the raw-coverage estimate, with old-pickle compatibility."""

    raw_depth = optional_float(depth_estimate.get("weighted_read_depth"))
    if raw_depth is not None:
        if not np.isfinite(haploid_depth) or float(haploid_depth) <= 0:
            return None
        return raw_depth / float(haploid_depth)
    return optional_float(depth_estimate.get("weighted_read_depth_N"))


def event_key_from_record(record: dict):
    event_key = record.get("event_key")
    if event_key is not None:
        return tuple(event_key)
    return tuple(sorted(int(index) for index in record["nclose_key"]))


def connected_components(event_keys, displayed_columns, path_usage):
    event_keys = set(event_keys)
    adjacency = {event_key: set() for event_key in event_keys}
    carriers = defaultdict(set)
    for column in displayed_columns:
        used = set(path_usage[column]) & event_keys
        for event_key in used:
            carriers[event_key].add(int(column))
            adjacency[event_key].update(used - {event_key})

    components = []
    unseen = set(event_keys)
    while unseen:
        seed = min(unseen, key=repr)
        stack = [seed]
        component = set()
        while stack:
            current = stack.pop()
            if current in component:
                continue
            component.add(current)
            stack.extend(adjacency[current] - component)
        unseen -= component
        components.append(component)
    components.sort(key=lambda component: min(map(repr, component)))
    return components, carriers


def classify_row(row: dict) -> str:
    if not row["countable"]:
        return "UNCOUNTABLE_QUERY_GAP"
    if row["overlaps_censat"] or min(row["mapq_a"], row["mapq_b"]) < 60:
        return "UNASSESSABLE_REPEAT_OR_LOWMAPQ"
    if row["junction_reads"] == 0:
        return "HIGH_QUALITY_ZERO"
    if row["current_depth_N"] <= 0:
        return "FILTERED_CURRENT_DEPTH"
    if row["read_depth_N"] is None:
        return "NO_READ_DEPTH_ESTIMATE"
    ratio = row["read_to_current_ratio"]
    if ROUGH_RATIO_LOW <= ratio <= ROUGH_RATIO_HIGH:
        return "AGREE_WITHIN_2X"
    return "DISAGREE_OVER_2X"


def safe_metric(values, default=None):
    return default if values is None or not np.isfinite(values) else float(values)


def metric_text(value) -> str:
    return "NA" if value is None else f"{float(value):.3f}"


def scatter_category(row: dict) -> str:
    status = row["status"]
    if status in {"AGREE_WITHIN_2X", "DISAGREE_OVER_2X"}:
        current_depth = row.get("current_depth_N")
        read_depth = row.get("read_depth_N")
        if (
            current_depth is not None
            and read_depth is not None
            and float(current_depth) > 0
            and abs(float(read_depth) / float(current_depth) - 1) > 0.25
        ):
            return "ASSESSABLE_DISCORDANT"
        return "ASSESSABLE"
    if status == "UNASSESSABLE_REPEAT_OR_LOWMAPQ":
        return "REPEAT_OR_LOWMAPQ"
    if status == "HIGH_QUALITY_ZERO":
        return "ZERO_SUPPORT"
    return "OTHER_PROBLEM"


def scatter_is_included(row: dict) -> bool:
    """Exclude an event when either endpoint overlaps CENSAT."""

    return not bool(row.get("overlaps_censat", False))


def proportionality_metrics(rows: list[dict]) -> dict:
    eligible = [
        row for row in rows
        if row["status"] in {"AGREE_WITHIN_2X", "DISAGREE_OVER_2X"}
        and row["current_depth_N"] > 0
        and row["read_depth_N"] is not None
    ]
    metrics = {
        "rough_ratio_lower": ROUGH_RATIO_LOW,
        "rough_ratio_upper": ROUGH_RATIO_HIGH,
        "all_event_count": len(rows),
        "high_quality_positive_event_count": len(eligible),
    }
    if not eligible:
        metrics.update({
            "slope_through_origin": None,
            "junction_reads_per_fitted_N_slope": None,
            "pearson_r": None,
            "spearman_rho": None,
            "median_read_to_current_ratio": None,
            "rough_agreement_fraction": None,
            "within_25_percent_fraction": None,
            "within_50_percent_fraction": None,
            "mean_absolute_error_N": None,
            "root_mean_square_error_N": None,
            "median_absolute_relative_error": None,
        })
        return metrics

    current = np.asarray([row["current_depth_N"] for row in eligible], dtype=float)
    read = np.asarray([row["read_depth_N"] for row in eligible], dtype=float)
    junction_reads = np.asarray(
        [row["junction_reads"] for row in eligible], dtype=float
    )
    ratios = read / current
    denominator = float(np.dot(current, current))
    slope = float(np.dot(current, read) / denominator) if denominator > 0 else None
    junction_slope = (
        float(np.dot(current, junction_reads) / denominator)
        if denominator > 0 else None
    )
    pearson = float(np.corrcoef(current, read)[0, 1]) if len(eligible) >= 2 else None
    spearman = float(spearmanr(current, read).statistic) if len(eligible) >= 2 else None
    metrics.update({
        "slope_through_origin": safe_metric(slope),
        "junction_reads_per_fitted_N_slope": safe_metric(junction_slope),
        "pearson_r": safe_metric(pearson),
        "spearman_rho": safe_metric(spearman),
        "median_read_to_current_ratio": float(np.median(ratios)),
        "rough_agreement_fraction": float(np.mean(
            (ratios >= ROUGH_RATIO_LOW) & (ratios <= ROUGH_RATIO_HIGH)
        )),
        "within_25_percent_fraction": float(np.mean(np.abs(ratios - 1) <= 0.25)),
        "within_50_percent_fraction": float(np.mean(np.abs(ratios - 1) <= 0.50)),
        "mean_absolute_error_N": float(np.mean(np.abs(read - current))),
        "root_mean_square_error_N": float(np.sqrt(np.mean((read - current) ** 2))),
        "median_absolute_relative_error": float(np.median(np.abs(ratios - 1))),
    })
    return metrics


def render_group(group_id: int, rows: list[dict], output_dir: str, pdf_pages) -> None:
    width = max(7.0, 1.45 * len(rows))
    figure, axis = plt.subplots(figsize=(width, 5.8))
    x = np.arange(len(rows), dtype=float)
    bar_width = 0.36
    current = np.asarray([row["current_depth_N"] for row in rows], dtype=float)
    read = np.asarray([
        np.nan if row["read_depth_N"] is None else row["read_depth_N"]
        for row in rows
    ], dtype=float)
    axis.bar(x - bar_width / 2, current, bar_width, label="Fitted cluster depth (N)")
    axis.bar(x + bar_width / 2, read, bar_width, label="Raw-read-derived depth (N)")
    axis.set_xticks(x)
    axis.set_xticklabels([
        f"{row['nclose_id']}\n{row['label']}"
        for row in rows
    ], rotation=35, ha="right")
    axis.set_ylabel("Depth (N)")
    axis.set_title(
        f"NClose co-carriage group {group_id:03d}: fitted vs raw-read-derived depth"
    )
    axis.grid(axis="y", alpha=0.25)
    axis.legend(loc="upper left")

    finite_values = np.concatenate((current[np.isfinite(current)], read[np.isfinite(read)]))
    ymax = max(float(np.max(finite_values)) if finite_values.size else 1.0, 0.25)
    axis.set_ylim(0, ymax * 1.35)
    for index, row in enumerate(rows):
        ratio_text = (
            "NA" if row["read_to_current_ratio"] is None
            else f"{row['read_to_current_ratio']:.2f}x"
        )
        bar_top = max(
            row["current_depth_N"],
            row["read_depth_N"] if row["read_depth_N"] is not None else 0,
        )
        axis.text(
            index, bar_top + ymax * 0.035,
            f"J={row['junction_reads']}\n{ratio_text}\n{SHORT_STATUS[row['status']]}",
            ha="center", va="bottom", fontsize=7,
            color=("#b2182b" if row["status"] == "DISAGREE_OVER_2X" else "black"),
        )
    figure.tight_layout()
    base = os.path.join(output_dir, f"group_{group_id:03d}")
    figure.savefig(f"{base}.png", dpi=180)
    figure.savefig(f"{base}.pdf")
    pdf_pages.savefig(figure)
    plt.close(figure)


def render_scatter(
    rows: list[dict], metrics: dict, prefix: str, haploid_depth: float
) -> None:
    positive_rows = [
        row for row in rows
        if row["current_depth_N"] > 0
    ]
    plotted_rows = [row for row in positive_rows if scatter_is_included(row)]
    excluded_censat_count = len(positive_rows) - len(plotted_rows)
    finite_rows = [
        row for row in plotted_rows if row["read_depth_N"] is not None
    ]
    missing_rows = [
        row for row in plotted_rows if row["read_depth_N"] is None
    ]
    if not finite_rows and not missing_rows:
        return

    finite_current = np.asarray(
        [row["current_depth_N"] * haploid_depth for row in finite_rows],
        dtype=float,
    )
    finite_read = np.asarray(
        [row["read_depth_N"] * haploid_depth for row in finite_rows],
        dtype=float,
    )
    missing_current = np.asarray(
        [row["current_depth_N"] * haploid_depth for row in missing_rows],
        dtype=float,
    )
    maxima = [1.0]
    if finite_current.size:
        maxima.extend((float(np.max(finite_current)), float(np.max(finite_read))))
    if missing_current.size:
        maxima.append(float(np.max(missing_current)))
    limit = max(maxima) * 1.12
    na_y = -0.075 * limit
    lower_limit = -0.14 * limit if missing_rows else 0.0

    figure, axis = plt.subplots(figsize=(8.0, 7.2))
    axis.plot([0, limit], [0, limit], "--", color="gray", label="y = x", zorder=1)

    for category, style in SCATTER_STYLES.items():
        category_rows = [
            row for row in finite_rows if scatter_category(row) == category
        ]
        if not category_rows:
            continue
        x_values = np.asarray([
            row["current_depth_N"] * haploid_depth for row in category_rows
        ], dtype=float)
        y_values = np.asarray([
            row["read_depth_N"] * haploid_depth for row in category_rows
        ], dtype=float)
        axis.scatter(
            x_values, y_values, s=62, marker=style["marker"],
            color=style["color"], edgecolors="black", linewidths=0.45,
            label=f"{style['label']} (n={len(category_rows)})", zorder=3,
        )
        for row, x_value, y_value in zip(category_rows, x_values, y_values):
            axis.annotate(
                row["nclose_id"].replace("SKYPE.nclose.", "N"),
                (x_value, y_value), xytext=(4, 4), textcoords="offset points",
                fontsize=7,
            )

    if missing_rows:
        axis.axhspan(lower_limit, 0, color="#f0f0f0", zorder=0)
        axis.scatter(
            missing_current, np.full(len(missing_rows), na_y),
            s=70, marker="X", color="#666666", edgecolors="black",
            linewidths=0.45,
            label=f"VAF depth unavailable (n={len(missing_rows)})", zorder=3,
        )
        for row, x_value in zip(missing_rows, missing_current):
            axis.annotate(
                row["nclose_id"].replace("SKYPE.nclose.", "N"),
                (x_value, na_y), xytext=(4, 4), textcoords="offset points",
                fontsize=7,
            )
        axis.text(
            0.015 * limit, na_y, "VAF depth NA", color="#555555",
            ha="left", va="center", fontsize=8,
        )

    assessable = [
        row for row in finite_rows
        if row["status"] in {"AGREE_WITHIN_2X", "DISAGREE_OVER_2X"}
    ]
    assessable_ratios = np.asarray([
        row["read_depth_N"] / row["current_depth_N"] for row in assessable
    ], dtype=float)
    within_25 = (
        float(np.mean(np.abs(assessable_ratios - 1) <= 0.25))
        if assessable_ratios.size else None
    )
    axis.set_xlim(0, limit)
    axis.set_ylim(lower_limit, limit)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("NClose-derived fitted depth (x)")
    axis.set_ylabel("VAF-derived read depth (x)")
    axis.set_title(
        "VAF-derived vs NClose-derived depth (non-CENSAT only)\n"
        f"shown={len(plotted_rows)}/{len(positive_rows)}, "
        f"CENSAT-excluded={excluded_censat_count}, "
        f"finite={len(finite_rows)}, "
        f"assessable={len(assessable)}, "
        f"assessable within 25%="
        f"{'NA' if within_25 is None else f'{within_25:.1%}'}"
    )
    axis.grid(alpha=0.25)
    axis.legend(loc="upper left", fontsize=8)
    figure.tight_layout()
    figure.savefig(os.path.join(prefix, SCATTER_PNG), dpi=180)
    figure.savefig(os.path.join(prefix, SCATTER_PDF))
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot fitted versus raw-read-derived NClose depth by path group"
    )
    parser.add_argument("prefix")
    parser.add_argument("depth_stat_path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    prefix = os.path.abspath(args.prefix)
    if not pipeline_mode_is_karyotype(load_pipeline_mode(prefix)):
        return

    result_path = os.path.join(prefix, COUNT_RESULT_PKL)
    if not os.path.isfile(result_path):
        raise FileNotFoundError(f"Missing cluster NClose read counts: {result_path}")
    with open(result_path, "rb") as handle:
        count_records = pickle.load(handle)
    with open(os.path.join(prefix, "contig_pat_vec_data.pkl"), "rb") as handle:
        paf_ans_list, _key_list, _int2key, _depth_list = pickle.load(handle)

    weights = np.load(os.path.join(prefix, "weight_cluster.npy"))
    path_usage = load_path_usage(prefix, expected_len=len(weights))
    catalog = load_event_catalog(prefix)
    event_by_key = {event["event_key"]: event for event in catalog}
    event_weights = calculate_event_weights(path_usage, weights)
    haploid_depth = median_depth(args.depth_stat_path) / 2.0
    displayed_columns = [
        column
        for column in range(len(paf_ans_list))
        if float(weights[column]) > TARGET_DEPTH_N * haploid_depth
    ]

    record_by_key = {event_key_from_record(record): record for record in count_records}
    components, carriers = connected_components(
        record_by_key, displayed_columns, path_usage
    )
    component_by_key = {
        event_key: group_id
        for group_id, component in enumerate(components, start=1)
        for event_key in component
    }

    rows = []
    for event_key, record in record_by_key.items():
        event = event_by_key.get(event_key)
        if event is None:
            continue
        depth_estimate = record.get("read_count_depth_estimate", {})
        read_depth_n = normalized_read_depth(depth_estimate, haploid_depth)
        current_depth_n = float(event_weights.get(event_key, 0.0)) / haploid_depth
        mapqs = tuple(int(value) for value in record.get("mapq", (0, 0)))
        if len(mapqs) != 2:
            mapqs = (0, 0)
        query_gap = record.get("nclose_query_gap")
        countable = query_gap is not None and int(query_gap) <= 5000
        ratio = (
            read_depth_n / current_depth_n
            if read_depth_n is not None and current_depth_n > 0 else None
        )
        row = {
            "group_id": component_by_key.get(event_key, 0),
            "nclose_id": record.get("nclose_id", ""),
            "event_key": repr(event_key),
            "label": event_label(event),
            "start_chr": event.get("start_chr", ""),
            "start_pos": event.get("start_pos", ""),
            "end_chr": event.get("end_chr", ""),
            "end_pos": event.get("end_pos", ""),
            "mapq_a": mapqs[0],
            "mapq_b": mapqs[1],
            "overlaps_censat": bool(record.get("overlaps_censat", False)),
            "query_gap": query_gap,
            "countable": countable,
            "normal_reads_a": int(record.get("read_counts", {}).get("d1", 0)),
            "junction_reads": int(record.get("read_counts", {}).get("d2", 0)),
            "normal_reads_b": int(record.get("read_counts", {}).get("d3", 0)),
            "local_depth_a": optional_float(depth_estimate.get("chr_a_local_depth")),
            "local_depth_b": optional_float(depth_estimate.get("chr_b_local_depth")),
            "current_depth_N": current_depth_n,
            "read_depth_N": read_depth_n,
            "read_to_current_ratio": ratio,
            "absolute_difference_N": (
                abs(read_depth_n - current_depth_n)
                if read_depth_n is not None else None
            ),
            "carrier_path_columns": ",".join(
                map(str, sorted(carriers.get(event_key, ())))
            ),
        }
        row["status"] = classify_row(row)
        rows.append(row)

    rows.sort(key=lambda row: (row["group_id"], row["nclose_id"]))
    metrics = proportionality_metrics(rows)
    render_scatter(rows, metrics, prefix, haploid_depth)

    output_dir = os.path.join(prefix, GROUP_DIR)
    os.makedirs(output_dir, exist_ok=True)
    for filename in os.listdir(output_dir):
        if filename.startswith("group_") and filename.endswith((".png", ".pdf")):
            os.remove(os.path.join(output_dir, filename))

    rows_by_group = defaultdict(list)
    for row in rows:
        rows_by_group[row["group_id"]].append(row)
    combined_pdf_path = os.path.join(prefix, GROUP_PDF)
    with PdfPages(combined_pdf_path) as pdf_pages:
        for group_id in sorted(rows_by_group):
            render_group(group_id, rows_by_group[group_id], output_dir, pdf_pages)

    fields = [
        "group_id", "nclose_id", "event_key", "label",
        "start_chr", "start_pos", "end_chr", "end_pos",
        "mapq_a", "mapq_b", "overlaps_censat", "query_gap", "countable",
        "normal_reads_a", "junction_reads", "normal_reads_b",
        "local_depth_a", "local_depth_b",
        "current_depth_N", "read_depth_N", "read_to_current_ratio",
        "absolute_difference_N", "status", "carrier_path_columns",
    ]
    with open(os.path.join(prefix, SUMMARY_TSV), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    metrics.update({
        "group_count": len(rows_by_group),
        "group_definition": "connected components of NCloses co-carried by displayed cluster paths",
        "summary_tsv": SUMMARY_TSV,
        "combined_group_pdf": GROUP_PDF,
        "scatter_png": SCATTER_PNG,
        "scatter_pdf": SCATTER_PDF,
    })
    with open(os.path.join(prefix, METRICS_JSON), "w") as handle:
        json.dump(metrics, handle, indent=2)


if __name__ == "__main__":
    main()
