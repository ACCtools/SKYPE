"""Count raw-read support for junctions shown in the cluster Virtual SKY.

This stage prepares the candidate schema consumed by ``03_Anal_bam.py`` and
then invokes that script in NClose-only mode.  Cluster-specific artifact names
keep these reporting counts separate from the optional stage-02 VAF filter.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import pickle
import subprocess
import sys

import numpy as np

from nclose_tracking import (
    load_event_catalog,
    load_filter_status,
    load_path_usage,
    nclose_event_id,
)
from skype_utils import load_pipeline_mode, pipeline_mode_is_karyotype


CTG_NAM = 0
CTG_STR = 2
CTG_END = 3
CTG_DIR = 4
CHR_NAM = 5
CHR_STR = 7
CHR_END = 8
CTG_CENSAT = 18

TARGET_DEPTH_N = 0.2
JULIA_BAM_THREAD_LIMIT = 32
CLUSTER_ARTIFACT_PREFIX = "cluster_"
PRECLUSTER_ARTIFACT_PREFIX = "precluster_"
NCLOSE_COUNT_CANDIDATE_PKL = "nclose_count_candidates.pkl"
NCLOSE_COUNT_RESULT_PKL = "nclose_count_result.pkl"
NCLOSE_COUNT_REPORT_TSV = "nclose_read_counts.tsv"


def read_ppc_data(path: str) -> list[list]:
    rows = []
    numeric_fields = (1, 2, 3, 6, 7, 8, 9)
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            row = line.rstrip("\r\n").split("\t")
            if len(row) <= CHR_END:
                raise ValueError(f"Malformed PPC PAF row {path}:{line_number}")
            for field in numeric_fields:
                row[field] = int(row[field])
            rows.append(row)
    return rows


def median_depth(path: str) -> float:
    opener = gzip.open if str(path).endswith(".gz") else open
    values = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 8 or fields[0] == "chrM":
                continue
            values.append(float(fields[7]))
    if not values:
        raise ValueError(f"No autosomal/sex-chromosome depth rows found in {path}")
    return float(np.median(values))


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


def breakend_coord(contig: list, side_index: int) -> int:
    nclose_at_start = side_index == 0
    is_forward = contig[CTG_DIR] == "+"
    return int(contig[CHR_STR if is_forward ^ nclose_at_start else CHR_END])


def expected_high_side(contig: list, side_index: int) -> str:
    nclose_at_start = side_index == 0
    is_forward = contig[CTG_DIR] == "+"
    return "right" if is_forward ^ nclose_at_start else "left"


def raw_endpoint(contig: list, side_index: int, contig_index: int) -> dict:
    nclose_at_start = side_index == 0
    is_forward = contig[CTG_DIR] == "+"
    return {
        "chrom": contig[CHR_NAM],
        "coord": breakend_coord(contig, side_index),
        "dir": contig[CTG_DIR],
        "match_forward": bool(is_forward == nclose_at_start),
        "contig_idx": int(contig_index),
        "ctg_name": contig[CTG_NAM],
        "ctg_st": int(contig[CTG_STR]),
        "ctg_nd": int(contig[CTG_END]),
        "ref_st": int(contig[CHR_STR]),
        "ref_nd": int(contig[CHR_END]),
        "expected_high_side": expected_high_side(contig, side_index),
    }


def canonical_raw_nclose_layout(contig_data: list[list], nclose_pair) -> dict:
    ordered_endpoints = [
        raw_endpoint(contig_data[int(contig_index)], side_index, int(contig_index))
        for side_index, contig_index in enumerate(nclose_pair)
    ]
    order = [0, 1]
    keys = [
        (chromosome_sort_key(ordered_endpoints[index]["chrom"]),
         ordered_endpoints[index]["coord"])
        for index in order
    ]
    if keys[1] < keys[0]:
        order = [1, 0]
    endpoints = [ordered_endpoints[index] for index in order]
    return {
        "nclose_key": tuple(sorted(map(int, nclose_pair))),
        "chroms": tuple(endpoint["chrom"] for endpoint in endpoints),
        "coords": tuple(endpoint["coord"] for endpoint in endpoints),
        "sides": tuple(endpoint["expected_high_side"] for endpoint in endpoints),
        "endpoints": tuple(endpoints),
        "ordered_endpoints": tuple(ordered_endpoints),
    }


def raw_count_node_pair(event: dict) -> tuple[int, int] | None:
    if event.get("kind") == "bnd":
        event_key = event.get("event_key", ())
        if len(event_key) == 2 and all(isinstance(index, int) for index in event_key):
            return tuple(sorted(map(int, event_key)))
    if event.get("graph_only") and event.get("kind") == "indel":
        type4_tuple = event.get("type4_tuple", ())
        if len(type4_tuple) == 2:
            return tuple(sorted(map(int, type4_tuple)))
    return None


def selected_junction_keys(
    prefix: str,
    depth_stat_path: str,
    event_catalog: list[dict],
    selection_stage: str = "cluster",
    target_depth_n: float = TARGET_DEPTH_N,
) -> tuple[set, list[int], float]:
    with open(os.path.join(prefix, "contig_pat_vec_data.pkl"), "rb") as handle:
        paf_ans_list, _key_list, _int2key, _depth_list = pickle.load(handle)
    weight_filename = "weight_cluster.npy" if selection_stage == "cluster" else "weight.npy"
    weights = np.load(os.path.join(prefix, weight_filename))
    path_usage = load_path_usage(prefix, expected_len=len(weights))
    mean_depth = median_depth(depth_stat_path)
    if selection_stage == "cluster":
        minimum_weight = float(target_depth_n) * mean_depth / 2.0
        selected_columns = [
            column_index
            for column_index in range(len(paf_ans_list))
            if float(weights[column_index]) > minimum_weight
        ]
    elif selection_stage == "base":
        filter_status = load_filter_status(prefix)
        if "base" not in filter_status.get("stages", {}):
            raise ValueError("Stage-23 base NClose filter status is missing")
        selected_columns = [
            int(column_index)
            for column_index in filter_status["stages"]["base"]["active_columns"]
            if int(column_index) < len(paf_ans_list)
        ]
    else:
        raise ValueError(f"Unsupported NClose selection stage: {selection_stage}")
    junction_keys = {
        event["event_key"]
        for event in event_catalog
        if raw_count_node_pair(event) is not None
    }
    displayed_keys = {
        event_key
        for column_index in selected_columns
        for event_key in path_usage[column_index]
        if event_key in junction_keys
    }
    return displayed_keys, selected_columns, mean_depth


def build_cluster_nclose_candidates(
    contig_data: list[list],
    event_catalog: list[dict],
    displayed_keys: set,
) -> list[dict]:
    candidates = []
    for catalog_index, event in enumerate(event_catalog):
        event_key = event["event_key"]
        if event_key not in displayed_keys:
            continue
        nclose_key = raw_count_node_pair(event)
        if nclose_key is None:
            continue
        if max(nclose_key) >= len(contig_data):
            raise IndexError(
                f"NClose {event_key} references PPC row outside 0..{len(contig_data) - 1}"
            )
        layout = canonical_raw_nclose_layout(contig_data, nclose_key)
        overlaps_censat = any(
            len(contig_data[index]) > CTG_CENSAT
            and contig_data[index][CTG_CENSAT] != "0"
            for index in nclose_key
        )
        candidates.append({
            "pair_id": len(candidates),
            "nclose_id": nclose_event_id(catalog_index + 1),
            "event_key": event_key,
            "event_kind": event.get("kind"),
            "nclose_key": nclose_key,
            "ctg_name": event.get("contig_name", layout["ordered_endpoints"][0]["ctg_name"]),
            "layout": layout,
            "overlaps_censat": overlaps_censat,
            "mapq": tuple(int(contig_data[index][9]) for index in nclose_key),
        })
    missing = displayed_keys - {candidate["event_key"] for candidate in candidates}
    if missing:
        raise ValueError(f"Cluster NClose candidates missing catalog events: {sorted(missing)!r}")
    return candidates


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Count raw-read support for NCloses displayed by cluster Virtual SKY"
    )
    parser.add_argument("ppc_paf_file_path")
    parser.add_argument("depth_stat_path")
    parser.add_argument("reference_fai_path")
    parser.add_argument("prefix")
    parser.add_argument("read_bam_loc")
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument(
        "--selection_stage", choices=("base", "cluster"), default="cluster"
    )
    parser.add_argument("--artifact_prefix", default=None)
    parser.add_argument("--progress", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    prefix = os.path.abspath(args.prefix)
    pipeline_mode = load_pipeline_mode(prefix)
    if not pipeline_mode_is_karyotype(pipeline_mode):
        logging.info("Cluster NClose raw-read counting skipped outside karyotype mode")
        return

    event_catalog = load_event_catalog(prefix)
    artifact_prefix = args.artifact_prefix or (
        PRECLUSTER_ARTIFACT_PREFIX
        if args.selection_stage == "base"
        else CLUSTER_ARTIFACT_PREFIX
    )
    if artifact_prefix and not all(
        char.isalnum() or char in {"_", "-", "."} for char in artifact_prefix
    ):
        raise ValueError("Artifact prefix contains unsupported characters")

    displayed_keys, displayed_columns, mean_depth = selected_junction_keys(
        prefix, args.depth_stat_path, event_catalog,
        selection_stage=args.selection_stage,
    )
    contig_data = read_ppc_data(args.ppc_paf_file_path)
    candidates = build_cluster_nclose_candidates(
        contig_data, event_catalog, displayed_keys
    )
    candidate_path = os.path.join(
        prefix, f"{artifact_prefix}{NCLOSE_COUNT_CANDIDATE_PKL}"
    )
    with open(candidate_path, "wb") as handle:
        pickle.dump(candidates, handle)
    for stale_name in (NCLOSE_COUNT_RESULT_PKL, NCLOSE_COUNT_REPORT_TSV):
        stale_path = os.path.join(prefix, f"{artifact_prefix}{stale_name}")
        if os.path.isfile(stale_path):
            os.remove(stale_path)

    logging.info(
        "%s junction raw-count candidates: %d events from %d selected paths "
        "(median depth %.6f)",
        args.selection_stage, len(candidates), len(displayed_columns), mean_depth,
    )
    julia_threads = max(1, min(int(args.thread), JULIA_BAM_THREAD_LIMIT))
    command = [
        sys.executable,
        "-X", f"juliacall-threads={julia_threads}",
        "-X", "juliacall-handle-signals=yes",
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "03_Anal_bam.py"),
        prefix,
        os.path.abspath(args.read_bam_loc),
        os.path.abspath(args.reference_fai_path),
        os.path.abspath(args.depth_stat_path),
        "--nclose_count_only",
        "--nclose_count_artifact_prefix", artifact_prefix,
    ]
    if args.progress:
        command.append("--progress")
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
