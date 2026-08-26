"""Breakend-graph construction and path search for native SKYPE stage 10.

The stage-01 handoff deliberately contains only the final contig-node snapshot,
path-ordered NClose pairs, and telomere connections.  Everything else used by
the graph is either derived here or loaded by the stage-10 CLI from fixed
reference inputs.
"""

from __future__ import annotations

import copy
import json
import logging
import os
import pickle
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path
from typing import Mapping, Sequence

import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

from skype_utils import (
    CHROMOSOME_COUNT,
    K,
    M,
    VCF_TYPE4_GRAPH_NODE_PREFIX,
)


# Final preprocessed-contig row schema.
CTG_NAM = 0
CTG_LEN = 1
CTG_STR = 2
CTG_END = 3
CTG_DIR = 4
CHR_NAM = 5
CHR_LEN = 6
CHR_STR = 7
CHR_END = 8
CTG_MAPQ = 9
CTG_TYP = 10
CTG_STRND = 11
CTG_ENDND = 12
CTG_TELCHR = 13
CTG_TELDIR = 14
CTG_TELCON = 15
CTG_RPTCHR = 16
CTG_RPTCASE = 17
CTG_CENSAT = 18
CTG_MAINFLOWDIR = 19
CTG_MAINFLOWCHR = 20
CTG_GLOBALIDX = 21

DIR_FOR = 1
DIR_BAK = 0
DIR_IN = 3
DIR_OUT = 2

DIR_CHANGE_LIMIT_ABS_MAX = 1
CHR_CHANGE_LIMIT_HARD_START = 2
HARD_GRAPH_NCLOSE_COUNT = 1000
FAIL_NCLOSE_COUNT = 10000

TOT_PATH_LIMIT = 3 * M
PAT_PATH_LIMIT = 10 * K
BND_OVERUSE_CNT = 2
CENSAT_VISIT_LIMIT = 2
PATH_MAJOR_COMPONENT = 3
PATH_COMPRESS_LIMIT = 50 * K
IGNORE_PATH_LIMIT = 50 * K
MIN_PATH_REF_LEN = 5 * M

TYPE4_INDEL_GRAPH_MIN_SPAN = 1 * M
TYPE4_INDEL_GRAPH_DEPTH_WINDOW = 500 * K
TYPE4_INDEL_GRAPH_DEPTH_DIFF_RATIO = 0.2
CENSAT_COMPRESSABLE_THRESHOLD = 1000 * K

LIMIT_COMBINATIONS_JSON = "limit_combinations.json"


@dataclass(frozen=True)
class Stage10Input:
    """The complete and intentionally small stage-01 to stage-10 contract."""

    contig_data: list
    nclose_nodes: dict
    telo_contig: dict


@dataclass(frozen=True)
class GraphBuildState:
    """Graph-only state derived from one immutable stage-01 handoff."""

    graph_nclose_nodes: Mapping
    selected_type4_indel_graph_edges: tuple
    type4_indel_zero_dim_edge_set: frozenset
    bnd_graph_adjacency: Mapping


@dataclass(frozen=True)
class PathSearchConfig:
    thread: int
    graph_depth: int = 4
    fixed_limit_combination: tuple[int, int] | None = None
    progress: bool = False
    verbose: bool = False
    raw_output_dir: str | None = None
    total_path_limit: int = TOT_PATH_LIMIT
    per_pair_path_limit: int = PAT_PATH_LIMIT


@dataclass
class PathSearchResult:
    chr_change_limit: int
    dir_change_limit: int
    path_items: list
    counts: list

    def as_path_dict(self) -> dict:
        return {key: paths for key, paths in self.path_items}


@dataclass(frozen=True)
class _PathWorkerContext:
    graph: nx.DiGraph
    contig_data: Sequence
    contig_data_size: int
    chr_rev_corr: Mapping[int, str]
    chr_len: Mapping[str, int]
    chr_change_limit: int
    dir_change_limit: int
    per_pair_path_limit: int
    verbose: bool
    raw_output_dir: str | None


_WORKER_CONTEXT: _PathWorkerContext | None = None


def load_stage10_input(path: str | os.PathLike) -> Stage10Input:
    """Load the named three-field handoff written by stage 01."""

    with open(path, "rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, Mapping):
        raise ValueError(f"Stage-10 input must be a named mapping: {path}")

    required_fields = {
        "contig_data",
        "nclose_nodes",
        "telo_contig",
    }
    missing = required_fields - set(payload)
    if missing:
        raise ValueError(
            f"Stage-10 input is missing fields {sorted(missing)}: {path}"
        )
    unexpected = set(payload) - required_fields
    if unexpected:
        raise ValueError(
            f"Stage-10 input has unexpected fields {sorted(unexpected)}: {path}"
        )

    return Stage10Input(
        contig_data=payload["contig_data"],
        nclose_nodes=payload["nclose_nodes"],
        telo_contig=payload["telo_contig"],
    )


def save_stage10_input(
    path: str | os.PathLike,
    contig_data,
    nclose_nodes,
    telo_contig,
) -> None:
    """Write the stage-01 handoff without adding hidden graph state."""

    with open(path, "wb") as handle:
        pickle.dump(
            {
                "contig_data": contig_data,
                "nclose_nodes": nclose_nodes,
                "telo_contig": telo_contig,
            },
            handle,
        )


def load_chr_lengths(path: str | os.PathLike) -> dict[str, int]:
    chr_len = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t")
            chr_len[fields[0]] = int(fields[1])
    return chr_len


def load_censat_intervals(path: str | os.PathLike) -> dict[str, list]:
    intervals = defaultdict(list)
    with open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t")
            interval = (int(fields[1]), int(fields[2]))
            if interval[1] - interval[0] > CENSAT_COMPRESSABLE_THRESHOLD:
                intervals[fields[0]].append(interval)
    return intervals


def load_depth_dataframe(path: str | os.PathLike) -> pd.DataFrame:
    return pd.read_csv(
        path,
        compression="gzip",
        comment="#",
        sep="\t",
        names=[
            "chr",
            "st",
            "nd",
            "length",
            "covsite",
            "totaldepth",
            "cov",
            "meandepth",
        ],
    ).query('chr != "chrM"')


def load_limit_combinations(path: str | os.PathLike) -> tuple[int, int]:
    try:
        with open(path, "rt", encoding="utf-8") as handle:
            data = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"could not read {path}: {exc}") from exc

    combination = data.get("limit_combinations") if isinstance(data, dict) else None
    if (
        not isinstance(combination, list)
        or len(combination) != 2
        or any(type(value) is not int for value in combination)
    ):
        raise ValueError(
            f"{path} must contain an integer pair named 'limit_combinations'"
        )

    chr_limit, dir_limit = combination
    if chr_limit < 1 or dir_limit not in (0, 1):
        raise ValueError(
            f"invalid limit_combinations in {path}: {combination!r}"
        )
    return chr_limit, dir_limit


def write_limit_combinations(
    path: str | os.PathLike,
    combination: tuple[int, int],
) -> None:
    temporary = f"{path}.tmp-{os.getpid()}"
    with open(temporary, "wt", encoding="utf-8") as handle:
        json.dump({"limit_combinations": list(combination)}, handle, indent=2)
        handle.write("\n")
    os.replace(temporary, path)


def chr_correlation_maker(contig_data: Sequence) -> tuple[dict, dict]:
    """Recreate the legacy chromosome-terminal index layout."""

    chr_corr = {}
    chr_rev_corr = {}
    contig_data_size = len(contig_data)
    for index in range(1, CHROMOSOME_COUNT):
        terminal_idx = contig_data_size + index - 1
        chr_corr[f"chr{index}f"] = terminal_idx
        chr_rev_corr[terminal_idx] = f"chr{index}f"
    chr_corr["chrXf"] = contig_data_size + CHROMOSOME_COUNT - 1
    chr_corr["chrYf"] = contig_data_size + CHROMOSOME_COUNT - 1
    chr_rev_corr[contig_data_size + CHROMOSOME_COUNT - 1] = "chrXf"

    for index in range(1, CHROMOSOME_COUNT):
        terminal_idx = contig_data_size + CHROMOSOME_COUNT + index - 1
        chr_corr[f"chr{index}b"] = terminal_idx
        chr_rev_corr[terminal_idx] = f"chr{index}b"
    chr_corr["chrXb"] = contig_data_size + 2 * CHROMOSOME_COUNT - 1
    chr_corr["chrYb"] = contig_data_size + 2 * CHROMOSOME_COUNT - 1
    chr_rev_corr[contig_data_size + 2 * CHROMOSOME_COUNT - 1] = "chrXb"
    return chr_corr, chr_rev_corr


def count_nclose_nodes(nclose_nodes: Mapping) -> int:
    return 2 * sum(len(pair_list) for pair_list in nclose_nodes.values())


def distance_checker(node_a: Sequence, node_b: Sequence) -> int:
    if max(int(node_a[CHR_STR]), int(node_b[CHR_STR])) < min(
        int(node_a[CHR_END]), int(node_b[CHR_END])
    ):
        return 0
    return min(
        abs(int(node_b[CHR_STR]) - int(node_a[CHR_END])),
        abs(int(node_b[CHR_END]) - int(node_a[CHR_STR])),
    )


def overlap_calculator(node_a: Sequence, node_b: Sequence) -> int:
    return min(
        abs(int(node_a[CHR_END]) - int(node_b[CHR_STR])),
        abs(int(node_b[CHR_END]) - int(node_a[CHR_STR])),
    )


def node_is_censat(node: Sequence) -> bool:
    return node[CTG_CENSAT] != "0"


def calculate_single_contig_ref_ratio(contig_data: Sequence) -> tuple[float, int]:
    total_ref_len = sum(
        int(node[CHR_END]) - int(node[CHR_STR]) for node in contig_data
    )
    if contig_data[0][CTG_DIR] == "+":
        estimated_ref_len = (
            int(contig_data[-1][CHR_END]) - int(contig_data[0][CHR_STR])
        )
    else:
        estimated_ref_len = (
            int(contig_data[0][CHR_END]) - int(contig_data[-1][CHR_STR])
        )
    return estimated_ref_len / total_ref_len, total_ref_len


def iter_nclose_owner_pairs(nclose_source):
    for contig_name, pair_list in nclose_source.items():
        for pair in pair_list:
            yield contig_name, tuple(pair)


def get_breakend_coord(contig: Sequence, side_idx: int) -> int:
    nclose_loc = side_idx == 0
    nclose_dir = contig[CTG_DIR] == "+"
    return int(contig[CHR_STR if nclose_dir ^ nclose_loc else CHR_END])


def canonical_nclose_snapshot(nclose_nodes: Mapping) -> dict:
    return {
        key: tuple(
            tuple(int(node_idx) for node_idx in pair) for pair in pair_list
        )
        for key, pair_list in nclose_nodes.items()
    }


def assert_vcf_nclose_has_no_indel_like(
    contig_data: Sequence,
    nclose_nodes: Mapping,
    context: str,
) -> None:
    offenders = []
    for event_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            if len(pair) != 2:
                raise AssertionError(
                    f"{context}: invalid canonical VCF nclose pair for "
                    f"{event_name}: {pair}"
                )
            node_a = contig_data[int(pair[0])]
            node_b = contig_data[int(pair[1])]
            if (
                node_a[CHR_NAM] == node_b[CHR_NAM]
                and node_a[CTG_DIR] == node_b[CTG_DIR]
            ):
                offenders.append(
                    (
                        event_name,
                        tuple(int(node_idx) for node_idx in pair),
                        node_a[CHR_NAM],
                        node_a[CTG_DIR],
                    )
                )

    if offenders:
        preview = ", ".join(
            f"{name}:{pair}:{chrom}{direction}/{chrom}{direction}"
            for name, pair, chrom, direction in offenders[:8]
        )
        raise AssertionError(
            f"{context}: canonical VCF nclose contains {len(offenders)} "
            f"Indel-like same-chrom/same-dir pair(s): {preview}"
        )


def initialize_bnd_graph(
    contig_data: Sequence,
    nclose_nodes: Mapping,
    telo_contig: Mapping,
    contig_data_size: int,
    chr_rev_corr: Mapping[int, str],
) -> dict:
    """Build the directed NClose/telomere adjacency used by path search."""

    telo_dir_map = {
        edge[1]: telo_name[-1]
        for telo_name, edge_list in telo_contig.items()
        for edge in edge_list
    }

    def is_telo_dir_consistent(telo_idx, approach_dir):
        telo_side = telo_dir_map.get(telo_idx)
        if telo_side is None:
            return True
        return (telo_side == "b" and approach_dir == "dec") or (
            telo_side == "f" and approach_dir == "inc"
        )

    adjacency = defaultdict(list)
    for telo_name, edge_list in telo_contig.items():
        for edge in edge_list:
            adjacency[telo_name].append(edge[:2])
            adjacency[(DIR_IN, edge[1])].append(telo_name)
            adjacency[telo_name].append((DIR_OUT, edge[1]))

    # Every pair under one owner is a separate branching NClose edge.
    for pair_list in nclose_nodes.values():
        for node_a_idx, node_b_idx in pair_list:
            adjacency[(DIR_FOR, node_a_idx)].append([DIR_FOR, node_b_idx])
            adjacency[(DIR_BAK, node_b_idx)].append([DIR_BAK, node_a_idx])

    owner_names = list(nclose_nodes)
    for owner_a_idx in range(len(owner_names)):
        for owner_b_idx in range(owner_a_idx + 1, len(owner_names)):
            for pair_a in nclose_nodes[owner_names[owner_a_idx]]:
                for pair_b in nclose_nodes[owner_names[owner_b_idx]]:
                    if len(pair_a) != 2 or len(pair_b) != 2:
                        raise ValueError(
                            f"Invalid NClose pair: {pair_a!r}, {pair_b!r}"
                        )
                    for endpoint_a in range(2):
                        for endpoint_b in range(2):
                            idx_a = pair_a[endpoint_a]
                            idx_b = pair_b[endpoint_b]
                            contig_a = contig_data[idx_a]
                            contig_b = contig_data[idx_b]
                            if contig_a[CHR_NAM] != contig_b[CHR_NAM]:
                                continue

                            end_a = ("b" if endpoint_a % 2 else "f") + contig_a[CTG_DIR]
                            end_b = ("b" if endpoint_b % 2 else "f") + contig_b[CTG_DIR]
                            a_st, a_end = contig_a[CHR_STR], contig_a[CHR_END]
                            b_st, b_end = contig_b[CHR_STR], contig_b[CHR_END]
                            overlaps = distance_checker(contig_a, contig_b) == 0

                            if end_a == "f+":
                                if end_b == "f-" and (a_st > b_st or overlaps):
                                    adjacency[(DIR_BAK, idx_a)].append([DIR_FOR, idx_b])
                                    adjacency[(DIR_BAK, idx_b)].append([DIR_FOR, idx_a])
                                elif end_b == "b+" and (a_st > b_end or overlaps):
                                    adjacency[(DIR_FOR, idx_b)].append([DIR_FOR, idx_a])
                                    adjacency[(DIR_BAK, idx_a)].append([DIR_BAK, idx_b])
                            elif end_a == "b+":
                                if end_b == "b-" and (b_end > a_end or overlaps):
                                    adjacency[(DIR_FOR, idx_a)].append([DIR_BAK, idx_b])
                                    adjacency[(DIR_FOR, idx_b)].append([DIR_BAK, idx_a])
                                elif end_b == "f+" and (b_st > a_end or overlaps):
                                    adjacency[(DIR_BAK, idx_b)].append([DIR_BAK, idx_a])
                                    adjacency[(DIR_FOR, idx_a)].append([DIR_FOR, idx_b])
                            elif end_a == "f-":
                                if end_b == "f+" and (b_st > a_st or overlaps):
                                    adjacency[(DIR_BAK, idx_a)].append([DIR_FOR, idx_b])
                                    adjacency[(DIR_BAK, idx_b)].append([DIR_FOR, idx_a])
                                elif end_b == "b-" and (b_end > a_st or overlaps):
                                    adjacency[(DIR_FOR, idx_b)].append([DIR_FOR, idx_a])
                                    adjacency[(DIR_BAK, idx_a)].append([DIR_BAK, idx_b])
                            elif end_a == "b-":
                                if end_b == "b+" and (a_end > b_end or overlaps):
                                    adjacency[(DIR_FOR, idx_a)].append([DIR_BAK, idx_b])
                                    adjacency[(DIR_FOR, idx_b)].append([DIR_BAK, idx_a])
                                elif end_b == "f-" and (a_end > b_st or overlaps):
                                    adjacency[(DIR_BAK, idx_b)].append([DIR_BAK, idx_a])
                                    adjacency[(DIR_FOR, idx_a)].append([DIR_FOR, idx_b])

    for telo_name, telo_edges in telo_contig.items():
        for telo_edge in telo_edges:
            telo_idx = telo_edge[1]
            telo_node = contig_data[telo_idx]
            for pair_list in nclose_nodes.values():
                for pair in pair_list:
                    if telo_idx in pair:
                        continue
                    for endpoint_idx, nclose_idx in enumerate(pair):
                        nclose_node = contig_data[nclose_idx]
                        if nclose_node[CHR_NAM] != telo_node[CHR_NAM]:
                            continue

                        nclose_end = (
                            ("b" if endpoint_idx % 2 else "f")
                            + nclose_node[CTG_DIR]
                        )
                        telo_front = telo_name[-1] == "f"
                        if telo_node[CTG_DIR] == "+":
                            telo_end = "f" if telo_front else "b"
                        else:
                            telo_end = "b" if telo_front else "f"
                        telo_end += telo_node[CTG_DIR]

                        dir1 = (
                            "inc"
                            if telo_node[CHR_END] <= nclose_node[CHR_END]
                            else "dec"
                        )
                        dir2 = (
                            "inc"
                            if telo_node[CHR_STR] <= nclose_node[CHR_STR]
                            else "dec"
                        )
                        if dir1 != dir2 or not is_telo_dir_consistent(telo_idx, dir1):
                            continue

                        if dir1 == "inc":
                            if telo_node[CTG_DIR] == "+" and nclose_end == "f+":
                                adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                            elif telo_node[CTG_DIR] == "-" and nclose_end == "b-":
                                adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                                adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                            elif telo_node[CTG_DIR] == "+" and nclose_end == "b-":
                                adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                            elif telo_node[CTG_DIR] == "-" and nclose_end == "f+":
                                adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                        else:
                            if telo_node[CTG_DIR] == "+" and nclose_end == "b+":
                                adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                                adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                            elif telo_node[CTG_DIR] == "-" and nclose_end == "f-":
                                adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                            elif telo_node[CTG_DIR] == "+" and nclose_end == "f-":
                                adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                            elif telo_node[CTG_DIR] == "-" and nclose_end == "b+":
                                adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])

    # Join the two telomere anchors belonging to each chromosome.
    for front_idx in range(
        contig_data_size,
        contig_data_size + CHROMOSOME_COUNT,
    ):
        back_idx = front_idx + CHROMOSOME_COUNT
        front_name = chr_rev_corr[front_idx]
        back_name = chr_rev_corr[back_idx]
        for node_a in telo_contig.get(front_name, ()):
            for node_b in telo_contig.get(back_name, ()):
                adjacency[(DIR_OUT, node_a[1])].append([DIR_IN, node_b[1]])
                adjacency[(DIR_OUT, node_b[1])].append([DIR_IN, node_a[1]])

    return adjacency


def collect_nclose_breakends(contig_data: Sequence, nclose_source: Mapping) -> dict:
    breakends_by_chrom = defaultdict(list)
    for contig_name, pair in iter_nclose_owner_pairs(nclose_source):
        pair_key = tuple(pair)
        for side_idx, contig_idx in enumerate(pair):
            contig = contig_data[contig_idx]
            breakends_by_chrom[contig[CHR_NAM]].append(
                (get_breakend_coord(contig, side_idx), pair_key, contig_name)
            )
    for chrom in breakends_by_chrom:
        breakends_by_chrom[chrom].sort(key=lambda value: value[0])
    return breakends_by_chrom


def get_graph_edge_weight(contig_data: Sequence, src_idx: int, dst_idx: int) -> int:
    contig_s = contig_data[src_idx]
    contig_e = contig_data[dst_idx]
    if contig_s[CTG_NAM] == contig_e[CTG_NAM]:
        if src_idx < dst_idx:
            idx_range = range(src_idx + 1, dst_idx)
        else:
            idx_range = range(dst_idx + 1, src_idx)
        return sum(
            int(contig_data[index][CHR_END])
            - int(contig_data[index][CHR_STR])
            for index in idx_range
        )

    source_len = int(contig_s[CHR_END]) - int(contig_s[CHR_STR])
    target_len = int(contig_e[CHR_END]) - int(contig_e[CHR_STR])
    if distance_checker(contig_s, contig_e) == 0:
        return source_len + target_len - overlap_calculator(contig_s, contig_e)
    return distance_checker(contig_s, contig_e) + source_len + target_len


def get_type4_breakend_span(contig_slice: Sequence) -> tuple[int, int]:
    return (
        get_breakend_coord(contig_slice[0], 0),
        get_breakend_coord(contig_slice[-1], 1),
    )


def get_type4_outer_ref_span(contig_slice: Sequence) -> tuple[int, int]:
    if contig_slice[0][CTG_DIR] == "+":
        return (
            int(contig_slice[0][CHR_STR]),
            int(contig_slice[-1][CHR_END]),
        )
    return (
        int(contig_slice[-1][CHR_STR]),
        int(contig_slice[0][CHR_END]),
    )


def point_overlaps_censat(censat_dict: Mapping, chrom: str, coord: int) -> bool:
    return any(
        max(start, coord) < min(end, coord + 1)
        for start, end in censat_dict.get(chrom, ())
    )


def span_overlaps_censat(
    censat_dict: Mapping,
    chrom: str,
    span_start: int,
    span_end: int,
) -> bool:
    return any(
        max(start, span_start) < min(end, span_end)
        for start, end in censat_dict.get(chrom, ())
    )


def make_depth_by_chrom(depth_df: pd.DataFrame) -> dict:
    return {
        chrom: chrom_df
        for chrom, chrom_df in depth_df.groupby("chr", sort=False)
    }


def weighted_depth_mean(
    depth_by_chrom: Mapping,
    chrom: str,
    start: int,
    end: int,
) -> float | None:
    if end <= start:
        return None
    chrom_df = depth_by_chrom.get(chrom)
    if chrom_df is None:
        return None
    mask = (chrom_df["nd"] > start) & (chrom_df["st"] < end)
    if not mask.any():
        return None

    rows = chrom_df.loc[mask]
    overlap_start = np.maximum(rows["st"].to_numpy(), start)
    overlap_end = np.minimum(rows["nd"].to_numpy(), end)
    weights = overlap_end - overlap_start
    if weights.sum() <= 0:
        return None
    return float(np.average(rows["meandepth"].to_numpy(), weights=weights))


def expected_depth_step(
    depth_by_chrom: Mapping,
    chrom: str,
    coord: int,
    expected_high_side: str,
    window: int = TYPE4_INDEL_GRAPH_DEPTH_WINDOW,
) -> float | None:
    left_mean = weighted_depth_mean(depth_by_chrom, chrom, coord - window, coord)
    right_mean = weighted_depth_mean(depth_by_chrom, chrom, coord, coord + window)
    if left_mean is None or right_mean is None:
        return None
    if expected_high_side == "left":
        return left_mean - right_mean
    if expected_high_side == "right":
        return right_mean - left_mean
    return None


def type4_indel_expected_high_sides(indel_kind: str) -> tuple[str, str]:
    if indel_kind == "deletion":
        return "left", "right"
    if indel_kind == "insertion":
        return "right", "left"
    raise ValueError(f"Unknown type4 indel kind: {indel_kind}")


def type4_indel_depth_supported(
    depth_by_chrom: Mapping,
    censat_dict: Mapping,
    chrom: str,
    span_start: int,
    span_end: int,
    indel_kind: str,
    normal_haploid_depth: float,
) -> tuple[bool, list, str]:
    if point_overlaps_censat(censat_dict, chrom, span_start) or point_overlaps_censat(
        censat_dict, chrom, span_end
    ):
        return False, [], "endpoint_censat"

    expected_high_sides = type4_indel_expected_high_sides(indel_kind)
    steps = [
        expected_depth_step(
            depth_by_chrom,
            chrom,
            span_start,
            expected_high_sides[0],
        ),
        expected_depth_step(
            depth_by_chrom,
            chrom,
            span_end,
            expected_high_sides[1],
        ),
    ]
    if any(step is None for step in steps):
        return False, steps, "missing_depth_window"

    min_step = TYPE4_INDEL_GRAPH_DEPTH_DIFF_RATIO * normal_haploid_depth
    if all(step >= min_step for step in steps):
        return True, steps, "pass"
    return False, steps, "depth_direction_or_size"


def collect_type4_indel_graph_candidates(
    contig_data: Sequence,
    depth_df: pd.DataFrame,
    censat_dict: Mapping,
    min_span: int = TYPE4_INDEL_GRAPH_MIN_SPAN,
) -> list:
    candidates = []
    depth_by_chrom = make_depth_by_chrom(depth_df)
    normal_haploid_depth = float(np.median(depth_df["meandepth"])) / 2
    start_idx = 0
    while start_idx < len(contig_data):
        end_idx = int(contig_data[start_idx][CTG_ENDND])
        contig_slice = contig_data[start_idx : end_idx + 1]
        contig_start = contig_data[start_idx]
        contig_end = contig_data[end_idx]

        if (
            int(contig_start[CTG_TYP]) == 4
            and contig_start[CHR_NAM] == contig_end[CHR_NAM]
            and contig_start[CTG_DIR] == contig_end[CTG_DIR]
        ):
            ratio, total_ref_len = calculate_single_contig_ref_ratio(contig_slice)
            ref_a, ref_b = get_type4_breakend_span(contig_slice)
            outer_ref_a, outer_ref_b = get_type4_outer_ref_span(contig_slice)
            span_start, span_end = sorted((ref_a, ref_b))
            outer_span_start, outer_span_end = sorted((outer_ref_a, outer_ref_b))
            span_len = span_end - span_start

            indel_kind = None
            if ratio < 0:
                indel_kind = "insertion"
            elif ratio > 0:
                indel_kind = "deletion"

            if indel_kind is not None and span_len >= min_span:
                censat_jump = span_overlaps_censat(
                    censat_dict,
                    contig_start[CHR_NAM],
                    span_start,
                    span_end,
                )
                # Preserve the legacy assembly rule: graph-only deletions are
                # retained only when their jump crosses CEN-SAT.
                if indel_kind == "deletion" and not censat_jump:
                    start_idx = end_idx + 1
                    continue

                depth_pass, depth_steps, depth_filter_reason = (
                    type4_indel_depth_supported(
                        depth_by_chrom,
                        censat_dict,
                        contig_start[CHR_NAM],
                        span_start,
                        span_end,
                        indel_kind,
                        normal_haploid_depth,
                    )
                )
                if not depth_pass:
                    start_idx = end_idx + 1
                    continue

                candidates.append(
                    {
                        "type4_tuple": (start_idx, end_idx),
                        "src": (DIR_FOR, start_idx),
                        "dst": (DIR_FOR, end_idx),
                        "type4_kind": f"raw_type4_{indel_kind}",
                        "indel_kind": indel_kind,
                        "chrom": contig_start[CHR_NAM],
                        "base_chrom": contig_start[CHR_NAM],
                        "base_st": span_start,
                        "base_nd": span_end,
                        "outer_base_st": outer_span_start,
                        "outer_base_nd": outer_span_end,
                        "span_st": span_start,
                        "span_nd": span_end,
                        "span_len": span_len,
                        "ratio": ratio,
                        "total_ref_len": total_ref_len,
                        "normal_haploid_depth": normal_haploid_depth,
                        "expected_high_sides": type4_indel_expected_high_sides(
                            indel_kind
                        ),
                        "depth_steps": tuple(depth_steps),
                        "depth_min_step": (
                            TYPE4_INDEL_GRAPH_DEPTH_DIFF_RATIO
                            * normal_haploid_depth
                        ),
                        "depth_filter_reason": depth_filter_reason,
                        "censat_jump": censat_jump,
                        "contig_name": contig_start[CTG_NAM],
                    }
                )

        start_idx = end_idx + 1
    return candidates


def select_type4_indel_graph_edges(
    contig_data: Sequence,
    nclose_nodes: Mapping,
    type4_indel_candidates: Sequence[Mapping],
) -> list:
    breakends_by_chrom = collect_nclose_breakends(contig_data, nclose_nodes)
    selected_edges = []

    for type4_idx, candidate in enumerate(type4_indel_candidates):
        start_idx, end_idx = candidate["type4_tuple"]
        chrom = candidate["chrom"]
        span_start = candidate["span_st"]
        span_end = candidate["span_nd"]
        inside_nclose_pairs = {
            pair_key
            for breakend, pair_key, _ in breakends_by_chrom.get(chrom, ())
            if span_start <= breakend <= span_end
        }
        if candidate["indel_kind"] == "insertion" and len(inside_nclose_pairs) < 2:
            continue

        for src, dst in (
            ((DIR_FOR, start_idx), (DIR_FOR, end_idx)),
            ((DIR_BAK, end_idx), (DIR_BAK, start_idx)),
        ):
            selected_edges.append(
                {
                    "type4_kind": candidate["type4_kind"],
                    "type4_idx": type4_idx,
                    "type4_tuple": tuple(candidate["type4_tuple"]),
                    "src": src,
                    "dst": dst,
                    "dist": get_graph_edge_weight(contig_data, src[1], dst[1]),
                    "indel_kind": candidate["indel_kind"],
                    "chrom": chrom,
                    "base_chrom": candidate["base_chrom"],
                    "base_st": candidate["base_st"],
                    "base_nd": candidate["base_nd"],
                    "outer_base_st": candidate["outer_base_st"],
                    "outer_base_nd": candidate["outer_base_nd"],
                    "span_st": span_start,
                    "span_nd": span_end,
                    "span_len": candidate["span_len"],
                    "inside_nclose_count": len(inside_nclose_pairs),
                    "ratio": candidate["ratio"],
                    "expected_high_sides": candidate["expected_high_sides"],
                    "depth_steps": candidate["depth_steps"],
                    "depth_min_step": candidate["depth_min_step"],
                    "censat_jump": candidate["censat_jump"],
                    "contig_name": candidate["contig_name"],
                }
            )
    return selected_edges


def augment_nclose_nodes_with_type4_indels(
    nclose_nodes: Mapping,
    selected_edges: Sequence[Mapping],
) -> defaultdict:
    augmented = defaultdict(list)
    existing_pairs = defaultdict(set)
    for key, pair_list in nclose_nodes.items():
        augmented[key].extend(tuple(pair) for pair in pair_list)
        existing_pairs[key].update(tuple(pair) for pair in pair_list)

    for edge in selected_edges:
        pair = tuple(edge["type4_tuple"])
        key = edge["contig_name"]
        if pair in existing_pairs[key]:
            continue
        augmented[key].append(pair)
        existing_pairs[key].add(pair)
    return augmented


def get_type4_indel_zero_dim_edge_set(selected_edges: Sequence[Mapping]) -> set:
    return {
        (tuple(edge["src"]), tuple(edge["dst"])) for edge in selected_edges
    }


def build_graph_state(
    stage_input: Stage10Input,
    *,
    add_indel_graph: bool = False,
    depth_df: pd.DataFrame | None = None,
    censat_intervals: Mapping | None = None,
    vcf_input: bool = False,
) -> GraphBuildState:
    """Build one graph state without mutating canonical stage-01 NCloses."""

    contig_data = stage_input.contig_data
    canonical_before = canonical_nclose_snapshot(stage_input.nclose_nodes)
    selected_edges = []
    graph_nclose_nodes = stage_input.nclose_nodes
    if add_indel_graph:
        if depth_df is None or censat_intervals is None:
            raise ValueError(
                "--add-indel-graph requires depth and CEN-SAT inputs"
            )
        candidates = collect_type4_indel_graph_candidates(
            contig_data,
            depth_df,
            censat_intervals,
        )
        selected_edges = select_type4_indel_graph_edges(
            contig_data,
            stage_input.nclose_nodes,
            candidates,
        )
        graph_nclose_nodes = augment_nclose_nodes_with_type4_indels(
            stage_input.nclose_nodes,
            selected_edges,
        )
        selected_kind_count = Counter(
            edge["indel_kind"] for edge in selected_edges
        )
        logging.info(
            "Added %d graph-only type4 indel edges (%s)",
            len(selected_edges),
            dict(selected_kind_count),
        )

    if canonical_nclose_snapshot(stage_input.nclose_nodes) != canonical_before:
        raise AssertionError(
            "Graph preparation mutated canonical nclose_nodes"
        )

    if vcf_input:
        assert_vcf_nclose_has_no_indel_like(
            contig_data,
            stage_input.nclose_nodes,
            "stage-10 canonical nclose",
        )
        if add_indel_graph:
            non_vcf_edges = [
                edge
                for edge in selected_edges
                if not str(edge.get("contig_name", "")).startswith(
                    VCF_TYPE4_GRAPH_NODE_PREFIX
                )
            ]
            if non_vcf_edges:
                raise AssertionError(
                    "VCF stage-10 selected graph-only edges outside VCF "
                    f"synthetic nodes: {non_vcf_edges[:2]}"
                )

    _, chr_rev_corr = chr_correlation_maker(contig_data)
    adjacency = initialize_bnd_graph(
        contig_data,
        graph_nclose_nodes,
        stage_input.telo_contig,
        len(contig_data),
        chr_rev_corr,
    )
    return GraphBuildState(
        graph_nclose_nodes=graph_nclose_nodes,
        selected_type4_indel_graph_edges=tuple(selected_edges),
        type4_indel_zero_dim_edge_set=frozenset(
            get_type4_indel_zero_dim_edge_set(selected_edges)
        ),
        bnd_graph_adjacency=adjacency,
    )


def build_dimension_graph(
    bnd_graph_adjacency: Mapping,
    contig_data: Sequence,
    zero_dim_edge_set: set | frozenset,
    chr_change_limit: int,
    dir_change_limit: int,
) -> nx.DiGraph:
    """Expand the base adjacency with chromosome/direction change counters."""

    graph = nx.DiGraph()
    for node in bnd_graph_adjacency:
        for chr_count in range(chr_change_limit + 1):
            for dir_count in range(dir_change_limit + 1):
                if isinstance(node, str):
                    graph.add_node((node, chr_count, dir_count))
                else:
                    graph.add_node(tuple(node) + (chr_count, dir_count))

    for node, edge_list in bnd_graph_adjacency.items():
        for edge in edge_list:
            if isinstance(edge, str):
                for chr_count in range(chr_change_limit + 1):
                    for dir_count in range(dir_change_limit + 1):
                        graph.add_edge(
                            tuple(node) + (chr_count, dir_count),
                            (edge, chr_count, dir_count),
                            weight=0,
                        )
                continue

            if isinstance(node, str):
                for chr_count in range(chr_change_limit + 1):
                    for dir_count in range(dir_change_limit + 1):
                        graph.add_edge(
                            (node, chr_count, dir_count),
                            tuple(edge) + (chr_count, dir_count),
                            weight=0,
                        )
                continue

            source_idx = int(node[1])
            target_idx = int(edge[1])
            source_contig = contig_data[source_idx]
            target_contig = contig_data[target_idx]
            weight = get_graph_edge_weight(
                contig_data,
                source_idx,
                target_idx,
            )
            zero_dimension = (
                (tuple(node), tuple(edge[:2])) in zero_dim_edge_set
            )

            if zero_dimension:
                for chr_count in range(chr_change_limit + 1):
                    for dir_count in range(dir_change_limit + 1):
                        graph.add_edge(
                            tuple(node) + (chr_count, dir_count),
                            tuple(edge) + (chr_count, dir_count),
                            weight=weight,
                        )
            elif source_contig[CHR_NAM] != target_contig[CHR_NAM]:
                for chr_count in range(chr_change_limit):
                    for dir_count in range(dir_change_limit + 1):
                        graph.add_edge(
                            tuple(node) + (chr_count, dir_count),
                            tuple(edge) + (chr_count + 1, dir_count),
                            weight=weight,
                        )
            elif (
                source_contig[CTG_DIR] != target_contig[CTG_DIR]
                and source_contig[CTG_NAM] == target_contig[CTG_NAM]
            ):
                for chr_count in range(chr_change_limit + 1):
                    for dir_count in range(dir_change_limit):
                        graph.add_edge(
                            tuple(node) + (chr_count, dir_count),
                            tuple(edge) + (chr_count, dir_count + 1),
                            weight=weight,
                        )
            else:
                for chr_count in range(chr_change_limit + 1):
                    for dir_count in range(dir_change_limit + 1):
                        graph.add_edge(
                            tuple(node) + (chr_count, dir_count),
                            tuple(edge) + (chr_count, dir_count),
                            weight=weight,
                        )
    return graph


def _to_graph_tool(nx_graph: nx.DiGraph):
    """Convert topology only, retaining Python node objects out of band."""

    import graph_tool.all as gt

    graph = gt.Graph(directed=True)
    node_to_vertex = {}
    vertex_to_node = {}
    for node in nx_graph.nodes:
        vertex = graph.add_vertex()
        node_to_vertex[node] = vertex
        vertex_to_node[int(vertex)] = node
    for source, target in nx_graph.edges:
        graph.add_edge(node_to_vertex[source], node_to_vertex[target])
    return gt, graph, node_to_vertex, vertex_to_node


def _write_debug_path(
    path: Sequence,
    chromosome_profile: Sequence,
    pair_count: int,
    context: _PathWorkerContext,
    source_name: str,
    target_name: str,
) -> None:
    if not context.verbose or context.raw_output_dir is None:
        return
    folder = Path(context.raw_output_dir) / f"{source_name}_{target_name}"
    folder.mkdir(parents=True, exist_ok=True)
    paf_path = folder / f"{pair_count}.paf"
    index_path = folder / f"{pair_count}.index.txt"
    with paf_path.open("wt", encoding="utf-8") as paf_handle, index_path.open(
        "wt", encoding="utf-8"
    ) as index_handle:
        for node in path:
            if isinstance(node[0], str):
                print(node, file=paf_handle)
                print(node, file=index_handle)
            else:
                paf_handle.write(
                    "\t".join(map(str, context.contig_data[node[1]])) + "\n"
                )
                index_handle.write("\t".join(map(str, node)) + "\n")
        print(chromosome_profile, file=paf_handle)
        print(chromosome_profile, file=index_handle)


def _path_is_acceptable(
    path: Sequence,
    source_name: str,
    target_name: str,
    path_compress: Mapping,
    context: _PathWorkerContext,
) -> tuple[bool, list, Counter]:
    """Apply all retained legacy path filters and profile compression."""

    path_len = len(path)
    if path_len < 4:
        return False, [], Counter()

    contig_data = context.contig_data
    path_counter = Counter()
    if (
        contig_data[path[1][1]][CTG_NAM]
        == contig_data[path[2][1]][CTG_NAM]
    ):
        path_counter[contig_data[path[1][1]][CHR_NAM]] += (
            int(contig_data[path[1][1]][CHR_END])
            - int(contig_data[path[1][1]][CHR_STR])
        )
    if (
        contig_data[path[-2][1]][CTG_NAM]
        == contig_data[path[-3][1]][CTG_NAM]
    ):
        path_counter[contig_data[path[-2][1]][CHR_NAM]] += (
            int(contig_data[path[-2][1]][CHR_END])
            - int(contig_data[path[-2][1]][CHR_STR])
        )

    contig_visits = Counter({path[-2][1]: 1})
    censat_visits = 0
    for node_offset in range(1, path_len - 2):
        node_idx = path[node_offset][1]
        next_idx = path[node_offset + 1][1]
        contig_visits[node_idx] += 1
        if contig_visits[node_idx] >= BND_OVERUSE_CNT:
            return False, [], Counter()

        current = contig_data[node_idx]
        following = contig_data[next_idx]
        current_len = int(current[CHR_END]) - int(current[CHR_STR])
        following_len = int(following[CHR_END]) - int(following[CHR_STR])
        if node_is_censat(current) and node_is_censat(following):
            censat_visits += 1
            if censat_visits >= CENSAT_VISIT_LIMIT:
                return False, [], Counter()

        if current[CTG_NAM] == following[CTG_NAM]:
            if node_idx < next_idx:
                intermediate = range(node_idx + 1, next_idx)
            else:
                intermediate = range(next_idx + 1, node_idx)
            for contig_idx in intermediate:
                node = contig_data[contig_idx]
                path_counter[node[CHR_NAM]] += (
                    int(node[CHR_END]) - int(node[CHR_STR])
                )
        elif distance_checker(current, following) == 0:
            path_counter[current[CHR_NAM]] += (
                current_len
                + following_len
                - overlap_calculator(current, following)
            )
        else:
            path_counter[current[CHR_NAM]] += (
                distance_checker(current, following)
                + current_len
                + following_len
            )

    total_path_ref_len = sum(path_counter.values())
    for chrom in [
        chrom
        for chrom, ref_len in path_counter.items()
        if ref_len < IGNORE_PATH_LIMIT
    ]:
        del path_counter[chrom]
    if not path_counter:
        return False, [], Counter()

    chromosome_profile = sorted(
        path_counter.items(),
        key=lambda value: -value[1],
    )
    longest = chromosome_profile[0]
    second = chromosome_profile[1] if len(chromosome_profile) > 1 else longest

    def rank_normalized(chrom):
        return "chrX" if chrom == "chrY" else chrom

    source_chrom = rank_normalized(source_name[:-1])
    target_chrom = rank_normalized(target_name[:-1])
    ranks = [
        rank
        for rank in range(min(PATH_MAJOR_COMPONENT, len(chromosome_profile)))
        if rank_normalized(chromosome_profile[rank][0])
        in (source_chrom, target_chrom)
    ]
    if len(ranks) == 2:
        has_major_terminals = True
    elif len(ranks) == 1:
        has_major_terminals = ranks[0] <= 1
    else:
        has_major_terminals = False

    if not has_major_terminals:
        return False, [], Counter()
    if (longest[1] + second[1]) / total_path_ref_len < 0.5:
        return False, [], Counter()
    if (
        context.chr_len[longest[0]] + context.chr_len[second[0]]
        < total_path_ref_len
    ):
        return False, [], Counter()
    if total_path_ref_len < MIN_PATH_REF_LEN:
        return False, [], Counter()

    profile_key = tuple(sorted(path_counter))
    for previous_counter in path_compress[profile_key]:
        if all(
            abs(previous_counter[chrom] - path_counter[chrom])
            <= PATH_COMPRESS_LIMIT
            for chrom in profile_key
        ):
            return False, [], Counter()
    return True, chromosome_profile, path_counter


def _run_terminal_pair(data):
    context = _WORKER_CONTEXT
    if context is None:
        raise RuntimeError("Stage-10 path worker was not initialised")

    source_idx, target_idx = data
    source_name = context.chr_rev_corr[source_idx]
    target_name = context.chr_rev_corr[target_idx]
    source = (source_name, 0, 0)
    graph_without_other_terminals = copy.deepcopy(context.graph)

    for terminal_idx in range(
        context.contig_data_size,
        context.contig_data_size + 2 * CHROMOSOME_COUNT,
    ):
        terminal_name = context.chr_rev_corr[terminal_idx]
        if terminal_name in (source_name, target_name):
            continue
        for chr_count in range(context.chr_change_limit + 1):
            for dir_count in range(context.dir_change_limit + 1):
                node = (terminal_name, chr_count, dir_count)
                if graph_without_other_terminals.has_node(node):
                    graph_without_other_terminals.remove_node(node)

    path_list = []
    path_compress = defaultdict(list)
    pair_count = 0
    for chr_count in range(context.chr_change_limit + 1):
        for dir_count in range(context.dir_change_limit + 1):
            candidate_graph = copy.deepcopy(graph_without_other_terminals)
            if source_name != target_name:
                for remove_chr in range(context.chr_change_limit + 1):
                    for remove_dir in range(context.dir_change_limit + 1):
                        if remove_chr or remove_dir:
                            candidate_graph.remove_nodes_from(
                                [(source_name, remove_chr, remove_dir)]
                            )
                        if remove_chr != chr_count or remove_dir != dir_count:
                            candidate_graph.remove_nodes_from(
                                [(target_name, remove_chr, remove_dir)]
                            )
            else:
                for remove_chr in range(context.chr_change_limit + 1):
                    for remove_dir in range(context.dir_change_limit + 1):
                        if (
                            (remove_chr != chr_count or remove_dir != dir_count)
                            and remove_chr + remove_dir != 0
                        ):
                            candidate_graph.remove_nodes_from(
                                [(target_name, remove_chr, remove_dir)]
                            )

            target = (target_name, chr_count, dir_count)
            if not candidate_graph.has_node(source) or not candidate_graph.has_node(
                target
            ):
                continue
            gt, tool_graph, node_to_vertex, vertex_to_node = _to_graph_tool(
                candidate_graph
            )
            source_vertex = node_to_vertex[source]
            target_vertex = node_to_vertex[target]
            if (
                gt.shortest_distance(
                    tool_graph,
                    source=source_vertex,
                    target=target_vertex,
                )
                == float("inf")
            ):
                continue

            for path_vertices in gt.all_paths(
                tool_graph,
                source_vertex,
                target_vertex,
            ):
                path = [vertex_to_node[int(vertex)] for vertex in path_vertices]
                accepted, chromosome_profile, profile_counter = _path_is_acceptable(
                    path,
                    source_name,
                    target_name,
                    path_compress,
                    context,
                )
                if not accepted:
                    continue

                profile_key = tuple(sorted(profile_counter))
                path_compress[profile_key].append(copy.deepcopy(profile_counter))
                pair_count += 1
                _write_debug_path(
                    path,
                    chromosome_profile,
                    pair_count,
                    context,
                    source_name,
                    target_name,
                )
                path_list.append((path, chromosome_profile))
                if pair_count >= context.per_pair_path_limit:
                    return (source_name, target_name), pair_count, path_list

    return (source_name, target_name), pair_count, path_list


def _init_path_worker(context: _PathWorkerContext) -> None:
    global _WORKER_CONTEXT
    _WORKER_CONTEXT = context


def _limit_combinations(
    config: PathSearchConfig,
    nclose_node_count: int,
) -> tuple[list[tuple[int, int]], int]:
    if config.fixed_limit_combination is not None:
        return [config.fixed_limit_combination], 0

    combinations = [
        (chr_limit, DIR_CHANGE_LIMIT_ABS_MAX)
        for chr_limit in range(config.graph_depth, 0, -1)
    ]
    combinations.append((1, 0))
    if nclose_node_count <= HARD_GRAPH_NCLOSE_COUNT:
        return combinations, 0

    hard_start = min(
        CHR_CHANGE_LIMIT_HARD_START,
        max(config.graph_depth, 1),
    )
    return combinations, combinations.index((hard_start, DIR_CHANGE_LIMIT_ABS_MAX))


def run_path_search(
    stage_input: Stage10Input,
    graph_state: GraphBuildState,
    chr_len: Mapping[str, int],
    config: PathSearchConfig,
) -> PathSearchResult | None:
    """Enumerate candidate paths with deterministic result serialization."""

    contig_data_size = len(stage_input.contig_data)
    _, chr_rev_corr = chr_correlation_maker(stage_input.contig_data)
    terminal_pairs = [
        (source, target)
        for source in range(
            contig_data_size,
            contig_data_size + 2 * CHROMOSOME_COUNT,
        )
        for target in range(
            source + 1,
            contig_data_size + 2 * CHROMOSOME_COUNT,
        )
    ]
    nclose_node_count = count_nclose_nodes(stage_input.nclose_nodes)
    combinations, combination_idx = _limit_combinations(
        config,
        nclose_node_count,
    )
    if config.fixed_limit_combination is not None:
        logging.info(
            "Using fixed limit_combinations: %s",
            config.fixed_limit_combination,
        )

    last_success = None
    while 0 <= combination_idx < len(combinations):
        chr_limit, dir_limit = combinations[combination_idx]
        dimension_graph = build_dimension_graph(
            graph_state.bnd_graph_adjacency,
            stage_input.contig_data,
            graph_state.type4_indel_zero_dim_edge_set,
            chr_limit,
            dir_limit,
        )
        worker_context = _PathWorkerContext(
            graph=dimension_graph,
            contig_data=stage_input.contig_data,
            contig_data_size=contig_data_size,
            chr_rev_corr=chr_rev_corr,
            chr_len=chr_len,
            chr_change_limit=chr_limit,
            dir_change_limit=dir_limit,
            per_pair_path_limit=config.per_pair_path_limit,
            verbose=config.verbose,
            raw_output_dir=config.raw_output_dir,
        )

        results = []
        total_path_count = 0
        with Pool(
            processes=max(int(config.thread), 1),
            initializer=_init_path_worker,
            initargs=(worker_context,),
        ) as pool:
            iterator = pool.imap_unordered(_run_terminal_pair, terminal_pairs)
            for result in tqdm(
                iterator,
                total=len(terminal_pairs),
                desc=f"Build breakend graph ({chr_limit},{dir_limit})",
                disable=not sys.stdout.isatty() and not config.progress,
            ):
                results.append(result)
                total_path_count += result[1]
                if total_path_count >= config.total_path_limit:
                    pool.terminate()
                    break

        success = total_path_count < config.total_path_limit
        if success:
            results.sort(key=lambda result: result[0])
            path_items = [
                (f"{terminal_pair[0]}_{terminal_pair[1]}", paths)
                for terminal_pair, _, paths in results
            ]
            counts = [
                (terminal_pair, count)
                for terminal_pair, count, _ in results
            ]
            last_success = PathSearchResult(
                chr_change_limit=chr_limit,
                dir_change_limit=dir_limit,
                path_items=path_items,
                counts=counts,
            )
            logging.info(
                "SUCCESS at %s, with %d paths",
                (chr_limit, dir_limit),
                total_path_count,
            )
            if total_path_count >= config.total_path_limit / 2:
                logging.info("Estimated next graph limit would exceed path budget")
                break
            combination_idx -= 1
            continue

        logging.info("FAIL at %s", (chr_limit, dir_limit))
        if config.fixed_limit_combination is not None or last_success is not None:
            break
        combination_idx += 1

    return last_success
