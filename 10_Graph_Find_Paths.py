#!/usr/bin/env python
"""Build the native SKYPE breakend graph and enumerate candidate paths."""

from __future__ import annotations

import argparse
import logging
import os
import pickle
import shutil
import sys
from pathlib import Path

from skype_options import load_stage_options, split_stage_options

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from breakend_graph import (  # noqa: E402
    FAIL_NCLOSE_COUNT,
    LIMIT_COMBINATIONS_JSON,
    PathSearchConfig,
    build_graph_state,
    count_nclose_nodes,
    load_censat_intervals,
    load_chr_lengths,
    load_depth_dataframe,
    load_limit_combinations,
    load_stage10_input,
    run_path_search,
    write_limit_combinations,
)
from skype_utils import (  # noqa: E402
    TYPE4_INDEL_GRAPH_EDGE_PKL,
    load_pipeline_input,
    pipeline_input_is_vcf,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build a breakend graph from the stage-01 NClose handoff"
    )
    parser.add_argument(
        "nclose_data_path",
        help="Named stage-01 pickle containing contig_data, nclose_nodes, and telo_contig",
    )
    parser.add_argument("reference_fai_path")
    parser.add_argument("prefix")
    parser.add_argument("-t", "--thread", type=int, required=True)
    parser.add_argument("-d", "--graph-depth", type=int, default=4)
    parser.add_argument("--progress", action="store_true")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Write debug PAF/index paths below PREFIX/00_raw",
    )
    parser.add_argument(
        "--add-indel-graph",
        "--add_indel_graph",
        dest="add_indel_graph",
        action="store_true",
        help="Add selected type4 indel edges without consuming graph dimensions",
    )
    parser.add_argument("--main-stat-path")
    parser.add_argument("--censat-bed-path")
    parser.add_argument(
        "--limit-combinations",
        help="Use exactly the integer pair in a limit_combinations.json file",
    )
    parser.add_argument(
        "--option_skype", "--option-skype", default=None,
        help="Replace saved stage-10 options when restarting directly at stage 10",
    )
    return parser


def _remove_stale_stage10_outputs(prefix: Path, verbose: bool) -> Path:
    for filename in (
        "path_data.pkl",
        "path_di_data.pkl",
        TYPE4_INDEL_GRAPH_EDGE_PKL,
        LIMIT_COMBINATIONS_JSON,
        "report.txt",
    ):
        path = prefix / filename
        if path.is_file():
            path.unlink()

    raw_output_dir = prefix / "00_raw"
    if verbose and raw_output_dir.is_dir():
        shutil.rmtree(raw_output_dir)
    return raw_output_dir


def main(argv=None) -> int:
    parser = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(argv)
    try:
        if args.option_skype is None:
            stage_options = load_stage_options(args.prefix, "10")
        else:
            stage01_options, stage_options = split_stage_options(args.option_skype)
            if stage01_options:
                parser.error("Preprocessing options require restarting at stage 01")
    except ValueError as exc:
        parser.error(str(exc))
    args = parser.parse_args(stage_options + argv)
    if args.thread < 1:
        parser.error("--thread must be positive")
    if args.graph_depth < 0:
        parser.error("--graph-depth cannot be negative")
    if args.add_indel_graph and (
        args.main_stat_path is None or args.censat_bed_path is None
    ):
        parser.error(
            "--add-indel-graph requires --main-stat-path and --censat-bed-path"
        )

    fixed_limit = None
    if args.limit_combinations is not None:
        try:
            fixed_limit = load_limit_combinations(args.limit_combinations)
        except ValueError as exc:
            parser.error(f"--limit-combinations: {exc}")
        if fixed_limit[0] > max(args.graph_depth, 1):
            parser.error(
                "--limit-combinations chromosome-change limit "
                f"{fixed_limit[0]} exceeds graph depth {args.graph_depth}"
            )

    prefix = Path(args.prefix)
    prefix.mkdir(parents=True, exist_ok=True)
    raw_output_dir = _remove_stale_stage10_outputs(prefix, args.verbose)

    logging.info("10_Graph_Find_Paths start")
    stage_input = load_stage10_input(args.nclose_data_path)
    nclose_node_count = count_nclose_nodes(stage_input.nclose_nodes)
    # Establish the inactive contract before any graph-size early exit.  An
    # old successful run must never leak type4 edges into stage 21.
    with (prefix / TYPE4_INDEL_GRAPH_EDGE_PKL).open("wb") as handle:
        pickle.dump([], handle)
    logging.info(
        "Stage-10 input: %d contig nodes, %d NClose nodes, %d telomere groups",
        len(stage_input.contig_data),
        nclose_node_count,
        len(stage_input.telo_contig),
    )
    if nclose_node_count > FAIL_NCLOSE_COUNT:
        logging.error(
            "NClose node count is too high (%d > %d); primary-only retry is disabled",
            nclose_node_count,
            FAIL_NCLOSE_COUNT,
        )
        return 1

    chr_len = load_chr_lengths(args.reference_fai_path)
    depth_df = None
    censat_intervals = None
    if args.add_indel_graph:
        depth_df = load_depth_dataframe(args.main_stat_path)
        censat_intervals = load_censat_intervals(args.censat_bed_path)

    pipeline_input = load_pipeline_input(prefix)
    graph_state = build_graph_state(
        stage_input,
        add_indel_graph=args.add_indel_graph,
        depth_df=depth_df,
        censat_intervals=censat_intervals,
        vcf_input=pipeline_input_is_vcf(pipeline_input),
    )
    with (prefix / TYPE4_INDEL_GRAPH_EDGE_PKL).open("wb") as handle:
        pickle.dump(list(graph_state.selected_type4_indel_graph_edges), handle)

    result = run_path_search(
        stage_input,
        graph_state,
        chr_len,
        PathSearchConfig(
            thread=args.thread,
            graph_depth=args.graph_depth,
            fixed_limit_combination=fixed_limit,
            progress=args.progress,
            verbose=args.verbose,
            raw_output_dir=str(raw_output_dir),
        ),
    )
    if result is None:
        logging.error(
            "Breakend graph is too divergent; primary-only retry is disabled"
        )
        return 1

    path_data = result.as_path_dict()
    with (prefix / "path_data.pkl").open("wb") as handle:
        pickle.dump(path_data, handle)
    write_limit_combinations(
        prefix / LIMIT_COMBINATIONS_JSON,
        (result.chr_change_limit, result.dir_change_limit),
    )

    total_path_count = sum(count for _, count in result.counts)
    with (prefix / "report.txt").open("wt", encoding="utf-8") as handle:
        print(Path(args.nclose_data_path).stem, file=handle)
        print(total_path_count, file=handle)
        for terminal_pair, count in sorted(
            result.counts,
            key=lambda value: (-value[1], value[0]),
        ):
            if count:
                print(terminal_pair[0], terminal_pair[1], count, file=handle)

    logging.info(
        "Final success settings: %s",
        (result.chr_change_limit, result.dir_change_limit),
    )
    logging.info("Total path count: %d", total_path_count)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
