#!/usr/bin/env python
"""Preprocess assembly/VCF input and write the stage-10 NClose handoff."""

from __future__ import annotations

import argparse
import logging

from nclose_preprocess import (
    NCLOSE_COUNT_DEFAULT_VAF_THRESHOLD,
    Stage01Config,
    parse_debug_nclose_endpoint,
    run_stage01,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Preprocess contig alignments, discover/filter NCloses, and write "
            "01_nclose_data.pkl"
        )
    )
    parser.add_argument(
        "paf_file_path",
        help=(
            "Compatibility primary PAF name (not read in assembly mode). "
            "In VCF mode this is the unitig PAF used for telomere anchors."
        ),
    )
    parser.add_argument("reference_fai_path")
    parser.add_argument("telomere_bed_path")
    parser.add_argument("repeat_bed_path")
    parser.add_argument("censat_bed_path")
    parser.add_argument("main_stat_path")
    parser.add_argument("prefix", help="Stage output directory")
    parser.add_argument("read_bam_path", help="Raw-read alignment BAM")
    parser.add_argument(
        "--alt",
        help=(
            "Required aligned unitig PAF in assembly mode; optional VCF INS "
            "alignment PAF in VCF mode"
        ),
    )
    parser.add_argument(
        "--original-paf-loc",
        "--original_paf_loc",
        dest="original_paf_loc",
        nargs="+",
        default=(),
        help=(
            "Legacy source PAF list; assembly stage 01 uses only its last "
            "value as the original unitig PAF"
        ),
    )
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument("--progress", action="store_true")
    parser.add_argument(
        "--exclude-nclose-list-loc",
        "--exclude_nclose_list_loc",
        dest="exclude_nclose_list_loc",
        default=None,
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--skip-bam-analysis",
        "--skip_bam_analysis",
        dest="skip_bam_analysis",
        action="store_true",
    )
    parser.add_argument(
        "--check-nclose-count",
        "--check_nclose_count",
        dest="check_nclose_count",
        action="store_true",
        help="Filter eligible NCloses using raw-read junction VAF",
    )
    parser.add_argument(
        "--nclose-count-vaf-threshold",
        "--nclose_count_vaf_threshold",
        dest="nclose_count_vaf_threshold",
        type=float,
        default=NCLOSE_COUNT_DEFAULT_VAF_THRESHOLD,
    )
    parser.add_argument(
        "--disable-alt-ctg-simple",
        "--disable_alt_ctg_simple",
        dest="disable_alt_ctg_simple",
        action="store_true",
        help="Disable primary-contig simple-alt rescue when that stage is enabled",
    )
    parser.add_argument(
        "--vcf-input",
        "--vcf_input",
        dest="vcf_input",
        default=None,
        help="Use retained VCF events instead of assembly-derived NCloses",
    )
    parser.add_argument(
        "--vcf-filter-pass",
        "--vcf_filter_pass",
        dest="vcf_filter_pass",
        nargs="+",
        metavar="FILTER",
        default=("PASS", "."),
    )
    parser.add_argument(
        "--debug-force-nclose",
        "--debug_force_nclose",
        dest="debug_force_nclose",
        action="append",
        nargs=2,
        type=parse_debug_nclose_endpoint,
        metavar=("ENDPOINT_A", "ENDPOINT_B"),
        default=None,
        help=(
            "Assembly-only debug injection of one synthetic NClose using two "
            "1-based chrom:pos:+/- endpoints; repeat the option for more pairs"
        ),
    )
    parser.add_argument("--test", action="store_true", help=argparse.SUPPRESS)
    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.thread < 1:
        parser.error("--thread must be positive")
    if args.vcf_input is None and args.alt is None:
        parser.error(
            "assembly mode requires --alt with the aligned unitig PAF; "
            "the positional primary PAF is not read"
        )
    if args.vcf_input is None and not args.original_paf_loc:
        parser.error(
            "assembly mode requires --original-paf-loc with the original "
            "unitig PAF as its last value"
        )
    if args.vcf_input is not None and args.debug_force_nclose:
        parser.error(
            "--debug-force-nclose is supported only in assembly mode"
        )

    logging.info("01_Preprocess_NClose start")
    config = Stage01Config(
        paf_file_path=args.paf_file_path,
        reference_fai_path=args.reference_fai_path,
        telomere_bed_path=args.telomere_bed_path,
        repeat_bed_path=args.repeat_bed_path,
        censat_bed_path=args.censat_bed_path,
        main_stat_path=args.main_stat_path,
        prefix=args.prefix,
        read_bam_path=args.read_bam_path,
        alt_path=args.alt,
        original_paf_paths=tuple(args.original_paf_loc),
        thread=args.thread,
        progress=args.progress,
        exclude_nclose_list_path=args.exclude_nclose_list_loc,
        skip_bam_analysis=args.skip_bam_analysis,
        check_nclose_count=args.check_nclose_count,
        nclose_count_vaf_threshold=args.nclose_count_vaf_threshold,
        disable_alt_ctg_simple=args.disable_alt_ctg_simple,
        vcf_input_path=args.vcf_input,
        vcf_filter_pass=tuple(args.vcf_filter_pass),
        debug_force_ncloses=tuple(
            tuple(endpoint_pair)
            for endpoint_pair in (args.debug_force_nclose or ())
        ),
    )
    result = run_stage01(config)
    logging.info(
        "01_Preprocess_NClose complete: %d contig node(s), %d NClose pair(s)",
        len(result.contig_data),
        sum(len(pairs) for pairs in result.nclose_nodes.values()),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
