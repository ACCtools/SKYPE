"""Native SKYPE stage orchestration for assembly and VCF inputs.

ACCtools prepares references and assembly alignments, then invokes this CLI.
Stage commands, completion checks, and partial restarts are owned here.
"""

import argparse
import logging
import os
import shlex
import subprocess
import shutil
import json
from pathlib import Path

import psutil

from skype_options import normalize_extra_args

MEM_SAFE_RATIO = 0.8


class SkypeArgumentError(ValueError):
    pass


def subprocess_print(args, **kwargs):
    print(*args)


def require_restart_inputs(prefix, stage, required_paths):
    missing = [
        os.path.join(prefix, path)
        for path in required_paths
        if not os.path.exists(os.path.join(prefix, path))
    ]
    if missing:
        raise SkypeArgumentError(
            f"Cannot restart at stage {stage}; missing prerequisite artifact(s): "
            + ", ".join(missing)
        )


def run_skype(CELL_LINE, PREFIX, ctg_paf, ctg_aln_paf, utg_paf, utg_aln_paf,
              depth_loc, thread, dep_folder, is_progress, skype_force, graph_depth,
              option_skype="", skype_start_at=0, print_args=False,
              benchmark_vcf_loc=None, vcf_ins_aln_paf=None,
              unitig_fasta=None, alignment_force=False, *,
              alignasm_ref, chr_fai, tel_bed, rpt_bed, rcs_bed, cyt_bed,
              ref_stat=None, raw_rescue_method=None, raw_rescue_options=None,
              reference_index_cache=None):
    # Execute the core SKYPE analysis scripts.
    dep_folder = os.path.abspath(dep_folder)
    reference_index_cache = os.path.abspath(
        reference_index_cache or os.path.join(dep_folder, 'reference_indexes'))
    skype_folder_loc = os.path.dirname(os.path.abspath(__file__))

    rescue_parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    rescue_parser.add_argument('--raw-rescue-method', '--raw_rescue_method',
                               choices=('off', 'read', 'olc'), default=None)
    rescue_parser.add_argument('--raw-rescue-options', default='')
    rescue_parser.add_argument('--depth-policy', choices=('legacy', 'all', 'nclose_l2', 'nclose_huber'), default='legacy')
    rescue_parser.add_argument('--robust-sigma-multiplier', type=float, default=3.)
    rescue_options, remaining_options = rescue_parser.parse_known_args(
        normalize_extra_args(option_skype)
    )
    raw_rescue_method = raw_rescue_method or rescue_options.raw_rescue_method
    previous_rescue = Path(PREFIX) / '24_raw_rescue' / 'summary.json'
    if raw_rescue_method is None and skype_start_at > 1 and previous_rescue.exists():
        raw_rescue_method = json.loads(previous_rescue.read_text())['method']
    raw_rescue_method = raw_rescue_method or ('off' if benchmark_vcf_loc else 'read')
    rescue_args = shlex.split(rescue_options.raw_rescue_options if raw_rescue_options is None else raw_rescue_options)
    option_skype = shlex.join(remaining_options)

    valid_start_stages = {0, 1, 10, 11, 21, 22, 23, 24, 31}
    if skype_start_at not in valid_start_stages:
        valid_text = ", ".join(map(str, sorted(valid_start_stages)))
        raise SkypeArgumentError(
            f"Native --skype_start_at must be one of {valid_text} "
            "for the native 01/10 route"
        )

    TEL_BED = tel_bed
    CHR_FAI = chr_fai
    RPT_BED = rpt_bed
    RCS_BED = rcs_bed
    CYT_BED = cyt_bed
    REF_STAT_LOC = ref_stat


    # python to bash variable
    MAIN_STAT_LOC = depth_loc
    MAIN_STAT_NORM_LOC = depth_loc if REF_STAT_LOC is None else depth_loc.replace('.win.stat.gz', '_normalized.win.stat.gz')
    PAF_LOC = ctg_aln_paf
    PAF_UTG_LOC = utg_aln_paf
    graph_paf_loc = PAF_UTG_LOC if benchmark_vcf_loc else PAF_LOC
    PPC_PAF_LOC = os.path.join(PREFIX, f"{os.path.basename(graph_paf_loc)}.ppc.paf")
    READ_BAM_LOC = os.path.join(os.path.dirname(depth_loc), f'{CELL_LINE}.bam')
    SORTED_READ_BAM_LOC = os.path.join(os.path.dirname(depth_loc), f'{CELL_LINE}.sorted.bam')
    if os.path.isfile(SORTED_READ_BAM_LOC):
        READ_BAM_LOC = SORTED_READ_BAM_LOC

    THREAD = str(thread)
    PROGRESS = ["--progress"] if is_progress else []

    subprocess_run = subprocess_print if print_args else subprocess.run
    EXTRA_SKYPE = ["--option_skype=" + shlex.join(normalize_extra_args(option_skype))]

    expected_outputs = [
        os.path.join(PREFIX, "total_cov.png"),
        os.path.join(PREFIX, "nclose_report.tsv"),
    ]
    if benchmark_vcf_loc:
        expected_outputs.append(os.path.join(PREFIX, "SV_benchmark_result.vcf"))
    else:
        expected_outputs.extend([
            os.path.join(PREFIX, "SV_call_result.vcf"),
            os.path.join(PREFIX, "SKYPE_result.bed"),
        ])

    if not print_args and skype_start_at > 0:
        if REF_STAT_LOC is not None and not os.path.isfile(MAIN_STAT_NORM_LOC):
            raise SkypeArgumentError(
                "Cannot restart SKYPE because the normalized depth file is "
                f"missing: {MAIN_STAT_NORM_LOC}"
            )
        restart_requirements = {
            10: [
                "01_nclose_data.pkl",
                "pipeline_mode.pkl",
                os.path.basename(PPC_PAF_LOC),
            ],
            11: [
                "path_data.pkl",
                "paf_file_path.pkl",
                os.path.basename(PPC_PAF_LOC),
            ],
            21: [
                "path_data.pkl",
                "11_ref_ratio_outliers",
                os.path.basename(PPC_PAF_LOC),
            ],
            22: ["path_data.pkl", "contig_pat_vec_data.pkl"],
            23: ["23_input.pkl"],
            24: ["01_nclose_data.pkl", "23_input.pkl", "B.npy", "predict_B.npy", "paf_file_path.pkl"],
            31: [
                "B.npy",
                "weight.npy",
                "predict_B.npy",
                "contig_pat_vec_data.pkl",
                "tot_loc_list.pkl",
                "23_input.pkl",
            ],
        }
        if skype_start_at in restart_requirements:
            require_restart_inputs(
                PREFIX,
                skype_start_at,
                restart_requirements[skype_start_at],
            )

    if not all(os.path.isfile(path) for path in expected_outputs) or skype_force or skype_start_at > 0:
        if skype_start_at <= 0 and REF_STAT_LOC is not None:
            subprocess_run([
                "python", os.path.join(skype_folder_loc, "00_depth_norm.py"),
                MAIN_STAT_LOC, REF_STAT_LOC, RCS_BED
            ], check=True)

        if skype_start_at <= 1:
            preprocess_cmd = [
                "python",
                os.path.join(skype_folder_loc, "01_Preprocess_NClose.py"),
                graph_paf_loc,
                CHR_FAI,
                TEL_BED,
                RPT_BED,
                RCS_BED,
                MAIN_STAT_NORM_LOC,
                PREFIX,
                READ_BAM_LOC,
                "-t",
                THREAD,
            ]
            if benchmark_vcf_loc:
                preprocess_cmd.extend([
                    "--vcf_input",
                    os.path.abspath(benchmark_vcf_loc),
                ])
                if vcf_ins_aln_paf is not None:
                    preprocess_cmd.extend([
                        "--alt",
                        os.path.abspath(vcf_ins_aln_paf),
                    ])
            else:
                if unitig_fasta is None:
                    raise SkypeArgumentError("Assembly preprocessing requires the effective unitig FASTA")
                censat_dir = os.path.splitext(utg_paf)[0] + ".censat_endpoints"
                prepare_cmd = [
                    "python", os.path.join(skype_folder_loc, "censat_endpoints.py"),
                    "--aln-paf", PAF_UTG_LOC,
                    "--raw-paf", utg_paf,
                    "--fasta", unitig_fasta,
                    "--reference", alignasm_ref,
                    "--reference-index-cache", reference_index_cache,
                    "--bed", RCS_BED,
                    "--outdir", censat_dir,
                    "-t", THREAD,
                ]
                if alignment_force:
                    prepare_cmd.append("--force")
                subprocess_run(prepare_cmd, check=True)
                preprocess_cmd.extend([
                    "--censat-endpoints-dir", censat_dir,
                    "--alt",
                    PAF_UTG_LOC,
                    "--original_paf_loc",
                    ctg_paf,
                    utg_paf,
                ])
            subprocess_run(preprocess_cmd + EXTRA_SKYPE + PROGRESS, check=True)
            if not print_args:
                rescue_dir = Path(PREFIX) / '24_raw_rescue'
                if rescue_dir.exists():
                    shutil.rmtree(rescue_dir)

        def run_fit(start_at):
            if start_at <= 10:
                graph_cmd = [
                    "python",
                    os.path.join(skype_folder_loc, "10_Graph_Find_Paths.py"),
                    os.path.join(PREFIX, "01_nclose_data.pkl"),
                    CHR_FAI,
                    PREFIX,
                    "-t",
                    THREAD,
                    "-d",
                    str(graph_depth),
                    "--main-stat-path",
                    MAIN_STAT_NORM_LOC,
                    "--censat-bed-path",
                    RCS_BED,
                ]
                if skype_start_at == 10 and option_skype:
                    graph_cmd.extend(EXTRA_SKYPE)
                subprocess_run(graph_cmd + PROGRESS, check=True)

            if start_at <= 11:
                subprocess_run([
                    "python", os.path.join(skype_folder_loc, "11_Ref_Outlier_Contig_Modify.py"),
                    CHR_FAI, PPC_PAF_LOC, PREFIX,
                ], check=True)

            if start_at <= 21:
                free_mem_gb = psutil.virtual_memory().available * MEM_SAFE_RATIO / (1024 ** 3)
                thread_lim = int(free_mem_gb / 6)

                subprocess_run([
                    "python", os.path.join(skype_folder_loc, "21_run_depth.py"),
                    PPC_PAF_LOC, PREFIX,
                    "--pandepth_loc", os.path.join(dep_folder, 'PanDepth', 'bin', 'pandepth'),
                    "-t", str(max(min(thread_lim, thread), 1))
                ] + PROGRESS, check=True)

            if start_at <= 22:
                subprocess_run([
                    "python", os.path.join(skype_folder_loc, "22_save_matrix.py"),
                    RCS_BED, MAIN_STAT_NORM_LOC,
                    PREFIX, "-t", THREAD,
                    '--depth-policy', rescue_options.depth_policy,
                    '--robust-sigma-multiplier', str(rescue_options.robust_sigma_multiplier),
                ] + PROGRESS, check=True)

            if start_at <= 23:
                # NNLS is faster with one BLAS thread on these depth matrices.
                subprocess_run([
                    "python", "23_run_nnls.py", os.path.abspath(PREFIX), "-t", "1"
                ], check=True, cwd=skype_folder_loc)

        run_fit(skype_start_at)

        if skype_start_at <= 24 and raw_rescue_method != 'off':
            rescue_dir = Path(PREFIX) / '24_raw_rescue'
            rescue_cmd = [
                'python', os.path.join(skype_folder_loc, '24_raw_nclose_rescue.py'),
                os.path.abspath(PREFIX), '--method', raw_rescue_method,
                '--bam', READ_BAM_LOC, '--reference', alignasm_ref,
                '--reference-index-cache', reference_index_cache,
                '--censat-bed', RCS_BED, '--repeat-bed', RPT_BED,
                '--ppc-paf', PPC_PAF_LOC, '-t', THREAD,
            ]
            subprocess_run(rescue_cmd + rescue_args, check=True)
            if not print_args:
                summary = json.loads((rescue_dir / 'summary.json').read_text())
                round_path = rescue_dir / 'round.json'
                round_state = json.loads(round_path.read_text()) if round_path.exists() else {}
                if summary['added_count'] and not round_state.get('refit_complete'):
                    augmented = Path(summary['augmented_dir'])
                    before = rescue_dir / 'before'
                    # Installing the same prepared snapshot again is safe after
                    # an interrupted copy. Discovery never runs a second time.
                    for filename in summary['install_files']:
                        target = Path(PREFIX) / filename
                        if target.exists() and not (before / filename).exists():
                            shutil.copyfile(target, before / filename)
                        shutil.copyfile(augmented / filename, target)
                    round_path.write_text(json.dumps(dict(applied=True, refit_complete=False, rounds=1), indent=2)+'\n')
                    logging.info('Stage 24 added %d NClose(s); rerunning stages 10--23 once', summary['added_count'])
                    run_fit(10)
                    round_path.write_text(json.dumps(dict(applied=True, refit_complete=True, rounds=1), indent=2)+'\n')

        if skype_start_at <= 31:
            subprocess_run([
                "python", os.path.join(skype_folder_loc, "31_depth_analysis.py"),
                RCS_BED, PPC_PAF_LOC, MAIN_STAT_NORM_LOC,
                TEL_BED, CHR_FAI, CYT_BED, PREFIX, "-t", THREAD
            ] + PROGRESS, check=True)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    for option, dest in (
        ("cell-line", "CELL_LINE"),
        ("prefix", "PREFIX"),
        ("ctg-paf", "ctg_paf"),
        ("ctg-aln-paf", "ctg_aln_paf"),
        ("utg-paf", "utg_paf"),
        ("utg-aln-paf", "utg_aln_paf"),
        ("depth-loc", "depth_loc"),
        ("dep-folder", "dep_folder"),
        ("alignasm-ref", "alignasm_ref"),
        ("chr-fai", "chr_fai"),
        ("tel-bed", "tel_bed"),
        ("rpt-bed", "rpt_bed"),
        ("rcs-bed", "rcs_bed"),
        ("cyt-bed", "cyt_bed"),
    ):
        parser.add_argument("--" + option, dest=dest, required=True)
    for option in ("ref-stat", "benchmark-vcf-loc", "vcf-ins-aln-paf", "unitig-fasta"):
        parser.add_argument("--" + option)
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument("-d", "--graph-depth", type=int, default=4)
    parser.add_argument("--progress", dest="is_progress", action="store_true")
    parser.add_argument("--skype-force", "--skype_force", action="store_true")
    parser.add_argument("--alignment-force", action="store_true")
    parser.add_argument("--option-skype", "--option_skype", default="")
    parser.add_argument("--skype-start-at", "--skype_start_at", type=int, default=0)
    parser.add_argument("--print-args", "--print_args", action="store_true")
    parser.add_argument('--raw-rescue-method', choices=('off', 'read', 'olc'))
    parser.add_argument('--raw-rescue-options')
    parser.add_argument('--reference-index-cache')
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    run_skype(**vars(args))


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s:%(message)s",
        datefmt="%m/%d/%Y %I:%M:%S %p",
    )
    try:
        main()
    except (FileNotFoundError, SkypeArgumentError) as exc:
        logging.error("%s", exc)
        raise SystemExit(1) from None
