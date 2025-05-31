import os
import sys
import logging
import argparse
import pickle as pkl
import pandas as pd
import numpy as np
import bisect
from collections import defaultdict 

from juliacall import Main as jl

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("03_anal_bam start")

K = 1000
M = 1000000

AC_WA_RATIO_LIMIT = 5
FLANK_LENGTH = 1 * M

BREAKEND_DEPTH_RATIO = 0.2
ABS_MAX_COVERAGE_RATIO = 3
CENSAT_COMPRESSABLE_THRESHOLD = 1000 * K

def import_censat_repeat_data(file_path : str) -> dict :
    fai_file = open(file_path, "r")
    repeat_data = defaultdict(list)
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        ref_data = (int(temp_list[1]), int(temp_list[2]))
        if abs(ref_data[1] - ref_data[0]) > CENSAT_COMPRESSABLE_THRESHOLD:
            repeat_data[temp_list[0]].append(ref_data)
    fai_file.close()
    return repeat_data


def preprocess_breakends(nclose_cord_list, df, repeat_censat_data) -> list:
    # Todo : kill outliers respect to meandepth
    # Remove breakends from centromeres
    meandepth = np.median(df['meandepth'])
    using_nclose_list = []
    df_by_chr = {}
    for chrom, subdf in df.groupby("chr"):
        df_by_chr[chrom] = subdf.reset_index(drop=True)

    ends_map = {
        chrom: [iv[1] for iv in intervals]
        for chrom, intervals in repeat_censat_data.items()
    }

    for nclose_cord in nclose_cord_list:
        using_nclose = True
        key = nclose_cord[-1]

        for i in (0, 3):
            chrom = nclose_cord[i]
            coord = nclose_cord[i + 1]

            # If chromosome not present in df, skip to next coordinate
            if chrom not in df_by_chr:
                continue
            df_chr = df_by_chr[chrom]

            intervals = repeat_censat_data.get(chrom, [])
            ends = ends_map.get(chrom, [])
            idx = bisect.bisect_left(ends, coord)
            iv_start = 0
            iv_end = 0

            if idx < len(intervals):
                iv_start, iv_end = intervals[idx]
            if iv_start <= coord and coord <= iv_end:
                pre_fail_key_dict[key] = 'CENSAT_NCLOSE'
                using_nclose = False
                break

            df_bin = df_chr[(df_chr["st"] <= coord) & (coord < df_chr["nd"])]
            assert(df_bin.empty == False)
            if np.median(df_bin["meandepth"]) > ABS_MAX_COVERAGE_RATIO * meandepth:
                pre_fail_key_dict[key] = 'HIGH_DEPTH_REGION_NCLOSE'
                using_nclose = False
                break
            
        if using_nclose:
            using_nclose_list.append(nclose_cord)
             

    return using_nclose_list


def postprocess_breakends(df, pre_cord_nclose_list, ac_nclose_cnt_dict) -> dict:
    meandepth = np.median(df['meandepth'])
    # Pre‐group by chromosome for faster lookup
    df_by_chr = {chrom: subdf.reset_index(drop=True) for chrom, subdf in df.groupby("chr")}
    chr_max_end = {chrom: subdf["nd"].max() for chrom, subdf in df_by_chr.items()}

    
    using_ac_nclose_cnt_dict = {}
    for bd in pre_cord_nclose_list:
        key = bd[-1]
        use_ac = True
        if key in ac_nclose_cnt_dict.keys():
            ac_cnt = ac_nclose_cnt_dict[key]
            if ac_cnt > ABS_MAX_COVERAGE_RATIO * meandepth:
                use_ac = False
            for i in (0, 3):
                chrom = bd[i]
                coord = bd[i + 1]

                chrom_max = chr_max_end[chrom]

                assert(chrom in df_by_chr)
                df_chr = df_by_chr[chrom]

                df_bin = df_chr[(df_chr["st"] <= coord) & (coord < df_chr["nd"])]
                assert(not df_bin.empty)

                low_bound       = max(0, coord - FLANK_LENGTH)
                left_window_end = coord
                right_window_start = coord
                high_bound      = min(coord + FLANK_LENGTH, chrom_max)

                df_left = df_chr[(df_chr["nd"] > low_bound) & (df_chr["nd"] <= left_window_end)].copy()
                df_right = df_chr[(df_chr["st"] >= right_window_start) & (df_chr["st"] < high_bound)].copy()

                left_mean = df_left["meandepth"].mean() if not df_left.empty else 0.0
                right_mean = df_right["meandepth"].mean() if not df_right.empty else 0.0
                both_end_depth_dict[key].extend([round(left_mean, 2), round(right_mean, 2)])

                mean_diff = abs(left_mean - right_mean)
                if max(ac_cnt, 0) * BREAKEND_DEPTH_RATIO > mean_diff:
                    pos_fail_key_dict[key] = 'FAIL_DEPTH_NCLOSE'
                    use_ac = False
            
            if use_ac:
                using_ac_nclose_cnt_dict[key] = ac_nclose_cnt_dict[key] 

    return using_ac_nclose_cnt_dict


parser = argparse.ArgumentParser(description="SKYPE depth analysis")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("read_bam_loc", 
                    help="read_aln_loc")

parser.add_argument("censat_bed_path", 
                    help="Path to the censat repeat information file.")

parser.add_argument("reference_fai_path", 
                    help="Path to the chromosome information file.")

parser.add_argument("main_stat_path", 
                    help="Path to the main stat file.")

parser.add_argument("--progress", 
                    help="Show progress bar", action='store_true')

args = parser.parse_args()

PREFIX = args.prefix
read_bam_loc = args.read_bam_loc
reference_fai_path = args.reference_fai_path
main_stat_loc = args.main_stat_path
CENSAT_PATH = args.censat_bed_path
repeat_censat_data = import_censat_repeat_data(CENSAT_PATH)

with open(f"{PREFIX}/nclose_cord_list.pkl", "rb") as f:
    nclose_cord_list, nclose_idx_corr, total_nclose_cord_list_contig_name = pkl.load(f)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

pos_fail_key_dict = dict()
pre_fail_key_dict = dict()
both_end_depth_dict = defaultdict(list)

pre_nclose_cord_list = preprocess_breakends(nclose_cord_list, df, repeat_censat_data)

jl.nclose_cord_vec = jl.Vector[jl.Vector[jl.Any]]()
for l in pre_nclose_cord_list:
    jl.push_b(jl.nclose_cord_vec, jl.Vector[jl.Any](l))

is_progress_bar = sys.stdout.isatty() or args.progress
jl.include(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'anal_bam.jl'))
ac_nclose_cnt_list, wa_nclose_cnt_list = jl.anal_bam(read_bam_loc, reference_fai_path, jl.nclose_cord_vec, is_progress_bar)

ac_nclose_cnt_dict = dict(list(ac_nclose_cnt_list))
wa_nclose_cnt_dict = dict(list(wa_nclose_cnt_list))

post_ac_nclose_cnt_dict = postprocess_breakends(df, pre_nclose_cord_list, ac_nclose_cnt_dict)

nclose2cov = dict()
nclose_cov_target = set()
for k, ac_v in post_ac_nclose_cnt_dict.items():
    if k in wa_nclose_cnt_dict:
        wa_v = wa_nclose_cnt_dict[k]
        if ac_v / wa_v >= AC_WA_RATIO_LIMIT:
            nclose2cov[nclose_idx_corr[k]] = ac_v
            nclose_cov_target.add(k)
    else:
        nclose2cov[nclose_idx_corr[k]] = ac_v
        nclose_cov_target.add(k)

with open(f"{PREFIX}/nclose2cov.pkl", "wb") as f:
    pkl.dump(nclose2cov, f)

with open(f"{PREFIX}/nclose_cov_report.tsv", "wt") as f2:
    print(*['NCLOSE_CHR1', 'NCLOSE_CORD1', 'NCLOSE_DIR1', 'NCLOSE_CHR2', 'NCLOSE_CORD2', 'NCLOSE_DIR2', 'NCLOSE_ID',
            'ORIGIN_CONTIG_NAME', 'ACCEPT_COUNT', 'FAIL_COUNT', 'COV_TARGET', 'PREPROCESS_FAIL_CODE',
            'NCLOSE_DEPTH_LEFT1', 'NCLOSE_DEPTH_RIGHT1', 'NCLOSE_DEPTH_LEFT2', 'NCLOSE_DEPTH_RIGHT2', 'POSTPROCESS_FAIL_CODE'], sep="\t", file=f2)
    
    for l in total_nclose_cord_list_contig_name:
        k = l[-2]

        print_list = l + [ac_nclose_cnt_dict.get(k, -1), wa_nclose_cnt_dict.get(k, -1), 'COV' if k in nclose_cov_target else '*']
        print_list.append(pre_fail_key_dict.get(k, '*'))
        print_list.extend(both_end_depth_dict.get(k, ['*'] * 4))
        print_list.append(pos_fail_key_dict.get(k, '*'))

        print(*print_list, sep="\t", file=f2)

with open(f"{PREFIX}/bam_nclose_cnt.pkl", "wb") as f:
    pkl.dump((ac_nclose_cnt_dict, wa_nclose_cnt_dict), f)
