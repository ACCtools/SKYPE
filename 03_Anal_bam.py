import os
import sys
import logging
import argparse
import pickle as pkl
import pandas as pd
import numpy as np
import bisect

from collections import defaultdict, Counter

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

DIFF_COMPARE_RAITO = 5
SIM_COMPARE_RAITO = 1.2
HARD_COMPARE_RAITO = 1.1

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

def preprocess_breakends(nclose_cord_list, df, repeat_censat_data) -> tuple:
    pre_fail_key_dict = dict()
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

    not_using_key = set()

    for nclose_cord in nclose_cord_list:
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
                not_using_key.add(key)
                break

            df_bin = df_chr[(df_chr["st"] <= coord) & (coord < df_chr["nd"])]
            if not df_bin.empty and np.median(df_bin["meandepth"]) > ABS_MAX_COVERAGE_RATIO * meandepth:
                pre_fail_key_dict[key] = 'HIGH_DEPTH_REGION_NCLOSE'
                not_using_key.add(key)
                break
            
    
    using_nclose_list = []
    for nclose_cord in nclose_cord_list:
        key = nclose_cord[-1]
        if key not in not_using_key:
            using_nclose_list.append(nclose_cord)

    return using_nclose_list, pre_fail_key_dict


def postprocess_breakends(df, pre_cord_nclose_list) -> dict:
    both_end_depth_dict = defaultdict(list)
    # Pre‐group by chromosome for faster lookup
    df_by_chr = {chrom: subdf.reset_index(drop=True) for chrom, subdf in df.groupby("chr")}
    chr_max_end = {chrom: subdf["nd"].max() for chrom, subdf in df_by_chr.items()}

    for bd in pre_cord_nclose_list:
        key = bd[-1]
        both_dict_check = key in both_end_depth_dict
        
        for i in (0, 3):
            chrom = bd[i]
            coord = bd[i + 1]

            chrom_max = chr_max_end[chrom]

            assert(chrom in df_by_chr)
            df_chr = df_by_chr[chrom]

            df_bin = df_chr[(df_chr["st"] <= coord) & (coord < df_chr["nd"])]

            if not df_bin.empty:
                low_bound       = max(0, coord - FLANK_LENGTH)
                left_window_end = coord
                right_window_start = coord
                high_bound      = min(coord + FLANK_LENGTH, chrom_max)

                df_left = df_chr[(df_chr["nd"] > low_bound) & (df_chr["nd"] <= left_window_end)].copy()
                df_right = df_chr[(df_chr["st"] >= right_window_start) & (df_chr["st"] < high_bound)].copy()

                left_mean = df_left["meandepth"].mean() if not df_left.empty else 0.0
                right_mean = df_right["meandepth"].mean() if not df_right.empty else 0.0
                if not both_dict_check:
                    both_end_depth_dict[key].extend([round(left_mean, 2), round(right_mean, 2)])

    return both_end_depth_dict

def similar_check(v1, v2, ratio=SIM_COMPARE_RAITO):
    assert(v1 >= 0 and v2 >= 0)
    mi, ma = sorted([v1, v2])
    return False if mi == 0 else (ma / mi <= ratio)

def over_check(v1, v2):
    assert(v1 >= 0 and v2 >= 0)
    
    if v1 == 0:
        return False
    elif v2 == 0:
        return True
    else:
        return v1 / v2 >= DIFF_COMPARE_RAITO

parser = argparse.ArgumentParser(description="SKYPE depth analysis")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("read_bam_loc", 
                    help="Raw read alignment bam location")

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

with open(f"{PREFIX}/03_anal_bam_input.pkl", "rb") as f:
    nclose_cord_list, nclose_idx_corr, total_nclose_cord_list_contig_name, \
    total_dir_data, transloc_k_set, nclose_nodes = pkl.load(f)

task_dict = defaultdict(list)
nclose2cov = dict()
repeat_censat_data = import_censat_repeat_data(CENSAT_PATH)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

pre_nclose_cord_list, pre_fail_key_dict = preprocess_breakends(nclose_cord_list, df, repeat_censat_data)

task_cnt = Counter()
task_dict = dict()
run_k_set = set()

jl.nclose_cord_vec = jl.Vector[jl.Vector[jl.Any]]()
for l in pre_nclose_cord_list:
    run_k_set.add(l[-1])
    jl.push_b(jl.nclose_cord_vec, jl.Vector[jl.Any](l))

is_progress_bar = sys.stdout.isatty() or args.progress
jl.include(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'anal_bam.jl'))

nclose_cnt_list, wa_nclose_cnt_list = jl.anal_bam(read_bam_loc, reference_fai_path, jl.nclose_cord_vec, is_progress_bar)

dir_order_list = [(True, False), (False, True), (True, True), (False, False)]

nclose_cnt_dict = dict()
for i, tl in zip(dir_order_list, nclose_cnt_list):
    nclose_cnt_dict[i] = dict(list(tl))

wa_nclose_cnt_dict = dict(list(wa_nclose_cnt_list))

with open(f"{PREFIX}/bam_nclose_cnt.pkl", "wb") as f:
    pkl.dump((nclose_cnt_dict, wa_nclose_cnt_dict), f)

# with open(f"{PREFIX}/bam_nclose_cnt.pkl", "rb") as f:
#     nclose_cnt_dict, wa_nclose_cnt_dict = pkl.load(f)

both_end_depth_dict = postprocess_breakends(df, pre_nclose_cord_list)

for k in run_k_set:
    l1, r1, l2, r2 = both_end_depth_dict[k]
    wa_v = wa_nclose_cnt_dict.get(k, 0)
    
    is_zero = wa_v == 0
    for d in dir_order_list:
        v = nclose_cnt_dict[d].get(k, 0)
        if v > 0:
            is_zero = False
            break
    
    if is_zero \
       and similar_check(l1, r1, ratio=HARD_COMPARE_RAITO) \
       and similar_check(l2, r2, ratio=HARD_COMPARE_RAITO):
        
        ctg_data = total_dir_data[k][(True, False)]

        key = ctg_data[0]
        nclose_tuple = (ctg_data[1], ctg_data[2])
        nclose_nodes[key].remove(nclose_tuple)

        if nclose_nodes[key] == []:
            del nclose_nodes[key]

        task_dict[k] = 3
        task_cnt[task_dict[k]] += 1

    else:
        for d in total_dir_data[k]:
            if d != (True, False):
                v = nclose_cnt_dict[d].get(k, 0)
                if over_check(v, wa_v):
                    line = total_dir_data[k][d]

                    nclose_nodes[line[0]].append(tuple(sorted((int(line[1]), int(line[2])))))

                    task_dict[k] = 1
                    task_cnt[task_dict[k]] += 1
                    
        if k in transloc_k_set:
            ac_v = nclose_cnt_dict[(True, False)].get(k, 0)
            rac_v = nclose_cnt_dict[(False, True)].get(k, 0)

            if over_check(ac_v, wa_v) and over_check(rac_v, wa_v):
                if similar_check(ac_v, rac_v) and similar_check(l1, r1) and similar_check(l2, r2):
                    if min([ac_v, rac_v]) < min([l1, r1, l2, r2]):
                        rev_ctg_data = total_dir_data[k][(False, True)]

                        # Coverage constraint limit
                        nclose2cov[nclose_idx_corr[k]] = ac_v
                        nclose2cov[(rev_ctg_data[1], rev_ctg_data[2])] = rac_v

                        task_dict[k] = 2
                        task_cnt[task_dict[k]] += 2

with open(f"{PREFIX}/03_anal_bam_output.pkl", "wb") as f:
    pkl.dump((nclose_nodes, task_cnt), f)

with open(f"{PREFIX}/nclose_cov_report.tsv", "wt") as f2:
    print(*['NCLOSE_CHR1', 'NCLOSE_CORD1', 'NCLOSE_DIR1', 'NCLOSE_CHR2', 'NCLOSE_CORD2', 'NCLOSE_DIR2', 'NCLOSE_ID',
            'CONTIG_NAME', 'ACCEPT_COUNT', 'REV_COUNT', 'CORD1_REV_COUNT', 'CORD2_REV_COUNT', 'FAIL_COUNT', 'NCLOSE_TYPE', 'PREPROCESS_FAIL_CODE',
            'NCLOSE_DEPTH_LEFT1', 'NCLOSE_DEPTH_RIGHT1', 'NCLOSE_DEPTH_LEFT2', 'NCLOSE_DEPTH_RIGHT2'], sep="\t", file=f2)
    
    for l in total_nclose_cord_list_contig_name:
        k = l[-2]

        print_list = l
        if k in run_k_set:
            nclose_type = '*'
            if k in task_dict:
                if task_dict[k] == 1:
                    nclose_type = 'ADD_NCLOSE'
                elif task_dict[k] == 2:
                    nclose_type = 'TRANSLOC_NCLOSE'
                elif task_dict[k] == 3:
                    nclose_type = 'DEL_NCLOSE'

            print_list.extend([nclose_cnt_dict[(True, False)].get(k, 0), nclose_cnt_dict[(False, True)].get(k, 0),
                               nclose_cnt_dict[(False, False)].get(k, 0), nclose_cnt_dict[(True, True)].get(k, 0),
                               wa_nclose_cnt_dict.get(k, 0), nclose_type])
        else:
            print_list.extend([-1, -1, -1, -1, -1, '*'])
        print_list.append(pre_fail_key_dict.get(k, '*'))
        print_list.extend(both_end_depth_dict.get(k, ['*'] * 4))

        print(*print_list, sep="\t", file=f2)
            
with open(f"{PREFIX}/nclose2cov.pkl", "wb") as f:
    pkl.dump(nclose2cov, f)