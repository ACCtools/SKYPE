import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *
from nclose_tracking import (
    calculate_event_weights,
    carrier_columns,
    load_filter_status,
    load_path_usage,
    record_filter_stage,
    save_filter_status,
)
from path_additivity import choose_removal, find_path_conflicts

from collections import defaultdict

import ast
import csv
import glob
import psutil
import logging
import argparse

import numpy as np
import pandas as pd
import pickle as pkl

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
CTG_CENSAT = 18
CTG_MAINFLOWDIR = 19
CTG_MAINFLOWCHR = 20
CTG_GLOBALIDX = 21

CLUSTER_START_DEPTH = 0.1
CLUSTER_FINAL_DEPTH = 0.2

TELOMERE_EXPANSION = 5 * K

TYPE4_CLUSTER_SIZE = 1 * M

BASE_ACCSUMABSMAX_RATIO = 1.5
CHROM_ERROR_FAIL_RATE = 2.0
PRECLUSTER_NCLOSE_COUNT_RESULT_PKL = 'precluster_nclose_count_result.pkl'
RAW_ZERO_DECISION_TSV = 'cluster_raw_zero_path_decisions.tsv'

VCF_FLANKING_LENGTH = 1*M
NCLOSE_SIM_COMPARE_RAITO = 1.2
NCLOSE_SIM_DIFF_THRESHOLD = 5

def get_relative_path(p):
    return tuple(p.split('/')[-3:])

def import_ppc_data(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        curr_contig = curr_contig.rstrip()
        temp_list = curr_contig.split("\t")
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END, \
                          CHR_LEN, CHR_STR, CHR_END, \
                          CTG_MAPQ, CTG_TYP, CTG_STRND, CTG_ENDND,]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    paf_file.close()
    return contig_data

def chr2int(x):
    chrXY2int = {'chrX' : 24, 'chrY' : 25}
    if x in chrXY2int:
        return chrXY2int[x]
    else:
        return int(x[3:])

def import_telo_data(file_path : str, chr_len : dict) -> dict :
    fai_file = open(file_path, "r")
    telo_data = []
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        int_induce_idx = [1, 2]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        if temp_list[1]>chr_len[temp_list[0]]/2:
            temp_list[1]-=TELOMERE_EXPANSION
            temp_list.append('b')
        else:
            temp_list.append('f')
            temp_list[2]+=TELOMERE_EXPANSION
        telo_data.append(tuple(temp_list))
    fai_file.close()
    return telo_data

def extract_telomere_connect_contig(telo_info_path : str) -> list:
    telomere_connect_contig = []
    with open(telo_info_path) as f:
        for curr_data in f:
            curr_data = curr_data.rstrip()
            temp_list = curr_data.split("\t")
            chr_info = temp_list[0]
            contig_id = ast.literal_eval(temp_list[1])
            telomere_connect_contig.append((chr_info, contig_id[1]))
    
    return telomere_connect_contig

def distance_checker(node_a : tuple, node_b : tuple) -> int :
    if max(int(node_a[0]), int(node_b[0])) < min(int(node_a[1]), int(node_b[1])):
        return 0   
    else:
        return min(abs(int(node_b[0]) - int(node_a[1])), abs(int(node_b[1]) - int(node_a[0])))

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    chr_data_file.close()
    return chr_len

def similar_check(v1, v2, ratio):
    if not (v1 >= 0 and v2 >= 0):
        logging.error(f"Invalid values for similarity check: v1={v1}, v2={v2}")
        return False
    mi, ma = sorted([v1, v2])
    return True if mi == 0 else (ma / mi <= ratio) or ma-mi < NCLOSE_SIM_DIFF_THRESHOLD

def exist_near_bnd(chrom, inside_st, inside_nd):
    df_chr = df[df['chr'] == chrom]

    def mean_depth(start, end):
        mask = (df_chr['nd'] > start) & (df_chr['st'] < end)
        return df_chr.loc[mask, 'meandepth'].mean()

    st_depth = mean_depth(inside_st - VCF_FLANKING_LENGTH, inside_st)
    nd_depth = mean_depth(inside_nd, inside_nd + VCF_FLANKING_LENGTH)
    if np.isnan(st_depth) or np.isnan(nd_depth):
        return True

    return not similar_check(st_depth, nd_depth, NCLOSE_SIM_COMPARE_RAITO)


def nclose_count_valid(st, ed):
    """
    Count "valid" sides:
      - censat side: always treated as valid
      - non-censat side: valid iff exist_near_bnd is True
    Returns 0, 1, or 2.
    """
    def side_valid(idx):
        if ppc_data[idx][CTG_CENSAT] != '0':
            return True
        return exist_near_bnd(
            ppc_data[idx][CHR_NAM],
            ppc_data[idx][CHR_STR],
            ppc_data[idx][CHR_END],
        )
    st_cnt = int(side_valid(st))
    ed_cnt = int(side_valid(ed))

    return st_cnt + ed_cnt

def is_filter_candidate(st, ed, tier_st, tier_ed):
    """
    Decide if nclose (st, ed) enters the filter set.

    Per-side tier verdict:
      Tier 1: filter if count_valid <= 1 (at least one side invalid)
      Tier 2: filter only if count_valid == 0 (both sides invalid)
      Tier 3: never filter
    nclose enters the set only when both sides' tiers vote filter (AND).
    """
    if ppc_data[st][CTG_NAM].startswith('virtual_censat_contig'):
        return False

    if tier_st >= 3 or tier_ed >= 3:
        return False

    count_valid = nclose_count_valid(st, ed)

    def tier_says_filter(tier):
        if tier == 1:
            return count_valid <= 1
        if tier == 2:
            return count_valid == 0
        return False

    return tier_says_filter(tier_st) and tier_says_filter(tier_ed)

def load_high_quality_zero_read_events(prefix, contig_data):
    result_path = os.path.join(prefix, PRECLUSTER_NCLOSE_COUNT_RESULT_PKL)
    if not os.path.isfile(result_path):
        raise FileNotFoundError(
            f'Missing pre-cluster NClose read counts: {result_path}. '
            'Run 25_cluster_nclose_read_count.py with --selection_stage base '
            'before stage 24.'
        )

    with open(result_path, 'rb') as f:
        records = pkl.load(f)

    selected = {}
    for record in records:
        query_gap = record.get('nclose_query_gap')
        if query_gap is None or int(query_gap) > 5 * K:
            continue
        if int(record.get('read_counts', {}).get('d2', 0)) != 0:
            continue
        if bool(record.get('overlaps_censat', False)):
            continue
        try:
            node_pair = tuple(sorted(int(index) for index in record['nclose_key']))
        except (KeyError, TypeError, ValueError):
            continue
        if len(node_pair) != 2 or max(node_pair) >= len(contig_data):
            continue
        mapqs = tuple(int(contig_data[index][CTG_MAPQ]) for index in node_pair)
        if mapqs != (60, 60):
            continue
        event_key = record.get('event_key')
        if event_key is None:
            event_key = node_pair
        else:
            event_key = tuple(event_key)
        selected[event_key] = {
            'event_key': event_key,
            'node_pair': node_pair,
            'nclose_id': record.get('nclose_id', ''),
            'contig_name': record.get('ctg_name', ''),
            'query_gap': int(query_gap),
            'mapq_a': mapqs[0],
            'mapq_b': mapqs[1],
            'normal_reads_a': int(record.get('read_counts', {}).get('d1', 0)),
            'junction_reads': 0,
            'normal_reads_b': int(record.get('read_counts', {}).get('d3', 0)),
        }
    return selected

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).

parser = argparse.ArgumentParser()

parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")

parser.add_argument("main_stat_loc", 
                    help="Cancer coverage location file")

parser.add_argument("telomere_bed_path", 
                    help="Path to the telomere information file.")

parser.add_argument("reference_fai_path", 
                    help="Path to the chromosome information file.")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("-t", "--thread", help="Number of threads", type=int)

args = parser.parse_args()

# t = """
# 24_cluster_weight.py /Data/hyunwoo/00_skype_run_data/COLO829/20_alignasm/COLO829.ctg.aln.paf.ppc.paf /Data/hyunwoo/00_skype_run_data/COLO829/01_depth/COLO829_normalized.win.stat.gz /hyunwoo/63_skype_test/deps/SKYPE/public_data/chm13v2.0_telomere.bed /hyunwoo/63_skype_test/deps/SKYPE/public_data/chm13v2.0.fa.fai /hyunwoo/63_skype_test/skype/COLO829_15_06_54 -t 16
# """
# args = parser.parse_args(t.strip().split()[1:])

PREFIX = args.prefix
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
main_stat_loc = args.main_stat_loc
TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
pipeline_mode_config = load_pipeline_mode(PREFIX)

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
TELO_CONNECT_NODES_INFO_PATH = PREFIX+"/telomere_connected_list.txt"
front_contig_path = RATIO_OUTLIER_FOLDER+"front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER+"back_jump/"
ecdna_contig_path = RATIO_OUTLIER_FOLDER+"ecdna/"
output_folder = f'{PREFIX}/21_pat_depth'
NCLOSE_FILE_PATH = f"{args.prefix}/nclose_nodes_index.txt"

use_julia_solver = pipeline_mode_is_karyotype(pipeline_mode_config)

if not use_julia_solver:
    exit(0)
logging.info("24_cluster_weight start")

THREAD = args.thread
core_num = psutil.cpu_count(logical=False)
if core_num is None:
    THREAD = THREAD
else:
    THREAD = min(int(THREAD), core_num)

NNLS_GRAM_MEM_RATIO = 0.8
BYTES_PER_GB = 1e9

def nnls_gram_mem_gb():
    return psutil.virtual_memory().available * NNLS_GRAM_MEM_RATIO / BYTES_PER_GB

ppc_data = import_ppc_data(PREPROCESSED_PAF_FILE_PATH)

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    paf_ans_list, key_list, int2key, _ = pkl.load(f)
paf_ans_dict = dict(paf_ans_list)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

meandepth = np.median(df['meandepth'])
N = meandepth / 2

chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)

# Julia remove small nclose
with open(f'{PREFIX}/23_input.pkl', 'rb') as f:
    chr_filt_st_list, path_nclose_set_dict, amplitude, bed_data, matrix_depth_meta = pkl.load(f)

ydf = df.query('chr == "chrY"')

bed_intervals = pd.IntervalIndex.from_tuples(bed_data['chrY'], closed='left')
y_mask = ydf.apply(
    lambda row: bed_intervals.overlaps(pd.Interval(row['st'], row['nd'], closed='left')),
    axis=1
)

correct_mask = y_mask.apply(any)
ydf_not_censat = ydf[~correct_mask]

chry_nz_len = len(ydf_not_censat.query('meandepth != 0'))
no_chrY = (chry_nz_len / len(ydf_not_censat)) < chrY_MINIMUM_RATIO

from juliacall import Main as jl

jl.include(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'nnls_solver.jl'))
jl.seval("using HDF5, LinearAlgebra")

A_jl, B_jl, b_start_ind = jl.load_nnls_array(f'{PREFIX}/matrix.h5')
A_fail_jl = jl.load_fail_array(f'{PREFIX}/matrix.h5')
B = np.asarray(B_jl)

b_start_ind = int(matrix_depth_meta.get("B_depth_start", b_start_ind))
b_depth_end = int(matrix_depth_meta.get("B_depth_end", b_start_ind + len(chr_filt_st_list)))

depth_success_slice = slice(b_start_ind, b_depth_end)
B_depth = B[depth_success_slice]
if len(B_depth) != len(chr_filt_st_list):
    raise ValueError(
        f"Depth row count mismatch: matrix has {len(B_depth)} clean depth rows, "
        f"23_input has {len(chr_filt_st_list)}"
    )

def make_jl_weight(weight):
    weight = np.maximum(np.asarray(weight), 0)
    return jl.Vector[jl.eltype(B_jl)](weight)

def run_nnls_map_gram_mem_gb(task_count):
    # run_nnls_map applies gram_mem_gb per spawned subproblem, so split the
    # 80% total budget across the maximum number of concurrently active tasks.
    active_tasks = min(max(int(task_count), 1), max(int(THREAD), 1))
    return nnls_gram_mem_gb() / active_tasks

# Get weight from 23_run_nnls.py
weight_base = np.load(f'{PREFIX}/weight.npy')
weight_base_jl = make_jl_weight(weight_base)
nclose_path_usage = load_path_usage(PREFIX, expected_len=len(weight_base))
nclose_filter_status = load_filter_status(PREFIX)
if 'base' not in nclose_filter_status['stages']:
    raise ValueError(
        "Missing stage-23 NClose filter provenance. Rerun SKYPE from stage 23."
    )

final_weights_fullsize = np.zeros(len(weight_base))
predict_suc_B_base = np.asarray(A_jl * weight_base_jl)

chrom_acc_sum_dict_base = defaultdict(int)
chrom_acc_sum_dict_max_base = defaultdict(int)
for i, (chrom, st) in enumerate(chr_filt_st_list):
    chrom_acc_sum_dict_base[chrom] += predict_suc_B_base[i] - B[i]
    if abs(chrom_acc_sum_dict_base[chrom]) > chrom_acc_sum_dict_max_base[chrom]:
        chrom_acc_sum_dict_max_base[chrom] = abs(chrom_acc_sum_dict_base[chrom])

nclose_total_weight_dict = defaultdict(
    float,
    calculate_event_weights(nclose_path_usage, weight_base),
)

with open(f'{PREFIX}/A_idx_list.pkl', 'rb') as f:
    A_idx_list = pkl.load(f)

def chrom_acc_sum_peak_stats(predict_suc_B):
    chrom_acc_sum_dict = defaultdict(float)
    chrom_peak_stats = {}
    for i, (chrom, st) in enumerate(chr_filt_st_list):
        chrom_acc_sum_dict[chrom] += float(predict_suc_B[i] - B[i])
        abs_acc_sum = abs(chrom_acc_sum_dict[chrom])
        if abs_acc_sum > chrom_peak_stats.get(chrom, {}).get('max_abs', -1):
            chrom_peak_stats[chrom] = {
                'max_abs': abs_acc_sum,
                'signed': chrom_acc_sum_dict[chrom],
                'idx': i,
                'coord': int(st),
            }
    return chrom_peak_stats

def get_fit_failures(
    predict_suc_B,
    fail_rate=CHROM_ERROR_FAIL_RATE,
    baseline_peak_stats=None,
):
    failures = []
    for chrom, stat in chrom_acc_sum_peak_stats(predict_suc_B).items():
        if chrom == 'chrY' and no_chrY:
            continue
        if baseline_peak_stats is None:
            acc_sum_max_base = chrom_acc_sum_dict_max_base[chrom]
        else:
            acc_sum_max_base = baseline_peak_stats.get(chrom, {}).get('max_abs', 0)
        # Preserve the previous 24_cluster_weight behavior: chromosomes with
        # a zero baseline envelope do not drive greedy rejection.
        if acc_sum_max_base == 0:
            continue
        ratio = stat['max_abs'] / acc_sum_max_base
        if ratio > fail_rate:
            failures.append({
                'chrom': chrom,
                'ratio': ratio,
                'coord': stat['coord'],
                'signed': stat['signed'],
                'max_abs': stat['max_abs'],
            })
    failures.sort(key=lambda x: x['ratio'], reverse=True)
    return failures

def solution_depth_error(solution):
    predict_depth = solution['predict_succ'][depth_success_slice]
    return float(np.linalg.norm(predict_depth - B_depth))

def solve_idx_list(final_idx_list, warm_fullsize):
    final_idx_array_jl = jl.Vector[jl.Int]([i + 1 for i in final_idx_list])

    A_final_jl = jl.view(A_jl, jl.Colon(), final_idx_array_jl)
    A_fail_final_jl = jl.view(A_fail_jl, jl.Colon(), final_idx_array_jl)

    solve_gram_mem_gb = nnls_gram_mem_gb()
    final_weight_jl = jl.nnls_solve(
        A_final_jl, B_jl, THREAD, False,
        gram_mem_gb=solve_gram_mem_gb,
        w_init=make_jl_weight(warm_fullsize[final_idx_list])
    )
    final_weight = np.asarray(final_weight_jl)

    final_weights_fullsize = np.zeros(len(weight_base), dtype=weight_base.dtype)
    for i, v in enumerate(final_idx_list):
        final_weights_fullsize[int(v)] = final_weight[i]

    predict_B_succ = np.asarray(A_final_jl * final_weight_jl)
    predict_B_fail = np.asarray(A_fail_final_jl * final_weight_jl)

    return {
        'idx_list': final_idx_list,
        'weight_jl': final_weight_jl,
        'weight': final_weight,
        'fullsize': final_weights_fullsize,
        'predict_succ': predict_B_succ,
        'predict_fail': predict_B_fail,
    }

# Pre-calculate max error rates for all possible candidate nclose pairs
nclose_max_error_dict = {}
pre_tar_nclose_list = []
pre_thread_data_jl = jl.Vector[jl.Vector[jl.Int]]()
base_A_idx_set = set(A_idx_list)

# Collect all potential candidate pairs using Tier 1 (most permissive — covers all tiers)
for k, v in nclose_total_weight_dict.items():
    if len(k) == 2:
        st, ed = k

        if v >= 0.01 * meandepth and ppc_data[st][CHR_NAM] != ppc_data[ed][CHR_NAM]:
            if is_filter_candidate(st, ed, 1, 1):
                pre_tar_nclose_list.append(k)

pre_nclose_notusing_idx_dict = defaultdict(list)
pre_solved_nclose_list = []
for nclose_pair in pre_tar_nclose_list:
    for path_k, event_usage in enumerate(nclose_path_usage):
        if path_k in base_A_idx_set and event_usage.get(nclose_pair, 0) == 0:
            pre_nclose_notusing_idx_dict[nclose_pair].append(path_k)
    if not pre_nclose_notusing_idx_dict[nclose_pair]:
        nclose_max_error_dict[nclose_pair] = float('inf')
        continue
    jl.push_b(pre_thread_data_jl, jl.Vector[jl.Int](pre_nclose_notusing_idx_dict[nclose_pair]))
    pre_solved_nclose_list.append(nclose_pair)

logging.info(f"Pre-calculating max error rate for pairs count : {len(pre_solved_nclose_list)}")
pre_map_gram_mem_gb = run_nnls_map_gram_mem_gb(len(pre_solved_nclose_list))

if pre_solved_nclose_list:
    pre_weight_list = jl.run_nnls_map(
        A_jl, B_jl, pre_thread_data_jl, False,
        gram_mem_gb=pre_map_gram_mem_gb,
        w_init=weight_base_jl
    )
else:
    pre_weight_list = []

# Store pre-calculated max error rates in dictionary
for nclose_pair, (weight_jl, predict_suc_B_jl) in zip(pre_solved_nclose_list, pre_weight_list):
    predict_suc_B = np.asarray(predict_suc_B_jl)

    chrom_acc_sum_dict = defaultdict(int)
    chrom_acc_sum_dict_max = defaultdict(int)
    for i, (chrom, st) in enumerate(chr_filt_st_list):
        chrom_acc_sum_dict[chrom] += predict_suc_B[i] - B[i]
        if abs(chrom_acc_sum_dict[chrom]) > chrom_acc_sum_dict_max[chrom]:
            chrom_acc_sum_dict_max[chrom] = abs(chrom_acc_sum_dict[chrom])

    max_error_rate = 0.0
    for chrom, acc_sum_max in chrom_acc_sum_dict_max.items():
        if chrom == 'chrY' and no_chrY:
            continue
        acc_sum_max_base = chrom_acc_sum_dict_max_base[chrom]
        # Calculate error rate only if base is not zero
        if acc_sum_max_base != 0:
            max_error_rate = max(max_error_rate, acc_sum_max / acc_sum_max_base)
    
    nclose_max_error_dict[nclose_pair] = max_error_rate

# Greedy filtering:
#   1. keep the depth-derived Tier-1 candidate pool,
#   2. order candidates by their single-removal error,
#   3. add one nclose at a time and refit,
#   4. if the trial fails, reject the current candidate.
nclose_notusing_idx_set_dict = {
    nclose_pair: set(pre_nclose_notusing_idx_dict[nclose_pair]) & base_A_idx_set
    for nclose_pair in pre_tar_nclose_list
}

def solve_removed_nclose(removed_nclose_set, warm_fullsize):
    using_idx_set = set(base_A_idx_set)
    for nclose_pair in removed_nclose_set:
        using_idx_set &= nclose_notusing_idx_set_dict[nclose_pair]
        if not using_idx_set:
            break
    if not using_idx_set:
        return None
    return solve_idx_list(sorted(using_idx_set), warm_fullsize)

greedy_candidate_list = sorted(
    {
        nclose_pair
        for nclose_pair in pre_tar_nclose_list
        if nclose_max_error_dict.get(nclose_pair, float('inf')) <= BASE_ACCSUMABSMAX_RATIO
    },
    key=lambda nclose_pair: (
        nclose_max_error_dict.get(nclose_pair, float('inf')),
        nclose_pair,
    )
)

logging.info(
    f'Greedy nclose filtering start: candidates={len(greedy_candidate_list)} '
    f'(precomputed pool={len(pre_tar_nclose_list)}, single-error <= {BASE_ACCSUMABSMAX_RATIO})'
)

accepted_nclose = set()
rejected_nclose = set()
filter_warm_fullsize = weight_base.copy()
current_solution = solve_removed_nclose(accepted_nclose, filter_warm_fullsize)
if current_solution is None:
    raise RuntimeError("No NNLS columns available before greedy nclose filtering.")
filter_warm_fullsize = current_solution['fullsize'].copy()

for nclose_pair in greedy_candidate_list:
    if nclose_pair in rejected_nclose:
        continue

    trial_nclose = set(accepted_nclose)
    trial_nclose.add(nclose_pair)
    trial_solution = solve_removed_nclose(trial_nclose, filter_warm_fullsize)
    if trial_solution is None:
        rejected_nclose.add(nclose_pair)
        logging.debug(f'{nclose_pair} : rejected because no NNLS columns remain')
        continue

    failures = get_fit_failures(trial_solution['predict_succ'])
    if not failures:
        accepted_nclose = trial_nclose
        current_solution = trial_solution
        filter_warm_fullsize = current_solution['fullsize'].copy()
        logging.debug(
            f'{nclose_pair} : greedy accepted '
            f'(single_error={nclose_max_error_dict.get(nclose_pair, float("inf")):.4f})'
        )
        continue

    worst = failures[0]
    rejected_nclose.add(nclose_pair)
    logging.debug(
        f'{nclose_pair} : greedy rejected '
        f'(worst={worst["chrom"]}:{worst["coord"]}, '
        f'ratio={worst["ratio"]:.3f})'
    )

not_essential_nclose = accepted_nclose
final_idx_list = current_solution['idx_list']
final_weight_jl = current_solution['weight_jl']
final_weight = current_solution['weight']
final_weights_fullsize = current_solution['fullsize']
predict_B_succ = current_solution['predict_succ']
predict_B_fail = current_solution['predict_fail']

logging.info(
    f'Filtered nclose count by greedy depth : {len(not_essential_nclose)} '
    f'(rejected={len(rejected_nclose)})'
)
for nclose_pair in sorted(not_essential_nclose):
    st, ed = nclose_pair
    st_row = ppc_data[st]
    ed_row = ppc_data[ed]
    logging.debug(
        f'  Filtered nclose: '
        f'{st_row[CTG_NAM]} ({st_row[CHR_NAM]}:{st_row[CHR_STR]}-{st_row[CHR_END]}) '
        f'<-> '
        f'{ed_row[CTG_NAM]} ({ed_row[CHR_NAM]}:{ed_row[CHR_STR]}-{ed_row[CHR_END]})'
    )

b_norm = np.linalg.norm(B_depth)

predict_B_succ_depth = predict_B_succ[depth_success_slice]
error = np.linalg.norm(predict_B_succ_depth - B_depth)
predict_B = np.concatenate((predict_B_succ_depth, predict_B_fail))

logging.info(f'Filter Error : {error:.4f}')
logging.info(f'Filter relative error : {error/b_norm:.4f}')

np.save(f'{PREFIX}/weight_filter.npy', final_weights_fullsize)
np.save(f'{PREFIX}/predict_B_filter.npy', predict_B)

filter_direct_reasons = {
    event_key: 'FILTERED_24_GREEDY_DEPTH'
    for event_key in not_essential_nclose
}
nclose_filter_status['stages'].pop('cluster', None)
record_filter_stage(
    nclose_filter_status,
    'filter',
    'base',
    nclose_path_usage,
    final_idx_list,
    filter_direct_reasons,
    'FILTERED_24_COFILTERED_PATH',
)
save_filter_status(PREFIX, nclose_filter_status)

weights_sorted_data = sorted(enumerate(final_weights_fullsize), key=lambda t:t[1], reverse=True)

chr_inf = max(chr_len.values())
chr_fb_len_dict = defaultdict(list)

telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)
telo_connected_node = extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH)

telo_dict = defaultdict(list)
for _ in telo_data:
    telo_dict[_[0]].append(_[1:])

telo_fb_dict = defaultdict(list)
for k, v in telo_dict.items():
    for i in v:
        telo_fb_dict[k+i[-1]].append([i[0], i[1]])


for chr_dir, node_id in telo_connected_node:
    telo_len = chr_inf
    for telo_bed in telo_fb_dict[chr_dir]:
        telo_len = min(telo_len, distance_checker(tuple(telo_bed), (ppc_data[node_id][CHR_STR], ppc_data[node_id][CHR_END])))
    chr_fb_len_dict[chr_dir].append((node_id, telo_len, chr_dir))

telo_len_data = []
for chr_dir, telo_len_list in chr_fb_len_dict.items():
    s_telo_len_list = sorted(telo_len_list, key=lambda t: t[1])
    telo_len_data.extend(filter(lambda t: t[1] > 0, s_telo_len_list[1:]))

need_label = defaultdict(list)
need_label_index = dict()
for node_id, telo_len, chr_dir in telo_len_data:
    need_label[chr_dir[:-1]].append((node_id, chr_dir[-1]))
    need_label_index[node_id] = (chr_dir, telo_len)


tot_loc_list = []
for loc, ll in paf_ans_list:
    tot_loc_list.append(loc)

fclen = len(glob.glob(front_contig_path+"*"))
bclen = len(glob.glob(back_contig_path+"*"))
eclen = len(glob.glob(ecdna_contig_path+"*"))

for i in range(1, fclen//4 + 1):
    bv_paf_loc = front_contig_path+f"{i}_base.paf"
    tot_loc_list.append(bv_paf_loc)

for i in range(1, bclen//4 + 1):
    bv_paf_loc = back_contig_path+f"{i}_base.paf"
    tot_loc_list.append(bv_paf_loc)

for i in range(1, eclen//2 + 1):
    ecdna_paf_loc = ecdna_contig_path+f"{i}.paf"
    tot_loc_list.append(ecdna_paf_loc)

with open(f'{PREFIX}/cen_fragment_data.pkl', 'rb') as f:
    cen_fragment_meta = pkl.load(f)

cen_fragment_list = sorted(cen_fragment_meta.items(), key=lambda kv: chr2int(kv[0]))
for chrom, info in cen_fragment_list:
    side = 'right' if info["dir"] else 'left'
    tot_loc_list.append(f'{PREFIX}/12_cent_fragment/{chrom}/{side}.fragment')

if len(tot_loc_list) != len(final_weights_fullsize):
    raise ValueError(
        "Stage-24 matrix-column metadata mismatch: "
        f"{len(tot_loc_list)} locations for {len(final_weights_fullsize)} weights."
    )
loc2weight = dict(zip(tot_loc_list, final_weights_fullsize))
filter_active_column_set = set(final_idx_list)
filter_active_locations = {
    tot_loc_list[col_idx]
    for col_idx in filter_active_column_set
}

# Build the initial stage-24 column pool without the former sequence-similarity
# clustering. Ordinary paths retain the existing 0.1N entry threshold, while
# pure chromosomes and the special type-4/centromere columns keep their prior
# inclusion rules.
initial_cluster_columns = set()
loc2column = {location: column for column, location in enumerate(tot_loc_list)}

# Add pure chromosomes even when their fitted weight is below 0.1N.
with open(f"{PREFIX}/tar_chr_data.pkl", "rb") as f:
    tar_chr_data = pkl.load(f)

loc_prefix = paf_ans_list[0][0].split('/')[:-3]

for path_tuple in tar_chr_data.values():
    path = '/'.join(loc_prefix + list(path_tuple))
    if path in filter_active_locations:
        initial_cluster_columns.add(loc2column[path])

for ind, w in weights_sorted_data:
    paf_loc = tot_loc_list[ind]
    if paf_loc.split('/')[-3] not in {'11_ref_ratio_outliers', '12_cent_fragment'}:
        if w > CLUSTER_START_DEPTH * N:
            initial_cluster_columns.add(ind)

for ncnt, (paf_loc, w) in enumerate(loc2weight.items()):
    if paf_loc.split('/')[-3] == '11_ref_ratio_outliers' and w > CLUSTER_FINAL_DEPTH * N:
        event_type = paf_loc.split('/')[-2]

        if event_type in {'front_jump', 'back_jump'}:
            with open(paf_loc, "r") as f:
                l = f.readline()
                l = l.rstrip()
                l = l.split("\t")
                pos1 = int(l[CHR_STR])
                pos2 = int(l[CHR_END])
                if abs(pos1-pos2) > TYPE4_CLUSTER_SIZE:
                    initial_cluster_columns.add(ncnt)
        elif event_type == 'ecdna':
            initial_cluster_columns.add(ncnt)
        else:
            assert(False)
    elif paf_loc.split('/')[-3] == '12_cent_fragment':
        initial_cluster_columns.add(ncnt)

initial_cluster_columns &= filter_active_column_set

# Raw-read veto runs before compositional A/B/AB conflict resolution.  A zero
# is actionable only when both PAF sides are non-CENSAT, MAPQ 60, and the
# 5-kb junction count was technically available.
raw_zero_event_records = load_high_quality_zero_read_events(PREFIX, ppc_data)
raw_zero_removed_columns = set()
raw_zero_applied_event_keys = set()
raw_zero_decision_rows = []
for event_key, record in sorted(
    raw_zero_event_records.items(), key=lambda item: repr(item[0])
):
    carrier_set = carrier_columns(nclose_path_usage, event_key)
    removed_carriers = carrier_set & initial_cluster_columns
    if removed_carriers:
        raw_zero_applied_event_keys.add(event_key)
        raw_zero_removed_columns.update(removed_carriers)
    raw_zero_decision_rows.append({
        **record,
        'event_key': repr(event_key),
        'carrier_columns': ','.join(map(str, sorted(carrier_set))),
        'removed_cluster_columns': ','.join(map(str, sorted(removed_carriers))),
        'removed_cluster_paths': ';'.join(
            '/'.join(get_relative_path(tot_loc_list[column]))
            for column in sorted(removed_carriers)
        ),
        'action': 'REMOVE_AND_REFIT' if removed_carriers else 'NO_CLUSTER_CARRIER',
    })

raw_zero_decision_fields = [
    'nclose_id', 'event_key', 'node_pair', 'contig_name', 'query_gap',
    'mapq_a', 'mapq_b',
    'normal_reads_a', 'junction_reads', 'normal_reads_b',
    'carrier_columns', 'removed_cluster_columns', 'removed_cluster_paths',
    'action',
]
with open(f'{PREFIX}/{RAW_ZERO_DECISION_TSV}', 'w', newline='') as f:
    writer = csv.DictWriter(
        f, fieldnames=raw_zero_decision_fields, delimiter='\t',
        extrasaction='ignore',
    )
    writer.writeheader()
    writer.writerows(raw_zero_decision_rows)

initial_cluster_columns -= raw_zero_removed_columns
logging.info(
    'High-quality zero-read filtering before A/B/AB: events=%d, '
    'removed_columns=%d',
    len(raw_zero_applied_event_keys), len(raw_zero_removed_columns),
)
if not initial_cluster_columns:
    raise RuntimeError("No NNLS columns remain for stage-24 cluster selection.")

current_solution = solve_idx_list(
    sorted(initial_cluster_columns),
    final_weights_fullsize,
)
if current_solution is None:
    raise RuntimeError("No NNLS solution for the initial stage-24 column pool.")

logging.info(
    f'Compositional path filtering start: initial_columns={len(initial_cluster_columns)}, '
    f'initial_error={solution_depth_error(current_solution):.4f}'
)

def is_ordinary_path_column(column):
    return tot_loc_list[column].split('/')[-3] not in {
        '11_ref_ratio_outliers', '12_cent_fragment'
    }

def path_conflict_candidate_columns(solution):
    return {
        int(column)
        for column in solution['idx_list']
        if is_ordinary_path_column(int(column))
        and solution['fullsize'][int(column)] > CLUSTER_START_DEPTH * N
    }

def evaluate_path_removal(columns, baseline_peak_stats):
    remaining = sorted(set(current_solution['idx_list']) - set(columns))
    if not remaining:
        return {
            'solution': None,
            'error': float('inf'),
            'failures': [{'chrom': 'all', 'ratio': float('inf')}],
            'feasible': False,
        }
    solution = solve_idx_list(remaining, current_solution['fullsize'])
    failures = get_fit_failures(
        solution['predict_succ'],
        baseline_peak_stats=baseline_peak_stats,
    )
    return {
        'solution': solution,
        'error': solution_depth_error(solution),
        'failures': failures,
        'feasible': not failures,
    }

def failure_summary(evaluation):
    if not evaluation['failures']:
        return ''
    worst = evaluation['failures'][0]
    return f"{worst['chrom']}:{worst.get('coord', -1)}:{worst['ratio']:.6g}"

path_conflict_decisions = []
path_conflict_removed_columns = set()
baseline_peak_stats = chrom_acc_sum_peak_stats(current_solution['predict_succ'])

while True:
    candidate_columns = path_conflict_candidate_columns(current_solution)
    conflicts = find_path_conflicts(
        candidate_columns,
        nclose_path_usage,
        current_solution['fullsize'],
    )
    if not conflicts:
        break

    # Both A/B/AB and A/BC/AB enter this one queue, so the globally smallest
    # current three-path weight sum is always resolved first.
    conflict = conflicts[0]
    eval_a = evaluate_path_removal((conflict.path_a,), baseline_peak_stats)
    eval_b = evaluate_path_removal((conflict.path_b,), baseline_peak_stats)
    eval_both = None
    if eval_a['feasible'] and eval_b['feasible']:
        eval_both = evaluate_path_removal(
            (conflict.path_a, conflict.path_b),
            baseline_peak_stats,
        )

    choice = choose_removal(
        conflict.path_a,
        conflict.path_b,
        a_feasible=eval_a['feasible'],
        b_feasible=eval_b['feasible'],
        a_error=eval_a['error'],
        b_error=eval_b['error'],
        both_feasible=(eval_both['feasible'] if eval_both is not None else None),
    )
    evaluation_by_removal = {
        (conflict.path_a,): eval_a,
        (conflict.path_b,): eval_b,
    }
    if eval_both is not None:
        evaluation_by_removal[
            tuple(sorted((conflict.path_a, conflict.path_b)))
        ] = eval_both
    selected_evaluation = evaluation_by_removal[choice.removed_columns]
    if selected_evaluation['solution'] is None:
        raise RuntimeError(
            f"Selected compositional path removal has no NNLS solution: {choice}"
        )

    decision_index = len(path_conflict_decisions) + 1
    path_conflict_decisions.append({
        'iteration': decision_index,
        'conflict_type': conflict.conflict_type,
        'path_a_role': conflict.role_a,
        'path_b_role': conflict.role_b,
        'path_ab_role': conflict.role_ab,
        'path_a_column': conflict.path_a,
        'path_b_column': conflict.path_b,
        'path_ab_column': conflict.path_ab,
        'path_a_location': '/'.join(get_relative_path(tot_loc_list[conflict.path_a])),
        'path_b_location': '/'.join(get_relative_path(tot_loc_list[conflict.path_b])),
        'path_ab_location': '/'.join(get_relative_path(tot_loc_list[conflict.path_ab])),
        'weight_a_N': conflict.weight_a / N,
        'weight_b_N': conflict.weight_b / N,
        'weight_ab_N': conflict.weight_ab / N,
        'weight_sum_N': conflict.total_weight / N,
        'remove_a_feasible': eval_a['feasible'],
        'remove_b_feasible': eval_b['feasible'],
        'remove_both_tested': eval_both is not None,
        'remove_both_feasible': eval_both['feasible'] if eval_both is not None else '',
        'remove_a_error': eval_a['error'],
        'remove_b_error': eval_b['error'],
        'remove_both_error': eval_both['error'] if eval_both is not None else '',
        'remove_a_worst_failure': failure_summary(eval_a),
        'remove_b_worst_failure': failure_summary(eval_b),
        'remove_both_worst_failure': (
            failure_summary(eval_both) if eval_both is not None else ''
        ),
        'selected_removed_columns': ','.join(map(str, choice.removed_columns)),
        'selection_reason': choice.reason,
        'selected_error': selected_evaluation['error'],
    })

    logging.info(
        f'Compositional path decision {decision_index} '
        f'[{conflict.conflict_type}]: '
        f'{conflict.role_a}={conflict.path_a}({conflict.weight_a/N:.3f}N), '
        f'{conflict.role_b}={conflict.path_b}({conflict.weight_b/N:.3f}N), '
        f'{conflict.role_ab}={conflict.path_ab}({conflict.weight_ab/N:.3f}N), '
        f'removed={choice.removed_columns}, reason={choice.reason}, '
        f'error={selected_evaluation["error"]:.4f}'
    )

    path_conflict_removed_columns.update(choice.removed_columns)
    current_solution = selected_evaluation['solution']
    # Each completed A/B/AB or A/BC/AB decision defines the comparison
    # envelope for the next decision. No later trial uses a stale global base.
    baseline_peak_stats = chrom_acc_sum_peak_stats(current_solution['predict_succ'])

decision_fields = [
    'iteration', 'conflict_type',
    'path_a_role', 'path_b_role', 'path_ab_role',
    'path_a_column', 'path_b_column', 'path_ab_column',
    'path_a_location', 'path_b_location', 'path_ab_location',
    'weight_a_N', 'weight_b_N', 'weight_ab_N', 'weight_sum_N',
    'remove_a_feasible', 'remove_b_feasible',
    'remove_both_tested', 'remove_both_feasible',
    'remove_a_error', 'remove_b_error', 'remove_both_error',
    'remove_a_worst_failure', 'remove_b_worst_failure',
    'remove_both_worst_failure',
    'selected_removed_columns', 'selection_reason', 'selected_error',
]
with open(f'{PREFIX}/cluster_additive_path_decisions.tsv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=decision_fields, delimiter='\t')
    writer.writeheader()
    writer.writerows(path_conflict_decisions)

remaining_path_conflicts = find_path_conflicts(
    path_conflict_candidate_columns(current_solution),
    nclose_path_usage,
    current_solution['fullsize'],
)
if remaining_path_conflicts:
    remaining_types = sorted({
        conflict.conflict_type for conflict in remaining_path_conflicts
    })
    raise RuntimeError(
        "Stage-24 compositional path filtering ended with unresolved "
        f"conflicts: {remaining_types}"
    )

logging.info(
    f'Compositional path filtering end: decisions={len(path_conflict_decisions)}, '
    f'removed_columns={len(path_conflict_removed_columns)}, '
    f'remaining_columns={len(current_solution["idx_list"])}'
)

cluster_column_list = list(current_solution['idx_list'])
final_idx_array_jl = jl.Vector[jl.Int]([i + 1 for i in cluster_column_list])
A_final_jl = jl.view(A_jl, jl.Colon(), final_idx_array_jl)
A_fail_final_jl = jl.view(A_fail_jl, jl.Colon(), final_idx_array_jl)
final_weight = current_solution['weight'].copy()
final_weights_fullsize = np.zeros(len(weight_base), dtype=weight_base.dtype)

cluster_weight_before_min = final_weight.copy()
cluster_min_mask = (
    (cluster_weight_before_min > 0)
    & (cluster_weight_before_min <= CLUSTER_START_DEPTH * N)
)
final_weight[cluster_weight_before_min <= 0] = 0
final_weight[cluster_min_mask] = 0
final_weight_jl = make_jl_weight(final_weight)

for i, v in enumerate(final_weight):
    final_weights_fullsize[cluster_column_list[i]] = v

predict_B_succ = np.asarray(A_final_jl * final_weight_jl)
predict_B_fail = np.asarray(A_fail_final_jl * final_weight_jl)

b_norm = np.linalg.norm(B_depth)

predict_B_succ_depth = predict_B_succ[depth_success_slice]
error = np.linalg.norm(predict_B_succ_depth - B_depth)
predict_B = np.concatenate((predict_B_succ_depth, predict_B_fail))

logging.info(f'Cluster error : {error:.4f}')
logging.info(f'Cluster relative error : {error/b_norm:.4f}')

np.save(f'{PREFIX}/weight_cluster.npy', final_weights_fullsize)
np.save(f'{PREFIX}/predict_B_cluster.npy', predict_B)

# A true NNLS zero remains a surviving 0/PASS result.  Only a positive
# coefficient removed by the explicit 0.1N floor is considered filtered.
cluster_thresholded_columns = {
    cluster_column_list[i]
    for i, is_thresholded in enumerate(cluster_min_mask)
    if is_thresholded
}
cluster_positive_columns = {
    cluster_column_list[i]
    for i, weight in enumerate(cluster_weight_before_min)
    if weight > CLUSTER_START_DEPTH * N
}
cluster_active_columns = set(cluster_column_list) - cluster_thresholded_columns
filter_active_columns = set(
    nclose_filter_status['stages']['filter']['active_columns']
)
selected_columns = set(cluster_column_list)
cluster_direct_reasons = {}
for event_key in nclose_filter_status['event_keys']:
    if event_key in nclose_filter_status['stages']['filter'].get('reasons', {}):
        continue
    if event_key in raw_zero_applied_event_keys:
        cluster_direct_reasons[event_key] = 'FILTERED_24_RAW_READ_ZERO'
        continue
    carriers = carrier_columns(nclose_path_usage, event_key) & filter_active_columns
    selected_carriers = carriers & selected_columns
    if not selected_carriers:
        cluster_direct_reasons[event_key] = 'FILTERED_24_CLUSTER_SELECTION'
    elif not (selected_carriers & cluster_positive_columns) and (
        selected_carriers & cluster_thresholded_columns
    ):
        cluster_direct_reasons[event_key] = 'FILTERED_24_CLUSTER_MIN_WEIGHT'

record_filter_stage(
    nclose_filter_status,
    'cluster',
    'filter',
    nclose_path_usage,
    cluster_active_columns,
    cluster_direct_reasons,
    'FILTERED_24_CLUSTER_COFILTERED_PATH',
    forced_reasons=cluster_direct_reasons,
)
save_filter_status(PREFIX, nclose_filter_status)
