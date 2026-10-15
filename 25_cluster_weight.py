"""Karyotype-mode filter and cluster fits (``_filter`` / ``_cluster``).

Ported from the former stage 24. The reason codes keep their FILTERED_24_*
names so existing NClose provenance stays comparable.
"""
import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *
from karyotype_mode import (
    MODE_KARYOTYPE,
    MODE_VARIANT,
    decide_karyotype_mode,
    describe_karyotype_mode,
    remove_suffixed_outputs,
    save_karyotype_mode,
)
from structure_nclose import load_or_build_structure_model
from nclose_tracking import (
    calculate_event_weights,
    carrier_columns,
    load_filter_status,
    load_path_usage,
    save_filter_status,
)

from collections import Counter, defaultdict

import re
import ast
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

NODE_NAME = 1
CHR_CHANGE_IDX = 2
DIR_CHANGE_IDX = 3

DIR_FOR = 1
DIR_BAK = 1

BND_TYPE = 0
CTG_IN_TYPE = 1
TEL_TYPE = 2

CLUSTER_START_DEPTH = 0.1
CLUSTER_FINAL_DEPTH = 0.2

JOIN_SAME_CHR_BASELINE = 0.8
JOIN_DIFF_CHR_BASELINE = 0.9

TELOMERE_EXPANSION = 5 * K
KARYOTYPE_SECTION_MINIMUM_LENGTH = 100*K

TYPE4_CLUSTER_SIZE = 1 * M

BASE_ACCSUMABSMAX_RATIO = 1.5
CHROM_ERROR_FAIL_RATE = 2.0

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

def import_index_path(file_path : str) -> list:
    file_path_list = file_path.split('/')
    key = file_path_list[-2]
    cnt = int(file_path_list[-1].split('.')[0]) - 1

    return path_list_dict[key][cnt][0]

def import_index_cnt(file_path : str) -> list:
    file_path_list = file_path.split('/')
    key = file_path_list[-2]
    cnt = int(file_path_list[-1].split('.')[0]) - 1

    return path_list_dict[key][cnt][1]

def get_final_paf_name_from_index(index_file_path):
    final_paf_list = index_file_path.split('/')
    cnt = final_paf_list[-1].split('.')[0]
    final_paf_list[-3] = '20_depth'
    final_paf_path = '/'.join(final_paf_list[:-1]) + f'/{cnt}.paf'

    return final_paf_path

def import_bed(bed_path : str) -> dict:
    bed_data = defaultdict(list)
    with open(bed_path, "r") as f:
        for curr_data in f:
            curr_data = curr_data.split("\t")
            bed_data[curr_data[0]].append((int(curr_data[1]), int(curr_data[2])))
    return bed_data

def build_tar_chr_data() -> dict:
    """Pick one reference-like telomere-to-telomere path per chromosome.

    Same selection as the former stage-22 ``tar_chr_data.pkl``: prefer the path
    between the outermost telomere nodes without any chromosome or direction
    change, otherwise the change-free path covering most of the chromosome.
    chrY is taken from a chrX path when such a path carries chrY sequence.
    """
    telo_connected_node_dict = defaultdict(list)
    for chr_info, contig_id in extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH):
        telo_connected_node_dict[chr_info].append(contig_id)

    telo_chr_list = sorted(set(n[:-1] for n in telo_connected_node_dict.keys()), key=chr2int)

    default_path_list = []
    for chr_name in telo_chr_list:
        if not (telo_connected_node_dict[chr_name + 'f'] and telo_connected_node_dict[chr_name + 'b']):
            logging.warning(f"Pure chromosome path skipped: {chr_name} lacks a telomere end")
            continue
        chr_def_paths = []
        for tar_chr_name, ctg_idx, func in [(chr_name + 'f', CHR_STR, min), (chr_name + 'b', CHR_END, max)]:
            tar_contig_data_list = [(i, ppc_data[i]) for i in telo_connected_node_dict[tar_chr_name]]
            tar_contig_data = func(tar_contig_data_list, key=lambda t: t[1][ctg_idx])
            chr_def_paths.append((tar_chr_name, tar_contig_data[0]))

        default_path_list.append(chr_def_paths)

    tar_chr_data = dict()
    for (chrf, chrfid), (chrb, chrbid) in default_path_list:
        tar_chr = chrf[:-1]
        key = f'{chrf}_{chrb}'
        n = len(path_list_dict.get(key, []))

        # chrY large indel
        if tar_chr == 'chrX':
            non_bnd_path_list = []
            for i in range(n):
                path_loc = f"{PREFIX}/00_raw/{key}/{i + 1}.index.txt"
                idx_data = import_index_path(path_loc)
                path_counter_list = import_index_cnt(path_loc)

                for ch, v in path_counter_list:
                    if ch == 'chrY':
                        non_bnd_path_list.append((path_loc, v, idx_data[-1][2]))

            if len(non_bnd_path_list) > 0:
                tar_path_loc = max(non_bnd_path_list, key=lambda t: (-t[2], t[1]))[0]
                tar_chr_data['chrY'] = get_relative_path(get_final_paf_name_from_index(tar_path_loc))

        tar_path_loc = None

        for i in range(n):
            path_loc = f"{PREFIX}/00_raw/{key}/{i + 1}.index.txt"
            idx_data = import_index_path(path_loc)
            if {idx_data[1][1], idx_data[-2][1]} == {chrfid, chrbid} and (idx_data[-1][1], idx_data[-1][2]) == (0, 0):
                tar_path_loc = path_loc
                break

        if tar_path_loc is None:
            non_bnd_path_list = []
            for i in range(n):
                path_loc = f"{PREFIX}/00_raw/{key}/{i + 1}.index.txt"
                idx_data = import_index_path(path_loc)
                if (idx_data[-1][1], idx_data[-1][2]) != (0, 0):
                    break

                path_counter_list = import_index_cnt(path_loc)

                for ch, v in path_counter_list:
                    if ch == tar_chr:
                        non_bnd_path_list.append((path_loc, v))

            if not non_bnd_path_list:
                logging.warning(f"Pure chromosome path skipped: no change-free {key} path")
                continue
            tar_path_loc = max(non_bnd_path_list, key=lambda t: t[1])[0]

        tar_chr_data[tar_chr] = get_relative_path(get_final_paf_name_from_index(tar_path_loc))

    return tar_chr_data

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


def parse_chromosome_labels(s):
    """
    Parse '...<f|b>_...<f|b>' into a canonical tuple:
      (left_label, left_is_f, right_label, right_is_f)

    Rules:
    - The two ends must end with 'f' or 'b' (assert if not).
    - Keep each '...' label string intact (e.g., 'chr12', 'scaf_007', etc.).
    - Canonicalize by sorting so the lexicographically smaller label comes first.
      If labels are equal, put 'f' (True) before 'b' (False).
    - When swapping due to sorting, directions stay attached to their original labels.
      (So 'chr12f_chr1b' becomes ('chr1', True, 'chr12', False).)
    """
    m = re.fullmatch(r'(.+?)([fb])_(.+?)([fb])', s)
    assert m is not None, "Input must match ...(f|b)_...(f|b) pattern"

    a_label, a_dir_ch, b_label, b_dir_ch = m.groups()
    a_is_f = (a_dir_ch == 'f')
    b_is_f = (b_dir_ch == 'f')

    # Canonical order by label; if same label, 'f' (True) first.
    if (a_label > b_label) or (a_label == b_label and not a_is_f and b_is_f):
        # Swap ends to enforce canonical order; keep directions with their labels.
        a_label, b_label = b_label, a_label
        a_is_f, b_is_f = b_is_f, a_is_f

    return (a_label, a_is_f, b_label, b_is_f)


def max_aligned_match_length(
    seq_a: list[tuple[tuple[str, str], int]],
    seq_b: list[tuple[tuple[str, str], int]],
) -> int:
    """
    Return the maximum total matched length after sliding two piecewise-constant
    label sequences along one axis. A and B are lists of ((chrom, strand), length).
    Only regions with exactly the same (chrom, strand) contribute to the score.

    Algorithm:
      1) Convert each sequence into absolute intervals [(start, end, label)].
      2) Consider candidate shifts = {a_ep - b_ep | a_ep in endpoints(A), b_ep in endpoints(B)}.
         (The overlap configuration only changes when an endpoint meets another.)
      3) For each shift, line-sweep over the two interval lists and accumulate
         overlap length where labels are equal.
      4) Return the maximum accumulated length across all shifts.

    Time complexity:
      Let n, m be #segments. Endpoints ~ (n+1), (m+1).
      Candidates O((n+1)*(m+1)); each evaluation O(n+m). Works well for tens~hundreds of segments.
    """
    # --- build absolute intervals: [(start, end, label)] and endpoint lists ---
    def build_intervals(seq):
        intervals = []
        endpoints = []
        pos = 0
        endpoints.append(pos)
        for (label, length) in seq:
            start = pos
            end = pos + length
            intervals.append((start, end, label))
            pos = end
            endpoints.append(pos)
        return intervals, endpoints

    A, A_ep = build_intervals(seq_a)
    B, B_ep = build_intervals(seq_b)

    if not A or not B:
        return 0

    # --- generate candidate shifts (all endpoint differences) ---
    # shift d means: compare A intervals with B intervals shifted by +d
    candidates = set()
    for a_e in A_ep:
        for b_e in B_ep:
            candidates.add(a_e - b_e)

    # --- overlap length for a given shift ---
    def match_length_for_shift(d: int) -> int:
        i, j = 0, 0
        total = 0
        # Two-pointer sweep over A and shifted-B
        while i < len(A) and j < len(B):
            a_s, a_e, a_lab = A[i]
            b_s, b_e, b_lab = B[j]
            b_s += d
            b_e += d

            # If no overlap, advance the one that ends earlier / starts later
            if a_e <= b_s:
                i += 1
                continue
            if b_e <= a_s:
                j += 1
                continue

            # Overlapping segment
            ov_s = a_s if a_s > b_s else b_s
            ov_e = a_e if a_e < b_e else b_e
            if ov_e > ov_s and a_lab == b_lab:
                total += (ov_e - ov_s)

            # Advance the interval that ends first
            if a_e <= b_e:
                i += 1
            else:
                j += 1
        return total

    best = 0
    # (Optional) small heuristic: iterate over sorted candidates for deterministic behavior
    for d in sorted(candidates):
        val = match_length_for_shift(d)
        if val > best:
            best = val

    return best

def should_join_by_baseline(
    seq_a: list[tuple[tuple[str, str], int]],
    seq_b: list[tuple[tuple[str, str], int]]
) -> bool:
    """
    Decide if two sequences should be joined based on:
      max_aligned_match_length(seq_a, seq_b) / max(total_len_a, total_len_b) >= JOIN_BASELINE

    Notes:
      - Returns False if both sequences have total length 0 (to avoid 0-division).
      - Assumes non-negative lengths.
      - Threshold is inclusive (>=).

    """

    total_a = sum(length for (_, length) in seq_a)
    total_b = sum(length for (_, length) in seq_b)

    denom = total_a if total_a >= total_b else total_b
    if denom == 0:
        return False

    seq_a_chr = (seq_a[0][0], seq_a[-1][0])
    seq_b_chr = (seq_b[0][0], seq_b[-1][0])

    score = max_aligned_match_length(seq_a, seq_b)
    return (score / denom) >= (JOIN_SAME_CHR_BASELINE if seq_a_chr == seq_b_chr else JOIN_DIFF_CHR_BASELINE)

def root_find(uf, path):
    if path not in uf:
        uf[path] = path
    elif uf[path] != path:
        uf[path] = root_find(uf, uf[path])
    return uf[path]

def get_karyotype_summary_relpath(non_type4_path_list : list):
    karyotypes_data_direction_include = {}
    karyotypes_nclose_count = {}
    
    for path_path in non_type4_path_list:
        pieces = []
        path = import_index_path(path_path)

        # padding for easier calculation
        if len(path[0]) < 4:
            path[0] = tuple([0] + list(path[0]))
        if len(path[-1]) < 4:
            path[-1] = tuple([0] + list(path[-1]))
        
        curr_ref = 0 if path[0][NODE_NAME][-1] =='f' else chr_len[path[0][NODE_NAME][:-1]]
        curr_incr = '+' if path[0][NODE_NAME][-1] =='f' else '-'
        curr_chr = [path[0][NODE_NAME][:-1], curr_incr]
        
        nclose_use_cnt = 0
        for i in range(1, len(path)-1):
            if path[i][CHR_CHANGE_IDX] > path[i-1][CHR_CHANGE_IDX] \
            or path[i][DIR_CHANGE_IDX] > path[i-1][DIR_CHANGE_IDX]:
                nclose_use_cnt += 1
                last_node = ppc_data[path[i-1][NODE_NAME]]
                curr_node = ppc_data[path[i][NODE_NAME]]
                # add last piece
                if curr_incr == '+':
                    pieces.append((tuple(curr_chr), last_node[CHR_END] - curr_ref))
                else:
                    pieces.append((tuple(curr_chr), curr_ref - last_node[CHR_STR]))
                
                # update info of new piece (starting ref, chromosome type, increment ..)
                if path[i][NODE_NAME] > path[i-1][NODE_NAME]:
                    curr_incr = curr_node[CTG_DIR]
                    curr_chr = [curr_node[CHR_NAM], curr_incr]
                    curr_ref = curr_node[CHR_STR] if curr_incr == '+' else curr_node[CHR_END]
                else:
                    curr_incr = '-' if curr_node[CTG_DIR] == '+' else '+'
                    curr_chr = [curr_node[CHR_NAM], curr_incr]
                    curr_ref = curr_node[CHR_STR] if curr_incr == '+' else curr_node[CHR_END]
        karyotypes_nclose_count[path_path] = nclose_use_cnt
        pieces.append((tuple(curr_chr), chr_len[curr_chr[0]] - curr_ref if curr_incr == '+' else curr_ref))
        karyotypes_data_direction_include[path_path] = pieces

    return karyotypes_data_direction_include

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

parser.add_argument("censat_bed_path",
                    help="Path to the censat repeat information file.")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("-t", "--thread", help="Number of threads", type=int)

parser.add_argument("--mode", choices=(MODE_KARYOTYPE, MODE_VARIANT),
                    default=MODE_KARYOTYPE,
                    help="Requested output mode. Karyotype mode is forced to variant "
                         "for VCF input or graphs above the NClose node limit.")

args = parser.parse_args()

PREFIX = args.prefix
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
main_stat_loc = args.main_stat_loc
TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
TELO_CONNECT_NODES_INFO_PATH = PREFIX+"/telomere_connected_list.txt"
front_contig_path = RATIO_OUTLIER_FOLDER+"front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER+"back_jump/"
ecdna_contig_path = RATIO_OUTLIER_FOLDER+"ecdna/"
output_folder = f'{PREFIX}/21_pat_depth'
NCLOSE_FILE_PATH = f"{args.prefix}/nclose_nodes_index.txt"

karyotype_mode_config = decide_karyotype_mode(
    args.mode,
    2 * sum(len(pairs) for pairs in load_nclose_nodes(PREFIX).values()),
    pipeline_input_is_vcf(load_pipeline_input(PREFIX)),
)
save_karyotype_mode(PREFIX, karyotype_mode_config)
logging.info(describe_karyotype_mode(karyotype_mode_config))

if not karyotype_mode_config["karyotype_mode"]:
    for removed_path in remove_suffixed_outputs(PREFIX):
        logging.info(f"Removed stale karyotype-mode output: {removed_path}")
    stale_status = load_filter_status(PREFIX)
    if {'filter', 'cluster'} & set(stale_status['stages']):
        stale_status['stages'].pop('filter', None)
        stale_status['stages'].pop('cluster', None)
        save_filter_status(PREFIX, stale_status)
    sys.exit(0)
logging.info("25_cluster_weight start")

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

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)

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
    matrix_depth_meta = pkl.load(f)
require_fixed_depth_matrix(matrix_depth_meta)
chr_filt_st_list = matrix_depth_meta["chr_filt_st_list"]
bed_data = import_bed(args.censat_bed_path)

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

def build_provenance_usage():
    """Path usage plus native structure constituents of every matrix column.

    Stage 31 counts the constituent NCloses of merged Type-4 and ecDNA columns,
    so filter/cluster provenance treats those columns as carriers as well.
    The greedy filter itself still removes columns by path usage only.
    """
    usage = [Counter(column_usage) for column_usage in nclose_path_usage]
    for structure in load_or_build_structure_model(PREFIX)["structures"]:
        column = structure["feature_index"]
        if column is None:
            continue
        counts = structure.get("source_nclose_counts", structure["nclose_counts"])
        for event_key, count in counts.items():
            usage[column][event_key] = max(usage[column][event_key], int(count))
    return usage

provenance_usage = build_provenance_usage()
nclose_filter_status = load_filter_status(PREFIX)
if 'base' not in nclose_filter_status['stages']:
    raise ValueError(
        "Missing stage-23 NClose filter provenance. Rerun SKYPE from stage 23."
    )
base_active_columns = set(nclose_filter_status['stages']['base']['active_columns'])
# Stage 31 also reports restored constituents outside the stage-23 catalog.
tracked_keys = list(dict.fromkeys(
    [*nclose_filter_status['event_keys'], *(key for usage in provenance_usage for key in usage)]
))


def record_structure_stage(stage, active_columns, reason_for):
    """Record each source's status in one fit from its native-structure carriers.

    A source passes while any carrier column stays in the fit, even if an
    earlier stage removed its own paths, matching the stage-31 accounting.
    """
    reasons = {}
    for key in tracked_keys:
        reason = reason_for(key, carrier_columns(provenance_usage, key))
        if reason != 'PASS':
            reasons[key] = reason
    nclose_filter_status['stages'][stage] = {
        'active_columns': sorted(active_columns), 'reasons': reasons,
    }
    save_filter_status(PREFIX, nclose_filter_status)


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

def get_fit_failures(predict_suc_B, fail_rate=CHROM_ERROR_FAIL_RATE):
    failures = []
    for chrom, stat in chrom_acc_sum_peak_stats(predict_suc_B).items():
        if chrom == 'chrY' and no_chrY:
            continue
        acc_sum_max_base = chrom_acc_sum_dict_max_base[chrom]
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
# Stage 23 fits every matrix column; explicitly excluded columns are zero
# columns outside the base provenance and must stay removed here.
base_A_idx_set = set(A_idx_list) & base_active_columns

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

filter_active_column_set = set(final_idx_list)
base_reasons = nclose_filter_status['stages']['base'].get('reasons', {})


def filter_reason(key, carriers):
    if carriers & filter_active_column_set:
        return 'PASS'
    if not carriers & base_active_columns:
        return base_reasons.get(key, 'FILTERED_02_COFILTERED_PATH')
    return 'FILTERED_24_GREEDY_DEPTH' if key in not_essential_nclose else 'FILTERED_24_COFILTERED_PATH'


nclose_filter_status['stages'].pop('cluster', None)
record_structure_stage('filter', filter_active_column_set, filter_reason)

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


# Stage 22 records the exact matrix-column order.
with open(f'{PREFIX}/tot_loc_list.pkl', 'rb') as f:
    tot_loc_list = pkl.load(f)

if len(tot_loc_list) != len(final_weights_fullsize):
    raise ValueError(
        "Stage-25 matrix-column metadata mismatch: "
        f"{len(tot_loc_list)} locations for {len(final_weights_fullsize)} weights."
    )
loc2weight = dict(zip(tot_loc_list, final_weights_fullsize))
filter_active_locations = {
    tot_loc_list[col_idx]
    for col_idx in filter_active_column_set
}

cluster_tar_path_list = []

# Add pure chromosome
tar_chr_data = build_tar_chr_data()

loc_prefix = paf_ans_list[0][0].split('/')[:-3]

for path_tuple in tar_chr_data.values():
    path = '/'.join(loc_prefix + list(path_tuple))
    if path in filter_active_locations:
        cluster_tar_path_list.append(path)

for ind, w in weights_sorted_data:
    paf_loc = tot_loc_list[ind]
    if paf_loc.split('/')[-3] not in {'11_ref_ratio_outliers', '12_cent_fragment'}:
        if w > CLUSTER_START_DEPTH * N:
            cluster_tar_path_list.append(paf_loc)

tot_loc_list2nclosecnt = dict()

for paf_loc in tot_loc_list:
    if paf_loc.split('/')[-3] not in {'11_ref_ratio_outliers', '12_cent_fragment'}:
        path = import_index_path(paf_loc)

        if len(path[0]) < 4:
            path[0] = tuple([0] + list(path[0])) # padding for easier calculation
        if len(path[-1]) < 4:
            path[-1] = tuple([0] + list(path[-1]))

        nclose_use_cnt = 0
        for i in range(1, len(path)-1):
            if path[i][CHR_CHANGE_IDX] > path[i-1][CHR_CHANGE_IDX] \
            or path[i][DIR_CHANGE_IDX] > path[i-1][DIR_CHANGE_IDX]:
                nclose_use_cnt += 1

        tot_loc_list2nclosecnt[paf_loc] = nclose_use_cnt
        
# karyotypes_data : key (rel_path) => value (karyotype value)
karyotypes_data = get_karyotype_summary_relpath(cluster_tar_path_list)

chr_set_merge = defaultdict(list)
rel_path2ncnt = dict((v, i) for i, v in enumerate(tot_loc_list))

for path, seq_list in karyotypes_data.items():
    matched_master = None

    # Compare only against master paths
    for master_path, slave_list in chr_set_merge.items():
        master_seq = karyotypes_data[master_path]
        if should_join_by_baseline(master_seq, seq_list):
            matched_master = master_path
            break

    if matched_master is not None:
        # Join to existing cluster as slave
        chr_set_merge[matched_master].append(path)
    else:
        # Create new cluster with this path as master
        chr_set_merge[path] = [path]

using_merge_ncnt_list = []
for path_list in chr_set_merge.values():
    merge_weight = sum(loc2weight[p] for p in path_list)
    if merge_weight > CLUSTER_FINAL_DEPTH * N:
        repr_path = min(path_list, key=lambda t : tot_loc_list2nclosecnt[t])
        using_merge_ncnt_list.append(rel_path2ncnt[repr_path])

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
                    using_merge_ncnt_list.append(ncnt)
        elif event_type == 'ecdna':
            using_merge_ncnt_list.append(ncnt)
        else:
            assert(False)
    elif paf_loc.split('/')[-3] == '12_cent_fragment':
        using_merge_ncnt_list.append(ncnt)

using_merge_ncnt_list = sorted(
    set(using_merge_ncnt_list) & filter_active_column_set
)
if not using_merge_ncnt_list:
    raise RuntimeError("No NNLS columns remain for stage-25 cluster selection.")

# Julia run partial NNLS by using_merge_ncnt_list

final_idx_array_jl = jl.Vector[jl.Int]([i + 1 for i in using_merge_ncnt_list])
A_final_jl = jl.view(A_jl, jl.Colon(), final_idx_array_jl)
A_fail_final_jl = jl.view(A_fail_jl, jl.Colon(), final_idx_array_jl)


solve_gram_mem_gb = nnls_gram_mem_gb()
final_weight_jl = jl.nnls_solve(
    A_final_jl, B_jl, THREAD, False,
    gram_mem_gb=solve_gram_mem_gb,
    w_init=make_jl_weight(final_weights_fullsize[using_merge_ncnt_list])
)
final_weight = np.asarray(final_weight_jl)
final_weights_fullsize = np.zeros(jl.size(A_jl, 2))

cluster_weight_before_min = final_weight.copy()
cluster_min_mask = (
    (cluster_weight_before_min > 0)
    & (cluster_weight_before_min <= CLUSTER_START_DEPTH * N)
)
final_weight[cluster_weight_before_min <= 0] = 0
final_weight[cluster_min_mask] = 0
final_weight_jl = jl.Vector[jl.eltype(final_weight_jl)](final_weight)

for i, v in enumerate(final_weight):
    final_weights_fullsize[using_merge_ncnt_list[i]] = v

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
    using_merge_ncnt_list[i]
    for i, is_thresholded in enumerate(cluster_min_mask)
    if is_thresholded
}
cluster_positive_columns = {
    using_merge_ncnt_list[i]
    for i, weight in enumerate(cluster_weight_before_min)
    if weight > CLUSTER_START_DEPTH * N
}
selected_columns = set(using_merge_ncnt_list)
cluster_active_columns = selected_columns - cluster_thresholded_columns
filter_reasons = nclose_filter_status['stages']['filter']['reasons']


def cluster_reason(key, carriers):
    # The floor cut wins even if an exact-zero carrier stays selected.
    selected = carriers & selected_columns
    if (selected & cluster_thresholded_columns) and not (selected & cluster_positive_columns):
        return 'FILTERED_24_CLUSTER_MIN_WEIGHT'
    if carriers & cluster_active_columns:
        return 'PASS'
    return filter_reasons.get(key, 'FILTERED_24_CLUSTER_SELECTION')


record_structure_stage('cluster', cluster_active_columns, cluster_reason)
