import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *

from collections import defaultdict
from juliacall import Main as jl

import ast
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
CTG_RPTCASE = 17
CTG_MAINFLOWDIR = 18
CTG_MAINFLOWCHR = 19

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


TYPE4_CLUSTER_SIZE = 10 * M

def get_relative_path(p):
    return tuple(p.split('/')[-3:])

def extract_nclose_node(nclose_path: str) -> list:
    nclose_list = []
    with open(nclose_path, "r") as f:
        for line in f:
            line = line.split()
            nclose_list.append((int(line[1]), int(line[2])))
    return nclose_list

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

def import_telo_data(file_path : str, chr_len : dict) -> dict :
    fai_file = open(file_path, "r")
    telo_data = [(0,1)]
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        int_induce_idx = [1, 2]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        if temp_list[0]!=telo_data[-1][0]:
            temp_list[2]+=TELOMERE_EXPANSION
            temp_list.append('f')
        else:
            if temp_list[1]>chr_len[temp_list[0]]/2:
                temp_list[1]-=TELOMERE_EXPANSION
                temp_list.append('b')
            else:
                temp_list.append('f')
                temp_list[2]+=TELOMERE_EXPANSION
        telo_data.append(tuple(temp_list))
    fai_file.close()
    return telo_data[1:]

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

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)

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
# 24_cluster_weight.py /home/hyunwoo/ACCtools-pipeline/90_skype_run/Caki-1/20_alignasm/Caki-1.ctg.aln.paf.ppc.paf /home/hyunwoo/ACCtools-pipeline/90_skype_run/Caki-1/01_depth/Caki-1_normalized.win.stat.gz public_data/chm13v2.0_telomere.bed public_data/chm13v2.0.fa.fai 30_skype_pipe/Caki-1_18_40_59 -t 128
# """
# args = parser.parse_args(t.strip().split()[1:])


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

with open(f"{PREFIX}/report.txt", 'r') as f:
    f.readline()
    path_cnt = int(f.readline().strip())

use_julia_solver = path_cnt <= HARD_PATH_COUNT_BASELINE

if not use_julia_solver:
    exit(0)
logging.info("24_cluster_weight start")

THREAD = args.thread
core_num = psutil.cpu_count(logical=False)
if core_num is None:
    THREAD = THREAD
else:
    THREAD = min(int(THREAD), core_num)

nclose_nodes = set(extract_nclose_node(NCLOSE_FILE_PATH))
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
weights = np.load(f'{PREFIX}/weight.npy')
weights_sorted_data = sorted(enumerate(weights), key=lambda t:t[1], reverse=True)


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

loc2weight = dict(zip(tot_loc_list, weights))

cluster_tar_path_list = []

# Add pure chromosome
with open(f"{PREFIX}/tar_chr_data.pkl", "rb") as f:
    tar_chr_data = pkl.load(f)

loc_prefix = paf_ans_list[0][0].split('/')[:-3]

for path_tuple in tar_chr_data.values():
    path = '/'.join(loc_prefix + list(path_tuple))
    cluster_tar_path_list.append(path)

for ind, w in weights_sorted_data:
    paf_loc = tot_loc_list[ind]
    if paf_loc.split('/')[-3] != '11_ref_ratio_outliers':
        if w > CLUSTER_START_DEPTH * N:
            cluster_tar_path_list.append(paf_loc)

tot_loc_list2nclosecnt = dict()

for paf_loc in tot_loc_list:
    if paf_loc.split('/')[-3] != '11_ref_ratio_outliers':
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
        with open(paf_loc, "r") as f:
            l = f.readline()
            l = l.rstrip()
            l = l.split("\t")
            chr_nam1 = l[CHR_NAM]
            chr_nam2 = l[CHR_NAM]
            pos1 = int(l[CHR_STR])
            pos2 = int(l[CHR_END])
            if abs(pos1-pos2) > TYPE4_CLUSTER_SIZE:
                using_merge_ncnt_list.append(ncnt)

using_merge_ncnt_list.sort()
using_merge_ncnt_arr = np.asarray(using_merge_ncnt_list)

jl.include(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'anderson_nnls.jl'))
jl.seval("using HDF5, LinearAlgebra")

A_jl, B_jl, b_start_ind = jl.load_nnls_array(f'{PREFIX}/matrix.h5')
A_fail_jl = jl.load_fail_array(f'{PREFIX}/matrix.h5')
B = np.asarray(B_jl)

final_idx_array_jl = jl.Vector[jl.Int]([i + 1 for i in using_merge_ncnt_list])
A_final_jl = jl.view(A_jl, jl.Colon(), final_idx_array_jl)
A_fail_final_jl = jl.view(A_fail_jl, jl.Colon(), final_idx_array_jl)


final_weight_jl = jl.nnls_solve(A_final_jl, B_jl, THREAD, False)
final_weight = np.asarray(final_weight_jl)
final_weights_fullsize = np.zeros(jl.size(A_jl, 2))

final_weight[final_weight <= CLUSTER_START_DEPTH * N] = 0
final_weight_jl = jl.Vector[jl.eltype(final_weight_jl)](final_weight)

for i, v in enumerate(final_weight):
    final_weights_fullsize[using_merge_ncnt_list[i]] = v

predict_B_succ = np.asarray(A_final_jl * final_weight_jl)
predict_B_fail = np.asarray(A_fail_final_jl * final_weight_jl)

b_norm = np.linalg.norm(B)

error = np.linalg.norm(predict_B_succ - B)
predict_B = np.concatenate((predict_B_succ, predict_B_fail))[b_start_ind:]

logging.info(f'Cluster error : {error:.4f}')
logging.info(f'Cluster relative error : {error/b_norm:.4f}')

np.save(f'{PREFIX}/weight_cluster.npy', final_weights_fullsize)
np.save(f'{PREFIX}/predict_B_cluster.npy', predict_B)
