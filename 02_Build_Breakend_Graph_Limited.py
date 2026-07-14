import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *
from parse_vcf import parse_vcf_events, select_vcf_type4_graph_events

import shutil
import argparse
import subprocess
import json

import pickle as pkl
import pandas as pd
import numpy as np

import ast

import copy
import bisect

import graph_tool.all as gt
import graph_tool
import networkx as nx

from collections import defaultdict, Counter

from multiprocessing import Pool
from tqdm import tqdm
from functools import partial
import logging

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).

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
INF = 1000000000
BUFFER = 10000000

CHUKJI_LIMIT = -1
BND_CONTIG_BOUND = 0.1
TYPE2_CONTIG_MINIMUM_LENGTH = 200*K
SUBTELOMERE_REPEAT_LENGTH = 0 # 100*K
TYPE2_FLANKING_LENGTH = 500 * K
TYPE2_SIM_COMPARE_RAITO = 1.3
TYPE34_BREAK_CHUKJI_LIMIT = 1*M
CHUKJI_FAIL_TYPE2_RESCUE_THRESHOLD = 2*K
CIRCUIT_ECDNA_LENGTH_LIMIT = 80*M

NCLOSE_SIM_COMPARE_RAITO = 1.2
NCLOSE_SIM_DIFF_THRESHOLD = 5

TOT_PATH_LIMIT = 3*M
PAT_PATH_LIMIT = 10*K

DIR_CHANGE_LIMIT_ABS_MAX = 1
CENSAT_VISIT_LIMIT = 2

CHR_CHANGE_LIMIT_HARD_START = 2
HARD_NCLOSE_COUNT = PIPELINE_MODE_NCLOSE_LIMIT
FAIL_NCLOSE_COUNT = 10000

BND_OVERUSE_CNT = 2
PATH_MAJOR_COMPONENT = 3
NCLOSE_COMPRESS_LIMIT = 100*K
NCLOSE_MERGE_LIMIT = 1*K
ALL_REPEAT_NCLOSE_COMPRESS_LIMIT = 500*K
SUBTELO_TIP_LIMIT = 500*K
OFFSET_DIR_GROUP_LIMIT = 100*K


PIPELINE_02_FILTER = True
ENABLE_CENSAT_NONCENSAT_OFFSET_DIR_FILTER = PIPELINE_02_FILTER

PATH_COMPRESS_LIMIT = 50*K
IGNORE_PATH_LIMIT = 50*K
NON_REPEAT_NOISE_RATIO=0.1

CONTIG_MINIMUM_SIZE = 100*K
BND_CONTIG_BOUND = 0.1
RPT_BND_CONTIG_BOUND = 0.2

MAPQ_BOUND = 60

TELOMERE_EXPANSION = 5 * K
TELOMERE_COMPRESS_RANGE = 100*K
CENSAT_COMPRESSABLE_THRESHOLD = 1000*K
VIRTUAL_TELOMERE_FLANK = 1 * K
VIRTUAL_TELOMERE_NODE_PREFIX = "virtual_telo_"

FORCE_TELOMERE_THRESHOLD = 10*K
TELOMERE_CLUSTER_THRESHOLD = 500*K
SUBTELOMERE_LENGTH = 500*K
MULTI_END_ALIGNMENT_WINDOW = 500*K

CENSAT_OUT_DIFF_RATIO = 0.30

MIN_FLANK_SIZE_BP = 1*M

WEAK_BREAKEND_CEN_RATIO_THRESHOLD = 1.3
BREAKEND_CEN_RATIO_THRESHOLD = 1.5

REPEAT_MERGE_GAP = 0

MAX_OVERLAP_SCORE = 3

SPLIT_CTG_LEN_LIMIT = 100 * K
TRUST_SPLIT_CTG_LEN_LIMIT = 20 * K

MIN_PATH_REF_LEN = 5 * M

TYPE2_CHUKJI_AS_TYPE4 = 5 * M
TYPE2_CONJOIN_COMPRESS_LIMIT = 1 * M
TYPE2_DIST_FLIP_THRESHOLD = 100 * K

TYPE4_INDEL_GRAPH_EDGE_PKL = 'type4_indel_graph_edges.pkl'
TYPE4_INDEL_GRAPH_MIN_SPAN = 1 * M
TYPE4_INDEL_GRAPH_DEPTH_WINDOW = 500 * K
TYPE4_INDEL_GRAPH_DEPTH_DIFF_RATIO = 0.2

RAW_TRANSLOCATION_CANDIDATE_PKL = 'raw_translocation_candidates.pkl'
RAW_TRANSLOCATION_RESULT_PKL = 'raw_translocation_result.pkl'
RAW_TRANSLOCATION_WINDOW = 5 * K
RAW_TRANSLOCATION_MIN_SAME_CHROM_SPAN = 10 * M

NCLOSE_COUNT_CANDIDATE_PKL = 'nclose_count_candidates.pkl'
NCLOSE_COUNT_RESULT_PKL = 'nclose_count_result.pkl'
NCLOSE_COUNT_DEFAULT_VAF_THRESHOLD = 0.1

ALT_SIMPLE_MIN_SEGMENT_LEN = 10 * K
ALT_SIMPLE_MAJOR_CHR_RATIO = 0.90
ALT_SIMPLE_EXISTING_NCLOSE_DIST = 10 * K

JULIA_BAM_THREAD_LIM = 32

def import_data(file_path : str) -> list :
    contig_data = []
    int_induce_idx = [1, 2, 3, 6, 7, 8, 9]
    idx = 0
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            a = curr_contig.split("\t")
            temp_list = a[:9]
            temp_list.append(a[11])
            for i in int_induce_idx:
                temp_list[i] = int(temp_list[i])
            temp_list.append(idx)
            contig_data.append(temp_list)
            idx+=1
    return contig_data


def terminal_alignment_anchors(row, window=MULTI_END_ALIGNMENT_WINDOW) -> set:
    """Return chromosome ends that fully contain this PAF alignment."""
    chrom = row[CHR_NAM]
    chrom_len = int(row[CHR_LEN])
    ref_st = int(row[CHR_STR])
    ref_nd = int(row[CHR_END])
    anchors = set()

    if ref_nd <= min(int(window), chrom_len):
        anchors.add((chrom, 'f'))
    if ref_st >= max(0, chrom_len - int(window)):
        anchors.add((chrom, 'b'))
    return anchors


def find_multi_end_aligned_contigs(contig_data, window=MULTI_END_ALIGNMENT_WINDOW):
    """Find contigs whose every PAF row lies within at least two chromosome ends."""
    row_indices_by_contig = defaultdict(list)
    anchors_by_contig = defaultdict(set)
    nonterminal_contigs = set()

    for row_idx, row in enumerate(contig_data):
        contig_name = row[CTG_NAM]
        row_indices_by_contig[contig_name].append(row_idx)
        anchors = terminal_alignment_anchors(row, window)
        if anchors:
            anchors_by_contig[contig_name].update(anchors)
        else:
            nonterminal_contigs.add(contig_name)

    excluded_contigs = {
        contig_name
        for contig_name, anchors in anchors_by_contig.items()
        if contig_name not in nonterminal_contigs and len(anchors) >= 2
    }
    excluded_row_indices = {
        row_idx
        for contig_name in excluded_contigs
        for row_idx in row_indices_by_contig[contig_name]
    }
    return excluded_contigs, excluded_row_indices


def import_data2(file_path : str) -> list :
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


def import_repeat_data_00(file_path : str) -> dict :
    fai_file = open(file_path, "r")
    repeat_data = defaultdict(list)
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        curr_patch = (int(temp_list[1]), int(temp_list[2]))
        if len(repeat_data[temp_list[0]])==0:
            repeat_data[temp_list[0]].append(curr_patch)
        else:
            if abs(curr_patch[0] - repeat_data[temp_list[0]][-1][1]) < REPEAT_MERGE_GAP:
                repeat_data[temp_list[0]][-1] = (repeat_data[temp_list[0]][-1][0], curr_patch[1])
            else:
                repeat_data[temp_list[0]].append(curr_patch)
    # for k, v in repeat_data.items():
    #     print(k, len(v))
    fai_file.close()
    return repeat_data

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

def is_stepwise_nonoverlapping(intervals: list) -> bool:
    n = len(intervals)
    if n < 2:
        return True

    for i in range(1, n):
        st_i, nd_i       = intervals[i]
        st_prev, nd_prev = intervals[i-1]

        if st_i <= st_prev:
            return False

        if nd_i <= nd_prev:
            return False

        if i >= 2:
            _, nd_prev2 = intervals[i-2]
            if st_i < nd_prev2:
                return False

    return True

def is_alt_contained_in_segments(segments: list, alt_segments: list) -> bool:
    n = len(segments)
    for st_b, nd_b in alt_segments:
        ok = False
        for i, (st_a, nd_a) in enumerate(segments):
            cond_start = True if i == 0 else (st_a <= st_b)
            cond_end   = True if i == n-1 else (nd_b <= nd_a)

            if cond_start and cond_end:
                ok = True
                break
        if not ok:
            return False
    return True

def max_overlap(intervals, target_intervals):
    events = []
    temp_inf = 2 * len(intervals) if target_intervals else 0

    for l in target_intervals:
        st, nd = l[1], l[2]
        events.append((st, +temp_inf))
        events.append((nd, -temp_inf))

    for l in intervals:
        st, nd = l[1], l[2]
        events.append((st, +1))
        events.append((nd, -1))

    events.sort(key=lambda x: (x[0], x[1]))

    curr = max_cnt = 0
    for _, delta in events:
        curr += delta
        max_cnt = max(max_cnt, curr)
    
    return max_cnt - temp_inf

def get_qry_cord_data(paf_path: str, get_ori_cord: bool = False) -> tuple:
    end_idx_dict = dict()
    paf_data_list = []

    with open(paf_path, "r") as paf_file:
        for line in paf_file:
            cols = line.rstrip("\n").split("\t")
            # [contig, start, end, mapq]
            tar_data = [
                cols[0],
                int(cols[2]),
                int(cols[3]),
                int(cols[11])
            ]

            if get_ori_cord:
                tar_col = None
                for col in cols:
                    if col.startswith('xi:A:'):
                        tar_col = col
                        break

                ori_type, ori_cord = tar_col.split(':')[-1].split('_')
                tar_data.append(int(ori_cord) if ori_type == 'P' else -1)

            paf_data_list.append(tar_data)

    if paf_data_list:
        current_contig = paf_data_list[0][0]
        start_idx = 0
        for i, entry in enumerate(paf_data_list):
            contig = entry[0]
            if contig != current_contig:
                end_idx_dict[current_contig] = (start_idx, i - 1)
                current_contig = contig
                start_idx = i

        end_idx_dict[current_contig] = (start_idx, len(paf_data_list) - 1)
    
    return paf_data_list, end_idx_dict

def div_repeat_paf(original_paf_path_list: list, aln_paf_path_list: list, contig_data: list) -> set:
    not_using_contig = set()
    
    ori_paf_data_list_data = []
    ori_end_idx_dict_data = []

    for original_paf_path in original_paf_path_list:
        ori_paf_data_list, ori_end_idx_dict = get_qry_cord_data(original_paf_path)

        ori_paf_data_list_data.append(ori_paf_data_list)
        ori_end_idx_dict_data.append(ori_end_idx_dict)

    aln_paf_data_list_data = []
    aln_end_idx_dict_data = []

    for aln_paf_path in aln_paf_path_list:
        aln_paf_data_list, aln_end_idx_dict = get_qry_cord_data(aln_paf_path, get_ori_cord=True)

        aln_paf_data_list_data.append(aln_paf_data_list)
        aln_end_idx_dict_data.append(aln_end_idx_dict)

    ppc_aln_end_idx_dict = dict()

    s = 0
    contig_data_size = len(contig_data)
    while s < contig_data_size:
        e = contig_data[s][CTG_ENDND]
        ori_ctg_name = contig_data[s][CTG_NAM]
        
        s_cl_ind, s_c_ind = map(int, contig_data[s][CTG_GLOBALIDX].split('.'))
        e_cl_ind, e_c_ind = map(int, contig_data[e][CTG_GLOBALIDX].split('.'))

        assert(s_cl_ind == e_cl_ind)

        if s_cl_ind < 2:
            assert(ori_ctg_name_data[s_cl_ind][s_c_ind] == ori_ctg_name_data[e_cl_ind][e_c_ind])
            ppc_aln_end_idx_dict[ori_ctg_name] = (s_cl_ind, s_c_ind, e_c_ind)
        
        s = e+1

    for ctg_name, (cl_ind, ppc_st, ppc_nd) in ppc_aln_end_idx_dict.items():
        ori_paf_data_list = ori_paf_data_list_data[cl_ind]
        ori_end_idx_dict = ori_end_idx_dict_data[cl_ind]

        aln_paf_data_list = aln_paf_data_list_data[cl_ind]
        aln_end_idx_dict = aln_end_idx_dict_data[cl_ind]

        ori_ctg_name = ori_ctg_name_data[cl_ind][ppc_st]

        ori_st, ori_nd = ori_end_idx_dict[ori_ctg_name]
        aln_st, aln_nd = aln_end_idx_dict[ori_ctg_name]
        assert(aln_st <= ppc_st and ppc_nd <= aln_nd)

        overlap_score = max_overlap(ori_paf_data_list[ori_st:ori_nd+1], [aln_paf_data_list[ppc_st], aln_paf_data_list[ppc_nd]])

        is_strong_paf = False

        ori_aln_ind = []
        for aln_ind in range(ppc_st, ppc_nd + 1):
            ori_ind = aln_paf_data_list[aln_ind][-1]

            if ori_ind != -1:
                ori_aln_ind.append(ori_ind)

        ori_aln_ind_set = set(ori_aln_ind)
        
        ori_aln_intervals = []
        ori_alt_intervals = []
        for ori_ind in range(ori_st, ori_nd + 1):
            _, st, nd, _ = ori_paf_data_list[ori_ind]

            if ori_ind in ori_aln_ind_set:
                ori_aln_intervals.append((st, nd))
            else:
                ori_alt_intervals.append((st, nd))
        
        if is_stepwise_nonoverlapping(ori_aln_intervals):
            is_strong_paf = is_alt_contained_in_segments(ori_aln_intervals, ori_alt_intervals)

        if overlap_score >= MAX_OVERLAP_SCORE and not is_strong_paf:
            not_using_contig.add(ctg_name)
        
    return not_using_contig

def get_overlap_bothend_score_dict(original_paf_path_list: list, aln_paf_path_list: list, contig_data: list) -> dict:
    ori_paf_data_list_data = []
    ori_end_idx_dict_data = []

    for original_paf_path in original_paf_path_list:
        ori_paf_data_list, ori_end_idx_dict = get_qry_cord_data(original_paf_path)

        ori_paf_data_list_data.append(ori_paf_data_list)
        ori_end_idx_dict_data.append(ori_end_idx_dict)

    aln_paf_data_list_data = []
    aln_end_idx_dict_data = []

    for aln_paf_path in aln_paf_path_list:
        aln_paf_data_list, aln_end_idx_dict = get_qry_cord_data(aln_paf_path, get_ori_cord=True)

        aln_paf_data_list_data.append(aln_paf_data_list)
        aln_end_idx_dict_data.append(aln_end_idx_dict)

    ppc_aln_end_idx_dict = dict()

    s = 0
    contig_data_size = len(contig_data)
    while s < contig_data_size:
        e = contig_data[s][CTG_ENDND]
        ori_ctg_name = contig_data[s][CTG_NAM]
        
        s_cl_ind, s_c_ind = map(int, contig_data[s][CTG_GLOBALIDX].split('.'))
        e_cl_ind, e_c_ind = map(int, contig_data[e][CTG_GLOBALIDX].split('.'))

        assert(s_cl_ind == e_cl_ind)

        if s_cl_ind < 2:
            assert(ori_ctg_name_data[s_cl_ind][s_c_ind] == ori_ctg_name_data[e_cl_ind][e_c_ind])
            ppc_aln_end_idx_dict[ori_ctg_name] = (s_cl_ind, s_c_ind, e_c_ind)
        s = e+1

    ctgname2overlap = dict()
    for ctg_name, (cl_ind, ppc_st, ppc_nd) in ppc_aln_end_idx_dict.items():
        ori_paf_data_list = ori_paf_data_list_data[cl_ind]
        ori_end_idx_dict = ori_end_idx_dict_data[cl_ind]

        aln_paf_data_list = aln_paf_data_list_data[cl_ind]
        aln_end_idx_dict = aln_end_idx_dict_data[cl_ind]

        ori_ctg_name = ori_ctg_name_data[cl_ind][ppc_st]

        ori_st, ori_nd = ori_end_idx_dict[ori_ctg_name]
        aln_st, aln_nd = aln_end_idx_dict[ori_ctg_name]
        assert(aln_st <= ppc_st and ppc_nd <= aln_nd)

        overlap_score = max_overlap(ori_paf_data_list[ori_st:ori_nd+1], [aln_paf_data_list[ppc_st], aln_paf_data_list[ppc_nd]])
        ctgname2overlap[ctg_name] = overlap_score
        
    return ctgname2overlap

def get_overlap_total_score_dict(original_paf_path_list: list) -> dict:
    ori_paf_data_list_data = []
    ori_end_idx_dict_data = []

    for original_paf_path in original_paf_path_list:
        ori_paf_data_list, ori_end_idx_dict = get_qry_cord_data(original_paf_path)

        ori_paf_data_list_data.append(ori_paf_data_list)
        ori_end_idx_dict_data.append(ori_end_idx_dict)

    ctgname2overlap = dict()
    for cl_ind, ori_end_idx_dict in enumerate(ori_end_idx_dict_data):
        for ori_ctg_name in ori_end_idx_dict.keys():
            ori_paf_data_list = ori_paf_data_list_data[cl_ind]
            ori_end_idx_dict = ori_end_idx_dict_data[cl_ind]

            ori_st, ori_nd = ori_end_idx_dict[ori_ctg_name]

            overlap_score = max_overlap(ori_paf_data_list[ori_st:ori_nd+1], [])
            ctgname2overlap[ori_ctg_name] = overlap_score
        
    return ctgname2overlap

def get_not_trust_contig_name(original_paf_path_list: list) -> set:
    not_using_contig = set()
    for file_path in original_paf_path_list:
        with open(file_path, "r") as paf_file:
            for curr_contig in paf_file:
                a = curr_contig.split("\t")
                flag = None
                for i in a:
                    if i.startswith("tp:A:"):
                        if i[5] == 'P':
                            flag = False
                        else:
                            flag = True
                        break
                assert(flag is not None)

                if flag:
                    not_using_contig.add(a[0])
                if int(a[11]) < 60:
                    not_using_contig.add(a[0])
        
    return not_using_contig

def import_repeat_data(file_path : str) -> dict :
    fai_file = open(file_path, "r")
    repeat_data = defaultdict(list)
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        ref_data = (int(temp_list[1]), int(temp_list[2]))
        if abs(ref_data[1] - ref_data[0]) > CENSAT_COMPRESSABLE_THRESHOLD:
            repeat_data[temp_list[0]].append(ref_data)
    fai_file.close()
    return repeat_data

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    chr_data_file.close()
    return chr_len

def extract(contig : list) -> list:
    return [contig[CTG_NAM], contig[CHR_STR], contig[CHR_END]]

def chr2int(x):
    chrXY2int = {'chrX' : 24, 'chrY' : 25}
    if x in chrXY2int:
        return chrXY2int[x]
    else:
        return int(x[3:])
    
def distance_checker(node_a : tuple, node_b : tuple) -> int :
    if max(int(node_a[CHR_STR]), int(node_b[CHR_STR])) < min(int(node_a[CHR_END]), int(node_b[CHR_END])):
        return 0   
    else:
        return min(abs(int(node_b[CHR_STR]) - int(node_a[CHR_END])), abs(int(node_b[CHR_END]) - int(node_a[CHR_STR])))

def distance_checker_tuple(node_a : tuple, node_b : tuple) -> int :
    if max(int(node_a[0]), int(node_b[0])) < min(int(node_a[1]), int(node_b[1])):
        return 0   
    else:
        return min(abs(int(node_b[0]) - int(node_a[1])), abs(int(node_b[1]) - int(node_a[0])))

def alt_simple_ref_len(node) -> int:
    return abs(int(node[CHR_END]) - int(node[CHR_STR]))

def alt_simple_chr_sort_key(chrom: str):
    try:
        return chr2int(chrom)
    except (ValueError, TypeError):
        return INF

def alt_simple_terminal_repeat(node, repeat_label: tuple, chr_len: dict) -> bool:
    if repeat_label[1] not in ("r", "rin"):
        return False
    curr_chr_len = chr_len.get(node[CHR_NAM])
    if curr_chr_len is None:
        return False
    return node[CHR_STR] <= TELOMERE_EXPANSION or node[CHR_END] >= curr_chr_len - TELOMERE_EXPANSION

def trim_alt_simple_terminal_indices(chunks: list, telo_labels: list, repeat_labels: list, chr_len: dict) -> list:
    def is_terminal_chunk(idx: int) -> bool:
        return telo_labels[idx][0] != '0' or alt_simple_terminal_repeat(chunks[idx], repeat_labels[idx], chr_len)

    left = 0
    right = len(chunks) - 1
    while left <= right and is_terminal_chunk(left):
        left += 1
    while right >= left and is_terminal_chunk(right):
        right -= 1
    return list(range(left, right + 1))

def select_alt_simple_major_chroms(indexed_chunks: list) -> set:
    chrom_len = Counter()
    total_len = 0
    for _, node in indexed_chunks:
        node_len = alt_simple_ref_len(node)
        if node_len <= ALT_SIMPLE_MIN_SEGMENT_LEN:
            continue
        chrom_len[node[CHR_NAM]] += node_len
        total_len += node_len

    if total_len == 0:
        return set()

    selected = set()
    acc_len = 0
    for chrom, chrom_span in sorted(chrom_len.items(), key=lambda x: (-x[1], alt_simple_chr_sort_key(x[0]), x[0])):
        selected.add(chrom)
        acc_len += chrom_span
        if acc_len / total_len >= ALT_SIMPLE_MAJOR_CHR_RATIO:
            break
    return selected

def find_alt_simple_transition_candidates(indexed_chunks: list, selected_chroms: set) -> tuple:
    same_dir_candidates = []
    diff_dir_transitions = []
    prev = None
    for chunk_idx, node in indexed_chunks:
        if node[CHR_NAM] not in selected_chroms:
            continue
        if alt_simple_ref_len(node) <= ALT_SIMPLE_MIN_SEGMENT_LEN:
            continue

        curr = (chunk_idx, node)
        if prev is None:
            prev = curr
            continue

        prev_node = prev[1]
        if prev_node[CHR_NAM] == node[CHR_NAM]:
            prev = curr
            continue

        if prev_node[CTG_DIR] == node[CTG_DIR]:
            same_dir_candidates.append((prev, curr))
        else:
            diff_dir_transitions.append((prev, curr))
        prev = curr
    return same_dir_candidates, diff_dir_transitions

def alt_simple_node_locus(node) -> tuple:
    return (node[CHR_NAM], int(node[CHR_STR]), int(node[CHR_END]))

def alt_simple_pair_loci(candidate: tuple) -> tuple:
    return (alt_simple_node_locus(candidate[0][1]), alt_simple_node_locus(candidate[1][1]))

def alt_simple_locus_close(locus_a: tuple, locus_b: tuple, max_dist: int) -> bool:
    if locus_a[0] != locus_b[0]:
        return False
    return distance_checker_tuple((locus_a[1], locus_a[2]), (locus_b[1], locus_b[2])) <= max_dist

def alt_simple_pair_close(pair_a: tuple, pair_b: tuple, max_dist: int) -> bool:
    return (
        alt_simple_locus_close(pair_a[0], pair_b[0], max_dist)
        and alt_simple_locus_close(pair_a[1], pair_b[1], max_dist)
    )

def alt_simple_candidate_near_existing(candidate: tuple, existing_pairs: list, max_dist: int) -> bool:
    cand_pair = alt_simple_pair_loci(candidate)
    cand_pair_rev = (cand_pair[1], cand_pair[0])
    for existing_pair in existing_pairs:
        if alt_simple_pair_close(cand_pair, existing_pair, max_dist):
            return True
        if alt_simple_pair_close(cand_pair_rev, existing_pair, max_dist):
            return True
    return False

def collect_existing_alt_simple_nclose_loci(contig_data: list) -> list:
    existing_pairs = []
    st = 0
    while st < len(contig_data):
        ed = int(contig_data[st][CTG_ENDND])
        curr_contig = contig_data[st:ed + 1]
        st = ed + 1
        if not curr_contig:
            continue
        if curr_contig[0][CTG_NAM].startswith("simple_ctg_alt_"):
            continue
        try:
            curr_type = int(curr_contig[0][CTG_TYP])
        except (TypeError, ValueError):
            continue
        if curr_type not in (1, 2):
            continue

        indexed_chunks = list(enumerate(curr_contig))
        long_chunks = [
            (idx, node)
            for idx, node in indexed_chunks
            if alt_simple_ref_len(node) > ALT_SIMPLE_MIN_SEGMENT_LEN
        ]
        if len(long_chunks) >= 2 and long_chunks[0][1][CHR_NAM] != long_chunks[-1][1][CHR_NAM]:
            existing_pairs.append((
                alt_simple_node_locus(long_chunks[0][1]),
                alt_simple_node_locus(long_chunks[-1][1]),
            ))

        selected_chroms = select_alt_simple_major_chroms(indexed_chunks)
        if len(selected_chroms) < 2:
            continue
        same_dir_candidates, _ = find_alt_simple_transition_candidates(indexed_chunks, selected_chroms)
        for candidate in same_dir_candidates:
            existing_pairs.append(alt_simple_pair_loci(candidate))

    return existing_pairs
    
def overlap_calculator(node_a : tuple, node_b : tuple) -> int :
    return min(abs(node_a[CHR_END] - node_b[CHR_STR]), abs(node_b[CHR_END] - node_a[CHR_STR]))

def inclusive_checker_tuple(tuple_a : tuple, tuple_b : tuple) -> bool :
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False

def chr_correlation_maker(contig_data):
    chr_corr = {}
    chr_rev_corr = {}
    contig_data_size = len(contig_data)
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'f'] = contig_data_size + i - 1
        chr_rev_corr[contig_data_size + i - 1] = 'chr'+str(i)+'f'
    chr_corr['chrXf'] = contig_data_size + CHROMOSOME_COUNT - 1
    chr_corr['chrYf'] = contig_data_size + CHROMOSOME_COUNT - 1
    chr_rev_corr[contig_data_size + CHROMOSOME_COUNT - 1] = 'chrXf'
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'b'] = contig_data_size + CHROMOSOME_COUNT + i - 1
        chr_rev_corr[contig_data_size + CHROMOSOME_COUNT + i - 1] = 'chr'+str(i)+'b'
    chr_corr['chrXb'] = contig_data_size + 2*CHROMOSOME_COUNT - 1
    chr_corr['chrYb'] = contig_data_size + 2*CHROMOSOME_COUNT - 1
    chr_rev_corr[contig_data_size + 2*CHROMOSOME_COUNT - 1] = 'chrXb'

    return chr_corr, chr_rev_corr

def extract_all_repeat_contig(contig_data : list, repeat_data : dict, ctg_index : int, baseline : float = 0) -> set:
    contig_data_size = len(contig_data)
    rpt_con = set()
    ends_map = {
        chrom: [iv[1] for iv in intervals]
        for chrom, intervals in repeat_data.items()
    }
    s = 0
    while s<contig_data_size:
        e = contig_data[s][CTG_ENDND]
        flag = True
        terminal_repeat = True
        total_ref_len = 0
        non_rpt_ref_len = 0
        for i in range(s, e+1):
            total_ref_len += contig_data[i][CHR_END] - contig_data[i][CHR_STR]
            if contig_data[i][ctg_index] == '0':
                non_rpt_ref_len += contig_data[i][CHR_END] - contig_data[i][CHR_STR]

        str_dir = contig_data[s][CTG_DIR]
        str_ref = contig_data[s][CHR_STR if str_dir == '+' else CHR_END]
        str_chr = contig_data[s][CHR_NAM]
        end_dir = contig_data[e][CTG_DIR]
        end_ref = contig_data[e][CHR_END if end_dir == '+' else CHR_STR]
        end_chr = contig_data[e][CHR_NAM]

        for ref, chrom in ((str_ref, str_chr), (end_ref, end_chr)):
            intervals = repeat_data.get(chrom, [])
            ends = ends_map.get(chrom, [])
            if not intervals:
                continue
            c_start, c_end = ref, ref
            idx = bisect.bisect_left(ends, c_start)
            if idx < len(intervals):
                iv_start, iv_end = intervals[idx]
                if iv_start > ref or ref > iv_end:
                    terminal_repeat = False

        if ctg_index == CTG_CENSAT:
            terminal_repeat = True
        if non_rpt_ref_len == 0 or baseline > non_rpt_ref_len/total_ref_len:
            if terminal_repeat:
                rpt_con.add(contig_data[s][CTG_NAM])
            # else:
            #     print(contig_data[s][CTG_NAM])
        s = e+1
    return rpt_con

def check_censat_contig(all_repeat_censat_con : set, ALIGNED_PAF_LOC_LIST : list, ORIGINAL_PAF_LOC_LIST : list, contig_data : list):
    div_repeat_paf_name = div_repeat_paf(ORIGINAL_PAF_LOC_LIST, ALIGNED_PAF_LOC_LIST, contig_data)
    return all_repeat_censat_con & div_repeat_paf_name
    
def extract_bnd_contig(contig_data : list) -> set:
    s = 0
    contig_data_size = len(contig_data)
    bnd_contig = set()
    while s<contig_data_size:
        e = contig_data[s][CTG_ENDND]
        if contig_data[s][CTG_TYP] in {1, 2}: # 2 넣기
            bnd_contig.add(contig_data[s][CTG_NAM])
        s = e+1
    return bnd_contig

def calculate_single_contig_ref_ratio(contig_data : list) -> tuple:
    total_ref_len = 0
    if contig_data[0][CTG_DIR] == '+':
        estimated_ref_len = contig_data[-1][CHR_END] - contig_data[0][CHR_STR]
    else:
        estimated_ref_len = contig_data[0][CHR_END] - contig_data[-1][CHR_STR]
    for node in contig_data:
        total_ref_len += node[CHR_END] - node[CHR_STR]
    return estimated_ref_len/total_ref_len, total_ref_len

# 00_contig_preprocessing
def inclusive_checker_node(contig_node : tuple, telomere_node : tuple) -> bool :
    if int(telomere_node[CHR_STR]) <= int(contig_node[CHR_STR]) and int(contig_node[CHR_END]) <= int(telomere_node[CHR_END]):
        return True
    else:
        return False
    
def telo_distance_checker(node: tuple, telo: tuple) -> int :
    return min(abs(telo[CHR_STR] - node[CHR_END]), abs(telo[CHR_END] - node[CHR_STR]))

def telo_distance_checker_cord(node_st, node_nd, telo_st, telo_nd) -> int :
    return min(abs(telo_st - node_nd), abs(telo_nd - node_st))

def label_node(contig_data : list, telo_data) -> list :
    label = []
    contig_data_size = len(contig_data)
    for i in range(contig_data_size):
        checker = 0
        for j in telo_data[contig_data[i][CHR_NAM]]:
            if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])) == 0:
                inclusive_label = ""
                if inclusive_checker_node(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])):
                    inclusive_label = "in"
                label.append((contig_data[i][CHR_NAM], j[2]+inclusive_label))
                checker = 1
                break
        if checker==0:
            label.append(('0', '0'))
    return label

def label_subtelo_node(contig_data : list, telo_data) -> list :
    label = []
    contig_data_size = len(contig_data)
    for i in range(contig_data_size):
        checker = 0
        for j in telo_data[contig_data[i][CHR_NAM]]:
            if j[2]=='f':
                if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1]+SUBTELOMERE_LENGTH)) == 0:
                    inclusive_label = ""
                    if inclusive_checker_node(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1]+SUBTELOMERE_LENGTH)):
                        inclusive_label = "in"
                    label.append((contig_data[i][CHR_NAM], j[2]+inclusive_label))
                    checker = 1
                    break
            else:
                if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0] - SUBTELOMERE_LENGTH, j[1])) == 0:
                    inclusive_label = ""
                    if inclusive_checker_node(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0] - SUBTELOMERE_LENGTH, j[1])):
                        inclusive_label = "in"
                    label.append((contig_data[i][CHR_NAM], j[2]+inclusive_label))
                    checker = 1
                    break
        if checker==0:
            label.append(('0', '0'))
    return label

def label_repeat_node(contig_data: list, repeat_data : dict, chr_len : dict) -> list:
    labels = []
    ends_map = {
        chrom: [iv[1] for iv in intervals]
        for chrom, intervals in repeat_data.items()
    }

    for contig in contig_data:
        chrom = contig[CHR_NAM]
        intervals = repeat_data.get(chrom, [])
        curr_chr_len = chr_len.get(chrom, 0)
        ends = ends_map.get(chrom, [])
        if not intervals:
            labels.append(('0','0'))
            continue
        c_start, c_end = contig[7], contig[8]
        idx = bisect.bisect_left(ends, c_start)
        interval_list = []
        if SUBTELOMERE_REPEAT_LENGTH > 0:
            interval_list.append((0, SUBTELOMERE_REPEAT_LENGTH))
            interval_list.append((curr_chr_len - SUBTELOMERE_REPEAT_LENGTH, curr_chr_len))
        if idx < len(intervals):
            interval_list.append(intervals[idx])
        r_chk = False
        rin_chk = False
        for iv in interval_list:
            iv_start, iv_end = iv[0], iv[1]
            if distance_checker(contig, (0,0,0,0,0,0,0, iv_start, iv_end)) == 0:
                r_chk = True
                rin_chk = True if inclusive_checker_node(contig, (0,0,0,0,0,0,0, iv_start, iv_end)) else False
        if rin_chk:
            labels.append((chrom, "rin"))
        elif r_chk:
            labels.append((chrom, "r"))
        else:
            labels.append(('0','0'))

    return labels

def preprocess_telo(contig_data : list, node_label : list) -> tuple :
    telo_preprocessed_contig = []
    telo_connect_info = {}
    semi_telomere_contig = []
    report_case = {'A':[], 'B':[], 'C':[], 'ESC':[], 'ALL_TELO_NON_ESC':[]}
    contig_data_size = len(contig_data)
    curr_contig_st = 0
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    for i in range(1, contig_data_size+1):
        if contig_data[i-1][CTG_NAM] != contig_data[i][CTG_NAM]:
            curr_contig_ed = i-1
            front_telo_bound = curr_contig_st
            end_telo_bound = curr_contig_ed
            while front_telo_bound<=curr_contig_ed and node_label[front_telo_bound][0] != '0' :
                front_telo_bound+=1
            while end_telo_bound>=curr_contig_st and node_label[end_telo_bound][0] != '0':
                end_telo_bound-=1
            st = curr_contig_st
            ed = curr_contig_ed
            all_telo = False
            escape = False
            if front_telo_bound == curr_contig_ed+1:
                all_telo = True
                if len(node_label[curr_contig_st][1])==1 \
                and (node_label[curr_contig_st][1] + contig_data[curr_contig_st][CTG_DIR] in ("b+", "f-")):
                    escape = True
                elif len(node_label[curr_contig_ed][1])==1 \
                and (node_label[curr_contig_ed][1] + contig_data[curr_contig_ed][CTG_DIR] in ("b-", "f+")):
                    escape = True
            else:
                escape = True
            if escape:
                if front_telo_bound > curr_contig_st:
                    front_telo_bound-=1
                    # Check if first telomere node is 'Nin' -> Else: escape
                    if len(node_label[curr_contig_st][1])==1 \
                    and (node_label[curr_contig_st][1] + contig_data[curr_contig_st][CTG_DIR] in ("b+", "f-")):
                        report_case['ESC'].append(curr_contig_st)
                        st = curr_contig_st
                    elif len(node_label[curr_contig_st][1])>1:
                        # Next node is connected with telomere node.
                        front_telo_bound+=1
                        if front_telo_bound <= curr_contig_ed:
                            dest = contig_data[front_telo_bound][CHR_NAM]
                            if contig_data[front_telo_bound][CTG_DIR] == '+':
                                dest += 'f'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['A'].append([dest, front_telo_bound])
                            else:
                                dest += 'b'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['A'].append([dest, front_telo_bound])
                        st = front_telo_bound
                    # If boundary node is not "Nin"
                    else:
                        if node_label[curr_contig_st][1] == 'f':
                            if contig_data[front_telo_bound][CTG_DIR]=='+':
                                # boundary node is connected with telomere node.
                                dest = contig_data[front_telo_bound][CHR_NAM]+'f'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['B'].append([dest, front_telo_bound])
                            else:
                                # Treat as "Nin" Case
                                front_telo_bound+=1
                                if front_telo_bound <= curr_contig_ed:
                                    dest = contig_data[front_telo_bound][CHR_NAM]
                                    if contig_data[front_telo_bound][CTG_DIR] == '+':
                                        dest += 'f'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                                    else:
                                        dest += 'b'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                        else:
                            if contig_data[front_telo_bound][CTG_DIR]=='-':
                                dest = contig_data[front_telo_bound][CHR_NAM]+'b'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['B'].append([dest, front_telo_bound])
                            else:
                                front_telo_bound+=1
                                if front_telo_bound <= curr_contig_ed:
                                    dest = contig_data[front_telo_bound][CHR_NAM]
                                    if contig_data[front_telo_bound][CTG_DIR] == '+':
                                        dest += 'f'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                                    else:
                                        dest += 'b'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                        st = front_telo_bound
                
                if end_telo_bound < curr_contig_ed:
                    end_telo_bound+=1
                    # Check if first telomere node is 'Nin' -> Else: escape
                    if len(node_label[curr_contig_ed][1])==1 \
                    and (node_label[curr_contig_ed][1] + contig_data[curr_contig_ed][CTG_DIR] in ("b-", "f+")):
                        report_case['ESC'].append(curr_contig_ed)
                        ed = curr_contig_ed
                    # If boundary node is 'Nin'
                    elif len(node_label[curr_contig_ed][1])>1:
                        # Next node is connected with telomere node.
                        end_telo_bound-=1
                        if end_telo_bound>=curr_contig_st:
                            dest = contig_data[end_telo_bound][CHR_NAM]
                            if contig_data[end_telo_bound][CTG_DIR] == '+':
                                dest += 'b'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['A'].append([dest, end_telo_bound])
                            else:
                                dest += 'f'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['A'].append([dest, end_telo_bound])
                        ed = end_telo_bound
                    # If boundary node is not "Nin"
                    else:
                        if node_label[curr_contig_ed][1] == 'b':
                            if contig_data[end_telo_bound][CTG_DIR]=='+':
                                # boundary node is connected with telomere node.
                                dest = contig_data[end_telo_bound][CHR_NAM]+'b'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['B'].append([dest, end_telo_bound])
                            else:
                                # Treat as "Nin" Case
                                end_telo_bound+=1
                                if end_telo_bound>=curr_contig_st:
                                    dest = contig_data[end_telo_bound][CHR_NAM]
                                    if contig_data[end_telo_bound][CTG_DIR] == '+':
                                        dest += 'b'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                                    else:
                                        dest += 'f'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                        else:
                            if contig_data[end_telo_bound][CTG_DIR]=='-':
                                dest = contig_data[end_telo_bound][CHR_NAM]+'f'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['B'].append([dest, end_telo_bound])
                            else:
                                end_telo_bound+=1
                                if end_telo_bound>=curr_contig_st:
                                    dest = contig_data[end_telo_bound][CHR_NAM]
                                    if contig_data[end_telo_bound][CTG_DIR] == '+':
                                        dest += 'b'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                                    else:
                                        dest += 'f'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                        ed = end_telo_bound
                if all_telo and st>ed:
                    if curr_contig_st <= st and st <= curr_contig_ed:
                        for j in range(curr_contig_st, st+1):
                            telo_preprocessed_contig.append(j)
                    elif curr_contig_st <= ed and ed <= curr_contig_ed:
                        for j in range(ed, curr_contig_ed+1):
                            telo_preprocessed_contig.append(j)
                else:
                    for j in range(st, ed+1):
                        telo_preprocessed_contig.append(j)
                
            else:
                report_case['ALL_TELO_NON_ESC'].append(curr_contig_st)
            # initialize
            curr_contig_st = curr_contig_ed+1
    contig_data = contig_data[0:-1]
    return telo_preprocessed_contig, report_case, telo_connect_info

def subtelo_cut(contig_data : list, node_label : list, subnode_label : list) -> list :
    telo_preprocessed_contig = []
    contig_data_size = len(contig_data)
    curr_contig_st = 0
    tcnt = 0
    idx = 0
    while curr_contig_st < contig_data_size:
        curr_contig_ed = contig_data[curr_contig_st][CTG_ENDND]
        front_telo_bound = curr_contig_st
        end_telo_bound = curr_contig_ed
        front_telcon = False
        end_telcon = False
        front_block_contain_telo = False
        end_block_contain_telo = False
        st = curr_contig_st
        ed = curr_contig_ed
        if contig_data[curr_contig_st][CTG_TELCON] != '0':
            front_telcon = True
        if contig_data[curr_contig_ed][CTG_TELCON] != '0':
            end_telcon = True
        front_dest = '0'
        end_dest = '0'
        if not front_telcon:
            while front_telo_bound<=curr_contig_ed and subnode_label[front_telo_bound][0] != '0' :
                if node_label[front_telo_bound][0] != '0':
                    front_block_contain_telo = True
                front_telo_bound+=1
        if not end_telcon:
            while end_telo_bound>=curr_contig_st and subnode_label[end_telo_bound][0] != '0':
                if node_label[end_telo_bound][0] != '0':
                    end_block_contain_telo = True
                end_telo_bound-=1
        for_cut = False
        bak_cut = False
        if front_telo_bound > curr_contig_st and front_block_contain_telo:
            st = front_telo_bound
            if contig_data[st][CTG_DIR]=='+':
                front_dest = contig_data[st][CHR_NAM] + 'f'
                for_cut = True
            else:
                front_dest = contig_data[st][CHR_NAM] + 'b'
                for_cut = True
        if end_telo_bound < curr_contig_ed and end_block_contain_telo:
            ed = end_telo_bound
            if contig_data[ed][CTG_DIR]=='+':
                end_dest = contig_data[ed][CHR_NAM] + 'b'
                bak_cut = True
            else:
                end_dest = contig_data[ed][CHR_NAM] + 'f'
                bak_cut = True
        if (for_cut or bak_cut) and st < ed:
            l = ed - st
            tcnt+=1
            if for_cut:
                temp_list = copy.deepcopy(contig_data[st])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                temp_list[CTG_TELCON] = front_dest
                telo_preprocessed_contig.append(temp_list)
            else:
                temp_list = copy.deepcopy(contig_data[st])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                telo_preprocessed_contig.append(temp_list)
            for i in range(st+1, ed):
                temp_list = copy.deepcopy(contig_data[i])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                telo_preprocessed_contig.append(temp_list)
            if bak_cut:
                temp_list = copy.deepcopy(contig_data[ed])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                temp_list[CTG_TELCON] = end_dest
                telo_preprocessed_contig.append(temp_list)
            else:
                temp_list = copy.deepcopy(contig_data[ed])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                telo_preprocessed_contig.append(temp_list)
            idx += l + 1
        
        curr_contig_st = curr_contig_ed+1
    return telo_preprocessed_contig


def initial_graph_build(contig_data : list, telo_data : dict) -> list :
    '''
    Initialize
    '''
    contig_data_size = len(contig_data)
    chr_corr, chr_rev_corr = chr_correlation_maker(contig_data)

    adjacency = [[[] for _ in range(contig_data_size+CHROMOSOME_COUNT*2)], [[] for _ in range(contig_data_size+CHROMOSOME_COUNT*2)]]
    '''
    Algorithm
    '''
    curr_contig_st = 0
    while curr_contig_st<contig_data_size:
        curr_contig_ed = contig_data[curr_contig_st][CTG_ENDND]
        if curr_contig_st != curr_contig_ed:
            if contig_data[curr_contig_st][CTG_TELCON] != '0':
                dest = chr_corr[contig_data[curr_contig_st][CTG_TELCON]]
                adjacency[DIR_FOR][dest].append([DIR_FOR, curr_contig_st, 0])
                adjacency[DIR_FOR][curr_contig_st].append([DIR_FOR, dest, 0])
                adjacency[DIR_BAK][curr_contig_st].append([DIR_FOR, dest, 0])
            if contig_data[curr_contig_ed][CTG_TELCON] != '0':
                dest = chr_corr[contig_data[curr_contig_ed][CTG_TELCON]]
                adjacency[DIR_FOR][curr_contig_ed].append([DIR_FOR, dest, 0])
                adjacency[DIR_BAK][curr_contig_ed].append([DIR_FOR, dest, 0])
                adjacency[DIR_FOR][dest].append([DIR_BAK, curr_contig_ed, 0])
        else:
            if contig_data[curr_contig_st][CTG_TELCON] != '0':
                if contig_data[curr_contig_st][CTG_TELCON][-1]=='f':
                    dest = chr_corr[contig_data[curr_contig_st][CTG_TELCON]]
                    adjacency[DIR_FOR][dest].append([DIR_FOR, curr_contig_st, 0])
                    adjacency[DIR_BAK][curr_contig_st].append([DIR_FOR, dest, 0])
                    adjacency[DIR_FOR][curr_contig_st].append([DIR_FOR, dest, 0])
                else:
                    dest = chr_corr[contig_data[curr_contig_ed][CTG_TELCON]]
                    adjacency[DIR_FOR][curr_contig_ed].append([DIR_FOR, dest, 0])
                    adjacency[DIR_BAK][curr_contig_ed].append([DIR_FOR, dest, 0])
                    adjacency[DIR_FOR][dest].append([DIR_BAK, curr_contig_ed, 0])

        curr_contig_st = curr_contig_ed + 1
    '''
    If telomere node has 0 connection, connect with closest node with same chromosome type.
    '''
    end_node_list = set()
    for i in contig_data:
        end_node_list.add((i[CTG_STRND], i[CTG_ENDND])) 
    
    now_telo_list = list(chr_corr.keys())

    if no_chrY:
        now_telo_list = [i for i in now_telo_list if 'chrY' not in i]

    for now_telo in now_telo_list:
        now_telo_chr = now_telo[:-1]
        curr_telo_set = set()
        i = chr_corr[now_telo]
        flag = False

        mini_telo_dist = INF
        for j in range(2):
            for connected_node in adjacency[j][i]:
                if contig_data[connected_node[1]][CHR_NAM] == now_telo_chr:
                    # 10K 이내면 없애기
                    if i < contig_data_size + CHROMOSOME_COUNT:
                        mini_telo_dist = min(mini_telo_dist, 0 if contig_data[connected_node[1]][CHR_STR] < telo_data[now_telo][0] else contig_data[connected_node[1]][CHR_STR] - telo_data[now_telo][0])
                        if contig_data[connected_node[1]][CHR_STR] < telo_data[now_telo][0]+FORCE_TELOMERE_THRESHOLD:
                            flag = True
                    else:
                        mini_telo_dist = min(mini_telo_dist, 0 if contig_data[connected_node[1]][CHR_END] > telo_data[now_telo][1] else telo_data[now_telo][1] - contig_data[connected_node[1]][CHR_END])
                        if contig_data[connected_node[1]][CHR_END] > telo_data[now_telo][1]-FORCE_TELOMERE_THRESHOLD:
                            flag = True
                    curr_telo_set.add(connected_node[1])
        
        if flag:
            continue

        telo_dist = mini_telo_dist
        telo_connect_node = INF
        telo_dir = 0
        temp_contig = (0, 0, 0, 0, 0, 0, 0, telo_data[now_telo][0], telo_data[now_telo][1])
        
        #우리가 연결하는 텔로미어 노드
        if now_telo[-1]=='f':
            for st, ed in end_node_list:
                # 텔로미어와 겹치지 않으며, 배제된 노드가 아니고, 중복이 아니며, + 방향이고, 거리가 갱신가능한 경우
                if contig_data[st][CHR_NAM] == now_telo_chr \
                and st not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[st], temp_contig) \
                and contig_data[st][CTG_DIR]=='+':
                    telo_connect_node = st
                    telo_dist = telo_distance_checker(contig_data[st], temp_contig)
                    telo_dir = '+'
                if contig_data[ed][CHR_NAM] == now_telo_chr \
                and ed not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[ed], temp_contig) \
                and contig_data[ed][CTG_DIR]=='-':
                    telo_connect_node = ed
                    telo_dist = telo_distance_checker(contig_data[ed], temp_contig)
                    telo_dir = '-'
            if telo_connect_node != INF:
                if telo_dir == '+':
                    adjacency[DIR_FOR][i].append([DIR_FOR, telo_connect_node, telo_dist])
                    adjacency[DIR_BAK][telo_connect_node].append([DIR_FOR, i, telo_dist])
                else:
                    adjacency[DIR_FOR][i].append([DIR_BAK, telo_connect_node, telo_dist])
                    adjacency[DIR_FOR][telo_connect_node].append([DIR_FOR, i, telo_dist])
        else:
            for st, ed in end_node_list:
                if contig_data[st][CHR_NAM] == now_telo_chr \
                and st not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[st], temp_contig) \
                and contig_data[st][CTG_DIR]=='-':
                    telo_connect_node = st
                    telo_dist = telo_distance_checker(contig_data[st], temp_contig)
                    telo_dir = '-'
                if contig_data[ed][CHR_NAM] == now_telo_chr \
                and ed not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[ed], temp_contig) \
                and contig_data[ed][CTG_DIR]=='+':
                    telo_connect_node = ed
                    telo_dist = telo_distance_checker(contig_data[ed], temp_contig)
                    telo_dir = '+'
            if telo_connect_node != INF:
                if telo_dir == '+':
                    adjacency[DIR_FOR][i].append([DIR_BAK, telo_connect_node, telo_dist])
                    adjacency[DIR_FOR][telo_connect_node].append([DIR_FOR, i, telo_dist])
                    adjacency[DIR_BAK][telo_connect_node].append([DIR_FOR, i, telo_dist])
                else:
                    adjacency[DIR_FOR][i].append([DIR_FOR, telo_connect_node, telo_dist])
                    adjacency[DIR_BAK][telo_connect_node].append([DIR_FOR, i, telo_dist])
                    adjacency[DIR_FOR][telo_connect_node].append([DIR_FOR, i, telo_dist])

    return adjacency

def edge_optimization(contig_data : list, contig_adjacency : list, telo_dict : dict,
                      asm2cov : dict, excluded_telomere_origins=None) -> tuple :

    contig_data_size = len(contig_data)
    excluded_telomere_origins = excluded_telomere_origins or set()
    excluded_candidate_nodes = set()
    chr_corr, chr_rev_corr = chr_correlation_maker(contig_data)
    contig_pair_nodes = defaultdict(list)
    telo_coverage = Counter()
    optimized_adjacency = [[[] for _ in range(contig_data_size + CHROMOSOME_COUNT*2)], [[] for _ in range(contig_data_size + CHROMOSOME_COUNT*2)]]
    for _ in range(2):
        for i in range(len(contig_data)+CHROMOSOME_COUNT*2):
            for edge in contig_adjacency[_][i]:
                if edge[1] >= contig_data_size or i >= contig_data_size:
                    optimized_adjacency[_][i].append(edge)
                    continue
                if contig_data[i][CTG_STRND]>contig_data[edge[1]][CTG_STRND]:
                    contig_pair_nodes[(contig_data[edge[1]][CTG_NAM], contig_data[i][CTG_NAM])] \
                    .append([edge[1], i])
                else:
                    contig_pair_nodes[(contig_data[i][CTG_NAM], contig_data[edge[1]][CTG_NAM])]\
                    .append([i, edge[1]])
    minmaxnode_contigpair = {}
    for pair in contig_pair_nodes:
        min_first_contig = INF
        max_first_contig = 0
        min_second_contig = INF
        max_second_contig = 0
        for edge in contig_pair_nodes[pair]:
            min_first_contig = min(min_first_contig, edge[0])
            max_first_contig = max(max_first_contig, edge[0])
            min_second_contig = min(min_second_contig, edge[1])
            max_second_contig = max(max_second_contig, edge[1])
            minmaxnode_contigpair[pair] = [min_first_contig, max_first_contig, min_second_contig, max_second_contig]
    
    contig_data_size = len(contig_data)
    for i in range(2):
        for j in range(contig_data_size):
            first_contig_name = contig_data[j][CTG_NAM]
            fcn_int = contig_data[j][CTG_STRND]
            for edge in contig_adjacency[i][j]:
                if edge[1] >= contig_data_size:
                    continue
                second_contig_name = contig_data[edge[1]][CTG_NAM]
                scn_int = contig_data[edge[1]][CTG_STRND]
                if fcn_int == scn_int:
                    optimized_adjacency[i][j].append(edge)
                elif fcn_int < scn_int:
                    if j in minmaxnode_contigpair[(first_contig_name, second_contig_name)][0:2]\
                    or edge[1] in minmaxnode_contigpair[(first_contig_name, second_contig_name)][2:4]:
                        optimized_adjacency[i][j].append(edge)
                else:
                    if edge[1] in minmaxnode_contigpair[(second_contig_name, first_contig_name)][0:2] \
                    or j in minmaxnode_contigpair[(second_contig_name, first_contig_name)][2:4]:
                        optimized_adjacency[i][j].append(edge)


    for i in range(contig_data_size, contig_data_size+2*CHROMOSOME_COUNT):
        telo_name = chr_rev_corr[i]
        telo_range = telo_dict[telo_name]
        for j in range(2):
            using_edge = []
            now_edge = [-1, 0, 0]
            curr_coverage = 0
            for edge in optimized_adjacency[j][i]:
                cl_ind, c_ind = map(int, contig_data[edge[1]][CTG_GLOBALIDX].split('.'))
                if (cl_ind, c_ind) in excluded_telomere_origins:
                    excluded_candidate_nodes.add(edge[1])
                    continue
                if cl_ind < 2:
                    name = ori_ctg_name_data[cl_ind][c_ind]
                else:
                    name = contig_data[edge[1]][CTG_NAM]
                    asm2cov[name] = -1
                if telo_name[-1]=='f':
                    if contig_data[edge[1]][CHR_STR]<=TELOMERE_CLUSTER_THRESHOLD:
                        if now_edge[0]<0:
                            curr_coverage += asm2cov[name]
                            now_edge = edge
                        else:
                            if contig_data[now_edge[1]][CHR_STR] > contig_data[edge[1]][CHR_STR]:
                                curr_coverage += asm2cov[name]
                                now_edge = edge
                            elif contig_data[now_edge[1]][CHR_STR] == contig_data[edge[1]][CHR_STR]:
                                if contig_data[now_edge[1]][CTG_LEN] > contig_data[edge[1]][CTG_LEN]:
                                    curr_coverage += asm2cov[name]
                                    now_edge = edge
                    else:
                        telo_compress_flag = False
                        for existing_edge in using_edge:
                            if distance_checker(contig_data[existing_edge[1]], contig_data[edge[1]]) < TELOMERE_COMPRESS_RANGE:
                                telo_compress_flag = True
                                telo_coverage[existing_edge[1]] += asm2cov[name]
                                break
                        if not telo_compress_flag:
                            telo_coverage[edge[1]] += asm2cov[name]
                            using_edge.append(edge)
                else:
                    if contig_data[edge[1]][CHR_END]>=telo_range[1]-TELOMERE_CLUSTER_THRESHOLD:
                        if now_edge[0]<0:
                            curr_coverage += asm2cov[name]
                            now_edge = edge
                        else:
                            if contig_data[now_edge[1]][CHR_END] < contig_data[edge[1]][CHR_END]:
                                curr_coverage += asm2cov[name]
                                now_edge = edge
                            elif contig_data[now_edge[1]][CHR_END] == contig_data[edge[1]][CHR_END]:
                                if contig_data[now_edge[1]][CTG_LEN] > contig_data[edge[1]][CTG_LEN]:
                                    curr_coverage += asm2cov[name]
                                    now_edge = edge
                    else:
                        telo_compress_flag = False
                        for existing_edge in using_edge:
                            if distance_checker(contig_data[existing_edge[1]], contig_data[edge[1]]) < TELOMERE_COMPRESS_RANGE:
                                telo_compress_flag = True
                                telo_coverage[existing_edge[1]] += asm2cov[name]
                                break
                        if not telo_compress_flag:
                            telo_coverage[edge[1]] += asm2cov[name]
                            using_edge.append(edge)
            if now_edge == [-1, 0, 0]:
                pass
            else:
                # telo_coverage[now_edge[1]] += curr_coverage # norm telomere coverage
                using_edge.append(now_edge)
            optimized_adjacency[j][i] = using_edge
        
    telo_connected_set = set()
    telo_connected_dict = dict()
    telo_connected_graph_dict = defaultdict(list)
    for j in range(2):
        for i in range(contig_data_size, contig_data_size + 2*(CHROMOSOME_COUNT)):
            for k in optimized_adjacency[j][i]:
                telo_connected_set.add(k[1])
                telo_connected_dict[k[1]] = chr_rev_corr[i]
                telo_connected_graph_dict[chr_rev_corr[i]].append(k)

    if excluded_candidate_nodes:
        logging.info(
            f"Skipped {len(excluded_candidate_nodes)} multi-end-aligned PAF nodes "
            "during telomere edge optimization"
        )

    return telo_connected_set, telo_connected_dict, telo_connected_graph_dict, telo_coverage

def _flip_ctg_dir(ctg_dir):
    return '-' if ctg_dir == '+' else '+'

def _cen_fragment_target_dir_from_meta(cen_fragment_meta, chrom):
    return '-' if cen_fragment_meta[chrom]['dir'] else '+'

def _normalized_telo_censat_dir(telo_name, ctg_dir):
    if telo_name[-1] == 'f':
        return ctg_dir
    return _flip_ctg_dir(ctg_dir)

def filter_telomere_connected_cen_fragment_mismatch(
    contig_data,
    telo_connected_graph_dict,
    cen_fragment_meta,
    stage_name,
):
    filtered_graph_dict = defaultdict(list)
    filtered_dict = {}
    filtered_set = set()
    removed = 0

    for telo_name, edge_list in telo_connected_graph_dict.items():
        for edge in edge_list:
            node_idx = edge[1]
            contig = contig_data[node_idx]
            chrom = contig[CHR_NAM]
            mismatch = False
            if contig[CTG_CENSAT] != '0' and chrom in cen_fragment_meta and telo_name[-1] in ('f', 'b'):
                norm_dir = _normalized_telo_censat_dir(telo_name, contig[CTG_DIR])
                mismatch = norm_dir != _cen_fragment_target_dir_from_meta(cen_fragment_meta, chrom)

            if mismatch:
                removed += 1
                continue

            filtered_graph_dict[telo_name].append(edge)
            filtered_dict[node_idx] = telo_name
            filtered_set.add(node_idx)

    logging.info(
        f"Removed {removed} {stage_name} telomere-connected censat nodes "
        f"with cen_fragment direction mismatch"
    )
    return filtered_set, filtered_dict, filtered_graph_dict

def calc_ratio(contig_data : list) -> dict:
    contig_data_size = len(contig_data)
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0))
    total_ref_len = 0
    curr_contig_first_fragment = contig_data[0]
    ref_qry_ratio = {}
    for i in range(1, contig_data_size+1):
        total_ref_len += contig_data[i-1][CHR_END] - contig_data[i-1][CHR_STR]
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CTG_DIR] == '+':
                estimated_ref_len = curr_contig_end_fragment[CHR_END] - curr_contig_first_fragment[CHR_STR]
            else:
                estimated_ref_len = curr_contig_first_fragment[CHR_END] - curr_contig_end_fragment[CHR_STR]
            try:
                ref_qry_ratio[curr_contig_name] = estimated_ref_len / total_ref_len
            except:
                if total_ref_len >=0:
                    ref_qry_ratio[curr_contig_name] = INF
                else:
                    ref_qry_ratio[curr_contig_name] = -INF
            total_ref_len = 0
            curr_contig_first_fragment = contig_data[i]
    contig_data = contig_data[:-1]
    return ref_qry_ratio            

def calc_chukji(contig_data : list) -> tuple:
    contig_data_size = len(contig_data)
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0))
    total_ref_len = 0
    curr_contig_first_fragment = contig_data[0]
    ref_qry_ratio = {}
    chukji = {}
    for i in range(1, contig_data_size+1):
        total_ref_len += contig_data[i-1][CHR_END] - contig_data[i-1][CHR_STR]
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CTG_DIR] == '+':
                estimated_ref_len = curr_contig_end_fragment[CHR_END] - curr_contig_first_fragment[CHR_STR]
            else:
                estimated_ref_len = curr_contig_first_fragment[CHR_END] - curr_contig_end_fragment[CHR_STR]
            try:
                ref_qry_ratio[curr_contig_name] = estimated_ref_len / total_ref_len
                chukji[curr_contig_name] = abs(estimated_ref_len)
            except:
                if total_ref_len >=0:
                    ref_qry_ratio[curr_contig_name] = INF
                else:
                    ref_qry_ratio[curr_contig_name] = -INF
            total_ref_len = 0
            curr_contig_first_fragment = contig_data[i]
    contig_data = contig_data[:-1]
    return ref_qry_ratio, chukji          

def find_mainflow(contig_data : list) -> dict:
    contig_data_size = len(contig_data)
    mainflow_dict = {}
    st = 0
    while st<contig_data_size:
        ed = contig_data[st][CTG_ENDND]
        ref_length_counter = Counter()
        for i in range(st, ed+1):
            ref_length_counter[(contig_data[i][CTG_DIR], contig_data[i][CHR_NAM])]\
            +=contig_data[i][CTG_END]-contig_data[i][CTG_STR]
        max_count = 0
        max_chr = (0, 0)
        total_len = 0
        for i in ref_length_counter:
            total_len += ref_length_counter[i]
            if ref_length_counter[i]>max_count:
                max_count = ref_length_counter[i]
                max_chr = i
        mainflow_dict[contig_data[st][CTG_NAM]] = max_chr
        st = ed+1
    return mainflow_dict

def pipeline_preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list, telo_connect_info : set) -> list :
    checker = 0
    contig_data_size = len(contig_data)-1
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    is_front_back_repeat = False
    chrM_flag = False
    idx = 0
    cnt = 0
    len_count = Counter()
    for i in range(1, contig_data_size+1):
        len_count[contig_data[i-1][CTG_NAM]]+=contig_data[i-1][CHR_END]-contig_data[i-1][CHR_STR]
        cnt+=1
        if (i-1) in telo_connect_info:
            is_telo = True
        if contig_data[i-1][CHR_NAM]=='chrM':
            chrM_flag = True
        # contig 넘어갈 때:
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            if repeat_label[i-1][0]!='0':
                is_front_back_repeat = True
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CHR_NAM] != curr_contig_end_fragment[CHR_NAM]:
                checker = 1
            elif curr_contig_first_fragment[CTG_DIR] != curr_contig_end_fragment[CTG_DIR]:
                checker = 2
            else:
                if is_front_back_repeat:
                    bound = RPT_BND_CONTIG_BOUND
                else:
                    bound = BND_CONTIG_BOUND
                if abs(ref_qry_ratio[curr_contig_name]-1)< bound:
                    if contig_data[i-1][CTG_LEN] > CONTIG_MINIMUM_SIZE or is_telo:
                        checker = 3
                    elif cnt==1:
                        checker = 3
                else:
                    checker = 4
            if chrM_flag:
                checker = 0
            if checker>0:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            # initialize
            curr_contig_first_fragment = contig_data[i]
            cnt = 0
            checker = 0
            is_telo = False
            is_front_back_repeat = False
            chrM_flag = False
            if i < contig_data_size and repeat_label[i][0]!='0':
                is_front_back_repeat = True
    contig_data = contig_data[:-1]
    return [using_contig_list, contig_type, contig_terminal_node, len_count]
    

def preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list, telo_connect_info : set) -> list :
    checker = 0
    contig_data_size = len(contig_data)-1
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    is_front_back_repeat = False
    chrM_flag = False
    first_idx = 0
    idx = 0
    cnt = 0
    len_count = Counter()
    for i in range(1, contig_data_size+1):
        len_count[contig_data[i-1][CTG_NAM]]+=contig_data[i-1][CHR_END]-contig_data[i-1][CHR_STR]
        cnt+=1
        if (i-1) in telo_connect_info:
            is_telo = True
        if contig_data[i-1][CHR_NAM]=='chrM':
            chrM_flag = True
        # contig 넘어갈 때:
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            splited_type4 = False
            curr_contig_name = contig_data[i-1][CTG_NAM]
            if repeat_label[i-1][0]!='0':
                is_front_back_repeat = True
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CHR_NAM] != curr_contig_end_fragment[CHR_NAM]:
                checker = 1
            elif curr_contig_first_fragment[CTG_DIR] != curr_contig_end_fragment[CTG_DIR]:
                checker = 2
            else:
                if is_front_back_repeat:
                    bound = RPT_BND_CONTIG_BOUND
                else:
                    bound = BND_CONTIG_BOUND
                if abs(ref_qry_ratio[curr_contig_name]-1)< bound:
                    if contig_data[i-1][CTG_LEN] > CONTIG_MINIMUM_SIZE or is_telo:
                        checker = 3
                    elif cnt==1:
                        checker = 3
                else:
                    checker = 4
            if chrM_flag:
                checker = 0
            if checker==3:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            # initialize
            curr_contig_first_fragment = contig_data[i]
            first_idx = i
            cnt = 0
            checker = 0
            is_telo = False
            is_front_back_repeat = False
            chrM_flag = False
            if i < contig_data_size and repeat_label[i][0]!='0':
                is_front_back_repeat = True
    contig_data = contig_data[:-1]
    return [using_contig_list, contig_type, contig_terminal_node, len_count]

def similar_check(v1, v2, ratio=TYPE2_SIM_COMPARE_RAITO):
    assert(v1 >= 0 and v2 >= 0)
    mi, ma = sorted([v1, v2])
    return False if mi == 0 else (ma / mi <= ratio) or ma-mi < NCLOSE_SIM_DIFF_THRESHOLD

def exist_near_bnd_point(chrom, inside_st):
    # subset of df for the given chromosome
    df_chr = df[df['chr'] == chrom]

    def mean_depth(start, end):
        """Return mean meandepth over windows overlapping [start, end)."""
        mask = (df_chr['nd'] > start) & (df_chr['st'] < end)
        return df_chr.loc[mask, 'meandepth'].mean()

    # for inside_st
    st_depth = mean_depth(inside_st - TYPE2_FLANKING_LENGTH, inside_st)
    nd_depth = mean_depth(inside_st, inside_st + TYPE2_FLANKING_LENGTH)

    # print(chrom, inside_st, inside_nd, not similar_check(st_depth, nd_depth))
    return not similar_check(st_depth, nd_depth)

def check_depth_near_bnd(st_data, nd_data):
    def mean_depth(start, end):
        """Return mean meandepth over windows overlapping [start, end)."""
        mask = (df_chr['nd'] > start) & (df_chr['st'] < end)
        return df_chr.loc[mask, 'meandepth'].mean()
    depth = []
    for data in [st_data, nd_data]:
        chrom = data[0]
        ref = data[1]
        df_chr = df[df['chr'] == chrom]
        depth.append(mean_depth(ref - TYPE2_FLANKING_LENGTH/2, ref + TYPE2_FLANKING_LENGTH/2))
    if similar_check(*depth, ratio = NCLOSE_SIM_COMPARE_RAITO):
        return True
    else:
        return False
        

def censat_overlap_check(censat_dict, chrom, inside_st, inside_nd):
    if chrom not in censat_dict.keys():
        return False
    for st, nd in censat_dict[chrom]:
        if not (nd < inside_st or inside_nd < st):
            return True
    return False

def alt_preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list, telo_connect_info : set, telo_dict : dict, censat_dict : dict) -> list :
    checker = 0
    contig_data_size = len(contig_data)-1
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    is_front_back_repeat = False
    chrM_flag = False
    len_count = Counter()
    idx = 0
    cnt = 0
    telo_node_count = 0
    all_mapq_0 = True
    
    telo_inside_dict = dict()
    for chr_name, data in telo_dict.items():
        for s, e, d in data:
            key = chr_name + d
            func = max if d == 'f' else min
            if key not in telo_inside_dict:
                telo_inside_dict[key] = (s, e)
            
            telo_inside_dict[key] = func(telo_inside_dict[key], (s, e), key=lambda t: t[0])
    
    telo_safe_dict = defaultdict(list)
    for i in range(1, contig_data_size+1):
        len_count[contig_data[i-1][CTG_NAM]]+=contig_data[i-1][CHR_END]-contig_data[i-1][CHR_STR]
        cnt+=1
        if (i-1) in telo_connect_info:
            telo_node_count += 1
            is_telo = True
        if contig_data[i-1][CHR_NAM] == 'chrM':
            chrM_flag = True
        if contig_data[i-1][CTG_MAPQ] != 0:
            all_mapq_0 = False
        # contig 넘어갈 때:
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            if repeat_label[i-1][0]!='0':
                is_front_back_repeat = True
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CHR_NAM] != curr_contig_end_fragment[CHR_NAM]:
                checker = 1
            elif curr_contig_first_fragment[CTG_DIR] != curr_contig_end_fragment[CTG_DIR]:
                checker = 2
            else:
                if is_telo:
                    checker = 5
                else:
                    if is_front_back_repeat:
                        bound = RPT_BND_CONTIG_BOUND
                    else:
                        bound = BND_CONTIG_BOUND
                    
                    if abs(ref_qry_ratio[curr_contig_name]-1) >= bound:
                        checker = 4
                    else:
                        checker = 3

            if chrM_flag:
                checker = 0
            elif checker == 2:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
                
            elif (checker>0 and checker != 3 and checker != 5):
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            elif is_telo:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+telo_node_count-1)
                idx+=telo_node_count

            if checker == 3 and cnt == 1:
                now_chr = contig_data[i-1][CHR_NAM]
                st_nd = telo_inside_dict[now_chr + 'f'][1]
                nd_st = telo_inside_dict[now_chr + 'b'][0]

                st, nd = contig_data[i-1][CHR_STR], contig_data[i-1][CHR_END]
                if st_nd < st and nd < nd_st:
                    s, n = telo_inside_dict[now_chr + 'f']
                    telo_safe_dict[now_chr + 'f'].append((telo_distance_checker_cord(st, nd, s, n), curr_contig_name, i-1))

                    s, n = telo_inside_dict[now_chr + 'b']
                    telo_safe_dict[now_chr + 'b'].append((telo_distance_checker_cord(st, nd, s, n), curr_contig_name, i-1))


            curr_contig_first_fragment = contig_data[i]
            all_mapq_0 = True
            cnt = 0
            checker = 0
            telo_node_count = 0
            is_telo = False
            is_front_back_repeat = False
            chrM_flag = False
            if i < contig_data_size and repeat_label[i][0]!='0':
                is_front_back_repeat = True
    contig_data = contig_data[:-1]

    using_type3_contig_data = set()
    for k, node_list in telo_safe_dict.items():
        dis, curr_contig_name, ind = min(node_list, key=lambda t: t[0])
        using_type3_contig_data.add((curr_contig_name, ind))
    using_type3_contig_data = list(using_type3_contig_data)
    using_type3_contig_list = []
    using_type3_contig_data.sort(key=lambda t: t[1])
    for curr_contig_name, _ in using_type3_contig_data:
        if curr_contig_name not in using_contig_list:
            using_type3_contig_list.append(curr_contig_name)

            cnt = 1
            contig_type[curr_contig_name] = 3
            contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
            idx+=cnt
    
    return [using_contig_list, using_type3_contig_list, contig_type, contig_terminal_node, len_count]

def preprocess_repeat(contig_data : list) -> list:
    repeat_preprocessed_contig = []
    contig_data_size = len(contig_data)
    curr_contig_st = 0
    idx = 0
    while curr_contig_st<contig_data_size:
        curr_contig_ed = contig_data[curr_contig_st][CTG_ENDND]
        front_telo_bound = curr_contig_st
        end_telo_bound = curr_contig_ed
        st = front_telo_bound
        ed = end_telo_bound
        if contig_data[curr_contig_st][CTG_TYP] != 3:
            while front_telo_bound <= curr_contig_ed and contig_data[front_telo_bound][CTG_TELCON] != '0':
                front_telo_bound+=1
            while end_telo_bound >= curr_contig_st and contig_data[end_telo_bound][CTG_TELCON] != '0':
                end_telo_bound-=1
            front_repeat_bound = front_telo_bound
            end_repeat_bound = end_telo_bound
            while front_repeat_bound<=curr_contig_ed \
                  and contig_data[front_repeat_bound][CTG_RPTCHR] != '0' \
                  and contig_data[front_repeat_bound][CTG_MAPQ] < MAPQ_BOUND:
                front_repeat_bound+=1
            if front_repeat_bound <= curr_contig_ed and front_repeat_bound > front_telo_bound:
                if contig_data[front_repeat_bound][CTG_RPTCHR] == '0':
                    front_repeat_bound-=1
            while end_repeat_bound>=curr_contig_st \
                  and contig_data[end_repeat_bound][CTG_RPTCHR] != '0' \
                  and contig_data[end_repeat_bound][CTG_MAPQ] < MAPQ_BOUND:
                end_repeat_bound-=1
            if end_repeat_bound >= curr_contig_st and end_repeat_bound < end_telo_bound:
                if contig_data[end_repeat_bound][CTG_RPTCHR] == '0':
                    end_repeat_bound+=1
            if front_repeat_bound >= front_telo_bound:
                st = front_repeat_bound
            if end_repeat_bound <= end_telo_bound:
                ed = end_repeat_bound
        repeat_temp_set = set()
        for i in range(curr_contig_st, front_telo_bound):
            repeat_temp_set.add(i)
        for i in range(st, ed+1):
            repeat_temp_set.add(i)
        for i in range(end_telo_bound+1, curr_contig_ed+1):
            repeat_temp_set.add(i)
        for i in sorted(list(repeat_temp_set)):
            repeat_preprocessed_contig.append(contig_data[i])
        for i in range(idx, len(repeat_preprocessed_contig)):
            repeat_preprocessed_contig[i][CTG_STRND] = idx
            repeat_preprocessed_contig[i][CTG_ENDND] = len(repeat_preprocessed_contig)-1
        
        idx = len(repeat_preprocessed_contig)
        curr_contig_st = curr_contig_ed+1
    return repeat_preprocessed_contig

def weighted_avg_meandepth(chrom_df, region_start, region_end):
    overlapping = chrom_df[(chrom_df['nd'] >= region_start) & (chrom_df['st'] <= region_end)]
    if overlapping.empty:
        return None
    total_weight = 0
    weighted_sum = 0
    for _, row in overlapping.iterrows():
        overlap_start = max(row['st'], region_start)
        overlap_end = min(row['nd'], region_end)
        if overlap_start <= overlap_end:
            length = overlap_end - overlap_start + 1
            weighted_sum += row['meandepth'] * length
            total_weight += length
    return weighted_sum / total_weight if total_weight > 0 else None


def find_breakend_centromere(
    repeat_censat_data : dict,
    chr_len : dict,
    df : pd.DataFrame,
    raw_nclose_nodes : dict = None,
    contig_data : list = None,
    log_context : str = "",
):
    results = []

    # Remove unused Y
    ydf = df.query('chr == "chrY"')
    ydepth = np.mean(ydf['meandepth'].to_numpy())

    depth = np.mean(df['meandepth'].to_numpy())

    sd, bd = sorted([depth, ydepth])
    yratio = bd / sd
    if yratio > 2:
        df = df.query('chr != "chrY"')

    # 각 염색체별 repeat 영역에 대해 flanking 영역 계산
    for chrom, intervals in repeat_censat_data.items():
        chrom_length = chr_len.get(chrom)
        if chrom_length is None:
            continue
        chrom_df = df[df['chr'] == chrom]
        if chrom_df.empty:
            continue
        
        for rep in intervals:
            rep_start_0, rep_end_0 = rep  # 0-indexed 좌표

            results.append({
                'chr': chrom,
                'repeat_start_0': rep_start_0,
                'repeat_end_0': rep_end_0,
            })

    meandepth = np.median(df['meandepth'])

    def _node_overlaps_region(node, region_start, region_end):
        return max(int(node[CHR_STR]), region_start) <= min(int(node[CHR_END]), region_end)

    def _normalized_pair_endpoint_dir(pair, endpoint_order):
        ctg_dir = contig_data[pair[endpoint_order]][CTG_DIR]
        if endpoint_order == 0:
            return ctg_dir
        return _flip_ctg_dir(ctg_dir)

    def _expected_intervention_dir(side, right_high):
        # For a right-high jump, chr1:31M+ as pair[1] and chr1:33M- as pair[0]
        # normalize to '-' on the right flank, while chr1:13M+ normalizes to '+'
        # on the left flank and can explain the left-side drop. Reverse the roles
        # for a left-high jump.
        if right_high:
            return '+' if side == 'left' else '-'
        return '-' if side == 'left' else '+'

    def _has_depth_compatible_raw_nclose(chrom, rep_start_0, rep_end_0, flank_bp, right_high):
        if raw_nclose_nodes is None or contig_data is None:
            return False

        left_start = max(1, (rep_start_0 + 1) - flank_bp)
        left_end = rep_start_0
        right_start = rep_end_0 + 1
        right_end = min(chr_len[chrom], rep_end_0 + flank_bp)

        for pair_list in raw_nclose_nodes.values():
            for pair in pair_list:
                for endpoint_order, node_idx in enumerate(pair):
                    node = contig_data[node_idx]
                    if node[CHR_NAM] != chrom:
                        continue

                    side = None
                    if left_end >= left_start and _node_overlaps_region(node, left_start, left_end):
                        side = 'left'
                    elif right_end >= right_start and _node_overlaps_region(node, right_start, right_end):
                        side = 'right'
                    if side is None:
                        continue

                    norm_dir = _normalized_pair_endpoint_dir(pair, endpoint_order)
                    if norm_dir == _expected_intervention_dir(side, right_high):
                        return True
        return False

    depth_diff_data = dict()
    depth_dir_data = dict()
    relaxed_depth_chroms = set()
    for chrom, intervals in repeat_censat_data.items():
        chrom_length = chr_len.get(chrom)
        if chrom_length is None:
            continue
        chrom_df = df[df['chr'] == chrom]
        if chrom_df.empty:
            continue
        
        rep_start_0, rep_end_0 = intervals[0]  # 0-indexed 좌표

        if rep_start_0 == 0 or rep_end_0 == chrom_length - 1:
            continue

        tmp_depth_data = []
        for FLANK_SIZE_BP in [1*M, 5*M, 10*M]:
            # 좌측 flanking: repeat의 1-indexed 시작은 rep_start_0 + 1
            if rep_start_0 > 0:
                left_flank_end = rep_start_0  # repeat 시작 전 마지막 base (1-indexed)
                left_flank_start = max(1, (rep_start_0 + 1) - FLANK_SIZE_BP)

                if left_flank_end < left_flank_start or (left_flank_end - left_flank_start + 1) < FLANK_SIZE_BP:
                    left_flank_start = None
                    left_flank_end = None
                    left_weighted = None
                else:
                    left_weighted = weighted_avg_meandepth(chrom_df, left_flank_start, left_flank_start + FLANK_SIZE_BP)
            else:
                left_flank_start = None
                left_flank_end = None
                left_weighted = None

            # 우측 flanking: repeat의 끝 이후 첫 base부터
            if rep_end_0 < chrom_length:
                right_flank_start = rep_end_0 + 1
                right_flank_end = min(chrom_length, rep_end_0 + FLANK_SIZE_BP)

                if right_flank_end < right_flank_start or (right_flank_end - right_flank_start + 1) < FLANK_SIZE_BP:
                    right_flank_start = None
                    right_flank_end = None
                    right_weighted = None
                else:
                    right_weighted = weighted_avg_meandepth(chrom_df, right_flank_end - FLANK_SIZE_BP, right_flank_end)
            else:
                right_flank_start = None
                right_flank_end = None
                right_weighted = None

            if not (right_weighted is None or left_weighted is None):
                tmp_depth_data.append(right_weighted - left_weighted)
            else:
                break

        if not tmp_depth_data or tmp_depth_data[0] == 0:
            continue

        first_sign = tmp_depth_data[0] > 0

        rules = [(1, 0.20), (2, 0.15), (3, 0.10)]

        for k, thr in rules:
            if len(tmp_depth_data) < k:
                continue

            head = tmp_depth_data[:k]
            if any(x == 0 for x in head):
                continue

            same_sign = all((x > 0) == first_sign for x in head)
            strong_enough = all(abs(x) >= abs(thr * meandepth) for x in head)
            any_strong = any(abs(x) >= abs(thr * meandepth) for x in head)
            localized = all(
                abs(x - tmp_depth_data[0]) <= abs(CENSAT_OUT_DIFF_RATIO * meandepth)
                for x in head[1:]
            )

            if same_sign and strong_enough and localized:
                depth_diff_data[chrom] = abs(tmp_depth_data[0])
                depth_dir_data[chrom] = first_sign
                break
            elif raw_nclose_nodes is not None and contig_data is not None and k > 1 and same_sign and any_strong and localized:
                # A compatible breakpoint outside the inner 5 Mb can still change
                # the 5 Mb average by altering the segment toward the censat edge.
                nclose_flank_bp = 10*M
                has_nclose_intervention = _has_depth_compatible_raw_nclose(
                    chrom,
                    rep_start_0,
                    rep_end_0,
                    nclose_flank_bp,
                    first_sign,
                )
                if not has_nclose_intervention:
                    depth_diff_data[chrom] = max(
                        abs(x) for x in head if abs(x) >= abs(thr * meandepth)
                    )
                    depth_dir_data[chrom] = first_sign
                    relaxed_depth_chroms.add(chrom)
                    break


    result_df = pd.DataFrame(results)

    cen_fragment_meta = {}
    for chrom, diff in depth_diff_data.items():
        intervals = repeat_censat_data[chrom]
        rep_start_0, rep_end_0 = intervals[0]
        mid_censat = (rep_start_0 + rep_end_0) // 2
        
        cen_fragment_meta[chrom] = {
            "dir": depth_dir_data[chrom],
            "mid": mid_censat,
            "depth_diff": diff,
            "chr_len": chr_len[chrom],
        }

    # censat(centromere) breakend 로 인식된 염색체 목록을 로그로 출력
    detected_summary = ', '.join(sorted(cen_fragment_meta, key=chr2int))
    log_prefix = f'{log_context} ' if log_context else ''
    logging.info(
        f'{log_prefix}Censat breakend chromosome: '
        f'{detected_summary if cen_fragment_meta else "(none)"}'
    )
    if relaxed_depth_chroms:
        relaxed_summary = ', '.join(sorted(relaxed_depth_chroms, key=chr2int))
        logging.info(
            f'{log_prefix}Censat breakend chromosome relaxed by raw nclose direction: '
            f'{relaxed_summary}'
        )

    def find_min_error_partition(depth_diff_data):
        """
        Partition keys into groups of size 1, 2, or 3 to minimize total error,
        while ensuring that in each group the dict key is the 'minuend' (largest
        by the rule below) and the value is the list of subtrahends.

        For size-2 groups: key is the one with the larger value (tie -> lex key).
        For size-3 groups: choose i that minimizes |v_i - (v_j + v_k)|; that i becomes key.
        """
        # Deterministic ordering
        all_keys = tuple(sorted(depth_diff_data.keys()))
        memo = {}

        def best_pair_mapping(k1, k2):
            v1, v2 = depth_diff_data[k1], depth_diff_data[k2]
            # key should be the one with larger value; tie -> lexicographic key
            if (v1 > v2) or (v1 == v2 and k1 < k2):
                return {k1: [k2]}, abs(v1 - v2)
            else:
                return {k2: [k1]}, abs(v1 - v2)

        def best_triple_mapping(a, b, c):
            v = {a: depth_diff_data[a], b: depth_diff_data[b], c: depth_diff_data[c]}
            # candidates: which one is the minuend
            candidates = []
            # a as key
            candidates.append( (a, abs(v[a] - (v[b] + v[c]))) )
            # b as key
            candidates.append( (b, abs(v[b] - (v[a] + v[c]))) )
            # c as key
            candidates.append( (c, abs(v[c] - (v[a] + v[b]))) )

            # pick by minimal error; tie-break by larger value, then lex key
            min_err = min(err for _, err in candidates)
            tied = [k for k, err in candidates if err == min_err]
            if len(tied) == 1:
                key = tied[0]
            else:
                # larger value first
                max_val = max(v[k] for k in tied)
                tied2 = [k for k in tied if v[k] == max_val]
                key = min(tied2) if len(tied2) > 1 else tied2[0]

            others = [x for x in (a, b, c) if x != key]
            return {key: others}, min_err

        def recurse(remaining):
            if not remaining:
                return {}, 0
            if remaining in memo:
                return memo[remaining]

            # we always pop the first, but mapping keys will be set by rules above
            first, *rest = remaining
            best_mapping = None
            best_err = float('inf')

            # Case 1: single
            err1 = abs(depth_diff_data[first])
            mapping1, err_sum1 = recurse(tuple(rest))
            total1 = err1 + err_sum1
            best_mapping = {first: []}
            best_mapping.update(mapping1)
            best_err = total1

            # Case 2: pairs
            for i, k2 in enumerate(rest):
                pair_map, pair_err = best_pair_mapping(first, k2)
                rem2 = rest[:i] + rest[i+1:]
                m2, e2 = recurse(tuple(rem2))
                total2 = pair_err + e2
                if total2 < best_err:
                    best_err = total2
                    # merge maps; ensure no conflicting keys
                    best_mapping = {}
                    best_mapping.update(pair_map)
                    best_mapping.update(m2)

            # Case 3: triples
            n = len(rest)
            for i in range(n):
                for j in range(i+1, n):
                    k2, k3 = rest[i], rest[j]
                    triple_map, triple_err = best_triple_mapping(first, k2, k3)
                    rem3 = tuple(rest[:i] + rest[i+1:j] + rest[j+1:])
                    m3, e3 = recurse(rem3)
                    total3 = triple_err + e3
                    if total3 < best_err:
                        best_err = total3
                        best_mapping = {}
                        best_mapping.update(triple_map)
                        best_mapping.update(m3)

            memo[remaining] = (best_mapping, best_err)
            return memo[remaining]

        return recurse(all_keys)

    cnt = 0
    vtg_list = []
    prefix = "virtual_censat_contig"

    # 가상 BND contig 생성 비활성화: centromere depth 차이는 22번 NNLS의
    # fragment chromosome column 으로 흡수한다. 페어링/노드 생성 코드는 보존.
    partition_dict = {}
    # partition_dict, error = find_min_error_partition(depth_diff_data)

    for k, v in partition_dict.items():
        connecting_pair = []
        if len(v) == 1:
            connecting_pair.append((k, v[0]))
        elif len(v) == 2:
            connecting_pair.append((k, v[0]))
            connecting_pair.append((k, v[1]))
        
        for chrom1, chrom2 in connecting_pair:
            row1 = result_df[result_df['chr'] == chrom1].iloc[0]
            row2 = result_df[result_df['chr'] == chrom2].iloc[0]
            if depth_dir_data[chrom1] and depth_dir_data[chrom2]: #right right
                cnt+=1
                N = int(CENSAT_COMPRESSABLE_THRESHOLD//2)
                mid_row1_censat = int(row1.repeat_start_0 + row1.repeat_end_0)//2
                mid_row2_censat = int(row2.repeat_start_0 + row2.repeat_end_0)//2

                temp_node1 = [f'{prefix}_{cnt}', N, 0, N//2, '-', row1.chr, 
                            chr_len[row1.chr], mid_row1_censat - N//2, mid_row1_censat, 
                            60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '-', row1.chr, f"2.{(cnt-1)*2}"]
                
                temp_node2 = [f'{prefix}_{cnt}', N, N//2, N, '+', row2.chr,
                            chr_len[row2.chr], mid_row2_censat, mid_row2_censat + N//2, 
                            60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '-', row1.chr, f"2.{(cnt-1)*2+1}"]
                
                vtg_list.append(temp_node1)
                vtg_list.append(temp_node2)
            elif not depth_dir_data[chrom1] and not depth_dir_data[chrom2]: #left left
                cnt+=1
                N = int(CENSAT_COMPRESSABLE_THRESHOLD//2)
                mid_row1_censat = int(row1.repeat_start_0 + row1.repeat_end_0)//2
                mid_row2_censat = int(row2.repeat_start_0 + row2.repeat_end_0)//2

                temp_node1 = [f'{prefix}_{cnt}', N, 0, N//2, '+', row1.chr, 
                            chr_len[row1.chr], mid_row1_censat - N//2, mid_row1_censat, 
                            60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '+', row1.chr, f"2.{(cnt-1)*2}"]
                
                temp_node2 = [f'{prefix}_{cnt}', N, N//2, N, '-', row2.chr,
                            chr_len[row2.chr], mid_row2_censat, mid_row2_censat + N//2, 
                            60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '+', row1.chr, f"2.{(cnt-1)*2+1}"]
                
                vtg_list.append(temp_node1)
                vtg_list.append(temp_node2)
            elif depth_dir_data[chrom1] and not depth_dir_data[chrom2]: #right left
                cnt+=1
                N = int(CENSAT_COMPRESSABLE_THRESHOLD//2)
                mid_row1_censat = int(row1.repeat_start_0 + row1.repeat_end_0)//2
                mid_row2_censat = int(row2.repeat_start_0 + row2.repeat_end_0)//2

                temp_node1 = [f'{prefix}_{cnt}', N, 0, N//2, '-', row1.chr, 
                        chr_len[row1.chr], mid_row1_censat - N//2, mid_row1_censat, 
                        60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '-', row1.chr, f"2.{(cnt-1)*2}"]
            
                temp_node2 = [f'{prefix}_{cnt}', N, N//2, N, '-', row2.chr,
                        chr_len[row2.chr], mid_row2_censat, mid_row2_censat + N//2, 
                        60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '-', row1.chr, f"2.{(cnt-1)*2+1}"]
                
                vtg_list.append(temp_node1)
                vtg_list.append(temp_node2)
            else: #left right
                cnt+=1
                N = int(CENSAT_COMPRESSABLE_THRESHOLD//2)
                mid_row1_censat = int(row1.repeat_start_0 + row1.repeat_end_0)//2
                mid_row2_censat = int(row2.repeat_start_0 + row2.repeat_end_0)//2

                temp_node1 = [f'{prefix}_{cnt}', N, 0, N//2, '+', row1.chr, 
                        chr_len[row1.chr], mid_row1_censat - N//2, mid_row1_censat, 
                        60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '+', row1.chr, f"2.{(cnt-1)*2}"]
            
                temp_node2 = [f'{prefix}_{cnt}', N, N//2, N, '+', row2.chr,
                        chr_len[row2.chr], mid_row2_censat, mid_row2_censat + N//2, 
                        60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '+', row1.chr, f"2.{(cnt-1)*2+1}"]
                
                vtg_list.append(temp_node1)
                vtg_list.append(temp_node2)

    return vtg_list, cen_fragment_meta

def break_double_telomere_contig(contig_data : list, telo_connected_set : set):
    s = 0
    vtg_list = []
    idx = 0
    count = 0
    tot_cur_ctg_cnt = 0
    contig_data_size = len(contig_data)
    while s<contig_data_size:
        e = contig_data[s][CTG_ENDND]
        cnt = 0
        for i in range(s+1, e):
            if contig_data[i][CTG_TELDIR][0] in ('f', 'b')\
            and s not in telo_connected_set \
            and e not in telo_connected_set:
                cnt+=1
        if cnt >= 2:
            st = s
            ed = e
            while st <= e and contig_data[st][CTG_TELDIR] == '0':
                st+=1
            while ed >= s and contig_data[ed][CTG_TELDIR] == '0':
                ed-=1
            if s <= st:
                tot_cur_ctg_cnt+=1
            for i in range(s, st+1):
                temp_list = copy.deepcopy(contig_data[i])
                temp_list[CTG_NAM] = f'telomere_middle_cut_contig_{tot_cur_ctg_cnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + st - s
                vtg_list.append(temp_list)
                count+=1
            idx += st-s+1
            if ed<=e:
                tot_cur_ctg_cnt+=1
            for i in range(ed, e+1):
                temp_list = copy.deepcopy(contig_data[i])
                temp_list[CTG_NAM] = f'telomere_middle_cut_contig_{tot_cur_ctg_cnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + e - ed
                vtg_list.append(temp_list)
                count+=1
        s = e+1
    return vtg_list

def delete_contig(contig_data : list, to_delete_contig_set : set) -> list:
    new_contig_data = []
    delete_count = 0
    for chunk in contig_data:
        if chunk[CTG_NAM] in to_delete_contig_set:
            delete_count += 1
        else:
            chunk[CTG_STRND] -= delete_count
            chunk[CTG_ENDND] -= delete_count
            new_contig_data.append(chunk)
    
    return new_contig_data

def break_type34_contig(contig_data : list):
    s = 0
    vtg_list = []
    broken_contig_set = set()
    
    while s < len(contig_data):
        e = contig_data[s][CTG_ENDND]
        curr_contig_name = contig_data[s][CTG_NAM]
        broken_chk = False  # reset per contig: avoid UnboundLocalError on non-type3 first contig and stale carry-over
        if contig_data[s][CTG_TYP] == 3:
            cnt = 0
            for i in range(s, e):
                front_chr = (contig_data[i][CHR_NAM], contig_data[i][CTG_DIR])
                back_chr = (contig_data[i+1][CHR_NAM], contig_data[i+1][CTG_DIR])
                ref_ratio, chukji = calculate_single_contig_ref_ratio(contig_data[i:i+2])
                if contig_data[i][CTG_RPTCHR] != '0' or contig_data[i+1][CTG_RPTCHR] != '0':
                    bnd_bound = RPT_BND_CONTIG_BOUND
                else:
                    bnd_bound = BND_CONTIG_BOUND
                if front_chr == back_chr and abs(1-ref_ratio) > bnd_bound and abs(chukji) > TYPE34_BREAK_CHUKJI_LIMIT:
                    broken_chk = True
                    cnt += 1
                    vtg = copy.deepcopy(contig_data[i:i+2])
                    for j in vtg:
                        j[CTG_NAM] = f'{curr_contig_name}_split34_{cnt}'
                        vtg_list.append(j)
        if broken_chk:
            broken_contig_set.add(curr_contig_name)
        s = e+1
    return vtg_list, broken_contig_set

def pass_pipeline(pre_contig_data, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, telo_ppc_passed, chr_len):
    if not telo_ppc_passed:
        if len(pre_contig_data)==0:
            return []
        contig_data = []

        for i in pre_contig_data:
            contig_data.append(i[:9] + [i[CTG_MAPQ], i[CTG_GLOBALIDX]])

        node_label = label_node(contig_data, telo_dict)

        repeat_label = label_repeat_node(contig_data, repeat_data, chr_len)

        telo_preprocessed_contig, report_case, telo_connect_info = preprocess_telo(contig_data, node_label)

        new_contig_data = []
        telcon_set = set()
        idxcnt = 0
        for i in telo_preprocessed_contig:
            temp_list = contig_data[i]
            if i in telo_connect_info:
                temp_list.append(telo_connect_info[i])
                telcon_set.add(idxcnt)
            else:
                temp_list.append("0")
            new_contig_data.append(temp_list)
            idxcnt+=1

        
        new_node_repeat_label = label_repeat_node(new_contig_data, repeat_data, chr_len)
        new_node_repeat_censat_label = label_repeat_node(new_contig_data, repeat_censat_data, chr_len)
        new_node_telo_label = label_node(new_contig_data, telo_dict)

        ref_qry_ratio = calc_ratio(new_contig_data)
        preprocess_result, \
        preprocess_contig_type, \
        preprocess_terminal_nodes, \
        len_counter = pipeline_preprocess_contig(new_contig_data, new_node_telo_label, ref_qry_ratio, new_node_repeat_label, telcon_set)
        preprocess_result = set(preprocess_result)
        new_contig_data = new_contig_data[0:-1]
        contig_data_size = len(new_contig_data)

        telo_ppc_contig = []
        cnt = 0
        for i in range(contig_data_size):
            if new_contig_data[i][CTG_NAM] in preprocess_result:
                temp_list = new_contig_data[i][:10]
                temp_list.append(preprocess_contig_type[new_contig_data[i][CTG_NAM]])
                temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][0])
                temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][1])
                temp_list.append(new_node_telo_label[i][0])
                temp_list.append(new_node_telo_label[i][1])
                temp_list.append(new_contig_data[i][11])
                temp_list.append(new_node_repeat_label[i][0])
                temp_list.append(new_node_repeat_label[i][1])
                temp_list.append(new_node_repeat_censat_label[i][1])
                temp_list.append(new_contig_data[i][10])
                telo_ppc_contig.append(temp_list)
                cnt+=1
        mainflow_dict = find_mainflow(telo_ppc_contig)
        telo_ppc_size = len(telo_ppc_contig)
        
        for i in range(telo_ppc_size):
            max_chr = mainflow_dict[telo_ppc_contig[i][CTG_NAM]]
            temp = [max_chr[0], max_chr[1], telo_ppc_contig[i][-1]]
            telo_ppc_contig[i] = telo_ppc_contig[i][:-1]
            telo_ppc_contig[i] += temp
    else:
        if len(pre_contig_data)==0:
            return []
        telo_ppc_contig = pre_contig_data
        mainflow_dict = find_mainflow(telo_ppc_contig)
        telo_ppc_size = len(telo_ppc_contig)
        
        for i in range(telo_ppc_size):
            max_chr = mainflow_dict[telo_ppc_contig[i][CTG_NAM]]
            telo_ppc_contig[i][CTG_MAINFLOWDIR] = max_chr[0]
            telo_ppc_contig[i][CTG_MAINFLOWCHR] = max_chr[1]

    

    final_contig = preprocess_repeat(telo_ppc_contig)

    final_contig_repeat_label = label_repeat_node(final_contig, repeat_data, chr_len)

    final_telo_node_label = label_node(final_contig, telo_dict)

    final_ref_qry_ratio = calc_ratio(final_contig)

    final_telo_connect = set()

    for v, i in enumerate(final_contig):
        try:
            if i[CTG_TELCON] != '0':
                final_telo_connect.add(v)
        except:
            pass

    final_using_contig, final_ctg_typ, final_preprocess_terminal_nodes, _ = pipeline_preprocess_contig(final_contig, final_telo_node_label, final_ref_qry_ratio, final_contig_repeat_label, final_telo_connect)

    final_contig = final_contig[:-1]

    real_final_contig = []

    final_using_contig = set(final_using_contig)
    for i in range(0, len(final_contig)):
        if final_contig[i][CTG_NAM] in final_using_contig:
            final_contig[i][CTG_TYP] = final_ctg_typ[final_contig[i][CTG_NAM]]
            final_contig[i][CTG_STRND] = final_preprocess_terminal_nodes[final_contig[i][CTG_NAM]][0]
            final_contig[i][CTG_ENDND] = final_preprocess_terminal_nodes[final_contig[i][CTG_NAM]][1]
            real_final_contig.append(final_contig[i])

    return real_final_contig



def extract_nclose_node(contig_data : list, bnd_contig : set, repeat_contig_name : set, \
                        censat_contig_name : set, repeat_censat_data : dict, ALIGNED_PAF_LOC_LIST : list, ORIGINAL_PAF_LOC_LIST : list, \
                        telo_set : set, telo_contig : dict, chr_len : dict, asm2cov : dict) -> tuple:
    s = 0
    fake_bnd = dict()
    contig_data_size = len(contig_data)
    nclose_compress = defaultdict(list)
    nclose_compress_track = defaultdict(list)
    nclose_start_compress = defaultdict(lambda : defaultdict(list))
    nclose_end_compress = defaultdict(lambda : defaultdict(list))
    censat_nclose_compress = set()
    nclose_dict = defaultdict(list)
    all_nclose_compress = defaultdict(list)
    telo_name_set = set()
    nclose_coverage = Counter()
    for i in telo_set:
        telo_name_set.add(contig_data[i][CTG_NAM])

    # [RECIPROCAL FIX] nclose 방향 정규화: '작은 염색체 노드를 먼저 읽는' 프레임의 (low_dir, high_dir).
    # all_nclose 출력의 get_corr_dir 와 동일 규칙이라 reverse-complement 대칭을 자동 보장한다
    # (예: chr12+ -> chr15+ 와 chr15- -> chr12- 는 같은 값으로 정규화됨).
    # balanced reciprocal translocation 두 산물(예: der(12) vs der(15))은 같은 breakpoint를
    # 공유해 좌표만으로는 구분되지 않으므로(distance_checker 가 겹침을 0으로 봄), 압축 병합
    # 판정에 이 정규화 방향을 함께 비교해 한쪽이 통째로 버려지는 것을 막는다.
    def nclose_canon_dir(low_node, high_node):
        is_for = low_node < high_node
        return (get_corr_dir(is_for, contig_data[low_node][CTG_DIR]),
                get_corr_dir(is_for, contig_data[high_node][CTG_DIR]))

    div_repeat_paf_name = div_repeat_paf(ORIGINAL_PAF_LOC_LIST, ALIGNED_PAF_LOC_LIST, contig_data)

    while s<contig_data_size:
        e = contig_data[s][CTG_ENDND]
        contig_s = contig_data[s]
        contig_e = contig_data[e]
        st = s
        ed = e

        if contig_data[s][CTG_NAM] in bnd_contig:
            if contig_data[s][CTG_TYP] in (1, 2):
                st_chr = [contig_data[s][CTG_DIR], contig_data[s][CHR_NAM]]
                ed_chr = [contig_data[e][CTG_DIR], contig_data[e][CHR_NAM]]
                # simple_ctg_alt: 1-chunk noise(다른 chr/dir 단일 chunk)는 무시하고 더 깊이 trim
                is_simple_ctg_alt = contig_data[s][CTG_NAM].startswith("simple_ctg_alt_")
                MAX_GAP = 1 if is_simple_ctg_alt else 0

                # Keep the two terminal trims from crossing on ABAB-like simple_ctg_alt paths.
                # Otherwise the 1-chunk gap tolerance can emit a reverse nclose pair.
                # forward trim
                gap = 0
                last_match_st = st  # st == s, st_chr 매칭 by construction
                st_cur = st + 1
                while st_cur < e:
                    if [contig_data[st_cur][CTG_DIR], contig_data[st_cur][CHR_NAM]] == st_chr:
                        last_match_st = st_cur
                        gap = 0
                    else:
                        gap += 1
                        if gap > MAX_GAP:
                            break
                    cut_ratio, _ = calculate_single_contig_ref_ratio(contig_data[s:last_match_st+1])
                    if abs(cut_ratio-1) > BND_CONTIG_BOUND:
                        break
                    st_cur += 1
                st = last_match_st

                # backward trim
                gap = 0
                last_match_ed = ed  # ed == e, ed_chr 매칭 by construction
                ed_cur = ed - 1
                while ed_cur > st:
                    if [contig_data[ed_cur][CTG_DIR], contig_data[ed_cur][CHR_NAM]] == ed_chr:
                        last_match_ed = ed_cur
                        gap = 0
                    else:
                        gap += 1
                        if gap > MAX_GAP:
                            break
                    cut_ratio, _ = calculate_single_contig_ref_ratio(contig_data[last_match_ed:e+1])
                    if abs(cut_ratio-1) > BND_CONTIG_BOUND:
                        break
                    ed_cur -= 1
                ed = last_match_ed

                if st + 1 < ed:
                    flags = [contig_data[ed-1][CTG_DIR], contig_data[ed-1][CHR_NAM]] == st_chr
                    flage = [contig_data[st+1][CTG_DIR], contig_data[st+1][CHR_NAM]] == ed_chr
                        
                    if flags and flage:
                        ratio_s, len_s = calculate_single_contig_ref_ratio(contig_data[s:ed])
                        ratio_e, len_e = calculate_single_contig_ref_ratio(contig_data[st+1:e+1])
                        if abs(ratio_s-1) < BND_CONTIG_BOUND and abs(ratio_e-1) < BND_CONTIG_BOUND:
                            if len_s > len_e:
                                fake_bnd[contig_s[CTG_NAM]] = (s, ed-1)
                                st = ed-1
                            else:
                                fake_bnd[contig_e[CTG_NAM]] = (st+1, e)
                                ed = st+1
                        elif abs(ratio_s-1) < BND_CONTIG_BOUND:
                            fake_bnd[contig_s[CTG_NAM]] = (s, ed-1)
                            st = ed-1
                        elif abs(ratio_e-1) < BND_CONTIG_BOUND:
                            fake_bnd[contig_e[CTG_NAM]] = (st+1, e)
                            ed = st+1
                        else:
                            if len_s > len_e:
                                st = ed-1
                            else:
                                ed = st+1
                    elif flags:
                        ratio_s, len_s = calculate_single_contig_ref_ratio(contig_data[s:ed])
                        if abs(ratio_s-1) < BND_CONTIG_BOUND:
                            fake_bnd[contig_data[s][CTG_NAM]] = (s, ed-1)
                        st = ed-1
                    elif flage:
                        ratio_e, len_e = calculate_single_contig_ref_ratio(contig_data[st+1:e+1])
                        if abs(ratio_e-1) < BND_CONTIG_BOUND:
                            fake_bnd[contig_data[s][CTG_NAM]] = (st+1, e)
                        ed = st+1

                all_nclose_compress[(st_chr[1], ed_chr[1])].append((st, ed))
                
                if contig_data[s][CTG_NAM] in div_repeat_paf_name \
                   and contig_data[s][CTG_NAM] in repeat_contig_name \
                   and st not in telo_set \
                   and ed not in telo_set \
                   and contig_data[s][CTG_NAM] not in telo_name_set:
                    s = e+1
                    continue
                flag = True
                # if st_chr[0] == ed_chr[0]:
                #     st_chr[0] = ed_chr[0] = '+'
                # else:
                #     if chr2int(st_chr[1]) < chr2int(ed_chr[1]):
                #         st_chr[0] = '+'
                #         ed_chr[0] = '+' #
                #     elif chr2int(st_chr[1]) > chr2int(ed_chr[1]):
                #         st_chr[0] = '+' #
                #         ed_chr[0] = '+'
                #     else:
                #         if (contig_data[st][CHR_STR], contig_data[st][CHR_END]) < (contig_data[ed][CHR_STR], contig_data[ed][CHR_END]):
                #             st_chr[0] = '+'
                #             ed_chr[0] = '+' #
                #         else:
                #             st_chr[0] = '+' #
                #             ed_chr[0] = '+'

                # if contig_data[st][CTG_CENSAT] != '0' and contig_data[ed][CTG_CENSAT] != '0' and contig_data[st][CTG_NAM] in div_repeat_paf_name:
                #     s = e+1
                #     continue
                st_chr[0] = ed_chr[0] = '='  

                st_chr = tuple(st_chr)
                ed_chr = tuple(ed_chr)

                max_chr = (contig_data[st][CTG_MAINFLOWDIR], contig_data[st][CTG_MAINFLOWCHR])
                combo = 0
                ref_combo = 0
                maxcombo=0
                max_ref_combo = 0
                st_idx = 0
                ed_idx = 0
                for i in range(st, ed+1):
                    if (contig_data[i][CTG_DIR], contig_data[i][CHR_NAM]) == max_chr:
                        combo+=1
                        ref_combo += contig_data[i][CHR_END] - contig_data[i][CHR_STR]
                    else:
                        combo=0
                        ref_combo = 0
                    if ref_combo>max_ref_combo:
                        maxcombo = combo
                        max_ref_combo = ref_combo
                        st_idx = i-maxcombo+1
                        ed_idx = i
                nclose_list = []
                # simple_ctg_alt contigs: skip mainflow split — keep exactly one nclose pair per contig
                is_simple_ctg_alt = contig_data[s][CTG_NAM].startswith("simple_ctg_alt_")
                if not is_simple_ctg_alt \
                    and max_chr not in ((contig_data[s][CTG_DIR], contig_data[s][CHR_NAM]), (contig_data[e][CTG_DIR], contig_data[e][CHR_NAM])) \
                    and st_idx != ed_idx:
                    st_chr = [contig_data[st][CTG_DIR], contig_data[st][CHR_NAM]]
                    ed_chr = [contig_data[st_idx][CTG_DIR], contig_data[st_idx][CHR_NAM]]
                    st_chr[0] = ed_chr[0] = '='
                    st_chr = tuple(st_chr)
                    ed_chr = tuple(ed_chr)
                    nclose_list.append([st, st_chr, st_idx, ed_chr])
                    st_chr = [contig_data[ed_idx][CTG_DIR], contig_data[ed_idx][CHR_NAM]]
                    ed_chr = [contig_data[ed][CTG_DIR], contig_data[ed][CHR_NAM]]
                    st_chr[0] = ed_chr[0] = '='
                    st_chr = tuple(st_chr)
                    ed_chr = tuple(ed_chr)
                    nclose_list.append([ed_idx, st_chr, ed, ed_chr])
                else:
                    nclose_list.append([st, st_chr, ed, ed_chr])

                for nclose in nclose_list:
                    flag = True  # reset per nclose: novel until proven duplicate (mainflow split can hold 2 ncloses)
                    st = nclose[0]
                    ed = nclose[2]
                    # Pass if 
                    # 1. terminal node is rin
                    # 2. AND both terminal have similar read depth
                    """
                    if contig_data[st][CHR_NAM] != contig_data[ed][CHR_NAM]:
                        if contig_data[st][CTG_RPTCASE] == 'rin' or contig_data[ed][CTG_RPTCASE] == 'rin':
                            st_data = (contig_data[st][CHR_NAM], contig_data[st][CHR_END])
                            ed_data = (contig_data[ed][CHR_NAM], contig_data[ed][CHR_STR])
                            if check_depth_near_bnd(st_data, ed_data):
                                continue
                    """
                    nclose_front_const = 0
                    nclose_back_const = 0
                    for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                        if inclusive_checker_tuple(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                            nclose_front_const = 1
                            if censat_ref_range[0] < K or censat_ref_range[1] > chr_len[contig_data[st][CHR_NAM]]-K:
                                nclose_front_const = 2
                    for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                        if inclusive_checker_tuple(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
                            nclose_back_const = 1
                            if censat_ref_range[0] < K or censat_ref_range[1] > chr_len[contig_data[ed][CHR_NAM]]-K:
                                nclose_back_const = 2
                    if nclose_front_const + nclose_back_const >= 3 \
                       and not (contig_data[st][CTG_MAPQ] == 60 and contig_data[ed][CTG_MAPQ] == 60):
                        continue
                    st_chr = nclose[1]
                    ed_chr = nclose[3]
                    upd_contig_name = contig_data[st][CTG_NAM]
                    cl_ind, c_ind = map(int, contig_data[st][CTG_GLOBALIDX].split('.'))
                    if cl_ind < 2:
                        cov_count_name = ori_ctg_name_data[cl_ind][c_ind]
                    else:
                        asm2cov[upd_contig_name] = -1
                        cov_count_name = upd_contig_name
                    is_curr_ctg_repeat = upd_contig_name in repeat_contig_name
                    if chr2int(st_chr[1]) < chr2int(ed_chr[1]): 
                        for i in nclose_compress[(st_chr, ed_chr)]:
                            if is_curr_ctg_repeat and i[0] in repeat_contig_name:
                                compress_limit = ALL_REPEAT_NCLOSE_COMPRESS_LIMIT
                            else:
                                compress_limit = NCLOSE_COMPRESS_LIMIT
                            dummy_list = [0,0,0,0,0,0,0,]
                            if distance_checker(contig_data[st], dummy_list+i[1]) < compress_limit \
                            and distance_checker(contig_data[ed], dummy_list + i[2]) < compress_limit \
                            and nclose_canon_dir(st, ed) == nclose_canon_dir(i[3], i[4]):
                                sorted_tar_tuple = tuple(sorted((i[3], i[4])))
                                nclose_coverage[sorted_tar_tuple] += asm2cov[cov_count_name]
                                nclose_compress_track[sorted_tar_tuple].append((st, ed))
                                flag = False
                                break
                        # passed first filtering
                        if flag:
                            censat_st_chr = [st_chr[0], 0]
                            censat_ed_chr = [ed_chr[0], 0]
                            if contig_s[CTG_CENSAT] != '0':
                                cnt = 0
                                for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                                    if inclusive_checker_tuple(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                                        censat_st_chr[1] = contig_data[st][CHR_NAM] + "." + str(cnt)
                                        break
                                    cnt+=1
                            if contig_e[CTG_CENSAT] != '0':
                                cnt = 0
                                for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                                    if inclusive_checker_tuple(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
                                        censat_ed_chr[1] = contig_data[ed][CHR_NAM] + "." + str(cnt)
                                        break
                                    cnt+=1
                            # at least one of nclose is not censat
                            if censat_st_chr[1] == 0 \
                            or censat_ed_chr[1] == 0 :
                                temp_list = [upd_contig_name, 
                                            [contig_data[st][CHR_STR], contig_data[st][CHR_END]], 
                                            [contig_data[ed][CHR_STR], contig_data[ed][CHR_END]], st, ed]
                                for i in nclose_compress[(st_chr, ed_chr)]:
                                    dummy_list = [0,0,0,0,0,0,0,]
                                    flag = False
                                    if distance_checker(contig_data[st], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                        flag = True
                                        nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                    if distance_checker(contig_data[ed], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                        flag = True
                                        nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                    if flag:
                                        break
                                nclose_compress[(st_chr, ed_chr)].append(temp_list)
                                nclose_coverage[(st, ed)] += asm2cov[cov_count_name]
                                nclose_dict[contig_data[s][CTG_NAM]].append((st, ed))
                            # both are censat
                            else:
                                key = (tuple(censat_st_chr), tuple(censat_ed_chr))
                                # if censat only one can survive
                                if key not in censat_nclose_compress:
                                    censat_nclose_compress.add(key)
                                    temp_list = [upd_contig_name, 
                                            [contig_data[st][CHR_STR], contig_data[st][CHR_END]], 
                                            [contig_data[ed][CHR_STR], contig_data[ed][CHR_END]], st, ed]
                                    for i in nclose_compress[(st_chr, ed_chr)]:
                                        dummy_list = [0,0,0,0,0,0,0,]
                                        flag = False
                                        if distance_checker(contig_data[st], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                            flag = True
                                            nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                        if distance_checker(contig_data[ed], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                            flag = True
                                            nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                        if flag:
                                            break
                                    nclose_compress[(st_chr, ed_chr)].append(temp_list)
                                    nclose_coverage[(st, ed)] += asm2cov[cov_count_name]
                                    nclose_dict[contig_data[s][CTG_NAM]].append((st, ed))
                    elif chr2int(st_chr[1]) > chr2int(ed_chr[1]):
                        for i in nclose_compress[(ed_chr, st_chr)]:
                            if is_curr_ctg_repeat and i[0] in repeat_contig_name:
                                compress_limit = ALL_REPEAT_NCLOSE_COMPRESS_LIMIT
                            else:
                                compress_limit = NCLOSE_COMPRESS_LIMIT
                            dummy_list = [0,0,0,0,0,0,0,]
                            if distance_checker(contig_data[st], dummy_list+i[2]) < compress_limit \
                            and distance_checker(contig_data[ed], dummy_list + i[1]) < compress_limit \
                            and nclose_canon_dir(ed, st) == nclose_canon_dir(i[3], i[4]):
                                sorted_tar_tuple = tuple(sorted((i[3], i[4])))
                                nclose_coverage[sorted_tar_tuple] += asm2cov[cov_count_name]
                                nclose_compress_track[sorted_tar_tuple].append((st, ed))
                                flag = False
                                break
                        if flag:
                            censat_st_chr = [st_chr[0], 0]
                            censat_ed_chr = [ed_chr[0], 0]
                            if contig_s[CTG_CENSAT] != '0':
                                cnt = 0
                                for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                                    if inclusive_checker_tuple(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                                        censat_st_chr[1] = contig_data[st][CHR_NAM] + "." + str(cnt)
                                        break
                                    cnt+=1
                            if contig_e[CTG_CENSAT] != '0':
                                cnt = 0
                                for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                                    if inclusive_checker_tuple(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
                                        censat_ed_chr[1] = contig_data[ed][CHR_NAM] + "." + str(cnt)
                                        break
                                    cnt+=1
                            if censat_st_chr[1] == 0 \
                            or censat_ed_chr[1] == 0 :
                                temp_list = [upd_contig_name,
                                        [contig_data[ed][CHR_STR], contig_data[ed][CHR_END]], 
                                        [contig_data[st][CHR_STR], contig_data[st][CHR_END]], ed, st]
                                for i in nclose_compress[(ed_chr, st_chr)]:
                                    dummy_list = [0,0,0,0,0,0,0,]
                                    flag = False
                                    if distance_checker(contig_data[ed], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                        nclose_start_compress[(ed_chr, st_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                        flag = True
                                    if distance_checker(contig_data[st], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                        nclose_end_compress[(ed_chr, st_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                        flag = True
                                    if flag:
                                        break
                                nclose_coverage[(st, ed)] += asm2cov[cov_count_name]
                                nclose_compress[(ed_chr, st_chr)].append(temp_list)
                                nclose_dict[contig_data[s][CTG_NAM]].append((st, ed))
                            else:
                                key = (tuple(censat_st_chr), tuple(censat_ed_chr))
                                if key not in censat_nclose_compress:
                                    censat_nclose_compress.add(key)
                                    temp_list = [upd_contig_name,
                                        [contig_data[ed][CHR_STR], contig_data[ed][CHR_END]], 
                                        [contig_data[st][CHR_STR], contig_data[st][CHR_END]], ed, st]
                                    for i in nclose_compress[(ed_chr, st_chr)]:
                                        dummy_list = [0,0,0,0,0,0,0,]
                                        flag = False
                                        if distance_checker(contig_data[ed], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                            nclose_start_compress[(ed_chr, st_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                            flag = True
                                        if distance_checker(contig_data[st], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                            nclose_end_compress[(ed_chr, st_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                            flag = True
                                        if flag:
                                            break
                                    nclose_coverage[(st, ed)] += asm2cov[cov_count_name]
                                    nclose_compress[(ed_chr, st_chr)].append(temp_list)
                                    nclose_dict[contig_data[s][CTG_NAM]].append((st, ed))
                    else:
                        if (contig_data[st][CHR_STR], contig_data[st][CHR_END]) <= (contig_data[ed][CHR_STR], contig_data[ed][CHR_END]):
                            for i in nclose_compress[(st_chr, ed_chr)]:
                                if is_curr_ctg_repeat and i[0] in repeat_contig_name:
                                    compress_limit = ALL_REPEAT_NCLOSE_COMPRESS_LIMIT
                                else:
                                    compress_limit = NCLOSE_COMPRESS_LIMIT
                                dummy_list = [0,0,0,0,0,0,0,]
                                if distance_checker(contig_data[st], dummy_list+i[1]) < compress_limit \
                                and distance_checker(contig_data[ed], dummy_list + i[2]) < compress_limit:
                                    sorted_tar_tuple = tuple(sorted((i[3], i[4])))
                                    nclose_coverage[sorted_tar_tuple] += asm2cov[cov_count_name]
                                    nclose_compress_track[sorted_tar_tuple].append((st, ed))
                                    flag = False
                                    break
                            if flag:
                                censat_st_chr = [st_chr[0], 0]
                                censat_ed_chr = [ed_chr[0], 0]
                                if contig_s[CTG_CENSAT] != '0':
                                    cnt = 0
                                    for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                                        if inclusive_checker_tuple(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                                            censat_st_chr[1] = contig_data[st][CHR_NAM] + "." + str(cnt)
                                            break
                                        cnt+=1
                                if contig_e[CTG_CENSAT] != '0':
                                    cnt = 0
                                    for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                                        if inclusive_checker_tuple(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
                                            censat_ed_chr[1] = contig_data[ed][CHR_NAM] + "." + str(cnt)
                                            break
                                        cnt+=1
                                if censat_st_chr[1] == 0 \
                                or censat_ed_chr[1] == 0 :
                                    temp_list = [upd_contig_name, 
                                            [contig_data[st][CHR_STR], contig_data[st][CHR_END]], 
                                            [contig_data[ed][CHR_STR], contig_data[ed][CHR_END]], st, ed]
                                    for i in nclose_compress[(st_chr, ed_chr)]:
                                        dummy_list = [0,0,0,0,0,0,0,]
                                        flag = False
                                        if i[3] < i[4]:
                                            if distance_checker(contig_data[st], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                                flag = True
                                                nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                            if distance_checker(contig_data[ed], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                                flag = True
                                                nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                            if flag:
                                                break
                                        else:
                                            if distance_checker(contig_data[st], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                                flag = True
                                                nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                            if distance_checker(contig_data[ed], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                                flag = True
                                                nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                            if flag:
                                                break
                                    nclose_coverage[(st, ed)] += asm2cov[cov_count_name]
                                    nclose_compress[(st_chr, ed_chr)].append(temp_list)
                                    nclose_dict[contig_data[s][CTG_NAM]].append((st, ed))      
                                else:
                                    key = (tuple(censat_st_chr), tuple(censat_ed_chr))
                                    if key not in censat_nclose_compress:
                                        censat_nclose_compress.add(key)
                                        temp_list = [upd_contig_name, 
                                            [contig_data[st][CHR_STR], contig_data[st][CHR_END]], 
                                            [contig_data[ed][CHR_STR], contig_data[ed][CHR_END]], st, ed]
                                        for i in nclose_compress[(st_chr, ed_chr)]:
                                            dummy_list = [0,0,0,0,0,0,0,]
                                            flag = False
                                            if i[3] < i[4]:
                                                if distance_checker(contig_data[st], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                                    flag = True
                                                    nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                if distance_checker(contig_data[ed], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                                    flag = True
                                                    nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                if flag:
                                                    break
                                            else:
                                                if distance_checker(contig_data[st], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                                    flag = True
                                                    nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                if distance_checker(contig_data[ed], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                                    flag = True
                                                    nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                if flag:
                                                    break
                                        nclose_coverage[(st, ed)] += asm2cov[cov_count_name]
                                        nclose_compress[(st_chr, ed_chr)].append(temp_list)
                                        nclose_dict[contig_data[s][CTG_NAM]].append((st, ed))  
                        else:
                            for i in nclose_compress[(ed_chr, st_chr)]:
                                if is_curr_ctg_repeat and i[0] in repeat_contig_name:
                                    compress_limit = ALL_REPEAT_NCLOSE_COMPRESS_LIMIT
                                else:
                                    compress_limit = NCLOSE_COMPRESS_LIMIT
                                dummy_list = [0,0,0,0,0,0,0,]
                                if distance_checker(contig_data[st], dummy_list+i[2]) < compress_limit \
                                and distance_checker(contig_data[ed], dummy_list + i[1]) < compress_limit:
                                    sorted_tar_tuple = tuple(sorted((i[3], i[4])))
                                    nclose_coverage[sorted_tar_tuple] += asm2cov[cov_count_name]
                                    nclose_compress_track[sorted_tar_tuple].append((st, ed))
                                    flag = False
                                    break
                            if flag:
                                censat_st_chr = [st_chr[0], 0]
                                censat_ed_chr = [ed_chr[0], 0]
                                if contig_s[CTG_CENSAT] != '0':
                                    cnt = 0
                                    for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                                        if inclusive_checker_tuple(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                                            censat_st_chr[1] = contig_data[st][CHR_NAM] + "." + str(cnt)
                                            break
                                        cnt+=1
                                if contig_e[CTG_CENSAT] != '0':
                                    cnt = 0
                                    for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                                        if inclusive_checker_tuple(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
                                            censat_ed_chr[1] = contig_data[ed][CHR_NAM] + "." + str(cnt)
                                            break
                                        cnt+=1
                                if censat_st_chr[1] == 0 \
                                or censat_ed_chr[1] == 0 :
                                    temp_list = [upd_contig_name,
                                            [contig_data[ed][CHR_STR], contig_data[ed][CHR_END]], 
                                            [contig_data[st][CHR_STR], contig_data[st][CHR_END]], ed, st]
                                    for i in nclose_compress[(ed_chr, st_chr)]:
                                        dummy_list = [0,0,0,0,0,0,0,]
                                        flag = False
                                        if i[3] > i[4]:
                                            if distance_checker(contig_data[st], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                                nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                flag = True
                                            if distance_checker(contig_data[ed], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                                nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                flag = True
                                        else:
                                            if distance_checker(contig_data[st], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                                nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                flag = True
                                            if distance_checker(contig_data[ed], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                                nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                flag = True
                                        if flag:
                                            break
                                    nclose_coverage[(st, ed)] += asm2cov[cov_count_name]
                                    nclose_compress[(ed_chr, st_chr)].append(temp_list)
                                    nclose_dict[contig_data[s][CTG_NAM]].append((st, ed))     
                                else:
                                    key = (tuple(censat_st_chr), tuple(censat_ed_chr))
                                    if key not in censat_nclose_compress:
                                        censat_nclose_compress.add(key)
                                        temp_list = [upd_contig_name,
                                            [contig_data[ed][CHR_STR], contig_data[ed][CHR_END]], 
                                            [contig_data[st][CHR_STR], contig_data[st][CHR_END]], ed, st]
                                        for i in nclose_compress[(ed_chr, st_chr)]:
                                            dummy_list = [0,0,0,0,0,0,0,]
                                            flag = False
                                            if i[3] > i[4]:
                                                if distance_checker(contig_data[st], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                                    nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                    flag = True
                                                if distance_checker(contig_data[ed], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                                    nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                    flag = True
                                            else:
                                                if distance_checker(contig_data[st], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                                    nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                    flag = True
                                                if distance_checker(contig_data[ed], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                                    nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                                    flag = True
                                            if flag:
                                                break
                                        nclose_coverage[(st, ed)] += asm2cov[cov_count_name]
                                        nclose_compress[(ed_chr, st_chr)].append(temp_list)
                                        nclose_dict[contig_data[s][CTG_NAM]].append((st, ed))
        s = e+1
    return nclose_dict, nclose_start_compress, nclose_end_compress, fake_bnd, all_nclose_compress, nclose_coverage, nclose_compress_track

def make_virtual_ord_ctg(contig_data, fake_bnd):
    contig_data_size = len(contig_data)
    s = 0
    virtual_ctg = []
    idx = 0
    while s < contig_data_size:
        e = contig_data[s][CTG_ENDND]
        curr_ctg = contig_data[s][CTG_NAM]
        if curr_ctg in fake_bnd.keys():
            vctg_rng = fake_bnd[curr_ctg]
            for i in range(vctg_rng[0], vctg_rng[1] +1):
                temp_list = list(contig_data[i])
                temp_list[CTG_NAM] = temp_list[CTG_NAM][:-1] + "l"
                temp_list[CTG_TYP] = 3
                temp_list[CTG_STRND] = idx + contig_data_size
                temp_list[CTG_ENDND] = idx + vctg_rng[1] - vctg_rng[0] + contig_data_size
                virtual_ctg.append(temp_list)
            idx+=vctg_rng[1] - vctg_rng[0]+1
        s = e+1
    return virtual_ctg

def extract_telomere_connect_contig(telo_info_path : str) -> dict:
    telomere_connect_contig = defaultdict(list)
    with open(telo_info_path) as f:
        for curr_data in f:
            curr_data = curr_data.rstrip()
            temp_list = curr_data.split("\t")
            chr_info = temp_list[0]
            contig_id = ast.literal_eval(temp_list[1])
            telomere_connect_contig[chr_info].append(contig_id)
    
    return telomere_connect_contig

def extract_telomere_connect_contig_bytuple(telo_info_path : str) -> list:
    telomere_connect_contig = []
    with open(telo_info_path) as f:
        for curr_data in f:
            curr_data = curr_data.rstrip()
            temp_list = curr_data.split("\t")
            chr_info = temp_list[0]
            contig_id = ast.literal_eval(temp_list[1])
            telomere_connect_contig.append((chr_info, contig_id[1]))
    
    return telomere_connect_contig


def telomere_name_sort_key(telo_name):
    telo_name = str(telo_name)
    side = telo_name[-1] if telo_name.endswith(('f', 'b')) else ''
    chrom = telo_name[:-1] if side else telo_name
    try:
        chrom_order = chr2int(chrom)
    except (TypeError, ValueError):
        chrom_order = INF
    return ({'f': 0, 'b': 1}.get(side, 2), chrom_order, chrom, telo_name)


def write_telomere_connected_outputs(prefix: str, telo_edges: list,
                                      contig_data: list) -> None:
    telo_edges_by_name = defaultdict(list)
    with open(f"{prefix}/telomere_connected_list.txt", "wt") as f:
        for telo_name, edge in telo_edges:
            edge = tuple(edge)
            telo_edges_by_name[telo_name].append(edge)
            print(telo_name, edge, sep="\t", file=f)

    with open(f"{prefix}/telomere_connected_list_readable.txt", "wt") as f:
        for telo_name in sorted(telo_edges_by_name, key=telomere_name_sort_key):
            edges = telo_edges_by_name[telo_name]
            print(telo_name, file=f)
            for edge in edges:
                print(tuple(contig_data[edge[1]]), file=f)
            print("", file=f)


def group_nclose_nodes_by_chrom(contig_data: list, nclose_nodes: dict) -> dict:
    nclose_type = defaultdict(list)
    for pair_list in nclose_nodes.values():
        for pair in pair_list:
            contig_a = contig_data[pair[0]]
            contig_b = contig_data[pair[1]]
            if chr2int(contig_a[CHR_NAM]) <= chr2int(contig_b[CHR_NAM]):
                chrom_pair = (contig_a[CHR_NAM], contig_b[CHR_NAM])
                normalized_pair = tuple(pair)
            else:
                chrom_pair = (contig_b[CHR_NAM], contig_a[CHR_NAM])
                normalized_pair = (pair[1], pair[0])
            nclose_type[chrom_pair].append(normalized_pair)
    return nclose_type


def write_nclose_nodes_list(path: str, nclose_type: dict, contig_data: list,
                            repeat_contig_names=()) -> int:
    node_count = 0
    with open(path, "wt") as f:
        for chrom_pair, pair_list in nclose_type.items():
            print(f"{chrom_pair[0]}, {chrom_pair[1]}, {len(pair_list)}", file=f)
            for pair in pair_list:
                node_count += 2
                contig_a = contig_data[pair[0]]
                contig_b = contig_data[pair[1]]
                is_for = pair[0] < pair[1]
                list_a = [
                    contig_a[CTG_NAM],
                    get_corr_dir(is_for, contig_a[CTG_DIR]),
                    contig_a[CHR_STR],
                    contig_a[CHR_END],
                ]
                list_b = [
                    contig_b[CTG_NAM],
                    get_corr_dir(is_for, contig_b[CTG_DIR]),
                    contig_b[CHR_STR],
                    contig_b[CHR_END],
                ]
                if contig_a[CTG_NAM] in repeat_contig_names:
                    print(list_a, list_b, "all_repeat", file=f)
                else:
                    print(list_a, list_b, file=f)
            print("", file=f)
    return node_count


def write_nclose_nodes_index(path: str, nclose_nodes: dict,
                             contig_data: list) -> int:
    node_count = 0
    with open(path, "wt") as f:
        for contig_name, pair_list in nclose_nodes.items():
            for pair in pair_list:
                node_count += 2
                print(
                    contig_name,
                    pair[0],
                    pair[1],
                    contig_data[pair[0]][CTG_TYP],
                    file=f,
                )
    return node_count


def write_virtual_ordinary_contig(path: str, contig_data: list) -> None:
    with open(path, "wt") as f:
        for contig in contig_data:
            for field in contig:
                print(field, end="\t", file=f)
            print("", file=f)
    
def initialize_bnd_graph(contig_data : list, nclose_nodes : dict, telo_contig : dict) -> dict:
    # Build telo direction lookup: telo_idx → 'f' or 'b'
    telo_dir_map = {}
    for telo_name in telo_contig:
        for edge in telo_contig[telo_name]:
            telo_dir_map[edge[1]] = telo_name[-1]

    def is_telo_dir_consistent(telo_idx, approach_dir):
        """
        Check if approach direction is consistent with telomere position.
        'b' telomere (chr end): consistent direction = "dec" (from right)
        'f' telomere (chr start): consistent direction = "inc" (from left)
        "mixed" is always inconsistent.
        Applies to both DIR_IN and DIR_OUT edges.
        """
        telo_side = telo_dir_map.get(telo_idx)
        if telo_side is None:
            return True
        return (telo_side == 'b' and approach_dir == "dec") or \
               (telo_side == 'f' and approach_dir == "inc")

    bnd_adjacency = defaultdict(list)
    for i in telo_contig:
        for j in telo_contig[i]:
            bnd_adjacency[i].append(j[:2])
            bnd_adjacency[(DIR_IN, j[1])].append(i)
            bnd_adjacency[i].append((DIR_OUT, j[1]))
    for j in nclose_nodes:
        for i in nclose_nodes[j]:
            bnd_adjacency[(DIR_FOR, i[0])].append([DIR_FOR, i[1]])
            bnd_adjacency[(DIR_BAK, i[1])].append([DIR_BAK, i[0]])

    nclose_nodes_key = list(nclose_nodes.keys())
    key_len = len(nclose_nodes)
    # nclose nodes connection.
    for i1 in range(key_len):
        for i2 in range(i1+1, key_len):
            for n1 in nclose_nodes[nclose_nodes_key[i1]]:
                for n2 in nclose_nodes[nclose_nodes_key[i2]]:
                    n1_len = len(n1)
                    assert(n1_len==2)
                    n2_len = len(n2)
                    assert(n2_len==2)
                    for i in range(n1_len):
                        for j in range(n2_len):
                            if contig_data[n1[i]][CHR_NAM]==contig_data[n2[j]][CHR_NAM]:
                                contig_i = contig_data[n1[i]]
                                contig_j = contig_data[n2[j]]
                                i_ind = 'b' if i%2 else 'f'
                                j_ind = 'b' if j%2 else 'f'
                                i_ind += contig_i[CTG_DIR]
                                j_ind += contig_j[CTG_DIR]
                                i_str = contig_i[CHR_STR]
                                j_str = contig_j[CHR_STR]
                                i_end = contig_i[CHR_END]
                                j_end = contig_j[CHR_END]
                                is_overlap = distance_checker(contig_i, contig_j)==0
                                if i_ind == 'f+':
                                    if j_ind == 'f-':
                                        if i_str > j_str or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_FOR, n2[j]])
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_FOR, n1[i]])
                                    elif j_ind == 'b+':
                                        if i_str > j_end or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_FOR, n1[i]])
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_BAK, n2[j]])
                                elif i_ind == 'b+':
                                    if j_ind == 'b-':
                                        if j_end > i_end or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_BAK, n2[j]])
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_BAK, n1[i]])
                                    elif j_ind == 'f+':
                                        if j_str > i_end or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_BAK, n1[i]])
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_FOR, n2[j]])
                                elif i_ind == 'f-':
                                    if j_ind == 'f+':
                                        if j_str > i_str or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_FOR, n2[j]])
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_FOR, n1[i]])
                                    elif j_ind == 'b-':
                                        if j_end > i_str or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_FOR, n1[i]])
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_BAK, n2[j]])
                                elif i_ind == 'b-':
                                    if j_ind == 'b+':
                                        if i_end > j_end or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_BAK, n2[j]])
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_BAK, n1[i]])
                                    elif j_ind == 'f-':
                                        if i_end > j_str or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_BAK, n1[i]])
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_FOR, n2[j]])

    # telo connected nodes - nclose nodes connection.
    for i in telo_contig:
        for j in telo_contig[i]:
            curr_contig = contig_data[j[1]]
            for nodes in nclose_nodes:
                for node in nclose_nodes[nodes]:
                    if j[1] in node:
                        continue
                    # curr_contig : telo-connected node
                    # index : j[1] and k
                    for k in range(len(node)):
                        nclose_contig = contig_data[node[k]]
                        telo_idx = j[1]
                        nclose_idx = node[k]
                        if nclose_contig[CHR_NAM] != curr_contig[CHR_NAM]:
                            continue
                        k_ind = 'b' if k%2 else 'f'
                        k_ind += nclose_contig[CTG_DIR]
                        k_str = nclose_contig[CHR_STR]
                        k_end = nclose_contig[CHR_END]
                        telo_fb = i[-1]=='f'
                        if curr_contig[CTG_DIR] == '+':
                            t_ind = 'f' if telo_fb else 'b'
                        else:
                            t_ind = 'b' if telo_fb else 'f'
                        t_ind += curr_contig[CTG_DIR]
                        t_str = curr_contig[CHR_STR]
                        t_end = curr_contig[CHR_END]
                        i_ind = k_ind
                        j_ind = t_ind
                        i_str = k_str
                        i_end = k_end
                        j_str = t_str
                        j_end = t_end
                        is_overlap = distance_checker(curr_contig, nclose_contig)==0
                        if curr_contig[CHR_END] <= nclose_contig[CHR_END]:
                            dir1 = "inc"
                        elif nclose_contig[CHR_END] < curr_contig[CHR_END]:
                            dir1 = "dec"
                        if curr_contig[CHR_STR] <= nclose_contig[CHR_STR]:
                            dir2 = "inc"
                        elif nclose_contig[CHR_STR] < curr_contig[CHR_STR]:
                            dir2 = "dec"
                        # if both end have consistency
                        if dir1 == dir2:
                            consistent = is_telo_dir_consistent(telo_idx, dir1)
                            if dir1 == "inc":
                                if curr_contig[CTG_DIR]=='+' and k_ind=='f+': # +
                                    if consistent:
                                        bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                        bnd_adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                                elif curr_contig[CTG_DIR]=='-' and k_ind=='b-':
                                    if consistent:
                                        bnd_adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                                        bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                elif curr_contig[CTG_DIR]=='+' and k_ind=='b-':
                                    if consistent:
                                        bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                        bnd_adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                                elif curr_contig[CTG_DIR]=='-' and k_ind=='f+': # +
                                    if consistent:
                                        bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                        bnd_adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                            else:
                                if curr_contig[CTG_DIR]=='+' and k_ind=='b+': # +
                                    if consistent:
                                        bnd_adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                                        bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                elif curr_contig[CTG_DIR]=='-' and k_ind=='f-':
                                    if consistent:
                                        bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                        bnd_adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                                elif curr_contig[CTG_DIR]=='+' and k_ind=='f-':
                                    if consistent:
                                        bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                        bnd_adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                                elif curr_contig[CTG_DIR]=='-' and k_ind=='b+': # +
                                    if consistent:
                                        bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                        bnd_adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                        # else: each end have different sign of increment
                        # -> one node is included in the other
                        # thus, we should determine by mean val.
                        else:
                            # mixed direction is always inconsistent → skip both DIR_IN and DIR_OUT
                            pass
                            

    
    for i in range(contig_data_size, contig_data_size + CHROMOSOME_COUNT):
        j = i + CHROMOSOME_COUNT
        for node_a in telo_contig[chr_rev_corr[i]]:
            for node_b in telo_contig[chr_rev_corr[j]]:
                ind = contig_data[node_a[1]][CTG_DIR] + contig_data[node_b[1]][CTG_DIR]
                if ind == "++":
                    bnd_adjacency[(DIR_OUT, node_a[1])].append([DIR_IN, node_b[1]])
                    bnd_adjacency[(DIR_OUT, node_b[1])].append([DIR_IN, node_a[1]])
                elif ind == "+-":
                    bnd_adjacency[(DIR_OUT, node_a[1])].append([DIR_IN, node_b[1]])
                    bnd_adjacency[(DIR_OUT, node_b[1])].append([DIR_IN, node_a[1]])
                elif ind=="-+":
                    bnd_adjacency[(DIR_OUT, node_a[1])].append([DIR_IN, node_b[1]])
                    bnd_adjacency[(DIR_OUT, node_b[1])].append([DIR_IN, node_a[1]])
                elif ind=="--":
                    bnd_adjacency[(DIR_OUT, node_a[1])].append([DIR_IN, node_b[1]])
                    bnd_adjacency[(DIR_OUT, node_b[1])].append([DIR_IN, node_a[1]])


    return bnd_adjacency

def initialize_inversion_only_graph(contig_data : list, nclose_nodes_original : dict) -> dict:
    bnd_adjacency = defaultdict(list)
    nclose_nodes = defaultdict(list)

    for j in nclose_nodes_original:
        for i in nclose_nodes_original[j]:
            if contig_data[i[0]][CHR_NAM] == contig_data[i[1]][CHR_NAM]:
                nclose_nodes[j].append(i)

    for j in nclose_nodes:
        for i in nclose_nodes[j]:
            bnd_adjacency[(DIR_FOR, i[0])].append([DIR_FOR, i[1]])
            bnd_adjacency[(DIR_BAK, i[1])].append([DIR_BAK, i[0]])
            if len(i)==4:
                bnd_adjacency[(DIR_FOR, i[2])].append([DIR_FOR, i[3]])
                bnd_adjacency[(DIR_BAK, i[3])].append([DIR_BAK, i[2]])
    nclose_nodes_key = list(nclose_nodes.keys())
    key_len = len(nclose_nodes)
    # nclose nodes connection.
    for i1 in range(key_len):
        for i2 in range(i1+1, key_len):
            for n1 in nclose_nodes[nclose_nodes_key[i1]]:
                for n2 in nclose_nodes[nclose_nodes_key[i2]]:
                    n1_len = len(n1)
                    assert(n1_len==2)
                    n2_len = len(n2)
                    assert(n2_len==2)
                    for i in range(n1_len):
                        for j in range(n2_len):
                            if contig_data[n1[i]][CHR_NAM]==contig_data[n2[j]][CHR_NAM]:
                                contig_i = contig_data[n1[i]]
                                contig_j = contig_data[n2[j]]
                                i_ind = 'b' if i%2 else 'f'
                                j_ind = 'b' if j%2 else 'f'
                                i_ind += contig_i[CTG_DIR]
                                j_ind += contig_j[CTG_DIR]
                                i_str = contig_i[CHR_STR]
                                j_str = contig_j[CHR_STR]
                                i_end = contig_i[CHR_END]
                                j_end = contig_j[CHR_END]
                                is_overlap = distance_checker(contig_i, contig_j)==0
                                if i_ind == 'f+':
                                    if j_ind == 'f-':
                                        if i_str > j_str or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_FOR, n2[j]])
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_FOR, n1[i]])
                                    elif j_ind == 'b+':
                                        if i_str > j_end or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_FOR, n1[i]])
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_BAK, n2[j]])
                                elif i_ind == 'b+':
                                    if j_ind == 'b-':
                                        if j_end > i_end or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_BAK, n2[j]])
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_BAK, n1[i]])
                                    elif j_ind == 'f+':
                                        if j_str > i_end or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_BAK, n1[i]])
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_FOR, n2[j]])
                                elif i_ind == 'f-':
                                    if j_ind == 'f+':
                                        if j_str > i_str or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_FOR, n2[j]])
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_FOR, n1[i]])
                                    elif j_ind == 'b-':
                                        if j_end > i_str or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_FOR, n1[i]])
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_BAK, n2[j]])
                                elif i_ind == 'b-':
                                    if j_ind == 'b+':
                                        if i_end > j_end or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_BAK, n2[j]])
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_BAK, n1[i]])
                                    elif j_ind == 'f-':
                                        if i_end > j_str or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_BAK, n1[i]])
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_FOR, n2[j]])

    return bnd_adjacency


def convert_all_nclose_comp_to_nclose_nodes(contig_data : list, all_nclose_comp : dict) -> dict:
    all_nclose_nodes = defaultdict(list)
    seen_nclose = set()

    for pairs in all_nclose_comp.values():
        for pair in pairs:
            pair = tuple(pair)
            if len(pair) != 2 or pair in seen_nclose:
                continue
            all_nclose_nodes[contig_data[pair[0]][CTG_NAM]].append(pair)
            seen_nclose.add(pair)

    return all_nclose_nodes


def build_ecdna_nclose_nodes(raw_nclose_nodes : dict, all_nclose_nodes : dict) -> dict:
    ecdna_nclose_nodes = defaultdict(list)
    seen_nclose = set()

    for nclose_source in (raw_nclose_nodes, all_nclose_nodes):
        for contig_name, pairs in nclose_source.items():
            for pair in pairs:
                pair = tuple(pair)
                if len(pair) != 2 or pair in seen_nclose:
                    continue
                ecdna_nclose_nodes[contig_name].append(pair)
                seen_nclose.add(pair)

    return ecdna_nclose_nodes


def make_inversion_nx_graph(bnd_graph_adjacency):
    G = nx.DiGraph()
    for i in bnd_graph_adjacency:
        G.add_node(i)
    for node in bnd_graph_adjacency:
        for edge in bnd_graph_adjacency[node]:
            if type(edge)==str:
                G.add_weighted_edges_from([(tuple(node), edge, 0)])
            elif type(node)==str:
                G.add_weighted_edges_from([(node, tuple(edge), 0)])
            else:
                d = 0
                contig_s = contig_data[node[1]]
                contig_e = contig_data[edge[1]]
                if contig_s[CTG_NAM] == contig_e[CTG_NAM]:
                    if node[1] < edge[1]:
                        d = 0
                        for i in range(node[1]+1, edge[1]):
                            d += contig_data[i][CHR_END] - contig_data[i][CHR_STR]
                    else:
                        d = 0
                        for i in range(edge[1]+1, node[1]):
                            d += contig_data[i][CHR_END] - contig_data[i][CHR_STR]
                            
                else:
                    ds = contig_s[CHR_END] - contig_s[CHR_STR]
                    de = contig_e[CHR_END] - contig_e[CHR_STR]
                    if distance_checker(contig_s, contig_e)==0:
                        d = ds + de - overlap_calculator(contig_s, contig_e)
                    else:  
                        d = distance_checker(contig_s, contig_e) + ds + de 
                G.add_weighted_edges_from([(tuple(node), tuple(edge), d)])
    return G

def circuit_length_calculator(circuit):
    circuit_len = 0
    for i in range(len(circuit)):
        curr_node = circuit[i]
        next_node = circuit[i+1] if i+1 < len(circuit) else circuit[0]
        circuit_len += abs(contig_data[curr_node][CHR_END] - contig_data[curr_node][CHR_STR])
        circuit_len += abs(contig_data[curr_node][CHR_END] - contig_data[next_node][CHR_STR])
    return circuit_len

def contig_preprocessing_00(PAF_FILE_PATH_ : list):

    original_node_count = 0
    excluded_telomere_origins = set()

    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
    telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)
    repeat_data = import_repeat_data_00(REPEAT_INFO_FILE_PATH)
    repeat_censat_data = import_censat_repeat_data(CENSAT_PATH)

    global df, no_chrY
    df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
    df = df.query('chr != "chrM"')

    cen_vtg_contig, cen_fragment_meta = find_breakend_centromere(
        repeat_censat_data,
        chr_len,
        df,
        log_context="Strict",
    )

    with open(f'{PREFIX}/cen_fragment_data.pkl', 'wb') as f:
        pkl.dump(cen_fragment_meta, f)

    telo_dict = defaultdict(list)
    for _ in telo_data:
        telo_dict[_[0]].append(_[1:])

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

    contig_data = import_data(PAF_FILE_PATH_[0])
    excluded_contigs, excluded_rows = find_multi_end_aligned_contigs(contig_data)
    excluded_telomere_origins.update((0, row_idx) for row_idx in excluded_rows)
    if excluded_contigs:
        logging.info(
            f"Detected {len(excluded_contigs)} multi-end-aligned contigs "
            f"({len(excluded_rows)} PAF rows) in {PAF_FILE_PATH_[0]}"
        )

    original_node_count += len(contig_data)

    node_label = label_node(contig_data, telo_dict)

    repeat_label = label_repeat_node(contig_data, repeat_data, chr_len)

    telo_preprocessed_contig, report_case, telo_connect_info = preprocess_telo(contig_data, node_label)
    excluded_telo_candidates = sum(
        row_idx in excluded_rows for row_idx in telo_connect_info
    )
    telo_connect_info = {
        row_idx: telo_name
        for row_idx, telo_name in telo_connect_info.items()
        if row_idx not in excluded_rows
    }
    if excluded_telo_candidates:
        logging.info(
            f"Excluded {excluded_telo_candidates} primary PAF telomere boundary candidates "
            "from multi-end-aligned contigs"
        )

    new_contig_data = []
    telcon_set = set()
    idxcnt = 0
    for i in telo_preprocessed_contig:
        temp_list = contig_data[i]
        if i in telo_connect_info:
            temp_list.append(telo_connect_info[i])
            telcon_set.add(idxcnt)
        else:
            temp_list.append("0")
        new_contig_data.append(temp_list)
        idxcnt+=1

    # with open("telo_preprocess_contig.txt", "wt") as f:
    #     for i in new_contig_data:
    #         for j in i:
    #             print(j, end="\t", file=f)
    #         print("", file=f)
    
    new_node_repeat_label = label_repeat_node(new_contig_data, repeat_data, chr_len)
    new_node_repeat_censat_label = label_repeat_node(new_contig_data, repeat_censat_data, chr_len)
    new_node_telo_label = label_node(new_contig_data, telo_dict)

    ref_qry_ratio = calc_ratio(new_contig_data)
    preprocess_result, \
    preprocess_contig_type, \
    preprocess_terminal_nodes, \
    len_counter = preprocess_contig(new_contig_data, new_node_telo_label, ref_qry_ratio, new_node_repeat_label, telcon_set)
    preprocess_result = set(preprocess_result)
    new_contig_data = new_contig_data[0:-1]
    contig_data_size = len(new_contig_data)
    telo_ppc_contig = []
    for i in range(contig_data_size):
        if new_contig_data[i][CTG_NAM] in preprocess_result:
            temp_list = new_contig_data[i][:10]
            temp_list.append(preprocess_contig_type[new_contig_data[i][CTG_NAM]])
            temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][0])
            temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][1])
            temp_list.append(new_node_telo_label[i][0])
            temp_list.append(new_node_telo_label[i][1])
            temp_list.append(new_contig_data[i][11])
            temp_list.append(new_node_repeat_label[i][0])
            temp_list.append(new_node_repeat_label[i][1])
            temp_list.append(new_node_repeat_censat_label[i][1])
            temp_list.append('0.'+str(new_contig_data[i][10]))
            telo_ppc_contig.append(temp_list)
    mainflow_dict = find_mainflow(telo_ppc_contig)
    telo_ppc_size = len(telo_ppc_contig)
    for i in range(telo_ppc_size):
        max_chr = mainflow_dict[telo_ppc_contig[i][CTG_NAM]]
        temp = [max_chr[0], max_chr[1], telo_ppc_contig[i][-1]]
        telo_ppc_contig[i] = telo_ppc_contig[i][:-1]
        telo_ppc_contig[i] += temp

    final_contig = preprocess_repeat(telo_ppc_contig)

    final_contig_repeat_label = label_repeat_node(final_contig, repeat_data, chr_len)

    final_telo_node_label = label_node(final_contig, telo_dict)

    final_ref_qry_ratio = calc_ratio(final_contig)

    final_telo_connect = set()

    for v, i in enumerate(final_contig):
        try:
            if i[CTG_TELCON] != '0':
                final_telo_connect.add(v)
        except:
            pass


    final_using_contig, final_ctg_typ, final_preprocess_terminal_nodes, _ = preprocess_contig(final_contig, final_telo_node_label, final_ref_qry_ratio, final_contig_repeat_label, final_telo_connect)

    final_contig = final_contig[:-1]

    real_final_contig = []

    final_using_contig = set(final_using_contig)
    for i in range(0, len(final_contig)):
        if final_contig[i][CTG_NAM] in final_using_contig:
            final_contig[i][CTG_TYP] = final_ctg_typ[final_contig[i][CTG_NAM]]
            final_contig[i][CTG_STRND] = final_preprocess_terminal_nodes[final_contig[i][CTG_NAM]][0]
            final_contig[i][CTG_ENDND] = final_preprocess_terminal_nodes[final_contig[i][CTG_NAM]][1]
            real_final_contig.append(final_contig[i])
    total_len = len(real_final_contig)

    # alt process
    if args.alt != None:
        contig_data = import_data(PAF_FILE_PATH_[1])
        excluded_contigs, excluded_rows = find_multi_end_aligned_contigs(contig_data)
        excluded_telomere_origins.update((1, row_idx) for row_idx in excluded_rows)
        if excluded_contigs:
            logging.info(
                f"Detected {len(excluded_contigs)} multi-end-aligned contigs "
                f"({len(excluded_rows)} PAF rows) in {PAF_FILE_PATH_[1]}"
            )
        original_node_count += len(contig_data)
        node_label = label_node(contig_data, telo_dict)
        telo_preprocessed_contig, report_case, telo_connect_info = preprocess_telo(contig_data, node_label)
        excluded_telo_candidates = sum(
            row_idx in excluded_rows for row_idx in telo_connect_info
        )
        telo_connect_info = {
            row_idx: telo_name
            for row_idx, telo_name in telo_connect_info.items()
            if row_idx not in excluded_rows
        }
        if excluded_telo_candidates:
            logging.info(
                f"Excluded {excluded_telo_candidates} alternative PAF telomere boundary candidates "
                "from multi-end-aligned contigs"
            )
        new_contig_data = []
        telcon_set = set()
        idxcnt = 0
        for i in telo_preprocessed_contig:
            temp_list = contig_data[i]
            if i in telo_connect_info:
                temp_list.append(telo_connect_info[i])
                telcon_set.add(idxcnt)
            else:
                temp_list.append("0")
            new_contig_data.append(temp_list)
            idxcnt+=1
        alt_telo_ppc_contig = []
        new_node_telo_label = label_node(new_contig_data, telo_dict)
        new_node_repeat_label = label_repeat_node(new_contig_data, repeat_data, chr_len)
        new_node_repeat_censat_label = label_repeat_node(new_contig_data, repeat_censat_data, chr_len)

        ref_qry_ratio = calc_ratio(new_contig_data)
        preprocess_result, \
        preprocess_type3_result, \
        preprocess_contig_type, \
        preprocess_terminal_nodes, \
        alt_len_counter = alt_preprocess_contig(new_contig_data, new_node_telo_label, ref_qry_ratio, new_node_repeat_label, telcon_set, telo_dict, repeat_censat_data)
        new_contig_data = new_contig_data[0:-1]
        contig_data_size = len(new_contig_data)
        for contig_name_list in [set(preprocess_result), set(preprocess_type3_result)]:
            for i in range(contig_data_size):
                if new_contig_data[i][CTG_NAM] in contig_name_list:
                    c_type = preprocess_contig_type[new_contig_data[i][CTG_NAM]]
                    if c_type == 5:
                        if i in telcon_set:
                            new_type = 3
                            temp_list = new_contig_data[i][:10]
                            temp_list.append(new_type)
                            temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][0])
                            temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][1])
                            temp_list.append(new_node_telo_label[i][0])
                            temp_list.append(new_node_telo_label[i][1])
                            temp_list.append(new_contig_data[i][11])
                            temp_list.append(new_node_repeat_label[i][0])
                            temp_list.append(new_node_repeat_label[i][1])
                            temp_list.append(new_node_repeat_censat_label[i][1])
                            temp_list.append('1.'+str(new_contig_data[i][10]))
                            alt_telo_ppc_contig.append(temp_list)
                    else:
                        temp_list = new_contig_data[i][:10]
                        temp_list.append(preprocess_contig_type[new_contig_data[i][CTG_NAM]])
                        temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][0])
                        temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][1])
                        temp_list.append(new_node_telo_label[i][0])
                        temp_list.append(new_node_telo_label[i][1])
                        temp_list.append(new_contig_data[i][11])
                        temp_list.append(new_node_repeat_label[i][0])
                        temp_list.append(new_node_repeat_label[i][1])
                        temp_list.append(new_node_repeat_censat_label[i][1])
                        temp_list.append('1.'+str(new_contig_data[i][10]))
                        alt_telo_ppc_contig.append(temp_list)
        # with open("alt_telo_ppc_contig.txt", "wt") as f:
        #     for i in alt_telo_ppc_contig:
        #         for j in i:
        #             print(j, end="\t", file=f)
        #         print("", file=f)
        alt_mainflow_dict = find_mainflow(alt_telo_ppc_contig)
        alt_telo_ppc_size = len(alt_telo_ppc_contig)
        for i in range(alt_telo_ppc_size):
            max_chr = alt_mainflow_dict[alt_telo_ppc_contig[i][CTG_NAM]]
            temp = [max_chr[0], max_chr[1], alt_telo_ppc_contig[i][-1]]
            alt_telo_ppc_contig[i] = alt_telo_ppc_contig[i][:-1]
            alt_telo_ppc_contig[i] += temp
        
        alt_final_contig = preprocess_repeat(alt_telo_ppc_contig)

        alt_final_repeat_node_label = label_repeat_node(alt_final_contig, repeat_data, chr_len)

        alt_final_telo_node_label = label_node(alt_final_contig, telo_dict)

        alt_final_telo_connect = set()

        for v, i in enumerate(alt_final_contig):
            try:
                if i[CTG_TELCON] != '0':
                    alt_final_telo_connect.add(v)
            except:
                pass

        alt_final_ref_qry_ratio = calc_ratio(alt_final_contig)

        alt_final_using_contig, alt_final_using_type3_contig, \
        alt_final_ctg_typ, alt_final_preprocess_terminal_nodes, _ = alt_preprocess_contig(alt_final_contig, alt_final_telo_node_label, alt_final_ref_qry_ratio, alt_final_repeat_node_label, alt_final_telo_connect, telo_dict, repeat_censat_data)
        
        alt_final_contig = alt_final_contig[:-1]
        
        bias = len(real_final_contig)

        real_alt_final_contig = []

        for alt_final_using_contig_ in [set(alt_final_using_contig), set(alt_final_using_type3_contig)]:
            for i in range(0, len(alt_final_contig)):
                if alt_final_contig[i][CTG_NAM] in alt_final_using_contig_:
                    if alt_final_ctg_typ[alt_final_contig[i][CTG_NAM]] == 5:
                        if alt_final_contig[i][CTG_TELCON] != '0':
                            alt_final_contig[i][CTG_TYP] = 3
                            alt_final_contig[i][CTG_STRND] = alt_final_preprocess_terminal_nodes[alt_final_contig[i][CTG_NAM]][0] + bias
                            alt_final_contig[i][CTG_ENDND] = alt_final_preprocess_terminal_nodes[alt_final_contig[i][CTG_NAM]][1] + bias
                            real_alt_final_contig.append(alt_final_contig[i])
                    else:
                        alt_final_contig[i][CTG_TYP] = alt_final_ctg_typ[alt_final_contig[i][CTG_NAM]]
                        alt_final_contig[i][CTG_STRND] = alt_final_preprocess_terminal_nodes[alt_final_contig[i][CTG_NAM]][0] + bias
                        alt_final_contig[i][CTG_ENDND] = alt_final_preprocess_terminal_nodes[alt_final_contig[i][CTG_NAM]][1] + bias
                        real_alt_final_contig.append(alt_final_contig[i])
        
        real_final_contig = real_final_contig + real_alt_final_contig
        total_len = len(real_final_contig)

    # with open("a.txt", "wt") as f:
    #     for i in real_final_contig:
    #         for j in i:
    #             print(j, end="\t", file=f)
    #         print("", file=f)
    
    telo_fb_dict = defaultdict(list)
    for k, v in telo_dict.items():
        for i in v:
            telo_fb_dict[k+i[-1]].append([i[0], i[1]])
    
    telo_bound_dict = {}
    for k, v in telo_fb_dict.items():
        if k[-1]=='f':
            telo_bound_dict[k] = min(v, key=lambda t:t[0])
        else:
            telo_bound_dict[k] = max(v, key=lambda t:t[1])

    telo_label = label_node(real_final_contig, telo_dict)
    subtelo_label = label_subtelo_node(real_final_contig, telo_dict)
    subtelo_ppc_node = subtelo_cut(real_final_contig, telo_label, subtelo_label)

    adjacency = initial_graph_build(real_final_contig, telo_bound_dict)

    telo_connected_node, telo_connected_dict, telo_connected_graph_dict, telo_coverage = edge_optimization(
        real_final_contig,
        adjacency,
        telo_bound_dict,
        asm2cov,
        excluded_telomere_origins,
    )
    telo_connected_node, telo_connected_dict, telo_connected_graph_dict = filter_telomere_connected_cen_fragment_mismatch(
        real_final_contig,
        telo_connected_graph_dict,
        cen_fragment_meta,
        "pre-break",
    )
    
    break_contig = break_double_telomere_contig(real_final_contig, telo_connected_node)

    type34_split_contig, broken_contig_set = break_type34_contig(real_final_contig)

    real_final_contig = delete_contig(real_final_contig, broken_contig_set)

    total_len = len(real_final_contig)

    if len(break_contig) > 0:
        final_break_contig = pass_pipeline(break_contig, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, False, chr_len)
    else:
        final_break_contig = []

    if len(cen_vtg_contig) > 0:
        final_cen_vtg_contig = pass_pipeline(cen_vtg_contig, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, False, chr_len)
    else:
        final_cen_vtg_contig = []

    if len(subtelo_ppc_node) > 0:
        final_subtelo_ppc_node = pass_pipeline(subtelo_ppc_node, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, True, chr_len)
    else:
        final_subtelo_ppc_node = []

    if len(type34_split_contig) > 0:
        final_type34_split_contig = pass_pipeline(type34_split_contig, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, False, chr_len)
    else:
        final_type34_split_contig = []

    for i in final_break_contig:
        temp_list = i
        temp_list[CTG_STRND] += total_len
        temp_list[CTG_ENDND] += total_len
        real_final_contig.append(temp_list)

    total_len += len(final_break_contig)

    for i in final_cen_vtg_contig:
        temp_list = i
        temp_list[CTG_STRND] += total_len
        temp_list[CTG_ENDND] += total_len
        real_final_contig.append(temp_list)

    total_len += len(final_cen_vtg_contig)

    for i in final_subtelo_ppc_node:
        temp_list = i
        temp_list[CTG_STRND] += total_len
        temp_list[CTG_ENDND] += total_len
        real_final_contig.append(temp_list)
    
    total_len += len(final_subtelo_ppc_node)

    for i in final_type34_split_contig:
        temp_list = i
        temp_list[CTG_STRND] += total_len
        temp_list[CTG_ENDND] += total_len
        real_final_contig.append(temp_list)

    total_len += len(final_type34_split_contig)

    overlap_low_split_contig = []
    if is_unitig_reduced == False:
        ctgname2overlap = get_overlap_total_score_dict(ORIGINAL_PAF_LOC_LIST)
        not_trust_contig_name = get_not_trust_contig_name(ORIGINAL_PAF_LOC_LIST)

        target_split_contig_nameset = set()
        for k, v in ctgname2overlap.items():
            if v <= 2:
                target_split_contig_nameset.add(k)

        s = 0
        r_l = len(real_final_contig)
        split_contig_counter = Counter()
        while s < r_l:
            contig_name = real_final_contig[s][CTG_NAM]
            e = real_final_contig[s][CTG_ENDND]
            if contig_name in target_split_contig_nameset \
               and real_final_contig[s][CTG_GLOBALIDX].split(".")[0]=='1' \
               and e-s >= 2:
                
                split_limit = SPLIT_CTG_LEN_LIMIT
                if contig_name not in not_trust_contig_name:
                    split_limit = TRUST_SPLIT_CTG_LEN_LIMIT
                
                chunk_start_idx = s
                current_idx = s + 1

                while current_idx <= e:
                    start_chr_info = (real_final_contig[chunk_start_idx][CTG_DIR], real_final_contig[chunk_start_idx][CHR_NAM])
                    current_chr_info = (real_final_contig[current_idx][CTG_DIR], real_final_contig[current_idx][CHR_NAM])
                    is_same_chromosome = (start_chr_info == current_chr_info)

                    is_ratio_valid = False
                    if is_same_chromosome:
                        chunk_slice = real_final_contig[chunk_start_idx : current_idx + 1]
                        curr_fragment_ratio, _ = calculate_single_contig_ref_ratio(chunk_slice)
                        if abs(curr_fragment_ratio - 1) < BND_CONTIG_BOUND:
                            is_ratio_valid = True

                    if not is_same_chromosome or not is_ratio_valid:
                        prev_element_idx = current_idx - 1
                        
                        if real_final_contig[current_idx][CTG_END] - real_final_contig[prev_element_idx][CTG_STR] >= split_limit:
                            split_contig_counter[contig_name] += 1
                            
                            temp_list = copy.deepcopy(real_final_contig[prev_element_idx : current_idx + 1])
                            
                            new_name = f"split_contig_{contig_name}_{split_contig_counter[contig_name]}"
                            temp_list[0][CTG_MAPQ] = temp_list[1][CTG_MAPQ] = 60
                            temp_list[0][CTG_NAM] = temp_list[1][CTG_NAM] = new_name
                            
                            overlap_low_split_contig += temp_list
                        
                        chunk_start_idx = current_idx
                        
                    current_idx += 1
            s = e+1

    if len(overlap_low_split_contig) > 0:
        final_low_split_contig = pass_pipeline(overlap_low_split_contig, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, False, chr_len)
    else:
        final_low_split_contig = []

    for i in final_low_split_contig:
        temp_list = i
        temp_list[CTG_STRND] += total_len
        temp_list[CTG_ENDND] += total_len
        real_final_contig.append(temp_list)
    
    total_len += len(final_low_split_contig)

    add_node_count = len(final_break_contig) + len(final_cen_vtg_contig) + len(final_low_split_contig)
    logging.info(f"Original PAF file length : {original_node_count}")
    logging.info(f"Final preprocessed PAF file length: {len(real_final_contig)}")
    logging.info(f"Number of virtual contigs added on preprocessing : {add_node_count}")

    # Simple ctg-as-alt (enabled by default; disabled with --disable_alt_ctg_simple):
    # primary contig PAF에서 아직 채택되지 않은 contig를 보되, telomere-like terminal chunk와
    # 10kb 이하 fragment를 빼고 90% major chromosome set에서 chr-change 후보를 뽑는다.
    # 같은 방향 chr-change만 기존 nclose와 비교해 10kb 이내면 중복으로 보고 제외하고,
    # 방향이 바뀌는 chr-change는 필터링하지 않고 그대로 유지한다.
    if not args.disable_alt_ctg_simple:
        primary_kept_set = set(final_using_contig)
        existing_names = {c[CTG_NAM] for c in real_final_contig}
        existing_nclose_loci = collect_existing_alt_simple_nclose_loci(real_final_contig)

        chunks_per_contig = defaultdict(list)
        # primary ctg PAF의 line idx를 함께 트래킹 → CTG_GLOBALIDX="0.{line_idx}"로 11번이 lookup 가능
        with open(PAF_FILE_PATH_[0], "rt") as f:
            for line_idx, line in enumerate(f):
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue
                ctg_name = cols[0]
                if ctg_name in primary_kept_set:
                    continue
                try:
                    chunk = [
                        ctg_name,        # 0  CTG_NAM
                        int(cols[1]),    # 1  CTG_LEN
                        int(cols[2]),    # 2  CTG_STR
                        int(cols[3]),    # 3  CTG_END
                        cols[4],         # 4  CTG_DIR
                        cols[5],         # 5  CHR_NAM
                        int(cols[6]),    # 6  CHR_LEN
                        int(cols[7]),    # 7  CHR_STR
                        int(cols[8]),    # 8  CHR_END
                        int(cols[11]),   # 9  CTG_MAPQ
                        line_idx,        # 10 line idx in primary PAF (used for CTG_GLOBALIDX)
                    ]
                except (IndexError, ValueError):
                    continue
                chunks_per_contig[ctg_name].append(chunk)

        simple_alt_count = 0
        simple_alt_skip_existing_count = 0
        simple_alt_keep_diff_dir_count = 0
        for ctg_name, chunks in sorted(chunks_per_contig.items()):
            if len(chunks) < 2:
                continue
            chunks.sort(key=lambda c: (c[CTG_STR], c[CTG_END]))

            # 기존 contig들과 동일한 telo / repeat / censat 라벨을 계산해서 채움
            # (label_node / label_repeat_node 는 chunk[CHR_NAM/CHR_STR/CHR_END] 만 봄)
            telo_labels = label_node(chunks, telo_dict)
            repeat_labels = label_repeat_node(chunks, repeat_data, chr_len)
            censat_labels = label_repeat_node(chunks, repeat_censat_data, chr_len)

            trimmed_indices = trim_alt_simple_terminal_indices(chunks, telo_labels, repeat_labels, chr_len)
            if len(trimmed_indices) < 2:
                continue

            indexed_chunks = [(chunk_i, chunks[chunk_i]) for chunk_i in trimmed_indices]
            selected_chroms = select_alt_simple_major_chroms(indexed_chunks)
            if len(selected_chroms) < 2:
                continue

            same_dir_candidates, diff_dir_transitions = find_alt_simple_transition_candidates(indexed_chunks, selected_chroms)
            candidates_to_add = []
            for candidate in same_dir_candidates:
                if alt_simple_candidate_near_existing(candidate, existing_nclose_loci, ALT_SIMPLE_EXISTING_NCLOSE_DIST):
                    simple_alt_skip_existing_count += 1
                    continue
                candidates_to_add.append(candidate)
            candidates_to_add.extend(diff_dir_transitions)
            simple_alt_keep_diff_dir_count += len(diff_dir_transitions)

            if not candidates_to_add:
                continue

            candidate_serial = 1
            syn_name_base = f"simple_ctg_alt_{ctg_name}"
            for candidate in candidates_to_add:
                syn_name = syn_name_base if candidate_serial == 1 else f"{syn_name_base}_nclose{candidate_serial}"
                while syn_name in existing_names:
                    candidate_serial += 1
                    syn_name = f"{syn_name_base}_nclose{candidate_serial}"
                existing_names.add(syn_name)
                candidate_serial += 1

                sv_type = 1 if candidate[0][1][CHR_NAM] != candidate[1][1][CHR_NAM] else 2
                s_idx = len(real_final_contig)
                e_idx = s_idx + 1
                for chunk_i, raw in candidate:
                    line_idx = raw[10]                       # tracked when parsing primary PAF
                    syn_chunk = list(raw[:10]) + [
                        sv_type,                             # 10 CTG_TYP
                        s_idx,                               # 11 CTG_STRND
                        e_idx,                               # 12 CTG_ENDND
                        telo_labels[chunk_i][0],             # 13 CTG_TELCHR
                        telo_labels[chunk_i][1],             # 14 CTG_TELDIR
                        '0',                                 # 15 CTG_TELCON (telcon_set 미참여)
                        repeat_labels[chunk_i][0],           # 16 CTG_RPTCHR
                        repeat_labels[chunk_i][1],           # 17 CTG_RPTCASE
                        censat_labels[chunk_i][1],           # 18 CTG_CENSAT
                        '0',                                 # 19 CTG_MAINFLOWDIR (find_mainflow가 덮어씀)
                        '0',                                 # 20 CTG_MAINFLOWCHR (find_mainflow가 덮어씀)
                        f'0.{line_idx}',                     # 21 CTG_GLOBALIDX (paf_file[0][line_idx]로 11번이 lookup)
                    ]
                    syn_chunk[CTG_NAM] = syn_name
                    real_final_contig.append(syn_chunk)

                if candidate[0][1][CTG_DIR] == candidate[1][1][CTG_DIR]:
                    existing_nclose_loci.append(alt_simple_pair_loci(candidate))
                simple_alt_count += 1

        if simple_alt_count > 0:
            logging.info(f"Number of simple ctg alt nclose candidates added : {simple_alt_count}")
        if simple_alt_skip_existing_count > 0:
            logging.info(f"Number of simple ctg alt candidates skipped as existing nclose : {simple_alt_skip_existing_count}")
        if simple_alt_keep_diff_dir_count > 0:
            logging.info(f"Number of simple ctg alt direction-changing transitions kept unfiltered : {simple_alt_keep_diff_dir_count}")

    adjacency = initial_graph_build(real_final_contig, telo_bound_dict)

    telo_connected_node, telo_connected_dict, telo_connected_graph_dict, telo_coverage = edge_optimization(
        real_final_contig,
        adjacency,
        telo_bound_dict,
        asm2cov,
        excluded_telomere_origins,
    )
    telo_connected_node, telo_connected_dict, telo_connected_graph_dict = filter_telomere_connected_cen_fragment_mismatch(
        real_final_contig,
        telo_connected_graph_dict,
        cen_fragment_meta,
        "final",
    )

    telo_edges = [
        (telo_name, tuple(edge))
        for telo_name, edges in telo_connected_graph_dict.items()
        for edge in edges
    ]
    virtual_telomere_nodes = add_missing_virtual_telomeres(
        real_final_contig,
        telo_edges,
        chr_len,
        telo_data,
        repeat_censat_data,
    )
    if virtual_telomere_nodes:
        logging.info(
            f"Added {virtual_telomere_nodes} missing chromosome-end virtual telomere nodes"
        )

    write_telomere_connected_outputs(PREFIX, telo_edges, real_final_contig)

    real_final_mainflow = find_mainflow(real_final_contig)
    for i in real_final_contig:
        try:
            i[CTG_MAINFLOWDIR] = real_final_mainflow[i[CTG_NAM]][0]
            i[CTG_MAINFLOWCHR] = real_final_mainflow[i[CTG_NAM]][1]
        except:
            pass


    with open(PREPROCESSED_PAF_FILE_PATH, "wt") as f:
        for i in real_final_contig:
            for j in i:
                print(j, end="\t", file=f)
            print("", file=f)

    return telo_coverage

def get_corr_dir(is_for : bool, dir_str : str) -> str:
    assert(dir_str == '+' or dir_str == '-')

    if is_for:
        return dir_str
    else:
        if dir_str == '+':
            return '-'
        else:
            return '+'

def get_left_right_centromere(repeat_censat_data : dict, chr_len : dict):
    left_start_cent = {}
    right_end_cent = {}
    for chrom, intervals in repeat_censat_data.items():
        chrom_length = chr_len[chrom]
        rep_start_0, rep_end_0 = intervals[0]  # 0-indexed 좌표
        if rep_start_0 == 0:
            left_start_cent[chrom] = rep_end_0
        elif rep_end_0 == chrom_length - 1:
            right_end_cent[chrom] = rep_start_0

    return left_start_cent, right_end_cent

def similar_centromere_nclose_cluster(nclose_dict : dict, contig_data : list, repeat_censat_data : dict, chr_len : dict):
    nclose_nodes = set()
    for j in nclose_dict:
        for k in nclose_dict[j]:
            nclose_nodes.add(k)

    left_cent, right_cent = get_left_right_centromere(repeat_censat_data, chr_len)
    centromere_nclose_master = set()
    slave_dict = dict()
    for (s, e) in nclose_nodes:
        contig_s = contig_data[s]
        contig_e = contig_data[e]
        if contig_s[CHR_NAM] in left_cent and contig_s[CHR_STR] <= left_cent[contig_s[CHR_NAM]] and contig_s[CTG_DIR] == '+':
            chk = True
            for (sm, em) in centromere_nclose_master:
                contig_sm = contig_data[sm]
                contig_em = contig_data[em]
                if contig_sm[CHR_NAM] in left_cent and contig_e[CHR_NAM] == contig_em[CHR_NAM] \
                        and contig_s[CTG_DIR] == contig_sm[CTG_DIR] and contig_e[CTG_DIR] == contig_em[CTG_DIR] \
                        and distance_checker(contig_e, contig_em) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
                elif contig_em[CHR_NAM] in left_cent and contig_e[CHR_NAM] == contig_sm[CHR_NAM] \
                        and contig_s[CTG_DIR] != contig_em[CTG_DIR] and contig_e[CTG_DIR] != contig_sm[CTG_DIR] \
                        and distance_checker(contig_e, contig_sm) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
            if chk:
                centromere_nclose_master.add((s, e))

        elif contig_e[CHR_NAM] in left_cent and contig_e[CHR_STR] <= left_cent[contig_e[CHR_NAM]] and contig_e[CTG_DIR] == '-':
            chk = True
            for (sm, em) in centromere_nclose_master:
                contig_sm = contig_data[sm]
                contig_em = contig_data[em]
                if contig_em[CHR_NAM] in left_cent and contig_s[CHR_NAM] == contig_sm[CHR_NAM] \
                        and contig_s[CTG_DIR] == contig_sm[CTG_DIR] and contig_e[CTG_DIR] == contig_em[CTG_DIR] \
                        and distance_checker(contig_s, contig_sm) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
                elif contig_sm[CHR_NAM] in left_cent and contig_s[CHR_NAM] == contig_em[CHR_NAM] \
                        and contig_s[CTG_DIR] != contig_em[CTG_DIR] and contig_e[CTG_DIR] != contig_sm[CTG_DIR] \
                        and distance_checker(contig_s, contig_em) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
            if chk:
                centromere_nclose_master.add((s, e))

        elif contig_s[CHR_NAM] in right_cent and contig_s[CHR_END] >= right_cent[contig_s[CHR_NAM]] and contig_s[CTG_DIR] == '-':
            chk = True
            for (sm, em) in centromere_nclose_master:
                contig_sm = contig_data[sm]
                contig_em = contig_data[em]
                if contig_sm[CHR_NAM] in right_cent and contig_e[CHR_NAM] == contig_em[CHR_NAM] \
                        and contig_s[CTG_DIR] == contig_sm[CTG_DIR] and contig_e[CTG_DIR] == contig_em[CTG_DIR] \
                        and distance_checker(contig_e, contig_em) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
                elif contig_em[CHR_NAM] in right_cent and contig_e[CHR_NAM] == contig_sm[CHR_NAM] \
                        and contig_s[CTG_DIR] != contig_em[CTG_DIR] and contig_e[CTG_DIR] != contig_sm[CTG_DIR] \
                        and distance_checker(contig_e, contig_sm) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
            if chk:
                centromere_nclose_master.add((s, e))
        elif contig_e[CHR_NAM] in right_cent and contig_e[CHR_END] >= right_cent[contig_e[CHR_NAM]] and contig_e[CTG_DIR] == '+':
            chk = True
            for (sm, em) in centromere_nclose_master:
                contig_sm = contig_data[sm]
                contig_em = contig_data[em]
                if contig_em[CHR_NAM] in right_cent and contig_s[CHR_NAM] == contig_sm[CHR_NAM] \
                        and contig_s[CTG_DIR] == contig_sm[CTG_DIR] and contig_e[CTG_DIR] == contig_em[CTG_DIR] \
                        and distance_checker(contig_s, contig_sm) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
                elif contig_sm[CHR_NAM] in right_cent and contig_s[CHR_NAM] == contig_em[CHR_NAM] \
                        and contig_s[CTG_DIR] != contig_em[CTG_DIR] and contig_e[CTG_DIR] != contig_sm[CTG_DIR] \
                        and distance_checker(contig_s, contig_em) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break

            if chk:
                centromere_nclose_master.add((s, e))

    return centromere_nclose_master, slave_dict

def conjoined_type4(contig_data, type2_nclose_node):
    """
    Detect type4 insertions/deletions by pairing type2 nodes within each chromosome.
    If a contig's internal span >= 100 kb, also test its flipped orientation in a nested loop.
    """
    conjoined_type4_ins = set()
    conjoined_type4_del = set()

    def flip_dir(ctg):
        flipped = list(ctg)
        flipped[CTG_DIR] = '+' if ctg[CTG_DIR] == '-' else '-'
        return flipped

    for chrom, type2_list in type2_nclose_node.items():
        L = len(type2_list)
        for i in range(L):
            for j in range(i + 1, L):
                t2n1 = type2_list[i]
                t2n2 = type2_list[j]

                c1fi, c1bi = t2n1
                c2fi, c2bi = t2n2

                c1f = contig_data[c1fi]
                c1b = contig_data[c1bi]
                c2f = contig_data[c2fi]
                c2b = contig_data[c2bi]

                # Internal spans; switch to CHR_END if that's your convention
                len_c1 = abs(c1f[CHR_STR] - c1b[CHR_STR])
                len_c2 = abs(c2f[CHR_STR] - c2b[CHR_STR])

                # Candidate orientations per contig (original, and flipped if span >= threshold)
                c1_flips = [(c1f, c1fi, c1b, c1bi)]
                c2_flips = [(c2f, c2fi, c2b, c2bi)]

                if len_c1 >= TYPE2_DIST_FLIP_THRESHOLD:
                    c1_flips.append((flip_dir(c1b), c1bi, flip_dir(c1f), c1fi))
                if len_c2 >= TYPE2_DIST_FLIP_THRESHOLD:
                    c2_flips.append((flip_dir(c2b), c2bi, flip_dir(c2f), c2fi))

                # Test all combinations (original/flip × original/flip)
                for c1f_use, c1fi_, c1b_use, c1bi_ in c1_flips:
                    for c2f_use, c2fi_, c2b_use, c2bi_ in c2_flips:
                        # Forward direction
                        if c1b_use[CTG_DIR] == c2f_use[CTG_DIR]:
                            dist_for = distance_checker(c1b_use, c2f_use)
                            if dist_for is not None and dist_for < TYPE2_CONJOIN_COMPRESS_LIMIT:
                                ratio, total_ref_len = calculate_single_contig_ref_ratio([c1f_use, c2b_use])
                                abs_estimated = total_ref_len * abs(ratio)
                                if abs(ratio - 1) > BND_CONTIG_BOUND and abs_estimated > TYPE2_CHUKJI_AS_TYPE4:
                                    # Insertion
                                    if ratio < 0:
                                        if contig_data[c1fi_][CHR_STR] > contig_data[c2bi_][CHR_STR]:
                                            conjoined_type4_ins.add((c1fi_, c1bi_, c2fi_, c2bi_))
                                        else:
                                            conjoined_type4_ins.add((c2bi_, c2fi_, c1bi_, c1fi_))
                                    # Deletion
                                    else:
                                        if contig_data[c1fi_][CHR_STR] < contig_data[c2bi_][CHR_STR]:
                                            conjoined_type4_del.add((c1fi_, c1bi_, c2fi_, c2bi_))
                                        else:
                                            conjoined_type4_del.add((c2bi_, c2fi_, c1bi_, c1fi_))


    conjoined_type4_ins = sorted(conjoined_type4_ins)
    conjoined_type4_del = sorted(conjoined_type4_del)

    return conjoined_type4_ins, conjoined_type4_del

def get_breakend_coord(contig, side_idx):
    nclose_loc = side_idx == 0
    nclose_dir = contig[CTG_DIR] == '+'
    return contig[CHR_STR if nclose_dir ^ nclose_loc else CHR_END]

def get_expected_high_side_from_contig(contig, side_idx):
    nclose_loc = side_idx == 0
    nclose_dir = contig[CTG_DIR] == '+'
    return 'right' if nclose_dir ^ nclose_loc else 'left'

def get_raw_endpoint(contig, side_idx, contig_idx):
    nclose_loc = side_idx == 0
    is_front_dir = contig[CTG_DIR] == '+'
    return {
        'chrom': contig[CHR_NAM],
        'coord': int(get_breakend_coord(contig, side_idx)),
        'dir': contig[CTG_DIR],
        'match_forward': bool(is_front_dir == nclose_loc),
        'contig_idx': int(contig_idx),
        'ctg_name': contig[CTG_NAM],
        'ctg_st': int(contig[CTG_STR]),
        'ctg_nd': int(contig[CTG_END]),
        'ref_st': int(contig[CHR_STR]),
        'ref_nd': int(contig[CHR_END]),
        'expected_high_side': get_expected_high_side_from_contig(contig, side_idx),
    }

def canonical_raw_nclose_layout(contig_data, nclose_pair):
    ordered_endpoints = [
        get_raw_endpoint(contig_data[contig_idx], side_idx, contig_idx)
        for side_idx, contig_idx in enumerate(nclose_pair)
    ]

    order = [0, 1]
    keys = [(chr2int(ordered_endpoints[i]['chrom']), ordered_endpoints[i]['coord']) for i in order]
    if keys[1] < keys[0]:
        order = [1, 0]

    endpoints = [ordered_endpoints[i] for i in order]
    return {
        'nclose_key': tuple(sorted(int(x) for x in nclose_pair)),
        'chroms': tuple(ep['chrom'] for ep in endpoints),
        'coords': tuple(ep['coord'] for ep in endpoints),
        'sides': tuple(ep['expected_high_side'] for ep in endpoints),
        'endpoints': tuple(endpoints),
        'ordered_endpoints': tuple(ordered_endpoints),
    }

def get_inner_boundary_interval(region_a, region_b):
    st_a, nd_a = int(region_a['ref_st']), int(region_a['ref_nd'])
    st_b, nd_b = int(region_b['ref_st']), int(region_b['ref_nd'])

    if (st_b, nd_b) < (st_a, nd_a):
        st_a, nd_a, st_b, nd_b = st_b, nd_b, st_a, nd_a

    if nd_a <= st_b:
        return nd_a, st_b

    return st_b, min(nd_a, nd_b)

def get_outer_reference_interval(region_a, region_b):
    return (
        min(int(region_a['ref_st']), int(region_b['ref_st'])),
        max(int(region_a['ref_nd']), int(region_b['ref_nd'])),
    )

def centered_reference_interval(chrom, inner_st, inner_nd, chr_len, size=RAW_TRANSLOCATION_WINDOW):
    center = int(round((int(inner_st) + int(inner_nd)) / 2))
    chrom_len = int(chr_len.get(chrom, center + size))
    st = center - size // 2
    nd = st + size

    if st < 0:
        st = 0
        nd = min(size, chrom_len)
    elif nd > chrom_len:
        nd = chrom_len
        st = max(0, nd - size)

    return int(st), int(nd)

def layout_has_split_contig_endpoint(layout):
    return any(
        str(endpoint['ctg_name']).startswith('split_contig_')
        for endpoint in layout['ordered_endpoints']
    )

def layout_endpoint_dirs_are_opposite(layout):
    endpoints = layout['ordered_endpoints']
    return len(endpoints) == 2 and endpoints[0]['dir'] != endpoints[1]['dir']

def is_large_same_chrom_raw_candidate(chrom_pair, layout_a, layout_b):
    if chrom_pair[0] != chrom_pair[1]:
        return False
    if not layout_endpoint_dirs_are_opposite(layout_a):
        return False
    if not layout_endpoint_dirs_are_opposite(layout_b):
        return False
    if layout_has_split_contig_endpoint(layout_a) or layout_has_split_contig_endpoint(layout_b):
        return False

    span_a = abs(int(layout_a['coords'][1]) - int(layout_a['coords'][0]))
    span_b = abs(int(layout_b['coords'][1]) - int(layout_b['coords'][0]))
    return max(span_a, span_b) >= RAW_TRANSLOCATION_MIN_SAME_CHROM_SPAN

def node_fully_inside_censat(node, repeat_censat_data):
    for censat_st, censat_nd in repeat_censat_data.get(node[CHR_NAM], []):
        if censat_st <= node[CHR_STR] and node[CHR_END] <= censat_nd:
            return True
    return False

def extract_raw_censat_type2_candidates(paf_path, repeat_censat_data):
    candidates = []
    if paf_path is None or not os.path.exists(paf_path):
        return candidates

    chunks_per_contig = defaultdict(list)
    with open(paf_path, "rt") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            try:
                chunk = [
                    cols[0],
                    int(cols[1]),
                    int(cols[2]),
                    int(cols[3]),
                    cols[4],
                    cols[5],
                    int(cols[6]),
                    int(cols[7]),
                    int(cols[8]),
                    int(cols[11]),
                ]
            except (IndexError, ValueError):
                continue
            chunks_per_contig[chunk[CTG_NAM]].append(chunk)

    for ctg_name, chunks in chunks_per_contig.items():
        if len(chunks) < 2:
            continue
        chunks.sort(key=lambda c: (c[CTG_STR], c[CTG_END]))
        if any(chunk[CTG_MAPQ] <= 0 for chunk in chunks):
            continue
        if any(not node_fully_inside_censat(chunk, repeat_censat_data) for chunk in chunks):
            continue

        start = chunks[0]
        end = chunks[-1]
        if start[CHR_NAM] != end[CHR_NAM]:
            continue
        if start[CTG_DIR] == end[CTG_DIR]:
            continue

        candidates.append({
            "ctg_name": ctg_name,
            "chrom": start[CHR_NAM],
            "pair": (tuple(start), tuple(end)),
        })

    return candidates

def connected_type2_endpoint_for_template(template_contig, template_side, type2_pair):
    type2_front, type2_back = type2_pair
    if template_side == 0:
        if template_contig[CTG_DIR] == '+':
            dist = distance_checker(type2_back, template_contig)
            if type2_back[CHR_END] < template_contig[CHR_STR] or dist == 0:
                return type2_back, dist
        else:
            dist = distance_checker(type2_front, template_contig)
            if type2_front[CHR_END] > template_contig[CHR_STR] or dist == 0:
                return type2_front, dist
    else:
        if template_contig[CTG_DIR] == '+':
            dist = distance_checker(type2_front, template_contig)
            if type2_front[CHR_STR] > template_contig[CHR_END] or dist == 0:
                return type2_front, dist
        else:
            dist = distance_checker(type2_back, template_contig)
            if type2_back[CHR_STR] < template_contig[CHR_END] or dist == 0:
                return type2_back, dist
    return None, None

def append_ppc_rows(ppc_path, rows):
    if not rows:
        return
    with open(ppc_path, "at") as f:
        for row in rows:
            for value in row:
                print(value, end="\t", file=f)
            print("", file=f)

def collect_missing_cen_fragment_dir_censat_noncensat(
    contig_data,
    nclose_nodes,
    cen_fragment_meta,
):
    candidates = []
    cent_fragment_chroms = set(cen_fragment_meta.keys())
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            if str(ctg_name).startswith('simple_ctg_alt_') \
               or contig_data[pair[0]][CTG_NAM].startswith('simple_ctg_alt_') \
               or contig_data[pair[1]][CTG_NAM].startswith('simple_ctg_alt_'):
                continue
            s_is_censat = contig_data[pair[0]][CTG_CENSAT] != '0'
            e_is_censat = contig_data[pair[1]][CTG_CENSAT] != '0'
            if s_is_censat == e_is_censat:
                continue
            if s_is_censat:
                cidx, nidx = pair[0], pair[1]
                censat_norm_dir = contig_data[cidx][CTG_DIR]
            else:
                cidx, nidx = pair[1], pair[0]
                censat_norm_dir = _flip_ctg_dir(contig_data[cidx][CTG_DIR])

            censat_chr = contig_data[cidx][CHR_NAM]
            if censat_chr not in cent_fragment_chroms:
                continue
            noncensat_chr = contig_data[nidx][CHR_NAM]
            noncensat_pos = (contig_data[nidx][CHR_STR] + contig_data[nidx][CHR_END]) // 2
            candidates.append({
                "ctg_name": ctg_name,
                "pair": tuple(pair),
                "cidx": cidx,
                "nidx": nidx,
                "censat_side": 0 if s_is_censat else 1,
                "censat_chr": censat_chr,
                "noncensat_chr": noncensat_chr,
                "noncensat_pos": noncensat_pos,
                "censat_norm_dir": censat_norm_dir,
            })

    grouped = defaultdict(list)
    for cand in candidates:
        grouped[(cand["censat_chr"], cand["noncensat_chr"])].append(cand)

    missing = []
    for (censat_chr, _), items in grouped.items():
        target_norm_dir = _cen_fragment_target_dir_from_meta(cen_fragment_meta, censat_chr)
        items.sort(key=lambda x: x["noncensat_pos"])
        i = 0
        while i < len(items):
            j = i + 1
            while j < len(items) and items[j]["noncensat_pos"] - items[j-1]["noncensat_pos"] < OFFSET_DIR_GROUP_LIMIT:
                j += 1
            sub = items[i:j]
            if all(it["censat_norm_dir"] != target_norm_dir for it in sub):
                missing.append(sub)
            i = j

    return missing

def filter_offset_direction_mismatched_censat_noncensat(
    contig_data,
    nclose_nodes,
    cen_fragment_meta,
):
    offset_filter_candidates = []
    cent_fragment_chroms = set(cen_fragment_meta.keys())
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            s_is_censat = contig_data[pair[0]][CTG_CENSAT] != '0'
            e_is_censat = contig_data[pair[1]][CTG_CENSAT] != '0'
            if s_is_censat == e_is_censat:
                continue
            if s_is_censat:
                cidx, nidx = pair[0], pair[1]
                censat_norm_dir = contig_data[cidx][CTG_DIR]
            else:
                cidx, nidx = pair[1], pair[0]
                censat_norm_dir = _flip_ctg_dir(contig_data[cidx][CTG_DIR])
            censat_chr = contig_data[cidx][CHR_NAM]
            if censat_chr not in cent_fragment_chroms:
                continue
            noncensat_chr = contig_data[nidx][CHR_NAM]
            noncensat_pos = (contig_data[nidx][CHR_STR] + contig_data[nidx][CHR_END]) // 2
            offset_filter_candidates.append(
                (ctg_name, tuple(pair), censat_chr, noncensat_chr, noncensat_pos, censat_norm_dir)
            )

    offset_group_map = defaultdict(list)
    for cand in offset_filter_candidates:
        offset_group_map[(cand[2], cand[3])].append(cand)

    offset_to_remove = set()
    for (censat_chr, _), items in offset_group_map.items():
        target_norm_dir = _cen_fragment_target_dir_from_meta(cen_fragment_meta, censat_chr)
        items.sort(key=lambda x: x[4])
        i = 0
        while i < len(items):
            j = i + 1
            while j < len(items) and items[j][4] - items[j-1][4] < OFFSET_DIR_GROUP_LIMIT:
                j += 1
            sub = items[i:j]
            dirs = {it[5] for it in sub}
            if len(dirs) > 1:
                for it in sub:
                    if it[5] != target_norm_dir:
                        offset_to_remove.add((it[0], it[1]))
            i = j

    offset_filtered_nclose_nodes = defaultdict(list)
    offset_removed = 0
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            if (ctg_name, tuple(pair)) in offset_to_remove:
                offset_removed += 1
                continue
            offset_filtered_nclose_nodes[ctg_name].append(pair)

    return offset_filtered_nclose_nodes, offset_removed

def add_nearest_combined_censat_noncensat_ncloses(
    contig_data,
    nclose_nodes,
    missing_candidate_groups,
    raw_censat_type2_candidates,
    repeat_censat_data,
    chr_len,
    cen_fragment_meta,
    ppc_path,
):
    if not missing_candidate_groups or not raw_censat_type2_candidates:
        return 0

    type2_by_chrom = defaultdict(list)
    for candidate in raw_censat_type2_candidates:
        type2_by_chrom[candidate["chrom"]].append(candidate)

    rows_to_append = []
    seen_signature = set()
    added = 0

    for group in missing_candidate_groups:
        best = None
        for cand in group:
            target_dir = _cen_fragment_target_dir_from_meta(cen_fragment_meta, cand["censat_chr"])
            censat_node = contig_data[cand["cidx"]]
            for type2 in type2_by_chrom.get(cand["censat_chr"], []):
                type2_endpoint, dist = connected_type2_endpoint_for_template(
                    censat_node,
                    cand["censat_side"],
                    type2["pair"],
                )
                if type2_endpoint is None:
                    continue
                if type2_endpoint[CTG_DIR] != target_dir:
                    continue
                best_key = (
                    dist,
                    abs(censat_node[CHR_STR] - type2_endpoint[CHR_STR]),
                    cand["noncensat_pos"],
                )
                if best is None or best_key < best[0]:
                    best = (best_key, cand, type2, type2_endpoint)

        if best is None:
            continue

        _, cand, type2, remaining_endpoint = best
        noncensat_node = contig_data[cand["nidx"]]
        # The synthetic contig is written as censat -> noncensat. If the
        # original nclose had censat on the second side, preserve the
        # noncensat endpoint in that reverse-complemented frame.
        noncensat_dir = noncensat_node[CTG_DIR]
        if cand["censat_side"] == 1:
            noncensat_dir = _flip_ctg_dir(noncensat_dir)
        signature = (
            remaining_endpoint[CTG_NAM],
            remaining_endpoint[CHR_NAM],
            remaining_endpoint[CHR_STR],
            remaining_endpoint[CHR_END],
            remaining_endpoint[CTG_DIR],
            cand["nidx"],
            noncensat_dir,
        )
        if signature in seen_signature:
            continue
        seen_signature.add(signature)

        new_name = f"combined_{type2['ctg_name']}_{noncensat_node[CTG_NAM]}_{added}"
        new_idx0 = len(contig_data)
        new_idx1 = new_idx0 + 1
        len0 = remaining_endpoint[CHR_END] - remaining_endpoint[CHR_STR]
        len1 = noncensat_node[CHR_END] - noncensat_node[CHR_STR]
        total_len = len0 + len1
        sv_type = 1 if remaining_endpoint[CHR_NAM] != noncensat_node[CHR_NAM] else 2
        censat_label = label_repeat_node([remaining_endpoint], repeat_censat_data, chr_len)[0][1]

        node0 = [
            new_name,
            total_len,
            0,
            len0,
            remaining_endpoint[CTG_DIR],
            remaining_endpoint[CHR_NAM],
            remaining_endpoint[CHR_LEN],
            remaining_endpoint[CHR_STR],
            remaining_endpoint[CHR_END],
            remaining_endpoint[CTG_MAPQ],
            sv_type,
            new_idx0,
            new_idx1,
            '0',
            '0',
            '0',
            '0',
            '0',
            censat_label,
            remaining_endpoint[CTG_DIR],
            remaining_endpoint[CHR_NAM],
            f'4.{2 * added}',
        ]
        node1 = [
            new_name,
            total_len,
            len0,
            total_len,
            noncensat_dir,
            noncensat_node[CHR_NAM],
            noncensat_node[CHR_LEN],
            noncensat_node[CHR_STR],
            noncensat_node[CHR_END],
            noncensat_node[CTG_MAPQ],
            sv_type,
            new_idx0,
            new_idx1,
            '0',
            '0',
            '0',
            noncensat_node[CTG_RPTCHR],
            noncensat_node[CTG_RPTCASE],
            noncensat_node[CTG_CENSAT],
            noncensat_dir,
            noncensat_node[CHR_NAM],
            f'4.{2 * added + 1}',
        ]

        contig_data.append(tuple(node0))
        contig_data.append(tuple(node1))
        rows_to_append.extend([node0, node1])
        nclose_nodes[new_name].append((new_idx0, new_idx1))
        added += 1

    append_ppc_rows(ppc_path, rows_to_append)
    return added

def build_raw_translocation_candidates(contig_data, nclose_nodes, chr_len,
                                       distance_threshold=RAW_TRANSLOCATION_WINDOW,
                                       candidate_filter=None):
    layout_by_key = {}
    for _, pair_list in nclose_nodes.items():
        for pair in pair_list:
            nclose_key = tuple(sorted(int(x) for x in pair))
            if len(nclose_key) != 2:
                continue
            layout_by_key[nclose_key] = canonical_raw_nclose_layout(contig_data, tuple(int(x) for x in pair))

    groups = defaultdict(list)
    for nclose_key, layout in layout_by_key.items():
        groups[layout['chroms']].append((nclose_key, layout))

    candidates = []
    seen_candidate_signatures = set()
    pair_id = 0
    for chrom_pair, items in groups.items():
        for i in range(len(items)):
            key_i, layout_i = items[i]
            for j in range(i + 1, len(items)):
                key_j, layout_j = items[j]
                if layout_i['sides'][0] == layout_j['sides'][0] or layout_i['sides'][1] == layout_j['sides'][1]:
                    continue
                if abs(layout_i['coords'][0] - layout_j['coords'][0]) > distance_threshold:
                    continue
                if abs(layout_i['coords'][1] - layout_j['coords'][1]) > distance_threshold:
                    continue

                if layout_i['sides'][0] == 'left':
                    key_a, layout_a = key_i, layout_i
                    key_b, layout_b = key_j, layout_j
                else:
                    key_a, layout_a = key_j, layout_j
                    key_b, layout_b = key_i, layout_i

                if candidate_filter is not None and not candidate_filter(chrom_pair, layout_a, layout_b):
                    continue

                candidate_signature = (
                    chrom_pair,
                    layout_a['coords'],
                    layout_a['sides'],
                    layout_b['coords'],
                    layout_b['sides'],
                )
                if candidate_signature in seen_candidate_signatures:
                    continue
                seen_candidate_signatures.add(candidate_signature)

                side_records = []
                for side_idx, chrom in enumerate(chrom_pair):
                    inner_st, inner_nd = get_inner_boundary_interval(
                        layout_a['endpoints'][side_idx],
                        layout_b['endpoints'][side_idx],
                    )
                    ref_st, ref_nd = centered_reference_interval(chrom, inner_st, inner_nd, chr_len)
                    path_drop_st, path_drop_nd = get_outer_reference_interval(
                        layout_a['endpoints'][side_idx],
                        layout_b['endpoints'][side_idx],
                    )
                    side_records.append({
                        'side_id': f'R{pair_id}.{side_idx}',
                        'chrom': chrom,
                        'inner_st': int(inner_st),
                        'inner_nd': int(inner_nd),
                        'ref_count_st': int(ref_st),
                        'ref_count_nd': int(ref_nd),
                        'path_drop_st': int(path_drop_st),
                        'path_drop_nd': int(path_drop_nd),
                        'no_spanning_rawread': False,
                        'no_spanning_utg': False,
                        'crossing_cols': [],
                        'accepted_forbid': False,
                    })

                candidates.append({
                    'pair_id': pair_id,
                    'nclose_key_a': key_a,
                    'nclose_key_b': key_b,
                    'chrom_pair': chrom_pair,
                    'layout_a': layout_a,
                    'layout_b': layout_b,
                    'split_junctions': [
                        {'count_name': 'd1', 'count_idx': 1, 'nclose_key': key_a, 'endpoints': layout_a['ordered_endpoints']},
                        {'count_name': 'd4', 'count_idx': 4, 'nclose_key': key_b, 'endpoints': layout_b['ordered_endpoints']},
                    ],
                    'span_junctions': [
                        {'count_name': 'd2', 'count_idx': 2, 'side_idx': 0, 'side_id': side_records[0]['side_id']},
                        {'count_name': 'd3', 'count_idx': 3, 'side_idx': 1, 'side_id': side_records[1]['side_id']},
                    ],
                    'side_records': side_records,
                    'pair_protected_from_1M_drop': False,
                })
                pair_id += 1

    return candidates

def raw_translocation_candidate_signature(candidate):
    return (
        candidate['chrom_pair'],
        candidate['layout_a']['coords'],
        candidate['layout_a']['sides'],
        candidate['layout_b']['coords'],
        candidate['layout_b']['sides'],
    )

def renumber_raw_translocation_candidates(candidates):
    for pair_id, candidate in enumerate(candidates):
        candidate['pair_id'] = pair_id
        for side_idx, side_record in enumerate(candidate['side_records']):
            side_record['side_id'] = f'R{pair_id}.{side_idx}'
        for span_junction in candidate['span_junctions']:
            side_idx = span_junction.get('side_idx')
            if side_idx is not None and 0 <= side_idx < len(candidate['side_records']):
                span_junction['side_id'] = candidate['side_records'][side_idx]['side_id']
    return candidates

def merge_raw_translocation_candidates(*candidate_lists):
    merged = []
    seen = set()
    for candidate_list in candidate_lists:
        for candidate in candidate_list:
            signature = raw_translocation_candidate_signature(candidate)
            if signature in seen:
                continue
            seen.add(signature)
            merged.append(candidate)
    return renumber_raw_translocation_candidates(merged)

def build_nclose_count_candidates(contig_data, nclose_nodes):
    candidates = []
    seen = set()
    pair_id = 0
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            nclose_key = tuple(sorted(int(x) for x in pair))
            if len(nclose_key) != 2:
                continue
            if nclose_key in seen:
                continue
            seen.add(nclose_key)
            layout = canonical_raw_nclose_layout(contig_data, tuple(int(x) for x in pair))
            overlaps_censat = any(contig_data[int(idx)][CTG_CENSAT] != '0' for idx in nclose_key)
            candidates.append({
                'pair_id': pair_id,
                'nclose_key': nclose_key,
                'ctg_name': ctg_name,
                'layout': layout,
                'overlaps_censat': overlaps_censat,
            })
            pair_id += 1
    return candidates

def collect_nclose_count_keep_pairs(prefix, vaf_threshold):
    result_path = f'{prefix}/{NCLOSE_COUNT_RESULT_PKL}'
    if not os.path.isfile(result_path):
        logging.warning(f'NClose raw-count result not found: {result_path}')
        return None, 0, 0

    with open(result_path, 'rb') as f:
        records = pkl.load(f)

    keep_pairs = set()
    for record in records:
        if 'keep_nclose' in record:
            if record['keep_nclose']:
                keep_pairs.add(tuple(sorted(int(x) for x in record['nclose_key'])))
            continue
        if record.get('filter_eligible') is False:
            keep_pairs.add(tuple(sorted(int(x) for x in record['nclose_key'])))
            continue
        vaf = record.get('vaf', {})
        chr_a_vaf = vaf.get('chr_a')
        chr_b_vaf = vaf.get('chr_b')
        if (
            (chr_a_vaf is not None and chr_a_vaf >= vaf_threshold) or
            (chr_b_vaf is not None and chr_b_vaf >= vaf_threshold)
        ):
            keep_pairs.add(tuple(sorted(int(x) for x in record['nclose_key'])))

    return keep_pairs, len(records), len(keep_pairs)

def filter_nclose_nodes_by_keep_pairs(nclose_nodes, keep_pairs):
    filtered = defaultdict(list)
    removed = 0
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            if tuple(sorted(int(x) for x in pair)) not in keep_pairs:
                removed += 1
                continue
            filtered[ctg_name].append(pair)
    return filtered, removed

def raw_false_value(value):
    return value is False or value == 0 or value == 'False' or value == 'false'

def raw_record_has_both_point_spans(record):
    no_span = record.get('raw_point_no_spanning')
    if isinstance(no_span, dict):
        return raw_false_value(no_span.get('point_a')) and raw_false_value(no_span.get('point_b'))

    side_records = record.get('side_records', [])
    if len(side_records) < 2:
        return False
    side_flags = [
        side.get('raw_point_no_spanning', side.get('no_spanning_rawread'))
        for side in side_records[:2]
    ]
    return raw_false_value(side_flags[0]) and raw_false_value(side_flags[1])

def raw_record_is_depth_balanced(record):
    if 'depth_balanced_translocation' not in record:
        return True
    return bool(record.get('depth_balanced_translocation'))

def collect_raw_virtual_inv_nclose_pairs(prefix):
    result_path = f'{prefix}/{RAW_TRANSLOCATION_RESULT_PKL}'
    if not os.path.isfile(result_path):
        return set()

    with open(result_path, 'rb') as f:
        records = pkl.load(f)

    pairs_to_remove = set()
    for record in records:
        if not raw_record_has_both_point_spans(record):
            continue
        if not raw_record_is_depth_balanced(record):
            continue
        for key_name in ('nclose_key_a', 'nclose_key_b'):
            key = record.get(key_name)
            if key is None:
                continue
            pairs_to_remove.add(tuple(sorted(int(x) for x in key)))

    return pairs_to_remove

def remove_nclose_pairs(nclose_nodes, pairs_to_remove):
    filtered = defaultdict(list)
    removed = 0
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            if tuple(sorted(int(x) for x in pair)) in pairs_to_remove:
                removed += 1
                continue
            filtered[ctg_name].append(pair)
    return filtered, removed

def collect_nclose_breakends(contig_data, nclose_nodes):
    breakends_by_chrom = defaultdict(list)
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            pair_key = tuple(pair)
            for side_idx, contig_idx in enumerate(pair):
                contig = contig_data[contig_idx]
                breakends_by_chrom[contig[CHR_NAM]].append((
                    get_breakend_coord(contig, side_idx),
                    pair_key,
                    ctg_name,
                ))

    for chrom in breakends_by_chrom:
        breakends_by_chrom[chrom].sort(key=lambda x: x[0])

    return breakends_by_chrom

def get_graph_edge_weight(contig_data, src_idx, dst_idx):
    contig_s = contig_data[src_idx]
    contig_e = contig_data[dst_idx]
    if contig_s[CTG_NAM] == contig_e[CTG_NAM]:
        if src_idx < dst_idx:
            idx_range = range(src_idx + 1, dst_idx)
        else:
            idx_range = range(dst_idx + 1, src_idx)
        return sum(
            contig_data[i][CHR_END] - contig_data[i][CHR_STR]
            for i in idx_range
        )

    ds = contig_s[CHR_END] - contig_s[CHR_STR]
    de = contig_e[CHR_END] - contig_e[CHR_STR]
    if distance_checker(contig_s, contig_e) == 0:
        return ds + de - overlap_calculator(contig_s, contig_e)
    return distance_checker(contig_s, contig_e) + ds + de

def get_type4_breakend_span(contig_slice):
    return (
        get_breakend_coord(contig_slice[0], 0),
        get_breakend_coord(contig_slice[-1], 1),
    )

def get_type4_outer_ref_span(contig_slice):
    if contig_slice[0][CTG_DIR] == '+':
        return contig_slice[0][CHR_STR], contig_slice[-1][CHR_END]
    return contig_slice[-1][CHR_STR], contig_slice[0][CHR_END]

def point_overlaps_censat(censat_dict, chrom, coord):
    for st, nd in censat_dict.get(chrom, []):
        if max(st, coord) < min(nd, coord + 1):
            return True
    return False

def span_overlaps_censat(censat_dict, chrom, span_st, span_nd):
    for st, nd in censat_dict.get(chrom, []):
        if max(st, span_st) < min(nd, span_nd):
            return True
    return False

def make_depth_by_chrom(depth_df):
    return {
        chrom: chrom_df
        for chrom, chrom_df in depth_df.groupby('chr', sort=False)
    }

def weighted_depth_mean(depth_by_chrom, chrom, start, end):
    if end <= start:
        return None
    chrom_df = depth_by_chrom.get(chrom)
    if chrom_df is None:
        return None
    mask = (chrom_df['nd'] > start) & (chrom_df['st'] < end)
    if not mask.any():
        return None

    rows = chrom_df.loc[mask]
    overlap_st = np.maximum(rows['st'].to_numpy(), start)
    overlap_nd = np.minimum(rows['nd'].to_numpy(), end)
    weights = overlap_nd - overlap_st
    if weights.sum() <= 0:
        return None
    return float(np.average(rows['meandepth'].to_numpy(), weights=weights))

def expected_depth_step(depth_by_chrom, chrom, coord, expected_high_side,
                        window=TYPE4_INDEL_GRAPH_DEPTH_WINDOW):
    left_mean = weighted_depth_mean(depth_by_chrom, chrom, coord - window, coord)
    right_mean = weighted_depth_mean(depth_by_chrom, chrom, coord, coord + window)
    if left_mean is None or right_mean is None:
        return None
    if expected_high_side == 'left':
        return left_mean - right_mean
    if expected_high_side == 'right':
        return right_mean - left_mean
    return None

def type4_indel_expected_high_sides(indel_kind):
    if indel_kind == 'deletion':
        return ('left', 'right')
    if indel_kind == 'insertion':
        return ('right', 'left')
    raise ValueError(f'Unknown type4 indel kind: {indel_kind}')

def type4_indel_depth_supported(depth_by_chrom, censat_dict, chrom, span_st, span_nd,
                                indel_kind, normal_haploid_depth):
    if point_overlaps_censat(censat_dict, chrom, span_st) or \
       point_overlaps_censat(censat_dict, chrom, span_nd):
        return False, None, 'endpoint_censat'

    expected_high_sides = type4_indel_expected_high_sides(indel_kind)
    steps = [
        expected_depth_step(depth_by_chrom, chrom, span_st, expected_high_sides[0]),
        expected_depth_step(depth_by_chrom, chrom, span_nd, expected_high_sides[1]),
    ]
    if any(step is None for step in steps):
        return False, steps, 'missing_depth_window'

    min_step = TYPE4_INDEL_GRAPH_DEPTH_DIFF_RATIO * normal_haploid_depth
    if all(step >= min_step for step in steps):
        return True, steps, 'pass'
    return False, steps, 'depth_direction_or_size'

def collect_type4_indel_graph_candidates(contig_data, depth_df, censat_dict,
                                         min_span=TYPE4_INDEL_GRAPH_MIN_SPAN):
    candidates = []
    depth_by_chrom = make_depth_by_chrom(depth_df)
    normal_haploid_depth = float(np.median(depth_df['meandepth'])) / 2
    s = 0
    contig_data_size = len(contig_data)
    while s < contig_data_size:
        e = contig_data[s][CTG_ENDND]
        contig_slice = contig_data[s:e + 1]
        contig_s = contig_data[s]
        contig_e = contig_data[e]

        if contig_s[CTG_TYP] == 4 \
            and contig_s[CHR_NAM] == contig_e[CHR_NAM] \
           and contig_s[CTG_DIR] == contig_e[CTG_DIR]:
            ratio, total_ref_len = calculate_single_contig_ref_ratio(contig_slice)
            ref_a, ref_b = get_type4_breakend_span(contig_slice)
            outer_ref_a, outer_ref_b = get_type4_outer_ref_span(contig_slice)
            span_st, span_nd = sorted((ref_a, ref_b))
            outer_span_st, outer_span_nd = sorted((outer_ref_a, outer_ref_b))
            span_len = span_nd - span_st

            indel_kind = None
            if ratio < 0:
                indel_kind = 'insertion'
            elif ratio > 0:
                indel_kind = 'deletion'

            if indel_kind is not None and span_len >= min_span:
                censat_jump = span_overlaps_censat(censat_dict, contig_s[CHR_NAM], span_st, span_nd)
                if indel_kind == 'deletion' and not censat_jump:
                    s = e + 1
                    continue

                depth_pass, depth_steps, depth_filter_reason = type4_indel_depth_supported(
                    depth_by_chrom, censat_dict, contig_s[CHR_NAM], span_st, span_nd,
                    indel_kind, normal_haploid_depth
                )
                if not depth_pass:
                    s = e + 1
                    continue

                candidates.append({
                    'type4_tuple': (s, e),
                    'src': (DIR_FOR, s),
                    'dst': (DIR_FOR, e),
                    'type4_kind': f'raw_type4_{indel_kind}',
                    'indel_kind': indel_kind,
                    'chrom': contig_s[CHR_NAM],
                    'base_chrom': contig_s[CHR_NAM],
                    'base_st': span_st,
                    'base_nd': span_nd,
                    'outer_base_st': outer_span_st,
                    'outer_base_nd': outer_span_nd,
                    'span_st': span_st,
                    'span_nd': span_nd,
                    'span_len': span_len,
                    'ratio': ratio,
                    'total_ref_len': total_ref_len,
                    'normal_haploid_depth': normal_haploid_depth,
                    'expected_high_sides': type4_indel_expected_high_sides(indel_kind),
                    'depth_steps': tuple(depth_steps),
                    'depth_min_step': TYPE4_INDEL_GRAPH_DEPTH_DIFF_RATIO * normal_haploid_depth,
                    'depth_filter_reason': depth_filter_reason,
                    'censat_jump': censat_jump,
                    'contig_name': contig_s[CTG_NAM],
                })

        s = e + 1

    return candidates

def select_type4_indel_graph_edges(contig_data, nclose_nodes, type4_indel_candidates):
    breakends_by_chrom = collect_nclose_breakends(contig_data, nclose_nodes)
    selected_edges = []

    for type4_idx, candidate in enumerate(type4_indel_candidates):
        s, e = candidate['type4_tuple']
        chrom = candidate['chrom']
        span_st = candidate['span_st']
        span_nd = candidate['span_nd']
        inside_nclose_pairs = {
            pair_key
            for bp, pair_key, _ in breakends_by_chrom.get(chrom, [])
            if span_st <= bp <= span_nd
        }
        if candidate['indel_kind'] == 'insertion' and len(inside_nclose_pairs) < 2:
            continue

        edge_specs = [
            ((DIR_FOR, s), (DIR_FOR, e)),
            ((DIR_BAK, e), (DIR_BAK, s)),
        ]
        for src, dst in edge_specs:
            selected_edges.append({
                'type4_kind': candidate['type4_kind'],
                'type4_idx': type4_idx,
                'type4_tuple': tuple(candidate['type4_tuple']),
                'src': src,
                'dst': dst,
                'dist': get_graph_edge_weight(contig_data, src[1], dst[1]),
                'indel_kind': candidate['indel_kind'],
                'chrom': chrom,
                'base_chrom': candidate['base_chrom'],
                'base_st': candidate['base_st'],
                'base_nd': candidate['base_nd'],
                'outer_base_st': candidate['outer_base_st'],
                'outer_base_nd': candidate['outer_base_nd'],
                'span_st': span_st,
                'span_nd': span_nd,
                'span_len': candidate['span_len'],
                'inside_nclose_count': len(inside_nclose_pairs),
                'ratio': candidate['ratio'],
                'expected_high_sides': candidate['expected_high_sides'],
                'depth_steps': candidate['depth_steps'],
                'depth_min_step': candidate['depth_min_step'],
                'censat_jump': candidate['censat_jump'],
                'contig_name': candidate['contig_name'],
            })

    return selected_edges

def augment_nclose_nodes_with_type4_indels(nclose_nodes, selected_edges):
    augmented_nclose_nodes = defaultdict(list)
    existing_pairs = defaultdict(set)
    for key, pair_list in nclose_nodes.items():
        augmented_nclose_nodes[key].extend(pair_list)
        existing_pairs[key].update(tuple(pair) for pair in pair_list)

    for edge in selected_edges:
        pair = tuple(edge['type4_tuple'])
        key = edge['contig_name']
        if pair in existing_pairs[key]:
            continue
        augmented_nclose_nodes[key].append(pair)
        existing_pairs[key].add(pair)

    return augmented_nclose_nodes

def get_type4_indel_zero_dim_edge_set(selected_edges):
    zero_dim_edge_set = set()
    for edge in selected_edges:
        src = tuple(edge['src'])
        dst = tuple(edge['dst'])
        zero_dim_edge_set.add((src, dst))

    return zero_dim_edge_set


VCF_SYNTHETIC_FLANK = 1 * K
VCF_SYNTHETIC_PAF_NAME = "vcf_synthetic.paf"
VCF_SKIPPED_RECORDS_TSV = "vcf_mode_skipped_records.tsv"
VCF_ORIENTATION_MISMATCH_TSV = "vcf_mode_orientation_mismatches.tsv"
VCF_TELOMERE_PAF_NODES_TSV = "vcf_telomere_paf_nodes.tsv"
LIMIT_COMBINATIONS_JSON = "limit_combinations.json"


def load_limit_combinations(path):
    try:
        with open(path, "rt", encoding="utf-8") as f:
            data = json.load(f)
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


def write_limit_combinations(path, combination):
    temporary = f"{path}.tmp-{os.getpid()}"
    with open(temporary, "wt", encoding="utf-8") as f:
        json.dump({"limit_combinations": list(combination)}, f, indent=2)
        f.write("\n")
    os.replace(temporary, path)


def load_vcf_ins_alt_alignment_spans(paf_path):
    best_row_by_query = {}
    best_key_by_query = {}
    rows_by_query = defaultdict(list)
    if not paf_path:
        return {}
    if not os.path.isfile(paf_path):
        logging.warning(f"VCF INS --alt PAF does not exist: {paf_path}")
        return {}

    with open(paf_path, "rt") as paf:
        for line in paf:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            try:
                qname = cols[0]
                chrom = cols[5]
                st, nd = sorted((int(cols[7]), int(cols[8])))
                matches = int(cols[9])
                block_len = int(cols[10])
                mapq = int(cols[11])
            except (IndexError, ValueError):
                continue
            if nd <= st:
                continue
            rows_by_query[qname].append(cols)
            # Keep every row for the INS type4 base PAF, but derive the
            # required single-chromosome representative coordinates from the
            # largest individual alignment block only.
            row_key = (block_len, mapq, matches, nd - st)
            if qname not in best_key_by_query or row_key > best_key_by_query[qname]:
                best_key_by_query[qname] = row_key
                best_row_by_query[qname] = {
                    "chrom": chrom,
                    "st": st,
                    "nd": nd,
                    "score": block_len,
                    "rows": 1,
                    "mapq": mapq,
                }

    best_by_query = {}
    for qname, best in best_row_by_query.items():
        if best["chrom"] is not None and best["nd"] > best["st"]:
            best = dict(best)
            best["paf_rows"] = rows_by_query.get(qname, [])
            best_by_query[qname] = best
    return best_by_query


def endpoint_interval(pos, chrom_len, path_dir, endpoint_role):
    pos = max(1, min(int(pos), int(chrom_len)))
    flank = min(VCF_SYNTHETIC_FLANK, max(1, int(chrom_len)))
    if endpoint_role == "exit":
        if path_dir == "+":
            ref_end = pos
            ref_st = max(0, ref_end - flank)
        else:
            ref_st = pos
            ref_end = min(int(chrom_len), ref_st + flank)
    else:
        if path_dir == "+":
            ref_st = pos
            ref_end = min(int(chrom_len), ref_st + flank)
        else:
            ref_end = pos
            ref_st = max(0, ref_end - flank)

    if ref_end <= ref_st:
        if ref_st <= 0:
            ref_st, ref_end = 0, 1
        else:
            ref_st, ref_end = ref_st - 1, ref_st
    return int(ref_st), int(ref_end)


def make_synthetic_vcf_node(contig_name, chrom, pos, path_dir, endpoint_role,
                            chr_len, ctg_typ, idx, global_idx):
    ref_st, ref_end = endpoint_interval(pos, chr_len[chrom], path_dir, endpoint_role)
    ref_len = max(1, ref_end - ref_st)
    return [
        contig_name,
        ref_len,
        0,
        ref_len,
        path_dir,
        chrom,
        int(chr_len[chrom]),
        ref_st,
        ref_end,
        60,
        ctg_typ,
        idx,
        idx,
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        path_dir,
        chrom,
        f"2.{global_idx}",
    ]


def make_synthetic_span_node(contig_name, chrom, st, nd, path_dir,
                             chr_len, ctg_typ, idx, global_idx):
    st, nd = sorted((int(st), int(nd)))
    if nd <= st:
        nd = st + 1
    st = max(0, min(st, int(chr_len[chrom]) - 1))
    nd = max(st + 1, min(nd, int(chr_len[chrom])))
    ref_len = nd - st
    return [
        contig_name,
        ref_len,
        0,
        ref_len,
        path_dir,
        chrom,
        int(chr_len[chrom]),
        st,
        nd,
        60,
        ctg_typ,
        idx,
        idx,
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        path_dir,
        chrom,
        f"2.{global_idx}",
    ]


def synthetic_node_to_paf_row(node):
    ref_len = max(1, int(node[CHR_END]) - int(node[CHR_STR]))
    return list(node[:9]) + [ref_len, ref_len, int(node[CTG_MAPQ]), "tp:A:P", f"cs:Z::{ref_len}"]


def telo_data_to_dict(telo_data):
    telo_dict = defaultdict(list)
    for entry in telo_data:
        telo_dict[entry[0]].append(entry[1:])
    return telo_dict


def telomere_label_to_name(label):
    if label[0] == "0" or label[1] == "0":
        return None
    side = str(label[1])[0]
    if side not in {"f", "b"}:
        return None
    return f"{label[0]}{side}"


def no_repeat_labels(contig_data):
    return [("0", "0") for _ in contig_data]


def telomere_bounds_by_name(telo_dict):
    telo_bounds = {}
    for chrom, intervals in telo_dict.items():
        intervals_by_side = defaultdict(list)
        for st, nd, side in intervals:
            intervals_by_side[side].append((int(st), int(nd)))
        for side, side_intervals in intervals_by_side.items():
            if side == "f":
                telo_bounds[f"{chrom}{side}"] = min(
                    side_intervals,
                    key=lambda interval: interval[0],
                )
            else:
                telo_bounds[f"{chrom}{side}"] = max(
                    side_intervals,
                    key=lambda interval: interval[1],
                )
    return telo_bounds


def is_terminal_telomere_candidate(row, telo_name, telo_bounds, chr_len):
    if telo_name.endswith("f"):
        return int(row[CHR_STR]) <= TELOMERE_CLUSTER_THRESHOLD

    telo_bound = telo_bounds.get(telo_name)
    if telo_bound is not None:
        telo_end = int(telo_bound[1])
    else:
        telo_end = int(chr_len[telo_name[:-1]])
    return int(row[CHR_END]) >= telo_end - TELOMERE_CLUSTER_THRESHOLD


def compress_paf_telomere_candidates(raw_rows, telo_connect_info, telo_dict,
                                     chr_len):
    """Cluster telomere anchors derived only from the assembly PAF."""
    telo_bounds = telomere_bounds_by_name(telo_dict)
    candidates_by_telo = defaultdict(list)
    for raw_idx, telo_name in telo_connect_info.items():
        candidates_by_telo[telo_name].append(raw_idx)

    compressed_candidates = []
    for telo_name, raw_indices in candidates_by_telo.items():
        nonterminal_indices = []
        terminal_idx = None
        for raw_idx in raw_indices:
            row = raw_rows[raw_idx]
            if is_terminal_telomere_candidate(
                row,
                telo_name,
                telo_bounds,
                chr_len,
            ):
                if terminal_idx is None:
                    terminal_idx = raw_idx
                    continue

                terminal_row = raw_rows[terminal_idx]
                if telo_name.endswith("f"):
                    is_better = (
                        int(row[CHR_STR]) < int(terminal_row[CHR_STR])
                        or (
                            int(row[CHR_STR]) == int(terminal_row[CHR_STR])
                            and int(row[CTG_LEN]) < int(terminal_row[CTG_LEN])
                        )
                    )
                else:
                    is_better = (
                        int(row[CHR_END]) > int(terminal_row[CHR_END])
                        or (
                            int(row[CHR_END]) == int(terminal_row[CHR_END])
                            and int(row[CTG_LEN]) < int(terminal_row[CTG_LEN])
                        )
                    )
                if is_better:
                    terminal_idx = raw_idx
                continue

            if any(
                distance_checker(raw_rows[kept_idx], row)
                < TELOMERE_COMPRESS_RANGE
                for kept_idx in nonterminal_indices
            ):
                continue
            nonterminal_indices.append(raw_idx)

        compressed_candidates.extend(
            (telo_name, raw_idx)
            for raw_idx in nonterminal_indices
        )
        if terminal_idx is not None:
            compressed_candidates.append((telo_name, terminal_idx))

    # break_double_telomere_contig prevents an original contig name from being
    # emitted for more than one final telomere edge.  This PAF-only path keeps
    # the strongest boundary anchor instead of manufacturing split contig paths.
    candidates_by_contig = defaultdict(list)
    for telo_name, raw_idx in compressed_candidates:
        candidates_by_contig[raw_rows[raw_idx][CTG_NAM]].append(
            (telo_name, raw_idx)
        )

    selected_candidates = set()
    duplicate_contigs_removed = 0
    for contig_candidates in candidates_by_contig.values():
        def candidate_rank(candidate):
            telo_name, raw_idx = candidate
            row = raw_rows[raw_idx]
            is_terminal = is_terminal_telomere_candidate(
                row,
                telo_name,
                telo_bounds,
                chr_len,
            )
            return (
                int(is_terminal),
                -int(row[CTG_MAPQ]),
                -(int(row[CHR_END]) - int(row[CHR_STR])),
                raw_idx,
                telo_name,
            )

        selected_candidates.add(min(contig_candidates, key=candidate_rank))
        duplicate_contigs_removed += len(contig_candidates) - 1

    return [
        candidate
        for candidate in compressed_candidates
        if candidate in selected_candidates
    ], duplicate_contigs_removed


def build_paf_telomere_nodes(paf_path, telo_data, repeat_censat_data,
                             chr_len, base_idx):
    """Build telomere nodes from PAF alignments without consulting VCF calls."""
    raw_rows = import_data(paf_path)
    if not raw_rows:
        return [], [], [], Counter()
    excluded_contigs, excluded_rows = find_multi_end_aligned_contigs(raw_rows)

    telo_dict = telo_data_to_dict(telo_data)
    raw_telo_labels = label_node(raw_rows, telo_dict)
    _, preprocess_report, telo_connect_info = preprocess_telo(
        raw_rows,
        raw_telo_labels,
    )
    # preprocess_telo appends a sentinel in-place while scanning contig groups.
    raw_rows.pop()

    boundary_candidate_count = len(telo_connect_info)
    excluded_boundary_candidate_count = sum(
        raw_idx in excluded_rows for raw_idx in telo_connect_info
    )
    telo_connect_info = {
        raw_idx: telo_name
        for raw_idx, telo_name in telo_connect_info.items()
        if raw_idx not in excluded_rows
    }
    if excluded_contigs:
        logging.info(
            f"Excluded {len(excluded_contigs)} multi-end-aligned contigs "
            f"({len(excluded_rows)} PAF rows, "
            f"{excluded_boundary_candidate_count} telomere boundary candidates) "
            f"from VCF-mode telomere PAF {paf_path}"
        )
    selected_candidates, duplicate_contigs_removed = \
        compress_paf_telomere_candidates(
            raw_rows,
            telo_connect_info,
            telo_dict,
            chr_len,
        )
    if not selected_candidates:
        return [], [], [], Counter({
            "telomere_paf_rows": len(raw_rows),
            "telomere_paf_contigs": 0,
            "telomere_paf_nodes": 0,
            "telomere_paf_edges": 0,
            "telomere_paf_boundary_candidates": boundary_candidate_count,
            "telomere_paf_duplicate_contigs_removed": duplicate_contigs_removed,
            "telomere_paf_multi_end_contigs_excluded": len(excluded_contigs),
            "telomere_paf_multi_end_rows_excluded": len(excluded_rows),
            "telomere_paf_multi_end_boundary_candidates_excluded": excluded_boundary_candidate_count,
        })

    selected_rows = [raw_rows[raw_idx] for _, raw_idx in selected_candidates]
    selected_censat_labels = label_repeat_node(
        selected_rows,
        repeat_censat_data,
        chr_len,
    )
    nodes = []
    edges = []
    report_rows = []
    for local_idx, (telo_name, raw_idx) in enumerate(selected_candidates):
        row = raw_rows[raw_idx]
        global_node_idx = base_idx + local_idx
        node = [
            row[CTG_NAM],
            int(row[CTG_LEN]),
            int(row[CTG_STR]),
            int(row[CTG_END]),
            row[CTG_DIR],
            row[CHR_NAM],
            int(row[CHR_LEN]),
            int(row[CHR_STR]),
            int(row[CHR_END]),
            int(row[CTG_MAPQ]),
            3,
            global_node_idx,
            global_node_idx,
            raw_telo_labels[raw_idx][0],
            raw_telo_labels[raw_idx][1],
            telo_name,
            "0",
            "0",
            selected_censat_labels[local_idx][1],
            row[CTG_DIR],
            row[CHR_NAM],
            f"1.{raw_idx}",
        ]
        nodes.append(node)
        edges.append((telo_name, (DIR_OUT, global_node_idx, 0)))
        report_rows.append((
            "preprocess_telo_boundary",
            telo_name,
            global_node_idx,
            row[CTG_NAM],
            row[CTG_STR],
            row[CTG_END],
            row[CHR_NAM],
            row[CHR_STR],
            row[CHR_END],
            row[CTG_DIR],
            raw_idx,
        ))

    metrics = Counter({
        "telomere_paf_rows": len(raw_rows),
        "telomere_paf_contigs": len(nodes),
        "telomere_paf_nodes": len(nodes),
        "telomere_paf_edges": len(edges),
        "telomere_paf_boundary_candidates": boundary_candidate_count,
        "telomere_paf_duplicate_contigs_removed": duplicate_contigs_removed,
        "telomere_paf_multi_end_contigs_excluded": len(excluded_contigs),
        "telomere_paf_multi_end_rows_excluded": len(excluded_rows),
        "telomere_paf_multi_end_boundary_candidates_excluded": excluded_boundary_candidate_count,
    })
    metrics.update({
        f"telomere_paf_preprocess_{case.lower()}": len(rows)
        for case, rows in preprocess_report.items()
    })
    return nodes, edges, report_rows, metrics


def annotate_synthetic_nodes(contig_data, telo_data, repeat_censat_data,
                             chr_len, annotate_telomeres=True):
    if not hasattr(telo_data, "items"):
        telo_data = telo_data_to_dict(telo_data)
    if annotate_telomeres:
        telo_labels = label_node(contig_data, telo_data)
    else:
        telo_labels = [("0", "0") for _ in contig_data]
    repeat_labels = no_repeat_labels(contig_data)
    censat_labels = label_repeat_node(contig_data, repeat_censat_data, chr_len)
    for idx, node in enumerate(contig_data):
        node[CTG_TELCHR] = telo_labels[idx][0]
        node[CTG_TELDIR] = telo_labels[idx][1]
        node[CTG_RPTCHR] = repeat_labels[idx][0]
        node[CTG_RPTCASE] = repeat_labels[idx][1]
        node[CTG_CENSAT] = censat_labels[idx][1]
        if (
            annotate_telomeres
            and node[CTG_TELCON] == "0"
            and str(node[CTG_NAM]).startswith(VIRTUAL_TELOMERE_NODE_PREFIX)
        ):
            node[CTG_TELCON] = str(node[CTG_NAM])[len(VIRTUAL_TELOMERE_NODE_PREFIX):]


def has_subtelomeric_telomere_node(contig_data, telo_edges, telo_name,
                                   chrom, chrom_len, side):
    for edge_telo_name, edge in telo_edges:
        if edge_telo_name != telo_name:
            continue
        node_idx = int(edge[1])
        if node_idx < 0 or node_idx >= len(contig_data):
            continue
        node = contig_data[node_idx]
        if node[CHR_NAM] != chrom:
            continue
        if side == "f":
            if int(node[CHR_STR]) <= SUBTELOMERE_LENGTH:
                return True
        elif int(chrom_len) - int(node[CHR_END]) <= SUBTELOMERE_LENGTH:
            return True
    return False


def add_missing_virtual_telomeres(contig_data, telo_edges, chr_len, telo_data,
                                  repeat_censat_data):
    virtual_nodes = []

    for chrom in [f"chr{i}" for i in range(1, 23)] + ["chrX"]:
        if chrom not in chr_len:
            continue
        chrom_len = int(chr_len[chrom])
        for side, path_dir, st, nd in (
            ("f", "+", 0, min(VIRTUAL_TELOMERE_FLANK, chrom_len)),
            ("b", "-", max(0, chrom_len - VIRTUAL_TELOMERE_FLANK), chrom_len),
        ):
            telo_name = f"{chrom}{side}"
            if has_subtelomeric_telomere_node(
                contig_data,
                telo_edges,
                telo_name,
                chrom,
                chrom_len,
                side,
            ):
                continue

            node_idx = len(contig_data) + len(virtual_nodes)
            node = make_synthetic_span_node(
                f"{VIRTUAL_TELOMERE_NODE_PREFIX}{telo_name}",
                chrom,
                st,
                nd,
                path_dir,
                chr_len,
                3,
                node_idx,
                node_idx,
            )
            node[CTG_TELCON] = telo_name
            virtual_nodes.append(node)
            telo_edges.append((telo_name, (DIR_OUT, node_idx, 0)))

    if virtual_nodes:
        annotate_synthetic_nodes(
            virtual_nodes,
            telo_data,
            repeat_censat_data,
            chr_len,
        )
        contig_data.extend(virtual_nodes)
    return len(virtual_nodes)


def add_vcf_nclose_pair(contig_data, nclose_nodes, event_name, chr_a, pos_a, dir_a,
                        chr_b, pos_b, dir_b, chr_len, global_idx_start):
    ctg_typ = 1 if chr_a != chr_b else 2
    s_idx = len(contig_data)
    e_idx = s_idx + 1
    node_a = make_synthetic_vcf_node(
        event_name, chr_a, pos_a, dir_a, "exit", chr_len, ctg_typ, s_idx, global_idx_start
    )
    node_b = make_synthetic_vcf_node(
        event_name, chr_b, pos_b, dir_b, "entry", chr_len, ctg_typ, e_idx, global_idx_start + 1
    )
    node_a[CTG_STRND] = node_b[CTG_STRND] = s_idx
    node_a[CTG_ENDND] = node_b[CTG_ENDND] = e_idx
    contig_data.extend([node_a, node_b])
    nclose_nodes[event_name].append((s_idx, e_idx))
    return global_idx_start + 2


def canonical_nclose_snapshot(nclose_nodes):
    return {
        key: tuple(tuple(int(node_idx) for node_idx in pair) for pair in pair_list)
        for key, pair_list in nclose_nodes.items()
    }


def assert_vcf_nclose_has_no_indel_like(contig_data, nclose_nodes, context):
    """Fail if canonical VCF nclose contains a same-chrom/same-dir pair."""

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


def add_vcf_type4_graph_pair(contig_data, event, event_index, chr_len,
                             global_idx_start):
    """Add graph-only breakpoint flanks for one non-INS VCF type4 event."""

    chrom = str(event["chrom"])
    st, nd = sorted((int(event["st"]), int(event["nd"])))
    event_type = str(event["event_type"])
    if event_type == "front_jump":
        exit_pos, entry_pos = st, nd
    elif event_type == "back_jump":
        exit_pos, entry_pos = nd, st
    else:
        raise ValueError(f"Invalid VCF type4 graph event_type: {event_type}")

    event_id = str(event.get("event_id", event.get("vcf_id", event_index)))
    contig_name = f"{VCF_TYPE4_GRAPH_NODE_PREFIX}{event_index}_{event_id}"
    s_idx = len(contig_data)
    e_idx = s_idx + 1
    node_a = make_synthetic_vcf_node(
        contig_name, chrom, exit_pos, "+", "exit", chr_len, 4,
        s_idx, global_idx_start,
    )
    node_b = make_synthetic_vcf_node(
        contig_name, chrom, entry_pos, "+", "entry", chr_len, 4,
        e_idx, global_idx_start + 1,
    )
    node_a[CTG_STRND] = node_b[CTG_STRND] = s_idx
    node_a[CTG_ENDND] = node_b[CTG_ENDND] = e_idx

    len_a = int(node_a[CTG_END]) - int(node_a[CTG_STR])
    len_b = int(node_b[CTG_END]) - int(node_b[CTG_STR])
    total_len = len_a + len_b
    node_a[CTG_LEN] = node_b[CTG_LEN] = total_len
    node_a[CTG_STR], node_a[CTG_END] = 0, len_a
    node_b[CTG_STR], node_b[CTG_END] = len_a, total_len
    contig_data.extend([node_a, node_b])
    return global_idx_start + 2


def build_vcf_mode_inputs():
    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
    vcf_ins_alt_alignments = load_vcf_ins_alt_alignment_spans(args.alt)
    parsed_vcf = parse_vcf_events(
        args.vcf_input,
        chr_len,
        pass_filters=args.vcf_filter_pass,
        ins_alt_alignments=vcf_ins_alt_alignments,
    )

    depth_df = pd.read_csv(
        main_stat_loc,
        compression="gzip",
        comment="#",
        sep="\t",
        names=["chr", "st", "nd", "length", "covsite", "totaldepth", "cov", "meandepth"],
    ).query('chr != "chrM"')

    repeat_censat_data = import_censat_repeat_data(CENSAT_PATH)
    telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)

    summary = dict(parsed_vcf.summary)
    skipped_records = list(parsed_vcf.skipped_records)
    orientation_mismatches = list(parsed_vcf.orientation_mismatches)
    contig_data = []
    nclose_nodes = defaultdict(list)
    # VCF indels are never added to the step-02 nclose graph.  The parser
    # places only >=100 kb indels here; step 11 consumes this handoff file.
    step11_vcf_type4_events = list(parsed_vcf.type4_events)
    global_idx = 0

    for spec in parsed_vcf.nclose_specs:
        global_idx = add_vcf_nclose_pair(
            contig_data,
            nclose_nodes,
            spec.event_name,
            spec.chrom_a,
            spec.pos_a,
            spec.dir_a,
            spec.chrom_b,
            spec.pos_b,
            spec.dir_b,
            chr_len,
            global_idx,
        )

    assert_vcf_nclose_has_no_indel_like(
        contig_data,
        nclose_nodes,
        "VCF parser handoff",
    )

    vcf_type4_graph_events = (
        select_vcf_type4_graph_events(step11_vcf_type4_events)
        if args.add_indel_graph
        else []
    )
    for event_index, event in enumerate(vcf_type4_graph_events):
        global_idx = add_vcf_type4_graph_pair(
            contig_data,
            event,
            event_index,
            chr_len,
            global_idx,
        )
    if args.add_indel_graph:
        invalid_graph_svtypes = sorted({
            str(event.get("svtype", "")).upper()
            for event in vcf_type4_graph_events
            if str(event.get("svtype", "")).upper() not in {"DEL", "DUP", "BND"}
        })
        if invalid_graph_svtypes:
            raise AssertionError(
                "VCF --add_indel_graph received ineligible SVTYPE(s): "
                + ", ".join(invalid_graph_svtypes)
            )
        graph_excluded_ins = sum(
            str(event.get("svtype", "")).upper() == "INS"
            for event in step11_vcf_type4_events
        )
        logging.info(
            f"VCF --add_indel_graph: prepared {len(vcf_type4_graph_events)} "
            f"DEL/DUP/BND Indel-like event(s); excluded {graph_excluded_ins} INS"
        )

    # Everything accumulated so far was created from VCF records.  Keep its
    # telomere labels at zero; only PAF/virtual nodes may carry telomere data.
    vcf_node_count = len(contig_data)
    telo_edges = []
    telomere_paf_report_rows = []
    telomere_paf_metrics = Counter()
    if args.paf_file_path is not None:
        # Telomere nodes come exclusively from the assembly PAF.  VCF records
        # above contribute nclose/type4 events but never telomere candidates.
        telomere_paf_nodes, telomere_paf_edges, telomere_paf_report_rows, telomere_paf_metrics = \
            build_paf_telomere_nodes(
                args.paf_file_path,
                telo_data,
                repeat_censat_data,
                chr_len,
                len(contig_data),
            )
        contig_data.extend(telomere_paf_nodes)
        telo_edges.extend(telomere_paf_edges)

    virtual_telomere_nodes = add_missing_virtual_telomeres(
        contig_data,
        telo_edges,
        chr_len,
        telo_data,
        repeat_censat_data,
    )
    if virtual_telomere_nodes:
        logging.info(
            f"Added {virtual_telomere_nodes} missing chromosome-end virtual telomere nodes"
        )

    annotate_synthetic_nodes(
        contig_data[:vcf_node_count],
        telo_data,
        repeat_censat_data,
        chr_len,
        annotate_telomeres=False,
    )
    contig_data = [tuple(row) for row in contig_data]

    with open(PREPROCESSED_PAF_FILE_PATH, "wt") as f:
        for row in contig_data:
            print("\t".join(map(str, row)), file=f)
    with open(PAF_FILE_PATH[0], "wt") as f:
        for row in contig_data:
            print("\t".join(map(str, synthetic_node_to_paf_row(row))), file=f)

    write_telomere_connected_outputs(PREFIX, telo_edges, contig_data)
    with open(f"{PREFIX}/{VCF_TELOMERE_PAF_NODES_TSV}", "wt") as f:
        print("kind\ttelomere\tnode_idx\tcontig\tquery_start\tquery_end\tchrom\tref_start\tref_end\tdir\tpaf_row_idx", file=f)
        for row in telomere_paf_report_rows:
            print("\t".join(map(str, row)), file=f)

    raw_nclose_nodes = nclose_nodes
    all_nclose_comp = defaultdict(list)
    for key, pair_list in nclose_nodes.items():
        for pair in pair_list:
            all_nclose_comp[key].append(tuple(pair))

    with open(f"{PREFIX}/conjoined_type4_ins_del.pkl", "wb") as f:
        pkl.dump(([], []), f)
    with open(f"{PREFIX}/{VCF_TYPE4_EVENTS_PKL}", "wb") as f:
        pkl.dump(step11_vcf_type4_events, f)
    with open(f"{PREFIX}/indel_exclude_idx_set.pkl", "wb") as f:
        pkl.dump(set(), f)

    _, cen_fragment_meta = find_breakend_centromere(
        repeat_censat_data,
        chr_len,
        depth_df,
        raw_nclose_nodes=raw_nclose_nodes,
        contig_data=contig_data,
        log_context="VCF mode",
    )
    with open(f"{PREFIX}/cen_fragment_data.pkl", "wb") as f:
        pkl.dump(cen_fragment_meta, f)

    chr_corr, chr_rev_corr = chr_correlation_maker(contig_data)
    telo_contig = extract_telomere_connect_contig(f"{PREFIX}/telomere_connected_list.txt")
    telo_set = {edge[1] for edge_list in telo_contig.values() for edge in edge_list}
    telo_node_count = len(telo_set)

    write_virtual_ordinary_contig(
        f"{PREFIX}/virtual_ordinary_contig.txt",
        [],
    )
    all_nclose_type = group_nclose_nodes_by_chrom(contig_data, all_nclose_comp)
    uncomp_node_count = write_nclose_nodes_list(
        f"{PREFIX}/all_nclose_nodes_list.txt",
        all_nclose_type,
        contig_data,
    )
    compressed_nclose_type = group_nclose_nodes_by_chrom(contig_data, nclose_nodes)
    write_nclose_nodes_list(
        f"{PREFIX}/compressed_nclose_nodes_list.txt",
        compressed_nclose_type,
        contig_data,
    )
    nclose_node_count = write_nclose_nodes_index(
        f"{PREFIX}/nclose_nodes_index.txt",
        nclose_nodes,
        contig_data,
    )

    summary.update({
        "nclose_pairs": sum(len(v) for v in nclose_nodes.values()),
        "synthetic_nodes": len(contig_data),
        "virtual_telomere_nodes": virtual_telomere_nodes,
        "type4_min_span": VCF_TYPE4_MIN_SPAN,
        "vcf_ins_alt_alignment_queries": len(vcf_ins_alt_alignments),
        "vcf_type4_graph_events": len(vcf_type4_graph_events),
        "vcf_type4_graph_nodes": 2 * len(vcf_type4_graph_events),
        "vcf_type4_graph_ins_excluded": (
            sum(
                str(event.get("svtype", "")).upper() == "INS"
                for event in step11_vcf_type4_events
            )
            if args.add_indel_graph
            else 0
        ),
    })
    summary.update(telomere_paf_metrics)
    with open(f"{PREFIX}/{VCF_MODE_SUMMARY_JSON}", "wt") as f:
        json.dump(dict(summary), f, indent=2, sort_keys=True)
    with open(f"{PREFIX}/{VCF_MODE_SUMMARY_TSV}", "wt") as f:
        print("metric\tvalue", file=f)
        for key, value in sorted(summary.items()):
            print(f"{key}\t{value}", file=f)
    with open(f"{PREFIX}/{VCF_SKIPPED_RECORDS_TSV}", "wt") as f:
        print("line_no\tid\tsvtype\treason\traw", file=f)
        for row in skipped_records:
            print("\t".join(map(str, row)), file=f)
    with open(f"{PREFIX}/{VCF_ORIENTATION_MISMATCH_TSV}", "wt") as f:
        print("line_no\tid\tmate_line_no\tmate_id\tdirs\tmate_dirs\talt\tmate_alt", file=f)
        for row in orientation_mismatches:
            print("\t".join(map(str, row)), file=f)

    no_chrY = False
    if "chrY" in set(depth_df["chr"]):
        ydf = depth_df.query('chr == "chrY"')
        if len(ydf) > 0:
            no_chrY = (len(ydf.query("meandepth != 0")) / len(ydf)) < chrY_MINIMUM_RATIO

    if VCF_TYPE4_MIN_SPAN % M == 0:
        vcf_indel_min_size_label = f"{VCF_TYPE4_MIN_SPAN // M}Mbp"
    elif VCF_TYPE4_MIN_SPAN % K == 0:
        vcf_indel_min_size_label = f"{VCF_TYPE4_MIN_SPAN // K}Kbp"
    else:
        vcf_indel_min_size_label = f"{VCF_TYPE4_MIN_SPAN}bp"

    ins_skip_reasons = []
    if summary["skipped_ins_small_size"] > 0:
        reason = f"SIZE <{vcf_indel_min_size_label}"
        if summary["skipped_ins_small_size"] != summary["skipped_ins"]:
            reason = f"{reason}: {summary['skipped_ins_small_size']}"
        ins_skip_reasons.append(reason)
    if summary["skipped_ins_no_alt_alignment"] > 0:
        reason = "no ALT alignment"
        if summary["skipped_ins_no_alt_alignment"] != summary["skipped_ins"]:
            reason = f"{reason}: {summary['skipped_ins_no_alt_alignment']}"
        ins_skip_reasons.append(reason)
    ins_skip_detail = f" ({', '.join(ins_skip_reasons)})" if ins_skip_reasons else ""

    logging.info(
        f"VCF mode: {summary['used_bnd_events']} BND, {summary['used_inv_events']} INV, "
        f"{summary['used_type4_events']} Indel events "
        f"({summary['used_bnd_type4_events']} BND Indel), "
        f"{summary['skipped_ins']} INS skipped{ins_skip_detail}"
    )

    return {
        "df": depth_df,
        "no_chrY": no_chrY,
        "repeat_censat_data": repeat_censat_data,
        "chr_len": chr_len,
        "contig_data": contig_data,
        "contig_data_size": len(contig_data),
        "chr_corr": chr_corr,
        "chr_rev_corr": chr_rev_corr,
        "telo_contig": telo_contig,
        "telo_node_count": telo_node_count,
        "telo_set": telo_set,
        "rpt_con": set(),
        "rpt_censat_con": set(),
        "bnd_contig": {key for key in nclose_nodes},
        "raw_nclose_nodes": raw_nclose_nodes,
        "nclose_nodes": nclose_nodes,
        "nclose_start_compress": defaultdict(dict),
        "nclose_end_compress": defaultdict(dict),
        "vctg_dict": {},
        "all_nclose_comp": all_nclose_comp,
        "nclose_coverage": Counter(),
        "nclose_compress_track": defaultdict(list),
        "st_compress": {},
        "ed_compress": {},
        "uncomp_node_count": uncomp_node_count,
        "nclose_node_count": nclose_node_count,
        "transloc_nclose_pair_count": sum(
            1
            for pair_list in nclose_nodes.values()
            for a, b in pair_list
            if contig_data[a][CHR_NAM] != contig_data[b][CHR_NAM]
        ),
        "indel_exclude_idx_set": set(),
        "telo_coverage": Counter(),
    }


def nclose_calc():
    repeat_censat_data = import_censat_repeat_data(CENSAT_PATH)
    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
    contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)

    TELO_CONNECT_NODES_INFO_PATH = PREFIX+"/telomere_connected_list.txt"

    os.makedirs(PREFIX, exist_ok=True)

    contig_data_size = len(contig_data)

    chr_corr, chr_rev_corr = chr_correlation_maker(contig_data)

    telo_contig = extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH)

    telo_node_count = 0
    telo_set = set()

    for i in telo_contig:
        for j in telo_contig[i]:
            telo_node_count+=1
            telo_set.add(j[1])

    repeat_data = import_repeat_data_00(REPEAT_INFO_FILE_PATH)
    repeat_censat_data = import_repeat_data(CENSAT_PATH)
    rpt_con = extract_all_repeat_contig(contig_data, repeat_data, CTG_RPTCASE, NON_REPEAT_NOISE_RATIO)
    rpt_censat_con = extract_all_repeat_contig(contig_data, repeat_censat_data, CTG_CENSAT)
    rpt_censat_con = check_censat_contig(rpt_censat_con, PAF_FILE_PATH, ORIGINAL_PAF_LOC_LIST, contig_data)

    bnd_contig = extract_bnd_contig(contig_data)


    # Type 1, 2, 4에 대해서 
    raw_nclose_nodes, nclose_start_compress, nclose_end_compress, vctg_dict, all_nclose_comp, nclose_coverage, nclose_compress_track = \
        extract_nclose_node(contig_data, bnd_contig, rpt_con, rpt_censat_con, repeat_censat_data, PAF_FILE_PATH, ORIGINAL_PAF_LOC_LIST, telo_set, telo_contig, chr_len, asm2cov)

    depth_df = df if 'df' in globals() else pd.read_csv(
        main_stat_loc,
        compression='gzip',
        comment='#',
        sep='\t',
        names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'],
    ).query('chr != "chrM"')
    _, cen_fragment_meta = find_breakend_centromere(
        repeat_censat_data,
        chr_len,
        depth_df,
        raw_nclose_nodes=raw_nclose_nodes,
        contig_data=contig_data,
        log_context="Raw-nclose adjusted",
    )
    with open(f'{PREFIX}/cen_fragment_data.pkl', 'wb') as f:
        pkl.dump(cen_fragment_meta, f)

    not_using_nclose_node = set()
    type1_nclose_node = []
    type2_nclose_node = defaultdict(list)

    nclose_nodes = raw_nclose_nodes
    for j in nclose_nodes:
        for s, e in nclose_nodes[j]:
            curr_contig_first_fragment = contig_data[s]
            curr_contig_end_fragment = contig_data[e]

            if curr_contig_first_fragment[CHR_NAM] == curr_contig_end_fragment[CHR_NAM]:
                type2_nclose_node[curr_contig_first_fragment[CHR_NAM]].append((s, e))
                if curr_contig_first_fragment[CTG_DIR] == '+':
                    inside_st, inside_nd = sorted([curr_contig_end_fragment[CHR_END], curr_contig_first_fragment[CHR_STR]])
                else:
                    inside_st, inside_nd = sorted([curr_contig_first_fragment[CHR_END], curr_contig_end_fragment[CHR_STR]])
                chukji = inside_nd - inside_st
                chukji_chrom = curr_contig_first_fragment[CHR_NAM]

                if (not exist_near_bnd_point(chukji_chrom, inside_st)) and (not exist_near_bnd_point(chukji_chrom, inside_nd)) and \
                   (not censat_overlap_check(repeat_censat_data, chukji_chrom, inside_st, inside_nd)):
                    not_using_nclose_node.add((s, e))
            else:
                type1_nclose_node.append((s, e))

    _, centromere_slave = similar_centromere_nclose_cluster(nclose_nodes, contig_data, repeat_censat_data, chr_len)

    with open(f"{PREFIX}/conjoined_type4_ins_del.pkl", "wb") as f:
        pkl.dump(conjoined_type4(contig_data, type2_nclose_node), f)

    saved_not_using_nclose_node = set()
    conjoined_nclose_node_set = set()

    for nunclose in not_using_nclose_node:
        chrom = contig_data[nunclose[0]][CHR_NAM]
        for t1nn in type1_nclose_node:
            if contig_data[t1nn[0]][CHR_NAM] != chrom and contig_data[t1nn[1]][CHR_NAM] != chrom:
                continue
            type2_contig_front = contig_data[nunclose[0]]
            type2_contig_back = contig_data[nunclose[1]]
            if contig_data[t1nn[0]][CHR_NAM] == chrom:
                template_contig = contig_data[t1nn[0]]
                if template_contig[CTG_DIR] == '+':
                    dist = distance_checker(type2_contig_back, template_contig)
                    if (type2_contig_back[CHR_END] < template_contig[CHR_STR] and dist < TYPE2_CONTIG_MINIMUM_LENGTH) or dist == 0:
                        saved_not_using_nclose_node.add(nunclose)
                        conjoined_nclose_node_set.add((t1nn[0], nunclose[1]))
                else:
                    dist = distance_checker(type2_contig_front, template_contig)
                    if (type2_contig_front[CHR_END] > template_contig[CHR_STR] and dist < TYPE2_CONTIG_MINIMUM_LENGTH) or dist == 0:
                        saved_not_using_nclose_node.add(nunclose)
                        conjoined_nclose_node_set.add((t1nn[0], nunclose[0]))
            elif contig_data[t1nn[1]][CHR_NAM] == chrom:
                template_contig = contig_data[t1nn[1]]  
                if template_contig[CTG_DIR] == '+':
                    dist = distance_checker(type2_contig_front, template_contig)
                    if (type2_contig_front[CHR_STR] > template_contig[CHR_END] and dist < TYPE2_CONTIG_MINIMUM_LENGTH) or dist == 0:
                        saved_not_using_nclose_node.add(nunclose)
                        conjoined_nclose_node_set.add((t1nn[1], nunclose[0]))
                else:
                    dist = distance_checker(type2_contig_back, template_contig)
                    if (type2_contig_back[CHR_STR] < template_contig[CHR_END] and dist < TYPE2_CONTIG_MINIMUM_LENGTH) or dist == 0:
                        saved_not_using_nclose_node.add(nunclose)
                        conjoined_nclose_node_set.add((t1nn[1], nunclose[1]))

    not_using_nclose_node |= set(centromere_slave.keys())

    for conjoined_nclose in conjoined_nclose_node_set:
        pass
        # do something

    virtual_ordinary_contig = make_virtual_ord_ctg(contig_data, vctg_dict)
    write_virtual_ordinary_contig(
        f"{PREFIX}/virtual_ordinary_contig.txt",
        virtual_ordinary_contig,
    )

    all_nclose_type = group_nclose_nodes_by_chrom(contig_data, all_nclose_comp)
    uncomp_node_count = write_nclose_nodes_list(
        f"{PREFIX}/all_nclose_nodes_list.txt",
        all_nclose_type,
        contig_data,
        rpt_con,
    )

    st_compress = dict()
    for ddict in nclose_start_compress.values():
        for ctg1, ctg2_list in ddict.items():
            for pair in nclose_nodes[ctg1]:
                ctg1_idx = 0
                contig_a = contig_data[pair[0]]
                contig_b = contig_data[pair[1]]
                if chr2int(contig_a[CHR_NAM]) <= chr2int(contig_b[CHR_NAM]):
                    st_compress[pair[0]] = pair[0]
                    ctg1_idx = pair[0]
                else:
                    st_compress[pair[1]] = pair[1]
                    ctg1_idx = pair[1]
            for ctg2 in ctg2_list:
                for pair in nclose_nodes[ctg2]:
                    contig_a = contig_data[pair[0]]
                    contig_b = contig_data[pair[1]]
                    if chr2int(contig_a[CHR_NAM]) <= chr2int(contig_b[CHR_NAM]):
                        st_compress[pair[0]] = ctg1_idx
                    else:
                        st_compress[pair[1]] = ctg1_idx

    ed_compress = dict()
    for ddict in nclose_end_compress.values():
        for ctg1, ctg2_list in ddict.items():
            for pair in nclose_nodes[ctg1]:
                ctg1_idx = 0
                contig_a = contig_data[pair[0]]
                contig_b = contig_data[pair[1]]
                if chr2int(contig_a[CHR_NAM]) <= chr2int(contig_b[CHR_NAM]):
                    ed_compress[pair[1]] = pair[1]
                    ctg1_idx = pair[1]
                else:
                    ed_compress[pair[0]] = pair[0]
                    ctg1_idx = pair[0]
            for ctg2 in ctg2_list:
                for pair in nclose_nodes[ctg2]:
                    contig_a = contig_data[pair[1]]
                    contig_b = contig_data[pair[0]]
                    if chr2int(contig_a[CHR_NAM]) <= chr2int(contig_b[CHR_NAM]):
                        ed_compress[pair[1]] = ctg1_idx
                    else:
                        ed_compress[pair[0]] = ctg1_idx
      
    final_nclose_nodes = defaultdict(list)
    for i in nclose_nodes:
        for j in nclose_nodes[i]:
            if j not in not_using_nclose_node or j in saved_not_using_nclose_node:
                final_nclose_nodes[i].append(j)
    
    nclose_nodes = final_nclose_nodes

    # censat 관련 nclose 필터:
    # - 양쪽 모두 censat: both endpoint MAPQ == 60 이면 keep, MAPQ < 60 이 있으면 drop
    #   단, 같은 염색체 censat +/- pair 는 drop
    #   단, 같은 두 censat locus 근처에서 보정 방향이 서로 반대인 pair 세트는 drop
    # - cen_fragment column 이 만들어진 chr 의 censat endpoint 는 아래에서
    #   depth 방향과 nclose 방향이 모두 match 해야만 보존한다.
    with open(f'{PREFIX}/cen_fragment_data.pkl', 'rb') as f:
        cen_fragment_meta_for_filter = pkl.load(f)
    cent_fragment_chroms = set(cen_fragment_meta_for_filter.keys())

    # chunk가 overlap하는 censat interval이 chr 경계와 정확히 맞닿은 경우만 chr 끝 판정
    # CHM13 v2.1 기준: chr13/14/15/21/22 (acrocentric) 의 censat이 s=0 에서 시작 -> match
    def _censat_at_chr_end(idx):
        chrom = contig_data[idx][CHR_NAM]
        cl = chr_len.get(chrom, 0)
        intervals = repeat_censat_data.get(chrom, [])
        chunk_str = contig_data[idx][CHR_STR]
        chunk_end = contig_data[idx][CHR_END]
        for s, e in intervals:
            if max(s, chunk_str) < min(e, chunk_end):
                if s == 0 or e == cl:
                    return True
        return False

    def _canonical_censat_pair_info(pair):
        a, b = pair
        contig_a = contig_data[a]
        contig_b = contig_data[b]
        a_key = (chr2int(contig_a[CHR_NAM]), contig_a[CHR_STR], contig_a[CHR_END])
        b_key = (chr2int(contig_b[CHR_NAM]), contig_b[CHR_STR], contig_b[CHR_END])
        if b_key < a_key:
            a, b = b, a
            contig_a, contig_b = contig_b, contig_a

        is_for = a < b
        return {
            "chroms": (contig_a[CHR_NAM], contig_b[CHR_NAM]),
            "idxs": (a, b),
            "dirs": (
                get_corr_dir(is_for, contig_a[CTG_DIR]),
                get_corr_dir(is_for, contig_b[CTG_DIR]),
            ),
        }

    def _opposite_dir_pair(dirs_a, dirs_b):
        return dirs_a[0] != dirs_b[0] and dirs_a[1] != dirs_b[1]

    def _flip_dir(ctg_dir):
        return '-' if ctg_dir == '+' else '+'

    def _cen_fragment_target_dir(chrom):
        return '-' if cen_fragment_meta_for_filter[chrom]['dir'] else '+'

    def _normalized_censat_endpoint_dirs(pair):
        endpoint_dirs = []
        for order_idx, node_idx in enumerate(pair):
            if contig_data[node_idx][CTG_CENSAT] == '0':
                continue
            ctg_dir = contig_data[node_idx][CTG_DIR]
            norm_dir = ctg_dir if order_idx == 0 else _flip_dir(ctg_dir)
            endpoint_dirs.append((node_idx, norm_dir))
        return endpoint_dirs

    censat_censat_mapq60_groups = defaultdict(list)
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            s_censat = contig_data[pair[0]][CTG_CENSAT] != '0'
            e_censat = contig_data[pair[1]][CTG_CENSAT] != '0'
            if not (s_censat and e_censat):
                continue
            if contig_data[pair[0]][CTG_MAPQ] < 60 or contig_data[pair[1]][CTG_MAPQ] < 60:
                continue
            pair_info = _canonical_censat_pair_info(pair)
            censat_censat_mapq60_groups[pair_info["chroms"]].append((ctg_name, tuple(pair), pair_info))

    opposite_censat_censat_pairs = set()
    for items in censat_censat_mapq60_groups.values():
        for i in range(len(items)):
            ctg_name_i, pair_i, info_i = items[i]
            a_i, b_i = info_i["idxs"]
            for j in range(i + 1, len(items)):
                ctg_name_j, pair_j, info_j = items[j]
                if not _opposite_dir_pair(info_i["dirs"], info_j["dirs"]):
                    continue
                a_j, b_j = info_j["idxs"]
                if distance_checker(contig_data[a_i], contig_data[a_j]) <= NCLOSE_COMPRESS_LIMIT \
                   and distance_checker(contig_data[b_i], contig_data[b_j]) <= NCLOSE_COMPRESS_LIMIT:
                    opposite_censat_censat_pairs.add((ctg_name_i, pair_i))
                    opposite_censat_censat_pairs.add((ctg_name_j, pair_j))

    censat_filtered_nclose_nodes = defaultdict(list)
    censat_removed_both = 0
    censat_removed_terminal_both = 0
    censat_removed_same_chrom_opposite_dir = 0
    censat_removed_opposite_dir = 0
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            s_censat = contig_data[pair[0]][CTG_CENSAT] != '0'
            e_censat = contig_data[pair[1]][CTG_CENSAT] != '0'
            if s_censat and e_censat:
                if _censat_at_chr_end(pair[0]) or _censat_at_chr_end(pair[1]):
                    censat_removed_terminal_both += 1
                    continue
                if contig_data[pair[0]][CHR_NAM] == contig_data[pair[1]][CHR_NAM] \
                   and contig_data[pair[0]][CTG_DIR] != contig_data[pair[1]][CTG_DIR]:
                    censat_removed_same_chrom_opposite_dir += 1
                    continue
                if (ctg_name, tuple(pair)) in opposite_censat_censat_pairs:
                    censat_removed_opposite_dir += 1
                    continue
                if contig_data[pair[0]][CTG_MAPQ] < 60 or contig_data[pair[1]][CTG_MAPQ] < 60:
                    censat_removed_both += 1
                    continue
            censat_filtered_nclose_nodes[ctg_name].append(pair)

    nclose_nodes = censat_filtered_nclose_nodes
    logging.info(
        f'Removed {censat_removed_both} censat-censat nclose where either endpoint MAPQ < 60, '
        f'{censat_removed_terminal_both} censat-censat nclose with a terminal-censat endpoint, '
        f'{censat_removed_same_chrom_opposite_dir} same-chromosome opposite-direction censat-censat nclose, '
        f'{censat_removed_opposite_dir} opposite-direction censat-censat nclose'
    )

    # censat-censat 의 cen_fragment 방향 mismatch 필터:
    # 두 endpoint 중 cen_fragment 로 잡힌 censat 은 모두 depth 방향과 match 해야만 보존한다.
    #   pair[0] censat: CTG_DIR 그대로
    #   pair[1] censat: pair 를 reverse-complement 해서 censat 이 앞에 온 것으로 보고 CTG_DIR flip
    # 정규화 후 '+' 는 앞쪽(low coord, p-arm) 높음(dir == False), '-' 는 뒤쪽(high coord, q-arm)
    # 높음(dir == True) 과 match 해야 한다.
    cen_fragment_dir_filtered_nclose_nodes = defaultdict(list)
    cen_fragment_dir_removed_censat_censat = 0
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            s_is_censat = contig_data[pair[0]][CTG_CENSAT] != '0'
            e_is_censat = contig_data[pair[1]][CTG_CENSAT] != '0'
            if not (s_is_censat and e_is_censat):
                cen_fragment_dir_filtered_nclose_nodes[ctg_name].append(pair)
                continue

            mismatch = False
            for node_idx, norm_dir in _normalized_censat_endpoint_dirs(pair):
                chrom = contig_data[node_idx][CHR_NAM]
                if chrom not in cent_fragment_chroms:
                    continue
                if norm_dir != _cen_fragment_target_dir(chrom):
                    mismatch = True
                    break

            if mismatch:
                cen_fragment_dir_removed_censat_censat += 1
                continue
            cen_fragment_dir_filtered_nclose_nodes[ctg_name].append(pair)
    nclose_nodes = cen_fragment_dir_filtered_nclose_nodes
    logging.info(
        f'Removed {cen_fragment_dir_removed_censat_censat} censat-censat '
        f'nclose pairs with cen_fragment direction mismatch'
    )

    # If a non-cent-fragment censat endpoint competes for the same non-censat
    # locus, prefer the synthetic simple_ctg_alt candidate in compressed output.
    simple_alt_candidates = defaultdict(list)
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            s_is_censat = contig_data[pair[0]][CTG_CENSAT] != '0'
            e_is_censat = contig_data[pair[1]][CTG_CENSAT] != '0'
            if s_is_censat == e_is_censat:
                continue
            if s_is_censat:
                cidx, nidx = pair[0], pair[1]
            else:
                cidx, nidx = pair[1], pair[0]

            censat_chr = contig_data[cidx][CHR_NAM]
            if censat_chr in cent_fragment_chroms:
                continue

            noncensat_chr = contig_data[nidx][CHR_NAM]
            noncensat_st = contig_data[nidx][CHR_STR]
            noncensat_nd = contig_data[nidx][CHR_END]
            is_simple_alt = (
                ctg_name.startswith('simple_ctg_alt_') or
                contig_data[pair[0]][CTG_NAM].startswith('simple_ctg_alt_') or
                contig_data[pair[1]][CTG_NAM].startswith('simple_ctg_alt_')
            )
            simple_alt_candidates[(censat_chr, noncensat_chr)].append(
                (ctg_name, pair, noncensat_st, noncensat_nd, is_simple_alt)
            )

    simple_alt_to_remove = set()
    for items in simple_alt_candidates.values():
        items.sort(key=lambda x: (x[2], x[3]))
        i = 0
        while i < len(items):
            j = i + 1
            group_end = items[i][3]
            while j < len(items) and items[j][2] < group_end:
                group_end = max(group_end, items[j][3])
                j += 1
            sub = items[i:j]
            if any(it[4] for it in sub):
                for it in sub:
                    if not it[4]:
                        simple_alt_to_remove.add((it[0], tuple(it[1])))
            i = j

    simple_alt_filtered_nclose_nodes = defaultdict(list)
    simple_alt_removed = 0
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            if (ctg_name, tuple(pair)) in simple_alt_to_remove:
                simple_alt_removed += 1
                continue
            simple_alt_filtered_nclose_nodes[ctg_name].append(pair)
    nclose_nodes = simple_alt_filtered_nclose_nodes
    logging.info(
        f'Removed {simple_alt_removed} non-simple censat-noncensat nclose pairs '
        f'overlapping a simple_ctg_alt non-censat locus'
    )

    if ENABLE_CENSAT_NONCENSAT_OFFSET_DIR_FILTER:
        missing_cen_fragment_dir_groups = collect_missing_cen_fragment_dir_censat_noncensat(
            contig_data,
            nclose_nodes,
            cen_fragment_meta_for_filter,
        )
        raw_censat_type2_candidates = extract_raw_censat_type2_candidates(
            args.alt,
            repeat_censat_data,
        )
        combined_added = add_nearest_combined_censat_noncensat_ncloses(
            contig_data,
            nclose_nodes,
            missing_cen_fragment_dir_groups,
            raw_censat_type2_candidates,
            repeat_censat_data,
            chr_len,
            cen_fragment_meta_for_filter,
            PREPROCESSED_PAF_FILE_PATH,
        )
        logging.info(
            f'Added {combined_added} nearest combined censat-noncensat nclose pairs '
            f'from {len(raw_censat_type2_candidates)} raw censat-internal type2 candidates '
            f'across {len(missing_cen_fragment_dir_groups)} target-missing groups'
        )
    else:
        combined_added = 0
        logging.info(
            'Censat-noncensat offset direction mismatch filter disabled; '
            'skip nearest combined censat-noncensat nclose construction'
        )
    if combined_added > 0:
        contig_data_size = len(contig_data)
        chr_corr, chr_rev_corr = chr_correlation_maker(contig_data)

    # Subtelomeric tip-orientation filter:
    # 양끝 SUBTELO_TIP_LIMIT tip을 동등 boundary로 가정. 두 chunk가 같은 telomere
    # orientation type ({f+, b-} = telo-inward 또는 {f-, b+} = telo-outward) 이면
    # 결과 path가 "정상 chromosome path + tip swap"으로 환원되어 정보 추가가 없음.
    # 다른 type끼리는 양쪽 본체가 다 path에 들어가는 fusion/inversion이라 보존.
    # 같은 chromosome head-to-tail (chr1b+ => chr1f+) 은 type이 달라 자동 보존됨.
    #
    # Telomere anchor가 실제 chromosome 끝이 아니라 중간 좌표에 잡힌 경우도 있다.
    # telomere_connected_list의 anchor와 겹치거나 SUBTELO_TIP_LIMIT 안에서 연결 가능한
    # chunk는 물리적 chr tip과 같은 방식으로 orientation을 부여한다.
    telo_anchor_by_chrom = defaultdict(list)
    telo_anchor_seen = set()

    def _add_telo_anchor(telo_name, idx):
        if idx < 0 or idx >= len(contig_data) or telo_name == '0':
            return
        key = (telo_name, idx)
        if key in telo_anchor_seen:
            return
        telo_anchor_seen.add(key)
        telo_anchor_by_chrom[contig_data[idx][CHR_NAM]].append((telo_name, idx))

    for telo_name, edge_list in telo_contig.items():
        for edge in edge_list:
            _add_telo_anchor(telo_name, edge[1])

    for idx, chunk in enumerate(contig_data):
        if chunk[CTG_TELCON] != '0':
            _add_telo_anchor(chunk[CTG_TELCON], idx)

    def _orientation_from_telo_side(telo_side, ctg_dir):
        if telo_side == 'f':
            return 'inward' if ctg_dir == '+' else 'outward'
        return 'outward' if ctg_dir == '+' else 'inward'

    def _subtelo_orientations(idx):
        chrom = contig_data[idx][CHR_NAM]
        chr_str = contig_data[idx][CHR_STR]
        chr_end = contig_data[idx][CHR_END]
        ctg_dir = contig_data[idx][CTG_DIR]
        cl = chr_len.get(chrom, 0)
        at_front = chr_str < SUBTELO_TIP_LIMIT
        at_back = chr_end > cl - SUBTELO_TIP_LIMIT
        orientations = set()
        if at_front:
            orientations.add(_orientation_from_telo_side('f', ctg_dir))
        if at_back:
            orientations.add(_orientation_from_telo_side('b', ctg_dir))

        for telo_name, anchor_idx in telo_anchor_by_chrom.get(chrom, []):
            if distance_checker(contig_data[idx], contig_data[anchor_idx]) < SUBTELO_TIP_LIMIT:
                orientations.add(_orientation_from_telo_side(telo_name[-1], ctg_dir))

        return orientations

    subtelo_filtered_nclose_nodes = defaultdict(list)
    subtelo_removed = 0
    for ctg_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            t_s = _subtelo_orientations(pair[0])
            t_e = _subtelo_orientations(pair[1])
            if t_s and not t_s.isdisjoint(t_e):
                subtelo_removed += 1
                continue
            subtelo_filtered_nclose_nodes[ctg_name].append(pair)
    nclose_nodes = subtelo_filtered_nclose_nodes
    logging.info(f"Removed {subtelo_removed} same-orientation subtelomeric tip nclose "
                 f"(both junctions within {SUBTELO_TIP_LIMIT//1000}kb of chromosome end or telomere anchor, "
                 f"same telo-in/out direction)")

    # Final offset direction mismatch filter:
    # Run once after synthetic combined nclose construction so the added
    # censat-noncensat pairs participate in the same direction arbitration.
    if ENABLE_CENSAT_NONCENSAT_OFFSET_DIR_FILTER:
        nclose_nodes, offset_removed = filter_offset_direction_mismatched_censat_noncensat(
            contig_data,
            nclose_nodes,
            cen_fragment_meta_for_filter,
        )
        logging.info(f'Removed {offset_removed} offset-direction-mismatched censat-noncensat nclose pairs')
    else:
        logging.info('Skipped offset-direction-mismatched censat-noncensat nclose filtering')

    def write_compressed_nclose_nodes_list(current_nclose_nodes):
        nclose_type = group_nclose_nodes_by_chrom(contig_data, current_nclose_nodes)

        transloc_nclose_pair_count = 0
        for chrom_pair, pair_list in nclose_type.items():
            if chrom_pair[0] != chrom_pair[1]:
                transloc_nclose_pair_count += len(pair_list)

            st_flag = (
                (('=', chrom_pair[0]), ('=', chrom_pair[1]))
                in nclose_start_compress
            )
            for pair in pair_list:
                contig_a = contig_data[pair[0]]
                if st_flag:
                    if contig_a[CTG_NAM] in nclose_start_compress[
                        (('=', chrom_pair[0]), ('=', chrom_pair[1]))
                    ]:
                        pass

                original_nclose = tuple(sorted(pair))
                assert(
                    original_nclose not in not_using_nclose_node
                    or original_nclose in saved_not_using_nclose_node
                )

        write_nclose_nodes_list(
            f"{PREFIX}/compressed_nclose_nodes_list.txt",
            nclose_type,
            contig_data,
            rpt_con,
        )
        return transloc_nclose_pair_count

    if args.check_nclose_count:
        if SKIP_BAM_ANAL:
            logging.warning("NClose raw-count VAF filter requested but BAM analysis is skipped")
        else:
            nclose_count_candidates = build_nclose_count_candidates(contig_data, nclose_nodes)
            with open(f"{PREFIX}/{NCLOSE_COUNT_CANDIDATE_PKL}", "wb") as f:
                pkl.dump(nclose_count_candidates, f)
            logging.info(
                f"NClose raw-count VAF filter candidates : {len(nclose_count_candidates)} "
                f"(threshold={args.nclose_count_vaf_threshold})"
            )

            if nclose_count_candidates:
                thread_lim = min(JULIA_BAM_THREAD_LIM, THREAD)
                PROGRESS = ['--progress'] if args.progress else []
                nclose_count_result = subprocess.run(
                    ['python', "-X", f"juliacall-threads={thread_lim}", "-X", "juliacall-handle-signals=yes",
                     os.path.join(os.path.dirname(os.path.abspath(__file__)), '03_Anal_bam.py'),
                     PREFIX, read_bam_loc, CHROMOSOME_INFO_FILE_PATH, main_stat_loc,
                     '--nclose_count_only',
                     '--nclose_count_vaf_threshold', str(args.nclose_count_vaf_threshold)] + PROGRESS
                )
                if nclose_count_result.returncode == 0:
                    keep_pairs, record_count, keep_count = collect_nclose_count_keep_pairs(
                        PREFIX, args.nclose_count_vaf_threshold
                    )
                    if keep_pairs is not None:
                        nclose_nodes, nclose_count_removed = filter_nclose_nodes_by_keep_pairs(
                            nclose_nodes, keep_pairs
                        )
                        logging.info(
                            f"NClose raw-count VAF filter : kept {keep_count}/{record_count} "
                            f"candidate nclose pairs, removed {nclose_count_removed}"
                        )
                else:
                    logging.warning(
                        f"NClose raw-count BAM analysis failed with exit code {nclose_count_result.returncode}; "
                        "skip nclose raw-count VAF filter"
                    )

    base_raw_translocation_candidates = build_raw_translocation_candidates(contig_data, nclose_nodes, chr_len)
    all_raw_candidate_nclose_nodes = convert_all_nclose_comp_to_nclose_nodes(contig_data, all_nclose_comp)
    raw_candidate_nclose_nodes = build_ecdna_nclose_nodes(raw_nclose_nodes, all_raw_candidate_nclose_nodes)
    logging.info(
        f"Raw-read translocation candidate input : "
        f"{sum(len(v) for v in raw_nclose_nodes.values())} raw nclose, "
        f"{sum(len(v) for v in all_raw_candidate_nclose_nodes.values())} all nclose, "
        f"{sum(len(v) for v in raw_candidate_nclose_nodes.values())} merged nclose"
    )
    large_same_chrom_raw_translocation_candidates = build_raw_translocation_candidates(
        contig_data,
        raw_candidate_nclose_nodes,
        chr_len,
        candidate_filter=is_large_same_chrom_raw_candidate,
    )
    raw_translocation_candidates = merge_raw_translocation_candidates(
        base_raw_translocation_candidates,
        large_same_chrom_raw_translocation_candidates,
    )
    with open(f"{PREFIX}/{RAW_TRANSLOCATION_CANDIDATE_PKL}", "wb") as f:
        pkl.dump(raw_translocation_candidates, f)
    logging.info(
        f"Raw-read translocation candidates : {len(raw_translocation_candidates)} "
        f"({len(base_raw_translocation_candidates)} final-nclose, "
        f"{len(large_same_chrom_raw_translocation_candidates)} large same-chrom)"
    )

    if SKIP_BAM_ANAL:
        logging.info("Raw-read translocation BAM analysis skipped")
        raw_virtual_inv_nclose_pairs = set()
    else:
        thread_lim = min(JULIA_BAM_THREAD_LIM, THREAD)
        PROGRESS = ['--progress'] if args.progress else []
        raw_bam_result = subprocess.run(
            ['python', "-X", f"juliacall-threads={thread_lim}", "-X", "juliacall-handle-signals=yes",
             os.path.join(os.path.dirname(os.path.abspath(__file__)), '03_Anal_bam.py'),
             PREFIX, read_bam_loc, CHROMOSOME_INFO_FILE_PATH, main_stat_loc,
             '--skip_nclose_count'] + PROGRESS
        )
        if raw_bam_result.returncode == 0:
            raw_virtual_inv_nclose_pairs = collect_raw_virtual_inv_nclose_pairs(PREFIX)
        else:
            logging.warning(
                f"Raw-read translocation BAM analysis failed with exit code {raw_bam_result.returncode}; "
                "skip raw virtual-inversion nclose removal"
            )
            raw_virtual_inv_nclose_pairs = set()

    if raw_virtual_inv_nclose_pairs:
        nclose_nodes, raw_virtual_inv_removed = remove_nclose_pairs(
            nclose_nodes, raw_virtual_inv_nclose_pairs
        )
        logging.info(
            f"Removed {raw_virtual_inv_removed} raw-read virtual-inversion nclose pairs "
            f"({len(raw_virtual_inv_nclose_pairs)} unique pairs) from final compressed nclose"
        )

    if args.exclude_nclose_list_loc is not None:
        with open(args.exclude_nclose_list_loc, 'r') as f:
            indel_exclude_idx_set = set()
            name_list = []

            for l in f:
                name: str = l.strip()
                if name.startswith('INDEL_INDEX_'):
                    indel_idx = int(name.removeprefix('INDEL_INDEX_'))
                    indel_exclude_idx_set.add(indel_idx)
                else:
                    if nclose_nodes.pop(name, None) is not None:
                        name_list.append(name)
            
            if name_list:
                logging.warning(f"Skipped contig : {', '.join(name_list)}")
            if len(indel_exclude_idx_set) > 0:
                logging.warning(f"Skipped indel index : {', '.join(map(str, indel_exclude_idx_set))}")
    else:
        indel_exclude_idx_set = set()

    with open(f'{PREFIX}/indel_exclude_idx_set.pkl', 'wb') as f:
        pkl.dump(indel_exclude_idx_set, f)

    transloc_nclose_pair_count = write_compressed_nclose_nodes_list(nclose_nodes)

    nclose_node_count = write_nclose_nodes_index(
        f"{PREFIX}/nclose_nodes_index.txt",
        nclose_nodes,
        contig_data,
    )

    logging.info(f"Uncompressed NClose node count : {uncomp_node_count}")    
    logging.info(f"NClose node count : {nclose_node_count}")
    logging.info(f"Telomere connected node count : {telo_node_count}")

    return locals()

def get_ori_ctg_name_data(PAF_FILE_PATH : list) -> list:
    ori_ctg_name_data = []

    for ori_paf_loc in PAF_FILE_PATH:
        ori_ctg_name_list = []
        with open(ori_paf_loc, 'r') as f:
            for paf_line in f:
                paf_line = paf_line.split('\t')
                ori_ctg_name_list.append(paf_line[CTG_NAM])

        ori_ctg_name_data.append(ori_ctg_name_list)

    return ori_ctg_name_data

def import_bed(bed_path: str) -> dict:
    bed_data_file = open(bed_path, "r")
    chr_len = defaultdict(list)
    for curr_data in bed_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]].append((int(curr_data[1]), int(curr_data[2])))
    bed_data_file.close()
    return chr_len

parser = argparse.ArgumentParser(description="Find breakend contigs with contig data and map data")

# 위치 인자 정의
parser.add_argument("paf_file_path", 
                    help="Path to the original PAF file. In VCF input mode, this is the reference-aligned PAF used for telomere/neotelomere anchors.")
parser.add_argument("reference_fai_path", 
                    help="Path to the chromosome information file.")
parser.add_argument("telomere_bed_path", 
                    help="Path to the telomere information file.")
parser.add_argument("repeat_bed_path", 
                    help="Path to the chromosome repeat information file.")                 
parser.add_argument("censat_bed_path", 
                    help="Path to the censat repeat information file.")
parser.add_argument("main_stat_path", 
                    help="Path to the main stat file.")
parser.add_argument("prefix", 
                    help="Pefix for pipeline")
parser.add_argument("read_bam_loc", 
                    help="Raw read alignment bam location")
parser.add_argument("--alt", 
                    help="Path to an alternative PAF file (optional). In VCF mode this is the VCF INS sequence alignment PAF.")
parser.add_argument("--original_paf_loc", nargs='+',
                    help="Original paf location to detect location (primary, alternative paf location)")
parser.add_argument("-t", "--thread", 
                    help="Number of thread", type=int)
parser.add_argument("-d", "--graph_depth", 
                    help="Depth of breakend graph", type=int, default=4)
parser.add_argument("--progress", 
                    help="Show progress bar", action='store_true')
parser.add_argument("--verbose", 
                    help="Enable index, paf output (Could be slow at HDD)", action='store_true')
parser.add_argument("--test", 
                    help=argparse.SUPPRESS, action='store_true')
parser.add_argument("--exclude_nclose_list_loc", 
                    help=argparse.SUPPRESS, default=None)                 
parser.add_argument("--skip_bam_analysis",
                    help="Skip bam analysis", action='store_true')
parser.add_argument("--check_nclose_count",
                    help="Count raw-read support for each nclose and keep only nclose with "
                         "chrA or chrB VAF above --nclose_count_vaf_threshold; "
                         "CEN-SAT-overlapping nclose are kept without VAF filtering.",
                    action='store_true')
parser.add_argument("--nclose_count_vaf_threshold",
                    help="Minimum single-side raw-read VAF to keep an nclose when "
                         "--check_nclose_count is enabled.",
                    type=float, default=NCLOSE_COUNT_DEFAULT_VAF_THRESHOLD)
parser.add_argument("--disable_alt_ctg_simple",
                    help="Disable the default primary-contig rescue that trims telomere-like "
                         "terminal chunks, ignores <=10kb fragments, selects chromosomes covering "
                         "90%% of the remaining span, skips same-direction chromosome-change "
                         "break candidates only when they are within 10kb of an existing nclose "
                         "candidate, and keeps direction-changing chromosome changes unfiltered.",
                    action='store_true')

parser.add_argument("--add_indel_graph",
                    dest="add_indel_graph",
                    help="Add selected type4 indel rescue edges to the breakend graph without increasing graph dimensions. "
                         "In VCF mode, DEL/DUP/Indel-like BND events are eligible and INS is excluded.",
                    action='store_true')
parser.add_argument("--vcf_input",
                    help="VCF input for benchmark mode; bypass PAF-derived nclose discovery and use VCF calls instead.",
                    default=None)
parser.add_argument("--vcf_filter_pass", nargs='+', metavar="FILTER",
                    help="Exact, case-sensitive VCF FILTER values to retain in VCF mode. "
                         "Supplying this option replaces the default PASS and . values.",
                    default=["PASS", "."])
parser.add_argument("--limit_combinations",
                    help="Path to a limit_combinations.json file. Use exactly that "
                         "graph-limit combination and fail instead of trying another "
                         "combination.",
                    default=None)

mode_group = parser.add_mutually_exclusive_group()
mode_group.add_argument("--karyotype_mode",
                        dest="pipeline_mode",
                        help="Use aggressive filtering for karyotype analysis (default).",
                        action="store_const",
                        const=PIPELINE_MODE_KARYOTYPE,
                        default=PIPELINE_MODE_KARYOTYPE)
mode_group.add_argument("--variant_mode",
                        dest="pipeline_mode",
                        help="Use depth-preserving mode for variant/VCF analysis.",
                        action="store_const",
                        const=PIPELINE_MODE_VARIANT)

args = parser.parse_args()

FIXED_LIMIT_COMBINATION = None
if args.limit_combinations is not None:
    try:
        FIXED_LIMIT_COMBINATION = load_limit_combinations(args.limit_combinations)
    except ValueError as exc:
        parser.error(f"--limit_combinations: {exc}")

# t = "02_Build_Breakend_Graph_Limited.py /home/hyunwoo/ACCtools-pipeline/90_skype_run/H1437/20_alignasm/H1437.ctg.aln.paf public_data/chm13v2.0.fa.fai public_data/chm13v2.0_telomere.bed public_data/chm13v2.0_repeat.m.bed public_data/chm13v2.0_censat_v2.1.m.bed /home/hyunwoo/ACCtools-pipeline/90_skype_run/H1437/01_depth/H1437_normalized.win.stat.gz 30_skype_pipe/H1437_23_00_00 /home/hyunwoo/ACCtools-pipeline/90_skype_run/H1437/01_depth/H1437.bam --alt /home/hyunwoo/ACCtools-pipeline/90_skype_run/H1437/20_alignasm/H1437.utg.aln.paf --original_paf_loc /home/hyunwoo/ACCtools-pipeline/90_skype_run/H1437/20_alignasm/H1437.ctg.paf /home/hyunwoo/ACCtools-pipeline/90_skype_run/H1437/20_alignasm/H1437.utg.paf --test --skip_bam_analysis -t 128"
# args = parser.parse_args(t.split()[1:])

PREFIX = args.prefix
REQUESTED_PIPELINE_MODE = PIPELINE_MODE_VARIANT if args.vcf_input is not None else args.pipeline_mode

initial_pipeline_mode_config = make_pipeline_mode_config(
    requested_mode=REQUESTED_PIPELINE_MODE,
    vcf_input=args.vcf_input is not None,
)
logging.info('SKYPE pipeline start')
logging.info(describe_pipeline_mode(initial_pipeline_mode_config))
logging.info("02_Build_Breakend_Graph start")

os.makedirs(PREFIX, exist_ok=True)

LIMIT_COMBINATIONS_OUTPUT_PATH = os.path.join(PREFIX, LIMIT_COMBINATIONS_JSON)
if os.path.isfile(LIMIT_COMBINATIONS_OUTPUT_PATH):
    os.remove(LIMIT_COMBINATIONS_OUTPUT_PATH)

PAF_FILE_PATH = []
if args.vcf_input is not None:
    PAF_FILE_PATH = [f"{PREFIX}/{VCF_SYNTHETIC_PAF_NAME}", args.paf_file_path]
elif args.alt is None:
    PAF_FILE_PATH = [args.paf_file_path]
else:
    PAF_FILE_PATH = [args.paf_file_path, args.alt]

PREPROCESSED_PAF_FILE_PATH = os.path.join(
    PREFIX, f"{os.path.basename(args.paf_file_path)}.ppc.paf"
)
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
REPEAT_INFO_FILE_PATH = args.repeat_bed_path
CENSAT_PATH = args.censat_bed_path
ORIGINAL_PAF_LOC_LIST_ = args.original_paf_loc
main_stat_loc = args.main_stat_path
read_bam_loc = args.read_bam_loc
PRINT_IDX_FILE = args.verbose
THREAD=args.thread
SKIP_BAM_ANAL=args.skip_bam_analysis
CHR_CHANGE_LIMIT_ABS_MAX = args.graph_depth

ORIGINAL_PAF_LOC_LIST = ORIGINAL_PAF_LOC_LIST_

if args.test:
    logging.warning("Test mode is enabled. This mode is for debugging purposes only. The results may be inaccurate and should not be trusted.")
    CHR_CHANGE_LIMIT_ABS_MAX = 2

if (
    FIXED_LIMIT_COMBINATION is not None
    and FIXED_LIMIT_COMBINATION[0] > max(CHR_CHANGE_LIMIT_ABS_MAX, 1)
):
    parser.error(
        "--limit_combinations chromosome-change limit "
        f"{FIXED_LIMIT_COMBINATION[0]} exceeds graph depth "
        f"{CHR_CHANGE_LIMIT_ABS_MAX}"
    )

if args.vcf_input is None:
    if args.vcf_filter_pass != ["PASS", "."]:
        logging.warning("--vcf_filter_pass is ignored without --vcf_input.")
    assert(len(PAF_FILE_PATH) == len(ORIGINAL_PAF_LOC_LIST))
else:
    ORIGINAL_PAF_LOC_LIST_ = []
    ORIGINAL_PAF_LOC_LIST = []

gfa_file_path = []
asm2cov = Counter()
# with open(f'{PREFIX}/asm2cov.pkl', 'rb') as f:
#     gfa_file_path, asm2cov = pkl.load(f)

is_unitig_reduced = False
ori_ctg_name_data = []

if args.vcf_input is not None:
    if args.original_paf_loc is not None:
        logging.warning("--original_paf_loc is ignored in VCF input mode.")
    globals().update(build_vcf_mode_inputs())
else:
    ori_ctg_name_data = get_ori_ctg_name_data(PAF_FILE_PATH)
    telo_coverage = contig_preprocessing_00(PAF_FILE_PATH)
    globals().update(nclose_calc())

telo_connected_node_tuple = extract_telomere_connect_contig_bytuple(PREFIX+"/telomere_connected_list.txt")
chr_fb_len_dict = defaultdict(list)

telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)
telo_dict = defaultdict(list)
for _ in telo_data:
    telo_dict[_[0]].append(_[1:])

telo_fb_dict = defaultdict(list)
for k, v in telo_dict.items():
    for i in v:
        telo_fb_dict[k+i[-1]].append([i[0], i[1]])

for chr_dir, node_id in telo_connected_node_tuple:
    telo_len = 1e9
    for telo_bed in telo_fb_dict[chr_dir]:
        telo_len = min(telo_len, distance_checker_tuple(tuple(telo_bed), (contig_data[node_id][CHR_STR], contig_data[node_id][CHR_END])))
    chr_fb_len_dict[chr_dir].append((node_id, telo_len, chr_dir))

telo_len_data = []
nonzero_telo_set = set()
for chr_dir, telo_len_list in chr_fb_len_dict.items():
    s_telo_len_list = sorted(telo_len_list, key=lambda t: t[1])
    telo_len_data.extend(filter(lambda t: t[1] > 0, s_telo_len_list[1:]))

for idx, len_value, chr in telo_len_data:
    if len_value > 0:
        nonzero_telo_set.add(idx)

if args.vcf_input is None and nclose_node_count > FAIL_NCLOSE_COUNT:
    logging.info("NClose node count is too high.")
    if len(PAF_FILE_PATH) == 1:
        logging.info("No method to reduce nclose node count.")
        sys.exit(1)
    else:
        is_unitig_reduced = True
        logging.info("Retrying with the primary PAF file.")

        PAF_FILE_PATH = [args.paf_file_path, args.paf_file_path]
        ORIGINAL_PAF_LOC_LIST = [ORIGINAL_PAF_LOC_LIST_[0], ORIGINAL_PAF_LOC_LIST_[0]]
        ori_ctg_name_data = get_ori_ctg_name_data(PAF_FILE_PATH)

        telo_coverage = contig_preprocessing_00(PAF_FILE_PATH)
        globals().update(nclose_calc())

        if nclose_node_count > FAIL_NCLOSE_COUNT:
            logging.info("No method to reduce nclose node count.")
            sys.exit(1)

selected_type4_indel_graph_edges = []
type4_indel_zero_dim_edge_set = set()
graph_nclose_nodes = nclose_nodes
canonical_nclose_before_indel_graph = canonical_nclose_snapshot(nclose_nodes)
if args.add_indel_graph:
    type4_indel_graph_candidates = collect_type4_indel_graph_candidates(
        contig_data, df, repeat_censat_data
    )
    selected_type4_indel_graph_edges = select_type4_indel_graph_edges(
        contig_data, nclose_nodes, type4_indel_graph_candidates
    )
    graph_nclose_nodes = augment_nclose_nodes_with_type4_indels(
        nclose_nodes, selected_type4_indel_graph_edges
    )
    type4_indel_zero_dim_edge_set = get_type4_indel_zero_dim_edge_set(
        selected_type4_indel_graph_edges
    )
    selected_type4_indel_kind_by_tuple = {
        edge['type4_tuple']: edge['indel_kind']
        for edge in selected_type4_indel_graph_edges
    }
    selected_type4_indel_kind_count = Counter(selected_type4_indel_kind_by_tuple.values())
    type4_indel_candidate_kind_count = Counter(
        candidate['indel_kind'] for candidate in type4_indel_graph_candidates
    )
    logging.info(
        f'Added {len(type4_indel_zero_dim_edge_set)} type4 indel graph edges '
        f'from {len(selected_type4_indel_kind_by_tuple)} selected type4 indel events '
        f'({dict(selected_type4_indel_kind_count)} selected; '
        f'{len(type4_indel_graph_candidates)} depth-supported candidates '
        f'{dict(type4_indel_candidate_kind_count)} >= {TYPE4_INDEL_GRAPH_MIN_SPAN} bp)'
    )

if canonical_nclose_snapshot(nclose_nodes) != canonical_nclose_before_indel_graph:
    raise AssertionError(
        "--add_indel_graph mutated canonical nclose_nodes; Indel-like edges "
        "must remain graph-only"
    )
if args.vcf_input is not None:
    assert_vcf_nclose_has_no_indel_like(
        contig_data,
        nclose_nodes,
        "post --add_indel_graph canonical nclose",
    )
    if args.add_indel_graph:
        non_vcf_graph_edges = [
            edge
            for edge in selected_type4_indel_graph_edges
            if not str(edge.get('contig_name', '')).startswith(
                VCF_TYPE4_GRAPH_NODE_PREFIX
            )
        ]
        if non_vcf_graph_edges:
            raise AssertionError(
                "VCF --add_indel_graph selected an edge outside graph-only "
                f"VCF nodes: {non_vcf_graph_edges[:2]}"
            )

bnd_graph_adjacency = initialize_bnd_graph(contig_data, graph_nclose_nodes, telo_contig)

with open(f'{PREFIX}/{TYPE4_INDEL_GRAPH_EDGE_PKL}', 'wb') as f:
    pkl.dump(selected_type4_indel_graph_edges, f)

save_loc = PREFIX + '/00_raw'
logging.info("Now saving results in folder : " + PREFIX)

def make_graph(CHR_CHANGE_LIMIT, DIR_CHANGE_LIMIT):
    G = nx.DiGraph()
    for i in bnd_graph_adjacency:
        for j in range(0, CHR_CHANGE_LIMIT+1):
            for k in range(0, DIR_CHANGE_LIMIT+1):
                if type(i)==str:
                    G.add_node((i, j, k))
                else:
                    G.add_node(tuple(list(i)+[j, k]))
    for node in bnd_graph_adjacency:
        for edge in bnd_graph_adjacency[node]:
            if type(edge)==str:
                for j in range(0, CHR_CHANGE_LIMIT+1):
                    for k in range(0, DIR_CHANGE_LIMIT+1):
                        node_limit = tuple(list(node)+[j, k])
                        edge_limit = (edge, j, k)
                        G.add_weighted_edges_from([(node_limit, edge_limit, 0)])
            elif type(node)==str:
                for j in range(0, CHR_CHANGE_LIMIT+1):
                    for k in range(0, DIR_CHANGE_LIMIT+1):
                        node_limit = (node, j, k)
                        edge_limit = tuple(list(edge) + [j, k])
                        G.add_weighted_edges_from([(node_limit, edge_limit, 0)])
            else:
                contig_s = contig_data[node[1]]
                contig_e = contig_data[edge[1]]
                if contig_s[CTG_NAM] == contig_e[CTG_NAM]:
                    if node[1] < edge[1]:
                        d = 0
                        for i in range(node[1]+1, edge[1]):
                            d += contig_data[i][CHR_END] - contig_data[i][CHR_STR]
                    else:
                        d = 0
                        for i in range(edge[1]+1, node[1]):
                            d += contig_data[i][CHR_END] - contig_data[i][CHR_STR]
                            
                else:
                    ds = contig_s[CHR_END] - contig_s[CHR_STR]
                    de = contig_e[CHR_END] - contig_e[CHR_STR]
                    if distance_checker(contig_s, contig_e)==0:
                        d = ds + de - overlap_calculator(contig_s, contig_e)
                    else:  
                        d = distance_checker(contig_s, contig_e) + ds + de 
                is_type4_indel_zero_dim_edge = (
                    (tuple(node), tuple(edge[:2])) in type4_indel_zero_dim_edge_set
                )
                if is_type4_indel_zero_dim_edge:
                    for j in range(0, CHR_CHANGE_LIMIT+1):
                        for k in range(0, DIR_CHANGE_LIMIT+1):
                            node_limit = tuple(list(node)+[j, k])
                            edge_limit = tuple(list(edge)+[j, k])
                            G.add_weighted_edges_from([(node_limit, edge_limit, d)])
                elif contig_s[CHR_NAM] != contig_e[CHR_NAM]:
                    for j in range(0, CHR_CHANGE_LIMIT):
                        for k in range(0, DIR_CHANGE_LIMIT+1):
                            node_limit = tuple(list(node)+[j, k])
                            edge_limit = tuple(list(edge)+[j+1, k])
                            G.add_weighted_edges_from([(node_limit, edge_limit, d)])
                elif contig_s[CTG_DIR] != contig_e[CTG_DIR] and contig_s[CTG_NAM] == contig_e[CTG_NAM]:
                    for j in range(0, CHR_CHANGE_LIMIT+1):
                        for k in range(0, DIR_CHANGE_LIMIT):
                            node_limit = tuple(list(node)+[j, k])
                            edge_limit = tuple(list(edge)+[j, k+1])
                            G.add_weighted_edges_from([(node_limit, edge_limit, d)])
                else:
                    for j in range(0, CHR_CHANGE_LIMIT+1):
                        for k in range(0, DIR_CHANGE_LIMIT+1):
                            node_limit = tuple(list(node)+[j, k])
                            edge_limit = tuple(list(edge)+[j, k])
                            G.add_weighted_edges_from([(node_limit, edge_limit, d)])

    return G

def get_prop_type(value, key=None):
    """
    Performs typing and value conversion for the graph_tool PropertyMap class.
    If a key is provided, it also ensures the key is in a format that can be
    used with the PropertyMap. Returns a tuple, (type name, value, key)
    """

    # Deal with the value
    if isinstance(value, bool):
        tname = 'bool'

    elif isinstance(value, int):
        tname = 'float'
        value = float(value)

    elif isinstance(value, float):
        tname = 'float'

    elif isinstance(value, dict):
        tname = 'object'

    else:
        tname = 'string'
        value = str(value)

    return tname, value, key

def nx2gt(nxG):
    """
    Converts a networkx graph to a graph-tool graph.
    """
    # Phase 0: Create a directed or undirected graph-tool Graph
    gtG = graph_tool.Graph(directed=nxG.is_directed())

    # Add the Graph properties as "internal properties"
    for key, value in nxG.graph.items():
        # Convert the value and key into a type for graph-tool
        tname, value, key = get_prop_type(value, key)

        prop = gtG.new_graph_property(tname) # Create the PropertyMap
        gtG.graph_properties[key] = prop     # Set the PropertyMap
        gtG.graph_properties[key] = value    # Set the actual value

    # Phase 1: Add the vertex and edge property maps
    # Go through all nodes and edges and add seen properties
    # Add the node properties first
    nprops = set() # cache keys to only add properties once
    for node, data in nxG.nodes(data=True):

        # Go through all the properties if not seen and add them.
        for key, val in data.items():
            if key in nprops: continue # Skip properties already added

            # Convert the value and key into a type for graph-tool
            tname, _, key  = get_prop_type(val, key)

            prop = gtG.new_vertex_property(tname) # Create the PropertyMap
            gtG.vertex_properties[key] = prop     # Set the PropertyMap

            # Add the key to the already seen properties
            nprops.add(key)

    # Also add the node id: in NetworkX a node can be any hashable type, but
    # in graph-tool node are defined as indices. So we capture any strings
    # in a special PropertyMap called 'id' -- modify as needed!
    gtG.vertex_properties['id'] = gtG.new_vertex_property('string')

    # Add the edge properties second
    eprops = set() # cache keys to only add properties once
    for src, dst, data in nxG.edges(data=True):

        # Go through all the edge properties if not seen and add them.
        for key, val in data.items():
            if key in eprops: continue # Skip properties already added

            # Convert the value and key into a type for graph-tool
            tname, _, key = get_prop_type(val, key)

            prop = gtG.new_edge_property(tname) # Create the PropertyMap
            gtG.edge_properties[key] = prop     # Set the PropertyMap

            # Add the key to the already seen properties
            eprops.add(key)

    # Phase 2: Actually add all the nodes and vertices with their properties
    # Add the nodes
    vertices = {} # vertex mapping for tracking edges later
    for node, data in nxG.nodes(data=True):

        # Create the vertex and annotate for our edges later
        v = gtG.add_vertex()
        vertices[node] = v

        # Set the vertex properties, not forgetting the id property
        data['id'] = str(node)
        for key, value in data.items():
            gtG.vp[key][v] = value # vp is short for vertex_properties

    # Add the edges
    for src, dst, data in nxG.edges(data=True):

        # Look up the vertex structs from our vertices mapping and add edge.
        e = gtG.add_edge(vertices[src], vertices[dst])

        # Add the edge properties
        for key, value in data.items():
            gtG.ep[key][e] = value # ep is short for edge_properties

    # Done, finally!
    return gtG

def find_vertex_by_id(gtG, node_id):
    for v in gtG.vertices():
        if gtG.vp['id'][v] == str(node_id):
            return v
    return None

def run_graph(data, nonzero_telo_set, CHR_CHANGE_LIMIT, DIR_CHANGE_LIMIT):
    path_list = []
    path_di_list = []

    i, j = data
    src = (chr_rev_corr[i], 0, 0)
    tar = chr_rev_corr[j]
    folder_name = f"{save_loc}/{src[0]}_{tar}"
    cnt = 0
    H = copy.deepcopy(G)

    # src와 tar에 해당하지 않는 노드 제거
    for k in range(contig_data_size, contig_data_size + 2 * CHROMOSOME_COUNT):
        if chr_rev_corr[k] not in (src[0], tar):
            for ii in range(0, CHR_CHANGE_LIMIT + 1):
                for jj in range(0, DIR_CHANGE_LIMIT + 1):
                    H.remove_node((chr_rev_corr[k], ii, jj))

    path_compress = defaultdict(list)
    for ii in range(0, CHR_CHANGE_LIMIT + 1):
        for jj in range(0, DIR_CHANGE_LIMIT + 1):
            I = copy.deepcopy(H)
            if src[0] != tar:
                for iii in range(0, CHR_CHANGE_LIMIT + 1):
                    for jjj in range(0, DIR_CHANGE_LIMIT + 1):
                        if iii or jjj:
                            I.remove_node((src[0], iii, jjj))
                        if iii != ii or jjj != jj:
                            I.remove_node((tar, iii, jjj))
            else:
                for iii in range(0, CHR_CHANGE_LIMIT + 1):
                    for jjj in range(0, DIR_CHANGE_LIMIT + 1):
                        if (iii != ii or jjj != jj) and (iii + jjj != 0):
                            if I.has_node((tar, iii, jjj)):
                                I.remove_node((tar, iii, jjj))

            # networkx 그래프 I를 graph-tool 그래프로 변환
            F = nx2gt(I)

            # graph-tool에서 'id' 속성을 이용해 source와 target vertex 찾기
            s = find_vertex_by_id(F, src)
            t = find_vertex_by_id(F, (tar, ii, jj))

            if gt.shortest_distance(F, source=s, target=t) == float('inf'):
                continue

            for path_vertices in gt.all_paths(F, s, t):
                path = [eval(F.vp['id'][v]) for v in path_vertices]

                path_len = len(path)
                if path_len == 1:
                    continue

                # path_counter 계산 (contig_data 관련 기존 로직)
                path_counter = Counter()
                if contig_data[path[1][1]][CTG_NAM] == contig_data[path[2][1]][CTG_NAM]:
                    path_counter[contig_data[path[1][1]][CHR_NAM]] += \
                        contig_data[path[1][1]][CHR_END] - contig_data[path[1][1]][CHR_STR]

                if contig_data[path[path_len - 2][1]][CTG_NAM] == contig_data[path[path_len - 3][1]][CTG_NAM]:
                    path_counter[contig_data[path[path_len - 2][1]][CHR_NAM]] += \
                        contig_data[path[path_len - 2][1]][CHR_END] - contig_data[path[path_len - 2][1]][CHR_STR]
                
                contig_set = Counter()
                contig_set[path[path_len - 2][1]] += 1

                nclose_st_cluster_set = set()
                nclose_ed_cluster_set = set()
                cluster_overlap_flag = False
                for node in range(1, path_len - 1):
                    if path[node][1] in st_compress and st_compress[path[node][1]] in nclose_st_cluster_set:
                        cluster_overlap_flag = True
                        break
                    elif path[node][1] in st_compress:
                        nclose_st_cluster_set.add(st_compress[path[node][1]])
                    if path[node][1] in ed_compress and ed_compress[path[node][1]] in nclose_ed_cluster_set:
                        cluster_overlap_flag = True
                        break
                    elif path[node][1] in ed_compress:
                        nclose_ed_cluster_set.add(ed_compress[path[node][1]])
                if cluster_overlap_flag:
                    continue

                censat_node_vis = 0
                contig_overuse_flag = False
                censat_overuse_flag = False

                for node in range(1, path_len - 2):
                    contig_set[path[node][1]] += 1
                    if contig_set[path[node][1]] >= BND_OVERUSE_CNT:
                        contig_overuse_flag = True
                        break
                    
                    contig_s = contig_data[path[node][1]]
                    contig_e = contig_data[path[node + 1][1]]
                    ds = contig_s[CHR_END] - contig_s[CHR_STR]
                    de = contig_e[CHR_END] - contig_e[CHR_STR]
                    if contig_s[CTG_CENSAT] != '0' and contig_e[CTG_CENSAT] != '0':
                        censat_node_vis += 1
                        if censat_node_vis >= CENSAT_VISIT_LIMIT:
                            censat_overuse_flag = True
                            break
                    if contig_s[CTG_NAM] == contig_e[CTG_NAM]:
                        if path[node][1] < path[node + 1][1]:
                            for same_contig in range(path[node][1] + 1, path[node + 1][1]):
                                path_counter[contig_data[same_contig][CHR_NAM]] += \
                                    contig_data[same_contig][CHR_END] - contig_data[same_contig][CHR_STR]
                        else:
                            for same_contig in range(path[node + 1][1] + 1, path[node][1]):
                                path_counter[contig_data[same_contig][CHR_NAM]] += \
                                    contig_data[same_contig][CHR_END] - contig_data[same_contig][CHR_STR]
                    else:
                        if distance_checker(contig_s, contig_e) == 0:
                            path_counter[contig_s[CHR_NAM]] += \
                                ds + de - overlap_calculator(contig_s, contig_e)
                        else:
                            ds = contig_s[CHR_END] - contig_s[CHR_STR]
                            de = contig_e[CHR_END] - contig_e[CHR_STR]
                            path_counter[contig_s[CHR_NAM]] += \
                                distance_checker(contig_s, contig_e) + ds + de
                if contig_overuse_flag or censat_overuse_flag:
                    continue

                total_path_ref_len = sum(path_counter.values())
                ig_k_list = [k for k, v in path_counter.items() if v < IGNORE_PATH_LIMIT]
                for k in ig_k_list:
                    del path_counter[k]
                if len(path_counter) == 0:
                    continue

                ack = sorted(path_counter.items(), key=lambda x: -x[1])
                longest_chr = ack[0]
                second_chr = ack[1] if len(ack) > 1 else ack[0]

                flag = True
                flagflag = False
                rank = []
                _rank_norm = lambda c: 'chrX' if c == 'chrY' else c
                for i in range(min(PATH_MAJOR_COMPONENT, len(ack))):
                    if _rank_norm(ack[i][0]) in (_rank_norm(src[0][:-1]), _rank_norm(tar[:-1])):
                        rank.append(i)
                if len(rank) == 2:
                    flagflag = True
                elif len(rank) == 1:
                    if rank[0] <= 1:
                        flagflag = True
                    else:
                        flagflag = False
                else:
                    flagflag = False
                flag = flagflag
                if (longest_chr[1] + second_chr[1]) / total_path_ref_len < 0.5:
                    flag = False
                if chr_len[longest_chr[0]] + chr_len[second_chr[0]] < total_path_ref_len:
                    flag = False
                if total_path_ref_len < MIN_PATH_REF_LEN:
                    flag = False
                if flag:
                    key = tuple(sorted(path_counter.keys()))
                    flagflag = True
                    for curr_path_counter in path_compress[key]:
                        check = all(
                            abs(curr_path_counter[chr_name] - path_counter[chr_name]) <= PATH_COMPRESS_LIMIT 
                            for chr_name in key
                        )
                        
                        if check:
                            flagflag = False
                            break
                    if flagflag:
                        path_compress[key].append(copy.deepcopy(path_counter))
                        cnt += 1

                        if PRINT_IDX_FILE:
                            os.makedirs(folder_name, exist_ok=True)
                            file_name = folder_name + f"/{cnt}.paf"
                            file_name2 = folder_name + f"/{cnt}.index.txt"
                            with open(file_name, "wt") as f, open(file_name2, "wt") as g:
                                for nodes in path:
                                    if isinstance(nodes[0], str):
                                        print(nodes, file=f)
                                        print(nodes, file=g)
                                    else:
                                        f.write("\t".join(map(str, contig_data[nodes[1]])) + "\n")
                                        g.write("\t".join(map(str, nodes)) + "\n")
                                print(ack, file=f)
                                print(ack, file=g)
                        neotelochk_1 = 0
                        neotelochk_2 = 0
                        double_neotelo = 0
                        if path[1][1] in nonzero_telo_set:
                            neotelochk_1 = 1
                        if path[-2][1] in nonzero_telo_set:
                            neotelochk_2 = 1
                        if neotelochk_1>0 and neotelochk_2>0:
                            double_neotelo = 1
                        path_list.append((path, ack))
                        path_di_list.append(ii+jj+neotelochk_1+neotelochk_2+double_neotelo) # Score edit

                        if cnt >= PAT_PATH_LIMIT:
                            return ((src[0], tar), cnt, path_list, path_di_list)

    return ((src[0], tar), cnt, path_list, path_di_list)

def init_worker(shared_graph):
    global G
    G = shared_graph

def run_graph_pipeline():
    tar_ind_list = []
    for i in range(contig_data_size, contig_data_size + 2*CHROMOSOME_COUNT):
        for j in range(i + 1, contig_data_size + 2*CHROMOSOME_COUNT):
            tar_ind_list.append((i, j))

    if FIXED_LIMIT_COMBINATION is not None:
        limit_combinations = [FIXED_LIMIT_COMBINATION]
        idx = 0
        logging.info(
            "Using fixed limit_combinations from "
            f"{args.limit_combinations}: {FIXED_LIMIT_COMBINATION}"
        )
    else:
        limit_combinations = [
            (c, DIR_CHANGE_LIMIT_ABS_MAX)
            for c in range(CHR_CHANGE_LIMIT_ABS_MAX, 0, -1)
        ]
        limit_combinations.append((1, 0))

        if nclose_node_count > HARD_NCLOSE_COUNT:
            hard_start_chr = min(CHR_CHANGE_LIMIT_HARD_START, CHR_CHANGE_LIMIT_ABS_MAX)
            idx = limit_combinations.index((hard_start_chr, DIR_CHANGE_LIMIT_ABS_MAX))
        else:
            idx = 0

    last_success = None

    while 0 <= idx < len(limit_combinations):
        CHR_CHANGE_LIMIT_PREFIX, DIR_CHANGE_LIMIT_PREFIX = limit_combinations[idx]
        path_list_dict_data = []
        path_di_list_dict_data = []
        tot_cnt = 0
        cnt_list = []

        G_obj = make_graph(
            CHR_CHANGE_LIMIT=CHR_CHANGE_LIMIT_PREFIX,
            DIR_CHANGE_LIMIT=DIR_CHANGE_LIMIT_PREFIX
        )

        with Pool(
            processes=THREAD,
            initializer=init_worker,
            initargs=(G_obj,)
        ) as pool:
            result_iterator = pool.imap_unordered(
                partial(
                    run_graph,
                    nonzero_telo_set=nonzero_telo_set,
                    CHR_CHANGE_LIMIT=CHR_CHANGE_LIMIT_PREFIX,
                    DIR_CHANGE_LIMIT=DIR_CHANGE_LIMIT_PREFIX
                ),
                tar_ind_list
            )

            for rs in tqdm(
                result_iterator,
                total=len(tar_ind_list),
                desc=f'Build breakend graph ({CHR_CHANGE_LIMIT_PREFIX},{DIR_CHANGE_LIMIT_PREFIX})',
                disable=not sys.stdout.isatty() and not args.progress
            ):
                cnt_list.append((rs[0], rs[1]))
                tot_cnt += rs[1]

                if tot_cnt >= TOT_PATH_LIMIT:
                    pool.terminate()
                    shutil.rmtree(save_loc, ignore_errors=True)
                    break

                path_list_dict_data.append((f'{rs[0][0]}_{rs[0][1]}', rs[2]))
                path_di_list_dict_data.append((f'{rs[0][0]}_{rs[0][1]}', rs[3]))

        del G_obj
        success = (tot_cnt < TOT_PATH_LIMIT)

        if success:
            last_success = {
                "chr":    CHR_CHANGE_LIMIT_PREFIX,
                "dir":    DIR_CHANGE_LIMIT_PREFIX,
                "paths":  path_list_dict_data,
                "paths_di": path_di_list_dict_data,
                "cnts":   cnt_list
            }
            
            path_count = sum(i[1] for i in cnt_list)
            logging.info(f'SUCCESS at {(CHR_CHANGE_LIMIT_PREFIX, DIR_CHANGE_LIMIT_PREFIX)}, with {path_count} paths')
            
            # Simple estimated path count 
            if path_count >= TOT_PATH_LIMIT / 2:
                logging.info('FAIL estimated with next stage.')
                break
            idx -= 1
        else:
            logging.info(f'FAIL at {(CHR_CHANGE_LIMIT_PREFIX, DIR_CHANGE_LIMIT_PREFIX)}')
            if FIXED_LIMIT_COMBINATION is not None:
                logging.info(
                    "Fixed limit_combinations failed; fallback is disabled."
                )
                break
            if last_success is not None:
                break

            idx += 1
    return last_success


last_success = run_graph_pipeline()
if not last_success:
    logging.info('Breakend graph is too divergent.')
    if (
        FIXED_LIMIT_COMBINATION is not None
        or args.vcf_input is not None
        or is_unitig_reduced
        or len(PAF_FILE_PATH) == 1
    ):
        logging.info('Breakend path failed.')
        sys.exit(1)
    else:
        is_unitig_reduced = True
        logging.info("Retrying with the primary PAF file.")

        PAF_FILE_PATH = [args.paf_file_path, args.paf_file_path]
        ORIGINAL_PAF_LOC_LIST = [ORIGINAL_PAF_LOC_LIST_[0], ORIGINAL_PAF_LOC_LIST_[0]]
        ori_ctg_name_data = get_ori_ctg_name_data(PAF_FILE_PATH)

        telo_coverage = contig_preprocessing_00(PAF_FILE_PATH)
        globals().update(nclose_calc())

        if nclose_node_count > FAIL_NCLOSE_COUNT:
            logging.info("NClose node count is too high.")
            logging.info("No method to reduce nclose node count.")
            sys.exit(1)
        
        last_success = run_graph_pipeline()
        if not last_success:
            logging.info('Breakend graph is too divergent.')
            logging.info('Breakend path failed.')
            sys.exit(1) 
        


CHR_CHANGE_LIMIT_PREFIX = last_success["chr"]
DIR_CHANGE_LIMIT_PREFIX = last_success["dir"]

path_list_dict_data, path_di_list_dict_data = last_success["paths"], last_success["paths_di"]
cnt_list = last_success["cnts"]

path_list_dict = dict()
for k, v in path_list_dict_data:
    path_list_dict[k] = v

path_di_list_dict = dict()
for k, v in path_di_list_dict_data:
    path_di_list_dict[k] = v

all_nclose_nodes = convert_all_nclose_comp_to_nclose_nodes(contig_data, all_nclose_comp)
ecdna_nclose_nodes = build_ecdna_nclose_nodes(raw_nclose_nodes, all_nclose_nodes)
logging.debug(
    f"ecDNA nclose input: {sum(len(v) for v in raw_nclose_nodes.values())} filtered nclose, "
    f"{sum(len(v) for v in all_nclose_nodes.values())} all-file nclose, "
    f"{sum(len(v) for v in ecdna_nclose_nodes.values())} merged nclose"
)

inversion_adjacency = initialize_inversion_only_graph(contig_data, ecdna_nclose_nodes)

inversion_graph = nx2gt(make_inversion_nx_graph(inversion_adjacency))

ecdna_circuit_candidate_set = set()

for circuit in gt.all_circuits(inversion_graph, max_length=4):
    if len(circuit) == 4:
        circuit_node_list = []
        for idx in circuit:
            v = inversion_graph.vertex(idx)
            node = int(inversion_graph.vp['id'][v][1:-1].split(", ")[-1])
            circuit_node_list.append(node)
        ecdna_circuit_candidate_set.add(tuple(sorted(circuit_node_list)))

ecdna_circuit_set = set()

for circuit in ecdna_circuit_candidate_set:
    if circuit_length_calculator(circuit) < CIRCUIT_ECDNA_LENGTH_LIMIT:
        ecdna_circuit_set.add(circuit)


# Todo : transform into vcf format
with open(f"{PREFIX}/ecdna_circuit_data.pkl", "wb") as f:
    pkl.dump((list(ecdna_circuit_set), ecdna_nclose_nodes), f)

logging.info(f'Final success settings: {(CHR_CHANGE_LIMIT_PREFIX, DIR_CHANGE_LIMIT_PREFIX)}')
write_limit_combinations(
    LIMIT_COMBINATIONS_OUTPUT_PATH,
    (CHR_CHANGE_LIMIT_PREFIX, DIR_CHANGE_LIMIT_PREFIX),
)

cancer_prefix = os.path.basename(PREPROCESSED_PAF_FILE_PATH).split('.')[0]

pipeline_mode_config = save_pipeline_mode(
    PREFIX,
    requested_mode=REQUESTED_PIPELINE_MODE,
    nclose_node_count=nclose_node_count,
    nclose_limit=HARD_NCLOSE_COUNT,
    vcf_input=args.vcf_input is not None,
    vcf_input_path=os.path.abspath(args.vcf_input) if args.vcf_input is not None else None,
)

with open(f'{PREFIX}/path_data.pkl', 'wb') as f:
    pkl.dump(path_list_dict, f)

with open(f'{PREFIX}/path_di_data.pkl', 'wb') as f:
    pkl.dump(path_di_list_dict, f)

with open(f'{PREFIX}/paf_file_path.pkl', 'wb') as f:
    pkl.dump(PAF_FILE_PATH, f)

with open(f'{PREFIX}/report.txt', 'w') as f:
    cnt = sum(i[1] for i in cnt_list)
    print(cancer_prefix, file=f)
    print(cnt, file=f)
    logging.info(f"Total path count : {cnt}")
    for (st, nd), c in sorted(cnt_list, key=lambda x:-x[1]):
        if c > 0:
            print(st, nd, c, file=f)

if args.vcf_input is not None:
    assert_vcf_nclose_has_no_indel_like(
        contig_data,
        nclose_nodes,
        "final nclose_chunk_data.pkl",
    )

with open(f'{PREFIX}/nclose_chunk_data.pkl', 'wb') as f:
    pkl.dump((nclose_nodes, st_compress, ed_compress), f)
