import sys
import argparse
import subprocess

import shutil
import os
import pickle as pkl

import pandas as pd
import numpy as np

import itertools
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

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info('SKYPE pipeline start')
logging.info("02_Build_Breakend_Graph start")

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
CHROMOSOME_COUNT = 23
K = 1000
M = 1000 * K
CHUKJI_LIMIT = -1
BND_CONTIG_BOUND = 0.1

TOT_PATH_LIMIT = 3*M
PAT_PATH_LIMIT = 10*K

CHR_CHANGE_LIMIT_ABS_MAX = 5
DIR_CHANGE_LIMIT_ABS_MAX = 1
CENSAT_VISIT_LIMIT = 2

CHR_CHANGE_LIMIT_HARD_START = 2
HARD_NCLOSE_COUNT = 1000
FAIL_NCLOSE_COUNT = 10000

BND_OVERUSE_CNT = 2
PATH_MAJOR_COMPONENT = 3
NCLOSE_COMPRESS_LIMIT = 100*K
NCLOSE_MERGE_LIMIT = 1*K
ALL_REPEAT_NCLOSE_COMPRESS_LIMIT = 500*K
PATH_COMPRESS_LIMIT = 50*K
IGNORE_PATH_LIMIT = 50*K
NON_REPEAT_NOISE_RATIO=0.1

CONTIG_MINIMUM_SIZE = 100*K
BND_CONTIG_BOUND = 0.1
RPT_BND_CONTIG_BOUND = 0.2
CHROMOSOME_COUNT = 23
MAPQ_BOUND = 60
TELOMERE_EXPANSION = 5 * K
TELOMERE_COMPRESS_RANGE = 100*K
CENSAT_COMPRESSABLE_THRESHOLD = 1000*K

FORCE_TELOMERE_THRESHOLD = 10*K
TELOMERE_CLUSTER_THRESHOLD = 500*K
SUBTELOMERE_LENGTH = 500*K


MIN_FLANK_SIZE_BP = 1*M

WEAK_BREAKEND_CEN_RATIO_THRESHOLD = 1.3
BREAKEND_CEN_RATIO_THRESHOLD = 1.5

REPEAT_MERGE_GAP = 0

MAX_OVERLAP_SCORE = 3

SPLIT_CTG_LEN_LIMIT = 100 * K

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


def import_repeat_data_00(file_path : str) -> dict :
    logging.info("Contig preprocessing start")
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
    
def overlap_calculator(node_a : tuple, node_b : tuple) -> int :
    return min(abs(node_a[CHR_END] - node_b[CHR_STR]), abs(node_b[CHR_END] - node_a[CHR_STR]))

def inclusive_checker_tuple(tuple_a : tuple, tuple_b : tuple) -> bool :
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False

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

def check_censat_contig(all_repeat_censat_con : set, ALIGNED_PAF_LOC_LIST : list, ORIGNAL_PAF_LOC_LIST : list, contig_data : list):
    div_repeat_paf_name = div_repeat_paf(ORIGNAL_PAF_LOC_LIST, ALIGNED_PAF_LOC_LIST, contig_data)
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

def label_repeat_node(contig_data: list, repeat_data) -> list:
    labels = []
    ends_map = {
        chrom: [iv[1] for iv in intervals]
        for chrom, intervals in repeat_data.items()
    }

    for contig in contig_data:
        chrom = contig[CHR_NAM]
        intervals = repeat_data.get(chrom, [])
        ends = ends_map.get(chrom, [])
        if not intervals:
            labels.append(('0','0'))
            continue
        c_start, c_end = contig[7], contig[8]
        idx = bisect.bisect_left(ends, c_start)
        if idx < len(intervals):
            iv_start, iv_end = intervals[idx]
            if distance_checker(contig, (0,0,0,0,0,0,0, iv_start, iv_end)) == 0:
                suffix = "in" if inclusive_checker_node(contig, (0,0,0,0,0,0,0, iv_start, iv_end)) else ""
                labels.append((chrom, "r" + suffix))
                continue
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
    for i in range(contig_data_size, contig_data_size + CHROMOSOME_COUNT*2):
        curr_telo_set = set()
        now_telo = chr_rev_corr[i]
        flag = False
        mini_telo_dist = INF
        for j in range(2):
            for _ in adjacency[j][i]:
                # 10K 이내면 없애기
                if i < contig_data_size + CHROMOSOME_COUNT:
                    mini_telo_dist = min(mini_telo_dist, 0 if contig_data[_[1]][CHR_STR] < telo_data[now_telo][0] else contig_data[_[1]][CHR_STR] - telo_data[now_telo][0])
                    if contig_data[_[1]][CHR_STR] < telo_data[now_telo][0]+FORCE_TELOMERE_THRESHOLD:
                        flag = True
                else:
                    mini_telo_dist = min(mini_telo_dist, 0 if contig_data[_[1]][CHR_END] > telo_data[now_telo][1] else telo_data[now_telo][1] - contig_data[_[1]][CHR_END])
                    if contig_data[_[1]][CHR_END] > telo_data[now_telo][1]-FORCE_TELOMERE_THRESHOLD:
                        flag = True
                curr_telo_set.add(_[1])
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
                if contig_data[st][CHR_NAM] == now_telo[0:-1] \
                and st not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[st], temp_contig) \
                and contig_data[st][CTG_DIR]=='+':
                    telo_connect_node = st
                    telo_dist = telo_distance_checker(contig_data[st], temp_contig)
                    telo_dir = '+'
                if contig_data[ed][CHR_NAM] == now_telo[0:-1] \
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
                if contig_data[st][CHR_NAM] == now_telo[0:-1] \
                and st not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[st], temp_contig) \
                and contig_data[st][CTG_DIR]=='-':
                    telo_connect_node = st
                    telo_dist = telo_distance_checker(contig_data[st], temp_contig)
                    telo_dir = '-'
                if contig_data[ed][CHR_NAM] == now_telo[0:-1] \
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

def edge_optimization(contig_data : list, contig_adjacency : list, telo_dict : dict, asm2cov : dict) -> tuple :
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

    return telo_connected_set, telo_connected_dict, telo_connected_graph_dict, telo_coverage

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
    chrM_flag = True
    idx = 0
    cnt = 0
    len_count = Counter()
    for i in range(1, contig_data_size+1):
        len_count[contig_data[i-1][CTG_NAM]]+=contig_data[i-1][CHR_END]-contig_data[i-1][CHR_STR]
        cnt+=1
        if (i-1) in telo_connect_info:
            is_telo = True
        if contig_data[i-1][CHR_NAM]!='chrM':
            chrM_flag = False
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
            chrM_flag = True
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
    chrM_flag = True
    idx = 0
    cnt = 0
    len_count = Counter()
    for i in range(1, contig_data_size+1):
        len_count[contig_data[i-1][CTG_NAM]]+=contig_data[i-1][CHR_END]-contig_data[i-1][CHR_STR]
        cnt+=1
        if (i-1) in telo_connect_info:
            is_telo = True
        if contig_data[i-1][CHR_NAM]!='chrM':
            chrM_flag = False
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
            if checker==3:
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
            chrM_flag = True
            if i < contig_data_size and repeat_label[i][0]!='0':
                is_front_back_repeat = True
    contig_data = contig_data[:-1]
    return [using_contig_list, contig_type, contig_terminal_node, len_count]

def alt_preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list, telo_connect_info : set, telo_dict : dict) -> list :
    checker = 0
    contig_data_size = len(contig_data)-1
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    is_front_back_repeat = False
    chrM_flag = True
    len_count = Counter()
    idx = 0
    cnt = 0
    telo_node_count = 0
    
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
        if contig_data[i-1][CHR_NAM] != 'chrM':
            chrM_flag = False
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
            if (checker>0 and checker != 3 and checker != 5):
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
            cnt = 0
            checker = 0
            telo_node_count = 0
            is_telo = False
            is_front_back_repeat = False
            chrM_flag = True
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


def find_breakend_centromere(repeat_censat_data : dict, chr_len : dict, df : pd.DataFrame):
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

            if rep_start_0 == 0 or rep_end_0 == chrom_length - 1:
                continue
            
            final_weighted_ratio = -1
            final_left_weighted = -1
            final_right_weighted = -1

            weak_cnt = 0
            for FLANK_SIZE_BP in [5 * M, 10 * M, 20 * M]:
                # 좌측 flanking: repeat의 1-indexed 시작은 rep_start_0 + 1
                if rep_start_0 > 0:
                    left_flank_end = rep_start_0  # repeat 시작 전 마지막 base (1-indexed)
                    left_flank_start = max(1, (rep_start_0 + 1) - FLANK_SIZE_BP)

                    if left_flank_end < left_flank_start or (left_flank_end - left_flank_start + 1) < MIN_FLANK_SIZE_BP:
                        left_flank_start = None
                        left_flank_end = None
                        left_weighted = None
                    else:
                        left_weighted = weighted_avg_meandepth(chrom_df, left_flank_start, left_flank_start + MIN_FLANK_SIZE_BP)
                else:
                    left_flank_start = None
                    left_flank_end = None
                    left_weighted = None

                # 우측 flanking: repeat의 끝 이후 첫 base부터
                if rep_end_0 < chrom_length:
                    right_flank_start = rep_end_0 + 1
                    right_flank_end = min(chrom_length, rep_end_0 + FLANK_SIZE_BP)

                    if right_flank_end < right_flank_start or (right_flank_end - right_flank_start + 1) < MIN_FLANK_SIZE_BP:
                        right_flank_start = None
                        right_flank_end = None
                        right_weighted = None
                    else:
                        right_weighted = weighted_avg_meandepth(chrom_df, right_flank_end - MIN_FLANK_SIZE_BP, right_flank_end)
                else:
                    right_flank_start = None
                    right_flank_end = None
                    right_weighted = None

                # 양쪽 flanking 영역 중 하나라도 존재하지 않으면 해당 repeat를 건너뜁니다.
                if left_flank_start is None or right_flank_start is None:
                    continue
                
                if left_weighted >= right_weighted:
                    weighted_ratio = left_weighted / right_weighted
                else:
                    weighted_ratio = right_weighted / left_weighted
                
                if final_weighted_ratio < weighted_ratio:
                    final_weighted_ratio = weighted_ratio
                    final_left_weighted = left_weighted
                    final_right_weighted = right_weighted

                if WEAK_BREAKEND_CEN_RATIO_THRESHOLD < weighted_ratio:
                    weak_cnt += 1

            if final_weighted_ratio == -1:
                continue
            
            baseline = WEAK_BREAKEND_CEN_RATIO_THRESHOLD if weak_cnt >= 2 else BREAKEND_CEN_RATIO_THRESHOLD 

            results.append({
                'chr': chrom,
                'repeat_start_0': rep_start_0,
                'repeat_end_0': rep_end_0,
                'left_weighted_meandepth': final_left_weighted,
                'right_weighted_meandepth': final_right_weighted,
                'weighted_ratio': final_weighted_ratio,
                'baseline' : baseline
            })

    result_df = pd.DataFrame(results)
    result_df = result_df[result_df['weighted_ratio']>=result_df['baseline']]

    censat_bnd_chr_list = sorted(set(result_df['chr']), key=lambda t : chr2int(t))
    logging.info(f'Breakend censat chr : {" ".join(censat_bnd_chr_list)}')

    right_df = result_df[result_df['right_weighted_meandepth'] > result_df['left_weighted_meandepth']]
    left_df = result_df[result_df['left_weighted_meandepth'] > result_df['right_weighted_meandepth']]
    vtg_list = []
    prefix = "virtual_censat_contig"
    cnt = 0
    for row1, row2 in itertools.combinations(right_df.itertuples(index = False), 2):
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


    for row1, row2 in itertools.combinations(left_df.itertuples(index = False), 2):
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

    for row1 in left_df.itertuples(index = False):
        for row2 in right_df.itertuples(index = False):
            cnt+=1
            N = int(CENSAT_COMPRESSABLE_THRESHOLD//2)
            mid_row1_censat = int(row1.repeat_start_0 + row1.repeat_end_0)//2
            mid_row2_censat = int(row2.repeat_start_0 + row2.repeat_end_0)//2

            temp_node1 = [f'{prefix}_{cnt}', N, 0, N//2, '+', row1.chr, 
                     chr_len[row1.chr], mid_row1_censat - N//2, mid_row1_censat, 
                     60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '+', row1.chr, f"2.{(cnt-1)*2}"]
        
            temp_node2 = [f'{prefix}_{cnt}', N, N//2, N, '+', row2.chr,
                     chr_len[row2.chr], mid_row2_censat, mid_row2_censat + N//2, 
                     60, 1, (cnt-1)*2, (cnt-1)*2+1, 0, 0, 0, 0, 0, 0, '-', row1.chr, f"2.{(cnt-1)*2+1}"]
            
            vtg_list.append(temp_node1)
            vtg_list.append(temp_node2)

    return vtg_list

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


def pass_pipeline(pre_contig_data, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, telo_ppc_passed):
    if not telo_ppc_passed:
        if len(pre_contig_data)==0:
            return []
        contig_data = []

        for i in pre_contig_data:
            contig_data.append(i[:9] + [i[CTG_MAPQ], i[CTG_GLOBALIDX]])

        node_label = label_node(contig_data, telo_dict)

        repeat_label = label_repeat_node(contig_data, repeat_data)

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

        
        new_node_repeat_label = label_repeat_node(new_contig_data, repeat_data)
        new_node_repeat_censat_label = label_repeat_node(new_contig_data, repeat_censat_data)
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

    final_contig_repeat_label = label_repeat_node(final_contig, repeat_data)

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
                        censat_contig_name : set, repeat_censat_data : dict, ALIGNED_PAF_LOC_LIST : list, ORIGNAL_PAF_LOC_LIST : list, \
                        telo_set : set, chr_len : dict, asm2cov : dict) -> tuple:
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

    div_repeat_paf_name = div_repeat_paf(ORIGNAL_PAF_LOC_LIST, ALIGNED_PAF_LOC_LIST, contig_data)

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
                while st <= e and [contig_data[st][CTG_DIR], contig_data[st][CHR_NAM]] == st_chr:
                    st+=1
                    cut_ratio, _ = calculate_single_contig_ref_ratio(contig_data[s:st+1])
                    if abs(cut_ratio-1) > BND_CONTIG_BOUND:
                        break
                st-=1
                while ed >= s and [contig_data[ed][CTG_DIR], contig_data[ed][CHR_NAM]] == ed_chr:
                    ed-=1
                    cut_ratio, _ = calculate_single_contig_ref_ratio(contig_data[ed:e+1])
                    if abs(cut_ratio-1) > BND_CONTIG_BOUND:
                        break
                ed+=1

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
                if max_chr not in ((contig_data[s][CTG_DIR], contig_data[s][CHR_NAM]), (contig_data[e][CTG_DIR], contig_data[e][CHR_NAM])) \
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
                    st = nclose[0]
                    ed = nclose[2]
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
                    if nclose_front_const + nclose_back_const >= 3:
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
                            and distance_checker(contig_data[ed], dummy_list + i[2]) < compress_limit:
                                nclose_coverage[(i[3], i[4])] += asm2cov[cov_count_name]
                                nclose_compress_track[(i[3], i[4])].append((st, ed))
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
                            and distance_checker(contig_data[ed], dummy_list + i[1]) < compress_limit:
                                nclose_coverage[(i[3], i[4])] += asm2cov[cov_count_name]
                                nclose_compress_track[(i[3], i[4])].append((st, ed))
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
                                    nclose_coverage[(i[3], i[4])] += asm2cov[cov_count_name]
                                    nclose_compress_track[(i[3], i[4])].append((st, ed))
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
                                    nclose_coverage[(i[3], i[4])] += asm2cov[cov_count_name]
                                    nclose_compress_track[(i[3], i[4])].append((st, ed))
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
    
def initialize_bnd_graph(contig_data : list, nclose_nodes : dict, telo_contig : dict) -> dict:
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
                            if dir1 == "inc":
                                if curr_contig[CTG_DIR]=='+' and k_ind=='f+': # +
                                    bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                    bnd_adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                                elif curr_contig[CTG_DIR]=='-' and k_ind=='b-':
                                    bnd_adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                                    bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                elif curr_contig[CTG_DIR]=='+' and k_ind=='b-':
                                    bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                    bnd_adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                                elif curr_contig[CTG_DIR]=='-' and k_ind=='f+': # +
                                    bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                    bnd_adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                            else:
                                if curr_contig[CTG_DIR]=='+' and k_ind=='b+': # +
                                    bnd_adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                                    bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                elif curr_contig[CTG_DIR]=='-' and k_ind=='f-':
                                    bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                    bnd_adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                                elif curr_contig[CTG_DIR]=='+' and k_ind=='f-':
                                    bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                    bnd_adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                                elif curr_contig[CTG_DIR]=='-' and k_ind=='b+': # +
                                    bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                    bnd_adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                        # else: each end have different sign of increment
                        # -> one node is included in the other
                        # thus, we should determine by mean val.
                        else:
                            if k_ind[0]=='f':
                                bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_FOR, nclose_idx])
                                bnd_adjacency[(DIR_BAK, nclose_idx)].append([DIR_IN, telo_idx])
                            else:
                                bnd_adjacency[(DIR_OUT, telo_idx)].append([DIR_BAK, nclose_idx])
                                bnd_adjacency[(DIR_FOR, nclose_idx)].append([DIR_IN, telo_idx])
                            

    
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

def contig_preprocessing_00(PAF_FILE_PATH_ : list):

    original_node_count = 0

    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
    telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)
    repeat_data = import_repeat_data_00(REPEAT_INFO_FILE_PATH)
    repeat_censat_data = import_censat_repeat_data(CENSAT_PATH)


    df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
    df = df.query('chr != "chrM"')

    cen_vtg_contig = find_breakend_centromere(repeat_censat_data, chr_len, df)

    telo_dict = defaultdict(list)
    for _ in telo_data:
        telo_dict[_[0]].append(_[1:])



    contig_data = import_data(PAF_FILE_PATH_[0])

    original_node_count += len(contig_data)

    node_label = label_node(contig_data, telo_dict)

    repeat_label = label_repeat_node(contig_data, repeat_data)

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

    # with open("telo_preprocess_contig.txt", "wt") as f:
    #     for i in new_contig_data:
    #         for j in i:
    #             print(j, end="\t", file=f)
    #         print("", file=f)
    
    new_node_repeat_label = label_repeat_node(new_contig_data, repeat_data)
    new_node_repeat_censat_label = label_repeat_node(new_contig_data, repeat_censat_data)
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
            temp_list.append('0.'+str(new_contig_data[i][10]))
            telo_ppc_contig.append(temp_list)
            cnt+=1
    
    mainflow_dict = find_mainflow(telo_ppc_contig)
    telo_ppc_size = len(telo_ppc_contig)
    for i in range(telo_ppc_size):
        max_chr = mainflow_dict[telo_ppc_contig[i][CTG_NAM]]
        temp = [max_chr[0], max_chr[1], telo_ppc_contig[i][-1]]
        telo_ppc_contig[i] = telo_ppc_contig[i][:-1]
        telo_ppc_contig[i] += temp

    final_contig = preprocess_repeat(telo_ppc_contig)

    final_contig_repeat_label = label_repeat_node(final_contig, repeat_data)

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
        original_node_count += len(contig_data)
        node_label = label_node(contig_data, telo_dict)
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
        alt_telo_ppc_contig = []
        new_node_telo_label = label_node(new_contig_data, telo_dict)
        new_node_repeat_label = label_repeat_node(new_contig_data, repeat_data)
        new_node_repeat_censat_label = label_repeat_node(new_contig_data, repeat_censat_data)

        ref_qry_ratio = calc_ratio(new_contig_data)
        preprocess_result, \
        preprocess_type3_result, \
        preprocess_contig_type, \
        preprocess_terminal_nodes, \
        alt_len_counter = alt_preprocess_contig(new_contig_data, new_node_telo_label, ref_qry_ratio, new_node_repeat_label, telcon_set, telo_dict)
        new_contig_data = new_contig_data[0:-1]
        contig_data_size = len(new_contig_data)
        bias = cnt
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

        alt_mainflow_dict = find_mainflow(alt_telo_ppc_contig)
        alt_telo_ppc_size = len(alt_telo_ppc_contig)
        for i in range(alt_telo_ppc_size):
            max_chr = alt_mainflow_dict[alt_telo_ppc_contig[i][CTG_NAM]]
            temp = [max_chr[0], max_chr[1], alt_telo_ppc_contig[i][-1]]
            alt_telo_ppc_contig[i] = alt_telo_ppc_contig[i][:-1]
            alt_telo_ppc_contig[i] += temp
        
        alt_final_contig = preprocess_repeat(alt_telo_ppc_contig)

        alt_final_repeat_node_label = label_repeat_node(alt_final_contig, repeat_data)

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
        alt_final_ctg_typ, alt_final_preprocess_terminal_nodes, _ = alt_preprocess_contig(alt_final_contig, alt_final_telo_node_label, alt_final_ref_qry_ratio, alt_final_repeat_node_label, alt_final_telo_connect, telo_dict)
        
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

    # with open("subtelo.txt", "wt") as f:
    #     for i in subtelo_ppc_node:
    #         for j in i:
    #             print(j, end="\t", file=f)
    #         print(file=f)

    adjacency = initial_graph_build(real_final_contig, telo_bound_dict)

    telo_connected_node, telo_connected_dict, _, telo_coverage = edge_optimization(real_final_contig, adjacency, telo_bound_dict, asm2cov)

    # for i in range(total_len):
    #     if i in telo_connected_node:
    #         real_final_contig[i][CTG_TELCON] = telo_connected_dict[i]
    
    break_contig = break_double_telomere_contig(real_final_contig, telo_connected_node)

    if len(break_contig) > 0:
        final_break_contig = pass_pipeline(break_contig, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, False)
    else:
        final_break_contig = []

    if len(cen_vtg_contig) > 0:
        final_cen_vtg_contig = pass_pipeline(cen_vtg_contig, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, False)
    else:
        final_cen_vtg_contig = []

    if len(subtelo_ppc_node) > 0:
        final_subtelo_ppc_node = pass_pipeline(subtelo_ppc_node, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, True)
    else:
        final_subtelo_ppc_node = []

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

    overlap_low_split_contig = []
    if is_unitig_reduced == False:
        ctgname2overlap = get_overlap_total_score_dict(ORIGNAL_PAF_LOC_LIST)
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
                        
                        if real_final_contig[current_idx][CTG_END] - real_final_contig[prev_element_idx][CTG_STR] >= SPLIT_CTG_LEN_LIMIT:
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
        final_low_split_contig = pass_pipeline(overlap_low_split_contig, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, False)
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

    adjacency = initial_graph_build(real_final_contig, telo_bound_dict)

    telo_connected_node, telo_connected_dict, telo_connected_graph_dict, telo_coverage = edge_optimization(real_final_contig, adjacency, telo_bound_dict, asm2cov)

    rev_telo_connected_dict = defaultdict(list)

    for i in telo_connected_dict:
        rev_telo_connected_dict[telo_connected_dict[i]].append(i)
    
    with open(f"{PREFIX}/telomere_connected_list_readable.txt", "wt") as f:
        for i in rev_telo_connected_dict:
            print(i, file=f)
            for j in rev_telo_connected_dict[i]:
                print(tuple(real_final_contig[j]), file=f)
            print("", file=f)
    
    with open(f"{PREFIX}/telomere_connected_list.txt", "wt") as f:
        for i in telo_connected_graph_dict:
            for j in telo_connected_graph_dict[i]:
                print(i, tuple(j), sep = "\t", file=f)

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

def nclose_calc():
    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
    contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)

    TELO_CONNECT_NODES_INFO_PATH = PREFIX+"/telomere_connected_list.txt"

    os.makedirs(PREFIX, exist_ok=True)

    contig_data_size = len(contig_data)

    chr_corr = {}
    chr_rev_corr = {}
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
    rpt_censat_con = check_censat_contig(rpt_censat_con, PAF_FILE_PATH, ORIGNAL_PAF_LOC_LIST, contig_data)

    bnd_contig = extract_bnd_contig(contig_data)


    # Type 1, 2, 4에 대해서 
    nclose_nodes, nclose_start_compress, nclose_end_compress, vctg_dict, all_nclose_comp, nclose_coverage, nclose_compress_track = \
        extract_nclose_node(contig_data, bnd_contig, rpt_con, rpt_censat_con, repeat_censat_data, PAF_FILE_PATH, ORIGNAL_PAF_LOC_LIST, telo_set, chr_len, asm2cov)

    virtual_ordinary_contig = make_virtual_ord_ctg(contig_data, vctg_dict)
    with open(f"{PREFIX}/virtual_ordinary_contig.txt", "wt") as f:
        for i in virtual_ordinary_contig:
            for j in i:
                print(j, end = "\t", file=f)
            print("", file = f)


    nclose_type = defaultdict(list)
    for j in nclose_nodes:
        for pair in nclose_nodes[j]:
            contig_a = contig_data[pair[0]]
            contig_b = contig_data[pair[1]]
            if chr2int(contig_a[CHR_NAM]) <= chr2int(contig_b[CHR_NAM]):
                nclose_type[(contig_a[CHR_NAM], contig_b[CHR_NAM])].append(pair)
            else:
                nclose_type[(contig_b[CHR_NAM], contig_a[CHR_NAM])].append((pair[1], pair[0]))

    all_nclose_type = defaultdict(list)

    for i in all_nclose_comp:
        pairs = all_nclose_comp[i]
        for pair in pairs:
            contig_a = contig_data[pair[0]]
            contig_b = contig_data[pair[1]]
            if chr2int(contig_a[CHR_NAM]) <= chr2int(contig_b[CHR_NAM]):
                all_nclose_type[(contig_a[CHR_NAM], contig_b[CHR_NAM])].append(pair)
            else:
                all_nclose_type[(contig_b[CHR_NAM], contig_a[CHR_NAM])].append((pair[1], pair[0]))

    uncomp_node_count = 0
    with open(f"{PREFIX}/all_nclose_nodes_list.txt", "wt") as f:
        for i in all_nclose_type:
            print(f"{i[0]}, {i[1]}, {len(all_nclose_type[i])}", file=f)
            for pair in all_nclose_type[i]:
                uncomp_node_count += 2
                contig_a = contig_data[pair[0]]
                contig_b = contig_data[pair[1]]

                is_for = pair[0] < pair[1]
                list_a = [contig_a[CTG_NAM], get_corr_dir(is_for, contig_a[CTG_DIR]), contig_a[CHR_STR], contig_a[CHR_END]]
                list_b = [contig_b[CTG_NAM], get_corr_dir(is_for, contig_b[CTG_DIR]), contig_b[CHR_STR], contig_b[CHR_END]]
                if contig_a[CTG_NAM] in rpt_con:
                    print(list_a, list_b, "all_repeat", file=f)
                else:
                    print(list_a, list_b, file=f)
            print("", file=f)

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

    transloc_nclose_pair_count = 0
    with open(f"{PREFIX}/compressed_nclose_nodes_list.txt", "wt") as f:
        for i in nclose_type:
            print(f"{i[0]}, {i[1]}, {len(nclose_type[i])}", file=f)
            if i[0] != i[1]:
                transloc_nclose_pair_count += len(nclose_type[i])

            st_flag = False
            ed_flag = False
            if (('=', i[0]), ('=', i[1])) in nclose_start_compress:
                st_flag = True
            if (('=', i[0]), ('=', i[1])) in nclose_end_compress:
                ed_flag = True
            for pair in nclose_type[i]:
                contig_a = contig_data[pair[0]]
                contig_b = contig_data[pair[1]]

                if st_flag:
                    if contig_a[CTG_NAM] in nclose_start_compress[(('=', i[0]), ('=', i[1]))]:
                        pass
                
                is_for = pair[0] < pair[1]
                list_a = [contig_a[CTG_NAM], get_corr_dir(is_for, contig_a[CTG_DIR]), contig_a[CHR_STR], contig_a[CHR_END]]
                list_b = [contig_b[CTG_NAM], get_corr_dir(is_for, contig_b[CTG_DIR]), contig_b[CHR_STR], contig_b[CHR_END]]
                if contig_a[CTG_NAM] in rpt_con:
                    print(list_a, list_b, "all_repeat", file=f)
                else:
                    print(list_a, list_b, file=f)
            print("", file=f)

    nclose_cord_list = []
    total_nclose_cord_list_contig_name = []
    nclose_idx_corr = []
    total_dir_data = {}
    transloc_k_set = set()

    trusted_nclose_count = 0
    for j in nclose_nodes:
        for i in nclose_nodes[j]:
            tf_set = set()

            s = i[0]
            e = i[1]
            curr_nclose_cord_list = []
            dir_data = defaultdict(dict)
            nclose_cord_list_contig_name = []
            nclose_maxcover_s = contig_data[s][CHR_END] if contig_data[s][CTG_DIR] == '+' else contig_data[s][CHR_STR]
            nclose_maxcover_e = contig_data[e][CHR_STR] if contig_data[e][CTG_DIR] == '+' else contig_data[e][CHR_END]
            start_chr = contig_data[s][CHR_NAM]
            end_chr = contig_data[e][CHR_NAM]
            temp_list = [start_chr, nclose_maxcover_s, contig_data[s][CTG_DIR],
                        end_chr, nclose_maxcover_e, contig_data[e][CTG_DIR],
                        trusted_nclose_count]

            template_dir = (contig_data[s][CTG_DIR], contig_data[e][CTG_DIR])
            tf_set.add((True, False))

            dir_data[trusted_nclose_count][(True, False)] = (contig_data[s][CTG_NAM], s, e, contig_data[s][CTG_TYP])
            curr_nclose_cord_list.append(temp_list)
            nclose_cord_list_contig_name.append(curr_nclose_cord_list[-1]+[contig_data[s][CTG_NAM]])

            for compressed_contig in nclose_compress_track[tuple(i)]:

                is_forward = None
                if contig_data[compressed_contig[0]][CHR_NAM] == start_chr and contig_data[compressed_contig[1]][CHR_NAM] == end_chr:
                    compress_s = compressed_contig[0]
                    compress_e = compressed_contig[1]
                    is_forward = True
                elif contig_data[compressed_contig[0]][CHR_NAM] == end_chr and contig_data[compressed_contig[1]][CHR_NAM] == start_chr:
                    compress_s = compressed_contig[1]
                    compress_e = compressed_contig[0]
                    is_forward = False
                else:
                    assert(False)

                nclose_corr_dir = (get_corr_dir(is_forward, contig_data[compress_s][CTG_DIR]),
                                   get_corr_dir(is_forward, contig_data[compress_e][CTG_DIR]))
                
                nclose_maxcover_s = contig_data[compress_s][CHR_END] if nclose_corr_dir[0] == '+' else contig_data[compress_s][CHR_STR]
                nclose_maxcover_e = contig_data[compress_e][CHR_STR] if nclose_corr_dir[1] == '+' else contig_data[compress_e][CHR_END]

                temp_list = [start_chr, nclose_maxcover_s, nclose_corr_dir[0],
                             end_chr, nclose_maxcover_e, nclose_corr_dir[1],
                             trusted_nclose_count]

                reference_dir = (True if nclose_corr_dir[0] == template_dir[0] else False, 
                                 False if nclose_corr_dir[1] == template_dir[1] else True)

                if reference_dir not in dir_data[trusted_nclose_count]:
                    dir_data[trusted_nclose_count][reference_dir] = (contig_data[compress_s][CTG_NAM], compress_s, compress_e, contig_data[compress_s][CTG_TYP])
                
                tf_set.add(reference_dir)
                curr_nclose_cord_list.append(temp_list)
                nclose_cord_list_contig_name.append(curr_nclose_cord_list[-1]+[contig_data[compress_e][CTG_NAM]])
            
            if len(tf_set) >= 1:
                if tf_set == {(True, False), (False, True)} and template_dir[0] == template_dir[1]:
                    transloc_k_set.add(trusted_nclose_count)
                
                assert(len(dir_data.values()) == 1)
                dd = list(dir_data.values())[0]

                assert(len(dd) >= 1)
                assert((True, False) in dd)

                nclose_cord_list += curr_nclose_cord_list
                total_nclose_cord_list_contig_name += nclose_cord_list_contig_name
                nclose_idx_corr.append(i)
                total_dir_data.update(dir_data)
                trusted_nclose_count += 1

    with open(f"{PREFIX}/03_anal_bam_input.pkl", "wb") as f:
        pkl.dump((nclose_cord_list, nclose_idx_corr, total_nclose_cord_list_contig_name,
                  total_dir_data, transloc_k_set, nclose_nodes), f)

    thread_lim = min(16, THREAD)
    PROGRESS = ['--progress'] if args.progress else []
    subprocess.run(['python', "-X", f"juliacall-threads={thread_lim}", "-X", "juliacall-handle-signals=yes", 
                    os.path.join(os.path.dirname(os.path.abspath(__file__)), '03_Anal_bam.py'),
                    PREFIX, read_bam_loc, CENSAT_PATH, CHROMOSOME_INFO_FILE_PATH, main_stat_loc] + PROGRESS)

    with open(f"{PREFIX}/03_anal_bam_output.pkl", "rb") as f:
        nclose_nodes, task_cnt = pkl.load(f)
    
    logging.info(f"Deleted NClose node with no coverage : {task_cnt[3]}")
    logging.info(f"Added directed NClose node with coverage : {task_cnt[1]}")
    logging.info(f"Translocation NClose node with coverage : {task_cnt[2]}")

    nclose_node_count = 0
    with open(f"{PREFIX}/nclose_nodes_index.txt", "wt") as f: 
        for j in nclose_nodes:
            for i in nclose_nodes[j]:
                nclose_node_count += 2
                print(j, i[0], i[1], contig_data[i[0]][CTG_TYP], file=f)

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
            
parser = argparse.ArgumentParser(description="Find breakend contigs with contig data and map data")

# 위치 인자 정의
parser.add_argument("paf_file_path", 
                    help="Path to the original PAF file.")
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
                    help="Path to an alternative PAF file (optional).")
parser.add_argument("--orignal_paf_loc", nargs='+',
                    help="Orignal paf location to detect location (primary, alternative paf location)")
parser.add_argument("-t", "--thread", 
                    help="Number of thread", type=int)
parser.add_argument("--progress", 
                    help="Show progress bar", action='store_true')
parser.add_argument("--verbose", 
                    help="Enable index, paf output (Could be slow at HDD)", action='store_true')

args = parser.parse_args()

# t = "02_Build_Breakend_Graph_Limited.py /home/hyunwoo/ACCtools-pipeline/90_skype_run/PC-3/20_alignasm/PC-3.ctg.aln.paf public_data/chm13v2.0.fa.fai public_data/chm13v2.0_telomere.bed public_data/chm13v2.0_repeat.m.bed public_data/chm13v2.0_censat_v2.1.m.bed /home/hyunwoo/ACCtools-pipeline/90_skype_run/PC-3/01_depth/PC-3.win.stat.gz 30_skype_pipe/PC-3_13_21_51 /home/hyunwoo/ACCtools-pipeline/90_skype_run/PC-3/01_depth/PC-3.bam --alt /home/hyunwoo/ACCtools-pipeline/90_skype_run/PC-3/20_alignasm/PC-3.utg.aln.paf --orignal_paf_loc /home/hyunwoo/ACCtools-pipeline/90_skype_run/PC-3/20_alignasm/PC-3.ctg.paf /home/hyunwoo/ACCtools-pipeline/90_skype_run/PC-3/20_alignasm/PC-3.utg.paf -t 128"
# args = parser.parse_args(t.split()[1:])

PREFIX = args.prefix

os.makedirs(PREFIX, exist_ok=True)

PAF_FILE_PATH = []
if args.alt is None:
    PAF_FILE_PATH = [args.paf_file_path]
else:
    PAF_FILE_PATH = [args.paf_file_path, args.alt]

PREPROCESSED_PAF_FILE_PATH = (PAF_FILE_PATH[0] +'.ppc.paf')
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
REPEAT_INFO_FILE_PATH = args.repeat_bed_path
CENSAT_PATH = args.censat_bed_path
ORIGNAL_PAF_LOC_LIST_ = args.orignal_paf_loc
main_stat_loc = args.main_stat_path
read_bam_loc = args.read_bam_loc
PRINT_IDX_FILE = args.verbose
THREAD=args.thread

ORIGNAL_PAF_LOC_LIST = ORIGNAL_PAF_LOC_LIST_


assert(len(PAF_FILE_PATH) == len(ORIGNAL_PAF_LOC_LIST))

gfa_file_path = []
asm2cov = Counter()
# with open(f'{PREFIX}/asm2cov.pkl', 'rb') as f:
#     gfa_file_path, asm2cov = pkl.load(f)

ori_ctg_name_data = get_ori_ctg_name_data(PAF_FILE_PATH)

is_unitig_reduced = False
telo_coverage = contig_preprocessing_00(PAF_FILE_PATH)
globals().update(nclose_calc())

if nclose_node_count > FAIL_NCLOSE_COUNT:
    logging.info("NClose node count is too high.")
    if len(PAF_FILE_PATH) == 1:
        logging.info("No method to reduce nclose node count.")
        sys.exit(1)
    else:
        is_unitig_reduced = True
        logging.info("Retrying with the primary PAF file.")

        PAF_FILE_PATH = [args.paf_file_path, args.paf_file_path]
        ORIGNAL_PAF_LOC_LIST = [ORIGNAL_PAF_LOC_LIST_[0], ORIGNAL_PAF_LOC_LIST_[0]]
        ori_ctg_name_data = get_ori_ctg_name_data(PAF_FILE_PATH)

        telo_coverage = contig_preprocessing_00(PAF_FILE_PATH)
        globals().update(nclose_calc())

        if nclose_node_count > FAIL_NCLOSE_COUNT:
            logging.info("No method to reduce nclose node count.")
            sys.exit(1)

bnd_graph_adjacency = initialize_bnd_graph(contig_data, nclose_nodes, telo_contig)

with open(f"{PREFIX}/bnd_only_graph.txt", "wt") as f: 
    for i in bnd_graph_adjacency:
        print(i, bnd_graph_adjacency[i], file=f)


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
                if contig_s[CHR_NAM] != contig_e[CHR_NAM]:
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

def run_graph(data, CHR_CHANGE_LIMIT, DIR_CHANGE_LIMIT):
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
                for i in range(min(PATH_MAJOR_COMPONENT, len(ack))):
                    if ack[i][0] in (src[0][:-1], tar[:-1]):
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
                        
                        path_list.append((path, ack))
                        path_di_list.append(ii)

                        if cnt >= PAT_PATH_LIMIT:
                            return ((src[0], tar), cnt, path_list, path_di_list)

    return ((src[0], tar), cnt, path_list, path_di_list)

def init_worker(shared_graph):
    global G
    G = shared_graph

def run_graph_pipeline():
    tar_ind_list = []
    for i in range(contig_data_size, contig_data_size + 2*CHROMOSOME_COUNT):
        for j in range(i, contig_data_size + 2*CHROMOSOME_COUNT):
            tar_ind_list.append((i, j))

    limit_combinations = [
        (c, DIR_CHANGE_LIMIT_ABS_MAX)
        for c in range(CHR_CHANGE_LIMIT_ABS_MAX, 0, -1)
    ]
    limit_combinations.append((1, 0))

    if nclose_node_count > HARD_NCLOSE_COUNT:
        idx = limit_combinations.index((CHR_CHANGE_LIMIT_HARD_START, DIR_CHANGE_LIMIT_ABS_MAX))
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
            if last_success is not None:
                break

            idx += 1
    return last_success

last_success = run_graph_pipeline()
if not last_success:
    logging.info('Breakend graph is too divergent.')
    if is_unitig_reduced or len(PAF_FILE_PATH) == 1:
        logging.info('Breakend path failed.')
        sys.exit(1)
    else:
        is_unitig_reduced = True
        logging.info("Retrying with the primary PAF file.")

        PAF_FILE_PATH = [args.paf_file_path, args.paf_file_path]
        ORIGNAL_PAF_LOC_LIST = [ORIGNAL_PAF_LOC_LIST_[0], ORIGNAL_PAF_LOC_LIST_[0]]
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

logging.info(f'Final success settings: {(CHR_CHANGE_LIMIT_PREFIX, DIR_CHANGE_LIMIT_PREFIX)}')

cancer_prefix = os.path.basename(PREPROCESSED_PAF_FILE_PATH).split('.')[0]

with open(f'{PREFIX}/path_data.pkl', 'wb') as f:
    pkl.dump(path_list_dict, f)

with open(f'{PREFIX}/path_di_data.pkl', 'wb') as f:
    pkl.dump(path_di_list_dict, f)

with open(f'{PREFIX}/paf_file_path.pkl', 'wb') as f:
    pkl.dump(PAF_FILE_PATH, f)

with open(f'{PREFIX}/report.txt', 'a') as f:
    cnt = sum(i[1] for i in cnt_list)
    print(cancer_prefix, file=f)
    print(cnt, file=f)
    logging.info(f"Total path count : {cnt}")
    for (st, nd), c in sorted(cnt_list, key=lambda x:-x[1]):
        if c > 0:
            print(st, nd, c, file=f)
