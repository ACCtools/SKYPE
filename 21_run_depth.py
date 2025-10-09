import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *

import re
import ast
import glob
import copy
import logging
import argparse
import subprocess

import numpy as np
import pickle as pkl
import networkx as nx

from tqdm import tqdm
from collections import defaultdict
from itertools import pairwise, groupby
from concurrent.futures import ProcessPoolExecutor, as_completed

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("21_run_depth_eff start")

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
CTG_MAINFLOWDIR = 19
CTG_MAINFLOWCHR = 20
CTG_GLOBALIDX = 21


BND_TYPE = 0
CTG_IN_TYPE = 1
TEL_TYPE = 2

DIR_FOR = 1
DIR_BAK = 0
DIR_IN = 3
DIR_OUT = 2
INF = 1000000000
BUFFER = 10000000


PANDEPTH_RETRY = 3
DEPTH_WINDOW = 100 * K
DEPTH_THREAD = 1

NCLOSE_COMPRESS_LIMIT = 50*K
VIRTUAL_CONTIG_PREFIX = "virtual_contig"

def import_data(file_path : list) -> list :
    contig_data = []
    int_induce_idx = [1, 2, 3, 6, 7, 8, 9, 10, 11]
    idx = 0
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            curr_contig = curr_contig.rstrip()
            temp_list = curr_contig.split("\t")
            for i in int_induce_idx:
                temp_list[i] = int(temp_list[i])
            contig_data.append(temp_list)
    return contig_data

def import_data2(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.split("\t")
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END,
                          CHR_LEN, CHR_STR, CHR_END,
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

    path = path_list_dict[key][cnt][0]
    return [tuple(list(n)[:2]) for n in path]

def distance_checker(node_a : tuple, node_b : tuple) -> int :
    # check if overlap
    if max(node_a[CHR_STR], node_b[CHR_STR]) < min(node_a[CHR_END], node_b[CHR_END]):
        return 0   
    else:
        # return distance of two ranges
        return min(abs(node_b[CHR_STR] - node_a[CHR_END]), abs(node_b[CHR_END] - node_a[CHR_STR]))

def telo_distance_checker(node: tuple, telo: tuple) -> int :
    return min(abs(telo[CHR_STR] - node[CHR_END]), abs(telo[CHR_END] - node[CHR_STR]))
    
def inclusive_checker(contig_node : tuple, telomere_node : tuple) -> bool :
    if int(telomere_node[CHR_STR]) <= int(contig_node[CHR_STR]) and int(contig_node[CHR_END]) <= int(telomere_node[CHR_END]):
        return True
    else:
        return False

def import_final_paf(final_paf_path : str) -> list:
    paf_data = []
    with open(final_paf_path, 'r') as f:
        for l in f:
            paf_data.append(l[:-1].split('\t'))

    return paf_data

def find_index(fill_index_path, final_paf_data, tar_index):
    index_i = 0
    for i, paf_line in enumerate(final_paf_data):
        if not paf_line[CTG_NAM].startswith(VIRTUAL_CONTIG_PREFIX):
            if fill_index_path[index_i] == tar_index:
                return i
            index_i += 1
    return None

def get_key_from_index_folder(key):
    index_key_data = []
    
    n = len(path_list_dict[key])
    for i in range(n):
        index_key_data.append(get_key_from_index_file(f"{PREFIX}/00_raw/{key}/{i + 1}.index.txt"))

    return key, index_key_data

def get_key_from_index_file(index_file_path):
    bnd_index_data = import_index_path(index_file_path)[1:-1]
    bnd_ctg_group = count_groups([contig_data[i[1]][CTG_NAM] for i in bnd_index_data])

    key_list = []
    is_contig_same = True
    for i, (ci, nci) in enumerate(pairwise(bnd_index_data)):  
        (ctg_dir, ctg_ind), (next_ctg_dir, next_ctg_ind) = ci, nci

        if i == 0 or i == len(bnd_index_data) - 2:
            if i == 0 and bnd_ctg_group[0] == 3:
                is_contig_same = False
            elif i == len(bnd_index_data) - 2 and bnd_ctg_group[-1] == 3:
                is_contig_same = False
            else:
                is_contig_same = True
        else:
            is_contig_same = True

        if i == 0 and bnd_ctg_group[0] == 2:
            key_list.append((TEL_TYPE, ctg_ind))

        if contig_data[ctg_ind][CTG_NAM] == contig_data[next_ctg_ind][CTG_NAM] and is_contig_same:
            key_list.append((CTG_IN_TYPE, (ci, nci)))
        else:
            key_list.append((BND_TYPE, (ci, nci)))
            
        if i == len(bnd_index_data) - 2 and bnd_ctg_group[-1] == 2:
            key_list.append((TEL_TYPE, ctg_ind))
    
    return index_file_path, key_list

def node_dir_transform(node):
    transform_node = list(copy.deepcopy(node))
    transform_node[CTG_STR] = node[CTG_LEN] - node[CTG_END]
    transform_node[CTG_END] = node[CTG_LEN] - node[CTG_STR]
    transform_node[CTG_DIR] = '+' if node[CTG_DIR] == '+' else '-'
    return transform_node

def form_virtual_contig(chr_len, chr_str, chr_end, chr_name, vcnt):
    virtual_contig = []
    N = chr_end - chr_str
    virtual_contig = [f"{VIRTUAL_CONTIG_PREFIX}_{vcnt}", N, 0, N]
    virtual_contig += ['+', chr_name, chr_len, chr_str, chr_end]
    virtual_contig += [N, N, 0]
    virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
    return virtual_contig

def form_normal_contig(node):
    paf_idx = int(node[CTG_GLOBALIDX].split('.')[0])
    if paf_idx < 2:
        node_idx = int(node[CTG_GLOBALIDX].split('.')[1])
        new_contig = paf_file[paf_idx][node_idx]
    else:
        N = node[CHR_END] - node[CHR_STR]
        new_contig = list(node[:9]) + [N, N, 0, 'tp:A:P', 'cs:Z:'+f":{N}"]
    return new_contig

def form_adjusted_contig(node, paf_data):
    new_node = list(copy.deepcopy(node))
    
    new_node[CHR_STR] = paf_data[0]
    new_node[CHR_END] = paf_data[1]
    new_node[CTG_STR] = paf_data[2]
    new_node[CTG_END] = paf_data[3]
    new_node[9] = paf_data[4]
    new_node[10] = paf_data[5]
    new_node[-1] = 'cs:Z:' + paf_data[6]
    return new_node

def cs_to_cigar(cs_tag: str) -> str:
   """
   cs:Z: 태그의 값(접두어 "cs:Z:" 제외)을 받아서 표준 CIGAR 문자열로 변환합니다.

   변환 규칙:
     - ":<number>" : 매칭(=, M) 연산. 이전 연산이 'M'이면 길이를 누적하고,
                    그렇지 않으면 이전 연산을 플러시한 후 'M'으로 시작합니다.
     - "+<seq>"    : 삽입(Insertion). 이전에 누적된 연산이 있으면 먼저 플러시하고,
                    삽입 길이는 해당 문자열의 길이로 계산합니다.
     - "-<seq>"    : 결실(Deletion). 위와 유사하게 처리합니다.
     - "*<base1><base2>": 치환(Mismatch). 연속되는 치환은 길이를 누적합니다.
   """
   pattern = re.compile(r"(:\d+|\*[a-z]{2}|[+\-][A-Za-z]+)")
   cigar = ""
   last_op = 'M'
   last_len = 0

   for match in pattern.finditer(cs_tag):
       part = match.group(0)
       op = part[0]
       if op == ":":
           # 예: ":10" → 매칭 10개
           length = int(part[1:])
           if last_op == 'M':
               last_len += length
           else:
               if last_len > 0:
                   cigar += f"{last_len}{last_op}"
               last_op = 'M'
               last_len = length
       elif op == "-":
           # 예: "-acgt" → deletion, 길이는 4
           length = len(part[1:])
           if last_len > 0:
               cigar += f"{last_len}{last_op}"
           cigar += f"{length}D"
           last_len = 0
           last_op = 'M'
       elif op == "+":
           # 예: "+ac" → insertion, 길이는 2
           length = len(part[1:])
           if last_len > 0:
               cigar += f"{last_len}{last_op}"
           cigar += f"{length}I"
           last_len = 0
           last_op = 'M'
       elif op == "*":
           # 예: "*at" : mismatch (치환)
           if last_op == 'X':
               last_len += 1
           else:
               if last_len > 0:
                   cigar += f"{last_len}{last_op}"
               last_op = 'X'
               last_len = 1
       else:
           # 정의되지 않은 연산은 무시
           continue

   if last_len > 0:
       cigar += f"{last_len}{last_op}"

   return cigar

def compute_aln_stats(cs_string):
    """
    cs_string을 파싱하여 매치 길이(mat_len)와 정렬 길이(aln_len)를 계산합니다.
    
    cs_string 내 연산자 종류:
      - :N    => N개의 매치 (매치와 정렬 길이에 N 추가)
      - *XY   => 치환 (정렬 길이에 1 추가)
      - +SEQ  => 삽입 (여기서는 기준(reference)에는 해당하지 않으므로 무시)
      - -SEQ  => 결실 (정렬 길이에 SEQ의 길이 추가)
    """
    mat_len = 0
    aln_len = 0

    # 정규표현식 패턴: 매치, 치환, 삽입, 결실을 그룹으로 구분
    pattern = r"(?::(\d+))|(?:\*([A-Za-z]{2}))|(?:\+([A-Za-z]+))|(?:\-([A-Za-z]+))"
    
    for match in re.finditer(pattern, cs_string):
        if match.group(1):  # ":" 매치: 숫자
            n = int(match.group(1))
            mat_len += n
            aln_len += n
        elif match.group(2):  # "*" 치환: 두 문자 (예: *ac)
            aln_len += 1
        elif match.group(3):  # "+" 삽입: query에만 있으므로 기준 정렬 길이에는 포함하지 않음
            aln_len += len(match.group(3))
        elif match.group(4):  # "-" 결실: 삭제된 서열 길이를 기준 정렬 길이에 더함
            aln_len += len(match.group(4))
            
    return mat_len, aln_len

def parse_cs(cs):
    """
    cs 문자열을 파싱하여 연산 리스트를 반환.
    각 연산은 (op, value) 튜플로 나타냅니다.
      - 매칭: (":", length)        ex: (":", 10)
      - 삽입: ("+", seq)            ex: ("+", "ACGT")
      - 결실: ("-", seq)             ex: ("-", "ACG")
      - 치환: ("*", (ref, qry))      ex: ("*", ("a", "t"))
    """
    pattern = re.compile(r"(:\d+|\*[a-z]{2}|[+\-][A-Za-z]+)")
    ops = []
    for token in pattern.findall(cs):
        if token.startswith(":"):
            length = int(token[1:])
            ops.append( (":", length) )
        elif token.startswith("+"):
            seq = token[1:]
            ops.append( ("+", seq) )
        elif token.startswith("-"):
            seq = token[1:]
            ops.append( ("-", seq) )
        elif token.startswith("*"):
            # token like "*at"
            ops.append( ("*", (token[1], token[2])) )
    return ops

def construct_cs(ops):
    """op 리스트를 다시 cs 문자열로 조합"""
    cs = ""
    for op, val in ops:
        if op == ":":
            cs += f":{val}"
        elif op in ("+", "-"):
            cs += f"{op}{val}"
        elif op == "*":
            cs += f"*{val[0]}{val[1]}"
    return cs

def trim_cs_right_extended(cs, keep_ref):
    """
    cs 문자열의 왼쪽 부분에서, 참조(consumed)가 keep_ref 만큼 소비될 때까지 자르고,
    그때 소비된 query의 길이도 함께 반환합니다.
    반환값: (새로운 cs 문자열, query_consumed)
    """
    ops = parse_cs(cs)
    new_ops = []
    ref_consumed = 0
    query_consumed = 0
    for op, val in ops:
        if op == ":":
            if ref_consumed + val < keep_ref:
                new_ops.append( (":", val) )
                ref_consumed += val
                query_consumed += val
            else:
                remaining = keep_ref - ref_consumed
                if remaining > 0:
                    new_ops.append( (":", remaining) )
                    ref_consumed += remaining
                    query_consumed += remaining
                break
        elif op == "-":
            n = len(val)
            if ref_consumed + n < keep_ref:
                new_ops.append( ("-", val) )
                ref_consumed += n
                # deletion does not consume query bases
            else:
                remaining = keep_ref - ref_consumed
                if remaining > 0:
                    new_ops.append( ("-", val[:remaining]) )
                    ref_consumed += remaining
                break
        elif op == "*":
            # substitution consumes 1 ref and 1 query
            if ref_consumed + 1 <= keep_ref:
                new_ops.append( ("*", val) )
                ref_consumed += 1
                query_consumed += 1
            else:
                break
        elif op == "+":
            # insertion: does not consume ref, only query
            new_ops.append( ("+", val) )
            query_consumed += len(val)
    return construct_cs(new_ops), query_consumed

def trim_cs_left_extended(cs, trim_ref):
    """
    cs 문자열의 왼쪽에서, 참조(consumed)가 trim_ref 만큼 소비되는 부분을 제거하고,
    제거된 부분에 해당하는 query의 길이(query_removed)도 함께 반환합니다.
    반환값: (새로운 cs 문자열, query_removed)
    """
    ops = parse_cs(cs)
    new_ops = []
    ref_skipped = 0
    query_skipped = 0
    i = 0
    while i < len(ops) and ref_skipped < trim_ref:
        op, val = ops[i]
        if op == ":":
            if ref_skipped + val <= trim_ref:
                ref_skipped += val
                query_skipped += val
                i += 1
            else:
                remaining = val - (trim_ref - ref_skipped)
                query_skipped += (trim_ref - ref_skipped)
                ref_skipped = trim_ref
                new_ops.append( (":", remaining) )
                i += 1
                break
        elif op == "-":
            n = len(val)
            if ref_skipped + n <= trim_ref:
                ref_skipped += n
                # deletion does not add to query consumption
                i += 1
            else:
                to_skip = trim_ref - ref_skipped
                remaining = n - to_skip
                ref_skipped = trim_ref
                new_ops.append( ("-", val[to_skip:]) )
                i += 1
                break
        elif op == "*":
            # substitution consumes 1 ref and 1 query
            if ref_skipped + 1 <= trim_ref:
                ref_skipped += 1
                query_skipped += 1
                i += 1
            else:
                i += 1
                break
        elif op == "+":
            # insertion: does not consume ref → skip it entirely and add its query length
            query_skipped += len(val)
            i += 1
    new_ops.extend(ops[i:])
    return construct_cs(new_ops), query_skipped

def adjust_paf_overlap(paf1_data, paf2_data):
    """
    두 개의 PAF 데이터를 입력받아 (각각  
    [ref_st, ref_nd, qry_st, qry_nd, mat_len, aln_len, cs_string]),
    두 정렬의 참조 영역이 반드시 겹친다고 가정합니다.
    겹치지 않으면 에러를 발생시킵니다.
    
    겹치는 영역의 중간(new_boundary)을 기준으로, paf1은 좌측 부분(참조: ref_st ~ new_boundary),
    paf2는 우측 부분(참조: new_boundary ~ ref_nd)를 남기도록 잘라냅니다.
    
    또한, aln_dir은 무조건 '+'라고 가정하여, 각 정렬의 query 좌표도 적절히 조정합니다.
    
    반환:
       new_paf1, new_paf2 (각각 [new_ref_st, new_ref_nd, new_qry_st, new_qry_nd, new_mat_len, new_aln_len, new_cs_string])
    """
    # unpack input data
    ref_st1, ref_nd1, qry_st1, qry_nd1, mat_len1, aln_len1, cs1 = paf1_data
    ref_st2, ref_nd2, qry_st2, qry_nd2, mat_len2, aln_len2, cs2 = paf2_data

    # 두 정렬의 참조 구간이 겹치는지 확인
    if ref_nd1 <= ref_st2 or ref_nd2 <= ref_st1:
        raise Exception("두 PAF의 reference 구간이 겹치지 않습니다.")
    
    # 겹치는 영역의 시작과 끝
    overlap_start = max(ref_st1, ref_st2)
    overlap_end   = min(ref_nd1, ref_nd2)
    # 중간점: 이 지점을 기준으로 나눔
    new_boundary = (overlap_start + overlap_end) // 2

    # paf1: 좌측 정렬
    # 유지할 참조 길이 = new_boundary - ref_st1
    keep_ref1 = new_boundary - ref_st1
    new_cs1, query_consumed1 = trim_cs_right_extended(cs1, keep_ref1)
    new_mat_len1, new_aln_len1 = compute_aln_stats(new_cs1)
    # query는 aln_dir이 '+'이므로, 좌측 부분은 qry_st1 그대로 시작하고,
    # new query end는 qry_st1 + query_consumed1
    new_qry1 = qry_st1 + query_consumed1
    new_paf1 = [ref_st1, new_boundary, qry_st1, new_qry1, new_mat_len1, new_aln_len1, new_cs1]

    # paf2: 우측 정렬
    # 제거할 참조 길이 = new_boundary - ref_st2 (왼쪽 부분 제거)
    trim_ref2 = new_boundary - ref_st2
    new_cs2, query_trimmed2 = trim_cs_left_extended(cs2, trim_ref2)
    new_mat_len2, new_aln_len2 = compute_aln_stats(new_cs2)
    # 우측 정렬은 query의 시작을 조정: new query start = qry_st2 + query_trimmed2, 끝은 그대로
    new_qry2 = qry_st2 + query_trimmed2
    new_paf2 = [new_boundary, ref_nd2, new_qry2, qry_nd2, new_mat_len2, new_aln_len2, new_cs2]

    return new_paf1, new_paf2

def process_raw_contig_list(full_connected_path, key_cnt):
    with open(f'{output_folder}/{key_cnt}.paf', 'w') as f:
        init_contig = form_normal_contig(contig_data[full_connected_path[0][1]])
        path_contig = [init_contig]
        full_connected_path_len = len(full_connected_path)
        vcnt = 0
        for i in range(1, full_connected_path_len):
            curr_contig = path_contig[-1]
            next_contig = contig_data[full_connected_path[i][1]]
            if curr_contig[CHR_NAM] == next_contig[CHR_NAM] and curr_contig[CTG_NAM] != next_contig[CTG_NAM] \
            or (curr_contig[CHR_NAM] == next_contig[CHR_NAM] and (full_connected_path[i-1][0] in (2, 3) or full_connected_path[i][0] in (2, 3))):
                dist = distance_checker(curr_contig, next_contig)
                if dist > 0:
                    if curr_contig[CHR_END] < next_contig[CHR_STR]:
                        if curr_contig[CTG_NAM] != next_contig[CTG_NAM]:
                            vcnt += 1
                            new_contig = form_virtual_contig(curr_contig[CHR_LEN], curr_contig[CHR_END], next_contig[CHR_STR], curr_contig[CHR_NAM], vcnt)
                            new_next_contig = form_normal_contig(next_contig)
                            path_contig.append(new_contig)
                            path_contig.append(new_next_contig)
                        else:
                            new_next_contig = form_normal_contig(next_contig)
                            path_contig.append(new_next_contig)
                    else:
                        if curr_contig[CTG_NAM] != next_contig[CTG_NAM]:
                            vcnt += 1
                            new_contig = form_virtual_contig(curr_contig[CHR_LEN], next_contig[CHR_END], curr_contig[CHR_STR], curr_contig[CHR_NAM], vcnt)
                            new_next_contig = form_normal_contig(next_contig)
                            path_contig.append(new_contig)
                            path_contig.append(new_next_contig)
                        else:
                            new_next_contig = form_normal_contig(next_contig)
                            path_contig.append(new_next_contig)
                elif curr_contig[CHR_END] == next_contig[CHR_STR] or curr_contig[CHR_STR] == next_contig[CHR_END]:
                    new_next_contig = form_normal_contig(next_contig)
                    path_contig.append(new_next_contig)
                else:
                    if curr_contig[CTG_DIR] == '-':
                        curr_contig_rc = node_dir_transform(curr_contig)
                    else:
                        curr_contig_rc = copy.deepcopy(curr_contig)
                    if next_contig[CTG_DIR] == '-':
                        next_contig_rc = node_dir_transform(next_contig)
                    else:
                        next_contig_rc = copy.deepcopy(next_contig)
                    curr_contig_data = [curr_contig_rc[CHR_STR], curr_contig_rc[CHR_END],
                                        curr_contig_rc[CTG_STR], curr_contig_rc[CTG_END],
                                        curr_contig[9], curr_contig[10],
                                        curr_contig[-1][5:]]

                    next_contig_origin = form_normal_contig(next_contig)
                    
                    next_contig_data = [next_contig_rc[CHR_STR], next_contig_rc[CHR_END],
                                        next_contig_rc[CTG_STR], next_contig_rc[CTG_END],
                                        next_contig_origin[9], next_contig_origin[10],
                                        next_contig_origin[-1][5:]]
                    # try:
                    #     assert(not (inclusive_checker_pair(tuple(curr_contig_data[0:2]), tuple(next_contig_data[0:2])) \
                    #            or inclusive_checker_pair(tuple(next_contig_data[0:2]), tuple(curr_contig_data[0:2]))))
                    # except:
                    #     print(tuple(curr_contig_data[0:2]), tuple(next_contig_data[0:2]))
                    #     print(index_file_path)
                    if full_connected_path[i-1][0] in {2, 3}:
                        is_index_inc = (next_contig[CTG_DIR]=='+' and full_connected_path[i][0]==1) or \
                                    (next_contig[CTG_DIR]=='-' and full_connected_path[i][0]==0)
                    elif full_connected_path[i][0] in {2, 3}:
                        is_index_inc = (curr_contig[CTG_DIR]=='+' and full_connected_path[i-1][0]==1) or \
                                (curr_contig[CTG_DIR]=='-' and full_connected_path[i-1][0]==0)
                    else:
                        assert(curr_contig[CTG_NAM] != next_contig[CTG_NAM])
                        is_index_inc = (curr_contig[CTG_DIR]=='+' and full_connected_path[i-1][0]==1) or \
                                (curr_contig[CTG_DIR]=='-' and full_connected_path[i-1][0]==0)

                    if is_index_inc: # 1324 2413

                        curr_paf, next_paf = adjust_paf_overlap(curr_contig_data, next_contig_data)
                        
                        curr_result = form_adjusted_contig(curr_contig, curr_paf)
                        if curr_contig[CTG_DIR]=='-':
                            curr_result = node_dir_transform(curr_result)
                        next_result = form_adjusted_contig(next_contig_origin, next_paf)
                        if next_contig[CTG_DIR]=='-':
                            next_result = node_dir_transform(next_result)
                        path_contig[-1] = curr_result
                        path_contig.append(next_result)
                    else:
                        next_paf, curr_paf = adjust_paf_overlap(next_contig_data, curr_contig_data)
                        
                        curr_result = form_adjusted_contig(curr_contig, curr_paf)
                        if curr_contig[CTG_DIR]=='-':
                            curr_result = node_dir_transform(curr_result)
                        next_result = form_adjusted_contig(next_contig_origin, next_paf)
                        if next_contig[CTG_DIR]=='-':
                            next_result = node_dir_transform(next_result)
                        path_contig[-1] = curr_result
                        path_contig.append(next_result)

                    
            else:
                next_contig_origin = form_normal_contig(next_contig)
                path_contig.append(next_contig_origin)
        cnt = 0
        for i in path_contig:
            cnt+=1
            cigar_str = 'cg:Z:' + cs_to_cigar(i[-1][5:])
            print("\t".join(list(map(str, i + [cigar_str]))), file=f)

def process_raw_contig_list_ecdna(full_connected_path):
    init_contig = form_normal_contig(contig_data[full_connected_path[0][1]])
    path_contig = [init_contig]
    full_connected_path_len = len(full_connected_path)
    vcnt = 0
    for i in range(1, full_connected_path_len):
        curr_contig = path_contig[-1]
        next_contig = contig_data[full_connected_path[i][1]]
        if curr_contig[CHR_NAM] == next_contig[CHR_NAM] and curr_contig[CTG_NAM] != next_contig[CTG_NAM] \
        or (curr_contig[CHR_NAM] == next_contig[CHR_NAM] and (full_connected_path[i-1][0] in (2, 3) or full_connected_path[i][0] in (2, 3))):
            dist = distance_checker(curr_contig, next_contig)
            if dist > 0:
                if curr_contig[CHR_END] < next_contig[CHR_STR]:
                    if curr_contig[CTG_NAM] != next_contig[CTG_NAM]:
                        vcnt += 1
                        new_contig = form_virtual_contig(curr_contig[CHR_LEN], curr_contig[CHR_END], next_contig[CHR_STR], curr_contig[CHR_NAM], vcnt)
                        new_next_contig = form_normal_contig(next_contig)
                        path_contig.append(new_contig)
                        path_contig.append(new_next_contig)
                    else:
                        new_next_contig = form_normal_contig(next_contig)
                        path_contig.append(new_next_contig)
                else:
                    if curr_contig[CTG_NAM] != next_contig[CTG_NAM]:
                        vcnt += 1
                        new_contig = form_virtual_contig(curr_contig[CHR_LEN], next_contig[CHR_END], curr_contig[CHR_STR], curr_contig[CHR_NAM], vcnt)
                        new_next_contig = form_normal_contig(next_contig)
                        path_contig.append(new_contig)
                        path_contig.append(new_next_contig)
                    else:
                        new_next_contig = form_normal_contig(next_contig)
                        path_contig.append(new_next_contig)
            elif curr_contig[CHR_END] == next_contig[CHR_STR] or curr_contig[CHR_STR] == next_contig[CHR_END]:
                new_next_contig = form_normal_contig(next_contig)
                path_contig.append(new_next_contig)
            else:
                if curr_contig[CTG_DIR] == '-':
                    curr_contig_rc = node_dir_transform(curr_contig)
                else:
                    curr_contig_rc = copy.deepcopy(curr_contig)
                if next_contig[CTG_DIR] == '-':
                    next_contig_rc = node_dir_transform(next_contig)
                else:
                    next_contig_rc = copy.deepcopy(next_contig)
                curr_contig_data = [curr_contig_rc[CHR_STR], curr_contig_rc[CHR_END],
                                    curr_contig_rc[CTG_STR], curr_contig_rc[CTG_END],
                                    curr_contig[9], curr_contig[10],
                                    curr_contig[-1][5:]]

                next_contig_origin = form_normal_contig(next_contig)
                
                next_contig_data = [next_contig_rc[CHR_STR], next_contig_rc[CHR_END],
                                    next_contig_rc[CTG_STR], next_contig_rc[CTG_END],
                                    next_contig_origin[9], next_contig_origin[10],
                                    next_contig_origin[-1][5:]]
                # try:
                #     assert(not (inclusive_checker_pair(tuple(curr_contig_data[0:2]), tuple(next_contig_data[0:2])) \
                #            or inclusive_checker_pair(tuple(next_contig_data[0:2]), tuple(curr_contig_data[0:2]))))
                # except:
                #     print(tuple(curr_contig_data[0:2]), tuple(next_contig_data[0:2]))
                #     print(index_file_path)
                if full_connected_path[i-1][0] in {2, 3}:
                    is_index_inc = (next_contig[CTG_DIR]=='+' and full_connected_path[i][0]==1) or \
                                (next_contig[CTG_DIR]=='-' and full_connected_path[i][0]==0)
                elif full_connected_path[i][0] in {2, 3}:
                    is_index_inc = (curr_contig[CTG_DIR]=='+' and full_connected_path[i-1][0]==1) or \
                            (curr_contig[CTG_DIR]=='-' and full_connected_path[i-1][0]==0)
                else:
                    assert(curr_contig[CTG_NAM] != next_contig[CTG_NAM])
                    is_index_inc = (curr_contig[CTG_DIR]=='+' and full_connected_path[i-1][0]==1) or \
                            (curr_contig[CTG_DIR]=='-' and full_connected_path[i-1][0]==0)

                if is_index_inc: # 1324 2413

                    curr_paf, next_paf = adjust_paf_overlap(curr_contig_data, next_contig_data)
                    
                    curr_result = form_adjusted_contig(curr_contig, curr_paf)
                    if curr_contig[CTG_DIR]=='-':
                        curr_result = node_dir_transform(curr_result)
                    next_result = form_adjusted_contig(next_contig_origin, next_paf)
                    if next_contig[CTG_DIR]=='-':
                        next_result = node_dir_transform(next_result)
                    path_contig[-1] = curr_result
                    path_contig.append(next_result)
                else:
                    next_paf, curr_paf = adjust_paf_overlap(next_contig_data, curr_contig_data)
                    
                    curr_result = form_adjusted_contig(curr_contig, curr_paf)
                    if curr_contig[CTG_DIR]=='-':
                        curr_result = node_dir_transform(curr_result)
                    next_result = form_adjusted_contig(next_contig_origin, next_paf)
                    if next_contig[CTG_DIR]=='-':
                        next_result = node_dir_transform(next_result)
                    path_contig[-1] = curr_result
                    path_contig.append(next_result)

                
        else:
            next_contig_origin = form_normal_contig(next_contig)
            path_contig.append(next_contig_origin)
    cnt = 0

    final_output_list = []
    for i in path_contig:
        cnt+=1
        cigar_str = 'cg:Z:' + cs_to_cigar(i[-1][5:])
        final_output_list.append("\t".join(list(map(str, i + [cigar_str]))))
    
    return final_output_list
    
def create_final_depth_paf(data):
    (key, key_cnt) = data
    key_type, key_val = key

    raw_contig_list = []
    if key_type == TEL_TYPE:
        # Single contig
        ctg_node_ind = key_val
        raw_contig_list.append((DIR_FOR, ctg_node_ind))
    elif key_type == CTG_IN_TYPE:
        (s_type, s), (e_type, e) = key_val
        assert(contig_data[s][CTG_NAM] == contig_data[e][CTG_NAM])
        if s<e:
            for i in range(s+1, e):
                raw_contig_list.append((DIR_FOR, i))
        else:
            for i in range(s-1, e, -1):
                raw_contig_list.append((DIR_BAK, i))
    elif key_type == BND_TYPE:
        (s_type, s), (e_type, e) = key_val
        
        if contig_data[s][CTG_NAM] != contig_data[e][CTG_NAM]:
            if nx.has_path(G, source=(DIR_OUT, s), target=(DIR_IN, e)):
                path = nx.shortest_path(G, source=(DIR_OUT, s), target=(DIR_IN, e), weight='weight')
                raw_contig_list.append((s_type, s))
                for node in path[1:-1]:
                    raw_contig_list.append((node[0], node[1]))
                raw_contig_list.append((e_type, e))
            else:
                raw_contig_list.append((s_type, s))
                raw_contig_list.append((e_type, e))
        else:
            if s<e:
                for i in range(s, e+1):
                    raw_contig_list.append((DIR_FOR, i))
            else:
                for i in range(s, e-1, -1):
                    raw_contig_list.append((DIR_BAK, i))
    else:
        # Never happen
        assert(False)

    if len(raw_contig_list) > 0:
        process_raw_contig_list(raw_contig_list, key_cnt)
    else:
        # Create empty file
        with open(f'{output_folder}/{key_cnt}.paf', 'w') as f:
            pass

def create_final_depth_paf_ecdna(ecdna_circuit, save_path):
    for idx, circuit in enumerate(ecdna_circuit):
        assert(len(circuit) == 4)
        data = []
        s1, e1, s2, e2 = circuit
        data.append((CTG_IN_TYPE, ((DIR_FOR, s1), (DIR_FOR, e1))))
        data.append((BND_TYPE, ((DIR_FOR, e1), (DIR_BAK, e2))))
        data.append((CTG_IN_TYPE, ((DIR_BAK, e2), (DIR_BAK, s2))))
        data.append((BND_TYPE, ((DIR_BAK, s2), (DIR_FOR, s1))))
        circuit_paf = []
        for (key_type, key_val) in data:
            raw_contig_list = []
            if key_type == TEL_TYPE:
                # Single contig
                ctg_node_ind = key_val
                raw_contig_list.append((DIR_FOR, ctg_node_ind))
            elif key_type == CTG_IN_TYPE:
                (s_type, s), (e_type, e) = key_val
                assert(contig_data[s][CTG_NAM] == contig_data[e][CTG_NAM])
                if s<e:
                    for i in range(s+1, e):
                        raw_contig_list.append((DIR_FOR, i))
                else:
                    for i in range(s-1, e, -1):
                        raw_contig_list.append((DIR_BAK, i))
            elif key_type == BND_TYPE:
                (s_type, s), (e_type, e) = key_val
                
                if contig_data[s][CTG_NAM] != contig_data[e][CTG_NAM]:
                    if nx.has_path(G, source=(DIR_OUT, s), target=(DIR_IN, e)):
                        path = nx.shortest_path(G, source=(DIR_OUT, s), target=(DIR_IN, e), weight='weight')
                        raw_contig_list.append((s_type, s))
                        for node in path[1:-1]:
                            raw_contig_list.append((node[0], node[1]))
                        raw_contig_list.append((e_type, e))
                    else:
                        raw_contig_list.append((s_type, s))
                        raw_contig_list.append((e_type, e))
                else:
                    if s<e:
                        for i in range(s, e+1):
                            raw_contig_list.append((DIR_FOR, i))
                    else:
                        for i in range(s, e-1, -1):
                            raw_contig_list.append((DIR_BAK, i))
            else:
                # Never happen
                assert(False)

            if len(raw_contig_list) > 0:
                circuit_paf += process_raw_contig_list_ecdna(raw_contig_list)
            else:
                # Create empty file
                    pass
        with open(f'{save_path}/{idx+1}.paf', 'wt') as f:
            for i in circuit_paf:
                print(i, file=f)

def create_final_depth_paf_type2(type2_ins_del, PREFIX):
    type2_ins, type2_del = type2_ins_del
    merge_save_path = f'{PREFIX}/11_ref_ratio_outliers/type2_ins'
    os.makedirs(merge_save_path, exist_ok=True)
    for idx, circuit in enumerate(list(type2_del) + list(type2_ins)):
        assert(len(circuit) == 4)
        data = []
        s1, e1, s2, e2 = circuit
        data.append((BND_TYPE, ((DIR_FOR, e1), (DIR_FOR, s2))))
        ins_paf = []
        for (key_type, key_val) in data:
            raw_contig_list = []
            if key_type == TEL_TYPE:
                # Single contig
                ctg_node_ind = key_val
                raw_contig_list.append((DIR_FOR, ctg_node_ind))
            elif key_type == CTG_IN_TYPE:
                (s_type, s), (e_type, e) = key_val
                assert(contig_data[s][CTG_NAM] == contig_data[e][CTG_NAM])
                if s<e:
                    for i in range(s+1, e):
                        raw_contig_list.append((DIR_FOR, i))
                else:
                    for i in range(s-1, e, -1):
                        raw_contig_list.append((DIR_BAK, i))
            elif key_type == BND_TYPE:
                (s_type, s), (e_type, e) = key_val
                
                if contig_data[s][CTG_NAM] != contig_data[e][CTG_NAM]:
                    if nx.has_path(G, source=(DIR_OUT, s), target=(DIR_IN, e)):
                        path = nx.shortest_path(G, source=(DIR_OUT, s), target=(DIR_IN, e), weight='weight')
                        raw_contig_list.append((s_type, s))
                        for node in path[1:-1]:
                            raw_contig_list.append((node[0], node[1]))
                        raw_contig_list.append((e_type, e))
                    else:
                        raw_contig_list.append((s_type, s))
                        raw_contig_list.append((e_type, e))
                else:
                    if s<e:
                        for i in range(s, e+1):
                            raw_contig_list.append((DIR_FOR, i))
                    else:
                        for i in range(s, e-1, -1):
                            raw_contig_list.append((DIR_BAK, i))
            else:
                # Never happen
                assert(False)

            if len(raw_contig_list) > 0:
                ins_paf += process_raw_contig_list_ecdna(raw_contig_list)
            else:
                # Create empty file
                    pass
        with open(f'{PREFIX}/11_ref_ratio_outliers/type2_ins/{idx+1}.paf', 'wt') as f:
            for i in ins_paf:
                print(i, file=f)


def rev_dir(d):
    if d == DIR_IN:
        return DIR_OUT
    elif d == DIR_OUT:
        return DIR_IN
    elif d == DIR_FOR:
        return DIR_BAK
    elif d == DIR_BAK:
        return DIR_FOR
    else:
        assert(False)

def rev_ind(ind):
    d, i = ind
    return (rev_dir(d), i)

def sorted_ind(f, b):
    if f[1] > b[1]:
        return False, (rev_ind(b), rev_ind(f))
    else:
        return True, (f, b)

def is_sorted_ind(key):
    if sorted_ind(key) == key:
        return True
    return False

def count_groups(lst):
    return [len(list(group)) for key, group in groupby(lst)]

def get_final_paf_name_from_index(index_file_path):
    final_paf_list = index_file_path.split('/')
    cnt = final_paf_list[-1].split('.')[0]
    final_paf_list[-3] = '20_depth'
    final_paf_path = '/'.join(final_paf_list[:-1]) + f'/{cnt}.paf'

    return final_paf_path

# 10_breakend_graph_build_paths
def initial_graph_build(contig_data : list) -> list :
    contig_data_size = len(contig_data)
    adjacency = [[[] for _ in range(contig_data_size)], [[] for _ in range(contig_data_size)]]
    curr_contig_st = 0
    while curr_contig_st<contig_data_size:
        curr_contig_ed = contig_data[curr_contig_st][CTG_ENDND]
        if contig_data[curr_contig_st][CTG_TYP] == 3:
            for j in range(curr_contig_st+1, curr_contig_ed+1):
                adjacency[DIR_FOR][j-1].append([DIR_FOR, j, 0])
                adjacency[DIR_BAK][j].append([DIR_BAK, j-1, 0])
        curr_contig_st = curr_contig_ed + 1
    return adjacency

def node_label(contig_data : list) -> list :
    contig_data_size = len(contig_data)
    label = [0 for _ in range(contig_data_size)]
    for i in range(0, contig_data_size) :
        contig_start = contig_data[contig_data[i][CTG_STRND]]
        contig_end = contig_data[contig_data[i][CTG_ENDND]]
        if contig_data[i][CTG_TYP] == 3 \
            and contig_data[i][CHR_NAM] == contig_start[CHR_NAM]:
            if i in (contig_data[i][CTG_STRND], contig_data[i][CTG_ENDND]):
                label[i] = "1-3"
            else:
                if contig_data[i][CTG_DIR] == contig_start[CTG_DIR]:
                    if contig_data[i][CTG_DIR] == '+':
                        if contig_start[CHR_END] <= contig_data[i][CHR_STR] \
                        and contig_data[i][CHR_END] <= contig_end[CHR_STR] :
                            label[i] = "1-1"
                    else:
                        if contig_start[CHR_STR] >= contig_data[i][CHR_END] \
                        and contig_data[i][CHR_STR] >= contig_end[CHR_END] :
                            label[i] = "1-2"
    return label

def graph_build(contig_adjacency : list, linkability_label : list) -> list :
    contig_data_size = len(contig_data)
    gap_contig_adjacency = [[[[], []] for __ in range(contig_data_size)] for _ in range(2)]
    for i in range(contig_data_size):
        if linkability_label[i] != 0:
            k = contig_data[i][CTG_ENDND]+1
            for j in range(k, contig_data_size):
                if linkability_label[j] != 0 \
                and contig_data[i][CHR_NAM] == contig_data[j][CHR_NAM] \
                and contig_data[i][CTG_NAM] != contig_data[j][CTG_NAM] :
                    if (inclusive_checker(contig_data[i], contig_data[j]) or inclusive_checker(contig_data[j], contig_data[i]))\
                    and (contig_data[i][CTG_ENDND] - contig_data[i][CTG_STRND] == 0 or contig_data[j][CTG_ENDND] - contig_data[j][CTG_STRND] == 0):
                        continue
                    dist = distance_checker(contig_data[i], contig_data[j])

                # 일단 dist가 0일 때부터. 즉, 겹칠 경우.
                    if dist == 0:
                        # i->j 방향: dir
                        if contig_data[i][CHR_END] <= contig_data[j][CHR_END]:
                            dir1 = "inc"
                        elif contig_data[j][CHR_END] < contig_data[i][CHR_END]:
                            dir1 = "dec"
                        if contig_data[i][CHR_STR] <= contig_data[j][CHR_STR]:
                            dir2 = "inc"
                        elif contig_data[j][CHR_STR] < contig_data[i][CHR_STR]:
                            dir2 = "dec"
                        # if both end have consistency
                        if dir1 == dir2:
                            if dir1 == "inc":
                                if contig_data[i][CTG_DIR]=='+' and contig_data[j][CTG_DIR]=='+':
                                    contig_adjacency[DIR_FOR][i].append([DIR_FOR, j, dist])
                                    contig_adjacency[DIR_BAK][j].append([DIR_BAK, i, dist])
                                elif contig_data[i][CTG_DIR]=='-' and contig_data[j][CTG_DIR]=='-':
                                    contig_adjacency[DIR_FOR][j].append([DIR_FOR, i, dist])
                                    contig_adjacency[DIR_BAK][i].append([DIR_BAK, j, dist])
                                elif contig_data[i][CTG_DIR]=='+' and contig_data[j][CTG_DIR]=='-':
                                    contig_adjacency[DIR_FOR][i].append([DIR_BAK, j, dist])
                                    contig_adjacency[DIR_FOR][j].append([DIR_BAK, i, dist])
                                elif contig_data[i][CTG_DIR]=='-' and contig_data[j][CTG_DIR]=='+':
                                    contig_adjacency[DIR_BAK][i].append([DIR_FOR, j, dist])
                                    contig_adjacency[DIR_BAK][j].append([DIR_FOR, i, dist])
                            else:
                                if contig_data[i][CTG_DIR]=='+' and contig_data[j][CTG_DIR]=='+':
                                    contig_adjacency[DIR_FOR][j].append([DIR_FOR, i, dist])
                                    contig_adjacency[DIR_BAK][i].append([DIR_BAK, j, dist])
                                elif contig_data[i][CTG_DIR]=='-' and contig_data[j][CTG_DIR]=='-':
                                    contig_adjacency[DIR_FOR][i].append([DIR_FOR, j, dist])
                                    contig_adjacency[DIR_BAK][j].append([DIR_BAK, i, dist])
                                elif contig_data[i][CTG_DIR]=='+' and contig_data[j][CTG_DIR]=='-':
                                    contig_adjacency[DIR_BAK][i].append([DIR_FOR, j, dist])
                                    contig_adjacency[DIR_BAK][j].append([DIR_FOR, i, dist])
                                elif contig_data[i][CTG_DIR]=='-' and contig_data[j][CTG_DIR]=='+':
                                    contig_adjacency[DIR_FOR][i].append([DIR_BAK, j, dist])
                                    contig_adjacency[DIR_FOR][j].append([DIR_BAK, i, dist])
                        # else: each end have different sign of increment
                        # -> one node is included in the other
                        # thus, we should determine by mean val.
                        else:
                            if contig_data[i][CTG_DIR] == contig_data[j][CTG_DIR]:
                                contig_adjacency[DIR_FOR][i].append([DIR_FOR, j, dist])
                                contig_adjacency[DIR_BAK][j].append([DIR_BAK, i, dist])
                                contig_adjacency[DIR_FOR][j].append([DIR_FOR, i, dist])
                                contig_adjacency[DIR_BAK][i].append([DIR_BAK, j, dist])
                            else:
                                contig_adjacency[DIR_FOR][i].append([DIR_BAK, j, dist])
                                contig_adjacency[DIR_FOR][j].append([DIR_BAK, i, dist])
                                contig_adjacency[DIR_BAK][i].append([DIR_FOR, j, dist])
                                contig_adjacency[DIR_BAK][j].append([DIR_FOR, i, dist])
                    else:
                        # i에서 j로 증가하는 경우.
                        if contig_data[i][CHR_END] < contig_data[j][CHR_STR]:
                            if contig_data[i][CTG_DIR]=='+' and contig_data[j][CTG_DIR]=='+':
                                if len(gap_contig_adjacency[DIR_FOR][i][0])==0:
                                    gap_contig_adjacency[DIR_FOR][i][0] = [DIR_FOR, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][i][0] = [DIR_FOR, j, dist]

                                if len(gap_contig_adjacency[DIR_BAK][j][0])==0:
                                    gap_contig_adjacency[DIR_BAK][j][0] = [DIR_BAK, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][j][0] = [DIR_BAK, i, dist]
                            
                            elif contig_data[i][CTG_DIR]=='-' and contig_data[j][CTG_DIR]=='-':
                                if len(gap_contig_adjacency[DIR_FOR][j][0])==0:
                                    gap_contig_adjacency[DIR_FOR][j][0] = [DIR_FOR, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][j][0] = [DIR_FOR, i, dist]

                                if len(gap_contig_adjacency[DIR_BAK][i][0])==0:
                                    gap_contig_adjacency[DIR_BAK][i][0] = [DIR_BAK, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][i][0] = [DIR_BAK, j, dist]
                            
                            elif contig_data[i][CTG_DIR]=='+' and contig_data[j][CTG_DIR]=='-':
                                if len(gap_contig_adjacency[DIR_FOR][i][0])==0:
                                    gap_contig_adjacency[DIR_FOR][i][0] = [DIR_BAK, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][i][0] = [DIR_BAK, j, dist]

                                if len(gap_contig_adjacency[DIR_FOR][j][0])==0:
                                    gap_contig_adjacency[DIR_FOR][j][0] = [DIR_BAK, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][j][0] = [DIR_BAK, i, dist]
                            
                            elif contig_data[i][CTG_DIR]=='-' and contig_data[j][CTG_DIR]=='+':
                                if len(gap_contig_adjacency[DIR_BAK][i][0])==0:
                                    gap_contig_adjacency[DIR_BAK][i][0] = [DIR_FOR, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][i][0] = [DIR_FOR, j, dist]

                                if len(gap_contig_adjacency[DIR_BAK][j][0])==0:
                                    gap_contig_adjacency[DIR_BAK][j][0] = [DIR_FOR, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][j][0] = [DIR_FOR, i, dist]
                        else:
                            if contig_data[j][CTG_DIR]=='+' and contig_data[i][CTG_DIR]=='+':
                                if len(gap_contig_adjacency[DIR_FOR][j][0])==0:
                                    gap_contig_adjacency[DIR_FOR][j][0] = [DIR_FOR, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][j][0] = [DIR_FOR, i, dist]

                                if len(gap_contig_adjacency[DIR_BAK][i][0])==0:
                                    gap_contig_adjacency[DIR_BAK][i][0] = [DIR_BAK, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][i][0] = [DIR_BAK, j, dist]
                            
                            elif contig_data[j][CTG_DIR]=='-' and contig_data[i][CTG_DIR]=='-':
                                if len(gap_contig_adjacency[DIR_FOR][i][0])==0:
                                    gap_contig_adjacency[DIR_FOR][i][0] = [DIR_FOR, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][i][0] = [DIR_FOR, j, dist]

                                if len(gap_contig_adjacency[DIR_BAK][j][0])==0:
                                    gap_contig_adjacency[DIR_BAK][j][0] = [DIR_BAK, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][j][0] = [DIR_BAK, i, dist]
                            
                            elif contig_data[j][CTG_DIR]=='+' and contig_data[i][CTG_DIR]=='-':
                                if len(gap_contig_adjacency[DIR_FOR][j][0])==0:
                                    gap_contig_adjacency[DIR_FOR][j][0] = [DIR_BAK, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][j][0] = [DIR_BAK, i, dist]

                                if len(gap_contig_adjacency[DIR_FOR][i][0])==0:
                                    gap_contig_adjacency[DIR_FOR][i][0] = [DIR_BAK, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][i][0] = [DIR_BAK, j, dist]
                            
                            elif contig_data[j][CTG_DIR]=='-' and contig_data[i][CTG_DIR]=='+':
                                if len(gap_contig_adjacency[DIR_BAK][j][0])==0:
                                    gap_contig_adjacency[DIR_BAK][j][0] = [DIR_FOR, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][j][0] = [DIR_FOR, i, dist]

                                if len(gap_contig_adjacency[DIR_BAK][i][0])==0:
                                    gap_contig_adjacency[DIR_BAK][i][0] = [DIR_FOR, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][i][0] = [DIR_FOR, j, dist]
    for i in range(2):
        for j in range(contig_data_size):
            for k in range(2):
                if len(gap_contig_adjacency[i][j][k])>0:
                    contig_adjacency[i][j].append(gap_contig_adjacency[i][j][k])
    return contig_adjacency

def edge_optimization(contig_data : list, contig_adjacency : list) -> list :
    contig_pair_nodes = defaultdict(list)
    contig_data_size = len(contig_data)
    optimized_adjacency = [[[] for _ in range(contig_data_size)], [[] for _ in range(contig_data_size)]]
    for _ in range(2):
        for i in range(contig_data_size):
            for edge in contig_adjacency[_][i]:
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
    for i in range(2):
        for j in range(contig_data_size):
            first_contig_name = contig_data[j][CTG_NAM]
            fcn_int = contig_data[j][CTG_STRND]
            for edge in contig_adjacency[i][j]:
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

    return optimized_adjacency

def find_using_node(contig_data, node_label):
    contig_data_size = len(contig_data)
    using_node = []
    for i in range(contig_data_size):
        if node_label[i] != 0:
            using_node.append(i)
    return using_node

def extract_bnd_contig(contig_data : list) -> set:
    s = 0
    contig_data_size = len(contig_data)
    bnd_contig = set()
    while s<contig_data_size:
        e = contig_data[s][CTG_ENDND]
        if contig_data[s][CTG_TYP] in set([1, 2]): # 2 넣기
            bnd_contig.add(contig_data[s][CTG_NAM])
        s = e+1
    return bnd_contig

def extract_telomere_connect_contig(telo_info_path : str) -> list:
    telomere_connect_contig = []
    with open(telo_info_path) as f:
        for curr_data in f:
            curr_data = curr_data.rstrip()
            temp_list = curr_data.split("\t")
            chr_info = temp_list[0]
            contig_id = ast.literal_eval(temp_list[1])
            telomere_connect_contig.append(contig_id[1])
    
    return telomere_connect_contig

def extract_nclose_node(nclose_path: str) -> list:
    nclose_list = []
    with open(nclose_path, "r") as f:
        for line in f:
            line = line.split()
            nclose_list.append(int(line[1]))
            nclose_list.append(int(line[2]))
    return nclose_list

def connect_nclose_telo(contig_data : list, using_node : list, type_3_graph : dict, nclose_nodes : list, telo_nodes : list) -> dict:
    # Telomere
    telo_nodes = telo_nodes + nclose_nodes
    full_bnd_graph = type_3_graph
    for telo_node_idx in telo_nodes:
        for type_3_idx in using_node:
            telo_node = contig_data[telo_node_idx]
            telo_node_len = telo_node[CTG_ENDND] - telo_node[CTG_STRND] + 1
            type_3 = contig_data[type_3_idx]
            type_3_dir = type_3[CTG_DIR]
            if telo_node_idx == type_3_idx and telo_node_len > 1:
                if telo_node_idx == telo_node[CTG_STRND]:
                    full_bnd_graph[(DIR_OUT, telo_node_idx)].append([DIR_FOR, telo_node_idx+1, 0])
                    full_bnd_graph[(DIR_BAK, telo_node_idx+1)].append([DIR_IN, telo_node_idx, 0])
                elif telo_node_idx == telo_node[CTG_ENDND]:
                    full_bnd_graph[(DIR_OUT, telo_node_idx)].append([DIR_BAK, telo_node_idx-1, 0])
                    full_bnd_graph[(DIR_FOR, telo_node_idx-1)].append([DIR_IN, telo_node_idx, 0])
            if telo_node[CTG_NAM] != type_3[CTG_NAM] \
            and telo_node[CHR_NAM] == type_3[CHR_NAM]:
                if telo_node[CHR_END] <= type_3[CHR_END]:
                    dir1 = "inc"
                elif type_3[CHR_END] < telo_node[CHR_END]:
                    dir1 = "dec"
                if telo_node[CHR_STR] <= type_3[CHR_STR]:
                    dir2 = "inc"
                elif type_3[CHR_STR] < telo_node[CHR_STR]:
                    dir2 = "dec"
                dist = distance_checker(telo_node, type_3)
                # if both end have consistency
                if dir1 == dir2:
                    if dir1 == "inc":
                        if type_3_dir=='+':
                            full_bnd_graph[(DIR_OUT, telo_node_idx)].append([DIR_FOR, type_3_idx, dist])
                            full_bnd_graph[(DIR_BAK, type_3_idx)].append([DIR_IN, telo_node_idx, dist])
                        elif type_3_dir=='-':
                            full_bnd_graph[(DIR_FOR, type_3_idx)].append([DIR_IN, telo_node_idx, dist])
                            full_bnd_graph[(DIR_OUT, telo_node_idx)].append([DIR_BAK, type_3_idx, dist])
                    else:
                        if type_3_dir=='+':
                            full_bnd_graph[(DIR_FOR, type_3_idx)].append([DIR_IN, telo_node_idx, dist])
                            full_bnd_graph[(DIR_OUT, telo_node_idx)].append([DIR_BAK, type_3_idx, dist])
                        elif type_3_dir=='-':
                            full_bnd_graph[(DIR_OUT, telo_node_idx)].append([DIR_FOR, type_3_idx, dist])
                            full_bnd_graph[(DIR_BAK, type_3_idx)].append([DIR_IN, telo_node_idx, dist])
                # else: each end have different sign of increment
                # -> one node is included in the other
                # thus, we should determine by mean val.
                else:
                    full_bnd_graph[(DIR_OUT, telo_node_idx)].append([DIR_FOR, type_3_idx, dist])
                    full_bnd_graph[(DIR_BAK, type_3_idx)].append([DIR_IN, telo_node_idx, dist])
                    full_bnd_graph[(DIR_OUT, telo_node_idx)].append([DIR_BAK, type_3_idx, dist])
                    full_bnd_graph[(DIR_FOR, type_3_idx)].append([DIR_IN, telo_node_idx, dist])
    
    return full_bnd_graph


parser = argparse.ArgumentParser(description="Find breakend contigs with contig data and map data")

# 위치 인자 정의
parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("--pandepth_loc", 
                    help="Pandepth binary location.", type=str, default="pandepth")

parser.add_argument("-t", "--thread", 
                    help="Number of thread", type=int)

parser.add_argument("--progress", 
                    help="Show progress bar", action='store_true')

args = parser.parse_args()

# t = "21_run_depth.py /home/hyunwoo/ACCtools-pipeline/90_skype_run/COLO829/20_alignasm/COLO829.ctg.aln.paf.ppc.paf 30_skype_pipe/COLO829_17_00_42/ --pandepth_loc ./PanDepth/bin/pandepth -t 128"
# args = parser.parse_args(t.split()[1:])

PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path

contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
contig_data_size = len(contig_data)

THREAD = args.thread
PREFIX = args.prefix
# 20_fill_path code

with open(f'{PREFIX}/paf_file_path.pkl', 'rb') as f:
    PAF_FILE_PATH = pkl.load(f)

paf_file = []
for paf_loc in PAF_FILE_PATH:
    paf_file.append(import_data(paf_loc))

# 10_breakend_graph_build_paths code

TELO_CONNECT_NODES_INFO_PATH = f"{args.prefix}/telomere_connected_list.txt"
NCLOSE_FILE_PATH = f"{args.prefix}/nclose_nodes_index.txt"

init_graph_adjacency = initial_graph_build(contig_data)
link_label = node_label(contig_data)
ord_graph_adjacency = graph_build(init_graph_adjacency, link_label)

with open(f"{PREFIX}/non_opt_bnd_connect_graph.txt", "wt") as f:
    for dir_ in range(2):
        for node in range(contig_data_size):
            print((dir_, node), ":", ord_graph_adjacency[dir_][node], file=f)

opt_ord_graph_adjacency = edge_optimization(contig_data, ord_graph_adjacency)

type_3_graph = defaultdict(list)
with open(f"{PREFIX}/opt_bnd_connect_graph.txt", "wt") as f:
    for dir_ in range(2):
        for node in range(contig_data_size):
            print((dir_, node), ":", opt_ord_graph_adjacency[dir_][node], file=f)
            for edge in opt_ord_graph_adjacency[dir_][node]:
                type_3_graph[(dir_, node)].append(edge)

using_node = find_using_node(contig_data, link_label)

bnd_contig = extract_bnd_contig(contig_data)
telo_contig = extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH)
nclose_nodes = extract_nclose_node(NCLOSE_FILE_PATH)

bnd_connected_graph = connect_nclose_telo(contig_data, using_node, type_3_graph, nclose_nodes, telo_contig)

os.makedirs(f'{PREFIX}/11_ref_ratio_outliers/ecdna', exist_ok=True)
with open(f'{PREFIX}/ecdna_circuit_data.pkl', 'rb') as f:
    ecdna_circuit, raw_nclose_nodes = pkl.load(file=f)

with open(f'{PREFIX}/conjoined_type4_ins_del.pkl', 'rb') as f:
    type4_ins_del = pkl.load(file=f)

raw_nclose_nodes_list = []
for vl in raw_nclose_nodes.values():
    for v in vl:
        raw_nclose_nodes_list.extend(list(v))

G = nx.DiGraph()

for node in telo_contig + raw_nclose_nodes_list:
    G.add_node((DIR_IN, node))
    G.add_node((DIR_OUT, node))

for node in bnd_connected_graph:
    G.add_node(node)
    for edge in bnd_connected_graph[node]:
        G.add_node(tuple(edge[:-1]))

for node in bnd_connected_graph:
    for edge in bnd_connected_graph[node]:
        G.add_weighted_edges_from([(node, tuple(edge[:-1]), edge[-1])])


create_final_depth_paf_ecdna(ecdna_circuit, f'{PREFIX}/11_ref_ratio_outliers/ecdna')

create_final_depth_paf_type2(type4_ins_del, PREFIX)

# 21_run_depth

G = nx.DiGraph()

for node in telo_contig + nclose_nodes:
    G.add_node((DIR_IN, node))
    G.add_node((DIR_OUT, node))

for node in bnd_connected_graph:
    G.add_node(node)
    for edge in bnd_connected_graph[node]:
        G.add_node(tuple(edge[:-1]))

for node in bnd_connected_graph:
    for edge in bnd_connected_graph[node]:
        G.add_weighted_edges_from([(node, tuple(edge[:-1]), edge[-1])])


key_cnt = 0
key2int = dict()
int2key = dict()

final_paf_vec_data = []

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
back_front_folder_path = glob.glob(RATIO_OUTLIER_FOLDER+"*")

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)
chr_chr_folder_path = path_list_dict.keys()

key_ord_list = []
index_data_list = []
with ProcessPoolExecutor(max_workers=THREAD) as executor:
    futures = [executor.submit(get_key_from_index_folder, folder_path) for folder_path in chr_chr_folder_path]

    for future in tqdm(as_completed(futures), total=len(futures), desc='Analyse index to split paf',
                       disable=not sys.stdout.isatty() and not args.progress):
        k, idd = future.result()

        key_ord_list.append(k)
        index_data_list.extend(idd)

for _, key_list in index_data_list:
    key_ind_data = []
    for key in key_list:
        if key not in key2int:
            key2int[key] = key_cnt
            int2key[key_cnt] = key
            key_cnt += 1

output_folder = f'{PREFIX}/21_pat_depth'
os.makedirs(output_folder, exist_ok=True)

with ProcessPoolExecutor(max_workers=THREAD) as executor:
    futures = [executor.submit(create_final_depth_paf, (key, key_cnt)) for key, key_cnt in key2int.items()]

    for future in tqdm(as_completed(futures), total=len(futures), desc='Split paf to run depth efficiently',
                       disable=not sys.stdout.isatty() and not args.progress):
        future.result()

paf_ans_list = []
for index_file_path, key_list in index_data_list:
    final_paf_path = get_final_paf_name_from_index(index_file_path)

    key_list = [key2int[k] for k in key_list]
    paf_ans_list.append((final_paf_path, key_list))

def get_paf_run(paf_loc):
    paf_base = os.path.splitext(paf_loc)[0]

    retry = 0
    while retry < PANDEPTH_RETRY:
        result = subprocess.run([args.pandepth_loc, '-w', str(int(DEPTH_WINDOW)),'-t', str(DEPTH_THREAD), '-i', paf_loc, '-o', paf_base],
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, capture_output=False)
        
        if result.returncode == 0:
            return
        else:
            retry += 1
            logging.info(f"{paf_loc} : Pandepth failed retry {retry} time")
    
    logging.info(f"{paf_loc} : Pandepth failed retry limit exceed")
    raise Exception("Pandepth failed")

with ProcessPoolExecutor(max_workers=THREAD) as executor:
    futures = []
    paf_paths = glob.glob(output_folder + "/*.paf")
    for paf_loc in paf_paths:
        futures.append(executor.submit(get_paf_run, paf_loc))

    for folder_path in back_front_folder_path:
        paf_paths = glob.glob(folder_path + "/*.paf")
        for paf_loc in paf_paths:
            futures.append(executor.submit(get_paf_run, paf_loc))
    
    # 제출된 작업들이 완료될 때까지 진행 상황을 tqdm으로 표시합니다.
    for future in tqdm(as_completed(futures), total=len(futures), desc='Run PanDepth for each path file',
                       disable=not sys.stdout.isatty() and not args.progress):
        future.result()

with open(f'{PREFIX}/path_di_data.pkl', 'rb') as f:
    path_di_list_dict = pkl.load(f)

dep_list = []
for k in key_ord_list:
    dep_list.extend(path_di_list_dict[k])

assert(len(dep_list) == len(paf_ans_list))
dep_sort_ind_list = np.argsort(-np.asarray(dep_list)).tolist()


dep_sort_list = []
paf_sort_ans_list = []

for sort_ind in dep_sort_ind_list:
    dep_sort_list.append(dep_list[sort_ind])
    paf_sort_ans_list.append(paf_ans_list[sort_ind])

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'wb') as f:
    pkl.dump((paf_sort_ans_list, list(key2int.values()), int2key, dep_sort_list), f)
