import re
import argparse
import os
import sys
import ast
import glob
import copy

from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor

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
CTG_GLOBALIDX=21

DIR_FOR = 1
DIR_BAK = 0
DIR_IN = 3
DIR_OUT = 2
INF = 1000000000
BUFFER = 10000000
CHROMOSOME_COUNT = 23
K = 1000

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
        curr_contig.rstrip()
        temp_list = curr_contig.split("\t")
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END, \
                          CHR_LEN, CHR_STR, CHR_END, \
                          CTG_MAPQ, CTG_TYP, CTG_STRND, CTG_ENDND,]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    return contig_data

def import_index_path(file_path : str) -> list:
    index_file = open(file_path, "r")
    index_data = []
    for curr_index in index_file:
        curr_index = curr_index.rstrip()
        index_data.append(ast.literal_eval(curr_index))
    index_file.close()
    return index_data

class PAFOverlapError(Exception):
    pass

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

# def compute_aln_stats(cs):
#     """
#     cs 문자열을 기반으로 매칭 길이(mat_len)와 정렬 길이(aln_len)을 계산.
#     - ":" 연산: ref와 query 모두 소비 (길이만큼)
#     - "+": 삽입 → aln_len에만 추가 (query 소비)
#     - "-": 결실 → aln_len에만 추가 (ref 소비)
#     - "*": 치환 → aln_len +1 (ref와 query 각각 1 소비)
#     """
#     pattern = re.compile(r"(:\d+|\*[a-z]{2}|[+\-][A-Za-z]+)")
#     mat, aln = 0, 0
#     for token in pattern.findall(cs):
#         if token.startswith(":"):
#             n = int(token[1:])
#             mat += n
#             aln += n
#         elif token.startswith("+"):
#             seq = token[1:]
#             aln += len(seq)
#         elif token.startswith("-"):
#             seq = token[1:]
#             aln += len(seq)
#         elif token.startswith("*"):
#             aln += 1
#     return mat, aln

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
        raise PAFOverlapError("두 PAF의 reference 구간이 겹치지 않습니다.")
    
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


def distance_checker(node_a : tuple, node_b : tuple) -> int :
    # check if overlap
    if max(node_a[CHR_STR], node_b[CHR_STR]) < min(node_a[CHR_END], node_b[CHR_END]):
        return 0   
    else:
        # return distance of two ranges
        return min(abs(node_b[CHR_STR] - node_a[CHR_END]), abs(node_b[CHR_END] - node_a[CHR_STR]))

def node_dir_transform(node):
    transform_node = list(copy.deepcopy(node))
    transform_node[CTG_STR] = node[CTG_LEN] - node[CTG_END]
    transform_node[CTG_END] = node[CTG_LEN] - node[CTG_STR]
    transform_node[CTG_DIR] = '+' if node[CTG_DIR] == '+' else '-'
    return transform_node

def form_virtual_contig(chr_len, chr_str, chr_end, chr_name, vcnt):
    virtual_contig = []
    N = chr_end - chr_str
    virtual_contig = [f"virtual_contig_{vcnt}", N, 0, N]
    virtual_contig += ['+', chr_name, chr_len, chr_str, chr_end]
    virtual_contig += [N, N, 0]
    virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
    return virtual_contig

def form_normal_contig(node, paf_file):
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
    new_node[13] = 'cs:Z:' + paf_data[6]
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

              
parser = argparse.ArgumentParser(description="Find breakend contigs with contig data and map data")

# 위치 인자 정의

parser.add_argument("paf_file_path", 
                    help="Path to the original PAF file.")

parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("--alt", 
                    help="Path to an alternative PAF file (optional).")

# t = "python 20_Fill_Path.py 20_acc_pipe/KKU-100.p/KKU-100.p.aln.paf 20_acc_pipe/KKU-100.p/KKU-100.p.aln.paf.ppc.paf 30_skype_pipe/KKU-100_06_27_56 --alt 20_acc_pipe/KKU-100.a/KKU-100.a.aln.paf"
# t = t.split()
args = parser.parse_args()
# args = parser.parse_args(t[2:])
PAF_FILE_PATH = []
paf_file = []
if args.alt is None:
    PAF_FILE_PATH = [args.paf_file_path]
    paf_file.append(import_data(PAF_FILE_PATH[0]))
else:
    PAF_FILE_PATH = [args.paf_file_path, args.alt]
    primary_paf = import_data(PAF_FILE_PATH[0])
    alt_paf = import_data(PAF_FILE_PATH[1])
    paf_file.append(primary_paf)
    paf_file.append(alt_paf)


PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path

contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
contig_data_size = len(contig_data)
cnt=0
PATH_FILE_FOLDER = f"{args.prefix}/10_fill/"
TYPE_4_PATH_FILE_FOLDER = f"{args.prefix}/11_ref_ratio_outliers/"

# fpath = "13_breakend_connect_graph/HuH-28.bnd_path.00:05:20/chr18f_chr18b"
# ifpath = "13_breakend_connect_graph/HuH-28.bnd_path.00:05:20/chr18f_chr18b/1.index.txt"

# fill_path((fpath, ifpath, 1, contig_data, paf_file))

# exit(1)


def fill_path(index_file_path):
    temp_file = index_file_path.split("/")[:-1]
    temp_file[2] = '20_depth'
    folder_path = "/".join(temp_file)
    cnt = int(index_file_path.split(".")[-3].split("/")[-1])
    os.makedirs(f"{folder_path}", exist_ok=True)
    with open(f"{folder_path}/{cnt}.paf", "wt") as f:
        full_connected_path = import_index_path(index_file_path)
        init_contig = form_normal_contig(contig_data[full_connected_path[0][1]], paf_file)
        path_contig = [init_contig]
        full_connected_path_len = len(full_connected_path)
        vcnt = 0
        for i in range(1, full_connected_path_len):
            if i==full_connected_path_len-1:
                t=1
            curr_contig = path_contig[-1]
            next_contig = contig_data[full_connected_path[i][1]]
            if curr_contig[CHR_NAM] == next_contig[CHR_NAM] and curr_contig[CTG_NAM] != next_contig[CTG_NAM] \
            or (curr_contig[CHR_NAM] == next_contig[CHR_NAM] and (full_connected_path[i-1][0] in (2, 3) or full_connected_path[i][0] in (2, 3))):
                dist = distance_checker(curr_contig, next_contig)
                if dist > 0:
                    vcnt += 1
                    if curr_contig[CHR_END] < next_contig[CHR_STR]:
                        new_contig = form_virtual_contig(curr_contig[CHR_LEN], curr_contig[CHR_END], next_contig[CHR_STR], curr_contig[CHR_NAM], vcnt)
                        new_next_contig = form_normal_contig(next_contig, paf_file)
                        path_contig.append(new_contig)
                        path_contig.append(new_next_contig)
                    else:
                        new_contig = form_virtual_contig(curr_contig[CHR_LEN], next_contig[CHR_END], curr_contig[CHR_STR], curr_contig[CHR_NAM], vcnt)
                        new_next_contig = form_normal_contig(next_contig, paf_file)
                        path_contig.append(new_contig)
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
                                        curr_contig[13][5:]]

                    next_contig_origin = form_normal_contig(next_contig, paf_file)
                    
                    next_contig_data = [next_contig_rc[CHR_STR], next_contig_rc[CHR_END],
                                        next_contig_rc[CTG_STR], next_contig_rc[CTG_END],
                                        next_contig_origin[9], next_contig_origin[10],
                                        next_contig_origin[13][5:]]

                    
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
                next_contig_origin = form_normal_contig(next_contig, paf_file)
                path_contig.append(next_contig_origin)
        for i in path_contig:
            cigar_str = 'cg:Z:' + cs_to_cigar(i[-1][5:])
            print("\t".join(list(map(str, i + [cigar_str]))), file=f)
    
    return


with ProcessPoolExecutor(max_workers=48) as executor:
    futures = []
    chr_chr_folder_path = glob.glob(PATH_FILE_FOLDER+"*")
    # print(chr_chr_folder_path)
    # chr_chr_folder_path = ["30_skype_pipe/Caki-1_03_52_36/10_fill/chr17f_chr17b"]
    for folder_path in chr_chr_folder_path:
        index_file_paths = glob.glob(folder_path + "/*index*")
        for index_file_path in index_file_paths:
            futures.append(executor.submit(fill_path, index_file_path))
        

# 제출된 작업들이 완료될 때까지 진행 상황을 tqdm으로 표시합니다.
for future in tqdm(as_completed(futures), total=len(futures), desc='Fill gap and modify path data', disable=not sys.stdout.isatty()):
    pass
