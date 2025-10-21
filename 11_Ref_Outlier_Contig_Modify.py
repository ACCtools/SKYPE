import os
import sys
import pickle as pkl
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *

import re
import logging
import argparse

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("11_Ref_Outlier_Contig_Modify start")

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

CHUKJI_LIMIT = 100*K

def import_origin_data(file_path : list) -> list :
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


def import_data(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.rstrip().split("\t")
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END, \
                          CHR_LEN, CHR_STR, CHR_END, \
                          CTG_MAPQ, CTG_TYP, CTG_STRND, CTG_ENDND,]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    paf_file.close()
    return contig_data

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    chr_data_file.close()
    return chr_len


def calculate_single_contig_ref_ratio(contig_data : list):
    total_ref_len = 0
    ref_st_ed = 0
    if contig_data[0][CTG_DIR] == '+':
        estimated_ref_len = contig_data[-1][CHR_END] - contig_data[0][CHR_STR]
        ref_st_ed = (contig_data[0][CHR_STR], contig_data[-1][CHR_END])
    else:
        estimated_ref_len = contig_data[0][CHR_END] - contig_data[-1][CHR_STR]
        ref_st_ed = (contig_data[-1][CHR_STR], contig_data[0][CHR_END])
    for node in contig_data:
        total_ref_len += node[CHR_END] - node[CHR_STR]
    return estimated_ref_len/total_ref_len, ref_st_ed

def distance_checker(node_a : tuple, node_b : tuple) -> int :
    if max(int(node_a[0]), int(node_b[0])) < min(int(node_a[1]), int(node_b[1])):
        return 0   
    else:
        return min(abs(int(node_b[0]) - int(node_a[1])), abs(int(node_b[1]) - int(node_a[0])))
    

def inclusive_checker(tuple_a : tuple, tuple_b : tuple) -> bool :
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False
    
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



parser = argparse.ArgumentParser(description="Find reference depth of reverse contig")

parser.add_argument("paf_file_path", 
                        help="Path to the original PAF file.")

parser.add_argument("reference_fai_path", 
                        help="Path to the chromosome information file.")

parser.add_argument("ppc_paf_file_path", 
                        help="Path to the preprocessed PAF file.")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("--alt", 
                        help="Path to an alternative PAF file (optional).")

args = parser.parse_args()

paf_file = []
if args.alt is None:
    paf_file = [import_origin_data(args.paf_file_path)]
else:
    paf_file = [import_origin_data(args.paf_file_path), import_origin_data(args.alt)]

PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
TYPE_4_VECTOR_PATH = f"{args.prefix}/11_ref_ratio_outliers"
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
os.makedirs(f"{TYPE_4_VECTOR_PATH}/front_jump", exist_ok=True)
os.makedirs(f"{TYPE_4_VECTOR_PATH}/back_jump", exist_ok=True)

chr_data = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
contig_data = import_data(PREPROCESSED_PAF_FILE_PATH)
contig_data_size = len(contig_data)

with open(f'{args.prefix}/conjoined_type4_ins_del.pkl', 'rb') as f:
    type4_ins, type4_del = pkl.load(file=f)
s = 0
cntfj = 0
cntbj = 0
while s<contig_data_size:
    e = contig_data[s][CTG_ENDND]

    if contig_data[s][CTG_TYP] == 4:
        chr_name = contig_data[s][CHR_NAM]
        chr_len = chr_data[chr_name]
        rat, ref_st_ed = calculate_single_contig_ref_ratio(contig_data[s:e+1])
        if abs(ref_st_ed[1]-ref_st_ed[0]) > CHUKJI_LIMIT or contig_data[s][CTG_LEN] > 2 * CHUKJI_LIMIT:
            if rat > 0:
                cntfj += 1
                with open(f"{TYPE_4_VECTOR_PATH}/front_jump/{cntfj}.paf", "wt") as f:
                    for i in range(s, e+1):
                        glob_paf_idx = int(contig_data[i][CTG_GLOBALIDX][0])
                        glob_idx = int(contig_data[i][CTG_GLOBALIDX][2:])
                        cigar_str = 'cg:Z:' + cs_to_cigar(paf_file[glob_paf_idx][glob_idx][-1][5:])
                        for j in paf_file[glob_paf_idx][glob_idx]:
                            print(j, end="\t", file=f)
                        print(cigar_str, end="\t", file=f)
                        print("", file=f)
                with open(f"{TYPE_4_VECTOR_PATH}/front_jump/{cntfj}_base.paf", "wt") as f:
                    N = ref_st_ed[1] - ref_st_ed[0]
                    virtual_contig = [f"base_contig_1", N, 0, N]
                    virtual_contig += ['+', chr_name, chr_len, ref_st_ed[0], ref_st_ed[1]]
                    virtual_contig += [N, N, 0]
                    virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
                    cigar_str = 'cg:Z:' + cs_to_cigar(virtual_contig[-1][5:])
                    virtual_contig += [cigar_str]
                    for j in virtual_contig:
                        print(j, end="\t", file=f)
                    print("", file=f)
            else:
                cntbj += 1
                with open(f"{TYPE_4_VECTOR_PATH}/back_jump/{cntbj}.paf", "wt") as f:
                    for i in range(s, e+1):
                        glob_paf_idx = int(contig_data[i][CTG_GLOBALIDX][0])
                        glob_idx = int(contig_data[i][CTG_GLOBALIDX][2:])
                        cigar_str = 'cg:Z:' + cs_to_cigar(paf_file[glob_paf_idx][glob_idx][-1][5:])
                        for j in paf_file[glob_paf_idx][glob_idx]:
                            print(j, end="\t", file=f)
                        print(cigar_str, end="\t", file=f)
                        print("", file=f)
                with open(f"{TYPE_4_VECTOR_PATH}/back_jump/{cntbj}_base.paf", "wt") as f:
                    N = ref_st_ed[0] - ref_st_ed[1]
                    virtual_contig = [f"base_contig_1", N, 0, N]
                    virtual_contig += ['+', chr_name, chr_len, ref_st_ed[1], ref_st_ed[0]]
                    virtual_contig += [N, N, 0]
                    virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
                    cigar_str = 'cg:Z:' + cs_to_cigar(virtual_contig[-1][5:])
                    virtual_contig += [cigar_str]
                    for j in virtual_contig:
                        print(j, end="\t", file=f)
                    print("", file=f)
            
    s = e+1


type2_indel_cnt = 0

for s1, e1, s2, e2 in list(type4_del) + list(type4_ins):
    type2_indel_cnt += 1
    chr_name = contig_data[s1][CHR_NAM]
    chr_len = chr_data[chr_name]
    rat, ref_st_ed = calculate_single_contig_ref_ratio([contig_data[s1], contig_data[e2]])
    if rat > 0:
        cntfj+=1
        with open(f"{TYPE_4_VECTOR_PATH}/front_jump/{cntfj}_type2_merge_{type2_indel_cnt}.paf", "wt") as f:
            for i in (s1, e2):
                glob_paf_idx = int(contig_data[i][CTG_GLOBALIDX][0])
                glob_idx = int(contig_data[i][CTG_GLOBALIDX][2:])
                cigar_str = 'cg:Z:' + cs_to_cigar(paf_file[glob_paf_idx][glob_idx][-1][5:])
                for j in paf_file[glob_paf_idx][glob_idx]:
                    print(j, end="\t", file=f)
                print(cigar_str, end="\t", file=f)
                print("", file=f)
            
        with open(f"{TYPE_4_VECTOR_PATH}/front_jump/{cntfj}_base.paf", "wt") as f:
                N = ref_st_ed[1] - ref_st_ed[0]
                virtual_contig = [f"base_contig_1", N, 0, N]
                virtual_contig += ['+', chr_name, chr_len, ref_st_ed[0], ref_st_ed[1]]
                virtual_contig += [N, N, 0]
                virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
                cigar_str = 'cg:Z:' + cs_to_cigar(virtual_contig[-1][5:])
                virtual_contig += [cigar_str]
                for j in virtual_contig:
                    print(j, end="\t", file=f)
                print("", file=f)
    else:
        cntbj+=1
        with open(f"{TYPE_4_VECTOR_PATH}/back_jump/{cntbj}_type2_merge_{type2_indel_cnt}.paf", "wt") as f:
            for i in (s1, e2):
                glob_paf_idx = int(contig_data[i][CTG_GLOBALIDX][0])
                glob_idx = int(contig_data[i][CTG_GLOBALIDX][2:])
                cigar_str = 'cg:Z:' + cs_to_cigar(paf_file[glob_paf_idx][glob_idx][-1][5:])
                for j in paf_file[glob_paf_idx][glob_idx]:
                    print(j, end="\t", file=f)
                print(cigar_str, end="\t", file=f)
                print("", file=f)
            
        with open(f"{TYPE_4_VECTOR_PATH}/back_jump/{cntbj}_base.paf", "wt") as f:
                N = ref_st_ed[0] - ref_st_ed[1]
                virtual_contig = [f"base_contig_1", N, 0, N]
                virtual_contig += ['+', chr_name, chr_len, ref_st_ed[1], ref_st_ed[0]]
                virtual_contig += [N, N, 0]
                virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
                cigar_str = 'cg:Z:' + cs_to_cigar(virtual_contig[-1][5:])
                virtual_contig += [cigar_str]
                for j in virtual_contig:
                    print(j, end="\t", file=f)
                print("", file=f)


logging.info(f"Forward-directed outlier contig count : {cntfj}")
logging.info(f"Backward-directed outlier contig count : {cntbj}")
logging.info(f"Total count : {cntfj + cntbj}")
