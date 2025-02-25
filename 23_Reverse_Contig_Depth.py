import re
import argparse
import pandas as pd
from collections import defaultdict
from collections import Counter

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
    return contig_data


def calculate_single_contig_ref_ratio(contig_data : list) -> float:
    total_ref_len = 0
    if contig_data[0][CTG_DIR] == '+':
        estimated_ref_len = contig_data[-1][CHR_END] - contig_data[0][CHR_STR]
    else:
        estimated_ref_len = contig_data[0][CHR_END] - contig_data[-1][CHR_STR]
    for node in contig_data:
        total_ref_len += node[CHR_END] - node[CHR_STR]
    return estimated_ref_len/total_ref_len, total_ref_len

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


parser = argparse.ArgumentParser(description="Find reference depth of reverse contig")

parser.add_argument("ppc_paf_file_path", 
                        help="Path to the preprocessed PAF file.")

args = parser.parse_args()

PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path

contig_data = import_data(PREPROCESSED_PAF_FILE_PATH)
contig_data_size = len(contig_data)

s = 0
while s<contig_data_size:
    e = contig_data[s][CTG_ENDND]
    if contig_data[s][CTG_TYP] == 4:
        if contig_data[s][CTG_DIR]=='+':
            ref_range = (contig_data[s][CHR_STR], contig_data[e][CHR_END])
        else:
            ref_range = 
        points = defaultdict(int)
        ranges = []
        for i in range(s, e+1):
            if contig_data[i][CHR_NAM] == contig_data[s][CHR_NAM]:
                ranges.append((contig_data[i][CHR_STR], contig_data[i][CHR_END]))
    
        # 구간의 시작점과 끝점 기록
        for start, end in ranges:
            points[start] += 1
            points[end] -= 1
        
        # 정렬된 키 기준으로 등장 횟수 계산
        sorted_keys = sorted(points.keys())
        result = []
        count = 0
        
        for i in range(len(sorted_keys) - 1):
            count += points[sorted_keys[i]]
            result.append(((sorted_keys[i], sorted_keys[i+1]), count))
        
        
            
            
                    

    s = e+1



