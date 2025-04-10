import sys
import argparse
from collections import defaultdict, Counter
import ast
import pickle as pkl

import graph_tool.all as gt
import graph_tool

import networkx as nx

import copy
import shutil
import os

from multiprocessing import Pool
from tqdm import tqdm
from functools import partial
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
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

CHR_CHANGE_LIMIT_PREFIX = 5
DIR_CHANGE_LIMIT = 1
CENSAT_VISIT_LIMIT = 2

BND_OVERUSE_CNT = 2
PATH_MAJOR_COMPONENT = 3
NCLOSE_COMPRESS_LIMIT = 100*K
NCLOSE_MERGE_LIMIT = 1*K
ALL_REPEAT_NCLOSE_COMPRESS_LIMIT = 500*K
PATH_COMPRESS_LIMIT = 50*K
IGNORE_PATH_LIMIT = 50*K
NON_REPEAT_NOISE_RATIO=0.1

CENSAT_COMPRESSABLE_THRESHOLD = 1000*K

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

def div_orignial_paf(file_path : str) -> list:
    not_using_contig = set()
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            a = curr_contig.split("\t")
            assert(a[16][:5]=='tp:A:')
            if a[16] != 'tp:A:P':
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
    
def overlap_calculator(node_a : tuple, node_b : tuple) -> int :
    return min(abs(node_a[CHR_END] - node_b[CHR_STR]), abs(node_b[CHR_END] - node_a[CHR_STR]))

def inclusive_checker(tuple_a : tuple, tuple_b : tuple) -> bool :
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False

def extract_all_repeat_contig(contig_data : list, ctg_index : int, baseline : float = 0) -> set:
    contig_data_size = len(contig_data)
    rpt_con = set()
    s = 0
    while s<contig_data_size:
        e = contig_data[s][CTG_ENDND]
        flag = True
        total_ref_len = 0
        non_rpt_ref_len = 0
        for i in range(s, e+1):
            total_ref_len += contig_data[i][CHR_END] - contig_data[i][CHR_STR]
            if contig_data[i][ctg_index] == '0':
                non_rpt_ref_len += contig_data[i][CHR_END] - contig_data[i][CHR_STR]
            # if (i == s or i == e) and contig_data[i][ctg_index] == 'r':
            #     flag = False
        if non_rpt_ref_len == 0 or baseline > non_rpt_ref_len/total_ref_len:
            rpt_con.add(contig_data[s][CTG_NAM])
        s = e+1
    return rpt_con

def check_censat_contig(all_repeat_censat_con : set, ORIGNAL_PAF_LOC_LIST : list):
    div_repeat_paf_name = set()
    for i in ORIGNAL_PAF_LOC_LIST:
        div_repeat_paf_name.update(div_orignial_paf(i))
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

def calculate_single_contig_ref_ratio(contig_data : list) -> float:
    total_ref_len = 0
    if contig_data[0][CTG_DIR] == '+':
        estimated_ref_len = contig_data[-1][CHR_END] - contig_data[0][CHR_STR]
    else:
        estimated_ref_len = contig_data[0][CHR_END] - contig_data[-1][CHR_STR]
    for node in contig_data:
        total_ref_len += node[CHR_END] - node[CHR_STR]
    return estimated_ref_len/total_ref_len, total_ref_len

def extract_nclose_node(contig_data : list, bnd_contig : set, repeat_contig_name : set, \
                        censat_contig_name : set, repeat_censat_data : dict, ORIGNAL_PAF_LOC_LIST : list, \
                        telo_set : set, chr_len : dict) -> dict:
    s = 0
    fake_bnd = {}
    contig_data_size = len(contig_data)
    nclose_compress = defaultdict(list)
    nclose_start_compress = defaultdict(lambda : defaultdict(list))
    nclose_end_compress = defaultdict(lambda : defaultdict(list))
    censat_nclose_compress = set()
    nclose_dict = defaultdict(list)
    all_nclose_compress = defaultdict(list)
    telo_name_set = set()
    for i in telo_set:
        telo_name_set.add(contig_data[i][CTG_NAM])

    div_repeat_paf_name = set()
    for i in ORIGNAL_PAF_LOC_LIST:
        div_repeat_paf_name.update(div_orignial_paf(i))

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
                while st < e and [contig_data[st][CTG_DIR], contig_data[st][CHR_NAM]] == st_chr:
                    st+=1
                st-=1
                while ed > s and [contig_data[ed][CTG_DIR], contig_data[ed][CHR_NAM]] == ed_chr:
                    ed-=1
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
                    st_chr = (contig_data[st][CTG_DIR], contig_data[st][CHR_NAM])
                    ed_chr = (contig_data[st_idx][CTG_DIR], contig_data[st_idx][CHR_NAM])
                    nclose_list.append([st, st_chr, st_idx, ed_chr])
                    st_chr = (contig_data[ed_idx][CTG_DIR], contig_data[ed_idx][CHR_NAM])
                    ed_chr = (contig_data[ed][CTG_DIR], contig_data[ed][CHR_NAM])
                    nclose_list.append([ed_idx, st_chr, ed, ed_chr])
                else:
                    nclose_list.append([st, st_chr, ed, ed_chr])
                for nclose in nclose_list:
                    st = nclose[0]
                    ed = nclose[2]
                    nclose_front_const = 0
                    nclose_back_const = 0
                    for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                        if inclusive_checker(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                            nclose_front_const = 1
                            if censat_ref_range[0] < K or censat_ref_range[1] > chr_len[contig_data[st][CHR_NAM]]-K:
                                nclose_front_const = 2
                    for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                        if inclusive_checker(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
                            nclose_back_const = 1
                            if censat_ref_range[0] < K or censat_ref_range[1] > chr_len[contig_data[ed][CHR_NAM]]-K:
                                nclose_back_const = 2
                    if nclose_front_const + nclose_back_const >= 3:
                        continue

                    st_chr = nclose[1]
                    ed_chr = nclose[3]
                    upd_contig_name = contig_data[st][CTG_NAM]  
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
                                flag = False
                                break
                        if flag:
                            censat_st_chr = [st_chr[0], 0]
                            censat_ed_chr = [ed_chr[0], 0]
                            if contig_s[CTG_CENSAT] != '0':
                                cnt = 0
                                for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                                    if inclusive_checker(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                                        censat_st_chr[1] = contig_data[st][CHR_NAM] + "." + str(cnt)
                                        break
                                    cnt+=1
                            if contig_e[CTG_CENSAT] != '0':
                                cnt = 0
                                for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                                    if inclusive_checker(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
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
                                    if distance_checker(contig_data[st], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                        flag = True
                                        nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                    if distance_checker(contig_data[ed], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                        flag = True
                                        nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                    if flag:
                                        break
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
                                        if distance_checker(contig_data[st], dummy_list+i[1]) < NCLOSE_MERGE_LIMIT:
                                            flag = True
                                            nclose_start_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                        if distance_checker(contig_data[ed], dummy_list+i[2]) < NCLOSE_MERGE_LIMIT:
                                            flag = True
                                            nclose_end_compress[(st_chr, ed_chr)][i[0]].append(contig_data[s][CTG_NAM])
                                        if flag:
                                            break
                                    nclose_compress[(st_chr, ed_chr)].append(temp_list)
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
                                flag = False
                                break
                        if flag:
                            censat_st_chr = [st_chr[0], 0]
                            censat_ed_chr = [ed_chr[0], 0]
                            if contig_s[CTG_CENSAT] != '0':
                                cnt = 0
                                for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                                    if inclusive_checker(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                                        censat_st_chr[1] = contig_data[st][CHR_NAM] + "." + str(cnt)
                                        break
                                    cnt+=1
                            if contig_e[CTG_CENSAT] != '0':
                                cnt = 0
                                for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                                    if inclusive_checker(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
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
                                    flag = False
                                    break
                            if flag:
                                censat_st_chr = [st_chr[0], 0]
                                censat_ed_chr = [ed_chr[0], 0]
                                if contig_s[CTG_CENSAT] != '0':
                                    cnt = 0
                                    for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                                        if inclusive_checker(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                                            censat_st_chr[1] = contig_data[st][CHR_NAM] + "." + str(cnt)
                                            break
                                        cnt+=1
                                if contig_e[CTG_CENSAT] != '0':
                                    cnt = 0
                                    for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                                        if inclusive_checker(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
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
                                    flag = False
                                    break
                            if flag:
                                censat_st_chr = [st_chr[0], 0]
                                censat_ed_chr = [ed_chr[0], 0]
                                if contig_s[CTG_CENSAT] != '0':
                                    cnt = 0
                                    for censat_ref_range in repeat_censat_data[contig_data[st][CHR_NAM]]:
                                        if inclusive_checker(censat_ref_range, (contig_data[st][CHR_STR], contig_data[st][CHR_END])):
                                            censat_st_chr[1] = contig_data[st][CHR_NAM] + "." + str(cnt)
                                            break
                                        cnt+=1
                                if contig_e[CTG_CENSAT] != '0':
                                    cnt = 0
                                    for censat_ref_range in repeat_censat_data[contig_data[ed][CHR_NAM]]:
                                        if inclusive_checker(censat_ref_range, (contig_data[ed][CHR_STR], contig_data[ed][CHR_END])):
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
                                        nclose_compress[(ed_chr, st_chr)].append(temp_list)
                                        nclose_dict[contig_data[s][CTG_NAM]].append((st, ed))
        s = e+1
    return nclose_dict, nclose_start_compress, nclose_end_compress, fake_bnd, all_nclose_compress

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



            
parser = argparse.ArgumentParser(description="Find breakend contigs with contig data and map data")

# 위치 인자 정의
parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")
parser.add_argument("reference_fai_path", 
                        help="Path to the chromosome information file.")
parser.add_argument("censat_bed_path", 
                    help="Path to the censat repeat information file.")
parser.add_argument("prefix", 
                    help="Pefix for pipeline")
parser.add_argument("--orignal_paf_loc", nargs='+',
                    help="Orignal paf location to detect location (primary, alternative paf location)")
parser.add_argument("-t", "--thread", 
                    help="Number of thread", type=int)
parser.add_argument("--progress", 
                    help="Show progress bar", action='store_true')
parser.add_argument("--verbose", 
                    help="Enable index, paf output (Counld be slow at HDD)", action='store_true')

args = parser.parse_args()

# t = "02_Build_Breakend_Graph_Limited.py 20_acc_pipe/SNU-1196.p/SNU-1196.p.aln.paf.ppc.paf public_data/chm13v2.0.fa.fai public_data/chm13v2.0_censat_v2.1.m.bed 30_skype_pipe/SNU-1196_14_07_00 --orignal_paf_loc 20_acc_pipe/SNU-1196.p/SNU-1196.p.paf 20_acc_pipe/SNU-1196.r/SNU-1196.r.paf -t 128"
# args = parser.parse_args(t.split()[1:])

PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
ORIGNAL_PAF_LOC_LIST = args.orignal_paf_loc
CENSAT_PATH = args.censat_bed_path

chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
PREFIX = args.prefix
PRINT_IDX_FILE = args.verbose

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

rpt_con = extract_all_repeat_contig(contig_data, CTG_RPTCASE, NON_REPEAT_NOISE_RATIO)
rpt_censat_con = extract_all_repeat_contig(contig_data, CTG_CENSAT)
rpt_censat_con = check_censat_contig(rpt_censat_con, ORIGNAL_PAF_LOC_LIST)
repeat_censat_data = import_repeat_data(CENSAT_PATH)

bnd_contig = extract_bnd_contig(contig_data)


# Type 1, 2, 4에 대해서 
nclose_nodes, nclose_start_compress, nclose_end_compress, vctg_dict, all_nclose_comp = \
    extract_nclose_node(contig_data, bnd_contig, rpt_con, rpt_censat_con, repeat_censat_data, ORIGNAL_PAF_LOC_LIST, telo_set, chr_len)

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

with open(f"{PREFIX}/all_nclose_nodes_list.txt", "wt") as f:
    for i in all_nclose_type:
        print(f"{i[0]}, {i[1]}, {len(all_nclose_type[i])}", file=f)
        for pair in all_nclose_type[i]:
            contig_a = contig_data[pair[0]]
            contig_b = contig_data[pair[1]]
            list_a = [contig_a[CTG_NAM], contig_a[CTG_DIR], contig_a[CHR_STR], contig_a[CHR_END]]
            list_b = [contig_b[CTG_NAM], contig_b[CTG_DIR], contig_b[CHR_STR], contig_b[CHR_END]]
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


with open(f"{PREFIX}/compressed_nclose_nodes_list.txt", "wt") as f:
    for i in nclose_type:
        print(f"{i[0]}, {i[1]}, {len(nclose_type[i])}", file=f)
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
            list_a = [contig_a[CTG_NAM], contig_a[CTG_DIR], contig_a[CHR_STR], contig_a[CHR_END]]
            list_b = [contig_b[CTG_NAM], contig_b[CTG_DIR], contig_b[CHR_STR], contig_b[CHR_END]]
            if contig_a[CTG_NAM] in rpt_con:
                print(list_a, list_b, "all_repeat", file=f)
            else:
                print(list_a, list_b, file=f)
        print("", file=f)


nclose_node_count = 0


bnd_graph_adjacency = initialize_bnd_graph(contig_data, nclose_nodes, telo_contig)

with open(f"{PREFIX}/nclose_nodes_index.txt", "wt") as f: 
    for j in nclose_nodes:
        for i in nclose_nodes[j]:
            nclose_node_count += 2
            print(j, i[0], i[1], contig_data[i[0]][CTG_TYP], file=f)

with open(f"{PREFIX}/bnd_only_graph.txt", "wt") as f: 
    for i in bnd_graph_adjacency:
        print(i, bnd_graph_adjacency[i], file=f)


save_loc = PREFIX + '/00_raw'
logging.info("Now saving results in folder : " + PREFIX)
logging.info(f"NClose node count : {nclose_node_count}")
logging.info(f"Telomere connected node count : {telo_node_count}")

def make_graph(CHR_CHANGE_LIMIT):
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

def run_graph(data, CHR_CHANGE_LIMIT):
    path_list = []

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

                        if cnt >= PAT_PATH_LIMIT:
                            return ((src[0], tar), cnt, path_list)

    return ((src[0], tar), cnt, path_list)


tar_ind_list = []
for i in range(contig_data_size, contig_data_size + 2*CHROMOSOME_COUNT):
    for j in range(i, contig_data_size + 2*CHROMOSOME_COUNT):
        tar_ind_list.append((i, j))

def init_worker(shared_graph):
    global G
    G = shared_graph

THREAD=args.thread

tot_cnt = TOT_PATH_LIMIT
cnt_list = []

path_list_dict_data = []
while tot_cnt >= TOT_PATH_LIMIT and CHR_CHANGE_LIMIT_PREFIX > 0:
    path_list_dict_data.clear()

    tot_cnt = 0
    cnt_list = []

    # 각 반복마다 새로운 G 객체 생성
    G_obj = make_graph(CHR_CHANGE_LIMIT=CHR_CHANGE_LIMIT_PREFIX)

    # Pool 생성 시 initializer에 G_obj를 전달 (initargs는 반드시 튜플로)
    with Pool(processes=THREAD, initializer=init_worker, initargs=(G_obj,)) as pool:
        # imap을 사용하여 결과를 순차적으로 받아옴
        result_iterator = pool.imap_unordered(partial(run_graph, CHR_CHANGE_LIMIT=CHR_CHANGE_LIMIT_PREFIX), tar_ind_list)
        for rs in tqdm(result_iterator,
                        total=len(tar_ind_list),
                        desc=f'Build breakend construct graph <{CHR_CHANGE_LIMIT_PREFIX}>',
                        disable=not sys.stdout.isatty() and not args.progress):
            
            cnt_list.append((rs[0], rs[1]))
            tot_cnt += rs[1]

            # 한계치를 넘으면 즉시 종료
            if tot_cnt >= TOT_PATH_LIMIT:
                pool.terminate()  # 작업 중인 프로세스들을 즉시 종료
                shutil.rmtree(save_loc, ignore_errors=True)
                break

            path_list_dict_data.append((f'{rs[0][0]}_{rs[0][1]}', rs[2]))
    
    if tot_cnt >= TOT_PATH_LIMIT:
        logging.info(f'CHR_CHANGE_LIMIT : {CHR_CHANGE_LIMIT_PREFIX} failed')
        CHR_CHANGE_LIMIT_PREFIX -= 1


if CHR_CHANGE_LIMIT_PREFIX == 0:
    logging.info('Cancer is too divergent.')
    logging.info('Breakend path failed')
    exit(1)

path_list_dict = dict()
for k, v in path_list_dict_data:
    path_list_dict[k] = v

logging.info(f'CHR_CHANGE_LIMIT : {CHR_CHANGE_LIMIT_PREFIX} success')

cancer_prefix = os.path.basename(PREPROCESSED_PAF_FILE_PATH).split('.')[0]

with open(f'{PREFIX}/path_data.pkl', 'wb') as f:
    pkl.dump(path_list_dict, f)

with open(f'{PREFIX}/report.txt', 'a') as f:
    cnt = sum(i[1] for i in cnt_list)
    print(cancer_prefix, file=f)
    print(cnt, file=f)
    logging.info(f"Total path count : {cnt}")
    for (st, nd), c in sorted(cnt_list, key=lambda x:-x[1]):
        if c > 0:
            print(st, nd, c, file=f)
