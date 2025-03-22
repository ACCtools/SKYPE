import argparse
import pandas as pd
import numpy as np
import copy
import itertools
import os
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p'
)
logging.info("00_Contig_Preprocessing start")

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
CTG_GLOBALIDX = 21

DIR_FOR = 1
DIR_BAK = 0

K = 1000
M = 1000 * K
INF = 10000000000
CONTIG_MINIMUM_SIZE = 100*K
BND_CONTIG_BOUND = 0.1
RPT_BND_CONTIG_BOUND = 0.2
CHROMOSOME_COUNT = 23
MAPQ_BOUND = 60
TELOMERE_EXPANSION = 5 * K
CENSAT_COMPRESSABLE_THRESHOLD = 1e6

FORCE_TELOMERE_THRESHOLD = 10*K
TELOMERE_CLUSTER_THRESHOLD = 500*K
SUBTELOMERE_LENGTH = 500*K

MIN_FLANK_SIZE_BP = 1*M

WEAK_BREAKEND_CEN_RATIO_THRESHOLD = 1.3
BREAKEND_CEN_RATIO_THRESHOLD = 1.5

def dbg() :
    print("hi")

def import_data(file_path : list) -> list :
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

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    chr_data_file.close()
    return chr_len

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

def import_repeat_data(file_path : str) -> dict :
    fai_file = open(file_path, "r")
    repeat_data = defaultdict(list)
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        repeat_data[temp_list[0]].append((int(temp_list[1]), int(temp_list[2])))
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

def distance_checker(node_a : tuple, node_b : tuple) -> int :
    if max(int(node_a[CHR_STR]), int(node_b[CHR_STR])) < min(int(node_a[CHR_END]), int(node_b[CHR_END])):
        return 0   
    else:
        return min(abs(int(node_b[CHR_STR]) - int(node_a[CHR_END])), abs(int(node_b[CHR_END]) - int(node_a[CHR_STR])))

def inclusive_checker(contig_node : tuple, telomere_node : tuple) -> bool :
    if int(telomere_node[CHR_STR]) <= int(contig_node[CHR_STR]) and int(contig_node[CHR_END]) <= int(telomere_node[CHR_END]):
        return True
    else:
        return False
    
def telo_distance_checker(node: tuple, telo: tuple) -> int :
    return min(abs(telo[CHR_STR] - node[CHR_END]), abs(telo[CHR_END] - node[CHR_STR]))
    
def label_node(contig_data : list, telo_data) -> list :
    label = []
    contig_data_size = len(contig_data)
    for i in range(contig_data_size):
        checker = 0
        for j in telo_data[contig_data[i][CHR_NAM]]:
            if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])) == 0:
                inclusive_label = ""
                if inclusive_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])):
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
                    if inclusive_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1]+SUBTELOMERE_LENGTH)):
                        inclusive_label = "in"
                    label.append((contig_data[i][CHR_NAM], j[2]+inclusive_label))
                    checker = 1
                    break
            else:
                if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0] - SUBTELOMERE_LENGTH, j[1])) == 0:
                    inclusive_label = ""
                    if inclusive_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0] - SUBTELOMERE_LENGTH, j[1])):
                        inclusive_label = "in"
                    label.append((contig_data[i][CHR_NAM], j[2]+inclusive_label))
                    checker = 1
                    break
        if checker==0:
            label.append(('0', '0'))
    return label

def label_repeat_node(contig_data : list, repeat_data) -> list :
    label = []
    contig_data_size = len(contig_data)
    for i in range(contig_data_size):
        try:
            k = repeat_data[contig_data[i][CHR_NAM]]
        except:
            label.append(('0', '0'))
            continue
        flag = False
        for j in k:
            if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])) == 0:
                inclusive_label = ""
                if inclusive_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])):
                    inclusive_label = "in"
                label.append((contig_data[i][CHR_NAM], "r" + inclusive_label))
                flag = True
                break
        if not flag:
            label.append(('0', '0'))
    return label

def preprocess_telo(contig_data : list, node_label : list) -> list :
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
            escape = False
            if front_telo_bound == curr_contig_ed+1:
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
                and contig_data[st][CTG_TELDIR] == '0' \
                and st not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[st], temp_contig) \
                and contig_data[st][CTG_DIR]=='+':
                    telo_connect_node = st
                    telo_dist = telo_distance_checker(contig_data[st], temp_contig)
                    telo_dir = '+'
                if contig_data[ed][CHR_NAM] == now_telo[0:-1] \
                and contig_data[ed][CTG_TELDIR] == '0' \
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
                and contig_data[st][CTG_TELDIR] == '0' \
                and st not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[st], temp_contig) \
                and contig_data[st][CTG_DIR]=='-':
                    telo_connect_node = st
                    telo_dist = telo_distance_checker(contig_data[st], temp_contig)
                    telo_dir = '-'
                if contig_data[ed][CHR_NAM] == now_telo[0:-1] \
                and contig_data[ed][CTG_TELDIR] == '0' \
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

def edge_optimization(contig_data : list, contig_adjacency : list, telo_dict : dict) -> list :
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
            for edge in optimized_adjacency[j][i]:
                if telo_name[-1]=='f':
                    if contig_data[edge[1]][CHR_STR]<=TELOMERE_CLUSTER_THRESHOLD:
                        if now_edge[0]<0:
                            now_edge = edge
                        else:
                            if contig_data[now_edge[1]][CHR_STR] > contig_data[edge[1]][CHR_STR]:
                                now_edge = edge
                            elif contig_data[now_edge[1]][CHR_STR] == contig_data[edge[1]][CHR_STR]:
                                if contig_data[now_edge[1]][CTG_LEN] > contig_data[edge[1]][CTG_LEN]:
                                    now_edge = edge
                    else:
                        using_edge.append(edge)
                else:
                    if contig_data[edge[1]][CHR_END]>=telo_range[1]-TELOMERE_CLUSTER_THRESHOLD:
                        if now_edge[0]<0:
                            now_edge = edge
                        else:
                            if contig_data[now_edge[1]][CHR_END] < contig_data[edge[1]][CHR_END]:
                                now_edge = edge
                            elif contig_data[now_edge[1]][CHR_END] == contig_data[edge[1]][CHR_END]:
                                if contig_data[now_edge[1]][CTG_LEN] == contig_data[edge[1]][CTG_LEN]:
                                    now_edge = edge
                    else:
                        using_edge.append(edge)
            if now_edge == [-1, 0, 0]:
                pass
            else:
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

    return telo_connected_set, telo_connected_dict, telo_connected_graph_dict

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

def calc_chukji(contig_data : list) -> dict:
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
    ed = contig_data[st][CTG_ENDND]
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
        if telo_label[i-1][0] != '0'  or ((i-1) in telo_connect_info):
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
        if telo_label[i-1][0] != '0'  or ((i-1) in telo_connect_info):
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

def alt_preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list, telo_connect_info : set) -> list :
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
    for i in range(1, contig_data_size+1):
        len_count[contig_data[i-1][CTG_NAM]]+=contig_data[i-1][CHR_END]-contig_data[i-1][CHR_STR]
        cnt+=1
        if telo_label[i-1][0] != '0' or ((i-1) in telo_connect_info):
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
                    checker = 3
                else:
                    if is_front_back_repeat:
                        bound = RPT_BND_CONTIG_BOUND
                    else:
                        bound = BND_CONTIG_BOUND
                    if abs(ref_qry_ratio[curr_contig_name]-1) >= bound:
                        checker = 4
            if chrM_flag:
                checker = 0   
            if (checker>0 and checker != 3) or is_telo:
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

def chr2int(x):
    chrXY2int = {'chrX' : 24, 'chrY' : 25}
    if x in chrXY2int:
        return chrXY2int[x]
    else:
        return int(x[3:])
    
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

def main():
    parser = argparse.ArgumentParser(description="Process file paths for telomere analysis.")

    # 위치 인자 정의
    parser.add_argument("paf_file_path", 
                        help="Path to the original PAF file.")
    parser.add_argument("telomere_bed_path", 
                        help="Path to the telomere information file.")
    parser.add_argument("reference_fai_path", 
                        help="Path to the chromosome information file.")
    parser.add_argument("repeat_bed_path", 
                        help="Path to the chromosome repeat information file.")
    parser.add_argument("censat_bed_path", 
                        help="Path to the censat repeat information file.")
    parser.add_argument("main_stat_path", 
                        help="Path to the main stat file.")
    parser.add_argument("prefix", 
                    help="Pefix for pipeline")
    parser.add_argument("--alt", 
                        help="Path to an alternative PAF file (optional).")
    parser.add_argument("--progress", 
                        help="Show progress bar", action='store_true')

    # 인자 파싱
    args = parser.parse_args()

    # t = "00_Contig_Preprocessing.py 20_acc_pipe/U2OS_telo.p/U2OS_telo.p.aln.paf public_data/chm13v2.0_telomere.bed public_data/chm13v2.0.fa.fai public_data/chm13v2.0_repeat.m.bed public_data/chm13v2.0_censat_v2.1.m.bed /home/hyunwoo/51g_cancer_denovo/51_depth_data/U2OS_telo.win.stat.gz 30_skype_pipe/U2OS_telo_23_19_30 --alt 20_acc_pipe/U2OS_telo.r/U2OS_telo.r.aln.paf".split()
    # args = parser.parse_args(t[1:])

    original_node_count = 0

    PREFIX = args.prefix

    os.makedirs(PREFIX, exist_ok=True)

    PAF_FILE_PATH = []
    if args.alt is None:
        PAF_FILE_PATH = [args.paf_file_path]
    else:
        PAF_FILE_PATH = [args.paf_file_path, args.alt]
    
    TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
    CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
    PREPROCESSED_PAF_FILE_PATH = (PAF_FILE_PATH[0] +'.ppc.paf')
    REPEAT_INFO_FILE_PATH = args.repeat_bed_path
    CENSAT_PATH = args.censat_bed_path
    main_stat_loc = args.main_stat_path

    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
    telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)
    repeat_data = import_repeat_data(REPEAT_INFO_FILE_PATH)
    repeat_censat_data = import_censat_repeat_data(CENSAT_PATH)


    df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
    df = df.query('chr != "chrM"')

    cen_vtg_contig = find_breakend_centromere(repeat_censat_data, chr_len, df)

    telo_dict = defaultdict(list)
    for _ in telo_data:
        telo_dict[_[0]].append(_[1:])



    contig_data = import_data(PAF_FILE_PATH[0])

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
        contig_data = import_data(PAF_FILE_PATH[1])
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
        preprocess_contig_type, \
        preprocess_terminal_nodes, \
        alt_len_counter = alt_preprocess_contig(new_contig_data, new_node_telo_label, ref_qry_ratio, new_node_repeat_label, telcon_set)
        new_contig_data = new_contig_data[0:-1]
        contig_data_size = len(new_contig_data)
        bias = cnt
        preprocess_result = set(preprocess_result)
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

        alt_final_using_contig, alt_final_ctg_typ, alt_final_preprocess_terminal_nodes, _ = alt_preprocess_contig(alt_final_contig, alt_final_telo_node_label, alt_final_ref_qry_ratio, alt_final_repeat_node_label, alt_final_telo_connect)
        
        alt_final_contig = alt_final_contig[:-1]
        
        bias = len(real_final_contig)

        real_alt_final_contig = []

        alt_final_using_contig = set(alt_final_using_contig)
        for i in range(0, len(alt_final_contig)):
            if alt_final_contig[i][CTG_NAM] in alt_final_using_contig:
                alt_final_contig[i][CTG_TYP] = alt_final_ctg_typ[alt_final_contig[i][CTG_NAM]]
                alt_final_contig[i][CTG_STRND] = alt_final_preprocess_terminal_nodes[alt_final_contig[i][CTG_NAM]][0] + bias
                alt_final_contig[i][CTG_ENDND] = alt_final_preprocess_terminal_nodes[alt_final_contig[i][CTG_NAM]][1] + bias
                real_alt_final_contig.append(alt_final_contig[i])
        
        real_final_contig = real_final_contig + real_alt_final_contig
        total_len = len(real_final_contig)
    
    
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

    with open("subtelo.txt", "wt") as f:
        for i in subtelo_ppc_node:
            for j in i:
                print(j, end="\t", file=f)
            print(file=f)

    adjacency = initial_graph_build(real_final_contig, telo_bound_dict)

    telo_connected_node, telo_connected_dict, _ = edge_optimization(real_final_contig, adjacency, telo_bound_dict)

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

    add_node_count = len(final_break_contig) + len(final_cen_vtg_contig)
    s = 0
    r_l = len(real_final_contig)
    logging.info(f"Original PAF file length : {original_node_count}")
    logging.info(f"Final preprocessed PAF file length: {r_l}")
    logging.info(f"Number of virtual contigs added on preprocessing : {add_node_count}")
    contig = set()
    while s<r_l:
        e=real_final_contig[s][CTG_ENDND]
        chkcensat = False
        chktelo=False
        for i in range(s, e+1):
            if real_final_contig[i][CTG_TELCHR] != '0':
                chkcensat = True
            elif real_final_contig[i][CTG_CENSAT] != '0':
                chktelo = True
        
        if chkcensat and chktelo:
            contig.add(contig_data[s][CTG_NAM])
        
        s = e+1
    
    # for i in contig:
    #     print(i)

    adjacency = initial_graph_build(real_final_contig, telo_bound_dict)

    telo_connected_node, telo_connected_dict, telo_connected_graph_dict = edge_optimization(real_final_contig, adjacency, telo_bound_dict)

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


    with open(PREPROCESSED_PAF_FILE_PATH, "wt") as f:
        for i in real_final_contig:
            for j in i:
                print(j, end="\t", file=f)
            print("", file=f)

if __name__ == "__main__":
    main()