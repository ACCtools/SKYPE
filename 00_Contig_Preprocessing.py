import os
import sys
import argparse
import pandas as pd
import copy

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

K = 1000
INF = 10000000000
CONTIG_MINIMUM_SIZE = 100*K
BND_CONTIG_BOUND = 0.1
RPT_BND_CONTIG_BOUND = 0.2
CHROMOSOME_COUNT = 23
MAPQ_BOUND = 60


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
            temp_list.append('f')
        else:
            if temp_list[1]>chr_len[temp_list[0]]/2:
                temp_list.append('b')
            else:
                temp_list.append('f')
        telo_data.append(tuple(temp_list))
    return telo_data[1:]

def import_repeat_data(file_path : str) -> dict :
    fai_file = open(file_path, "r")
    repeat_data = defaultdict(list)
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        repeat_data[temp_list[0]].append((int(temp_list[1]), int(temp_list[2])))
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
    report_case = {'A':[], 'B':[], 'C':[], 'ESC':[], 'ALL_TELO_NON_ESC':[]}
    contig_data_size = len(contig_data)
    curr_contig_st = 0
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    for i in range(1, contig_data_size+1):
        if contig_data[i-1][CTG_NAM] != contig_data[i][CTG_NAM]:
            curr_contig_ed = i-1
            front_telo_bound = curr_contig_st
            end_telo_bound = curr_contig_ed
            while node_label[front_telo_bound][0] != '0' and front_telo_bound<=curr_contig_ed:
                front_telo_bound+=1
            while node_label[end_telo_bound][0] != '0' and end_telo_bound>=curr_contig_st:
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
    total_contig_count = len(contig_data)
    mainflow_dict = {}
    st = 0
    ed = contig_data[st][CTG_ENDND]
    while st<total_contig_count:
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
    

def preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list) -> list :
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
        if telo_label[i-1][0] != '0':
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

def alt_preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list) -> list :
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
        if telo_label[i-1][0] != '0':
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
    
    parser.add_argument("--alt", 
                        help="Path to an alternative PAF file (optional).")

    # 인자 파싱
    args = parser.parse_args()

    # t = "python 00_Contig_Preprocessing.py 20_acc_pipe/OZ.p/OZ.p.aln.paf public_data/chm13v2.0_telomere.bed public_data/chm13v2.0.fa.fai public_data/chm13v2.0_repeat.m.bed public_data/chm13v2.0_censat_v2.1.m.bed --alt 20_acc_pipe/OZ.a/OZ.a.aln.paf".split()
    # args = parser.parse_args(t[2:])

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

    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
    telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)
    repeat_data = import_repeat_data(REPEAT_INFO_FILE_PATH)
    repeat_censat_data = import_repeat_data(CENSAT_PATH)

    telo_dict = defaultdict(list)
    for _ in telo_data:
        telo_dict[_[0]].append(_[1:])


    contig_data = import_data(PAF_FILE_PATH[0])

    node_label = label_node(contig_data, telo_dict)

    repeat_label = label_repeat_node(contig_data, repeat_data)

    telo_preprocessed_contig, report_case, telo_connect_info = preprocess_telo(contig_data, node_label)

    new_contig_data = []

    for i in telo_preprocessed_contig:
        temp_list = contig_data[i]
        if i in telo_connect_info:
            temp_list.append(telo_connect_info[i])
        else:
            temp_list.append("0")
        new_contig_data.append(temp_list)

    with open("telo_preprocess_contig.txt", "wt") as f:
        for i in new_contig_data:
            for j in i:
                print(j, end="\t", file=f)
            print("", file=f)
    
    new_node_repeat_label = label_repeat_node(new_contig_data, repeat_data)
    new_node_repeat_censat_label = label_repeat_node(new_contig_data, repeat_censat_data)
    new_node_telo_label = label_node(new_contig_data, telo_dict)

    ref_qry_ratio = calc_ratio(new_contig_data)
    preprocess_result, \
    preprocess_contig_type, \
    preprocess_terminal_nodes, \
    len_counter = preprocess_contig(new_contig_data, new_node_telo_label, ref_qry_ratio, new_node_repeat_label)
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

    final_using_contig, final_ctg_typ, final_preprocess_terminal_nodes, _ = preprocess_contig(final_contig, final_telo_node_label, final_ref_qry_ratio, final_contig_repeat_label)

    final_contig = final_contig[:-1]

    real_final_contig = []

    with open(PREPROCESSED_PAF_FILE_PATH, "wt") as f:
        final_using_contig = set(final_using_contig)
        for i in range(0, len(final_contig)):
            if final_contig[i][CTG_NAM] in final_using_contig:
                final_contig[i][CTG_TYP] = final_ctg_typ[final_contig[i][CTG_NAM]]
                final_contig[i][CTG_STRND] = final_preprocess_terminal_nodes[final_contig[i][CTG_NAM]][0]
                final_contig[i][CTG_ENDND] = final_preprocess_terminal_nodes[final_contig[i][CTG_NAM]][1]
                real_final_contig.append(final_contig[i])
                for j in final_contig[i]:
                    print(j, end="\t", file=f)
                print("", file=f)


    # alt process
    if args.alt != None:
        contig_data = import_data(PAF_FILE_PATH[1])
        node_label = label_node(contig_data, telo_dict)
        telo_preprocessed_contig, report_case, telo_connect_info = preprocess_telo(contig_data, node_label)
        new_contig_data = []
        for i in telo_preprocessed_contig:
            temp_list = contig_data[i]
            if i in telo_connect_info:
                temp_list.append(telo_connect_info[i])
            else:
                temp_list.append("0")
            new_contig_data.append(temp_list)
        alt_telo_ppc_contig = []
        new_node_telo_label = label_node(new_contig_data, telo_dict)
        new_node_repeat_label = label_repeat_node(new_contig_data, repeat_data)
        new_node_repeat_censat_label = label_repeat_node(new_contig_data, repeat_censat_data)

        ref_qry_ratio = calc_ratio(new_contig_data)
        preprocess_result, \
        preprocess_contig_type, \
        preprocess_terminal_nodes, \
        alt_len_counter = alt_preprocess_contig(new_contig_data, new_node_telo_label, ref_qry_ratio, new_node_repeat_label)
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

        alt_final_ref_qry_ratio = calc_ratio(alt_final_contig)

        alt_final_using_contig, alt_final_ctg_typ, alt_final_preprocess_terminal_nodes, _ = alt_preprocess_contig(alt_final_contig, alt_final_telo_node_label, alt_final_ref_qry_ratio, alt_final_repeat_node_label)
        
        alt_final_contig = alt_final_contig[:-1]
        
        bias = len(real_final_contig)

        real_alt_final_contig = []

        with open(PREPROCESSED_PAF_FILE_PATH, "a") as f:
            alt_final_using_contig = set(alt_final_using_contig)
            for i in range(0, len(alt_final_contig)):
                if alt_final_contig[i][CTG_NAM] in alt_final_using_contig:
                    alt_final_contig[i][CTG_TYP] = alt_final_ctg_typ[alt_final_contig[i][CTG_NAM]]
                    alt_final_contig[i][CTG_STRND] = alt_final_preprocess_terminal_nodes[alt_final_contig[i][CTG_NAM]][0] + bias
                    alt_final_contig[i][CTG_ENDND] = alt_final_preprocess_terminal_nodes[alt_final_contig[i][CTG_NAM]][1] + bias
                    for j in alt_final_contig[i]:
                        print(j, end="\t", file=f)
                    real_alt_final_contig.append(alt_final_contig[i])
                    print("", file=f)
        real_final_contig = real_final_contig + real_alt_final_contig
        real_final_ratio, chukji = calc_chukji(real_final_contig)
        real_final_contig = real_final_contig[:-1]
        real_final_ctg_typ = final_ctg_typ | alt_final_ctg_typ
        # with open(f"06_repeat_preprocess/{PAF_FILE_PATH[0].split("/")[-1].split(".")[0]}_real_final_ratio.txt", "wt") as f:
        #     a =[]
        #     for i in real_final_ratio:
        #         if real_final_ctg_typ[i] in (3, 4):
        #             a.append((i, real_final_ratio[i], chukji[i]))
        #     a = sorted(a, key=lambda x:x[1])
        #     for i in a:
        #         print(i[0], i[1], int(i[2]/1000), file=f)

        final_len_count = Counter()
        s = 0
        while s<len(real_final_contig):
            e = real_final_contig[s][CTG_ENDND]
            for i in range(s, e+1):
                final_len_count[real_final_contig[i][CTG_NAM]] += real_final_contig[i][CHR_END] - real_final_contig[i][CHR_STR]
            s = e+1
        a = len_counter + alt_len_counter
        ppc_bna_ratio = []
        
        # prefix = os.path.basename(PREPROCESSED_PAF_FILE_PATH).split('.')[0]
        # with open(f"06_repeat_preprocess/repeat_ppc_before_and_after.{prefix}.txt", "wt") as f:
        #     for i in a:
        #         ppc_bna_ratio.append((i, final_len_count[i]/a[i]))
        #     ppc_bna_ratio.sort(key=lambda x:x[1])
        #     for i in ppc_bna_ratio:
        #         if 0 < i[1] < 1:
        #             print(i[0],i[1], file=f)

if __name__ == "__main__":
    main()