import re
import argparse
import pandas as pd
from collections import defaultdict

CTG_NAM = 0
CTG_LEN = 1
CTG_STR = 2
CTG_END = 3
CTG_DIR = 4
CHR_NAM = 5
CHR_LEN = 6
CHR_STR = 7
CHR_END = 8
CTG_TYP = 9
CTG_STRND = 10
CTG_ENDND = 11
CTG_TELCHR = 12
CTG_TELDIR = 13
DIR_FOR = 1
DIR_BAK = 0
INF = 1000000000
BUFFER = 10000000
CHROMOSOME_COUNT = 23
K = 1000

def dbg():
    print("hi")

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

def import_data(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.rstrip().split("\t")
        int_induce_idx = [1, 2, 3, 6, 7, 8, 9, 10, 11,]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    return contig_data


def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    return chr_len

def distance_checker(node_a : tuple, node_b : tuple) -> int :
    # check if overlap
    if max(node_a[CHR_STR], node_b[CHR_STR]) < min(node_a[CHR_END], node_b[CHR_END]):
        return 0   
    else:
        # return distance of two ranges
        return min(abs(node_b[CHR_STR] - node_a[CHR_END]), abs(node_b[CHR_END] - node_a[CHR_STR]))

def telo_distance_checker(node: tuple, telo: tuple) -> int :
    if node[CTG_DIR] == '+':
        return min(abs(telo[CHR_STR] - node[CHR_END]), abs(telo[CHR_END] - node[CHR_STR]))
    else:
        return min(abs(telo[CHR_STR] - node[CHR_STR]), abs(telo[CHR_END] - node[CHR_END]))
    
def inclusive_checker(contig_node : tuple, telomere_node : tuple) -> bool :
    if int(telomere_node[CHR_STR]) <= int(contig_node[CHR_STR]) and int(contig_node[CHR_END]) <= int(telomere_node[CHR_END]):
        return True
    else:
        return False

def initial_graph_build(contig_data : list, telo_data : dict) -> list :
    '''
    Initialize
    '''
    contig_data_size = len(contig_data)
    report_case = {'A':[], 'B':[], 'C':[], 'ESC':[], 'ALL_TELO_NON_ESC':[]}
    chr_corr = {}
    chr_rev_corr = {}
    total_contig_count = len(contig_data)
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'f'] = total_contig_count + i - 1
        chr_rev_corr[total_contig_count + i - 1] = 'chr'+str(i)+'f'
    chr_corr['chrXf'] = total_contig_count + CHROMOSOME_COUNT - 1
    chr_corr['chrYf'] = total_contig_count + CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + CHROMOSOME_COUNT - 1] = 'chrXf'
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'b'] = total_contig_count + CHROMOSOME_COUNT + i - 1
        chr_rev_corr[total_contig_count + CHROMOSOME_COUNT + i - 1] = 'chr'+str(i)+'b'
    chr_corr['chrXb'] = total_contig_count + 2*CHROMOSOME_COUNT - 1
    chr_corr['chrYb'] = total_contig_count + 2*CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + 2*CHROMOSOME_COUNT - 1] = 'chrXb'

    adjacency = [[[] for _ in range(total_contig_count+CHROMOSOME_COUNT*2)], [[] for _ in range(total_contig_count+CHROMOSOME_COUNT*2)]]
    curr_contig_name = 0
    last_node = 0
    telo_detect = False
    st = 0
    ed = INF
    
    '''
    Algorithm
    '''
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, '0', '0', '0', '0', ))
    curr_contig_st = 0
    curr_contig_ed = contig_data[curr_contig_st][CTG_ENDND]
    while curr_contig_st<total_contig_count:
        # front of contig
        st = curr_contig_st
        ed = curr_contig_ed+1
        front_telo_bound = curr_contig_st
        end_telo_bound = curr_contig_ed
        while contig_data[front_telo_bound][CTG_TELCHR] != '0' and front_telo_bound<=curr_contig_ed:
            front_telo_bound+=1
        front_telo_bound-=1
        while contig_data[end_telo_bound][CTG_TELCHR] != '0' and end_telo_bound>=curr_contig_st:
            end_telo_bound-=1
        end_telo_bound+=1
        escape = False
        if front_telo_bound == curr_contig_ed:
            if len(contig_data[curr_contig_st][CTG_TELDIR])==1 \
            and (contig_data[curr_contig_st][CTG_TELDIR] + contig_data[curr_contig_st][CTG_DIR] in ("b+", "f-")):
                escape = True
            elif len(contig_data[curr_contig_ed][CTG_TELDIR])==1 \
            and (contig_data[curr_contig_ed][CTG_TELDIR] + contig_data[curr_contig_ed][CTG_DIR] in ("b-", "f+")):
                escape = True
        else:
            escape = True
        if escape:
            if front_telo_bound > curr_contig_st:
                # Check if first telomere node is 'Nin' -> Else: escape
                if len(contig_data[curr_contig_st][CTG_TELDIR])==1 \
                and (contig_data[curr_contig_st][CTG_TELDIR] + contig_data[curr_contig_st][CTG_DIR] in ("b+", "f-")):
                    report_case['ESC'].append(curr_contig_st)
                    st = curr_contig_st
                # If boundary node is 'Nin'
                elif len(contig_data[front_telo_bound][CTG_TELDIR])>1:
                    # Next node is connected with telomere node.
                    front_telo_bound+=1
                    dest = contig_data[front_telo_bound][CHR_NAM]
                    if contig_data[front_telo_bound][CTG_DIR] == '+':
                        dest += 'f'
                        adjacency[DIR_FOR][chr_corr[dest]].append([DIR_FOR, front_telo_bound, 0])
                        adjacency[DIR_FOR][front_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                        adjacency[DIR_BAK][chr_corr[dest]].append([DIR_BAK, front_telo_bound, 0])
                        adjacency[DIR_BAK][front_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                        report_case['A'].append([dest, front_telo_bound])
                    else:
                        dest += 'b'
                        adjacency[DIR_FOR][chr_corr[dest]].append([DIR_BAK, front_telo_bound, 0])
                        adjacency[DIR_FOR][front_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                        adjacency[DIR_BAK][chr_corr[dest]].append([DIR_FOR, front_telo_bound, 0])
                        adjacency[DIR_BAK][front_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                        report_case['A'].append([dest, front_telo_bound])
                    st = front_telo_bound
                # If boundary node is not "Nin"
                else:
                    if contig_data[front_telo_bound][CTG_TELDIR] == 'f':
                        if contig_data[front_telo_bound][CTG_DIR]=='+':
                            # boundary node is connected with telomere node.
                            dest = contig_data[front_telo_bound][CHR_NAM]+'f'
                            adjacency[DIR_FOR][chr_corr[dest]].append([DIR_FOR, front_telo_bound, 0])
                            adjacency[DIR_FOR][front_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                            adjacency[DIR_BAK][chr_corr[dest]].append([DIR_BAK, front_telo_bound, 0])
                            adjacency[DIR_BAK][front_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                            report_case['B'].append([dest, front_telo_bound])
                        else:
                            # Treat as "Nin" Case
                            front_telo_bound+=1
                            dest = contig_data[front_telo_bound][CHR_NAM]
                            if contig_data[front_telo_bound][CTG_DIR] == '+':
                                dest += 'f'
                                adjacency[DIR_FOR][chr_corr[dest]].append([DIR_FOR, front_telo_bound, 0])
                                adjacency[DIR_FOR][front_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                                adjacency[DIR_BAK][chr_corr[dest]].append([DIR_BAK, front_telo_bound, 0])
                                adjacency[DIR_BAK][front_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                                report_case['C'].append([dest, front_telo_bound])
                            else:
                                dest += 'b'
                                adjacency[DIR_FOR][chr_corr[dest]].append([DIR_BAK, front_telo_bound, 0])
                                adjacency[DIR_FOR][front_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                                adjacency[DIR_BAK][chr_corr[dest]].append([DIR_FOR, front_telo_bound, 0])
                                adjacency[DIR_BAK][front_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                                report_case['C'].append([dest, front_telo_bound])
                    else:
                        if contig_data[front_telo_bound][CTG_DIR]=='-':
                            dest = contig_data[front_telo_bound][CHR_NAM]+'b'
                            adjacency[DIR_FOR][chr_corr[dest]].append([DIR_BAK, front_telo_bound, 0])
                            adjacency[DIR_FOR][front_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                            adjacency[DIR_BAK][chr_corr[dest]].append([DIR_FOR, front_telo_bound, 0])
                            adjacency[DIR_BAK][front_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                            report_case['B'].append([dest, front_telo_bound])
                        else:
                            front_telo_bound+=1
                            dest = contig_data[front_telo_bound][CHR_NAM]
                            if contig_data[front_telo_bound][CTG_DIR] == '+':
                                dest += 'f'
                                adjacency[DIR_FOR][chr_corr[dest]].append([DIR_FOR, front_telo_bound, 0])
                                adjacency[DIR_FOR][front_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                                adjacency[DIR_BAK][chr_corr[dest]].append([DIR_BAK, front_telo_bound, 0])
                                adjacency[DIR_BAK][front_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                                report_case['C'].append([dest, front_telo_bound])
                            else:
                                dest += 'b'
                                adjacency[DIR_FOR][chr_corr[dest]].append([DIR_BAK, front_telo_bound, 0])
                                adjacency[DIR_FOR][front_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                                adjacency[DIR_BAK][chr_corr[dest]].append([DIR_FOR, front_telo_bound, 0])
                                adjacency[DIR_BAK][front_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                                report_case['C'].append([dest, front_telo_bound])
                    st = front_telo_bound

            if end_telo_bound < curr_contig_ed:
                # Check if first telomere node is 'Nin' -> Else: escape
                if len(contig_data[curr_contig_ed][CTG_TELDIR])==1 \
                and (contig_data[curr_contig_ed][CTG_TELDIR] + contig_data[curr_contig_ed][CTG_DIR] in ("b-", "f+")):
                    report_case['ESC'].append(curr_contig_ed)
                    ed = curr_contig_ed+1
                # If boundary node is 'Nin'
                elif len(contig_data[end_telo_bound][CTG_TELDIR])>1:
                    # Next node is connected with telomere node.
                    end_telo_bound-=1
                    dest = contig_data[end_telo_bound][CHR_NAM]
                    if contig_data[end_telo_bound][CTG_DIR] == '+':
                        dest += 'b'
                        adjacency[DIR_FOR][chr_corr[dest]].append([DIR_FOR, end_telo_bound, 0])
                        adjacency[DIR_FOR][end_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                        adjacency[DIR_BAK][chr_corr[dest]].append([DIR_BAK, end_telo_bound, 0])
                        adjacency[DIR_BAK][end_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                        report_case['A'].append([dest, end_telo_bound])
                    else:
                        dest += 'f'
                        adjacency[DIR_FOR][chr_corr[dest]].append([DIR_BAK, end_telo_bound, 0])
                        adjacency[DIR_FOR][end_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                        adjacency[DIR_BAK][chr_corr[dest]].append([DIR_FOR, end_telo_bound, 0])
                        adjacency[DIR_BAK][end_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                        report_case['A'].append([dest, end_telo_bound])
                    ed = end_telo_bound
                # If boundary node is not "Nin"
                else:
                    if contig_data[end_telo_bound][CTG_TELDIR] == 'b':
                        if contig_data[end_telo_bound][CTG_DIR]=='+':
                            # boundary node is connected with telomere node.
                            dest = contig_data[end_telo_bound][CHR_NAM]+'b'
                            adjacency[DIR_FOR][chr_corr[dest]].append([DIR_FOR, end_telo_bound, 0])
                            adjacency[DIR_FOR][end_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                            adjacency[DIR_BAK][chr_corr[dest]].append([DIR_BAK, end_telo_bound, 0])
                            adjacency[DIR_BAK][end_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                            report_case['B'].append([dest, end_telo_bound])
                        else:
                            # Treat as "Nin" Case
                            end_telo_bound+=1
                            dest = contig_data[end_telo_bound][CHR_NAM]
                            if contig_data[end_telo_bound][CTG_DIR] == '+':
                                dest += 'b'
                                adjacency[DIR_FOR][chr_corr[dest]].append([DIR_FOR, end_telo_bound, 0])
                                adjacency[DIR_FOR][end_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                                adjacency[DIR_BAK][chr_corr[dest]].append([DIR_BAK, end_telo_bound, 0])
                                adjacency[DIR_BAK][end_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                                report_case['C'].append([dest, end_telo_bound])
                            else:
                                dest += 'f'
                                adjacency[DIR_FOR][chr_corr[dest]].append([DIR_BAK, end_telo_bound, 0])
                                adjacency[DIR_FOR][end_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                                adjacency[DIR_BAK][chr_corr[dest]].append([DIR_FOR, end_telo_bound, 0])
                                adjacency[DIR_BAK][end_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                                report_case['C'].append([dest, end_telo_bound])
                    else:
                        if contig_data[end_telo_bound][CTG_DIR]=='-':
                            dest = contig_data[end_telo_bound][CHR_NAM]+'f'
                            adjacency[DIR_FOR][chr_corr[dest]].append([DIR_BAK, end_telo_bound, 0])
                            adjacency[DIR_FOR][end_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                            adjacency[DIR_BAK][chr_corr[dest]].append([DIR_FOR, end_telo_bound, 0])
                            adjacency[DIR_BAK][end_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                            report_case['B'].append([dest, end_telo_bound])
                        else:
                            end_telo_bound+=1
                            dest = contig_data[end_telo_bound][CHR_NAM]
                            if contig_data[end_telo_bound][CTG_DIR] == '+':
                                dest += 'b'
                                adjacency[DIR_FOR][chr_corr[dest]].append([DIR_FOR, end_telo_bound, 0])
                                adjacency[DIR_FOR][end_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                                adjacency[DIR_BAK][chr_corr[dest]].append([DIR_BAK, end_telo_bound, 0])
                                adjacency[DIR_BAK][end_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                                report_case['C'].append([dest, end_telo_bound])
                            else:
                                dest += 'f'
                                adjacency[DIR_FOR][chr_corr[dest]].append([DIR_BAK, end_telo_bound, 0])
                                adjacency[DIR_FOR][end_telo_bound].append([DIR_BAK, chr_corr[dest], 0])
                                adjacency[DIR_BAK][chr_corr[dest]].append([DIR_FOR, end_telo_bound, 0])
                                adjacency[DIR_BAK][end_telo_bound].append([DIR_FOR, chr_corr[dest], 0])
                                report_case['C'].append([dest, end_telo_bound])
                    ed = end_telo_bound
            for j in range(st+1, ed):
                adjacency[DIR_FOR][j-1].append([DIR_FOR, j, 0])
                adjacency[DIR_BAK][j].append([DIR_BAK, j-1, 0])
        else:
            report_case['ALL_TELO_NON_ESC'].append(curr_contig_st)
        curr_contig_st = curr_contig_ed + 1
        curr_contig_ed = contig_data[curr_contig_st][CTG_ENDND]
    '''
    If telomere node has 0 connection, connect with closest node with same chromosome type.
    '''
    end_node_list = set()
    contig_data = contig_data[:-1]
    for i in contig_data:
        end_node_list.add((i[CTG_STRND], i[CTG_ENDND])) 
    for i in range(contig_data_size, contig_data_size + CHROMOSOME_COUNT*2):
        if True:
            telo_dist = INF
            telo_connect_node = 0
            telo_dir = 0
            telo_sign = 0
            now_telo = chr_rev_corr[i]
            temp_contig = (0, 0, 0, 0, 0, 0, 0, telo_data[now_telo][0], telo_data[now_telo][1])
            for st, ed in end_node_list:
                if contig_data[st][CHR_NAM] == now_telo[0:-1] \
                    and telo_distance_checker(contig_data[st], temp_contig) < telo_dist\
                    and contig_data[st][CTG_TELCHR] == '0':
                    telo_dist = telo_distance_checker(contig_data[st], temp_contig)
                    telo_dir = 'f'
                    telo_sign = contig_data[st][CTG_DIR]
                    telo_connect_node = st
                if contig_data[ed][CHR_NAM] == now_telo[0:-1] \
                    and telo_distance_checker(contig_data[ed], temp_contig) < telo_dist\
                    and contig_data[ed][CTG_TELCHR] == '0':
                    telo_dist = telo_distance_checker(contig_data[ed], temp_contig)
                    telo_connect_node = ed
                    telo_sign = contig_data[st][CTG_DIR]
                    telo_dir = 'b'
            # 고칠 점: 최종 값만 하지 말고 모든 값에 대해 해봐야 됨
            if now_telo[-1]=='f':
                if telo_dir+telo_sign in {"f+", "b-"}:
                    for sign in (DIR_FOR, DIR_BAK):
                        for j in (telo_connect_node, i):
                            chk = 1
                            for k in adjacency[sign][j]:
                                if k[1] == telo_connect_node+i-j:
                                    chk = 0
                            if chk:
                                adjacency[sign][j].append([sign, telo_connect_node+i-j, telo_dist])
                else:
                    for sign in (DIR_FOR, DIR_BAK):
                        for j in (telo_connect_node, i):
                            chk = 1
                            for k in adjacency[sign][j]:
                                if k[1] == telo_connect_node+i-j:
                                    chk = 0
                            if chk:
                                adjacency[sign][j].append([1 - sign, telo_connect_node+i-j, telo_dist])
            else:
                if telo_dir+telo_sign in {"b+", "f-"}:
                    for sign in (DIR_FOR, DIR_BAK):
                        for j in (telo_connect_node, i):
                            chk = 1
                            for k in adjacency[sign][j]:
                                if k[1] == telo_connect_node+i-j:
                                    chk = 0
                            if chk:
                                adjacency[sign][j].append([sign, telo_connect_node+i-j, telo_dist])
                else:
                    for sign in (DIR_FOR, DIR_BAK):
                        for j in (telo_connect_node, i):
                            chk = 1
                            for k in adjacency[sign][j]:
                                if k[1] == telo_connect_node+i-j:
                                    chk = 0
                            if chk:
                                adjacency[sign][j].append([1 - sign, telo_connect_node+i-j, telo_dist])

    return [adjacency, report_case]

def node_label(contig_data : list) -> list :
    total_contig_count = len(contig_data)
    label = [0 for _ in range(total_contig_count)]
    for i in range(0, total_contig_count) :
        contig_start = contig_data[contig_data[i][CTG_STRND]]
        contig_end = contig_data[contig_data[i][CTG_ENDND]]

        # Contig's start and end node is labeled as "linkable"

        if contig_data[i] == contig_start or contig_data[i] == contig_end:
            label[i] = "4"

        # TYPE 3 : START AND END HAS SAME CHROMOSOME TYPE & SAME DIRECTION
        # Condition: If node is formed with same type of chromosome & type = 3
        elif contig_data[i][CTG_TYP] == 3 \
             and contig_data[i][CHR_NAM] == contig_start[CHR_NAM] \
             and contig_data[i][CTG_DIR] == contig_start[CTG_DIR]:
            if contig_data[i][CTG_DIR] == '+':
                if contig_start[CHR_END] <= contig_data[i][CHR_STR] \
                and contig_data[i][CHR_END] <= contig_end[CHR_STR] :
                    label[i] = "1-1"
            else:
                if contig_start[CHR_STR] >= contig_data[i][CHR_END] \
                and contig_data[i][CHR_STR] >= contig_end[CHR_END] :
                    label[i] = "1-2"

        # TYPE 2 : START AND END HAS SAME CHROMOSOME TYPE & OPPOSITE DIRECTION
        elif contig_data[i][CTG_TYP] == 2 and contig_data[i][CHR_NAM] == contig_start[CHR_NAM]:
            template_node = contig_start
            curr_dir = contig_data[i][CTG_DIR]
            # CASE 1 : Current node has same direction with end node : start has to be predicted
            if(contig_data[i][CTG_DIR] != contig_start[CTG_DIR]):
                template_node = contig_end
                # CASE 1-1 : Current node's direction is positive
                if curr_dir == '+' : 
                    predicted_start = contig_end[CHR_END] \
                                    - (contig_data[i][CHR_END] - contig_data[i][CHR_STR]) \
                                    - BUFFER
                    if predicted_start <= contig_data[i][CHR_STR] \
                    and contig_data[i][CHR_END] <= contig_end[CHR_STR] :
                        label[i] = "2-1-1"
                # CASE 1-2 : Current node's direction is negative
                # In this case, curr.ed > curr.st > ctg_end.ed > ctg_end.st
                # So, curr.st > ctg_end.ed, curr.ed < predicted_start, which derived by ctg_end.st
                elif curr_dir == '-':
                    predicted_start = contig_end[CHR_STR] \
                                    + (contig_data[i][CHR_END] - contig_data[i][CHR_STR]) \
                                    + BUFFER
                    if predicted_start >= contig_data[i][CHR_END] \
                       and contig_data[i][CHR_STR] >= contig_end[CHR_END] : 
                        label[i] = "2-1-2"
            # CASE 2 : Current node has same direction with start node : end has to be predicted
            else:
                # CASE 2-1 : Current node's direction is positive
                if curr_dir == '+' : 
                    predicted_end = contig_start[CHR_STR] \
                                    + (contig_data[i][CHR_END] - contig_data[i][CHR_STR]) \
                                    + BUFFER
                    if contig_start[CHR_END] <= contig_data[i][CHR_STR] \
                    and contig_data[i][CHR_END] <= predicted_end :
                        label[i] = "2-2-1"
                # CASE 2-2 : Current node's direction is negative
                # In this case, ctg_str.ed > ctg_str.st > curr.ed > curr.st > pred_end
                # So, ctg_str.st > curr.ed && curr.st > pred_end, pred_end is derived from ctg_str.ed
                elif curr_dir == '-':
                    predicted_end = contig_start[CHR_END] \
                                    - (contig_data[i][CHR_END] - contig_data[i][CHR_STR]) \
                                    - BUFFER
                    if contig_start[CHR_STR] >= contig_data[i][CHR_END] \
                    and contig_data[i][CHR_STR] >= predicted_end : 
                        label[i] = "2-2-2"
        
        # TYPE 1 : START AND END HAS DIFFERENT CHROMOSOME TYPE
        elif contig_data[i][CTG_TYP] == 1:
            curr_dir = contig_data[i][CTG_DIR]
            # To be labeled as linkable, its chromosome type 
            # AND direction has to be same with terminal nodes
            # Case 1 : End chromosome has same type & direction : Start has to be predicted
            if contig_data[i][CHR_NAM] == contig_end[CHR_NAM] \
            and curr_dir == contig_end[CTG_DIR] :
                if curr_dir == '+':
                    predicted_start = contig_end[CHR_END] \
                                    - (contig_data[i][CHR_END] - contig_data[i][CHR_STR]) \
                                    - BUFFER
                    if predicted_start <= contig_data[i][CHR_STR] \
                    and contig_data[i][CHR_END] <= contig_end[CHR_STR] :
                        label[i] = "3-1-1"
                elif curr_dir == '-':
                    predicted_start = contig_end[CHR_STR] \
                                    + (contig_data[i][CHR_END] - contig_data[i][CHR_STR]) \
                                    + BUFFER
                    if predicted_start >= contig_data[i][CHR_END] \
                       and contig_data[i][CHR_STR] >= contig_end[CHR_END] : 
                        label[i] = "3-1-2"
            # CASE 2 : Current node has same direction with start node : end has to be predicted
            elif contig_data[i][CHR_NAM] == contig_start[CHR_NAM] \
            and curr_dir == contig_start[CTG_DIR] :
                # CASE 2-1 : Current node's direction is positive
                if curr_dir == '+' : 
                    predicted_end = contig_start[CHR_STR] \
                                    + (contig_data[i][CHR_END] - contig_data[i][CHR_STR]) \
                                    + BUFFER
                    if contig_start[CHR_END] <= contig_data[i][CHR_STR] \
                    and contig_data[i][CHR_END] <= predicted_end :
                        label[i] = "3-2-1"
                # CASE 2-2 : Current node's direction is negative
                # In this case, ctg_str.ed > ctg_str.st > curr.ed > curr.st > pred_end
                # So, ctg_str.st > curr.ed && curr.st > pred_end, pred_end is derived from ctg_str.ed
                elif curr_dir == '-':
                    predicted_end = contig_start[CHR_END] \
                                    - (contig_data[i][CHR_END] - contig_data[i][CHR_STR]) \
                                    - BUFFER
                    if contig_start[CHR_STR] >= contig_data[i][CHR_END] \
                    and contig_data[i][CHR_STR] >= predicted_end : 
                        label[i] = "3-2-2"
                
    return label


def graph_build(contig_data : list, contig_adjacency : list, linkability_label : list) -> list :
    total_contig_count = len(contig_data)
    gap_contig_adjacency = [[[[], []] for __ in range(total_contig_count)] for _ in range(2)]
    for i in range(total_contig_count):
        if linkability_label[i] != 0:
            k = contig_data[i][CTG_ENDND]+1
            for j in range(k, total_contig_count):
                if linkability_label[j] != 0 \
                and contig_data[i][CHR_NAM] == contig_data[j][CHR_NAM] \
                and contig_data[i][CTG_NAM] != contig_data[j][CTG_NAM] :
                    if (inclusive_checker(contig_data[i], contig_data[j]) or inclusive_checker(contig_data[j], contig_data[i]))\
                    and (contig_data[i][CTG_ENDND] - contig_data[i][CTG_STRND] == 0 or contig_data[j][CTG_ENDND] - contig_data[j][CTG_STRND] == 0):
                        continue
                    dist = distance_checker(contig_data[i], contig_data[j])
                    # CASE 1 : Both contigs have same direction
                    if contig_data[i][CTG_DIR] == contig_data[j][CTG_DIR]:
                        # 일단 dist가 0일 때부터.
                        if dist == 0:
                            contig_adjacency[DIR_FOR][i].append([DIR_FOR, j, dist])
                            contig_adjacency[DIR_FOR][j].append([DIR_FOR, i, dist])
                            contig_adjacency[DIR_BAK][i].append([DIR_BAK, j, dist])
                            contig_adjacency[DIR_BAK][j].append([DIR_BAK, i, dist])
                        else:
                            # i에서 j로 증가하는 경우.
                            if contig_data[i][CHR_END] < contig_data[j][CHR_STR]:
                                if len(gap_contig_adjacency[DIR_FOR][i][0])==0:
                                    gap_contig_adjacency[DIR_FOR][i][0] = [DIR_FOR, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][i][0] = [DIR_FOR, j, dist]
                                
                                if len(gap_contig_adjacency[DIR_BAK][i][0])==0:
                                    gap_contig_adjacency[DIR_BAK][i][0] = [DIR_BAK, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][i][0] = [DIR_BAK, j, dist]
                                
                                if len(gap_contig_adjacency[DIR_FOR][j][1])==0:
                                    gap_contig_adjacency[DIR_FOR][j][1] = [DIR_FOR, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][j][1][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][j][1] = [DIR_FOR, i, dist]
                                
                                if len(gap_contig_adjacency[DIR_BAK][j][1])==0:
                                    gap_contig_adjacency[DIR_BAK][j][1] = [DIR_BAK, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][j][1][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][j][1] = [DIR_BAK, i, dist]
                            # j에서 i로 증가하는 경우.
                            else:
                                if len(gap_contig_adjacency[DIR_FOR][i][1])==0:
                                    gap_contig_adjacency[DIR_FOR][i][1] = [DIR_FOR, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][i][1][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][i][1] = [DIR_FOR, j, dist]
                                
                                if len(gap_contig_adjacency[DIR_BAK][i][1])==0:
                                    gap_contig_adjacency[DIR_BAK][i][1] = [DIR_BAK, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][i][1][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][i][1] = [DIR_BAK, j, dist]
                                
                                if len(gap_contig_adjacency[DIR_FOR][j][0])==0:
                                    gap_contig_adjacency[DIR_FOR][j][0] = [DIR_FOR, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][j][0] = [DIR_FOR, i, dist]
                                
                                if len(gap_contig_adjacency[DIR_BAK][j][0])==0:
                                    gap_contig_adjacency[DIR_BAK][j][0] = [DIR_BAK, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][j][0] = [DIR_BAK, i, dist]
                    # Case two: two contigs have opposite direction
                    else:
                        if dist == 0:
                            contig_adjacency[DIR_FOR][i].append([DIR_BAK, j, dist])
                            contig_adjacency[DIR_FOR][j].append([DIR_BAK, i, dist])
                            contig_adjacency[DIR_BAK][i].append([DIR_FOR, j, dist])
                            contig_adjacency[DIR_BAK][j].append([DIR_FOR, i, dist])
                        else:
                            # i에서 j로 증가하는 경우.
                            if contig_data[i][CHR_END] < contig_data[j][CHR_STR]:
                                if len(gap_contig_adjacency[DIR_FOR][i][0])==0:
                                    gap_contig_adjacency[DIR_FOR][i][0] = [DIR_BAK, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][i][0] = [DIR_BAK, j, dist]
                                
                                if len(gap_contig_adjacency[DIR_BAK][i][0])==0:
                                    gap_contig_adjacency[DIR_BAK][i][0] = [DIR_FOR, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][i][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][i][0] = [DIR_FOR, j, dist]
                                
                                if len(gap_contig_adjacency[DIR_FOR][j][1])==0:
                                    gap_contig_adjacency[DIR_FOR][j][1] = [DIR_BAK, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][j][1][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][j][1] = [DIR_BAK, i, dist]
                                
                                if len(gap_contig_adjacency[DIR_BAK][j][1])==0:
                                    gap_contig_adjacency[DIR_BAK][j][1] = [DIR_FOR, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][j][1][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][j][1] = [DIR_FOR, i, dist]
                            # j에서 i로 증가하는 경우.
                            else:
                                if len(gap_contig_adjacency[DIR_FOR][i][1])==0:
                                    gap_contig_adjacency[DIR_FOR][i][1] = [DIR_BAK, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][i][1][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][i][1] = [DIR_BAK, j, dist]
                                
                                if len(gap_contig_adjacency[DIR_BAK][i][1])==0:
                                    gap_contig_adjacency[DIR_BAK][i][1] = [DIR_FOR, j, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][i][1][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][i][1] = [DIR_FOR, j, dist]
                                
                                if len(gap_contig_adjacency[DIR_FOR][j][0])==0:
                                    gap_contig_adjacency[DIR_FOR][j][0] = [DIR_BAK, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_FOR][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_FOR][j][0] = [DIR_BAK, i, dist]
                                
                                if len(gap_contig_adjacency[DIR_BAK][j][0])==0:
                                    gap_contig_adjacency[DIR_BAK][j][0] = [DIR_FOR, i, dist]
                                else:
                                    if gap_contig_adjacency[DIR_BAK][j][0][2]>dist:
                                        gap_contig_adjacency[DIR_BAK][j][0] = [DIR_FOR, i, dist]
    for i in range(2):
        for j in range(total_contig_count):
            for k in range(2):
                if len(gap_contig_adjacency[i][j][k])>0:
                    contig_adjacency[i][j].append(gap_contig_adjacency[i][j][k])
    return contig_adjacency
                        
def edge_optimization(contig_data : list, contig_adjacency : list, telo_dict : dict) -> list :
    chr_corr = {}
    chr_rev_corr = {}
    total_contig_count = len(contig_data)
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'f'] = total_contig_count + i - 1
        chr_rev_corr[total_contig_count + i - 1] = 'chr'+str(i)+'f'
    chr_corr['chrXf'] = total_contig_count + CHROMOSOME_COUNT - 1
    chr_corr['chrYf'] = total_contig_count + CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + CHROMOSOME_COUNT - 1] = 'chrXf'
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'b'] = total_contig_count + CHROMOSOME_COUNT + i - 1
        chr_rev_corr[total_contig_count + CHROMOSOME_COUNT + i - 1] = 'chr'+str(i)+'b'
    chr_corr['chrXb'] = total_contig_count + 2*CHROMOSOME_COUNT - 1
    chr_corr['chrYb'] = total_contig_count + 2*CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + 2*CHROMOSOME_COUNT - 1] = 'chrXb'
    contig_pair_nodes = defaultdict(list)
    optimized_adjacency = [[[] for _ in range(total_contig_count + CHROMOSOME_COUNT*2)], [[] for _ in range(total_contig_count + CHROMOSOME_COUNT*2)]]
    for _ in range(2):
        for i in range(len(contig_data)+CHROMOSOME_COUNT*2):
            for edge in contig_adjacency[_][i]:
                if edge[1] >= total_contig_count or i >= total_contig_count:
                    optimized_adjacency[_][i].append(edge)
                    continue
                a = (ctg_name_to_int(contig_data[i][CTG_NAM]))
                b = (ctg_name_to_int(contig_data[edge[1]][CTG_NAM]))
                if int(a)>int(b):
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
    total_contig_count = len(contig_data)
    for i in range(2):
        for j in range(total_contig_count):
            first_contig_name = contig_data[j][CTG_NAM]
            fcn_int = int(ctg_name_to_int(first_contig_name))
            for edge in contig_adjacency[i][j]:
                if edge[1] >= total_contig_count:
                    continue
                second_contig_name = contig_data[edge[1]][CTG_NAM]
                scn_int = int(ctg_name_to_int(second_contig_name))
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
    for i in range(total_contig_count, total_contig_count+2*CHROMOSOME_COUNT):
        telo_name = chr_rev_corr[i]
        telo_range = telo_dict[telo_name]
        for j in range(2):
            using_edge = []
            now_edge = [-1, 0, 0]
            for edge in optimized_adjacency[j][i]:
                if telo_name[-1]=='f':
                    if contig_data[edge[1]][CHR_END]<=10*K:
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
                    if contig_data[edge[1]][CHR_STR]>=telo_range[1]-10*K:
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


    return optimized_adjacency

def ctg_name_to_int(ctg_name):
    return re.search(r'\d+$', ctg_name[:-1]).group()

# Test area
    
def main():
    parser = argparse.ArgumentParser(description="Build graph using preprocessed paf")
    
    # 위치 인자 정의
    parser.add_argument("ppc_paf_file_path", 
                        help="Path to the preprocessed PAF file.")
    parser.add_argument("telomere_bed_path", 
                        help="Path to the telomere information file.")
    parser.add_argument("reference_fai_path", 
                        help="Path to the chromosome information file.")

    args = parser.parse_args()

    # All preprocessed contigs are numbered from 0 to contig_count-1

    PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
    TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
    CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path

    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
    telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)

    telo_dict = {}
    for i in telo_data:
        telo_dict[i[0]+i[-1]] = i[1:3]

    contig_data = import_data(PREPROCESSED_PAF_FILE_PATH)


    contig_adjacency, report_case = initial_graph_build(contig_data, telo_dict)

    for _ in report_case:
        print(_, report_case[_])

    contig_data = contig_data[0:-1]

    label = node_label(contig_data)

    a = pd.Series(label)
    print(a.value_counts())

    contig_adjacency = graph_build(contig_data, contig_adjacency, label)

    optimized_adjacency = edge_optimization(contig_data, contig_adjacency, telo_dict)

    chr_corr = {}
    chr_rev_corr = {}
    total_contig_count = len(contig_data)
    print(total_contig_count)
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'f'] = total_contig_count + i - 1
        chr_rev_corr[total_contig_count + i - 1] = 'chr'+str(i)+'f'
    chr_corr['chrXf'] = total_contig_count + CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + CHROMOSOME_COUNT - 1] = 'chrXf'
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'b'] = total_contig_count + CHROMOSOME_COUNT + i - 1
        chr_rev_corr[total_contig_count + CHROMOSOME_COUNT + i - 1] = 'chr'+str(i)+'b'
    chr_corr['chrXb'] = total_contig_count + 2*CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + 2*CHROMOSOME_COUNT - 1] = 'chrXb'

    cnt = 0
    print("Number of Nodes:", len(contig_data))
    with open(PREPROCESSED_PAF_FILE_PATH + ".np.graph.txt", "wt") as f:
        for i in range(2):
            for j in range(len(contig_data)+CHROMOSOME_COUNT*2):
                f.write(f"({i}, {j}) : ")
                for k in contig_adjacency[i][j]:
                    f.write(f"{k}, ")
                    cnt+=1
                f.write("\n")
    print("Number of Edges:", cnt)

    cnt = 0
    with open(PREPROCESSED_PAF_FILE_PATH + ".op.graph.txt", "wt") as g:
        for i in range(2):
            for j in range(len(contig_data)+CHROMOSOME_COUNT*2):
                g.write(f"({i}, {j}) : ")
                for k in optimized_adjacency[i][j]:
                    g.write(f"{k}, ")
                    cnt+=1
                g.write("\n")
    
    with open(PREPROCESSED_PAF_FILE_PATH + ".telomere_connected_node.txt", "wt") as h:
        for i in range(2):
            for j in range(len(contig_data), len(contig_data)+CHROMOSOME_COUNT*2):
                h.write(f"({i}, {chr_rev_corr[j]}) : \n")
                for k in optimized_adjacency[i][j]:
                    h.write(f"CONNECTED_NODE_DIR: {k[0]}, DIST: {k[2]}, NODE_NUM : {k[1]}, NODE_DATA: {contig_data[k[1]]}\n")


if __name__ == "__main__":
    main()