import re
import sys
import argparse
import os
from collections import defaultdict
import ast
import glob
import networkx as nx

from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("10_Breakend_Graph_Build_Paths start")

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
CTG_MAINFLOWDIR = 18
CTG_MAINFLOWCHR = 19


DIR_FOR = 1
DIR_BAK = 0
DIR_IN = 3
DIR_OUT = 2
INF = 1000000000
BUFFER = 10000000
CHROMOSOME_COUNT = 23
K = 1000

NCLOSE_COMPRESS_LIMIT = 50*K

def import_data(file_path : str) -> dict :
    graph_file = open(file_path, "r")
    graph_adjacency = {}
    cnt = 0
    for curr_edge in graph_file:
        l, r = curr_edge.split(":")
        r.lstrip()
        r.rstrip(',')
        r = ast.literal_eval('['+r+']')
        if cnt==0:
            cnt+=1
        l = ast.literal_eval(l)
        graph_adjacency[l] = r
    graph_file.close()
    return graph_adjacency

def import_data2(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.split("\t")
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END, \
                          CHR_LEN, CHR_STR, CHR_END, \
                          CTG_MAPQ, CTG_TYP, CTG_STRND, CTG_ENDND,]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    paf_file.close()
    return contig_data

def import_index_path(file_path : str) -> list:
    index_file = open(file_path, "r")
    index_data = []
    for curr_index in index_file:
        curr_index.rstrip()
        if curr_index[0] == '(':
            index_data.append(ast.literal_eval(curr_index))
        elif curr_index[0] != '[':
            temp_list = curr_index.split("\t")
            index_data.append(tuple((int(temp_list[0]), int(temp_list[1]))))
    index_file.close()
    return index_data

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
    
def extract_telomere_connect_contig(contig_data : list, graph_adjacency : dict) -> list:
    contig_data_size = len(contig_data)
    telomere_connect_contig = []
    for i in range(contig_data_size, contig_data_size+2*CHROMOSOME_COUNT):
        for j in graph_adjacency[(DIR_FOR, i)]:
            telomere_connect_contig.append(j[1])
    
    return telomere_connect_contig

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


def extract_nclose_node(nclose_path: str) -> list:
    nclose_list = []
    with open(nclose_path, "r") as f:
        for line in f:
            line = line.split()
            nclose_list.append(int(line[1]))
            nclose_list.append(int(line[2]))
    return nclose_list

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

def find_using_node(contig_data, node_label):
    contig_data_size = len(contig_data)
    using_node = []
    for i in range(contig_data_size):
        if node_label[i] != 0:
            using_node.append(i)
    return using_node


def graph_build(contig_data : list, contig_adjacency : list, linkability_label : list) -> list :
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

def ctg_name_to_int(ctg_name):
    return re.search(r'\d+$', ctg_name[:-1]).group()

def edge_optimization(contig_data : list, contig_adjacency : list) -> list :
    contig_pair_nodes = defaultdict(list)
    contig_data_size = len(contig_data)
    optimized_adjacency = [[[] for _ in range(contig_data_size)], [[] for _ in range(contig_data_size)]]
    for _ in range(2):
        for i in range(contig_data_size):
            for edge in contig_adjacency[_][i]:
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
    for i in range(2):
        for j in range(contig_data_size):
            first_contig_name = contig_data[j][CTG_NAM]
            fcn_int = int(ctg_name_to_int(first_contig_name))
            for edge in contig_adjacency[i][j]:
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

    return optimized_adjacency

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

def connect_path_folder(folder_path):
    index_file_paths = glob.glob(folder_path + "/*index*")
    for index_file_path in index_file_paths:
        chr_folder = index_file_path.split('/')[-2]
        connect_path(f'{prefix}/10_fill/{chr_folder}', index_file_path)

def connect_path(folder_path, index_file_path):
    os.makedirs(f"{folder_path}/", exist_ok=True)
    cnt = int(index_file_path.split(".")[-3].split("/")[-1])
    with open(f"{folder_path}/{cnt}.paf", "wt") as g, \
         open(f"{folder_path}/{cnt}.index.txt", "wt") as h:
        bnd_index_path = import_index_path(index_file_path)
        str_telo = bnd_index_path[0]
        end_telo = bnd_index_path[-1]
        bnd_index_path = bnd_index_path[1:-1]
        bnd_index_path_size = len(bnd_index_path)
        print(str_telo, file=g)
        for i in range(1, bnd_index_path_size):
            s = bnd_index_path[i-1][1]
            e = bnd_index_path[i][1]
            print(contig_data[s], file=g)
            print(tuple(bnd_index_path[i-1][0:2]), file=h)
            if contig_data[s][CTG_NAM] != contig_data[e][CTG_NAM]:
                if nx.has_path(G_obj, source=(DIR_OUT, s), target=(DIR_IN, e)):
                    path = nx.shortest_path(G_obj, source=(DIR_OUT, s), target=(DIR_IN, e), weight='weight')
                    for node in path[1:-1]:
                        assert(contig_data[node[1]][CTG_TYP] == 3)
                        print(contig_data[node[1]], file=g)
                        print(tuple(node[0:2]), file=h)
            else:
                if s<e:
                    for i in range(s+1, e):
                        print(contig_data[i], file=g)
                        print(tuple((DIR_FOR, i)), file=h)
                else:
                    for i in range(s-1, e, -1):
                        print(contig_data[i], file=g)
                        print(tuple((DIR_BAK, i)), file=h)
        last_path_idx = bnd_index_path[-1][1]
        print(contig_data[last_path_idx], file=g)
        print(tuple((DIR_IN, last_path_idx)), file=h)
        print(end_telo, file=g)
                    
def pool_init(contig_data_, G, prefix_):
    global G_obj, contig_data, prefix

    G_obj = G
    contig_data = contig_data_
    prefix = prefix_

def main():
    parser = argparse.ArgumentParser(description="Find breakend contigs with contig data and map data")
    
    # 위치 인자 정의
    parser.add_argument("ppc_paf_file_path", 
                        help="Path to the preprocessed PAF file.")
    
    parser.add_argument("graph_file_txt", 
                    help="Path to the graph text file.")

    parser.add_argument("prefix", 
                    help="Pefix for pipeline")
    
    parser.add_argument("-t", "--thread", 
                        help="Number of thread", type=int)
    
    parser.add_argument("--progress", 
                        help="Show progress bar", action='store_true')
    
    args = parser.parse_args()

    graph_data = args.graph_file_txt

    PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path

    BREAKEND_GRAPH_PATH_FILE_PREFIX = args.prefix

    NCLOSE_FILE_PATH = f"{args.prefix}/nclose_nodes_index.txt"

    VIRTUAL_TYPE_3_PATH = f"{args.prefix}/virtual_ordinary_contig.txt"

    THREAD = args.thread
    
    contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
    virtual_type_3_data = import_data2(VIRTUAL_TYPE_3_PATH)
    #contig_data += virtual_type_3_data
    contig_data_size = len(contig_data)
    
    init_graph_adjacency = initial_graph_build(contig_data)
    link_label = node_label(contig_data)
    using_node = find_using_node(contig_data, link_label)
    ord_graph_adjacency = graph_build(contig_data, init_graph_adjacency, link_label)
    with open(f"{BREAKEND_GRAPH_PATH_FILE_PREFIX}/non_opt_bnd_connect_graph.txt", "wt") as f:
        for dir in range(2):
            for node in range(contig_data_size):
                print((dir, node), ":", ord_graph_adjacency[dir][node], file=f)
    opt_ord_graph_adjacency = edge_optimization(contig_data, ord_graph_adjacency)
    type_3_graph = defaultdict(list)
    with open(f"{BREAKEND_GRAPH_PATH_FILE_PREFIX}/opt_bnd_connect_graph.txt", "wt") as f:
        for dir in range(2):
            for node in range(contig_data_size):
                print((dir, node), ":", opt_ord_graph_adjacency[dir][node], file=f)
                for edge in opt_ord_graph_adjacency[dir][node]:
                    type_3_graph[(dir, node)].append(edge)

    graph_adjacency = import_data(graph_data)

    bnd_contig = extract_bnd_contig(contig_data)

    telo_contig = extract_telomere_connect_contig(contig_data, graph_adjacency)

    nclose_nodes = extract_nclose_node(NCLOSE_FILE_PATH)

    bnd_connected_graph = connect_nclose_telo(contig_data, using_node, type_3_graph, nclose_nodes, telo_contig)
    node_cnt = 0
    edge_cnt = 0
    with open(f"{BREAKEND_GRAPH_PATH_FILE_PREFIX}/full_bnd_connect_graph.txt", "wt") as f:
        for _ in bnd_connected_graph:
            node_cnt += 1
            print(_, ":", end=" ", file=f)
            for __ in bnd_connected_graph[_]:
                edge_cnt +=1
                print(__, end=" ", file=f)
            print("", file=f)
    #print(f"Node : {node_cnt}, Edge: {edge_cnt}")

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

    chr_chr_folder_path = glob.glob(BREAKEND_GRAPH_PATH_FILE_PREFIX+"/00_raw/*")

    with ProcessPoolExecutor(max_workers=THREAD, initializer=pool_init, initargs=(contig_data, G, args.prefix)) as executor:
        futures = []
        for folder_path in chr_chr_folder_path:
            futures.append(executor.submit(connect_path_folder, folder_path))
    
        # 제출된 작업들이 완료될 때까지 진행 상황을 tqdm으로 표시합니다.
        for future in tqdm(as_completed(futures), total=len(futures), desc='Fill breakend graph', disable=not sys.stdout.isatty() and not args.progress):
            future.result()

if __name__ == "__main__":
    main()