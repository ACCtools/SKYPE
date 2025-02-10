import re
import argparse
import pandas as pd
from collections import defaultdict
import ast
import matplotlib.pyplot as plt

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

    return graph_adjacency

def import_data2(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.split("\t")
        int_induce_idx = [1, 2, 3, 6, 7, 8, 9, 10, 11]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    return contig_data


CHROMOSOME_COUNT = 23
PREPROCESSED_PAF_FILE_PATH = '20_acc_pipe/U2OS_telo.r/U2OS_telo.r.aln.paf.ppc.paf'
graph_data = '20_acc_pipe/U2OS_telo.r/U2OS_telo.r.aln.paf.ppc.paf.op.graph.txt'

contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
contig_data_size = len(contig_data)
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

graph_adjacency = import_data(graph_data)

breakend_adjacency = {}
BOUNDARY = 10000*K
cnt = 0
a = []

breakend_list = set()
cnt = 0

def extract(contig : list) -> list:
    return [contig[CTG_NAM], contig[CHR_STR], contig[CHR_END]]


def chr2int(x):
    chrXY2int = {'chrX' : 24, 'chrY' : 25}
    if x in chrXY2int:
        return chrXY2int[x]
    else:
        return int(x[3:])


with open("endpoint_different_contigs.txt", "wt") as f:
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    data_dict = defaultdict(list)
    for i in range(contig_data_size):
        if contig_data[i][CTG_TYP] ==1 and contig_data[i][CTG_NAM] != contig_data[i+1][CTG_NAM]:
            curr_ctg = contig_data[i]
            str_node = contig_data[curr_ctg[CTG_STRND]]
            end_node = contig_data[curr_ctg[CTG_ENDND]]
            if len(str_node[CTG_TELDIR])>1 or len(end_node[CTG_TELDIR])>1:
                continue
            cnt+=1
            str_nam = str_node[CHR_NAM]
            end_nam = end_node[CHR_NAM]
            if chr2int(str_nam) > chr2int(end_nam):
                data_dict[(end_nam, str_nam)].append(i)
            else:
                data_dict[(str_nam, end_nam)].append(i)
    for _ in data_dict:
        f.write(f"{_[0]},{_[1]},{len(data_dict[_])}\n")
        ctgs = []
        for i in data_dict[_]:
            if contig_data[contig_data[i][CTG_STRND]][CHR_NAM] < contig_data[contig_data[i][CTG_ENDND]][CHR_NAM]:
                str_ctg = contig_data[contig_data[i][CTG_STRND]]
                end_ctg = contig_data[contig_data[i][CTG_ENDND]]
            else:
                str_ctg = contig_data[contig_data[i][CTG_ENDND]]
                end_ctg = contig_data[contig_data[i][CTG_STRND]]
            str_ctg = extract(str_ctg)
            end_ctg = extract(end_ctg)
            ctgs.append([str_ctg, end_ctg])
        ctgs = sorted(ctgs, key=lambda t : t[0][1])
        for i in ctgs:
            print(i[0], i[1], sep = " ", file=f)
            
            
            
        f.write("\n")
    f.write(str(cnt))
    contig_data = contig_data[:-1]
with open("breakends.txt", "wt") as f:
    f.write("Breakend 1 : chromosome type changes\n")
    # Breakend 1 : chromosome type changes
    for i in range(1, contig_data_size):
        if contig_data[i][CTG_NAM]==contig_data[i-1][CTG_NAM] \
        and contig_data[i][CHR_NAM] != contig_data[i-1][CHR_NAM]:
            breakend_list.add(i)
            cnt+=1
            f.write(f'{contig_data[i-1][CHR_NAM]},{contig_data[i][CHR_NAM]}\n')
            f.write(str(contig_data[i-1])+"\n")
            f.write(str(contig_data[i])+"\n\n")
    a.append(cnt)
    cnt=0
    # Breakend 2 : Direction changes
    f.write("Breakend 2 : Direction changes\n")
    for i in range(contig_data_size):
        if contig_data[i][CTG_NAM]==contig_data[i-1][CTG_NAM] \
        and contig_data[i][CHR_NAM] == contig_data[i-1][CHR_NAM] \
        and contig_data[i][CTG_DIR] != contig_data[i-1][CTG_DIR]:
            breakend_list.add(i)
            cnt+=1
            f.write(str(contig_data[i-1])+"\n")
            f.write(str(contig_data[i])+"\n\n")
    a.append(cnt)
    cnt=0
    # Breakend 3 : Big changes
    f.write("Breakend 3 : Big changes\n")
    for i in range(contig_data_size):
        if contig_data[i][CTG_NAM]==contig_data[i-1][CTG_NAM] \
        and contig_data[i][CHR_NAM] == contig_data[i-1][CHR_NAM] \
        and abs(contig_data[i][CHR_STR] - contig_data[i-1][CHR_END]) > BOUNDARY:
            breakend_list.add(i)
            cnt+=1
            f.write(str(contig_data[i-1])+"\n")
            f.write(str(contig_data[i])+"\n\n")
    a.append(cnt)
    cnt=0
    f.write("Breakend 4 : Node connected with telomere nodes\n")
    # Breakend 4 : Node connected with telomere nodes
    for j in range(2):
        for i in range(contig_data_size, contig_data_size+2*CHROMOSOME_COUNT):
            try:
                for _ in graph_adjacency[(j, i)]:
                    f.write(str(chr_rev_corr[i])+"\n")
                    f.write(str(contig_data[_[1]])+"\n\n")
                    cnt+=1
            except:
                pass
    a.append(cnt)
    cnt=0
    print(sum(a), a)

