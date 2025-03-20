import sys
import subprocess
import argparse
import os
from collections import defaultdict
import ast
import glob
import pickle as pkl
from itertools import pairwise, groupby

from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging

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
CTG_MAINFLOWDIR = 18
CTG_MAINFLOWCHR = 19

BND_TYPE = 0
CTG_IN_TYPE = 1
TEL_TYPE = 2

DIR_FOR = 1
DIR_BAK = 0
DIR_IN = 3
DIR_OUT = 2
INF = 1000000000
BUFFER = 10000000
CHROMOSOME_COUNT = 23
K = 1000

DEPTH_WINDOW=100 * K
DEPTH_THREAD=1

NCLOSE_COMPRESS_LIMIT = 50*K
VIRTUAL_CONTIG_PREFIX = "virtual_contig"

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
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END,
                          CHR_LEN, CHR_STR, CHR_END,
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

def import_fill_index_path(file_path : str) -> list:
    index_file = open(file_path, "r")
    index_data = []
    for curr_index in index_file:
        curr_index = curr_index.rstrip()
        index_data.append(ast.literal_eval(curr_index))
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

def compare_paf_list(l1, l2):
    if len(l1) != len(l2):
        return False
    
    for p1, p2 in zip(l1, l2):
        if p1 != p2:
            if p1[CTG_NAM].startswith(VIRTUAL_CONTIG_PREFIX) == p2[CTG_NAM].startswith(VIRTUAL_CONTIG_PREFIX):
                continue
            return False
    return True

def compare_paf_list_test(l1, l2):
    if len(l1) != len(l2):
        print(len(l1), len(l2))
        return False
    
    for p1, p2 in zip(l1, l2):
        if p1 != p2:
            if p1[CTG_NAM].startswith(VIRTUAL_CONTIG_PREFIX) == p2[CTG_NAM].startswith(VIRTUAL_CONTIG_PREFIX):
                continue
            for i1, i2 in zip(p1, p2):
                if i1 != i2:
                    print(i1, i2)
            return False
    return True

def path_printer(index_file_path):
    final_paf_list = index_file_path.split('/')
    cnt = final_paf_list[-1].split('.')[0]
    
    final_paf_list[-3] = '10_fill'
    fill_index_path_loc = '/'.join(final_paf_list)

    final_paf_list[-3] = '20_depth'
    final_paf_path = '/'.join(final_paf_list[:-1]) + f'/{cnt}.paf'

    print(os.path.abspath(final_paf_path))
    print(os.path.abspath(fill_index_path_loc))
    print(os.path.abspath(index_file_path))

def split_final_paf(index_file_path):
    final_paf_list = index_file_path.split('/')
    cnt = final_paf_list[-1].split('.')[0]
    
    final_paf_list[-3] = '10_fill'
    fill_index_path_loc = '/'.join(final_paf_list)
    fill_index_path = import_fill_index_path(fill_index_path_loc)

    final_paf_list[-3] = '20_depth'
    final_paf_path = '/'.join(final_paf_list[:-1]) + f'/{cnt}.paf'

    bnd_index_data = import_index_path(index_file_path)[1:-1]
    final_paf_data = import_final_paf(final_paf_path)

    bnd_ctg_group = count_groups([contig_data[i[1]][CTG_NAM] for i in bnd_index_data])
    # print(final_paf_path)
    # print(fill_index_path_loc)
    # print(index_file_path)

    index_list = []
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
            ii = find_index(fill_index_path, final_paf_data, ci)
            index_list.append((ii, ii, TEL_TYPE,
                               ctg_ind))

        if contig_data[ctg_ind][CTG_NAM] == contig_data[next_ctg_ind][CTG_NAM] and is_contig_same:
            ii1 = find_index(fill_index_path, final_paf_data, ci)
            ii2 = find_index(fill_index_path, final_paf_data, nci)

            index_list.append((ii1 + 1, ii2 - 1, CTG_IN_TYPE,
                               (ci, nci)))
        else:
            ii1 = find_index(fill_index_path, final_paf_data, ci)
            ii2 = find_index(fill_index_path, final_paf_data, nci)

            index_list.append((ii1, ii2, BND_TYPE,
                               (ci, nci)))
            
            
        if i == len(bnd_index_data) - 2 and bnd_ctg_group[-1] == 2:
            ii = find_index(fill_index_path, final_paf_data, nci)
            index_list.append((ii, ii, TEL_TYPE,
                               ctg_ind))

    # print(*index_list)
    for (si1, si2, st, si), (ni1, ni2, nt, si) in pairwise(index_list):
        assert(si2 + 1 == ni1)
    
    assert(index_list[0][0] == 0)
    assert(index_list[-1][1] == len(final_paf_data) - 1)

    key_list = []
    for i1, i2, data_type, data_key in index_list:
        key = (data_type, data_key)
        if key not in index_data:
            index_data[key] = final_paf_data[i1:i2+1]
            # index_data_loc[key] = index_file_path
        else:
            assert(compare_paf_list(index_data[key], final_paf_data[i1:i2+1]))

        key_list.append(key)
    
    return final_paf_path, key_list

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


parser = argparse.ArgumentParser(description="Find breakend contigs with contig data and map data")

# 위치 인자 정의
parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")

parser.add_argument("prefix", 
                help="Pefix for pipeline")

parser.add_argument("-t", "--thread", 
                help="Number of thread", type=int)

parser.add_argument("--progress", 
                    help="Show progress bar", action='store_true')

args = parser.parse_args()

PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path

THREAD = args.thread
PREFIX = args.prefix

PATH_FILE_FOLDER = f"{PREFIX}/20_depth/"
RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
back_front_folder_path = glob.glob(RATIO_OUTLIER_FOLDER+"*")

contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)

bnd_dir_dict = defaultdict(set)

tot_ind_set = [set() for _ in range(3)]
tot_ind_data = [[] for _ in range(3)]

norm_dir_type = {DIR_FOR : True, DIR_BAK : True, DIR_IN : False, DIR_OUT : False}

for k, vs in bnd_dir_dict.items():
    type_set = set((norm_dir_type[i], norm_dir_type[j]) for i, j in vs)
    assert(len(type_set) == len(vs))

index_data = dict()

chr_chr_folder_path = list(glob.glob(PREFIX+"/00_raw/*"))
final_paf_vec_data = []

# tqdm 전에꺼랑 같이 desc 적고
for folder_path in tqdm(chr_chr_folder_path, desc='Parse and verify final paf data',
                        disable=not sys.stdout.isatty() and not args.progress):
    index_file_paths = glob.glob(folder_path + "/*index*")
    for index_file_path in index_file_paths:
        final_paf_vec_data.append(split_final_paf(index_file_path))

key2int = dict()

output_folder = f'{PREFIX}/21_pat_depth'
os.makedirs(output_folder, exist_ok=True)

for i, (key, paf_list) in enumerate(index_data.items()):
    key2int[key] = i
    with open(f'{output_folder}/{i}.paf', 'w') as f:
        for paf_line in paf_list:
            print(*paf_line, sep='\t', file=f)

paf_ans_list = []
for final_paf_path, key_list in final_paf_vec_data:
    key_list = [key2int[k] for k in key_list]
    paf_ans_list.append((final_paf_path, key_list))

def get_paf_run(paf_loc):
    paf_base = os.path.splitext(paf_loc)[0]
    result = subprocess.run(['./PanDepth/bin/pandepth', '-w', str(int(DEPTH_WINDOW)),'-t', str(DEPTH_THREAD), '-i', paf_loc, '-o', paf_base],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, capture_output=False)
    
    if result.returncode != 0:
        raise Exception("Pandepth failed!")

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

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'wb') as f:
    pkl.dump((paf_ans_list, list(key2int.values())), f)
