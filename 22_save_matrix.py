import numpy as np
import pandas as pd
import pickle as pkl

import os
import ast
import sys
import h5py
import logging
import argparse
import collections
import glob

from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("22_save_matrix start")


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

ABS_MAX_COVERAGE_RATIO = 3
K = 1000
M = K * 1000
MAX_PATH_CNT = 100

CHROMOSOME_COUNT = 23
DIR_FOR = 1
TELOMERE_EXPANSION = 5 * K


def import_index_path(file_path : str) -> list:
    file_path_list = file_path.split('/')
    key = file_path_list[-2]
    cnt = int(file_path_list[-1].split('.')[0]) - 1

    return path_list_dict[key][cnt][0]

def import_index_cnt(file_path : str) -> list:
    file_path_list = file_path.split('/')
    key = file_path_list[-2]
    cnt = int(file_path_list[-1].split('.')[0]) - 1

    return path_list_dict[key][cnt][1]

def extract_groups(lst):
    if not lst:
        return []

    result = []
    seen = set()
    current = lst[0]
    result.append(current)
    seen.add(current)

    for num in lst[1:]:
        if num != current:
            # 새로운 숫자가 등장했는데 이미 이전에 등장한 적이 있다면 에러 처리
            if num in seen:
                raise ValueError(f"Error: {num}가 연속된 구간 이후에 다시 등장합니다.")
            result.append(num)
            seen.add(num)
            current = num
    return result


def find_chr_len(file_path: str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    chr_data_file.close()
    return chr_len


def chr_correlation_maker(contig_data):
    chr_corr = {}
    chr_rev_corr = {}
    contig_data_size = len(contig_data)
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr' + str(i) + 'f'] = contig_data_size + i - 1
        chr_rev_corr[contig_data_size + i - 1] = 'chr' + str(i) + 'f'
    chr_corr['chrXf'] = contig_data_size + CHROMOSOME_COUNT - 1
    chr_corr['chrYf'] = contig_data_size + CHROMOSOME_COUNT - 1
    chr_rev_corr[contig_data_size + CHROMOSOME_COUNT - 1] = 'chrXf'
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr' + str(i) + 'b'] = contig_data_size + CHROMOSOME_COUNT + i - 1
        chr_rev_corr[contig_data_size + CHROMOSOME_COUNT + i - 1] = 'chr' + str(i) + 'b'
    chr_corr['chrXb'] = contig_data_size + 2 * CHROMOSOME_COUNT - 1
    chr_corr['chrYb'] = contig_data_size + 2 * CHROMOSOME_COUNT - 1
    chr_rev_corr[contig_data_size + 2 * CHROMOSOME_COUNT - 1] = 'chrXb'
    return chr_corr, chr_rev_corr


def import_telo_data(file_path: str, chr_len: dict) -> dict:
    fai_file = open(file_path, "r")
    telo_data = [(0, 1)]
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        int_induce_idx = [1, 2]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        if temp_list[0] != telo_data[-1][0]:
            temp_list[2] += TELOMERE_EXPANSION
            temp_list.append('f')
        else:
            if temp_list[1] > chr_len[temp_list[0]] / 2:
                temp_list[1] -= TELOMERE_EXPANSION
                temp_list.append('b')
            else:
                temp_list.append('f')
                temp_list[2] += TELOMERE_EXPANSION
        telo_data.append(tuple(temp_list))
    fai_file.close()
    return telo_data[1:]


def distance_checker(node_a: tuple, node_b: tuple) -> int:
    if max(int(node_a[0]), int(node_b[0])) < min(int(node_a[1]), int(node_b[1])):
        return 0
    else:
        return min(abs(int(node_b[0]) - int(node_a[1])), abs(int(node_b[1]) - int(node_a[0])))


def chr2int(x):
    chrXY2int = {'chrX': 24, 'chrY': 25}
    if x in chrXY2int:
        return chrXY2int[x]
    else:
        return int(x[3:])


def extract_telomere_connect_contig(telo_info_path: str) -> list:
    telomere_connect_contig = []
    with open(telo_info_path) as f:
        for curr_data in f:
            curr_data = curr_data.rstrip()
            temp_list = curr_data.split("\t")
            chr_info = temp_list[0]
            contig_id = ast.literal_eval(temp_list[1])
            telomere_connect_contig.append((chr_info, contig_id[1]))

    return telomere_connect_contig


def rebin_dataframe(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """
    Group rows in the DataFrame into bins spanning n consecutive units.
    The unit is determined from the 'length' of the first row in each chromosome group.

    For each group:
      - The new start (st) is taken from the first row.
      - The new end (nd) is taken from the last row in the group.
      - For a complete group (n rows), values are summed.
      - For an incomplete group (fewer than n rows), the new bin's length is the sum of the available lengths,
        and new covsite and totaldepth are computed as the weighted average:
            new_value = sum( value_i * length_i ) / sum(length_i)
      - New coverage (cov) is computed as (new_covsite / new_length) * 100.
      - New mean depth (meandepth) is computed as new_totaldepth / new_length.

    Parameters:
        df (pd.DataFrame): Input DataFrame with columns ['chr', 'st', 'nd', 'length',
                          'covsite', 'totaldepth', 'cov', 'meandepth'].
        n (int): Number of consecutive units (rows) to combine into each bin.

    Returns:
        pd.DataFrame: New DataFrame with binned rows.
    """
    new_rows = []

    # Process each chromosome separately.
    for chrom, sub_df in df.groupby('chr'):
        # Sort rows by starting position.
        sub_df = sub_df.sort_values('st').reset_index(drop=True)

        # Group rows in chunks of size n.
        for i in range(0, len(sub_df), n):
            chunk = sub_df.iloc[i:i + n]
            new_st = chunk['st'].iloc[0]
            new_nd = chunk['nd'].iloc[-1]

            # Use the actual sum of lengths in the chunk.
            sum_length = chunk['length'].sum()
            new_meandepth = np.sum(chunk['totaldepth']) / sum_length

            new_rows.append({
                'chr': chrom,
                'st': new_st,
                'nd': new_nd,
                'meandepth': new_meandepth
            })

    return pd.DataFrame(new_rows)


def rebin_dataframe_B(df: pd.DataFrame, n: int) -> np.array:
    new_rows = dict()

    # Process each chromosome separately.
    for chrom, sub_df in df.groupby('chr'):
        chr_mean_list = []
        sub_df = sub_df.sort_values('st').reset_index(drop=True)

        # Group rows in chunks of size n.
        for i in range(0, len(sub_df), n):
            chunk = sub_df.iloc[i:i + n]
            new_st = chunk['st'].iloc[0]
            new_nd = chunk['nd'].iloc[-1]

            # Use the actual sum of lengths in the chunk.
            sum_length = chunk['length'].sum()
            new_meandepth = np.sum(chunk['totaldepth']) / sum_length

            chr_mean_list.extend([new_meandepth for _ in range(len(chunk))])

        new_rows[chrom] = np.asarray(chr_mean_list)

    ans_B = []
    for c in chr_order_list:
        ans_B.append(new_rows[c])

    return np.hstack(ans_B)


def import_data2(file_path: str) -> list:
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.split("\t")
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END, \
                          CHR_LEN, CHR_STR, CHR_END, \
                          CTG_MAPQ, CTG_TYP, CTG_STRND, CTG_ENDND, ]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    paf_file.close()
    return contig_data


def import_bed(bed_path: str) -> dict:
    bed_data_file = open(bed_path, "r")
    chr_len = collections.defaultdict(list)
    for curr_data in bed_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]].append((int(curr_data[1]), int(curr_data[2])))
    bed_data_file.close()
    return chr_len


def inclusive_checker(tuple_a: tuple, tuple_b: tuple) -> bool:
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False


def get_final_paf_name_from_index(index_file_path):
    final_paf_list = index_file_path.split('/')
    cnt = final_paf_list[-1].split('.')[0]
    final_paf_list[-3] = '20_depth'
    final_paf_path = '/'.join(final_paf_list[:-1]) + f'/{cnt}.paf'

    return final_paf_path


def np_safe_divide(a, b):
    out = np.zeros_like(a, dtype=np.float64)

    out[(b == 0) & (a > 0)] = np.inf
    out[(b == 0) & (a < 0)] = -np.inf

    return np.divide(a, b, out=out, where=(b != 0))

def get_relative_path(p):
    return tuple(p.split('/')[-3:])

parser = argparse.ArgumentParser(description="SKYPE depth analysis")

parser.add_argument("censat_bed_path",
                    help="Path to the censat repeat information file.")

parser.add_argument("ppc_paf_file_path",
                    help="Path to the preprocessed PAF file.")

parser.add_argument("main_stat_loc",
                    help="Cancer coverage location file")

parser.add_argument("telomere_bed_path",
                    help="Path to the telomere information file.")

parser.add_argument("reference_fai_path",
                    help="Path to the chromosome information file.")

parser.add_argument("reference_cytobands_path",
                    help="Path to the cytoband information file.")

parser.add_argument("prefix",
                    help="Pefix for pipeline")

parser.add_argument("-t", "--thread",
                    help="Number of thread", type=int)

parser.add_argument("--progress",
                    help="Show progress bar", action='store_true')

args = parser.parse_args()

# t = "22_save_matrix.py public_data/chm13v2.0_censat_v2.1.m.bed 20_acc_pipe/Caki-1.p/Caki-1.p.aln.paf.ppc.paf /home/hyunwoo/51g_cancer_denovo/51_depth_data/Caki-1.win.stat.gz public_data/chm13v2.0_telomere.bed public_data/chm13v2.0.fa.fai public_data/chm13v2.0_cytobands_allchrs.bed 30_skype_pipe/Caki-1_03_00_34 -t 64"

# args = parser.parse_args(t.split()[1:])

bed_data = import_bed(args.censat_bed_path)

PREFIX = args.prefix
THREAD = args.thread
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
main_stat_loc = args.main_stat_loc
TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
front_contig_path = RATIO_OUTLIER_FOLDER + "front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER + "back_jump/"
TELO_CONNECT_NODES_INFO_PATH = PREFIX + "/telomere_connected_list.txt"

contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
telo_connected_node = extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH)

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)

telo_connected_node_dict = collections.defaultdict(list)
for chr_info, contig_id in telo_connected_node:
    telo_connected_node_dict[chr_info].append(contig_id)

telo_chr_list = list(set(n[:-1] for n in telo_connected_node_dict.keys()))

default_path_list = []
for chr_name in telo_chr_list:
    chr_def_paths = []
    for tar_chr_name, ctg_idx, func in [(chr_name + 'f', CHR_STR, min), (chr_name + 'b', CHR_END, max)]:
        tar_contig_data_list = [(i, contig_data[i]) for i in telo_connected_node_dict[tar_chr_name]]
        tar_contig_data = func(tar_contig_data_list, key=lambda t: t[1][ctg_idx])
        chr_def_paths.append((tar_chr_name, tar_contig_data[0]))

    default_path_list.append(chr_def_paths)

tar_chr_data = dict()
tar_def_path_list = []

for (chrf, chrfid), (chrb, chrbid) in default_path_list:
    tar_chr = chrf[:-1]

    key = f'{chrf}_{chrb}'
    n = len(path_list_dict[key])
    tar_path_loc = None

    for i in range(n):
        path_loc = f"{PREFIX}/00_raw/{key}/{i + 1}.inedx.txt"
        idx_data = import_index_path(path_loc)
        if {idx_data[1][1], idx_data[-2][1]} == {chrfid, chrbid} and (idx_data[-1][1], idx_data[-1][2]) == (0, 0):
            tar_path_loc = path_loc
            break

    if tar_path_loc is None:
        non_bnd_path_list = []
        for i in range(n):
            path_loc = f"{PREFIX}/00_raw/{key}/{i + 1}.inedx.txt"
            idx_data = import_index_path(path_loc)
            if (idx_data[-1][1], idx_data[-1][2]) != (0, 0):
                break

            path_counter_list = import_index_cnt(path_loc)

            for ch, v in path_counter_list:
                if ch == tar_chr:
                    non_bnd_path_list.append((path_loc, v))

        tar_path_loc = max(non_bnd_path_list, key=lambda t: t[1])[0]

    tar_final_path_loc = get_relative_path(get_final_paf_name_from_index(tar_path_loc))
    tar_chr_data[tar_chr] = tar_final_path_loc
    tar_def_path_list.append(tar_final_path_loc)

tar_def_path_set = set(tar_def_path_list)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t',
                 names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

ydf = df.query('chr == "chrY"')
chry_nz_len = len(ydf.query('meandepth != 0'))
no_chrY = (chry_nz_len / len(ydf)) < 0.5

meandepth = np.median(df['meandepth'])
chr_order_list = extract_groups(list(df['chr']))

chr_filt_idx_dict = collections.defaultdict(list)
chr_no_filt_idx_dict = collections.defaultdict(list)

chr_filt_st_list = []
chr_no_filt_st_list = []

chr_filt_idx_list = []
chr_no_filt_idx_list = []

for ind, l in enumerate(df.itertuples(index=False)):
    flag = True
    if l.meandepth > ABS_MAX_COVERAGE_RATIO * meandepth:
        flag = False
    for i in bed_data:
        if l.chr == i:
            for j in bed_data[i]:
                if inclusive_checker(j, (l.st, l.nd)):
                    flag = False
                    break
        if not flag:
            break

    if flag:
        chr_filt_st_list.append((l.chr, l.st))
        chr_filt_idx_list.append(ind)
    else:
        chr_no_filt_st_list.append((l.chr, l.st))
        chr_no_filt_idx_list.append(ind)

filter_len = len(chr_filt_st_list)

for i, (chr_, st_) in enumerate(chr_filt_st_list):
    chr_filt_idx_dict[chr_].append(i)

for i, (chr_, st_) in enumerate(chr_no_filt_st_list):
    chr_no_filt_idx_dict[chr_].append(i + filter_len)


def get_vec_from_stat_loc(stat_loc_):
    df = pd.read_csv(stat_loc_, compression='gzip', comment='#', sep='\t',
                     names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
    df = df.query('chr != "chrM"')

    chr_st_data = dict()
    for l in df.itertuples(index=False):
        chr_st_data[(l.chr, l.st)] = l.meandepth

    v = []
    for cs in chr_filt_st_list:
        if cs not in chr_st_data:
            v.append(0)
        else:
            v.append(chr_st_data[cs])

    for cs in chr_no_filt_st_list:
        if cs not in chr_st_data:
            v.append(0)
        else:
            v.append(chr_st_data[cs])

    return np.asarray(v, dtype=np.float32)


def get_vec_from_ki(ki):
    stat_loc_ = f'{output_folder}/{ki}.win.stat.gz'
    df = pd.read_csv(stat_loc_, compression='gzip', comment='#', sep='\t',
                     names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
    df = df.query('chr != "chrM"')

    chr_st_data = dict()
    for l in df.itertuples(index=False):
        chr_st_data[(l.chr, l.st)] = l.meandepth

    v = []
    for cs in chr_filt_st_list:
        if cs not in chr_st_data:
            v.append(0)
        else:
            v.append(chr_st_data[cs])

    for cs in chr_no_filt_st_list:
        if cs not in chr_st_data:
            v.append(0)
        else:
            v.append(chr_st_data[cs])

    return ki, np.asarray(v, dtype=np.float32)


PATH_FILE_FOLDER = f"{PREFIX}/20_depth"
chr_chr_folder_path = sorted(glob.glob(PATH_FILE_FOLDER + "/*"))

chr_st_data = dict()
for l in df.itertuples(index=False):
    if l.chr == 'chrY' and no_chrY:
        chr_st_data[(l.chr, l.st)] = 0
    else:
        chr_st_data[(l.chr, l.st)] = l.meandepth

v = []
for cs in chr_filt_st_list:
    if cs not in chr_st_data:
        v.append(0)
    else:
        v.append(chr_st_data[cs])

for cs in chr_no_filt_st_list:
    if cs not in chr_st_data:
        v.append(0)
    else:
        v.append(chr_st_data[cs])

B = np.asarray(v, dtype=np.float32)

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    paf_ans_list, key_list, int2key = pkl.load(f)

vec_dict = [None] * (len(key_list) + 1)
output_folder = f'{PREFIX}/21_pat_depth'
with ProcessPoolExecutor(max_workers=THREAD) as executor:
    futures = []
    for ki in key_list:
        futures.append(executor.submit(get_vec_from_ki, ki))
    for future in tqdm(as_completed(futures), total=len(futures), desc='Parse depth for each seperated paths\' gz file',
                       disable=not sys.stdout.isatty() and not args.progress):
        i, v = future.result()
        vec_dict[i] = v

fclen = len(glob.glob(front_contig_path + "*"))
bclen = len(glob.glob(back_contig_path + "*"))

m = np.shape(B)[0]
n = len(paf_ans_list) + fclen // 4 + bclen // 4
ncnt = 0

A = np.empty((n, m), dtype=np.float32, order='C')

filter_vec_list = []
tot_loc_list = []
bv_loc_list = []

tmp_v = np.zeros(m, dtype=np.float32)

tar_def_path_ind_dict = dict()
for path, key_int_list in tqdm(paf_ans_list, desc='Recover depth from seperated paths',
                 disable=not sys.stdout.isatty() and not args.progress):
    ki = key_int_list[0]
    np.copyto(tmp_v, vec_dict[ki])

    for ki in key_int_list[1:]:
        tmp_v += vec_dict[ki]
    A[ncnt, :] = tmp_v

    path = get_relative_path(path)
    if path in tar_def_path_set:
        tar_def_path_ind_dict[path] = ncnt

    ncnt += 1

for_dir_data = []
for i in tqdm(range(1, fclen // 4 + 1), desc='Parse coverage from forward-directed outlier contig gz files',
              disable=not sys.stdout.isatty() and not args.progress):
    bv_paf_loc = front_contig_path + f"{i}_base.paf"
    ov_loc = front_contig_path + f"{i}.win.stat.gz"
    bv_loc = front_contig_path + f"{i}_base.win.stat.gz"
    bv_loc_list.append(bv_paf_loc)
    ov = get_vec_from_stat_loc(ov_loc)
    bv = get_vec_from_stat_loc(bv_loc)

    vec = ov - bv

    with open(bv_paf_loc, 'r') as f:
        l = f.readline().strip()
        ll = l.split('\t')
        tar_chr = ll[CHR_NAM]

    tar_def_loc = tar_chr_data[tar_chr]
    temp_ncnt = tar_def_path_ind_dict[tar_def_loc]

    ref_vec = A[temp_ncnt, :][chr_filt_idx_dict[tar_chr]]
    ref_loc_vec = vec[chr_filt_idx_dict[tar_chr]]
    ref_loc_vec[ref_vec == 0] = 0

    vec_div = np_safe_divide(ref_loc_vec, ref_vec)
    ref_weight = -np.min(vec_div)

    if ref_weight > 5:
        ref_weight = 0
        vec = np.zeros_like(vec)
    elif ref_weight > 0:
        vec[chr_filt_idx_dict[tar_chr]] = ref_loc_vec
        vec += A[temp_ncnt, :] * ref_weight
    else:
        ref_weight = 0

    A[ncnt, :] = vec
    for_dir_data.append([temp_ncnt, ncnt, ref_weight])

    ncnt += 1

with open(f'{PREFIX}/for_dir_data.pkl', 'wb') as f:
    pkl.dump(for_dir_data, f)

for i in tqdm(range(1, bclen // 4 + 1), desc='Parse coverage from backward-directed outlier contig gz files',
              disable=not sys.stdout.isatty() and not args.progress):
    bv_paf_loc = back_contig_path + f"{i}_base.paf"
    ov_loc = back_contig_path + f"{i}.win.stat.gz"
    bv_loc = back_contig_path + f"{i}_base.win.stat.gz"
    bv_loc_list.append(bv_paf_loc)
    ov = get_vec_from_stat_loc(ov_loc)
    bv = get_vec_from_stat_loc(bv_loc)

    A[ncnt, :] = ov + bv
    ncnt += 1

with h5py.File(f'{PREFIX}/matrix.h5', 'w') as hf:
    dset_A = hf.create_dataset('A', shape=A[:, :filter_len].shape, dtype=A.dtype)
    dset_B = hf.create_dataset('B', shape=B[:filter_len].shape, dtype=B.dtype)
    dset_A_fail = hf.create_dataset('A_fail', shape=A[:, filter_len:].shape, dtype=A.dtype)
    dset_B_fail = hf.create_dataset('B_fail', shape=B[filter_len:].shape, dtype=B.dtype)

    dset_A.write_direct(A, source_sel=np.s_[:, :filter_len])
    dset_B.write_direct(B, source_sel=np.s_[:filter_len])
    dset_A_fail.write_direct(A, source_sel=np.s_[:, filter_len:])
    dset_B_fail.write_direct(B, source_sel=np.s_[filter_len:])
