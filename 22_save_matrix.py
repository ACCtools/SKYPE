import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *
from skype_output_files import (
    discover_ecdna_depth_inputs,
    discover_jump_depth_inputs,
)
from raw_translocation_prior import (
    RAW_TRANSLOCATION_REPORT_TSV,
    calculate_length_normalized_prior_scales,
    load_normal_prior_excluded_chromosomes,
    select_normal_prior_chromosome_columns,
)
from nclose_tracking import (
    bnd_event_keys,
    compressed_bnd_event_keys,
    count_ecdna_circuit_events,
    count_index_path_events,
    event_catalog_by_key,
    indel_event_key,
    initialise_filter_status,
    load_event_catalog,
    load_type4_edge_event_map,
    save_filter_status,
    save_path_usage,
)

import numpy as np
import pandas as pd
import pickle as pkl

import ast
import h5py
import logging
import argparse
import itertools
import collections
import glob

from collections import Counter, defaultdict

from tqdm import tqdm
from scipy.optimize import nnls
from scipy.signal import butter, filtfilt
from concurrent.futures import ProcessPoolExecutor, as_completed

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).
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
CTG_CENSAT = 18
CTG_MAINFLOWDIR = 19
CTG_MAINFLOWCHR = 20
CTG_GLOBALIDX = 21

ABS_MAX_COVERAGE_RATIO = 3
MAX_PATH_CNT = 100
NORMAL_PRIOR_STRENGTH = 0.00

DIR_FOR = 1
TELOMERE_EXPANSION = 5 * K



def highpass_filter(data, cutoff, fs, order=3):
    """
    Butterworth high-pass filter를 이용하여 저주파 성분을 제거합니다.
    
    Parameters:
    - data: 1D array, 입력 신호 (CN 값)
    - cutoff: float, 컷오프 주파수
    - fs: float, 샘플링 주파수 (데이터 포인트 간 간격)
    - order: int, 필터 차수
    
    Returns:
    - 필터링된 신호 (노이즈 성분 추출)
    """
    nyq = 0.5 * fs  # Nyquist 주파수
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return filtfilt(b, a, data)

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


def import_telo_data(file_path: str, chr_len: dict) -> dict:
    fai_file = open(file_path, "r")
    telo_data = []
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        int_induce_idx = [1, 2]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        if temp_list[1] > chr_len[temp_list[0]] / 2:
            temp_list[1] -= TELOMERE_EXPANSION
            temp_list.append('b')
        else:
            temp_list.append('f')
            temp_list[2] += TELOMERE_EXPANSION
        telo_data.append(tuple(temp_list))
    fai_file.close()
    return telo_data


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

normal_prior_strength = float(NORMAL_PRIOR_STRENGTH)
if normal_prior_strength < 0:
    raise ValueError("NORMAL_PRIOR_STRENGTH must be non-negative")

# t = "22_save_matrix.py public_data/chm13v2.0_censat_v2.1.m.bed /home/hyunwoo/ACCtools-pipeline/90_skype_run/COLO829/20_alignasm/COLO829.ctg.aln.paf.ppc.paf /home/hyunwoo/ACCtools-pipeline/90_skype_run/COLO829/01_depth/COLO829_normalized.win.stat.gz public_data/chm13v2.0_telomere.bed public_data/chm13v2.0.fa.fai public_data/chm13v2.0_cytobands_allchrs.bed 30_skype_pipe/COLO829_14_12_42 -t 4"
# args = parser.parse_args(t.split()[1:])

bed_data = import_bed(args.censat_bed_path)

PREFIX = args.prefix
THREAD = args.thread
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
main_stat_loc = args.main_stat_loc
TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
pipeline_mode_config = load_pipeline_mode(PREFIX)

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
front_contig_path = RATIO_OUTLIER_FOLDER + "front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER + "back_jump/"
ecdna_contig_path = RATIO_OUTLIER_FOLDER + "ecdna/"
type2_ins_contig_path = RATIO_OUTLIER_FOLDER + "type2_ins/"

TELO_CONNECT_NODES_INFO_PATH = PREFIX + "/telomere_connected_list.txt"

if pipeline_mode_is_variant(pipeline_mode_config):
    logging.info("Normal chromosome prior disabled: variant mode")
    normal_prior_strength = 0.0

normal_prior_excluded_chroms = set()
if normal_prior_strength > 0:
    raw_translocation_report_path = os.path.join(
        PREFIX, RAW_TRANSLOCATION_REPORT_TSV
    )
    normal_prior_excluded_chroms = load_normal_prior_excluded_chromosomes(
        raw_translocation_report_path
    )
    if normal_prior_excluded_chroms:
        logging.info(
            "Normal chromosome prior raw-read exclusions: %s",
            ",".join(sorted(normal_prior_excluded_chroms)),
        )

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t',
                 names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
telo_connected_node = extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH)

with open(f'{PREFIX}/nclose_chunk_data.pkl', 'rb') as f:
    nclose_nodes_pkl, _, _ = pkl.load(f)

nclose_set = set()
for vl in nclose_nodes_pkl.values():
    for v in vl:
        nclose_set.add(tuple(sorted(v)))

nclose_event_catalog = load_event_catalog(PREFIX)
nclose_event_by_key = event_catalog_by_key(nclose_event_catalog)
nclose_bnd_keys = bnd_event_keys(nclose_event_catalog)
compressed_nclose_bnd_keys = compressed_bnd_event_keys(nclose_event_catalog)
if nclose_set != compressed_nclose_bnd_keys:
    raise ValueError(
        "NClose BND catalog does not match nclose_chunk_data.pkl. "
        "Rerun SKYPE from stage 02."
    )
type4_edge_to_event = load_type4_edge_event_map(PREFIX, nclose_event_catalog)

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

chrY_type4_skip = False
for (chrf, chrfid), (chrb, chrbid) in default_path_list:
    tar_chr = chrf[:-1]
    key = f'{chrf}_{chrb}'
    n = len(path_list_dict[key])

    # chrY large indel
    if tar_chr == 'chrX':
        non_bnd_path_list = []
        for i in range(n):
            path_loc = f"{PREFIX}/00_raw/{key}/{i + 1}.index.txt"
            idx_data = import_index_path(path_loc)
            path_counter_list = import_index_cnt(path_loc)

            for ch, v in path_counter_list:
                if ch == 'chrY':
                    non_bnd_path_list.append((path_loc, v, idx_data[-1][2]))
        
        if len(non_bnd_path_list) > 0:
            tar_path_loc = max(non_bnd_path_list, key=lambda t: (-t[2], t[1]))[0]

            tar_final_path_loc = get_relative_path(get_final_paf_name_from_index(tar_path_loc))
            tar_chr_data['chrY'] = tar_final_path_loc
            tar_def_path_list.append(tar_final_path_loc)
        else:
            chrY_type4_skip = True

    tar_path_loc = None

    for i in range(n):
        path_loc = f"{PREFIX}/00_raw/{key}/{i + 1}.index.txt"
        idx_data = import_index_path(path_loc)
        if {idx_data[1][1], idx_data[-2][1]} == {chrfid, chrbid} and (idx_data[-1][1], idx_data[-1][2]) == (0, 0):
            tar_path_loc = path_loc
            break

    if tar_path_loc is None:
        non_bnd_path_list = []
        for i in range(n):
            path_loc = f"{PREFIX}/00_raw/{key}/{i + 1}.index.txt"
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

ydf = df.query('chr == "chrY"')

bed_intervals = pd.IntervalIndex.from_tuples(bed_data['chrY'], closed='left')
y_mask = ydf.apply(
    lambda row: bed_intervals.overlaps(pd.Interval(row['st'], row['nd'], closed='left')),
    axis=1
)

correct_mask = y_mask.apply(any)
ydf_not_censat = ydf[~correct_mask]

chry_nz_len = len(ydf_not_censat.query('meandepth != 0'))
no_chrY = (chry_nz_len / len(ydf_not_censat)) < chrY_MINIMUM_RATIO

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

grouped_data = defaultdict(lambda: {"positions": [], "values": []})
for i, (chrom, pos) in enumerate(chr_filt_st_list + chr_no_filt_st_list):
    grouped_data[chrom]["positions"].append(pos)
    grouped_data[chrom]["values"].append(B[i])

# 필터링 파라미터 (데이터 특성에 따라 조정)
cutoff_frequency = 0.1  # 컷오프 주파수 (예시 값)
fs = 1.0                 # 샘플링 주파수 (데이터 포인트 간 간격이 1이라고 가정)
order = 3                # 필터 차수

# 각 염색체별로 필터링 수행
filtered_values_list = []
for chrom, data in grouped_data.items():
    positions = np.array(data["positions"])
    values = np.array(data["values"])
    
    # 만약 위치가 정렬되어 있지 않다면 정렬 (예시에서는 이미 정렬된 상태)
    sort_idx = np.argsort(positions)
    positions = positions[sort_idx]
    values = values[sort_idx]
    
    # 고주파 필터 적용: 기본 CN 트렌드를 제거하여 노이즈 성분 추출
    filtered_values = highpass_filter(values, cutoff=cutoff_frequency, fs=fs, order=order)

    filtered_values_list.append(filtered_values)

noise_array = np.hstack(filtered_values_list)
amplitude = np.std(noise_array)

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    paf_ans_list, key_list, int2key, dep_list = pkl.load(f)

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

front_depth_inputs = discover_jump_depth_inputs(
    front_contig_path, type2_ins_contig_path
)
back_depth_inputs = discover_jump_depth_inputs(
    back_contig_path, type2_ins_contig_path
)
ecdna_depth_inputs = discover_ecdna_depth_inputs(ecdna_contig_path)

with open(f'{PREFIX}/ecdna_circuit_data.pkl', 'rb') as f:
    ecdna_circuits, _ = pkl.load(f)
if len(ecdna_circuits) != len(ecdna_depth_inputs):
    raise ValueError(
        "ecDNA circuit/depth-input count mismatch: "
        f"{len(ecdna_circuits)} != {len(ecdna_depth_inputs)}. "
        "Rerun SKYPE from stage 21."
    )

with open(f'{PREFIX}/cen_fragment_data.pkl', 'rb') as f:
    cen_fragment_meta = pkl.load(f)
cen_fragment_list = sorted(cen_fragment_meta.items(), key=lambda kv: chr2int(kv[0]))

m = np.shape(B)[0]
n = (
    len(paf_ans_list)
    + len(front_depth_inputs)
    + len(back_depth_inputs)
    + len(ecdna_depth_inputs)
    + len(cen_fragment_list)
)
ncnt = 0

tar_def_path_ind_dict = {}
for path_idx, (path, _) in enumerate(paf_ans_list):
    path_rel = get_relative_path(path)
    if path_rel in tar_def_path_set:
        tar_def_path_ind_dict[path_rel] = path_idx

init_cols = [
    tar_def_path_ind_dict[path_rel]
    for path_rel in tar_chr_data.values()
    if path_rel in tar_def_path_ind_dict
]

missing_init_paths = [
    path_rel for path_rel in tar_chr_data.values()
    if path_rel not in tar_def_path_ind_dict
]
if missing_init_paths:
    logging.warning(
        f"Default chromosome paths missing from NNLS columns: {missing_init_paths}"
    )

prior_chrom_cols = select_normal_prior_chromosome_columns(
    tar_chr_data,
    tar_def_path_ind_dict,
    normal_prior_excluded_chroms,
)
prior_chroms = [chrom for chrom, _ in prior_chrom_cols]
prior_cols = [col_idx for _, col_idx in prior_chrom_cols]
prior_row_capacity = len(prior_cols) if normal_prior_strength > 0 and prior_cols else 0

use_julia_solver = pipeline_mode_is_karyotype(pipeline_mode_config)

shape = (n, m + prior_row_capacity)
A_arr = np.empty(shape, dtype=np.float32, order='C')

if prior_row_capacity:
    A_arr[:, m:] = 0.0

prior_col_start = m
fm = filter_len

filter_vec_list = []
tot_loc_list = []

tmp_v = np.zeros(m, dtype=np.float32)

ncnt = 0
# Matrix column index -> canonical event tags.
# ordinary nclose: (left_contig_idx, right_contig_idx), sorted tuple of ints
# step-11 INDEL: (event_type, event_idx, type2_merge_idx)
# ecDNA columns carry their two compressed BND keys; cent/ecDNA identity tags
# remain only in path_nclose_dict_set for existing downstream path handling.
path_nclose_dict_set = defaultdict(set)
path_nclose_usage = [Counter() for _ in range(n)]
explicit_filter_reasons = {}
for path, key_int_list in tqdm(paf_ans_list, desc='Recover depth from separated paths',
                                    disable=not sys.stdout.isatty() and not args.progress):
    ki = key_int_list[0]
    np.copyto(tmp_v, vec_dict[ki])
    for ki in key_int_list[1:]:
        tmp_v += vec_dict[ki]

    idx_path = import_index_path(path)

    usage = count_index_path_events(
        idx_path, compressed_nclose_bnd_keys, type4_edge_to_event
    )
    path_nclose_usage[ncnt].update(usage)
    path_nclose_dict_set[ncnt] = set(usage)

    A_arr[ncnt, :m] = tmp_v
    ncnt += 1

with open(f'{PREFIX}/indel_exclude_idx_set.pkl', 'rb') as f:
    indel_exclude_idx_set = pkl.load( f)
tv_empty = np.zeros(len(chr_filt_st_list) + len(chr_no_filt_st_list), dtype=np.float32)

indel_idx = 0
# Process forward-directed outlier contigs
for i, ov_loc, bv_loc, type2_ins_idx, type2_ins_loc in tqdm(
              front_depth_inputs,
              desc='Parse coverage from forward-directed outlier contig gz files',
              disable=not sys.stdout.isatty() and not args.progress):
    ov = get_vec_from_stat_loc(ov_loc)
    if type2_ins_loc is not None:
        ov += get_vec_from_stat_loc(type2_ins_loc)
    bv = get_vec_from_stat_loc(bv_loc)
    
    event_key = indel_event_key('front_jump', i, type2_ins_idx)
    if event_key not in nclose_event_by_key:
        raise ValueError(f"Missing step-11 INDEL catalog event: {event_key}")
    if indel_idx in indel_exclude_idx_set:
        tv = tv_empty
        explicit_filter_reasons[event_key] = 'FILTERED_02_EXCLUDE_LIST'
    else:
        tv = ov - bv
    
    A_arr[ncnt, :m] = tv
        
    path_nclose_usage[ncnt][event_key] += 1
    path_nclose_dict_set[ncnt].add(event_key)
    ncnt += 1
    indel_idx += 1
# Process backward-directed outlier contigs
for i, ov_loc, bv_loc, type2_ins_idx, type2_ins_loc in tqdm(
              back_depth_inputs,
              desc='Parse coverage from backward-directed outlier contig gz files',
              disable=not sys.stdout.isatty() and not args.progress):
    ov = get_vec_from_stat_loc(ov_loc)
    if type2_ins_loc is not None:
        ov += get_vec_from_stat_loc(type2_ins_loc)
    bv = get_vec_from_stat_loc(bv_loc)

    event_key = indel_event_key('back_jump', i, type2_ins_idx)
    if event_key not in nclose_event_by_key:
        raise ValueError(f"Missing step-11 INDEL catalog event: {event_key}")
    if indel_idx in indel_exclude_idx_set:
        tv = tv_empty
        explicit_filter_reasons[event_key] = 'FILTERED_02_EXCLUDE_LIST'
    else:
        tv = ov + bv
    
    A_arr[ncnt, :m] = tv
        
    path_nclose_usage[ncnt][event_key] += 1
    path_nclose_dict_set[ncnt].add(event_key)
    ncnt += 1
    indel_idx += 1

for i, ov_loc in ecdna_depth_inputs:
    ov = get_vec_from_stat_loc(ov_loc)
    type2_ins_idx = -1

    A_arr[ncnt, :m] = ov
        
    ecdna_usage = count_ecdna_circuit_events(
        ecdna_circuits[i - 1], nclose_bnd_keys
    )
    path_nclose_usage[ncnt].update(ecdna_usage)
    path_nclose_dict_set[ncnt].update(ecdna_usage)
    path_nclose_dict_set[ncnt].add((ov_loc.split('/')[-2], i, type2_ins_idx))
    ncnt += 1

# Centromere depth-difference fragment chromosomes
# 02번에서 detect 된 (chrom, dir, mid_censat) 기반으로 fragment 1개 = column 1개.
# 단위 indicator vector(0/1)로 채우면 NNLS weight ≈ 추가 copy 수.
all_st_list = chr_filt_st_list + chr_no_filt_st_list
cent_fragment_cols = []
for chrom, info in cen_fragment_list:
    mid = info["mid"]
    chr_length = info["chr_len"]
    if info["dir"]:
        st_lo, st_hi = mid, chr_length
    else:
        st_lo, st_hi = 0, mid

    tv = np.zeros(m, dtype=np.float32)
    for idx, (c, st) in enumerate(all_st_list):
        if c == chrom and st_lo <= st < st_hi:
            tv[idx] = 1.0

    A_arr[ncnt, :m] = tv

    path_nclose_dict_set[ncnt].add(('cent_fragment', chrom, info["dir"]))
    cent_fragment_cols.append(ncnt)
    ncnt += 1

dep_list.extend([0] * (
    len(front_depth_inputs)
    + len(back_depth_inputs)
    + len(ecdna_depth_inputs)
    + len(cen_fragment_list)
))

assert(len(dep_list) == n)
assert ncnt == n

initial_active_columns = set(range(n))
for event_key in explicit_filter_reasons:
    for col_idx, usage in enumerate(path_nclose_usage):
        if usage.get(event_key, 0) > 0:
            initial_active_columns.discard(col_idx)

nclose_filter_status = initialise_filter_status(
    nclose_event_catalog,
    path_nclose_usage,
    initial_active_columns,
    explicit_filter_reasons,
)
save_path_usage(PREFIX, path_nclose_usage)
save_filter_status(PREFIX, nclose_filter_status)
logging.info(
    "NClose path usage : %d events across %d/%d matrix columns",
    len({key for usage in path_nclose_usage for key in usage}),
    sum(bool(usage) for usage in path_nclose_usage),
    len(path_nclose_usage),
)

for (i1, i2) in itertools.pairwise(dep_list):
    assert(i1 >= i2)

if init_cols:
    # Estimate default-chromosome targets with all cent_fragment columns present
    # as auxiliary explanatory variables.  Raw-read-excluded default paths
    # remain auxiliary only and do not receive normal-chromosome prior rows.
    prior_fit_cols = init_cols + cent_fragment_cols
    A_pri = A_arr[prior_fit_cols, :fm].T
    fit_w_pri = nnls(A_pri, B[:filter_len])[0].astype(np.float32)
    init_target_by_col = dict(zip(init_cols, fit_w_pri[:len(init_cols)]))
    w_pri = np.asarray(
        [init_target_by_col[col_idx] for col_idx in prior_cols],
        dtype=np.float32,
    )
else:
    prior_fit_cols = []
    fit_w_pri = np.asarray([], dtype=np.float32)
    w_pri = np.asarray([], dtype=np.float32)

normal_prior_row_count = 0
normal_prior_scale = 0.0
normal_prior_scales = np.asarray([], dtype=np.float32)
normal_prior_base_row_count = len(init_cols)
normal_prior_base_col_norm_squares = np.asarray([], dtype=np.float64)
normal_prior_col_norm_squares = np.asarray([], dtype=np.float64)
normal_prior_length_normalization_factor = 0.0
normal_prior_total_scale_squared = 0.0
B_prior = np.asarray([], dtype=B.dtype)

if normal_prior_strength > 0 and prior_cols:
    if len(init_cols) != len(set(init_cols)):
        raise ValueError("Default chromosome prior columns must be unique")
    if len(prior_cols) != len(set(prior_cols)):
        raise ValueError("Eligible normal chromosome prior columns must be unique")

    # The prior rows encode scale_j * w_j ~= scale_j * w_pri_j.  Normalize
    # scale_j^2 by each default column's clean-depth squared norm so the prior
    # has the same relative curvature across chromosome lengths.  The full
    # pre-exclusion init_cols set defines the denominator: excluding one
    # chromosome removes its prior mass without strengthening the survivors.
    base_depth_vectors = A_arr[np.asarray(init_cols), :fm]
    normal_prior_base_col_norm_squares = np.sum(
        np.square(base_depth_vectors, dtype=np.float64),
        axis=1,
        dtype=np.float64,
    )
    base_norm_square_by_col = dict(
        zip(init_cols, normal_prior_base_col_norm_squares)
    )
    normal_prior_col_norm_squares = np.asarray(
        [base_norm_square_by_col[col_idx] for col_idx in prior_cols],
        dtype=np.float64,
    )
    (
        normal_prior_scales_list,
        normal_prior_scale,
        normal_prior_length_normalization_factor,
    ) = calculate_length_normalized_prior_scales(
        normal_prior_base_col_norm_squares,
        normal_prior_col_norm_squares,
        fm,
        normal_prior_strength,
    )
    normal_prior_scales = np.asarray(
        normal_prior_scales_list,
        dtype=np.float32,
    )
    normal_prior_total_scale_squared = float(
        np.sum(
            np.square(normal_prior_scales, dtype=np.float64),
            dtype=np.float64,
        )
    )

    normal_prior_row_count = len(prior_cols)
    for row_idx, (col_idx, prior_scale) in enumerate(
        zip(prior_cols, normal_prior_scales)
    ):
        A_arr[col_idx, prior_col_start + row_idx] = prior_scale

    B_prior = (normal_prior_scales * w_pri).astype(B.dtype, copy=False)
    logging.info(
        "Normal chromosome prior enabled: "
        f"strength={normal_prior_strength}, rows={normal_prior_row_count}, "
        f"base_rows={normal_prior_base_row_count}, "
        f"default_rows={len(init_cols)}, excluded_rows={len(init_cols) - len(prior_cols)}, "
        f"cent_fragment_rows={len(cent_fragment_cols)}, "
        f"reference_scale={normal_prior_scale:.6g}, "
        f"scale_min={float(np.min(normal_prior_scales)):.6g}, "
        f"scale_median={float(np.median(normal_prior_scales)):.6g}, "
        f"scale_max={float(np.max(normal_prior_scales)):.6g}, "
        f"active_total_scale_squared={normal_prior_total_scale_squared:.6g}, "
        f"base_total_scale_squared={normal_prior_strength * fm:.6g}, "
        f"target_weight_median={float(np.median(w_pri)):.6g}"
    )
else:
    logging.info("Normal chromosome prior disabled.")

A_width = fm + normal_prior_row_count

with h5py.File(f'{PREFIX}/matrix.h5', 'w') as hf:
    dset_A = hf.create_dataset('A', shape=(n, A_width), dtype=A_arr.dtype)
    dset_A.write_direct(A_arr, source_sel=np.s_[:, :fm], dest_sel=np.s_[:, :fm])
    if normal_prior_row_count:
        dset_A.write_direct(
            A_arr,
            source_sel=np.s_[:, prior_col_start:prior_col_start + normal_prior_row_count],
            dest_sel=np.s_[:, fm:A_width],
        )

    dset_A_fail = hf.create_dataset('A_fail', shape=(n, m - fm), dtype=A_arr.dtype)
    dset_A_fail.write_direct(A_arr, source_sel=np.s_[:, fm:m])


    dset_B = hf.create_dataset('B', shape=(A_width,), dtype=B.dtype)
    dset_B.write_direct(B, source_sel=np.s_[:fm], dest_sel=np.s_[:fm])
    if normal_prior_row_count:
        dset_B.write_direct(B_prior, dest_sel=np.s_[fm:A_width])

    dset_B_fail = hf.create_dataset('B_fail', shape=B[fm:].shape, dtype=B.dtype)
    dset_B_fail.write_direct(B, source_sel=np.s_[fm:])

    hf.create_dataset('B_depth_start', data=0)
    hf.create_dataset('B_depth_end', data=fm)
    hf.create_dataset('normal_prior_strength', data=normal_prior_strength)
    hf.create_dataset('normal_prior_row_count', data=normal_prior_row_count)
    hf.create_dataset('normal_prior_base_row_count', data=normal_prior_base_row_count)
    # Backward-compatible scalar summary (RMS over the fixed base set).
    # The actual augmented rows use normal_prior_scales.
    hf.create_dataset('normal_prior_scale', data=normal_prior_scale)
    hf.create_dataset('normal_prior_scales', data=normal_prior_scales)
    hf.create_dataset(
        'normal_prior_total_scale_squared',
        data=normal_prior_total_scale_squared,
    )
    hf.create_dataset(
        'normal_prior_length_normalization_factor',
        data=normal_prior_length_normalization_factor,
    )

with open(f"{PREFIX}/23_input.pkl", "wb") as f:
    pkl.dump((
        chr_filt_st_list,
        path_nclose_dict_set,
        amplitude,
        bed_data,
        {
            "B_depth_start": 0,
            "B_depth_end": fm,
        },
    ), f)

with open(f"{PREFIX}/tar_chr_data.pkl", "wb") as f:
    pkl.dump(tar_chr_data, f)

with open(f"{PREFIX}/normal_prior_data.pkl", "wb") as f:
    pkl.dump({
        "strength": normal_prior_strength,
        "normalization": "column_norm_fixed_base_total",
        # Backward-compatible scalar summary; use "scales" for actual rows.
        "scale": normal_prior_scale,
        "scales": normal_prior_scales,
        "total_scale_squared": normal_prior_total_scale_squared,
        "row_count": normal_prior_row_count,
        "base_row_count": normal_prior_base_row_count,
        "base_col_norm_squares": normal_prior_base_col_norm_squares,
        "prior_col_norm_squares": normal_prior_col_norm_squares,
        "length_normalization_factor": normal_prior_length_normalization_factor,
        "prior_chroms": prior_chroms,
        "prior_cols": prior_cols,
        "prior_fit_cols": prior_fit_cols,
        "init_cols": init_cols,
        "cent_fragment_cols": cent_fragment_cols,
        "excluded_chroms": sorted(normal_prior_excluded_chroms),
        "targets": w_pri,
        "fit_targets": fit_w_pri,
        "tar_chr_data": tar_chr_data,
    }, f)

np.save(f'{PREFIX}/B.npy', B)
