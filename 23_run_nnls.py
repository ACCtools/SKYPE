import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *

import h5py
import glob
import psutil
import logging
import warnings
import argparse

import pickle as pkl
import numpy as np
import pandas as pd

from collections import defaultdict, Counter
from scipy.stats import mannwhitneyu, ttest_ind

from skglm import GeneralizedLinearEstimator
from skglm.datafits import Quadratic
from skglm.penalties import PositiveConstraint
from skglm.solvers import AndersonCD

from h5py._hl import selections as sel

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

FILTER_P_VALUE = 0.01

CHROM_ERROR_FAIL_RATE = 2.0

def normalize_path_nclose_tag(tag):
    if isinstance(tag, (tuple, list)):
        if len(tag) == 2:
            return tuple(sorted(tag))
        if len(tag) == 3 and tag[0] in {'front_jump', 'back_jump', 'ecdna'}:
            event_type, event_idx, type2_merge_idx = tag
            return (event_type, int(event_idx), int(type2_merge_idx))
        if len(tag) == 3 and tag[0] == 'cent_fragment':
            return tuple(tag)
        
    raise ValueError(f'Unexpected nclose tag format: {tag!r}')

def normalize_path_nclose_set_dict(path_nclose_set_dict):
    return {
        k: {normalize_path_nclose_tag(tag) for tag in path_nclose_list}
        for k, path_nclose_list in path_nclose_set_dict.items()
    }

def import_ppc_contig_data(file_path: str) -> list:
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

# Hypothesis test nclose filtering

def get_contig_pair(i, cd):
    nclose_loc = i == 0
    nclose_dir = '+' == cd[CTG_DIR]

    chr_name = cd[CHR_NAM]
    chr_cord = cd[CHR_STR if nclose_dir ^ nclose_loc else CHR_END]

    return chr_name, chr_cord

def get_chr_cord_data(contig_data, nclose_nodes, type4_nclose_nodes):
    chr_cord_dict = defaultdict(list)
    ncloseidx2cord = defaultdict(list)
    ncloseidx2path_tag = {}
    
    for ctg_name, nclose_idx_pair_list in nclose_nodes.items():
        for j, nclose_idx_pair in enumerate(nclose_idx_pair_list):
            nclose_key = (ctg_name, j)
            ncloseidx2path_tag[nclose_key] = normalize_path_nclose_tag(nclose_idx_pair)
            for i, idx in enumerate(nclose_idx_pair):
                chr_name, chr_cord = get_contig_pair(i, contig_data[idx])

                chr_cord_dict[chr_name].append((chr_cord, nclose_key))
                ncloseidx2cord[nclose_key].append((chr_name, chr_cord))
                

    for type4_nclose_idx, (c1, i1, c2, i2) in type4_nclose_nodes.items():
        type4_nclose_idx = normalize_path_nclose_tag(type4_nclose_idx)

        chr_cord_dict[c1].append((i1, type4_nclose_idx))
        chr_cord_dict[c2].append((i2, type4_nclose_idx))
        ncloseidx2cord[type4_nclose_idx].append((c1, i1))
        ncloseidx2cord[type4_nclose_idx].append((c2, i2))
        ncloseidx2path_tag[type4_nclose_idx] = type4_nclose_idx


    for l in chr_cord_dict.values():
        l.sort(key=lambda t: t[0])

    return chr_cord_dict, ncloseidx2cord, ncloseidx2path_tag

def group_close_coordinates(data_list, distance_threshold=50 * K):
    """
    Groups a list of tuples by their coordinate values.
    
    This function assumes the input data_list is already sorted by the 
    coordinate (the first element of the tuple). Tuples are grouped together 
    if the distance between consecutive coordinates is less than or equal 
    to the distance_threshold (default is 100k).
    
    Args:
        data_list: A list of tuples formatted as (cord, ctg_name, nclose_idx).
        distance_threshold: Maximum distance between consecutive coordinates to be grouped.
        
    Returns:
        A list of grouped tuples formatted as (group_start_cord, [(ctg_name, nclose_idx), ...]).
    """
    if not data_list:
        return []

    result = []
    
    # Initialize the first group using the first element in the list
    first_cord = data_list[0][0]
    last_cord = data_list[0][0]
    current_group_items = [data_list[0][1]]

    for i in range(1, len(data_list)):
        cord, nclose_data = data_list[i]

        # Check if the current coordinate is within the threshold from the previous one
        if cord - last_cord <= distance_threshold:
            current_group_items.append(nclose_data)
        else:
            # Save the completed group and start a new one
            result.append((first_cord, current_group_items))
            first_cord = cord
            current_group_items = [nclose_data]
            
        # Update the last seen coordinate for the next iteration
        last_cord = cord

    # Append the final group after the loop finishes
    result.append((first_cord, current_group_items))

    return result

def merge_regions_iteratively(df_cov, df_exclude,
                              target_chr, chr_cord_dict,
                              p_value):
    """
    Iteratively merges adjacent regions based on t-test p-values.
    It prevents merging only if the specific split point falls within any excluded areas.
    """
    split_points = group_close_coordinates(chr_cord_dict[target_chr])
    out_split_points = []
    
    # Filter dataframes for the specific chromosome to speed up operations
    df_cov_chr = df_cov[df_cov['chr'] == target_chr].copy()
    df_excl_chr = df_exclude[df_exclude['chr'] == target_chr]
    
    # Filter out the excluded regions from the coverage data upfront
    # Bins that overlap with any excluded region are completely removed
    for _, excl_row in df_excl_chr.iterrows():
        df_cov_chr = df_cov_chr[
            (df_cov_chr['nd'] <= excl_row['st']) | 
            (df_cov_chr['st'] >= excl_row['nd'])
        ]
        
    # Copy the split points to safely modify the active list during iteration
    active_splits = split_points.copy()
    
    # Determine the absolute maximum boundary from the remaining valid coverage data
    max_nd = df_cov_chr['nd'].max() if not df_cov_chr.empty else 0
    
    while True:
        merged_in_this_round = False
        
        for i in range(len(active_splits)):
            current_split = active_splits[i][0]
            
            # Check if the specific split point itself falls within any excluded region
            # Skip the merge process if the exact point is located inside an excluded area
            is_split_point_excluded = (
                (df_excl_chr['st'] <= current_split) & 
                (df_excl_chr['nd'] >= current_split)
            ).any()
            
            if is_split_point_excluded:
                continue
            
            # Define boundaries for the left and right regions
            left_start = active_splits[i-1][0] if i > 0 else 0
            left_end = current_split
            
            right_start = current_split
            right_end = active_splits[i+1][0] if i < len(active_splits) - 1 else max_nd
            
            # Extract coverage values for left and right regions based on bin coordinates
            # Since df_cov_chr has already been filtered, excluded regions naturally won't be included
            left_cov = list(df_cov_chr[
                (df_cov_chr['nd'] > left_start) & 
                (df_cov_chr['st'] < left_end)
            ]['meandepth'])
            
            right_cov = list(df_cov_chr[
                (df_cov_chr['nd'] > right_start) & 
                (df_cov_chr['st'] < right_end)
            ]['meandepth'])
            
            # Ensure both regions have at least 2 data points to perform the statistical test
            # This naturally avoids statistical calculation errors
            if len(left_cov) >= 5 and len(right_cov) >= 5:
                # stat, p_val = ttest_ind(left_cov, right_cov)
                stat, p_val = mannwhitneyu(left_cov, right_cov)

                # If distributions are not significantly different, merge them
                if p_val >= p_value and (not np.isnan(p_val)):
                    out_split_points.append(active_splits[i])
                    active_splits.pop(i)
                    merged_in_this_round = True
                    break # Break the inner loop to restart the while loop with updated regions
                    
        # If no merges occurred in the entire pass, the regions are stable
        if not merged_in_this_round:
            break
            
    return active_splits, out_split_points

def get_censet_count(cdf, c1, i1, c2, i2):
    cnt = 0
    for chr_name, chr_cord in [(c1, i1), (c2, i2)]:
        mask = (cdf['chr'] == chr_name) & (cdf['st'] <= chr_cord) & (cdf['nd'] >= chr_cord)
        if mask.any():
            cnt += 1

    return cnt

def filter_nclose_by_test(contig_data, nclose_nodes, df, cdf):
    type4_tot_loc_list = []
    for bv_paf_loc in glob.glob(f'{RATIO_OUTLIER_FOLDER}/front_jump/*_base.paf'):
        type4_tot_loc_list.append(bv_paf_loc)

    for bv_paf_loc in glob.glob(f'{RATIO_OUTLIER_FOLDER}/back_jump/*_base.paf'):
        type4_tot_loc_list.append(bv_paf_loc)

    type4_nclose_nodes = dict()
    for paf_loc in type4_tot_loc_list:
        event_type = paf_loc.split('/')[-2]
        event_idx = int(paf_loc.split('/')[-1].split('_')[0])

        type2_merge_paf_idx = -1
        if not os.path.isfile(f'{RATIO_OUTLIER_FOLDER}/{event_type}/{event_idx}.paf'):
            type2_merge_paf_loc = glob.glob(f'{RATIO_OUTLIER_FOLDER}/{event_type}/{event_idx}_type2_merge_*.paf')[0]
            type2_merge_paf_idx = int(type2_merge_paf_loc.split('/')[-1].split('.')[0].split('_')[-1])

        with open(paf_loc, "r") as f:
            l = f.readline().rstrip().split("\t")
            chr_nam1 = l[CHR_NAM]
            chr_nam2 = l[CHR_NAM]
            pos1 = int(l[CHR_STR])
            pos2 = int(l[CHR_END])

            type4_tag = normalize_path_nclose_tag((event_type, event_idx, type2_merge_paf_idx))
            type4_nclose_nodes[type4_tag] = (chr_nam1, pos1, chr_nam2, pos2)

    chr_cord_dict, ncloseidx2cord, ncloseidx2path_tag = get_chr_cord_data(contig_data, nclose_nodes, type4_nclose_nodes)

    in_cnt = defaultdict(lambda : [0, 0])
    for chr_key in chr_cord_dict:
        if chr_key in {'chr1', 'chr18'}:
            t = 1

        in_split, out_split = merge_regions_iteratively(df, cdf, chr_key, chr_cord_dict, p_value=FILTER_P_VALUE)

        for ctg_idx, nclose_idx_list in in_split:
            for nclose_idx in nclose_idx_list:
                in_cnt[nclose_idx][0] += 1

        for ctg_idx, nclose_idx_list in out_split:
            for nclose_idx in nclose_idx_list:
                in_cnt[nclose_idx][1] += 1

    return in_cnt, ncloseidx2cord, ncloseidx2path_tag

def _group_consecutive(sorted_idx_list):
    """
    Group a sorted list of indices into (start, end) ranges
    where end is inclusive.
    e.g. [0,1,2,5,6,9] -> [(0,2), (5,6), (9,9)]
    """
    if not sorted_idx_list:
        return []
    ranges = []
    start = sorted_idx_list[0]
    prev = start
    for idx in sorted_idx_list[1:]:
        if idx == prev + 1:
            prev = idx
        else:
            ranges.append((start, prev))
            start = idx
            prev = idx
    ranges.append((start, prev))
    return ranges


def read_selected_columns_direct(dA, A, A_idx_list):
    """
    Read specific columns from an HDF5 dataset directly into a slice of a
    pre-allocated NumPy array using the low-level HDF5 API.

    A_idx_list must be sorted. Consecutive indices are batched into
    contiguous slice reads for much better I/O performance.
    """
    ranges = _group_consecutive(A_idx_list)

    dest_col = 0
    for start, end in ranges:
        ncols = end - start + 1

        # Source: contiguous column slice in the file
        source_sel = sel.select(dA.shape, np.s_[:, start:end + 1], dA)
        fspace = source_sel.id

        # Destination: corresponding slice in the output array
        dest_sel = sel.select(A.shape, np.s_[:, dest_col:dest_col + ncols])

        for mspace in dest_sel.broadcast(source_sel.array_shape):
            dA.id.read(mspace, fspace, A, dxpl=dA._dxpl)

        dest_col += ncols

warnings.simplefilter(action='ignore', category=FutureWarning)

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).
logging.info("23_run_nnls start")

parser = argparse.ArgumentParser(description="SKYPE depth analysis")
parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")

parser.add_argument("prefix", help="Prefix for pipeline")

parser.add_argument("main_stat_path", help="Path to the main depth statistics file")

parser.add_argument("censat_bed_path", 
                    help="Path to the censat repeat information file.")

parser.add_argument("-t", "--thread", help="Number of threads", type=int)

# args = parser.parse_args()

t = "23_run_nnls.py /Data/hyunwoo/00_skype_run_data/HG008-T/20_alignasm/HG008-T.ctg.aln.paf.ppc.paf /Data/hyunwoo/00_skype_run_data/HG008-T/30_skype /Data/hyunwoo/00_skype_run_data/HG008-T/01_depth/HG008-T_normalized.win.stat.gz /hyunwoo/63_skype_test/deps/SKYPE/public_data/chm13v2.0_censat_v2.1.m.bed -t 72"
args = parser.parse_args(t.split()[1:])

PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
PREFIX = args.prefix
main_stat_loc = args.main_stat_path
THREAD = args.thread
CENSAT_PATH = args.censat_bed_path

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"

with open(f"{PREFIX}/report.txt", 'r') as f:
    f.readline()
    path_cnt = int(f.readline().strip())

use_julia_solver = path_cnt <= HARD_PATH_COUNT_BASELINE

core_num = psutil.cpu_count(logical=False)
if core_num is None:
    THREAD = THREAD
else:
    THREAD = min(int(THREAD), core_num)

contig_data = import_ppc_contig_data(PREPROCESSED_PAF_FILE_PATH)

with open(f'{PREFIX}/nclose_chunk_data.pkl', 'rb') as f:
    nclose_nodes, _, _ = pkl.load(f)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

cdf = pd.read_csv(CENSAT_PATH, sep='\t', names=['chr', 'st', 'nd'])

with open(f'{PREFIX}/23_input.pkl', 'rb') as f:
    chr_filt_st_list, path_nclose_set_dict, amplitude, bed_data = pkl.load(f)
path_nclose_set_dict = normalize_path_nclose_set_dict(path_nclose_set_dict)

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

with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
    dA = f["A"]
    A = np.empty(dA.shape, dtype=dA.dtype)
    dA.read_direct(A)

    dB = f["B"]
    B = np.empty(dB.shape, dtype=dB.dtype)
    dB.read_direct(B)

    dAf = f["A_fail"]
    A_fail = np.empty(dAf.shape, dtype=dAf.dtype)
    dAf.read_direct(A_fail)

    b_start_ind = int(f["B_depth_start"][()])

nnls = GeneralizedLinearEstimator(
    datafit=Quadratic(),
    penalty=PositiveConstraint(),
    solver=AndersonCD(fit_intercept=False)
)

nnls.fit(A, B)
weights_fullsize = nnls.coef_
predict_suc_B_base = A.dot(weights_fullsize)

chrom_acc_sum_dict_base = defaultdict(int)
chrom_acc_sum_dict_max_base = defaultdict(int)
for i, (chrom, st) in enumerate(chr_filt_st_list):
    chrom_acc_sum_dict_base[chrom] += predict_suc_B_base[i] - B[i]
    if abs(chrom_acc_sum_dict_base[chrom]) > chrom_acc_sum_dict_max_base[chrom]:
        chrom_acc_sum_dict_max_base[chrom] = abs(chrom_acc_sum_dict_base[chrom])

# Loop start
prev_fail_chrom_list = None
fail_chrom_list = []

while True:
    in_cnt, ncloseidx2cord, ncloseidx2path_tag = filter_nclose_by_test(contig_data, nclose_nodes, df, cdf)

    not_using_nclose_set = set()

    for nclose_key, nclose_cnt in in_cnt.items():
        if nclose_key == ('utg008527l', 0):
            q = 1
        (c1, i1), (c2, i2) = ncloseidx2cord[nclose_key]
        censat_count = get_censet_count(cdf, c1, i1, c2, i2)
        if (not (censat_count == 2 or nclose_cnt[0] - censat_count >= 1) and
            c1 not in fail_chrom_list and c2 not in fail_chrom_list):
            not_using_nclose_set.add(ncloseidx2path_tag[nclose_key])

    fail_chrom_text = ','.join(fail_chrom_list) if fail_chrom_list else 'none'
    logging.info(f'Hyp. test filtering nclose count : {len(not_using_nclose_set)} (fail chroms: {fail_chrom_text})')
    logging.debug(
        f'Hyp. test filtering nclose breakdown : '
        f'{dict(Counter("ordinary" if len(tag) == 2 else tag[0] for tag in not_using_nclose_set))}'
    )

    A_idx_list = []
    for k, path_nclose_list in path_nclose_set_dict.items():
        if not path_nclose_list or len(not_using_nclose_set & set(path_nclose_list)) == 0:
            A_idx_list.append(k)

    A_idx_list.sort()

    with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
        dA = f["A"]
        dAf = f["A_fail"]

        read_selected_columns_direct(dA, A, A_idx_list)
        read_selected_columns_direct(dAf, A_fail, A_idx_list)

    A_new = A[:, :len(A_idx_list)]
    A_fail_new = A_fail[:, :len(A_idx_list)]

    nnls.fit(A_new, B)
    filter_weights = nnls.coef_

    predict_suc_B = A_new.dot(filter_weights)
    predict_fal_B = A_fail_new.dot(filter_weights)

    chrom_acc_sum_dict = defaultdict(int)
    chrom_acc_sum_dict_max = defaultdict(int)
    for i, (chrom, st) in enumerate(chr_filt_st_list):
        chrom_acc_sum_dict[chrom] += predict_suc_B[i] - predict_suc_B_base[i]
        if abs(chrom_acc_sum_dict[chrom]) > chrom_acc_sum_dict_max[chrom]:
            chrom_acc_sum_dict_max[chrom] = abs(chrom_acc_sum_dict[chrom])

    new_fail_chrom_list = []
    for chrom, acc_sum_max in chrom_acc_sum_dict_max.items():
        if chrom == 'chrY' and no_chrY:
            continue
        acc_sum_max_base = chrom_acc_sum_dict_max_base[chrom]
        chrom_error_rate = acc_sum_max / acc_sum_max_base
        if chrom_error_rate > CHROM_ERROR_FAIL_RATE:
            new_fail_chrom_list.append(chrom)

    # 새로 발견된 fail chrom을 누적
    fail_chrom_list = sorted(set(fail_chrom_list + new_fail_chrom_list))

    # Break: 새로 추가된 fail chrom이 없으면 종료
    if not new_fail_chrom_list or (prev_fail_chrom_list is not None and fail_chrom_list == prev_fail_chrom_list):
        break

    prev_fail_chrom_list = fail_chrom_list[:]

# === Cent_fragment greedy off-test ===
# nclose 필터링 직후를 fail 기준점(고정 baseline)으로 잡고, weight 작은 cent_fragment
# 컬럼부터 누적 제거 시도. PASS 시 컬럼 누적 제거(=offset 누적)되며, fail 비교 대상
# baseline 은 항상 nclose 필터링 직후로 고정 — 누적 drift 도 base 대비 한도 안이면 OK.
cent_fragment_col2chrom = {}
for k, tags in path_nclose_set_dict.items():
    for tag in tags:
        if isinstance(tag, tuple) and len(tag) > 0 and tag[0] == 'cent_fragment':
            cent_fragment_col2chrom[k] = tag[1]
            break

def _compute_acc_max(predict_a, predict_b):
    acc = defaultdict(float)
    acc_max = defaultdict(float)
    for i, (chrom, _) in enumerate(chr_filt_st_list):
        acc[chrom] += predict_a[i] - predict_b[i]
        if abs(acc[chrom]) > acc_max[chrom]:
            acc_max[chrom] = abs(acc[chrom])
    return acc_max

# 고정 baseline (nclose-filter 직후) — 누적 제거 후 비교 대상
predict_suc_B_post = predict_suc_B.copy()
cf_acc_max_base_fixed = _compute_acc_max(predict_suc_B_post, B)

# 누적 상태 (PASS마다 갱신)
A_idx_list_curr = sorted(A_idx_list)
predict_suc_B_curr = predict_suc_B.copy()
predict_fal_B_curr = predict_fal_B.copy()
filter_weights_curr = filter_weights.copy()

logging.info(f'CF greedy start: pool={len(cent_fragment_col2chrom)} cent_fragment cols')

removed_cf = []
while True:
    break
    pos_of = {col: i for i, col in enumerate(A_idx_list_curr)}
    pool = [(filter_weights_curr[pos_of[col]], col, chrom)
            for col, chrom in cent_fragment_col2chrom.items() if col in pos_of]
    if not pool:
        break
    pool.sort(key=lambda x: x[0])

    progressed = False
    for w, col, chrom in pool:
        A_idx_test = sorted(set(A_idx_list_curr) - {col})
        with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
            read_selected_columns_direct(f["A"], A, A_idx_test)
            read_selected_columns_direct(f["A_fail"], A_fail, A_idx_test)
        
        A_test = A[:, :len(A_idx_test)]
        A_fail_test = A_fail[:, :len(A_idx_test)]

        nnls.fit(A_test, B)
        weights_test = nnls.coef_
        predict_suc_B_test = A_test.dot(weights_test)
        predict_fal_B_test = A_fail_test.dot(weights_test)

        # numerator: 누적 제거된 test prediction vs FIXED baseline (nclose-filter 직후)
        acc_max_diff = _compute_acc_max(predict_suc_B_test, predict_suc_B_post)

        ok = True
        worst_chrom = None
        worst_ratio = 0.0
        for chk, max_diff in acc_max_diff.items():
            if chk == 'chrY' and no_chrY:
                continue
            base_max = cf_acc_max_base_fixed.get(chk, 0.0)
            if base_max == 0:
                continue
            ratio = max_diff / base_max
            if ratio > CHROM_ERROR_FAIL_RATE:
                ok = False
                if ratio > worst_ratio:
                    worst_ratio = ratio
                    worst_chrom = chk

        if ok:
            logging.debug(f'CF PASS: drop col {col} chrom={chrom} w={w:.4f}')
            A_idx_list_curr = A_idx_test
            predict_suc_B_curr = predict_suc_B_test
            predict_fal_B_curr = predict_fal_B_test
            filter_weights_curr = weights_test
            removed_cf.append((col, chrom))
            progressed = True
            break
        else:
            logging.debug(f'CF FAIL: keep col {col} chrom={chrom} w={w:.4f} '
                         f'(worst chrom={worst_chrom} ratio={worst_ratio:.3f})')

    if not progressed:
        break

logging.info(f'CF greedy done. Removed {len(removed_cf)} cols: {removed_cf}')

A_idx_list = A_idx_list_curr
filter_weights = filter_weights_curr
predict_suc_B = predict_suc_B_curr
predict_fal_B = predict_fal_B_curr

# === BP step / depth ratio 기반 nclose 추가 필터링 ===
# predict_suc_B와 B 각각에서, 정상 nclose (type 4 제외) 의 두 BP 점 양옆 ±500kb
# 절대 step 을 (그 nclose 의 count*weight 합) 으로 나눈 ratio 가 모두 10% 미만이면 drop set에 추가.
# 최종 제거 대상은 predict_suc_B 기반 drop set과 B 기반 drop set의 union.
# BP 가 censat 영역에 있거나 ±500kb window 못 만들면 false 처리 (=제거 안 함).
# cent_fragment greedy off-test 가 끝나 weight 가 정리된 상태에서 적용 (cent_fragment 의 잔여 영향 배제).
BP_STEP_WIN = 500 * K
BP_STEP_RATIO_THRESHOLD = 0.25

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    bp_path_list_dict = pkl.load(f)
with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    bp_paf_sort_ans_list, _, _, _ = pkl.load(f)

bp_nclose_set = set()
for ctg, pairs in nclose_nodes.items():
    for pair in pairs:
        bp_nclose_set.add(tuple(pair))

# 각 path column 의 nclose 사용 Counter (31번 path_nclose_usage 와 동일)
bp_path_nclose_count = []
for col_idx in range(len(bp_paf_sort_ans_list)):
    paf_loc = bp_paf_sort_ans_list[col_idx][0]
    key = paf_loc.split('/')[-2]
    cnt = int(paf_loc.split('/')[-1].split('.')[0]) - 1
    idx_path = bp_path_list_dict[key][cnt][0]
    ctr = Counter()
    s = 1
    while s < len(idx_path) - 2:
        cand = tuple(sorted([idx_path[s][1], idx_path[s+1][1]]))
        if cand in bp_nclose_set:
            ctr[cand] += 1
            s += 2
        else:
            s += 1
    bp_path_nclose_count.append(ctr)

# 각 nclose 의 두 BP (chr, cord)
bp_nclose_cords = {}
for ctg, pairs in nclose_nodes.items():
    for pair in pairs:
        bp_nclose_cords[tuple(pair)] = [get_contig_pair(i, contig_data[idx]) for i, idx in enumerate(pair)]

# 현재 filter_weights 를 전체 column 길이로 unfold
bp_weight_full = np.zeros(len(bp_paf_sort_ans_list), dtype=np.float64)
for j, col in enumerate(A_idx_list):
    if col < len(bp_paf_sort_ans_list):
        bp_weight_full[col] = filter_weights[j]

# nclose 별 depth = Σ (count × weight)
bp_nclose_depth = defaultdict(float)
for col_idx, ctr in enumerate(bp_path_nclose_count):
    for nc, v in ctr.items():
        bp_nclose_depth[nc] += v * bp_weight_full[col_idx]

# chr -> sorted [(bin_idx, st), ...]
bp_chr_to_bins = defaultdict(list)
for i, (chrom, st) in enumerate(chr_filt_st_list):
    bp_chr_to_bins[chrom].append((i, st))
for chrom in bp_chr_to_bins:
    bp_chr_to_bins[chrom].sort(key=lambda x: x[1])

def bp_in_censat(chrom, cord):
    mask = (cdf['chr'] == chrom) & (cdf['st'] <= cord) & (cdf['nd'] >= cord)
    return bool(mask.any())

def bp_signal_passes_ratio(signal_B, left_bins, right_bins, depth_val):
    left_rows = [b_start_ind + i for i in left_bins]
    right_rows = [b_start_ind + i for i in right_bins]
    step_abs = abs(signal_B[right_rows].mean() - signal_B[left_rows].mean())
    return (step_abs / depth_val) < BP_STEP_RATIO_THRESHOLD

def bp_passes_ratio(chrom, bp, depth_val, signal_B):
    if bp_in_censat(chrom, bp):
        return False
    bins = bp_chr_to_bins.get(chrom, [])
    left = [i for i, st in bins if bp - BP_STEP_WIN <= st < bp]
    right = [i for i, st in bins if bp <= st < bp + BP_STEP_WIN]
    if not left or not right:
        return False
    return bp_signal_passes_ratio(signal_B, left, right, depth_val)

def bp_find_nclose_to_drop(signal_B):
    nclose_to_drop = set()
    for nc, cords in bp_nclose_cords.items():
        if nc == (1921, 1922):
            t = 1
        depth = bp_nclose_depth.get(nc, 0.0)
        if depth < 1e-6:
            continue
        if all(bp_passes_ratio(chrom, pos, depth, signal_B) for chrom, pos in cords):
            nclose_to_drop.add(nc)
    return nclose_to_drop

bp_nclose_to_drop_predict = bp_find_nclose_to_drop(predict_suc_B)
bp_nclose_to_drop_observed = bp_find_nclose_to_drop(B)
bp_nclose_to_drop = bp_nclose_to_drop_predict | bp_nclose_to_drop_observed

logging.info(
    f'BP step/depth filter : nclose to remove = {len(bp_nclose_to_drop)} '
    f'(predict={len(bp_nclose_to_drop_predict)}, observed={len(bp_nclose_to_drop_observed)})'
)

bp_drop_cols = set()
for col_idx, ctr in enumerate(bp_path_nclose_count):
    if any(nc in bp_nclose_to_drop for nc in ctr):
        bp_drop_cols.add(col_idx)

if bp_drop_cols:
    bp_new_A_idx_list = sorted(c for c in A_idx_list if c not in bp_drop_cols)
    with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
        dA = f["A"]
        dAf = f["A_fail"]
        read_selected_columns_direct(dA, A, bp_new_A_idx_list)
        read_selected_columns_direct(dAf, A_fail, bp_new_A_idx_list)
    
    A_new = A[:, :len(bp_new_A_idx_list)]
    A_fail_new = A_fail[:, :len(bp_new_A_idx_list)]
    nnls.fit(A_new, B)
    filter_weights = nnls.coef_
    predict_suc_B = A_new.dot(filter_weights)
    predict_fal_B = A_fail_new.dot(filter_weights)
    A_idx_list = bp_new_A_idx_list
    logging.info(f'BP step/depth filter : columns kept {len(A_idx_list)} after refit')
else:
    logging.info('BP step/depth filter : no nclose removed, skip refit')

final_weights_fullsize = np.zeros_like(weights_fullsize)
final_weights_fullsize[A_idx_list] = filter_weights

b_norm = np.linalg.norm(B)

error = np.linalg.norm(predict_suc_B - B)
predict_B = np.concatenate((predict_suc_B, predict_fal_B))[b_start_ind:]

logging.info(f'Error : {error:.4f}')
logging.info(f'Relative error : {error/b_norm:.4f}')

np.save(f'{PREFIX}/weight.npy', final_weights_fullsize)
np.save(f'{PREFIX}/predict_B.npy', predict_B)

with open(f'{PREFIX}/A_idx_list.pkl', 'wb') as f:
    pkl.dump(A_idx_list, f)
