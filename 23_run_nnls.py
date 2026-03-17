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

from collections import defaultdict
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

CHROM_ERROR_FAIL_RATE = 1.0

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
    
    for ctg_name, nclose_idx_pair_list in nclose_nodes.items():
        for j, nclose_idx_pair in enumerate(nclose_idx_pair_list):
            for i, idx in enumerate(nclose_idx_pair):
                chr_name, chr_cord = get_contig_pair(i, contig_data[idx])

                chr_cord_dict[chr_name].append((chr_cord, (ctg_name, j)))
                ncloseidx2cord[(ctg_name, j)].append((chr_name, chr_cord))
                

    for type4_nclose_idx, (c1, i1, c2, i2) in type4_nclose_nodes.items():
        chr_cord_dict[c1].append((i1, type4_nclose_idx))
        chr_cord_dict[c2].append((i2, type4_nclose_idx))

        ncloseidx2cord[type4_nclose_idx].append((c1, i1))
        ncloseidx2cord[type4_nclose_idx].append((c2, i2))


    for l in chr_cord_dict.values():
        l.sort(key=lambda t: t[0])

    return chr_cord_dict, ncloseidx2cord

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

            type4_nclose_nodes[(event_type, event_idx, type2_merge_paf_idx)] = (chr_nam1, pos1, chr_nam2, pos2)

    chr_cord_dict, ncloseidx2cord = get_chr_cord_data(contig_data, nclose_nodes, type4_nclose_nodes)

    in_cnt = defaultdict(lambda : [0, 0])
    for chr_key in chr_cord_dict:
        in_split, out_split = merge_regions_iteratively(df, cdf, chr_key, chr_cord_dict, p_value=FILTER_P_VALUE)

        for ctg_idx, nclose_idx_list in in_split:
            for nclose_idx in nclose_idx_list:
                in_cnt[nclose_idx][0] += 1

        for ctg_idx, nclose_idx_list in out_split:
            for nclose_idx in nclose_idx_list:
                in_cnt[nclose_idx][1] += 1

    return in_cnt, ncloseidx2cord

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

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("23_run_nnls start")

parser = argparse.ArgumentParser(description="SKYPE depth analysis")
parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")

parser.add_argument("prefix", help="Prefix for pipeline")

parser.add_argument("main_stat_path", help="Path to the main depth statistics file")

parser.add_argument("censat_bed_path", 
                    help="Path to the censat repeat information file.")

parser.add_argument("-t", "--thread", help="Number of threads", type=int)

args = parser.parse_args()

# t = "23_run_nnls.py /home/hyunwoo/ACCtools-pipeline/90_skype_run/Caki-1/20_alignasm/Caki-1.ctg.aln.paf.ppc.paf 30_skype_pipe/Caki-1_13_47_01 /home/hyunwoo/ACCtools-pipeline/90_skype_run/Caki-1/01_depth/Caki-1_normalized.win.stat.gz public_data/chm13v2.0_censat_v2.1.m.bed -t 128"
# args = parser.parse_args(t.split()[1:])

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
    in_cnt, ncloseidx2cord = filter_nclose_by_test(contig_data, nclose_nodes, df, cdf)

    not_using_nclose_set = set()

    for nclose_key, nclose_cnt in in_cnt.items():
        (c1, i1), (c2, i2) = ncloseidx2cord[nclose_key]
        censat_count = get_censet_count(cdf, c1, i1, c2, i2)
        if (not (censat_count == 2 or nclose_cnt[0] - censat_count >= 1) and
            c1 not in fail_chrom_list and c2 not in fail_chrom_list):
            not_using_nclose_set.add(nclose_key)

    logging.info(f'Hyp. test filtering nclose count : {len(not_using_nclose_set)}')

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

final_weights_fullsize = np.zeros_like(weights_fullsize)
final_weights_fullsize[A_idx_list] = filter_weights

b_norm = np.linalg.norm(B)

error = np.linalg.norm(predict_suc_B - B)
predict_B = np.concatenate((predict_suc_B, predict_fal_B))[b_start_ind:]

logging.info(f'Error : {error:.4f}')
logging.info(f'Relative error : {error/b_norm:.4f}')

np.save(f'{PREFIX}/weight.npy', final_weights_fullsize)
np.save(f'{PREFIX}/predict_B.npy', predict_B)

with open('A_idx_list.pkl', 'wb') as f:
    pkl.dump(A_idx_list, f)

