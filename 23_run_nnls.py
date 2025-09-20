import os
import psutil
import logging
import warnings
import argparse

import pickle as pkl
import numpy as np
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

K = 1000
M = 1000 * K

DEPTH_VECTOR_WINDOW = 100 * K
DEPTH_ERROR_ALLOW_LEN = 500 * K

BASE_ACCSUMABSMAX_RATIO = 0.1

VCF_FLANKING_LENGTH = 1*M
NCLOSE_SIM_COMPARE_RAITO = 1.2
NCLOSE_SIM_DIFF_THRESHOLD = 5

HARD_PATH_COUNT_BASELINE = 100 * K

def similar_check(v1, v2, ratio):
    try:
        assert(v1 >= 0 and v2 >= 0)
    except:
        logging.error(f"Invalid values for similarity check: v1={v1}, v2={v2}")
    mi, ma = sorted([v1, v2])
    return True if mi == 0 else (ma / mi <= ratio) or ma-mi < NCLOSE_SIM_DIFF_THRESHOLD

def check_near_bnd(chrom, inside_st, inside_nd):
    # subset of df for the given chromosome
    df_chr = df[df['chr'] == chrom]

    def mean_depth(start, end):
        """Return mean meandepth over windows overlapping [start, end)."""
        mask = (df_chr['nd'] > start) & (df_chr['st'] < end)
        return df_chr.loc[mask, 'meandepth'].mean()

    # for inside_st
    st_depth = mean_depth(inside_st - VCF_FLANKING_LENGTH, inside_st)
    nd_depth = mean_depth(inside_nd, inside_nd + VCF_FLANKING_LENGTH)
    if np.isnan(st_depth) or np.isnan(nd_depth):
        return True
    
    return similar_check(st_depth, nd_depth, NCLOSE_SIM_COMPARE_RAITO)

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
parser.add_argument("-t", "--thread", help="Number of threads", type=int)
args = parser.parse_args()

PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
PREFIX = args.prefix
main_stat_loc = args.main_stat_path
THREAD = args.thread

core_num = psutil.cpu_count(logical=False)
if core_num is None:
    THREAD = THREAD
else:
    THREAD = min(int(THREAD), core_num)

ppc_contig_data = import_ppc_contig_data(PREPROCESSED_PAF_FILE_PATH)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t',
                 names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')
meandepth = df['meandepth'].median()

ydf = df.query('chr == "chrY"')
chry_nz_len = len(ydf.query('meandepth != 0'))
no_chrY = (chry_nz_len / len(ydf)) < 0.5

with open(f'{PREFIX}/23_input.pkl', 'rb') as f:
#     dep_list, init_cols, w_pri, using_ncnt_array, \
#     div_nclose_set, each_nclose_notusing_path_dict, \
    chr_filt_st_list, path_nclose_set_dict, amplitude = pkl.load(f)

with open(f"{PREFIX}/report.txt", 'r') as f:
    f.readline()
    path_cnt = int(f.readline().strip())

use_julia_solver = path_cnt <= HARD_PATH_COUNT_BASELINE

if use_julia_solver:
    from juliacall import Main as jl

    # Julia solver
    jl.include(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'anderson_nnls.jl'))
    jl.seval("using HDF5, LinearAlgebra")

    A_jl, B_jl, b_start_ind = jl.load_nnls_array(f'{PREFIX}/matrix.h5')
    A_fail_jl = jl.load_fail_array(f'{PREFIX}/matrix.h5')

    B = np.asarray(B_jl)

    weight_base_jl = jl.nnls_solve(A_jl, B_jl, THREAD, False)
    weight_base = np.asarray(weight_base_jl)

    final_weights_fullsize = np.zeros(len(weight_base))
    predict_suc_B_base = np.asarray(A_jl * weight_base_jl)

    chrom_acc_sum_dict_base = defaultdict(int)
    chrom_acc_sum_dict_max_base = defaultdict(int)
    for i, (chrom, st) in enumerate(chr_filt_st_list):
        chrom_acc_sum_dict_base[chrom] += predict_suc_B_base[i] - B[i]
        if abs(chrom_acc_sum_dict_base[chrom]) > chrom_acc_sum_dict_max_base[chrom]:
            chrom_acc_sum_dict_max_base[chrom] = abs(chrom_acc_sum_dict_base[chrom])

    # print("="*50)
    # print(f"Base acc sum max")
    # for chrom, acc_sum_max in chrom_acc_sum_dict_max_base.items():
    #     if chrom == 'chrY' and no_chrY:
    #         continue
    #     acc_sum_max_base = chrom_acc_sum_dict_max_base[chrom]
    #     print(f"{chrom} : {acc_sum_max:.4f}")
    # print("="*50)

    nclose_total_weight_dict = defaultdict(float)
    for i, v in enumerate(weight_base):
        for j in path_nclose_set_dict[i]:
            nclose_total_weight_dict[j] += v

    div_nclose_set = set()
    for k, v in nclose_total_weight_dict.items():
        st, ed = k
        # depth check, type 1 only
        if v > 0.2 * meandepth \
        and ppc_contig_data[st][CHR_NAM] != ppc_contig_data[ed][CHR_NAM]:
            # Check if at least one of endpoint have non-diverse depth & not a virtual censat contig
            if (check_near_bnd(ppc_contig_data[st][CHR_NAM], ppc_contig_data[st][CHR_STR], ppc_contig_data[st][CHR_END]) or \
            check_near_bnd(ppc_contig_data[ed][CHR_NAM], ppc_contig_data[ed][CHR_STR], ppc_contig_data[ed][CHR_END])) and \
            not ppc_contig_data[st][CTG_NAM].startswith('virtual_censat_contig'):
                div_nclose_set.add(k)
            # Check if at least one of endpoint is censat contig
            if (ppc_contig_data[st][CTG_CENSAT] != '0' or ppc_contig_data[ed][CTG_CENSAT] != '0') and not ppc_contig_data[st][CTG_NAM].startswith('virtual_censat_contig'):
                div_nclose_set.add(k)


    int2nclose = dict()

    nclose_notusing_idx_dict = defaultdict(list)
    for nclose_pair in div_nclose_set:
        for k, path_nclose_usage in path_nclose_set_dict.items():
            if (nclose_pair not in path_nclose_usage):
                nclose_notusing_idx_dict[nclose_pair].append(k)

    tar_nclose_list = []
    thread_data_jl = jl.Vector[jl.Vector[jl.Int]]()
    for k, v in nclose_notusing_idx_dict.items():
        jl.push_b(thread_data_jl, jl.Vector[jl.Int](v))
        tar_nclose_list.append(k)

    weight_list = jl.run_nnls_map(A_jl, B_jl, thread_data_jl)

    not_essential_nclose = set()

    logging.info(f"Testing nclose pairs count : {len(tar_nclose_list)}")

    for nclose_pair, (weight_jl, predict_suc_B_jl) in zip(tar_nclose_list, weight_list):
        predict_suc_B = np.asarray(predict_suc_B_jl)
        predict_diff = predict_suc_B - predict_suc_B_base

        chrom_acc_sum_dict = defaultdict(int)
        chrom_acc_sum_dict_max = defaultdict(int)
        for i, (chrom, st) in enumerate(chr_filt_st_list):
            chrom_acc_sum_dict[chrom] += predict_suc_B[i] - B[i]
            if abs(chrom_acc_sum_dict[chrom]) > chrom_acc_sum_dict_max[chrom]:
                chrom_acc_sum_dict_max[chrom] = abs(chrom_acc_sum_dict[chrom])

        # print("***********************************")
        # print(nclose_pair)

        max_error_rate = 0
        for chrom, acc_sum_max in chrom_acc_sum_dict_max.items():
            if chrom == 'chrY' and no_chrY:
                continue
            acc_sum_max_base = chrom_acc_sum_dict_max_base[chrom]
            # print(f"{chrom} : {acc_sum_max:.4f} (error rate: {abs(acc_sum_max - acc_sum_max_base) / acc_sum_max_base:.4f})")
            max_error_rate = max(max_error_rate, abs(acc_sum_max - acc_sum_max_base) / acc_sum_max_base)
        
        # print(f"Max base ratio : {max_error_rate:.4f}")
        
        # print("***********************************")
        # print(np.max(np.abs(predict_diff)))
        # for i in range(10):
        #     i = i + 1
        #     print(f"{round(i / 20, 3)}", np.sum(np.abs(predict_diff) > i / 20 * meandepth))

        # if np.sum(np.abs(predict_diff) > WINDOW_DEPTH_ERROR_THRESHOLD * meandepth) < DEPTH_ERROR_ALLOW_LEN / DEPTH_VECTOR_WINDOW:
        #     print("Difference : This nclose pair is not essential.")

        if max_error_rate < BASE_ACCSUMABSMAX_RATIO:
            logging.debug(f"{nclose_pair} : Nclose removed")
            not_essential_nclose.add(nclose_pair)

    logging.info(f'Filtered nclose count by depth : {len(not_essential_nclose)}')
    # logging.info(f'{not_essential_nclose}')

    using_idx_set = set(range(len(weight_base)))
    for nclose_pair in not_essential_nclose:
        using_idx_set &= set(nclose_notusing_idx_dict[nclose_pair])

    final_idx_list = sorted(list(using_idx_set))
    final_idx_array_jl = jl.Vector[jl.Int]([i + 1 for i in final_idx_list]) # 1-index

    A_final_jl = jl.view(A_jl, jl.Colon(), final_idx_array_jl)
    A_fail_final_jl = jl.view(A_fail_jl, jl.Colon(), final_idx_array_jl)

    final_weight_jl = jl.nnls_solve(A_final_jl, B_jl, THREAD, False)
    final_weight = np.asarray(final_weight_jl)

    for i, v in enumerate(final_idx_list):
        final_weights_fullsize[int(v)] = final_weight[i]

    predict_B_succ = np.asarray(A_final_jl * final_weight_jl)
    predict_B_fail = np.asarray(A_fail_final_jl * final_weight_jl)

else:
    import h5py

    from skglm import GeneralizedLinearEstimator
    from skglm.datafits import Quadratic
    from skglm.penalties import PositiveConstraint
    from skglm.solvers import AndersonCD

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
    weights = nnls.coef_
    final_weights_fullsize = weights

    predict_B_succ = A.dot(weights)
    predict_B_fail = A_fail.dot(weights)

b_norm = np.linalg.norm(B)

error = np.linalg.norm(predict_B_succ - B)
predict_B = np.concatenate((predict_B_succ, predict_B_fail))[b_start_ind:]

logging.info(f'Error : {error:.4f}')
logging.info(f'Relative error : {error/b_norm:.4f}')

np.save(f'{PREFIX}/weight.npy', final_weights_fullsize)
np.save(f'{PREFIX}/predict_B.npy', predict_B)
