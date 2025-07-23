import h5py
import logging
import warnings
import argparse
import pickle as pkl

from collections import defaultdict

from networkx import add_path, edge_betweenness_centrality_subset
import numpy as np
import pandas as pd

from skglm import GeneralizedLinearEstimator
from skglm.datafits import Quadratic
from skglm.penalties import PositiveConstraint
from skglm.solvers import AndersonCD
from threadpoolctl import threadpool_limits


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
DEPTH_ERROR_ALLOW_LEN = 1*M

ACC_SUM_ERROR_THRESHOLD = 0.5

VCF_FLANKING_LENGTH = 1*M
NCLOSE_SIM_COMPARE_RAITO = 1.2
NCLOSE_SIM_DIFF_THRESHOLD = 5

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

ppc_contig_data = import_ppc_contig_data(PREPROCESSED_PAF_FILE_PATH)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t',
                 names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

meandepth = df['meandepth'].median()

path_count = 0

with open(f'{PREFIX}/23_input.pkl', 'rb') as f:
#     dep_list, init_cols, w_pri, using_ncnt_array, \
#     div_nclose_set, each_nclose_notusing_path_dict, \
    path_nclose_set_dict, amplitude = pkl.load(f)

with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
    dA = f["A"]
    A = np.empty(dA.shape, dtype=dA.dtype)
    At = np.empty(dA.shape, dtype=dA.dtype)
    path_count = dA.shape[1]
    dA.read_direct(A)

    dB = f["B"]
    B = np.empty(dB.shape, dtype=dB.dtype)
    dB.read_direct(B)

    dAf = f["A_fail"]
    A_fail = np.empty(dAf.shape, dtype=dAf.dtype)
    A_failt = np.empty(dAf.shape, dtype=dAf.dtype)
    dAf.read_direct(A_fail)

    b_start_ind = int(f["B_depth_start"][()])

nnls = GeneralizedLinearEstimator(
    datafit=Quadratic(),
    penalty=PositiveConstraint(),
    solver=AndersonCD(fit_intercept=False)
)

with threadpool_limits(limits=THREAD):
    nnls.fit(A, B)

weights_base = nnls.coef_.copy()

final_weights_base = np.zeros(len(weights_base))

predict_suc_B_base = A.dot(weights_base)

nclose_total_weight_dict = defaultdict(float)

for i, v in enumerate(weights_base):
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

using_ncnt_array = []

for i, path_nclose_usage in path_nclose_set_dict.items():
    using_path = True
    for nclose_pair in path_nclose_usage:
        if nclose_pair in div_nclose_set:
            using_path = False
            break
    if using_path:
        using_ncnt_array.append(i)

using_ncnt_set = set(using_ncnt_array)

each_nclose_notusing_path_dict = defaultdict(list)
for nclose_pair in div_nclose_set:
    for i, path_nclose_usage in path_nclose_set_dict.items():
        if (nclose_pair not in path_nclose_usage) and (i not in using_ncnt_set):
            each_nclose_notusing_path_dict[nclose_pair].append(i)

using_ncnt_array = np.array(sorted(using_ncnt_array))

print(f"Div nclose set length : {len(div_nclose_set)}")
print(amplitude)

cn = len(using_ncnt_array)
At[:, :cn] = A[:, using_ncnt_array]
Atl = At[:, :cn]    
A_failt[:, :cn] = A_fail[:, using_ncnt_array]
A_failtl = A_failt[:, :cn]

print(Atl.shape, A_failtl.shape)

not_essential_nclose = set()

for nclose in list(div_nclose_set)[::-1]:

    print(f"Current nclose : {nclose}")
    var_ncnt_nclose_array = each_nclose_notusing_path_dict[nclose]
    vn = len(var_ncnt_nclose_array)

    assert(len(set(var_ncnt_nclose_array) & set(using_ncnt_array)) == 0)
    At[:, cn:cn+vn] = A[:, var_ncnt_nclose_array]
    Atl = At[:, :cn+vn]

    with threadpool_limits(limits=THREAD):
        nnls.fit(Atl, B)

    weights = nnls.coef_.copy()

    predict_suc_B = Atl.dot(weights)

    predict_diff = predict_suc_B - predict_suc_B_base

    print(np.max(np.abs(predict_diff)))
    for i in [1, 2, 3, 4, 5, 6, 7, 8, 9]:
        print(f"0.{i}", np.sum(np.abs(predict_diff) > i / 10 * amplitude))

    if np.sum(np.abs(predict_diff) > ACC_SUM_ERROR_THRESHOLD * meandepth) < DEPTH_ERROR_ALLOW_LEN / DEPTH_VECTOR_WINDOW:
        not_essential_nclose.add(nclose)


logging.info(f'Filtered nclose count : {len(not_essential_nclose)}')
logging.info(f'{not_essential_nclose}')

At[:, :cn] = A[:, using_ncnt_array]
Atl = At[:, :cn]    

add_path_set = set()
for k, s in each_nclose_notusing_path_dict.items():
    add_path_set |= set(s)
# For every non-essential nclose, we will add the paths that are not using the nclose
# -> Add intersection of all paths that not using non-essential nclose
for i in not_essential_nclose:
    add_path_set &= set(each_nclose_notusing_path_dict[i])

add_path_array = np.array(sorted(list(add_path_set)))
avn = len(add_path_array)


if avn > 0:
    At[:, cn:cn+avn] = A[:, add_path_array]
    A_failt[:, cn:cn+avn] = A_fail[:, add_path_array]

    Atl = At[:, :cn+avn]
    A_failtl = A_failt[:, :cn+avn]


with threadpool_limits(limits=THREAD):
    nnls.fit(Atl, B)

final_weights = nnls.coef_.copy()

final_index = np.concatenate((using_ncnt_array, add_path_array))
for i, v in enumerate(final_index):
    final_weights_base[int(v)] = final_weights[i]

predict_B_succ = Atl.dot(final_weights)
predict_B_fail = A_failtl.dot(final_weights)
b_norm = np.linalg.norm(B)
error = np.linalg.norm(predict_B_succ - B)
predict_B = np.concatenate((predict_B_succ, predict_B_fail))[b_start_ind:]

logging.info(f'Error : {error:.4f}')
logging.info(f'Relative error : {error/b_norm:.4f}')

np.save(f'{PREFIX}/weight.npy',  final_weights_base)
np.save(f'{PREFIX}/predict_B.npy', predict_B)
