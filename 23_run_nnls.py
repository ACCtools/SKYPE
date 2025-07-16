import h5py
import logging
import warnings
import argparse

import numpy as np
import pickle as pkl

from scipy.sparse import issparse

from skglm import GeneralizedLinearEstimator
from skglm.datafits import Quadratic, QuadraticSVC
from skglm.penalties import PositiveConstraint
from skglm.solvers import AndersonCD
from skglm.utils.jit_compilation import compiled_clone

from sklearn.utils.validation import check_array

from threadpoolctl import threadpool_limits

PATH_WEIGHT_INIT_RATIO = 0.05

warnings.simplefilter(action='ignore', category=FutureWarning)

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("23_run_nnls start")

parser = argparse.ArgumentParser(description="SKYPE depth analysis")
parser.add_argument("prefix", help="Prefix for pipeline")
parser.add_argument("-t", "--thread", help="Number of threads", type=int)
args = parser.parse_args()

PREFIX = args.prefix
THREAD = args.thread

with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
    dA = f["A"]
    mx, my = dA.shape
    A = np.empty((mx + 1, my), dtype=dA.dtype)
    dA.read_direct(A, dest_sel=np.s_[:-1, :])

    dB = f["B"]

    B = np.zeros(dB.shape[0] + 1, dtype=dB.dtype)
    dB.read_direct(B, dest_sel=np.s_[:-1])

    dAf = f["A_fail"]
    A_fail = np.empty(dAf.shape, dtype=dAf.dtype)
    dAf.read_direct(A_fail)

    b_start_ind = int(f["B_depth_start"][()])


with open(f"{PREFIX}/23_input.pkl", "rb") as f:
    dep_list, init_cols, w_pri = pkl.load(f)
dep_np = np.array(dep_list, dtype=np.float32)

b_rms = np.sqrt(np.mean(B[:-1] ** 2))
PATH_WEIGHT = b_rms * PATH_WEIGHT_INIT_RATIO
A[-1, :] = dep_np * PATH_WEIGHT

nnls = GeneralizedLinearEstimator(
    datafit=Quadratic(),
    penalty=PositiveConstraint(),
    solver=AndersonCD(fit_intercept=False)
)

with threadpool_limits(limits=THREAD):
    nnls.fit(A, B)

weights = nnls.coef_.copy()

predict_suc_B = A[:-1, :].dot(weights)
b_norm = np.linalg.norm(B[:-1])
error = np.linalg.norm(predict_suc_B - B[:-1])

predict_B = np.concatenate([predict_suc_B, A_fail.dot(weights)])[b_start_ind:]
logging.info(f'Error : {error:.4f}')
logging.info(f'Relative error : {error/b_norm:.4f}')
np.save(f'{PREFIX}/weight.npy', weights)
np.save(f'{PREFIX}/predict_B.npy', predict_B)
