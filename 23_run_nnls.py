import h5py
import logging
import warnings
import argparse

import numpy as np
import pickle as pkl

from skglm import GeneralizedLinearEstimator
from skglm.datafits import Quadratic
from skglm.penalties import PositiveConstraint
from skglm.solvers import AndersonCD
from threadpoolctl import threadpool_limits

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
    A = np.empty(dA.shape, dtype=dA.dtype)
    dA.read_direct(A)

    dB = f["B"]
    B = np.empty(dB.shape, dtype=dB.dtype)
    dB.read_direct(B)

    dAf = f["A_fail"]
    A_fail = np.empty(dAf.shape, dtype=dAf.dtype)
    dAf.read_direct(A_fail)

    b_start_ind = int(f["B_depth_start"][()])


with open(f"{PREFIX}/23_input.pkl", "rb") as f:
    dep_list, init_cols, w_pri = pkl.load(f)

n = len(dep_list)
x0 = np.zeros(n, dtype=np.float32)
x0[init_cols] = w_pri

max_order = dep_list[0]
for order in range(0, max_order + 1):
    tar_ind = None
    for i, v in enumerate(dep_list):
        if v <= order:
            tar_ind = i
            break

    assert(tar_ind is not None)

    logging.info(f"Now order : {order}")
    nnls = GeneralizedLinearEstimator(
        datafit=Quadratic(),
        penalty=PositiveConstraint(),
        solver=AndersonCD(fit_intercept=False, warm_start=True)
    )

    nnls.coef_ = x0[tar_ind:]
    with threadpool_limits(limits=THREAD):
        nnls.fit(A[:, tar_ind:], B)

    if max_order == order:
        weights = nnls.coef_
    else:
        x0[tar_ind:] = nnls.coef_

predict_suc_B = A.dot(weights)
error = np.linalg.norm(predict_suc_B - B)
b_norm = np.linalg.norm(B)

predict_B = np.concatenate([predict_suc_B, A_fail.dot(weights)])[b_start_ind:]
logging.info(f'Error : {error:.4f}')
logging.info(f'Relative error : {error/b_norm:.4f}')
np.save(f'{PREFIX}/weight.npy', weights)
np.save(f'{PREFIX}/predict_B.npy', predict_B)
