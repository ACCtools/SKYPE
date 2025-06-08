from concurrent.futures import thread
import h5py
import logging
import argparse

import numpy as np

from skglm import GeneralizedLinearEstimator
from skglm.datafits import Quadratic
from skglm.penalties import PositiveConstraint
from skglm.solvers import AndersonCD
from threadpoolctl import threadpool_limits

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


nnls = GeneralizedLinearEstimator(
    datafit=Quadratic(),
    penalty=PositiveConstraint(),
    solver=AndersonCD(fit_intercept=False)
)

with threadpool_limits(limits=THREAD):
    nnls.fit(A, B)

weights = nnls.coef_

predict_suc_B = A.dot(weights)
error   = np.linalg.norm(predict_suc_B - B)
b_norm  = np.linalg.norm(B)

predict_B = np.concatenate([predict_suc_B, A_fail.dot(weights)])[b_start_ind:]

logging.info(f'Error       : {error:.4f}')
logging.info(f'Relative err: {error/b_norm:.4f}')

np.save(f'{PREFIX}/weight.npy', weights)
np.save(f'{PREFIX}/predict_B.npy', predict_B)
