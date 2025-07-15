import h5py
import logging
import warnings
import argparse

import numpy as np
import pickle as pkl

from scipy.sparse import issparse

from skglm import GeneralizedLinearEstimator
from skglm.datafits import Quadratic, Logistic, QuadraticSVC
from skglm.penalties import PositiveConstraint, WeightedL1
from skglm.solvers import AndersonCD
from skglm.utils.jit_compilation import compiled_clone

from sklearn.utils.validation import check_array

from threadpoolctl import threadpool_limits

def _glm_fit_custom(X, y, model):
    solver, datafit, penalty = model.solver, model.datafit, model.penalty

    is_classif = False
    fit_intercept = solver.fit_intercept

    check_X_params = dict(
        dtype=[np.float64, np.float32], order='F',
        accept_sparse='csc', copy=fit_intercept)
    check_y_params = dict(ensure_2d=False, order='F')

    X, y = model._validate_data(
        X, y, validate_separately=(check_X_params, check_y_params))
    X = check_array(X, 'csc', dtype=[np.float64, np.float32],
                    order='F', copy=False, accept_large_sparse=False)
    y = check_array(y, 'csc', dtype=X.dtype.type, order='F', copy=False,
                    ensure_2d=False)

    if y.ndim == 2 and y.shape[1] == 1:
        warnings.warn("DataConversionWarning('A column-vector y"
                      " was passed when a 1d array was expected")
        y = y[:, 0]

    if not hasattr(model, "n_features_in_"):
        model.n_features_in_ = X.shape[1]

    n_samples = X.shape[0]
    if n_samples != y.shape[0]:
        raise ValueError("X and y have inconsistent dimensions (%d != %d)"
                         % (n_samples, y.shape[0]))

    # if not model.warm_start or not hasattr(model, "coef_"):
    if not solver.warm_start or not hasattr(model, "coef_"):
        model.coef_ = None

    X_ = X
    n_samples, n_features = X_.shape

    penalty_jit = compiled_clone(penalty)
    datafit_jit = compiled_clone(datafit, to_float32=X.dtype == np.float32)
    if issparse(X):
        datafit_jit.initialize_sparse(X_.data, X_.indptr, X_.indices, y)
    else:
        datafit_jit.initialize(X_, y)

    # if model.warm_start and hasattr(model, 'coef_') and model.coef_ is not None:
    if solver.warm_start and hasattr(model, 'coef_') and model.coef_ is not None:
        if isinstance(datafit, QuadraticSVC):
            w = model.dual_coef_[0, :].copy()
        elif is_classif:
            w = model.coef_[0, :].copy()
        else:
            w = model.coef_.copy()
        if fit_intercept:
            w = np.hstack([w, model.intercept_])
        Xw = X_ @ w[:w.shape[0] - fit_intercept] + fit_intercept * w[-1]
    else:
        # TODO this should be solver.get_init() do delegate the work
        if y.ndim == 1:
            w = np.zeros(n_features + fit_intercept, dtype=X_.dtype)
            Xw = np.zeros(n_samples, dtype=X_.dtype)
        else:  # multitask
            w = np.zeros((n_features + fit_intercept, y.shape[1]), dtype=X_.dtype)
            Xw = np.zeros(y.shape, dtype=X_.dtype)

    coefs, p_obj, kkt = solver.solve(X_, y, datafit_jit, penalty_jit, w, Xw, run_checks=False)
    model.coef_, model.stop_crit_ = coefs[:n_features], kkt
    if y.ndim == 1:
        model.intercept_ = coefs[-1] if fit_intercept else 0.
    else:
        model.intercept_ = coefs[-1, :] if fit_intercept else np.zeros(
            y.shape[1])

    model.n_iter_ = len(p_obj)

    return model

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

    nnls.coef_ = x0[tar_ind:].copy()
    with threadpool_limits(limits=THREAD):
        _glm_fit_custom(A[:, tar_ind:], B, nnls)

    if max_order == order:
        weights = nnls.coef_.copy()
    else:
        x0[tar_ind:] = nnls.coef_.copy()

predict_suc_B = A.dot(weights)
error = np.linalg.norm(predict_suc_B - B)
b_norm = np.linalg.norm(B)

predict_B = np.concatenate([predict_suc_B, A_fail.dot(weights)])[b_start_ind:]
logging.info(f'Error : {error:.4f}')
logging.info(f'Relative error : {error/b_norm:.4f}')
np.save(f'{PREFIX}/weight.npy', weights)
np.save(f'{PREFIX}/predict_B.npy', predict_B)
