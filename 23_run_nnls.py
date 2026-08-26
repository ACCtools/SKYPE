"""Fit every stage-22 depth feature with one raw non-negative solve."""

import argparse
import logging
import os
import pickle
import sys

import h5py
import numpy as np
from adelie.solver import bvls

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from denoised_relative_error import (  # noqa: E402
    TV_LAMBDA_OVER_NOISE_SIGMA,
    calculate_denoised_relative_error,
)
from nclose_tracking import (  # noqa: E402
    load_filter_status,
    load_path_usage,
    record_filter_stage,
    save_filter_status,
)
from skype_utils import DEPTH_ONLY_MATRIX_CONTRACT, LOG_LEVEL  # noqa: E402,F401


MATRIX_CONTRACT = DEPTH_ONLY_MATRIX_CONTRACT


def fit_raw_nnls(matrix, target):
    """Run the pipeline's single non-negative least-squares solve."""

    matrix = np.asarray(matrix)
    target = np.asarray(target, dtype=matrix.dtype)
    if matrix.ndim != 2:
        raise ValueError(f"NNLS matrix must be 2-dimensional: {matrix.shape}")
    if target.shape != (matrix.shape[0],):
        raise ValueError(
            f"NNLS target has shape {target.shape}, expected {(matrix.shape[0],)}"
        )
    if matrix.shape[1] == 0:
        raise ValueError("NNLS matrix has no feature columns")

    lower = np.zeros(matrix.shape[1], dtype=matrix.dtype)
    upper = np.full(matrix.shape[1], np.finfo(matrix.dtype).max, dtype=matrix.dtype)
    state = bvls(
        matrix,
        target,
        lower,
        upper,
        n_threads=1,
    )
    return np.asarray(state.beta)


def _contract_text(value):
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def load_depth_only_matrix(prefix):
    """Load and validate the stage-22 feature-major depth matrix."""

    matrix_path = os.path.join(prefix, "matrix.h5")
    with h5py.File(matrix_path, "r") as handle:
        contract = _contract_text(handle.attrs.get("matrix_contract", ""))
        if contract != MATRIX_CONTRACT:
            raise ValueError(
                "Stage-23 requires a depth-only stage-22 matrix; "
                "rerun the pipeline from stage 22."
            )

        required = {"A", "A_fail", "B", "B_fail"}
        missing = sorted(required - set(handle.keys()))
        if missing:
            raise ValueError(f"Stage-22 matrix is missing datasets: {missing}")

        feature_depth = handle["A"][:]
        feature_fail = handle["A_fail"][:]
        target_depth = handle["B"][:]
        target_fail = handle["B_fail"][:]

    if feature_depth.ndim != 2 or feature_fail.ndim != 2:
        raise ValueError("Stage-22 A/A_fail datasets must be 2-dimensional")
    if feature_depth.shape[0] != feature_fail.shape[0]:
        raise ValueError(
            "Stage-22 A/A_fail feature count mismatch: "
            f"{feature_depth.shape[0]} != {feature_fail.shape[0]}"
        )
    if feature_depth.shape[1] != len(target_depth):
        raise ValueError(
            "Stage-22 clean-depth row mismatch: "
            f"{feature_depth.shape[1]} != {len(target_depth)}"
        )
    if feature_fail.shape[1] != len(target_fail):
        raise ValueError(
            "Stage-22 failed-depth row mismatch: "
            f"{feature_fail.shape[1]} != {len(target_fail)}"
        )
    return feature_depth, feature_fail, target_depth, target_fail


def build_parser():
    parser = argparse.ArgumentParser(description="SKYPE raw depth NNLS")
    parser.add_argument("prefix")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    prefix = args.prefix
    logging.info("23_run_nnls start")

    with open(os.path.join(prefix, "23_input.pkl"), "rb") as handle:
        matrix_meta = pickle.load(handle)
    if matrix_meta.get("matrix_contract") != MATRIX_CONTRACT:
        raise ValueError(
            "Stage-23 requires depth-only stage-22 metadata; "
            "rerun the pipeline from stage 22."
        )

    chr_filt_st_list = matrix_meta.get("chr_filt_st_list", [])
    feature_depth, feature_fail, target_depth, _ = load_depth_only_matrix(prefix)
    if (
        int(matrix_meta.get("B_depth_start", -1)) != 0
        or int(matrix_meta.get("B_depth_end", -1)) != len(target_depth)
    ):
        raise ValueError("Stage-22 depth metadata does not match the NNLS target")
    if len(chr_filt_st_list) != len(target_depth):
        raise ValueError(
            "Stage-22 clean-depth coordinates do not match the NNLS target: "
            f"{len(chr_filt_st_list)} != {len(target_depth)}"
        )

    # This is the only optimizer invocation in the native SKYPE pipeline.
    solver_matrix = feature_depth.T
    weights = fit_raw_nnls(solver_matrix, target_depth)
    predict_depth = solver_matrix.dot(weights)
    predict_fail = feature_fail.T.dot(weights)
    predict_all = np.concatenate((predict_depth, predict_fail))

    target_norm = np.linalg.norm(target_depth)
    error = np.linalg.norm(predict_depth - target_depth)
    relative_error = error / target_norm if target_norm else float("nan")
    denoised_error, denoised_relative_error, _ = calculate_denoised_relative_error(
        chr_filt_st_list,
        target_depth,
        predict_depth,
        TV_LAMBDA_OVER_NOISE_SIGMA,
    )

    logging.info("Raw NNLS feature count : %d", len(weights))
    logging.info("Error : %.4f", error)
    logging.info("Relative error : %.4f", relative_error)
    logging.info(
        "Denoised target TV lambda/noise sigma : %g",
        TV_LAMBDA_OVER_NOISE_SIGMA,
    )
    logging.info("Denoised error : %.4f", denoised_error)
    logging.info("Denoised relative error : %.4f", denoised_relative_error)

    np.save(os.path.join(prefix, "weight.npy"), weights)
    np.save(os.path.join(prefix, "predict_B.npy"), predict_all)
    active_columns = list(range(len(weights)))
    with open(os.path.join(prefix, "A_idx_list.pkl"), "wb") as handle:
        pickle.dump(active_columns, handle)

    path_usage = load_path_usage(prefix, expected_len=len(weights))
    status = load_filter_status(prefix)
    status["stages"].pop("base", None)
    status["stages"].pop("filter", None)
    status["stages"].pop("cluster", None)
    provenance_active_columns = status["stages"]["initial"]["active_columns"]
    record_filter_stage(
        status,
        stage="base",
        previous_stage="initial",
        path_usage=path_usage,
        active_columns=provenance_active_columns,
        direct_reasons={},
        cofiltered_reason="FILTERED_02_COFILTERED_PATH",
    )
    save_filter_status(prefix, status)


if __name__ == "__main__":
    main()
