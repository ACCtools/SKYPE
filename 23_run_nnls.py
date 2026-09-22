"""Fit stage-22 depth features with raw NNLS or the optional robust policy."""

import argparse
import json
import logging
import os
import pickle
import sys

import h5py
import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from denoised_relative_error import (  # noqa: E402
    TV_LAMBDA_OVER_NOISE_SIGMA,
    calculate_denoised_relative_error,
)
from depth_nnls import predict_depth_from_weights, solve_depth_nnls  # noqa: E402
from nclose_tracking import (  # noqa: E402
    load_filter_status,
    load_path_usage,
    record_filter_stage,
    save_filter_status,
)
from skype_utils import DEPTH_ONLY_MATRIX_CONTRACT, LOG_LEVEL  # noqa: E402,F401


MATRIX_CONTRACT = DEPTH_ONLY_MATRIX_CONTRACT


def fit_raw_nnls(matrix, target, *, return_diagnostics=False):
    """Run the Float64 SciPy working-set fit with an approximate Adelie fallback."""
    return solve_depth_nnls(matrix, target, return_diagnostics=return_diagnostics)


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


def _positive_thread_count(value):
    count = int(value)
    if count < 1:
        raise argparse.ArgumentTypeError("threads must be a positive integer")
    return count


def build_parser():
    parser = argparse.ArgumentParser(description="SKYPE raw depth NNLS")
    parser.add_argument("prefix")
    parser.add_argument(
        "-t", "--thread", "--threads", type=_positive_thread_count, default=1,
        help="BLAS thread limit for fitting and prediction (default: 1)",
    )
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

    # One solve per pass; stage 24 may request one additional graph/fit pass.
    solver_matrix = feature_depth.T
    logging.info("NNLS BLAS thread limit : %d", args.thread)
    with threadpool_limits(limits=args.thread, user_api="blas"):
        high_mask = np.zeros(len(target_depth), dtype=bool)
        if matrix_meta.get('depth_policy') == 'nclose_huber':
            high_mask[matrix_meta.get('high_depth_rows', [])] = True
        if high_mask.any():
            from high_depth import solve_robust
            initial_path = os.path.join(prefix, 'weight.npy')
            initial = np.load(initial_path) if os.path.isfile(initial_path) else None
            if initial is not None and initial.shape != (solver_matrix.shape[1],):
                initial = None
            weights, solver_diagnostics = solve_robust(
                solver_matrix, np.asarray(target_depth, dtype=float), high_mask,
                np.asarray(matrix_meta['high_depth_tau']), initial=initial)
        else:
            weights, solver_diagnostics = fit_raw_nnls(
                solver_matrix, target_depth, return_diagnostics=True
            )
            if matrix_meta.get('depth_policy') == 'nclose_huber':
                solver_diagnostics['high_depth_fallback'] = 'no_high_depth_rows'
        predict_depth = predict_depth_from_weights(solver_matrix, weights)
        predict_fail = predict_depth_from_weights(feature_fail.T, weights)
        solver_diagnostics["blas_threads_requested"] = args.thread
        solver_diagnostics["blas_threadpools"] = [
            {key: pool[key] for key in ("filepath", "internal_api", "num_threads")}
            for pool in threadpool_info() if pool["user_api"] == "blas"
        ]
    predict_all = np.concatenate((predict_depth, predict_fail))
    if matrix_meta.get('depth_policy') == 'nclose_huber':
        np.save(os.path.join(prefix, 'unexplained_high_depth.npy'),
                np.maximum(target_depth[high_mask] - predict_depth[high_mask]
                           - np.asarray(matrix_meta.get('high_depth_tau', [])), 0))

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
    logging.info("NNLS solver : %s", solver_diagnostics["solver"])
    if solver_diagnostics["kkt_scaled_max"] is None:
        logging.info("NNLS convergence : solver-default; strict KKT not checked")
        if solver_diagnostics["solver_tolerance"] is not None:
            logging.info("NNLS solver tolerance : %g", solver_diagnostics["solver_tolerance"])
    else:
        logging.info("NNLS aggregated-matrix scaled KKT violation : %.3g", solver_diagnostics["kkt_scaled_max"])
    logging.info(
        "NNLS strategy : %s; working rounds : %d; largest subproblem : %d columns",
        solver_diagnostics["strategy"], solver_diagnostics["working_set_rounds"],
        solver_diagnostics["max_working_columns"],
    )
    if solver_diagnostics["working_set_fallback"] is not None:
        logging.info("NNLS working-set fallback : %s", solver_diagnostics["working_set_fallback"])
    logging.info(
        "NNLS equivalent depth rows : %d -> %d",
        solver_diagnostics["observation_rows"], solver_diagnostics["aggregated_rows"],
    )
    logging.info("Error : %.4f", error)
    logging.info("Relative error : %.4f", relative_error)
    logging.info(
        "Denoised target TV lambda/noise sigma : %g",
        TV_LAMBDA_OVER_NOISE_SIGMA,
    )
    logging.info("Denoised error : %.4f", denoised_error)
    logging.info("Denoised relative error : %.4f", denoised_relative_error)

    np.save(os.path.join(prefix, "weight.npy"), weights)
    with open(os.path.join(prefix, "nnls_diagnostics.json"), "w") as handle:
        json.dump(solver_diagnostics, handle, indent=2)
        handle.write("\n")
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
