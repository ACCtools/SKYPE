"""Depth NNLS with a SciPy working set and an approximate Adelie fallback."""

import numpy as np
from scipy.optimize import nnls


COLUMN_BLOCK_SIZE = 256
WORKING_SET_TOLERANCE = 1e-12
WORKING_SET_BATCH_SIZE = 32
FULL_SOLVE_COLUMN_THRESHOLD = 2048
ADELIE_TOLERANCE = 1e-7
ADELIE_MAX_ITERS = 100000


def predict_depth_from_weights(matrix, weights):
    """Accumulate nonzero contributions in Float64 without a full matrix cast."""
    prediction = np.zeros(matrix.shape[0], dtype=np.float64)
    indices = np.flatnonzero(weights)
    for start in range(0, len(indices), COLUMN_BLOCK_SIZE):
        selected = indices[start:start + COLUMN_BLOCK_SIZE]
        block = np.asarray(matrix[:, selected], dtype=np.float64)
        prediction += block @ weights[selected]
    return prediction


def _full_nnls(matrix, target):
    return nnls(matrix, target, maxiter=max(3 * matrix.shape[1], 10000))[0]


def _adelie_nnls(matrix, target):
    """Solve all aggregated columns with Adelie's default accuracy."""
    # Import only when needed so successful SciPy fits avoid Adelie startup.
    from adelie.matrix import dense
    from adelie.solver import bvls

    matrix = np.asfortranarray(matrix, dtype=np.float64)
    p = matrix.shape[1]
    state = bvls(
        dense(matrix, method="naive", n_threads=1),
        np.asarray(target, dtype=np.float64),
        np.zeros(p, dtype=np.float64),
        np.full(p, np.inf, dtype=np.float64),
        tol=ADELIE_TOLERANCE, max_iters=ADELIE_MAX_ITERS, n_threads=1,
    )
    # Adelie can report an error on the state without raising an exception.
    if state.error:
        raise RuntimeError(f"Adelie NNLS fallback failed: {state.error}")
    weights = np.array(state.beta, dtype=np.float64, copy=True)
    if weights.shape != (p,):
        raise RuntimeError(f"Adelie NNLS returned invalid coefficient shape: {weights.shape}")
    if not np.isfinite(weights).all() or np.any(weights < 0):
        raise RuntimeError("Adelie NNLS returned non-finite or negative coefficients")
    return weights


def _solve_working_set(matrix, target, target_norm, *, batch_size=WORKING_SET_BATCH_SIZE,
                       max_rounds=128):
    """Grow SciPy subproblems toward full KKT convergence, falling back to Adelie.

    Admitted columns remain eligible in every subsequent solve, even when
    their coefficient becomes zero. Keeping them prevents repeated eviction
    and re-entry of nearly interchangeable depth profiles. No reference fit
    or previous graph is needed to initialize the working set.
    """
    p = matrix.shape[1]
    norms = np.sqrt(np.einsum("ij,ij->j", matrix, matrix))
    scale = norms * (target_norm if target_norm else 1.0)
    weights = np.zeros(p, dtype=np.float64)
    admitted = np.zeros(p, dtype=bool)
    residual = -target.copy()
    target_energy = float(target @ target)
    max_columns = 0
    rounds = 0
    reason = "working-set round limit"
    for iteration in range(max_rounds + 1):
        gradient = matrix.T @ residual
        violations = np.where(weights > 0, np.abs(gradient), np.maximum(-gradient, 0))
        violations = np.divide(violations, scale, out=np.zeros(p), where=scale > 0)
        if not np.isfinite(violations).all():
            reason = "non-finite working-set optimality residual"
            break
        max_violation = float(np.max(violations, initial=0))
        if max_violation <= WORKING_SET_TOLERANCE:
            return weights, dict(strategy="working_set", working_set_rounds=rounds,
                max_working_columns=max_columns, working_set_fallback=None,
                kkt_scaled_max=max_violation)
        if iteration == max_rounds:
            break
        candidates = np.flatnonzero((~admitted) & (gradient < 0) &
                                    (violations > WORKING_SET_TOLERANCE))
        if not len(candidates):
            reason = "subproblem left an admitted-column optimality violation"
            break
        ranked = np.argsort(-violations[candidates], kind="stable")[:batch_size]
        admitted[candidates[ranked]] = True
        indices = np.flatnonzero(admitted)
        max_columns = max(max_columns, len(indices))
        subproblem = matrix[:, indices]
        try:
            values = _full_nnls(subproblem, target)
        except RuntimeError:
            reason = "subproblem exceeded its iteration limit"
            break
        if not np.isfinite(values).all() or np.any(values < 0):
            reason = "subproblem returned invalid coefficients"
            break
        new_residual = subproblem @ values - target
        old_sse = float(residual @ residual)
        new_sse = float(new_residual @ new_residual)
        allowance = 64 * np.finfo(float).eps * max(target_energy, old_sse, 1.)
        if not np.isfinite(new_sse) or new_sse > old_sse + allowance:
            reason = "subproblem increased the objective"
            break
        weights[:] = 0
        weights[indices] = values
        residual = new_residual
        rounds += 1
    return _adelie_nnls(matrix, target), dict(strategy="adelie_fallback",
        working_set_rounds=rounds, max_working_columns=max_columns,
        working_set_fallback=reason)


def aggregate_identical_depth_rows(matrix, target):
    """Combine identical consecutive rows without changing the minimizer.

    For a run of k identical feature rows a, replace its observations with
    sqrt(k)*a and sum(b)/sqrt(k). The removed within-run SSE is independent
    of the coefficients. Runs need not respect chromosome boundaries: only
    exact equality of every model feature permits aggregation.
    """
    changed = np.zeros(matrix.shape[0] - 1, dtype=bool)
    for start in range(0, matrix.shape[1], COLUMN_BLOCK_SIZE):
        block = matrix[:, start:start + COLUMN_BLOCK_SIZE]
        if not np.isfinite(block).all():
            raise ValueError("NNLS matrix contains non-finite values")
        changed |= np.any(block[1:] != block[:-1], axis=1)
    starts = np.r_[0, np.flatnonzero(changed) + 1]
    counts = np.diff(np.r_[starts, len(target)])
    root_counts = np.sqrt(counts)
    reduced = np.empty((len(starts), matrix.shape[1]), dtype=np.float64, order="F")
    for start in range(0, matrix.shape[1], COLUMN_BLOCK_SIZE):
        stop = start + COLUMN_BLOCK_SIZE
        reduced[:, start:stop] = matrix[starts, start:stop] * root_counts[:, None]
    reduced_target = np.add.reduceat(target, starts) / root_counts
    return reduced, reduced_target


def _fit_diagnostics(matrix, target, weights):
    """Compute residual statistics without a gradient over all columns."""
    if not np.isfinite(weights).all() or np.any(weights < 0):
        raise RuntimeError("NNLS solver returned non-finite or negative coefficients")
    prediction = predict_depth_from_weights(matrix, weights)
    residual = prediction - target
    if not np.isfinite(residual).all():
        raise RuntimeError("NNLS residual contains non-finite values")
    target_norm = np.linalg.norm(target)
    return residual, {
        "sse": float(residual @ residual),
        "relative_error": float(np.linalg.norm(residual) / target_norm) if target_norm else None,
        "min_weight": float(np.min(weights)),
        "positive_weights": int(np.count_nonzero(weights > 0)),
    }


def certify_depth_nnls(matrix, target, weights, *, target_norm=None):
    """Recompute scaled KKT violations on the supplied matrix.

    The certificate is invariant to positive rescaling of target or columns:
    max_j violation_j / (||A_j|| * ||b||). At positive coefficients violation
    is |A_j.T @ (Aw-b)|; at zero it is max(-A_j.T @ (Aw-b), 0).
    Aggregated problems pass the original target norm to preserve this scale.
    Residual statistics always describe the supplied matrix and target.
    """
    residual, diagnostics = _fit_diagnostics(matrix, target, weights)
    if target_norm is None:
        target_norm = np.linalg.norm(target)
    max_violation = 0.0
    for start in range(0, matrix.shape[1], COLUMN_BLOCK_SIZE):
        stop = start + COLUMN_BLOCK_SIZE
        block = np.asarray(matrix[:, start:stop], dtype=np.float64)
        gradient = block.T @ residual
        norms = np.linalg.norm(block, axis=0)
        violation = np.where(weights[start:stop] > 0, np.abs(gradient), np.maximum(-gradient, 0))
        denominator = norms * (target_norm if target_norm else 1.0)
        normalized = np.divide(violation, denominator, out=np.zeros_like(violation), where=denominator > 0)
        if not np.isfinite(normalized).all():
            raise RuntimeError("NNLS optimality check contains non-finite values")
        max_violation = max(max_violation, float(np.max(normalized, initial=0)))
    diagnostics["kkt_scaled_max"] = max_violation
    return diagnostics


def solve_depth_nnls(matrix, target, *, return_diagnostics=False):
    """Fit all columns without a separate KKT check after the solver returns."""
    matrix = np.asarray(matrix)
    target = np.asarray(target, dtype=np.float64)
    if matrix.ndim != 2:
        raise ValueError(f"NNLS matrix must be 2-dimensional: {matrix.shape}")
    if target.shape != (matrix.shape[0],):
        raise ValueError(f"NNLS target has shape {target.shape}, expected {(matrix.shape[0],)}")
    if not matrix.shape[0] or not matrix.shape[1]:
        raise ValueError("NNLS matrix must have observation rows and feature columns")
    if not np.isfinite(target).all():
        raise ValueError("NNLS target contains non-finite values")
    reduced, reduced_target = aggregate_identical_depth_rows(matrix, target)
    target_norm = np.linalg.norm(target)
    if matrix.shape[1] <= FULL_SOLVE_COLUMN_THRESHOLD:
        weights = _full_nnls(reduced, reduced_target)
        strategy = dict(strategy="full", working_set_rounds=0,
                        max_working_columns=matrix.shape[1], working_set_fallback=None)
    else:
        weights, strategy = _solve_working_set(reduced, reduced_target, target_norm)
    uses_adelie = strategy["strategy"] == "adelie_fallback"
    uses_working_set = strategy["strategy"] == "working_set"
    if not uses_working_set:
        strategy["kkt_scaled_max"] = None
    _, diagnostics = _fit_diagnostics(matrix, target, weights)
    diagnostics.update(
        solver="adelie.solver.bvls" if uses_adelie else "scipy.optimize.nnls",
        dtype="float64", kkt_tolerance=WORKING_SET_TOLERANCE if uses_working_set else None,
        kkt_matrix="aggregated" if uses_working_set else None,
        convergence_check="aggregated_kkt" if uses_working_set else "solver_default",
        solver_tolerance=ADELIE_TOLERANCE if uses_adelie else None,
        solver_max_iters=ADELIE_MAX_ITERS if uses_adelie else None,
        observation_rows=matrix.shape[0], aggregated_rows=reduced.shape[0],
        feature_columns=matrix.shape[1], converged=True if uses_working_set else None,
        **strategy,
    )
    return (weights, diagnostics) if return_diagnostics else weights
