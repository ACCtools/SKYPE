"""Optional NClose-gated high-depth rows and nonnegative robust fitting."""
from pathlib import Path
import pickle
import logging

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from native_vcf_bnd import source_junctions
from depth_nnls import (
    aggregate_identical_depth_rows, predict_depth_from_weights, solve_depth_nnls,
    _full_nnls,
)
from denoised_relative_error import estimate_noise_sigma_by_chromosome


HIGH_DEPTH_COLUMNS = {
    'chrom': 'object', 'start': 'int64', 'end': 'int64', 'depth': 'float64',
    'nclose_gate': 'bool', 'nclose_ids': 'object', 'nclose_overlap_bp': 'int64',
    'noise_sigma': 'float64', 'direct_nclose_gate': 'bool',
    'block_start': 'int64', 'block_end': 'int64', 'block_nclose_ids': 'object',
}


def empty_high_depth_table():
    """Keep a usable schema and boolean masks even when no candidate exists."""
    return pd.DataFrame({name: pd.Series(dtype=dtype)
                         for name, dtype in HIGH_DEPTH_COLUMNS.items()})


def nclose_intervals(prefix):
    """Use candidate membership only, never a fitted weight or raw read count.

    Stage 22 saves the current topology before selecting high-depth rows.
    Source geometry and membership do not depend on fitted coefficients.
    """
    root = Path(prefix)
    with (root / 'structure_nclose_model.pkl').open('rb') as f:
        model = pickle.load(f)
    with (root / '01_nclose_data.pkl').open('rb') as f:
        nodes = pickle.load(f)['contig_data']
    represented = {key for s in model['structures'] if s.get('feature_index') is not None
                   for key, count in s.get('source_nclose_counts', {}).items() if count > 0}
    intervals = {}
    for key, event in model['nclose_sources'].items():
        if event['kind'] != 'bnd' or key not in represented:
            continue
        node_ids = {i for j in source_junctions(event, nodes) for i in j['node_pair']}
        for i in node_ids:
            node = nodes[i]
            if int(node[8]) > int(node[7]):
                nid = model['ncloses'][model['nclose_aliases'][key]]['nclose_id']
                intervals.setdefault(node[5], []).append((int(node[7]), int(node[8]), nid))
    return intervals


def high_depth_gate(prefix, frame, old_clean, censat_intervals, *, use_nclose=True):
    threshold = 3 * float(frame.meandepth.median())
    candidates = [row for row in frame.itertuples(index=False)
                  if row.meandepth > threshold
                  and not any(a <= row.st and row.nd <= b
                              for a, b in censat_intervals.get(row.chr, []))]
    if not candidates:
        return empty_high_depth_table()
    intervals = nclose_intervals(prefix) if use_nclose else {}
    depth_by_coord = {(row.chr, int(row.st)): float(row.meandepth) for row in frame.itertuples(index=False)}
    noise = estimate_noise_sigma_by_chromosome(old_clean, [depth_by_coord[key] for key in old_clean])
    rows = []
    for row in candidates:
        overlap = [(max(a, row.st-1), min(b, row.nd), nid) for a, b, nid in intervals.get(row.chr, [])
                   if max(a, row.st-1) < min(b, row.nd)]
        ids = sorted({nid for a, b, nid in overlap})
        spans = sorted((a, b) for a, b, nid in overlap)
        merged = []
        for a, b in spans:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        rows.append(dict(chrom=row.chr, start=int(row.st), end=int(row.nd), depth=float(row.meandepth),
                         nclose_gate=bool(ids), nclose_ids=';'.join(ids),
                         nclose_overlap_bp=sum(b-a for a, b in merged),
                         noise_sigma=float(noise[row.chr])))
    result = pd.DataFrame(rows)
    result['direct_nclose_gate'] = result.nclose_gate
    # A boundary NClose may support a long amplified segment. Propagate its
    # eligibility only within a contiguous run of already-high depth bins,
    # not along a whole chromosome-spanning model path.
    for chrom, group in result.groupby('chrom', sort=False):
        blocks = []
        for row in group.sort_values('start').itertuples():
            if not blocks or row.start != result.loc[blocks[-1][-1], 'end'] + 1:
                blocks.append([])
            blocks[-1].append(row.Index)
        for indices in blocks:
            ids = sorted({nid for text in result.loc[indices, 'nclose_ids'] for nid in text.split(';') if nid})
            result.loc[indices, 'block_start'] = int(result.loc[indices, 'start'].min())
            result.loc[indices, 'block_end'] = int(result.loc[indices, 'end'].max())
            result.loc[indices, 'block_nclose_ids'] = ';'.join(ids)
            result.loc[indices, 'nclose_gate'] = bool(ids)
    return result


def robust_loss_gradient(matrix, target, x, high_mask, tau):
    residual = matrix @ x - target
    derivative = residual.copy()
    derivative[high_mask] = np.maximum(residual[high_mask], -tau)
    loss = 0.5 * float(residual[~high_mask] @ residual[~high_mask])
    rh = residual[high_mask]
    loss += float(np.where(rh >= -tau, 0.5*rh*rh, -tau*rh-0.5*tau*tau).sum())
    return loss, matrix.T @ derivative


def _warm_weighted_nnls(a, b, seed, tolerance=1e-12):
    admitted = seed > 0
    norms = np.linalg.norm(a, axis=0)
    for iteration in range(128):
        cols = np.flatnonzero(admitted)
        x = np.zeros(a.shape[1])
        if len(cols):
            x[cols] = _full_nnls(a[:, cols], b)
        g = a.T @ (a @ x - b)
        violation = np.where(x > 0, np.abs(g), np.maximum(-g, 0))
        v = np.divide(violation, norms, out=np.zeros_like(g), where=norms > 0)
        if float(v.max()) <= tolerance:
            return x
        candidates = np.flatnonzero((~admitted) & (g < 0) & (v > tolerance))
        if not len(candidates):
            raise RuntimeError(f'Weighted NNLS inner KKT failed on admitted columns: {v.max()}')
        admitted[candidates[np.argsort(-v[candidates])[:32]]] = True
    raise RuntimeError('Weighted NNLS inner iteration limit')


def solve_robust(matrix, target, high_mask, tau, initial=None, tolerance=1e-10):
    """Convex one-sided Huber; clean rows retain their exact quadratic loss.

    Normalize columns/target for conditioning. KKT is recomputed for this
    objective; it is not the ordinary NNLS KKT certificate.
    """
    high_mask = np.asarray(high_mask, dtype=bool)
    if high_mask.shape != np.shape(target):
        raise ValueError('High-depth mask must match the target vector')
    if not high_mask.any():
        weights, diagnostics = solve_depth_nnls(matrix, target, return_diagnostics=True)
        diagnostics.update(high_depth_fallback='no_high_depth_rows',
                           unexplained_high_depth_sum=0., saturated_high_depth_bins=0)
        return weights, diagnostics
    tau = np.asarray(tau, dtype=float)
    if tau.shape != (int(high_mask.sum()),) or not np.isfinite(tau).all() or np.any(tau <= 0):
        raise ValueError('High-depth thresholds must be finite, positive, and match high-depth rows')
    clean_matrix, clean_target = aggregate_identical_depth_rows(matrix[~high_mask], target[~high_mask])
    a = np.vstack([clean_matrix, np.asarray(matrix[high_mask], dtype=float)])
    b = np.r_[clean_target, target[high_mask]]
    high = np.arange(len(b)) >= len(clean_target)
    column_norm = np.linalg.norm(a, axis=0)
    valid = column_norm > 0
    b_norm = float(np.linalg.norm(target)) or 1.
    a = np.asfortranarray(a[:, valid] / column_norm[valid])
    b = b / b_norm
    scaled_tau = np.asarray(tau, dtype=float) / b_norm
    x0 = np.zeros(int(valid.sum())) if initial is None else np.asarray(initial)[valid] * column_norm[valid] / b_norm
    history = []
    def objective(x):
        return robust_loss_gradient(a, b, x, high, scaled_tau)
    def callback(x):
        loss, g = objective(x)
        v = np.where(x > 0, np.abs(g), np.maximum(-g, 0))
        history.append((len(history)+1, loss, float(v.max())))
        if len(history) % 200 == 0:
            logging.info('Experimental robust iteration %d, scaled KKT %.3g', len(history), v.max())
    result = minimize(objective, x0, jac=True, method='L-BFGS-B', bounds=[(0., None)]*len(x0),
        callback=callback, options=dict(gtol=tolerance/10, ftol=0., maxiter=20000, maxls=50, maxcor=50))
    x = result.x
    loss, gradient = objective(x)
    violation = np.where(x > 0, np.abs(gradient), np.maximum(-gradient, 0))
    maximum = float(violation.max())
    polish_history = []
    for outer in range(100):
        if maximum <= tolerance:
            break
        residual = a @ x - b
        row_weight = np.ones(len(b))
        row_weight[high] = np.minimum(1., scaled_tau/np.maximum(-residual[high], scaled_tau))
        root_weight = np.sqrt(row_weight)
        trial = _warm_weighted_nnls(a*root_weight[:, None], b*root_weight, x)
        # One-sided loss needs a line search if an underfit row crosses into
        # overprediction; do not assume symmetric-Huber IRLS majorization.
        step = 1.
        allowance = 64*np.finfo(float).eps*max(1.,abs(loss))
        for backtrack in range(40):
            candidate = (1-step)*x+step*trial
            candidate_loss, candidate_gradient = objective(candidate)
            if candidate_loss <= loss+allowance:
                break
            step *= .5
        else:
            raise RuntimeError('Robust weighted-NNLS polish could not decrease objective')
        x, loss, gradient = candidate, candidate_loss, candidate_gradient
        violation = np.where(x > 0, np.abs(gradient), np.maximum(-gradient, 0))
        maximum = float(violation.max())
        polish_history.append(dict(iteration=outer+1,kkt_scaled_max=maximum,step=step))
        logging.info('Experimental robust polish %d, scaled KKT %.3g', outer+1, maximum)
    weights = np.zeros(matrix.shape[1]); weights[valid] = x*b_norm/column_norm[valid]
    predicted = predict_depth_from_weights(matrix, weights)
    unexplained = np.maximum(target[high_mask] - predicted[high_mask] - tau, 0)
    diagnostics = dict(solver='scipy.optimize.minimize + scipy.optimize.nnls', strategy='one_sided_huber',
        working_set_rounds=int(result.nit), max_working_columns=int(valid.sum()), working_set_fallback=None,
        kkt_scaled_max=maximum, kkt_tolerance=tolerance, kkt_matrix='clean exact aggregation plus individual high-depth rows',
        convergence_check='one-sided-Huber KKT', converged=bool(maximum <= tolerance),
        solver_message=str(result.message), observation_rows=int(matrix.shape[0]), aggregated_rows=int(len(b)),
        feature_columns=int(matrix.shape[1]), robust_objective_scaled=loss,
        unexplained_high_depth_sum=float(unexplained.sum()), saturated_high_depth_bins=int((unexplained>0).sum()),
        sse=float(np.sum((predicted-target)**2)), relative_error=float(np.linalg.norm(predicted-target)/np.linalg.norm(target)),
        min_weight=float(weights.min()), positive_weights=int(np.count_nonzero(weights)),
        history=history, weighted_nnls_polish=polish_history)
    if not np.isfinite(maximum) or maximum > tolerance:
        raise RuntimeError(f'Robust fit failed its own KKT check: {maximum}, {result.message}')
    return weights, diagnostics
