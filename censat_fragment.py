"""Detect a depth discontinuity across CEN-SAT after separating flank changes.

Depth-stat bins are 1-based inclusive; internal intervals are 0-based half-open.
Only supported calls enter the legacy fragment dictionary. Every other result
is retained in the diagnostic report. Defaults are exploratory, not calibrated
on labelled centromere rearrangements.
"""

from dataclasses import asdict, dataclass
import json
from pathlib import Path

import numpy as np

from denoised_relative_error import estimate_noise_sigma_by_chromosome


@dataclass(frozen=True)
class CensatFragmentPolicy:
    flank_bp: int = 10_000_000
    edge_guard_bp: int = 100_000
    min_segment_bins: int = 5
    compare_bp: int = 2_000_000
    penalty_scale: float = 3.0
    nclose_penalty_ratio: float = 0.75
    min_effect_fraction: float = 0.10
    min_effect_sigma: float = 1.5
    bootstrap_block_bins: int = 5
    bootstrap_repeats: int = 500


def nclose_breakend_positions(contig_data, nclose_nodes):
    """Use junction-facing coordinates, not whole endpoint alignment spans."""
    result = {}
    if contig_data is None or nclose_nodes is None:
        return result
    for pairs in nclose_nodes.values():
        for pair in pairs:
            for order, index in enumerate(pair):
                row = contig_data[index]
                position = row[7 if ((row[4] == '+') ^ (order == 0)) else 8]
                result.setdefault(str(row[5]), set()).add(int(position))
    return {chrom: sorted(positions) for chrom, positions in result.items()}


def _segment(values, starts, sigma, breakends, policy):
    """Penalized piecewise-constant fit; raw NCloses only discount split cost.

    Clamp isolated spikes around a five-bin median before fitting. Persistent
    steps survive this operation. All bin boundaries remain eligible, including
    those without an assembly junction. A hinted split still needs depth gain.
    """
    count = len(values)
    if count < 2 * policy.min_segment_bins:
        return [(0, count)]
    local = np.array([
        np.median(values[max(0, i - 2):min(count, i + 3)])
        for i in range(count)
    ])
    fitted_values = np.clip(values, local - 4 * sigma, local + 4 * sigma)
    sums = np.r_[0.0, np.cumsum(fitted_values)]
    squares = np.r_[0.0, np.cumsum(fitted_values ** 2)]
    penalty = policy.penalty_scale * sigma ** 2 * np.log(max(count, 2))
    step = float(np.median(np.diff(starts)))
    hinted = {int(np.argmin(abs(starts - pos))) for pos in breakends
              if starts[0] - step <= pos <= starts[-1] + step}
    cost = np.full(count + 1, np.inf)
    previous = np.zeros(count + 1, dtype=int)
    cost[0] = 0.0
    for end in range(policy.min_segment_bins, count + 1):
        candidates = np.r_[0, np.arange(policy.min_segment_bins,
                                        end - policy.min_segment_bins + 1)]
        lengths = end - candidates
        residual = (squares[end] - squares[candidates]
                    - (sums[end] - sums[candidates]) ** 2 / lengths)
        split_cost = np.array([
            0.0 if start == 0 else penalty * (
                policy.nclose_penalty_ratio if int(start) in hinted else 1.0
            ) for start in candidates
        ])
        scores = cost[candidates] + np.maximum(residual, 0) + split_cost
        best = int(np.argmin(scores))
        cost[end], previous[end] = scores[best], candidates[best]
    segments = []
    end = count
    while end:
        start = int(previous[end])
        segments.append((start, end))
        end = start
    return list(reversed(segments))


def _adjacent_flank(rows, censat, side, all_censat, sigma, breakends, policy):
    start, end = censat
    if side == 'left':
        lo, hi = max(0, start - policy.flank_bp), start - policy.edge_guard_bp
    else:
        lo, hi = end + policy.edge_guard_bp, end + policy.flank_bp
    # Whole bins only: do not leak censat bases into flank summaries.
    selected = [(a, b, value) for a, b, value in rows
                if a >= lo and b <= hi and np.isfinite(value) and value > 0
                and not any(a < d and c < b for c, d in all_censat)]
    if not selected:
        return {'status': 'insufficient_flank', 'segments': []}, np.array([])
    starts, ends, values = np.asarray(selected, dtype=float).T
    step = float(np.median(ends - starts))
    # A masked/missing stretch is a boundary, not permission to skip outward.
    if (hi - ends[-1] if side == 'left' else starts[0] - lo) > step:
        return {'status': 'edge_gap', 'segments': []}, np.array([])
    gaps = np.flatnonzero(starts[1:] - ends[:-1] > 1) + 1
    if len(gaps):
        keep = slice(int(gaps[-1]), None) if side == 'left' else slice(None, int(gaps[0]))
        starts, ends, values = starts[keep], ends[keep], values[keep]
    segments = _segment(values, starts, sigma, breakends, policy)
    records = [dict(start=int(starts[a]), end=int(ends[b - 1]),
                    bins=b - a, depth=float(np.median(values[a:b])))
               for a, b in segments]
    a, b = segments[-1] if side == 'left' else segments[0]
    indices = np.arange(a, b)
    if side == 'left':
        indices = indices[ends[b - 1] - starts[indices] <= policy.compare_bp][::-1]
    else:
        indices = indices[ends[indices] - starts[a] <= policy.compare_bp]
    nearby = values[indices]
    record = {'status': 'ok' if len(nearby) >= policy.min_segment_bins else 'short_plateau',
              'segments': records, 'bins': len(nearby), 'bin_bp': int(step),
              'start': int(min(starts[indices])) if len(indices) else None,
              'end': int(max(ends[indices])) if len(indices) else None,
              'depth': float(np.median(nearby)) if len(nearby) else None}
    return record, nearby


def _bootstrap_differences(left, right, policy, rng):
    def resample(values):
        # Moving blocks retain short-range depth correlation. The analytic
        # noise floor below handles plateaus too short to vary under blocks.
        block = min(policy.bootstrap_block_bins, len(values))
        count = int(np.ceil(len(values) / block))
        starts = rng.integers(len(values), size=(policy.bootstrap_repeats, count))
        indices = (starts[..., None] + np.arange(block)) % len(values)
        indices = indices.reshape(policy.bootstrap_repeats, -1)[:, :len(values)]
        return np.median(values[indices], axis=1)
    return resample(right) - resample(left)


def detect_censat_fragments(depth_df, censat_by_chrom, chromosome_lengths,
                            breakends=None, policy=None):
    """Return (supported legacy metadata, diagnostics for all CEN-SAT loci)."""
    policy = policy or CensatFragmentPolicy()
    breakends = breakends or {}
    rows_by_chrom, coordinates, noise_depth = {}, [], []
    for row in depth_df.itertuples(index=False):
        chrom = str(row.chr)
        if chrom == 'chrM':
            continue
        start, end, value = int(row.st) - 1, int(row.nd), float(row.meandepth)
        rows_by_chrom.setdefault(chrom, []).append((start, end, value))
        if np.isfinite(value) and value > 0 and not any(
            start < b and a < end for a, b in censat_by_chrom.get(chrom, ())
        ):
            coordinates.append((chrom, start))
            noise_depth.append(value)
    for rows in rows_by_chrom.values():
        rows.sort()
    # Sort the noise input too, so input file ordering cannot merge chromosomes.
    ordered = sorted(zip(coordinates, noise_depth))
    coordinates = [c for c, _ in ordered]
    noise_depth = np.array([v for _, v in ordered])
    median_depth = float(np.median(noise_depth)) if len(noise_depth) else 0.0
    noise = estimate_noise_sigma_by_chromosome(coordinates, noise_depth)
    accepted, diagnostics = {}, []
    for chrom, intervals in censat_by_chrom.items():
        for locus, (start, end) in enumerate(intervals):
            sigma = float(noise.get(chrom, noise['__global__']))
            minimum = max(policy.min_effect_fraction * median_depth,
                          policy.min_effect_sigma * sigma)
            result = dict(chrom=chrom, censat_start=int(start), censat_end=int(end),
                          locus=locus, status='uncertain', reason='', noise_sigma=sigma,
                          minimum_effect=minimum, median_depth=median_depth)
            diagnostics.append(result)
            length = chromosome_lengths.get(chrom, 0)
            if start <= 0 or end >= length - 1:
                result.update(status='not_evaluable', reason='terminal_censat')
                continue
            if chrom not in rows_by_chrom:
                result.update(status='not_evaluable', reason='missing_depth')
                continue
            sides = [_adjacent_flank(rows_by_chrom[chrom], (start, end), side,
                                    intervals, sigma, breakends.get(chrom, ()), policy)
                     for side in ('left', 'right')]
            (left_info, left), (right_info, right) = sides
            result.update(left=left_info, right=right_info,
                          nearby_breakends=[int(p) for p in breakends.get(chrom, ())
                                            if start - policy.flank_bp <= p <= end + policy.flank_bp])
            if left_info['status'] != 'ok' or right_info['status'] != 'ok':
                result['reason'] = 'insufficient_adjacent_plateau'
                continue
            difference = float(np.median(right) - np.median(left))
            boot = _bootstrap_differences(left, right, policy, np.random.default_rng(0))
            low, high = np.quantile(boot, [.025, .975])
            floor = 1.96 * 1.2533 * sigma * np.sqrt(1 / len(left) + 1 / len(right))
            low, high = min(float(low), difference - floor), max(float(high), difference + floor)
            window_differences = []
            for bp in (500_000, 1_000_000, policy.compare_bp):
                nl = min(len(left), max(policy.min_segment_bins, bp // left_info['bin_bp']))
                nr = min(len(right), max(policy.min_segment_bins, bp // right_info['bin_bp']))
                window_differences.append(float(np.median(right[:nr]) - np.median(left[:nl])))
            stable = all(value * difference > 0 and abs(value) >= minimum / 2
                         for value in window_differences)
            result.update(depth_diff=difference, ci95=[low, high],
                          window_differences=window_differences, stable=stable)
            if low <= 0 <= high:
                result.update(status='unsupported', reason='no_resolved_edge_difference')
            elif abs(difference) < minimum:
                result['reason'] = 'small_edge_difference'
            elif not stable:
                result['reason'] = 'window_sensitive_difference'
            else:
                result.update(status='supported', reason='stable_adjacent_depth_step')
    # The downstream model has one fragment per chromosome. Conflicting or
    # multiple supported loci cannot be represented by silently choosing one.
    for chrom in censat_by_chrom:
        calls = [r for r in diagnostics if r['chrom'] == chrom and r['status'] == 'supported']
        if len(calls) > 1:
            for result in calls:
                result.update(status='uncertain', reason='multiple_supported_censat_loci')
        elif calls:
            result = calls[0]
            accepted[chrom] = dict(dir=result['depth_diff'] > 0,
                                  mid=(result['censat_start'] + result['censat_end']) // 2,
                                  depth_diff=abs(result['depth_diff']),
                                  chr_len=int(chromosome_lengths[chrom]))
    return accepted, dict(method='adjacent_plateau_v1', policy=asdict(policy), loci=diagnostics)


def write_censat_diagnostics(prefix, report):
    prefix = Path(prefix)
    prefix.mkdir(parents=True, exist_ok=True)
    (prefix / 'cen_fragment_diagnostics.json').write_text(
        json.dumps(report, indent=2, allow_nan=False) + '\n')
    columns = ('chrom', 'censat_start', 'censat_end', 'status', 'reason',
               'depth_diff', 'minimum_effect', 'noise_sigma', 'ci95')
    with (prefix / 'cen_fragment_diagnostics.tsv').open('w') as handle:
        handle.write('\t'.join(columns) + '\n')
        for result in report['loci']:
            handle.write('\t'.join(str(result.get(key, '')) for key in columns) + '\n')
