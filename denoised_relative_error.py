"""Denoise SKYPE depth targets before reporting their relative fit error.

The public entry point is :func:`calculate_denoised_relative_error`.  It uses
only NumPy so that it can run in the existing ``skype`` environment without an
additional dependency.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence

import numpy as np


DEPTH_BIN_STEP = 100_000
NORMAL_MAD_QUANTILE = 0.6744897501960817
MIN_CHROM_NOISE_DIFFERENCES = 20


# Why TV, and where does lambda / noise_sigma = 3 come from?
# ---------------------------------------------------------------------------
# This is an empirical CASTLE HiFi calibration, not a constant prescribed by
# the total-variation paper.  Condat (2013, doi:10.1109/LSP.2013.2278339)
# supplies the exact O(n) solver used below for a *given* lambda; that paper
# explicitly treats selection of lambda as a separate problem.
#
# Calibration data and target leakage control:
#   * H1437, H2009, HCC1937, and HCC1954 CASTLE HiFi depth were each aligned to
#     hs1 and hg38.  The two reference alignments used their own assembly-
#     specific chr_filt_st_list, giving eight equally weighted depth tracks.
#   * There was no clean/noiseless depth truth.  Comparing a denoised track
#     directly with its noisy input would therefore select the identity
#     transform.  Instead, deterministic 10-fold masked prediction was used:
#     every eligible bin was hidden once, filled only from visible neighbours
#     in the same contiguous run, denoised, and scored at the hidden bin.
#   * The selection score was RRMSE = sqrt(sum((prediction - observed)^2) /
#     sum(observed^2)), averaged with equal weight across all eight tracks.
#     Forty-four raw/median/Gaussian/Savitzky-Golay/TV/Haar configurations were
#     compared.  The TV grid was c={0.25, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4,
#     4.5, 5, 6, 8, 12}.  TV c=3 had the lowest mean masked RRMSE (0.098588),
#     versus 0.106095 for the raw/interpolation baseline.  Nearby TV values
#     formed a shallow but ordered optimum: c=3.5 gave 0.098615, c=2.5 gave
#     0.098679, and c=4 gave 0.098711.  Leave-one-sample-out selection returned
#     c=3 for three held-out samples and c=3.5 when HCC1954 was held out.
#
# Why total variation was retained:
#   * Copy-number depth is approximately piecewise constant with real, sharp
#     breakpoints.  The L1 penalty on adjacent differences removes small local
#     fluctuations while preserving steps better than convolutional smoothing
#     such as a Gaussian filter.
#   * A dimensionless multiplier is more portable than a fixed raw lambda.
#     For each chromosome, noise_sigma is estimated robustly from adjacent
#     positive-depth differences: MAD(diff) / 0.67449 / sqrt(2).  Dividing by
#     sqrt(2) converts the scale of differences of two independent noisy bins
#     back to a per-bin noise scale.  Chromosomes with fewer than 20 usable
#     differences fall back to the track-wide estimate.
#   * Denoising is restarted at chromosome boundaries and whenever successive
#     starts are not 100 kb apart.  It can therefore never smooth across a
#     centromere/satellite/high-depth gap removed by chr_filt_st_list.
#
# Thus c=3 means lambda = 3 * chromosome_noise_sigma in
#   0.5 * ||x - y||_2^2 + lambda * sum_i |x[i+1] - x[i]|.
# It is the common setting selected for these CASTLE HiFi 100-kb tracks.  It
# should be recalibrated if the assay, bin width, normalization, or filtering
# policy changes materially.
TV_LAMBDA_OVER_NOISE_SIGMA = 3.0


def _depth_vector(values: Sequence[float] | np.ndarray, name: str) -> np.ndarray:
    result = np.asarray(values, dtype=float)
    if result.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    if not np.isfinite(result).all():
        raise ValueError(f"{name} contains a non-finite value")
    return result


def _contiguous_run_slices(
    coordinates: Sequence[tuple[str, int]],
) -> list[tuple[str, slice]]:
    """Split coordinates without crossing a chromosome or filtered 100-kb gap."""
    if len(coordinates) == 0:
        return []

    runs: list[tuple[str, slice]] = []
    run_start = 0
    previous_chrom = str(coordinates[0][0])
    previous_start = int(coordinates[0][1])
    for index, (chrom_value, start_value) in enumerate(coordinates[1:], start=1):
        chrom = str(chrom_value)
        start = int(start_value)
        if chrom != previous_chrom or start - previous_start != DEPTH_BIN_STEP:
            runs.append((previous_chrom, slice(run_start, index)))
            run_start = index
        previous_chrom = chrom
        previous_start = start
    runs.append((previous_chrom, slice(run_start, len(coordinates))))
    return runs


def _robust_normal_sigma(values: np.ndarray) -> float:
    if values.size == 0:
        return float("nan")
    center = np.median(values)
    return float(np.median(np.abs(values - center)) / NORMAL_MAD_QUANTILE)


def estimate_noise_sigma_by_chromosome(
    coordinates: Sequence[tuple[str, int]],
    observed_depth: Sequence[float] | np.ndarray,
) -> dict[str, float]:
    """Estimate per-bin noise from adjacent positive bins in contiguous runs."""
    depth = _depth_vector(observed_depth, "observed_depth")
    if len(coordinates) != depth.size:
        raise ValueError(
            "coordinate/depth length mismatch: "
            f"{len(coordinates)} != {depth.size}"
        )

    differences_by_chromosome: dict[str, list[float]] = defaultdict(list)
    all_differences: list[float] = []
    for chrom, run_slice in _contiguous_run_slices(coordinates):
        run = depth[run_slice]
        if run.size < 2:
            continue
        # A zero-depth bin can be a mapping dropout or a true deep deletion;
        # neither is a sample from symmetric measurement noise.  Only adjacent
        # pairs for which both depths are positive contribute to sigma.
        positive_pair = (run[:-1] > 0) & (run[1:] > 0)
        differences = np.diff(run)[positive_pair]
        differences_by_chromosome[chrom].extend(differences.tolist())
        all_differences.extend(differences.tolist())

    global_differences = np.asarray(all_differences, dtype=float)
    global_sigma = _robust_normal_sigma(global_differences) / np.sqrt(2.0)
    if not np.isfinite(global_sigma) or global_sigma <= 0:
        # A constant/all-zero track has no estimable adjacent noise.  One depth
        # unit is a conservative, deterministic fallback that keeps TV valid.
        global_sigma = 1.0

    result: dict[str, float] = {"__global__": float(global_sigma)}
    for chrom in dict.fromkeys(str(chrom) for chrom, _ in coordinates):
        differences = np.asarray(
            differences_by_chromosome.get(chrom, []), dtype=float
        )
        sigma = _robust_normal_sigma(differences) / np.sqrt(2.0)
        if (
            differences.size < MIN_CHROM_NOISE_DIFFERENCES
            or not np.isfinite(sigma)
            or sigma <= 0
        ):
            sigma = global_sigma
        result[chrom] = float(sigma)
    return result


def tv1d_denoise(observed: Sequence[float] | np.ndarray, lam: float) -> np.ndarray:
    """Return argmin_x 0.5*||x-y||^2 + lam*sum(abs(diff(x)))."""
    y = _depth_vector(observed, "observed")
    n = y.size
    if n <= 1 or lam <= 0:
        return y.copy()

    # Direct O(n) 1-D TV algorithm from Condat (2013).  This implementation is
    # kept local so the production SKYPE environment needs no denoising package.
    output = np.empty(n, dtype=float)
    k = 0
    k0 = 0
    umin = lam
    umax = -lam
    vmin = y[0] - lam
    vmax = y[0] + lam
    kplus = 0
    kminus = 0
    two_lam = 2.0 * lam
    minus_lam = -lam

    while True:
        while k == n - 1:
            if umin < 0.0:
                output[k0 : kminus + 1] = vmin
                k0 = kminus + 1
                k = k0
                kminus = k0
                vmin = y[k0]
                umin = lam
                umax = vmin + lam - vmax
            elif umax > 0.0:
                output[k0 : kplus + 1] = vmax
                k0 = kplus + 1
                k = k0
                kplus = k0
                vmax = y[k0]
                umax = minus_lam
                umin = vmax - lam - vmin
            else:
                vmin += umin / (k - k0 + 1)
                output[k0 : k + 1] = vmin
                return output

        umin += y[k + 1] - vmin
        if umin < minus_lam:
            output[k0 : kminus + 1] = vmin
            k0 = kminus + 1
            k = k0
            kplus = k0
            kminus = k0
            vmin = y[k0]
            vmax = vmin + two_lam
            umin = lam
            umax = minus_lam
            continue

        umax += y[k + 1] - vmax
        if umax > lam:
            output[k0 : kplus + 1] = vmax
            k0 = kplus + 1
            k = k0
            kplus = k0
            kminus = k0
            vmax = y[k0]
            vmin = vmax - two_lam
            umin = lam
            umax = minus_lam
            continue

        k += 1
        if umin >= lam:
            kminus = k
            vmin += (umin - lam) / (kminus - k0 + 1)
            umin = lam
        if umax <= minus_lam:
            kplus = k
            vmax += (umax + lam) / (kplus - k0 + 1)
            umax = minus_lam


def denoise_depth_target(
    coordinates: Sequence[tuple[str, int]],
    observed_depth: Sequence[float] | np.ndarray,
    lambda_over_noise_sigma: float = TV_LAMBDA_OVER_NOISE_SIGMA,
) -> np.ndarray:
    """TV-denoise an observed target independently within every genomic run."""
    depth = _depth_vector(observed_depth, "observed_depth")
    if len(coordinates) != depth.size:
        raise ValueError(
            "coordinate/depth length mismatch: "
            f"{len(coordinates)} != {depth.size}"
        )
    if not np.isfinite(lambda_over_noise_sigma) or lambda_over_noise_sigma <= 0:
        raise ValueError("lambda_over_noise_sigma must be finite and positive")

    noise = estimate_noise_sigma_by_chromosome(coordinates, depth)
    denoised = np.empty(depth.size, dtype=float)
    for chrom, run_slice in _contiguous_run_slices(coordinates):
        lam = lambda_over_noise_sigma * noise[chrom]
        denoised[run_slice] = tv1d_denoise(depth[run_slice], lam)
    # Depth cannot be negative.  TV of non-negative input is theoretically
    # non-negative, but clipping also protects against floating-point residue.
    return np.maximum(denoised, 0.0)


def calculate_denoised_relative_error(
    coordinates: Sequence[tuple[str, int]],
    observed_depth: Sequence[float] | np.ndarray,
    predicted_depth: Sequence[float] | np.ndarray,
    lambda_over_noise_sigma: float = TV_LAMBDA_OVER_NOISE_SIGMA,
) -> tuple[float, float, np.ndarray]:
    """Return absolute error, relative error, and the TV-denoised target.

    The reported relative metric is
    ``||predicted - denoised_observed||_2 / ||denoised_observed||_2``.
    The SKYPE prediction is not denoised.
    """
    predicted = _depth_vector(predicted_depth, "predicted_depth")
    denoised = denoise_depth_target(
        coordinates, observed_depth, lambda_over_noise_sigma
    )
    if predicted.size != denoised.size:
        raise ValueError(
            "prediction/depth length mismatch: "
            f"{predicted.size} != {denoised.size}"
        )

    denominator = float(np.linalg.norm(denoised))
    if denominator == 0:
        raise ValueError("cannot calculate relative error for a zero-norm target")
    absolute_error = float(np.linalg.norm(predicted - denoised))
    return absolute_error, absolute_error / denominator, denoised
