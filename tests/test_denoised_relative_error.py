from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


SKYPE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKYPE_ROOT))

from denoised_relative_error import (  # noqa: E402
    calculate_denoised_relative_error,
    denoise_depth_target,
    estimate_noise_sigma_by_chromosome,
    tv1d_denoise,
)


class DenoisedRelativeErrorTests(unittest.TestCase):
    def test_tv_preserves_constant_signal(self) -> None:
        observed = np.full(8, 7.5)
        np.testing.assert_allclose(tv1d_denoise(observed, 3.0), observed)

    def test_gap_and_chromosome_boundaries_are_not_smoothed_across(self) -> None:
        coordinates = [
            ("chr1", 0),
            ("chr1", 100_000),
            ("chr1", 400_000),
            ("chr1", 500_000),
            ("chr2", 0),
            ("chr2", 100_000),
        ]
        observed = np.asarray([1.0, 1.0, 20.0, 20.0, 50.0, 50.0])

        denoised = denoise_depth_target(
            coordinates, observed, lambda_over_noise_sigma=3.0
        )

        np.testing.assert_allclose(denoised, observed)

    def test_sparse_chromosome_uses_global_noise_fallback(self) -> None:
        coordinates = [("chr1", index * 100_000) for index in range(30)]
        coordinates.extend(("chr2", index * 100_000) for index in range(3))
        observed = np.asarray(
            [10.0 + (-1.0) ** index for index in range(30)] + [5.0, 6.0, 5.0]
        )

        noise = estimate_noise_sigma_by_chromosome(coordinates, observed)

        self.assertGreater(noise["__global__"], 0)
        self.assertEqual(noise["chr2"], noise["__global__"])

    def test_metric_uses_denoised_target_norm_as_denominator(self) -> None:
        coordinates = [("chr1", index * 100_000) for index in range(6)]
        observed = np.asarray([10.0, 9.0, 11.0, 10.0, 9.0, 11.0])
        predicted = np.asarray([10.0, 10.0, 10.0, 10.0, 10.0, 10.0])

        absolute, relative, denoised = calculate_denoised_relative_error(
            coordinates, observed, predicted
        )

        expected_absolute = np.linalg.norm(predicted - denoised)
        self.assertAlmostEqual(absolute, expected_absolute)
        self.assertAlmostEqual(relative, expected_absolute / np.linalg.norm(denoised))

    def test_length_mismatch_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "length mismatch"):
            calculate_denoised_relative_error(
                [("chr1", 0)], np.asarray([1.0]), np.asarray([1.0, 2.0])
            )


if __name__ == "__main__":
    unittest.main()
