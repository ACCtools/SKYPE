from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from raw_translocation_prior import (
    calculate_length_normalized_prior_scales,
    load_normal_prior_excluded_chromosomes,
    select_normal_prior_chromosome_columns,
)


class RawTranslocationPriorTests(unittest.TestCase):
    def test_excludes_both_chromosomes_when_either_side_has_no_span(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            report_path = Path(temporary) / "raw_translocation_read_counts.tsv"
            report_path.write_text(
                "chrom_a\tchrom_b\tpoint_a_no_spanning_rawread\t"
                "point_b_no_spanning_rawread\n"
                "chr3\tchr11\tTrue\tFalse\n"
                "chr5\tchr5\tFalse\tTrue\n"
                "chr12\tchr15\t1\t0\n"
            )

            excluded = load_normal_prior_excluded_chromosomes(report_path)

            self.assertEqual(
                excluded,
                {"chr3", "chr11", "chr5", "chr12", "chr15"},
            )

    def test_missing_report_has_no_exclusions(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            report_path = Path(temporary) / "missing.tsv"

            self.assertEqual(
                load_normal_prior_excluded_chromosomes(report_path),
                set(),
            )

    def test_rejects_report_without_no_span_columns(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            report_path = Path(temporary) / "raw_translocation_read_counts.tsv"
            report_path.write_text("chrom_a\tchrom_b\nchr1\tchr2\n")

            with self.assertRaisesRegex(ValueError, "missing required columns"):
                load_normal_prior_excluded_chromosomes(report_path)

    def test_excluded_chromosomes_do_not_receive_prior_columns(self) -> None:
        chromosome_paths = {
            "chr1": ("20_depth", "chr1f_chr1b", "1.paf"),
            "chr2": ("20_depth", "chr2f_chr2b", "1.paf"),
            "chr3": ("20_depth", "chr3f_chr3b", "1.paf"),
        }
        path_columns = {
            chromosome_paths["chr1"]: 7,
            chromosome_paths["chr2"]: 11,
        }

        selected = select_normal_prior_chromosome_columns(
            chromosome_paths,
            path_columns,
            {"chr2"},
        )

        self.assertEqual(selected, [("chr1", 7)])

    def test_length_normalization_equalizes_relative_curvature(self) -> None:
        base_norm_squares = [2000.0, 500.0, 250.0]

        scales, reference_scale, factor = (
            calculate_length_normalized_prior_scales(
                base_norm_squares,
                base_norm_squares,
                depth_row_count=2750,
                strength=0.01,
            )
        )

        relative_strengths = [
            scale * scale / norm_square
            for scale, norm_square in zip(scales, base_norm_squares)
        ]
        self.assertTrue(
            all(
                abs(value - relative_strengths[0]) < 1e-12
                for value in relative_strengths
            )
        )
        self.assertAlmostEqual(sum(scale * scale for scale in scales), 27.5)
        self.assertAlmostEqual(reference_scale ** 2, 27.5 / 3)
        self.assertAlmostEqual(factor, 0.01)

    def test_exclusion_does_not_strengthen_surviving_prior(self) -> None:
        base_norm_squares = [2000.0, 500.0, 250.0]
        all_scales, _, _ = calculate_length_normalized_prior_scales(
            base_norm_squares,
            base_norm_squares,
            depth_row_count=2750,
            strength=0.01,
        )
        surviving_scales, _, _ = calculate_length_normalized_prior_scales(
            base_norm_squares,
            [2000.0, 250.0],
            depth_row_count=2750,
            strength=0.01,
        )

        self.assertAlmostEqual(surviving_scales[0], all_scales[0])
        self.assertAlmostEqual(surviving_scales[1], all_scales[2])
        self.assertLess(
            sum(scale * scale for scale in surviving_scales),
            sum(scale * scale for scale in all_scales),
        )


if __name__ == "__main__":
    unittest.main()
