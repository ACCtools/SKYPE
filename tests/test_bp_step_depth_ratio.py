from __future__ import annotations

import collections
import tempfile
import unittest

import numpy as np
import vcfpy

from bp_step_depth_ratio import (
    BP_STEP_DEPTH_RATIO_B,
    BP_STEP_DEPTH_RATIO_PREDICT_B,
    BreakpointStepDepthRatio,
    add_ratio_info,
    bnd_expected_high_side,
)


class BreakpointStepDepthRatioTests(unittest.TestCase):
    def setUp(self):
        coordinates = [
            ("chr1", 100_000),
            ("chr1", 300_000),
            ("chr1", 500_000),
            ("chr1", 700_000),
            ("chr1", 900_000),
            ("chr2", 100_000),
            ("chr2", 300_000),
        ]
        self.calculator = BreakpointStepDepthRatio(
            coordinates,
            observed_depth=np.array([10, 14, 20, 22, 24, 3, 5], dtype=float),
            predicted_depth=np.array([8, 12, 18, 26, 30, 4, 6], dtype=float),
            censat_intervals={"chr1": [(850_000, 950_000)]},
        )

    def test_signed_ratio_uses_500kb_window_and_raw_event_depth(self):
        ratios = self.calculator.pair(
            [
                ("chr1", 500_000, "left"),
                ("chr1", 700_000, "right"),
            ],
            event_depth=4,
        )

        self.assertAlmostEqual(ratios[BP_STEP_DEPTH_RATIO_B][0], -2.5)
        self.assertAlmostEqual(ratios[BP_STEP_DEPTH_RATIO_B][1], 1.5)
        self.assertAlmostEqual(
            ratios[BP_STEP_DEPTH_RATIO_PREDICT_B][0], -11.0 / 3.0
        )
        self.assertAlmostEqual(
            ratios[BP_STEP_DEPTH_RATIO_PREDICT_B][1], 3.25
        )

    def test_missing_windows_censat_and_zero_depth_are_missing(self):
        self.assertEqual(
            self.calculator.endpoint("chr2", 100_000, "left", 4),
            (None, None),
        )
        self.assertEqual(
            self.calculator.endpoint("chr1", 900_000, "left", 4),
            (None, None),
        )
        self.assertEqual(
            self.calculator.endpoint("chr1", 500_000, "left", 0),
            (None, None),
        )

    def test_expected_sides_match_bnd_roles(self):
        self.assertEqual(bnd_expected_high_side("+", "exit"), "left")
        self.assertEqual(bnd_expected_high_side("-", "exit"), "right")
        self.assertEqual(bnd_expected_high_side("+", "entry"), "right")
        self.assertEqual(bnd_expected_high_side("-", "entry"), "left")

    def test_second_bnd_mate_reverses_arrays(self):
        ratio_info = {
            BP_STEP_DEPTH_RATIO_B: [1.0, 2.0],
            BP_STEP_DEPTH_RATIO_PREDICT_B: [3.0, 4.0],
        }
        first = collections.OrderedDict()
        second = collections.OrderedDict()
        add_ratio_info(first, ratio_info)
        add_ratio_info(second, ratio_info, reverse=True)

        self.assertEqual(first[BP_STEP_DEPTH_RATIO_B], [1.0, 2.0])
        self.assertEqual(second[BP_STEP_DEPTH_RATIO_B], [2.0, 1.0])
        self.assertEqual(second[BP_STEP_DEPTH_RATIO_PREDICT_B], [4.0, 3.0])

    def test_vcfpy_serializes_partial_missing_float_array(self):
        with tempfile.NamedTemporaryFile(suffix=".vcf") as output:
            header = vcfpy.Header(lines=[vcfpy.HeaderLine("fileformat", "VCFv4.3")])
            header.add_info_line(collections.OrderedDict([
                ("ID", BP_STEP_DEPTH_RATIO_B),
                ("Number", 2),
                ("Type", "Float"),
                ("Description", "test"),
            ]))
            with vcfpy.Writer.from_path(output.name, header) as writer:
                writer.write_record(vcfpy.Record(
                    CHROM="chr1",
                    POS=1,
                    ID=["test"],
                    REF="N",
                    ALT=[vcfpy.SymbolicAllele("DEL")],
                    QUAL=60,
                    FILTER=[],
                    INFO={BP_STEP_DEPTH_RATIO_B: [None, 0.5]},
                ))
            with open(output.name) as handle:
                self.assertIn("BP_STEP_DEPTH_RATIO_B=.,0.5", handle.read())


if __name__ == "__main__":
    unittest.main()
