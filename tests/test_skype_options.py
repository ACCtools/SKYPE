import contextlib
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from skype_options import load_stage_options, save_stage_options, split_stage_options


def load_cli(filename):
    path = Path(__file__).resolve().parents[1] / filename
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SkypeOptionTests(unittest.TestCase):
    def test_stage01_applies_options_and_resets_the_downstream_handoff(self):
        stage01 = load_cli("01_Preprocess_NClose.py")
        with tempfile.TemporaryDirectory() as prefix:
            argv = ["paf", "fai", "tel", "repeat", "censat", "depth", prefix,
                    "bam", "--vcf_input", "input.vcf"]
            result = SimpleNamespace(contig_data=[], nclose_nodes={})
            with patch.object(stage01, "run_stage01", return_value=result) as run:
                stage01.main(argv + ["--option_skype=--skip_bam_analysis --verbose"])
                self.assertTrue(run.call_args.args[0].skip_bam_analysis)
                self.assertEqual(load_stage_options(prefix, "10"), ["--verbose"])
                stage01.main(argv)
                self.assertFalse(run.call_args.args[0].skip_bam_analysis)
                self.assertEqual(load_stage_options(prefix, "10"), [])

    def test_stage10_reuses_or_replaces_saved_options_before_loading_data(self):
        stage10 = load_cli("10_Graph_Find_Paths.py")
        with tempfile.TemporaryDirectory() as prefix:
            save_stage_options(prefix, [], ["--verbose"])
            argv = ["handoff", "fai", prefix, "-t", "1"]
            for extra, expected_verbose in [([], True), (["--option_skype="], False)]:
                with self.subTest(extra=extra), patch.object(
                    stage10, "_remove_stale_stage10_outputs"
                ) as cleanup, patch.object(
                    stage10, "load_stage10_input", side_effect=RuntimeError("loaded")
                ), self.assertRaisesRegex(RuntimeError, "loaded"):
                    stage10.main(argv + extra)
                cleanup.assert_called_once_with(Path(prefix), expected_verbose)
            # A graph-only restart must not silently ignore preprocessing changes.
            with contextlib.redirect_stderr(io.StringIO()), patch.object(
                stage10, "load_stage10_input"
            ) as load, self.assertRaises(SystemExit):
                stage10.main(argv + ["--option_skype=--skip_bam_analysis"])
            load.assert_not_called()

    def test_split_options_are_partitioned_by_stage(self):
        stage01, stage10 = split_stage_options(
            "--check_nclose_count --nclose_count_vaf_threshold 0.2 "
            "--vcf_filter_pass PASS . "
            "--debug-force-nclose chr1:123:+ chr2:456:- "
            "--debug_force_nclose chr3:789:- chr4:1000:+ "
            "--add_indel_graph "
            "--limit_combinations limits.json"
        )
        self.assertEqual(
            stage01,
            [
                "--check-nclose-count",
                "--nclose-count-vaf-threshold",
                "0.2",
                "--vcf-filter-pass",
                "PASS",
                ".",
                "--debug-force-nclose",
                "chr1:123:+",
                "chr2:456:-",
                "--debug-force-nclose",
                "chr3:789:-",
                "chr4:1000:+",
            ],
        )
        self.assertEqual(
            stage10,
            [
                "--add-indel-graph",
                "--limit-combinations",
                "limits.json",
            ],
        )

    def test_split_options_reject_unknown_owned_and_duplicate_flags(self):
        for value, message in (
            ("--unknown", "Unknown"),
            ("--alt replacement.paf", "controlled"),
            ("--verbose --verbose", "Duplicate"),
        ):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, message):
                    split_stage_options(value)

    def test_debug_force_nclose_requires_exactly_two_values(self):
        for value in (
            "--debug-force-nclose chr1:123:+",
            "--debug-force-nclose=chr1:123:+ chr2:456:-",
        ):
            with self.subTest(value=value):
                with self.assertRaisesRegex(ValueError, "requires 2"):
                    split_stage_options(value)
