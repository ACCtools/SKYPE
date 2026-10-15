import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pipeline
from karyotype_mode import (
    KARYOTYPE_MODE_JSON,
    decide_karyotype_mode,
    load_karyotype_mode,
    remove_suffixed_outputs,
    save_karyotype_mode,
)
from virtual_sky_plotting import render_karyotype_diagram


SKYPE_ROOT = Path(__file__).resolve().parents[1]


class KaryotypeModeDecisionTests(unittest.TestCase):
    def test_karyotype_is_kept_up_to_the_nclose_node_limit(self):
        config = decide_karyotype_mode("karyotype", 1000, vcf_input=False)
        self.assertTrue(config["karyotype_mode"])
        self.assertFalse(config["forced_variant"])

    def test_large_graphs_and_vcf_input_force_variant_mode(self):
        for count, vcf_input in ((1002, False), (10, True)):
            with self.subTest(count=count, vcf_input=vcf_input):
                config = decide_karyotype_mode("karyotype", count, vcf_input)
                self.assertEqual(config["mode"], "variant")
                self.assertFalse(config["karyotype_mode"])
                self.assertTrue(config["forced_variant"])
                self.assertEqual(config["requested_mode"], "karyotype")

    def test_requested_variant_is_not_reported_as_forced(self):
        config = decide_karyotype_mode("variant", 10, vcf_input=False)
        self.assertEqual(config["mode"], "variant")
        self.assertFalse(config["forced_variant"])

    def test_missing_decision_means_unsuffixed_outputs_only(self):
        with tempfile.TemporaryDirectory() as prefix:
            with self.assertLogs(level="WARNING"):
                self.assertFalse(load_karyotype_mode(prefix)["karyotype_mode"])
            save_karyotype_mode(prefix, decide_karyotype_mode("karyotype", 10, False))
            self.assertTrue(load_karyotype_mode(prefix)["karyotype_mode"])

    def test_variant_cleanup_removes_only_suffixed_outputs(self):
        with tempfile.TemporaryDirectory() as prefix:
            names = ["weight.npy", "SKYPE_result.bed", "weight_filter.npy",
                     "SKYPE_result_cluster.bed", "SV_call_result_filter.nclose_bnds.tsv",
                     "karyotype_cluster.txt"]
            for name in names:
                Path(prefix, name).write_text("")
            removed = {Path(path).name for path in remove_suffixed_outputs(prefix)}
            self.assertEqual(removed, set(names[2:]))
            self.assertEqual(sorted(p.name for p in Path(prefix).iterdir()),
                             ["SKYPE_result.bed", "weight.npy"])


class CentromereFragmentRowTests(unittest.TestCase):
    def test_fragment_display_adds_isochromosome_row(self):
        with tempfile.TemporaryDirectory() as prefix:
            result = render_karyotype_diagram(
                output_prefix=prefix,
                cell_line="TEST",
                chromosome_lengths={"chr1": 1000, "chr12": 800},
                grouped_norm_data={},
                display_indel={},
                virtual_inv_display=[],
                maxh=1000,
                path_depth_n={},
                path_karyotype={},
                fig_prefix="_cluster",
                fragment_display=[("chr12", "right", 400, 800, 1.234)],
            )
            rows = Path(result["karyotype_path"]).read_text().splitlines()
            self.assertEqual(Path(result["karyotype_path"]).name, "karyotype_cluster.txt")
            self.assertTrue(Path(prefix, "virtual_sky_cluster.png").is_file())
        self.assertEqual(rows, ["karyotype\tdepth", "i(12)(q10)\t1.23N"])


class PipelineStageTests(unittest.TestCase):
    def run_commands(self, prefix, **kwargs):
        calls = []

        def record(args, **call_kwargs):
            calls.append((list(args), call_kwargs))

        options = dict(
            CELL_LINE="CELL", PREFIX=prefix, ctg_paf="ctg.paf", ctg_aln_paf="ctg.aln.paf",
            utg_paf="utg.paf", utg_aln_paf="utg.aln.paf", depth_loc="/depth/CELL.win.stat.gz",
            thread=4, dep_folder="/deps", is_progress=False, skype_force=False, graph_depth=2,
            print_args=True, unitig_fasta="utg.fa", alignasm_ref="ref.fa", chr_fai="ref.fai",
            tel_bed="tel.bed", rpt_bed="rpt.bed", rcs_bed="rcs.bed", cyt_bed="cyt.bed",
            raw_rescue_method="off",
        )
        options.update(kwargs)
        with patch.object(pipeline, "subprocess_print", side_effect=record):
            pipeline.run_skype(**options)
        return calls

    @staticmethod
    def script_names(calls):
        names = []
        for args, _ in calls:
            script = next((arg for arg in args if arg.endswith(".py")), None)
            names.append(Path(script).name if script else args[0])
        return names

    @staticmethod
    def command(calls, script):
        return next((args, kwargs) for args, kwargs in calls
                    if any(arg.endswith(script) for arg in args))

    def test_default_run_restores_filter_cluster_and_virtual_sky_stages(self):
        with tempfile.TemporaryDirectory() as prefix:
            calls = self.run_commands(prefix)
        self.assertEqual(self.script_names(calls)[-5:], [
            "22_save_matrix.py", "23_run_nnls.py", "25_cluster_weight.py",
            "30_virtual_sky.py", "31_depth_analysis.py",
        ])
        args, kwargs = self.command(calls, "25_cluster_weight.py")
        self.assertEqual(kwargs["cwd"], str(SKYPE_ROOT))
        self.assertEqual(args[:5], ["python", "-X", "juliacall-threads=4",
                                    "-X", "juliacall-handle-signals=yes"])
        self.assertEqual(args[-2:], ["--mode", "karyotype"])
        self.assertTrue(all(Path(arg).is_absolute() for arg in args[6:12]))
        args, _ = self.command(calls, "30_virtual_sky.py")
        self.assertEqual(args[-2:], [prefix, "CELL"])

    def test_variant_option_is_consumed_by_the_orchestrator(self):
        with tempfile.TemporaryDirectory() as prefix:
            calls = self.run_commands(prefix, option_skype="--variant_mode --verbose")
        args, _ = self.command(calls, "25_cluster_weight.py")
        self.assertEqual(args[-2:], ["--mode", "variant"])
        preprocess, _ = self.command(calls, "01_Preprocess_NClose.py")
        forwarded = [arg for arg in preprocess if arg.startswith("--option_skype=")]
        self.assertEqual(forwarded, ["--option_skype=--verbose"])

    def test_removed_depth_options_are_rejected_before_every_restart(self):
        with tempfile.TemporaryDirectory() as prefix:
            for start in (0, 22, 25, 31):
                for option in ("--depth-policy legacy", "--depth-policy all",
                               "--depth-policy nclose_l2", "--depth-policy=nclose_huber",
                               "--robust-sigma-multiplier=3"):
                    with self.subTest(start=start, option=option):
                        with self.assertRaisesRegex(pipeline.SkypeArgumentError,
                                                    "Unknown --option_skype argument"):
                            self.run_commands(prefix, skype_start_at=start, option_skype=option)

    def test_restarts_reuse_the_requested_mode(self):
        with tempfile.TemporaryDirectory() as prefix:
            Path(prefix, KARYOTYPE_MODE_JSON).write_text(json.dumps(
                decide_karyotype_mode("variant", 10, False)))
            calls = self.run_commands(prefix, skype_start_at=25)
            self.assertEqual(self.script_names(calls), [
                "25_cluster_weight.py", "30_virtual_sky.py", "31_depth_analysis.py"])
            self.assertEqual(self.command(calls, "25_cluster_weight.py")[0][-1], "variant")
            calls = self.run_commands(prefix, skype_start_at=30, option_skype="--karyotype-mode")
            self.assertEqual(self.script_names(calls), [
                "30_virtual_sky.py", "31_depth_analysis.py"])

    def test_restart_at_25_requires_the_fit_matrix_in_karyotype_mode(self):
        with tempfile.TemporaryDirectory() as prefix:
            for name in ("nclose_nodes.pkl", "nclose_filter_status.pkl", "nclose_path_usage.pkl",
                         "weight.npy", "A_idx_list.pkl", "tot_loc_list.pkl", "23_input.pkl"):
                Path(prefix, name).write_text("")
            with self.assertRaisesRegex(pipeline.SkypeArgumentError, "matrix.h5"):
                self.run_commands(prefix, skype_start_at=25, print_args=False)

    def test_unknown_restart_stage_is_rejected(self):
        with tempfile.TemporaryDirectory() as prefix:
            with self.assertRaisesRegex(pipeline.SkypeArgumentError, "25, 30, 31"):
                self.run_commands(prefix, skype_start_at=26)


if __name__ == "__main__":
    unittest.main()
