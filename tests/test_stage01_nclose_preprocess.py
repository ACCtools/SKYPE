from __future__ import annotations

import importlib.util
import json
import pickle
import subprocess
import sys
import tempfile
import unittest
from collections import Counter, defaultdict
from dataclasses import fields
from pathlib import Path
from unittest.mock import patch

import networkx as nx
import pandas as pd

import nclose_preprocess as np
from nclose_candidate import NCloseCandidate, NCloseRejection
from skype_utils import load_pipeline_input


SKYPE_ROOT = Path(__file__).resolve().parents[1]


def make_node(
    name,
    chrom,
    ref_start,
    ref_end,
    *,
    direction="+",
    index=0,
    terminal_range=None,
    contig_type=1,
):
    ref_len = ref_end - ref_start
    terminal_start, terminal_end = terminal_range or (index, index)
    return (
        name,
        ref_len,
        0,
        ref_len,
        direction,
        chrom,
        5_000_000,
        ref_start,
        ref_end,
        60,
        contig_type,
        terminal_start,
        terminal_end,
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        direction,
        chrom,
        f"0.{index}",
    )


def make_config(prefix):
    return np.Stage01Config(
        paf_file_path="primary.aln.paf",
        reference_fai_path="reference.fai",
        telomere_bed_path="telomere.bed",
        repeat_bed_path="repeat.bed",
        censat_bed_path="censat.bed",
        main_stat_path="depth.win.stat.gz",
        prefix=str(prefix),
        read_bam_path="reads.bam",
        alt_path="unitig.aln.paf",
        original_paf_paths=("primary.paf", "unitig.paf"),
        skip_bam_analysis=True,
    )


class Stage01ContractTests(unittest.TestCase):
    def test_unitig_builder_has_one_small_result_contract(self):
        self.assertEqual(
            [field.name for field in fields(np.NCloseCandidateBuildResult)],
            ["candidates", "raw_nodes", "virtual_contigs", "all_nodes"],
        )
        contig_data = [
            make_node(
                "event",
                "chr1",
                100,
                200,
                index=0,
                terminal_range=(0, 1),
            ),
            make_node(
                "event",
                "chr2",
                300,
                400,
                direction="-",
                index=1,
                terminal_range=(0, 1),
            ),
        ]
        result = np.build_unitig_nclose_candidates(
            np.NCloseCandidateBuildContext(
                contig_data=contig_data,
                bnd_contig={"event"},
                repeat_contig_names=set(),
                repeat_censat_data=defaultdict(list),
                paf_file_paths=("primary.aln.paf",),
                original_paf_paths=("primary.paf",),
                telo_set=set(),
                telo_contig={},
                chr_len={"chr1": 5_000_000, "chr2": 5_000_000},
                original_contig_names=[["event"]],
            )
        )

        self.assertEqual(
            [candidate.identity for candidate in result.candidates],
            [("event", (0, 1))],
        )
        self.assertEqual(dict(result.raw_nodes), {"event": [(0, 1)]})
        self.assertEqual(result.virtual_contigs, {})

    def test_pipeline_stages_can_be_inserted_at_named_boundaries(self):
        pipeline = np.default_stage01_pipeline()
        contig_stage = np.ContigPipelineStage(
            "future_unitig_split",
            lambda _context, state: state,
        )
        nclose_stage = np.make_nclose_filter_stage(
            "future_nclose_filter",
            lambda _context, candidates: (list(candidates), []),
        )

        pipeline = pipeline.with_contig_stage(
            contig_stage,
            after="telomere_children",
        ).with_nclose_stage(
            nclose_stage,
            after="censat_pair",
        )

        self.assertEqual(
            [stage.name for stage in pipeline.contig_stages],
            ["telomere_children", "future_unitig_split"],
        )
        nclose_names = [stage.name for stage in pipeline.nclose_stages]
        self.assertEqual(
            nclose_names[nclose_names.index("censat_pair") + 1],
            "future_nclose_filter",
        )

    def test_filter_stage_rejects_non_rejection_provenance(self):
        candidate = NCloseCandidate("event", (0, 1), origin="test")
        stage = np.make_nclose_filter_stage(
            "invalid_filter",
            lambda _context, candidates: (list(candidates), ["not-a-rejection"]),
        )
        state = np.NClosePipelineState(
            candidates=[candidate],
            contig_data_size=2,
            chr_corr={},
            chr_rev_corr={},
        )
        with self.assertRaisesRegex(TypeError, "non-NCloseRejection"):
            stage.run(None, state)

    def test_assembly_source_uses_only_unitig_paths(self):
        source = np.resolve_pregraph_source(
            np.PregraphSourceMode.CONFIGURED_PAF,
            "missing-primary.aln.paf",
            "unitig.aln.paf",
            ("missing-primary.paf", "unitig.paf"),
            "output",
        )
        self.assertEqual(source.paf_file_paths, ("unitig.aln.paf",))
        self.assertEqual(source.original_paf_paths, ("unitig.paf",))
        self.assertEqual(source.secondary_candidate_paf, "unitig.aln.paf")
        with self.assertRaisesRegex(ValueError, "requires --alt"):
            np.resolve_pregraph_source(
                np.PregraphSourceMode.CONFIGURED_PAF,
                "primary.aln.paf",
                None,
                ("primary.paf",),
                "output",
            )

    def test_stage01_cli_does_not_accept_graph_search_options(self):
        path = SKYPE_ROOT / "01_Preprocess_NClose.py"
        spec = importlib.util.spec_from_file_location("stage01_cli_test", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        destinations = {action.dest for action in module.build_parser()._actions}
        self.assertTrue(
            {"thread", "progress", "vcf_input", "check_nclose_count"}
            <= destinations
        )
        self.assertTrue(
            {"graph_depth", "verbose", "add_indel_graph", "limit_combinations"}
            .isdisjoint(destinations)
        )

    def test_legacy_stage02_writes_the_new_downstream_nclose_contract(self):
        source = (SKYPE_ROOT / "02_Build_Breakend_Graph_Limited.py").read_text(
            encoding="utf-8"
        )
        self.assertIn("save_nclose_nodes(PREFIX, nclose_nodes)", source)
        self.assertNotIn("nclose_chunk_data.pkl", source)

    def test_persistence_writes_exact_handoff_and_lightweight_nclose_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            legacy_path = prefix / "nclose_chunk_data.pkl"
            legacy_path.write_bytes(b"stale")
            config = make_config(prefix)
            source = np.resolve_pregraph_source(
                np.PregraphSourceMode.CONFIGURED_PAF,
                config.paf_file_path,
                config.alt_path,
                config.original_paf_paths,
                config.prefix,
            )
            contig_data = [
                make_node("event", "chr1", 100, 200, index=0),
                make_node("event", "chr2", 300, 400, index=1),
            ]
            nclose_nodes = defaultdict(list, {"event": [(0, 1)]})
            result = np.NCloseBuildResult(
                df=None,
                no_chrY=False,
                repeat_censat_data={},
                chr_len={"chr1": 5_000_000, "chr2": 5_000_000},
                contig_data=contig_data,
                contig_data_size=2,
                chr_corr={},
                chr_rev_corr={},
                telo_contig={},
                telo_node_count=0,
                telo_set=set(),
                rpt_con=set(),
                bnd_contig={"event"},
                raw_nclose_nodes=nclose_nodes,
                nclose_nodes=nclose_nodes,
                vctg_dict={},
                all_nclose_comp={"event": [(0, 1)]},
                uncomp_node_count=2,
                nclose_node_count=2,
                transloc_nclose_pair_count=1,
                indel_exclude_idx_set=set(),
                telo_coverage=Counter(),
                contig_stage_records=(
                    {
                        "name": "telomere_children",
                        "before_nodes": 2,
                        "after_nodes": 2,
                    },
                ),
                nclose_stage_records=(
                    {
                        "name": "initial_rejections",
                        "before_candidates": 2,
                        "after_candidates": 1,
                    },
                ),
                nclose_rejections=(
                    NCloseRejection(
                        NCloseCandidate("removed", (0, 1), origin="test"),
                        "initial_rejections",
                        "synthetic_test",
                    ),
                ),
            )

            handoff_path = np.persist_stage01_outputs(config, source, result)

            with open(handoff_path, "rb") as handle:
                handoff = pickle.load(handle)
            with (prefix / "nclose_nodes.pkl").open("rb") as handle:
                persisted_nclose = pickle.load(handle)
            with (prefix / "ecdna_circuit_data.pkl").open("rb") as handle:
                circuits, _ = pickle.load(handle)

            self.assertEqual(
                set(handoff),
                {"contig_data", "nclose_nodes", "telo_contig"},
            )
            self.assertEqual(handoff["contig_data"], contig_data)
            self.assertEqual(handoff["nclose_nodes"], nclose_nodes)
            self.assertEqual(persisted_nclose, nclose_nodes)
            self.assertEqual(circuits, [])
            self.assertFalse(legacy_path.exists())
            self.assertFalse(load_pipeline_input(prefix)["vcf_input"])
            with (prefix / "paf_file_path.pkl").open("rb") as handle:
                self.assertEqual(pickle.load(handle), ["unitig.aln.paf"])
            summary = json.loads(
                (prefix / "stage01_nclose_summary.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(summary["ignored_primary_paf"], "primary.aln.paf")
            self.assertEqual(summary["effective_unitig_paf"], "unitig.aln.paf")
            self.assertEqual(summary["recorded_rejection_count"], 1)
            rejection_lines = (
                prefix / "stage01_nclose_rejections.tsv"
            ).read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(rejection_lines), 2)
            self.assertIn("initial_rejections", rejection_lines[1])
            self.assertIn("synthetic_test", rejection_lines[1])

    def test_ecdna_circuits_are_unique_sorted_four_node_cycles(self):
        contig_data = [
            make_node("cycle", "chr1", index * 1_000, index * 1_000 + 100, index=index)
            for index in range(4)
        ]
        cycle_graph = nx.DiGraph()
        cycle_nodes = [(np.DIR_FOR, index) for index in range(4)]
        cycle_graph.add_edges_from(
            zip(cycle_nodes, cycle_nodes[1:] + cycle_nodes[:1])
        )

        with patch.object(np, "initialize_inversion_only_graph", return_value={}), patch.object(
            np,
            "make_inversion_nx_graph",
            return_value=cycle_graph,
        ):
            circuits, ecdna_nodes = np.build_ecdna_circuits(
                contig_data,
                {},
                {},
            )

        self.assertEqual(circuits, [(0, 1, 2, 3)])
        self.assertEqual(dict(ecdna_nodes), {})

    def test_fake_unitig_builder_runs_all_postprocessing_and_persistence(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            config = make_config(prefix)
            contig_data = [
                make_node("event", "chr1", 1_000_000, 1_100_000, index=0),
                make_node("event", "chr2", 2_000_000, 2_100_000, index=1),
            ]
            ppc_path = prefix / "primary.aln.paf.ppc.paf"
            ppc_path.write_text(
                "".join(
                    "\t".join(map(str, node)) + "\n" for node in contig_data
                ),
                encoding="utf-8",
            )
            (prefix / "reference.fai").write_text(
                "chr1\t5000000\nchr2\t5000000\nchrY\t5000000\n",
                encoding="utf-8",
            )
            for filename in (
                "telomere.bed",
                "repeat.bed",
                "censat.bed",
                "telomere_connected_list.txt",
            ):
                (prefix / filename).write_text("", encoding="utf-8")

            config = np.Stage01Config(
                paf_file_path="primary.aln.paf",
                reference_fai_path=str(prefix / "reference.fai"),
                telomere_bed_path=str(prefix / "telomere.bed"),
                repeat_bed_path=str(prefix / "repeat.bed"),
                censat_bed_path=str(prefix / "censat.bed"),
                main_stat_path="unused",
                prefix=str(prefix),
                read_bam_path="unused",
                alt_path="unitig.aln.paf",
                original_paf_paths=("primary.paf", "unitig.paf"),
                skip_bam_analysis=True,
            )
            source = np.resolve_pregraph_source(
                np.PregraphSourceMode.CONFIGURED_PAF,
                config.paf_file_path,
                config.alt_path,
                config.original_paf_paths,
                config.prefix,
            )
            context = np.PregraphBuildContext(
                source=source,
                ori_ctg_name_data=[],
                prefix=str(prefix),
                preprocessed_paf_path=str(ppc_path),
                reference_fai_path=config.reference_fai_path,
                telomere_bed_path=config.telomere_bed_path,
                repeat_bed_path=config.repeat_bed_path,
                censat_bed_path=config.censat_bed_path,
                main_stat_path=config.main_stat_path,
                asm2cov={},
                disable_alt_ctg_simple=False,
            )
            depth_df = pd.DataFrame(
                {
                    "chr": ["chr1", "chr2", "chrY"],
                    "st": [0, 0, 0],
                    "nd": [100_000, 100_000, 100_000],
                    "meandepth": [20.0, 20.0, 20.0],
                }
            )

            def fake_builder(_context):
                raw_nodes = defaultdict(list, {"event": [(0, 1)]})
                return np.NCloseCandidateBuildResult(
                    candidates=(
                        NCloseCandidate("event", (0, 1), origin="test"),
                    ),
                    raw_nodes=raw_nodes,
                    virtual_contigs={},
                    all_nodes={("chr1", "chr2"): [(0, 1)]},
                )

            np._activate_stage01_runtime(config)
            result = np.nclose_calc(
                context,
                np.ContigPreprocessResult(Counter(), depth_df, False),
                candidate_builder=fake_builder,
            )
            np.persist_stage01_outputs(config, source, result)

            self.assertEqual(dict(result.nclose_nodes), {"event": [(0, 1)]})
            self.assertTrue((prefix / "compressed_nclose_nodes_list.txt").is_file())
            self.assertTrue((prefix / "nclose_event_catalog.pkl").is_file())
            self.assertTrue((prefix / "ecdna_circuit_data.pkl").is_file())
            self.assertTrue((prefix / "nclose_nodes.pkl").is_file())
            self.assertTrue((prefix / "01_nclose_data.pkl").is_file())

    def test_synthetic_assembly_writes_complete_stage01_state(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix = root / "output"
            prefix.mkdir()
            for filename in (
                "01_nclose_data.pkl",
                "nclose_nodes.pkl",
                "ecdna_circuit_data.pkl",
            ):
                (prefix / filename).write_bytes(b"stale")
            fai_path = SKYPE_ROOT / "public_data/chm13v2.0.fa.fai"
            chrom_lengths = {
                fields[0]: int(fields[1])
                for fields in (
                    line.split("\t")
                    for line in fai_path.read_text(encoding="utf-8").splitlines()
                )
            }
            primary_path = root / "missing-primary.aln.paf"
            unitig_path = root / "unitig.aln.paf"
            unitig_path.write_text(
                "ctg1\t400000\t0\t200000\t+\tchr1\t"
                f"{chrom_lengths['chr1']}\t1000000\t1200000\t190000\t"
                "200000\t60\ttp:A:P\tcs:Z::200000\n"
                "ctg1\t400000\t200000\t400000\t+\tchr2\t"
                f"{chrom_lengths['chr2']}\t2000000\t2200000\t190000\t"
                "200000\t60\ttp:A:P\tcs:Z::200000\n",
                encoding="utf-8",
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(SKYPE_ROOT / "01_Preprocess_NClose.py"),
                    str(primary_path),
                    str(fai_path),
                    str(SKYPE_ROOT / "public_data/chm13v2.0_telomere.bed"),
                    str(SKYPE_ROOT / "public_data/chm13v2.0_repeat.m.bed"),
                    str(SKYPE_ROOT / "public_data/chm13v2.0_censat_v2.1.m.bed"),
                    str(SKYPE_ROOT / "public_data/CHM13.win.stat.gz"),
                    str(prefix),
                    str(root / "unused.bam"),
                    "--alt",
                    str(unitig_path),
                    "--original-paf-loc",
                    str(root / "missing-primary.paf"),
                    str(unitig_path),
                    "--skip-bam-analysis",
                    "-t",
                    "1",
                ],
                cwd=SKYPE_ROOT,
                capture_output=True,
                text=True,
            )

            self.assertEqual(
                result.returncode,
                0,
                msg=result.stdout + "\n" + result.stderr,
            )
            ppc_path = prefix / "missing-primary.aln.paf.ppc.paf"
            self.assertTrue(ppc_path.is_file())
            self.assertIn("ctg1", ppc_path.read_text(encoding="utf-8"))
            self.assertTrue((prefix / "telomere_connected_list.txt").is_file())
            self.assertTrue((prefix / "01_nclose_data.pkl").is_file())
            self.assertTrue((prefix / "nclose_nodes.pkl").is_file())
            self.assertTrue((prefix / "ecdna_circuit_data.pkl").is_file())
            self.assertTrue((prefix / "stage01_nclose_summary.json").is_file())
            self.assertTrue((prefix / "stage01_nclose_rejections.tsv").is_file())
            with (prefix / "paf_file_path.pkl").open("rb") as handle:
                self.assertEqual(pickle.load(handle), [str(unitig_path)])
            stage_summary = json.loads(
                (prefix / "stage01_nclose_summary.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                [stage["name"] for stage in stage_summary["contig_stages"]],
                ["telomere_children"],
            )
            self.assertEqual(
                [stage["name"] for stage in stage_summary["nclose_stages"]],
                [
                    "initial_rejections",
                    "censat_pair",
                    "censat_fragment_direction",
                    "simple_alt_preference",
                    "combined_censat_noncensat",
                    "subtelomeric_orientation",
                    "offset_direction",
                    "raw_count_vaf",
                    "raw_translocation_artifact",
                    "raw_virtual_inversion",
                    "user_exclusion",
                ],
            )
            handoff = load_pipeline_input(prefix)
            self.assertFalse(handoff["vcf_input"])
            with (prefix / "01_nclose_data.pkl").open("rb") as handle:
                stage10_input = pickle.load(handle)
            self.assertEqual(
                set(stage10_input),
                {"contig_data", "nclose_nodes", "telo_contig"},
            )
            self.assertIsInstance(stage10_input["nclose_nodes"], dict)

            limit_path = root / "limit_combinations.json"
            limit_path.write_text(
                json.dumps({"limit_combinations": [1, 0]}),
                encoding="utf-8",
            )
            graph_result = subprocess.run(
                [
                    sys.executable,
                    str(SKYPE_ROOT / "10_Graph_Find_Paths.py"),
                    str(prefix / "01_nclose_data.pkl"),
                    str(fai_path),
                    str(prefix),
                    "-t",
                    "1",
                    "-d",
                    "0",
                    "--limit-combinations",
                    str(limit_path),
                ],
                cwd=SKYPE_ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                graph_result.returncode,
                0,
                msg=graph_result.stdout + "\n" + graph_result.stderr,
            )
            self.assertTrue((prefix / "path_data.pkl").is_file())
            self.assertTrue((prefix / "limit_combinations.json").is_file())
            self.assertFalse((prefix / "path_di_data.pkl").exists())


if __name__ == "__main__":
    unittest.main()
