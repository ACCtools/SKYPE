from __future__ import annotations

import ast
import json
import pickle
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import breakend_graph as bg


SKYPE_ROOT = Path(__file__).resolve().parents[1]
STAGE10_PATH = SKYPE_ROOT / "10_Graph_Find_Paths.py"
STAGE21_PATH = SKYPE_ROOT / "21_run_depth.py"


def load_stage10_namespace():
    namespace = {
        "__name__": "stage10_test",
        "__file__": str(STAGE10_PATH),
    }
    tree = ast.parse(STAGE10_PATH.read_text(encoding="utf-8"))
    exec(compile(tree, str(STAGE10_PATH), "exec"), namespace)
    return namespace


def make_node(
    name,
    chrom,
    ref_start,
    ref_end,
    *,
    direction="+",
    contig_type=1,
    start_idx=0,
    end_idx=0,
):
    ref_len = ref_end - ref_start
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
        start_idx,
        end_idx,
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        direction,
        chrom,
        "0.0",
    )


class Stage10HandoffTests(unittest.TestCase):
    def test_named_handoff_round_trip_has_exact_public_values(self):
        contig_data = [make_node("a", "chr1", 0, 10)]
        nclose_nodes = {"utg001": [(0, 0)]}
        telo_contig = {"chr1f": [(bg.DIR_FOR, 0, 0)]}

        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "01_nclose_data.pkl"
            bg.save_stage10_input(
                path,
                contig_data,
                nclose_nodes,
                telo_contig,
            )
            with path.open("rb") as handle:
                raw = pickle.load(handle)
            observed = bg.load_stage10_input(path)

        self.assertEqual(
            set(raw),
            {"contig_data", "nclose_nodes", "telo_contig"},
        )
        self.assertEqual(observed.contig_data, contig_data)
        self.assertEqual(observed.nclose_nodes, nclose_nodes)
        self.assertEqual(observed.telo_contig, telo_contig)

    def test_handoff_rejects_hidden_extra_state(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "01_nclose_data.pkl"
            with path.open("wb") as handle:
                pickle.dump(
                    {
                        "contig_data": [],
                        "nclose_nodes": {},
                        "telo_contig": {},
                        "removed_compression_state": {},
                    },
                    handle,
                )
            with self.assertRaisesRegex(ValueError, "unexpected fields"):
                bg.load_stage10_input(path)

    def test_one_owner_can_branch_from_a_shared_endpoint(self):
        contig_data = [
            make_node("utg001", "chr1", 0, 10),
            make_node("utg001", "chr2", 0, 10),
            make_node("utg001", "chr3", 0, 10),
        ]
        nclose_nodes = {"utg001": [(0, 1), (0, 2)]}
        _, reverse_terminals = bg.chr_correlation_maker(contig_data)

        adjacency = bg.initialize_bnd_graph(
            contig_data,
            nclose_nodes,
            {},
            len(contig_data),
            reverse_terminals,
        )

        self.assertEqual(
            adjacency[(bg.DIR_FOR, 0)],
            [[bg.DIR_FOR, 1], [bg.DIR_FOR, 2]],
        )
        self.assertIn([bg.DIR_BAK, 0], adjacency[(bg.DIR_BAK, 1)])
        self.assertIn([bg.DIR_BAK, 0], adjacency[(bg.DIR_BAK, 2)])

    def test_telomere_adjacency_has_entry_exit_and_chromosome_bridge(self):
        contig_data = [
            make_node("front", "chr1", 0, 10),
            make_node("back", "chr1", 20, 30),
        ]
        telo_contig = {
            "chr1f": [(bg.DIR_FOR, 0, 0)],
            "chr1b": [(bg.DIR_BAK, 1, 0)],
        }
        _, reverse_terminals = bg.chr_correlation_maker(contig_data)
        adjacency = bg.initialize_bnd_graph(
            contig_data,
            {},
            telo_contig,
            len(contig_data),
            reverse_terminals,
        )

        self.assertEqual(
            adjacency["chr1f"][:2],
            [(bg.DIR_FOR, 0), (bg.DIR_OUT, 0)],
        )
        self.assertIn("chr1f", adjacency[(bg.DIR_IN, 0)])
        self.assertIn([bg.DIR_IN, 1], adjacency[(bg.DIR_OUT, 0)])
        self.assertIn([bg.DIR_IN, 0], adjacency[(bg.DIR_OUT, 1)])


class DimensionGraphTests(unittest.TestCase):
    def test_cross_chromosome_edge_consumes_chromosome_dimension(self):
        contig_data = [
            make_node("a", "chr1", 0, 10),
            make_node("b", "chr2", 0, 10),
        ]
        graph = bg.build_dimension_graph(
            {(bg.DIR_FOR, 0): [(bg.DIR_FOR, 1)]},
            contig_data,
            set(),
            1,
            1,
        )
        self.assertTrue(
            graph.has_edge(
                (bg.DIR_FOR, 0, 0, 0),
                (bg.DIR_FOR, 1, 1, 0),
            )
        )
        self.assertFalse(
            graph.has_edge(
                (bg.DIR_FOR, 0, 0, 0),
                (bg.DIR_FOR, 1, 0, 0),
            )
        )

    def test_same_contig_direction_change_consumes_direction_dimension(self):
        contig_data = [
            make_node("a", "chr1", 0, 10, direction="+"),
            make_node("a", "chr1", 20, 30, direction="-"),
        ]
        graph = bg.build_dimension_graph(
            {(bg.DIR_FOR, 0): [(bg.DIR_FOR, 1)]},
            contig_data,
            set(),
            1,
            1,
        )
        self.assertTrue(
            graph.has_edge(
                (bg.DIR_FOR, 0, 0, 0),
                (bg.DIR_FOR, 1, 0, 1),
            )
        )

    def test_type4_edge_consumes_no_dimension(self):
        contig_data = [
            make_node("type4", "chr1", 0, 10),
            make_node("type4", "chr2", 20, 30),
        ]
        source = (bg.DIR_FOR, 0)
        target = (bg.DIR_FOR, 1)
        graph = bg.build_dimension_graph(
            {source: [target]},
            contig_data,
            {(source, target)},
            1,
            1,
        )
        self.assertTrue(
            graph.has_edge(
                (bg.DIR_FOR, 0, 0, 0),
                (bg.DIR_FOR, 1, 0, 0),
            )
        )


class GraphBuildStateTests(unittest.TestCase):
    def test_depth_supported_type4_adds_both_graph_only_directions(self):
        contig_data = [
            make_node(
                "vcf_type4_graph_0",
                "chr1",
                1_000_000,
                2_000_000,
                contig_type=4,
                start_idx=0,
                end_idx=1,
            ),
            make_node(
                "vcf_type4_graph_0",
                "chr1",
                4_000_000,
                5_000_000,
                contig_type=4,
                start_idx=0,
                end_idx=1,
            ),
        ]
        depth_df = bg.pd.DataFrame(
            {
                "chr": ["chr1"] * 4,
                "st": [1_500_000, 2_000_000, 3_500_000, 4_000_000],
                "nd": [2_000_000, 2_500_000, 4_000_000, 4_500_000],
                "meandepth": [40.0, 10.0, 10.0, 40.0],
            }
        )
        stage_input = bg.Stage10Input(contig_data, {}, {})

        state = bg.build_graph_state(
            stage_input,
            add_indel_graph=True,
            depth_df=depth_df,
            censat_intervals={"chr1": [(2_500_000, 3_500_000)]},
            vcf_input=True,
        )

        self.assertEqual(stage_input.nclose_nodes, {})
        self.assertEqual(
            state.graph_nclose_nodes["vcf_type4_graph_0"],
            [(0, 1)],
        )
        self.assertEqual(len(state.selected_type4_indel_graph_edges), 2)
        self.assertEqual(
            state.type4_indel_zero_dim_edge_set,
            {
                ((bg.DIR_FOR, 0), (bg.DIR_FOR, 1)),
                ((bg.DIR_BAK, 1), (bg.DIR_BAK, 0)),
            },
        )

    def test_graph_only_type4_does_not_mutate_canonical_nclose(self):
        contig_data = [
            make_node("canonical", "chr1", 0, 10),
            make_node("canonical", "chr2", 20, 30),
            make_node("vcf_type4_graph_0", "chr1", 40, 50, contig_type=4),
            make_node("vcf_type4_graph_0", "chr1", 60, 70, contig_type=4),
        ]
        canonical = {"canonical": [(0, 1)]}
        stage_input = bg.Stage10Input(contig_data, canonical, {})
        selected = [
            {
                "type4_tuple": (2, 3),
                "indel_kind": "deletion",
                "contig_name": "vcf_type4_graph_0",
                "src": (bg.DIR_FOR, 2),
                "dst": (bg.DIR_FOR, 3),
            }
        ]

        with patch.object(
            bg,
            "collect_type4_indel_graph_candidates",
            return_value=[{"type4_tuple": (2, 3)}],
        ), patch.object(
            bg,
            "select_type4_indel_graph_edges",
            return_value=selected,
        ):
            state = bg.build_graph_state(
                stage_input,
                add_indel_graph=True,
                depth_df=object(),
                censat_intervals={},
                vcf_input=True,
            )

        self.assertEqual(canonical, {"canonical": [(0, 1)]})
        self.assertNotIn("vcf_type4_graph_0", canonical)
        self.assertEqual(
            state.graph_nclose_nodes["vcf_type4_graph_0"],
            [(2, 3)],
        )
        self.assertEqual(
            state.type4_indel_zero_dim_edge_set,
            {((bg.DIR_FOR, 2), (bg.DIR_FOR, 3))},
        )

    def test_successive_inputs_do_not_share_graph_state(self):
        first = bg.Stage10Input(
            [
                make_node("first", "chr1", 0, 10),
                make_node("first", "chr2", 20, 30),
            ],
            {"first": [(0, 1)]},
            {},
        )
        second = bg.Stage10Input(
            [
                make_node("second", "chr3", 0, 10),
                make_node("second", "chr4", 20, 30),
            ],
            {"second": [(0, 1)]},
            {},
        )
        first_state = bg.build_graph_state(first)
        second_state = bg.build_graph_state(second)

        self.assertIn("first", first_state.graph_nclose_nodes)
        self.assertNotIn("first", second_state.graph_nclose_nodes)
        self.assertIn("second", second_state.graph_nclose_nodes)


class Stage10CliTests(unittest.TestCase):
    def test_minimal_graph_writes_paths_without_path_score_artifact(self):
        # Two telomere anchors 5 Mb apart form one accepted chr1f->chr1b path.
        contig_data = [
            make_node("front", "chr1", 0, 1_000_000),
            make_node("back", "chr1", 4_000_000, 5_000_000),
        ]
        telo_contig = {
            "chr1f": [(bg.DIR_FOR, 0, 0)],
            "chr1b": [(bg.DIR_BAK, 1, 0)],
        }

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "01_nclose_data.pkl"
            prefix = root / "output"
            fai_path = root / "reference.fai"
            limit_path = root / "limits.json"
            bg.save_stage10_input(input_path, contig_data, {}, telo_contig)
            fai_path.write_text("chr1\t5000000\n", encoding="utf-8")
            limit_path.write_text(
                json.dumps({"limit_combinations": [1, 0]}),
                encoding="utf-8",
            )

            namespace = load_stage10_namespace()
            return_code = namespace["main"](
                [
                    str(input_path),
                    str(fai_path),
                    str(prefix),
                    "-t",
                    "1",
                    "-d",
                    "1",
                    "--limit-combinations",
                    str(limit_path),
                ]
            )

            with (prefix / "path_data.pkl").open("rb") as handle:
                path_data = pickle.load(handle)

            self.assertEqual(return_code, 0)
            self.assertFalse((prefix / "path_di_data.pkl").exists())
            self.assertIn("chr1f_chr1b", path_data)
            self.assertEqual(len(path_data["chr1f_chr1b"]), 1)

    def test_hard_nclose_failure_replaces_stale_type4_edges(self):
        contig_data = [make_node("a", "chr1", 0, 10)]
        too_many_pairs = [(0, 0)] * (bg.FAIL_NCLOSE_COUNT // 2 + 1)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            input_path = root / "01_nclose_data.pkl"
            prefix = root / "output"
            prefix.mkdir()
            edge_path = prefix / "type4_indel_graph_edges.pkl"
            bg.save_stage10_input(
                input_path,
                contig_data,
                {"owner": too_many_pairs},
                {},
            )
            with edge_path.open("wb") as handle:
                pickle.dump([{"stale": True}], handle)

            return_code = load_stage10_namespace()["main"](
                [
                    str(input_path),
                    str(root / "unused.fai"),
                    str(prefix),
                    "-t",
                    "1",
                ]
            )
            with edge_path.open("rb") as handle:
                edges = pickle.load(handle)

        self.assertEqual(return_code, 1)
        self.assertEqual(edges, [])


class GraphLimitFallbackTests(unittest.TestCase):
    @staticmethod
    def run_with_counts(count_by_limit, *, fixed_limit=None):
        attempts = []

        class FakePool:
            def __init__(self, *, processes, initializer, initargs):
                del processes
                self.context = initargs[0]
                attempts.append(
                    (
                        self.context.chr_change_limit,
                        self.context.dir_change_limit,
                    )
                )
                initializer(*initargs)

            def __enter__(self):
                return self

            def __exit__(self, exc_type, exc_value, traceback):
                return False

            def terminate(self):
                pass

            def imap_unordered(self, function, terminal_pairs):
                del function
                count = count_by_limit.get(attempts[-1], 0)
                for pair_index, (source, target) in enumerate(terminal_pairs):
                    pair = (
                        self.context.chr_rev_corr[source],
                        self.context.chr_rev_corr[target],
                    )
                    yield pair, count if pair_index == 0 else 0, []

        stage_input = bg.Stage10Input([], {}, {})
        graph_state = bg.GraphBuildState(
            graph_nclose_nodes={},
            selected_type4_indel_graph_edges=(),
            type4_indel_zero_dim_edge_set=frozenset(),
            bnd_graph_adjacency=defaultdict(list),
        )
        config = bg.PathSearchConfig(
            thread=1,
            graph_depth=2,
            fixed_limit_combination=fixed_limit,
            total_path_limit=10,
        )
        with patch.object(bg, "Pool", FakePool), patch.object(
            bg,
            "build_dimension_graph",
            return_value=bg.nx.DiGraph(),
        ):
            result = bg.run_path_search(stage_input, graph_state, {}, config)
        return result, attempts

    def test_automatic_limits_fall_back_to_last_success(self):
        result, attempts = self.run_with_counts(
            {
                (2, 1): 10,
                (1, 1): 1,
                (1, 0): 0,
            }
        )

        self.assertEqual(
            attempts,
            [(2, 1), (1, 1), (2, 1)],
        )
        self.assertIsNotNone(result)
        self.assertEqual(
            (result.chr_change_limit, result.dir_change_limit),
            (1, 1),
        )

    def test_all_automatic_limits_can_fail(self):
        result, attempts = self.run_with_counts(
            {
                (2, 1): 10,
                (1, 1): 10,
                (1, 0): 10,
            }
        )

        self.assertIsNone(result)
        self.assertEqual(attempts, [(2, 1), (1, 1), (1, 0)])

    def test_failed_fixed_limit_does_not_fall_back(self):
        result, attempts = self.run_with_counts(
            {(2, 1): 10},
            fixed_limit=(2, 1),
        )

        self.assertIsNone(result)
        self.assertEqual(attempts, [(2, 1)])


class Stage21PathScoreContractTests(unittest.TestCase):
    def test_stage21_no_longer_reads_or_sorts_by_path_di(self):
        source = STAGE21_PATH.read_text(encoding="utf-8")
        self.assertNotIn("path_di_data.pkl", source)
        self.assertNotIn("dep_sort_ind_list", source)
        self.assertIn("dep_list = [0] * len(paf_ans_list)", source)

    def test_stage10_core_has_no_removed_cluster_or_score_state(self):
        source = (SKYPE_ROOT / "breakend_graph.py").read_text(encoding="utf-8")
        self.assertNotIn("st_compress", source)
        self.assertNotIn("ed_compress", source)
        self.assertNotIn("nonzero_telo_set", source)
        self.assertNotIn("path_di_list", source)
        self.assertNotIn("path_di_data.pkl", source)


if __name__ == "__main__":
    unittest.main()
