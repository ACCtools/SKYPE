from __future__ import annotations

import ast
import logging
import pickle
import tempfile
import unittest
from collections import Counter, defaultdict
from dataclasses import dataclass, fields
from pathlib import Path
from types import SimpleNamespace


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "02_Build_Breakend_Graph_Limited.py"
)


def load_graph_attempt_helpers():
    """Load graph-attempt definitions without running stage 02's CLI."""

    definition_names = {
        "PregraphState",
        "GraphAttemptState",
        "canonical_nclose_snapshot",
        "augment_nclose_nodes_with_type4_indels",
        "get_type4_indel_zero_dim_edge_set",
        "build_nonzero_telo_set",
        "build_graph_attempt_state",
    }
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    definitions = [
        node
        for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef))
        and node.name in definition_names
    ]
    module = ast.Module(body=definitions, type_ignores=[])
    namespace = {
        "dataclass": dataclass,
        "NCloseSourceConfig": object,
        "NCloseBuildResult": object,
        "Counter": Counter,
        "defaultdict": defaultdict,
        "logging": logging,
        "pkl": pickle,
        "CHR_STR": 0,
        "CHR_END": 1,
        "TYPE4_INDEL_GRAPH_MIN_SPAN": 1_000_000,
        "TYPE4_INDEL_GRAPH_EDGE_PKL": "type4_indel_graph_edges.pkl",
        "TELOMERE_INFO_FILE_PATH": "/unused/telomeres.bed",
        "VCF_TYPE4_GRAPH_NODE_PREFIX": "vcf_type4_graph_",
    }
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


HELPERS = load_graph_attempt_helpers()
GraphAttemptState = HELPERS["GraphAttemptState"]
build_graph_attempt_state = HELPERS["build_graph_attempt_state"]


def interval_distance(first, second):
    if first[1] < second[0]:
        return second[0] - first[1]
    if second[1] < first[0]:
        return first[0] - second[1]
    return 0


def make_pregraph(label, node_offset):
    contig_data = [
        (10_000 + idx * 10, 10_005 + idx * 10, f"unused-{idx}")
        for idx in range(node_offset + 4)
    ]
    contig_data[node_offset] = (0, 10, f"{label}-nearest")
    contig_data[node_offset + 1] = (20, 30, f"{label}-penalized")
    contig_data[node_offset + 2] = (40, 50, f"{label}-type4-a")
    contig_data[node_offset + 3] = (60, 70, f"{label}-type4-b")

    canonical_nodes = {
        f"canonical-{label}": [(node_offset, node_offset + 1)]
    }
    result = SimpleNamespace(
        contig_data=contig_data,
        chr_len={"chr1": 1_000_000},
        nclose_nodes=canonical_nodes,
        df={"type4_pair": (node_offset + 2, node_offset + 3)},
        repeat_censat_data={},
        telo_contig={f"chr1f-{label}": [("in", node_offset)]},
        contig_data_size=len(contig_data),
        chr_rev_corr={len(contig_data): "chr1f"},
    )
    return SimpleNamespace(
        label=label,
        connected_nodes=[
            ("chr1f", node_offset),
            ("chr1f", node_offset + 1),
        ],
        nclose=result,
    )


class GraphAttemptStateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.active_pregraph = None
        self.initializer_calls = []

        def activate_pregraph_state(pregraph):
            self.active_pregraph = pregraph

        def connected_nodes(_path):
            return list(self.active_pregraph.connected_nodes)

        def import_telo_data(_path, _chr_len):
            return [("chr1", 0, 10, "f")]

        def collect_type4(_contig_data, depth_data, _censat_data):
            pair = tuple(depth_data["type4_pair"])
            return [{"type4_tuple": pair, "indel_kind": "deletion"}]

        def select_type4(_contig_data, _nclose_nodes, candidates):
            pair = tuple(candidates[0]["type4_tuple"])
            return [
                {
                    "type4_tuple": pair,
                    "indel_kind": "deletion",
                    "contig_name": f"type4-{self.active_pregraph.label}",
                    "src": ("for", pair[0]),
                    "dst": ("for", pair[1]),
                }
            ]

        def initialize_bnd_graph(
            contig_data,
            graph_nodes,
            telo_contig,
            contig_data_size,
            chr_rev_corr,
        ):
            node_indices = {
                node_idx
                for pair_list in graph_nodes.values()
                for pair in pair_list
                for node_idx in pair
            }
            self.initializer_calls.append(
                (
                    contig_data,
                    graph_nodes,
                    telo_contig,
                    contig_data_size,
                    chr_rev_corr,
                )
            )
            return {f"adjacency-{self.active_pregraph.label}": node_indices}

        HELPERS.update(
            {
                "activate_pregraph_state": activate_pregraph_state,
                "extract_telomere_connect_contig_bytuple": connected_nodes,
                "import_telo_data": import_telo_data,
                "distance_checker_tuple": interval_distance,
                "collect_type4_indel_graph_candidates": collect_type4,
                "select_type4_indel_graph_edges": select_type4,
                "initialize_bnd_graph": initialize_bnd_graph,
                "args": SimpleNamespace(add_indel_graph=True, vcf_input=None),
            }
        )

    def test_state_has_explicit_attempt_contract(self) -> None:
        self.assertEqual(
            [field.name for field in fields(GraphAttemptState)],
            [
                "pregraph",
                "nonzero_telo_set",
                "selected_type4_indel_graph_edges",
                "type4_indel_zero_dim_edge_set",
                "graph_nclose_nodes",
                "bnd_graph_adjacency",
            ],
        )

    def test_successive_pregraphs_rebuild_all_node_indexed_graph_state(self) -> None:
        first_pregraph = make_pregraph("first", 0)
        second_pregraph = make_pregraph("second", 4)
        first_canonical_before = pickle.dumps(
            first_pregraph.nclose.nclose_nodes
        )
        second_canonical_before = pickle.dumps(
            second_pregraph.nclose.nclose_nodes
        )

        with tempfile.TemporaryDirectory() as temporary:
            HELPERS["PREFIX"] = temporary
            first = build_graph_attempt_state(first_pregraph)
            second = build_graph_attempt_state(second_pregraph)

            with (Path(temporary) / "type4_indel_graph_edges.pkl").open(
                "rb"
            ) as handle:
                persisted_edges = pickle.load(handle)

        self.assertIs(first.pregraph, first_pregraph)
        self.assertIs(second.pregraph, second_pregraph)
        self.assertEqual(first.nonzero_telo_set, {1})
        self.assertEqual(second.nonzero_telo_set, {5})
        self.assertTrue(
            first.nonzero_telo_set.isdisjoint(second.nonzero_telo_set)
        )

        first_zero_dim = first.type4_indel_zero_dim_edge_set
        second_zero_dim = second.type4_indel_zero_dim_edge_set
        self.assertEqual(first_zero_dim, {(('for', 2), ('for', 3))})
        self.assertEqual(second_zero_dim, {(('for', 6), ('for', 7))})
        self.assertTrue(first_zero_dim.isdisjoint(second_zero_dim))

        first_adjacency_nodes = first.bnd_graph_adjacency["adjacency-first"]
        second_adjacency_nodes = second.bnd_graph_adjacency["adjacency-second"]
        self.assertEqual(first_adjacency_nodes, {0, 1, 2, 3})
        self.assertEqual(second_adjacency_nodes, {4, 5, 6, 7})
        self.assertTrue(first_adjacency_nodes.isdisjoint(second_adjacency_nodes))

        self.assertEqual(
            pickle.dumps(first_pregraph.nclose.nclose_nodes),
            first_canonical_before,
        )
        self.assertEqual(
            pickle.dumps(second_pregraph.nclose.nclose_nodes),
            second_canonical_before,
        )
        self.assertNotIn("type4-first", first_pregraph.nclose.nclose_nodes)
        self.assertNotIn("type4-second", second_pregraph.nclose.nclose_nodes)
        self.assertIn("type4-first", first.graph_nclose_nodes)
        self.assertIn("type4-second", second.graph_nclose_nodes)

        self.assertEqual(len(self.initializer_calls), 2)
        self.assertIs(self.initializer_calls[0][0], first_pregraph.nclose.contig_data)
        self.assertIs(self.initializer_calls[1][0], second_pregraph.nclose.contig_data)
        self.assertEqual(
            persisted_edges,
            second.selected_type4_indel_graph_edges,
        )


class GraphAttemptOrchestrationTests(unittest.TestCase):
    def test_node_count_retry_finishes_before_graph_attempt_build(self) -> None:
        tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
        node_retry = next(
            node
            for node in tree.body
            if isinstance(node, ast.If)
            and any(
                isinstance(child, ast.Name)
                and child.id == "FAIL_NCLOSE_COUNT"
                for child in ast.walk(node.test)
            )
        )
        retry_build = next(
            node
            for node in ast.walk(node_retry)
            if isinstance(node, ast.Assign)
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == "pregraph_state"
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id == "rebuild_primary_only_pregraph_state"
        )
        graph_build = next(
            node
            for node in tree.body
            if isinstance(node, ast.Assign)
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == "graph_attempt_state"
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id == "build_and_activate_graph_attempt_state"
        )

        self.assertLess(retry_build.lineno, graph_build.lineno)
        self.assertLess(node_retry.end_lineno, graph_build.lineno)
        self.assertEqual(len(graph_build.value.args), 1)
        self.assertIsInstance(graph_build.value.args[0], ast.Name)
        self.assertEqual(graph_build.value.args[0].id, "pregraph_state")

    def test_divergence_retry_rebuilds_attempt_before_second_graph_run(self) -> None:
        tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
        fallback = next(
            node
            for node in tree.body
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.UnaryOp)
            and isinstance(node.test.op, ast.Not)
            and isinstance(node.test.operand, ast.Name)
            and node.test.operand.id == "last_success"
        )

        assignments = [
            node
            for node in ast.walk(fallback)
            if isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
        ]

        def assignment(target_name, call_name):
            return next(
                node
                for node in assignments
                if node.targets[0].id == target_name
                and node.value.func.id == call_name
            )

        rebuild = assignment(
            "pregraph_state", "rebuild_primary_only_pregraph_state"
        )
        graph_rebuild = assignment(
            "graph_attempt_state", "build_and_activate_graph_attempt_state"
        )
        second_run = assignment("last_success", "run_graph_pipeline")

        self.assertLess(rebuild.lineno, graph_rebuild.lineno)
        self.assertLess(graph_rebuild.lineno, second_run.lineno)
        self.assertEqual(len(graph_rebuild.value.args), 1)
        self.assertIsInstance(graph_rebuild.value.args[0], ast.Name)
        self.assertEqual(graph_rebuild.value.args[0].id, "pregraph_state")


if __name__ == "__main__":
    unittest.main()
