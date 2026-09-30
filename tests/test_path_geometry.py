from __future__ import annotations

import ast
from collections import defaultdict
import copy
import itertools
import json
import logging
import os
from pathlib import Path
import pickle
import re
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import networkx as nx
import numpy as np

import breakend_graph as bg
import nclose_preprocess as pre
from nclose_tracking import discover_step11_indel_events
from path_geometry import conjoined_outer_geometry, indel_depth_vector, reference_connection_allowed, walked_strand


ROOT = Path(__file__).resolve().parents[1]


def script_helpers(filename):
    """Execute real stage helpers without launching the script's pipeline."""
    path = ROOT / filename
    tree = ast.parse(path.read_text())
    namespace = dict(
        copy=copy, re=re, os=os, nx=nx, np=np, logging=logging,
        defaultdict=defaultdict, pairwise=itertools.pairwise, groupby=itertools.groupby,
        reference_connection_allowed=reference_connection_allowed,
        conjoined_outer_geometry=conjoined_outer_geometry,
        walked_strand=walked_strand,
    )
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id.isupper():
                    try:
                        namespace[target.id] = ast.literal_eval(node.value)
                    except (ValueError, TypeError):
                        pass
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef)]
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)
    return namespace


def rows_and_paf(specs):
    rows, paf = [], []
    for index, (name, strand, start, end, kind) in enumerate(specs):
        members = [i for i, spec in enumerate(specs) if spec[0] == name]
        length = end - start
        query_start = sum(specs[i][3] - specs[i][2] for i in members if i < index)
        query_len = sum(specs[i][3] - specs[i][2] for i in members)
        base = [name, query_len, query_start, query_start + length, strand, "chr1", 200_000_000, start, end]
        rows.append(base + [60, kind, min(members), max(members), "0", "0", "0", "0", "0", "0", strand, "chr1", f"0.{index}"])
        paf.append(base + [length, length, 60, "tp:A:P", f"cs:Z::{length}"])
    return rows, paf


def bind_depth(specs):
    namespace = script_helpers("21_run_depth.py")
    rows, paf = rows_and_paf(specs)
    namespace.update(contig_data=rows, paf_file=[paf])
    return namespace


def depth_graph(namespace, endpoints):
    rows = namespace['contig_data']
    labels = namespace['node_label'](rows)
    adjacency = namespace['graph_build'](namespace['initial_graph_build'](rows), labels)
    adjacency = namespace['edge_optimization'](rows, adjacency)
    typed = defaultdict(list)
    for direction in range(2):
        for index in range(len(rows)):
            typed[direction, index].extend(adjacency[direction][index])
    adjacency = namespace['connect_nclose_telo'](
        rows, namespace['find_using_node'](rows, labels), typed, endpoints, [],
    )
    graph = nx.DiGraph()
    for index in endpoints:
        graph.add_nodes_from([(2, index), (3, index)])
    for source, edges in adjacency.items():
        for direction, index, distance in edges:
            graph.add_edge(source, (direction, index), weight=distance)
    namespace['G'] = graph
    return graph


def paf_spans(lines):
    return [(fields[5], int(fields[7]), int(fields[8])) for fields in (line.split('\t') for line in lines)]


def conjoin(rows, pairs):
    with patch.object(pre, 'estimate_type2_global_noise_sigma', return_value=0.0), patch.object(
        pre, 'depth_df_to_by_chrom', return_value={}
    ), patch.object(pre, 'breakpoint_is_depth_balanced', return_value=False):
        return pre.conjoined_type4(rows, {'chr1': pairs}, None, None, {})


class ReferencePathTests(unittest.TestCase):
    SPECS = [
        ('a', '+', 300_000, 400_000, 1), ('a', '-', 700_000, 800_000, 1),
        ('b', '-', 500_000, 600_000, 1), ('b', '+', 100_000, 200_000, 1),
        ('cover', '+', 0, 500_000, 3),
    ]

    def test_containing_bridge_uses_the_walked_strand_in_both_directions(self):
        for cover_strand in ('+', '-'):
            specs = self.SPECS[:-1] + [('cover', cover_strand, 0, 500_000, 3)]
            ns = bind_depth(specs)
            depth_graph(ns, [0, 1, 2, 3])
            for key in (((1, 3), (1, 0)), ((0, 0), (0, 3)), ((2, 3), (3, 0)), ((2, 0), (3, 3))):
                with self.subTest(strand=cover_strand, key=key):
                    raw = ns['bnd_raw_contig_list'](key)
                    output, skipped = ns['process_raw_contig_list'](raw)
                    intervals = sorted(paf_spans(output))
                    self.assertEqual(intervals[0][1], 100_000)
                    self.assertEqual(intervals[-1][2], 400_000)
                    self.assertEqual(sum(end - start for _, start, end in intervals), 300_000)
                    self.assertTrue(all(a[2] == b[1] for a, b in itertools.pairwise(intervals)))
                    self.assertEqual(skipped, 0)

    def test_overlap_with_opposite_walked_strands_is_rejected(self):
        ns = bind_depth(self.SPECS)
        with self.assertRaisesRegex(AssertionError, 'against the path strand'):
            ns['process_raw_contig_list']([(1, 3), (0, 4), (1, 0)])

    def test_no_valid_graph_route_falls_back_to_a_forward_reference_bridge(self):
        ns = bind_depth(self.SPECS)
        graph = nx.DiGraph()
        graph.add_weighted_edges_from([((2, 3), (0, 4), 0), ((0, 4), (3, 0), 0)])
        ns['G'] = graph
        self.assertEqual(ns['bnd_raw_contig_list'](((1, 3), (1, 0))), [(1, 3), (1, 0)])
        output, _ = ns['process_raw_contig_list']([(1, 3), (1, 0)])
        self.assertEqual(paf_spans(output), [('chr1', 100_000, 200_000), ('chr1', 200_000, 300_000), ('chr1', 300_000, 400_000)])

    def test_internal_unitig_inversion_is_preserved(self):
        ns = bind_depth(self.SPECS)
        output, _ = ns['process_raw_contig_list']([(1, 0), (1, 1)])
        self.assertEqual(len(output), 2)
        self.assertEqual([row.split('\t')[4] for row in output], ['+', '-'])


class TelomerePathTests(unittest.TestCase):
    def test_stage10_cli_path_keeps_both_terminal_anchors_in_stage21(self):
        for with_nclose in (False, True):
            with self.subTest(with_nclose=with_nclose), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                ns = bind_depth([('a', '+', 0, 10_000_000, 3), ('a', '+', 10_000_000, 20_000_000, 3)])
                rows = ns['contig_data']
                if with_nclose:
                    rows[1][5] = ns['paf_file'][0][1][5] = 'chr2'
                last_terminal = 'chr2b' if with_nclose else 'chr1b'
                terminals = {'chr1f': [(1, 0, 0)], last_terminal: [(0, 1, 0)]}
                bg.save_stage10_input(root / 'input.pkl', rows, {'a': [(0, 1)]} if with_nclose else {}, terminals)
                (root / 'ref.fai').write_text('chr1\t200000000\nchr2\t200000000\n')
                (root / 'limit.json').write_text(json.dumps({'limit_combinations': [1, 0]}))
                process = subprocess.run([
                    sys.executable, str(ROOT / '10_Graph_Find_Paths.py'), str(root / 'input.pkl'),
                    str(root / 'ref.fai'), str(root / 'output'), '-t', '1', '-d', '1',
                    '--limit-combinations', str(root / 'limit.json'),
                ], text=True, capture_output=True, timeout=60)
                self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
                with (root / 'output/path_data.pkl').open('rb') as handle:
                    paths = pickle.load(handle)
                name = f'chr1f_{last_terminal}'
                self.assertEqual(len(paths[name]), 1)
                ns.update(path_list_dict=paths, output_folder=tmp, G=nx.DiGraph())
                _, keys = ns['get_key_from_index_file'](f'{tmp}/{name}/1.index.txt')
                self.assertEqual([key for key in keys if key[0] == 2], [(2, 0), (2, 1)])
                output = []
                for index, key in enumerate(keys):
                    ns['create_final_depth_paf']((key, index))
                    output.extend((root / f'{index}.paf').read_text().splitlines())
                self.assertEqual(paf_spans(output), [('chr1', 0, 10_000_000), ('chr2' if with_nclose else 'chr1', 10_000_000, 20_000_000)])

    def test_self_nclose_terminal_has_a_reverse_exit(self):
        rows, _ = rows_and_paf([('u', '+', 0, 10_000_000, 1), ('u', '+', 10_000_000, 20_000_000, 1), ('back', '+', 90_000_000, 100_000_000, 3)])
        rows[1][5] = rows[2][5] = 'chr2'
        _, reverse = bg.chr_correlation_maker(rows)
        adjacency = bg.initialize_bnd_graph(rows, {'u': [(0, 1)]}, {'chr1f': [(1, 0, 0)], 'chr2b': [(0, 2, 0)]}, len(rows), reverse)
        graph = nx.DiGraph()
        for source, targets in adjacency.items():
            for target in targets:
                graph.add_edge(source, target if isinstance(target, str) else tuple(target))
        self.assertTrue(nx.has_path(graph, 'chr1f', 'chr2b'))
        self.assertTrue(nx.has_path(graph, 'chr2b', 'chr1f'))

    def test_reversed_reference_telomere_anchors_are_not_joined(self):
        rows, _ = rows_and_paf([('f', '+', 100, 200, 3), ('b', '+', 0, 50, 3)])
        _, reverse = bg.chr_correlation_maker(rows)
        adjacency = bg.initialize_bnd_graph(rows, {}, {'chr1f': [(1, 0, 0)], 'chr1b': [(0, 1, 0)]}, len(rows), reverse)
        self.assertNotIn([3, 1], adjacency[2, 0])
        self.assertNotIn([3, 0], adjacency[2, 1])


class EcdnaLengthTests(unittest.TestCase):
    SPECS = [('a', '+', 1_000_000, 1_010_000, 1), ('a', '-', 51_000_000, 51_010_000, 1), ('b', '-', 50_980_000, 50_990_000, 1), ('b', '+', 980_000, 990_000, 1)]

    def test_short_circle_survives_long_reference_jumps_and_matches_stage21(self):
        ns = bind_depth(self.SPECS)
        rows = ns['contig_data']
        circuits = pre.find_ecdna_circuits(rows, {'a': [(0, 1)], 'b': [(2, 3)]})
        self.assertEqual(circuits, [(0, 1, 2, 3)])
        self.assertEqual(pre.circuit_length_calculator(circuits[0], rows), 60_000)
        self.assertEqual(pre.circuit_length_calculator(circuits[0][::-1], rows), 60_000)
        depth_graph(ns, [0, 1, 2, 3])
        with tempfile.TemporaryDirectory() as tmp:
            ns['create_final_depth_paf_ecdna'](circuits, tmp)
            output = (Path(tmp) / '1.paf').read_text().splitlines()
        self.assertEqual(sum(end - start for _, start, end in paf_spans(output)), 60_000)
        with patch.object(pre, 'CIRCUIT_ECDNA_LENGTH_LIMIT', 60_000):
            self.assertEqual(pre.find_ecdna_circuits(rows, {'a': [(0, 1)], 'b': [(2, 3)]}), [])

    def test_internal_pieces_count_but_nclose_reference_jump_does_not(self):
        specs = [self.SPECS[0], ('a', '+', 80_000_000, 80_005_000, 1), *self.SPECS[1:]]
        rows, _ = rows_and_paf(specs)
        self.assertEqual(pre.circuit_length_calculator((0, 2, 3, 4), rows), 65_000)


class ConjoinedGeometryTests(unittest.TestCase):
    REVERSED = [('a', '+', 160_000, 170_000, 1), ('a', '-', 80_000, 90_000, 1), ('b', '-', 100_000, 110_000, 1), ('b', '+', 150_000, 160_000, 1)]
    OVERLAPPING_DUP = [('a', '+', 100_000, 120_000, 1), ('a', '-', 140_000, 150_000, 1), ('b', '-', 130_000, 140_000, 1), ('b', '+', 110_000, 130_000, 1)]

    def test_short_candidates_are_independent_of_owner_order(self):
        rows, _ = rows_and_paf(self.REVERSED)
        self.assertEqual(conjoin(rows, [(0, 1), (2, 3)]), ([], [(1, 0, 3, 2)]))
        self.assertEqual(conjoin(rows, [(2, 3), (0, 1)]), ([], [(1, 0, 3, 2)]))

    def test_overlapping_duplication_is_classified_by_breakends(self):
        rows, _ = rows_and_paf(self.OVERLAPPING_DUP)
        ins, dels = conjoin(rows, [(0, 1), (2, 3)])
        self.assertIn((0, 1, 2, 3), ins)
        self.assertEqual(dels, [])
        self.assertEqual(conjoin(rows, [(2, 3), (0, 1)]), (ins, dels))

    def write_compound(self, tmp, specs, circuit, event_type):
        rows, paf = rows_and_paf(specs)
        ns = script_helpers('11_Ref_Outlier_Contig_Modify.py')
        ns.update(contig_data=rows, paf_file=[paf], chr_data={'chr1': 200_000_000})
        root = Path(tmp) / '11_ref_ratio_outliers'
        (root / event_type).mkdir(parents=True)
        ns['write_conjoined_depth_pafs'](circuit, 1, 1, str(root))
        circuits = ([circuit], []) if event_type == 'back_jump' else ([], [circuit])
        with (Path(tmp) / 'conjoined_type4_ins_del.pkl').open('wb') as handle:
            pickle.dump(circuits, handle)
        return discover_step11_indel_events(tmp)[0]

    def test_catalog_uses_rc_traversal_without_rewriting_alignment_strands(self):
        with tempfile.TemporaryDirectory() as tmp:
            event = self.write_compound(tmp, self.REVERSED, (1, 0, 3, 2), 'front_jump')
            self.assertEqual((event['start_pos'], event['start_dir'], event['end_pos'], event['end_dir']), (90_000, '+', 100_000, '+'))
            self.assertEqual((event['st'], event['nd']), (90_000, 100_000))
            self.assertEqual([line.split('\t')[4] for line in Path(event['primary_paf']).read_text().splitlines()], ['-', '-'])

    def test_dup_label_and_span_change_without_changing_depth_delta(self):
        with tempfile.TemporaryDirectory() as tmp:
            event = self.write_compound(tmp, self.OVERLAPPING_DUP, (0, 1, 2, 3), 'back_jump')
            self.assertEqual((event['indel_kind'], event['st'], event['nd']), ('insertion', 110_000, 120_000))
            self.assertEqual(event['depth_base_sign'], -1)
            ns = bind_depth(self.OVERLAPPING_DUP)
            depth_graph(ns, [0, 1, 2, 3])
            ns['create_final_depth_paf_type2'](([(0, 1, 2, 3)], []), tmp)
            primary = Path(event['primary_paf']).read_text().splitlines()
            inner = (Path(tmp) / '11_ref_ratio_outliers/type2_ins/1.paf').read_text().splitlines()
            base = Path(event['base_paf']).read_text().splitlines()
            def coverage(lines):
                return np.array([sum(start <= pos < end for _, start, end in paf_spans(lines)) for pos in range(100_000, 150_000, 10_000)])
            delta = indel_depth_vector(coverage(primary + inner), coverage(base), event)
            np.testing.assert_array_equal(delta, [0, 1, 0, 1, 1])

    def test_nonoverlapping_dup_keeps_positive_baseline_and_ordinary_indels_keep_signs(self):
        first, last = rows_and_paf([('a', '+', 200, 300, 1), ('b', '+', 0, 100, 1)])[0]
        geometry = conjoined_outer_geometry(first, last, (0, 1, 2, 3))
        self.assertEqual((geometry['st'], geometry['nd']), (0, 300))
        self.assertEqual((geometry['depth_base_st'], geometry['depth_base_nd'], geometry['depth_base_sign']), (100, 200, 1))
        self.assertEqual(indel_depth_vector(2, 3, {'event_type': 'front_jump'}), -1)
        self.assertEqual(indel_depth_vector(2, 3, {'event_type': 'back_jump'}), 5)


if __name__ == '__main__':
    unittest.main()
