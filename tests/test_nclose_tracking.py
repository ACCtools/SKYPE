import csv
import os
import pickle
import tempfile
import unittest
from collections import Counter, OrderedDict

import numpy as np

from nclose_tracking import (
    REPORT_COLUMNS,
    append_ecdna_events_to_catalog,
    bed_visible_ecdna_indices,
    bed_visible_ecdna_indices_across_stages,
    build_bnd_event_catalog,
    build_nclose_report_rows,
    calculate_event_weights,
    compressed_bnd_event_keys,
    count_ecdna_circuit_events,
    count_index_path_events,
    discover_step11_indel_events,
    ecdna_circuit_event_keys,
    format_nclose_ids,
    initialise_filter_status,
    load_type4_edge_event_map,
    nclose_event_id_by_key,
    record_filter_stage,
    save_event_catalog,
    write_nclose_report,
)


def paf_row(name, strand, chrom, chrom_start, chrom_end):
    return [
        name, "1000", "0", "1000", strand, chrom, "250000000",
        str(chrom_start), str(chrom_end),
    ]


def write_paf(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wt") as handle:
        for row in rows:
            handle.write("\t".join(map(str, row)) + "\n")


def make_event(key, pos):
    return {
        "event_key": key,
        "kind": "bnd",
        "start_chr": "chr1",
        "start_pos": pos,
        "start_dir": "+",
        "end_chr": "chr2",
        "end_pos": pos + 1,
        "end_dir": "-",
    }


class NCloseCatalogTests(unittest.TestCase):
    def test_bnd_catalog_matches_compressed_txt_order_and_direction(self):
        contig_data = [
            paf_row("utg050585l", "-", "chr1", 189967205, 189986475),
            paf_row("utg050585l", "+", "chr13", 94503081, 94547399),
            paf_row("second", "+", "chr2", 20, 30),
            paf_row("second", "-", "chr2", 80, 90),
        ]
        nclose_type = OrderedDict([
            (("chr1", "chr13"), [(0, 1)]),
            (("chr2", "chr2"), [(3, 2)]),
        ])

        catalog = build_bnd_event_catalog(nclose_type, contig_data)

        self.assertEqual([event["event_key"] for event in catalog], [(0, 1), (2, 3)])
        self.assertEqual(
            (
                catalog[0]["start_chr"], catalog[0]["start_pos"], catalog[0]["start_dir"],
                catalog[0]["end_chr"], catalog[0]["end_pos"], catalog[0]["end_dir"],
            ),
            ("chr1", 189967205, "-", "chr13", 94503081, "+"),
        )
        # The writer prints pair (3, 2); both directions therefore flip while
        # the event key remains canonical for path matching.
        self.assertEqual(
            (catalog[1]["start_pos"], catalog[1]["start_dir"],
             catalog[1]["end_pos"], catalog[1]["end_dir"]),
            (80, "+", 20, "-"),
        )

    def test_native_and_empty_vcf_indel_breakends(self):
        with tempfile.TemporaryDirectory() as prefix:
            front = os.path.join(prefix, "11_ref_ratio_outliers", "front_jump")
            back = os.path.join(prefix, "11_ref_ratio_outliers", "back_jump")

            native_rows = [
                paf_row("native", "+", "chr4", 100, 200),
                paf_row("native", "-", "chr4", 400, 500),
            ]
            write_paf(os.path.join(front, "1.paf"), native_rows)
            write_paf(os.path.join(front, "1_base.paf"), native_rows)

            write_paf(os.path.join(front, "2.paf"), [])
            write_paf(
                os.path.join(front, "2_base.paf"),
                [paf_row("empty-del", "+", "chr5", 300, 700)],
            )
            write_paf(os.path.join(back, "1.paf"), [])
            write_paf(
                os.path.join(back, "1_base.paf"),
                [paf_row("empty-dup", "+", "chr6", 900, 1200)],
            )
            reverse_native_rows = [
                paf_row("reverse-native", "-", "chr7", 1000, 1100),
                paf_row("reverse-native", "+", "chr7", 1500, 1600),
            ]
            write_paf(os.path.join(back, "2.paf"), reverse_native_rows)
            write_paf(
                os.path.join(back, "2_base.paf"),
                [paf_row("reverse-base", "+", "chr7", 1600, 1000)],
            )

            metadata = {
                ("front_jump", 2): {"chrom": "chr5", "st": 300, "nd": 700},
                ("back_jump", 1): {"chrom": "chr6", "st": 900, "nd": 1200},
            }
            with open(os.path.join(prefix, "vcf_type4_outlier_index.pkl"), "wb") as handle:
                pickle.dump(metadata, handle)

            events = discover_step11_indel_events(prefix)

        self.assertEqual(
            [event["event_key"] for event in events],
            [
                ("front_jump", 1, -1),
                ("front_jump", 2, -1),
                ("back_jump", 1, -1),
                ("back_jump", 2, -1),
            ],
        )
        self.assertEqual(
            (events[0]["start_chr"], events[0]["start_pos"], events[0]["start_dir"],
             events[0]["end_chr"], events[0]["end_pos"], events[0]["end_dir"]),
            ("chr4", 200, "+", "chr4", 500, "-"),
        )
        self.assertEqual(
            (events[1]["start_pos"], events[1]["start_dir"],
             events[1]["end_pos"], events[1]["end_dir"]),
            (300, "+", 700, "+"),
        )
        self.assertEqual(
            (events[2]["start_pos"], events[2]["start_dir"],
             events[2]["end_pos"], events[2]["end_dir"]),
            (1200, "+", 900, "+"),
        )
        self.assertEqual(
            (events[3]["start_pos"], events[3]["start_dir"],
             events[3]["end_pos"], events[3]["end_dir"],
             events[3]["st"], events[3]["nd"]),
            (1000, "-", 1500, "+", 1000, 1600),
        )
        self.assertEqual([event["report_index"] for event in events], [0, 1, 2, 3])

    def test_ecdna_inversions_get_stable_report_ids(self):
        contig_data = [
            paf_row(f"ctg{i // 2}", "+" if i % 2 == 0 else "-", "chr7", i * 10, i * 10 + 5)
            for i in range(6)
        ]
        compressed = build_bnd_event_catalog(
            OrderedDict([(("chr7", "chr7"), [(0, 1)])]),
            contig_data,
        )
        circuits = [(0, 1, 2, 3), (2, 3, 4, 5)]
        with tempfile.TemporaryDirectory() as prefix:
            save_event_catalog(prefix, compressed)
            catalog = append_ecdna_events_to_catalog(
                prefix, circuits, contig_data
            )

        self.assertEqual(
            [event["event_key"] for event in catalog],
            [(0, 1), (2, 3), (4, 5)],
        )
        self.assertEqual(compressed_bnd_event_keys(catalog), {(0, 1)})
        event_ids = nclose_event_id_by_key(catalog)
        self.assertEqual(event_ids[(0, 1)], "SKYPE.nclose.1")
        self.assertEqual(
            format_nclose_ids(
                ecdna_circuit_event_keys(circuits[1]), event_ids
            ),
            "SKYPE.nclose.2;SKYPE.nclose.3",
        )


class NCloseUsageTests(unittest.TestCase):
    def test_bed_visible_ecdna_gate_is_strict(self):
        self.assertEqual(
            bed_visible_ecdna_indices(
                np.array([1.0, 1.01, 2.0]),
                {7: 0, 8: 1, 9: 2},
                raw_weight_cutoff=1.0,
            ),
            [8, 9],
        )

    def test_bed_visible_ecdna_indices_include_all_bed_stages(self):
        self.assertEqual(
            bed_visible_ecdna_indices_across_stages(
                {
                    "base": np.array([1.01, 0.0, 0.0, 0.0]),
                    "filter": np.array([0.0, 1.01, 0.0, 0.0]),
                    "cluster": np.array([0.0, 0.0, 1.01, 1.0]),
                },
                {7: 0, 8: 1, 9: 2, 10: 3},
                raw_weight_cutoff=1.0,
            ),
            [7, 8, 9],
        )

    def test_all_path_types_and_weight_multiplicity(self):
        bnd_key = (5, 6)
        amplicon_bnd_key = (1, 2)
        second_amplicon_bnd_key = (3, 4)
        indel_key = ("front_jump", 1, -1)
        path = [
            (0, "start"), (0, 5), (0, 6), (0, 5), (0, 6),
            (0, "end"), (0, "terminal"),
        ]
        usage = count_index_path_events(
            path,
            {bnd_key},
            {((0, 6), (0, 5)): indel_key},
        )
        self.assertEqual(usage, Counter({bnd_key: 2, indel_key: 1}))

        # Each Amplicon consists of the two inversion NCloses (s1,e1) and
        # (s2,e2), and its matrix-column weight contributes once to each.
        amplicon_usage = count_ecdna_circuit_events(
            [1, 2, 3, 4],
            {bnd_key, amplicon_bnd_key, second_amplicon_bnd_key},
        )
        self.assertEqual(
            amplicon_usage,
            Counter({second_amplicon_bnd_key: 1, amplicon_bnd_key: 1}),
        )

        path_usage = [usage, Counter({indel_key: 1}), amplicon_usage, Counter()]
        raw = calculate_event_weights(path_usage, np.array([3.0, 5.0, 7.0, 100.0]))
        amplicon_component = calculate_event_weights(
            path_usage,
            np.array([3.0, 5.0, 7.0, 100.0]),
            active_columns={2},
        )
        self.assertEqual(raw[bnd_key], 6.0)
        self.assertEqual(raw[amplicon_bnd_key], 7.0)
        self.assertEqual(raw[second_amplicon_bnd_key], 7.0)
        self.assertEqual(raw[indel_key], 8.0)
        self.assertEqual(raw[bnd_key] / 2.0, 3.0)
        self.assertNotIn(bnd_key, amplicon_component)
        self.assertEqual(amplicon_component[amplicon_bnd_key], 7.0)
        self.assertEqual(amplicon_component[second_amplicon_bnd_key], 7.0)

    def test_type4_graph_edge_maps_to_emitted_indel(self):
        indel_key = ("front_jump", 1, -1)
        catalog = [{
            "event_key": indel_key,
            "kind": "indel",
            "event_type": "front_jump",
            "start_chr": "chr7", "start_pos": 100, "start_dir": "+",
            "end_chr": "chr7", "end_pos": 200, "end_dir": "+",
            "chrom": "chr7", "st": 100, "nd": 200,
        }]
        edge = {
            "src": (0, 11), "dst": (0, 12),
            "indel_kind": "deletion",
            "base_chrom": "chr7", "base_st": 102, "base_nd": 198,
        }
        with tempfile.TemporaryDirectory() as prefix:
            with open(os.path.join(prefix, "type4_indel_graph_edges.pkl"), "wb") as handle:
                pickle.dump([edge], handle)
            mapping = load_type4_edge_event_map(prefix, catalog)
        self.assertEqual(mapping, {((0, 11), (0, 12)): indel_key})


class NCloseStatusAndReportTests(unittest.TestCase):
    def setUp(self):
        self.keys = [(i, i + 100) for i in range(7)]
        self.catalog = [make_event(key, 1000 + i * 10) for i, key in enumerate(self.keys)]
        self.usage = [Counter({key: 1}) for key in self.keys]

    def build_status(self):
        status = initialise_filter_status(
            self.catalog, self.usage, range(len(self.usage))
        )
        record_filter_stage(
            status,
            "base",
            "initial",
            self.usage,
            {0, 3, 4, 5, 6},
            {self.keys[1]: "FILTERED_23_RECIPROCAL_PATH"},
            "FILTERED_23_COFILTERED_PATH",
        )
        record_filter_stage(
            status,
            "filter",
            "base",
            self.usage,
            {0, 4, 5, 6},
            {self.keys[3]: "FILTERED_24_GREEDY_DEPTH"},
            "FILTERED_24_COFILTERED_PATH",
        )
        record_filter_stage(
            status,
            "cluster",
            "filter",
            self.usage,
            {0, 6},
            {
                self.keys[4]: "FILTERED_24_CLUSTER_SELECTION",
                self.keys[5]: "FILTERED_24_CLUSTER_MIN_WEIGHT",
            },
            "FILTERED_24_CLUSTER_COFILTERED_PATH",
        )
        return status

    def test_stage02_explicit_no_path_and_cofiltered_reasons(self):
        catalog = self.catalog[:3]
        usage = [Counter({self.keys[0]: 1, self.keys[2]: 1})]
        status = initialise_filter_status(
            catalog,
            usage,
            active_columns=set(),
            explicit_reasons={self.keys[0]: "FILTERED_02_EXCLUDE_LIST"},
        )
        reasons = status["stages"]["initial"]["reasons"]
        self.assertEqual(reasons[self.keys[0]], "FILTERED_02_EXCLUDE_LIST")
        self.assertEqual(reasons[self.keys[1]], "FILTERED_02_NO_ELIGIBLE_PATH")
        self.assertEqual(reasons[self.keys[2]], "FILTERED_02_COFILTERED_PATH")

    def test_first_filter_reason_propagation_and_karyotype_values(self):
        status = self.build_status()
        weights = {
            "base": np.array([2, 2, 2, 2, 2, 2, 0], dtype=float),
            "filter": np.array([2, 0, 0, 0, 2, 0.1, 0], dtype=float),
            "cluster": np.array([2, 0, 0, 0, 0, 0, 0], dtype=float),
        }
        rows = build_nclose_report_rows(
            self.catalog, self.usage, status, 2.0, weights, True
        )

        self.assertEqual(rows[0]["nclose_cn"], "1")
        self.assertEqual(rows[0]["nclose_cluster_reason"], "PASS")
        self.assertEqual(
            (rows[1]["nclose_cn"], rows[1]["nclose_cn_reason"]),
            ("-1", "FILTERED_23_RECIPROCAL_PATH"),
        )
        self.assertEqual(rows[1]["nclose_filter_reason"], "FILTERED_23_RECIPROCAL_PATH")
        self.assertEqual(rows[2]["nclose_cn_reason"], "FILTERED_23_COFILTERED_PATH")
        self.assertEqual(rows[3]["nclose_filter_reason"], "FILTERED_24_GREEDY_DEPTH")
        self.assertEqual(rows[4]["nclose_cluster_reason"], "FILTERED_24_CLUSTER_SELECTION")
        self.assertEqual(rows[5]["nclose_cluster_reason"], "FILTERED_24_CLUSTER_MIN_WEIGHT")
        self.assertEqual(
            (rows[6]["nclose_cluster"], rows[6]["nclose_cluster_reason"]),
            ("0", "PASS"),
        )

    def test_variant_na_and_fixed_tsv_schema(self):
        status = self.build_status()
        weights = {"base": np.ones(len(self.usage), dtype=float)}
        rows = build_nclose_report_rows(
            self.catalog, self.usage, status, 2.0, weights, False
        )
        self.assertEqual(rows[0]["nclose_filter"], "NA")
        self.assertEqual(rows[0]["nclose_filter_reason"], "NOT_RUN_VARIANT_MODE")
        self.assertEqual(rows[0]["nclose_cluster"], "NA")
        self.assertEqual(rows[0]["nclose_id"], "SKYPE.nclose.1")

        with tempfile.TemporaryDirectory() as prefix:
            output = write_nclose_report(
                prefix, self.catalog, self.usage, status, 2.0, weights, False
            )
            with open(output, "rt") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                output_rows = list(reader)
        self.assertEqual(tuple(reader.fieldnames), REPORT_COLUMNS)
        self.assertEqual(reader.fieldnames[0], "nclose_id")
        self.assertEqual(len(output_rows), len(self.catalog))

    def test_cluster_min_reason_overrides_an_exact_zero_carrier(self):
        # One event can occur in both a positive-small column and an exact-zero
        # column.  The zero carrier must not hide the fact that a positive
        # cluster coefficient was explicitly cut by the 0.1N rule.
        key = self.keys[0]
        usage = [Counter({key: 1}), Counter({key: 1})]
        catalog = [self.catalog[0]]
        status = initialise_filter_status(catalog, usage, {0, 1})
        record_filter_stage(
            status, "base", "initial", usage, {0, 1}, {}, "BASE_COFILTERED"
        )
        record_filter_stage(
            status, "filter", "base", usage, {0, 1}, {}, "FILTER_COFILTERED"
        )
        record_filter_stage(
            status,
            "cluster",
            "filter",
            usage,
            {1},
            {key: "FILTERED_24_CLUSTER_MIN_WEIGHT"},
            "CLUSTER_COFILTERED",
            forced_reasons={key: "FILTERED_24_CLUSTER_MIN_WEIGHT"},
        )
        rows = build_nclose_report_rows(
            catalog,
            usage,
            status,
            2.0,
            {
                "base": np.array([1.0, 0.0]),
                "filter": np.array([1.0, 0.0]),
                "cluster": np.array([0.0, 0.0]),
            },
            True,
        )
        self.assertEqual(rows[0]["nclose_cluster"], "-1")
        self.assertEqual(
            rows[0]["nclose_cluster_reason"],
            "FILTERED_24_CLUSTER_MIN_WEIGHT",
        )


if __name__ == "__main__":
    unittest.main()
