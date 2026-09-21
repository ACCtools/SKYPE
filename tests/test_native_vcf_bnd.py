import collections
import csv
from pathlib import Path
import pickle
import tempfile
import unittest

import vcfpy

from native_bnd import node_pair_endpoints
from native_structure_output import display_events, write_native_vcf
from native_vcf_bnd import bnd_calls, nclose_bnd_memberships, source_junctions
from nclose_tracking import _build_bnd_event
from structure_nclose import MODEL_VERSION, StructureWeights, add_virtual_structure
from test_native_bnd import export_namespace, node


def chain(name, *alignments):
    return [node(name, i * 100, (i + 1) * 100, strand, chrom, start, end)
            for i, (chrom, start, end, strand) in enumerate(alignments)]


ABC = (("chr1", 100, 200, "+"), ("chr2", 300, 400, "+"), ("chr3", 500, 600, "+"))


def make_context(chains, structures, weights):
    nodes, pairs, events = [], [], {}
    for index, alignments in enumerate(chains, 1):
        pair = (len(nodes), len(nodes) + len(alignments) - 1)
        nodes.extend(alignments)
        pairs.append(pair)
        event = _build_bnd_event(pair, nodes, "compressed_nclose_nodes_list.txt")
        endpoints, names = node_pair_endpoints(nodes, pair)
        event.update(nclose_id=f"SKYPE.nclose.{index}", endpoints=endpoints, contig_names=names)
        events[pair] = event
    records = [dict(structure_id=f"S{i}", feature_index=i, kind=kind,
                    nclose_counts={pairs[index]: count for index, count in uses.items()},
                    telomere_counts={}, source=f"feature:{i}", legacy_ids=[])
               for i, (kind, uses) in enumerate(structures)]
    model = dict(version=MODEL_VERSION, column_count=len(records), next_nclose_id=len(events) + 1,
                 ncloses=events, structures=records)
    return StructureWeights(model, weights, 10.), nodes, pairs


class NativeVcfProjectionTests(unittest.TestCase):
    def test_vcf_split_leaves_model_nclose_totals_and_display_unchanged(self):
        context, nodes, (pair,) = make_context([chain("unitig", *ABC)], [("PATH", {0: 1})], [20.])
        display = display_events(context, nodes, {})
        before = pickle.dumps((context.__dict__, nodes))
        calls = bnd_calls(context, nodes)
        self.assertEqual([call["endpoints"] for call in calls], [
            (("chr1", 200, "L"), ("chr2", 300, "R")),
            (("chr2", 400, "L"), ("chr3", 500, "R")),
        ])
        self.assertEqual([call["weight"] for call in calls], [2., 2.])
        self.assertEqual(context.totals[pair], 20.)
        self.assertEqual(display_events(context, nodes, {}), display)
        self.assertEqual(pickle.dumps((context.__dict__, nodes)), before)

    def test_censat_and_telomere_children_are_split(self):
        for name in ("utg_censat", "subtelomere_cut_contig_1", "telomere_middle_cut_contig_2"):
            with self.subTest(name=name):
                context, nodes, _ = make_context([chain(name, *ABC)], [("PATH", {0: 1})], [10.])
                self.assertEqual(len(bnd_calls(context, nodes)), 2)
                self.assertEqual([j["mode"] for js in nclose_bnd_memberships(context, nodes).values() for j in js],
                                 ["ADJACENT_ALIGNMENT", "ADJACENT_ALIGNMENT"])

    def test_read_and_olc_rescue_split_including_restored_sources(self):
        for method in ("read", "olc"):
            for named in (True, False):
                with self.subTest(method=method, named=named):
                    name = f"raw_rescue_{method}_R0001" if named else "legacy_rescue"
                    context, nodes, (pair,) = make_context([chain(name, *ABC)], [("AMP", {0: 1})], [10.])
                    if not named:
                        context.model["nclose_sources"][pair]["rescue_method"] = method
                    before = pickle.dumps((context.__dict__, nodes))
                    calls = bnd_calls(context, nodes)
                    self.assertEqual([call["endpoints"] for call in calls], [
                        (("chr1", 200, "L"), ("chr2", 300, "R")),
                        (("chr2", 400, "L"), ("chr3", 500, "R")),
                    ])
                    self.assertEqual([call["weight"] for call in calls], [1., 1.])
                    self.assertEqual([j["mode"] for j in nclose_bnd_memberships(context, nodes)[pair]],
                                     ["ADJACENT_ALIGNMENT", "ADJACENT_ALIGNMENT"])
                    self.assertEqual(context.totals[pair], 10.)
                    self.assertEqual(pickle.dumps((context.__dict__, nodes)), before)

    def test_raw_and_unitig_exact_junctions_merge_source_weights(self):
        for method in ("read", "olc"):
            with self.subTest(method=method):
                context, nodes, (unitig, raw) = make_context(
                    [chain("unitig", *ABC), chain(f"raw_rescue_{method}_R1", *ABC)],
                    [("PATH", {0: 1}), ("PATH", {1: 1})], [20., 50.])
                self.assertEqual(context.model["nclose_aliases"][raw], unitig)
                self.assertEqual(context.totals[unitig], 70.)
                calls = bnd_calls(context, nodes)
                self.assertEqual(len(calls), 2)
                for call in calls:
                    self.assertEqual(call["keys"], {unitig, raw})
                    self.assertEqual(call["weight"], 7.)
                    self.assertEqual({row["source_nclose_key"]: row["contribution_N"]
                                      for row in call["rows"]}, {unitig: 2., raw: 5.})

    def test_raw_and_unitig_exact_outer_aliases_keep_separate_geometry_and_weight(self):
        context, nodes, (unitig, raw) = make_context(
            [chain("unitig", *ABC), chain("raw_rescue_olc_R1", ABC[0], ABC[-1])],
            [("PATH", {0: 1}), ("PATH", {1: 1})], [20., 50.])
        self.assertEqual(context.model["nclose_aliases"][raw], unitig)
        self.assertEqual(context.totals[unitig], 70.)
        calls = bnd_calls(context, nodes)
        self.assertEqual(len(calls), 3)
        outer, = [call for call in calls if call["keys"] == {raw}]
        self.assertEqual(outer["weight"], 5.)
        self.assertEqual([call["weight"] for call in calls if call["keys"] == {unitig}], [2., 2.])

    def test_same_outer_endpoints_with_different_interiors_do_not_share_contributions(self):
        other = (ABC[0], ("chr4", 700, 800, "+"), ABC[2])
        context, nodes, (a, b) = make_context([chain("a", *ABC), chain("b", *other)],
                                            [("PATH", {0: 1}), ("PATH", {1: 1})], [30., 70.])
        self.assertEqual(context.model["nclose_aliases"][b], a)
        calls = bnd_calls(context, nodes)
        self.assertEqual(len(calls), 4)
        self.assertEqual(sorted(call["weight"] for call in calls), [3., 3., 7., 7.])
        self.assertTrue(all(len(call["keys"]) == 1 for call in calls))

    def test_shared_internal_bnd_merges_only_exact_both_coordinates_and_sides(self):
        other = (ABC[0], ABC[1], ("chr4", 700, 800, "+"))
        chains = [chain("a", *ABC), chain("b", *other),
                  chain("near", ("chr1", 100, 201, "+"), ABC[1]),
                  chain("different_side", ("chr1", 200, 300, "-"), ABC[1]),
                  chain("other_mate", ABC[0], ("chr2", 301, 401, "+"))]
        context, nodes, pairs = make_context(chains, [("PATH", {i: 1}) for i in range(5)], [10.] * 5)
        calls = bnd_calls(context, nodes)
        shared, = [call for call in calls if len(call["keys"]) == 2]
        self.assertEqual(shared["keys"], set(pairs[:2]))
        self.assertEqual(shared["weight"], 2.)
        self.assertEqual(len(calls), 6)
        self.assertEqual(len({call["endpoints"] for call in calls}), 6)

    def test_reverse_complement_chain_merges_exact_junctions(self):
        reverse = tuple((chrom, start, end, "-") for chrom, start, end, _ in reversed(ABC))
        context, nodes, pairs = make_context([chain("forward", *ABC), chain("reverse", *reverse)],
                                            [("PATH", {0: 1}), ("PATH", {1: 1})], [10., 20.])
        calls = bnd_calls(context, nodes)
        self.assertEqual(len(calls), 2)
        self.assertTrue(all(call["weight"] == 3. and call["keys"] == set(pairs) for call in calls))

    def test_compounds_use_constituent_chains_without_cross_unitig_bridges(self):
        other = (("chr4", 100, 200, "-"), ("chr5", 300, 400, "+"), ("chr6", 500, 600, "+"))
        context, nodes, _ = make_context([chain("a", *ABC), chain("b", *other)],
                                        [("AMP", {0: 1, 1: 1}), ("MERGE_TYPE4", {0: 1, 1: 1})], [10., 20.])
        calls = bnd_calls(context, nodes)
        self.assertEqual(len(calls), 4)
        self.assertTrue(all(call["weight"] == 3. for call in calls))
        self.assertTrue(all(call["info"]["SVCLASS"] == ["AMPLICON", "MERGED_TYPE4"] for call in calls))
        self.assertFalse(any({ep[0] for ep in call["endpoints"]} == {"chr3", "chr4"} for call in calls))

    def test_virtual_inversion_restored_source_is_split_with_its_own_weight(self):
        context, nodes, (original,) = make_context([chain("a", *ABC)], [("PATH", {0: 1})], [30.])
        nodes.extend(chain("restored", ("chr4", 100, 200, "+"), ("chr5", 300, 400, "+"), ("chr6", 500, 600, "+")))
        record = dict(pair_id=7, nclose_key_a=original, nclose_key_b=(3, 5), layout_a={},
                      layout_b=dict(ordered_endpoints=(
                          dict(chrom="chr4", coord=200, dir="+", ctg_name="restored"),
                          dict(chrom="chr6", coord=500, dir="+", ctg_name="restored"))))
        add_virtual_structure(context.model, record, 20., [], nodes)
        context = StructureWeights(context.model, [30.], 10.)
        calls = bnd_calls(context, nodes)
        self.assertEqual(len(calls), 4)
        for call in calls:
            self.assertEqual(call["info"]["VIRTUAL_WEIGHT"], 2.)
            self.assertEqual(call["info"]["MODEL_WEIGHT"], 3. if original in call["keys"] else 0.)

    def test_repeated_junctions_multiply_source_nclose_occurrences(self):
        repeated = (ABC[0], ABC[1], ABC[0], ABC[1], ABC[2])
        context, nodes, _ = make_context([chain("repeat", *repeated)], [("PATH", {0: 3})], [20.])
        calls = bnd_calls(context, nodes)
        repeated_call, = [call for call in calls if call["weight"] == 12.]
        self.assertEqual(repeated_call["info"]["PARENT_MULTIPLICITY"], [6])
        row, = repeated_call["rows"]
        self.assertEqual((row["nclose_occurrence_count"], row["bnd_occurrences_per_nclose"], row["junction_indices"]),
                         (3, 2, (1, 3)))
        self.assertEqual(sorted(call["weight"] for call in calls), [6., 6., 12.])

    def test_only_selected_nclose_interval_is_split(self):
        nodes = chain("a", ("chr7", 10, 110, "+"), *ABC, ("chr8", 10, 110, "+"))
        event = _build_bnd_event((1, 3), nodes, "AMP")
        event["endpoints"] = node_pair_endpoints(nodes, (1, 3))[0]
        self.assertEqual([j["node_pair"] for j in source_junctions(event, nodes)], [(1, 2), (2, 3)])

    def test_normal_continuations_are_skipped_but_reference_discontinuities_remain(self):
        for strand in ("+", "-"):
            with self.subTest(strand=strand):
                spans = [(100, 200), (230, 330), (400, 500)]
                if strand == "-":
                    spans = [(600 - end, 600 - start) for start, end in spans]
                nodes = [node("a", i * 130, i * 130 + 100, strand, "chr1", start, end)
                         for i, (start, end) in enumerate(spans)]
                event = _build_bnd_event((0, 2), nodes, "AMP")
                event["endpoints"] = node_pair_endpoints(nodes, (0, 2))[0]
                self.assertEqual([j["node_pair"] for j in source_junctions(event, nodes)], [(1, 2)])
                for method in ("read", "olc"):
                    with self.subTest(method=method):
                        for alignment in nodes:
                            alignment[0] = f"raw_rescue_{method}_R1"
                        event["rescue_method"] = method
                        self.assertEqual([j["node_pair"] for j in source_junctions(event, nodes)], [(1, 2)])

    def test_entirely_normal_chain_does_not_fall_back_to_outer_bnd(self):
        normal = (("chr1", 100, 200, "+"), ("chr1", 200, 300, "+"), ("chr1", 300, 400, "+"))
        context, nodes, (pair,) = make_context([chain("normal", *normal)], [("PATH", {0: 1})], [20.])
        self.assertEqual(bnd_calls(context, nodes), [])
        self.assertEqual(context.totals[pair], 20.)

    def test_synthetic_pair_preserves_explicit_geometry_without_alignment_evidence(self):
        nodes = chain("debug_forced_nclose_1", ("chr1", 100, 200, "+"), ("chr1", 200, 300, "+"))
        event = _build_bnd_event((0, 1), nodes, "compressed_nclose_nodes_list.txt")
        event["endpoints"] = node_pair_endpoints(nodes, (0, 1))[0]
        junction, = source_junctions(event, nodes)
        self.assertEqual(junction["mode"], "SYNTHETIC_OUTER")
        self.assertEqual(junction["endpoints"], event["endpoints"])

    def test_threshold_applies_after_shared_primitive_bnd_sum(self):
        other = (ABC[0], ABC[1], ("chr4", 700, 800, "+"))
        context, nodes, _ = make_context([chain("a", *ABC), chain("b", *other)],
                                        [("PATH", {0: 1}), ("PATH", {1: 1})], [.6, .6])
        self.assertEqual(context.visible_ncloses(), [])
        call, = bnd_calls(context, nodes)
        self.assertEqual(call["weight"], .12)
        self.assertEqual(call["endpoints"], (("chr1", 200, "L"), ("chr2", 300, "R")))

    def test_serialized_rescue_mapping_weights_mates_and_depth_ratio_use_primitive_bnds(self):
        other = (ABC[0], ABC[1], ("chr4", 700, 800, "+"))
        context, nodes, _ = make_context([chain("raw_rescue_read_R1", *ABC), chain("raw_rescue_olc_R2", *other)],
                                        [("PATH", {0: 1}), ("PATH", {1: 1})], [.6, .6])
        ns = export_namespace()
        ratio_requests = []
        class Ratios:
            def pair(self, endpoints, weight):
                ratio_requests.append((endpoints, weight))
                return {ns["BP_STEP_DEPTH_RATIO_B"]: [None, None],
                        ns["BP_STEP_DEPTH_RATIO_PREDICT_B"]: [None, None]}
        before = pickle.dumps((context.__dict__, nodes))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calls.vcf"
            result = write_native_vcf(context, {f"chr{i}": 1000 for i in range(1, 5)}, path, Ratios(),
                                      ns["build_vcf_header"], ns["write_bnd_vcf_pair"], ns["write_symbolic_vcf_record"],
                                      nodes=nodes)
            with vcfpy.Reader.from_path(path) as reader:
                records = list(reader)
            with path.with_suffix(".nclose_bnds.tsv").open() as handle:
                links = list(csv.DictReader(handle, delimiter="\t"))
            with path.with_suffix(".bnd_weights.tsv").open() as handle:
                contributions = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(result, (1, 2))
        self.assertEqual(len(links), 4)
        self.assertEqual({link["mode"] for link in links}, {"ADJACENT_ALIGNMENT"})
        self.assertEqual(collections.Counter(link["bnd_id"] for link in links), {"SKYPE.BND.1": 2, ".": 2})
        self.assertEqual({link["source_nclose_key"] for link in links}, {"0:2", "3:5"})
        self.assertAlmostEqual(sum(float(row["contribution_N"]) for row in contributions), .12)
        for row in contributions:
            matching = [link for link in links if link["bnd_id"] == row["bnd_id"]
                        and link["source_nclose_key"] == row["source_nclose_key"]
                        and link["junction_index"] in row["junction_indices"].split(";")]
            self.assertEqual(len(matching), int(row["bnd_occurrences_per_nclose"]))
        a, b = records
        self.assertEqual((a.CHROM, a.POS, b.CHROM, b.POS), ("chr1", 200, "chr2", 301))
        self.assertEqual(a.INFO["MATEID"], b.ID[0])
        self.assertEqual(b.INFO["MATEID"], a.ID[0])
        self.assertEqual(a.ALT[0].mate_pos, b.POS)
        self.assertEqual(a.INFO["NCLOSE_KEYS"], ["0:2", "3:5"])
        self.assertAlmostEqual(a.INFO["WEIGHT"], .12)
        self.assertAlmostEqual(b.INFO["WEIGHT"], .12)
        self.assertEqual(ratio_requests[0][0], [("chr1", 200, "left"), ("chr2", 300, "right")])
        self.assertAlmostEqual(ratio_requests[0][1], 1.2)
        self.assertEqual(pickle.dumps((context.__dict__, nodes)), before)


if __name__ == "__main__":
    unittest.main()
