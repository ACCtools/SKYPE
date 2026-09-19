import collections
import csv
from pathlib import Path
import tempfile
import unittest

import vcfpy

from native_bnd import node_pair_endpoints
from native_structure_output import bnd_calls, display_events, write_native_bed, write_native_vcf
from nclose_tracking import _build_bnd_event
from structure_nclose import (
    StructureWeights, add_bnd_metadata, add_virtual_structure, count_path_members,
    set_path_splits, write_structure_reports,
)
from test_native_bnd import export_namespace, node


PAIR = (0, 1)
OTHER = (2, 3)
TYPE4 = ("front_jump", 1, -1)


def fixture(structures):
    nodes = [node("u1", 0, 100, "+", "chr1", 100, 200),
             node("u1", 100, 200, "-", "chr1", 400, 500),
             node("u2", 0, 100, "-", "chr1", 210, 310),
             node("u2", 100, 200, "+", "chr1", 410, 510)]
    ncloses = {}
    for index, key in enumerate((PAIR, OTHER), 1):
        event = _build_bnd_event(key, nodes, "compressed_nclose_nodes_list.txt")
        event.update(nclose_id=f"SKYPE.nclose.{index}")
        event["endpoints"], event["contig_names"] = node_pair_endpoints(nodes, key)
        ncloses[key] = event
    ncloses[TYPE4] = dict(kind="indel", event_key=TYPE4, nclose_id="SKYPE.nclose.3",
                         event_type="front_jump", type2_merge_idx=-1, chrom="chr1", st=700, nd=800,
                         start_chr="chr1", end_chr="chr1", start_pos=700, end_pos=800,
                         start_dir="+", end_dir="+", source="INDEL_INDEX_0")
    for column, s in enumerate(structures):
        s.update(structure_id=f"S{column}", feature_index=column,
                 source=f"feature:{column}", legacy_ids=[], telomere_counts={})
    return nodes, dict(version=1, ncloses=ncloses, structures=structures,
                       column_count=len(structures), next_nclose_id=4, path_splits=[])


class StructureWeightTests(unittest.TestCase):
    def test_path_and_independent_type4_use_own_coefficients(self):
        _, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1, TYPE4: 1}),
                            dict(kind="TYPE4", nclose_counts={TYPE4: 1})])
        context = StructureWeights(model, [20., 10.], 10.)
        self.assertEqual(context.totals[TYPE4], 30.)
        self.assertEqual(context.totals[PAIR], 20.)
        self.assertEqual([r["structure_weight"] for r in context.by_nclose[TYPE4]], [20., 10.])

    def test_path_positions_deduplicate_aliases_but_preserve_repeat_traversals(self):
        path = [(0, "tel"), (1, 0), (1, 1), (1, 0), (1, 1), (0, "end")]
        edge_map = {((1, 0), (1, 1)): "merged"}
        self.assertEqual(count_path_members(path, {PAIR}, edge_map,
                                           {"merged": collections.Counter({PAIR: 1, OTHER: 1})}),
                         {PAIR: 2, OTHER: 2})

    def test_direct_and_compound_exact_alias_at_one_position_is_one_occurrence(self):
        path = [(0, "tel"), (1, 0), (1, 1), (1, 0), (1, 1), (0, "end")]
        counts = count_path_members(path, {PAIR}, {((1, 0), (1, 1)): "merged"},
                                   {"merged": collections.Counter({OTHER: 1})},
                                   {PAIR: PAIR, OTHER: PAIR})
        self.assertEqual(counts, {PAIR: 2})
        # Two physical uses in a compound still count twice per position.
        counts = count_path_members(path, {PAIR}, {((1, 0), (1, 1)): "merged"},
                                   {"merged": collections.Counter({PAIR: 1, OTHER: 1})},
                                   {PAIR: PAIR, OTHER: PAIR})
        self.assertEqual(counts, {PAIR: 2, OTHER: 2})

    def test_virtual_new_source_reuses_existing_fitted_nclose(self):
        nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1})])
        nodes.extend([node("virtual", 0, 100, "+", "chr1", 100, 200),
                      node("virtual", 100, 200, "-", "chr1", 400, 500)])
        layout = dict(ordered_endpoints=(dict(chrom="chr1", coord=200, dir="+", ctg_name="virtual"),
                                        dict(chrom="chr1", coord=500, dir="-", ctg_name="virtual")))
        record = dict(pair_id=1, nclose_key_a=(4, 5), nclose_key_b=OTHER, layout_a=layout, layout_b={})
        add_virtual_structure(model, record, 5., [], nodes)
        context = StructureWeights(model, [10.], 10.)
        call = next(c for c in bnd_calls(context) if PAIR in c["keys"])
        self.assertEqual(call["info"]["NCLOSE_IDS"], ["SKYPE.nclose.1"])
        self.assertEqual(call["weight"], 1.5)
        self.assertEqual(call["info"]["MODEL_WEIGHT"], 1.)
        self.assertEqual(call["info"]["VIRTUAL_WEIGHT"], .5)

    def test_compounds_and_repeated_constituents_add_once_per_structure(self):
        _, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1}),
                            dict(kind="AMP", nclose_counts={PAIR: 2, OTHER: 1}),
                            dict(kind="MERGE_TYPE4", nclose_counts={PAIR: 1, OTHER: 1})])
        context = StructureWeights(model, [10., 20., 30.], 10.)
        self.assertEqual(context.totals[PAIR], 80.)
        self.assertEqual(context.totals[OTHER], 50.)
        call = next(c for c in bnd_calls(context) if PAIR in c["keys"])
        self.assertEqual(call["info"]["PARENT_MULTIPLICITY"], [1, 2, 1])
        self.assertEqual(call["weight"], 8.)

    def test_virtual_is_additive_and_two_display_arms_are_one_structure(self):
        nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1})])
        record = dict(pair_id=1, nclose_key_a=PAIR, nclose_key_b=OTHER,
                      layout_a={}, layout_b={})
        add_virtual_structure(model, record, 8.6, [("chr1", 0, 200), ("chr2", 0, 300)], nodes)
        context = StructureWeights(model, [12.], 10.)
        call = next(c for c in bnd_calls(context) if PAIR in c["keys"])
        self.assertAlmostEqual(call["weight"], 2.06)
        self.assertEqual(call["info"]["MODEL_WEIGHT"], 1.2)
        self.assertAlmostEqual(call["info"]["VIRTUAL_WEIGHT"], .86)
        self.assertEqual(len(context.cn_lists()[0]), 2)

    def test_threshold_applies_after_sum_in_all_nclose_views(self):
        nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1}),
                                dict(kind="PATH", nclose_counts={PAIR: 1})])
        context = StructureWeights(model, [.6, .6], 10.)
        self.assertEqual(len(bnd_calls(context)), 1)
        self.assertEqual(bnd_calls(context)[0]["weight"], .12)
        self.assertEqual(len(display_events(context, nodes, {})), 1)
        self.assertEqual(context.cn_lists()[0], [.12])

    def test_exact_geometry_reuses_identity_and_preserves_both_occurrences(self):
        _, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1, OTHER: 1})])
        model["ncloses"][OTHER]["endpoints"] = model["ncloses"][PAIR]["endpoints"]
        context = StructureWeights(model, [10.], 10.)
        call, = bnd_calls(context)
        self.assertEqual(call["weight"], 2.)
        self.assertEqual(call["info"]["NCLOSE_IDS"], ["SKYPE.nclose.1"])
        self.assertEqual(call["info"]["NCLOSE_KEYS"], ["0:1", "2:3"])
        self.assertEqual(call["info"]["PARENT_MULTIPLICITY"], [2])
        self.assertEqual(context.totals[PAIR], 20.)
        self.assertNotIn(OTHER, context.ncloses)
        self.assertEqual(context.cn_lists()[0], [2.])
        self.assertEqual(model["nclose_aliases"][OTHER], PAIR)

    def test_restored_amp_reuses_exact_reverse_order_geometry_without_allocating_id(self):
        nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1})])
        nodes.extend([node("amp", 0, 100, "+", "chr1", 450, 500),
                      node("amp", 100, 200, "-", "chr1", 150, 200)])
        before = model["next_nclose_id"]
        canonical = add_bnd_metadata(model, (4, 5), nodes, "AMP")
        self.assertEqual(canonical, PAIR)
        self.assertEqual(model["next_nclose_id"], before)
        self.assertEqual(model["nclose_sources"][(4, 5)]["contig_names"], ("amp",))

    def test_nearby_coordinate_or_different_retained_side_is_not_reused(self):
        nodes, model = fixture([])
        nodes.extend([node("near", 0, 100, "+", "chr1", 100, 201),
                      node("near", 100, 200, "-", "chr1", 400, 500),
                      node("side", 0, 100, "-", "chr1", 200, 300),
                      node("side", 100, 200, "+", "chr1", 500, 600)])
        self.assertEqual(add_bnd_metadata(model, (4, 5), nodes, "AMP"), (4, 5))
        self.assertEqual(add_bnd_metadata(model, (6, 7), nodes, "AMP"), (6, 7))
        self.assertNotEqual(model["ncloses"][(4, 5)]["nclose_id"], model["ncloses"][PAIR]["nclose_id"])

    def test_compressed_identity_is_preferred_to_earlier_amp_catalog_entry(self):
        _, model = fixture([dict(kind="AMP", nclose_counts={PAIR: 1})])
        model["ncloses"][PAIR]["source"] = "ecdna_circuit_data.pkl"
        model["ncloses"][OTHER]["endpoints"] = tuple(reversed(model["ncloses"][PAIR]["endpoints"]))
        context = StructureWeights(model, [10.], 10.)
        self.assertEqual(model["nclose_aliases"][PAIR], OTHER)
        self.assertEqual(bnd_calls(context)[0]["info"]["NCLOSE_IDS"], ["SKYPE.nclose.2"])

    def test_split_removes_only_its_source_among_exact_aliases(self):
        _, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1, OTHER: 1}),
                            dict(kind="AMP", nclose_counts={OTHER: 1})])
        model["ncloses"][OTHER]["endpoints"] = model["ncloses"][PAIR]["endpoints"]
        key = (1, 1, "chr1", 200, "+", "chr2", 300, "+", "u1")
        set_path_splits(model, {key: {0: 1}}, {1: {0: 1}}, {1: PAIR})
        context = StructureWeights(model, [10., 20.], 10.)
        self.assertEqual(context.totals[PAIR], 40.)
        self.assertEqual(sorted(context.projected_totals.values()), [10., 30.])
        calls = bnd_calls(context)
        self.assertEqual(sorted(c["weight"] for c in calls), [1., 3.])
        self.assertTrue(all(c["info"]["NCLOSE_IDS"] == ["SKYPE.nclose.1"] for c in calls))

    def test_alias_sums_pass_threshold_once_in_report_bed_and_vcf(self):
        nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1}),
                                dict(kind="PATH", nclose_counts={OTHER: 1})])
        model["ncloses"][OTHER]["endpoints"] = model["ncloses"][PAIR]["endpoints"]
        context = StructureWeights(model, [.6, .6], 10.)
        self.assertEqual(bnd_calls(context)[0]["weight"], .12)
        self.assertEqual(len(display_events(context, nodes, {})), 1)
        self.assertEqual(context.cn_lists()[0], [.12])
        with tempfile.TemporaryDirectory() as tmp:
            write_structure_reports(tmp, context, {})
            def read(name):
                with (Path(tmp)/name).open() as handle:
                    return list(csv.DictReader(handle, delimiter="\t"))
            reports = read("nclose_report.tsv")
            sources = read("nclose_sources.tsv")
            usage = read("structure_nclose_usage.tsv")
        self.assertEqual([r["nclose_id"] for r in reports], ["SKYPE.nclose.1", "SKYPE.nclose.3"])
        alias = next(r for r in sources if r["source_nclose_key"] == "2:3")
        self.assertEqual(alias["source_nclose_id"], "SKYPE.nclose.2")
        self.assertEqual(alias["nclose_id"], "SKYPE.nclose.1")
        self.assertEqual(alias["is_alias"], "1")
        self.assertEqual(usage[1]["source_nclose_keys"], "2:3")
        self.assertEqual(usage[1]["source_occurrence_counts"], "1")

    def test_split_transfers_only_corresponding_path_occurrences(self):
        _, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 2}),
                            dict(kind="AMP", nclose_counts={PAIR: 1})])
        key = (1, 1, "chr1", 200, "+", "chr2", 300, "+", "u1")
        set_path_splits(model, {key: {0: 1}}, {1: {0: 1}}, {1: PAIR})
        context = StructureWeights(model, [10., 20.], 10.)
        self.assertEqual(context.totals[PAIR], 40.)
        self.assertEqual(sorted(context.projected_totals.values()), [10., 30.])
        self.assertEqual(context.cn_lists()[0], [4.])

    def test_amp_bed_has_structure_unit_and_nclose_total(self):
        nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1}),
                                dict(kind="AMP", nclose_counts={PAIR: 1, OTHER: 1},
                                     span=("chr1", 100, 510))])
        context = StructureWeights(model, [10., 20.], 10.)
        events = display_events(context, nodes, {})
        amp, = [e for e in events if e["kind"] == "AMP"]
        self.assertEqual(amp["weight_N"], 2.)
        self.assertEqual(amp["weight_scope"], "STRUCTURE")
        with tempfile.TemporaryDirectory() as tmp:
            write_native_bed(tmp, events)
            with open(Path(tmp)/"SKYPE_result.bed") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([float(r["weight (N)"]) for r in rows if r["type"] == "Amplicon"], [2.])

    def test_preprocessing_exclusion_does_not_hide_other_structure_support(self):
        _, model = fixture([dict(kind="AMP", nclose_counts={PAIR: 1})])
        context = StructureWeights(model, [10.], 10.)
        history = {"stages": {"base": {"reasons": {PAIR: "FILTERED_02_NO_ELIGIBLE_PATH"}}}}
        with tempfile.TemporaryDirectory() as tmp:
            write_structure_reports(tmp, context, history)
            with open(Path(tmp)/"nclose_report.tsv") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(float(rows[0]["nclose_cn"]), 1.)
        self.assertEqual(rows[0]["preprocessing_reason"], "FILTERED_02_NO_ELIGIBLE_PATH")
        self.assertEqual(float(rows[1]["nclose_cn"]), 0.)

    def test_native_serialization_and_contribution_sum(self):
        nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1, TYPE4: 1}),
                                dict(kind="AMP", nclose_counts={PAIR: 1, OTHER: 1})])
        key = (1, 1, "chr1", 200, "+", "chr2", 300, "+", "u1")
        set_path_splits(model, {key: {0: 1}}, {1: {0: 1}}, {1: PAIR})
        context = StructureWeights(model, [10., 20.], 10.)
        ns = export_namespace()
        class Ratios:
            def pair(self, *args):
                return {ns["BP_STEP_DEPTH_RATIO_B"]: [None, None],
                        ns["BP_STEP_DEPTH_RATIO_PREDICT_B"]: [None, None]}
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/"calls.vcf"
            write_native_vcf(context, {"chr1": 1000, "chr2": 1000}, path, Ratios(),
                             ns["build_vcf_header"], ns["write_bnd_vcf_pair"], ns["write_symbolic_vcf_record"])
            with vcfpy.Reader.from_path(path) as reader:
                records = list(reader)
            self.assertEqual(len(records), 7)
            with path.with_suffix(".bnd_weights.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
        sums = collections.defaultdict(float)
        for row in rows:
            sums[row["bnd_id"]] += float(row["contribution_N"])
        by_id = {r.ID[0]: r for r in records}
        for record in records:
            if record.INFO["SVTYPE"] == "BND":
                self.assertAlmostEqual(record.INFO["WEIGHT"], sums[record.ID[0].rsplit("_", 1)[0]])
                mate = by_id[record.INFO["MATEID"]]
                self.assertEqual(mate.INFO["MATEID"], record.ID[0])
                self.assertEqual(record.ALT[0].mate_pos, mate.POS)


if __name__ == "__main__":
    unittest.main()
