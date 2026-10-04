"""Molecule support belongs to ordered junctions, not just outer anchors."""
import importlib.util
from pathlib import Path
import sys
from types import SimpleNamespace
import unittest
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SPEC = importlib.util.spec_from_file_location("raw_rescue_support", ROOT / "24_raw_nclose_rescue.py")
rescue = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rescue)
ARGS = SimpleNamespace(min_sv_size=1000, min_mapq=20, min_anchor=500,
                       min_alignment_identity=.9, max_query_overlap=300,
                       breakpoint_cluster=200, min_support=3)


def anchor(chrom, ref, query, length=2000):
    return dict(chrom=chrom, start0=ref, end0=ref+length, strand="+",
                qstart=query, qend=query+length, qlen=10000,
                mapq=60, identity=1., cigar=f"{length}M")


def event(name, chain):
    return dict(rescue.pair_event(chain[0], chain[-1]), name=name, chain=chain,
                minimum_anchor=2000, minimum_mapq=60)


def reverse_chain(chain):
    return [dict(row, strand="-" if row["strand"] == "+" else "+",
                 qstart=10000-row["qend"], qend=10000-row["qstart"])
            for row in chain[::-1]]


class RawRescueSupportTests(unittest.TestCase):
    def test_resolved_internal_template_survives_collinear_outer_anchors(self):
        query = dict(chrom="chr1", start0=9000, end0=15000, query_id="Q001")
        for last_start in [12000, 14000]:
            chain = [anchor("chr1", 10000, 0), anchor("chr2", 30000, 2000),
                     anchor("chr1", last_start, 4000)]
            for oriented in [chain, reverse_chain(chain)]:
                found, reason = rescue.find_outer_event("template", oriented, query, ARGS)
                self.assertIsNotNone(found, reason)
                self.assertEqual(len(rescue.chain_junctions(found['chain'])), 2)

    def test_unmapped_insertion_and_reference_continuation_stay_excluded(self):
        query = dict(chrom="chr1", start0=9000, end0=20000, query_id="Q001")
        for last_ref, last_query in [(12000, 4000), (14000, 4000)]:
            chain = [anchor("chr1", 10000, 0), anchor("chr1", last_ref, last_query)]
            found, _ = rescue.find_outer_event("negative", chain, query, ARGS)
            self.assertIsNone(found)
        small_changes = [anchor("chr1", 10000, 0), anchor("chr1", 12010, 2000),
                         anchor("chr1", 14000, 4000)]
        self.assertIsNone(rescue.find_outer_event("small", small_changes, query, ARGS)[0])

    def test_compound_reference_template_is_not_zero_span_standalone_indel(self):
        chain = [anchor("chr1", 10000, 0), anchor("chr1", 300000, 2000, 200000),
                 anchor("chr1", 12000, 202000)]
        candidate = dict(event("compound", chain), status="supported")
        rescue.filter_small_indel_candidates([candidate])
        self.assertEqual(candidate['status'], 'supported')
        self.assertEqual(rescue.rescue_handoff_type(chain), 2)

    def test_compression_keeps_different_internal_templates(self):
        chain = [anchor("chr1", 10000, 0), anchor("chr2", 30000, 2000),
                 anchor("chr5", 10000, 4000)]
        nodes = [["original", r['qlen'], r['qstart'], r['qend'], r['strand'],
                  r['chrom'], 1000000, r['start0'], r['end0'], r['mapq']] for r in chain]
        original = event("original", chain)
        changed = event("changed", [chain[0], anchor("chr3", 30000, 2000), chain[-1]])
        self.assertTrue(rescue.same_event(original, changed, 200))
        self.assertFalse(rescue.compatible_source_chain(changed, nodes, (0, 2), 200))
        reverse = event("reverse", reverse_chain(chain))
        self.assertTrue(rescue.compatible_source_chain(reverse, nodes, (0, 2), 200))

    def test_simple_indel_size_gate_uses_primitive_span(self):
        for span, expected in [(99999, 'below_min_indel_span'), (100000, 'supported')]:
            chain = [anchor("chr1", 10000, 0), anchor("chr1", 12000+span, 2000)]
            candidate = dict(event("simple", chain), status="supported")
            rescue.filter_small_indel_candidates([candidate])
            self.assertEqual(candidate['status'], expected)
            self.assertEqual(rescue.rescue_handoff_type(chain), 4)

    def test_cn_control_removes_only_prediction_condition(self):
        args = rescue.build_parser().parse_args(["unused"])
        coords = [("chr1", i*100000+1) for i in range(80)]
        noise = np.random.default_rng(81).normal(0., 1., 80)
        observed = np.r_[np.full(40, 30.), np.full(40, 60.)] + noise
        predicted = np.r_[np.full(40, 30.), np.full(40, 60.)]
        residual = rescue.detect_depth(coords, observed, predicted, args)
        cn_only = rescue.detect_depth(coords, observed, predicted, args, filter_prediction=False)
        self.assertEqual(len(residual[0]), 0)
        self.assertEqual(len(cn_only[0]), 1)
        self.assertEqual(len(residual[2]), len(cn_only[2]))
        for left, right in zip(residual[2], cn_only[2]):
            for field in ['representative0', 'observed_step', 'z2', 'hom', 'pred_range']:
                self.assertEqual(left[field], right[field])
        unexplained = np.full_like(observed, 30.)
        self.assertEqual(rescue.detect_depth(coords, observed, unexplained, args)[:3],
                         rescue.detect_depth(coords, observed, unexplained, args, filter_prediction=False)[:3])

    def test_same_outer_anchors_cannot_pool_different_interiors(self):
        events = [event(f"read{i}", [anchor("chr1", 10000, 0),
                                    anchor(chrom, 20000, 2000),
                                    anchor("chr5", 30000, 4000)])
                  for i, chrom in enumerate(["chr2", "chr3", "chr4"])]
        groups = rescue.cluster_events(events, 200)
        self.assertEqual(sorted(len(g) for g in groups), [1, 1, 1])

    def test_reverse_complement_and_reference_fragmentation_keep_support(self):
        chain = [anchor("chr1", 10000, 0), anchor("chr2", 20000, 2000),
                 anchor("chr5", 30000, 4000)]
        fragmented = [anchor("chr1", 10000, 0, 1000),
                      anchor("chr1", 11000, 1000, 1000), *chain[1:]]
        events = [event("forward", chain), event("reverse", reverse_chain(chain)),
                  event("fragmented", fragmented)]
        groups = rescue.cluster_events(events, 200)
        self.assertEqual([len(g) for g in groups], [3])

    def test_olc_checks_internal_junction_in_addition_to_terminals(self):
        chain = [anchor(chrom, 10000, 2000*i)
                 for i, chrom in enumerate(["chr1", "chr2", "chr3", "chr4"])]
        reads = {**{f"left{i}": chain[:2] for i in range(3)},
                 **{f"right{i}": chain[-2:] for i in range(3)},
                 "middle": chain[1:3]}
        support = rescue.junction_support(chain, reads, list(reads), ARGS)
        self.assertEqual(list(map(len, support)), [3, 1, 3])
        candidate = event("assembled", chain)
        rescue.annotate_junction_evidence(candidate, reads, list(reads), ARGS)
        self.assertEqual(candidate["minimum_junction_support"], 1)
        self.assertEqual(candidate["status"], "insufficient_raw_support")

    def test_each_molecule_counts_once_and_linkage_is_separate(self):
        chain = [anchor(chrom, 10000, 2000*i)
                 for i, chrom in enumerate(["chr1", "chr2", "chr3"])]
        reads = {**{f"left{i}": chain[:2] for i in range(3)},
                 **{f"right{i}": chain[1:] for i in range(3)}}
        candidate = event("assembled", chain)
        rescue.annotate_junction_evidence(candidate, reads, list(reads)*2, ARGS)
        self.assertEqual(candidate["junction_support_counts"], [3, 3])
        self.assertEqual(candidate["full_chain_support_reads"], [])
        self.assertEqual(candidate["status"], "supported")
        reads["spanning"] = reverse_chain(chain)
        rescue.annotate_junction_evidence(candidate, reads, list(reads), ARGS)
        self.assertEqual(candidate["full_chain_support_reads"], ["spanning"])


if __name__ == "__main__":
    unittest.main()
