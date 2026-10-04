"""Physical ports, assembly context and historical input-binding boundaries."""
from collections import Counter
import copy
import hashlib
from itertools import product
import json
from pathlib import Path
import pickle
import tempfile
from types import SimpleNamespace
import unittest

import networkx as nx

from terminal_evidence import (
    BINDING_SCHEMA, begin_component_capture, component_inputs_snapshot, capture_component, display_class, project_source_boundary,
    sha256, validate_upstream_binding, write_component_context, write_terminal_evidence,
)
from terminal_repeat import describe_source, repeat_tracts, reverse_complement
from test_depth_paf import load_bnd_raw_contig_list
from test_touching_reference_depth import load_renderer, row as aligned_row, traversal


def piece(qstart, qend, start=10000, chrom="chr1", strand="+", index=0, cs=None):
    return dict(query="query", qlen=10000, qstart=qstart, qend=qend, strand=strand,
                chrom=chrom, target_length=100000, start=start, end=start + qend - qstart,
                mapq=60, row_index=index, cs=cs or f":{qend-qstart}")


def signature(path):
    return dict(path=str(Path(path).resolve()), sha256=sha256(path))


class TerminalSequenceTests(unittest.TestCase):
    def test_direct_repeat_reverse_complement_and_internal_flank(self):
        sequence = "ACGT" * 500 + "TTAGGG" * 200 + "ACGT" * 250
        host, distal = piece(0, 2000), piece(3200, 4200, start=40000, index=1)
        result = describe_source(sequence, 2000, 1, [host, distal], host)
        self.assertEqual(result["context"], "direct_repeat_extension")
        self.assertTrue(result["internal_repeat_with_farther_aligned_flank"])
        n = len(sequence)
        mirrored = [dict(p, qstart=n-p["qend"], qend=n-p["qstart"], strand="-") for p in [host, distal]]
        rc_result = describe_source(reverse_complement(sequence), n-2000, -1, mirrored, mirrored[0])
        self.assertEqual(rc_result["context"], result["context"])
        self.assertTrue(rc_result["internal_repeat_with_farther_aligned_flank"])
        self.assertEqual(rc_result["adjacent_tracts"][0]["outward_coverage"], result["adjacent_tracts"][0]["outward_coverage"])
        self.assertEqual(display_class("resolved", "assembly_evidence_bound", result["context"],
                                      "no_qualifying_outward_reference_tract"), "observed_nonreference_telomeric_extension")

    def test_donor_route_keeps_order_and_short_pieces(self):
        sequence = "ACGT" * 875 + "CCCTAA" * 200
        host = piece(0, 2000)
        short = piece(2000, 2200, start=50000, chrom="chr2", index=1)
        donor = piece(2200, 3500, start=70000, chrom="chr3", index=2)
        result = describe_source(sequence, 2000, 1, [host, short, donor], host)
        self.assertEqual(result["context"], "repeat_reached_through_ordered_aligned_source_pieces")
        route = result["ordered_repeat_routes"][0]
        self.assertEqual([p["row_index"] for p in route["ordered_source_pieces"]], [0, 1, 2])
        self.assertEqual(len(route["ordered_transitions"]), 2)

    def test_reference_fragmentation_does_not_invent_a_donor(self):
        sequence = "ACGT" * 875 + "CCCTAA" * 200
        whole = piece(0, 3500)
        first, last = piece(0, 2000), piece(2000, 3500, start=12000, index=1)
        unsplit = describe_source(sequence, 2000, 1, [whole], whole)
        split = describe_source(sequence, 2000, 1, [first, last], first)
        self.assertEqual(split["context"], "repeat_reached_through_reference_compatible_source_pieces")
        self.assertTrue(all(t["reference_continuation_compatible"] for t in split["ordered_repeat_routes"][0]["ordered_transitions"]))
        for result in (unsplit, split):
            self.assertNotEqual(display_class("resolved", "assembly_evidence_bound", result["context"],
                                              "no_qualifying_outward_reference_tract"), "donor_route_repeat")

    def test_smoothed_flank_overlap_is_not_canonical_microhomology(self):
        sequence = "CCCTAA" * 200 + "ACGT" * 7 + "CCCTGA" + "ACGT" * 500
        host = piece(1200, len(sequence))
        result = describe_source(sequence, 1200, -1, [host], host)
        tract = next(t for t in result["adjacent_tracts"] if t["tract"]["primary"])
        overlap = tract["aligned_side_overlap_coverage"]
        self.assertGreater(overlap["span_bp"], overlap["canonical_bases"])
        self.assertGreater(tract["outward_coverage"]["canonical_bases"], 1100)

    def test_absent_outward_sequence_has_unknown_precedence(self):
        host = piece(0, 2000)
        result = describe_source("ACGT" * 500, 2000, 1, [host], host)
        self.assertEqual(result["context"], "no_outward_sequence_observed")
        self.assertEqual(display_class("resolved", "assembly_evidence_bound", result["context"],
                                      "no_qualifying_outward_reference_tract"), "unavailable_sequence_or_ambiguous_physical_port")

    def test_whole_tract_canonical_bases_inside_host_do_not_qualify_extension(self):
        # Scientific-review witness: all canonical bases are on the aligned side.
        sequence = "A" * 1000 + "CCCTAA" * 10 + "CCCTGA" * 60 + "A" * 50
        host = piece(0, 1070, start=1000)
        result = describe_source(sequence, 1070, 1, [host], host)
        tract = next(row for row in result["adjacent_tracts"] if row["tract"]["primary"])
        self.assertGreaterEqual(tract["tract"]["canonical_bases"], 60)
        self.assertEqual(tract["outward_coverage"]["canonical_bases"], 0)
        self.assertFalse(tract["outward_primary"])
        self.assertEqual(result["context"], "repeat_context_outward_unqualified")
        self.assertEqual(display_class("resolved", "assembly_evidence_bound", result["context"],
                                      "no_qualifying_outward_reference_tract"), "unavailable_sequence_or_ambiguous_physical_port")

    def test_low_confidence_separator_remains_ambiguous(self):
        sequence = "ACGT" * 600 + "CCCTAA" * 200
        host = piece(0, 2000)
        uncertain = dict(piece(2000, 2400, chrom="chr2", index=1), mapq=0)
        result = describe_source(sequence, 2000, 1, [host, uncertain], host)
        self.assertEqual(result["context"], "repeat_context_ambiguous")
        self.assertEqual(len(result["adjacent_tracts"][0]["low_confidence_intervening_pieces"]), 1)

    def test_insertion_and_deletion_boundary_projections(self):
        for strand in ("+", "-"):
            host = piece(100, 1102, strand=strand, cs=":400+aa:600")
            host["end"] = host["start"] + 1000
            candidates, status = project_source_boundary(host, host["start"] + 400)
            self.assertEqual(len(candidates), 2)
            self.assertEqual(status, "ambiguous_source_boundary")
            deletion = dict(host, qend=1100, end=11002, cs=":400-aa:600")
            self.assertEqual(project_source_boundary(deletion, 10401)[1], "inside_reference_gap")


class TerminalCaptureTests(unittest.TestCase):
    def test_actual_renderer_capture_keeps_zero_span_and_all_strands(self):
        for a, b, walked in product(("+", "-"), repeat=3):
            rows = [aligned_row("first", 10000, 11000, a), aligned_row("last", 11000, 12000, b)]
            raw = [(traversal(a, walked), 0), (traversal(b, walked), 1)]
            renderer = load_renderer(rows)
            unchanged = renderer(raw)
            captured = {}
            self.assertEqual(renderer(raw, terminal_context=captured), unchanged)
            if walked == "-":
                self.assertEqual(unchanged, ([], 2))
                self.assertIsNone(captured["start"])
                self.assertIsNone(captured["end"])
            else:
                self.assertIsNotNone(captured["start"])
                self.assertIsNotNone(captured["end"])
            nodes = renderer.__globals__["contig_data"]
            value = capture_component(0, (0, raw), raw, nodes, captured)
            self.assertEqual(value["endpoints"]["start"]["raw_query_traversal"], raw[0][0])
            self.assertEqual(value["endpoints"]["end"]["host_depth_row_emitted"], walked == "+")

    def test_actual_strandless_resolver_captures_resolved_walk(self):
        rows = [aligned_row("first", 10000, 11000, "+"), aligned_row("last", 11000, 12000, "-")]
        for first, last in ((2, 0), (1, 3), (2, 3)):
            raw = load_bnd_raw_contig_list(rows, nx.DiGraph())(((first, 0), (last, 1)))
            self.assertEqual([i for _, i in raw], [0, 1])
            self.assertTrue(all(direction in (0, 1) for direction, _ in raw))


class TerminalReportTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.prefix = self.root / "result"
        self.prefix.mkdir()
        self.fasta, self.reference = self.root / "source.fa", self.root / "reference.fa"
        sequence = "CCCTAA" * 200 + "ACGT" * 1000 + "TTAGGG" * 200
        self.fasta.write_text(">query\n" + sequence + "\n")
        self.reference.write_text(">chr1\n" + "ACGT" * 25000 + "\n")
        self.paf, self.raw, self.alt = self.root / "source.aln.paf", self.root / "source.paf", self.root / "source.alt.paf"
        self.paf.write_text("query\t6400\t1200\t5200\t+\tchr1\t100000\t10000\t14000\t4000\t4000\t60\tcs:Z::4000\n")
        self.raw.write_bytes(self.paf.read_bytes())
        self.alt.write_text("")
        self.nodes = [["query", 6400, 1200, 5200, "+", "chr1", 100000, 10000, 14000,
                       60, 3, 0, 0, "chr1", "f", "chr1f", "0", "0", "0", "+", "chr1", "0.0"],
                      ["virtual_contig_chr2", 1000, 0, 1000, "+", "chr2", 100000, 0, 1000,
                       0, 3, 1, 1, "chr2", "f", "chr2f", "0", "0", "0", "+", "chr2", "2.0"]]
        source = str(self.prefix / "20_depth/chr1f_chr1b/1.paf")
        self.structure = dict(structure_id="SKYPE.STRUCTURE.1", feature_index=0, kind="PATH", source=source,
                              telomere_counts=Counter({0: 2}), raw_weight=0.0, weight_N=0.0)
        self.context = SimpleNamespace(structures={self.structure["structure_id"]: self.structure})
        path = [("chr1f", 0, 0), (2, 0), (3, 0), ("chr1b", 0, 1)]
        for name, value in {
            "path_data.pkl": {"chr1f_chr1b": [(path,)]},
            "contig_pat_vec_data.pkl": ([(source, [0])], [0], {0: (2, 0)}, [0]),
            "paf_file_path.pkl": [str(self.paf)],
            "01_nclose_data.pkl": {"contig_data": self.nodes},
        }.items():
            (self.prefix / name).write_bytes(pickle.dumps(value))
        self.ppc = self.prefix / "source.ppc.paf"
        self.ppc.write_text("preserved historical node bytes\n")
        (self.prefix / "telomere_connected_list.txt").write_text("chr1f\t(1, 0, 0)\nchr1b\t(0, 0, 0)\nchr2f\t(1, 1, 0)\n")
        record = capture_component(0, (2, 0), [(1, 0)], self.nodes,
                                   {"start": self.paf.read_text().strip(), "end": self.paf.read_text().strip()})
        self.write_capture(record)
        self.index, self.mapper = self.root / "reference.mmi", self.root / "minimap2"
        self.index.write_bytes(b"test-only reference index")
        self.mapper.write_bytes(b"test-only minimap2")
        inputs = dict(fasta=signature(self.fasta), reference=signature(self.reference))
        outputs = dict(primary_paf=signature(self.raw), alternate_paf=signature(self.alt))
        index_binding = dict(schema="SKYPE.reference_index_source.v1", status="complete", reference=signature(self.reference),
            index=signature(self.index), preset="asm20", minimap2=dict(signature(self.mapper), version="fixture"),
            generation_output=signature(self.index),
            generation_command=[str(self.mapper), "-x", "asm20", "-d", str(self.index), str(self.reference)])
        self.binding = self.root / "source.paf.source_binding.json"
        self.binding.write_text(json.dumps(dict(schema=BINDING_SCHEMA, status="complete",
            producer_stage="raw_and_alternate_alignment_generation", inputs_before=inputs, inputs_after=inputs,
            outputs=outputs, outputs_at_generation=outputs, reference_index_binding=index_binding,
            generation_commands=[[str(self.mapper), str(self.index), str(self.fasta), "-o", str(self.raw)]])))
        Path(str(self.paf) + ".alignasm.json").write_text(json.dumps(dict(inputs=outputs, output=signature(self.paf))))

    def write_capture(self, record):
        write_component_context(self.prefix, [record], self.ppc,
                                expected_inputs=component_inputs_snapshot(self.prefix, self.ppc))

    def report(self, binding=True):
        return write_terminal_evidence(self.prefix, self.context, self.nodes,
            source_fasta=str(self.fasta), reference_fasta=str(self.reference), raw_paf=str(self.raw),
            source_binding=str(self.binding) if binding else str(self.root / "no_binding.json"))

    def test_all_ends_zero_weights_graph_only_and_byte_invariance(self):
        # Both uses of one zero-weight node are separate physical occurrences.
        for name in ("matrix.npz", "weight.npy", "predict_B.npy", "SV_call_result.vcf", "SV_call_result.bed"):
            (self.prefix / name).write_bytes(b"preserve every byte\x00\n")
        before = {p.name: sha256(p) for p in self.prefix.iterdir() if p.is_file()}
        objects_before = pickle.dumps((self.context.structures, self.nodes))
        result = self.report()
        self.assertEqual(pickle.dumps((self.context.structures, self.nodes)), objects_before)
        self.assertEqual(result["counts"]["terminal_occurrences"], 2)
        self.assertEqual(result["counts"]["graph_hosts"], 3)
        self.assertEqual(result["counts"]["graph_only_hosts"], 1)
        self.assertEqual(result["counts"]["physical_source_contexts"], 2)
        self.assertEqual({r["display_class"] for r in result["contexts"]}, {"observed_nonreference_telomeric_extension"})
        self.assertTrue(all(sha256(self.prefix / name) == digest for name, digest in before.items()))
        self.assertFalse(Path(str(self.fasta) + ".fai").exists())
        self.assertFalse(Path(str(self.reference) + ".fai").exists())

    def test_legacy_readable_sequence_does_not_become_bound(self):
        result = self.report(binding=False)
        for row in result["contexts"]:
            self.assertEqual(row["assembly_evidence_state"], "assembly_evidence_unverified")
            self.assertEqual(row["source_repeat_context"], "direct_repeat_extension")
            self.assertEqual(row["display_class"], "unavailable_sequence_or_ambiguous_physical_port")

    def test_preprocessing_query_alias_keeps_original_sequence_provenance(self):
        self.nodes[0][0] = "telomere_middle_cut_contig_1"
        record = capture_component(0, (2, 0), [(1, 0)], self.nodes,
            {"start": self.paf.read_text().strip(), "end": self.paf.read_text().strip()})
        self.write_capture(record)
        result = self.report()
        for row in result["contexts"]:
            self.assertEqual(row["graph_query_label"], "telomere_middle_cut_contig_1")
            self.assertEqual(row["source_query"], "query")
            self.assertTrue(row["source_query_alias_resolved"])
            self.assertEqual(row["assembly_evidence_state"], "assembly_evidence_bound")

    def test_alias_direction_offset_and_length_mismatch_remain_unknown(self):
        original = copy.deepcopy(self.nodes[0])
        for field, value in ((4, "-"), (2, 1201), (1, 6399)):
            with self.subTest(field=field):
                self.nodes[0] = list(original)
                self.nodes[0][0] = "an_alias_with_no_parseable_name_convention"
                self.nodes[0][field] = value
                record = capture_component(0, (2, 0), [(1, 0)], self.nodes,
                    {"start": self.paf.read_text().strip(), "end": self.paf.read_text().strip()})
                self.write_capture(record)
                result = self.report()
                self.assertTrue(all(row["assembly_evidence_state"] in (
                    "source_query_frame_mismatch", "source_alignment_identity_mismatch") for row in result["contexts"]))
                self.assertTrue(all(row["display_class"] == "unavailable_sequence_or_ambiguous_physical_port" for row in result["contexts"]))

    def test_content_mismatch_same_length_and_cached_alignasm_is_insufficient(self):
        original = self.fasta.read_bytes()
        self.fasta.write_bytes(original.replace(b"ACGT", b"TGCA", 1))
        result = self.report()
        self.assertTrue(all(row["source_binding_status"] == "source_or_reference_content_mismatch" for row in result["contexts"]))
        self.assertTrue(all(row["display_class"] == "unavailable_sequence_or_ambiguous_physical_port" for row in result["contexts"]))

    def test_old_reference_index_signature_cannot_certify_fresh_raw_mapping(self):
        data = json.loads(self.binding.read_text())
        del data["reference_index_binding"]
        self.binding.write_text(json.dumps(data))
        result = self.report()
        self.assertTrue(all(row["source_binding_status"] == "unverified_reference_index_source_binding" for row in result["contexts"]))

    def test_missing_stale_or_wrong_raw_source_has_explicit_unknown_state(self):
        self.raw.write_text("different raw input\n")
        result = self.report()
        self.assertTrue(all(row["source_binding_status"] == "raw_or_alternate_content_mismatch" for row in result["contexts"]))
        self.fasta.unlink()
        result = self.report()
        self.assertTrue(all(row["assembly_evidence_state"] == "unavailable" for row in result["contexts"]))
        self.assertTrue(all("source_evidence" not in row for row in result["contexts"]))

    def test_missing_capture_does_not_guess_from_strandless_graph_ports(self):
        (self.prefix / "terminal_component_context.json").unlink()
        result = self.report()
        self.assertEqual(result["counts"]["terminal_occurrences"], 2)
        self.assertTrue(all(row["physical_endpoint"] is None for row in result["contexts"]))
        self.assertTrue(all(row["display_class"] == "unavailable_sequence_or_ambiguous_physical_port" for row in result["contexts"]))

    def test_clipped_reference_end_is_not_a_complete_negative_comparator(self):
        # Keep source binding valid, but move the physical/source locus to 20 bp
        # from the chromosome start before these fixture alignments are bound.
        self.nodes[0][7:9] = [20, 4020]
        self.paf.write_text(self.paf.read_text().replace("\t10000\t14000\t", "\t20\t4020\t"))
        self.raw.write_bytes(self.paf.read_bytes())
        record = capture_component(0, (2, 0), [(1, 0)], self.nodes,
            {"start": self.paf.read_text().strip(), "end": self.paf.read_text().strip()})
        self.write_capture(record)
        binding = json.loads(self.binding.read_text())
        binding["outputs"]["primary_paf"] = signature(self.raw)
        binding["outputs_at_generation"] = binding["outputs"]
        self.binding.write_text(json.dumps(binding))
        Path(str(self.paf) + ".alignasm.json").write_text(json.dumps(dict(inputs=binding["outputs"], output=signature(self.paf))))
        result = self.report()
        first = next(row for row in result["contexts"] if row["physical_endpoint"][2] == "R")
        self.assertEqual(first["assembly_evidence_state"], "assembly_evidence_bound")
        self.assertEqual(first["reference_comparison_status"], "reference_end_truncated_comparison")
        self.assertNotEqual(first["display_class"], "observed_nonreference_telomeric_extension")

    def test_source_PAF_mutation_during_component_capture_is_not_published(self):
        before = begin_component_capture(self.prefix, self.ppc)
        self.assertFalse((self.prefix / "terminal_component_context.json").exists())
        self.paf.write_text(self.paf.read_text() + "changed input\n")
        with self.assertRaisesRegex(ValueError, "inputs changed during stage 21"):
            write_component_context(self.prefix, [], self.ppc, expected_inputs=before)
        self.assertFalse((self.prefix / "terminal_component_context.json").exists())

    def test_current_output_hash_does_not_certify_changed_generation_output(self):
        data = json.loads(self.binding.read_text())
        data["outputs_at_generation"]["primary_paf"]["sha256"] = "0" * 64
        self.binding.write_text(json.dumps(data))
        result = self.report()
        self.assertTrue(all(row["source_binding_status"] == "alignment_output_generation_binding_mismatch" for row in result["contexts"]))


if __name__ == "__main__":
    unittest.main()
