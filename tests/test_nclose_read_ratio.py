import math
import os
import tempfile
import unittest
import unittest.mock
from pathlib import Path

import nclose_read_ratio as rr


def row(name, ctg_start, ctg_end, direction, chrom, ref_start, ref_end, mapq=60):
    data = [None] * 22
    data[rr.CTG_NAM] = name
    data[rr.CTG_LEN] = 100_000
    data[rr.CTG_STR] = ctg_start
    data[rr.CTG_END] = ctg_end
    data[rr.CTG_DIR] = direction
    data[rr.CHR_NAM] = chrom
    data[rr.CHR_STR] = ref_start
    data[rr.CHR_END] = ref_end
    data[rr.CTG_MAPQ] = mapq
    return data


def aln(target, ts, te, score, strand="+", qs=0, qe=None, qlen=20_000):
    qe = qlen if qe is None else qe
    return rr.Alignment(target, strand, qlen, qs, qe, 20_000, ts, te, 0, te - ts, score)


class BreakpointGeometryTests(unittest.TestCase):
    def test_end_a_leaves_forward_anchor_at_reference_end(self):
        self.assertEqual(rr.anchor_breakpoint(row("u", 0, 500, "+", "chr1", 100, 600), "A"), (600, "left"))

    def test_end_a_leaves_reverse_anchor_at_reference_start(self):
        self.assertEqual(rr.anchor_breakpoint(row("u", 0, 500, "-", "chr1", 100, 600), "A"), (100, "right"))

    def test_end_b_enters_forward_anchor_at_reference_start(self):
        self.assertEqual(rr.anchor_breakpoint(row("u", 0, 500, "+", "chr1", 100, 600), "B"), (100, "right"))

    def test_end_b_enters_reverse_anchor_at_reference_end(self):
        self.assertEqual(rr.anchor_breakpoint(row("u", 0, 500, "-", "chr1", 100, 600), "B"), (600, "left"))


def make_site():
    site = rr.Site(0, "chr1", 10_000, "left", 1000, 1000)
    site.ref_target = "REF|chr1|5000|15000"
    site.ref_offset = 5000
    site.local_targets = {"ALT|0"}
    site.ctx_lo, site.ctx_hi = 0, 1_000_000
    return site


ALT = rr.MemberEnd(0, ("u", (0, 1)), ("u", (0, 1)), "A", 0, "chr1", 10_000, "left", 60,
                   alt_seq="N" * 10_000, alt_bp=5000, cont_len=5000, cls="C0")


class ReadCallTests(unittest.TestCase):
    def call(self, alns):
        return rr.call_read(make_site(), alns, {"ALT|0": ALT}, {"ALT|0": "C0"}, margin=50)

    def test_reference_read_crossing_both_sides(self):
        alns = [aln("REF|chr1|5000|15000", 1000, 9000, 8000), aln("ALT|0", 1000, 5000, 4000)]
        self.assertEqual(self.call(alns), ("REF", None))

    def test_junction_read_matches_alt_continuation(self):
        alns = [aln("REF|chr1|5000|15000", 1000, 5000, 4000), aln("ALT|0", 1000, 9000, 8000)]
        self.assertEqual(self.call(alns), ("ALT", "C0"))

    def test_small_margin_is_ambiguous(self):
        alns = [aln("REF|chr1|5000|15000", 1000, 9000, 8000), aln("ALT|0", 1000, 9000, 7990)]
        self.assertEqual(self.call(alns)[0], "AMB")

    def test_read_leaving_reference_into_unknown_sequence_is_other(self):
        alns = [aln("REF|chr1|5000|15000", 1000, 5000, 4000, qs=0, qe=4000, qlen=9000)]
        self.assertEqual(self.call(alns), ("OTHER", None))

    def test_read_ending_at_breakpoint_is_not_informative(self):
        alns = [aln("REF|chr1|5000|15000", 1000, 5000, 4000, qs=0, qe=4000, qlen=4100)]
        self.assertEqual(self.call(alns), ("NA", None))

    def test_other_junction_elsewhere_on_the_read_does_not_compete(self):
        # 30 kb read: this junction at read 5000, another NClose's junction near read 20000.
        alns = [aln("REF|chr1|5000|15000", 0, 5000, 5000, qs=0, qe=5000, qlen=30_000),
                aln("ALT|0", 0, 10_000, 9_990, qs=0, qe=10_000, qlen=30_000),
                aln("ALT|9", 0, 10_000, 9_995, qs=15_000, qe=25_000, qlen=30_000)]
        self.assertEqual(self.call(alns), ("ALT", "C0"))

    def test_longer_reference_window_does_not_mask_a_junction_read(self):
        # Merged reference window: the retained side aligns over 12 kb, the ALT sequence holds 5 kb of it.
        site = make_site()
        site.ref_target, site.ref_offset = "REF|chr1|0|20000", 12_000
        alns = [aln("REF|chr1|0|20000", 0, 12_000, 12_000, qs=0, qe=12_000, qlen=20_000),
                aln("ALT|0", 0, 10_000, 9_990, qs=7_000, qe=17_000, qlen=20_000)]
        self.assertEqual(rr.call_read(site, alns, {"ALT|0": ALT}, {"ALT|0": "C0"}, margin=50), ("ALT", "C0"))

    def test_sequence_carrying_the_same_junction_does_not_compete(self):
        # The other end's ALT of the same junction (inserted sequence between anchors) spans more of the read.
        alns = [aln("REF|chr1|5000|15000", 1000, 5000, 4000, qs=0, qe=4000, qlen=12_000),
                aln("ALT|0", 1000, 9000, 7990, qs=0, qe=8000, qlen=12_000),
                aln("ALT|7", 0, 11_000, 10_990, qs=0, qe=11_000, qlen=12_000)]
        self.assertEqual(self.call(alns)[0], "AMB")
        got = rr.call_read(make_site(), alns, {"ALT|0": ALT}, {"ALT|0": "C0", "ALT|7": "C7"}, margin=50,
                           carriers={"C0": {"ALT|7"}})
        self.assertEqual(got, ("ALT", "C0"))

    def test_carriers_need_the_junction_with_its_flanks(self):
        a = rr.MemberEnd(0, ("n", (0, 1)), ("n", (0, 1)), "B", 1, "chr17", 0, "left", 60,
                         alt_seq="N" * 6127, alt_bp=1127, cont_len=5000, cls="C97")
        b = rr.MemberEnd(1, ("n", (0, 1)), ("n", (0, 1)), "A", 0, "chr4", 0, "right", 60,
                         alt_seq="N" * 10_000, alt_bp=5000, cont_len=5000, cls="C96")
        c = rr.MemberEnd(2, ("m", (2, 3)), ("m", (2, 3)), "B", 3, "chr4", 0, "left", 60,
                         alt_seq="N" * 8036, alt_bp=5000, cont_len=3036, cls="C95")
        alns = {"ALT|0": [rr.Alignment("ALT|1", "-", 6127, 0, 6127, 10_000, 2091, 8218, 6127, 6127, 6127),
                          rr.Alignment("ALT|2", "+", 6127, 0, 1296, 8036, 6736, 8036, 1294, 1300, 1254)]}
        self.assertEqual(rr.junction_carriers([a, b, c], alns), {"C97": {"ALT|1"}})

    def test_reach_past_homology_cannot_be_bypassed_through_another_template(self):
        # This end needs 3 kb past its breakpoint (2 kb homology); the other end's template of the
        # same class needs only 1 kb but is not at this breakpoint.
        site = make_site()
        site.req_cont = 3000
        here = rr.MemberEnd(0, ("u", (0, 1)), ("u", (0, 1)), "A", 0, "chr1", 10_000, "left", 60,
                            alt_seq="N" * 10_000, alt_bp=5000, cont_len=5000, cls="C0", req_cont=3000)
        there = rr.MemberEnd(1, ("u", (0, 1)), ("u", (0, 1)), "B", 1, "chr2", 50_000, "right", 60,
                             alt_seq="N" * 10_000, alt_bp=5000, cont_len=5000, cls="C0", req_cont=1000)
        alns = [aln("REF|chr1|5000|15000", 1000, 11_200, 10_190, qs=0, qe=10_200, qlen=10_200),
                aln("ALT|0", 1000, 6200, 5190, qs=0, qe=5200, qlen=10_200),
                aln("ALT|1", 3800, 9000, 5190, qs=0, qe=5200, qlen=10_200)]
        info = {"ALT|0": here, "ALT|1": there}
        classes = {"ALT|0": "C0", "ALT|1": "C0"}
        self.assertEqual(rr.call_read(site, alns, info, classes, margin=50), ("REF", None))
        self.assertEqual(rr.call_read(site, alns[1:], info, classes, margin=50), ("NA", None))
        # A read reaching past the homology on this end's own template is a junction read.
        reach = [aln("ALT|0", 1000, 8200, 7190, qs=0, qe=7200, qlen=7200)]
        self.assertEqual(rr.call_read(site, reach, info, classes, margin=50), ("ALT", "C0"))

    def test_competing_explanation_of_the_same_read_interval_is_ambiguous(self):
        alns = [aln("REF|chr1|5000|15000", 1000, 5000, 4000, qs=0, qe=4000, qlen=9000),
                aln("ALT|0", 1000, 9000, 7990, qs=0, qe=8000, qlen=9000),
                aln("REF|chr7|0|10000", 0, 8000, 7980, qs=0, qe=8000, qlen=9000)]
        self.assertEqual(self.call(alns)[0], "AMB")


class ReadIndelTests(unittest.TestCase):
    @staticmethod
    def read_aln(target, strand, qlen, qs, qe, ts, te, score, cigar):
        a = rr.Alignment(target, strand, qlen, qs, qe, 30_000, ts, te, 0, te - ts, score)
        a.indels = rr.read_indels(a, cigar)
        return a

    def test_read_alignment_is_cut_at_a_structural_indel_at_the_breakpoint(self):
        whole = self.read_aln("REF|chr1|0|30000", "+", 20_000, 0, 20_000, 2000, 27_000, 20_000,
                              "10000=5000D10000=")
        pieces = [(2000, 12_000, 0, 10_000, 10_000), (17_000, 27_000, 10_000, 20_000, 10_000)]
        for bp in (12_000, 17_000, 14_000):  # either end of the deletion, or inside it
            self.assertEqual([(a.ts, a.te, a.qs, a.qe, a.score) for a in rr.cut_at(whole, bp)], pieces)
        self.assertEqual(rr.cut_at(whole, 25_000), [whole])
        reverse = self.read_aln("REF|chr1|0|30000", "-", 23_000, 0, 23_000, 2000, 22_000, 20_000,
                                "10000=3000I10000=")
        self.assertEqual([(a.ts, a.te, a.qs, a.qe) for a in rr.cut_at(reverse, 12_000)],
                         [(2000, 12_000, 13_000, 23_000), (12_000, 22_000, 0, 10_000)])
        # Sequencing-error sized indels are not structural.
        small = self.read_aln("REF|chr1|0|30000", "+", 24_960, 0, 24_960, 2000, 27_000, 24_900,
                              "10000=40D14960=")
        self.assertEqual(rr.cut_at(small, 12_000), [small])

    def test_germline_indel_near_the_breakpoint_keeps_a_reference_read(self):
        site = make_site()
        near = self.read_aln(site.ref_target, "+", 9500, 0, 9500, 0, 10_000, 9400, "4000=500D5500=")
        self.assertEqual(rr._query_position(near, site.ref_offset), 4500)
        self.assertEqual(rr.call_read(site, [near], {"ALT|0": ALT}, {"ALT|0": "C0"}, margin=50), ("REF", None))
        at = self.read_aln(site.ref_target, "+", 9500, 0, 9500, 0, 10_000, 9400, "5000=500D4500=")
        self.assertEqual(rr.call_read(site, [at], {"ALT|0": ALT}, {"ALT|0": "C0"}, margin=50), ("OTHER", None))

    def test_read_carrying_a_deletion_at_the_breakpoint_is_not_a_reference_read(self):
        # Merged reference window holding both breakpoints of a 5 kb deletion; the read carries it.
        site = make_site()
        site.ref_target, site.ref_offset = "REF|chr1|0|30000", 12_000
        with tempfile.TemporaryDirectory() as tmp:
            paf = Path(tmp) / "reads.paf"
            paf.write_text(
                "read\t20000\t0\t20000\t+\tREF|chr1|0|30000\t30000\t2000\t27000\t20000\t25000\t60\t"
                "AS:i:19970\tcg:Z:10000=5000D10000=\n"
                "read\t20000\t5000\t15000\t+\tALT|0\t10000\t0\t10000\t10000\t10000\t60\t"
                "AS:i:10000\tcg:Z:10000=\n")
            alns = rr.parse_paf(paf)["read"]
        self.assertEqual(rr.call_read(site, alns, {"ALT|0": ALT}, {"ALT|0": "C0"}, margin=50), ("ALT", "C0"))
        # Without its ALT sequence the read leaves the reference at the breakpoint.
        no_alt = [a for a in alns if a.target != "ALT|0"]
        self.assertEqual(rr.call_read(site, no_alt, {"ALT|0": ALT}, {"ALT|0": "C0"}, margin=50), ("OTHER", None))


def ref_aln(qs, qe, ts, te, cigar, strand="+"):
    return rr.Alignment("REF|chr1|5000|15000", strand, 10_000, qs, qe, 10_000, ts, te, 0, qe - qs,
                        qe - qs, cigar)


class JunctionHomologyTests(unittest.TestCase):
    def setUp(self):
        self.site = make_site()
        self.item = rr.MemberEnd(0, ("u", (0, 1)), ("u", (0, 1)), "A", 0, "chr1", 10_000, "left", 60,
                                 alt_seq="N" * 10_000, alt_bp=5000, cont_len=5000, cls="C0", site=0)

    def homology(self, alns):
        return rr.junction_homology(self.item, self.site, {"ALT|0": alns})

    def test_identity_is_measured_past_the_breakpoint_on_either_strand(self):
        fwd = ref_aln(0, 10_000, 0, 10_000, "5000=10X4990=")
        rev = ref_aln(0, 10_000, 0, 10_000, "4990=10X5000=", strand="-")
        self.assertEqual(rr._identity_beyond(fwd, 5000), (4990, 5000))
        self.assertEqual(rr._identity_beyond(rev, 5000), (4990, 5000))

    def test_continuation_copying_the_reference_is_homologous(self):
        self.assertEqual(self.homology([ref_aln(0, 7000, 0, 7000, "7000=")]), 2000)

    def test_diverged_continuation_still_discriminates(self):
        self.assertEqual(self.homology([ref_aln(0, 10_000, 0, 10_000, "5000=250X4750=")]), 0)

    def test_junction_copy_elsewhere_within_read_reach_is_homologous(self):
        # A repeat copy 1 kb away holds the whole ALT: reads from it cross the junction without carrying it.
        self.assertEqual(self.homology([ref_aln(0, 10_000, 1000, 11_000, "10000=")]), 5000)

    def test_junction_copy_in_the_read_reach_context_is_homologous(self):
        self.site.ctx_target = "CTX|chr1|0|60000"
        copy = rr.Alignment("CTX|chr1|0|60000", "+", 10_000, 0, 10_000, 60_000, 40_000, 50_000,
                            9950, 10_000, 9800, "4000=10X5990=")
        self.assertEqual(self.homology([copy]), 5000)

    def test_deletion_at_the_breakpoint_is_not_homology_on_either_strand(self):
        # The reference holds 5 kb between the flank and the continuation: not the contiguous junction.
        for strand in "+-":
            deletion = ref_aln(0, 10_000, 0, 15_000, "5000=5000D5000=", strand=strand)
            self.assertEqual(rr._identity_beyond(deletion, 5000), (5000, 10_000))
            self.assertEqual(self.homology([deletion]), 0)

    def test_copy_with_an_indel_inside_the_required_flank_is_homology(self):
        # The copy lacks 500 bases 0.5-1 kb before the breakpoint.  A read from it still passes as
        # ALT (a germline-sized indel away from the breakpoint does not cut its alignment), so the
        # copy shares the junction as far as read calling is concerned.
        self.assertEqual(self.homology([ref_aln(0, 10_000, 0, 9500, "4000=500I5500=")]), 5000)

    def test_duplicated_unique_sequence_is_not_homology(self):
        # Tandem duplication of unique sequence: the ALT spans the reference only with the copy as an insertion.
        self.assertEqual(self.homology([ref_aln(0, 10_000, 0, 9000, "5000=1000I4000=")]), 0)

    def test_copy_sharing_less_than_the_required_retained_flank_is_not_homology(self):
        # Reads from this copy diverge from the ALT 313 bp before the junction (1 kb must match).
        self.assertEqual(self.homology([ref_aln(4687, 10_000, 1000, 6313, "5313=")]), 0)
        self.assertEqual(self.homology([ref_aln(4100, 10_000, 1000, 6900, "5900=")]), 5000)

    def test_retained_flank_alone_is_not_homology(self):
        # The ALT's own anchor ends at the breakpoint; the continuation aligns elsewhere.
        self.assertEqual(self.homology([ref_aln(0, 5100, 0, 5100, "5100="),
                                        ref_aln(5000, 10_000, 20_000, 25_000, "5000=")]), 0)

    def test_partial_homology_moves_the_required_reach(self):
        rr.apply_junction_homology([self.item], [self.site], {"ALT|0": [ref_aln(0, 7000, 0, 7000, "7000=")]})
        self.assertEqual((self.item.homology, self.site.req_cont, self.item.req_cont), (2000, 3000, 3000))
        self.assertEqual(self.item.reason, "")

    def test_fully_homologous_continuation_is_not_testable(self):
        rr.apply_junction_homology([self.item], [self.site],
                                   {"ALT|0": [ref_aln(0, 10_000, 0, 10_000, "10000=")]})
        self.assertEqual(self.item.reason, "homologous_junction")


class _Reference:
    def __init__(self, length):
        self.length = length

    def get_reference_length(self, chrom):
        return self.length

    def fetch(self, chrom, start, end):
        return "A" * (end - start)


class _Source:
    """Reference pieces named by their row; oriented_reference returns them as is."""

    def oriented_reference(self, row):
        return row[rr.CHR_NAM] * (int(row[rr.CTG_END]) - int(row[rr.CTG_STR]))


class ReferenceAltTests(unittest.TestCase):
    def rows(self, gap):
        return [row("u", 0, 6000, "+", "A", 0, 6000), row("u", 6000 + gap, 12_000 + gap, "+", "B", 0, 6000)]

    def test_adjacent_pieces_make_the_junction(self):
        seq, ret, cont = rr._reference_alt(_Source(), self.rows(0), (0, 1), "A")
        self.assertEqual((seq, ret, cont), ("A" * 5000 + "B" * 5000, 5000, 5000))

    def test_unaligned_insertion_between_pieces_is_not_invented_away(self):
        self.assertEqual(rr._reference_alt(_Source(), self.rows(20_000), (0, 1), "A"), ("", 0, 0))
        self.assertEqual(rr._reference_alt(_Source(), self.rows(20_000), (0, 1), "B"), ("", 0, 0))


class PresetTests(unittest.TestCase):
    def test_combined_minimap2_flags(self):
        for cl, preset in (("minimap2 -ax map-ont ref.fa r.fq", "map-ont"), ("minimap2 -t 8 -x lr:hq a b", "lr:hq"),
                           ("minimap2 -a -x map-pb a b", "map-pb"), ("samtools sort -o x.bam", None)):
            header = {"PG": [{"ID": "x", "CL": cl}]}
            with unittest.mock.patch.object(rr.pysam, "AlignmentFile") as bam:
                bam.return_value.__enter__.return_value.header.to_dict.return_value = header
                self.assertEqual(rr.detect_preset("x.bam"), preset or "map-hifi")


class ContextTargetTests(unittest.TestCase):
    def test_context_follows_the_reach_of_the_fetched_reads(self):
        near = rr.Site(0, "chr1", 50_000, "left", 1000, 1000)
        near.reach_lo, near.reach_hi = 30_000, 75_000
        far = rr.Site(1, "chr1", 500_000, "left", 1000, 1000)
        far.reach_lo, far.reach_hi = 300_000, 900_000
        bare = rr.Site(2, "chr1", 2_000_000, "left", 1000, 1000)
        targets = rr.build_context_targets([near, far, bare], _Reference(3_000_000))
        self.assertEqual(near.ctx_target, "CTX|chr1|25000|80000")
        # Ultra-long reads: the whole reach, however far.
        self.assertEqual(far.ctx_target, "CTX|chr1|295000|905000")
        self.assertEqual(bare.ctx_target, f"CTX|chr1|{2_000_000 - rr.FLANK}|{2_000_000 + rr.FLANK}")
        self.assertEqual(len(targets), 3)

    def test_merged_context_keeps_each_sites_own_reach(self):
        a = rr.Site(0, "chr1", 100_000, "left", 1000, 1000)
        a.reach_lo, a.reach_hi = 95_000, 105_000
        b = rr.Site(1, "chr1", 120_000, "left", 1000, 1000)
        b.reach_lo, b.reach_hi = 100_000, 150_000
        rr.build_context_targets([a, b], _Reference(3_000_000))
        self.assertEqual(a.ctx_target, b.ctx_target)
        copy = rr.Alignment(a.ctx_target, "+", 10_000, 0, 10_000, 65_000, 50_000, 60_000, 10_000, 10_000,
                            10_000, "10000=")  # reference 140-150 kb: only b's reads reach it
        self.assertFalse(rr._within_reach(copy, a))
        self.assertTrue(rr._within_reach(copy, b))


class LocalCopyNumberTests(unittest.TestCase):
    def windows(self):
        # 100 kb windows on chr1 at depth 40 (2N with N = 20) everywhere.
        return {"chr1": [(i * 100_000 + 1, (i + 1) * 100_000, 40.0) for i in range(20)]}

    def site(self, site_id, bp, retained, records):
        site = rr.Site(site_id, "chr1", bp, retained, 1000, 1000)
        site.mapq_records, site.unique_fraction = records, 1.0
        return site

    def test_segment_bounded_by_another_breakpoint_uses_depth_at_the_breakpoint(self):
        # A 1.1 kb segment between two breakpoints carries extra copies the windows cannot see.
        sites = [self.site(i, 300_000 + 400_000 * i, "left", 40) for i in range(5)]
        sites += [self.site(5, 1_550_000, "right", 100), self.site(6, 1_551_100, "left", 100)]
        per_copy = rr.assign_local_cn(sites, self.windows(), 20.0)
        self.assertAlmostEqual(per_copy, 20.0)
        self.assertEqual([s.cn_source for s in sites], ["windows"] * 5 + ["point", "point"])
        self.assertAlmostEqual(sites[5].cn_local, 5.0)

    def test_other_side_of_the_same_breakpoint_does_not_bound(self):
        sites = [self.site(0, 1_000_004, "left", 40), self.site(1, 1_000_000, "right", 40)]
        rr.assign_local_cn(sites, self.windows(), 20.0)
        self.assertEqual([s.cn_source for s in sites], ["windows", "windows"])
        self.assertAlmostEqual(sites[0].cn_local, 2.0)


class ShortSegmentTests(unittest.TestCase):
    def test_short_retained_segment_is_assigned_but_not_estimated(self):
        nclose = ("u", (0, 1))
        sites = [measured_site(0, 2.0, 20, {"C0": 10}), measured_site(1, 2.0, 20, {})]
        ends = [member_end(0, nclose, "A", 0, "C0"),
                rr.MemberEnd(1, nclose, nclose, "B", 1, "chr1", 0, "left", 60, alt_seq="A", alt_bp=1100,
                             cont_len=5000, cls="C1", site=1, reason="short_retained_segment")]
        results, _ = rr.summarize(ends, sites)
        self.assertAlmostEqual(results[nclose]["j"], 10 / 30 * 2.0)


def measured_site(site_id, cn, ref, alt_counts, other=0, unique=1.0, names=None):
    """A measured site; ``names`` (class -> read names) gives reads seen at other sites too."""
    site = rr.Site(site_id, "chr1", 10_000 * (site_id + 1), "left", 1000, 1000)
    site.cn_local = cn
    site.unique_fraction = unique
    informative = ref + sum(alt_counts.values()) + other
    site.counts = {"REF": ref, "ALT": sum(alt_counts.values()), "OTHER": other, "AMB": 0,
                   "informative": informative}
    site.alt_counts = dict(alt_counts)
    site.alt_names = {c: set((names or {}).get(c) or (f"s{site_id}.{c}.{i}" for i in range(n)))
                      for c, n in alt_counts.items()}
    return site


def member_end(end_id, nclose, end, site, cls, mapq=60):
    return rr.MemberEnd(end_id, nclose, nclose, end, 0, "chr1", 0, "left", mapq,
                        alt_seq="A", alt_bp=5000, cont_len=5000, cls=cls, site=site)


class SummaryTests(unittest.TestCase):
    def test_single_read_side_branch_is_removed_before_the_fit(self):
        nclose = ("utg", (0, 1))
        sites = [measured_site(0, 1.0, 17, {"C0": 1}, unique=0.0),
                 measured_site(1, 1.6, 26, {"C0": 1})]
        ends = [member_end(0, nclose, "A", 0, "C0", mapq=1), member_end(1, nclose, "B", 1, "C0")]
        results, reasons = rr.summarize(ends, sites)
        est = results[nclose]
        self.assertEqual(reasons[0], "anchor_not_unique")
        self.assertEqual(est["status"], "TESTABLE")
        self.assertAlmostEqual(est["j"], 1 / 27 * 1.6)
        self.assertEqual(rr.prefit_decision(est), "read_ratio_low")

    def test_shared_class_does_not_borrow_another_ncloses_site(self):
        weak, strong = ("a", (0, 1)), ("b", (2, 3))
        sites = [measured_site(0, 2.0, 5, {}), measured_site(1, 2.0, 40, {})]
        ends = [member_end(0, weak, "A", 0, "C0"), member_end(1, strong, "A", 0, "C0"),
                member_end(2, strong, "B", 1, "C0")]
        results, _ = rr.summarize(ends, sites)
        self.assertAlmostEqual(results[weak]["r_eff"], 2.5)
        self.assertEqual(rr.prefit_decision(results[weak]), "")
        self.assertEqual(rr.prefit_decision(results[strong]), "read_ratio_low")
        self.assertEqual(results[weak]["shared"], "C0")

    def test_another_ncloses_unique_anchor_does_not_place_this_one(self):
        ambiguous, unique = ("a", (0, 1)), ("b", (2, 3))
        sites = [measured_site(0, 2.0, 40, {}), measured_site(1, 2.0, 40, {})]
        ends = [member_end(0, ambiguous, "A", 0, "C0", mapq=0), member_end(1, ambiguous, "B", 1, "C1", mapq=0),
                member_end(2, unique, "A", 0, "C2"), member_end(3, unique, "B", 1, "C3")]
        results, _ = rr.summarize(ends, sites)
        self.assertEqual(results[ambiguous]["status"], "NOT_TESTABLE")
        self.assertEqual(results[ambiguous]["reason"], "anchor_not_unique")
        self.assertEqual(rr.prefit_decision(results[ambiguous]), "")
        self.assertEqual(rr.prefit_decision(results[unique]), "read_ratio_low")

    def test_untestable_ncloses_are_kept(self):
        nclose = ("u", (0, 1))
        sites = [measured_site(0, 2.0, 30, {}, unique=0.2)]
        ends = [member_end(0, nclose, "A", 0, "C0"), member_end(1, nclose, "B", 0, "C0")]
        results, _ = rr.summarize(ends, sites)
        self.assertEqual(results[nclose]["status"], "NOT_TESTABLE")
        self.assertEqual(rr.prefit_decision(results[nclose]), "")

    def test_compressed_distinct_junctions_add_up(self):
        nclose = ("rep", (0, 1))
        sites = [measured_site(0, 2.0, 20, {"C0": 10}), measured_site(1, 2.0, 20, {"C1": 10})]
        ends = [member_end(0, nclose, "A", 0, "C0"),
                rr.MemberEnd(1, nclose, ("dup", (2, 3)), "A", 0, "chr1", 0, "left", 60,
                             alt_seq="A", alt_bp=5000, cont_len=5000, cls="C1", site=1)]
        results, _ = rr.summarize(ends, sites)
        self.assertAlmostEqual(results[nclose]["j"], 2 * (10 / 30 * 2.0))
        self.assertEqual(results[nclose]["groups"], 2)

    def placed_end(self, end_id, member, end, chrom, bp, site, cls):
        return rr.MemberEnd(end_id, ("rep", (0, 1)), member, end, 0, chrom, bp, "left", 60,
                            alt_seq="A", alt_bp=5000, cont_len=5000, cls=cls, site=site)

    def alternatives(self, left_cn, right_cn, swap):
        # Two compound alternatives share the start junction class L and end in distinct classes R1, R2.
        rep, alt = ("rep", (0, 1)), ("alt", (2, 3))
        # 40 informative reads at CN_local 2: 20 reads per N.
        sites = [measured_site(0, 2.0, 40 - round(left_cn * 20), {"L": round(left_cn * 20)}),
                 measured_site(1, 2.0, 40 - round(right_cn * 20), {"R1": round(right_cn * 20)}),
                 measured_site(2, 2.0, 40 - round(right_cn * 20), {"R2": round(right_cn * 20)})]
        alt_ends = ("B", "A") if swap else ("A", "B")
        ends = [self.placed_end(0, rep, "A", "chr1", 100_000, 0, "L"),
                self.placed_end(1, rep, "B", "chr2", 500_000, 1, "R1"),
                self.placed_end(2, alt, alt_ends[0], "chr1", 100_000, 0, "L"),
                self.placed_end(3, alt, alt_ends[1], "chr2", 501_000, 2, "R2")]
        est = rr.summarize(ends, sites)[0][rep]
        return est

    def test_alternatives_sharing_one_junction_add_up_to_its_copies(self):
        est = self.alternatives(0.30, 0.15, swap=False)
        self.assertAlmostEqual(est["j"], 0.30)
        self.assertEqual(rr.prefit_decision(est), "")

    def test_shared_junction_bounds_the_alternatives_whatever_the_end_labels(self):
        # End labels of compressed members need not match the representative's.
        for swap in (False, True):
            self.assertAlmostEqual(self.alternatives(0.20, 0.15, swap)["j"], 0.20)

    def test_class_shared_across_opposite_ends_bounds_both_members_once(self):
        rep, alt = ("rep", (0, 1)), ("alt", (2, 3))
        sites = [measured_site(i, 2.0, 37, {c: 3}) for i, c in enumerate(("C", "D", "E"))]
        ends = [self.placed_end(0, rep, "A", "chr1", 100_000, 0, "C"),
                self.placed_end(1, rep, "B", "chr2", 500_000, 1, "D"),
                self.placed_end(2, alt, "A", "chr1", 101_000, 2, "E"),
                self.placed_end(3, alt, "B", "chr2", 500_500, 0, "C")]
        est = rr.summarize(ends, sites)[0][rep]
        self.assertAlmostEqual(est["j"], 0.15)
        self.assertEqual(rr.prefit_decision(est), "read_ratio_low")

    def test_tied_estimates_report_the_better_sampled_one(self):
        nclose = ("u", (0, 1))
        sites = [measured_site(0, 1.0, 10, {}), measured_site(1, 1.0, 30, {})]
        ends = [member_end(0, nclose, "A", 0, "C0"), member_end(1, nclose, "B", 1, "C1")]
        results, _ = rr.summarize(ends, sites)
        self.assertAlmostEqual(results[nclose]["r_eff"], 30.0)


class RescueTemplateTests(unittest.TestCase):
    def test_olc_owners_use_their_assembled_unitig(self):
        with tempfile.TemporaryDirectory() as tmp:
            rescue = Path(tmp) / "24_raw_rescue"
            rescue.mkdir()
            (rescue / "raw_candidates.json").write_text(
                '[{"owner": "raw_rescue_read_R1", "name": "read1"},'
                ' {"owner": "raw_rescue_olc_R2", "name": "Q001_utg1"}]')
            (rescue / "raw_realign.fa").write_text(">read1\nacgt\n")
            (rescue / "olc_realign.fa").write_text(">Q001_utg1\nggcc\n")
            templates = rr.load_rescue_templates(tmp)
        self.assertEqual(templates, {"raw_rescue_read_R1": "ACGT", "raw_rescue_olc_R2": "GGCC"})


class ReportWeightTests(unittest.TestCase):
    def test_source_aliases_carry_the_canonical_weight(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp)
            (prefix / "01_nclose_data.pkl").write_bytes(b"")
            os.utime(prefix / "01_nclose_data.pkl", (1, 1))
            (prefix / "nclose_report.tsv").write_text(
                "nclose_id\tnclose_cn\tkind\tnclose_key\n"
                "SKYPE.nclose.1\t1.5\tbnd\t1:0\n")
            (prefix / "nclose_sources.tsv").write_text(
                "source_nclose_key\tsource_nclose_id\tnclose_id\tcanonical_nclose_key\tis_alias\n"
                "0:1\tSKYPE.nclose.1\tSKYPE.nclose.1\t0:1\t0\n"
                "3:2\t.\tSKYPE.nclose.1\t0:1\t1\n")
            weights, per_stage = rr._report_weights(prefix)
        self.assertEqual(weights[(0, 1)], 1.5)
        self.assertEqual(weights[(2, 3)], 1.5)
        self.assertEqual(per_stage[(2, 3)], {"raw": 1.5})

    def test_scale_converts_every_reported_weight(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp)
            (prefix / "01_nclose_data.pkl").write_bytes(b"")
            os.utime(prefix / "01_nclose_data.pkl", (1, 1))
            (prefix / "nclose_report.tsv").write_text(
                "nclose_id\tnclose_cn\tkind\tnclose_key\n"
                "SKYPE.nclose.1\t1.5\tbnd\t1:0\n")
            (prefix / "nclose_report_cluster.tsv").write_text(
                "nclose_id\tnclose_cn\tkind\tnclose_key\n"
                "SKYPE.nclose.1\t3.0\tbnd\t1:0\n")
            (prefix / "nclose_sources.tsv").write_text(
                "source_nclose_key\tsource_nclose_id\tnclose_id\tcanonical_nclose_key\tis_alias\n"
                "3:2\t.\tSKYPE.nclose.1\t0:1\t1\n")
            weights, per_stage = rr._report_weights(prefix, 2 / 3)
        self.assertAlmostEqual(weights[(0, 1)], 2.0)
        self.assertAlmostEqual(per_stage[(0, 1)]["raw"], 1.0)
        self.assertAlmostEqual(per_stage[(2, 3)]["raw"], 1.0)


class ReportUnitTests(unittest.TestCase):
    """Stage 31's NClose CN (copy unit) is converted into median depth / 2."""

    DEPTHS = [10.0, 20.0, 30.0, 40.0, 50.0]   # median 30: median depth / 2 = 15

    def write_stat(self, prefix):
        import gzip
        path = prefix / "main.win.stat.gz"
        with gzip.open(path, "wt") as handle:
            handle.write("#chr\tst\tnd\tlength\tcovsite\ttotaldepth\tcov\tmeandepth\n")
            for i, depth in enumerate(self.DEPTHS):
                st = i * 100_000
                handle.write(f"chr1\t{st}\t{st + 100_000}\t100000\t100000\t0\t1\t{depth}\n")
            handle.write("chrM\t0\t16569\t16569\t16569\t0\t1\t9999.0\n")
        return path

    def test_baf_unit_is_converted_to_median_half(self):
        import depth_scale
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp)
            stat = self.write_stat(prefix)
            depth_scale.write(str(prefix), str(stat), self.DEPTHS, {}, {"N": 10.0})
            ratio, unit, source = rr.report_unit_ratio(prefix, stat)
        self.assertEqual((unit, source), (10.0, "baf"))
        self.assertAlmostEqual(ratio, 10.0 / 15.0)
        # A weight of 30 depth reported by stage 31 as 3.0 (in 1N = 10) is 2.0 here.
        self.assertAlmostEqual(3.0 * ratio, 30.0 / 15.0)

    def test_without_depth_scale_json_the_report_is_already_in_median_half(self):
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp)
            stat = self.write_stat(prefix)
            self.assertEqual(rr.report_unit_ratio(prefix, stat), (1.0, 15.0, "median"))

    def test_json_of_another_depth_table_falls_back_like_stage_31(self):
        import depth_scale
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp)
            stat = self.write_stat(prefix)
            depth_scale.write(str(prefix), str(stat), [1.0, 2.0, 3.0], {}, {"N": 10.0})
            with self.assertLogs(level="WARNING"):
                ratio, _, source = rr.report_unit_ratio(prefix, stat)
        self.assertEqual((ratio, source), (1.0, "median"))


class PosthocTests(unittest.TestCase):
    def test_poisson_lower_tail(self):
        self.assertAlmostEqual(rr.poisson_cdf(1, 7.5), math.exp(-7.5) * 8.5)

    def test_lower_tail_stays_accurate_at_large_means(self):
        # exp(-mu) underflows near mu = 745; the tail itself is far from zero here.
        self.assertAlmostEqual(rr.poisson_cdf(1200, 1260), 0.0464, places=3)
        est = dict(status="TESTABLE", j=6.0, r_eff=200.0, alt_reads=1200, exposure=200.0)
        self.assertEqual(rr.posthoc_decision(est, 6.3, 0.1)[0], "")

    def test_overweight_requires_effect_size_and_tail(self):
        est = dict(status="TESTABLE", j=0.06, r_eff=17.0, alt_reads=1, exposure=17.0)
        self.assertEqual(rr.posthoc_decision(est, 0.5)[0], "read_ratio_overweight")
        self.assertEqual(rr.posthoc_decision(est, 0.2)[0], "")
        self.assertEqual(rr.posthoc_decision(dict(est, j=0.9), 1.0)[0], "")

    def test_reads_per_copy_uncertainty_widens_the_tail(self):
        self.assertAlmostEqual(rr.marginal_poisson_cdf(3, 7.5, 0.0), rr.poisson_cdf(3, 7.5))
        self.assertGreater(rr.marginal_poisson_cdf(47, 74.0, 0.35), 50 * rr.poisson_cdf(47, 74.0))

    def test_sharp_tail_is_integrated_accurately(self):
        # Dense integration of the same model gives 0.010834: not below the 0.01 threshold.
        self.assertAlmostEqual(rr.marginal_poisson_cdf(99, 1000.0, 1.0), 0.010834, places=5)
        est = dict(status="TESTABLE", j=1.98, r_eff=50.0, alt_reads=99, exposure=50.0)
        self.assertEqual(rr.posthoc_decision(est, 20.0, 1.0)[0], "")

    def test_real_junction_with_some_excess_weight_is_kept(self):
        # Many ALT reads (J well above the display floor) and w within the scaling uncertainty.
        est = dict(status="TESTABLE", j=1.35, r_eff=35.3, alt_reads=47, exposure=35.3)
        self.assertEqual(rr.posthoc_decision(est, 2.1, 0.0)[0], "read_ratio_overweight")
        self.assertEqual(rr.posthoc_decision(est, 2.1, 0.35)[0], "")

    def test_weight_split_over_member_groups_is_tested_on_whole_reads(self):
        # Two compressed members with distinct junctions: 0/40 reads and 5/42 reads at CN 2.
        a, b = ("a", (0, 1)), ("b", (2, 3))
        sites = [measured_site(0, 2.0, 40, {}), measured_site(1, 2.0, 37, {"C1": 5})]
        ends = [rr.MemberEnd(0, a, a, "A", 0, "chr1", 0, "left", 60, alt_seq="A", alt_bp=5000,
                             cont_len=5000, cls="C0", site=0),
                rr.MemberEnd(1, a, b, "A", 0, "chr1", 0, "left", 60, alt_seq="A", alt_bp=5000,
                             cont_len=5000, cls="C1", site=1)]
        est = rr.summarize(ends, sites)[0][a]
        self.assertEqual((est["alt_reads"], est["exposure"]), (5, 20.0))
        self.assertEqual(rr.prefit_decision(est), "")
        # All of w in the less-sampled group predicts 12 reads; 5 are seen.
        self.assertEqual(rr.posthoc_decision(est, 0.6, 0.1)[0], "")

    def test_reads_seen_at_both_breakpoints_count_once(self):
        # The same 11 split reads are fetched at both ends of a simple junction.
        nclose = ("u", (0, 1))
        reads = [f"r{i}" for i in range(11)]
        sites = [measured_site(0, 2.0, 29, {"C0": 11}, names={"C0": reads}),
                 measured_site(1, 2.0, 29, {"C0": 11}, names={"C0": reads})]
        ends = [member_end(0, nclose, "A", 0, "C0"), member_end(1, nclose, "B", 1, "C0")]
        est = rr.summarize(ends, sites)[0][nclose]
        self.assertAlmostEqual(est["j"], 0.55)
        self.assertEqual((est["alt_reads"], est["exposure"]), (11, 20.0))
        # 11 molecules against 20 expected, not 22 observations against 40.
        self.assertEqual(rr.posthoc_decision(est, 1.0, 0.1)[0], "")

    def test_ends_disagreeing_on_reads_per_copy_use_the_lower(self):
        # H2009 chr13-chr22: the same 83 split reads at both ends; the chr22 short-arm end
        # holds 2.5x the line's reads per copy (26 vs 65), a calibration error at that end.
        nclose = ("u", (0, 1))
        reads = [f"r{i}" for i in range(83)]
        sites = [measured_site(0, 4.262, 28, {"C0": 83}, names={"C0": reads}),
                 measured_site(1, 3.723, 158, {"C0": 83}, names={"C0": reads})]
        ends = [member_end(0, nclose, "A", 0, "C0"), member_end(1, nclose, "B", 1, "C0")]
        est = rr.summarize(ends, sites)[0][nclose]
        self.assertEqual(est["alt_reads"], 83)
        self.assertAlmostEqual(est["exposure"], 111 / 4.262)
        # w 2.76 predicts 72 reads at the lower end, 179 at the higher; 83 are seen.
        self.assertEqual(rr.posthoc_decision(est, 2.76, 0.14)[0], "")

    def test_members_chained_by_classes_are_tested(self):
        # Members (L1,R1), (L1,R2), (L2,R2): no class is common to all three.
        nclose = ("u", (0, 1))
        left, right = [f"a{i}" for i in range(6)], [f"b{i}" for i in range(6)]
        sites = [measured_site(0, 2.0, 28, {"L1": 6, "L2": 6}, names={"L1": left, "L2": right}),
                 measured_site(1, 2.0, 28, {"R1": 6, "R2": 6}, names={"R1": left[:3] + right[:3],
                                                                      "R2": left[3:] + right[3:]})]
        ends = []
        for k, (lc, rc) in enumerate((("L1", "R1"), ("L1", "R2"), ("L2", "R2"))):
            member = ("u", (2 * k, 2 * k + 1))
            for end, site, cls in (("A", 0, lc), ("B", 1, rc)):
                ends.append(rr.MemberEnd(len(ends), nclose, member, end, 0, "chr1", 0, "left", 60,
                                         alt_seq="A", alt_bp=5000, cont_len=5000, cls=cls, site=site))
        est = rr.summarize(ends, sites)[0][nclose]
        self.assertAlmostEqual(est["j"], 0.6)
        self.assertEqual((est["status"], est["alt_reads"], est["exposure"]), ("TESTABLE", 12, 20.0))
        self.assertEqual(rr.posthoc_decision(est, 2.0, 0.1)[0], "read_ratio_overweight")

    def test_member_without_a_tested_end_of_its_own_is_not_tested_after_the_fit(self):
        # Member b shares a's allele but its own anchor is not uniquely placed.
        a, b = ("a", (0, 1)), ("b", (2, 3))
        sites = [measured_site(0, 2.0, 40, {})]
        ends = [rr.MemberEnd(0, a, a, "A", 0, "chr1", 0, "left", 60, alt_seq="A", alt_bp=5000,
                             cont_len=5000, cls="C0", site=0),
                rr.MemberEnd(1, a, b, "A", 0, "chr1", 0, "left", 0, alt_seq="A", alt_bp=5000,
                             cont_len=5000, cls="C0", site=0)]
        est = rr.summarize(ends, sites)[0][a]
        self.assertEqual(est["status"], "TESTABLE")
        self.assertTrue(math.isnan(est["exposure"]))
        self.assertEqual(rr.posthoc_decision(est, 2.0, 0.1)[0], "")

    def test_absent_junction_is_removed_despite_scaling_uncertainty(self):
        est = dict(status="TESTABLE", j=0.0, r_eff=12.4, alt_reads=0, exposure=12.4)
        self.assertEqual(rr.posthoc_decision(est, 0.74, 0.35)[0], "read_ratio_overweight")

    def test_untestable_is_never_removed(self):
        est = dict(status="NOT_TESTABLE", j=float("nan"), r_eff=float("nan"))
        self.assertEqual(rr.posthoc_decision(est, 3.0)[0], "")


if __name__ == "__main__":
    unittest.main()
