"""New-telomere rules, on unitig layouts observed in U2OS and HCC1954 (hs1)."""
from __future__ import annotations

import unittest
from collections import defaultdict
from pathlib import Path

import nclose_preprocess as npp
from nclose_candidate import NCloseCandidate

PUBLIC = Path(__file__).resolve().parents[1] / "public_data"
CHR_LEN = npp.find_chr_len(str(PUBLIC / "chm13v2.0.fa.fai"))
TELO_DICT = defaultdict(list)
for _row in npp.import_telo_data(str(PUBLIC / "chm13v2.0_telomere.bed"), CHR_LEN):
    TELO_DICT[_row[0]].append(_row[1:])


def paf_rows(name, *pieces):
    """PAF-like rows (import_data layout) in query order: (chrom, start, end, dir, mapq)."""
    rows, query = [], 0
    for index, (chrom, start, end, direction, mapq) in enumerate(pieces):
        rows.append([name, 0, query, query + end - start, direction, chrom,
                     CHR_LEN[chrom], start, end, mapq, index])
        query += end - start
    for row in rows:
        row[1] = query
    return rows


def graph_rows(name, *pieces, telcon=()):
    """Stage-01 node rows for one contig; telcon maps piece index to a TELCON name."""
    telcon = dict(telcon)
    rows = []
    last = len(pieces) - 1
    for index, (chrom, start, end, direction, mapq) in enumerate(pieces):
        rows.append([name, 0, 0, end - start, direction, chrom, CHR_LEN[chrom], start, end,
                     mapq, 1, 0, last, "0", "0", telcon.get(index, "0"), "0", "0", "0",
                     direction, chrom, f"0.{index}"])
    return rows


def telomere_connections(rows):
    labels = npp.label_node(rows, TELO_DICT)
    _, _, info = npp.preprocess_telo([list(row) for row in rows], labels)
    return {(rows[idx][npp.CHR_NAM], rows[idx][npp.CHR_STR]): name for idx, name in info.items()}


class NewTelomereHostTests(unittest.TestCase):
    def test_telomere_capped_unique_host_is_kept(self):
        # U2OS utg024970l: telomere repeats + chr6 from 224,853 (reads: 10 telomeric clips).
        rows = paf_rows(
            "utg024970l",
            ("chr3", 201101722, 201105939, "-", 0),
            ("chr9", 0, 2270, "+", 60),
            ("chr4", 1, 3046, "+", 1),
            ("chr6", 224853, 281045, "+", 60),
        )
        self.assertEqual(telomere_connections(rows), {("chr6", 224853): "chr6f"})

    def test_multi_mapping_host_is_rejected(self):
        # U2OS utg089255l: chr2q14 interstitial telomeric sequence, MAPQ 0, no split reads.
        rows = paf_rows(
            "utg089255l",
            ("chr2", 114016003, 114017310, "+", 0),
            ("chr3", 201090946, 201104563, "+", 1),
            ("chr19", 61704425, 61705921, "+", 60),
            ("chr21", 45086593, 45087510, "+", 1),
        )
        self.assertNotIn(("chr2", 114016003), telomere_connections(rows))
        labels = npp.label_node(rows, TELO_DICT)
        self.assertEqual(npp.new_telomere_rejection_reason(rows, labels, 0, "chr2b"), "host_mapq")

    def test_sub_kilobase_host_is_rejected(self):
        # U2OS utg085665l: 330 bp of chr7 between chr16 and telomere repeats.
        rows = paf_rows(
            "utg085665l",
            ("chr16", 14591582, 14604979, "+", 60),
            ("chr7", 27251170, 27251500, "+", 60),
            ("chr21", 45088088, 45088786, "+", 0),
            ("chr9", 304, 3652, "-", 0),
        )
        self.assertEqual(telomere_connections(rows), {})

    def test_host_inserted_into_another_chromosome_end_is_rejected(self):
        # U2OS utg011541l: 1.5 kb of chrX inside the chr20p end (split reads on both
        # sides go to chr20p); chrX itself is not truncated.
        rows = paf_rows(
            "utg011541l",
            ("chr9", 14732, 132237, "-", 60),
            ("chr20", 3148, 81664, "-", 60),
            ("chrX", 72458078, 72459546, "-", 60),
            ("chr7", 1025, 2683, "-", 0),
        )
        self.assertEqual(telomere_connections(rows), {})

    def test_reference_end_connection_is_unaffected(self):
        rows = paf_rows(
            "ref_end",
            ("chr7", 0, 3000, "+", 0),
            ("chr7", 3000, 3800, "+", 0),
        )
        labels = npp.label_node(rows, TELO_DICT)
        self.assertIsNone(npp.new_telomere_rejection_reason(rows, labels, 1, "chr7f"))


class SubteloCutDirectionTests(unittest.TestCase):
    def cut(self, rows):
        return npp.subtelo_cut(
            rows, npp.label_node(rows, TELO_DICT), npp.label_subtelo_node(rows, TELO_DICT),
        )

    def test_fusion_toward_partner_telomere_is_not_cut(self):
        # HCC1954 utg001852l: chr8 then chr5 from 1,904 inward (69 split reads to chr5p).
        rows = graph_rows(
            "utg001852l",
            ("chr8", 107225450, 107249407, "-", 60),
            ("chr8", 106402703, 106410551, "+", 60),
            ("chr5", 1904, 39836, "+", 60),
        )
        self.assertEqual(self.cut(rows), [])

    def test_captured_terminal_segment_is_cut(self):
        rows = graph_rows(
            "capture",
            ("chr8", 107225450, 107249407, "-", 60),
            ("chr8", 106402703, 106410551, "+", 60),
            ("chr5", 1904, 39836, "-", 60),
        )
        cut = self.cut(rows)
        self.assertEqual([row[npp.CTG_TELCON] for row in cut], ["0", "chr8b"])


class SubtelomericOrientationTests(unittest.TestCase):
    def filter(self, rows, telo_contig=None):
        candidate = NCloseCandidate(rows[0][npp.CTG_NAM], (0, 1))
        kept, _ = npp.apply_subtelomeric_orientation_filter(
            [candidate], rows, telo_contig or {}, CHR_LEN,
        )
        return kept

    def test_internal_new_telomere_does_not_make_endpoint_terminal(self):
        # U2OS utg037835l: chrX 86.82 Mb -> chr8q end; a new chrX telomere lies 12.5 kb away.
        rows = graph_rows(
            "utg037835l",
            ("chrX", 86815582, 86828132, "+", 60),
            ("chr8", 146239312, 146256507, "-", 60),
        ) + graph_rows("utg070929l", ("chrX", 86815582, 86820878, "+", 60), telcon={0: "chrXf"})
        self.assertEqual(len(self.filter(rows)), 1)

    def test_reference_subtelomere_pair_is_still_rejected(self):
        # U2OS utg036495l: chr6p end <-> chr20q end, same telomere orientation.
        rows = graph_rows(
            "utg036495l",
            ("chr6", 6298, 40000, "-", 60),
            ("chr20", 66103776, 66150000, "+", 60),
        )
        self.assertEqual(self.filter(rows), [])


if __name__ == "__main__":
    unittest.main()
