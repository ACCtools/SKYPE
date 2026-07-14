from __future__ import annotations

import collections
import io
import unittest

import vcfpy

from skype_vcf_writer import write_vcf_record_with_fallback


class SkypeVcfWriterTests(unittest.TestCase):
    def test_incompatible_info_value_falls_back_to_raw_record(self) -> None:
        raw_record = (
            "chr1\t10\trecord1\tA\tT\t60\tPASS\t"
            "LINKED_BY;SKYPE_CN=9;SKYPE_STATUS=OLD\tGT\t0/1"
        )
        input_vcf = io.StringIO(
            "##fileformat=VCFv4.3\n"
            '##INFO=<ID=LINKED_BY,Number=.,Type=String,Description="links">\n'
            '##INFO=<ID=SKYPE_CN,Number=1,Type=Float,Description="copy number">\n'
            '##INFO=<ID=SKYPE_STATUS,Number=1,Type=String,Description="status">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
            f"{raw_record}\n"
        )
        reader = vcfpy.Reader.from_stream(input_vcf)
        record = next(reader)
        output = io.StringIO()
        writer = vcfpy.Writer.from_stream(output, reader.header)
        updates = collections.OrderedDict([
            ("SKYPE_CN", "0"),
            ("SKYPE_STATUS", "TEST"),
        ])

        error = write_vcf_record_with_fallback(
            writer,
            record,
            raw_record,
            updates,
        )

        self.assertIsInstance(error, ValueError)
        output_record = output.getvalue().splitlines()[-1]
        self.assertEqual(
            output_record,
            "chr1\t10\trecord1\tA\tT\t60\tPASS\t"
            "LINKED_BY;SKYPE_CN=0;SKYPE_STATUS=TEST\tGT\t0/1",
        )

    def test_compatible_record_uses_vcfpy_writer(self) -> None:
        raw_record = "chr1\t10\trecord1\tA\tT\t60\tPASS\tCOUNT=2"
        input_vcf = io.StringIO(
            "##fileformat=VCFv4.3\n"
            '##INFO=<ID=COUNT,Number=1,Type=Integer,Description="count">\n'
            '##INFO=<ID=SKYPE_CN,Number=1,Type=Float,Description="copy number">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            f"{raw_record}\n"
        )
        reader = vcfpy.Reader.from_stream(input_vcf)
        record = next(reader)
        output = io.StringIO()
        writer = vcfpy.Writer.from_stream(output, reader.header)

        error = write_vcf_record_with_fallback(
            writer,
            record,
            raw_record,
            collections.OrderedDict([("SKYPE_CN", "1.5")]),
        )

        self.assertIsNone(error)
        self.assertEqual(
            output.getvalue().splitlines()[-1],
            "chr1\t10\trecord1\tA\tT\t60\tPASS\tCOUNT=2;SKYPE_CN=1.5",
        )


if __name__ == "__main__":
    unittest.main()
