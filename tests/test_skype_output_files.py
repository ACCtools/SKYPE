from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from skype_output_files import (
    discover_ecdna_depth_inputs,
    discover_jump_depth_inputs,
)


class SkypeOutputFileTests(unittest.TestCase):
    def test_jump_discovery_supports_ordinary_and_type2_events(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            jump = root / "front_jump"
            type2 = root / "type2_ins"
            jump.mkdir()
            type2.mkdir()

            for name in (
                "1.win.stat.gz",
                "1_base.win.stat.gz",
                "2_type2_merge_7.win.stat.gz",
                "2_base.win.stat.gz",
            ):
                (jump / name).touch()
            (type2 / "7.win.stat.gz").touch()

            records = discover_jump_depth_inputs(str(jump), str(type2))

            self.assertEqual([record[0] for record in records], [1, 2])
            self.assertEqual([record[3] for record in records], [-1, 7])
            self.assertIsNone(records[0][4])
            self.assertEqual(records[1][4], str(type2 / "7.win.stat.gz"))

    def test_jump_discovery_rejects_mixed_run_primary_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            jump = root / "front_jump"
            type2 = root / "type2_ins"
            jump.mkdir()
            type2.mkdir()

            for name in (
                "1.win.stat.gz",
                "1_type2_merge_7.win.stat.gz",
                "1_base.win.stat.gz",
            ):
                (jump / name).touch()
            (type2 / "7.win.stat.gz").touch()

            with self.assertRaisesRegex(ValueError, "different runs were mixed"):
                discover_jump_depth_inputs(str(jump), str(type2))

    def test_ecdna_discovery_rejects_non_contiguous_indices(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            ecdna = Path(temporary)
            (ecdna / "1.win.stat.gz").touch()
            (ecdna / "3.win.stat.gz").touch()

            with self.assertRaisesRegex(ValueError, "Non-contiguous"):
                discover_ecdna_depth_inputs(str(ecdna))


if __name__ == "__main__":
    unittest.main()
