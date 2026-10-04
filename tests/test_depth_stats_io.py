import gzip
import os
from pathlib import Path
import stat
import sys
import tempfile
import threading
import unittest
from unittest import mock

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from depth_stats_io import atomic_gzip_text


class AtomicDepthTableTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.folder = Path(self.directory.name)
        self.output = self.folder / "sample_normalized.win.stat.gz"
        self.output.write_bytes(gzip.compress(b"previous complete table\n"))

    def read(self):
        with gzip.open(self.output, "rt") as handle:
            return handle.read()

    def test_csv_payload_matches_existing_pandas_serialization(self):
        frame = pd.DataFrame({"chr": ["chr1", "chrX", "chrY"], "start": [1, 100001, 200001],
                              "depth": [0.0, 51.23456789, -0.0001]})
        legacy = self.folder / "legacy.gz"
        frame.to_csv(legacy, sep="\t", index=False, header=False, compression="gzip")
        with atomic_gzip_text(self.output) as output:
            frame.to_csv(output, sep="\t", index=False, header=False)
        self.assertEqual(gzip.decompress(legacy.read_bytes()), gzip.decompress(self.output.read_bytes()))

    def test_failed_writer_preserves_previous_bytes(self):
        before = self.output.read_bytes()
        with self.assertRaisesRegex(ValueError, "synthetic interrupted write"):
            with atomic_gzip_text(self.output) as output:
                output.write("partial new table\n")
                output.flush()
                raise ValueError("synthetic interrupted write")
        self.assertEqual(before, self.output.read_bytes())
        self.assertEqual([], list(self.folder.glob(".*.tmp")))

    def test_reader_and_second_writer_never_see_first_partial_output(self):
        staged, release = threading.Event(), threading.Event()
        errors = []

        def writer():
            try:
                with atomic_gzip_text(self.output) as output:
                    output.write("first writer beginning\n")
                    output.flush()
                    staged.set()
                    if not release.wait(10):
                        raise RuntimeError("Test release was not signaled")
                    output.write("first writer ending\n")
            except BaseException as exc:
                errors.append(exc)

        thread = threading.Thread(target=writer)
        thread.start()
        try:
            self.assertTrue(staged.wait(10))
            self.assertEqual("previous complete table\n", self.read())
            with atomic_gzip_text(self.output) as output:
                output.write("second complete table\n")
            self.assertEqual("second complete table\n", self.read())
        finally:
            release.set()
            thread.join(10)
        self.assertFalse(thread.is_alive())
        self.assertEqual([], errors)
        self.assertEqual("first writer beginning\nfirst writer ending\n", self.read())

    def test_new_destination_is_absent_until_complete(self):
        self.output.unlink()
        with atomic_gzip_text(self.output) as output:
            output.write("first completed table\n")
            output.flush()
            self.assertFalse(self.output.exists())
        self.assertEqual("first completed table\n", self.read())

    def test_existing_symlink_and_permissions_are_preserved(self):
        target = self.folder / "actual.gz"
        self.output.rename(target)
        target.chmod(0o640)
        self.output.symlink_to(target.name)
        with atomic_gzip_text(self.output) as output:
            output.write("new complete table\n")
        self.assertTrue(self.output.is_symlink())
        self.assertEqual(0o640, stat.S_IMODE(target.stat().st_mode))
        self.assertEqual("new complete table\n", self.read())

    def test_failed_publication_preserves_old_output(self):
        before = self.output.read_bytes()
        with mock.patch("depth_stats_io.os.replace", side_effect=OSError("publication failed")):
            with self.assertRaisesRegex(OSError, "publication failed"):
                with atomic_gzip_text(self.output) as output:
                    output.write("complete but unpublished table\n")
        self.assertEqual(before, self.output.read_bytes())
        self.assertEqual([], list(self.folder.glob(".*.tmp")))


if __name__ == "__main__":
    unittest.main()
