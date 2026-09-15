"""Shared minimap2 reference indexes for SKYPE and ACCtools-pipeline."""

import fcntl
import hashlib
import json
import logging
import os
from pathlib import Path
import subprocess
import tempfile

DEFAULT_CACHE_DIR = Path(__file__).resolve().parent.parent / "reference_indexes"


def file_signature(path):
    stat = os.stat(path)
    return {"path": os.path.abspath(path), "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns}


def ensure_reference_index(reference, preset, thread=1, cache_dir=None):
    """Share minimap2 indexes across samples and run directories."""
    reference = os.path.realpath(reference)
    manifest = {
        "reference": file_signature(reference),
        "preset": preset,
        "minimap2_version": subprocess.check_output(
            ["minimap2", "--version"], text=True
        ).strip(),
    }
    key = hashlib.sha256(json.dumps(manifest, sort_keys=True).encode()).hexdigest()[:20]
    cache_dir = os.path.abspath(cache_dir or DEFAULT_CACHE_DIR)
    os.makedirs(cache_dir, exist_ok=True)
    index = os.path.join(cache_dir, f"{os.path.basename(reference)}.{key}.mmi")
    # Concurrent runs share the build; publish only a completed index.
    with open(index + ".lock", "a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if os.path.isfile(index) and os.path.getsize(index) > 0:
            logging.info("Reusing reference index (%s): %s", preset, index)
            return index
        logging.info("Building reference index (%s): %s", preset, index)
        with tempfile.TemporaryDirectory(prefix="index-", dir=cache_dir) as tmp:
            temporary_index = os.path.join(tmp, "reference.mmi")
            with open(index + ".log", "w") as log:
                subprocess.run(
                    ["minimap2", "-x", preset, "-t", str(thread),
                     "-d", temporary_index, reference],
                    stdout=log, stderr=subprocess.STDOUT, check=True,
                )
            with open(index + ".json", "w") as handle:
                json.dump(manifest, handle, indent=2, sort_keys=True)
            os.replace(temporary_index, index)
        logging.info("Reference index ready: %s", index)
    return index
