"""Publish complete compressed depth tables without exposing partial writers."""
from contextlib import contextmanager
import gzip
import io
import os
from pathlib import Path
import stat
import uuid


@contextmanager
def atomic_gzip_text(destination):
    """Yield a UTF-8 gzip writer, then replace the destination after close.

    Readers with the old file open keep that complete version. Concurrent
    writers each publish a complete file; this does not choose which inputs
    a shared output name should represent or attest their generation.
    Existing destination symlinks and target permission bits are preserved.
    """
    target = Path(destination).resolve()
    temporary = target.with_name(f".{target.name}.{uuid.uuid4().hex}.tmp")
    created = False
    try:
        fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
        created = True
        with os.fdopen(fd, "wb") as raw:
            try:
                mode = stat.S_IMODE(target.stat().st_mode)
            except FileNotFoundError:
                mode = None
            if mode is not None:
                os.fchmod(raw.fileno(), mode)
            # Name the final output in the gzip header, not the random staging
            # path. Pandas' ordinary gzip output uses the same default codec.
            with gzip.GzipFile(filename=target.name, mode="wb", fileobj=raw) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text:
                    yield text
            raw.flush()
            os.fsync(raw.fileno())
        os.replace(temporary, target)
    finally:
        if created:
            temporary.unlink(missing_ok=True)
