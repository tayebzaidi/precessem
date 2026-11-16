"""Read a dated harness log whichever form retention has left it in.

The weekly sweep gzips a cited log at 7 days (`docs/testing/retention-policy.md`
§1), so a name that was plain when its test landed becomes ``<name>.gz`` under
it. `OPS-73`: a test that reads a record goes through here, and never
decompresses one on disk.
"""

from __future__ import annotations

import glob
import gzip
from pathlib import Path


def resolve_log_path(path) -> Path:
    """The file on disk: ``path`` if it exists, else ``path`` + ``.gz``."""
    plain = Path(path)
    gz = Path(f"{plain}.gz")
    if plain.exists():
        return plain
    if gz.exists():
        return gz
    raise FileNotFoundError(f"neither {plain} nor {gz} exists")


def read_log_text(path, errors="replace") -> str:
    """A harness log's text, from the plain file or its ``.gz``."""
    found = resolve_log_path(path)
    if found.suffix == ".gz":
        with gzip.open(found, "rt", errors=errors) as fh:
            return fh.read()
    return found.read_text(errors=errors)


def glob_log_paths(pattern) -> list:
    """Sorted plain names of the logs matching ``pattern``, in either form."""
    pattern = str(pattern)
    found = glob.glob(pattern) + glob.glob(f"{pattern}.gz")
    return sorted({p.removesuffix(".gz") for p in found})
