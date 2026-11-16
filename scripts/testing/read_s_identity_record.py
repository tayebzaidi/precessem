#!/usr/bin/env python3
"""`OPS-77`: every recorded ``||S_study - S_direct||/||S_direct||`` reading against its band.

Stdlib only; run from the repo root.  Reads each log through a ``.gz``-safe route, the stdlib
twin of ``tests.log_records.read_log_text``, and prints log:line, value and multiple of the
band.  The band is the one the gate printed on the same line, so nothing here restates it.  A
log's ``[WF-16]`` line is the four-port fixture's, printed apart as a neighbour; the sweep-vs-point
reads of ``mri:38`` are the floor when both sides share a route.
"""

from __future__ import annotations

import glob
import gzip
import re
import sys
from pathlib import Path

LOGS = Path("docs/testing/logs")
PATTERNS = ("*_EX-167*.log", "*_WF-15-s10*.log", "*_WF-15-s11*.log", "*_OPS-77-mri37-*.log")
FLOOR_LOG = "20261102T060640Z_EX-170.log"
READ = re.compile(r"\|\|S_study - S_direct\|\|/\|\|S_direct\|\|\s+([0-9.]+e[+-][0-9]+)")
BAND = re.compile(r"IDENTITY_BAND\s+([0-9.]*e[+-]?[0-9]+)")
FLOOR = re.compile(r"\[\(S\)\] (\d+ MHz): sweep vs point rel ([0-9.]+e[+-][0-9]+)")


def _text(path):
    """A harness log's text, from the plain file or its ``.gz``."""
    if path.exists():
        return path.read_text(errors="replace")
    with gzip.open(f"{path}.gz", "rt", errors="replace") as fh:
        return fh.read()


def _names():
    found = []
    for pat in PATTERNS:
        found += glob.glob(str(LOGS / pat)) + glob.glob(str(LOGS / f"{pat}.gz"))
    return sorted({Path(p.removesuffix(".gz")).name for p in found})


def main():
    names = _names()
    print(f"RECORD {len(names)} logs matched {PATTERNS}")
    samples, carrying = [], 0
    for name in names:
        hit = False
        for n, line in enumerate(_text(LOGS / name).splitlines(), 1):
            m = READ.search(line)
            if not m:
                continue
            hit = True
            value = float(m.group(1))
            if "[WF-16]" in line:
                print(f"NEIGHBOUR {name}:{n} four-port (I) {value:.3e}, another fixture")
                continue
            b = BAND.search(line)
            if b is None:
                sys.exit(f"{name}:{n}: an (S) line with no printed band")
            band = float(b.group(1))
            samples.append(value)
            print(f"SAMPLE {name}:{n} {value:.3e} = {value / band:.4f} x band {band:g}, headroom "
                  f"{band / value:.1f} x")
        carrying += hit
    for n, line in enumerate(_text(LOGS / FLOOR_LOG).splitlines(), 1):
        m = FLOOR.search(line)
        if m:
            print(f"FLOOR {FLOOR_LOG}:{n} {m.group(1)}: {float(m.group(2)):.3e}")
    if not samples:
        sys.exit("no (S) sample found")
    print(f"RESULT logs_carrying={carrying} samples={len(samples)} min={min(samples):.3e} "
          f"max={max(samples):.3e} over_1e-13={sum(v > 1e-13 for v in samples)}")


if __name__ == "__main__":
    main()
