"""`ANS-9` step 1d (§9 item 410): ``_write_comparison`` keeps the weekly's Verdict section.

Render-only, on scratch copies of the committed ``metrics.json`` / ``COMPARISON.md``:

(K1)     ``D1_RAGGED_COLUMN`` equals the committed degree-1 values exactly — the
         top-level keys while ``metrics.json`` still carries the 2026-10-04 column,
         ``columns.d1_ragged`` once step 1d's degree-2 column has been written;
(R-keep) the tail from the first ``## Verdict`` line is byte-identical to the
         committed file's tail (heading count 1);
(N-tail) with the Verdict section removed from the scratch copy the output carries
         no ``## Verdict`` line: the tail is read, not hard-coded.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans9_coil_driven_gel_heating_10MHz"
SCRIPT = CASE / "09_coil_driven_gel_heating_10MHz.py"
MARK = b"## Verdict"


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ans9():
    return _load("ans9_writer_ut", SCRIPT)


def _count(data: bytes) -> int:
    return sum(1 for line in data.split(b"\n") if line.startswith(MARK))


def _tail(data: bytes) -> bytes:
    i = data.find(b"\n" + MARK)
    assert i >= 0
    return data[i + 1:]


def test_k1_d1_column_matches_committed(ans9):
    m = json.loads((CASE / "metrics.json").read_text(encoding="utf-8"))
    d1 = ans9.D1_RAGGED_COLUMN
    if "columns" in m:
        src = m["columns"]["d1_ragged"]
        print(f"(K1) source columns.d1_ragged, {len(src)} keys")
        assert src == d1
    else:
        src = m
        shared = [k for k in d1 if k in src]
        print(f"(K1) source top-level (pre-1d), {len(shared)} shared keys; em_element {src['em_element']}")
        assert src["em_element"] == "N1curl degree 1"
        assert len(shared) >= 20
        for k in shared:
            assert src[k] == d1[k], k


def test_r_keep_verdict_tail_byte_identical(ans9, tmp_path):
    shutil.copy(CASE / "COMPARISON.md", tmp_path / "COMPARISON.md")
    m = json.loads((CASE / "metrics.json").read_text(encoding="utf-8"))
    committed = (CASE / "COMPARISON.md").read_bytes()
    before = committed
    ans9._write_comparison(tmp_path, m)
    out = (tmp_path / "COMPARISON.md").read_bytes()
    t_in, t_out = _tail(committed), _tail(out)
    print(f"(R-keep) tail {len(t_in)} B in, {len(t_out)} B out; '## Verdict' count in {_count(committed)}, "
          f"out {_count(out)}")
    assert _count(committed) == 1
    assert _count(out) == 1
    assert t_out == t_in
    assert (CASE / "COMPARISON.md").read_bytes() == before  # the repo file is untouched
    # both columns present: primary and the carried degree-1 header
    head = out[: len(out) - len(t_out)].decode("utf-8")
    assert "Ours (degree 2, exact, h_p 0.0025)" in head
    assert "Ours (degree 1, ragged; 2026-10-04)" in head
    assert "| 4 | theta probe D at 1254 s, K |" in head


def test_n_tail_verdict_is_read_not_hardcoded(ans9, tmp_path):
    committed = (CASE / "COMPARISON.md").read_bytes()
    stripped = committed[: committed.find(b"\n" + MARK) + 1]
    assert _count(stripped) == 0
    (tmp_path / "COMPARISON.md").write_bytes(stripped)
    m = json.loads((CASE / "metrics.json").read_text(encoding="utf-8"))
    ans9._write_comparison(tmp_path, m)
    out = (tmp_path / "COMPARISON.md").read_bytes()
    print(f"(N-tail) '## Verdict' count after rendering a Verdict-less copy: {_count(out)}")
    assert _count(out) == 0
