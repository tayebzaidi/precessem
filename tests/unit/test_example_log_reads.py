"""`OPS-75`: the two example scripts read a dated log in either form.

`OPS-73` taught the tests to read a log retention has gzipped; these two reads
are the route it left, running the example, and an example cannot import
`tests/log_records.py`. Fixtures are short synthetic text under ``tmp_path``.

(H′) on a ``.gz``-only copy, ``node_log_readings`` returns the plain copy's dict
     and ``parse_h_table_logs`` its record, under the **listed** name (the
     numbers file records it, `OPS-73`'s (X) red); each script's ``_read_log``
     equals ``read_log_text`` on that path.
(A′) with neither file present, each raises `FileNotFoundError`.
Control: the read as it stood on `main`, ``open(path)``, raises on a ``.gz``-only
path, so (H′) is not true by construction.
"""

from __future__ import annotations

import gzip
import importlib.util
import json
import sys
from pathlib import Path

import pytest

from tests.log_records import read_log_text

REPO = Path(__file__).resolve().parents[2]
ANS4 = (REPO / "examples" / "ansys_benchmarks" / "ans4_birdcage_four_port_10_64_128MHz"
        / "04_birdcage_four_port_10_64_128MHz.py")
ANS8 = (REPO / "examples" / "ansys_benchmarks" / "ans8_f_human_quadrature_mode_content_64MHz"
        / "08_f_human_quadrature_mode_content.py")
NODE = "h9p-node-1.0000-d2.json"
NODE_TEXT = (
    "[WF-10 step2] [ccw] w = (1.000+0.000j, 0.000-1.000j); (P_w) P_src 5.000000000e+00 "
    "sumP_sheet 6.000000000e+00 P_phantom 7.000000000e+00 P_rest 0.000000000e+00 "
    "P_surf 8.000000000e+00 W\n"
    "[WF-10 step2] [ccw] (Lin) currents rel 1.000e-16 (band); Y_in current route vs "
    "_active_y_in 1.184e-15 (band)\n"
    f"[WF-10 step2] wrote /workspace/logs/wf10/{NODE}\n"
    "tail ≈ with a non-ASCII byte\n"
)
_S_ROW = "S_{k}k = +1.0e-01+2.0e-01j +3.0e-01-4.0e-01j +5.0e-01+6.0e-01j -7.0e-01+8.0e-01j"
H_TEXT = ("[ANS-4 step5] degree 2, -n 8, cells 1234\n"
          "[ANS-4 step5] 64 MHz kappa-off S (10 digits):\n"
          + "".join(_S_ROW.format(k=k) + "\n" for k in range(1, 5)))


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def ans4():
    return _load("ops75_ans4_ut", ANS4)


@pytest.fixture(scope="module")
def ans8():
    return _load("ops75_ans8_ut", ANS8)


def _pair(tmp_path: Path, text: str):
    """The same name as a plain copy and a `.gz`-only copy, in two directories."""
    plain = tmp_path / "plain" / "20261102T000000Z_X.log"
    gz_only = tmp_path / "gz" / plain.name
    plain.parent.mkdir()
    gz_only.parent.mkdir()
    plain.write_text(text)
    with gzip.open(f"{gz_only}.gz", "wb") as fh:
        fh.write(text.encode())
    return plain, gz_only


def test_h_prime_node_log_readings(ans8, tmp_path):
    plain, gz_only = _pair(tmp_path, NODE_TEXT)
    got, want = ans8.node_log_readings(gz_only, NODE), ans8.node_log_readings(plain, NODE)
    print(f"\n[OPS-75] (H') node readings {got}")
    assert want["p_w"]["p_src"] == 5.0 and len(want["lin"]) == 1
    assert got == want
    assert ans8._read_log(gz_only) == read_log_text(gz_only) == NODE_TEXT


def test_h_prime_h_table_keeps_the_listed_name(ans4, tmp_path):
    plain, gz_only = _pair(tmp_path, H_TEXT)
    got = ans4.parse_h_table_logs([("r", gz_only)])["rungs"][0]
    want = ans4.parse_h_table_logs([("r", plain)])["rungs"][0]
    print(f"\n[OPS-75] (H') h rung log {got['log']!r}, cells {got['cells']}")
    assert want["cells"] == 1234 and len(want["freqs"]["64"]["rows"]) == 4
    assert got == want
    assert got["log"] == plain.name  # never the resolved `.gz` name
    assert ans4._read_log(gz_only) == read_log_text(gz_only) == H_TEXT


def test_h_table_reads_a_path_listed_as_gz(ans4, tmp_path):
    """`H_TABLE_LOGS`' ×1 rung names its `.gz` outright; that route stays."""
    plain, gz_only = _pair(tmp_path, H_TEXT)
    listed = Path(f"{gz_only}.gz")
    got = ans4.parse_h_table_logs([("r", listed)])["rungs"][0]
    assert got["log"] == listed.name and got["cells"] == 1234
    assert got["freqs"] == ans4.parse_h_table_logs([("r", plain)])["rungs"][0]["freqs"]


def test_h_table_logs_read_directly_equal_the_numbers_file(ans4):
    """The disclosed route on the real logs, no scratch copy: two of the three are
    listed plain and exist only as `.gz`, and the numbers file records those names."""
    paths = [(r, REPO / p) for r, p in ans4.H_TABLE_LOGS]
    on_disk = [p.exists() for _, p in paths]
    print(f"\n[OPS-75] H_TABLE_LOGS listed name on disk: {on_disk}")
    regen = json.dumps(ans4.parse_h_table_logs(paths), indent=1, sort_keys=True) + "\n"
    assert regen == ans4.H_TABLE_JSON.read_text()


def test_a_prime_both_absent_refused(ans4, ans8, tmp_path):
    absent = tmp_path / "absent.log"
    with pytest.raises(FileNotFoundError):
        ans8.node_log_readings(absent, NODE)
    with pytest.raises(FileNotFoundError):
        ans4.parse_h_table_logs([("r", absent)])


def test_nc_bare_open_raises_on_gz_only(tmp_path):
    _, gz_only = _pair(tmp_path, NODE_TEXT)
    with pytest.raises(FileNotFoundError):
        open(gz_only)
