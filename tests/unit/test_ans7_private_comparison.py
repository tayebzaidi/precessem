"""`ANS-7` v2 (§9 item 68, chain T10): the private-mode comparison writer.

Synthetic AED JSON on the key layout of ``aed/ans7_hfss_pyaedt.py:44–61``
(names only; every value is a sentinel), in a scratch copy of the case.
Item 60's four anchors (`tests/unit/test_ans6_private_comparison.py`):

1. every sentinel lands in the private table, in its source's row;
2. the tracked files (``metrics.json``, ``COMPARISON.md``) are byte-identical
   after the write, and no sentinel reaches them;
3. ``FEM_EM_ANS_SKIP_PRIVATE=1`` writes nothing;
4. no ``aed_results/`` ⇒ nothing written.

Negative control: a *leaking* variant (the private text appended to the
tracked ``COMPARISON.md``) is caught by the same sentinel scan anchor 2 runs.
"""

from __future__ import annotations

import importlib.util
import json
import shutil
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans7_tuned_copper_birdcage_driven_sweep"
SCRIPT = CASE / "07_tuned_copper_birdcage_driven_sweep.py"
TRACKED = ("metrics.json", "COMPARISON.md")

# Sentinels: values that occur nowhere in our metrics; formatted as the writer does.
F0 = 61.2345678e6          # -> "61.234568"
R_IN = 4.4444321e3         # -> "4.444432e+03"
Q = 777.7777               # -> "777.7777"
K_SHARE = 0.123456         # -> "0.123456"
SENTINELS = ("61.234568", "4.444432e+03", "777.7777", "0.123456")


@pytest.fixture(scope="module")
def ans7():
    spec = importlib.util.spec_from_file_location("ans7_private_ut", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ans7_private_ut"] = mod
    spec.loader.exec_module(mod)
    return mod


def _synthetic_aed() -> dict:
    rows = {"f0_hz": F0, "r_in_ohm": R_IN, "q_loaded": Q, "dbdf_s_per_hz": 1.0,
            "bracket_hz": [F0 - 5e4, F0 + 5e4], "crossings_up_hz": [F0],
            "crossings_down_hz": []}
    return {
        "_provenance": "synthetic", "_meta": {}, "_reading": "",
        "design_first_order": {
            "FineWindow": {"freq_hz": [F0], "Zin": [{"re": 1.0, "im": 0.0}],
                           "Z4": [], "S4": [], "zin_z_vs_s_max_rel": 0.0,
                           "rows": rows, "k1_share_at_f0": K_SHARE},
            "elements": 0, "passes": 0, "delta_s": 0.0, "solve_time": 0.0,
            "basis_order": "first_order", "model": "synthetic",
        },
        "_route_comparison": {"cu_first_order|FineWindow": {
            "df0_rel": 0.0, "dq_rel": 0.0, "dr_in_rel": 0.0, "max_dzin_rel": 0.0}},
        "_column_difference": {"first_order|B|FineWindow": {
            "df0_hz": 0.0, "dq": 0.0, "dr_in_ohm": 0.0}},
    }


@pytest.fixture
def scratch(tmp_path, ans7, monkeypatch):
    case = tmp_path / "case"
    case.mkdir()
    for name in TRACKED:
        shutil.copy(CASE / name, case / name)
    aed = case / "aed_results" / "ans7_aed_results.json"
    aed.parent.mkdir()
    aed.write_text(json.dumps(_synthetic_aed()))
    monkeypatch.setattr(ans7, "CASE_DIR", case)
    monkeypatch.setattr(ans7, "AED_RESULTS_PATH", aed)
    monkeypatch.delenv(ans7.SKIP_PRIVATE_ENV, raising=False)
    return case


def _leaks(case: Path) -> list:
    return [(n, s) for n in TRACKED for s in SENTINELS if s in (case / n).read_text()]


def test_private_writer_sentinels_and_tracked_bytes(ans7, scratch):
    before = {n: (scratch / n).read_bytes() for n in TRACKED}
    assert _leaks(scratch) == [], "a sentinel already occurs in the tracked copy"
    metrics = json.loads((scratch / "metrics.json").read_text())
    note = ans7._write_private_if_available(metrics)
    private = scratch / "COMPARISON_private.md"
    assert private.exists(), note
    text = private.read_text()
    # Anchor 1: every sentinel in the AED row of the private table.
    aed_rows = [ln for ln in text.splitlines() if ln.startswith("| AED design_first_order")]
    assert len(aed_rows) == 1
    for s in SENTINELS:
        assert s in aed_rows[0], f"sentinel {s} not in the AED row: {aed_rows[0]}"
    # Ours rows present for both routes, both columns (Route A written by item 68).
    for route in ("A", "B"):
        for col in ("cu", "pec"):
            assert f"| ours {col} (route {route}) |" in text
    # Anchor 2: tracked files byte-identical, no sentinel in them.
    for n in TRACKED:
        assert (scratch / n).read_bytes() == before[n], f"{n} changed"
    assert _leaks(scratch) == []


def test_leaking_variant_is_caught(ans7, scratch):
    """Negative control: the private text written into the tracked COMPARISON.md."""
    metrics = json.loads((scratch / "metrics.json").read_text())
    aed = ans7._load_aed_results()
    leaked = (scratch / "COMPARISON.md").read_text() + ans7._private_comparison_text(metrics, aed)
    (scratch / "COMPARISON.md").write_text(leaked)
    found = _leaks(scratch)
    assert {s for _n, s in found} == set(SENTINELS), found


def test_skip_private_knob(ans7, scratch, monkeypatch):
    monkeypatch.setenv(ans7.SKIP_PRIVATE_ENV, "1")
    note = ans7._write_private_if_available(json.loads((scratch / "metrics.json").read_text()))
    assert "skipped" in note
    assert not (scratch / "COMPARISON_private.md").exists()


def test_no_aed_results_writes_nothing(ans7, scratch, monkeypatch):
    monkeypatch.setattr(ans7, "AED_RESULTS_PATH", scratch / "absent.json")
    note = ans7._write_private_if_available(json.loads((scratch / "metrics.json").read_text()))
    assert "no aed_results" in note
    assert not (scratch / "COMPARISON_private.md").exists()
