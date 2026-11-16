"""`ANS-8` step 1 (§9 item 141): the private-mode comparison writer.

A synthetic AED JSON (one design carrying a 32×32 ``"S"`` on a shuffled port
order, every value synthetic) in a scratch copy of the case; its rows reach
the private table **through** ``ring_mode_rows`` — the same function our
side's stage R uses.  `ANS-7` item 68's four anchors:

1. every sentinel (AED's own C_mode, C_k2 and drive-(ii) share, as the writer
   formats them) lands in the AED row of the private table;
2. the tracked files (``metrics.json``, ``COMPARISON.md``) are byte-identical
   after the write and no sentinel reaches them;
3. ``FEM_EM_ANS_SKIP_PRIVATE=1`` writes nothing;
4. no ``aed_results/`` ⇒ nothing written.

Negative control (the leak): the private text appended to the tracked
``COMPARISON.md`` is caught by the same sentinel scan anchor 2 runs.

SPEC v2 (§9 item 161): the writer goes through ``ring_mode_rows_v2``.
(P5) the v2 sentinels land in the AED row, none in tracked files, and AED's 2b
runs at AED's own row-1 C_mode (bitwise), not ours; (P6, cross-wiring
control) AED's 2b recomputed at our synthetic C_mode differs by > 1e-6;
(P7) a pre-existing private file with a foreign first line is kept
byte-identical.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from fem_em_solver.ports.circuit import z_to_s
from fem_em_solver.ports.ring_modes import N_SLOTS, ring_mode_rows

REPO = Path(__file__).resolve().parents[2]
CASE = REPO / "examples" / "ansys_benchmarks" / "ans8_f_human_quadrature_mode_content_64MHz"
SCRIPT = CASE / "08_f_human_quadrature_mode_content.py"
TRACKED = ("metrics.json", "COMPARISON.md")
F = 64.0e6
W = 2 * np.pi * F


@pytest.fixture(scope="module")
def ans8():
    spec = importlib.util.spec_from_file_location("ans8_private_ut", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ans8_private_ut"] = mod
    spec.loader.exec_module(mod)
    return mod


def _z(c_mode, x2):
    j = np.arange(N_SLOTS)
    f = np.exp(-2j * np.pi * np.outer(j, j) / N_SLOTS) / np.sqrt(N_SLOTS)
    w = np.kron(np.eye(2), f)
    zh = np.zeros((32, 32), dtype=np.complex128)
    for k in range(N_SLOTS):
        kk = min(k, N_SLOTS - k)
        x = 1.0 / (W * c_mode) if kk == 1 else (x2 if kk == 2 else 40.0 + 9.0 * kk)
        amb, apb = (1.7 + 0.1 * kk) + 1j * x, 2.0 - 1j * (120.0 + 30.0 * kk)
        a, b = (apb + amb) / 2, (apb - amb) / 2
        zh[np.ix_([k, N_SLOTS + k], [k, N_SLOTS + k])] = [[a, b], [b, a]]
    return w.conj().T @ zh @ w


NAMES = [f"P{n}" for n in range(17, 49)]


OURS_C_MODE = 1.13e-10


def _ours_metrics(ans8):
    s = z_to_s(_z(OURS_C_MODE, 31.0), 50.0)
    pm = [(j, r) for r in ("bottom", "top") for j in range(N_SLOTS)]
    rows = ring_mode_rows(s, port_map=pm, freq_hz=F, c_f=ans8.SPEC_C_F,
                          feeds=[(0, "bottom"), (4, "bottom")], drives=ans8.SPEC_DRIVES)
    v2 = ans8.rows_v2(s, pm, NAMES, F)
    return {"cells": 1, "kappa": 0.0, "freq_hz": F, "rows": rows, "rows_v2": v2, "aed": None}


def _aed_json():
    s = z_to_s(_z(1.0987654e-10, 47.0), 50.0)
    perm = np.random.default_rng(11).permutation(32)
    s = s[np.ix_(perm, perm)]
    return {"_provenance": "synthetic",
            "design_zero_order": {"freq_hz": F, "port_names": [NAMES[p] for p in perm],
                                  "S": {"re": s.real.tolist(), "im": s.imag.tolist()}}}


def _sentinels(ans8):
    blk = _aed_json()["design_zero_order"]
    r = ans8.aed_rows(blk, F)["fixed_c"]
    return (f"{r['c_mode_f'] * 1e12:.6f}", f"{r['c_k2_f'] * 1e12:.6f}",
            f"{r['drives']['ccw']['rings']['bottom']['share_k1']:.6f}",
            f"{r['drives']['ccw']['rings']['bottom']['leak_k0']:.6f}")


def _v2_sentinels(ans8):
    """AED's v2 rows as the writer formats them: 2b ccw share / leak at AED's
    own C_mode, the 2c opposite share at the SPEC C, and 4a f₀ [MHz]."""
    v2 = ans8.aed_rows(_aed_json()["design_zero_order"], F)
    own = v2["own_c_mode"]["drives"]["ccw"]["rings"]["bottom"]
    opp = v2["controls"]["fixed_c"]["opposite"]["drives"]["opposite"]["rings"]["bottom"]
    return (f"{own['share_k1']:.6f}", f"{own['leak_k0']:.6f}", f"{opp['share_k1']:.6f}",
            f"{v2['f0_one_point_hz'] / 1e6:.6f}")


@pytest.fixture
def scratch(tmp_path, ans8, monkeypatch):
    case = tmp_path / "case"
    case.mkdir()
    ans8.write_outputs(_ours_metrics(ans8), case)
    aed = case / "aed_results" / "ans8_aed_results.json"
    aed.parent.mkdir()
    aed.write_text(json.dumps(_aed_json()))
    monkeypatch.setattr(ans8, "CASE_DIR", case)
    monkeypatch.setattr(ans8, "AED_RESULTS_PATH", aed)
    monkeypatch.delenv(ans8.SKIP_PRIVATE_ENV, raising=False)
    return case


def _leaks(case, sentinels):
    return [(n, s) for n in TRACKED for s in sentinels if s in (case / n).read_text()]


def test_private_writer_sentinels_and_tracked_bytes(ans8, scratch):
    sent = _sentinels(ans8)
    # The AED C_mode came through ring_mode_rows on the shuffled order.
    assert abs(float(sent[0]) - 109.87654) < 1e-6, sent
    before = {n: (scratch / n).read_bytes() for n in TRACKED}
    assert _leaks(scratch, sent) == [], "a sentinel already occurs in the tracked copy"
    note = ans8._write_private_if_available(json.loads((scratch / "metrics.json").read_text()))
    private = scratch / "COMPARISON_private.md"
    assert private.exists(), note
    aed_rows = [ln for ln in private.read_text().splitlines() if ln.startswith("| AED design_zero_order")]
    assert len(aed_rows) == 1
    for s in sent:
        assert s in aed_rows[0], f"sentinel {s} not in the AED row: {aed_rows[0]}"
    for n in TRACKED:
        assert (scratch / n).read_bytes() == before[n], f"{n} changed"
    assert _leaks(scratch, sent) == []


def test_leaking_variant_is_caught(ans8, scratch):
    sent = _sentinels(ans8)
    metrics = json.loads((scratch / "metrics.json").read_text())
    leaked = (scratch / "COMPARISON.md").read_text() + ans8._private_comparison_text(
        metrics, ans8._load_aed_results())
    (scratch / "COMPARISON.md").write_text(leaked)
    assert {s for _n, s in _leaks(scratch, sent)} == set(sent)


def test_skip_private_knob(ans8, scratch, monkeypatch):
    monkeypatch.setenv(ans8.SKIP_PRIVATE_ENV, "1")
    note = ans8._write_private_if_available(json.loads((scratch / "metrics.json").read_text()))
    assert "skipped" in note
    assert not (scratch / "COMPARISON_private.md").exists()


def test_no_aed_results_writes_nothing(ans8, scratch, monkeypatch):
    monkeypatch.setattr(ans8, "AED_RESULTS_PATH", scratch / "absent.json")
    note = ans8._write_private_if_available(json.loads((scratch / "metrics.json").read_text()))
    assert "no aed_results" in note
    assert not (scratch / "COMPARISON_private.md").exists()


# ---- SPEC v2 (§9 item 161): P5–P7 ----
def test_p5_v2_sentinels_at_aed_own_c_mode(ans8, scratch):
    v2 = ans8.aed_rows(_aed_json()["design_zero_order"], F)
    # AED 2b runs at AED's own row-1 C_mode (109.87654 pF synthetic), bitwise — not ours.
    assert v2["own_c_mode"]["c_f"] == v2["fixed_c"]["c_mode_f"]
    assert abs(v2["own_c_mode"]["c_f"] * 1e12 - 109.87654) < 1e-6
    assert abs(v2["own_c_mode"]["c_f"] / OURS_C_MODE - 1.0) > 1e-2
    sent = _v2_sentinels(ans8)
    print(f"P5 v2 sentinels {sent}; AED own c_f {v2['own_c_mode']['c_f']:.12e}")
    before = {n: (scratch / n).read_bytes() for n in TRACKED}
    assert _leaks(scratch, sent) == [], "a v2 sentinel already occurs in the tracked copy"
    note = ans8._write_private_if_available(json.loads((scratch / "metrics.json").read_text()))
    private = scratch / "COMPARISON_private.md"
    assert private.exists(), note
    row = [ln for ln in private.read_text().splitlines() if ln.startswith("| AED design_zero_order")]
    assert len(row) == 1
    for s in sent:
        assert s in row[0], f"v2 sentinel {s} not in the AED row"
    for n in TRACKED:
        assert (scratch / n).read_bytes() == before[n], f"{n} changed"
    assert _leaks(scratch, sent) == []


def test_p6_cross_wiring_control(ans8):
    blk = _aed_json()["design_zero_order"]
    v2 = ans8.aed_rows(blk, F)
    s = ans8._parse_complex_matrix(blk["S"])
    pm = ans8._ring_modes().ring_port_map_from_names(blk["port_names"])
    feeds = [pm[blk["port_names"].index(n)] for n in ans8.SPEC_FEED_NAMES]
    crossed = ring_mode_rows(s, port_map=pm, freq_hz=F, c_f=OURS_C_MODE, feeds=feeds,
                             drives=ans8.SPEC_DRIVES)
    a = v2["own_c_mode"]["drives"]["ccw"]["rings"]["bottom"]["share_k1"]
    b = crossed["drives"]["ccw"]["rings"]["bottom"]["share_k1"]
    print(f"P6 AED 2b ccw share at its own C {a:.9f} vs at our C_mode {b:.9f}: |d| {abs(a - b):.3e}")
    assert abs(a - b) > 1e-6


def test_p7_clobber_guard(ans8, scratch):
    private = scratch / "COMPARISON_private.md"
    foreign = "# Operator's own AED table\n| x | 1.234 |\n".encode()
    private.write_bytes(foreign)
    note = ans8._write_private_if_available(json.loads((scratch / "metrics.json").read_text()))
    print(f"P7 {note}")
    assert "kept" in note
    assert private.read_bytes() == foreign
