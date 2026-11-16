"""`PORT-23` step 2 -- scikit-rf reproduces `ans:7`'s Route B (§9 item 155).

On the 51 stored grid 4×4 S of ``ans:7`` (``logs/ans7/B1.json`` 62.0–64.4 MHz,
``B2.json`` 64.5–67.0 MHz; Cu and PEC columns; gitignored, written 2026-09-29),
the circuit layer is run twice: once through the raw numpy Route B of the
``ans:7`` script (``z_from_s`` -> ``route_b_z``), once through `scikit-rf`
(``network_from_records`` -> ``terminate`` with P2–P4 in ``Z_C = 1/(jωC_F)`` ->
``input_impedance``).  C₁ is added identically (``add_c1``) on both sides and
``parallel_rows`` reduces each to f₀ / R_in / Q_loaded.

Gates (asserted):
  (P0) the stored S, through the script's own S route (``route_b_s``, the route
       that wrote the stored ``Z_in``), reproduces the stored ``Z_in`` to 1e-12
       relative and ``DENSE_ROWS`` (metrics.json, imported) within ``BANDS``;
  (R1) skrf ``Z_bare`` vs ``route_b_z(z_from_s(S))[1]`` at all 51 points, both
       columns, to ``S_VS_Z_ROUTE_RTOL``;
  (R2) rows from both sides agree to 1e-9 relative, same ``bracket_hz``;
  (R3) the skrf-side rows match ``DENSE_ROWS`` within ``BANDS``.
Negative controls (*predicted*, printed only; rule (e)): (N1) the 50 Ω S labelled
z0 = 75 Ω; (N2) a port map swapping network ports P1 ↔ P3 (drive moved to the
opposite leg).  Either reading inside 1e-9 means the gate is blind.

Only pure functions of the ``ans:7`` script are used -- never ``_stage_*``,
``_assemble``, ``_write_outputs`` or ``main`` (they rewrite tracked files).
Stored S, numpy and skrf only; no field solve.  Scope: identities on stored
matrices -- no absolute Z_in / f₀ / Q claim (``ANS-7``, §2.2).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

skrf = pytest.importorskip("skrf")

from fem_em_solver.ports.circuit_skrf import (  # noqa: E402
    input_impedance,
    network_from_records,
    terminate,
)
from tests.validation.test_ans7_sparse_route_b import BANDS, DENSE_ROWS  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
SCRIPT = (REPO / "examples" / "ansys_benchmarks" / "ans7_tuned_copper_birdcage_driven_sweep"
          / "07_tuned_copper_birdcage_driven_sweep.py")
STAGE_DIR = REPO / "logs" / "ans7"
COLUMNS = ("cu", "pec")
ROWS = ("f0_hz", "r_in_ohm", "q_loaded")
Z0 = 50.0

STORED_Z_IN_RTOL = 1.0e-12   # (P0) reproduction of the stored Z_in (brief, §9 item 155)
ROW_AGREE_RTOL = 1.0e-9      # (R2) skrf rows vs raw rows (brief)
BLIND_FLOOR = 1.0e-9         # (N1)/(N2): a reading inside this means the gate is blind
N1_Z0 = 75.0
N2_PORT_MAP = (2, 1, 0, 3)   # logical P1 (drive) wired to network P3 and vice versa


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_ANS7 = _load(SCRIPT, "ans7_for_port23")
z_from_s = _ANS7.z_from_s
route_b_z = _ANS7.route_b_z
route_b_s = _ANS7.route_b_s
add_c1 = _ANS7.add_c1
parallel_rows = _ANS7.parallel_rows
C_F = _ANS7.C_F
S_VS_Z_ROUTE_RTOL = _ANS7.S_VS_Z_ROUTE_RTOL
_unmat = _ANS7._unmat


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _rel(a, b):
    return abs(a - b) / abs(b)


def _row_rel(rows, ref):
    return {k: (_rel(rows[k], ref[k]) if k in rows else np.inf) for k in ROWS}


@pytest.fixture(scope="module")
def data():
    pts = []
    for name in ("B1.json", "B2.json"):
        p = STAGE_DIR / name
        if not p.exists():
            pytest.fail(f"BLOCKED: {p} absent -- unblock by an ans:7 B1/B2 re-run (§9 item 155)")
        doc = json.loads(p.read_text())
        _say(f"[port23] {name}: {len(doc['points'])} points, mpi_ranks {doc.get('mpi_ranks')}")
        pts.extend(doc["points"])
    pts.sort(key=lambda p: p["f_hz"])
    f_hz = np.array([p["f_hz"] for p in pts], dtype=float)
    assert len(pts) == 51 and np.all(np.diff(f_hz) > 0)
    out = {"f_hz": f_hz}
    for col in COLUMNS:
        out[col] = {
            "s": np.array([_unmat(p[col]["s"]) for p in pts]),
            "z_in_stored": np.array([complex(p[col]["Z_in"]["re"], p[col]["Z_in"]["im"])
                                     for p in pts]),
        }
    _say(f"[port23] skrf {skrf.__version__}; C_F {C_F:.6e} F; z0 {Z0} ohm; "
         f"S_VS_Z_ROUTE_RTOL {S_VS_Z_ROUTE_RTOL:.0e}")
    return out


def _skrf_z_bare(s, f_hz, z0=Z0, port_map=None):
    net = network_from_records(s, z0, f_hz)
    zc = 1.0 / (1j * 2.0 * np.pi * f_hz * C_F)
    s11 = terminate(net, {1: zc, 2: zc, 3: zc}, port_map=port_map)
    assert s11.shape == (len(f_hz), 1, 1)
    return np.array([input_impedance(v, z0) for v in s11[:, 0, 0]])


def _rows_from_bare(f_hz, z_bare):
    omega = 2.0 * np.pi * f_hz
    return parallel_rows(f_hz, add_c1(np.asarray(z_bare), omega))


@complex_only
def test_p0_stored_route_b_reproduces_dense(data):
    f_hz = data["f_hz"]
    fails = []
    for col in COLUMNS:
        s = data[col]["s"]
        z_in = np.array([route_b_s(s[k], f_hz[k])[0] for k in range(len(f_hz))])
        dev = np.max(np.abs(z_in - data[col]["z_in_stored"]) / np.abs(data[col]["z_in_stored"]))
        z_in_z = np.array([route_b_z(z_from_s(s[k]), f_hz[k])[0] for k in range(len(f_hz))])
        dev_z = np.max(np.abs(z_in_z - data[col]["z_in_stored"])
                       / np.abs(data[col]["z_in_stored"]))
        rows = parallel_rows(f_hz, z_in)
        rel = _row_rel(rows, DENSE_ROWS[col])
        _say(f"[port23] (P0) {col}: S-route Z_in vs stored max rel {dev:.3e} "
             f"(<= {STORED_Z_IN_RTOL:.0e}); Z-route vs stored {dev_z:.3e} (printed); "
             f"rows vs DENSE_ROWS f0 {rel['f0_hz']:.3e} R_in {rel['r_in_ohm']:.3e} "
             f"Q {rel['q_loaded']:.3e} (BANDS {BANDS}); up-crossings "
             f"{rows['crossings_up_hz']}")
        if not dev <= STORED_Z_IN_RTOL:
            fails.append((col, "stored Z_in", dev))
        for k in ROWS:
            if not rel[k] <= BANDS[k]:
                fails.append((col, k, rel[k]))
    assert not fails, f"(P0) BLOCKED: {fails}"


@complex_only
def test_r1_r2_r3_skrf_reproduces_route_b(data):
    f_hz = data["f_hz"]
    fails = []
    for col in COLUMNS:
        s = data[col]["s"]
        zb_raw = np.array([route_b_z(z_from_s(s[k]), f_hz[k])[1] for k in range(len(f_hz))])
        zb_skrf = _skrf_z_bare(s, f_hz)
        r1 = np.abs(zb_skrf - zb_raw) / np.abs(zb_raw)
        _say(f"[port23] (R1) {col}: skrf Z_bare vs route_b_z max rel {r1.max():.3e} at "
             f"{f_hz[int(np.argmax(r1))] / 1e6:.1f} MHz (<= {S_VS_Z_ROUTE_RTOL:.0e}, "
             f"predicted ~1e-13), 51 points")
        if not np.all(r1 <= S_VS_Z_ROUTE_RTOL):
            fails.append((col, "R1", float(r1.max())))

        rows_raw = _rows_from_bare(f_hz, zb_raw)
        rows_skrf = _rows_from_bare(f_hz, zb_skrf)
        r2 = _row_rel(rows_skrf, rows_raw)
        same_bracket = rows_skrf.get("bracket_hz") == rows_raw.get("bracket_hz")
        _say(f"[port23] (R2) {col}: skrf vs raw rows f0 {r2['f0_hz']:.3e} R_in "
             f"{r2['r_in_ohm']:.3e} Q {r2['q_loaded']:.3e} (<= {ROW_AGREE_RTOL:.0e}); "
             f"bracket skrf {rows_skrf.get('bracket_hz')} raw {rows_raw.get('bracket_hz')} "
             f"same {same_bracket}")
        for k in ROWS:
            if not r2[k] <= ROW_AGREE_RTOL:
                fails.append((col, "R2", k, r2[k]))
        if not same_bracket:
            fails.append((col, "R2", "bracket"))

        r3 = _row_rel(rows_skrf, DENSE_ROWS[col])
        _say(f"[port23] (R3) {col}: skrf rows f0 {rows_skrf.get('f0_hz')!r} Hz, R_in "
             f"{rows_skrf.get('r_in_ohm')!r} ohm, Q {rows_skrf.get('q_loaded')!r}; vs "
             f"DENSE_ROWS f0 {r3['f0_hz']:.3e} R_in {r3['r_in_ohm']:.3e} Q "
             f"{r3['q_loaded']:.3e} (BANDS {BANDS})")
        for k in ROWS:
            if not r3[k] <= BANDS[k]:
                fails.append((col, "R3", k, r3[k]))

        # Negative controls -- predicted, printed beside, never asserted (rule (e)).
        zb_n1 = _skrf_z_bare(s, f_hz, z0=N1_Z0)
        rows_n1 = _rows_from_bare(f_hz, zb_n1)
        n1_zb = float(np.max(np.abs(zb_n1 - zb_raw) / np.abs(zb_raw)))
        n1_rows = _row_rel(rows_n1, rows_skrf)
        n1 = max(n1_rows.values())
        _say(f"[port23] (N1, predicted) {col}: S labelled z0 = {N1_Z0:g} ohm -> up-crossings "
             f"{rows_n1['crossings_up_hz']} (predicted none in 62–67 MHz; k = 1 near "
             f"{64 / np.sqrt(1.5):.1f} MHz); Z_bare shift max rel {n1_zb:.3e}; row shift "
             f"{n1:.3e} (blind if <= {BLIND_FLOOR:.0e})")
        zb_n2 = _skrf_z_bare(s, f_hz, port_map=N2_PORT_MAP)
        rows_n2 = _rows_from_bare(f_hz, zb_n2)
        n2_zb = float(np.max(np.abs(zb_n2 - zb_raw) / np.abs(zb_raw)))
        n2_rows = _row_rel(rows_n2, rows_skrf)
        n2 = max(n2_rows.values())
        _say(f"[port23] (N2, predicted) {col}: port_map {N2_PORT_MAP} (P1<->P3) -> Z_bare "
             f"shift max rel {n2_zb:.3e}; row shift f0 {n2_rows['f0_hz']:.3e} R_in "
             f"{n2_rows['r_in_ohm']:.3e} Q {n2_rows['q_loaded']:.3e}, max {n2:.3e} "
             f"(predicted >> 1e-9; blind if <= {BLIND_FLOOR:.0e})")
        if n1 <= BLIND_FLOOR or n2 <= BLIND_FLOOR:
            _say(f"[port23] GATE BLIND ({col}): N1 {n1:.3e}, N2 {n2:.3e} -- stop and report")
    assert not fails, f"PORT-23 step 2: {fails}"
