"""`ANS-8` H30a (§9 item 576, §10 chain H, provisional): the v3 record route.

SPEC v3 rows 4–5 per accepted watt, the line-cut − partition gap and X_s /
ΣX_s / P_acc read from a node record (``h9p-node-*.json``) and its harness
log(s) — no builder, no field, no solve.  The record's ``s_node`` is the
current route's S, so under the drive w ``V_port = ½ V_s (𝟙 + S) w`` and
``I = (w − V_port) / Z_th``; the formulas are ``_stage_v3``'s with these V / I.

* ``test_ans8_h30a_coding`` (ungated): (C1) ``record_port_state`` on item 369's
  synthetic two-feed network, (C2) ``node_log_readings`` on a synthetic log,
  (C3) the knob parser.
* ``test_ans8_h30a_degree1_reproduction`` (``FEM_EM_ANS8_H30A=1``): the degree-1
  record against item 386's ``v3-rows.json`` at ``REUSE_REPRODUCTION_RTOL``.
* ``test_ans8_h30a_degree2_rehearsal`` (``FEM_EM_ANS8_H30A=1``): entry 35's
  record (degree 2 at 118.2437 pF, detuned) and entry 42's (degree 2 at the
  root, no lattice in its log) — printed, provenance / algebra asserted.
* ``test_ans8_h30_degree2_root`` (``FEM_EM_ANS8_H30=1``): H30's own window on
  entry 41's record — lands skipped; H30 runs it.

Rows 3 / 6 and CS_s are in no record: printed "not in record", never filled.
Not claimed: H30's reading, any degree-2 root row as a result, the
known-issues 2026-10-11 retirement, any AED value or verdict.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

from tests.log_records import glob_log_paths, read_log_text
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_drive_superposition import UNIT_WEIGHT_CONTROL_RTOL
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
from tests.validation.test_wf10_f_human_copper_esr_root import (
    JSON_DIR,
    Y0_CODING_BAND,
    _active_y_in,
    _current_route_s,
    _source_voltage,
    _thevenin_feed_impedance,
    _z_rc,
)
from tests.validation.test_wf10_f_human_tuned_b1_sar import (
    C_IN_F_HUMAN_D2_PF,
    C_IN_F_HUMAN_PF,
    _node_json_path,
)
from tests.validation.test_wf10_lattice_d2 import LOG_DIR, parse_table
from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR

CASE = (
    Path(__file__).resolve().parents[2]
    / "examples"
    / "ansys_benchmarks"
    / "ans8_f_human_quadrature_mode_content_64MHz"
)
SCRIPT = CASE / "08_f_human_quadrature_mode_content.py"
KNOB = "FEM_EM_ANS8_H30A"
H30_KNOB = "FEM_EM_ANS8_H30"
V3_ROWS_JSON = "/workspace/logs/ans8/v3-rows.json"
H30_JSON = "/workspace/logs/ans8/v3-rows-d2-115.5209.json"

# Logs (tracked, docs/testing/logs/).
D1_NODE_LOG = "20261018T094815Z_MAT-15-step6b-W3-regate.log"  # wrote :6299
D1_LATTICE_LOG = "20261013T111006Z_WF-10-step4b-W1-regate.log"  # c* :6390
R35_LOG = "20261016T070008Z_WF-10-step4-H20.log"  # node + lattice
R42_LOG = "20261019T070008Z_MAT-15-step6-P22-d2root.log"  # node only
H30_LOG_GLOB = "20261022T*_WF-10-step4-H28.log"

# Predictions (brief, drafter's unlogged stdlib arithmetic) — printed, never asserted.
PREDICTED = {
    "d1": {
        "p_acc": 3.057629563251e-02,
        "x": (1.329324349e-04, 1.213651853e-04),
        "sum_x_over_p_acc": 8.316822e-03,
        "row5": (0.268112, 0.019926, 0.703644, 0.991683),
        "row4": 0.1587143,
        "nc361": 4.9700e-02,
        "esr_blind": 0.71196,
    },
    "R35": {
        "p_acc": 1.907593125981e-02,
        "x": (2.050749253e-04, 1.899495613e-04),
        "sum_x_over_p_acc": 2.070801e-02,
        "row5": (0.250469, 0.019414, 0.709409, 0.979292),
        "row4": 0.150620,
        "nc361": 9.7512e-02,
        "esr_blind": 0.73012,
    },
    "R42": {
        "p_acc": 2.910591023476e-02,
        "x": (3.262132656e-04, 2.994809404e-04),
        "sum_x_over_p_acc": 2.149715e-02,
        "row5": (0.260954, 0.020083, 0.697465, 0.978503),
        "row4": None,
        "nc361": 6.2316e-02,
        "esr_blind": 0.71896,
    },
}
#: For H30 (printed beside R42 only): the degree-2 root row-4 window and point.
PREDICTED_H30_ROW4 = ((0.145, 0.170), 0.159)
PREDICTED_H30_SUM_X_OVER_P = 2.15e-2


def _mod():
    spec = importlib.util.spec_from_file_location("ans8_h30a", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ans8_h30a"] = mod
    spec.loader.exec_module(mod)
    return mod


def _enabled(knob=KNOB) -> bool:
    raw = os.environ.get(knob, "").strip()
    if not raw:
        return False
    if raw != "1":
        raise ValueError(f"{knob} must be unset or '1', got {raw!r}")
    return True


def _need(knob=KNOB):
    if not _enabled(knob):
        pytest.skip(f"{knob} unset (reads gitignored node records under {JSON_DIR})")


def _log(name):
    return os.path.join(LOG_DIR, name)


def _node_log_readings(mod, name, node_basename):
    """`mod.node_log_readings` on a log in either form: it opens a plain path, so
    it reads a scratch copy of the text (`OPS-73`), never a decompressed record."""
    with tempfile.TemporaryDirectory() as d:
        scratch = os.path.join(d, name)
        with open(scratch, "w") as fh:
            fh.write(read_log_text(_log(name)))
        return mod.node_log_readings(scratch, node_basename)


def _load(path):
    assert os.path.exists(path), f"{path} absent: item 576 BLOCKED; never re-solve"
    with open(path) as fh:
        return json.load(fh)


def _rel(a, b):
    return abs(a / b - 1.0)


def _check_record(mod, tag, rec, json_path, node_log, lattice_log, expect):
    """(Prov) / (Rec-Y) / (Lin) / (P_acc-1) / (Prov-L) / NC 361 on one record.
    Returns ``(rows, lattice_or_None, log_readings)``."""
    rd = _node_log_readings(mod, node_log, os.path.basename(json_path))
    print(
        f"[H30a {tag}] {node_log}: last wrote at :{rd['wrote_line']} ({rd['n_wrote']} wrote "
        f"line(s)) -> {rd['wrote_path']}"
    )
    # (Prov)
    assert rd["wrote_path"] == json_path, (rd["wrote_path"], json_path)
    got = (rec.get("degree"), rec["c_in_pf"], rec.get("node_knob", None))
    print(f"[H30a {tag}] (Prov) degree / c_in_pf / node_knob {got} (expected {expect})")
    assert got == expect, (got, expect)
    part = rec["partitions"]["ccw"]
    prov = {
        k: _rel(rd["p_w"][k if k != "p_sheets" else "p_sheets"], part[k])
        for k in ("p_src", "p_sheets", "p_phantom", "p_surf")
    }
    print(
        f"[H30a {tag}] (Prov) P_w line vs JSON partitions[ccw]: "
        + ", ".join(f"{k} {v:.3e}" for k, v in prov.items())
        + f" (band {REUSE_REPRODUCTION_RTOL:g} ASSERTED; PREDICTED <= 6e-10)"
    )
    assert max(prov.values()) <= REUSE_REPRODUCTION_RTOL, prov
    # (Rec-Y)
    s_node = np.array([[complex(a, b) for a, b in row] for row in rec["s_node"]])
    w = np.array([complex(a, b) for a, b in rec["w"]["ccw"]])
    y_rec = [complex(a, b) for a, b in rec["y_in"]["ccw"]]
    y_now = _active_y_in(s_node, w, float(REFERENCE_IMPEDANCE_OHM))
    r_y = max(abs(a - b) for a, b in zip(y_now, y_rec)) / max(abs(b) for b in y_rec)
    print(
        f"[H30a {tag}] (Rec-Y) y_in[ccw] vs _active_y_in(s_node, w_ccw): {r_y:.3e} (band "
        f"{Y0_CODING_BAND:g} ASSERTED; PREDICTED 0); w_ccw[1] = {w[1]!r}"
    )
    assert r_y <= Y0_CODING_BAND
    # (Lin) — the window's own field-side gate, re-asserted from its log.
    assert rd["lin"], f"{node_log}: no (Lin) line"
    for sense, r_lin, r_yl in rd["lin"]:
        print(
            f"[H30a {tag}] (Lin, field side, from the log) [{sense}] currents {r_lin:.3e} (band "
            f"{UNIT_WEIGHT_CONTROL_RTOL:g}); Y_in route {r_yl:.3e} (band {Y0_CODING_BAND:g})"
        )
        assert r_lin <= UNIT_WEIGHT_CONTROL_RTOL and r_yl <= Y0_CODING_BAND
    rows = mod.v3_record_rows(rec)
    # (P_acc-1), record route: the S <-> V / I algebra.
    print(
        f"[H30a {tag}] (P_acc-1, record route: S <-> V/I algebra) line cut "
        f"{rows['power_w']['p_acc_spec']:.12e} W vs wave form {rows['power_w']['p_acc_wave']:.12e}"
        f" W: rel {rows['r_acc1']:.3e} (band {DISCRETE_IDENTITY_RTOL:g} ASSERTED; PREDICTED <= 2e-15)"
    )
    assert rows["r_acc1"] <= DISCRETE_IDENTITY_RTOL
    nc = rows["negative_controls"]
    bar = BLIND_FACTOR * DISCRETE_IDENTITY_RTOL
    print(
        f"[H30a {tag}] NC item 361's P_acc in (P_acc-1) (ASSERTED >= {bar:g}): miss "
        f"{nc['item361_p_acc_in_p_acc_1']:.4e} (PREDICTED {PREDICTED[tag]['nc361']:.4e}); "
        f"ESR-blind (printed) {nc['esr_blind']:.5f} (PREDICTED {PREDICTED[tag]['esr_blind']}); "
        f"feed-ESR-blind {nc['feed_esr_blind']:.4e}"
    )
    assert nc["item361_p_acc_in_p_acc_1"] >= bar
    lat = None
    if lattice_log is not None:
        lat = parse_table(_log(lattice_log))
        r_l = _rel(lat["p_phantom_ccw"], part["p_phantom"])
        print(
            f"[H30a {tag}] (Prov-L) {lattice_log}: {len(lat['labels'])} rows, c* "
            f"{lat['cstar_label']} SAR10g {lat['cstar_sar']:.9e} W/kg; P_phantom ccw "
            f"{lat['p_phantom_ccw']:.9e} vs JSON {part['p_phantom']:.9e}: rel {r_l:.3e} (band "
            f"{REUSE_REPRODUCTION_RTOL:g} ASSERTED)"
        )
        assert not lat["reduced"]
        assert r_l <= REUSE_REPRODUCTION_RTOL
        rows = mod.v3_record_rows(rec, sar10g_cstar_w_per_kg=lat["cstar_sar"])
    return rows, lat, rd


def _print_rows(tag, rows, label=""):
    p = rows["power_w"]
    pr = PREDICTED.get(tag)
    print(
        f"[H30a {tag}] rows{label}: P_acc,SPEC {p['p_acc_spec']:.12e} W"
        + (f" (PREDICTED {pr['p_acc']:.12e})" if pr else "")
        + f"; P_src {p['p_src']:.9e}, P_RC,30 {p['p_rc30']:.9e}, P_ESR,feeds {p['p_esr_feeds']:.9e},"
        f" P_phantom {p['p_phantom']:.9e}, P_surf {p['p_surf']:.9e}, P_rest {p['p_rest']:.9e} W"
    )
    f, fp = rows["row5"], rows["row5_partition"]
    print(
        f"[H30a {tag}] row 5 (of P_acc): phantom {f['phantom']:.6f} / copper {f['copper']:.6f} / "
        f"ESR {f['esr']:.6f}, sum {f['sum']:.6f}"
        + (f" (PREDICTED {pr['row5']})" if pr else "")
        + f"; of the partition: {fp['phantom']:.6f} / {fp['copper']:.6f} / {fp['esr']:.6f}, sum "
        f"{fp['sum']:.15f}"
    )
    for k, (pid, m) in enumerate(rows["gap"]["mechanism"].items()):
        print(
            f"[H30a {tag}] X_s {pid}: {m['x_w']:.9e} W = P_sheet,field {m['p_sheet_field_w']:.9e} − "
            f"lumped Norton {m['lumped_norton_w']:.9e} W"
            + (f" (PREDICTED {pr['x'][k]:.9e})" if pr else "")
        )
    gap = rows["gap"]["line_cut_minus_partition_w"]
    print(
        f"[H30a {tag}] gap (line cut − partition) {gap:.9e} W, rel {rows['r_acc2']:.4e}; ΣX_s "
        f"{rows['sum_x_w']:.9e} W; ΣX_s / P_acc {rows['sum_x_over_p_acc']:.6e}"
        + (f" (PREDICTED {pr['sum_x_over_p_acc']:.6e})" if pr else "")
        + " — printed; never read against the weekly's threshold here"
    )
    print(
        f"[H30a {tag}] (printed) P_src vs Norton lumped {rows['norton_lumped_w']:.12e} W: rel "
        f"{_rel(p['p_src'], rows['norton_lumped_w']):.3e}; ΣX_s vs gap rel "
        f"{_rel(rows['sum_x_w'], gap):.3e}; ΣX_s vs p_361 − terminal rel "
        f"{_rel(rows['sum_x_w'], p['p_acc_item361'] - rows['terminal_w']):.3e} (PREDICTED <= 1e-12)"
    )
    if "row4" in rows:
        r4 = rows["row4"]
        print(
            f"[H30a {tag}] row 4: SAR10g(c*) {r4['sar10g_w_per_kg']:.9e} W/kg; per W of P_acc "
            f"{r4['sar10g_per_w']:.7f} W/kg/W"
            + (f" (PREDICTED {pr['row4']})" if pr else "")
            + f"; per W of the partition {r4['sar10g_per_w_partition']:.7f}"
        )
    else:
        print(f"[H30a {tag}] row 4: not in this log (no 10 g lattice)")
    print(f"[H30a {tag}] not in record: {rows['not_in_record']} (need a field; never filled)")


# ------------------------------------------------------------------------------
def _network():
    """Item 369's synthetic two-feed network (``tests/unit/test_ans8_v3_rows.py``,
    seed 369), copied, not imported."""
    rng = np.random.default_rng(369)
    a = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    z_net = 0.5 * (a + a.T) * 3.0 + np.diag([4.0 + 2.0j, 5.0 - 1.0j])
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    z_rc = _z_rc(118.2437e-12)
    z_th = _thevenin_feed_impedance(z_rc, z0)
    v_s = _source_voltage(z_rc, z0)
    m = z_net + z_th * np.eye(2)
    curs = [np.linalg.solve(m, e) for e in np.eye(2, dtype=np.complex128)]
    s2 = _current_route_s(curs, [0, 1], z_th, v_s)
    return z_th, v_s, curs, s2


def test_ans8_h30a_coding(monkeypatch, tmp_path):
    mod = _mod()
    # (C1)
    z_th, v_s, curs, s2 = _network()
    for w in ((1.0, -1.0j), (1.0, 1.0j), (1.0, 0.0), (0.3 + 0.2j, -0.7j)):
        w = np.asarray(w, dtype=np.complex128)
        cur = w[0] * curs[0] + w[1] * curs[1]
        v_true = w - z_th * cur
        v, i = mod.record_port_state(s2, w, z_th, v_s)
        rv = np.max(np.abs(v - v_true)) / np.max(np.abs(v_true))
        ri = np.max(np.abs(i - cur)) / np.max(np.abs(cur))
        print(f"(C1) w={w}: V rel {rv:.3e}, I rel {ri:.3e} (band {Y0_CODING_BAND:g} ASSERTED)")
        assert rv <= Y0_CODING_BAND and ri <= Y0_CODING_BAND
    # (C2) a synthetic log in the window's formats: two writes, P_w before each.
    pw = (
        "[WF-10 step2] [ccw] w = (1.000+0.000j, 0.000-1.000j); (P_w) P_src {a:.9e} sumP_sheet "
        "{b:.9e} P_phantom {c:.9e} P_rest 0.000000000e+00 P_surf {d:.9e} W; rel 4.854e-13 (band "
        "1e-06, ASSERTED; PREDICTED <= 1e-10); P_surf-blind miss 7.0397e-03 (PREDICTED ≈ 0.007, printed)"
    )
    lin = (
        "[WF-10 step2] [{s}] (Lin) currents rel {r:.3e} (band 1e-12 UNIT_WEIGHT_CONTROL_RTOL, "
        "ASSERTED; PREDICTED <= 1e-13); Y_in current route vs _active_y_in 1.184e-15 (band 1e-12, "
        "ASSERTED)"
    )
    wrote = "[WF-10 step2] wrote /workspace/logs/wf10/h9p-node-1.0000-d2.json"
    text = (
        "\n".join(
            [
                pw.format(a=1.0, b=2.0, c=3.0, d=4.0),
                lin.format(s="ccw", r=1e-16),
                wrote,
                pw.format(a=5.0, b=6.0, c=7.0, d=8.0),
                lin.format(s="cw", r=2e-16),
                wrote,
                pw.format(a=9.0, b=9.0, c=9.0, d=9.0),
                "tail",
            ]
        )
        + "\n"
    )
    path = tmp_path / "synthetic.log"
    path.write_text(text)
    rd = mod.node_log_readings(str(path), "h9p-node-1.0000-d2.json")
    print(f"(C2) {rd}")
    assert rd["wrote_line"] == 6 and rd["n_wrote"] == 2
    assert rd["wrote_path"] == "/workspace/logs/wf10/h9p-node-1.0000-d2.json"
    assert rd["p_w"] == {
        "p_src": 5.0,
        "p_sheets": 6.0,
        "p_phantom": 7.0,
        "p_rest": 0.0,
        "p_surf": 8.0,
    }
    assert rd["lin"] == [("ccw", 1e-16, 1.184e-15), ("cw", 2e-16, 1.184e-15)]
    with pytest.raises(AssertionError):
        mod.node_log_readings(str(path), "h9p-node-2.0000-d2.json")
    # (C3) the knob parser.
    monkeypatch.delenv(KNOB, raising=False)
    assert _enabled() is False
    monkeypatch.setenv(KNOB, "1")
    assert _enabled() is True
    monkeypatch.setenv(KNOB, "yes")
    with pytest.raises(ValueError):
        _enabled()
    print("(C3) knob: unset -> skip, '1' -> run, 'yes' -> ValueError")


def test_ans8_h30a_degree1_reproduction():
    _need()
    from mpi4py import MPI

    assert MPI.COMM_WORLD.size == 1, "record route: -n 1"
    mod = _mod()
    json_path = _node_json_path(None)
    rec = _load(json_path)
    ref = _load(V3_ROWS_JSON)
    rows, lat, _rd = _check_record(
        mod, "d1", rec, json_path, D1_NODE_LOG, D1_LATTICE_LOG, (1, C_IN_F_HUMAN_PF, "")
    )
    _print_rows("d1", rows)
    # (Rep-S)
    r_s = _rel(lat["cstar_sar"], mod.SAR10G_CSTAR_H9P_RECORD_W_PER_KG)
    print(
        f"[H30a d1] (Rep-S) c* {lat['cstar_label']} (record {mod.V3_CSTAR_LABEL_RECORD}); SAR10g "
        f"rel {r_s:.3e} (band {REUSE_REPRODUCTION_RTOL:g} ASSERTED)"
    )
    assert lat["cstar_label"] == mod.V3_CSTAR_LABEL_RECORD
    assert r_s <= REUSE_REPRODUCTION_RTOL
    # (Rep-d1) against item 386's v3-rows.json.
    pairs = {}
    for k in (
        "p_acc_spec",
        "p_acc_wave",
        "p_acc_item361",
        "p_src",
        "p_rc30",
        "p_esr_feeds",
        "p_phantom",
        "p_surf",
    ):
        pairs[f"power_w.{k}"] = (rows["power_w"][k], ref["power_w"][k])
    for blk in ("row5", "row5_partition"):
        for k in ("phantom", "copper", "esr"):
            pairs[f"{blk}.{k}"] = (rows[blk][k], ref[blk][k])
    for k in ("line_cut_minus_partition_w", "terminal_w", "norton_lumped_w"):
        pairs[f"gap.{k}"] = (rows["gap"][k], ref["gap"][k])
    for pid, m in rows["gap"]["mechanism"].items():
        for k in ("x_w", "lumped_norton_w", "p_sheet_field_w"):
            pairs[f"gap.mechanism.{pid}.{k}"] = (m[k], ref["gap"]["mechanism"][pid][k])
    pairs["row4.sar10g_per_w"] = (rows["row4"]["sar10g_per_w"], ref["row4"]["sar10g_per_w"])
    pairs["per_partition.row4_sar10g_per_w"] = (
        rows["row4"]["sar10g_per_w_partition"],
        ref["per_partition"]["row4_sar10g_per_w"],
    )
    worst = 0.0
    for name, (a, b) in pairs.items():
        r = _rel(a, b)
        worst = max(worst, r)
        print(
            f"[H30a d1] (Rep-d1) {name:<38s} record route {a:+.12e} vs v3-rows {b:+.12e}: rel {r:.3e}"
        )
    print(
        f"[H30a d1] (Rep-d1) worst {worst:.3e} over {len(pairs)} keys (band "
        f"{REUSE_REPRODUCTION_RTOL:g} ASSERTED; PREDICTED <= 1e-9)"
    )
    assert worst <= REUSE_REPRODUCTION_RTOL
    # NC (asserted, floor imported): feed-label swap on X_s.
    bar = BLIND_FACTOR * REUSE_REPRODUCTION_RTOL
    mech, rmech = rows["gap"]["mechanism"], ref["gap"]["mechanism"]
    f0, f1 = rows["feeds"]
    swap = min(_rel(mech[f1]["x_w"], rmech[f0]["x_w"]), _rel(mech[f0]["x_w"], rmech[f1]["x_w"]))
    print(
        f"[H30a d1] NC feed-label swap: X_s({f0}) <-> X_s({f1}) misses (Rep-d1) by {swap:.4e} "
        f"(ASSERTED >= {bar:g}; PREDICTED |X21/X17 − 1| 8.70e-02)"
    )
    assert swap >= bar
    # NC (asserted, floor imported): wrong-record provenance, d1 JSON vs entry 35's P_w.
    rd35 = _node_log_readings(mod, R35_LOG, os.path.basename(_node_json_path(2)))
    wrong = _rel(rd35["p_w"]["p_src"], rec["partitions"]["ccw"]["p_src"])
    print(
        f"[H30a d1] NC wrong record: d1 JSON P_src vs entry 35's last P_w line: miss {wrong:.4e} "
        f"(ASSERTED >= {bar:g}; PREDICTED 3.18e-02)"
    )
    assert wrong >= bar


def test_ans8_h30a_degree2_rehearsal():
    _need()
    from mpi4py import MPI

    assert MPI.COMM_WORLD.size == 1, "record route: -n 1"
    mod = _mod()
    bar = BLIND_FACTOR * REUSE_REPRODUCTION_RTOL
    # (R35) entry 35: degree 2 at the degree-1 node (detuned), log = node + lattice.
    p35 = _node_json_path(2)
    rec35 = _load(p35)
    rows35, _lat, rd35 = _check_record(
        mod, "R35", rec35, p35, R35_LOG, R35_LOG, (2, C_IN_F_HUMAN_PF, None)
    )
    _print_rows("R35", rows35, " (degree 2 at 118.2437 pF, detuned — a labelled record)")
    # (R42) entry 42: degree 2 at the root, no lattice.
    p42 = _node_json_path(2, c_pf=C_IN_F_HUMAN_D2_PF)
    rec42 = _load(p42)
    rows42, _l, _rd = _check_record(
        mod, "R42", rec42, p42, R42_LOG, None, (2, C_IN_F_HUMAN_D2_PF, "d2root")
    )
    _print_rows("R42", rows42, " (degree 2 at the root, entry 42's record)")
    wrong = _rel(rd35["p_w"]["p_src"], rec42["partitions"]["ccw"]["p_src"])
    print(
        f"[H30a R42] NC wrong record: R42 JSON P_src vs entry 35's last P_w line: miss "
        f"{wrong:.4e} (ASSERTED >= {bar:g}; PREDICTED 3.36e-03)"
    )
    assert wrong >= bar
    lo, hi = PREDICTED_H30_ROW4[0]
    print(
        f"[H30a R42] FOR H30 (printed, never asserted, never H30's reading): entry 41 should "
        f"reproduce R42 at {REUSE_REPRODUCTION_RTOL:g}; ΣX_s / P_acc predicted ≈ "
        f"{PREDICTED_H30_SUM_X_OVER_P:g} (R42 here {rows42['sum_x_over_p_acc']:.6e}); row 4 at the "
        f"root predicted in ({lo}, {hi}), point ≈ {PREDICTED_H30_ROW4[1]}"
    )


def test_ans8_h30_degree2_root():
    """H30's window (§10 H30) on entry 41's record — lands skipped; H30 runs it."""
    _need(H30_KNOB)
    from mpi4py import MPI

    assert MPI.COMM_WORLD.size == 1, "record route: -n 1"
    mod = _mod()
    found = glob_log_paths(_log(H30_LOG_GLOB))
    assert len(found) == 1, f"expected exactly one {H30_LOG_GLOB}, found {found}"
    log = os.path.basename(found[0])
    path = _node_json_path(2, c_pf=C_IN_F_HUMAN_D2_PF)
    rec = _load(path)
    assert rec.get("node_knob") == "d2root" and rec.get("degree") == 2, (
        rec.get("node_knob"),
        rec.get("degree"),
    )
    rows, _lat, _rd = _check_record(
        mod, "R42", rec, path, log, log, (2, C_IN_F_HUMAN_D2_PF, "d2root")
    )
    _print_rows("H30", rows, " (degree 2 at the root, entry 41's record)")
    ref = _load(V3_ROWS_JSON)
    print(
        f"[H30] degree-1 column (v3-rows.json): P_acc {ref['power_w']['p_acc_spec']:.12e} W; row 5 "
        f"{ref['row5']['phantom']:.6f} / {ref['row5']['copper']:.6f} / {ref['row5']['esr']:.6f}; "
        f"row 4 {ref['row4']['sar10g_per_w']:.7f} W/kg/W; ΣX_s / P_acc "
        f"{ref['gap']['rel_line_cut']:.6e}"
    )
    out = {k: v for k, v in rows.items() if k not in ("v_port", "i_feed")}
    out["v_port"] = [[complex(v).real, complex(v).imag] for v in rows["v_port"]]
    out["i_feed"] = [[complex(v).real, complex(v).imag] for v in rows["i_feed"]]
    out["node_log"] = log
    assert os.path.basename(H30_JSON) != mod.V3_JSON_NAME
    os.makedirs(os.path.dirname(H30_JSON), exist_ok=True)
    with open(H30_JSON, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(f"[H30] wrote {H30_JSON}")
