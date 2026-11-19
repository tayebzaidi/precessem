"""`ANS-8` step 2 runnable half (§9 item 369, §10 H21): SPEC v3 row algebra.

(V1) P_acc algebra: on a synthetic 2-feed network driven through the
Thevenin pair (w, Z_th), ½ Re Σ V_port conj(I_line) equals the wave form
(‖v_s‖² − ‖S v_s‖²) / (8 z0), v_s = w·``_source_voltage``, S from
``_current_route_s``, to ≤ ``Y0_CODING_BAND``.  Negative control (asserted):
item 361's P_acc (no feed-ESR term) misses the wave form by > 1e-6.
(V2) row-5 fractions of a synthetic partition with P_acc = Σ terms sum to 1.
(V3) ``comparison_text`` with a synthetic ``rows_v3``: the six rows land in
the table with every AED cell empty, the v1 / v2 text is an unchanged prefix,
and no planted AED sentinel reaches the text.
(V4, §9 item 386) ``sheet_variance_split`` recovers a synthetic sheet
field's transverse term and drive-component variance to ≤ 1e-12.
"""

from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

from tests.validation.test_wf10_f_human_copper_esr_root import (
    Y0_CODING_BAND,
    _current_route_s,
    _source_voltage,
    _thevenin_feed_impedance,
    _z_rc,
)

CASE = (
    Path(__file__).resolve().parents[2]
    / "examples"
    / "ansys_benchmarks"
    / "ans8_f_human_quadrature_mode_content_64MHz"
)
SCRIPT = CASE / "08_f_human_quadrature_mode_content.py"
Z0 = 50.0


def _mod():
    spec = importlib.util.spec_from_file_location("ans8_v3_ut", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ans8_v3_ut"] = mod
    spec.loader.exec_module(mod)
    return mod


def _network():
    rng = np.random.default_rng(369)
    a = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    z_net = 0.5 * (a + a.T) * 3.0 + np.diag([4.0 + 2.0j, 5.0 - 1.0j])  # reciprocal, lossy
    z_rc = _z_rc(118.2437e-12)
    z_th = _thevenin_feed_impedance(z_rc, Z0)
    v_s = _source_voltage(z_rc, Z0)
    m = z_net + z_th * np.eye(2)
    curs = [np.linalg.solve(m, e) for e in np.eye(2, dtype=np.complex128)]
    s2 = _current_route_s(curs, [0, 1], z_th, v_s)
    return z_net, z_rc, z_th, v_s, curs, s2


def test_v1_p_acc_line_cut_equals_wave_form():
    mod = _mod()
    z_net, z_rc, z_th, v_s, curs, s2 = _network()
    for w in ((1.0, -1.0j), (1.0, 1.0j), (1.0, 0.0), (0.3 + 0.2j, -0.7j)):
        w = np.asarray(w, dtype=np.complex128)
        cur = w[0] * curs[0] + w[1] * curs[1]
        p_line, v, _il = mod.p_acc_line(w, cur, z_th, z_rc)
        p_wave = mod.p_acc_wave(s2, w, v_s, Z0)
        rel = abs(p_line - p_wave) / abs(p_wave)
        # item 361's P_acc: power into the coil terminals only (no feed-gap ESR).
        p_361 = float(0.5 * np.real(np.sum(v * np.conj(cur))))
        miss = abs(p_361 - p_wave) / abs(p_wave)
        print(
            f"V1 w={w}: line {p_line:.15e} wave {p_wave:.15e} rel {rel:.3e} (band "
            f"{Y0_CODING_BAND:g}); item-361 form miss {miss:.3e} (> 1e-6 ASSERTED)"
        )
        assert p_wave > 0.0
        assert rel <= Y0_CODING_BAND
        assert miss > 1e-6


def test_v2_row5_fractions_sum_to_one():
    mod = _mod()
    terms = dict(p_rc30=2.1e-2, p_esr_feeds=1.27e-3, p_phantom=8.2e-3, p_surf=6.1e-4, p_rest=3.3e-9)
    p_acc = sum(terms.values())
    fr = mod.row5_fractions(p_acc, **terms)
    print(f"V2 fractions {fr}; |sum - 1| {abs(fr['sum'] - 1.0):.3e}")
    assert abs(fr["sum"] - 1.0) <= 4.0 * np.finfo(float).eps
    assert abs(fr["esr"] - (terms["p_rc30"] + terms["p_esr_feeds"]) / p_acc) == 0.0


def _rows_v3():
    sp = {"share_k1": 0.99, "leak_k0": 1e-3, "fractions": [0.0] * 16, "abs_k_folds": [0.0] * 9}
    return {
        "c_in_pf": 118.2437,
        "r_ohm": 0.2,
        "cells": 369387,
        "feeds": ["P17", "P21"],
        "row1_y_in_ccw_s": {"P17": [1.0e-3, 1.0e-5], "P21": [1.1e-3, -2.0e-4]},
        "row2_field": {"bottom": dict(sp), "top": dict(sp)},
        "row3": {n: {"b1p_per_sqrt_w": 1.6e-6} for n in ("O", "P1", "P2", "P3", "P4")},
        "row4": {"c_star_m": [0.044, 0.0088, 0.0], "sar10g_per_w": 0.16},
        "row5": {"phantom": 0.27, "copper": 0.02, "esr": 0.71, "rest": 0.0, "sum": 1.0},
        "row6_purity_o": 21.5,
        "power_w": {"p_acc_spec": 3.03e-2},
    }


def test_v3_comparison_text():
    mod = _mod()
    base = json.loads((CASE / "metrics.json").read_text())
    base.pop("rows_v3", None)
    before = mod.comparison_text(base)
    sentinel = 7.777777e-07
    m = copy.deepcopy(base)
    m["rows_v3"] = _rows_v3()
    m["aed"] = {"v3": {"row3_O": sentinel, "row4": sentinel}}
    after = mod.comparison_text(m)
    assert after.startswith(before.rstrip("\n")), "v1 / v2 lines must be an unchanged prefix"
    added = after[len(before.rstrip("\n")) :].splitlines()
    rows = [ln for ln in added if ln.startswith("| (")]
    labels = {ln.split(")")[0] + ")" for ln in rows}
    print(f"V3: {len(rows)} v3 table rows, labels {sorted(labels)}")
    assert labels == {f"| ({k})" for k in range(1, 7)}
    assert all(ln.endswith("| |") for ln in rows), "every AED cell empty"
    for tok in ("7.777777e-07", "7.78e-07", "7.777777e-7", "0.000000777777"):
        assert tok not in after
    assert "lossless (step 1)" in after


def test_v4_sheet_variance_split_recovers_both_terms():
    """(V4, §9 item 386 amendment 3) the discrete analogue of the split: on a
    sampled sheet field E = f·ĥ + g·p̂ (p̂ ⟂ ĥ, both in the sheet plane) with
    quadrature weights, ``sheet_variance_split`` recovers the transverse term
    ½Re(Y)Σw|g|² and the drive variance ½Re(Y)(Σw|f|² − |Σw f|²/A) to ≤ 1e-12."""
    mod = _mod()
    rng = np.random.default_rng(386)
    n = 400
    wts = rng.uniform(0.5, 1.5, size=n) * 1e-6  # quadrature weights [m²]
    f = (2.0 + 0.3 * rng.normal(size=n)) + 1j * (0.5 + 0.2 * rng.normal(size=n))
    g = 0.1 * rng.normal(size=n) + 0.05j * rng.normal(size=n)
    h_hat = np.array([0.0, 0.0, 1.0])
    p_hat = np.array([np.cos(0.3), np.sin(0.3), 0.0])
    n_hat = np.cross(h_hat, p_hat)
    e = f[:, None] * h_hat + g[:, None] * p_hat + 0.7 * rng.normal(size=n)[:, None] * n_hat
    nxe = np.cross(n_hat, e)
    e_h = e @ h_hat
    area = wts.sum()
    j = np.sum(wts * e_h)
    l_tan = np.sum(wts * np.sum(np.abs(nxe) ** 2, axis=1))
    t_drv = np.sum(wts * np.abs(e_h) ** 2)
    y_s = 1.0 / (0.4 + 0.1j)
    cs, trans, var = mod.sheet_variance_split(area, j, l_tan, t_drv, y_s)
    half = 0.5 * np.real(y_s)
    trans_ref = half * np.sum(wts * np.abs(g) ** 2)
    var_ref = half * (np.sum(wts * np.abs(f) ** 2) - abs(np.sum(wts * f)) ** 2 / area)
    r_t = abs(trans / trans_ref - 1.0)
    r_v = abs(var / var_ref - 1.0)
    r_cs = abs(cs / (trans_ref + var_ref) - 1.0)
    print(
        f"V4: transverse {trans:.15e} (ref {trans_ref:.15e}, rel {r_t:.3e}); variance "
        f"{var:.15e} (ref {var_ref:.15e}, rel {r_v:.3e}); CS rel {r_cs:.3e} (band 1e-12)"
    )
    assert r_t <= 1e-12 and r_v <= 1e-12 and r_cs <= 1e-12
    # Uniform drive field, no transverse part: both terms vanish (CS = 0).
    cs0, t0, v0 = mod.sheet_variance_split(area, 2.0 * area, 4.0 * area, 4.0 * area, y_s)
    assert abs(cs0) <= 1e-12 * half * 4.0 * area and abs(t0) == 0.0
