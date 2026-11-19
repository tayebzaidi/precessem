"""`ANS-9` step 1c (§9 item 396; §10 chain P23): far-side point SAR and θ₄₀ refined.

docs/planning/chunks/ANS-9.md, "§9 item 396".  Our values only: no AED number is
read or written here.

P21 (a) + (b)'s d2x window (371a / 371b) re-run through route B (`GEO-27` step 2's
builder: ``tests/mesh/test_birdcage_port_sheets._build(True, phantom_resolution=h_p,
leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)`` → ``build_four_port_sweep(built=…,
build_only=True)`` → ``run_n_port_sparameter_sweep(keep_fields=True)``, one
factorisation, P1's field with the other three ports on their 50 Ω sheets), on the
exact-rectangle F-small fixture at 10 MHz, κ off, at degree 2 (d2x) and degree 1 on
the same mesh (d1x).  Knob ``FEM_EM_ANS9_1C_RUNG`` ∈ {unset, ``x1``, ``h005``,
``h0025``}; unset ⇒ the gated tests skip; anything else ⇒ ``ValueError``.

Asserted (every band imported): (P1W), (G2), (i), (R1′), (E), (R3′) at every rung;
x1 only (F0′), (X1c), (F0-ctl), (d6′); h0025 only (R).  Printed, never asserted:
the order effects, the h_p moves against P21's comparands, (d6′) on the refined
rungs (ruling Q5), (Live-h), ρ(h) beside 2.119196 (``_branch`` is never called),
and at h0025 the three-way split against ``ans:9``'s degree-1 ragged column.
"""

from __future__ import annotations

import gc
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import _build
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_S_P1_RECORDS,
    ANS4_EXACT_STRIP_S_P2_X1,
)
from tests.validation.test_ans4_resolution_ladder import STEP3A_RECORD_RTOL
from tests.validation.test_ans9_far_side_sar import (
    PROBE_NAMES,
    _column_readout,
    _rel,
    _thermal_readout,
)
from tests.validation.test_birdcage_sar_mass_averaged import (
    EXACT_IDENTITY_RTOL,
    _global_tag_cell_count,
)
from tests.validation.test_geo27_one_gram_rung_degree2 import RUNGS
from tests.validation.test_imp3_palace_twin import _sheet_areas
from tests.validation.test_imp3_rod_birdcage import _peak_rss_gib
from tests.validation.test_imp3_tip_ladder_1g import I14_OFF_CELLS, I14_OFF_DOFS
from tests.validation.test_imp3_tip_ratio_1g import NOMINAL_STRIP_AREA
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, _robin_tags
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, S2_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import D_CONTROL_FLOOR, _tag_submesh
from tests.validation.test_mat19_probe_coil_phantom import ROW2_MEAN_K, ROW4_POINT_K
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD, STRIP_BAND
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

TAG = "[ANS-9 1c]"
KNOB_NAME = "FEM_EM_ANS9_1C_RUNG"
RUNG_MAP = {"x1": None, **RUNGS}

# --- Printed comparands (never bands), each with its log line -------------------
# 371a / 371b, item 192's fixture through `_four_port_rung`, x1 phantom (537 tag-3):
# docs/testing/logs/20261011T185145Z_ANS-9-step1b-W2b-d2x-thermal.log
P21_D2X_SAR = {
    "A": 3.696839972e-01,
    "B": 1.848801809e01,
    "C": 7.033595117e-03,  # :1978-1982
    "D": 3.064054166e00,
    "E": 2.584895127e00,
}
P21_D1X_SAR = {
    "A": 5.858834833e-01,
    "B": 1.771392340e01,
    "C": 3.028258304e-01,  # :2035-2039
    "D": 3.094049586e00,
    "E": 2.283028850e00,
}
P21_D2X_THETA40 = {
    "A": 1.070257,
    "B": 3.166010,
    "C": 0.929158,
    "D": 0.723341,  # :2046
    "E": 1.195038,
}
P21_D1X_THETA40 = {
    "A": 1.050658,
    "B": 3.189674,
    "C": 0.914637,
    "D": 0.620227,  # :2052
    "E": 1.203214,
}
P21_ROW2_D2X = {10: 0.327896, 20: 0.633639, 30: 0.920643, 40: 1.190592}  # :2043
P21_ROW2_D1X = {10: 0.327423, 20: 0.632311, 30: 0.918342, 40: 1.187285}  # :2049
P21_DC_D2X, P21_DA_D2X = 435.631297, 8.288306  # :1983
P21_DC_D1X, P21_DA_D1X = 10.217258, 5.280998  # :2040
# ans:9's own degree-1 ragged column on its 0.0025 m phantom (199 920 cells):
# docs/testing/logs/20261011T184301Z_ANS-9-step1b-W1-d1.log:2013-2018
ANS9_D1_SAR = {
    "A": 4.642360792e-01,
    "B": 1.717658032e01,
    "C": 1.289384239e-02,
    "D": 2.650515686e00,
    "E": 2.551981669e00,
}
ANS9_D1_DC = 205.564455
P21_RHO = 2.119196
# docs/testing/logs/20261012T034731Z_GEO-27-step2-W1-h005.log:1942
GEO27_H005_CELLS = (130427, 8498)

# --- Predictions (written before W0; printed, never asserted) -------------------
PRED_CELLS = {"x1": (117206, 537), "h005": (130427, 8498), "h0025": (201072, 58960)}
PRED_KAPPA = {"x1": "1.34", "h005": "between", "h0025": "1.0-1.15"}
PRED_MEAN_ORDER = (1.0028 - 0.005, 1.0028 + 0.005)
PRED_F0 = 1.0e-10
PRED_F0_CTL = 1.02e-2
PRED_D6 = 1.000


def _knob():
    raw = os.environ.get(KNOB_NAME)
    if raw is None or raw == "":
        return None
    if raw not in RUNG_MAP:
        raise ValueError(f"{KNOB_NAME}={raw!r}: expected one of {sorted(RUNG_MAP)}")
    return raw


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


needs_rung = pytest.mark.skipif(_knob() is None, reason=f"{KNOB_NAME} unset (x1 | h005 | h0025)")


def test_ans9_1c_coding(monkeypatch):
    """W0, ungated: knob parser, rung map, probe names."""
    for ok in ("x1", "h005", "h0025"):
        monkeypatch.setenv(KNOB_NAME, ok)
        assert _knob() == ok
    for bad in ("h001", "d2x", "X1"):
        monkeypatch.setenv(KNOB_NAME, bad)
        with pytest.raises(ValueError):
            _knob()
    monkeypatch.delenv(KNOB_NAME, raising=False)
    assert _knob() is None
    assert RUNG_MAP["x1"] is None
    assert RUNG_MAP["h005"] == 0.005
    assert RUNG_MAP["h0025"] == 0.0025
    assert PROBE_NAMES == ("A", "B", "C", "D", "E")


def _sweep(case, degree):
    return run_n_port_sparameter_sweep(
        case["problem"],
        case["port_defs"],
        lumped_sheet_ports=case["specs"],
        lumped_sheet_facet_tags=case["facet_tags"],
        keep_fields=True,
        degree=degree,
    )


def _recip(s):
    return float(np.linalg.norm(s - s.T) / np.linalg.norm(s))


@pytest.fixture(scope="module")
def rung():
    mode = _knob()
    comm = MPI.COMM_WORLD
    h_p = RUNG_MAP[mode]
    t_body = time.perf_counter()
    out = {"mode": mode}

    # 1-2. Build (sizing to _build only), degree 2.
    t0 = time.perf_counter()
    built = _build(True, phantom_resolution=h_p, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
    case = build_four_port_sweep(built=built, build_only=True, degree=2)
    comm.Barrier()
    t_build = time.perf_counter() - t0
    msh, ct = case["mesh"], case["cell_tags"]
    cells = int(case["cells"])
    tag3 = _global_tag_cell_count(msh, ct, PHANTOM_CELL_TAG, comm)
    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    n_e = int(msh.topology.index_map(1).size_global)
    n_f = int(msh.topology.index_map(2).size_global)
    areas = _sheet_areas(msh, case["facet_tags"], comm)
    out.update(cells=cells, tag3=tag3, n_e=n_e, n_f=n_f, areas=areas)
    _say(
        f"\n{TAG} rung {mode}: h_p = {h_p} m, {cells} cells ({tag3} tag-3; predicted "
        f"{PRED_CELLS[mode][0]} / {PRED_CELLS[mode][1]}), E = {n_e}, F = {n_f}, build {t_build:.2f} s"
    )
    if mode == "h005":
        _say(f"{TAG} cells / tag-3 vs GEO27_H005_CELLS {GEO27_H005_CELLS} (printed)")
    _say(
        f"{TAG} (G2) N_i areas "
        + ", ".join(f"{k} {a:.12e}" for k, a in areas.items())
        + f" vs {NOMINAL_STRIP_AREA:.6e} (STRIP_BAND {STRIP_BAND:.0e})"
    )

    # 3. d2x sweep, one factorisation.
    comm.Barrier()
    t0 = time.perf_counter()
    r2 = _sweep(case, 2)
    comm.Barrier()
    t_s2 = time.perf_counter() - t0
    s2 = np.asarray(r2.s_matrix, dtype=np.complex128)
    out["port_ids"] = tuple(r2.port_ids)
    out["s_d2x"] = s2
    v = r2.fields["P1"].e_complex.function_space
    out["dofs2"] = int(v.dofmap.index_map.size_global * v.dofmap.index_map_bs)
    out["recip2"] = _recip(s2)
    _say(
        f"{TAG} d2x sweep {t_s2:.2f} s, port_ids {out['port_ids']}, d2 dofs {out['dofs2']}; "
        f"(i) ||S-S^T||/||S|| = {out['recip2']:.3e} (band {RECIPROCITY_BAND:.0e})"
    )

    # 4-5. d2x readout + thermal half.
    f2 = r2.fields["P1"]
    out["d2x"] = _column_readout(f"{mode}-d2x", f2.e_complex, f2.sigma_field, ct, comm)
    t0 = time.perf_counter()
    sub, parent_ids = _tag_submesh(msh, ct, PHANTOM_CELL_TAG)
    ft = _robin_tags(sub)
    out["th_d2x"] = _thermal_readout(
        f"{mode}-d2x",
        msh,
        ct,
        f2.e_complex,
        f2.sigma_field,
        out["d2x"]["G"],
        out["d2x"]["pts"],
        sub,
        parent_ids,
        ft,
        comm,
        control=True,
    )
    t_th = time.perf_counter() - t0
    del f2, r2
    gc.collect()

    # 6. d1x on the same mesh.
    case1 = build_four_port_sweep(
        reuse={
            k: case[k] for k in ("mesh", "cell_tags", "facet_tags", "sheets", "halves", "cells")
        },
        build_only=True,
        degree=1,
    )
    comm.Barrier()
    t0 = time.perf_counter()
    r1 = _sweep(case1, 1)
    comm.Barrier()
    t_s1 = time.perf_counter() - t0
    s1 = np.asarray(r1.s_matrix, dtype=np.complex128)
    out["s_d1x"] = s1
    out["recip1"] = _recip(s1)
    _say(f"{TAG} d1x sweep {t_s1:.2f} s; (i) ||S-S^T||/||S|| = {out['recip1']:.3e}")
    f1 = r1.fields["P1"]
    out["d1x"] = _column_readout(
        f"{mode}-d1x", f1.e_complex, f1.sigma_field, case1["cell_tags"], comm
    )
    t0 = time.perf_counter()
    out["th_d1x"] = _thermal_readout(
        f"{mode}-d1x",
        msh,
        ct,
        f1.e_complex,
        f1.sigma_field,
        out["d1x"]["G"],
        out["d1x"]["pts"],
        sub,
        parent_ids,
        ft,
        comm,
    )
    t_th += time.perf_counter() - t0
    del f1, r1
    gc.collect()

    # 7. Prints.
    d2, d1 = out["d2x"], out["d1x"]
    t2, t1 = out["th_d2x"], out["th_d1x"]
    _say(
        f"{TAG} [{mode}] order effect point SAR d2x/d1x: "
        + "  ".join(f"{k}={d2['sar'][k] / d1['sar'][k]:.6f}" for k in PROBE_NAMES)
        + "  (x1 P21: D=0.990)"
    )
    _say(
        f"{TAG} [{mode}] order effect theta40 d2x/d1x: "
        + "  ".join(f"{k}={t2['row4'][k] / t1['row4'][k]:.6f}" for k in PROBE_NAMES)
        + f"; <theta>40 {t2['row2'][40] / t1['row2'][40]:.6f} (predicted in "
        f"[{PRED_MEAN_ORDER[0]:.4f}, {PRED_MEAN_ORDER[1]:.4f}]; x1 P21 theta_D 1.166252)"
    )
    _say(
        f"{TAG} [{mode}] theta40 D/C: d2x {t2['row4']['D'] / t2['row4']['C']:.6f}, "
        f"d1x {t1['row4']['D'] / t1['row4']['C']:.6f};  D/A: d2x {t2['row4']['D'] / t2['row4']['A']:.6f}, "
        f"d1x {t1['row4']['D'] / t1['row4']['A']:.6f}"
    )
    _say(
        f"{TAG} [{mode}] h_p move SAR vs P21 (value/P21 - 1): d2x "
        + "  ".join(f"{k}={d2['sar'][k] / P21_D2X_SAR[k] - 1:+.6e}" for k in PROBE_NAMES)
        + ";  d1x "
        + "  ".join(f"{k}={d1['sar'][k] / P21_D1X_SAR[k] - 1:+.6e}" for k in PROBE_NAMES)
    )
    _say(
        f"{TAG} [{mode}] h_p move theta40 vs P21: d2x "
        + "  ".join(f"{k}={t2['row4'][k] / P21_D2X_THETA40[k] - 1:+.6e}" for k in PROBE_NAMES)
        + ";  d1x "
        + "  ".join(f"{k}={t1['row4'][k] / P21_D1X_THETA40[k] - 1:+.6e}" for k in PROBE_NAMES)
    )
    _say(
        f"{TAG} [{mode}] h_p move row2 <theta>40 vs P21: d2x {t2['row2'][40] / P21_ROW2_D2X[40] - 1:+.6e}, "
        f"d1x {t1['row2'][40] / P21_ROW2_D1X[40] - 1:+.6e}"
    )
    _say(
        f"{TAG} [{mode}] h_p move D/C: d2x {d2['DC'] / P21_DC_D2X - 1:+.6e}, d1x {d1['DC'] / P21_DC_D1X - 1:+.6e};"
        f"  D/A: d2x {d2['DA'] / P21_DA_D2X - 1:+.6e}, d1x {d1['DA'] / P21_DA_D1X - 1:+.6e}"
    )
    if mode == "h0025":
        _say(
            f"{TAG} [h0025] three-way split, point SAR: sheet step d1x/ans9-d1 "
            + "  ".join(f"{k}={d1['sar'][k] / ANS9_D1_SAR[k]:.6f}" for k in PROBE_NAMES)
            + ";  order step d2x/d1x "
            + "  ".join(f"{k}={d2['sar'][k] / d1['sar'][k]:.6f}" for k in PROBE_NAMES)
        )
        _say(
            f"{TAG} [h0025] three-way split, theta40: sheet step d1x/ROW4 "
            + "  ".join(f"{k}={t1['row4'][k] / ROW4_POINT_K[k]:.6f}" for k in PROBE_NAMES)
            + ";  order step d2x/d1x "
            + "  ".join(f"{k}={t2['row4'][k] / t1['row4'][k]:.6f}" for k in PROBE_NAMES)
            + f";  row2 d1x/ans9 {t1['row2'][40] / ROW2_MEAN_K[40]:.6f}"
        )
    rho = d2["DC"] / ANS9_D1_DC
    _say(
        f"{TAG} [{mode}] rho(h) = (D/C)_d2x / {ANS9_D1_DC} = {rho:.6f} (P21: {P21_RHO}; printed only, "
        f"no branch read)"
    )
    if mode != "x1":
        live = _rel(s2, ANS4_EXACT_STRIP_S_P2_X1["10"])
        out["live"] = live
        _say(f"{TAG} [{mode}] (Live-h) ||S(h_p) - S_x1||/||S_x1|| = {live:.3e} (predicted > 1e-6)")
        _say(
            f"{TAG} [{mode}] (d6', printed per ruling Q5) {t2.get('ctl_d', float('nan')):.6f} "
            f"(predicted {PRED_D6:.3f})"
        )
    body = time.perf_counter() - t_body
    out["body"] = body
    rss = _peak_rss_gib()
    _say(
        f"{TAG} PRICE rung {mode} h_p {h_p} cells {cells} tag3 {tag3} E {n_e} F {n_f} d2dofs {out['dofs2']} "
        f"build {t_build:.2f} s d2sweep {t_s2:.2f} s d1sweep {t_s1:.2f} s thermal {t_th:.2f} s "
        f"body {body:.2f} s peakRSS {rss:.2f} GiB ranks {comm.size}"
    )
    return out


@complex_only
@needs_rung
def test_p1w_and_g2_and_reciprocity(rung):
    for lab in ("d2x", "d1x"):
        rd = rung[lab]
        assert all(rd["valid"]), (lab, rd["valid"])
        assert rd["p1w"] <= EXACT_IDENTITY_RTOL, (lab, rd["p1w"])
    assert len(rung["areas"]) == 4, rung["areas"]
    for k, a in rung["areas"].items():
        assert abs(a / NOMINAL_STRIP_AREA - 1.0) <= STRIP_BAND, (k, a)
    assert rung["recip2"] <= RECIPROCITY_BAND, rung["recip2"]
    assert rung["recip1"] <= RECIPROCITY_BAND, rung["recip1"]


@complex_only
@needs_rung
def test_thermal_anchors(rung):
    for lab in ("th_d2x", "th_d1x"):
        t = rung[lab]
        assert all(t["valid4"]), (lab, t["valid4"])
        assert t["r1"] <= B_BAND, (lab, t["r1"])
        for run, v in t["E"].items():
            assert v <= E_BAND, (lab, run, v)
        assert t["r3"] <= S2_BAND, (lab, t["r3"])


@complex_only
@needs_rung
def test_x1_route_tie(rung):
    if rung["mode"] != "x1":
        pytest.skip("(F0'), (X1c), (F0-ctl), (d6') asserted at x1 only")
    r2 = _rel(rung["s_d2x"], ANS4_EXACT_STRIP_S_P2_X1["10"])
    r1 = _rel(rung["s_d1x"], ANS4_EXACT_STRIP_S_P1_RECORDS["10"])
    ctl = _rel(rung["s_d2x"], ANS4_EXACT_STRIP_S_P1_RECORDS["10"])
    cell_rel = abs(rung["cells"] / ANS4_EXACT_STRIP_CELL_RECORD - 1.0)
    _say(
        f"\n{TAG} (F0') d2x vs P2 x1 record {r2:.3e}; d1x vs P1 record {r1:.3e} "
        f"(band {STEP3A_RECORD_RTOL:.0e}, predicted <= {PRED_F0:.0e}); (F0-ctl) {ctl:.3e} "
        f"(asserted > band, predicted {PRED_F0_CTL:.2e}); (X1c) cells {rung['cells']} vs "
        f"{ANS4_EXACT_STRIP_CELL_RECORD} rel {cell_rel:.3e} (band {CELL_COUNT_BAND}); "
        f"(d6') {rung['th_d2x']['ctl_d']:.6f} (asserted > {D_CONTROL_FLOOR:.0e})"
    )
    assert cell_rel <= CELL_COUNT_BAND, (rung["cells"], ANS4_EXACT_STRIP_CELL_RECORD)
    assert r2 <= STEP3A_RECORD_RTOL, r2
    assert r1 <= STEP3A_RECORD_RTOL, r1
    assert ctl > STEP3A_RECORD_RTOL, ctl
    assert rung["th_d2x"]["ctl_d"] > D_CONTROL_FLOOR, rung["th_d2x"]["ctl_d"]


@complex_only
@needs_rung
def test_h0025_reaches_i14_fixture(rung):
    if rung["mode"] != "h0025":
        pytest.skip("(R) is the h0025 window's")
    _say(
        f"\n{TAG} (R) cells {rung['cells']} vs I14_OFF_CELLS {I14_OFF_CELLS}; d2 dofs "
        f"{rung['dofs2']} vs I14_OFF_DOFS {I14_OFF_DOFS}"
    )
    assert abs(rung["cells"] / I14_OFF_CELLS - 1.0) <= CELL_COUNT_BAND, rung["cells"]
    assert rung["dofs2"] == I14_OFF_DOFS, rung["dofs2"]
