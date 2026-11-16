"""`ANS-9` step 1b (a) — local SAR at the SPEC monitor points A–E on two EM columns.

§9 item 371a (docs/planning/chunks/ANS-9.md, "§9 item 371a"; §10 chain P21 (a)).
Our values only: no AED number is read or written here.

Columns (opt-in knob ``FEM_EM_ANS9_1B_COLUMN`` ∈ {unset, ``d1``, ``d2x``}):

* ``d1``  — ``ans:9``'s own EM solve, ``_run_heated(return_fields=True)``
  (ANS-2's 10 MHz ragged-sheet fixture, degree 1), drive k = 0 (P1), G = 1/P₁.
* ``d2x`` — item 192's exact-rectangle rung at degree 2, κ off, P1's field from
  ``_solve_driven``; plus ``d1x`` (same mesh, degree 1, printed decomposition).

Point SAR = G·σ_saline·|E|²/(2ρ), the N1curl field evaluated on the parent mesh
with ``owner="lowest_original_cell"``.  Anchors: (P1W) independent ∫ G·½σ|E|²
over the phantom = 1 W within ``EXACT_IDENTITY_RTOL``; (ORD, d1) SAR_B > SAR_D >
SAR_C; (F0) the rung's 10 MHz S against item 192's records within
``STEP3A_RECORD_RTOL`` and (F0-ctl) the d2x S misses the degree-1 record by more
than that band.  Printed: SAR at A–E, D/C, D/A, and in W2 the weekly's decision
ratio ρ = (D/C)_d2x / (D/C)_d1 with its branch (pre-registered, unchanged).
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.sar import mean_sar

from tests.complex_mode import complex_only
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_S_P1_RECORDS,
    ANS4_EXACT_STRIP_S_P2_X1,
)
from tests.validation.test_ans4_resolution_ladder import STEP3A_RECORD_RTOL
from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3, _solve_driven
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_lossy_sphere_fullwave import SALINE_SIGMA
from dolfinx import default_scalar_type, fem

from tests.validation.test_mat10_sar_to_heat_sphere import (
    B_BAND,
    ROBIN_TAG,
    _allsum,
    _robin_tags,
    _say,
    _to_sub,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND, RHO_C, S2_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    D_CONTROL_FLOOR,
    DT_1,
    N_MAIN,
    T_INF,
    T_STAR,
    _be_one,
    _cell_average_dg0_tag,
    _e_prime,
    _run_heated,
    _tag_submesh,
    _transient,
)
from tests.validation.test_mat19_probe_coil_phantom import (
    ROW2_MEAN_K,
    ROW4_POINT_K,
    _load_spec_probes,
)
from tests.validation.test_port_birdcage_leg_offset_sweep import _four_port_rung
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_gap_voltage_impedance import (
    FREQUENCY_HZ,
    SIGMA_WIRE_S_PER_M,
)
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

COLUMN_ENV = "FEM_EM_ANS9_1B_COLUMN"
COLUMNS = ("d1", "d2x")
D1_JSON = Path("/workspace/logs/ans9_1b/d1.json")
ANS2_METRICS = (
    Path(__file__).resolve().parents[2]
    / "examples/ansys_benchmarks/ans2_birdcage_coil_driven_sar_10MHz/metrics.json"
)
PROBE_NAMES = ("A", "B", "C", "D", "E")

# The weekly's pre-registered rule (2026-10-10, §10 P21), "within 1.5×" read
# by the 2026-10-11 03:00 review as 1/1.5 ≤ ρ ≤ 1.5.
BRANCH_RAGGED = "degree-1 ragged source"
BRANCH_FAR_LEG = "far-leg current model"
BRANCH_NONE = "no branch"


def _column():
    raw = os.environ.get(COLUMN_ENV, "").strip()
    if not raw:
        return None
    if raw not in COLUMNS:
        raise ValueError(f"{COLUMN_ENV}={raw!r}: only unset, 'd1' or 'd2x'")
    return raw


def _branch(rho):
    if rho >= 3.0:
        return BRANCH_RAGGED
    if 1.0 / 1.5 <= rho <= 1.5:
        return BRANCH_FAR_LEG
    return BRANCH_NONE


needs_column = pytest.mark.skipif(_column() is None, reason=f"{COLUMN_ENV} unset (opt-in window)")


def _point_sar(e, pts, G, sigma):
    vals, ok = evaluate_vector_field_parallel(e, pts, owner="lowest_original_cell")
    e2 = np.sum(np.abs(np.asarray(vals)) ** 2, axis=1)
    return G * sigma * e2 / (2.0 * PHANTOM_RHO_KG_PER_M3), np.asarray(ok, dtype=bool)


def _column_readout(label, e, sigma_field, cell_tags, comm):
    """G, (P1W) and the point SAR at A–E for one solved P1 field."""
    P = float(mean_sar(e, sigma=sigma_field, rho=PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
                       comm=comm, subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"])
    G = 1.0 / P
    dxp = ufl.Measure("dx", domain=e.function_space.mesh, subdomain_data=cell_tags)(PHANTOM_CELL_TAG)
    p_vol = float(np.real(_allsum(comm, G * 0.5 * SALINE_SIGMA * ufl.inner(e, e) * dxp)))
    p1w = abs(p_vol - 1.0)
    spec = _load_spec_probes()
    pts = np.array([[float(v) for v in spec[k]] for k in PROBE_NAMES], dtype=np.float64)
    sar, ok = _point_sar(e, pts, G, SALINE_SIGMA)
    s = dict(zip(PROBE_NAMES, (float(x) for x in sar)))
    _say(f"\nANS9-1b [{label}] P_1 = {P:.12e} W, G = {G:.12e} 1/W; (P1W) int_phantom G*sigma|E|^2/2 "
         f"= {p_vol:.15f} W, |.-1| = {p1w:.3e} (band {EXACT_IDENTITY_RTOL:.0e})")
    for k in PROBE_NAMES:
        _say(f"ANS9-1b [{label}] SAR_{k} at {tuple(pts[PROBE_NAMES.index(k)])} = {s[k]:.9e} W/kg "
             f"(1 W in phantom; valid {bool(ok[PROBE_NAMES.index(k)])})")
    dc, da = s["D"] / s["C"], s["D"] / s["A"]
    _say(f"ANS9-1b [{label}] D/C = {dc:.6f}  D/A = {da:.6f}")
    return {"P": P, "G": G, "p1w": p1w, "sar": s, "valid": [bool(v) for v in ok],
            "DC": dc, "DA": da, "pts": pts}


def _thermal_readout(label, msh, cell_tags, e, sigma_field, G, pts, sub, parent_ids, ft, comm,
                     control=False):
    """§9 item 371b: ``ans:9``'s thermal half (step 1a schedule) on one P1 source."""
    Wp, q_par = _cell_average_dg0_tag(msh, cell_tags, PHANTOM_CELL_TAG,
                                      G * 0.5 * sigma_field * ufl.inner(e, e))
    q_h = _to_sub(Wp, q_par, sub, parent_ids)
    one = fem.Constant(sub, default_scalar_type(1.0))
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    area = float(np.real(_allsum(comm, one * ds)))
    int_q = float(np.real(_allsum(comm, q_h * ufl.dx(domain=sub))))
    r1 = abs(int_q - 1.0)
    q_s = float(np.real(_allsum(comm, q_h * ds)))
    kappa = (q_s / area) / (int_q / v_sub)
    t1 = time.perf_counter()
    res = _transient(sub, ft, q_h, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
                     startup_half_steps=2, snapshot_steps=range(N_MAIN + 1))
    t_cn = time.perf_counter() - t1
    mean = {n: float(np.real(_allsum(comm, (res.snapshots[n] - T_INF) * ufl.dx(domain=sub)))) / v_sub
            for n in (10, 20, 30, 40)}
    ep, stor, conv = _e_prime(sub, ds, res, 1.0)
    res1 = _be_one(sub, ft, q_h)
    thbar1 = float(np.real(_allsum(comm, (res1.temperature - T_INF) * ufl.dx(domain=sub)))) / v_sub
    ratio = RHO_C * thbar1 / (DT_1 * int_q / v_sub)
    r3 = abs(ratio - 1.0)
    vals, okT = evaluate_vector_field_parallel(res.snapshots[N_MAIN], pts)
    th = dict(zip(PROBE_NAMES, (float(x) for x in np.real(np.asarray(vals)[:, 0]) - T_INF)))
    E = {"cn": res.max_balance_residual, "be": res1.max_balance_residual}
    _say(f"ANS9-1b-th [{label}] (R1') int_sub q_h = {int_q:.15f} W, |.-1| = {r1:.3e} "
         f"(bound {B_BAND:.0e}); V_sub = {v_sub:.9e} m^3, oint ds = {area:.6e} m^2, kappa_s = {kappa:.4f}")
    _say(f"ANS9-1b-th [{label}] row 2 <theta> [K]: "
         + "  ".join(f"n={n} {mean[n]:.6f} (ans:9 d1 {ROW2_MEAN_K[n]:.6f})" for n in (10, 20, 30, 40))
         + f";  CN solve {t_cn:.1f} s")
    _say(f"ANS9-1b-th [{label}] (E) max_balance CN {E['cn']:.3e}, BE {E['be']:.3e} (bound {E_BAND:.0e}); "
         f"(E', printed) stored {stor:.9e} + surface {conv:.9e} vs 1 W: rel {ep:.3e}")
    _say(f"ANS9-1b-th [{label}] row 3 (R3') rho c theta1/(dt_1 int q/V) = {ratio:.12f}, |.-1| = {r3:.4e} "
         f"(bound {S2_BAND:.0e}; predicted 5e-6..2e-5)")
    _say(f"ANS9-1b-th [{label}] row 4 theta40 at A-E [K]: "
         + "  ".join(f"{k}={th[k]:.6f} (ans:9 d1 {ROW4_POINT_K[k]:.6f})" for k in PROBE_NAMES)
         + f"; valid {[bool(v) for v in okT]}; D/C = {th['D'] / th['C']:.6f}  D/A = {th['D'] / th['A']:.6f}")
    out = {"r1": r1, "r3": r3, "E": E, "row2": mean, "row4": th, "valid4": [bool(v) for v in okT]}
    if control:
        q_two = fem.Function(q_h.function_space)
        q_two.x.array[:] = 2.0 * q_h.x.array
        q_two.x.scatter_forward()
        res_d = _transient(sub, ft, q_two, dt=T_STAR / N_MAIN, n_steps=N_MAIN, theta=0.5,
                           startup_half_steps=2, snapshot_steps=(N_MAIN - 1, N_MAIN))
        ep_d, _, _ = _e_prime(sub, ds, res_d, 1.0)
        out["ctl_d"] = ep_d
        out["E"]["ctl"] = res_d.max_balance_residual
        _say(f"ANS9-1b-th [{label}] control (d6') doubled source: (E') vs 1 W = {ep_d:.6f} "
             f"(predicted 1.000, asserted > {D_CONTROL_FLOOR:.0e}); its (E) {res_d.max_balance_residual:.3e}")
    return out


@pytest.fixture(scope="module")
def column():
    col = _column()
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    res = {"column": col}
    if col == "d1":
        out = _run_heated(return_fields=True)
        f = out["solves"][out["by_k"][0]]["fields"]
        rd = _column_readout("d1", f.e_complex, f.sigma_field, out["cell_tags"], comm)
        res["d1"] = rd
        # (Blind-sigma), printed only.
        blind, _ = _point_sar(f.e_complex, rd["pts"], rd["G"], SIGMA_WIRE_S_PER_M)
        _say("ANS9-1b [d1] (Blind-sigma) SAR with sigma_wire at A-E: "
             + "  ".join(f"{k}={b:.4e}" for k, b in zip(PROBE_NAMES, blind))
             + f"  ratio to saline {blind[0] / rd['sar']['A']:.6e} "
             f"(sigma ratio {SIGMA_WIRE_S_PER_M / SALINE_SIGMA:.6e})")
        # ANS-2's r = 0.015 m named points, our N1curl route vs its CG1 metrics.json.
        m = json.loads(ANS2_METRICS.read_text())
        npts = np.asarray(m["named_points_m"], dtype=np.float64)
        ours, _ = _point_sar(f.e_complex, npts, rd["G"], SALINE_SIGMA)
        p_ans2 = float(m["phantom_power"]["0"]["dissipated_power_w"])
        theirs = np.asarray(m["pointwise"]["0"]["sar_w_per_kg"], dtype=float) / p_ans2
        for p, a, b in zip(npts, ours, theirs):
            _say(f"ANS9-1b [d1] ANS-2 named point {tuple(p)}: N1curl owner-rule {a:.6e} W/kg vs "
                 f"ANS-2 metrics.json ({m.get('pointwise_rule', 'cg1')}) {b:.6e} W/kg (both 1 W; printed)")
        # θ40 at A–E through the default owner on the submesh (item 218's route).
        T = out["T_star"][0]
        vals, okT = evaluate_vector_field_parallel(T, rd["pts"])
        th = np.real(np.asarray(vals)[:, 0]) - T_INF
        _say("ANS9-1b [d1] theta40 at A-E [K]: "
             + "  ".join(f"{k}={t:.6f} (rec {ROW4_POINT_K[k]:.6f})" for k, t in zip(PROBE_NAMES, th))
             + f"; valid {[bool(v) for v in okT]}; row 2 mean rec at 40 min {ROW2_MEAN_K[40]}")
        if comm.rank == 0:
            D1_JSON.parent.mkdir(parents=True, exist_ok=True)
            D1_JSON.write_text(json.dumps({"sar": rd["sar"], "DC": rd["DC"], "DA": rd["DA"],
                                           "G": rd["G"]}, indent=2))
    elif col == "d2x":
        rung = _four_port_rung("ANS-9 1b exact 10 MHz", np.zeros(4), frequency_hz=FREQUENCY_HZ,
                               degree=2, leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
        res["s_d2x"] = np.asarray(rung["s"], dtype=np.complex128)
        _say(f"\nANS9-1b [d2x] rung cells {int(rung['cells'])} (record 117 206), degree {rung['degree']}")
        f2 = _solve_driven(rung, "P1")["fields"]
        res["d2x"] = _column_readout("d2x", f2.e_complex, f2.sigma_field, rung["cell_tags"], comm)
        r1 = _four_port_rung("ANS-9 1b exact 10 MHz d1x", np.zeros(4), frequency_hz=FREQUENCY_HZ,
                             degree=1, reuse=rung)
        res["s_d1x"] = np.asarray(r1["s"], dtype=np.complex128)
        f1 = _solve_driven(r1, "P1")["fields"]
        res["d1x"] = _column_readout("d1x", f1.e_complex, f1.sigma_field, r1["cell_tags"], comm)
        # §9 item 371b: the thermal half on both sources, one tag-3 submesh of d2x's mesh.
        t_th = time.perf_counter()
        msh, ct = rung["mesh"], rung["cell_tags"]
        sub, parent_ids = _tag_submesh(msh, ct, PHANTOM_CELL_TAG)
        ft = _robin_tags(sub)
        _say(f"ANS9-1b-th tag-{PHANTOM_CELL_TAG} submesh cells {sub.topology.index_map(3).size_global} "
             f"(predicted 500-1500)")
        res["th_d2x"] = _thermal_readout("d2x", msh, ct, f2.e_complex, f2.sigma_field,
                                         res["d2x"]["G"], res["d2x"]["pts"], sub, parent_ids, ft,
                                         comm, control=True)
        res["th_d1x"] = _thermal_readout("d1x", msh, ct, f1.e_complex, f1.sigma_field,
                                         res["d1x"]["G"], res["d1x"]["pts"], sub, parent_ids, ft, comm)
        _say("ANS9-1b-th order effect theta40 d2x/d1x (same thermal mesh): "
             + "  ".join(f"{k}={res['th_d2x']['row4'][k] / res['th_d1x']['row4'][k]:.6f}"
                         for k in PROBE_NAMES)
             + "; <theta>40 " + f"{res['th_d2x']['row2'][40] / res['th_d1x']['row2'][40]:.6f}"
             + f"; thermal {time.perf_counter() - t_th:.1f} s")
        if D1_JSON.exists():
            d1 = json.loads(D1_JSON.read_text())
            rho = res["d2x"]["DC"] / float(d1["DC"])
            res["rho"] = rho
            _say(f"ANS9-1b decision: (D/C)_d1 = {float(d1['DC']):.6f}  (D/C)_d1x = "
                 f"{res['d1x']['DC']:.6f}  (D/C)_d2x = {res['d2x']['DC']:.6f};  rho = "
                 f"(D/C)_d2x/(D/C)_d1 = {rho:.6f} => branch: {_branch(rho)}")
        else:
            _say("ANS9-1b decision: W1 JSON absent")
    _say(f"ANS9-1b column {col} elapsed {time.perf_counter() - t0:.1f} s")
    return res


def test_ans9_1b_coding(monkeypatch):
    """W0: knob parser, the SPEC's A–E, and the pre-registered branch function."""
    for raw, want in (("", None), ("d1", "d1"), (" d2x ", "d2x")):
        monkeypatch.setenv(COLUMN_ENV, raw)
        assert _column() == want
    for bad in ("d1x", "D1", "2"):
        monkeypatch.setenv(COLUMN_ENV, bad)
        with pytest.raises(ValueError):
            _column()
    spec = _load_spec_probes()
    assert tuple(spec) == PROBE_NAMES
    assert np.allclose(np.asarray(spec["D"], dtype=float), (-0.020, 0.0, 0.0), rtol=0, atol=0)
    assert np.allclose(np.asarray(spec["C"], dtype=float), (0.0, 0.020, 0.0), rtol=0, atol=0)
    assert _branch(3.0) == BRANCH_RAGGED
    assert _branch(1.0) == BRANCH_FAR_LEG
    assert _branch(1.0 / 1.5) == BRANCH_FAR_LEG
    assert _branch(2.0) == BRANCH_NONE
    assert _branch(0.5) == BRANCH_NONE


@complex_only
@needs_column
def test_p1w_phantom_power_is_one_watt(column):
    for lab in ("d1", "d2x", "d1x"):
        if lab in column:
            rd = column[lab]
            assert all(rd["valid"]), f"{lab}: invalid point evaluation {rd['valid']}"
            assert rd["p1w"] <= EXACT_IDENTITY_RTOL, f"{lab}: (P1W) {rd['p1w']:.3e}"


@complex_only
@needs_column
def test_ord_d1_far_side_between_driven_and_adjacent(column):
    if column["column"] != "d1":
        pytest.skip("(ORD) is the d1 column's")
    s = column["d1"]["sar"]
    assert s["B"] > s["D"] > s["C"], s


def _rel(s, rec):
    rec = np.asarray(rec, dtype=np.complex128)
    return float(np.linalg.norm(s - rec) / np.linalg.norm(rec))


@complex_only
@needs_column
def test_f0_fixture_is_item_192(column):
    if column["column"] != "d2x":
        pytest.skip("(F0) is the d2x window's")
    r2 = _rel(column["s_d2x"], ANS4_EXACT_STRIP_S_P2_X1["10"])
    r1 = _rel(column["s_d1x"], ANS4_EXACT_STRIP_S_P1_RECORDS["10"])
    ctl = _rel(column["s_d2x"], ANS4_EXACT_STRIP_S_P1_RECORDS["10"])
    _say(f"\nANS9-1b (F0) d2x vs P2 x1 record {r2:.3e}; d1x vs P1 record {r1:.3e} "
         f"(band {STEP3A_RECORD_RTOL:.0e}); (F0-ctl) d2x vs P1 record {ctl:.3e} "
         f"(asserted > band, predicted ≈ 2e-2)")
    assert r2 <= STEP3A_RECORD_RTOL
    assert r1 <= STEP3A_RECORD_RTOL
    assert ctl > STEP3A_RECORD_RTOL


@complex_only
@needs_column
def test_371b_thermal_half_on_degree2_exact_source(column):
    """§9 item 371b anchors: (R1'), (E), (R3') per column; (d6') doubled-source control."""
    if column["column"] != "d2x":
        pytest.skip("371b's thermal half is the d2x window's")
    for lab in ("th_d2x", "th_d1x"):
        t = column[lab]
        assert all(t["valid4"]), f"{lab}: invalid row-4 evaluation {t['valid4']}"
        assert t["r1"] <= B_BAND, f"{lab}: (R1') {t['r1']:.3e}"
        for run, v in t["E"].items():
            assert v <= E_BAND, f"{lab}: (E) {run} {v:.3e}"
        assert t["r3"] <= S2_BAND, f"{lab}: (R3') {t['r3']:.3e}"
    assert column["th_d2x"]["ctl_d"] > D_CONTROL_FLOOR, column["th_d2x"]["ctl_d"]
