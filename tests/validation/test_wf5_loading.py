"""`WF-5` step 2 (§9 item 142, chain T15): loading — phantom σ 0.5 → 0 at fixed C.

The tuned F-small birdcage on `ANS-7`'s Route A (step 1's fixture, tracked κ)
is solved twice per column: loaded (σ = ``th14.SALINE_SIGMA``) and unloaded
(σ = 0, ε_r = 78 kept).  By first-order perturbation removing σ moves ω̃ along
the imaginary axis, so Δf₀ is second order and is printed only; the asserted
content is the loss partition:

* (A0) the loaded secant f₀ / Q reproduce the tracked ``metrics.json`` rows;
* (A1) P_phantom,U = 0.0 exactly (the knob reached the solve);
* (A2) unloaded |P_loss,U − (p_src − p_sheet(P1))| ≤ ``DISCRETE_IDENTITY_RTOL``·p_src;
* (A3) Cu unloaded Q two ways inside B_Q,U (frozen from W1);
* (A4) |Δ(1/Q)_slope − Π|/Π ≤ B_Δ, Π = 1/Q_field,L − 1/Q_field,U, B_Δ by the
  triangle inequality from (A3) and step 1's B_Q.

Negative control (binary): PEC σ = 0 is lossless, 1/Q_slope,U ≤ 1e-6.
Never divides by G on the PEC unloaded column — 1/Q = 2G/(f₀ dB/df) throughout.
Tests only.  Degree 1, complex mode.
"""

from __future__ import annotations

import contextlib
import json

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.validation.test_wf5_q_two_ways import (
    ANS7_DIR,
    ANS7_SCRIPT,
    B_Q_FROZEN,
    BAND_FLOOR,
    COLUMNS,
    CUBIC_MHZ,
    RECORD_RTOL,
    _budget,
    _load,
    _y_in,
    ceil_1sf,
)

UNLOADED_MHZ = tuple(round(63.7 + 0.1 * k, 1) for k in range(7))  # 63.7 … 64.3
CONTROL_INV_Q_MAX = 1.0e-6
# *Predicted, printed only* (§9 item 142).
PREDICTED_E_SEC_U = 1.0e-3
PREDICTED_DF0_ABS = 1.0e-3
PREDICTED_PI_PH_REL = 2.0e-2
STEP1_PI_PH = {
    "cu": 6.489692e-04,
    "pec": 6.493197e-04,
}  # 20261001T141511Z_WF-5-step1-W2.log:971,978
# Frozen from W1 (``20261001T184835Z_WF-5-step2-W1.log:962``): Cu e_sec,U
# 5.7992e-4 ⇒ ceil_1sf(3 × 5.7992e-4) = 2e-3 (predicted 3e-3).  Never loosened.
B_Q_U_FROZEN = 2.0e-3


@contextlib.contextmanager
def _unloaded(th14):
    """σ_phantom → 0 on every rank; ``_problem`` reads the global at call time."""
    saved = th14.SALINE_SIGMA
    th14.SALINE_SIGMA = 0.0
    try:
        yield
    finally:
        th14.SALINE_SIGMA = saved


@pytest.fixture(scope="module")
def ans7():
    return _load(ANS7_SCRIPT, "ans7_driven_sweep_wf5")


@pytest.fixture(scope="module")
def fixture_build(ans7):
    """`test_wf5_q_two_ways.py:99–117` verbatim (κ from the tracked metrics.json)."""
    from tests.validation import test_th14_birdcage_copper as th14

    comm = MPI.COMM_WORLD
    tracked = json.loads((ANS7_DIR / "metrics.json").read_text())
    kappa = float(tracked["kappa"])
    msh, cell_tags, facet_tags, _diag, _t_mesh = th14._sheets_build(True, as_hole=True)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    outer_tags, _census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_k, port_defs_k, specs_k = th14._build_ports(
        msh, cell_tags, comm, width_correction_kappa=kappa
    )
    geo_k = (msh, cell_tags, facet_tags, outer_tags, tags_k, port_defs_k)
    if comm.rank == 0:
        print(
            f"[wf5-2] hole mesh {int(msh.topology.index_map(tdim).size_global)} cells; "
            f"kappa (tracked) = {kappa:.9e}; SALINE_SIGMA (loaded) = {th14.SALINE_SIGMA}",
            flush=True,
        )
    return {"th14": th14, "geo": geo_k, "specs": specs_k, "tracked": tracked, "comm": comm}


def _sweep(ans7, fb, grid_mhz, col):
    y = []
    for f_mhz in grid_mhz:
        s11, _ = ans7._route_a_solve(fb["th14"], fb["geo"], fb["specs"], f_mhz * 1e6, col)
        y.append(_y_in(ans7, s11, f_mhz * 1e6))
    return np.array(y)


@pytest.fixture(scope="module")
def sweep_unloaded(ans7, fixture_build):
    th14 = fixture_build["th14"]
    out = {}
    with _unloaded(th14):
        for col in COLUMNS:
            out[col] = _sweep(ans7, fixture_build, UNLOADED_MHZ, col)
    return out


@pytest.fixture(scope="module")
def sweep_loaded(ans7, fixture_build):
    return {col: _sweep(ans7, fixture_build, CUBIC_MHZ, col) for col in COLUMNS}


def _rows_on_grid(f_mhz, y):
    """The − → + crossing of Im Y_in on a grid: secant on the bracketing pair and
    a cubic through the 4 points centred on it.  Returns 1/Q = 2G/(f₀ dB/df),
    never Q (G may be ≈ 0 on a lossless column)."""
    f_mhz = np.asarray(f_mhz, dtype=float)
    y = np.asarray(y)
    b, g = y.imag, y.real
    up = [k for k in range(len(b) - 1) if b[k] < 0 <= b[k + 1]]
    assert (
        len(up) == 1
    ), f"expected exactly one - -> + crossing of Im Y on {list(f_mhz)} MHz, got {up}"
    k = up[0]
    f_hz = f_mhz * 1e6
    t = -b[k] / (b[k + 1] - b[k])
    f0 = f_hz[k] + t * (f_hz[k + 1] - f_hz[k])
    g0 = g[k] + t * (g[k + 1] - g[k])
    dbdf = (b[k + 1] - b[k]) / (f_hz[k + 1] - f_hz[k])
    inv_q_sec = 2.0 * g0 / (f0 * dbdf)
    lo = min(max(k - 1, 0), len(b) - 4)
    idx = slice(lo, lo + 4)
    fc = f_mhz[k]
    x = f_mhz[idx] - fc
    pb = np.polyfit(x, b[idx], 3)
    pg = np.polyfit(x, g[idx], 3)
    xk1 = f_mhz[k + 1] - fc
    roots = [
        r.real for r in np.roots(pb) if abs(r.imag) < 1e-12 and -1e-12 <= r.real <= xk1 + 1e-12
    ]
    assert len(roots) == 1, f"cubic B: expected one root in the bracket, got {roots}"
    x0 = roots[0]
    f0_c = (fc + x0) * 1e6
    dbdf_c = np.polyval(np.polyder(pb), x0) / 1.0e6
    inv_q_cub = 2.0 * np.polyval(pg, x0) / (f0_c * dbdf_c)
    return {
        "f0_hz": float(f0),
        "g0": float(g0),
        "dbdf": float(dbdf),
        "inv_q_sec": float(inv_q_sec),
        "f0_cubic_hz": float(f0_c),
        "inv_q_cubic": float(inv_q_cub),
        "bracket_mhz": (float(f_mhz[k]), float(f_mhz[k + 1])),
    }


def _e_sec(rows):
    """|Q_sec − Q_cub|/Q_cub in the 1/Q form."""
    return abs(rows["inv_q_cubic"] - rows["inv_q_sec"]) / rows["inv_q_sec"]


@complex_only
def test_unloaded_band_candidate(fixture_build, sweep_unloaded):
    """B_Q,U from the Cu unloaded secant's error against the cubic; PEC σ=0 control."""
    comm = fixture_build["comm"]
    rank0 = comm.rank == 0
    rows = {col: _rows_on_grid(UNLOADED_MHZ, sweep_unloaded[col]) for col in COLUMNS}
    for col in COLUMNS:
        r = rows[col]
        if rank0:
            print(
                f"[wf5-2] unloaded {col}: bracket {r['bracket_mhz']} MHz; secant f0 "
                f"{r['f0_hz']:.6f} Hz, G0 {r['g0']:.9e} S, dB/df {r['dbdf']:.9e} S/Hz, "
                f"1/Q_sec {r['inv_q_sec']:.9e}; cubic f0 {r['f0_cubic_hz']:.6f} Hz, "
                f"1/Q_cub {r['inv_q_cubic']:.9e}",
                flush=True,
            )
    e_sec_u = _e_sec(rows["cu"])
    b_q_u = ceil_1sf(max(3.0 * e_sec_u, BAND_FLOOR))
    ctrl = rows["pec"]["inv_q_sec"]
    if rank0:
        print(
            f"[wf5-2] cu unloaded Q_sec {1.0 / rows['cu']['inv_q_sec']:.6f}, Q_cub "
            f"{1.0 / rows['cu']['inv_q_cubic']:.6f}; e_sec,U {e_sec_u:.4e} (PREDICTED "
            f"{PREDICTED_E_SEC_U:.0e}); B_Q,U candidate = ceil_1sf(max(3*{e_sec_u:.4e}, "
            f"{BAND_FLOOR:.0e})) = {b_q_u:.1e}; frozen B_Q_U_FROZEN = {B_Q_U_FROZEN}",
            flush=True,
        )
        print(
            f"[wf5-2] control: PEC sigma=0 1/Q_slope,U = {ctrl:+.6e} (ASSERTED |.| <= "
            f"{CONTROL_INV_Q_MAX:.0e}; PREDICTED <= 1e-9; dead knob reads 6.50e-4)",
            flush=True,
        )
    assert abs(ctrl) <= CONTROL_INV_Q_MAX, f"PEC sigma=0 control: 1/Q_slope,U = {ctrl}"
    if B_Q_U_FROZEN is not None:
        assert b_q_u <= B_Q_U_FROZEN, f"B_Q,U candidate {b_q_u} exceeds the frozen {B_Q_U_FROZEN}"


def _w_tot(b):
    return b["w_m"] + b["w_e"] + sum(b["w_sheet"].values()) + b["w_c1"] + b["w_surf"]


@complex_only
def test_loading_partition(ans7, fixture_build, sweep_loaded, sweep_unloaded):
    """A0–A4 and the control; Δf₀, Π_ph, shares printed."""
    assert B_Q_U_FROZEN is not None, "B_Q,U not frozen yet (run W1 first)"
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL

    fb = fixture_build
    th14 = fb["th14"]
    comm = fb["comm"]
    rank0 = comm.rank == 0
    tracked = fb["tracked"]["rows_route_a"]
    fails = []
    for col in COLUMNS:
        rl = _rows_on_grid(CUBIC_MHZ, sweep_loaded[col])
        ru = _rows_on_grid(UNLOADED_MHZ, sweep_unloaded[col])
        # (A0)
        q_l_slope = 1.0 / rl["inv_q_sec"]
        df0_rec = abs(rl["f0_hz"] - tracked[col]["f0_hz"]) / tracked[col]["f0_hz"]
        dq_rec = abs(q_l_slope - tracked[col]["q_loaded"]) / tracked[col]["q_loaded"]
        bl = _budget(ans7, fb, rl["f0_hz"], col)
        with _unloaded(th14):
            bu = _budget(ans7, fb, ru["f0_hz"], col)
        wl, wu = _w_tot(bl), _w_tot(bu)
        pl = bl["p_phantom"] + bl["p_surf"]
        pu = bu["p_phantom"] + bu["p_surf"]
        inv_qf_l = pl / (bl["omega"] * wl)
        inv_qf_u = pu / (bu["omega"] * wu)
        # (A2)
        id_u = abs(pu - (bu["p_src"] - bu["p_sheet1"])) / bu["p_src"]
        # (A3), Cu only (PEC unloaded Q is infinite)
        miss_u = (
            abs(1.0 / inv_qf_u - 1.0 / ru["inv_q_sec"]) * ru["inv_q_sec"] if col == "cu" else None
        )
        # (A4): |1/Q_s − 1/Q_f| = (|Q_f − Q_s|/Q_s)/Q_f ≤ B·(1/Q_f) per state.
        pi = inv_qf_l - inv_qf_u
        d_slope = rl["inv_q_sec"] - ru["inv_q_sec"]
        b_delta = (B_Q_FROZEN * inv_qf_l + B_Q_U_FROZEN * inv_qf_u) / pi
        miss_d = abs(d_slope - pi) / pi
        pi_ph = bl["p_phantom"] / (bl["omega"] * wl)
        df0 = (ru["f0_hz"] - rl["f0_hz"]) / rl["f0_hz"]
        if rank0:
            print(
                f"[wf5-2] {col}: (A0) loaded secant f0 {rl['f0_hz']:.6f} Hz (record "
                f"{tracked[col]['f0_hz']:.6f}, rel {df0_rec:.3e}); Q_slope,L {q_l_slope:.10f} "
                f"(record {tracked[col]['q_loaded']:.10f}, rel {dq_rec:.3e}; ASSERTED <= "
                f"{RECORD_RTOL:.0e})",
                flush=True,
            )
            print(
                f"[wf5-2] {col}: unloaded f0 {ru['f0_hz']:.6f} Hz; Delta f0/f0 = {df0:+.6e} "
                f"(PREDICTED |.| <~ {PREDICTED_DF0_ABS:.0e}, sign not predicted)",
                flush=True,
            )
            print(
                f"[wf5-2] {col}: (A1) P_phantom,U = {bu['p_phantom']!r} W (ASSERTED == 0.0); "
                f"P_surf,U {bu['p_surf']:.9e} W; p_src,U {bu['p_src']:.9e} W, p_sheet(P1),U "
                f"{bu['p_sheet1']:.9e} W; (A2) |P_loss,U - (p_src - p_sheet)|/p_src = {id_u:.3e} "
                f"(ASSERTED <= {DISCRETE_IDENTITY_RTOL:.0e})",
                flush=True,
            )
            print(
                f"[wf5-2] {col}: W_tot L {wl:.9e} J, U {wu:.9e} J; shares L: W_m "
                f"{bl['w_m'] / wl:.6f}, W_e {bl['w_e'] / wl:.6f}, W_sheets "
                f"{sum(bl['w_sheet'].values()) / wl:.6f}, W_C1 {bl['w_c1'] / wl:.6f}, W_surf "
                f"{bl['w_surf'] / wl:.4e}; shares U: W_m {bu['w_m'] / wu:.6f}, W_e "
                f"{bu['w_e'] / wu:.6f}, W_sheets {sum(bu['w_sheet'].values()) / wu:.6f}, W_C1 "
                f"{bu['w_c1'] / wu:.6f}, W_surf {bu['w_surf'] / wu:.4e}",
                flush=True,
            )
            print(
                f"[wf5-2] {col}: surface 1/Q share L {bl['p_surf'] / (bl['omega'] * wl):.6e}, "
                f"U {bu['p_surf'] / (bu['omega'] * wu):.6e}; 1/Q_field L {inv_qf_l:.9e}, U "
                f"{inv_qf_u:.9e}; 1/Q_slope L {rl['inv_q_sec']:.9e}, U {ru['inv_q_sec']:.9e}",
                flush=True,
            )
            if col == "cu":
                print(
                    f"[wf5-2] cu: Q_U slope {1.0 / ru['inv_q_sec']:.6f}, field "
                    f"{1.0 / inv_qf_u:.6f} (PREDICTED ~ 2.3e3); (A3) |Q_field,U - Q_slope,U|"
                    f"/Q_slope,U = {miss_u:.4e} (ASSERTED <= B_Q,U {B_Q_U_FROZEN:.0e})",
                    flush=True,
                )
            print(
                f"[wf5-2] {col}: (A4) Delta(1/Q)_slope {d_slope:.9e}, Pi {pi:.9e}; "
                f"|Delta - Pi|/Pi = {miss_d:.4e} (ASSERTED <= B_Delta {b_delta:.4e})",
                flush=True,
            )
            print(
                f"[wf5-2] {col}: Pi_ph = P_ph,L/(w0 W_L) {pi_ph:.6e} (step 1 "
                f"{STEP1_PI_PH[col]:.6e}); |Delta - Pi_ph|/Pi_ph = "
                f"{abs(d_slope - pi_ph) / pi_ph:.4e} (PREDICTED <= {PREDICTED_PI_PH_REL:.0e})",
                flush=True,
            )
            if col == "pec":
                print(
                    f"[wf5-2] control: PEC sigma=0 1/Q_slope,U = {ru['inv_q_sec']:+.6e} "
                    f"(ASSERTED |.| <= {CONTROL_INV_Q_MAX:.0e})",
                    flush=True,
                )
        if not (df0_rec <= RECORD_RTOL and dq_rec <= RECORD_RTOL):
            fails.append((col, "A0", df0_rec, dq_rec))
        if not bu["p_phantom"] == 0.0:
            fails.append((col, "A1", bu["p_phantom"]))
        if not id_u <= DISCRETE_IDENTITY_RTOL:
            fails.append((col, "A2", id_u))
        if col == "cu" and not miss_u <= B_Q_U_FROZEN:
            fails.append((col, "A3", miss_u))
        if not miss_d <= b_delta:
            fails.append((col, "A4", miss_d, b_delta))
        if col == "pec" and not abs(ru["inv_q_sec"]) <= CONTROL_INV_Q_MAX:
            fails.append((col, "control", ru["inv_q_sec"]))
    assert not fails, f"WF-5 step 2 gates: {fails}"
