"""`MAT-15` step 3 (§9 item 236, chain P5) -- heating per accepted watt on the
tuned copper F-small birdcage at 64.0 MHz, degree 2, exact rectangles, κ off.

Fixture: the ×1 exact-strip hole mesh (``_exact_geo(th14, comm,
GATED_WIDTH_FRACTION)``, 80 798 cells), copper by Leontovich, Z_C(C_mode,Cu) on
P2–P4 and C₁ across P1 (circuit layer, lossless), one drive of P1 at 64.0 MHz
through ``_route_a_solve(..., degree=2)``.  C_mode,Cu is read from T21's
gitignored ``logs/ans7_exact/tune_p2.json`` and never recomputed.  **This
module does not claim the coil is tuned to 64.000 MHz** (known-issues
2026-10-04, T21): every asserted anchor is an identity that holds for any C.

Normalisation G = 1 W / P_acc, P_acc = P_src − ΣP_sheet (one accepted watt).
Thermal side: step 1's (`test_mat15_coil_driven_heating_phantom.py`) loop
body for one drive, on the tag-3 submesh.

Anchors (bands imported, none restated): (F) fixture tie; (Φ) |sin φ₁| ≤ 1e-9;
(W) S11 against T21's 64.0 MHz S11_A ≤ ``NODE_S_RTOL``; (P) accepted-power
identity ≤ ``DISCRETE_IDENTITY_RTOL``; (B) ≤ ``B_BAND``; (E), (E′) ≤ ``E_BAND``;
(S2) ≤ ``S2_BAND``; (M1)/(M2) σ_y mirror 10 g ball means ≤ ``C4_COVARIANCE_BAND``.
Negative controls (rule (e)): (d) doubled source, asserted only > ``D_CONTROL_FLOOR``;
(f′) off-mirror pair θ̄(φ₁+90°) vs θ̄(φ₁), asserted only > ``F_CONTROL_FLOOR``;
(P-ctl) printed.

Opt-in knob ``FEM_EM_MAT15_STEP3=1``; unset ⇒ every test skips; any other
value raises.  Writes nothing outside the harness log.
"""

from __future__ import annotations

import dataclasses
import json
import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.core.cavity import surface_resistance_ohm
from fem_em_solver.post.sar import averaging_ball_radius, mean_sar
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import _sheet_azimuth_deg
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans7_exact_strip import (
    GATED_WIDTH_FRACTION,
    OUT_DIR,
    _assert_tie,
    _exact_geo,
    _fixture_tie,
    _route_a_solve,
)
from tests.validation.test_ans7_sparse_route_b import _ANS7, NODE_S_RTOL, _p_inc_w
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    PHANTOM_RHO_KG_PER_M3,
)
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_birdcage_sar_mass_averaged import CENTRE_RADIUS_M
from tests.validation.test_mass_averaged_sar_standard_masses import (
    QUADRATURE_DEGREE,
    TEN_GRAM_KG,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    B_BAND,
    ROBIN_TAG,
    _allsum,
    _robin_tags,
    _say,
    _to_sub,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    E_BAND,
    RHO_C,
    S2_BAND,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    CAD_AREA_M2,
    D_CONTROL_FLOOR,
    DT_1,
    F_CONTROL_FLOOR,
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
    T_R,
    T_STAR,
    A,
    _be_one,
    _cell_average_dg0_tag,
    _cn,
    _e_prime,
    _max_imag,
    _tag_submesh,
)

KNOB = "FEM_EM_MAT15_STEP3"
DRIVE_HZ = 64.0e6
PHI_BAND = 1.0e-9  # (Φ) precondition: |sin φ₁|
V_CYL_M3 = 2.26195e-4  # πa²·2L, printed only
V_SUB_STEP1_RECORD_M3 = 2.259978941112e-04  # ANS-2 rung, printed only
P1_THETA_PER_DEP_W = 1.153  # step 1 mean θ(t*) per deposited watt, printed only
KAPPA_S_STEP1 = 1.65  # printed only


def _enabled():
    raw = os.environ.get(KNOB, "").strip()
    if not raw:
        return False
    if raw != "1":
        raise ValueError(f"{KNOB} must be unset or '1', got {raw!r}")
    return True


def _need():
    if not _enabled():
        pytest.skip(f"{KNOB} unset (heavy -n 8 window, gitignored T21 inputs)")


def _read_t21():
    p_tune = OUT_DIR / "tune_p2.json"
    p_in = OUT_DIR / "tune_inmodel_p2.json"
    assert p_tune.exists() and p_in.exists(), f"T21 inputs absent: {p_tune}, {p_in}"
    tune = json.loads(p_tune.read_text())
    inm = json.loads(p_in.read_text())
    assert tune["degree"] == 2 and inm["degree"] == 2, (tune["degree"], inm["degree"])
    assert abs(float(tune["f_hz"]) - 64.0e6) < 1.0, tune["f_hz"]
    assert int(tune["cells"]) == 80798 and int(inm["cells"]) == 80798
    assert list(inm["tune_mhz"]) == [63.9, 64.0, 64.1], inm["tune_mhz"]
    c_mode = float(tune["c_mode"]["cu"])
    assert float(inm["c_mode"]) == c_mode, (inm["c_mode"], c_mode)
    s11_ref = complex(*inm["s11_a"][1])
    return c_mode, s11_ref, inm


def _run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    c_mode, s11_ref, inm = _read_t21()
    _say(
        f"\nMAT15s3 C_mode,Cu = {c_mode:.15e} F (read from tune_p2.json, not recomputed); "
        f"S11_ref(64.0 MHz) = {s11_ref.real:+.15e}{s11_ref.imag:+.15e}j; in-model zero "
        f"{inm['rows_a']['f0_hz']:.3f} Hz (NOT a tuning claim; drive at {DRIVE_HZ:.1f} Hz)"
    )
    tie = _fixture_tie()
    geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
    msh, cell_tags, facet_tags, _outer, tags_f, _pdefs = geo
    t_mesh = time.perf_counter() - t0
    phi1 = _sheet_azimuth_deg(msh, tags_f, SHEET_IFACE + 1, comm)
    sin_phi = abs(np.sin(np.radians(phi1)))
    _say(
        f"MAT15s3 (F) cells {cells}; (Phi) phi_1 = {phi1:.12f} deg, |sin phi_1| = {sin_phi:.3e} "
        f"(bound {PHI_BAND:.0e}); mesh+tie {t_mesh:.1f} s"
    )
    out = {"cells": cells, "tie": tie, "phi1": phi1, "sin_phi": sin_phi}

    t1 = time.perf_counter()
    s11, fields = _route_a_solve(
        th14, geo, specs, DRIVE_HZ, "cu", c_f=c_mode, return_fields=True, degree=2
    )
    comm.Barrier()
    t_em = time.perf_counter() - t1
    out["W"] = abs(s11 - s11_ref) / abs(s11_ref)
    _say(
        f"MAT15s3 (W) S11 = {s11.real:+.15e}{s11.imag:+.15e}j  rel to T21 {out['W']:.3e} "
        f"(bound {NODE_S_RTOL:.0e}, predicted <= 1e-12); degree-2 drive {t_em:.1f} s"
    )

    # (P) -- TH-14's (b) block on the Z_C-replaced specs.
    omega = 2.0 * np.pi * DRIVE_HZ
    z_c = _ANS7.z_cap(omega, c_mode)
    sp_list = [
        (
            dataclasses.replace(s, port_impedance_ohm=z_c)
            if i in _ANS7.TUNING_TERMINATED_INDICES
            else s
        )
        for i, s in enumerate(specs)
    ]
    e = fields.e_complex
    sigma_f = fields.sigma_field
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)

    def _vol(measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    t2 = time.perf_counter()
    acct = {"mesh": msh, "facet_tags": tags_f}
    sheet_objs = {sp.port_id: sp.sheet(driven=(sp.port_id == th14.DRIVEN)) for sp in sp_list}
    p_sheet = {
        k: float(_sheet_field_dissipation_w(acct, sh, e, omega)) for k, sh in sheet_objs.items()
    }
    p_sheets = sum(p_sheet.values())
    p_src = float(_source_power_w(acct, sheet_objs[th14.DRIVEN], e, omega))
    p_total = _vol(dx)
    p_phantom = _vol(dx(th14.PHANTOM_CELL_TAG))
    p_ph_sar = float(
        np.real(
            mean_sar(
                e,
                sigma=sigma_f,
                rho=PHANTOM_RHO_KG_PER_M3,
                cell_tags=cell_tags,
                comm=comm,
                subdomain_ids=th14.PHANTOM_CELL_TAG,
            )["dissipated_power_w"]
        )
    )
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, th14.COPPER_SIGMA)
    p_surf = th14._surface_loss_w(msh, facet_tags, z_s, e, comm)
    p_non = p_total - p_phantom
    out["P"] = abs(p_src - p_sheets - p_phantom - p_non - p_surf) / abs(p_src)
    p_ctl = abs(p_src - p_phantom - p_non - p_surf) / abs(p_src)
    p_acc = p_src - p_sheets
    G = 1.0 / p_acc
    eta = p_phantom / p_acc
    l_s = 1.0 - abs(s11) ** 2
    l_f = p_acc / _p_inc_w()
    out.update(G=G, eta=eta, p_acc=p_acc)
    _say(
        f"MAT15s3 (P) P_src = {p_src:.15e} W; sheets "
        + "  ".join(f"{k}:{v:.6e}" for k, v in p_sheet.items())
        + f" (sum {p_sheets:.15e}); P_total {p_total:.15e}; P_phantom {p_phantom:.15e}; "
        f"P_non-phantom {p_non:.15e}; P_surf {p_surf:.15e}"
    )
    _say(
        f"MAT15s3 (P) identity rel {out['P']:.3e} (bound {DISCRETE_IDENTITY_RTOL:.0e}, predicted "
        f"<= 1e-10); (P-ctl) sheets dropped {p_ctl:.6f} (predicted 0.998, printed)"
    )
    _say(
        f"MAT15s3 P_phantom mean_sar vs _vol rel {abs(p_ph_sar / p_phantom - 1.0):.3e} "
        f"(predicted <= 1e-12); P_acc = {p_acc:.15e} W, G = 1/P_acc = {G:.12e} 1/W"
    )
    _say(
        f"MAT15s3 eta = P_phantom/P_acc = {eta:.9f} (predicted [0.52, 0.62], estimate 0.582); "
        f"P_surf/P_acc = {p_surf / p_acc:.9f} (1 - eta = {1.0 - eta:.9f}); P_non/P_acc = "
        f"{p_non / p_acc:.3e}"
    )
    _say(
        f"MAT15s3 L_S = 1-|S11|^2 = {l_s:.9e}  L_F = P_acc/P_inc = {l_f:.9e}  |diff| = "
        f"{abs(l_s - l_f):.3e} (predicted <= 1e-4, printed)"
    )
    _say(
        f"MAT15s3 detuning estimate: (308.55 Hz / 32.9 kHz)^2 = {(308.55 / 32.9e3) ** 2:.2e} "
        f"(printed)"
    )

    # Thermal -- step 1's per-drive loop body, one drive.
    sub, parent_ids = _tag_submesh(msh, cell_tags, th14.PHANTOM_CELL_TAG)
    ft = _robin_tags(sub)
    n_sub = sub.topology.index_map(3).size_global
    one = fem.Constant(sub, default_scalar_type(1.0))
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    area = float(np.real(_allsum(comm, one * ds)))
    _say(
        f"MAT15s3 sub cells = {n_sub} (predicted [300, 3000]); V_sub = {v_sub:.12e} m^3 vs "
        f"pi a^2 2L {V_CYL_M3:.5e} ({v_sub / V_CYL_M3 - 1.0:+.4e}, predicted [-3e-2, 0]); step 1 "
        f"record {V_SUB_STEP1_RECORD_M3:.12e} (printed); oint ds = {area:.9e} m^2 vs CAD "
        f"{CAD_AREA_M2} ({area / CAD_AREA_M2 - 1.0:+.4e}, printed)"
    )
    _say(
        f"MAT15s3 t_R = {T_R:.3f} s  t* = {T_STAR:.3f} s  dt = {T_STAR / N_MAIN:.4f} s  "
        f"dt_1 = {DT_1:.6e} s  Bi = {H_GEL * A / K_GEL:.3f}"
    )
    out["n_sub"] = n_sub

    Wp, q_par = _cell_average_dg0_tag(
        msh, cell_tags, th14.PHANTOM_CELL_TAG, G * 0.5 * sigma_f * ufl.inner(e, e)
    )
    q_h = _to_sub(Wp, q_par, sub, parent_ids)
    dep = G * p_phantom
    int_q = float(np.real(_allsum(comm, q_h * ufl.dx(domain=sub))))
    out["B"] = abs(int_q - dep) / abs(dep)
    q_s = float(np.real(_allsum(comm, q_h * ds)))
    kappa = (q_s / area) / (int_q / v_sub)

    t3 = time.perf_counter()
    res = _cn(sub, ft, q_h)
    t_cn = time.perf_counter() - t3
    ep, stor, conv = _e_prime(sub, ds, res, dep)
    out["Ep"] = ep

    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    centre_c = fem.Constant(sub, np.zeros(3, dtype=default_scalar_type))
    xs = ufl.SpatialCoordinate(sub)
    off = xs - centre_c
    inside = ufl.conditional(ufl.lt(ufl.real(ufl.dot(off, off)), a10**2), 1.0, 0.0)
    dxq = ufl.dx(domain=sub, metadata={"quadrature_degree": int(QUADRATURE_DEGREE)})
    V1 = fem.functionspace(sub, ("Lagrange", 1))
    Ws = fem.functionspace(sub, ("DG", 0))
    th_ball = fem.Function(V1)
    q_ball = fem.Function(Ws)
    f_vol = fem.form(inside * dxq)
    f_th = fem.form(th_ball * inside * dxq)
    f_q = fem.form(q_ball * inside * dxq)

    def _red(form):
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    centres = {}
    for m in range(8):
        ang = np.radians(phi1 + 45.0 * m)
        centres[m] = np.array([CENTRE_RADIUS_M * np.cos(ang), CENTRE_RADIUS_M * np.sin(ang), 0.0])
    ball_vol = {}
    v_ball = 4.0 / 3.0 * np.pi * a10**3
    for m in range(8):
        centre_c.value[:] = centres[m]
        ball_vol[m] = _red(f_vol)
    _say(
        "MAT15s3 ball volumes [m^3] (4/3 pi a10^3 = %.9e): " % v_ball
        + "  ".join(f"m{m}:{ball_vol[m]:.9e}" for m in range(8))
    )

    def _ball_mean(form, m):
        centre_c.value[:] = centres[m]
        return _red(form) / ball_vol[m]

    th_ball.x.array[:] = res.temperature.x.array - T_INF
    th_ball.x.scatter_forward()
    mean_th = float(np.real(_allsum(comm, th_ball * ufl.dx(domain=sub)))) / v_sub
    q_ball.x.array[:] = q_h.x.array
    q_ball.x.scatter_forward()
    th = {m: _ball_mean(f_th, m) for m in range(8)}
    qb = {m: _ball_mean(f_q, m) for m in range(8)}

    res1 = _be_one(sub, ft, q_h)
    thbar1 = float(np.real(_allsum(comm, (res1.temperature - T_INF) * ufl.dx(domain=sub)))) / v_sub
    ratio = RHO_C * thbar1 / (DT_1 * dep / v_sub)
    leak = float(np.real(_allsum(comm, H_GEL * (res1.temperature - T_INF) * ds))) / dep
    estimate = H_GEL * DT_1 * area / (RHO_C * v_sub)
    out["S2"] = abs(ratio - 1.0)
    out["E"] = {"cn": res.max_balance_residual, "be": res1.max_balance_residual}
    tmax_imag = max(_max_imag(res.temperature), _max_imag(res1.temperature))

    _say(
        f"MAT15s3 (B) int_sub q_h = {int_q:.15e}  G*P_phantom = {dep:.15e}  rel {out['B']:.3e} "
        f"(bound {B_BAND:.0e}, predicted <= 1e-13);  kappa_s = {kappa:.4f} (predicted [1.5, 6]; "
        f"step 1 {KAPPA_S_STEP1}, printed)"
    )
    _say(
        f"MAT15s3 CN N={N_MAIN}: mean theta(t*) per accepted watt = {mean_th:.6f} K (predicted "
        f"[0.55, 0.75]; eta x step 1's {P1_THETA_PER_DEP_W} = {eta * P1_THETA_PER_DEP_W:.4f}, "
        f"printed)  (E) max_balance = {res.max_balance_residual:.3e}  solve {t_cn:.1f} s"
    )
    _say(
        f"MAT15s3 (E') stored {stor:.12e} W + surface {conv:.12e} W vs deposited G*P_phantom "
        f"{dep:.12e} W: rel {ep:.3e} (bound {E_BAND:.0e}, predicted <= 1e-12)"
    )
    _say(
        f"MAT15s3 (S2) one BE step: ratio = {ratio:.12f}  |ratio-1| = {out['S2']:.4e} (bound "
        f"{S2_BAND:.0e}, predicted [5e-6, 2.1e-5]); deficit {leak:.4e}  estimate "
        f"{estimate:.3e}*kappa_s = {estimate * kappa:.4e}  balance {res1.max_balance_residual:.3e}"
    )
    for m in range(8):
        _say(
            f"MAT15s3   10 g ball m={m} (phi_1{45 * m:+d} deg): theta = {th[m]:.9e} K  "
            f"q_h = {qb[m]:.9e} W/m^3"
        )
    peak = max(th, key=th.get)
    _say(
        f"MAT15s3 peak 10 g theta per accepted watt = {th[peak]:.9e} K at azimuth "
        f"{(phi1 + 45.0 * peak) % 360.0:.3f} deg (m={peak}, printed)"
    )

    out["M1"] = abs(th[2] - th[6]) / abs(th[2])
    out["M2"] = abs(th[1] - th[7]) / abs(th[1])
    q_m1 = abs(qb[2] - qb[6]) / abs(qb[2])
    q_m2 = abs(qb[1] - qb[7]) / abs(qb[1])
    c2 = abs(th[0] - th[4]) / abs(th[0])
    out["fprime"] = abs(th[2] - th[0]) / abs(th[0])
    _say(
        f"MAT15s3 (M1) theta(phi_1+90) vs theta(phi_1-90): {out['M1'] * 100:.4f}% (bound "
        f"{C4_COVARIANCE_BAND * 100:.1f}%, predicted <= 1%); q_h pair {q_m1 * 100:.4f}% "
        f"(printed, predicted <= 3%)"
    )
    _say(
        f"MAT15s3 (M2) theta(phi_1+45) vs theta(phi_1-45): {out['M2'] * 100:.4f}% (bound "
        f"{C4_COVARIANCE_BAND * 100:.1f}%, predicted <= 1%); q_h pair {q_m2 * 100:.4f}% "
        f"(printed, predicted <= 3%)"
    )
    _say(
        f"MAT15s3 C2 pair theta(phi_1) vs theta(phi_1+180): {c2 * 100:.4f}% (not a symmetry; "
        f"printed, predicted <= 2%)"
    )
    _say(
        f"MAT15s3 control (f') theta(phi_1+90) vs theta(phi_1): {out['fprime'] * 100:.4f}% "
        f"(asserted > {F_CONTROL_FLOOR * 100:.1f}%, predicted 40-80%)"
    )

    q_two = fem.Function(q_h.function_space)
    q_two.x.array[:] = 2.0 * q_h.x.array
    q_two.x.scatter_forward()
    res_d = _cn(sub, ft, q_two)
    ep_d, _sd, _cd = _e_prime(sub, ds, res_d, dep)
    out["ctl_d"] = ep_d
    out["E"]["ctl_d"] = res_d.max_balance_residual
    out["imag"] = tmax_imag
    _say(
        f"MAT15s3 control (d) doubled source: (E') vs G*P_phantom = {ep_d:.6f} (predicted "
        f"1.000, asserted > {D_CONTROL_FLOOR:.0e}); its (E) = {res_d.max_balance_residual:.3e}"
    )
    worst = max(out["E"], key=out["E"].get)
    _say(
        f"MAT15s3 (E) worst balance over {len(out['E'])} runs = {out['E'][worst]:.3e} on "
        f"{worst} (bound {E_BAND:.0e}); max|Im theta| = {tmax_imag:.3e} (predicted 0)"
    )
    _say(
        f"MAT15s3 elapsed: mesh+tie {t_mesh:.1f} s, EM {t_em:.1f} s, post+thermal "
        f"{time.perf_counter() - t2:.1f} s, total {time.perf_counter() - t0:.1f} s"
    )
    return out


@pytest.fixture(scope="module")
def tuned():
    _need()
    return _run()


@complex_only
def test_f_fixture_tie(tuned):
    _assert_tie(tuned["cells"], tuned["tie"])


@complex_only
def test_phi_driven_leg_on_mirror_plane(tuned):
    assert tuned["sin_phi"] <= PHI_BAND, (tuned["phi1"], tuned["sin_phi"])


@complex_only
def test_w_source_wiring_reproduces_t21(tuned):
    assert tuned["W"] <= NODE_S_RTOL, f"(W) {tuned['W']:.3e} > {NODE_S_RTOL:.0e}"


@complex_only
def test_p_accepted_power_identity(tuned):
    assert tuned["P"] <= DISCRETE_IDENTITY_RTOL, f"(P) {tuned['P']:.3e}"


@complex_only
def test_b_source_transfer(tuned):
    assert tuned["B"] <= B_BAND, f"(B) {tuned['B']:.3e} > {B_BAND:.0e}"


@complex_only
def test_e_every_run_balances(tuned):
    for name, r in tuned["E"].items():
        assert r <= E_BAND, f"(E) run {name}: {r:.3e} > {E_BAND:.0e}"


@complex_only
def test_e_prime_energy_identity_from_em_side(tuned):
    assert tuned["Ep"] <= E_BAND, f"(E') {tuned['Ep']:.3e} > {E_BAND:.0e}"


@complex_only
def test_s2_calorimetric_slope(tuned):
    assert tuned["S2"] <= S2_BAND, f"(S2) {tuned['S2']:.4e} > {S2_BAND:.0e}"


@complex_only
@pytest.mark.parametrize("pair", ["M1", "M2"])
def test_m_sigma_y_mirror_ten_gram_theta(tuned, pair):
    r = tuned[pair]
    assert r <= C4_COVARIANCE_BAND, f"({pair}) {r * 100:.4f}% > {C4_COVARIANCE_BAND * 100:.1f}%"


@complex_only
def test_control_d_doubled_source_misses(tuned):
    assert tuned["ctl_d"] > D_CONTROL_FLOOR, f"control (d) {tuned['ctl_d']:.3e} did not miss"


@complex_only
def test_control_f_prime_off_mirror_pair_misses(tuned):
    r = tuned["fprime"]
    assert r > F_CONTROL_FLOOR, f"control (f') {r * 100:.4f}% inside {F_CONTROL_FLOOR * 100:.1f}%"
