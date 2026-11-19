"""`MAT-15` step 1 — coil-driven SAR heats the gel phantom (§9 item 130, chain P1).

Fixture: `ANS-2`'s 10 MHz loaded F-small birdcage on the 0.0025 m phantom rung
(``build_four_port_sweep(phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG)``,
199 920 / 58 866 cells), the four single-leg drives through ``_solve_driven``,
C4 ordering as ``test_birdcage_sar_mass_averaged.py:239–247``.  Per drive
P_k = ``mean_sar(…, subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"]`` and
the heat source is q_h = G·(DG0 cell average of ½σ|E|² on tag 3), transferred
to the tag-3 submesh.  **One** G = 1 W / P₁ for all four drives — a per-drive G
would hide a C4 power asymmetry.

Thermal side: the tag-3 submesh, Robin on every exterior facet, gel k 0.6,
ρc 4.18e6, h 10 (Bi = 0.5 at a = ``PHANTOM_RADIUS``), T∞ = T(0) = 20, every
error on θ = T − 20; t_R = a²ρc/k = 6270 s, t* = 0.2 t_R, CN N = 40 with the
Rannacher start (`MAT-12`'s ``_cn`` shape); Δt₁ = 2.5e-6 t_R.

Anchors (asserted): (B) ∫_sub q_h dV = G·P_k ≤ ``B_BAND`` and V_sub against the
POST-8 record ≤ 1e-10; (E) every run's ``max_balance_residual`` ≤ ``E_BAND``;
(E′) per drive, from snapshots N−1 / N, ρc∫Δθ/Δt + ∮hθ^(n+½) = G·P_k **from the
EM side** ≤ 1e-9; (S2) one BE step Δt₁ from θ = 0, the volume-mean slope
|ρc·θ̄¹/(Δt₁·G·P_k/V_sub) − 1| ≤ ``S2_BAND``; (C4) 10 g ball means of θ at t* at
the four rotated centres, four cyclic pairs ≤ ``C4_COVARIANCE_BAND``.
Negative controls (rule (e)): (d) doubled source — (E′) against G·P₁ asserted
only > 1e-3 (predicted 1.000); (f) far-side ball θ̄(c_{k+2}; k) against
θ̄(c_k; k) asserted > 5 % (predicted 40–80 %).

Scope: one rung, degree-1 EM, P1 thermal, no perfusion, constant σ, a 1 W
normalisation.  Not: a closed form for the coil-driven θ, convergence, an
absolute temperature, `ANS-9`, tissue, ASTM F2182.
"""

from __future__ import annotations

import time

import dolfinx
import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.post.sar import averaging_ball_radius, mean_sar
from fem_em_solver.thermal import solve_transient_conduction
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
from tests.mesh.test_birdcage_port_tags import PHANTOM_RADIUS
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    PHANTOM_RHO_KG_PER_M3,
    _solve_driven,
)
from tests.validation.test_birdcage_b1_quadrature import (
    QUADRATURE_STEP_DEG,
    _port_index,
)
from tests.validation.test_birdcage_sar_1g_rung import PHANTOM_RESOLUTION_1G_RUNG
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
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

# Gel constants (pinned by the 2026-10-01 03:00 review, §7 MAT-15 row).
K_GEL = 0.6  # W/m/K
H_GEL = 10.0  # W/m²/K; Bi = hR/k = 0.5 at R = PHANTOM_RADIUS
T_INF = 20.0  # = T(0); errors are on θ = T − 20
A = PHANTOM_RADIUS
T_R = A**2 * RHO_C / K_GEL  # 6270 s
T_STAR = 0.2 * T_R
N_MAIN = 40
DT_1 = 2.5e-6 * T_R  # the calorimetric step

# Records on this rung (`20260927T112115Z_POST-8-step1-n2.log:1878`).  V_sub is
# the line's mass 2.259978941112e-01 kg / ρ = 1000 (13 digits) — the volume
# column printed beside it (`:1879`, 2.259978941e-04) carries only 10 digits, a
# rounding of 2e-10 relative, which is above the 1e-10 band itself; the two are
# the same quantity (`:1880`: mass vs ρ·V_phantom 1.98e-13).
V_SUB_RECORD_M3 = 2.259978941112e-04
P1_RECORD_W = 5.503204728588e-08
V_BAND = 1.0e-10
CAD_AREA_M2 = 0.020735  # 2πa·0.08 + 2πa², printed only

# Negative controls (rule (e)): asserted for sign only.
D_CONTROL_FLOOR = 1.0e-3  # (d) predicted 1.000
F_CONTROL_FLOOR = C4_COVARIANCE_BAND  # (f) predicted 40–80 %


def _tag_submesh(msh, cell_tags, tag):
    """Cell tag ``tag`` of the parent as its own mesh (MAT-10's, tag-generic)."""
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_parent = imap.size_local + imap.num_ghosts
    parent_values = np.full(n_parent, -1, dtype=np.int32)
    parent_values[cell_tags.indices] = cell_tags.values
    keep = np.flatnonzero(parent_values == tag).astype(np.int32)
    sub, emap, _vmap, _gmap = dmesh.create_submesh(msh, tdim, keep)
    simap = sub.topology.index_map(tdim)
    n_sub = simap.size_local + simap.num_ghosts
    parent_ids = np.asarray(
        emap.sub_topology_to_topology(np.arange(n_sub, dtype=np.int32), False),
        dtype=np.int32,
    )
    return sub, parent_ids


def _cell_average_dg0_tag(msh, cell_tags, tag, integrand):
    """Parent DG0 cell average of ``integrand`` over ``tag`` (MAT-10's, tag-generic)."""
    W = fem.functionspace(msh, ("DG", 0))
    w = ufl.TestFunction(W)
    dxt = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)(tag)
    one = fem.Constant(msh, default_scalar_type(1.0))

    def _vec(form):
        b = fem.assemble_vector(fem.form(form))
        b.scatter_reverse(dolfinx.la.InsertMode.add)
        b.scatter_forward()
        return b.array.copy()

    num = _vec(ufl.inner(integrand, w) * dxt)
    vol = _vec(ufl.inner(one, w) * dxt)
    avg = np.zeros_like(num)
    mask = np.abs(vol) > 0.0
    avg[mask] = num[mask] / vol[mask]
    return W, avg


def _transient(sub, ft, q, *, dt, n_steps, theta, startup_half_steps=0, snapshot_steps=()):
    return solve_transient_conduction(
        sub,
        ft,
        k=K_GEL,
        rho_c=RHO_C,
        q=q,
        h=H_GEL,
        t_inf=T_INF,
        t_init=T_INF,
        robin_tag=ROBIN_TAG,
        dt=dt,
        n_steps=n_steps,
        theta=theta,
        startup_half_steps=startup_half_steps,
        degree=1,
        snapshot_steps=snapshot_steps,
    )


def _cn(sub, ft, q):
    """CN, N = 40 to t*, two BE half-steps first; snapshots N−1, N for (E′)."""
    return _transient(
        sub,
        ft,
        q,
        dt=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        snapshot_steps=(N_MAIN - 1, N_MAIN),
    )


def _be_one(sub, ft, q):
    return _transient(sub, ft, q, dt=DT_1, n_steps=1, theta=1.0)


def _e_prime(sub, ds, res, deposited):
    """ρc∫(θᴺ−θᴺ⁻¹)/Δt + ∮hθ^(n+½) against the EM-side deposited power."""
    comm = sub.comm
    tn, tm = res.snapshots[N_MAIN], res.snapshots[N_MAIN - 1]
    dt = T_STAR / N_MAIN
    stor = float(np.real(_allsum(comm, RHO_C * (tn - tm) / dt * ufl.dx(domain=sub))))
    conv = float(np.real(_allsum(comm, H_GEL * (0.5 * (tn + tm) - T_INF) * ds)))
    return abs(stor + conv - deposited) / abs(deposited), stor, conv


def _max_imag(T):
    n_own = T.function_space.dofmap.index_map.size_local
    loc = float(np.max(np.abs(np.imag(T.x.array[:n_own])))) if n_own else 0.0
    return T.function_space.mesh.comm.allreduce(loc, op=MPI.MAX)


def _run_heated(*, return_fields=False):
    """Step 1's fixture body, lifted to module level (EX-79, rule (a)).

    ``return_fields=True`` additionally returns the mesh objects and per-drive
    fields (``msh``, ``cell_tags``, ``sub``, ``ft``, ``by_k``, ``centres``,
    ``q_h`` and ``T_star`` keyed by k, plus ``sweep`` and ``solves`` — ANS-9
    step 1a, additive); prints and the other keys are unchanged.
    """
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()

    sweep = build_four_port_sweep(phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG)
    msh, cell_tags = sweep["mesh"], sweep["cell_tags"]
    t_mesh = time.perf_counter() - t0

    azimuths = {f"P{s['tag'] - SHEET_IFACE}": float(s["azimuth_deg"]) for s in sweep["sheets"]}
    order = sorted(azimuths)
    indices = {
        pid: _port_index(azimuths[pid], azimuths["P1"], QUADRATURE_STEP_DEG) for pid in order
    }
    by_k = {indices[pid]: pid for pid in order}
    solves = {pid: _solve_driven(sweep, pid) for pid in order}
    t_em = time.perf_counter() - t0

    P = {}
    for k in range(4):
        f = solves[by_k[k]]["fields"]
        P[k] = float(
            mean_sar(
                f.e_complex,
                sigma=f.sigma_field,
                rho=PHANTOM_RHO_KG_PER_M3,
                cell_tags=cell_tags,
                comm=comm,
                subdomain_ids=PHANTOM_CELL_TAG,
            )["dissipated_power_w"]
        )
    G = 1.0 / P[0]
    _say(
        f"\nMAT15 mesh {int(sweep['cells'])} cells, f = {sweep['problem'].frequency_hz:.3e} Hz; "
        f"mesh+rung {t_mesh:.1f} s, EM total {t_em:.1f} s; drive order "
        + ", ".join(f"k={k}:{by_k[k]}@{azimuths[by_k[k]]:.1f}deg" for k in range(4))
    )
    _say(
        f"MAT15 G = 1 W / P_1 = {G:.12e} 1/W;  P_k [W] = "
        + "  ".join(f"{P[k]:.12e}" for k in range(4))
        + f"  (P_1 record {P1_RECORD_W:.12e}, rel {abs(P[0] / P1_RECORD_W - 1.0):.3e}, printed)"
    )

    sub, parent_ids = _tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)
    ft = _robin_tags(sub)
    n_sub = sub.topology.index_map(3).size_global
    one = fem.Constant(sub, default_scalar_type(1.0))
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    area = float(np.real(_allsum(comm, one * ds)))
    err_V = abs(v_sub / V_SUB_RECORD_M3 - 1.0)
    _say(
        f"MAT15 sub cells = {n_sub} (record 58866); V_sub = {v_sub:.12e} m^3 vs record "
        f"{V_SUB_RECORD_M3:.12e} rel {err_V:.3e} (bound {V_BAND:.0e}); oint ds = {area:.9e} m^2 "
        f"vs CAD {CAD_AREA_M2} ({area / CAD_AREA_M2 - 1.0:+.4e}, printed)"
    )
    _say(
        f"MAT15 t_R = {T_R:.3f} s  t* = {T_STAR:.3f} s  dt = {T_STAR / N_MAIN:.4f} s  "
        f"dt_1 = {DT_1:.6e} s  Bi = {H_GEL * A / K_GEL:.3f}"
    )

    # Ball-mean forms, compiled once each (centre as a Constant).
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

    centres = {
        k: np.array(
            [
                CENTRE_RADIUS_M * np.cos(np.radians(azimuths[by_k[k]])),
                CENTRE_RADIUS_M * np.sin(np.radians(azimuths[by_k[k]])),
                0.0,
            ]
        )
        for k in range(4)
    }
    ball_vol = {}
    for j in range(4):
        centre_c.value[:] = centres[j]
        ball_vol[j] = _red(f_vol)

    def _ball_mean(form, j):
        centre_c.value[:] = centres[j]
        return _red(form) / ball_vol[j]

    out = {
        "P": P,
        "G": G,
        "v_sub": v_sub,
        "err_V": err_V,
        "B": {},
        "E": {},
        "Ep": {},
        "S2": {},
        "th_own": {},
        "th_far": {},
        "q_own": {},
    }
    tmax_imag = 0.0
    q_first = None
    q_h_by_k, t_star_by_k = {}, {}
    for k in range(4):
        f = solves[by_k[k]]["fields"]
        Wp, q_par = _cell_average_dg0_tag(
            msh,
            cell_tags,
            PHANTOM_CELL_TAG,
            G * 0.5 * f.sigma_field * ufl.inner(f.e_complex, f.e_complex),
        )
        q_h = _to_sub(Wp, q_par, sub, parent_ids)
        dep = G * P[k]
        int_q = float(np.real(_allsum(comm, q_h * ufl.dx(domain=sub))))
        out["B"][k] = abs(int_q - dep) / abs(dep)
        q_s = float(np.real(_allsum(comm, q_h * ds)))
        kappa = (q_s / area) / (int_q / v_sub)

        t1 = time.perf_counter()
        res = _cn(sub, ft, q_h)
        t_cn = time.perf_counter() - t1
        ep, stor, conv = _e_prime(sub, ds, res, dep)
        out["Ep"][k] = ep
        th_ball.x.array[:] = res.temperature.x.array - T_INF
        th_ball.x.scatter_forward()
        mean_th = float(np.real(_allsum(comm, th_ball * ufl.dx(domain=sub)))) / v_sub
        out["th_own"][k] = _ball_mean(f_th, k)
        out["th_far"][k] = _ball_mean(f_th, (k + 2) % 4)
        q_ball.x.array[:] = q_h.x.array
        q_ball.x.scatter_forward()
        out["q_own"][k] = _ball_mean(f_q, k)

        res1 = _be_one(sub, ft, q_h)
        thbar1 = (
            float(np.real(_allsum(comm, (res1.temperature - T_INF) * ufl.dx(domain=sub)))) / v_sub
        )
        ratio = RHO_C * thbar1 / (DT_1 * dep / v_sub)
        leak = float(np.real(_allsum(comm, H_GEL * (res1.temperature - T_INF) * ds))) / dep
        estimate = H_GEL * DT_1 * area / (RHO_C * v_sub)
        out["S2"][k] = abs(ratio - 1.0)
        out["E"][f"cn{k}"] = res.max_balance_residual
        out["E"][f"be{k}"] = res1.max_balance_residual
        tmax_imag = max(tmax_imag, _max_imag(res.temperature), _max_imag(res1.temperature))

        _say(
            f"MAT15 drive k={k} ({by_k[k]}): (B) int_sub q_h = {int_q:.15e}  G*P_k = {dep:.15e}  "
            f"rel {out['B'][k]:.3e} (bound {B_BAND:.0e});  kappa_s = {kappa:.4f} (predicted 2-5)"
        )
        _say(
            f"MAT15   CN N={N_MAIN}: mean theta(t*) = {mean_th:.6f} K (predicted ≈ 1.2)  "
            f"(E) max_balance = {res.max_balance_residual:.3e}  steps {len(res.balance_residuals)}  "
            f"solve {t_cn:.1f} s"
        )
        _say(
            f"MAT15   (E') stored {stor:.12e} W + surface {conv:.12e} W vs deposited G*P_k "
            f"{dep:.12e} W: rel {ep:.3e} (bound {E_BAND:.0e})"
        )
        _say(
            f"MAT15   (S2) one BE step: rho c mean(theta^1)/(dt_1 G P_k/V) = {ratio:.12f}  "
            f"|ratio-1| = {abs(ratio - 1.0):.4e} (bound {S2_BAND:.0e}); deficit oint h theta^1/(G P_k) "
            f"= {leak:.4e}  estimate {estimate:.3e}*kappa_s = {estimate * kappa:.4e}  "
            f"balance {res1.max_balance_residual:.3e}"
        )
        _say(
            f"MAT15   10 g ball means at t*: theta(c_k; k) = {out['th_own'][k]:.9e} K  "
            f"theta(c_k+2; k) = {out['th_far'][k]:.9e} K;  q_h(c_k; k) = {out['q_own'][k]:.9e} W/m^3"
        )
        if k == 0:
            q_first = q_h
        q_h_by_k[k] = q_h
        t_star_by_k[k] = res.temperature

    # Negative control (d): doubled source on drive k = 0.
    q_two = fem.Function(q_first.function_space)
    q_two.x.array[:] = 2.0 * q_first.x.array
    q_two.x.scatter_forward()
    res_d = _cn(sub, ft, q_two)
    ep_d, stor_d, conv_d = _e_prime(sub, ds, res_d, G * P[0])
    out["ctl_d"] = ep_d
    out["ctl_d_E"] = res_d.max_balance_residual
    _say(
        f"MAT15 control (d) doubled source, k=0: (E') vs G*P_1 = {ep_d:.6f} (predicted 1.000, "
        f"asserted > {D_CONTROL_FLOOR:.0e});  its internal (E) = {res_d.max_balance_residual:.3e} "
        f"(predicted ≤ 1e-12 — blind to a mis-scaled source)"
    )

    out["pairs"] = {
        k: abs(out["th_own"][(k + 1) % 4] - out["th_own"][k]) / out["th_own"][k] for k in range(4)
    }
    out["q_pairs"] = {
        k: abs(out["q_own"][(k + 1) % 4] - out["q_own"][k]) / out["q_own"][k] for k in range(4)
    }
    out["far"] = {k: abs(out["th_far"][k] - out["th_own"][k]) / out["th_own"][k] for k in range(4)}
    for k in range(4):
        _say(
            f"MAT15 (C4) k={k}->{(k + 1) % 4}: theta 10 g pair {out['pairs'][k] * 100:.4f}% "
            f"(bound {C4_COVARIANCE_BAND * 100:.1f}%, predicted ≤ 0.1%);  q_h 10 g pair "
            f"{out['q_pairs'][k] * 100:.4f}% (printed);  control (f) far-side "
            f"{out['far'][k] * 100:.4f}% (asserted > {F_CONTROL_FLOOR * 100:.1f}%, predicted 40-80%)"
        )
    worst = max(out["E"], key=out["E"].get)
    _say(
        f"MAT15 (E) worst balance over {len(out['E'])} runs = {out['E'][worst]:.3e} on {worst} "
        f"(bound {E_BAND:.0e}); max|Im theta| = {tmax_imag:.3e} (predicted 0)"
    )
    _say(f"MAT15 elapsed: EM {t_em:.1f} s, total {time.perf_counter() - t0:.1f} s")
    if return_fields:
        out.update(
            msh=msh,
            cell_tags=cell_tags,
            sub=sub,
            ft=ft,
            by_k=by_k,
            centres=centres,
            q_h=q_h_by_k,
            T_star=t_star_by_k,
            sweep=sweep,
            solves=solves,
        )
    return out


@pytest.fixture(scope="module")
def heated():
    return _run_heated()


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_b_source_transfer_matches_em_power(heated, k):
    assert heated["B"][k] <= B_BAND, f"(B) drive {k}: {heated['B'][k]:.3e} > {B_BAND:.0e}"


@complex_only
def test_b_submesh_volume_matches_record(heated):
    assert heated["err_V"] <= V_BAND, f"V_sub rel {heated['err_V']:.3e} > {V_BAND:.0e}"


@complex_only
def test_e_every_run_balances(heated):
    for name, r in heated["E"].items():
        assert r <= E_BAND, f"(E) run {name}: {r:.3e} > {E_BAND:.0e}"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_e_prime_stored_plus_surface_equals_em_deposited(heated, k):
    assert heated["Ep"][k] <= E_BAND, f"(E') drive {k}: {heated['Ep'][k]:.3e} > {E_BAND:.0e}"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_s2_first_step_slope_is_mean_sar_over_c(heated, k):
    assert heated["S2"][k] <= S2_BAND, f"(S2) drive {k}: {heated['S2'][k]:.4e} > {S2_BAND:.0e}"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_c4_ten_gram_theta_follows_the_drive(heated, k):
    r = heated["pairs"][k]
    assert r <= C4_COVARIANCE_BAND, f"(C4) pair {k}: {r * 100:.4f}% > 5 %"


@complex_only
def test_control_d_doubled_source_misses_the_balance(heated):
    assert (
        heated["ctl_d"] > D_CONTROL_FLOOR
    ), f"control (d) doubled source reads (E') {heated['ctl_d']:.3e} — did not miss"


@complex_only
@pytest.mark.parametrize("k", [0, 1, 2, 3])
def test_control_f_far_side_ball_does_not_agree(heated, k):
    r = heated["far"][k]
    assert r > F_CONTROL_FLOOR, f"control (f) drive {k}: far-side {r * 100:.4f}% inside 5 %"
