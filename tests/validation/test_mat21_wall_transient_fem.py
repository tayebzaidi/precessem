"""`MAT-21` step 20c — the walled box's Crank–Nicolson transient against the
exact two-layer slab transient (§9 item 682, §10 P48c; brief
`docs/planning/chunks/MAT-21.md` "§9 item 682").

Geometry: item 670's **1D-reduced** walled box (`test_mat21_wall_fem`'s
``_walled_box`` with the wall on ±z only, Robin on z = ±(C + d) only, the four
side faces natural Neumann), rung m = 4 (52, 32, 12) = 119 808 cells.  There
the 1D two-layer slab series is exact for any A, B.  Mesh helpers are imported
from item 670's module (never edited); k is item 670's DG0 field, ρc a local
DG0 field in ``_materials``' pattern (gel ``RHO_C``, wall ``RHO_C_W``), the
source Q·χ_gel switched on at t = 0.

Solver: ``solve_transient_conduction`` θ = 0.5 with two backward-Euler start-up
half-steps, Δt = t*/N_MAIN, 8·N_MAIN steps, snapshots at t* = 0.2 T_C (step
N_MAIN) and t₂ = 8 t* = 1.6 T_C (step 8 N_MAIN; near item 673's printed (ii)
maximum at 1.62 T_C).

Comparand: item 677's `two_layer_slab_transient_temperature` (400 modes) in z,
interpolated into P2 in chunks; errors on θ = T − T∞.

Asserted (imported gates, unmoved): gel-interior (|x|<0.9A, |y|<0.9B,
|z|<0.9C — item 670's (S1) interior selection) rel L2 ≤ ``INTERIOR_BAND`` at
both snapshots; the max per-step balance ≤ ``E_BAND``; tag volumes ≤ item
670's ``_vol_bar`` (rule (f) bar).  Negative control (asserted): the
massless-wall series (ρ_wc_w = 1e-9 × nominal, which item 673/677's (K0) proved
equal to the bare slab at h_eff) misses the FEM at t₂ by more than 2× the FEM's
own interior error at t₂.

``K_W`` / ``RHO_C_W`` are nominal PMMA values, citation owed (dashboard W3):
the identity tested here (FEM = the exact series for the same coefficients)
does not depend on their values.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from dolfinx import fem
from dolfinx import mesh as dmesh

from fem_em_solver.thermal import solve_transient_conduction
from fem_em_solver.utils.analytical import (
    two_layer_slab_eigenvalues,
    two_layer_slab_transient_temperature,
)

from tests.validation.test_mat10_sar_to_heat_sphere import _say
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, N_MAIN
from tests.validation.test_mat15_finite_cylinder_closed_form import H, K, Q, RHO_C, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import _rel
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX, T_C
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, RHO_C_W
from tests.validation.test_mat21_wall_fem import (
    BASE,
    GEL,
    ROBIN,
    WALL,
    _centre,
    _materials,
    _robin_facets,
    _vol_bar,
    _volumes,
    _walled_box,
)
from tests.validation.test_mat21_wall_transient_closed_form import N_SERIES, T_STAR

RUNG = 4
N2 = 8 * N_MAIN
T2 = N2 * (T_STAR / N_MAIN)
CTL_FACTOR = 2.0
MASSLESS = 1e-9  # (K0)'s ρ_wc_w scale
INTERP_CHUNK = 5000
SLAB = dict(half_thickness=C_BOX, wall=D_WALL, k=K, k_w=K_W, rho_c=RHO_C)


def _lam(rho_c_w):
    return two_layer_slab_eigenvalues(h=H, n=N_SERIES, rho_c_w=rho_c_w, **SLAB)


def _series(z, t, rho_c_w, lam):
    return two_layer_slab_transient_temperature(z, t, q=Q, h=H, t_inf=0.0, rho_c_w=rho_c_w,
                                                n_terms=N_SERIES, eigenvalues=lam, **SLAB)


def _exact_t(msh, t, rho_c_w, lam):
    V2 = fem.functionspace(msh, ("Lagrange", 2))
    g = fem.Function(V2)

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, INTERP_CHUNK):
            out[s:s + INTERP_CHUNK] = _series(x[2, s:s + INTERP_CHUNK], t, rho_c_w, lam)
        return out

    g.interpolate(ev)
    return g


def _rho_c_field(msh):
    """ρc as DG0 in item 670's ``_materials`` pattern (copied, not edited there)."""
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    mids = dmesh.compute_midpoints(msh, tdim, np.arange(n, dtype=np.int32))
    gel = ((np.abs(mids[:, 0]) < A_BOX) & (np.abs(mids[:, 1]) < B_BOX)
           & (np.abs(mids[:, 2]) < C_BOX))
    D0 = fem.functionspace(msh, ("DG", 0))
    dofs = D0.dofmap.list[:n, 0]
    rcf = fem.Function(D0)
    rcf.x.array[dofs] = np.where(gel, RHO_C, RHO_C_W)
    return rcf


@pytest.fixture(scope="module")
def walled():
    t0 = time.perf_counter()
    # Constant identity: item 670's K, H (step 0's source) are MAT-15's K_GEL, H_GEL.
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    msh, ns, nw = _walled_box(RUNG, (False, False, True))
    n = msh.topology.index_map(3).size_global
    chi, kf, inter, ct = _materials(msh)
    rcf = _rho_c_field(msh)
    ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
    vg, vw = _volumes(msh, ct)
    v_gel = 8.0 * A_BOX * B_BOX * C_BOX
    v_wall = 4.0 * A_BOX * B_BOX * 2.0 * D_WALL
    out = {"n": n, "dvg": abs(vg / v_gel - 1.0), "dvw": abs(vw / v_wall - 1.0)}
    _say(f"\nMAT21-20c [W1] 1D-reduced walled box m = {RUNG} cells {tuple(ns)} N_w = {nw} "
         f"(BASE {BASE}, tags gel {GEL} / wall {WALL}): n = {n}; V_gel rel dev {out['dvg']:.3e}, "
         f"V_wall rel dev {out['dvw']:.3e} (bar {_vol_bar(n):.3e}); k {K} / {K_W}, rho c "
         f"{RHO_C} / {RHO_C_W}, h {H}, q {Q}, d {D_WALL}; T_C {T_C:.3f} s, t* {T_STAR:.3f} s, "
         f"t2 {T2:.3f} s = {T2 / T_C:.4f} T_C; dt {T_STAR / N_MAIN:.4f} s, {N2} CN steps")
    tl = time.perf_counter()
    res = solve_transient_conduction(
        msh, ft, k=kf, rho_c=rcf, q=Q * chi, h=H_GEL, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN, theta=0.5, startup_half_steps=2, dt=T_STAR / N_MAIN, n_steps=N2,
        snapshot_steps=(N_MAIN, N2))
    _say(f"MAT21-20c [W1] CN solve {time.perf_counter() - tl:.1f} s; source power "
         f"{res.source_power:.9e} W (Q V_gel {Q * v_gel:.9e}); max balance "
         f"{res.max_balance_residual:.3e} over {len(res.balance_residuals)} (half-)steps")
    out["bal"] = res.max_balance_residual
    lam_w = _lam(RHO_C_W)
    lam_0 = _lam(MASSLESS * RHO_C_W)
    out["snap"] = {}
    for step, t in ((N_MAIN, T_STAR), (N2, T2)):
        th = res.snapshots[step] - T_INF
        ex = _exact_t(msh, t, RHO_C_W, lam_w)
        ex0 = _exact_t(msh, t, MASSLESS * RHO_C_W, lam_0)
        e_int = _rel(msh, th, ex, inter)
        e_whole = _rel(msh, th, ex)
        miss0 = _rel(msh, th, ex0, inter)
        c_h = _centre(res.snapshots[step]) - T_INF
        c_ex = float(_series(0.0, t, RHO_C_W, lam_w)[0])
        c_0 = float(_series(0.0, t, MASSLESS * RHO_C_W, lam_0)[0])
        out["snap"][step] = {"e_int": e_int, "e_whole": e_whole, "miss0": miss0,
                             "c_h": c_h, "c_ex": c_ex, "c_0": c_0}
        _say(f"MAT21-20c [W1] step {step} t = {t:.3f} s ({t / T_C:.4f} T_C): interior rel L2 vs "
             f"two-layer series {e_int:.6e} (bound {INTERIOR_BAND:.0e}); whole {e_whole:.6e} "
             f"(printed); massless-wall series vs FEM, interior {miss0:.6e}")
        _say(f"MAT21-20c [W1] step {step} centre theta: FEM {c_h:.9f} K, two-layer series "
             f"{c_ex:.9f} K (rel {c_h / c_ex - 1.0:+.3e}), massless-wall series {c_0:.9f} K "
             f"(series diff {c_0 / c_ex - 1.0:+.4e})")
    _say(f"MAT21-20c [W1] elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_volume_identity(walled):
    assert walled["dvg"] <= _vol_bar(walled["n"]), walled
    assert walled["dvw"] <= _vol_bar(walled["n"]), walled


def test_balance(walled):
    assert walled["bal"] <= E_BAND, walled["bal"]


def test_interior_both_snapshots(walled):
    errs = {s: v["e_int"] for s, v in walled["snap"].items()}
    assert all(e <= INTERIOR_BAND for e in errs.values()), errs


def test_control_massless_wall_misses_at_t2(walled):
    s = walled["snap"][N2]
    _say(f"MAT21-20c (control) at t2: massless-wall series misses the FEM by {s['miss0']:.6e} "
         f"(interior rel L2) vs {CTL_FACTOR:.0f} x FEM interior error {CTL_FACTOR * s['e_int']:.6e} "
         f"(separation {s['miss0'] / s['e_int']:.2f}x)")
    assert s["miss0"] > CTL_FACTOR * s["e_int"], s
