"""Example (`EX-81`): the phantom's own cylinder against its closed form — uniform
heating, Robin walls, FEM beside the double eigen-series, steady and at t*.

`mat:4` / `6` / `7` / `8` compare heat *balls* with a one-index series and `mat:9`
heats this very cylinder (a = 0.03 m, 2L = 0.08 m, the 0.0025 m phantom rung, 58 866
cells) under a coil source with no comparand.  Here the source is uniform q and the
comparand is `MAT-15` step 2(a)'s two-index eigen-series, so the FEM error can be
read pointwise, in L2, and mapped onto the faceted lateral wall and the corners.

Everything — the fixture body, every band, every record, the series — comes from
``tests/validation/test_mat15_finite_cylinder_fem.py`` (`MAT-15` step 2(b)) and its
imports: ``_run_transient(return_fields=True)`` is that gate's `transient` fixture
body lifted to module level (rule (a) lift, disclosed; the transient nodes were
re-gated green in the same slot).  Nothing is restated.

Asserted, each printed beside its record: cells = ``FINEST_CELLS``; (T) CN interior,
(E) balance and (inf) BE-vs-steady against ``20261002T093610Z_MAT-15-step2b-W5-transient.log:1875-1883``;
the steady interior and the Gauss balance (composed from the gate's ``_exact`` /
``_rel``) against ``20261002T093152Z_MAT-15-step2b-W4-steady.log:7174-7176``;
negative control (c2) Dirichlet radial roots misses ``INTERIOR_BAND``.
Scope: the 0.0025 rung, P1, uniform q, Robin walls, steady and CN at t*.  Not claimed:
the LSQ rate, EM, perfusion, degree 2, `ANS-9`.

Run through the example runner (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:10 -n 2 -t 480

Output: ``examples/materials/paraview_output/materials_10_finite_cylinder_robin_series_combined.xdmf``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.setup_figure import write_setup_figure

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.validation import test_mat12_transient_sar_heating_sphere as m12  # noqa: E402
from tests.validation import test_mat13_transient_perfusion_sphere as m13  # noqa: E402
from tests.validation import test_mat15_finite_cylinder_fem as m15b  # noqa: E402
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    A,
    H_GEL,
    T_INF,
    T_STAR,
    T_R,
)
from tests.validation.test_mat15_finite_cylinder_closed_form import KW, L, Q  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_fem import (  # noqa: E402
    E_BAND,
    FINEST_CELLS,
    INF_BAND,
    INTERIOR_BAND,
    RUNGS,
    S2_GAUSS_BAND,
    V_CAD,
    _exact,
    _rel,
    _run_transient,
    _steady_ex,
    _trans_ex,
)
from fem_em_solver.utils.analytical import robin_finite_cylinder_steady_temperature  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_10_finite_cylinder_robin_series"
W4 = "20261002T093152Z_MAT-15-step2b-W4-steady.log"
W5 = "20261002T093610Z_MAT-15-step2b-W5-transient.log"

RADIAL_FRAC = (0.0, 0.25, 0.5, 0.75, 0.9)  # r / a at z = 0
AXIAL_FRAC = (0.25, 0.5, 0.75, 0.9)  # z / L at r = 0


def _at(T, pts):
    """theta_h = T - T_inf at points, real part, NaN where no rank owns the point."""
    vals, ok = evaluate_vector_field_parallel(T, pts)
    out = []
    for i in range(pts.shape[0]):
        out.append(float(np.real(np.ravel(vals[i])[0])) - T_INF if bool(ok[i]) else float("nan"))
    return np.array(out)


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    structural = {
        "m15b.INTERIOR_BAND is m13.INTERIOR_BAND": m15b.INTERIOR_BAND is m13.INTERIOR_BAND,
        "m15b.E_BAND is m12.E_BAND": m15b.E_BAND is m12.E_BAND,
        "m15b.INF_BAND is m12.INF_BAND": m15b.INF_BAND is m12.INF_BAND,
    }
    for name, ok in structural.items():
        assert ok, f"structural import check failed: {name}"
    if rank0:
        print("=" * 72)
        print("EX-81 — the phantom's cylinder vs its double eigen-series (MAT-15 step 2(b), re-executed)")
        print("=" * 72)
        print("[structural] " + "; ".join(f"{k}: {v}" for k, v in structural.items()))
        print(f"[records] transient anchors beside {W5}:1875-1883; steady beside {W4}:7174-7176",
              flush=True)

    r = _run_transient(return_fields=True)
    sub, ft, v_sub, st, res, ex_t, chi = (
        r["sub"], r["ft"], r["v_sub"], r["st"], r["res"], r["ex_t"], r["chi"])

    # EX-57 setup figure on the thermal submesh; no-op unless FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        sub,
        None,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={0: "gel phantom (thermal domain; Robin on every face)"},
        slice_normal=(0.0, 1.0, 0.0),
        title="mat:10 -- gel phantom cylinder, uniform heating, Robin walls (MAT-15 2b)",
        comm=comm,
    )

    # ---- steady anchors composed from the gate's own helpers ----------------
    ex_s = _exact(sub, _steady_ex())
    th_s = st.temperature - T_INF
    steady_interior = _rel(sub, th_s, ex_s, chi)
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    conv = float(np.real(_allsum(comm, H_GEL * th_s * ds)))
    gauss = abs(conv - Q * v_sub) / (Q * v_sub)

    cells = r["n"]
    anchors = [
        ("cells", f"{cells}", f"{FINEST_CELLS}", cells == FINEST_CELLS),
        ("(T) CN interior", f"{r['T']:.6e}", f"3.494532e-04 <= {INTERIOR_BAND:.0e}",
         r["T"] <= INTERIOR_BAND),
        ("(E) worst balance", f"{max(r['E'].values()):.3e}", f"7.234e-14 <= {E_BAND:.0e}",
         max(r["E"].values()) <= E_BAND),
        ("(inf) BE vs steady", f"{r['inf']:.3e}", f"2.204e-14 <= {INF_BAND:.0e}",
         r["inf"] <= INF_BAND),
        ("max|Im theta|", f"{r['imag']:.3e}", "0", r["imag"] == 0.0),
        ("steady interior", f"{steady_interior:.6e}", f"7.373685e-04 <= {INTERIOR_BAND:.0e}",
         steady_interior <= INTERIOR_BAND),
        ("steady Gauss", f"{gauss:.3e}", f"1.326e-14 <= {S2_GAUSS_BAND:.0e}",
         gauss <= S2_GAUSS_BAND),
    ]
    ctl = r["c2"]
    if rank0:
        print()
        for name, got, rec, ok in anchors:
            print(f"[anchor] {name}: {got}  (record {rec})  {'ok' if ok else 'MISS'}")
        print(f"[control] (c2) Dirichlet radial roots interior miss = {ctl:.5e} > "
              f"{INTERIOR_BAND:.0e} ({ctl / INTERIOR_BAND:.1f}x the band; W5 record 7.7946e-01)",
              flush=True)
    for name, got, rec, ok in anchors:
        assert ok, f"anchor {name}: {got} vs record {rec}"
    assert f"{r['T']:.6e}" == "3.494532e-04", f"(T) digits off W5: {r['T']:.6e}"
    assert f"{steady_interior:.6e}" == "7.373685e-04", f"steady digits off W4: {steady_interior:.6e}"
    assert f"{ctl:.4e}" == "7.7946e-01", f"(c2) digits off W5: {ctl:.4e}"
    assert ctl > INTERIOR_BAND, f"control (c2) did not miss: {ctl:.4e}"

    # ---- printed (no band) ---------------------------------------------------
    vdef = v_sub / V_CAD - 1.0
    pts0 = np.zeros((1, 3))
    th_c = float(_at(st.temperature, pts0)[0])
    th_c_ex = float(robin_finite_cylinder_steady_temperature(0.0, 0.0, t_inf=0.0, **KW))
    th_ct = float(_at(res.temperature, pts0)[0])
    fn_t = _trans_ex(T_STAR)[0]
    fn_s = _steady_ex()[0]
    ex_ct = float(fn_t(np.array([0.0]), np.array([0.0]))[0])
    if rank0:
        print(f"\n[printed] V_sub deficit vs V_CAD = {vdef:+.4e} (W4 -8.6995e-4)")
        print(f"[printed] theta_h(0) = {th_c:.6f} K beside W4 2.980951 K and series {th_c_ex:.6f} K "
              f"(rel {th_c / th_c_ex - 1.0:+.3e})")
        print(f"[printed] theta_h(0, t*)/theta_h(0, inf) = {th_ct / th_c:.6f} beside series "
              f"{ex_ct / th_c_ex:.6f}  (t*/t_R = {T_STAR / T_R:.3f})")

    def profile(label, pts, coord):
        h_s = _at(st.temperature, pts)
        h_t = _at(res.temperature, pts)
        rr = np.hypot(pts[:, 0], pts[:, 1])
        zz = pts[:, 2]
        s_s = fn_s(rr, zz)
        s_t = fn_t(rr, zz)
        if rank0:
            print(f"[printed] {label}: steady theta_h | series | rel ; t* theta_h | series | rel")
            for i, c in enumerate(coord):
                print(f"  {c:5.2f}  {h_s[i]:.6f} | {s_s[i]:.6f} | {h_s[i] / s_s[i] - 1.0:+.3e}   "
                      f"{h_t[i]:.6f} | {s_t[i]:.6f} | {h_t[i] / s_t[i] - 1.0:+.3e}")

    radial = np.array([[f * A, 0.0, 0.0] for f in RADIAL_FRAC])
    axial = np.array([[0.0, 0.0, f * L] for f in AXIAL_FRAC])
    profile("radial profile at z = 0, r/a", radial, RADIAL_FRAC)
    profile("axial profile at r = 0, z/L", axial, AXIAL_FRAC)

    # ---- combined XDMF on the thermal submesh --------------------------------
    V1 = st.temperature.function_space
    theta_h_steady = fem.Function(V1, name="theta_h_steady")
    theta_h_steady.x.array[:] = st.temperature.x.array - T_INF
    theta_series_steady = fem.Function(V1, name="theta_series_steady")
    theta_series_steady.interpolate(ex_s)
    theta_err_steady = fem.Function(V1, name="theta_err_steady")
    theta_err_steady.x.array[:] = theta_h_steady.x.array - theta_series_steady.x.array
    theta_h_tstar = fem.Function(V1, name="theta_h_tstar")
    theta_h_tstar.x.array[:] = res.temperature.x.array - T_INF
    theta_series_tstar = fem.Function(V1, name="theta_series_tstar")
    theta_series_tstar.interpolate(ex_t)
    theta_err_tstar = fem.Function(V1, name="theta_err_tstar")
    theta_err_tstar.x.array[:] = theta_h_tstar.x.array - theta_series_tstar.x.array
    chi.name = "interior_mask"
    fields = {}
    for f in (theta_h_steady, theta_series_steady, theta_err_steady,
              theta_h_tstar, theta_series_tstar, theta_err_tstar, chi):
        f.x.scatter_forward()
        fields[f.name] = f

    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None, fields, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    if rank0:
        print(f"\n[output] {xdmf_path.name}: P1 theta_h_steady, theta_series_steady, theta_err_steady, "
              f"theta_h_tstar, theta_series_tstar, theta_err_tstar (K, t* = {T_STAR:.1f} s) and DG0 "
              f"interior_mask (r < 0.9 a, |z| < 0.9 L)")
        print(f"[scope] rung h_p = {RUNGS[-1]}, P1, uniform q = {Q:g} W/m^3, Robin walls")
        print("EX-81 anchors and control: PASS")


if __name__ == "__main__":
    main()
