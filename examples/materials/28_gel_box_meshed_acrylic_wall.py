"""Example (`EX-159`): the F2182 gel box with a meshed acrylic wall -- the wall as
cells, converging at order 2 to the two-layer slab, and the bare slab missing.

`mat:12` has a spherical wall on the ball; `mat:14` puts a wall into h_eff with
no wall cells. This example meshes a 6 mm acrylic wall on the *box*: DolfinX
``create_box`` over the outer box, every rank's geometry remapped so grid planes
fall exactly at the gel faces, wall cells carrying a DG0 k_w and no source.
It runs `MAT-21` step 20 / P48b on the 1D-reduced geometry (wall on +-z, Robin on
z = +-(C + d) only, Neumann sides), which is exactly item 660's two-layer slab.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Mesh builder,
materials, Robin tags, volumes, the exact two-layer interpolant, the solve
wrapper and the rung table are imported by name from
``tests/validation/test_mat21_wall_fem.py`` (`w`); the bands from the modules
that gate imports them from. The gate's ``reduced`` fixture is a pytest fixture
and is not called: its loop is written out here with every call it makes, plus
``np.real`` on each integral (the group sources the complex build). No gate
edit, no ``src/`` edit.

Asserted: (V) gel and wall volumes against analytic at the rule-(f) bar;
(S1) whole-domain rel L2 LSQ rate >= ``RATE_BAR`` and finest gel-interior rel L2
<= ``INTERIOR_BAND``; (E) balance <= ``E_BAND`` on every rung; (c-h) the bare
Robin slab at h (no wall) misses the FEM centre by > ``INTERIOR_BAND``.
Printed only: the six-wall box at m = 4 beside the bare box at h_eff and the
triple series ("lateral wall conduction is not tested"), the wall's relative
resistance eps_w, cells and seconds.

Run it::

    ./run_examples.sh -e mat:28 -n 2 -t 300

Output: ``examples/materials/paraview_output/
materials_28_gel_box_meshed_acrylic_wall_combined.xdmf`` (theta on the m = 4
reduced box with the DG0 k field and the error against the two-layer series).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import solve_steady_conduction  # noqa: E402
from fem_em_solver.utils.analytical import (  # noqa: E402
    robin_box_steady_temperature,
    wall_effective_heat_transfer,
)

from tests.validation.test_mat10_sar_to_heat_sphere import _allsum  # noqa: E402
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_closed_form import T_INF  # noqa: E402
from tests.validation.test_mat15_finite_cylinder_fem import RATE_BAR, _rel  # noqa: E402
from tests.validation.test_mat21_box_closed_form import (  # noqa: E402
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    H,
    K,
    Q,
    SERIES_TERMS,
)
from tests.validation.test_mat21_wall_closed_form import D_WALL, K_W, _bare_series  # noqa: E402
from tests.validation.test_mat21_wall_fem import (  # noqa: E402
    BASE,
    GEL,
    ROBIN,
    RUNGS,
    W2_RUNG,
    WALL,
    _centre,
    _exact_z,
    _materials,
    _robin_facets,
    _slab_theta,
    _solve,
    _vol_bar,
    _volumes,
    _walled_box,
)
from dolfinx import mesh as dmesh  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_28_gel_box_meshed_acrylic_wall"

#: Gate records, real build -- `20261027T200231Z_MAT-21-step20-W1-fem.log:53-57`,
#: `:67`, and `20261027T200257Z_MAT-21-step20-W2-fem.log:52-53`.  Printed beside
#: this run, a difference reported and never asserted.
RECORD = {
    "rate": 1.9896,
    "interior": 1.190636e-3,
    "balance": 2.7e-14,
    "miss": 2.2967e-1,
    "six_wall": 14.5339,
    "bare_heff": 14.6001,
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example expects the complex DolfinX build the runner's `mat:` group "
            "sources (source /usr/local/bin/dolfinx-complex-mode)"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-159 -- gel box with a meshed acrylic wall (MAT-21 step 20 / P48b)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: RUNGS {RUNGS}, RATE_BAR {RATE_BAR}, INTERIOR_BAND "
               f"{INTERIOR_BAND:.0e}, E_BAND {E_BAND:.0e}; (A, B, C) = {HALF}, d = {D_WALL}, "
               f"k = {K}, k_w = {K_W}, h = {H}, q = {Q}")

    v_gel = 8.0 * A_BOX * B_BOX * C_BOX
    v_wall = 4.0 * A_BOX * B_BOX * 2.0 * D_WALL
    rungs = []
    keep = {}
    for m in RUNGS:
        tl = time.perf_counter()
        msh, ns, nw = _walled_box(m, (False, False, True))
        n = msh.topology.index_map(3).size_global
        chi, kf, inter, ct = _materials(msh)
        if m == W2_RUNG:
            write_setup_figure(
                msh,
                ct,
                FIGURE_DIR / f"{BASENAME}_setup.png",
                region_names={GEL: "gel (k = 0.6 W/m/K, uniform q)",
                              WALL: "acrylic wall, 6 mm (k_w, no source)"},
                translucent_tags=(GEL,),
                slice_normal=(0.0, 1.0, 0.0),
                slice_origin=(0.0, 0.0, 0.0),
                title="mat:28 -- 1D-reduced gel box, wall on +-z, Robin on z = +-(C + d), m = 4",
                comm=comm,
            )
        ft = _robin_facets(msh, (A_BOX, B_BOX, C_BOX + D_WALL), (2,))
        vg, vw = _volumes(msh, ct)
        st = _solve(msh, ft, kf, chi)
        th = st.temperature
        ex = _exact_z(msh)
        e = _rel(msh, th, ex)
        ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN)
        conv = float(np.real(_allsum(msh.comm, H * th * ds)))
        bal = abs(conv - Q * v_gel) / (Q * v_gel)
        hm = ((vg + vw) / n) ** (1.0 / 3.0)
        rec = {"m": m, "n": n, "h": hm, "e": e, "bal": bal,
               "dvg": abs(vg / v_gel - 1.0), "dvw": abs(vw / v_wall - 1.0)}
        _say(comm, f"[W1] m = {m} cells {tuple(ns)} N_w = {nw}: n = {n}; V_gel rel dev "
                   f"{rec['dvg']:.3e}, V_wall rel dev {rec['dvw']:.3e} (bar {_vol_bar(n):.3e}); "
                   f"h = {hm:.6e} m; whole rel L2 vs two-layer slab {e:.6e}; balance {bal:.3e}; "
                   f"{time.perf_counter() - tl:.1f} s")
        rungs.append(rec)
        if m == W2_RUNG:
            keep = {"msh": msh, "th": th, "kf": kf, "ct": ct, "ex": ex}
        if m == RUNGS[-1]:
            interior = _rel(msh, th, ex, inter)
            centre = _centre(th)
            centre_ex = float(np.real(_slab_theta(0.0)))
    hs = np.array([r["h"] for r in rungs])
    es = np.array([r["e"] for r in rungs])
    slope = float(np.polyfit(np.log(hs), np.log(es), 1)[0])
    pw = [float(np.log(es[i] / es[i + 1]) / np.log(hs[i] / hs[i + 1])) for i in range(len(hs) - 1)]
    _say(comm, f"[S1] LSQ slope {slope:.4f} (record {RECORD['rate']}, bar {RATE_BAR}); pairwise "
               + ", ".join(f"{p:.4f}" for p in pw))
    _say(comm, f"[S1] finest interior rel L2 {interior:.6e} (record {RECORD['interior']:.6e}, "
               f"bound {INTERIOR_BAND:.0e}); centre theta_h {centre:.9f} K vs two-layer "
               f"{centre_ex:.9f} K (rel {centre / centre_ex - 1.0:+.3e})")
    worst_bal = max(r["bal"] for r in rungs)
    _say(comm, f"[E] worst balance {worst_bal:.3e} (record <= {RECORD['balance']:.1e}, "
               f"bound {E_BAND:.0e})")

    # ---- control: the bare slab at h --------------------------------------------------
    bare = float(np.real(_bare_series(0.0, H)[0])) - T_INF
    miss = abs(centre - bare) / abs(bare)
    _say(comm, f"[c-h] FEM centre {centre:.9f} K vs bare Robin slab at h (no wall) {bare:.9f} K: "
               f"rel miss {miss:.4e} (record {RECORD['miss']:.4e}; asserted > {INTERIOR_BAND:.0e})")

    # ---- six-wall box, one rung (print only) -------------------------------------------
    t6 = time.perf_counter()
    m = W2_RUNG
    msh6, ns6, nw6 = _walled_box(m, (True, True, True))
    n6 = msh6.topology.index_map(3).size_global
    chi6, kf6, _i6, ct6 = _materials(msh6)
    outer = tuple(L + D_WALL for L in HALF)
    ft6 = _robin_facets(msh6, outer, (0, 1, 2))
    vg6, vw6 = _volumes(msh6, ct6)
    v_wall6 = 8.0 * np.prod(outer) - v_gel
    st6 = _solve(msh6, ft6, kf6, chi6)
    ds6 = ufl.Measure("ds", domain=msh6, subdomain_data=ft6)(ROBIN)
    conv6 = float(np.real(_allsum(msh6.comm, H * st6.temperature * ds6)))
    bal6 = abs(conv6 - Q * v_gel) / (Q * v_gel)
    th_wall = _centre(st6.temperature)
    h_eff = wall_effective_heat_transfer(H, D_WALL, K_W)
    ng = [b * m for b in BASE]
    mb = dmesh.create_box(comm, [-np.array(HALF), np.array(HALF)], ng,
                          cell_type=dmesh.CellType.tetrahedron)
    ftb = _robin_facets(mb, HALF, (0, 1, 2))
    stb = solve_steady_conduction(mb, ftb, k=K, q=Q, h=h_eff, t_inf=0.0, robin_tag=ROBIN, degree=1)
    th_bare = _centre(stb.temperature)
    kw = dict(half_widths=HALF, k=K, q=Q, t_inf=0.0, n_terms=SERIES_TERMS)
    ser_eff = float(np.real(robin_box_steady_temperature(0.0, 0.0, 0.0, h=h_eff, **kw)))
    ser_h = float(np.real(robin_box_steady_temperature(0.0, 0.0, 0.0, h=H, **kw)))
    dvg6 = abs(vg6 / v_gel - 1.0)
    dvw6 = abs(vw6 / v_wall6 - 1.0)
    eps_w = th_bare / th_wall - 1.0
    _say(comm, f"\n[W2] six-wall box m = {m} cells {tuple(ns6)} N_w = {nw6}: n = {n6}; V_gel rel dev "
               f"{dvg6:.3e}, V_wall rel dev {dvw6:.3e} (bar {_vol_bar(n6):.3e}); balance {bal6:.3e}")
    _say(comm, f"[W2] (print; lateral wall conduction is not tested) centre rise: six-wall meshed FEM "
               f"{th_wall:.6f} K (record {RECORD['six_wall']}); bare gel box FEM at h_eff {h_eff:.6f} "
               f"{th_bare:.6f} K (record {RECORD['bare_heff']}); triple series at h_eff {ser_eff:.6f} K, "
               f"at h {ser_h:.6f} K; {time.perf_counter() - t6:.1f} s")
    _say(comm, f"[W2] eps_w (h_eff-only miss, bare-at-h_eff / six-wall - 1) = {eps_w:+.4e}; "
               f"wall Biot h*d/k_w = {H * D_WALL / K_W:.4e}")

    # ---- output: m = 4 reduced box ----------------------------------------------------
    msh4, th4, kf4, ct4, ex4 = (keep[k] for k in ("msh", "th", "kf", "ct", "ex"))
    V = th4.function_space
    th_out = fem.Function(V, name="theta_h")
    th_out.x.array[:] = th4.x.array - T_INF
    th_out.x.scatter_forward()
    th_ser = fem.Function(V, name="theta_two_layer_series")
    err = fem.Function(V, name="theta_error")
    ex_p1 = fem.Function(V)
    ex_p1.interpolate(lambda x: np.real(_slab_theta(x[2])))
    th_ser.x.array[:] = ex_p1.x.array
    th_ser.x.scatter_forward()
    err.x.array[:] = th_out.x.array - ex_p1.x.array
    err.x.scatter_forward()
    kf4.name = "k_field_DG0"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh4, ct4,
        {"theta_h": th_out, "theta_two_layer_series": th_ser, "theta_error": err,
         "k_field_DG0": kf4},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- m = 4 reduced box: theta_h, two-layer series, error, DG0 k")

    # ---- assertions (the gate's) ------------------------------------------------------
    for r in rungs:
        assert r["dvg"] <= _vol_bar(r["n"]) and r["dvw"] <= _vol_bar(r["n"]), f"(V) {r}"
    assert dvg6 <= _vol_bar(n6) and dvw6 <= _vol_bar(n6), f"(V) six-wall {dvg6}, {dvw6}"
    assert slope >= RATE_BAR, f"(S1) LSQ slope {slope:.4f}"
    assert interior <= INTERIOR_BAND, f"(S1) interior {interior:.3e}"
    assert worst_bal <= E_BAND, f"(E) balance {worst_bal:.3e}"
    assert bal6 <= E_BAND, f"(E) six-wall balance {bal6:.3e}"
    assert miss > INTERIOR_BAND, f"(c-h) bare slab at h does not miss: {miss:.3e}"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
