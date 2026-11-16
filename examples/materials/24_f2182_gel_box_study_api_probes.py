"""Example (`EX-130`, `mat:24`): the F2182-sized gel box through the package's own route --
``ThermalBoxGeometry`` -> ``Study.temperature`` -> ``Study.probe`` -- with a probe table at the
six named F2182 positions (P0 centre, P1 / P2 the +-y walls, P3 / P4 the +-x walls, P5 the top face).

The physics is `mat:19`'s uniform steady box (q = Q, Robin on all six faces, m = 4 rung, P1).  The
delivery mechanism is new to the examples: no example calls ``ThermalBoxGeometry`` or ``Study.probe``.
`WF-15` step 5 (``Study.temperature`` on the box) and step 6 (``Study.probe`` on the box) are the
gates; `MAT-21` step 10 (P37) names the six positions.

Everything -- every band, every helper, every position -- is imported from the gate modules
(``tests.package.test_wf15_study_thermal_box`` as ``s5``, ``tests.package.test_wf15_study_probe_thermal``
as ``s6``, ``tests.validation.test_mat21_box_probe_points_fem`` as ``p37`` and the modules they import
from).  Records live only in the example-local ``REC`` dict with their ``log:line``, printed beside the
readings and never asserted by equality.

Asserted (imported bands only): (K) build / ranks / divisions; (S1-m4) through ``Study.temperature``
the cell and dof counts, whole-box L2 at ``.4e`` equal to step 5's record, interior <= ``INTERIOR_BAND``,
Gauss <= ``S2_GAUSS_BAND``; (I-rep) the theta inside ``Study.probe(...).thermal`` equals
``Study.temperature``'s <= ``L0_BAND``; (A-box) the whole-box ball <= ``B_BAND`` on value and volume;
(R-box) the five step-6 balls <= ``INTERIOR_BAND`` of the exact ball mean; (V) the six positions are
valid on the owner route; (M) every inset centre keeps its ball inside the box.  Controls: (c-wall) h
doubled misses (backed by step 6's record, asserted > ``INTERIOR_BAND``, < 1); (c-half) a ball centred
on the +y wall holds only half its nominal volume (predicted, asserted > ``INTERIOR_BAND``, < 1).
Printed only: the six point values against the uniform series, the six inset ball means against the
exact ball mean, (c-mv).

Not claimed: transient or duty probes on the box, the quadratic source, the YAML route, a real probe's
volume or self-heating, an h-rate, EM, a coil, a container wall or perfusion, compliance or an
external check.

Run through the example runner (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:24 -n 2 -t 180

Output: ``examples/materials/paraview_output/materials_24_f2182_gel_box_study_api_probes_combined.xdmf``
and ``..._probe_table.csv``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver import RobinCooling, Study, ThermalBoxGeometry, UniformHeatSource
from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.evaluation import point_owner_cells
from fem_em_solver.post.setup_figure import write_setup_figure
from fem_em_solver.thermal.probe import ProbeBall
from fem_em_solver.utils.analytical import robin_box_steady_temperature

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import tests.package.test_wf15_study_probe_thermal as s6  # noqa: E402
import tests.package.test_wf15_study_thermal_box as s5  # noqa: E402
import tests.validation.test_mat21_box_probe_points_fem as p37  # noqa: E402
from tests.complex_mode import IS_COMPLEX  # noqa: E402
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, T_INF  # noqa: E402
from tests.validation.test_mat21_box_closed_form import HALF, SERIES_TERMS, H, K, Q  # noqa: E402
from tests.validation.test_mat21_box_fem import BASE, _exact, _interior_chi, _steady_ex  # noqa: E402

L0_BAND = s5.L0_BAND
S2_GAUSS_BAND = s5.S2_GAUSS_BAND
_on, _theta_owned_rel, _rel, _max_imag = s5._on, s5._theta_owned_rel, s5._rel, s5._max_imag
A_BOXPROBE, BOX_PROBES, WHOLE_BOX = s6.A_BOXPROBE, s6.BOX_PROBES, s6.WHOLE_BOX
BOX_DIVISIONS, S5_C_WALL = s6.BOX_DIVISIONS, s6.S5_C_WALL
_series_at = s6._series_at
PROBE_PTS, PROBE_LABELS, _probe_vals = p37.PROBE_PTS, p37.PROBE_LABELS, p37._probe_vals

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_24_f2182_gel_box_study_api_probes"

S5W1 = "20261018T112154Z_WF-15-s5-W1-steady.log"
S6W2 = "20261018T123416Z_WF-15-s6-W2-box.log"

# Records, printed beside readings only; never asserted by equality.
REC = {
    "i_dir": 5.918e-16,                       # S5W1 :67
    "s5_int": 3.2236e-03, "s5_gauss": 2.708e-13,  # S5W1 :69
    "s5_cwall": 3.8048e-01,                   # S5W1 :70
    "a_val": 2.495e-13, "a_vol": 2.509e-13,   # S6W2 :63
    "rbox_max": 2.7066e-03,                   # S6W2 :64-69
    "cmv": 6.0750e-02,                        # S6W2 :71 (K)
    "s6_cwall": 3.5712e-01,                   # S6W2 :72
}


def fmt(v: float) -> str:
    """theta at <= 6 significant digits (leak guard matches at 8)."""
    v = float(v)
    return f"{v:.4f}" if abs(v) >= 10.0 else f"{v:.5e}"


class _Clock:
    """Per-phase wall clock, printed on rank 0 with flush, prefix EX-130."""

    def __init__(self):
        self.t = time.perf_counter()
        self.t0 = self.t

    def lap(self, phase):
        now = time.perf_counter()
        if MPI.COMM_WORLD.rank == 0:
            print(f"EX-130 clock {phase}: {now - self.t:.1f} s", flush=True)
        self.t = now


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    def say(msg=""):
        if rank0:
            print(msg, flush=True)

    # ---- (K) ----------------------------------------------------------------------------
    assert IS_COMPLEX, "this example needs the complex DolfinX build (the runner's `mat:` group sources it)"
    assert comm.size == 2, f"comm.size {comm.size}"
    assert BOX_DIVISIONS == tuple(4 * b for b in BASE), (BOX_DIVISIONS, BASE)
    assert K == K_GEL and H == H_GEL, (K, K_GEL, H, H_GEL)
    say("=" * 72)
    say("EX-130 -- the F2182-sized gel box through ThermalBoxGeometry -> Study.temperature -> Study.probe")
    say("=" * 72)
    say(f"[records] beside {S5W1} and {S6W2}; divisions {BOX_DIVISIONS}, ranks {comm.size}, "
        f"SERIES_TERMS {SERIES_TERMS}")
    say(f"[bands] imported: INTERIOR_BAND {INTERIOR_BAND:.0e}, S2_GAUSS_BAND {S2_GAUSS_BAND:.0e}, "
        f"L0_BAND {L0_BAND:.0e}, B_BAND {B_BAND:.0e}")
    clk = _Clock()

    geo = ThermalBoxGeometry(HALF, BOX_DIVISIONS)
    cool = RobinCooling(K_GEL, H_GEL, T_INF)

    # ---- (S1-m4) Study.temperature ------------------------------------------------------
    res = Study(geo, None, UniformHeatSource(Q)).temperature(cool)
    clk.lap("Study.temperature steady m=4")
    msh = res.submesh
    ft = res.facet_tags
    tdim = msh.topology.dim
    imap_c = msh.topology.index_map(tdim)
    n_cells_all = imap_c.size_local + imap_c.num_ghosts
    cell_tags = dmesh.meshtags(
        msh, tdim, np.arange(n_cells_all, dtype=np.int32), np.full(n_cells_all, 1, dtype=np.int32))
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel (650 x 420 x 90 mm)"},
        title="mat:24 -- F2182-sized gel box via ThermalBoxGeometry, uniform q, Robin on all six "
              "faces: z = 0 slice",
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        clip_normal=(0.0, 0.0, -1.0),
        comm=comm,
    )
    clk.lap("setup figure (no-op unless FEM_EM_SETUP_FIGURES=1)")

    cells = int(imap_c.size_global)
    p1_dofs = int(res.temperature.function_space.dofmap.index_map.size_global)
    th_h = res.temperature - T_INF
    ex, npts = _exact(msh, _steady_ex())
    clk.lap(f"steady series interpolation ({npts} P2 points)")
    whole = _rel(msh, th_h, ex)
    chi = _interior_chi(msh)
    interior = _rel(msh, th_h, ex, chi)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN_TAG)
    conv = float(np.real(_allsum(comm, H_GEL * th_h * ds)))
    one = fem.Constant(msh, default_scalar_type(1.0))
    vh = float(np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    gauss = abs(conv - Q * vh) / (Q * vh)
    imag = _max_imag(res.temperature)
    say(f"[anchor] (S1-m4) cells {cells} (record {s5.STEP1_M4_CELLS}); P1 dofs {p1_dofs} "
        f"(record {s5.STEP1_M4_P1_DOFS}); whole L2 {whole:.6e} -> {whole:.4e} (record "
        f"{s5.STEP1_M4_WHOLE}); interior {interior:.6e} <= {INTERIOR_BAND:.0e} (record "
        f"{REC['s5_int']:.4e}, {S5W1}:69); Gauss {gauss:.3e} <= {S2_GAUSS_BAND:.0e} (record "
        f"{REC['s5_gauss']:.3e}); max imag {imag:.1e}")

    # ---- the six points: validity (V) and the point table ---------------------------------
    pts = np.asarray(PROBE_PTS, dtype=float)
    theta_pt, ok_pt = _probe_vals(res.temperature, pts)
    owner_orig, _ = point_owner_cells(msh, pts)
    owner_orig = np.asarray(owner_orig, dtype=np.int64)
    ser_pt = np.asarray(robin_box_steady_temperature(
        pts[:, 0], pts[:, 1], pts[:, 2], half_widths=HALF, k=K, q=Q, h=H, t_inf=0.0,
        n_terms=SERIES_TERMS), dtype=float)
    e_p = np.abs(theta_pt - ser_pt) / np.abs(ser_pt)
    clk.lap("six point values + series")
    say()
    for i in range(len(PROBE_PTS)):
        say(f"[points] {PROBE_LABELS[i]} {PROBE_PTS[i]}: theta_h {fmt(theta_pt[i])} K, series "
            f"{fmt(ser_pt[i])} K, e_p {e_p[i]:.3e}, owner {owner_orig[i]}")
    say(f"[points] max e_p {np.max(e_p):.3e} (predicted ~3e-3, range [1e-4, 1e-2]; in range: "
        f"{1e-4 <= float(np.max(e_p)) <= 1e-2}; printed only)")

    # ---- inset centres (M) -----------------------------------------------------------------
    a = A_BOXPROBE
    normals = [(0, 0, 0), (0, 1, 0), (0, -1, 0), (1, 0, 0), (-1, 0, 0), (0, 0, 1)]
    inset_centres = [tuple(float(p[k]) - a * n[k] for k in range(3))
                     for p, n in zip(PROBE_PTS, normals)]
    for c in inset_centres:
        for x, hw in zip(c, HALF):
            assert abs(x) + a <= hw + 1e-12, ("(M)", c, hw)
    inset_probes = [ProbeBall(c, a) for c in inset_centres]
    half_probe = ProbeBall(PROBE_PTS[1], a)  # centred exactly on the +y wall
    say(f"[anchor] (M) six inset centres, radius {a}: each |c_i| + a <= hw_i + 1e-12 "
        f"(centres {[tuple(round(v, 6) for v in c) for c in inset_centres]})")

    # ---- Study.probe: one call, 13 balls ---------------------------------------------------
    balls = [WHOLE_BOX, *BOX_PROBES, *inset_probes, half_probe]
    t0 = time.perf_counter()
    pres = Study(geo, None, UniformHeatSource(Q)).probe(balls, cool)
    say(f"[time] Study.probe, {len(balls)} balls: {time.perf_counter() - t0:.1f} s")
    clk.lap("Study.probe 13 balls")
    vals, vols, noms = pres.values, pres.volumes, pres.nominal_volumes

    # (I-rep): the re-solve inside probe(...) against Study.temperature's theta
    t_rep = _on(msh, pres.thermal.temperature)
    d_rep = _theta_owned_rel(t_rep, res.temperature, comm)
    say(f"[anchor] (I-rep) theta in Study.probe(...).thermal vs Study.temperature: {d_rep:.3e} "
        f"<= {L0_BAND:.0e} (predicted 0 to 1e-15; record (I-dir) {REC['i_dir']:.3e}, {S5W1}:67)")

    # (A-box)
    one_p = fem.Constant(pres.theta.function_space.mesh, default_scalar_type(1.0))
    pmsh = pres.theta.function_space.mesh
    integral = float(np.real(_allsum(comm, pres.theta * ufl.dx(domain=pmsh))))
    v_h = float(np.real(_allsum(comm, one_p * ufl.dx(domain=pmsh))))
    mean_h = integral / v_h
    a_val = abs(float(vals[0]) - mean_h) / abs(mean_h)
    a_vol = abs(float(vols[0]) / v_h - 1.0)
    say(f"[anchor] (A-box) |P - mean| / |mean| = {a_val:.3e}, |V_probe / V_h - 1| = {a_vol:.3e} "
        f"<= {B_BAND:.0e} (record {REC['a_val']:.3e} / {REC['a_vol']:.3e}, {S6W2}:63)")

    # (R-box)
    mv = Q * a**2 / (10.0 * K_GEL)
    r_rel, point_c = [], []
    for j, c in enumerate(s6.BOX_CENTRES):
        pt = _series_at(c)
        point_c.append(pt)
        cf = pt - mv
        p = float(vals[1 + j])
        r = abs(p / cf - 1.0)
        r_rel.append(r)
        say(f"[anchor] (R-box) c{j} {c}: P {fmt(p)}, series {fmt(pt)}, exact ball mean {fmt(cf)}, "
            f"rel {r:.4e} <= {INTERIOR_BAND:.0e}; volume {float(vols[1 + j]):.6e} / nominal "
            f"{float(noms[1 + j]):.6e}")
    say(f"[anchor] (R-box) max rel {max(r_rel):.4e} (record {REC['rbox_max']:.4e}, {S6W2}:64-69; "
        f"printed only, build differs)")

    # inset ball means, printed only
    say()
    inset_rows = []
    for j, c in enumerate(inset_centres):
        k = 6 + j
        cf = _series_at(c) - mv
        p = float(vals[k])
        r = abs(p / cf - 1.0)
        inset_rows.append((p, cf, r, float(vols[k]) / float(noms[k])))
        say(f"[inset] {PROBE_LABELS[j]} centre {tuple(round(v, 6) for v in c)}: ball mean {fmt(p)} K, "
            f"exact ball mean {fmt(cf)} K, rel {r:.4e}, V/V_nom {float(vols[k]) / float(noms[k]):.6f}; "
            f"theta_h at the wall point {fmt(theta_pt[j])} K")
    max_inset = max(r[2] for r in inset_rows)
    say(f"[inset] max rel {max_inset:.4e} (predicted <= 5e-3, range [1e-4, 1e-2]; in range: "
        f"{1e-4 <= max_inset <= 1e-2}; printed only, not asserted)")

    # (c-half)
    k = 6 + len(inset_centres)
    c_half = abs(float(vols[k]) / float(noms[k]) - 1.0)
    say(f"[control] (c-half) ball centred on P1 (+y wall): |V/V_nom - 1| {c_half:.4f} > "
        f"{INTERIOR_BAND:.0e} (predicted 0.5, range [0.45, 0.55]; in range: "
        f"{0.45 <= c_half <= 0.55}; ceiling 1); its mean {fmt(float(vals[k]))} K beside "
        f"theta_h(P1) {fmt(theta_pt[1])} K")

    # (c-mv)
    point0 = point_c[0]
    say(f"[control] (c-mv) point series at c0 {fmt(point0)} K vs ball mean {fmt(point0 - mv)} K: "
        f"Q a^2/(10K) = {mv:.4e} K ({mv / point0:.3e} relative; record {REC['cmv']:.4e} K, "
        f"{S6W2}:71; printed only)")

    # (c-wall)
    res2 = Study(geo, None, UniformHeatSource(Q)).probe(
        BOX_PROBES[:1], RobinCooling(K_GEL, 2.0 * H_GEL, T_INF))
    clk.lap("Study.probe h doubled")
    c_wall = abs(float(res2.values[0]) / (point0 - mv) - 1.0)
    say(f"[control] (c-wall) h doubled at c0: {c_wall:.4e} > {INTERIOR_BAND:.0e} (backed, record "
        f"{REC['s6_cwall']:.4e}, {S6W2}:72; step 5 {S5_C_WALL:.4e}, {S5W1}:70; ceiling 1)")

    say(f"[print] max|Im theta| fraction in probe {pres.max_imag_ratio:.1e}")

    # ---- assertions ------------------------------------------------------------------------
    assert cells == s5.STEP1_M4_CELLS, f"(S1-m4) cells {cells}"
    assert p1_dofs == s5.STEP1_M4_P1_DOFS, f"(S1-m4) dofs {p1_dofs}"
    assert f"{whole:.4e}" == s5.STEP1_M4_WHOLE, f"(S1-m4) whole {whole:.4e}"
    assert interior <= INTERIOR_BAND, f"(S1-m4) interior {interior:.4e}"
    assert gauss <= S2_GAUSS_BAND, f"(S1-m4) Gauss {gauss:.3e}"
    assert d_rep <= L0_BAND, f"(I-rep) {d_rep:.3e}"
    assert a_val <= B_BAND, f"(A-box) value {a_val:.3e}"
    assert a_vol <= B_BAND, f"(A-box) volume {a_vol:.3e}"
    for j, r in enumerate(r_rel):
        assert r <= INTERIOR_BAND, f"(R-box) c{j} {r:.4e}"
    assert bool(np.all(ok_pt)), f"(V) ok {ok_pt}"
    assert bool(np.all(owner_orig >= 0)), f"(V) owners {owner_orig}"
    assert bool(np.all(np.isfinite(theta_pt))), "(V) finite"
    assert c_wall > INTERIOR_BAND, f"(c-wall) did not miss: {c_wall:.4e}"
    assert c_wall < 1.0, f"(c-wall) above ceiling: {c_wall:.4e}"
    assert c_half > INTERIOR_BAND, f"(c-half) did not miss: {c_half:.4e}"
    assert c_half < 1.0, f"(c-half) above ceiling: {c_half:.4e}"

    # ---- outputs: CSV, then XDMF, then adopt_host_ownership ---------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    if rank0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_probe_table.csv"
        with open(csv_path, "w") as fh:
            fh.write("label,x_m,y_m,z_m,theta_h_K,theta_series_K,e_p,cx_m,cy_m,cz_m,"
                     "ball_mean_K,exact_ball_mean_K,ball_rel,volume_ratio\n")
            for i in range(len(PROBE_PTS)):
                p, cf, r, vr = inset_rows[i]
                fh.write(f"{PROBE_LABELS[i]},{PROBE_PTS[i][0]:.9e},{PROBE_PTS[i][1]:.9e},"
                         f"{PROBE_PTS[i][2]:.9e},{theta_pt[i]:.5e},{ser_pt[i]:.5e},{e_p[i]:.5e},"
                         f"{inset_centres[i][0]:.9e},{inset_centres[i][1]:.9e},"
                         f"{inset_centres[i][2]:.9e},{p:.5e},{cf:.5e},{r:.5e},{vr:.6f}\n")
    comm.Barrier()

    V1 = fem.functionspace(msh, ("Lagrange", 1))

    def p1(name, src, offset=0.0):
        f = fem.Function(V1, name=name)
        f.interpolate(src)
        f.x.array[:] = np.real(f.x.array) - offset
        f.x.scatter_forward()
        return f

    f_th = p1("theta_study_steady", res.temperature, T_INF)
    f_se = p1("theta_series_steady", ex)
    f_d = fem.Function(V1, name="theta_minus_series")
    f_d.x.array[:] = np.real(f_th.x.array) - np.real(f_se.x.array)
    f_d.x.scatter_forward()
    f_w = fem.Function(V1, name="theta_wall_2h")
    f_w.x.array[:] = np.real(_on(msh, res2.thermal.temperature).x.array) - T_INF
    f_w.x.scatter_forward()
    f_b = fem.Function(V1, name="probe_balls")
    xd = np.real(V1.tabulate_dof_coordinates())
    ind = np.zeros(xd.shape[0])
    for j, c in enumerate(inset_centres):
        ind[np.sum((xd - np.asarray(c)) ** 2, axis=1) <= a * a] = j + 1.0
    f_b.x.array[:] = ind
    f_b.x.scatter_forward()
    ft.name = "robin_facets"
    fields = {f.name: f for f in (f_th, f_se, f_d, f_w, f_b)}
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, fields, comm=comm, facet_tags=ft)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    clk.lap("outputs")
    if rank0:
        print(f"\n[output] {xdmf_path.name}: CG1 theta_study_steady, theta_series_steady, "
              "theta_minus_series, theta_wall_2h (h doubled), probe_balls (j + 1 inside inset ball j), "
              "Robin facets")
        print(f"[output] {BASENAME}_probe_table.csv ({len(PROBE_PTS)} rows)")
        print(f"[time] total {time.perf_counter() - clk.t0:.1f} s")
        print("[scope] steady uniform box at m = 4, P1, through the Study; ball probes and point "
              "values at the six F2182 positions; -n 2, complex build; no transient, no EM, no "
              "coil, no h-rate, no compliance, no external check claimed")
        print("All assertions hold")


if __name__ == "__main__":
    main()
