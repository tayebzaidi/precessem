"""Example (`EX-78`): where the tuned birdcage's 1/Q goes -- Q from the field's
stored energy against the Y_in slope, loaded and unloaded, and the phantom's
share.

`WF-5` steps 1-2 gate the tuned F-small copper birdcage on `ANS-7`'s Route A
hole mesh (80 181 cells, degree 1, tracked kappa, f0 ~ 63.98 MHz) two ways:
Q from the susceptance slope of the input admittance (Foster), and Q = w0 W / P
from one solve's stored energy over its loss. This example runs the Cu column
loaded (phantom sigma = `th14.SALINE_SIGMA`) and unloaded (sigma = 0, eps_r
kept), and shows the loss partition on the mesh: |E| loaded vs unloaded, the
phantom loss density, and the copper loss per adjacent cell.

**It asserts, it does not merely render.** Every record, band and helper is
*imported* from the modules that gate the capability
(``tests/validation/test_wf5_q_two_ways.py``, ``test_wf5_loading.py``,
``test_birdcage_power_identity.py``), never restated. The fixture is composed
from imported helpers exactly as ``test_wf5_loading.py:74-93`` does; the fields
are captured by wrapping ``ans7._route_a_solve`` around each ``_budget`` call
(the gate's own pattern, ``test_wf5_q_two_ways.py:186-206``), restored in
``finally`` on every rank.

Asserted: global cells == tracked ``n_cells``; (A0) loaded secant f0 / Q reproduce
``rows_route_a.cu`` <= ``RECORD_RTOL``; loaded |Q_field - Q_slope|/Q_slope <=
``B_Q_FROZEN`` and the loaded identity <= ``DISCRETE_IDENTITY_RTOL``; (A1)
P_phantom,U == 0.0; (A2) the unloaded identity; (A3) unloaded Q two ways <=
``B_Q_U_FROZEN``; (A4) |Delta(1/Q)_slope - Pi|/Pi <= B_Delta; the map sums
reproduce P_phantom,L and P_surf (L and U) <= ``DISCRETE_IDENTITY_RTOL`` (DG0
partition of unity: a consistency check, not physics). Negative control
(asserted only as "> band"): the budget with W_sheets dropped misses all three
Q bands. Printed, not asserted: energy shares, the 1/Q partition, Pi_ph beside
``STEP1_PI_PH``, Delta f0/f0 (second order).

Scope: Cu column, default hole mesh, degree 1, tracked kappa, ``-n 4``. No PEC
column, no phantom resolution, no absolute Q, no tuning, no AED claim.

Run it through the example runner (the ``mri:`` group sources the complex
build automatically)::

    ./scripts/run_examples.sh -e mri:5 -n 4 -t 300

Output lands in ``examples/mri/paraview_output``: open
``mri_05_birdcage_q_loaded_unloaded_combined.xdmf``.
"""

from __future__ import annotations

import contextlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

from fem_em_solver.core.cavity import surface_resistance_ohm
from fem_em_solver.io.paraview_utils import (
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure

# The gated fixtures live in tests; the runner puts only ``src`` on PYTHONPATH.
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.validation import test_th14_birdcage_copper as th14  # noqa: E402
from tests.validation.test_birdcage_power_identity import (  # noqa: E402
    DISCRETE_IDENTITY_RTOL,
)
from tests.validation.test_wf5_loading import (  # noqa: E402
    B_Q_U_FROZEN,
    STEP1_PI_PH,
    UNLOADED_MHZ,
    _rows_on_grid,
    _sweep,
    _unloaded,
    _w_tot,
)
from tests.validation.test_wf5_q_two_ways import (  # noqa: E402
    ANS7_DIR,
    ANS7_SCRIPT,
    B_Q_FROZEN,
    CUBIC_MHZ,
    RECORD_RTOL,
    _budget,
    _load,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "mri_05_birdcage_q_loaded_unloaded"
PORTS18_SCRIPT = _REPO_ROOT / "examples" / "ports" / "18_birdcage_copper_leontovich.py"


@contextlib.contextmanager
def _capture_fields(ans7, store):
    """Stash the ``fields`` ``_route_a_solve`` returns (called with
    ``return_fields=True`` by ``_budget``); restore the original on every rank."""
    original = ans7._route_a_solve

    def wrapper(*args, **kwargs):
        out = original(*args, **kwargs)
        if kwargs.get("return_fields"):
            store["fields"] = out[1]
        return out

    ans7._route_a_solve = wrapper
    try:
        yield
    finally:
        ans7._route_a_solve = original


def _rss_gib(comm) -> float:
    """Summed peak RSS over ranks, GiB (collective)."""
    kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return float(comm.allreduce(kib, op=MPI.SUM)) / 1024.0**2


def _e_magnitude_cg1(msh, e, name):
    v = fem.functionspace(msh, ("Lagrange", 1))
    f = fem.Function(v, name=name)
    f.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), v.element.interpolation_points)
    )
    f.x.scatter_forward()
    return f


def _phantom_loss_density(msh, fields, comm):
    """Cell-averaged phantom loss density 0.5 sigma |E|^2, W/m^3 (DG0), and the
    sum of (density x cell volume) over the mesh (allreduced), W."""
    from dolfinx import la

    dg0 = fem.functionspace(msh, ("DG", 0))
    v = ufl.TestFunction(dg0)
    dx = ufl.Measure("dx", domain=msh)
    e = fields.e_complex
    # complex mode: the real DG0 test function enters conjugated.
    b = fem.assemble_vector(
        fem.form(0.5 * ufl.real(fields.sigma_field) * ufl.inner(e, e) * ufl.conj(v) * dx)
    )
    b.scatter_reverse(la.InsertMode.add)
    b.scatter_forward()
    vol = fem.assemble_vector(fem.form(ufl.conj(v) * dx))
    vol.scatter_reverse(la.InsertMode.add)
    vol.scatter_forward()
    n_loc = dg0.dofmap.index_map.size_local
    dens = fem.Function(dg0, name="phantom_loss_density_W_per_m3")
    dens.x.array[:] = np.real(b.array) / np.real(vol.array)
    dens.x.scatter_forward()
    total = float(comm.allreduce(np.sum(dens.x.array[:n_loc] * np.real(vol.array[:n_loc])),
                                 op=MPI.SUM))
    return dens, total


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mri:` "
            "group sources it automatically"
        )

    ans7 = _load(ANS7_SCRIPT, "ans7_driven_sweep_wf5")
    ports18 = _load(PORTS18_SCRIPT, "ports18_copper_leontovich_ex78")

    # ---- fixture: test_wf5_loading.py:74-93, tracked kappa ---------------------
    tracked_all = json.loads((ANS7_DIR / "metrics.json").read_text())
    kappa = float(tracked_all["kappa"])
    n_cells_tracked = int(tracked_all["n_cells"])
    tracked = tracked_all["rows_route_a"]["cu"]
    msh, cell_tags, facet_tags, _diag, t_mesh = th14._sheets_build(True, as_hole=True)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    n_cells = int(msh.topology.index_map(tdim).size_global)
    outer_tags, _census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_k, port_defs_k, specs_k = th14._build_ports(
        msh, cell_tags, comm, width_correction_kappa=kappa
    )
    geo_k = (msh, cell_tags, facet_tags, outer_tags, tags_k, port_defs_k)
    fb = {"th14": th14, "geo": geo_k, "specs": specs_k, "tracked": tracked_all, "comm": comm}

    region_names = {2: "air", th14.PHANTOM_CELL_TAG: "phantom"}
    for spec in specs_k:
        region_names[int(spec.facet_tag)] = f"port {spec.port_id}"
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=(2,),
        translucent_tags=(th14.PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="mri:5 -- tuned copper birdcage F-small hole mesh, phantom, four gap ports",
        comm=comm,
    )
    comm.Barrier()
    if rank0:
        print("=" * 78)
        print("EX-78 -- where the tuned birdcage's 1/Q goes: Q two ways, loaded / unloaded")
        print("=" * 78)
        print(
            f"\n[fixture] hole mesh {n_cells} cells (tracked {n_cells_tracked}), mesh "
            f"{t_mesh:.1f} s; kappa (tracked) = {kappa:.9e}; SALINE_SIGMA (loaded) = "
            f"{th14.SALINE_SIGMA}; degree 1, Cu column",
            flush=True,
        )
    assert n_cells == n_cells_tracked, (
        f"cell count {n_cells} != tracked {n_cells_tracked}: the fixture is not the record's"
    )

    # ---- Y_in sweeps (collective) ---------------------------------------------
    t_sw = time.perf_counter()
    y_l = _sweep(ans7, fb, CUBIC_MHZ, "cu")
    with _unloaded(th14):
        y_u = _sweep(ans7, fb, UNLOADED_MHZ, "cu")
    rl = _rows_on_grid(CUBIC_MHZ, y_l)
    ru = _rows_on_grid(UNLOADED_MHZ, y_u)
    sweep_seconds = time.perf_counter() - t_sw

    # ---- budgets with captured fields (collective) ----------------------------
    t_bud = time.perf_counter()
    store_l, store_u = {}, {}
    with _capture_fields(ans7, store_l):
        bl = _budget(ans7, fb, rl["f0_hz"], "cu")
    with _unloaded(th14):
        with _capture_fields(ans7, store_u):
            bu = _budget(ans7, fb, ru["f0_hz"], "cu")
    budget_seconds = time.perf_counter() - t_bud

    # ---- composition (cited: test_wf5_loading.py:209-230, q_two_ways.py:283-292)
    q_l_slope = 1.0 / rl["inv_q_sec"]
    q_u_slope = 1.0 / ru["inv_q_sec"]
    df0_rec = abs(rl["f0_hz"] - tracked["f0_hz"]) / tracked["f0_hz"]
    dq_rec = abs(q_l_slope - tracked["q_loaded"]) / tracked["q_loaded"]
    wl, wu = _w_tot(bl), _w_tot(bu)
    pl = bl["p_phantom"] + bl["p_surf"]
    pu = bu["p_phantom"] + bu["p_surf"]
    om_l, om_u = bl["omega"], bu["omega"]
    q_field_l = om_l * wl / pl
    q_field_u = om_u * wu / pu
    inv_qf_l, inv_qf_u = 1.0 / q_field_l, 1.0 / q_field_u
    miss_l = abs(q_field_l - q_l_slope) / q_l_slope
    id_l = abs(pl - (bl["p_src"] - bl["p_sheet1"])) / pl
    id_u = abs(pu - (bu["p_src"] - bu["p_sheet1"])) / bu["p_src"]
    miss_u = abs(q_field_u - q_u_slope) / q_u_slope
    pi = inv_qf_l - inv_qf_u
    d_slope = rl["inv_q_sec"] - ru["inv_q_sec"]
    b_delta = (B_Q_FROZEN * inv_qf_l + B_Q_U_FROZEN * inv_qf_u) / pi
    miss_d = abs(d_slope - pi) / pi
    pi_ph = bl["p_phantom"] / (om_l * wl)
    pi_cu = bl["p_surf"] / (om_l * wl)
    df0 = (ru["f0_hz"] - rl["f0_hz"]) / rl["f0_hz"]
    ws_l = sum(bl["w_sheet"].values())
    ws_u = sum(bu["w_sheet"].values())
    # negative control: W_sheets dropped
    q_ns_l = om_l * (wl - ws_l) / pl
    q_ns_u = om_u * (wu - ws_u) / pu
    miss_ns_l = abs(q_ns_l - q_l_slope) / q_l_slope
    miss_ns_u = abs(q_ns_u - q_u_slope) / q_u_slope
    pi_ns = 1.0 / q_ns_l - 1.0 / q_ns_u
    miss_ns_d = abs(d_slope - pi_ns) / pi_ns

    # ---- maps (collective) ------------------------------------------------------
    e_l, e_u = store_l["fields"].e_complex, store_u["fields"].e_complex
    z_s_l = (1.0 + 1.0j) * surface_resistance_ohm(om_l, th14.COPPER_SIGMA)
    z_s_u = (1.0 + 1.0j) * surface_resistance_ohm(om_u, th14.COPPER_SIGMA)
    cu_l, cu_l_total = ports18._coil_surface_loss_per_cell_w(msh, facet_tags, z_s_l, e_l, comm)
    cu_u, cu_u_total = ports18._coil_surface_loss_per_cell_w(msh, facet_tags, z_s_u, e_u, comm)
    cu_l.name = "coil_surface_loss_per_cell_W_loaded"
    cu_u.name = "coil_surface_loss_per_cell_W_unloaded"
    dens_l, ph_total = _phantom_loss_density(msh, store_l["fields"], comm)
    e_mag_l = _e_magnitude_cg1(msh, e_l, "E_magnitude_loaded")
    e_mag_u = _e_magnitude_cg1(msh, e_u, "E_magnitude_unloaded")
    map_ph = abs(ph_total - bl["p_phantom"]) / bl["p_phantom"]
    map_cu_l = abs(cu_l_total - bl["p_surf"]) / bl["p_surf"]
    map_cu_u = abs(cu_u_total - bu["p_surf"]) / bu["p_surf"]
    rss = _rss_gib(comm)

    if rank0:
        print(
            f"\n[sweeps] loaded {CUBIC_MHZ} + unloaded {UNLOADED_MHZ} MHz, Cu: "
            f"{sweep_seconds:.1f} s; budgets {budget_seconds:.1f} s; RSS(sum) {rss:.2f} GiB",
            flush=True,
        )
        print(
            f"[A0] loaded secant f0 {rl['f0_hz']:.6f} Hz (record {tracked['f0_hz']:.6f}, rel "
            f"{df0_rec:.3e}); Q_slope,L {q_l_slope:.10f} (record {tracked['q_loaded']:.10f}, "
            f"rel {dq_rec:.3e}; ASSERTED <= {RECORD_RTOL:.0e})",
            flush=True,
        )
        print(
            f"[energy] loaded W_tot {wl:.9e} J: W_m {bl['w_m'] / wl:.6f}, W_e {bl['w_e'] / wl:.6f}, "
            f"W_sheets {ws_l / wl:.6f}, W_C1 {bl['w_c1'] / wl:.6f}, W_surf {bl['w_surf'] / wl:.4e}",
            flush=True,
        )
        print(
            f"[energy] unloaded W_tot {wu:.9e} J: W_m {bu['w_m'] / wu:.6f}, W_e "
            f"{bu['w_e'] / wu:.6f}, W_sheets {ws_u / wu:.6f}, W_C1 {bu['w_c1'] / wu:.6f}, "
            f"W_surf {bu['w_surf'] / wu:.4e}",
            flush=True,
        )
        print(
            f"[Q_L] Q_field,L {q_field_l:.6f} vs Q_slope,L {q_l_slope:.6f}: miss {miss_l:.4e} "
            f"(ASSERTED <= B_Q {B_Q_FROZEN:.0e}; WF-5 record 920.010244 vs 919.037377, "
            f"1.0586e-3); loaded identity {id_l:.3e} (ASSERTED <= {DISCRETE_IDENTITY_RTOL:.0e})",
            flush=True,
        )
        print(
            f"[A1] P_phantom,U = {bu['p_phantom']!r} W (ASSERTED == 0.0); [A2] unloaded "
            f"identity {id_u:.3e} (ASSERTED <= {DISCRETE_IDENTITY_RTOL:.0e})",
            flush=True,
        )
        print(
            f"[A3] unloaded f0 {ru['f0_hz']:.6f} Hz; Q_field,U {q_field_u:.6f} vs Q_slope,U "
            f"{q_u_slope:.6f}: miss {miss_u:.4e} (ASSERTED <= B_Q,U {B_Q_U_FROZEN:.0e}; WF-5 "
            f"record 2283.280107 vs 2279.420027, 1.6934e-3)",
            flush=True,
        )
        print(
            f"[A4] Delta(1/Q)_slope {d_slope:.9e}, Pi {pi:.9e}: |Delta - Pi|/Pi {miss_d:.4e} "
            f"(ASSERTED <= B_Delta {b_delta:.4e}; WF-5 record 6.3012e-4 <= 6.3743e-3)",
            flush=True,
        )
        tot = inv_qf_l
        print(
            f"[1/Q] loaded 1/Q_field {inv_qf_l:.9e} (record 1.086944418e-3): phantom {pi_ph:.6e} "
            f"(record 6.489692e-4; STEP1_PI_PH {STEP1_PI_PH['cu']:.6e}) = {pi_ph / tot:.1%}, "
            f"copper {pi_cu:.6e} (record 4.379752e-4) = {pi_cu / tot:.1%}; unloaded 1/Q_field "
            f"{inv_qf_u:.9e} (all copper)",
            flush=True,
        )
        print(
            f"[loading] Delta f0/f0 = {df0:+.6e} (printed; second order, WF-5 record "
            f"-4.523479e-6); |Delta - Pi_ph|/Pi_ph = {abs(d_slope - pi_ph) / pi_ph:.4e} (printed)",
            flush=True,
        )
        print(
            f"[maps] sum(phantom density x vol) {ph_total:.9e} W vs P_phantom,L "
            f"{bl['p_phantom']:.9e} W: rel {map_ph:.3e}; sum(copper per cell) L {cu_l_total:.9e} "
            f"vs P_surf,L {bl['p_surf']:.9e}: rel {map_cu_l:.3e}; U {cu_u_total:.9e} vs "
            f"{bu['p_surf']:.9e}: rel {map_cu_u:.3e} (ASSERTED <= {DISCRETE_IDENTITY_RTOL:.0e}; "
            f"DG0 partition of unity, consistency not physics)",
            flush=True,
        )
        print(
            f"[control] W_sheets omitted: loaded miss {miss_ns_l:.4e} (ASSERTED > "
            f"{B_Q_FROZEN:.0e}; WF-5 record 0.24536); unloaded miss {miss_ns_u:.4e} (ASSERTED > "
            f"{B_Q_U_FROZEN:.0e}; PREDICTED 0.2449 in [0.240, 0.250]); Pi' vs Delta(1/Q)_slope "
            f"{miss_ns_d:.4e} (ASSERTED > B_Delta {b_delta:.4e}; PREDICTED 0.2457 in "
            f"[0.240, 0.250])",
            flush=True,
        )

    # ---- ParaView ---------------------------------------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cell_tags,
        {
            "E_magnitude_loaded": e_mag_l,
            "E_magnitude_unloaded": e_mag_u,
            "phantom_loss_density_W_per_m3": dens_l,
            "coil_surface_loss_per_cell_W_loaded": cu_l,
            "coil_surface_loss_per_cell_W_unloaded": cu_u,
        },
        comm=comm,
        facet_tags=facet_tags,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    if rank0:
        print(f"\n[paraview] wrote {xdmf_path}", flush=True)

    fails = []
    if not (df0_rec <= RECORD_RTOL and dq_rec <= RECORD_RTOL):
        fails.append(("A0", df0_rec, dq_rec))
    if not miss_l <= B_Q_FROZEN:
        fails.append(("Q_L", miss_l))
    if not id_l <= DISCRETE_IDENTITY_RTOL:
        fails.append(("loaded identity", id_l))
    if not bu["p_phantom"] == 0.0:
        fails.append(("A1", bu["p_phantom"]))
    if not id_u <= DISCRETE_IDENTITY_RTOL:
        fails.append(("A2", id_u))
    if not miss_u <= B_Q_U_FROZEN:
        fails.append(("A3", miss_u))
    if not miss_d <= b_delta:
        fails.append(("A4", miss_d, b_delta))
    for name, val in (("map phantom", map_ph), ("map copper L", map_cu_l),
                      ("map copper U", map_cu_u)):
        if not val <= DISCRETE_IDENTITY_RTOL:
            fails.append((name, val))
    if not miss_ns_l > B_Q_FROZEN:
        fails.append(("control loaded inside band", miss_ns_l))
    if not miss_ns_u > B_Q_U_FROZEN:
        fails.append(("control unloaded inside band", miss_ns_u))
    if not miss_ns_d > b_delta:
        fails.append(("control Pi' inside band", miss_ns_d))
    assert not fails, f"EX-78 gates: {fails}"

    if rank0:
        print(f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.",
              flush=True)


if __name__ == "__main__":
    main()
