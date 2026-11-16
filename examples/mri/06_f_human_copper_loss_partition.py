"""Example (`EX-80`): the F-human 32x32 on copper and its loss partition.

`WF-12` step 2 (`tests/validation/test_wf12_f_human_copper_matrix.py`) drives
all 32 ring-gap ports of the F-human hole mesh (369 387 cells, degree 1,
64 MHz, 50 ohm) with `TH-14`'s Leontovich copper wall and one PEC column.
``_run_step2`` is that gate's own body lifted to module level (rule (a) lift,
EX-79's precedent) with an optional ``field_hook``; nothing is restated here.

The example re-executes the gate through ``_run_step2(field_hook=...)``, asserts
the printed records beside ``20261002T094751Z_WF-12-step2-W2.log:6352-6414``,
and, on the P17 drive (Cu and PEC), builds per-cell maps inside the hook:
|E|, the phantom loss density and the copper loss per adjacent cell.  Their
owned sums are asserted against the gate's P_phantom / P_surf (example-side
identities, true by construction: a DG0 partition of unity).

Negative control (asserted): on every Cu drive the energy identity with P_surf
dropped misses by >= BLIND_FACTOR x DISCRETE_IDENTITY_RTOL; on P17 also with
the map-sum P_phantom.

Scope: step 2's fixture, degree 1, 64 MHz, 50 ohm, one PEC column.  Not
|B1+| or a hotspot, not C_mode, not a PEC 32x32, not degree 2, not loss
against AED, not absolute loss.

Run through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:6 -n 8 -t 480

Output: ``examples/mri/paraview_output/mri_06_f_human_copper_loss_partition_combined.xdmf``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem, la
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.setup_figure import write_setup_figure

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.validation import test_wf7_f_human_b1_quadrature as wf7  # noqa: E402
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_birdcage_ring_column import (  # noqa: E402
    AZIMUTH_MATCH_DEG,
    OPPOSITE_SPREAD_BAND,
)
from tests.validation.test_port_birdcage_ring_matrix import (  # noqa: E402
    MATRIX_PASSIVITY_CEILING,
    N_SYMMETRY_CLASSES,
    REUSE_REPRODUCTION_RTOL,
)
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND  # noqa: E402
from tests.validation.test_th14_birdcage_copper import (  # noqa: E402
    BIRDCAGE_CONDUCTOR_SURFACE_TAG,
)
from tests.validation.test_wf12_f_human_copper_matrix import (  # noqa: E402
    BLIND_FACTOR,
    DRIVEN,
    HOLE_CELL_RECORD,
    P17_P_SRC_RECORD_W,
    P17_P_SURF_RECORD_W,
    PEC_SURF_RATIO_CEILING,
    PREDICTED_COIL_SHARE,
    SIGMA800_COIL_SHARE,
    _run_step2,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "mri_06_f_human_copper_loss_partition"
W2 = "20261002T094751Z_WF-12-step2-W2.log"


def _dg0_vec(form, dg0):
    b = fem.assemble_vector(fem.form(form))
    b.scatter_reverse(la.InsertMode.add)
    b.scatter_forward()
    return np.real(b.array)


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mri:` "
            "group sources it automatically"
        )

    if rank0:
        print("=" * 72)
        print("EX-80 -- F-human 32x32 on copper and its loss partition (WF-12 step 2, re-executed)")
        print("=" * 72)
        print(f"[records] every anchor below is beside {W2}:6352-6414", flush=True)

    store = {"maps": {}, "sums": {}}

    # Stash the gate's ``built`` (facet and cell tags) as the gate builds it: the
    # gate imports ``_build_f_human`` inside the call, so this wrapper is seen.
    original_build = wf7._build_f_human

    def _capturing_build(*args, **kwargs):
        built = original_build(*args, **kwargs)
        store["built"] = built
        return built

    def hook(kind, pid, e, sigma_f, z_s):
        if pid != DRIVEN:
            return
        built = store["built"]
        msh = e.function_space.mesh
        cell_tags = built["ctx"]["cell_tags"]
        facets = built["facets"]
        dg0 = fem.functionspace(msh, ("DG", 0))
        v = ufl.TestFunction(dg0)
        n_loc = dg0.dofmap.index_map.size_local
        dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
        vol = _dg0_vec(ufl.conj(v) * ufl.dx(domain=msh), dg0)

        def named(name, values):
            f = fem.Function(dg0, name=name)
            f.x.array[:] = values
            f.x.scatter_forward()
            return f

        e_sq = fem.Function(dg0, name=f"e_sq_{kind}")
        e_sq.interpolate(fem.Expression(ufl.inner(e, e), dg0.element.interpolation_points))
        mag = np.sqrt(np.maximum(np.real(e_sq.x.array), 0.0))
        loss = _dg0_vec(0.5 * sigma_f * ufl.inner(e, e) * ufl.conj(v) * dx(PHANTOM_CELL_TAG), dg0)
        maps = {
            f"E_magnitude_{kind}_{pid}": named(f"E_magnitude_{kind}_{pid}", mag),
            f"phantom_loss_density_{kind}_{pid}": named(
                f"phantom_loss_density_{kind}_{pid}", loss / vol),
        }
        store["sums"][f"phantom_{kind}"] = float(comm.allreduce(np.sum(loss[:n_loc]), op=MPI.SUM))
        if kind == "cu":
            ds_c = ufl.Measure("ds", domain=msh, subdomain_data=facets,
                               subdomain_id=BIRDCAGE_CONDUCTOR_SURFACE_TAG)
            nrm = ufl.FacetNormal(msh)
            re_y = fem.Constant(msh, default_scalar_type(0.5 * np.real(1.0 / complex(z_s))))
            surf = _dg0_vec(re_y * ufl.inner(ufl.cross(nrm, e), ufl.cross(nrm, e))
                            * ufl.conj(v) * ds_c, dg0)
            name = f"coil_surface_loss_per_cell_W_cu_{pid}"
            maps[name] = named(name, surf)
            store["sums"]["surf_cu"] = float(comm.allreduce(np.sum(surf[:n_loc]), op=MPI.SUM))
        store["maps"].update(maps)

    wf7._build_f_human = _capturing_build
    try:
        r = _run_step2(field_hook=hook)
    finally:
        wf7._build_f_human = original_build

    msh, cell_tags = r["msh"], r["cell_tags"]
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            CONDUCTOR_CELL_TAG: "copper conductor",
            2: "air",
            PHANTOM_CELL_TAG: "phantom",
            **{PORT_LOWER + int(s.port_id[1:]): f"port {s.port_id}" for s in r["specs"]},
            **{PORT_UPPER + int(s.port_id[1:]): f"port {s.port_id}" for s in r["specs"]},
        },
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="mri:6 -- F-human birdcage on copper, 32 ring-gap ports, phantom",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- anchors, asserted at the imported bands ----------------------------
    cols, m, pec = r["cols"], r["matrix"], r["pec"]
    first = cols[DRIVEN]
    cell_rel = r["cells"] / HOLE_CELL_RECORD - 1.0
    col_min, col_max = float(m["col_norms"].min()), float(m["col_norms"].max())
    worst_spread = m["spreads"][m["worst"]][0]
    sums = store["sums"]
    p_src, p_sh = first["p_src"], first["p_sheets"]
    map_cu_ph = abs(sums["phantom_cu"] / first["p_phantom"] - 1.0)
    map_cu_surf = abs(sums["surf_cu"] / first["p_surf"] - 1.0)
    map_pec_ph = abs(sums["phantom_pec"] / pec["d"]["p_phantom"] - 1.0)
    blind_map = abs(p_src - p_sh - sums["phantom_cu"] - first["p_rest"]) / abs(p_src)
    if rank0:
        print(f"\n[anchor] (A0) cells {r['cells']} rel {cell_rel:.3e} <= {CELL_COUNT_BAND}  (W2: 369387, 0)")
        print(f"[anchor] (R) P_src rel {r['r_src']:.3e} and P_surf rel {r['r_surf']:.3e} "
              f"<= {REUSE_REPRODUCTION_RTOL:g}  (W2: 2.111e-11, 1.507e-11; records "
              f"{P17_P_SRC_RECORD_W:.9e}, {P17_P_SURF_RECORD_W:.9e})")
        print(f"[anchor] (P) worst rel {r['worst_rel']:.3e} <= {DISCRETE_IDENTITY_RTOL:g} over "
              f"{len(cols)} drives, {r['kinds'].count('solve')} factorising + "
              f"{r['kinds'].count('held')} held  (W2: 5.117e-13)")
        print(f"[anchor] (ii) reciprocity {m['ratio']:.6e} <= {RECIPROCITY_BAND:g}  (W2: round-off "
              f"level, 9.166882e-15)")
        print(f"[anchor] (iii) sigma_max {m['sigma_max']:.12f} <= {MATRIX_PASSIVITY_CEILING}; "
              f"column sum|S|^2 {col_min:.6f} .. {col_max:.6f}  (W2: 0.999993339099; 0.993097 .. 0.993680)")
        print(f"[anchor] (iv') classes {m['n_classes']} (expected {N_SYMMETRY_CLASSES}); residual "
              f"{m['worst_residual_deg']:.3e} deg < {AZIMUTH_MATCH_DEG:g}; worst spread "
              f"{worst_spread:.4e} <= {OPPOSITE_SPREAD_BAND}  (W2: 18, 5.684e-14, 3.2748e-02)")
        print(f"[anchor] (PEC) identity w/o P_surf rel {pec['rel_pec']:.3e} <= {DISCRETE_IDENTITY_RTOL:g}; "
              f"P_surf(PEC) = {pec['surf_ratio']:.3e} x Cu <= {PEC_SURF_RATIO_CEILING:g}  "
              f"(W2: round-off level, 2.78e-23)")
        print("[classes] same/other ring, steps, n, mean|S|, spread %")
        for (same, steps), (sp, cnt, mean) in sorted(m["spreads"].items(),
                                                     key=lambda kv: (not kv[0][0], kv[0][1])):
            print(f"   {'same ' if same else 'other'} {steps:2d}  n={cnt:4d}  mean {mean:.9e}  "
                  f"spread {sp * 100:8.4f}%")
        print(f"[map] Cu phantom map sum rel {map_cu_ph:.3e}, Cu coil-loss map sum rel "
              f"{map_cu_surf:.3e}, PEC phantom map sum rel {map_pec_ph:.3e} "
              f"<= {DISCRETE_IDENTITY_RTOL:g}")
        print(f"[control] P_surf-blind: min over {len(cols)} Cu drives {r['min_blind'] / DISCRETE_IDENTITY_RTOL:.1f}x "
              f"the band (bar {BLIND_FACTOR:.0f}x; W2 123.3x); P17 with map P_phantom "
              f"{blind_map / DISCRETE_IDENTITY_RTOL:.1f}x")
        print("\n[printed] P17 partition, fractions of P_src: "
              f"sheets {p_sh / p_src:.6f}, phantom {first['p_phantom'] / p_src:.6e}, "
              f"rest {first['p_rest'] / p_src:.6e}, copper wall {first['p_surf'] / p_src:.6e}")
        print(f"[printed] coil share P_surf/(P_src - sum P_sheet) min {min(r['shares']):.4f} max "
              f"{max(r['shares']):.4f} (record 0.3699-0.3737; PREDICTED {PREDICTED_COIL_SHARE}); "
              f"sigma = 800 {SIGMA800_COIL_SHARE} -- never a share of P_src")
        print(f"[printed] P_phantom Cu {first['p_phantom']:.9e} W, PEC {pec['d']['p_phantom']:.9e} W")
        print(f"[printed] Cu != PEC: delta sum|S_i,17|^2 {pec['dsum']:.4e}, max|dS| "
              f"{pec['max_ds']:.4e}  (W2: 3.3903e-04, 2.2094e-04)", flush=True)

    assert abs(cell_rel) <= CELL_COUNT_BAND and r["cells"] == HOLE_CELL_RECORD, r["cells"]
    assert r["r_src"] <= REUSE_REPRODUCTION_RTOL and r["r_surf"] <= REUSE_REPRODUCTION_RTOL
    assert r["worst_rel"] <= DISCRETE_IDENTITY_RTOL, r["worst_rel"]
    assert r["kinds"].count("solve") == 1 and r["kinds"].count("held") == 31, r["kinds"]
    assert m["ratio"] <= RECIPROCITY_BAND, m["ratio"]
    assert m["sigma_max"] <= MATRIX_PASSIVITY_CEILING and col_max <= MATRIX_PASSIVITY_CEILING
    assert m["n_classes"] == N_SYMMETRY_CLASSES, m["n_classes"]
    assert m["worst_residual_deg"] < AZIMUTH_MATCH_DEG, m["worst_residual_deg"]
    assert worst_spread <= OPPOSITE_SPREAD_BAND, worst_spread
    assert pec["rel_pec"] <= DISCRETE_IDENTITY_RTOL and pec["surf_ratio"] <= PEC_SURF_RATIO_CEILING
    # reproduction of the W2 digits (deterministic quantities only)
    assert abs(m["sigma_max"] - 0.999993339099) <= 1e-9, m["sigma_max"]
    assert round(col_min, 6) == 0.993097 and round(col_max, 6) == 0.993680, (col_min, col_max)
    assert abs(worst_spread - 3.2748e-02) <= 5e-7, worst_spread
    assert abs(pec["dsum"] - 3.3903e-04) <= 5e-9 and abs(pec["max_ds"] - 2.2094e-04) <= 5e-9
    assert abs(r["min_blind"] / DISCRETE_IDENTITY_RTOL - 123.3) <= 0.05, r["min_blind"]
    assert abs(min(r["shares"]) - 0.3699) <= 5e-5 and abs(max(r["shares"]) - 0.3737) <= 5e-5
    # example-side map-sum identities and the P_surf-blind negative control
    for name, val in (("Cu phantom", map_cu_ph), ("Cu coil loss", map_cu_surf),
                      ("PEC phantom", map_pec_ph)):
        assert val <= DISCRETE_IDENTITY_RTOL, f"map sum {name}: {val:.3e}"
    assert r["min_blind"] >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, r["min_blind"]
    assert blind_map >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, blind_map

    # ---- combined XDMF -------------------------------------------------------
    order = ["E_magnitude_cu_P17", "phantom_loss_density_cu_P17",
             "coil_surface_loss_per_cell_W_cu_P17", "E_magnitude_pec_P17",
             "phantom_loss_density_pec_P17"]
    fields = {k: store["maps"][k] for k in order}
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags, fields, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    if rank0:
        print(f"\n[output] {xdmf_path.name}: DG0 " + ", ".join(order) + ", CellTags")
        print("EX-80 anchors, map sums and control: PASS")


if __name__ == "__main__":
    main()
