"""Example (`EX-98`): the container wall on the coil-driven gel phantom, four drives.

`mat:12` has a conforming wall on a ball under an imposed field; `mat:9` has the
coil-driven phantom with no wall.  This example puts the wall on the coil-driven
phantom the way `MAT-18` step 3 did: **a series resistance in the Robin coefficient,
not geometry** -- h_eff = 1 / (1/h + t_w/k2), planar, one scalar on the lateral faces
and the caps of `MAT-15` P1's tag-3 submesh, which has no wall cells and no wall heat
capacity.

Everything -- the fixture body, every band, every record -- comes from
``tests/validation/test_mat18_wall_coil_driven.py`` (`MAT-18` step 3):
``_run_walled(return_fields=True)`` is that gate's fixture body lifted to a
module-level function (rule (a) lift, `EX-98`; the gate re-ran green after it).
Nothing is restated; the step-3 record values appear only inside print strings.

Asserted: (R), (E') per drive, (E) every wall run, (L0) helper equivalence, (C4) the four
cyclic 10 g theta pairs, controls (d) doubled source and (W-ctl) the wall is felt.
Printed, not asserted: the wall effect (the +4.63 % move) with a ``.6f`` string-equality
verdict against the step-3 record, and the (h-swap) prediction.

Scope: 1 W normalisation (no absolute temperature), degree 1, one rung, untuned 10 MHz,
sigma = 800 S/m F-small birdcage; no conforming wall, no wall heat storage, no contact
resistance, no convergence claim.

Run through the example runner (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:14 -n 2 -t 590

Output: ``examples/materials/paraview_output/materials_14_coil_driven_phantom_in_container_wall_combined.xdmf``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.setup_figure import write_setup_figure

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.mesh.test_birdcage_port_tags import LEG_COUNT  # noqa: E402
from tests.validation import test_birdcage_b1_plus_map as b1map  # noqa: E402
from tests.validation import test_mat10_sar_to_heat_sphere as m10  # noqa: E402
from tests.validation import test_mat12_transient_sar_heating_sphere as m12  # noqa: E402
from tests.validation import test_mat15_coil_driven_heating_phantom as m15  # noqa: E402
from tests.validation import test_mat18_two_layer_sphere_closed_form as m18  # noqa: E402
from tests.validation import test_mat18_wall_coil_driven as m18s3  # noqa: E402
from tests.validation.test_mat18_wall_coil_driven import (  # noqa: E402
    B_BAND,
    C4_COVARIANCE_BAND,
    D_CONTROL_FLOOR,
    E_BAND,
    H_GEL,
    R_BAND,
    ROBIN_TAG,
    T_INF,
    T_STAR,
    W_CONTROL_FLOOR,
    _e_prime_h,
    _run_walled,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import PHANTOM_CELL_TAG  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_14_coil_driven_phantom_in_container_wall"
W1 = "20261006T200649Z_MAT-18-step3-W1.log"
# Printed-only records (W1:1960-1972); never asserted against.
RECORD_RATIO = ("1.046316", "1.046314", "1.046316", "1.046314")
RECORD_MEAN_W = ("1.206730553", "1.205637378", "1.206047101", "1.205345805")
RECORD_MEAN_P1 = ("1.153313650", "1.152271459", "1.152660000", "1.151992598")


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
        "m18s3.E_BAND is m12.E_BAND": m18s3.E_BAND is m12.E_BAND,
        "m18s3.B_BAND is m10.B_BAND": m18s3.B_BAND is m10.B_BAND,
        "m18s3.C4_COVARIANCE_BAND is b1map.C4_COVARIANCE_BAND":
            m18s3.C4_COVARIANCE_BAND is b1map.C4_COVARIANCE_BAND,
        "m18s3.D_CONTROL_FLOOR is m15.D_CONTROL_FLOOR":
            m18s3.D_CONTROL_FLOOR is m15.D_CONTROL_FLOOR,
        "m18s3._run_heated is m15._run_heated": m18s3._run_heated is m15._run_heated,
        "m18s3.K2 is m18.K2": m18s3.K2 is m18.K2,
    }
    for name, ok in structural.items():
        assert ok, f"structural import check failed: {name}"
    if rank0:
        print("=" * 72)
        print("EX-98 -- the container wall on the coil-driven gel phantom "
              "(MAT-18 step 3, re-executed)")
        print("=" * 72)
        print("[structural] " + "; ".join(f"{k}: {v}" for k, v in structural.items()))
        print(f"[records] every anchor below is beside {W1}", flush=True)

    r = _run_walled(return_fields=True)
    assert r["R"] <= R_BAND, f"(R) {r['R']:.3e} > {R_BAND:.0e}"
    base = r["base"]
    msh, cell_tags, sub, ft = base["msh"], base["cell_tags"], base["sub"], base["ft"]
    h_eff, h_lat = r["h_eff"], r["h_lat"]

    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            CONDUCTOR_CELL_TAG: "conductor",
            2: "air",
            PHANTOM_CELL_TAG: (f"phantom (gel; wall as series resistance, Robin h_eff "
                               f"{h_eff:.4g} W/m^2/K -- no wall cells)"),
            **{PORT_LOWER + i: f"leg L{i} port lower" for i in range(1, LEG_COUNT + 1)},
            **{PORT_UPPER + i: f"leg L{i} port upper" for i in range(1, LEG_COUNT + 1)},
        },
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title=("mat:14 -- container wall on the coil-driven gel phantom: series "
               "resistance h_eff, four drives (MAT-18 step 3)"),
        comm=comm,
    )

    # ---- anchors and controls, asserted at the imported bands ------------------
    worst_Ep = max(r["Ep"].values())
    worst_E = max(r["E"].values())
    worst_L0 = max(r["L0"].values())
    pair_pct = [r["pairs"][k] * 100 for k in range(4)]
    qpair_pct = [r["q_pairs"][k] * 100 for k in range(4)]
    for k in range(4):
        assert r["Ep"][k] <= E_BAND, f"(E') drive {k}: {r['Ep'][k]:.3e}"
        assert r["L0"][k] <= B_BAND, f"(L0) drive {k}: {r['L0'][k]:.3e}"
        assert r["pairs"][k] <= C4_COVARIANCE_BAND, f"(C4) pair {k}: {pair_pct[k]:.4f}%"
        assert r["W"][k] > W_CONTROL_FLOOR, f"control (W-ctl) drive {k}: {r['W'][k]:.3e}"
    for name, val in r["E"].items():
        assert val <= E_BAND, f"(E) run {name}: {val:.3e}"
    assert r["ctl_d"] > D_CONTROL_FLOOR, f"control (d) did not miss: {r['ctl_d']:.3e}"

    # (h-swap): predicted, printed, never asserted -- the balance sees which h was used.
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    G, P = base["G"], base["P"]
    hswap, _, _ = _e_prime_h(sub, ds, r["wall_runs"][0], G * P[0], H_GEL)

    ratio_str = [f"{r['ratio'][k]:.6f}" for k in range(4)]
    verdict = tuple(ratio_str) == RECORD_RATIO
    kpk = max(r["th_own"], key=r["th_own"].get)

    if rank0:
        print(f"\n[wall] h_eff (planar, used) = {h_eff:.9f} W/m^2/K (W1: 6.551724138); "
              f"lateral cylindrical {h_lat:.9f} (printed, not used; W1: 8.304028064)")
        print(f"[anchor] (R) {r['R']:.3e} <= {R_BAND:.0e}  (W1: 4.985e-16)")
        print("[anchor] (E') per drive " + ", ".join(f"{r['Ep'][k]:.3e}" for k in range(4))
              + f" <= {E_BAND:.0e}  (W1: 6.661e-15, 2.267e-14, 1.133e-14, 9.337e-15)")
        print(f"[anchor] (E) worst over {len(r['E'])} wall runs {worst_E:.3e} <= {E_BAND:.0e}  "
              f"(W1: 4.554e-14)")
        print("[anchor] (L0) helper equivalence (t_w = 0 through the wall path is P1's code "
              "path; not an independent wall -> 0 limit) per drive "
              + ", ".join(f"{r['L0'][k]:.3e}" for k in range(4))
              + f" <= {B_BAND:.0e}  (W1: 0, 6.661e-16, 2.220e-16, 2.220e-16)")
        print("[anchor] (C4) wall theta 10 g pairs % = " + ", ".join(f"{p:.4f}" for p in pair_pct)
              + f" <= {C4_COVARIANCE_BAND * 100:.1f}  (W1: 0.0402, 0.0210, 0.0292, 0.0484)")
        print("[print]  q_h pairs % = " + ", ".join(f"{p:.4f}" for p in qpair_pct))
        print(f"[control] (d) doubled-source (E') with wall = {r['ctl_d']:.6f} > "
              f"{D_CONTROL_FLOOR:.0e}  (W1: 1.000000)")
        print("[control] (W-ctl) wall is felt, |ratio - 1| per drive "
              + ", ".join(f"{r['W'][k]:.4e}" for k in range(4))
              + f" > {W_CONTROL_FLOOR:.0e}  (W1: 4.631e-2 on all four)")
        print(f"[control] (h-swap) k=0 wall run with the bare gel h in the surface term: "
              f"(E') = {hswap:.4e} (PREDICTED 8.697e-02, printed, never asserted)")
        for k in range(4):
            print(f"[effect] k={k} <theta>(t*) wall {r['mean_w'][k]!r} K  P1 {r['mean_p1'][k]!r} K "
                  f"(W1: {RECORD_MEAN_W[k]} / {RECORD_MEAN_P1[k]});  ratio {ratio_str[k]} = "
                  f"{100 * (r['ratio'][k] - 1.0):+.4f} % (W1 ratio {RECORD_RATIO[k]})")
        print(f"[effect] .6f string-equality verdict vs the step-3 record ratios: {verdict}")
        print(f"[print]  peak 10 g theta with wall {r['th_own'][kpk]:.9e} K on k={kpk}, P1 "
              f"{r['th_p1'][kpk]:.9e} K, ratio {r['th_own'][kpk] / r['th_p1'][kpk]:.6f} "
              f"(W1: 1.029608 on k=0)")
        print(f"[print]  max|Im theta| = {r['max_imag']:.3e} (W1: 0);  global phantom sub-cells "
              f"{sub.topology.index_map(3).size_global}", flush=True)

    # ---- combined XDMF on the phantom submesh -----------------------------------
    fields = {}
    for k in range(4):
        wall_T = r["wall_runs"][k].temperature
        V = wall_T.function_space
        tw = fem.Function(V, name=f"theta_wall_k{k}")
        tw.x.array[:] = wall_T.x.array - T_INF
        tw.x.scatter_forward()
        tb = fem.Function(base["T_star"][k].function_space, name=f"theta_bare_k{k}")
        tb.x.array[:] = base["T_star"][k].x.array - T_INF
        tb.x.scatter_forward()
        dt_ = fem.Function(V, name=f"dtheta_k{k}")
        dt_.x.array[:] = tw.x.array - tb.x.array
        dt_.x.scatter_forward()
        q = base["q_h"][k]
        q.name = f"q_h_k{k}"
        fields[tw.name], fields[tb.name], fields[dt_.name], fields[q.name] = tw, tb, dt_, q
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None, fields, comm=comm, facet_tags=ft
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    if rank0:
        print(f"\n[output] {xdmf_path.name}: P1 theta_wall_k0..k3, theta_bare_k0..k3, "
              f"dtheta_k0..k3 (T - {T_INF:g} K at t* = {T_STAR:.1f} s) and DG0 q_h_k0..k3 "
              f"(W/m^3 at G = 1 W/P_1); Robin facets as their own block")
        print("EX-98 anchors and controls: PASS")


if __name__ == "__main__":
    main()
