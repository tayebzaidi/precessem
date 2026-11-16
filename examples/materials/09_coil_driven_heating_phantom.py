"""Example (`EX-79`): coil-driven SAR heats the gel phantom — four drives, four hot spots.

`mat:4` / `mat:7` heat a ball under an *imposed* uniform field and `mri:3` shows
the four single-leg drives' SAR without heat.  This example closes the loop on a
coil: `ANS-2`'s 10 MHz loaded F-small birdcage on the 0.0025 m phantom rung
(199 920 cells, 58 866 in the phantom), the four single-leg drives through
``_solve_driven``, the dissipated power ½σ|E|² as a DG0 heat source on the tag-3
submesh at **one** normalisation G = 1 W / P_1, and a Robin-cooled gel
transient (k 0.6, ρc 4.18e6, h 10, T_inf = 20) to t* = 0.2 t_R = 1254 s.

Everything — the fixture body, every band, every record — comes from
``tests/validation/test_mat15_coil_driven_heating_phantom.py`` (`MAT-15` step 1):
``_run_heated(return_fields=True)`` is that gate's own fixture body lifted to a
module-level function (rule (a) lift, EX-79).  Nothing is restated.

Asserted, each printed beside ``20261001T123309Z_MAT-15-step1-W1.log:1924-1954``:
(B) source transfer, V_sub, (E) every run, (E') per drive, (S2), (C4) the four
cyclic 10 g theta pairs; negative controls (d) doubled source and (f) far-side
ball.  Scope: 1 W normalisation (no absolute temperature), degree 1, one rung,
no perfusion, no tissue, no quadrature drive, no convergence claim.

Run through the example runner (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:9 -n 2 -t 480

Output: ``examples/materials/paraview_output/materials_09_coil_driven_heating_phantom_combined.xdmf``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
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
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    B_BAND,
    C4_COVARIANCE_BAND,
    D_CONTROL_FLOOR,
    E_BAND,
    F_CONTROL_FLOOR,
    PHANTOM_CELL_TAG,
    S2_BAND,
    T_INF,
    T_STAR,
    V_BAND,
    _run_heated,
)
from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_09_coil_driven_heating_phantom"
W1 = "20261001T123309Z_MAT-15-step1-W1.log"


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    # Structural checks: the bands this example gates on are the gate modules' objects.
    structural = {
        "m15.E_BAND is m12.E_BAND": m15.E_BAND is m12.E_BAND,
        "m15.S2_BAND is m12.S2_BAND": m15.S2_BAND is m12.S2_BAND,
        "m15.B_BAND is m10.B_BAND": m15.B_BAND is m10.B_BAND,
        "m15.C4_COVARIANCE_BAND is b1map.C4_COVARIANCE_BAND":
            m15.C4_COVARIANCE_BAND is b1map.C4_COVARIANCE_BAND,
        "m15._solve_driven is b1map._solve_driven": m15._solve_driven is b1map._solve_driven,
    }
    for name, ok in structural.items():
        assert ok, f"structural import check failed: {name}"
    if rank0:
        print("=" * 72)
        print("EX-79 — coil-driven SAR heats the gel phantom (MAT-15 step 1, re-executed)")
        print("=" * 72)
        print("[structural] " + "; ".join(f"{k}: {v}" for k, v in structural.items()))
        print(f"[records] every anchor below is beside {W1}:1924-1954", flush=True)

    r = _run_heated(return_fields=True)
    msh, cell_tags, sub = r["msh"], r["cell_tags"], r["sub"]

    # EX-57 setup figure.  The fixture builds the mesh inside ``_run_heated``, so the
    # render follows the call; it touches no printed number.  No-op unless
    # FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            CONDUCTOR_CELL_TAG: "conductor",
            2: "air",
            PHANTOM_CELL_TAG: "phantom (gel)",
            **{PORT_LOWER + i: f"leg L{i} port lower" for i in range(1, LEG_COUNT + 1)},
            **{PORT_UPPER + i: f"leg L{i} port upper" for i in range(1, LEG_COUNT + 1)},
        },
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="mat:9 -- coil-driven SAR heats the gel phantom, four drives (MAT-15 step 1)",
        comm=comm,
    )

    # ---- anchors, asserted at the imported bands ---------------------------
    worst_B = max(r["B"].values())
    worst_Ep = max(r["Ep"].values())
    worst_S2 = max(r["S2"].values())
    worst_E = max(r["E"].values())
    pair_pct = [r["pairs"][k] * 100 for k in range(4)]
    far_pct = [r["far"][k] * 100 for k in range(4)]
    if rank0:
        print(f"\n[anchor] (B) worst rel over 4 drives {worst_B:.3e} <= {B_BAND:.0e}  "
              f"(W1: <= 1.7e-14)")
        print(f"[anchor] V_sub rel {r['err_V']:.3e} <= {V_BAND:.0e}  (W1: 3.447e-13)")
        print(f"[anchor] (E) worst over {len(r['E'])} runs {worst_E:.3e} <= {E_BAND:.0e}  "
              f"(W1: 8.479e-12)")
        print(f"[anchor] (E') worst over 4 drives {worst_Ep:.3e} <= {E_BAND:.0e}")
        print(f"[anchor] (S2) worst |ratio-1| {worst_S2:.4e} <= {S2_BAND:.0e}  (W1: 5.7543e-06)")
        print("[anchor] (C4) theta 10 g pairs % = " + ", ".join(f"{p:.4f}" for p in pair_pct)
              + f" <= {C4_COVARIANCE_BAND * 100:.1f}  (W1: 0.0403, 0.0202, 0.0287, 0.0487)")
        print(f"[control] (d) doubled-source (E') = {r['ctl_d']:.6f} > {D_CONTROL_FLOOR:.0e}  "
              f"(W1: 1.000000)")
        print("[control] (f) far-side ball miss % = " + ", ".join(f"{p:.4f}" for p in far_pct)
              + f" > {F_CONTROL_FLOOR * 100:.1f}  (W1: 75.2462-75.2815)", flush=True)

    for k in range(4):
        assert r["B"][k] <= B_BAND, f"(B) drive {k}: {r['B'][k]:.3e}"
        assert r["Ep"][k] <= E_BAND, f"(E') drive {k}: {r['Ep'][k]:.3e}"
        assert r["S2"][k] <= S2_BAND, f"(S2) drive {k}: {r['S2'][k]:.4e}"
        assert r["pairs"][k] <= C4_COVARIANCE_BAND, f"(C4) pair {k}: {pair_pct[k]:.4f}%"
        assert r["far"][k] > F_CONTROL_FLOOR, f"control (f) drive {k}: {far_pct[k]:.4f}%"
    assert r["err_V"] <= V_BAND, f"V_sub rel {r['err_V']:.3e}"
    for name, val in r["E"].items():
        assert val <= E_BAND, f"(E) run {name}: {val:.3e}"
    assert r["ctl_d"] > D_CONTROL_FLOOR, f"control (d) did not miss: {r['ctl_d']:.3e}"

    # ---- combined XDMF on the phantom submesh -------------------------------
    fields = {}
    for k in range(4):
        T = r["T_star"][k]
        th = fem.Function(T.function_space, name=f"theta_k{k}")
        th.x.array[:] = T.x.array - T_INF
        th.x.scatter_forward()
        fields[f"theta_k{k}"] = th
        q = r["q_h"][k]
        q.name = f"q_h_k{k}"  # the XDMF dataset name is the Function's, not the dict key
        fields[f"q_h_k{k}"] = q
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None, fields, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    if rank0:
        print(f"\n[output] {xdmf_path.name}: P1 theta_k0..k3 (T - {T_INF:g} K at t* = "
              f"{T_STAR:.1f} s) and DG0 q_h_k0..k3 (W/m^3 at G = 1 W/P_1)")
        print("EX-79 anchors and controls: PASS")


if __name__ == "__main__":
    main()
