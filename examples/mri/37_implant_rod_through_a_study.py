"""Example (`EX-167`): a PEC implant rod in the coil-driven F-small phantom through a ``Study``,
against the same ``Study`` with the rod removed.

`WF-15` step 10 (item 709) put the rod on a public surface: ``FourPortBirdcageGeometry(implant_rod=rod)``
(YAML twin: an ``implant:`` section).  No example drives it.  ``imp:3``, ``imp:10`` and ``imp:15`` call
the generators directly.  This one runs

* the rod ``Study`` (degree 1, 64 MHz, four-port ``LumpedPortDrive``, exact-rectangle ports, kappa off);
* the same ``Study`` with ``implant_rod=None``;
* the direct ``_build_rod`` + ``build_four_port_sweep`` as the comparand.

Gates, collected first and then asserted (every band imported by name):

* (Cells) the rod Study's cells equal ``ROD_EXACT_CELL_RECORD`` and the direct build's; the off Study's
  are within ``CELL_COUNT_BAND`` of ``ANS4_EXACT_STRIP_CELL_RECORD``;
* (S) the rod Study's S equals the direct S within ``IDENTITY_BAND``;
* (i) reciprocity within ``RECIPROCITY_BAND``; (ii) sigma_max <= 1 + ``PASSIVITY_SIGMA_TOLERANCE`` and
  column power <= 1;
* (E) max |E| on the rod's PEC wall (tag 402) <= ``PEC_WALL_MAX_E``; wall facets > 0 with the rod, 0 off;
* (c-off), the negative control: the off Study misses the rod S by more than ``IDENTITY_BAND``.

Printed only: the tip 1 g Q per drive (no record; the off arm's phantom cannot hold a 1 g ball, so it is
not a tip ratio) and (M).

Run it through the example runner::

    ./scripts/run_examples.sh -e mri:37 -n 2 -t 560

Not claimed: a tip ratio, a coil-driven tip SAR accuracy, a sweep, degree 2, a tissue model.
"""

from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_implant_rod import ROD, _build_rod  # noqa: E402
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND  # noqa: E402
from tests.mesh.test_step_import_four_port import POOLED_SEPARATION_FLOOR  # noqa: E402
from tests.package.test_wf15_study_implant_rod import ROD_EXACT_CELL_RECORD  # noqa: E402
from tests.validation.test_birdcage_power_identity import STEP2D_FREQUENCY_HZ  # noqa: E402
from tests.validation.test_imp3_palace_twin import _sweep  # noqa: E402
from tests.validation.test_imp3_rod_birdcage import (  # noqa: E402
    PEC_WALL_MAX_E,
    PORT_IDS,
    TIPS,
    _matrix_quantities,
    _wall_dofs,
    _wall_max_e,
)
from tests.validation.test_imp3_rod_heating_phantom import _parent_wall_count  # noqa: E402
from tests.validation.test_lossy_sphere_fullwave import (  # noqa: E402
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_ops61_exact_strip_ports import ANS4_EXACT_STRIP_CELL_RECORD  # noqa: E402
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    PASSIVITY_SIGMA_TOLERANCE,
    build_four_port_sweep,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND  # noqa: E402
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION  # noqa: E402
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_37_implant_rod_through_a_study"
RECORD_WIDTH = 2

# Records of the gate runs, printed beside the readings (never asserted by equality).
REC_WF15_LOG = "20261101T104952Z_WF-15-s10-W1.log:5588-5600"
REC_OFF_MISS = 1.462391e-04          # 20261101T134820Z_WF-17-s25-W1-r2.log (item 717)


def _p0(comm, text):
    if comm.rank == 0:
        print(text, flush=True)


def _study(geometry):
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, STEP2D_FREQUENCY_HZ)
    return Study(geometry, phantom, LumpedPortDrive(), frequency_hz=STEP2D_FREQUENCY_HZ, degree=1)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    assert comm.size == RECORD_WIDTH, f"the record holds at -n {RECORD_WIDTH}, got -n {comm.size}"
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    comm.Barrier()

    # ---- the two Studies and the direct comparand (one window) ------------------------------------
    # The result does not keep the build's facet tags, which (E) needs; a wrapper on `build`
    # captures them (example-side, as in item 709's fixture; no src/ edit).
    captured = {}
    real_build = FourPortBirdcageGeometry.build

    def _capturing_build(self, comm_=MPI.COMM_WORLD):
        out = real_build(self, comm_)
        captured["rod" if self.implant_rod else "off"] = out
        return out

    centres = list(TIPS.values())
    FourPortBirdcageGeometry.build = _capturing_build
    try:
        t0 = time.perf_counter()
        rod_study = _study(FourPortBirdcageGeometry(implant_rod=ROD))
        rod_sar = rod_study.sar(mass_g=1.0, centres=centres, comm=comm)   # one field-keeping solve
        rod = rod_study.fields_result
        t_rod = time.perf_counter() - t0

        t0 = time.perf_counter()
        built = _build_rod(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
        case = build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True)
        direct = _sweep(case, keep=True)
        t_direct = time.perf_counter() - t0

        t0 = time.perf_counter()
        off_study = _study(FourPortBirdcageGeometry())
        off_sar = off_study.sar(mass_g=1.0, centres=centres, comm=comm)
        off = off_study.fields_result
        t_off = time.perf_counter() - t0
    finally:
        FourPortBirdcageGeometry.build = real_build
    _p0(comm, f"[time] -n {comm.size}: rod Study {t_rod:.1f} s (mesh {rod.mesh_time:.1f}, sweep "
              f"{rod.sweep_time:.1f}, incl. its 1 g reads); direct {t_direct:.1f} s; off Study {t_off:.1f} s")

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) ---------------------------------------
    msh, cell_tags = rod.mesh, rod.cell_tags
    local_tags = np.unique(cell_tags.values)
    all_tags = sorted({int(x) for part in comm.allgather(local_tags.tolist()) for x in part})
    phantom_tag = rod_study.geometry.phantom_tag
    other = tuple(x for x in all_tags if x != phantom_tag)
    names = {phantom_tag: "gel phantom (rod hole)"}
    for x in other:
        names[x] = "air / copper regions"
    write_setup_figure(
        msh, cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=other, translucent_tags=(phantom_tag,),
        slice_normal=(0.0, 1.0, 0.0), slice_origin=(0.0, 0.0, 0.0),
        title="mri:37 -- F-small phantom with the PEC rod hole (degree 1)",
        comm=comm)
    gc.collect()
    comm.Barrier()

    # ---- the gates: collect everything, then assert ----------------------------------------------
    mq = _matrix_quantities(rod)
    s_st, s_di = np.asarray(rod.s_matrix), np.asarray(direct.s_matrix)
    rel_s = float(np.linalg.norm(s_st - s_di) / np.linalg.norm(s_di))
    s_off = np.asarray(off.s_matrix)
    miss = float(np.linalg.norm(s_off - s_st) / np.linalg.norm(s_st))
    gates = {
        "cells_record": rod.cells == ROD_EXACT_CELL_RECORD,
        "cells_direct": rod.cells == case["cells"],
        "S_direct": rel_s <= IDENTITY_BAND,
        "i": mq["recip"] <= RECIPROCITY_BAND,
        "ii_sigma": float(mq["sigma"][0]) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE,
        "ii_colpow": bool(np.all(mq["colpow"] <= 1.0)),
    }
    _p0(comm, f"[(Cells)] rod Study {rod.cells}, direct {case['cells']}, record {ROD_EXACT_CELL_RECORD} "
              f"(ASSERTED equal)")
    _p0(comm, f"[(S)] ||S_study - S_direct||/||S_direct|| {rel_s:.3e} (ASSERTED <= IDENTITY_BAND "
              f"{IDENTITY_BAND:.0e}); record {REC_WF15_LOG}")
    _p0(comm, f"[(i)] {mq['recip']:.3e} (ASSERTED <= RECIPROCITY_BAND {RECIPROCITY_BAND:.0e}); "
              f"[(ii)] sigma_max {float(mq['sigma'][0]):.12f} (ASSERTED <= 1 + "
              f"{PASSIVITY_SIGMA_TOLERANCE:g}), column power max {float(np.max(mq['colpow'])):.7f} "
              f"(ASSERTED <= 1)")
    _p0(comm, "[(M)] sigma_y pairs " + ", ".join(f"{v:.3e}" for v in mq["mirror"])
        + f" (printed; band {ADJACENT_SPREAD_BAND})")

    built_rod = captured["rod"]
    v_space = rod.fields["P1"].e_complex.function_space
    owned, n_wall = _wall_dofs(rod.mesh, built_rod["mesh_facet_tags"], v_space)
    wall = {p: _wall_max_e(rod.fields[p].e_complex, owned, comm) for p in PORT_IDS}
    gates["E"] = max(wall.values()) <= PEC_WALL_MAX_E
    gates["wall_dofs"] = n_wall > 0
    _p0(comm, f"[(E)] rod wall (tag 402) dofs {n_wall}; max|E| " +
        ", ".join(f"{p} {v:.3e}" for p, v in wall.items()) +
        f" (ASSERTED <= PEC_WALL_MAX_E {PEC_WALL_MAX_E:g})")

    wall_rod = _parent_wall_count(rod.mesh, built_rod["mesh_facet_tags"])
    wall_off = _parent_wall_count(off.mesh, captured["off"]["mesh_facet_tags"])
    gates["off_cells"] = abs(off.cells / ANS4_EXACT_STRIP_CELL_RECORD - 1.0) <= CELL_COUNT_BAND
    gates["wall_rod"] = wall_rod > 0
    gates["wall_off"] = wall_off == 0
    _p0(comm, f"[(Cells)] off Study {off.cells} (record {ANS4_EXACT_STRIP_CELL_RECORD}, band "
              f"{CELL_COUNT_BAND}); tag-402 owned facets rod {wall_rod}, off {wall_off} (ASSERTED > 0 / 0)")

    gates["c_off"] = miss > IDENTITY_BAND
    _p0(comm, f"[(c-off)] the off Study misses the rod S by {miss:.6e} = {miss / IDENTITY_BAND:.3e} x "
              f"IDENTITY_BAND {IDENTITY_BAND:.0e} (ASSERTED >); item 717's miss {REC_OFF_MISS:.6e}; "
              f"POOLED_SEPARATION_FLOOR {POOLED_SEPARATION_FLOOR:g} (printed)")

    tips = list(TIPS)
    for pid in PORT_IDS:
        row = []
        for i, tip in enumerate(tips):
            a = float(rod_sar[pid][i]["averaged_sar_w_per_kg"])
            b = float(off_sar[pid][i]["averaged_sar_w_per_kg"])
            q = a / b if b != 0.0 else float("nan")
            row.append(f"{tip} {q:.4f} (rod {a:.4e}, off {b:.4e} W/kg)")
        _p0(comm, f"[tip 1 g Q] {pid} (a reading, no record; NOT a tip ratio): " + "; ".join(row))

    red = sorted(comm.allreduce([k for k, v in gates.items() if not v], op=MPI.SUM))
    _p0(comm, f"[gates] red: {red or 'none'} of {len(gates)}")

    # ---- ParaView: |E| of the P1 drive with the rod, sigma, cell tags ------------------------------
    e = rod.fields["P1"].e_complex
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    e_abs = fem.Function(p1_space, name="E_abs")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1_space.element.interpolation_points))
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()
    dg0 = fem.functionspace(msh, ("DG", 0))
    sigma_out = fem.Function(dg0, name="sigma")
    sigma_out.interpolate(rod.fields["P1"].sigma_field)
    sigma_out.x.array[:] = np.real(sigma_out.x.array)
    sigma_out.x.scatter_forward()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, cell_tags,
        {"E_abs": e_abs, "sigma": sigma_out}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _p0(comm, f"[paraview] wrote {xdmf_path} (|E| of the P1 drive with the rod, sigma, cell tags)")

    assert not red, red
    _p0(comm, f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
