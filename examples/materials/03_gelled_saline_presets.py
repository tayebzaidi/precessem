"""Example (`EX-65`): the low / mid / high gelled-saline presets side by side on
the lossy sphere, through the ``phantom_material`` path.

`MAT-1` (re-anchored 2026-09-30) gates the three ``GelledSalinePhantomMaterial``
presets, passed as ``phantom_material`` with ``phantom_tag = SPHERE_TAG``,
against the lossy-sphere closed form: mean SAR and interior Im/Re E_z each
within `MAT-4`'s 10 % band. `mat:2` sets the sphere by ``material_map`` at two
conductivities; `mri:1` hand-types a phantom and asserts nothing. No example
calls ``from_preset`` or anchors the ``phantom_material`` path every
coil-and-phantom workflow uses. This one solves all three presets on one mesh
and shows the counter-intuitive reading: on an **imposed field** the more
conductive gel absorbs *less* (SAR_high / SAR_low = 0.57 where a sigma-blind
field says 2.58).

**Import, never restate.** ``_solve_case``, ``_closed_case``, ``MAT4_BAND``,
``_PREDICTED`` and the fixture constants are imported from
``tests/materials/test_phantom_material_model.py``. The one change to that
module is additive and disclosed: the mesh build is lifted verbatim into a
module-level ``build_preset_sphere_mesh(comm)`` that the gate calls too. The
gate module is re-run green in the same slot (`EX-65-mat1-regate`).

**Anchors (asserted):** per preset, mean SAR and interior Im/Re E_z within
``MAT4_BAND`` of ``_closed_case``. The gate's printed record (SAR misses
3.523 / 3.548 / 3.557 %, Im/Re 0.404 / 0.835 / 1.499 %, ratio 0.5693) is printed
beside the measurement; a difference is reported, never asserted.

**Negative control (asserted only as "misses ``MAT4_BAND``"):** each preset's FEM
SAR against each *other* preset's closed form: six misses. The predictions
(from the closed forms and the gate's own printed misses) are printed beside
the measured values, not asserted.

Scope: imposed uniform field, 64 MHz, degree 1, one mesh rung. No coil, no mass
averaging, no accuracy claim against measured gels.

Run it through the example runner (the ``mat:`` group sources the complex
build automatically)::

    ./run_examples.sh -e mat:3 -n 2 -t 300
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

import dolfinx
from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial  # noqa: E402
from fem_em_solver.post.sar import point_sar  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.materials.test_phantom_material_model import (  # noqa: E402
    E0,
    FREQUENCY_HZ,
    MAT4_BAND,
    RHO_KG_M3,
    SPHERE_TAG,
    _PREDICTED,
    _closed_case,
    _solve_case,
    build_preset_sphere_mesh,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_03_gelled_saline_presets_combined"
TABLE_NAME = "materials_03_gelled_saline_presets_table.csv"
FIGURE_BASENAME = "materials_03_gelled_saline_presets"

NAMES = ("low", "mid", "high")

#: Plan figures (`20260930T003229Z_MAT-1-reanchor-W1.log:496-499`): printed
#: beside the measurement, never asserted.
PLAN_SAR_ERR_PCT = {"low": 3.523, "mid": 3.548, "high": 3.557}
PLAN_IM_RE_ERR_PCT = {"low": 0.404, "mid": 0.835, "high": 1.499}
PLAN_RATIO = (0.5693, 0.5692, 2.5788)  # FEM, closed, sigma-blind


def _sar_dg0(run, sigma, msh, cell_tags, name):
    """Pointwise sigma|E|^2/(2 rho) at owned cell centroids, zero off the
    sphere tag. `point_sar` is collective over one point list identical on
    every rank (`mat:2` `_build_dg0_sar_fields`)."""
    comm = msh.comm
    n_local = msh.topology.index_map(msh.topology.dim).size_local
    cells = np.arange(n_local, dtype=np.int32)
    centroids = dolfinx.mesh.compute_midpoints(msh, msh.topology.dim, cells)
    gathered = comm.allgather(centroids)
    global_points = np.vstack(gathered) if len(gathered) > 1 else np.atleast_2d(gathered[0])
    offset = int(sum(chunk.shape[0] for chunk in gathered[: comm.rank]))
    sar_global = point_sar(
        run["fields"].e_real,
        run["fields"].e_imag,
        global_points,
        sigma=sigma,
        rho=RHO_KG_M3,
        comm=comm,
    )
    values = sar_global[offset : offset + n_local]
    dg0 = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(dg0, name=name)
    f.x.array[:] = 0.0
    sphere_cells = cell_tags.indices[cell_tags.values == SPHERE_TAG]
    for cell in sphere_cells:
        if cell < n_local:
            dof = dg0.dofmap.cell_dofs(int(cell))
            f.x.array[dof] = values[int(cell)]
    f.x.scatter_forward()
    return f


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    if comm.rank == 0:
        print("=" * 72)
        print("EX-65 -- the low / mid / high gelled-saline presets on the lossy sphere")
        print("=" * 72, flush=True)

    msh, cell_tags, _ = build_preset_sphere_mesh(comm)
    n_global = msh.topology.index_map(msh.topology.dim).size_global

    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{FIGURE_BASENAME}_setup.png",
        region_names={SPHERE_TAG: "sphere (gelled saline preset)", 2: "air"},
        hide_tags=(2,),
        translucent_tags=(SPHERE_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:3 -- gelled-saline presets on the lossy sphere (MAT-1 fixture)",
        comm=comm,
    )

    solve_started = time.perf_counter()
    rows = {}
    for name in NAMES:
        mat = GelledSalinePhantomMaterial.from_preset(name, frequency_hz=FREQUENCY_HZ)
        run = _solve_case(msh, cell_tags, mat.sigma, mat.epsilon_r, phantom=mat)
        c = _closed_case(mat.sigma, mat.epsilon_r)
        im_re = run["ez_mean"].imag / run["ez_mean"].real
        rows[name] = {
            "mat": mat,
            "closed": c,
            "run": run,
            "sar_err": abs(run["mean_sar"] - c["sar_w_per_kg"]) / c["sar_w_per_kg"],
            "im_re": im_re,
            "im_re_err": abs(im_re - c["im_re"]) / abs(c["im_re"]),
            "sar_blind": mat.sigma * (3.0 * E0 / (mat.epsilon_r + 2.0)) ** 2 / (2.0 * RHO_KG_M3),
        }
    solve_seconds = time.perf_counter() - solve_started

    ratio_fem = rows["high"]["run"]["mean_sar"] / rows["low"]["run"]["mean_sar"]
    ratio_closed = rows["high"]["closed"]["sar_w_per_kg"] / rows["low"]["closed"]["sar_w_per_kg"]
    ratio_blind = rows["high"]["sar_blind"] / rows["low"]["sar_blind"]

    # Negative control: each preset's FEM SAR against each OTHER preset's closed form.
    control = {}
    for a in NAMES:
        for b in NAMES:
            if a == b:
                continue
            fem_sar = rows[a]["run"]["mean_sar"]
            closed_b = rows[b]["closed"]["sar_w_per_kg"]
            measured = abs(fem_sar - closed_b) / closed_b
            # prediction: gate-record miss of `a`, closed forms from the gate's _PREDICTED
            pred = abs(
                (1.0 + PLAN_SAR_ERR_PCT[a] / 100.0) * _PREDICTED[a][0] / _PREDICTED[b][0] - 1.0
            )
            control[(a, b)] = (measured, pred)

    if comm.rank == 0:
        print(f"\n[mesh] {n_global} cells, f = {FREQUENCY_HZ / 1e6:.3f} MHz, degree 1, phantom_tag = {SPHERE_TAG}")
        print(f"[solve] {solve_seconds:.1f} s for three preset solves", flush=True)
        for name in NAMES:
            r = rows[name]
            p = _PREDICTED[name]
            print(
                f"  preset {name:4s} sigma = {r['mat'].sigma:.2f}, eps_r = {r['mat'].epsilon_r:.0f}: "
                f"|k_in|R = {r['closed']['k_in_r']:.4f} (pred {p[2]}), "
                f"closed SAR = {r['closed']['sar_w_per_kg']:.4e} (pred {p[0]:.3g}), "
                f"FEM SAR = {r['run']['mean_sar']:.4e} err {r['sar_err']:.3%} "
                f"(record {PLAN_SAR_ERR_PCT[name]:.3f}%); "
                f"Im/Re closed {r['closed']['im_re']:.4f} (pred {p[1]}), FEM {r['im_re']:.4f} "
                f"err {r['im_re_err']:.3%} (record {PLAN_IM_RE_ERR_PCT[name]:.3f}%)"
            )
        print(
            f"  ratio SAR_high/SAR_low: FEM {ratio_fem:.4f} (record {PLAN_RATIO[0]}), "
            f"closed {ratio_closed:.4f} (record {PLAN_RATIO[1]}), "
            f"sigma-blind {ratio_blind:.4f} (record {PLAN_RATIO[2]}) -- printed, not asserted"
        )
        for (a, b), (measured, pred) in control.items():
            print(
                f"  [control] FEM SAR({a}) vs closed({b}): miss {measured:.3%} "
                f"(predicted {pred:.3%}, band {MAT4_BAND:.0%})"
            )
        print(flush=True)

    # ---- anchors, asserted --------------------------------------------------
    for name in NAMES:
        r = rows[name]
        assert r["sar_err"] < MAT4_BAND, (
            f"preset {name}: mean SAR {r['run']['mean_sar']:.4e} vs closed "
            f"{r['closed']['sar_w_per_kg']:.4e} ({r['sar_err']:.2%}) -- over MAT4_BAND"
        )
        assert r["im_re_err"] < MAT4_BAND, (
            f"preset {name}: Im/Re E_z {r['im_re']:.4f} vs closed "
            f"{r['closed']['im_re']:.4f} ({r['im_re_err']:.2%}) -- over MAT4_BAND"
        )

    # ---- negative control, asserted only as "misses MAT4_BAND" ---------------
    for (a, b), (measured, _pred) in control.items():
        assert measured > MAT4_BAND, (
            f"control: FEM SAR({a}) vs closed({b}) misses by only {measured:.2%} "
            f"-- inside MAT4_BAND, the band cannot tell the presets apart"
        )

    # ---- outputs -------------------------------------------------------------
    fields = {
        f"SAR_{name}": _sar_dg0(rows[name]["run"], rows[name]["mat"].sigma, msh, cell_tags, f"SAR_{name}")
        for name in NAMES
    }
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    write_xdmf_with_tags(OUTPUT_DIR / BASENAME, msh, cell_tags, fields, comm=comm)
    if comm.rank == 0:
        with open(OUTPUT_DIR / TABLE_NAME, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(
                ["preset", "sigma_S_per_m", "eps_r", "closed_sar_w_per_kg", "fem_sar_w_per_kg",
                 "sar_err", "closed_im_re", "fem_im_re", "im_re_err"]
            )
            for name in NAMES:
                r = rows[name]
                w.writerow(
                    [name, r["mat"].sigma, r["mat"].epsilon_r, r["closed"]["sar_w_per_kg"],
                     r["run"]["mean_sar"], r["sar_err"], r["closed"]["im_re"], r["im_re"], r["im_re_err"]]
                )
    adopt_host_ownership(OUTPUT_DIR, comm)

    elapsed = time.perf_counter() - started
    if comm.rank == 0:
        print(
            f"[write] {BASENAME}.xdmf: SAR_low, SAR_mid, SAR_high (DG0, zero off the sphere tag), CellTags"
            f"\n[write] {TABLE_NAME}"
            f"\n\nEX-65 complete in {elapsed:.1f} s. All anchors and the six-miss control held.",
            flush=True,
        )


if __name__ == "__main__":
    main()
