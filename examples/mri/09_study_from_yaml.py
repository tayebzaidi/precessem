"""Example (`EX-93`): the same studies from YAML -- ``femem run`` on two files gives
the Python route's outputs (`WF-1`).

`WF-1`'s ``fem_em_solver.cli`` reads a study from a YAML file and runs it. This
example drives ``cli.main(["run", <yaml>, "--output-dir", <dir>])`` in-process on
two YAML files that sit beside it -- `MAT-4` step 1's lossy sphere and `WF-16`'s
F-small four-port birdcage, both at degree 1 -- and reads the outputs against the
Python route's records.

Every constant, band and helper is imported from the gate modules (the ANS-1
rule); the example restates none.

Asserted:

* (K)        each YAML is spec-equal (``float.hex``) to its ``tests/package/wf1/``
             original and to the Python-built ``Study`` (``_spec_diff == []``);
* (R-sphere) ``run`` returns 0 on every rank and ``summary.json``'s mean SAR is
             within ``L0_BAND`` of ``MAT4_DEG1_RECORD``; degree, ranks and cells match;
* (R-S)      the F-small ``summary.json`` S matrix is within
             ``LEG_D0_REPRODUCTION_BAND`` of ``P1 = _record(1)``.

Negative controls (both asserted):

* (X)  the unknown-key file raises ``StudyConfigError`` naming
       ``geometry.sphere_raduis``; ``run`` returns 2 on every rank and writes no
       output directory;
* (c3) the S comparator is not blind: ``_rel(S_yaml, _record(2))`` exceeds
       ``LEG_D0_REPRODUCTION_BAND`` (no solve).

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:9 -n 2 -t 300

Not claimed: a stable YAML schema, degree 2, the installed-wheel ``femem`` entry
point, Touchstone accuracy, or anything `WF-1` / `WF-15` / `WF-16` do not claim.
"""

from __future__ import annotations

import json
import shutil
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

from fem_em_solver import StudyConfigError, load_study  # noqa: E402
from fem_em_solver import cli  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.package.test_wf15_study_lossy_sphere import MAT4_DEG1_RECORD, _study  # noqa: E402
from tests.package.test_wf16_study_four_port import _record, _rel  # noqa: E402
from tests.package.test_wf1_study_yaml import (  # noqa: E402
    BIRDCAGE_FILE,
    MALFORMED,
    SPHERE_FILE,
    WF1,
    _py_birdcage,
    _rc_all,
    _spec_diff,
)
from tests.validation.test_ans4_exact_strip_records import (  # noqa: E402
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND  # noqa: E402
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    LEG_D0_REPRODUCTION_BAND,
)

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "paraview_output"
FIGURE_DIR = HERE / "figures"
BASENAME = "mri_09_study_from_yaml"
SPHERE_YAML = HERE / "09_study_from_yaml_sphere.yaml"
FOUR_PORT_YAML = HERE / "09_study_from_yaml_four_port.yaml"
MALFORMED_NAME = "malformed_unknown_key.yaml"
AIR_TAG = 2


def _fresh(path: Path, comm) -> None:
    if comm.rank == 0 and path.exists():
        shutil.rmtree(path)
    comm.Barrier()


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"record width -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}"
    )
    if comm.rank == 0:
        OUTPUT_DIR.mkdir(exist_ok=True)
    comm.Barrier()

    # ---- (K): the YAML specs equal the originals and the Python studies ----
    study = load_study(SPHERE_YAML)
    d_sphere_orig = _spec_diff(study, load_study(SPHERE_FILE))
    d_cage_orig = _spec_diff(load_study(FOUR_PORT_YAML), load_study(BIRDCAGE_FILE))
    d_sphere = _spec_diff(study, _study(degree=1))
    d_cage = _spec_diff(load_study(FOUR_PORT_YAML), _py_birdcage())

    # ---- (X): the unknown-key file ------------------------------------------
    bad = WF1 / MALFORMED_NAME
    expected_field = MALFORMED[MALFORMED_NAME]
    try:
        load_study(bad)
        x_message = None
    except StudyConfigError as exc:
        x_message = str(exc)
        x_field = exc.field
    bad_out = OUTPUT_DIR / f"{BASENAME}_malformed"
    _fresh(bad_out, comm)
    rc_bad = cli.main(["run", str(bad), "--output-dir", str(bad_out)])
    bad_rc = _rc_all(rc_bad, comm)
    bad_exists = comm.allreduce(bool(bad_out.exists()), op=MPI.LOR)

    # ---- (R-sphere): femem run on the sphere YAML ---------------------------
    sphere_out = OUTPUT_DIR / f"{BASENAME}_sphere"
    _fresh(sphere_out, comm)
    t0 = time.perf_counter()
    rc_s = cli.main(["run", str(SPHERE_YAML), "--output-dir", str(sphere_out)])
    t_sphere = time.perf_counter() - t0
    sphere_rc = _rc_all(rc_s, comm)
    comm.Barrier()
    with open(sphere_out / "summary.json") as fh:
        s_summary = json.load(fh)
    sar_yaml = float(s_summary["mean_sar"]["mean_sar_w_per_kg"])
    rel_sphere = abs(sar_yaml - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD

    # ---- (R-S): femem run on the four-port YAML -----------------------------
    cage_out = OUTPUT_DIR / f"{BASENAME}_four_port"
    _fresh(cage_out, comm)
    t0 = time.perf_counter()
    rc_c = cli.main(["run", str(FOUR_PORT_YAML), "--output-dir", str(cage_out)])
    t_cage = time.perf_counter() - t0
    cage_rc = _rc_all(rc_c, comm)
    comm.Barrier()
    with open(cage_out / "summary.json") as fh:
        c_summary = json.load(fh)
    s_yaml = np.array([[complex(*p) for p in row] for row in c_summary["s_matrix"]])
    p1 = _record(1)
    rel_s = _rel(s_yaml, p1)
    c3 = _rel(s_yaml, _record(2))

    # ---- the loaded sphere Study, solved once for the XDMF + figure ---------
    t0 = time.perf_counter()
    res = study.solve(comm)
    solve_seconds = time.perf_counter() - t0
    msh = res.mesh
    n_cells = int(msh.topology.index_map(msh.topology.dim).size_global)
    sigma_value = study.material.sigma

    write_setup_figure(
        msh,
        res.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            study.geometry.sphere_tag: f"gel sphere (sigma {sigma_value:g} S/m)",
            AIR_TAG: "air box (imposed-field walls)",
        },
        hide_tags=(AIR_TAG,),
        translucent_tags=(study.geometry.sphere_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:9 -- the sphere Study loaded from YAML (WF-1)",
        comm=comm,
    )

    e = res.fields.e_complex
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    e_abs = fem.Function(p1_space, name="E_abs")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1_space.element.interpolation_points)
    )
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()

    dg0 = fem.functionspace(msh, ("DG", 0))
    sar_field = fem.Function(dg0, name="SAR")
    sar_field.interpolate(
        fem.Expression(
            res.fields.sigma_field * ufl.real(ufl.inner(e, e)) / (2.0 * study.density_kg_m3),
            dg0.element.interpolation_points,
        )
    )
    sar_field.x.array[:] = np.real(sar_field.x.array)
    sar_field.x.scatter_forward()

    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        res.cell_tags,
        {"E_abs": e_abs, "SAR": sar_field},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- prints (rank 0, no collective inside) ------------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-93 -- the same studies from YAML (WF-1)")
        print("=" * 72)
        print(
            f"[(K)] sphere vs original {d_sphere_orig}; F-small vs original {d_cage_orig}; "
            f"sphere vs Python {d_sphere}; F-small vs Python {d_cage}"
        )
        print(
            f"[(R-sphere)] run rc (min, max) {sphere_rc}, {t_sphere:.1f} s; SAR = {sar_yaml!r} "
            f"({sar_yaml:.10e}); record {MAT4_DEG1_RECORD:.10e}; rel {rel_sphere:.3e} "
            f"(band L0_BAND {L0_BAND:g})"
        )
        print(
            f"[(R-sphere)] summary degree {s_summary['degree']}, ranks {s_summary['ranks']}, "
            f"cells {s_summary['cells']} / solved {n_cells} (printed; version-tagged 74020); "
            f"loaded-Study solve {solve_seconds:.1f} s"
        )
        print(
            f"[(R-S)] run rc (min, max) {cage_rc}, {t_cage:.1f} s; ||S_yaml - P1||/||P1|| "
            f"{rel_s:.3e} (band LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:g}); "
            f"cells {c_summary['cells']} (printed; version-tagged 117206)"
        )
        print(f"[(R-S)] touchstone {cage_out / c_summary['touchstone']}")
        print(f"[(X)] field {x_field!r} (expected {expected_field!r}); message: {x_message}")
        print(f"[(X)] run rc (min, max) {bad_rc}; output dir exists {bad_exists}")
        print(
            f"[(c3)] ||S_yaml - _record(2)||/||_record(2)|| = {c3:.4f} "
            f"(must exceed band {LEG_D0_REPRODUCTION_BAND:g})"
        )
        print(f"[paraview] wrote {xdmf_path}")

    assert d_sphere_orig == [], d_sphere_orig
    assert d_cage_orig == [], d_cage_orig
    assert d_sphere == [], d_sphere
    assert d_cage == [], d_cage
    assert sphere_rc == (0, 0) and cage_rc == (0, 0)
    assert s_summary["degree"] == 1 and s_summary["ranks"] == 2
    assert s_summary["cells"] == n_cells
    assert rel_sphere <= L0_BAND, f"(R-sphere) rel {rel_sphere:.3e} > {L0_BAND:g}"
    assert c_summary["degree"] == 1 and c_summary["ranks"] == 2
    assert rel_s <= LEG_D0_REPRODUCTION_BAND, f"(R-S) rel {rel_s:.3e}"
    assert x_message is not None and x_field == expected_field and expected_field in x_message
    assert bad_rc == (2, 2) and not bad_exists, (bad_rc, bad_exists)
    assert c3 > LEG_D0_REPRODUCTION_BAND, f"(c3) blind: {c3:.4f}"

    if comm.rank == 0:
        print(
            f"\nAll anchors and controls hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
