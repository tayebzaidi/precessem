"""Example (`EX-92`): four-port S from one ``Study`` -- the F-small birdcage on exact
port rectangles, with its Touchstone file.

`WF-16`'s ``Study`` composes ``FourPortBirdcageGeometry.f_small()`` (exact CAD port
rectangles, kappa off, operator ruling 2026-10-02), a gelled-saline phantom and a
``LumpedPortDrive``, solves the four drives at 64 MHz and returns S, Z and a
Touchstone file.  This example builds it at degree 1 on ``-n 2`` -- the width and
degree of item 192's registered record -- and holds S, the file read back through
``scikit-rf`` and `PORT-11`'s identities to the gate modules' own bands.

Every constant, band and helper is imported from the gate modules (the ANS-1
rule); the example states no sigma, epsilon_r, frequency, geometry dimension or
impedance of its own.

Asserted:

* (W) ``comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH`` before any S comparison;
* (R) ``_rel(S, P1) <= LEG_D0_REPRODUCTION_BAND`` (P1 = the degree-1 record);
* (G) `PORT-11`'s reciprocity / passivity / class-spread / pooled-separation gates
  on the Study's own S and Z;
* (T) the ``.s4p`` read back through ``skrf.Network``: max |S_skrf - S| <=
  ``ANCHOR_TOL``, frequency, 4 ports, z0 from ``LumpedPortDrive``, and
  ``_rel(S_file, P1) <= LEG_D0_REPRODUCTION_BAND``;
* (X) a ``keep_fields=True`` sweep on the Study's own problem reproduces the
  Study's S inside ``LEG_D0_REPRODUCTION_BAND``, so the field written to XDMF
  belongs to the gated S.

Negative controls (both asserted, no solve):

* (c1) the record comparator separates degrees: ``_rel(S, _record(2))`` exceeds the
  band;
* (c2) it sees a port relabelling: S with P1 and P2 swapped, against P1, exceeds
  the band (a cyclic rotation reads far smaller, by C4 symmetry; printed).

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:8 -n 2 -t 300

Not claimed: degree 2 (the Study default, held at ``-n 8`` by `WF-16`), tuning or
resonance, an absolute S11 / Z_in, a frequency sweep, anything `WF-16` does not
claim.
"""

from __future__ import annotations

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

from fem_em_solver import (  # noqa: E402
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    Study,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial  # noqa: E402
from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL, FREQ_ATOL_HZ  # noqa: E402
from tests.package.test_wf16_study_four_port import _record, _rel  # noqa: E402
from tests.validation.test_ans4_exact_strip_records import (  # noqa: E402
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
)
from tests.validation.test_lossy_sphere_fullwave import (  # noqa: E402
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_ops61_exact_strip_ports import (  # noqa: E402
    ANS4_EXACT_STRIP_CELL_RECORD,
)
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    _circulant_classes,
    _class_spread,
)
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    ADJACENT_SPREAD_BAND,
)
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_08_study_four_port_birdcage"
DEGREE = 1  # the record's degree (item 192); the Study default is 2.


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- (W) the record holds at its own width (OPS-41) --------------------
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 record is asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}; "
        f"got -n {comm.size}"
    )

    # ---- the Study ---------------------------------------------------------
    geometry = FourPortBirdcageGeometry.f_small()
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    drive = LumpedPortDrive()
    study = Study(geometry, phantom, drive, frequency_hz=FREQUENCY_64_HZ, degree=DEGREE)
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_study = time.perf_counter() - t0
    s = res.s_matrix
    z = res.z_matrix

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    msh = res.mesh
    local_tags = np.unique(res.cell_tags.values)
    all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist()) for t in part})
    named = {
        geometry.conductor_tag: f"copper-proxy conductor (sigma {geometry.conductor_sigma_s_per_m:g} S/m)",
        geometry.phantom_tag: "gel phantom",
    }
    air_tags = tuple(t for t in all_tags if t not in named)
    region_names = dict(named)
    for t in air_tags:
        region_names[t] = "air (PEC outer wall)"
    write_setup_figure(
        msh,
        res.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=air_tags,
        translucent_tags=(geometry.phantom_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:8 -- Study on the F-small four-port birdcage (WF-16, degree 1)",
        comm=comm,
    )

    # ---- (R) S against item 192's degree-1 record --------------------------
    p1 = _record(1)
    rel_r = _rel(s, p1)

    # ---- (G) PORT-11's gates on the Study's own S and Z --------------------
    recip = float(np.linalg.norm(s - s.T) / np.linalg.norm(s))
    sigma_sv = np.linalg.svd(s, compute_uv=False)
    sigma_max = float(np.max(sigma_sv))
    power = np.sum(np.abs(s) ** 2, axis=0)
    max_power = float(np.max(power))
    classes = _circulant_classes(z)
    spreads = {name: _class_spread(v) for name, v in classes.items()}
    pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
    worst = max(spreads.values())
    separation = pooled / worst if worst > 0.0 else np.inf

    # ---- (T) the Touchstone file, read back through scikit-rf --------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    s4p_path = OUTPUT_DIR / f"{BASENAME}.s4p"
    t0 = time.perf_counter()
    written = res.write_touchstone(s4p_path)
    t_ts = time.perf_counter() - t0
    readings = None
    if comm.rank == 0:
        import skrf

        ntw = skrf.Network(str(written))
        sk = np.asarray(ntw.s[0], dtype=np.complex128)
        readings = {
            "d_skrf": float(np.max(np.abs(sk - s))),
            "df": abs(float(ntw.f[0]) - FREQUENCY_64_HZ),
            "nports": int(ntw.nports),
            "z0": np.asarray(ntw.z0).real.ravel().tolist(),
            "rel_file": _rel(sk, p1),
        }
    readings = comm.bcast(readings, root=0)
    z0_expected = float(drive.port_impedance_ohm)

    # ---- negative controls (no solve) --------------------------------------
    c1 = _rel(s, _record(2))
    perm = [1, 0, 2, 3]
    c2 = _rel(s[np.ix_(perm, perm)], p1)
    rot = [1, 2, 3, 0]
    c_rot = _rel(s[np.ix_(rot, rot)], p1)

    # ---- (X) the P1-drive field, from the Study's own problem and mesh -----
    t0 = time.perf_counter()
    sweep = run_n_port_sparameter_sweep(
        res.problem,
        res.port_definitions,
        lumped_sheet_ports=res.specs,
        lumped_sheet_facet_tags=res.facet_tags,
        degree=res.degree,
        keep_fields=True,
    )
    t_sweep = time.perf_counter() - t0
    s_x = np.asarray(sweep.s_matrix, dtype=np.complex128)
    rel_x = _rel(s_x, s)
    fields = sweep.fields["P1"]
    e = fields.e_complex

    # Every collective on every rank.  The mesh is affine, so the N1curl Function
    # is interpolated, not written: |E| into P1, sigma into DG0.
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    e_abs = fem.Function(p1_space, name="E_abs")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1_space.element.interpolation_points)
    )
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()
    dg0 = fem.functionspace(msh, ("DG", 0))
    sigma_out = fem.Function(dg0, name="sigma")
    sigma_out.interpolate(fields.sigma_field)
    sigma_out.x.array[:] = np.real(sigma_out.x.array)
    sigma_out.x.scatter_forward()

    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        res.cell_tags,
        {"E_abs": e_abs, "sigma": sigma_out},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- prints (rank 0, no collective inside) -----------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-92 -- four-port S from one Study on the F-small birdcage (WF-16)")
        print("=" * 72)
        print(
            f"[study] degree {study.degree}, -n {comm.size}, f = {study.frequency_hz:.6e} Hz, "
            f"z0 = {drive.z0_ohm:g} ohm, terminated port impedance {z0_expected:g} ohm"
        )
        print(
            f"[mesh] {res.cells} cells (printed; version-tagged "
            f"{ANS4_EXACT_STRIP_CELL_RECORD}); Study {t_study:.2f} s (mesh "
            f"{res.mesh_time:.2f} s, sweep {res.sweep_time:.2f} s); Touchstone {t_ts:.3f} s; "
            f"keep_fields sweep {t_sweep:.2f} s"
        )
        print(f"[(W)] comm.size {comm.size} == record width {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}")
        print(
            f"[(R)] ||S - P1||/||P1|| {rel_r:.3e} (band LEG_D0_REPRODUCTION_BAND "
            f"{LEG_D0_REPRODUCTION_BAND:g})"
        )
        print(
            f"[(G)] reciprocity {recip:.3e} (band {RECIPROCITY_BAND:g}); sigma_max "
            f"{sigma_max:.12f}; max column power {max_power:.12f} (<= 1 + "
            f"{PASSIVITY_SIGMA_TOLERANCE:g}); spreads "
            + ", ".join(f"{k} {v * 100:.4f}%" for k, v in spreads.items())
            + f" (band {ADJACENT_SPREAD_BAND:g}); separation {separation:.2f}x "
            f"(floor {POOLED_SEPARATION_FLOOR:g})"
        )
        print(
            f"[(T)] max|S_skrf - S| {readings['d_skrf']:.3e} (ANCHOR_TOL {ANCHOR_TOL:g}); "
            f"|df| {readings['df']:.3e} Hz (FREQ_ATOL_HZ {FREQ_ATOL_HZ:g}); nports "
            f"{readings['nports']}; z0 {readings['z0']}; "
            f"||S_file - P1||/||P1|| {readings['rel_file']:.3e}"
        )
        print(
            f"[(c1)] ||S - degree-2 record||/||.|| {c1:.4e} "
            f"({c1 / LEG_D0_REPRODUCTION_BAND:.1e}x the band; must exceed it)"
        )
        print(
            f"[(c2)] P1<->P2 swap vs P1 {c2:.4f} (must exceed {LEG_D0_REPRODUCTION_BAND:g}); "
            f"printed: cyclic rotation {c_rot:.3e} -- C4 symmetry makes a rotation nearly "
            f"invisible, so the swap, not the rotation, is the control"
        )
        print(
            f"[(X)] keep_fields sweep S vs Study S {rel_x:.3e} (band "
            f"{LEG_D0_REPRODUCTION_BAND:g}; printed, not a record)"
        )
        print(f"[paraview] wrote {xdmf_path}; Touchstone {written}")

    assert readings["d_skrf"] <= ANCHOR_TOL, f"(T) |S_skrf - S| {readings['d_skrf']:.3e}"
    assert readings["df"] <= FREQ_ATOL_HZ, f"(T) frequency off by {readings['df']:.3e} Hz"
    assert readings["nports"] == 4, f"(T) nports {readings['nports']}"
    assert all(v == z0_expected for v in readings["z0"]), f"(T) z0 {readings['z0']}"
    assert readings["rel_file"] <= LEG_D0_REPRODUCTION_BAND, (
        f"(T) file vs P1 {readings['rel_file']:.3e}"
    )
    assert rel_r <= LEG_D0_REPRODUCTION_BAND, f"(R) rel {rel_r:.3e}"
    assert recip <= RECIPROCITY_BAND, f"(G) reciprocity {recip:.3e}"
    assert sigma_max <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, f"(G) sigma_max {sigma_max!r}"
    for v in power:
        assert v <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, f"(G) column power {v!r}"
    for name, v in spreads.items():
        assert v <= ADJACENT_SPREAD_BAND, f"(G) spread {name} {v:.3e}"
    assert separation >= POOLED_SEPARATION_FLOOR, f"(G) separation {separation:.2f}"
    assert c1 > LEG_D0_REPRODUCTION_BAND, f"(c1) blind: {c1:.3e}"
    assert c2 > LEG_D0_REPRODUCTION_BAND, f"(c2) blind: {c2:.3e}"
    assert rel_x <= LEG_D0_REPRODUCTION_BAND, f"(X) sweep S vs Study S {rel_x:.3e}"

    if comm.rank == 0:
        print(
            f"\nAll anchors and controls hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
