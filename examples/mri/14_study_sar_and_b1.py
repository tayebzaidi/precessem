"""Example (`EX-106`): 10 g SAR and |B1+| from one ``Study`` -- one field-keeping solve,
two results, on the F-small four-port birdcage.

`WF-16` steps 4 and 5 put ``Study.sar(mass_g=10.0)`` and ``Study.b1_plus()`` /
``Study.b1_plus_at()`` on the exact-rectangle four-port birdcage (κ off, operator
ruling 2026-10-02) at 64 MHz.  This example builds the same ``Study`` at degree 1 on
``-n 2`` -- the width and degree of item 192's registered record -- asks it for
``|B1+|`` first (that call runs the one field-keeping solve and caches it on
``study.fields_result``), then for the 10 g mass-averaged SAR on the four single
drives, and shows that the second result reuses the first one's solve.

Every band, constant and helper is imported from the gate modules (the ANS-1 rule);
the example states no σ, εr, frequency, radius, band or tolerance of its own.  The
gate modules export no SAR or |B1+| record constant, so the record digits live in the
example-local ``REC`` dict below, each with its gate-log ``file:line``, and are
asserted at the printed precision (the printed digit rounds to the recorded digit).

Asserted:

* (W) ``comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH`` before anything else;
* (R) ``_rel(S, P1) <= LEG_D0_REPRODUCTION_BAND`` (S of the field-keeping solve);
* (One) after the SAR calls, ``study.fields_result is fr`` and ``study.result is None``;
* (ii) the 10 g C4 quadrant identity at ``C4_COVARIANCE_BAND``;
* (i) the whole-phantom ball against ``mean_sar`` at ``EXACT_IDENTITY_RTOL``;
* (iv) containment of the 10 g ball and the whole-phantom ball;
* (B1) the |B1+| C4 covariance identity at ``C4_COVARIANCE_BAND``;
* (B3) every sample point valid on every image;
* (T) ``PHANTOM_CELL_TAG == study.geometry.phantom_tag``;
* every printed record digit equal to its cited gate-log digit (``REC``).

Negative controls (both asserted in sign, no solve; the size is printed only):

* (c1) the mis-paired far-side ball against the right one, above the band at every k;
* (cB) P2 read at the counter-rotated points against P1, above the band.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:14 -n 2 -t 300

Not claimed: an absolute SAR, a C95.3 / IEC compliance figure or limit, a spatial-peak
search, 1 g SAR, |B1+| homogeneity or CV, an absolute or closed-form |B1+|, a
quadrature drive, tuning or resonance, degree 2, h-convergence, copper, human scale,
anything `WF-16` does not claim.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

import dolfinx
from dolfinx import default_scalar_type, fem

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
from fem_em_solver.post.sar import averaging_ball_radius, mean_sar  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.study import SAR_QUADRATURE_DEGREE  # noqa: E402
from fem_em_solver.utils.instrumentation import report_peak_rss  # noqa: E402

from tests.mesh.test_birdcage_implant_rod import _cell_values_with_ghosts  # noqa: E402
from tests.package.test_wf16_study_four_port import _record, _rel  # noqa: E402
from tests.validation.test_ans4_exact_strip_records import (  # noqa: E402
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
)
from tests.validation.test_birdcage_b1_plus_map import (  # noqa: E402
    C4_COVARIANCE_BAND,
    PHANTOM_CELL_TAG,
    _relative_l2,
    _rotate_z,
    _sample_points,
)
from tests.validation.test_birdcage_b1_quadrature import (  # noqa: E402
    QUADRATURE_STEP_DEG,
    _port_index,
)
from tests.validation.test_birdcage_sar_mass_averaged import (  # noqa: E402
    CENTRE_RADIUS_M,
    CONTROL_CEILING,
    CONTROL_PREDICTED,
    EXACT_IDENTITY_RTOL,
    WHOLE_PHANTOM_BALL_RADIUS_M,
)
from tests.validation.test_lossy_sphere_fullwave import (  # noqa: E402
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_mass_averaged_sar_standard_masses import TEN_GRAM_KG  # noqa: E402
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    LEG_D0_REPRODUCTION_BAND,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_14_study_sar_and_b1"
DEGREE = 1  # the record's degree (item 192); the Study default is 2.
SAR = "averaged_sar_w_per_kg"
FIGURE_BALL_TAGS = (901, 902, 903, 904)  # figure-only; never passed to a solve.

_S4 = "docs/testing/logs/20261008T184215Z_WF-16-s4-W1-deg1.log"
_S5 = "docs/testing/logs/20261008T185141Z_WF-16-s5-W1-deg1.log"

# Record digits of the WF-16 step 4 / 5 degree-1 gate windows (-n 2).  Each entry is
# asserted at the printed precision by ``_Rec.check`` -- never re-recorded from here.
REC = {
    # _S4:1881
    "cells": 117206,
    "phantom_cells": 537,
    "a10_mm": 13.3650,
    "a10_over_h": 0.8910,
    # _S4:66 (iv)
    "reach_mm": 28.3650,
    # _S4:1893-1896 (ii) diagonal, W/kg and kg
    "sar_diag": [1.980612231e-05, 1.963377184e-05, 2.032653726e-05, 2.073991516e-05],
    "mass_diag": [1.000012094e-02, 1.001130889e-02, 1.000069359e-02, 9.934035589e-03],
    "mass_kernel": 1.000000000e-02,
    # _S4:1897 cyclic pairs, percent
    "pairs_pct": [0.8702, 3.5284, 2.0337, 4.5024],
    # _S4:1898 max over four evaluated centres (P4)
    "sar_max": 2.073991516e-05,
    # _S4:1888 (i): ball power W, mass kg (equal to the tagged mean_sar values)
    "whole_power_w": 1.585977541981e-06,
    "whole_mass_kg": 2.195294683148e-01,
    # _S4:1906 (c1) percent
    "c1_pct": [88.1162, 87.9366, 88.6600, 88.9315],
    # _S5:1877
    "delta_deg": 90.000000,
    "n_points": 51,
    # _S5:1880 (B3) valid points per image
    "b3_counts": 51,
    # _S5:1884 (B1) percent: P2@+D, P4@-D, P3@2D
    "b1_pct": [1.8695, 1.9494, 1.6405],
    # _S5:1885 |B1+| at the phantom centre, T: P1..P4
    "centre_t": [1.667281487e-08, 1.700738851e-08, 1.665270745e-08, 1.674331925e-08],
    # _S5:1892 (cB) percent
    "cb_pct": 35.0040,
}


class _Rec:
    """Collects print-precision comparisons against ``REC``; asserted once at the end."""

    def __init__(self) -> None:
        self.rows: list[tuple[str, str, str, bool]] = []

    def check(self, label: str, value, rec, fmt: str) -> bool:
        got, want = format(value, fmt), format(rec, fmt)
        ok = got == want
        self.rows.append((label, got, want, ok))
        return ok

    def failures(self) -> list[tuple[str, str, str, bool]]:
        return [r for r in self.rows if not r[3]]


def _global_count(comm, local: int) -> int:
    return int(comm.allreduce(int(local), op=MPI.SUM))


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    rec = _Rec()

    # ---- (W) the record holds at its own width (OPS-41) --------------------
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 record is asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}; "
        f"got -n {comm.size}"
    )
    if np.dtype(default_scalar_type).kind != "c":
        raise RuntimeError(
            "mri:14 needs the complex DolfinX build (source /usr/local/bin/dolfinx-complex-mode)"
        )

    # ---- the Study ---------------------------------------------------------
    geometry = FourPortBirdcageGeometry.f_small()
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    study = Study(
        geometry, phantom, LumpedPortDrive(), frequency_hz=FREQUENCY_64_HZ, degree=DEGREE
    )
    assert PHANTOM_CELL_TAG == geometry.phantom_tag, "(T) sample-set tag is not the Study's"

    # ---- step 1: ask for |B1+| -- this runs the one field-keeping solve -----
    t0 = time.perf_counter()
    maps = study.b1_plus(comm=comm)
    t_first = time.perf_counter() - t0
    fr = study.fields_result
    assert fr is not None and set(maps) == {"P1", "P2", "P3", "P4"}

    # ---- step 2: sheets, azimuths, drive indices, centres ------------------
    sheets = {f"P{s['tag'] - geometry.sheet_iface}": s for s in fr.sheets}
    az = {pid: float(s["azimuth_deg"]) for pid, s in sheets.items()}
    delta_deg = (az["P2"] - az["P1"]) % 360.0
    delta = np.radians(delta_deg)
    idx = {pid: _port_index(az[pid], az["P1"], QUADRATURE_STEP_DEG) for pid in az}
    by_k = {k: pid for pid, k in idx.items()}
    assert sorted(by_k) == [0, 1, 2, 3]
    centres = {
        k: (
            CENTRE_RADIUS_M * float(np.cos(np.radians(az[by_k[k]]))),
            CENTRE_RADIUS_M * float(np.sin(np.radians(az[by_k[k]]))),
            0.0,
        )
        for k in range(4)
    }
    rho = study.density_kg_m3
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=rho)

    # ---- the setup figure: figure-only ball tags, never passed to a solve ---
    msh = fr.mesh
    tdim = msh.topology.dim
    n_owned = int(msh.topology.index_map(tdim).size_local)
    n_total = n_owned + int(msh.topology.index_map(tdim).num_ghosts)
    local_tags = np.unique(fr.cell_tags.values)
    all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist()) for t in part})
    assert not set(FIGURE_BALL_TAGS) & set(all_tags), (
        f"figure-only tags {FIGURE_BALL_TAGS} collide with the mesh tags {all_tags}"
    )
    mids_all = dolfinx.mesh.compute_midpoints(msh, tdim, np.arange(n_total, dtype=np.int32))
    ball_k = np.full(n_total, -1, dtype=np.int64)
    for k in range(4):
        c = np.asarray(centres[k], dtype=np.float64)
        d2 = np.sum((mids_all - c) ** 2, axis=1)
        ball_k[(d2 <= a10**2) & (ball_k < 0)] = k
    ball_owned_counts = [
        _global_count(comm, int(np.count_nonzero(ball_k[:n_owned] == k))) for k in range(4)
    ]

    fig_values = _cell_values_with_ghosts(msh, fr.cell_tags)[:n_owned].astype(np.int32).copy()
    for k in range(4):
        fig_values[ball_k[:n_owned] == k] = FIGURE_BALL_TAGS[k]
    fig_tags = dolfinx.mesh.meshtags(
        msh, tdim, np.arange(n_owned, dtype=np.int32), fig_values
    )
    named = {
        geometry.conductor_tag: f"copper-proxy conductor (sigma {geometry.conductor_sigma_s_per_m:g} S/m)",
        geometry.phantom_tag: "gel phantom",
    }
    air_tags = tuple(t for t in all_tags if t not in named)
    region_names = dict(named)
    for t in air_tags:
        region_names[t] = "air (PEC outer wall)"
    for k in range(4):
        region_names[FIGURE_BALL_TAGS[k]] = f"10 g ball at {by_k[k]} centre"
    write_setup_figure(
        msh,
        fig_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=air_tags,
        translucent_tags=(geometry.phantom_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:14 -- 10 g balls in the gel phantom, F-small four-port birdcage (WF-16)",
        comm=comm,
    )

    # ---- step 3: |B1+| at the sample set, five images, and at the origin ----
    points = _sample_points({"mesh": fr.mesh, "cell_tags": fr.cell_tags})
    images = {
        "P1@0": ("P1", points),
        "P2@+D": ("P2", _rotate_z(points, delta)),
        "P4@-D": ("P4", _rotate_z(points, -delta)),
        "P3@2D": ("P3", _rotate_z(points, 2.0 * delta)),
        "P2@-D (control)": ("P2", _rotate_z(points, -delta)),
    }
    t0 = time.perf_counter()
    read = {
        label: study.b1_plus_at(pts, drives=[pid], comm=comm)[pid]
        for label, (pid, pts) in images.items()
    }
    centre = {
        pid: float(study.b1_plus_at(np.zeros((1, 3)), drives=[pid], comm=comm)[pid][0][0])
        for pid in ("P1", "P2", "P3", "P4")
    }
    t_eval = time.perf_counter() - t0

    # ---- step 4: 10 g SAR through the same Study (9 operator calls) ---------
    t0 = time.perf_counter()
    diag = {
        k: study.sar(mass_g=10.0, centres=[centres[k]], drives=[by_k[k]], comm=comm)[by_k[k]][0]
        for k in range(4)
    }
    control = {
        k: study.sar(
            mass_g=10.0, centres=[centres[(k + 2) % 4]], drives=[by_k[k]], comm=comm
        )[by_k[k]][0]
        for k in range(4)
    }
    whole = study.sar(
        radius_m=WHOLE_PHANTOM_BALL_RADIUS_M,
        centres=[(0.0, 0.0, 0.0)],
        drives=["P1"],
        comm=comm,
    )["P1"][0]
    t_ops = time.perf_counter() - t0
    same_solve = study.fields_result is fr
    result_none = study.result is None

    # ---- step 5: tagged mean_sar for (i) ------------------------------------
    f1 = fr.fields["P1"]
    tagged = mean_sar(
        f1.e_complex,
        sigma=f1.sigma_field,
        rho=rho,
        cell_tags=fr.cell_tags,
        subdomain_ids=geometry.phantom_tag,
        comm=comm,
    )

    # ---- (R) ----------------------------------------------------------------
    rel_r = _rel(fr.s_matrix, _record(1))

    # ---- (iv) containment ---------------------------------------------------
    circum = float(np.hypot(geometry.phantom_radius, 0.5 * geometry.phantom_height))
    nearest_conductor = min(
        geometry.ring_radius - 0.5 * geometry.leg_width,
        geometry.ring_radius - geometry.ring_minor_radius,
    )
    reach = CENTRE_RADIUS_M + a10
    contained = reach < geometry.phantom_radius
    whole_contained = circum < WHOLE_PHANTOM_BALL_RADIUS_M < nearest_conductor

    # ---- (i) whole-phantom coverage ----------------------------------------
    d_p = abs(whole["dissipated_power_w"] / tagged["dissipated_power_w"] - 1.0)
    d_m = abs(whole["mass_kg"] / (rho * tagged["volume_m3"]) - 1.0)

    # ---- (ii) the 10 g C4 quadrant identity ---------------------------------
    pairs = {k: abs(diag[(k + 1) % 4][SAR] - diag[k][SAR]) / diag[k][SAR] for k in range(4)}
    kernel_mass = rho * 4.0 / 3.0 * np.pi * a10**3
    k_max = max(range(4), key=lambda k: diag[k][SAR])

    # ---- (c1) mis-paired far-side ball --------------------------------------
    miss_c1 = {k: abs(control[k][SAR] - diag[k][SAR]) / diag[k][SAR] for k in range(4)}

    # ---- (B1) / (B3) / (cB) -------------------------------------------------
    n_pts = int(points.shape[0])
    counts = {label: int(v.sum()) for label, (_, v) in read.items()}
    mask = np.logical_and.reduce([v for _, v in read.values()])
    ref = read["P1@0"][0]
    b1_table = {
        label: _relative_l2(vals, ref, mask)
        for label, (vals, _) in read.items()
        if label not in ("P1@0", "P2@-D (control)")
    }
    miss_cb = _relative_l2(read["P2@-D (control)"][0], ref, mask)

    # ---- fields for ParaView (every collective on every rank) ---------------
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    dg0 = fem.functionspace(msh, ("DG", 0))
    out_fields = {}
    for pid in ("P1", "P2", "P3", "P4"):
        v = maps[pid].x.array.reshape(-1, 3)
        fn = fem.Function(p1_space, name=f"B1plus_abs_{pid}")
        assert fn.x.array.shape[0] == v.shape[0], "CG1 vector and P1 scalar node layouts differ"
        fn.x.array[:] = np.abs(v[:, 0] + 1j * v[:, 1]) / 2.0
        fn.x.scatter_forward()
        out_fields[f"B1plus_abs_{pid}"] = fn
    e = f1.e_complex
    e_abs = fem.Function(p1_space, name="E_abs_P1")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1_space.element.interpolation_points)
    )
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()
    out_fields["E_abs_P1"] = e_abs
    pdens = fem.Function(dg0, name="power_density_P1")
    pdens.interpolate(
        fem.Expression(
            0.5 * ufl.real(f1.sigma_field * ufl.inner(e, e)), dg0.element.interpolation_points
        )
    )
    pdens.x.array[:] = np.real(pdens.x.array)
    pdens.x.scatter_forward()
    out_fields["power_density_P1"] = pdens
    sigma_out = fem.Function(dg0, name="sigma")
    sigma_out.interpolate(f1.sigma_field)
    sigma_out.x.array[:] = np.real(sigma_out.x.array)
    sigma_out.x.scatter_forward()
    out_fields["sigma"] = sigma_out
    mask_fn = fem.Function(dg0, name="ball10g_mask")
    mask_fn.x.array[:] = 0.0
    dofs = dg0.dofmap.list.reshape(-1)
    inside = ball_k >= 0
    mask_fn.x.array[dofs[inside]] = (ball_k[inside] + 1).astype(np.float64)
    mask_fn.x.scatter_forward()
    out_fields["ball10g_mask"] = mask_fn

    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, fr.cell_tags, out_fields, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    phantom_cells = _global_count(
        comm,
        int(np.count_nonzero(fr.cell_tags.values[fr.cell_tags.indices < n_owned] == geometry.phantom_tag)),
    )
    peak = report_peak_rss(comm, label="[mri:14] summed peak ru_maxrss")

    # ---- record digits against REC (print precision) ------------------------
    rec.check("cells", fr.cells, REC["cells"], "d")
    rec.check("phantom cells", phantom_cells, REC["phantom_cells"], "d")
    rec.check("a_10g mm", a10 * 1e3, REC["a10_mm"], ".4f")
    rec.check("a_10g / h", a10 / geometry.resolution, REC["a10_over_h"], ".4f")
    rec.check("r0 + a_10g mm", reach * 1e3, REC["reach_mm"], ".4f")
    for k in range(4):
        rec.check(f"SAR_10g k={k}", diag[k][SAR], REC["sar_diag"][k], ".9e")
        rec.check(f"mass k={k}", diag[k]["mass_kg"], REC["mass_diag"][k], ".9e")
        rec.check(f"pair k={k}", pairs[k] * 100, REC["pairs_pct"][k], ".4f")
        rec.check(f"c1 k={k}", miss_c1[k] * 100, REC["c1_pct"][k], ".4f")
    rec.check("kernel mass", kernel_mass, REC["mass_kernel"], ".9e")
    rec.check("SAR max", diag[k_max][SAR], REC["sar_max"], ".9e")
    rec.check("ball power", whole["dissipated_power_w"], REC["whole_power_w"], ".12e")
    rec.check("tagged power", tagged["dissipated_power_w"], REC["whole_power_w"], ".12e")
    rec.check("ball mass", whole["mass_kg"], REC["whole_mass_kg"], ".12e")
    rec.check("rho*V", rho * tagged["volume_m3"], REC["whole_mass_kg"], ".12e")
    rec.check("drive rotation deg", delta_deg, REC["delta_deg"], ".6f")
    rec.check("sample points", n_pts, REC["n_points"], "d")
    for label, c in counts.items():
        rec.check(f"B3 {label}", c, REC["b3_counts"], "d")
    for j, label in enumerate(b1_table):
        rec.check(f"B1 {label}", b1_table[label] * 100, REC["b1_pct"][j], ".4f")
    for j, pid in enumerate(("P1", "P2", "P3", "P4")):
        rec.check(f"|B1+| centre {pid}", centre[pid], REC["centre_t"][j], ".9e")
    rec.check("cB", miss_cb * 100, REC["cb_pct"], ".4f")

    # ---- prints (rank 0, no collective inside) ------------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-106 -- 10 g SAR and |B1+| from one Study on the F-small birdcage (WF-16)")
        print("=" * 72)
        print(
            f"[study] degree {study.degree}, -n {comm.size}, f = {study.frequency_hz:.6e} Hz, "
            f"quadrature_degree {SAR_QUADRATURE_DEGREE}"
        )
        print(
            f"[mesh] {fr.cells} cells (printed; version-tagged {REC['cells']}), "
            f"{phantom_cells} phantom (tag {geometry.phantom_tag}) cells; resolution "
            f"{geometry.resolution} m, a_10g = {a10 * 1e3:.4f} mm, a_10g / h = "
            f"{a10 / geometry.resolution:.4f}"
        )
        print(
            f"[time] first call b1_plus() (the one field-keeping solve + 4 CG1 projections) "
            f"{t_first:.2f} s (mesh {fr.mesh_time:.2f} s, sweep {fr.sweep_time:.2f} s); "
            f"|B1+| evaluations {t_eval:.2f} s; 9 operator calls {t_ops:.2f} s"
        )
        print(f"[(W)] comm.size {comm.size} == record width {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}")
        print(
            f"[(R)] ||S_keep - P1||/||P1|| {rel_r:.3e} (band LEG_D0_REPRODUCTION_BAND "
            f"{LEG_D0_REPRODUCTION_BAND:g}; record 9.557e-11 at {_S4}:1884)"
        )
        print(
            f"[(One)] study.fields_result is fr: {same_solve}; study.result is None: "
            f"{result_none} (solve() never called; the SAR calls reused the b1_plus() solve)"
        )
        print(
            f"[(iv)] r0 + a_10g = {reach * 1e3:.4f} mm < {geometry.phantom_radius * 1e3:.1f} mm; "
            f"circumradius {circum:.4f} m < {WHOLE_PHANTOM_BALL_RADIUS_M} m < conductor "
            f"{nearest_conductor:.4f} m (gate {_S4}:66)"
        )
        print(
            f"[(i)] ball power {whole['dissipated_power_w']:.12e} W, tagged "
            f"{tagged['dissipated_power_w']:.12e} W: rel {d_p:.3e}; mass {whole['mass_kg']:.12e} kg "
            f"vs rho*V {rho * tagged['volume_m3']:.12e} kg: rel {d_m:.3e} "
            f"(ASSERTED <= {EXACT_IDENTITY_RTOL:.0e}; gate {_S4}:1888)"
        )
        lines = "\n".join(
            f"    k={k} ({by_k[k]} at c_{k}): SAR_10g = {diag[k][SAR]:.9e} W/kg, mass "
            f"{diag[k]['mass_kg']:.9e} kg vs rho*4/3*pi*a^3 {kernel_mass:.9e} kg"
            for k in range(4)
        )
        print(
            f"[(ii)] 10 g diagonal (gate {_S4}:1893-1896):\n{lines}\n"
            "    cyclic pairs |SAR(c_k+1;k+1) - SAR(c_k;k)| / SAR(c_k;k): "
            + ", ".join(f"{k}->{(k + 1) % 4} {pairs[k] * 100:.4f}%" for k in range(4))
            + f" (ASSERTED <= {C4_COVARIANCE_BAND * 100:.1f}%; gate {_S4}:1897)\n"
            f"    10 g max over 4 evaluated centres at 1 V drive, NOT a spatial search, "
            f"NOT a C95.3 compliance figure: {diag[k_max][SAR]:.9e} W/kg ({by_k[k_max]})"
        )
        print(
            "[(c1)] mis-paired far-side ball |SAR(c_k+2;k) - SAR(c_k;k)| / SAR(c_k;k): "
            + ", ".join(f"k={k} {miss_c1[k] * 100:.4f}%" for k in range(4))
            + f" (sign ASSERTED > {C4_COVARIANCE_BAND * 100:.1f}%; size PREDICTED "
            f"≈ {CONTROL_PREDICTED * 100:.0f}%, ceiling {CONTROL_CEILING * 100:.0f}%, "
            f"not asserted; gate {_S4}:1906)"
        )
        print(
            f"[(B3)] valid points per image (of {n_pts}): {counts} (ASSERTED full; "
            f"drive rotation {delta_deg:.6f} deg; gate {_S5}:1877,1880)"
        )
        print(
            "[(B1)] ||B1+(image) - B1+(P1)||_2 / ||B1+(P1)||_2: "
            + ", ".join(f"{k} {v * 100:.4f}%" for k, v in b1_table.items())
            + f" (ASSERTED <= {C4_COVARIANCE_BAND * 100:.1f}%; gate {_S5}:1884)\n"
            "    |B1+| at the phantom centre per 1 V drive (a record, not absolute, not "
            "homogeneity): "
            + ", ".join(f"{p} {v:.9e} T" for p, v in centre.items())
        )
        print(
            f"[(cB)] control P2@-D vs P1: {miss_cb * 100:.4f}% (sign ASSERTED > "
            f"{C4_COVARIANCE_BAND * 100:.1f}%; size printed, not asserted; gate {_S5}:1892)"
        )
        print(
            f"[ball mask] owned cells per 10 g ball (picture only, not the operator's "
            f"quadrature ball): "
            + ", ".join(f"{by_k[k]} {ball_owned_counts[k]}" for k in range(4))
        )
        print(f"[mem] summed peak {peak / 2**30:.4f} GiB")
        bad = rec.failures()
        print(
            f"[REC] {len(rec.rows)} record digits compared at print precision; "
            f"{len(bad)} differ" + "".join(f"\n    DIFF {r[0]}: got {r[1]}, record {r[2]}" for r in bad)
        )
        print(f"[paraview] wrote {xdmf_path}")

    # ---- assertions ---------------------------------------------------------
    assert rel_r <= LEG_D0_REPRODUCTION_BAND, f"(R) rel {rel_r:.3e}"
    assert same_solve, "(One) the SAR calls re-solved: study.fields_result changed"
    assert result_none, "(One) study.result was written: solve() ran"
    assert contained, f"(iv) r0 + a_10g {reach:.5f} m reaches the phantom wall"
    assert whole_contained, "(iv) whole-phantom ball not between circumradius and conductor"
    assert d_p <= EXACT_IDENTITY_RTOL, f"(i) power rel {d_p:.3e}"
    assert d_m <= EXACT_IDENTITY_RTOL, f"(i) mass rel {d_m:.3e}"
    for k in range(4):
        assert pairs[k] <= C4_COVARIANCE_BAND, f"(ii) pair {k}->{(k + 1) % 4} {pairs[k]:.4%}"
    for k in range(4):
        assert miss_c1[k] > C4_COVARIANCE_BAND, f"(c1) blind at k={k}: {miss_c1[k]:.3e}"
    assert n_pts > 0
    for label, c in counts.items():
        assert c == n_pts, f"(B3) {label}: {c} of {n_pts}"
    for label, v in b1_table.items():
        assert v <= C4_COVARIANCE_BAND, f"(B1) {label} {v:.4%}"
    assert miss_cb > C4_COVARIANCE_BAND, f"(cB) blind: {miss_cb:.3e}"
    bad = rec.failures()
    assert not bad, f"record digits differ from the gate logs: {bad}"

    if comm.rank == 0:
        print(
            f"\nAll anchors and controls hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
