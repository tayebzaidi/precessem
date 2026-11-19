"""`POST-8` step 3a — the averaged Row 1 on `ans:2`'s named points, width-invariant.

PROJECT_PLAN §9 item 188; brief ``docs/planning/chunks/POST-8.md`` (§9 item
188 section).  Width is ``comm.size``; the module runs at ``-n 2`` (writes the
reference npz), then ``-n 4`` and ``-n 1`` (load it) on one commit.

Construction: item 177's ``case`` fixture (``test_post8_owner_rule_ans2.py``):
``build_four_port_sweep`` on the 0.0025 m 1 g rung, four single-drive
``_solve_driven`` solves P1…P4.  `ans:2`'s script is loaded by importlib only for
``_pointwise_sar``, ``_phantom_power`` and its constants — **never** ``main()``,
which writes tracked and private files.

Three reading families, each 4 × 4 (drive P1…P4 × point +x, +y, −x, −y):

* **raw** — the owner-rule SAR σ|E|²/(2ρ) through
  ``evaluate_vector_field_parallel(..., owner="lowest_original_cell")``;
* **cell** — the owner cell's average of σ|E|²/(2ρ), by two routes, both exact
  because |E|² is quadratic for degree-1 N1curl on an affine tet:
  (i) ``e_complex`` at the 4-point degree-2 rule on the owning rank, shared by
  buffer ``Allreduce`` SUM; (ii) ∫σ inner(E,E) χ/(2ρ) / ∫χ over a DG0 indicator
  χ of the owner cell at ``quadrature_degree`` 2, both integrals reduced;
* **1 g** — ``mass_averaged_sar`` centred on each point, radius
  ``averaging_ball_radius(ONE_GRAM_KG, ρ)``, ``QUADRATURE_DEGREE``, ρ from
  ``build_density_field`` as ``_build_mass_averaged`` builds it.

Asserted: (A0) cells and P1 phantom power; (A2) raw vs ``STEP1_N1CURL_SAR`` at
``A2_RTOL``; (R) the two cell routes agree ≤ 1e-12 and route (ii)'s
|Im| ≤ 1e-14·Re; (M) every 1 g ball inside the phantom; (W) at ``-n 4`` /
``-n 1`` owner indices and cell count equal the ``-n 2`` npz exactly and each
of the 48 values is within 1e-9 relative of its ``-n 2`` entry; (G) a planted
1 + 1e-8 perturbation per family is flagged by the (W) comparator, and only it.
Negative control (rule (e)): C1b asserted at ``-n 2`` / ``-n 4``; at ``-n 1``
printed only.  Printed, not asserted: cell/raw (predicted in [0.8, 1.25]) and
1 g/raw ratios, and the (W) per-family maxima beside their predictions
(raw / cell ≤ 5e-10, 1 g ≤ 1e-10).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.post.evaluation import (
    evaluate_vector_field_parallel,
    point_owner_cells,
)
from fem_em_solver.post.sar import (
    averaging_ball_radius,
    build_density_field,
    mass_averaged_sar,
)
from tests.complex_mode import complex_only
from tests.mesh.helpers import global_cell_tag_set
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT, PHANTOM_RADIUS
from tests.validation.test_birdcage_b1_plus_map import _solve_driven
from tests.validation.test_birdcage_b1_quadrature import QUADRATURE_STEP_DEG, _port_index
from tests.validation.test_birdcage_sar_1g_rung import (
    ONE_GRAM_RUNG_CELL_RECORD,
    PHANTOM_RESOLUTION_1G_RUNG,
)
from tests.validation.test_mass_averaged_sar_standard_masses import (
    ONE_GRAM_KG,
    QUADRATURE_DEGREE,
)
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_post8_owner_rule_ans2 import (
    A2_RTOL,
    ANS2_PATH,
    C1B_MIN_SEPARATION,
    PROBE_PATH,
    STEP1_N1CURL_SAR,
)

pytestmark = complex_only

#: (W) — the weekly's cross-width bound on each of the 48 values.
WIDTH_RTOL = 1.0e-9
#: (R) — the two exact cell-average routes; round-off only.
ROUTE_RTOL = 1.0e-12
#: (R) — route (ii)'s imaginary part relative to its real part.
ROUTE_IMAG_RTOL = 1.0e-14
#: (G) — the planted perturbation (10× the (W) bound).
PLANT = 1.0e-8
#: Printed predictions (not asserted).
PREDICTED_RAW_CELL = 5.0e-10
PREDICTED_ONE_GRAM = 1.0e-10
PREDICTED_CELL_RAW_RATIO = (0.8, 1.25)

#: The 4-point degree-2 tetrahedral rule (equal weights 1/4).
QUAD_A = 0.5854101966249685
QUAD_B = 0.1381966011250105

NPZ_DIR = Path("/workspace/logs/post8")
NPZ_PATH = NPZ_DIR / "step3a-n2.npz"
FAMILIES = ("raw", "cell", "one_gram")


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def case():
    """Item 177's ``case`` fixture, as built there."""
    comm = MPI.COMM_WORLD
    ans2 = _load(ANS2_PATH, "post8_step3a_ans2_module")
    probe = _load(PROBE_PATH, "post8_step3a_probe_module")
    assert not ans2.RESOLUTION_OVERRIDDEN

    sweep = build_four_port_sweep(phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG)
    azimuths = {f"P{s['tag'] - SHEET_IFACE}": float(s["azimuth_deg"]) for s in sweep["sheets"]}
    order = sorted(azimuths)
    indices = {
        pid: _port_index(azimuths[pid], azimuths["P1"], QUADRATURE_STEP_DEG) for pid in order
    }
    solves = {pid: _solve_driven(sweep, pid) for pid in order}
    built = {
        "mesh": sweep["mesh"],
        "cell_tags": sweep["cell_tags"],
        "cells": int(sweep["cells"]),
        "solves": solves,
    }
    r = ans2.CENTRE_RADIUS_M
    points = np.array(
        [[+r, 0.0, 0.0], [0.0, +r, 0.0], [-r, 0.0, 0.0], [0.0, -r, 0.0]],
        dtype=np.float64,
    )
    if comm.rank == 0:
        print(
            f"\n[POST-8 step 3a] -n {comm.size}: cells {built['cells']}; drives {order}; "
            f"port indices {indices}; points r = {r} m",
            flush=True,
        )
    return {"ans2": ans2, "probe": probe, "built": built, "points": points, "comm": comm}


def _sar_point(sigma, rho, e):
    return float(sigma) * float(np.sum(np.abs(e) ** 2)) / (2.0 * float(rho))


@pytest.fixture(scope="module")
def readings(case):
    comm, built, ans2 = case["comm"], case["built"], case["ans2"]
    points = case["points"]
    msh = built["mesh"]
    sigma = float(ans2.SALINE_SIGMA)
    rho = float(ans2.PHANTOM_RHO_KG_PER_M3)

    owner_original, owner_local = point_owner_cells(msh, points, comm)
    owner_original = np.asarray(owner_original, dtype=np.int64)

    # ---- raw (owner rule) ----
    raw = np.zeros((4, 4))
    for k in range(4):
        e = built["solves"][f"P{k + 1}"]["fields"].e_complex
        vals, valid = evaluate_vector_field_parallel(e, points, comm, owner="lowest_original_cell")
        assert bool(np.all(valid))
        vals = np.asarray(vals)
        for i in range(4):
            raw[k, i] = _sar_point(sigma, rho, vals[i])

    # ---- cell route (i): 4-point degree-2 rule on the owning rank ----
    gx = np.asarray(msh.geometry.x)
    gdofs = np.asarray(msh.geometry.dofmap)
    bary = np.full((4, 4), QUAD_B)
    np.fill_diagonal(bary, QUAD_A)
    cell_i_local = np.zeros((4, 4))
    for i in range(4):
        c = int(owner_local[i])
        if c < 0:
            continue
        verts = gx[gdofs[c]][:, :3]
        qpts = np.ascontiguousarray(bary @ verts)
        cells_arg = np.full(4, c, dtype=np.int32)
        for k in range(4):
            e = built["solves"][f"P{k + 1}"]["fields"].e_complex
            ev = np.asarray(e.eval(qpts, cells_arg)).reshape(4, 3)
            cell_i_local[k, i] = np.mean([_sar_point(sigma, rho, ev[q]) for q in range(4)])
    cell_i = np.zeros_like(cell_i_local)
    comm.Allreduce(cell_i_local, cell_i, op=MPI.SUM)

    # ---- cell route (ii): DG0 indicator, degree-2 quadrature, reduced ----
    V0 = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(V0)
    dx2 = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": 2})
    num_forms = [
        fem.form(
            sigma
            * ufl.inner(
                built["solves"][f"P{k + 1}"]["fields"].e_complex,
                built["solves"][f"P{k + 1}"]["fields"].e_complex,
            )
            * chi
            / (2.0 * rho)
            * dx2
        )
        for k in range(4)
    ]
    den_form = fem.form(chi * dx2)
    cell_ii = np.zeros((4, 4))
    cell_ii_imag = np.zeros((4, 4))
    for i in range(4):
        chi.x.array[:] = 0.0
        c = int(owner_local[i])
        if c >= 0:
            chi.x.array[V0.dofmap.cell_dofs(c)] = 1.0
        chi.x.scatter_forward()
        den = comm.allreduce(complex(fem.assemble_scalar(den_form)), op=MPI.SUM)
        assert den.real > 0.0
        for k in range(4):
            num = comm.allreduce(complex(fem.assemble_scalar(num_forms[k])), op=MPI.SUM)
            val = num / den
            cell_ii[k, i] = val.real
            cell_ii_imag[k, i] = val.imag

    # ---- 1 g ball on each point ----
    cell_tags = built["cell_tags"]
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh,
        rho,
        cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG},
    )
    a1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=rho)
    one_gram = np.zeros((4, 4))
    for k in range(4):
        fields = built["solves"][f"P{k + 1}"]["fields"]
        for i in range(4):
            res = mass_averaged_sar(
                fields.e_complex,
                sigma=fields.sigma_field,
                rho=rho_field,
                center=tuple(float(v) for v in points[i]),
                radius=a1g,
                comm=comm,
                quadrature_degree=QUADRATURE_DEGREE,
            )
            one_gram[k, i] = float(res["averaged_sar_w_per_kg"])

    # ---- C1b input: ans:2's CG1 route (default owner=None) ----
    pointwise = ans2._pointwise_sar(built, points, comm)
    cg1 = np.array([pointwise[k]["sar_w_per_kg"] for k in range(4)])

    out = {
        "raw": raw,
        "cell": cell_ii,
        "cell_quad": cell_i,
        "cell_imag": cell_ii_imag,
        "one_gram": one_gram,
        "owner_original": owner_original,
        "cells": int(built["cells"]),
        "a1g": a1g,
        "cg1": cg1,
        "rho": rho,
        "sigma": sigma,
    }

    if comm.rank == 0:
        print(f"\n[POST-8 3a] -n {comm.size}: owner original_cell_index {owner_original.tolist()}")
        print(f"[POST-8 3a] 1 g radius a1g = {a1g:.17e} m (rho {rho}, sigma {sigma})")
        for name, arr in (
            ("raw", raw),
            ("cell(i)", cell_i),
            ("cell(ii)", cell_ii),
            ("1g", one_gram),
            ("CG1", cg1),
        ):
            for k in range(4):
                print(f"  {name:8s} P{k + 1}: " + "  ".join(f"{v:.17e}" for v in arr[k]))
        print("[POST-8 3a] ratios [drive][point]: cell/raw | 1g/raw")
        for k in range(4):
            print(
                f"  P{k + 1}: "
                + "  ".join(f"{c / s:.6f}" for c, s in zip(cell_ii[k], raw[k]))
                + "  |  "
                + "  ".join(f"{g / s:.6f}" for g, s in zip(one_gram[k], raw[k]))
            )
        ratio = cell_ii / raw
        print(
            f"[POST-8 3a] cell/raw range [{ratio.min():.4f}, {ratio.max():.4f}] "
            f"(predicted within {PREDICTED_CELL_RAW_RATIO}, printed only)"
        )

    if comm.size == 2:
        if comm.rank == 0:
            NPZ_DIR.mkdir(parents=True, exist_ok=True)
            np.savez(
                NPZ_PATH,
                raw=raw,
                cell=cell_ii,
                one_gram=one_gram,
                owner_original=owner_original,
                cells=np.int64(built["cells"]),
            )
            print(f"[POST-8 3a] -n 2: wrote {NPZ_PATH}", flush=True)
        comm.Barrier()
    return out


def _flag(measured: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """The (W) comparator: True where |measured / reference − 1| > WIDTH_RTOL."""
    return np.abs(measured / reference - 1.0) > WIDTH_RTOL


def _load_reference():
    if not NPZ_PATH.exists():
        pytest.fail(f"{NPZ_PATH} missing: run the -n 2 window first (never recompute)")
    with np.load(NPZ_PATH) as data:
        return {k: np.array(data[k]) for k in data.files}


def test_a0_cells_and_power(case, readings):
    comm, built, ans2, probe = case["comm"], case["built"], case["ans2"], case["probe"]
    powers = ans2._phantom_power(built)
    p1 = powers[0]["dissipated_power_w"]
    rel = abs(p1 / probe.TRACKED_P1_POWER_W - 1.0)
    if comm.rank == 0:
        print(
            f"\n[POST-8 3a A0] -n {comm.size}: cells {built['cells']} vs "
            f"{ONE_GRAM_RUNG_CELL_RECORD}; P1 phantom power {p1:.17e} W rel {rel:.3e} "
            f"(<= {probe.POWER_RTOL:g})"
        )
    assert built["cells"] == ONE_GRAM_RUNG_CELL_RECORD
    assert rel <= probe.POWER_RTOL


def test_a2_raw_gate(case, readings):
    comm = case["comm"]
    rel = np.abs(readings["raw"] / STEP1_N1CURL_SAR - 1.0)
    if comm.rank == 0:
        print(f"\n[POST-8 3a A2] -n {comm.size}: max rel {rel.max():.3e} (<= {A2_RTOL:g})")
    assert float(rel.max()) <= A2_RTOL, rel


def test_r_cell_routes_agree(case, readings):
    comm = case["comm"]
    cell_ii, cell_i = readings["cell"], readings["cell_quad"]
    rel = np.abs(cell_i / cell_ii - 1.0)
    imag = np.abs(readings["cell_imag"]) / cell_ii
    if comm.rank == 0:
        print(
            f"\n[POST-8 3a R] -n {comm.size}: max |route(i)/route(ii) - 1| {rel.max():.3e} "
            f"(<= {ROUTE_RTOL:g}); max |Im|/Re route (ii) {imag.max():.3e} "
            f"(<= {ROUTE_IMAG_RTOL:g})"
        )
    assert float(rel.max()) <= ROUTE_RTOL, rel
    assert float(imag.max()) <= ROUTE_IMAG_RTOL, imag


def test_m_balls_inside_phantom(case, readings):
    comm, points = case["comm"], case["points"]
    a1g = readings["a1g"]
    radial = np.linalg.norm(points[:, :2], axis=1) + a1g
    if comm.rank == 0:
        print(
            f"\n[POST-8 3a M] -n {comm.size}: max |c| + a1g = {radial.max():.6e} m "
            f"(< PHANTOM_RADIUS {PHANTOM_RADIUS}); a1g {a1g:.6e} (< H/2 {PHANTOM_HEIGHT / 2})"
        )
    assert np.all(radial < PHANTOM_RADIUS)
    assert np.all(np.abs(points[:, 2]) + a1g < PHANTOM_HEIGHT / 2.0)
    assert a1g < PHANTOM_HEIGHT / 2.0


def test_w_width_invariance(case, readings):
    comm = case["comm"]
    width = comm.size
    if width == 2:
        if comm.rank == 0:
            print(f"\n[POST-8 3a W] -n 2 is the reference width ({NPZ_PATH} written)")
        return
    ref = _load_reference()
    assert np.array_equal(readings["owner_original"], ref["owner_original"])
    assert int(readings["cells"]) == int(ref["cells"])
    maxima = {}
    for fam in FAMILIES:
        rel = np.abs(readings[fam] / ref[fam] - 1.0)
        maxima[fam] = float(rel.max())
        if comm.rank == 0:
            pred = PREDICTED_ONE_GRAM if fam == "one_gram" else PREDICTED_RAW_CELL
            print(
                f"\n[POST-8 3a W] -n {width} vs -n 2, {fam}: max rel {rel.max():.3e} "
                f"(<= {WIDTH_RTOL:g}; predicted <= {pred:g}, printed)"
            )
            for k in range(4):
                print(f"  P{k + 1}: " + "  ".join(f"{v:.3e}" for v in rel[k]))
    for fam in FAMILIES:
        assert maxima[fam] <= WIDTH_RTOL, (fam, maxima[fam])


def test_g_planted_control(case, readings):
    comm = case["comm"]
    width = comm.size
    ref = {fam: readings[fam] for fam in FAMILIES} if width == 2 else _load_reference()
    plants = {"raw": (0, 0), "cell": (1, 2), "one_gram": (3, 1)}
    for fam in FAMILIES:
        planted = np.array(readings[fam], copy=True)
        planted[plants[fam]] *= 1.0 + PLANT
        flags = _flag(planted, ref[fam])
        expected = np.zeros((4, 4), dtype=bool)
        expected[plants[fam]] = True
        if comm.rank == 0:
            print(
                f"\n[POST-8 3a G] -n {width}, {fam}: planted {plants[fam]} x(1+{PLANT:g}); "
                f"flagged {np.argwhere(flags).tolist()}"
            )
        assert np.array_equal(flags, expected), (fam, np.argwhere(flags).tolist())


def test_c1b_negative_control(case, readings):
    comm, probe = case["comm"], case["probe"]
    width = comm.size
    row = readings["cg1"][0]
    seps = {}
    for w in (2, 4):
        ref = np.array(probe.LOGGED_P1_ROW[w])
        seps[w] = float(np.min(np.abs(row / ref - 1.0)))
    if comm.rank == 0:
        print(f"\n[POST-8 3a C1b] -n {width}: CG1 P1 row " + "  ".join(f"{v:.12e}" for v in row))
        for w in (2, 4):
            ref = np.array(probe.LOGGED_P1_ROW[w])
            print(
                f"  vs logged -n {w}: rel "
                + "  ".join(f"{a / b - 1.0:+.4e}" for a, b in zip(row, ref))
                + f"  min |rel| {seps[w]:.4f}"
            )
    if width in (2, 4):
        other = 4 if width == 2 else 2
        if comm.rank == 0:
            print(
                f"[POST-8 3a C1b] asserted: min |rel| vs -n {other} {seps[other]:.4f} "
                f"(>= {C1B_MIN_SEPARATION})"
            )
        assert seps[other] >= C1B_MIN_SEPARATION
    elif comm.rank == 0:
        print("[POST-8 3a C1b] -n 1: printed only (no logged -n 1 row)")
