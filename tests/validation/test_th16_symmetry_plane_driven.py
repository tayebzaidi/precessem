"""Validation: symmetry planes on the driven path, with port rescaling (TH-16 step 2).

One saline box, one interior lumped sheet (ĥ = ŷ, 50 Ω, 1 V), 64 MHz, degree 1.
The driven path pins only the ``pec_facet_tags`` faces; every other exterior
facet is natural, i.e. PMC (Jin §8.7.5). The full box is symmetric under
x → 0.48 − x (a plane containing ĥ: PMC, E_t even) and y → 0.40 − y (a plane
normal to ĥ: PEC, E_t odd), so on a **mirror-exact** mesh the half and quarter
models are the full discrete problem restricted to the even / odd subspace,
and the identities below hold to solver round-off:

(0) structure — the quarter Q = create_box([0, .24]×[0, .20]×[0, .48],
    (6, 5, 12)) mirrored in numpy (plane vertices reused, no crack);
(1) the full field is mirror-symmetric at 64 interior points;
(2) PMC half: Z_h = 2Z_f, V_h = V_f, I_h = I_f/2, E_h = E_f;
(3) PEC half: Z_h = Z_f/2, V_h = V_f/2, I_h = I_f;
(4) quarter: Z_q = Z_f; P = ½Re(VI*) is ½, ½, ¼ of P_f; Re Z_f > 0.

Band: `LEG_D0_REPRODUCTION_BAND` (1e-9), imported.

Negative controls: NC1 — the PMC half with its cut *pinned* (the wrong
condition), |Z/(2Z_f) − 1| predicted ≥ 10 % in W1 and asserted > 10 % in W2
under ``FEM_EM_TH16_STEP2_CONTROL=1``; NC2 — a directly meshed full box
(create_box (12, 10, 12), not mirror-symmetric) against the PMC half, printed
only (the asymmetric-mesh floor step 3 will face).

Not claimed: a mesh cutter, the birdcage, multi-port S, κ on cut sheets, any
memory figure, the eigen path.
"""

import os
import time

import basix.ufl
import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem
from fem_em_solver.ports.definitions import PortDefinition
from fem_em_solver.ports.lumped import LumpedSheetPortSpec, run_lumped_sheet_port_case
from fem_em_solver.ports.symmetry import full_impedance, halve_sheet_spec
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from tests.complex_mode import complex_only
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

pytestmark = complex_only

QUARTER_EXTENT = (0.24, 0.20, 0.48)
QUARTER_DIVISIONS = (6, 5, 12)
SPACING = 0.04
X_CUT = 0.24  # PMC plane (contains ĥ)
Y_CUT = 0.20  # PEC plane (normal to ĥ)
SHEET_Z = 0.24
SHEET_X = (0.16, 0.32)
SHEET_Y = (0.12, 0.28)
Z_PORT = 50.0
V_DRIVE = 1.0
TOL = 1e-8

# Face tags; the sheet tag is interior.
X0, X1, Y0, Y1, Z0, Z1 = 1, 2, 3, 4, 5, 6
SHEET_TAG = 10
ALL_FACES = (X0, X1, Y0, Y1, Z0, Z1)
NO_PMC_CUT = (X0, Y0, Y1, Z0, Z1)  # every face except x = 0.24

EXPECTED = {  # (vertices, cells, owned sheet facets)
    "quarter": (546, 2160, 8),
    "pmc_half": (1001, 4320, 16),
    "pec_half": (1014, 4320, 16),
    "full": (1859, 8640, 32),
}
EXPECTED_PMC_CUT_FACETS = {"quarter": 120, "pmc_half": 240}
AREA_RTOL = 1e-12
NC1_FLOOR = 0.10
CONTROL = os.environ.get("FEM_EM_TH16_STEP2_CONTROL") == "1"

SALINE = HomogeneousMaterial(sigma=SALINE_SIGMA, epsilon_r=SALINE_EPSILON_R, mu_r=1.0)

# 64 fixed interior points of the quarter, off every mesh plane.
_g = np.meshgrid(
    0.0137 + 0.0561 * np.arange(4),
    0.0113 + 0.0467 * np.arange(4),
    0.0291 + 0.1379 * np.arange(4),
    indexing="ij",
)
POINTS = np.ascontiguousarray(np.stack([a.ravel() for a in _g], axis=1))
FULL_X, FULL_Y = 2 * X_CUT, 2 * Y_CUT


def _mx(p):
    q = p.copy()
    q[:, 0] = FULL_X - q[:, 0]
    return q


def _my(p):
    q = p.copy()
    q[:, 1] = FULL_Y - q[:, 1]
    return q


# Reflection of a vector: E(mirror p) = S·E(p). PMC (x) even: flip E_x.
# PEC (y) odd: −R_y E = (−E_x, E_y, −E_z).
S_X = np.array([-1.0, 1.0, 1.0])
S_Y = np.array([-1.0, 1.0, -1.0])


def _quarter_arrays():
    q = dmesh.create_box(
        MPI.COMM_SELF,
        [np.zeros(3), np.array(QUARTER_EXTENT)],
        QUARTER_DIVISIONS,
        dmesh.CellType.tetrahedron,
    )
    x = np.asarray(q.geometry.x, dtype=np.float64)[:, :3].copy()
    cells = np.asarray(q.geometry.dofmap, dtype=np.int64).reshape(-1, 4).copy()
    return x, cells


def _mirror(x, cells, axis, plane):
    """Mirror (x, cells) across ``plane`` on ``axis``; plane vertices are reused."""
    xm = x.copy()
    xm[:, axis] = 2.0 * plane - xm[:, axis]
    both_x = np.vstack([x, xm])
    both_c = np.vstack([cells, cells[:, [1, 0, 2, 3]] + x.shape[0]])
    key = np.rint(both_x / SPACING).astype(np.int64)
    uniq, inverse = np.unique(key, axis=0, return_inverse=True)
    return uniq.astype(np.float64) * SPACING, inverse.reshape(-1)[both_c]


def _lattice(x, cells):
    key = np.rint(x / SPACING).astype(np.int64)
    uniq, inverse = np.unique(key, axis=0, return_inverse=True)
    return uniq.astype(np.float64) * SPACING, inverse.reshape(-1)[cells]


def _distribute(x, cells, comm):
    if comm.rank != 0:
        x = np.empty((0, 3), dtype=np.float64)
        cells = np.empty((0, 4), dtype=np.int64)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    partitioner = dmesh.create_cell_partitioner(dmesh.GhostMode.shared_facet, 2)
    return dmesh.create_mesh(
        comm,
        cells=np.ascontiguousarray(cells),
        x=np.ascontiguousarray(x),
        e=domain,
        partitioner=partitioner,
    )


def _tag(msh, extent):
    fdim = 2
    msh.topology.create_connectivity(fdim, 3)
    idx, val = [], []
    for tag, axis, v in (
        (X0, 0, 0.0),
        (X1, 0, extent[0]),
        (Y0, 1, 0.0),
        (Y1, 1, extent[1]),
        (Z0, 2, 0.0),
        (Z1, 2, extent[2]),
    ):
        f = dmesh.locate_entities_boundary(msh, fdim, lambda p, a=axis, c=v: np.abs(p[a] - c) < TOL)
        idx.append(f)
        val.append(np.full(f.size, tag, dtype=np.int32))
    sheet = dmesh.locate_entities(
        msh,
        fdim,
        lambda p: (np.abs(p[2] - SHEET_Z) < TOL)
        & (p[0] > SHEET_X[0] - TOL)
        & (p[0] < SHEET_X[1] + TOL)
        & (p[1] > SHEET_Y[0] - TOL)
        & (p[1] < SHEET_Y[1] + TOL),
    )
    idx.append(sheet)
    val.append(np.full(sheet.size, SHEET_TAG, dtype=np.int32))
    idx = np.concatenate(idx).astype(np.int32)
    val = np.concatenate(val).astype(np.int32)
    order = np.argsort(idx, kind="stable")
    idx, val = idx[order], val[order]
    return dmesh.meshtags(msh, fdim, idx, val)


def _owned_count(msh, tags, tag):
    n_local = msh.topology.index_map(2).size_local
    f = tags.find(tag)
    return msh.comm.allreduce(int(np.sum(f < n_local)), op=MPI.SUM)


def _sheet_area(msh, tags):
    dS = ufl.Measure("dS", domain=msh, subdomain_data=tags, subdomain_id=SHEET_TAG)
    one = fem.Constant(msh, default_scalar_type(1.0))
    local = fem.assemble_scalar(fem.form(one * dS))
    return float(np.real(msh.comm.allreduce(local, op=MPI.SUM)))


def _cut_plane_exterior(msh):
    """Owned exterior facets of ``msh`` lying on x = 0.24 or y = 0.20 (global)."""
    msh.topology.create_connectivity(2, 3)
    ext = dmesh.exterior_facet_indices(msh.topology)
    n_local = msh.topology.index_map(2).size_local
    ext = ext[ext < n_local]
    mid = dmesh.compute_midpoints(msh, 2, ext.astype(np.int32)) if ext.size else np.zeros((0, 3))
    on = (np.abs(mid[:, 0] - X_CUT) < TOL) | (np.abs(mid[:, 1] - Y_CUT) < TOL)
    return msh.comm.allreduce(int(np.sum(on)), op=MPI.SUM)


def _solve(msh, tags, spec, pec_tags, points, return_fields=False):
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_64_HZ,
        material=SALINE,
        boundary_condition="pec_zero_tangential_a",
        facet_tags=tags,
        pec_facet_tags=pec_tags,
    )
    port = PortDefinition(
        port_id="P1",
        positive_tag=SHEET_TAG,
        negative_tag=SHEET_TAG + 1,
        orientation="sheet_plus_y",
        z0_ohm=Z_PORT,
    )
    result, fields = run_lumped_sheet_port_case(
        problem,
        [port],
        [spec],
        facet_tags=tags,
        driven_port_id="P1",
        degree=1,
        verbose=False,
        return_fields=True,
    )
    r = result.responses["P1"]
    v, i = complex(r.voltage_v), complex(r.current_a)
    values, valid = evaluate_vector_field_parallel(fields.e_complex, points)
    out = {
        "V": v,
        "I": i,
        "Z": v / i,
        "P": 0.5 * (v * np.conj(i)).real,
        "E": values,
        "valid": valid,
    }
    if return_fields:
        out["fields"] = fields
    return out


def _say(comm, msg):
    if comm.rank == 0:
        print(f"\n[TH-16 step 2] {msg}", flush=True)


def build_cases():
    """The four mirror-exact meshes, their specs and pinned faces (EX-66 lift).

    Returns ``(cases, full_spec)``; ``cases[name] = (x, cells, extent, spec,
    pec_tags)`` for ``quarter``, ``pmc_half``, ``pec_half``, ``full``.
    """
    xq, cq = _quarter_arrays()
    xq, cq = _lattice(xq, cq)
    x_pmc, c_pmc = _mirror(xq, cq, 1, Y_CUT)  # cut stays x = 0.24
    x_pec, c_pec = _mirror(xq, cq, 0, X_CUT)  # cut stays y = 0.20
    x_full, c_full = _mirror(x_pmc, c_pmc, 0, X_CUT)

    full_spec = LumpedSheetPortSpec(
        port_id="P1",
        facet_tag=SHEET_TAG,
        port_impedance_ohm=Z_PORT + 0.0j,
        gap_height_m=SHEET_Y[1] - SHEET_Y[0],
        sheet_width_m=SHEET_X[1] - SHEET_X[0],
        drive_direction=(0.0, 1.0, 0.0),
        drive_voltage_v=V_DRIVE + 0.0j,
        interior=True,
    )
    pmc_spec = halve_sheet_spec(full_spec, "pmc", plane_normal=(1.0, 0.0, 0.0))
    pec_spec = halve_sheet_spec(full_spec, "pec", plane_normal=(0.0, 1.0, 0.0))
    q_spec = halve_sheet_spec(pmc_spec, "pec", plane_normal=(0.0, 1.0, 0.0))

    cases = {
        "quarter": (xq, cq, (X_CUT, Y_CUT, 0.48), q_spec, NO_PMC_CUT),
        "pmc_half": (x_pmc, c_pmc, (X_CUT, FULL_Y, 0.48), pmc_spec, NO_PMC_CUT),
        "pec_half": (x_pec, c_pec, (FULL_X, Y_CUT, 0.48), pec_spec, ALL_FACES),
        "full": (x_full, c_full, (FULL_X, FULL_Y, 0.48), full_spec, ALL_FACES),
    }
    return cases, full_spec


@pytest.fixture(scope="module")
def run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    cases, full_spec = build_cases()
    full_points = np.ascontiguousarray(
        np.vstack([POINTS, _mx(POINTS), _my(POINTS), _mx(_my(POINTS))])
    )
    out = {"structure": {}, "solve": {}}
    for name, (x, c, extent, spec, pec_tags) in cases.items():
        msh = _distribute(x, c, comm)
        tags = _tag(msh, extent)
        s = {
            "vertices": msh.topology.index_map(0).size_global,
            "cells": msh.topology.index_map(3).size_global,
            "sheet_facets": _owned_count(msh, tags, SHEET_TAG),
            "pmc_cut_facets": _owned_count(msh, tags, X1),
            "area": _sheet_area(msh, tags),
            "wh": spec.sheet_width_m * spec.gap_height_m,
        }
        if name == "full":
            s["cut_exterior"] = _cut_plane_exterior(msh)
        out["structure"][name] = s
        _say(comm, f"{name}: {s}")
        pts = full_points if name == "full" else POINTS
        out["solve"][name] = _solve(msh, tags, spec, pec_tags, pts)
        r = out["solve"][name]
        _say(
            comm,
            f"{name}: V={r['V']:.12e} I={r['I']:.12e} Z={r['Z']:.12e} "
            f"P={r['P']:.12e} valid={int(np.sum(r['valid']))}/{len(pts)}",
        )
        if name == "pmc_half":
            # NC1: same mesh and spec, the x = 0.24 cut pinned (wrong condition).
            out["nc1"] = _solve(msh, tags, spec, ALL_FACES, POINTS)

    # NC2: a directly meshed (not mirror-symmetric) full box.
    box = dmesh.create_box(
        comm,
        [np.zeros(3), np.array([FULL_X, FULL_Y, 0.48])],
        (12, 10, 12),
        dmesh.CellType.tetrahedron,
        ghost_mode=dmesh.GhostMode.shared_facet,
    )
    btags = _tag(box, (FULL_X, FULL_Y, 0.48))
    out["nc2"] = _solve(box, btags, full_spec, ALL_FACES, POINTS)

    zf = out["solve"]["full"]["Z"]
    out["nc1_metric"] = abs(out["nc1"]["Z"] / (2.0 * zf) - 1.0)
    out["nc2_metric"] = abs(out["solve"]["pmc_half"]["Z"] / (2.0 * out["nc2"]["Z"]) - 1.0)
    _say(
        comm,
        f"NC1 (PMC half, cut pinned) |Z/(2Z_f) - 1| = {out['nc1_metric']:.6e} "
        f"(Z = {out['nc1']['Z']:.12e}); predicted >= {NC1_FLOOR:.0%}: "
        f"{'HOLDS' if out['nc1_metric'] >= NC1_FLOOR else 'FAILS'}; "
        f"asserted: {CONTROL}",
    )
    _say(
        comm,
        f"NC2 (create_box (12,10,12) full vs PMC half) |Z_h/(2Z_box) - 1| = "
        f"{out['nc2_metric']:.6e} (Z_box = {out['nc2']['Z']:.12e}) — printed only",
    )
    out["elapsed"] = time.perf_counter() - t0
    _say(comm, f"fixture elapsed {out['elapsed']:.1f} s")
    return out


def _rel(a, b):
    return abs(a / b - 1.0)


def test_structure(run):
    """(0) Mirror construction: counts exact, no crack, sheet area w·h."""
    st = run["structure"]
    for name, (nv, nc, ns) in EXPECTED.items():
        assert st[name]["vertices"] == nv, name
        assert st[name]["cells"] == nc, name
        assert st[name]["sheet_facets"] == ns, name
        assert abs(st[name]["area"] / st[name]["wh"] - 1.0) <= AREA_RTOL, name
    for name, n in EXPECTED_PMC_CUT_FACETS.items():
        assert st[name]["pmc_cut_facets"] == n, name
    assert st["full"]["cut_exterior"] == 0


def test_full_field_mirror_symmetric(run):
    """(1) The full field is mirror-symmetric at 64 interior points."""
    r = run["solve"]["full"]
    assert r["valid"].all()
    e = r["E"]
    n = POINTS.shape[0]
    e0, ex, ey, exy = e[:n], e[n : 2 * n], e[2 * n : 3 * n], e[3 * n :]
    scale = np.max(np.linalg.norm(e0, axis=1))
    resid = (
        max(
            np.max(np.linalg.norm(ex - S_X * e0, axis=1)),
            np.max(np.linalg.norm(ey - S_Y * e0, axis=1)),
            np.max(np.linalg.norm(exy - S_X * S_Y * e0, axis=1)),
        )
        / scale
    )
    _say(MPI.COMM_WORLD, f"(1) mirror residual {resid:.3e} (max|E| = {scale:.6e})")
    assert resid <= LEG_D0_REPRODUCTION_BAND


def test_pmc_half(run):
    """(2) PMC half: Z_h = 2Z_f, V_h = V_f, I_h = I_f/2, E_h = E_f, P_h = P_f/2."""
    f, h = run["solve"]["full"], run["solve"]["pmc_half"]
    n = POINTS.shape[0]
    assert h["valid"].all()
    e_rel = np.max(np.linalg.norm(h["E"] - f["E"][:n], axis=1)) / np.max(
        np.linalg.norm(f["E"][:n], axis=1)
    )
    ratios = {
        "Z": _rel(full_impedance(h["Z"], "pmc"), f["Z"]),
        "V": _rel(h["V"], f["V"]),
        "I": _rel(h["I"], f["I"] / 2.0),
        "P": _rel(h["P"], f["P"] / 2.0),
        "E": e_rel,
    }
    _say(MPI.COMM_WORLD, f"(2) PMC half rel errors {ratios}")
    for k, v in ratios.items():
        assert v <= LEG_D0_REPRODUCTION_BAND, k


def test_pec_half(run):
    """(3) PEC half: Z_h = Z_f/2, V_h = V_f/2, I_h = I_f, P_h = P_f/2."""
    f, h = run["solve"]["full"], run["solve"]["pec_half"]
    ratios = {
        "Z": _rel(full_impedance(h["Z"], "pec"), f["Z"]),
        "V": _rel(h["V"], f["V"] / 2.0),
        "I": _rel(h["I"], f["I"]),
        "P": _rel(h["P"], f["P"] / 2.0),
    }
    _say(MPI.COMM_WORLD, f"(3) PEC half rel errors {ratios}")
    for k, v in ratios.items():
        assert v <= LEG_D0_REPRODUCTION_BAND, k


def test_quarter(run):
    """(4) Quarter: Z_q = Z_f, P_q = P_f/4, Re Z_f > 0."""
    f, q = run["solve"]["full"], run["solve"]["quarter"]
    ratios = {"Z": _rel(q["Z"], f["Z"]), "P": _rel(q["P"], f["P"] / 4.0)}
    _say(MPI.COMM_WORLD, f"(4) quarter rel errors {ratios}; Re Z_f = {f['Z'].real:.9e}")
    for k, v in ratios.items():
        assert v <= LEG_D0_REPRODUCTION_BAND, k
    assert f["Z"].real > 0.0


def test_nc1_pinned_pmc_cut(run):
    """NC1: pinning the PMC cut breaks Z_h = 2Z_f. Asserted only under the flag (W2)."""
    if not CONTROL:
        pytest.skip(
            f"NC1 predicted in W1 (printed {run['nc1_metric']:.6e}); asserted in W2 "
            "under FEM_EM_TH16_STEP2_CONTROL=1"
        )
    assert run["nc1_metric"] > NC1_FLOOR
