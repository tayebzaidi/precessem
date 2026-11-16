"""GEO-38 C6a (§9 item 667, §10 chain C step 6): the ASTM-box-loaded
F-clinical-1.5T module — I36b's arm R (the rod fixture with tag 403 mapped to
the gel) as the box-loaded coil, build-only.

Brief: ``docs/planning/chunks/GEO-38.md`` "§9 item 667".

Selected by ``FEM_EM_GEO38_C6`` only (never by ``-k``):

* unset — the gated test skips; the coding tests run (no gmsh);
* ``build`` — I36b's route (`_route_drive`) into C4a's fixture: (A0) / census /
  (U1) / (D-rod) asserted inside, then (A0m) on the gel twin and the lattice's
  containment on the real mesh (tag 3 at every point), RSS printed;
* ``1`` — `_run_box` (§9 item 694, ``xl`` entry 52; composed, never run in a
  slot): arm R alone, one factorisation and a held column on P33 / P41 at the
  node, r = 0.2 Ω on all 64 gaps, then `_box_readouts` under ccw = (1, −1j) with
  (1, +1j) printed beside it at c\\*.  Gates (imported bands): (A0) / census /
  (U1) inside the fixture, (A0m), (Z), (R1), (P) + NEG on both columns,
  P_surf > 0, (P_RC), coverage (i) power and mass, the lattice in tag 3.  No
  mirror gate, no class gate;
* anything else — ``ValueError``.

I36b's node knob ``FEM_EM_GEO36_S8_DRIVE_C_PF`` is required under ``1`` and
refused under ``build`` (I36b's `_node_knob`, reused unchanged).

C34 (§9 item 820): ``FEM_EM_GEO38_C6_DEGREE`` — unset ⇒ degree 1, ``2`` ⇒ 2,
else ``ValueError``; refused under ``build``; routed into C4a's fixture and arm
R's solver as I36b's ``FEM_EM_GEO36_S8_DRIVE_DEGREE`` is, the JSON gaining a
``-d2`` suffix.  The foreign guard refuses I36's knobs, so ``xxl`` entry 62 is
two pytest calls in one window (K-62).

**The lattice.**  18 named points, §9 item 693 (Review ruling 2026-10-28 18:00,
provisional): x offsets {−Lx/4, 0, +Lx/4} from the box centre (Lx = 0.42 m) ×
``C6A_LATTICE_DZ_M`` on z (Lz = 0.65 m), at mid-depth y.

**The gel is tags (3, 403) together** (Review ruling 2026-10-30 18:00): in arm R
tag 403 carries gel σ, so the density and the tagged mean count it with tag 3;
the coverage ball carries no conductor limit (its radius passes the rungs, and
on a Leontovich hole mesh off the gel σ = ρ = 0).  K6's control is that defect.

Scope: the module and its coding tests.  Not claimed: any solve, |B₁⁺| or SAR on
the coil, tuning or compliance, the weekly's choice of sense or lattice.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

C6_ENV = "FEM_EM_GEO38_C6"
#: C34 (§9 item 820): the degree knob of ``1``, read only by `test_geo38_c6a_box_loaded`.
C6_DEGREE_ENV = "FEM_EM_GEO38_C6_DEGREE"
#: The ruled lattice size: §9 item 693 (Review ruling 2026-10-28 18:00), docs/planning/chunks/GEO-38.md.
BRIEF_LATTICE_N = 18
#: Six named z offsets from the box centre, none at z = 0 (GEO-38.md "Review ruling 2026-10-28 18:00").
C6A_LATTICE_DZ_M = (-0.25, -0.15, -0.05, 0.05, 0.15, 0.25)
#: Variables a C6a window must never read (I35 / I36 / C4 / step-1b / stop-rule knobs).
FOREIGN_ENVS = (
    "FEM_EM_GEO36_S8_DRIVE",
    "FEM_EM_GEO36_S8_DRIVE_DEGREE",
    "FEM_EM_GEO36_S8_COLUMN",
    "FEM_EM_GEO36_ROD_BOX",
    "FEM_EM_GEO38_STEP1B",
    "FEM_EM_GEO38_C4",
    "FEM_EM_GEO38_C4_NODES_PF",
    "FEM_EM_GEO38_C4_DEGREE",
    "FEM_EM_GMSH_STOP_RULE",
)
# Predictions written first (brief), printed only.
PREDICTED_WINDOW_S = (470.0, 520.0)
#: entry 52's predicted window under ``1``, printed only.
PREDICTED_BOX_WINDOW_S = (850.0, 1000.0)
#: The JSON of ``1`` (gitignored logs), one per node.
C6_JSON_DIR = "/workspace/logs/geo38"


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-38 C6a] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Knobs.
# ---------------------------------------------------------------------------
def _c6_mode():
    raw = os.environ.get(C6_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("build", "1"):
        return raw
    raise ValueError(f"{C6_ENV} must be unset, 'build' or '1', got {raw!r}")


def _knobs():
    """Foreign guard, then ``(mode, c_pf)`` with I36b's node-knob semantics."""
    from tests.validation import test_geo36_rod_f_clinical_drive as i36

    bad = [k for k in FOREIGN_ENVS if os.environ.get(k) is not None]
    if bad:
        raise ValueError(f"C6a never reads these variables; unset them: {bad}")
    mode = _c6_mode()
    return mode, i36._node_knob(mode)


def _c6_degree(mode):
    """C34 (§9 item 820): `_degree_knob_c4`'s semantics under ``C6_DEGREE_ENV`` —
    unset / empty ⇒ ``None`` (degree 1), ``"2"`` ⇒ 2, anything else ⇒ ValueError —
    refused under ``build`` (I36b's message: the box mesh's (Deg-2) is already
    asserted, `20261027T170559Z_GEO-36-I36b-W1-build.log:11989–11990`).  Separate
    from `_knobs`, so K1's equalities are unchanged."""
    raw = os.environ.get(C6_DEGREE_ENV, "").strip()
    if raw == "":
        degree = None
    elif raw == "2":
        degree = 2
    else:
        raise ValueError(f"{C6_DEGREE_ENV}={raw!r}: only unset/empty (degree 1) or '2' are accepted")
    if mode == "build" and degree is not None:
        raise ValueError(f"{C6_DEGREE_ENV} is refused under {C6_ENV}=build: "
                         f"the build proof covers both degrees (Deg-2)")
    return degree


# ---------------------------------------------------------------------------
# The box, the rod and the lattice (pure python).
# ---------------------------------------------------------------------------
def _box():
    """``(centre, half)`` of the ASTM box from the generator's own convention
    (`io/mesh.py` `phantom_box`: ``[-Lx/2, Lx/2] × [-Ly/2 + off, Ly/2 + off] ×
    [-Lz/2, Lz/2]``) and `_box_kwargs()` (`test_geo38_f_clinical_box.py:54–65`)."""
    from tests.validation.test_geo38_f_clinical_box import _box_kwargs

    box = _box_kwargs()["phantom_box"]
    lx, ly, lz = (float(v) for v in box["size_m"])
    return np.array([0.0, float(box["table_offset_m"]), 0.0]), 0.5 * np.array([lx, ly, lz])


def _lattice(heights=C6A_LATTICE_DZ_M):
    """The ruled lattice (§9 item 693): x offsets {−Lx/4, 0, +Lx/4} from the box
    centre × ``C6A_LATTICE_DZ_M`` on z, at mid-depth y.  ``heights`` admits K3's
    printed attempt-1 control.  ``[(label, (x, y, z))]``."""
    centre, half = _box()
    quarter = 0.5 * half[0]  # Lx / 4
    out = []
    for dx in (-quarter, 0.0, quarter):
        for dz in heights:
            p = centre + np.array([dx, 0.0, dz])
            out.append((f"x{dx * 1e3:+.0f}mm dz{dz * 1e3:+.0f}mm", tuple(float(v) for v in p)))
    return out


def _wall_containment(c, a):
    """WF-10's `_containment` (`test_wf10_f_human_tuned_sar10g.py:95`) with the box's
    six walls in place of R / H: ``<= 0`` ⇒ the ball of radius ``a`` is ≥ 1 mm clear
    of every wall."""
    from tests.validation.test_wf7_f_human_sar10g_hotspot import CONTAINMENT_MARGIN_M

    centre, half = _box()
    d = np.abs(np.asarray(c, dtype=float) - centre) + a - (half - CONTAINMENT_MARGIN_M)
    return float(d.max())


def _rod_axis_distance(c):
    """Distance from ``c`` to the rod's axis segment (centre ± ROD_LENGTH/2 on ẑ);
    the rod (tag 403) is the capsule of radius ROD_RADIUS about that segment."""
    from tests.mesh.test_birdcage_implant_rod import ROD_LENGTH
    from tests.validation import test_geo36_rod_f_clinical as i35

    rc = np.asarray(i35._rod_centre(), dtype=float)
    d = np.asarray(c, dtype=float) - rc
    dz = max(0.0, abs(d[2]) - 0.5 * ROD_LENGTH)
    return float(np.sqrt(d[0] ** 2 + d[1] ** 2 + dz ** 2))


def _a10():
    from fem_em_solver.post.sar import averaging_ball_radius
    from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_mass_averaged_sar_standard_masses import TEN_GRAM_KG

    return averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)


def _lattice_report(heights=C6A_LATTICE_DZ_M):
    """``[(label, point, wall, rod_dist, in_gel)]``; ``in_gel`` ⇔ the 10 g ball is
    ≥ 1 mm clear of the rod capsule (rod distance ≥ ROD_RADIUS + a10 + margin) and
    inside the box."""
    from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
    from tests.validation.test_wf7_f_human_sar10g_hotspot import CONTAINMENT_MARGIN_M

    a = _a10()
    rows = []
    for label, c in _lattice(heights):
        wall = _wall_containment(c, a)
        rd = _rod_axis_distance(c)
        rows.append((label, c, wall, rd,
                     bool(rd >= ROD_RADIUS + a + CONTAINMENT_MARGIN_M and wall <= a)))
    return rows


# ---------------------------------------------------------------------------
# C6b (§9 item 819, §10 C6b): the box-wide 10 g reader, knob-gated OFF by default.
# ---------------------------------------------------------------------------
#: Unset (or empty) ⇒ off; ``1`` ⇒ on; anything else ⇒ ValueError.  Read inside
#: `_run_box`, so C21 inherits it; ``xl`` entries 58 / 59 leave it unset (Review
#: ruling 2026-11-05 03:00, premise 4).
C6B_GRID_ENV = "FEM_EM_GEO38_C6_GRID"
#: The z and y steps.  The x step is Lx/8 (52.5 mm on the real box), not 50 mm: a
#: 50 mm grid anchored at the box centre cannot contain the lattice's x = ±Lx/4 =
#: ±105 mm (Review ruling 2026-11-05 03:00, premise 3; control (G-x)).
C6B_GRID_DZ_M = 0.05
C6B_GRID_DY_M = 0.05
#: Grid points are matched to lattice points by coordinates, to this tolerance (m).
C6B_MATCH_M = 1.0e-12


def _grid_knob():
    raw = os.environ.get(C6B_GRID_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{C6B_GRID_ENV} must be unset or '1', got {raw!r}")


def _grid_indices(half, step):
    """Indices n with |n·step| strictly inside ``half`` (1e-9 relative slack, so
    0.21 / 0.0525 = 4 gives |n| ≤ 3)."""
    n = int(np.ceil(half / step - 1.0e-9)) - 1
    return range(-n, n + 1)


def _grid_report(*, x_step=None, y_indices=None):
    """Every grid candidate as ``[(label, point, wall, rod_dist, kept)]``
    (`_lattice_report`'s shape): centres = box centre + (i·Lx/8, j·0.05, k·0.05),
    centres strictly inside the box.  ``kept`` ⇔ wall ≤ 0 (strict form) and the
    lattice's own rod clearance (rod distance ≥ ROD_RADIUS + a10 + margin).
    `_box()` is read at call time (K8's stub patches it).  ``x_step`` /
    ``y_indices`` admit controls (G-x) / (G-c)."""
    from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
    from tests.validation.test_wf7_f_human_sar10g_hotspot import CONTAINMENT_MARGIN_M

    centre, half = _box()
    dx = float(2.0 * half[0] / 8.0) if x_step is None else float(x_step)
    js = _grid_indices(half[1], C6B_GRID_DY_M) if y_indices is None else y_indices
    a = _a10()
    bar = ROD_RADIUS + a + CONTAINMENT_MARGIN_M
    rows = []
    for i in _grid_indices(half[0], dx):
        for j in js:
            for k in _grid_indices(half[2], C6B_GRID_DZ_M):
                off = np.array([i * dx, j * C6B_GRID_DY_M, k * C6B_GRID_DZ_M])
                c = tuple(float(v) for v in centre + off)
                wall = _wall_containment(c, a)
                rd = _rod_axis_distance(c)
                rows.append((f"x{off[0] * 1e3:+.1f}mm y{off[1] * 1e3:+.1f}mm dz{off[2] * 1e3:+.1f}mm",
                             c, wall, rd, bool(wall <= 0.0 and rd >= bar)))
    return rows


def _match_lattice(grid_points, lattice_points):
    """``[lattice index or None]`` per grid point, by coordinates (≤ C6B_MATCH_M)."""
    lat = np.asarray(lattice_points, dtype=float).reshape(-1, 3)
    out = []
    for p in np.asarray(grid_points, dtype=float).reshape(-1, 3):
        d = np.max(np.abs(lat - p), axis=1) if lat.size else np.zeros(0)
        hit = np.flatnonzero(d <= C6B_MATCH_M)
        out.append(int(hit[0]) if hit.size else None)
    return out


# ---------------------------------------------------------------------------
# The build branch (real mesh).
# ---------------------------------------------------------------------------
def _point_tags(msh, cell_tags, points, comm):
    """Global set of cell tags owning each point (rank-safe: gathered)."""
    from dolfinx import geometry

    from tests.mesh.test_birdcage_implant_rod import _cell_values_with_ghosts

    tagv = _cell_values_with_ghosts(msh, cell_tags)
    pts = np.ascontiguousarray(np.asarray(points, dtype=np.float64))
    tree = geometry.bb_tree(msh, msh.topology.dim)
    cand = geometry.compute_collisions_points(tree, pts)
    hit = geometry.compute_colliding_cells(msh, cand, pts)
    local = [set(int(tagv[c]) for c in hit.links(i)) for i in range(len(pts))]
    gathered = comm.allgather(local)
    return [set().union(*(g[i] for g in gathered)) for i in range(len(pts))]


def _run_build(comm, monkeypatch):
    import dataclasses

    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo36_rod_f_clinical_drive as i36
    from tests.validation import test_geo38_f_clinical_column as col
    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    t_window = time.perf_counter()
    _say(comm, f"\nstart: -n {comm.size}, BUILD-ONLY; I36b's route, arm R map; PREDICTED window "
               f"{PREDICTED_WINDOW_S[0]:.0f}-{PREDICTED_WINDOW_S[1]:.0f} s (printed)")
    i36._route_drive(monkeypatch)
    monkeypatch.delenv(c4.C4_DEGREE_ENV, raising=False)
    fixture = c4._build_c4_fixture(comm)  # (A0) / census / (U1) / (D-rod) / (F) asserted inside
    ctx = fixture["built"]["ctx"]
    msh = ctx["msh"]
    problem_m = ctx["solver"].problem
    problem_r = dataclasses.replace(
        problem_m, material_map=i36._gel_twin_map(dict(problem_m.material_map)))
    census_m = col._owned_tag_counts(msh, problem_m.cell_tags, comm)
    census_r = col._owned_tag_counts(msh, problem_r.cell_tags, comm)
    a0m = i36._a0m_failures(census_m, census_r)
    rows = _lattice_report()
    tags = _point_tags(msh, ctx["cell_tags"], [r[1] for r in rows], comm)
    for (label, c, wall, rd, _g), t in zip(rows, tags):
        _say(comm, f"lattice {label} {c}: wall {wall:+.4e} m, rod-axis {rd:.4e} m, tags {sorted(t)}")
    bad = [r[0] for r, t in zip(rows, tags) if r[2] > 0.0 or t != {i35.PHANTOM_CELL_TAG}]
    _m, per, tot = _rss_gib(comm)
    _say(comm, f"(A0m) {a0m or 'none'}; lattice n {len(rows)} (brief {BRIEF_LATTICE_N}), failures "
               f"{bad or 'none'}; RSS per rank {', '.join(f'{r:.2f}' for r in (per or []))} GiB, "
               f"summed {tot:.3f} GiB; window {time.perf_counter() - t_window:.1f} s")
    assert not a0m, a0m
    assert len(rows) == BRIEF_LATTICE_N, len(rows)
    assert not bad, bad


# ---------------------------------------------------------------------------
# The ``1`` branch (§9 item 694): the readouts as a function of fields only, so
# K6 runs them on a small mesh, then the body.
# ---------------------------------------------------------------------------
def _coverage_radius():
    """The box's half-diagonal + ``COVERAGE_BALL_PAD_M`` from `_box()`."""
    from tests.validation.test_wf7_f_human_sar10g_hotspot import COVERAGE_BALL_PAD_M

    _centre, half = _box()
    return float(np.linalg.norm(half) + COVERAGE_BALL_PAD_M)


def _gel_density(msh, cell_tags):
    """ρ_gel on tags (3, 403), 0 on every other tag (ruling 2026-10-30 18:00, defect 3)."""
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from fem_em_solver.post.sar import build_density_field
    from tests.mesh.helpers import global_cell_tag_set
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3

    return build_density_field(
        msh, PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(global_cell_tag_set(msh, cell_tags))
                     if int(t) not in (i35.PHANTOM_CELL_TAG, BIRDCAGE_IMPLANT_CELL_TAG)})


def _coverage_identity(e, sigma_f, rho_field, cell_tags, ids, centre, radius, comm):
    """Coverage (i), `test_wf10_f_human_tuned_sar10g.py:418–426`: the ball on
    ``rho_field`` against ``mean_sar`` on ``ids``.  ``(whole, tagged, power_miss,
    mass_miss)``, reduced."""
    from fem_em_solver.post.sar import mass_averaged_sar, mean_sar
    from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_mass_averaged_sar_standard_masses import QUADRATURE_DEGREE

    whole = mass_averaged_sar(e, sigma=sigma_f, rho=rho_field,
                              center=tuple(float(v) for v in centre), radius=radius,
                              comm=comm, quadrature_degree=QUADRATURE_DEGREE)
    tagged = mean_sar(e, sigma=sigma_f, rho=PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
                      subdomain_ids=ids, comm=comm)
    power_miss = abs(whole["dissipated_power_w"] / tagged["dissipated_power_w"] - 1.0)
    mass_miss = abs(whole["mass_kg"] / (PHANTOM_RHO_KG_PER_M3 * tagged["volume_m3"]) - 1.0)
    return whole, tagged, float(power_miss), float(mass_miss)


def _box_sample_points(msh, cell_tags, comm):
    """`_sample_points`'s gather / lexsort / stride (`test_birdcage_b1_plus_map.py:246`)
    on tag-3 centroids, no cylinder cut (its ``SAMPLE_RADIUS_M`` is F-small's);
    identical on every rank."""
    import dolfinx

    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation.test_birdcage_b1_plus_map import MAX_SAMPLE_POINTS

    tdim = msh.topology.dim
    owned = int(msh.topology.index_map(tdim).size_local)
    tagged = cell_tags.find(i35.PHANTOM_CELL_TAG)
    local_cells = np.asarray(tagged[tagged < owned], dtype=np.int32)
    mids = (dolfinx.mesh.compute_midpoints(msh, tdim, local_cells) if local_cells.size
            else np.zeros((0, 3), dtype=np.float64))
    gathered = comm.allgather(np.ascontiguousarray(mids, dtype=np.float64))
    points = np.concatenate(gathered, axis=0)
    points = points[np.lexsort((points[:, 2], points[:, 1], points[:, 0]))]
    if points.shape[0] > MAX_SAMPLE_POINTS:
        points = points[::int(np.ceil(points.shape[0] / MAX_SAMPLE_POINTS))]
    return points


def _box_readouts(msh, cell_tags, e_ccw, e_cw, sigma_f, omega, rho_field, rows, sample_points, comm,
                  *, coverage=None, grid_rows=None):
    """All of C6's post-processing, numbers only.  ``rows`` = `_lattice_report`'s shape
    (label and point read).  ``coverage`` = ``(centre, radius)``, default the box's
    (keyword-only: K6 passes its own mesh's half-diagonal + 1 mm)."""
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from fem_em_solver.post.faraday import magnetic_flux_density_from_e, project_to_cg1
    from fem_em_solver.post.sar import mass_averaged_sar
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation.test_birdcage_b1_quadrature import _cv, _read_senses
    from tests.validation.test_mass_averaged_sar_standard_masses import QUADRATURE_DEGREE

    assert e_ccw.function_space.mesh is msh and e_cw.function_space.mesh is msh
    labels = [r[0] for r in rows]
    lattice = np.ascontiguousarray(np.array([r[1] for r in rows], dtype=np.float64))
    senses = {}
    for sense, e in (("ccw", e_ccw), ("cw", e_cw)):
        b_cg1 = project_to_cg1(magnetic_flux_density_from_e(e, omega, name=f"B_{sense}"),
                               name=f"B_{sense}_cg1")
        plus, minus, valid = _read_senses(b_cg1, sample_points)
        plus_l, minus_l, valid_l = _read_senses(b_cg1, lattice)
        mean_plus = float(np.mean(plus[valid]))
        mean_minus = float(np.mean(minus[valid]))
        by_z = {}
        for z, v in zip(lattice[:, 2], plus_l):
            by_z.setdefault(float(z), []).append(float(v))
        z_means = np.array([np.mean(by_z[z]) for z in sorted(by_z)])
        senses[sense] = {
            "b_plus": plus, "b_minus": minus, "valid": valid, "n_valid": int(np.sum(valid)),
            "mean_b_plus": mean_plus, "mean_b_minus": mean_minus,
            "sense_ratio": mean_plus / mean_minus, "cv": _cv(plus[valid]),
            "lattice_b_plus": plus_l, "lattice_b_minus": minus_l,
            "lattice_valid": bool(np.all(valid_l)),
            "lattice_z_spread": float(np.ptp(z_means) / np.mean(z_means)),
        }

    def _ball(e, c):
        return mass_averaged_sar(e, sigma=sigma_f, rho=rho_field,
                                 center=tuple(float(v) for v in c), radius=_a10(), comm=comm,
                                 quadrature_degree=QUADRATURE_DEGREE)

    balls = [_ball(e_ccw, c) for c in lattice]
    sar = np.array([b["averaged_sar_w_per_kg"] for b in balls])
    star = int(np.argmax(sar))
    c_star = lattice[star]
    cw_star = _ball(e_cw, c_star)
    if coverage is None:
        coverage = (_box()[0], _coverage_radius())
    whole, tagged, power_miss, mass_miss = _coverage_identity(
        e_ccw, sigma_f, rho_field, cell_tags, (i35.PHANTOM_CELL_TAG, BIRDCAGE_IMPLANT_CELL_TAG),
        coverage[0], coverage[1], comm)
    result = {
        "senses": senses, "labels": labels, "sar10g": sar,
        "ball_mass_kg": np.array([b["mass_kg"] for b in balls]),
        "star": star, "star_label": labels[star], "c_star": tuple(float(v) for v in c_star),
        "sar_star": float(sar[star]), "sar_star_cw": float(cw_star["averaged_sar_w_per_kg"]),
        "whole": whole, "tagged": tagged, "power_miss": power_miss, "mass_miss": mass_miss,
    }
    if grid_rows is None:
        return result
    # C6b: the lattice values reused at their grid indices (matched by coordinates),
    # one ccw ball per new point, one cw ball at g*.
    t_grid = time.perf_counter()
    g_labels = [r[0] for r in grid_rows]
    g_points = np.ascontiguousarray(np.array([r[1] for r in grid_rows], dtype=np.float64))
    match = _match_lattice(g_points, lattice)
    # §9 item 844: each new ball's wall time (each centre is a cold JIT compile,
    # `post/sar.py:257`), rank-max reduced; printed, never asserted.
    ball_t = []

    def _timed_ball(p):
        t0 = time.perf_counter()
        v = float(_ball(e_ccw, p)["averaged_sar_w_per_kg"])
        ball_t.append(time.perf_counter() - t0)
        return v

    g_sar = np.array([float(sar[mk]) if mk is not None else _timed_ball(p)
                      for p, mk in zip(g_points, match)])
    if ball_t:
        bt = np.empty(len(ball_t))
        comm.Allreduce(np.array(ball_t, dtype=np.float64), bt, op=MPI.MAX)
        if comm.rank == 0:
            print(f"(C6b) ball times n {len(bt)} first {bt[0]:.2f} median {np.median(bt):.2f} "
                  f"max {np.max(bt):.2f} sum {np.sum(bt):.2f}", flush=True)
    g_star = int(np.argmax(g_sar))
    g_cw = _ball(e_cw, g_points[g_star])
    t_grid = comm.allreduce(time.perf_counter() - t_grid, op=MPI.MAX)
    result["grid"] = {
        "n": len(grid_rows), "labels": g_labels, "points": g_points, "sar10g": g_sar,
        "star": g_star, "c_star": tuple(float(v) for v in g_points[g_star]),
        "sar_star": float(g_sar[g_star]), "sar_star_cw": float(g_cw["averaged_sar_w_per_kg"]),
        "lattice_match": match, "n_new": int(sum(mk is None for mk in match)), "t_s": float(t_grid),
    }
    return result


def _sense_p_src(env, specs, feeds, e, w):
    """Field-partition P_src under weights ``w`` on the two feeds, WF-10's route
    (`test_wf10_f_human_tuned_b1_sar.py:576–591`)."""
    import dataclasses

    from tests.validation.test_birdcage_power_identity import _source_power_w

    by_id = {sp.port_id: sp for sp in specs}
    return float(sum(
        _source_power_w(env["acct"],
                        dataclasses.replace(by_id[pid].sheet(driven=True), source_voltage_v=complex(wk)),
                        e, env["omega"])
        for pid, wk in zip(feeds, w)))


def _p_acc_spec(s2, w, v_s, z0):
    """The lumped line cut (ruling 2026-10-12): |V_s|²/(8 z0)·(|a|² − |b|²), b = S a,
    H33's `p_acc_wave` under two incident waves."""
    a = np.asarray(w, dtype=np.complex128)
    b = np.asarray(s2, dtype=np.complex128) @ a
    return float(abs(v_s) ** 2 / (8.0 * z0) * (np.sum(np.abs(a) ** 2) - np.sum(np.abs(b) ** 2)))


def _run_box(comm, monkeypatch, c_pf, *, degree=None):
    """Entry 52's body: arm R of I36b's `_run_drive`
    (`test_geo36_rod_f_clinical_drive.py:472–630`), then `_box_readouts`.  Every
    reading printed first; asserts last, collected, so one red hides no other.
    ``degree`` (C34): ``None`` / 1 or 2, routed as I36b's `_run_drive` does."""
    if degree is None and os.environ.get(C6_DEGREE_ENV) is not None:
        raise ValueError(f"{C6_DEGREE_ENV} is read only by C6a's gated test; unset it "
                         f"(a C21 / P61 window never routes a degree)")
    deg = degree or 1
    if deg not in (1, 2):
        raise ValueError(f"degree {degree!r}: only None / 1 or 2")
    import dataclasses
    import gc
    import json

    from fem_em_solver.core import TimeHarmonicSolver
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo36_rod_f_clinical_drive as i36
    from tests.validation import test_geo38_f_clinical_column as col
    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_geo38_f_clinical_matrix import BLIND_FACTOR
    from tests.validation.test_port_birdcage_ring_column import (
        AZIMUTH_MATCH_DEG,
        _reciprocity_ratio,
    )
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf10_f_human_tuned_b1_sar import _sense_weights
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    t_window = time.perf_counter()
    grid_on = _grid_knob()  # C6b: read here, before gmsh, so C21 / P61 inherit it
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    _say(comm, f"\nstart: -n {comm.size}, C6 body; node C {c_pf:.4f} pF, degree {deg}, r = {c4.ESR_OHM} "
               f"Ohm on all 64 gaps, kappa off; arm R only (tag 403 -> gel); PREDICTED window "
               f"{PREDICTED_BOX_WINDOW_S[0]:.0f}-{PREDICTED_BOX_WINDOW_S[1]:.0f} s (printed)")
    i36._route_drive(monkeypatch)
    if deg == 2:  # I36b's `_run_drive` (`test_geo36_rod_f_clinical_drive.py:504–510`)
        monkeypatch.setenv(c4.C4_DEGREE_ENV, "2")  # C4a's fixture reads its own name
    else:
        monkeypatch.delenv(c4.C4_DEGREE_ENV, raising=False)
    fixture = c4._build_c4_fixture(comm)  # (A0) / census / (U1) / (D-rod) / (F) asserted inside
    monkeypatch.delenv(c4.C4_DEGREE_ENV, raising=False)
    assert fixture["degree"] == deg, (fixture["degree"], deg)
    if deg == 2:  # `_apply_degree_c4` asserted 2E + 2F inside; the rod mesh's record here
        n2 = fixture["counts"]["n_dofs"]
        _say(comm, f"(Deg) degree 2 global unknowns {n2} vs record "
                   f"{i35.F_CLINICAL_ROD_UNKNOWNS_RECORD[1]} (EXACT, ASSERTED)")
        assert n2 == i35.F_CLINICAL_ROD_UNKNOWNS_RECORD[1], n2
    feeds, ids = fixture["feeds"], fixture["ids"]
    ctx = fixture["built"]["ctx"]
    env = i36._drive_env(fixture, comm)
    frequency_hz = env["omega"] / (2.0 * np.pi)
    others = sorted((pid for pid in ids if pid not in feeds), key=lambda p: int(p[1:]))
    order = list(feeds) + others

    # ---- the senses are valid only on an azimuth-increasing pair -------------
    az_feeds = [float(fixture["geo"][pid][1]) for pid in feeds]
    d_az = (az_feeds[1] - az_feeds[0]) % 360.0
    _say(comm, f"feeds {feeds}: az {az_feeds[0]:.9f} / {az_feeds[1]:.9f} deg, az1 - az0 = {d_az:.9f} "
               f"deg (+90 within {AZIMUTH_MATCH_DEG:.0e}, ASSERTED)")
    assert abs(d_az - 90.0) < AZIMUTH_MATCH_DEG, d_az

    # ---- (Z) the ESR specs at the node ---------------------------------------
    c_f = c_pf * 1e-12
    base_specs = list(ctx["specs"])
    base_ids = [sp.port_id for sp in base_specs]
    specs, z_rc, z_th = c4._esr_specs(base_specs, feeds, c_f, frequency_hz, z0)
    z_fails, im_worst, th_worst = c4._z_failures(specs, base_ids, feeds, c_f, frequency_hz, z0)
    v_s = h17._source_voltage(z_rc, z0)
    _say(comm, f"(Z) node {c_pf:.4f} pF: Z_RC {z_rc:.9e} Ohm, Z_th {z_th:.9e} Ohm, V_s {v_s:.9e} V; "
               f"Im rel {im_worst:.1e}, Z_th rel {th_worst:.1e}; failures {z_fails or 'none'} (ASSERTED)")

    # ---- arm R: the gel twin, (A0m), the lattice on the real mesh -------------
    msh, cell_tags = env["msh"], env["cell_tags"]
    problem_m = ctx["solver"].problem
    problem_r = dataclasses.replace(
        problem_m, material_map=i36._gel_twin_map(dict(problem_m.material_map)))
    census_m = col._owned_tag_counts(msh, problem_m.cell_tags, comm)
    census_r = col._owned_tag_counts(msh, problem_r.cell_tags, comm)
    a0m_fails = i36._a0m_failures(census_m, census_r)
    rows = _lattice_report()
    tags = _point_tags(msh, ctx["cell_tags"], [r[1] for r in rows], comm)
    lattice_fails = [r[0] for r, t in zip(rows, tags) if r[2] > 0.0 or t != {i35.PHANTOM_CELL_TAG}]
    _say(comm, f"(A0m) census M {census_m} vs R {census_r}: failures {a0m_fails or 'none'} (ASSERTED "
               f"exact); lattice n {len(rows)}, not in tag 3 {lattice_fails or 'none'} (ASSERTED)")
    grid_rows, grid_fails = None, []
    if grid_on:
        grid_all = _grid_report()
        grid_rows = [r for r in grid_all if r[4]]
        gel_tags = {i35.PHANTOM_CELL_TAG, BIRDCAGE_IMPLANT_CELL_TAG}
        gtags = _point_tags(msh, ctx["cell_tags"], [r[1] for r in grid_rows], comm)
        grid_fails = [r[0] for r, t in zip(grid_rows, gtags) if not t or not t <= gel_tags]
        excluded = [f"{r[0]} (wall {r[2]:+.4e} m, rod-axis {r[3]:.4e} m)" for r in grid_all if not r[4]]
        _say(comm, f"(C6b) grid n {len(grid_rows)} of {len(grid_all)} candidates kept; excluded "
                   f"{excluded or 'none'}; kept points not in tags {sorted(gel_tags)} "
                   f"{grid_fails or 'none'} (ASSERTED)")
    del ctx["solver"]  # a clinical factorisation is ≈ 40 GiB; the arm-M operator is never needed
    gc.collect()
    solver_r = TimeHarmonicSolver(problem_r, degree=deg)

    # ---- one factorisation, one held column ----------------------------------
    cur, kinds, times, parts, es = [], [], [], [], []
    sigma_f = None
    for n, pid in enumerate(feeds):
        d = i36._drive(env, solver_r, specs, pid, held=(n > 0))
        kinds.append(d["kind"])
        times.append(d["t_solve"] + d["t_post"])
        cur.append(np.array([d["currents"][p] for p in order], dtype=np.complex128))
        parts.append({k: v for k, v in d.items() if k.startswith("p_")})
        es.append(d["e"])
        sigma_f = d["sigma"]
        _say(comm, f"column {pid}: {d['kind']}, solve {d['t_solve']:.2f} s, post {d['t_post']:.2f} s")
        del d
    s2 = h17._current_route_s(cur, [0, 1], z_th, v_s)
    recip = _reciprocity_ratio(s2)
    del solver_r
    gc.collect()

    fails = []
    _say(comm, f"kinds {kinds}; S2 " + "; ".join(", ".join(f"{complex(v):.10e}" for v in row) for row in s2)
               + f"; (R1) ||S-S^T||_F/||S||_F {recip:.3e} (band {RECIPROCITY_BAND:.0e}, ASSERTED)")
    if kinds != ["solve", "held"]:
        fails.append(f"kinds {kinds}")
    if not recip <= RECIPROCITY_BAND:
        fails.append(f"(R1) {recip:.3e}")
    for j, p in enumerate(parts):
        rel, blind = i36._partition(p)
        _say(comm, f"(P) column {feeds[j]}: P_src {p['p_src']:.9e} sumP_sheet {p['p_sheets']:.9e} "
                   f"P_phantom {p['p_phantom']:.9e} P_rod {p['p_rod']:.9e} P_rest {p['p_rest']:.9e} "
                   f"P_surf {p['p_surf']:.9e} W; rel {rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, "
                   f"ASSERTED); NEG blind {blind:.4e} (>= {BLIND_FACTOR:g}x band, ASSERTED)")
        if not rel <= DISCRETE_IDENTITY_RTOL:
            fails.append(f"(P) {feeds[j]} {rel:.3e}")
        if not blind >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL:
            fails.append(f"NEG {feeds[j]} {blind:.3e}")
        if not p["p_surf"] > 0.0:
            fails.append(f"P_surf {feeds[j]} {p['p_surf']!r}")
    rel_prc, _b = i36._partition(parts[0])
    p_rc, p_feeds, dz = i36._p_rc_split(parts[0], feeds)
    _say(comm, f"(P_RC) column {feeds[0]} at the node: rel {rel_prc:.3e} (band "
               f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED); 62 RC {p_rc:.9e} W, 2 feeds {p_feeds:.9e} W; "
               f"(Dz) without the RC sheets {dz:.4e} (printed)")

    # ---- the two senses, one linear combination each --------------------------
    weights = _sense_weights()
    fields = {}
    for sense in ("ccw", "cw"):
        f = es[0].copy()
        f.x.array[:] = weights[sense][0] * es[0].x.array + weights[sense][1] * es[1].x.array
        f.x.scatter_forward()
        fields[sense] = f
    y_in = {sense: h17._active_y_in(s2, weights[sense], z0) for sense in fields}
    y_series = i36._y_in_1p1j(s2, z0)
    for sense, ys in y_in.items():
        _say(comm, f"Y_in {sense} {tuple(weights[sense])}: {feeds[0]} {ys[0].real:+.9e} "
                   f"{ys[0].imag:+.9e}j S, {feeds[1]} {ys[1].real:+.9e} {ys[1].imag:+.9e}j S (Im "
                   f"printed: tuned on C4's cylinder, never a red)")
    _say(comm, f"C4's root series is {feeds[0]} under (1,+1j), here the cw sense: "
               f"{y_series[0].real:+.9e} {y_series[0].imag:+.9e}j S")

    # ---- readouts ------------------------------------------------------------
    t_post = time.perf_counter()
    rho_field = _gel_density(msh, cell_tags)
    sample_points = _box_sample_points(msh, cell_tags, comm)
    out = _box_readouts(msh, cell_tags, fields["ccw"], fields["cw"], sigma_f, env["omega"], rho_field,
                        rows, sample_points, comm, grid_rows=grid_rows)
    p_tot, p_gel, p_rod = i36._vol_powers(env, sigma_f, fields["ccw"])
    p_src = _sense_p_src(env, specs, feeds, fields["ccw"], weights["ccw"])
    p_acc = _p_acc_spec(s2, weights["ccw"], v_s, z0)
    t_post = comm.allreduce(time.perf_counter() - t_post, op=MPI.MAX)
    ccw, cw = out["senses"]["ccw"], out["senses"]["cw"]
    p_box = out["tagged"]["dissipated_power_w"]
    p_cover = out["whole"]["dissipated_power_w"]
    x, y, z = out["c_star"]
    _say(comm, f"mean |B1+|_ccw {ccw['mean_b_plus']:.9e} T over {ccw['n_valid']} / "
               f"{len(sample_points)} tag-3 centroids, mean |B1-| {ccw['mean_b_minus']:.9e} T, (Sense) "
               f"{ccw['sense_ratio']:.6f} (printed); CV {ccw['cv']:.6f} (printed, no homogeneity "
               f"claim); (1,+1j) mean |B1+| {cw['mean_b_plus']:.9e} T; lattice z spread "
               f"{ccw['lattice_z_spread']:.6f}, all valid {ccw['lattice_valid']}")
    for k, (label, s) in enumerate(zip(out["labels"], out["sar10g"])):
        _say(comm, f"lattice {label}: |B1+|_ccw {ccw['lattice_b_plus'][k]:.9e} T, SAR10g {s:.9e} "
                   f"W/kg, mass {out['ball_mass_kg'][k]:.6e} kg")
    _say(comm, f"c* {out['star_label']}: SAR10g {out['sar_star']:.9e} W/kg (ccw), (1,+1j) ball "
               f"{out['sar_star_cw']:.9e} W/kg; r {np.hypot(x, y):.6f} m, dz {z:+.6f} m, az "
               f"{np.degrees(np.arctan2(y, x)):.3f} deg")
    grid = out.get("grid")
    if grid is not None:
        gx, gy, gz = grid["c_star"]
        _say(comm, f"(C6b) grid n {grid['n']} ({grid['n_new']} new balls, "
                   f"{grid['n'] - grid['n_new']} lattice values reused); g* {grid['labels'][grid['star']]} "
                   f"({gx:.6f}, {gy:.6f}, {gz:.6f}) m: SAR10g {grid['sar_star']:.9e} W/kg (ccw), (1,+1j) "
                   f"ball {grid['sar_star_cw']:.9e} W/kg; beside c* {out['star_label']} "
                   f"({x:.6f}, {y:.6f}, {z:.6f}) m: SAR10g {out['sar_star']:.9e} W/kg; grid max >= "
                   f"lattice max (ASSERTED); g*/c* {grid['sar_star'] / out['sar_star']:.6f}")
        if grid["n_new"] != grid["n"] - len(rows):
            fails.append(f"(C6b) lattice not a grid subset: {grid['n'] - grid['n_new']} of {len(rows)} matched")
        if not grid["sar_star"] >= out["sar_star"]:
            fails.append(f"(C6b) grid max {grid['sar_star']:.9e} < lattice max {out['sar_star']:.9e}")
    if grid_fails:
        fails.append(f"(C6b) grid points outside the gel {grid_fails}")
    _say(comm, f"(i) ball radius {_coverage_radius():.6f} m on tags (3, 403): power {p_cover:.12e} vs "
               f"{p_box:.12e} W, miss {out['power_miss']:.3e}; mass {out['whole']['mass_kg']:.9e} kg, "
               f"miss {out['mass_miss']:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED)")
    _say(comm, f"(Sigma_P) P_gel + P_tag403 {p_gel + p_rod:.12e} W (P_total {p_tot:.12e}) vs mean_sar "
               f"{p_box:.12e} W, rel {abs(p_gel + p_rod - p_box) / abs(p_box):.3e} (printed)")
    _say(comm, f"normalisers: P_box {p_box:.9e} W, P_src (field partition, ccw) {p_src:.9e} W, P_acc "
               f"(lumped line cut, ruling 2026-10-12) {p_acc:.9e} W; P_box/P_src {p_box / p_src:.6e}, "
               f"P_box/P_acc {p_box / p_acc:.6e}; SAR10g(c*) per box watt "
               f"{out['sar_star'] / p_box:.6e}, per P_acc {out['sar_star'] / p_acc:.6e}, per P_src "
               f"{out['sar_star'] / p_src:.6e} W/kg/W; mean |B1+|_ccw / sqrt(W): box "
               f"{ccw['mean_b_plus'] / np.sqrt(p_box):.6e}, P_acc "
               f"{ccw['mean_b_plus'] / np.sqrt(p_acc):.6e} T")
    for name, miss in (("power", out["power_miss"]), ("mass", out["mass_miss"])):
        if not miss <= EXACT_IDENTITY_RTOL:
            fails.append(f"(i) {name} {miss:.3e}")
    if comm.rank == 0:
        os.makedirs(C6_JSON_DIR, exist_ok=True)
        suffix = "-d2" if deg == 2 else ""  # `_i36_json_path`'s pattern; no new key
        path = os.path.join(C6_JSON_DIR, f"c6-node-{c_pf:.4f}{suffix}.json")
        rec = {
            "c_pf": c_pf, "feeds": list(feeds), "az_feeds": az_feeds,
            "s2": [[[complex(v).real, complex(v).imag] for v in row] for row in s2],
            "y_in": {s: [[y_.real, y_.imag] for y_ in ys] for s, ys in y_in.items()},
            "mean_b_plus": {s: out["senses"][s]["mean_b_plus"] for s in out["senses"]},
            "mean_b_minus": {s: out["senses"][s]["mean_b_minus"] for s in out["senses"]},
            "lattice": [{"label": lb, "c": list(r[1]), "b_plus_ccw": float(bp), "sar10g": float(s)}
                        for lb, r, bp, s in zip(out["labels"], rows, ccw["lattice_b_plus"],
                                                out["sar10g"])],
            "c_star": list(out["c_star"]), "star_label": out["star_label"],
            "sar_star": out["sar_star"], "sar_star_cw": out["sar_star_cw"],
            "p_box": p_box, "p_src": p_src, "p_acc": p_acc, "parts": [
                {k: v for k, v in p.items() if k != "p_sheet_by_id"} for p in parts],
        }
        if grid is not None:  # the key exists only with the knob on (control (G-off))
            rec["grid"] = {
                "n": grid["n"], "labels": list(grid["labels"]),
                "points": [list(map(float, p)) for p in grid["points"]],
                "sar10g": [float(v) for v in grid["sar10g"]], "star": grid["star"],
                "star_label": grid["labels"][grid["star"]], "c_star": list(grid["c_star"]),
                "sar_star": grid["sar_star"], "sar_star_cw": grid["sar_star_cw"],
                "lattice_match": grid["lattice_match"],
            }
        with open(path, "w") as fh:
            json.dump(rec, fh, indent=1)
        print(f"[GEO-38 C6a] wrote {path}", flush=True)
    comm.Barrier()
    _m, per, tot = _rss_gib(comm)
    _say(comm, f"price: build {fixture['t_build']:.1f} s; columns {', '.join(f'{t:.1f}' for t in times)} "
               f"s; readouts {t_post:.1f} s; RSS per rank {', '.join(f'{r:.2f}' for r in (per or []))} "
               f"GiB, summed {tot:.3f} GiB; window {time.perf_counter() - t_window:.1f} s; tag "
               f"{BIRDCAGE_IMPLANT_CELL_TAG} counted as gel"
               + (f"; grid {grid['t_s']:.1f} s" if grid is not None else ""))
    assert not z_fails, z_fails
    assert not a0m_fails, a0m_fails
    assert not lattice_fails, lattice_fails
    assert rel_prc <= DISCRETE_IDENTITY_RTOL, f"(P_RC) {rel_prc:.3e} — never re-band"
    assert not fails, fails


@complex_only
def test_geo38_c6a_box_loaded(monkeypatch):
    mode, _c_pf = _knobs()  # node knob missing under 1 / set under build => ValueError, before gmsh
    _degree = _c6_degree(mode)  # C34: set under build => ValueError, before gmsh
    if mode is None:
        pytest.skip(f"{C6_ENV} unset: GEO-38 C6a is selected by environment only")
    if mode == "1":
        _run_box(MPI.COMM_WORLD, monkeypatch, _c_pf, degree=_degree)
        return
    _run_build(MPI.COMM_WORLD, monkeypatch)


# ---------------------------------------------------------------------------
# Coding tests (no gmsh); skip whenever C6_ENV is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(os.environ.get(C6_ENV, "").strip() != "",
                             reason=f"{C6_ENV} set: coding tests skip in build windows")


@_coding
def test_c6a_k1_knob_parsing(monkeypatch):
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    for k in (C6_ENV, S8_DRIVE_C_PF_ENV) + FOREIGN_ENVS:
        monkeypatch.delenv(k, raising=False)
    assert _knobs() == (None, None)
    monkeypatch.setenv(C6_ENV, "build")
    assert _knobs() == ("build", None)
    monkeypatch.setenv(S8_DRIVE_C_PF_ENV, "193.5178")
    with pytest.raises(ValueError, match="refused"):
        _knobs()
    monkeypatch.setenv(C6_ENV, "1")
    assert _knobs() == ("1", 193.5178)
    for bad in ("2", "x"):
        monkeypatch.setenv(C6_ENV, bad)
        with pytest.raises(ValueError):
            _c6_mode()
    monkeypatch.setenv(C6_ENV, "build")
    monkeypatch.delenv(S8_DRIVE_C_PF_ENV)
    for k in FOREIGN_ENVS:
        with monkeypatch.context() as m:
            m.setenv(k, "1")
            with pytest.raises(ValueError):
                _knobs()
    print(f"\n[GEO-38 C6a K1] mode unset/build/1 accepted, 2/x refused; node knob refused under "
          f"build, read under 1; {len(FOREIGN_ENVS)} foreign variables refused (ASSERTED)")


@_coding
def test_c6a_k2_one_without_node_stops_before_gmsh(monkeypatch):
    from fem_em_solver.io.mesh import MeshGenerator
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    def _sentinel(*a, **k):
        raise AssertionError("gmsh reached")

    with monkeypatch.context() as m:
        for k in FOREIGN_ENVS + (S8_DRIVE_C_PF_ENV,):
            m.delenv(k, raising=False)
        m.setenv(C6_ENV, "1")
        m.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(_sentinel))
        with pytest.raises(ValueError, match=S8_DRIVE_C_PF_ENV):
            _knobs()
    print("\n[GEO-38 C6a K2] C6=1 without the node knob raised ValueError before gmsh (ASSERTED)")


@_coding
def test_c6a_k3_lattice_on_box_constants():
    """The ruled lattice (§9 item 693) on the box constants, no mesh: count = 18
    and every 10 g ball ≥ 1 mm clear of every wall and of the rod capsule.  Every
    row printed before the asserts; the attempt-1 lattice (``LATTICE_DZ_M``, nine
    points) printed through the same arithmetic as the control, never asserted."""
    from tests.mesh.test_birdcage_implant_rod import ROD_LENGTH, ROD_RADIUS
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation.test_wf7_f_human_sar10g_hotspot import (
        CONTAINMENT_MARGIN_M,
        LATTICE_DZ_M,
    )

    centre, half = _box()
    a = _a10()
    rows = _lattice_report()
    print(f"\n[GEO-38 C6a K3] box centre {tuple(centre)} m, half extents {tuple(half)} m "
          f"(Lx, Ly, Lz = {tuple(2 * half)}); rod centre {i35._rod_centre()}, ROD_LENGTH "
          f"{ROD_LENGTH}, ROD_RADIUS {ROD_RADIUS}; a10 {a:.6e} m; C6A_LATTICE_DZ_M "
          f"{C6A_LATTICE_DZ_M} ({len(C6A_LATTICE_DZ_M)} heights); rod clearance bar "
          f"{ROD_RADIUS + a + CONTAINMENT_MARGIN_M:.4e} m")
    for label, c, wall, rd, g in _lattice_report(LATTICE_DZ_M):
        print(f"[GEO-38 C6a K3] control (attempt-1 lattice, PRINTED) {label}: wall {wall:+.4e} m, "
              f"rod-axis {rd:.4e} m, in gel {g}")
    for label, c, wall, rd, g in rows:
        print(f"[GEO-38 C6a K3] {label} {tuple(round(v, 6) for v in c)}: wall {wall:+.4e} m "
              f"(<= 0 contained), min wall distance {a - wall + 1e-3:.4e} m, rod-axis {rd:.4e} m, "
              f"in gel {g}")
    out_gel = [r[0] for r in rows if not r[4]]
    out_wall = [r[0] for r in rows if r[2] > 0.0]
    print(f"[GEO-38 C6a K3] n {len(rows)} vs brief {BRIEF_LATTICE_N}; wall failures "
          f"{out_wall or 'none'}; not in gel {out_gel or 'none'}")
    assert len(rows) == BRIEF_LATTICE_N, f"lattice has {len(rows)} points, brief says {BRIEF_LATTICE_N}"
    assert not out_wall, out_wall
    assert not out_gel, out_gel


@_coding
def test_c6a_k4_coverage_radius_reaches_every_corner():
    """(K4): half-diagonal + 1 mm, every corner of the box inside (ruling 2026-10-30
    18:00, defect 2); pure python."""
    import itertools

    centre, half = _box()
    radius = _coverage_radius()
    expected = float(np.sqrt(0.21 ** 2 + 0.045 ** 2 + 0.325 ** 2)) + 1.0e-3
    corners = [centre + half * np.array(s) for s in itertools.product((-1.0, 1.0), repeat=3)]
    far = max(float(np.linalg.norm(c - centre)) for c in corners)
    print(f"\n[GEO-38 C6a K4] radius {radius:.12f} m, expected {expected:.12f} m, farthest corner "
          f"{far:.12f} m (ASSERTED)")
    assert abs(radius - expected) <= 1.0e-12
    assert round(radius, 6) == 0.390551
    assert far < radius


@_coding
def test_c6a_k5_sense_weights():
    """(K5): ccw = (1, −1j), cw = (1, +1j) (ruling 2026-10-30 18:00, defect 1)."""
    from tests.validation.test_wf10_f_human_tuned_b1_sar import _sense_weights

    w = _sense_weights()
    d = [max(abs(np.asarray(w[s]) - np.array(t))) for s, t in (("ccw", (1.0, -1.0j)), ("cw", (1.0, 1.0j)))]
    print(f"\n[GEO-38 C6a K5] ccw {w['ccw']}, cw {w['cw']}; misses {d} (<= 1e-15, ASSERTED)")
    assert max(d) <= 1.0e-15, d


def _k6_fixture(comm):
    """K6's box: tags 3 and 403 (centroid max-norm < 0.03 m), E = c·½ x̂ × r."""
    import dolfinx
    from dolfinx import fem

    half = np.array([0.1, 0.05, 0.15])
    msh = dolfinx.mesh.create_box(comm, [(-half).tolist(), half.tolist()], [10, 6, 14],
                                  dolfinx.mesh.CellType.tetrahedron)
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    cells = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
    mids = dolfinx.mesh.compute_midpoints(msh, tdim, cells)
    values = np.where(np.max(np.abs(mids), axis=1) < 0.03, 403, 3).astype(np.int32)
    cell_tags = dolfinx.mesh.meshtags(msh, tdim, cells, values)
    c = 1.0 + 0.5j
    v = fem.functionspace(msh, ("N1curl", 1))

    def _field(factor):
        e = fem.Function(v)
        e.interpolate(lambda x: factor * 0.5 * c * np.vstack([0.0 * x[0], -x[2], x[1]]))
        e.x.scatter_forward()
        return e

    sigma_f = fem.Function(fem.functionspace(msh, ("DG", 0)))
    sigma_f.x.array[:] = 0.5
    return msh, cell_tags, _field(1.0), _field(2.0), sigma_f, c, half


@complex_only
@_coding
def test_c6a_k6_readouts_on_a_constant_curl_field():
    """(K6): the readouts' closed form on a small mesh.  (K6-B) |B₁^±| = |c|/(2ω) and
    ``e_cw`` twice that; (K6-i) coverage on tags (3, 403); (K6-c) the negative control:
    the identity on tag 3 alone with 403's density zeroed in only one side misses."""
    from fem_em_solver.post.sar import build_density_field
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import (
        EXACT_IDENTITY_RTOL,
        PHANTOM_RHO_KG_PER_M3,
    )
    from tests.validation.test_wf7_f_human_sar10g_hotspot import COVERAGE_BALL_PAD_M

    comm = MPI.COMM_WORLD
    msh, cell_tags, e_ccw, e_cw, sigma_f, c, half = _k6_fixture(comm)
    omega = 2.0 * np.pi * 64.0e6
    centre = np.zeros(3)
    radius = float(np.linalg.norm(half) + COVERAGE_BALL_PAD_M)
    rho_field = _gel_density(msh, cell_tags)
    points = _box_sample_points(msh, cell_tags, comm)
    assert np.array_equal(points, comm.bcast(points, root=0))
    rows = [(f"k6-{k}", (0.0, 0.0, dz), None, None, True)
            for k, dz in enumerate((0.0, 0.05, 0.1, -0.05))]
    out = _box_readouts(msh, cell_tags, e_ccw, e_cw, sigma_f, omega, rho_field, rows, points, comm,
                        coverage=(centre, radius))
    ccw, cw = out["senses"]["ccw"], out["senses"]["cw"]
    b0 = abs(c) / (2.0 * omega)
    miss_b = {
        "ccw+": float(np.max(np.abs(ccw["b_plus"][ccw["valid"]] / b0 - 1.0))),
        "ccw-": float(np.max(np.abs(ccw["b_minus"][ccw["valid"]] / b0 - 1.0))),
        "cw+": float(np.max(np.abs(cw["b_plus"][cw["valid"]] / (2.0 * b0) - 1.0))),
    }
    r_star = abs(out["sar_star_cw"] / out["sar_star"] - 4.0) / 4.0

    # (K6-c): one side keeps tag 403, the other does not (defect 3, both ways round).
    rho_403_zero = build_density_field(msh, PHANTOM_RHO_KG_PER_M3, cell_tags=cell_tags,
                                       density_map={403: 0.0})
    _w, tag3, p_miss_c, m_miss_zero = _coverage_identity(e_ccw, sigma_f, rho_403_zero, cell_tags, 3,
                                                         centre, radius, comm)
    _w, _t, _p, m_miss_c = _coverage_identity(e_ccw, sigma_f, rho_field, cell_tags, 3, centre, radius,
                                              comm)
    v_all = out["tagged"]["volume_m3"]
    v_403 = v_all - tag3["volume_m3"]
    if comm.rank == 0:
        print(f"\n[GEO-38 C6a K6] {len(points)} tag-3 centroids, {ccw['n_valid']} valid; b0 = "
              f"|c|/(2 omega) {b0:.9e} T; (K6-B) misses {miss_b} (band {DISCRETE_IDENTITY_RTOL:g}, "
              f"ASSERTED); (Sense) {ccw['sense_ratio']:.12f}, CV {ccw['cv']:.3e}, lattice z spread "
              f"{ccw['lattice_z_spread']:.3e} (printed)")
        print(f"[GEO-38 C6a K6] c* {out['star_label']} {out['c_star']}, SAR10g {out['sar10g'].tolist()} "
              f"W/kg; cw/ccw ball at c* rel to 4 {r_star:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED)")
        print(f"[GEO-38 C6a K6] (K6-i) radius {radius:.9f} m: power miss {out['power_miss']:.3e}, mass "
              f"miss {out['mass_miss']:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED)")
        print(f"[GEO-38 C6a K6] (K6-c) tag 3 alone: ball rho(403) = 0, power miss {p_miss_c:.6e} (mass "
              f"then matches, {m_miss_zero:.3e}); ball rho(403) = gel, mass miss {m_miss_c:.6e}; V403/V "
              f"{v_403 / v_all:.6e}, V403/V3 {v_403 / tag3['volume_m3']:.6e} (> band, ASSERTED)",
              flush=True)
    assert len(points) > 0 and ccw["n_valid"] == len(points) and cw["n_valid"] == len(points)
    assert float(np.min(np.max(np.abs(points), axis=1))) >= 0.03  # tag 403's cells are not sampled
    assert max(miss_b.values()) <= DISCRETE_IDENTITY_RTOL, miss_b
    assert out["star_label"] == "k6-2", out["star_label"]  # |E| grows with |z| on the axis
    assert r_star <= EXACT_IDENTITY_RTOL, r_star
    assert out["power_miss"] <= EXACT_IDENTITY_RTOL, out["power_miss"]
    assert out["mass_miss"] <= EXACT_IDENTITY_RTOL, out["mass_miss"]
    assert p_miss_c > EXACT_IDENTITY_RTOL, p_miss_c
    assert m_miss_c > EXACT_IDENTITY_RTOL, m_miss_c


class _FixtureReached(Exception):
    pass


@_coding
def test_c6a_k7_one_with_node_reaches_run_box(monkeypatch):
    """(K7): ``1`` with the node knob reaches `_run_box`, stopping at C4a's fixture
    (a sentinel, so no gmsh)."""
    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation.test_geo36_rod_f_clinical_drive import S8_DRIVE_C_PF_ENV

    def _sentinel(_comm):
        raise _FixtureReached

    with monkeypatch.context() as m:
        for k in FOREIGN_ENVS:
            m.delenv(k, raising=False)
        m.setenv(C6_ENV, "1")
        m.setenv(S8_DRIVE_C_PF_ENV, "193.5996")
        m.setattr(c4, "_build_c4_fixture", _sentinel)
        with pytest.raises(_FixtureReached):
            test_geo38_c6a_box_loaded(m)
    print("\n[GEO-38 C6a K7] C6=1 with the node knob reached _run_box's fixture call (ASSERTED)")


#: K8's node (the 193.5996 pF of ``xl`` entry 52's command).
K8_NODE_PF = 193.5996


def _k8_stubs(m, fx, json_dir, *, p_src_p2=1.0):
    """§9 item 701: everything `_run_box` reads around its own lines, stubbed on K6's
    mesh.  Column 1 is E, column 2 is 0.5j E; the powers close TH-14's identity
    exactly unless ``p_src_p2`` breaks column 2's.  Returns ``(lines, seen, z_sym)``."""
    import dataclasses
    import sys

    import ufl
    from dolfinx import fem

    import fem_em_solver.core as core
    from fem_em_solver.core import HomogeneousMaterial
    from tests.validation import test_birdcage_power_identity as power_identity
    from tests.validation import test_geo36_rod_f_clinical_drive as i36
    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation import test_wf10_f_human_copper_esr_root as h17

    @dataclasses.dataclass
    class _Sheet:
        port_id: str
        source_voltage_v: complex = 1.0

        def sheet(self, driven=False):
            return self

    @dataclasses.dataclass
    class _Problem:
        material_map: dict
        cell_tags: object

    comm = MPI.COMM_WORLD
    msh, cell_tags, e, _e2, sigma_f, _c, half = fx
    omega = 2.0 * np.pi * 64.0e6
    e_2 = fem.Function(e.function_space)
    e_2.x.array[:] = 0.5j * e.x.array
    e_2.x.scatter_forward()
    columns = (e, e_2)

    feeds, order = ["P1", "P2"], ["P1", "P2", "P3", "P4"]
    specs = [_Sheet(pid) for pid in order]
    gel = HomogeneousMaterial(sigma=0.5, epsilon_r=78.0, mu_r=1.0)
    rod = HomogeneousMaterial(sigma=1.0e6, epsilon_r=78.0, mu_r=1.0)
    problem = _Problem(material_map={3: gel, 403: rod}, cell_tags=cell_tags)
    ctx = {"specs": specs, "solver": type("_Solver", (), {"problem": problem})(), "msh": msh,
           "cell_tags": cell_tags}
    fixture = {"degree": 1, "feeds": feeds, "ids": order,
               "geo": {pid: (0.35, 90.0 * k) for k, pid in enumerate(order)},
               "built": {"ctx": ctx}, "t_build": 0.0}
    env = {"comm": comm, "msh": msh, "cell_tags": cell_tags, "omega": omega,
           "dx": ufl.Measure("dx", domain=msh, subdomain_data=cell_tags), "acct": {"mesh": msh}}

    # K5's reciprocal synthetic circuit (`_k8_arm`, `test_geo36_rod_f_clinical_drive.py:993–1001`).
    z0 = 50.0
    rng = np.random.default_rng(656)
    a = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    z_sym = a + a.T + np.eye(4) * (20.0 + 5.0j)
    z_rc = h17._z_rc(K8_NODE_PF * 1e-12, 64.0e6)
    z_th = h17._thevenin_feed_impedance(z_rc, z0)
    loads = h17._parallel_loads(4, [0, 1], z_rc, z_th)

    def _drive(_env, _solver, _specs, driven_id, *, held):
        j = feeds.index(driven_id)
        v = np.zeros(4, dtype=np.complex128)
        v[j] = 1.0
        cur = h17._parallel_currents(z_sym, loads, v)
        return {"currents": dict(zip(order, cur)), "kind": "held" if held else "solve",
                "t_solve": 0.0, "t_post": 0.0, "e": columns[j].copy(), "sigma": sigma_f,
                "p_src": p_src_p2 if j == 1 else 1.0, "p_sheets": 0.5, "p_phantom": 0.2,
                "p_rest": 0.2, "p_surf": 0.1, "p_rod": 0.0,
                "p_sheet_by_id": {"P1": 0.1, "P2": 0.1, "P3": 0.15, "P4": 0.15}}

    seen = []

    def _source_power_w(_acct, sheet, _e, _omega):
        seen.append((sheet.port_id, complex(sheet.source_voltage_v)))
        return 1.0

    lines = []
    module = sys.modules[_run_box.__module__]
    real_say = module._say

    def _say_recorded(comm_, msg):
        lines.append(msg)
        real_say(comm_, msg)

    k6_rows = [(f"k8-{k}", (0.0, 0.0, dz), -1.0, None, True)
               for k, dz in enumerate((0.05, 0.1, -0.05, -0.1))]  # K6's rows clear of tag 403
    radius = float(np.linalg.norm(half) + 1.0e-3)

    m.setattr(i36, "_route_drive", lambda _monkeypatch: None)
    m.setattr(c4, "_build_c4_fixture", lambda _comm: fixture)
    m.setattr(i36, "_drive_env", lambda _fixture, _comm: env)
    m.setattr(c4, "_esr_specs", lambda *_a: (specs, z_rc, z_th))
    m.setattr(c4, "_z_failures", lambda *_a: ([], 0.0, 0.0))
    m.setattr(core, "TimeHarmonicSolver", lambda *_a, **_k: object())
    m.setattr(i36, "_drive", _drive)
    m.setattr(power_identity, "_source_power_w", _source_power_w)
    m.setattr(module, "_lattice_report", lambda: k6_rows)
    m.setattr(module, "_box", lambda: (np.zeros(3), half))
    m.setattr(module, "_coverage_radius", lambda: radius)
    m.setattr(module, "C6_JSON_DIR", json_dir)
    m.setattr(module, "_say", _say_recorded)
    return lines, seen, z_sym


def _k8_json(comm, directory):
    """The node's JSON, parsed on rank 0 and identical on every rank."""
    import json

    rec = None
    if comm.rank == 0:
        with open(os.path.join(directory, f"c6-node-{K8_NODE_PF:.4f}.json")) as fh:
            rec = json.load(fh)
    return comm.bcast(rec, root=0)


@complex_only
@_coding
def test_c6a_k8_run_box_end_to_end_on_stubs(monkeypatch, tmp_path):
    """(K8), §9 item 701: `_run_box`'s every line on K6's mesh, around stubbed FEM.
    Column 2 = 0.5j × column 1, so ccw = (1, −1j) is 1.5 E and cw = (1, +1j) is 0.5 E:
    (K8-B) |B₁⁺| ratio 3, (K8-S) SAR ratio 9, (K8-Ps) the P_src recorder, (K8-Pa)
    `_p_acc_spec`, (K8-J) the JSON, (K8-L) the last prints.  Control: column 2's P_src
    perturbed by 1e-3 raises naming (P) P2, and the JSON is still written."""
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_wf10_f_human_tuned_b1_sar import _sense_weights

    comm = MPI.COMM_WORLD
    fx = _k6_fixture(comm)
    z0 = 50.0
    root = comm.bcast(str(tmp_path), root=0)
    good_dir, neg_dir = os.path.join(root, "good"), os.path.join(root, "neg")

    with monkeypatch.context() as m:
        lines, seen, z_sym = _k8_stubs(m, fx, good_dir)
        _run_box(comm, m, K8_NODE_PF)
    rec = _k8_json(comm, good_dir)

    with monkeypatch.context() as m:
        _lines, _seen, _z = _k8_stubs(m, fx, neg_dir, p_src_p2=1.0 + 1.0e-3)
        with pytest.raises(AssertionError, match=r"\(P\) P2") as neg:
            _run_box(comm, m, K8_NODE_PF)
    rec_neg = _k8_json(comm, neg_dir)

    r_b = abs(rec["mean_b_plus"]["ccw"] / rec["mean_b_plus"]["cw"] / 3.0 - 1.0)
    r_s = abs(rec["sar_star"] / rec["sar_star_cw"] / 9.0 - 1.0)
    w = _sense_weights()["ccw"]
    expected_seen = [("P1", complex(w[0])), ("P2", complex(w[1]))]
    z_rc = h17._z_rc(K8_NODE_PF * 1e-12, 64.0e6)
    s_route = h17._circuit_current_route_s(z_sym, [0, 1], z_rc, z_rc, z0)
    s_json = np.array([[complex(*v) for v in row] for row in rec["s2"]])
    r_s2 = float(np.linalg.norm(s_json - s_route) / np.linalg.norm(s_route))

    v_s = 2.0 + 1.0j
    p_zero = _p_acc_spec(np.zeros((2, 2)), w, v_s, z0)
    r_pa = abs(p_zero / (abs(v_s) ** 2 * 2.0 / (8.0 * z0)) - 1.0)
    unitary = np.array([[0.6, 0.8j], [0.8j, 0.6]])
    p_unitary = abs(_p_acc_spec(unitary, w, v_s, z0))

    keys = {"c_pf", "feeds", "az_feeds", "s2", "y_in", "mean_b_plus", "mean_b_minus", "lattice",
            "c_star", "star_label", "sar_star", "sar_star_cw", "p_box", "p_src", "p_acc", "parts"}
    heads = ("(i) ball radius", "normalisers:", "price:")
    missing_heads = [h for h in heads if not any(ln.startswith(h) for ln in lines)]
    if comm.rank == 0:
        print(f"\n[GEO-38 C6a K8] (K8-B) ccw/cw miss from 3 {r_b:.3e} (band "
              f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED); (K8-S) miss from 9 {r_s:.3e} (band "
              f"{EXACT_IDENTITY_RTOL:g}, ASSERTED)")
        print(f"[GEO-38 C6a K8] (K8-Ps) seen {seen}, JSON p_src {rec['p_src']!r} (ASSERTED); "
              f"(K8-Pa) S=0 miss {r_pa:.3e}, unitary {p_unitary:.3e} (ASSERTED)")
        print(f"[GEO-38 C6a K8] (K8-J) keys {sorted(rec)}; s2 vs circuit route {r_s2:.3e} (band "
              f"{EXACT_IDENTITY_RTOL:g}, ASSERTED); (K8-L) missing {missing_heads or 'none'}")
        print(f"[GEO-38 C6a K8] control: raised {str(neg.value).splitlines()[0]}; JSON still "
              f"written, parts[1] p_src {rec_neg['parts'][1]['p_src']!r} (ASSERTED)", flush=True)
    assert r_b <= DISCRETE_IDENTITY_RTOL, r_b
    assert r_s <= EXACT_IDENTITY_RTOL, r_s
    assert seen == expected_seen, seen
    assert rec["p_src"] == 2.0, rec["p_src"]
    assert r_pa <= 1.0e-12, r_pa
    assert p_unitary <= 1.0e-12, p_unitary
    assert set(rec) == keys, sorted(rec)
    assert r_s2 <= EXACT_IDENTITY_RTOL, r_s2
    assert not missing_heads, missing_heads
    assert rec_neg["parts"][1]["p_src"] == 1.0 + 1.0e-3, rec_neg["parts"][1]
    assert "(P) P2" in str(neg.value)


# ---------------------------------------------------------------------------
# C34 (§9 item 820): the degree-2 knob, routed on stubs.
# ---------------------------------------------------------------------------
#: K8's exact JSON key set (`test_c6a_k8_run_box_end_to_end_on_stubs`); degree 2 adds none.
K8_JSON_KEYS = frozenset({"c_pf", "feeds", "az_feeds", "s2", "y_in", "mean_b_plus", "mean_b_minus",
                          "lattice", "c_star", "star_label", "sar_star", "sar_star_cw", "p_box",
                          "p_src", "p_acc", "parts"})


def _c34_recorders(m, *, degree):
    """C34: wrap the stubbed `c4._build_c4_fixture` (reporting ``degree``, with the
    rod mesh's degree-2 count) and `core.TimeHarmonicSolver`; record C4a's degree
    variable at the fixture call and every solver's ``degree``.  Installed after
    `_k8_stubs` / `_k9_stubs`, whose fixtures hard-code degree 1 (never edited)."""
    import fem_em_solver.core as core
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo38_f_clinical_esr_root as c4

    stub_fixture, stub_solver = c4._build_c4_fixture, core.TimeHarmonicSolver
    env_seen, degrees = [], []

    def _fixture(comm_):
        env_seen.append(os.environ.get(c4.C4_DEGREE_ENV))
        fx = dict(stub_fixture(comm_))
        fx["degree"] = degree
        fx["counts"] = {"n_dofs": i35.F_CLINICAL_ROD_UNKNOWNS_RECORD[1]}
        return fx

    def _solver(*a, **k):
        degrees.append(k.get("degree"))
        return stub_solver(*a, **k)

    m.setattr(c4, "_build_c4_fixture", _fixture)
    m.setattr(core, "TimeHarmonicSolver", _solver)
    return env_seen, degrees


@_coding
def test_c34_kdeg_knob_parsing(monkeypatch):
    """(K-deg): unset and ``2`` accepted under ``None`` / ``1``; ``1``, ``x`` refused;
    any value under ``build`` refused; a stray knob through `_run_box(degree=None)`
    refused before any fixture."""
    monkeypatch.delenv(C6_DEGREE_ENV, raising=False)
    assert _c6_degree(None) is None and _c6_degree("1") is None and _c6_degree("build") is None
    monkeypatch.setenv(C6_DEGREE_ENV, "")
    assert _c6_degree("1") is None
    monkeypatch.setenv(C6_DEGREE_ENV, "2")
    assert _c6_degree(None) == 2 and _c6_degree("1") == 2
    for bad in ("1", "x", "3"):
        monkeypatch.setenv(C6_DEGREE_ENV, bad)
        with pytest.raises(ValueError, match="only unset"):
            _c6_degree("1")
    for v in ("2", "x"):
        monkeypatch.setenv(C6_DEGREE_ENV, v)
        with pytest.raises(ValueError):
            _c6_degree("build")
    monkeypatch.setenv(C6_DEGREE_ENV, "2")
    with pytest.raises(ValueError, match="refused under"):
        _c6_degree("build")
    with pytest.raises(ValueError, match="read only by C6a"):
        _run_box(MPI.COMM_WORLD, monkeypatch, K8_NODE_PF)
    with pytest.raises(ValueError, match="only None / 1 or 2"):
        _run_box(MPI.COMM_WORLD, monkeypatch, K8_NODE_PF, degree=3)
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-38 C34 K-deg] {C6_DEGREE_ENV}: unset / '' / 2 accepted, 1 / x / 3 refused; "
              f"2 / x under build refused; a stray knob through _run_box(degree=None) refused "
              f"(ASSERTED)", flush=True)


# ---------------------------------------------------------------------------
# C6b coding tests (§9 item 819).
# ---------------------------------------------------------------------------
@_coding
def test_c6b_k0_grid_knob(monkeypatch):
    """The knob: unset / empty ⇒ off, ``1`` ⇒ on, anything else ⇒ ValueError."""
    monkeypatch.delenv(C6B_GRID_ENV, raising=False)
    assert _grid_knob() is False
    monkeypatch.setenv(C6B_GRID_ENV, "")
    assert _grid_knob() is False
    monkeypatch.setenv(C6B_GRID_ENV, "1")
    assert _grid_knob() is True
    for bad in ("0", "2", "on", "true", "x"):
        monkeypatch.setenv(C6B_GRID_ENV, bad)
        with pytest.raises(ValueError):
            _grid_knob()
    print(f"\n[GEO-38 C6b K0] {C6B_GRID_ENV}: unset/empty off, 1 on, 0/2/on/true/x refused (ASSERTED)")


def _g1_check(tag):
    """(G1) + (G-c) + (G-x) on whatever `_box()` / `_rod_axis_distance` currently read."""
    from tests.mesh.test_birdcage_implant_rod import ROD_RADIUS
    from tests.validation.test_wf7_f_human_sar10g_hotspot import CONTAINMENT_MARGIN_M

    centre, half = _box()
    a = _a10()
    bar = ROD_RADIUS + a + CONTAINMENT_MARGIN_M
    grid_all = _grid_report()
    kept = [r for r in grid_all if r[4]]
    excluded = [r for r in grid_all if not r[4]]
    lattice = [r[1] for r in _lattice_report()]
    match = _match_lattice([r[1] for r in kept], lattice)
    hit = sorted(mk for mk in match if mk is not None)
    lat_miss = max(min(float(np.max(np.abs(np.asarray(p) - np.asarray(r[1])))) for r in kept)
                   for p in lattice)
    wall_max = max(r[2] for r in kept)
    rod_min = min(r[3] for r in kept)
    # (G-c): the y = ±50 mm layer.
    layer = _grid_report(y_indices=(-1, 1))
    layer_kept = [r for r in layer if r[4]]
    layer_walls = sorted({round(r[2], 12) for r in layer})
    # (G-x): a uniform 50 mm x grid anchored at the centre.
    gx = _grid_report(x_step=0.05)
    gx_match = _match_lattice([r[1] for r in gx if r[4]], lattice)
    gx_missed = len(lattice) - len({mk for mk in gx_match if mk is not None})
    print(f"\n[GEO-38 C6b G1{tag}] box centre {tuple(centre)} half {tuple(half)}; x step "
          f"{2 * half[0] / 8 * 1e3:.4f} mm, y/z {C6B_GRID_DY_M * 1e3:.1f}/{C6B_GRID_DZ_M * 1e3:.1f} mm; "
          f"candidates {len(grid_all)}, kept {len(kept)}; excluded "
          f"{[(r[0], round(r[2], 6), round(r[3], 6)) for r in excluded]}")
    print(f"[GEO-38 C6b G1{tag}] lattice {len(lattice)} points, matched {len(set(hit))} (worst nearest "
          f"{lat_miss:.3e} m, tol {C6B_MATCH_M:g}); kept wall max {wall_max:+.4e} m (<= 0), rod-axis min "
          f"{rod_min:.4e} m (bar {bar:.4e})")
    print(f"[GEO-38 C6b G1{tag}] (G-c) y = +/-50 mm layer: {len(layer)} candidates, {len(layer_kept)} "
          f"kept, walls {layer_walls}; (G-x) 50 mm x grid: {len(gx)} candidates, lattice points "
          f"missed {gx_missed} of {len(lattice)}")
    return {"n_all": len(grid_all), "n_kept": len(kept), "excluded": excluded, "hit": hit,
            "lat_miss": lat_miss, "wall_max": wall_max, "rod_min": rod_min, "bar": bar,
            "layer": layer, "layer_kept": layer_kept, "gx_missed": gx_missed, "kept": kept}


def _g1_asserts(g):
    assert g["n_all"] == 91, g["n_all"]
    assert g["n_kept"] == 90, g["n_kept"]
    assert len(g["excluded"]) == 1
    ex = g["excluded"][0]
    assert ex[0] == "x+105.0mm y+0.0mm dz+0.0mm" and ex[3] == 0.0, ex
    assert g["hit"] == list(range(BRIEF_LATTICE_N)), g["hit"]  # the weekly's subset gate
    assert g["lat_miss"] <= C6B_MATCH_M, g["lat_miss"]
    assert g["wall_max"] <= 0.0, g["wall_max"]
    assert g["rod_min"] >= g["bar"] and round(g["bar"], 6) == 1.5365e-02, (g["rod_min"], g["bar"])
    # (G-c): wholly rejected, wall +1.9365e-02 m on each of the 182.
    assert len(g["layer"]) == 182 and not g["layer_kept"], (len(g["layer"]), len(g["layer_kept"]))
    assert all(abs(r[2] - 1.9365e-02) <= 5.0e-7 and r[2] > 0.0 for r in g["layer"])
    # (G-x): the x step is forced.
    assert g["gx_missed"] == 12, g["gx_missed"]


@_coding
def test_c6b_g1_grid_on_box_constants():
    """(G1) pure python on the real box: 91 candidates, 90 kept, the 18 lattice points
    a subset to 1e-12; controls (G-c) and (G-x)."""
    _g1_asserts(_g1_check(""))


@complex_only
@_coding
def test_c34_k8_d2_run_box_on_stubs(monkeypatch, tmp_path):
    """(K8-d2): `_run_box(…, degree=2)` on K8's stubs — C4a's variable reads "2" at
    the fixture call and is unset after; arm R's solver gets ``degree=2``; the JSON
    lands at the ``-d2`` name with K8's exact key set; (K8-B) 3 and (K8-S) 9
    unchanged.  Control (d2-c): a fixture reporting degree 1 under ``degree=2``
    raises at the degree assert, before any drive."""
    import json

    from tests.validation import test_geo36_rod_f_clinical_drive as i36
    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    comm = MPI.COMM_WORLD
    fx = _k6_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    good_dir, neg_dir = os.path.join(root, "d2"), os.path.join(root, "d2c")

    with monkeypatch.context() as m:
        m.delenv(C6_DEGREE_ENV, raising=False)
        m.delenv(c4.C4_DEGREE_ENV, raising=False)
        lines, _seen, _z = _k8_stubs(m, fx, good_dir)
        env_seen, degrees = _c34_recorders(m, degree=2)
        _run_box(comm, m, K8_NODE_PF, degree=2)
        env_after = os.environ.get(c4.C4_DEGREE_ENV)

    with monkeypatch.context() as m:
        m.delenv(C6_DEGREE_ENV, raising=False)
        _k8_stubs(m, fx, neg_dir)
        _env_c, degrees_c = _c34_recorders(m, degree=1)
        drives = []
        real_drive = i36._drive

        def _counted(*a, **k):
            drives.append(a[3])
            return real_drive(*a, **k)

        m.setattr(i36, "_drive", _counted)
        with pytest.raises(AssertionError) as d2c:
            _run_box(comm, m, K8_NODE_PF, degree=2)

    rec, d1_exists, neg_files = None, None, None
    if comm.rank == 0:
        with open(os.path.join(good_dir, f"c6-node-{K8_NODE_PF:.4f}-d2.json")) as fh:
            rec = json.load(fh)
        d1_exists = os.path.exists(os.path.join(good_dir, f"c6-node-{K8_NODE_PF:.4f}.json"))
        neg_files = sorted(os.listdir(neg_dir)) if os.path.isdir(neg_dir) else []
    rec, d1_exists, neg_files = comm.bcast((rec, d1_exists, neg_files), root=0)
    r_b = abs(rec["mean_b_plus"]["ccw"] / rec["mean_b_plus"]["cw"] / 3.0 - 1.0)
    r_s = abs(rec["sar_star"] / rec["sar_star_cw"] / 9.0 - 1.0)
    deg_lines = [ln for ln in lines if ln.startswith("(Deg) degree 2")]
    start_d2 = any(ln.lstrip().startswith("start:") and "degree 2," in ln for ln in lines)
    if comm.rank == 0:
        print(f"\n[GEO-38 C34 K8-d2] C4 variable at the fixture call {env_seen}, after {env_after!r}; "
              f"solver degrees {degrees}; JSON -d2 keys {sorted(rec) == sorted(K8_JSON_KEYS)}, "
              f"degree-1 name written {d1_exists}; start line degree 2 {start_d2}; "
              f"{deg_lines[0] if deg_lines else '(Deg) line MISSING'} (ASSERTED)")
        print(f"[GEO-38 C34 K8-d2] (K8-B) miss from 3 {r_b:.3e} (band {DISCRETE_IDENTITY_RTOL:g}); "
              f"(K8-S) miss from 9 {r_s:.3e} (band {EXACT_IDENTITY_RTOL:g}) (ASSERTED)")
        print(f"[GEO-38 C34 d2-c] degree 1 under degree=2 raised {str(d2c.value).splitlines()[0]}; "
              f"drives {drives}, solvers {degrees_c}, files {neg_files} (ASSERTED)", flush=True)
    assert env_seen == ["2"], env_seen
    assert env_after is None, env_after
    assert degrees == [2], degrees
    assert set(rec) == K8_JSON_KEYS, sorted(rec)
    assert not d1_exists
    assert start_d2
    assert len(deg_lines) == 1, deg_lines
    assert r_b <= DISCRETE_IDENTITY_RTOL, r_b
    assert r_s <= EXACT_IDENTITY_RTOL, r_s
    assert "(1, 2)" in str(d2c.value), str(d2c.value)
    assert drives == [] and degrees_c == [] and neg_files == [], (drives, degrees_c, neg_files)


@_coding
def test_c34_k62_entry62_is_two_calls(monkeypatch):
    """(K-62): C6a's foreign guard refuses I36's knobs in C6a's own environment, so
    ``xxl`` entry 62 runs C6a and I36 as two pytest calls, each with its own prefix.
    No guard is changed."""
    from tests.validation import test_geo36_rod_f_clinical_drive as i36

    for k in FOREIGN_ENVS + (i36.S8_DRIVE_C_PF_ENV, C6_DEGREE_ENV):
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv(C6_ENV, "1")
    monkeypatch.setenv(i36.S8_DRIVE_C_PF_ENV, f"{K8_NODE_PF:.4f}")
    monkeypatch.setenv(C6_DEGREE_ENV, "2")
    assert _knobs() == ("1", K8_NODE_PF) and _c6_degree("1") == 2  # C6a's call alone: accepted
    monkeypatch.setenv(i36.S8_DRIVE_ENV, "1")
    with pytest.raises(ValueError) as one:
        _knobs()
    monkeypatch.setenv(i36.S8_DRIVE_DEGREE_ENV, "2")
    with pytest.raises(ValueError) as two:
        _knobs()
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-38 C34 K-62] C6=1 + C6_DEGREE=2 + node accepted alone; with "
              f"{i36.S8_DRIVE_ENV}=1 refused: {one.value}; with {i36.S8_DRIVE_DEGREE_ENV}=2 too: "
              f"{two.value} — entry 62 is two pytest calls in one window (ASSERTED)", flush=True)
    assert i36.S8_DRIVE_ENV in str(one.value)
    assert i36.S8_DRIVE_ENV in str(two.value) and i36.S8_DRIVE_DEGREE_ENV in str(two.value)


def test_c34_deg_c_degree1_count_is_not_degree2():
    """(Deg-c): `_check_unknowns` accepts 2E + 2F = 8 327 394 at degree 2 and refuses
    the degree-1 count there; U₂/U₁ = 5.1828 is the ceiling, asserted > 5 on integers."""
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation.test_geo38_f_clinical_matrix import _check_unknowns

    u1, u2 = i35.F_CLINICAL_ROD_UNKNOWNS_RECORD
    n_edges, n_faces = 1_606_744, 2_556_953  # I36b's build window, :11989–11990
    assert _check_unknowns(u2, 2, n_edges, n_faces) == u2
    with pytest.raises(AssertionError):
        _check_unknowns(u1, 2, n_edges, n_faces)
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-38 C34 Deg-c] 2E + 2F = {2 * n_edges + 2 * n_faces} = record {u2}; U1 {u1} at "
              f"degree 2 refused; U2/U1 {u2 / u1:.4f} > 5 (ASSERTED on integers)", flush=True)
    assert u1 == n_edges
    assert u2 > 5 * u1


@complex_only
@_coding
def test_c6b_k8_grid_run_box_on_stubs(monkeypatch, tmp_path):
    """(K8-grid) `_run_box` on `_k8_stubs` with the knob on: 35 kept on K6's box, K8's
    axis rows a grid subset, grid max ≥ lattice max, g* at |dz| = 0.10 m, cw/ccw at g*
    = 1/9; the JSON keys K8's ∪ {"grid"}.  Control (G-off): the knob unset ⇒ no "grid"
    key, ``set(rec) == keys``."""
    comm = MPI.COMM_WORLD
    fx = _k6_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    on_dir, off_dir = os.path.join(root, "on"), os.path.join(root, "off")
    keys = {"c_pf", "feeds", "az_feeds", "s2", "y_in", "mean_b_plus", "mean_b_minus", "lattice",
            "c_star", "star_label", "sar_star", "sar_star_cw", "p_box", "p_src", "p_acc", "parts"}
    _k8_grid_check(comm, monkeypatch, fx, on_dir, off_dir, keys, "", lambda m: None,
                   lambda m, d: _k8_stubs(m, fx, d), K8_NODE_PF, _k8_json)


def _k8_grid_check(comm, monkeypatch, fx, on_dir, off_dir, keys, tag, patches, stubs, node, read_json):
    """Shared by (K8-grid) and C21's (K8-grid-3T): off run (G-off), then on run."""
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    module = __import__("sys").modules[_run_box.__module__]
    with monkeypatch.context() as m:
        m.delenv(C6B_GRID_ENV, raising=False)
        patches(m)
        off_lines = stubs(m, off_dir)[0]
        module._run_box(comm, m, node)
    rec_off = read_json(comm, off_dir)

    with monkeypatch.context() as m:
        m.setenv(C6B_GRID_ENV, "1")
        patches(m)
        lines = stubs(m, on_dir)[0]
        module._run_box(comm, m, node)
    rec = read_json(comm, on_dir)

    g = rec["grid"]
    lat_max = max(p["sar10g"] for p in rec["lattice"])
    star_dz = abs(g["c_star"][2])
    r_cw = abs(g["sar_star_cw"] / g["sar_star"] * 9.0 - 1.0)
    axis_hit = sorted(mk for mk in g["lattice_match"] if mk is not None)
    heads = ("(i) ball radius", "normalisers:", "price:", "(C6b) grid n")
    missing_heads = [h for h in heads if not any(ln.startswith(h) for ln in lines)]
    price = [ln for ln in lines if ln.startswith("price:")]
    off_c6b = [ln for ln in off_lines if ln.startswith("(C6b)")]
    if comm.rank == 0:
        print(f"\n[GEO-38 C6b K8-grid{tag}] n {g['n']} (35 expected); lattice rows matched {axis_hit}; g* "
              f"{g['star_label']} {g['c_star']}: SAR10g {g['sar_star']:.9e} vs lattice max {lat_max:.9e} "
              f"W/kg (>=, ASSERTED); |dz| {star_dz:.6f} m; cw/ccw x 9 miss {r_cw:.3e} (band "
              f"{EXACT_IDENTITY_RTOL:g}, ASSERTED)")
        print(f"[GEO-38 C6b K8-grid{tag}] keys {sorted(rec)}; price line {price}; missing heads "
              f"{missing_heads or 'none'}")
        print(f"[GEO-38 C6b K8-grid{tag}] control (G-off): keys {sorted(rec_off)} == K8's set "
              f"{set(rec_off) == keys}; (C6b) lines {len(off_c6b)} (ASSERTED)", flush=True)
    assert g["n"] == 35, g["n"]
    assert axis_hit == [0, 1, 2, 3], axis_hit
    assert g["sar_star"] >= lat_max, (g["sar_star"], lat_max)
    assert abs(star_dz - 0.10) <= 1.0e-12, g["c_star"]
    assert r_cw <= EXACT_IDENTITY_RTOL, r_cw
    assert set(rec) == keys | {"grid"}, sorted(rec)
    assert not missing_heads, missing_heads
    assert len(price) == 1 and "; grid " in price[0], price
    # (G-off)
    assert "grid" not in rec_off and set(rec_off) == keys, sorted(rec_off)
    assert not off_c6b, off_c6b
    return rec, rec_off
