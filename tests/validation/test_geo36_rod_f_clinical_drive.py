"""GEO-36 I36a (§9 item 593, §10 chain I36a): the rod fixture's two-port drive
module — C4a's two-port build route (`test_geo38_f_clinical_esr_root.py`) run
unchanged on I33's rod-in-ASTM-box F-clinical mesh by composition (I35's
no-edit route), build-only.

Brief: ``docs/planning/chunks/GEO-36.md`` "§9 item 593" and its review ruling
2026-10-22 03:00 (binds over the brief): heavy, daytime, slot-first; the ``1``
branch stops before any build (I36 is the 10-28 weekly's); the rod material on
the degree-2 problem is asserted as an exact map identity (M2).

Selected by ``FEM_EM_GEO36_S8_DRIVE`` only (never by ``-k``):

* unset — the gated test skips; the coding tests (K0), (K1), (K2), (N-route),
  (N-U) run (no gmsh);
* ``build`` — C4a's ``_run_build`` on the rod fixture through `_route_drive`
  (four monkeypatches, no environment variable): I35's (D-rod) / (M) / (N-M) /
  (C) / (U1) inside the fixture, C4a's (f) / (A0) / (iv) / (S) / (C) / (W) /
  (Az) / (U1) / (F) / (Z) / (Deg-2) against the rod records, then (M2);
* ``1`` — reserved for I36: ``NotImplementedError`` before any build;
* anything else — ``ValueError``.

I36b (§9 item 656, ``docs/planning/chunks/GEO-36.md`` "§9 item 656"): the
``1`` branch is now I36's body (`_run_drive`), superseding the "reserved"
line above.  It requires the node knob ``FEM_EM_GEO36_S8_DRIVE_C_PF`` (one C in
pF, ≤ 4 decimals; missing ⇒ ``ValueError`` before any gmsh; refused under
``build``).  Arm M: the rod mesh (C4a's fixture through `_route_drive`), tag
403 at σ_wire, r = 0.2 Ω ESR on all 64 gaps at C, feeds P33 / P41, one
factorisation (solve + held column).  Arm R: the same mesh with tag 403 mapped
to the gel (I41's ``replace`` route), a second factorisation, no build.
Gates (imported bands): (A0) / census / (U1) / (D-rod) inside the fixture,
(M2), (A0m), (Z), (P) + NEG on both arms, (P_RC), (R1).  Printed only: tip
1 g ratio Q, η_rod(T1), the rod's power share, Im Y_in.  ``1`` runs on ``xl``
only (entry 48); item 656 ran ``build``.

``FEM_EM_GEO36_S8_DRIVE_DEGREE``: unset ⇒ 1; ``2`` ⇒ 2; else ``ValueError``
(C4a's ``_degree_knob_c4`` semantics under its own name).  Refused under
``build``: the build proof covers both degrees by construction (Deg-2).

Never reuses I35's ``_route_column`` (it sets ``FEM_EM_GEO38_STEP1B``).  Module
objects are imported inside function bodies only, never ``test_*`` names.

Scope: the module, its coding tests and a build-only proof.  Not claimed: any
solve, drive, S entry, Y_in, root or tuning statement; a tip 1 g ratio, SAR or
heating; the no-rod twin (I36: tag 403 → the gel's material, never
``problem_gel_only``); I36's gates.  The ASTM box has no feed-bisector mirror
(its mirror planes are the coil-frame xz / yz planes; the bisector is at
50.625° and the rod at 0° breaks yz as well) — named for the weekly.
"""

from __future__ import annotations

import math
import os
import time

import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

S8_DRIVE_ENV = "FEM_EM_GEO36_S8_DRIVE"
S8_DRIVE_DEGREE_ENV = "FEM_EM_GEO36_S8_DRIVE_DEGREE"
#: I36b's node knob (one C in pF): required under ``1``, refused under ``build``.
S8_DRIVE_C_PF_ENV = "FEM_EM_GEO36_S8_DRIVE_C_PF"
#: I36's JSON (gitignored logs), one per node and degree.
I36_JSON_DIR = "/workspace/logs/geo36"
#: Variables an I36a window must never read (the I35 / C4 / step-1b / stop-rule knobs).
FOREIGN_ENVS = (
    "FEM_EM_GEO36_S8_COLUMN",
    "FEM_EM_GEO36_ROD_BOX",
    "FEM_EM_GEO38_STEP1B",
    "FEM_EM_GEO38_C4",
    "FEM_EM_GEO38_C4_NODES_PF",
    "FEM_EM_GEO38_C4_DEGREE",
    "FEM_EM_GMSH_STOP_RULE",
)
#: entry 40's 1.5 T r1 records, read back in (K1) / (N-U) (not new targets).
_R1_CELL_RECORD_EXPECTED = 1_147_366
_R1_UNKNOWNS_RECORD_EXPECTED = (1_525_193, 7_885_674)

# Predictions written first (brief, printed only, never asserted).
PREDICTED_ROD_AZ_DEG = 0.0
PREDICTED_BISECTOR_DEG = 50.625
PREDICTED_T_BUILD_S = (410.0, 560.0)
PREDICTED_WINDOW_S = 467.0
PREDICTED_RSS_GIB = 2.1


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-36 I36a] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Knobs.
# ---------------------------------------------------------------------------
def _guard_foreign():
    bad = [k for k in FOREIGN_ENVS if os.environ.get(k) is not None]
    if bad:
        raise ValueError(f"I36a never reads these variables; unset them: {bad}")


def _drive_mode():
    raw = os.environ.get(S8_DRIVE_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("build", "1"):
        return raw
    raise ValueError(f"{S8_DRIVE_ENV} must be unset, 'build' or '1', got {raw!r}")


def _drive_degree():
    """Unset / empty ⇒ ``None`` (degree 1); ``"2"`` ⇒ 2; anything else ⇒ ValueError."""
    raw = os.environ.get(S8_DRIVE_DEGREE_ENV, "").strip()
    if raw == "":
        return None
    if raw == "2":
        return 2
    raise ValueError(
        f"{S8_DRIVE_DEGREE_ENV}={raw!r}: only unset/empty (degree 1) or '2' are accepted"
    )


def _knobs():
    """Guard, then ``(mode, degree)``; a set degree knob under ``build`` is refused."""
    _guard_foreign()
    mode = _drive_mode()
    degree = _drive_degree()
    if mode == "build" and degree is not None:
        raise ValueError(
            f"{S8_DRIVE_DEGREE_ENV} is refused under {S8_DRIVE_ENV}=build: "
            f"the build proof covers both degrees (Deg-2)"
        )
    return mode, degree


# ---------------------------------------------------------------------------
# The route: C4a's build on the rod fixture (four patches, no env variable).
# ---------------------------------------------------------------------------
def _route_drive(monkeypatch):
    """Patch C4a's call-time names onto I35's rod fixture.  Returns the list the
    `_apply_degree_c4` recorder fills with ``(degree, material_map copy)``; the
    solver itself is never held (one solver, C4a)."""
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo38_f_clinical_build as s1a
    from tests.validation import test_geo38_f_clinical_column as col
    from tests.validation import test_geo38_f_clinical_esr_root as c4

    records = []
    original = c4._apply_degree_c4

    def _recording_apply_degree(ctx, comm):
        solver, degree, counts = original(ctx, comm)
        records.append((degree, dict(solver.problem.material_map)))
        return solver, degree, counts

    monkeypatch.setattr(col, "_build_f_clinical", i35._build_rod_f_clinical)
    monkeypatch.setattr(s1a, "F_CLINICAL_R1_CELL_RECORD", i35.F_CLINICAL_ROD_CELL_RECORD)
    monkeypatch.setattr(s1a, "F_CLINICAL_R1_UNKNOWNS_RECORD", i35.F_CLINICAL_ROD_UNKNOWNS_RECORD)
    monkeypatch.setattr(c4, "_apply_degree_c4", _recording_apply_degree)
    return records


def _m2_failures(records):
    """(M2): degrees [1, 2]; tag 403 σ = SIGMA_WIRE with the gel's ε_r / μ_r; tag 3
    σ = SALINE_SIGMA, ε_r = SALINE_EPSILON_R.  Exact equality of imported constants."""
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
    from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M

    fails = []
    degrees = [d for d, _m in records]
    if degrees != [1, 2]:
        fails.append(f"degrees {degrees} != [1, 2]")
    for d, mmap in records:
        rod = mmap.get(BIRDCAGE_IMPLANT_CELL_TAG)
        gel = mmap.get(i35.PHANTOM_CELL_TAG)
        if rod is None or gel is None:
            fails.append(f"degree {d}: tag 403 {rod!r} / tag 3 {gel!r} missing")
            continue
        if rod.sigma != SIGMA_WIRE_S_PER_M:
            fails.append(f"degree {d}: tag 403 sigma {rod.sigma!r} != {SIGMA_WIRE_S_PER_M}")
        if rod.epsilon_r != gel.epsilon_r or rod.mu_r != gel.mu_r:
            fails.append(
                f"degree {d}: tag 403 eps/mu {rod.epsilon_r!r}/{rod.mu_r!r} != gel "
                f"{gel.epsilon_r!r}/{gel.mu_r!r}"
            )
        if gel.sigma != SALINE_SIGMA or gel.epsilon_r != SALINE_EPSILON_R:
            fails.append(
                f"degree {d}: tag 3 sigma/eps {gel.sigma!r}/{gel.epsilon_r!r} != "
                f"{SALINE_SIGMA}/{SALINE_EPSILON_R}"
            )
    return fails


# ---------------------------------------------------------------------------
# I36b: the node knob and the body's pure-python pieces (exercised by the
# coding tests on small arrays).
# ---------------------------------------------------------------------------
def _node_knob(mode):
    """``S8_DRIVE_C_PF_ENV``: required under ``1`` (missing ⇒ ValueError before any
    gmsh), refused under ``build``, not read when the mode is unset.  One float
    > 0 that round-trips at 4 decimals (C4a's `_nodes_pf` naming rule)."""
    import numpy as np

    raw = os.environ.get(S8_DRIVE_C_PF_ENV, "").strip()
    if mode is None:
        return None
    if mode == "build":
        if raw != "":
            raise ValueError(f"{S8_DRIVE_C_PF_ENV} is refused under {S8_DRIVE_ENV}=build")
        return None
    if raw == "":
        raise ValueError(f"{S8_DRIVE_C_PF_ENV} is required under {S8_DRIVE_ENV}=1 (the node C, pF)")
    try:
        c_pf = float(raw)
    except ValueError as exc:
        raise ValueError(f"{S8_DRIVE_C_PF_ENV}={raw!r}: not a float") from exc
    if not np.isfinite(c_pf) or c_pf <= 0.0:
        raise ValueError(f"{S8_DRIVE_C_PF_ENV}={raw!r}: the node must be finite and > 0")
    if abs(float(f"{c_pf:.4f}") - c_pf) > 1e-12 * max(1.0, abs(c_pf)):
        raise ValueError(f"{S8_DRIVE_C_PF_ENV}={raw!r}: at most 4 decimals")
    return c_pf


def _i36_json_path(c_pf, degree):
    """``i36-node-<pF:.4f>.json`` at degree 1, a ``-d2`` sibling at degree 2."""
    if degree in (None, 1):
        return os.path.join(I36_JSON_DIR, f"i36-node-{c_pf:.4f}.json")
    if degree == 2:
        return os.path.join(I36_JSON_DIR, f"i36-node-{c_pf:.4f}-d2.json")
    raise ValueError(degree)


def _partition(p):
    """C4a's `_partition` closure (`test_geo38_f_clinical_esr_root.py:568–572`),
    module level: ``(rel, blind)`` of TH-14's identity P_src = ΣP_sheet + P_vol +
    P_surf, ``blind`` with P_surf dropped (NEG)."""
    p_vol = p["p_phantom"] + p["p_rest"]
    rel = abs(p["p_src"] - p["p_sheets"] - p_vol - p["p_surf"]) / abs(p["p_src"])
    blind = abs(p["p_src"] - p["p_sheets"] - p_vol) / abs(p["p_src"])
    return rel, blind


def _p_rc_split(p, feeds):
    """C4a's (P_RC) split (`:609–612`): ``(p_rc, p_feeds, dz)``, ``dz`` the
    identity's miss with the 62 RC sheets dropped (printed)."""
    p_rc = sum(v for pid, v in p["p_sheet_by_id"].items() if pid not in feeds)
    p_feeds = sum(v for pid, v in p["p_sheet_by_id"].items() if pid in feeds)
    p_vol = p["p_phantom"] + p["p_rest"]
    dz = abs(p["p_src"] - (p["p_sheets"] - p_rc) - p_vol - p["p_surf"]) / abs(p["p_src"])
    return p_rc, p_feeds, dz


def _gel_twin_map(mmap):
    """Arm R's map (I41's `_arm_with_403`, `test_imp3_rod_heating_matched.py:113–118`):
    ``{**mmap, 403: mmap[3]}``; the input is not mutated."""
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from tests.validation import test_geo36_rod_f_clinical as i35

    gel = mmap[i35.PHANTOM_CELL_TAG]
    out = {**mmap, BIRDCAGE_IMPLANT_CELL_TAG: gel}
    assert out[BIRDCAGE_IMPLANT_CELL_TAG] is out[i35.PHANTOM_CELL_TAG]
    return out


def _a0m_failures(census_m, census_r):
    """(A0m): arm R's tag-3 + tag-403 census equals arm M's, exactly."""
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from tests.validation import test_geo36_rod_f_clinical as i35

    fails = []
    for tag in (i35.PHANTOM_CELL_TAG, BIRDCAGE_IMPLANT_CELL_TAG):
        a, b = census_m.get(tag), census_r.get(tag)
        if a is None or b is None or a != b or a <= 0:
            fails.append(f"tag {tag}: arm M {a!r} vs arm R {b!r}")
    return fails


def _probe_points(centre):
    """``(tips, t1)`` in F-clinical's frame: tips at centre ± ROD_LENGTH/2 on ẑ (the
    brief's literal; F-small's `TIPS` sit at ±(L/2 + r), printed beside as
    ``apex``), T1 = centre + (F-small's T1 − F-small's ROD_CENTRE).  Never `TIPS`
    itself (F-small coordinates, `test_imp3_rod_birdcage.py:140`)."""
    import numpy as np

    from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE, ROD_LENGTH, ROD_RADIUS
    from tests.validation.test_imp3_palace_fields import PNAMES, PTS

    c = np.asarray(centre, dtype=float)
    half = 0.5 * ROD_LENGTH
    tips = {
        "top": tuple(c + (0.0, 0.0, half)),
        "bottom": tuple(c - (0.0, 0.0, half)),
        "top_apex": tuple(c + (0.0, 0.0, half + ROD_RADIUS)),
        "bottom_apex": tuple(c - (0.0, 0.0, half + ROD_RADIUS)),
    }
    t1 = c + (
        np.asarray(PTS[PNAMES.index("T1")], dtype=float) - np.asarray(ROD_CENTRE, dtype=float)
    )
    return tips, tuple(t1)


def _y_in_1p1j(s2, z0):
    """Active Y_in at both feeds under the (1, +1j) drive (H17's `_active_y_in`)."""
    from tests.validation import test_wf10_f_human_copper_esr_root as h17

    return h17._active_y_in(s2, h17.SENSES["(1,+1j)"], z0)


# ---------------------------------------------------------------------------
# I36b: the body's FEM pieces.  C4a's `_drive` / `_pair` closures
# (`test_geo38_f_clinical_esr_root.py:497–538`, `:551–563`) lifted to module
# level with their free variables passed in ``env``.
# ---------------------------------------------------------------------------
def _drive_env(fixture, comm):
    """The free variables of C4a's `_drive` (`:487–493`), built once per mesh."""
    import ufl

    from fem_em_solver.core.cavity import surface_resistance_ohm
    from tests.validation import test_geo38_f_clinical_matrix as c3
    from tests.validation.test_th14_birdcage_copper import _leontovich_term

    built = fixture["built"]
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, c3.COPPER_SIGMA)
    return {
        "comm": comm,
        "msh": msh,
        "tags_f": tags_f,
        "omega": omega,
        "cell_tags": cell_tags,
        "facets": built["facets"],
        "z_s": z_s,
        "wall": _leontovich_term(msh, built["facets"], z_s, omega),
        "dx": ufl.Measure("dx", domain=msh, subdomain_data=cell_tags),
        "acct": {"mesh": msh, "facet_tags": tags_f},
    }


def _drive(env, solver, specs, driven_id, *, held):
    """C4a's `_drive` (`:497–538`) with the partition always on, tag 403 split
    out of P_rest (``p_rod``, printed), and the field kept (``e``, ``sigma``)."""
    import numpy as np
    import ufl
    from dolfinx import fem

    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.validation.test_birdcage_power_identity import (
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_th14_birdcage_copper import _surface_loss_w

    comm, msh, tags_f, omega = env["comm"], env["msh"], env["tags_f"], env["omega"]
    sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in specs]
    dsheet = next(s for s in sheets if s.port_id == driven_id)
    linear = [
        lambda test, _s=dsheet: lumped_port_linear_term(
            msh, tags_f, _s, test, omega_rad_per_s=omega
        )
    ]
    comm.Barrier()
    t0 = time.perf_counter()
    if held:
        fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
    else:
        bil = [env["wall"]] + [
            (
                lambda trial, test, _s=s: lumped_port_bilinear_term(
                    msh, tags_f, _s, trial, test, omega_rad_per_s=omega
                )
            )
            for s in sheets
        ]
        fields = solver.solve(
            current_density=None,
            project_source=False,
            extra_bilinear_terms=bil,
            extra_linear_terms=linear,
        )
    comm.Barrier()
    t_solve = time.perf_counter() - t0
    t1 = time.perf_counter()
    e = fields.e_complex
    currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
    out = {"currents": currents, "kind": "held" if held else "solve", "t_solve": t_solve}
    sigma_f = fields.sigma_field
    dx = env["dx"]

    # F821: pyflakes sees the `del` at the end of this function; the closure runs before it.
    def _vol(measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)  # noqa: F821
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    by_id = {s.port_id: float(_sheet_field_dissipation_w(env["acct"], s, e, omega)) for s in sheets}
    out["p_sheet_by_id"] = by_id
    out["p_sheets"] = float(sum(by_id.values()))
    out["p_src"] = float(_source_power_w(env["acct"], dsheet, e, omega))
    p_total = _vol(dx)
    out["p_phantom"] = _vol(dx(PHANTOM_CELL_TAG))
    out["p_rod"] = _vol(dx(BIRDCAGE_IMPLANT_CELL_TAG))
    out["p_rest"] = p_total - out["p_phantom"]
    out["p_surf"] = _surface_loss_w(msh, env["facets"], env["z_s"], e, comm)
    out["e"] = e.copy()
    out["sigma"] = sigma_f
    del fields, e
    out["t_post"] = time.perf_counter() - t1
    return out


def _vol_powers(env, sigma_f, e):
    """``(P_total, P_gel, P_rod)`` of ½σ|E|², reduced (rank-safe)."""
    import numpy as np
    import ufl
    from dolfinx import fem

    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

    comm, dx = env["comm"], env["dx"]

    def _vol(measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    return _vol(dx), _vol(dx(PHANTOM_CELL_TAG)), _vol(dx(BIRDCAGE_IMPLANT_CELL_TAG))


def _arm(env, solver, label, specs, feeds, order, z_th, v_s, z0, rho_field, tips, t1):
    """One arm: C4a's `_pair` (`:551–563`, solve + held), S (current route), the
    partition on both columns, then the (1, +1j) superposition read at the tips
    (1 g ball, `_arm_d`'s `mass_averaged_sar` call, `test_imp3_tip_ratio_1g.py:207–213`)
    and T1 (`evaluate_vector_field_parallel`, owner ``lowest_original_cell``).
    Returns numbers only; the fields are released."""
    import gc

    import numpy as np

    from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
    from fem_em_solver.post.sar import averaging_ball_radius, mass_averaged_sar
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_mass_averaged_sar_standard_masses import (
        ONE_GRAM_KG,
        QUADRATURE_DEGREE,
    )

    comm = env["comm"]
    cur, kinds, times, parts, es = [], [], [], [], []
    sigma_f = None
    for n_drive, pid in enumerate(feeds):
        d = _drive(env, solver, specs, pid, held=(n_drive > 0))
        kinds.append(d["kind"])
        times.append(d["t_solve"] + d["t_post"])
        cur.append(np.array([d["currents"][p] for p in order], dtype=np.complex128))
        parts.append({k: v for k, v in d.items() if k.startswith("p_")})
        es.append(d["e"])
        sigma_f = d["sigma"]
        del d
    s2 = h17._current_route_s(cur, [0, 1], z_th, v_s)
    w = h17.SENSES["(1,+1j)"]
    e_sup = es[0].copy()
    e_sup.x.array[:] = w[0] * es[0].x.array + w[1] * es[1].x.array
    e_sup.x.scatter_forward()
    p_tot, p_gel, p_rod = _vol_powers(env, sigma_f, e_sup)
    a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    sar, mass = {}, {}
    for name, centre in tips.items():
        o = mass_averaged_sar(
            e_sup,
            sigma=sigma_f,
            rho=rho_field,
            center=centre,
            radius=a_1g,
            quadrature_degree=QUADRATURE_DEGREE,
        )
        sar[name] = float(o["averaged_sar_w_per_kg"])
        mass[name] = float(o["mass_kg"])
    pts = np.ascontiguousarray(np.array([t1], dtype=np.float64))
    vals, mask = evaluate_vector_field_parallel(e_sup, pts, owner="lowest_original_cell")
    vals = comm.bcast(vals if comm.rank == 0 else None, root=0)
    mask = comm.bcast(mask if comm.rank == 0 else None, root=0)
    e_t1 = np.asarray(vals).reshape(1, 3)[0]
    out = {
        "label": label,
        "s2": s2,
        "kinds": kinds,
        "times": times,
        "parts": parts,
        "p_sup": {"total": p_tot, "gel": p_gel, "rod": p_rod},
        "sar": sar,
        "mass": mass,
        "a_1g": a_1g,
        "e_t1": e_t1,
        "t1_valid": bool(np.all(mask)),
        "y_in": _y_in_1p1j(s2, z0),
    }
    del es, e_sup, sigma_f
    gc.collect()
    return out


def _run_drive(comm, monkeypatch, degree, c_pf):
    """I36's body (the ``1`` branch): arm M then arm R on one mesh, gates asserted
    at the end (every reading printed first)."""
    import dataclasses
    import gc
    import json

    import numpy as np

    from fem_em_solver.core import TimeHarmonicSolver
    from fem_em_solver.post.sar import build_density_field
    from tests.mesh.helpers import global_cell_tag_set
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo38_f_clinical_column as col
    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_geo38_f_clinical_matrix import BLIND_FACTOR
    from tests.validation.test_imp3_tip_ladder_1g import ENTRY43_READOUT_INTERVAL
    from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    t_window = time.perf_counter()
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    deg = degree or 1
    _say(
        comm,
        f"\nstart: -n {comm.size}, I36 body; node C {c_pf:.4f} pF, degree {deg}, "
        f"r = {c4.ESR_OHM} Ohm on all 64 gaps, kappa off; arm M (rod, tag 403 sigma_wire) "
        f"then arm R (tag 403 -> gel, same mesh)",
    )
    records = _route_drive(monkeypatch)
    if deg == 2:
        monkeypatch.setenv(c4.C4_DEGREE_ENV, "2")  # C4a's fixture reads its own name
    else:
        monkeypatch.delenv(c4.C4_DEGREE_ENV, raising=False)
    fixture = c4._build_c4_fixture(comm)  # (A0) / census / (U1) / (D-rod) / (F) asserted inside
    monkeypatch.delenv(c4.C4_DEGREE_ENV, raising=False)
    assert fixture["degree"] == deg, (fixture["degree"], deg)
    assert len(records) == 1 and records[0][0] == deg, [r[0] for r in records]
    # (M2): the map is degree-independent; the degree is asserted on the line above.
    m2_fails = _m2_failures([(1, records[0][1]), (2, records[0][1])])
    feeds, ids = fixture["feeds"], fixture["ids"]
    ctx = fixture["built"]["ctx"]
    env = _drive_env(fixture, comm)
    frequency_hz = env["omega"] / (2.0 * np.pi)
    others = sorted((pid for pid in ids if pid not in feeds), key=lambda p: int(p[1:]))
    order = list(feeds) + others

    # ---- (Z) the ESR specs at the node -----------------------------------
    c_f = c_pf * 1e-12
    base_specs = list(ctx["specs"])
    base_ids = [sp.port_id for sp in base_specs]
    specs, z_rc, z_th = c4._esr_specs(base_specs, feeds, c_f, frequency_hz, z0)
    z_fails, im_worst, th_worst = c4._z_failures(specs, base_ids, feeds, c_f, frequency_hz, z0)
    v_s = h17._source_voltage(z_rc, z0)
    _say(
        comm,
        f"(M2) failures {m2_fails or 'none'} (ASSERTED exact); (Z) node {c_pf:.4f} pF: Z_RC "
        f"{z_rc:.9e} Ohm, Z_th {z_th:.9e} Ohm, V_s {v_s:.9e} V; Im rel {im_worst:.1e}, Z_th "
        f"rel {th_worst:.1e}; failures {z_fails or 'none'} (ASSERTED)",
    )

    # ---- shared readers ----------------------------------------------------
    msh, cell_tags = env["msh"], env["cell_tags"]
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh,
        PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != i35.PHANTOM_CELL_TAG},
    )
    tips, t1 = _probe_points(i35._rod_centre())

    # ---- arm M -------------------------------------------------------------
    solver_m = ctx["solver"]
    problem_m = solver_m.problem
    census_m = col._owned_tag_counts(msh, problem_m.cell_tags, comm)
    arm_m = _arm(env, solver_m, "M (rod)", specs, feeds, order, z_th, v_s, z0, rho_field, tips, t1)
    del solver_m
    del ctx["solver"]
    gc.collect()

    # ---- arm R: same mesh, tag 403 -> gel, a second factorisation -----------
    problem_r = dataclasses.replace(
        problem_m, material_map=_gel_twin_map(dict(problem_m.material_map))
    )
    census_r = col._owned_tag_counts(msh, problem_r.cell_tags, comm)
    a0m_fails = _a0m_failures(census_m, census_r)
    solver_r = TimeHarmonicSolver(problem_r, degree=deg)
    arm_r = _arm(
        env, solver_r, "R (gel twin)", specs, feeds, order, z_th, v_s, z0, rho_field, tips, t1
    )
    del solver_r
    gc.collect()

    # ---- readings ----------------------------------------------------------
    fails = []
    for arm in (arm_m, arm_r):
        recip = _reciprocity_ratio(arm["s2"])
        arm["recip"] = recip
        _say(
            comm,
            f"arm {arm['label']}: kinds {arm['kinds']}, times "
            f"{', '.join(f'{t:.2f}' for t in arm['times'])} s; S2 "
            + "; ".join(", ".join(f"{complex(v):.10e}" for v in row) for row in arm["s2"])
            + f"; (R1) ||S-S^T||_F/||S||_F {recip:.3e} (band {RECIPROCITY_BAND:.0e}, ASSERTED)",
        )
        if arm["kinds"] != ["solve", "held"]:
            fails.append(f"{arm['label']} kinds {arm['kinds']}")
        if not recip <= RECIPROCITY_BAND:
            fails.append(f"{arm['label']} (R1) {recip:.3e}")
        for j, p in enumerate(arm["parts"]):
            rel, blind = _partition(p)
            _say(
                comm,
                f"  (P) column {feeds[j]}: P_src {p['p_src']:.9e} sumP_sheet {p['p_sheets']:.9e} "
                f"P_phantom {p['p_phantom']:.9e} P_rod {p['p_rod']:.9e} P_rest {p['p_rest']:.9e} "
                f"P_surf {p['p_surf']:.9e} W; rel {rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, "
                f"ASSERTED); NEG blind {blind:.4e} (>= {BLIND_FACTOR:g}x band, ASSERTED)",
            )
            if not rel <= DISCRETE_IDENTITY_RTOL:
                fails.append(f"{arm['label']} (P) {feeds[j]} {rel:.3e}")
            if not blind >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL:
                fails.append(f"{arm['label']} NEG {feeds[j]} {blind:.3e}")
            if not p["p_surf"] > 0.0:
                fails.append(f"{arm['label']} P_surf {p['p_surf']!r}")
        ya, yb = arm["y_in"]
        ps = arm["p_sup"]
        _say(
            comm,
            f"  Y_in (1,+1j): {feeds[0]} {ya.real:+.9e} {ya.imag:+.9e}j S, {feeds[1]} "
            f"{yb.real:+.9e} {yb.imag:+.9e}j S (Im printed: the coil was tuned on C4's "
            f"cylinder; never a red); (1,+1j) P_vol {ps['total']:.9e} W, gel {ps['gel']:.9e}, "
            f"rod {ps['rod']:.9e} (rod share {ps['rod'] / ps['total']:.6e}, printed)",
        )
    pm = arm_m["parts"][0]
    rel_prc, _b = _partition(pm)
    p_rc, p_feeds, dz = _p_rc_split(pm, feeds)
    _say(
        comm,
        f"(P_RC) arm M column {feeds[0]} at the node: rel {rel_prc:.3e} (band "
        f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED); 62 RC {p_rc:.9e} W, 2 feeds {p_feeds:.9e} W; "
        f"(Dz) without the RC sheets {dz:.4e} (printed)",
    )
    _say(
        comm,
        f"(A0m) census M {census_m} vs R {census_r}: failures {a0m_fails or 'none'} (ASSERTED "
        f"exact)",
    )
    q = {n: arm_m["sar"][n] / arm_r["sar"][n] for n in tips}
    eta_t1 = float(np.linalg.norm(arm_m["e_t1"]) / np.linalg.norm(arm_r["e_t1"]))
    _say(
        comm,
        "tip 1 g ratio Q = SAR_1g(M)/SAR_1g(R), (1,+1j), degree "
        f"{deg}: "
        + ", ".join(
            f"{n} {q[n]:.6f} (M {arm_m['sar'][n]:.6e}, R "
            f"{arm_r['sar'][n]:.6e} W/kg, mass {arm_m['mass'][n]:.3e} kg)"
            for n in tips
        )
        + f" beside ENTRY43_READOUT_INTERVAL {ENTRY43_READOUT_INTERVAL} (printed); "
        f"eta_rod(T1 = {tuple(round(x, 6) for x in t1)}) = |E_M|/|E_R| {eta_t1:.6f} (valid "
        f"{arm_m['t1_valid'] and arm_r['t1_valid']}; printed)",
    )
    if comm.rank == 0:
        os.makedirs(I36_JSON_DIR, exist_ok=True)
        path = _i36_json_path(c_pf, deg)
        rec = {
            "c_pf": c_pf,
            "degree": deg,
            "feeds": list(feeds),
            "q": q,
            "eta_t1": eta_t1,
            "census_m": {str(k): v for k, v in census_m.items()},
        }
        for arm, key in ((arm_m, "M"), (arm_r, "R")):
            rec[key] = {
                "s2": [[[complex(v).real, complex(v).imag] for v in row] for row in arm["s2"]],
                "y_in": [[y.real, y.imag] for y in arm["y_in"]],
                "recip": arm["recip"],
                "parts": [
                    {k: v for k, v in p.items() if k != "p_sheet_by_id"} for p in arm["parts"]
                ],
                "p_sup": arm["p_sup"],
                "sar": arm["sar"],
                "times": arm["times"],
            }
        with open(path, "w") as fh:
            json.dump(rec, fh, indent=1)
        print(f"[GEO-36 I36a] wrote {path}", flush=True)
    comm.Barrier()
    _m, per, tot = _rss_gib(comm)
    _say(
        comm,
        f"price: build {fixture['t_build']:.1f} s; arm M {sum(arm_m['times']):.1f} s, arm R "
        f"{sum(arm_r['times']):.1f} s; RSS per rank {', '.join(f'{r:.2f}' for r in (per or []))} "
        f"GiB, summed {tot:.3f} GiB; window {time.perf_counter() - t_window:.1f} s",
    )
    assert not m2_fails, m2_fails
    assert not z_fails, z_fails
    assert not a0m_fails, a0m_fails
    assert rel_prc <= DISCRETE_IDENTITY_RTOL, f"(P_RC) {rel_prc:.3e} — never re-band"
    assert not fails, fails


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
@complex_only
def test_geo36_i36_rod_drive(monkeypatch):
    mode, _degree = _knobs()
    c_pf = _node_knob(mode)  # I36b: missing under 1 / set under build => ValueError, before gmsh
    if mode is None:
        pytest.skip(f"{S8_DRIVE_ENV} unset: GEO-36 I36a is selected by environment only")
    if mode == "1":
        _run_drive(MPI.COMM_WORLD, monkeypatch, _degree, c_pf)
        return

    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()
    _say(
        comm,
        f"\nstart: -n {comm.size}, BUILD-ONLY; C4a's _run_build on the rod fixture "
        f"(records cells {i35.F_CLINICAL_ROD_CELL_RECORD}, U "
        f"{i35.F_CLINICAL_ROD_UNKNOWNS_RECORD}); PREDICTED t_build "
        f"{PREDICTED_T_BUILD_S[0]:.0f}-{PREDICTED_T_BUILD_S[1]:.0f} s (printed)",
    )
    records = _route_drive(monkeypatch)
    c4._run_build(comm, monkeypatch)

    m2_fails = _m2_failures(records)
    _say(
        comm,
        "(M2) recorded _apply_degree_c4 calls: "
        + "; ".join(f"degree {d}: tag 403 {m.get(403)!r}, tag 3 {m.get(3)!r}" for d, m in records)
        + f"; failures {m2_fails or 'none'} (ASSERTED exact)",
    )
    x, y, _z = i35._rod_centre()
    rod_az = math.degrees(math.atan2(y, x))
    _say(
        comm,
        f"rod azimuth atan2(y, x) of {i35._rod_centre()} = {rod_az:.3f} deg (PREDICTED "
        f"{PREDICTED_ROD_AZ_DEG}); bisector: see C4a's (F) line above (PREDICTED "
        f"{PREDICTED_BISECTOR_DEG}); rod - bisector {rod_az - PREDICTED_BISECTOR_DEG:.3f} "
        f"deg vs the predicted bisector (printed)",
    )
    _m, per, tot = _rss_gib(comm)
    _say(
        comm,
        f"BUILD-ONLY return: window {time.perf_counter() - t_window:.1f} s (PREDICTED "
        f"≈{PREDICTED_WINDOW_S:.0f} s); RSS per rank "
        f"{', '.join(f'{r:.2f}' for r in (per or []))} GiB, summed {tot:.3f} GiB "
        f"(PREDICTED ≈{PREDICTED_RSS_GIB}; printed)",
    )
    assert not m2_fails, m2_fails


# ---------------------------------------------------------------------------
# Coding tests (no gmsh); skip whenever S8_DRIVE_ENV is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(
    os.environ.get(S8_DRIVE_ENV, "").strip() != "",
    reason=f"{S8_DRIVE_ENV} set: coding tests skip in build windows",
)


class _Routed(Exception):
    def __init__(self, kwargs):
        super().__init__("routed")
        self.kwargs = kwargs


def _install_recorder(m):
    from fem_em_solver.io.mesh import MeshGenerator

    def _rec(*args, **kwargs):
        raise _Routed(kwargs)

    m.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(_rec))


@_coding
def test_i36a_k0_knob_parsing(monkeypatch):
    """(K0): drive / degree knobs, degree under build refused, foreign guard."""
    for k in (S8_DRIVE_ENV, S8_DRIVE_DEGREE_ENV) + FOREIGN_ENVS:
        monkeypatch.delenv(k, raising=False)
    assert _knobs() == (None, None)
    for ok in ("build", "1"):
        monkeypatch.setenv(S8_DRIVE_ENV, ok)
        assert _drive_mode() == ok
    for bad in ("2", "x"):
        monkeypatch.setenv(S8_DRIVE_ENV, bad)
        with pytest.raises(ValueError):
            _drive_mode()
    monkeypatch.delenv(S8_DRIVE_ENV)
    assert _drive_degree() is None
    monkeypatch.setenv(S8_DRIVE_DEGREE_ENV, "2")
    assert _drive_degree() == 2
    assert _knobs() == (None, 2)
    monkeypatch.setenv(S8_DRIVE_ENV, "1")
    assert _knobs() == ("1", 2)
    monkeypatch.setenv(S8_DRIVE_ENV, "build")
    with pytest.raises(ValueError):
        _knobs()
    for bad in ("1", "x"):
        monkeypatch.setenv(S8_DRIVE_DEGREE_ENV, bad)
        with pytest.raises(ValueError):
            _drive_degree()
    monkeypatch.delenv(S8_DRIVE_DEGREE_ENV)
    assert _knobs() == ("build", None)
    for k in FOREIGN_ENVS:
        with monkeypatch.context() as m:
            m.setenv(k, "1")
            with pytest.raises(ValueError):
                _knobs()
    print(
        "\n[GEO-36 I36a K0] drive unset/build/1 accepted, 2/x refused; degree unset/2 "
        f"accepted, 1/x refused; degree under build refused; {len(FOREIGN_ENVS)} foreign "
        "variables refused"
    )


@_coding
def test_i36a_k1_routing_tripwire(monkeypatch):
    """(K1): routed C4a fixture reaches gmsh with the kept rod and the box; all
    four patches restored on exit; no env variable set."""
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo38_f_clinical_build as s1a
    from tests.validation import test_geo38_f_clinical_column as col
    from tests.validation import test_geo38_f_clinical_esr_root as c4

    assert os.environ.get("FEM_EM_GEO38_STEP1B") is None
    orig_build, orig_apply = col._build_f_clinical, c4._apply_degree_c4
    with monkeypatch.context() as m:
        _route_drive(m)
        _install_recorder(m)
        assert col._build_f_clinical is i35._build_rod_f_clinical
        assert s1a.F_CLINICAL_R1_CELL_RECORD == i35.F_CLINICAL_ROD_CELL_RECORD
        assert s1a.F_CLINICAL_R1_UNKNOWNS_RECORD == i35.F_CLINICAL_ROD_UNKNOWNS_RECORD
        assert c4._apply_degree_c4 is not orig_apply
        with pytest.raises(_Routed) as exc:
            c4._build_c4_fixture(MPI.COMM_WORLD)
        kw = exc.value.kwargs
    rod = kw.get("implant_rod")
    print(
        f"\n[GEO-36 I36a K1] routed kwargs: implant_rod {rod}; phantom_box {kw.get('phantom_box')}"
    )
    assert rod is not None and rod["keep_cells"] is True, rod
    assert tuple(rod["centre"]) == tuple(i35._rod_centre()), rod
    assert kw.get("phantom_box") is not None, sorted(kw)
    assert col._build_f_clinical is orig_build
    assert c4._apply_degree_c4 is orig_apply
    assert s1a.F_CLINICAL_R1_CELL_RECORD == _R1_CELL_RECORD_EXPECTED
    assert s1a.F_CLINICAL_R1_UNKNOWNS_RECORD == _R1_UNKNOWNS_RECORD_EXPECTED
    assert os.environ.get("FEM_EM_GEO38_STEP1B") is None


@_coding
def test_i36a_k2_stub_stops_before_gmsh(monkeypatch):
    """(K2), I36b meaning: S8_DRIVE=1 **without** the node knob raises ValueError,
    never the gmsh sentinel (still a stop-before-gmsh test)."""
    with monkeypatch.context() as m:
        for k in FOREIGN_ENVS + (S8_DRIVE_DEGREE_ENV, S8_DRIVE_C_PF_ENV):
            m.delenv(k, raising=False)
        m.setenv(S8_DRIVE_ENV, "1")
        _install_recorder(m)
        with pytest.raises(ValueError, match=S8_DRIVE_C_PF_ENV):
            test_geo36_i36_rod_drive(m)
    print(
        "\n[GEO-36 I36a K2] S8_DRIVE=1 without the node knob raised ValueError before gmsh "
        "(ASSERTED)"
    )


@_coding
def test_i36b_k3_node_knob_parsing(monkeypatch):
    """(K3): the node knob — required under 1, refused under build, unread when
    unset; float > 0 with ≤ 4 decimals; build + knob stops before gmsh."""
    monkeypatch.delenv(S8_DRIVE_C_PF_ENV, raising=False)
    assert _node_knob(None) is None
    assert _node_knob("build") is None
    with pytest.raises(ValueError, match="required"):
        _node_knob("1")
    for ok, want in (("193.5178", 193.5178), ("115.5", 115.5), (" 100 ", 100.0)):
        monkeypatch.setenv(S8_DRIVE_C_PF_ENV, ok)
        assert _node_knob("1") == want
        assert _node_knob(None) is None
        with pytest.raises(ValueError, match="refused"):
            _node_knob("build")
    for bad in ("193.51781", "x", "0", "-5", "nan", "inf", "1,2"):
        monkeypatch.setenv(S8_DRIVE_C_PF_ENV, bad)
        with pytest.raises(ValueError):
            _node_knob("1")
    assert _i36_json_path(193.5178, None).endswith("i36-node-193.5178.json")
    assert _i36_json_path(193.5178, 2).endswith("i36-node-193.5178-d2.json")
    with pytest.raises(ValueError):
        _i36_json_path(1.0, 3)
    with monkeypatch.context() as m:
        for k in FOREIGN_ENVS + (S8_DRIVE_DEGREE_ENV,):
            m.delenv(k, raising=False)
        m.setenv(S8_DRIVE_ENV, "build")
        m.setenv(S8_DRIVE_C_PF_ENV, "193.5178")
        _install_recorder(m)
        with pytest.raises(ValueError, match="refused"):
            test_geo36_i36_rod_drive(m)
    print(
        "\n[GEO-36 I36b K3] node knob: required under 1, refused under build (gated test "
        "stops before gmsh), unread when unset; 3 accepted, 7 refused (ASSERTED)"
    )


@_coding
def test_i36b_k4_partition_and_p_rc_arithmetic():
    """(K4): `_partition` / `_p_rc_split` on a synthetic balance with known misses."""
    feeds = ["P33", "P41"]
    by_id = {"P33": 0.5, "P41": 0.25, "P34": 1.0, "P35": 2.0}
    p = {
        "p_src": 10.0,
        "p_sheets": 3.75,
        "p_sheet_by_id": by_id,
        "p_phantom": 4.0,
        "p_rest": 1.25,
        "p_surf": 1.0,
    }
    rel, blind = _partition(p)
    assert rel == 0.0, rel
    assert abs(blind - 0.1) < 1e-15, blind
    p_rc, p_feeds, dz = _p_rc_split(p, feeds)
    assert (p_rc, p_feeds) == (3.0, 0.75)
    assert abs(dz - 0.3) < 1e-15, dz
    q = dict(p, p_src=10.0 * (1 + 2e-6))
    rel_q, _ = _partition(q)
    assert abs(rel_q - 2e-6 / (1 + 2e-6)) < 1e-15, rel_q
    print(
        f"\n[GEO-36 I36b K4] partition rel 0 / blind {blind:.3e} / dz {dz:.3e} / planted "
        f"{rel_q:.6e} (ASSERTED)"
    )


@_coding
def test_i36b_k5_r1_and_current_route_on_synthetic_2x2():
    """(K5): H17's current route on a reciprocal 4-port circuit (unit Thevenin
    drives through Z_th, RC loads) gives a symmetric 2×2 — (R1) ratio ≤ 1e-12 —
    and a planted non-reciprocal Z fails RECIPROCITY_BAND; Y_in is read."""
    import numpy as np

    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    z0 = 50.0
    rng = np.random.default_rng(656)
    a = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    z_sym = a + a.T + np.eye(4) * (20.0 + 5.0j)
    z_rc = h17._z_rc(193.5178e-12, 64.0e6)
    z_th = h17._thevenin_feed_impedance(z_rc, z0)
    v_s = h17._source_voltage(z_rc, z0)

    def _s(zz):
        loads = h17._parallel_loads(4, [0, 1], z_rc, z_th)
        cur = []
        for j in (0, 1):
            e = np.zeros(4, dtype=np.complex128)
            e[j] = 1.0
            cur.append(h17._parallel_currents(zz, loads, e))
        return h17._current_route_s(cur, [0, 1], z_th, v_s)

    s_sym = _s(z_sym)
    r_sym = _reciprocity_ratio(s_sym)
    z_bad = z_sym.copy()
    z_bad[0, 1] *= 1.5
    r_bad = _reciprocity_ratio(_s(z_bad))
    y = _y_in_1p1j(s_sym, z0)
    print(
        f"\n[GEO-36 I36b K5] (R1) reciprocal {r_sym:.3e} <= 1e-12; planted {r_bad:.3e} > "
        f"{RECIPROCITY_BAND:.0e}; Y_in {y[0]:.6e}, {y[1]:.6e} (ASSERTED)"
    )
    assert r_sym <= 1e-12, r_sym
    assert r_bad > RECIPROCITY_BAND, r_bad
    assert len(y) == 2 and all(np.isfinite(complex(v)) for v in y)
    b = s_sym @ np.array(h17.SENSES["(1,+1j)"])
    assert abs(y[0] - (1 - b[0]) / (z0 * (1 + b[0]))) < 1e-15 * max(1.0, abs(y[0]))


@_coding
def test_i36b_k6_probe_points_in_f_clinical_frame():
    """(K6): tips at centre ± ROD_LENGTH/2 on ẑ, apexes at ±(L/2 + r), T1 at
    F-small's offset from its rod centre — never F-small's coordinates."""
    import numpy as np

    from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE, ROD_LENGTH, ROD_RADIUS
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation.test_imp3_rod_birdcage import TIPS

    c = np.asarray(i35._rod_centre())
    tips, t1 = _probe_points(tuple(c))
    assert np.allclose(np.asarray(tips["top"]) - c, (0, 0, ROD_LENGTH / 2), atol=1e-15, rtol=0)
    assert np.allclose(np.asarray(tips["bottom"]) - c, (0, 0, -ROD_LENGTH / 2), atol=1e-15, rtol=0)
    apex = np.asarray(tips["top_apex"]) - c
    assert np.allclose(apex, np.asarray(TIPS["top"]) - np.asarray(ROD_CENTRE), atol=1e-15, rtol=0)
    assert abs(apex[2] - (ROD_LENGTH / 2 + ROD_RADIUS)) < 1e-15
    d = np.asarray(t1) - c
    assert abs(d[2] - 0.017531) < 1e-12 and abs(d[0] - 0.000023) < 1e-12, d
    assert not np.allclose(np.asarray(tips["top"]), np.asarray(TIPS["top"]))
    print(
        f"\n[GEO-36 I36b K6] centre {tuple(c)}; top {tips['top']}; apex {tips['top_apex']}; "
        f"T1 {t1} (ASSERTED)"
    )


@_coding
def test_i36b_k7_gel_twin_map_and_a0m():
    """(K7): arm R's map is `{**mmap, 403: mmap[3]}` (input untouched, 403 is 3);
    (A0m) passes on equal censuses and fails on a planted one-cell difference."""
    from fem_em_solver.core import HomogeneousMaterial

    gel = HomogeneousMaterial(sigma=0.5, epsilon_r=78.0, mu_r=1.0)
    rod = HomogeneousMaterial(sigma=1.0e6, epsilon_r=78.0, mu_r=1.0)
    mmap = {3: gel, 403: rod, 2: HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0)}
    twin = _gel_twin_map(mmap)
    assert twin[403] is gel and mmap[403] is rod and twin[2] is mmap[2]
    assert set(twin) == set(mmap)
    census = {3: 1000, 403: 12, 2: 5}
    assert _a0m_failures(census, dict(census)) == []
    assert _a0m_failures(census, {**census, 403: 11})
    assert _a0m_failures(census, {3: 1000})
    print(
        "\n[GEO-36 I36b K7] gel twin map 403 -> gel object; (A0m) equal passes, planted "
        "differences fail (ASSERTED)"
    )


@_coding
def test_i36a_n_route_unpatched_reaches_15t_fixture(monkeypatch):
    """(N-route), asserted binary: C4a's fixture without `_route_drive` reaches
    the 1.5 T stand-in kwargs, with no rod."""
    from tests.validation import test_geo38_f_clinical_esr_root as c4

    with monkeypatch.context() as m:
        _install_recorder(m)
        with pytest.raises(_Routed) as exc:
            c4._build_c4_fixture(MPI.COMM_WORLD)
        kw = exc.value.kwargs
    print(f"\n[GEO-36 I36a N-route] unpatched kwargs implant_rod: {kw.get('implant_rod')}")
    assert kw.get("implant_rod") is None, kw.get("implant_rod")


@_coding
def test_i36a_n_u_15t_unknowns_miss_rod_mesh():
    """(N-U), asserted by construction: the 1.5 T unknown records miss the rod mesh,
    so an unbound record patch turns the exact (U1) / (Deg-2) red."""
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo38_f_clinical_build as s1a

    assert s1a.F_CLINICAL_R1_UNKNOWNS_RECORD == _R1_UNKNOWNS_RECORD_EXPECTED
    rod_u = i35.F_CLINICAL_ROD_UNKNOWNS_RECORD
    miss = abs(rod_u[0] / s1a.F_CLINICAL_R1_UNKNOWNS_RECORD[0] - 1.0)
    print(
        f"\n[GEO-36 I36a N-U] |{rod_u[0]}/{s1a.F_CLINICAL_R1_UNKNOWNS_RECORD[0]} - 1| = "
        f"{miss:.3e} (PREDICTED 5.347e-02); U2 {rod_u[1]} vs "
        f"{s1a.F_CLINICAL_R1_UNKNOWNS_RECORD[1]}"
    )
    assert rod_u[0] != s1a.F_CLINICAL_R1_UNKNOWNS_RECORD[0]
    assert rod_u[1] != s1a.F_CLINICAL_R1_UNKNOWNS_RECORD[1]


def _k8_arm(monkeypatch, column_fields):
    """(K8)'s stubbed `_arm` on the small mesh: ``column_fields(E_u) -> (e₁, e₂)``.
    Returns ``(arm, msh, cell_tags, E_u, sigma_f, rho, E0, tips, t1, a_1g)``."""
    import sys

    import dolfinx
    import numpy as np
    import ufl
    from dolfinx import fem

    from fem_em_solver.post.sar import averaging_ball_radius, build_density_field
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_mass_averaged_sar_standard_masses import ONE_GRAM_KG

    comm = MPI.COMM_WORLD
    msh = dolfinx.mesh.create_box(
        comm, [[-0.05] * 3, [0.05] * 3], [12, 12, 12], dolfinx.mesh.CellType.tetrahedron
    )
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    cells = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
    mids = dolfinx.mesh.compute_midpoints(msh, tdim, cells)
    values = np.where(np.max(np.abs(mids), axis=1) < 0.012, 403, 3).astype(np.int32)
    cell_tags = dolfinx.mesh.meshtags(msh, tdim, cells, values)
    E0 = 3.0
    e_u = fem.Function(fem.functionspace(msh, ("N1curl", 1)))
    e_u.interpolate(lambda x: np.vstack([E0 * np.ones_like(x[0]), 0.0 * x[0], 0.0 * x[0]]))
    e_u.x.scatter_forward()
    sigma_f = fem.Function(fem.functionspace(msh, ("DG", 0)))
    sigma_f.x.array[:] = 0.5
    rho_field = build_density_field(msh, PHANTOM_RHO_KG_PER_M3)

    # K5's reciprocal synthetic circuit (`test_i36b_k5_r1_and_current_route_on_synthetic_2x2`).
    z0 = 50.0
    rng = np.random.default_rng(656)
    a = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    z_sym = a + a.T + np.eye(4) * (20.0 + 5.0j)
    z_rc = h17._z_rc(193.5178e-12, 64.0e6)
    z_th = h17._thevenin_feed_impedance(z_rc, z0)
    v_s = h17._source_voltage(z_rc, z0)
    loads = h17._parallel_loads(4, [0, 1], z_rc, z_th)
    feeds, order = ["P1", "P2"], ["P1", "P2", "P3", "P4"]
    columns = column_fields(e_u)

    def _stub(_env, _solver, _specs, driven_id, *, held):
        j = feeds.index(driven_id)
        v = np.zeros(4, dtype=np.complex128)
        v[j] = 1.0
        cur = h17._parallel_currents(z_sym, loads, v)
        return {
            "currents": dict(zip(order, cur)),
            "kind": "held" if held else "solve",
            "t_solve": 0.0,
            "t_post": 0.0,
            "p_src": 1.0,
            "e": columns[j].copy(),
            "sigma": sigma_f,
        }

    monkeypatch.setattr(sys.modules[_arm.__module__], "_drive", _stub)
    env = {"comm": comm, "dx": ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)}
    tips = {"a": (0.03, -0.03, 0.0), "b": (-0.03, 0.03, 0.0)}
    t1 = (0.03, 0.03, 0.0)
    arm = _arm(env, None, "K8", [], feeds, order, z_th, v_s, z0, rho_field, tips, t1)
    a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    return arm, msh, cell_tags, e_u, sigma_f, PHANTOM_RHO_KG_PER_M3, E0, tips, t1, a_1g


@complex_only
@_coding
def test_i36c_k8_arm_glue_on_small_mesh(monkeypatch):
    """(K8), §9 item 695: `_arm` with `_drive` stubbed on a small mesh.  (1, +1j) on
    e₁ = e₂ = E_u gives e_sup = (1 + j) E_u, so |e_sup|² = 2 E₀²: (K8-P) the volume
    powers, (K8-S) the 1 g balls, (K8-E) the T1 read, (K8-k) kinds and reciprocity;
    the single-column control (e₂ = 0) misses (K8-S)'s closed form by 0.5."""
    import numpy as np
    import ufl
    from dolfinx import fem

    from fem_em_solver.post.sar import mass_averaged_sar
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_mass_averaged_sar_standard_masses import QUADRATURE_DEGREE
    from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    arm, msh, cell_tags, e_u, sigma_f, rho, e0, tips, t1, a_1g = _k8_arm(
        monkeypatch, lambda e: (e, e)
    )
    comm = msh.comm
    sigma = 0.5

    def _volume(measure):
        one = fem.Constant(msh, np.complex128(1.0))
        return float(
            np.real(comm.allreduce(fem.assemble_scalar(fem.form(one * measure)), op=MPI.SUM))
        )

    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    v_mesh = _volume(dx)
    v_403 = _volume(dx(403))
    p_expected = 0.5 * sigma * 2.0 * e0**2 * v_mesh
    sar_expected = sigma * e0**2 / rho
    r_total = abs(arm["p_sup"]["total"] / p_expected - 1.0)
    r_rod = abs((arm["p_sup"]["rod"] / arm["p_sup"]["total"]) / (v_403 / v_mesh) - 1.0)
    r_sar = {n: abs(arm["sar"][n] / sar_expected - 1.0) for n in tips}
    r_mass = {}
    for n, c in tips.items():
        ball = mass_averaged_sar(
            e_u,
            sigma=sigma_f,
            rho=rho,
            center=c,
            radius=a_1g,
            comm=comm,
            quadrature_degree=QUADRATURE_DEGREE,
        )
        r_mass[n] = abs(arm["mass"][n] / (rho * ball["volume_m3"]) - 1.0)
    e_ref = np.array([(1.0 + 1.0j) * e0, 0.0, 0.0])
    r_e = float(np.linalg.norm(arm["e_t1"] - e_ref) / np.linalg.norm(e_ref))
    recip = _reciprocity_ratio(arm["s2"])

    with monkeypatch.context() as m:
        control, *_rest = _k8_arm(m, lambda e: (e, fem.Function(e.function_space)))
    miss_c = {n: abs(control["sar"][n] / sar_expected - 1.0) for n in tips}
    if comm.rank == 0:
        print(
            f"\n[GEO-36 I36c K8] a_1g {a_1g:.6e} m; V_mesh {v_mesh:.9e}, V403 {v_403:.9e} m^3; "
            f"(K8-P) total {r_total:.3e}, rod share {r_rod:.3e} (band {EXACT_IDENTITY_RTOL:g}, "
            f"ASSERTED); (K8-S) {r_sar}, mass {r_mass} (ASSERTED)"
        )
        print(
            f"[GEO-36 I36c K8] (K8-E) e_t1 {arm['e_t1']} miss {r_e:.3e} (band "
            f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED), valid {arm['t1_valid']}; (K8-k) kinds "
            f"{arm['kinds']}, (R1) {recip:.3e} (band {RECIPROCITY_BAND:.0e}, ASSERTED); Y_in "
            f"{arm['y_in'][0]:.6e}, {arm['y_in'][1]:.6e}"
        )
        print(
            f"[GEO-36 I36c K8] control e2 = 0: (K8-S) misses {miss_c} (0.5 within 1e-10, > band, "
            f"ASSERTED)",
            flush=True,
        )
    assert r_total <= EXACT_IDENTITY_RTOL, r_total
    assert r_rod <= EXACT_IDENTITY_RTOL, r_rod
    assert max(r_sar.values()) <= EXACT_IDENTITY_RTOL, r_sar
    assert max(r_mass.values()) <= EXACT_IDENTITY_RTOL, r_mass
    assert min(arm["mass"].values()) > 0.0, arm["mass"]
    assert r_e <= DISCRETE_IDENTITY_RTOL, r_e
    assert arm["t1_valid"] is True
    assert arm["kinds"] == ["solve", "held"], arm["kinds"]
    assert recip <= RECIPROCITY_BAND, recip
    assert all(v > EXACT_IDENTITY_RTOL for v in miss_c.values()), miss_c
    assert max(abs(v - 0.5) for v in miss_c.values()) <= 1.0e-10, miss_c


#: K9's node (the 193.5996 pF of ``xl`` entry 48's command).
K9_NODE_PF = 193.5996


def _k9_stubs(m, fx, json_dir, *, p_src_r2=1.0):
    """§9 item 702: everything `_run_drive` reads around its own lines, stubbed on
    `_k8_arm`'s mesh; `_arm` runs real.  Arm M's columns are E_u, arm R's 0.5 E_u, told
    apart by the solver token; the powers close TH-14's identity exactly unless
    ``p_src_r2`` breaks arm R's column 2.  Returns ``(lines, z_sym)``."""
    import dataclasses
    import sys
    from types import SimpleNamespace

    import numpy as np
    import ufl
    from dolfinx import fem

    import fem_em_solver.core as core
    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
    from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M

    @dataclasses.dataclass
    class _Problem:
        material_map: dict
        cell_tags: object

    comm = MPI.COMM_WORLD
    msh, cell_tags, e_u, sigma_f, tips, t1 = fx
    e_half = fem.Function(e_u.function_space)
    e_half.x.array[:] = 0.5 * e_u.x.array
    e_half.x.scatter_forward()

    gel = SimpleNamespace(sigma=SALINE_SIGMA, epsilon_r=SALINE_EPSILON_R, mu_r=1.0)
    rod = SimpleNamespace(sigma=SIGMA_WIRE_S_PER_M, epsilon_r=SALINE_EPSILON_R, mu_r=1.0)
    mmap = {3: gel, 403: rod}
    feeds, order = ["P1", "P2"], ["P1", "P2", "P3", "P4"]
    specs = [SimpleNamespace(port_id=pid) for pid in order]
    problem = _Problem(material_map=mmap, cell_tags=cell_tags)
    ctx = {
        "specs": specs,
        "solver": SimpleNamespace(problem=problem),
        "msh": msh,
        "cell_tags": cell_tags,
    }
    fixture = {"degree": 1, "feeds": feeds, "ids": order, "built": {"ctx": ctx}, "t_build": 0.0}
    env = {
        "comm": comm,
        "msh": msh,
        "cell_tags": cell_tags,
        "omega": 2.0 * np.pi * 64.0e6,
        "dx": ufl.Measure("dx", domain=msh, subdomain_data=cell_tags),
        "acct": {"mesh": msh},
    }

    # K5's reciprocal synthetic circuit, as `_k8_arm`, at this node.
    z0 = 50.0
    rng = np.random.default_rng(656)
    a = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    z_sym = a + a.T + np.eye(4) * (20.0 + 5.0j)
    z_rc = h17._z_rc(K9_NODE_PF * 1e-12, 64.0e6)
    z_th = h17._thevenin_feed_impedance(z_rc, z0)
    loads = h17._parallel_loads(4, [0, 1], z_rc, z_th)
    token_r = object()

    def _stub(_env, solver, _specs, driven_id, *, held):
        j = feeds.index(driven_id)
        arm_r = solver is token_r
        v = np.zeros(4, dtype=np.complex128)
        v[j] = 1.0
        cur = h17._parallel_currents(z_sym, loads, v)
        return {
            "currents": dict(zip(order, cur)),
            "kind": "held" if held else "solve",
            "t_solve": 0.0,
            "t_post": 0.0,
            "e": (e_half if arm_r else e_u).copy(),
            "sigma": sigma_f,
            "p_src": p_src_r2 if arm_r and j == 1 else 1.0,
            "p_sheets": 0.5,
            "p_phantom": 0.2,
            "p_rest": 0.2,
            "p_surf": 0.1,
            "p_rod": 0.0,
            "p_sheet_by_id": {"P1": 0.1, "P2": 0.1, "P3": 0.15, "P4": 0.15},
        }

    lines = []
    module = sys.modules[_arm.__module__]
    real_say = module._say

    def _say_recorded(comm_, msg):
        lines.append(msg)
        real_say(comm_, msg)

    m.setattr(module, "_route_drive", lambda _monkeypatch: [(1, mmap)])
    m.setattr(c4, "_build_c4_fixture", lambda _comm: fixture)
    m.setattr(module, "_drive_env", lambda _fixture, _comm: env)
    m.setattr(c4, "_esr_specs", lambda *_a: (specs, z_rc, z_th))
    m.setattr(c4, "_z_failures", lambda *_a: ([], 0.0, 0.0))
    m.setattr(core, "TimeHarmonicSolver", lambda *_a, **_k: token_r)
    m.setattr(module, "_drive", _stub)
    m.setattr(module, "_probe_points", lambda _centre: (tips, t1))
    m.setattr(module, "I36_JSON_DIR", json_dir)
    m.setattr(module, "_say", _say_recorded)
    return lines, z_sym


def _k9_fixture(comm):
    """`_k8_arm`'s mesh, field and probes (copied: `_k8_arm` ends in a call to `_arm`)."""
    import dolfinx
    import numpy as np
    from dolfinx import fem

    msh = dolfinx.mesh.create_box(
        comm, [[-0.05] * 3, [0.05] * 3], [12, 12, 12], dolfinx.mesh.CellType.tetrahedron
    )
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    cells = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
    mids = dolfinx.mesh.compute_midpoints(msh, tdim, cells)
    values = np.where(np.max(np.abs(mids), axis=1) < 0.012, 403, 3).astype(np.int32)
    cell_tags = dolfinx.mesh.meshtags(msh, tdim, cells, values)
    e_u = fem.Function(fem.functionspace(msh, ("N1curl", 1)))
    e_u.interpolate(lambda x: np.vstack([3.0 * np.ones_like(x[0]), 0.0 * x[0], 0.0 * x[0]]))
    e_u.x.scatter_forward()
    sigma_f = fem.Function(fem.functionspace(msh, ("DG", 0)))
    sigma_f.x.array[:] = 0.5
    tips = {"a": (0.03, -0.03, 0.0), "b": (-0.03, 0.03, 0.0)}
    return msh, cell_tags, e_u, sigma_f, tips, (0.03, 0.03, 0.0)


def _k9_json(comm, directory):
    """The node's JSON, parsed on rank 0 and identical on every rank."""
    import json

    rec = None
    if comm.rank == 0:
        with open(os.path.join(directory, f"i36-node-{K9_NODE_PF:.4f}.json")) as fh:
            rec = json.load(fh)
    return comm.bcast(rec, root=0)


@complex_only
@_coding
def test_i36d_k9_run_drive_end_to_end_on_stubs(monkeypatch, tmp_path):
    """(K9), §9 item 702: `_run_drive`'s every line on `_k8_arm`'s mesh, `_arm` real and
    `_drive` stubbed.  Arm R's field is half arm M's on both columns: (K9-Q) every
    tip's Q = 4, (K9-E) η = 2, (K9-A) the census, (K9-J) the JSON, (K9-L) the last
    prints.  Control: arm R's column 2 P_src perturbed by 1e-3 raises naming it, and
    the JSON is still written."""
    import numpy as np

    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    comm = MPI.COMM_WORLD
    fx = _k9_fixture(comm)
    msh, cell_tags = fx[0], fx[1]
    root = comm.bcast(str(tmp_path), root=0)
    good_dir, neg_dir = os.path.join(root, "good"), os.path.join(root, "neg")

    with monkeypatch.context() as m:
        lines, z_sym = _k9_stubs(m, fx, good_dir)
        _run_drive(comm, m, None, K9_NODE_PF)
    rec = _k9_json(comm, good_dir)

    with monkeypatch.context() as m:
        _k9_stubs(m, fx, neg_dir, p_src_r2=1.0 + 1.0e-3)
        with pytest.raises(AssertionError, match=r"R \(gel twin\) \(P\) P2") as neg:
            _run_drive(comm, m, None, K9_NODE_PF)
    rec_neg = _k9_json(comm, neg_dir)

    r_q = {n: abs(v / 4.0 - 1.0) for n, v in rec["q"].items()}
    r_e = abs(rec["eta_t1"] / 2.0 - 1.0)

    def _s(block):
        return np.array([[complex(*v) for v in row] for row in rec[block]["s2"]])

    r_mr = float(np.linalg.norm(_s("M") - _s("R")) / np.linalg.norm(_s("M")))
    z_rc = h17._z_rc(K9_NODE_PF * 1e-12, 64.0e6)
    s_route = h17._circuit_current_route_s(z_sym, [0, 1], z_rc, z_rc, 50.0)
    r_route = float(np.linalg.norm(_s("M") - s_route) / np.linalg.norm(s_route))

    size_local = msh.topology.index_map(3).size_local
    owned = np.asarray(cell_tags.indices) < size_local
    tag_values = np.asarray(cell_tags.values)[owned]
    census = {
        str(t): int(comm.allreduce(int(np.count_nonzero(tag_values == t)), op=MPI.SUM))
        for t in (3, 403)
    }
    a0m_clean = any(
        ln.startswith("(A0m)") and ln.rstrip().endswith("failures none (ASSERTED exact)")
        for ln in lines
    )
    heads = ("(P_RC)", "(A0m)", "tip 1 g ratio", "price:")
    missing_heads = [h for h in heads if not any(ln.startswith(h) for ln in lines)]
    if comm.rank == 0:
        print(
            f"\n[GEO-36 I36d K9] (K9-Q) misses from 4 {r_q} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED); "
            f"(K9-E) miss from 2 {r_e:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED)"
        )
        print(
            f"[GEO-36 I36d K9] (K9-A) census_m {rec['census_m']} vs the mesh's {census}, (A0m) clean "
            f"{a0m_clean} (ASSERTED); (K9-J) keys {sorted(rec)}; s2 M vs R {r_mr:.3e}, vs the circuit "
            f"route {r_route:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED)"
        )
        print(
            f"[GEO-36 I36d K9] (K9-L) missing {missing_heads or 'none'}; control: raised "
            f"{str(neg.value).splitlines()[0]}; JSON still written, R parts[1] p_src "
            f"{rec_neg['R']['parts'][1]['p_src']!r} (ASSERTED)",
            flush=True,
        )
    assert len(r_q) == len(fx[4]) and max(r_q.values()) <= EXACT_IDENTITY_RTOL, r_q
    assert r_e <= DISCRETE_IDENTITY_RTOL, r_e
    assert rec["census_m"] == census, rec["census_m"]
    assert sum(census.values()) == msh.topology.index_map(3).size_global
    assert census["403"] > 0
    assert a0m_clean
    assert {"M", "R"} <= set(rec)
    assert r_mr <= EXACT_IDENTITY_RTOL, r_mr
    assert r_route <= EXACT_IDENTITY_RTOL, r_route
    assert not missing_heads, missing_heads
    assert rec_neg["R"]["parts"][1]["p_src"] == 1.0 + 1.0e-3, rec_neg["R"]["parts"][1]
    assert rec_neg["M"]["parts"][1]["p_src"] == 1.0, rec_neg["M"]["parts"][1]


@complex_only
@_coding
def test_c34_k9_d2_run_drive_on_stubs(monkeypatch, tmp_path):
    """(K9-d2), GEO-38 C34 (§9 item 820): `_run_drive(…, 2, …)` on K9's stubs, with
    `_route_drive` → ``[(2, mmap)]`` and the fixture at degree 2 — C4a's variable
    reads "2" at the fixture call; arm R's solver gets ``degree=2`` (`:553`); the JSON
    lands at ``i36-node-193.5996-d2.json`` with ``"degree": 2``; (K9-Q) 4 and (K9-E) 2
    unchanged.  Control (d2-c): a fixture reporting degree 1 under degree 2 raises at
    the degree assert, before any drive.  No knob is added here."""
    import json
    import sys

    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_geo38_f_clinical_box_loaded import _c34_recorders

    comm = MPI.COMM_WORLD
    module = sys.modules[_run_drive.__module__]
    fx = _k9_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    good_dir, neg_dir = os.path.join(root, "d2"), os.path.join(root, "d2c")

    def _install(m, json_dir, degree):
        _k9_stubs(m, fx, json_dir)
        mmap = c4._build_c4_fixture(comm)["built"]["ctx"]["solver"].problem.material_map
        m.setattr(module, "_route_drive", lambda _monkeypatch: [(2, mmap)])
        return _c34_recorders(m, degree=degree)

    with monkeypatch.context() as m:
        m.delenv(c4.C4_DEGREE_ENV, raising=False)
        env_seen, degrees = _install(m, good_dir, 2)
        _run_drive(comm, m, 2, K9_NODE_PF)
        env_after = os.environ.get(c4.C4_DEGREE_ENV)

    with monkeypatch.context() as m:
        _env_c, degrees_c = _install(m, neg_dir, 1)
        drives = []
        real_drive = module._drive

        def _counted(*a, **k):
            drives.append(a[3])
            return real_drive(*a, **k)

        m.setattr(module, "_drive", _counted)
        with pytest.raises(AssertionError) as d2c:
            _run_drive(comm, m, 2, K9_NODE_PF)

    rec, d1_exists, neg_files = None, None, None
    if comm.rank == 0:
        with open(os.path.join(good_dir, f"i36-node-{K9_NODE_PF:.4f}-d2.json")) as fh:
            rec = json.load(fh)
        d1_exists = os.path.exists(os.path.join(good_dir, f"i36-node-{K9_NODE_PF:.4f}.json"))
        neg_files = sorted(os.listdir(neg_dir)) if os.path.isdir(neg_dir) else []
    rec, d1_exists, neg_files = comm.bcast((rec, d1_exists, neg_files), root=0)
    r_q = {n: abs(v / 4.0 - 1.0) for n, v in rec["q"].items()}
    r_e = abs(rec["eta_t1"] / 2.0 - 1.0)
    if comm.rank == 0:
        print(
            f"\n[GEO-36 C34 K9-d2] C4 variable at the fixture call {env_seen}, after {env_after!r}; "
            f"solver degrees {degrees} (arm R's); JSON -d2 degree {rec['degree']!r}, degree-1 name "
            f"written {d1_exists} (ASSERTED)"
        )
        print(
            f"[GEO-36 C34 K9-d2] (K9-Q) misses from 4 {r_q} (band {EXACT_IDENTITY_RTOL:g}); (K9-E) "
            f"miss from 2 {r_e:.3e} (band {DISCRETE_IDENTITY_RTOL:g}) (ASSERTED)"
        )
        print(
            f"[GEO-36 C34 d2-c] degree 1 under degree 2 raised {str(d2c.value).splitlines()[0]}; "
            f"drives {drives}, solvers {degrees_c}, files {neg_files} (ASSERTED)",
            flush=True,
        )
    assert env_seen == ["2"], env_seen
    assert env_after is None, env_after
    assert degrees == [2], degrees
    assert rec["degree"] == 2, rec["degree"]
    assert not d1_exists
    assert len(r_q) == len(fx[4]) and max(r_q.values()) <= EXACT_IDENTITY_RTOL, r_q
    assert r_e <= DISCRETE_IDENTITY_RTOL, r_e
    assert "(1, 2)" in str(d2c.value), str(d2c.value)
    assert drives == [] and degrees_c == [] and neg_files == [], (drives, degrees_c, neg_files)
