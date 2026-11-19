"""`GEO-38` C30 (§9 item 790, §10 chain C): I31's kept aspect-5 disc inside
F-clinical-3T's ASTM box — CAD and mesh identities only, no solve.

A twin module, so I31's (`test_imp6_disc_f_clinical.py`, `IMP-6` step 4) stays
byte-identical: its build body runs unchanged under C20's 3 T route,
``m.setattr(s2, "_f_clinical_kwargs", t3._f_clinical_3t_kwargs)``
(`test_geo38_f_clinical_3t_rod_box.py:106–110`). I31's ``_kwargs()`` calls
``_box_kwargs`` by name, and that function reads ``_f_clinical_kwargs`` from
s2's namespace at call time, so the swap reaches it; patching ``s1a`` would not.
The disc (a = 6 mm, c = 1.2 mm, centre (0.015, 0, 0), axis x̂) is I31's — a
modelling choice, no ASTM-site claim. Brief: docs/planning/chunks/GEO-38.md,
"§9 item 790".

The build test is selected by environment only (``FEM_EM_IMP6_C30=1``), never
by ``-k``; it skips when unset. The coding tests are pure and always run.

Anchors: I31's (D0)–(D10) at the 3 T can, asserted inside I31's body (its
(D4) reads ``air_domain_cad_mass_m3``, C20's route; its (D9) is census 64; its
gmsh stop rule is armed at 450 s). Added here, on the captured build:
  (T3) shield Ø 0.743, length 0.70, rung Ø 0.704 at IDENTITY_BAND — the 3 T
       sentinel; the 1.5 T rung Ø 0.710 is a red (coding control below);
  (Bore) the box strictly inside the bore, by C20's ``_bore_gap``;
  U₁ from a function-space build = E; the cell count and (U₁, U₂) registered
  as version-tagged records, binding once registered.
Negative control (asserted, coding): without the 3 T swap, I31's ``_kwargs()``
reaches the builder with the 1.5 T rung Ø, which (T3) rejects.
Scope: build-only. Not claimed: a solve, a tuned coil, a price, a 3 T gel
value, an ASTM placement.
"""

from __future__ import annotations

import os
import resource

import pytest
from mpi4py import MPI

import tests.validation.test_geo38_f_clinical_3t as t3
import tests.validation.test_geo38_f_clinical_box as s2
import tests.validation.test_imp6_disc_f_clinical as i31
from tests.validation.gmsh_stop_rule import GMSH_STOP_ENV
from tests.validation.test_geo38_f_clinical_3t_rod_box import _bore_gap
from tests.validation.test_geo38_f_clinical_build import IDENTITY_BAND, LEG_COUNT

C30_ENV = "FEM_EM_IMP6_C30"
#: Variables a C30 window must never read (I31's own knob is set by the module).
FOREIGN_ENVS = (
    i31.STEP4_ENV,
    "FEM_EM_GEO38_C20",
    "FEM_EM_GEO38_C21",
    "FEM_EM_GEO38_C8",
    "FEM_EM_GEO38_C11",
    "FEM_EM_GEO36_ROD_BOX",
)

#: I31's 1.5 T counts, printed beside the 3 T ones only
#: (`20261018T201127Z_IMP-6-step4-a2-W1-n4.log.gz:11910`, `:11924`).
I31_CELLS_15T = 1_198_350
I31_U2_15T = 8_208_970

#: Disc and box on the 3 T table, version-tagged ((1*) licence), default image,
#: 2026-11-03: `20261103T180215Z_GEO-38-C30-W1-build.log:11946`, `:11949` (-n 2).
F_CLINICAL_3T_DISC_CELL_RECORD = 1_110_290
# (U1 = E, U2 = 2E + 2F): E 1 487 097, F 2 347 021.
F_CLINICAL_3T_DISC_UNKNOWNS_RECORD = (1_487_097, 7_668_236)


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(f"[GEO-38 C30] {msg}", flush=True)


def _guard_foreign():
    bad = [k for k in FOREIGN_ENVS if os.environ.get(k) is not None]
    if bad:
        raise ValueError(f"C30 never reads these variables; unset them: {bad}")


def _mode():
    raw = os.environ.get(C30_ENV, "").strip()
    if raw == "":
        return None
    if raw == "1":
        return raw
    raise ValueError(f"{C30_ENV} must be unset or '1', got {raw!r}")


def _route_3t(monkeypatch):
    """C20's route: the 3 T table in s2's namespace, where `_box_kwargs` reads it."""
    monkeypatch.setattr(s2, "_f_clinical_kwargs", t3._f_clinical_3t_kwargs)


def _capture_builder(monkeypatch):
    """Wrap `birdcage_port_domain` so the build I31's body runs is kept here."""
    from fem_em_solver.io.mesh import MeshGenerator

    orig = MeshGenerator.birdcage_port_domain
    seen = {}

    def _capture(*a, **k):
        seen["kw"] = dict(k)
        out = orig(*a, **k)
        seen["out"] = out
        return out

    monkeypatch.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(_capture))
    return seen


def _t3_devs(rung_d, shield_d=None, shield_l=None):
    devs = {"rung": abs(rung_d / t3.RUNG_DIAMETER_3T_M - 1.0)}
    if shield_d is not None:
        devs["shield_d"] = abs(shield_d / t3.SHIELD_DIAMETER_3T_M - 1.0)
        devs["shield_l"] = abs(shield_l / t3.SHIELD_LENGTH_3T_M - 1.0)
    return devs


# --- Coding tests: pure, no gmsh. ---------------------------------------------


def test_c30_k0_knob_and_foreign_guard(monkeypatch):
    for k in FOREIGN_ENVS + (C30_ENV,):
        monkeypatch.delenv(k, raising=False)
    assert _mode() is None
    _guard_foreign()
    monkeypatch.setenv(C30_ENV, "1")
    assert _mode() == "1"
    for bad in ("2", "3t", "build"):
        monkeypatch.setenv(C30_ENV, bad)
        with pytest.raises(ValueError):
            _mode()
    for k in FOREIGN_ENVS:
        with monkeypatch.context() as m:
            m.setenv(k, "1")
            with pytest.raises(ValueError):
                _guard_foreign()


def test_c30_k1_disc_kwargs_on_the_3t_table(monkeypatch):
    kw15 = i31._kwargs()
    with monkeypatch.context() as m:
        _route_3t(m)
        kw = i31._kwargs()
    for key, val in t3.F_CLINICAL_3T_OVERRIDES.items():
        assert kw[key] == val, key
    assert kw["phantom_box"] == kw15["phantom_box"]
    assert kw["implant_spheroid"] == kw15["implant_spheroid"]
    assert kw["implant_spheroid"]["keep_cells"] is True
    assert "implant_rod" not in kw
    # The sentinel and its control: 3 T rung Ø 0.704 passes, 1.5 T's misses.
    dev3 = _t3_devs(2.0 * kw["ring_radius"])["rung"]
    dev15 = _t3_devs(2.0 * kw15["ring_radius"])["rung"]
    _say(
        f"(T3-k) rung Ø 3 T {2.0 * kw['ring_radius']:.12f} dev {dev3:.3e}; "
        f"1.5 T {2.0 * kw15['ring_radius']:.12f} dev {dev15:.3e} (IDENTITY_BAND {IDENTITY_BAND:.0e})"
    )
    assert dev3 <= IDENTITY_BAND, dev3
    assert dev15 > IDENTITY_BAND, dev15
    # The swap is scoped: I31's `_kwargs()` reads the 1.5 T table again afterwards.
    assert i31._kwargs()["ring_radius"] == kw15["ring_radius"]


def test_c30_k2_bore(monkeypatch):
    with monkeypatch.context() as m:
        _route_3t(m)
        gap, corner = _bore_gap(i31._kwargs())
    _say(f"(Bore-k) largest box-corner radius {corner:.9f} m, gap {gap:.9f} m (ASSERTED > 0)")
    assert gap > 0.0, (gap, corner)


def test_c30_k3_3t_kwargs_reach_the_builder(monkeypatch):
    """The 3 T disc kwargs pass the refusals and reach the builder (tripwire)."""
    from fem_em_solver.io.mesh import MeshGenerator

    _route_3t(monkeypatch)
    seen = _capture_builder(monkeypatch)

    def _tripwire(*a, **k):
        raise AssertionError("reached the gmsh builder")

    monkeypatch.setattr(MeshGenerator, "_build_birdcage_port_model", staticmethod(_tripwire))
    if MPI.COMM_WORLD.rank == 0:
        expected = dict(expected_exception=AssertionError, match="reached the gmsh builder")
    else:
        expected = dict(expected_exception=RuntimeError, match="generation failed on rank 0")
    with pytest.raises(**expected):
        MeshGenerator.birdcage_port_domain(
            comm=MPI.COMM_WORLD, return_diagnostics=True, **i31._kwargs()
        )
    assert _t3_devs(2.0 * seen["kw"]["ring_radius"])["rung"] <= IDENTITY_BAND
    assert seen["kw"]["implant_spheroid"]["keep_cells"] is True


# --- The build (environment-selected). ----------------------------------------


def test_geo38_c30_disc_in_3t_box(monkeypatch):
    if _mode() is None:
        pytest.skip(f"{C30_ENV} unset: GEO-38 C30 is selected by environment only")
    _guard_foreign()
    from dolfinx import fem

    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND

    comm = MPI.COMM_WORLD
    _route_3t(monkeypatch)
    monkeypatch.setenv(i31.STEP4_ENV, "1")
    monkeypatch.setenv(GMSH_STOP_ENV, "1")
    seen = _capture_builder(monkeypatch)
    gap, corner = _bore_gap(i31._kwargs())
    _say(
        f"-n {comm.size}; 3 T route on I31's body; (Bore) corner {corner:.9f} m vs rung radius "
        f"{t3.RUNG_DIAMETER_3T_M / 2} m, gap {gap:.9f} m; /proc/loadavg {t3._loadavg()}"
    )

    # I31's (D0)–(D10), unchanged, at the 3 T can.
    i31.test_imp6_step4_disc_in_astm_box(monkeypatch)
    _say(f"I31 body returned; /proc/loadavg {t3._loadavg()}")

    kw = seen["kw"]
    msh, _cells, _facets, diag = seen["out"]
    n_cells = int(msh.topology.index_map(msh.topology.dim).size_global)
    shield_d = 2.0 * float(diag["shield_radius_m"])
    shield_l = float(diag["shield_length_m"])
    rung_d = 2.0 * float(kw["ring_radius"])
    devs = _t3_devs(rung_d, shield_d, shield_l)
    gap_b, corner_b = _bore_gap(kw)

    msh.topology.create_entities(1)
    msh.topology.create_entities(2)
    E = comm.allreduce(msh.topology.index_map(1).size_local, op=MPI.SUM)
    F = comm.allreduce(msh.topology.index_map(2).size_local, op=MPI.SUM)
    V1 = fem.functionspace(msh, ("N1curl", 1))
    U1 = int(V1.dofmap.index_map.size_global * V1.dofmap.index_map_bs)
    U2 = 2 * E + 2 * F
    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2, op=MPI.SUM)
    miss_15t = abs(n_cells / I31_CELLS_15T - 1.0)
    _say(
        f"{n_cells} cells; beside F_CLINICAL_3T_BOX_CELL_RECORD "
        f"{t3.F_CLINICAL_3T_BOX_CELL_RECORD}: disc adds {n_cells - t3.F_CLINICAL_3T_BOX_CELL_RECORD}; "
        f"I31 1.5 T {I31_CELLS_15T} (|ratio - 1| {miss_15t:.3e}, {miss_15t / CELL_COUNT_BAND:.2f}x "
        f"CELL_COUNT_BAND, printed); gmsh {diag.get('mesh_wall_time_s', float('nan')):.2f} s; "
        f"summed peak RSS {rss:.3f} GiB\n"
        f"[GEO-38 C30] (T3) shield Ø {shield_d:.12f} (vs {t3.SHIELD_DIAMETER_3T_M}), length "
        f"{shield_l:.12f} (vs {t3.SHIELD_LENGTH_3T_M}), rung Ø {rung_d:.12f} (vs "
        f"{t3.RUNG_DIAMETER_3T_M}); devs {devs} (IDENTITY_BAND {IDENTITY_BAND:.0e})\n"
        f"[GEO-38 C30] (Bore) built kwargs corner {corner_b:.9f} m, gap {gap_b:.9f} m\n"
        f"[GEO-38 C30] unknowns E {E}, F {F}; U1 {U1} (function space), U2 = 2E + 2F = {U2} "
        f"(I31 1.5 T U2 {I31_U2_15T}); records {F_CLINICAL_3T_DISC_CELL_RECORD}, "
        f"{F_CLINICAL_3T_DISC_UNKNOWNS_RECORD}"
    )

    # (T3): the 3 T sentinel.
    for key, dev in devs.items():
        assert dev <= IDENTITY_BAND, (key, dev)
    # (Bore)
    assert gap_b > 0.0, (gap_b, corner_b)
    assert n_cells < i31.F_CLINICAL_DISC_CELL_GUARD, n_cells
    assert len(diag["ring_sheet_gmsh_triangles"]) == 2 * LEG_COUNT == 64
    assert U1 == E, (U1, E)
    # The records bind once registered.
    if F_CLINICAL_3T_DISC_CELL_RECORD is not None:
        assert abs(n_cells / F_CLINICAL_3T_DISC_CELL_RECORD - 1.0) <= CELL_COUNT_BAND, n_cells
        assert (U1, U2) == F_CLINICAL_3T_DISC_UNKNOWNS_RECORD, (U1, U2)
