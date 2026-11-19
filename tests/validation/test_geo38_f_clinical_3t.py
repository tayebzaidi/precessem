"""`GEO-38` C8 (§9 item 518, §10 chain C) — the F-clinical-3T twin through the
1.5 T generator: the 1.5 T test bodies run unchanged on a 3 T kwargs table.

Selected by environment only (``FEM_EM_GEO38_C8`` in {``r1``, ``r2``, ``box``}),
never by ``-k``. Brief: docs/planning/chunks/GEO-38.md, "§9 item 518" and
"Review ruling 2026-10-18 03:00".

The geometry knob is the kwargs provider: each windowed test swaps
``_f_clinical_kwargs`` in the 1.5 T modules' own namespaces for the 3 T table
(``monkeypatch``) and calls the 1.5 T body, so every C1 (a) / C2 anchor is
asserted at its own band; no 1.5 T module and no ``src/`` file is edited.

Rung length (review ruling (b), 2026-10-18 03:00): the weekly's 0.58 m is the
rung (leg) span, i.e. ``coil_length`` (the body's (L) ties the cavity wall's
max z to coil_length/2); ``leg_spacing`` keeps the 1.5 T absolute leg overhang,
0.58 − 2 × 0.0321428571 m. The chosen parameter is printed.

Added anchors after each body returns (``IDENTITY_BAND``): shield Ø 0.743 m and
length 0.70 m from the recorded ``diag``; rung Ø 0.704 m and rung span 0.58 m
from the table; the leg-box corner equals its closed form at 1e-12.
Coding tests (K0)–(K3), (SR0): table derivation, pre-gmsh refusals, and the
gmsh stop-rule watchdog firing inside ``generate``.
Scope: geometry and build price only — no solve, no port gate, no tuning.
"""

from __future__ import annotations

import math
import os
import threading
import time

import pytest
from mpi4py import MPI

import tests.validation.test_geo38_f_clinical_box as s2
import tests.validation.test_geo38_f_clinical_build as s1a

C8_ENV = "FEM_EM_GEO38_C8"
C8_SELECTORS = ("r1", "r2", "box")
GMSH_STOP_S = 450.0
STOP_RULE_RC = 86

# Nominal F-clinical-3T geometry (§10 table, `GEO-38.md:31–34`).
RUNG_DIAMETER_3T_M = 0.704
SHIELD_DIAMETER_3T_M = 0.743
SHIELD_LENGTH_3T_M = 0.70
RUNG_LENGTH_3T_M = 0.58
FREQUENCY_3T_HZ = 128e6  # enters no assertion (build-only)

# Literal copy of the 1.5 T table as of `d94320db`; (K0) asserts the live
# table still equals it, so the derived 3 T table cannot move silently.
F_CLINICAL_15T_FROZEN = dict(
    leg_count=32,
    ring_radius=0.355,
    leg_spacing=0.45,
    coil_length=0.45 + 2.0 * 0.0321428571,
    leg_width=0.012,
    ring_minor_radius=0.004,
    ring_gap_length=0.0171428571,
    port_box_size=(0.004, 0.004, 0.0214285714),
    phantom_radius=0.10969,
    phantom_height=0.65,
    air_padding=(0.03 / 0.07) * 0.355,
    shield_radius=0.3725,
    as_hole=True,
    emit_port_sheets=True,
    ring_sheet_orientation="longitudinal",
    resolution=0.015,
)

LEG_OVERHANG_M = 0.0321428571  # the 1.5 T absolute leg overhang past each ring
_COIL_LENGTH_3T = RUNG_LENGTH_3T_M
_LEG_SPACING_3T = RUNG_LENGTH_3T_M - 2.0 * LEG_OVERHANG_M
_Z_NEED_3T = max(
    0.5 * _COIL_LENGTH_3T,
    0.5 * _LEG_SPACING_3T + F_CLINICAL_15T_FROZEN["ring_minor_radius"],
    0.5 * F_CLINICAL_15T_FROZEN["phantom_height"],
)
F_CLINICAL_3T_OVERRIDES = dict(
    ring_radius=0.5 * RUNG_DIAMETER_3T_M,
    leg_spacing=_LEG_SPACING_3T,
    coil_length=_COIL_LENGTH_3T,
    shield_radius=0.5 * SHIELD_DIAMETER_3T_M,
    # The can spans ±z_extent, z_extent = max(...) + air_padding
    # (`io/mesh.py:5636`); computed, never typed.
    air_padding=0.5 * SHIELD_LENGTH_3T_M - _Z_NEED_3T,
)
F_CLINICAL_3T_KWARGS = dict(s1a.F_CLINICAL_15T_KWARGS, **F_CLINICAL_3T_OVERRIDES)

# Closed form of the refusal's leg-box corner (`io/mesh.py:4815–4828`), brief.
LEG_BOX_CORNER_3T_M = 0.363814913104
LEG_BOX_CORNER_ABS_TOL = 1.0e-12

# Printed predictions (brief, "Predictions, written first"), never asserted.
PREDICTED_CELLS_3T = {"r1": (1.33e6, (1.15e6, 1.55e6)), "r2": (1.06e6, (0.90e6, 1.25e6))}
PREDICTED_U1_3T = {"r1": 1.77e6, "r2": 1.41e6}
PREDICTED_U2_3T = {"r1": 9.1e6, "r2": 7.2e6}
PREDICTED_BUILD_S_3T = {"r1": (400.0, (340.0, 480.0)), "r2": (300.0, (250.0, 380.0))}
PREDICTED_BOX_BUILD_S_3T = (400.0, (350.0, 480.0))

# Version-tagged records ((1*) licence), registered after the last window from
# the measured values (default image, 2026-10-18, item 518); the windowed path
# never reads them. r1 `20261018T170415Z_GEO-38-C8-W1-r1.log` (-n 2), r2
# `20261018T170953Z_GEO-38-C8-W2-r2.log` (-n 2), box
# `20261018T171406Z_GEO-38-C8-W3-box-n4.log` (-n 4).
F_CLINICAL_3T_R1_CELL_RECORD = 1_057_895
F_CLINICAL_3T_R2_CELL_RECORD = 808_126
# (U1 = E, U2 = 2E + 2F)
F_CLINICAL_3T_R1_UNKNOWNS_RECORD = (1_425_765, 7_335_988)
F_CLINICAL_3T_R2_UNKNOWNS_RECORD = (1_062_903, 5_518_344)
F_CLINICAL_3T_BOX_CELL_RECORD = 1_061_517

_RECORD: dict = {}


def _selector():
    sel = os.environ.get(C8_ENV, "").strip()
    if sel and sel not in C8_SELECTORS:
        raise ValueError(f"{C8_ENV}={sel!r} not in {C8_SELECTORS}")
    return sel


def _f_clinical_3t_kwargs(rung):
    return dict(F_CLINICAL_3T_KWARGS, conductor_resolution=s1a.CONDUCTOR_RESOLUTION[rung])


def _loadavg():
    try:
        with open("/proc/loadavg") as fh:
            return fh.read().strip()
    except OSError:
        return "n/a"


class _GmshStopRule:
    """Watchdog over gmsh ``generate`` → ``optimize`` (the interval
    `io/mesh.py:6917–6920` reports as ``mesh wall time``), rank 0 only.

    A ``threading.Timer(stop_s)`` starts on entry to ``generate`` and is
    cancelled when ``optimize`` returns (or on exit). The default action prints
    the stop line and aborts the MPI job with ``STOP_RULE_RC``.
    """

    def __init__(self, monkeypatch, stop_s, label, action=None):
        self.mp = monkeypatch
        self.stop_s = float(stop_s)
        self.label = label
        self.action = action or self._abort
        self.timer = None
        self.fired = False
        self.t_test0 = time.perf_counter()
        self.t_gen = None

    def _abort(self):
        comm = MPI.COMM_WORLD
        print(
            f"[GEO-38 C8] GMSH STOP RULE: gmsh > {self.stop_s:.0f} s (3T, {self.label}, -n "
            f"{comm.size}; generate entered at t={self.t_gen - self.t_test0:.1f} s after test "
            f"start; /proc/loadavg {_loadavg()}) — price recorded as a lower bound; the build "
            f"goes to the weekly's xl list",
            flush=True,
        )
        comm.Abort(STOP_RULE_RC)

    def _fire(self):
        self.fired = True
        self.action()

    def _cancel(self):
        if self.timer is not None:
            self.timer.cancel()
            self.timer.join()
            self.timer = None

    def __enter__(self):
        import gmsh

        if MPI.COMM_WORLD.rank != 0:
            return self
        orig_gen = gmsh.model.mesh.generate
        orig_opt = gmsh.model.mesh.optimize

        def generate(*a, **k):
            if self.timer is None:
                self.t_gen = time.perf_counter()
                self.timer = threading.Timer(self.stop_s, self._fire)
                self.timer.daemon = True
                self.timer.start()
            return orig_gen(*a, **k)

        def optimize(*a, **k):
            try:
                return orig_opt(*a, **k)
            finally:
                self._cancel()

        self.mp.setattr(gmsh.model.mesh, "generate", staticmethod(generate))
        self.mp.setattr(gmsh.model.mesh, "optimize", staticmethod(optimize))
        return self

    def __exit__(self, *exc):
        self._cancel()
        return False


def _install_recording_domain(monkeypatch, label):
    from fem_em_solver.io.mesh import MeshGenerator

    orig = MeshGenerator.birdcage_port_domain
    comm = MPI.COMM_WORLD

    def wrapper(*a, **k):
        t0 = time.perf_counter()
        if comm.rank == 0:
            print(
                f"[GEO-38 C8] build start t={time.time():.1f} ({label}, -n {comm.size}; "
                f"/proc/loadavg {_loadavg()})",
                flush=True,
            )
        with _GmshStopRule(monkeypatch, GMSH_STOP_S, label):
            out = orig(*a, **k)
        msh, diag = out[0], out[3]
        _RECORD["diag"] = diag
        _RECORD["kw"] = dict(k)
        _RECORD["cells"] = int(msh.topology.index_map(msh.topology.dim).size_global)
        if comm.rank == 0:
            print(
                f"[GEO-38 C8] build done {time.perf_counter() - t0:.1f} s; gmsh "
                f"{diag.get('mesh_wall_time_s', float('nan')):.1f} s; {_RECORD['cells']} cells; "
                f"/proc/loadavg {_loadavg()}",
                flush=True,
            )
        return out

    monkeypatch.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(wrapper))


def _assert_3t_table(label, box=False):
    diag, kw = _RECORD["diag"], _RECORD["kw"]
    band = s1a.IDENTITY_BAND
    shield_d = 2.0 * float(diag["shield_radius_m"])
    shield_l = float(diag["shield_length_m"])
    rung_d = 2.0 * kw["ring_radius"]
    rung_l = kw["coil_length"]
    corner = float(diag["leg_box_corner_radius_m"])
    z_e = 0.5 * shield_l
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[GEO-38 C8] 3T table ({label}): shield Ø {shield_d:.12f} (vs {SHIELD_DIAMETER_3T_M}), "
            f"length {shield_l:.12f} (vs {SHIELD_LENGTH_3T_M}); rung Ø {rung_d:.12f} (vs "
            f"{RUNG_DIAMETER_3T_M}); rung (leg) span = coil_length {rung_l:.12f} (vs "
            f"{RUNG_LENGTH_3T_M}; ruling (b): the chosen parameter is coil_length, leg_spacing "
            f"{kw['leg_spacing']:.10f} keeps the 1.5 T overhang {LEG_OVERHANG_M}); leg-box corner "
            f"{corner:.12f} vs closed form {LEG_BOX_CORNER_3T_M} (|diff| "
            f"{abs(corner - LEG_BOX_CORNER_3T_M):.3e}); clearance R_s - corner "
            f"{(0.5 * shield_d - corner) * 1e3:.4f} mm; cap clearances: phantom ends "
            f"{(z_e - 0.5 * kw.get('phantom_height', 0.65)) * 1e3:.2f} mm (ruling (c), named for "
            f"the weekly), leg ends {(z_e - 0.5 * rung_l) * 1e3:.2f} mm, ring tubes "
            f"{(z_e - 0.5 * kw['leg_spacing'] - kw['ring_minor_radius']) * 1e3:.2f} mm",
            flush=True,
        )
    assert abs(shield_d / SHIELD_DIAMETER_3T_M - 1.0) <= band, shield_d
    assert abs(shield_l / SHIELD_LENGTH_3T_M - 1.0) <= band, shield_l
    assert abs(rung_d / RUNG_DIAMETER_3T_M - 1.0) <= band, rung_d
    assert abs(rung_l / RUNG_LENGTH_3T_M - 1.0) <= band, rung_l
    assert kw["leg_count"] == s1a.LEG_COUNT == 32
    assert len(s1a._ring_ports()) == 64
    assert abs(corner - LEG_BOX_CORNER_3T_M) <= LEG_BOX_CORNER_ABS_TOL, corner
    if box:
        assert tuple(diag["phantom_box"]["size_m"]) == (0.42, 0.09, 0.65), diag["phantom_box"]


def test_geo38_c8_3t_identities(monkeypatch):
    rung = _selector()
    if rung not in ("r1", "r2"):
        pytest.skip(f"{C8_ENV} not in {{r1, r2}}: C8 identities are selected by environment only")
    monkeypatch.setenv(s1a.STEP1A_ENV, rung)
    monkeypatch.setattr(s1a, "_f_clinical_kwargs", lambda r: _f_clinical_3t_kwargs(r))
    monkeypatch.setattr(s1a, "PREDICTED_CELLS", PREDICTED_CELLS_3T)
    monkeypatch.setattr(s1a, "PREDICTED_BUILD_S", PREDICTED_BUILD_S_3T)
    monkeypatch.setattr(s1a, "PREDICTED_U1", PREDICTED_U1_3T)
    monkeypatch.setattr(s1a, "PREDICTED_U2", PREDICTED_U2_3T)
    _install_recording_domain(monkeypatch, rung)
    s1a.test_geo38_step1a_f_clinical_identities()
    _assert_3t_table(rung)
    if MPI.COMM_WORLD.rank == 0:
        rho = _RECORD["cells"] / s1a.F_HUMAN_N1
        print(
            f"[GEO-38 C8] brief-rho projection (rho = N1/369387 on cells = {rho:.4f}; printed "
            f"only): degree-2 xl t ≈ t_b + {287.1 * rho ** 0.84:.1f} s, memory "
            f"{89.31 * rho ** 0.84:.1f} GiB (x rho^(4/3): {89.31 * rho ** (4 / 3):.1f} GiB)",
            flush=True,
        )


def test_geo38_c8_3t_box_gates(monkeypatch):
    if _selector() != "box":
        pytest.skip(f"{C8_ENV} != box: C8 box gates are selected by environment only")
    monkeypatch.setenv(s2.STEP2_ENV, "1")
    monkeypatch.setattr(s2, "_f_clinical_kwargs", lambda r: _f_clinical_3t_kwargs(r))
    _install_recording_domain(monkeypatch, "box")
    if MPI.COMM_WORLD.rank == 0:
        pb, br = PREDICTED_BOX_BUILD_S_3T
        print(
            f"[GEO-38 C8] box PREDICTED build ≈{pb} s {br}; cells ≈ W1 + 0.16 %; the (A0) line "
            f"below compares with the 1.5 T record and is meaningless at 3 T",
            flush=True,
        )
    s2.test_geo38_step2_astm_box_mesh_gates()
    _assert_3t_table("box", box=True)


# --- coding tests (skip whenever a window selector is set) -------------------

_coding = pytest.mark.skipif(
    bool(os.environ.get(C8_ENV, "").strip()),
    reason=f"{C8_ENV} set: windowed run, coding tests skipped",
)


@_coding
def test_geo38_c8_k0_table_derivation(monkeypatch):
    assert s1a.F_CLINICAL_15T_KWARGS == F_CLINICAL_15T_FROZEN
    for key, val in s1a.F_CLINICAL_15T_KWARGS.items():
        if key not in F_CLINICAL_3T_OVERRIDES:
            assert F_CLINICAL_3T_KWARGS[key] == val, key
    kw = F_CLINICAL_3T_KWARGS
    z_e = (
        max(
            0.5 * kw["coil_length"],
            0.5 * kw["leg_spacing"] + kw["ring_minor_radius"],
            0.5 * kw["phantom_height"],
        )
        + kw["air_padding"]
    )
    R, g, w = kw["ring_radius"], kw["ring_gap_length"], 0.010
    alpha = g / (2.0 * R)
    chord = 2.0 * R * math.tan(alpha)
    print(
        f"\n[GEO-38 C8] (K0) 3T table {F_CLINICAL_3T_OVERRIDES}; z_e {z_e!r} (== 0.35: "
        f"{z_e == 0.35}); alpha {alpha:.9f} rad; chord {chord:.9e} m; sheet area "
        f"{chord * w:.9e} m^2; can volume {math.pi * kw['shield_radius'] ** 2 * 2 * z_e:.6f} m^3",
        flush=True,
    )
    assert abs(z_e / 0.35 - 1.0) <= 1.0e-15, z_e
    assert abs(kw["coil_length"] / RUNG_LENGTH_3T_M - 1.0) <= s1a.IDENTITY_BAND
    assert abs(chord * w / 1.714624621e-4 - 1.0) <= 1.0e-9
    monkeypatch.setenv(C8_ENV, "r3")
    with pytest.raises(ValueError):
        _selector()


@_coding
def test_geo38_c8_k1_3t_table_reaches_builder(monkeypatch):
    with pytest.raises(AssertionError, match="reached the gmsh builder"):
        s1a._refuse_before_gmsh(monkeypatch, **_f_clinical_3t_kwargs("r2"))


@_coding
def test_geo38_c8_k2_attempt1_box_refused(monkeypatch):
    kw = _f_clinical_3t_kwargs("r2")
    kw["port_box_size"] = s1a.ATTEMPT1_PORT_BOX_SIZE
    with pytest.raises(ValueError, match="leg-box corner") as exc:
        s1a._refuse_before_gmsh(monkeypatch, **kw)
    print(f"\n[GEO-38 C8] (K2) refusal: {exc.value}", flush=True)


@_coding
def test_geo38_c8_k3_shield_threshold(monkeypatch):
    kw = _f_clinical_3t_kwargs("r2")
    kw["shield_radius"] = 0.3645
    with pytest.raises(ValueError) as exc:
        s1a._refuse_before_gmsh(monkeypatch, **kw)
    print(f"\n[GEO-38 C8] (K3) R_s 0.3645 refused: {exc.value}", flush=True)
    kw["shield_radius"] = 0.3650
    with pytest.raises(AssertionError, match="reached the gmsh builder"):
        s1a._refuse_before_gmsh(monkeypatch, **kw)


@_coding
def test_geo38_c8_sr0_watchdog_fires_inside_gmsh(monkeypatch):
    import gmsh

    if MPI.COMM_WORLD.rank != 0:
        pytest.skip("SR0 runs on rank 0")
    baseline = threading.active_count()
    seen = {}

    def record():
        seen["fired_in_generate"] = seen.get("in_generate", False)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)

        def build():
            gmsh.clear()
            gmsh.model.occ.addBox(0, 0, 0, 1, 1, 1)
            gmsh.model.occ.synchronize()
            gmsh.option.setNumber("Mesh.MeshSizeMax", 0.03)

        def run(stop_s, label):
            # A fresh MonkeyPatch per rule: a stacked earlier rule's wrapper
            # would otherwise re-arm its own timer (W0 first run).
            with (
                pytest.MonkeyPatch.context() as mp,
                _GmshStopRule(mp, stop_s, label, action=record),
            ):

                def flagged(*a, **k):
                    seen["in_generate"] = True
                    try:
                        return gmsh.model.mesh.generate(*a, **k)
                    finally:
                        seen["in_generate"] = False

                t0 = time.perf_counter()
                flagged(3)
                dt = time.perf_counter() - t0
                gmsh.model.mesh.optimize("Netgen")
            return dt

        build()
        dt = run(0.25, "SR0")
        print(
            f"\n[GEO-38 C8] (SR0) unit cube generate(3) {dt:.2f} s; fired in generate: "
            f"{seen.get('fired_in_generate')}",
            flush=True,
        )
        assert dt >= 0.25, f"inconclusive: generate returned in {dt:.3f} s < 0.25 s"
        assert seen.get("fired_in_generate") is True
        seen.clear()
        build()
        run(600.0, "SR0-long")
        assert "fired_in_generate" not in seen
        assert threading.active_count() == baseline, (threading.active_count(), baseline)
    finally:
        gmsh.finalize()
