"""GEO-38 C12 (§9 item 595, §10 chain C12): the F-clinical-3T twin's ESR
root module — C4a's module (``test_geo38_f_clinical_esr_root.py``) composed
with C8's geometry knob at 128 MHz, the way C11
(``test_geo38_f_clinical_3t_matrix.py``) composed C3a's module.

Composition, never an edit (review ruling 2026-10-22 03:00 (1)): C4a's
``_run_build`` / ``_run_nodes`` run unchanged; ``_c12_patches`` applies C11's
seven patches by calling ``c11._c11_patches`` (no copy), then three more on
module globals C4a reads at call time —

8.  ``c4.FREQUENCY_HZ_C4``      -> ``t3.FREQUENCY_3T_HZ`` (C4a's (f-C4) check)
9.  ``c4._c4_node_json_path``   -> ``_c12_node_json_path`` (1.5 T node JSONs safe)
10. ``c4._apply_degree_c4``     -> a recorder that calls the original and stores
    ``(degree, problem.frequency_hz, problem.material_map)`` — never the solver

— all ten restored at teardown by ``monkeypatch``.

Selected by ``FEM_EM_GEO38_C12`` only (never by ``-k``): unset ⇒ the gated
test skips and the coding tests (K0), (K1), (Z128), (N-f), (N-json) run;
``build`` ⇒ build-only (C4a's fixture anchors against C8's 3 T records through
the patches, (Z) at 128 MHz, (Deg-2), plus (f), (T3), (G) and the file guard
here); ``1`` ⇒ the 3 T node run (the window the 10-28 weekly may commission;
never run by item 595 — it assert-stops before any build on the absent C11
full-branch npz, ruling (4)); anything else ⇒ ``ValueError``.
``FEM_EM_GEO38_C12_NODES_PF`` / ``FEM_EM_GEO38_C12_DEGREE`` are copied into
C4a's knobs under ``1`` (C4a's parsers validate them) and refused under
``build``. A C12 window never reads a C4 or C11 variable: those (and C3, C8,
any STEP*) are refused on entry to the gated test.

The gel at 128 MHz: (G) asserts the problem the solver holds carries the
module constants ``SALINE_SIGMA`` / ``SALINE_EPSILON_R`` — the 2026-10-21
weekly's *provisional* ruling (W3), not a claim about F2182's nominal.

Imports are module objects inside function bodies, never ``test_*`` names.

Scope: the module, its coding tests and a build-only proof. Not claimed: any
Y_in, root, C_in, C_mode or tuning statement at 3 T; any S entry.
"""

from __future__ import annotations

import glob
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

C12_ENV = "FEM_EM_GEO38_C12"
C12_NODES_ENV = "FEM_EM_GEO38_C12_NODES_PF"
C12_DEGREE_ENV = "FEM_EM_GEO38_C12_DEGREE"
FOREIGN_ENVS = (
    "FEM_EM_GEO38_C4",
    "FEM_EM_GEO38_C4_NODES_PF",
    "FEM_EM_GEO38_C4_DEGREE",
    "FEM_EM_GEO38_C11",
    "FEM_EM_GEO38_C11_DRIVES",
    "FEM_EM_GEO38_C3",
    "FEM_EM_GEO38_C8",
)
FOREIGN_PREFIX = "FEM_EM_GEO38_STEP"

# Printed only (brief "Predictions, written first"; never asserted).
PREDICTED_T_BUILD_S = (316.0, (290.0, 560.0))
PREDICTED_WINDOW_S = 318.0
PREDICTED_RSS_SUM_GIB = 2.0
PREDICTED_F = 2_242_229
PREDICTED_Z_RC = complex(0.2, -12.43397993)
C11_PRINTED_3T_COLUMN_S = 318.0  # C11's printed 3 T column (projection input only)
ENTRY40_15T_COLUMN_S = 199.0  # entry 40's measured 1.5 T column (beside it, printed)

_REC_DEG: list = []


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-38 C12] {msg}", flush=True)


def _modules():
    import tests.validation.test_geo38_f_clinical_3t as t3
    import tests.validation.test_geo38_f_clinical_3t_matrix as c11
    import tests.validation.test_geo38_f_clinical_build as s1a
    import tests.validation.test_geo38_f_clinical_column as col
    import tests.validation.test_geo38_f_clinical_esr_root as c4
    import tests.validation.test_geo38_f_clinical_matrix as c3
    import tests.validation.test_wf7_f_human_b1_quadrature as wf7

    return {"c4": c4, "c11": c11, "c3": c3, "col": col, "s1a": s1a, "t3": t3, "wf7": wf7}


# ---------------------------------------------------------------------------
# Knobs.
# ---------------------------------------------------------------------------
def _mode():
    raw = os.environ.get(C12_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("build", "1"):
        return raw
    raise ValueError(f"{C12_ENV} must be unset, 'build' or '1', got {raw!r}")


def _guard_foreign():
    """A C12 window never reads a C4 / C11 / C3 / C8 / STEP* variable."""
    bad = [k for k in FOREIGN_ENVS if k in os.environ]
    bad += sorted(k for k in os.environ if k.startswith(FOREIGN_PREFIX))
    if bad:
        raise ValueError(f"foreign GEO-38 variables set in a C12 window: {bad}")


def _copy_node_knobs(monkeypatch, mode):
    """``build``: either C12 node/degree knob set ⇒ ValueError.  ``1``: copy
    them into C4a's knobs so C4a's own parsers validate them."""
    c4 = _modules()["c4"]
    vals = {k: os.environ.get(k, "").strip() for k in (C12_NODES_ENV, C12_DEGREE_ENV)}
    if mode == "build":
        bad = [k for k, v in vals.items() if v != ""]
        if bad:
            raise ValueError(f"{bad} set under {C12_ENV}=build")
        return
    if mode != "1":
        raise ValueError(mode)
    for src, dst in ((C12_NODES_ENV, c4.C4_NODES_ENV), (C12_DEGREE_ENV, c4.C4_DEGREE_ENV)):
        if vals[src] == "":
            monkeypatch.delenv(dst, raising=False)
        else:
            monkeypatch.setenv(dst, vals[src])


def _c12_node_json_path(c_pf, degree):
    """``c12-node-<pF:.4f>-128mhz.json`` at degree 1 (``degree is None``), a
    ``-128mhz-d2.json`` sibling at degree 2; anything else raises."""
    c4 = _modules()["c4"]
    if degree is None:
        return os.path.join(c4.JSON_DIR, f"c12-node-{c_pf:.4f}-128mhz.json")
    if degree == 2:
        return os.path.join(c4.JSON_DIR, f"c12-node-{c_pf:.4f}-128mhz-d2.json")
    raise ValueError(degree)


def _c12_patches(monkeypatch):
    """C11's seven by call, then patches 8-10 on C4a's call-time globals."""
    m = _modules()
    c4, c11, t3 = m["c4"], m["c11"], m["t3"]
    c11._c11_patches(monkeypatch)
    _REC_DEG.clear()
    monkeypatch.setattr(c4, "FREQUENCY_HZ_C4", t3.FREQUENCY_3T_HZ)
    monkeypatch.setattr(c4, "_c4_node_json_path", _c12_node_json_path)
    orig = c4._apply_degree_c4

    def degree_recorder(ctx, comm):
        out = orig(ctx, comm)
        problem = out[0].problem
        _REC_DEG.append(
            (int(out[1]), float(problem.frequency_hz), dict(problem.material_map or {}))
        )
        return out

    monkeypatch.setattr(c4, "_apply_degree_c4", degree_recorder)


def _record_files():
    """``{path: mtime_ns}`` of every c4-node-* / c12-node-* under JSON_DIR."""
    c4 = _modules()["c4"]
    out = {}
    for pat in ("c4-node-*", "c12-node-*"):
        for p in glob.glob(os.path.join(c4.JSON_DIR, pat)):
            out[p] = os.stat(p).st_mtime_ns
    return out


def _path_state(p):
    return (p.exists(), p.stat().st_mtime_ns if p.exists() else None)


# ---------------------------------------------------------------------------
# Coding tests; skip whenever FEM_EM_GEO38_C12 is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(
    os.environ.get(C12_ENV, "").strip() != "", reason=f"{C12_ENV} set: coding tests skip in windows"
)


@_coding
def test_c12_k0_knob_parsing(monkeypatch):
    c4 = _modules()["c4"]
    for k in (C12_ENV, C12_NODES_ENV, C12_DEGREE_ENV, c4.C4_NODES_ENV, c4.C4_DEGREE_ENV):
        monkeypatch.delenv(k, raising=False)
    assert _mode() is None
    for good in ("build", "1"):
        monkeypatch.setenv(C12_ENV, good)
        assert _mode() == good
    for bad in ("2", "x"):
        monkeypatch.setenv(C12_ENV, bad)
        with pytest.raises(ValueError):
            _mode()
    monkeypatch.delenv(C12_ENV, raising=False)

    # build: either node knob set ⇒ ValueError; neither ⇒ accepted.
    _copy_node_knobs(monkeypatch, "build")
    for k, v in ((C12_NODES_ENV, "60.25,58.5"), (C12_DEGREE_ENV, "2")):
        monkeypatch.setenv(k, v)
        with pytest.raises(ValueError):
            _copy_node_knobs(monkeypatch, "build")
        monkeypatch.delenv(k)

    # 1: copied into C4a's knobs; C4a's parsers then validate.
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv(C12_NODES_ENV, "60.25,58.5")
        mp.setenv(C12_DEGREE_ENV, "2")
        _copy_node_knobs(mp, "1")
        assert os.environ[c4.C4_NODES_ENV] == "60.25,58.5"
        assert c4._nodes_pf() == [58.5, 60.25]
        assert c4._degree_knob_c4() == 2
        mp.setenv(C12_NODES_ENV, "58.50001,58.50002")
        mp.delenv(C12_DEGREE_ENV)
        _copy_node_knobs(mp, "1")
        assert c4.C4_DEGREE_ENV not in os.environ
        assert c4._degree_knob_c4() is None
        with pytest.raises(ValueError):
            c4._nodes_pf()
    assert c4.C4_NODES_ENV not in os.environ and c4.C4_DEGREE_ENV not in os.environ

    # The foreign-variable guard.
    _guard_foreign()
    for k in FOREIGN_ENVS + ("FEM_EM_GEO38_STEP1B", "FEM_EM_GEO38_STEP2"):
        with pytest.MonkeyPatch.context() as mp:
            mp.setenv(k, "1")
            with pytest.raises(ValueError):
                _guard_foreign()
    _guard_foreign()
    print(
        f"\n[GEO-38 C12 K0] C12 unset/build/1 accepted, 2/x refused; node/degree knobs refused "
        f"under build, copied under 1 (C4a parsers: [58.5, 60.25], collision refused); "
        f"{len(FOREIGN_ENVS) + 2} foreign variables refused",
        flush=True,
    )


class _Sentinel(Exception):
    pass


@_coding
def test_c12_k1_routing_tripwire():
    """(K1): under the ten patches C4a's fixture reaches gmsh with the 3 T
    table, 128 MHz, C11's npz and the C12 JSON name; on exit all are restored."""
    from fem_em_solver.io.mesh import MeshGenerator

    m = _modules()
    c4, c11, c3, s1a, wf7 = m["c4"], m["c11"], m["c3"], m["s1a"], m["wf7"]
    orig_save, orig_json, orig_deg = c3.SAVE_PATH, c4._c4_node_json_path, c4._apply_degree_c4
    orig_domain = MeshGenerator.birdcage_port_domain
    seen = {}

    def recorder(*a, **k):
        seen.update(k)
        raise _Sentinel

    with pytest.MonkeyPatch.context() as mp:
        _c12_patches(mp)
        assert c4.FREQUENCY_HZ_C4 == 128e6
        assert wf7.FREQUENCY_HZ == 128e6
        assert c3.SAVE_PATH == c11.C11_SAVE_PATH
        name = os.path.basename(c4._c4_node_json_path(100.0, None))
        assert name == "c12-node-100.0000-128mhz.json", name
        mp.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(recorder))
        with pytest.raises(_Sentinel):
            c4._build_c4_fixture(MPI.COMM_WORLD)
    print(
        f"\n[GEO-38 C12 K1] tripwire: ring_radius {seen.get('ring_radius')}, shield_radius "
        f"{seen.get('shield_radius')}; JSON name {name}",
        flush=True,
    )
    assert abs(seen["ring_radius"] - 0.352) <= 1e-15, seen["ring_radius"]
    assert abs(seen["shield_radius"] - 0.3715) <= 1e-15, seen["shield_radius"]
    assert c4.FREQUENCY_HZ_C4 == 64e6
    assert wf7.FREQUENCY_HZ == 64e6
    assert s1a._f_clinical_kwargs("r1")["ring_radius"] == 0.355
    assert c3.SAVE_PATH == orig_save
    assert c4._c4_node_json_path is orig_json
    assert c4._apply_degree_c4 is orig_deg
    assert MeshGenerator.birdcage_port_domain is orig_domain


@_coding
def test_c12_z128_specs_and_n_z():
    """(Z128): C4a's (Z) on synthetic 64-port specs at 128 MHz; (N-Z): the same
    specs checked against 64 MHz fail with worst Im rel exactly 0.5."""
    from fem_em_solver.ports.lumped import LumpedSheetPortSpec

    m = _modules()
    c4, t3 = m["c4"], m["t3"]
    base = [
        LumpedSheetPortSpec(
            port_id=f"P{i}",
            facet_tag=1000 + i,
            port_impedance_ohm=50.0,
            gap_height_m=0.01,
            sheet_width_m=0.01,
            drive_direction=(0.0, 0.0, 1.0),
        )
        for i in range(33, 97)
    ]
    ids = [sp.port_id for sp in base]
    feeds = ["P33", "P41"]
    c_f = c4.SPEC_WIRING_C_PF * 1e-12
    f3 = t3.FREQUENCY_3T_HZ
    specs, z_rc, z_th = c4._esr_specs(base, feeds, c_f, f3, 50.0)
    fails, im_w, th_w = c4._z_failures(specs, ids, feeds, c_f, f3, 50.0)
    x_ref = -1.0 / (2.0 * np.pi * 128e6 * 100e-12)
    fails64, im64, th64 = c4._z_failures(specs, ids, feeds, c_f, 64e6, 50.0)
    print(
        f"\n[GEO-38 C12 Z128] Z_RC {z_rc:.9e} (Im ref {x_ref:.9e}), Z_th {z_th:.9e}, fails "
        f"{fails or 'none'}; N-Z vs 64 MHz: worst Im rel {im64:.15f}, Z_th rel {th64:.3e}, "
        f"{len(fails64)} failures",
        flush=True,
    )
    assert f3 == 128e6
    assert not fails, fails
    assert z_rc.real == c4.ESR_OHM == 0.2
    assert abs(z_rc.imag / x_ref - 1.0) <= 1e-12, (z_rc, x_ref)
    assert abs(z_rc.imag - PREDICTED_Z_RC.imag) <= 1e-8, z_rc
    assert fails64
    assert abs(im64 - 0.5) <= 1e-12, im64
    assert th64 > 1e-15
    assert any("Z_th" in s for s in fails64), fails64


@_coding
def test_c12_n_f_frequency_patch_load_bearing():
    """(N-f), asserted by construction: C11's patches alone leave C4a's global
    at 64 MHz and its (f-C4) check at ω = 2π·128e6 reads exactly 1.0."""
    m = _modules()
    c4, t3 = m["c4"], m["t3"]
    omega = 2.0 * np.pi * t3.FREQUENCY_3T_HZ
    with pytest.MonkeyPatch.context() as mp:
        m["c11"]._c11_patches(mp)
        assert c4.FREQUENCY_HZ_C4 == 64e6
        f_rel_c11 = omega / (2.0 * np.pi * c4.FREQUENCY_HZ_C4) - 1.0
    with pytest.MonkeyPatch.context() as mp:
        _c12_patches(mp)
        f_rel_c12 = omega / (2.0 * np.pi * c4.FREQUENCY_HZ_C4) - 1.0
    print(
        f"\n[GEO-38 C12 N-f] C11 patches only: f-C4 rel {f_rel_c11!r} (ASSERTED 1.0 > 1e-12); "
        f"C12 patches: {f_rel_c12!r}",
        flush=True,
    )
    assert abs(f_rel_c11 - 1.0) <= 1e-15, f_rel_c11
    assert abs(f_rel_c11) > 1e-12
    assert abs(f_rel_c12) <= 1e-12, f_rel_c12


@_coding
def test_c12_n_json_name_patch_load_bearing():
    """(N-json): without patch 9 a 3 T node at 100 pF writes the 1.5 T name;
    with it the name differs from every c4-node-* name."""
    c4 = _modules()["c4"]
    unpatched = os.path.basename(c4._c4_node_json_path(100.0, None))
    assert unpatched == "c4-node-100.0000.json", unpatched
    c4_names = {c4._c4_node_json_path(c, d) for c in (58.5, 100.0, 120.1234) for d in (None, 2)}
    c12_names = {_c12_node_json_path(c, d) for c in (58.5, 100.0, 120.1234) for d in (None, 2)}
    with pytest.raises(ValueError):
        _c12_node_json_path(58.5, 1)
    print(f"\n[GEO-38 C12 N-json] unpatched {unpatched}; C12 names {sorted(c12_names)}", flush=True)
    assert len(c12_names) == 6 and not (c12_names & c4_names)
    for p in c12_names:
        assert os.path.dirname(p) == c4.JSON_DIR
        assert os.path.basename(p).startswith("c12-node-") and "-128mhz" in p


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
@complex_only
def test_geo38_c12_3t_esr_root(monkeypatch):
    _guard_foreign()
    mode = _mode()
    if mode is None:
        pytest.skip(f"{C12_ENV} unset: GEO-38 C12 is selected by environment only")
    _copy_node_knobs(monkeypatch, mode)
    m = _modules()
    c4, c11, c3, t3 = m["c4"], m["c11"], m["c3"], m["t3"]
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

    comm = MPI.COMM_WORLD
    c3_save = c3.SAVE_PATH
    c3_state = _path_state(c3_save)
    c11_state = _path_state(c11.C11_SAVE_PATH)
    files_before = _record_files()
    tb, tb_range = PREDICTED_T_BUILD_S
    _say(
        comm,
        f"PREDICTED (printed, never asserted): (A0) rel 0; (Az) worst ≈5.7e-14 deg; t_build "
        f"≈{tb} s {tb_range} at loadavg <= 4; window ≈{PREDICTED_WINDOW_S} s; summed RSS "
        f"≈{PREDICTED_RSS_SUM_GIB} GiB after degree 2; F = {PREDICTED_F}; (F) P33 / P41 at "
        f"5.625 / 95.625 deg, bisector 50.625; (Z) Z_RC {PREDICTED_Z_RC}; C4a's (f-C4) label "
        f"reads 'rel to 64 MHz' (string literal) — C12's (f) line is the reading of record",
    )
    _say(
        comm,
        f"start: mode {mode}, -n {comm.size}; C3 npz {c3_save} (exists, mtime) {c3_state}; "
        f"C11 npz {c11.C11_SAVE_PATH} exists {c11_state[0]}; record JSONs before "
        f"{len(files_before)}",
    )

    _c12_patches(monkeypatch)
    t0 = time.perf_counter()
    if mode == "build":
        c4._run_build(comm, monkeypatch)
    else:
        c4._run_nodes(comm)
    t_run = time.perf_counter() - t0

    # (f) the frequency reached the fixture.
    omega = c11._REC["omega"]
    omega_ref = 2.0 * np.pi * t3.FREQUENCY_3T_HZ
    f_rel = abs(omega / omega_ref - 1.0)
    _say(
        comm,
        f"(f) recorded omega {omega.real:.12e} vs 2*pi*128e6 {omega_ref:.12e}: rel {f_rel:.3e} "
        f"(band {c11.OMEGA_REL_BAND:.0e}, ASSERTED; the reading of record for f)",
    )
    assert f_rel <= c11.OMEGA_REL_BAND, (omega, omega_ref)

    # (T3) the 3 T table reached the builder.
    s1a = m["s1a"]
    diag, kw = c11._REC["diag"], c11._REC["kw"]
    band = s1a.IDENTITY_BAND
    shield_d = 2.0 * float(diag["shield_radius_m"])
    shield_l = float(diag["shield_length_m"])
    rung_d = 2.0 * kw["ring_radius"]
    _say(
        comm,
        f"(T3) shield Ø {shield_d:.12f} (vs {t3.SHIELD_DIAMETER_3T_M}), length {shield_l:.12f} "
        f"(vs {t3.SHIELD_LENGTH_3T_M}), rung Ø {rung_d:.12f} (vs {t3.RUNG_DIAMETER_3T_M}) "
        f"(IDENTITY_BAND {band:.0e}, ASSERTED)",
    )
    assert abs(shield_d / t3.SHIELD_DIAMETER_3T_M - 1.0) <= band, shield_d
    assert abs(shield_l / t3.SHIELD_LENGTH_3T_M - 1.0) <= band, shield_l
    assert abs(rung_d / t3.RUNG_DIAMETER_3T_M - 1.0) <= band, rung_d

    # (G) the problem the 3 T root would solve (exact equality of imported constants).
    degrees = [d for d, _f, _mm in _REC_DEG]
    for d, f_hz, mm in _REC_DEG:
        gel = mm.get(PHANTOM_CELL_TAG)
        _say(
            comm,
            f"(G) degree {d}: problem.frequency_hz {f_hz!r} (== {t3.FREQUENCY_3T_HZ!r}); "
            f"material_map tags {sorted(mm)}; tag {PHANTOM_CELL_TAG} sigma "
            f"{getattr(gel, 'sigma', None)} eps_r {getattr(gel, 'epsilon_r', None)} mu_r "
            f"{getattr(gel, 'mu_r', None)} (== SALINE_SIGMA {SALINE_SIGMA}, SALINE_EPSILON_R "
            f"{SALINE_EPSILON_R}, 1; provisional gel ruling W3); tag "
            f"{BIRDCAGE_IMPLANT_CELL_TAG} absent {BIRDCAGE_IMPLANT_CELL_TAG not in mm} (ASSERTED)",
        )
    expected_degrees = [1, 2] if mode == "build" else [int(c4._degree_knob_c4() or 1)]
    assert degrees == expected_degrees, degrees
    for d, f_hz, mm in _REC_DEG:
        assert f_hz == t3.FREQUENCY_3T_HZ, (d, f_hz)
        gel = mm[PHANTOM_CELL_TAG]
        assert (gel.sigma, gel.epsilon_r, gel.mu_r) == (SALINE_SIGMA, SALINE_EPSILON_R, 1.0), gel
        assert BIRDCAGE_IMPLANT_CELL_TAG not in mm, sorted(mm)

    if mode != "build":
        return

    # (File guard) build-only wrote and touched no record file.
    c3_after = _path_state(c3_save)
    files_after = _record_files()
    _say(
        comm,
        f"(File guard) C11 npz {c11.C11_SAVE_PATH} exists {c11.C11_SAVE_PATH.exists()} (before "
        f"{c11_state[0]}; ASSERTED False); C3 npz (exists, mtime) {c3_after} == before "
        f"{c3_state} (ASSERTED); c4-node-* / c12-node-* files {len(files_after)} == before "
        f"{len(files_before)}, mtimes unchanged (ASSERTED)",
    )
    assert not c11.C11_SAVE_PATH.exists(), c11.C11_SAVE_PATH
    assert c3_after == c3_state, (c3_after, c3_state)
    assert files_after == files_before, sorted(set(files_after.items()) ^ set(files_before.items()))

    _say(
        comm,
        f"BUILD-ONLY return: c4._run_build {t_run:.1f} s at -n {comm.size}. PRICE INPUT "
        f"(projection, printed only): build {t_run:.1f} s + C11's printed 3 T column "
        f"{C11_PRINTED_3T_COLUMN_S} s = {t_run + C11_PRINTED_3T_COLUMN_S:.1f} s per degree-1 "
        f"node window floor, beside entry 40's measured 1.5 T column {ENTRY40_15T_COLUMN_S} s",
    )
