"""GEO-38 C11 (§9 item 566, §10 chain C11): the F-clinical-3T twin's 64×64
copper sweep module — C3a's module (``test_geo38_f_clinical_matrix.py``)
composed with C8's geometry knob (``test_geo38_f_clinical_3t.py``) at 128 MHz.

Composition, never an edit (review ruling 2026-10-20 03:00 (2)): C3a's
``_run_c3`` runs unchanged; ``_c11_patches`` rebinds the module globals it and
the column module read at call time —

1. ``s1a._f_clinical_kwargs``              -> ``t3._f_clinical_3t_kwargs``
2. ``s1a.F_CLINICAL_R1_CELL_RECORD``       -> ``t3.F_CLINICAL_3T_R1_CELL_RECORD``
3. ``s1a.F_CLINICAL_R1_UNKNOWNS_RECORD``   -> ``t3.F_CLINICAL_3T_R1_UNKNOWNS_RECORD``
4. ``wf7.FREQUENCY_HZ``                    -> ``t3.FREQUENCY_3T_HZ`` (128 MHz)
5. ``c3.SAVE_PATH``                        -> ``C11_SAVE_PATH`` (entry 44's npz is safe)
6. a recording wrapper on ``col._context_from_mesh`` (stores ``ctx["omega"]``)
7. a recording wrapper on ``MeshGenerator.birdcage_port_domain`` (stores
   ``kw`` and ``diag``; arms nothing — the 3 T r1 build is priced, so C8's
   450 s gmsh stop rule is not reused)
8. ``c3._degree_knob_c3``                  -> ``_degree_knob_c11`` (§9 item 869, C5e)
9. ``c3.SAVE_PATH_D2``                     -> ``C11_SAVE_PATH_D2`` (entry 54's npz is safe)

— all restored at teardown by ``monkeypatch``.

Selected by ``FEM_EM_GEO38_C11`` only (never by ``-k``): unset ⇒ the gated test
skips and the coding tests (K0), (K1), (N1)–(N3) run; ``build`` ⇒ build-only
(C3a's (A0), (iv), outer census, (S), (C), (W), (Az), (U1) against C8's 3 T
records, plus (f) and (T3) here); ``1`` ⇒ the full 64-drive sweep at 128 MHz
(a future window the weekly may commission; never run by item 566); anything
else ⇒ ``ValueError``. ``FEM_EM_GEO38_C11_DRIVES`` (1..64) is parsed as C3a's
``_n_drives`` so a C11 window never reads a C3 variable.
``FEM_EM_GEO38_C3_DEGREE`` set in a selected window ⇒ ``ValueError`` before gmsh
(§9 item 706, C17): the degree-2 3 T path is C5e's.
``FEM_EM_GEO38_C11_DEGREE`` (§9 item 869, C5e; C3's ``_degree_knob_c3`` body):
unset / empty ⇒ degree 1, ``2`` ⇒ C3's ``_apply_degree_c3`` swaps in the
degree-2 solver and asserts (Deg) against the 3 T record 2E + 2F, and
``_save_path()`` ⇒ ``C11_SAVE_PATH_D2``, through patches 8 and 9.  C3's own knob
stays refused first, so a stray 1.5 T knob never routes a 3 T build.

Imports are module objects inside function bodies, never ``test_*`` names, so
C3a's, C8's and WF-7's gated tests are never collected from here.

Not set here: the gel's 128 MHz σ / ε_r (the full branch would use the
frequency-independent ``SALINE_SIGMA`` / ``SALINE_EPSILON_R``; the weekly sets
them with a cited source). Build-only never reads them.  Nor, since C5e: any
solve or ``_sweep`` at 3 T degree 2, or a 3 T degree-2 price.

Scope: the module, its coding tests and a build-only proof. No S entry,
reciprocity, passivity, class reading, 3 T column, resonance or tuning claim.
C5e adds the degree knob, its save path and the degree-2 count, build-only.
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

C11_ENV = "FEM_EM_GEO38_C11"
C11_DRIVES_ENV = "FEM_EM_GEO38_C11_DRIVES"
C3_DEGREE_ENV = "FEM_EM_GEO38_C3_DEGREE"
C11_DEGREE_ENV = "FEM_EM_GEO38_C11_DEGREE"
C11_SAVE_PATH = Path("logs/geo38/c11-cu-s64-d1-128mhz.npz")
C11_SAVE_PATH_D2 = Path("logs/geo38/c11-cu-s64-d2-128mhz.npz")
N_PORTS_C11 = 64
OMEGA_REL_BAND = 1.0e-12  # H17's own ω check (`test_wf10_f_human_copper_esr_root.py:118`)

# Printed predictions (brief, "Predictions, written first"), never asserted.
PREDICTED_T_BUILD_S = (313.0, (290.0, 560.0))
PREDICTED_CHECKS_S = (3.0, 5.0)
PREDICTED_WINDOW_S = 325.0
PREDICTED_RSS_SUM_GIB = 2.1
PREDICTED_AZ_RESIDUAL_DEG = 1e-13
F_HUMAN_N1 = 369_387
F_HUMAN_U1 = 500_568
F_HUMAN_COL_S = 39.2
ENTRY40_COLUMN_S = 364.0

_REC: dict = {}


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-38 C11] {msg}", flush=True)


def _mode():
    raw = os.environ.get(C11_ENV, "").strip()
    if raw == "":
        return None
    if raw in ("build", "1"):
        return raw
    raise ValueError(f"{C11_ENV} must be unset, 'build' or '1', got {raw!r}")


def _n_drives():
    raw = os.environ.get(C11_DRIVES_ENV, "").strip()
    n = int(raw) if raw else N_PORTS_C11
    if not 1 <= n <= N_PORTS_C11:
        raise ValueError(f"{C11_DRIVES_ENV} must be in 1..{N_PORTS_C11}, got {raw!r}")
    return n


def _degree_knob_c11():
    """C3's ``_degree_knob_c3`` body on ``C11_DEGREE_ENV``: unset / empty ⇒
    ``None`` (degree 1); ``"2"`` ⇒ 2; anything else ⇒ ValueError."""
    raw = os.environ.get(C11_DEGREE_ENV, "").strip()
    if raw == "":
        return None
    if raw == "2":
        return 2
    raise ValueError(f"{C11_DEGREE_ENV}={raw!r}: only unset/empty (degree 1) or '2' are accepted")


def _refuse_c3_degree():
    """C3's degree knob is refused in a C11 window (§9 item 706, C17; reworded
    by item 869, C5e): it would route ``_apply_degree_c3`` and ``_save_path``
    through a 3 T build; the 3 T degree-2 path is ``C11_DEGREE_ENV``."""
    raw = os.environ.get(C3_DEGREE_ENV, "").strip()
    if raw != "":
        raise ValueError(f"{C3_DEGREE_ENV}={raw!r}: refused in a C11 window; the 3 T degree-2 "
                         f"path is {C11_DEGREE_ENV} (C5e)")


def _modules():
    import tests.validation.test_geo38_f_clinical_3t as t3
    import tests.validation.test_geo38_f_clinical_build as s1a
    import tests.validation.test_geo38_f_clinical_column as col
    import tests.validation.test_geo38_f_clinical_matrix as c3
    import tests.validation.test_wf7_f_human_b1_quadrature as wf7

    return c3, col, s1a, t3, wf7


def _c11_patches(monkeypatch):
    """The geometry knob at 128 MHz: rebind the globals read at call time."""
    from fem_em_solver.io.mesh import MeshGenerator

    c3, col, s1a, t3, wf7 = _modules()
    _REC.clear()
    monkeypatch.setattr(s1a, "_f_clinical_kwargs", t3._f_clinical_3t_kwargs)
    monkeypatch.setattr(s1a, "F_CLINICAL_R1_CELL_RECORD", t3.F_CLINICAL_3T_R1_CELL_RECORD)
    monkeypatch.setattr(s1a, "F_CLINICAL_R1_UNKNOWNS_RECORD", t3.F_CLINICAL_3T_R1_UNKNOWNS_RECORD)
    monkeypatch.setattr(wf7, "FREQUENCY_HZ", t3.FREQUENCY_3T_HZ)
    monkeypatch.setattr(c3, "SAVE_PATH", C11_SAVE_PATH)
    monkeypatch.setattr(c3, "_degree_knob_c3", _degree_knob_c11)
    monkeypatch.setattr(c3, "SAVE_PATH_D2", C11_SAVE_PATH_D2)

    orig_ctx = col._context_from_mesh

    def ctx_recorder(*a, **k):
        out = orig_ctx(*a, **k)
        _REC["omega"] = complex(out[0]["omega"])
        return out

    monkeypatch.setattr(col, "_context_from_mesh", ctx_recorder)

    orig_domain = MeshGenerator.birdcage_port_domain

    def domain_recorder(*a, **k):
        out = orig_domain(*a, **k)
        _REC["kw"] = dict(k)
        _REC["diag"] = out[3]
        return out

    monkeypatch.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(domain_recorder))


class _Tripwire(Exception):
    pass


def _tripwire_kwargs(monkeypatch):
    """Replace the builder with a recorder that raises before gmsh; return the
    kwargs ``col._build_f_clinical`` reached it with."""
    from fem_em_solver.io.mesh import MeshGenerator

    _c3, col, _s1a, _t3, _wf7 = _modules()
    seen = {}

    def recorder(*a, **k):
        seen.update(k)
        raise _Tripwire

    monkeypatch.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(recorder))
    with pytest.raises(_Tripwire):
        col._build_f_clinical(MPI.COMM_WORLD)
    return seen


def _price_lines(comm, cells, u1):
    rho3 = cells / F_HUMAN_N1
    col_s = F_HUMAN_COL_S * (u1 / F_HUMAN_U1) ** 2
    _say(comm, f"PRICE INPUT (printed only): rho3 = {cells} / {F_HUMAN_N1} = {rho3:.4f}; conservative "
               f"column {F_HUMAN_COL_S}*(U1/{F_HUMAN_U1})^2 = {col_s:.1f} s (entry 40's 1.5 T column "
               f"{ENTRY40_COLUMN_S:.0f} s; entry 40 recalibrates both)")


# ---------------------------------------------------------------------------
# Coding tests; skip whenever FEM_EM_GEO38_C11 is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(os.environ.get(C11_ENV, "").strip() != "",
                             reason=f"{C11_ENV} set: coding tests skip in windows")


@_coding
def test_c11_k0_knob_parsing(monkeypatch):
    monkeypatch.delenv(C11_ENV, raising=False)
    assert _mode() is None
    for good in ("build", "1"):
        monkeypatch.setenv(C11_ENV, good)
        assert _mode() == good
    for bad in ("2", "x"):
        monkeypatch.setenv(C11_ENV, bad)
        with pytest.raises(ValueError):
            _mode()
    monkeypatch.delenv(C11_DRIVES_ENV, raising=False)
    assert _n_drives() == 64
    for good in ("1", "7", "64"):
        monkeypatch.setenv(C11_DRIVES_ENV, good)
        assert _n_drives() == int(good)
    for bad in ("0", "65"):
        monkeypatch.setenv(C11_DRIVES_ENV, bad)
        with pytest.raises(ValueError):
            _n_drives()
    # A C3 variable never steers a C11 window.
    monkeypatch.setenv("FEM_EM_GEO38_C3_DRIVES", "5")
    monkeypatch.delenv(C11_DRIVES_ENV, raising=False)
    assert _n_drives() == 64


@_coding
def test_c11_k1_routing_tripwire():
    """(K1): under the patches the column module's builder reaches gmsh with
    the 3 T table; on exit the 1.5 T table and 64 MHz are back."""
    _c3, _col, s1a, _t3, wf7 = _modules()
    with pytest.MonkeyPatch.context() as mp:
        _c11_patches(mp)
        assert wf7.FREQUENCY_HZ == 128e6
        kw = _tripwire_kwargs(mp)
    print(f"\n[GEO-38 C11 K1] tripwire kwargs: ring_radius {kw.get('ring_radius')}, shield_radius "
          f"{kw.get('shield_radius')}, coil_length {kw.get('coil_length')}", flush=True)
    assert abs(kw["ring_radius"] - 0.352) <= 1e-15, kw["ring_radius"]
    assert abs(kw["shield_radius"] - 0.3715) <= 1e-15, kw["shield_radius"]
    assert s1a._f_clinical_kwargs("r1")["ring_radius"] == 0.355
    assert wf7.FREQUENCY_HZ == 64e6


@_coding
def test_c11_n1_15t_record_misses_3t_mesh():
    """(N1), asserted by construction: a patch that silently failed would turn
    (A0) red in W1 — the 1.5 T record is 7.8 bands off the 3 T one."""
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND

    _c3, _col, s1a, t3, _wf7 = _modules()
    rel = abs(t3.F_CLINICAL_3T_R1_CELL_RECORD / s1a.F_CLINICAL_R1_CELL_RECORD - 1.0)
    print(f"\n[GEO-38 C11 N1] |{t3.F_CLINICAL_3T_R1_CELL_RECORD}/{s1a.F_CLINICAL_R1_CELL_RECORD} - 1| = "
          f"{rel:.4e} = {rel / CELL_COUNT_BAND:.2f}x CELL_COUNT_BAND", flush=True)
    assert abs(rel - 7.80e-2) <= 5e-4, rel
    assert rel > 7.0 * CELL_COUNT_BAND, rel
    assert t3.F_CLINICAL_3T_R1_UNKNOWNS_RECORD[0] != s1a.F_CLINICAL_R1_UNKNOWNS_RECORD[0]


@_coding
def test_c11_n2_unpatched_tripwire_is_15t():
    """(N2), K1's twin: without the patches the same call reaches 0.355 m."""
    with pytest.MonkeyPatch.context() as mp:
        kw = _tripwire_kwargs(mp)
    print(f"\n[GEO-38 C11 N2] unpatched tripwire ring_radius {kw.get('ring_radius')}", flush=True)
    assert kw["ring_radius"] == 0.355, kw["ring_radius"]


@_coding
def test_c11_n3_save_path_guard():
    c3, _col, _s1a, _t3, _wf7 = _modules()
    original = c3.SAVE_PATH
    assert C11_SAVE_PATH != original
    with pytest.MonkeyPatch.context() as mp:
        _c11_patches(mp)
        assert c3.SAVE_PATH == C11_SAVE_PATH
    assert c3.SAVE_PATH == original
    print(f"\n[GEO-38 C11 N3] C11 saves to {C11_SAVE_PATH}, never C3's {original}", flush=True)


@_coding
def test_c11_k3t_refuses_c3_degree(monkeypatch):
    """(K-3T): the helper refuses the knob; unset and empty pass; a selected
    window refuses it before ``_run_c3``, so before gmsh."""
    c3, _col, _s1a, _t3, _wf7 = _modules()
    assert C3_DEGREE_ENV == c3.C3_DEGREE_ENV
    monkeypatch.delenv(C3_DEGREE_ENV, raising=False)
    _refuse_c3_degree()
    monkeypatch.setenv(C3_DEGREE_ENV, "")
    _refuse_c3_degree()
    monkeypatch.setenv(C3_DEGREE_ENV, "2")
    with pytest.raises(ValueError, match=r"FEM_EM_GEO38_C3_DEGREE.*refused.*FEM_EM_GEO38_C11_DEGREE \(C5e\)"):
        _refuse_c3_degree()

    def never(*a, **k):
        raise AssertionError("_run_c3 reached past the refusal")

    monkeypatch.setattr(c3, "_run_c3", never)
    monkeypatch.setenv(C11_ENV, "build")
    with pytest.raises(ValueError, match="C5e"):
        test_geo38_c11_3t_matrix(monkeypatch)
    print("\n[GEO-38 C11 K-3T] knob refused in the helper and in the selected path", flush=True)


@_coding
def test_c11_k0c_degree_knob_parsing(monkeypatch):
    """(K0c, §9 item 869): C3's knob grammar on the C11 variable."""
    monkeypatch.delenv(C11_DEGREE_ENV, raising=False)
    assert _degree_knob_c11() is None
    monkeypatch.setenv(C11_DEGREE_ENV, "")
    assert _degree_knob_c11() is None
    monkeypatch.setenv(C11_DEGREE_ENV, "2")
    assert _degree_knob_c11() == 2
    for bad in ("1", "x"):
        monkeypatch.setenv(C11_DEGREE_ENV, bad)
        with pytest.raises(ValueError):
            _degree_knob_c11()


@_coding
def test_c11_k5_degree2_save_path_route(monkeypatch):
    """(K5, §9 item 869): under the patches the C11 knob routes C3's
    ``_save_path`` to C11's own degree-2 file; unset keeps ``C11_SAVE_PATH``;
    both C3 names are restored on exit."""
    c3, _col, _s1a, _t3, _wf7 = _modules()
    orig_knob, orig_d1, orig_d2 = c3._degree_knob_c3, c3.SAVE_PATH, c3.SAVE_PATH_D2
    monkeypatch.delenv(C3_DEGREE_ENV, raising=False)
    with pytest.MonkeyPatch.context() as mp:
        _c11_patches(mp)
        mp.setenv(C11_DEGREE_ENV, "2")
        d2 = c3._save_path()
        assert d2 == C11_SAVE_PATH_D2, d2
        assert c3._degree_knob_c3() == 2
        assert d2 not in (orig_d1, orig_d2, C11_SAVE_PATH), d2
        mp.delenv(C11_DEGREE_ENV)
        assert c3._save_path() == C11_SAVE_PATH
    assert c3._degree_knob_c3 is orig_knob
    assert c3.SAVE_PATH == orig_d1 and c3.SAVE_PATH_D2 == orig_d2
    print(f"\n[GEO-38 C11 K5] C11 d2 {C11_SAVE_PATH_D2}; C11 d1 {C11_SAVE_PATH}; C3 d1 {orig_d1}, "
          f"d2 {orig_d2}; restored on exit", flush=True)


@_coding
def test_c11_n4_deg_refuses_degree_1_count():
    """(N4, §9 item 869; asserted, integers only): (Deg) against the 3 T 2E + 2F
    refuses the degree-1 size; the separation is U2 / U1 = 5.145."""
    c3, _col, _s1a, t3, _wf7 = _modules()
    u1, u2 = t3.F_CLINICAL_3T_R1_UNKNOWNS_RECORD
    n_edges, n_faces = 1_425_765, 2_242_229
    assert n_edges == u1
    assert c3._check_unknowns(u2, 2, n_edges, n_faces) == u2 == 7_335_988
    with pytest.raises(AssertionError):
        c3._check_unknowns(1_425_765, 2, 1_425_765, 2_242_229)
    print(f"\n[GEO-38 C11 N4] U2 / U1 = {u2 / u1:.4f} (> 5, ASSERTED)", flush=True)
    assert u2 / u1 > 5.0, (u1, u2)


@_coding
def test_c11_n5_15t_degree2_record_misses_3t_count():
    """(N5, §9 item 869): a silently failed patch 3 would assert the 1.5 T
    degree-2 count on the 3 T mesh; the two records are 7.49e-02 apart."""
    _c3, _col, s1a, t3, _wf7 = _modules()
    u2_15, u2_3 = s1a.F_CLINICAL_R1_UNKNOWNS_RECORD[1], t3.F_CLINICAL_3T_R1_UNKNOWNS_RECORD[1]
    rel = abs(u2_15 / u2_3 - 1.0)
    print(f"\n[GEO-38 C11 N5] |{u2_15}/{u2_3} - 1| = {rel:.4e} (> 7e-02, ASSERTED)", flush=True)
    assert u2_15 == 7_885_674 and u2_3 == 7_335_988
    assert u2_15 != u2_3
    assert rel > 7.0e-2, rel


@_coding
def test_c11_n6_c11_knob_does_not_leak_into_15t(monkeypatch):
    """(N6, §9 item 869; item 705's flag (2) in reverse): without the patches
    the C11 knob is invisible to C3's route."""
    c3, _col, _s1a, _t3, _wf7 = _modules()
    monkeypatch.delenv(C3_DEGREE_ENV, raising=False)
    monkeypatch.setenv(C11_DEGREE_ENV, "2")
    assert c3._degree_knob_c3() is None
    assert c3._save_path() == c3.SAVE_PATH
    print(f"\n[GEO-38 C11 N6] {C11_DEGREE_ENV}=2 unpatched: C3 knob None, save path {c3.SAVE_PATH}",
          flush=True)


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
@complex_only
def test_geo38_c11_3t_matrix(monkeypatch):
    mode = _mode()
    if mode is None:
        pytest.skip(f"{C11_ENV} unset: GEO-38 C11 is selected by environment only")
    _refuse_c3_degree()
    c3, _col, s1a, t3, _wf7 = _modules()
    comm = MPI.COMM_WORLD
    c3_save = c3.SAVE_PATH
    c3_state = (c3_save.exists(), c3_save.stat().st_mtime_ns if c3_save.exists() else None)
    c3_save_d2 = c3.SAVE_PATH_D2  # read before _c11_patches rebinds it (C5e)
    tb, tb_range = PREDICTED_T_BUILD_S
    _say(comm, f"PREDICTED (printed, never asserted): (A0) rel 0; t_build ≈{tb} s {tb_range} at loadavg "
               f"<= 4; checks ≈{PREDICTED_CHECKS_S[0]}-{PREDICTED_CHECKS_S[1]} s; window ≈"
               f"{PREDICTED_WINDOW_S} s; summed RSS ≈{PREDICTED_RSS_SUM_GIB} GiB; (Az) worst "
               f"residual ≈{PREDICTED_AZ_RESIDUAL_DEG:.0e} deg")

    _c11_patches(monkeypatch)
    if mode == "build":
        save_path = c3._save_path()  # after the patches: C11's d1 or d2 path (C5e)
        guarded = (C11_SAVE_PATH, c3_save_d2)
        guarded_before = [comm.bcast(c3._size_mtime(p) if comm.rank == 0 else None, root=0)
                          for p in guarded]
        out, c11_state = c3._build_keeps_save_state(comm, save_path,
                                                    lambda: c3._run_c3(comm, build_only=True))
        guarded_after = [comm.bcast(c3._size_mtime(p) if comm.rank == 0 else None, root=0)
                         for p in guarded]
    else:
        out = c3._run_c3(comm, build_only=False, n_drives=_n_drives())

    # (f) the frequency reached the fixture.
    omega = _REC["omega"]
    omega_ref = 2.0 * np.pi * t3.FREQUENCY_3T_HZ
    f_rel = abs(omega / omega_ref - 1.0)
    _say(comm, f"(f) recorded omega {omega.real:.12e} vs 2*pi*128e6 {omega_ref:.12e}: rel {f_rel:.3e} "
               f"(band {OMEGA_REL_BAND:.0e}, ASSERTED)")
    assert f_rel <= OMEGA_REL_BAND, (omega, omega_ref)

    # (T3) the 3 T table reached the builder.
    diag, kw = _REC["diag"], _REC["kw"]
    band = s1a.IDENTITY_BAND
    shield_d = 2.0 * float(diag["shield_radius_m"])
    shield_l = float(diag["shield_length_m"])
    rung_d = 2.0 * kw["ring_radius"]
    _say(comm, f"(T3) shield Ø {shield_d:.12f} (vs {t3.SHIELD_DIAMETER_3T_M}), length {shield_l:.12f} "
               f"(vs {t3.SHIELD_LENGTH_3T_M}), rung Ø {rung_d:.12f} (vs {t3.RUNG_DIAMETER_3T_M}) "
               f"(IDENTITY_BAND {band:.0e}, ASSERTED)")
    assert abs(shield_d / t3.SHIELD_DIAMETER_3T_M - 1.0) <= band, shield_d
    assert abs(shield_l / t3.SHIELD_LENGTH_3T_M - 1.0) <= band, shield_l
    assert abs(rung_d / t3.RUNG_DIAMETER_3T_M - 1.0) <= band, rung_d

    if mode != "build":
        return
    _price_lines(comm, out["cells"], out["u1"])
    c3_after = (c3_save.exists(), c3_save.stat().st_mtime_ns if c3_save.exists() else None)
    _say(comm, f"save paths after build-only: (S) C11 {save_path} (size, mtime_ns) or absence "
               f"unchanged: {c11_state} (ASSERTED); C3 {c3_save} (exists, mtime) "
               f"{c3_after} == before {c3_state} (ASSERTED)")
    for p, before, after in zip(guarded, guarded_before, guarded_after):
        _say(comm, f"(S) {p} (size, mtime_ns) {after} == before {before} (ASSERTED)")
    assert c3_after == c3_state, (c3_after, c3_state)
    assert guarded_after == guarded_before, (guarded, guarded_after, guarded_before)
    if save_path == C11_SAVE_PATH_D2:
        assert c11_state is None, (save_path, c11_state)  # absent before, and so after
