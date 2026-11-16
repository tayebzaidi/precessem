"""GEO-36 step 11 / §10 C32 (§9 item 822): I36-3T's body — I36b's arm M / arm R
route (`test_geo36_rod_f_clinical_drive.py`, `_run_drive`) composed with C21's
3 T patch set (`test_geo38_f_clinical_3t_box_loaded.py`, `_c21_patches`) on
C20's 3 T rod mesh, stub-run (no solve, no gmsh).

Brief: ``docs/planning/chunks/GEO-36.md`` "§9 item 822" and the "Review ruling
2026-11-05 03:00" paragraphs above it and in ``GEO-38.md``.

Composition, never an edit: I36b, C21, C20 and ``src`` are byte-identical.
`_run_drive_3t` applies `c21._c21_patches` (C20's mesh route, C20's records
read by `i36._route_drive` at call time, 128 MHz on ``wf7`` / ``c4``) and then
calls `i36._run_drive(comm, monkeypatch, degree, c_pf)` — degree 1 unless
``FEM_EM_GEO36_S11_DRIVE_DEGREE=2`` (C52, §9 item 929).
Disclosed: `_c21_patches`' patch 4 wraps ``box._box_readouts``, which I36's
route never calls, so it is inert here; and I36b's Y_in print still reads "the
coil was tuned on C4's cylinder", a stale label at 3 T (not edited).  The tag-403
σ_wire cells, the node knob and the two-port drive on C16's feeds P33 / P41
(``c4._feeds_32``, the same at both field strengths) are I36b's, unchanged.

Selected by ``FEM_EM_GEO36_S11_DRIVE`` only (never by ``-k``): unset ⇒ the gated
test skips and the coding tests run; ``1`` ⇒ `_run_drive_3t` (node knob
``FEM_EM_GEO36_S8_DRIVE_C_PF`` required, read through ``i36._node_knob``);
anything else ⇒ ``ValueError``.  ``FEM_EM_GEO36_S8_DRIVE_DEGREE`` is in the
foreign list, so a set 1.5 T degree knob is refused by the foreign guard; the
3 T degree is routed only by C52's ``FEM_EM_GEO36_S11_DRIVE_DEGREE`` (unset /
empty ⇒ degree 1, ``2`` ⇒ 2, else ``ValueError``).  At degree 2 a pass-through
wrapper on ``c4._build_c4_fixture`` asserts (Deg) n_dofs = C20's U₂ record.

Print-only beside the 3 T values: entry 48's 1.5 T tip Q 4.341688–4.460568
(``20261102T080008Z_GEO-36-step9-I36.log:12501``) and η_rod(T1) 3.642.

Scope: the body and its stub proof.  Not claimed: any 128 MHz tip Q or rod
reading (``xl`` entry 61), the second factorisation, a published comparand.
Box containment is not asserted on this route.
"""

from __future__ import annotations

import os

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

S11_DRIVE_ENV = "FEM_EM_GEO36_S11_DRIVE"
#: C52 (§9 item 929): the 3 T degree knob (`i36._drive_degree`'s semantics).
S11_DRIVE_DEGREE_ENV = "FEM_EM_GEO36_S11_DRIVE_DEGREE"
#: The 3 T node (``xl`` entry 49 (C16): C_in(F-clinical-3T, d1) = 44.0615 pF).
K_NODE_3T_PF = 44.0615
#: Entry 48's 1.5 T readings (print-only comparands, never asserted).
ENTRY48_Q_RANGE = (4.341688, 4.460568)
ENTRY48_ETA_T1 = 3.642
#: Added to I36b's FOREIGN_ENVS (brief: I36b's drive / degree knobs, C6 / C21 /
#: C20 and the C11 / C12 knobs).
FOREIGN_ENVS_C32_EXTRA = (
    "FEM_EM_GEO36_S8_DRIVE",
    "FEM_EM_GEO36_S8_DRIVE_DEGREE",
    "FEM_EM_GEO38_C6",
    "FEM_EM_GEO38_C21",
    "FEM_EM_GEO38_C21_DEGREE",
    "FEM_EM_GEO38_C20",
    "FEM_EM_GEO38_C11",
    "FEM_EM_GEO38_C11_DRIVES",
    "FEM_EM_GEO38_C12",
    "FEM_EM_GEO38_C12_NODES_PF",
    "FEM_EM_GEO38_C12_DEGREE",
)
#: The foreign guard's message head (K1 matches it).
FOREIGN_MSG = "C32 never reads these variables"


def _modules():
    import tests.validation.test_geo36_rod_f_clinical as i35
    import tests.validation.test_geo36_rod_f_clinical_drive as i36
    import tests.validation.test_geo38_f_clinical_3t as t3
    import tests.validation.test_geo38_f_clinical_3t_box_loaded as c21
    import tests.validation.test_geo38_f_clinical_3t_rod_box as c20
    import tests.validation.test_geo38_f_clinical_box as s2
    import tests.validation.test_geo38_f_clinical_build as s1a
    import tests.validation.test_geo38_f_clinical_esr_root as c4
    import tests.validation.test_wf7_f_human_b1_quadrature as wf7

    return {"i35": i35, "i36": i36, "t3": t3, "c21": c21, "c20": c20, "s2": s2, "s1a": s1a,
            "c4": c4, "wf7": wf7}


def _foreign_envs():
    return tuple(_modules()["i36"].FOREIGN_ENVS) + FOREIGN_ENVS_C32_EXTRA


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-36 C32] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Knobs.
# ---------------------------------------------------------------------------
def _guard_foreign():
    bad = [k for k in _foreign_envs() if os.environ.get(k) is not None]
    if bad:
        raise ValueError(f"{FOREIGN_MSG}; unset them: {bad}")


def _s11_mode():
    raw = os.environ.get(S11_DRIVE_ENV, "").strip()
    if raw == "":
        return None
    if raw == "1":
        return raw
    raise ValueError(f"{S11_DRIVE_ENV} must be unset or '1', got {raw!r}")


def _knobs():
    """Foreign guard, then ``(mode, c_pf)`` with I36b's node-knob semantics."""
    _guard_foreign()
    mode = _s11_mode()
    return mode, _modules()["i36"]._node_knob(mode)


def _s11_degree():
    """C52: unset / empty ⇒ ``None`` (degree 1); ``"2"`` ⇒ 2; anything else ⇒
    ValueError (`i36._drive_degree`'s semantics under ``S11_DRIVE_DEGREE_ENV``)."""
    raw = os.environ.get(S11_DRIVE_DEGREE_ENV, "").strip()
    if raw == "":
        return None
    if raw == "2":
        return 2
    raise ValueError(f"{S11_DRIVE_DEGREE_ENV}={raw!r}: only unset/empty (degree 1) or '2' are accepted")


# ---------------------------------------------------------------------------
# The body.
# ---------------------------------------------------------------------------
def _run_drive_3t(comm, monkeypatch, c_pf, *, degree=None):
    """C21's 3 T patch set, then I36b's `_run_drive` at ``degree`` (None / 1 or 2;
    arm M then arm R; every gate of I36b's asserted inside).  At degree 2 a
    pass-through wrapper on ``c4._build_c4_fixture`` prints and asserts (Deg)
    n_dofs = C20's U₂ record (C6a's `:566–570` line, which I36b's `_run_drive`
    lacks).  On return, entry 48's 1.5 T Q / η_rod(T1) are printed beside this
    node's JSON values (print-only)."""
    import json

    deg = degree or 1
    if deg not in (1, 2):
        raise ValueError(f"degree {degree!r}: only None / 1 or 2")
    m = _modules()
    i36, t3, c4, wf7, i35 = m["i36"], m["t3"], m["c4"], m["wf7"], m["i35"]
    m["c21"]._c21_patches(monkeypatch)
    _say(comm, f"3 T patch set (C21's): rung Ø {t3.RUNG_DIAMETER_3T_M} m, f {wf7.FREQUENCY_HZ:.6e} Hz "
               f"(C4 {c4.FREQUENCY_HZ_C4:.6e}); records cells {i35.F_CLINICAL_ROD_CELL_RECORD}, U "
               f"{i35.F_CLINICAL_ROD_UNKNOWNS_RECORD}; node {c_pf:.4f} pF, degree {deg}; patch 4 "
               f"(box._box_readouts) is inert on this route (I36 never calls it)")
    if deg == 2:
        inner = c4._build_c4_fixture

        def _fixture_with_deg(comm_):
            fx = inner(comm_)
            n2 = fx["counts"]["n_dofs"]
            _say(comm_, f"(Deg) degree 2 global unknowns {n2} vs record "
                        f"{i35.F_CLINICAL_ROD_UNKNOWNS_RECORD[1]} (EXACT, ASSERTED)")
            assert n2 == i35.F_CLINICAL_ROD_UNKNOWNS_RECORD[1], (n2, i35.F_CLINICAL_ROD_UNKNOWNS_RECORD)
            return fx

        monkeypatch.setattr(c4, "_build_c4_fixture", _fixture_with_deg)
    i36._run_drive(comm, monkeypatch, degree, c_pf)
    rec = None
    if comm.rank == 0:
        path = i36._i36_json_path(c_pf, deg)
        with open(path) as fh:
            rec = json.load(fh)
    rec = comm.bcast(rec, root=0)
    _say(comm, "3 T tip 1 g ratio Q: " + ", ".join(f"{n} {v:.6f}" for n, v in rec["q"].items())
               + f"; eta_rod(T1) {rec['eta_t1']:.6f} — beside entry 48's 1.5 T Q "
               f"{ENTRY48_Q_RANGE[0]}-{ENTRY48_Q_RANGE[1]} and eta_rod(T1) {ENTRY48_ETA_T1} (printed, "
               f"never asserted)")
    return rec


@complex_only
def test_geo36_c32_rod_drive_3t(monkeypatch):
    mode, c_pf = _knobs()  # node knob missing under 1 => ValueError, before gmsh
    if mode is None:
        pytest.skip(f"{S11_DRIVE_ENV} unset: GEO-36 C32 is selected by environment only")
    _run_drive_3t(MPI.COMM_WORLD, monkeypatch, c_pf, degree=_s11_degree())


# ---------------------------------------------------------------------------
# Coding tests (no gmsh); skip whenever S11_DRIVE_ENV is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(os.environ.get(S11_DRIVE_ENV, "").strip() != "",
                             reason=f"{S11_DRIVE_ENV} set: coding tests skip in windows")


@_coding
def test_c32_k1_knob_parsing(monkeypatch):
    """(K1): mode unset / 1 accepted, else ValueError; node knob required under 1
    and read through I36b's `_node_knob`; every foreign variable — the degree knob
    included — is refused by the foreign guard (that is what raises for it)."""
    i36 = _modules()["i36"]
    foreign = _foreign_envs()
    for k in (S11_DRIVE_ENV, i36.S8_DRIVE_C_PF_ENV) + foreign:
        monkeypatch.delenv(k, raising=False)
    assert _knobs() == (None, None)
    monkeypatch.setenv(S11_DRIVE_ENV, "1")
    with pytest.raises(ValueError, match="required"):
        _knobs()
    monkeypatch.setenv(i36.S8_DRIVE_C_PF_ENV, f"{K_NODE_3T_PF}")
    assert _knobs() == ("1", K_NODE_3T_PF)
    for bad in ("2", "build", "x"):
        monkeypatch.setenv(S11_DRIVE_ENV, bad)
        with pytest.raises(ValueError, match=S11_DRIVE_ENV):
            _s11_mode()
    monkeypatch.setenv(S11_DRIVE_ENV, "1")
    for k in foreign:
        with monkeypatch.context() as m:
            m.setenv(k, "1")
            with pytest.raises(ValueError, match=FOREIGN_MSG) as exc:
                _knobs()
            assert k in str(exc.value), (k, str(exc.value))
    # The degree knob: refused by the foreign guard at either value I36b accepts.
    for val in ("2", ""):
        with monkeypatch.context() as m:
            m.setenv(i36.S8_DRIVE_DEGREE_ENV, val)
            with pytest.raises(ValueError, match=FOREIGN_MSG) as exc:
                _knobs()
            assert i36.S8_DRIVE_DEGREE_ENV in str(exc.value)
    assert set(i36.FOREIGN_ENVS) <= set(foreign)
    assert i36.S8_DRIVE_DEGREE_ENV in foreign and i36.S8_DRIVE_ENV in foreign
    assert i36.S8_DRIVE_C_PF_ENV not in foreign
    print(f"\n[GEO-36 C32 K1] mode unset/1 accepted, 2/build/x refused; node knob required under 1, "
          f"read {K_NODE_3T_PF}; {len(foreign)} foreign variables refused by the guard (I36b's "
          f"{len(i36.FOREIGN_ENVS)} + {FOREIGN_ENVS_C32_EXTRA}); degree knob '2' / '' refused by the "
          f"foreign guard (ASSERTED)")


@_coding
def test_c32_k2_one_without_node_stops_before_gmsh(monkeypatch):
    """(K2): ``1`` without the node knob raises ValueError naming it, never gmsh."""
    from fem_em_solver.io.mesh import MeshGenerator

    i36 = _modules()["i36"]

    def _sentinel(*a, **k):
        raise AssertionError("gmsh reached")

    with monkeypatch.context() as m:
        for k in _foreign_envs() + (i36.S8_DRIVE_C_PF_ENV,):
            m.delenv(k, raising=False)
        m.setenv(S11_DRIVE_ENV, "1")
        m.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(_sentinel))
        with pytest.raises(ValueError, match=i36.S8_DRIVE_C_PF_ENV):
            test_geo36_c32_rod_drive_3t(m)
    print("\n[GEO-36 C32 K2] S11=1 without the node knob raised ValueError before gmsh (ASSERTED)")


@_coding
def test_c32_k7_one_with_node_reaches_builder_with_3t_kwargs(monkeypatch):
    """(K7-3T): ``1`` + node 44.0615 runs `_run_drive_3t` into C21's
    `_builder_sentinel`: rung Ø 0.704, C20's rod with ``keep_cells=True``, f 128 MHz
    on both globals, C20's records; (Bore) gap > 0 under the patch.  Negative
    (asserted, `20261103T103420Z_GEO-38-C21-W0.log:83`): I36b's `_run_drive` without
    C21's patches reaches rung Ø 0.710 / 64 MHz / 1 217 025 cells."""
    mods = _modules()
    t3, c20, c21, i35, i36 = mods["t3"], mods["c20"], mods["c21"], mods["i35"], mods["i36"]
    kw15 = i35._kwargs()  # the 1.5 T rod kwargs, unpatched
    kw3 = c20._kwargs_3t()
    comm = MPI.COMM_WORLD

    seen = []
    with monkeypatch.context() as m:
        for k in _foreign_envs():
            m.delenv(k, raising=False)
        m.setenv(S11_DRIVE_ENV, "1")
        m.setenv(i36.S8_DRIVE_C_PF_ENV, f"{K_NODE_3T_PF}")
        c21._builder_sentinel(m, seen)
        with pytest.raises(c21._BuilderReached):
            test_geo36_c32_rod_drive_3t(m)
        # (Bore) under the patch set (still applied inside this context).
        gap, corner = c20._bore_gap(i35._kwargs())

    neg = []
    with monkeypatch.context() as m:
        for k in _foreign_envs():
            m.delenv(k, raising=False)
        c21._builder_sentinel(m, neg)
        with pytest.raises(c21._BuilderReached):
            i36._run_drive(comm, m, None, K_NODE_3T_PF)

    assert len(seen) == 1 and len(neg) == 1, (len(seen), len(neg))
    got, got_neg = seen[0], neg[0]
    kw, kwn = got["kw"], got_neg["kw"]
    rung_d, rung_d_neg = 2.0 * kw["ring_radius"], 2.0 * kwn["ring_radius"]
    rod = kw.get("implant_rod")
    print(f"\n[GEO-36 C32 K7] S11=1 node {K_NODE_3T_PF:.4f}: builder reached with rung Ø {rung_d:.12f} m "
          f"(3 T {t3.RUNG_DIAMETER_3T_M}); rod {rod}; f wf7 {got['f_wf7']:.6e} / C4 {got['f_c4']:.6e} "
          f"Hz; records cells {got['cells']}, U {got['unknowns']}; (Bore) corner {corner:.9f} m, gap "
          f"{gap:.9f} m (> 0) (ASSERTED)")
    print(f"[GEO-36 C32 K7] negative (I36b's _run_drive, no C21 patches): rung Ø {rung_d_neg:.12f} m, "
          f"f wf7 {got_neg['f_wf7']:.6e} / C4 {got_neg['f_c4']:.6e} Hz, records cells {got_neg['cells']} "
          f"(ASSERTED: 0.710 / 64 MHz / 1 217 025)", flush=True)
    # (K7-3T)
    assert abs(rung_d / t3.RUNG_DIAMETER_3T_M - 1.0) <= 1e-15 and t3.RUNG_DIAMETER_3T_M == 0.704, rung_d
    assert rod is not None and rod == kw3["implant_rod"] and rod["keep_cells"] is True, rod
    assert got["f_wf7"] == t3.FREQUENCY_3T_HZ == 128e6 and got["f_c4"] == 128e6, got
    assert got["cells"] == c20.F_CLINICAL_3T_ROD_CELL_RECORD == 1_128_430
    assert got["unknowns"] == c20.F_CLINICAL_3T_ROD_UNKNOWNS_RECORD == (1_508_256, 7_783_028)
    assert gap > 0.0, (gap, corner)
    # Negative: categorical separation.
    assert abs(rung_d_neg - 0.710) <= 1e-12 and kwn["ring_radius"] == kw15["ring_radius"], rung_d_neg
    assert got_neg["f_wf7"] == 64.0e6 and got_neg["f_c4"] == 64.0e6, got_neg
    assert got_neg["cells"] == 1_217_025, got_neg["cells"]
    # Teardown restored the 1.5 T namespace.
    assert i35.F_CLINICAL_ROD_CELL_RECORD == 1_217_025
    assert mods["wf7"].FREQUENCY_HZ == 64.0e6


def _k9_stubs_3t(m, fx, json_dir, *, p_src_r2=1.0):
    """I36d's `_k9_stubs` (`drive.py:1097`), then the 3 T overrides: env ω =
    2π·128e6, the ESR specs at 44.0615 pF and 128 MHz, and `_drive` on the 128 MHz
    loads (the same seed-656 reciprocal circuit).  Arm R (the solver token) gets
    0.5 E_u.  Returns ``(lines, z_sym, z_rc)``."""
    from types import SimpleNamespace

    from dolfinx import fem

    import fem_em_solver.core as core
    from tests.validation import test_wf10_f_human_copper_esr_root as h17

    mods = _modules()
    i36, c4, t3 = mods["i36"], mods["c4"], mods["t3"]
    comm = MPI.COMM_WORLD
    lines, z_sym = i36._k9_stubs(m, fx, json_dir, p_src_r2=p_src_r2)
    env = dict(i36._drive_env(None, comm))
    env["omega"] = 2.0 * np.pi * t3.FREQUENCY_3T_HZ
    specs, _zrc, _zth = c4._esr_specs()
    z0 = 50.0
    z_rc = h17._z_rc(K_NODE_3T_PF * 1e-12, t3.FREQUENCY_3T_HZ)
    z_th = h17._thevenin_feed_impedance(z_rc, z0)
    loads = h17._parallel_loads(4, [0, 1], z_rc, z_th)
    _msh, _ct, e_u, sigma_f, _tips, _t1 = fx
    e_half = fem.Function(e_u.function_space)
    e_half.x.array[:] = 0.5 * e_u.x.array
    e_half.x.scatter_forward()
    feeds, order = ["P1", "P2"], ["P1", "P2", "P3", "P4"]
    token_r = SimpleNamespace(name="arm R token (3 T)")

    def _drive(_env, solver, _specs, driven_id, *, held):
        j = feeds.index(driven_id)
        arm_r = solver is token_r
        v = np.zeros(4, dtype=np.complex128)
        v[j] = 1.0
        cur = h17._parallel_currents(z_sym, loads, v)
        return {"currents": dict(zip(order, cur)), "kind": "held" if held else "solve",
                "t_solve": 0.0, "t_post": 0.0, "e": (e_half if arm_r else e_u).copy(),
                "sigma": sigma_f, "p_src": p_src_r2 if arm_r and j == 1 else 1.0, "p_sheets": 0.5,
                "p_phantom": 0.2, "p_rest": 0.2, "p_surf": 0.1, "p_rod": 0.0,
                "p_sheet_by_id": {"P1": 0.1, "P2": 0.1, "P3": 0.15, "P4": 0.15}}

    m.setattr(i36, "_drive_env", lambda _fixture, _comm: env)
    m.setattr(c4, "_esr_specs", lambda *_a: (specs, z_rc, z_th))
    m.setattr(core, "TimeHarmonicSolver", lambda *_a, **_k: token_r)
    m.setattr(i36, "_drive", _drive)
    return lines, z_sym, z_rc


def _k9_json_3t(comm, directory):
    import json

    rec = None
    if comm.rank == 0:
        with open(os.path.join(directory, f"i36-node-{K_NODE_3T_PF:.4f}.json")) as fh:
            rec = json.load(fh)
    return comm.bcast(rec, root=0)


@complex_only
@_coding
def test_c32_k9_run_drive_3t_end_to_end_on_stubs(monkeypatch, tmp_path):
    """(K9-3T): `_run_drive_3t` on I36d's `_k9_stubs` with the 3 T overrides.
    (K9-Q) q = 4, (K9-E) η = 2, (K9-A) census, (K9-J) JSON at
    ``i36-node-44.0615.json`` with M's and R's s2 equal to each other and to the
    128 MHz circuit route, (K9-L) heads.  Controls (asserted): arm R column 2's
    P_src × (1 + 1e-3) raises naming ``R (gel twin) (P) P2`` with the JSON still
    written (I36d: 9.990e-04 vs 1e-6); the 64 MHz route misses s2 by > the exact
    band (`20261103T103420Z_GEO-38-C21-W0.log:149`: 2.515e-01)."""
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    mods = _modules()
    i36 = mods["i36"]
    comm = MPI.COMM_WORLD
    fx = i36._k9_fixture(comm)
    msh, cell_tags = fx[0], fx[1]
    root = comm.bcast(str(tmp_path), root=0)
    good_dir, neg_dir = os.path.join(root, "good"), os.path.join(root, "neg")

    with monkeypatch.context() as m:
        lines, z_sym, z_rc = _k9_stubs_3t(m, fx, good_dir)
        ret = _run_drive_3t(comm, m, K_NODE_3T_PF)
        f_seen = (mods["wf7"].FREQUENCY_HZ, mods["c4"].FREQUENCY_HZ_C4)
    rec = _k9_json_3t(comm, good_dir)

    with monkeypatch.context() as m:
        _k9_stubs_3t(m, fx, neg_dir, p_src_r2=1.0 + 1.0e-3)
        with pytest.raises(AssertionError, match=r"R \(gel twin\) \(P\) P2") as neg:
            _run_drive_3t(comm, m, K_NODE_3T_PF)
    rec_neg = _k9_json_3t(comm, neg_dir)

    r_q = {n: abs(v / 4.0 - 1.0) for n, v in rec["q"].items()}
    r_e = abs(rec["eta_t1"] / 2.0 - 1.0)

    def _s(block):
        return np.array([[complex(*v) for v in row] for row in rec[block]["s2"]])

    r_mr = float(np.linalg.norm(_s("M") - _s("R")) / np.linalg.norm(_s("M")))
    s_route = h17._circuit_current_route_s(z_sym, [0, 1], z_rc, z_rc, 50.0)
    r_route = float(np.linalg.norm(_s("M") - s_route) / np.linalg.norm(s_route))
    z_rc_64 = h17._z_rc(K_NODE_3T_PF * 1e-12, 64.0e6)
    s_64 = h17._circuit_current_route_s(z_sym, [0, 1], z_rc_64, z_rc_64, 50.0)
    r_64 = float(np.linalg.norm(_s("M") - s_64) / np.linalg.norm(s_route))

    size_local = msh.topology.index_map(3).size_local
    owned = np.asarray(cell_tags.indices) < size_local
    tag_values = np.asarray(cell_tags.values)[owned]
    census = {str(t): int(comm.allreduce(int(np.count_nonzero(tag_values == t)), op=MPI.SUM))
              for t in (3, 403)}
    a0m_clean = any(ln.startswith("(A0m)") and ln.rstrip().endswith("failures none (ASSERTED exact)")
                    for ln in lines)
    heads = ("(P_RC)", "(A0m)", "tip 1 g ratio", "price:")
    missing_heads = [h for h in heads if not any(ln.startswith(h) for ln in lines)]
    neg_rel = abs(rec_neg["R"]["parts"][1]["p_src"] - 0.5 - 0.4 - 0.1) / rec_neg["R"]["parts"][1]["p_src"]
    if comm.rank == 0:
        print(f"\n[GEO-36 C32 K9] f during the run wf7 / C4 {f_seen}; (K9-Q) misses from 4 {r_q} (band "
              f"{EXACT_IDENTITY_RTOL:g}, ASSERTED); (K9-E) miss from 2 {r_e:.3e} (band "
              f"{DISCRETE_IDENTITY_RTOL:g}, ASSERTED)")
        print(f"[GEO-36 C32 K9] (K9-A) census_m {rec['census_m']} vs the mesh's {census}, (A0m) clean "
              f"{a0m_clean} (ASSERTED); (K9-J) keys {sorted(rec)}, c_pf {rec['c_pf']!r}; s2 M vs R "
              f"{r_mr:.3e}, vs the 128 MHz circuit route {r_route:.3e} (band {EXACT_IDENTITY_RTOL:g}, "
              f"ASSERTED); vs the 64 MHz route {r_64:.3e} (> band, ASSERTED; C21 W0 read 2.515e-01)")
        print(f"[GEO-36 C32 K9] (K9-L) missing {missing_heads or 'none'}; control: raised "
              f"{str(neg.value).splitlines()[0]}; JSON still written, R parts[1] p_src "
              f"{rec_neg['R']['parts'][1]['p_src']!r}, rel {neg_rel:.3e} (I36d 9.990e-04; ASSERTED)",
              flush=True)
    assert f_seen == (128e6, 128e6), f_seen
    assert ret["q"] == rec["q"] and ret["eta_t1"] == rec["eta_t1"]
    assert len(r_q) == len(fx[4]) and max(r_q.values()) <= EXACT_IDENTITY_RTOL, r_q
    assert r_e <= DISCRETE_IDENTITY_RTOL, r_e
    assert rec["census_m"] == census, rec["census_m"]
    assert sum(census.values()) == msh.topology.index_map(3).size_global
    assert census["403"] > 0
    assert a0m_clean
    assert {"M", "R"} <= set(rec) and rec["c_pf"] == K_NODE_3T_PF and rec["degree"] == 1
    assert r_mr <= EXACT_IDENTITY_RTOL, r_mr
    assert r_route <= EXACT_IDENTITY_RTOL, r_route
    assert r_64 > EXACT_IDENTITY_RTOL, r_64
    assert not missing_heads, missing_heads
    assert rec_neg["R"]["parts"][1]["p_src"] == 1.0 + 1.0e-3, rec_neg["R"]["parts"][1]
    assert rec_neg["M"]["parts"][1]["p_src"] == 1.0, rec_neg["M"]["parts"][1]


# ---------------------------------------------------------------------------
# C52 (§9 item 929): the 3 T degree knob, routed on stubs.
# ---------------------------------------------------------------------------
@_coding
def test_c52_k1_s11_degree_parsing(monkeypatch):
    """(K1-S11deg): unset / '' ⇒ None, '2' ⇒ 2, '1' / 'x' / '3' refused naming the
    knob; ``FEM_EM_GEO38_C21_DEGREE`` is refused by the foreign guard, naming it.
    (K-66), C32's half: `_knobs()` with ``FEM_EM_GEO38_C21`` or its degree knob set
    raises naming it, so entry 66 is two pytest calls."""
    i36 = _modules()["i36"]
    for k in (S11_DRIVE_ENV, S11_DRIVE_DEGREE_ENV, i36.S8_DRIVE_C_PF_ENV) + _foreign_envs():
        monkeypatch.delenv(k, raising=False)
    assert _s11_degree() is None
    monkeypatch.setenv(S11_DRIVE_DEGREE_ENV, "")
    assert _s11_degree() is None
    monkeypatch.setenv(S11_DRIVE_DEGREE_ENV, "2")
    assert _s11_degree() == 2
    for bad in ("1", "x", "3"):
        monkeypatch.setenv(S11_DRIVE_DEGREE_ENV, bad)
        with pytest.raises(ValueError, match=S11_DRIVE_DEGREE_ENV):
            _s11_degree()
    monkeypatch.setenv(S11_DRIVE_DEGREE_ENV, "2")
    monkeypatch.setenv(S11_DRIVE_ENV, "1")
    monkeypatch.setenv(i36.S8_DRIVE_C_PF_ENV, f"{K_NODE_3T_PF:.4f}")
    alone = (_knobs(), _s11_degree())  # C32's call alone: accepted
    msgs = {}
    for k, v in (("FEM_EM_GEO38_C21_DEGREE", "2"), ("FEM_EM_GEO38_C21", "1")):
        with monkeypatch.context() as m:
            m.setenv(k, v)
            with pytest.raises(ValueError, match=FOREIGN_MSG) as exc:
                _knobs()
            msgs[k] = str(exc.value)
    if MPI.COMM_WORLD.rank == 0:
        print(f"\n[GEO-36 C52 K1-S11deg] {S11_DRIVE_DEGREE_ENV}: unset / '' / 2 accepted, 1 / x / 3 "
              f"refused; S11=1 + degree 2 + node alone {alone}; C21 knobs refused: {msgs} — entry 66 "
              f"is two pytest calls (ASSERTED)", flush=True)
    assert alone == (("1", K_NODE_3T_PF), 2), alone
    assert all(k in msg for k, msg in msgs.items()), msgs
    assert "FEM_EM_GEO38_C21_DEGREE" in FOREIGN_ENVS_C32_EXTRA
    assert S11_DRIVE_DEGREE_ENV not in _foreign_envs()


def _k9_d2_3t(comm, m, fx, json_dir, *, report_degree=2, n_dofs=None):
    """C34's `_install` on the 3 T stubs: `_k9_stubs_3t`, `_route_drive` → ``[(2,
    mmap)]`` through ``sys.modules[_run_drive.__module__]``, `_c34_recorders`; an
    optional override of the reported ``n_dofs`` (N-Deg)."""
    import sys

    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation.test_geo38_f_clinical_box_loaded import _c34_recorders

    i36 = _modules()["i36"]
    module = sys.modules[i36._run_drive.__module__]
    _lines, z_sym, z_rc = _k9_stubs_3t(m, fx, json_dir)
    mmap = c4._build_c4_fixture(comm)["built"]["ctx"]["solver"].problem.material_map
    m.setattr(module, "_route_drive", lambda _monkeypatch: [(2, mmap)])
    env_seen, degrees = _c34_recorders(m, degree=report_degree)
    if n_dofs is not None:
        rec_fx = c4._build_c4_fixture

        def _fx(comm_):
            out = dict(rec_fx(comm_))
            out["counts"] = {"n_dofs": n_dofs}
            return out

        m.setattr(c4, "_build_c4_fixture", _fx)
    return (z_sym, z_rc), env_seen, degrees, module


@complex_only
@_coding
def test_c52_k9_3t_d2_run_drive_on_stubs(monkeypatch, tmp_path):
    """(K9-3T-d2): `_run_drive_3t(degree=2)` on the 3 T K9 stubs — C4a's variable
    '2' at the fixture call; arm R's solver degree [2]; ``i36-node-44.0615-d2.json``
    with ``"degree": 2``; (K9-Q) 4, (K9-E) 2, s2 against the 128 MHz route; the
    wrapper's (Deg) prints 7783028.  Controls: (d2-c-3T) a fixture reporting degree 1
    raises at the degree assert before any drive, solver or file; (N-Deg) a fixture
    reporting 8327394 makes the wrapper raise, before any drive."""
    import json

    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    i36 = _modules()["i36"]
    comm = MPI.COMM_WORLD
    fx = i36._k9_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    dirs = {k: os.path.join(root, k) for k in ("d2", "d2c", "ndeg")}

    with monkeypatch.context() as m:
        m.delenv(c4.C4_DEGREE_ENV, raising=False)
        (z_sym, z_rc), env_seen, degrees, module = _k9_d2_3t(comm, m, fx, dirs["d2"])
        say_lines = []
        real_say = globals()["_say"]
        m.setitem(globals(), "_say", lambda c_, msg: (say_lines.append(msg), real_say(c_, msg)))
        ret = _run_drive_3t(comm, m, K_NODE_3T_PF, degree=2)
        env_after = os.environ.get(c4.C4_DEGREE_ENV)
    drives = {"d2c": [], "ndeg": []}
    excs = {}
    for key, kw in (("d2c", {"report_degree": 1}), ("ndeg", {"n_dofs": 8_327_394})):
        with monkeypatch.context() as m:
            _l, _e, degs, module = _k9_d2_3t(comm, m, fx, dirs[key], **kw)
            real_drive = module._drive

            def _counted(*a, _key=key, _real=real_drive, **k):
                drives[_key].append(a[3])
                return _real(*a, **k)

            m.setattr(module, "_drive", _counted)
            with pytest.raises(AssertionError) as exc:
                _run_drive_3t(comm, m, K_NODE_3T_PF, degree=2)
            excs[key] = (str(exc.value), list(degs))

    rec, d1_exists, files = None, None, None
    if comm.rank == 0:
        with open(os.path.join(dirs["d2"], f"i36-node-{K_NODE_3T_PF:.4f}-d2.json")) as fh:
            rec = json.load(fh)
        d1_exists = os.path.exists(os.path.join(dirs["d2"], f"i36-node-{K_NODE_3T_PF:.4f}.json"))
        files = {k: (sorted(os.listdir(dirs[k])) if os.path.isdir(dirs[k]) else []) for k in ("d2c", "ndeg")}
    rec, d1_exists, files = comm.bcast((rec, d1_exists, files), root=0)
    r_q = {n: abs(v / 4.0 - 1.0) for n, v in rec["q"].items()}
    r_e = abs(rec["eta_t1"] / 2.0 - 1.0)
    s_m = np.array([[complex(*v) for v in row] for row in rec["M"]["s2"]])
    s_r = np.array([[complex(*v) for v in row] for row in rec["R"]["s2"]])
    r_mr = float(np.linalg.norm(s_m - s_r) / np.linalg.norm(s_m))
    s_route = h17._circuit_current_route_s(z_sym, [0, 1], z_rc, z_rc, 50.0)
    r_route = float(np.linalg.norm(s_m - s_route) / np.linalg.norm(s_route))
    deg_lines = [ln for ln in say_lines if ln.startswith("(Deg) degree 2")]
    n_deg = int(deg_lines[0].split("unknowns ")[1].split()[0]) if deg_lines else None
    if comm.rank == 0:
        print(f"\n[GEO-36 C52 K9-3T-d2] C4 variable at the fixture call {env_seen}, after {env_after!r}; "
              f"solver degrees {degrees} (arm R's); JSON -d2 degree {rec['degree']!r}, degree-1 name "
              f"written {d1_exists}; {deg_lines[0] if deg_lines else '(Deg) MISSING'} (ASSERTED)")
        print(f"[GEO-36 C52 K9-3T-d2] (K9-Q) misses from 4 {r_q}; (K9-E) miss from 2 {r_e:.3e}; s2 M vs R "
              f"{r_mr:.3e}, vs the 128 MHz circuit route {r_route:.3e} (ASSERTED)")
        print(f"[GEO-36 C52 d2-c-3T] degree 1 under degree=2 raised {excs['d2c'][0].splitlines()[0]}; "
              f"drives {drives['d2c']}, solvers {excs['d2c'][1]}, files {files['d2c']} (ASSERTED)")
        print(f"[GEO-36 C52 N-Deg] n_dofs 8327394 under degree 2 raised "
              f"{excs['ndeg'][0].splitlines()[0]}; drives {drives['ndeg']}, files {files['ndeg']} "
              f"(ASSERTED)", flush=True)
    assert env_seen == ["2"] and env_after is None, (env_seen, env_after)
    assert degrees == [2], degrees
    assert rec["degree"] == 2 and not d1_exists, (rec["degree"], d1_exists)
    assert ret["q"] == rec["q"]
    assert max(r_q.values()) <= EXACT_IDENTITY_RTOL, r_q
    assert r_e <= DISCRETE_IDENTITY_RTOL, r_e
    assert r_mr <= EXACT_IDENTITY_RTOL and r_route <= EXACT_IDENTITY_RTOL, (r_mr, r_route)
    assert len(deg_lines) == 1 and n_deg == 7_783_028, deg_lines
    assert "(1, 2)" in excs["d2c"][0], excs["d2c"][0]
    assert drives["d2c"] == [] and excs["d2c"][1] == [] and files["d2c"] == [], (drives, excs, files)
    assert "8327394" in excs["ndeg"][0] and "7783028" in excs["ndeg"][0], excs["ndeg"][0]
    assert drives["ndeg"] == [] and files["ndeg"] == [], (drives, files)
