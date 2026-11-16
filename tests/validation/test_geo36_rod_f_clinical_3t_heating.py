"""GEO-36 step 12 / §10 I85 (§9 item 932): I32-3T's body — I32b's `_run_heat`
(`test_geo36_rod_f_clinical_heating.py`) composed with C21's 3 T patch set
(`test_geo38_f_clinical_3t_box_loaded.py`, `_c21_patches`), exactly as C32 runs
I36b's `_run_drive` under it (`test_geo36_rod_f_clinical_3t_drive.py`), stub-proved
on C32's 3 T stub with I32b's literals (no solve, no gmsh).

Brief: ``docs/planning/chunks/GEO-36.md`` "§9 item 932".

Composition, never an edit: the heating module (``xl`` entry 63's), C32's (entry
61's), I36b's, C21's, C20's and ``src`` are byte-identical.  The thermal half
needs no 3 T override: `_run_heat` takes ω from `_drive_env` (C4's fixture reads
the ``wf7`` / ``c4`` globals that C21's patch 3 rebinds) and the rod centre from
``i35._rod_centre()`` (reads ``s2._f_clinical_kwargs`` at call time; patch 1).

Selected by ``FEM_EM_GEO36_S12_HEAT`` only: unset ⇒ the gated test skips and the
coding tests run; ``1`` ⇒ `_run_heat_3t` (node knob ``FEM_EM_GEO36_S8_DRIVE_C_PF``
required, read through ``i36._node_knob``); anything else ⇒ ``ValueError``.
Degree 1 only: the S8 drive / degree / heat knobs and C32's S11 knob are in the
foreign list (19 names).

Disclosed: `_columns_e_sup` (used by K12-11) hard-codes the 64 MHz node for its own
s2; K12-11 reads only its fields (not edited).

Print-only: entry 55's 1.5 T tip θ ratio rod / gel twin 1.088964 / 1.089826
(``20261106T080008Z_GEO-36-step10-I32.log:12521``).

Scope: the 3 T body and its stub proof.  Not claimed: any 128 MHz θ or tip figure
(``xl`` entry 67), entry 61's EM record, a band, a published comparand.
"""

from __future__ import annotations

import os
import re
import sys

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

S12_HEAT_ENV = "FEM_EM_GEO36_S12_HEAT"
#: Entry 55's 1.5 T tip θ(t*) ratio rod / gel twin (printed, never asserted).
ENTRY55_TIP_RATIO = ("1.088964", "1.089826")
#: The foreign guard's message head (K1-12 matches it).
FOREIGN_MSG = "I85 (GEO-36 step 12) never reads these variables"
#: I32b's stub literals (`20261108T181211Z_GEO-36-I32b-L-W0b.log:103–130`):
#: (line head, tokens).  Residual and timing tokens are excluded (brief).
T_REP = (
    ("rod share of P_vol", ("1.041667e-02", "P_vol 4.500000000e-03 W", "gel 4.453125000e-03",
                            "rod 4.687500000e-05")),
    ("(S0) census rod {3: 10260, 403: 108}", ("merged tag 3 10368",)),
    ("[gel twin] sense (1,+1j): P_gel", ("P_gel 1.125000000000e-03 W", "G 8.888888888889e+02",
                                         "G*P 1.000000000000e+00", "V_sub 9.999999999999e-04",
                                         "exterior 1728 = outer 1728 + rod wall 0")),
    ("[gel twin] F2182-shaped probe a", ("theta(900 s) +2.121079410e-01",)),
    ("[gel twin] F2182-shaped probe b", ("theta(900 s) +2.121079410e-01",)),
    ("[gel twin] 1 g ball (a 6.2035 mm) at a,", ("theta(t*) 6.062014386e-01",)),
    ("[gel twin] 1 g ball (a 6.2035 mm) at b,", ("theta(t*) 6.062135760e-01",)),
    ("[gel twin] tip ball a,", ("R_Pi dt1/10 8.0980e-07",)),
    ("[gel twin] tip ball b,", ("R_Pi dt1/10 8.0969e-07",)),
    ("[gel twin] probe interior", ("theta(t*) +6.070276390e-01",)),
    ("[gel twin] probe face", ("theta(t*) +5.217085743e-01",)),
    ("[rod] sense (1,+1j): P_gel", ("G*P 3.958333333333e+00", "V_sub 9.895833333333e-04",
                                    "rod wall 144")),
    ("[rod] F2182-shaped probe a", ("theta(900 s) +8.484317694e-01",)),
    ("[rod] F2182-shaped probe b", ("theta(900 s) +8.484317694e-01",)),
    ("[rod] 1 g ball (a 6.2035 mm) at a,", ("theta(t*) 2.424561997e+00",)),
    ("[rod] 1 g ball (a 6.2035 mm) at b,", ("theta(t*) 2.424610520e+00",)),
    ("[rod] tip ball a,", ("R_Pi dt1/10 8.0963e-07",)),
    ("[rod] tip ball b,", ("R_Pi dt1/10 8.1058e-07",)),
    ("[rod] probe interior", ("theta(t*) +2.427815119e+00",)),
    ("[rod] probe face", ("theta(t*) +2.086793306e+00",)),
    ("tip theta(t*) ratio rod / gel twin", ("a 3.999598, b 3.999598",)),
    ("gates:", ("failures none",)),
)
#: C32's stub (M2) / (Z) tokens at 44.0615 pF / 128 MHz (`20261105T133319Z_GEO-36-C32-W0.log:81`).
Z_3T = ("Z_RC 2.000000000e-01-2.821960199e+01j", "Z_th 1.215760095e+01-2.127285736e+01j",
        "V_s 1.012556726e+00+1.771729051e+00j")
#: I32b's K11 arm r R_Π dt1/10 (`L-W0b.log:159–160`) and K11-c (`:162`), by string.
K11_R_PI_R_STR = {"a": "6.855244e-07", "interior": "8.521914e-07"}
K11C_STR = "1.5793e-03"
#: K10-R's six heads (`test_geo36_rod_f_clinical_heating.py:900`).
K10R_HEADS = ("(M2)", "(S0)", "rod share", "tip theta(t*) ratio", "price:", "gates: failures none")


def _modules():
    import tests.validation.test_geo36_rod_f_clinical as i35
    import tests.validation.test_geo36_rod_f_clinical_3t_drive as c32
    import tests.validation.test_geo36_rod_f_clinical_drive as i36
    import tests.validation.test_geo36_rod_f_clinical_heating as heat
    import tests.validation.test_geo38_f_clinical_3t as t3
    import tests.validation.test_geo38_f_clinical_3t_box_loaded as c21
    import tests.validation.test_geo38_f_clinical_3t_rod_box as c20
    import tests.validation.test_geo38_f_clinical_esr_root as c4
    import tests.validation.test_wf7_f_human_b1_quadrature as wf7

    # Pytest's double import: the heating module that `_run_heat` resolves names in
    # (through `_s0_failures`, a name no stub here patches: `_heat_columns` is wrapped).
    heat = sys.modules[heat._s0_failures.__module__]
    return {"i35": i35, "c32": c32, "i36": i36, "heat": heat, "t3": t3, "c21": c21, "c20": c20,
            "c4": c4, "wf7": wf7}


def _foreign_envs():
    m = _modules()
    return tuple(m["c32"]._foreign_envs()) + (m["heat"].HEAT_ENV, m["c32"].S11_DRIVE_ENV)


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-36 I85] {msg}", flush=True)


# ---------------------------------------------------------------------------
# Knobs.
# ---------------------------------------------------------------------------
def _guard_foreign():
    bad = [k for k in _foreign_envs() if os.environ.get(k) is not None]
    if bad:
        raise ValueError(f"{FOREIGN_MSG}; unset them: {bad}")


def _s12_mode():
    raw = os.environ.get(S12_HEAT_ENV, "").strip()
    if raw == "":
        return None
    if raw == "1":
        return raw
    raise ValueError(f"{S12_HEAT_ENV} must be unset or '1', got {raw!r}")


def _knobs():
    _guard_foreign()
    mode = _s12_mode()
    return mode, _modules()["i36"]._node_knob(mode)


# ---------------------------------------------------------------------------
# The body.
# ---------------------------------------------------------------------------
def _run_heat_3t(comm, monkeypatch, c_pf):
    """C21's 3 T patch set, then I32b's `_run_heat` at degree 1 (every gate of
    I32b's asserted inside)."""
    m = _modules()
    t3, c4, wf7, i35 = m["t3"], m["c4"], m["wf7"], m["i35"]
    m["c21"]._c21_patches(monkeypatch)
    _say(comm, f"3 T patch set (C21's): rung Ø {t3.RUNG_DIAMETER_3T_M} m, f {wf7.FREQUENCY_HZ:.6e} Hz "
               f"(C4 {c4.FREQUENCY_HZ_C4:.6e}); records cells {i35.F_CLINICAL_ROD_CELL_RECORD}, U "
               f"{i35.F_CLINICAL_ROD_UNKNOWNS_RECORD}; node {c_pf:.4f} pF, degree 1; patch 4 "
               f"(box._box_readouts) is inert on this route (I32 never calls it)")
    return m["heat"]._run_heat(comm, monkeypatch, c_pf)


@complex_only
def test_geo36_i32_3t_rod_heating(monkeypatch):
    mode, c_pf = _knobs()  # node knob missing under 1 => ValueError, before gmsh
    if mode is None:
        pytest.skip(f"{S12_HEAT_ENV} unset: GEO-36 I85 is selected by environment only")
    _run_heat_3t(MPI.COMM_WORLD, monkeypatch, c_pf)


# ---------------------------------------------------------------------------
# Coding tests (no gmsh); skip whenever S12_HEAT_ENV is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(os.environ.get(S12_HEAT_ENV, "").strip() != "",
                             reason=f"{S12_HEAT_ENV} set: coding tests skip in windows")


def _literal(comm):
    return comm.size == _modules()["heat"].OPS83_RECORD_RANKS


def _lit_word(comm):
    rr = _modules()["heat"].OPS83_RECORD_RANKS
    return "ASSERTED" if _literal(comm) else f"printed at -n {comm.size} (record is -n {rr})"


@_coding
def test_i85_k1_12_knob_parsing(monkeypatch):
    """(K1-12): unset / 1 accepted, 2 / build / x refused; each of the 19 foreign
    names refused by this module's guard; the node knob required under 1."""
    mods = _modules()
    i36, c32, heat = mods["i36"], mods["c32"], mods["heat"]
    foreign = _foreign_envs()
    for k in (S12_HEAT_ENV, i36.S8_DRIVE_C_PF_ENV) + foreign:
        monkeypatch.delenv(k, raising=False)
    assert _knobs() == (None, None)
    monkeypatch.setenv(S12_HEAT_ENV, "1")
    with pytest.raises(ValueError, match="required"):
        _knobs()
    monkeypatch.setenv(i36.S8_DRIVE_C_PF_ENV, f"{c32.K_NODE_3T_PF}")
    assert _knobs() == ("1", c32.K_NODE_3T_PF)
    for bad in ("2", "build", "x"):
        monkeypatch.setenv(S12_HEAT_ENV, bad)
        with pytest.raises(ValueError, match=S12_HEAT_ENV):
            _s12_mode()
    monkeypatch.setenv(S12_HEAT_ENV, "1")
    for k in foreign:
        with monkeypatch.context() as m:
            m.setenv(k, "1")
            with pytest.raises(ValueError, match=re.escape(FOREIGN_MSG)) as exc:
                _knobs()
            assert k in str(exc.value), (k, str(exc.value))
    assert len(foreign) == len(set(foreign)) == 20, foreign
    assert "FEM_EM_GEO38_C21_DEGREE" in foreign, foreign
    assert {heat.HEAT_ENV, c32.S11_DRIVE_ENV, i36.S8_DRIVE_ENV, i36.S8_DRIVE_DEGREE_ENV} <= set(foreign)
    assert i36.S8_DRIVE_C_PF_ENV not in foreign and S12_HEAT_ENV not in foreign
    print(f"\n[GEO-36 I85 K1-12] mode unset/1 accepted, 2/build/x refused; node knob required under 1, "
          f"read {c32.K_NODE_3T_PF}; {len(foreign)} foreign variables refused by the guard ({foreign}) "
          f"(ASSERTED)", flush=True)


@_coding
def test_i85_k2_12_one_without_node_stops_before_gmsh(monkeypatch):
    """(K2-12): ``1`` without the node knob raises ValueError naming it, never gmsh."""
    from fem_em_solver.io.mesh import MeshGenerator

    i36 = _modules()["i36"]

    def _sentinel(*a, **k):
        raise AssertionError("gmsh reached")

    with monkeypatch.context() as m:
        for k in _foreign_envs() + (i36.S8_DRIVE_C_PF_ENV,):
            m.delenv(k, raising=False)
        m.setenv(S12_HEAT_ENV, "1")
        m.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(_sentinel))
        with pytest.raises(ValueError, match=i36.S8_DRIVE_C_PF_ENV):
            test_geo36_i32_3t_rod_heating(m)
    print("\n[GEO-36 I85 K2-12] S12=1 without the node knob raised ValueError naming "
          f"{i36.S8_DRIVE_C_PF_ENV} before gmsh (ASSERTED)", flush=True)


@_coding
def test_i85_k7_12_one_with_node_reaches_builder_with_3t_kwargs(monkeypatch):
    """(K7-12): ``1`` + node 44.0615 runs the gated test into C21's
    `_builder_sentinel` with C32's 3 T values; inside the patched context the rod
    centre is C20's.  (N-K7): `heat._run_heat` without C21's patches reaches the
    builder at 0.710 m / 64 MHz / 1 217 025 cells (`…C32-W0.log:68`)."""
    mods = _modules()
    t3, c20, c21, c32, i35, i36, heat = (mods[k] for k in ("t3", "c20", "c21", "c32", "i35", "i36",
                                                            "heat"))
    kw3 = c20._kwargs_3t()
    comm = MPI.COMM_WORLD
    centre_unpatched = i35._rod_centre()

    seen = []
    with monkeypatch.context() as m:
        for k in _foreign_envs():
            m.delenv(k, raising=False)
        m.setenv(S12_HEAT_ENV, "1")
        m.setenv(i36.S8_DRIVE_C_PF_ENV, f"{c32.K_NODE_3T_PF}")
        c21._builder_sentinel(m, seen)
        with pytest.raises(c21._BuilderReached):
            test_geo36_i32_3t_rod_heating(m)
        centre_patched = i35._rod_centre()
        centre_3t = c20._rod_centre_3t()

    neg = []
    with monkeypatch.context() as m:
        for k in _foreign_envs():
            m.delenv(k, raising=False)
        c21._builder_sentinel(m, neg)
        with pytest.raises(c21._BuilderReached):
            heat._run_heat(comm, m, c32.K_NODE_3T_PF)

    assert len(seen) == 1 and len(neg) == 1, (len(seen), len(neg))
    got, got_neg = seen[0], neg[0]
    kw, kwn = got["kw"], got_neg["kw"]
    rung_d, rung_d_neg = 2.0 * kw["ring_radius"], 2.0 * kwn["ring_radius"]
    rod = kw.get("implant_rod")
    print(f"\n[GEO-36 I85 K7-12] S12=1 node {c32.K_NODE_3T_PF:.4f}: builder reached with rung Ø "
          f"{rung_d:.12f} m (3 T {t3.RUNG_DIAMETER_3T_M}); rod {rod}; f wf7 {got['f_wf7']:.6e} / C4 "
          f"{got['f_c4']:.6e} Hz; records cells {got['cells']}, U {got['unknowns']}; rod centre in the "
          f"patched context {centre_patched} vs C20's {centre_3t} (ASSERTED), unpatched {centre_unpatched} "
          f"(printed; the same ASTM box)")
    print(f"[GEO-36 I85 K7-12] (N-K7) heat._run_heat without C21's patches: rung Ø {rung_d_neg:.12f} m, "
          f"f wf7 {got_neg['f_wf7']:.6e} / C4 {got_neg['f_c4']:.6e} Hz, records cells {got_neg['cells']} "
          f"(ASSERTED: 0.710 / 64 MHz / 1 217 025)", flush=True)
    assert abs(rung_d / t3.RUNG_DIAMETER_3T_M - 1.0) <= 1e-15 and t3.RUNG_DIAMETER_3T_M == 0.704, rung_d
    assert rod is not None and rod == kw3["implant_rod"] and rod["keep_cells"] is True, rod
    assert got["f_wf7"] == t3.FREQUENCY_3T_HZ == 128e6 and got["f_c4"] == 128e6, got
    assert got["cells"] == c20.F_CLINICAL_3T_ROD_CELL_RECORD == 1_128_430
    assert got["unknowns"] == c20.F_CLINICAL_3T_ROD_UNKNOWNS_RECORD == (1_508_256, 7_783_028)
    assert tuple(centre_patched) == tuple(centre_3t), (centre_patched, centre_3t)
    # (N-K7): categorical separation.
    assert abs(rung_d_neg - 0.710) <= 1e-12, rung_d_neg
    assert got_neg["f_wf7"] == 64.0e6 and got_neg["f_c4"] == 64.0e6, got_neg
    assert got_neg["cells"] == 1_217_025, got_neg["cells"]
    assert i35.F_CLINICAL_ROD_CELL_RECORD == 1_217_025
    assert mods["wf7"].FREQUENCY_HZ == 64.0e6


def _stubs_3t(m, fx, json_dir, *, p_src_r2=1.0):
    """The heating module's `_stubs` with its `_k9_stubs` replaced by C32's
    `_k9_stubs_3t` (by import), and a recorder on `_heat_columns` keeping each
    call's s2 and ω.  Returns ``(lines, keep, cols)``."""
    mods = _modules()
    heat, c32 = mods["heat"], mods["c32"]
    this = heat
    keep = {}

    def _wrapper(m_, fx_, jd_):
        lines_, z_sym, z_rc = c32._k9_stubs_3t(m_, fx_, jd_, p_src_r2=p_src_r2)
        keep.update(z_sym=z_sym, z_rc=z_rc)
        return lines_, z_sym

    m.setattr(this, "_k9_stubs", _wrapper)
    lines = heat._stubs(m, fx, json_dir)
    real_cols = this._heat_columns
    cols = []

    def _recorded(env, *a, **k):
        out = real_cols(env, *a, **k)
        cols.append({"s2": np.array(out["s2"], dtype=np.complex128),
                     "f": env["omega"] / (2.0 * np.pi)})
        return out

    m.setattr(this, "_heat_columns", _recorded)
    return lines, keep, cols


def _find(lines, head, tokens):
    hits = [ln for ln in lines if ln.lstrip("\n").startswith(head)]
    return [t for t in tokens if not any(t in ln for ln in hits)] if hits else list(tokens)


@complex_only
@_coding
def test_i85_k12r_run_heat_3t_end_to_end_on_stubs(monkeypatch, tmp_path):
    """(K12-R): `_run_heat_3t` on C32's 3 T stubs at node 44.0615 — green on its
    own gates; (Z-3T), (s2-3T), (T-rep), K10-R's heads, the sense on every θ line.
    Controls: (N-64) the 64 MHz route misses s2; (N-P) P_src × (1 + 1e-3) raises
    naming ``gel twin (P) P2 9.990e-04``."""
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL

    mods = _modules()
    i36, c32, heat = mods["i36"], mods["c32"], mods["heat"]
    comm = MPI.COMM_WORLD
    fx = i36._k9_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    good_dir, neg_dir = os.path.join(root, "good"), os.path.join(root, "neg")

    with monkeypatch.context() as m:
        lines, keep, cols = _stubs_3t(m, fx, good_dir)
        out = _run_heat_3t(comm, m, c32.K_NODE_3T_PF)
        f_seen = (mods["wf7"].FREQUENCY_HZ, mods["c4"].FREQUENCY_HZ_C4)

    with monkeypatch.context() as m:
        _stubs_3t(m, fx, neg_dir, p_src_r2=1.0 + 1.0e-3)
        with pytest.raises(AssertionError) as neg:
            _run_heat_3t(comm, m, c32.K_NODE_3T_PF)
    neg_text = str(neg.value)
    np_label = "gel twin (P) P2 9.990e-04"

    z_sym, z_rc = keep["z_sym"], keep["z_rc"]
    s_rod, s_twin = cols[0]["s2"], cols[1]["s2"]
    s_route = h17._circuit_current_route_s(z_sym, [0, 1], z_rc, z_rc, 50.0)
    r_rt = float(np.linalg.norm(s_rod - s_twin) / np.linalg.norm(s_rod))
    r_route = float(np.linalg.norm(s_rod - s_route) / np.linalg.norm(s_route))
    z_rc_64 = h17._z_rc(c32.K_NODE_3T_PF * 1e-12, 64.0e6)
    s_64 = h17._circuit_current_route_s(z_sym, [0, 1], z_rc_64, z_rc_64, 50.0)
    r_64 = float(np.linalg.norm(s_rod - s_64) / np.linalg.norm(s_route))
    f_cols = tuple(c["f"] for c in cols)

    literal = _literal(comm)
    z_miss = _find(lines, "(M2) / (Z) node 44.0615 pF", Z_3T)
    t_miss = {head: miss for head, toks in T_REP if (miss := _find(lines, head, toks))}
    heads_missing = [h for h in K10R_HEADS if not any(ln.lstrip("\n").startswith(h) for ln in lines)]
    sense_free = [ln for ln in lines if "theta(" in ln and heat.SENSE not in ln]
    ratio_line = [ln for ln in lines if ln.startswith("tip theta(t*) ratio")]
    if comm.rank == 0:
        print(f"\n[GEO-36 I85 K12-R] f during the run wf7 / C4 {f_seen}, per _heat_columns call {f_cols}; "
              f"(s2-3T) rod vs twin {r_rt:.3e}, rod vs the 128 MHz circuit route {r_route:.3e} (band "
              f"{EXACT_IDENTITY_RTOL:g}, ASSERTED; C32 read 0.000e+00 on both); (N-64) vs the 64 MHz route "
              f"{r_64:.3e} (> band, ASSERTED; C32 read 2.515e-01)")
        print(f"[GEO-36 I85 K12-R] (Z-3T) tokens missing {z_miss or 'none'}; (T-rep) I32b's literals "
              f"missing {t_miss or 'none'} ({len(T_REP)} lines; {_lit_word(comm)}); K10-R heads missing "
              f"{heads_missing or 'none'}, theta lines without {heat.SENSE}: {len(sense_free)} (ASSERTED)")
        print(f"[GEO-36 I85 K12-R] 3 T stub {ratio_line[0] if ratio_line else '(no ratio line)'}; beside "
              f"1.5 T, entry 55, before I32b's gate: a {ENTRY55_TIP_RATIO[0]}, b {ENTRY55_TIP_RATIO[1]} "
              f"(printed, never asserted)")
        print(f"[GEO-36 I85 K12-R] (N-P) control P_src x (1 + 1e-3): raised {neg_text.splitlines()[0]}; "
              f"contains {np_label!r}: {np_label in neg_text} (ceiling 1e-3/1.001 = "
              f"{1e-3 / 1.001:.3e}; C32 / I36d 9.990e-04; ASSERTED)", flush=True)
    assert set(out) == {"rod", "twin", "ratios"}
    assert f_seen == (128e6, 128e6), f_seen
    assert len(cols) == 2 and f_cols == (128e6, 128e6), f_cols
    assert r_rt <= EXACT_IDENTITY_RTOL, r_rt
    assert r_route <= EXACT_IDENTITY_RTOL, r_route
    assert r_64 > EXACT_IDENTITY_RTOL, r_64
    assert not heads_missing, heads_missing
    assert not sense_free, sense_free
    assert np_label in neg_text, neg_text
    if literal:
        assert not z_miss, z_miss
        assert not t_miss, t_miss


@complex_only
@_coding
def test_i85_k12_11_projected_tip_on_composed_route(monkeypatch, tmp_path):
    """(K12-11): K11's set-up (`heating.py:946–1018`) with C32's 3 T stubs and
    C21's patches live: (K11-Π), arm r's R_Π dt1/10 by string, (K11-rep), |pq_q|,
    the K11-u identity, (E) / (S2) on both arms, (B) on u.  (N-c): the retired DG0
    tip_s2[a] on q_r > ``S2_BAND``, reading 1.5793e-03 by string."""
    from tests.validation.test_ops83_s2tip_form_on_stub import PROJECTION_ID_TOL, _checkerboard

    mods = _modules()
    i36, heat, c21 = mods["i36"], mods["heat"], mods["c21"]
    comm = MPI.COMM_WORLD
    fx = i36._k9_fixture(comm)
    msh, cell_tags = fx[0], fx[1]
    root = comm.bcast(str(tmp_path), root=0)
    this = heat
    tips = {"a": fx[4]["a"], "interior": heat.STUB_PROBES["interior"]}
    with monkeypatch.context() as m:
        c21._c21_patches(m)
        m.setattr(this, "_k9_stubs",
                  lambda m_, fx_, jd_: mods["c32"]._k9_stubs_3t(m_, fx_, jd_)[:2])
        e_sup, sigma_f, _ = heat._columns_e_sup(m, fx, root)
        f_live = (mods["wf7"].FREQUENCY_HZ, this._drive_env(None, comm)["omega"] / (2.0 * np.pi))
        merged = heat._merge_rod_into_phantom(msh, cell_tags)
        arm_u = heat._thermal_arm("K11_3T u", msh, merged, e_sup, sigma_f, None, heat.STUB_CENTRE, tips,
                                  heat.STUB_PROBES, merged=True, tip_projection=True)
        real_source = this._heat_source

        def _rough(*a):
            return _checkerboard(real_source(*a))[0]

        m.setattr(this, "_heat_source", _rough)
        arm_r = heat._thermal_arm("K11_3T r", msh, merged, e_sup, sigma_f, None, heat.STUB_CENTRE, tips,
                                  heat.STUB_PROBES, merged=True, tip_projection=True)
    S2, E, B = heat.S2_BAND, heat.E_BAND, heat.B_BAND
    fails, lit_fails = [], []
    for arm in (arm_u, arm_r):
        if not arm["E"] <= E:
            fails.append(f"{arm['label']} (E) {arm['E']:.3e}")
        if not arm["S2"] <= S2:
            fails.append(f"{arm['label']} (S2) {arm['S2']:.3e}")
    if not arm_u["B"] <= B:
        fails.append(f"K11_3T u (B) {arm_u['B']:.3e}")
    ru, rr = arm_u["tip_rows"]["dt1/10"], arm_r["tip_rows"]["dt1/10"]
    for n in tips:
        if not arm_r["tip_s2_pi"][n] <= S2:
            fails.append(f"(K11-Pi) {n} {arm_r['tip_s2_pi'][n]:.4e}")
        if f"{rr[n]['R_Pi']:.6e}" != K11_R_PI_R_STR[n]:
            lit_fails.append(f"R_Pi(r) dt1/10 {n} {rr[n]['R_Pi']:.6e} vs {K11_R_PI_R_STR[n]}")
        if not abs(rr[n]["R_Pi"] - heat.OPS83_R_PI_R[n]) <= heat.REP_TOL_4E:
            lit_fails.append(f"(K11-rep) R_Pi {n} {rr[n]['R_Pi']:.6e} vs {heat.OPS83_R_PI_R[n]:.4e}")
        if not abs(abs(rr[n]["pq_q"]) - heat.OPS83_PQ_Q_R[n]) <= heat.REP_TOL_3E:
            fails.append(f"(K11-rep) pq_q {n} {rr[n]['pq_q']:.6e} vs {heat.OPS83_PQ_Q_R[n]:.3e}")
        if not abs(ru[n]["R_Pi"] - ru[n]["R_DG0"]) <= PROJECTION_ID_TOL:
            fails.append(f"(K11-u) R_Pi - R_DG0 {n} {ru[n]['R_Pi'] - ru[n]['R_DG0']:.3e}")
    if not abs(ru["a"]["R_Pi"] - heat.OPS83_R_PI_U_TIP) <= heat.REP_TOL_4E:
        lit_fails.append(f"(K11-u) R_Pi tip {ru['a']['R_Pi']:.6e} vs {heat.OPS83_R_PI_U_TIP:.4e}")
    neg = arm_r["tip_s2"]["a"]
    if f"{neg:.4e}" != K11C_STR:
        lit_fails.append(f"(N-c) {neg:.4e} vs {K11C_STR}")
    if comm.rank == 0:
        for arm in (arm_u, arm_r):
            print(f"\n[GEO-36 I85 K12-11] [{arm['label']}] (B) {arm['B']:.3e}; (E) {arm['E']:.3e} (band "
                  f"{E:.0e}, ASSERTED); global (S2) {arm['S2']:.4e} (band {S2:.0e}, ASSERTED)", flush=True)
            for n in tips:
                b = arm["tip_rows"]["dt1/10"][n]
                print(f"[GEO-36 I85 K12-11] [{arm['label']}] ball {n}: R_Pi dt1/10 {b['R_Pi']:.6e}; R_DG0 "
                      f"dt1/10 {b['R_DG0']:.6e}; <Pq>/<q> - 1 {b['pq_q']:+.6e}", flush=True)
        rep = ", ".join(f"{n} |diff| {abs(rr[n]['R_Pi'] - heat.OPS83_R_PI_R[n]):.3e}" for n in tips)
        print(f"[GEO-36 I85 K12-11] f live wf7 / env {f_live}; (K11-rep) {rep} (I32b read 4.402e-12 / "
              f"1.436e-12; u tip {abs(ru['a']['R_Pi'] - heat.OPS83_R_PI_U_TIP):.3e}, I32b 3.954e-12); "
              f"literal strings {_lit_word(comm)}: misses {lit_fails or 'none'}", flush=True)
        print(f"[GEO-36 I85 K12-11] (N-c) control: retired DG0 tip_s2[a] on q_r {neg:.4e} = "
              f"{neg / S2:.2f} x S2_BAND (> band, ASSERTED; I32b read 1.5793e-03); failures "
              f"{fails or 'none'}", flush=True)
    assert f_live == (128e6, 128e6), f_live
    assert not fails, fails
    assert neg > S2, neg
    if _literal(comm):
        assert not lit_fails, lit_fails
