"""`ANS-4` H33a (§9 item 624) — the lumped-sheet excess against the sheet width.

Brief: ``docs/planning/chunks/ANS-4.md`` § "§9 items 624 / 625".  F-small
exact-rectangle four-port, κ off, σ = 800 S/m wire, 64 MHz, P1 driven alone,
at three sheet widths ``f · port_dx`` with ``f`` ∈ {0.125, 0.25, 0.5} (w₀/4,
w₀/2, w₀ — the daily's 2026-10-25 re-pose).  Every sheet is the exact CAD
rectangle: ``_four_port_rung(sheet_width_fraction=f)`` sets both the strip
fraction and the narrowing fraction to ``f``; **(Wd)** proves the meshed area.

**The diagnostic, ported.**  Per gated sheet ``s`` on the driven solution:
``A_s = ∫1``, ``J_s = ∫E·ĥ``, ``L_s = ∫|n̂×E|²``, ``T_s = ∫|E·ĥ|²`` (facet
quadrature pinned), ``CS_s`` from ``sheet_variance_split`` — **copied** from
``examples/ansys_benchmarks/ans8_f_human_quadrature_mode_content_64MHz/
08_f_human_quadrature_mode_content.py:444–458`` (that script imports F-human
fixtures at call time and is not a package; the arithmetic is reproduced
verbatim) — and ``X_s = P_sheet,field − ½Re(1/Z_p)|V_s|²`` with
``P_sheet,field`` from ``test_birdcage_power_identity._sheet_field_dissipation_w``
and the gap voltage ``V_s = V_src,s − I_s Z_p`` from the package's own
``sheet_terminal_current`` (the convention of ``ports/lumped.py:572``), so
**(CS)** ``CS_s/X_s = 1`` compares two independent code paths.  Mechanism
loop ported from the same script's ``:666–709``.

Anchors (asserted, ``DISCRETE_IDENTITY_RTOL``): (CS) on all four sheets per
width; (P_w) ``P_src,exact = P_vol + ΣP_sheet,field`` (the analog of the
script's ``:581``) with its conductor-blind control missing (``:582``);
(P_acc-1) line cut ``½Re(V₁Ī₁)`` = wave form ``|V_src|²/(8 z0)(1 − |S₁₁|²)``
with ``S₁₁`` from the rung's own sweep (``:634``); (R0) at ``f = 0.5`` the
keyword path's S equals ``ANS4_EXACT_STRIP_S_P1_RECORDS["64"]`` at
``LEG_D0_REPRODUCTION_BAND`` (``-n 2`` only); (Wd) sheet area = ``f ·
port_dx · LEG_GAP_LENGTH`` at ``STRIP_BAND``.  Negative control: a planted
×1.01 on V moves ``CS/X`` off 1 by ≥ ``BLIND_FACTOR × DISCRETE_IDENTITY_RTOL``.

**Printed, never asserted:** X_s/P_acc per sheet, ΣX_s/P_acc, the LSQ width
exponent, σ_max(S) per width.  Non-monotone is a record, not a failure.

Knobs: ``FEM_EM_ANS4_H33=1`` runs the solve tests; ``FEM_EM_ANS4_H33_DEGREE``
∈ {1, 2}, unset ⇒ 1, anything else raises.

Run (complex build, ``-n 2``)::

    scripts/testing/run_and_log.sh ANS-4-H33a-W1 "docker compose exec -T fem-em-solver \\
      bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && \\
      PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_ANS4_H33=1 \\
      timeout -k 30 590 mpiexec -n 2 python3 -m pytest tests/environment \\
      tests/validation/test_ans4_h33_sheet_width_excess.py -v -s --tb=short'"
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

H33_ENV = "FEM_EM_ANS4_H33"
DEGREE_ENV = "FEM_EM_ANS4_H33_DEGREE"
# Item 625 (H33b): one width per window at degree 2.  Unset => all three.
WIDTH_ENV = "FEM_EM_ANS4_H33_WIDTH"
WIDTH_FRACTIONS = (0.125, 0.25, 0.5)
RECORD_FRACTION = 0.5
FREQUENCY_HZ_64 = 64.0e6
FREQ_KEY = "64"
DRIVEN = "P1"
PORT_IDS = ("P1", "P2", "P3", "P4")
# Planted negative control: ×1.01 on V before forming X_s.
PLANTED_V_FACTOR = 1.01


def _enabled() -> bool:
    raw = os.environ.get(H33_ENV, "").strip()
    if raw in ("", "0"):
        return False
    if raw == "1":
        return True
    raise ValueError(f"{H33_ENV}={raw!r}: only unset, 0 or 1")


def _degree() -> int:
    raw = os.environ.get(DEGREE_ENV, "").strip()
    if raw == "":
        return 1
    if raw in ("1", "2"):
        return int(raw)
    raise ValueError(f"{DEGREE_ENV}={raw!r}: only unset, 1 or 2")


def _widths() -> tuple:
    raw = os.environ.get(WIDTH_ENV, "").strip()
    if raw == "":
        return WIDTH_FRACTIONS
    for f in WIDTH_FRACTIONS:
        if raw == repr(f):
            return (f,)
    raise ValueError(f"{WIDTH_ENV}={raw!r}: only unset or one of {WIDTH_FRACTIONS}")


def sheet_variance_split(area, j_drive, l_tan, t_drive, y_s):
    """Copied verbatim from the ans8 script ``:444–458`` (see module docstring).

    Returns ``(cs, transverse, variance)`` with ``cs = ½ Re(Y_s)(L − |J|²/A)``.
    """
    g = 0.5 * float(np.real(y_s))
    a = float(np.real(area))
    mean_sq = abs(complex(j_drive)) ** 2 / a
    transverse = g * (float(np.real(l_tan)) - float(np.real(t_drive)))
    variance = g * (float(np.real(t_drive)) - mean_sq)
    return g * (float(np.real(l_tan)) - mean_sq), transverse, variance


def excess_w(p_field, v_gap, z_p):
    """``X_s = P_sheet,field − ½ Re(1/Z_p)|V_s|²``."""
    return float(p_field) - 0.5 * float(np.real(1.0 / complex(z_p))) * abs(complex(v_gap)) ** 2


def gap_voltage(v_src, current, z_p):
    """``V = V_src − I·Z_p`` — the package's convention (``ports/lumped.py:572``)."""
    return complex(v_src) - complex(current) * complex(z_p)


def p_acc_wave(v_src, s11, z0):
    """Available power ``|V_src|²/(8 Re z0)`` times ``1 − |S₁₁|²``."""
    return abs(complex(v_src)) ** 2 / (8.0 * float(np.real(z0))) * (1.0 - abs(complex(s11)) ** 2)


def width_exponent(fractions, values):
    """LSQ slope of ``log|value|`` on ``log f``; NaN unless all values share one sign."""
    v = np.asarray(values, dtype=float)
    if not (np.all(v > 0) or np.all(v < 0)):
        return float("nan")
    slope, _ = np.polyfit(np.log(np.asarray(fractions, dtype=float)), np.log(np.abs(v)), 1)
    return float(slope)


def is_monotone(values) -> bool:
    d = np.diff(np.asarray(values, dtype=float))
    return bool(np.all(d > 0) or np.all(d < 0))


# --------------------------------------------------------------------- W0 tests


def test_knobs_parse(monkeypatch):
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    assert _degree() == 1
    monkeypatch.setenv(DEGREE_ENV, "2")
    assert _degree() == 2
    monkeypatch.setenv(DEGREE_ENV, "3")
    with pytest.raises(ValueError):
        _degree()
    monkeypatch.setenv(H33_ENV, "yes")
    with pytest.raises(ValueError):
        _enabled()
    monkeypatch.delenv(H33_ENV, raising=False)
    assert _enabled() is False
    monkeypatch.delenv(WIDTH_ENV, raising=False)
    assert _widths() == WIDTH_FRACTIONS
    monkeypatch.setenv(WIDTH_ENV, "0.25")
    assert _widths() == (0.25,)
    monkeypatch.setenv(WIDTH_ENV, "0.3")
    with pytest.raises(ValueError):
        _widths()


def test_cs_equals_x_on_a_synthetic_sheet():
    """On a sheet with ``R = Z w/h`` and ``V = −J/w`` the split is the excess, exactly."""
    z_p, w, h = 50.0, 3.5e-3, 8.0e-3
    r = z_p * w / h
    y = 1.0 / r
    area = w * h
    j = (0.3 - 0.7j) * area / h * 1.3
    l_tan = 1.7 * abs(j) ** 2 / area
    t_drive = 1.2 * abs(j) ** 2 / area
    cs, tr, var = sheet_variance_split(area, j, l_tan, t_drive, y)
    assert abs(tr + var - cs) <= 1e-15 * abs(cs)
    p_field = 0.5 * y * l_tan
    v_gap = -j / w
    x = excess_w(p_field, v_gap, z_p)
    assert abs(cs / x - 1.0) <= 1e-12
    # the planted ×1.01 on V moves it far off
    x_bad = excess_w(p_field, PLANTED_V_FACTOR * v_gap, z_p)
    assert abs(cs / x_bad - 1.0) >= 1e-3


def test_line_cut_equals_wave_form_on_a_lumped_load():
    z0, v_src, z_in = 50.0, 1.0 + 0.0j, 23.0 - 41.0j
    i = v_src / (z0 + z_in)
    v = gap_voltage(v_src, i, z0)
    line = 0.5 * float(np.real(v * np.conj(i)))
    s11 = (z_in - z0) / (z_in + z0)
    assert abs(line / p_acc_wave(v_src, s11, z0) - 1.0) <= 1e-12


def test_width_exponent_and_monotone():
    f = np.array(WIDTH_FRACTIONS)
    assert abs(width_exponent(f, 3.0 * f**-1.5) + 1.5) <= 1e-12
    assert np.isnan(width_exponent(f, [1.0, -1.0, 2.0]))
    assert is_monotone([1, 2, 3]) and not is_monotone([1, 3, 2])


def test_keyword_rejects_bad_values():
    from tests.validation.test_port_birdcage_leg_offset_sweep import _four_port_rung

    with pytest.raises(ValueError):
        _four_port_rung("bad", np.zeros(4), sheet_width_fraction=1.0)
    with pytest.raises(ValueError):
        _four_port_rung(
            "bad", np.zeros(4), sheet_width_fraction=0.25, leg_sheet_strip_fraction=0.25
        )


# ------------------------------------------------------------------ W1 (solve)

needs_h33 = pytest.mark.skipif(not _enabled(), reason=f"{H33_ENV} unset")


def _width_case(f, degree, comm):
    import ufl
    from dolfinx import default_scalar_type, fem

    from tests.validation.test_birdcage_b1_plus_map import _solve_driven
    from tests.validation.test_birdcage_power_identity import (
        _exact_shares,
        _reduce,
    )
    from tests.validation.test_port_birdcage_four_port import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_leg_offset_sweep import _four_port_rung

    t0 = time.perf_counter()
    rung = _four_port_rung(
        f"ANS-4 H33 f={f}",
        np.zeros(4),
        frequency_hz=FREQUENCY_HZ_64,
        degree=degree,
        sheet_width_fraction=f,
    )
    sweep = {
        "mesh": rung["mesh"],
        "cell_tags": rung["cell_tags"],
        "problem": rung["problem"],
        "facet_tags": rung["facet_tags"],
        "specs": rung["specs"],
        "degree": int(degree),
    }
    solved = _solve_driven(sweep, DRIVEN)
    exact = _exact_shares(sweep, solved)
    msh, tags_f = sweep["mesh"], sweep["facet_tags"]
    e = solved["fields"].e_complex
    e_r = e("+")
    n_hat = ufl.FacetNormal(msh)("+")
    one = fem.Constant(msh, default_scalar_type(1.0))
    qd = {"quadrature_degree": 2 * int(degree)}
    specs = {sp.port_id: sp for sp in sweep["specs"]}
    per = {}
    for pid in PORT_IDS:
        sh = exact["sheets"][pid]
        assert sh.interior, pid
        ds_s = ufl.Measure(
            "dS",
            domain=msh,
            subdomain_data=tags_f,
            subdomain_id=(int(sh.facet_tag),),
            metadata=qd,
        )
        h_hat = ufl.as_vector([default_scalar_type(c) for c in sh.unit_drive()])
        e_h = ufl.dot(e_r, h_hat)
        nxe = ufl.cross(n_hat, e_r)
        a_s = _reduce(comm, one * ds_s)
        j_s = _reduce(comm, e_h * ds_s)
        l_s = _reduce(comm, ufl.inner(nxe, nxe) * ds_s)
        t_s = _reduce(comm, ufl.inner(e_h, e_h) * ds_s)
        y_s = 1.0 / complex(sh.sheet_resistivity)
        cs, trans, var = sheet_variance_split(a_s, j_s, l_s, t_s, y_s)
        z_p = complex(specs[pid].port_impedance_ohm)
        v = gap_voltage(sh.source_voltage_v, solved["currents"][pid], z_p)
        p_field = float(exact["sheet_field"][pid])
        x = excess_w(p_field, v, z_p)
        x_bad = excess_w(p_field, PLANTED_V_FACTOR * v, z_p)
        per[pid] = {
            "area": float(np.real(a_s)),
            "j": complex(j_s),
            "l": float(np.real(l_s)),
            "t": float(np.real(t_s)),
            "cs": cs,
            "trans": trans,
            "var": var,
            "x": x,
            "x_bad": x_bad,
            "p_field": p_field,
            "v": v,
            "z_p": z_p,
            "v_from_j": abs(complex(j_s)) / (float(np.real(a_s)) / float(sh.gap_height_m)),
        }
    v1, i1 = per[DRIVEN]["v"], complex(solved["currents"][DRIVEN])
    v_src = complex(exact["sheets"][DRIVEN].source_voltage_v)
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    p_acc_line = 0.5 * float(np.real(v1 * np.conj(i1)))
    p_wave = p_acc_wave(v_src, rung["s"][0, 0], z0)
    z_p1 = per[DRIVEN]["z_p"]
    p_norton = 0.5 * float(np.real(v1 * np.conj(v_src / z_p1)))
    return {
        "f": f,
        "rung": rung,
        "exact": exact,
        "per": per,
        "p_acc_line": p_acc_line,
        "p_acc_wave": p_wave,
        "p_norton": p_norton,
        "v_src": v_src,
        "time_s": time.perf_counter() - t0,
    }


@pytest.fixture(scope="module")
def h33_cases():
    comm = MPI.COMM_WORLD
    d = _degree()
    out = {}
    for f in _widths():
        out[f] = _width_case(f, d, comm)
        if comm.rank == 0:
            c = out[f]
            print(
                f"\n[ANS-4 H33] degree {d} f = {f}: {c['rung']['cells']} cells, "
                f"{c['time_s']:.1f} s (rung + P1 re-solve + diagnostic)",
                flush=True,
            )
    return {"degree": d, "cases": out}


@needs_h33
def test_wd_sheet_area_is_the_cad_rectangle(h33_cases):
    from tests.mesh.test_birdcage_leg_gaps import LEG_GAP_LENGTH
    from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND

    for f, c in h33_cases["cases"].items():
        dx = c["rung"]["sheet_analytic"] / LEG_GAP_LENGTH
        target = f * dx * LEG_GAP_LENGTH
        for s in c["rung"]["sheets"]:
            rel = abs(s["area"] / target - 1.0)
            if MPI.COMM_WORLD.rank == 0:
                print(
                    f"[ANS-4 H33] (Wd) f={f} sheet {s['tag']}: area {s['area']:.12e} vs "
                    f"f*dx*g {target:.12e} rel {rel:.3e} (band {STRIP_BAND:g}); "
                    f"h {s['h']:.9e} w {s['w']:.9e}",
                    flush=True,
                )
            assert rel <= STRIP_BAND, (f, s["tag"], rel)


@needs_h33
def test_cs_identity_and_planted_control(h33_cases):
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR

    for f, c in h33_cases["cases"].items():
        for pid, r in c["per"].items():
            rel = abs(r["cs"] / r["x"] - 1.0)
            bad = abs(r["cs"] / r["x_bad"] - 1.0)
            if MPI.COMM_WORLD.rank == 0:
                print(
                    f"[ANS-4 H33] f={f} {pid}: A_s {r['area']:.9e} J_s {r['j']:.9e} "
                    f"L_s {r['l']:.9e} T_s {r['t']:.9e}; CS_s {r['cs']:.9e} W, X_s "
                    f"{r['x']:.9e} W, CS/X {r['cs'] / r['x']:.12f} (rel {rel:.3e}, band "
                    f"{DISCRETE_IDENTITY_RTOL:g}); transverse {r['trans']:.9e} + variance "
                    f"{r['var']:.9e}; |V| {abs(r['v']):.9e} (|J|/w {r['v_from_j']:.9e}); "
                    f"NEG x1.01 V miss {bad:.4e}",
                    flush=True,
                )
            assert rel <= DISCRETE_IDENTITY_RTOL, (f, pid, rel)
            assert bad >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, (f, pid, bad)


@needs_h33
def test_power_identities(h33_cases):
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR

    for f, c in h33_cases["cases"].items():
        ex = c["exact"]
        rel_w = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
        blind = abs(ex["p_src"] - ex["phantom"] - ex["sheet_field_total"]) / abs(ex["p_src"])
        rel_acc = abs(c["p_acc_line"] - c["p_acc_wave"]) / abs(c["p_acc_wave"])
        rel_n = abs(ex["p_src"] / c["p_norton"] - 1.0)
        if MPI.COMM_WORLD.rank == 0:
            print(
                f"[ANS-4 H33] f={f} (P_w) P_src {ex['p_src']:.12e} = P_vol {ex['p_vol']:.12e} "
                f"+ ΣP_sheet {ex['sheet_field_total']:.12e}: rel {rel_w:.3e}; NEG conductor-"
                f"blind miss {blind:.4e}; P_src vs Norton lumped {c['p_norton']:.12e} rel "
                f"{rel_n:.3e} (printed)",
                flush=True,
            )
            print(
                f"[ANS-4 H33] f={f} (P_acc-1) line {c['p_acc_line']:.12e} W vs wave "
                f"{c['p_acc_wave']:.12e} W: rel {rel_acc:.3e}",
                flush=True,
            )
        assert rel_w <= DISCRETE_IDENTITY_RTOL, (f, rel_w)
        assert blind >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, (f, blind)
        assert rel_acc <= DISCRETE_IDENTITY_RTOL, (f, rel_acc)


@needs_h33
def test_r0_record_width_reproduces_the_record(h33_cases):
    from tests.validation.test_ans4_exact_strip_records import (
        ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
        ANS4_EXACT_STRIP_S_P1_RECORDS,
    )
    from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

    if h33_cases["degree"] != 1:
        pytest.skip("(R0) is degree 1")
    if RECORD_FRACTION not in h33_cases["cases"]:
        pytest.skip("(R0) needs the record width in this window")
    if MPI.COMM_WORLD.size != ANS4_EXACT_STRIP_RECORD_RANK_WIDTH:
        pytest.skip(f"record asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}")
    s = h33_cases["cases"][RECORD_FRACTION]["rung"]["s"]
    rec = np.asarray(ANS4_EXACT_STRIP_S_P1_RECORDS[FREQ_KEY], dtype=np.complex128)
    rel = float(np.linalg.norm(s - rec) / np.linalg.norm(rec))
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[ANS-4 H33] (R0) f={RECORD_FRACTION} ||S - S_rec||/||S_rec|| {rel:.3e} "
            f"(band {LEG_D0_REPRODUCTION_BAND:g})",
            flush=True,
        )
    assert rel <= LEG_D0_REPRODUCTION_BAND


@needs_h33
def test_print_the_width_table(h33_cases):
    """The record: printed, no band.  Non-monotone is a measured negative."""
    d = h33_cases["degree"]
    fr = list(h33_cases["cases"])
    sums = []
    for f in fr:
        c = h33_cases["cases"][f]
        p_acc = c["p_acc_wave"]
        xs = {pid: r["x"] / p_acc for pid, r in c["per"].items()}
        sums.append(sum(xs.values()))
        if MPI.COMM_WORLD.rank == 0:
            print(
                f"[ANS-4 H33] TABLE degree {d} f={f}: P_acc {p_acc:.9e} W; "
                + " ".join(f"X_{pid}/P_acc {v:+.6e}" for pid, v in xs.items())
                + f"; ΣX_s/P_acc {sums[-1]:+.6e}; sigma_max(S) "
                f"{float(c['rung']['sigma'][0]):.9f}",
                flush=True,
            )
    if len(fr) < 2:
        if MPI.COMM_WORLD.rank == 0:
            print(
                f"[ANS-4 H33] degree {d}: single-width window ({WIDTH_ENV}); exponent "
                f"and monotonicity read across windows",
                flush=True,
            )
        return
    expo = width_exponent(fr, sums)
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[ANS-4 H33] degree {d}: ΣX_s/P_acc monotone in width: {is_monotone(sums)}; "
            f"LSQ width exponent {expo:+.4f}",
            flush=True,
        )
