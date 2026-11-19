"""`GEO-4` step 2 (§9 item 185): the box-term padding heuristic, gated.

``MeshGenerator.two_torus_box_padding(τ)`` turns `PORT-1` decision (4)'s fit
``deficit(W) = D∞ + C·W^(−p)`` (fitted on step 3b-xi's rungs at 0.08 / 0.10 /
0.12 m) into a padding whose PEC-box term ``|C|·W^(−p)`` is ≤ τ pp.  The gate
is physical, not a restatement of the formula: the two-torus pair is re-solved
at the heuristic's paddings for τ = 6 (0.11 m, inside the fit's range) and
τ = 4 (0.14 m, **outside** it) and the *measured* ``|deficit − D∞|`` must be
≤ τ.

Assertions:

* **(A0)** the rung constants equal the probe's; the fit equals the probe's
  output bitwise and lies within half an ulp of the printed D∞ / C / p.
* **(A1)** τ = 6 → 0.11 m, τ = 4 → 0.14 m; ``B(padding) ≤ τ < B(padding − q)``;
  ``B(unrounded) = τ`` to 1e-12.
* **(P0)** premise: the 0.08 m solve sits within step 2f's pin and within
  0.05 pp of 3b-xi's −8.0324 pp (predicted band: no 0.11-image run of this
  fixture exists before this one).
* **(G)** gate: every deficit < 0; ``|deficit|`` strictly decreasing over
  0.08 → 0.11 → 0.14; ``|deficit(W_τ) − D∞| ≤ τ``.  The p = 3 fit's
  predictions are printed beside the readings: at 0.14 m that fit gives
  ``|B| = 4.37 > 4``, so 0.14 m tells the two exponents apart.
* **(C1)** negative control, asserted: the un-sized production padding 0.08 m
  has box term ≈ 9.726 pp > both τ (predicted, and re-measured by P0).
* **(C2)** negative control: the dipolar-exponent heuristic (p = 3, the
  probe's pinned least-squares C) returns 0.10 m for τ = 4, whose recorded box
  term 6.719 pp > 4 (asserted); for τ = 6 it returns 0.09 m at ≈ 8.0 pp
  (predicted, printed).

Scope: two-torus only, an effective-range fit, not a converged box-free
deficit; the coil+phantom 0.35 rule is untouched; `PORT-1`'s ✅ is untouched
and the 0.14 m reading is a fourth-rung record for it.

Cost: heavy tier, three meshes + one driven column each (3b-xi timings).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
from mpi4py import MPI

from fem_em_solver.io.mesh import (
    PORT1_PADDING_DEFICITS_PP,
    PORT1_PADDING_RUNGS_M,
    MeshGenerator,
    _free_exponent_fit,
    _padding_for_box_term,
)
from tests.complex_mode import complex_only
from tests.validation.test_port_box_padding_sweep import (
    STEP2F_DEFICIT_AT_008,
    STEP2F_DEFICIT_TOLERANCE,
)
from tests.validation.test_port_reaction_impedance import (
    MAJOR_RADIUS,
    OMEGA,
    SEPARATION,
    TORUS_TAGS,
    _solve_reaction_z,
    mutual_inductance,
)

PROBE_PATH = Path(__file__).resolve().parents[2] / "scripts" / "probes" / "port1_padding_fit.py"

# The fit as printed in `20260812T213337Z_PORT-1-dec4-fit.log` ([fit] block),
# with the half-ulp of each printed digit.
PRINTED_D_INF, HALF_ULP_D_INF = 1.6934, 5.0e-5
PRINTED_C, HALF_ULP_C = -1.478719e-01, 5.0e-8
PRINTED_P, HALF_ULP_P = 1.6574, 5.0e-5
# The p = 3 diagnostic's C as printed in the same log.
PRINTED_C_PIN, HALF_ULP_C_PIN = -3.410151e-3, 5.0e-10

QUANTUM_M = 0.01
PRODUCTION_PADDING_M = 0.08
TAUS_PP = (6.0, 4.0)
EXPECTED_PADDING_M = {6.0: 0.11, 4.0: 0.14}
# Brief values (§9 item 185), given there as "≈".  A half-ulp pin on those
# digits was tried first and failed in `20261003T034956Z_GEO-4-smoke.log`:
# measured unrounded 0.1070673 (brief ≈ 0.10706, |Δ| 7.3e-6) and |B(0.14)|
# 3.846889 pp (brief 3.8467, |Δ| 1.9e-4) — the brief's last digits are
# truncation slips, not a different fit (A0 pins the fit bitwise).  These are
# therefore pinned to 1e-4 relative, the "≈" they were stated with; the
# binding A1 assertions are the inequality and B(unrounded) = τ to 1e-12.
EXPECTED_REL_TOL = 1.0e-4
EXPECTED_UNROUNDED_M = {6.0: 0.10706, 4.0: 0.13674}
EXPECTED_B_PP = {6.0: 5.737, 4.0: 3.8467}

# (P0): 3b-xi's 0.08 m rung, and the predicted band — a third of the smaller
# of the gate margins (0.263 / 0.153 pp under the fit).
P0_RUNG_PP = -8.0324
P0_BAND_PP = 0.05

# (C2): 3b-xi's recorded 0.10 m box term under the free fit, 6.719 pp.
C2_RECORDED_BOX_TERM_AT_010_PP = 6.719

GATE_PADDINGS = (0.08, 0.11, 0.14)


def _load_probe():
    name = "geo4_port1_padding_fit_probe"
    spec = importlib.util.spec_from_file_location(name, PROBE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _root_print(msg: str) -> None:
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _fit():
    d_inf, p = _free_exponent_fit(PORT1_PADDING_RUNGS_M, PORT1_PADDING_DEFICITS_PP)
    c = (PORT1_PADDING_DEFICITS_PP[0] - d_inf) * PORT1_PADDING_RUNGS_M[0] ** p
    return d_inf, c, p


def _dipolar_fit():
    """The probe's p = 3 least-squares fit, same arithmetic as its diagnostic."""
    x = [w ** (-3.0) for w in PORT1_PADDING_RUNGS_M]
    y = PORT1_PADDING_DEFICITS_PP
    n = len(x)
    sx, sy = sum(x), sum(y)
    sxx = sum(v * v for v in x)
    sxy = sum(a * b for a, b in zip(x, y))
    det = n * sxx - sx * sx
    c_pin = (n * sxy - sx * sy) / det
    d_pin = (sy - c_pin * sx) / n
    return d_pin, c_pin


# ---------------------------------------------------------------- arithmetic


def test_a0_fit_is_the_probes_fit():
    probe = _load_probe()
    assert tuple(PORT1_PADDING_RUNGS_M) == tuple(probe.WIDTHS_M)
    assert tuple(PORT1_PADDING_DEFICITS_PP) == tuple(probe.DEFICITS_PP)

    d_probe, p_probe = probe.power_law_extrapolation(probe.WIDTHS_M, probe.DEFICITS_PP)
    c_probe = (probe.DEFICITS_PP[0] - d_probe) * probe.WIDTHS_M[0] ** p_probe
    d_inf, c, p = _fit()
    out = MeshGenerator.two_torus_box_padding(6.0)
    _root_print(
        f"[GEO-4 A0] D_inf = {d_inf!r} pp, C = {c!r}, p = {p!r} "
        f"(probe: {d_probe!r}, {c_probe!r}, {p_probe!r})"
    )
    assert d_inf == d_probe and p == p_probe and c == c_probe
    assert out["d_inf_pp"] == d_probe and out["p"] == p_probe and out["c_pp"] == c_probe
    assert abs(d_inf - PRINTED_D_INF) <= HALF_ULP_D_INF
    assert abs(c - PRINTED_C) <= HALF_ULP_C
    assert abs(p - PRINTED_P) <= HALF_ULP_P

    d_pin, c_pin = _dipolar_fit()
    _root_print(f"[GEO-4 A0] dipolar p=3 fit: D_inf = {d_pin:+.6f} pp, C = {c_pin:+.6e}")
    assert abs(c_pin - PRINTED_C_PIN) <= HALF_ULP_C_PIN


@pytest.mark.parametrize("tau", TAUS_PP)
def test_a1_padding_for_tolerance(tau):
    out = MeshGenerator.two_torus_box_padding(tau, quantum_m=QUANTUM_M)
    c, p = out["c_pp"], out["p"]

    def box(w):
        return abs(c) * w ** (-p)

    pad, unr = out["padding_m"], out["unrounded_padding_m"]
    _root_print(
        f"[GEO-4 A1] tau = {tau} pp: padding {pad} m (unrounded {unr:.6f}), "
        f"predicted |B| = {out['predicted_box_term_pp']:.6f} pp; "
        f"|B(padding - q)| = {box(pad - QUANTUM_M):.6f} pp"
    )
    assert pad == EXPECTED_PADDING_M[tau]
    u_ref = EXPECTED_UNROUNDED_M[tau]
    b_ref = EXPECTED_B_PP[tau]
    assert abs(unr - u_ref) <= EXPECTED_REL_TOL * u_ref
    assert abs(out["predicted_box_term_pp"] - b_ref) <= EXPECTED_REL_TOL * b_ref
    assert box(pad) <= tau < box(pad - QUANTUM_M)
    assert abs(box(unr) - tau) <= 1.0e-12


def test_a1_rejects_bad_tolerance():
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError):
            MeshGenerator.two_torus_box_padding(bad)
    assert _free_exponent_fit((0.08, 0.10, 0.12), (-3.0, -5.0, -4.0)) == (None, None)


def test_c1_unsized_production_padding_exceeds_both_tolerances():
    """(C1) asserted on the fit: 0.08 m has |B| ≈ 9.726 pp, 1.62×τ=6, 2.43×τ=4."""
    d_inf, c, p = _fit()
    b = abs(c) * PRODUCTION_PADDING_M ** (-p)
    _root_print(
        f"[GEO-4 C1] production padding {PRODUCTION_PADDING_M} m: |B| = {b:.4f} pp "
        f"= {b / 6.0:.3f} x tau 6, {b / 4.0:.3f} x tau 4"
    )
    assert abs(b - 9.726) <= 5.0e-4
    for tau in TAUS_PP:
        assert b > tau
        assert MeshGenerator.two_torus_box_padding(tau)["padding_m"] > PRODUCTION_PADDING_M


def test_c2_dipolar_exponent_heuristic_undersizes():
    """(C2) p = 3 heuristic: τ = 4 → 0.10 m, recorded |B| 6.719 > 4 (asserted)."""
    d_inf, c, p = _fit()
    _d_pin, c_pin = _dipolar_fit()
    pad4, unr4 = _padding_for_box_term(c_pin, 3.0, 4.0, QUANTUM_M)
    recorded_b_010 = abs(PORT1_PADDING_DEFICITS_PP[1] - d_inf)
    pad6, unr6 = _padding_for_box_term(c_pin, 3.0, 6.0, QUANTUM_M)
    pred6 = abs(c) * pad6 ** (-p)
    _root_print(
        f"[GEO-4 C2] p=3 heuristic tau 4: {unr4:.4f} -> {pad4} m, recorded |B(0.10)| "
        f"= {recorded_b_010:.4f} pp (asserted > 4)"
    )
    _root_print(
        f"[GEO-4 C2] p=3 heuristic tau 6: {unr6:.4f} -> {pad6} m, free-fit |B| "
        f"= {pred6:.4f} pp (predicted, printed only)"
    )
    assert pad4 == 0.10
    assert abs(unr4 - 0.0948) <= 5.0e-5
    assert abs(recorded_b_010 - C2_RECORDED_BOX_TERM_AT_010_PP) <= 5.0e-4
    assert recorded_b_010 > 4.0


# ---------------------------------------------------------------- solves


@pytest.fixture(scope="module")
def gate_deficits_pp():
    """Signed ``(Im Z₂₁/ωM₁₂ − 1)·100`` at 0.08 / 0.11 / 0.14 m, one column each."""
    omega_m12 = OMEGA * mutual_inductance(MAJOR_RADIUS, MAJOR_RADIUS, SEPARATION)
    out = {}
    for padding in GATE_PADDINGS:
        z = _solve_reaction_z(
            SEPARATION, (TORUS_TAGS[0],), f"GEO-4 step 2 pad {padding}", air_padding=padding
        )
        out[padding] = 100.0 * (z[1, 0].imag / omega_m12 - 1.0)
        _root_print(
            f"[GEO-4] padding {padding} m: Im Z21 = {z[1, 0].imag:+.6e} Ohm, "
            f"omega*M12 = {omega_m12:+.6e} Ohm, deficit {out[padding]:+.4f} pp"
        )
    return out


@complex_only
def test_p0_premise_reproduces_3bxi_rung(gate_deficits_pp):
    d = gate_deficits_pp[0.08]
    delta = abs(d - P0_RUNG_PP)
    _root_print(
        f"[GEO-4 P0] 0.08 m deficit {d:+.4f} pp vs 3b-xi {P0_RUNG_PP:+.4f} pp: "
        f"|delta| = {delta:.4f} pp (band {P0_BAND_PP})"
    )
    assert abs(d / 100.0 - STEP2F_DEFICIT_AT_008) < STEP2F_DEFICIT_TOLERANCE
    assert delta <= P0_BAND_PP, (
        f"P0 off: {d:+.4f} pp vs {P0_RUNG_PP:+.4f} pp — the fit's data moved; "
        "stop before the gate"
    )


@complex_only
def test_c1_measured_unsized_box_term(gate_deficits_pp):
    d_inf, _c, _p = _fit()
    b = abs(gate_deficits_pp[0.08] - d_inf)
    _root_print(f"[GEO-4 C1] measured |deficit(0.08) - D_inf| = {b:.4f} pp > tau 6 and 4")
    assert b > max(TAUS_PP)


@complex_only
def test_g_signs_and_monotone(gate_deficits_pp):
    vals = [gate_deficits_pp[w] for w in GATE_PADDINGS]
    _root_print(
        "[GEO-4 G] deficits "
        + ", ".join(f"{w} m: {v:+.4f} pp" for w, v in zip(GATE_PADDINGS, vals))
    )
    assert all(v < 0.0 for v in vals)
    mags = [abs(v) for v in vals]
    assert all(a > b for a, b in zip(mags, mags[1:]))


@complex_only
@pytest.mark.parametrize("tau", TAUS_PP)
def test_g_measured_box_term_within_tolerance(gate_deficits_pp, tau):
    d_inf, c, p = _fit()
    out = MeshGenerator.two_torus_box_padding(tau, quantum_m=QUANTUM_M)
    w = out["padding_m"]
    assert w in gate_deficits_pp
    measured = gate_deficits_pp[w]
    excess = abs(measured - d_inf)
    d_pin, c_pin = _dipolar_fit()
    pred_free = d_inf + c * w ** (-p)
    pred_p3 = d_pin + c_pin * w ** (-3.0)
    _root_print(
        f"[GEO-4 G] tau {tau} pp at {w} m: deficit {measured:+.4f} pp, "
        f"|deficit - D_inf| = {excess:.4f} pp, margin {tau - excess:+.4f} pp "
        f"(fit margin {tau - out['predicted_box_term_pp']:+.4f}); free-fit prediction "
        f"{pred_free:+.4f} pp, p=3 prediction {pred_p3:+.4f} pp "
        f"(|p=3 prediction - D_inf| = {abs(pred_p3 - d_inf):.4f} pp)"
    )
    assert excess <= tau, f"tau = {tau}: measured box term {excess:.4f} pp at {w} m exceeds tau"
