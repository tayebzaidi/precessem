"""`PORT-3` (§9 item 154): the four-port network gates as one callable.

``fem_em_solver.ports.network_gates.four_port_network_gates`` packages the
`PORT-9` / `PORT-11` birdcage 4x4 gates.  Stored matrices only, no solve:

- (A0) each default band equals the test constant it was lifted from;
- (A1) on both stored ``S`` (10 MHz ``LEG_D_S_MATRIX_10MHZ``, 64 MHz
  ``PORT11_S_64MHZ_RECORD`` — a printed control in its home module, input data
  here) every reading equals (``==``) the in-test expressions recomputed with
  the test modules' own ``_circulant_classes`` / ``_class_spread``;
- (A2) spreads, pooled spread, sigma_max and max column power reproduce the
  logged records to half a unit in the last printed digit (§9 item 178
  ruling; the separation is printed against its record, not asserted) and
  every verdict passes (reciprocity on stored ``S`` reads 0.0 exactly: the
  printed ``S`` is symmetric — expected, not a finding);
- (A3) the derived ``Z`` matches ``LEG_D0_Z_COLUMN`` and the logged 64 MHz
  ``Z`` to ``FREQUENCY_CONTROL_BAND``;
- (N1)-(N4) planted controls, each moving exactly one gate.
"""

from __future__ import annotations

import numpy as np
import pytest

from fem_em_solver.ports import C4_SPREAD_BAND as GATE_C4_SPREAD_BAND
from fem_em_solver.ports import PASSIVITY_SIGMA_TOLERANCE as GATE_PASSIVITY_SIGMA_TOLERANCE
from fem_em_solver.ports import POOLED_SEPARATION_FLOOR as GATE_POOLED_SEPARATION_FLOOR
from fem_em_solver.ports import RECIPROCITY_BAND as GATE_RECIPROCITY_BAND
from fem_em_solver.ports import four_port_network_gates, terminated_z_from_s
from tests.validation.test_ans4_resolution_ladder import PORT11_S_64MHZ_RECORD
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_Z_COLUMN,
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    _circulant_classes,
    _class_spread,
)
from tests.validation.test_port_birdcage_larmor_gate import (
    FREQUENCY_CONTROL_BAND,
    LEG_D_S_MATRIX_10MHZ,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

Z0_OHM = 50.0

# The 64 MHz terminated Z logged beside PORT11_S_64MHZ_RECORD,
# `20260826T110434Z_PORT-11-step2.log.gz:9259-9262` (rows i, columns k = port k
# driven).  Version-tagged: 116 085-cell record mesh, `c4_congruent_sheets`
# off, power-wave route, 0.11 image.
PORT11_Z_64MHZ_RECORD = np.array(
    [
        [
            +2.647082952e01 + 4.646185233e01j,
            +1.877248178e01 + 6.880903941e-01j,
            +1.429887332e01 - 4.749644673e00j,
            +1.877380478e01 + 6.929877997e-01j,
        ],
        [
            +1.877079735e01 + 6.864775531e-01j,
            +2.647408805e01 + 4.647247309e01j,
            +1.877301423e01 + 6.893105431e-01j,
            +1.429704807e01 - 4.748654779e00j,
        ],
        [
            +1.429428638e01 - 4.749063864e00j,
            +1.876908728e01 + 6.896008822e-01j,
            +2.649109802e01 + 4.648038833e01j,
            +1.876465515e01 + 6.882647066e-01j,
        ],
        [
            +1.877383419e01 + 6.947656906e-01j,
            +1.429928713e01 - 4.746587865e00j,
            +1.877029497e01 + 6.913645278e-01j,
            +2.647483350e01 + 4.645438954e01j,
        ],
    ],
    dtype=np.complex128,
)

# (A2) records, stored as the strings the logs printed.  64 MHz:
# `20260826T110434Z_PORT-11-step2.log.gz:9269-9271` (sigma, column power,
# spreads) and `:13920-13928` (separation).  10 MHz: `LEG_D_S_MATRIX_10MHZ`'s
# gate prints (spreads 0.0553 / 0.0353 / 0.0214 %, separation 166.6766x,
# sigma_max 0.999992805), as cited by the §9 item 154 brief.  Per the §9 item
# 178 ruling (2026-10-02 10:30 review) each asserted reading must lie within
# half a unit in the last printed digit of its record; the separation is
# printed against its record, not asserted.
RECORDS = {
    "64MHz": {
        "spreads_pct": {"self": "0.0573", "adjacent": "0.0599", "opposite": "0.0370"},
        "pooled_pct": "40.1838",
        "separation": "671.0527",
        "sigma_max": "0.999721388",
        "column_power_max": "0.804704664",
    },
    "10MHz": {
        "spreads_pct": {"self": "0.0553", "adjacent": "0.0353", "opposite": "0.0214"},
        "separation": "166.6766",
        "sigma_max": "0.999992805",
    },
}

STORED_S = {"10MHz": LEG_D_S_MATRIX_10MHZ, "64MHz": PORT11_S_64MHZ_RECORD}


def _half_ulp(text):
    """Half a unit in the last printed digit of the record string ``text``."""
    decimals = len(text.split(".", 1)[1]) if "." in text else 0
    return 0.5 * 10.0 ** (-decimals)


def _within_half_ulp(label, reading, text):
    gap = abs(reading - float(text))
    tol = _half_ulp(text)
    print(
        f"    {label}: reading {reading:.12g} vs record {text}  "
        f"gap {gap:.3e} <= half-ulp {tol:.0e}  margin {tol - gap:.3e}"
    )
    return gap <= tol


def test_a0_default_bands_equal_their_test_constants():
    pairs = [
        ("RECIPROCITY_BAND", GATE_RECIPROCITY_BAND, RECIPROCITY_BAND),
        ("PASSIVITY_SIGMA_TOLERANCE", GATE_PASSIVITY_SIGMA_TOLERANCE, PASSIVITY_SIGMA_TOLERANCE),
        ("C4_SPREAD_BAND / ADJACENT_SPREAD_BAND", GATE_C4_SPREAD_BAND, ADJACENT_SPREAD_BAND),
        ("POOLED_SEPARATION_FLOOR", GATE_POOLED_SEPARATION_FLOOR, POOLED_SEPARATION_FLOOR),
    ]
    for name, mine, theirs in pairs:
        print(f"[PORT-3 A0] {name}: callable {mine!r} vs test {theirs!r}")
        assert mine == theirs, name


@pytest.mark.parametrize("key", ["10MHz", "64MHz"])
def test_a1_readings_equal_in_test_expressions_bitwise(key):
    s_matrix = np.asarray(STORED_S[key], dtype=np.complex128)
    gates = four_port_network_gates(s_matrix, z0_ohm=Z0_OHM)

    n = s_matrix.shape[0]
    eye = np.eye(n, dtype=np.complex128)
    z_matrix = np.empty_like(s_matrix)
    for j in range(n):
        z_matrix[:, j] = Z0_OHM * (eye[:, j] + s_matrix[:, j]) / (1.0 - s_matrix[j, j])

    classes = _circulant_classes(z_matrix)
    spreads = {name: _class_spread(v) for name, v in classes.items()}
    sigma = np.linalg.svd(s_matrix, compute_uv=False)
    column_power = np.sum(np.abs(s_matrix) ** 2, axis=0)
    reciprocity = float(np.linalg.norm(s_matrix - s_matrix.T) / np.linalg.norm(s_matrix))
    z_reciprocity = float(np.linalg.norm(z_matrix - z_matrix.T) / np.linalg.norm(z_matrix))
    pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
    worst = max(spreads.values())
    separation = pooled / worst if worst > 0.0 else np.inf

    print(f"\n[PORT-3 A1] {key}: callable vs in-test (== required)")
    checks = [
        ("reciprocity", gates.reciprocity, reciprocity),
        ("z_reciprocity", gates.z_reciprocity, z_reciprocity),
        ("sigma_max", gates.sigma_max, float(np.max(sigma))),
        ("column_power_max", gates.column_power_max, float(np.max(column_power))),
        ("pooled", gates.pooled, pooled),
        ("worst_spread", gates.worst_spread, worst),
        ("separation", gates.separation, separation),
    ] + [(f"spread[{k}]", gates.spreads[k], spreads[k]) for k in spreads]
    for name, mine, theirs in checks:
        print(f"    {name}: {mine!r} vs {theirs!r}  " f"{'EQUAL' if mine == theirs else 'DIFFER'}")
    for name, mine, theirs in checks:
        assert mine == theirs, name
    assert np.array_equal(gates.sigma, sigma)
    assert np.array_equal(gates.column_power, column_power)
    assert np.array_equal(terminated_z_from_s(s_matrix, Z0_OHM), z_matrix)


@pytest.mark.parametrize("key", ["10MHz", "64MHz"])
def test_a2_readings_reproduce_the_records(key):
    gates = four_port_network_gates(STORED_S[key], z0_ohm=Z0_OHM)
    rec = RECORDS[key]
    print(f"\n[PORT-3 A2] {key}:")
    print(
        f"    reciprocity ||S - S^T||/||S|| = {gates.reciprocity!r} "
        "(stored S is symmetric: 0.0 expected)"
    )
    print(f"    z_reciprocity = {gates.z_reciprocity:.9e} (reported)")
    results = {}
    for name, text in rec["spreads_pct"].items():
        results[f"spread {name} %"] = _within_half_ulp(
            f"spread {name} %", gates.spreads[name] * 100.0, text
        )
    if "pooled_pct" in rec:
        results["pooled %"] = _within_half_ulp("pooled %", gates.pooled * 100.0, rec["pooled_pct"])
    else:
        print(f"    pooled: {gates.pooled * 100.0:.6f}% (reported)")
    results["sigma_max"] = _within_half_ulp("sigma_max", gates.sigma_max, rec["sigma_max"])
    if "column_power_max" in rec:
        results["column_power_max"] = _within_half_ulp(
            "column_power_max", gates.column_power_max, rec["column_power_max"]
        )
    else:
        print(f"    column_power_max: {gates.column_power_max:.12f} (reported)")
    sep_text = rec["separation"]
    sep_gap = abs(gates.separation - float(sep_text))
    print(
        f"    separation: {gates.separation:.7f}x vs record {sep_text}x  "
        f"deviation {sep_gap / _half_ulp(sep_text):.2f} half-ulps, "
        f"{sep_gap / float(sep_text):.3e} relative  "
        "-- printed, not asserted (item 178 ruling)"
    )
    print(
        f"    verdicts: reciprocal={gates.reciprocal} passive={gates.passive} "
        f"c4={gates.c4} resolved={gates.resolved} passed={gates.passed}"
    )
    assert gates.reciprocity == 0.0
    for name, ok in results.items():
        assert ok, name
    assert gates.reciprocal and gates.passive and gates.c4 and gates.resolved
    assert gates.passed


def _relative(derived, record):
    return np.abs(derived - record) / np.abs(record)


def test_a3_derived_z_matches_the_logged_z():
    z10 = terminated_z_from_s(LEG_D_S_MATRIX_10MHZ, Z0_OHM)
    z64 = terminated_z_from_s(PORT11_S_64MHZ_RECORD, Z0_OHM)
    dev10 = _relative(z10[:, 0], LEG_D0_Z_COLUMN)
    dev64 = _relative(z64, PORT11_Z_64MHZ_RECORD)
    print(f"\n[PORT-3 A3] band {FREQUENCY_CONTROL_BAND:.0e} relative")
    for i in range(4):
        print(
            f"    10 MHz Z_{i + 1}1: derived {z10[i, 0]:.9e} vs "
            f"{LEG_D0_Z_COLUMN[i]:.9e}  rel {dev10[i]:.3e}"
        )
    for i in range(4):
        for k in range(4):
            print(
                f"    64 MHz Z_{i + 1}{k + 1}: derived {z64[i, k]:.9e} vs "
                f"{PORT11_Z_64MHZ_RECORD[i, k]:.9e}  rel {dev64[i, k]:.3e}"
            )
    print(f"    max rel: 10 MHz {np.max(dev10):.3e}, 64 MHz {np.max(dev64):.3e}")
    assert float(np.max(dev10)) <= FREQUENCY_CONTROL_BAND
    assert float(np.max(dev64)) <= FREQUENCY_CONTROL_BAND


def _verdicts(g):
    return {"reciprocal": g.reciprocal, "passive": g.passive, "c4": g.c4, "resolved": g.resolved}


def _planted(label, s, z, expect_false):
    gates = four_port_network_gates(s, z=z, z0_ohm=Z0_OHM)
    verdicts = _verdicts(gates)
    print(
        f"\n[PORT-3 {label}] reciprocity {gates.reciprocity:.6e}  sigma_max "
        f"{gates.sigma_max:.12f}  column_power_max {gates.column_power_max:.9f}  "
        f"spreads "
        + ", ".join(f"{k} {v * 100:.4f}%" for k, v in gates.spreads.items())
        + f"  pooled {gates.pooled * 100:.4f}%  separation {gates.separation:.4f}x"
        f"\n    verdicts {verdicts}  passed={gates.passed}  "
        f"(predicted: only {expect_false} False)"
    )
    for name, value in verdicts.items():
        assert value is (name != expect_false), (label, name, value)
    assert gates.passed is False
    return gates


# Planted controls on the 64 MHz record (pooled 40%: N3's 1% self spread keeps
# the separation above its floor, so only the C4 gate can move).
def _base():
    s = np.array(PORT11_S_64MHZ_RECORD, dtype=np.complex128)
    return s, terminated_z_from_s(s, Z0_OHM)


N1_PHASE_RAD = np.deg2rad(1.0)


def test_n1_planted_nonreciprocal_phase_fails_only_reciprocity():
    """(N1), §9 item 178: S' = D S D^H, D = diag(e^{j phi}, 1, 1, 1), phi = 1 deg.

    A lossless non-reciprocal phase shifter in port 1's feed.  D is unitary, so
    sigma(S') = sigma(S) and |S'_ik| = |S_ik| (passivity and column powers
    cannot move); the reciprocity residual is 2|sin phi| K with
    K = sqrt(2 sum_{k=1..3} |S_0k|^2) / ||S||_F = 0.538108 on the 64 MHz
    record, predicted 1.8783e-2.  Built elementwise (diagonal untouched
    bitwise) and passed the UNPLANTED derived z, so the Z-class gates cannot
    move either.
    """
    s0, z = _base()
    s = s0.copy()
    for k in range(1, 4):
        s[0, k] *= np.exp(1j * N1_PHASE_RAD)
        s[k, 0] *= np.exp(-1j * N1_PHASE_RAD)
    assert np.array_equal(np.diag(s), np.diag(s0))
    gates = _planted("N1", s, z, "reciprocal")
    base = four_port_network_gates(s0, z=z, z0_ohm=Z0_OHM)
    d_sigma = gates.sigma_max - base.sigma_max
    d_col = float(np.max(np.abs(gates.column_power - base.column_power)))
    k_factor = np.sqrt(2.0 * np.sum(np.abs(s0[0, 1:]) ** 2)) / np.linalg.norm(s0)
    predicted = 2.0 * abs(np.sin(N1_PHASE_RAD)) * k_factor
    print(
        f"    N1: reciprocity {gates.reciprocity:.6e} vs predicted "
        f"2|sin phi|K = {predicted:.6e} (K = {k_factor:.6f}; predicted, "
        f"not asserted); {gates.reciprocity / GATE_RECIPROCITY_BAND:.2f}x band"
    )
    print(
        f"    N1: sigma_max(S') - sigma_max(S) = {d_sigma:.3e}; "
        f"max |delta column power| = {d_col:.3e}"
    )
    print(
        f"    N1: spreads/pooled/separation unchanged: "
        f"{gates.spreads == base.spreads and gates.pooled == base.pooled and gates.separation == base.separation}"
    )


def test_n2_planted_gain_fails_only_passivity():
    s, z = _base()
    sigma_max = float(np.max(np.linalg.svd(s, compute_uv=False)))
    s = s * ((1.0 + 1.0e-6) / sigma_max)
    _planted("N2", s, z, "passive")


def test_n3_planted_self_term_fails_only_c4():
    s, z = _base()
    z[0, 0] *= 1.01
    gates = _planted("N3", s, z, "c4")
    assert abs(gates.spreads["self"] - 0.01) < 0.002


def test_n4_planted_adjacent_equal_opposite_fails_only_resolved():
    s, z = _base()
    classes = _circulant_classes(z)
    opposite_mean = complex(np.mean(classes["opposite"]))
    n = z.shape[0]
    for row in range(n):
        for col in range(n):
            if (row - col) % n in (1, n - 1):
                z[row, col] = opposite_mean
    _planted("N4", s, z, "resolved")
