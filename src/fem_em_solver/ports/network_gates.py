"""Four-port network gates as one callable (`PORT-3`, §9 item 154).

The `PORT-9` / `PORT-11` birdcage 4x4 gates — reciprocity, passivity, C4
circulance of ``Z`` and the pooled-separation control that keeps the C4 gate
from passing on noise — packaged so a stored or freshly swept ``S`` (and
optionally ``Z``) can be judged without re-deriving the expressions.  The
expressions are those of
``tests/validation/test_port_birdcage_leg_offset_sweep.py`` (the sweep's
readings) and ``tests/validation/test_port_birdcage_four_port.py`` (the gate
verdicts); the band defaults are the bands those modules pre-stated.

Stored matrices only: nothing here solves.  The `PORT-9` / `PORT-11` test
modules are not rewired onto this callable.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# Gate (i), `PORT-9` step 2c's band (``test_port_lumped_sheet_sweep``).
RECIPROCITY_BAND = 1.0e-3
# Gate (ii), pre-stated 2026-08-16 (``test_port_birdcage_four_port``).
PASSIVITY_SIGMA_TOLERANCE = 1.0e-9
# Gate (iii'), leg (c)'s adjacent-spread band (``test_port_birdcage_lumped_column``).
C4_SPREAD_BAND = 0.005
# Gate (iii')'s negative control (2): pooled / worst intra-class spread floor.
POOLED_SEPARATION_FLOOR = 10.0


def terminated_z_from_s(s, z0):
    """``Z`` from a matched-load ``S``: ``Z[:, j] = z0 (e_j + S[:, j]) / (1 - S_jj)``.

    The matched-load route's ``z_matrix``: with every other port terminated in
    ``z0``, driving port ``j`` gives ``V_i = z0 (δ_ij + S_ij) a_j`` and
    ``I_j = (1 - S_jj) a_j``, so column ``j`` of the terminated impedance is
    their ratio.  Not the open-circuit ``z0 (I + S)(I - S)^-1``; it need not be
    symmetric even when ``S`` is.
    """
    s = np.asarray(s, dtype=np.complex128)
    n = s.shape[0]
    z = np.empty_like(s)
    eye = np.eye(n, dtype=np.complex128)
    for j in range(n):
        z[:, j] = z0 * (eye[:, j] + s[:, j]) / (1.0 - s[j, j])
    return z


def class_spread(values):
    """Maximum pairwise disagreement in a class, relative to its mean magnitude.

    The class generalisation of leg (c)/(d0)'s ``|Z₂₁ − Z₄₁|/|Z₂₁|``: the widest
    complex separation any two members of the class show, divided by the class's
    own mean magnitude.  On a perfectly C4-symmetric layout every circulant class
    is constant and this is zero.
    """
    entries = np.asarray(values, dtype=np.complex128)
    scale = float(np.mean(np.abs(entries)))
    gap = 0.0
    for a in range(entries.size):
        for b in range(a + 1, entries.size):
            gap = max(gap, float(abs(entries[a] - entries[b])))
    return gap / scale if scale > 0.0 else np.inf


def circulant_classes(z_matrix):
    """The three C4 classes of a 4-port circulant ``Z``: self, adjacent, opposite.

    Indices are taken modulo ``LEG_COUNT`` on the physical layout: port ``i+1``
    and ``i-1`` are the driven leg's azimuthal neighbours (90 deg away), ``i+2``
    is the leg across the coil (180 deg).  Every off-diagonal entry of the 4x4
    lands in exactly one class, and their union is the whole off-diagonal.
    """
    n = z_matrix.shape[0]
    self_terms, adjacent, opposite = [], [], []
    for row in range(n):
        for col in range(n):
            step = (row - col) % n
            if step == 0:
                self_terms.append(z_matrix[row, col])
            elif step == n // 2:
                opposite.append(z_matrix[row, col])
            else:
                adjacent.append(z_matrix[row, col])
    return {
        "self": np.array(self_terms, dtype=np.complex128),
        "adjacent": np.array(adjacent, dtype=np.complex128),
        "opposite": np.array(opposite, dtype=np.complex128),
    }


@dataclass(frozen=True)
class FourPortNetworkGates:
    """Readings and verdicts of the four-port network gates."""

    reciprocity: float
    z_reciprocity: float
    sigma: np.ndarray
    sigma_max: float
    column_power: np.ndarray
    column_power_max: float
    spreads: dict
    worst_spread: float
    pooled: float
    separation: float
    reciprocal: bool
    passive: bool
    c4: bool
    resolved: bool
    passed: bool


def four_port_network_gates(
    s,
    *,
    z=None,
    z0_ohm=50.0,
    reciprocity_band=RECIPROCITY_BAND,
    passivity_sigma_tolerance=PASSIVITY_SIGMA_TOLERANCE,
    c4_spread_band=C4_SPREAD_BAND,
    pooled_separation_floor=POOLED_SEPARATION_FLOOR,
):
    """Judge a four-port ``S`` (and its ``Z``) against the `PORT-9` gates.

    ``z=None`` derives ``Z`` from ``S`` by :func:`terminated_z_from_s` at
    ``z0_ohm``.  Gate (i): ``‖S − Sᵀ‖/‖S‖ ≤ reciprocity_band``.  Gate (ii):
    ``σ_max(S)`` and every column power sum ``≤ 1 + passivity_sigma_tolerance``.
    Gate (iii'): every circulant class of ``Z`` spreads ``≤ c4_spread_band``.
    Control (2) (``resolved``): pooled off-diagonal spread over the worst
    intra-class spread ``≥ pooled_separation_floor``.
    """
    s_matrix = np.asarray(s, dtype=np.complex128)
    if z is None:
        z_matrix = terminated_z_from_s(s_matrix, z0_ohm)
    else:
        z_matrix = np.asarray(z, dtype=np.complex128)

    classes = circulant_classes(z_matrix)
    spreads = {n: class_spread(v) for n, v in classes.items()}
    sigma = np.linalg.svd(s_matrix, compute_uv=False)
    column_power = np.sum(np.abs(s_matrix) ** 2, axis=0)
    reciprocity = float(
        np.linalg.norm(s_matrix - s_matrix.T) / np.linalg.norm(s_matrix)
    )
    z_reciprocity = float(
        np.linalg.norm(z_matrix - z_matrix.T) / np.linalg.norm(z_matrix)
    )
    pooled = class_spread(
        np.concatenate([classes["adjacent"], classes["opposite"]])
    )
    worst = max(spreads.values())
    separation = pooled / worst if worst > 0.0 else np.inf
    sigma_max = float(np.max(sigma))
    power_max = float(np.max(column_power))

    reciprocal = bool(reciprocity <= reciprocity_band)
    passive = bool(
        sigma_max <= 1.0 + passivity_sigma_tolerance
        and all(v <= 1.0 + passivity_sigma_tolerance for v in column_power)
    )
    c4 = bool(all(v <= c4_spread_band for v in spreads.values()))
    resolved = bool(separation >= pooled_separation_floor)
    return FourPortNetworkGates(
        reciprocity=reciprocity,
        z_reciprocity=z_reciprocity,
        sigma=sigma,
        sigma_max=sigma_max,
        column_power=column_power,
        column_power_max=power_max,
        spreads=spreads,
        worst_spread=worst,
        pooled=pooled,
        separation=separation,
        reciprocal=reciprocal,
        passive=passive,
        c4=c4,
        resolved=resolved,
        passed=reciprocal and passive and c4 and resolved,
    )
