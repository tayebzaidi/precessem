"""Layered (foil-on-backing) surface impedance by the transmission-line recursion (`TH-18`).

A planar stack seen from the field side is a cascade of transmission-line
sections terminated by a backing impedance ``Z_b`` (Jin §3.4.2, eqs.
(3.87)–(3.89); the surface-impedance boundary itself is Jin §1.5.3, eqs.
(1.54)–(1.56)). One layer of characteristic impedance ``Z_c``, propagation
constant ``γ`` and thickness ``t`` maps the impedance ``Z_L`` below it to

    Z_in = Z_c (Z_L + Z_c tanh γt) / (Z_c + Z_L tanh γt).

Conductor layers use the **good-conductor** form, bit for bit the one
`TH-14`'s Leontovich wall uses (``core.cavity.surface_resistance_ohm``):

    R_s = √(ωμ₀/(2σ)),  Z_c = (1 + j) R_s,  γ = (1 + j)/δ,  δ = √(2/(ωμ₀σ)).

The exact √(jωμ/(σ + jωε)) differs from this by ≈ ωε₀/σ (≈ 1.6e-6 at
σ = 1e4 S/m, 291 MHz); a thick layer must reproduce the Leontovich Z_s of
`TH-14` exactly, so the displacement term is deliberately dropped.

Dielectric layers are lossless: ``Z_c = η₀√(μ_r/ε_r)``, ``γ = jω√(μ_r ε_r)/c``.

Convention e^{+jωt}: ``X = Im Z_s > 0`` is inductive. Every function takes the
angular frequency ω [rad/s], not Hz (the ``skin_depth`` helpers elsewhere in
the package take Hz).

numpy and ``utils.constants`` only — no DolfinX import.
"""

from dataclasses import dataclass
from typing import Sequence, Tuple, Union

import numpy as np

from ..utils.constants import C_0, ETA_0, MU_0


@dataclass(frozen=True)
class ConductorLayer:
    """A good-conductor layer of conductivity ``sigma_s_per_m`` and ``thickness_m``."""

    sigma_s_per_m: float
    thickness_m: float


@dataclass(frozen=True)
class DielectricLayer:
    """A lossless dielectric layer (``eps_r``, ``mu_r``) of ``thickness_m``."""

    eps_r: float
    thickness_m: float
    mu_r: float = 1.0


Layer = Union[ConductorLayer, DielectricLayer]


def conductor_surface_resistance_ohm(omega_rad_s: float, sigma_s_per_m: float) -> float:
    """``R_s = √(ωμ₀/(2σ))`` — the same expression as ``core.cavity.surface_resistance_ohm``."""

    return float(np.sqrt(omega_rad_s * MU_0 / (2.0 * sigma_s_per_m)))


def conductor_skin_depth_m(omega_rad_s: float, sigma_s_per_m: float) -> float:
    """``δ = √(2/(ωμ₀σ))`` [m], from ω (rad/s)."""

    return float(np.sqrt(2.0 / (omega_rad_s * MU_0 * sigma_s_per_m)))


def layer_line_constants(omega_rad_s: float, layer: Layer) -> Tuple[complex, complex]:
    """``(Z_c, γ)`` of one layer at ω."""

    if isinstance(layer, ConductorLayer):
        r_s = conductor_surface_resistance_ohm(omega_rad_s, layer.sigma_s_per_m)
        delta = conductor_skin_depth_m(omega_rad_s, layer.sigma_s_per_m)
        return (1.0 + 1.0j) * r_s, (1.0 + 1.0j) / delta
    if isinstance(layer, DielectricLayer):
        z_c = ETA_0 * np.sqrt(layer.mu_r / layer.eps_r)
        gamma = 1.0j * omega_rad_s * np.sqrt(layer.mu_r * layer.eps_r) / C_0
        return complex(z_c), complex(gamma)
    raise TypeError(f"unknown layer type {type(layer).__name__}")


def transform_impedance(z_c: complex, gamma_t: complex, z_load: complex) -> complex:
    """One transmission-line section: ``Z_c (Z_L + Z_c tanh γt)/(Z_c + Z_L tanh γt)``."""

    th = np.tanh(gamma_t)
    return complex(z_c * (z_load + z_c * th) / (z_c + z_load * th))


def layered_surface_impedance(
    omega_rad_s: float,
    layers: Sequence[Layer],
    backing_impedance_ohm: complex = 0.0,
) -> complex:
    """Input surface impedance ``Z_s = R + jX`` [Ω] of a planar stack.

    ``layers`` are listed **from the field side**; the recursion starts at the
    backing (``backing_impedance_ohm``; 0 is a PEC backing, η₀ a matched
    free-space backing) and works outward through ``reversed(layers)``.
    """

    z = complex(backing_impedance_ohm)
    for layer in reversed(list(layers)):
        z_c, gamma = layer_line_constants(omega_rad_s, layer)
        z = transform_impedance(z_c, gamma * layer.thickness_m, z)
    return z
