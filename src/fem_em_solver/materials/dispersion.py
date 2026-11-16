"""Frequency-dispersive permittivity: Debye and Cole–Cole ε̂(ω) (`MAT-3` step 1).

The multi-pole Cole–Cole form used for biological tissue by Gabriel, Lau &
Gabriel (1996), *Phys. Med. Biol.* 41:2271 (their eq. 3, four poles) is

    ε̂(ω) = ε∞ + Σₙ Δεₙ / (1 + (jωτₙ)^(1−αₙ)) + σ_i / (jωε₀)

with ``Δεₙ = ε_s,n − ε∞`` the pole strength, ``τₙ`` its relaxation time,
``αₙ ∈ [0, 1)`` the Cole–Cole broadening (Cole & Cole 1941, *J. Chem. Phys.*
9:341 — the depressed-arc locus) and ``σ_i`` the static ionic conductivity.
``αₙ = 0`` is the Debye pole ``Δε / (1 + jωτ)``; Jin, *The Finite Element
Method in Electromagnetics* 3rd ed., ch. 12 §12.4.1 writes the same Debye
permittivity (as an equation image, so the sign was checked here, not
copied: Jin's frequency-domain chapters use ``e^{+jωt}``, under which a
passive pole must have ``Im ε̂ < 0``).

Convention: the project's ``e^{+jωt}`` (`TH-1` formulation note), so

    ε̂ = ε′ − j ε″,   ε″ > 0 for a lossy (passive) medium,

and the effective conductivity the solver sees is ``σ_eff = ω ε₀ ε″`` with
``εᵣ = ε′`` — exactly the inverse of
:func:`fem_em_solver.utils.analytical.complex_permittivity`.

The fractional power is evaluated in polar form,
``(jωτ)^(1−α) = (ωτ)^(1−α) · e^{jπ(1−α)/2}``, never through a complex ``**``
on ``1j*ω*τ``, so the branch cannot drift.

This module is numpy-only.  ε̂(ω) is evaluated here, at the solve frequency,
and handed to the FEM as real scalars through :meth:`as_homogeneous`; no
dispersion enters UFL.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np

from ..utils.constants import EPSILON_0

Pole = Tuple[float, float, float]  # (Δε, τ [s], α)


@dataclass(frozen=True)
class ColeColeMaterial:
    """Multi-pole Cole–Cole dielectric with static ionic conductivity.

    Parameters
    ----------
    epsilon_inf:
        High-frequency relative permittivity ``ε∞``.
    poles:
        Sequence of ``(Δε, τ, α)``: pole strength, relaxation time in seconds,
        Cole–Cole broadening ``0 ≤ α < 1`` (``α = 0`` is Debye).
    sigma_ionic:
        Static ionic conductivity ``σ_i`` in S/m.
    mu_r:
        Relative permeability, carried through to :meth:`as_homogeneous`.
    """

    epsilon_inf: float
    poles: Tuple[Pole, ...] = ()
    sigma_ionic: float = 0.0
    mu_r: float = 1.0

    def __post_init__(self) -> None:
        poles = tuple(tuple(float(v) for v in p) for p in self.poles)
        for delta, tau, alpha in poles:
            if not (np.isfinite(delta) and delta >= 0.0):
                raise ValueError(f"pole strength Δε must be finite and >= 0, got {delta!r}")
            if not (np.isfinite(tau) and tau > 0.0):
                raise ValueError(f"relaxation time τ must be finite and > 0, got {tau!r}")
            if not (0.0 <= alpha < 1.0):
                raise ValueError(f"Cole–Cole α must satisfy 0 <= α < 1, got {alpha!r}")
        if not (np.isfinite(self.epsilon_inf) and self.epsilon_inf > 0.0):
            raise ValueError(f"epsilon_inf must be finite and > 0, got {self.epsilon_inf!r}")
        if not (np.isfinite(self.sigma_ionic) and self.sigma_ionic >= 0.0):
            raise ValueError(f"sigma_ionic must be finite and >= 0, got {self.sigma_ionic!r}")
        object.__setattr__(self, "poles", poles)

    @classmethod
    def debye(
        cls,
        epsilon_inf: float,
        poles: Sequence[Tuple[float, float]] = (),
        sigma_ionic: float = 0.0,
        mu_r: float = 1.0,
    ) -> "ColeColeMaterial":
        """Debye material: ``poles`` are ``(Δε, τ)`` pairs, ``α = 0`` for all."""
        return cls(
            epsilon_inf=epsilon_inf,
            poles=tuple((float(d), float(t), 0.0) for d, t in poles),
            sigma_ionic=sigma_ionic,
            mu_r=mu_r,
        )

    @property
    def epsilon_static(self) -> float:
        """``ε_s = ε∞ + Σ Δεₙ`` (the ω → 0 limit of the dielectric part)."""
        return float(self.epsilon_inf + sum(p[0] for p in self.poles))

    def complex_permittivity(self, frequency_hz):
        """ε̂(f), vectorised over ``frequency_hz`` (scalar in → complex out)."""
        f = np.asarray(frequency_hz, dtype=float)
        if np.any(f <= 0.0):
            raise ValueError("frequency_hz must be > 0")
        omega = 2.0 * np.pi * f
        eps = np.full(f.shape, complex(self.epsilon_inf), dtype=complex)
        for delta, tau, alpha in self.poles:
            power = 1.0 - alpha
            # (jωτ)^(1−α) in polar form: |·| = (ωτ)^(1−α), arg = π(1−α)/2.
            z = (omega * tau) ** power * np.exp(1j * 0.5 * np.pi * power)
            eps = eps + delta / (1.0 + z)
        # σ/(jωε₀) = −jσ/(ωε₀) under e^{+jωt}.
        eps = eps - 1j * self.sigma_ionic / (omega * EPSILON_0)
        if eps.ndim == 0:
            return complex(eps)
        return eps

    def sigma_epsilon_at(self, frequency_hz):
        """``(σ_eff, εᵣ) = (ωε₀·(−Im ε̂), Re ε̂)`` at ``frequency_hz``."""
        eps = self.complex_permittivity(frequency_hz)
        omega = 2.0 * np.pi * np.asarray(frequency_hz, dtype=float)
        sigma_eff = omega * EPSILON_0 * (-np.imag(eps))
        eps_r = np.real(eps)
        if np.ndim(sigma_eff) == 0:
            return float(sigma_eff), float(eps_r)
        return sigma_eff, eps_r

    def as_homogeneous(self, frequency_hz: float):
        """The frozen-frequency :class:`HomogeneousMaterial` the solver consumes."""
        # Imported here: core/time_harmonic.py imports ``materials`` before it
        # defines HomogeneousMaterial, so a module-level import is circular.
        from ..core.time_harmonic import HomogeneousMaterial

        sigma_eff, eps_r = self.sigma_epsilon_at(float(frequency_hz))
        return HomogeneousMaterial(sigma=sigma_eff, epsilon_r=eps_r, mu_r=self.mu_r)


__all__ = ["ColeColeMaterial"]
