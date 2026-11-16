"""Temperature-dependent electrical conductivity (`MAT-5` step 1).

``σ(T) = σ₀ [1 + α (T − T₀)]`` — the linear model used for saline and gel
phantoms over a few tens of kelvin (α of order 0.02 /K for saline; the value
is an input here, never a material claim).

Numpy evaluation carries the positivity guard: ``.sigma(T)`` raises
``ValueError`` wherever ``1 + α(T − T₀) ≤ 0``.  The UFL factor
``.factor_ufl(T)`` is the bare polynomial ``1 + α(T − T₀)`` — no
``conditional``, ``max_value`` or comparison enters a form (`MAT-5` trap).
Nothing here is wired into ``materials/phantom.py`` (`MAT-1`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LinearTemperatureConductivity:
    """``σ(T) = σ_ref [1 + α (T − T_ref)]``.

    Parameters
    ----------
    sigma_ref : float
        Conductivity at ``t_ref`` [S/m], must be positive.
    alpha_per_k : float
        Linear temperature coefficient α [1/K].
    t_ref : float
        Reference temperature T₀ (same units as the ``T`` passed in).
    """

    sigma_ref: float
    alpha_per_k: float
    t_ref: float

    def __post_init__(self):
        if not self.sigma_ref > 0.0:
            raise ValueError(f"sigma_ref must be positive, got {self.sigma_ref!r}")

    def factor(self, T):
        """Numpy ``1 + α(T − T₀)`` (no guard)."""
        T = np.asarray(T, dtype=float)
        return 1.0 + self.alpha_per_k * (T - self.t_ref)

    def sigma(self, T):
        """Vectorised σ(T) [S/m]; raises where ``1 + α(T − T₀) ≤ 0``."""
        f = self.factor(T)
        if np.any(f <= 0.0):
            raise ValueError(
                "1 + alpha*(T - t_ref) <= 0: linear conductivity model is "
                f"non-positive at T = {np.asarray(T)[np.asarray(f <= 0.0)]!r}"
            )
        return self.sigma_ref * f

    def factor_ufl(self, T):
        """UFL expression ``1 + α(T − T₀)`` (bare polynomial, no guard)."""
        return 1.0 + self.alpha_per_k * (T - self.t_ref)

    def sigma_ufl(self, T):
        """UFL expression ``σ₀[1 + α(T − T₀)]``."""
        return self.sigma_ref * self.factor_ufl(T)
