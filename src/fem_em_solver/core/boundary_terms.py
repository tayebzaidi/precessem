"""Surface terms for the complex curl-curl form (`TH-5` step 1a).

Convention: ``e^{+jωt}``; the weak form the terms add to is

    ∫ (∇×u)·(∇×v̄) dx − k₀² ∫ u·v̄ dx + γ ∮ (n×u)·(n×v̄) ds = …

which follows from ∫(∇×∇×E)·v̄ = ∫(∇×E)·(∇×v̄) + ∮(n×∇×E)·v̄ ds and the
boundary condition n×∇×E = −γ n×(n×E). For an outgoing plane wave
E = ẑ e^{−jk₀x} on the face n = x̂ this holds with γ = +jk₀, the first-order
ABC (Jin 3e §9.2 (9.44) at s = 1, §9.3 (9.54)). A conjugated γ turns the
boundary into a source (the trap named in the `TH-5` entry).

The Leontovich impedance condition n×∇×E = −(jωμ₀/Z_s) n×(n×E) is the same
term with γ = jωμ₀/Z_s; at Z_s = η₀ it is the first-order ABC exactly.

Nothing here is exported from ``core/__init__.py`` (additive by design).
Step 3 appends a polynomial-graded PEC-backed PML slab (Jin 3e §9.6) below.
"""

from __future__ import annotations

import numpy as np
import ufl

from fem_em_solver.utils.constants import C_0, MU_0


def impedance_boundary_term(msh, facet_tags, tag, surface_impedance_ohm, omega_rad_s):
    """Return ``f(trial, test)`` = jωμ₀/Z_s ∮ (n×u)·(n×v̄) ds(tag)."""
    gamma = 1j * float(omega_rad_s) * MU_0 / complex(surface_impedance_ohm)
    return _tangential_mass_term(msh, facet_tags, tag, gamma)


def first_order_abc_term(msh, facet_tags, tag, omega_rad_s, cos_theta=1.0):
    """Return ``f(trial, test)`` = j k₀ cosθ ∮ (n×u)·(n×v̄) ds(tag); vacuum only.

    ``cos_theta = 1`` is the first-order ABC; ``cos_theta = cos θ`` is Jin's
    exact condition for a plane wave at a known angle θ (reflectionless there).
    """
    k0 = float(omega_rad_s) / C_0
    gamma = 1j * k0 * float(cos_theta)
    return _tangential_mass_term(msh, facet_tags, tag, gamma)


def _tangential_mass_term(msh, facet_tags, tag, gamma):
    n = ufl.FacetNormal(msh)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=facet_tags)

    def form(trial, test):
        return gamma * ufl.inner(ufl.cross(n, trial), ufl.cross(n, test)) * ds(tag)

    return form


def abc_first_order_reflection_te(cos_theta):
    """Plane-wave TE reflection of the first-order ABC, R = (cosθ − 1)/(cosθ + 1).

    Jin 3e (9.48) / the first-order Padé (9.4): 0 at normal incidence, −1/3 at
    60°, → −1 at grazing incidence.
    """
    c = np.asarray(cos_theta, dtype=float)
    return (c - 1.0) / (c + 1.0)


# ------------------------------------------------------------- TH-5 step 3 ---


def pml_slab_term(msh, x_start, thickness, s_max, order, omega_rad_s):
    """Return ``f(trial, test)``: the correction turning the vacuum form into a PML.

    Jin 3e §9.6.2 (9.152)–(9.154) anisotropic absorber in the slab
    ``x_start < x < x_start + thickness`` (normal x̂, outgoing +x), in the FE
    form of §9.6.3 (9.157):

        ∫ [(Λ⁻¹ − I)∇×u·∇×v̄ − k₀²(Λ − I)u·v̄] dx,

    Λ = diag(1/s, s, s), s = 1 − j s_max χ ξ^m, ξ = (x − x_start)/thickness
    (the (9.165) profile). χ is a DG0 indicator set in numpy at cell midpoints
    (no UFL comparison); outside the slab s = 1 and the correction vanishes.
    Vacuum only. Under e^{+jωt}, ``s_max > 0`` absorbs; ``s_max < 0`` is the
    conjugated, amplifying stretch (the γ trap's analogue), kept callable for
    the negative control. Quadrature degree pinned at 6.
    """
    from dolfinx import fem

    k0 = float(omega_rad_s) / C_0
    x_lo, x_hi = float(x_start), float(x_start) + float(thickness)
    dg0 = fem.functionspace(msh, ("DG", 0))
    chi = fem.Function(dg0)
    dtype = chi.x.array.dtype
    chi.interpolate(lambda x: ((x[0] > x_lo) & (x[0] < x_hi)).astype(dtype))
    x = ufl.SpatialCoordinate(msh)
    xi = (x[0] - x_lo) / float(thickness)
    s = 1.0 - 1j * float(s_max) * chi * xi ** int(order)
    lam = ufl.as_matrix([[1.0 / s, 0, 0], [0, s, 0], [0, 0, s]])
    lam_inv = ufl.as_matrix([[s, 0, 0], [0, 1.0 / s, 0], [0, 0, 1.0 / s]])
    eye = ufl.Identity(3)
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": 6})

    def form(trial, test):
        return (
            ufl.inner((lam_inv - eye) * ufl.curl(trial), ufl.curl(test))
            - k0 * k0 * ufl.inner((lam - eye) * trial, test)
        ) * dx

    return form


def pml_reflection_pec_backed(s_max, thickness, order, k0, cos_theta=1.0):
    """Continuum Γ_c of a PEC-backed graded PML at the front face, e^{+jωt}.

        Γ_c = −exp(−2j k₀ cosθ d) · exp(−2 k₀ cosθ s_max d / (m + 1)).

    Jin 3e §9.6.1 (9.133) (constant stretch, |R| = exp(−2 k₀ s″ L cosθ)) with
    s″L replaced by its integral over the §9.6.4 (9.165) profile
    s″(ξ) = s_max ξ^m — the drafter's derivation, not a numbered equation.
    σ form (s″ = σ/(ωε₀)): |R| = exp(−2 η₀ σ_max d cosθ / (m + 1)).
    """
    c = float(cos_theta)
    d = float(thickness)
    return complex(-np.exp(-2j * k0 * c * d) * np.exp(-2.0 * k0 * c * s_max * d / (order + 1)))


def pml_s_max_for_reflection(r, thickness, order, k0, cos_theta=1.0):
    """s_max giving design |Γ_c| = r: −ln r · (m + 1) / (2 k₀ cosθ d)."""
    return float(-np.log(r) * (order + 1) / (2.0 * k0 * float(cos_theta) * float(thickness)))
