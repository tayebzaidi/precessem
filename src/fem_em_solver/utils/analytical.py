"""Analytical solutions for electromagnetic validation cases."""

import numpy as np
from typing import Union, Callable

from .constants import MU_0


class AnalyticalSolutions:
    """Collection of analytical solutions for EM validation."""
    
    @staticmethod
    def straight_wire_magnetic_field(
        points: np.ndarray,
        current: float,
        wire_position: np.ndarray = None
    ) -> np.ndarray:
        """Analytical B-field for infinite straight wire.
        
        For a wire along the z-axis carrying current I:
            B_φ = μ₀I / (2πr)
        
        In Cartesian coordinates:
            B_x = -μ₀I * y / (2πr²)
            B_y =  μ₀I * x / (2πr²)
            B_z = 0
        
        Parameters
        ----------
        points : np.ndarray
            Array of shape (n, 3) with evaluation points
        current : float
            Current in wire [A]
        wire_position : np.ndarray, optional
            (x, y) position of wire. Default is (0, 0)
            
        Returns
        -------
        np.ndarray
            B-field at points, shape (n, 3)
        """
        if wire_position is None:
            wire_position = np.array([0.0, 0.0])
        
        mu_0 = 4 * np.pi * 1e-7  # H/m
        
        # Extract coordinates relative to wire
        x = points[:, 0] - wire_position[0]
        y = points[:, 1] - wire_position[1]
        
        # Cylindrical radius
        r = np.sqrt(x**2 + y**2)
        
        # Avoid division by zero at wire location
        r = np.maximum(r, 1e-10)
        
        # B-field magnitude (azimuthal direction)
        B_phi = mu_0 * current / (2 * np.pi * r)
        
        # Convert to Cartesian
        # B_φ direction is (-y/r, x/r, 0)
        B = np.zeros_like(points)
        B[:, 0] = -B_phi * y / r  # B_x
        B[:, 1] =  B_phi * x / r  # B_y
        B[:, 2] = 0.0              # B_z
        
        return B
    
    @staticmethod
    def straight_wire_vector_potential(
        points: np.ndarray,
        current: float,
        wire_position: np.ndarray = None,
        A_ref: float = 0.0,
        wire_radius: float = None
    ) -> np.ndarray:
        """Analytical A-field for infinite straight wire.

        The magnetic vector potential for a wire along z-axis:
            A_z = -μ₀I / (2π) * ln(r/r_ref)

        Choosing r_ref such that A=0 at some reference point.

        With ``wire_radius = a`` the uniform-current-density conductor form is
        used instead, gauged so that A_z(a) = 0::

            A_z(r) = -μ₀I/(2π) · ln(r/a)          r ≥ a
            A_z(r) =  μ₀I/(4π) · (1 - r²/a²)      r ≤ a

        Both branches have the same curl as the filament solution outside the
        conductor, and the interior branch is *finite on the axis*. That
        matters when this field is imposed as Dirichlet data (MAG-13): the
        domain end caps of ``straight_wire_domain`` cross r = 0, where the
        filament form diverges logarithmically.

        Parameters
        ----------
        points : np.ndarray
            Array of shape (n, 3) with evaluation points
        current : float
            Current in wire [A]
        wire_position : np.ndarray, optional
            (x, y) position of wire
        A_ref : float
            Reference potential (gauge choice)
        wire_radius : float, optional
            Conductor radius [m]. If given, use the finite-conductor form
            above; if None (default) the pure filament form ``-μ₀I/(2π)·ln r``
            is returned unchanged.

        Returns
        -------
        np.ndarray
            A-field at points, shape (n, 3)
        """
        if wire_position is None:
            wire_position = np.array([0.0, 0.0])

        mu_0 = 4 * np.pi * 1e-7

        x = points[:, 0] - wire_position[0]
        y = points[:, 1] - wire_position[1]
        r = np.sqrt(x**2 + y**2)

        # A has only z-component
        A = np.zeros_like(points)

        if wire_radius is None:
            # Avoid log(0)
            r = np.maximum(r, 1e-10)
            A[:, 2] = -mu_0 * current / (2 * np.pi) * np.log(r) + A_ref
            return A

        if wire_radius <= 0:
            raise ValueError(f"wire_radius must be positive, got {wire_radius!r}")

        outside = r >= wire_radius
        a_z = np.empty_like(r)
        a_z[outside] = (
            -mu_0 * current / (2 * np.pi) * np.log(r[outside] / wire_radius)
        )
        a_z[~outside] = (
            mu_0 * current / (4 * np.pi) * (1.0 - (r[~outside] / wire_radius) ** 2)
        )
        A[:, 2] = a_z + A_ref
        return A

    @staticmethod
    def circular_loop_vector_potential(
        points: np.ndarray,
        current: float,
        radius: float,
        loop_center: float = 0.0
    ) -> np.ndarray:
        """Off-axis vector potential of a circular current loop (Jackson 5.37).

        For a loop of radius ``a`` in the plane z = z0, carrying current I
        counter-clockwise seen from +z, the only non-zero component is
        azimuthal::

            A_φ(ρ,z) = (μ₀I/π)·√(a/ρ)·[(1 − k²/2)·K(k) − E(k)] / k
            k² = 4aρ / ((a+ρ)² + (z−z0)²)

        Two conventions matter here. ``scipy.special.ellipk``/``ellipe`` take
        the *parameter* ``m = k²``, not the modulus ``k`` — passing ``k`` is a
        silent factor error. And ρ → 0 needs its own branch: A_φ = 0 on the
        axis by symmetry, while the formula there is 0/0.

        Parameters
        ----------
        points : np.ndarray
            Array of shape (n, 3) with evaluation points
        current : float
            Loop current [A]
        radius : float
            Loop radius ``a`` [m]
        loop_center : float
            z-position of the loop plane [m]

        Returns
        -------
        np.ndarray
            A-field at points in Cartesian components, shape (n, 3)
        """
        from scipy.special import ellipe, ellipk

        if radius <= 0:
            raise ValueError(f"radius must be positive, got {radius!r}")

        mu_0 = 4 * np.pi * 1e-7

        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2] - loop_center
        rho = np.sqrt(x**2 + y**2)

        A = np.zeros_like(points)
        # On-axis (and at the wire itself, where the potential diverges) the
        # azimuthal direction is undefined; A_φ = 0 is the symmetry answer.
        on_axis = rho < 1e-14 * radius
        active = ~on_axis
        if not np.any(active):
            return A

        rho_a = rho[active]
        z_a = z[active]
        m = 4.0 * radius * rho_a / ((radius + rho_a) ** 2 + z_a**2)  # m = k²
        k = np.sqrt(m)
        a_phi = (
            mu_0 * current / np.pi
            * np.sqrt(radius / rho_a)
            * ((1.0 - m / 2.0) * ellipk(m) - ellipe(m))
            / k
        )

        # φ̂ = (-y, x, 0)/ρ
        A[active, 0] = -a_phi * y[active] / rho_a
        A[active, 1] = a_phi * x[active] / rho_a
        return A

    @staticmethod
    def circular_loop_magnetic_field(
        points: np.ndarray,
        current: float,
        radius: float,
        loop_center: float = 0.0
    ) -> np.ndarray:
        """Off-axis B-field of a circular current loop (elliptic integrals).

        For a loop of radius ``a`` in the plane z = z0, carrying current I
        counter-clockwise seen from +z, with α² = (a−ρ)² + (z−z0)²,
        β² = (a+ρ)² + (z−z0)² and parameter m = k² = 4aρ/β²::

            B_ρ = μ₀I (z−z0) / (2π ρ α² β) · [(a² + ρ² + (z−z0)²)·E(m) − α²·K(m)]
            B_z = μ₀I / (2π α² β) · [(a² − ρ² − (z−z0)²)·E(m) + α²·K(m)]

        As with :meth:`circular_loop_vector_potential`, scipy's ``ellipk`` /
        ``ellipe`` take the *parameter* m = k², not the modulus k. The ρ → 0
        limit reduces to :meth:`circular_loop_magnetic_field_on_axis` and is
        handled on a separate branch (B_ρ = 0 by symmetry). The field diverges
        on the filament itself (α → 0); α² is clamped there so evaluation stays
        finite for visualization, but values near the wire are not meaningful.

        Parameters
        ----------
        points : np.ndarray
            Array of shape (n, 3) with evaluation points
        current : float
            Loop current [A]
        radius : float
            Loop radius ``a`` [m]
        loop_center : float
            z-position of the loop plane [m]

        Returns
        -------
        np.ndarray
            B-field at points in Cartesian components, shape (n, 3)
        """
        from scipy.special import ellipe, ellipk

        if radius <= 0:
            raise ValueError(f"radius must be positive, got {radius!r}")

        mu_0 = 4 * np.pi * 1e-7

        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2] - loop_center
        rho = np.sqrt(x**2 + y**2)

        B = np.zeros_like(points, dtype=np.float64)

        on_axis = rho < 1e-14 * radius
        # On-axis branch: B_ρ = 0 by symmetry, closed-form B_z.
        B[on_axis, 2] = AnalyticalSolutions.circular_loop_magnetic_field_on_axis(
            points[on_axis, 2], current, radius, loop_center
        )

        active = ~on_axis
        if not np.any(active):
            return B

        rho_a = rho[active]
        z_a = z[active]

        alpha2 = (radius - rho_a) ** 2 + z_a**2
        beta = np.sqrt((radius + rho_a) ** 2 + z_a**2)
        # Clamp the filament singularity so evaluation stays finite.
        alpha2 = np.maximum(alpha2, (1e-6 * radius) ** 2)
        m = 4.0 * radius * rho_a / beta**2
        m = np.minimum(m, 1.0 - 1e-12)

        K = ellipk(m)
        E = ellipe(m)
        pref = mu_0 * current / (2 * np.pi * alpha2 * beta)
        s2 = radius**2 + rho_a**2 + z_a**2

        B_rho = pref * (z_a / rho_a) * (s2 * E - alpha2 * K)
        B_z = pref * ((radius**2 - rho_a**2 - z_a**2) * E + alpha2 * K)

        # ρ̂ = (x, y, 0)/ρ
        B[active, 0] = B_rho * x[active] / rho_a
        B[active, 1] = B_rho * y[active] / rho_a
        B[active, 2] = B_z
        return B

    @staticmethod
    def circular_loop_magnetic_field_on_axis(
        z: np.ndarray,
        current: float,
        radius: float,
        loop_center: float = 0.0
    ) -> np.ndarray:
        """B-field on axis of circular current loop.
        
        For a loop of radius a in the xy-plane, centered at z=z0:
            B_z(z) = μ₀Ia² / (2(a² + (z-z0)²)^(3/2))
        
        Parameters
        ----------
        z : np.ndarray
            Positions along z-axis [m]
        current : float
            Current in loop [A]
        radius : float
            Loop radius [m]
        loop_center : float
            z-position of loop center [m]
            
        Returns
        -------
        np.ndarray
            B_z at positions z
        """
        mu_0 = 4 * np.pi * 1e-7
        
        dz = z - loop_center
        denom = 2 * (radius**2 + dz**2)**(3/2)
        
        B_z = mu_0 * current * radius**2 / denom
        
        return B_z
    
    @staticmethod
    def helmholtz_coil_field_on_axis(
        z: np.ndarray,
        current: float,
        radius: float,
        separation: float = None
    ) -> np.ndarray:
        """B-field on axis of Helmholtz coil.
        
        Helmholtz coil: two identical loops separated by distance = radius,
        carrying current in same direction. This configuration gives
        maximally uniform field in the center.
        
        Parameters
        ----------
        z : np.ndarray
            Positions along z-axis [m]
        current : float
            Current in each loop [A]
        radius : float
            Loop radius [m]
        separation : float, optional
            Distance between loops. Default is radius (Helmholtz condition)
            
        Returns
        -------
        np.ndarray
            B_z at positions z
        """
        if separation is None:
            separation = radius  # Helmholtz condition
        
        # Center the configuration at z=0
        z1 = -separation / 2
        z2 = separation / 2
        
        B1 = AnalyticalSolutions.circular_loop_magnetic_field_on_axis(
            z, current, radius, z1
        )
        B2 = AnalyticalSolutions.circular_loop_magnetic_field_on_axis(
            z, current, radius, z2
        )
        
        return B1 + B2


_ARC_GAUSS_ORDER = 64
_ARC_GAUSS_NODES, _ARC_GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(
    _ARC_GAUSS_ORDER
)


def _finite_segment_kernel(
    points: np.ndarray, start: np.ndarray, end: np.ndarray
) -> np.ndarray:
    """B per unit current of one straight filament from ``start`` to ``end``.

    The closed form is the standard finite-segment Biot--Savart result written
    without any trigonometry (Jin 3e §2.2 in vector form)::

        B = (mu_0 I / 4 pi) * (|r1| + |r2|)
            / ( |r1||r2| (|r1||r2| + r1.r2) ) * (r1 x r2)

    with ``r1 = r - start`` and ``r2 = r - end``.  The denominator vanishes
    only on the segment itself (where ``r1 x r2 = 0`` *and* the scalar factor
    diverges); on the *extension* of the segment it stays positive and the
    cross product is zero, so those points return an exact zero.
    """
    r1 = points - start
    r2 = points - end
    n1 = np.linalg.norm(r1, axis=1)
    n2 = np.linalg.norm(r2, axis=1)
    denom = n1 * n2 * (n1 * n2 + np.einsum("ij,ij->i", r1, r2))
    if np.any(denom <= 0.0):
        raise ValueError(
            "birdcage_filament_field: an evaluation point lies on a leg "
            "filament (the Biot-Savart integrand is singular there)"
        )
    scale = (MU_0 / (4.0 * np.pi)) * (n1 + n2) / denom
    return scale[:, None] * np.cross(r1, r2)


def _arc_kernel(
    points: np.ndarray,
    radius: float,
    z_plane: float,
    phi_start: float,
    phi_end: float,
) -> np.ndarray:
    """B per unit current of one circular arc, by Gauss-Legendre quadrature.

    ``_ARC_GAUSS_ORDER`` nodes per arc: the integrand is analytic on any point
    set that keeps a finite distance from the ring, so the rule is spectrally
    accurate and the residual is at machine level for the fixtures used here.
    """
    half = 0.5 * (phi_end - phi_start)
    mid = 0.5 * (phi_end + phi_start)
    phi = mid + half * _ARC_GAUSS_NODES
    weights = half * _ARC_GAUSS_WEIGHTS

    source = np.stack(
        [radius * np.cos(phi), radius * np.sin(phi), np.full_like(phi, z_plane)],
        axis=1,
    )
    # dl/dphi = radius * phi_hat
    tangent = radius * np.stack(
        [-np.sin(phi), np.cos(phi), np.zeros_like(phi)], axis=1
    )

    delta = points[:, None, :] - source[None, :, :]
    dist = np.linalg.norm(delta, axis=2)
    if np.any(dist <= 0.0):
        raise ValueError(
            "birdcage_filament_field: an evaluation point lies on a ring arc "
            "(the Biot-Savart integrand is singular there)"
        )
    integrand = np.cross(np.broadcast_to(tangent, delta.shape), delta)
    integrand /= (dist**3)[:, :, None]
    return (MU_0 / (4.0 * np.pi)) * np.einsum("nqi,q->ni", integrand, weights)


def birdcage_filament_field(
    points: np.ndarray,
    *,
    ring_radius: float,
    coil_length: float,
    leg_currents,
    ring_currents=None,
) -> np.ndarray:
    """Biot-Savart B of an ideal filament birdcage: N legs plus two end rings.

    `WF-6` step 4a.  This is the closed-form anchor the coil-driven ``|B1+|``
    work needs: pure numpy, no mesh, no FEM (keep it that way -- importing
    dolfinx here would move the unit tier off the smoke tier).

    Geometry.  ``N = len(leg_currents)`` straight filaments of length
    ``coil_length`` sit at radius ``ring_radius`` and azimuths
    ``theta_n = 2 pi n / N``, each running from ``z = -coil_length/2`` to
    ``z = +coil_length/2``; ``leg_currents[n]`` is that leg's current in the
    **+z** direction.  Two end rings of the same radius close the circuit at
    ``z = +-coil_length/2``, each split into ``N`` arcs, arc ``n`` running from
    ``theta_n`` to ``theta_{n+1}`` in the **+phi** direction.

    Ring currents.  With ``ring_currents=None`` they come from Kirchhoff's
    current law at every leg-ring node.  At top node ``n`` the leg delivers
    ``I_n`` and the arcs carry it away::

        I_n + J_{n-1} - J_n = 0   =>   J_n = J_{n-1} + I_n

    which closes around the ring only if ``sum(I_n) = 0`` (raised otherwise),
    and then fixes ``J`` up to one constant.  The unique zero-mean solution is
    ``J = cumsum(I) - mean(cumsum(I))``.  The bottom ring's node equation is
    the same with the leg current leaving instead of arriving, so its zero-mean
    solution is exactly ``-J``: the two rings are opposite.

    Passing ``ring_currents`` explicitly overrides the Kirchhoff solve (and
    with it the zero-sum requirement, since the caller is then no longer
    claiming a closed circuit): shape ``(N,)`` gives the top-ring arc currents
    with the bottom ring at ``-``that; shape ``(2, N)`` gives top and bottom
    independently.  ``ring_currents=np.zeros(N)`` is the open-circuit negative
    control -- legs with no return path.

    Parameters
    ----------
    points : np.ndarray
        ``(n, 3)`` evaluation points [m].
    ring_radius : float
        Leg / ring radius ``R`` [m].
    coil_length : float
        Leg length ``L`` [m]; the rings sit at ``z = +-L/2``.
    leg_currents : array_like
        ``(N,)`` leg currents [A], +z positive.  Complex is allowed (a phased
        drive); the returned field is then complex.
    ring_currents : array_like, optional
        ``(N,)`` or ``(2, N)`` arc currents [A], +phi positive.  Default
        ``None`` solves Kirchhoff as above.

    Returns
    -------
    np.ndarray
        ``(n, 3)`` magnetic flux density [T], complex if any current is.
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n, 3)")
    if ring_radius <= 0.0:
        raise ValueError("ring_radius must be > 0")
    if coil_length <= 0.0:
        raise ValueError("coil_length must be > 0")

    legs = np.asarray(leg_currents)
    if legs.ndim != 1 or legs.size == 0:
        raise ValueError("leg_currents must be a non-empty 1-D array")
    n_legs = legs.size

    if ring_currents is None:
        scale = np.max(np.abs(legs))
        residual = np.abs(np.sum(legs))
        if scale > 0.0 and residual > 1e-12 * scale:
            raise ValueError(
                "leg_currents must sum to zero for the ring currents to close "
                f"around the ring (sum = {np.sum(legs)!r}); pass ring_currents "
                "explicitly to model an open circuit"
            )
        cumulative = np.cumsum(legs)
        top = cumulative - np.mean(cumulative)
        bottom = -top
    else:
        rings = np.asarray(ring_currents)
        if rings.shape == (n_legs,):
            top = rings
            bottom = -rings
        elif rings.shape == (2, n_legs):
            top, bottom = rings[0], rings[1]
        else:
            raise ValueError(
                "ring_currents must have shape (N,) or (2, N) with "
                f"N = {n_legs}; got {rings.shape}"
            )

    dtype = np.result_type(float, legs.dtype, np.asarray(top).dtype,
                           np.asarray(bottom).dtype)
    field = np.zeros((points.shape[0], 3), dtype=dtype)

    theta = 2.0 * np.pi * np.arange(n_legs) / n_legs
    half_length = 0.5 * coil_length

    for n in range(n_legs):
        if legs[n] == 0:
            continue
        foot = np.array(
            [ring_radius * np.cos(theta[n]), ring_radius * np.sin(theta[n]),
             -half_length]
        )
        head = foot + np.array([0.0, 0.0, coil_length])
        field += legs[n] * _finite_segment_kernel(points, foot, head)

    for n in range(n_legs):
        phi_start = theta[n]
        phi_end = theta[n] + 2.0 * np.pi / n_legs
        for z_plane, current in ((half_length, top[n]), (-half_length, bottom[n])):
            if current == 0:
                continue
            field += current * _arc_kernel(
                points, ring_radius, z_plane, phi_start, phi_end
            )

    return field


def dielectric_sphere_quasistatic_field(*, radius, epsilon_c, e0=1.0):
    """Quasi-static field of a (lossy) dielectric sphere in a uniform ``E₀ ẑ``.

    Exterior: uniform + dipole with ``β = (ε_c − 1)/(ε_c + 2)``; interior:
    uniform ``E_in = 3E₀/(ε_c + 2)``.  Returns a dolfinx-convention callable
    ``x (3, n) -> (3, n)`` of ``PETSc.ScalarType``.  The arithmetic is copied
    in the same operation order from ``tests/validation/test_lossy_sphere_sar.py``
    ``_exact_exterior_numpy`` (`MAT-4` step 1's wall data; `WF-15` (F) asserts
    the two bitwise), so the package does not import ``tests/``.
    """
    from petsc4py import PETSc

    beta = (epsilon_c - 1.0) / (epsilon_c + 2.0)
    e_in = 3.0 * e0 / (epsilon_c + 2.0)
    R = radius

    def field(x: np.ndarray) -> np.ndarray:
        r = np.sqrt(x[0] ** 2 + x[1] ** 2 + x[2] ** 2)
        r_safe = np.where(r > 0.0, r, 1.0)
        z = x[2]
        common = beta * e0 * R**3 / r_safe**5
        ex_out = common * 3.0 * z * x[0]
        ey_out = common * 3.0 * z * x[1]
        ez_out = e0 + common * (3.0 * z * z - r_safe**2)

        inside = r < R
        ex = np.where(inside, 0.0, ex_out)
        ey = np.where(inside, 0.0, ey_out)
        ez = np.where(inside, e_in, ez_out)
        return np.array([ex, ey, ez], dtype=PETSc.ScalarType)

    return field


def complex_permittivity(
    epsilon_r: float, sigma: float, frequency_hz: float, conjugate: bool = False
) -> complex:
    """``ε_c = εᵣ − j·σ/(ωε₀)`` in the project's ``e^{+jωt}`` convention.

    ``conjugate=True`` returns the wrong-sign value ``εᵣ + j·σ/(ωε₀)`` and
    exists only so a negative control can pass it deliberately (`TH-1`
    formulation note: the sign convention is part of the spec).
    """
    from .constants import EPSILON_0

    omega = 2.0 * np.pi * frequency_hz
    loss = sigma / (omega * EPSILON_0)
    return complex(epsilon_r, loss if conjugate else -loss)


class LossySphereSeries:
    """Full-wave series solution for a lossy dielectric sphere in a plane wave.

    ``TH-10`` step 1 anchor.  The sphere of radius ``a`` and complex relative
    permittivity ``ε_c`` sits in free space; the incident wave travels along
    ``+z`` and is polarised along ``x``.  In the project's ``e^{+jωt}``
    convention (``TH-1`` formulation note) that incident field is

        E_inc(r) = E₀ x̂ · e^{−j k₀ z},      ε_c = εᵣ − j·σ/(ωε₀)

    and the interior/scattered fields are the classical Mie series.  Textbook
    Mie theory (Bohren & Huffman, *Absorption and Scattering of Light by Small
    Particles*, ch. 4, eqs. 4.37/4.40/4.45/4.50/4.53) is written in the
    ``e^{−iωt}`` convention, where the lossy permittivity is
    ``εᵣ + i·σ/(ωε₀)``.  **This class evaluates the textbook series with the
    conjugated permittivity and conjugates the resulting field**, which is the
    convention import the ``TH-1`` note demands; the conjugated-convention
    evaluation is available as ``conjugate_convention=True`` purely as a
    negative control and is *not* physical here.  Special-function definitions
    and the Wronskian follow Jin, *The FEM in Electromagnetics* 3rd ed.,
    App. E.2 (spherical Bessel/Hankel functions, eqs. E.24–E.31); the FEM-side
    dielectric-sphere scattering fixture this anchors is Jin §9.4 (Fig. 9.11).

    Only ``numpy``/``scipy`` — no mesh, no DolfinX.

    Notes
    -----
    ``scipy.special.spherical_jn`` rejects complex arguments, so the interior
    radial functions go through ``scipy.special.jv(n+½, z)``, which does not:
    ``j_n(z) = √(π/2z)·J_{n+½}(z)``.  The exterior argument ``k₀r`` is real.
    """

    def __init__(
        self,
        radius: float,
        epsilon_c: complex,
        frequency_hz: float,
        e0: float = 1.0,
        n_terms: int = None,
        conjugate_convention: bool = False,
    ):
        from .constants import EPSILON_0, MU_0

        if radius <= 0.0:
            raise ValueError("radius must be positive")
        if frequency_hz <= 0.0:
            raise ValueError("frequency_hz must be positive")

        self.radius = float(radius)
        self.frequency_hz = float(frequency_hz)
        self.e0 = float(e0)
        self.conjugate_convention = bool(conjugate_convention)

        # e^{+jωt} spec value, kept for reporting.
        self.epsilon_c = complex(epsilon_c)
        # Textbook (e^{-iωt}) permittivity.  The negative control skips the
        # conjugation, i.e. evaluates the series as if ε_c = εᵣ + jσ/(ωε₀).
        self.epsilon_c_textbook = (
            self.epsilon_c if conjugate_convention else np.conj(self.epsilon_c)
        )

        omega = 2.0 * np.pi * self.frequency_hz
        self.k0 = omega * np.sqrt(MU_0 * EPSILON_0)
        # Relative refractive index m = k_in/k₀ = √ε_c, principal branch with
        # Re m > 0 (the physical root for a passive medium in this convention).
        m = np.sqrt(complex(self.epsilon_c_textbook))
        self.m = m if m.real >= 0.0 else -m

        self.size_parameter = self.k0 * self.radius
        self.n_terms = int(n_terms) if n_terms is not None else self._wiscombe_n()
        if self.n_terms < 1:
            raise ValueError("n_terms must be >= 1")

        self._coefficients = self._mie_coefficients()

    # ------------------------------------------------------------------ setup

    def _wiscombe_n(self) -> int:
        """Wiscombe's truncation ``x + 4x^{1/3} + 2`` on the larger of the two
        size parameters ``k₀a`` and ``|m|k₀a``, floored at 4 terms."""
        x = max(self.size_parameter, abs(self.m) * self.size_parameter)
        return int(max(4, np.ceil(x + 4.0 * x ** (1.0 / 3.0) + 2.0)))

    @staticmethod
    def _sph_jn(n: np.ndarray, z):
        """``j_n(z)`` for complex ``z`` via ``J_{n+½}``."""
        from scipy.special import jv

        z = np.asarray(z, dtype=complex)
        return np.sqrt(np.pi / (2.0 * z)) * jv(n + 0.5, z)

    @staticmethod
    def _sph_h1(n: np.ndarray, z):
        """``h_n^{(1)}(z) = √(π/2z)·H^{(1)}_{n+½}(z)`` (Jin App. E, eq. E.28)."""
        from scipy.special import hankel1

        z = np.asarray(z, dtype=complex)
        return np.sqrt(np.pi / (2.0 * z)) * hankel1(n + 0.5, z)

    @classmethod
    def _riccati(cls, n_max: int, z, kind: str):
        """Riccati-Bessel ``ψ_n(z) = z·j_n(z)`` / ``ξ_n(z) = z·h_n^{(1)}(z)``
        and its derivative, for ``n = 1…n_max``.

        The derivative uses the standard downward relation
        ``ψ_n'(z) = ψ_{n-1}(z) − n·ψ_n(z)/z`` (Jin App. E, eq. E.30 applied to
        the Riccati form), which needs the ``n = 0`` value as well.
        """
        z = np.asarray(z, dtype=complex)
        orders = np.arange(0, n_max + 1)
        shape = (n_max + 1,) + z.shape
        zz = np.broadcast_to(z, shape)
        nn = orders.reshape((-1,) + (1,) * z.ndim)
        radial = cls._sph_jn(nn, zz) if kind == "j" else cls._sph_h1(nn, zz)
        psi = zz * radial
        d_psi = psi[:-1] - nn[1:] * psi[1:] / zz[1:]
        return psi[1:], d_psi  # both indexed n = 1…n_max

    def _mie_coefficients(self) -> dict:
        """``a_n, b_n`` (scattered) and ``c_n, d_n`` (internal), B&H eq. 4.53.

        For ``m = 1`` the numerators reduce to the Wronskian
        ``ψ_nξ_n' − ξ_nψ_n' = i`` and ``c_n = d_n = 1`` exactly — that identity
        is what the empty-limit gate exercises.
        """
        n_max = self.n_terms
        x = np.array([self.size_parameter], dtype=complex)
        m = self.m

        psi_x, dpsi_x = self._riccati(n_max, x, "j")
        xi_x, dxi_x = self._riccati(n_max, x, "h1")
        psi_mx, dpsi_mx = self._riccati(n_max, m * x, "j")

        psi_x, dpsi_x = psi_x[:, 0], dpsi_x[:, 0]
        xi_x, dxi_x = xi_x[:, 0], dxi_x[:, 0]
        psi_mx, dpsi_mx = psi_mx[:, 0], dpsi_mx[:, 0]

        a_n = (m * psi_mx * dpsi_x - psi_x * dpsi_mx) / (
            m * psi_mx * dxi_x - xi_x * dpsi_mx
        )
        b_n = (psi_mx * dpsi_x - m * psi_x * dpsi_mx) / (
            psi_mx * dxi_x - m * xi_x * dpsi_mx
        )
        wronskian = m * (psi_x * dxi_x - xi_x * dpsi_x)
        c_n = wronskian / (psi_mx * dxi_x - m * xi_x * dpsi_mx)
        d_n = wronskian / (m * psi_mx * dxi_x - xi_x * dpsi_mx)
        return {"a": a_n, "b": b_n, "c": c_n, "d": d_n}

    @property
    def coefficients(self) -> dict:
        return self._coefficients

    def last_term_bound(self) -> float:
        """``|E_N|·max(|c_N|,|d_N|)·|j_N(m k₀a)|/E₀`` — the dropped tail.

        The radial factor is what actually kills the tail for a small sphere
        (``j_N(z) ~ z^N/(2N+1)!!``), so leaving it out reports a bound of order
        one for the empty limit and reads as alarming when the true tail is
        ``1e-9``.  Printed rather than assumed: a sweep in ``n_terms`` is the
        honest check, and this is the cheap companion to it.
        """
        n = self.n_terms
        e_n = (2.0 * n + 1.0) / (n * (n + 1.0))
        radial = abs(
            complex(
                self._sph_jn(
                    np.array([n]), np.array([self.m * self.k0 * self.radius])
                )[0]
            )
        )
        return float(
            e_n
            * max(abs(self._coefficients["c"][-1]), abs(self._coefficients["d"][-1]))
            * radial
        )

    # ------------------------------------------------------------------ fields

    @staticmethod
    def _angular_functions(n_max: int, mu: np.ndarray):
        """``π_n(cosθ) = P_n¹/sinθ`` and ``τ_n = dP_n¹/dθ`` by upward recurrence
        (B&H eq. 4.47); both are finite on the axis."""
        pi_n = np.zeros((n_max + 1,) + mu.shape)
        tau_n = np.zeros((n_max + 1,) + mu.shape)
        pi_n[1] = 1.0
        tau_n[1] = mu
        for n in range(2, n_max + 1):
            pi_n[n] = ((2 * n - 1) / (n - 1)) * mu * pi_n[n - 1] - (n / (n - 1)) * pi_n[
                n - 2
            ]
            tau_n[n] = n * mu * pi_n[n] - (n + 1) * pi_n[n - 1]
        return pi_n[1:], tau_n[1:]

    def _spherical_frame(self, points: np.ndarray):
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 3:
            raise ValueError("points must have shape (n, 3)")
        r = np.linalg.norm(points, axis=1)
        if np.any(r < 1e-14):
            raise ValueError(
                "the series is evaluated in spherical coordinates; r = 0 is a "
                "removable singularity that this implementation does not take "
                "the limit of — offset the probe point"
            )
        mu = points[:, 2] / r
        mu = np.clip(mu, -1.0, 1.0)
        sin_theta = np.sqrt(np.maximum(0.0, 1.0 - mu**2))
        phi = np.arctan2(points[:, 1], points[:, 0])
        cos_phi, sin_phi = np.cos(phi), np.sin(phi)
        e_r = np.stack(
            [sin_theta * cos_phi, sin_theta * sin_phi, mu], axis=1
        )
        e_theta = np.stack(
            [mu * cos_phi, mu * sin_phi, -sin_theta], axis=1
        )
        e_phi = np.stack([-sin_phi, cos_phi, np.zeros_like(phi)], axis=1)
        return r, mu, sin_theta, cos_phi, sin_phi, e_r, e_theta, e_phi

    def _series(self, points: np.ndarray, which: str) -> np.ndarray:
        """Evaluate one of ``incident|internal|scattered`` at ``points``.

        Returns the textbook (``e^{-iωt}``) field; the public methods conjugate.
        """
        n_max = self.n_terms
        r, mu, sin_theta, cos_phi, sin_phi, e_r, e_theta, e_phi = (
            self._spherical_frame(points)
        )
        pi_n, tau_n = self._angular_functions(n_max, mu)
        n = np.arange(1, n_max + 1).reshape(-1, 1)

        if which == "internal":
            rho = self.m * self.k0 * r
            psi, d_psi = self._riccati(n_max, rho, "j")
            coeff_m, coeff_n = self._coefficients["c"], self._coefficients["d"]
        elif which == "incident":
            rho = (self.k0 * r).astype(complex)
            psi, d_psi = self._riccati(n_max, rho, "j")
            coeff_m = np.ones(n_max, dtype=complex)
            coeff_n = np.ones(n_max, dtype=complex)
        elif which == "scattered":
            rho = (self.k0 * r).astype(complex)
            psi, d_psi = self._riccati(n_max, rho, "h1")
            coeff_m, coeff_n = self._coefficients["b"], self._coefficients["a"]
        else:
            raise ValueError(f"unknown series {which!r}")

        z_n = psi / rho  # the radial function itself
        d_n = d_psi / rho  # [ρ z_n(ρ)]'/ρ
        e_n = 1j**n.astype(float) * self.e0 * (2.0 * n + 1.0) / (n * (n + 1.0))
        coeff_m = coeff_m.reshape(-1, 1)
        coeff_n = coeff_n.reshape(-1, 1)

        # M_o1n and N_e1n, B&H eq. 4.50.
        m_theta = cos_phi * pi_n * z_n
        m_phi = -sin_phi * tau_n * z_n
        n_r = cos_phi * n * (n + 1.0) * sin_theta * pi_n * z_n / rho
        n_theta = cos_phi * tau_n * d_n
        n_phi = -sin_phi * pi_n * d_n

        if which == "scattered":
            # E_s = Σ E_n (i a_n N^(3) − b_n M^(3));  coeff_m = b, coeff_n = a.
            w_m, w_n = -coeff_m, 1j * coeff_n
        else:
            # E = Σ E_n (c_n M^(1) − i d_n N^(1)); incident is c = d = 1.
            w_m, w_n = coeff_m, -1j * coeff_n

        comp_r = np.sum(e_n * w_n * n_r, axis=0)
        comp_theta = np.sum(e_n * (w_m * m_theta + w_n * n_theta), axis=0)
        comp_phi = np.sum(e_n * (w_m * m_phi + w_n * n_phi), axis=0)

        return (
            comp_r[:, None] * e_r
            + comp_theta[:, None] * e_theta
            + comp_phi[:, None] * e_phi
        )

    def incident_field(self, points: np.ndarray) -> np.ndarray:
        """Series form of the incident plane wave, shape (n, 3), complex."""
        return np.conj(self._series(points, "incident"))

    def incident_field_closed_form(self, points: np.ndarray) -> np.ndarray:
        """``E₀ x̂ e^{−j k₀ z}`` — what the incident series must reproduce."""
        points = np.asarray(points, dtype=float)
        field = np.zeros((points.shape[0], 3), dtype=complex)
        field[:, 0] = self.e0 * np.exp(-1j * self.k0 * points[:, 2])
        return field

    def internal_field(self, points: np.ndarray) -> np.ndarray:
        """Interior field (valid for ``r < a``), shape (n, 3), complex."""
        return np.conj(self._series(points, "internal"))

    def scattered_field(self, points: np.ndarray) -> np.ndarray:
        """Scattered field (valid for ``r > a``), shape (n, 3), complex."""
        return np.conj(self._series(points, "scattered"))

    def total_field(self, points: np.ndarray) -> np.ndarray:
        """Piecewise total field: interior series inside, incident + scattered
        outside.  This is the callable a later ``TH-10`` step drives the box
        wall with, exactly as ``TH-8`` drives its box with the quasi-static
        closed form."""
        points = np.asarray(points, dtype=float)
        r = np.linalg.norm(points, axis=1)
        inside = r < self.radius
        field = np.zeros((points.shape[0], 3), dtype=complex)
        if np.any(inside):
            field[inside] = self.internal_field(points[inside])
        outside = ~inside
        if np.any(outside):
            field[outside] = self.incident_field(points[outside]) + (
                self.scattered_field(points[outside])
            )
        return field

    def quasistatic_internal_field(self) -> complex:
        """``3E₀/(ε_c + 2)`` — the ``TH-8`` closed form continued onto the
        imaginary axis of ``ε_c``.  The full-wave interior field must approach
        this (uniform, x̂-directed) value as ``|m|k₀a → 0``."""
        return 3.0 * self.e0 / (self.epsilon_c + 2.0)


class ErrorMetrics:
    """Error metrics for comparing numerical and analytical solutions."""
    
    @staticmethod
    def l2_error(numerical: np.ndarray, analytical: np.ndarray) -> float:
        """Compute L2 norm of error.
        
        ||u_num - u_ana||₂ = sqrt(∫ |u_num - u_ana|² dx)
        
        For discrete points:
            L2_error = sqrt(sum(|u_num - u_ana|²))
        """
        diff = numerical - analytical
        return np.sqrt(np.sum(diff**2))
    
    @staticmethod
    def l2_relative_error(numerical: np.ndarray, analytical: np.ndarray) -> float:
        """Compute relative L2 error."""
        l2_err = ErrorMetrics.l2_error(numerical, analytical)
        l2_ana = np.sqrt(np.sum(analytical**2))
        
        if l2_ana < 1e-15:
            return l2_err
        
        return l2_err / l2_ana
    
    @staticmethod
    def max_error(numerical: np.ndarray, analytical: np.ndarray) -> float:
        """Compute maximum absolute error (L∞ norm)."""
        return np.max(np.abs(numerical - analytical))
    
    @staticmethod
    def max_relative_error(numerical: np.ndarray, analytical: np.ndarray) -> float:
        """Compute maximum relative error."""
        abs_err = np.abs(numerical - analytical)
        max_abs_ana = np.max(np.abs(analytical))
        
        if max_abs_ana < 1e-15:
            return np.max(abs_err)
        
        return np.max(abs_err / (np.abs(analytical) + 1e-15))


def robin_sphere_steady_temperature(
    r,
    *,
    radius: float,
    k: float,
    q: float,
    h: float,
    t_inf: float = 0.0,
):
    """Steady temperature in a uniformly heated, convectively cooled ball (`MAT-9`).

    Solves ``-k ∇²T = Q`` in ``r < R`` with the Robin (Newton cooling)
    condition ``-k ∂T/∂r = h (T - T∞)`` on ``r = R``.

    Derivation.  Spherical symmetry reduces the equation to
    ``-(k/r²) d/dr (r² dT/dr) = Q``.  Integrating once, ``r² dT/dr =
    -Q r³/(3k) + C₁``, and regularity at ``r = 0`` forces ``C₁ = 0``, so
    ``dT/dr = -Q r/(3k)`` and ``T = C₂ - Q r²/(6k)``.  At ``r = R`` the heat
    flux ``-k dT/dr = QR/3`` (which is also ``∫Q dV / 4πR²`` — the global
    balance) must equal ``h (T(R) - T∞)``, hence ``T(R) = T∞ + QR/(3h)`` and

    .. math::

        T(r) = T_∞ + \\frac{Q R}{3 h} + \\frac{Q (R^2 - r^2)}{6 k}.

    This is the solid-sphere-with-uniform-generation result of Incropera et
    al., *Fundamentals of Heat and Mass Transfer*, §3.5 / Appendix C
    (``T(r) = Q R²/(6k) (1 - r²/R²) + T_s``) combined with the surface energy
    balance ``T_s = T∞ + QR/(3h)`` (``q'' = QR/3``); Carslaw & Jaeger,
    *Conduction of Heat in Solids*, ch. IX (sphere with heat production,
    radiation boundary) gives the same steady limit.  Equation numbers are
    NOT pinned here: no copy of either text was available in the `MAT-9`
    step-1 slot (2026-09-28), and the derivation above is self-contained;
    pinning them is left open rather than guessed.

    Parameters
    ----------
    r : array_like
        Radial coordinate(s) [m], ``0 ≤ r ≤ R``.
    radius, k, q, h, t_inf
        ``R`` [m], conductivity [W/m/K], volumetric source [W/m³], film
        coefficient [W/m²/K], ambient temperature.

    Returns
    -------
    ndarray
        ``T(r)`` in the units of ``t_inf``.
    """
    r = np.asarray(r, dtype=float)
    R = float(radius)
    return t_inf + q * R / (3.0 * h) + q * (R**2 - r**2) / (6.0 * k)


def robin_sphere_linear_source_coefficients(*, radius, k, q0, alpha, h):
    """``(m, D, A, λ₁)`` of :func:`robin_sphere_linear_source_temperature`.

    ``m = √(αq₀/k)``, ``D = mR cos mR + (Bi − 1) sin mR``,
    ``A = Bi·R/(αD)`` (so ``θ = A sin(mr)/r − 1/α``) and ``λ₁`` the first
    root of `MAT-9`'s eigen-equation.  Raises ``ValueError`` if ``α ≤ 0`` or
    ``mR ≥ λ₁`` (thermal runaway: no bounded steady state).
    """
    R = float(radius)
    if alpha <= 0.0 or q0 <= 0.0 or k <= 0.0 or h <= 0.0:
        raise ValueError("radius, k, q0, alpha, h must be positive")
    bi = h * R / k
    m = np.sqrt(alpha * q0 / k)
    lam1 = float(robin_sphere_eigenvalues(bi, 1)[0])
    if m * R >= lam1:
        raise ValueError(
            f"m*R = {m * R!r} >= lambda_1 = {lam1!r}: thermal runaway, "
            "no steady state"
        )
    D = m * R * np.cos(m * R) + (bi - 1.0) * np.sin(m * R)
    A = bi * R / (alpha * D)
    return float(m), float(D), float(A), lam1


def robin_sphere_linear_source_temperature(
    r,
    *,
    radius: float,
    k: float,
    q0: float,
    alpha: float,
    h: float,
    t_inf: float = 0.0,
):
    """Robin-cooled ball with a source linear in temperature (`MAT-5` step 1).

    Solves ``-k ∇²T = q₀[1 + α(T − T∞)]`` in ``r < R`` with
    ``-k ∂T/∂r = h(T − T∞)`` on ``r = R`` — the fixed-\\|E\\| Joule source
    ``½σ(T)|E|²`` with ``σ(T) = σ₀[1 + α(T − T₀)]`` and ``T₀ = T∞``.

    Derivation.  With ``θ = T − T∞`` and ``m² = αq₀/k`` the equation is the
    spherical Helmholtz problem ``θ'' + (2/r)θ' + m²θ = −q₀/k``.  The bounded
    solution is ``θ = A sin(mr)/r − 1/α`` (particular part ``−1/α``;
    ``cos(mr)/r`` is excluded by regularity).  With
    ``θ'(r) = A (mr cos mr − sin mr)/r²`` the Robin condition
    ``−kθ'(R) = hθ(R)``, multiplied by ``R²/k``, reads
    ``−A(mR cos mR − sin mR) = Bi·A sin mR − Bi·R/α`` (``Bi = hR/k``), so
    ``A = Bi·R/(αD)`` with ``D = mR cos mR + (Bi − 1) sin mR`` and

    .. math::

        θ(r) = \\frac{1}{α}\\left[\\frac{Bi\\,R \\sin(mr)}{r\\,D} - 1\\right].

    ``D = 0`` is `MAT-9`'s eigen-equation ``(1 − Bi) sin λ = λ cos λ``, so
    the steady state exists only for ``mR < λ₁`` (runaway at ``mR = λ₁``).
    As ``α → 0`` it reduces to :func:`robin_sphere_steady_temperature`.
    This is the steady case of Carslaw & Jaeger, *Conduction of Heat in
    Solids*, ch. IX (heat production a linear function of temperature);
    the equation number is NOT pinned — no copy of the text was available in
    the `MAT-5` step-1 slot (2026-09-30), and the self-contained derivation
    above stands in for it (as `MAT-9` item 66 did).

    ``sin(mr)/r → m`` at ``r = 0``.  Raises ``ValueError`` if
    ``mR ≥ robin_sphere_eigenvalues(Bi, 1)[0]``.
    """
    r = np.asarray(r, dtype=float)
    m, D, A, _ = robin_sphere_linear_source_coefficients(
        radius=radius, k=k, q0=q0, alpha=alpha, h=h
    )
    safe = np.where(r > 0.0, r, 1.0)
    s = np.where(r > 0.0, np.sin(m * safe) / safe, m)
    return t_inf + A * s - 1.0 / alpha


def _robin_sphere_perfusion_scaled(radius, k, q, w, h):
    """``(m, x, Bi, Dn)`` with ``Dn = 2 e^{−x} D`` (overflow-free), ``x = mR``."""
    R = float(radius)
    if w <= 0.0 or q <= 0.0 or k <= 0.0 or h <= 0.0 or R <= 0.0:
        raise ValueError("radius, k, q, w, h must be positive")
    bi = h * R / k
    m = float(np.sqrt(w / k))
    x = m * R
    e2 = np.exp(-2.0 * x)
    dn = x * (1.0 + e2) + (bi - 1.0) * (-np.expm1(-2.0 * x))
    return m, x, bi, float(dn)


def robin_sphere_perfusion_coefficients(*, radius, k, q, w, h):
    """``(m, D, A, θ̄)`` of :func:`robin_sphere_perfusion_temperature` (`MAT-11`).

    ``m = √(w/k)``, ``D = mR cosh mR + (Bi − 1) sinh mR`` (positive for every
    ``mR > 0``, ``Bi > 0``: a sink has no runaway), ``A = Bi·mR/D`` and the
    volume mean ``θ̄ = (q/w)[1 − 3A(x cosh x − sinh x)/x³]``, ``x = mR``.
    ``θ̄`` is evaluated in a scaled form that does not overflow; ``D`` itself
    overflows to ``inf`` (and ``A`` to 0) past ``mR ≈ 710``.  Raises
    ``ValueError`` for ``w ≤ 0``.
    """
    m, x, bi, dn = _robin_sphere_perfusion_scaled(radius, k, q, w, h)
    with np.errstate(over="ignore"):
        D = x * np.cosh(x) + (bi - 1.0) * np.sinh(x)
    A = bi * x / D
    e2 = np.exp(-2.0 * x)
    # 3A(x cosh x − sinh x)/x³ = 3 Bi [x(1 + e^{−2x}) − (1 − e^{−2x})] / (x² Dn)
    bracket = 3.0 * bi * (x * (1.0 + e2) + np.expm1(-2.0 * x)) / (x * x * dn)
    return m, float(D), float(A), float(q / w * (1.0 - bracket))


def robin_sphere_perfusion_temperature(
    r,
    *,
    radius: float,
    k: float,
    q: float,
    w: float,
    h: float,
    t_inf: float = 0.0,
):
    r"""Robin-cooled, uniformly heated ball with a Pennes perfusion sink (`MAT-11`).

    Solves ``-k ∇²T + w (T − T∞) = q`` in ``r < R`` with
    ``-k ∂T/∂r = h (T − T∞)`` on ``r = R``: the Pennes (1948) bioheat equation
    at steady state with the arterial temperature equal to the ambient one and
    ``w = ω_b ρ_b c_b`` [W/m³/K].

    Derivation.  With ``θ = T − T∞`` and ``m² = w/k`` the equation is the
    modified spherical Helmholtz problem ``θ'' + (2/r)θ' − m²θ = −q/k``.  The
    particular part is ``q/w``; the bounded homogeneous part is
    ``sinh(mr)/r`` (``cosh(mr)/r`` is excluded by regularity), so
    ``θ = q/w − C sinh(mr)/r``.  With ``θ' = −C (mr cosh mr − sinh mr)/r²``
    the Robin condition ``−kθ'(R) = hθ(R)``, multiplied by ``R²/k``, reads
    ``C (mR cosh mR − sinh mR) = Bi·R·q/w − Bi·C sinh mR`` (``Bi = hR/k``),
    hence ``C = Bi·R·(q/w)/D`` with ``D = mR cosh mR + (Bi − 1) sinh mR`` and

    .. math::

        θ(r) = \frac{q}{w}\left[1 - A\,\frac{\sinh(mr)}{mr}\right],
        \qquad A = \frac{Bi\,mR}{D}.

    It is :func:`robin_sphere_linear_source_temperature` continued to
    ``α = −w/q`` (``sin → sinh``); ``D > 0`` always, so there is no runaway.
    As ``w → 0`` it reduces to :func:`robin_sphere_steady_temperature` (the
    bracket then loses ``eps/(mR)²`` digits); as ``mR → ∞`` the interior tends
    to ``q/w`` and ``wθ(R)/q → (mR − 1)/(mR + Bi − 1)``.  The global balance is
    ``h·4πR²·θ(R) + w·V·θ̄ = q·V``.

    This is the steady sphere with heat production and a linear loss term of
    Carslaw & Jaeger, *Conduction of Heat in Solids*, ch. IX, and the bioheat
    section of Incropera et al.; equation numbers are NOT pinned — no copy of
    either text was available in the `MAT-11` slot (2026-09-30), and the
    self-contained derivation above stands in (as for `MAT-9` and `MAT-5`).

    The ratio ``sinh(mr)/D`` is evaluated as
    ``e^{−m(R−r)}(1 − e^{−2mr})/Dn`` so nothing overflows at large ``mR``.
    Raises ``ValueError`` for ``w ≤ 0``.
    """
    r = np.asarray(r, dtype=float)
    R = float(radius)
    m, x, bi, dn = _robin_sphere_perfusion_scaled(radius, k, q, w, h)
    safe = np.where(r > 0.0, r, 1.0)
    # (1 − e^{−2mr})/r → 2m at r = 0.
    s = np.where(r > 0.0, -np.expm1(-2.0 * m * safe) / safe, 2.0 * m)
    ratio = bi * R * np.exp(-m * (R - r)) * s / dn
    return t_inf + q / w * (1.0 - ratio)


def robin_sphere_eigenvalues(biot: float, n: int):
    """First ``n`` positive roots of ``(1 − Bi) sin λ − λ cos λ = 0`` (`MAT-9`).

    Equivalent to ``λ cot λ = 1 − Bi``: the radial eigenvalues of the
    Robin-cooled ball (Carslaw & Jaeger, *Conduction of Heat in Solids*,
    ch. IX; Incropera et al., *Fundamentals of Heat and Mass Transfer*,
    §5.6, where the same equation is written ``1 − ζₙ cot ζₙ = Bi``).  For
    ``Bi > 0`` exactly one root lies in each ``((n − 1)π, nπ)``; it is
    bracketed there and found with ``scipy.optimize.brentq``.  At ``Bi = 1``
    the roots are ``(n − ½)π``.
    """
    from scipy.optimize import brentq

    if biot <= 0.0:
        raise ValueError(f"biot must be positive, got {biot!r}")

    def f(lam):
        return (1.0 - biot) * np.sin(lam) - lam * np.cos(lam)

    roots = np.empty(int(n))
    for i in range(1, int(n) + 1):
        lo = (i - 1) * np.pi + 1e-12
        hi = i * np.pi - 1e-12
        roots[i - 1] = brentq(f, lo, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps,
                              maxiter=500)
    return roots


def robin_sphere_transient_coefficients(lam, *, radius, k, q, h):
    """``(Cₙ, Bₙ)`` of :func:`robin_sphere_transient_temperature`.

    ``Cₙ`` projects the constant 1 and ``Bₙ`` projects the steady excess
    ``θ_ss`` onto ``Xₙ(ρ) = sin(λₙρ)/(λₙρ)`` with weight ``ρ²``.
    """
    lam = np.asarray(lam, dtype=float)
    C = 4.0 * (np.sin(lam) - lam * np.cos(lam)) / (2.0 * lam - np.sin(2.0 * lam))
    B = q * radius**2 / k * C / lam**2
    return C, B


def robin_sphere_transient_temperature(
    r,
    t,
    *,
    radius: float,
    k: float,
    rho_c: float,
    q: float,
    h: float,
    t_inf: float,
    t_init: float,
    n_terms: int = 400,
):
    """Transient temperature in a Robin-cooled ball with a step source (`MAT-9`).

    Solves ``ρc ∂ₜT = k∇²T + Q`` in ``r < R``, ``−k ∂ᵣT = h(T − T∞)`` on
    ``r = R``, ``T(r, 0) = T_init`` (uniform), ``Q`` switched on at ``t = 0``.

    Derivation.  Let ``θ = T − T∞`` and ``θ_ss`` the steady excess of
    :func:`robin_sphere_steady_temperature`.  ``u = θ − θ_ss`` satisfies the
    source-free equation with the homogeneous Robin condition and
    ``u(r, 0) = θᵢ − θ_ss(r)``, ``θᵢ = T_init − T∞``.  Separation in
    ``ρ = r/R`` and Fourier number ``Fo = αt/R²`` (``α = k/ρc``) gives the
    radial eigenfunctions ``Xₙ = sin(λₙρ)/(λₙρ)``; the Robin condition
    ``−k X'(1)/R = h X(1)`` reduces to ``(1 − Bi) sin λ − λ cos λ = 0``,
    ``Bi = hR/k`` (:func:`robin_sphere_eigenvalues`).  The ``Xₙ`` are
    orthogonal with weight ``ρ²`` and

        ∫₀¹ρ²Xₙ dρ = (sin λₙ − λₙ cos λₙ)/λₙ³,
        ∫₀¹ρ²Xₙ² dρ = (2λₙ − sin 2λₙ)/(4λₙ³),

    so the projection of 1 is ``Cₙ = 4(sin λₙ − λₙ cos λₙ)/(2λₙ − sin 2λₙ)``.
    For ``θ_ss`` Green's second identity (both ``θ_ss`` and ``Xₙ`` satisfy the
    homogeneous Robin condition, ``−∇²θ_ss = Q/k``, ``−∇²Xₙ = λₙ²Xₙ/R²``) gives
    ``Bₙ = (QR²/k) Cₙ/λₙ²``.  Hence

    .. math::

        θ(r, t) = θ_{ss}(r) + \\sum_n (θ_i C_n − B_n)\\,
                  e^{-λ_n^2 α t / R^2}\\,\\frac{\\sin(λ_n r/R)}{λ_n r/R}.

    This is the sphere-with-radiation-boundary and heat-production solution
    of Carslaw & Jaeger, *Conduction of Heat in Solids*, ch. IX, and the
    one-term / series sphere solution of Incropera et al., *Fundamentals of
    Heat and Mass Transfer*, §5.6 (``Cₙ`` as written there) superposed on
    the steady generation solution.  Equation numbers are not pinned
    (neither text is in ``docs/references/``); the derivation above is
    self-contained.

    Returns ``T(r, t)`` for broadcastable ``r`` [m] and ``t`` [s].
    """
    r = np.asarray(r, dtype=float)
    t = np.asarray(t, dtype=float)
    R = float(radius)
    lam = robin_sphere_eigenvalues(h * R / k, n_terms)
    C, B = robin_sphere_transient_coefficients(lam, radius=R, k=k, q=q, h=h)
    theta_i = t_init - t_inf
    alpha = k / rho_c
    rho, tt = np.broadcast_arrays(r / R, t)
    rho = rho[..., None]
    tt = tt[..., None]
    x = lam * rho
    with np.errstate(invalid="ignore", divide="ignore"):
        X = np.where(x == 0.0, 1.0, np.sin(x) / np.where(x == 0.0, 1.0, x))
    series = np.sum((theta_i * C - B) * np.exp(-lam**2 * alpha * tt / R**2) * X, axis=-1)
    theta_ss = robin_sphere_steady_temperature(
        rho[..., 0] * R, radius=R, k=k, q=q, h=h, t_inf=0.0
    )
    return t_inf + theta_ss + series


def robin_sphere_perfusion_transient_coefficients(lam, *, radius, k, q, w, h):
    """``(Cₙ, Bₙ^w)`` of :func:`robin_sphere_perfusion_transient_temperature`.

    ``Bₙ^w`` projects the *perfused* steady excess θ_ss^w onto ``Xₙ``:
    Green's second identity on θ_ss^w (``−∇²θ + m²θ = q/k``) and ``Xₙ``
    (``−∇²Xₙ = λₙ²Xₙ/R²``), both under the homogeneous Robin condition, gives
    ``Bₙ^w = (qR²/k) Cₙ/(λₙ² + x²)``, ``x² = wR²/k``; ``w = 0`` is `MAT-9`'s ``Bₙ``.
    """
    lam = np.asarray(lam, dtype=float)
    C, _ = robin_sphere_transient_coefficients(lam, radius=radius, k=k, q=q, h=h)
    x2 = w * radius**2 / k
    return C, q * radius**2 / k * C / (lam**2 + x2)


def robin_sphere_perfusion_transient_temperature(
    r,
    t,
    *,
    radius: float,
    k: float,
    rho_c: float,
    q: float,
    w: float,
    h: float,
    t_inf: float,
    t_init: float,
    n_terms: int = 400,
):
    """Transient Pennes ball: step source, perfusion sink, Robin cooling (`MAT-13`).

    Solves ``ρc ∂ₜT = k∇²T − w(T − T∞) + q`` in ``r < R``,
    ``−k ∂ᵣT = h(T − T∞)`` on ``r = R``, ``T(r, 0) = T_init``, arterial
    temperature equal to ``T∞``.  ``u = θ − θ_ss^w`` (θ_ss^w from
    :func:`robin_sphere_perfusion_temperature`) is source-free,
    ``ρc ∂ₜu = k∇²u − wu``, so on `MAT-9`'s ``Xₙ = sin(λₙρ)/(λₙρ)`` mode ``n``
    decays at ``(λₙ² + x²) k/(ρcR²)``, ``x = mR``::

        θ(r, t) = θ_ss^w(r) + Σ (θᵢCₙ − Bₙ^w) e^{−(λₙ² + x²)Fo} Xₙ(r/R),

    ``Fo = kt/(ρcR²)``, ``Bₙ^w`` from
    :func:`robin_sphere_perfusion_transient_coefficients`.  No textbook
    equation is pinned (``docs/references/`` holds Jin only); the derivation
    above and `MAT-13`'s numpy anchors (quad, the ``w → 0`` limit, the
    time-convolution identity) stand in.  Raises ``ValueError`` for ``w ≤ 0``.
    """
    if w <= 0.0:
        raise ValueError(f"w must be positive, got {w!r}")
    r = np.asarray(r, dtype=float)
    t = np.asarray(t, dtype=float)
    R = float(radius)
    lam = robin_sphere_eigenvalues(h * R / k, n_terms)
    C, B = robin_sphere_perfusion_transient_coefficients(lam, radius=R, k=k, q=q, w=w, h=h)
    theta_i = t_init - t_inf
    x2 = w * R**2 / k
    alpha = k / rho_c
    rho, tt = np.broadcast_arrays(r / R, t)
    rho = rho[..., None]
    tt = tt[..., None]
    x = lam * rho
    with np.errstate(invalid="ignore", divide="ignore"):
        X = np.where(x == 0.0, 1.0, np.sin(x) / np.where(x == 0.0, 1.0, x))
    series = np.sum(
        (theta_i * C - B) * np.exp(-(lam**2 + x2) * alpha * tt / R**2) * X, axis=-1
    )
    theta_ss = robin_sphere_perfusion_temperature(
        rho[..., 0] * R, radius=R, k=k, q=q, w=w, h=h, t_inf=0.0
    )
    return t_inf + theta_ss + series


# ---------------------------------------------------------------------------
# IMP-1 step 3(a): exact finite-k l = 1 TM field of a PEC sphere in a lossy
# host (§9 item 135).  e^{+jωt}; k = k₀√ε_c on the principal root (Im k < 0).
# Spherical Bessel/Hankel definitions: Jin, *The FEM in Electromagnetics* 3rd
# ed., App. E.2 (E.24–E.31).  The outgoing radial function is h₁⁽²⁾, which is
# outgoing in e^{+jωt}; the regular part uses j₁ of the *complex* argument via
# J_{3/2} (``spherical_jn`` rejects complex input and the j₁ closed form
# cancels at |z| < 0.5).  Never build the regular part as conj(outgoing): that
# identity holds only for real k.
# ---------------------------------------------------------------------------


def _l1_sph_j1(z):
    """``(j₁(z), (z j₁)′(z))`` for complex ``z`` via ``J_{1/2}``, ``J_{3/2}``."""
    from scipy.special import jv

    z = np.asarray(z, dtype=complex)
    pref = np.sqrt(np.pi / (2.0 * z))
    j0 = pref * jv(0.5, z)
    j1 = pref * jv(1.5, z)
    return j1, z * j0 - j1


def _l1_sph_h2(z):
    """``(h₁⁽²⁾(z), (z h₁⁽²⁾)′(z))``; ``h₁⁽²⁾ = e^{−jz}(j/z² − 1/z)``."""
    z = np.asarray(z, dtype=complex)
    ph = np.exp(-1j * z)
    h = ph * (1j / z**2 - 1.0 / z)
    dzh = ph * (1.0 / z + 1j - 1j / z**2)
    return h, dzh


def _l1_tm_bracket(points, k, kind):
    """``2cosθ f(x)/x r̂ − sinθ (x f)′/x θ̂`` with x = k r, in the axis-safe
    Cartesian form ``cosθ [2f − (xf)′]/x r̂ + (xf)′/x ẑ``."""
    pts = np.asarray(points, dtype=float)
    r = np.linalg.norm(pts, axis=1)
    r_hat = pts / r[:, None]
    cos_t = r_hat[:, 2]
    x = k * r
    f, dxf = _l1_sph_j1(x) if kind == "j" else _l1_sph_h2(x)
    radial = cos_t * (2.0 * f - dxf) / x
    out = radial[:, None] * r_hat
    out[:, 2] = out[:, 2] + dxf / x
    return out


def pec_sphere_l1_tm_parts(points, *, k, radius):
    """Unmasked ``(E_reg, E_out)`` at ``points`` (N, 3), E₀ = 1.

    ``E_reg = (3/2)[…j₁…]`` → ẑ as k → 0; ``E_out = −j(ka)³[…h₁⁽²⁾…]`` → the
    a³ static dipole ``(a/r)³(3cosθ r̂ − ẑ)`` as k → 0.
    """
    k = complex(k)
    e_reg = 1.5 * _l1_tm_bracket(points, k, "j")
    e_out = -1j * (k * radius) ** 3 * _l1_tm_bracket(points, k, "h")
    return e_reg, e_out


def pec_sphere_l1_tm_beta(k, radius):
    """``β_k = −(3/2)(x j₁)′ / [−j x³ (x h₁⁽²⁾)′]`` at x = ka, so that
    ``n × (E_reg + β_k E_out) = 0`` on r = a; β_k → 1 as k → 0."""
    x = complex(k) * float(radius)
    _, dxj = _l1_sph_j1(x)
    _, dxh = _l1_sph_h2(x)
    return complex(-1.5 * dxj / (-1j * x**3 * dxh))


def pec_sphere_l1_tm_field(points, *, k, radius, e0=1.0, beta=None):
    """Exact l = 1 TM field ``E₀(E_reg + β_k E_out)`` of a PEC sphere of radius
    ``radius`` in a homogeneous (lossy) host of wavenumber ``k`` under a
    uniform ẑ-directed drive; (N, 3) complex, zero for r < a.  ``beta``
    overrides β_k (negative controls only).  Points on r = a evaluate the
    exterior field (the mask is r < a(1 − 1e-12))."""
    pts = np.asarray(points, dtype=float)
    if beta is None:
        beta = pec_sphere_l1_tm_beta(k, radius)
    e_reg, e_out = pec_sphere_l1_tm_parts(pts, k=k, radius=radius)
    field = e0 * (e_reg + beta * e_out)
    r = np.linalg.norm(pts, axis=1)
    field[r < radius * (1.0 - 1e-12)] = 0.0
    return field


# ---------------------------------------------------------------------------
# IMP-2 step 0: PEC prolate spheroid in a uniform axial field (quasi-static,
# §9 item 136).  Sources: Landau & Lifshitz ECM §4; Stratton §3.27; Osborn,
# Phys. Rev. 67, 351 (1945) — cited, not verified in-repo (Jin has no
# spheroid).  The identities the test asserts are derived here.
#
# Coordinates.  Semi-axes c (along z) > a, e = √(1 − a²/c²), focal half-
# distance f = c e.  With r₁, r₂ the distances to the foci (0, 0, +f) and
# (0, 0, −f): ξ = (r₂ + r₁)/(2f) ≥ 1, η = (r₂ − r₁)/(2f) ∈ [−1, 1], z = fξη;
# the body surface is ξ = ξ₀ = 1/e.
#
# Laplace separates as Σ [A P_n(ξ) + B Q_n(ξ)] P_n(η); the uniform drive
# φ₀ = −E₀ z = −E₀ f P₁(ξ) P₁(η) excites n = 1 only, and the exterior solution
# decaying at infinity is Q₁(ξ) = ξ artanh(1/ξ) − 1 (∼ 1/(3ξ²)).  φ = 0 on
# ξ = ξ₀ fixes
#     φ = −E₀ f η [ξ − ξ₀ Q₁(ξ)/Q₁(ξ₀)].
# Q₁ satisfies (1 − ξ²)y″ − 2ξy′ + 2y = 0 with
#     Q₁′ = artanh(1/ξ) − ξ/(ξ² − 1),    Q₁″ = 2/(ξ² − 1)²,
# and the Wronskian ξQ₁′ − Q₁ = −1/(ξ² − 1).
# Tip field (η = 1, z = fξ): E_tip/E₀ = 1 − ξ₀Q₁′(ξ₀)/Q₁(ξ₀)
#     = (Q₁ − ξ₀Q₁′)/Q₁ = 1/[(ξ₀² − 1)Q₁(ξ₀)] ≡ 1/L_z,
# and with ξ₀ = 1/e, L_z = (ξ₀² − 1)Q₁(ξ₀) = (1 − e²)/e³ (artanh e − e)
#     = (1 − e²) Σ_{k≥0} e^{2k}/(2k + 3)  →  1/3 − 2e²/15  (1/L_z → 3 + 6e²/5).
# Far field: φ_s = φ + E₀z = E₀ f η ξ₀ Q₁(ξ)/Q₁(ξ₀) → E₀ cosθ f³ξ₀/(3Q₁(ξ₀) r²)
#     = p cosθ / r²,   p = E₀ a²c/(3 L_z)   (= E₀a³ for the sphere);
# outside any sphere enclosing the focal segment the P₁ projection
# R²·(3/2)∫ φ_s μ dμ equals p exactly (higher multipoles are orthogonal).
# Cancellation traps: ξ artanh(1/ξ) − 1 and artanh e − e use series branches
# for ξ > 10 and e < 0.05.
# ---------------------------------------------------------------------------

_PROLATE_XI_SERIES = 10.0
_PROLATE_E_SERIES = 0.05
_PROLATE_NTERMS = 40


def prolate_spheroid_depolarisation(e, *, branch=None):
    """Axial depolarisation factor ``L_z(e) = (1 − e²)/e³ (artanh e − e)`` of a
    prolate spheroid of eccentricity ``e``; the transverse factors are
    ``L_x = L_y = (1 − L_z)/2``.  For ``e < 0.05`` the series
    ``(1 − e²) Σ e^{2k}/(2k + 3)`` (no cancellation); ``branch`` ∈ {"closed",
    "series"} forces one form (tests only).  L_z(0) = 1/3; the tip field of a
    PEC spheroid is ``E₀/L_z``."""
    e = float(e)
    use_series = (e < _PROLATE_E_SERIES) if branch is None else (branch == "series")
    if use_series:
        k = np.arange(_PROLATE_NTERMS)
        return float((1.0 - e * e) * np.sum(e ** (2 * k) / (2 * k + 3)))
    return float((1.0 - e * e) / e**3 * (np.arctanh(e) - e))


def prolate_legendre_q1(xi, *, branch=None):
    """``(Q₁(ξ), Q₁′(ξ), Q₁″(ξ))`` for ξ > 1 (arrays).  Closed form via
    ``artanh(1/ξ) = ½ log1p(2/(ξ − 1))``; for ξ > 10 the series
    ``Q₁ = Σ 1/((2k + 3) ξ^{2k+2})``, ``Q₁′ = −Σ (2k + 2)/((2k + 3) ξ^{2k+3})``.
    ``branch`` ∈ {"closed", "series"} forces one form (tests only)."""
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    q = np.empty_like(xi)
    dq = np.empty_like(xi)
    if branch is None:
        ser = xi > _PROLATE_XI_SERIES
    else:
        ser = np.full(xi.shape, branch == "series")
    xm1 = xi - 1.0
    xx1 = xm1 * (xi + 1.0)
    if np.any(~ser):
        x = xi[~ser]
        at = 0.5 * np.log1p(2.0 / xm1[~ser])
        q[~ser] = x * at - 1.0
        dq[~ser] = at - x / xx1[~ser]
    if np.any(ser):
        x = xi[ser]
        k = np.arange(_PROLATE_NTERMS)[:, None]
        inv2 = 1.0 / (x * x)
        pw = inv2[None, :] ** (k + 1)
        q[ser] = np.sum(pw / (2 * k + 3), axis=0)
        dq[ser] = -np.sum((2 * k + 2) * pw / ((2 * k + 3) * x[None, :]), axis=0)
    d2q = 2.0 / (xx1 * xx1)
    return q, dq, d2q


def prolate_spheroidal_coordinates(points, focal):
    """``(ξ, η)`` of Cartesian ``points`` (N, 3) for foci at z = ±``focal``:
    ξ = (r₂ + r₁)/(2f), η = (r₂ − r₁)/(2f), r₁ to (0,0,+f), r₂ to (0,0,−f)."""
    pts = np.asarray(points, dtype=float)
    rho2 = pts[:, 0] ** 2 + pts[:, 1] ** 2
    r1 = np.sqrt(rho2 + (pts[:, 2] - focal) ** 2)
    r2 = np.sqrt(rho2 + (pts[:, 2] + focal) ** 2)
    xi = (r2 + r1) / (2.0 * focal)
    eta = np.clip((r2 - r1) / (2.0 * focal), -1.0, 1.0)
    return np.maximum(xi, 1.0), eta


def pec_prolate_spheroid_potential(
    points, *, semi_major, semi_minor, e0=1.0, scattered_only=False
):
    """Quasi-static potential ``φ = −E₀ f η [ξ − ξ₀ Q₁(ξ)/Q₁(ξ₀)]`` outside a
    grounded PEC prolate spheroid (semi-axis ``semi_major`` = c along z,
    ``semi_minor`` = a) in the uniform field E₀ẑ (φ₀ = −E₀z); zero inside.
    ``scattered_only`` returns φ + E₀z (the induced part).  Dipole moment
    ``E₀ a²c/(3 L_z)``; tip field ``E₀/L_z``."""
    c = float(semi_major)
    a = float(semi_minor)
    if not c > a > 0.0:
        raise ValueError("prolate spheroid needs semi_major > semi_minor > 0")
    e = np.sqrt(1.0 - (a / c) ** 2)
    f = c * e
    xi0 = 1.0 / e
    pts = np.asarray(points, dtype=float)
    xi, eta = prolate_spheroidal_coordinates(pts, f)
    q0 = prolate_legendre_q1(xi0)[0][0]
    q = prolate_legendre_q1(xi)[0]
    inside = (pts[:, 0] ** 2 + pts[:, 1] ** 2) / a**2 + pts[:, 2] ** 2 / c**2 < 1.0 - 1e-12
    phi_s = e0 * f * eta * xi0 * q / q0
    out = phi_s if scattered_only else phi_s - e0 * f * eta * xi
    out = np.where(inside, 0.0, out)
    return out


def pec_prolate_spheroid_field(
    points, *, semi_major, semi_minor, e0=1.0, interior="zero"
):
    """Quasi-static field ``E = E₀ẑ − ∇φ_s`` outside the grounded PEC prolate
    spheroid of :func:`pec_prolate_spheroid_potential` (`IMP-2` step 1, §9
    item 144), as an (N, 3) array.

    With ``φ_s = E₀ f ξ₀ η Q₁(ξ)/Q₁(ξ₀)``,
    ``∇φ_s = E₀ f ξ₀/Q₁(ξ₀) · [Q₁(ξ) ∇η + η Q₁′(ξ) ∇ξ]``, where
    ``∇ξ = (r̂₁ + r̂₂)/(2f)`` and ``∇η = (r̂₂ − r̂₁)/(2f)`` (r̂ᵢ the unit vectors
    from the foci (0, 0, ±f) to the point).  ``interior="zero"`` returns 0
    inside the body (the Dirichlet trace); ``"continue"`` evaluates the
    exterior expression everywhere off the focal segment (quadrature points on
    the chords of a curved cavity wall)."""
    if interior not in ("zero", "continue"):
        raise ValueError(f"interior must be 'zero' or 'continue', got {interior!r}")
    c = float(semi_major)
    a = float(semi_minor)
    if not c > a > 0.0:
        raise ValueError("prolate spheroid needs semi_major > semi_minor > 0")
    e = np.sqrt(1.0 - (a / c) ** 2)
    f = c * e
    xi0 = 1.0 / e
    pts = np.asarray(points, dtype=float)
    d1 = pts - np.array([0.0, 0.0, f])
    d2 = pts + np.array([0.0, 0.0, f])
    r1 = np.linalg.norm(d1, axis=1)
    r2 = np.linalg.norm(d2, axis=1)
    u1 = d1 / r1[:, None]
    u2 = d2 / r2[:, None]
    xi = np.maximum((r2 + r1) / (2.0 * f), 1.0)
    eta = np.clip((r2 - r1) / (2.0 * f), -1.0, 1.0)
    grad_xi = (u1 + u2) / (2.0 * f)
    grad_eta = (u2 - u1) / (2.0 * f)
    q0 = prolate_legendre_q1(xi0)[0][0]
    q, dq, _ = prolate_legendre_q1(xi)
    scale = e0 * f * xi0 / q0
    grad_phi_s = scale * (q[:, None] * grad_eta + (eta * dq)[:, None] * grad_xi)
    out = -grad_phi_s
    out[:, 2] += e0
    if interior == "zero":
        inside = (pts[:, 0] ** 2 + pts[:, 1] ** 2) / a**2 + pts[:, 2] ** 2 / c**2 < 1.0 - 1e-12
        out[inside] = 0.0
    return out


# ---------------------------------------------------------------------------
# IMP-6 step 0: PEC oblate spheroid in a uniform field, axial and in-plane
# (quasi-static, §9 item 403).  Sources (cited, not verified in-repo — Jin has
# no spheroid): Landau & Lifshitz ECM §4; Stratton §3.27; Osborn, Phys. Rev.
# 67, 351 (1945).  The in-repo cross-check is the ellipsoid integral
#     L_i = (a₁a₂a₃/2) ∫₀^∞ ds / ((a_i² + s) √((a₁² + s)(a₂² + s)(a₃² + s))),
# evaluated by quadrature in the test.
#
# Geometry.  Semi-axes a along x and y (``semi_major``), c along z
# (``semi_minor``), a > c.  e = √(1 − c²/a²), focal-ring radius f = a e,
# ξ₀ = c/f = √(1 − e²)/e.  Oblate spheroidal coordinates
#     ρ = f √((ξ² + 1)(1 − η²)),   z = f ξ η,   ξ ≥ 0, η ∈ [−1, 1];
# the body surface is ξ = ξ₀ (f²(ξ₀² + 1) = a², f ξ₀ = c).  ξ of a point is
# the positive root of F = ρ²/(ξ² + 1) + z²/ξ² − f² = 0, i.e. of
#     f² ξ⁴ + (f² − r²) ξ² − z² = 0      (r² = ρ² + z²),
# taken in the cancellation-free form (B = r² − f²):
#     ξ² = (B + √(B² + 4f²z²))/(2f²)   for B ≥ 0,
#     ξ² = 2z² / (√(B² + 4f²z²) − B)   for B < 0.
#
# Separation.  Laplace separates as R(ξ) S(η) e^{imφ} with the radial
# equation of the imaginary-argument Legendre functions.
#  * Axial drive φ₀ = −E₀ z = −E₀ f ξ η (n = 1, m = 0; regular solution ξ).
#    The decaying partner is q₁(ξ) = ξ arccot ξ − 1 (∼ −1/(3ξ²)).  Since
#    η q₁(ξ) = (z/f) h(ξ) with h = q₁/ξ = arccot ξ − 1/ξ,
#        φ = −E₀ z [1 − h(ξ)/h(ξ₀)]        (φ = 0 on ξ = ξ₀).
#  * In-plane drive φ₀ = −E₀ x = −E₀ f √(ξ² + 1) √(1 − η²) cos φ (n = m = 1;
#    regular solution √(1 + ξ²)).  Reduction of order,
#    K = √(1 + ξ²) ∫_ξ^∞ dt / ((1 + t²) · (1 + t²)) · const, gives the
#    decaying K = √(1 + ξ²) [arccot ξ − ξ/(1 + ξ²)] with Wronskian
#    W(√(1+ξ²), K) = −2/(1 + ξ²).  With g = K/√(1 + ξ²),
#        φ = −E₀ x [1 − g(ξ)/g(ξ₀)].
#  Both drives: φ = −E₀ x_i [1 − D_i(ξ)/D_i(ξ₀)], x_i ∈ {z, x}, D ∈ {h, g},
#        h′ = 1/(ξ²(1 + ξ²)),    g′ = −2/(1 + ξ²)².
#
# Field.  E = E₀ x̂_i − ∇φ_s, φ_s = E₀ x_i D(ξ)/D(ξ₀),
#     ∇φ_s = (E₀/D₀) [D(ξ) x̂_i + x_i D′(ξ) ∇ξ],
#     ∇ξ = −∇F/F_ξ,  ∇F = (2x/(ξ² + 1), 2y/(ξ² + 1), 2z/ξ²),
#     F_ξ = −2ξρ²/(ξ² + 1)² − 2z²/ξ³.
#
# Depolarisation factors (each from its own drive's surface field).
#  * Pole (0, 0, c), axial drive: ∂ξ/∂z = ξ₀/c there, D = D₀, so
#        E_z = −E₀ ξ₀ h′(ξ₀)/h₀ = −E₀ / (ξ₀(1 + ξ₀²) h₀) ≡ E₀/L_z,
#        L_z = (1 + ξ₀²)(1 − ξ₀ arccot ξ₀) = (1/e²)[1 − √(1 − e²) arcsin(e)/e]
#    (arccot ξ₀ = arcsin e, 1 + ξ₀² = 1/e²).
#  * Rim (a, 0, 0), in-plane drive: ∂ξ/∂x = (1 + ξ₀²)/(ξ₀ a) there, so
#        E_x = −E₀ a g′(ξ₀)(1 + ξ₀²)/(g₀ ξ₀ a) = 2E₀/(ξ₀(1 + ξ₀²) g₀) ≡ E₀/L_x,
#        L_x = L_y = (ξ₀/2)[(1 + ξ₀²) arccot ξ₀ − ξ₀]
#    — computed from g, independently of L_z, so 2L_x + L_z = 1 is a test.
#  Sphere limit e → 0: L_z = 1/3 + 2e²/15 + …, L_x = 1/3 − e²/15 + ….
#  Disc limit ξ₀ → 0: L_x = πξ₀/4 − ξ₀² + (π/4)ξ₀³ + …, so p_x/a³ → 4/(3π).
#
# Far field.  h ∼ −1/(3ξ³), g ∼ 2/(3ξ³), ξ ∼ r/f, so φ_s → p_i x_i / r³ with
#     p_i = E₀ a²c / (3 L_i)    (= E₀a³ for the sphere);
# outside any sphere enclosing the focal disc the l = 1 projection of φ_s on
# a sphere of radius R equals p_i exactly.
#
# Cancellation traps.  arccot ξ − 1/ξ and arccot ξ − ξ/(1 + ξ²) cancel for
# large ξ; for ξ > 10 the series
#     h = Σ_{k≥1} (−1)^k ξ^{−(2k+1)}/(2k + 1),
#     g = Σ_{k≥1} (−1)^{k+1} (2k/(2k + 1)) ξ^{−(2k+1)}
# (40 terms).  arccot ξ is evaluated as arctan(1/ξ), ξ > 0 only (never inside
# the focal disc).
# ---------------------------------------------------------------------------

_OBLATE_XI_SERIES = 10.0
_OBLATE_NTERMS = 40


def oblate_spheroid_d_functions(xi, *, branch=None):
    """``(h, h′, g, g′)`` of the oblate exterior solutions for ξ > 0 (arrays):
    ``h = arccot ξ − 1/ξ`` (axial drive), ``g = arccot ξ − ξ/(1 + ξ²)``
    (in-plane drive), ``h′ = 1/(ξ²(1 + ξ²))``, ``g′ = −2/(1 + ξ²)²``.  For
    ξ > 10 ``h`` and ``g`` come from their 40-term series; ``branch`` ∈
    {"closed", "series"} forces one form (tests only).  The derivatives have
    no cancellation and are always closed-form."""
    xi = np.atleast_1d(np.asarray(xi, dtype=float))
    h = np.empty_like(xi)
    g = np.empty_like(xi)
    if branch is None:
        ser = xi > _OBLATE_XI_SERIES
    else:
        ser = np.full(xi.shape, branch == "series")
    if np.any(~ser):
        x = xi[~ser]
        acot = np.arctan(1.0 / x)
        h[~ser] = acot - 1.0 / x
        g[~ser] = acot - x / (1.0 + x * x)
    if np.any(ser):
        x = xi[ser]
        k = np.arange(1, _OBLATE_NTERMS + 1)[:, None]
        pw = (1.0 / x)[None, :] ** (2 * k + 1)
        sign = np.where(k % 2 == 0, 1.0, -1.0)
        h[ser] = np.sum(sign * pw / (2 * k + 1), axis=0)
        g[ser] = np.sum(-sign * (2.0 * k / (2 * k + 1)) * pw, axis=0)
    one = 1.0 + xi * xi
    dh = 1.0 / (xi * xi * one)
    dg = -2.0 / (one * one)
    return h, dh, g, dg


def oblate_spheroid_depolarisation(e, *, branch=None):
    """``(L_x, L_z)`` of an oblate spheroid of eccentricity ``e`` (0 < e < 1;
    semi-axes a, a, c with c/a = √(1 − e²)), from ξ₀ = √(1 − e²)/e:
    ``L_z = −ξ₀(1 + ξ₀²) h(ξ₀) = (1 + ξ₀²)(1 − ξ₀ arccot ξ₀)`` and,
    independently, ``L_x = ξ₀(1 + ξ₀²) g(ξ₀)/2``.  ``L_y = L_x``.  The PEC
    rim field under E₀x̂ is ``E₀/L_x`` and the pole field under E₀ẑ is
    ``E₀/L_z``.  ``branch`` is passed to :func:`oblate_spheroid_d_functions`."""
    e = float(e)
    if not 0.0 < e < 1.0:
        raise ValueError("oblate spheroid needs 0 < e < 1")
    xi0 = np.sqrt((1.0 - e) * (1.0 + e)) / e
    return _oblate_l_from_xi0(xi0, branch=branch)


def _oblate_l_from_xi0(xi0, *, branch=None):
    h, _, g, _ = oblate_spheroid_d_functions(xi0, branch=branch)
    w = xi0 * (1.0 + xi0 * xi0)
    return float(0.5 * w * g[0]), float(-w * h[0])


def oblate_spheroidal_xi(points, focal):
    """ξ ≥ 0 of Cartesian ``points`` (N, 3) for the focal ring of radius
    ``focal`` (oblate coordinates ρ = f√((ξ²+1)(1−η²)), z = fξη), from the
    stable root of f²ξ⁴ + (f² − r²)ξ² − z² = 0."""
    pts = np.asarray(points, dtype=float)
    f2 = float(focal) ** 2
    z2 = pts[:, 2] ** 2
    b = np.sum(pts * pts, axis=1) - f2
    disc = np.sqrt(b * b + 4.0 * f2 * z2)
    xi2 = np.empty_like(b)
    pos = b >= 0.0
    xi2[pos] = (b[pos] + disc[pos]) / (2.0 * f2)
    neg = ~pos
    den = disc[neg] - b[neg]
    xi2[neg] = 2.0 * z2[neg] / den
    return np.sqrt(xi2)


def _oblate_setup(semi_major, semi_minor, drive):
    a = float(semi_major)
    c = float(semi_minor)
    if not a > c > 0.0:
        raise ValueError("oblate spheroid needs semi_major > semi_minor > 0")
    if drive not in ("axial", "inplane"):
        raise ValueError(f"drive must be 'axial' or 'inplane', got {drive!r}")
    f = np.sqrt((a - c) * (a + c))
    xi0 = c / f
    return a, c, f, xi0, (2 if drive == "axial" else 0)


def _oblate_inside(pts, a, c):
    return (pts[:, 0] ** 2 + pts[:, 1] ** 2) / a**2 + pts[:, 2] ** 2 / c**2 < 1.0 - 1e-12


def pec_oblate_spheroid_potential(
    points, *, semi_major, semi_minor, drive="axial", e0=1.0, scattered_only=False
):
    """Quasi-static potential ``φ = −E₀ x_i [1 − D(ξ)/D(ξ₀)]`` outside a
    grounded PEC oblate spheroid (equatorial semi-axis ``semi_major`` = a in
    x and y, ``semi_minor`` = c along z) in the uniform field E₀x̂_i:
    ``drive="axial"`` (x_i = z, D = h) or ``"inplane"`` (x_i = x, D = g).
    Zero inside.  ``scattered_only`` returns ``φ_s = φ + E₀x_i``.  Dipole
    moment ``E₀a²c/(3L_i)``."""
    a, c, f, xi0, axis = _oblate_setup(semi_major, semi_minor, drive)
    pts = np.asarray(points, dtype=float)
    xi = oblate_spheroidal_xi(pts, f)
    h, _, g, _ = oblate_spheroid_d_functions(xi)
    h0, _, g0, _ = oblate_spheroid_d_functions(xi0)
    ratio = h / h0[0] if axis == 2 else g / g0[0]
    xi_ = pts[:, axis]
    phi_s = e0 * xi_ * ratio
    out = phi_s if scattered_only else phi_s - e0 * xi_
    return np.where(_oblate_inside(pts, a, c), 0.0, out)


def pec_oblate_spheroid_field(
    points, *, semi_major, semi_minor, drive="axial", e0=1.0, interior="zero"
):
    """Quasi-static field ``E = E₀x̂_i − ∇φ_s`` outside the grounded PEC oblate
    spheroid of :func:`pec_oblate_spheroid_potential` (`IMP-6` step 0, §9
    item 403), as an (N, 3) array, with
    ``∇φ_s = (E₀/D₀)[D(ξ) x̂_i + x_i D′(ξ) ∇ξ]`` and ``∇ξ = −∇F/F_ξ``.
    ``interior="zero"`` returns 0 inside the body (the Dirichlet trace);
    ``"continue"`` evaluates the exterior expression everywhere off the focal
    disc (same contract as :func:`pec_prolate_spheroid_field`)."""
    if interior not in ("zero", "continue"):
        raise ValueError(f"interior must be 'zero' or 'continue', got {interior!r}")
    a, c, f, xi0, axis = _oblate_setup(semi_major, semi_minor, drive)
    pts = np.asarray(points, dtype=float)
    xi = oblate_spheroidal_xi(pts, f)
    h, dh, g, dg = oblate_spheroid_d_functions(xi)
    h0, _, g0, _ = oblate_spheroid_d_functions(xi0)
    if axis == 2:
        d, dd, d0 = h, dh, h0[0]
    else:
        d, dd, d0 = g, dg, g0[0]
    rho2 = pts[:, 0] ** 2 + pts[:, 1] ** 2
    z2 = pts[:, 2] ** 2
    xp1 = xi * xi + 1.0
    grad_f = np.column_stack(
        [2.0 * pts[:, 0] / xp1, 2.0 * pts[:, 1] / xp1, 2.0 * pts[:, 2] / (xi * xi)]
    )
    f_xi = -2.0 * xi * rho2 / (xp1 * xp1) - 2.0 * z2 / xi**3
    grad_xi = -grad_f / f_xi[:, None]
    grad_phi_s = (e0 / d0) * (pts[:, axis] * dd)[:, None] * grad_xi
    grad_phi_s[:, axis] += (e0 / d0) * d
    out = -grad_phi_s
    out[:, axis] += e0
    if interior == "zero":
        out[_oblate_inside(pts, a, c)] = 0.0
    return out


# ---------------------------------------------------------------------------
# `MAT-15` step 2(a): the finite Robin cylinder with a uniform source
# ---------------------------------------------------------------------------


def robin_slab_eigenvalues(biot: float, n: int):
    """First ``n`` roots of ``γ sin γ = Bi cos γ`` (`MAT-15` step 2(a)).

    Axial eigenvalues of the symmetric slab ``|z| ≤ L`` with the Robin
    condition ``−k∂_zθ = hθ`` at ``z = ±L`` (``Bi = hL/k``); the
    eigenfunctions are ``cos(γₙ z/L)``.  For ``Bi > 0`` exactly one root lies
    in ``((n − 1)π, (n − ½)π)`` (``γ tan γ`` runs 0 → ∞ there); found with
    ``brentq``.
    """
    from scipy.optimize import brentq

    if biot <= 0.0:
        raise ValueError(f"biot must be positive, got {biot!r}")

    def f(g):
        return g * np.sin(g) - biot * np.cos(g)

    roots = np.empty(int(n))
    for i in range(1, int(n) + 1):
        lo = (i - 1) * np.pi
        hi = (i - 0.5) * np.pi
        roots[i - 1] = brentq(f, lo, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps,
                              maxiter=500)
    return roots


def robin_cylinder_eigenvalues(biot: float, n: int):
    """First ``n`` roots of ``β J₁(β) = Bi J₀(β)`` (`MAT-15` step 2(a)).

    Radial eigenvalues of the solid cylinder ``r ≤ a`` with ``−k∂_rθ = hθ``
    at ``r = a`` (``Bi = ha/k``); eigenfunctions ``J₀(βₘ r/a)``.  One root lies
    in each ``(j₁,ₘ₋₁, j₀,ₘ)`` (``j₁,₀ = 0``): there ``f = βJ₁ − BiJ₀`` goes
    from ``−Bi J₀(j₁,ₘ₋₁)`` to ``j₀,ₘ J₁(j₀,ₘ)``, of opposite signs.
    """
    from scipy.optimize import brentq
    from scipy.special import j0, j1, jn_zeros

    if biot <= 0.0:
        raise ValueError(f"biot must be positive, got {biot!r}")
    n = int(n)
    z0 = jn_zeros(0, n)
    z1 = np.concatenate([[0.0], jn_zeros(1, n - 1)]) if n > 1 else np.array([0.0])

    def f(b):
        return b * j1(b) - biot * j0(b)

    roots = np.empty(n)
    for i in range(n):
        roots[i] = brentq(f, z1[i], z0[i], xtol=1e-15, rtol=4 * np.finfo(float).eps,
                          maxiter=500)
    return roots


def _robin_cylinder_coefficients(beta):
    """``Cₘ = 2J₁(β)/[β(J₀² + J₁²)]`` — the projection of 1 onto J₀(βr/a)
    with weight r, ∫₀ᵃ J₀(βr/a)² r dr = (a²/2)(J₀² + J₁²) at any β."""
    from scipy.special import j0, j1

    beta = np.asarray(beta, dtype=float)
    a0, a1 = j0(beta), j1(beta)
    return 2.0 * a1 / (beta * (a0**2 + a1**2))


def _robin_slab_coefficients(gamma):
    """``Cₙ = 2 sin γ/(γ + sin γ cos γ)`` — the projection of 1 onto cos(γz/L)."""
    gamma = np.asarray(gamma, dtype=float)
    return 2.0 * np.sin(gamma) / (gamma + np.sin(gamma) * np.cos(gamma))


def robin_finite_cylinder_steady_temperature(r, z, *, radius, half_height, k, q, h,
                                             t_inf=0.0, n_terms=400,
                                             return_gradient=False):
    """Steady temperature of the finite cylinder ``r ≤ a, |z| ≤ L`` with a
    uniform source ``q`` and Robin walls ``−k∂ₙT = h(T − T∞)`` everywhere
    (`MAT-15` step 2(a)).

    Derivation (Carslaw & Jaeger ch. III / VIII style; derived here, not pinned
    to a textbook equation).  Write ``θ = T − T∞ = θ_slab(z) + w``, where
    ``θ_slab = (q/2k)(L² − z²) + qL/h`` solves ``kθ'' + q = 0`` with the Robin
    condition at ``z = ±L``.  Then ``w`` is harmonic, satisfies the homogeneous
    Robin condition at ``z = ±L`` and ``−k∂_r w − h w = h θ_slab`` at ``r = a``.
    With ``γₙ tan γₙ = Bi_z = hL/k`` and ``θ_slab = Σ sₙ cos(γₙz/L)``,
    ``sₙ = (qL²/k)Cₙ/γₙ²``, ``Cₙ = 2 sin γₙ/(γₙ + sin γₙ cos γₙ)``,
    ``w = Σ Dₙ I₀(γₙr/L) cos(γₙz/L)`` with
    ``Dₙ = −h sₙ/[k(γₙ/L) I₁(γₙa/L) + h I₀(γₙa/L)]``.  I₀ ratios are formed
    through ``i0e`` / ``i1e`` × ``e^{γ(r−a)/L}`` (no overflow).

    ``return_gradient`` also returns ``(∂_r T, ∂_z T)``.
    """
    from scipy.special import i0e, i1e

    a = float(radius)
    L = float(half_height)
    r = np.asarray(r, dtype=float)
    z = np.asarray(z, dtype=float)
    gam = robin_slab_eigenvalues(h * L / k, n_terms)
    cn = _robin_slab_coefficients(gam)
    sn = (q * L**2 / k) * cn / gam**2
    xa = gam * a / L
    denom = k * (gam / L) * i1e(xa) + h * i0e(xa)   # scaled by e^{−xa}
    dn_s = -h * sn / denom                           # Dₙ e^{xa}
    rr = r[..., None]
    zz = z[..., None]
    xr = gam * rr / L
    ex = np.exp(xr - xa)
    cz = np.cos(gam * zz / L)
    theta = (q / (2.0 * k)) * (L**2 - z**2) + q * L / h + np.sum(dn_s * i0e(xr) * ex * cz, axis=-1)
    if not return_gradient:
        return t_inf + theta
    dr = np.sum(dn_s * (gam / L) * i1e(xr) * ex * cz, axis=-1)
    dz = -(q / k) * z + np.sum(-dn_s * i0e(xr) * ex * (gam / L) * np.sin(gam * zz / L), axis=-1)
    return t_inf + theta, dr, dz


def robin_finite_cylinder_transient_temperature(r, z, t, *, radius, half_height, k, q, h,
                                                rho_c, t_inf=0.0, t_init=None,
                                                n_radial=200, n_axial=200, n_terms=400,
                                                radial_roots=None, return_gradient=False):
    """Transient temperature of the finite Robin cylinder, uniform source ``q``
    switched on at ``t = 0`` from the uniform ``T(0) = t_init`` (default T∞)
    (`MAT-15` step 2(a)).

    ``θ = θ_ss + Σₘ Σₙ CₘCₙ [θ₀ − q/(kλ²)] e^{−αλ²t} J₀(βₘr/a) cos(γₙz/L)``,
    ``λ² = βₘ²/a² + γₙ²/L²``, ``α = k/ρc``, ``θ₀ = t_init − T∞``: each product
    eigenfunction satisfies both Robin conditions, and at ``t = 0`` the sum
    cancels ``θ_ss`` (whose eigen-expansion has coefficients ``qCₘCₙ/(kλ²)``)
    and leaves θ₀.  ``radial_roots`` overrides the βₘ (negative control only;
    Cₘ is then the same projection formula at those roots).  The sum
    converges slowly at ``t = 0`` — never assert it there.
    """
    from scipy.special import j0, j1

    a = float(radius)
    L = float(half_height)
    if t_init is None:
        t_init = t_inf
    theta0 = float(t_init) - float(t_inf)
    alpha = k / rho_c
    r = np.asarray(r, dtype=float)
    z = np.asarray(z, dtype=float)
    beta = (robin_cylinder_eigenvalues(h * a / k, n_radial) if radial_roots is None
            else np.asarray(radial_roots, dtype=float))
    gam = robin_slab_eigenvalues(h * L / k, n_axial)
    cm = _robin_cylinder_coefficients(beta)
    cn = _robin_slab_coefficients(gam)
    lam2 = (beta[:, None] / a) ** 2 + (gam[None, :] / L) ** 2
    amp = cm[:, None] * cn[None, :] * (theta0 - q / (k * lam2)) * np.exp(-alpha * lam2 * float(t))
    rr = r[..., None]
    zz = z[..., None]
    jr = j0(beta * rr / a)
    cz = np.cos(gam * zz / L)
    s = np.einsum("...m,mn,...n->...", jr, amp, cz)
    ss = robin_finite_cylinder_steady_temperature(
        r, z, radius=a, half_height=L, k=k, q=q, h=h, t_inf=t_inf, n_terms=n_terms,
        return_gradient=return_gradient)
    if not return_gradient:
        return ss + s
    tss, dr_ss, dz_ss = ss
    dr = np.einsum("...m,mn,...n->...", -(beta / a) * j1(beta * rr / a), amp, cz)
    dz = np.einsum("...m,mn,...n->...", jr, amp, -(gam / L) * np.sin(gam * zz / L))
    return tss + s, dr_ss + dr, dz_ss + dz


# ---------------------------------------------------------------------------
# MAT-21 step 0: the rectangular Robin box with a uniform source (closed form)
# ---------------------------------------------------------------------------

_ROBIN_BOX_CHUNK = 2**22   # doubles in the (P × N_y·N_z) intermediate (32 MiB)


def _robin_box_axes(half_widths, k, h, n_terms, roots):
    """Per-axis (L, γ, C) for the box; ``roots[d]`` overrides axis d's γ."""
    if roots is None:
        roots = (None, None, None)
    out = []
    for L, n, rt in zip(half_widths, n_terms, roots):
        L = float(L)
        gam = (robin_slab_eigenvalues(h * L / k, int(n)) if rt is None
               else np.asarray(rt, dtype=float))
        out.append((L, gam, _robin_slab_coefficients(gam)))
    return out


def _robin_box_basis(s, L, gam, deriv):
    """``cos(γ s/L)`` (deriv False) or its s-derivative, shape (P, N)."""
    arg = np.outer(s, gam / L)
    if deriv:
        return -(gam / L)[None, :] * np.sin(arg)
    return np.cos(arg)


def _robin_box_sum(W, x, y, z, axes, deriv=(False, False, False)):
    """``Σ_lmn W_lmn X_l(x) Y_m(y) Z_n(z)`` contracted separably by BLAS.

    ``R = X @ W.reshape(N_x, N_y N_z)`` (P × N_y N_z), times Y and reduced
    over m, then dotted with Z.  Points are processed in chunks so that the
    (P × N_y N_z) intermediate stays ≤ 2²² doubles.
    """
    (Lx, gx, _), (Ly, gy, _), (Lz, gz, _) = axes
    nx, ny, nz = W.shape
    W2 = W.reshape(nx, ny * nz)
    x = np.ravel(x)
    y = np.ravel(y)
    z = np.ravel(z)
    p = x.size
    out = np.empty(p)
    step = max(1, _ROBIN_BOX_CHUNK // max(1, ny * nz))
    for i0 in range(0, p, step):
        sl = slice(i0, min(p, i0 + step))
        X = _robin_box_basis(x[sl], Lx, gx, deriv[0])
        Y = _robin_box_basis(y[sl], Ly, gy, deriv[1])
        Z = _robin_box_basis(z[sl], Lz, gz, deriv[2])
        R = (X @ W2).reshape(-1, ny, nz)
        R = np.einsum("pmn,pm->pn", R, Y, optimize=True)
        out[sl] = np.einsum("pn,pn->p", R, Z)
    return out


def _robin_box_eval(W, x, y, z, axes, t_inf, return_gradient):
    x, y, z = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                                  np.asarray(z, dtype=float))
    shape = x.shape
    T = t_inf + _robin_box_sum(W, x, y, z, axes).reshape(shape)
    if not return_gradient:
        return T
    grads = tuple(
        _robin_box_sum(W, x, y, z, axes, deriv=tuple(i == d for i in range(3))).reshape(shape)
        for d in range(3))
    return T, grads


def _robin_box_lambda(axes):
    (Lx, gx, _), (Ly, gy, _), (Lz, gz, _) = axes
    return ((gx / Lx) ** 2)[:, None, None] + ((gy / Ly) ** 2)[None, :, None] \
        + ((gz / Lz) ** 2)[None, None, :]


def _robin_box_c3(axes):
    (_, _, cx), (_, _, cy), (_, _, cz) = axes
    return cx[:, None, None] * cy[None, :, None] * cz[None, None, :]


def robin_box_steady_temperature(x, y, z, *, half_widths, k, q, h, t_inf=0.0,
                                 n_terms=(40, 40, 40), roots=None, return_gradient=False):
    """Steady temperature of the box ``|x| ≤ A, |y| ≤ B, |z| ≤ C`` with a
    uniform source ``q`` and Robin walls ``−k∂ₙT = h(T − T∞)`` on all six
    faces (`MAT-21` step 0).

    Derivation (Carslaw & Jaeger's product solution of three slab problems;
    a pointer, not a pinned equation — derived here).  For each axis d with
    half-width L_d and ``Bi_d = hL_d/k`` the roots γ of ``γ sin γ = Bi_d cos γ``
    make ``X(s) = cos(γ s/L_d)`` satisfy ``−k X'(±L_d)·(±1) = h X(±L_d)``
    (``kγ/L sin γ = h cos γ``).  ``C = 2 sin γ/(γ + sin γ cos γ)`` projects 1
    onto X (``_robin_slab_coefficients``).  The products X_l Y_m Z_n are
    orthogonal on the box and each satisfies every face's Robin condition, so
    ``1 = Σ C_l C_m C_n X_l Y_m Z_n`` and, with ``∇²(XYZ) = −Λ XYZ``,
    ``Λ_lmn = (γ_l/A)² + (γ_m/B)² + (γ_n/C)²``, ``kΔθ + q = 0`` gives
    ``θ_ss = (q/k) Σ C_l C_m C_n X_l Y_m Z_n / Λ_lmn``.

    Truncated identities (exact for the truncated sum θ_N): Gauss,
    ``h∮θ_N dA = −k∫Δθ_N dV = q∫P_N 1 dV = qV Π_d S_d``,
    ``S_d = Σ_{n≤N_d} C_n sin γ_n/γ_n``.  The truncated sum does not solve
    the PDE pointwise (``kΔθ_N + q = q(1 − Π P_N 1)``).

    ``roots`` is a 3-tuple of ``None`` or a per-axis root array (negative
    control only; that axis takes its length from the array, and C is the same
    projection formula at those roots).  ``return_gradient`` returns
    ``(T, (∂ₓT, ∂ᵧT, ∂_zT))``.
    """
    axes = _robin_box_axes(half_widths, k, h, n_terms, roots)
    W = (q / k) * _robin_box_c3(axes) / _robin_box_lambda(axes)
    return _robin_box_eval(W, x, y, z, axes, float(t_inf), return_gradient)


def robin_box_transient_temperature(x, y, z, t, *, half_widths, k, q, h, rho_c, t_inf=0.0,
                                    t_init=None, n_terms=(40, 40, 40), roots=None,
                                    return_gradient=False):
    """Transient temperature of the Robin box, uniform ``q`` switched on at
    ``t = 0`` from the uniform ``T(0) = t_init`` (default T∞) (`MAT-21` step 0).

    ``θ(t) = Σ C_l C_m C_n X Y Z [(q/(kΛ))(−expm1(−αΛt)) + θ₀ e^{−αΛt}]``,
    ``α = k/ρc``, ``θ₀ = t_init − T∞``: each term is the modal ODE
    ``θ' = −αΛθ + (q/ρc)C_lC_mC_n`` from the projected start.  The term set is
    the steady function's, so ``t → ∞`` returns θ_ss to round-off, and from
    θ₀ = 0 the −expm1 form has no small-t cancellation.  Duhamel (exact at the
    same truncation): ``θ(t) = (q/ρc)∫₀ᵗ U_x U_y U_z ds``,
    ``U_d(s, ·) = Σ C_n X_n e^{−α(γ_n/L_d)² s}``.
    """
    if t_init is None:
        t_init = t_inf
    theta0 = float(t_init) - float(t_inf)
    alpha = k / rho_c
    axes = _robin_box_axes(half_widths, k, h, n_terms, roots)
    lam = _robin_box_lambda(axes)
    at = alpha * lam * float(t)
    W = _robin_box_c3(axes) * ((q / (k * lam)) * (-np.expm1(-at)) + theta0 * np.exp(-at))
    return _robin_box_eval(W, x, y, z, axes, float(t_inf), return_gradient)


def robin_slab_cosine_coefficients(gamma, j):
    """``D⁽ʲ⁾ₙ`` — the projection of ``cos(jπs/L)`` onto ``cos(γₙ s/L)`` on
    ``|s| ≤ L`` (`MAT-21` step 6).

    ``∫cos(as)cos(bs)ds = 2[a sin(aL)cos(bL) − b cos(aL)sin(bL)]/(a² − b²)``
    with ``a = jπ/L``, ``b = γ/L`` gives ``2L(−1)^{j+1} γ sin γ/((jπ)² − γ²)``;
    dividing by ``∫cos²(γs/L)ds = L(γ + sin γ cos γ)/γ``,
    ``D⁽ʲ⁾ = (−1)^{j+1} C γ²/((jπ)² − γ²)`` with C the projection of 1.  For
    Bi > 0, γₙ ∈ ((n − 1)π, (n − ½)π) never equals jπ, so there is no pole.
    ``j = 0`` returns C itself.
    """
    j = int(j)
    c = _robin_slab_coefficients(gamma)
    if j == 0:
        return c
    gamma = np.asarray(gamma, dtype=float)
    sign = -1.0 if j % 2 == 0 else 1.0
    return sign * c * gamma**2 / ((j * np.pi) ** 2 - gamma**2)


def _robin_box_separable_axes(half_widths, k, h, n_terms, profiles):
    axes = _robin_box_axes(half_widths, k, h, n_terms, None)
    axes_p = []
    for (L, gam, c), (a, b, j) in zip(axes, profiles):
        if float(a) == 1.0 and float(b) == 0.0:
            p = c
        else:
            p = float(a) * c + float(b) * robin_slab_cosine_coefficients(gam, j)
        axes_p.append((L, gam, p))
    return axes, axes_p


def robin_box_separable_steady_temperature(x, y, z, *, half_widths, k, q, h, profiles,
                                           t_inf=0.0, n_terms=(40, 40, 40),
                                           return_gradient=False):
    """Steady Robin box (as ``robin_box_steady_temperature``) under the
    separable source ``q·q_x(x)q_y(y)q_z(z)``, ``q_d(s) = a_d + b_d cos(j_d π s/L_d)``
    (`MAT-21` step 6; Carslaw & Jaeger's product solution as a pointer).

    The source projects onto the product basis with per-axis coefficient
    ``P_d = a_d C + b_d D⁽ʲᵈ⁾`` (``robin_slab_cosine_coefficients``), so
    ``θ_ss = (q/k) Σ P_l P_m P_n X_l Y_m Z_n/Λ_lmn``.  Truncated Gauss:
    ``h∮θ_N dA = qV Π_d S′_d``, ``S′_d = Σ P_n sin γₙ/γₙ`` → mean of q_d.
    An axis with ``(1.0, 0.0, ·)`` keeps C itself, so the uniform profile is
    bitwise ``robin_box_steady_temperature`` (same association order).
    """
    axes, axes_p = _robin_box_separable_axes(half_widths, k, h, n_terms, profiles)
    W = (q / k) * _robin_box_c3(axes_p) / _robin_box_lambda(axes)
    return _robin_box_eval(W, x, y, z, axes, float(t_inf), return_gradient)


def robin_box_separable_transient_temperature(x, y, z, t, *, half_widths, k, q, h, rho_c,
                                              profiles, t_inf=0.0, n_terms=(40, 40, 40),
                                              return_gradient=False):
    """Transient separable-source Robin box from θ₀ = 0 (`MAT-21` step 6):
    ``θ = Σ P_l P_m P_n XYZ (q/(kΛ))(−expm1(−αΛt))``; Duhamel
    ``θ = (q/ρc)∫₀ᵗ U′_x U′_y U′_z ds``, ``U′_d = Σ P_n X_n e^{−α(γₙ/L_d)²s}``.
    The uniform profile is bitwise ``robin_box_transient_temperature`` at
    ``t_init = t_inf``."""
    alpha = k / rho_c
    axes, axes_p = _robin_box_separable_axes(half_widths, k, h, n_terms, profiles)
    lam = _robin_box_lambda(axes)
    at = alpha * lam * float(t)
    W = _robin_box_c3(axes_p) * ((q / (k * lam)) * (-np.expm1(-at)))
    return _robin_box_eval(W, x, y, z, axes, float(t_inf), return_gradient)


def robin_slab_quadratic_coefficients(gamma):
    """``Eₙ`` — the projection of ``(s/L)²`` onto ``cos(γₙ s/L)`` on ``|s| ≤ L``
    (`MAT-21` step 9).

    ``∫_{−L}^{L}(s/L)² cos(γs/L) ds = 2L[sin γ/γ + 2 cos γ/γ² − 2 sin γ/γ³]``;
    dividing by ``∫cos²(γs/L)ds = L(γ + sin γ cos γ)/γ`` gives
    ``E = 2[sin γ + 2 cos γ/γ − 2 sin γ/γ²]/(γ + sin γ cos γ)``, valid at any
    γ > 0 (no division by sin γ).  At a Robin root (``cot γ = γ/Bi``) this
    equals ``C(1 + 2/Bi − 2/γ²)``; that identity is a test, not used here.
    """
    g = np.asarray(gamma, dtype=float)
    s = np.sin(g)
    c = np.cos(g)
    return 2.0 * (s + 2.0 * c / g - 2.0 * s / g**2) / (g + s * c)


def _robin_box_quadratic_tensor(half_widths, k, h, n_terms, weights):
    w0, wx, wy, wz = (float(w) for w in weights)
    axes = _robin_box_axes(half_widths, k, h, n_terms, None)
    ax_e = [(L, gam, robin_slab_quadratic_coefficients(gam)) for (L, gam, _) in axes]
    axes_ex = [ax_e[0], axes[1], axes[2]]
    axes_ey = [axes[0], ax_e[1], axes[2]]
    axes_ez = [axes[0], axes[1], ax_e[2]]
    T = (w0 * _robin_box_c3(axes) + wx * _robin_box_c3(axes_ex)
         + wy * _robin_box_c3(axes_ey) + wz * _robin_box_c3(axes_ez))
    return axes, T


def robin_box_quadratic_steady_temperature(x, y, z, *, half_widths, k, q, h, weights,
                                           t_inf=0.0, n_terms=(40, 40, 40),
                                           return_gradient=False):
    """Steady Robin box (as ``robin_box_steady_temperature``) under the
    quadratic source ``q·[w₀ + w_x(x/A)² + w_y(y/B)² + w_z(z/C)²]``,
    ``weights = (w₀, w_x, w_y, w_z)`` (`MAT-21` step 9; Carslaw & Jaeger's
    product solution as a pointer, derived here).

    Each source term is separable: ``(s/L_d)²`` projects onto axis d's basis
    with ``Eₙ`` (``robin_slab_quadratic_coefficients``) and 1 onto the others
    with ``Cₙ``, so the coefficient tensor is
    ``T = w₀ C⊗C⊗C + w_x E⊗C⊗C + w_y C⊗E⊗C + w_z C⊗C⊗E`` and
    ``θ_ss = (q/k) Σ T_lmn X_l Y_m Z_n/Λ_lmn``.  Truncated Gauss:
    ``h∮θ_N dA = qV[w₀ Π S_d + Σ_d w_d S″_d Π_{e≠d} S_e]``,
    ``S″_d = Σ Eₙ sin γₙ/γₙ`` → 1/3 (the mean of (s/L)²).  At
    ``weights = (1, 0, 0, 0)`` this is bitwise ``robin_box_steady_temperature``
    (``1.0·c3 == c3``, ``c3 + 0.0 == c3``, same association of W).
    """
    axes, T = _robin_box_quadratic_tensor(half_widths, k, h, n_terms, weights)
    W = (q / k) * T / _robin_box_lambda(axes)
    return _robin_box_eval(W, x, y, z, axes, float(t_inf), return_gradient)


def robin_box_quadratic_transient_temperature(x, y, z, t, *, half_widths, k, q, h, rho_c,
                                              weights, t_inf=0.0, n_terms=(40, 40, 40),
                                              return_gradient=False):
    """Transient quadratic-source Robin box from θ₀ = 0 (`MAT-21` step 9):
    ``θ = Σ T_lmn XYZ (q/(kΛ))(−expm1(−αΛt))``, ``α = k/ρc``.  Duhamel
    ``θ = (q/ρc)∫₀ᵗ[w₀U_xU_yU_z + w_xU″_xU_yU_z + w_yU_xU″_yU_z + w_zU_xU_yU″_z] ds``,
    ``U″_d = Σ Eₙ Xₙ e^{−α(γₙ/L_d)²s}``.  At unit weights it is bitwise
    ``robin_box_transient_temperature`` at ``t_init = t_inf``."""
    alpha = k / rho_c
    axes, T = _robin_box_quadratic_tensor(half_widths, k, h, n_terms, weights)
    lam = _robin_box_lambda(axes)
    at = alpha * lam * float(t)
    W = T * ((q / (k * lam)) * (-np.expm1(-at)))
    return _robin_box_eval(W, x, y, z, axes, float(t_inf), return_gradient)


def robin_slab_odd_eigenvalues(biot: float, n: int):
    """First ``n`` roots of ``γ cos γ = −Bi sin γ`` (`MAT-21` step 19).

    Odd eigenfunctions ``sin(γₙ s/L)`` of the symmetric slab ``|s| ≤ L`` with
    ``−k∂ₙθ = hθ`` at ``s = ±L`` (``Bi = hL/k``): ``−k(γ/L)cos γ = h sin γ``.
    Carslaw & Jaeger §3.10 (the slab with radiation at both faces) give the
    even and odd families; the odd one is needed for an initial state odd in s.
    ``f = γ cos γ + Bi sin γ`` is ``Bi(−1)^{n+1}`` at ``(n − ½)π`` and
    ``nπ(−1)ⁿ`` at ``nπ`` (opposite signs), and ``tan γ = −γ/Bi < 0`` puts
    exactly one root in each ``((n − ½)π, nπ)``; bracketed by ``brentq``,
    never seeded.
    """
    from scipy.optimize import brentq

    if biot <= 0.0:
        raise ValueError(f"biot must be positive, got {biot!r}")

    def f(g):
        return g * np.cos(g) + biot * np.sin(g)

    roots = np.empty(int(n))
    for i in range(1, int(n) + 1):
        lo = (i - 0.5) * np.pi
        hi = i * np.pi
        roots[i - 1] = brentq(f, lo, hi, xtol=1e-15, rtol=4 * np.finfo(float).eps,
                              maxiter=500)
    return roots


def robin_slab_odd_linear_coefficients(gamma):
    """``bₙ`` — the projection of ``s/L`` onto ``sin(γₙ s/L)`` on ``|s| ≤ L``
    (`MAT-21` step 19).

    ``∫_{−L}^{L} s sin(γs/L) ds = 2L²(sin γ − γ cos γ)/γ²`` and
    ``∫ sin²(γs/L) ds = L(γ − sin γ cos γ)/γ`` give
    ``b = 2(sin γ − γ cos γ)/(γ(γ − sin γ cos γ))``, valid at any γ > 0 (at the
    Dirichlet roots γ = mπ it is the Fourier sine coefficient 2(−1)^{m+1}/(mπ)).
    """
    g = np.asarray(gamma, dtype=float)
    s = np.sin(g)
    c = np.cos(g)
    return 2.0 * (s - g * c) / (g * (g - s * c))


_ROBIN_GRAD_CHUNK = 2**18   # points per chunk of the per-axis (P × N) factors


def _robin_axis_factor(s, L, gam, coef, rate, odd, deriv):
    """``Σₙ coefₙ e^{−rateₙ} φₙ(s)`` (φ = cos or sin of γs/L) or its s-derivative."""
    arg = np.outer(s, gam / L)
    w = coef * np.exp(-rate)
    if odd:
        B = (gam / L)[None, :] * np.cos(arg) if deriv else np.sin(arg)
    else:
        B = -(gam / L)[None, :] * np.sin(arg) if deriv else np.cos(arg)
    return B @ w


def robin_box_initial_gradient_temperature(x, y, z, t, *, half_widths, k, h, rho_c, g,
                                           t_inf=0.0, theta_uniform=0.0, n_terms=40,
                                           roots_z=None, return_gradient=False):
    """Temperature of the source-free Robin box from the non-uniform initial
    state ``T(0) = T∞ + ḡ + g·z`` (``ḡ = theta_uniform``) (`MAT-21` step 19).

    Returns the temperature ``T`` (``T∞`` included, as every Robin-box helper
    does; ``T − T∞`` is the rise θ).  The initial state is uniform in x and y,
    so the x and y factors carry only the projection of 1 onto the even modes
    (``Cₙ``, ``_robin_slab_coefficients``); in z, ``ḡ`` projects onto the even
    modes (``Cₙ``) and ``g·z`` onto the odd modes ``sin(γₙ z/C)``,
    ``γ cos γ = −Bi sin γ`` (``robin_slab_odd_eigenvalues``; Carslaw & Jaeger
    §3.10), with coefficient ``g·C·bₙ`` (``robin_slab_odd_linear_coefficients``).
    Every product mode decays as ``e^{−αΛt}``, ``Λ = Σ_d (γ_d/L_d)²``,
    ``α = k/ρc``, and ``e^{−αΛt}`` is itself a product over axes, so

        ``θ = U_x(x,t) U_y(y,t) [ḡ U_z(z,t) + g C V_z(z,t)]``,
        ``U_d = Σ Cₙ cos(γₙ s/L_d) e^{−α(γₙ/L_d)²t}``,
        ``V_z = Σ bₙ sin(γₙ' z/C) e^{−α(γₙ'/C)²t}``

    — the triple sum factorised exactly.  ``n_terms`` is an int or a per-axis
    3-tuple (the z entry serves both z families).  ``roots_z`` overrides the
    odd z roots (negative control only; ``bₙ`` is the same formula there).
    ``return_gradient`` returns ``(T, (∂ₓT, ∂ᵧT, ∂_zT))``.
    """
    if np.ndim(n_terms) == 0:
        n_terms = (int(n_terms),) * 3
    alpha = k / rho_c
    t = float(t)
    g = float(g)
    gbar = float(theta_uniform)
    axes = _robin_box_axes(half_widths, k, h, n_terms, None)
    Lz = float(half_widths[2])
    gz_odd = (robin_slab_odd_eigenvalues(h * Lz / k, int(n_terms[2])) if roots_z is None
              else np.asarray(roots_z, dtype=float))
    bz = robin_slab_odd_linear_coefficients(gz_odd)
    rate_odd = alpha * (gz_odd / Lz) ** 2 * t
    rates = [alpha * (gam / L) ** 2 * t for (L, gam, _) in axes]

    x, y, z = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float),
                                  np.asarray(z, dtype=float))
    shape = x.shape
    xr, yr, zr = x.ravel(), y.ravel(), z.ravel()
    p = xr.size
    T = np.empty(p)
    G = [np.empty(p) for _ in range(3)] if return_gradient else None
    (Lx, gx, cx), (Ly, gy, cy), (_, gze, cz) = axes
    for i0 in range(0, p, _ROBIN_GRAD_CHUNK):
        sl = slice(i0, min(p, i0 + _ROBIN_GRAD_CHUNK))
        ux = _robin_axis_factor(xr[sl], Lx, gx, cx, rates[0], False, False)
        uy = _robin_axis_factor(yr[sl], Ly, gy, cy, rates[1], False, False)
        wz = gbar * _robin_axis_factor(zr[sl], Lz, gze, cz, rates[2], False, False) \
            + g * Lz * _robin_axis_factor(zr[sl], Lz, gz_odd, bz, rate_odd, True, False)
        T[sl] = ux * uy * wz
        if return_gradient:
            dux = _robin_axis_factor(xr[sl], Lx, gx, cx, rates[0], False, True)
            duy = _robin_axis_factor(yr[sl], Ly, gy, cy, rates[1], False, True)
            dwz = gbar * _robin_axis_factor(zr[sl], Lz, gze, cz, rates[2], False, True) \
                + g * Lz * _robin_axis_factor(zr[sl], Lz, gz_odd, bz, rate_odd, True, True)
            G[0][sl] = dux * uy * wz
            G[1][sl] = ux * duy * wz
            G[2][sl] = ux * uy * dwz
    T = float(t_inf) + T.reshape(shape)
    if not return_gradient:
        return T
    return T, tuple(gd.reshape(shape) for gd in G)


# ---------------------------------------------------------------------------
# MAT-21 step 20: a container wall on the slab, steady, in closed form
# ---------------------------------------------------------------------------


def wall_effective_heat_transfer(h, wall, k_w):
    """``h_eff = (1/h + d/k_w)⁻¹`` — a planar, source-free wall of thickness
    ``d`` and conductivity ``k_w`` in series with the film ``h`` (`MAT-21`
    step 20).  Algebraic: in steady 1D the wall carries a constant flux, so
    its conduction resistance d/k_w adds to the film's 1/h."""
    return 1.0 / (1.0 / float(h) + float(wall) / float(k_w))


def two_layer_slab_steady_temperature(z, *, half_thickness, wall, k, k_w, q, h, t_inf=0.0,
                                      return_gradient=False):
    """Steady temperature of the symmetric 1D two-layer slab (`MAT-21` step 20).

    Gel ``|z| ≤ C`` (conductivity k, uniform source q) inside a source-free
    wall ``C < |z| ≤ C + d`` (conductivity k_w), Robin ``−k_w ∂ₙT = h(T − T∞)``
    at ``|z| = C + d``.  Derived here (no reference equation pinned): the flux
    through ``|z| = s ≤ C`` is ``q s`` by symmetry, so the gel is the parabola
    ``T = T_C + q(C² − z²)/(2k)``; the wall carries the constant flux ``qC``
    and is linear, ``T = T_C − (qC/k_w)(|z| − C)``; the film gives
    ``T(C + d) − T∞ = qC/h``.  Hence ``T_C = T∞ + qC(1/h + d/k_w) = T∞ +
    qC/h_eff`` — on the gel side, exactly the bare Robin slab at ``h_eff``
    (`wall_effective_heat_transfer`).  Temperature and flux are continuous at
    ``|z| = C`` by construction.

    Points with ``|z| ≤ C`` take the gel branch, ``|z| > C`` the wall branch;
    ``return_gradient`` returns ``(T, ∂_zT)``.
    """
    z = np.asarray(z, dtype=float)
    C = float(half_thickness)
    d = float(wall)
    q = float(q)
    t_c = float(t_inf) + q * C * (1.0 / float(h) + d / float(k_w))
    a = np.abs(z)
    gel = a <= C
    T = np.where(gel, t_c + q * (C * C - z * z) / (2.0 * float(k)),
                 t_c - (q * C / float(k_w)) * (a - C))
    if not return_gradient:
        return T
    dT = np.where(gel, -q * z / float(k), -(q * C / float(k_w)) * np.sign(z))
    return T, dT


# ---------------------------------------------------------------------------
# MAT-21 step 20b: the two-layer slab transient (gel + wall), closed form
# ---------------------------------------------------------------------------


def _sinc_x(x):
    """``sin(x)/x`` (1 at 0)."""
    return np.sinc(np.asarray(x, dtype=float) / np.pi)


def _two_layer_slab_mu(x, k, k_w, rho_c, rho_c_w):
    x = np.asarray(x, dtype=float)
    return x * np.sqrt(float(rho_c) / float(k)), x * np.sqrt(float(rho_c_w) / float(k_w))


def _two_layer_slab_F(x, *, half_thickness, wall, k, k_w, rho_c, rho_c_w, h,
                      return_scale=False):
    """Robin residual of the even mode at ``|z| = C + d`` as a function of
    ``x = √λ``: ``F = k_w X′(C + d) + h X(C + d)``.  Entire in λ, ``F(0) = h``."""
    C, d = float(half_thickness), float(wall)
    m1, m2 = _two_layer_slab_mu(x, k, k_w, rho_c, rho_c_w)
    c1, s1 = np.cos(m1 * C), np.sin(m1 * C)
    c2, s2 = np.cos(m2 * d), np.sin(m2 * d)
    sd = d * _sinc_x(m2 * d)                     # sin(μ₂d)/μ₂
    t1 = h * c1 * c2
    t2 = -float(k) * m1 * s1 * c2
    t3 = -float(k_w) * m2 * m2 * c1 * sd
    t4 = -h * (float(k) / float(k_w)) * m1 * s1 * sd
    F = t1 + t2 + t3 + t4
    if return_scale:
        return F, np.abs(t1) + np.abs(t2) + np.abs(t3) + np.abs(t4)
    return F


def two_layer_slab_eigenvalues(*, half_thickness, wall, k, k_w, rho_c, rho_c_w, h, n,
                               refine=64):
    r"""First ``n`` decay rates ``λₙ`` [1/s] of the symmetric two-layer Robin slab
    (`MAT-21` step 20b): even modes of ``(k X′)′ + ρc λ X = 0`` on
    ``0 ≤ z ≤ C + d`` with ``X′(0) = 0``, T and flux continuous at ``z = C``,
    ``k_w X′ + h X = 0`` at ``z = C + d``.

    Bracketing is on the λ axis through the monotone map ``x = √λ`` (both
    layer wavenumbers are ∝ x, so neither μ₁ nor μ₂ sets the ladder): ``F`` is
    scanned at ``π/(refine·Λ)``, ``Λ = C√(ρc/k) + d√(ρ_wc_w/k_w)`` (the
    asymptotic root spacing in x is π/Λ), each sign change refined with
    ``brentq`` and polished to the |F|-minimising double within ±16 ulp —
    the `two_layer_sphere_eigenvalues` precedent in slab geometry.
    """
    from scipy.optimize import brentq

    kw = dict(half_thickness=half_thickness, wall=wall, k=k, k_w=k_w, rho_c=rho_c,
              rho_c_w=rho_c_w, h=h)
    lam_len = (float(half_thickness) * np.sqrt(rho_c / k)
               + float(wall) * np.sqrt(rho_c_w / k_w))
    step = np.pi / (float(refine) * lam_len)
    n = int(n)
    roots = []
    x_lo = 0.0
    while len(roots) < n:
        need = n - len(roots)
        x_hi = x_lo + (need + 8) * 1.25 * np.pi / lam_len
        xs = np.arange(x_lo, x_hi + step, step)
        fs = _two_layer_slab_F(xs, **kw)
        idx = np.nonzero(np.sign(fs[:-1]) * np.sign(fs[1:]) < 0)[0]
        for i in idx:
            x0 = brentq(lambda x: float(_two_layer_slab_F(x, **kw)), xs[i], xs[i + 1],
                        xtol=1e-300, rtol=4 * np.finfo(float).eps, maxiter=500)
            cand = x0 + np.spacing(x0) * np.arange(-16, 17)
            roots.append(float(cand[np.argmin(np.abs(_two_layer_slab_F(cand, **kw)))]))
            if len(roots) == n:
                break
        x_lo = xs[-1]
    x = np.asarray(roots[:n])
    return x * x


def two_layer_slab_mode(z, lam, *, half_thickness, wall, k, k_w, rho_c, rho_c_w, deriv=0):
    """Even mode ``Xₙ(z)`` (columns over ``lam``) or its ``deriv``-th z-derivative
    (0, 1, 2).  Gel ``cos(μ₁z)``; wall, with ``s = |z| − C``,
    ``cos(μ₁C) cos(μ₂s) − (kμ₁/k_w) sin(μ₁C) sin(μ₂s)/μ₂`` (T and flux
    continuous at ``|z| = C`` by construction)."""
    C = float(half_thickness)
    z = np.atleast_1d(np.asarray(z, dtype=float))[:, None]
    lam = np.atleast_1d(np.asarray(lam, dtype=float))[None, :]
    m1, m2 = _two_layer_slab_mu(np.sqrt(lam), k, k_w, rho_c, rho_c_w)
    a = np.abs(z)
    gel = a <= C
    s = np.where(gel, 0.0, a - C)
    A = np.cos(m1 * C)
    Bt = -(float(k) / float(k_w)) * m1 * np.sin(m1 * C)    # X′ at s = 0⁺
    cg, sg = np.cos(m1 * z), np.sin(m1 * z)
    cw, sw = np.cos(m2 * s), np.sin(m2 * s)
    S = s * _sinc_x(m2 * s)                                 # sin(μ₂s)/μ₂
    if deriv == 0:
        return np.where(gel, cg, A * cw + Bt * S)
    if deriv == 1:
        return np.where(gel, -m1 * sg, np.sign(z) * (-A * m2 * m2 * S + Bt * cw))
    if deriv == 2:
        return np.where(gel, -m1 * m1 * cg, -m2 * m2 * (A * cw + Bt * S))
    raise ValueError("deriv must be 0, 1 or 2")


def _sin2_over(x, d):
    """``(2x − sin 2x)/(4x³)`` × d³ with x = μ₂d (Taylor through x⁸ below 0.1)."""
    x = np.asarray(x, dtype=float)
    small = np.abs(x) < 0.1
    xs = np.where(small, 1.0, x)
    big = (2 * xs - np.sin(2 * xs)) / (4 * xs**3)
    ser = 1.0 / 3.0 - x**2 / 15.0 + 2.0 * x**4 / 315.0 - x**6 / 2835.0 + x**8 / 77962.5
    return d**3 * np.where(small, ser, big)


def two_layer_slab_mode_integrals(lam, *, half_thickness, wall, k, k_w, rho_c, rho_c_w):
    """Closed-form half-slab integrals over ``0 ≤ z ≤ C + d``, per mode:
    ``N = ∫ρc X²``, ``I = ∫ρc X``, ``G = ∫_gel X = sin(μ₁C)/μ₁``."""
    C, d = float(half_thickness), float(wall)
    lam = np.asarray(lam, dtype=float)
    m1, m2 = _two_layer_slab_mu(np.sqrt(lam), k, k_w, rho_c, rho_c_w)
    A = np.cos(m1 * C)
    Bt = -(float(k) / float(k_w)) * m1 * np.sin(m1 * C)
    x2 = m2 * d
    G = C * _sinc_x(m1 * C)
    n_gel = rho_c * (0.5 * C + 0.5 * C * _sinc_x(2 * m1 * C))
    # wall: X = A cos(μ₂s) + B̃ sin(μ₂s)/μ₂
    cc = 0.5 * d + 0.5 * d * _sinc_x(2 * x2)                # ∫cos²
    ss = _sin2_over(x2, d)                                  # ∫sin²/μ₂²
    cs = 0.5 * d * d * _sinc_x(x2) ** 2                     # ∫cos·sin/μ₂ = sin²(μ₂d)/(2μ₂²)
    n_wall = rho_c_w * (A * A * cc + Bt * Bt * ss + 2.0 * A * Bt * cs)
    i_wall = rho_c_w * (A * d * _sinc_x(x2) + Bt * 0.5 * d * d * _sinc_x(0.5 * x2) ** 2)
    return n_gel + n_wall, rho_c * G + i_wall, G


def two_layer_slab_transient_temperature(z, t, *, half_thickness, wall, k, k_w, rho_c,
                                         rho_c_w, q, h, t_inf=0.0, n_terms=400,
                                         eigenvalues=None):
    r"""Transient of the symmetric 1D two-layer slab (`MAT-21` step 20b).

    Gel ``|z| < C`` (k, ρc, uniform q switched on at t = 0), source-free wall
    ``C < |z| < C + d`` (k_w, ρ_wc_w), Robin h at ``|z| = C + d``, θ(0) = 0:
    ``θ = θ_ss(z) − Σ aₙ Xₙ(z) e^{−λₙt}`` with θ_ss
    `two_layer_slab_steady_temperature`, Xₙ `two_layer_slab_mode`, λₙ
    `two_layer_slab_eigenvalues`.

    Derivation (here; composite-slab textbook form: citation owed).  The
    operator ``−(1/ρc)(kX′)′`` with the interface and Robin conditions is
    self-adjoint in the ρc-weighted inner product, so the modes are
    ρc-orthogonal and ``aₙ = ∫ρc θ_ss Xₙ / ∫ρc Xₙ²``.  The numerator is taken
    in closed form by Green's identity: ``ρcλₙXₙ = −(kXₙ′)′`` and
    ``−(kθ_ss′)′ = q χ_gel`` with the same homogeneous boundary conditions give
    ``∫ρc θ_ss Xₙ = (q/λₙ) ∫_gel Xₙ = (q/λₙ) sin(μ₁C)/μ₁``; the denominator
    is the closed form of `two_layer_slab_mode_integrals`.  Both integrals are
    over the half-slab 0 ≤ z ≤ C + d.  ``z`` is 1D, ``t`` a scalar.
    """
    kw = dict(half_thickness=half_thickness, wall=wall, k=k, k_w=k_w, rho_c=rho_c,
              rho_c_w=rho_c_w)
    lam = (two_layer_slab_eigenvalues(h=h, n=n_terms, **kw) if eigenvalues is None
           else np.asarray(eigenvalues, dtype=float))
    a = two_layer_slab_transient_coefficients(lam, q=q, **kw)
    z = np.atleast_1d(np.asarray(z, dtype=float))
    t = np.asarray(t, dtype=float)
    ss = two_layer_slab_steady_temperature(z, half_thickness=half_thickness, wall=wall, k=k,
                                           k_w=k_w, q=q, h=h, t_inf=t_inf)
    X = two_layer_slab_mode(z, lam, **kw)
    return ss - X @ (a * np.exp(-lam * t))


def two_layer_slab_transient_coefficients(lam, *, half_thickness, wall, k, k_w, rho_c,
                                          rho_c_w, q):
    """``aₙ = (q/λₙ) Gₙ / Nₙ`` (Green's identity; see the transient docstring)."""
    N, _, G = two_layer_slab_mode_integrals(lam, half_thickness=half_thickness, wall=wall,
                                            k=k, k_w=k_w, rho_c=rho_c, rho_c_w=rho_c_w)
    return float(q) * G / (np.asarray(lam, dtype=float) * N)


def _two_layer_slab_duty_state(t, lam, a, *, period, duty, side="right"):
    """Source state ``s`` and modal amplitudes ``Aₙ`` of the duty-cycled
    two-layer slab at time ``t`` (scalar), from rest.  Phases are
    ``[jP, jP + DP]`` on and ``[jP + DP, (j + 1)P]`` off; each ``Aₙ`` obeys
    ``Ȧₙ = λₙ(aₙ s − Aₙ)`` and is advanced exactly phase by phase,
    ``A ← A e^{−λΔt} + s aₙ (1 − e^{−λΔt})`` (``_one_minus_exp_neg``, no
    cancellation at small ``λΔt``).  At a switch instant ``side="right"``
    takes the phase that starts there (``χ = 1`` while ``0 ≤ t mod P < DP``),
    ``side="left"`` the phase that ends there.  ``A`` is continuous; only ``s``
    differs between the two sides."""
    lam = np.asarray(lam, dtype=float)
    a = np.asarray(a, dtype=float)
    t = float(t)
    P = float(period)
    t_on = float(duty) * P
    A = np.zeros_like(lam)

    def adv(A, dt, on):
        if dt <= 0.0:
            return A
        e = np.exp(-lam * dt)
        om = _one_minus_exp_neg(lam * dt).real
        return A * e + on * a * om

    j = 0
    while True:
        t0 = j * P
        phases = [(t0, t0 + t_on, 1.0)]
        if t_on < P:
            phases.append((t0 + t_on, t0 + P, 0.0))
        for p0, p1, on in phases:
            past = t > p1 or (t == p1 and side == "right")
            if past:
                A = adv(A, p1 - p0, on)
                continue
            return on, adv(A, t - p0, on)
        j += 1


def two_layer_slab_duty_cycle_temperature(z, t, *, half_thickness, wall, k, k_w, rho_c,
                                          rho_c_w, q, h, period, duty, t_inf=0.0,
                                          n_terms=400, eigenvalues=None, side="right"):
    r"""Duty-cycled transient of the symmetric 1D two-layer slab (`MAT-21`
    step 21a, §9 item 687).

    The `two_layer_slab_transient_temperature` geometry (gel |z| < C with
    uniform q, source-free wall C < |z| < C + d, Robin h at |z| = C + d),
    θ(0) = 0, with the source on while ``0 ≤ (t mod P) < D·P``.

    Derivation (the same modal ODE as `MAT-17` step 1, here on the two-layer
    modes).  With θ_ss `two_layer_slab_steady_temperature`, Xₙ
    `two_layer_slab_mode`, λₙ `two_layer_slab_eigenvalues` and
    aₙ `two_layer_slab_transient_coefficients` (θ_ss = Σ aₙXₙ, ρc-weighted
    projection), write θ = Σ Aₙ(t) Xₙ; projecting ρc ∂ₜθ = (kθ′)′ + q s(t) χ_gel
    on the ρc-orthogonal modes gives ``Ȧₙ = λₙ(aₙ s(t) − Aₙ)``, Aₙ(0) = 0,
    advanced exactly phase by phase (`_two_layer_slab_duty_state`; the
    `pec_sphere_implant_duty_cycle_from_modes` pattern).  The slowly
    converging Σ aₙXₙ is replaced by the closed-form θ_ss in the source term:
    ``θ = s(t) θ_ss(z) − Σ (aₙ s(t) − Aₙ) Xₙ(z)``.  Duty 1 reduces to
    `two_layer_slab_transient_temperature`.  ``z`` is 1D, ``t`` a scalar;
    ``side`` selects the left / right limit at a switch instant (the two
    differ by ``Δs·(θ_ss − Σ aₙXₙ)``, the series' truncation residual).
    Raises ``ValueError`` for ``duty ∉ (0, 1]`` or ``period ≤ 0``.
    """
    duty = float(duty)
    if not (0.0 < duty <= 1.0):
        raise ValueError(f"duty must lie in (0, 1], got {duty}")
    if not float(period) > 0.0:
        raise ValueError(f"period must be positive, got {period}")
    if side not in ("left", "right"):
        raise ValueError("side must be 'left' or 'right'")
    kw = dict(half_thickness=half_thickness, wall=wall, k=k, k_w=k_w, rho_c=rho_c,
              rho_c_w=rho_c_w)
    lam = (two_layer_slab_eigenvalues(h=h, n=n_terms, **kw) if eigenvalues is None
           else np.asarray(eigenvalues, dtype=float))
    a = two_layer_slab_transient_coefficients(lam, q=q, **kw)
    z = np.atleast_1d(np.asarray(z, dtype=float))
    s, A = _two_layer_slab_duty_state(t, lam, a, period=period, duty=duty, side=side)
    ss = two_layer_slab_steady_temperature(z, half_thickness=half_thickness, wall=wall, k=k,
                                           k_w=k_w, q=q, h=h, t_inf=0.0)
    X = two_layer_slab_mode(z, lam, **kw)
    return float(t_inf) + s * ss - X @ (a * s - A)


# ---------------------------------------------------------------------------
# MAT-17 step 1: duty-cycled source on the Robin-cooled ball (closed form)
# ---------------------------------------------------------------------------


def _one_minus_exp_neg(z):
    """``1 − e^{−z}`` for complex ``z`` without cancellation at small ``|z|``."""
    z = np.asarray(z, dtype=complex)
    small = np.abs(z) < 0.1
    zs = np.where(small, z, 0.0)
    ser = np.zeros_like(zs)
    term = np.ones_like(zs)
    for n in range(1, 14):
        term = term * (-zs) / n
        ser = ser - term
    return np.where(small, ser, 1.0 - np.exp(-np.where(small, 0.0, z)))


def _robin_sphere_complex_w(r, w, *, radius, k, q, h):
    """Steady excess ``Θ`` of ``−k∇²Θ + wΘ = q`` with Robin cooling, complex ``w``.

    `MAT-11`'s perfused steady state continued to complex ``w``:
    ``Θ = (q/w)[1 − Bi·R·S(r)/Dₙ]`` with ``γ = √(w/k)`` (principal root,
    ``Re γ ≥ 0``), ``S(r) = e^{−γR} sinh(γr)/r = e^{−γ(R−r)}(1 − e^{−2γr})/(2r)``
    (``S(0) = γe^{−γR}``) and ``Dₙ = e^{−γR}[γR cosh γR + (Bi − 1) sinh γR]
    = [γR(1 + e^{−2γR}) + (Bi − 1)(1 − e^{−2γR})]/2``.  The scaled form never
    overflows.  At ``w = iωρc`` it is the response ``Θ e^{iωt}`` to a source
    ``q e^{iωt}`` of ``ρc ∂ₜθ = k∇²θ + q``.  Returns complex ``Θ``.
    """
    r = np.asarray(r, dtype=float)
    R = float(radius)
    w = np.asarray(w, dtype=complex)
    bi = h * R / k
    g = np.sqrt(w / k)
    e2 = np.exp(-2.0 * g * R)
    dn = 0.5 * (g * R * (1.0 + e2) + (bi - 1.0) * _one_minus_exp_neg(2.0 * g * R))
    safe = np.where(r > 0.0, r, 1.0)
    s = np.where(
        r > 0.0,
        np.exp(-g * (R - safe)) * _one_minus_exp_neg(2.0 * g * safe) / (2.0 * safe),
        g * np.exp(-g * R),
    )
    return q / w * (1.0 - bi * R * s / dn)


def robin_sphere_periodic_harmonic(r, omega, *, radius, k, rho_c, q, h):
    """Complex harmonic ``Θ(r; ω)`` at ``w = iωρc`` (convention ``e^{+iωt}``)."""
    return _robin_sphere_complex_w(r, 1j * omega * rho_c, radius=radius, k=k, q=q, h=h)


def _robin_sphere_complex_w_dR(*, radius, k, q, h, w):
    """``Θ′(R)`` of :func:`_robin_sphere_complex_w` (analytic)."""
    R = float(radius)
    w = np.asarray(w, dtype=complex)
    bi = h * R / k
    g = np.sqrt(w / k)
    e2 = np.exp(-2.0 * g * R)
    om = _one_minus_exp_neg(2.0 * g * R)
    dn = 0.5 * (g * R * (1.0 + e2) + (bi - 1.0) * om)
    dS = g * (1.0 + e2) / (2.0 * R) - om / (2.0 * R * R)
    return -q / w * bi * R * dS / dn


def _duty_cycle_modal(
    r, t, *, radius, k, rho_c, q, h, period, duty, n_terms=400, periodic=False,
    lam=None, geometric=True, want="theta",
):
    """Modal duty-cycle excess (see :func:`robin_sphere_duty_cycle_temperature`).

    ``lam`` overrides the eigenvalues (negative control: Dirichlet roots, with
    ``Bₙ`` re-projected by the same formula); ``geometric=False`` drops the
    geometric sum (``g_{j,0} = 0`` every cycle).  ``want='dR'`` returns
    ``∂θ/∂r`` at ``r`` from analytic derivatives.
    """
    R = float(radius)
    P = float(period)
    D = float(duty)
    if not (0.0 < D <= 1.0):
        raise ValueError(f"duty must lie in (0, 1], got {duty!r}")
    if P <= 0.0:
        raise ValueError("period must be positive")
    if lam is None:
        lam = robin_sphere_eigenvalues(h * R / k, n_terms)
    lam = np.asarray(lam, dtype=float)
    _, B = robin_sphere_transient_coefficients(lam, radius=R, k=k, q=q, h=h)
    s = lam**2 * k / (rho_c * R**2)
    om_P = -np.expm1(-s * P)
    om_DP = -np.expm1(-s * D * P)
    G = np.exp(-s * (1.0 - D) * P) * om_DP / om_P

    r = np.asarray(r, dtype=float)
    t = np.asarray(t, dtype=float)
    rr, tt = np.broadcast_arrays(r, t)
    j = np.floor(tt / P)
    tau = tt - j * P
    on = tau < D * P
    rho = rr[..., None] / R
    x = lam * rho
    xs = np.where(x == 0.0, 1.0, x)
    if want == "theta":
        X = np.where(x == 0.0, 1.0, np.sin(xs) / xs)
        ss = robin_sphere_steady_temperature(rr, radius=R, k=k, q=q, h=h, t_inf=0.0)
    else:
        X = np.where(x == 0.0, 0.0, (xs * np.cos(xs) - np.sin(xs)) / (xs * xs)) * lam / R
        ss = -q * rr / (3.0 * k)
    jj = j[..., None]
    if periodic:
        g0 = np.broadcast_to(G, jj.shape[:-1] + G.shape)
    elif geometric:
        g0 = G * (-np.expm1(-jj * s * P))
    else:
        g0 = np.zeros(jj.shape[:-1] + G.shape)
    ta = tau[..., None]
    on_series = np.sum(B * X * (1.0 - g0) * np.exp(-s * ta), axis=-1)
    g_end = 1.0 - (1.0 - g0) * np.exp(-s * D * P)
    off_series = np.sum(B * X * g_end * np.exp(-s * np.maximum(ta - D * P, 0.0)), axis=-1)
    return np.where(on, ss - on_series, off_series)


def robin_sphere_duty_cycle_temperature(
    r,
    t,
    *,
    radius: float,
    k: float,
    rho_c: float,
    q: float,
    h: float,
    t_inf: float,
    period: float,
    duty: float,
    n_terms: int = 400,
    periodic: bool = False,
):
    r"""Robin-cooled ball under a duty-cycled uniform source (`MAT-17` step 1).

    Solves ``ρc ∂ₜT = k∇²T + Q·χ(t)`` in ``r < R``, ``−k∂ᵣT = h(T − T∞)`` on
    ``r = R``, ``T(r, 0) = T∞``, with ``χ = 1`` while ``0 ≤ (t mod P) < DP``
    and 0 otherwise (``periodic=True`` returns the periodic steady state).

    Derivation.  With ``θ = T − T∞`` and `MAT-9`'s eigenpairs
    (:func:`robin_sphere_eigenvalues`, ``Xₙ = sin(λₙρ)/(λₙρ)``) and the
    projection ``θ_ss = Σ BₙXₙ`` (:func:`robin_sphere_transient_coefficients`),
    write ``θ = Σ BₙXₙ gₙ(t)``.  Projecting the equation gives
    ``gₙ′ = sₙ(χ − gₙ)``, ``sₙ = λₙ²k/(ρcR²)``.  In cycle ``j``
    (``t = jP + τ``) from a start value ``g_{j,0}``:
    on, ``g = 1 − (1 − g_{j,0})e^{−sτ}``; off, ``g = g_end e^{−s(τ − DP)}`` with
    ``g_end = 1 − (1 − g_{j,0})e^{−sDP}``.  So
    ``g_{j+1,0} = e^{−s(1−D)P}[1 − (1 − g_{j,0})e^{−sDP}]``, an affine map with
    multiplier ``e^{−sP}`` and fixed point
    ``G = e^{−s(1−D)P}(1 − e^{−sDP})/(1 − e^{−sP})``; from ``g_{0,0} = 0``,
    ``g_{j,0} = G(1 − e^{−jsP})`` (the geometric sum).  In the on phase the
    sum ``Σ BₙXₙ`` (which converges only as ``1/N²``) is replaced by the
    closed-form ``θ_ss`` (:func:`robin_sphere_steady_temperature`), leaving
    ``θ = θ_ss − Σ BₙXₙ(1 − g_{j,0})e^{−sₙτ}``; the off phase keeps the
    decaying factor.  Duty 1 reduces to
    :func:`robin_sphere_transient_temperature`.  Equivalent routes (used as
    the test's anchors): Duhamel superposition of step responses over the
    switch instants (Carslaw & Jaeger, *Conduction of Heat in Solids*, ch. I)
    and the Fourier series built on :func:`robin_sphere_periodic_harmonic`.
    Equation numbers are not pinned (no copy of the texts in the repo); the
    derivation is self-contained.  Never evaluate within a few ``1/s_N`` of a
    switch with the default truncation.  Raises ``ValueError`` for
    ``duty ∉ (0, 1]``.
    """
    return t_inf + _duty_cycle_modal(
        r, t, radius=radius, k=k, rho_c=rho_c, q=q, h=h, period=period,
        duty=duty, n_terms=n_terms, periodic=periodic,
    )


# ---------------------------------------------------------------------------
# MAT-18 step 1 — the two-layer (core + shell) Robin sphere in closed form.
# ---------------------------------------------------------------------------


def _two_layer_sphere_steady(r, *, a, b, k1, k2, q1, q2, h, layer=None, want="theta"):
    """θ = T − T∞ (``want="theta"``) or ∂ᵣθ (``want="dr"``) of the steady
    two-layer sphere; ``layer`` forces the core ("core") or shell ("shell")
    formula regardless of ``r`` (interface checks)."""
    r = np.asarray(r, dtype=float)
    theta_b = (q1 * a**3 + q2 * (b**3 - a**3)) / (3.0 * h * b**2)
    safe = np.where(r > 0.0, r, 1.0)
    if want == "theta":
        shell = (theta_b + (q1 - q2) * a**3 * (1.0 / safe - 1.0 / b) / (3.0 * k2)
                 + q2 * (b**2 - r**2) / (6.0 * k2))
        theta_a = (theta_b + (q1 - q2) * a**3 * (1.0 / a - 1.0 / b) / (3.0 * k2)
                   + q2 * (b**2 - a**2) / (6.0 * k2))
        core = theta_a + q1 * (a**2 - r**2) / (6.0 * k1)
    elif want == "dr":
        shell = -((q1 - q2) * a**3 / (3.0 * k2 * safe**2) + q2 * r / (3.0 * k2))
        core = -q1 * r / (3.0 * k1)
    else:
        raise ValueError(f"want must be 'theta' or 'dr', got {want!r}")
    if layer == "core":
        return core
    if layer == "shell":
        return shell
    return np.where(r <= a, core, shell)


def two_layer_sphere_steady_temperature(r, *, a, b, k1, k2, q1, q2=0.0, h, t_inf):
    r"""Steady temperature of a two-layer Robin-cooled sphere (`MAT-18` step 1).

    Core ``r < a`` (conductivity ``k₁``, uniform source ``q₁``), shell
    ``a < r < b`` (``k₂``, ``q₂``), perfect thermal contact at ``r = a``
    (``θ`` and ``k∂ᵣθ`` continuous) and ``−k₂∂ᵣT = h(T − T∞)`` on ``r = b``.

    Derivation.  In each layer ``−(kᵢ/r²)(r²θ′)′ = qᵢ``.  The heat crossing
    the sphere of radius ``r`` equals the heat produced inside it:
    ``−k r²θ′ = q₁a³/3 + q₂(r³ − a³)/3`` in the shell and ``q₁r³/3`` in the
    core (regularity at 0).  At ``r = b`` the Robin condition gives
    ``θ(b) = [q₁a³ + q₂(b³ − a³)]/(3hb²)``; integrating
    ``θ′ = −(q₁ − q₂)a³/(3k₂r²) − q₂r/(3k₂)`` inward from ``b``,

    .. math::

        θ(r) = θ(b) + \frac{(q_1 − q_2)a^3}{3k_2}\Bigl(\frac1r − \frac1b\Bigr)
               + \frac{q_2(b^2 − r^2)}{6k_2}, \qquad a ≤ r ≤ b,

    and in the core ``θ = θ(a) + q₁(a² − r²)/(6k₁)``.  With ``q₂ = 0`` this is
    the spherical-wall resistance ``(1/r₁ − 1/r₂)/(4πk)`` in series with the
    film ``1/(4πb²h)`` (Incropera et al., ch. 3, radial systems) on top of
    the solid-sphere generation profile; Carslaw & Jaeger, *Conduction of
    Heat in Solids*, 2e, ch. IX (composite sphere).  Equation numbers are not
    pinned (no copy of either text in the repo); the derivation above is
    self-contained.  With ``k₂ = k₁``, ``q₂ = q₁`` it reduces to
    :func:`robin_sphere_steady_temperature` with ``radius = b``.
    """
    return t_inf + _two_layer_sphere_steady(r, a=a, b=b, k1=k1, k2=k2, q1=q1, q2=q2, h=h)


def _two_layer_sphere_shell_start(mu1, a, k1, k2, flux_jump=True):
    """``(u₂(a), u₂′(a))`` for the core mode ``u = sin(μ₁r)/μ₁``."""
    z = mu1 * a
    ua = np.sin(z) / mu1
    du1 = np.cos(z)
    if not flux_jump:  # negative control (c1): u₂′(a) = u₁′(a)
        return ua, du1
    return ua, (k1 * du1 - (k1 - k2) * ua / a) / k2


def _two_layer_sphere_F(x, *, a, b, k1, k2, rho_c1, rho_c2, h, flux_jump=True,
                        want="F"):
    """Eigen-function ``F = k₂u₂′(b) − (k₂/b − h)u₂(b)`` at ``√s = x``.

    ``want="scale"`` returns the rounding scale of ``F``: the sum of the
    magnitudes of its constituent terms *before* cancellation,
    ``k₂(|Aμ₂ sin θ| + |B̄ cos θ|) + |k₂/b − h|(|A cos θ| + |B̄ sin θ/μ₂|)``
    with ``θ = μ₂(b − a)``, ``A = u₂(a)`` and
    ``B̄ = (k₁|cos μ₁a| + |k₁ − k₂||sin μ₁a|/(μ₁a))/k₂ ≥ |u₂′(a)|``.  (`MAT-18`
    W1, 2026-10-03: the post-cancellation ``k₂|u₂′(b)| + |k₂/b − h||u₂(b)|``
    read 7.8e-11 over 400 roots because ``u₂′(b)`` is ≈ 4e-3 at the roots while
    its two constituents are O(1); the roots themselves were at round-off by
    every other check — Nₙ 1.2e-15, Gram 1.7e-14, Robin residual 7e-17.)
    ``want="scale_post"`` keeps that W1 normalisation for the record.  Entire
    in ``s`` (``sin(μr)/μ`` and ``cos(μr)`` are even in ``μ``); ``F(0) = hb``.
    """
    x = np.asarray(x, dtype=float)
    pos = x > 0.0
    mu1 = np.where(pos, x * np.sqrt(rho_c1 / k1), 1.0)
    mu2 = np.where(pos, x * np.sqrt(rho_c2 / k2), 1.0)
    A, B = _two_layer_sphere_shell_start(mu1, a, k1, k2, flux_jump)
    A = np.where(pos, A, a)
    B = np.where(pos, B, 1.0)
    th = np.where(pos, mu2 * (b - a), 0.0)
    sinc_l = np.where(pos, np.sin(th) / mu2, b - a)
    ub = A * np.cos(th) + B * sinc_l
    dub = -A * np.where(pos, mu2, 0.0) * np.sin(th) + B * np.cos(th)
    if want == "scale_post":
        return k2 * np.abs(dub) + abs(k2 / b - h) * np.abs(ub)
    if want == "scale":
        z = mu1 * a
        bbar = np.where(pos, (k1 * np.abs(np.cos(z))
                              + abs(k1 - k2) * np.abs(np.sin(z)) / z) / k2, 1.0)
        if not flux_jump:
            bbar = np.where(pos, np.abs(np.cos(z)), 1.0)
        return (k2 * (np.abs(A) * np.where(pos, mu2, 0.0) * np.abs(np.sin(th))
                      + bbar * np.abs(np.cos(th)))
                + abs(k2 / b - h) * (np.abs(A) * np.abs(np.cos(th))
                                     + bbar * np.abs(sinc_l)))
    return k2 * dub - (k2 / b - h) * ub


def two_layer_sphere_eigenvalues(*, a, b, k1, k2, rho_c1, rho_c2, h, n,
                                 _flux_jump=True):
    r"""First ``n`` decay rates ``sₙ`` [1/s] of the two-layer Robin sphere (`MAT-18`).

    Modes ``ψ = u/r`` of ``−∇·(k∇ψ) = sρcψ``; ``μᵢ = √(s/αᵢ)``, ``αᵢ = kᵢ/ρcᵢ``.
    Core ``u = sin(μ₁r)/μ₁`` (regular at 0; the ``1/μ₁`` normalisation makes
    ``F`` entire in ``s`` with ``F(0) = hb ≠ 0``); shell
    ``u₂ = u₂(a)cos μ₂(r − a) + u₂′(a) sin μ₂(r − a)/μ₂`` with ``u₂(a) = u₁(a)``
    (θ continuous) and ``u₂′(a) = [k₁u₁′(a) − (k₁ − k₂)u₁(a)/a]/k₂`` (from
    ``k₁ψ₁′ = k₂ψ₂′``, ``ψ′ = (u′r − u)/r²``).  The Robin condition
    ``−k₂ψ′(b) = hψ(b)`` times ``b²`` is
    ``F(s) = k₂u₂′(b) − (k₂/b − h)u₂(b) = 0``.

    Composite roots are not one per ``π`` interval, so ``F`` is scanned in
    ``x = √s`` at a step ``π/(32Λ)``, ``Λ = a/√α₁ + (b − a)/√α₂`` (the
    asymptotic root spacing is ``π/Λ``), and each sign change is refined
    with ``brentq``.  Nothing divides by ``F`` and no ``cot`` is used.
    Carslaw & Jaeger 2e, ch. IX (composite sphere); equation numbers not
    pinned.  ``_flux_jump=False`` is the (c1) negative control only.
    """
    from scipy.optimize import brentq

    kw = dict(a=a, b=b, k1=k1, k2=k2, rho_c1=rho_c1, rho_c2=rho_c2, h=h,
              flux_jump=_flux_jump)
    lam_len = a * np.sqrt(rho_c1 / k1) + (b - a) * np.sqrt(rho_c2 / k2)
    step = np.pi / (32.0 * lam_len)
    n = int(n)
    roots = []
    x_lo = 0.0
    while len(roots) < n:
        need = n - len(roots)
        x_hi = x_lo + (need + 8) * 1.25 * np.pi / lam_len
        xs = np.arange(x_lo, x_hi + step, step)
        fs = _two_layer_sphere_F(xs, **kw)
        idx = np.nonzero(np.sign(fs[:-1]) * np.sign(fs[1:]) < 0)[0]
        for i in idx:
            x0 = brentq(lambda x: float(_two_layer_sphere_F(x, **kw)),
                        xs[i], xs[i + 1], xtol=1e-300,
                        rtol=4 * np.finfo(float).eps, maxiter=500)
            # brentq stops at rtol = 4ε; polish to the |F|-minimising double
            # within ±16 ulp (the argument μ₁a reaches ≈ 1e3 at n = 400).
            cand = x0 + np.spacing(x0) * np.arange(-16, 17)
            roots.append(float(cand[np.argmin(np.abs(_two_layer_sphere_F(cand, **kw)))]))
            if len(roots) == n:
                break
        x_lo = xs[-1]
    x = np.asarray(roots[:n])
    return x * x


def _two_layer_sphere_modes(r, s, *, a, b, k1, k2, rho_c1, rho_c2, layer=None,
                            flux_jump=True):
    """``(ψ, ∂ᵣψ)`` of the modes at ``s`` (last axis) and radii ``r``."""
    r = np.asarray(r, dtype=float)[..., None]
    s = np.asarray(s, dtype=float)
    mu1 = np.sqrt(s * rho_c1 / k1)
    mu2 = np.sqrt(s * rho_c2 / k2)
    A, B = _two_layer_sphere_shell_start(mu1, a, k1, k2, flux_jump)
    safe = np.where(r > 0.0, r, 1.0)
    z = mu1 * r
    zs = np.where(z > 0.0, z, 1.0)
    psi_c = np.where(z > 0.0, np.sin(zs) / zs, 1.0)
    dpsi_c = np.where(z > 0.0, mu1 * (zs * np.cos(zs) - np.sin(zs)) / zs**2, 0.0)
    th = mu2 * (r - a)
    u2 = A * np.cos(th) + B * np.sin(th) / mu2
    du2 = -A * mu2 * np.sin(th) + B * np.cos(th)
    psi_s = u2 / safe
    dpsi_s = (du2 * safe - u2) / safe**2
    if layer == "core":
        return psi_c, dpsi_c
    if layer == "shell":
        return psi_s, dpsi_s
    core = r <= a
    return np.where(core, psi_c, psi_s), np.where(core, dpsi_c, dpsi_s)


def _two_layer_sphere_mode_integrals(s, *, a, b, k1, k2, rho_c1, rho_c2,
                                     flux_jump=True):
    """``(Nₙ, I₁, I₂)``: ``Nₙ = ∫₀ᵇρc u² dr``, ``I₁ = ∫₀ᵃ u r dr``,
    ``I₂ = ∫ₐᵇ u r dr`` (so ``∫ρcψ r²dr = ρc₁I₁ + ρc₂I₂``), all analytic."""
    s = np.asarray(s, dtype=float)
    mu1 = np.sqrt(s * rho_c1 / k1)
    mu2 = np.sqrt(s * rho_c2 / k2)
    A, B = _two_layer_sphere_shell_start(mu1, a, k1, k2, flux_jump)
    c = B / mu2
    z = mu1 * a
    T = mu2 * (b - a)
    n_core = rho_c1 / mu1**2 * (a / 2.0 - np.sin(2.0 * z) / (4.0 * mu1))
    n_shell = rho_c2 / mu2 * (A**2 * (T / 2.0 + np.sin(2.0 * T) / 4.0)
                              + A * c * np.sin(T) ** 2
                              + c**2 * (T / 2.0 - np.sin(2.0 * T) / 4.0))
    I1 = (np.sin(z) - z * np.cos(z)) / mu1**3
    int_cos = b * np.sin(T) / mu2 + (np.cos(T) - 1.0) / mu2**2
    int_sin = -b * np.cos(T) / mu2 + a / mu2 + np.sin(T) / mu2**2
    I2 = A * int_cos + c * int_sin
    return n_core + n_shell, I1, I2


def _two_layer_sphere_series(r, t, *, a, b, k1, k2, rho_c1, rho_c2, q1, q2, h,
                             theta_i, n_terms=400, s=None, layer=None,
                             want="theta", flux_jump=True):
    """``θ − θ_ss`` (``want="theta"``) or its ∂ᵣ (``"dr"``), shape bcast(r, t)."""
    mk = dict(a=a, b=b, k1=k1, k2=k2, rho_c1=rho_c1, rho_c2=rho_c2)
    if s is None:
        s = two_layer_sphere_eigenvalues(**mk, h=h, n=n_terms, _flux_jump=flux_jump)
    N, I1, I2 = _two_layer_sphere_mode_integrals(s, **mk, flux_jump=flux_jump)
    cn = (q1 * I1 + q2 * I2) / (s * N)          # θ_ss projection (Green form)
    pn = (rho_c1 * I1 + rho_c2 * I2) / N        # projection of the constant 1
    dn = theta_i * pn - cn
    rr, tt = np.broadcast_arrays(np.asarray(r, dtype=float), np.asarray(t, dtype=float))
    psi, dpsi = _two_layer_sphere_modes(rr, s, **mk, layer=layer, flux_jump=flux_jump)
    f = psi if want == "theta" else dpsi
    return np.sum(dn * np.exp(-s * tt[..., None]) * f, axis=-1)


def two_layer_sphere_transient_temperature(
    r, t, *, a, b, k1, k2, rho_c1, rho_c2, q1, q2=0.0, h, t_inf, t_init,
    n_terms=400,
):
    r"""Transient two-layer Robin sphere with a step source (`MAT-18` step 1).

    Solves ``ρcᵢ∂ₜT = kᵢ∇²T + qᵢ`` in the core ``r < a`` and the shell
    ``a < r < b``, perfect contact at ``a``, ``−k₂∂ᵣT = h(T − T∞)`` at ``b``,
    ``T(r, 0) = T_init`` uniform, sources switched on at ``t = 0``.

    Derivation.  With ``θ = T − T∞``, ``θ_ss`` from
    :func:`two_layer_sphere_steady_temperature` and ``θᵢ = T_init − T∞``,
    ``w = θ − θ_ss`` is source-free with homogeneous interface and Robin
    conditions.  The modes ``ψₙ`` and rates ``sₙ`` of
    :func:`two_layer_sphere_eigenvalues` are orthogonal under the piecewise
    weight ``ρc r²`` with ``Nₙ = ∫ρc ψₙ²r²dr = ∫ρc uₙ²dr`` (closed form per
    layer).  Green's second identity on ``−∇·(k∇·)`` (both ``θ_ss`` and ``ψₙ``
    satisfy the interface and Robin conditions, so all surface terms cancel)
    gives the projection of ``θ_ss``,
    ``cₙ = ∫ρcθ_ssψₙr²dr/Nₙ = ∫qψₙr²dr/(sₙNₙ)``, and the constant projects to
    ``pₙ = ∫ρcψₙr²dr/Nₙ``.  Hence

    .. math::

        θ(r, t) = θ_{ss}(r) + \sum_n (θ_i p_n − c_n)\,ψ_n(r)\,e^{−s_n t}.

    ``∫u r dr`` per layer is elementary (``_two_layer_sphere_mode_integrals``).
    Carslaw & Jaeger 2e, ch. IX (composite sphere); equation numbers not
    pinned (no copy in the repo), the derivation is self-contained.  The
    series converges as ``1/N²`` at ``t = 0``: never evaluate it there.
    With equal properties and ``q₂ = q₁`` it is
    :func:`robin_sphere_transient_temperature` at ``radius = b``.
    """
    ss = _two_layer_sphere_steady(r, a=a, b=b, k1=k1, k2=k2, q1=q1, q2=q2, h=h)
    ser = _two_layer_sphere_series(
        r, t, a=a, b=b, k1=k1, k2=k2, rho_c1=rho_c1, rho_c2=rho_c2, q1=q1,
        q2=q2, h=h, theta_i=t_init - t_inf, n_terms=n_terms,
    )
    return t_inf + ss + ser


# ---------------------------------------------------------------------------
# IMP-4 step 0 (§9 item 202): PEC circular cylinder (infinite, radius a) in a
# lossy host at finite k, single-mode drives, closed form.
#
# Convention e^{+jωt}; H^(2)_n(kρ) is the outgoing radial function (Im k < 0,
# so it decays outward). Expansion in H^(2)_n(kρ) e^{jnφ}: Jin 3e §11.2.1,
# eq (11.56); Bessel's equation Jin (E.1); Wronskian Jin (E.19) written for
# H^(2): J_n H_n^(2)' − J_n' H_n^(2) = −2j/(πz). Series as in Harrington,
# Time-Harmonic EM Fields §5-8. Host μ = μ0, so ε_c = k²/(ω²μ0).
#
# TM_z, n = 0 (axial E):
#   E_z(ρ) = E0 [J0(kρ) − J0(ka) H0(kρ)/H0(ka)]          (E_z(a) = 0, PEC)
#   Faraday ∇×E = −jωμ0 H with E = E_z(ρ) ẑ gives ∇×E = −∂ρE_z φ̂, so
#   H_φ = ∂ρE_z/(jωμ0). At ρ = a, ∂ρE_z = E0 k [J0'H0 − J0 H0']/H0(ka)
#   = E0 k (2j/(πka))/H0(ka) by the Wronskian, so the wire current is
#   I = 2πa H_φ(a) = 4 E0 / (ωμ0 H0^(2)(ka)).
#
# TE_z, n = 1 (transverse E):
#   H_z = H0 f(ρ) cos φ,  f = J1(kρ) − [J1'(ka)/H1'(ka)] H1(kρ)
#   Ampère ∇×H = jωε_c E with ∇×(H_z ẑ) = (∂φH_z/ρ) ρ̂ − ∂ρH_z φ̂:
#     E_ρ = −H0 f sin φ /(jωε_c ρ),  E_φ = −H0 f' cos φ /(jωε_c).
#   With H0 = −2jωε_c E0/k the incident part (J1(kρ) ≈ kρ/2) gives
#   E_ρ = E0 sin φ, E_φ = E0 cos φ, i.e. E0 ŷ on the axis. Note
#   −H0/(jωε_c) = 2E0/k, so ω cancels out of E. PEC: E_φ(a) = 0 ⇔ f'(a) = 0.
#   Static limit (k → 0): f ≈ kρ/2 + k a²/(2ρ), so
#   E_ρ = E0(1 + a²/ρ²) sin φ, E_φ = E0(1 − a²/ρ²) cos φ.
# ---------------------------------------------------------------------------


def _pec_cylinder_tm0_radial(rho, *, k, radius, hankel_kind=2):
    """Return (f, f', f'') of the TM n = 0 radial function, f(a) = 0.

    f = J0(kρ) − J0(ka) H0(kρ)/H0(ka); derivatives are the analytic
    chain-rule ones from scipy's recurrence derivatives (jvp / h2vp).
    ``hankel_kind=1`` is the wrong-convention negative control only.
    """
    from scipy import special

    hank = special.hankel2 if hankel_kind == 2 else special.hankel1
    hankp = special.h2vp if hankel_kind == 2 else special.h1vp
    rho = np.asarray(rho, dtype=float)
    k = complex(k)
    z, za = k * rho, k * radius
    coef = special.jv(0, za) / hank(0, za)
    f = special.jv(0, z) - coef * hank(0, z)
    fp = k * (special.jvp(0, z, 1) - coef * hankp(0, z, 1))
    fpp = k**2 * (special.jvp(0, z, 2) - coef * hankp(0, z, 2))
    return f, fp, fpp


def _pec_cylinder_te1_radial(rho, *, k, radius, condition="neumann"):
    """Return (f, f', f'') of the TE n = 1 radial function.

    ``condition="neumann"``: f'(a) = 0 (PEC, E_φ(a) = 0).
    ``condition="dirichlet"``: f(a) = 0 — the wrong-condition negative
    control only.
    """
    from scipy import special

    rho = np.asarray(rho, dtype=float)
    k = complex(k)
    z, za = k * rho, k * radius
    if condition == "neumann":
        coef = special.jvp(1, za, 1) / special.h2vp(1, za, 1)
    elif condition == "dirichlet":
        coef = special.jv(1, za) / special.hankel2(1, za)
    else:
        raise ValueError(f"unknown condition {condition!r}")
    f = special.jv(1, z) - coef * special.hankel2(1, z)
    fp = k * (special.jvp(1, z, 1) - coef * special.h2vp(1, z, 1))
    fpp = k**2 * (special.jvp(1, z, 2) - coef * special.h2vp(1, z, 2))
    return f, fp, fpp


def _cylinder_polar(points):
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 3:
        raise ValueError(f"points must have shape (N, 3), got {pts.shape}")
    x, y = pts[:, 0], pts[:, 1]
    rho = np.hypot(x, y)
    return pts, rho


def pec_cylinder_tm0_field(points, *, k, radius, e0=1.0, interior="zero",
                           _hankel_kind=2):
    """Total E (complex, (N, 3) Cartesian) of the TM_z n = 0 drive on an
    infinite PEC cylinder of radius ``radius`` along z (derivation above).

    E = E0 [J0(kρ) − J0(ka) H0^(2)(kρ)/H0^(2)(ka)] ẑ, e^{+jωt}.
    ``interior="zero"`` returns 0 for ρ < a; ``"continue"`` evaluates the
    formula there.
    """
    if interior not in ("zero", "continue"):
        raise ValueError(f"interior must be 'zero' or 'continue', got {interior!r}")
    pts, rho = _cylinder_polar(points)
    out = np.zeros((pts.shape[0], 3), dtype=complex)
    f, _, _ = _pec_cylinder_tm0_radial(rho, k=k, radius=radius,
                                       hankel_kind=_hankel_kind)
    out[:, 2] = e0 * f
    if interior == "zero":
        out[rho < radius] = 0.0
    return out


def pec_cylinder_te1_field(points, *, k, radius, omega, e0=1.0,
                           interior="zero", _condition="neumann"):
    """Total E (complex, (N, 3) Cartesian) of the TE_z n = 1 drive on an
    infinite PEC cylinder (derivation above): incident E0 ŷ on the axis.

    H_z = H0 f(ρ) cos φ with H0 = −2jωε_c E0/k, ε_c = k²/(ω²μ0);
    E_ρ = −H0 f sin φ/(jωε_c ρ), E_φ = −H0 f' cos φ/(jωε_c).
    The ρ = 0 axis is guarded against the 1/ρ division and returned as 0;
    callers never sample it.
    """
    if interior not in ("zero", "continue"):
        raise ValueError(f"interior must be 'zero' or 'continue', got {interior!r}")
    mu_0 = 4 * np.pi * 1e-7
    k = complex(k)
    eps_c = k**2 / (omega**2 * mu_0)
    h0 = -2j * omega * eps_c * e0 / k
    pref = -h0 / (1j * omega * eps_c)  # = 2 E0 / k
    pts, rho = _cylinder_polar(points)
    on_axis = rho == 0.0
    rho_safe = np.where(on_axis, 1.0, rho)
    cos_p = np.where(on_axis, 1.0, pts[:, 0] / rho_safe)
    sin_p = np.where(on_axis, 0.0, pts[:, 1] / rho_safe)
    f, fp, _ = _pec_cylinder_te1_radial(rho_safe, k=k, radius=radius,
                                        condition=_condition)
    e_rho = pref * f * sin_p / rho_safe
    e_phi = pref * fp * cos_p
    out = np.zeros((pts.shape[0], 3), dtype=complex)
    out[:, 0] = e_rho * cos_p - e_phi * sin_p
    out[:, 1] = e_rho * sin_p + e_phi * cos_p
    out[on_axis] = 0.0
    if interior == "zero":
        out[rho < radius] = 0.0
    return out


def pec_cylinder_tm0_current(*, k, radius, omega, e0=1.0):
    """Total axial wire current of the TM_z n = 0 drive:
    I = 2πa H_φ(a) = 4 E0 / (ωμ0 H0^(2)(ka)) (Wronskian, derivation above).
    ``radius`` enters only through ka."""
    from scipy import special

    mu_0 = 4 * np.pi * 1e-7
    return 4.0 * e0 / (omega * mu_0 * special.hankel2(0, complex(k) * radius))


# ---------------------------------------------------------------------------
# IMP-5 step 0 (§9 item 210): dielectric-coated PEC circular cylinder
# (infinite, PEC radius a, coating a ≤ ρ < b) in a lossy host at finite k,
# single-mode drives, closed form.
#
# Convention e^{+jωt}; H^(2)_n outgoing in the host (Im k_host < 0). The
# coating (k1 = k_coat) excludes ρ = 0, so it carries J_n / Y_n (regular
# there); the host (k2 = k_host) carries J_n + (coef) H^(2)_n. Both regions
# have μ = μ0, so ε_c = k²/(ω²μ0) per region and ε_c,coat/ε_c,host = k1²/k2²
# =: r. Jin 3e: z-invariant scalar reductions with piecewise εr, μr (1.44),
# (1.45) give the scalar interface conditions —
#   TM (E = E_z ẑ):  E_z and μr⁻¹ ∂ρE_z continuous;
#   TE (H = H_z ẑ):  H_z and ε_c⁻¹ ∂ρH_z continuous —
# which are the independent pair of the field conditions (1.46)–(1.49)
# (tangential E, tangential H; normal D and B follow). PEC: n̂ × E = 0 (1.52).
# Outgoing expansion (11.53)/(11.56); Bessel's equation (E.1); Wronskians
# (E.19): J_n Y_n' − J_n' Y_n = 2/(πz), J_n H_n^(2)' − J_n' H_n^(2) = −2j/(πz).
# Series: the layered-cylinder eigenfunction solution (Harrington, Time-
# Harmonic EM Fields ch. 5 §5-8; Balanis, Adv. Eng. EM ch. 11 — neither held
# in the repo nor re-verified; every coefficient is derived here).
#
# TM_z, n = 0:
#   coating  E_z = E0 c g(ρ),  g = J0(k1ρ)Y0(k1a) − Y0(k1ρ)J0(k1a), g(a) = 0;
#   host     E_z = E0 [J0(k2ρ) + s H0(k2ρ)].
#   At ρ = b (E_z and ∂ρE_z continuous, μ equal):
#     c g(b)  − s H0(k2b)      = J0(k2b)
#     c g'(b) − s k2 H0'(k2b)  = k2 J0'(k2b)
#   H_φ = ∂ρE_z/(jωμ0) (Faraday, as in IMP-4). g'(a) = k1[J0'Y0 − Y0'J0](k1a)
#   = −2/(πa) by (E.19), so the wire current is
#     I = 2πa H_φ(a) = 2πa E0 c (−2/(πa))/(jωμ0) = −4 E0 c/(jωμ0).
#
# TE_z, n = 1 (incident E0 ŷ on the axis, H0 = −2jωε_c,host E0/k2 as IMP-4):
#   coating  H_z = H0 c h(ρ) cos φ, h = J1(k1ρ)Y1'(k1a) − Y1(k1ρ)J1'(k1a),
#            h'(a) = 0 ⇔ E_φ(a) = 0 (PEC);
#   host     H_z = H0 [J1(k2ρ) + t H1(k2ρ)] cos φ.
#   At ρ = b (H_z and ε_c⁻¹∂ρH_z continuous; multiply the second by ε_c,coat):
#     c h(b)  − t H1(k2b)        = J1(k2b)
#     c h'(b) − r t k2 H1'(k2b)  = r k2 J1'(k2b)
#   E = ∇×(H_z ẑ)/(jωε_c) with the region's own ε_c:
#     E_ρ = −H0 F sin φ/(jωε_c ρ),  E_φ = −H0 F' cos φ/(jωε_c),
#   F the region's radial function. −H0/(jωε_c,host) = 2E0/k2 and
#   −H0/(jωε_c,coat) = (2E0/k2)(k2²/k1²): ω cancels out of E.
#
# coat_radius == radius is legal: g(a) = 0, g'(a) = −2/(πa) ≠ 0, h'(a) = 0,
# h(a) = 2/(πk1a) ≠ 0 (both by (E.19)), so the 2 × 2 systems stay
# non-singular and reduce to IMP-4's coefficients s0 = −J0/H0 (k2a),
# t0 = −J1'/H1' (k2a).
#
# Private ``condition="swapped"`` applies the other polarisation's flux
# condition (TM: ε-weighted, second row weighted by r; TE: unweighted) and
# exists for the negative controls only.
# ---------------------------------------------------------------------------


def _coated_cylinder_validate(radius, coat_radius, drive, condition, region=None):
    if not (0.0 < radius <= coat_radius):
        raise ValueError(
            f"need 0 < radius <= coat_radius, got radius={radius!r}, "
            f"coat_radius={coat_radius!r}")
    if drive not in ("tm0", "te1"):
        raise ValueError(f"drive must be 'tm0' or 'te1', got {drive!r}")
    if condition not in ("exact", "swapped"):
        raise ValueError(f"condition must be 'exact' or 'swapped', got {condition!r}")
    if region not in (None, "coat", "host"):
        raise ValueError(f"region must be None, 'coat' or 'host', got {region!r}")


def _coated_cylinder_coefficients(drive, *, k_host, k_coat, radius, coat_radius,
                                  condition="exact"):
    """Return (c, s) (TM n = 0) or (c, t) (TE n = 1): the 2 × 2 interface
    solve at ρ = coat_radius (derivation above)."""
    from scipy import special

    _coated_cylinder_validate(radius, coat_radius, drive, condition)
    k1, k2 = complex(k_coat), complex(k_host)
    r = k1**2 / k2**2
    n = 0 if drive == "tm0" else 1
    za, zb1, zb2 = k1 * radius, k1 * coat_radius, k2 * coat_radius
    if drive == "tm0":
        ja, ya = special.jv(0, za), special.yv(0, za)
        u_b = special.jv(0, zb1) * ya - special.yv(0, zb1) * ja
        up_b = k1 * (special.jvp(0, zb1, 1) * ya - special.yvp(0, zb1, 1) * ja)
        w = 1.0 if condition == "exact" else r
    else:
        jpa, ypa = special.jvp(1, za, 1), special.yvp(1, za, 1)
        u_b = special.jv(1, zb1) * ypa - special.yv(1, zb1) * jpa
        up_b = k1 * (special.jvp(1, zb1, 1) * ypa - special.yvp(1, zb1, 1) * jpa)
        w = r if condition == "exact" else 1.0
    mat = np.array([[u_b, -special.hankel2(n, zb2)],
                    [up_b, -w * k2 * special.h2vp(n, zb2, 1)]], dtype=complex)
    rhs = np.array([special.jv(n, zb2), w * k2 * special.jvp(n, zb2, 1)],
                   dtype=complex)
    c, s = np.linalg.solve(mat, rhs)
    return complex(c), complex(s)


def _coated_cylinder_radial(rho, drive, region, *, k_host, k_coat, radius,
                            coat_radius, condition="exact"):
    """Return (F, F', F'') of the region's radial function: E_z = E0 F (TM
    n = 0) or H_z = H0 F cos φ (TE n = 1). ``region`` is "coat" (c·g or c·h
    with k_coat) or "host" (J_n + coef·H_n^(2) with k_host); derivatives are
    scipy's recurrence derivatives with the analytic chain rule."""
    from scipy import special

    _coated_cylinder_validate(radius, coat_radius, drive, condition, region)
    if region is None:
        raise ValueError("region must be 'coat' or 'host' here")
    c, s = _coated_cylinder_coefficients(drive, k_host=k_host, k_coat=k_coat,
                                         radius=radius, coat_radius=coat_radius,
                                         condition=condition)
    rho = np.asarray(rho, dtype=float)
    n = 0 if drive == "tm0" else 1
    if region == "coat":
        k = complex(k_coat)
        z, za = k * rho, k * radius
        if drive == "tm0":
            wj, wy = special.yv(0, za), special.jv(0, za)
        else:
            wj, wy = special.yvp(1, za, 1), special.jvp(1, za, 1)
        f = c * (special.jv(n, z) * wj - special.yv(n, z) * wy)
        fp = c * k * (special.jvp(n, z, 1) * wj - special.yvp(n, z, 1) * wy)
        fpp = c * k**2 * (special.jvp(n, z, 2) * wj - special.yvp(n, z, 2) * wy)
    else:
        k = complex(k_host)
        z = k * rho
        f = special.jv(n, z) + s * special.hankel2(n, z)
        fp = k * (special.jvp(n, z, 1) + s * special.h2vp(n, z, 1))
        fpp = k**2 * (special.jvp(n, z, 2) + s * special.h2vp(n, z, 2))
    return f, fp, fpp


def _coated_region_masks(rho, radius, coat_radius, interior, region):
    if interior not in ("zero", "continue"):
        raise ValueError(f"interior must be 'zero' or 'continue', got {interior!r}")
    if region == "coat":
        return np.ones(rho.shape, bool), np.zeros(rho.shape, bool), None
    if region == "host":
        return np.zeros(rho.shape, bool), np.ones(rho.shape, bool), None
    host = rho >= coat_radius
    zero = (rho < radius) if interior == "zero" else None
    return ~host, host, zero


def coated_pec_cylinder_tm0_field(points, *, k_host, k_coat, radius, coat_radius,
                                  e0=1.0, interior="zero", region=None,
                                  _condition="exact"):
    """Total E (complex, (N, 3) Cartesian) of the TM_z n = 0 drive on an
    infinite PEC cylinder of radius ``radius`` with a dielectric coating to
    ``coat_radius`` in a host of wavenumber ``k_host`` (derivation above).

    ``region=None`` selects by ρ: ρ < radius → ``interior`` ("zero" returns
    0, "continue" evaluates the coating formula), radius ≤ ρ < coat_radius
    → coating, ρ ≥ coat_radius → host. ``region="coat"`` / ``"host"`` forces
    that region's formula at every point.
    """
    _coated_cylinder_validate(radius, coat_radius, "tm0", _condition, region)
    pts, rho = _cylinder_polar(points)
    coat, host, zero = _coated_region_masks(rho, radius, coat_radius, interior,
                                            region)
    out = np.zeros((pts.shape[0], 3), dtype=complex)
    kw = dict(k_host=k_host, k_coat=k_coat, radius=radius,
              coat_radius=coat_radius, condition=_condition)
    for mask, reg in ((coat, "coat"), (host, "host")):
        if np.any(mask):
            f, _, _ = _coated_cylinder_radial(rho[mask], "tm0", reg, **kw)
            out[mask, 2] = e0 * f
    if zero is not None:
        out[zero] = 0.0
    return out


def coated_pec_cylinder_te1_field(points, *, k_host, k_coat, radius, coat_radius,
                                  omega, e0=1.0, interior="zero", region=None,
                                  _condition="exact"):
    """Total E (complex, (N, 3) Cartesian) of the TE_z n = 1 drive on the
    coated PEC cylinder (derivation above): incident E0 ŷ on the axis.

    Region selection as :func:`coated_pec_cylinder_tm0_field`. E uses the
    region's own ε_c = k²/(ω²μ0); ω cancels. The ρ = 0 axis is guarded
    against the 1/ρ division and returned as 0; callers never sample it.
    """
    _coated_cylinder_validate(radius, coat_radius, "te1", _condition, region)
    mu_0 = 4 * np.pi * 1e-7
    k1, k2 = complex(k_coat), complex(k_host)
    eps_host = k2**2 / (omega**2 * mu_0)
    eps_coat = k1**2 / (omega**2 * mu_0)
    h0 = -2j * omega * eps_host * e0 / k2
    pts, rho = _cylinder_polar(points)
    on_axis = rho == 0.0
    rho_safe = np.where(on_axis, 1.0, rho)
    cos_p = np.where(on_axis, 1.0, pts[:, 0] / rho_safe)
    sin_p = np.where(on_axis, 0.0, pts[:, 1] / rho_safe)
    coat, host, zero = _coated_region_masks(rho, radius, coat_radius, interior,
                                            region)
    e_rho = np.zeros(pts.shape[0], dtype=complex)
    e_phi = np.zeros(pts.shape[0], dtype=complex)
    kw = dict(k_host=k_host, k_coat=k_coat, radius=radius,
              coat_radius=coat_radius, condition=_condition)
    for mask, reg, eps in ((coat, "coat", eps_coat), (host, "host", eps_host)):
        if np.any(mask):
            pref = -h0 / (1j * omega * eps)
            f, fp, _ = _coated_cylinder_radial(rho_safe[mask], "te1", reg, **kw)
            e_rho[mask] = pref * f * sin_p[mask] / rho_safe[mask]
            e_phi[mask] = pref * fp * cos_p[mask]
    out = np.zeros((pts.shape[0], 3), dtype=complex)
    out[:, 0] = e_rho * cos_p - e_phi * sin_p
    out[:, 1] = e_rho * sin_p + e_phi * cos_p
    out[on_axis] = 0.0
    if zero is not None:
        out[zero] = 0.0
    return out


def coated_pec_cylinder_tm0_current(*, k_host, k_coat, radius, coat_radius, omega,
                                    e0=1.0):
    """Total axial wire current of the coated TM_z n = 0 drive:
    I = 2πa H_φ(a) = −4 E0 c/(jωμ0) (Wronskian, derivation above)."""
    mu_0 = 4 * np.pi * 1e-7
    c, _ = _coated_cylinder_coefficients("tm0", k_host=k_host, k_coat=k_coat,
                                         radius=radius, coat_radius=coat_radius)
    return -4.0 * e0 * c / (1j * omega * mu_0)


# ---------------------------------------------------------------------------
# IMP-7 step 0 (§9 item 579): steady heating around a PEC-sphere implant by its
# own quasi-static dipole SAR in a Robin-cooled gel shell.  Numpy only.
# ---------------------------------------------------------------------------


def pec_sphere_implant_source(points, *, radius, q_far, beta=1.0):
    r"""Heat source ``q = q_far |E|²/E₀²`` of the quasi-static dipole field
    around a sphere of radius ``a`` (`IMP-7` step 0).

    Outside the sphere ``E = E₀ẑ + βE₀a³(3(ẑ·r̂)r̂ − ẑ)/r³`` (``β = 1`` for
    PEC, ``(ε − 1)/(ε + 2)`` for a dielectric; `IMP-1`, `TH-15`).  With
    ``μ = z/r`` and ``u = (a/r)³``, ``E_r = E₀μ(1 + 2βu)``,
    ``E_ϑ = −E₀ sinϑ (1 − βu)``, so

    .. math::

        |E|^2/E_0^2 = (1 + 2βu)^2μ^2 + (1 − βu)^2(1 − μ^2)
                    = 1 + 2β^2u^2 + (4βu + 2β^2u^2)P_2(μ),

    using ``μ² = (1 + 2P₂)/3`` (Jin 3e App. E §E.3, eqs (E.34)–(E.36)).  At
    ``β = 1`` it is 9 at the poles of ``r = a`` and 0 on its equator.
    ``points`` has shape ``(N, 3)``; the source is zero for ``r < a`` (the
    metal carries no SAR).
    """
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    r2 = np.sum(pts * pts, axis=1)
    r_safe = np.sqrt(np.where(r2 > 0.0, r2, 1.0))
    u = (radius / r_safe) ** 3
    mu = pts[:, 2] / r_safe
    p2 = 0.5 * (3.0 * mu * mu - 1.0)
    q = q_far * (1.0 + 2.0 * beta**2 * u**2
                 + (4.0 * beta * u + 2.0 * beta**2 * u**2) * p2)
    return np.where(r2 < radius * radius, 0.0, q)


def _pec_sphere_implant_terms(coeffs, l):
    """``[(c, n)]`` with ``R_l(r) = Σ c rⁿ`` in the shell (power laws)."""
    a, k2 = coeffs["radius"], coeffs["k_gel"]
    q, beta = coeffs["q_far"], coeffs["beta"]
    if l == 0:
        return [(coeffs["A0"], 0), (coeffs["B0"], -1),
                (-q / (6.0 * k2), 2), (-beta**2 * q * a**6 / (6.0 * k2), -4)]
    if l == 2:
        return [(coeffs["D"], 2), (coeffs["E"], -3),
                (coeffs["C1"], -1), (coeffs["C2"], -4)]
    raise ValueError(f"l must be 0 or 2, got {l!r}")


def pec_sphere_implant_coefficients(*, radius, outer_radius, k_core, k_gel, q_far,
                                    h, beta=1.0):
    r"""Coefficients of the exact steady temperature ``θ = T − T∞ =
    R₀(r) + R₂(r)P₂(cos ϑ)`` around a sphere implant (`IMP-7` step 0).

    Geometry: core ``r < a`` (conductivity ``k₁ = k_core``, no source), gel
    shell ``a < r < b`` (``k₂ = k_gel``, source
    :func:`pec_sphere_implant_source`), perfect contact at ``r = a``,
    ``−k₂∂ᵣθ = hθ`` on ``r = b``.  Carslaw & Jaeger, *Conduction of Heat in
    Solids* 2e, ch. IX (the sphere) for the method; equation numbers not
    pinned, the derivation below is self-contained.

    In the shell ``−k₂Δθ = q``.  For a single harmonic ``Rₗ(r)Pₗ`` the radial
    operator acts on powers as ``Δₗ rⁿ = [n(n + 1) − l(l + 1)] rⁿ⁻²``.

    *l = 0.*  ``q₀ = q_far(1 + 2β²a⁶r⁻⁶)``.  Particular parts: ``−q_far r²/(6k₂)``
    (``n = 2``: 6) and ``−β²q_far a⁶/(6k₂r⁴)`` (``n = −4``: 12), so
    ``R₀ = A₀ + B₀/r − q_far r²/(6k₂) − β²q_far a⁶/(6k₂r⁴)``.  The core has no
    source and is regular, so its ``l = 0`` part is a constant and carries no
    flux; hence ``R₀′(a⁺) = 0`` for every ``k₁`` (including 0 and ∞):
    ``B₀ = (2β² − 1)q_far a³/(3k₂)``.  All the heat leaves through ``r = b``:
    ``P = (4πq_far/3)(b³ − a³ + 2β²a³ − 2β²a⁶/b³) = 4πb²hθ₀(b)``, so
    ``θ₀(b) = q_far(b³ − a³ + 2β²a³ − 2β²a⁶/b³)/(3hb²)`` and
    ``A₀ = θ₀(b) − B₀/b + q_far b²/(6k₂) + β²q_far a⁶/(6k₂b⁴)``.  The core
    constant is ``θ₀_core = R₀(a⁺) = A₀ + B₀/a − (1 + β²)q_far a²/(6k₂)``
    (equal to ``A₀`` only at ``β = ±1``).

    *l = 2.*  ``q₂ = q_far(4βa³r⁻³ + 2β²a⁶r⁻⁶)``.  Particular parts
    ``C₁/r`` (``n = −1``: −6) and ``C₂/r⁴`` (``n = −4``: 6):
    ``C₁ = 2βq_far a³/(3k₂)``, ``C₂ = −β²q_far a⁶/(3k₂)``.  Homogeneous parts
    ``Dr² + Er⁻³``.  The core is ``F r²`` (regular).  Conditions:
    ``R₂(a⁻) = R₂(a⁺)``, ``k₁R₂′(a⁻) = k₂R₂′(a⁺)`` and
    ``k₂R₂′(b) + hR₂(b) = 0``.  Eliminating ``F = R₂(a⁺)/a²`` turns the first
    two into ``k₂R₂′(a⁺) − (2k₁/a)R₂(a⁺) = 0``; with the Robin row this is a
    2 × 2 system for ``(D, E)``, solved by Cramer's rule.  Exact branches:
    ``k_core = 0`` (adiabatic hole: ``R₂′(a⁺) = 0``, core excluded, ``F`` NaN)
    and ``k_core = np.inf`` (isothermal core: ``R₂(a⁺) = 0``, ``F = 0``).

    Returns a dict with ``A0, B0, C1, C2, D, E, F, theta0_core, theta0_b,
    power`` and the inputs ``radius, outer_radius, q_far, k_gel, k_core, h,
    beta`` (read by :func:`pec_sphere_implant_radial`).
    """
    a, b = float(radius), float(outer_radius)
    k2, q, hh, be = float(k_gel), float(q_far), float(h), float(beta)
    k1 = float(k_core)
    if not (k1 >= 0.0):
        raise ValueError(f"k_core must be >= 0 (or np.inf), got {k_core!r}")
    shape = b**3 - a**3 + 2.0 * be**2 * a**3 - 2.0 * be**2 * a**6 / b**3
    B0 = (2.0 * be**2 - 1.0) * q * a**3 / (3.0 * k2)
    theta0_b = q * shape / (3.0 * hh * b**2)
    A0 = theta0_b - B0 / b + q * b**2 / (6.0 * k2) + be**2 * q * a**6 / (6.0 * k2 * b**4)
    theta0_core = A0 + B0 / a - (1.0 + be**2) * q * a**2 / (6.0 * k2)
    C1 = 2.0 * be * q * a**3 / (3.0 * k2)
    C2 = -(be**2) * q * a**6 / (3.0 * k2)

    def p(r):
        return C1 / r + C2 / r**4

    def dp(r):
        return -C1 / r**2 - 4.0 * C2 / r**5

    a21, a22 = 2.0 * k2 * b + hh * b**2, -3.0 * k2 / b**4 + hh / b**3
    r2 = -(k2 * dp(b) + hh * p(b))
    if k1 == 0.0:
        a11, a12, r1 = 2.0 * a, -3.0 / a**4, -dp(a)
    elif np.isinf(k1):
        a11, a12, r1 = a**2, a**-3, -p(a)
    else:
        a11 = 2.0 * k2 * a - 2.0 * k1 * a
        a12 = -3.0 * k2 / a**4 - 2.0 * k1 / a**4
        r1 = -(k2 * dp(a) - (2.0 * k1 / a) * p(a))
    det = a11 * a22 - a12 * a21
    D = (r1 * a22 - a12 * r2) / det
    E = (a11 * r2 - r1 * a21) / det
    r2a = D * a**2 + E / a**3 + p(a)
    if k1 == 0.0:
        F = float("nan")
    elif np.isinf(k1):
        F = 0.0
    else:
        F = r2a / a**2
    return {
        "radius": a, "outer_radius": b, "q_far": q, "k_gel": k2, "k_core": k1,
        "h": hh, "beta": be,
        "A0": A0, "B0": B0, "C1": C1, "C2": C2, "D": D, "E": E, "F": F,
        "theta0_core": theta0_core, "theta0_b": theta0_b,
        "power": 4.0 * np.pi * q * shape / 3.0,
    }


def pec_sphere_implant_radial(r, coeffs, l, *, derivative=0):
    """Shell radial function ``R_l`` (``l`` ∈ {0, 2}) or its first / second
    derivative, as pure power laws from
    :func:`pec_sphere_implant_coefficients`.  No r-mask: it continues the
    shell formula to any ``r > 0`` (`IMP-7` step 1 evaluates it slightly
    inside ``r = a``)."""
    if derivative not in (0, 1, 2):
        raise ValueError(f"derivative must be 0, 1 or 2, got {derivative!r}")
    r = np.asarray(r, dtype=float)
    out = np.zeros_like(r)
    for c, n in _pec_sphere_implant_terms(coeffs, l):
        fac = 1.0 if derivative == 0 else (n if derivative == 1 else n * (n - 1))
        if fac != 0:
            out = out + c * fac * r ** (n - derivative)
    return out


def pec_sphere_implant_steady_temperature(points, *, radius, outer_radius, k_core,
                                          k_gel, q_far, h, t_inf, beta=1.0):
    """Exact steady temperature ``T`` at ``points`` (shape ``(N, 3)``) around a
    sphere implant heated by its own dipole SAR (`IMP-7` step 0; derivation in
    :func:`pec_sphere_implant_coefficients`).  Shell points (``r ≥ a``):
    ``T∞ + R₀(r) + R₂(r)P₂(z/r)``.  Core points: ``T∞ + θ₀_core + F r²P₂``,
    or NaN on the adiabatic branch (``k_core = 0``: the core is not part of
    the domain)."""
    c = pec_sphere_implant_coefficients(radius=radius, outer_radius=outer_radius,
                                        k_core=k_core, k_gel=k_gel, q_far=q_far,
                                        h=h, beta=beta)
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    r2 = np.sum(pts * pts, axis=1)
    shell = r2 >= radius * radius
    r2_s = np.where(shell, r2, radius * radius)
    r_s = np.sqrt(r2_s)
    r2p2 = 0.5 * (3.0 * pts[:, 2] ** 2 - r2)  # r² P₂(μ), regular at r = 0
    th_shell = (pec_sphere_implant_radial(r_s, c, 0)
                + pec_sphere_implant_radial(r_s, c, 2) * r2p2 / r2_s)
    if c["k_core"] == 0.0:
        th_core = np.full_like(r2, np.nan)
    else:
        th_core = c["theta0_core"] + c["F"] * r2p2
    return t_inf + np.where(shell, th_shell, th_core)


# ---------------------------------------------------------------------------
# IMP-7 step 4 (§9 item 658, §10 I45): the transient of the step 0 problem on
# the adiabatic hole (Neumann at r = a, Robin at r = b), started from T∞.
# Hollow-sphere eigen-expansion, Carslaw & Jaeger 2e ch. IX; spherical Bessel
# identities, Jin 3e App. E.
# ---------------------------------------------------------------------------


def _imp7_sph_bessel(l, x, derivative=0):
    """``(j_l(x), y_l(x))`` or their first / second x-derivatives.  The second
    derivative comes from the recurrence ``f_l′ = f_{l−1} − (l + 1)f_l/x``
    differentiated once (``f_0″ = −f_1′``), not from the Bessel ODE, so the
    (PDE) identity in `IMP-7` step 4 is not circular."""
    from scipy.special import spherical_jn, spherical_yn

    x = np.asarray(x, dtype=float)
    if derivative in (0, 1):
        d = bool(derivative)
        return spherical_jn(l, x, derivative=d), spherical_yn(l, x, derivative=d)
    if derivative != 2:
        raise ValueError(f"derivative must be 0, 1 or 2, got {derivative!r}")
    if l == 0:
        return (-spherical_jn(1, x, derivative=True),
                -spherical_yn(1, x, derivative=True))
    out = []
    for f in (spherical_jn, spherical_yn):
        out.append(f(l - 1, x, derivative=True)
                   - (l + 1) / x * f(l, x, derivative=True)
                   + (l + 1) / x**2 * f(l, x))
    return tuple(out)


def pec_sphere_implant_mode(r, lam, l, *, radius, inner="neumann", derivative=0):
    r"""Radial eigenfunction of the hollow sphere (`IMP-7` step 4).

    ``inner="neumann"``: ``R(r) = j_l(λr)y_l′(λa) − y_l(λr)j_l′(λa)``
    (``R′(a) = 0``).  ``inner="dirichlet"``: ``R(r) = j_l(λr)y_l(λa) −
    y_l(λr)j_l(λa)`` (``R(a) = 0``; the step 4 control).  ``derivative`` ∈
    {0, 1, 2} returns the r-derivative.  Unnormalised: callers scale by the
    mode's maximum on (a, b) before any residual (y_l blows up as λa → 0)."""
    if inner not in ("neumann", "dirichlet"):
        raise ValueError(f"inner must be 'neumann' or 'dirichlet', got {inner!r}")
    ja, ya = _imp7_sph_bessel(l, lam * radius, 1 if inner == "neumann" else 0)
    jr, yr = _imp7_sph_bessel(l, lam * np.asarray(r, dtype=float), derivative)
    return lam**derivative * (jr * ya - yr * ja)


def _imp7_robin_det(lam, l, radius, outer_radius, k, h, inner):
    """``k R′(b) + h R(b)`` divided by the (positive) size of the a-side
    coefficients — same roots, tame magnitude."""
    ja, ya = _imp7_sph_bessel(l, lam * radius, 1 if inner == "neumann" else 0)
    scale = np.hypot(ja, ya)
    return (k * pec_sphere_implant_mode(outer_radius, lam, l, radius=radius,
                                        inner=inner, derivative=1)
            + h * pec_sphere_implant_mode(outer_radius, lam, l, radius=radius,
                                          inner=inner)) / scale


def pec_sphere_implant_eigenvalues(l, n, *, radius, outer_radius, k, h,
                                   inner="neumann"):
    r"""First ``n`` roots ``λ`` of ``k R′(b) + h R(b) = 0`` for the hollow-sphere
    mode :func:`pec_sphere_implant_mode` (`IMP-7` step 4).

    Sign changes of the scaled determinant are bracketed on a uniform grid of
    spacing ``π/(b − a)/16`` (finer than the brief's ``π/(b − a)/8``), then
    refined with ``brentq``.  The determinant is smooth for ``λ > 0`` (no
    poles), so every sign change is a root; a bracket with two roots would be
    silent, which the step 4 test excludes by asserting the gaps.  With
    ``h > 0`` there is no zero eigenvalue.
    """
    from scipy.optimize import brentq

    a, b = float(radius), float(outer_radius)
    span = b - a
    dl = np.pi / span / 16.0
    roots = []
    lo = 1e-3 * dl
    f_lo = _imp7_robin_det(lo, l, a, b, k, h, inner)
    while len(roots) < n:
        hi = lo + dl
        f_hi = _imp7_robin_det(hi, l, a, b, k, h, inner)
        if f_lo == 0.0:
            roots.append(lo)
        elif np.sign(f_lo) != np.sign(f_hi) and f_hi != 0.0:
            roots.append(brentq(_imp7_robin_det, lo, hi,
                                args=(l, a, b, k, h, inner),
                                xtol=1e-15 * hi, rtol=4.0 * np.finfo(float).eps,
                                maxiter=500))
        lo, f_lo = hi, f_hi
    return np.asarray(roots[:n])


def _imp7_gl_nodes(a, b, n_panels, order=16):
    x, w = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(a, b, n_panels + 1)
    mid = 0.5 * (edges[1:] + edges[:-1])[:, None]
    half = 0.5 * (edges[1:] - edges[:-1])[:, None]
    return (mid + half * x).ravel(), (half * w).ravel()


def pec_sphere_implant_source_radial(r, l, *, radius, q_far, beta=1.0):
    """``q_l(r)`` of :func:`pec_sphere_implant_source` = ``q₀ + q₂P₂``:
    ``q₀ = q_far(1 + 2β²u²)``, ``q₂ = q_far(4βu + 2β²u²)``, ``u = (a/r)³``."""
    u = (radius / np.asarray(r, dtype=float)) ** 3
    if l == 0:
        return q_far * (1.0 + 2.0 * beta**2 * u**2)
    if l == 2:
        return q_far * (4.0 * beta * u + 2.0 * beta**2 * u**2)
    raise ValueError(f"l must be 0 or 2, got {l!r}")


def pec_sphere_implant_transient_modes(*, radius, outer_radius, k_gel, q_far, h,
                                       beta=1.0, n_terms=200, inner="neumann",
                                       t_init_offset=0.0):
    r"""Eigenvalues ``λ_{l,n}``, coefficients ``c_{l,n}`` and r²-norms of the
    hollow-sphere expansion (`IMP-7` step 4), for ``l`` ∈ {0, 2}.

    ``inner="neumann"`` (the adiabatic hole):
    ``c = ⟨θ_ss,l − θ_init,l, R⟩_{r²}/⟨R, R⟩_{r²}`` with ``θ_ss,l`` the exact
    steady radial part (:func:`pec_sphere_implant_radial`, adiabatic branch)
    and ``θ_init − T∞ = t_init_offset`` (l = 0 only; 0 when ``θ_init = T∞``).
    ``inner="dirichlet"`` (the control; no closed steady form here):
    ``c = ⟨q_l, R⟩_{r²}/(k λ²⟨R, R⟩_{r²})``, the projection of the Dirichlet
    steady state.  Quadrature: composite Gauss–Legendre, ``max(64, 2n)``
    panels × 16 nodes on (a, b).
    """
    a, b = float(radius), float(outer_radius)
    rq, wq = _imp7_gl_nodes(a, b, max(64, 2 * int(n_terms)))
    wr = wq * rq * rq
    coeffs = pec_sphere_implant_coefficients(radius=a, outer_radius=b, k_core=0.0,
                                             k_gel=k_gel, q_far=q_far, h=h,
                                             beta=beta)
    out = {"inner": inner, "coeffs": coeffs, "radius": a, "outer_radius": b,
           "k_gel": float(k_gel)}
    for l in (0, 2):
        lam = pec_sphere_implant_eigenvalues(l, n_terms, radius=a, outer_radius=b,
                                             k=k_gel, h=h, inner=inner)
        R = np.array([pec_sphere_implant_mode(rq, lm, l, radius=a, inner=inner)
                      for lm in lam])
        scale = np.max(np.abs(R), axis=1)
        R = R / scale[:, None]
        norms = R**2 @ wr
        if inner == "neumann":
            f = pec_sphere_implant_radial(rq, coeffs, l)
            if l == 0:
                f = f - t_init_offset
            c = (R @ (f * wr)) / norms
        else:
            ql = pec_sphere_implant_source_radial(rq, l, radius=a, q_far=q_far,
                                                  beta=beta)
            c = (R @ (ql * wr)) / (k_gel * lam**2 * norms)
        out[l] = {"lam": lam, "scale": scale, "norms": norms, "c": c}
    return out


def _imp7_points_rmu(points):
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    r = np.sqrt(np.sum(pts * pts, axis=1))
    mu = pts[:, 2] / r
    return r, 0.5 * (3.0 * mu * mu - 1.0)


def _imp7_mode_matrix(r, modes, l):
    m = modes[l]
    return np.array([pec_sphere_implant_mode(r, lm, l, radius=modes["radius"],
                                             inner=modes["inner"]) / s
                     for lm, s in zip(m["lam"], m["scale"])])


def pec_sphere_implant_transient_from_modes(points, t, modes, *, rho_c, t_inf,
                                            ls=(0, 2)):
    r"""``T(r, ϑ, t)`` from a :func:`pec_sphere_implant_transient_modes` dict.
    Neumann: ``θ = θ_ss − Σ_l P_l Σ_n c R e^{−κλ²t}`` (θ_ss exact).
    Dirichlet: ``θ = Σ_l P_l Σ_n c R (1 − e^{−κλ²t})`` (steady part as its
    series).  ``ls`` selects the harmonics (``(0,)`` is the isotropic-source
    control)."""
    r, p2 = _imp7_points_rmu(points)
    kappa = modes["k_gel"] / rho_c
    theta = np.zeros_like(r)
    for l in ls:
        m = modes[l]
        R = _imp7_mode_matrix(r, modes, l)
        dec = np.exp(-kappa * m["lam"] ** 2 * float(t))
        if modes["inner"] == "neumann":
            part = (pec_sphere_implant_radial(r, modes["coeffs"], l)
                    - (m["c"] * dec) @ R)
        else:
            part = (m["c"] * (1.0 - dec)) @ R
        theta = theta + (part if l == 0 else part * p2)
    return t_inf + theta


def pec_sphere_implant_transient_temperature(points, t, *, radius, outer_radius,
                                             k_gel, rho_c, q_far, h, t_inf,
                                             beta=1.0, n_terms=200):
    r"""Exact-series transient temperature around the adiabatic hole
    (`IMP-7` step 4, §10 I45): ``ρc ∂ₜθ = k Δθ + q`` in ``a < r < b`` with
    ``q`` :func:`pec_sphere_implant_source`, ``∂ᵣθ = 0`` at ``r = a``,
    ``−k∂ᵣθ = hθ`` at ``r = b``, ``θ(t = 0) = 0`` (``T = T∞``).

    ``θ = θ_ss − Σ_{l∈{0,2}} P_l(cos ϑ) Σ_n c_{l,n} R_{l,n}(r) exp(−κλ²_{l,n}t)``,
    ``κ = k/ρc``; ``θ_ss`` from :func:`pec_sphere_implant_steady_temperature`
    (``k_core = 0``), ``R_{l,n}`` :func:`pec_sphere_implant_mode`, ``λ``
    :func:`pec_sphere_implant_eigenvalues`, ``c``
    :func:`pec_sphere_implant_transient_modes`.  Truncated at ``n_terms`` per
    harmonic; (IC) / (Slope) converge with ``n_terms`` only."""
    modes = pec_sphere_implant_transient_modes(
        radius=radius, outer_radius=outer_radius, k_gel=k_gel, q_far=q_far, h=h,
        beta=beta, n_terms=n_terms)
    return pec_sphere_implant_transient_from_modes(points, t, modes, rho_c=rho_c,
                                                   t_inf=t_inf)


def pec_sphere_implant_duty_cycle_from_modes(points, t, modes, *, rho_c, t_inf,
                                             period, duty=0.5):
    r"""Direct modal duty-cycled transient (source on for ``duty·period`` at the
    start of each period), Neumann modes only.  Each modal amplitude ``A_n``
    obeys ``Ȧ = μ(c s(t) − A)``, ``μ = κλ²``, ``A(0) = 0``, advanced exactly
    phase by phase (the `MAT-14` ``_duty_cycle_modal`` pattern); then
    ``θ = s(t)θ_ss − Σ (c s(t) − A)R``."""
    if modes["inner"] != "neumann":
        raise ValueError("duty cycle is defined on the Neumann (adiabatic) modes")
    r, p2 = _imp7_points_rmu(points)
    kappa = modes["k_gel"] / rho_c
    t_on = duty * period
    edges = []
    s = 0.0
    while s < t:
        edges.append((s, min(s + t_on, t), 1.0))
        if s + t_on < t:
            edges.append((s + t_on, min(s + period, t), 0.0))
        s += period
    s_now = edges[-1][2]
    theta = np.zeros_like(r)
    for l in (0, 2):
        m = modes[l]
        mu_ = kappa * m["lam"] ** 2
        A = np.zeros_like(m["c"])
        for t0, t1, on in edges:
            e = np.exp(-mu_ * (t1 - t0))
            A = on * m["c"] + (A - on * m["c"]) * e
        R = _imp7_mode_matrix(r, modes, l)
        part = (s_now * pec_sphere_implant_radial(r, modes["coeffs"], l)
                - (s_now * m["c"] - A) @ R)
        theta = theta + (part if l == 0 else part * p2)
    return t_inf + theta
