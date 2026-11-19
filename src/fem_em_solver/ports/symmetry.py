"""Port bookkeeping across symmetry planes (`TH-16` step 2).

A lumped sheet cut by a symmetry plane keeps its constitutive law — the
sheet resistivity ``R = Z_p·w/h`` (:func:`~fem_em_solver.ports.lumped.
sheet_resistivity_ohm_per_square`) and the impressed field ``E_src = V/h`` —
while the *terminal* quantities of the retained piece rescale:

* ``"pmc"`` — the plane **contains** the drive direction ĥ (a magnetic wall,
  E_t even). It halves the sheet's width: ``w/2`` and ``Z_p×2``, so R is
  unchanged; V is unchanged and the retained piece carries half the current,
  so the half-model impedance is twice the full one.
* ``"pec"`` — the plane is **normal** to ĥ (an electric wall, E_t odd). It
  halves the terminal-to-terminal height: ``h/2``, ``Z_p/2`` and ``V/2``, so R
  and ``E_src`` are unchanged; the current is unchanged and the half-model
  impedance is half the full one.

:func:`full_impedance` maps a reduced-model impedance back to the full one.
Two cuts compose: ``pmc`` then ``pec`` gives the quarter, whose impedance
equals the full model's.
"""

from __future__ import annotations

import dataclasses

import numpy as np

from .lumped import LumpedSheetPortSpec

_CUTS = ("pmc", "pec")
_ORIENTATION_TOL = 1e-12


def _check_cut(cut: str) -> str:
    if cut not in _CUTS:
        raise ValueError(f"cut must be one of {_CUTS}, got {cut!r}")
    return cut


def halve_sheet_spec(spec: LumpedSheetPortSpec, cut: str, *, plane_normal) -> LumpedSheetPortSpec:
    """The spec of the half of ``spec`` retained on one side of a symmetry plane.

    ``plane_normal`` is the cut plane's normal; it is used only to check that
    the cut kind agrees with the drive direction ĥ — ``"pmc"`` needs ĥ in the
    plane (ĥ·n̂ = 0), ``"pec"`` needs ĥ along n̂ (|ĥ·n̂| = 1). Anything else
    raises: a sheet cut obliquely has no half-model bookkeeping.
    """
    _check_cut(cut)
    if spec.has_path:
        raise ValueError(
            f"port '{spec.port_id}': halving a spec with a path integral is not "
            "supported (the path would need cutting too)"
        )
    n = np.asarray(plane_normal, dtype=float)
    h = np.asarray(spec.drive_direction, dtype=float)
    if n.shape != (3,) or float(np.linalg.norm(n)) <= 0.0:
        raise ValueError(f"plane_normal must be a non-zero 3-vector, got {plane_normal!r}")
    cos = abs(float(np.dot(n / np.linalg.norm(n), h / np.linalg.norm(h))))
    if cut == "pmc":
        if cos > _ORIENTATION_TOL:
            raise ValueError(
                f"port '{spec.port_id}': a PMC cut must contain the drive direction "
                f"(|ĥ·n̂| = {cos:.3e}, expected 0)"
            )
        return dataclasses.replace(
            spec,
            sheet_width_m=float(spec.sheet_width_m) / 2.0,
            port_impedance_ohm=complex(spec.port_impedance_ohm) * 2.0,
        )
    if abs(cos - 1.0) > _ORIENTATION_TOL:
        raise ValueError(
            f"port '{spec.port_id}': a PEC cut must be normal to the drive direction "
            f"(|ĥ·n̂| = {cos:.3e}, expected 1)"
        )
    return dataclasses.replace(
        spec,
        gap_height_m=float(spec.gap_height_m) / 2.0,
        port_impedance_ohm=complex(spec.port_impedance_ohm) / 2.0,
        drive_voltage_v=complex(spec.drive_voltage_v) / 2.0,
    )


def full_impedance(z: complex, cut: str) -> complex:
    """Full-model impedance from a model reduced by one ``cut``: z/2 (PMC), 2z (PEC)."""
    _check_cut(cut)
    return complex(z) / 2.0 if cut == "pmc" else 2.0 * complex(z)
