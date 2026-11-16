"""`WF-1` step 1 (§9 item 256): a thin YAML layer over ``Study`` (schema v0).

``load_study(path)`` parses a YAML study file into the same ``Study`` object a
script builds in Python; ``run_study_file(path)`` solves it and writes a JSON
summary (plus a Touchstone file for a port study).  No solver code lives here
and nothing is imported from ``tests/``.

Schema v0 (``femem: 1``).  Units are SI and given by the key name; a unit
string is a type error.

* ``femem``: int, must be 1 (required).
* ``study``: ``frequency_hz`` (number, required), ``degree`` (int in {1, 2},
  optional — omitted means the ``Study``'s own default, which is 2),
  ``density_kg_m3`` (number, optional — omitted means the ``Study``'s default),
  ``sweep_hz`` (a strictly increasing list of numbers, optional, `WF-1` step 21:
  the solid ``four_port_birdcage`` with ``lumped_port_drive`` only).
* ``geometry``: ``kind`` (required):

  - ``sphere_in_box``: ``sphere_radius``, ``box_half_width``,
    ``resolution_sphere``, ``resolution_far`` (numbers, required),
    ``sphere_tag`` (int, optional);
  - ``four_port_birdcage``: ``preset: f_small`` (optional) plus any
    ``FourPortBirdcageGeometry`` field as an override.  There is no κ key
    (exact rectangles, κ off: operator ruling 2026-10-02) — an unknown key.
    ``as_hole`` (a YAML ``true`` / ``false``, `WF-1` step 4) is the only
    ``bool`` field settable from YAML; every other ``bool`` field is refused.
  - ``f_human_birdcage`` (`WF-1` step 7a, §9 item 390a): any
    ``FHumanBirdcageGeometry`` field as an override, no ``preset`` key;
    ``as_hole`` is required and must be ``true`` (the copper hole route only).

* ``implant`` (optional, `WF-15` step 10, §9 item 709): on a
  ``four_port_birdcage`` only, ``rod: {radius_m, length_m, centre_m, axis}``,
  all required, to ``FourPortBirdcageGeometry(implant_rod=...)``: a rod hole
  with a PEC wall.  Refused on every other geometry and with ``load_result``;
  beside ``capacitor_sheets`` on the hole route or with a ``conductor:``
  section (field ``implant``; admitted on the PEC four-port, `WF-15` step 13,
  §9 item 802) and at any resolved degree but 1, an omitted ``study.degree``
  included (field ``study.degree``; `WF-15` step 12, §9 item 729).

* ``material``: either ``preset`` (``low`` / ``mid`` / ``high``, through the
  ``Study``'s own string route) or ``sigma`` + ``epsilon_r`` (+ ``mu_r``).  The
  material frequency is the study's; there is no frequency key.
* ``excitation``: ``kind`` (required): ``quasi_static_sphere_field`` (``e0``)
  or ``lumped_port_drive`` (``z0_ohm``, ``port_impedance_ohm``,
  ``drive_voltage_v`` — a number or a ``[re, im]`` pair), or, on an
  ``f_human_birdcage`` only (`WF-1` step 7a, §9 item 390a),
  ``two_port_quadrature`` (``c_in_f`` required; ``esr_ohm``, ``feeds`` — two
  distinct ring-port ids — and ``z0_ohm`` optional, omitted means the
  ``TwoPortQuadratureDrive`` default); it needs a ``conductor:`` section.
  An ``f_human_birdcage`` also admits ``lumped_port_drive`` (`WF-1` step 12,
  §9 item 484: the 32-port ring route) without the quadrature keys
  (``c_in_f`` / ``esr_ohm`` / ``feeds`` are refused at ``excitation.kind``);
  its ``thermal:``, ``outputs.sar`` and ``outputs.b1_plus`` are refused.
  On a ``four_port_birdcage`` only (`WF-1` step 13, §9 item 496),
  ``capacitor_sheets`` (``c_f`` and ``sheets`` required — a list of distinct
  port ids ``P1``…``P{leg_count}`` that leaves at least one port to drive;
  ``esr_ohm`` optional, omitted means the ``CapacitorSheets`` default) builds
  ``CapacitorSheets(c_f, esr_ohm, on=sheets)``; the kept ports take
  ``LumpedPortDrive()``'s defaults.  The key is ``sheets``, not ``on``: YAML
  1.1 reads a bare ``on`` key as ``true``, and that key is refused at
  ``excitation.on``.  Its ``tune:``, ``outputs.sar`` and ``outputs.b1_plus``
  are refused.
* ``thermal_box`` (`WF-1` step 15, §9 item 521b): ``geometry: {kind:
  thermal_box, half_widths_m: [a, b, c], divisions: [i, j, k]}`` with the
  ``thermal_ball`` sections below, except that ``thermal.transient`` is
  optional (a steady box) and ``geometry.wall`` is refused.
  ``excitation.source: {kind: uniform | transverse_quadratic}`` (optional,
  `WF-1` step 16, §9 item 528) selects ``TransverseQuadraticHeatSource(q)``
  (mean q, shape 3(y² + z²)/(B² + C²)) in place of ``UniformHeatSource(q)``.
  ``thermal.dose: {points_m: [[x, y, z], …], window_s: w}`` (`WF-1` step 30,
  §9 item 878; needs ``thermal.transient``, a uniform source, no duty, no
  probes; points inside the gel box) adds a ``DoseReadout``: the summary's
  ``dose`` carries the point histories, `MAT-20`'s ``cem43`` and the
  trailing-window peak per point.
* ``thermal_ball`` (`WF-1` step 14, §9 item 500): a thermal-only Study, no
  field solve.  Sections ``geometry`` (``kind: thermal_ball``, ``radius``,
  ``resolution``, optional ``wall: {thickness_m, k_w_per_m_k,
  rho_c_j_per_m3_k}``), ``excitation`` (``kind: uniform_heat_source``,
  ``q_w_per_m3``) and ``thermal`` (the ``sphere_in_box`` keys plus a required
  ``transient: {dt_s, n_steps, rho_c_j_per_m3_k, theta, startup_half_steps}``
  and an optional ``duty: {period_s, duty}``) are required; ``study``,
  ``material``, ``conductor`` and ``tune`` are refused, as is
  ``thermal.source_scale``; ``outputs`` admits ``directory`` and ``probes``
  (`WF-1` steps 17 / 19).  Builds
  ``Study(ThermalBallGeometry(...), None, UniformHeatSource(q))``; the solve
  legs are gated in the real build.
* ``conductor`` (optional, `WF-1` step 4, §9 item 343): ``kind: leontovich``
  (required) and ``sigma_s_per_m`` (number > 0, optional — omitted means
  ``LeontovichWall``'s own default).  Builds ``Study(..., conductor=
  LeontovichWall(...))``; refused on a ``sphere_in_box`` geometry and on a
  ``four_port_birdcage`` without ``geometry.as_hole: true``; admitted on an
  ``f_human_birdcage`` (`WF-1` step 7a, §9 item 390a).
* ``thermal`` (optional, `WF-1` step 3, §9 item 295): ``k_w_per_m_k``,
  ``h_w_per_m2_k``, ``t_inf_c`` (numbers, required), ``source_scale``
  (number), ``degree`` (int) — omitted means ``RobinCooling``'s own default.
  Drives ``Study.temperature``; refused on a ``four_port_birdcage`` and an
  ``f_human_birdcage`` geometry.  On a ``sphere_in_box`` (`WF-1` step 29,
  §9 item 873): ``transient: {t_end_s, n_steps, rho_c_j_per_m3_k}``
  (required; ``theta``, ``startup_half_steps`` optional; dt = t_end_s /
  n_steps) drives ``Study.temperature(cooling, transient=...)``;
  ``thermal.duty``, ``thermal.initial`` and ``outputs.probes`` beside it are
  refused with the API's reasons.
* ``outputs`` (optional): ``directory``; on a ``four_port_birdcage`` only
  (`WF-1` step 5, §9 item 344): ``sar: {mass_g, centres}`` (``centres`` a
  list of ``[x, y, z]`` or ``{ring_radius_m: r}`` — one centre per port at
  ``(r cos az, r sin az, 0)``, ``az`` off the solved sheets) drives
  ``Study.sar`` (every drive at every centre) and ``b1_plus: {points}`` drives
  ``Study.b1_plus_at``; both land in the JSON summary.  They reuse the cached
  field-keeping solve, a second solve beside ``Study.solve`` (disclosed).
  On an ``f_human_birdcage`` (`WF-1` step 7a, §9 item 390a) ``sar`` and
  ``b1_plus`` are admitted per sense and reuse the tuned solve; the
  ``{ring_radius_m}`` centre form is refused there (literal centres only).
  On a ``sphere_in_box`` with a ``thermal:`` section only (`WF-1` step 6, §9
  item 345): ``probes: [{centre: [x, y, z], radius: r > 0}, …]`` drives
  ``Study.probe(probes, cooling)`` (``levels`` not exposed — the API default)
  and lands in the summary as ``probes``.  The probe values come from the
  thermal solve inside ``Study.probe``, a second (cheap) thermal solve beside
  the one the ``temperature`` block reports; the EM solve is reused (disclosed).
  On a ``four_port_birdcage`` (`WF-1` step 20, §9 item 564): ``probes:
  {points: [[x, y, z], …]}``, each strictly inside the phantom, drives
  ``post.sar.point_sar`` (owner-rule default) on the cached kept fields and
  lands in the summary as ``probes`` (point SAR and ``|E|`` per drive).
* ``tune`` (optional, `WF-1` step 11, §9 item 442): ``k`` (int, 1 or 2;
  omitted means ``Study.tune``'s own default).  Top-level in the file, carried
  in the parsed outputs under ``"tune"``; defined on a ``four_port_birdcage``
  and (`WF-1` step 12, §9 item 484) on an ``f_human_birdcage`` with
  ``excitation.kind: lumped_port_drive`` (the 32-port ring read).  Drives ``Study.tune`` on the solved Study (no second solve) and lands
  in the summary as ``tune``.  ``at_hz`` (number, `WF-1` step 27, §9 item 730)
  selects the nearest frequency of a multi-frequency ``load_result`` file, and
  the block's ``frequency_hz`` is that selected frequency; anywhere else it is
  refused by field ``tune.at_hz``, and a multi-frequency file with no ``at_hz``
  is refused by the same field (`OPS-76`, §9 item 754).
* ``load_result`` (optional, `WF-1` step 24, §9 item 674): a path to a
  single-frequency Touchstone file (relative paths resolve against the YAML
  file's directory).  Admitted only on a ``four_port_birdcage``, (`WF-1`
  step 25, §9 item 697) an ``f_clinical_birdcage`` or (`WF-1` step 26, §9 item
  714) an ``f_human_birdcage`` with ``excitation.kind: lumped_port_drive`` and
  a ``tune`` section; refused at parse with ``conductor`` (except on an
  ``f_human_birdcage``, where it documents the saved column and nothing is
  solved), ``thermal``, ``study.sweep_hz``, ``outputs.sar``,
  ``outputs.b1_plus`` or ``outputs.probes`` (and by ``--backend palace``).
  ``run_study_file`` then calls ``Study.load_result`` + ``Study.tune`` with no
  mesh and no solve, and writes the same ``tune`` summary block.

Validation order is fixed and the first error raises: the version, unknown
top-level keys, missing sections; then study → geometry → material →
excitation → conductor → thermal → outputs → tune, and within a section unknown keys, missing keys, then
types and values.  Every rank parses the same file, so a ``StudyConfigError``
is raised identically on every rank with no collective before it.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional

from mpi4py import MPI

from .materials import GelledSalinePhantomMaterial
from .study import (
    CapacitorSheets,
    ContainerWall,
    DoseReadout,
    DutyCycle,
    FClinicalBirdcageGeometry,
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    QuasiStaticSphereField,
    RobinCooling,
    SphereInBoxGeometry,
    Study,
    ThermalBallGeometry,
    ThermalBoxGeometry,
    TransientSpec,
    TransverseQuadraticHeatSource,
    TwoPortQuadratureDrive,
    UniformHeatSource,
    _f_clinical_tier,
)

__all__ = ["StudyConfigError", "load_study", "run_study_file"]

SCHEMA_VERSION = 1
_TOP_KEYS = (
    "femem", "study", "geometry", "material", "excitation", "conductor", "thermal", "outputs", "tune",
    "load_result",  # `WF-1` step 24 (§9 item 674, U99)
    "implant",  # `WF-15` step 10 (§9 item 709, U97′)
)
_CONDUCTOR_KINDS = ("leontovich",)
# The only bool geometry field settable from YAML (`WF-1` step 4, §9 item 343).
_SETTABLE_BOOLS = ("as_hole",)
_REQUIRED_SECTIONS = ("study", "geometry", "material", "excitation")
_YAML_FLOAT_HINT = (
    "YAML reads `64.0e6` as a string; write `6.4e+07` or `64000000.0`"
)
_SPHERE_NUMBERS = ("sphere_radius", "box_half_width", "resolution_sphere", "resolution_far")
_EXCITATION_FOR = {
    "sphere_in_box": "quasi_static_sphere_field",
    "four_port_birdcage": "lumped_port_drive",
    # `WF-1` step 7a (§9 item 390a): the tuned copper F-human, hole route only.
    "f_human_birdcage": "two_port_quadrature",
    # `WF-13` step 12 (§9 item 663): F-clinical-1.5T, estimate only.
    "f_clinical_birdcage": "lumped_port_drive",
}


class StudyConfigError(ValueError):
    """A study file that does not parse to a ``Study``; ``.field`` names where."""

    def __init__(self, source, field: str, reason: str):
        self.source = str(source)
        self.field = field
        self.reason = reason
        super().__init__(f"{self.source}: {field}: {reason}")


# ---------------------------------------------------------------------------
# value checks
# ---------------------------------------------------------------------------


def _number(value, field_name, source) -> float:
    # bool is an int subclass: refuse it before the int test.
    if isinstance(value, bool):
        raise StudyConfigError(source, field_name, f"expected a number, got bool {value!r}")
    if isinstance(value, str):
        raise StudyConfigError(
            source,
            field_name,
            f"expected a number, got the string {value!r} ({_YAML_FLOAT_HINT}; "
            "units are given by the key name, never as a string)",
        )
    if not isinstance(value, (int, float)):
        raise StudyConfigError(
            source, field_name, f"expected a number, got {type(value).__name__}"
        )
    return float(value)


def _integer(value, field_name, source) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise StudyConfigError(
            source, field_name, f"expected an integer, got {type(value).__name__} {value!r}"
        )
    return int(value)


def _string(value, field_name, source) -> str:
    if not isinstance(value, str):
        raise StudyConfigError(
            source, field_name, f"expected a string, got {type(value).__name__} {value!r}"
        )
    return value


def _section(mapping, name, source) -> dict:
    value = mapping[name]
    if not isinstance(value, dict):
        raise StudyConfigError(
            source, name, f"expected a mapping, got {type(value).__name__}"
        )
    return value


def _unknown(section: dict, allowed, prefix, source):
    for key in section:
        if key not in allowed:
            raise StudyConfigError(
                source,
                f"{prefix}.{key}",
                f"unknown key (allowed: {', '.join(sorted(allowed))})",
            )


def _missing(section: dict, required, prefix, source):
    for key in required:
        if key not in section:
            raise StudyConfigError(source, f"{prefix}.{key}", "required key is missing")


def _kind(section: dict, kinds, prefix, source) -> str:
    if "kind" not in section:
        raise StudyConfigError(source, f"{prefix}.kind", "required key is missing")
    kind = _string(section["kind"], f"{prefix}.kind", source)
    if kind not in kinds:
        raise StudyConfigError(
            source, f"{prefix}.kind", f"unknown kind {kind!r} (allowed: {', '.join(kinds)})"
        )
    return kind


# ---------------------------------------------------------------------------
# sections
# ---------------------------------------------------------------------------


def _parse_study(sec, source) -> dict:
    p = "study"
    _unknown(sec, ("frequency_hz", "degree", "density_kg_m3", "sweep_hz"), p, source)
    _missing(sec, ("frequency_hz",), p, source)
    kwargs = {"frequency_hz": _number(sec["frequency_hz"], f"{p}.frequency_hz", source)}
    if kwargs["frequency_hz"] <= 0.0:
        raise StudyConfigError(source, f"{p}.frequency_hz", "must be > 0")
    # An omitted degree / density is NOT filled in here: the Study's own
    # default applies (degree 2).  A loader default would drift silently if
    # the Study's ever changed; (Deg) in the WF-1 test checks the omitted case.
    if "degree" in sec:
        degree = _integer(sec["degree"], f"{p}.degree", source)
        if degree not in (1, 2):
            raise StudyConfigError(source, f"{p}.degree", f"must be 1 or 2, got {degree}")
        kwargs["degree"] = degree
    if "density_kg_m3" in sec:
        kwargs["density_kg_m3"] = _number(sec["density_kg_m3"], f"{p}.density_kg_m3", source)
    return kwargs


def _parse_geometry(sec, source):
    p = "geometry"
    kind = _kind(sec, tuple(_EXCITATION_FOR), p, source)
    if kind == "sphere_in_box":
        _unknown(sec, ("kind",) + _SPHERE_NUMBERS + ("sphere_tag",), p, source)
        _missing(sec, _SPHERE_NUMBERS, p, source)
        kwargs = {k: _number(sec[k], f"{p}.{k}", source) for k in _SPHERE_NUMBERS}
        if "sphere_tag" in sec:
            kwargs["sphere_tag"] = _integer(sec["sphere_tag"], f"{p}.sphere_tag", source)
        return kind, SphereInBoxGeometry(**kwargs)

    if kind == "f_clinical_birdcage":
        # `WF-13` step 12 (§9 item 663): the 1.5 T r1 fixture, no overrides but
        # (`WF-13` step 15, §9 item 782) ``field_strength_t``, 1.5 or 3.0 only.
        _unknown(sec, ("kind", "field_strength_t"), p, source)
        if "field_strength_t" not in sec:
            return kind, FClinicalBirdcageGeometry()
        b0 = _number(sec["field_strength_t"], f"{p}.field_strength_t", source)
        if b0 not in (1.5, 3.0):
            raise StudyConfigError(source, f"{p}.field_strength_t",
                                   f"must be 1.5 or 3.0, got {b0}")
        return kind, FClinicalBirdcageGeometry(field_strength_t=float(b0))

    if kind == "f_human_birdcage":
        # `WF-1` step 7a (§9 item 390a): copied from the four-port loop below,
        # not factored out of it (a shared helper would edit that path).
        defaults = FHumanBirdcageGeometry()
        names = [f.name for f in dataclasses.fields(FHumanBirdcageGeometry)]
        _unknown(sec, ["kind"] + names, p, source)
        _missing(sec, ("as_hole",), p, source)
        overrides = {}
        for name in names:
            if name not in sec:
                continue
            default = getattr(defaults, name)
            where = f"{p}.{name}"
            value = sec[name]
            if isinstance(default, bool):
                if name not in _SETTABLE_BOOLS:
                    raise StudyConfigError(source, where, "not settable from YAML")
                if not isinstance(value, bool):
                    raise StudyConfigError(
                        source,
                        where,
                        f"expected true or false, got {type(value).__name__} {value!r}",
                    )
                overrides[name] = value
                continue
            if isinstance(default, int):
                overrides[name] = _integer(value, where, source)
            elif isinstance(default, float):
                overrides[name] = _number(value, where, source)
            else:  # pragma: no cover - every field is int / float / bool today
                raise StudyConfigError(source, where, "not settable from YAML")
        if overrides["as_hole"] is not True:
            raise StudyConfigError(
                source,
                f"{p}.as_hole",
                "must be true (the YAML f_human_birdcage route is the copper hole route only)",
            )
        return kind, dataclasses.replace(defaults, **overrides)

    defaults = FourPortBirdcageGeometry.f_small()
    names = [f.name for f in dataclasses.fields(FourPortBirdcageGeometry)]
    _unknown(sec, ["kind", "preset"] + names, p, source)
    if "preset" in sec:
        preset = _string(sec["preset"], f"{p}.preset", source)
        if preset != "f_small":
            raise StudyConfigError(source, f"{p}.preset", f"unknown preset {preset!r} (allowed: f_small)")
    overrides = {}
    for name in names:
        if name not in sec:
            continue
        default = getattr(defaults, name)
        where = f"{p}.{name}"
        value = sec[name]
        if isinstance(default, bool):
            if name not in _SETTABLE_BOOLS:
                raise StudyConfigError(source, where, "not settable from YAML")
            # YAML `true` / `false` only: an int (`1`) is refused, not coerced.
            if not isinstance(value, bool):
                raise StudyConfigError(
                    source,
                    where,
                    f"expected true or false, got {type(value).__name__} {value!r}",
                )
            overrides[name] = value
            continue
        if isinstance(default, int):
            overrides[name] = _integer(value, where, source)
        elif isinstance(default, float):
            overrides[name] = _number(value, where, source)
        elif isinstance(default, tuple):
            if not isinstance(value, (list, tuple)) or len(value) != len(default):
                raise StudyConfigError(
                    source, where, f"expected a list of {len(default)} numbers"
                )
            # A YAML list is a list; the dataclass default is a tuple.
            overrides[name] = tuple(
                _number(v, f"{where}[{i}]", source) for i, v in enumerate(value)
            )
        else:  # pragma: no cover - every field is int / float / tuple today
            raise StudyConfigError(source, where, "not settable from YAML")
    return kind, dataclasses.replace(defaults, **overrides)


def _clinical_rod(sec):
    """The ``implant:`` section as entry 48's rod mapping, or ``None`` for
    anything else (`WF-13` step 16); never raises, so the refusal stays as is."""
    from .study import F_CLINICAL_ROD_ENTRY48

    if not isinstance(sec, dict) or set(sec) != {"rod"} or not isinstance(sec["rod"], dict):
        return None
    rod = sec["rod"]
    if set(rod) != {"radius_m", "length_m", "centre_m", "axis"}:
        return None
    try:
        got = {"radius": float(rod["radius_m"]), "length": float(rod["length_m"]),
               "centre": tuple(float(v) for v in rod["centre_m"]),
               "axis": tuple(float(v) for v in rod["axis"])}
    except (TypeError, ValueError):
        return None
    return got if got == F_CLINICAL_ROD_ENTRY48 else None


def _parse_implant(sec, geometry_kind, geometry, source):
    """``implant: {rod: {...}}`` to the generator's ``implant_rod`` mapping
    (`WF-15` step 10, §9 item 709): a ``four_port_birdcage`` only."""
    p = "implant"
    if geometry_kind == "f_clinical_birdcage":
        # `WF-13` step 16 (§9 item 783): entry 48's rod at 1.5 T only, estimate-only;
        # every other implant here keeps the message below unchanged.
        rod = _clinical_rod(sec)
        if rod is not None and geometry == FClinicalBirdcageGeometry():
            return dataclasses.replace(geometry, implant_rod=rod)
    if geometry_kind != "four_port_birdcage":
        raise StudyConfigError(
            source, p, f"an implant is defined only on a four_port_birdcage, not {geometry_kind!r}"
        )
    _unknown(sec, ("rod",), p, source)
    _missing(sec, ("rod",), p, source)
    p = "implant.rod"
    rod = _section(sec, "rod", source)
    keys = ("radius_m", "length_m", "centre_m", "axis")
    _unknown(rod, keys, p, source)
    _missing(rod, keys, p, source)
    implant_rod = {
        "radius": _number(rod["radius_m"], f"{p}.radius_m", source),
        "length": _number(rod["length_m"], f"{p}.length_m", source),
        "centre": _triple(rod["centre_m"], f"{p}.centre_m", source),
        "axis": _triple(rod["axis"], f"{p}.axis", source),
    }
    return dataclasses.replace(geometry, implant_rod=implant_rod)


def _parse_material(sec, frequency_hz, source):
    p = "material"
    _unknown(sec, ("preset", "sigma", "epsilon_r", "mu_r"), p, source)
    explicit = [k for k in ("sigma", "epsilon_r", "mu_r") if k in sec]
    if "preset" in sec:
        if explicit:
            raise StudyConfigError(
                source, p, f"give either `preset` or `sigma` + `epsilon_r`, not both ({explicit})"
            )
        preset = _string(sec["preset"], f"{p}.preset", source)
        if preset not in GelledSalinePhantomMaterial.preset_names():
            raise StudyConfigError(
                source,
                f"{p}.preset",
                f"unknown preset {preset!r} (allowed: "
                f"{', '.join(GelledSalinePhantomMaterial.preset_names())})",
            )
        return preset  # the Study's own string route
    _missing(sec, ("sigma", "epsilon_r"), p, source)
    sigma = _number(sec["sigma"], f"{p}.sigma", source)
    epsilon_r = _number(sec["epsilon_r"], f"{p}.epsilon_r", source)
    if "mu_r" in sec:
        mu_r = _number(sec["mu_r"], f"{p}.mu_r", source)
        return GelledSalinePhantomMaterial(sigma, epsilon_r, frequency_hz, mu_r)
    return GelledSalinePhantomMaterial(sigma, epsilon_r, frequency_hz)


def _parse_excitation(sec, geometry_kind, geometry, source):
    p = "excitation"
    kind = _kind(sec, tuple(_EXCITATION_FOR.values()) + ("capacitor_sheets",), p, source)
    if kind != _EXCITATION_FOR[geometry_kind] and not (
        kind == "lumped_port_drive" and geometry_kind == "f_human_birdcage"  # `WF-1` step 12
    ) and not (
        kind == "capacitor_sheets" and geometry_kind == "four_port_birdcage"  # `WF-1` step 13
    ):
        raise StudyConfigError(
            source,
            f"{p}.kind",
            f"{kind!r} does not drive a {geometry_kind!r} geometry "
            f"(use {_EXCITATION_FOR[geometry_kind]!r})",
        )
    if kind == "quasi_static_sphere_field":
        _unknown(sec, ("kind", "e0"), p, source)
        kwargs = {}
        if "e0" in sec:
            kwargs["e0"] = _number(sec["e0"], f"{p}.e0", source)
        return QuasiStaticSphereField(**kwargs)
    if kind == "two_port_quadrature":  # `WF-1` step 7a (§9 item 390a)
        return _parse_two_port_quadrature(sec, geometry, source)
    if kind == "capacitor_sheets":  # `WF-1` step 13 (§9 item 496)
        return _parse_capacitor_sheets(sec, geometry, source)
    if geometry_kind == "f_human_birdcage":  # `WF-1` step 12 (§9 item 484), E1b
        stray = [key for key in ("c_in_f", "esr_ohm", "feeds") if key in sec]
        if stray:
            raise StudyConfigError(
                source,
                f"{p}.kind",
                "lumped_port_drive takes no c_in_f / esr_ohm / feeds "
                f"(two_port_quadrature's keys; found {', '.join(stray)})",
            )
    _unknown(sec, ("kind", "z0_ohm", "port_impedance_ohm", "drive_voltage_v"), p, source)
    kwargs = {}
    for key in ("z0_ohm", "port_impedance_ohm"):
        if key in sec:
            kwargs[key] = _number(sec[key], f"{p}.{key}", source)
    if "drive_voltage_v" in sec:
        value = sec["drive_voltage_v"]
        where = f"{p}.drive_voltage_v"
        if isinstance(value, (list, tuple)):
            if len(value) != 2:
                raise StudyConfigError(source, where, "expected a number or a [re, im] pair")
            kwargs["drive_voltage_v"] = complex(
                _number(value[0], f"{where}[0]", source), _number(value[1], f"{where}[1]", source)
            )
        else:
            kwargs["drive_voltage_v"] = complex(_number(value, where, source))
    return LumpedPortDrive(**kwargs)


def _parse_two_port_quadrature(sec, geometry, source) -> TwoPortQuadratureDrive:
    """``two_port_quadrature`` → ``TwoPortQuadratureDrive`` (`WF-1` step 7a, §9 item 390a).

    Every value is checked here, so the dataclass's own ``ValueError`` never
    escapes bare; omitted ``esr_ohm`` / ``feeds`` / ``z0_ohm`` are NOT filled
    in (the dataclass defaults apply).  A feed must be a ring port of the
    parsed geometry (the tuned solve's own refusal, raised here at parse).
    """
    import math

    p = "excitation"
    _unknown(sec, ("kind", "c_in_f", "esr_ohm", "feeds", "z0_ohm"), p, source)
    _missing(sec, ("c_in_f",), p, source)
    c_in_f = _number(sec["c_in_f"], f"{p}.c_in_f", source)
    if not (math.isfinite(c_in_f) and c_in_f > 0.0):
        raise StudyConfigError(source, f"{p}.c_in_f", f"must be finite and > 0, got {c_in_f!r}")
    kwargs = {"c_in_f": c_in_f}
    if "esr_ohm" in sec:
        esr = _number(sec["esr_ohm"], f"{p}.esr_ohm", source)
        if not (math.isfinite(esr) and esr >= 0.0):
            raise StudyConfigError(source, f"{p}.esr_ohm", f"must be finite and >= 0, got {esr!r}")
        kwargs["esr_ohm"] = esr
    if "feeds" in sec:
        value = sec["feeds"]
        where = f"{p}.feeds"
        if not isinstance(value, (list, tuple)) or len(value) != 2:
            raise StudyConfigError(source, where, f"expected a list of two port ids, got {value!r}")
        feeds = tuple(_string(v, f"{where}[{i}]", source) for i, v in enumerate(value))
        if feeds[0] == feeds[1]:
            raise StudyConfigError(source, where, f"the two feeds must differ, got {list(feeds)!r}")
        valid = {f"P{i}" for i in geometry.ring_ordinals()}
        for i, fid in enumerate(feeds):
            if fid not in valid:
                raise StudyConfigError(
                    source, f"{where}[{i}]", f"{fid!r} is not a ring port of this geometry"
                )
        kwargs["feeds"] = feeds
    if "z0_ohm" in sec:
        z0 = _number(sec["z0_ohm"], f"{p}.z0_ohm", source)
        if not (math.isfinite(z0) and z0 > 0.0):
            raise StudyConfigError(source, f"{p}.z0_ohm", f"must be finite and > 0, got {z0!r}")
        kwargs["z0_ohm"] = z0
    return TwoPortQuadratureDrive(**kwargs)


def _parse_capacitor_sheets(sec, geometry, source) -> CapacitorSheets:
    """``capacitor_sheets`` → ``CapacitorSheets`` (`WF-1` step 13, §9 item 496).

    Every value is checked here, so the dataclass's own ``ValueError`` never
    escapes bare; an omitted ``esr_ohm`` is NOT filled in (the dataclass
    default applies).  ``sheets`` maps to ``CapacitorSheets.on`` (a bare
    ``on`` key loads as YAML 1.1 ``true`` and is refused at ``excitation.on``);
    the kept ports take ``LumpedPortDrive()``'s defaults.
    """
    import math

    p = "excitation"
    if True in sec:
        raise StudyConfigError(
            source, f"{p}.on", "YAML 1.1 reads a bare `on` key as true; write `sheets: [...]`"
        )
    _unknown(sec, ("kind", "c_f", "esr_ohm", "sheets"), p, source)
    _missing(sec, ("c_f", "sheets"), p, source)
    c_f = _number(sec["c_f"], f"{p}.c_f", source)
    if not (math.isfinite(c_f) and c_f > 0.0):
        raise StudyConfigError(source, f"{p}.c_f", f"must be finite and > 0, got {c_f!r}")
    kwargs = {"c_f": c_f}
    if "esr_ohm" in sec:
        esr = _number(sec["esr_ohm"], f"{p}.esr_ohm", source)
        if not (math.isfinite(esr) and esr >= 0.0):
            raise StudyConfigError(source, f"{p}.esr_ohm", f"must be finite and >= 0, got {esr!r}")
        kwargs["esr_ohm"] = esr
    value = sec["sheets"]
    where = f"{p}.sheets"
    if not isinstance(value, (list, tuple)):
        raise StudyConfigError(source, where, f"expected a list of port ids, got {value!r}")
    if not value:
        raise StudyConfigError(source, where, "must name at least one port sheet")
    valid = {f"P{i + 1}" for i in range(int(geometry.leg_count))}
    sheets = []
    for i, v in enumerate(value):
        sid = _string(v, f"{where}[{i}]", source)
        if sid not in valid:
            raise StudyConfigError(
                source, f"{where}[{i}]", f"{sid!r} is not a port of this geometry "
                f"(allowed: {', '.join(sorted(valid, key=lambda s: int(s[1:])))})"
            )
        sheets.append(sid)
    if len(set(sheets)) != len(sheets):
        raise StudyConfigError(source, where, f"names a port sheet twice: {sheets!r}")
    if len(sheets) >= len(valid):
        raise StudyConfigError(source, where, "leaves no port to drive")
    return CapacitorSheets(on=tuple(sheets), **kwargs)


def _parse_conductor(sec, geometry_kind, geometry, source) -> LeontovichWall:
    """``conductor:`` → ``LeontovichWall`` (`WF-1` step 4, §9 item 343).

    Refused at ``conductor`` on a sphere and without ``geometry.as_hole: true``
    (``_conductor_setup``'s own ``ValueError``, raised here at parse, before
    any mesh).  An omitted ``sigma_s_per_m`` is NOT filled in: the dataclass
    default applies, as for ``study.degree``.
    """
    p = "conductor"
    if geometry_kind not in ("four_port_birdcage", "f_human_birdcage"):
        raise StudyConfigError(
            source, p, f"a conductor is defined only on a four_port_birdcage, not {geometry_kind!r}"
        )
    if not getattr(geometry, "as_hole", False):
        raise StudyConfigError(
            source, p, "a leontovich conductor needs `geometry.as_hole: true` (the hole route)"
        )
    _unknown(sec, ("kind", "sigma_s_per_m"), p, source)
    _kind(sec, _CONDUCTOR_KINDS, p, source)
    if "sigma_s_per_m" in sec:
        sigma = _number(sec["sigma_s_per_m"], f"{p}.sigma_s_per_m", source)
        if not sigma > 0.0:
            raise StudyConfigError(source, f"{p}.sigma_s_per_m", f"must be > 0, got {sigma!r}")
        return LeontovichWall(sigma_s_per_m=sigma)
    return LeontovichWall()


# `WF-1` step 14 (§9 item 500): a thermal-only ball (no frequency, no material).
_BALL_REQUIRED = ("geometry", "excitation", "thermal")
_BALL_REFUSED = {
    "study": "a thermal_ball study has no frequency / EM degree / density",
    "material": "a thermal_ball study has no EM material: k is `thermal.k_w_per_m_k`, "
    "rho_c is `thermal.transient.rho_c_j_per_m3_k`",
    "conductor": "a thermal_ball study has no conductor",
    "tune": "a thermal_ball study has no tune route",
    "implant": "a thermal study has no implant",  # `WF-15` step 10
}
_WALL_KEYS = ("thickness_m", "k_w_per_m_k", "rho_c_j_per_m3_k")


_THERMAL_REQUIRED = ("k_w_per_m_k", "h_w_per_m2_k", "t_inf_c")


def _parse_thermal(sec, geometry_kind, source, implant=False) -> RobinCooling:
    """``thermal:`` → ``RobinCooling`` (`WF-1` step 3, §9 item 295).

    Every value is checked here, so ``RobinCooling``'s own ``ValueError`` never
    escapes bare; omitted ``source_scale`` / ``degree`` are NOT filled in (the
    dataclass defaults apply, as for ``study.degree``).  On a
    ``four_port_birdcage`` with ``implant:`` (`WF-15` step 14, §9 items 831 /
    884) it is a ``_RodThermal``.
    """
    p = "thermal"
    if geometry_kind == "four_port_birdcage" and implant:
        return _parse_thermal_rod(sec, source)
    if geometry_kind == "four_port_birdcage":
        raise StudyConfigError(
            source, p, "temperature is not defined on a port Study (four_port_birdcage)"
        )
    elif geometry_kind == "f_human_birdcage":  # `WF-1` step 10a (§9 item 431a)
        return _parse_thermal_f_human(sec, source)
    elif geometry_kind == "sphere_in_box" and any(  # `WF-1` step 29 (§9 item 873)
            k in sec for k in ("transient", "duty", "initial")):
        return _parse_thermal_sphere_transient(sec, source)
    _unknown(sec, _THERMAL_REQUIRED + ("source_scale", "degree"), p, source)
    _missing(sec, _THERMAL_REQUIRED, p, source)
    kwargs = {k: _number(sec[k], f"{p}.{k}", source) for k in _THERMAL_REQUIRED}
    for key in ("k_w_per_m_k", "h_w_per_m2_k"):
        if not kwargs[key] > 0.0:
            raise StudyConfigError(source, f"{p}.{key}", f"must be > 0, got {kwargs[key]!r}")
    if "source_scale" in sec:
        scale = _number(sec["source_scale"], f"{p}.source_scale", source)
        if not scale > 0.0:
            raise StudyConfigError(source, f"{p}.source_scale", f"must be > 0, got {scale!r}")
        kwargs["source_scale"] = scale
    if "degree" in sec:
        degree = _integer(sec["degree"], f"{p}.degree", source)
        if degree < 1:
            raise StudyConfigError(source, f"{p}.degree", f"must be >= 1, got {degree}")
        kwargs["degree"] = degree
    return RobinCooling(**kwargs)


@dataclass(frozen=True)
class _FHumanThermal:
    """A parsed F-human ``thermal:`` section (`WF-1` step 10a, §9 item 431a).

    ``window_s`` / ``sense`` are ``None`` when omitted: the API default
    applies, never a loader default.
    """

    cooling: RobinCooling
    transient: Any
    balls: dict
    window_s: Optional[float] = None
    sense: Optional[str] = None


def _parse_thermal_f_human(sec, source) -> _FHumanThermal:
    """``thermal:`` on an ``f_human_birdcage``: a transient Robin-cooled
    temperature of the phantom, 10 g ball means and a running maximum.

    Order as everywhere: unknown keys, missing, types, values.  A section
    without ``transient`` (or ``balls``) is refused at ``thermal``.
    """
    import math

    from .study import TransientSpec

    p = "thermal"
    allowed = _THERMAL_REQUIRED + ("source_scale", "degree", "transient", "balls",
                                   "window_s", "sense")
    _unknown(sec, allowed, p, source)
    if "transient" not in sec:
        raise StudyConfigError(
            source, p, "f_human_birdcage needs `thermal.transient` (a transient Study)"
        )
    if "balls" not in sec:
        raise StudyConfigError(
            source, p, "f_human_birdcage needs `thermal.balls` (10 g ball centres)"
        )
    base = {k: sec[k] for k in sec if k not in ("transient", "balls", "window_s", "sense")}
    cooling = _parse_thermal(base, "sphere_in_box", source)

    tp = f"{p}.transient"
    tsec = sec["transient"]
    if not isinstance(tsec, dict):
        raise StudyConfigError(source, tp, f"expected a mapping, got {type(tsec).__name__}")
    t_req = ("dt_s", "n_steps", "rho_c_j_per_m3_k")
    _unknown(tsec, t_req + ("theta", "startup_half_steps"), tp, source)
    _missing(tsec, t_req, tp, source)
    dt = _number(tsec["dt_s"], f"{tp}.dt_s", source)
    n_steps = _integer(tsec["n_steps"], f"{tp}.n_steps", source)
    rho_c = _number(tsec["rho_c_j_per_m3_k"], f"{tp}.rho_c_j_per_m3_k", source)
    tkw = {}
    if "theta" in tsec:
        tkw["theta"] = _number(tsec["theta"], f"{tp}.theta", source)
    if "startup_half_steps" in tsec:
        tkw["startup_half_steps"] = _integer(
            tsec["startup_half_steps"], f"{tp}.startup_half_steps", source
        )
    if not (math.isfinite(dt) and dt > 0.0):
        raise StudyConfigError(source, f"{tp}.dt_s", f"must be finite and > 0, got {dt!r}")
    if n_steps < 1:
        raise StudyConfigError(source, f"{tp}.n_steps", f"must be >= 1, got {n_steps}")
    if not (math.isfinite(rho_c) and rho_c > 0.0):
        raise StudyConfigError(
            source, f"{tp}.rho_c_j_per_m3_k", f"must be finite and > 0, got {rho_c!r}"
        )
    if "theta" in tkw and not 0.0 <= tkw["theta"] <= 1.0:
        raise StudyConfigError(source, f"{tp}.theta", f"must be in [0, 1], got {tkw['theta']!r}")
    if "startup_half_steps" in tkw and (
        tkw["startup_half_steps"] < 0 or tkw["startup_half_steps"] % 2
    ):
        raise StudyConfigError(
            source, f"{tp}.startup_half_steps",
            f"must be a non-negative even integer, got {tkw['startup_half_steps']}",
        )
    transient = TransientSpec(dt_s=dt, n_steps=n_steps, rho_c_j_per_m3_k=rho_c, **tkw)

    bp = f"{p}.balls"
    bsec = sec["balls"]
    if not isinstance(bsec, dict):
        raise StudyConfigError(source, bp, f"expected a mapping, got {type(bsec).__name__}")
    _unknown(bsec, ("mass_g", "centres"), bp, source)
    _missing(bsec, ("mass_g", "centres"), bp, source)
    mass_g = _number(bsec["mass_g"], f"{bp}.mass_g", source)
    if not (math.isfinite(mass_g) and mass_g > 0.0):
        raise StudyConfigError(source, f"{bp}.mass_g", f"must be finite and > 0, got {mass_g!r}")
    balls = {"mass_g": mass_g, "centres": _triples(bsec["centres"], f"{bp}.centres", source)}

    window_s = None
    if "window_s" in sec:
        window_s = _number(sec["window_s"], f"{p}.window_s", source)
        if not (math.isfinite(window_s) and window_s > 0.0):
            raise StudyConfigError(
                source, f"{p}.window_s", f"must be finite and > 0, got {window_s!r}"
            )
    sense = None
    if "sense" in sec:
        sense = _string(sec["sense"], f"{p}.sense", source)
        if sense not in ("ccw", "cw"):
            raise StudyConfigError(source, f"{p}.sense", f"expected ccw or cw, got {sense!r}")
    return _FHumanThermal(cooling, transient, balls, window_s, sense)


@dataclass(frozen=True)
class _RodThermal:
    """A parsed rod four-port ``thermal:`` section (`WF-15` step 14, §9 items 831 / 884)."""

    cooling: RobinCooling
    transient: Any
    balls: dict
    drive: str


def _parse_thermal_rod(sec, source) -> _RodThermal:
    """``thermal:`` on a ``four_port_birdcage`` with ``implant:``: ``{k, h, t_inf,
    source_scale?, degree?, transient, balls: {mass_g, centres}, drive}``.

    The ``transient`` / ``balls`` blocks are ``_parse_thermal_f_human``'s, called
    (that path is not edited) after this section's own key checks; ``drive`` names
    one port at 1 V (``Study.temperature(drive=...)``).
    """
    p = "thermal"
    _unknown(sec, _THERMAL_REQUIRED + ("source_scale", "degree", "transient", "balls", "drive"),
             p, source)
    for key, what in (("transient", "a transient Study"), ("balls", "1 g ball centres"),
                      ("drive", "one port name")):
        if key not in sec:
            raise StudyConfigError(
                source, p, f"the implant-rod four_port_birdcage needs `thermal.{key}` ({what})"
            )
    drive = _string(sec["drive"], f"{p}.drive", source)
    if not drive:
        raise StudyConfigError(source, f"{p}.drive", "must name one port")
    parsed = _parse_thermal_f_human({k: v for k, v in sec.items() if k != "drive"}, source)
    return _RodThermal(parsed.cooling, parsed.transient, parsed.balls, drive)


@dataclass(frozen=True)
class _BallThermal:
    """A parsed ``thermal_ball`` ``thermal:`` section (`WF-1` step 14, §9 item 500)."""

    cooling: RobinCooling
    transient: Any
    duty: Optional[Any] = None


def _positive_finite(value, where, source) -> float:
    import math

    x = _number(value, where, source)
    if not (math.isfinite(x) and x > 0.0):
        raise StudyConfigError(source, where, f"must be finite and > 0, got {x!r}")
    return x


_SPHERE_DUTY_REFUSAL = (
    "duty is not available from YAML on the EM sphere (WF-1 step 29 parses transient "
    "only; the API's DutyCycle is WF-15 step 17)"
)

_SPHERE_TRANSIENT_DT_REFUSAL = (
    "the sphere transient takes t_end_s (dt = t_end_s / n_steps), not dt_s"
)


def _parse_thermal_sphere_transient(sec, source) -> _BallThermal:
    """A ``sphere_in_box`` ``thermal:`` with ``transient`` (`WF-1` step 29, §9 item 873).

    `MAT-12`'s SAR-driven transient on the EM sphere (U131's
    ``Study.temperature(cooling, transient=...)``).  ``duty`` and ``initial``
    are refused with the API's own reasons (``Study.temperature``); the
    ``transient`` block is ``_parse_thermal_box_mapping``'s (copied, gradient
    lines dropped), with ``t_end_s`` in place of ``dt_s`` (dt = t_end_s / n_steps).
    """
    import math

    p = "thermal"
    if "duty" in sec:
        # YAML-owned since §9 item 891: the API accepts duty on the sphere since
        # U135 (WF-15 step 17); no sphere ``thermal.duty`` key is commissioned.
        raise StudyConfigError(source, f"{p}.duty", _SPHERE_DUTY_REFUSAL)
    if "initial" in sec:
        raise StudyConfigError(
            source, f"{p}.initial",
            "initial_gradient_k_per_m is defined on a ThermalBoxGeometry only",
        )
    # Copied from _parse_thermal_box_mapping's transient block (gradient lines
    # dropped), not factored out of it (a shared helper would edit that path).
    tp = f"{p}.transient"
    tsec = sec["transient"]
    if not isinstance(tsec, dict):
        raise StudyConfigError(source, tp, f"expected a mapping, got {type(tsec).__name__}")
    t_req = ("t_end_s", "n_steps", "rho_c_j_per_m3_k")
    if "dt_s" in tsec:
        # §9 item 892 (18:00 ruling on 891's red): refused at the block, before
        # _unknown, so the landed s14 row 322 keeps its field.
        raise StudyConfigError(source, tp, _SPHERE_TRANSIENT_DT_REFUSAL)
    _unknown(tsec, t_req + ("theta", "startup_half_steps"), tp, source)
    _missing(tsec, t_req, tp, source)
    t_end = _number(tsec["t_end_s"], f"{tp}.t_end_s", source)
    n_steps = _integer(tsec["n_steps"], f"{tp}.n_steps", source)
    rho_c = _number(tsec["rho_c_j_per_m3_k"], f"{tp}.rho_c_j_per_m3_k", source)
    tkw = {}
    if "theta" in tsec:
        tkw["theta"] = _number(tsec["theta"], f"{tp}.theta", source)
    if "startup_half_steps" in tsec:
        tkw["startup_half_steps"] = _integer(
            tsec["startup_half_steps"], f"{tp}.startup_half_steps", source
        )
    if not (math.isfinite(t_end) and t_end > 0.0):
        raise StudyConfigError(source, f"{tp}.t_end_s", f"must be finite and > 0, got {t_end!r}")
    if n_steps < 1:
        raise StudyConfigError(source, f"{tp}.n_steps", f"must be >= 1, got {n_steps}")
    if not (math.isfinite(rho_c) and rho_c > 0.0):
        raise StudyConfigError(
            source, f"{tp}.rho_c_j_per_m3_k", f"must be finite and > 0, got {rho_c!r}"
        )
    if "theta" in tkw and not 0.0 <= tkw["theta"] <= 1.0:
        raise StudyConfigError(
            source, f"{tp}.theta", f"must be in [0, 1], got {tkw['theta']!r}"
        )
    if "startup_half_steps" in tkw and (
        tkw["startup_half_steps"] < 0 or tkw["startup_half_steps"] % 2
    ):
        raise StudyConfigError(
            source, f"{tp}.startup_half_steps",
            f"must be a non-negative even integer, got {tkw['startup_half_steps']}",
        )
    transient = TransientSpec(
        dt_s=t_end / n_steps, n_steps=n_steps, rho_c_j_per_m3_k=rho_c, **tkw
    )
    base = {k: sec[k] for k in sec if k != "transient"}
    cooling = _parse_thermal(base, "sphere_in_box", source)
    return _BallThermal(cooling, transient, None)


def _parse_thermal_ball_mapping(mapping, source):
    """A ``thermal_ball`` file -> ``(Study, "thermal_ball", outputs, _BallThermal)``.

    `WF-1` step 14 (§9 item 500).  Order: required sections, refused
    sections, geometry, excitation, thermal, outputs; first error raises.
    """
    import math

    for name in _BALL_REQUIRED:
        if name not in mapping:
            raise StudyConfigError(source, name, "required section is missing")
    for name, reason in _BALL_REFUSED.items():
        if name in mapping:
            raise StudyConfigError(source, name, reason)

    gp = "geometry"
    gsec = _section(mapping, "geometry", source)
    _unknown(gsec, ("kind", "radius", "resolution", "wall"), gp, source)
    _missing(gsec, ("radius", "resolution"), gp, source)
    radius = _positive_finite(gsec["radius"], f"{gp}.radius", source)
    resolution = _positive_finite(gsec["resolution"], f"{gp}.resolution", source)
    wall = None
    if "wall" in gsec:
        wp = f"{gp}.wall"
        wsec = gsec["wall"]
        if not isinstance(wsec, dict):
            raise StudyConfigError(source, wp, f"expected a mapping, got {type(wsec).__name__}")
        _unknown(wsec, _WALL_KEYS, wp, source)
        _missing(wsec, _WALL_KEYS, wp, source)
        wall = ContainerWall(
            **{k: _positive_finite(wsec[k], f"{wp}.{k}", source) for k in _WALL_KEYS}
        )
    geometry = (
        ThermalBallGeometry(radius, resolution)
        if wall is None
        else ThermalBallGeometry(radius, resolution, wall=wall)
    )

    ep = "excitation"
    esec = _section(mapping, "excitation", source)
    _kind(esec, ("uniform_heat_source",), ep, source)
    _unknown(esec, ("kind", "q_w_per_m3"), ep, source)
    _missing(esec, ("q_w_per_m3",), ep, source)
    q = _positive_finite(esec["q_w_per_m3"], f"{ep}.q_w_per_m3", source)

    p = "thermal"
    sec = _section(mapping, "thermal", source)
    _unknown(sec, _THERMAL_REQUIRED + ("degree", "source_scale", "transient", "duty"), p, source)
    if "source_scale" in sec:
        raise StudyConfigError(
            source, f"{p}.source_scale",
            "applies to an EM source; scale `excitation.q_w_per_m3` instead",
        )
    if "transient" not in sec:
        raise StudyConfigError(source, p, "a thermal_ball Study needs `thermal.transient`")
    base = {k: sec[k] for k in sec if k not in ("transient", "duty")}
    cooling = _parse_thermal(base, "sphere_in_box", source)

    # Copied from _parse_thermal_f_human's transient block, not factored out
    # of it (a shared helper would edit that path).
    tp = f"{p}.transient"
    tsec = sec["transient"]
    if not isinstance(tsec, dict):
        raise StudyConfigError(source, tp, f"expected a mapping, got {type(tsec).__name__}")
    t_req = ("dt_s", "n_steps", "rho_c_j_per_m3_k")
    _unknown(tsec, t_req + ("theta", "startup_half_steps"), tp, source)
    _missing(tsec, t_req, tp, source)
    dt = _number(tsec["dt_s"], f"{tp}.dt_s", source)
    n_steps = _integer(tsec["n_steps"], f"{tp}.n_steps", source)
    rho_c = _number(tsec["rho_c_j_per_m3_k"], f"{tp}.rho_c_j_per_m3_k", source)
    tkw = {}
    if "theta" in tsec:
        tkw["theta"] = _number(tsec["theta"], f"{tp}.theta", source)
    if "startup_half_steps" in tsec:
        tkw["startup_half_steps"] = _integer(
            tsec["startup_half_steps"], f"{tp}.startup_half_steps", source
        )
    if not (math.isfinite(dt) and dt > 0.0):
        raise StudyConfigError(source, f"{tp}.dt_s", f"must be finite and > 0, got {dt!r}")
    if n_steps < 1:
        raise StudyConfigError(source, f"{tp}.n_steps", f"must be >= 1, got {n_steps}")
    if not (math.isfinite(rho_c) and rho_c > 0.0):
        raise StudyConfigError(
            source, f"{tp}.rho_c_j_per_m3_k", f"must be finite and > 0, got {rho_c!r}"
        )
    if "theta" in tkw and not 0.0 <= tkw["theta"] <= 1.0:
        raise StudyConfigError(source, f"{tp}.theta", f"must be in [0, 1], got {tkw['theta']!r}")
    if "startup_half_steps" in tkw and (
        tkw["startup_half_steps"] < 0 or tkw["startup_half_steps"] % 2
    ):
        raise StudyConfigError(
            source, f"{tp}.startup_half_steps",
            f"must be a non-negative even integer, got {tkw['startup_half_steps']}",
        )
    transient = TransientSpec(dt_s=dt, n_steps=n_steps, rho_c_j_per_m3_k=rho_c, **tkw)

    duty = None
    if "duty" in sec:
        dp = f"{p}.duty"
        dsec = sec["duty"]
        if not isinstance(dsec, dict):
            raise StudyConfigError(source, dp, f"expected a mapping, got {type(dsec).__name__}")
        _unknown(dsec, ("period_s", "duty"), dp, source)
        _missing(dsec, ("period_s", "duty"), dp, source)
        period = _positive_finite(dsec["period_s"], f"{dp}.period_s", source)
        d = _number(dsec["duty"], f"{dp}.duty", source)
        if not 0.0 < d <= 1.0:
            raise StudyConfigError(source, f"{dp}.duty", f"must be in (0, 1], got {d!r}")
        duty = DutyCycle(period_s=period, duty=d)

    outputs = {}
    if "outputs" in mapping:
        osec = _section(mapping, "outputs", source)
        _unknown(osec, ("directory", "probes"), "outputs", source)
        if "directory" in osec:
            outputs["directory"] = _string(osec["directory"], "outputs.directory", source)
        if "probes" in osec:  # `WF-1` step 19 (§9 item 545)
            if wall is not None:  # review ruling 2026-10-19 03:00, Q3: no gel-plus-wall probe record
                raise StudyConfigError(
                    source, "outputs.probes", "probes are not admitted beside `geometry.wall`"
                )
            outputs["probes"] = _parse_box_probes(osec["probes"], source)

    study = Study(geometry, None, UniformHeatSource(q))
    return study, "thermal_ball", outputs, _BallThermal(cooling, transient, duty)


def _parse_thermal_box_mapping(mapping, source):
    """A ``thermal_box`` file -> ``(Study, "thermal_box", outputs, _BallThermal)``.

    `WF-1` step 15 (§9 item 521b).  ``_parse_thermal_ball_mapping``'s order
    and refusals; ``thermal.transient`` is optional (a steady box).
    """
    import math

    for name in _BALL_REQUIRED:
        if name not in mapping:
            raise StudyConfigError(source, name, "required section is missing")
    for name, reason in _BALL_REFUSED.items():
        if name in mapping:
            raise StudyConfigError(source, name, reason.replace("thermal_ball", "thermal_box"))

    gp = "geometry"
    gsec = _section(mapping, "geometry", source)
    _unknown(gsec, ("kind", "half_widths_m", "divisions", "wall"), gp, source)
    if "wall" in gsec:
        # `WF-1` step 28 (§9 item 830, review ruling 2026-11-05 03:00): refused at the
        # field before its keys are read when there is no transient (the pinned row,
        # test_wf1_yaml_thermal_box.py's variant table).
        tsec_raw = mapping.get("thermal")
        if not (isinstance(tsec_raw, dict) and "transient" in tsec_raw):
            raise StudyConfigError(
                source, f"{gp}.wall", "a walled box needs `thermal.transient` (no steady record)"
            )
    _missing(gsec, ("half_widths_m", "divisions"), gp, source)
    hp, dvp = f"{gp}.half_widths_m", f"{gp}.divisions"
    hraw, draw = gsec["half_widths_m"], gsec["divisions"]
    if not isinstance(hraw, (list, tuple)) or len(hraw) != 3:
        raise StudyConfigError(source, hp, f"expected three numbers, got {hraw!r}")
    half = tuple(_positive_finite(v, hp, source) for v in hraw)
    if not isinstance(draw, (list, tuple)) or len(draw) != 3:
        raise StudyConfigError(source, dvp, f"expected three integers, got {draw!r}")
    divisions = tuple(_integer(v, dvp, source) for v in draw)
    if any(n < 1 for n in divisions):
        raise StudyConfigError(source, dvp, f"must be integers >= 1, got {list(divisions)}")
    box_wall, box_nw = None, None
    if "wall" in gsec:  # `WF-1` step 28: _WALL_KEYS plus `divisions` (N_w cells across)
        wp = f"{gp}.wall"
        wsec = gsec["wall"]
        if not isinstance(wsec, dict):
            raise StudyConfigError(source, wp, f"expected a mapping, got {type(wsec).__name__}")
        _unknown(wsec, _WALL_KEYS + ("divisions",), wp, source)
        _missing(wsec, _WALL_KEYS + ("divisions",), wp, source)
        box_wall = ContainerWall(
            **{k: _positive_finite(wsec[k], f"{wp}.{k}", source) for k in _WALL_KEYS}
        )
        box_nw = _integer(wsec["divisions"], f"{wp}.divisions", source)
        if box_nw < 1:
            raise StudyConfigError(source, f"{wp}.divisions", f"must be >= 1, got {box_nw}")
    geometry = (
        ThermalBoxGeometry(half, divisions)
        if box_wall is None
        else ThermalBoxGeometry(half, divisions, wall=box_wall, wall_divisions=box_nw)
    )

    ep = "excitation"
    esec = _section(mapping, "excitation", source)
    _kind(esec, ("uniform_heat_source",), ep, source)
    _unknown(esec, ("kind", "q_w_per_m3", "source"), ep, source)
    _missing(esec, ("q_w_per_m3",), ep, source)
    q = _positive_finite(esec["q_w_per_m3"], f"{ep}.q_w_per_m3", source)
    # `WF-1` step 16 (§9 item 528): the optional `excitation.source` shape.
    source_kind = "uniform"
    if "source" in esec and box_wall is not None:  # `WF-1` step 28
        raise StudyConfigError(
            source, f"{ep}.source", "a source shape is not admitted beside `geometry.wall`"
        )
    if "source" in esec:
        sp = f"{ep}.source"
        ssec = esec["source"]
        if not isinstance(ssec, dict):
            raise StudyConfigError(source, sp, f"expected a mapping, got {type(ssec).__name__}")
        _unknown(ssec, ("kind",), sp, source)
        source_kind = _kind(ssec, ("uniform", "transverse_quadratic"), sp, source)
    heat_source = (
        TransverseQuadraticHeatSource(q) if source_kind == "transverse_quadratic"
        else UniformHeatSource(q)
    )

    p = "thermal"
    sec = _section(mapping, "thermal", source)
    _unknown(sec, _THERMAL_REQUIRED + ("degree", "source_scale", "transient", "duty", "initial",
                                       "dose"),
             p, source)
    if "source_scale" in sec:
        raise StudyConfigError(
            source, f"{p}.source_scale",
            "applies to an EM source; scale `excitation.q_w_per_m3` instead",
        )
    if "duty" in sec and "transient" not in sec:
        raise StudyConfigError(source, f"{p}.duty", "needs `thermal.transient`")
    if "dose" in sec:  # `WF-1` step 30 (§9 item 878)
        if "transient" not in sec:
            raise StudyConfigError(source, f"{p}.dose", "needs `thermal.transient`")
        if "duty" in sec:
            raise StudyConfigError(source, f"{p}.dose", "a dose is not admitted beside `thermal.duty`")
        if source_kind != "uniform":
            raise StudyConfigError(
                source, f"{p}.dose", "a dose is admitted with a uniform source only"
            )
    if "duty" in sec and box_wall is not None:  # `WF-1` step 28
        raise StudyConfigError(source, f"{p}.duty", "a duty is not admitted beside `geometry.wall`")
    # `WF-1` step 28 (§9 item 830): `thermal.initial: {gradient_k_per_m: g}`, θ₀ = g·z.
    gradient = None
    if "initial" in sec:
        ip = f"{p}.initial"
        if "transient" not in sec:
            raise StudyConfigError(source, ip, "needs `thermal.transient`")
        isec = sec["initial"]
        if not isinstance(isec, dict):
            raise StudyConfigError(source, ip, f"expected a mapping, got {type(isec).__name__}")
        _unknown(isec, ("gradient_k_per_m",), ip, source)
        _missing(isec, ("gradient_k_per_m",), ip, source)
        gradient = _number(isec["gradient_k_per_m"], f"{ip}.gradient_k_per_m", source)
        if not math.isfinite(gradient):
            raise StudyConfigError(
                source, f"{ip}.gradient_k_per_m", f"must be finite, got {gradient!r}"
            )
        if source_kind != "uniform":
            raise StudyConfigError(
                source, ip, "an initial gradient is admitted with a uniform source only"
            )
    base = {k: sec[k] for k in sec if k not in ("transient", "duty", "initial", "dose")}
    cooling = _parse_thermal(base, "sphere_in_box", source)

    # Copied from _parse_thermal_ball_mapping's transient / duty blocks, not
    # factored out of it (a shared helper would edit that path).
    transient = None
    if "transient" in sec:
        tp = f"{p}.transient"
        tsec = sec["transient"]
        if not isinstance(tsec, dict):
            raise StudyConfigError(source, tp, f"expected a mapping, got {type(tsec).__name__}")
        t_req = ("dt_s", "n_steps", "rho_c_j_per_m3_k")
        _unknown(tsec, t_req + ("theta", "startup_half_steps"), tp, source)
        _missing(tsec, t_req, tp, source)
        dt = _number(tsec["dt_s"], f"{tp}.dt_s", source)
        n_steps = _integer(tsec["n_steps"], f"{tp}.n_steps", source)
        rho_c = _number(tsec["rho_c_j_per_m3_k"], f"{tp}.rho_c_j_per_m3_k", source)
        tkw = {}
        if "theta" in tsec:
            tkw["theta"] = _number(tsec["theta"], f"{tp}.theta", source)
        if "startup_half_steps" in tsec:
            tkw["startup_half_steps"] = _integer(
                tsec["startup_half_steps"], f"{tp}.startup_half_steps", source
            )
        if not (math.isfinite(dt) and dt > 0.0):
            raise StudyConfigError(source, f"{tp}.dt_s", f"must be finite and > 0, got {dt!r}")
        if n_steps < 1:
            raise StudyConfigError(source, f"{tp}.n_steps", f"must be >= 1, got {n_steps}")
        if not (math.isfinite(rho_c) and rho_c > 0.0):
            raise StudyConfigError(
                source, f"{tp}.rho_c_j_per_m3_k", f"must be finite and > 0, got {rho_c!r}"
            )
        if "theta" in tkw and not 0.0 <= tkw["theta"] <= 1.0:
            raise StudyConfigError(
                source, f"{tp}.theta", f"must be in [0, 1], got {tkw['theta']!r}"
            )
        if "startup_half_steps" in tkw and (
            tkw["startup_half_steps"] < 0 or tkw["startup_half_steps"] % 2
        ):
            raise StudyConfigError(
                source, f"{tp}.startup_half_steps",
                f"must be a non-negative even integer, got {tkw['startup_half_steps']}",
            )
        if gradient is not None:
            tkw["initial_gradient_k_per_m"] = gradient
        transient = TransientSpec(dt_s=dt, n_steps=n_steps, rho_c_j_per_m3_k=rho_c, **tkw)

    duty = None
    if "duty" in sec:
        dp = f"{p}.duty"
        dsec = sec["duty"]
        if not isinstance(dsec, dict):
            raise StudyConfigError(source, dp, f"expected a mapping, got {type(dsec).__name__}")
        _unknown(dsec, ("period_s", "duty"), dp, source)
        _missing(dsec, ("period_s", "duty"), dp, source)
        period = _positive_finite(dsec["period_s"], f"{dp}.period_s", source)
        d = _number(dsec["duty"], f"{dp}.duty", source)
        if not 0.0 < d <= 1.0:
            raise StudyConfigError(source, f"{dp}.duty", f"must be in (0, 1], got {d!r}")
        duty = DutyCycle(period_s=period, duty=d)

    # `WF-1` step 30 (§9 item 878): `thermal.dose: {points_m: [[x, y, z], …], window_s: w}`,
    # carried in ``outputs["dose"]`` (as ``tune`` is); ``_BallThermal`` is unedited.
    dose = None
    if "dose" in sec:
        op = f"{p}.dose"
        osec_d = sec["dose"]
        if not isinstance(osec_d, dict):
            raise StudyConfigError(source, op, f"expected a mapping, got {type(osec_d).__name__}")
        _unknown(osec_d, ("points_m", "window_s"), op, source)
        _missing(osec_d, ("points_m",), op, source)
        pts = _triples(osec_d["points_m"], f"{op}.points_m", source)
        for i, pt in enumerate(pts):
            if not all(math.isfinite(c) for c in pt) or any(
                    abs(c) > h for c, h in zip(pt, half)):
                raise StudyConfigError(
                    source, f"{op}.points_m[{i}]",
                    f"must be a finite point inside the gel box [-a, a] x [-b, b] x [-c, c] "
                    f"= {half!r}, got {list(pt)!r}",
                )
        dkw = {}
        if "window_s" in osec_d:
            dkw["window_s"] = _positive_finite(osec_d["window_s"], f"{op}.window_s", source)
        dose = DoseReadout(tuple(tuple(pt) for pt in pts), **dkw)

    outputs = {}
    if "outputs" in mapping:
        osec = _section(mapping, "outputs", source)
        _unknown(osec, ("directory", "probes"), "outputs", source)
        if "directory" in osec:
            outputs["directory"] = _string(osec["directory"], "outputs.directory", source)
        if "probes" in osec:  # `WF-1` step 17 (§9 item 539)
            if dose is not None:  # `WF-1` step 30: a dose readout is not a probe
                raise StudyConfigError(
                    source, "outputs.probes", "probes are not admitted beside `thermal.dose`"
                )
            if box_wall is not None:  # `WF-1` step 28: no probe beside a box wall
                raise StudyConfigError(
                    source, "outputs.probes", "probes are not admitted beside `geometry.wall`"
                )
            outputs["probes"] = _parse_box_probes(osec["probes"], source)
    if dose is not None:
        outputs["dose"] = dose

    study = Study(geometry, None, heat_source)
    return study, "thermal_box", outputs, _BallThermal(cooling, transient, duty)


def _triple(value, where, source) -> tuple:
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise StudyConfigError(source, where, f"expected an [x, y, z] triple, got {value!r}")
    return tuple(_number(v, f"{where}[{i}]", source) for i, v in enumerate(value))


def _triples(value, where, source) -> list:
    if not isinstance(value, (list, tuple)) or len(value) == 0:
        raise StudyConfigError(source, where, "expected a non-empty list of [x, y, z] triples")
    return [_triple(v, f"{where}[{i}]", source) for i, v in enumerate(value)]


def _parse_sar(sec, source) -> dict:
    """``outputs.sar`` (`WF-1` step 5, §9 item 344): ``mass_g`` + ``centres``.

    ``centres`` is a list of ``[x, y, z]`` triples or ``{ring_radius_m: r}``,
    one centre per port at ``(r cos az, r sin az, 0)`` with ``az`` read off the
    solved Study's own sheets (an API convenience, ruling 2026-10-08 (3)).
    """
    p = "outputs.sar"
    if not isinstance(sec, dict):
        raise StudyConfigError(source, p, f"expected a mapping, got {type(sec).__name__}")
    _unknown(sec, ("mass_g", "centres"), p, source)
    _missing(sec, ("mass_g", "centres"), p, source)
    mass_g = _number(sec["mass_g"], f"{p}.mass_g", source)
    if not mass_g > 0.0:
        raise StudyConfigError(source, f"{p}.mass_g", f"must be > 0, got {mass_g!r}")
    centres = sec["centres"]
    if isinstance(centres, dict):
        _unknown(centres, ("ring_radius_m",), f"{p}.centres", source)
        _missing(centres, ("ring_radius_m",), f"{p}.centres", source)
        r = _number(centres["ring_radius_m"], f"{p}.centres.ring_radius_m", source)
        if not r > 0.0:
            raise StudyConfigError(source, f"{p}.centres.ring_radius_m", f"must be > 0, got {r!r}")
        return {"mass_g": mass_g, "centres": {"ring_radius_m": r}}
    return {"mass_g": mass_g, "centres": _triples(centres, f"{p}.centres", source)}


def _parse_b1_plus(sec, source) -> dict:
    """``outputs.b1_plus`` (`WF-1` step 5, §9 item 344): ``points``."""
    p = "outputs.b1_plus"
    if not isinstance(sec, dict):
        raise StudyConfigError(source, p, f"expected a mapping, got {type(sec).__name__}")
    _unknown(sec, ("points",), p, source)
    _missing(sec, ("points",), p, source)
    return {"points": _triples(sec["points"], f"{p}.points", source)}


def _parse_probes(value, geometry_kind, has_thermal, source) -> tuple:
    """``outputs.probes`` (`WF-1` step 6, §9 item 345) → a tuple of ``ProbeBall``.

    Needs a ``sphere_in_box`` and a ``thermal:`` section; every radius is
    checked here, before ``ProbeBall`` is built, so its ``ValueError`` never
    escapes bare.
    """
    from .thermal.probe import ProbeBall

    p = "outputs.probes"
    if geometry_kind != "sphere_in_box":
        raise StudyConfigError(
            source, p, f"probes are defined only on a sphere_in_box, not {geometry_kind!r}"
        )
    if not has_thermal:
        raise StudyConfigError(
            source, p, "probes read the steady temperature: a `thermal:` section is required"
        )
    if not isinstance(value, (list, tuple)) or len(value) == 0:
        raise StudyConfigError(source, p, "expected a non-empty list of {centre, radius} mappings")
    parsed = []
    for i, entry in enumerate(value):
        where = f"{p}[{i}]"
        if not isinstance(entry, dict):
            raise StudyConfigError(
                source, where, f"expected a mapping, got {type(entry).__name__}"
            )
        _unknown(entry, ("centre", "radius"), where, source)
        _missing(entry, ("centre", "radius"), where, source)
        centre = _triple(entry["centre"], f"{where}.centre", source)
        radius = _number(entry["radius"], f"{where}.radius", source)
        if not radius > 0.0:
            raise StudyConfigError(source, f"{where}.radius", f"must be > 0, got {radius!r}")
        parsed.append((centre, radius))
    return tuple(ProbeBall(c, r) for c, r in parsed)


def _parse_box_probes(value, source) -> tuple:
    """``outputs.probes`` on a ``thermal_box`` file (`WF-1` step 17, §9 item 539).

    ``_parse_probes``'s list-and-entry block, copied rather than factored (a
    shared helper would edit the sphere's path).  No geometry or ``thermal:``
    check: a box file always carries ``thermal:``.
    """
    from .thermal.probe import ProbeBall

    p = "outputs.probes"
    if not isinstance(value, (list, tuple)) or len(value) == 0:
        raise StudyConfigError(source, p, "expected a non-empty list of {centre, radius} mappings")
    parsed = []
    for i, entry in enumerate(value):
        where = f"{p}[{i}]"
        if not isinstance(entry, dict):
            raise StudyConfigError(
                source, where, f"expected a mapping, got {type(entry).__name__}"
            )
        _unknown(entry, ("centre", "radius"), where, source)
        _missing(entry, ("centre", "radius"), where, source)
        centre = _triple(entry["centre"], f"{where}.centre", source)
        radius = _number(entry["radius"], f"{where}.radius", source)
        if not radius > 0.0:
            raise StudyConfigError(source, f"{where}.radius", f"must be > 0, got {radius!r}")
        parsed.append((centre, radius))
    return tuple(ProbeBall(c, r) for c, r in parsed)


def _parse_em_probes(value, source) -> tuple:
    """``outputs.probes`` on a ``four_port_birdcage`` (`WF-1` step 20, §9 item 564).

    ``{points: [[x, y, z], …]}`` (the ``b1_plus`` shape): EM probes are
    points, thermal probes are balls.  The in-phantom check runs in
    ``_parse_mapping``, where the geometry is in hand.
    """
    p = "outputs.probes"
    if not isinstance(value, dict):
        raise StudyConfigError(
            source, p, f"expected a mapping {{points: [...]}}, got {type(value).__name__}"
        )
    _unknown(value, ("points",), p, source)
    _missing(value, ("points",), p, source)
    return tuple(_triples(value["points"], f"{p}.points", source))


def _parse_outputs(sec, geometry_kind, source, has_thermal=False) -> dict:
    p = "outputs"
    _unknown(sec, ("directory", "sar", "b1_plus", "probes"), p, source)
    out = {}
    if "directory" in sec:
        out["directory"] = _string(sec["directory"], f"{p}.directory", source)
    for key, parse in (("sar", _parse_sar), ("b1_plus", _parse_b1_plus)):
        if key in sec:
            if geometry_kind not in ("four_port_birdcage", "f_human_birdcage"):
                raise StudyConfigError(
                    source, f"{p}.{key}",
                    f"{key} is defined only on a four_port_birdcage, not {geometry_kind!r}",
                )
            out[key] = parse(sec[key], source)
            if (
                geometry_kind == "f_human_birdcage"
                and key == "sar"
                and isinstance(out[key]["centres"], dict)
            ):  # `WF-1` step 7a (§9 item 390a): that form reads four-port sheets
                raise StudyConfigError(
                    source,
                    f"{p}.sar.centres",
                    "the {ring_radius_m} form is not defined on an f_human_birdcage; "
                    "give a list of [x, y, z] centres",
                )
    if "probes" in sec:
        out["probes"] = _parse_em_probes(sec["probes"], source) if geometry_kind == "four_port_birdcage" else _parse_probes(sec["probes"], geometry_kind, has_thermal, source)  # `WF-1` step 20
    return out


def _parse_tune(sec, geometry_kind, source, excitation=None, has_load_result=False) -> dict:
    """``tune:`` → ``{}``, ``{"k": k}``, plus ``"at_hz"`` (`WF-1` step 11, §9 item 442).

    Admitted on a ``four_port_birdcage`` and, since `WF-1` step 12 (§9 item
    484), on an ``f_human_birdcage`` driven by ``lumped_port_drive``.  Since
    `WF-1` step 25 (§9 item 697) also on an ``f_clinical_birdcage`` (always
    ``lumped_port_drive``): without ``load_result:`` the run itself is refused.
    ``at_hz`` (`WF-1` step 27, §9 item 730) needs ``load_result:``; the file's
    own frequencies are checked in ``_run_load_result``.
    """
    p = "tune"
    if not (
        geometry_kind == "four_port_birdcage"
        or (geometry_kind in ("f_human_birdcage", "f_clinical_birdcage")
            and isinstance(excitation, LumpedPortDrive))
    ):
        raise StudyConfigError(
            source,
            p,
            "tune is defined on a four_port_birdcage, or an f_human_birdcage or "
            f"f_clinical_birdcage with excitation.kind lumped_port_drive, only, not {geometry_kind}",
        )
    _unknown(sec, ("k", "at_hz"), p, source)
    out = {}
    if "k" in sec:
        k = _integer(sec["k"], f"{p}.k", source)
        if k not in (1, 2):
            raise StudyConfigError(source, f"{p}.k", f"must be 1 or 2, got {k}")
        out["k"] = k
    if "at_hz" in sec:  # `WF-1` step 27 (§9 item 730)
        at_hz = _number(sec["at_hz"], f"{p}.at_hz", source)
        if not math.isfinite(at_hz) or at_hz <= 0.0:
            raise StudyConfigError(source, f"{p}.at_hz", f"must be finite and > 0, got {at_hz!r}")
        if not has_load_result:
            raise StudyConfigError(source, f"{p}.at_hz", "at_hz needs load_result")
        out["at_hz"] = at_hz
    return out


def _parse_sweep(sec, kind, geometry, excitation, conductor, outputs, has_tune, source) -> tuple:
    """``study.sweep_hz`` → a float tuple (`WF-1` step 21, §9 item 572).

    Admitted on the solid ``four_port_birdcage`` with ``lumped_port_drive``
    only: ``Study.solve_sweep`` has no conductor route, and no SAR / B₁⁺ /
    probe / tune output per frequency.  Refusals name ``study.sweep_hz``.
    """
    p = "study.sweep_hz"
    value = sec["sweep_hz"]
    if not isinstance(value, list):
        raise StudyConfigError(
            source, p, f"expected a list of frequencies in Hz, got {type(value).__name__} {value!r}"
        )
    freqs = []
    for i, v in enumerate(value):
        f = _number(v, f"{p}[{i}]", source)
        if not math.isfinite(f) or f <= 0.0:
            raise StudyConfigError(source, f"{p}[{i}]", f"must be finite and > 0, got {f!r}")
        freqs.append(f)
    if not freqs:
        raise StudyConfigError(source, p, "needs at least one frequency")
    for a, b in zip(freqs, freqs[1:]):
        if not b > a:
            raise StudyConfigError(
                source, p, f"must be strictly increasing (Touchstone needs sorted points): {a!r} then {b!r}"
            )
    if kind != "four_port_birdcage":
        raise StudyConfigError(source, p, f"a sweep is defined on a four_port_birdcage only, not {kind}")
    if not isinstance(excitation, LumpedPortDrive):
        raise StudyConfigError(
            source, p, "a sweep needs excitation.kind lumped_port_drive (Study.solve_sweep refuses others)"
        )
    if conductor is not None:
        raise StudyConfigError(
            source, p, "a sweep with a conductor is not supported (Study.solve_sweep has no conductor route)"
        )
    if getattr(geometry, "as_hole", False):
        raise StudyConfigError(source, p, "a sweep on geometry.as_hole: true is not supported")
    for key in ("sar", "b1_plus", "probes"):
        if key in outputs:
            raise StudyConfigError(source, p, f"a sweep has no outputs.{key} route")
    if has_tune:
        raise StudyConfigError(source, p, "a sweep has no tune route")
    return tuple(freqs)


def _refuse_with_load_result(mapping, source) -> None:
    """`WF-1` step 24 (§9 item 674): the fields a loaded S matrix cannot serve,
    refused by name before any other section is parsed (no mesh)."""
    why = "a load_result Study has an S matrix only (no mesh, no fields)"
    geo = mapping.get("geometry")
    f_human = isinstance(geo, dict) and geo.get("kind") == "f_human_birdcage"
    for name in ("conductor", "thermal", "implant"):
        if name == "conductor" and f_human:  # `WF-1` step 26 (§9 item 714): nothing is solved
            continue
        if name in mapping:
            raise StudyConfigError(source, name, f"{name} is not supported with load_result: {why}")
    if isinstance(mapping.get("study"), dict) and "sweep_hz" in mapping["study"]:
        raise StudyConfigError(
            source, "study.sweep_hz", f"study.sweep_hz is not supported with load_result: {why}"
        )
    outs = mapping.get("outputs")
    if isinstance(outs, dict):
        for key in ("sar", "b1_plus", "probes"):
            if key in outs:
                raise StudyConfigError(
                    source, f"outputs.{key}", f"outputs.{key} is not supported with load_result: {why}"
                )
    if "tune" not in mapping:
        raise StudyConfigError(source, "load_result", "load_result needs a tune section")


def _parse_load_result(value, kind, excitation, source) -> Path:
    """``load_result: <path>`` -> an absolute ``Path`` (`WF-1` step 24, §9 item 674).

    The file is not read here (``run_study_file`` reads it through
    ``Study.load_result``); a relative path resolves against the YAML file's
    directory."""
    p = "load_result"
    # `WF-1` step 25 (§9 item 697), step 26 (§9 item 714)
    if kind not in ("four_port_birdcage", "f_clinical_birdcage", "f_human_birdcage"):
        raise StudyConfigError(
            source, p, "load_result is defined on a four_port_birdcage, an f_clinical_birdcage "
            f"or an f_human_birdcage only, not {kind}"
        )
    if not isinstance(excitation, LumpedPortDrive):
        raise StudyConfigError(source, p, "load_result needs excitation.kind lumped_port_drive")
    path = Path(_string(value, p, source))
    if not path.is_absolute() and isinstance(source, Path):
        path = Path(source).resolve().parent / path
    return path


def _parse_mapping(mapping, source):
    """The private mapping entry point: ``(Study, kind, outputs, cooling)``.

    ``cooling`` is a ``RobinCooling`` or ``None`` (no ``thermal`` section);
    on an ``f_human_birdcage`` it is an ``_FHumanThermal`` (`WF-1` step 10a).
    """
    if not isinstance(mapping, dict):
        raise StudyConfigError(
            source, "<document>", f"expected a mapping at top level, got {type(mapping).__name__}"
        )
    if "femem" not in mapping:
        raise StudyConfigError(source, "femem", f"required key is missing (schema version {SCHEMA_VERSION})")
    version = _integer(mapping["femem"], "femem", source)
    if version != SCHEMA_VERSION:
        raise StudyConfigError(source, "femem", f"unsupported schema version {version} (this is {SCHEMA_VERSION})")
    for key in mapping:
        if key not in _TOP_KEYS:
            raise StudyConfigError(
                source, str(key), f"unknown top-level key (allowed: {', '.join(_TOP_KEYS)})"
            )
    geo = mapping.get("geometry")  # `WF-1` step 14 (§9 item 500)
    if "load_result" in mapping and isinstance(geo, dict) and geo.get("kind") in (
        "thermal_ball", "thermal_box"
    ):  # `WF-1` step 24 (§9 item 674)
        raise StudyConfigError(
            source, "load_result", "load_result is defined on a four_port_birdcage, an "
            f"f_clinical_birdcage or an f_human_birdcage only, not {geo.get('kind')}"
        )
    if isinstance(geo, dict) and geo.get("kind") == "thermal_ball":
        return _parse_thermal_ball_mapping(mapping, source)
    if isinstance(geo, dict) and geo.get("kind") == "thermal_box":  # `WF-1` step 15 (§9 item 521b)
        return _parse_thermal_box_mapping(mapping, source)
    for name in _REQUIRED_SECTIONS:
        if name not in mapping:
            raise StudyConfigError(source, name, "required section is missing")
    study_kwargs = _parse_study(_section(mapping, "study", source), source)
    if "load_result" in mapping:  # `WF-1` step 24 (§9 item 674): refused fields first
        _refuse_with_load_result(mapping, source)
    kind, geometry = _parse_geometry(_section(mapping, "geometry", source), source)
    if "implant" in mapping:  # `WF-15` step 10 (§9 item 709)
        geometry = _parse_implant(_section(mapping, "implant", source), kind, geometry, source)
    material = _parse_material(
        _section(mapping, "material", source), study_kwargs["frequency_hz"], source
    )
    excitation = _parse_excitation(
        _section(mapping, "excitation", source), kind, geometry, source
    )
    if "implant" in mapping:  # `WF-15` step 12 (§9 item 729): `Study`'s rod refusals, by field
        # `WF-15` step 13 (§9 item 802): admitted on the gated path only -- the PEC
        # four-port (no ``as_hole``, no ``conductor:``), at degree 1 (`study.degree` below).
        if isinstance(excitation, CapacitorSheets) and (
                getattr(geometry, "as_hole", False) or "conductor" in mapping):
            raise StudyConfigError(
                source, "implant", "implant with capacitor_sheets is not gated (WF-15 step 12)"
            )
        # Omitted means the Study's own default, degree 2 (see `_parse_study`).
        if study_kwargs.get("degree") != 1:
            raise StudyConfigError(
                source, "study.degree",
                f"implant needs degree 1, got {study_kwargs.get('degree', 'the default')} "
                "(WF-15 step 12)",
            )
    conductor = (
        _parse_conductor(_section(mapping, "conductor", source), kind, geometry, source)
        if "conductor" in mapping
        else None
    )
    if kind == "f_human_birdcage" and conductor is None:  # `WF-1` step 7a (§9 item 390a)
        raise StudyConfigError(
            source,
            "conductor",
            "required section is missing: two_port_quadrature needs a leontovich conductor",
        )
    cooling = (
        _parse_thermal(_section(mapping, "thermal", source), kind, source,
                       implant="implant" in mapping)
        if "thermal" in mapping
        else None
    )
    if isinstance(cooling, _RodThermal) and "sweep_hz" in mapping["study"]:  # `WF-15` step 14
        raise StudyConfigError(
            source, "thermal", "thermal is not defined beside study.sweep_hz (one frequency)"
        )
    ring = kind == "f_human_birdcage" and isinstance(excitation, LumpedPortDrive)  # `WF-1` step 12
    if ring and cooling is not None:  # E4: the per-sense routes need two_port_quadrature
        raise StudyConfigError(
            source, "thermal", "lumped_port_drive on an f_human_birdcage has no thermal route "
            "(use two_port_quadrature)"
        )
    outputs = (
        _parse_outputs(
            _section(mapping, "outputs", source), kind, source, has_thermal=cooling is not None
        )
        if "outputs" in mapping
        else {}
    )
    if isinstance(cooling, _BallThermal) and "probes" in outputs:  # `WF-1` step 29 (§9 item 873)
        raise StudyConfigError(
            source, "outputs.probes", "transient probes are not defined on the EM sphere "
            "(WF-15 step 15 gates temperature only)"
        )
    for key in ("sar", "b1_plus"):  # E4 (`WF-1` step 12, §9 item 484)
        if ring and key in outputs:
            raise StudyConfigError(
                source, f"outputs.{key}", "lumped_port_drive on an f_human_birdcage has no "
                f"{key} route (use two_port_quadrature)"
            )
    caps = isinstance(excitation, CapacitorSheets)  # `WF-1` step 13 (§9 item 496)
    for key in ("sar", "b1_plus"):
        if caps and key in outputs:
            raise StudyConfigError(
                source, f"outputs.{key}", f"capacitor_sheets has no {key} route "
                "(Study refuses a CapacitorSheets Study)"
            )
    if kind == "four_port_birdcage" and "probes" in outputs:  # `WF-1` step 20 (§9 item 564)
        if caps:
            raise StudyConfigError(
                source, "outputs.probes", "capacitor_sheets has no probes route "
                "(Study refuses a CapacitorSheets Study)"
            )
        for i, (x, y, z) in enumerate(outputs["probes"]):  # sigma and rho are the phantom's
            if not (math.hypot(x, y) < geometry.phantom_radius
                    and abs(z) < geometry.phantom_height / 2.0):
                raise StudyConfigError(
                    source, f"outputs.probes.points[{i}]",
                    f"({x!r}, {y!r}, {z!r}) is outside the phantom (r < "
                    f"{geometry.phantom_radius!r}, |z| < {geometry.phantom_height / 2.0!r})",
                )
    if caps and "tune" in mapping:
        raise StudyConfigError(
            source, "tune", "capacitor_sheets has no tune route (Study.tune refuses a "
            "CapacitorSheets Study; c_f is an input)"
        )
    if "tune" in mapping:  # `WF-1` step 11 (§9 item 442): rides in outputs (4-tuple kept)
        outputs = {
            **outputs,
            "tune": _parse_tune(
                _section(mapping, "tune", source), kind, source, excitation,
                has_load_result="load_result" in mapping,
            ),
        }
    if "sweep_hz" in mapping["study"]:  # `WF-1` step 21 (§9 item 572): rides in outputs
        outputs = {
            **outputs,
            "sweep_hz": _parse_sweep(
                mapping["study"], kind, geometry, excitation, conductor, outputs,
                "tune" in mapping, source,
            ),
        }
    if "load_result" in mapping:  # `WF-1` step 24 (§9 item 674): rides in outputs
        outputs = {
            **outputs,
            "load_result": _parse_load_result(mapping["load_result"], kind, excitation, source),
        }
    frequency_hz = study_kwargs.pop("frequency_hz")
    study = Study(
        geometry, material, excitation, frequency_hz=frequency_hz, conductor=conductor,
        **study_kwargs,
    )
    return study, kind, outputs, cooling


def _study_from_mapping(mapping, source="<mapping>") -> Study:
    return _parse_mapping(mapping, source)[0]


# ---------------------------------------------------------------------------
# YAML reading
# ---------------------------------------------------------------------------


def _read_mapping(path: Path):
    import yaml  # inside: `import fem_em_solver` must not need PyYAML

    class _StrictLoader(yaml.SafeLoader):
        """SafeLoader that refuses a duplicate mapping key (YAML 1.1 keeps
        the last one silently), naming the key by its dotted path."""

        def construct_mapping(self, node, deep=False):
            seen = set()
            for key_node, _value in node.value:
                key = self.construct_object(key_node, deep=True)
                if key in seen:
                    raise _DuplicateKey(node, key)
                seen.add(key)
            return super().construct_mapping(node, deep=deep)

    class _DuplicateKey(Exception):
        def __init__(self, node, key):
            self.node, self.key = node, key

    def _path_of(root, target, prefix=()):
        """Dotted path of mapping node ``target`` under ``root``."""
        if root is target:
            return prefix
        if isinstance(root, yaml.MappingNode):
            for k, v in root.value:
                found = _path_of(v, target, prefix + (str(k.value),))
                if found is not None:
                    return found
        elif isinstance(root, yaml.SequenceNode):
            for i, v in enumerate(root.value):
                found = _path_of(v, target, prefix + (f"[{i}]",))
                if found is not None:
                    return found
        return None

    text = path.read_text()
    loader = _StrictLoader(text)
    try:
        try:
            root = loader.get_single_node()
            if root is None:
                return None
            return loader.construct_document(root)
        except _DuplicateKey as dup:
            parts = _path_of(root, dup.node) or ()
            raise StudyConfigError(
                path, ".".join(parts + (str(dup.key),)), "duplicate key (YAML would keep the last silently)"
            ) from None
        except yaml.YAMLError as exc:
            raise StudyConfigError(path, "<yaml>", f"not valid YAML: {exc}") from None
    finally:
        loader.dispose()


def _load(path):
    path = Path(path)
    if not path.is_file():
        raise StudyConfigError(path, "<file>", "no such file")
    mapping = _read_mapping(path)
    return _parse_mapping(mapping, path)


def load_study(path) -> Study:
    """Parse a schema-v0 YAML study file into a ``Study`` (no mesh, no solve)."""
    return _load(path)[0]


# ---------------------------------------------------------------------------
# running
# ---------------------------------------------------------------------------


@dataclass
class StudyRun:
    """What ``run_study_file`` returns (not exported)."""

    study: Study
    result: Any
    summary: dict
    paths: dict = field(default_factory=dict)
    thermal: Any = None  # ThermalStudyResult when the file has `thermal:`
    probes: Any = None  # ProbeStudyResult when the file has `outputs.probes`


def _pairs(matrix):
    return [[[float(z.real), float(z.imag)] for z in row] for row in matrix]


def _temperature_summary(thermal, comm) -> dict:
    """Powers and θ max = max over owned dofs of Re(T − t_inf), ``MPI.MAX``-reduced.

    Taken in NumPy (ordering does not compile on complex UFL operands, `OPS-22`).
    """
    import numpy as np

    T = thermal.temperature
    n_own = T.function_space.dofmap.index_map.size_local * T.function_space.dofmap.index_map_bs
    theta = np.real(np.asarray(T.x.array[:n_own])) - float(thermal.cooling.t_inf_c)
    local = float(np.max(theta)) if theta.size else -np.inf
    return {
        "source_power_w": float(np.real(thermal.source_power)),
        "convected_power_w": float(np.real(thermal.convected_power)),
        "theta_max_k": float(comm.allreduce(local, op=MPI.MAX)),
    }


def _tuned_temperature_summary(thermal) -> dict:
    """The tuned F-human transient temperature as plain floats (`WF-1` step 10a)."""
    return {
        "sense": str(thermal.sense),
        "times_s": [float(t) for t in thermal.times_s],
        "volume_mean_k": [float(v) for v in thermal.volume_mean_k],
        "balls": [
            {
                "centre": [float(x) for x in b["centre"]],
                "radius_m": float(b["radius_m"]),
                "volume_m3": float(b["volume_m3"]),
                "theta_k": [float(v) for v in b["theta_k"]],
            }
            for b in thermal.balls
        ],
        "six_minute": {
            "window_s": float(thermal.six_minute["window_s"]),
            "volume_mean": {
                "peak_k": float(thermal.six_minute["volume_mean"]["peak_k"]),
                "t_max_s": float(thermal.six_minute["volume_mean"]["t_max_s"]),
            },
            "balls": [
                {"peak_k": float(b["peak_k"]), "t_max_s": float(b["t_max_s"])}
                for b in thermal.six_minute["balls"]
            ],
        },
        "source_power_w": float(thermal.source_power_w),
        "max_balance_residual": float(thermal.max_balance_residual),
        "source_scale": float(thermal.cooling.source_scale),
    }


def _sar_summary(study, spec, comm) -> dict:
    """``Study.sar`` at every centre for every drive (collective; `WF-1` step 5).

    ``{ring_radius_m: r}`` expands to one centre per sheet, in the solved
    sheets' order, ``(r cos az, r sin az, 0)`` with ``az`` in radians.
    """
    import numpy as np

    from .post.sar import averaging_ball_radius

    centres = spec["centres"]
    if isinstance(centres, dict):
        res = study._solve_port_fields(comm)  # cached: the same field set sar uses
        r = centres["ring_radius_m"]
        centres = [
            (
                r * float(np.cos(np.radians(float(s["azimuth_deg"])))),
                r * float(np.sin(np.radians(float(s["azimuth_deg"])))),
                0.0,
            )
            for s in res.sheets
        ]
    values = study.sar(mass_g=spec["mass_g"], centres=centres, comm=comm)
    return {
        "mass_g": float(spec["mass_g"]),
        "radius_m": float(
            averaging_ball_radius(mass_kg=spec["mass_g"] / 1000.0, rho=study.density_kg_m3)
        ),
        "centres": [[float(x) for x in c] for c in centres],
        "values": {
            pid: [{k: float(v) for k, v in d.items()} for d in per_centre]
            for pid, per_centre in values.items()
        },
    }


def _b1_plus_summary(study, spec, comm) -> dict:
    """``Study.b1_plus_at`` on the parsed points (collective; `WF-1` step 5)."""
    points = spec["points"]
    read = study.b1_plus_at(points, comm=comm)
    return {
        "points": [[float(x) for x in p] for p in points],
        "values": {pid: [float(v) for v in vals] for pid, (vals, _) in read.items()},
        "valid": {pid: [bool(b) for b in ok] for pid, (_, ok) in read.items()},
    }


def _em_probe_summary(study, points, comm) -> dict:
    """Point SAR / ``|E|`` per drive on the parsed points (collective; `WF-1` step 20).

    ``post.sar.point_sar`` under its owner-rule default on the cached kept
    fields, fed the split ``e_real`` / ``e_imag`` (never the complex field
    twice); ``|E| = sqrt(2 rho SAR / sigma)`` with the phantom's sigma and rho.
    """
    import numpy as np

    from .post.sar import point_sar

    res = study._solve_port_fields(comm)  # cached: the same field set sar uses
    pts = np.asarray(points, dtype=np.float64).reshape(-1, 3)
    sigma = float(study.material.sigma)
    rho = float(study.density_kg_m3)
    values = {}
    for pid, f in res.fields.items():
        sar = point_sar(f.e_real, f.e_imag, pts, sigma=sigma, rho=rho, comm=res.comm)
        e_abs = np.sqrt(2.0 * rho * sar / sigma)
        values[str(pid)] = {
            "sar_w_per_kg": [float(v) for v in sar],
            "e_abs_v_per_m": [float(v) for v in e_abs],
        }
    return {
        "points": [[float(x) for x in p] for p in points],
        "owner": "lowest_original_cell",
        "sigma_s_per_m": sigma,
        "rho_kg_m3": rho,
        "values": values,
    }


def _tune_summary(study, spec, comm) -> dict:
    """``Study.tune`` on the solved Study as plain floats (`WF-1` step 11)."""
    t = study.tune(**spec, comm=comm)
    return {
        "k": int(t.k),
        "frequency_hz": float(t.frequency_hz),
        "c_mode_farad": float(t.c_mode_farad),
        "z_k": [[float(complex(z).real), float(complex(z).imag)] for z in t.z_k],
        "off_diagonal_residual": float(t.off_diagonal_residual),
    }


def _run_thermal_ball(study, kind, cooling, path, out, digest, started, comm) -> StudyRun:
    """A ``thermal_ball`` file: ``Study.temperature`` only (`WF-1` step 14, §9 item 500)."""
    import numpy as np

    thermal = study.temperature(
        cooling.cooling, transient=cooling.transient, duty=cooling.duty, comm=comm
    )
    msh = thermal.mesh
    tr, duty, wall = cooling.transient, cooling.duty, study.geometry.wall
    # theta max: _temperature_summary's lines, copied (a transient result has
    # no convected_power).
    T = thermal.temperature
    n_own = T.function_space.dofmap.index_map.size_local * T.function_space.dofmap.index_map_bs
    theta = np.real(np.asarray(T.x.array[:n_own])) - float(cooling.cooling.t_inf_c)
    local = float(np.max(theta)) if theta.size else -np.inf
    summary = {
        "femem": SCHEMA_VERSION,
        "source": str(path),
        "sha256": digest,
        "family": kind,
        "ranks": int(comm.size),
        "cells": int(msh.topology.index_map(msh.topology.dim).size_global),
        "temperature": {
            "dt_s": float(tr.dt_s),
            "n_steps": int(tr.n_steps),
            "theta": float(tr.theta),
            "startup_half_steps": int(tr.startup_half_steps),
            "rho_c_j_per_m3_k": float(tr.rho_c_j_per_m3_k),
            "thermal_degree": int(cooling.cooling.degree),
            "q_w_per_m3": float(study.excitation.q_w_per_m3),
            "duty": None
            if duty is None
            else {"period_s": float(duty.period_s), "duty": float(duty.duty)},
            "wall": None
            if wall is None
            else {
                "thickness_m": float(wall.thickness_m),
                "k_w_per_m_k": float(wall.k_w_per_m_k),
                "rho_c_j_per_m3_k": float(wall.rho_c_j_per_m3_k),
            },
            "source_multipliers": [float(m) for m in thermal.source_multipliers],
            "source_power_w": float(np.real(thermal.source_power)),
            "max_balance_residual": float(thermal.max_balance_residual),
            "theta_max_k": float(comm.allreduce(local, op=MPI.MAX)),
        },
    }
    summary["elapsed_s"] = float(time.perf_counter() - started)
    summary_path = out / "summary.json"
    if comm.rank == 0:
        out.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as fh:
            json.dump(summary, fh, indent=1)
    comm.Barrier()
    return StudyRun(
        study=study, result=None, summary=summary, paths={"summary": summary_path},
        thermal=thermal,
    )


def _run_rod_thermal(study, kind, cooling, path, out, digest, started, comm) -> StudyRun:
    """The implant-rod four-port ``thermal:`` (`WF-15` step 14, §9 items 831 / 884):
    ``Study.temperature(drive=...)`` only -- one port-fields solve, no S summary --
    writing the per-ball θ series and the energy balance."""
    thermal = study.temperature(
        cooling.cooling, transient=cooling.transient, balls=cooling.balls,
        drive=cooling.drive, comm=comm,
    )
    tr = cooling.transient
    summary = {
        "femem": SCHEMA_VERSION,
        "source": str(path),
        "sha256": digest,
        "family": kind,
        "degree": int(study.degree),
        "frequency_hz": float(study.frequency_hz),
        "ranks": int(comm.size),
        "cells": int(thermal.cells),
        "temperature": {
            "drive": str(thermal.sense),
            "dt_s": float(tr.dt_s),
            "n_steps": int(tr.n_steps),
            "theta": float(tr.theta),
            "startup_half_steps": int(tr.startup_half_steps),
            "rho_c_j_per_m3_k": float(tr.rho_c_j_per_m3_k),
            "thermal_degree": int(cooling.cooling.degree),
            "source_scale": float(cooling.cooling.source_scale),
            "times_s": [float(t) for t in thermal.times_s],
            "volume_mean_k": [float(v) for v in thermal.volume_mean_k],
            "balls": [
                {
                    "centre": [float(x) for x in b["centre"]],
                    "radius_m": float(b["radius_m"]),
                    "volume_m3": float(b["volume_m3"]),
                    "theta_k": [float(v) for v in b["theta_k"]],
                }
                for b in thermal.balls
            ],
            "source_power_w": float(thermal.source_power_w),
            "max_balance_residual": float(thermal.max_balance_residual),
        },
    }
    summary["elapsed_s"] = float(time.perf_counter() - started)
    summary_path = out / "summary.json"
    if comm.rank == 0:
        out.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as fh:
            json.dump(summary, fh, indent=1)
    comm.Barrier()
    return StudyRun(
        study=study, result=None, summary=summary, paths={"summary": summary_path},
        thermal=thermal,
    )


def _run_thermal_box(study, kind, cooling, path, out, digest, started, comm) -> StudyRun:
    """A ``thermal_box`` file: ``Study.temperature`` only (`WF-1` step 15, §9 item 521b)."""
    import numpy as np

    thermal = study.temperature(
        cooling.cooling, transient=cooling.transient, duty=cooling.duty, comm=comm
    )
    tr, duty = cooling.transient, cooling.duty
    msh = thermal.submesh if tr is None else thermal.mesh
    # theta max: _temperature_summary's lines, copied.
    T = thermal.temperature
    n_own = T.function_space.dofmap.index_map.size_local * T.function_space.dofmap.index_map_bs
    theta = np.real(np.asarray(T.x.array[:n_own])) - float(cooling.cooling.t_inf_c)
    local = float(np.max(theta)) if theta.size else -np.inf
    temp = {
        "half_widths_m": [float(v) for v in study.geometry.half_widths_m],
        "divisions": [int(v) for v in study.geometry.divisions],
        "thermal_degree": int(cooling.cooling.degree),
        "q_w_per_m3": float(study.excitation.q_w_per_m3),
        "transient": None
        if tr is None
        else {
            "dt_s": float(tr.dt_s),
            "n_steps": int(tr.n_steps),
            "theta": float(tr.theta),
            "startup_half_steps": int(tr.startup_half_steps),
            "rho_c_j_per_m3_k": float(tr.rho_c_j_per_m3_k),
        },
        "duty": None
        if duty is None
        else {"period_s": float(duty.period_s), "duty": float(duty.duty)},
        "source_power_w": float(np.real(thermal.source_power)),
    }
    if tr is None:
        temp["convected_power_w"] = float(np.real(thermal.convected_power))
    else:
        temp["source_multipliers"] = [float(m) for m in thermal.source_multipliers]
        temp["max_balance_residual"] = float(thermal.max_balance_residual)
    temp["theta_max_k"] = float(comm.allreduce(local, op=MPI.MAX))
    temp["source"] = (  # `WF-1` step 16 (§9 item 528)
        "transverse_quadratic" if isinstance(study.excitation, TransverseQuadraticHeatSource)
        else "uniform"
    )
    w = study.geometry.wall  # `WF-1` step 28 (§9 item 830): keys only when set
    if w is not None:
        temp["wall"] = {
            "thickness_m": float(w.thickness_m),
            "k_w_per_m_k": float(w.k_w_per_m_k),
            "rho_c_j_per_m3_k": float(w.rho_c_j_per_m3_k),
            "divisions": int(study.geometry.wall_divisions),
        }
    if tr is not None and tr.initial_gradient_k_per_m is not None:
        temp["initial_gradient_k_per_m"] = float(tr.initial_gradient_k_per_m)
    summary = {
        "femem": SCHEMA_VERSION,
        "source": str(path),
        "sha256": digest,
        "family": kind,
        "ranks": int(comm.size),
        "cells": int(msh.topology.index_map(msh.topology.dim).size_global),
        "temperature": temp,
    }
    summary["elapsed_s"] = float(time.perf_counter() - started)
    summary_path = out / "summary.json"
    if comm.rank == 0:
        out.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as fh:
            json.dump(summary, fh, indent=1)
    comm.Barrier()
    return StudyRun(
        study=study, result=None, summary=summary, paths={"summary": summary_path},
        thermal=thermal,
    )


def _run_thermal_box_dose(study, kind, cooling, outputs, path, out, digest, started,
                          comm) -> StudyRun:
    """A ``thermal_box`` file with ``thermal.dose`` (`WF-1` step 30, §9 item 878).

    ``_run_thermal_box``'s lines (`WF-1` step 15's span, copied, that function
    unedited) with ``dose=outputs["dose"]`` passed to ``Study.temperature`` and
    ``summary["dose"]`` = the result's ``dose`` (`MAT-20`'s ``cem43`` and
    trailing-window mean per point; a readout, no limit implied).
    """
    import numpy as np

    thermal = study.temperature(
        cooling.cooling, transient=cooling.transient, duty=cooling.duty, dose=outputs["dose"],
        comm=comm,
    )
    tr, duty = cooling.transient, cooling.duty
    msh = thermal.submesh if tr is None else thermal.mesh
    # theta max: _temperature_summary's lines, copied.
    T = thermal.temperature
    n_own = T.function_space.dofmap.index_map.size_local * T.function_space.dofmap.index_map_bs
    theta = np.real(np.asarray(T.x.array[:n_own])) - float(cooling.cooling.t_inf_c)
    local = float(np.max(theta)) if theta.size else -np.inf
    temp = {
        "half_widths_m": [float(v) for v in study.geometry.half_widths_m],
        "divisions": [int(v) for v in study.geometry.divisions],
        "thermal_degree": int(cooling.cooling.degree),
        "q_w_per_m3": float(study.excitation.q_w_per_m3),
        "transient": None
        if tr is None
        else {
            "dt_s": float(tr.dt_s),
            "n_steps": int(tr.n_steps),
            "theta": float(tr.theta),
            "startup_half_steps": int(tr.startup_half_steps),
            "rho_c_j_per_m3_k": float(tr.rho_c_j_per_m3_k),
        },
        "duty": None
        if duty is None
        else {"period_s": float(duty.period_s), "duty": float(duty.duty)},
        "source_power_w": float(np.real(thermal.source_power)),
    }
    if tr is None:
        temp["convected_power_w"] = float(np.real(thermal.convected_power))
    else:
        temp["source_multipliers"] = [float(m) for m in thermal.source_multipliers]
        temp["max_balance_residual"] = float(thermal.max_balance_residual)
    temp["theta_max_k"] = float(comm.allreduce(local, op=MPI.MAX))
    temp["source"] = (  # `WF-1` step 16 (§9 item 528)
        "transverse_quadratic" if isinstance(study.excitation, TransverseQuadraticHeatSource)
        else "uniform"
    )
    w = study.geometry.wall  # `WF-1` step 28 (§9 item 830): keys only when set
    if w is not None:
        temp["wall"] = {
            "thickness_m": float(w.thickness_m),
            "k_w_per_m_k": float(w.k_w_per_m_k),
            "rho_c_j_per_m3_k": float(w.rho_c_j_per_m3_k),
            "divisions": int(study.geometry.wall_divisions),
        }
    if tr is not None and tr.initial_gradient_k_per_m is not None:
        temp["initial_gradient_k_per_m"] = float(tr.initial_gradient_k_per_m)
    summary = {
        "femem": SCHEMA_VERSION,
        "source": str(path),
        "sha256": digest,
        "family": kind,
        "ranks": int(comm.size),
        "cells": int(msh.topology.index_map(msh.topology.dim).size_global),
        "temperature": temp,
    }
    summary["dose"] = thermal.dose
    summary["elapsed_s"] = float(time.perf_counter() - started)
    summary_path = out / "summary.json"
    if comm.rank == 0:
        out.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as fh:
            json.dump(summary, fh, indent=1)
    comm.Barrier()
    return StudyRun(
        study=study, result=None, summary=summary, paths={"summary": summary_path},
        thermal=thermal,
    )


def _run_thermal_box_probes(study, kind, cooling, outputs, path, out, digest, started,
                            comm) -> StudyRun:
    """A ``thermal_box`` file with ``outputs.probes`` (`WF-1` step 17, §9 item 539).

    ``_run_thermal_box`` unedited, then ``Study.probe`` (collective; it
    re-solves the temperature, as the sphere route does).  ``summary.json`` is
    written twice; the second write, with ``probes``, is the record.
    """
    run = _run_thermal_box(study, kind, cooling, path, out, digest, started, comm)
    probes = study.probe(
        outputs["probes"], cooling.cooling, transient=cooling.transient, duty=cooling.duty
    )  # collective: every rank enters, owner or not
    summary = run.summary
    summary["probes"] = {
        "centres": [[float(x) for x in pb.centre] for pb in probes.probes],
        "radii": [float(pb.radius) for pb in probes.probes],
        "values_k": [float(v) for v in probes.values],
        "volumes_m3": [float(v) for v in probes.volumes],
    }
    summary["elapsed_s"] = float(time.perf_counter() - started)
    if comm.rank == 0:
        out.mkdir(parents=True, exist_ok=True)
        with open(run.paths["summary"], "w") as fh:
            json.dump(summary, fh, indent=1)
    comm.Barrier()
    run.probes = probes
    return run


def _run_thermal_ball_probes(study, kind, cooling, outputs, path, out, digest, started,
                             comm) -> StudyRun:
    """A ``thermal_ball`` file with ``outputs.probes`` (`WF-1` step 19, §9 item 545).

    ``_run_thermal_box_probes`` copied with ``_run_thermal_ball`` in place of
    ``_run_thermal_box``: ``_run_thermal_ball`` unedited, then ``Study.probe``
    (collective; it re-solves the temperature).  ``summary.json`` is written
    twice; the second write, with ``probes``, is the record.
    """
    run = _run_thermal_ball(study, kind, cooling, path, out, digest, started, comm)
    probes = study.probe(
        outputs["probes"], cooling.cooling, transient=cooling.transient, duty=cooling.duty
    )  # collective: every rank enters, owner or not
    summary = run.summary
    summary["probes"] = {
        "centres": [[float(x) for x in pb.centre] for pb in probes.probes],
        "radii": [float(pb.radius) for pb in probes.probes],
        "values_k": [float(v) for v in probes.values],
        "volumes_m3": [float(v) for v in probes.volumes],
    }
    summary["elapsed_s"] = float(time.perf_counter() - started)
    if comm.rank == 0:
        out.mkdir(parents=True, exist_ok=True)
        with open(run.paths["summary"], "w") as fh:
            json.dump(summary, fh, indent=1)
    comm.Barrier()
    run.probes = probes
    return run


# `ARCH-1` step 1d (§9 item 590, U76): the only ``outputs`` key the Palace
# route honours; every other key is refused by name (an allow-list, so a key
# added later is refused rather than silently dropped).
_PALACE_ROUTED_OUTPUTS = ("directory",)
# Display names for keys that are not ``outputs.*`` in the YAML.
_PALACE_OUTPUT_FIELD = {"tune": "tune", "sweep_hz": "study.sweep_hz", "load_result": "load_result"}


def _palace_input_refusal(study, path, comm):
    """`ARCH-1` step 1e (§9 item 601, U79): the two inputs ``solve_palace``
    would otherwise reject by ``ValueError`` traceback, refused by name.

    Reads only the geometry and ``comm.size`` (no collective), so every rank
    raises before any mesh.  The degree needs no check here: the loader
    already refuses anything but 1 or 2 at ``study.degree``.
    """
    if study.geometry.as_hole:
        raise StudyConfigError(path, "--backend", "palace: geometry.as_hole is not supported")
    if comm.size != 1:
        raise StudyConfigError(
            path, "--backend",
            f"palace: a {comm.size}-rank launch is not supported "
            "(the case writer is serial; run at mpiexec -n 1)")


def _run_palace(study, kind, outputs, path, out, digest, started, comm,
                *, cooling=None) -> StudyRun:
    """`ARCH-1` step 1a (§9 item 517): write / read a Palace case; S only."""
    if kind != "four_port_birdcage":
        raise StudyConfigError(path, "--backend", f"palace supports four_port_birdcage only, got {kind}")
    if cooling is not None:  # `ARCH-1` step 1d (§9 item 590)
        raise StudyConfigError(path, "--backend", "palace: thermal is not supported")
    for key in sorted(outputs):  # `ARCH-1` step 1d (§9 item 590): allow-list
        if key not in _PALACE_ROUTED_OUTPUTS:
            raise StudyConfigError(
                path, "--backend",
                f"palace: {_PALACE_OUTPUT_FIELD.get(key, 'outputs.' + key)} is not supported")
    if study.conductor is not None:
        raise StudyConfigError(path, "--backend", "palace: a conductor is not supported")
    if isinstance(study.excitation, CapacitorSheets):
        raise StudyConfigError(path, "--backend", "palace: capacitor_sheets is not supported")
    _palace_input_refusal(study, path, comm)
    result = study.solve(comm, backend="palace", case_dir=out / "palace_case")
    summary = {
        "femem": SCHEMA_VERSION,
        "source": str(path),
        "sha256": digest,
        "family": kind,
        "degree": int(study.degree),
        "frequency_hz": float(study.frequency_hz),
        "ranks": int(comm.size),
        "backend": "palace",
        "cells": int(result.cells),
        "palace": {
            "case_dir": str(result.case_dir),
            "order": int(result.order),
            "sha256": dict(result.sha256),
            "palace_total_s": result.palace_total_s,
        },
        "s_matrix": _pairs(result.s_matrix),
    }
    summary["elapsed_s"] = float(time.perf_counter() - started)
    summary_path = out / "summary.json"
    out.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=1)
    return StudyRun(study=study, result=result, summary=summary,
                    paths={"summary": summary_path})


def _run_four_port_sweep(study, kind, outputs, path, out, digest, started, comm) -> StudyRun:
    """`WF-1` step 21 (§9 item 572): ``Study.solve_sweep`` on one mesh, one ``.s4p``."""
    sw = study.solve_sweep(outputs["sweep_hz"], comm)
    touchstone = out / f"{path.stem}.s4p"
    sw.write_touchstone(touchstone)  # collective: rank 0 writes, Barrier
    summary = {
        "femem": SCHEMA_VERSION,
        "source": str(path),
        "sha256": digest,
        "family": kind,
        "degree": int(study.degree),
        "frequency_hz": float(study.frequency_hz),
        "sweep_hz": [float(f) for f in sw.frequencies_hz],
        "ranks": int(comm.size),
        "cells": int(sw.results[0].cells),
        "s_matrices": [_pairs(r.s_matrix) for r in sw.results],
        "z_matrices": [_pairs(r.z_matrix) for r in sw.results],
        "touchstone": touchstone.name,
        "mesh_time_s": float(sw.mesh_time),
        "sweep_times_s": [float(t) for t in sw.sweep_times],
    }
    summary["elapsed_s"] = float(time.perf_counter() - started)
    summary_path = out / "summary.json"
    if comm.rank == 0:
        out.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as fh:
            json.dump(summary, fh, indent=1)
    comm.Barrier()
    return StudyRun(study=study, result=sw, summary=summary,
                    paths={"summary": summary_path, "touchstone": touchstone})


def _run_load_result(study, kind, outputs, path, out, digest, started, comm) -> StudyRun:
    """`WF-1` step 24 (§9 item 674): ``Study.load_result`` then ``Study.tune``.

    No mesh and no solve; the ``tune`` block is ``_tune_summary``'s, as on the
    solve route.  No Touchstone is written (the input is one)."""
    loaded = study.load_result(outputs["load_result"])
    if len(loaded.frequencies_hz) >= 2 and "at_hz" not in outputs["tune"]:  # `OPS-76`, item 754
        raise StudyConfigError(path, "tune.at_hz", "a multi-frequency load_result needs at_hz")
    if "at_hz" in outputs["tune"]:  # `WF-1` step 27 (§9 item 730): every rank read the file
        at_hz, freqs = outputs["tune"]["at_hz"], loaded.frequencies_hz
        if len(freqs) < 2:
            raise StudyConfigError(path, "tune.at_hz", "at_hz needs a multi-frequency load_result")
        if not min(freqs) <= at_hz <= max(freqs):
            raise StudyConfigError(path, "tune.at_hz", f"{at_hz!r} Hz is outside the file's "
                                   f"span [{min(freqs)!r}, {max(freqs)!r}] Hz")
    summary = {
        "femem": SCHEMA_VERSION,
        "source": str(path),
        "sha256": digest,
        "family": kind,
        "degree": int(study.degree),
        "frequency_hz": float(study.frequency_hz),
        "ranks": int(comm.size),
        "load_result": str(outputs["load_result"]),
        "port_ids": [str(p) for p in loaded.port_ids],
        "s_matrix": _pairs(loaded.s_matrix),
        "tune": _tune_summary(study, outputs["tune"], comm),
    }
    summary["elapsed_s"] = float(time.perf_counter() - started)
    summary_path = out / "summary.json"
    if comm.rank == 0:
        out.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as fh:
            json.dump(summary, fh, indent=1)
    comm.Barrier()
    return StudyRun(study=study, result=loaded, summary=summary, paths={"summary": summary_path})


def _run_study_file_palace(path, *, output_dir=None,
                           comm: MPI.Intracomm = MPI.COMM_WORLD) -> StudyRun:
    """``femem run --backend palace``: ``run_study_file``'s prelude, then the
    Palace write / read stages. Private, so the generated API page of
    ``run_study_file`` is unchanged (`WF-17` (G))."""
    path = Path(path).resolve()
    started = time.perf_counter()
    study, kind, outputs, cooling = _load(path)
    if kind == "f_clinical_birdcage":  # `WF-13` step 12 (§9 item 663): before any mesh
        raise StudyConfigError(
            path, "geometry.kind", f"f_clinical_birdcage solve is {_f_clinical_tier(study.degree)} "
            "(use `femem run --estimate`)"
        )
    out = Path(
        output_dir
        if output_dir is not None
        else outputs.get("directory", os.path.join("femem_out", path.stem))
    )
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return _run_palace(study, kind, outputs, path, out, digest, started, comm,
                       cooling=cooling)


def run_study_file(path, *, output_dir=None, comm: MPI.Intracomm = MPI.COMM_WORLD) -> StudyRun:
    """``load_study``, ``Study.solve``, then the outputs.

    Output directory: the argument, else ``outputs.directory``, else
    ``./femem_out/<stem>``.  Rank 0 writes the JSON summary into that
    directory after every reduction, then ``Barrier``.  A port study also
    writes ``<stem>.s<N>p``.  An ``f_human_birdcage`` ``two_port_quadrature``
    study (`WF-1` step 7a, §9 item 390a) writes no Touchstone: its node S is
    not a reference-impedance network, so the summary carries ``s_node`` /
    ``y_in`` and no ``z_matrix``.  An ``f_human_birdcage`` ``lumped_port_drive``
    study (`WF-1` step 12, §9 item 484) writes ``<stem>.s32p``.  A
    ``capacitor_sheets`` study (`WF-1` step 13, §9 item 496) writes no
    Touchstone: the summary carries the reduced ``s_matrix`` over
    ``kept_ports`` and a ``capacitors`` block, and no ``z_matrix``.
    A ``thermal_ball`` study (`WF-1` step 14, §9 item 500) runs
    ``Study.temperature`` with no field solve and writes no Touchstone.
    A ``thermal_box`` study (`WF-1` step 15, §9 item 521b) does the same, steady or transient.
    A ``thermal_box`` study with ``outputs.probes`` (`WF-1` step 17, §9 item
    539) then runs ``Study.probe`` and adds ``probes`` to the summary; so
    does a ``thermal_ball`` study with ``outputs.probes`` (`WF-1` step 19, §9 item 545).
    A ``four_port_birdcage`` study with ``outputs.probes`` (`WF-1` step 20, §9
    item 564) adds point SAR / ``|E|`` per drive to the summary as ``probes``.
    A ``four_port_birdcage`` study with ``study.sweep_hz`` (`WF-1` step 21, §9
    item 572) runs ``Study.solve_sweep`` instead: one mesh, one ``<stem>.s4p``.
    """
    path = Path(path).resolve()
    started = time.perf_counter()
    study, kind, outputs, cooling = _load(path)
    # `WF-13` step 12 (§9 item 663): before any mesh; a saved matrix needs no solve
    # (`WF-1` step 25, §9 item 697).
    if kind == "f_clinical_birdcage" and "load_result" not in outputs:
        raise StudyConfigError(
            path, "geometry.kind", f"f_clinical_birdcage solve is {_f_clinical_tier(study.degree)} "
            "(use `femem run --estimate`)"
        )
    out = Path(
        output_dir
        if output_dir is not None
        else outputs.get("directory", os.path.join("femem_out", path.stem))
    )
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if isinstance(study.geometry, ThermalBoxGeometry) and "probes" in outputs:  # §9 item 539
        return _run_thermal_box_probes(
            study, kind, cooling, outputs, path, out, digest, started, comm
        )
    if isinstance(study.geometry, ThermalBoxGeometry) and "dose" in outputs:  # §9 item 878
        return _run_thermal_box_dose(
            study, kind, cooling, outputs, path, out, digest, started, comm
        )
    if isinstance(study.geometry, ThermalBoxGeometry):  # `WF-1` step 15 (§9 item 521b)
        return _run_thermal_box(study, kind, cooling, path, out, digest, started, comm)
    if isinstance(study.geometry, ThermalBallGeometry) and "probes" in outputs:  # §9 item 545
        return _run_thermal_ball_probes(
            study, kind, cooling, outputs, path, out, digest, started, comm
        )
    if isinstance(study.geometry, ThermalBallGeometry):  # `WF-1` step 14 (§9 item 500)
        return _run_thermal_ball(study, kind, cooling, path, out, digest, started, comm)
    if isinstance(cooling, _RodThermal):  # `WF-15` step 14 (§9 items 831 / 884)
        return _run_rod_thermal(study, kind, cooling, path, out, digest, started, comm)
    if "sweep_hz" in outputs:  # WF-1 step 21 (§9 item 572)
        return _run_four_port_sweep(study, kind, outputs, path, out, digest, started, comm)
    if "load_result" in outputs:  # `WF-1` step 24 (§9 item 674): no mesh, no solve
        return _run_load_result(study, kind, outputs, path, out, digest, started, comm)

    result = study.solve(comm)
    summary = {
        "femem": SCHEMA_VERSION,
        "source": str(path),
        "sha256": digest,
        "family": kind,
        "degree": int(study.degree),
        "frequency_hz": float(study.frequency_hz),
        "ranks": int(comm.size),
    }
    if study.conductor is not None:  # `WF-1` step 4: absent otherwise
        summary["conductor"] = {
            "kind": "leontovich",
            "sigma_s_per_m": float(study.conductor.sigma_s_per_m),
        }
    paths = {}
    thermal = None
    probes = None
    if isinstance(study.excitation, CapacitorSheets):  # `WF-1` step 13 (§9 item 496)
        zc = complex(result.capacitor_impedance_ohm)
        summary["cells"] = int(result.cells)
        summary["s_matrix"] = _pairs(result.s_matrix)
        summary["kept_ports"] = [str(pd.port_id) for pd in result.port_definitions]
        summary["capacitors"] = {
            "sheets": [str(s) for s in result.capacitor_port_ids],
            "c_f": float(study.excitation.c_f),
            "esr_ohm": float(study.excitation.esr_ohm),
            "impedance_ohm": [float(zc.real), float(zc.imag)],
        }
    elif isinstance(study.geometry, FourPortBirdcageGeometry):
        summary["cells"] = int(result.cells)
        n = int(result.s_matrix.shape[0])
        touchstone = out / f"{path.stem}.s{n}p"
        result.write_touchstone(touchstone)  # collective: rank 0 writes, Barrier
        summary["s_matrix"] = _pairs(result.s_matrix)
        summary["z_matrix"] = _pairs(result.z_matrix)
        summary["touchstone"] = touchstone.name
        paths["touchstone"] = touchstone
        if "sar" in outputs:
            summary["sar"] = _sar_summary(study, outputs["sar"], comm)
        if "b1_plus" in outputs:
            summary["b1_plus"] = _b1_plus_summary(study, outputs["b1_plus"], comm)
        if "probes" in outputs:  # WF-1 step 20 (§9 item 564)
            summary["probes"] = _em_probe_summary(study, outputs["probes"], comm)
        if "tune" in outputs:  # `WF-1` step 11 (§9 item 442): reuses study.result
            summary["tune"] = _tune_summary(study, outputs["tune"], comm)
    elif isinstance(study.geometry, FHumanBirdcageGeometry) and isinstance(
        study.excitation, LumpedPortDrive
    ):  # `WF-1` step 12 (§9 item 484), E5: the 32-port ring route
        summary["cells"] = int(result.cells)
        n = int(result.s_matrix.shape[0])
        touchstone = out / f"{path.stem}.s{n}p"
        result.write_touchstone(touchstone)  # collective: rank 0 writes, Barrier
        summary["s_matrix"] = _pairs(result.s_matrix)
        summary["z_matrix"] = _pairs(result.z_matrix)
        summary["touchstone"] = touchstone.name
        paths["touchstone"] = touchstone
        if "tune" in outputs:
            summary["tune"] = _tune_summary(study, outputs["tune"], comm)
    elif isinstance(study.geometry, FHumanBirdcageGeometry):  # `WF-1` step 7a (§9 item 390a)
        drive = study.excitation

        def _ri(z):
            z = complex(z)
            return [float(z.real), float(z.imag)]

        summary["cells"] = int(result.cells)
        summary["excitation"] = {
            "kind": "two_port_quadrature",
            "c_in_f": float(drive.c_in_f),
            "esr_ohm": float(drive.esr_ohm),
            "feeds": [str(f) for f in drive.feeds],
            "z0_ohm": float(drive.z0_ohm),
            "gap_impedance_ohm": _ri(result.gap_impedance_ohm),
            "feed_impedance_ohm": _ri(result.feed_impedance_ohm),
            "source_voltage_v": _ri(result.source_voltage_v),
        }
        summary["az_feeds"] = [float(a) for a in result.az_feeds]
        summary["az_mirror"] = float(result.az_mirror)
        summary["s_node"] = _pairs(result.s_node)
        summary["y_in"] = {
            str(sense): [_ri(z) for z in vals] for sense, vals in result.y_in.items()
        }
        if "sar" in outputs:
            summary["sar"] = _sar_summary(study, outputs["sar"], comm)
        if "b1_plus" in outputs:
            summary["b1_plus"] = _b1_plus_summary(study, outputs["b1_plus"], comm)
        if cooling is not None:  # `WF-1` step 10a (§9 item 431a): an _FHumanThermal
            thermal = study.temperature(
                cooling.cooling, transient=cooling.transient, balls=cooling.balls,
                window_s=cooling.window_s, sense=cooling.sense, comm=comm,
            )  # reuses study.result (no second solve)
            summary["temperature"] = _tuned_temperature_summary(thermal)
    else:
        msh = result.mesh
        summary["cells"] = int(msh.topology.index_map(msh.topology.dim).size_global)
        sar = study.mean_sar()  # collective reduction inside
        summary["mean_sar"] = {k: float(v) for k, v in sar.items()}
        if isinstance(cooling, _BallThermal):  # `WF-1` step 29 (§9 item 873)
            import numpy as np

            tr = cooling.transient
            thermal = study.temperature(
                cooling.cooling, transient=tr, comm=comm
            )  # reuses study.result (no second solve)
            # _run_thermal_box's transient lines, copied.
            T = thermal.temperature
            n_own = (T.function_space.dofmap.index_map.size_local
                     * T.function_space.dofmap.index_map_bs)
            theta = np.real(np.asarray(T.x.array[:n_own])) - float(cooling.cooling.t_inf_c)
            local = float(np.max(theta)) if theta.size else -np.inf
            summary["temperature"] = {
                "transient": {
                    "t_end_s": float(tr.dt_s * tr.n_steps),
                    "dt_s": float(tr.dt_s),
                    "n_steps": int(tr.n_steps),
                    "theta": float(tr.theta),
                    "startup_half_steps": int(tr.startup_half_steps),
                    "rho_c_j_per_m3_k": float(tr.rho_c_j_per_m3_k),
                },
                "source_power_w": float(np.real(thermal.source_power)),
                "source_multipliers": [float(m) for m in thermal.source_multipliers],
                "max_balance_residual": float(thermal.max_balance_residual),
                "theta_max_k": float(comm.allreduce(local, op=MPI.MAX)),
            }
        elif cooling is not None:
            thermal = study.temperature(cooling)  # reuses study.result
            summary["temperature"] = _temperature_summary(thermal, comm)
            if "probes" in outputs:
                # Study.probe reuses study.result (no EM re-solve) and re-solves
                # the thermal problem; the values come from that θ (disclosed).
                probes = study.probe(outputs["probes"], cooling)  # collective
                summary["probes"] = {
                    "centres": [[float(x) for x in pb.centre] for pb in probes.probes],
                    "radii": [float(pb.radius) for pb in probes.probes],
                    "values_k": [float(v) for v in probes.values],
                    "volumes_m3": [float(v) for v in probes.volumes],
                }
    summary["elapsed_s"] = float(time.perf_counter() - started)
    summary_path = out / "summary.json"
    paths["summary"] = summary_path
    if comm.rank == 0:
        out.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w") as fh:
            json.dump(summary, fh, indent=1)
    comm.Barrier()
    return StudyRun(
        study=study, result=result, summary=summary, paths=paths, thermal=thermal, probes=probes
    )
