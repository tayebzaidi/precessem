"""`Study` v0 (`WF-15`): one scriptable object composing an existing solve.

A ``Study`` composes five things that already exist in the package:

* a geometry generator (``SphereInBoxGeometry`` wraps
  ``MeshGenerator.sphere_in_box_domain``);
* a ``MAT-1`` material (``GelledSalinePhantomMaterial``, or a preset name);
* an imposed field (``QuasiStaticSphereField``: the quasi-static dielectric
  sphere field on the box wall, ``utils.analytical``);
* a ``TimeHarmonicSolver`` solve, **degree 2 by default** (operator order
  rule 2026-10-01; the solver's own default is 1);
* a SAR result (``post.sar.mean_sar``).

No solver code lives here, and nothing is imported from ``tests/``.  v0 is one
fixture family (sphere in a box, imposed field); ports, coils and YAML are
out of scope (`WF-16` / `WF-1`).
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field as _dc_field
from typing import Any, Optional, Union

import numpy as np
from mpi4py import MPI

from .core import HomogeneousMaterial, TimeHarmonicProblem, TimeHarmonicSolver
from .io.mesh import MeshGenerator
from .materials import ColeColeMaterial, GelledSalinePhantomMaterial
from .post.sar import mean_sar
from .utils.analytical import dielectric_sphere_quasistatic_field
from .utils.constants import EPSILON_0

__all__ = [
    "SphereInBoxGeometry",
    "QuasiStaticSphereField",
    "Study",
    "StudyResult",
    "FourPortBirdcageGeometry",
    "LumpedPortDrive",
    "CapacitorSheets",
    "TwoPortQuadratureDrive",
    "PortStudySweepResult",
    "RobinCooling",
    "ThermalBoxGeometry",
    "ThermalBallGeometry",
    "UniformHeatSource",
    "TransverseQuadraticHeatSource",
    "TransientSpec",
    "DutyCycle",
    "ContainerWall",
]


@dataclass(frozen=True)
class SphereInBoxGeometry:
    """A sphere (tag ``sphere_tag``) centred in a cubic air box."""

    sphere_radius: float
    box_half_width: float
    resolution_sphere: float
    resolution_far: float
    sphere_tag: int = 1

    def build(self, comm: MPI.Intracomm = MPI.COMM_WORLD):
        """Return ``(mesh, cell_tags, facet_tags)`` from the generator."""
        return MeshGenerator.sphere_in_box_domain(
            sphere_radius=self.sphere_radius,
            box_half_width=self.box_half_width,
            resolution_sphere=self.resolution_sphere,
            resolution_far=self.resolution_far,
            comm=comm,
        )


@dataclass(frozen=True)
class QuasiStaticSphereField:
    """Imposed quasi-static field of a dielectric sphere in a uniform ``e0 ẑ``."""

    e0: float = 1.0

    def dirichlet_field(self, *, radius: float, epsilon_c: complex):
        return dielectric_sphere_quasistatic_field(
            radius=radius, epsilon_c=epsilon_c, e0=self.e0
        )


@dataclass
class StudyResult:
    """What ``Study.solve`` returns."""

    mesh: Any
    cell_tags: Any
    facet_tags: Any
    fields: Any
    solver: TimeHarmonicSolver
    degree: int
    comm: Any


class Study:
    """Geometry + material + imposed field + solve + SAR, composed."""

    # `WF-15` step 8: the original ``ColeColeMaterial`` when one was resolved
    # at the Study's frequency; ``None`` on every other path.
    dispersive_material = None

    def __init__(
        self,
        geometry: SphereInBoxGeometry,
        material: Union[GelledSalinePhantomMaterial, str],
        excitation: QuasiStaticSphereField,
        *,
        frequency_hz: Optional[float] = None,
        degree: int = 2,
        density_kg_m3: float = 1000.0,
        conductor=None,
    ):
        if isinstance(geometry, ThermalBoxGeometry):
            # `WF-15` step 5 (§9 item 521a): a thermal-only box Study.
            if material is not None:
                raise TypeError("a ThermalBoxGeometry Study takes material=None")
            if not isinstance(excitation, (UniformHeatSource, TransverseQuadraticHeatSource)):
                raise TypeError(
                    f"a ThermalBoxGeometry Study needs a UniformHeatSource, "
                    f"got {type(excitation).__name__}"
                )
            self.geometry = geometry
            self.material = None
            self.excitation = excitation
            self.frequency_hz = None
            self.degree = int(degree)
            self.density_kg_m3 = float(density_kg_m3)
            self.conductor = None
            self.result = None
            return
        if isinstance(geometry, ThermalBallGeometry):
            # `WF-15` step 3 (§9 item 330): a thermal-only Study -- no material,
            # no field, no EM solve; the heat source is imposed.
            if material is not None:
                raise TypeError("a ThermalBallGeometry Study takes material=None")
            if not isinstance(excitation, UniformHeatSource):
                raise TypeError(
                    f"a ThermalBallGeometry Study needs a UniformHeatSource, "
                    f"got {type(excitation).__name__}"
                )
            self.geometry = geometry
            self.material = None
            self.excitation = excitation
            self.frequency_hz = None
            self.degree = int(degree)
            self.density_kg_m3 = float(density_kg_m3)
            self.conductor = None
            self.result = None
            return
        if frequency_hz is None:
            raise TypeError("frequency_hz is required for an EM Study")
        if isinstance(material, ColeColeMaterial):
            # `WF-15` step 8 (§9 item 653): a Cole–Cole medium is resolved at
            # the Study's frequency to the real (σ_eff, εᵣ) pair the solver
            # consumes; only the imposed-field sphere takes one.
            if not isinstance(geometry, SphereInBoxGeometry):
                raise TypeError(
                    f"a ColeColeMaterial material needs a SphereInBoxGeometry, "
                    f"got {type(geometry).__name__}"
                )
            self.dispersive_material = material
            sigma_eff, eps_r = material.sigma_epsilon_at(frequency_hz)
            material = GelledSalinePhantomMaterial(
                sigma=sigma_eff, epsilon_r=eps_r, frequency_hz=frequency_hz, mu_r=material.mu_r
            )
        if isinstance(material, str):
            material = GelledSalinePhantomMaterial.from_preset(
                material, frequency_hz=frequency_hz
            )
        if not isinstance(material, GelledSalinePhantomMaterial):
            raise TypeError(
                f"material must be a GelledSalinePhantomMaterial or a preset name, "
                f"got {type(material).__name__}"
            )
        self.geometry = geometry
        self.material = material
        self.excitation = excitation
        self.frequency_hz = float(frequency_hz)
        self.degree = int(degree)
        self.density_kg_m3 = float(density_kg_m3)
        # `WF-16` step 7 (§9 item 324): ``None`` -- the default -- is the
        # unchanged path; a ``LeontovichWall`` needs a hole-route port geometry.
        self.conductor = conductor
        self.result: Optional[StudyResult] = None

    @property
    def omega(self) -> float:
        _refuse_thermal_only_attr(self, "frequency_hz", "omega")
        return 2.0 * np.pi * self.frequency_hz

    def epsilon_c(self) -> complex:
        """``ε_c = εᵣ − jσ/(ωε₀)``, the expression the solver assembles."""
        _refuse_thermal_only_attr(self, "material", "epsilon_c")
        return complex(
            self.material.epsilon_r, -self.material.sigma / (self.omega * EPSILON_0)
        )

    def dirichlet_field(self):
        _refuse_thermal_only_attr(self, "excitation", "dirichlet_field")
        return self.excitation.dirichlet_field(
            radius=self.geometry.sphere_radius, epsilon_c=self.epsilon_c()
        )

    def solve(self, comm: MPI.Intracomm = MPI.COMM_WORLD, *, backend="dolfinx",
              case_dir=None, palace_ranks=8) -> StudyResult:
        _refuse_f_clinical(self, "solve")  # `WF-13` step 12 (§9 item 663), before any mesh
        _refuse_rod_route(self, "solve")  # `WF-15` step 12 (§9 item 729), before any mesh
        if backend != "dolfinx":
            if backend != "palace":
                raise ValueError(f"backend {backend!r} not in ('dolfinx', 'palace')")
            from .backends.palace import solve_palace

            return solve_palace(self, comm, case_dir=case_dir, palace_ranks=palace_ranks)
        if isinstance(self.geometry, ThermalBoxGeometry):
            raise TypeError("a thermal-only Study has no field solve; use temperature()")
        if isinstance(self.geometry, ThermalBallGeometry):
            raise TypeError("a thermal-only Study has no field solve; use temperature()")
        if isinstance(self.geometry, (FourPortBirdcageGeometry, FHumanBirdcageGeometry)):
            return self._solve_ports(comm)
        msh, cell_tags, facet_tags = self.geometry.build(comm)
        problem = TimeHarmonicProblem(
            mesh=msh,
            frequency_hz=self.frequency_hz,
            material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0),
            cell_tags=cell_tags,
            phantom_material=self.material,
            phantom_tag=self.geometry.sphere_tag,
            boundary_condition="pec_zero_tangential_a",
            dirichlet_e_field=self.dirichlet_field(),
        )
        solver = TimeHarmonicSolver(problem, degree=self.degree)
        fields = solver.solve()
        self.result = StudyResult(
            mesh=msh,
            cell_tags=cell_tags,
            facet_tags=facet_tags,
            fields=fields,
            solver=solver,
            degree=self.degree,
            comm=comm,
        )
        return self.result

    def mean_sar(self) -> dict:
        """``post.sar.mean_sar`` over the sphere, its dict unchanged."""
        _refuse_f_human(self, "mean_sar")
        if isinstance(self.geometry, FourPortBirdcageGeometry):
            raise TypeError("mean_sar is not defined on a port Study (WF-16)")
        _refuse_loaded(self, "mean_sar")  # WF-16 step 14 (§9 item 661)
        if self.result is None:
            self.solve()
        res = self.result
        return mean_sar(
            res.fields.e_complex,
            sigma=res.fields.sigma_field,
            rho=self.density_kg_m3,
            cell_tags=res.cell_tags,
            subdomain_ids=self.geometry.sphere_tag,
            comm=res.comm,
        )


# ---------------------------------------------------------------------------
# `WF-16` step 1: the F-small four-port birdcage on exact port rectangles.
# Additive; nothing above changed except one dispatch branch in ``solve`` and
# one refusal in ``mean_sar``.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class FourPortBirdcageGeometry:
    """Item 192's F-small gapped 4-leg birdcage, sheets on exact CAD rectangles.

    No κ and no ragged route (operator ruling 2026-10-02): ``strip_fraction``
    must be a number; it is both the generator's CAD strip fraction and the
    midpoint filter's fraction, as in the direct call.
    """

    leg_count: int = 4
    ring_radius: float = 0.07
    leg_width: float = 0.012
    leg_spacing: float = 0.11
    coil_length: float = 0.14
    ring_minor_radius: float = 0.004
    phantom_radius: float = 0.03
    phantom_height: float = 0.08
    port_box_size: tuple = (0.010, 0.008, 0.010)
    leg_gap_length: float = 0.008
    air_padding: float = 0.03
    resolution: float = 0.015
    conductor_resolution: float = 0.4 * 0.004
    strip_fraction: float = 0.5
    sheet_iface: int = 210
    port_lower: int = 100
    port_upper: int = 200
    conductor_tag: int = 1
    phantom_tag: int = 3
    conductor_sigma_s_per_m: float = 8.0e2
    # `WF-16` step 7 (§9 item 324), additive: ``False`` forwards nothing, so the
    # generator call is byte-identical; ``True`` cuts the coil out of the air
    # (`TH-15` step 3a's hole route), for a ``LeontovichWall`` or PEC conductor.
    as_hole: bool = False
    # `WF-15` step 10 (§9 item 709), additive: ``None`` forwards nothing, so the
    # generator call is byte-identical; a mapping is the generator's
    # ``implant_rod`` (a rod hole with a PEC wall, tag 402).  A dict, so out of
    # the hash; ``==`` still compares it.
    implant_rod: Optional[dict] = _dc_field(default=None, hash=False)

    def __post_init__(self):
        if self.strip_fraction is None or not isinstance(
            self.strip_fraction, (float, int)
        ):
            raise ValueError(
                "strip_fraction must be a number: ports are exact CAD rectangles "
                "(no ragged sheet, operator ruling 2026-10-02)"
            )

    @classmethod
    def f_small(cls) -> "FourPortBirdcageGeometry":
        return cls()

    def build(self, comm: MPI.Intracomm = MPI.COMM_WORLD) -> dict:
        """Mesh, cell tags, narrowed sheet facet tags and one record per port."""
        import time

        from .ports._sheet_geometry import (
            _facet_group_area,
            _narrowed_radial,
            _port_frame,
            _projected_extents,
            _sheet_areas,
            _sheet_azimuth_deg,
            _sheet_facet_count,
        )

        offsets = np.zeros(self.leg_count)
        _hole_extra = {"as_hole": True} if self.as_hole else {}
        _implant_extra = {"implant_rod": dict(self.implant_rod)} if self.implant_rod else {}
        started = time.perf_counter()
        msh, cell_tags, _facet_tags, diag = MeshGenerator.birdcage_port_domain(
            leg_count=self.leg_count,
            ring_radius=self.ring_radius,
            leg_width=self.leg_width,
            leg_spacing=self.leg_spacing,
            coil_length=self.coil_length,
            ring_minor_radius=self.ring_minor_radius,
            phantom_radius=self.phantom_radius,
            phantom_height=self.phantom_height,
            port_box_size=self.port_box_size,
            leg_gap_length=self.leg_gap_length,
            leg_azimuth_offsets_rad=offsets,
            emit_port_sheets=True,
            air_padding=self.air_padding,
            resolution=self.resolution,
            conductor_resolution=self.conductor_resolution,
            c4_congruent_sheets=False,
            comm=comm,
            return_diagnostics=True,
            leg_sheet_strip_fraction=float(self.strip_fraction),
            **_hole_extra,
            **_implant_extra,
        )
        t_mesh = time.perf_counter() - started
        tdim = msh.topology.dim
        ncells = int(msh.topology.index_map(tdim).size_global)
        msh.topology.create_connectivity(tdim - 1, tdim)
        msh.topology.create_entity_permutations()

        ports_idx = list(range(1, self.leg_count + 1))
        tags_f, full_areas = _sheet_areas(
            msh,
            cell_tags,
            ports_idx,
            comm,
            self.sheet_iface,
            self.port_lower,
            self.port_upper,
        )
        dx, _dy, dz = diag["port_box_size_m"]
        sheet_analytic = float(dx) * float(dz)

        geometry = {}
        for i in ports_idx:
            tag = self.sheet_iface + i
            azimuth = _sheet_azimuth_deg(msh, tags_f, tag, comm)
            radial, azimuthal, axial = _port_frame(azimuth)
            spans, centres = _projected_extents(
                msh, tags_f, tag, comm, [radial, azimuthal, axial]
            )
            geometry[i] = {
                "tag": tag,
                "azimuth_deg": float(azimuth),
                "radial": radial,
                "w_full": float(spans[0]),
                "out_of_plane_full": float(spans[1]),
                "h_full": float(spans[2]),
                "centre_radial": float(centres[0]),
                "area_full": float(full_areas[i]),
                "area_ratio_full": float(full_areas[i] / sheet_analytic),
            }

        for i in ports_idx:
            g = geometry[i]
            tags_f = _narrowed_radial(
                msh,
                tags_f,
                g["tag"],
                float(self.strip_fraction),
                g["centre_radial"],
                g["radial"],
                0.5 * g["w_full"],
            )

        sheets = []
        for i in ports_idx:
            g = geometry[i]
            n_facets = _sheet_facet_count(msh, tags_f, g["tag"], comm)
            if n_facets <= 0:
                raise RuntimeError(f"sheet {g['tag']}: no owned facets anywhere")
            area = _facet_group_area(msh, tags_f, g["tag"], comm)
            radial, azimuthal, axial = _port_frame(g["azimuth_deg"])
            spans, _centres = _projected_extents(
                msh, tags_f, g["tag"], comm, [radial, azimuthal, axial]
            )
            sheets.append(
                {
                    **g,
                    "facets": int(n_facets),
                    "area": float(area),
                    "h": float(spans[2]),
                    "w": float(area / spans[2]),
                    "w_bbox": float(spans[0]),
                    "out_of_plane": float(spans[1]),
                }
            )
        return {
            "mesh": msh,
            "cell_tags": cell_tags,
            "facet_tags": tags_f,
            "sheets": sheets,
            "cells": ncells,
            "sheet_analytic": sheet_analytic,
            "mesh_time": float(t_mesh),
            "mesh_facet_tags": _facet_tags,
        }


# ---------------------------------------------------------------------------
# `WF-16` step 8 (§9 item 355, chain U27): the F-human 32-port ring-gap
# birdcage (`WF-7` step 1's fixture) through the port Study.  Additive; the
# sheets are the full CAD ring-gap interfaces (no narrowing, κ never applied).
# ---------------------------------------------------------------------------

# `mesh:9`'s base rung (the `GEO-25` ladder's base), scaled by
# ``s = ring_radius / F_HUMAN_BASE_RING_RADIUS`` exactly as `GEO-25`'s
# ``_params`` scales it.  Conductor cross-section and sizing are not scaled.
F_HUMAN_BASE_RING_RADIUS = 0.07
F_HUMAN_BASE_LEG_SPACING = 0.11
F_HUMAN_BASE_COIL_LENGTH = 0.14
F_HUMAN_BASE_PHANTOM_RADIUS = 0.03
F_HUMAN_BASE_PHANTOM_HEIGHT = 0.08
F_HUMAN_BASE_PORT_BOX_SIZE = (0.010, 0.008, 0.010)
F_HUMAN_BASE_RING_GAP_LENGTH = 0.008
F_HUMAN_BASE_AIR_PADDING = 0.03


@dataclass(frozen=True)
class FHumanBirdcageGeometry:
    """`WF-7` step 1's F-human birdcage: 16 legs, 32 longitudinal ring-gap ports.

    ``emit_leg_ports=False`` (default) is `GEO-35` step 2a's leg-box-free
    fixture; ``True`` adds no key, so the generator dict is the boxed
    pre-switch build byte for byte.
    """

    leg_count: int = 16
    ring_radius: float = 0.15
    leg_width: float = 0.012
    ring_minor_radius: float = 0.004
    resolution: float = 0.015
    conductor_resolution: float = 0.4 * 0.004
    sheet_iface: int = 210
    port_lower: int = 100
    port_upper: int = 200
    conductor_tag: int = 1
    phantom_tag: int = 3
    conductor_sigma_s_per_m: float = 8.0e2
    emit_leg_ports: bool = False
    #: `WF-16` step 9 (§9 item 382): ``True`` is the copper hole route (the coil
    #: a cavity, leg port boxes kept); ``False`` -- the default -- is unchanged.
    as_hole: bool = False

    def generator_kwargs(self) -> dict:
        """The ``birdcage_port_domain`` kwargs (pure, no mesh)."""
        s = self.ring_radius / F_HUMAN_BASE_RING_RADIUS
        kwargs = dict(
            leg_count=self.leg_count,
            ring_radius=F_HUMAN_BASE_RING_RADIUS * s,
            leg_width=self.leg_width,
            leg_spacing=F_HUMAN_BASE_LEG_SPACING * s,
            coil_length=F_HUMAN_BASE_COIL_LENGTH * s,
            ring_minor_radius=self.ring_minor_radius,
            phantom_radius=F_HUMAN_BASE_PHANTOM_RADIUS * s,
            phantom_height=F_HUMAN_BASE_PHANTOM_HEIGHT * s,
            port_box_size=tuple(v * s for v in F_HUMAN_BASE_PORT_BOX_SIZE),
            leg_gap_length=None,
            ring_gap_length=F_HUMAN_BASE_RING_GAP_LENGTH * s,
            emit_port_sheets=True,
            air_padding=F_HUMAN_BASE_AIR_PADDING * s,
            resolution=self.resolution,
            conductor_resolution=self.conductor_resolution,
        )
        kwargs["ring_sheet_orientation"] = "longitudinal"
        if self.as_hole:
            # The hole route keeps the leg port boxes (`GEO-35` step 1 (K)).
            kwargs["as_hole"] = True
        elif not self.emit_leg_ports:
            kwargs["emit_leg_ports"] = False
        return kwargs

    def ring_ordinals(self) -> list:
        return list(range(self.leg_count + 1, 3 * self.leg_count + 1))

    def drive_direction(self, ordinal: int) -> tuple:
        """φ̂ at the ring gap's centre azimuth, from the ordinal alone."""
        first = self.leg_count + 1
        j = (ordinal - first) % self.leg_count
        phi_c = 2.0 * np.pi * j / self.leg_count + np.pi / self.leg_count
        return (float(-np.sin(phi_c)), float(np.cos(phi_c)), 0.0)

    def build(self, comm: MPI.Intracomm = MPI.COMM_WORLD) -> dict:
        """Mesh, cell tags, full ring-sheet facet tags and one record per port."""
        import time

        from .ports._sheet_geometry import _sheet_areas, _sheet_facet_count

        started = time.perf_counter()
        msh, cell_tags, _facet_tags, diag = MeshGenerator.birdcage_port_domain(
            comm=comm, return_diagnostics=True, **self.generator_kwargs()
        )
        t_mesh = time.perf_counter() - started
        tdim = msh.topology.dim
        ncells = int(msh.topology.index_map(tdim).size_global)
        ords = self.ring_ordinals()
        tags_f, areas = _sheet_areas(
            msh, cell_tags, ords, comm, self.sheet_iface, self.port_lower,
            self.port_upper,
        )
        msh.topology.create_connectivity(tdim - 1, tdim)
        msh.topology.create_entity_permutations()
        chord = float(diag["ring_port_layout"]["ring_port_gap_chord_m"])
        sheets = []
        for i in ords:
            tag = self.sheet_iface + i
            sheets.append(
                {
                    "tag": tag,
                    "ordinal": i,
                    "facets": _sheet_facet_count(msh, tags_f, tag, comm),
                    "area": float(areas[i]),
                    "h": chord,
                    "w": float(areas[i]) / chord,
                    "drive_direction": self.drive_direction(i),
                    "orientation": f"ring_gap_{i}_phi_hat",
                }
            )
        return {
            "mesh": msh,
            "cell_tags": cell_tags,
            "facet_tags": tags_f,
            "sheets": sheets,
            "cells": ncells,
            "chord": chord,
            "mesh_time": float(t_mesh),
            "mesh_facet_tags": _facet_tags,
        }


def _refuse_f_human(study, method: str) -> None:
    """`WF-16` step 8: the F-human Study is the S route only (before any mesh);
    estimate() is priced, `WF-13` steps 10-11."""
    if isinstance(study.geometry, FHumanBirdcageGeometry):
        raise TypeError(f"{method} is not defined on an FHumanBirdcageGeometry Study "
                        "(WF-16 step 8: S-parameters only; estimate() is priced, "
                        "WF-13 steps 10-11)")


# ---------------------------------------------------------------------------
# `WF-13` step 12 (§9 item 663, chain U96): the F-clinical-1.5T birdcage
# (`GEO-38`'s r1 fixture) through ``Study.estimate`` only.  Additive: the kind
# carries no builder; ``estimate`` reads `GEO-38`'s record constants and the
# measured ``PRICE_TABLE_CLINICAL`` row; ``solve`` is refused before any mesh.
# ---------------------------------------------------------------------------

#: Copied from ``tests/validation/test_geo38_f_clinical_build.py:86``
#: (``F_CLINICAL_R1_CELL_RECORD``) and ``:89`` (``F_CLINICAL_R1_UNKNOWNS_RECORD``,
#: ``(U1 = E, U2 = 2E + 2F)``); ``src`` may not import from ``tests/``, so the
#: new test asserts equality to the test constants.
F_CLINICAL_R1_CELL_RECORD = 1_147_366
F_CLINICAL_R1_UNKNOWNS_RECORD = (1_525_193, 7_885_674)
#: `WF-13` step 15 (§9 item 782, §10 U118): the F-clinical-3T r1 twin, copied
#: from ``tests/validation/test_geo38_f_clinical_3t.py:105``
#: (``F_CLINICAL_3T_R1_CELL_RECORD``) and ``:108``
#: (``F_CLINICAL_3T_R1_UNKNOWNS_RECORD``); the new test asserts equality.
F_CLINICAL_3T_R1_CELL_RECORD = 1_057_895
F_CLINICAL_3T_R1_UNKNOWNS_RECORD = (1_425_765, 7_335_988)
#: `WF-13` step 16 (§9 item 783, §10 U119): entry 48's rod in the ASTM box
#: (`GEO-36` I36): ``tests/mesh/test_birdcage_implant_rod.py:97``'s ``ROD``
#: moved to ``_rod_centre()`` (``tests/validation/test_geo36_rod_f_clinical.py:88``:
#: x = 0.25 Lx, Lx = 2 B_BOX = 0.42 m; y = the table offset, 0).  Records copied
#: from ``test_geo36_rod_f_clinical.py:410`` / ``:412``; the new test asserts equality.
F_CLINICAL_ROD_ENTRY48 = {"radius": 1.0e-3, "length": 0.030,
                          "centre": (0.25 * (2.0 * 0.21), 0.0, 0.0), "axis": (0.0, 0.0, 1.0)}
F_CLINICAL_ROD_CELL_RECORD = 1_217_025
F_CLINICAL_ROD_UNKNOWNS_RECORD = (1_606_744, 8_327_394)


@dataclass(frozen=True)
class FClinicalBirdcageGeometry:
    """`GEO-38`'s F-clinical-1.5T birdcage (32 legs, 64 ring-gap sheets, r1 sizing).

    The 1.5 T fixture only: the fields name ``F_CLINICAL_15T_KWARGS``
    (``tests/validation/test_geo38_f_clinical_build.py:44``) at the r1
    conductor resolution and are informational -- there is no ``build``; a
    solve is xl tier only.  ``Study.estimate`` reads record constants.  A saved
    64-port matrix loads and tunes (`WF-16` step 16); a solve stays xl tier.
    """

    leg_count: int = 32
    ring_radius: float = 0.355
    shield_radius: float = 0.3725
    resolution: float = 0.015
    conductor_resolution: float = 1.6e-3
    field_strength_t: float = 1.5
    #: `WF-13` step 16 (§9 item 783, §10 U119): estimate-only.  ``None`` or
    #: exactly entry 48's rod (``F_CLINICAL_ROD_ENTRY48``) at 1.5 T; any other
    #: value is refused.  A dict, so out of the hash (as on the four-port).
    implant_rod: Optional[dict] = _dc_field(default=None, hash=False)

    def __post_init__(self):
        rod = self.implant_rod
        if rod is None:
            return
        if self.field_strength_t != 1.5 or rod != F_CLINICAL_ROD_ENTRY48:
            raise ValueError(
                "FClinicalBirdcageGeometry.implant_rod accepts only entry 48's rod "
                f"{F_CLINICAL_ROD_ENTRY48!r} at field_strength_t=1.5 (WF-13 step 16), "
                f"got {rod!r} at {self.field_strength_t}")


def _f_clinical_tier(degree: int) -> str:
    """The tier an F-clinical solve runs on, for the refusals (`WF-13` step 13).

    Degree 2 is the 10-31 weekly's: a root grid of ≤ 3 nodes fits ``xl``, a
    column set needs ``xxl``.  Every text carries "xl tier only".
    """
    if int(degree) == 2:
        return ("xl tier only at <= 3 root-grid nodes at degree 2, "
                "and xxl tier for a column set")
    return "xl tier only"


def _refuse_f_clinical(study, method: str) -> None:
    """`WF-13` step 12: an F-clinical Study is priced only (before any mesh)."""
    if isinstance(study.geometry, FClinicalBirdcageGeometry):
        raise ValueError(f"{method} on an f_clinical_birdcage (FClinicalBirdcageGeometry) "
                         f"Study is {_f_clinical_tier(study.degree)} "
                         "(WF-13 step 12: estimate() is priced)")


def _refuse_rod_route(study, method: str) -> None:
    """`WF-15` step 12 (§9 item 729): a rod solves at degree 1 with lumped ports
    only (before any mesh).  Neither other route is known wrong; no step gates them."""
    geometry = study.geometry
    if isinstance(geometry, FClinicalBirdcageGeometry) and geometry.implant_rod:
        # `WF-13` step 16 (§9 item 783): the clinical rod is estimate-only, by tier.
        _refuse_f_clinical(study, method)
    if not isinstance(geometry, FourPortBirdcageGeometry) or not geometry.implant_rod:
        return
    # `WF-15` step 13 (§9 item 802): lifted on exactly the path that step gates --
    # ``solve`` on the PEC four-port (no ``as_hole``, no conductor) at degree 1.  Every
    # other entry point and the copper hole route keep this TypeError.
    if isinstance(study.excitation, CapacitorSheets) and (
            method != "solve" or geometry.as_hole
            or getattr(study, "conductor", None) is not None):
        raise TypeError(f"{method}: implant_rod with CapacitorSheets is not gated (WF-15 step 12); "
                        "only Study.solve on the PEC four-port at degree 1 is (WF-15 step 13)")
    # The degree refusal below still applies to that path at degree != 1.
    if study.degree != 1:
        raise ValueError(f"{method}: implant_rod at degree {study.degree} is not gated; "
                         "it needs degree 1 (WF-15 step 12)")


def _refuse_thermal_only(study, method: str) -> None:
    """`WF-13` step 9: a thermal-only Study has no EM solve to price (before any mesh)."""
    if isinstance(study.geometry, ThermalBoxGeometry):
        raise TypeError(f"{method} is not defined on a thermal-only Study "
                        "(ThermalBoxGeometry has no EM solve; use temperature())")
    if isinstance(study.geometry, ThermalBallGeometry):
        raise TypeError(f"{method} is not defined on a thermal-only Study "
                        "(ThermalBallGeometry has no EM solve; use temperature())")


def _refuse_thermal_only_attr(study, field: str, method: str) -> None:
    """`WF-15` step 6 (§9 item 523): an EM helper on a thermal-only Study
    raises ``StudyConfigError`` naming ``field`` (before any mesh)."""
    if study.frequency_hz is None:
        # Function-local: ``study_yaml`` imports this module at load time.
        from .study_yaml import StudyConfigError

        raise StudyConfigError(
            "<Study>", field,
            f"{method} is not defined on a thermal-only Study "
            "(no EM material / frequency / imposed field)",
        )


@dataclass(frozen=True)
class LumpedPortDrive:
    """All-port lumped-sheet drive/sweep: one ``LumpedSheetPortSpec`` per sheet."""

    z0_ohm: float = 50.0
    port_impedance_ohm: float = 50.0
    drive_voltage_v: complex = 1.0 + 0.0j

    def specs_and_definitions(self, geometry: FourPortBirdcageGeometry, sheets):
        from .ports.definitions import PortDefinition
        from .ports.lumped import LumpedSheetPortSpec

        specs = []
        for s in sheets:
            specs.append(
                LumpedSheetPortSpec(
                    port_id=f"P{s['tag'] - geometry.sheet_iface}",
                    facet_tag=int(s["tag"]),
                    port_impedance_ohm=self.port_impedance_ohm,
                    gap_height_m=s["h"],
                    sheet_width_m=s["w"],
                    drive_direction=s.get("drive_direction", (0.0, 0.0, 1.0)),
                    drive_voltage_v=self.drive_voltage_v,
                    interior=True,
                )
            )
        port_defs = []
        for s in sheets:
            port_defs.append(
                PortDefinition(
                    port_id=f"P{s['tag'] - geometry.sheet_iface}",
                    positive_tag=int(s["tag"]),
                    negative_tag=geometry.conductor_tag,
                    orientation=s.get("orientation", "leg_gap_axial_plus_z"),
                    z0_ohm=self.z0_ohm,
                )
            )
        return specs, port_defs


@dataclass
class PortStudyResult(StudyResult):
    """``StudyResult`` plus the S-parameter sweep of a port Study (`WF-16`)."""

    problem: Any = None
    sweep: Any = None
    s_matrix: Any = None
    z_matrix: Any = None
    port_definitions: Any = None
    specs: Any = None
    sheets: Any = None
    cells: int = 0
    mesh_time: float = 0.0
    sweep_time: float = 0.0
    conductor_terms: Any = None
    conductor_facet_tags: Any = None

    def write_touchstone(self, path):
        """Rank 0 writes (no CSV companion), ``Barrier``; the path on every rank."""
        from pathlib import Path

        from .ports.touchstone import export_touchstone

        path = Path(path)
        if self.comm.rank == 0:
            path.parent.mkdir(parents=True, exist_ok=True)
            export_touchstone(
                [self.sweep],
                self.port_definitions,
                output_path=path,
                write_csv_companion=False,
            )
        self.comm.Barrier()
        return path


def _study_solve_ports(self, comm: MPI.Intracomm = MPI.COMM_WORLD) -> PortStudyResult:
    import time

    from .ports.sparameters import run_n_port_sparameter_sweep

    if not isinstance(self.excitation, LumpedPortDrive):
        hook = getattr(self.excitation, "_solve_port_study", None)
        if hook is not None:
            return hook(self, comm)
        raise TypeError("a FourPortBirdcageGeometry Study needs a LumpedPortDrive")
    if (isinstance(self.geometry, FHumanBirdcageGeometry) and self.conductor is not None
            and not self.geometry.as_hole):
        raise TypeError("the solid FHumanBirdcageGeometry is the sigma = 800 route only; "
                        "a conductor needs as_hole=True (WF-16 step 8 / 12)")
    geom = self.geometry
    built = geom.build(comm)
    msh = built["mesh"]
    cond_kwargs, cond_terms = _conductor_setup(self, built, comm)
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=self.frequency_hz,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=built["cell_tags"],
        # `WF-16` step 7: the hole route has no conductor cells to map.
        material_map=({} if geom.as_hole else {
            geom.conductor_tag: HomogeneousMaterial(
                sigma=geom.conductor_sigma_s_per_m, epsilon_r=1.0, mu_r=1.0
            ),
        }),
        phantom_material=self.material,
        phantom_tag=geom.phantom_tag,
        boundary_condition="pec_zero_tangential_a",
        **cond_kwargs,
    )
    specs, port_defs = self.excitation.specs_and_definitions(geom, built["sheets"])
    comm.Barrier()
    t0 = time.perf_counter()
    sweep = run_n_port_sparameter_sweep(
        problem,
        port_defs,
        lumped_sheet_ports=specs,
        lumped_sheet_facet_tags=built["facet_tags"],
        degree=self.degree,
        **({"extra_bilinear_terms": cond_terms} if cond_terms else {}),
    )
    comm.Barrier()
    t_sweep = time.perf_counter() - t0
    self.result = PortStudyResult(
        mesh=msh,
        cell_tags=built["cell_tags"],
        facet_tags=built["facet_tags"],
        fields=None,
        solver=None,
        degree=self.degree,
        comm=comm,
        problem=problem,
        sweep=sweep,
        s_matrix=np.asarray(sweep.s_matrix, dtype=np.complex128),
        z_matrix=np.asarray(sweep.z_matrix, dtype=np.complex128),
        port_definitions=port_defs,
        specs=specs,
        sheets=built["sheets"],
        cells=built["cells"],
        mesh_time=built["mesh_time"],
        sweep_time=float(t_sweep),
        conductor_terms=cond_terms,
        conductor_facet_tags=built["mesh_facet_tags"] if cond_terms else None,
    )
    return self.result


Study._solve_ports = _study_solve_ports


# ---------------------------------------------------------------------------
# `WF-16` step 10 (§9 item 373, chain U37): capacitor sheets with ESR.
# Insertions only, plus the three-line hook inside `_study_solve_ports`'s
# refusal branch (a `LumpedPortDrive` never enters it).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CapacitorSheets:
    """Series-RC terminations on the named port sheets; the rest are lumped ports.

    Every sheet whose port id is in ``on`` carries
    ``series_rlc_impedance(f, r_ohm=esr_ohm, l_h=0, c_f=c_f)``; every other
    sheet is a lumped port with ``drive``'s reference impedance, impedance and
    voltage, driven in turn.  ``Study.solve`` returns the S-matrix of the
    reduced network over the kept ports (ascending port order).  A record of
    the circuit reduction, not a tuning result: ``c_f`` is an input.
    """

    c_f: float
    esr_ohm: float = 0.0
    on: tuple = ()
    drive: LumpedPortDrive = _dc_field(default_factory=LumpedPortDrive)

    def __post_init__(self):
        c = float(self.c_f)
        if not np.isfinite(c) or c <= 0.0:
            raise ValueError("c_f must be finite and positive")
        r = float(self.esr_ohm)
        if not np.isfinite(r) or r < 0.0:
            raise ValueError("esr_ohm must be finite and non-negative")
        if isinstance(self.on, str):
            raise ValueError("on is a tuple of port ids, e.g. ('P2',)")
        on = tuple(self.on)
        if not on:
            raise ValueError("on must name at least one port sheet")
        if len(set(on)) != len(on):
            raise ValueError("on names a port sheet twice")
        object.__setattr__(self, "on", on)
        if not isinstance(self.drive, LumpedPortDrive):
            raise ValueError("drive must be a LumpedPortDrive")

    def impedance_ohm(self, frequency_hz: float) -> complex:
        from .ports.lumped import series_rlc_impedance

        return series_rlc_impedance(frequency_hz, r_ohm=self.esr_ohm, l_h=0.0, c_f=self.c_f)

    def _solve_port_study(self, study, comm):
        return _study_solve_capacitor_sheets(study, comm)


@dataclass
class CapacitorStudyResult(PortStudyResult):
    """``PortStudyResult`` over the kept ports, plus the capacitor sheets."""

    all_port_definitions: Any = None
    capacitor_port_ids: Any = None
    capacitor_impedance_ohm: Any = None
    solve_times: Any = None

    def write_touchstone(self, path):
        raise NotImplementedError("a reduced-network Touchstone is out of scope (WF-16 step 10)")


def _study_solve_capacitor_sheets(self, comm: MPI.Intracomm = MPI.COMM_WORLD):
    import dataclasses
    import time

    from .ports.lumped import run_lumped_sheet_port_case
    from .ports.sparameters import _power_waves

    exc = self.excitation
    if isinstance(self.geometry, FHumanBirdcageGeometry):
        raise TypeError("CapacitorSheets is defined on the F-small four-port Study only "
                        "(WF-16 step 10)")
    if not isinstance(self.geometry, FourPortBirdcageGeometry):
        raise TypeError("CapacitorSheets needs a FourPortBirdcageGeometry")
    geom = self.geometry
    valid = tuple(f"P{i + 1}" for i in range(int(geom.leg_count)))
    bad = [p for p in exc.on if p not in valid]
    if bad:
        raise ValueError(f"on names {bad}, outside {valid}")
    if len(exc.on) >= len(valid):
        raise ValueError("on leaves no port to drive")
    # The problem lines of `_study_solve_ports`, copied (not factored out).
    built = geom.build(comm)
    msh = built["mesh"]
    cond_kwargs, cond_terms = _conductor_setup(self, built, comm)
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=self.frequency_hz,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=built["cell_tags"],
        material_map=({} if geom.as_hole else {
            geom.conductor_tag: HomogeneousMaterial(
                sigma=geom.conductor_sigma_s_per_m, epsilon_r=1.0, mu_r=1.0
            ),
        }),
        phantom_material=self.material,
        phantom_tag=geom.phantom_tag,
        boundary_condition="pec_zero_tangential_a",
        **cond_kwargs,
    )
    specs, port_defs = exc.drive.specs_and_definitions(geom, built["sheets"])
    z_cap = exc.impedance_ohm(self.frequency_hz)
    specs_t = [dataclasses.replace(s, port_impedance_ohm=z_cap) if s.port_id in exc.on else s
               for s in specs]
    kept = sorted((p for p in port_defs if p.port_id not in exc.on),
                  key=lambda p: int(p.port_id[1:]))
    z0s = {float(p.z0_ohm) for p in kept}
    if len(z0s) != 1:
        raise ValueError("the kept ports need one common z0")
    z0 = z0s.pop()
    from .core import TimeHarmonicSolver

    shared = TimeHarmonicSolver(problem, degree=self.degree)
    k = len(kept)
    s_red = np.zeros((k, k), dtype=np.complex128)
    times = []
    comm.Barrier()
    t0 = time.perf_counter()
    for col, drv in enumerate(kept):
        t1 = time.perf_counter()
        res = run_lumped_sheet_port_case(
            problem, port_defs, specs_t, facet_tags=built["facet_tags"],
            driven_port_id=drv.port_id, verbose=False, degree=self.degree,
            solver=shared, extra_bilinear_terms=(cond_terms if cond_terms else None),
        )
        d = res.responses[drv.port_id]
        a_d, _ = _power_waves(d.voltage_v, d.current_a, z0)
        if a_d == 0:
            raise RuntimeError("zero incident wave on the driven port")
        for row, rcv in enumerate(kept):
            rr = res.responses[rcv.port_id]
            _, b = _power_waves(rr.voltage_v, rr.current_a, z0)
            s_red[row, col] = b / a_d
        times.append(time.perf_counter() - t1)
    comm.Barrier()
    t_sweep = time.perf_counter() - t0
    self.result = CapacitorStudyResult(
        mesh=msh,
        cell_tags=built["cell_tags"],
        facet_tags=built["facet_tags"],
        fields=None,
        solver=None,
        degree=self.degree,
        comm=comm,
        problem=problem,
        sweep=None,
        s_matrix=s_red,
        z_matrix=None,
        port_definitions=kept,
        specs=specs_t,
        sheets=built["sheets"],
        cells=built["cells"],
        mesh_time=built["mesh_time"],
        sweep_time=float(t_sweep),
        conductor_terms=cond_terms,
        conductor_facet_tags=built["mesh_facet_tags"] if cond_terms else None,
        all_port_definitions=port_defs,
        capacitor_port_ids=tuple(exc.on),
        capacitor_impedance_ohm=z_cap,
        solve_times=tuple(times),
    )
    return self.result


# ---------------------------------------------------------------------------
# `WF-16` step 2: a discrete frequency sweep through the port Study.
# Insertions only.  One mesh, one factorisation per frequency.  The gel is
# constant-property across the sweep: σ, εᵣ and μᵣ are copied and only the
# frequency tag follows the problem (there is no dispersion model here).
# ---------------------------------------------------------------------------


@dataclass
class PortStudySweepResult:
    """What ``Study.solve_sweep`` returns: one ``PortStudyResult`` per frequency."""

    frequencies_hz: tuple
    results: list
    mesh_time: float
    sweep_times: tuple
    comm: Any

    @property
    def s_matrices(self) -> np.ndarray:
        return np.asarray([r.s_matrix for r in self.results], dtype=np.complex128)

    def write_touchstone(self, path):
        """One multi-frequency file; rank 0 writes (no CSV companion), ``Barrier``."""
        from pathlib import Path

        from .ports.touchstone import export_touchstone

        path = Path(path)
        if self.comm.rank == 0:
            path.parent.mkdir(parents=True, exist_ok=True)
            export_touchstone(
                [r.sweep for r in self.results],
                self.results[0].port_definitions,
                output_path=path,
                write_csv_companion=False,
            )
        self.comm.Barrier()
        return path


# `WF-1` step 23 (§9 item 589): the YAML refusal's reason
# (``study_yaml.py``, field ``study.sweep_hz``), byte for byte.
_SWEEP_CONDUCTOR_REASON = (
    "a sweep with a conductor is not supported (Study.solve_sweep has no conductor route)"
)


def _validate_sweep_frequencies(frequencies_hz) -> tuple:
    freqs = tuple(float(f) for f in frequencies_hz)
    if len(freqs) == 0:
        raise ValueError("solve_sweep needs at least one frequency")
    for f in freqs:
        if not np.isfinite(f) or f <= 0.0:
            raise ValueError(f"sweep frequency must be finite and positive, got {f!r}")
    for a, b in zip(freqs, freqs[1:]):
        if not b > a:
            raise ValueError(
                "sweep frequencies must be strictly increasing (result order is input "
                f"order; Touchstone needs sorted points): {a!r} then {b!r}"
            )
    return freqs


def _study_solve_sweep(
    self, frequencies_hz, comm: MPI.Intracomm = MPI.COMM_WORLD, *, fields_at=None
) -> PortStudySweepResult:
    """Solve the port Study at each frequency on one mesh (`WF-16` step 2).

    ``self.frequency_hz``, ``self.material`` and ``self.result`` are not
    changed; the sweep is also stored as ``self.sweep_result``.  The phantom is
    constant-property: ``dataclasses.replace(self.material, frequency_hz=f)``.
    A Study with a ``conductor`` is refused (``TypeError`` naming ``conductor``,
    before the frequency list is checked or anything is built; `WF-1` step 23).

    ``fields_at`` (`WF-16` step 18, §9 item 780): one frequency of
    ``frequencies_hz`` (``float`` equality) whose ``results[i].fields`` keeps the
    drive fields (``{port_id: TimeHarmonicFields}``, taken inside that
    frequency's own sweep call, ``keep_fields=True``); every other result keeps
    ``fields=None``.  A value not in ``frequencies_hz`` raises ``ValueError``
    naming ``fields_at``, before anything is built.  ``None`` (the default) is
    the old contract: every ``fields`` is ``None``.
    """
    import dataclasses
    import time

    from .ports.sparameters import run_n_port_sparameter_sweep

    _refuse_rod_route(self, "solve_sweep")  # `WF-15` step 12 (§9 item 729)
    if not isinstance(self.geometry, FourPortBirdcageGeometry):
        raise TypeError("solve_sweep is defined on a port Study only (WF-16)")
    if not isinstance(self.excitation, LumpedPortDrive):
        raise TypeError("a FourPortBirdcageGeometry Study needs a LumpedPortDrive")
    if self.conductor is not None:  # WF-1 step 23 (§9 item 589)
        raise TypeError(f"conductor: {_SWEEP_CONDUCTOR_REASON}")
    freqs = _validate_sweep_frequencies(frequencies_hz)
    keep_at = None
    if fields_at is not None:  # WF-16 step 18 (§9 item 780), before any mesh
        try:
            keep_at = float(fields_at)
        except (TypeError, ValueError):
            raise ValueError(f"fields_at must be one of frequencies_hz, got {fields_at!r}")
        if not any(float(f) == keep_at for f in freqs):
            raise ValueError(
                f"fields_at must be one of frequencies_hz {list(freqs)!r}, got {fields_at!r}"
            )
    geom = self.geometry
    built = geom.build(comm)
    msh = built["mesh"]
    specs, port_defs = self.excitation.specs_and_definitions(geom, built["sheets"])
    results = []
    times = []
    for f in freqs:
        problem = TimeHarmonicProblem(
            mesh=msh,
            frequency_hz=f,
            material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
            cell_tags=built["cell_tags"],
            material_map={
                geom.conductor_tag: HomogeneousMaterial(
                    sigma=geom.conductor_sigma_s_per_m, epsilon_r=1.0, mu_r=1.0
                ),
            },
            phantom_material=dataclasses.replace(self.material, frequency_hz=f),
            phantom_tag=geom.phantom_tag,
            boundary_condition="pec_zero_tangential_a",
        )
        keep = keep_at is not None and float(f) == keep_at
        comm.Barrier()
        t0 = time.perf_counter()
        sweep = run_n_port_sparameter_sweep(
            problem,
            port_defs,
            lumped_sheet_ports=specs,
            lumped_sheet_facet_tags=built["facet_tags"],
            degree=self.degree,
            **({"keep_fields": True} if keep else {}),
        )
        comm.Barrier()
        t_sweep = time.perf_counter() - t0
        times.append(float(t_sweep))
        results.append(
            PortStudyResult(
                mesh=msh,
                cell_tags=built["cell_tags"],
                facet_tags=built["facet_tags"],
                fields=dict(sweep.fields) if keep else None,
                solver=None,
                degree=self.degree,
                comm=comm,
                problem=problem,
                sweep=sweep,
                s_matrix=np.asarray(sweep.s_matrix, dtype=np.complex128),
                z_matrix=np.asarray(sweep.z_matrix, dtype=np.complex128),
                port_definitions=port_defs,
                specs=specs,
                sheets=built["sheets"],
                cells=built["cells"],
                mesh_time=built["mesh_time"],
                sweep_time=float(t_sweep),
            )
        )
    self.sweep_result = PortStudySweepResult(
        frequencies_hz=freqs,
        results=results,
        mesh_time=built["mesh_time"],
        sweep_times=tuple(times),
        comm=comm,
    )
    return self.sweep_result


Study.solve_sweep = _study_solve_sweep


# ---------------------------------------------------------------------------
# `WF-16` step 6(a): the mode capacitance of a port Study from its own S.
# Additive; no existing line changed.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class TuneResult:
    """Lossless modal tuning read of a four-port Study (``Study.tune``), or the
    32-port F-human ring read (``z_k`` = the ring-antisymmetric eigenvalue per
    k = 0…N/2)."""

    k: int
    frequency_hz: float
    z_k: Any
    c_mode_farad: float
    off_diagonal_residual: float


def _four_port_dft(n: int = 4) -> np.ndarray:
    # Copied from `TH-17`'s four-port DFT helper (the unitary n-point DFT);
    # ``src`` never imports ``tests``.
    idx = np.arange(n)
    return np.exp(-2j * np.pi * np.outer(idx, idx) / n) / np.sqrt(n)


def _study_tune(self, k: int = 1, *, comm: MPI.Intracomm = MPI.COMM_WORLD,
                n_slots=None, at_hz=None) -> TuneResult:
    """Mode-``k`` capacitance ``C = 1/(omega Im z_k)`` from the Study's 4x4 S.

    ``Z = s_to_z(S, z0)``; ``z_k = diag(F Z F^H)`` with ``F`` the unitary 4-point
    DFT; ``omega`` is the Study's own ``2 pi frequency_hz``.  ``k`` is 1 or 2;
    ``Im z_k <= 0`` (always so for ``k = 0``) raises ``ValueError``.  Defined on a
    ``FourPortBirdcageGeometry`` Study with a ``LumpedPortDrive`` only; any other
    Study raises ``TypeError`` before meshing.  Reuses ``self.result`` and solves
    once if it is absent.  The read is lossless: a uniform series resistance on
    the four capacitors shifts ``Re z_k`` only, so ``C`` does not depend on it
    and no resistance argument is taken.

    ``n_slots`` (`WF-16` step 13): the F-human ring read only (default
    ``geometry.leg_count``); on any other route an explicit value raises
    ``TypeError`` before meshing.  An F-clinical Study takes the same ring read
    (`WF-16` step 16): loaded it needs no solve, unloaded it raises the solve's
    "xl tier only" ``ValueError``.

    ``at_hz`` (`WF-16` step 17): required on a multi-frequency loaded four-port
    result, where it selects the file's nearest frequency (no interpolation;
    outside the file's span raises ``ValueError``) and ``omega`` is that
    frequency's, not the Study's; on any other result an explicit value raises
    ``TypeError`` before meshing.
    """
    from .ports.circuit import s_to_z

    if (isinstance(self.geometry, (FHumanBirdcageGeometry, FClinicalBirdcageGeometry))
            and isinstance(self.excitation, LumpedPortDrive)):
        if at_hz is not None:  # WF-16 step 17 (§9 item 722)
            raise TypeError("at_hz: defined on a multi-frequency four-port result only")
        return _study_tune_ring(self, k, comm=comm, n_slots=n_slots)
    if n_slots is not None:  # WF-16 step 13 (§9 item 588)
        raise TypeError("n_slots: defined on the f_human_birdcage ring read only")
    _refuse_f_human(self, "tune")
    if not isinstance(self.geometry, FourPortBirdcageGeometry):
        raise TypeError("tune is defined on a four-port Study only (WF-16)")
    if not isinstance(self.excitation, LumpedPortDrive):
        raise TypeError("tune needs a LumpedPortDrive Study (WF-16)")
    if isinstance(k, bool) or not isinstance(k, (int, np.integer)):
        raise TypeError(f"k must be an int, got {type(k).__name__}")
    k = int(k)
    if k not in (0, 1, 2, 3):
        raise ValueError(f"k must be a mode index of the four-port, got {k}")
    loaded_multi = (isinstance(self.result, LoadedPortResult)  # WF-16 step 17 (§9 item 722)
                    and len(self.result.frequencies_hz) > 1)
    if at_hz is not None and not loaded_multi:
        raise TypeError("at_hz: defined on a multi-frequency loaded result only")
    if loaded_multi and at_hz is None:
        raise ValueError(f"tune on {self.result.source}: a multi-frequency result needs at_hz; "
                         f"it holds {len(self.result.frequencies_hz)} frequencies")
    if self.result is None:
        if k not in (1, 2):
            raise ValueError(f"tune is defined for k = 1 or 2, got k = {k}")
        self.solve(comm)
    frequency_hz = float(self.frequency_hz)
    if loaded_multi:
        freqs = np.asarray(self.result.frequencies_hz, dtype=np.float64)
        at = float(at_hz)
        if not freqs.min() <= at <= freqs.max():
            raise ValueError(f"at_hz {at!r} Hz is outside the file's span "
                             f"[{freqs.min()!r}, {freqs.max()!r}] Hz")
        col = int(np.argmin(np.abs(freqs - at)))
        s = np.asarray(self.result.s_matrices[col], dtype=np.complex128)
        frequency_hz = float(freqs[col])
        omega = 2.0 * np.pi * frequency_hz
    else:
        s = np.asarray(self.result.s_matrix, dtype=np.complex128)
        omega = self.omega
    z = s_to_z(s, float(self.excitation.z0_ohm))
    F = _four_port_dft(s.shape[0])
    m = F @ z @ F.conj().T
    zk = np.diag(m).copy()
    off = m - np.diag(zk)
    resid = float(np.linalg.norm(off) / np.linalg.norm(m))
    if k not in (1, 2):
        raise ValueError(f"tune is defined for k = 1 or 2, got k = {k} (z_k = {zk[k]!r})")
    if not zk[k].imag > 0.0:
        raise ValueError(f"Im z_{k} = {zk[k].imag!r} <= 0: no capacitive tuning")
    c = 1.0 / (omega * float(zk[k].imag))
    return TuneResult(k=k, frequency_hz=frequency_hz, z_k=zk,
                      c_mode_farad=float(c), off_diagonal_residual=resid)


def _study_tune_ring(self, k: int = 1, *, comm: MPI.Intracomm = MPI.COMM_WORLD,
                     n_slots=None) -> TuneResult:
    """`WF-16` step 12 (§9 item 472): the 32-port F-human ring read.

    ``Z = s_to_z(S, z0)`` in canonical (bottom 0…15, top 0…15) order from the
    port names ``P17…P48``; ``z_k`` = the ring-antisymmetric eigenvalue of the
    k-block of ``W Z Wᴴ`` (`ports.ring_modes`), ``C = 1/(omega Im z_k)``.
    Lossless, as the four-port read.

    ``n_slots`` (`WF-16` step 13, §9 item 588): ``None`` ⇒ the geometry's
    ``leg_count``; any other value must equal it (``ValueError`` naming
    ``n_slots``, before any solve).  An N-leg ring names its gaps
    P(N+1)…P(3N); ``z_k`` holds k = 0…N/2.
    """
    from .ports.circuit import s_to_z
    from .ports.ring_modes import (_canonical_permutation, _check_n_slots,
                                   antisymmetric_eigenpair, k_blocks,
                                   ring_port_map_from_names)

    if isinstance(k, bool) or not isinstance(k, (int, np.integer)):
        raise TypeError(f"k must be an int, got {type(k).__name__}")
    k = int(k)
    if k not in (1, 2):
        raise ValueError(f"tune is defined for k = 1 or 2, got k = {k}")
    n = self.geometry.leg_count if n_slots is None else n_slots
    n = _check_n_slots(n)
    if n != self.geometry.leg_count:
        raise ValueError(f"n_slots: {n} != the geometry's leg_count {self.geometry.leg_count}")
    if self.result is None:
        self.solve(comm)
    if isinstance(self.result, LoadedPortResult):  # WF-16 step 15 (§9 item 680)
        names = list(self.result.port_ids)
    else:
        names = [pd.port_id for pd in self.result.port_definitions]
    perm = _canonical_permutation(ring_port_map_from_names(names, n_slots=n), n_slots=n)
    s = np.asarray(self.result.s_matrix, dtype=np.complex128)
    z = s_to_z(s, float(self.excitation.z0_ohm))
    blocks, off = k_blocks(z[np.ix_(perm, perm)], n_slots=n)
    lam = np.array([antisymmetric_eigenpair(blocks[m])[0] for m in range(n // 2 + 1)])
    if not lam[k].imag > 0.0:
        raise ValueError(f"Im z_{k} = {lam[k].imag!r} <= 0: no capacitive tuning")
    return TuneResult(k=k, frequency_hz=float(self.frequency_hz), z_k=lam,
                      c_mode_farad=float(1.0 / (self.omega * lam[k].imag)),
                      off_diagonal_residual=float(off))


Study.tune = _study_tune


# ---------------------------------------------------------------------------
# `WF-16` step 14 (§9 item 661, §10 U94): a saved four-port result, loaded.
# Additive.  ``PortStudyResult`` requires every ``StudyResult`` field (mesh,
# cell tags, ...), so a loaded matrix is its own small class; ``tune`` reads
# only ``s_matrix`` and runs on it unchanged, every method that needs fields
# refuses it (``_refuse_loaded``).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class LoadedPortResult:
    """A port S matrix (four-port, the 32-port F-human ring since `WF-16`
    step 15, the 64-port F-clinical ring since step 16) read back from a Touchstone file (``Study.load_result``).

    No mesh and no fields: ``Study.tune`` reads it; ``sar`` / ``b1_plus`` /
    ``b1_plus_at`` / ``temperature`` / ``mean_sar`` refuse it.

    A multi-frequency four-port file (`WF-16` step 17) also fills
    ``frequencies_hz`` and ``s_matrices``; ``s_matrix`` and ``frequency_hz`` then
    hold the column nearest the Study's frequency."""

    s_matrix: Any
    frequency_hz: float
    port_ids: tuple
    source: str
    frequencies_hz: tuple = ()
    s_matrices: Any = None


def _refuse_loaded(study, method: str) -> None:
    """`WF-16` step 14: a field method on a loaded result raises (before any mesh)."""
    if isinstance(getattr(study, "result", None), LoadedPortResult):
        raise TypeError(f"{method} needs solved fields; this Study's result came from "
                        "load_result (an S matrix only) -- call solve() first")


#: Relative match of the file's single frequency to ``Study.frequency_hz``.
LOAD_RESULT_FREQUENCY_RTOL = 1e-12


def _study_load_result(self, path) -> LoadedPortResult:
    """Read a Touchstone file into ``self.result``; no solve.

    Defined on a ``FourPortBirdcageGeometry``, (`WF-16` step 15)
    ``FHumanBirdcageGeometry`` or (step 16) ``FClinicalBirdcageGeometry`` Study
    with a ``LumpedPortDrive`` only (``TypeError`` otherwise).  Refused (``ValueError``, nothing stored) unless
    the file holds exactly one frequency equal to ``frequency_hz`` to
    ``LOAD_RESULT_FREQUENCY_RTOL`` relative and the Study's port count:
    ``leg_count`` for the four-port, ``2 * leg_count`` (the ring-gap sheets)
    for the F-human and F-clinical.  Every rank reads the file.  ``tune`` then
    runs on the loaded matrix without solving; the ring read takes its port
    names from the file's ``port_ids``.

    `WF-16` step 17: on a ``FourPortBirdcageGeometry`` a multi-frequency file
    also loads, provided ``frequency_hz`` lies inside its span; ``tune`` then
    needs ``at_hz``.  A multi-frequency ring file stays refused.
    """
    from pathlib import Path

    from .ports.touchstone import load_touchstone

    if not isinstance(self.geometry, (FourPortBirdcageGeometry, FHumanBirdcageGeometry,
                                      FClinicalBirdcageGeometry)):
        raise TypeError("load_result is defined on a FourPortBirdcageGeometry, "
                        "FHumanBirdcageGeometry or FClinicalBirdcageGeometry Study only "
                        "(WF-16 steps 14-16)")
    if not isinstance(self.excitation, LumpedPortDrive):
        raise TypeError("load_result needs a LumpedPortDrive Study (WF-16 step 14)")
    path = Path(path)
    freqs, s_matrices, port_ids = load_touchstone(path)
    freqs = np.asarray(freqs, dtype=np.float64)
    multi = freqs.size > 1 and isinstance(self.geometry, FourPortBirdcageGeometry)
    if freqs.shape != (1,) and not multi:
        raise ValueError(f"load_result reads exactly one frequency; {path.name} holds "
                         f"{freqs.size}")
    if multi:  # WF-16 step 17 (§9 item 722)
        if not freqs.min() <= self.frequency_hz <= freqs.max():
            raise ValueError(f"load_result: the Study's frequency_hz {self.frequency_hz!r} is "
                             f"outside {path.name}'s span [{freqs.min()!r}, {freqs.max()!r}] Hz")
        col = int(np.argmin(np.abs(freqs - self.frequency_hz)))
    else:
        col = 0
    f = float(freqs[col])
    if not multi and (abs(f - self.frequency_hz)
                      > LOAD_RESULT_FREQUENCY_RTOL * abs(self.frequency_hz)):
        raise ValueError(f"load_result: file frequency {f!r} Hz != the Study's "
                         f"frequency_hz {self.frequency_hz!r}")
    n = int(self.geometry.leg_count)
    if isinstance(self.geometry, (FHumanBirdcageGeometry, FClinicalBirdcageGeometry)):  # WF-16 steps 15-16
        n = 2 * n
    if len(port_ids) != n or s_matrices.shape[1:] != (n, n):
        raise ValueError(f"load_result: {len(port_ids)} ports in {path.name}, the Study "
                         f"defines {n}")
    self.result = LoadedPortResult(
        s_matrix=np.asarray(s_matrices[col], dtype=np.complex128).copy(),
        frequency_hz=f,
        port_ids=tuple(port_ids),
        source=str(path),
        frequencies_hz=tuple(float(x) for x in freqs) if multi else (),
        s_matrices=np.asarray(s_matrices, dtype=np.complex128).copy() if multi else None,
    )
    return self.result


Study.load_result = _study_load_result


# ---------------------------------------------------------------------------
# `WF-15` step 2 (§9 item 260): steady temperature through the Study.
# Additive.  The four private helpers are `MAT-10`'s
# (tests/validation/test_mat10_sar_to_heat_sphere.py:74–130), copied verbatim
# in operation order and parameterised only by the tag (``src`` never imports
# ``tests``).
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class RobinCooling:
    """Steady P-``degree`` conduction with Robin cooling on the sphere surface.

    ``source_scale`` multiplies q = ½σ|E|² only (for an imposed amplitude E′,
    (E′/e0)²); it never touches the EM solve.  ``degree`` is the thermal
    Lagrange degree, not the Study's EM degree.
    """

    k_w_per_m_k: float
    h_w_per_m2_k: float
    t_inf_c: float
    source_scale: float = 1.0
    degree: int = 1

    def __post_init__(self):
        for name in ("k_w_per_m_k", "h_w_per_m2_k", "source_scale"):
            if not float(getattr(self, name)) > 0.0:
                raise ValueError(f"{name} must be positive, got {getattr(self, name)!r}")
        if int(self.degree) < 1:
            raise ValueError(f"degree must be >= 1, got {self.degree!r}")


@dataclass
class ThermalStudyResult:
    """What ``Study.temperature`` returns (not exported)."""

    temperature: Any
    submesh: Any
    parent_cells: Any
    facet_tags: Any
    heat_source: Any
    source_power: float
    convected_power: float
    cooling: RobinCooling


_ROBIN_TAG = 1


def _tag_submesh(msh, cell_tags, tag):
    """Cell tag ``tag`` of the parent as its own mesh (`MAT-10` :74–89)."""
    from dolfinx import mesh as dmesh

    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_parent = imap.size_local + imap.num_ghosts
    parent_values = np.full(n_parent, -1, dtype=np.int32)
    parent_values[cell_tags.indices] = cell_tags.values
    keep = np.flatnonzero(parent_values == tag).astype(np.int32)
    sub, emap, _vmap, _gmap = dmesh.create_submesh(msh, tdim, keep)
    simap = sub.topology.index_map(tdim)
    n_sub = simap.size_local + simap.num_ghosts
    parent_ids = np.asarray(
        emap.sub_topology_to_topology(np.arange(n_sub, dtype=np.int32), False),
        dtype=np.int32,
    )
    return sub, parent_ids


def _exterior_robin_tags(sub):
    """Every exterior facet of the submesh -> ``_ROBIN_TAG`` (`MAT-10` :92–98)."""
    from dolfinx import mesh as dmesh

    tdim = sub.topology.dim
    fdim = tdim - 1
    sub.topology.create_connectivity(fdim, tdim)
    ext = np.sort(np.asarray(dmesh.exterior_facet_indices(sub.topology), dtype=np.int32))
    return dmesh.meshtags(sub, fdim, ext, np.full(ext.size, _ROBIN_TAG, dtype=np.int32))


def _dg0_cell_average(msh, cell_tags, integrand, tag):
    """Parent DG0 cell average of ``integrand`` over ``tag`` (`MAT-10` :101–119)."""
    import dolfinx
    import ufl
    from dolfinx import default_scalar_type, fem

    W = fem.functionspace(msh, ("DG", 0))
    w = ufl.TestFunction(W)
    dx1 = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)(tag)
    one = fem.Constant(msh, default_scalar_type(1.0))

    def _vec(form):
        b = fem.assemble_vector(fem.form(form))
        b.scatter_reverse(dolfinx.la.InsertMode.add)
        b.scatter_forward()
        return b.array.copy()

    num = _vec(ufl.inner(integrand, w) * dx1)
    vol = _vec(ufl.inner(one, w) * dx1)
    avg = np.zeros_like(num)
    mask = np.abs(vol) > 0.0
    avg[mask] = num[mask] / vol[mask]
    return W, avg


def _dg0_to_submesh(W_parent, avg_parent, sub, parent_ids):
    """Map a parent DG0 array to a sub DG0 Function (`MAT-10` :122–130)."""
    from dolfinx import fem

    Ws = fem.functionspace(sub, ("DG", 0))
    q = fem.Function(Ws)
    pdofs = W_parent.dofmap.list[parent_ids, 0]
    sdofs = Ws.dofmap.list[np.arange(parent_ids.size), 0]
    q.x.array[sdofs] = avg_parent[pdofs]
    q.x.scatter_forward()
    return q


def _study_temperature(self, cooling: RobinCooling, *, transient=None, duty=None,
                       balls=None, window_s=None, sense=None, dose=None,
                       drive: Optional[str] = None,
                       comm: MPI.Intracomm = MPI.COMM_WORLD):
    """Steady temperature of the sphere heated by ``source_scale``·½σ|E|².

    Reuses ``self.result`` (solves once if absent, as ``mean_sar`` does).

    On a thermal-only Study (``ThermalBallGeometry``, `WF-15` step 3) a
    ``TransientSpec`` is required and the result is a
    ``TransientStudyResult``; ``duty`` (a ``DutyCycle``) becomes the
    transient solver's ``source_schedule``.

    On the EM sphere (``SphereInBoxGeometry``) a ``TransientSpec`` gives
    `MAT-12`'s SAR-driven θ-scheme transient on the same submesh, Robin facets
    and DG0 source as the steady branch, returned as a ``TransientStudyResult``
    (`WF-15` step 15). EM degree 1 and the default degree 2 are gated there
    (`WF-15` steps 15–16; degree 1 reproduces `MAT-12`'s records, degree 2
    is a registered record against the same closed form). A ``DutyCycle``
    is gated on the EM sphere at EM degree 1 (`WF-15` step 17); ``duty``
    without ``transient`` and ``initial_gradient_k_per_m`` are refused there.

    On a ``TwoPortQuadratureDrive`` Study (`WF-1` step 10a) the excitation's
    hook gives a transient temperature of the phantom with 10 g ball means
    (``balls``), the running maximum over ``window_s`` and the ``sense``
    field's source; ``balls`` / ``window_s`` / ``sense`` are refused elsewhere.
    ``duty`` is refused there (`WF-15` step 19).

    ``dose`` (a ``DoseReadout``, `WF-1` step 30, §9 item 878): on a transient
    ``ThermalBoxGeometry`` with a uniform source and no duty, the temperature
    history at the readout's points with `MAT-20`'s ``cem43`` and trailing
    window mean per point, on ``TransientStudyResult.dose``; refused elsewhere.

    ``drive`` (a port name, `WF-15` step 14, §9 items 831 / 884): on the
    implant-rod ``FourPortBirdcageGeometry`` the transient temperature of the
    phantom under that one port at 1 V, cooled on the outer boundary only (the
    rod wall adiabatic), with 1 g ball means; refused elsewhere.
    """
    hook = getattr(self.excitation, "_temperature_port_study", None)
    if dose is not None and (hook is not None
                             or not isinstance(self.geometry, ThermalBoxGeometry)):
        raise ValueError("dose is defined on a transient ThermalBoxGeometry only (WF-1 step 30)")
    rod_route = (isinstance(self.geometry, FourPortBirdcageGeometry)
                 and bool(self.geometry.implant_rod))
    if drive is not None and (hook is not None or not rod_route):
        raise ValueError("drive is defined on the implant-rod four-port Study only (WF-15 step 14)")
    if hook is None and rod_route:
        # `WF-15` step 14 (§9 item 831): degree / caps refusals before any mesh.
        _refuse_rod_route(self, "temperature")
        if duty is not None or window_s is not None or sense is not None:
            raise ValueError("duty / window_s / sense are not defined on the rod four-port Study")
        return _rod_port_temperature(self, cooling, transient=transient, balls=balls,
                                     drive=drive, comm=comm)
    if hook is not None:
        if duty is not None:
            raise TypeError("duty is not defined on a tuned two-port Study (WF-15 step 19)")
        kw = {} if sense is None else {"sense": sense}
        return hook(self, cooling, transient=transient, balls=balls, window_s=window_s,
                    comm=comm, **kw)
    if balls is not None or window_s is not None or sense is not None:
        raise ValueError("balls / window_s / sense are defined on a tuned two-port Study only")
    import ufl

    from .thermal.steady import solve_steady_conduction

    if isinstance(self.geometry, ThermalBoxGeometry):
        if dose is not None:  # `WF-1` step 30 (§9 item 878)
            return _study_box_temperature(self, cooling, transient, duty, dose=dose)
        return _study_box_temperature(self, cooling, transient, duty)
    if isinstance(self.geometry, ThermalBallGeometry):
        return _study_transient_temperature(self, cooling, transient, duty)
    _refuse_f_human(self, "temperature")
    if (transient is not None or duty is not None) and not isinstance(
            self.geometry, (SphereInBoxGeometry, FourPortBirdcageGeometry)):
        raise ValueError("transient / duty are defined on a thermal-only Study in v0")
    if isinstance(self.geometry, FourPortBirdcageGeometry):
        raise TypeError("temperature is not defined on a port Study (WF-16)")
    if not isinstance(cooling, RobinCooling):
        raise TypeError(f"cooling must be a RobinCooling, got {type(cooling).__name__}")
    _refuse_loaded(self, "temperature")  # WF-16 step 14 (§9 item 661)
    # `WF-15` step 15 (§9 item 862): `MAT-12`'s SAR-driven transient on the EM
    # sphere — degree-1 EM is what is gated; an initial gradient is not.  A
    # ``DutyCycle`` is gated at EM degree 1 (`WF-15` step 17, §9 item 877).
    if transient is None and duty is not None:
        raise ValueError("duty needs a transient")
    if duty is not None and not isinstance(duty, DutyCycle):
        raise TypeError(f"duty must be a DutyCycle, got {type(duty).__name__}")
    if transient is not None:
        if not isinstance(transient, TransientSpec):
            raise TypeError(f"transient must be a TransientSpec, got {type(transient).__name__}")
        if transient.initial_gradient_k_per_m is not None:
            raise ValueError("initial_gradient_k_per_m is defined on a ThermalBoxGeometry only")
    if self.result is None:
        self.solve()
    res = self.result
    tag = self.geometry.sphere_tag
    sub, parent_ids = _tag_submesh(res.mesh, res.cell_tags, tag)
    ft = _exterior_robin_tags(sub)
    E, sigma = res.fields.e_complex, res.fields.sigma_field
    q_expr = cooling.source_scale * 0.5 * sigma * ufl.inner(E, E)
    Wp, q_par = _dg0_cell_average(res.mesh, res.cell_tags, q_expr, tag)
    q_h = _dg0_to_submesh(Wp, q_par, sub, parent_ids)
    if transient is not None:
        return _sphere_transient_temperature(cooling, transient, sub, ft, q_h, duty=duty)
    th = solve_steady_conduction(
        sub,
        ft,
        k=cooling.k_w_per_m_k,
        q=q_h,
        h=cooling.h_w_per_m2_k,
        t_inf=cooling.t_inf_c,
        robin_tag=_ROBIN_TAG,
        degree=cooling.degree,
    )
    return ThermalStudyResult(
        temperature=th.temperature,
        submesh=sub,
        parent_cells=parent_ids,
        facet_tags=ft,
        heat_source=q_h,
        source_power=th.source_power,
        convected_power=th.convected_power,
        cooling=cooling,
    )


def _sphere_transient_temperature(cooling, transient, sub, ft, q_h, duty=None):
    """`WF-15` step 15: the ball route's transient call on the EM sphere's submesh."""
    from .thermal.transient import solve_transient_conduction

    t_inf = float(cooling.t_inf_c)
    t_init = t_inf if transient.t_init_c is None else float(transient.t_init_c)
    s = solve_transient_conduction(
        sub,
        ft,
        k=float(cooling.k_w_per_m_k),
        rho_c=float(transient.rho_c_j_per_m3_k),
        q=q_h,
        h=float(cooling.h_w_per_m2_k),
        t_inf=t_inf,
        t_init=t_init,
        robin_tag=_ROBIN_TAG,
        dt=float(transient.dt_s),
        n_steps=int(transient.n_steps),
        theta=float(transient.theta),
        startup_half_steps=int(transient.startup_half_steps),
        degree=int(cooling.degree),
        source_schedule=duty,
    )
    return TransientStudyResult(
        temperature=s.temperature,
        mesh=sub,
        cell_tags=None,
        facet_tags=ft,
        k=float(cooling.k_w_per_m_k),
        rho_c=float(transient.rho_c_j_per_m3_k),
        q=q_h,
        source_multipliers=list(s.source_multipliers),
        source_power=s.source_power,
        max_balance_residual=s.max_balance_residual,
        solver_result=s,
        transient=transient,
        duty=duty,
        cooling=cooling,
    )


Study.temperature = _study_temperature


# ---------------------------------------------------------------------------
# `WF-15` step 3 (§9 item 330, chain P16): ``Study.temperature`` with a duty
# cycle and a container wall, on a thermal-only ball with an imposed uniform
# source (no EM solve).  The records reproduced are `MAT-17` step 2 (duty)
# and `MAT-18` step 2(b) (wall); the solver is ``solve_transient_conduction``
# unchanged, the mesh ``MeshGenerator.ball_domain`` / ``two_layer_ball_domain``.
# ---------------------------------------------------------------------------


def _positive(obj, *names):
    for name in names:
        v = getattr(obj, name)
        if not float(v) > 0.0:
            raise ValueError(f"{name} must be positive, got {v!r}")


@dataclass(frozen=True)
class ContainerWall:
    """A shell ``radius < r < radius + thickness`` around the ball, no source."""

    thickness_m: float
    k_w_per_m_k: float
    rho_c_j_per_m3_k: float

    def __post_init__(self):
        _positive(self, "thickness_m", "k_w_per_m_k", "rho_c_j_per_m3_k")


@dataclass(frozen=True)
class ThermalBallGeometry:
    """A ball of ``radius`` (cell tag 1), optionally inside a ``ContainerWall``
    (cell tag 2); facet tag 1 is the outer surface (Robin cooling)."""

    radius: float
    resolution: float
    wall: Optional[ContainerWall] = None

    def __post_init__(self):
        _positive(self, "radius", "resolution")
        if self.wall is not None and not isinstance(self.wall, ContainerWall):
            raise TypeError(f"wall must be a ContainerWall, got {type(self.wall).__name__}")

    @property
    def outer_radius(self) -> float:
        return self.radius if self.wall is None else self.radius + self.wall.thickness_m

    def build(self, comm: MPI.Intracomm = MPI.COMM_WORLD):
        """Return ``(mesh, cell_tags, facet_tags)`` from the generator."""
        if self.wall is None:
            return MeshGenerator.ball_domain(
                radius=self.radius, resolution=self.resolution, comm=comm
            )
        return MeshGenerator.two_layer_ball_domain(
            inner_radius=self.radius,
            outer_radius=self.outer_radius,
            resolution=self.resolution,
            comm=comm,
        )


@dataclass(frozen=True)
class ThermalBoxGeometry:
    """`WF-15` step 5 (§9 item 521a): the structured ``create_box`` box
    [-A, A] x [-B, B] x [-C, C], ``divisions`` hexahedra split into
    tetrahedra (`MAT-21`'s ``_rung``); cell tag 1 everywhere, facet tag 1
    (Robin cooling) on all six faces.

    `WF-1` step 28 (§9 item 830): with ``wall`` (a ``ContainerWall``) and
    ``wall_divisions`` N_w >= 1, the box is `MAT-21` P48b's 1D-reduced walled
    box — the wall on ±z only (cell tag 2, ``thickness_m`` thick, N_w cells
    across), the gel cell tag 1, and facet tag 1 (Robin) on z = ±(C + d) only;
    the four side faces are natural Neumann."""

    half_widths_m: tuple
    divisions: tuple
    wall: Optional[ContainerWall] = None
    wall_divisions: Optional[int] = None

    def __post_init__(self):
        hw = tuple(self.half_widths_m)
        dv = tuple(self.divisions)
        if len(hw) != 3 or not all(np.isfinite(float(v)) and float(v) > 0.0 for v in hw):
            raise ValueError(f"half_widths_m must be three finite positive floats, got {hw!r}")
        if len(dv) != 3 or not all(int(v) == v and int(v) >= 1 for v in dv):
            raise ValueError(f"divisions must be three integers >= 1, got {dv!r}")
        object.__setattr__(self, "half_widths_m", tuple(float(v) for v in hw))
        object.__setattr__(self, "divisions", tuple(int(v) for v in dv))
        if self.wall is None:
            if self.wall_divisions is not None:
                raise ValueError("wall_divisions is defined only with a wall")
        else:
            if not isinstance(self.wall, ContainerWall):
                raise TypeError(f"wall must be a ContainerWall, got {type(self.wall).__name__}")
            nw = self.wall_divisions
            if nw is None or isinstance(nw, bool) or int(nw) != nw or int(nw) < 1:
                raise ValueError(f"wall_divisions must be an integer >= 1 with a wall, got {nw!r}")
            object.__setattr__(self, "wall_divisions", int(nw))

    def build(self, comm: MPI.Intracomm = MPI.COMM_WORLD):
        """Return ``(mesh, cell_tags, facet_tags)``."""
        from dolfinx import mesh as dmesh

        a, b, c = self.half_widths_m
        if self.wall is not None:
            return self._build_walled(comm)
        msh = dmesh.create_box(comm, [np.array([-a, -b, -c]), np.array([a, b, c])],
                               list(self.divisions), cell_type=dmesh.CellType.tetrahedron)
        imap = msh.topology.index_map(3)
        cells = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
        cell_tags = dmesh.meshtags(msh, 3, cells, np.ones(cells.size, dtype=np.int32))
        return msh, cell_tags, _exterior_robin_tags(msh)

    def _build_walled(self, comm):
        """`MAT-21` P48b's walled box (``tests/validation/test_mat21_wall_fem.py:91–104``,
        ``_walled_box(m, (False, False, True))``), its cell midpoint test (``:107–124``)
        and its ``_robin_facets(..., (2,))`` marker (``:128–139``), copied."""
        from dolfinx import mesh as dmesh

        a, b, c = self.half_widths_m
        d = float(self.wall.thickness_m)
        nw = int(self.wall_divisions)
        nx, ny, nz = self.divisions
        outer = (a, b, c + d)
        msh = dmesh.create_box(comm, [np.array([-o for o in outer]), np.array(outer)],
                               [nx, ny, nz + 2 * nw], cell_type=dmesh.CellType.tetrahedron)
        x = msh.geometry.x  # owned + ghost nodes on this rank
        x[:, 2] = _box_wall_remap(x[:, 2].copy(), c, d, nz, nw)
        tdim = msh.topology.dim
        imap = msh.topology.index_map(tdim)
        n = imap.size_local + imap.num_ghosts
        cells = np.arange(n, dtype=np.int32)
        mids = dmesh.compute_midpoints(msh, tdim, cells)
        gel = (np.abs(mids[:, 0]) < a) & (np.abs(mids[:, 1]) < b) & (np.abs(mids[:, 2]) < c)
        cell_tags = dmesh.meshtags(msh, tdim, cells, np.where(gel, 1, 2).astype(np.int32))
        fdim = tdim - 1
        msh.topology.create_connectivity(fdim, tdim)

        def mark(xx):
            return np.isclose(np.abs(xx[2]), outer[2], rtol=0.0, atol=1e-12)

        f = np.sort(dmesh.locate_entities_boundary(msh, fdim, mark).astype(np.int32))
        ft = dmesh.meshtags(msh, fdim, f, np.full(f.size, _ROBIN_TAG, dtype=np.int32))
        return msh, cell_tags, ft


def _box_wall_remap(s, L, d, n_gel, n_w):
    """Uniform grid on [−(L+d), L+d] (n_gel + 2 n_w cells) → planes at ±L.

    Copied from ``tests/validation/test_mat21_wall_fem.py:78–88`` (`MAT-21` step 20,
    §9 item 670); ``src`` imports nothing from ``tests``."""
    if n_w == 0:
        return s
    tot = n_gel + 2 * n_w
    u = (s + (L + d)) / (2.0 * (L + d)) * tot  # grid-index coordinate
    out = np.where(
        u <= n_w, -(L + d) + u * (d / n_w),
        np.where(u <= n_w + n_gel, -L + (u - n_w) * (2.0 * L / n_gel),
                 L + (u - n_w - n_gel) * (d / n_w)))
    return out


def _study_box_temperature(self, cooling, transient, duty, dose=None):
    """Steady (``transient=None``) or transient temperature of the box."""
    if dose is not None:  # `WF-1` step 30 (§9 item 878)
        if not isinstance(dose, DoseReadout):
            raise TypeError(f"dose must be a DoseReadout, got {type(dose).__name__}")
        if transient is None:
            raise ValueError("dose needs a transient (WF-1 step 30)")
        if isinstance(self.excitation, TransverseQuadraticHeatSource):
            raise ValueError("dose is defined with a UniformHeatSource only (WF-1 step 30)")
        if duty is not None:
            raise ValueError("dose is not admitted beside a duty (WF-1 step 30)")
    if self.geometry.wall is not None:  # `WF-1` step 28 (§9 item 830)
        if isinstance(self.excitation, TransverseQuadraticHeatSource):
            raise ValueError("a walled box admits a UniformHeatSource only")
        if transient is None:
            raise ValueError("a walled box needs a transient (no steady walled-box record)")
        if duty is not None:
            raise ValueError("duty is not admitted beside a box wall")
    if isinstance(self.excitation, TransverseQuadraticHeatSource):
        # `WF-1` step 16 / `WF-15` step 7 (§9 item 528).
        return _study_box_quadratic_temperature(self, cooling, transient, duty)
    if transient is None and duty is not None:
        raise ValueError("duty needs a transient")
    if transient is not None:
        if dose is not None:  # `WF-1` step 30 (§9 item 878)
            return _study_transient_temperature(self, cooling, transient, duty, dose=dose)
        return _study_transient_temperature(self, cooling, transient, duty)
    from .thermal.steady import solve_steady_conduction

    if not isinstance(cooling, RobinCooling):
        raise TypeError(f"cooling must be a RobinCooling, got {type(cooling).__name__}")
    if float(cooling.source_scale) != 1.0:
        raise ValueError("source_scale applies to an EM source; scale UniformHeatSource.q instead")
    msh, _cell_tags, ft = self.geometry.build(MPI.COMM_WORLD)
    q = float(self.excitation.q_w_per_m3)
    th = solve_steady_conduction(
        msh, ft, k=cooling.k_w_per_m_k, q=q, h=cooling.h_w_per_m2_k, t_inf=cooling.t_inf_c,
        robin_tag=_ROBIN_TAG, degree=cooling.degree,
    )
    return ThermalStudyResult(
        temperature=th.temperature,
        submesh=msh,
        parent_cells=None,
        facet_tags=ft,
        heat_source=q,
        source_power=th.source_power,
        convected_power=th.convected_power,
        cooling=cooling,
    )


@dataclass(frozen=True)
class UniformHeatSource:
    """A uniform volumetric heat source ``q`` in the ball (W/m³); the wall has none."""

    q_w_per_m3: float

    def __post_init__(self):
        _positive(self, "q_w_per_m3")


@dataclass(frozen=True)
class TransverseQuadraticHeatSource:
    """`WF-15` step 7 (§9 item 528): q(x) = 3 q (y² + z²)/(B² + C²) on a
    `ThermalBoxGeometry` with half-widths (A, B, C); its box mean is q
    (`MAT-21` step 9's transverse eddy-SAR shape).  Box Studies only."""

    q_w_per_m3: float

    def __post_init__(self):
        _positive(self, "q_w_per_m3")


def _study_box_quadratic_temperature(self, cooling, transient, duty):
    """The box under ``TransverseQuadraticHeatSource``: steady or transient,
    the source a UFL expression in P36's operation order (`MAT-21` step 9,
    ``tests/validation/test_mat21_box_quadratic_fem.py:141–142``)."""
    import ufl

    from .thermal.steady import solve_steady_conduction
    from .thermal.transient import solve_transient_conduction

    if transient is None and duty is not None:
        raise ValueError("duty needs a transient")
    if not isinstance(cooling, RobinCooling):
        raise TypeError(f"cooling must be a RobinCooling, got {type(cooling).__name__}")
    if transient is not None and not isinstance(transient, TransientSpec):
        raise TypeError("a thermal-only Study needs transient=TransientSpec(...)")
    if transient is not None and transient.initial_gradient_k_per_m is not None:
        raise ValueError("initial_gradient_k_per_m is defined on a uniform-source box only")
    if duty is not None and not isinstance(duty, DutyCycle):
        raise TypeError(f"duty must be a DutyCycle, got {type(duty).__name__}")
    if float(cooling.source_scale) != 1.0:
        raise ValueError("source_scale applies to an EM source; scale the heat source's q instead")
    geo = self.geometry
    msh, cell_tags, facet_tags = geo.build(MPI.COMM_WORLD)
    q0 = float(self.excitation.q_w_per_m3)
    _a, B, C = geo.half_widths_m
    xs = ufl.SpatialCoordinate(msh)
    q = q0 * 3.0 * (xs[1] * xs[1] + xs[2] * xs[2]) / (B**2 + C**2)
    if transient is None:
        th = solve_steady_conduction(
            msh, facet_tags, k=cooling.k_w_per_m_k, q=q, h=cooling.h_w_per_m2_k,
            t_inf=cooling.t_inf_c, robin_tag=_ROBIN_TAG, degree=cooling.degree,
        )
        return ThermalStudyResult(
            temperature=th.temperature,
            submesh=msh,
            parent_cells=None,
            facet_tags=facet_tags,
            heat_source=q,
            source_power=th.source_power,
            convected_power=th.convected_power,
            cooling=cooling,
        )
    k, rho_c = float(cooling.k_w_per_m_k), float(transient.rho_c_j_per_m3_k)
    t_inf = float(cooling.t_inf_c)
    t_init = t_inf if transient.t_init_c is None else float(transient.t_init_c)
    s = solve_transient_conduction(
        msh,
        facet_tags,
        k=k,
        rho_c=rho_c,
        q=q,
        h=float(cooling.h_w_per_m2_k),
        t_inf=t_inf,
        t_init=t_init,
        robin_tag=_ROBIN_TAG,
        dt=float(transient.dt_s),
        n_steps=int(transient.n_steps),
        theta=float(transient.theta),
        startup_half_steps=int(transient.startup_half_steps),
        degree=int(cooling.degree),
        source_schedule=duty,
    )
    return TransientStudyResult(
        temperature=s.temperature,
        mesh=msh,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        k=k,
        rho_c=rho_c,
        q=q,
        source_multipliers=list(s.source_multipliers),
        source_power=s.source_power,
        max_balance_residual=s.max_balance_residual,
        solver_result=s,
        transient=transient,
        duty=duty,
        cooling=cooling,
    )


@dataclass(frozen=True)
class TransientSpec:
    """θ-scheme time stepping from ``T(0) = t_init_c`` (default: ``t_inf_c``).

    ``rho_c_j_per_m3_k`` is the ball's volumetric heat capacity; ``theta`` 0.5
    with ``startup_half_steps`` 2 is Crank–Nicolson with a Rannacher start.

    ``initial_gradient_k_per_m`` g (`WF-1` step 28, §9 item 830): T(0) =
    t_init + g·z on every cell (`MAT-21` P47); a ``ThermalBoxGeometry`` with a
    uniform source only, refused elsewhere.
    """

    dt_s: float
    n_steps: int
    rho_c_j_per_m3_k: float
    theta: float = 0.5
    startup_half_steps: int = 2
    t_init_c: Optional[float] = None
    initial_gradient_k_per_m: Optional[float] = None

    def __post_init__(self):
        _positive(self, "dt_s", "rho_c_j_per_m3_k")
        g = self.initial_gradient_k_per_m
        if g is not None and not (isinstance(g, (int, float)) and not isinstance(g, bool)
                                  and np.isfinite(float(g))):
            raise ValueError(f"initial_gradient_k_per_m must be a finite number, got {g!r}")
        if int(self.n_steps) < 1:
            raise ValueError(f"n_steps must be >= 1, got {self.n_steps!r}")
        if not 0.0 <= float(self.theta) <= 1.0:
            raise ValueError(f"theta must be in [0, 1], got {self.theta!r}")
        if int(self.startup_half_steps) < 0 or int(self.startup_half_steps) % 2:
            raise ValueError("startup_half_steps must be a non-negative even integer")


@dataclass(frozen=True)
class DoseReadout:
    """`WF-1` step 30 (§9 item 878): a point-history readout of `MAT-20`'s
    thermal dose on a transient ``ThermalBoxGeometry`` Study.

    ``points_m`` is a non-empty tuple of finite (x, y, z) points inside the gel
    box; at each, ``Study.temperature(…, dose=…)`` returns ``cem43`` of the
    temperature history on the full-step grid and the trailing ``window_s``
    running average (``thermal.dose.cem43`` / ``running_average``).  Readout
    only: no dose limit or threshold is implied.

    ``window_s`` defaults to 360.0 s, ``running_average``'s own default
    (``thermal/dose.py:108``) and ``_SIX_MINUTE_WINDOW_S`` below; written as a
    literal because that name is defined after this class.
    """

    points_m: tuple
    window_s: float = 360.0

    def __post_init__(self):
        pts = self.points_m
        if not isinstance(pts, tuple) or len(pts) == 0:
            raise ValueError(f"points_m must be a non-empty tuple of (x, y, z), got {pts!r}")
        for p in pts:
            if (not isinstance(p, tuple) or len(p) != 3
                    or not all(isinstance(c, (int, float)) and not isinstance(c, bool)
                               and np.isfinite(float(c)) for c in p)):
                raise ValueError(f"each dose point must be a finite (x, y, z) tuple, got {p!r}")
        w = self.window_s
        if not (isinstance(w, (int, float)) and not isinstance(w, bool)
                and np.isfinite(float(w)) and float(w) > 0.0):
            raise ValueError(f"window_s must be finite and > 0, got {w!r}")
        # Plain floats (an ``np.float64`` coordinate is stored as its float value,
        # same bits), so the repr does not depend on the caller's array type.
        object.__setattr__(self, "points_m", tuple(tuple(float(c) for c in p) for p in pts))
        object.__setattr__(self, "window_s", float(w))


@dataclass(frozen=True)
class DutyCycle:
    """Source on while ``0 <= (t mod period) < duty * period`` (`MAT-17`)."""

    period_s: float
    duty: float

    def __post_init__(self):
        _positive(self, "period_s")
        if not 0.0 < float(self.duty) <= 1.0:
            raise ValueError(f"duty must be in (0, 1], got {self.duty!r}")

    def __call__(self, t: float) -> float:
        return 1.0 if (t % self.period_s) < self.duty * self.period_s else 0.0


@dataclass
class TransientStudyResult:
    """What ``Study.temperature`` returns on a thermal-only Study, or the EM sphere
    with ``transient`` (`WF-15` step 15; ``cell_tags`` is None there; EM degree 1
    and the default degree 2 gated, `WF-15` steps 15–16; a ``DutyCycle`` is
    gated on the EM sphere at EM degree 1, `WF-15` step 17) (not exported).

    ``dose`` (`WF-1` step 30, `MAT-20`): with a ``DoseReadout`` on a box, plain
    floats — ``points_m``, ``times_s``, ``temperature_c`` (n_t × n_pts),
    ``cem43_s`` and ``six_minute`` {``window_s``, ``peak_c``, ``t_max_s``} per
    point; otherwise None."""

    temperature: Any
    mesh: Any
    cell_tags: Any
    facet_tags: Any
    k: Any
    rho_c: Any
    q: Any
    source_multipliers: list
    source_power: float
    max_balance_residual: float
    solver_result: Any
    transient: TransientSpec
    duty: Optional[DutyCycle]
    cooling: RobinCooling
    dose: Optional[dict] = None


def _dg0_two_region(msh, cell_tags, v_core, v_shell):
    """DG0 coefficient: ``v_core`` on tag 1, ``v_shell`` on tag 2, ghost-updated."""
    from dolfinx import fem

    W = fem.functionspace(msh, ("DG", 0))
    f = fem.Function(W)
    f.x.array[:] = np.nan
    f.x.array[cell_tags.find(1)] = v_core
    f.x.array[cell_tags.find(2)] = v_shell
    f.x.scatter_forward()
    return f


def _study_transient_temperature(self, cooling, transient, duty,
                                 dose=None) -> TransientStudyResult:
    from .thermal.transient import solve_transient_conduction

    if dose is not None:  # `WF-1` step 30 (§9 item 878): before any mesh
        if not isinstance(self.geometry, ThermalBoxGeometry):
            raise ValueError("dose is defined on a transient ThermalBoxGeometry only (WF-1 step 30)")
        half = [float(v) for v in self.geometry.half_widths_m]
        for i, p in enumerate(dose.points_m):
            if any(abs(float(c)) > h for c, h in zip(p, half)):
                raise ValueError(f"dose point {i} {p!r} lies outside the gel box "
                                 f"[-A, A] x [-B, B] x [-C, C] = {tuple(half)!r}")

    if not isinstance(cooling, RobinCooling):
        raise TypeError(f"cooling must be a RobinCooling, got {type(cooling).__name__}")
    if not isinstance(transient, TransientSpec):
        raise TypeError("a thermal-only Study needs transient=TransientSpec(...)")
    if duty is not None and not isinstance(duty, DutyCycle):
        raise TypeError(f"duty must be a DutyCycle, got {type(duty).__name__}")
    if float(cooling.source_scale) != 1.0:
        raise ValueError("source_scale applies to an EM source; scale UniformHeatSource.q instead")
    geo = self.geometry
    g = transient.initial_gradient_k_per_m
    if g is not None and not isinstance(geo, ThermalBoxGeometry):
        raise ValueError("initial_gradient_k_per_m is defined on a ThermalBoxGeometry only")
    msh, cell_tags, facet_tags = geo.build(MPI.COMM_WORLD)
    q0 = float(self.excitation.q_w_per_m3)
    if geo.wall is None:
        # Scalars: the unchanged ``fem.Constant`` path of `MAT-17` step 2.
        k, rho_c, q = float(cooling.k_w_per_m_k), float(transient.rho_c_j_per_m3_k), q0
    else:
        w = geo.wall
        k = _dg0_two_region(msh, cell_tags, float(cooling.k_w_per_m_k), float(w.k_w_per_m_k))
        rho_c = _dg0_two_region(
            msh, cell_tags, float(transient.rho_c_j_per_m3_k), float(w.rho_c_j_per_m3_k)
        )
        q = _dg0_two_region(msh, cell_tags, q0, 0.0)
    t_inf = float(cooling.t_inf_c)
    t_init = t_inf if transient.t_init_c is None else float(transient.t_init_c)
    extra = {}
    if g is not None:  # `WF-1` step 28: P63's ``T_INF + G_GRAD * x[2]`` on every cell
        gg = float(g)
        extra["t_init_expr"] = lambda x: t_init + gg * x[2]
    if dose is not None:  # `WF-1` step 30: every full step, read below then cleared
        extra["snapshot_steps"] = tuple(range(int(transient.n_steps) + 1))
    s = solve_transient_conduction(
        msh,
        facet_tags,
        k=k,
        rho_c=rho_c,
        q=q,
        h=float(cooling.h_w_per_m2_k),
        t_inf=t_inf,
        t_init=t_init,
        robin_tag=_ROBIN_TAG,
        dt=float(transient.dt_s),
        n_steps=int(transient.n_steps),
        theta=float(transient.theta),
        startup_half_steps=int(transient.startup_half_steps),
        degree=int(cooling.degree),
        source_schedule=duty,
        **extra,
    )
    if dose is None:
        return TransientStudyResult(
            temperature=s.temperature,
            mesh=msh,
            cell_tags=cell_tags,
            facet_tags=facet_tags,
            k=k,
            rho_c=rho_c,
            q=q,
            source_multipliers=list(s.source_multipliers),
            source_power=s.source_power,
            max_balance_residual=s.max_balance_residual,
            solver_result=s,
            transient=transient,
            duty=duty,
            cooling=cooling,
        )
    return TransientStudyResult(
        temperature=s.temperature,
        mesh=msh,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        k=k,
        rho_c=rho_c,
        q=q,
        source_multipliers=list(s.source_multipliers),
        source_power=s.source_power,
        max_balance_residual=s.max_balance_residual,
        solver_result=s,
        transient=transient,
        duty=duty,
        cooling=cooling,
        dose=_dose_readout(s, transient, dose),
    )


def _dose_readout(s, transient, dose) -> dict:
    """`WF-1` step 30 (§9 item 878): each snapshot read on one fixed batch (the
    dose points, in order) by ``evaluate_vector_field_parallel`` (collective),
    then cleared — `MAT-21` P68's route
    (``tests/validation/test_mat21_wall_pre_equilibration_duty25_dose_fem.py:193–200``);
    per point `MAT-20`'s ``cem43`` and ``running_average`` on the full-step grid."""
    from .post.evaluation import evaluate_vector_field_parallel
    from .thermal.dose import cem43, running_average

    n_steps = int(transient.n_steps)
    pts = np.array([[float(c) for c in p] for p in dose.points_m], dtype=float)
    hist = np.empty((n_steps + 1, pts.shape[0]))
    for n in range(n_steps + 1):
        vals, ok = evaluate_vector_field_parallel(s.snapshots[n], pts)
        for j in range(pts.shape[0]):
            if not bool(ok[j]):
                s.snapshots.clear()
                raise ValueError(f"dose point {j} {dose.points_m[j]!r} was not found on the mesh")
            hist[n, j] = float(np.real(np.ravel(vals[j])[0]))
    s.snapshots.clear()
    times = float(transient.dt_s) * np.arange(n_steps + 1)
    window = float(dose.window_s)
    cem, peak, t_max = [], [], []
    for j in range(pts.shape[0]):
        cem.append(float(cem43(times, hist[:, j])))
        _means, mx, t_mx = running_average(times, hist[:, j], window=window)
        peak.append(float(mx))
        t_max.append(float(t_mx))
    return {
        "points_m": [[float(c) for c in p] for p in dose.points_m],
        "times_s": [float(t) for t in times],
        "temperature_c": [[float(v) for v in row] for row in hist],
        "cem43_s": cem,
        "six_minute": {"window_s": window, "peak_c": peak, "t_max_s": t_max},
    }


# ---------------------------------------------------------------------------
# `WF-15` step 4 (§9 item 326, chain U30): ``Study.probe`` — `MAT-19`'s
# ball-average virtual temperature probe through the API.  Additive: one call
# path, the public ``thermal.probe.ball_average``; no copy of the algorithm.
# ---------------------------------------------------------------------------


@dataclass
class ProbeStudyResult:
    """What ``Study.probe`` returns (not exported): one reading per probe."""

    probes: tuple
    values: np.ndarray
    volumes: np.ndarray
    nominal_volumes: np.ndarray
    max_imag_ratio: float
    theta: Any
    thermal: ThermalStudyResult


def _study_probe(self, probes, cooling: RobinCooling, *, levels: int = 3, transient=None,
                 duty=None) -> ProbeStudyResult:
    """Ball-averaged θ = T − ``t_inf_c`` over each ``ProbeBall`` (collective).

    Reuses ``self.temperature(cooling)`` (which solves once if ``self.result``
    is absent); the probe is ``thermal.probe.ball_average`` on the P1 submesh
    temperature, owned cells only (`MAT-19`'s ownership contract).
    """
    from dolfinx import fem

    from .thermal.probe import ProbeBall, ball_average

    probes = tuple(probes)
    if not probes:
        raise ValueError("probe: at least one ProbeBall is required")
    if isinstance(self.geometry, ThermalBoxGeometry) and self.geometry.wall is not None:
        raise ValueError("probe is not defined on a walled box (WF-1 step 28)")
    for p in probes:
        if not isinstance(p, ProbeBall):
            raise TypeError(f"probes must be ProbeBall instances, got {type(p).__name__}")
    if (transient is not None or duty is not None) and isinstance(
            self.geometry, SphereInBoxGeometry):
        raise ValueError("transient probes are not defined on the EM sphere "
                         "(WF-15 step 15 gates temperature only)")
    if transient is None and duty is None:
        th = self.temperature(cooling)
    else:
        # `WF-15` step 6 (§9 item 523): the thermal-only transient pass-through.
        th = self.temperature(cooling, transient=transient, duty=duty)
    T = th.temperature
    theta = fem.Function(T.function_space)
    theta.x.array[:] = T.x.array - cooling.t_inf_c
    theta.x.scatter_forward()
    values, volumes = [], []
    for p in probes:
        v, vol = ball_average(theta, p, levels=levels)
        values.append(v)
        volumes.append(vol)
    n_own = T.function_space.dofmap.index_map.size_local
    arr = np.asarray(theta.x.array[:n_own])
    comm = T.function_space.mesh.comm
    im = comm.allreduce(float(np.max(np.abs(np.imag(arr)), initial=0.0)), op=MPI.MAX)
    re = comm.allreduce(float(np.max(np.abs(np.real(arr)), initial=0.0)), op=MPI.MAX)
    return ProbeStudyResult(
        probes=probes,
        values=np.asarray(values, dtype=np.float64),
        volumes=np.asarray(volumes, dtype=np.float64),
        nominal_volumes=np.asarray([4.0 * np.pi * p.radius ** 3 / 3.0 for p in probes]),
        max_imag_ratio=im / re if re > 0.0 else 0.0,
        theta=theta,
        thermal=th,
    )


Study.probe = _study_probe


# ---------------------------------------------------------------------------
# `WF-13` step 2 (§9 items 268 / 287): ``Study.estimate()`` — a resource estimate
# before a solve.  Additive: builds the mesh and the N1curl space only (no
# form, no assembly, no solve) and prices the solve from measured rows.
# ---------------------------------------------------------------------------

#: Measured price rows (ours only; never an external number).  ``None`` = not
#: printed by the source log.  Memory is summed ``ru_maxrss`` over ranks on
#: every row (GiB); time is the recorded figure of each row's own definition
#: (heterogeneous on purpose, §9 item 287 — never redefined here).  One row
#: per fixture and mesh: duplicates of a fixture are circular in a
#: leave-one-out fit.  Sources (docs/testing/logs/ unless a ledger):
#:   P1 deg 1: 20261006T035826Z_WF-15-s2-W1-deg1.log (item 260 Outcome);
#:      N 90 209 from 20261006T185702Z_WF-13-s2-W1.log:204
#:   P2 20261005T170537Z_WF-15-r2-W2a-deg2-study.log:158 — Study solve + mesh,
#:      1 RHS, -n 2
#:   P3 deg 1: 20261003T093408Z_ANS-4-item192-W1-p1.log:1899 (time)
#:   P3 N and memory: 20261015T201542Z_WF-4-s1a-Wm.log:1874 (memory), :1876 (N)
#:      (Study route, the same fixture; time stays item 192's rung) — `WF-13`
#:      step 6, §9 item 471a.  Time and memory come from two windows of one
#:      fixture (the P1 precedent).  P3's peak is below P1's (another route's
#:      window), so the two-row degree-1 memory exponent is negative: the
#:      degree-1 rows price their own fixtures, not a trend.
#:   T49 deg 1 (the tune YAML, F-small exact-rectangle PEC hole, -n 2):
#:      time 20261014T110326Z_WF-1-s11-W1.log:1830 (item 442's
#:      run_study_file 29.668 s); N and memory
#:      20261016T035745Z_WF-13-s7-Wm.log:1788 (N), :1786 (memory) — `WF-13`
#:      step 7, §9 item 479.  Time and memory from two windows of one fixture
#:      (the P1 / P3 precedent).
#:   F16 20261005T203015Z_WF-16-W2-deg2.log:3850 (four driven solves, -n 8),
#:      :3871 (memory); N from 20261006T185702Z_WF-13-s2-W1.log:1989
#:   Cu7 20261004T202814Z_ANS-7-item228-W1-nodes.log:2159 — 231.3 s / 4 nodes
#:   Cu7h 20261004T111415Z_ANS-7-item214-W1-probe.log:2135 — one node, -n 12
#:   A6 20261001T070008Z_ANS-6-step2d.log:954 (N), :957 (PEC sweep, memory)
#:   S4c 20261006T070008Z_ANS-4-step4c.log:4525 (ZMUMPS N; the ``≈1,661,075``
#:      at :4863 is an estimate, not used), :4862 (memory), :4863 (four drives)
#:   P4 20261004T183338Z_ANS-4-item225-W3-p2-cf0p5-128-koff.log:2083 (cells,
#:      four drives), :2102 (memory); N 1 658 512 = size_global*bs read by
#:      20261007T033308Z_WF-13-s2-r-W1b-p4-unknowns.log:1918 (cells 259 524)
#:   P5 20261005T070009Z_ANS-4-step6b.log:2169 (ZMUMPS N), :2506 (summed
#:      ru_maxrss; the cgroup ``memory.peak`` 90.76 GiB is xl-ledger.md:51),
#:      :3707 (first sweep)
#:   P6 xxl-ledger.md:21 (20260926T070009Z_WF-7-step0d.log.gz:10427, :11903;
#:      xxl entry 9 — xl entry 29 carries the same price and has not run):
#:      factorising drive 1 of one factorisation + 31 back-substitutions.
#:   FH1 deg 1 (the copper Leontovich F-human hole, 369 387 cells, kappa off,
#:      tests/package/wf1/f_human_hole_copper_tune.yaml, -n 8):
#:      20261015T213534Z_WF-16-s12-W1.log:7343 (sweep 93.4 s — the EM solve,
#:      mesh excluded; 1 factorisation + 31 held), :7344 (memory); N from
#:      20261011T140313Z_WF-10-step4a-W1.log:6215 — `WF-13` step 10, §9 item
#:      562.  Time and memory from one window, N from another window of the
#:      same fixture (the P1 / P3 / T49 precedent).  T49's seconds include its
#:      mesh; FH1's do not (disclosed).
#:   FH2 deg 2 (the same copper F-human hole, 369 387 cells, kappa off,
#:      -n 16; `xl` entry 41): 20261022T070008Z_WF-10-step4-H28.log:6855
#:      (node pair 298.9 s — the EM solve, build excluded: assembly + 1
#:      factorisation + 1 held; summed ru_maxrss 88.58 GiB), N :6613
#:      (2E + 2F) — `WF-13` step 11, §9 item 585.  FH1 carries 31 held
#:      back-substitutions, FH2 one (disclosed).  The cgroup memory.peak
#:      87.12 GiB is :14681 (the P5 precedent).  Entries 35 / 42, H26 and
#:      entry 41's second-module pair are duplicates of this fixture: printed
#:      by the test, never fitted (item 212's precedent).
#: Item 212's two rows (P4's mesh at 10 / 64 MHz; 36.12 GiB / 186.93 s at
#: -n 8, 38.40 GiB / 183.12 s at -n 12) duplicate P4's fixture and are not
#: fitted (printed by the test as a rank-sensitivity reading).
PRICE_TABLE = (
    dict(row="P1", degree=1, ranks=2, cells=74020, unknowns=90209, seconds=13.7, peak_gib=2.39),
    dict(row="P3", degree=1, ranks=2, cells=117206, unknowns=140460, seconds=29.77, peak_gib=1.8185),
    dict(row="T49", degree=1, ranks=2, cells=80798, unknowns=111946, seconds=29.668, peak_gib=1.5296),
    dict(row="FH1", degree=1, ranks=8, cells=369387, unknowns=500568, seconds=93.4, peak_gib=9.7652),
    dict(row="P2", degree=2, ranks=2, cells=74020, unknowns=482158, seconds=282.4, peak_gib=20.41),
    dict(row="F16", degree=2, ranks=8, cells=117206, unknowns=753062, seconds=70.11, peak_gib=19.37),
    dict(row="Cu7", degree=2, ranks=8, cells=80798, unknowns=570332, seconds=57.83, peak_gib=14.72),
    dict(row="Cu7h", degree=2, ranks=12, cells=140270, unknowns=968016, seconds=210.1, peak_gib=32.83),
    dict(row="A6", degree=2, ranks=16, cells=226642, unknowns=1536270, seconds=283.2, peak_gib=82.34),
    dict(row="S4c", degree=2, ranks=16, cells=259543, unknowns=1658680, seconds=102.0, peak_gib=42.73),
    dict(row="P4", degree=2, ranks=8, cells=259524, unknowns=1658512, seconds=138.01, peak_gib=36.28),
    dict(row="P5", degree=2, ranks=16, cells=282198, unknowns=1822094, seconds=283.86, peak_gib=86.30),
    dict(row="P6", degree=2, ranks=16, cells=507266, unknowns=3256418, seconds=322.95, peak_gib=104.781),
    dict(row="FH2", degree=2, ranks=16, cells=369387, unknowns=2570024, seconds=298.9, peak_gib=88.58),
)

#: `WF-13` step 12 (§9 item 663): the F-clinical-1.5T r1 price, **outside**
#: ``PRICE_TABLE`` (a row there would move every per-degree fit).  Measured,
#: `xl` entries 44 / 46 at -n 16: ``20261026T070008Z_GEO-38-C4.log:12142``
#: (cells 1 147 366 = record, t_build 357.7 s), ``:14156`` (build 357.7 s,
#: Z50 pair 203.3 s -- the EM solve, build excluded, the FH1 / FH2 convention),
#: ``:14211`` (cgroup memory.peak 52 689 350 656 B = 49.07 GiB; the summed RSS
#: 64.89 GiB double-counts shared pages, xl-ledger.md:69).  Read directly, never
#: fitted; the last row of a degree wins.
#:
#: `WF-13` step 13 (§9 item 715): FC2, degree 2, measured, one `xxl` window
#: (entry 50, -n 16): ``20261031T070008Z_GEO-38-C5a.log:12595`` (build 341.8 s,
#: Z50 pair 2892.1 s, FC1's convention, build excluded), ``:12650`` (cgroup
#: memory.peak 507 256 434 688 B = 472.42 GiB; the summed RSS 525.24 GiB
#: double-counts shared pages).  Factorisation driver 2 525 -- 2 700 s, measured
#: (``:12285``, ``:12506``).  MUMPS INFOG(19) varies 7.5 % node to node, so this
#: is a single window, not a band.  Tier: `xl` at ≤ 3 root-grid nodes, `xxl`
#: for a column set.  No 64 x 64 price, thermal price or 3 T price is implied.
PRICE_TABLE_CLINICAL = (
    dict(row="FC1", degree=1, ranks=16, cells=1147366, unknowns=1525193, seconds=203.3,
         build_seconds=357.7, peak_gib=49.07),
    dict(row="FC2", degree=2, ranks=16, cells=1147366, unknowns=7885674, seconds=2892.1,
         build_seconds=341.8, peak_gib=472.42),
)

#: `WF-13` step 15 (§9 item 782, §10 U118): the F-clinical-3T r1 twin's price,
#: in its **own** table, keyed by ``field_strength_t == 3.0``, so the "last row
#: of a degree wins" lookup cannot let a 3 T row shadow FC1 (and the arity pins
#: on ``PRICE_TABLE_CLINICAL`` in steps 12 / 13's modules stay unedited).
#: FC3, degree 1, measured, `xl` entry 49 (-n 16):
#: ``20261101T080008Z_GEO-38-C16.log:14125`` (build 327.6 s, **Z50 pair 135.1 s**
#: -- FC1's convention, one EM solve, build excluded; the node pairs 146.3 --
#: 173.0 s are not the price), ``:14183`` (cgroup memory.peak 42 297 978 880 B =
#: 39.39 GiB; the summed RSS 49.10 GiB double-counts shared pages).  Cells and
#: unknowns are `GEO-38` C8's 3 T r1 records.  FC4 (3 T, degree 2) has **no**
#: row: ``not_priced``, with U_2 = 7 335 988 from the record.  Read directly,
#: never fitted.
PRICE_TABLE_CLINICAL_3T = (
    dict(row="FC3", degree=1, ranks=16, cells=1057895, unknowns=1425765, seconds=135.1,
         build_seconds=327.6, peak_gib=39.39),
)

#: `WF-13` step 16 (§9 item 783, §10 U119): the rod-in-clinical-box price, in
#: its **own** table (the arity pins on ``PRICE_TABLE_CLINICAL`` stay unedited,
#: as item 782's ``PRICE_TABLE_CLINICAL_3T``).  FC-ROD1, 1.5 T, degree 1,
#: measured, `xl` entry 48 (-n 8): ``20261102T080008Z_GEO-36-step9-I36.log:12503``
#: (build 461.7 s, **arm M 263.5 s** -- FC1's convention, one EM solve, build
#: excluded; arm M is the rod-bearing solve, arm R 262.3 s the reference),
#: ``:12626`` / ``:12627`` (cgroup memory.peak 45 519 204 352 B = 42.39 GiB; the
#: summed RSS 46.048 GiB is not the price).  Cells / unknowns are I33's records.
PRICE_TABLE_CLINICAL_ROD = (
    dict(row="FC-ROD1", degree=1, ranks=8, cells=1217025, unknowns=1606744, seconds=263.5,
         build_seconds=461.7, peak_gib=42.39),
)


def _n1curl_unknowns(k: int, edges: int, faces: int, cells: int) -> int:
    """First-kind N1curl degree ``k`` on tetrahedra (`WF-13` step 14, §9 item 731).

    ``U_k = k E + k(k - 1) F + k(k - 1)(k - 2) T / 2`` (Nédélec 1980; local
    dimension ``k(k + 2)(k + 3)/2``).
    """
    return k * edges + k * (k - 1) * faces + k * (k - 1) * (k - 2) * cells // 2


def _f_clinical_unknowns(degree: int, record=None, cells=None) -> int:
    """``U_k`` from the record: ``U_1 = E``, ``U_2 = 2E + 2F``; no mesh.

    ``record`` / ``cells`` default to the 1.5 T r1 fixture; the 3 T twin passes
    its own (`WF-13` step 15, §9 item 782).
    """
    k = int(degree)
    if k < 1:
        raise ValueError(f"degree must be >= 1, got {degree!r}")
    u1, u2 = F_CLINICAL_R1_UNKNOWNS_RECORD if record is None else record
    cells = F_CLINICAL_R1_CELL_RECORD if cells is None else cells
    return _n1curl_unknowns(k, u1, (u2 - 2 * u1) // 2, cells)


def _estimate_f_clinical(study):
    """``Study.estimate`` on an ``FClinicalBirdcageGeometry``: no mesh, no fit."""
    if study.geometry == FClinicalBirdcageGeometry(implant_rod=F_CLINICAL_ROD_ENTRY48):
        # `WF-13` step 16 (§9 item 783): entry 48's rod, its own records and table.
        cells = F_CLINICAL_ROD_CELL_RECORD
        unknowns = _f_clinical_unknowns(study.degree, F_CLINICAL_ROD_UNKNOWNS_RECORD, cells)
        table = PRICE_TABLE_CLINICAL_ROD
    elif study.geometry == FClinicalBirdcageGeometry(field_strength_t=3.0):
        # `WF-13` step 15 (§9 item 782): the 3 T twin, its own records and table.
        cells = F_CLINICAL_3T_R1_CELL_RECORD
        unknowns = _f_clinical_unknowns(study.degree, F_CLINICAL_3T_R1_UNKNOWNS_RECORD, cells)
        table = PRICE_TABLE_CLINICAL_3T
    elif study.geometry == FClinicalBirdcageGeometry():
        cells = F_CLINICAL_R1_CELL_RECORD
        unknowns = _f_clinical_unknowns(study.degree)
        table = PRICE_TABLE_CLINICAL
    else:
        raise ValueError("only the F-clinical-1.5T r1 fixture is priced "
                         "(FClinicalBirdcageGeometry() defaults; WF-13 step 12), "
                         "and its 3 T twin (field_strength_t=3.0; WF-13 step 15)")
    rows = [r for r in table if r["degree"] == int(study.degree)]
    if not rows:
        return ResourceEstimate(
            cells=cells, unknowns=unknowns, degree=study.degree, ranks=16,
            predicted_peak_gib=float("nan"), predicted_seconds=float("nan"),
            rows_used=(),
            not_priced=(f"degree {study.degree} (no measured F-clinical row)",
                        "thermal transient"),
        )
    row = rows[-1]
    return ResourceEstimate(
        cells=cells, unknowns=unknowns, degree=study.degree, ranks=int(row["ranks"]),
        predicted_peak_gib=float(row["peak_gib"]), predicted_seconds=float(row["seconds"]),
        rows_used=(row["row"],), not_priced=("thermal transient",),
    )


@dataclass(frozen=True)
class ResourceEstimate:
    """What ``Study.estimate`` returns (not exported, as ``PortStudyResult``)."""

    cells: int
    unknowns: int
    degree: int
    ranks: int
    predicted_peak_gib: float
    predicted_seconds: float
    rows_used: tuple
    not_priced: tuple = ()


def _fit_power_law(xs, ys):
    """Least-squares ``log y = log a + b log x``; returns ``(a, b)``."""
    lx = np.log(np.asarray(xs, dtype=float))
    ly = np.log(np.asarray(ys, dtype=float))
    b, loga = np.polyfit(lx, ly, 1)
    return float(np.exp(loga)), float(b)


def _predict(unknowns: int, degree: int, ranks: int, n_rhs: int = 1, table=PRICE_TABLE):
    """Price a solve from ``table``: ``mem = a·N^b`` and ``time·ranks = c·N^e``.

    Fitted per degree on the rows that carry unknowns and the quantity.
    ``n_rhs`` back-substitutions are not priced in v0 (one factorisation
    dominates).  Returns ``(peak_gib, seconds, rows_used)``.
    """
    rows = [r for r in table if r["degree"] == int(degree) and r["unknowns"]]
    mem_rows = [r for r in rows if r["peak_gib"]]
    t_rows = [r for r in rows if r["seconds"]]
    if len(mem_rows) < 2 or len(t_rows) < 2:
        raise ValueError(f"fewer than two priced rows at degree {degree}")
    a, b = _fit_power_law([r["unknowns"] for r in mem_rows], [r["peak_gib"] for r in mem_rows])
    c, e = _fit_power_law(
        [r["unknowns"] for r in t_rows], [r["seconds"] * r["ranks"] for r in t_rows]
    )
    peak = a * float(unknowns) ** b
    seconds = c * float(unknowns) ** e / float(ranks)
    return float(peak), float(seconds), tuple(r["row"] for r in rows)


def _study_estimate(self, comm: MPI.Intracomm = MPI.COMM_WORLD, *, n_rhs: Optional[int] = None):
    """Cells, exact N1curl unknowns at ``self.degree``, and a priced peak/time.

    Builds the mesh and the function space only; ``self.result`` is untouched.
    """
    if isinstance(self.geometry, FClinicalBirdcageGeometry):  # `WF-13` step 12 (§9 item 663)
        return _estimate_f_clinical(self)
    from dolfinx import fem

    _refuse_thermal_only(self, "estimate")
    not_priced = ()
    if isinstance(self.geometry, FHumanBirdcageGeometry):
        msh = self.geometry.build(comm)["mesh"]
        V = fem.functionspace(msh, ("N1curl", self.degree))
        n_rhs = 1 if n_rhs is None else n_rhs
        not_priced = ("thermal transient",)
    elif isinstance(self.geometry, FourPortBirdcageGeometry):
        msh = self.geometry.build(comm)["mesh"]
        V = fem.functionspace(msh, ("N1curl", self.degree))
        n_rhs = 4 if n_rhs is None else n_rhs
    else:
        msh, cell_tags, _ = self.geometry.build(comm)
        problem = TimeHarmonicProblem(
            mesh=msh,
            frequency_hz=self.frequency_hz,
            material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0),
            cell_tags=cell_tags,
            phantom_material=self.material,
            phantom_tag=self.geometry.sphere_tag,
            boundary_condition="pec_zero_tangential_a",
            dirichlet_e_field=self.dirichlet_field(),
        )
        V = TimeHarmonicSolver(problem, degree=self.degree).function_space()
        n_rhs = 1 if n_rhs is None else n_rhs
    cells = int(msh.topology.index_map(msh.topology.dim).size_global)
    unknowns = int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)
    ranks = int(comm.size)
    peak, seconds, used = _predict(unknowns, self.degree, ranks, n_rhs)
    return ResourceEstimate(
        cells=cells,
        unknowns=unknowns,
        degree=self.degree,
        ranks=ranks,
        predicted_peak_gib=peak,
        predicted_seconds=seconds,
        rows_used=used,
        not_priced=not_priced,
    )


Study.estimate = _study_estimate


# ---------------------------------------------------------------------------
# `WF-16` step 4 (§9 item 321, chain U23): ``Study.sar`` — the C95.3
# mass-averaging operator (``post.sar.mass_averaged_sar``) on the port Study's
# kept drive fields.  Insertions only: ``_study_solve_ports``, ``solve``,
# ``solve_sweep`` and ``mean_sar`` are untouched; the field-keeping solve
# repeats ``_study_solve_ports``'s problem lines (copied, not factored).
# ---------------------------------------------------------------------------

#: `MAT-4` step 3's pinned quadrature degree (``QUADRATURE_DEGREE`` in
#: tests/validation/test_mass_averaged_sar_standard_masses.py, re-exported to
#: test_birdcage_sar_mass_averaged.py); the step-4 test asserts equality by
#: import.  Never unpinned: a ``SpatialCoordinate`` integral at the automatic
#: degree can compile for > 9 min (POST-5).
SAR_QUADRATURE_DEGREE = 16


def _study_solve_port_fields(
    self, comm: MPI.Intracomm = MPI.COMM_WORLD
) -> PortStudyResult:
    """``_study_solve_ports`` with ``keep_fields=True``; cached on ``self.fields_result``.

    ``fields`` of the returned ``PortStudyResult`` is the sweep's
    ``{port_id: TimeHarmonicFields}``.  ``self.result`` is not touched.
    """
    import time

    from .ports.sparameters import run_n_port_sparameter_sweep

    # `WF-15` step 12 (§9 item 729), ahead of the cache: sar and b1_plus reach it here.
    _refuse_rod_route(self, "the field-keeping port solve")
    cached = getattr(self, "fields_result", None)
    if cached is not None:
        return cached
    if not isinstance(self.geometry, FourPortBirdcageGeometry):
        raise TypeError("the field-keeping port solve is defined on a port Study only")
    if not isinstance(self.excitation, LumpedPortDrive):
        raise TypeError("a FourPortBirdcageGeometry Study needs a LumpedPortDrive")
    geom = self.geometry
    built = geom.build(comm)
    msh = built["mesh"]
    # `WF-16` step 11 (§9 item 388): the conductor set-up, the hole-route
    # `material_map`, `**cond_kwargs`, the conditional `extra_bilinear_terms`
    # and the two `conductor_*` result fields are copied, in order, from
    # `_study_solve_ports` (main at e122cead) -- copied, not factored out.
    cond_kwargs, cond_terms = _conductor_setup(self, built, comm)
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=self.frequency_hz,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=built["cell_tags"],
        # `WF-16` step 7: the hole route has no conductor cells to map.
        material_map=({} if geom.as_hole else {
            geom.conductor_tag: HomogeneousMaterial(
                sigma=geom.conductor_sigma_s_per_m, epsilon_r=1.0, mu_r=1.0
            ),
        }),
        phantom_material=self.material,
        phantom_tag=geom.phantom_tag,
        boundary_condition="pec_zero_tangential_a",
        **cond_kwargs,
    )
    specs, port_defs = self.excitation.specs_and_definitions(geom, built["sheets"])
    comm.Barrier()
    t0 = time.perf_counter()
    sweep = run_n_port_sparameter_sweep(
        problem,
        port_defs,
        lumped_sheet_ports=specs,
        lumped_sheet_facet_tags=built["facet_tags"],
        degree=self.degree,
        keep_fields=True,
        **({"extra_bilinear_terms": cond_terms} if cond_terms else {}),
    )
    comm.Barrier()
    t_sweep = time.perf_counter() - t0
    self.fields_result = PortStudyResult(
        mesh=msh,
        cell_tags=built["cell_tags"],
        facet_tags=built["facet_tags"],
        fields=dict(sweep.fields),
        solver=None,
        degree=self.degree,
        comm=comm,
        problem=problem,
        sweep=sweep,
        s_matrix=np.asarray(sweep.s_matrix, dtype=np.complex128),
        z_matrix=np.asarray(sweep.z_matrix, dtype=np.complex128),
        port_definitions=port_defs,
        specs=specs,
        sheets=built["sheets"],
        cells=built["cells"],
        mesh_time=built["mesh_time"],
        sweep_time=float(t_sweep),
        conductor_terms=cond_terms,
        conductor_facet_tags=built["mesh_facet_tags"] if cond_terms else None,
    )
    return self.fields_result


Study._solve_port_fields = _study_solve_port_fields


def _global_tags(cell_tags, comm) -> list:
    """The MPI-reduced set of cell-tag values (``cell_tags.values`` is rank-local)."""
    local = {int(v) for v in np.unique(np.asarray(cell_tags.values))}
    merged = set()
    for part in comm.allgather(sorted(local)):
        merged.update(part)
    return sorted(merged)


def _study_density_field(self, res):
    """ρ = ``density_kg_m3`` on the phantom tag, 0.0 on every other tag."""
    from .post.sar import build_density_field

    tags = _global_tags(res.cell_tags, res.comm)
    return build_density_field(
        res.mesh,
        self.density_kg_m3,
        cell_tags=res.cell_tags,
        density_map={t: 0.0 for t in tags if t != int(self.geometry.phantom_tag)},
    )


Study._density_field = _study_density_field


class _MassDefault:
    """``mass_g``'s default: 10.0 g unless ``radius_m`` is given."""

    def __repr__(self) -> str:
        return "10.0"


_MASS_DEFAULT = _MassDefault()


def _study_sar(
    self,
    mass_g: Any = _MASS_DEFAULT,
    *,
    centres,
    drives=None,
    radius_m: Optional[float] = None,
    quadrature_degree: int = SAR_QUADRATURE_DEGREE,
    comm: MPI.Intracomm = MPI.COMM_WORLD,
) -> dict:
    """Mass-averaged SAR (``post.sar.mass_averaged_sar``) per drive per centre.

    Port Study only.  Exactly one of ``mass_g`` / ``radius_m``: ``mass_g``
    defaults to 10.0 g when ``radius_m`` is not given; an explicit ``mass_g``
    together with ``radius_m``, or ``mass_g=None`` without it, is refused.  The ball radius for ``mass_g`` is
    ``averaging_ball_radius(mass_g/1000, density_kg_m3)``.  Returns
    ``{port_id: [operator dict per centre]}``, each dict unchanged.  A symmetry
    / operator record per 1 V drive — no compliance claim.
    """
    from .post.sar import averaging_ball_radius, mass_averaged_sar

    _refuse_loaded(self, "sar")  # WF-16 step 14 (§9 item 661)
    if not isinstance(self.geometry, FourPortBirdcageGeometry):
        hook = getattr(self.excitation, "_sar_port_study", None)
        if hook is not None:
            return hook(self, mass_g, centres=centres, drives=drives, radius_m=radius_m,
                        quadrature_degree=quadrature_degree, comm=comm)
        raise TypeError("sar is defined on a port Study only (WF-16); use mean_sar")
    if mass_g is _MASS_DEFAULT:
        mass_g = 10.0 if radius_m is None else None
    if (mass_g is None) == (radius_m is None):
        raise ValueError("give exactly one of mass_g and radius_m")
    if radius_m is None:
        if not np.isfinite(float(mass_g)) or float(mass_g) <= 0.0:
            raise ValueError(f"mass_g must be finite and positive, got {mass_g!r}")
        radius = averaging_ball_radius(
            mass_kg=float(mass_g) / 1000.0, rho=self.density_kg_m3
        )
    else:
        radius = float(radius_m)
        if not np.isfinite(radius) or radius <= 0.0:
            raise ValueError(f"radius_m must be finite and positive, got {radius_m!r}")
    res = self._solve_port_fields(comm)
    rho = getattr(self, "_rho_field", None)
    if rho is None:
        rho = self._density_field(res)
        self._rho_field = rho
    ids = list(res.fields) if drives is None else list(drives)
    out = {}
    for pid in ids:
        f = res.fields[pid]
        out[pid] = [
            mass_averaged_sar(
                f.e_complex,
                sigma=f.sigma_field,
                rho=rho,
                center=tuple(float(x) for x in c),
                radius=radius,
                comm=res.comm,
                quadrature_degree=int(quadrature_degree),
            )
            for c in centres
        ]
    return out


Study.sar = _study_sar


def _study_b1_plus(
    self, *, drives=None, comm: MPI.Intracomm = MPI.COMM_WORLD
) -> dict:
    """The gated CG1 ``B`` estimator per drive (`WF-6` step 1d), on the kept field.

    Port Study only.  Reuses ``self.fields_result`` (``_solve_port_fields``) and
    returns ``{port_id: CG1 vector Function}`` =
    ``project_to_cg1(magnetic_flux_density_from_e(e, ω))`` — never the DG0 curl.
    ``|B₁⁺| = |B_x + jB_y|/2`` is formed *after* point evaluation
    (:meth:`Study.b1_plus_at`), the order `WF-6` gated: ``|·|`` is not linear,
    and ``post.faraday.b1_plus`` is a DG0 operator, so it is not applied to the
    projection.  Cached per drive on ``self._b1_cg1``.
    """
    from .post.faraday import magnetic_flux_density_from_e, project_to_cg1

    _refuse_loaded(self, "b1_plus")  # WF-16 step 14 (§9 item 661); b1_plus_at calls this
    if not isinstance(self.geometry, FourPortBirdcageGeometry):
        hook = getattr(self.excitation, "_b1_plus_port_study", None)
        if hook is not None:
            return hook(self, drives=drives, comm=comm)
        raise TypeError("b1_plus is defined on a port Study only (WF-16)")
    res = self._solve_port_fields(comm)
    cache = getattr(self, "_b1_cg1", None)
    if cache is None:
        cache = {}
        self._b1_cg1 = cache
    ids = list(res.fields) if drives is None else list(drives)
    for pid in ids:
        if pid not in cache:
            cache[pid] = project_to_cg1(
                magnetic_flux_density_from_e(res.fields[pid].e_complex, self.omega)
            )
    return {pid: cache[pid] for pid in ids}


Study.b1_plus = _study_b1_plus


def _study_b1_plus_at(
    self, points, *, drives=None, comm: MPI.Intracomm = MPI.COMM_WORLD
) -> dict:
    """``|B₁⁺|`` at ``points`` per drive: ``{port_id: (values, valid)}``.

    Evaluates :meth:`Study.b1_plus`'s CG1 vector through
    ``post.evaluation.evaluate_vector_field_parallel`` (never ``f.eval``) and
    forms ``|B_x + jB_y|/2`` from the evaluated complex vector.  ``points`` must
    be identical on every rank.  Per 1 V drive; a record, not an absolute field.
    """
    from .post.evaluation import evaluate_vector_field_parallel

    maps = self.b1_plus(drives=drives, comm=comm)
    pts = np.ascontiguousarray(np.asarray(points, dtype=np.float64).reshape(-1, 3))
    out = {}
    for pid, fn in maps.items():
        values, valid = evaluate_vector_field_parallel(fn, pts)
        values = np.asarray(values).reshape(-1, 3)
        out[pid] = (
            np.abs(values[:, 0] + 1j * values[:, 1]) / 2.0,
            np.asarray(valid, dtype=bool),
        )
    return out


Study.b1_plus_at = _study_b1_plus_at


# ---------------------------------------------------------------------------
# `WF-16` step 7 (§9 item 324, chain U26): the port Study on copper.
# Insertions only.  A ``LeontovichWall`` on the hole route's conductor wall
# (tag ``BIRDCAGE_CONDUCTOR_SURFACE_TAG``, 401), the outer box pinned as its own
# group (``OUTER_BOX_TAG``, 499) -- `TH-14` step 2's fixture through the API.
# ``conductor=None`` adds no keyword and no term: the default is unchanged.
# ---------------------------------------------------------------------------

#: The gated port-strip fraction (`tests/validation/test_port_lumped_narrowed_sheet`).
GATED_WIDTH_FRACTION = 0.5
#: Copper conductivity, S/m (`tests/validation/test_th14_birdcage_copper`).
COPPER_SIGMA = 5.8e7
#: The outer-box facet group the copper Study pins (exterior minus tag 401).
OUTER_BOX_TAG = 499


@dataclass(frozen=True)
class LeontovichWall:
    """Surface-impedance conductor: ``Z_s = (1 + j) R_s(omega, sigma)`` on tag 401 (Jin 5.8.3)."""

    sigma_s_per_m: float = COPPER_SIGMA

    def surface_impedance_ohm(self, omega_rad_s: float) -> complex:
        from .core.cavity import surface_resistance_ohm

        return (1.0 + 1.0j) * surface_resistance_ohm(omega_rad_s, self.sigma_s_per_m)


def _outer_box_tags(msh, facet_tags):
    import dolfinx

    from .io.mesh import BIRDCAGE_CONDUCTOR_SURFACE_TAG

    fdim = msh.topology.dim - 1
    msh.topology.create_connectivity(fdim, msh.topology.dim)
    exterior = np.asarray(dolfinx.mesh.exterior_facet_indices(msh.topology), dtype=np.int32)
    cavity = np.asarray(facet_tags.find(BIRDCAGE_CONDUCTOR_SURFACE_TAG), dtype=np.int32)
    outer = np.setdiff1d(exterior, cavity).astype(np.int32)
    return dolfinx.mesh.meshtags(
        msh, fdim, outer, np.full(outer.size, OUTER_BOX_TAG, dtype=np.int32)
    )


def _conductor_setup(self, built, comm):
    """``(problem keywords, extra bilinear terms)`` for ``self.conductor``."""
    cond = getattr(self, "conductor", None)
    if cond is None:
        return {}, None
    if not isinstance(cond, LeontovichWall):
        raise TypeError(
            f"conductor must be None or a LeontovichWall, got {type(cond).__name__}"
        )
    if not getattr(self.geometry, "as_hole", False):
        raise ValueError(
            "a LeontovichWall conductor needs FourPortBirdcageGeometry(as_hole=True)"
        )
    from .core.boundary_terms import impedance_boundary_term
    from .io.mesh import BIRDCAGE_CONDUCTOR_SURFACE_TAG

    msh = built["mesh"]
    outer = _outer_box_tags(msh, built["mesh_facet_tags"])
    term = impedance_boundary_term(
        msh,
        built["mesh_facet_tags"],
        BIRDCAGE_CONDUCTOR_SURFACE_TAG,
        cond.surface_impedance_ohm(self.omega),
        self.omega,
    )
    return {"facet_tags": outer, "pec_facet_tags": (OUTER_BOX_TAG,)}, [term]


# ---------------------------------------------------------------------------
# `WF-16` step 9 part A (§9 item 382, chain U36): the tuned copper F-human
# through the ``Study``.  Insertions only (beside the ``as_hole`` field on
# ``FHumanBirdcageGeometry``); reached through item 373's
# ``_solve_port_study`` hook inside ``_study_solve_ports``' refusal branch.
# ---------------------------------------------------------------------------


# Copied, with provenance, from `main` at this item's base (e163ca93):
# tests/validation/test_wf10_f_human_copper_esr_root.py lines 76–103
# (`_thevenin_feed_impedance`, `_source_voltage`, `_current_route_s`) and
# lines 121–131 (`_y_from_gamma`, `_active_y_in`).
def _thevenin_feed_impedance(z_c, z0):
    """Z_th of a z0 source in parallel with the gap's Z_C."""
    return z0 * z_c / (z0 + z_c)


def _source_voltage(z_c, z0):
    """The z0 source's open-circuit V_s for V_th = 1 V across Z_C."""
    return (z0 + z_c) / z_c


def _current_route_s(currents_by_drive, feeds_idx, z_th, v_s):
    """``S`` (len(feeds)²) from port currents: column j from the drive of feed j."""
    n = len(feeds_idx)
    s = np.zeros((n, n), dtype=np.complex128)
    for j in range(n):
        cur = currents_by_drive[j]
        for i, idx in enumerate(feeds_idx):
            d = 1.0 if i == j else 0.0
            v = d - z_th * cur[idx]
            s[i, j] = 2.0 * v / v_s - d
    return s


def _y_from_gamma(gamma, z0):
    return (1.0 - gamma) / (z0 * (1.0 + gamma))


def _active_y_in(s2, a, z0):
    """``[Y_in at feed 0, Y_in at feed 1]`` under incident waves ``a`` (2-vector)."""
    a = np.asarray(a, dtype=np.complex128)
    b = np.asarray(s2, dtype=np.complex128) @ a
    return [complex(_y_from_gamma(b[i] / a[i], z0)) for i in range(2)]


def _sheet_azimuth_deg(msh, facet_tags, tag, comm) -> float:
    """Azimuth of a sheet's bounding-box centre [deg], reduced over ranks.

    Copied from `tests/mesh/test_birdcage_leg_offset.py` lines 174–225
    (`_sheet_nodes` + `_sheet_azimuth_deg`) at e163ca93.
    """
    import dolfinx

    fdim = msh.topology.dim - 1
    facets = np.asarray(facet_tags.find(tag), dtype=np.int32)
    if facets.size:
        nodes = dolfinx.cpp.mesh.entities_to_geometry(msh._cpp_object, fdim, facets, False)
        pts = msh.geometry.x[nodes.reshape(-1)]
        lo, hi = pts.min(axis=0), pts.max(axis=0)
    else:
        lo, hi = np.full(3, np.inf), np.full(3, -np.inf)
    lo = np.min(np.array(comm.allgather(lo)), axis=0)
    hi = np.max(np.array(comm.allgather(hi)), axis=0)
    centre = 0.5 * (lo + hi)
    return float(np.degrees(np.arctan2(centre[1], centre[0])) % 360.0)


@dataclass(frozen=True)
class TwoPortQuadratureDrive:
    """The two-port 90° feed of the tuned F-human with series-RC gaps.

    The non-feed ring gaps carry ``Z_RC = series_rlc_impedance(f, r_ohm=esr_ohm,
    l_h=0, c_f=c_in_f)``; each feed is a Thevenin source across its own gap
    (``Z_th = z0 Z_RC/(z0 + Z_RC)``, ``V_s = (z0 + Z_RC)/Z_RC``).  ``Study.solve``
    returns the feeds' node S, the active input admittance per sense and the two
    superposed sense fields (``ccw`` lags with azimuth).  ``c_in_f`` is an input,
    not a tuning result.
    """

    c_in_f: float
    esr_ohm: float = 0.0
    feeds: tuple = ("P17", "P21")
    z0_ohm: float = 50.0

    def __post_init__(self):
        c = float(self.c_in_f)
        if not np.isfinite(c) or c <= 0.0:
            raise ValueError("c_in_f must be finite and positive")
        r = float(self.esr_ohm)
        if not np.isfinite(r) or r < 0.0:
            raise ValueError("esr_ohm must be finite and non-negative")
        z0 = float(self.z0_ohm)
        if not np.isfinite(z0) or z0 <= 0.0:
            raise ValueError("z0_ohm must be finite and positive")
        if isinstance(self.feeds, str):
            raise ValueError("feeds is a pair of port ids, e.g. ('P17', 'P21')")
        feeds = tuple(self.feeds)
        if len(feeds) != 2 or feeds[0] == feeds[1]:
            raise ValueError("feeds must name two distinct port sheets")
        object.__setattr__(self, "feeds", feeds)

    def gap_impedance_ohm(self, frequency_hz: float) -> complex:
        from .ports.lumped import series_rlc_impedance

        return series_rlc_impedance(frequency_hz, r_ohm=self.esr_ohm, l_h=0.0,
                                    c_f=self.c_in_f)

    def feed_impedance_ohm(self, frequency_hz: float) -> complex:
        return _thevenin_feed_impedance(self.gap_impedance_ohm(frequency_hz), self.z0_ohm)

    def source_voltage_v(self, frequency_hz: float) -> complex:
        return _source_voltage(self.gap_impedance_ohm(frequency_hz), self.z0_ohm)

    def _solve_port_study(self, study, comm):
        return _study_solve_tuned_quadrature(study, comm)

    def _sar_port_study(self, study, mass_g, **kwargs):
        return _tuned_quadrature_sar(study, mass_g, **kwargs)

    def _b1_plus_port_study(self, study, **kwargs):
        return _tuned_quadrature_b1_plus(study, **kwargs)

    def _temperature_port_study(self, study, cooling, *, transient, balls, window_s,
                                sense="ccw", comm=MPI.COMM_WORLD):
        return _tuned_quadrature_temperature(
            study, cooling, transient=transient, balls=balls, window_s=window_s,
            sense=sense, comm=comm,
        )


@dataclass
class _SenseField:
    """One superposed sense field, in the shape ``sar`` / ``b1_plus`` read."""

    e_complex: Any
    sigma_field: Any


@dataclass
class TunedQuadratureResult:
    """What ``Study.solve`` returns for a ``TwoPortQuadratureDrive`` (not exported)."""

    mesh: Any
    cell_tags: Any
    facet_tags: Any
    problem: Any
    specs: Any
    port_definitions: Any
    feeds: tuple
    az_feeds: tuple
    az_mirror: float
    weights: dict
    s_node: Any
    y_in: dict
    currents: Any
    fields: dict
    cells: int
    times: dict
    conductor_terms: Any
    conductor_facet_tags: Any
    comm: Any
    degree: int
    gap_impedance_ohm: complex = 0j
    feed_impedance_ohm: complex = 0j
    source_voltage_v: complex = 0j
    sheets: Any = None


def _study_solve_tuned_quadrature(self, comm: MPI.Intracomm = MPI.COMM_WORLD):
    """The node pair and the two sense superpositions.

    The node lines are lifted, with provenance, from `main` at this item's base
    (e163ca93): tests/validation/test_wf10_f_human_tuned_b1_sar.py lines
    364–413 (`_drive`), 438–453 (`_pair`), 489–495 (the node pair) and
    525–536 (the superposition).  The fixture is the Study's own build.
    """
    import dataclasses
    import time

    from dolfinx import fem

    from .ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )

    exc = self.excitation
    geom = self.geometry
    if not isinstance(geom, FHumanBirdcageGeometry):
        raise TypeError("TwoPortQuadratureDrive needs an FHumanBirdcageGeometry")
    if not geom.as_hole:
        raise ValueError("TwoPortQuadratureDrive needs FHumanBirdcageGeometry(as_hole=True)")
    if not isinstance(self.conductor, LeontovichWall):
        raise TypeError("TwoPortQuadratureDrive needs conductor=LeontovichWall()")
    valid = {f"P{i}" for i in geom.ring_ordinals()}
    bad = [p for p in exc.feeds if p not in valid]
    if bad:
        raise ValueError(f"feeds name {bad}, not ring ports of this geometry")

    times = {}
    t0 = time.perf_counter()
    built = geom.build(comm)
    times["build"] = time.perf_counter() - t0
    msh = built["mesh"]
    tags_f = built["facet_tags"]
    cond_kwargs, cond_terms = _conductor_setup(self, built, comm)
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=self.frequency_hz,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=built["cell_tags"],
        material_map={},
        phantom_material=self.material,
        phantom_tag=geom.phantom_tag,
        boundary_condition="pec_zero_tangential_a",
        **cond_kwargs,
    )
    omega = self.omega
    z0 = float(exc.z0_ohm)
    specs, port_defs = LumpedPortDrive(z0_ohm=z0).specs_and_definitions(geom, built["sheets"])
    z_rc = exc.gap_impedance_ohm(self.frequency_hz)
    z_th = _thevenin_feed_impedance(z_rc, z0)
    v_s = _source_voltage(z_rc, z0)
    feeds = tuple(exc.feeds)
    specs_x = [dataclasses.replace(sp, port_impedance_ohm=(z_th if sp.port_id in feeds else z_rc))
               for sp in specs]
    order = list(feeds) + sorted((sp.port_id for sp in specs_x if sp.port_id not in feeds),
                                 key=lambda p: int(p[1:]))
    feeds_idx = [0, 1]

    az_feeds = tuple(_sheet_azimuth_deg(msh, tags_f, geom.sheet_iface + int(p[1:]), comm)
                     for p in feeds)
    d_az = (az_feeds[1] - az_feeds[0]) % 360.0
    az_mirror = az_feeds[0] + 0.5 * d_az
    d_rad = np.radians(d_az)
    weights = {"ccw": (1.0 + 0j, complex(np.exp(-1j * d_rad))),
               "cw": (1.0 + 0j, complex(np.exp(+1j * d_rad)))}

    solver = TimeHarmonicSolver(problem, degree=self.degree)
    wall = list(cond_terms)
    keep = {}
    cur = []
    t_pair = []
    for n_drive, pid in enumerate(feeds):
        sheets = [sp.sheet(driven=(sp.port_id == pid)) for sp in specs_x]
        dsheet = next(s for s in sheets if s.port_id == pid)
        linear = [lambda test, _s=dsheet: lumped_port_linear_term(
            msh, tags_f, _s, test, omega_rad_per_s=omega)]
        comm.Barrier()
        t1 = time.perf_counter()
        if n_drive > 0:
            fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
        else:
            bil = wall + [
                (lambda trial, test, _s=s: lumped_port_bilinear_term(
                    msh, tags_f, _s, trial, test, omega_rad_per_s=omega))
                for s in sheets
            ]
            fields = solver.solve(current_density=None, project_source=False,
                                  extra_bilinear_terms=bil, extra_linear_terms=linear)
        comm.Barrier()
        t_pair.append(time.perf_counter() - t1)
        e = fields.e_complex
        by_id = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
        cur.append(np.array([by_id[p] for p in order], dtype=np.complex128))
        if "sigma" not in keep:
            keep["sigma"] = fields.sigma_field
        cp = fem.Function(e.function_space, name=f"E_{pid}")
        cp.x.array[:] = e.x.array
        cp.x.scatter_forward()
        keep[pid] = cp
        del fields, e
    times["node_pair"] = float(sum(t_pair))
    s_node = _current_route_s(cur, feeds_idx, z_th, v_s)

    e_a, e_b = keep.pop(feeds[0]), keep.pop(feeds[1])
    # Shape, mesh and element checks -- never object identity (`POST-6` step 1).
    if not (e_a.function_space.mesh is e_b.function_space.mesh is msh
            and e_a.x.array.shape == e_b.x.array.shape
            and str(e_a.function_space.ufl_element()) == str(e_b.function_space.ufl_element())):
        raise RuntimeError("the two feed fields do not share one space")
    t1 = time.perf_counter()
    sense_fields, y_in, by_sense = {}, {}, {}
    for sense, (w_a, w_b) in weights.items():
        ef = fem.Function(e_a.function_space, name=f"E_{sense}")
        ef.x.array[:] = w_a * e_a.x.array + w_b * e_b.x.array
        ef.x.scatter_forward()
        sense_fields[sense] = _SenseField(e_complex=ef, sigma_field=keep["sigma"])
        by_sense[sense] = w_a * cur[0] + w_b * cur[1]
        y_in[sense] = _active_y_in(s_node, (w_a, w_b), z0)
    del e_a, e_b
    times["superposition"] = time.perf_counter() - t1
    self.result = TunedQuadratureResult(
        mesh=msh,
        cell_tags=built["cell_tags"],
        facet_tags=tags_f,
        problem=problem,
        specs=specs_x,
        port_definitions=port_defs,
        feeds=feeds,
        az_feeds=az_feeds,
        az_mirror=float(az_mirror),
        weights=weights,
        s_node=s_node,
        y_in=y_in,
        currents={"order": order, "by_drive": cur, "by_sense": by_sense},
        fields=sense_fields,
        cells=built["cells"],
        times=times,
        conductor_terms=cond_terms,
        conductor_facet_tags=built["mesh_facet_tags"],
        comm=comm,
        degree=self.degree,
        gap_impedance_ohm=complex(z_rc),
        feed_impedance_ohm=complex(z_th),
        source_voltage_v=complex(v_s),
        sheets=built["sheets"],
    )
    return self.result


# ---------------------------------------------------------------------------
# `WF-16` step 9 part B (§9 item 382b, chain U36): ``Study.sar`` and
# ``Study.b1_plus`` / ``b1_plus_at`` on the tuned Study.  Reached through the
# ``_sar_port_study`` / ``_b1_plus_port_study`` hooks inside the existing
# geometry-refusal branches of ``_study_sar`` / ``_study_b1_plus`` (no
# pre-existing route enters those branches with a hook-bearing excitation).
# ---------------------------------------------------------------------------

#: The senses ``TwoPortQuadratureDrive`` superposes (``ccw`` lags with azimuth).
_TUNED_SENSES = ("ccw", "cw")


def _tuned_result(study, comm):
    res = getattr(study, "result", None)
    if isinstance(res, TunedQuadratureResult):
        return res
    return _study_solve_tuned_quadrature(study, comm)


def _tuned_senses(drives):
    ids = list(_TUNED_SENSES) if drives is None else list(drives)
    for s in ids:
        if s not in _TUNED_SENSES:
            raise ValueError(f"drives names senses {_TUNED_SENSES}, got {s!r}")
    return ids


def _tuned_quadrature_sar(
    study,
    mass_g: Any = _MASS_DEFAULT,
    *,
    centres,
    drives=None,
    radius_m: Optional[float] = None,
    quadrature_degree: int = SAR_QUADRATURE_DEGREE,
    comm: MPI.Intracomm = MPI.COMM_WORLD,
) -> dict:
    """``Study.sar`` on a ``TwoPortQuadratureDrive`` Study: ``{sense: [dict per centre]}``.

    The ``mass_g`` / ``radius_m`` validation is ``_study_sar``'s (copied); ρ is
    ``Study._density_field`` (phantom tag only); every ``mass_averaged_sar``
    call is collective on every rank.  A record per sense field, no compliance
    claim.
    """
    from .post.sar import averaging_ball_radius, mass_averaged_sar

    if mass_g is _MASS_DEFAULT:
        mass_g = 10.0 if radius_m is None else None
    if (mass_g is None) == (radius_m is None):
        raise ValueError("give exactly one of mass_g and radius_m")
    if radius_m is None:
        if not np.isfinite(float(mass_g)) or float(mass_g) <= 0.0:
            raise ValueError(f"mass_g must be finite and positive, got {mass_g!r}")
        radius = averaging_ball_radius(
            mass_kg=float(mass_g) / 1000.0, rho=study.density_kg_m3
        )
    else:
        radius = float(radius_m)
        if not np.isfinite(radius) or radius <= 0.0:
            raise ValueError(f"radius_m must be finite and positive, got {radius_m!r}")
    ids = _tuned_senses(drives)
    res = _tuned_result(study, comm)
    rho = getattr(study, "_rho_field", None)
    if rho is None:
        rho = study._density_field(res)
        study._rho_field = rho
    out = {}
    for sense in ids:
        f = res.fields[sense]
        out[sense] = [
            mass_averaged_sar(
                f.e_complex,
                sigma=f.sigma_field,
                rho=rho,
                center=tuple(float(x) for x in c),
                radius=radius,
                comm=res.comm,
                quadrature_degree=int(quadrature_degree),
            )
            for c in centres
        ]
    return out


def _tuned_quadrature_b1_plus(
    study, *, drives=None, comm: MPI.Intracomm = MPI.COMM_WORLD
) -> dict:
    """``Study.b1_plus`` on a ``TwoPortQuadratureDrive`` Study: ``{sense: CG1 B}``.

    ``project_to_cg1(magnetic_flux_density_from_e(e_s, ω))`` per sense field,
    cached on ``study._b1_cg1``; ``Study.b1_plus_at`` is unchanged.
    """
    from .post.faraday import magnetic_flux_density_from_e, project_to_cg1

    ids = _tuned_senses(drives)
    res = _tuned_result(study, comm)
    cache = getattr(study, "_b1_cg1", None)
    if cache is None:
        cache = {}
        study._b1_cg1 = cache
    for sense in ids:
        if sense not in cache:
            cache[sense] = project_to_cg1(
                magnetic_flux_density_from_e(res.fields[sense].e_complex, study.omega)
            )
    return {sense: cache[sense] for sense in ids}


# ---------------------------------------------------------------------------
# `WF-1` step 10a (§9 item 431a, chain U46): ``Study.temperature`` on the tuned
# Study, reached through the ``_temperature_port_study`` hook at the top of
# ``_study_temperature``.  The route is `MAT-15` step 5's (P17), lifted with
# provenance from `main` at 59809ba7,
# tests/validation/test_mat15_f_human_tuned_heating.py: the source (lines
# 346–351: the DG0 average over the phantom tag, the scale applied after the
# average, ruling 361-4), the conduction call (step 1's ``_transient``,
# tests/validation/test_mat15_coil_driven_heating_phantom.py:145–151) and the
# 10 g ball-mean forms and per-snapshot series (lines 426–455 and 539–559).
# Never imported from ``tests``.
# ---------------------------------------------------------------------------

#: ``running_average``'s default window (s), used when ``window_s`` is omitted.
_SIX_MINUTE_WINDOW_S = 360.0


@dataclass
class TunedTemperatureResult:
    """What ``Study.temperature`` returns on a tuned two-port Study (not exported).

    ``balls`` is a list of ``{centre, radius_m, volume_m3, theta_k}`` mappings;
    ``six_minute`` holds ``window_s``, ``volume_mean`` ``{peak_k, t_max_s}``
    and one such pair per ball.  θ = T − t_inf, a record per ``source_scale``;
    no compliance claim.
    """

    times_s: list
    volume_mean_k: list
    balls: list
    six_minute: dict
    source_power_w: float
    max_balance_residual: float
    cells: int
    n_sub: int
    sense: str
    cooling: Any
    transient: Any


def _tuned_quadrature_temperature(study, cooling, *, transient, balls, window_s=None,
                                  sense="ccw", comm: MPI.Intracomm = MPI.COMM_WORLD):
    """Transient Robin-cooled temperature of the tuned Study's phantom.

    Every argument is checked before any mesh; every assembly and reduction
    is collective on every rank (a rank may own no phantom cell).
    """
    if not isinstance(cooling, RobinCooling):
        raise TypeError(f"cooling must be a RobinCooling, got {type(cooling).__name__}")
    if not isinstance(transient, TransientSpec):
        raise TypeError("a tuned two-port Study needs a TransientSpec (transient=)")
    if transient.initial_gradient_k_per_m is not None:  # `WF-1` step 28 (§9 item 830)
        raise ValueError("initial_gradient_k_per_m is defined on a ThermalBoxGeometry only")
    if sense not in _TUNED_SENSES:
        raise ValueError(f"sense must be one of {_TUNED_SENSES}, got {sense!r}")
    if not isinstance(balls, dict) or set(balls) != {"mass_g", "centres"}:
        raise ValueError("balls must be a mapping with mass_g and centres")
    mass_g = float(balls["mass_g"])
    if not np.isfinite(mass_g) or mass_g <= 0.0:
        raise ValueError(f"balls mass_g must be finite and positive, got {mass_g!r}")
    centres = [tuple(float(x) for x in c) for c in balls["centres"]]
    if not centres or any(len(c) != 3 for c in centres):
        raise ValueError("balls centres must be a non-empty list of [x, y, z]")
    window = _SIX_MINUTE_WINDOW_S if window_s is None else float(window_s)
    if not np.isfinite(window) or window <= 0.0:
        raise ValueError(f"window_s must be finite and positive, got {window_s!r}")

    res = _tuned_result(study, comm)  # the kept field: no second solve
    return _phantom_transient_core(
        study, res.mesh, res.cell_tags, res.comm, res.fields[sense], cooling, transient,
        mass_g, centres, window, _exterior_robin_tags, int(res.cells), str(sense))


def _phantom_transient_core(study, msh, cell_tags, rcomm, f, cooling, transient, mass_g,
                            centres, window, robin_fn, cells, sense):
    """The tuned route's source / conduction / ball-mean core on one field ``f``.

    ``robin_fn(sub)`` returns the Robin facet tags (``_ROBIN_TAG``); the tuned
    two-port Study passes ``_exterior_robin_tags``, the rod four-port
    ``_rod_outer_robin_tags`` (`WF-15` step 14, §9 item 831).
    """
    import ufl
    from dolfinx import default_scalar_type, fem

    from .post.sar import averaging_ball_radius
    from .thermal.dose import running_average
    from .thermal.transient import solve_transient_conduction

    tag = study.geometry.phantom_tag
    sub, parent_ids = _tag_submesh(msh, cell_tags, tag)
    ft = robin_fn(sub)
    Wp, avg = _dg0_cell_average(
        msh, cell_tags, 0.5 * f.sigma_field * ufl.inner(f.e_complex, f.e_complex), tag
    )
    avg *= cooling.source_scale  # after the average (ruling 361-4)
    q = _dg0_to_submesh(Wp, avg, sub, parent_ids)

    def _red(form):
        return float(np.real(rcomm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = _red(fem.form(one * ufl.dx(domain=sub)))
    source_power = _red(fem.form(q * ufl.dx(domain=sub)))
    n_sub = int(sub.topology.index_map(sub.topology.dim).size_global)

    n_steps = int(transient.n_steps)
    dt = float(transient.dt_s)
    t_inf = float(cooling.t_inf_c)
    t_init = t_inf if transient.t_init_c is None else float(transient.t_init_c)
    steps = tuple(range(n_steps + 1))
    run = solve_transient_conduction(
        sub, ft, k=cooling.k_w_per_m_k, rho_c=transient.rho_c_j_per_m3_k, q=q,
        h=cooling.h_w_per_m2_k, t_inf=t_inf, t_init=t_init, robin_tag=_ROBIN_TAG,
        dt=dt, n_steps=n_steps, theta=float(transient.theta),
        startup_half_steps=int(transient.startup_half_steps), degree=int(cooling.degree),
        snapshot_steps=steps,
    )
    snaps = run.snapshots

    radius = float(averaging_ball_radius(mass_kg=mass_g * 1e-3, rho=study.density_kg_m3))
    centre_c = fem.Constant(sub, np.zeros(3, dtype=default_scalar_type))
    xs = ufl.SpatialCoordinate(sub)
    off = xs - centre_c
    inside = ufl.conditional(ufl.lt(ufl.real(ufl.dot(off, off)), radius**2), 1.0, 0.0)
    dxq = ufl.dx(domain=sub, metadata={"quadrature_degree": int(SAR_QUADRATURE_DEGREE)})
    th = fem.Function(snaps[0].function_space)
    f_vol = fem.form(inside * dxq)
    f_th = fem.form(th * inside * dxq)
    f_mean = fem.form(th * ufl.dx(domain=sub))

    vols = []
    for c in centres:
        centre_c.value[:] = np.array(c)
        vols.append(_red(f_vol))
    mean_series = []
    ball_series = [[] for _ in centres]
    for n in steps:
        th.x.array[:] = snaps[n].x.array - t_inf
        th.x.scatter_forward()
        mean_series.append(_red(f_mean) / v_sub)
        for i, c in enumerate(centres):
            centre_c.value[:] = np.array(c)
            ball_series[i].append(_red(f_th) / vols[i])
    times = [n * dt for n in steps]

    def _six(series):
        _means, peak, t_max = running_average(times, np.asarray(series), window=window)
        return {"peak_k": float(peak), "t_max_s": float(t_max)}

    return TunedTemperatureResult(
        times_s=times,
        volume_mean_k=mean_series,
        balls=[
            {"centre": list(c), "radius_m": radius, "volume_m3": vols[i],
             "theta_k": ball_series[i]}
            for i, c in enumerate(centres)
        ],
        # ``window`` None (the rod four-port, `WF-15` step 14): no running maximum.
        six_minute={"window_s": None} if window is None else {
            "window_s": window,
            "volume_mean": _six(mean_series),
            "balls": [_six(s) for s in ball_series],
        },
        source_power_w=source_power,
        max_balance_residual=float(run.max_balance_residual),
        cells=int(cells),
        n_sub=n_sub,
        sense=str(sense),
        cooling=cooling,
        transient=transient,
    )


# ---------------------------------------------------------------------------
# `WF-15` step 14 (§9 item 831, §10 U127): ``Study.temperature(drive=...)`` on
# the PEC implant-rod four-port, on `IMP-3` I41's route.  The Robin selection
# is I39's ``_outer_robin_tags`` (copied from `main` at 2f95023f,
# tests/validation/test_imp3_rod_heating_phantom.py:115–140): exterior facets
# of the phantom submesh farther than ``_ROD_WALL_SELECT_M`` from the rod's
# capsule axis segment are cooled, the rod wall (adiabatic) is not.  Never
# ``_exterior_robin_tags`` here, which would cool the rod wall.
# ---------------------------------------------------------------------------

#: I39's ``WALL_SELECT_M``: wall at r = 1 mm, nearest outer phantom surface >= 15 mm.
_ROD_WALL_SELECT_M = 0.005


def _rod_outer_robin_counts(sub, rod):
    """``(outer facets, {exterior, outer, wall} owned counts reduced)`` on ``sub``."""
    from dolfinx import mesh as dmesh

    tdim = sub.topology.dim
    fdim = tdim - 1
    sub.topology.create_connectivity(fdim, tdim)
    ext = np.sort(np.asarray(dmesh.exterior_facet_indices(sub.topology), dtype=np.int32))
    c = np.asarray(rod["centre"], dtype=float)
    a = np.asarray(rod["axis"], dtype=float)
    a = a / np.linalg.norm(a)
    half = 0.5 * float(rod["length"])  # cylindrical section between the cap centres
    if ext.size:
        x = dmesh.compute_midpoints(sub, fdim, ext) - c
        s = np.clip(x @ a, -half, half)
        d = np.linalg.norm(x - s[:, None] * a[None, :], axis=1)
        keep = d > _ROD_WALL_SELECT_M
    else:
        keep = np.zeros(0, dtype=bool)
    outer = ext[keep].astype(np.int32)
    wall = ext[~keep]
    n_owned = sub.topology.index_map(fdim).size_local
    counts = {
        "exterior": sub.comm.allreduce(int(np.count_nonzero(ext < n_owned)), op=MPI.SUM),
        "outer": sub.comm.allreduce(int(np.count_nonzero(outer < n_owned)), op=MPI.SUM),
        "wall": sub.comm.allreduce(int(np.count_nonzero(wall < n_owned)), op=MPI.SUM),
    }
    return outer, counts


def _rod_outer_robin_tags(sub, rod):
    from dolfinx import mesh as dmesh

    outer, _counts = _rod_outer_robin_counts(sub, rod)
    return dmesh.meshtags(sub, sub.topology.dim - 1, outer,
                          np.full(outer.size, _ROBIN_TAG, dtype=np.int32))


def _rod_port_temperature(study, cooling, *, transient, balls, drive,
                          comm: MPI.Intracomm = MPI.COMM_WORLD):
    """Transient Robin-cooled temperature of the rod four-port's phantom, one port drive.

    The field is ``study._solve_port_fields(comm).fields[drive]`` (one port at
    1 V, the page's convention); balls are 1 g by default via ``mass_g``.  A
    record per ``source_scale``; no compliance claim.
    """
    if not isinstance(cooling, RobinCooling):
        raise TypeError(f"cooling must be a RobinCooling, got {type(cooling).__name__}")
    if not isinstance(transient, TransientSpec):
        raise TypeError("the rod four-port Study needs a TransientSpec (transient=)")
    if transient.initial_gradient_k_per_m is not None:  # `WF-1` step 28 (§9 item 830)
        raise ValueError("initial_gradient_k_per_m is defined on a ThermalBoxGeometry only")
    if not isinstance(drive, str) or not drive:
        raise ValueError(f"drive must name one port, got {drive!r}")
    if not isinstance(balls, dict) or set(balls) != {"mass_g", "centres"}:
        raise ValueError("balls must be a mapping with mass_g and centres")
    mass_g = float(balls["mass_g"])
    if not np.isfinite(mass_g) or mass_g <= 0.0:
        raise ValueError(f"balls mass_g must be finite and positive, got {mass_g!r}")
    centres = [tuple(float(x) for x in c) for c in balls["centres"]]
    if not centres or any(len(c) != 3 for c in centres):
        raise ValueError("balls centres must be a non-empty list of [x, y, z]")
    res = study._solve_port_fields(comm)
    if drive not in res.fields:
        raise ValueError(f"drive must be one of {sorted(res.fields)}, got {drive!r}")
    rod = dict(study.geometry.implant_rod)
    return _phantom_transient_core(
        study, res.mesh, res.cell_tags, res.comm, res.fields[drive], cooling, transient,
        mass_g, centres, None, lambda sub: _rod_outer_robin_tags(sub, rod),
        int(res.cells), drive)
