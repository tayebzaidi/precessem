"""`TH-17` — eigenmodes of the loaded F-small birdcage with the PEC coil and
`PORT-15` step 3's tuned capacitor sheets.

Fixture: `TH-15` step 3's hole mesh (``_sheets_build(True, as_hole=True)`` —
the coil cut from the air, its cavity wall exterior and therefore PEC under the
``pec_facet_tags=None`` reasoning that fixture documents), the phantom at
``SALINE_SIGMA`` / ``SALINE_EPSILON_R``, and the four gap sheets carrying
``Z_p = 1/(jωC)`` at `PORT-15` step 3's ``C_tuned`` (**imported** from
``tests.validation.test_port_circuit_layer_field`` — ``select_c_tuned`` on
``tuning_sweep``, never restated here).  No source, no termination: the whole
problem is the one non-Hermitian generalised pencil (G1) of
:func:`fem_em_solver.core.cavity.solve_general_mesh_modes`.

Step 1 lands the **cost probe** first (§5.1, §9 item 2): this mesh has never
carried an eigensolve — `TH-9` and `TH-14` were boxes.  ``TH17_PROBE=1`` runs
mesh + one shift-invert solve at the 64 MHz target with ``nev = 6`` and prints
cells / dofs / wall / ``ru_maxrss`` per rank.  The item's **pre-registered STOP
rule**: > 15 min or > 40 GiB at ``-n 4`` ⇒ the step re-prices to an XL cost
probe and nothing else runs.
"""

from __future__ import annotations

import os
import resource
import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from dolfinx.fem.petsc import assemble_matrix
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core.cavity import (
    complex_relative_permittivity,
    solve_general_mesh_modes,
)
from fem_em_solver.ports.lumped import (
    LumpedPortSheet,
    lumped_port_bilinear_term,
    series_rlc_impedance,
)
from fem_em_solver.utils.constants import C_0
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import (
    _projected_extents,
    _sheet_areas,
    _sheet_azimuth_deg,
)
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
from tests.mesh.test_birdcage_port_sheets import _build as _sheets_build
from tests.mesh.test_birdcage_port_tags import LEG_COUNT
from tests.mesh.test_two_torus_port_facets import _facet_group_area
from tests.mesh.test_two_torus_port_sheet import _sheet_facet_count
from tests.validation.test_port_birdcage_leg_offset_sweep import (
    _narrowed_radial,
    _port_frame,
)
from tests.validation.test_port_birdcage_lumped_column import (
    PHANTOM_CELL_TAG,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port_circuit_layer_field import (
    S_64MHZ_EPS0_RECORD,
    STEP3_REGISTERED_FREQUENCY_HZ,
    select_c_tuned,
    tuning_sweep,
)
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

PROBE_ENV = "TH17_PROBE"
#: Step 1b's two diagnostics (§9 item 10 (i) and (ii)) — env-gated exactly as the
#: probe is, so the default collection of this module is unchanged.
STEP1B_ENV = "TH17_STEP1B"
#: (ii) the eight shift-invert targets [Hz], the item's list verbatim.
SCAN_TARGETS_HZ = (1e6, 4e6, 16e6, 32e6, 64e6, 128e6, 256e6, 512e6)
#: (ii) only converged eigenvalues above this are interesting — below it is the
#: N1curl gradient (null-space) cluster step 1 part (i) found (`Re f` ≈ 2 kHz).
SCAN_PRINT_FLOOR_HZ = 1.0e6
#: `PORT-15` step 3's registered tuning frequency is the eigen target.
TARGET_FREQUENCY_HZ = STEP3_REGISTERED_FREQUENCY_HZ
NEV = 6


def c_tuned_farad() -> float:
    """`PORT-15` step 3's ``C_tuned`` — imported, recomputed, never restated."""

    _, _, roots = tuning_sweep(S_64MHZ_EPS0_RECORD, STEP3_REGISTERED_FREQUENCY_HZ)
    tuned = select_c_tuned(roots)
    assert tuned is not None, "PORT-15 step 3 tuning found no zero of Im Z_in"
    return float(tuned["c_f"])


def build_hole_fixture():
    """`TH-15` step 3's hole mesh + its four narrowed gap sheets' geometry."""

    comm = MPI.COMM_WORLD
    ports_idx = list(range(1, LEG_COUNT + 1))
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags, _diag, t_mesh = _sheets_build(True, as_hole=True)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()

    tags_f, _full_areas = _sheet_areas(msh, cell_tags, ports_idx, comm)
    geometry = {}
    for i in ports_idx:
        tag = SHEET_IFACE + i
        azimuth = _sheet_azimuth_deg(msh, tags_f, tag, comm)
        radial, azimuthal, axial = _port_frame(azimuth)
        spans, centres = _projected_extents(msh, tags_f, tag, comm, [radial, azimuthal, axial])
        geometry[i] = {
            "tag": tag,
            "azimuth_deg": float(azimuth),
            "radial": radial,
            "w_full": float(spans[0]),
            "centre_radial": float(centres[0]),
        }
    for i in ports_idx:
        g = geometry[i]
        tags_f = _narrowed_radial(
            msh,
            tags_f,
            g["tag"],
            GATED_WIDTH_FRACTION,
            g["centre_radial"],
            g["radial"],
            0.5 * g["w_full"],
        )
    sheets = []
    for i in ports_idx:
        g = geometry[i]
        n_facets = _sheet_facet_count(msh, tags_f, g["tag"], comm)
        assert n_facets > 0, f"sheet {g['tag']}: no owned facets anywhere"
        area = _facet_group_area(msh, tags_f, g["tag"], comm)
        radial, azimuthal, axial = _port_frame(g["azimuth_deg"])
        spans, _c = _projected_extents(msh, tags_f, g["tag"], comm, [radial, azimuthal, axial])
        sheets.append(
            {
                **g,
                "facets": int(n_facets),
                "area": float(area),
                "h": float(spans[2]),
                "w": float(area / spans[2]),
            }
        )
    return {
        "mesh": msh,
        "cell_tags": cell_tags,
        "facet_tags": facet_tags,
        "sheet_facet_tags": tags_f,
        "sheets": sheets,
        "n_cells": int(msh.topology.index_map(tdim).size_global),
        "t_mesh_s": float(t_mesh),
        "t_fixture_s": float(time.perf_counter() - t0),
    }


def permittivity_field(msh, cell_tags, omega_lin_rad_s, sigma_s_per_m=None):
    """DG0 ``ε_r,c``: the phantom linearised at ``ω_lin``, air elsewhere.

    ``sigma_s_per_m`` (default ``None`` ⇒ ``SALINE_SIGMA``, byte-identical to
    every earlier caller) is the phantom conductivity; gate (e) passes 0.0.
    """

    sigma = SALINE_SIGMA if sigma_s_per_m is None else float(sigma_s_per_m)
    Q = fem.functionspace(msh, ("DG", 0))
    eps = fem.Function(Q, name="epsilon_r_complex")
    eps.x.array[:] = 1.0 + 0.0j
    phantom_cells = np.asarray(cell_tags.find(PHANTOM_CELL_TAG), dtype=np.int32)
    if phantom_cells.size:
        dofs = fem.locate_dofs_topological(Q, msh.topology.dim, phantom_cells)
        eps.x.array[dofs] = complex_relative_permittivity(SALINE_EPSILON_R, sigma, omega_lin_rad_s)
    eps.x.scatter_forward()
    return eps


def _capacitor_mass_forms(
    msh, sheet_facet_tags, sheets, c_f, trial, test, *, omega_ref_rad_s, live=None, z_ref_hz=None
):
    """``−(c²/ω_ref²)·L1(ω_ref)`` per sheet — the (L1) law, moved onto ``k₀²``.

    ``live`` (default: all) selects which sheets carry the capacitor; the
    shuffled-sheet negative control passes a subset of tags.

    ``z_ref_hz`` (default ``None`` ⇒ ``TARGET_FREQUENCY_HZ``, byte-identical to
    every earlier caller) is the frequency ``Z_p = 1/(jωC)`` is evaluated at.
    The mass is ``∝ C`` only when it equals ``ω_ref/2π``; otherwise it scales by
    ``ω_z/ω_ref`` — step 2's fixed-point re-solve at ``ω_lin ≠ 2π·64 MHz``
    passes ``ω_ref/2π`` so the four capacitors stay ``C``.
    """

    z_hz = TARGET_FREQUENCY_HZ if z_ref_hz is None else float(z_ref_hz)
    z_p = series_rlc_impedance(z_hz, c_f=float(c_f))
    scale = -(C_0**2) / (float(omega_ref_rad_s) ** 2)
    forms = []
    for s in sheets:
        if live is not None and s["tag"] not in live:
            continue
        sheet = LumpedPortSheet(
            port_id=f"P{s['tag'] - SHEET_IFACE}",
            facet_tag=int(s["tag"]),
            port_impedance_ohm=z_p,
            gap_height_m=s["h"],
            sheet_width_m=s["w"],
            drive_direction=(0.0, 0.0, 1.0),
            source_voltage_v=0.0 + 0.0j,
            interior=True,
        )
        forms.append(
            scale
            * lumped_port_bilinear_term(
                msh,
                sheet_facet_tags,
                sheet,
                trial,
                test,
                omega_rad_per_s=float(omega_ref_rad_s),
            )
        )
    return forms


def solve_modes(
    fx,
    c_f,
    *,
    omega_lin_rad_s,
    degree=1,
    nev=NEV,
    live=None,
    target_hz=TARGET_FREQUENCY_HZ,
    z_ref_hz=None,
    sigma_s_per_m=None,
    return_vectors=False,
):
    """One (G1) solve; returns ``(eigenvalues, n_converged, n_dofs, wall_s)``.

    ``return_vectors`` (default ``False``: unchanged; `EX-76`) appends one
    element, ``{"vectors", "V", "eps", "forms", "k0_sq"}``: the eigenvectors in
    the solver's order plus what re-assembles the pencil.

    ``target_hz`` (default: the registered 64 MHz target, so step 1's probe is
    unchanged) is the shift-invert target; step 1b (ii) sweeps it.
    ``z_ref_hz`` is passed through to :func:`_capacitor_mass_forms` (default
    ``None``: unchanged).  ``sigma_s_per_m`` is passed through to
    :func:`permittivity_field` (default ``None`` ⇒ ``SALINE_SIGMA``: unchanged).
    """

    msh = fx["mesh"]
    V = fem.functionspace(msh, ("N1curl", degree))
    eps = permittivity_field(msh, fx["cell_tags"], omega_lin_rad_s, sigma_s_per_m=sigma_s_per_m)
    trial, test = ufl.TrialFunction(V), ufl.TestFunction(V)
    forms = _capacitor_mass_forms(
        msh,
        fx["sheet_facet_tags"],
        fx["sheets"],
        c_f,
        trial,
        test,
        omega_ref_rad_s=omega_lin_rad_s,
        live=live,
        z_ref_hz=z_ref_hz,
    )
    k0_sq = (2.0 * np.pi * float(target_hz) / C_0) ** 2
    t0 = time.perf_counter()
    if return_vectors:
        values, n_converged, n_bc, vectors = solve_general_mesh_modes(
            V,
            eps,
            target_k0_sq=k0_sq,
            sheet_mass_forms=forms,
            nev=nev,
            comm=msh.comm,
            return_vectors=True,
        )
    else:
        values, n_converged, n_bc = solve_general_mesh_modes(
            V,
            eps,
            target_k0_sq=k0_sq,
            sheet_mass_forms=forms,
            nev=nev,
            comm=msh.comm,
        )
    wall = float(time.perf_counter() - t0)
    n_dofs = int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)
    if return_vectors:
        return (
            values,
            int(n_converged),
            n_dofs,
            wall,
            int(n_bc),
            {"vectors": vectors, "V": V, "eps": eps, "forms": forms, "k0_sq": k0_sq},
        )
    return values, int(n_converged), n_dofs, wall, int(n_bc)


def frequency_of(lam: complex) -> complex:
    """``f = c√λ / 2π`` [Hz], complex."""

    return C_0 * np.sqrt(complex(lam)) / (2.0 * np.pi)


@complex_only
@pytest.mark.skipif(os.environ.get(PROBE_ENV) != "1", reason=f"set {PROBE_ENV}=1")
def test_cost_probe():
    """§9 item 2 (i): the unmeasured cost of one (G1) solve on the hole mesh.

    Prints cells / dofs / wall / ``ru_maxrss`` per rank and evaluates the
    item's pre-registered STOP rule (> 15 min or > 40 GiB at ``-n 4``).
    Asserts nothing but that an eigenpair converged — this is a probe.
    """

    comm = MPI.COMM_WORLD
    t_all = time.perf_counter()
    c_f = c_tuned_farad()
    fx = build_hole_fixture()
    omega_lin = 2.0 * np.pi * TARGET_FREQUENCY_HZ
    values, n_converged, n_dofs, wall, n_bc = solve_modes(fx, c_f, omega_lin_rad_s=omega_lin)
    total = float(time.perf_counter() - t_all)

    rss_gib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024.0**2)
    all_rss = comm.allgather(float(rss_gib))
    summed = float(sum(all_rss))

    if comm.rank == 0:
        print(
            f"\n[TH-17 probe] ranks = {comm.size}, cells = {fx['n_cells']}, "
            f"dofs (N1curl deg 1) = {n_dofs}, constrained (rank-local) = {n_bc}"
        )
        print(
            f"[TH-17 probe] mesh {fx['t_mesh_s']:.2f} s, fixture {fx['t_fixture_s']:.2f} s, "
            f"eigensolve {wall:.2f} s, probe total {total:.2f} s"
        )
        print(
            "[TH-17 probe] ru_maxrss per rank [GiB]: "
            + ", ".join(f"{v:.3f}" for v in all_rss)
            + f"; summed {summed:.3f}, max {max(all_rss):.3f}"
        )
        print(
            f"[TH-17 probe] C_tuned = {c_f:.15e} F (imported, PORT-15 step 3); "
            f"target f = {TARGET_FREQUENCY_HZ:.6e} Hz"
        )
        print(f"[TH-17 probe] n_converged = {n_converged} of nev = {NEV}")
        for i, lam in enumerate(values[:NEV]):
            f_c = frequency_of(lam)
            q = abs(f_c.real / (2.0 * f_c.imag)) if f_c.imag != 0 else float("inf")
            print(
                f"    lam[{i}] = {lam.real:+.9e} {lam.imag:+.9e}j  "
                f"f = {f_c.real:.9e} {f_c.imag:+.9e}j Hz  Q = {q:.6e}"
            )
        stop_time = total > 15.0 * 60.0
        stop_mem = summed > 40.0
        print(
            f"[TH-17 probe] STOP rule (> 15 min or > 40 GiB at -n 4): "
            f"time {total / 60.0:.3f} min -> {'STOP' if stop_time else 'ok'}; "
            f"summed rss {summed:.3f} GiB -> {'STOP' if stop_mem else 'ok'}"
        )
        print(f"[TH-17 probe] VERDICT: {'STOP' if (stop_time or stop_mem) else 'PROCEED'}")
    assert n_converged > 0, "no eigenpair converged near the 64 MHz target"


# ---------------------------------------------------------------------------
# Step 1b (§9 item 10) — the two diagnostics, printed and never gated.
#
# Step 1 part (i) found only the N1curl gradient cluster at the 64 MHz target
# (`Re f` ~ 2 kHz, `f085e51`).  Two things have to be ruled out before the
# formulation is blamed: that the (L1) sheet forms never reached the pencil at
# all (facet tags / restriction — diagnostic (i) below), and that the loaded
# modes simply sit somewhere other than 64 MHz (diagnostic (ii)).  Neither
# asserts a physical number: the item's asserted anchor is the closed-form LC
# loop control (iii), which is a different fixture.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def step1b_fixture():
    """One hole-mesh build shared by both step 1b diagnostics."""

    return build_hole_fixture()


def _sheet_dof_indicator(V, fx):
    """A ``V``-Function that is 1 on the four gap sheets' dofs and 0 elsewhere.

    Returns ``(indicator, n_sheet_dofs_global_owned)``; the owned count is
    reduced across ranks (``locate_dofs_topological`` is rank-local).
    """

    msh = V.mesh
    tdim = msh.topology.dim
    tags = fx["sheet_facet_tags"]
    pieces = [np.asarray(tags.find(int(s["tag"])), dtype=np.int32) for s in fx["sheets"]]
    facets = np.unique(np.concatenate(pieces)) if pieces else np.empty(0, dtype=np.int32)
    dofs = fem.locate_dofs_topological(V, tdim - 1, facets)
    ind = fem.Function(V)
    ind.x.array[:] = 0.0
    ind.x.array[dofs] = 1.0
    ind.x.scatter_forward()
    n_owned_local = int(
        np.count_nonzero(dofs < V.dofmap.index_map.size_local * V.dofmap.index_map_bs)
    )
    return ind, int(msh.comm.allreduce(n_owned_local, op=MPI.SUM))


def _quadratic_form(mat, vec):
    """``|xᴴ M x|`` — the mass a matrix carries on the indicator's support."""

    out = mat.createVecLeft()
    mat.mult(vec, out)
    value = complex(vec.dot(out))
    out.destroy()
    return abs(value), value


@complex_only
@pytest.mark.skipif(os.environ.get(STEP1B_ENV) != "1", reason=f"set {STEP1B_ENV}=1")
def test_step1b_sheet_mass_norms(step1b_fixture):
    """(i) ``‖B_sheet‖`` against ``‖B‖`` on the sheet dofs — printed, not gated.

    A zero (or a ratio at round-off) says the (L1) forms never reached the
    pencil — the facet tags or the interior restriction — and is the first thing
    the item rules out.  The only assertion is that the sheet term is nonzero and
    finite, which is the premise every later step of `TH-17` rests on; the ratio
    itself is a diagnostic reading and carries no band.
    """

    fx = step1b_fixture
    comm = fx["mesh"].comm
    c_f = c_tuned_farad()
    omega_lin = 2.0 * np.pi * TARGET_FREQUENCY_HZ
    V = fem.functionspace(fx["mesh"], ("N1curl", 1))
    eps = permittivity_field(fx["mesh"], fx["cell_tags"], omega_lin)
    trial, test = ufl.TrialFunction(V), ufl.TestFunction(V)

    volume_ufl = ufl.inner(eps * trial, test) * ufl.dx
    forms = _capacitor_mass_forms(
        fx["mesh"],
        fx["sheet_facet_tags"],
        fx["sheets"],
        c_f,
        trial,
        test,
        omega_ref_rad_s=omega_lin,
    )
    assert forms, "no capacitor sheet forms were built at all"
    sheet_ufl = forms[0]
    for form in forms[1:]:
        sheet_ufl = sheet_ufl + form

    b_vol = assemble_matrix(fem.form(volume_ufl))
    b_vol.assemble()
    b_sheet = assemble_matrix(fem.form(sheet_ufl))
    b_sheet.assemble()

    n_vol = float(b_vol.norm(PETSc.NormType.FROBENIUS))
    n_sheet = float(b_sheet.norm(PETSc.NormType.FROBENIUS))

    ind, n_sheet_dofs = _sheet_dof_indicator(V, fx)
    x = ind.x.petsc_vec
    q_vol, q_vol_c = _quadratic_form(b_vol, x)
    q_sheet, q_sheet_c = _quadratic_form(b_sheet, x)

    n_dofs = int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)
    if comm.rank == 0:
        print(
            f"\n[TH-17 1b(i)] cells = {fx['n_cells']}, N1curl deg-1 dofs = {n_dofs}, "
            f"sheet dofs (owned, reduced) = {n_sheet_dofs}, "
            f"C_tuned = {c_f:.9e} F"
        )
        for s in fx["sheets"]:
            print(
                f"    sheet tag {s['tag']}: facets {s['facets']}, area {s['area']:.6e} m^2, "
                f"h {s['h']:.6e} m, w {s['w']:.6e} m, h/w {s['h'] / s['w']:.6e}"
            )
        print(
            f"[TH-17 1b(i)] Frobenius: ||B_volume|| = {n_vol:.9e}, "
            f"||B_sheet|| = {n_sheet:.9e}, ratio = {n_sheet / n_vol:.9e}"
        )
        print(
            f"[TH-17 1b(i)] on the sheet-dof indicator: x^H B_volume x = {q_vol_c:.9e}, "
            f"x^H B_sheet x = {q_sheet_c:.9e}, ratio |sheet/volume| = "
            f"{(q_sheet / q_vol if q_vol else float('inf')):.9e}"
        )
    b_vol.destroy()
    b_sheet.destroy()

    assert (
        np.isfinite(n_sheet) and n_sheet > 0.0
    ), f"||B_sheet|| = {n_sheet!r}: the (L1) sheet forms never reached the pencil"
    assert (
        np.isfinite(q_sheet) and q_sheet > 0.0
    ), f"x^H B_sheet x = {q_sheet!r} on the sheet dofs: the restriction is empty"


@complex_only
@pytest.mark.skipif(os.environ.get(STEP1B_ENV) != "1", reason=f"set {STEP1B_ENV}=1")
def test_step1b_target_scan(step1b_fixture):
    """(ii) shift-invert at eight targets; every converged `Re f` > 1 MHz printed.

    Printed, never gated (the item): this locates where the capacitor-loaded
    modes of the F-small hole mesh actually sit, so a later step can say whether
    the circuit layer's `C_tuned` puts the *eigen* resonance at 64 MHz.  The one
    assertion is that at least one target converged something — a scan that
    converges nothing is a solver report, not a spectrum.
    """

    fx = step1b_fixture
    comm = fx["mesh"].comm
    c_f = c_tuned_farad()
    omega_lin = 2.0 * np.pi * TARGET_FREQUENCY_HZ
    total_converged = 0
    above_floor = []
    for target in SCAN_TARGETS_HZ:
        values, n_converged, n_dofs, wall, _n_bc = solve_modes(
            fx, c_f, omega_lin_rad_s=omega_lin, target_hz=target
        )
        total_converged += int(n_converged)
        k0_sq = (2.0 * np.pi * float(target) / C_0) ** 2
        if comm.rank == 0:
            print(
                f"\n[TH-17 1b(ii)] target {target / 1e6:.0f} MHz "
                f"(k0^2 = {k0_sq:.6e}), nev = {NEV}, converged = {n_converged}, "
                f"{wall:.2f} s, dofs = {n_dofs}"
            )
        for lam in values:
            f_c = frequency_of(lam)
            if abs(f_c.real) <= SCAN_PRINT_FLOOR_HZ:
                continue
            above_floor.append((float(target), complex(lam), complex(f_c)))
            q = abs(f_c.real / (2.0 * f_c.imag)) if f_c.imag != 0 else float("inf")
            if comm.rank == 0:
                print(
                    f"    lam = {lam.real:+.9e} {lam.imag:+.9e}j  "
                    f"f = {f_c.real:.9e} {f_c.imag:+.9e}j Hz  Q = {q:.6e}  "
                    f"|lam - target| = {abs(lam - k0_sq):.6e}"
                )
        if comm.rank == 0 and not any(t == float(target) for t, _l, _f in above_floor):
            print(
                f"    (nothing above {SCAN_PRINT_FLOOR_HZ / 1e6:.0f} MHz at this target "
                f"— gradient cluster only)"
            )
    if comm.rank == 0:
        print(
            f"\n[TH-17 1b(ii)] summary: {len(above_floor)} converged eigenvalue(s) "
            f"with Re f > {SCAN_PRINT_FLOOR_HZ / 1e6:.0f} MHz across "
            f"{len(SCAN_TARGETS_HZ)} targets"
        )
        for t, lam, f_c in above_floor:
            print(
                f"    target {t / 1e6:.0f} MHz -> Re f = {f_c.real:.9e} Hz, "
                f"lam = {lam.real:+.9e} {lam.imag:+.9e}j"
            )
    assert total_converged > 0, "no eigenpair converged at any of the eight targets"


# ---------------------------------------------------------------------------
# Step 1c (§9 item 24) — the closed-form LC control on the two-torus hole
# fixture.  One gapped PEC loop (torus 1's hole, its sheet 211 carrying the
# capacitor mass form exactly as `_capacitor_mass_forms` builds it) is an LC
# circuit: `f_LC = 1/(2π√(L·C))`.  `L` is this fixture's own driven reading,
# `L_fem = Im Z_in/ω` of sheet 211 at the target, so both sides of the 5 % are
# this mesh's numbers and the box's effect on `L` cancels.  Torus 2's gap is
# left open in both solves (no sheet form in the pencil; ≥ 1e9 Ω driven).
#
# Note on the facet tag: the item text names "facet tag 201"; on the
# ``as_hole=True`` route the gap<->conductor groups 201/202 do not exist
# (`io/mesh.py` docstring, `as_hole`) — torus 1's port **sheet** is 211, the
# group every lumped-sheet caller of this fixture uses (``SHEET_FACET_TAGS``).
# ---------------------------------------------------------------------------

STEP1C_ENV = "TH17_STEP1C"
#: Open-gap termination for torus 2's sheet in the driven solve (the item: ≥ 1e9 Ω).
STEP1C_OPEN_OHM = 1.0e9
#: The driven P1 sheet's source impedance — any finite value; Z_in is V/I.
STEP1C_SOURCE_OHM = 50.0
#: The two driven frequencies (both printed; the target's `L` is the one used).
STEP1C_L_FREQUENCIES_HZ = (10.0e6, TARGET_FREQUENCY_HZ)
#: Anchor (a): the eigenvalue nearest f_LC within 5 % (item 24, registered).
STEP1C_FLC_BAND = 0.05
#: Anchor (b): at 2·C the same mode moves by 1/√2 within 2 % (item 24).
STEP1C_SQRT2_BAND = 0.02
#: Negative control: no eigenvalue with Re f in f_LC·[0.8, 1.2] without the form.
STEP1C_CONTROL_WINDOW = (0.8, 1.2)


def build_two_torus_hole_fixture():
    """`TH-15` step 2a's two-torus PEC-hole mesh + sheet 211/212 geometry.

    Every constant is imported (`test_port_reaction_impedance` for the torus
    and box, `test_port_gap_voltage_impedance` for the gap), never restated.
    """

    from fem_em_solver.io.mesh import MeshGenerator
    from tests.mesh.test_two_torus_port_sheet import SHEET_FACET_TAGS, _sheet_extents
    from tests.validation.test_port_gap_voltage_impedance import (
        GAP_ANGLE,
        GAP_ARC_RESOLUTION,
        GAP_BURIAL,
        GAP_OVERHANG,
    )
    from tests.validation.test_port_reaction_impedance import (
        AIR_PADDING,
        H_FAR,
        H_WIRE,
        MAJOR_RADIUS,
        MINOR_RADIUS,
        SEPARATION,
    )

    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags = MeshGenerator.two_torus_domain(
        separation=SEPARATION,
        major_radius=MAJOR_RADIUS,
        minor_radius=MINOR_RADIUS,
        resolution=H_FAR,
        air_padding=AIR_PADDING,
        wire_resolution=H_WIRE,
        far_resolution=H_FAR,
        port_gap=True,
        gap_angle=GAP_ANGLE,
        gap_burial=GAP_BURIAL,
        gap_overhang=GAP_OVERHANG,
        gap_arc_resolution=GAP_ARC_RESOLUTION,
        emit_port_sheet=True,
        as_hole=True,
        comm=comm,
    )
    t_mesh = time.perf_counter() - t0
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    sheets = []
    for tag in SHEET_FACET_TAGS:
        tag = int(tag)
        n_facets = _sheet_facet_count(msh, facet_tags, tag, comm)
        assert n_facets > 0, f"sheet {tag}: no owned facets anywhere"
        area = _facet_group_area(msh, facet_tags, tag, comm)
        extents = _sheet_extents(msh, facet_tags, tag, comm)
        h = float(extents[1])
        sheets.append(
            {"tag": tag, "facets": int(n_facets), "area": float(area), "h": h, "w": float(area / h)}
        )
    return {
        "mesh": msh,
        "cell_tags": cell_tags,
        "facet_tags": facet_tags,
        "sheets": sheets,
        "n_cells": int(comm.allreduce(msh.topology.index_map(tdim).size_local, op=MPI.SUM)),
        "t_mesh_s": float(t_mesh),
    }


def _driven_loop_impedance(fx, frequency_hz):
    """P1 driven, P2 open: ``(Z_in = V1/I1, Z21/Z11, wall_s)`` off the solved field."""

    from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem
    from fem_em_solver.ports.definitions import PortDefinition
    from fem_em_solver.ports.lumped import LumpedSheetPortSpec, run_lumped_sheet_port_case

    msh = fx["mesh"]
    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=float(frequency_hz),
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=fx["cell_tags"],
        boundary_condition="pec_zero_tangential_a",
        facet_tags=fx["facet_tags"],
        pec_facet_tags=None,
    )
    ports, specs = [], []
    for k, s in enumerate(fx["sheets"]):
        pid = f"P{k + 1}"
        ports.append(
            PortDefinition(
                port_id=pid,
                positive_tag=int(s["tag"]),
                negative_tag=int(s["tag"]),
                orientation="gap_azimuthal_plus_y",
                z0_ohm=STEP1C_SOURCE_OHM,
            )
        )
        specs.append(
            LumpedSheetPortSpec(
                port_id=pid,
                facet_tag=int(s["tag"]),
                port_impedance_ohm=STEP1C_SOURCE_OHM if k == 0 else STEP1C_OPEN_OHM,
                gap_height_m=s["h"],
                sheet_width_m=s["w"],
                drive_direction=(0.0, 1.0, 0.0),
                drive_voltage_v=1.0 + 0.0j,
                interior=True,
            )
        )
    t0 = time.perf_counter()
    res = run_lumped_sheet_port_case(
        problem,
        ports,
        specs,
        facet_tags=fx["facet_tags"],
        driven_port_id="P1",
        verbose=False,
    )
    wall = float(time.perf_counter() - t0)
    r1, r2 = res.responses["P1"], res.responses["P2"]
    z11 = complex(r1.voltage_v) / complex(r1.current_a)
    z21 = complex(r2.voltage_v) / complex(r1.current_a)
    return z11, z21 / z11, wall


def _solve_loop_modes(fx, c_f, *, target_hz, live):
    """One (G1) solve with torus 1's sheet (211) carrying ``C``; torus 2 open."""

    msh = fx["mesh"]
    V = fem.functionspace(msh, ("N1curl", 1))
    Q = fem.functionspace(msh, ("DG", 0))
    eps = fem.Function(Q, name="epsilon_r_air")
    eps.x.array[:] = 1.0 + 0.0j
    trial, test = ufl.TrialFunction(V), ufl.TestFunction(V)
    loop_sheet = [fx["sheets"][0]]
    forms = _capacitor_mass_forms(
        msh,
        fx["facet_tags"],
        loop_sheet,
        c_f,
        trial,
        test,
        omega_ref_rad_s=2.0 * np.pi * TARGET_FREQUENCY_HZ,
        live=live,
    )
    k0_sq = (2.0 * np.pi * float(target_hz) / C_0) ** 2
    t0 = time.perf_counter()
    values, n_converged, _n_bc = solve_general_mesh_modes(
        V,
        eps,
        target_k0_sq=k0_sq,
        sheet_mass_forms=forms,
        nev=NEV,
        comm=msh.comm,
    )
    wall = float(time.perf_counter() - t0)
    n_dofs = int(V.dofmap.index_map.size_global * V.dofmap.index_map_bs)
    freqs = [frequency_of(lam) for lam in values]
    return freqs, int(n_converged), n_dofs, wall, len(forms)


def _nearest(freqs, f_hz):
    return min(freqs, key=lambda f: abs(f.real - f_hz))


@complex_only
@pytest.mark.skipif(os.environ.get(STEP1C_ENV) != "1", reason=f"set {STEP1C_ENV}=1")
def test_step1c_lc_loop_closed_form():
    """(iii) the pencil-with-sheet-mass reproduces ``f_LC = 1/(2π√(LC))``.

    Anchors (asserted): (a) nearest eigenvalue to ``f_LC`` within 5 %;
    (b) at ``2·C`` the same mode moves by ``1/√2`` within 2 %.  Negative
    control (asserted): with the capacitor form removed (``live=()``) no
    eigenvalue has ``Re f`` in ``f_LC·[0.8, 1.2]``.
    """

    from tests.validation.test_port_reaction_impedance import MAJOR_RADIUS, MINOR_RADIUS
    from tests.validation.test_port_self_impedance_energy import grover_loop_inductance

    comm = MPI.COMM_WORLD
    t_all = time.perf_counter()
    fx = build_two_torus_hole_fixture()
    rank0 = comm.rank == 0
    if rank0:
        print(
            f"\n[TH-17 1c] two-torus hole fixture: {fx['n_cells']} cells, "
            f"mesh {fx['t_mesh_s']:.1f} s, ranks {comm.size}"
        )
        for s in fx["sheets"]:
            print(
                f"    sheet {s['tag']}: facets {s['facets']}, area {s['area']:.9e} m^2, "
                f"h {s['h']:.9e} m, w = A/h {s['w']:.9e} m, h/w {s['h'] / s['w']:.9e}"
            )

    l_fem = {}
    for f_hz in STEP1C_L_FREQUENCIES_HZ:
        z11, ratio21, wall = _driven_loop_impedance(fx, f_hz)
        omega = 2.0 * np.pi * f_hz
        l_fem[f_hz] = z11.imag / omega
        if rank0:
            print(
                f"[TH-17 1c] driven f = {f_hz / 1e6:.3f} MHz: Z_in = "
                f"{z11.real:+.9e} {z11.imag:+.9e}j Ohm, L_fem = Im Z_in/w = "
                f"{l_fem[f_hz]:.9e} H, Z21/Z11 = {ratio21.real:+.6e} "
                f"{ratio21.imag:+.6e}j (|.| = {abs(ratio21):.6e}), {wall:.1f} s"
            )
    l_use = l_fem[TARGET_FREQUENCY_HZ]
    assert l_use > 0.0, f"L_fem = {l_use!r} H at the target: not an inductive loop"
    l_grover = grover_loop_inductance(MAJOR_RADIUS, MINOR_RADIUS)
    omega_t = 2.0 * np.pi * TARGET_FREQUENCY_HZ
    c_f = 1.0 / (omega_t**2 * l_use)
    f_lc = 1.0 / (2.0 * np.pi * np.sqrt(l_use * c_f))
    if rank0:
        print(
            f"[TH-17 1c] Grover (printed, never gated) = {l_grover:.9e} H; "
            f"L_fem(target)/Grover = {l_use / l_grover:.6f}, "
            f"L_fem(10 MHz)/Grover = {l_fem[STEP1C_L_FREQUENCIES_HZ[0]] / l_grover:.6f}"
        )
        print(f"[TH-17 1c] C = 1/(w_t^2 L_fem) = {c_f:.9e} F; f_LC = {f_lc:.9e} Hz")

    freqs_c, n_c, n_dofs, wall_c, n_forms = _solve_loop_modes(fx, c_f, target_hz=f_lc, live=None)
    f_c = _nearest(freqs_c, f_lc)
    f_lc2 = f_lc / np.sqrt(2.0)
    freqs_2c, n_2c, _d, wall_2c, _n = _solve_loop_modes(fx, 2.0 * c_f, target_hz=f_lc2, live=None)
    f_2c = _nearest(freqs_2c, f_lc2)
    freqs_0, n_0, _d, wall_0, n_forms0 = _solve_loop_modes(fx, c_f, target_hz=f_lc, live=())
    lo, hi = (w * f_lc for w in STEP1C_CONTROL_WINDOW)
    in_window = [f for f in freqs_0 if lo <= f.real <= hi]

    err_a = abs(f_c.real - f_lc) / f_lc
    ratio_b = f_2c.real / f_c.real
    err_b = abs(ratio_b * np.sqrt(2.0) - 1.0)
    if rank0:
        print(
            f"[TH-17 1c] N1curl deg-1 dofs = {n_dofs}; sheet forms in pencil = "
            f"{n_forms}, control = {n_forms0}"
        )
        for label, freqs, n, wall in (
            ("C", freqs_c, n_c, wall_c),
            ("2C", freqs_2c, n_2c, wall_2c),
            ("no-C control", freqs_0, n_0, wall_0),
        ):
            print(f"[TH-17 1c] {label}: converged {n} of nev {NEV}, {wall:.1f} s")
            for f in freqs:
                print(f"    f = {f.real:.9e} {f.imag:+.9e}j Hz")
        print(
            f"[TH-17 1c] (a) nearest Re f = {f_c.real:.9e} Hz (Im {f_c.imag:+.6e}) vs "
            f"f_LC = {f_lc:.9e} Hz: rel err {err_a:.6e} (band {STEP1C_FLC_BAND})"
        )
        print(
            f"[TH-17 1c] (b) Re f(2C) = {f_2c.real:.9e} Hz; f(2C)/f(C) = {ratio_b:.9e} "
            f"vs 1/sqrt2 = {1.0 / np.sqrt(2.0):.9e}: rel err {err_b:.6e} "
            f"(band {STEP1C_SQRT2_BAND})"
        )
        above = sorted((f for f in freqs_0 if f.real > SCAN_PRINT_FLOOR_HZ), key=lambda f: f.real)
        first = f"{above[0].real:.9e} Hz" if above else "none (gradient cluster only)"
        print(
            f"[TH-17 1c] control: {len(in_window)} eigenvalue(s) in "
            f"f_LC*{STEP1C_CONTROL_WINDOW}; first converged Re f above "
            f"{SCAN_PRINT_FLOOR_HZ / 1e6:.0f} MHz = {first}"
        )
        print(f"[TH-17 1c] total {time.perf_counter() - t_all:.1f} s")

    assert n_forms == 1 and n_forms0 == 0, "the capacitor form did not reach/leave the pencil"
    assert n_c > 0 and n_2c > 0, "no eigenpair converged near f_LC or f_LC/sqrt2"
    assert (
        err_a <= STEP1C_FLC_BAND
    ), f"(a) nearest Re f {f_c.real:.6e} Hz vs f_LC {f_lc:.6e} Hz: {err_a:.3e} > 5 %"
    assert (
        err_b <= STEP1C_SQRT2_BAND
    ), f"(b) f(2C)/f(C) = {ratio_b:.6e} vs 1/sqrt2: {err_b:.3e} > 2 %"
    assert not in_window, (
        f"negative control: {len(in_window)} eigenvalue(s) in f_LC*[0.8,1.2] " f"with no capacitor"
    )


# ---------------------------------------------------------------------------
# Step 2 (§9 item 51, chain T5a) — the four-capacitor spectrum against the
# circuit layer.  Same fixture as 1b (ii) (hole mesh, uncorrected sheet widths,
# ``C_tuned`` on all four gaps).  The driven side runs `TH-15`'s four-port
# matched sweep at each f (every port a 50 Ω sheet, the exact power-wave route),
# converts S to the open-circuit Z, and diagonalises it with the four-port DFT:
# with the same capacitor on every port the structure resonates where
# ``det(Z + Z_C·I) = 0``, i.e. at a zero of ``z_k(f) + 1/(jωC)`` for one k —
# ``Im z_k = 1/(ωC)``.  ``z_k`` is the diagonal of ``F·Z·Fᴴ`` (never a class
# average, which would make (iii) true by construction).
#
# Linearisation: the pencil fixes the phantom's ε_c at 64 MHz, so every driven
# solve here carries ``σ·f/(64 MHz)`` (``ε_c = ε_r − jσ/(ωε₀)`` then matches the
# pencil's, `core/time_harmonic.py`).  The physical-σ zero and one eigen
# re-solve at ``ω_lin = 2π·Re f₁`` (the fixed-point shift) are printed only.
# ---------------------------------------------------------------------------

STEP2_ENV = "FEM_EM_TH17_STEP2"
#: The bracket grid: 21 points, 100–200 MHz in 5 MHz steps (the item).
STEP2_GRID_HZ = tuple(100.0e6 + 5.0e6 * i for i in range(21))
#: ≤ 6 secant steps per root (the item).
STEP2_MAX_SECANT = 6
#: 1b (ii)'s recorded k = 1 pair (`20260920T093810Z_TH-17.log:1769–1770`) —
#: printed beside the live pair, never gated (the live pair is the comparand).
STEP2_RECORDED_PAIR_HZ = (1.295936260e08, 1.296011466e08)
#: The eigen target for the live pair (1b (ii)'s 128 MHz target, which found it).
STEP2_EIGEN_TARGET_HZ = 128.0e6
#: Printed prediction for the 0.5·C control (rule (e)): ≈ √2 × 129.6 MHz, from
#: 1c's 0.70905 ratio on another fixture (`20260921T110359Z_TH-17.log:572`).
STEP2_HALF_C_PREDICTED_HZ = float(np.sqrt(2.0) * 129.6e6)


def _four_port_dft():
    n = LEG_COUNT
    idx = np.arange(n)
    return np.exp(-2j * np.pi * np.outer(idx, idx) / n) / np.sqrt(n)


def _driven_hole_z(fx, f_hz, *, sigma_s_per_m):
    """Open-circuit 4×4 Z of the hole fixture at ``f_hz`` (one matched sweep).

    Problem and specs exactly as ``test_th15_birdcage_pec_hole._hole_rung``
    builds them (every sheet a 50 Ω termination, P_k driven at 1 V in turn);
    only the phantom σ is the caller's (linearised or physical).  S is the
    exact matched-drive power-wave matrix, so ``Z = z0 (I+S)(I−S)⁻¹``.
    """
    from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicProblem
    from fem_em_solver.ports.definitions import PortDefinition
    from fem_em_solver.ports.lumped import LumpedSheetPortSpec
    from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
    from tests.validation.test_port_birdcage_four_port import (
        TERMINATED_PORT_IMPEDANCE_OHM,
    )
    from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

    problem = TimeHarmonicProblem(
        mesh=fx["mesh"],
        frequency_hz=float(f_hz),
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=fx["cell_tags"],
        material_map={
            PHANTOM_CELL_TAG: HomogeneousMaterial(
                sigma=float(sigma_s_per_m), epsilon_r=SALINE_EPSILON_R, mu_r=1.0
            ),
        },
        boundary_condition="pec_zero_tangential_a",
        facet_tags=fx["facet_tags"],
        pec_facet_tags=None,
    )
    port_defs, specs = [], []
    for s in fx["sheets"]:
        pid = f"P{s['tag'] - SHEET_IFACE}"
        port_defs.append(
            PortDefinition(
                port_id=pid,
                positive_tag=int(s["tag"]),
                negative_tag=CONDUCTOR_CELL_TAG,
                orientation="leg_gap_axial_plus_z",
                z0_ohm=REFERENCE_IMPEDANCE_OHM,
            )
        )
        specs.append(
            LumpedSheetPortSpec(
                port_id=pid,
                facet_tag=int(s["tag"]),
                port_impedance_ohm=TERMINATED_PORT_IMPEDANCE_OHM,
                gap_height_m=s["h"],
                sheet_width_m=s["w"],
                drive_direction=(0.0, 0.0, 1.0),
                drive_voltage_v=1.0 + 0.0j,
                interior=True,
            )
        )
    result = run_n_port_sparameter_sweep(
        problem,
        port_defs,
        lumped_sheet_ports=specs,
        lumped_sheet_facet_tags=fx["sheet_facet_tags"],
        degree=1,
    )
    s = np.asarray(result.s_matrix, dtype=np.complex128)
    eye = np.eye(s.shape[0])
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    return z0 * (eye + s) @ np.linalg.inv(eye - s)


def _modal_z(z):
    """``diag(F·Z·Fᴴ)`` and the off-diagonal residual (C4 diagonalisation check)."""
    F = _four_port_dft()
    m = F @ z @ F.conj().T
    off = m - np.diag(np.diag(m))
    return np.diag(m).copy(), float(np.linalg.norm(off) / np.linalg.norm(m))


def _loop_z(zk, f_hz, c_f):
    """``z_k + 1/(jωC)`` — zero at a resonance of the four-capacitor structure."""
    return complex(zk) + 1.0 / (1j * 2.0 * np.pi * float(f_hz) * float(c_f))


@complex_only
@pytest.mark.skipif(os.environ.get(STEP2_ENV) != "1", reason=f"set {STEP2_ENV}=1")
def test_step2_circuit_eigen_identity():
    """§9 item 51: circuit-layer zeros of the four-capacitor network vs the pencil.

    Asserted: (i) ``|f₁,circ − Re f_eig|/Re f_eig ≤ κ`` for both pair members
    (κ = `PORT-14`'s ``STEP2D_C_OVER_TERMINAL_64MHZ_P1``, imported); (ii) the
    eigen pair's split ≤ ``ADJACENT_SPREAD_BAND``; (iii) the z₁ / z₃ zeros'
    split ≤ the same band; negative control: at 0.5·C_tuned (same Z(f), no
    extra solves) the k = 1 zero misses ``Re f_eig`` by more than κ.
    """
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
    from tests.validation.test_port_circuit_layer_field import TUNING_IM_Z_RTOL
    from tests.validation.test_port_lumped_rlc_termination import (
        STEP2D_C_OVER_TERMINAL_64MHZ_P1 as KAPPA,
    )

    comm = MPI.COMM_WORLD
    t_all = time.perf_counter()
    c_f = c_tuned_farad()
    fx = build_hole_fixture()
    omega_64 = 2.0 * np.pi * TARGET_FREQUENCY_HZ
    counter = {"sweeps": 0}

    def say(msg):
        if comm.rank == 0:
            print(msg, flush=True)

    say(
        f"\n[TH-17 step2] cells = {fx['n_cells']}, mesh {fx['t_mesh_s']:.2f} s, "
        f"C_tuned = {c_f:.9e} F, kappa = {KAPPA:.6e}, "
        f"spread band = {ADJACENT_SPREAD_BAND}, rtol = {TUNING_IM_Z_RTOL:.1e}"
    )
    for s in fx["sheets"]:
        say(
            f"    sheet {s['tag']}: azimuth {s['azimuth_deg']:.3f} deg, facets {s['facets']}, "
            f"h {s['h']:.9e} m, w {s['w']:.9e} m, h/w {s['h'] / s['w']:.6f}"
        )

    # --- eigen side: the live pair at 1b (ii)'s linearisation -----------------
    values, n_conv, n_dofs, wall_e, _ = solve_modes(
        fx, c_f, omega_lin_rad_s=omega_64, target_hz=STEP2_EIGEN_TARGET_HZ
    )
    phys = sorted(
        (frequency_of(v) for v in values if frequency_of(v).real > SCAN_PRINT_FLOOR_HZ),
        key=lambda f: f.real,
    )
    say(
        f"[TH-17 step2] eigen (omega_lin = 2pi*64 MHz): dofs {n_dofs}, converged "
        f"{n_conv}, {wall_e:.2f} s"
    )
    for f in phys:
        say(f"    f = {f.real:.9e} {f.imag:+.9e}j Hz")
    assert len(phys) >= 2, "fewer than two physical eigenvalues at the 128 MHz target"
    pair = [phys[0].real, phys[1].real]
    f_eig = 0.5 * (pair[0] + pair[1])
    split_eig = abs(pair[1] - pair[0]) / f_eig
    say(
        f"[TH-17 step2] eigen pair {pair[0]:.9e} / {pair[1]:.9e} Hz (recorded 1b(ii) "
        f"{STEP2_RECORDED_PAIR_HZ[0]:.9e} / {STEP2_RECORDED_PAIR_HZ[1]:.9e}), "
        f"split {split_eig:.3e}"
    )

    # --- circuit side: the bracket grid (linearised σ) -----------------------
    def sigma_lin(f_hz):
        return SALINE_SIGMA * float(f_hz) / TARGET_FREQUENCY_HZ

    cache = {}

    def modal(f_hz, *, physical=False):
        key = (round(float(f_hz), 3), physical)
        if key not in cache:
            t0 = time.perf_counter()
            sig = SALINE_SIGMA if physical else sigma_lin(f_hz)
            z = _driven_hole_z(fx, f_hz, sigma_s_per_m=sig)
            zk, off = _modal_z(z)
            counter["sweeps"] += 1
            cache[key] = (zk, off, time.perf_counter() - t0)
        return cache[key]

    grid = np.array(STEP2_GRID_HZ)
    zk_grid = []
    for f in grid:
        zk, off, dt = modal(f)
        zk_grid.append(zk)
        say(
            f"    grid {f / 1e6:6.1f} MHz ({dt:.2f} s): off-diag {off:.2e}; Im z_k - 1/wC = "
            + ", ".join(f"{_loop_z(zk[k], f, c_f).imag:+.6e}" for k in range(LEG_COUNT))
        )
    zk_grid = np.array(zk_grid)

    def brackets(k, c_use):
        g = np.array([_loop_z(zk_grid[i, k], grid[i], c_use).imag for i in range(len(grid))])
        return [
            (i, g[i], g[i + 1]) for i in range(len(grid) - 1) if np.sign(g[i]) != np.sign(g[i + 1])
        ]

    def g_of(k, f, physical):
        return _loop_z(modal(f, physical=physical)[0][k], f, c_f)

    def secant(k, fa, fb, *, physical=False):
        ga, gb = g_of(k, fa, physical).imag, g_of(k, fb, physical).imag
        for _ in range(STEP2_MAX_SECANT):
            if gb == ga:
                break
            fc = fb - gb * (fb - fa) / (gb - ga)
            fa, ga = fb, gb
            fb = float(fc)
            zl = g_of(k, fb, physical)
            gb = zl.imag
            if abs(zl.imag) / abs(zl) <= TUNING_IM_Z_RTOL:
                break
        zl = g_of(k, fb, physical)
        rel = abs(zl.imag) / abs(zl)
        return float(fb), float(rel), bool(rel <= TUNING_IM_Z_RTOL)

    zeros = {}
    for k in range(LEG_COUNT):
        for i, _ga, _gb in brackets(k, c_f):
            f_root, rel, is_zero = secant(k, grid[i], grid[i + 1])
            say(
                f"[TH-17 step2] k = {k}: bracket [{grid[i] / 1e6:.0f}, {grid[i + 1] / 1e6:.0f}] "
                f"MHz -> f = {f_root:.9e} Hz, |Im|/|z| = {rel:.3e} "
                f"({'ZERO' if is_zero else 'pole/unconverged'})"
            )
            if is_zero:
                zeros.setdefault(k, []).append(f_root)

    def nearest_zero(k):
        cands = zeros.get(k, [])
        return min(cands, key=lambda f: abs(f - f_eig)) if cands else None

    f1, f3 = nearest_zero(1), nearest_zero(3)
    say(f"[TH-17 step2] z_1 zero {f1}, z_3 zero {f3} (Hz); all zeros {zeros}")
    assert f1 is not None and f3 is not None, "no classified zero for k = 1 or k = 3"
    f_circ = 0.5 * (f1 + f3)
    errs = [abs(f_circ - p) / p for p in pair]
    split_circ = abs(f1 - f3) / f_circ

    # --- negative control: 0.5·C_tuned on the same grid Z(f), no extra solves -
    half = 0.5 * c_f
    nc = [
        float(grid[i] - ga * (grid[i + 1] - grid[i]) / (gb - ga)) for i, ga, gb in brackets(1, half)
    ]
    nc_near = min(nc, key=lambda f: abs(f - f_eig)) if nc else None
    nc_err = abs(nc_near - f_eig) / f_eig if nc_near is not None else float("inf")
    say(
        f"[TH-17 step2] NEG 0.5*C: k=1 sign changes (grid-interpolated, zero or pole) "
        f"{nc}; nearest {nc_near}, |f - Re f_eig|/Re f_eig = {nc_err:.4e} "
        f"(= {nc_err / KAPPA:.2f} kappa); PREDICTED (printed, rule (e)) "
        f"~{STEP2_HALF_C_PREDICTED_HZ:.4e} Hz"
    )

    # --- printed only: physical-σ zero, fixed-point eigen re-solve ------------
    f1_phys, rel_p, ok_p = secant(1, f1, f1 * (1.0 + 1.0e-3), physical=True)
    say(
        f"[TH-17 step2] physical-sigma k=1 zero {f1_phys:.9e} Hz (|Im|/|z| {rel_p:.3e}, "
        f"{'ZERO' if ok_p else 'unconverged'}); shift vs linearised "
        f"{(f1_phys - f1) / f1:+.4e} (printed, not gated)"
    )
    values2, n2, _d, wall2, _ = solve_modes(
        fx,
        c_f,
        omega_lin_rad_s=2.0 * np.pi * pair[0],
        target_hz=STEP2_EIGEN_TARGET_HZ,
        z_ref_hz=pair[0],
    )
    phys2 = sorted(
        (frequency_of(v) for v in values2 if frequency_of(v).real > SCAN_PRINT_FLOOR_HZ),
        key=lambda f: f.real,
    )
    if phys2:
        say(
            f"[TH-17 step2] fixed-point re-solve (omega_lin = 2pi*{pair[0]:.6e}): "
            f"{wall2:.2f} s, lowest "
            + ", ".join(f"{f.real:.9e} {f.imag:+.9e}j" for f in phys2[:2])
            + f"; shift {(phys2[0].real - pair[0]) / pair[0]:+.4e} (printed, not gated)"
        )
    else:
        say(f"[TH-17 step2] fixed-point re-solve: no physical eigenvalue ({n2} converged)")

    total = time.perf_counter() - t_all
    say(
        f"[TH-17 step2] (i) f_circ = {f_circ:.9e} Hz (z1 {f1:.9e}, z3 {f3:.9e}) vs pair: "
        f"errs {errs[0]:.4e} / {errs[1]:.4e} (kappa {KAPPA:.4e}); (ii) eig split "
        f"{split_eig:.3e}; (iii) circ split {split_circ:.3e}; NEG {nc_err:.4e}; "
        f"{counter['sweeps']} sweeps, total {total:.1f} s"
    )

    assert split_eig <= ADJACENT_SPREAD_BAND, f"(ii) eigen split {split_eig:.3e}"
    assert split_circ <= ADJACENT_SPREAD_BAND, f"(iii) circuit split {split_circ:.3e}"
    assert nc_err > KAPPA, f"negative control: 0.5*C zero within kappa ({nc_err:.3e})"
    for p, e in zip(pair, errs):
        assert (
            e <= KAPPA
        ), f"(i) circuit zero {f_circ:.9e} vs eigen {p:.9e}: {e:.4e} > kappa {KAPPA:.4e}"


# ---------------------------------------------------------------------------
# Step 3 (§9 item 59, chain T5b): tune by the eigen route.  ``C_mode`` is the
# capacitance that zeroes the k = 1 loop reactance at 64 MHz,
# ``C_mode = 1/(ω₆₄·Im z₁(64 MHz))`` (z₁ ≡ z₃ by reciprocity: averaged, no C4
# claim).  At 64 MHz the linearised σ·f/64 MHz is the physical σ.  The pencil
# at ``C_mode`` (ω_lin = 2π·64 MHz, ``z_ref_hz=None`` so the capacitor mass is
# ∝ C) must put the lowest physical pair within κ of 64 MHz — a first reading
# (1c read circuit ↔ eigen at 2.053 % on the two-torus).  Mode ordering, the
# k = 0 / 2 circuit sign changes and the next eigenvalue are printed, never
# asserted.
# ---------------------------------------------------------------------------

STEP3_ENV = "FEM_EM_TH17_STEP3"
#: Circuit-side grid at ``C_mode``: 60–100 MHz in 5 MHz steps (the item).
STEP3_GRID_HZ = tuple(60.0e6 + 5.0e6 * i for i in range(9))
#: Printed prediction (item 59): L_eff = L₀/(1 − (f/523 MHz)²) fit to step 2.
STEP3_C_MODE_PREDICTED_F = 6.63e-11


@complex_only
@pytest.mark.skipif(os.environ.get(STEP3_ENV) != "1", reason=f"set {STEP3_ENV}=1")
def test_step3_eigen_route_tuning():
    """§9 item 59: the eigen pencil at ``C_mode`` lands the k = 1 pair on 64 MHz.

    Asserted: ``Im z₁(64 MHz) > 0``; (a) both members of the lowest physical
    pair within κ of 64 MHz (κ = ``STEP2D_C_OVER_TERMINAL_64MHZ_P1``, imported);
    (b) pair split ≤ ``ADJACENT_SPREAD_BAND``; negative control: at 2·C_mode
    (target 64/√2 MHz) the pair misses 64 MHz by more than κ.
    """
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
    from tests.validation.test_port_circuit_layer_field import TUNING_IM_Z_RTOL
    from tests.validation.test_port_lumped_rlc_termination import (
        STEP2D_C_OVER_TERMINAL_64MHZ_P1 as KAPPA,
    )

    comm = MPI.COMM_WORLD
    t_all = time.perf_counter()
    fx = build_hole_fixture()
    f64 = TARGET_FREQUENCY_HZ
    omega_64 = 2.0 * np.pi * f64
    counter = {"sweeps": 0}

    def say(msg):
        if comm.rank == 0:
            print(msg, flush=True)

    say(
        f"\n[TH-17 step3] cells = {fx['n_cells']}, mesh {fx['t_mesh_s']:.2f} s, "
        f"kappa = {KAPPA:.6e}, spread band = {ADJACENT_SPREAD_BAND}"
    )

    cache = {}

    def modal(f_hz):
        key = round(float(f_hz), 3)
        if key not in cache:
            t0 = time.perf_counter()
            z = _driven_hole_z(fx, f_hz, sigma_s_per_m=SALINE_SIGMA * float(f_hz) / f64)
            zk, off = _modal_z(z)
            counter["sweeps"] += 1
            cache[key] = (zk, off, time.perf_counter() - t0)
        return cache[key]

    def physical_sorted(vals):
        return sorted(
            (frequency_of(v) for v in vals if frequency_of(v).real > SCAN_PRINT_FLOOR_HZ),
            key=lambda f: f.real,
        )

    # --- C_mode from one sweep at 64 MHz -------------------------------------
    zk64, off64, dt64 = modal(f64)
    z1 = complex(0.5 * (zk64[1] + zk64[3]))
    say(
        f"[TH-17 step3] 64 MHz sweep ({dt64:.2f} s): off-diag {off64:.2e}; z_k = "
        + ", ".join(f"{complex(z):.6e}" for z in zk64)
        + f"; z1 (avg z1,z3) = {z1:.9e}, |z1-z3|/|z1| {abs(zk64[1] - zk64[3]) / abs(z1):.2e}"
    )
    assert z1.imag > 0.0, f"Im z1(64 MHz) = {z1.imag:.6e} <= 0: no C_mode"
    c_mode = 1.0 / (omega_64 * z1.imag)
    say(
        f"[TH-17 step3] C_mode = {c_mode:.9e} F (PREDICTED ~{STEP3_C_MODE_PREDICTED_F:.3e}, "
        f"pure 1/f^2 6.32e-11; printed); c_tuned_farad() = {c_tuned_farad():.6e} (not used)"
    )

    # --- eigen solve at C_mode -----------------------------------------------
    values, n_conv, n_dofs, wall_e, _ = solve_modes(
        fx, c_mode, omega_lin_rad_s=omega_64, target_hz=f64
    )
    phys = physical_sorted(values)
    say(
        f"[TH-17 step3] eigen at C_mode (target 64 MHz): dofs {n_dofs}, converged "
        f"{n_conv}, {wall_e:.2f} s"
    )
    for f in phys:
        say(f"    f = {f.real:.9e} {f.imag:+.9e}j Hz")
    assert len(phys) >= 2, "fewer than two physical eigenvalues near 64 MHz"
    pair = [phys[0].real, phys[1].real]
    f_pair = 0.5 * (pair[0] + pair[1])
    errs = [abs(p - f64) / f64 for p in pair]
    split = abs(pair[1] - pair[0]) / f_pair
    nxt = f"{phys[2].real:.9e} Hz" if len(phys) > 2 else "none converged"
    say(
        f"[TH-17 step3] pair {pair[0]:.9e} / {pair[1]:.9e} Hz: |Re f - 64|/64 = "
        f"{errs[0]:.4e} / {errs[1]:.4e} ({errs[0] / KAPPA:.3f} / {errs[1] / KAPPA:.3f} kappa); "
        f"split {split:.3e}; next above the pair: {nxt}"
    )

    # --- circuit side at C_mode (printed; ordering not asserted) -------------
    grid = np.array(STEP3_GRID_HZ)
    zk_grid = []
    for f in grid:
        zk, off, dt = modal(f)
        zk_grid.append(zk)
        say(
            f"    grid {f / 1e6:6.1f} MHz ({dt:.2f} s): off-diag {off:.2e}; Im z_k - 1/wC = "
            + ", ".join(f"{_loop_z(zk[k], f, c_mode).imag:+.6e}" for k in range(LEG_COUNT))
        )
    zk_grid = np.array(zk_grid)

    def g_grid(k, i):
        return _loop_z(zk_grid[i, k], grid[i], c_mode).imag

    def brackets(k):
        return [
            i for i in range(len(grid) - 1) if np.sign(g_grid(k, i)) != np.sign(g_grid(k, i + 1))
        ]

    def g_of(k, f):
        return _loop_z(modal(f)[0][k], f, c_mode)

    def secant(k, fa, fb):
        ga, gb = g_of(k, fa).imag, g_of(k, fb).imag
        for _ in range(STEP2_MAX_SECANT):
            if gb == ga:
                break
            fc = fb - gb * (fb - fa) / (gb - ga)
            fa, ga = fb, gb
            fb = float(fc)
            zl = g_of(k, fb)
            gb = zl.imag
            if abs(zl.imag) / abs(zl) <= TUNING_IM_Z_RTOL:
                break
        zl = g_of(k, fb)
        rel = abs(zl.imag) / abs(zl)
        return float(fb), float(rel), bool(rel <= TUNING_IM_Z_RTOL)

    zeros = {}
    for k in (1, 0, 2):
        for i in brackets(k):
            if k == 1:
                f_root, rel, is_zero = secant(k, grid[i], grid[i + 1])
                tag = "ZERO" if is_zero else "pole/unconverged"
            else:  # k = 0 / 2: grid-interpolated sign change only (no extra sweeps)
                ga, gb = g_grid(k, i), g_grid(k, i + 1)
                f_root = float(grid[i] - ga * (grid[i + 1] - grid[i]) / (gb - ga))
                rel, tag = float("nan"), "grid-interpolated, zero or pole"
            say(
                f"[TH-17 step3] circuit k = {k}: bracket [{grid[i] / 1e6:.0f}, "
                f"{grid[i + 1] / 1e6:.0f}] MHz -> f = {f_root:.9e} Hz, |Im|/|z| = {rel:.3e} "
                f"({tag})"
            )
            zeros.setdefault(k, []).append(f_root)
    say(f"[TH-17 step3] circuit sign changes at C_mode, 60-100 MHz (printed): {zeros}")

    # --- negative control: 2·C_mode, target 64/√2 MHz -------------------------
    values_n, _n, _d, wall_n, _ = solve_modes(
        fx, 2.0 * c_mode, omega_lin_rad_s=omega_64, target_hz=f64 / np.sqrt(2.0)
    )
    phys_n = physical_sorted(values_n)
    for f in phys_n:
        say(f"    NEG f = {f.real:.9e} {f.imag:+.9e}j Hz")
    assert len(phys_n) >= 2, "negative control: fewer than two physical eigenvalues"
    pair_n = [phys_n[0].real, phys_n[1].real]
    miss = min(abs(p - f64) / f64 for p in pair_n)
    ratio = (0.5 * (pair_n[0] + pair_n[1]) / f_pair) / (1.0 / np.sqrt(2.0))
    say(
        f"[TH-17 step3] NEG 2*C_mode: pair {pair_n[0]:.9e} / {pair_n[1]:.9e} Hz "
        f"({wall_n:.2f} s), miss {miss:.4e} (= {miss / KAPPA:.2f} kappa; PREDICTED ~27.4 "
        f"kappa, ceiling 27.6); (f_neg/f_pair)/(1/sqrt2) = {ratio:.6f} (printed)"
    )

    # --- printed only: fixed-point re-solve ----------------------------------
    values2, n2, _d, wall2, _ = solve_modes(
        fx, c_mode, omega_lin_rad_s=2.0 * np.pi * pair[0], target_hz=f64, z_ref_hz=pair[0]
    )
    phys2 = physical_sorted(values2)
    if phys2:
        say(
            f"[TH-17 step3] fixed-point re-solve (omega_lin = 2pi*{pair[0]:.6e}): "
            f"{wall2:.2f} s, lowest "
            + ", ".join(f"{f.real:.9e} {f.imag:+.9e}j" for f in phys2[:2])
            + f"; shift {(phys2[0].real - pair[0]) / pair[0]:+.4e} (printed, not gated)"
        )
    else:
        say(f"[TH-17 step3] fixed-point re-solve: no physical eigenvalue ({n2} converged)")

    total = time.perf_counter() - t_all
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2
    say(
        f"[TH-17 step3] SUMMARY Im z1 = {z1.imag:.6e} ohm, C_mode = {c_mode:.6e} F; (a) errs "
        f"{errs[0]:.4e} / {errs[1]:.4e} vs kappa {KAPPA:.4e}; (b) split {split:.3e}; "
        f"NEG miss {miss:.4e} ({miss / KAPPA:.2f} kappa); {counter['sweeps']} sweeps, "
        f"total {total:.1f} s, rank-0 peak RSS {peak:.2f} GiB"
    )

    assert split <= ADJACENT_SPREAD_BAND, f"(b) pair split {split:.3e}"
    assert miss > KAPPA, f"negative control: 2*C_mode pair within kappa ({miss:.3e})"
    for p, e in zip(pair, errs):
        assert e <= KAPPA, f"(a) eigen {p:.9e} Hz: |Re f - 64|/64 = {e:.4e} > kappa {KAPPA:.4e}"


# ---------------------------------------------------------------------------
# Gate (e) — §9 item 117: the σ → 0 control on the tuned pencil.  PEC coil,
# phantom at σ = 0 with ε_r kept: A (curl-curl), B_ε and B_sheet (Z_p = 1/(jωC)
# ⇒ jωμ₀/R real) are all real symmetric, so every λ of (G1) is real and the
# GNHEP read must return Im λ at round-off.  Asserted bound: TH-14's PEC
# control ``CONTROL_IM_RE_BOUND`` (imported).  Negative control (L): the lossy
# pencil's pair reads Im λ ≫ bound (measured 6.56e-4, TH-17-step3 log).
# ---------------------------------------------------------------------------

GATE_E_ENV = "FEM_EM_TH17_GATE_E"
#: σ for the printed-only Im(xᴴBx) comparison against σ = 0 (the item).
GATE_E_LOSSY_PROBE_SIGMA = 0.5
#: Step 3's (L) pair and Q, printed beside the live values (never gated).
GATE_E_STEP3_PAIR_HZ = (6.432991449e07, 6.433349959e07)
GATE_E_STEP3_Q_PREDICTED = (1523.7, 1528.8)


def _physical_pair_lambdas(values):
    """The two lowest λ with ``Re f > SCAN_PRINT_FLOOR_HZ``, by ``Re f``."""

    phys = sorted(
        (complex(v) for v in values if frequency_of(v).real > SCAN_PRINT_FLOOR_HZ),
        key=lambda v: frequency_of(v).real,
    )
    return phys


def _im_xbx_ratio(fx, c_f, *, omega_lin_rad_s, sigma_s_per_m, target_hz):
    """``|Im(xᴴBx)|/|xᴴBx|`` for one ``setRandom`` x on (G1)'s B (global dot)."""
    from fem_em_solver.core.cavity import general_mesh_pencil

    msh = fx["mesh"]
    V = fem.functionspace(msh, ("N1curl", 1))
    eps = permittivity_field(msh, fx["cell_tags"], omega_lin_rad_s, sigma_s_per_m=sigma_s_per_m)
    trial, test = ufl.TrialFunction(V), ufl.TestFunction(V)
    forms = _capacitor_mass_forms(
        msh,
        fx["sheet_facet_tags"],
        fx["sheets"],
        c_f,
        trial,
        test,
        omega_ref_rad_s=omega_lin_rad_s,
    )
    k0_sq = (2.0 * np.pi * float(target_hz) / C_0) ** 2
    A, B, _n = general_mesh_pencil(V, eps, bc_diagonal=1.0e4 * k0_sq, sheet_mass_forms=forms)
    x = B.createVecRight()
    rnd = PETSc.Random().create(comm=msh.comm)
    rnd.setFromOptions()
    rnd.setSeed(117)
    x.setRandom(rnd)
    y = B.createVecLeft()
    B.mult(x, y)
    xbx = complex(y.dot(x))  # VecDot(y, x) = xᴴ y = xᴴ B x, global
    for obj in (A, B, x, y, rnd):
        obj.destroy()
    return abs(xbx.imag) / abs(xbx), xbx


@complex_only
@pytest.mark.skipif(os.environ.get(GATE_E_ENV) != "1", reason=f"set {GATE_E_ENV}=1")
def test_gate_e_lossless_phantom_control():
    """§9 item 117 — `TH-17` gate (e): σ → 0 ⇒ Hermitian pencil ⇒ Im λ at round-off.

    Asserted: (H1) σ = 0 at ``C_mode`` and (H2) σ = 0 at ``c_tuned_farad()``,
    each member of the physical pair ``|Im λ|/Re λ ≤ CONTROL_IM_RE_BOUND``;
    negative control (L) σ = ``SALINE_SIGMA`` at ``C_mode``: both members
    ``> CONTROL_IM_RE_BOUND`` with ``Im f > 0``.  Printed: the Re f shifts,
    (L)'s Q, the (L)/(H1) ratio, and ``|Im(xᴴBx)|/|xᴴBx|`` at σ = 0 / 0.5.
    """
    from tests.validation.test_cavity_leontovich_q import CONTROL_IM_RE_BOUND

    comm = MPI.COMM_WORLD
    t_all = time.perf_counter()
    fx = build_hole_fixture()
    f64 = TARGET_FREQUENCY_HZ
    omega_64 = 2.0 * np.pi * f64

    def say(msg):
        if comm.rank == 0:
            print(msg, flush=True)

    say(
        f"\n[TH-17 gate e] cells = {fx['n_cells']}, mesh {fx['t_mesh_s']:.2f} s, "
        f"bound {CONTROL_IM_RE_BOUND:.0e}"
    )

    # --- C_mode from one sweep at 64 MHz, as step 3 ---------------------------
    t0 = time.perf_counter()
    z = _driven_hole_z(fx, f64, sigma_s_per_m=SALINE_SIGMA * float(f64) / f64)
    zk64, off64 = _modal_z(z)
    z1 = complex(0.5 * (zk64[1] + zk64[3]))
    say(
        f"[TH-17 gate e] 64 MHz sweep ({time.perf_counter() - t0:.2f} s): off-diag "
        f"{off64:.2e}; z1 = {z1:.9e}"
    )
    assert z1.imag > 0.0, f"Im z1(64 MHz) = {z1.imag:.6e} <= 0: no C_mode"
    c_mode = 1.0 / (omega_64 * z1.imag)
    c_tuned = c_tuned_farad()
    say(f"[TH-17 gate e] C_mode = {c_mode:.9e} F; c_tuned_farad() = {c_tuned:.9e} F")

    cases = {
        "L": dict(c_f=c_mode, sigma=SALINE_SIGMA, target_hz=f64),
        "H1": dict(c_f=c_mode, sigma=0.0, target_hz=f64),
        "H2": dict(c_f=c_tuned, sigma=0.0, target_hz=STEP2_EIGEN_TARGET_HZ),
    }
    pairs = {}
    for name, cs in cases.items():
        values, n_conv, n_dofs, wall, _ = solve_modes(
            fx,
            cs["c_f"],
            omega_lin_rad_s=omega_64,
            target_hz=cs["target_hz"],
            z_ref_hz=None,
            sigma_s_per_m=cs["sigma"],
        )
        say(
            f"[TH-17 gate e] ({name}) sigma = {cs['sigma']}, C = {cs['c_f']:.9e} F, target "
            f"{cs['target_hz']:.6e} Hz: dofs {n_dofs}, converged {n_conv}, {wall:.2f} s"
        )
        for v in sorted((complex(v) for v in values), key=lambda v: v.real):
            f = frequency_of(v)
            say(
                f"    ({name}) lambda = {v.real:.9e} {v.imag:+.9e}j; f = {f.real:.9e} "
                f"{f.imag:+.9e}j Hz; |Im l|/Re l = {abs(v.imag) / abs(v.real):.3e}"
            )
        phys = _physical_pair_lambdas(values)
        assert len(phys) >= 2, f"({name}) fewer than two physical eigenvalues"
        pair = phys[:2]
        fr = [frequency_of(v) for v in pair]
        ratios = [abs(v.imag) / v.real for v in pair]
        split = abs(fr[1].real - fr[0].real) / (0.5 * (fr[0].real + fr[1].real))
        say(
            f"[TH-17 gate e] ({name}) pair {fr[0].real:.9e} / {fr[1].real:.9e} Hz, split "
            f"{split:.3e}; Im f {fr[0].imag:+.3e} / {fr[1].imag:+.3e} Hz; |Im l|/Re l "
            f"{ratios[0]:.3e} / {ratios[1]:.3e}"
        )
        pairs[name] = (pair, fr, ratios)

    # --- printed records ------------------------------------------------------
    _l_pair, l_f, l_r = pairs["L"]
    q_l = [f.real / (2.0 * abs(f.imag)) for f in l_f]
    say(
        f"[TH-17 gate e] (L) Q = {q_l[0]:.1f} / {q_l[1]:.1f} (PREDICTED "
        f"{GATE_E_STEP3_Q_PREDICTED[0]} / {GATE_E_STEP3_Q_PREDICTED[1]}); pair vs step 3 "
        f"{GATE_E_STEP3_PAIR_HZ[0]:.9e} / {GATE_E_STEP3_PAIR_HZ[1]:.9e}: rel "
        + " / ".join(f"{(f.real - r) / r:+.3e}" for f, r in zip(l_f, GATE_E_STEP3_PAIR_HZ))
    )
    h1_f, h1_r = pairs["H1"][1], pairs["H1"][2]
    h2_f, h2_r = pairs["H2"][1], pairs["H2"][2]
    say(
        "[TH-17 gate e] Re f shift (H1) - (L), relative: "
        + " / ".join(f"{(h.real - l.real) / l.real:+.4e}" for h, l in zip(h1_f, l_f))
    )
    say(
        "[TH-17 gate e] (H2) vs STEP2_RECORDED_PAIR_HZ (lossy, 128 MHz tuned; no prior): "
        + " / ".join(f"{(h.real - r) / r:+.4e}" for h, r in zip(h2_f, STEP2_RECORDED_PAIR_HZ))
    )
    say(
        "[TH-17 gate e] (L)/(H1) |Im l|/Re l ratio (PREDICTED >= 1e9, printed): "
        + " / ".join(f"{lr / max(hr, 1e-300):.3e}" for lr, hr in zip(l_r, h1_r))
    )

    rb = {}
    for sig in (0.0, GATE_E_LOSSY_PROBE_SIGMA):
        rb[sig], xbx = _im_xbx_ratio(
            fx, c_mode, omega_lin_rad_s=omega_64, sigma_s_per_m=sig, target_hz=f64
        )
        say(
            f"[TH-17 gate e] |Im(x^H B x)|/|x^H B x| at sigma = {sig}: {rb[sig]:.3e} "
            f"(x^H B x = {xbx:.6e})"
        )

    total = time.perf_counter() - t_all
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2
    say(
        f"[TH-17 gate e] SUMMARY (H1) {h1_r[0]:.3e} / {h1_r[1]:.3e}; (H2) {h2_r[0]:.3e} / "
        f"{h2_r[1]:.3e}; (L) {l_r[0]:.3e} / {l_r[1]:.3e} (bound {CONTROL_IM_RE_BOUND:.0e}); "
        f"Im(xBx) sigma 0 {rb[0.0]:.3e}, sigma {GATE_E_LOSSY_PROBE_SIGMA} "
        f"{rb[GATE_E_LOSSY_PROBE_SIGMA]:.3e}; total {total:.1f} s, rank-0 peak RSS "
        f"{peak:.2f} GiB"
    )

    # --- asserts: negative control first, then the gate -----------------------
    for f, r in zip(l_f, l_r):
        assert r > CONTROL_IM_RE_BOUND, f"(L) |Im l|/Re l = {r:.3e} <= bound: sigma never reached"
        assert f.imag > 0.0, f"(L) Im f = {f.imag:.3e} <= 0"
    for name in ("H1", "H2"):
        for f, r in zip(pairs[name][1], pairs[name][2]):
            assert (
                r <= CONTROL_IM_RE_BOUND
            ), f"({name}) f = {f.real:.9e} Hz: |Im l|/Re l = {r:.3e} > {CONTROL_IM_RE_BOUND:.0e}"


# ---------------------------------------------------------------------------
# Gate (c) and (d′) — §9 item 131 (chain T5c, row amended by rule (j)).  (c) is
# step 2's identity (i) on the k = 2 class: the circuit zero of
# ``Im(z₂ + 1/(jωC))`` (z_k = diag(F·Z·Fᴴ) from the S-derived network Z,
# linearised σ·f/64 MHz exactly as steps 2 / 3) against the pencil's lowest
# physical eigenvalue above the k = 1 pair, checked a singlet, at C_mode (M)
# and at ``c_tuned_farad()`` (T); k = 0 only if a classified zero falls in the
# leg's grid.  (d′) at 2·C_mode: ``|(f̄_pair(2C)/f̄_pair(C))·√2 − 1| ≤ κ``.
# Negative controls: (N1) the k = 2 circuit zero against the k = 1 pair, (N2)
# the C_mode k = 2 circuit zero against the pencil's k = 2 at 2·C_mode.  The
# ``modal`` / ``secant`` closures of steps 2 / 3 are copied into
# ``_CircuitSide`` below; steps 2 / 3 are untouched.
# ---------------------------------------------------------------------------

GATE_C_ENV = "FEM_EM_TH17_GATE_C"
#: Step 2's printed k = 2 circuit zero at ``c_tuned_farad()``
#: (`20260927T094905Z_TH-17-step2.log`) — leg (T) must reproduce it to 1e-9.
GATE_C_STEP2_K2_ZERO_HZ = 1.581275645e08
GATE_C_STEP2_K2_REPRO_RTOL = 1.0e-9


class _CircuitSide:
    """Steps 2 / 3's ``modal`` (cached, linearised σ) and ``secant`` closures."""

    def __init__(self, fx, c_f, grid_hz, say, tag):
        self.fx, self.c_f, self.say, self.tag = fx, float(c_f), say, tag
        self.grid = np.array(grid_hz)
        self.cache = {}
        self.sweeps = 0
        self.zk_grid = None

    def modal(self, f_hz):
        key = round(float(f_hz), 3)
        if key not in self.cache:
            t0 = time.perf_counter()
            sig = SALINE_SIGMA * float(f_hz) / TARGET_FREQUENCY_HZ
            z = _driven_hole_z(self.fx, f_hz, sigma_s_per_m=sig)
            zk, off = _modal_z(z)
            self.sweeps += 1
            self.cache[key] = (zk, off, time.perf_counter() - t0)
        return self.cache[key]

    def scan(self):
        rows = []
        for f in self.grid:
            zk, off, dt = self.modal(f)
            rows.append(zk)
            self.say(
                f"    [{self.tag}] grid {f / 1e6:6.1f} MHz ({dt:.2f} s): off-diag "
                f"{off:.2e}; Im z_k - 1/wC = "
                + ", ".join(f"{_loop_z(zk[k], f, self.c_f).imag:+.6e}" for k in range(LEG_COUNT))
            )
        self.zk_grid = np.array(rows)

    def brackets(self, k):
        g = [
            _loop_z(self.zk_grid[i, k], self.grid[i], self.c_f).imag for i in range(len(self.grid))
        ]
        return [i for i in range(len(self.grid) - 1) if np.sign(g[i]) != np.sign(g[i + 1])]

    def g_of(self, k, f):
        return _loop_z(self.modal(f)[0][k], f, self.c_f)

    def secant(self, k, fa, fb, rtol):
        ga, gb = self.g_of(k, fa).imag, self.g_of(k, fb).imag
        for _ in range(STEP2_MAX_SECANT):
            if gb == ga:
                break
            fc = fb - gb * (fb - fa) / (gb - ga)
            fa, ga = fb, gb
            fb = float(fc)
            zl = self.g_of(k, fb)
            gb = zl.imag
            if abs(zl.imag) / abs(zl) <= rtol:
                break
        zl = self.g_of(k, fb)
        rel = abs(zl.imag) / abs(zl)
        return float(fb), float(rel), bool(rel <= rtol)

    def zeros(self, k, rtol):
        out = []
        for i in self.brackets(k):
            f_root, rel, ok = self.secant(k, self.grid[i], self.grid[i + 1], rtol)
            self.say(
                f"[{self.tag}] circuit k = {k}: bracket [{self.grid[i] / 1e6:.0f}, "
                f"{self.grid[i + 1] / 1e6:.0f}] MHz -> f = {f_root:.9e} Hz, |Im|/|z| = "
                f"{rel:.3e} ({'ZERO' if ok else 'pole/unconverged'})"
            )
            if ok:
                out.append(f_root)
        return out


def _physical_freqs(values):
    return sorted(
        (frequency_of(v) for v in values if frequency_of(v).real > SCAN_PRINT_FLOOR_HZ),
        key=lambda f: f.real,
    )


def _gate_c_leg(fx, c_f, *, grid_hz, eigen_target_hz, say, tag):
    """One leg of gate (c): returns a dict of the readings (asserts by the caller)."""
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
    from tests.validation.test_port_circuit_layer_field import TUNING_IM_Z_RTOL

    omega_64 = 2.0 * np.pi * TARGET_FREQUENCY_HZ
    # --- eigen side: the k = 1 pair and the k = 2 candidate above it ----------
    values, n_conv, n_dofs, wall, _ = solve_modes(
        fx, c_f, omega_lin_rad_s=omega_64, target_hz=eigen_target_hz
    )
    phys = _physical_freqs(values)
    say(
        f"[{tag}] eigen at C = {c_f:.9e} F (target {eigen_target_hz:.6e} Hz): dofs {n_dofs}, "
        f"converged {n_conv}, {wall:.2f} s"
    )
    for f in phys:
        say(f"    [{tag}] f = {f.real:.9e} {f.imag:+.9e}j Hz")
    assert len(phys) >= 3, f"[{tag}] fewer than three physical eigenvalues (no k = 2 candidate)"
    pair = [phys[0].real, phys[1].real]
    f_pair = 0.5 * (pair[0] + pair[1])
    f2 = phys[2].real
    # --- singlet check: shift-invert at Re f₂; the second-nearest is > band away
    values_s, n_s, _d, wall_s, _ = solve_modes(fx, c_f, omega_lin_rad_s=omega_64, target_hz=f2)
    phys_s = _physical_freqs(values_s)
    by_dist = sorted(phys_s, key=lambda f: abs(f.real - f2))
    say(
        f"[{tag}] singlet solve at target Re f2 = {f2:.9e} Hz ({wall_s:.2f} s, converged "
        f"{n_s}): " + ", ".join(f"{f.real:.9e}" for f in phys_s)
    )
    assert by_dist, f"[{tag}] singlet solve found no physical eigenvalue"
    self_rel = abs(by_dist[0].real - f2) / f2
    nb = by_dist[1].real if len(by_dist) > 1 else None
    nb_gap = abs(nb - f2) / f2 if nb is not None else float("inf")
    say(
        f"[{tag}] k = 2 candidate {f2:.9e} Hz: re-found at {by_dist[0].real:.9e} "
        f"(rel {self_rel:.2e}); nearest other {nb} Hz, gap {nb_gap:.4e} "
        f"(band {ADJACENT_SPREAD_BAND})"
    )
    # --- circuit side --------------------------------------------------------
    cs = _CircuitSide(fx, c_f, grid_hz, say, tag)
    cs.scan()
    z1 = cs.zeros(1, TUNING_IM_Z_RTOL)
    z2 = cs.zeros(2, TUNING_IM_Z_RTOL)
    z0 = cs.zeros(0, TUNING_IM_Z_RTOL)
    say(f"[{tag}] k = 0 sign changes in grid at indices {cs.brackets(0)}; classified zeros {z0}")
    # --- printed only: fixed-point re-solve at ω_lin = 2π·Re f₂ ---------------
    values_f, _nf, _d, wall_f, _ = solve_modes(
        fx, c_f, omega_lin_rad_s=2.0 * np.pi * f2, target_hz=f2, z_ref_hz=f2
    )
    phys_f = _physical_freqs(values_f)
    if phys_f:
        ff = min(phys_f, key=lambda f: abs(f.real - f2))
        say(
            f"[{tag}] fixed-point re-solve (omega_lin = 2pi*{f2:.6e}): {wall_f:.2f} s, "
            f"nearest {ff.real:.9e} {ff.imag:+.9e}j; shift {(ff.real - f2) / f2:+.4e} "
            f"(printed, not gated)"
        )
    else:
        say(f"[{tag}] fixed-point re-solve: no physical eigenvalue")
    return dict(
        phys=phys,
        pair=pair,
        f_pair=f_pair,
        f2=f2,
        self_rel=self_rel,
        nb_gap=nb_gap,
        z1=z1,
        z2=z2,
        z0=z0,
        sweeps=cs.sweeps,
        all_eig=[f.real for f in phys + phys_s],
    )


def _assert_gate_c_leg(r, kappa, tag, say):
    """Asserted on both legs: singlet, N1, f₁ < f₂ both sides, k = 0 if any, (c)."""
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND

    f2 = r["f2"]
    assert r["z2"], f"[{tag}] no classified k = 2 circuit zero in the grid"
    assert r["z1"], f"[{tag}] no classified k = 1 circuit zero in the grid"
    fc2 = min(r["z2"], key=lambda f: abs(f - f2))
    fc1 = min(r["z1"], key=lambda f: abs(f - r["f_pair"]))
    err_c = abs(fc2 - f2) / f2
    n1 = abs(fc2 - r["f_pair"]) / r["f_pair"]
    say(
        f"[{tag}] (c) circuit k = 2 zero {fc2:.9e} Hz vs Re f2 {f2:.9e} Hz: {err_c:.4e} "
        f"(= {err_c / kappa:.3f} kappa); f1 circ {fc1:.9e} vs f2 circ {fc2:.9e}; eigen pair "
        f"{r['pair'][0]:.9e} / {r['pair'][1]:.9e} vs {f2:.9e}; (N1) k2 zero vs k1 pair "
        f"{n1:.4e} (= {n1 / kappa:.2f} kappa)"
    )
    assert (
        r["self_rel"] <= 1.0e-6
    ), f"[{tag}] singlet solve did not re-find f2 ({r['self_rel']:.2e})"
    assert r["nb_gap"] > ADJACENT_SPREAD_BAND, f"[{tag}] k = 2 not a singlet ({r['nb_gap']:.3e})"
    assert n1 > kappa, f"[{tag}] (N1) class swap within kappa ({n1:.3e})"
    assert max(r["pair"]) < f2, f"[{tag}] eigen side: f1 >= f2"
    assert fc1 < fc2, f"[{tag}] circuit side: f1 {fc1:.9e} >= f2 {fc2:.9e}"
    for f0 in r["z0"]:  # pre-registered: the identity on a k = 0 zero at the same band
        fe = min(r["all_eig"], key=lambda f: abs(f - f0))
        e0 = abs(f0 - fe) / fe
        say(f"[{tag}] k = 0 zero {f0:.9e} vs nearest eigen {fe:.9e}: {e0:.4e}")
        assert e0 <= kappa, f"[{tag}] (c) k = 0: {e0:.4e} > kappa {kappa:.4e}"
    assert err_c <= kappa, f"[{tag}] (c) k = 2: {err_c:.4e} > kappa {kappa:.4e}"
    return fc2, err_c, n1


@complex_only
@pytest.mark.skipif(os.environ.get(GATE_C_ENV) != "1", reason=f"set {GATE_C_ENV}=1")
def test_gate_c_k2_identity_c_mode():
    """§9 item 131, leg (M): gate (c) on the k = 2 class at C_mode, (d′), N1, N2."""
    from tests.validation.test_port_lumped_rlc_termination import (
        STEP2D_C_OVER_TERMINAL_64MHZ_P1 as KAPPA,
    )

    comm = MPI.COMM_WORLD
    t_all = time.perf_counter()
    fx = build_hole_fixture()
    f64 = TARGET_FREQUENCY_HZ
    omega_64 = 2.0 * np.pi * f64
    tag = "TH-17 gate c M"

    def say(msg):
        if comm.rank == 0:
            print(msg, flush=True)

    say(f"\n[{tag}] cells = {fx['n_cells']}, mesh {fx['t_mesh_s']:.2f} s, kappa {KAPPA:.6e}")
    t0 = time.perf_counter()
    z = _driven_hole_z(fx, f64, sigma_s_per_m=SALINE_SIGMA)
    zk64, off64 = _modal_z(z)
    z1 = complex(0.5 * (zk64[1] + zk64[3]))
    assert z1.imag > 0.0, f"Im z1(64 MHz) = {z1.imag:.6e} <= 0: no C_mode"
    c_mode = 1.0 / (omega_64 * z1.imag)
    say(
        f"[{tag}] 64 MHz sweep ({time.perf_counter() - t0:.2f} s): off-diag {off64:.2e}; "
        f"z1 = {z1:.9e}; C_mode = {c_mode:.9e} F"
    )

    r = _gate_c_leg(fx, c_mode, grid_hz=STEP3_GRID_HZ, eigen_target_hz=f64, say=say, tag=tag)
    fc2, err_c, n1 = _assert_gate_c_leg(r, KAPPA, tag, say)

    # --- (d′) and (N2): the pencil at 2·C_mode, target 64/√2 MHz ---------------
    values_n, _n, _d, wall_n, _ = solve_modes(
        fx, 2.0 * c_mode, omega_lin_rad_s=omega_64, target_hz=f64 / np.sqrt(2.0)
    )
    phys_n = _physical_freqs(values_n)
    for f in phys_n:
        say(f"    [{tag}] 2C f = {f.real:.9e} {f.imag:+.9e}j Hz")
    assert len(phys_n) >= 3, "2*C_mode: fewer than three physical eigenvalues"
    f_pair_n = 0.5 * (phys_n[0].real + phys_n[1].real)
    f2_n = phys_n[2].real
    d_prime = abs((f_pair_n / r["f_pair"]) * np.sqrt(2.0) - 1.0)
    n2 = abs(fc2 - f2_n) / f2_n
    k2_ratio = (f2_n / r["f2"]) * np.sqrt(2.0)
    total = time.perf_counter() - t_all
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2
    say(
        f"[{tag}] SUMMARY (c-M) {err_c:.4e} ({err_c / KAPPA:.3f} kappa); (d') {d_prime:.4e} "
        f"({d_prime / KAPPA:.3f} kappa; pair(2C) {f_pair_n:.9e}, pair(C) {r['f_pair']:.9e}); "
        f"(N1) {n1:.4e} ({n1 / KAPPA:.2f} kappa); (N2) circuit k2 {fc2:.9e} vs pencil k2 at "
        f"2C {f2_n:.9e}: {n2:.4e} ({n2 / KAPPA:.2f} kappa); k = 2 sqrt2 ratio {k2_ratio:.6f} "
        f"(printed); {r['sweeps'] + 1} sweeps, 2C solve {wall_n:.2f} s, total {total:.1f} s, "
        f"rank-0 peak RSS {peak:.2f} GiB"
    )
    assert n2 > KAPPA, f"(N2) C mismatch within kappa ({n2:.3e})"
    assert d_prime <= KAPPA, f"(d') {d_prime:.4e} > kappa {KAPPA:.4e}"


@complex_only
@pytest.mark.skipif(os.environ.get(GATE_C_ENV) != "1", reason=f"set {GATE_C_ENV}=1")
def test_gate_c_k2_identity_c_tuned():
    """§9 item 131, leg (T): gate (c) on the k = 2 class at ``c_tuned_farad()``, N1."""
    from tests.validation.test_port_lumped_rlc_termination import (
        STEP2D_C_OVER_TERMINAL_64MHZ_P1 as KAPPA,
    )

    comm = MPI.COMM_WORLD
    t_all = time.perf_counter()
    fx = build_hole_fixture()
    c_f = c_tuned_farad()
    tag = "TH-17 gate c T"

    def say(msg):
        if comm.rank == 0:
            print(msg, flush=True)

    say(
        f"\n[{tag}] cells = {fx['n_cells']}, mesh {fx['t_mesh_s']:.2f} s, kappa {KAPPA:.6e}, "
        f"C_tuned {c_f:.9e} F"
    )
    r = _gate_c_leg(
        fx, c_f, grid_hz=STEP2_GRID_HZ, eigen_target_hz=STEP2_EIGEN_TARGET_HZ, say=say, tag=tag
    )
    fc2, err_c, n1 = _assert_gate_c_leg(r, KAPPA, tag, say)
    repro = abs(fc2 - GATE_C_STEP2_K2_ZERO_HZ) / GATE_C_STEP2_K2_ZERO_HZ
    total = time.perf_counter() - t_all
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0**2
    say(
        f"[{tag}] SUMMARY (c-T) {err_c:.4e} ({err_c / KAPPA:.3f} kappa); k2 zero {fc2:.9e} vs "
        f"step 2 {GATE_C_STEP2_K2_ZERO_HZ:.9e}: {repro:.2e}; (N1) {n1:.4e} "
        f"({n1 / KAPPA:.2f} kappa); {r['sweeps']} sweeps, total {total:.1f} s, rank-0 peak "
        f"RSS {peak:.2f} GiB"
    )
    assert (
        repro <= GATE_C_STEP2_K2_REPRO_RTOL
    ), f"k = 2 zero does not reproduce step 2 ({repro:.2e})"
