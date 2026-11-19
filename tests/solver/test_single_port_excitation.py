"""Tests for single-port excitation (chunk E3/D2) and `PORT-4`.

**`PORT-4` — with unequal terminations, drive and read-back stay adjoint, so Y
is symmetric.**  Replaces the monkeypatched finiteness placeholder (the
``PORT-0`` heuristic; known-issues entry 6).  On the lumped-sheet route a
port's termination is ``LumpedSheetPortSpec.port_impedance_ohm`` (sheet
R = Z_p·w/h).  For any termination set, Y[i, j] = I_i(drive j)/V_src is
symmetric: the sheet term is complex-symmetric for any R, the source is
(V_src/(R_j h_j))∫ĥ·v and the read-back is (1/(R_i h_i))∫E·ĥ — the same
functional (``ports/lumped.py``).  At all-50 Ω a read-back scaled by the wrong
port's R is invisible; P2 = 200 Ω / P4 = 25 Ω exposes it.

**Premise (asserted first).**  Cells = ``STEP1_CELL_RECORD`` (116 085) and
the 50 Ω sweep reproduces ``S_10MHZ_EPS0_RECORD`` at ``RECORD_RTOL``.

**Band (asserted).**  ‖Y − Yᵀ‖_F/‖Y‖_F ≤ 1e-12 on the unequal sweep
(predicted ≈ 1e-14: the fixture's 50 Ω reciprocity floor is 2.152e-14).

**Negative control (asserted).**  A termination-blind read-back
Y_b = diag(Z_p/z0)·Y has per-pair asymmetry |1 − Z_p,j/Z_p,i| (normalised
by the larger entry): 0.875 on (P2, P4), 0.75 on (P2, P1/P3), 0.5 on
(P4, P1/P3).  Asserted: the worst pair ≥ 0.5 (exactly 0.875 if the anchor
holds; ceiling 1).

Records (printed, never gated): the sweep's own ``s_matrix`` reciprocity
under these terminations (it is not an S-matrix there), and
‖Y − (s_to_z(S₅₀, 50) + diag Z_p)⁻¹‖/‖Y‖ against `PORT-14`'s 200 Ω floor
7.249519e-04.  ``z_matrix`` is never read as an open-circuit Z.

Scope: one fixture, 10 MHz, degree 1, lumped route, resistive terminations.

Run (complex build, ``-n 2`` only)::

    scripts/testing/run_and_log.sh PORT-4-W1 "docker compose exec -T fem-em-solver \\
      bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && \\
      PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 \\
      mpiexec -n 2 python3 -m pytest tests/environment \\
      tests/solver/test_single_port_excitation.py -v -s --tb=short'"
"""

from __future__ import annotations

import dataclasses
import time

import numpy as np
import pytest
from dolfinx import fem
from dolfinx.mesh import create_unit_cube, meshtags
from mpi4py import MPI

from fem_em_solver.core import HomogeneousMaterial, TimeHarmonicFields, TimeHarmonicProblem
from fem_em_solver.ports import (
    PortDefinition,
    run_n_port_sparameter_sweep,
    run_single_port_excitation_case,
)
from fem_em_solver.ports.circuit import s_to_z
from tests.complex_mode import complex_only
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_circuit_layer_field import (
    RECORD_RANK_WIDTH,
    RECORD_RTOL,
    S_10MHZ_EPS0_RECORD,
    _residual,
)
from tests.validation.test_port_lumped_rlc_termination import (
    REDUCTION_FLOOR_F_SMALL,
    STEP1_CELL_RECORD,
)

UNEQUAL_TERMINATIONS_OHM = {"P2": 200.0, "P4": 25.0}
Y_SYMMETRY_BAND = 1.0e-12
PREDICTED_Y_ASYMMETRY = 1.0e-14  # the fixture's 50 Ω reciprocity floor is 2.152e-14
BLIND_CONTROL_FLOOR = 0.5  # asserted
BLIND_CONTROL_PREDICTED = 0.875  # |1 − 25/200| on (P2, P4); printed


def _build_test_problem():
    comm = MPI.COMM_WORLD
    mesh = create_unit_cube(comm, 2, 1, 1)

    tdim = mesh.topology.dim
    n_local_cells = mesh.topology.index_map(tdim).size_local
    cell_indices = np.arange(n_local_cells, dtype=np.int32)

    # Four terminal regions: P1+ (11), P1- (12), P2+ (21), P2- (22)
    tags = np.array([11, 12, 21, 22], dtype=np.int32)
    cell_values = tags[cell_indices % len(tags)]
    cell_tags = meshtags(mesh, tdim, cell_indices, cell_values)

    problem = TimeHarmonicProblem(
        mesh=mesh,
        frequency_hz=127.74e6,
        material=HomogeneousMaterial(sigma=0.2, epsilon_r=5.0, mu_r=1.0),
        cell_tags=cell_tags,
        facet_tags=None,
    )
    return problem


@complex_only
def test_unequal_terminations_keep_y_symmetric():
    comm = MPI.COMM_WORLD
    assert (
        comm.size == RECORD_RANK_WIDTH
    ), f"records are -n {RECORD_RANK_WIDTH} records; run at that width"

    b = build_four_port_sweep(build_only=True)
    ncells = int(b["cells"])
    specs = list(b["specs"])
    ids = [s.port_id for s in specs]
    assert ids == ["P1", "P2", "P3", "P4"], ids
    z0 = 50.0
    for s in specs:
        assert complex(s.port_impedance_ohm) == z0
    unequal_specs = [
        (
            dataclasses.replace(s, port_impedance_ohm=UNEQUAL_TERMINATIONS_OHM[s.port_id])
            if s.port_id in UNEQUAL_TERMINATIONS_OHM
            else s
        )
        for s in specs
    ]
    zp = np.array([complex(s.port_impedance_ohm) for s in unequal_specs])

    def _sweep(port_specs):
        comm.Barrier()
        t0 = time.perf_counter()
        r = run_n_port_sparameter_sweep(
            b["problem"],
            b["port_defs"],
            lumped_sheet_ports=port_specs,
            lumped_sheet_facet_tags=b["facet_tags"],
            degree=b["degree"],
        )
        comm.Barrier()
        y = np.array(
            [
                [
                    r.excitation_results[pj.port_id].responses[pi.port_id].current_a
                    / complex(pj.drive_voltage_v)
                    for pj in port_specs
                ]
                for pi in port_specs
            ],
            dtype=np.complex128,
        )
        return np.asarray(r.s_matrix, dtype=np.complex128), y, time.perf_counter() - t0

    def _asym(m):
        return float(np.linalg.norm(m - m.T) / np.linalg.norm(m))

    s50, y50, t50 = _sweep(specs)

    # --- premise (asserted first) -----------------------------------------
    premise = _residual(s50, np.asarray(S_10MHZ_EPS0_RECORD))
    if comm.rank == 0:
        print(
            f"\n[PORT-4] premise: {ncells} cells (record {STEP1_CELL_RECORD}), "
            f"mesh {b['mesh_time']:.2f} s; 50 Ω sweep {t50:.2f} s at -n {comm.size}; "
            f"‖S − S_10MHZ_EPS0_RECORD‖_F/‖rec‖_F = {premise:.3e} "
            f"(ASSERTED ≤ {RECORD_RTOL:.0e})",
            flush=True,
        )
    assert ncells == STEP1_CELL_RECORD
    assert premise <= RECORD_RTOL

    # --- the unequal-termination sweep on the same mesh --------------------
    s_u, y_u, t_u = _sweep(unequal_specs)
    anchor = _asym(y_u)

    y_b = np.diag(zp / z0) @ y_u
    pair_asym = {
        (ids[i], ids[j]): float(abs(y_b[i, j] - y_b[j, i]) / max(abs(y_b[i, j]), abs(y_b[j, i])))
        for i in range(4)
        for j in range(i + 1, 4)
    }
    worst_pair = max(pair_asym, key=pair_asym.get)
    control = pair_asym[worst_pair]

    y_pred = np.linalg.inv(s_to_z(s50, z0) + np.diag(zp))

    if comm.rank == 0:
        print(
            f"[PORT-4] unequal sweep (P2 200 Ω, P4 25 Ω) {t_u:.2f} s\n"
            f"[PORT-4] anchor ‖Y − Yᵀ‖_F/‖Y‖_F = {anchor:.3e} (ASSERTED ≤ "
            f"{Y_SYMMETRY_BAND:.0e}; predicted ≈ {PREDICTED_Y_ASYMMETRY:.0e}); "
            f"50 Ω sweep {_asym(y50):.3e} (printed)",
            flush=True,
        )
        for i in range(4):
            for j in range(i + 1, 4):
                print(
                    f"    Y[{ids[i]},{ids[j]}]/Y[{ids[j]},{ids[i]}] = "
                    f"{y_u[i, j] / y_u[j, i]:.12f}   blind asym "
                    f"{pair_asym[(ids[i], ids[j])]:.6f}",
                    flush=True,
                )
        print(
            f"[PORT-4] negative control: termination-blind worst pair {worst_pair} = "
            f"{control:.6f} (ASSERTED ≥ {BLIND_CONTROL_FLOOR}; predicted "
            f"{BLIND_CONTROL_PREDICTED}, ceiling 1); ‖Y_b − Y_bᵀ‖_F/‖Y_b‖_F = "
            f"{_asym(y_b):.4f} (printed)\n"
            f"[PORT-4] record (a): sweep s_matrix asymmetry under these terminations "
            f"= {_asym(s_u):.3e} (not an S-matrix there; printed)\n"
            f"[PORT-4] record (b): ‖Y − (s_to_z(S₅₀, 50) + diag Z_p)⁻¹‖_F/‖Y‖_F = "
            f"{np.linalg.norm(y_u - y_pred) / np.linalg.norm(y_u):.6e}, sign-flipped "
            f"{np.linalg.norm(y_u + y_pred) / np.linalg.norm(y_u):.6e} (printed; "
            f"`PORT-14` 200 Ω floor {REDUCTION_FLOOR_F_SMALL['R = 200 Ohm']:.6e})",
            flush=True,
        )

    assert control >= BLIND_CONTROL_FLOOR
    assert anchor <= Y_SYMMETRY_BAND


def test_single_port_excitation_rejects_missing_required_tags():
    problem = _build_test_problem()
    ports = [PortDefinition(port_id="P1", positive_tag=11, negative_tag=99, orientation="cw")]

    with pytest.raises(ValueError, match=r"missing required port tags"):
        run_single_port_excitation_case(
            problem,
            ports,
            driven_port_id="P1",
        )


def test_single_port_excitation_rejects_driven_index_mismatch(monkeypatch):
    problem = _build_test_problem()

    from fem_em_solver.ports import excitation as excitation_module

    class DummyTimeHarmonicSolver:
        def __init__(self, problem, degree=1):
            self.problem = problem

        def solve(self, **_kwargs):
            dg_vec = fem.functionspace(self.problem.mesh, ("DG", 1, (3,)))
            dg0 = fem.functionspace(self.problem.mesh, ("DG", 0))
            e_real = fem.Function(dg_vec, name="E_real_dummy")
            e_imag = fem.Function(dg_vec, name="E_imag_dummy")
            sigma = fem.Function(dg0, name="sigma_dummy")
            epsilon_r = fem.Function(dg0, name="epsilon_dummy")
            return TimeHarmonicFields(
                e_real=e_real,
                e_imag=e_imag,
                frequency_hz=self.problem.frequency_hz,
                sigma_field=sigma,
                epsilon_r_field=epsilon_r,
            )

    monkeypatch.setattr(excitation_module, "TimeHarmonicSolver", DummyTimeHarmonicSolver)

    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw"),
    ]

    with pytest.raises(ValueError, match=r"driven_port_id/index mismatch"):
        run_single_port_excitation_case(
            problem,
            ports,
            driven_port_id="P1",
            driven_port_index=1,
        )


def test_single_port_excitation_rejects_invalid_passive_termination_map(monkeypatch):
    problem = _build_test_problem()

    from fem_em_solver.ports import excitation as excitation_module

    class DummyTimeHarmonicSolver:
        def __init__(self, problem, degree=1):
            self.problem = problem

        def solve(self, **_kwargs):
            dg_vec = fem.functionspace(self.problem.mesh, ("DG", 1, (3,)))
            dg0 = fem.functionspace(self.problem.mesh, ("DG", 0))
            e_real = fem.Function(dg_vec, name="E_real_dummy")
            e_imag = fem.Function(dg_vec, name="E_imag_dummy")
            sigma = fem.Function(dg0, name="sigma_dummy")
            epsilon_r = fem.Function(dg0, name="epsilon_dummy")
            return TimeHarmonicFields(
                e_real=e_real,
                e_imag=e_imag,
                frequency_hz=self.problem.frequency_hz,
                sigma_field=sigma,
                epsilon_r_field=epsilon_r,
            )

    monkeypatch.setattr(excitation_module, "TimeHarmonicSolver", DummyTimeHarmonicSolver)

    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw"),
    ]

    with pytest.raises(ValueError, match=r"must not include driven port"):
        run_single_port_excitation_case(
            problem,
            ports,
            driven_port_id="P1",
            passive_port_terminations_ohm={"P1": 25.0},
        )

    with pytest.raises(ValueError, match=r"passive termination must be finite and positive"):
        run_single_port_excitation_case(
            problem,
            ports,
            driven_port_id="P1",
            passive_port_terminations_ohm={"P2": 0.0},
        )
