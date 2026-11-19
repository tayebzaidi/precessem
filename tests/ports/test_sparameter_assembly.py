"""Tests for N-port sweep and S-parameter matrix assembly (chunk E4/D3)."""

from __future__ import annotations

import numpy as np
import pytest

from fem_em_solver.ports import (
    PortDefinition,
    run_n_port_sparameter_sweep,
    summarize_sparameter_sanity,
)


class _DummyComm:
    rank = 0


class _DummyMesh:
    comm = _DummyComm()


class _DummyProblem:
    def __init__(self, frequency_hz: float = 127.74e6):
        self.frequency_hz = frequency_hz
        self.mesh = _DummyMesh()


# --- `PORT-5` (§9 item 98): the assembler and the sanity metrics on typed S. ---
# Every gated S this repo has solved is reciprocal, so the `[recv, drive]`
# convention of `_assemble_sparameter_matrix` had never met S != S^T.  An ideal
# circulator does.  S0 is typed here and never read back from the code.

Z0_OHM = 50.0
THETA_LINE = 0.7
# Incident amplitude of the fake's driven port (non-zero: the assembler guards it).
A0 = 0.3 - 0.2j

S_LINE = np.array(
    [[0.0, np.exp(-1j * THETA_LINE)], [np.exp(-1j * THETA_LINE), 0.0]], dtype=np.complex128
)
S_CIRC = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=np.complex128)
S_FLAT = 0.7 * np.ones((2, 2), dtype=np.complex128)


def _s_rand() -> np.ndarray:
    rng = np.random.default_rng(20260930)
    m = rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))
    return 0.8 * m / np.linalg.svd(m, compute_uv=False).max()


S_RAND = _s_rand()
ASSEMBLY_CASES = {"line": S_LINE, "circulator": S_CIRC, "random": S_RAND}
ASSEMBLY_ATOL = 1.0e-13
METRIC_ATOL = 1.0e-12


def _ports(n: int) -> list[PortDefinition]:
    return [
        PortDefinition(
            port_id=f"P{k + 1}",
            positive_tag=10 * (k + 1) + 1,
            negative_tag=10 * (k + 1) + 2,
            orientation="cw",
            z0_ohm=Z0_OHM,
        )
        for k in range(n)
    ]


def _sweep_known_s(monkeypatch, s0: np.ndarray):
    """Run the sweep over a fake whose port states are exactly those of ``s0``.

    Drive j: a = A0·e_j, b = A0·s0[:, j], V = √z0 (a + b), I = (a − b)/√z0.
    """
    from fem_em_solver.ports import sparameters as sparam_module

    ports = _ports(s0.shape[0])

    def _fake_single_port_case(problem, ports, *, driven_port_id, **_kwargs):
        from fem_em_solver.ports import (
            PortSolveContext,
            PortVoltageCurrentEstimate,
            SinglePortExcitationResult,
        )

        drive_idx = [p.port_id for p in ports].index(driven_port_id)
        responses = {}
        solve_context = {}
        for idx, port in enumerate(ports):
            is_driven = idx == drive_idx
            a = A0 if is_driven else 0.0 + 0.0j
            b = A0 * complex(s0[idx, drive_idx])
            responses[port.port_id] = PortVoltageCurrentEstimate(
                port_id=port.port_id,
                voltage_v=np.sqrt(port.z0_ohm) * (a + b),
                current_a=(a - b) / np.sqrt(port.z0_ohm),
                is_driven=is_driven,
                termination_ohm=port.z0_ohm,
            )
            ring = abs(idx - drive_idx)
            solve_context[port.port_id] = PortSolveContext(
                port_id=port.port_id,
                port_index=idx,
                driven_port_id=driven_port_id,
                driven_port_index=drive_idx,
                is_driven=is_driven,
                wrapped_ring_distance=min(ring, len(ports) - ring),
                coupling_factor=1.0 if is_driven else 0.0,
                termination_ohm=port.z0_ohm,
            )
        return SinglePortExcitationResult(
            driven_port_id=driven_port_id,
            frequency_hz=problem.frequency_hz,
            responses=responses,
            solve_context=solve_context,
        )

    monkeypatch.setattr(sparam_module, "run_placeholder_port_coupling_case", _fake_single_port_case)
    return run_n_port_sparameter_sweep(_DummyProblem(), ports, drive_voltage_v=1.0 + 0.0j)


@pytest.mark.parametrize("name", list(ASSEMBLY_CASES))
def test_n_port_sweep_assembles_the_known_s_matrix(monkeypatch, name):
    """Anchor A: the sweep returns the S whose port states the fake supplied."""
    s0 = ASSEMBLY_CASES[name]
    result = _sweep_known_s(monkeypatch, s0)
    err = float(np.max(np.abs(result.s_matrix - s0)))
    err_t = float(np.max(np.abs(result.s_matrix - s0.T)))
    print(
        f"\nPORT5_ASSEMBLY {name}: n={s0.shape[0]} max|S - S0|={err:.3e} "
        f"max|S - S0^T|={err_t:.6f} (negative-control comparand)"
    )
    assert result.port_ids == tuple(f"P{k + 1}" for k in range(s0.shape[0]))
    assert result.s_matrix.shape == s0.shape
    assert err <= ASSEMBLY_ATOL


def test_transposed_and_hermitian_comparands_are_seen_only_off_reciprocal_data(monkeypatch):
    """(n1) S_circ^T misses by exactly 1; (n2) S_line^H misses by 2 sin θ."""
    circ = _sweep_known_s(monkeypatch, S_CIRC).s_matrix
    line = _sweep_known_s(monkeypatch, S_LINE).s_matrix
    n1_circ = float(np.max(np.abs(circ - S_CIRC.T)))
    n1_line = float(np.max(np.abs(line - S_LINE.T)))
    n2_line = float(np.max(np.abs(line - S_LINE.conj().T)))
    reciprocity_line = summarize_sparameter_sanity(line).reciprocity_max_abs_delta
    print(
        f"\nPORT5_CONTROL n1 transposed comparand: circulator {n1_circ:.12f} (exactly 1), "
        f"line {n1_line:.3e} (reciprocal data cannot see a transposed assembly)"
    )
    print(
        f"PORT5_CONTROL n2 Hermitian comparand on the line: {n2_line:.12f} "
        f"(2 sin theta = {2.0 * np.sin(THETA_LINE):.12f}) where the reciprocity metric reads "
        f"{reciprocity_line:.3e}"
    )
    assert abs(n1_circ - 1.0) <= ASSEMBLY_ATOL
    assert n1_line <= ASSEMBLY_ATOL
    assert abs(n2_line - 2.0 * np.sin(THETA_LINE)) <= METRIC_ATOL
    assert abs(n2_line - 1.288435) <= 5e-7
    assert reciprocity_line <= ASSEMBLY_ATOL


def _metrics(report):
    return (
        report.reciprocity_max_abs_delta,
        report.reciprocity_max_rel_delta,
        report.passivity_max_sigma,
        report.passivity_max_column_power_sum,
    )


def _kinds(report):
    """Which warnings fired: (reciprocity abs, reciprocity rel, sigma, column power)."""
    text = report.warnings
    return (
        any("reciprocity abs" in w for w in text),
        any("reciprocity rel" in w for w in text),
        any("passivity sigma" in w for w in text),
        any("passivity column power" in w for w in text),
    )


@pytest.mark.parametrize(
    "name, s_matrix, want, kinds",
    [
        ("line", S_LINE, (0.0, 0.0, 1.0, 1.0), (False, False, False, False)),
        ("circulator", S_CIRC, (1.0, 1.0, 1.0, 1.0), (True, True, False, False)),
        ("line x 1.02", 1.02 * S_LINE, (0.0, 0.0, 1.02, 1.02**2), (False, False, False, False)),
        ("line x 1.03", 1.03 * S_LINE, (0.0, 0.0, 1.03, 1.03**2), (False, False, False, True)),
        ("line x 1.06", 1.06 * S_LINE, (0.0, 0.0, 1.06, 1.06**2), (False, False, True, True)),
        ("flat 0.7", S_FLAT, (0.0, 0.0, 1.4, 0.98), (False, False, True, False)),
    ],
)
def test_sanity_metrics_equal_their_closed_forms(name, s_matrix, want, kinds):
    """Anchor B: (abs Δ, rel Δ, σ_max, max column power) and which warnings fire."""
    report = summarize_sparameter_sanity(s_matrix)
    got = _metrics(report)
    print(
        f"\nPORT5_METRICS {name}: got={tuple(f'{g:.15g}' for g in got)} want={want} "
        f"warnings={_kinds(report)}"
    )
    assert max(abs(g - w) for g, w in zip(got, want)) <= METRIC_ATOL
    assert _kinds(report) == kinds
    assert len(report.warnings) == sum(kinds)


def test_sweep_sanity_report_equals_the_direct_call_on_the_circulator(monkeypatch):
    result = _sweep_known_s(monkeypatch, S_CIRC)
    direct = summarize_sparameter_sanity(S_CIRC)
    diff = max(abs(a - b) for a, b in zip(_metrics(result.sanity_report), _metrics(direct)))
    print(f"\nPORT5_SWEEP_REPORT circulator: max metric difference sweep vs direct = {diff:.3e}")
    assert diff <= METRIC_ATOL
    assert _kinds(result.sanity_report) == _kinds(direct) == (True, True, False, False)


def test_n_port_sweep_rejects_zero_incident_drive(monkeypatch):
    from fem_em_solver.ports import sparameters as sparam_module

    ports = [
        PortDefinition(port_id="P1", positive_tag=11, negative_tag=12, orientation="cw"),
        PortDefinition(port_id="P2", positive_tag=21, negative_tag=22, orientation="cw"),
    ]

    def _fake_zero_incident(problem, ports, *, driven_port_id, **_kwargs):
        from fem_em_solver.ports import (
            PortSolveContext,
            PortVoltageCurrentEstimate,
            SinglePortExcitationResult,
        )

        responses = {}
        solve_context = {}
        drive_idx = [p.port_id for p in ports].index(driven_port_id)
        for idx, port in enumerate(ports):
            is_driven = port.port_id == driven_port_id
            if is_driven:
                voltage = 1.0 + 0.0j
                current = -(1.0 / port.z0_ohm) + 0.0j  # forces a = (V + Z0*I)/(2*sqrt(Z0)) = 0
            else:
                voltage = 0.0 + 0.0j
                current = 0.0 + 0.0j

            responses[port.port_id] = PortVoltageCurrentEstimate(
                port_id=port.port_id,
                voltage_v=voltage,
                current_a=current,
                is_driven=is_driven,
                termination_ohm=port.z0_ohm,
            )
            solve_context[port.port_id] = PortSolveContext(
                port_id=port.port_id,
                port_index=idx,
                driven_port_id=driven_port_id,
                driven_port_index=drive_idx,
                is_driven=is_driven,
                wrapped_ring_distance=min(abs(idx - drive_idx), len(ports) - abs(idx - drive_idx)),
                coupling_factor=1.0 if is_driven else 0.0,
                termination_ohm=port.z0_ohm,
            )

        return SinglePortExcitationResult(
            driven_port_id=driven_port_id,
            frequency_hz=problem.frequency_hz,
            responses=responses,
            solve_context=solve_context,
        )

    monkeypatch.setattr(sparam_module, "run_placeholder_port_coupling_case", _fake_zero_incident)

    with pytest.raises(ValueError, match="incident wave.*is zero"):
        run_n_port_sparameter_sweep(_DummyProblem(), ports)


def test_sparameter_sanity_metrics_report_low_reciprocity_delta_for_symmetric_matrix():
    s_matrix = np.array(
        [
            [0.05 + 0.01j, 0.10 - 0.02j, 0.08 + 0.00j],
            [0.10 - 0.02j, 0.04 + 0.00j, 0.07 + 0.01j],
            [0.08 + 0.00j, 0.07 + 0.01j, 0.03 - 0.01j],
        ],
        dtype=np.complex128,
    )

    report = summarize_sparameter_sanity(s_matrix)

    assert report.reciprocity_max_abs_delta <= 1e-12
    assert report.reciprocity_max_rel_delta <= 1e-12
    assert report.passivity_max_sigma <= 1.0
    assert report.passivity_max_column_power_sum <= 1.0
    assert report.warnings == ()


def test_sparameter_sanity_metrics_emit_warnings_for_non_reciprocal_or_non_passive_matrix():
    s_matrix = np.array(
        [
            [0.10 + 0.00j, 0.70 + 0.00j],
            [0.20 + 0.00j, 1.20 + 0.00j],
        ],
        dtype=np.complex128,
    )

    report = summarize_sparameter_sanity(s_matrix)

    assert report.reciprocity_max_abs_delta > 0.05
    assert report.passivity_max_sigma > 1.05
    assert report.passivity_max_column_power_sum > 1.05
    assert any("reciprocity" in warning for warning in report.warnings)
    assert any("passivity" in warning for warning in report.warnings)
