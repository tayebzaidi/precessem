"""Validation checks for phantom |E|/|B| consistency diagnostics (chunk C5)."""

from __future__ import annotations

import functools
import math

import pytest

from fem_em_solver.post import compute_field_consistency_diagnostics
from tests.complex_mode import complex_only


def test_field_consistency_metrics_are_finite_and_warning_oriented():
    """Consistency indicators should be finite and warning-free for balanced stats."""
    e_stats = {"min": 1.0, "max": 3.0, "mean": 2.0}
    b_stats = {"min": 0.5, "max": 1.5, "mean": 1.0}

    diagnostics = compute_field_consistency_diagnostics(
        e_stats, b_stats, reference_speed_m_per_s=1.0
    )

    assert math.isfinite(diagnostics["e_to_b_mean_ratio"])
    assert math.isfinite(diagnostics["e_to_b_max_ratio"])
    assert math.isfinite(diagnostics["e_span_ratio"])
    assert math.isfinite(diagnostics["b_span_ratio"])
    assert math.isfinite(diagnostics["mean_balance_rel_diff"])

    assert diagnostics["e_to_b_mean_ratio"] == 2.0
    assert diagnostics["e_to_b_max_ratio"] == 2.0
    assert diagnostics["e_span_ratio"] == (3.0 - 1.0) / 3.0
    assert diagnostics["b_span_ratio"] == (1.5 - 0.5) / 1.5
    assert diagnostics["mean_balance_rel_diff"] == 0.5
    assert diagnostics["warnings"] == []


def test_field_consistency_metrics_emit_actionable_warnings_for_extreme_imbalance():
    """High span and severe |E|/|B| mean imbalance should emit warning strings."""
    e_stats = {"min": 0.0, "max": 10.0, "mean": 10.0}
    b_stats = {"min": 0.0, "max": 1.0, "mean": 0.01}

    diagnostics = compute_field_consistency_diagnostics(
        e_stats, b_stats, reference_speed_m_per_s=1.0
    )

    assert diagnostics["e_span_ratio"] == 1.0
    assert diagnostics["b_span_ratio"] == 1.0
    assert (
        diagnostics["mean_balance_rel_diff"]
        > diagnostics["thresholds"]["max_mean_balance_rel_diff_warn"]
    )

    warnings = " | ".join(diagnostics["warnings"])
    assert "span ratio" in warnings
    assert "strongly imbalanced" in warnings


# ---------------------------------------------------------------------------
# `POST-2` (2026-09-29, §9 item 92): the diagnostics on a field with a closed
# form.  `TH-6`'s lossy plane wave E = ẑ e^{−jkx} has ∇×E = jk e^{−jkx} ŷ, so
# |E|/|B| = ω/|k| at every point and the ratio of the means is exactly ω/|k|;
# both magnitudes fall as e^{−αx}, so each span ratio is 1 − e^{−αΔx_c} with
# Δx_c the (allreduced) centroid x-extent.  Imports are by name only — never a
# ``test_*`` function, or pytest would re-collect `TH-6` here.
# ---------------------------------------------------------------------------

POST2_BAND = 0.05  # `TH-6`'s MVP 5 %, PROJECT_PLAN §7 `POST-2` Done-when.


@functools.lru_cache(maxsize=None)
def _post2_read_diagnostics(n: int):
    import numpy as np
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    from fem_em_solver.post import compute_tagged_vector_magnitude_stats
    from fem_em_solver.post.faraday import magnetic_flux_density_from_e
    from tests.validation.test_lossy_plane_wave import (
        FREQUENCY_HZ,
        SIGMA,
        _solve_plane_wave,
    )

    comm = MPI.COMM_WORLD
    _rel, alpha_fit, beta_fit, ncells, msh, fields = _solve_plane_wave(n, SIGMA, return_fields=True)
    omega = 2.0 * np.pi * FREQUENCY_HZ
    b_fn = magnetic_flux_density_from_e(fields.e_complex, omega)

    tdim = msh.topology.dim
    cmap = msh.topology.index_map(tdim)
    n_local = cmap.size_local + cmap.num_ghosts
    cells = np.arange(n_local, dtype=np.int32)
    tags = dmesh.meshtags(msh, tdim, cells, np.ones(n_local, dtype=np.int32))

    e_stats = compute_tagged_vector_magnitude_stats(fields.e_complex, tags, 1, comm=comm)
    b_stats = compute_tagged_vector_magnitude_stats(b_fn, tags, 1, comm=comm)
    diagnostics = compute_field_consistency_diagnostics(e_stats, b_stats)

    # Centroid x-extent over owned cells, MIN/MAX-allreduced (rank-local otherwise).
    owned = np.arange(cmap.size_local, dtype=np.int32)
    x = msh.geometry.x
    dofmap = msh.geometry.dofmap
    cx = np.mean(x[dofmap[owned], 0], axis=1) if owned.size else np.array([])
    x_min = comm.allreduce(float(cx.min()) if cx.size else np.inf, op=MPI.MIN)
    x_max = comm.allreduce(float(cx.max()) if cx.size else -np.inf, op=MPI.MAX)
    return diagnostics, e_stats, b_stats, x_max - x_min, alpha_fit, beta_fit, ncells


@complex_only
@pytest.mark.integration
def test_post2_diagnostics_read_the_lossy_plane_wave_closed_form():
    """`e_to_b_mean_ratio` = ω/|k| and spans = 1 − e^{−αΔx_c} within 5 %."""
    import numpy as np
    from mpi4py import MPI

    from fem_em_solver.utils.constants import EPSILON_0, MU_0
    from tests.validation.test_lossy_plane_wave import (
        EPSILON_R,
        FREQUENCY_HZ,
        SIGMA,
        _analytic_alpha_beta,
        _wavenumber,
    )

    rank0 = MPI.COMM_WORLD.rank == 0
    omega = 2.0 * np.pi * FREQUENCY_HZ
    k = _wavenumber(SIGMA)
    alpha_cf, beta_cf = _analytic_alpha_beta(SIGMA)
    alpha = -k.imag
    ratio_cf = omega / abs(k)
    c0 = 1.0 / np.sqrt(MU_0 * EPSILON_0)
    ratio_sigma_blind = c0 / np.sqrt(EPSILON_R)  # negative-control comparand

    results = {}
    for n in (24, 12):
        diag, e_stats, b_stats, dx_c, a_fit, b_fit, ncells = _post2_read_diagnostics(n)
        span_cf = 1.0 - np.exp(-alpha * dx_c)
        mean_miss = diag["e_to_b_mean_ratio"] / ratio_cf - 1.0
        max_miss = diag["e_to_b_max_ratio"] / ratio_cf - 1.0
        e_span_miss = diag["e_span_ratio"] / span_cf - 1.0
        b_span_miss = diag["b_span_ratio"] / span_cf - 1.0
        balance_cf = 1.0 - abs(k) / omega
        control_miss = diag["e_to_b_mean_ratio"] / ratio_sigma_blind - 1.0
        results[n] = (mean_miss, e_span_miss, b_span_miss, control_miss)
        if rank0:
            print(
                f"\nPOST-2 n={n}: cells={ncells} samples E/B={e_stats['count']}/{b_stats['count']}"
                f"\n  premise fit: alpha={a_fit:.5f} (cf {alpha_cf:.5f}, {abs(a_fit/alpha_cf-1)*100:.3f} %),"
                f" beta={b_fit:.5f} (cf {beta_cf:.5f}, {abs(b_fit/beta_cf-1)*100:.3f} %)"
                f"\n  mean |E|={e_stats['mean']:.6e} V/m, mean |B|={b_stats['mean']:.6e} T"
                f"\n  e_to_b_mean_ratio={diag['e_to_b_mean_ratio']:.6e} vs omega/|k|={ratio_cf:.6e}"
                f"  miss {mean_miss*100:+.4f} %"
                f"\n  e_to_b_max_ratio={diag['e_to_b_max_ratio']:.6e}  miss {max_miss*100:+.4f} % (printed)"
                f"\n  dx_c={dx_c:.6f}  span_cf={span_cf:.6f}"
                f"\n  e_span_ratio={diag['e_span_ratio']:.6f} miss {e_span_miss*100:+.4f} %"
                f"\n  b_span_ratio={diag['b_span_ratio']:.6f} miss {b_span_miss*100:+.4f} %"
                f"\n  mean_balance_rel_diff={diag['mean_balance_rel_diff']!r} with no reference speed"
                f" (the old cross-unit reading was 1-|k|/omega={balance_cf:.12f})"
                f"\n  'strongly imbalanced' warning present: "
                f"{any('strongly imbalanced' in w for w in diag['warnings'])}"
                f"\n  negative control vs c/sqrt(78)={ratio_sigma_blind:.6e}: miss {control_miss*100:+.4f} %"
            )

    mean_miss, e_span_miss, b_span_miss, control_miss = results[24]
    assert abs(mean_miss) < POST2_BAND, f"mean-ratio miss {mean_miss:.4%} at 24^3"
    assert abs(e_span_miss) < POST2_BAND, f"E span miss {e_span_miss:.4%} at 24^3"
    assert abs(b_span_miss) < POST2_BAND, f"B span miss {b_span_miss:.4%} at 24^3"
    # Comparand swap (predicted 21.21 %): the band resolves σ.
    assert abs(control_miss) > POST2_BAND, f"sigma-blind comparand miss only {control_miss:.4%}"


@complex_only
@pytest.mark.integration
def test_post2_balance_is_dimensionless_against_a_reference_speed():
    """`mean_balance_rel_diff` at v = ω/|k| is |r − v|/max(r, v); ``None`` without v."""
    import numpy as np
    from mpi4py import MPI

    from fem_em_solver.utils.constants import EPSILON_0, MU_0
    from tests.validation.test_lossy_plane_wave import (
        EPSILON_R,
        FREQUENCY_HZ,
        SIGMA,
        _wavenumber,
    )

    rank0 = MPI.COMM_WORLD.rank == 0
    omega = 2.0 * np.pi * FREQUENCY_HZ
    v_cf = omega / abs(_wavenumber(SIGMA))
    v_blind = 1.0 / np.sqrt(MU_0 * EPSILON_0) / np.sqrt(EPSILON_R)

    def imbalanced(d):
        return any("strongly imbalanced" in w for w in d["warnings"])

    readings = {}
    for n in (24, 12):
        diag_none, e_stats, b_stats, *_rest = _post2_read_diagnostics(n)  # cached: no new solve
        at_cf = compute_field_consistency_diagnostics(
            e_stats, b_stats, reference_speed_m_per_s=v_cf
        )
        at_blind = compute_field_consistency_diagnostics(
            e_stats, b_stats, reference_speed_m_per_s=v_blind
        )
        at_one = compute_field_consistency_diagnostics(
            e_stats, b_stats, reference_speed_m_per_s=1.0
        )
        tight_cf = compute_field_consistency_diagnostics(
            e_stats,
            b_stats,
            reference_speed_m_per_s=v_cf,
            max_mean_balance_rel_diff_warn=POST2_BAND,
        )
        tight_blind = compute_field_consistency_diagnostics(
            e_stats,
            b_stats,
            reference_speed_m_per_s=v_blind,
            max_mean_balance_rel_diff_warn=POST2_BAND,
        )
        r = at_cf["e_to_b_mean_ratio"]
        identity = abs(r - v_cf) / max(r, v_cf)
        readings[n] = (diag_none, at_cf, at_blind, at_one, tight_cf, tight_blind, identity)
        if rank0:
            print(
                f"\nPOST-2 balance n={n}: v=omega/|k|={v_cf:.6e} m/s, r={r:.6e} m/s"
                f"\n  balance at omega/|k| = {at_cf['mean_balance_rel_diff']:.6e}"
                f" vs |r-v|/max(r,v) = {identity:.6e}; warning: {imbalanced(at_cf)}"
                f"\n  no speed: mean_balance_rel_diff={diag_none['mean_balance_rel_diff']!r};"
                f" warning: {imbalanced(diag_none)}"
                f"\n  control v=c/sqrt(78)={v_blind:.6e}: balance = {at_blind['mean_balance_rel_diff']:.6e}"
                f" (predicted 21.20 %); warning at 0.95: {imbalanced(at_blind)};"
                f" at {POST2_BAND}: {imbalanced(tight_blind)} (omega/|k| at {POST2_BAND}: {imbalanced(tight_cf)})"
                f"\n  v=1.0 (the old cross-unit reading): balance = {at_one['mean_balance_rel_diff']:.12f};"
                f" warning: {imbalanced(at_one)}"
            )

    diag_none, at_cf, at_blind, at_one, tight_cf, tight_blind, identity = readings[24]
    # Anchor: the dimensionless balance is the ratio miss already gated by POST-2.
    assert abs(at_cf["mean_balance_rel_diff"] / identity - 1.0) <= 1e-9
    assert at_cf["mean_balance_rel_diff"] < POST2_BAND
    assert not imbalanced(at_cf)
    # No speed: the key is present, None, and nothing is called imbalanced.
    assert "mean_balance_rel_diff" in diag_none and diag_none["mean_balance_rel_diff"] is None
    assert not imbalanced(diag_none)
    # Negative control: the σ-blind speed misses the band; 0.95 stays silent.
    assert at_blind["mean_balance_rel_diff"] > POST2_BAND
    assert not imbalanced(at_blind)
    assert imbalanced(tight_blind) and not imbalanced(tight_cf)
    # v = 1.0 reproduces the old cross-unit reading and its warning.
    assert at_one["mean_balance_rel_diff"] > 0.95 and imbalanced(at_one)
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="reference_speed_m_per_s"):
            compute_field_consistency_diagnostics(
                {"min": 1.0, "max": 3.0, "mean": 2.0},
                {"min": 0.5, "max": 1.5, "mean": 1.0},
                reference_speed_m_per_s=bad,
            )
