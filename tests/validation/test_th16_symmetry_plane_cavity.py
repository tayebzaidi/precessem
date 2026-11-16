"""Validation: per-face PEC / PMC on the cavity eigen path (TH-16 step 1).

A symmetry plane is a perfect electric (odd E_t) or perfect magnetic (even
E_t) wall. In the E-field curl-curl form the PMC wall n × ∇×E = 0 is the
*natural* condition (Jin §8.7.5): no ``ds`` term, only the PEC facets pinned.
Gated three ways at `TH-9`'s imported 1 % band, degree 2:

(i)   one PMC face (z+): per-axis k = (p + ½)π/d on the PMC axis, closed form
      from :func:`analytic_box_cavity_modes`;
(ii)  partition: the all-PEC 1.0 × 0.8 × 1.2 box's spectrum is the union of
      its magnetic-wall half (i) and electric-wall half (`TH-9`'s box)
      (Jin §7.1 Figs 7.6–7.7);
(iii) two PMC faces (x+, y+) — the birdcage's x = 0 / y = 0 cuts in box form;
(iv)  flag-off: ``pmc_faces=()`` reproduces ``solve_pec_cavity_modes``.

Not claimed: a mesh cutter, port half-impedance bookkeeping, or the
quarter-birdcage gate (later TH-16 steps).
"""

import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.core.cavity import (
    analytic_box_cavity_modes,
    solve_box_cavity_modes,
    solve_pec_cavity_modes,
)

from tests.validation.test_cavity_resonances import (  # TH-9's 1 % band
    FREQUENCY_TOLERANCE_PCT,
)

EDGES = (1.0, 0.8, 0.6)
DIVISIONS = (6, 5, 4)
FULL_EDGES = (1.0, 0.8, 1.2)
FULL_DIVISIONS = (6, 5, 8)

# Closed-form values quoted in PROJECT_PLAN §9 item 65 [MHz] — cross-checked
# here against ``analytic_box_cavity_modes`` so the helper cannot drift.
QUOTED_Z_PMC = (195.12, 225.19, 270.52, 270.52, 324.78, 374.95)
QUOTED_PEC_HALF = (239.95, 291.35, 312.28, 346.40, 346.40, 353.53)
QUOTED_XY_PMC = (119.98, 243.58, 277.14, 277.14, 290.88, 348.92, 348.92, 359.93)

# Owned facets on the PMC planes of a (6, 5, 4) box, two triangles per square.
EXPECTED_PMC_FACETS_Z = 6 * 5 * 2  # 60
EXPECTED_PMC_FACETS_XY = 5 * 4 * 2 + 6 * 4 * 2  # 88

FLAG_OFF_RTOL = 1e-10

# Gradient (null-space) eigenvalues returned by the shift-invert window are
# counted and discarded by the solver, never allowed into the physical band.
# First draft asserted ``null_mode_count == 0`` (TH-9's premise at n_modes = 4);
# measured 2026-09-29 (20260929T033604Z_TH-16-step1-W1.log) the window at
# nev = n_modes + 4 returns 2 / 3 / 4 of them for (i) / (ii) / (iii) — (ii) is
# all-PEC, so this is the window width, not the PMC route. Every physical mode
# still matched mode-for-mode (worst 0.0508 %). The honest check is TH-9's own
# null-cluster one: whatever null eigenvalues come back sit at machine zero.
NULL_CLUSTER_RATIO = 1e-8
# Negative control: all-PEC half vs (i)'s closed form. Closed-form gap 22.98 %
# is the ceiling; asserted strictly above 10 %.
NEGATIVE_CONTROL_FLOOR_PCT = 10.0


def _rel_pct(numeric, exact):
    return 100.0 * np.abs(np.asarray(numeric) - np.asarray(exact)) / np.asarray(exact)


def _report(tag, spectrum, analytic, elapsed):
    err = _rel_pct(spectrum.frequencies_hz, analytic)
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[TH-16] {tag}: cells={spectrum.n_cells} dofs={spectrum.n_dofs} "
            f"converged={spectrum.n_converged} null_discarded={spectrum.null_mode_count} "
            f"null_max/lam0={spectrum.diagnostics.get('null_max_abs_over_lam0', 0):.2e} "
            f"pmc_facets_owned={spectrum.diagnostics.get('n_pmc_facets_owned')} "
            f"elapsed={elapsed:.2f}s"
        )
        for f_num, f_exact, e in zip(spectrum.frequencies_hz, analytic, err):
            print(f"  {f_num/1e6:10.4f} MHz vs {f_exact/1e6:10.4f} MHz  ({e:.4f}%)")
        print(f"  max err {err.max():.4f}% (band {FREQUENCY_TOLERANCE_PCT}%)")
    return err


def test_closed_form_partition_multiset():
    """Numpy only: closed form (i) ∪ PEC half = full PEC box, as multisets."""

    n = 40
    full = analytic_box_cavity_modes(FULL_EDGES, (), n)
    half_pmc = analytic_box_cavity_modes(EDGES, ("z+",), n)
    half_pec = analytic_box_cavity_modes(EDGES, (), n)
    union = np.sort(np.concatenate([half_pmc, half_pec]))
    # Compare below the largest value both sides are complete up to.
    cap = min(half_pmc[-1], half_pec[-1], full[-1]) * (1 - 1e-12)
    lhs = union[union < cap]
    rhs = full[full < cap]
    assert lhs.size == rhs.size and lhs.size >= 20
    np.testing.assert_allclose(lhs, rhs, rtol=1e-12)

    # The quoted item-65 figures are what the helper produces.
    np.testing.assert_allclose(
        analytic_box_cavity_modes(EDGES, ("z+",), 6) / 1e6, QUOTED_Z_PMC, atol=0.006
    )
    np.testing.assert_allclose(
        analytic_box_cavity_modes(EDGES, (), 6) / 1e6, QUOTED_PEC_HALF, atol=0.006
    )
    np.testing.assert_allclose(
        analytic_box_cavity_modes(EDGES, ("x+", "y+"), 8) / 1e6,
        QUOTED_XY_PMC,
        atol=0.006,
    )


@pytest.mark.validation
def test_one_pmc_face_matches_closed_form():
    """(i) z+ PMC: lowest 5 within 1 %; 60 owned PMC facets (allreduced)."""

    t0 = time.perf_counter()
    spectrum = solve_box_cavity_modes(
        EDGES, DIVISIONS, degree=2, n_modes=5, pmc_faces=("z+",)
    )
    elapsed = time.perf_counter() - t0
    analytic = analytic_box_cavity_modes(EDGES, ("z+",), 5)
    err = _report("(i) z+ PMC", spectrum, analytic, elapsed)

    assert spectrum.diagnostics["n_pmc_facets_owned"] == EXPECTED_PMC_FACETS_Z
    assert spectrum.diagnostics["null_max_abs_over_lam0"] < NULL_CLUSTER_RATIO
    assert np.all(err < FREQUENCY_TOLERANCE_PCT), err

    # Negative control: the all-PEC half against (i)'s closed form.
    pec = solve_pec_cavity_modes(EDGES, DIVISIONS, degree=2, n_modes=4)
    gap = float(_rel_pct(pec.frequencies_hz[0], analytic[0]))
    closed_gap = float(
        _rel_pct(analytic_box_cavity_modes(EDGES, (), 1)[0], analytic[0])
    )
    wrong_face = analytic_box_cavity_modes(EDGES, ("x+",), 1)[0]
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"  [neg] all-PEC mode 1 {pec.frequencies_hz[0]/1e6:.4f} MHz vs (i) "
            f"{analytic[0]/1e6:.4f} MHz: gap {gap:.2f}% (closed-form gap "
            f"{closed_gap:.2f}%, floor {NEGATIVE_CONTROL_FLOOR_PCT}%)"
        )
        print(
            f"  [predicted, printed] x+ for z+: mode 1 at {wrong_face/1e6:.2f} MHz "
            f"({float(_rel_pct(wrong_face, analytic[0])):.2f}% from (i))"
        )
    assert gap > NEGATIVE_CONTROL_FLOOR_PCT


@pytest.mark.validation
def test_partition_full_box_is_union_of_halves():
    """(ii) All-PEC 1.0 × 0.8 × 1.2: lowest 11 vs sorted((i) ∪ PEC half)."""

    n = 11
    t0 = time.perf_counter()
    full = solve_box_cavity_modes(FULL_EDGES, FULL_DIVISIONS, degree=2, n_modes=n)
    elapsed = time.perf_counter() - t0
    union = np.sort(
        np.concatenate(
            [
                analytic_box_cavity_modes(EDGES, ("z+",), n),
                analytic_box_cavity_modes(EDGES, (), n),
            ]
        )
    )[:n]
    err = _report("(ii) full box vs union of halves", full, union, elapsed)
    assert full.diagnostics["n_pmc_facets_owned"] == 0
    assert full.diagnostics["null_max_abs_over_lam0"] < NULL_CLUSTER_RATIO
    assert np.all(err < FREQUENCY_TOLERANCE_PCT), err


@pytest.mark.validation
def test_two_pmc_faces_match_closed_form():
    """(iii) x+, y+ PMC: lowest 7 within 1 %; 88 owned PMC facets."""

    t0 = time.perf_counter()
    spectrum = solve_box_cavity_modes(
        EDGES, DIVISIONS, degree=2, n_modes=7, pmc_faces=("x+", "y+")
    )
    elapsed = time.perf_counter() - t0
    analytic = analytic_box_cavity_modes(EDGES, ("x+", "y+"), 7)
    err = _report("(iii) x+, y+ PMC", spectrum, analytic, elapsed)
    assert spectrum.diagnostics["n_pmc_facets_owned"] == EXPECTED_PMC_FACETS_XY
    assert spectrum.diagnostics["null_max_abs_over_lam0"] < NULL_CLUSTER_RATIO
    assert np.all(err < FREQUENCY_TOLERANCE_PCT), err


@pytest.mark.validation
def test_flag_off_reproduces_pec_solver():
    """(iv) ``pmc_faces=()`` equals ``solve_pec_cavity_modes`` to 1e-10 rel."""

    t0 = time.perf_counter()
    new = solve_box_cavity_modes(EDGES, DIVISIONS, degree=2, n_modes=4, pmc_faces=())
    old = solve_pec_cavity_modes(EDGES, DIVISIONS, degree=2, n_modes=4)
    elapsed = time.perf_counter() - t0
    rel = np.abs(new.frequencies_hz - old.frequencies_hz) / old.frequencies_hz
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[TH-16] (iv) flag-off max rel {rel.max():.3e} "
            f"(bound {FLAG_OFF_RTOL:g}); elapsed={elapsed:.2f}s"
        )
    assert new.diagnostics["n_pmc_facets_owned"] == 0
    assert np.all(rel <= FLAG_OFF_RTOL), rel
