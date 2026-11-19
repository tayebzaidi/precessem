"""`PORT-18` step 2 (§9 item 246) — the sheet's two one-sided current readings
on the exact-rectangle strips, for a solved field.

``sheet_terminal_current`` (``src/fem_em_solver/ports/lumped.py:338-364``)
reads ``I⁺ = (1/(R h)) ∫_N (E('+') + E_src)·ĥ dS`` through ``_restrict``, which
hard-codes ``'+'`` on an interior sheet.  `PORT-18` step 1 (2026-09-10) measured
the sheet normal as azimuthal, so ``ĥ = ẑ`` lies *in* the sheet and ``E·ĥ`` is
a tangential component N1curl keeps single-valued: ``I⁺ = I⁻`` is then an
identity up to round-off whichever cell DOLFINx calls ``'+'``.  Step 1 read
that only on the CAD face ``F_i`` and only for interpolated fields; this module
asserts it on the gated strip ``N_i`` of the exact-rectangle ×1 fixture
(``leg_sheet_strip_fraction = GATED_WIDTH_FRACTION``, κ off) for a P1-driven
degree-1 10 MHz solve.  It tests no collapse of a ragged-sheet discrepancy —
the 09-10 record holds none.

**Bands imported, never restated**: ``EXACTNESS_BAND`` and ``PREMISE_BAND``
from ``scripts/probes/port18_sheet_plus_side_census.py`` by that probe's own
``sys.path`` insert idiom.

**Anchors (asserted, every sheet):** (P0) fixture identity — cells, ``N_i``
area vs the SPEC rectangle, and (at the record width ``-n 2`` only) ``S`` vs
``ANS4_EXACT_STRIP_S_P1_RECORDS["10"]``; (P1) premise ``∫(n('+')·ĥ)² dS / A_i
≤ PREMISE_BAND``; (A) ``E_lin = (−z/ρ₀, 0, 1 + x/ρ₀)`` reads ``R^± = 1 +
x̄_i/ρ₀`` within ``EXACTNESS_BAND`` (relative denominator ``max(|target|, 1)``,
ratified 2026-09-10 10:30); (B0) the test's ``'+'`` mirror reproduces the
production ``I⁺``; (B) ``|I⁺ − I⁻|/|I⁺|`` and the jump metric
``J_i[E] = sqrt(∫|(E('+') − E('−'))·ĥ|² / ∫|E('+')·ĥ|²)`` within
``EXACTNESS_BAND`` (also ``J_i`` of the smooth interpolant).

**Negative control (asserted — arithmetic):** the planted DG0 field
``g = (χ_L − χ_U) ẑ`` jumps sign across every ``N_i`` facet, so ``J_i[g] = 2``
exactly; a jump metric blind to a discontinuity reads 0 and fails.

**Scope:** degree 1, 10 MHz, P1 drive, ×1 exact mesh.  No band, record, κ or
``src/`` line moves; a pass licenses no ``S₁₁`` / ``Z_in`` / ``f₀`` figure.

Run (complex build)::

    scripts/testing/run_and_log.sh PORT-18-step2-W1 "docker compose exec -T fem-em-solver \\
      bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && \\
      PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 300 \\
      mpiexec -n 2 python3 -m pytest tests/environment \\
      tests/validation/test_port18_exact_strip_one_sided.py -s -v --tb=short'"
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.ports import lumped
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
from tests.mesh.test_birdcage_port_tags import RING_RADIUS
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
    ANS4_EXACT_STRIP_S_P1_RECORDS,
    SPEC_AREA_M2,
)
from tests.validation.test_ans4_resolution_ladder import STEP3A_RECORD_RTOL
from tests.validation.test_birdcage_b1_plus_map import _solve_driven
from tests.validation.test_ops61_exact_strip_ports import (
    ANS4_EXACT_STRIP_CELL_RECORD,
    STRIP_BAND,
)
from tests.validation.test_port_birdcage_leg_offset_sweep import _four_port_rung
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

# The probe's own path-insert idiom (port18_sheet_plus_side_census.py:71-77).
_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
if str(_REPO_ROOT / "scripts" / "probes") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "scripts" / "probes"))

import port18_sheet_plus_side_census as port18_probe  # noqa: E402

EXACTNESS_BAND = port18_probe.EXACTNESS_BAND
PREMISE_BAND = port18_probe.PREMISE_BAND

PORTS = (1, 2, 3, 4)
DRIVEN = "P1"
# Step 1's smooth-interpolant four-sheet spread on F_i at x1
# (`20260910T093526Z_PORT-18.log:7012`) — printed comparand, GEO-32's lead.
STEP1_SMOOTH_SPREAD_F_X1 = 7.223213e-06
# The negative control's closed form: g('+')·ĥ = -g('-')·ĥ = ±1 on every facet.
CONTROL_JUMP = 2.0


def _rank0(text):
    if MPI.COMM_WORLD.rank == 0:
        print(text, flush=True)


class _PerSheet:
    """One compiled form per quantity, four assemblies selected by Constants
    (the probe's ``_PerSheet`` idiom, complex-valued, every number allreduced)."""

    def __init__(self, msh, facet_tags, comm):
        self.msh = msh
        self.facet_tags = facet_tags
        self.comm = comm
        self.w = [fem.Constant(msh, default_scalar_type(0.0)) for _ in PORTS]

    def eval(self, builder, q=None):
        metadata = {} if q is None else {"quadrature_degree": int(q)}
        dS = ufl.Measure("dS", domain=self.msh, subdomain_data=self.facet_tags, metadata=metadata)
        integrand = None
        for k, i in enumerate(PORTS):
            term = self.w[k] * builder(i) * dS(SHEET_IFACE + i)
            integrand = term if integrand is None else integrand + term
        form = fem.form(integrand)
        out = []
        for k in range(len(PORTS)):
            for j, c in enumerate(self.w):
                c.value = 1.0 if j == k else 0.0
            local = fem.assemble_scalar(form)
            out.append(complex(self.comm.allreduce(local, op=MPI.SUM)))
        return np.asarray(out, dtype=np.complex128)


def _one_sided_current(msh, facet_tags, sheet, e_field, comm, side):
    """``lumped.py:356-364`` verbatim with the restriction side as a parameter."""
    sheet.validate()
    ds_sheet = lumped._sheet_measure(msh, facet_tags, sheet)
    unit = sheet.unit_drive()
    h_hat = ufl.as_vector([default_scalar_type(c) for c in unit])
    e_src = complex(sheet.source_voltage_v) / sheet.gap_height_m
    r = lambda f: f(side)  # noqa: E731
    integrand = (
        ufl.inner(r(e_field), h_hat) + default_scalar_type(e_src)
    ) / sheet.sheet_resistivity
    local = fem.assemble_scalar(fem.form(integrand * ds_sheet))
    return complex(comm.allreduce(local, op=MPI.SUM)) / sheet.gap_height_m


def _jump_ratio(num, den, label):
    """``sqrt(num/den)`` in numpy after the reduction; imaginary parts printed.

    The numerator's imaginary part is held to ``1e-12`` of the *denominator's*
    real part (the quantity's own scale): the numerator of an exact identity is
    ≈ 0, where a bound relative to itself is undefined.
    """
    out = []
    for i, (a, b) in zip(PORTS, zip(num, den)):
        assert abs(a.imag) <= 1e-12 * abs(
            b.real
        ), f"{label} P{i}: numerator imag {a.imag:.3e} vs real {a.real:.3e}"
        assert abs(b.imag) <= 1e-12 * abs(
            b.real
        ), f"{label} P{i}: denominator imag {b.imag:.3e} vs real {b.real:.3e}"
        out.append(float(np.sqrt(max(a.real, 0.0) / b.real)))
    return np.asarray(out)


@pytest.fixture(scope="module")
def case():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    rung = _four_port_rung(
        "PORT-18 exact 10 MHz",
        np.zeros(4),
        degree=1,
        leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
    )
    t_rung = time.perf_counter() - t0
    msh = rung["mesh"]
    cell_tags = rung["cell_tags"]
    tags_f = rung["facet_tags"]
    tdim = msh.topology.dim

    t1 = time.perf_counter()
    solved = _solve_driven(rung, DRIVEN)
    t_solve = time.perf_counter() - t1
    E = solved["fields"].e_complex
    sheets = [s.sheet(driven=(s.port_id == DRIVEN)) for s in rung["specs"]]
    assert [s.port_id for s in sheets] == [f"P{i}" for i in PORTS]
    h_vec = sheets[0].unit_drive()
    for s in sheets:
        assert np.allclose(s.unit_drive(), h_vec)
    h_hat = ufl.as_vector([default_scalar_type(c) for c in h_vec])

    t2 = time.perf_counter()
    # Production I+ (from _solve_driven), the '+' mirror (B0) and the '-' mirror.
    i_prod = np.asarray([solved["currents"][s.port_id] for s in sheets])
    i_plus = np.asarray([_one_sided_current(msh, tags_f, s, E, comm, "+") for s in sheets])
    i_minus = np.asarray([_one_sided_current(msh, tags_f, s, E, comm, "-") for s in sheets])

    ps = _PerSheet(msh, tags_f, comm)
    x = ufl.SpatialCoordinate(msh)
    n = ufl.FacetNormal(msh)
    one = fem.Constant(msh, default_scalar_type(1.0))
    area = ps.eval(lambda i: one, 2).real
    premise = ps.eval(lambda i: ufl.dot(n("+"), h_hat) * ufl.dot(n("+"), h_hat), 2)
    premise = premise.real / area
    xb = ps.eval(lambda i: x("+")[0], 2).real / area
    yb = ps.eval(lambda i: x("+")[1], 2).real / area
    zb = ps.eval(lambda i: x("+")[2], 2).real / area
    rho0 = float(np.mean(np.hypot(xb, yb)))
    zbar = float(np.mean(zb))
    gap_h = float(np.mean([s.gap_height_m for s in sheets]))
    lam = 4.0 * gap_h

    def _jump(field):
        d = lambda: ufl.dot(field("+") - field("-"), h_hat)  # noqa: E731
        p = lambda: ufl.dot(field("+"), h_hat)  # noqa: E731
        num = ps.eval(lambda i: ufl.inner(d(), d()))
        den = ps.eval(lambda i: ufl.inner(p(), p()))
        return num, den

    # (B) jump metric on the solved field.
    jE_num, jE_den = _jump(E)
    # Interpolated fields (the probe's definitions, :243-254, :276-281).
    V = E.function_space
    Es = fem.Function(V)
    Es.interpolate(
        lambda X: np.vstack(
            (
                np.zeros(X.shape[1]),
                np.zeros(X.shape[1]),
                np.cos(2.0 * np.pi * (X[2] - zbar) / lam)
                * (1.0 + (X[0] ** 2 + X[1] ** 2) / rho0**2),
            )
        )
    )
    Es.x.scatter_forward()
    jS_num, jS_den = _jump(Es)
    rs_plus = ps.eval(lambda i: ufl.dot(Es("+"), h_hat), 2) / area
    El = fem.Function(V)
    El.interpolate(lambda X: np.vstack((-X[2] / rho0, np.zeros(X.shape[1]), 1.0 + X[0] / rho0)))
    El.x.scatter_forward()
    lin_plus = ps.eval(lambda i: ufl.dot(El("+"), h_hat), 2) / area
    lin_minus = ps.eval(lambda i: ufl.dot(El("-"), h_hat), 2) / area

    # Planted control: g = (chi_L - chi_U) z-hat, DG0 on owned + ghost cells.
    n_owned = int(msh.topology.index_map(tdim).size_local)
    Q = fem.functionspace(msh, ("DG", 0))
    chiL = {i: port18_probe._indicator(Q, cell_tags, PORT_LOWER + i, n_owned) for i in PORTS}
    chiU = {i: port18_probe._indicator(Q, cell_tags, PORT_UPPER + i, n_owned) for i in PORTS}
    ez = ufl.as_vector(
        [default_scalar_type(0.0), default_scalar_type(0.0), default_scalar_type(1.0)]
    )

    def _gnum(i):
        d = ufl.dot((chiL[i] - chiU[i])("+") * ez - (chiL[i] - chiU[i])("-") * ez, h_hat)
        return ufl.inner(d, d)

    def _gden(i):
        p = ufl.dot((chiL[i] - chiU[i])("+") * ez, h_hat)
        return ufl.inner(p, p)

    jG_num = ps.eval(_gnum)
    jG_den = ps.eval(_gden)
    f_plus = ps.eval(lambda i: chiL[i]("+"), 2).real / area
    g_mean_plus = ps.eval(lambda i: ufl.dot((chiL[i] - chiU[i])("+") * ez, h_hat), 2).real / area
    g_mean_minus = ps.eval(lambda i: ufl.dot((chiL[i] - chiU[i])("-") * ez, h_hat), 2).real / area
    t_read = time.perf_counter() - t2

    jE = _jump_ratio(jE_num, jE_den, "J[E]")
    jS = _jump_ratio(jS_num, jS_den, "J[smooth]")
    jG = _jump_ratio(jG_num, jG_den, "J[g]")
    target = 1.0 + xb / rho0
    denom = np.maximum(np.abs(target), 1.0)
    lin_plus_rel = np.abs(lin_plus - target) / denom
    lin_minus_rel = np.abs(lin_minus - target) / denom
    b0_rel = np.abs(i_plus - i_prod) / np.abs(i_prod)
    b_rel = np.abs(i_plus - i_minus) / np.abs(i_plus)
    s_rec = np.asarray(ANS4_EXACT_STRIP_S_P1_RECORDS["10"], dtype=np.complex128)
    s_rel = float(np.linalg.norm(rung["s"] - s_rec) / np.linalg.norm(s_rec))
    rs_real = rs_plus.real
    smooth_spread = float((rs_real.max() - rs_real.min()) / abs(rs_real.mean()))

    _rank0(
        f"\n[PORT-18 step2] -n {comm.size}; cells {rung['cells']} (record "
        f"{ANS4_EXACT_STRIP_CELL_RECORD}); rung {t_rung:.2f} s (box-speed gate 60 s at "
        f"-n 2), sweep {rung['sweep_time']:.2f} s, P1 solve {solved['solve_time']:.2f} s "
        f"(+ read {t_solve - solved['solve_time']:.2f} s), one-sided readings {t_read:.2f} s"
    )
    _rank0(
        f"[PORT-18 step2] ||S - S_rec(10)||/||S_rec|| = {s_rel:.3e} (band "
        f"{STEP3A_RECORD_RTOL:.0e}, asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH} only)"
    )
    _rank0(
        f"[PORT-18 step2] h_hat {tuple(float(c) for c in h_vec)}; rho0 {rho0:.12e} "
        f"(RING_RADIUS {RING_RADIUS}); zbar {zbar:.12e}; h {gap_h:.12e}; lambda {lam:.12e}; "
        f"EXACTNESS_BAND {EXACTNESS_BAND:.0e}, PREMISE_BAND {PREMISE_BAND:.0e} (imported)"
    )
    for k, (i, s) in enumerate(zip(PORTS, rung["sheets"])):
        _rank0(
            f"  P{i}: facets {s['facets']}, N_i area (rung) {s['area']:.12e} rel SPEC "
            f"{s['area'] / SPEC_AREA_M2 - 1:+.3e}, assembled A {area[k]:.12e}; premise "
            f"{premise[k]:.3e}; xbar {xb[k]:+.9e} ybar {yb[k]:+.9e} zbar {zb[k]:.9e}"
        )
        _rank0(
            f"     I+ prod {i_prod[k]:+.10e}  |I+| {abs(i_prod[k]):.10e}\n"
            f"     I+ mirror {i_plus[k]:+.10e}  (B0 rel {b0_rel[k]:.3e})\n"
            f"     I- mirror {i_minus[k]:+.10e}  |I+ - I-|/|I+| {b_rel[k]:.3e} "
            f"(predicted 1e-15 .. 1e-13)"
        )
        _rank0(
            f"     J[E] {jE[k]:.3e} (num imag {jE_num[k].imag:+.2e}); J[smooth] "
            f"{jS[k]:.3e} (num imag {jS_num[k].imag:+.2e}); J[g] {jG[k]:.15f} "
            f"(|J-2| {abs(jG[k] - CONTROL_JUMP):.3e}); f_i {f_plus[k]:.6f}; g one-sided "
            f"means + {g_mean_plus[k]:+.6f} - {g_mean_minus[k]:+.6f} "
            f"(2f-1 = {2 * f_plus[k] - 1:+.6f})"
        )
        _rank0(
            f"     E_lin R+ {lin_plus[k].real:+.15e} R- {lin_minus[k].real:+.15e} "
            f"target {target[k]:+.15e}; rel {lin_plus_rel[k]:.3e} / {lin_minus_rel[k]:.3e}; "
            f"smooth R+ {rs_plus[k].real:.12e}"
        )
    _rank0(
        f"[PORT-18 step2] smooth-interpolant four-sheet spread on N_i {smooth_spread:.6e} "
        f"(step 1 on F_i x1: {STEP1_SMOOTH_SPREAD_F_X1:.6e}; printed, GEO-32's lead)"
    )
    return {
        "comm_size": comm.size,
        "cells": int(rung["cells"]),
        "areas_rung": np.asarray([s["area"] for s in rung["sheets"]]),
        "s_rel": s_rel,
        "premise": premise,
        "lin_plus_rel": lin_plus_rel,
        "lin_minus_rel": lin_minus_rel,
        "b0_rel": b0_rel,
        "b_rel": b_rel,
        "jE": jE,
        "jS": jS,
        "jG": jG,
    }


@complex_only
def test_p0_fixture_identity(case):
    ratio = case["cells"] / ANS4_EXACT_STRIP_CELL_RECORD
    assert abs(ratio - 1.0) <= CELL_COUNT_BAND, f"cells {case['cells']}"
    for i, a in zip(PORTS, case["areas_rung"]):
        assert abs(a / SPEC_AREA_M2 - 1.0) <= STRIP_BAND, f"P{i} area {a:.12e}"
    if case["comm_size"] == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH:
        assert case["s_rel"] <= STEP3A_RECORD_RTOL, f"S vs record {case['s_rel']:.3e}"


@complex_only
def test_p1_premise_h_lies_in_the_sheet(case):
    for i, v in zip(PORTS, case["premise"]):
        assert v <= PREMISE_BAND, f"P{i}: (n('+')·h)^2 mean {v:.3e}"


@complex_only
def test_a_e_lin_closed_form(case):
    for i, a, b in zip(PORTS, case["lin_plus_rel"], case["lin_minus_rel"]):
        assert a <= EXACTNESS_BAND, f"P{i}: E_lin R+ rel {a:.3e}"
        assert b <= EXACTNESS_BAND, f"P{i}: E_lin R- rel {b:.3e}"


@complex_only
def test_b0_mirror_reproduces_production(case):
    for i, v in zip(PORTS, case["b0_rel"]):
        assert v <= EXACTNESS_BAND, f"P{i}: '+' mirror vs production {v:.3e}"


@complex_only
def test_b_one_sided_currents_agree(case):
    for i, v, j, js in zip(PORTS, case["b_rel"], case["jE"], case["jS"]):
        assert v <= EXACTNESS_BAND, f"P{i}: |I+ - I-|/|I+| {v:.3e}"
        assert j <= EXACTNESS_BAND, f"P{i}: J[E] {j:.3e}"
        assert js <= EXACTNESS_BAND, f"P{i}: J[smooth] {js:.3e}"


@complex_only
def test_negative_control_planted_dg0_jump_reads_two(case):
    for i, j in zip(PORTS, case["jG"]):
        assert abs(j - CONTROL_JUMP) <= EXACTNESS_BAND, f"P{i}: J[g] {j:.15f}"
