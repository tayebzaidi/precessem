"""Validation test: layered (foil-on-backing) surface impedance on the TH-14 cavity (TH-18 step 1).

`TH-14`'s fixture — the 1.0 × 0.8 × 0.6 m box, (9, 7, 6), degree 2, σ = 1e4 S/m,
TE₁₀₁ at ω₀ ≈ 2π·291.35 MHz — with the Leontovich ``(1 + j)R_s`` replaced by
the layered input impedance of ``core.layered_impedance`` (Jin §3.4.2
(3.87)–(3.89)) through ``solve_impedance_wall_cavity_mode(surface_impedance_fn=)``.
The layer thickness t is fixed in metres from δ(ω₀); δ follows ω_lin, t does not.

Perturbation carried to Z_s = R + jX: Q = Q_c·R_s/R, Δf/f = −(X/R_s)/(2Q_c).

Solves: ``pec`` (Dirichlet control), ``leo`` (TH-14 Leontovich), ``t5`` (5δ on
PEC), ``thin`` (δ/2 on PEC), ``air`` (δ/2 on η₀, printed only).

Asserted (PROJECT_PLAN §9 item 82, 2026-09-29 10:30 review):
  (P)   premise: Q_leo = 8.017679e+02 (`20260929T033804Z_TH-16-step1-W2-ruleC.log:103`) to 1e-5;
  (i)   |Q_t5/Q_leo − 0.999874426| ≤ 1e-5 (re-derived from Z_c·tanh(5(1+j)));
  (ii)  Q_thin/Q_leo = 6.242716 within 0.5 %, and Q_thin = Q_c·6.242716 = 5004.69
        within the imported `Q_TOLERANCE` (5 %);
  (iii) (f_thin − f_pec)/(f_leo − f_pec) = 0.967980 within 2 % — sees X and its sign;
  negative control (rule (e)): bulk Leontovich against the thin-layer
        prediction, 1 − Q_leo/5004.69 (≈ 0.8398) > 0.5.
Printed only: ``air`` Q (hand-computed ≈ Q_c/2.011 ≈ 399).
"""

import time

import numpy as np
import pytest
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core.cavity import (
    pozar_te10l_q_conductor,
    solve_impedance_wall_cavity_mode,
)
from fem_em_solver.core.layered_impedance import (
    ConductorLayer,
    conductor_skin_depth_m,
    conductor_surface_resistance_ohm,
    layered_surface_impedance,
)
from fem_em_solver.utils.constants import C_0, ETA_0
from tests.validation.test_cavity_leontovich_q import Q_TOLERANCE

pytestmark = pytest.mark.skipif(
    not np.issubdtype(PETSc.ScalarType, np.complexfloating),
    reason="surface-impedance cavity needs the complex DolfinX build",
)

EDGES = (1.0, 0.8, 0.6)
FINE = (9, 7, 6)
DEGREE = 2
SIGMA = 1.0e4

# Pre-registered (§9 item 82).
Q_LEO_RECORD = 8.017679e02  # ruleC.log:103
PREMISE_TOLERANCE = 1e-5
T5_RATIO = 0.999874426
T5_TOLERANCE = 1e-5
THIN_RATIO = 6.242716
THIN_RATIO_TOLERANCE = 0.005
Q_THIN_PREDICTED = 5004.69  # Q_c · 6.242716
SHIFT_RATIO = 0.967980
SHIFT_TOLERANCE = 0.02
NEGATIVE_CONTROL_FLOOR = 0.5


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _omega0():
    a, _b, d = EDGES
    return float(np.pi * np.sqrt((1.0 / a) ** 2 + (1.0 / d) ** 2) * C_0)


@pytest.fixture(scope="module")
def solves():
    comm = MPI.COMM_WORLD
    omega0 = _omega0()
    delta0 = conductor_skin_depth_m(omega0, SIGMA)
    t5 = 5.0 * delta0
    t_half = 0.5 * delta0

    def stack(thickness, backing):
        return lambda w: layered_surface_impedance(w, [ConductorLayer(SIGMA, thickness)], backing)

    cases = (
        ("pec", dict(sigma_s_per_m=None)),
        ("leo", dict(sigma_s_per_m=SIGMA)),
        ("t5", dict(surface_impedance_fn=stack(t5, 0.0))),
        ("thin", dict(surface_impedance_fn=stack(t_half, 0.0))),
        ("air", dict(surface_impedance_fn=stack(t_half, ETA_0))),
    )
    out = {}
    t0 = time.perf_counter()
    _say(
        f"[TH-18] omega0={omega0:.9e} rad/s delta0={delta0:.9e} m t5={t5:.9e} m t_half={t_half:.9e} m"
    )
    for key, kw in cases:
        t = time.perf_counter()
        out[key] = solve_impedance_wall_cavity_mode(
            edges=EDGES, divisions=FINE, degree=DEGREE, comm=comm, **kw
        )
        r = out[key]
        _say(
            f"[TH-18] {key:5s} cells={r.n_cells} dofs={r.n_dofs} nconv={r.n_converged} "
            f"Z_s={r.surface_impedance_ohm} lambda={r.eigenvalue:.10e} "
            f"f={r.omega_rad_s / (2 * np.pi) / 1e6:.6f} MHz Q={r.q:.6e} "
            f"({time.perf_counter() - t:.1f}s)"
        )
    r_s = conductor_surface_resistance_ohm(omega0, SIGMA)
    for key in ("t5", "thin", "air"):
        z = out[key].surface_impedance_ohm
        _say(f"[TH-18] {key}: R/R_s={z.real / r_s:.10f} X/R_s={z.imag / r_s:.10f}")
    for key in ("leo", "t5", "thin", "air"):
        # Printed only (e^{jωt}: Im ω > 0 is lossy).
        _say(f"[TH-18] Im omega {key}={out[key].omega_rad_s.imag:+.6e} rad/s")
    _say(f"[TH-18] all solves elapsed={time.perf_counter() - t0:.1f}s")
    return out


@pytest.mark.validation
def test_premise_leo_reproduces_th14_record(solves):
    q = solves["leo"].q
    rel = abs(q / Q_LEO_RECORD - 1.0)
    _say(
        f"[TH-18] (P) Q_leo={q:.7e} vs record {Q_LEO_RECORD:.6e} rel={rel:.3e} (bound {PREMISE_TOLERANCE:.0e})"
    )
    assert rel <= PREMISE_TOLERANCE


@pytest.mark.validation
def test_i_five_skin_depths_is_leontovich(solves):
    ratio = solves["t5"].q / solves["leo"].q
    _say(
        f"[TH-18] (i) Q_t5/Q_leo={ratio:.10f} vs {T5_RATIO} |diff|={abs(ratio - T5_RATIO):.3e} "
        f"(bound {T5_TOLERANCE:.0e})"
    )
    assert abs(ratio - T5_RATIO) <= T5_TOLERANCE


@pytest.mark.validation
def test_ii_half_skin_depth_q(solves):
    q_c = pozar_te10l_q_conductor(EDGES, 1, SIGMA)
    q_thin = solves["thin"].q
    ratio = q_thin / solves["leo"].q
    miss_ratio = ratio / THIN_RATIO - 1.0
    miss_q = q_thin / Q_THIN_PREDICTED - 1.0
    _say(
        f"[TH-18] (ii) Q_thin/Q_leo={ratio:.7f} vs {THIN_RATIO} miss={100 * miss_ratio:+.4f}% "
        f"(band {100 * THIN_RATIO_TOLERANCE:.1f}%); Q_thin={q_thin:.6e} vs Q_c*{THIN_RATIO}="
        f"{Q_THIN_PREDICTED} (Q_c={q_c:.6e}, Q_c*ratio={q_c * THIN_RATIO:.4f}) miss={100 * miss_q:+.4f}% "
        f"(band {100 * Q_TOLERANCE:.0f}%)"
    )
    assert abs(miss_ratio) <= THIN_RATIO_TOLERANCE
    assert abs(miss_q) <= Q_TOLERANCE


@pytest.mark.validation
def test_iii_half_skin_depth_frequency_shift(solves):
    f = {k: solves[k].omega_rad_s.real / (2 * np.pi) for k in ("pec", "leo", "thin", "air")}
    ratio = (f["thin"] - f["pec"]) / (f["leo"] - f["pec"])
    miss = ratio / SHIFT_RATIO - 1.0
    _say(
        f"[TH-18] (iii) f_pec={f['pec'] / 1e6:.6f} f_leo={f['leo'] / 1e6:.6f} "
        f"f_thin={f['thin'] / 1e6:.6f} MHz; (f_thin-f_pec)/(f_leo-f_pec)={ratio:.7f} vs {SHIFT_RATIO} "
        f"miss={100 * miss:+.4f}% (band {100 * SHIFT_TOLERANCE:.0f}%)"
    )
    q_c = pozar_te10l_q_conductor(EDGES, 1, SIGMA)
    _say(
        f"[TH-18] printed: air (delta/2 on eta0) Q={solves['air'].q:.6e} Q_c/Q_air={q_c / solves['air'].q:.4f} "
        f"(hand-computed ~ Q_c/2.011 ~ 399); f_air={f['air'] / 1e6:.6f} MHz"
    )
    assert abs(miss) <= SHIFT_TOLERANCE


@pytest.mark.validation
def test_negative_control_bulk_leontovich_misses_thin_layer(solves):
    gap = 1.0 - solves["leo"].q / Q_THIN_PREDICTED
    _say(
        f"[TH-18] negative control: 1 - Q_leo/{Q_THIN_PREDICTED} = {gap:.4f} (expected 0.8398; floor {NEGATIVE_CONTROL_FLOOR})"
    )
    assert gap > NEGATIVE_CONTROL_FLOOR
