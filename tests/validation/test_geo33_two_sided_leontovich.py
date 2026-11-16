"""GEO-33 step 0b: the two-sided Leontovich term on step 0a's duplicated lips.

The cracked box (``core.cracked_box``, y = 0.4 plane duplicated) is two
separate 1.0 x 0.4 x 0.6 half-cavities; each lip is an exterior facet with its
own outward normal, so TH-14's ``ds`` term acts on both sides of the sheet.
For TE101 k does not depend on b, so R_s is common and Pozar gives

    Q_sep / Q0 = [0.4 / 1.7888] / [0.8 / 2.7616] = 0.7719141.

Anchors (asserted; TH-14's 5 % band and TH-9's 1 % band, imported):
(Q0) merged mesh, walls only, Q within Q_TOLERANCE of Pozar (801.6849);
(Qs) cracked mesh, Q_sep within Q_TOLERANCE of Pozar((1.0, 0.4, 0.6));
(R)  Q_sep / Q0 within Q_TOLERANCE of 0.771914;
(F)  Re f within FREQUENCY_TOLERANCE_PCT of 291.3459 MHz on all three solves;
(D)  Im omega > 0 on all three solves.
Negative control (rule (e)): the merged mesh with the single-valued interior
sheet on y = 0.4 (R = Z_s/2). TE101 has E_t = 0 there, so Q_sheet/Q0 = 1.000:
asserted within Q_TOLERANCE of 1 and outside (R)'s band, and no larger miss
from 0.771914 than the ceiling 1/0.771914 - 1 (a sheet only adds loss).
"""
import time

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.core.cavity import (
    pozar_te10l_q_conductor,
    solve_impedance_mesh_mode,
)
from fem_em_solver.core.cracked_box import LIP_LOWER, LIP_UPPER, cracked_box_mesh
from fem_em_solver.utils.constants import C_0
from tests.validation.test_cavity_leontovich_q import Q_TOLERANCE
from tests.validation.test_cavity_resonances import FREQUENCY_TOLERANCE_PCT

pytestmark = pytest.mark.skipif(
    not np.issubdtype(PETSc.ScalarType, np.complexfloating),
    reason="surface-impedance cavity needs the complex DolfinX build",
)

EDGES = (1.0, 0.8, 0.6)
HALF_EDGES = (1.0, 0.4, 0.6)
SIGMA_GATE = 1.0e4
SIGMA_SCALING = 1.0e6
DEGREE = 2
RATIO_LITERAL = 0.771914  # Pozar, two half-cavities (03:00 review arithmetic)
RATIO_LITERAL_TOL = 1e-6
F_TE101_MHZ = 291.3459
SEPTUM_TAG = 21
K0_SQ = float(np.pi**2 * ((1.0 / EDGES[0]) ** 2 + (1.0 / EDGES[2]) ** 2))


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(f"[GEO-33 0b] {msg}", flush=True)


def _area(msh, tags, tag):
    ds = ufl.Measure("ds", domain=msh, subdomain_data=tags)
    one = fem.Constant(msh, PETSc.ScalarType(1.0))
    local = fem.assemble_scalar(fem.form(one * ds(tag)))
    return float(np.real(msh.comm.allreduce(local, op=MPI.SUM)))


def _check_mode(label, res, i=0):
    omega = res.omega(i)
    f_mhz = omega.real / (2.0 * np.pi) / 1e6
    err_pct = 100.0 * abs(f_mhz / F_TE101_MHZ - 1.0)
    _say(f"{label}[{i}]: lambda={complex(res.eigenvalues[i])!r}, Re f={f_mhz:.6f} MHz "
         f"({err_pct:.4f} % vs {F_TE101_MHZ}), Im omega={omega.imag:.6e}, Q={res.q(i):.6e}")
    assert err_pct < FREQUENCY_TOLERANCE_PCT, (label, f_mhz)  # (F)
    assert omega.imag > 0.0, (label, omega)  # (D)
    return res.q(i)


def test_pozar_ratio_literal():
    q0 = pozar_te10l_q_conductor(EDGES, 1, SIGMA_GATE)
    qs = pozar_te10l_q_conductor(HALF_EDGES, 1, SIGMA_GATE)
    ratio = qs / q0
    _say(f"Pozar Q0={q0:.6f}, Q_sep={qs:.6f}, ratio={ratio:.8f}; "
         f"1/Q_sep - 1/Q0 = {1 / qs - 1 / q0:.6e}")
    assert abs(ratio - RATIO_LITERAL) < RATIO_LITERAL_TOL, ratio


def test_two_sided_lips_against_pozar_and_sheet_control():
    t0 = time.perf_counter()
    merged, mtags, spread_m, septum = cracked_box_mesh(merge=True, septum_tag=SEPTUM_TAG)
    cracked, ctags, spread_c = cracked_box_mesh(merge=False)
    assert merged is not None and cracked is not None, (spread_m, spread_c)
    lo, up = _area(cracked, ctags, LIP_LOWER), _area(cracked, ctags, LIP_UPPER)
    n_septum = merged.comm.allreduce(
        int(np.sum(septum.indices < merged.topology.index_map(2).size_local)), op=MPI.SUM
    )
    _say(f"lip areas lower={lo!r} upper={up!r}; merged septum owned facets {n_septum}")

    q0_pozar = pozar_te10l_q_conductor(EDGES, 1, SIGMA_GATE)
    qs_pozar = pozar_te10l_q_conductor(HALF_EDGES, 1, SIGMA_GATE)

    # (Q0) merged mesh, walls only.
    r0 = solve_impedance_mesh_mode(merged, sigma_s_per_m=SIGMA_GATE, k0_sq=K0_SQ, degree=DEGREE)
    _say(f"merged walls: dofs {r0.n_dofs}, converged {r0.n_converged}")
    q0 = _check_mode("Q0", r0)
    e0 = abs(q0 / q0_pozar - 1.0)
    _say(f"(Q0) Q0={q0:.6e} vs Pozar {q0_pozar:.6e}: {100 * e0:.4f} %")

    # (Qs) cracked mesh: both members of the degenerate pair.
    rs = solve_impedance_mesh_mode(cracked, sigma_s_per_m=SIGMA_GATE, k0_sq=K0_SQ, degree=DEGREE)
    _say(f"cracked: dofs {rs.n_dofs}, converged {rs.n_converged}")
    qs = _check_mode("Qs", rs, 0)
    if rs.n_converged >= 2:
        qs_b = _check_mode("Qs", rs, 1)
        _say(f"pair member 2 Q={qs_b:.6e} ({100 * abs(qs_b / qs_pozar - 1):.4f} % vs Pozar)")
    else:
        _say("pair member 2 unconverged")
    es = abs(qs / qs_pozar - 1.0)
    ratio = qs / q0
    er = abs(ratio / RATIO_LITERAL - 1.0)
    _say(f"(Qs) Q_sep={qs:.6e} vs Pozar {qs_pozar:.6e}: {100 * es:.4f} %")
    _say(f"(R) Q_sep/Q0={ratio:.6f} vs {RATIO_LITERAL}: {100 * er:.4f} %; "
         f"1/Q_sep - 1/Q0 = {1 / qs - 1 / q0:.6e} (predicted 3.686e-4)")

    # Negative control: merged mesh + single-valued septum sheet.
    rc = solve_impedance_mesh_mode(
        merged, sigma_s_per_m=SIGMA_GATE, k0_sq=K0_SQ, degree=DEGREE,
        sheet=(septum, SEPTUM_TAG),
    )
    qc = _check_mode("control", rc)
    rc_ratio = qc / q0
    miss = abs(rc_ratio / RATIO_LITERAL - 1.0)
    ceiling = 1.0 / RATIO_LITERAL - 1.0
    _say(f"control Q_sheet/Q0={rc_ratio:.8f} (|delta| {abs(rc_ratio - 1):.3e}); miss from "
         f"{RATIO_LITERAL}: {100 * miss:.3f} % (ceiling {100 * ceiling:.3f} %, "
         f"{miss / Q_TOLERANCE:.2f}x band)")

    # Printed: the ratio at sigma = 1e6.
    q0_hi = solve_impedance_mesh_mode(
        merged, sigma_s_per_m=SIGMA_SCALING, k0_sq=K0_SQ, degree=DEGREE).q(0)
    qs_hi = solve_impedance_mesh_mode(
        cracked, sigma_s_per_m=SIGMA_SCALING, k0_sq=K0_SQ, degree=DEGREE).q(0)
    _say(f"sigma=1e6: Q0={q0_hi:.6e}, Q_sep={qs_hi:.6e}, ratio={qs_hi / q0_hi:.6f} "
         f"(predicted 0.771914)")
    _say(f"elapsed {time.perf_counter() - t0:.1f} s, -n {MPI.COMM_WORLD.size}")

    assert e0 < Q_TOLERANCE, (q0, q0_pozar)  # (Q0)
    assert es < Q_TOLERANCE, (qs, qs_pozar)  # (Qs)
    assert er < Q_TOLERANCE, ratio  # (R)
    assert abs(rc_ratio - 1.0) < Q_TOLERANCE, rc_ratio  # control reads 1.000
    assert miss > Q_TOLERANCE, rc_ratio  # control outside (R)'s band
    assert miss <= ceiling * (1.0 + 1e-6), (miss, ceiling)  # no larger factor
