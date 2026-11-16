"""Example (`EX-143`): the first-order ABC at oblique incidence -- a TE10 guide at 60 degrees.

A rectangular PEC guide of width a = lambda0 / (2 sin 60 deg) is driven in
TE10 at wavelength 1 m.  Its mode is a plane-wave pair travelling at 60 degrees
to the guide axis (beta = k0 cos 60 deg = pi rad/m).  The far end (facet tag 2)
is terminated three ways: Jin's angle-matched condition gamma = j k0 cos(theta)
(`first_order_abc_term(..., cos_theta=cos 60 deg)`), the plain first-order ABC
(Jin 3e 9.2 (9.44), s = 1, normal-incidence form) and PEC.  Gamma is read by the
gate's 41-point fit over x in FIT_O and set against Jin's
R(theta) = (cos(theta) - 1) / (cos(theta) + 1) = -1/3 at 60 degrees ((9.48)).

**It asserts, and it does not re-implement** (the `ANS-1` rule).  `K0`, `LAMBDA0`,
`THETA_O`, `A_O`, `L_O`, `BETA_O`, `FIT_O`, `O_BAND` and `_te10` are imported by
name from `tests/validation/test_th5_abc_reflection.py` (`TH-5` step 1a), and
`abc_first_order_reflection_te` from `fem_em_solver.core.boundary_terms`.  Two
disclosed rule (a) lifts in the gate module (the `O_BAND` hoist and the
additive `return_fields` keyword); no `src/` edit.  Asserted: (O1) |Gamma_exact|
<= O_BAND; (O2) |Gamma_ABC - R| <= O_BAND; (R) |R - (-1/3)| <= 1e-15;
(c-PEC) |Gamma_PEC - R| > O_BAND.  Records are printed beside the readings and
never asserted for equality.

Run it::

    ./run_examples.sh -e th:18 -n 2 -t 600

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_18_abc_oblique_incidence_abc_combined.xdmf`` (CG1 ``E_real``,
``E_imag``; facet tags in the same file).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.core.boundary_terms import abc_first_order_reflection_te  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_th5_abc_reflection import (  # noqa: E402
    A_O,
    BETA_O,
    FIT_O,
    K0,
    L_O,
    LAMBDA0,
    O_BAND,
    THETA_O,
    _te10,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_18_abc_oblique_incidence"

#: Gate records -- printed beside this run, a difference reported and never asserted.
#: W1 = 20260930T094044Z_TH-5-step1a2-W1.log
RECORD = {
    "cells": 66240,                                       # W1:123
    "shape": (120, 46, 2),                                # W1:123
    "a": 0.57735,                                         # W1:123
    "beta": 3.141593,                                     # W1:123
    "gamma_exact": 9.895163e-05 + 6.079623e-05j,          # W1:124
    "abs_gamma_exact": 1.161362e-04,                      # W1:124
    "gamma_abc": -3.332454e-01 + 2.468819e-04j,           # W1:125
    "abc_miss": 2.620831e-04,                             # W1:125
    "gamma_pec": -1.000000 + 0.000659j,                   # W1:129
    "pec_over_abc": 2543.7,                               # W1:129
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _cg1(msh, name):
    v = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    return v, fem.Function(v, name=name)


def _export(msh, ft, e_h, comm):
    """CG1 interpolants of Re/Im E (qualitative picture only)."""
    v_vec, e_real = _cg1(msh, "E_real")
    e_real.interpolate(fem.Expression(ufl.real(e_h), v_vec.element.interpolation_points))
    _, e_imag = _cg1(msh, "E_imag")
    e_imag.interpolate(fem.Expression(ufl.imag(e_h), v_vec.element.interpolation_points))
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_abc_combined", msh, None,
        {"E_real": e_real, "E_imag": e_imag},
        comm=comm, facet_tags=ft,
    )
    return path


def _region_tags(msh):
    """Cell tags for the figure: 1 guide interior, 2 the fit window FIT_O in x."""
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    idx = np.arange(n_owned, dtype=np.int32)
    mids = dmesh.compute_midpoints(msh, tdim, idx)
    vals = np.ones(n_owned, dtype=np.int32)
    vals[(mids[:, 0] > FIT_O[0]) & (mids[:, 0] < FIT_O[1])] = 2
    return dmesh.meshtags(msh, tdim, idx, vals)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "This example needs the complex DolfinX build: "
            "source /usr/local/bin/dolfinx-complex-mode (the runner does this "
            "automatically for the `th:` group)."
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-143 -- the first-order ABC at oblique incidence: TE10 at 60 degrees")
    _say(comm, "=" * 78)
    cos_t = float(np.cos(THETA_O))
    r_ref = float(abc_first_order_reflection_te(cos_t))
    _say(comm, f"\n[gate] O_BAND {O_BAND}, K0 {K0:.6f}, LAMBDA0 {LAMBDA0}, theta 60 deg, "
               f"L_O {L_O} m, FIT_O {FIT_O} (imported from the TH-5 step 1a gate module)")
    _say(comm, f"[closed form] R(60 deg) = (cos - 1)/(cos + 1) = {r_ref:.17f}; "
               f"|R - (-1/3)| = {abs(r_ref + 1.0 / 3.0):.3e}")

    # ---- the three solves, degree 1, lambda0/80 (the gate's rung) --------------
    t = time.perf_counter()
    g_exact, nc, shape, _ = _te10(comm, 80, 1, "exact", cos_theta=cos_t)
    _say(comm, f"\n[exact] mesh {shape} (record {RECORD['shape']}), {nc} cells "
               f"(record {RECORD['cells']}), a = {A_O:.5f} m (record {RECORD['a']}), "
               f"beta = {BETA_O:.6f} (record {RECORD['beta']}) ({time.perf_counter() - t:.1f} s)")
    _say(comm, f"[exact] angle-matched: Gamma_exact = {g_exact:.6e} "
               f"(record {RECORD['gamma_exact']:.6e})")

    t = time.perf_counter()
    g_abc, _, _, _, e_abc, msh, ft = _te10(comm, 80, 1, "abc", return_fields=True)
    _say(comm, f"[abc] first-order ABC: Gamma_ABC = {g_abc:.6e} "
               f"(record {RECORD['gamma_abc']:.6e}) ({time.perf_counter() - t:.1f} s)")

    write_setup_figure(
        msh, _region_tags(msh), FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "guide interior (PEC y, z walls; ABC at x = 1.5 m)",
                      2: "fit window 0.25 < x < 1.25 m"},
        translucent_tags=(1,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.5 * float(2.0 * LAMBDA0 / 80)),
        title="th:18 -- TE10 guide a = 0.577 m (60 deg), drive x = 0, "
              "first-order ABC at x = 1.5 m; slice z = mid",
        comm=comm,
    )
    path = _export(msh, ft, e_abc, comm)
    _say(comm, f"[output] {path}")
    del e_abc, msh, ft

    t = time.perf_counter()
    g_pec, _, _, _ = _te10(comm, 80, 1, "pec")
    _say(comm, f"[pec] Gamma_PEC = {g_pec:.6e} (record {RECORD['gamma_pec']:.6e}) "
               f"({time.perf_counter() - t:.1f} s)")

    # ---- gates ----------------------------------------------------------------
    o1 = abs(g_exact)
    o2 = abs(g_abc - r_ref)
    pm = abs(g_pec - r_ref)
    r_id = abs(r_ref - (-1.0 / 3.0))
    _say(comm, "")
    _say(comm, f"(O1) |Gamma_exact| = {o1:.6e} (record {RECORD['abs_gamma_exact']:.6e}) "
               f"<= O_BAND {O_BAND}: ASSERTED")
    _say(comm, f"(O2) |Gamma_ABC - R| = {o2:.6e} (record {RECORD['abc_miss']:.6e}) "
               f"<= O_BAND {O_BAND}: ASSERTED")
    _say(comm, f"(R)  |R - (-1/3)| = {r_id:.3e} <= 1e-15: ASSERTED")
    _say(comm, f"(c-PEC) |Gamma_PEC - R| = {pm:.4f} > O_BAND {O_BAND} "
               f"(x{pm / O_BAND:.0f}): ASSERTED")
    _say(comm, f"[printed] |Gamma_PEC - R| / |Gamma_ABC - R| = {pm / o2:.1f} "
               f"(record {RECORD['pec_over_abc']})")

    # ---- R(theta) curve (numpy, printed only) ---------------------------------
    _say(comm, "\n[R(theta) curve] Jin (9.48), first-order ABC, TE (closed form, printed only)")
    for deg in (0, 30, 45, 60, 75, 85):
        r = float(abc_first_order_reflection_te(float(np.cos(np.deg2rad(deg)))))
        _say(comm, f"  theta = {deg:2d} deg: R = {r:+.6f}")

    assert o1 <= O_BAND, g_exact
    assert o2 <= O_BAND, (g_abc, r_ref)
    assert r_id <= 1e-15, r_ref
    assert pm > O_BAND, g_pec

    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
