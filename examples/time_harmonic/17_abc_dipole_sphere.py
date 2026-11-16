"""Example (`EX-142`): where the absorbing boundary works -- a dipole in a sphere at k0 R = pi.

A Hertzian dipole's outgoing field E_out is imposed (Dirichlet) in the inner
ball r <= R_DRIVE = 0.15 m of an air ball of radius R_BALL = 0.5 m, wavelength
1 m, so k0 R = pi.  The outer sphere (facet tag 1, faceted) is terminated three
ways: the exact wave impedance (Z_s = eta0 f/g), the first-order ABC (Jin 3e
9.2 (9.44), s = 1), and PEC.  The reflection coefficient Gamma is read by a
volume least-squares projection of E onto {E_out, conj(E_out)} over
0.22 < r < 0.42 m and set against the closed form
Gamma_sph(x) = -exp(-2jx) / (2x^2 - 1 + 2jx) at x = k0 R.

**It asserts, and it does not re-implement** (the `ANS-1` rule).  `K0`, `R_BALL`,
`R_DRIVE`, `PROJ_R`, `BAND`, `BAND_CEILING`, `_gamma_sphere`, `_e_out_np`,
`_ball`, `_solve_ball` and `_DipoleProjection` are imported by name from
`tests/validation/test_th5_dipole_sphere.py` (`TH-5` step 1b); no gate module
or `src/` edit.  Asserted exactly as the gate's own four asserts: (B) BAND <=
BAND_CEILING; (E) |Gamma_E| <= BAND and (A) |Gamma_A - Gamma_sph(pi)| <= BAND at
lambda0/40; (c-PEC) |Gamma_PEC - Gamma_sph(pi)| > BAND at lambda0/20.  Records
are printed beside the readings and never asserted for equality.

Run it::

    ./run_examples.sh -e th:17 -n 2 -t 600

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_17_abc_dipole_sphere_abc_combined.xdmf`` (CG1 ``E_real``,
``E_imag``, ``E_out_real``, ``E_out_imag``; facet tags in the same file).
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

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_th5_dipole_sphere import (  # noqa: E402
    BAND,
    BAND_CEILING,
    K0,
    PROJ_R,
    R_BALL,
    R_DRIVE,
    _ball,
    _DipoleProjection,
    _e_out_np,
    _gamma_sphere,
    _solve_ball,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_17_abc_dipole_sphere"

#: Gate records -- printed beside this run, a difference reported and never asserted.
#: W2 = 20260930T185052Z_TH-5-step1b-W2.log
RECORD = {
    "cells_040": 151498,                                  # W2:120
    "cells_020": 19553,                                   # W2:112
    "gamma_e": 1.278291e-03 - 1.223502e-03j,              # W2:120
    "gamma_a": -4.673816e-02 + 1.472716e-02j,             # W2:122
    "gamma_pec": 7.825549e-01 - 6.237844e-01j,            # W2:116
    "abs_gamma_e": 1.7695e-03,                            # W2:124
    "abs_gamma_a_miss": 1.8336e-03,                       # W2:124
    "pec_miss": 1.05,                                     # W2:126 (about 175 x BAND)
    "abs_gamma_sph_pi": 0.0505957066,                     # W2:89
    "pec_over_abc": 21.61,                                # W2:126 (same rung, lambda0/20)
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _cg1(msh, name):
    v = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    return v, fem.Function(v, name=name)


def _export(msh, ft, e_h, comm):
    """CG1 interpolants of Re/Im E and of the imposed E_out (qualitative picture only)."""
    v_vec, e_real = _cg1(msh, "E_real")
    e_real.interpolate(fem.Expression(ufl.real(e_h), v_vec.element.interpolation_points))
    _, e_imag = _cg1(msh, "E_imag")
    e_imag.interpolate(fem.Expression(ufl.imag(e_h), v_vec.element.interpolation_points))
    e_out_c = fem.Function(v_vec)
    e_out_c.interpolate(_e_out_np(K0))
    _, e_out_re = _cg1(msh, "E_out_real")
    e_out_re.x.array[:] = e_out_c.x.array.real
    _, e_out_im = _cg1(msh, "E_out_imag")
    e_out_im.x.array[:] = e_out_c.x.array.imag
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_abc_combined", msh, None,
        {"E_real": e_real, "E_imag": e_imag,
         "E_out_real": e_out_re, "E_out_imag": e_out_im},
        comm=comm, facet_tags=ft,
    )
    return path


def _region_tags(msh):
    """Cell tags for the figure: 1 air, 2 imposed-drive ball, 3 projection shell."""
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    idx = np.arange(n_owned, dtype=np.int32)
    mids = dmesh.compute_midpoints(msh, tdim, idx)
    rr = np.linalg.norm(mids, axis=1)
    vals = np.ones(n_owned, dtype=np.int32)
    vals[rr <= R_DRIVE] = 2
    vals[(rr > PROJ_R[0]) & (rr < PROJ_R[1])] = 3
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
    _say(comm, "EX-142 -- where the absorbing boundary works: a dipole in a sphere, k0 R = pi")
    _say(comm, "=" * 78)
    g_sph = _gamma_sphere(np.pi)
    _say(comm, f"\n[gate] BAND {BAND}, BAND_CEILING {BAND_CEILING:.4e}, k0 {K0:.6f}, "
               f"R_BALL {R_BALL} m, R_DRIVE {R_DRIVE} m, PROJ_R {PROJ_R} "
               "(all imported from the TH-5 step 1b gate module)")
    _say(comm, f"[closed form] Gamma_sph(pi) = {g_sph:.10f}, |Gamma_sph(pi)| = {abs(g_sph):.10f} "
               f"(record {RECORD['abs_gamma_sph_pi']}; 1/sqrt(1 + 4 pi^4) = "
               f"{1.0 / np.sqrt(1.0 + 4.0 * np.pi**4):.10f})")

    # ---- rung lambda0/40: exact + ABC ----------------------------------------
    t = time.perf_counter()
    msh, ft, n40 = _ball(comm, 0.025)
    t_mesh = time.perf_counter() - t
    _say(comm, f"\n[lambda0/40] {n40} cells (record {RECORD['cells_040']}), mesh {t_mesh:.1f} s")

    write_setup_figure(
        msh, _region_tags(msh), FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "air (first-order ABC on the outer sphere)",
                      2: "imposed dipole drive r <= 0.15 m",
                      3: "projection shell 0.22 < r < 0.42 m"},
        translucent_tags=(1,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="th:17 -- air ball R = 0.5 m (k0 R = pi), imposed dipole core, "
              "first-order ABC; slice y = 0",
        comm=comm,
    )

    t = time.perf_counter()
    g_e = _DipoleProjection(_solve_ball(msh, ft, K0, "exact"), msh, K0).gamma()
    _say(comm, f"[lambda0/40] exact impedance: Gamma_E = {g_e:.6e} "
               f"(record {RECORD['gamma_e']:.6e}) ({time.perf_counter() - t:.1f} s)")
    t = time.perf_counter()
    e_abc = _solve_ball(msh, ft, K0, "abc")
    g_a = _DipoleProjection(e_abc, msh, K0).gamma()
    _say(comm, f"[lambda0/40] first-order ABC: Gamma_A = {g_a:.6e} "
               f"(record {RECORD['gamma_a']:.6e}) ({time.perf_counter() - t:.1f} s)")
    path = _export(msh, ft, e_abc, comm)
    _say(comm, f"[output] {path}")
    del e_abc, msh, ft

    # ---- control rung lambda0/20: PEC ----------------------------------------
    t = time.perf_counter()
    msh20, ft20, n20 = _ball(comm, 0.05)
    g_p = _DipoleProjection(_solve_ball(msh20, ft20, K0, "pec"), msh20, K0).gamma()
    _say(comm, f"[lambda0/20] {n20} cells (record {RECORD['cells_020']}), PEC: "
               f"Gamma_PEC = {g_p:.6e} (record {RECORD['gamma_pec']:.6e}) "
               f"({time.perf_counter() - t:.1f} s)")

    # ---- gates, the gate module's own four asserts ---------------------------
    e_ = abs(g_e)
    a_ = abs(g_a - g_sph)
    p_ = abs(g_p - g_sph)
    _say(comm, "")
    _say(comm, f"(B) BAND {BAND} <= BAND_CEILING {BAND_CEILING:.4e}: ASSERTED")
    _say(comm, f"(E) |Gamma_E| = {e_:.4e} (record {RECORD['abs_gamma_e']:.4e}) "
               f"<= BAND {BAND}: ASSERTED")
    _say(comm, f"(A) |Gamma_A - Gamma_sph(pi)| = {a_:.4e} (record "
               f"{RECORD['abs_gamma_a_miss']:.4e}) <= BAND {BAND}: ASSERTED")
    _say(comm, f"(c-PEC) |Gamma_PEC - Gamma_sph(pi)| = {p_:.4f} (record about "
               f"{RECORD['pec_miss']}) > BAND {BAND} (x{p_ / BAND:.0f}): ASSERTED")
    _say(comm, f"[printed] |Gamma_PEC|/|Gamma_A| (mixed rungs: lambda0/20 PEC over lambda0/40 ABC) = "
               f"{abs(g_p) / abs(g_a):.2f} (gate's same-rung record {RECORD['pec_over_abc']})")

    # ---- size curve (numpy, printed only) ------------------------------------
    _say(comm, "\n[size curve] |Gamma_sph(x)| at x = k0 R (closed form, printed only)")
    for x in (0.04, 0.5, 1.0, np.pi, 2 * np.pi, 3 * np.pi):
        gm = abs(_gamma_sphere(x))
        asym = f"{1.0 / (2 * x * x):.6f}" if x >= 1.0 else "n/a (not the x >> 1 regime)"
        _say(comm, f"  x = {x:8.4f}: |Gamma_sph| = {gm:.6f}   1/(2x^2) = {asym}")

    assert BAND <= BAND_CEILING, BAND
    assert e_ <= BAND, g_e
    assert a_ <= BAND, g_a
    assert p_ > BAND, g_p

    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
