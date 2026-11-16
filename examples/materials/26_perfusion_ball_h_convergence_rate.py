"""Example (`EX-149`): the P1 error on the perfused Robin ball falls at rate
about 2 in h against the Pennes closed form.

`mat:6` (`EX-69`) shows one rung's error on this ball. This example runs the
whole `MAT-11` gate (ii) ladder, 0.02 / 0.01 / 0.005 m, and prints the fitted
convergence rate in h of the whole-ball relative L2 error on theta = T - T_inf,
beside the interior error and the power balance on every rung.

**It asserts, and it does not re-implement** (the `ANS-1` rule). The ladder,
bands, constants and helpers are imported by name from
``tests/validation/test_mat11_perfusion_sphere.py`` (its ``_ball``, ``LADDER``
and ``_interior_indicator`` are `MAT-9`'s, re-exported there). No gate edit.

Asserted, each as the gate writes it:

* (iii) balance <= ``BALANCE_BAND`` on every rung;
* (ii) least-squares whole-ball rate >= ``RATE_FLOOR``;
* (i) interior relative L2 at 0.01 m <= ``INTERIOR_BAND``;
* (c-w0) the gate's control (a): w = 0 at 0.01 m misses ``INTERIOR_BAND``
  against the perfused closed form, and its ``perfused_power`` is exactly 0.0.

Printed only: the w = 0 whole-ball error on all three rungs and its fitted
rate (predicted a plateau, rate about 0; never asserted).

Run it::

    ./run_examples.sh -e mat:26 -n 2 -t 420

Output (finest rung): ``examples/materials/paraview_output/
materials_26_perfusion_ball_h_convergence_rate_finest_combined.xdmf`` carrying
P1 ``theta_h``, P1 ``theta_closed_form`` and DG0 ``interior``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import solve_steady_conduction  # noqa: E402

from tests.validation.test_mat11_perfusion_sphere import (  # noqa: E402
    BALANCE_BAND,
    H,
    INTERIOR_BAND,
    K,
    LADDER,
    Q,
    RATE_FLOOR,
    ROBIN_TAG,
    T_INF,
    V_BALL,
    W,
    _ball,
    _balance,
    _centre_surface,
    _errors,
    _interior_indicator,
    _max_abs_imag,
    _theta,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_26_perfusion_ball_h_convergence_rate"

#: Gate records, complex build -- `20261001T004409Z_MAT-11-cx-W2.log:117-130`.
#: Printed beside this run, a difference reported and never asserted.
RECORD = {
    "cells": (2571, 19463, 151548),
    "balance": (3.227e-15, 1.703e-15, 3.820e-15),
    "rate": 2.060,
    "pairwise": (2.068, 2.052),
    "interior": 1.8285e-03,
    "ctl": 1.1055e00,
    "centre": 5.5394,
    "surface": 3.2668,
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _rate(h, e):
    return float(np.polyfit(np.log(h), np.log(e), 1)[0])


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-149 -- P1 h-convergence rate on the perfused Robin ball (MAT-11 gate (ii))")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] imported from MAT-11: ladder {LADDER}, interior band {INTERIOR_BAND:.0e}, "
               f"balance band {BALANCE_BAND:.0e}, rate floor {RATE_FLOOR}; w = {W:g} W/m^3/K")

    rows = []
    finest = None
    for i, res in enumerate(LADDER):
        t0 = time.perf_counter()
        msh, cell_tags, ft = _ball(res)
        if res == 0.01:
            write_setup_figure(
                msh,
                cell_tags,
                FIGURE_DIR / f"{BASENAME}_setup.png",
                region_names={1: "gel ball (Pennes sink w = 240)"},
                translucent_tags=(1,),
                slice_normal=(0.0, 0.0, 1.0),
                slice_origin=(0.0, 0.0, 0.0),
                title="mat:26 -- 10 cm Robin ball with perfusion, rung 0.01 m of the h ladder",
                comm=comm,
            )
        sol = solve_steady_conduction(
            msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1,
            perfusion_w_per_m3_k=W,
        )
        sol0 = solve_steady_conduction(
            msh, ft, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1,
        )
        e_in, e_all, vol, n = _errors(msh, sol.temperature, W)
        _e0_in, e0_all, _v, _n = _errors(msh, sol0.temperature, W)
        bal = _balance(sol)
        h_eff = (vol / n) ** (1.0 / 3.0)
        centre, surf = _centre_surface(msh, sol.temperature)
        rows.append((res, n, h_eff, e_in, e_all, bal, e0_all))
        if res == 0.01:
            imag = _max_abs_imag(sol.temperature)
            _say(comm, f"[rung] max|Im theta_h| at 0.01 m = {imag:.3e} (predicted 0.0; "
                       f"ScalarType={np.dtype(default_scalar_type).name})")
        _say(comm, f"[rung] res={res:.3f} cells={n} (record {RECORD['cells'][i]}) "
                   f"h_eff={h_eff:.5e} interior_rel_L2={e_in:.4e} whole_rel_L2={e_all:.4e}")
        _say(comm, f"[rung]   balance {bal:.3e} (record {RECORD['balance'][i]:.3e}) "
                   f"ASSERTED <= {BALANCE_BAND:.0e}; source {sol.source_power:.6f} W, "
                   f"convected {sol.convected_power:.6f} W, perfused {sol.perfused_power:.6f} W "
                   f"({sol.perfused_power / sol.source_power:.4f} of the source)")
        _say(comm, f"[rung]   theta_h centre {centre:.4f} K (closed form {RECORD['centre']:.4f}), "
                   f"surface mean {surf:.4f} K (closed form {RECORD['surface']:.4f}); "
                   f"Vh/V = {vol / V_BALL:.6f}; w = 0 whole_rel_L2 = {e0_all:.4e}; "
                   f"{time.perf_counter() - t0:.1f} s")
        if i == len(LADDER) - 1:
            finest = (msh, cell_tags, sol)

    h = np.array([r[2] for r in rows])
    e = np.array([r[4] for r in rows])
    e0 = np.array([r[6] for r in rows])
    pair = [float(np.log(e[i] / e[i + 1]) / np.log(h[i] / h[i + 1])) for i in range(len(rows) - 1)]
    rate = _rate(h, e)
    rate0 = _rate(h, e0)
    interior_001 = rows[LADDER.index(0.01)][3]
    _say(comm, f"\n[rate] pairwise {LADDER[0]}->{LADDER[1]}: {pair[0]:.3f}, "
               f"{LADDER[1]}->{LADDER[2]}: {pair[1]:.3f} (records {RECORD['pairwise'][0]:.3f} / "
               f"{RECORD['pairwise'][1]:.3f}; printed)")
    _say(comm, f"[rate] (ii) least-squares whole-ball P1 L2 rate = {rate:.3f} "
               f"(record {RECORD['rate']:.3f}, diff {rate - RECORD['rate']:+.1e}) "
               f"ASSERTED >= {RATE_FLOOR}")
    _say(comm, f"[rate] (i) interior rel L2 at 0.01 m = {interior_001:.4e} "
               f"(record {RECORD['interior']:.4e}) ASSERTED <= {INTERIOR_BAND:.0e}")

    # ---- the control: the gate's own control (a) ---------------------------------
    msh_c, _ct, ft_c = _ball(0.01)
    sol_a = solve_steady_conduction(
        msh_c, ft_c, k=K, q=Q, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    e_ctl, _, _, _ = _errors(msh_c, sol_a.temperature, W)
    _say(comm, f"\n[control] (c-w0) w = 0 at 0.01 m against the perfused closed form: "
               f"interior_rel_L2 = {e_ctl:.4e} (record {RECORD['ctl']:.4e}, "
               f"{e_ctl / INTERIOR_BAND:.0f}x the band) ASSERTED > {INTERIOR_BAND:.0e}; "
               f"perfused_power = {sol_a.perfused_power!r} ASSERTED == 0.0")
    _say(comm, f"[control] w = 0 whole-ball error on the ladder: "
               + ", ".join(f"{x:.4e}" for x in e0)
               + f"; fitted rate {rate0:.3f} (PREDICTION only: a plateau, rate about 0; "
                 f"perfused rate {rate:.3f}; never asserted)")

    # ---- assertions ---------------------------------------------------------------
    for r in rows:
        assert r[5] <= BALANCE_BAND, f"(iii) balance {r[5]:.3e} at res={r[0]}"
    assert rate >= RATE_FLOOR, f"(ii) P1 rate {rate:.3f}"
    assert interior_001 <= INTERIOR_BAND, f"(i) interior {interior_001:.4e}"
    assert e_ctl > INTERIOR_BAND, f"(c-w0) control {e_ctl:.4e} did not miss"
    assert sol_a.perfused_power == 0.0, f"(c-w0) perfused_power {sol_a.perfused_power!r}"

    # ---- output: finest rung, P1 theta_h, P1 closed form, DG0 interior -------------
    msh_f, _tags_f, sol_f = finest
    Tf = sol_f.temperature
    V = Tf.function_space
    th_h = fem.Function(V, name="theta_h")
    th_h.x.array[:] = Tf.x.array - T_INF
    th_h.x.scatter_forward()
    th_cf = fem.Function(V, name="theta_closed_form")
    xyz = V.tabulate_dof_coordinates()
    th_cf.x.array[:] = np.asarray(_theta(np.linalg.norm(xyz, axis=1), W), dtype=float)
    th_cf.x.scatter_forward()
    chi = _interior_indicator(msh_f)
    chi.name = "interior"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_finest_combined", msh_f, None,
        {"theta_h": th_h, "theta_closed_form": th_cf, "interior": chi}, comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- finest rung: P1 theta_h, P1 theta_closed_form, DG0 interior")

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
