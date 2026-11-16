"""Example (`EX-74`): how long a gel ball takes to settle -- the time axis.

`mat:6` solves `MAT-9`'s unperfused transient to Fo 0.2 and `mat:4` asserts the
same dimensionless problem; every example stops there, where the ball is about
half-way to steady. This one runs the same 10 cm gel ball (k = 0.6, Q = 2000
W/m^3, rho c = 4.18e6, Robin h = 10 W/m^2/K, Bi = 5/3, T_inf = T(0) = 20 C,
`MAT-9`'s fixture) out to 3 t_R and shows the settling time: the slowest mode
lambda_1 = 1.907 sets tau_1 = t_R / lambda_1^2 (about 5.3 h), and the 1 % settle
time t_99 is about a day.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Everything is
imported from ``tests/validation/test_mat9_transient_conduction_sphere.py``
(`MAT-9`) and the complex-safe norm twins of
``test_mat13_transient_perfusion_sphere.py`` (`MAT-13`); `MAT-9`'s own
``_interior_error`` / ``_rel_diff`` call ``float()`` on a complex value and are
not used. No gate edit.

Asserted: (i) CN + Rannacher, N = 40, 0.01 m, against the series, <=
``INTERIOR_BAND`` and the formatted error equals ``MAT9_RECORD``; (iii) every
run <= ``BALANCE_BAND``; structural: `MAT-13`'s by-value fixture constants equal
`MAT-9`'s. Negative control, asserted only as a band miss: the quasi-static
(steady) field standing in for T(t*) misses ``INTERIOR_BAND``.

Printed only, never gated: (iv) long backward Euler against the steady solve;
the backward-Euler and CN least-squares rates beside the records 0.984 / 2.001;
the settle curve (volume-mean theta, Fo 0.05 ... 3, CN N = 120) against the
series, the one-mode asymptote e^(-lambda_1^2 Fo) and t_99 (series and FEM).

Run it::

    ./run_examples.sh -e mat:8 -n 2 -t 180

Output (on the ball mesh): ``examples/materials/paraview_output/
materials_08_gel_ball_settling_time_settle_combined.xdmf`` (P1 ``theta_heating``
and ``theta_cooling`` at steps 0, 4, ..., 120) and ``..._settle_curve.csv``.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI
from scipy.integrate import simpson
from scipy.optimize import brentq

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_time_series,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal import (  # noqa: E402
    solve_steady_conduction,
    solve_transient_conduction,
)
from fem_em_solver.utils.analytical import (  # noqa: E402
    robin_sphere_eigenvalues,
    robin_sphere_steady_temperature,
)

from tests.validation import test_mat9_transient_conduction_sphere as mat9  # noqa: E402
from tests.validation import test_mat13_transient_perfusion_sphere as mat13  # noqa: E402
from tests.validation.test_mat9_transient_conduction_sphere import (  # noqa: E402
    BI,
    H,
    K,
    LADDER_N,
    Q,
    R,
    RHO_C,
    ROBIN_TAG,
    T_END,
    T_INF,
    T_R,
    _allsum,
    _ball,
    _rates,
    _series,
    _transient,
)
from tests.validation.test_mat13_transient_perfusion_sphere import (  # noqa: E402
    BALANCE_BAND,
    INTERIOR_BAND,
    MAT9_RECORD,
    _exact_theta,
    _interior_error,
    _rel_diff,
    _series_0,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_08_gel_ball_settling_time"
RES = 0.01
N_SETTLE = 120
FO_END = 3.0
FRAMES = tuple(range(0, N_SETTLE + 1, 4))
FO_PRINT = (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0)

#: Gate records, printed beside this run and never asserted:
#: `20260929T170037Z_MAT-9-step2b.log:115,124,130` (control, BE rate, CN rate).
RECORD = {"ctl": 0.94421, "be_rate": 0.984, "cn_rate": 2.001, "e_i": MAT9_RECORD}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _theta_of(T_h, name):
    """theta = T - T_inf as a new function in T_h's space."""
    f = fem.Function(T_h.function_space, name=name)
    f.x.array[:] = T_h.x.array - T_INF
    f.x.scatter_forward()
    return f


def _cool_of(T_ss, T_h, name):
    """theta_ss - theta_heat: exact for the linear theta-scheme (see main)."""
    f = fem.Function(T_h.function_space, name=name)
    f.x.array[:] = T_ss.x.array - T_h.x.array
    f.x.scatter_forward()
    return f


def _mean(comm, msh, f, vol):
    return float(np.real(_allsum(comm, f * ufl.dx(domain=msh)))) / vol


_RGRID = np.linspace(0.0, R, 2001)


def _series_mean(t):
    """Volume-mean theta of the series: 3/R^3 int theta r^2 dr."""
    th = _series(_RGRID, t) - T_INF
    return float(3.0 / R**3 * simpson(th * _RGRID**2, x=_RGRID))


def _first_below(ts, deficits, level):
    """First time the (monotone) deficit crosses `level`, log-linear."""
    for i in range(1, len(ts)):
        if deficits[i] <= level < deficits[i - 1]:
            a, b = np.log(deficits[i - 1]), np.log(deficits[i])
            return float(ts[i - 1] + (np.log(level) - a) / (b - a) * (ts[i] - ts[i - 1]))
    return float("nan")


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    # Structural check: MAT-13 restates MAT-9's fixture by value.
    assert (mat13.R, mat13.K, mat13.Q, mat13.H, mat13.T_INF) == (R, K, Q, H, T_INF)
    assert mat13.RHO_C == RHO_C and mat13.T_R == T_R and mat13.T_END == T_END
    assert mat13.BI == BI and mat13.ROBIN_TAG == ROBIN_TAG
    assert tuple(mat13.LADDER_N) == tuple(LADDER_N) and mat13.RES == RES
    assert mat13.QDEG == mat9.QDEG

    _say(comm, "=" * 78)
    _say(comm, "EX-74 -- how long a gel ball takes to settle: the time axis")
    _say(comm, "=" * 78)
    lam1 = float(robin_sphere_eigenvalues(BI, 1)[0])
    tau1 = T_R / lam1**2
    _say(comm, f"\n[gate] bands imported: interior {INTERIOR_BAND:.0e}, balance "
               f"{BALANCE_BAND:.0e}; Bi = {BI:.4f}, t_R = {T_R:.1f} s, t* = {T_END:.1f} s "
               f"(Fo 0.2), lambda_1 = {lam1:.6f}, tau_1 = t_R/lambda_1^2 = {tau1:.1f} s "
               f"= {tau1 / 3600:.2f} h")

    msh, cell_tags, ft = _ball(RES)
    n_cells = msh.topology.index_map(3).size_global
    _say(comm, f"[mesh] MAT-9's _ball({RES}): {n_cells} cells")

    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel ball (Robin h = 10, Bi = 5/3)"},
        translucent_tags=(1,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:8 -- 10 cm Robin gel ball (MAT-9 fixture), settling time",
        comm=comm,
    )
    one = fem.Constant(msh, default_scalar_type(1.0))
    vol = float(np.real(_allsum(comm, one * ufl.dx(domain=msh))))
    balances: dict[str, float] = {}

    # ---- (i) CN + Rannacher N = 40 at t*, against the series ---------------------
    s_i = _transient(msh, ft, 40, 0.5, 2)
    balances["(i) CN N=40"] = s_i.max_balance_residual
    theta_ex = _exact_theta(msh, _series_0)
    e_i = _interior_error(msh, s_i.temperature, theta_ex)
    _say(comm, f"\n[transient] (i) CN + Rannacher N = 40 at t*: interior rel L2 on theta "
               f"{e_i:.4e} (record {MAT9_RECORD}); string "
               f"{'equal' if f'{e_i:.4e}' == MAT9_RECORD else 'DIFFERENT'}  "
               f"ASSERTED <= {INTERIOR_BAND:.0e} and == record")

    # ---- negative control: the steady field stands in for T(t*) ---------------------
    st = solve_steady_conduction(msh, ft, k=K, q=Q, h=H, t_inf=T_INF,
                                 robin_tag=ROBIN_TAG, degree=1)
    e_ctl = _interior_error(msh, st.temperature, theta_ex)
    _say(comm, f"[control] quasi-static (steady) field for theta(t*): interior rel L2 "
               f"{e_ctl:.5f} (record {RECORD['ctl']:.5f}); ASSERTED > {INTERIOR_BAND:.0e}")

    # ---- (iv) long BE against the steady solve (printed) -----------------------------
    s_inf = solve_transient_conduction(
        msh, ft, k=K, rho_c=RHO_C, q=Q, h=H, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=5.0 * T_R, n_steps=20, theta=1.0,
    )
    balances["(iv) BE 20x5t_R"] = s_inf.max_balance_residual
    d_inf = _rel_diff(msh, s_inf.temperature, st.temperature)
    _say(comm, f"[transient] (iv) BE 20 x 5 t_R vs steady: rel diff {d_inf:.3e} (printed; "
               f"MAT-9 asserts it at 1e-8)")

    # ---- the Delta t ladder: BE and CN rates (printed beside the records) -------------
    rates = {}
    for label, theta, startup in (("BE", 1.0, 0), ("CN+Rannacher", 0.5, 2)):
        sols = []
        for n in LADDER_N:
            r = _transient(msh, ft, n, theta, startup)
            balances[f"ladder {label} N={n}"] = r.max_balance_residual
            sols.append(r.temperature)
        d = [_rel_diff(msh, sols[i], sols[i + 1]) for i in range(len(sols) - 1)]
        lsq, pair = _rates(d)
        rates[label] = lsq
        rec = RECORD["be_rate"] if label == "BE" else RECORD["cn_rate"]
        _say(comm, f"[ladder] {label}: d_N (N = 10, 20, 40, 80) = {[f'{x:.4e}' for x in d]}, "
                   f"LSQ rate {lsq:.3f} (record {rec:.3f}, diff {lsq - rec:+.3f}), "
                   f"pairwise {[f'{p:.3f}' for p in pair]}  (printed)")

    # ---- the settle curve: CN N = 120 out to 3 t_R -------------------------------------
    dt = FO_END * T_R / N_SETTLE
    steps_all = tuple(range(0, N_SETTLE + 1))
    sres = solve_transient_conduction(
        msh, ft, k=K, rho_c=RHO_C, q=Q, h=H, t_inf=T_INF, t_init=T_INF,
        robin_tag=ROBIN_TAG, dt=dt, n_steps=N_SETTLE, theta=0.5,
        startup_half_steps=2, degree=1, snapshot_steps=steps_all,
    )
    balances["settle CN N=120"] = sres.max_balance_residual
    worst = max(balances, key=balances.get)
    _say(comm, f"\n[transient] (iii) worst per-step balance {balances[worst]:.3e} on run "
               f"'{worst}' of {len(balances)}  ASSERTED <= {BALANCE_BAND:.0e}")

    steady_fem = _mean(comm, msh, _theta_of(st.temperature, "t"), vol)
    steady_ser = float(3.0 / R**3 * simpson(
        robin_sphere_steady_temperature(_RGRID, radius=R, k=K, q=Q, h=H, t_inf=T_INF)
        * _RGRID**2 - T_INF * _RGRID**2, x=_RGRID))
    fem_mean = np.array([_mean(comm, msh, _theta_of(sres.snapshots[i], "t"), vol)
                         for i in steps_all])
    t_all = np.array(steps_all) * dt
    fo_all = t_all / T_R
    _say(comm, f"[settle] steady volume-mean theta: FEM {steady_fem:.6f} K, closed form "
               f"{steady_ser:.6f} K (rel diff {abs(steady_fem / steady_ser - 1):.2e}); "
               f"CN N = {N_SETTLE}, dt = {dt:.1f} s = {dt / T_R:.4f} t_R")
    _say(comm, "[settle]    Fo     t [h]   FEM mean theta [K]   series mean [K]   "
               "FEM def   series def   one-mode def")

    ser_def_end = steady_ser - _series_mean(FO_END * T_R)
    amp = ser_def_end * np.exp(lam1**2 * FO_END)  # one-mode amplitude from Fo 3
    rows = []
    for fo in FO_PRINT:
        i = int(round(fo * N_SETTLE / FO_END))
        assert abs(fo_all[i] - fo) < 1e-9
        sm = _series_mean(fo * T_R)
        fdef = (steady_fem - fem_mean[i]) / steady_fem
        sdef = (steady_ser - sm) / steady_ser
        odef = amp / steady_ser * np.exp(-lam1**2 * fo)
        rows.append((fo, fem_mean[i], sm, fdef, sdef, odef))
        _say(comm, f"[settle] {fo:5.2f}  {fo * T_R / 3600:8.3f}   {fem_mean[i]:.6e}     "
                   f"{sm:.6e}    {fdef:.3e}   {sdef:.3e}   {odef:.3e}")
    slope_fem = np.log(((steady_fem - fem_mean[-1]) / (steady_fem - fem_mean[-21]))) / (
        fo_all[-1] - fo_all[-21])
    _say(comm, f"[settle] late-time log-slope of the FEM deficit, Fo 2.5 -> 3: {slope_fem:.4f} "
               f"(one-mode -lambda_1^2 = {-lam1**2:.4f}; printed)")

    # ---- t_99 (printed) ---------------------------------------------------------------
    fem_def = (steady_fem - fem_mean) / steady_fem
    t99_fem = _first_below(t_all, fem_def, 0.01)
    t99_ser = brentq(lambda t: (steady_ser - _series_mean(t)) / steady_ser - 0.01,
                     0.2 * T_R, 3.0 * T_R, xtol=1.0)
    t99_one = tau1 * np.log(amp / steady_ser / 0.01)
    _say(comm, f"[t99] 1 % settle time of the volume-mean theta: series {t99_ser:.0f} s "
               f"= {t99_ser / 3600:.2f} h = {t99_ser / 86400:.3f} d (Fo {t99_ser / T_R:.3f}); "
               f"FEM {t99_fem:.0f} s = {t99_fem / 3600:.2f} h (rel diff to series "
               f"{t99_fem / t99_ser - 1:+.2e}); one-mode {t99_one:.0f} s "
               f"(rel diff {t99_one / t99_ser - 1:+.2e})  (printed, not gated)")
    _say(comm, "[cooling] the cooling curve theta_ss,h - theta_heat,h is the exact cooling "
               "from the discrete steady state for the linear theta-scheme (q off, the "
               "discrete steady state is a fixed point of the scheme with q on); it is "
               "written, not solved separately.")

    # ---- assertions ---------------------------------------------------------------------------
    assert e_i <= INTERIOR_BAND, f"(i) {e_i:.4e}"
    assert f"{e_i:.4e}" == MAT9_RECORD, f"(i) string {e_i:.4e} != {MAT9_RECORD}"
    assert balances[worst] <= BALANCE_BAND, f"(iii) {worst} {balances[worst]:.3e}"
    assert e_ctl > INTERIOR_BAND, f"control {e_ctl:.4e} did not miss {INTERIOR_BAND}"

    # ---- outputs (on the ball mesh; CG1 exports are interpolants, nothing asserted) --------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    frames = [
        (step * dt, {
            "theta_heating": _theta_of(sres.snapshots[step], "theta_heating"),
            "theta_cooling": _cool_of(st.temperature, sres.snapshots[step], "theta_cooling"),
        })
        for step in FRAMES
    ]
    series_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_settle_combined", msh, None, frames, comm=comm,
    )
    table = OUTPUT_DIR / f"{BASENAME}_settle_curve.csv"
    if comm.rank == 0:
        with table.open("w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["step", "Fo", "t_s", "mean_theta_heating_K", "mean_theta_cooling_K",
                         "steady_mean_K", "series_mean_theta_heating_K"])
            for i in steps_all:
                wr.writerow([i, f"{fo_all[i]:.6f}", f"{t_all[i]:.3f}", f"{fem_mean[i]:.9e}",
                             f"{steady_fem - fem_mean[i]:.9e}", f"{steady_fem:.9e}",
                             f"{_series_mean(t_all[i]):.9e}"])
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {series_path} -- P1 theta_heating, theta_cooling at steps "
               f"{FRAMES[0]}, {FRAMES[1]}, ..., {FRAMES[-1]}")
    _say(comm, f"[output] {table}")

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
