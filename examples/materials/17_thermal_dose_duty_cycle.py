"""Example (`EX-104`): thermal dose and the six-minute average on closed-form histories.

`mat:11` (`EX-85`) switches a uniform source 0/1 on `MAT-9`'s 10 cm Robin gel
ball and reads the temperature; it explicitly claims nothing about a dose or a
six-minute average. This example asks what two post-processors report on the same
history -- `MAT-20` step 1's ``cem43`` (thermal dose, exact on the piecewise-linear
interpolant) and step 2's ``running_average`` (trailing window mean) -- and checks
each against a closed form. The history is `MAT-17` case A's periodic closed form
(P = 1200 s, D = 0.5), linear ramps, and `MAT-9`'s step response. **No FEM history
is used anywhere**; the mesh exists for the setup figure and as the support of a
closed-form field map ("closed form evaluated at vertex radii -- not an FEM
history"). Nothing is assembled or solved.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Every helper and
band is imported from ``tests/validation/test_mat20_cem43_closed_form.py`` (`s1`),
``test_mat20_running_average.py`` (`s2`) and ``test_mat17_duty_cycle_closed_form.py``
(`m17cf`); the schedule ``_chi`` from ``test_mat17_duty_cycle_fem.py`` and the ball
from ``test_mat9_transient_conduction_sphere.py``. No gate edit, no ``src/`` edit.
The brentq / quad settings of (A2) and (B2) are copied from the gate with citations
(settings, not bands).

Part A recomputes the gated operator identities and prints them beside the gate
logs' records: (A1) three ramps x three samplings, (A2) `MAT-9`'s step response
re-based to 37 C at N = 2^18, (B1) integer-period window means, (B2) the window
maximum's phase. Controls (asserted): (Rc) the planted R swap, (Wc) a
non-integer-period window. Part B adds: (S1) the power's six-minute maximum is 1,
(S2) its whole-period mean is D, (J) Jensen at the centre (strict inequality,
behind the asserted premise min T > 43 C), (S1c) the power's whole-period-mean
control. Printed only: (Sc), (Pc), (Jc), the six-minute temperature excess, the
doses.

Run it::

    ./scripts/run_examples.sh -e mat:17 -n 2 -t 180

Output (under the gitignored ``paraview_output/``): ``materials_17_thermal_dose_duty_cycle_map_combined.xdmf``
(closed-form P1 map), ``..._period_combined.xdmf`` (17-frame series) and
``..._centre_history.csv``.
"""

from __future__ import annotations

import csv
import math
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI
from scipy.integrate import quad
from scipy.optimize import brentq

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_time_series,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal.dose import cem43, running_average  # noqa: E402
from fem_em_solver.utils.analytical import (  # noqa: E402
    robin_sphere_duty_cycle_temperature,
    robin_sphere_steady_temperature,
)

from tests.validation import test_mat17_duty_cycle_closed_form as m17cf  # noqa: E402
from tests.validation import test_mat20_cem43_closed_form as s1  # noqa: E402
from tests.validation import test_mat20_running_average as s2  # noqa: E402
from tests.validation.test_mat17_duty_cycle_fem import _chi  # noqa: E402
from tests.validation.test_mat9_transient_conduction_sphere import _ball  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_17_thermal_dose_duty_cycle"
RES = 0.02
LEFT = 1.0 - 1e-12  # mat:11's left-limit idiom at a switch
M_MAP = 2**10  # map samples per period
LABEL = "closed form evaluated at vertex radii -- not an FEM history"

#: Gate records (printed beside; strings, never re-derived here).
#: S1 = 20261008T171434Z_MAT-20-step1.log ; S2 = 20261008T171843Z_MAT-20-step2.log
REC_A1 = {  # name: (ref, {n: (op, rect_sum)})
    "rise 41->45 / 600 s": ("7.506522634625388e+02", {
        2: ("7.506522634625388e+02", "3.7500000000e+01"),
        7: ("7.506522634625389e+02", "5.7240835460e+02"),
        1000: ("7.506522634625388e+02", "7.4947029984e+02")}),
    "shallow 42.5->43.5 / 60 s": ("5.749545675161140e+01", {
        2: ("5.749545675161138e+01", "3.0000000000e+01"),
        7: ("5.749545675161139e+01", "5.3060441491e+01"),
        1000: ("5.749545675161138e+01", "5.7468002593e+01")}),
    "fall 45->40 / 900 s": ("9.068690858587969e+02", {
        2: ("9.068690858587970e+02", "3.6000000000e+03"),
        7: ("9.068690858587967e+02", "1.2317912341e+03"),
        1000: ("9.068690858587963e+02", "9.0866476831e+02")}),
}
REC_RC = {"rise 41->45 / 600 s": "1.3784", "shallow 42.5->43.5 / 60 s": "0.1937",
          "fall 45->40 / 900 s": "1.3982"}
REC_SC = {"rise 41->45 / 600 s": "1.2973", "shallow 42.5->43.5 / 60 s": "0.1292",
          "fall 45->40 / 900 s": "0.1096"}
REC_A2 = {"t43": "14087.0956163533", "quad": "1.708516120835631e+05",
          "op131072": "1.708516120774260e+05", "op262144": "1.708516120820288e+05",
          "rel": "8.980e-12", "pred": "8.980e-12",
          "T0": "37.000000027996", "T1": "47.237564602670"}
REC_B1 = {(1, 4096): ("5.959e-15", 4097), (1, 8192): ("8.139e-15", 8193),
          (2, 4096): ("1.032e-14", 4097), (2, 8192): ("9.011e-15", 8193)}
REC_WC = {(1, 4096): ("1.957e-03", "3.915e-03"), (1, 8192): ("1.957e-03", "3.915e-03"),
          (2, 4096): ("1.174e-03", "2.349e-03"), (2, 8192): ("1.174e-03", "2.349e-03")}
REC_B2 = {4096: ("779.999309", "779.882812", "-0.398"),
          8192: ("779.999309", "780.029297", "+0.205")}
REC_PC = {4096: "779.882812", 8192: "779.882812"}

_MISS: list[str] = []  # hard-record string mismatches (asserted empty at the end)


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _eq(got: str, rec: str, hard: bool = True) -> str:
    if got == rec:
        return "equal"
    if hard:
        _MISS.append(f"{got} != {rec}")
    return "DIFFERENT"


def _field(V, name, values):
    f = fem.Function(V, name=name)
    f.x.array[:] = values
    f.x.scatter_forward()
    return f


def _theta_block(rv, t):
    """theta(r_i, t_j) for P-periodic case A, shape (len(t), len(rv)), in <= 8192-pair blocks."""
    out = np.empty((t.size, rv.size))
    per = max(1, 8192 // t.size)
    for a in range(0, rv.size, per):
        r_blk = rv[a:a + per]
        rr, tt = np.meshgrid(r_blk, t, indexing="ij")  # (nr, nt)
        out[:, a:a + per] = s2_theta_r(rr, tt).T
    return out


def s2_theta_r(rr, tt):
    """s2._theta generalised to a radius array (same call, same constants)."""
    return robin_sphere_duty_cycle_temperature(
        rr, tt, radius=s2.R, k=s2.K, rho_c=s2.RHO_C, q=s2.Q, h=s2.H, t_inf=0.0,
        period=s2.P, duty=s2.D, periodic=True,
    )


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    # ---- structural checks: one ball, three gate modules ----------------------------
    assert (s2.P, s2.D) == m17cf.CASES["A"]
    assert s2.TH_SS0 == m17cf.TH_SS0
    assert (s2.R, s2.K, s2.RHO_C, s2.Q, s2.H) == (m17cf.R, m17cf.K, m17cf.RHO_C, m17cf.Q, m17cf.H)
    assert (s1.PHYS["radius"], s1.PHYS["k"], s1.PHYS["rho_c"], s1.PHYS["q"], s1.PHYS["h"]) == (
        m17cf.R, m17cf.K, m17cf.RHO_C, m17cf.Q, m17cf.H)
    P, D, TH0 = s2.P, s2.D, s2.TH_SS0
    BAR1, BAR2 = s1.BAR, s2.BAR

    _say(comm, "=" * 78)
    _say(comm, "EX-104 -- thermal dose and the six-minute average on closed-form histories")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] structural checks pass; bands imported: s1.BAR {BAR1:g}, s2.BAR {BAR2:g}; "
               f"case A P = {P:g} s, D = {D:g}, theta_ss(0) = {TH0:.6f} K; base {s1.T_BASE:g} C")
    tau1 = m17cf.TAU1
    _say(comm, f"[time] tau_1 = {tau1:.6e} s (predicted 1.915501e4); P/tau_1 = {P / tau1:.4f} "
               f"(predicted 0.0626)")

    # ---- mesh + setup figure (collective) --------------------------------------------------
    msh, cell_tags, _ft = _ball(RES)
    n_cells = msh.topology.index_map(3).size_global
    _say(comm, f"[mesh] MAT-9's _ball({RES}): {n_cells} cells (predicted 2571); {LABEL}")
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel ball R = 0.1 m (Robin h = 10; closed-form histories only, no FEM solve)"},
        translucent_tags=(1,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:17 -- 10 cm Robin gel ball (MAT-9 / MAT-17 fixture): CEM43 and six-minute average "
              "from closed forms",
        comm=comm,
    )

    # =============================== Part A ======================================================
    _say(comm, "\n--- Part A: the gated operator identities, recomputed ---")
    for name, (a, b_end, tau) in s1.RAMPS.items():
        ref = s1.ramp_dose(a, b_end, tau)
        ref_rec, per_n = REC_A1[name]
        for n in s1.SAMPLINGS:
            t, T = s1._ramp(a, b_end, tau, n)
            got = cem43(t, T)
            err = s1._rel(got, ref)
            rs = s1._rect_sum(t, T)
            op_rec, rs_rec = per_n[n]
            e_ref = _eq(f"{ref:.15e}", ref_rec)
            e_rs = _eq(f"{rs:.10e}", rs_rec)
            e_op = _eq(f"{got:.15e}", op_rec, hard=False)
            _say(comm, f"[A1] {name:28s} n={n:5d} ref={ref:.15e} ({e_ref}) op={got:.15e} "
                       f"(rec {op_rec}, {e_op}) rel={err:.3e} rect_sum={rs:.10e} ({e_rs})  "
                       f"ASSERTED rel <= {BAR1:g}")
            assert err <= BAR1, (name, n, err)

    # (A2): the gate's grid, brentq / quad settings copied from
    # test_mat20_cem43_closed_form.py:142-149 (settings, not bands).
    T0, T1 = s1._T_centre([0.0, s1.T_END])
    e0, e1_ = _eq(f"{T0:.12f}", REC_A2["T0"]), _eq(f"{T1:.12f}", REC_A2["T1"])
    _say(comm, f"[A2] bracket T(0)={T0:.12f} ({e0}) T(T_END={s1.T_END:g})={T1:.12f} ({e1_}) "
               f"(S1 log :113)")
    assert T0 < 43.0 < T1, "brentq bracket"
    t43 = brentq(lambda s: s1._T_centre(s)[0] - 43.0, 0.0, s1.T_END, xtol=1e-12,
                 rtol=4 * np.finfo(float).eps)

    def f_dose(s):
        T = s1._T_centre(s)[0]
        return (0.5 if T >= 43.0 else 0.25) ** (43.0 - T)

    i1, _ = quad(f_dose, 0.0, t43, epsabs=1e-13, epsrel=1e-13, limit=2000)
    i2, _ = quad(f_dose, t43, s1.T_END, epsabs=1e-13, epsrel=1e-13, limit=2000)
    ref_a2 = i1 + i2
    _say(comm, f"[A2] t43 = {t43:.10f} s ({_eq(f'{t43:.10f}', REC_A2['t43'])}, rec {REC_A2['t43']}); "
               f"quad = {ref_a2:.15e} equivalent s ({_eq(f'{ref_a2:.15e}', REC_A2['quad'])}, "
               f"rec {REC_A2['quad']})")
    k_max = 19
    t_fine = np.linspace(0.0, s1.T_END, 2**k_max + 1)
    T_fine = s1._T_centre(t_fine)
    d17 = cem43(t_fine[::4], T_fine[::4])
    d18 = cem43(t_fine[::2], T_fine[::2])
    err18 = s1._rel(d18, ref_a2)
    pred18 = abs(d18 - d17) / 3.0 / abs(d18)
    _say(comm, f"[A2] N = 131072: op = {d17:.15e} equivalent s ({_eq(f'{d17:.15e}', REC_A2['op131072'])}, "
               f"rec {REC_A2['op131072']}, S1 :122)")
    _say(comm, f"[A2] N = 262144: op = {d18:.15e} equivalent s ({_eq(f'{d18:.15e}', REC_A2['op262144'])}, "
               f"rec {REC_A2['op262144']}); rel {err18:.3e} (rec {REC_A2['rel']}, "
               f"{_eq(f'{err18:.3e}', REC_A2['rel'])}) pred {pred18:.3e} (rec {REC_A2['pred']}, "
               f"{_eq(f'{pred18:.3e}', REC_A2['pred'])})  ASSERTED rel <= {BAR1:g} "
               f"(N = 2^18 is the gate's recorded selection; no selection loop)")
    assert err18 <= BAR1, err18
    dose_cont_min = d18 / 60.0

    # (B1)
    target = D * TH0
    for k in (1, 2):
        for M in s2.M_SEQ:
            t, v = s2._series(M, k + 1)
            means, _, _ = running_average(t, v, window=k * P)
            ok = ~np.isnan(means)
            res = float(np.max(np.abs(means[ok] - target))) / target
            rec_res, rec_n = REC_B1[(k, M)]
            _say(comm, f"[B1] k={k} M={M}: {int(ok.sum())} windows (rec {rec_n}), max |mean - D*theta_ss|"
                       f"/(D*theta_ss) = {res:.3e} (rec {rec_res}, {_eq(f'{res:.3e}', rec_res)})  "
                       f"ASSERTED <= {BAR2:g}")
            assert int(ok.sum()) == rec_n
            assert res <= BAR2, (k, M, res)
            # (Wc) asserted: W = kP + P/2
            Wc = k * P + 0.5 * P
            t2, v2 = s2._series(M, k + 2)
            mc, _, _ = running_average(t2, v2, window=Wc)
            okc = ~np.isnan(mc)
            miss = float(np.max(np.abs(mc[okc] - target))) / target
            pred = 0.5 * (v.max() - v.min()) * (0.5 * P / Wc) / target
            rec_miss, rec_pred = REC_WC[(k, M)]
            _say(comm, f"[Wc] k={k} M={M}: W = {Wc:g} s miss {miss:.3e} (rec {rec_miss}, "
                       f"{_eq(f'{miss:.3e}', rec_miss)}); O(ripple) ceiling {pred:.3e} (rec {rec_pred}, "
                       f"{_eq(f'{pred:.3e}', rec_pred)})  ASSERTED miss > {BAR2:g}")
            assert miss > BAR2, (k, M, miss)

    # (B2): bracket and brentq copied from test_mat20_running_average.py:106-109.
    g = lambda s: float(s2._theta(np.array([s]))[0] - s2._theta(np.array([s - s2.W_B2]))[0])  # noqa: E731
    a_b, b_b = D * P + 1.0, D * P + s2.W_B2 - 1.0
    assert g(a_b) > 0.0 > g(b_b), (g(a_b), g(b_b))
    root = brentq(g, a_b, b_b, xtol=1e-12, rtol=4 * np.finfo(float).eps)
    for M in s2.M_SEQ:
        dt = P / M
        t, v = s2._series(M, 2)
        means, _, _ = running_average(t, v, window=s2.W_B2)
        sel = t >= P
        i = int(np.nanargmax(np.where(sel, means, -np.inf)))
        t_max = t[i] - P
        off = t_max - root
        rr_, rt_, ro_ = REC_B2[M]
        _say(comm, f"[B2] M={M}: brentq root {root:.6f} s ({_eq(f'{root:.6f}', rr_)}, rec {rr_}), "
                   f"t_max {t_max:.6f} s ({_eq(f'{t_max:.6f}', rt_)}, rec {rt_}), offset "
                   f"{off / dt:+.3f} samples ({_eq(f'{off / dt:+.3f}', ro_)}, rec {ro_})  "
                   f"ASSERTED |offset| <= dt = {dt:.6f} s")
        assert abs(off) <= dt
        mp, _, _ = running_average(t, v, window=s2.W_B2 - dt)
        ip = int(np.nanargmax(np.where(sel, mp, -np.inf)))
        _say(comm, f"[Pc] M={M}: W - dt -> t_max {t[ip] - P:.6f} s (rec {REC_PC[M]}, "
                   f"{_eq(f'{t[ip] - P:.6f}', REC_PC[M], hard=False)}; printed only)")

    # (Rc) asserted, (Sc) printed
    for name, (a, b_end, tau) in s1.RAMPS.items():
        ref = s1.ramp_dose(a, b_end, tau)
        misses = []
        for n in s1.SAMPLINGS:
            t, T = s1._ramp(a, b_end, tau, n)
            miss = s1._rel(cem43(t, T, r_above=0.25, r_below=0.5), ref)
            misses.append(miss)
            assert miss > BAR1, (name, n, miss)
        e_rc = [_eq(f"{m:.4f}", REC_RC[name]) for m in misses]
        r_ = 0.5 if a >= 43.0 else 0.25
        sc = s1._rel(s1._side(a, b_end, tau, r_), ref)
        _say(comm, f"[Rc] {name:28s} swap miss at n = 2 / 7 / 1000: "
                   f"{misses[0]:.4f} / {misses[1]:.4f} / {misses[2]:.4f} (rec {REC_RC[name]}, {e_rc}); "
                   f"ASSERTED > {BAR1:g}; (Sc) no split {sc:.4f} (rec {REC_SC[name]}; printed only)")

    # =============================== Part B ======================================================
    _say(comm, "\n--- Part B: the gel-phantom exposure story ---")
    M0 = s2.M_SEQ[0]
    t, th = s2._series(M0, 2)  # centre, two periods, switches on the grid
    sec = t >= P  # second period
    peak = float(m17cf._modal(0.0, D * P * LEFT, P, D, periodic=True))
    trough = float(m17cf._modal(0.0, P * LEFT, P, D, periodic=True))
    _say(comm, f"[centre] periodic peak / trough / mean: {peak:.6f} / {trough:.6f} / {target:.6f} K "
               f"(predicted 6.182881 / 6.039341 / 6.111111); as C on the {s1.T_BASE:g} base: "
               f"{s1.T_BASE + peak:.6f} / {s1.T_BASE + trough:.6f} / {s1.T_BASE + target:.6f}")

    # (S1), (S2), (S1c): the power
    chi = _chi(P, D)
    x = np.fromiter((chi(float(s)) for s in t), dtype=float, count=t.size)
    pm, pmx, _ = running_average(t, x, window=s2.W_B2)
    s1_dev = abs(pmx - 1.0)
    _say(comm, f"[S1] six-minute (W = {s2.W_B2:g} s) maximum of the power schedule chi = {pmx:.15f}; "
               f"|max - 1| = {s1_dev:.3e}  ASSERTED <= {BAR2:g}")
    assert s1_dev <= BAR2, s1_dev
    pw, _, _ = running_average(t, x, window=P)
    okw = ~np.isnan(pw)
    s2_dev = float(np.max(np.abs(pw[okw] - D))) / D
    _say(comm, f"[S2] whole-period (W = P) window mean of chi: {int(okw.sum())} windows, max |mean - D|/D = "
               f"{s2_dev:.3e}  ASSERTED <= {BAR2:g}")
    assert s2_dev <= BAR2, s2_dev
    _say(comm, f"[power] six-minute max {pmx:.6f} vs period mean {float(np.nanmean(pw[okw])):.6f}: the "
               f"six-minute average reads {pmx / D:.3f} x the period average (predicted 2.000)")
    pc, _, _ = running_average(t, x, window=1.5 * P)
    okc = ~np.isnan(pc)
    s1c = float(np.max(np.abs(pc[okc] - D))) / D
    _say(comm, f"[S1c] control: W = 1.5 P on chi: max |mean - D|/D = {s1c:.6f} (predicted ceiling 1/3 = "
               f"{1 / 3:.6f}; no factor asserted)  ASSERTED > {BAR2:g}")
    assert s1c > BAR2, s1c

    # six-minute temperature at the centre (printed)
    mth, _, _ = running_average(t, th, window=s2.W_B2)
    i = int(np.nanargmax(np.where(sec, mth, -np.inf)))
    th6 = float(mth[i])
    _say(comm, f"[temperature] six-minute max of theta(0) = {th6:.4f} K (predicted 6.1614 [6.150, 6.172]) at "
               f"t = {t[i] - P:.3f} s (predicted ~780); {100 * (th6 / target - 1):+.2f} % over D*theta_ss(0) "
               f"(predicted +0.82 [+0.6, +1.0] %)")

    # (J) Jensen at the centre
    T_c = s1.T_BASE + th[sec]
    t_c = t[sec]
    t_min = float(T_c.min())
    _say(comm, f"[J] premise: min T over the periodic second period = {t_min:.6f} C (predicted 43.039341)  "
               f"ASSERTED > 43")
    assert t_min > 43.0, t_min
    Tbar = s1.T_BASE + target
    dose_duty = cem43(t_c, T_c)
    dose_stand = cem43(np.array([0.0, P]), np.array([Tbar, Tbar]))
    margin = dose_duty / dose_stand - 1.0
    _say(comm, f"[J] centre, 37 C base, per period: duty {dose_duty:.2f} equivalent s = {dose_duty / 60:.3f} min "
               f"(predicted 1296.6 [1296.3, 1297.2] s); stand-in {dose_stand:.2f} s = {dose_stand / 60:.3f} min "
               f"(predicted 1296.07 s = 21.601 min; arithmetic P*2^(D*theta_ss(0) - 6) = "
               f"{P * 2.0 ** (target - 6.0):.2f} s)")
    _say(comm, f"[J] Jensen margin duty / stand-in - 1 = {margin:+.3e} (predicted [2e-4, 8e-4]); dose per "
               f"minute of elapsed time {dose_duty / P:.4f} equivalent min/min (predicted 1.080)  "
               f"ASSERTED duty > stand-in")
    assert dose_duty > dose_stand, (dose_duty, dose_stand)

    # surface and bench gel (printed)
    th_s = s2_theta_r(np.full(t.size, s2.R), t)
    dose_surf = cem43(t[sec], s1.T_BASE + th_s[sec])
    _say(comm, f"[surface] r = R, 37 C base, per period: {dose_surf:.3f} equivalent s (predicted 30 [20, 45]); "
               f"D*theta_ss(R) = "
               f"{D * float(robin_sphere_steady_temperature(s2.R, radius=s2.R, k=s2.K, q=s2.Q, h=s2.H)):.4f} K")
    dose_bench = cem43(t_c, m17cf.T_INF + th[sec])
    _say(comm, f"[bench] centre, {m17cf.T_INF:g} C base, per period: {dose_bench:.3e} equivalent s "
               f"(predicted 8e-8 [1e-8, 1e-6]); CEM43 is a tissue metric, the 37 C base is the gate's "
               f"re-basing, not a claim about any phantom")
    _say(comm, f"[continuous] heating from rest (A2's history) over {s1.T_END / 3600:g} h: {d18:.6e} equivalent s "
               f"= {dose_cont_min:.1f} min (predicted 1.708516e5 s = 2847.5 min); the ten-hour exposure is "
               f"the gate's T_END, not a protocol")

    # ---- the closed-form field map (not an FEM history) -------------------------------------
    V1 = fem.functionspace(msh, ("Lagrange", 1))
    n_loc = V1.dofmap.index_map.size_local
    n_tot = n_loc + V1.dofmap.index_map.num_ghosts
    coords = V1.tabulate_dof_coordinates()[:n_loc]
    rv = np.minimum(np.linalg.norm(coords, axis=1), s2.R)
    t_map = np.arange(2 * M_MAP + 1) * (P / M_MAP)
    sec_m = t_map >= P
    t0 = time.perf_counter()
    th_map = _theta_block(rv, t_map)  # (nt, n_loc)
    wall_map = time.perf_counter() - t0
    th_ss_r = np.asarray(robin_sphere_steady_temperature(rv, radius=s2.R, k=s2.K, q=s2.Q, h=s2.H))
    mean_k = D * th_ss_r
    wm, _, _ = running_average(t_map, th_map, window=s2.W_B2)
    wm_max = np.max(np.where(sec_m[:, None], wm, -np.inf), axis=0)
    dose_d = cem43(t_map[sec_m], s1.T_BASE + th_map[sec_m]) / 60.0  # equivalent min per period
    Tbar_v = s1.T_BASE + mean_k
    dose_s = cem43(np.array([0.0, P]), np.stack([Tbar_v, Tbar_v])) / 60.0
    jensen = dose_d / dose_s - 1.0
    above = (Tbar_v >= 43.0).astype(float)
    six_pct = 100.0 * (wm_max / mean_k - 1.0)
    jmin_loc = float(jensen.min()) if n_loc else float("inf")
    r_at_loc = float(rv[int(np.argmin(jensen))]) if n_loc else float("nan")
    jmin, r_jmin = min(comm.allgather((jmin_loc, r_at_loc)))
    jmax = comm.allreduce(float(jensen.max()) if n_loc else -np.inf, op=MPI.MAX)
    ncell_above = comm.allreduce(int(above.sum()), op=MPI.SUM)
    _say(comm, f"[map] {LABEL}: {n_loc} owned P1 dofs on rank 0 (map {wall_map:.1f} s); "
               f"vertices with mean T >= 43 C: {ncell_above}; r = 0.0200 m is D*theta_ss = 6 K")
    _say(comm, f"[Jc] minimum per-vertex Jensen excess (duty / stand-in - 1) = {jmin:+.3e} at r = {r_jmin:.4f} m "
               f"(predicted negative, in [-3e-2, -1e-4], near r = 0.0200 m); maximum {jmax:+.3e}; printed only")

    # ---- centre history CSV (second period, M = M_SEQ[0]) ---------------------------------------
    mth_c = mth[sec]
    x_c = x[sec]
    pm_c = pm[sec]
    seg = np.array([cem43(t_c[j:j + 2], T_c[j:j + 2]) for j in range(t_c.size - 1)])
    cum = np.concatenate([[0.0], np.cumsum(seg)]) / 60.0
    _say(comm, f"[csv] cumulative CEM43 at the end of the second period {cum[-1] * 60:.2f} equivalent s "
               f"(operator over the whole period {dose_duty:.2f}); {t_c.size} rows")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    if comm.rank == 0:
        with (OUTPUT_DIR / f"{BASENAME}_centre_history.csv").open("w", newline="") as fh:
            wr = csv.writer(fh)
            wr.writerow(["t_s", "theta_K", "T37_C", "chi", "chi_six_min_mean",
                         "theta_six_min_mean_K", "cumulative_cem43_min"])
            for j in range(t_c.size):
                wr.writerow([f"{t_c[j]:.6f}", f"{th[sec][j]:.9e}", f"{T_c[j]:.9e}", f"{x_c[j]:.1f}",
                             f"{pm_c[j]:.9e}", f"{mth_c[j]:.9e}", f"{cum[j]:.9e}"])

    # ---- XDMF outputs (collective; closed-form P1 maps, nothing asserted) ---------------------------
    frames = []
    for j in range(17):
        idx = M_MAP + 64 * j  # t = P + 75 j s, on the grid
        assert abs(t_map[idx] - (P + 75.0 * j)) <= 1e-9
        frames.append((float(t_map[idx]), {
            "T_C": _field(V1, "T_C", np.concatenate([s1.T_BASE + th_map[idx],
                                                     np.zeros(n_tot - n_loc)])),
            "T_six_min_mean_C": _field(V1, "T_six_min_mean_C", np.concatenate(
                [s1.T_BASE + wm[idx], np.zeros(n_tot - n_loc)])),
        }))
    series_path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_period_combined", msh, None, frames, comm=comm,
    )

    def _pad(a):
        return np.concatenate([a, np.zeros(n_tot - n_loc)])

    fields = {
        "theta_mean_K": _field(V1, "theta_mean_K", _pad(mean_k)),
        "T_six_min_max_C": _field(V1, "T_six_min_max_C", _pad(s1.T_BASE + wm_max)),
        "six_min_excess_pct": _field(V1, "six_min_excess_pct", _pad(six_pct)),
        "cem43_duty_per_period_min": _field(V1, "cem43_duty_per_period_min", _pad(dose_d)),
        "cem43_stand_in_per_period_min": _field(V1, "cem43_stand_in_per_period_min", _pad(dose_s)),
        "jensen_excess": _field(V1, "jensen_excess", _pad(jensen)),
        "above_43_mean": _field(V1, "above_43_mean", _pad(above)),
    }
    map_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_map_combined", msh, None, fields, comm=comm,
    )
    comm.Barrier()
    if comm.rank == 0:  # provenance label inside the XDMF
        for pth in (map_path, series_path):
            txt = Path(pth).read_text(encoding="utf-8")
            assert "<Domain>" in txt, pth
            note = (f'<Information Name="provenance" Value="{LABEL}; case A periodic closed form, '
                    f'37 C base, no FEM solve"/>')
            Path(pth).write_text(txt.replace("<Domain>", "<Domain>\n    " + note, 1), encoding="utf-8")
    comm.Barrier()
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {map_path} -- 7 P1 fields, {LABEL}")
    _say(comm, f"[output] {series_path} -- T_C and T_six_min_mean_C, 17 frames, t = P + 75 j s")
    _say(comm, f"[output] {OUTPUT_DIR / (BASENAME + '_centre_history.csv')} ({t_c.size} rows)")

    # ---- every hard record reproduced to its printed digits -------------------------------------------
    assert not _MISS, f"example/test divergence: {_MISS[:6]}"
    _say(comm, f"\n[records] every hard record (A1 ref / rect_sum, A2, B1, Wc, B2, Rc) equals its gate-log digits")
    _say(comm, f"[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
