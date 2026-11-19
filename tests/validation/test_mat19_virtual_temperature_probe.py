"""`MAT-19` step 1 — the virtual temperature probe on `MAT-12`'s sphere (§9 item 243).

Brief: ``docs/planning/chunks/MAT-19.md`` (§9 item 243).  `MAT-12`'s fixture
(EM solve of `MAT-4`'s lossy sphere, tag-1 submesh, Robin h = 100, CN +
Rannacher N = 40 to t* = 0.2 t_R) with every full step recorded; three probes
of radius a = 2 mm (centres 0, 3.5 mm on x, 6.5 mm off-axis) read the
ball-averaged temperature through ``thermal.probe`` — an owned-cell integral,
no point owner (the submesh owner rule raises, known-issues).

Comparand: the exact ball mean of the Robin eigen-series,
θ̄(d, a, t) = qR/(3h) + q(R² − d² − 3a²/5)/(6k)
             − Σ Bₙ e^{−λₙ²αt/R²} Xₙ(d)·3[sin κₙa − κₙa cos κₙa]/(κₙa)³.

Asserted: (N1) θ̄ vs radial shell-cap quadrature of the series ≤ 1e-10·θ_ss(0);
(M) |cⱼ| + a ≤ 0.9R; (A0) 19 496 submesh cells, snapshot 0 reads 0.0;
(T) q_cf series within ``P_BAND`` 1 % of θ̄ at 3 probes × 40 times;
(W) at ``-n 4`` / ``-n 1`` the q_cf series is within 1e-9 relative of the
``-n 2`` npz and ``cells`` equal; (C1) the 1 mm inward-displaced probe misses
1 % against the nominal comparand (max over n).  Printed: (G), q_h series,
(C2), the displaced probe against its own comparand, imag ratio.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal import ProbeBall, probe_series
from fem_em_solver.thermal.probe import _lattice_points
from fem_em_solver.utils.analytical import (
    robin_sphere_eigenvalues,
    robin_sphere_transient_coefficients,
    robin_sphere_transient_temperature,
)
from fem_em_solver.utils.constants import EPSILON_0
from tests.complex_mode import complex_only
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
    _solve_lossy_sphere,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    A_BAND,
    P_BAND,
    G,
    H,
    K,
    R,
    _cell_average_dg0,
    _interior_indicator,
    _interior_rel,
    _robin_tags,
    _say,
    _sphere_submesh,
    _to_sub,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    N_MAIN,
    RHO_C,
    T_STAR,
    _cn,
    _exact_theta,
)

A_PROBE = 2.0e-3  # 2 h_in = 0.2 R, a fixture parameter, not a real probe size
CENTRES = [
    (0.0, 0.0, 0.0),
    (3.5e-3, 0.0, 0.0),
    tuple(6.5e-3 * np.ones(3) / np.sqrt(3.0)),
]
DISPLACED = tuple((5.5 / 6.5) * np.asarray(CENTRES[2]))  # (C1): 1 mm radially inward
N_SUB = 19496  # MAT-12-W1.log:159
N1_BAND = 1.0e-10  # ·θ_ss(0); MAT-17 step 1 (P3), test_mat17_duty_cycle_closed_form.py:206–207
WIDTH_BAND = 1.0e-9  # MAT-19 row, PROJECT_PLAN.md:2512
N_TERMS = 400
DT = T_STAR / N_MAIN
NPZ = "/workspace/logs/mat19/step1-n2.npz"
_CACHE: dict = {}


def _q_cf():
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    return G * RHO_KG_M3 * cf["sar_w_per_kg"]


def _point(r, t, q):
    return robin_sphere_transient_temperature(
        np.atleast_1d(np.asarray(r, dtype=float)),
        t,
        radius=R,
        k=K,
        rho_c=RHO_C,
        q=q,
        h=H,
        t_inf=0.0,
        t_init=0.0,
        n_terms=N_TERMS,
    )


def _j1fac(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    small = x < 1.0e-4
    xs = x[~small]
    out[~small] = 3.0 * (np.sin(xs) - xs * np.cos(xs)) / xs**3
    out[small] = 1.0 - x[small] ** 2 / 10.0
    return out


def _ball_mean(d, a, t, q):
    """Exact mean of the series over B(c, a), |c| = d (Helmholtz mean value)."""
    if "lam" not in _CACHE:
        _CACHE["lam"] = robin_sphere_eigenvalues(H * R / K, N_TERMS)
    lam = _CACHE["lam"]
    _, B = robin_sphere_transient_coefficients(lam, radius=R, k=K, q=q, h=H)
    kap = lam / R
    alpha = K / RHO_C
    X = np.ones_like(kap) if d == 0.0 else np.sin(kap * d) / (kap * d)
    steady = q * R / (3.0 * H) + q * (R**2 - d**2 - 3.0 * a**2 / 5.0) / (6.0 * K)
    return steady - float(np.sum(B * np.exp(-(lam**2) * alpha * t / R**2) * X * _j1fac(kap * a)))


def _ball_quad(d, a, t, q, n=400):
    """Radial Gauss–Legendre with the shell-cap weight, normalised by 4πa³/3."""
    xg, wg = np.polynomial.legendre.leggauss(n)

    def seg(lo, hi):
        r = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
        return r, 0.5 * (hi - lo) * wg

    rs, ws = [], []
    if d == 0.0:
        r, w = seg(0.0, a)
        rs.append(r)
        ws.append(w * 4.0 * np.pi * r**2)
    else:
        if d < a:
            r, w = seg(0.0, a - d)
            rs.append(r)
            ws.append(w * 4.0 * np.pi * r**2)
        r, w = seg(abs(d - a), d + a)
        rs.append(r)
        ws.append(w * np.pi * r * (a**2 - (r - d) ** 2) / d)
    r = np.concatenate(rs)
    w = np.concatenate(ws)
    return float(np.sum(w * _point(r, t, q))) / (4.0 * np.pi * a**3 / 3.0)


def _theta_ss0(q):
    return q * R / (3.0 * H) + q * R**2 / (6.0 * K)


@complex_only
def test_mat19_n1_ball_mean_against_quadrature():
    q = _q_cf()
    ss0 = _theta_ss0(q)
    worst, worst_lim = 0.0, 0.0
    for c in CENTRES + [DISPLACED]:
        d = float(np.linalg.norm(c))
        for n in range(1, N_MAIN + 1):
            t = n * DT
            worst = max(worst, abs(_ball_mean(d, A_PROBE, t, q) - _ball_quad(d, A_PROBE, t, q)))
            worst_lim = max(
                worst_lim, abs(_ball_mean(d, 1.0e-12, t, q) - float(_point(d, t, q)[0]))
            )
    _say(
        f"MAT19 (N1) closed-form ball mean vs shell-cap quadrature: max |diff|/theta_ss(0) = "
        f"{worst / ss0:.3e}; a->0 vs point {worst_lim / ss0:.3e}  (band {N1_BAND:.0e}, "
        f"predicted <= 1e-13); theta_ss(0) = {ss0:.7f} K"
    )
    assert worst / ss0 <= N1_BAND
    assert worst_lim / ss0 <= N1_BAND


def _build_run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    em = _solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)
    msh, cell_tags = em["mesh"], em["cell_tags"]
    E, sigma = em["e_complex"], em["sigma_field"]
    t_em = time.perf_counter() - t0
    q_cf = _q_cf()

    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    ft = _robin_tags(sub)
    n_sub = sub.topology.index_map(3).size_global
    Wp, q_par = _cell_average_dg0(msh, cell_tags, G * 0.5 * sigma * ufl.inner(E, E))
    q_h = _to_sub(Wp, q_par, sub, parent_ids)

    t1 = time.perf_counter()
    res_cf = _cn(sub, ft, q_cf, snapshot_steps=range(N_MAIN + 1))
    res_h = _cn(sub, ft, q_h, snapshot_steps=range(N_MAIN + 1))
    t_cn = time.perf_counter() - t1
    chi = _interior_indicator(sub)
    theta = _exact_theta(sub, T_STAR, q_cf)
    err_P = _interior_rel(sub, chi, res_cf.temperature - theta, theta)

    probes = [ProbeBall(c, A_PROBE) for c in CENTRES] + [ProbeBall(DISPLACED, A_PROBE)]
    t2 = time.perf_counter()
    s_cf = probe_series(res_cf.snapshots, DT, probes)
    s_h = probe_series(res_h.snapshots, DT, probes)
    t_probe = time.perf_counter() - t2

    # (G): the same lattice, closed form at the lattice points, at t*.
    V = res_cf.temperature.function_space
    owned = np.arange(sub.topology.index_map(3).size_local, dtype=np.int32)
    geo = []
    for p in probes:
        pts, w = _lattice_points(V, owned, p)
        r = np.sqrt(np.sum(pts**2, axis=1))
        val = 0.0
        for s in range(0, r.size, 20000):
            val += float(np.sum(w[s : s + 20000] * _point(r[s : s + 20000], T_STAR, q_cf)))
        tot = comm.allreduce(np.array([val, float(np.sum(w))]), op=MPI.SUM)
        geo.append(tot[0] / tot[1])
    # Contributing cells per probe (global), for (W)'s "cells".
    cells = [n_sub]
    for p in probes:
        pts, _w = _lattice_points(V, owned, p)
        cells.append(int(pts.shape[0]))
    cells = np.array([cells[0]] + list(comm.allreduce(np.array(cells[1:]), op=MPI.SUM)))
    return dict(
        q_cf=q_cf,
        n_sub=n_sub,
        err_P=err_P,
        s_cf=s_cf,
        s_h=s_h,
        geo=np.array(geo),
        cells=cells,
        t_em=t_em,
        t_cn=t_cn,
        t_probe=t_probe,
        t_total=time.perf_counter() - t0,
        comm=comm,
        sub=sub,
        snap_cf=res_cf.snapshots,
        snap_h=res_h.snapshots,
        probes=probes,
    )


@pytest.fixture(scope="module")
def run():
    return _build_run()


def _check_probe_series(run):
    comm = run["comm"]
    q = run["q_cf"]
    s_cf, s_h = run["s_cf"], run["s_h"]
    _say(
        f"MAT19 q_cf = {q:.6f} W/m^3  dt = {DT:.4f} s  t* = {T_STAR:.4f} s  "
        f"(P_t) interior rel L2 at t* = {run['err_P']:.4e} (MAT-12 1.6845e-03)"
    )
    _say(
        f"MAT19 (A0) submesh cells = {run['n_sub']} (record {N_SUB}); probe V_h/V - 1 = "
        + ", ".join(f"{v:+.3e}" for v in s_cf.volumes / s_cf.nominal_volumes - 1.0)
        + " (predicted |.| <= 2e-3); lattice points inside per probe = "
        + str(list(run["cells"][1:]))
    )
    m_ok = all(np.linalg.norm(c) + A_PROBE <= 0.9 * R for c in CENTRES + [DISPLACED])
    assert list(s_cf.steps) == list(range(N_MAIN + 1))

    d = [float(np.linalg.norm(c)) for c in CENTRES]
    ref = np.array([[_ball_mean(dj, A_PROBE, n * DT, q) for dj in d] for n in range(N_MAIN + 1)])
    pts = np.array([[float(_point(dj, n * DT, q)[0]) for dj in d] for n in range(N_MAIN + 1)])
    P = s_cf.values[:, :3]
    relT = np.abs(P[1:] / ref[1:] - 1.0)
    _say("MAT19 (T) q_cf probe vs ball-mean closed form, rel error (n: c0 c1 c2):")
    for n in range(1, N_MAIN + 1):
        _say(
            f"MAT19   n={n:2d} t={n * DT:9.4f}  P = {P[n, 0]:.6e} {P[n, 1]:.6e} {P[n, 2]:.6e}"
            f"  ref = {ref[n, 0]:.6e} {ref[n, 1]:.6e} {ref[n, 2]:.6e}  rel = "
            f"{relT[n - 1, 0]:.3e} {relT[n - 1, 1]:.3e} {relT[n - 1, 2]:.3e}"
        )
    worstT = float(np.max(relT))
    _say(
        f"MAT19 (T) max over 3 probes x 40 times = {worstT:.4e}  (band {P_BAND:.0e}, "
        f"predicted <= 3e-3)"
    )
    for j in range(3):
        g = run["geo"][j] / ref[N_MAIN, j] - 1.0
        tt = P[N_MAIN, j] / ref[N_MAIN, j] - 1.0
        _say(
            f"MAT19 (G) probe {j} at t*: geometry {g:+.4e}  total {tt:+.4e}  FEM (T - G) "
            f"{tt - g:+.4e}  (G predicted |.| <= 1e-3)"
        )
    relH = np.abs(s_h.values[1:, :3] / ref[1:] - 1.0)
    _say(
        "MAT19 q_h series vs closed form at t*: "
        + " ".join(f"{s_h.values[N_MAIN, j] / ref[N_MAIN, j] - 1:+.4e}" for j in range(3))
        + f"  max over (j, n) {np.max(relH):.4e}  (printed beside A_BAND {A_BAND}; predicted "
        f"≈ +3.4–3.6 % at t*)"
    )
    _say(f"MAT19 max_imag_ratio: q_cf {s_cf.max_imag_ratio:.3e}  q_h {s_h.max_imag_ratio:.3e}")

    # (C1) displaced probe against the nominal comparand.
    Pd = s_cf.values[:, 3]
    missC1 = np.abs(Pd[1:] / ref[1:, 2] - 1.0)
    dd = float(np.linalg.norm(DISPLACED))
    own = np.array([_ball_mean(dd, A_PROBE, n * DT, q) for n in range(1, N_MAIN + 1)])
    ownrel = np.abs(Pd[1:] / own - 1.0)
    _say(
        "MAT19 (C1) displaced probe vs nominal theta_bar(6.5 mm): miss at n = 1/5/6/10/20/40 = "
        + " ".join(f"{missC1[n - 1]:.3e}" for n in (1, 5, 6, 10, 20, 40))
        + " (predicted 1.28e-4 8.15e-3 1.089e-2 2.116e-2 3.887e-2 5.466e-2); max "
        f"{np.max(missC1):.4e} (asserted > {P_BAND:.0e}); vs own theta_bar(5.5 mm) max "
        f"{np.max(ownrel):.4e} (printed)"
    )
    # (C2) point comparand, printed.
    c2 = P[N_MAIN] / pts[N_MAIN] - 1.0
    c2n20 = P[20, 2] / pts[20, 2] - 1.0
    _say(
        f"MAT19 (C2) ball series vs point closed form at t*: {c2[0]:+.4e} {c2[1]:+.4e} "
        f"{c2[2]:+.4e} (predicted +0.732 % +0.872 % +1.275 %); c2 at n=20 {c2n20:+.4e} "
        f"(predicted 1.14 %)"
    )

    # Width artefact / (W).
    w_cf = w_h = None
    if comm.size == 2:
        if comm.rank == 0:
            os.makedirs(os.path.dirname(NPZ), exist_ok=True)
            np.savez(
                NPZ,
                series_cf=s_cf.values,
                series_h=s_h.values,
                volumes=s_cf.volumes,
                cells=run["cells"],
            )
            print(f"MAT19 (W) wrote {NPZ}", flush=True)
        comm.Barrier()
    else:
        if not os.path.exists(NPZ):
            pytest.fail(f"(W) {NPZ} missing: run W1 (-n 2) first; never recompute")
        z = np.load(NPZ)
        cells_eq = bool(np.array_equal(z["cells"], run["cells"]))
        w_cf = float(np.max(np.abs(s_cf.values[1:] / z["series_cf"][1:] - 1.0)))
        w_h = float(np.max(np.abs(s_h.values[1:] / z["series_h"][1:] - 1.0)))
        w_vol = float(np.max(np.abs(s_cf.volumes / z["volumes"] - 1.0)))
        _say(
            f"MAT19 (W) -n {comm.size} vs -n 2: cells equal {cells_eq} ({list(run['cells'])} vs "
            f"{list(z['cells'])}); series_cf max rel {w_cf:.3e} (band {WIDTH_BAND:.0e}, predicted "
            f"<= 1e-12); series_h max rel {w_h:.3e} (printed, predicted <= 3e-10); volumes "
            f"{w_vol:.3e}"
        )
    _say(
        f"MAT19 elapsed: EM {run['t_em']:.1f} s, two CN {run['t_cn']:.1f} s, probes "
        f"{run['t_probe']:.2f} s, fixture total {run['t_total']:.1f} s"
    )

    assert m_ok, "(M) a probe ball leaves r < 0.9 R"
    assert run["n_sub"] == N_SUB, f"(A0) submesh cells {run['n_sub']} != {N_SUB}"
    assert np.all(s_cf.values[0] == 0.0) and np.all(s_h.values[0] == 0.0), "(A0) snapshot 0"
    assert worstT <= P_BAND, f"(T) {worstT:.4e} > 1 %"
    assert np.max(missC1) > P_BAND, f"(C1) displaced probe did not miss: {np.max(missC1):.4e}"
    if w_cf is not None:
        assert cells_eq, "(W) cells differ across widths"
        assert w_cf <= WIDTH_BAND, f"(W) series_cf {w_cf:.3e} > 1e-9"
    return dict(
        ref=ref,
        pts=pts,
        P=P,
        relT=relT,
        worstT=worstT,
        relH=relH,
        missC1=missC1,
        ownrel=ownrel,
        c2=c2,
        c2n20=c2n20,
        w_cf=w_cf,
        w_h=w_h,
    )


@complex_only
def test_mat19_probe_series(run):
    _check_probe_series(run)
