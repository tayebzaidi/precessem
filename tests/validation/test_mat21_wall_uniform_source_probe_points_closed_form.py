"""`MAT-21` step 42 — the uniform-source walled box at entry 59's six probe points, closed form
(§9 item 876, §10 P75; brief `docs/planning/chunks/MAT-21.md` "§9 item 876").

Numpy.  Sets no band and moves none.  Three forms at P37's six points ``PROBE_PTS``:

- **W_src(t)** — P48a″'s 1D-reduced walled sourced rise (P64's ``_step``), at the point's
  height z ∈ {0, C}; P1–P4 (side-wall points the 1D reduction cannot see) carry the
  centre value, labelled "1D: z = 0".  History = P64's ``_uniform_w((0.0,))`` (row 0 = 0).
- **W_drift(t)** — P67's capacitive-wall drift at t_w = 0 (``_drift_cap_z``), unscaled
  (it carries no watt); history = P70's ``_cap_hists()["0"]``.
- **U3(t)** — step 0's six-face unwalled uniform-source series
  (``robin_box_transient_temperature``, ``SERIES_TERMS``); history row 0 = 0.

Per watt: P_dep = Q·8ABC = 49.14 W.  Asserted (collected, then asserted — item 668): the
preconditions; (Sup) q = 0 / g = 0 limits, homogeneity θ(2Q, 2g) = 2θ(Q, g), per-watt
invariance of U3, at ``B1_BAND``; (Lim) U3(1000 T_C) against the steady series at
``B1_BAND`` and P37's six steady values within 5e-7 K; (Bq) drift N vs 2N at ``B1_BAND``;
MAT-20's (Q-d) with the R-swap control and (Win-d) on U3/P, W_src/P and W_src/P + W_drift;
the two reproductions; controls (c-wall) and (c-drift).  Where a limit's reference is
exactly 0 (the odd drift at z = 0) the error is the absolute difference in K (``_err``).

Not claimed: entry 59's figures (P76), a non-uniform SAR, a 3D walled box, a dose band, a
compliance figure, a cited PMMA ρc.
"""

from __future__ import annotations

import numpy as np

from fem_em_solver.utils.analytical import (
    robin_box_steady_temperature,
    robin_box_transient_temperature,
    two_layer_slab_eigenvalues,
    two_layer_slab_transient_temperature,
)
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL
from tests.validation.test_mat15_coil_driven_heating_phantom import T_INF as T_INF_59
from tests.validation.test_mat21_box_closed_form import (
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_closed_form import SERIES_TERMS, T_C, T_INF, H, K, Q
from tests.validation.test_mat21_box_dose_closed_form import AVG_BAR, DOSE_BAR, T_GRID, W_B2
from tests.validation.test_mat21_box_dt_duty125_dose_closed_form import _q_check_g, _win_check_g
from tests.validation.test_mat21_box_fem import T_STAR as T_STAR_BOX
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS
from tests.validation.test_mat21_initial_gradient_closed_form import G_GRAD, N_TERMS
from tests.validation.test_mat21_wall_capacity_drift_closed_form import _drift_cap, _drift_cap_z
from tests.validation.test_mat21_wall_capacity_duty25_dose_closed_form import _cap_hists
from tests.validation.test_mat21_wall_pre_equilibration_duty25_dose_closed_form import (
    Z,
    _step,
    _uniform_w,
)
from tests.validation.test_mat21_wall_transient_closed_form import N_SERIES, T_STAR, WALL

P_DEP = Q * 8.0 * A_BOX * B_BOX * C_BOX
FEM_P_DEP = 49.14  # `20261105T225057Z_MAT-21-P69-W1.log:53`, source_power 4.914000000e+01 W
TIMES = (("t*", T_STAR), ("360 s", float(W_B2)))
# P37's steady record (`20261018T111612Z_MAT-21-step10-W2a-points-n2.log.gz:90–95`), K.
P37_STEADY = (12.003947, 6.236654, 6.236654, 6.150392, 6.150392, 8.733966)
P37_TOL = 5e-7
REPRO_SRC_P0 = (1.327311508, 5e-10)  # `20261028T021107Z_MAT-21-step20c-W1-fem.log.gz:55`
REPRO_DRIFT_P5 = (1.819866e-01, 5e-8)  # `…P67-W1.log:68`
ZC = {0.0: 2, float(C_BOX): 0}  # height -> column of Z
_LAM = {}
_CACHE = {}


def _say(msg):
    print(msg, flush=True)


def _err(a, b):
    """Relative error; the absolute difference (K) where the reference is exactly 0."""
    return abs(a - b) / abs(b) if b != 0.0 else abs(a - b)


def _col(i):
    return ZC[float(PROBE_PTS[i][2])]


def _lam(n=N_SERIES):
    if n not in _LAM:
        _LAM[n] = two_layer_slab_eigenvalues(h=H, n=n, **WALL)
    return _LAM[n]


def _src_q(t, q, n=N_SERIES):
    return np.asarray(
        two_layer_slab_transient_temperature(
            Z, t, q=q, h=H, t_inf=0.0, n_terms=n, eigenvalues=_lam(n), **WALL
        ),
        dtype=float,
    )


def _u3(t, q=Q, n_terms=SERIES_TERMS):
    P = np.asarray(PROBE_PTS, dtype=float)
    return np.asarray(
        robin_box_transient_temperature(
            P[:, 0],
            P[:, 1],
            P[:, 2],
            t,
            half_widths=HALF,
            k=K,
            q=q,
            h=H,
            rho_c=RHO_C_0,
            t_inf=0.0,
            n_terms=n_terms,
        ),
        dtype=float,
    )


def _hists():
    """Per-watt / unscaled (41, 6) histories on T_GRID, read as one batch each."""
    if not _CACHE:
        u = np.zeros((T_GRID.size, len(PROBE_PTS)))
        for n in range(1, T_GRID.size):
            u[n] = _u3(T_GRID[n])
        ws = _uniform_w((0.0,))
        wd = _cap_hists()["0"]
        cols = [_col(i) for i in range(len(PROBE_PTS))]
        _CACHE["U3"] = u / P_DEP
        _CACHE["Wsrc"] = ws[:, cols] / P_DEP
        _CACHE["Wsd"] = ws[:, cols] / P_DEP + wd[:, cols]
        _CACHE["Wsrc_raw"] = ws
    return _CACHE


def _theta(j, t, q=Q, g=G_GRAD):
    return float(_src_q(t, q)[j]) + _drift_cap((0.0, 0.0, Z[j]), t, g=g)


# ------------------------------------------------------------------ tests
def test_preconditions_and_height_map():
    gates = []
    zs = [float(p[2]) for p in PROBE_PTS]
    gates.append((("T_STAR", T_STAR, T_STAR_BOX), abs(T_STAR - T_STAR_BOX) <= 1e-9, 0))
    gates.append((("K", K, K_GEL), K == K_GEL, 0))
    gates.append((("H", H, H_GEL), H == H_GEL, 0))
    gates.append((("T_INF", T_INF, T_INF_59), T_INF == T_INF_59, 0))
    pr = abs(P_DEP / FEM_P_DEP - 1.0)
    gates.append((("P_dep", P_DEP), pr <= 1e-12, pr))
    gates.append((("z map", zs), all(z in (0.0, float(C_BOX)) for z in zs), 0))
    gates.append((("Z", Z.tolist()), Z[0] == C_BOX and Z[2] == 0.0, 0))
    gates.append((("row0",), bool(np.all(_hists()["Wsrc_raw"][0] == 0.0)), 0))
    _say(
        f"MAT21-42 T_STAR {T_STAR:.6f} s vs T_STAR_BOX {T_STAR_BOX:.6f} s; K {K} K_GEL {K_GEL}; "
        f"H {H} H_GEL {H_GEL}; T_INF {T_INF} T_INF_59 {T_INF_59}; P_dep {P_DEP:.12e} W vs FEM "
        f"{FEM_P_DEP} W, |ratio-1| {pr:.3e} (bound 1e-12); heights {zs}; Z {Z.tolist()}; "
        f"T_GRID {T_GRID.size} rows to {T_GRID[-1]:.4f} s"
    )
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_sup_lim_bq_and_reproductions():
    gates = []
    w = {"src-eq": 0.0, "q0": 0.0, "g0": 0.0, "hom": 0.0, "pw": 0.0, "Bq": 0.0}
    for tl, t in TIMES:
        s = _step(t)
        sq = _src_q(t, Q)
        u1, u2 = _u3(t) / P_DEP, _u3(t, q=2 * Q) / (2 * P_DEP)
        w["pw"] = max(w["pw"], float(np.max(np.abs(u2 - u1) / np.abs(u1))))
        for zl, j in (("z=C", 0), ("z=0", 2)):
            p = (0.0, 0.0, Z[j])
            d = _drift_cap(p, t)
            w["src-eq"] = max(w["src-eq"], _err(sq[j], s[j]))
            w["q0"] = max(w["q0"], _err(float(_src_q(t, 0.0)[j]) + d, d))
            w["g0"] = max(w["g0"], _err(float(sq[j]) + _drift_cap(p, t, g=0.0), float(s[j])))
            w["hom"] = max(w["hom"], _err(_theta(j, t, 2 * Q, 2 * G_GRAD), 2 * _theta(j, t)))
            dz, dz2 = float(_drift_cap_z(Z, t)[j]), float(_drift_cap_z(Z, t, n=2 * N_TERMS)[j])
            w["Bq"] = max(w["Bq"], _err(dz, dz2))
            _say(f"MAT21-42 (Bq) {tl} {zl}: drift N {dz:+.12e} K, 2N {dz2:+.12e} K")
    for k_, v in w.items():
        _say(
            f"MAT21-42 ({k_}) max over both heights, t* and 360 s: {v:.3e} (band {B1_BAND:g}, "
            f"ASSERTED)"
        )
        gates.append(((k_,), v <= B1_BAND, v))
    # (Lim)
    ss = np.asarray(
        robin_box_steady_temperature(
            *np.asarray(PROBE_PTS, dtype=float).T,
            half_widths=HALF,
            k=K,
            q=Q,
            h=H,
            t_inf=0.0,
            n_terms=SERIES_TERMS,
        ),
        dtype=float,
    )
    ul = _u3(1000.0 * T_C)
    lim = float(np.max(np.abs(ul - ss) / np.abs(ss)))
    _say(
        f"MAT21-42 (Lim) U3(1000 T_C) vs steady series max rel {lim:.3e} (band {B1_BAND:g}, "
        f"ASSERTED; predicted 0)"
    )
    gates.append((("Lim",), lim <= B1_BAND, lim))
    for i, lab in enumerate(PROBE_LABELS):
        dv = abs(ul[i] - P37_STEADY[i])
        _say(
            f"MAT21-42 (Lim) {lab}: U3(1000 T_C) {ul[i]:.9f} K, steady {ss[i]:.9f} K, P37 "
            f"{P37_STEADY[i]:.6f} K, |diff| {dv:.3e} (bound {P37_TOL:g})"
        )
        gates.append(((lab, "P37"), dv <= P37_TOL, dv))
    # reproductions
    ws0 = float(_step(T_STAR)[2])
    wd5 = float(_drift_cap_z(Z, T_STAR)[0])
    _say(
        f"MAT21-42 (repro) W_src P0 t* {ws0:.12f} K vs {REPRO_SRC_P0[0]} |diff| "
        f"{abs(ws0 - REPRO_SRC_P0[0]):.3e} (bound {REPRO_SRC_P0[1]:g}); W_drift P5 t* "
        f"{wd5:+.9e} K vs {REPRO_DRIFT_P5[0]:+.6e} |diff| {abs(wd5 - REPRO_DRIFT_P5[0]):.3e} "
        f"(bound {REPRO_DRIFT_P5[1]:g})"
    )
    gates.append((("repro-src",), abs(ws0 - REPRO_SRC_P0[0]) <= REPRO_SRC_P0[1], ws0))
    gates.append((("repro-drift",), abs(wd5 - REPRO_DRIFT_P5[0]) <= REPRO_DRIFT_P5[1], wd5))
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_q_d_win_d_on_per_watt_histories():
    h = _hists()
    gates = []
    worst = {"q": 0.0, "w": 0.0, "rc": np.inf}
    for key in ("U3", "Wsrc", "Wsd"):
        for i, lab in enumerate(PROBE_LABELS):
            th = h[key][:, i]
            _d, _quad, rel, _swap, miss = _q_check_g(T_INF + th, T_GRID)
            mx, t_mx, lo, hi, _ref, wrel = _win_check_g(th, T_GRID)
            gates.append(((key, lab, "Q"), rel <= DOSE_BAR, rel))
            gates.append(((key, lab, "swap"), miss > DOSE_BAR, miss))
            gates.append(
                (
                    (key, lab, "win-loc"),
                    lo <= mx <= hi and T_GRID[0] + W_B2 <= t_mx <= T_GRID[-1],
                    t_mx,
                )
            )
            gates.append(((key, lab, "Win"), wrel <= AVG_BAR, wrel))
            worst["q"] = max(worst["q"], rel)
            worst["w"] = max(worst["w"], wrel)
            worst["rc"] = min(worst["rc"], miss)
    _say(
        f"MAT21-42 (Q-d) worst rel {worst['q']:.3e} (bar {DOSE_BAR:g}, ASSERTED); (Win-d) worst "
        f"rel {worst['w']:.3e} (bar {AVG_BAR:g}, ASSERTED); R swap least miss {worst['rc']:.4e} "
        f"(asserted > {DOSE_BAR:g}); 3 histories x 6 points"
    )
    bad = [g for g in gates if not g[1]]
    assert not bad, bad


def test_controls_c_wall_c_drift():
    u0 = float(_u3(T_STAR)[0])
    ws0 = float(_step(T_STAR)[2])
    cw = abs(ws0 / u0 - 1.0)
    ws5 = float(_step(T_STAR)[0]) / P_DEP
    wd5 = float(_drift_cap_z(Z, T_STAR)[0])
    cd = abs((ws5 + wd5) / ws5 - 1.0)
    d0 = float(_drift_cap_z(Z, T_STAR)[2])
    _say(
        f"MAT21-42 (c-wall) P0 t*: W_src {ws0:.9e} K vs U3 {u0:.9e} K, |W_src/U3-1| {cw:.4e} "
        f"(ASSERTED > {B1_BAND:g}; no log predicts it)"
    )
    _say(
        f"MAT21-42 (c-drift) P5 t*: W_src/P {ws5:.9e} K/W, + W_drift {ws5 + wd5:.9e}, rel miss "
        f"{cd:.4e} (ASSERTED > {B1_BAND:g}); drift at P0 {d0!r} K (0 by odd symmetry)"
    )
    assert cw > B1_BAND and cd > B1_BAND, (cw, cd)


def test_record_printed():
    """The record (no band); asserts only that it was formed."""
    h = _hists()
    _say(f"MAT21-42 P_dep {P_DEP:.6f} W (Q {Q:.6e} W/m^3 x 8ABC)")
    for tl, t in TIMES:
        u = _u3(t) / P_DEP
        s = _step(t)
        dr = _drift_cap_z(Z, t)
        _say(
            f"MAT21-42 TABLE {tl} point | U3/P K/W | W_src/P K/W | W_drift K | W_src/P+W_drift "
            f"| W_src/U3"
        )
        for i, lab in enumerate(PROBE_LABELS):
            j = _col(i)
            tag = " (1D: z = 0)" if i in (1, 2, 3, 4) else ""
            ws = float(s[j]) / P_DEP
            _say(
                f"MAT21-42 TABLE {tl} {lab}{tag} | {u[i]:.7g} | {ws:.7g} | {float(dr[j]):+.7g} | "
                f"{ws + float(dr[j]):.7g} | {ws / u[i]:.7g}"
            )
        if tl == "t*":
            hot = int(np.argmax(u))
            _say(f"MAT21-42 hottest by U3 at t*: {PROBE_LABELS[hot]} (predicted P0 centre)")
            assert u[hot] > 0.0
    _say(
        "MAT21-42 DOSE point | cem43/W U3 | W_src | W_src+drift (s) | max 360 s mean/W U3 | "
        "W_src | W_src+drift (K)"
    )
    for i, lab in enumerate(PROBE_LABELS):
        tag = " (1D: z = 0)" if i in (1, 2, 3, 4) else ""
        ds = [_q_check_g(T_INF + h[k][:, i], T_GRID)[0] for k in ("U3", "Wsrc", "Wsd")]
        ms = [_win_check_g(h[k][:, i], T_GRID)[0] for k in ("U3", "Wsrc", "Wsd")]
        _say(
            f"MAT21-42 DOSE {lab}{tag} | {ds[0]:.7g} | {ds[1]:.7g} | {ds[2]:.7g} | {ms[0]:.7g} | "
            f"{ms[1]:.7g} | {ms[2]:.7g}"
        )
    # (Bq) printed only
    bq_s = max(
        _err(float(a), float(b))
        for t in (T_STAR, float(W_B2))
        for a, b in zip(_src_q(t, Q), _src_q(t, Q, n=2 * N_SERIES))
    )
    _say(
        f"MAT21-42 (Bq-src) N_SERIES {N_SERIES} vs 2N max rel {bq_s:.3e} (predicted <= 1e-12, "
        f"printed only)"
    )
    u40, u80 = _u3(T_STAR), _u3(T_STAR, n_terms=(80, 80, 80))
    rel3 = np.abs(u40 / u80 - 1.0)
    _say(
        "MAT21-42 (Bq-3D) (40,40,40) vs (80,80,80) at t*, rel per point: "
        + ", ".join(f"{lab.split()[0]} {r:.3e}" for lab, r in zip(PROBE_LABELS, rel3))
        + " (printed; steady tail 1.4e-07 to 3.0e-06)"
    )
