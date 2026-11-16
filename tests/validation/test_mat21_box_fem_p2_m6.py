"""`MAT-21` step 5 — step 1's steady anchors at Lagrange degree 2 against a P3
comparand on the m = 6 rung (§9 item 455, §10 chain P30, daily 2026-10-14 —
provisional; brief `docs/planning/chunks/MAT-21.md`).

A sibling of step 3's ``test_mat21_box_fem_p2.py``: imports its degree-free
helpers, edits no landed module, sets no band and moves none.  Steady only —
the transient and duty legs at m = 6 P2 are priced in the brief and excluded.
m = 2–4 are NOT re-solved: their whole-box errors and h_eff are step 3's logged
values (``STEP3_WHOLE`` / ``STEP3_HEFF``) and feed only the printed rate.

Asserted: (P0) step-0 constants = FEM's; (Pdeg) global P2 dofs on m = 6 =
380 725; (S1-int) m = 6 interior ≤ ``INTERIOR_BAND``; (S2) Gauss ≤
``S2_GAUSS_BAND``.  Controls (separations only): (c1) 2h comparand interior >
``INTERIOR_BAND``; (c-deg) P1 interior > P2 interior on the same rung.
Printed only: whole-box error, the four-rung P2 rate, pairwise rates, the
comparand-floor ratios, theta_h(0), RSS, clocks.
"""

from __future__ import annotations

import time

import numpy as np
import pytest
import ufl
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.thermal import solve_steady_conduction

from tests.complex_mode import complex_only
from tests.validation.test_imp6_oblate_fem import _rss_gib
from tests.validation.test_mat10_sar_to_heat_sphere import ROBIN_TAG, _allsum, _say
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, T_INF
from tests.validation.test_mat15_finite_cylinder_fem import (
    RATE_BAR,
    S2_GAUSS_BAND,
    _max_imag,
)
from tests.validation.test_mat21_box_closed_form import H, K, Q
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0
from tests.validation.test_mat21_box_fem import V_BOX, _interior_chi, _steady_ex
from tests.validation.test_mat21_box_fem_p2 import (
    CMP_DEGREE,
    COMPARAND_FLOOR,
    DEGREE,
    QDEG_P,
    _dofs_pred,
    _exact_p,
    _rel_p,
    _rung_p,
)

M6 = 6
# Step 3 W1 (20261014T095542Z_MAT-21-step3-W1-p2-steady.log :94, :101, :105),
# copied to every printed digit; they feed the printed rate only, never an assertion.
STEP3_WHOLE = {2: 1.884027e-4, 3: 5.851175e-5, 4: 2.489307e-5}
STEP3_HEFF = {2: 1.350103e-2, 3: 9.000686e-3, 4: 6.750514e-3}


class _Clock5:
    """Per-phase wall clock, printed on rank 0 with flush, prefix MAT21-5."""

    def __init__(self, tag):
        self.tag = tag
        self.t = time.perf_counter()

    def lap(self, phase):
        now = time.perf_counter()
        _say(f"MAT21-5 [{self.tag}] clock {phase}: {now - self.t:.1f} s")
        self.t = now


@pytest.fixture(scope="module")
def steady_m6():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    clk = _Clock5("steady")
    out = {}
    _say(f"\nMAT21-5 DEGREE = {DEGREE}, comparand P{CMP_DEGREE}, quadrature {QDEG_P}; "
         f"rung m = {M6}; comparand floor = {COMPARAND_FLOOR:.4e} (printed)")

    # 1. solve, then RSS while the factor is alive.
    msh, ft, ns, n, v_h, st = _rung_p(M6, clk, DEGREE)
    out["rss"] = _rss_gib(comm)
    _say(f"MAT21-5 RSS after P2 solve (sum over ranks) = {out['rss']:.3f} GiB (predicted 6-30)")
    out["ndof"] = int(st.temperature.function_space.dofmap.index_map.size_global)
    out["imag"] = _max_imag(st.temperature)

    # 2. Gauss.
    ds = ufl.Measure("ds", domain=msh, subdomain_data=ft)(ROBIN_TAG)
    th_h = st.temperature - T_INF
    conv = float(np.real(_allsum(comm, H_GEL * th_h * ds)))
    out["gauss"] = abs(conv - Q * v_h) / (Q * v_h)
    clk.lap("Gauss")

    # 3. comparand, whole / interior, centre.
    ex, npts = _exact_p(msh, _steady_ex())
    clk.lap(f"steady series interpolation m={M6} ({npts} P3 points)")
    chi = _interior_chi(msh)
    out["whole"] = _rel_p(msh, th_h, ex)
    out["interior"] = _rel_p(msh, th_h, ex, chi)
    h_eff = (v_h / n) ** (1.0 / 3.0)
    out["h_eff"] = h_eff
    try:
        vals, ok = evaluate_vector_field_parallel(st.temperature, np.zeros((1, 3)))
        th_c = float(np.real(np.ravel(vals[0])[0])) - T_INF if bool(ok[0]) else float("nan")
    except Exception as exc:  # printed quantity only
        _say(f"MAT21-5 centre evaluation failed: {exc!r}")
        th_c = float("nan")
    th_c_ex = float(np.real(_steady_ex()(np.zeros(1), np.zeros(1), np.zeros(1))[0]))
    _say(f"MAT21-5 rung m = {M6} {tuple(ns)}: cells {n} (6*prod N = {6 * int(np.prod(ns))}); "
         f"P2 dofs {out['ndof']} (predicted {_dofs_pred(M6)}); P3 points {npts} "
         f"(predicted {_dofs_pred(M6, CMP_DEGREE)}, 1260775); V_h = {v_h:.15e} vs 8ABC "
         f"{V_BOX:.15e} ({v_h / V_BOX - 1.0:+.3e}); h_eff = {h_eff:.6e} m (predicted 4.500343e-3)")
    _say(f"MAT21-5 (S1-int) interior rel L2 = {out['interior']:.6e} (bound {INTERIOR_BAND:.0e}, "
         f"predicted 5.0e-6, range 4e-6 to 1.2e-5; step 3 P2 m = 4 read 1.669359e-5, different mask); "
         f"/floor = {out['interior'] / COMPARAND_FLOOR:.4f}")
    _say(f"MAT21-5 whole rel L2 = {out['whole']:.6e} (printed; predicted 7.5e-6, range 6.5e-6 "
         f"to 1.2e-5); /floor = {out['whole'] / COMPARAND_FLOOR:.4f}")
    _say(f"MAT21-5 (S2) Gauss = {out['gauss']:.3e} (bound {S2_GAUSS_BAND:.0e}, predicted <= 5e-12)")
    _say(f"MAT21-5 theta_h(0) = {th_c:.6f} K vs 40-term series {th_c_ex:.6f} K "
         f"(rel {th_c / th_c_ex - 1.0:+.3e}) and 120-term record 12.003957 K "
         f"(predicted 12.003957-12.003961; P2 m = 4 read 12.003960)")
    clk.lap("readings")

    # 4. (c1) comparand at 2h.
    ex_c1, _ = _exact_p(msh, _steady_ex(h=2.0 * H))
    out["c1"] = _rel_p(msh, th_h, ex_c1, chi)
    _say(f"MAT21-5 control (c1) comparand at 2h: interior miss = {out['c1']:.5f} "
         f"(predicted 0.5973, range 0.5960-0.5990; step 1 P1 m = 6 read 0.59734 same mask; "
         f"P2 m = 4 read 0.60389, different mask; asserted > {INTERIOR_BAND:.0e})")
    clk.lap("control (c1)")

    # 5. (c-deg) degree 1 on the same mesh and ft, same P3 comparand and chi.
    st1 = solve_steady_conduction(msh, ft, k=K_GEL, q=Q, h=H_GEL, t_inf=T_INF,
                                  robin_tag=ROBIN_TAG, degree=1)
    clk.lap(f"(c-deg) steady solve m={M6} degree=1")
    th1 = st1.temperature - T_INF
    out["e1_int"] = _rel_p(msh, th1, ex, chi)
    out["e1_whole"] = _rel_p(msh, th1, ex)
    ndof1 = int(st1.temperature.function_space.dofmap.index_map.size_global)
    _say(f"MAT21-5 control (c-deg) m = {M6}: P1 dofs {ndof1} (predicted {_dofs_pred(M6, 1)}); "
         f"interior P1 {out['e1_int']:.6e} (predicted 1.41e-3 to 1.43e-3; step 1 1.417055e-3) vs "
         f"P2 {out['interior']:.6e}, ratio e1/e2 = {out['e1_int'] / out['interior']:.2f} "
         f"(predicted 280, range 120-360; asserted > 1); whole P1 {out['e1_whole']:.6e} "
         f"(predicted 1.540e-3 to 1.555e-3; step 1 1.547121e-3)")
    clk.lap("control (c-deg)")

    # 6. printed rate.
    ms = (2, 3, 4)
    hs = np.array([STEP3_HEFF[m] for m in ms] + [h_eff])
    es = np.array([STEP3_WHOLE[m] for m in ms] + [out["whole"]])
    out["slope"] = float(np.polyfit(np.log(hs), np.log(es), 1)[0])
    out["pairwise"] = [float(np.log(es[i] / es[i + 1]) / np.log(hs[i] / hs[i + 1]))
                       for i in range(len(hs) - 1)]
    _say(f"MAT21-5 P2 rate m in (2,3,4,6) (m = 2-4 from step 3 W1, printed only): LSQ slope = "
         f"{out['slope']:.4f} (≈ 3; RATE_BAR {RATE_BAR} named, not asserted; predicted 2.85-3.0); "
         f"pairwise " + ", ".join(f"{p:.4f}" for p in out["pairwise"])
         + " (4 -> 6 predicted 2.9, range 1.8-3.3; < 2.5 reads as comparand floor)")
    _say(f"MAT21-5 max|Im theta| = {out['imag']:.3e} (predicted 0); "
         f"elapsed {time.perf_counter() - t0:.1f} s")
    return out


def test_precondition_constants():
    assert K == K_GEL
    assert H == H_GEL
    assert RHO_C_0 == RHO_C


@complex_only
def test_precondition_degree_dofs_m6(steady_m6):
    assert steady_m6["ndof"] == _dofs_pred(M6) == 380725, (
        f"(Pdeg) {steady_m6['ndof']} != {_dofs_pred(M6)}")


@complex_only
def test_s1_p2_m6_interior(steady_m6):
    assert steady_m6["interior"] <= INTERIOR_BAND, f"(S1-int) {steady_m6['interior']:.6e}"


@complex_only
def test_s2_p2_m6_gauss_balance(steady_m6):
    assert steady_m6["gauss"] <= S2_GAUSS_BAND, f"(S2) {steady_m6['gauss']:.3e}"


@complex_only
def test_control_c1_p2_m6_doubled_h_misses(steady_m6):
    assert steady_m6["c1"] > INTERIOR_BAND, f"(c1) {steady_m6['c1']:.4e}"


@complex_only
def test_control_cdeg_p1_m6_worse(steady_m6):
    assert steady_m6["e1_int"] > steady_m6["interior"], (
        f"(c-deg) e1 {steady_m6['e1_int']:.6e} <= e2 {steady_m6['interior']:.6e}")
