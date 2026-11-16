"""`IMP-7` step 3 (§9 item 647, §10 I42): P2 on curved (P2-geometry) cells,
metal core only.

Brief: docs/planning/chunks/IMP-7.md, "§9 item 647 — 2026-10-25".  The
question I40 left: its affine P2 L2 rate is ≈ 2 (LSQ 1.988), not 3.  The
polyhedral-floor reading (chord depth ≈ h²/(8a) on the affine r = a / r = b
surfaces) predicts the rate rises toward 3 on second-order geometry and the
error falls well below the affine P2's; "something else" predicts ≈ 2 again.

Asserted (bands imported from `IMP-7`'s module by name, order-free):
(S3c) balance ≤ S3_BAND every curved rung; (S2c) interior rel L2 ≤ S2_BAND on
the finest curved rung run; (c-β0c) the β = 0 source on curved LADDER[1]
misses S2_BAND.  Printed only, no band: (Geo) volumes against 4/3·π·a³ /
4/3·π·b³ on both meshes, (Rate) curved rel L2 and pairwise rate, the
curved / affine ratio, the l = 2 readout.

The curved mesh never enters `m7._ball` (its cache is keyed by res alone);
this module owns `_ball_curved`, cached by (res, 2).  Real build.
"""

from __future__ import annotations

import inspect
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator
from tests.validation import test_imp7_implant_heating_fem as m7

# Names imported from `IMP-7`'s module by name (checked at execution below).
_M7_NAMES = ("LADDER", "DEG_ROUTE", "S2_BAND", "S3_BAND", "A", "B", "K_CORE",
             "K_GEL", "_dg0", "_source_dg2", "_solve", "_p2", "_theta_metal",
             "_q_shell", "_errors", "_allsum", "_indicator", "_volume",
             "_l2_readout", "SERIES_READ", "_ball")
LADDER = m7.LADDER
DEG_ROUTE = m7.DEG_ROUTE
S2_BAND = m7.S2_BAND
S3_BAND = m7.S3_BAND
A = m7.A
B = m7.B
K_CORE = m7.K_CORE
K_GEL = m7.K_GEL

V_CORE = 4.0 / 3.0 * np.pi * A**3
V_BALL = 4.0 / 3.0 * np.pi * B**3
VOL_QMETA = {"quadrature_degree": 4}  # detJ is degree 3 on P2 tets: exact
QCHECK = {"quadrature_degree": 8}  # adequacy check on DEG_ROUTE[2]'s qdeg 6

# I40 affine P2, metal core, `-n 2` (20261025T123436Z_IMP-7-step2-W1.log:61-63):
# (whole, interior) rel L2 per LADDER rung.  Width differs (-n 4 here):
# compare rates and ratios, never digits.
AFFINE_P2 = {LADDER[0]: (1.5827e-3, 1.6370e-3),
             LADDER[1]: (8.1490e-4, 8.4803e-4),
             LADDER[2]: (4.0631e-4, 4.2200e-4)}
AFFINE_P2_PAIRWISE = (1.937, 2.040)  # W1:66-67
AFFINE_P2_READ = (1.2224e-2, 0.362259, 0.343898)  # W1:65 (finest affine rung)
AFFINE_BETA0 = 5.9510e-2  # W1:70 (affine LADDER[1])
# Curved P2 pairwise L2 rate bar — the 2026-10-26 weekly's (§10 I44; §9 item
# 657): 0.85 × the order (3), the family's convention (`S1_RATE` 1.7 for P1,
# test_imp7_implant_heating_fem.py:58).  Records: I42 3.036 / 3.105, I43
# 3.234 / 3.129; the affine floor 1.937–2.043 (I40) misses it.
RATE_BAR_P2_CURVED = 2.6

_CURVED: dict = {}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _check_names():
    missing = [n for n in _M7_NAMES if not hasattr(m7, n)]
    assert not missing, f"m7 lacks {missing}"
    assert not any(n.startswith("test_") for n in _M7_NAMES)


def _ball_curved(res):
    key = (res, 2)
    if key not in _CURVED:
        t0 = time.perf_counter()
        out = MeshGenerator.two_layer_ball_domain(
            inner_radius=A, outer_radius=B, resolution=res, geometry_order=2
        )
        _CURVED[key] = (*out, time.perf_counter() - t0)
    return _CURVED[key]


def _cmap_degree(msh):
    cm = msh.geometry.cmap
    d = getattr(cm, "degree")
    return int(d() if callable(d) else d)


def _volumes(msh, ct):
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=VOL_QMETA)
    vball = m7._allsum(comm, fem.Constant(msh, default_scalar_type(1.0)) * dx)
    vcore = m7._allsum(comm, m7._dg0(msh, ct, 1.0, 0.0) * dx)
    return vcore / V_CORE - 1.0, vball / V_BALL - 1.0, vball


def _indicator_count(msh):
    chi = m7._indicator(msh, lambda m: np.sum(m**2, axis=1) < 0.81 * B**2)
    nloc = msh.topology.index_map(3).size_local
    return int(round(msh.comm.allreduce(float(np.sum(chi.x.array[:nloc])), op=MPI.SUM)))


def _curved_metal_rung(res, q_fn=None):
    """`m7._metal_rung`'s body (:335–351) with the curved mesh explicit and
    rt = DEG_ROUTE[2] (DG3 source, P3 comparand, quadrature 6)."""
    q_fn = m7._q_shell if q_fn is None else q_fn
    rt = DEG_ROUTE[2]
    msh, ct, ft, t_build = _ball_curved(res)
    k_f = m7._dg0(msh, ct, K_CORE, K_GEL)
    q_f = m7._source_dg2(msh, ct.find(2), q_fn, dg=rt["q_dg"])
    sol, t_solve = m7._solve(msh, ft, k_f, q_f, degree=2)
    theta = m7._p2(msh, m7._theta_metal, p=rt["comparand"])
    chi_sh = m7._dg0(msh, ct, 0.0, 1.0)
    e_all, e_in, e_sh = m7._errors(msh, chi_sh, sol.temperature, theta, B,
                                   qmeta=rt["qmeta"])
    e_all8, e_in8, _ = m7._errors(msh, chi_sh, sol.temperature, theta, B,
                                  qmeta=QCHECK)
    nglob = msh.topology.index_map(3).size_global
    _, _, vol = _volumes(msh, ct)
    h_eff = (vol / nglob) ** (1.0 / 3.0)
    bal = abs(sol.convected_power - sol.source_power) / sol.source_power
    return dict(msh=msh, ct=ct, sol=sol, t_build=t_build, t_solve=t_solve,
                nglob=nglob, h_eff=h_eff, e_all=e_all, e_in=e_in,
                e_all8=e_all8, e_in8=e_in8, bal=bal)


def _geo_print(res):
    mc, ctc, _, _ = _ball_curved(res)
    ma, cta, _, _ = m7._ball(res)  # affine (flat) mesh into the flat cache only
    dc_core, dc_ball, _ = _volumes(mc, ctc)
    da_core, da_ball, _ = _volumes(ma, cta)
    _say(f"IMP7-3 (Geo) res={res} cells curved={mc.topology.index_map(3).size_global} "
         f"affine={ma.topology.index_map(3).size_global}: "
         f"Vcore/(4/3πa³)-1 curved={dc_core:.4e} affine={da_core:.4e}; "
         f"Vball/(4/3πb³)-1 curved={dc_ball:.4e} affine={da_ball:.4e}; "
         f"interior-indicator cells curved={_indicator_count(mc)} "
         f"affine={_indicator_count(ma)} (printed)")


def _rung_report(res, r):
    aw, ai = AFFINE_P2[res]
    _say(f"IMP7-3 [metal curved P2] res={res} cells={r['nglob']} h_eff={r['h_eff']:.5e} "
         f"build={r['t_build']:.1f}s solve={r['t_solve']:.1f}s "
         f"whole_rel_L2={r['e_all']:.4e} interior_rel_L2={r['e_in']:.4e} "
         f"(qdeg8 check: whole {r['e_all8']:.4e} interior {r['e_in8']:.4e}) "
         f"balance={r['bal']:.3e}")
    _say(f"IMP7-3 [metal] res={res} curved/affine ratio (I40 -n 2): "
         f"whole {r['e_all'] / aw:.4f} interior {r['e_in'] / ai:.4f}")
    assert r["bal"] <= S3_BAND, f"(S3c) balance {r['bal']:.3e} > {S3_BAND} at res={res}"


def _pairwise(r0, r1, res0, res1, i):
    rw = np.log(r0["e_all"] / r1["e_all"]) / np.log(r0["h_eff"] / r1["h_eff"])
    ri = np.log(r0["e_in"] / r1["e_in"]) / np.log(r0["h_eff"] / r1["h_eff"])
    _say(f"IMP7-3 (Rate) [metal curved P2] pairwise {res0}->{res1}: whole {rw:.3f} "
         f"interior {ri:.3f} (affine P2 I40 whole {AFFINE_P2_PAIRWISE[i]:.3f}) "
         f"— band whole >= RATE_BAR_P2_CURVED {RATE_BAR_P2_CURVED}")
    assert rw >= RATE_BAR_P2_CURVED, (
        f"(Rate-c) [metal] curved P2 pairwise {res0}->{res1} {rw:.3f} < {RATE_BAR_P2_CURVED}")
    return rw


def _readout(r, res):
    _say(f"IMP7-3 l=2 readout on curved res={res}; affine P2 (I40, finest affine rung): "
         f"R2_h={AFFINE_P2_READ[0]:.4e} pole={AFFINE_P2_READ[1]:.6f} "
         f"equator={AFFINE_P2_READ[2]:.6f} K")
    m7._l2_readout(r["sol"].temperature, f"metal curved P2 res={res}",
                   m7.SERIES_READ["metal"])


# ---- tests --------------------------------------------------------------------


def test_imp7_3_w0_coding():
    """W0: names, the keyword's contract, the coarse curved ball."""
    _check_names()
    sig = inspect.signature(MeshGenerator.two_layer_ball_domain)
    assert sig.parameters["geometry_order"].default == 1
    assert list(sig.parameters)[:5] == [
        "inner_radius", "outer_radius", "resolution", "comm", "rank"]
    with pytest.raises(ValueError):
        MeshGenerator.two_layer_ball_domain(
            inner_radius=A, outer_radius=B, resolution=4e-3, geometry_order=3)
    t0 = time.perf_counter()
    mc, ctc, ftc = MeshGenerator.two_layer_ball_domain(
        inner_radius=A, outer_radius=B, resolution=4e-3, geometry_order=2)
    t_c = time.perf_counter() - t0
    ma, cta, fta = MeshGenerator.two_layer_ball_domain(
        inner_radius=A, outer_radius=B, resolution=4e-3)
    dc, da = _cmap_degree(mc), _cmap_degree(ma)
    nc = mc.topology.index_map(3).size_global
    na = ma.topology.index_map(3).size_global
    vc = _volumes(mc, ctc)
    va = _volumes(ma, cta)
    _say(f"IMP7-3 W0 res=4e-3: cmap degree curved={dc} affine={da}; cells "
         f"curved={nc} affine={na}; build curved={t_c:.1f}s; "
         f"Vcore-1 curved={vc[0]:.4e} affine={va[0]:.4e}; "
         f"Vball-1 curved={vc[1]:.4e} affine={va[1]:.4e}")
    assert dc == 2, f"curved cmap degree {dc}"
    assert da == 1, f"default cmap degree {da}"
    assert nc == na, f"setOrder(2) changed the cell count: {nc} vs {na}"
    # The curved core volume must beat the affine chord deficit (O(h²) -> O(h⁴)).
    assert abs(vc[0]) < abs(va[0]), f"curved Vcore {vc[0]:.3e} vs affine {va[0]:.3e}"


def test_imp7_3c_metal_record_control():
    """Record control (no mesh): the affine floor I40 measured misses the curved bar."""
    _say(f"IMP7-3c record control [metal]: affine P2 pairwise (I40) {AFFINE_P2_PAIRWISE} "
         f"max {max(AFFINE_P2_PAIRWISE):.3f} vs RATE_BAR_P2_CURVED {RATE_BAR_P2_CURVED}")
    assert max(AFFINE_P2_PAIRWISE) < RATE_BAR_P2_CURVED, (
        "record control: the affine floor I40 measured misses the curved bar")


def test_imp7_3_metal_curved_w1():
    """W1: curved LADDER[0], LADDER[1] at P2, plus (c-β0c) on LADDER[1]."""
    _check_names()
    rows = {}
    for res in LADDER[:2]:
        _geo_print(res)
        rows[res] = _curved_metal_rung(res)
        _rung_report(res, rows[res])
    _pairwise(rows[LADDER[0]], rows[LADDER[1]], LADDER[0], LADDER[1], 0)
    _readout(rows[LADDER[1]], LADDER[1])
    fin = rows[LADDER[1]]
    _say(f"IMP7-3 (S2c) curved interior rel L2 at res={LADDER[1]} = {fin['e_in']:.4e} "
         f"(band <= S2_BAND {S2_BAND})")
    assert fin["e_in"] <= S2_BAND, f"(S2c) interior {fin['e_in']:.4e} > {S2_BAND}"
    c = _curved_metal_rung(LADDER[1], q_fn=lambda x: m7._q_shell(x, beta=0.0))
    _say(f"IMP7-3 (c-beta0c) curved res={LADDER[1]}: interior miss {c['e_in']:.4e} "
         f"(affine P2 I40 {AFFINE_BETA0:.4e}; asserted > {S2_BAND}); "
         f"whole {c['e_all']:.4e}; balance {c['bal']:.3e}; solve={c['t_solve']:.1f}s")
    assert c["bal"] <= S3_BAND, f"(S3c) c-beta0c balance {c['bal']:.3e}"
    assert c["e_in"] > S2_BAND, f"(c-beta0c) control did not fire: {c['e_in']:.4e}"


def test_imp7_3_metal_curved_w2_finest():
    """W2 (optional): curved LADDER[1] and LADDER[2]; (S2c) on the finest."""
    _check_names()
    rows = {}
    for res in LADDER[1:]:
        if res == LADDER[2]:
            _geo_print(res)
        rows[res] = _curved_metal_rung(res)
        _rung_report(res, rows[res])
    _pairwise(rows[LADDER[1]], rows[LADDER[2]], LADDER[1], LADDER[2], 1)
    _readout(rows[LADDER[2]], LADDER[2])
    fin = rows[LADDER[2]]
    _say(f"IMP7-3 (S2c) curved finest interior rel L2 = {fin['e_in']:.4e} "
         f"(band <= S2_BAND {S2_BAND})")
    assert fin["e_in"] <= S2_BAND, f"(S2c) interior {fin['e_in']:.4e} > {S2_BAND}"
