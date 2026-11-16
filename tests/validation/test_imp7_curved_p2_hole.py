"""`IMP-7` step 3b (§9 item 651, §10 I43): P2 on curved (P2-geometry) cells,
adiabatic hole only.

Brief: docs/planning/chunks/IMP-7.md, "§9 item 651 — 2026-10-26".  I42
(step 3) confirmed the polyhedral-floor reading on the metal core (curved P2
pairwise rate ≈ 3 against affine ≈ 2).  The hole's inner surface carries the
natural (zero-flux) condition on r = a, not a conductivity jump; this module
measures whether the floor reading holds on that second boundary type.

Asserted (bands imported from `IMP-7`'s module by name, order-free):
(G-hole-c) census / vertex radii / Robin area on the curved tag-2 sub-mesh,
and its coordinate element is degree 2; (S3h) balance ≤ S3_BAND every curved
rung and the control; (S2h) interior rel L2 ≤ S2_BAND on the finest curved
rung run; (c-β0h) the β = 0 source on curved LADDER[1] misses S2_BAND; no NaN
in the hole comparand.  Printed only, no band: shell volumes (curved and
affine), r = a / r = b facet counts and areas on both meshes, curved rel L2
and pairwise rate, curved / affine ratio against I40, the l = 2 readout.

The curved mesh never enters `m7._ball` (its cache is keyed by res alone);
`_ball_curved` (step 3's module) caches it by (res, 2).  Real build.
"""

from __future__ import annotations

import time

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.io.mesh import MeshGenerator
from tests.validation import test_imp7_implant_heating_fem as m7
from tests.validation.test_imp7_curved_p2 import (RATE_BAR_P2_CURVED, _ball_curved,
                                                  _cmap_degree)

# Names read from `IMP-7`'s module by name (checked at execution below),
# including everything `m7._hole_rung` (:393–432) reads (AREA_A among them).
_M7_NAMES = ("LADDER", "CONTROL_RES", "DEG_ROUTE", "S2_BAND", "S3_BAND", "A", "B",
             "K_GEL", "ROBIN_TAG", "AREA_B", "AREA_A", "G_AREA_BAND", "SERIES_READ",
             "_shell_submesh", "_robin_outer", "_allsum", "_source_dg2", "_solve",
             "_coeffs", "_theta_shell", "_p2", "_q_shell", "_errors", "_volume",
             "_vertex_radii", "_rmin_rmax", "_l2_readout", "_ball")
LADDER = m7.LADDER
CONTROL_RES = m7.CONTROL_RES
DEG_ROUTE = m7.DEG_ROUTE
S2_BAND = m7.S2_BAND
S3_BAND = m7.S3_BAND
A = m7.A
B = m7.B
K_GEL = m7.K_GEL
ROBIN_TAG = m7.ROBIN_TAG
AREA_B = m7.AREA_B
AREA_A = m7.AREA_A
G_AREA_BAND = m7.G_AREA_BAND

V_SHELL = 4.0 / 3.0 * np.pi * (B**3 - A**3)
GEO_QMETA = {"quadrature_degree": 4}  # detJ degree 3 (cells) / 2 (facets) on P2: exact

# I40 affine P2 hole, `-n 2` (20261025T123915Z_IMP-7-step2-W2.log:61-70).
# Width differs (-n 4 here): compare rates and ratios, never digits.
AFFINE_P2 = {LADDER[0]: (1.3382e-3, 1.3695e-3),
             LADDER[1]: (6.8842e-4, 7.0965e-4),
             LADDER[2]: (3.4273e-4, 3.5242e-4)}
AFFINE_P2_PAIRWISE = (1.936, 2.043)  # W2:66-67
AFFINE_P2_READ = (4.4953e-2, 0.394957, 0.327547)  # W2:65 (finest affine rung)
AFFINE_BETA0 = 6.2240e-2  # W2:70 (affine LADDER[1])


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _check_names():
    missing = [n for n in _M7_NAMES if not hasattr(m7, n)]
    assert not missing, f"m7 lacks {missing}"
    assert not any(n.startswith("test_") for n in _M7_NAMES)


def _shell_geo(sub, sft):
    """(volume, area r=b, area r=a, n facets r=b, n facets r=a), all reduced."""
    comm = sub.comm
    one = fem.Constant(sub, default_scalar_type(1.0))
    vol = m7._allsum(comm, one * ufl.dx(domain=sub, metadata=GEO_QMETA))
    ds = ufl.Measure("ds", domain=sub, subdomain_data=sft, metadata=GEO_QMETA)
    area_b = m7._allsum(comm, one * ds(ROBIN_TAG))
    area_a = m7._allsum(comm, one * ds(2))
    nf_own = sub.topology.index_map(sub.topology.dim - 1).size_local
    idx, vals = sft.indices, sft.values
    own = idx < nf_own
    nb = comm.allreduce(int(np.sum(vals[own] == ROBIN_TAG)), op=MPI.SUM)
    na = comm.allreduce(int(np.sum(vals[own] == 2)), op=MPI.SUM)
    return vol, area_b, area_a, nb, na


def _affine_print(res):
    """Affine (flat) shell geometry, printed only; the flat mesh enters only
    the flat cache."""
    msh, ct, _, _ = m7._ball(res)
    sub = m7._shell_submesh(msh, ct)
    sft = m7._robin_outer(sub, A, B)
    vol, ab, aa, nb, na = _shell_geo(sub, sft)
    _say(f"IMP7-3b (Geo) [hole affine] res={res} sub_cells="
         f"{sub.topology.index_map(3).size_global} cmap={_cmap_degree(sub)} "
         f"Vshell/(4/3π(b³-a³))-1={vol / V_SHELL - 1:.4e} "
         f"facets r=b {nb} area/4πb²-1={ab / AREA_B - 1:.4e}; "
         f"facets r=a {na} area/4πa²-1={aa / AREA_A - 1:.4e} (printed)")


def _curved_hole_rung(res, q_fn=None, check=False):
    """`m7._hole_rung`'s body (:393–432) with the mesh from `_ball_curved(res)`
    and rt = DEG_ROUTE[2] (DG3 source, P3 comparand, quadrature 6)."""
    q_fn = m7._q_shell if q_fn is None else q_fn
    comm = MPI.COMM_WORLD
    rt = DEG_ROUTE[2]
    msh, ct, ft, t_build = _ball_curved(res)
    c2 = ct.find(2)
    n2 = comm.allreduce(int(np.sum(c2 < msh.topology.index_map(3).size_local)), op=MPI.SUM)
    sub = m7._shell_submesh(msh, ct)
    deg = _cmap_degree(sub)
    assert deg == 2, f"curved hole sub-mesh cmap degree {deg} (parent {_cmap_degree(msh)})"
    sft = m7._robin_outer(sub, A, B)
    nsub = sub.topology.index_map(3).size_global
    vol, area_b, area_a, nb, na = _shell_geo(sub, sft)
    if check:
        srmin, srmax = m7._rmin_rmax(
            comm, m7._vertex_radii(sub, np.arange(sub.topology.index_map(3).size_local,
                                                  dtype=np.int32)))
        _say(
            f"IMP7-3b G-hole-c res={res} sub_cells={nsub} parent_tag2={n2} cmap={deg} "
            f"shell_rmin/a-1={srmin / A - 1:.3e} shell_rmax/b-1={srmax / B - 1:.3e} "
            f"Vshell/(4/3π(b³-a³))-1={vol / V_SHELL - 1:.4e} "
            f"facets r=b {nb} robin_area/4πb²-1={area_b / AREA_B - 1:.4e}; "
            f"facets r=a {na} r=a_area/4πa²-1={area_a / AREA_A - 1:.4e} (printed)"
        )
        assert nsub == n2, f"(G-hole-c) hole cells {nsub} != parent tag-2 {n2}"
        assert srmin >= A * (1 - 1e-9), f"(G-hole-c) hole vertex at r={srmin!r} < a"
        assert srmax <= B * (1 + 1e-9), f"(G-hole-c) hole vertex at r={srmax!r} > b"
        assert abs(area_b / AREA_B - 1) <= G_AREA_BAND, (
            f"(G-hole-c) Robin area {area_b / AREA_B - 1:.3e}")
    q_f = m7._source_dg2(sub, None, q_fn, dg=rt["q_dg"])
    sol, t_solve = m7._solve(sub, sft, K_GEL, q_f, degree=2)
    c0 = m7._coeffs(0.0)
    theta = m7._p2(sub, lambda x: m7._theta_shell(x, c0), p=rt["comparand"])
    nan_loc = int(np.count_nonzero(~np.isfinite(theta.x.array)))
    nan = comm.allreduce(nan_loc, op=MPI.SUM)
    assert nan == 0, f"hole comparand holds {nan} NaN dofs"
    e_all, e_in, _ = m7._errors(sub, None, sol.temperature, theta, B, qmeta=rt["qmeta"])
    h_eff = (vol / nsub) ** (1.0 / 3.0)
    bal = abs(sol.convected_power - sol.source_power) / sol.source_power
    return dict(sol=sol, t_build=t_build, t_solve=t_solve, nsub=nsub, h_eff=h_eff,
                e_all=e_all, e_in=e_in, bal=bal)


def _rung_report(res, r):
    aw, ai = AFFINE_P2[res]
    _say(f"IMP7-3b [hole curved P2] res={res} cells={r['nsub']} h_eff={r['h_eff']:.5e} "
         f"build={r['t_build']:.1f}s solve={r['t_solve']:.1f}s "
         f"whole_rel_L2={r['e_all']:.4e} interior_rel_L2={r['e_in']:.4e} "
         f"balance={r['bal']:.3e}")
    _say(f"IMP7-3b [hole] res={res} curved/affine ratio (I40 -n 2: whole {aw:.4e} "
         f"interior {ai:.4e}): whole {r['e_all'] / aw:.4f} interior {r['e_in'] / ai:.4f}")
    assert r["bal"] <= S3_BAND, f"(S3h) balance {r['bal']:.3e} > {S3_BAND} at res={res}"


def _pairwise(r0, r1, res0, res1, i):
    rw = np.log(r0["e_all"] / r1["e_all"]) / np.log(r0["h_eff"] / r1["h_eff"])
    ri = np.log(r0["e_in"] / r1["e_in"]) / np.log(r0["h_eff"] / r1["h_eff"])
    _say(f"IMP7-3b (Rate) [hole curved P2] pairwise {res0}->{res1}: whole {rw:.3f} "
         f"interior {ri:.3f} (affine P2 I40 whole {AFFINE_P2_PAIRWISE[i]:.3f}) "
         f"— band whole >= RATE_BAR_P2_CURVED {RATE_BAR_P2_CURVED}")
    assert rw >= RATE_BAR_P2_CURVED, (
        f"(Rate-c) [hole] curved P2 pairwise {res0}->{res1} {rw:.3f} < {RATE_BAR_P2_CURVED}")
    return rw


def _readout(r, res):
    _say(f"IMP7-3b l=2 readout on curved res={res}; affine P2 (I40, finest affine rung): "
         f"R2_h={AFFINE_P2_READ[0]:.4e} pole={AFFINE_P2_READ[1]:.6f} "
         f"equator={AFFINE_P2_READ[2]:.6f} K")
    m7._l2_readout(r["sol"].temperature, f"hole curved P2 res={res}",
                   m7.SERIES_READ["hole"])


# ---- tests --------------------------------------------------------------------


def test_imp7_3b_w0_coding():
    """W0: names, the curved tag-2 sub-mesh's cmap degree, shell volumes."""
    _check_names()
    mc, ctc, _ = MeshGenerator.two_layer_ball_domain(
        inner_radius=A, outer_radius=B, resolution=4e-3, geometry_order=2)
    ma, cta, _ = MeshGenerator.two_layer_ball_domain(
        inner_radius=A, outer_radius=B, resolution=4e-3)
    subc = m7._shell_submesh(mc, ctc)
    suba = m7._shell_submesh(ma, cta)
    dc, da = _cmap_degree(subc), _cmap_degree(suba)
    gc = _shell_geo(subc, m7._robin_outer(subc, A, B))
    ga = _shell_geo(suba, m7._robin_outer(suba, A, B))
    _say(f"IMP7-3b W0 res=4e-3: sub cmap degree curved={dc} affine={da}; sub cells "
         f"curved={subc.topology.index_map(3).size_global} "
         f"affine={suba.topology.index_map(3).size_global}; "
         f"Vshell-1 curved={gc[0] / V_SHELL - 1:.4e} affine={ga[0] / V_SHELL - 1:.4e}; "
         f"facets r=b curved={gc[3]} affine={ga[3]}, r=a curved={gc[4]} affine={ga[4]}")
    assert dc == 2, f"curved hole sub-mesh cmap degree {dc}"
    assert da == 1, f"affine hole sub-mesh cmap degree {da}"


def test_imp7_3c_hole_record_control():
    """Record control (no mesh): the affine floor I40 measured misses the curved bar."""
    _say(f"IMP7-3c record control [hole]: affine P2 pairwise (I40) {AFFINE_P2_PAIRWISE} "
         f"max {max(AFFINE_P2_PAIRWISE):.3f} vs RATE_BAR_P2_CURVED {RATE_BAR_P2_CURVED}")
    assert max(AFFINE_P2_PAIRWISE) < RATE_BAR_P2_CURVED, (
        "record control: the affine floor I40 measured misses the curved bar")


def test_imp7_3b_hole_curved_w1():
    """W1: curved LADDER[0], LADDER[1] at P2, plus (c-β0h) on LADDER[1]."""
    _check_names()
    rows = {}
    for res in LADDER[:2]:
        _affine_print(res)
        rows[res] = _curved_hole_rung(res, check=True)
        _rung_report(res, rows[res])
    _pairwise(rows[LADDER[0]], rows[LADDER[1]], LADDER[0], LADDER[1], 0)
    _readout(rows[LADDER[1]], LADDER[1])
    fin = rows[LADDER[1]]
    _say(f"IMP7-3b (S2h) curved interior rel L2 at res={LADDER[1]} = {fin['e_in']:.4e} "
         f"(band <= S2_BAND {S2_BAND})")
    assert fin["e_in"] <= S2_BAND, f"(S2h) interior {fin['e_in']:.4e} > {S2_BAND}"
    assert CONTROL_RES == LADDER[1]
    c = _curved_hole_rung(CONTROL_RES, q_fn=lambda x: m7._q_shell(x, beta=0.0))
    _say(f"IMP7-3b (c-beta0h) curved res={CONTROL_RES}: interior miss {c['e_in']:.4e} "
         f"(affine P2 I40 {AFFINE_BETA0:.4e}; asserted > {S2_BAND}); "
         f"whole {c['e_all']:.4e}; balance {c['bal']:.3e}; solve={c['t_solve']:.1f}s")
    assert c["bal"] <= S3_BAND, f"(S3h) c-beta0h balance {c['bal']:.3e}"
    assert c["e_in"] > S2_BAND, f"(c-beta0h) control did not fire: {c['e_in']:.4e}"


def test_imp7_3b_hole_curved_w2_finest():
    """W2 (optional): curved LADDER[1] and LADDER[2]; (S2h) on the finest."""
    _check_names()
    rows = {}
    for res in LADDER[1:]:
        if res == LADDER[2]:
            _affine_print(res)
        rows[res] = _curved_hole_rung(res, check=True)
        _rung_report(res, rows[res])
    _pairwise(rows[LADDER[1]], rows[LADDER[2]], LADDER[1], LADDER[2], 1)
    _readout(rows[LADDER[2]], LADDER[2])
    fin = rows[LADDER[2]]
    _say(f"IMP7-3b (S2h) curved finest interior rel L2 = {fin['e_in']:.4e} "
         f"(band <= S2_BAND {S2_BAND})")
    assert fin["e_in"] <= S2_BAND, f"(S2h) interior {fin['e_in']:.4e} > {S2_BAND}"
