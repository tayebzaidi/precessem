"""`WF-10` step 3b (§10 H18 part 2, §9 item 367b) — step 2's gate set at the
refined-phantom root.

Item 367a's fixture (``build_refined_fixture``: the copper F-human hole mesh
with the phantom at 0.0075 m, κ off, r = 0.2 Ω on all 32 gaps, the two-port
90° feed) at its in-model root ``C_IN_REFINED_PF``, with item 358's sense
superposition and ``|B₁⁺|`` gates and item 359's 10 g gates lifted from the
audited commits (2026-10-11 03:00 review ruling (2)).

Brief: ``docs/planning/chunks/WF-10.md`` "§9 item 367b" and the 2026-10-11
03:00 review rulings.

Tests:

* ``test_wf10_step3_lattice_coding`` (L0) — numpy only, always runs.
* ``test_wf10_step3_refined_b1`` — selected by ``FEM_EM_WF10_STEP3_B1=1``.
* ``test_wf10_step3_refined_sar10g`` — selected by ``FEM_EM_WF10_STEP3_SAR=1``.

Both gated tests are heavy, ``-n 8``, never ``-k``.

Scope: degree 1, 64 MHz, κ off, one copper hole fixture, one phantom rung at
0.0075 m; c* is a lattice argmax; no h-statement, no tuning / absolute SAR /
compliance claim, not degree 2, not AED.
"""

from __future__ import annotations

import dataclasses
import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.validation.test_birdcage_b1_quadrature import _mirror_xy
from tests.validation.test_wf7_f_human_sar10g_hotspot import (
    CONTAINMENT_MARGIN_M,
    LATTICE_R_M,
)
from tests.validation.test_wf10_f_human_copper_esr_root import (
    ESR_OHM,
    JSON_DIR,
    Y0_CODING_BAND,
    _active_y_in,
    _thevenin_feed_impedance,
    _z_rc,
)
from tests.validation.test_wf10_f_human_refined_reroot import ROOT_JSON, build_refined_fixture
from tests.validation.test_wf10_f_human_tuned_b1_sar import _sense_weights, _signed_k_shares
from tests.validation.test_wf10_f_human_tuned_sar10g import (
    AZ17_RECORD_DEG,
    NOMINAL_HALF_H_M,
    NOMINAL_R_M,
    _containment,
    _lattice,
    _plane_distance,
)

STEP3_B1_ENV = "FEM_EM_WF10_STEP3_B1"
STEP3_SAR_ENV = "FEM_EM_WF10_STEP3_SAR"
#: C_in′ on the 0.0075 m phantom, rounded to 4 decimals, typed from item 367a's
#: W1 log (`20261011T124038Z_WF-10-step3a-W1.log:6318`: "C_in′ … = 118.231230 pF").
C_IN_REFINED_PF = 118.2312
NODE_JSON = os.path.join(JSON_DIR, f"h18-node-{C_IN_REFINED_PF:.4f}.json")

# 0.015 m readings (items 358 / 359), printed beside, never asserted.
ITEM358_MIRROR = 1.9647e-03  # 20261010T170356Z_WF-10-step2-W1.log:6335
ITEM358_MISPAIRED = 0.9548
ITEM358_C16_WORST = 0.0080
ITEM359_SPREAD = {0.015: 0.167, 0.030: 0.210, 0.045: 0.186}  # …step2s-W1.log:6375
ITEM359_SAR_PER_W = 5.607e-02
ITEM359_PARTNER = 0.0738
ITEM359_PEAK_TO_MEAN = 1.309
ITEM359_MIRROR_STAR = 0.006413
ITEM359_MIRROR_L2 = 0.012789
ITEM359_CONTROL_MISS = 36.75
BALL_S_REFINED_SOLID = 2.67  # 20261010T111003Z_GEO-35-s2bi-W2.log:10831
BALL_S_DEFAULT_HOLE = 3.44
# Predictions written first (brief), printed only.
PREDICTED_SPREAD = (0.03, 0.10)
PREDICTED_SAR_PER_W_BAND = 0.15
PREDICTED_PARTNER = 0.05
PREDICTED_BLIND = 7.0e-3
PREDICTED_SREP = 1.0e-10


def test_wf10_step3_lattice_coding():
    """(L0), 359's assertions re-run through the imported helpers."""
    from fem_em_solver.post.sar import averaging_ball_radius
    from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_mass_averaged_sar_standard_masses import TEN_GRAM_KG

    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    for az17 in (AZ17_RECORD_DEG, 0.0, 33.7):
        lat = _lattice(az17, 0.0)
        assert len(lat) == 18
        pts = np.array([p for *_x, p in lat])
        az_m = az17 + 45.0
        img = _mirror_xy(pts, az_m)
        rows = pts[:12].reshape(4, 3, 3)
        irow = img[:12].reshape(4, 3, 3)
        d_pair = max(
            np.max(np.abs(irow[0] - rows[1])),
            np.max(np.abs(irow[1] - rows[0])),
            np.max(np.abs(irow[2] - rows[3])),
            np.max(np.abs(irow[3] - rows[2])),
        )
        dist = _plane_distance(pts, az_m)
        bar = min(LATTICE_R_M) * np.sin(np.radians(45.0))
        worst = max(
            _containment(c, (0.0, 0.0, 0.0), NOMINAL_R_M, NOMINAL_HALF_H_M, a10)
            for c in np.vstack([pts, img])
        )
        print(
            f"[WF-10 step3b L0] az17 {az17}: 18 points; mirror pairing {d_pair:.3e} m "
            f"(<= 1e-12 ASSERTED); min plane distance {dist.min() * 1e3:.4f} mm (>= "
            f"{bar * 1e3:.4f} mm ASSERTED); worst containment {worst * 1e3:.4f} mm (<= 0 "
            f"ASSERTED, a10 {a10 * 1e3:.4f} mm, margin {CONTAINMENT_MARGIN_M})"
        )
        assert d_pair <= 1e-12, d_pair
        assert dist.min() >= bar - 1e-12, dist.min()
        assert worst <= 0.0, worst
    lat = _lattice(AZ17_RECORD_DEG, 0.0)
    assert sum(1 for x in lat if x[2] == 0.0) == 12


def build_refined_quadrature(comm):
    """367a's fixture, the node pair at ``C_IN_REFINED_PF`` with fields kept, and
    the two sense superpositions.  Never a ``test_*``.

    The superposition / partition / (Lin) block and the JSON record are lifted
    from 63bccf11:tests/validation/test_wf10_f_human_tuned_b1_sar.py lines
    414–474 (``build_tuned_quadrature``, item 358), with the node constant and
    the ``h18-`` JSON name changed.
    """
    import ufl
    from dolfinx import fem

    from fem_em_solver.ports.lumped import sheet_terminal_current
    from tests.validation.test_birdcage_power_identity import (
        DISCRETE_IDENTITY_RTOL,
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_th14_birdcage_copper import _surface_loss_w
    from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ

    z0 = float(REFERENCE_IMPEDANCE_OHM)

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-10 step3b] {msg}", flush=True)

    root = None
    if comm.rank == 0:
        with open(ROOT_JSON) as fh:
            root = json.load(fh)
    root = comm.bcast(root, root=0)
    say(
        f"C_IN_REFINED_PF {C_IN_REFINED_PF} pF beside h18-root.json full precision "
        f"{root['c_in_prime_pf']!r} pF (bracket {root['bracket_pf']}); diff "
        f"{C_IN_REFINED_PF - root['c_in_prime_pf']:+.3e} pF"
    )
    assert abs(C_IN_REFINED_PF - root["c_in_prime_pf"]) <= 5e-5, root["c_in_prime_pf"]

    fx = build_refined_fixture(comm)
    ctx, times, order, feeds = fx["ctx"], fx["times"], fx["order"], fx["feeds"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    acct = {"mesh": msh, "facet_tags": tags_f}
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    c_f = C_IN_REFINED_PF * 1e-12
    z_rc = _z_rc(c_f, FREQUENCY_HZ)
    z_th = _thevenin_feed_impedance(z_rc, z0)

    r = fx["node"](C_IN_REFINED_PF, keep_fields=True)
    times["node_pair"] = r["t_pair"]
    s_node, cur_node, specs_n, pn = r["s_node"], r["currents"], r["specs"], r["partition"]
    p_rc_node, blind_node = r["rel_p_rc"], r["blind"]
    anchors = {"p_rc_node": p_rc_node, "blind_node": blind_node, "dz": r["dz"]}
    e17, e21 = r["e17"], r["e21"]
    sigma = fx["sigma"]
    assert e17.function_space.mesh is e21.function_space.mesh is msh
    assert e17.x.array.shape == e21.x.array.shape
    assert str(e17.function_space.ufl_element()) == str(e21.function_space.ufl_element())
    assert (
        e17.function_space.dofmap.index_map.size_local
        == e21.function_space.dofmap.index_map.size_local
    )

    def _vol(sigma_f, e, measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    def _partition(p):
        p_vol = p["p_phantom"] + p["p_rest"]
        rel = abs(p["p_src"] - p["p_sheets"] - p_vol - p["p_surf"]) / abs(p["p_src"])
        blind = abs(p["p_src"] - p["p_sheets"] - p_vol) / abs(p["p_src"])
        return rel, blind

    # ---- the two senses: superposition, partition, (Lin), Y_in (358 :414–463) --
    w = _sense_weights()
    specs_by = {sp.port_id: sp for sp in specs_n}
    e_s, cur_s, y_in, parts, lin = {}, {}, {}, {}, {}
    t0 = time.perf_counter()
    for sense in ("ccw", "cw"):
        w17, w21 = w[sense]
        ef = fem.Function(e17.function_space, name=f"E_{sense}")
        ef.x.array[:] = w17 * e17.x.array + w21 * e21.x.array
        ef.x.scatter_forward()
        e_s[sense] = ef
        wf = {feeds[0]: w17, feeds[1]: w21}
        cur_s[sense] = w17 * cur_node[0] + w21 * cur_node[1]
        sheets = {
            pid: (
                dataclasses.replace(
                    specs_by[pid].sheet(driven=True), source_voltage_v=complex(wf[pid])
                )
                if pid in feeds
                else specs_by[pid].sheet(driven=False)
            )
            for pid in order
        }
        cur_read = np.array(
            [sheet_terminal_current(msh, tags_f, sheets[pid], ef, comm) for pid in order],
            dtype=np.complex128,
        )
        r_lin = float(np.max(np.abs(cur_read - cur_s[sense])) / np.max(np.abs(cur_s[sense])))
        y_cur = [
            complex(1.0 / z_rc + cur_read[i] / (wf[pid] - z_th * cur_read[i]))
            for i, pid in enumerate(feeds)
        ]
        y_ref = _active_y_in(s_node, (w17, w21), z0)
        r_y = float(max(abs(a - b) for a, b in zip(y_cur, y_ref)) / max(abs(b) for b in y_ref))
        y_in[sense] = y_ref
        p = {
            "p_src": float(sum(_source_power_w(acct, sheets[pid], ef, omega) for pid in feeds)),
            "p_sheet_by_id": {
                pid: float(_sheet_field_dissipation_w(acct, sheets[pid], ef, omega))
                for pid in order
            },
        }
        p["p_sheets"] = float(sum(p["p_sheet_by_id"].values()))
        p_total = _vol(sigma, ef, dx)
        p["p_phantom"] = _vol(sigma, ef, dx(PHANTOM_CELL_TAG))
        p["p_rest"] = p_total - p["p_phantom"]
        p["p_surf"] = _surface_loss_w(msh, fx["built"]["facets"], fx["z_s"], ef, comm)
        rel_w, blind_w = _partition(p)
        p["rel"], p["blind"] = rel_w, blind_w
        parts[sense] = p
        lin[sense] = (r_lin, r_y)
        say(
            f"[{sense}] w = ({w17:.3f}, {w21:.3f}); (P_w) P_src {p['p_src']:.9e} sumP_sheet "
            f"{p['p_sheets']:.9e} P_phantom {p['p_phantom']:.9e} P_rest {p['p_rest']:.9e} P_surf "
            f"{p['p_surf']:.9e} W; rel {rel_w:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); "
            f"P_surf-blind miss {blind_w:.4e} (PREDICTED ≈ {PREDICTED_BLIND:g}, printed)"
        )
        say(
            f"[{sense}] (Lin) currents rel {r_lin:.3e} (band UNIT_WEIGHT_CONTROL_RTOL, ASSERTED); "
            f"Y_in current route vs _active_y_in {r_y:.3e} (band {Y0_CODING_BAND:g}, ASSERTED)"
        )
    times["superposition"] = time.perf_counter() - t0
    anchors["lin"] = lin
    anchors["p_w"] = {s: (parts[s]["rel"], parts[s]["blind"]) for s in parts}

    rec = {
        "c_f": c_f,
        "c_in_pf": C_IN_REFINED_PF,
        "r_ohm": ESR_OHM,
        "freq_hz": float(FREQUENCY_HZ),
        "order": order,
        "feeds": feeds,
        "az_feeds": list(fx["az_feeds"]),
        "az_mirror": fx["az_mirror"],
        "s_node": [[[complex(v).real, complex(v).imag] for v in row] for row in s_node],
        "w": {s: [[complex(x).real, complex(x).imag] for x in w[s]] for s in w},
        "y_in": {s: [[y.real, y.imag] for y in ys] for s, ys in y_in.items()},
        "partitions": {s: parts[s] for s in parts},
        "partition_p17_node": dict(pn),
        "p_rc_node": p_rc_node,
        "blind_node": blind_node,
        "times_s": dict(times),
        "cells": int(fx["cells"]),
        "phantom_cells": int(fx["phantom_cells"]),
    }
    if comm.rank == 0:
        os.makedirs(JSON_DIR, exist_ok=True)
        with open(NODE_JSON, "w") as fh:
            json.dump(rec, fh, indent=1)
        print(f"[WF-10 step3b] wrote {NODE_JSON}", flush=True)
    comm.Barrier()

    return {
        "built": fx["built"],
        "ctx": ctx,
        "order": order,
        "slots": fx["slots"],
        "feeds": feeds,
        "az_feeds": fx["az_feeds"],
        "az_mirror": fx["az_mirror"],
        "w": w,
        "e": e_s,
        "e17": e17,
        "e21": e21,
        "sigma": sigma,
        "s_node": s_node,
        "cur": cur_s,
        "y_in": y_in,
        "parts": parts,
        "p_rc_node": p_rc_node,
        "blind_node": blind_node,
        "times": times,
        "cells": fx["cells"],
        "phantom_cells": fx["phantom_cells"],
        "anchors": anchors,
        "say": say,
        "omega": omega,
        "t_window": fx["t_window"],
    }


@complex_only
def test_wf10_step3_refined_b1():
    """358's ``|B₁⁺|`` sequence (63bccf11 :493–630), lifted; no Z50 block."""
    if os.environ.get(STEP3_B1_ENV, "").strip() != "1":
        pytest.skip(f"{STEP3_B1_ENV} unset: heavy -n 8 F-human window, environment-selected")

    from fem_em_solver.post.faraday import magnetic_flux_density_from_e, project_to_cg1
    from fem_em_solver.post.sar import mean_sar
    from tests.validation.test_birdcage_b1_plus_map import (
        C4_COVARIANCE_BAND,
        PHANTOM_RHO_KG_PER_M3,
        _relative_l2,
        _rotate_z,
        _sample_points,
    )
    from tests.validation.test_birdcage_b1_quadrature import _read_senses
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_drive_superposition import UNIT_WEIGHT_CONTROL_RTOL
    from tests.validation.test_wf10_f_human_tuned_circuit import _k1_shares
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD
    out = build_refined_quadrature(comm)
    say, omega, ctx = out["say"], out["omega"], out["ctx"]
    order, slots, feeds = out["order"], out["slots"], out["feeds"]
    times = out["times"]

    for sense, label in (("cw", "(1,+1j)"), ("ccw", "(1,-1j)")):
        for i, pid in enumerate(feeds):
            y = out["y_in"][sense][i]
            say(f"Y_in a = {label} [{sense}] {pid}: {y.real:+.9e} {y.imag:+.9e}j S (reading)")
    yc, yw = out["y_in"]["ccw"], out["y_in"]["cw"]
    say(
        f"M-paired Y_in spreads: |Y17(ccw) - Y21(cw)|/|Y17(ccw)| {abs(yc[0] - yw[1]) / abs(yc[0]):.4e}; "
        f"|Y21(ccw) - Y17(cw)|/|Y21(ccw)| {abs(yc[1] - yw[0]) / abs(yc[1]):.4e} (reading)"
    )
    for sense in ("ccw", "cw"):
        sh = _signed_k_shares(out["cur"][sense], order, slots)
        un = _k1_shares(out["cur"][sense], order, slots)
        say(
            f"[{sense}] signed k shares: ring+ lead {sh[1.0]['lead']:.6f} lag {sh[1.0]['lag']:.6f}; "
            f"ring- lead {sh[-1.0]['lead']:.6f} lag {sh[-1.0]['lag']:.6f}; unsigned +-1 ring+ "
            f"{un[1.0]:.6f} ring- {un[-1.0]:.6f} (reading)"
        )

    comm.Barrier()
    t0 = time.perf_counter()
    cg1 = {}
    for sense in ("ccw", "cw"):
        b_dg = magnetic_flux_density_from_e(out["e"][sense], omega, name=f"B_{sense}")
        cg1[sense] = project_to_cg1(b_dg, name=f"B_{sense}_cg1")
        del b_dg
    comm.Barrier()
    t_proj = time.perf_counter() - t0
    sweep = {"mesh": ctx["msh"], "cell_tags": ctx["cell_tags"]}
    points = _sample_points(sweep)
    az_m = out["az_mirror"]
    az17 = out["az_feeds"][0]
    comm.Barrier()
    t0 = time.perf_counter()
    ccw_x = _read_senses(cg1["ccw"], points)
    cw_x = _read_senses(cg1["cw"], points)
    cw_mx = _read_senses(cg1["cw"], _mirror_xy(points, az_m))
    cw_m17 = _read_senses(cg1["cw"], _mirror_xy(points, az17))
    rot = {
        m: _read_senses(cg1["ccw"], _rotate_z(points, np.radians(m * 22.5))) for m in range(1, 16)
    }
    comm.Barrier()
    t_eval = time.perf_counter() - t0
    mask = np.logical_and.reduce(
        [ccw_x[2], cw_x[2], cw_mx[2], cw_m17[2]] + [r[2] for r in rot.values()]
    )
    n_pts, n_mask = int(points.shape[0]), int(np.count_nonzero(mask))
    plus_x, minus_x = ccw_x[0], ccw_x[1]
    mirror = _relative_l2(cw_mx[1], plus_x, mask)
    mispaired = _relative_l2(cw_mx[0], plus_x, mask)
    mispaired_pred = _relative_l2(minus_x, plus_x, mask)
    wrong_plane = _relative_l2(cw_m17[1], plus_x, mask)
    c16 = {m: _relative_l2(rot[m][0], plus_x, mask) for m in rot}
    worst_m = max(c16, key=c16.get)
    sense_ratio = float(np.mean(plus_x[mask]) / np.mean(minus_x[mask]))
    cw_ratio = float(np.mean(cw_x[1][mask]) / np.mean(cw_x[0][mask]))
    mean_b = float(np.mean(plus_x[mask]))
    cv_b = float(np.std(plus_x[mask]) / mean_b)
    p_src_ccw = out["parts"]["ccw"]["p_src"]
    say(
        f"B: two CG1 projections {t_proj:.2f} s; {n_pts} sample points, {n_mask} valid on every "
        f"image; evaluations {t_eval:.2f} s"
    )
    say(
        f"(Sense) mean|B1+|_ccw / mean|B1-|_ccw = {sense_ratio:.4f} (> 1 ASSERTED); cw: "
        f"mean|B1-|/mean|B1+| = {cw_ratio:.4f} (printed)"
    )
    say(
        f"(M_B) feed-bisector mirror (az {az_m:.6f} deg): |B1-|_cw(Mx) vs |B1+|_ccw(x) rel L2 "
        f"{mirror:.4e} (band {C4_COVARIANCE_BAND:g}, ASSERTED; 0.015 m phantom {ITEM358_MIRROR:.4e})"
    )
    say(
        f"mis-paired control {mispaired:.4e} beside |B1-|_ccw vs |B1+|_ccw {mispaired_pred:.4e} "
        f"(0.015 m phantom {ITEM358_MISPAIRED:.4f}; printed)"
    )
    say(f"wrong-plane mirror (az17 {az17:.6f} deg): {wrong_plane:.4e} (printed, reading)")
    say(
        "C16 readings: "
        + ", ".join(f"m={m}: {v:.4e}" for m, v in c16.items())
        + f"; worst m={worst_m} {c16[worst_m]:.4e} (0.015 m phantom {ITEM358_C16_WORST:.4f}; "
        f"printed, never asserted)"
    )
    say(
        f"mean|B1+|_ccw {mean_b:.6e} T, CV {cv_b:.4e} (reading); mean|B1+|_ccw/sqrt(P_src) "
        f"{mean_b / np.sqrt(p_src_ccw):.6e} T/sqrt(W)"
    )

    sar = mean_sar(
        out["e"]["ccw"],
        sigma=out["sigma"],
        rho=PHANTOM_RHO_KG_PER_M3,
        cell_tags=ctx["cell_tags"],
        subdomain_ids=PHANTOM_CELL_TAG,
        comm=comm,
    )
    p_ph = out["parts"]["ccw"]["p_phantom"]
    r_sigma = abs(sar["dissipated_power_w"] - p_ph) / abs(p_ph)
    say(
        f"(Sigma_P) P_phantom {p_ph:.12e} W vs mean_sar dissipated {sar['dissipated_power_w']:.12e} "
        f"W, rel {r_sigma:.3e} (band {EXACT_IDENTITY_RTOL:g}, ASSERTED)"
    )
    say(
        f"SAR (ccw): whole-phantom mean {sar['mean_sar_w_per_kg']:.6e} W/kg, per W of P_src "
        f"{sar['mean_sar_w_per_kg'] / p_src_ccw:.6e} W/kg/W; P_phantom/P_src {p_ph / p_src_ccw:.6e}"
    )

    _m, _rss_per, rss_sum = _rss_gib(comm)
    say(
        f"price: build {times['build']:.1f} s, node pair {times['node_pair']:.1f} s, "
        f"superposition + partitions {times['superposition']:.1f} s, B projections {t_proj:.1f} s "
        f"+ evaluations {t_eval:.1f} s; RSS summed {rss_sum:.2f} GiB; window elapsed "
        f"{time.perf_counter() - out['t_window']:.1f} s at -n {comm.size}"
    )

    for sense in ("ccw", "cw"):
        r_lin, r_y = out["anchors"]["lin"][sense]
        rel_w, _b = out["anchors"]["p_w"][sense]
        assert rel_w <= DISCRETE_IDENTITY_RTOL, f"(P_w) {sense} {rel_w:.3e} — never re-band"
        assert r_lin <= UNIT_WEIGHT_CONTROL_RTOL, f"(Lin) currents {sense} {r_lin:.3e}"
        assert r_y <= Y0_CODING_BAND, f"(Lin) Y_in {sense} {r_y:.3e}"
    assert r_sigma <= EXACT_IDENTITY_RTOL, f"(Sigma_P) {r_sigma:.3e}"
    assert n_mask > 0, "no sample point valid on every image"
    assert sense_ratio > 1.0, f"(Sense) red: B1- dominant under ccw ({sense_ratio:.4f})"
    assert mirror <= C4_COVARIANCE_BAND, f"(M_B) mirror {mirror:.4e} — never re-band"


@complex_only
def test_wf10_step3_refined_sar10g():
    """359's 10 g sequence (7ef44c35 :128–362), lifted; (S_rep) against W2's JSON."""
    if os.environ.get(STEP3_SAR_ENV, "").strip() != "1":
        pytest.skip(f"{STEP3_SAR_ENV} unset: heavy -n 8 F-human window, environment-selected")

    from fem_em_solver.post.sar import (
        averaging_ball_radius,
        build_density_field,
        mass_averaged_sar,
        mean_sar,
    )
    from tests.mesh.helpers import global_cell_tag_set
    from tests.validation.test_birdcage_b1_plus_map import (
        C4_COVARIANCE_BAND,
        PHANTOM_RHO_KG_PER_M3,
        _relative_l2,
    )
    from tests.validation.test_birdcage_f_human_rung import F_HUMAN_RING_RADIUS
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_mass_averaged_sar_standard_masses import (
        QUADRATURE_DEGREE,
        TEN_GRAM_KG,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_port_drive_superposition import UNIT_WEIGHT_CONTROL_RTOL
    from tests.validation.test_wf7_f_human_sar10g_hotspot import (
        CONDUCTOR_CLEARANCE_M,
        COVERAGE_BALL_PAD_M,
        FIRST_BALL_STOP_S,
        _phantom_geometry,
    )
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD

    # ---- W2's record, read BEFORE the builder (which rewrites the file) -------
    assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: W3 needs W2's record"
    rec = None
    if comm.rank == 0:
        with open(NODE_JSON) as fh:
            rec = json.load(fh)
    rec = comm.bcast(rec, root=0)
    assert "c_in_pf" in rec, "the JSON is module A's node record, not W2's"
    assert abs(rec["c_in_pf"] - C_IN_REFINED_PF) == 0.0, rec["c_in_pf"]
    s_rec = np.array(
        [[complex(a, b) for a, b in row] for row in rec["s_node"]], dtype=np.complex128
    )

    out = build_refined_quadrature(comm)
    say, ctx, times = out["say"], out["ctx"], out["times"]
    msh, cell_tags, sigma = ctx["msh"], ctx["cell_tags"], out["sigma"]

    s_node = np.asarray(out["s_node"], dtype=np.complex128)
    r_srep = _residual(s_node, s_rec)
    say(
        f"(S_rep) s_node vs W2's JSON at {C_IN_REFINED_PF} pF: residual {r_srep:.6e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED; PREDICTED <= {PREDICTED_SREP:g}); JSON "
        f"az_feeds {rec['az_feeds']} vs run {list(out['az_feeds'])}"
    )

    centre, radius, half_height, lo, hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    z_c = float(centre[2])
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    az17 = float(out["az_feeds"][0])
    az_m = float(out["az_mirror"])
    say(
        f"phantom bbox lo {lo} hi {hi}: centre ({centre[0]:.6e}, {centre[1]:.6e}, "
        f"{centre[2]:.6e}) m, R {radius * 1e3:.4f} mm, H/2 {half_height * 1e3:.4f} mm; 10 g "
        f"a = {a10 * 1e3:.4f} mm; az17 {az17:.9f} deg, mirror plane {az_m:.9f} deg"
    )
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh,
        PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG},
    )

    fields = {"ccw": out["e"]["ccw"], "cw": out["e"]["cw"], "P17": out["e17"], "P21": out["e21"]}
    containment = []
    walls = []

    def ball(key, c, label):
        c = tuple(float(v) for v in c)
        containment.append((label, c, _containment(c, centre, radius, half_height, a10)))
        t = time.perf_counter()
        r = mass_averaged_sar(
            fields[key],
            sigma=sigma,
            rho=rho_field,
            center=c,
            radius=a10,
            comm=comm,
            quadrature_degree=QUADRATURE_DEGREE,
        )
        r["wall_s"] = comm.allreduce(time.perf_counter() - t, op=MPI.MAX)
        walls.append(r["wall_s"])
        say(
            f"  ball {key:<3s} {label:<30s} at ({c[0] * 1e3:8.3f}, {c[1] * 1e3:8.3f}, "
            f"{c[2] * 1e3:8.3f}) mm: SAR10g {r['averaged_sar_w_per_kg']:.9e} W/kg, mass "
            f"{r['mass_kg'] * 1e3:.6f} g, {r['wall_s']:.2f} s"
        )
        return r

    comm.Barrier()
    t_sar = time.perf_counter()
    lattice = _lattice(az17, z_c)
    ccw = {}
    reduced = False
    for idx, (label, _r, dz, _d, c) in enumerate(lattice):
        if reduced and dz != 0.0:
            continue
        ccw[label] = ball("ccw", c, label)
        if idx == 0 and ccw[label]["wall_s"] > FIRST_BALL_STOP_S:
            reduced = True
            say(
                f"REDUCED LATTICE: first ball {ccw[label]['wall_s']:.2f} s > "
                f"{FIRST_BALL_STOP_S:.0f} s — z = z_c rows only (pre-registered, not a red)"
            )
    coords = {label: np.array(c) for label, *_x, c in lattice}
    used = [label for label, *_x in lattice if label in ccw]
    key_star = max(used, key=lambda k: ccw[k]["averaged_sar_w_per_kg"])
    sar_star = ccw[key_star]["averaged_sar_w_per_kg"]
    c_star = coords[key_star]
    cw = {}
    for label in used:
        mc = _mirror_xy(coords[label][None, :], az_m)[0]
        cw[label] = ball("cw", mc, f"M({label})")
    mc_star = _mirror_xy(c_star[None, :], az_m)[0]
    ctl17_c = ball("P17", c_star, "control P17 c*")["averaged_sar_w_per_kg"]
    ctl17_m = ball("P17", mc_star, "control P17 Mc*")["averaged_sar_w_per_kg"]
    ctl21_m = ball("P21", mc_star, "control P21 Mc*")["averaged_sar_w_per_kg"]
    coverage_radius = float(np.hypot(radius, half_height) + COVERAGE_BALL_PAD_M)
    conductor_limit = F_HUMAN_RING_RADIUS - CONDUCTOR_CLEARANCE_M
    t = time.perf_counter()
    whole = mass_averaged_sar(
        fields["ccw"],
        sigma=sigma,
        rho=rho_field,
        center=tuple(float(v) for v in centre),
        radius=coverage_radius,
        comm=comm,
        quadrature_degree=QUADRATURE_DEGREE,
    )
    tagged = mean_sar(
        fields["ccw"],
        sigma=sigma,
        rho=PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        subdomain_ids=PHANTOM_CELL_TAG,
        comm=comm,
    )
    t_cov = comm.allreduce(time.perf_counter() - t, op=MPI.MAX)
    t_sar = comm.allreduce(time.perf_counter() - t_sar, op=MPI.MAX)
    power_miss = abs(whole["dissipated_power_w"] / tagged["dissipated_power_w"] - 1.0)
    mass_miss = abs(whole["mass_kg"] / (PHANTOM_RHO_KG_PER_M3 * tagged["volume_m3"]) - 1.0)
    p_ph = out["parts"]["ccw"]["p_phantom"]
    p_src = out["parts"]["ccw"]["p_src"]
    r_sigma = abs(tagged["dissipated_power_w"] - p_ph) / abs(p_ph)

    sar_c = np.array([ccw[k]["averaged_sar_w_per_kg"] for k in used])
    sar_m = np.array([cw[k]["averaged_sar_w_per_kg"] for k in used])
    mask = np.ones(len(used), dtype=bool)
    mirror_star = abs(cw[key_star]["averaged_sar_w_per_kg"] - sar_star) / sar_star
    mirror_l2 = _relative_l2(sar_m, sar_c, mask)
    ctl_miss = abs(ctl17_m - ctl17_c) / ctl17_c
    ctl_partner = abs(ctl21_m - ctl17_c) / ctl17_c
    spreads = {}
    for r in LATTICE_R_M:
        vals = [
            ccw[lab]["averaged_sar_w_per_kg"]
            for lab, rr, dz, _d, _c in lattice
            if rr == r and dz == 0.0
        ]
        spreads[r] = (max(vals) - min(vals)) / np.mean(vals)
    worst_cont = max(x[2] for x in containment)
    _m, _rss_per, rss_sum = _rss_gib(comm)
    sar_per_w = sar_star / p_src

    if comm.rank == 0:
        print(
            f"[WF-10 step3b] ccw 10 g table ({len(used)} points"
            f"{', REDUCED' if reduced else ''}); cw at the feed-bisector image:",
            flush=True,
        )
        for k in used:
            d = abs(cw[k]["averaged_sar_w_per_kg"] / ccw[k]["averaged_sar_w_per_kg"] - 1.0)
            print(
                f"[WF-10 step3b]   {k:<26s} ccw {ccw[k]['averaged_sar_w_per_kg']:.9e} cw(M) "
                f"{cw[k]['averaged_sar_w_per_kg']:.9e} W/kg, rel {d * 100:.4f}%"
                f"{'   <- c*' if k == key_star else ''}",
                flush=True,
            )
        print(
            f"[WF-10 step3b] c* = {key_star} ({c_star[0] * 1e3:.3f}, {c_star[1] * 1e3:.3f}, "
            f"{c_star[2] * 1e3:.3f}) mm, SAR10g(c*) {sar_star:.9e} W/kg (reading)",
            flush=True,
        )
        print(
            f"[WF-10 step3b] (iii-a) mirror at c*: {mirror_star * 100:.4f}% (band "
            f"{C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED; 0.015 m phantom "
            f"{ITEM359_MIRROR_STAR * 100:.4f}%)",
            flush=True,
        )
        print(
            f"[WF-10 step3b] (iii-b) mirror lattice rel L2 over {len(used)} points: "
            f"{mirror_l2 * 100:.4f}% (band {C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED; 0.015 m "
            f"phantom {ITEM359_MIRROR_L2 * 100:.4f}%)",
            flush=True,
        )
        print(
            "[WF-10 step3b] azimuthal spread of the four z_c rows (RECORD; PREDICTED "
            f"{PREDICTED_SPREAD[0] * 100:g}-{PREDICTED_SPREAD[1] * 100:g}%): "
            + ", ".join(
                f"r={r * 1e3:.0f} mm {v * 100:.4f}% (0.015 m: "
                f"{ITEM359_SPREAD.get(round(float(r), 3), float('nan')) * 100:.1f}%)"
                for r, v in spreads.items()
            ),
            flush=True,
        )
        print(
            f"[WF-10 step3b] control (printed): P17 SAR10g(c*) {ctl17_c:.6e}, P17 SAR10g(Mc*) "
            f"{ctl17_m:.6e} W/kg, miss {ctl_miss * 100:.4f}% (0.015 m {ITEM359_CONTROL_MISS:g}%); "
            f"M-partner P21(Mc*) {ctl21_m:.6e} vs P17(c*): {ctl_partner * 100:.4f}% (0.015 m "
            f"{ITEM359_PARTNER * 100:.2f}%; PREDICTED < {PREDICTED_PARTNER * 100:g}%)",
            flush=True,
        )
        print(
            f"[WF-10 step3b] (i) coverage ball {coverage_radius * 1e3:.4f} mm (< conductor "
            f"limit {conductor_limit * 1e3:.1f} mm ASSERTED): power "
            f"{whole['dissipated_power_w']:.12e} vs tagged {tagged['dissipated_power_w']:.12e} W, "
            f"rel {power_miss:.3e}; mass {whole['mass_kg']:.12e} vs rho*V "
            f"{PHANTOM_RHO_KG_PER_M3 * tagged['volume_m3']:.12e} kg, rel {mass_miss:.3e} (band "
            f"{EXACT_IDENTITY_RTOL:g}, ASSERTED); {t_cov:.2f} s",
            flush=True,
        )
        print(
            f"[WF-10 step3b] (Sigma_P) P_phantom {p_ph:.12e} vs mean_sar dissipated "
            f"{tagged['dissipated_power_w']:.12e} W, rel {r_sigma:.3e} (band "
            f"{EXACT_IDENTITY_RTOL:g}, ASSERTED)",
            flush=True,
        )
        for sense in ("ccw", "cw"):
            rel_w, blind_w = out["anchors"]["p_w"][sense]
            r_lin, r_y = out["anchors"]["lin"][sense]
            print(
                f"[WF-10 step3b] [{sense}] (P_w) {rel_w:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, "
                f"ASSERTED), blind {blind_w:.4e}; (Lin) {r_lin:.3e} / Y {r_y:.3e} (ASSERTED)",
                flush=True,
            )
        print(
            f"[WF-10 step3b] (P_RC) {out['p_rc_node']:.3e} (ASSERTED in the node); node blind "
            f"{out['blind_node']:.4e} (printed)",
            flush=True,
        )
        print(
            f"[WF-10 step3b] (iv) worst containment {worst_cont * 1e3:.4f} mm over "
            f"{len(containment)} balls (<= 0 ASSERTED, measured bbox)",
            flush=True,
        )
        print(
            f"[WF-10 step3b] SAR per watt (ccw): SAR10g(c*)/P_src {sar_per_w:.6e} W/kg/W beside "
            f"0.015 m {ITEM359_SAR_PER_W:.4e} ({sar_per_w / ITEM359_SAR_PER_W - 1.0:+.4%}; "
            f"PREDICTED within +-{PREDICTED_SAR_PER_W_BAND * 100:g}%), /P_phantom "
            f"{sar_star / p_ph:.6e}; whole-phantom mean {tagged['mean_sar_w_per_kg']:.6e} W/kg, "
            f"/P_src {tagged['mean_sar_w_per_kg'] / p_src:.6e}; peak-to-mean "
            f"{sar_star / tagged['mean_sar_w_per_kg']:.4f} (0.015 m {ITEM359_PEAK_TO_MEAN}); "
            f"P_phantom/P_src {p_ph / p_src:.6e}",
            flush=True,
        )
        print(
            f"[WF-10 step3b] PRICE: {len(walls)} 10 g balls + coverage, per ball mean "
            f"{np.mean(walls):.2f} s max {np.max(walls):.2f} s (refined-solid "
            f"{BALL_S_REFINED_SOLID} s, default hole {BALL_S_DEFAULT_HOLE} s); SAR phase "
            f"{t_sar:.1f} s; build {times['build']:.1f} s, node pair {times['node_pair']:.1f} s, "
            f"superposition {times['superposition']:.1f} s; RSS summed {rss_sum:.2f} GiB; window "
            f"{time.perf_counter() - out['t_window']:.1f} s at -n {comm.size}",
            flush=True,
        )

    assert r_srep <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_srep:.6e}"
    for sense in ("ccw", "cw"):
        rel_w, _b = out["anchors"]["p_w"][sense]
        r_lin, r_y = out["anchors"]["lin"][sense]
        assert rel_w <= DISCRETE_IDENTITY_RTOL, f"(P_w) {sense} {rel_w:.3e} — never re-band"
        assert r_lin <= UNIT_WEIGHT_CONTROL_RTOL, f"(Lin) {sense} {r_lin:.3e}"
        assert r_y <= Y0_CODING_BAND, f"(Lin) Y_in {sense} {r_y:.3e}"
    assert r_sigma <= EXACT_IDENTITY_RTOL, f"(Sigma_P) {r_sigma:.3e}"
    assert coverage_radius < conductor_limit, (coverage_radius, conductor_limit)
    assert power_miss <= EXACT_IDENTITY_RTOL, f"(i) power {power_miss:.3e}"
    assert mass_miss <= EXACT_IDENTITY_RTOL, f"(i) mass {mass_miss:.3e}"
    assert worst_cont <= 0.0, f"(iv) containment {worst_cont:.4e} m"
    assert mirror_star <= C4_COVARIANCE_BAND, f"(iii-a) {mirror_star:.4e} — never re-band"
    assert mirror_l2 <= C4_COVARIANCE_BAND, f"(iii-b) {mirror_l2:.4e} — never re-band"
