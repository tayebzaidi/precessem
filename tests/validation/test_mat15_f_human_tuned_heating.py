"""`MAT-15` step 5 (§9 item 361, §10 P17, provisional) -- step 3's thermal
identities on the coil-driven SAR of H9′'s tuned copper F-human.

Fixture: item 358's ``build_tuned_quadrature(comm, with_z50=False)`` -- the
copper F-human hole mesh (369 387 cells), degree-1 EM, κ off, r = 0.2 Ω on all
32 gaps, C_in,F-human = 118.2437 pF, P17 / P21 driven through Z_th, two senses
ccw = (1, −1j) and cw = (1, +1j) as superpositions of one node pair, 64 MHz.
Thermal: step 1's P1 transient (`test_mat15_coil_driven_heating_phantom.py`)
on the cylindrical 0.015 m phantom's tag-3 submesh, step 1's gel constants and
time grid **in seconds** (t* = 1 254 s; not rescaled to F-human's t_R).

Normalisation (ruling 361-1): one G = 1 W / P_acc,ccw for every thermal run,
P_acc,s := P_src,s − P_sheet,s(P17) − P_sheet,s(P21) -- the power the two
Thevenin feeds deliver across their sheets into the coil network.  Each feed
gap's own r + C sits inside Z_th, so its ESR loss is on the source side of this
cut; that term is printed beside P_acc, as are per-W-of-P_src figures.

Anchors (every band imported, none restated): (S_rep) ≤ REUSE_REPRODUCTION_RTOL;
builder (A0) / (P_RC) / node NEG re-asserted; (P_w), (Lin) both senses; (Σ_P)
≤ EXACT_IDENTITY_RTOL; (V) ≤ V_BAND against 359's logged volume; (B) per field
≤ B_BAND; (E) every run ≤ E_BAND; (E′) ccw / cw ≤ E_BAND; (S2) ccw / cw ≤
S2_BAND; (M_θ-a) / (M_θ-b) feed-bisector mirror of 10 g θ̄ ≤ C4_COVARIANCE_BAND;
(iv) containment ≤ 0.  Negative controls (rule (e)): (d) doubled source,
asserted only > D_CONTROL_FLOOR; (f′) single-feed P17 at c*_θ vs Mc*_θ,
asserted only > F_CONTROL_FLOOR.  Printed only: M-partner P21(Mc*), mis-paired
ccw(Mc*), P_surf-blind misses.

Not claimed: ASTM box phantom, (8e) closed, tissue / perfusion / σ(T), any
absolute temperature or compliance figure, a closed form, h-convergence,
degree 2, f₀ / Q / match, C16 invariance, ANS-9.

Opt-in knob ``FEM_EM_MAT15_STEP5=1``; unset ⇒ the gated tests skip; any other
value raises.  Writes nothing outside the harness log (the builder rewrites
358's node JSON, which is read first).

Additive knob ``FEM_EM_MAT15_STEP5_SERIES=1`` (`MAT-20` step 3, §9 item 380),
honoured only with ``FEM_EM_MAT15_STEP5=1``; any other value raises.  When on,
after every step-5 print the ccw CN run is repeated with every step kept and
rank 0 writes the volume-mean, 10 g ball-mean (at c*_θ) and 5 mm probe θ(t)
series to the gitignored ``logs/mat15_s5/theta-series-ccw.json``.  Unset ⇒
every print and number above is unchanged and nothing is written.

Step 6a (§9 item 413, §10 P22a): ``FEM_EM_MAT15_STEP5_DEGREE`` -- unset /
empty ⇒ the degree-1 path above, untouched; ``"2"`` ⇒ degree-2 EM, honoured only
together with H20a's ``FEM_EM_WF10_DEGREE=2`` (the builder's switch); a
mismatch either way, or any other value (``"1"`` included), raises.  At degree 2
the series goes to ``logs/mat15_s5/theta-series-ccw-d2.json`` and (S_rep) is a
degree-1 record comparand, printed only (``STEP5_DEGREE2_GATES``).  Two more
tests: ``test_mat15_step6a_degree_knob_coding`` (K0, numpy only) and
``test_mat15_step6a_build_only`` (K1, ``FEM_EM_MAT15_STEP6A=1``, build-only).

Step 6b (§9 item 515, §10 P22b): the node follows H28's ``FEM_EM_WF10_NODE``
(``_node_pf()``, read after the step-6a guard): ``d2root`` (honoured only with
both degree knobs at ``"2"``) moves Z_RC / Z_th, the (S_rep) comparand path
(``h9p-node-115.5209-d2.json``), the prints and the series payload to
``C_IN_F_HUMAN_D2_PF``; unset ⇒ every number above bitwise.  Adds the printed-only
Y_in and SAR10g(c*_theta) lines, ``test_mat15_step6b_node_coding`` (K0b, numpy only)
and ``test_mat15_step6b_build_only`` (K1b, ``FEM_EM_MAT15_STEP6B=1``, build-only).
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.post.sar import averaging_ball_radius, mean_sar
from fem_em_solver.thermal.probe import ProbeBall, probe_series
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    PHANTOM_RHO_KG_PER_M3,
    _relative_l2,
)
from tests.validation.test_birdcage_b1_quadrature import _mirror_xy
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_mass_averaged_sar_standard_masses import (
    QUADRATURE_DEGREE,
    TEN_GRAM_KG,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    B_BAND,
    ROBIN_TAG,
    _allsum,
    _robin_tags,
    _say,
    _to_sub,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import (
    E_BAND,
    RHO_C,
    S2_BAND,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    D_CONTROL_FLOOR,
    DT_1,
    F_CONTROL_FLOOR,
    H_GEL,
    K_GEL,
    N_MAIN,
    T_INF,
    T_R,
    T_STAR,
    V_BAND,
    _be_one,
    _cell_average_dg0_tag,
    _cn,
    _e_prime,
    _max_imag,
    _tag_submesh,
    _transient,
)
from tests.validation.test_mat19_probe_coil_phantom import A_PROBE
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_circuit_layer_field import _residual
from tests.validation.test_port_drive_superposition import UNIT_WEIGHT_CONTROL_RTOL
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE, FREQUENCY_HZ
from tests.validation.test_wf7_f_human_sar10g_hotspot import LATTICE_R_M, _phantom_geometry
from tests.validation.test_wf10_f_human_copper_esr_root import (
    ESR_OHM,
    Y0_CODING_BAND,
    _thevenin_feed_impedance,
    _z_rc,
)
from tests.validation.test_wf10_f_human_tuned_b1_sar import (
    C_IN_F_HUMAN_D2_PF,
    C_IN_F_HUMAN_PF,
    DEGREE_ENV,
    E_RECORD,
    F_RECORD,
    N_DOFS_D2_RECORD,
    NODE_ENV,
    NODE_JSON,
    PREDICTED_RSS_D2_GIB,
    _apply_degree,
    _degree_knob,
    _node_json_path,
    _node_pf,
    build_tuned_quadrature,
)
from tests.validation.test_wf10_f_human_tuned_sar10g import (
    AZ17_RECORD_DEG,
    _containment,
    _lattice,
    _plane_distance,
)
from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR, HOLE_CELL_RECORD
from tests.validation.test_wf12_f_human_hole import _rss_gib

KNOB = "FEM_EM_MAT15_STEP5"
SERIES_KNOB = "FEM_EM_MAT15_STEP5_SERIES"
SERIES_JSON = "logs/mat15_s5/theta-series-ccw.json"

# ---- step 6a (§9 item 413, §10 P22a): the degree knob -----------------------------
DEGREE_KNOB = "FEM_EM_MAT15_STEP5_DEGREE"
STEP6A_KNOB = "FEM_EM_MAT15_STEP6A"
#: At degree 2: which step-5 gates are degree-1 record comparands (printed, not
#: asserted) and which stay asserted.  (V) is the tag-3 mesh volume against 359's
#: logged value (ruling 361-3) -- degree-independent, so it stays asserted.
STEP5_DEGREE2_GATES = {
    "printed": ("(S_rep)",),
    "asserted": (
        "(A0)",
        "(P_RC)",
        "NEG node",
        "(P_w)",
        "(Lin)",
        "(Sigma_P)",
        "(V)",
        "(B)",
        "(E)",
        "(E')",
        "(S2)",
        "(M_theta-a)",
        "(M_theta-b)",
        "(iv)",
        "(d)",
        "(f')",
    ),
}


def _step5_degree():
    """``FEM_EM_MAT15_STEP5_DEGREE`` read at call time: unset / empty ⇒ ``None``
    (degree-1 path), ``"2"`` ⇒ 2, anything else (``"1"`` included) ⇒ ValueError.
    Then guarded: must agree with H20a's ``FEM_EM_WF10_DEGREE`` (both unset, or both
    2) -- one degree, two switches thrown together."""
    raw = os.environ.get(DEGREE_KNOB, "").strip()
    if raw == "":
        mine = None
    elif raw == "2":
        mine = 2
    else:
        raise ValueError(f"{DEGREE_KNOB}={raw!r}: only unset/empty (degree 1) or '2' are accepted")
    theirs = _degree_knob()
    if mine != theirs:
        raise ValueError(
            f"{DEGREE_KNOB}={mine!r} disagrees with {DEGREE_ENV}={theirs!r}: set "
            f"both to '2' or leave both unset"
        )
    return mine


def _series_json_path(degree):
    """``SERIES_JSON`` on the degree-1 path; a ``-d2`` sibling at degree 2, so the
    degree-1 series `MAT-20` step 3 / ``mat:18`` read is never overwritten."""
    if degree is None:
        return SERIES_JSON
    if degree == 2:
        return SERIES_JSON[: -len(".json")] + "-d2.json"
    raise ValueError(degree)


# ---- step 6b (§9 item 515, §10 P22b): the node knob ---------------------------------
STEP6B_KNOB = "FEM_EM_MAT15_STEP6B"


def _step6b_node():
    """``(degree, c_pf)``: the step-6a guard first, then H28's ``_node_pf()``."""
    degree = _step5_degree()
    c_pf = _node_pf()
    if degree is None:
        assert c_pf is C_IN_F_HUMAN_PF, c_pf
    return degree, c_pf


def _s_rep_path(degree, c_pf):
    """The (S_rep) comparand path: ``NODE_JSON`` at degree 1; the builder's own
    ``-d2`` record name for the node at degree 2; anything else ⇒ ValueError."""
    if degree is None:
        return NODE_JSON
    if degree == 2:
        return _node_json_path(2, c_pf=c_pf)
    raise ValueError(degree)


def _node_impedances(c_pf, z0):
    """``(z_rc, z_th)`` at the node ``c_pf`` [pF]: the series RC and its Thevenin."""
    z_rc = _z_rc(c_pf * 1e-12, FREQUENCY_HZ)
    return z_rc, _thevenin_feed_impedance(z_rc, z0)


#: Record comparand (printed and asserted at step 1's V_BAND, ruling 361-3, not a
#: band): the meshed tag-3 volume item 359's coverage line printed as mass / ρ,
#: ρ = 1000 (`20261010T171705Z_WF-10-step2s-W1.log:6377`, 13 digits).
V_PHANTOM_H9P_RECORD_M3 = 2.210693836095e-03
#: Nominal phantom of this fixture (`…step2s-W1.log:6313`), printed only.
NOMINAL_R_M = 0.0642857
NOMINAL_HALF_H_M = 0.0857143
V_CAD_M3 = np.pi * NOMINAL_R_M**2 * 2.0 * NOMINAL_HALF_H_M  # 2.225676e-3, printed only
AREA_CAD_M2 = (
    2.0 * np.pi * NOMINAL_R_M * 2.0 * NOMINAL_HALF_H_M + 2.0 * np.pi * NOMINAL_R_M**2
)  # 9.521e-2, printed only
T_R_H = NOMINAL_R_M**2 * RHO_C / K_GEL  # ≈ 28 791 s, printed only
BI_H = H_GEL * NOMINAL_R_M / K_GEL  # ≈ 1.07, printed only

# Predictions written first (brief, §9 item 361), printed only.
PRED_SUB_CELLS = (2000, 8000)
PRED_V_SUB_VS_CAD = -6.73e-3
PRED_AREA_DEFICIT = (-2e-2, 0.0)
PRED_IDENTITY = 1e-12
PRED_SIGMA_P = 1e-14
PRED_B = 1e-13
PRED_ETA = (0.22, 0.34)
PRED_PACC_SENSE = 1e-2
PRED_KAPPA = (1.2, 2.5)
PRED_S2 = (1.9e-6, 4.0e-6)
PRED_MEAN_THETA_PER_ACC_W = (0.025, 0.048)
PRED_PEAK_TO_MEAN = (1.1, 1.35)
PRED_MIRROR = (0.2e-2, 3e-2)
PRED_Q_MIRROR = (0.5e-2, 3e-2)
PRED_FPRIME_MIN = 1.0
PRED_PARTNER = 0.10
PRED_MISPAIRED = (0.03, 0.25)
PRED_BLIND = 7.0e-3


def _enabled():
    raw = os.environ.get(KNOB, "").strip()
    if not raw:
        return False
    if raw != "1":
        raise ValueError(f"{KNOB} must be unset or '1', got {raw!r}")
    return True


def _series_enabled():
    """``FEM_EM_MAT15_STEP5_SERIES``: unset ⇒ False; ``"1"`` ⇒ on only with the step-5
    knob on; anything else raises (`MAT-20` step 3, additive)."""
    raw = os.environ.get(SERIES_KNOB, "").strip()
    if not raw:
        return False
    if raw != "1":
        raise ValueError(f"{SERIES_KNOB} must be unset or '1', got {raw!r}")
    return _enabled()


def _need():
    if not _enabled():
        pytest.skip(f"{KNOB} unset (heavy -n 8 F-human window, gitignored 358/WF-12 inputs)")


def test_mat15_s5_coding(monkeypatch):
    """Ungated (W0): knob parser, feed-bisector mirror coding, time-grid identity."""
    monkeypatch.delenv(KNOB, raising=False)
    assert _enabled() is False
    monkeypatch.setenv(KNOB, "1")
    assert _enabled() is True
    for bad in ("0", "yes", "2"):
        monkeypatch.setenv(KNOB, bad)
        with pytest.raises(ValueError):
            _enabled()
    # MAT-20 step 3's series knob: honoured only with the step-5 knob on.
    monkeypatch.delenv(KNOB, raising=False)
    monkeypatch.delenv(SERIES_KNOB, raising=False)
    assert _series_enabled() is False
    monkeypatch.setenv(SERIES_KNOB, "1")
    assert _series_enabled() is False
    monkeypatch.setenv(KNOB, "1")
    assert _series_enabled() is True
    monkeypatch.delenv(SERIES_KNOB, raising=False)
    assert _series_enabled() is False
    for bad in ("0", "yes", "2"):
        monkeypatch.setenv(SERIES_KNOB, bad)
        with pytest.raises(ValueError):
            _series_enabled()
    monkeypatch.delenv(SERIES_KNOB, raising=False)
    monkeypatch.delenv(KNOB, raising=False)
    lat = _lattice(AZ17_RECORD_DEG, 0.0)
    assert len(lat) == 18
    pts = np.array([p for *_x, p in lat])
    az_m = AZ17_RECORD_DEG + 45.0
    img = _mirror_xy(pts, az_m)
    invol = float(np.max(np.abs(_mirror_xy(img, az_m) - pts)))
    dist = _plane_distance(pts, az_m)
    bar = min(LATTICE_R_M) * np.sin(np.radians(45.0))
    print(
        f"[MAT-15 s5 coding] involution {invol:.3e} m (<= 1e-12 ASSERTED); min plane distance "
        f"{dist.min() * 1e3:.4f} mm (>= {bar * 1e3:.4f} mm ASSERTED); T_STAR = {T_STAR!r} = "
        f"0.2*T_R {0.2 * T_R!r}; DT_1 = {DT_1!r} = 2.5e-6*T_R {2.5e-6 * T_R!r}; "
        f"t*/t_R,H = {T_STAR / T_R_H:.4f} (t_R,H {T_R_H:.1f} s); Bi_H = {BI_H:.4f}"
    )
    assert invol <= 1e-12, invol
    assert dist.min() >= bar - 1e-12, dist.min()
    assert T_STAR == 0.2 * T_R
    assert DT_1 == 2.5e-6 * T_R


def _run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    step5_degree, c_pf = _step6b_node()  # step 6a guard first, then the step-6b node
    _say(
        f"MAT15s5 start: effective EM degree {1 if step5_degree is None else step5_degree} "
        f"({DEGREE_KNOB} {step5_degree if step5_degree is not None else 'unset'}, {DEGREE_ENV} "
        f"{_degree_knob() if _degree_knob() is not None else 'unset'}; fixture DEGREE {DEGREE})"
    )
    _say(
        f"MAT15s5 start: effective node C_in = {c_pf} pF ({NODE_ENV} "
        f"{os.environ.get(NODE_ENV, '').strip() or 'unset'!r}; degree-1 node {C_IN_F_HUMAN_PF} pF, "
        f"degree-2 secant root {C_IN_F_HUMAN_D2_PF} pF)"
    )
    assert int(DEGREE) == 1, DEGREE

    # ---- 1. 358's record, read BEFORE the builder (which rewrites it) -------------
    s_rec = None
    if step5_degree is None:
        assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: item 361 stops (BLOCKED)"
        rec = None
        if comm.rank == 0:
            with open(NODE_JSON) as fh:
                rec = json.load(fh)
        rec = comm.bcast(rec, root=0)
        assert abs(rec["c_in_pf"] - C_IN_F_HUMAN_PF) == 0.0, rec["c_in_pf"]
        s_rec = np.array(
            [[complex(a, b) for a, b in row] for row in rec["s_node"]], dtype=np.complex128
        )
    else:
        node_json_d2 = _s_rep_path(step5_degree, c_pf)  # step 6b: this node's record only
        rec = None
        if comm.rank == 0 and os.path.exists(node_json_d2):
            with open(node_json_d2) as fh:
                rec = json.load(fh)
        rec = comm.bcast(rec, root=0)
        if rec is not None:
            s_rec = np.array(
                [[complex(a, b) for a, b in row] for row in rec["s_node"]], dtype=np.complex128
            )
        else:
            _say(
                f"MAT15s5 (S_rep) degree {step5_degree}: no degree-2 node record "
                f"({node_json_d2} absent); printed only (STEP5_DEGREE2_GATES)"
            )

    # ---- 2. builder ----------------------------------------------------------------
    out = build_tuned_quadrature(comm, with_z50=False)
    t_build = time.perf_counter() - t0
    assert out["c_in_pf"] == c_pf, (out["c_in_pf"], c_pf)
    pfs = comm.allgather(c_pf)
    assert len(set(pfs)) == 1, pfs
    ctx, times = out["ctx"], out["times"]
    msh, cell_tags, sigma = ctx["msh"], ctx["cell_tags"], out["sigma"]
    feeds = list(out["feeds"])
    res = {
        "cells": int(out["cells"]),
        "anchors": out["anchors"],
        "p_rc_node": out["p_rc_node"],
        "blind_node": out["blind_node"],
    }
    s_node = np.asarray(out["s_node"], dtype=np.complex128)
    res["degree"] = step5_degree
    if step5_degree is None:
        res["S_rep"] = _residual(s_node, s_rec)
        _say(
            f"MAT15s5 (S_rep) s_node vs 358's JSON: residual {res['S_rep']:.6e} (band "
            f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED; predicted <= {PRED_IDENTITY:g}); (A0) cells "
            f"{res['cells']} vs {HOLE_CELL_RECORD}; (P_RC) {res['p_rc_node']:.3e}; node NEG "
            f"{res['blind_node']:.4e} = {res['blind_node'] / DISCRETE_IDENTITY_RTOL:.1f}x band "
            f"(bar {BLIND_FACTOR:.0f}x); builder {t_build:.1f} s"
        )
    else:
        res["S_rep"] = None if s_rec is None else _residual(s_node, s_rec)
        s_rep_txt = "n/a (no degree-2 record)" if s_rec is None else f"{res['S_rep']:.6e}"
        _say(
            f"MAT15s5 (S_rep) degree {step5_degree}: s_node vs the -d2 node record read before "
            f"the builder: residual {s_rep_txt} (PRINTED only, STEP5_DEGREE2_GATES); (A0) cells "
            f"{res['cells']} vs {HOLE_CELL_RECORD}; (P_RC) {res['p_rc_node']:.3e}; node NEG "
            f"{res['blind_node']:.4e} = {res['blind_node'] / DISCRETE_IDENTITY_RTOL:.1f}x band "
            f"(bar {BLIND_FACTOR:.0f}x); builder {t_build:.1f} s"
        )

    # ---- 3. power per sense ------------------------------------------------------------
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    z_rc, z_th = _node_impedances(c_pf, z0)  # step 6b: at the builder's node
    p_acc, esr, pvi = {}, {}, {}
    for s in ("ccw", "cw"):
        p = out["parts"][s]
        p_acc[s] = p["p_src"] - sum(p["p_sheet_by_id"][pid] for pid in feeds)
        cur = out["cur"][s]
        ws = out["w"][s]
        v_port = [complex(ws[i]) - z_th * complex(cur[i]) for i in range(2)]
        pvi[s] = 0.5 * sum((v_port[i] * np.conj(complex(cur[i]))).real for i in range(2))
        esr[s] = sum(0.5 * ESR_OHM * abs(v_port[i] / z_rc) ** 2 for i in range(2))
    G = 1.0 / p_acc["ccw"]
    eta = {s: out["parts"][s]["p_phantom"] / p_acc[s] for s in p_acc}
    res.update(G=G, eta=eta, p_acc=p_acc)
    for s in ("ccw", "cw"):
        p = out["parts"][s]
        rel_w, blind_w = out["anchors"]["p_w"][s]
        _say(
            f"MAT15s5 [{s}] P_src {p['p_src']:.12e} W; P_sheet(P17) "
            f"{p['p_sheet_by_id'][feeds[0]]:.6e} P_sheet(P21) {p['p_sheet_by_id'][feeds[1]]:.6e}; "
            f"P_acc {p_acc[s]:.12e} W (1/2 Re sum V_port I* {pvi[s]:.12e}, rel "
            f"{abs(pvi[s] / p_acc[s] - 1.0):.3e}, printed); feed-gap ESR term {esr[s]:.6e} W, P_acc "
            f"+ ESR {p_acc[s] + esr[s]:.12e} W; P_phantom {p['p_phantom']:.12e} P_rest "
            f"{p['p_rest']:.6e} P_surf {p['p_surf']:.6e} W; eta {eta[s]:.6f} (predicted "
            f"{PRED_ETA}); P_phantom/P_src {p['p_phantom'] / p['p_src']:.6f}; (P_w) {rel_w:.3e} "
            f"(band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED), blind {blind_w:.4e} (predicted "
            f"{PRED_BLIND:g}, printed); (Lin) {out['anchors']['lin'][s][0]:.3e} / Y "
            f"{out['anchors']['lin'][s][1]:.3e}"
        )
    _say(
        f"MAT15s5 G = 1/P_acc,ccw = {G:.12e} 1/W (one G, every run); P_acc,cw/P_acc,ccw - 1 = "
        f"{p_acc['cw'] / p_acc['ccw'] - 1.0:+.4e} (predicted |.| <= {PRED_PACC_SENSE:g})"
    )
    # step 6b (Y): Y_in at the node, printed only, never asserted.
    for s, label in (("cw", "(1,+1j)"), ("ccw", "(1,-1j)")):
        for i, pid in enumerate(feeds):
            y = complex(out["y_in"][s][i])
            root = "  <- the P17 (1,+1j) root reading" if (s == "cw" and i == 0) else ""
            _say(
                f"MAT15s5 (Y) Y_in a = {label} [{s}] {pid} at C_in {c_pf} pF: {y.real:+.9e} "
                f"{y.imag:+.9e}j S (printed; |Im| > 1e-3 S => not at the root){root}"
            )

    tagged = mean_sar(
        out["e"]["ccw"],
        sigma=sigma,
        rho=PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        subdomain_ids=PHANTOM_CELL_TAG,
        comm=comm,
    )
    p_ph_ccw = out["parts"]["ccw"]["p_phantom"]
    res["Sigma_P"] = abs(float(np.real(tagged["dissipated_power_w"])) - p_ph_ccw) / abs(p_ph_ccw)
    _say(
        f"MAT15s5 (Sigma_P) P_phantom,ccw {p_ph_ccw:.15e} vs mean_sar "
        f"{float(np.real(tagged['dissipated_power_w'])):.15e} W, rel {res['Sigma_P']:.3e} (band "
        f"{EXACT_IDENTITY_RTOL:g}, ASSERTED; predicted <= {PRED_SIGMA_P:g})"
    )

    # ---- 4. thermal ---------------------------------------------------------------------
    t_th = time.perf_counter()
    fields = {"ccw": out["e"]["ccw"], "cw": out["e"]["cw"], "P17": out["e17"], "P21": out["e21"]}
    dxp = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)(PHANTOM_CELL_TAG)
    p_ph = {}
    for f, e in fields.items():
        p_ph[f] = float(np.real(_allsum(comm, 0.5 * sigma * ufl.inner(e, e) * dxp)))
    for s in ("ccw", "cw"):
        _say(
            f"MAT15s5 P_phantom,{s} _vol {p_ph[s]:.15e} vs partition "
            f"{out['parts'][s]['p_phantom']:.15e} rel "
            f"{abs(p_ph[s] / out['parts'][s]['p_phantom'] - 1.0):.3e} (printed)"
        )
    _say(f"MAT15s5 P_phantom P17 {p_ph['P17']:.12e} P21 {p_ph['P21']:.12e} W (_vol)")

    sub, parent_ids = _tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)
    ft = _robin_tags(sub)
    simap = sub.topology.index_map(3)
    n_sub = simap.size_global
    counts = comm.gather((simap.size_local, simap.num_ghosts), root=0)
    one = fem.Constant(sub, default_scalar_type(1.0))
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    area = float(np.real(_allsum(comm, one * ds)))
    res["V"] = abs(v_sub / V_PHANTOM_H9P_RECORD_M3 - 1.0)
    res["n_sub"] = n_sub
    _say(
        f"MAT15s5 sub cells {n_sub} (predicted {PRED_SUB_CELLS}); per-rank (owned, ghost) "
        f"{counts}"
    )
    _say(
        f"MAT15s5 (V) V_sub {v_sub:.15e} m^3 vs record {V_PHANTOM_H9P_RECORD_M3:.12e}: rel "
        f"{res['V']:.3e} (band {V_BAND:g}, ASSERTED; predicted <= 1e-13); vs CAD {V_CAD_M3:.6e} "
        f"{v_sub / V_CAD_M3 - 1.0:+.4e} (predicted {PRED_V_SUB_VS_CAD:+.2e}); oint ds {area:.9e} "
        f"m^2 vs CAD {AREA_CAD_M2:.6e} {area / AREA_CAD_M2 - 1.0:+.4e} (predicted "
        f"{PRED_AREA_DEFICIT})"
    )
    _say(
        f"MAT15s5 grid: t* = {T_STAR:.3f} s, N = {N_MAIN}, dt_1 = {DT_1:.6e} s (step 1's, "
        f"seconds); t_R,H = {T_R_H:.1f} s, t*/t_R,H = {T_STAR / T_R_H:.4f}; Bi_H = {BI_H:.4f}"
    )

    q, B, kappa = {}, {}, {}
    for f, e in fields.items():
        W, avg = _cell_average_dg0_tag(
            msh, cell_tags, PHANTOM_CELL_TAG, 0.5 * sigma * ufl.inner(e, e)
        )
        avg *= G  # G applied after the form (ruling 361-4): one compiled form for all four
        q[f] = _to_sub(W, avg, sub, parent_ids)
        del W, avg
        dep = G * p_ph[f]
        int_q = float(np.real(_allsum(comm, q[f] * ufl.dx(domain=sub))))
        q_s = float(np.real(_allsum(comm, q[f] * ds)))
        B[f] = abs(int_q - dep) / abs(dep)
        kappa[f] = (q_s / area) / (int_q / v_sub)
        _say(
            f"MAT15s5 (B) [{f}] int q = {int_q:.15e} vs G*P_phantom {dep:.15e}: rel "
            f"{B[f]:.3e} (band {B_BAND:g}, ASSERTED; predicted <= {PRED_B:g}); kappa_s "
            f"{kappa[f]:.4f} (predicted {PRED_KAPPA})"
        )
    res["B"] = B
    res["kappa"] = kappa

    E, Ep, theta, mean_th = {}, {}, {}, {}
    imag = 0.0
    for f in fields:
        t1 = time.perf_counter()
        r = _cn(sub, ft, q[f])
        t_cn = time.perf_counter() - t1
        dep = G * p_ph[f]
        ep, stor, conv = _e_prime(sub, ds, r, dep)
        E[f"cn_{f}"] = r.max_balance_residual
        Ep[f] = ep
        th = fem.Function(r.temperature.function_space)
        th.x.array[:] = r.temperature.x.array - T_INF
        th.x.scatter_forward()
        theta[f] = th
        mean_th[f] = float(np.real(_allsum(comm, th * ufl.dx(domain=sub)))) / v_sub
        imag = max(imag, _max_imag(r.temperature))
        _say(
            f"MAT15s5 CN [{f}] N={N_MAIN}: mean theta(t*) {mean_th[f]:.9e} K per accepted W "
            f"(G of ccw); (E) {r.max_balance_residual:.3e}; (E') stored {stor:.12e} + surface "
            f"{conv:.12e} vs G*P_phantom {dep:.12e}: rel {ep:.3e} (band {E_BAND:g}; asserted on "
            f"ccw/cw); solve {t_cn:.2f} s"
        )
    res["Ep"] = Ep

    S2 = {}
    for s in ("ccw", "cw"):
        dep = G * p_ph[s]
        r1 = _be_one(sub, ft, q[s])
        E[f"be_{s}"] = r1.max_balance_residual
        imag = max(imag, _max_imag(r1.temperature))
        thbar1 = (
            float(np.real(_allsum(comm, (r1.temperature - T_INF) * ufl.dx(domain=sub)))) / v_sub
        )
        ratio = RHO_C * thbar1 / (DT_1 * dep / v_sub)
        leak = float(np.real(_allsum(comm, H_GEL * (r1.temperature - T_INF) * ds))) / dep
        est = H_GEL * DT_1 * area / (RHO_C * v_sub)
        S2[s] = abs(ratio - 1.0)
        _say(
            f"MAT15s5 (S2) [{s}] ratio {ratio:.12f} |ratio-1| {S2[s]:.4e} (band {S2_BAND:g}, "
            f"ASSERTED; predicted {PRED_S2}); deficit {leak:.4e}, estimate {est:.3e}*kappa_s = "
            f"{est * kappa[s]:.4e}; (E) {r1.max_balance_residual:.3e}"
        )
    res["S2"] = S2

    q_two = fem.Function(q["ccw"].function_space)
    q_two.x.array[:] = 2.0 * q["ccw"].x.array
    q_two.x.scatter_forward()
    r_d = _cn(sub, ft, q_two)
    ep_d, _sd, _cd = _e_prime(sub, ds, r_d, G * p_ph["ccw"])
    E["ctl_d"] = r_d.max_balance_residual
    imag = max(imag, _max_imag(r_d.temperature))
    res["ctl_d"] = ep_d
    res["E"] = E
    res["imag"] = imag
    _say(
        f"MAT15s5 control (d) doubled source: (E') vs G*P_phantom,ccw = {ep_d:.6f} (predicted "
        f"1.000, asserted > {D_CONTROL_FLOOR:g}); its (E) {r_d.max_balance_residual:.3e}"
    )
    worst = max(E, key=E.get)
    _say(
        f"MAT15s5 (E) worst over {len(E)} runs {E[worst]:.3e} on {worst} (band {E_BAND:g}, "
        f"ASSERTED); max|Im theta| = {imag:.3e} (predicted 0)"
    )
    t_thermal = time.perf_counter() - t_th

    # ---- 5. ball means (step 3's forms, copied) ----------------------------------------
    # Copied from `main` at base 9c38d4e6, tests/validation/test_mat15_tuned_copper_heating.py
    # lines 251–282 (the ball-mean block inside `_run`), adapted: one θ / q holder
    # per call, centres from item 359's lattice instead of the eight ring points.
    t_b = time.perf_counter()
    centre, radius, half_height, lo, hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    z_c = float(centre[2])
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    centre_c = fem.Constant(sub, np.zeros(3, dtype=default_scalar_type))
    xs = ufl.SpatialCoordinate(sub)
    off = xs - centre_c
    inside = ufl.conditional(ufl.lt(ufl.real(ufl.dot(off, off)), a10**2), 1.0, 0.0)
    dxq = ufl.dx(domain=sub, metadata={"quadrature_degree": int(QUADRATURE_DEGREE)})
    V1 = theta["ccw"].function_space
    Ws = q["ccw"].function_space
    th_ball = fem.Function(V1)
    q_ball = fem.Function(Ws)
    f_vol = fem.form(inside * dxq)
    f_th = fem.form(th_ball * inside * dxq)
    f_q = fem.form(q_ball * inside * dxq)
    v_ball = 4.0 / 3.0 * np.pi * a10**3

    def _red(form):
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    containment = []

    def _ball(field, c, label):
        c = np.array([float(np.real(v)) for v in c])
        containment.append((label, _containment(c, centre, radius, half_height, a10)))
        centre_c.value[:] = c
        vol = _red(f_vol)
        th_ball.x.array[:] = theta[field].x.array
        th_ball.x.scatter_forward()
        q_ball.x.array[:] = q[field].x.array
        q_ball.x.scatter_forward()
        tb = _red(f_th) / vol
        qb = _red(f_q) / vol
        _say(
            f"MAT15s5   ball {field:<3s} {label:<32s} at ({c[0] * 1e3:8.3f}, {c[1] * 1e3:8.3f}, "
            f"{c[2] * 1e3:8.3f}) mm: theta {tb:.9e} K  q {qb:.9e} W/m^3  vol {vol:.6e} "
            f"(4/3 pi a^3 {v_ball:.6e}, {vol / v_ball - 1.0:+.3e})"
        )
        return tb, qb

    az17 = float(out["az_feeds"][0])
    az_m = float(out["az_mirror"])
    _say(
        f"MAT15s5 phantom bbox centre ({centre[0]:.6e}, {centre[1]:.6e}, {centre[2]:.6e}) m, R "
        f"{radius * 1e3:.4f} mm, H/2 {half_height * 1e3:.4f} mm; a10 {a10 * 1e3:.4f} mm; az17 "
        f"{az17:.9f} deg, mirror plane {az_m:.9f} deg"
    )
    lattice = _lattice(az17, z_c)
    labels = [lab for lab, *_x in lattice]
    coords = {lab: np.array(c) for lab, *_x, c in lattice}
    ccw = {lab: _ball("ccw", coords[lab], lab) for lab in labels}
    cw = {
        lab: _ball("cw", _mirror_xy(coords[lab][None, :], az_m)[0], f"M({lab})") for lab in labels
    }
    key_star = max(labels, key=lambda k: ccw[k][0])
    c_star = coords[key_star]
    mc_star = _mirror_xy(c_star[None, :], az_m)[0]
    th17_c = _ball("P17", c_star, "control P17 c*")[0]
    th17_m = _ball("P17", mc_star, "control P17 Mc*")[0]
    th21_m = _ball("P21", mc_star, "control P21 Mc*")[0]
    thccw_m = _ball("ccw", mc_star, "mis-paired ccw Mc*")[0]

    th_c = np.array([ccw[k][0] for k in labels])
    th_m = np.array([cw[k][0] for k in labels])
    q_c = np.array([ccw[k][1] for k in labels])
    q_m = np.array([cw[k][1] for k in labels])
    mask = np.ones(len(labels), dtype=bool)
    th_star = ccw[key_star][0]
    res["M_a"] = abs(cw[key_star][0] - th_star) / th_star
    res["M_b"] = _relative_l2(th_m, th_c, mask)
    q_ma = abs(cw[key_star][1] - ccw[key_star][1]) / ccw[key_star][1]
    q_mb = _relative_l2(q_m, q_c, mask)
    res["fprime"] = abs(th17_m - th17_c) / th17_c
    partner = abs(th21_m - th17_c) / th17_c
    mispaired = abs(thccw_m - th_star) / th_star
    res["cont"] = max(x[1] for x in containment)
    t_balls = time.perf_counter() - t_b

    for k in labels:
        _say(
            f"MAT15s5   {k:<26s} theta ccw {ccw[k][0]:.9e} cw(M) {cw[k][0]:.9e} K rel "
            f"{abs(cw[k][0] / ccw[k][0] - 1.0) * 100:.4f}%; q ccw {ccw[k][1]:.6e} cw(M) "
            f"{cw[k][1]:.6e} rel {abs(cw[k][1] / ccw[k][1] - 1.0) * 100:.4f}%"
            f"{'   <- c*_theta' if k == key_star else ''}"
        )
    _say(
        f"MAT15s5 c*_theta = {key_star} ({c_star[0] * 1e3:.3f}, {c_star[1] * 1e3:.3f}, "
        f"{c_star[2] * 1e3:.3f}) mm; peak 10 g theta-bar per accepted W {th_star:.9e} K "
        f"(predicted at r = 45 mm, z_c, az17); peak-to-mean {th_star / mean_th['ccw']:.4f} "
        f"(predicted {PRED_PEAK_TO_MEAN})"
    )
    _say(
        f"MAT15s5 (M_theta-a) mirror at c*_theta {res['M_a'] * 100:.4f}% (band "
        f"{C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED; predicted {PRED_MIRROR}); q-ball pair "
        f"{q_ma * 100:.4f}% (printed, predicted {PRED_Q_MIRROR})"
    )
    _say(
        f"MAT15s5 (M_theta-b) lattice rel L2 over {len(labels)} points {res['M_b'] * 100:.4f}% "
        f"(band {C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED; predicted {PRED_MIRROR}); q-ball rel "
        f"L2 {q_mb * 100:.4f}% (printed)"
    )
    _say(
        f"MAT15s5 control (f') P17 theta(c*) {th17_c:.6e} vs P17 theta(Mc*) {th17_m:.6e} K: "
        f"miss {res['fprime'] * 100:.4f}% (asserted > {F_CONTROL_FLOOR * 100:.1f}%, predicted >= "
        f"{PRED_FPRIME_MIN * 100:.0f}%)"
    )
    _say(
        f"MAT15s5 M-partner P21(Mc*) {th21_m:.6e} vs P17(c*): {partner * 100:.4f}% (printed, "
        f"predicted <= {PRED_PARTNER * 100:.0f}%); mis-paired ccw(Mc*) vs ccw(c*): "
        f"{mispaired * 100:.4f}% (printed, predicted {PRED_MISPAIRED})"
    )
    _say(
        f"MAT15s5 (iv) worst containment {res['cont'] * 1e3:.4f} mm over {len(containment)} "
        f"balls (<= 0 ASSERTED, measured bbox)"
    )

    # ---- readings per watt ----------------------------------------------------------------
    p_src = out["parts"]["ccw"]["p_src"]
    to_src = p_acc["ccw"] / p_src
    m = mean_th["ccw"]
    _say(
        f"MAT15s5 READING mean theta(t*) [ccw]: {m:.6e} K per accepted W (predicted "
        f"{PRED_MEAN_THETA_PER_ACC_W}); {m * to_src:.6e} K per W of P_src; "
        f"{m / eta['ccw']:.6e} K per deposited W; eta x t*/(rho c V_sub) = "
        f"{eta['ccw'] * T_STAR / (RHO_C * v_sub):.6e} K (no-loss bound); cw mean "
        f"{mean_th['cw']:.6e} K"
    )
    _say(
        f"MAT15s5 READING peak 10 g theta-bar [ccw]: {th_star:.6e} K per accepted W; "
        f"{th_star * to_src:.6e} K per W of P_src; {th_star / eta['ccw']:.6e} K per deposited W"
    )
    sar10 = ccw[key_star][1] / PHANTOM_RHO_KG_PER_M3
    _say(
        f"MAT15s5 READING SAR10g(c*_theta) [ccw] at C_in {c_pf} pF: {sar10:.6e} W/kg per accepted "
        f"W; {sar10 * to_src:.6e} W/kg per W of P_src; {sar10 / eta['ccw']:.6e} W/kg per deposited "
        f"W (q-ball mean / rho, printed)"
    )
    _m, _per, rss = _rss_gib(comm)
    _say(
        f"MAT15s5 PRICE: build {times['build']:.1f} s, node pair {times['node_pair']:.1f} s, "
        f"superposition {times['superposition']:.1f} s; builder total {t_build:.1f} s; thermal "
        f"{t_thermal:.1f} s; balls {t_balls:.1f} s ({len(containment)}); RSS summed {rss:.2f} "
        f"GiB; total {time.perf_counter() - t0:.1f} s at -n {comm.size}"
    )

    # ---- 6. MAT-20 step 3 (§9 item 380): θ(t) series, additive knob ----------------------
    if _series_enabled():
        t_s = time.perf_counter()
        dt = T_STAR / N_MAIN
        steps = tuple(range(N_MAIN + 1))
        r_s = _transient(
            sub,
            ft,
            q["ccw"],
            dt=dt,
            n_steps=N_MAIN,
            theta=0.5,
            startup_half_steps=2,
            snapshot_steps=steps,
        )
        snaps = r_s.snapshots
        assert sorted(snaps) == list(steps), sorted(snaps)
        th_s = fem.Function(snaps[0].function_space)
        f_mean = fem.form(th_s * ufl.dx(domain=sub))
        centre_c.value[:] = np.array([float(np.real(v)) for v in c_star])
        vol_star = _red(f_vol)
        s_mean, s_ball = [], []
        for n in steps:
            th_s.x.array[:] = snaps[n].x.array - T_INF
            th_s.x.scatter_forward()
            th_ball.x.array[:] = th_s.x.array
            th_ball.x.scatter_forward()
            s_mean.append(_red(f_mean) / v_sub)
            s_ball.append(_red(f_th) / vol_star)
        pser = probe_series(snaps, dt, [ProbeBall(tuple(c_star), A_PROBE)])
        s_probe = [float(v) - T_INF for v in pser.values[:, 0]]
        times_s = [n * dt for n in steps]
        _say(
            f"MAT15s5 SERIES ccw CN N={N_MAIN}, every step kept ({len(steps)} snapshots, "
            f"dt {dt:.6e} s); n=N: mean theta {s_mean[-1]:.9e} vs module {mean_th['ccw']:.9e} "
            f"(rel {abs(s_mean[-1] / mean_th['ccw'] - 1.0):.3e}); ball theta-bar(c*) "
            f"{s_ball[-1]:.9e} vs module {th_star:.9e} (rel "
            f"{abs(s_ball[-1] / th_star - 1.0):.3e}); probe a={A_PROBE:g} m {s_probe[-1]:.9e} "
            f"(vol {pser.volumes[0]:.6e} vs nominal {pser.nominal_volumes[0]:.6e}; max Im/Re "
            f"{pser.max_imag_ratio:.3e}); printed, predicted equal to rounding"
        )
        for n in steps:
            _say(
                f"MAT15s5 SERIES n={n:2d} t={times_s[n]:9.3f} s: mean {s_mean[n]:.12e} ball "
                f"{s_ball[n]:.12e} probe {s_probe[n]:.12e} K"
            )
        payload = {
            "source": "MAT-15 step 5 fixture, ccw CN every step kept (MAT-20 step 3, item 380)",
            "times_s": times_s,
            "dt_s": dt,
            "t_inf_c": T_INF,
            "n_main": N_MAIN,
            "G_per_w": G,
            "p_acc_w": {s: float(p_acc[s]) for s in p_acc},
            "p_vi_w": {s: float(pvi[s]) for s in pvi},
            "p_src_w": float(out["parts"]["ccw"]["p_src"]),
            "esr_w": {s: float(esr[s]) for s in esr},
            "p_phantom_ccw_w": float(out["parts"]["ccw"]["p_phantom"]),
            "c_star_label": key_star,
            "c_star_m": [float(np.real(v)) for v in c_star],
            "a10_m": float(a10),
            "a_probe_m": float(A_PROBE),
            "ball_volume_m3": vol_star,
            "probe_volume_m3": float(pser.volumes[0]),
            "theta_mean_k": s_mean,
            "theta_ball10g_k": s_ball,
            "theta_probe_k": s_probe,
            "theta_mean_tstar_module_k": float(mean_th["ccw"]),
            "theta_ball10g_tstar_module_k": float(th_star),
            "cells": int(out["cells"]),
            "n_sub": int(n_sub),
            "ranks": int(comm.size),
            "normalisation": "per accepted W, ruling 361-1: G = 1 / P_acc,ccw",
            "c_in_pf": float(c_pf),
            "node_knob": os.environ.get(NODE_ENV, "").strip(),
            "degree": 1 if step5_degree is None else int(step5_degree),
        }
        series_json = _series_json_path(step5_degree)
        if comm.rank == 0:
            os.makedirs(os.path.dirname(series_json), exist_ok=True)
            with open(series_json, "w") as fh:
                json.dump(payload, fh, indent=1)
        comm.Barrier()
        _say(
            f"MAT15s5 SERIES wrote {series_json} (rank 0); series block "
            f"{time.perf_counter() - t_s:.1f} s"
        )
    return res


@pytest.fixture(scope="module")
def heated():
    _need()
    return _run()


@complex_only
def test_s_rep_source_wiring(heated):
    if heated.get("degree") == 2:
        pytest.skip(
            "(S_rep): degree-1 record comparand, printed only at degree 2 " "(STEP5_DEGREE2_GATES)"
        )
    assert heated["S_rep"] <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {heated['S_rep']:.6e}"


@complex_only
def test_builder_anchors(heated):
    assert abs(heated["cells"] / HOLE_CELL_RECORD - 1.0) <= CELL_COUNT_BAND, heated["cells"]
    assert heated["p_rc_node"] <= DISCRETE_IDENTITY_RTOL, heated["p_rc_node"]
    assert heated["blind_node"] >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, heated["blind_node"]


@complex_only
@pytest.mark.parametrize("sense", ["ccw", "cw"])
def test_p_w_and_lin(heated, sense):
    rel_w, _b = heated["anchors"]["p_w"][sense]
    r_lin, r_y = heated["anchors"]["lin"][sense]
    assert rel_w <= DISCRETE_IDENTITY_RTOL, f"(P_w) {sense} {rel_w:.3e}"
    assert r_lin <= UNIT_WEIGHT_CONTROL_RTOL, f"(Lin) {sense} {r_lin:.3e}"
    assert r_y <= Y0_CODING_BAND, f"(Lin) Y {sense} {r_y:.3e}"


@complex_only
def test_sigma_p(heated):
    assert heated["Sigma_P"] <= EXACT_IDENTITY_RTOL, f"(Sigma_P) {heated['Sigma_P']:.3e}"


@complex_only
def test_v_fixture_tie(heated):
    assert heated["V"] <= V_BAND, f"(V) {heated['V']:.3e}"


@complex_only
@pytest.mark.parametrize("field", ["ccw", "cw", "P17", "P21"])
def test_b_source_transfer(heated, field):
    assert heated["B"][field] <= B_BAND, f"(B) {field} {heated['B'][field]:.3e}"


@complex_only
def test_e_every_run_balances(heated):
    for name, r in heated["E"].items():
        assert r <= E_BAND, f"(E) {name} {r:.3e}"


@complex_only
@pytest.mark.parametrize("sense", ["ccw", "cw"])
def test_e_prime_energy_identity(heated, sense):
    assert heated["Ep"][sense] <= E_BAND, f"(E') {sense} {heated['Ep'][sense]:.3e}"


@complex_only
@pytest.mark.parametrize("sense", ["ccw", "cw"])
def test_s2_calorimetric_slope(heated, sense):
    assert heated["S2"][sense] <= S2_BAND, f"(S2) {sense} {heated['S2'][sense]:.4e}"


@complex_only
def test_m_theta_a_mirror_at_c_star(heated):
    assert heated["M_a"] <= C4_COVARIANCE_BAND, f"(M_theta-a) {heated['M_a']:.4e}"


@complex_only
def test_m_theta_b_mirror_on_lattice(heated):
    assert heated["M_b"] <= C4_COVARIANCE_BAND, f"(M_theta-b) {heated['M_b']:.4e}"


@complex_only
def test_iv_containment(heated):
    assert heated["cont"] <= 0.0, f"(iv) {heated['cont']:.4e} m"


@complex_only
def test_control_d_doubled_source_misses(heated):
    assert heated["ctl_d"] > D_CONTROL_FLOOR, f"control (d) {heated['ctl_d']:.3e} did not miss"


@complex_only
def test_control_f_prime_single_feed_misses(heated):
    r = heated["fprime"]
    assert r > F_CONTROL_FLOOR, f"control (f') {r * 100:.4f}% inside {F_CONTROL_FLOOR * 100:.1f}%"


def test_mat15_step6a_degree_knob_coding(monkeypatch):
    """(K0) step 6a: the parser, the guard, the series path, the gate table."""
    monkeypatch.delenv(DEGREE_KNOB, raising=False)
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    monkeypatch.delenv(NODE_ENV, raising=False)
    assert _step5_degree() is None
    for raw in ("", "  "):
        monkeypatch.setenv(DEGREE_KNOB, raw)
        assert _step5_degree() is None, raw
    monkeypatch.setenv(DEGREE_ENV, "2")
    for raw in ("2", " 2 "):
        monkeypatch.setenv(DEGREE_KNOB, raw)
        assert _step5_degree() == 2, raw
    # Negative control (asserted): every rejected value raises.
    for raw in ("1", "3", "x", "2.0"):
        monkeypatch.setenv(DEGREE_KNOB, raw)
        with pytest.raises(ValueError):
            _step5_degree()
    # The guard: both mismatches raise.
    monkeypatch.setenv(DEGREE_KNOB, "2")
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    with pytest.raises(ValueError, match=DEGREE_ENV):
        _step5_degree()
    monkeypatch.delenv(DEGREE_KNOB, raising=False)
    monkeypatch.setenv(DEGREE_ENV, "2")
    with pytest.raises(ValueError, match=DEGREE_KNOB):
        _step5_degree()
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    p1, p2 = _series_json_path(None), _series_json_path(2)
    print(f"[MAT-15 step6a K0] series JSON: unset -> {p1}; degree 2 -> {p2}")
    assert p1 == SERIES_JSON
    assert (
        p2.endswith("-d2.json")
        and os.path.dirname(p2) == os.path.dirname(SERIES_JSON)
        and p2 != SERIES_JSON
    )
    with pytest.raises(ValueError):
        _series_json_path(1)
    printed, asserted = STEP5_DEGREE2_GATES["printed"], STEP5_DEGREE2_GATES["asserted"]
    print(f"[MAT-15 step6a K0] degree-2 gates: printed {printed}; asserted {asserted}")
    assert "(S_rep)" in printed
    assert not set(printed) & set(asserted)
    for gate in (
        "(V)",
        "(B)",
        "(E)",
        "(E')",
        "(S2)",
        "(M_theta-a)",
        "(M_theta-b)",
        "(iv)",
        "(d)",
        "(f')",
        "(Sigma_P)",
    ):
        assert gate in asserted, gate
    # The records' own arithmetic (exact integers, H20a's records).
    assert N_DOFS_D2_RECORD == 2 * E_RECORD + 2 * F_RECORD == 2_570_024


@complex_only
def test_mat15_step6a_build_only(monkeypatch):
    """(K1) step 6a: one build, both paths' unknown counts; no form compiled."""
    if os.environ.get(STEP6A_KNOB, "").strip() != "1":
        pytest.skip(f"{STEP6A_KNOB} unset: build-only -n 2 F-human window, environment-selected")
    comm = MPI.COMM_WORLD

    def say(msg):
        if comm.rank == 0:
            print(f"[MAT-15 step6a K1] {msg}", flush=True)

    monkeypatch.delenv(DEGREE_KNOB, raising=False)
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    monkeypatch.delenv(NODE_ENV, raising=False)
    assert _step5_degree() is None
    t0 = time.perf_counter()
    out = build_tuned_quadrature(comm, with_z50=False, build_only=True)
    t_build = time.perf_counter() - t0
    ctx = out["ctx"]
    n1, e, f = out["n_dofs"], out["n_edges"], out["n_faces"]
    same = ctx["solver"] is out["solver_built"]
    say(
        f"(A0) hole cells {out['cells']} vs record {HOLE_CELL_RECORD} (band {CELL_COUNT_BAND}, "
        f"ASSERTED in the builder and here)"
    )
    say(
        f"(Deg-1) both knobs unset: degree {out['degree']}, solver is _build_f_human's: {same}; "
        f"unknowns {n1} vs E {e} (record E {E_RECORD}; ASSERTED exact); F {f} (record "
        f"{F_RECORD}); build {t_build:.1f} s"
    )
    assert abs(out["cells"] / HOLE_CELL_RECORD - 1.0) <= CELL_COUNT_BAND, out["cells"]
    assert out["degree"] == 1 and ctx["solver"].degree == 1
    assert same, "(Deg-1) the unset path rebuilt the solver: the knob leaks"
    assert n1 == e, (n1, e)
    out.pop("solver_built")

    monkeypatch.setenv(DEGREE_KNOB, "2")
    monkeypatch.setenv(DEGREE_ENV, "2")
    assert _step5_degree() == 2
    t0 = time.perf_counter()
    solver2, deg2, counts2 = _apply_degree(ctx, comm)
    t_d2 = time.perf_counter() - t0
    n2 = counts2["n_dofs"]
    assert counts2["n_edges"] == e and counts2["n_faces"] == f
    _m, _per, rss_sum = _rss_gib(comm)
    say(
        f"(Deg-2) both knobs 2: degree {deg2} (solver.degree {solver2.degree}); unknowns {n2} "
        f"vs 2E+2F = {2 * e + 2 * f} (record {N_DOFS_D2_RECORD}; ASSERTED exact); degree-2 "
        f"space {t_d2:.1f} s"
    )
    say(
        f"NEG (asserted, rule (e)): n(2) {n2} != n(1) {n1}; n(2)/n(1) = {n2 / n1:.6f} vs "
        f"(2E+2F)/E = {(2 * e + 2 * f) / e:.6f} (exact integer identity ASSERTED)"
    )
    say(
        f"RSS after the degree-2 space, summed {rss_sum:.2f} GiB (PREDICTED <= "
        f"{PREDICTED_RSS_D2_GIB:g}; H20a 1.10 GiB; printed)"
    )
    assert deg2 == 2 and solver2.degree == 2 and ctx["solver"] is solver2
    assert n2 == 2 * e + 2 * f, (n2, e, f)
    assert n2 != n1
    assert n2 * e == n1 * (2 * e + 2 * f), (n1, n2, e, f)


def _closed_form_z_rc(c_pf):
    return complex(ESR_OHM, -1.0 / (2.0 * np.pi * float(FREQUENCY_HZ) * c_pf * 1e-12))


def test_mat15_step6b_node_coding(monkeypatch):
    """(K0b) step 6b: node knob, guard order, (S_rep) path, the (N2) closed form."""
    for name in (NODE_ENV, DEGREE_ENV, DEGREE_KNOB):
        monkeypatch.delenv(name, raising=False)
    deg, c = _step6b_node()
    assert deg is None and c is C_IN_F_HUMAN_PF
    monkeypatch.setenv(NODE_ENV, "d2root")
    monkeypatch.setenv(DEGREE_ENV, "2")
    monkeypatch.setenv(DEGREE_KNOB, "2")
    deg, c = _step6b_node()
    assert (deg, c) == (2, 115.5209) and c is C_IN_F_HUMAN_D2_PF
    # d2root with both degree knobs unset: _node_pf's guard names DEGREE_ENV.
    monkeypatch.delenv(DEGREE_ENV, raising=False)
    monkeypatch.delenv(DEGREE_KNOB, raising=False)
    with pytest.raises(ValueError, match=DEGREE_ENV):
        _step6b_node()
    # d2root with only FEM_EM_WF10_DEGREE=2: the step-6a guard fires first.
    monkeypatch.setenv(DEGREE_ENV, "2")
    with pytest.raises(ValueError, match=DEGREE_KNOB):
        _step6b_node()
    monkeypatch.setenv(DEGREE_KNOB, "2")
    for bad in ("115.5209", "x", "d2"):
        monkeypatch.setenv(NODE_ENV, bad)
        with pytest.raises(ValueError):
            _step6b_node()
    for name in (NODE_ENV, DEGREE_ENV, DEGREE_KNOB):
        monkeypatch.delenv(name, raising=False)
    p_none = _s_rep_path(None, C_IN_F_HUMAN_PF)
    p_d2_old = _s_rep_path(2, C_IN_F_HUMAN_PF)
    p_d2_root = _s_rep_path(2, C_IN_F_HUMAN_D2_PF)
    print(
        f"[MAT-15 step6b K0b] (S_rep) paths: unset {p_none}; degree 2 at {C_IN_F_HUMAN_PF} "
        f"{p_d2_old}; degree 2 at {C_IN_F_HUMAN_D2_PF} {p_d2_root}"
    )
    assert p_none == NODE_JSON
    assert p_d2_old == _node_json_path(2)
    assert p_d2_root.endswith("h9p-node-115.5209-d2.json")
    assert p_d2_root != p_d2_old and p_d2_root != p_none
    with pytest.raises(ValueError):
        _s_rep_path(1, C_IN_F_HUMAN_PF)
    assert _series_json_path(None) == SERIES_JSON
    # (N2) closed form at both nodes.
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    for cc in (C_IN_F_HUMAN_PF, C_IN_F_HUMAN_D2_PF):
        z_rc, z_th = _node_impedances(cc, z0)
        zc = _closed_form_z_rc(cc)
        r1 = abs(z_rc - zc) / abs(zc)
        r2 = abs(z_th - z0 * zc / (z0 + zc)) / abs(z0 * zc / (z0 + zc))
        print(
            f"[MAT-15 step6b K0b] (N2) C {cc} pF: Z_RC {z_rc:.9e} rel {r1:.2e}; Z_th rel "
            f"{r2:.2e} (<= 1e-12 ASSERTED)"
        )
        assert r1 <= 1e-12 and r2 <= 1e-12, (cc, r1, r2)
    # NEG (asserted, closed form): a knob-blind :360 reads Z_RC at 118.2437 pF.
    ratio = (
        _node_impedances(C_IN_F_HUMAN_D2_PF, z0)[0].imag
        / _node_impedances(C_IN_F_HUMAN_PF, z0)[0].imag
    )
    print(
        f"[MAT-15 step6b K0b] NEG Im Z_RC(root)/Im Z_RC(118.2437) = {ratio:.6f} "
        f"(= {C_IN_F_HUMAN_PF / C_IN_F_HUMAN_D2_PF:.6f}; a blind module gives 1)"
    )
    assert abs(ratio - C_IN_F_HUMAN_PF / C_IN_F_HUMAN_D2_PF) <= 1e-12
    assert abs(ratio - 1.0) >= 0.02, ratio


@complex_only
def test_mat15_step6b_build_only():
    """(K1b) step 6b: this module sees ``FEM_EM_WF10_NODE=d2root`` through
    ``mpiexec`` at degree 2; build-only; reads the window's environment."""
    if os.environ.get(STEP6B_KNOB, "").strip() != "1":
        pytest.skip(f"{STEP6B_KNOB} unset: build-only -n 2 F-human window, environment-selected")
    comm = MPI.COMM_WORLD

    def say(msg):
        if comm.rank == 0:
            print(f"[MAT-15 step6b K1b] {msg}", flush=True)

    assert _step6b_node() == (2, C_IN_F_HUMAN_D2_PF), "window environment wrong (wiring red)"
    json_dir = os.path.dirname(_node_json_path(2))
    files = [
        os.path.join(json_dir, f"h9p-node-{c:.4f}-d2.json")
        for c in (C_IN_F_HUMAN_D2_PF, C_IN_F_HUMAN_PF)
    ]

    def _stat(p):
        return (os.path.exists(p), os.path.getmtime(p) if os.path.exists(p) else None)

    before = [_stat(p) for p in files]
    t0 = time.perf_counter()
    out = build_tuned_quadrature(comm, with_z50=False, build_only=True)
    t_build = time.perf_counter() - t0
    after = [_stat(p) for p in files]
    n2, e, f = out["n_dofs"], out["n_edges"], out["n_faces"]
    pfs = comm.allgather(_node_pf())
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    z_rc = _node_impedances(out["c_in_pf"], z0)[0]
    zc = _closed_form_z_rc(C_IN_F_HUMAN_D2_PF)
    r_n2 = abs(z_rc - zc) / abs(zc)
    ratio = z_rc.imag / _node_impedances(C_IN_F_HUMAN_PF, z0)[0].imag
    _m, _per, rss_sum = _rss_gib(comm)
    say(
        f"(A0) cells {out['cells']} vs {HOLE_CELL_RECORD}; (Deg-2) degree {out['degree']}; "
        f"unknowns {n2} vs 2E+2F = {2 * e + 2 * f} (record {N_DOFS_D2_RECORD}; ASSERTED exact)"
    )
    say(
        f"node: c_in_pf {out['c_in_pf']} (want {C_IN_F_HUMAN_D2_PF}, ASSERTED exact); node_json "
        f"{os.path.basename(out['node_json'])}; _s_rep_path {_s_rep_path(2, _node_pf())} "
        f"(== builder's ASSERTED); per-rank _node_pf {pfs}"
    )
    say(
        f"(N2) Z_RC at node {z_rc:.9e} vs closed form rel {r_n2:.2e}; NEG Im ratio vs 118.2437 "
        f"{ratio:.6f} (blind module 1)"
    )
    say(f"-d2 records before {before} after {after} (unchanged ASSERTED)")
    say(
        f"build {t_build:.1f} s; RSS summed {rss_sum:.2f} GiB (PREDICTED <= "
        f"{PREDICTED_RSS_D2_GIB:g}; H28 W1 113.1 s / 1.10 GiB; printed)"
    )
    assert abs(out["cells"] / HOLE_CELL_RECORD - 1.0) <= CELL_COUNT_BAND, out["cells"]
    assert out["degree"] == 2
    assert n2 == 2 * e + 2 * f == N_DOFS_D2_RECORD, (n2, e, f)
    assert out["c_in_pf"] == C_IN_F_HUMAN_D2_PF and out["c_in_pf"] != C_IN_F_HUMAN_PF
    assert os.path.basename(out["node_json"]) == "h9p-node-115.5209-d2.json", out["node_json"]
    assert _s_rep_path(2, _node_pf()) == out["node_json"]
    assert _s_rep_path(2, _node_pf()) != _node_json_path(2)
    assert len(set(pfs)) == 1 and pfs[0] == C_IN_F_HUMAN_D2_PF, pfs
    assert r_n2 <= 1e-12, r_n2
    assert abs(ratio - C_IN_F_HUMAN_PF / C_IN_F_HUMAN_D2_PF) <= 1e-12, ratio
    assert after == before, (before, after)
