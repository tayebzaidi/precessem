"""`WF-10` step 5 (§9 item 422, §10 H24): the 10 g lattice under `EX-75`'s
four-port 0 / 90 / 180 / 270 degree drive, gated by C4.

The four-port loading (feeds at slots 0 / 4 / 8 / 12 of the feed ring at Z_th,
the other 28 gaps at Z_RC) maps onto itself under a 90° rotation, and the
``ccw`` weights pick up a common phase, so |E|² — hence SAR10g — is C4-invariant
in CAD; only the mesh breaks it.  359's ``_lattice`` puts the z_c rows at
az₁₇ + 0 / 90 / 180 / 270°, so each radius's four z_c points are one C4 orbit.

Asserted: (C4₄) max over m = 1, 2, 3 of ``_relative_l2`` (row m vs row 0, over
the three radii) ≤ ``C4_COVARIANCE_BAND``; the two-port ccw control on the
same 12 points > band (backed by item 359's logged values); the builder's
anchors, (S_rep), the two-port (P_w) / (Lin); the four-port construction
identities (X4-K) / (X4-P) / (X4-L) and feed azimuths; (Σ_P)₄, (i)₄ coverage;
(iv) containment.  Printed only: the per-radius spread (the record) and its
verdict line against ⅓ of 359's 16.7 / 21.0 / 18.6 %.

Scope: degree 1, 64 MHz, κ off, one copper hole fixture, the 0.015 m phantom,
r = 0.2 Ω on all 32 gaps, C at the two-port root (the four-port drive runs
untuned, as in `EX-75`).  No tuning, f₀, Q, match, S₁₁ / Z_in, absolute SAR,
compliance, h-statement, degree 2 or AED claim.
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
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    _relative_l2,
    _rotate_z,
)
from tests.validation.test_birdcage_b1_quadrature import quadrature_phase_weights
from tests.validation.test_wf7_f_human_sar10g_hotspot import LATTICE_R_M
from tests.validation.test_wf10_f_human_tuned_b1_sar import (
    C_IN_F_HUMAN_PF,
    NODE_JSON,
    _degree_knob,
    _signed_k_shares,
)
from tests.validation.test_wf10_f_human_tuned_sar10g import _containment, _lattice

STEP5_SAR_ENV = "FEM_EM_WF10_STEP5_SAR"

#: Item 359's logged ccw SAR10g on the 12 z_c points (W/kg), rows az₁₇ + 0 / 90
#: / 180 / 270°, columns r = 15 / 30 / 45 mm.  Coding test (L4)(d) only.
#: Source: docs/testing/logs/20261010T171705Z_WF-10-step2s-W1.log.
TWO_PORT_ZC_RECORD = np.array(
    [
        [
            6.830007305e-04,  # :6354  r=15 az17+0
            2.230610345e-03,  # :6355  r=30 az17+0
            4.852895959e-03,
        ],  # :6356  r=45 az17+0
        [
            5.781188640e-04,  # :6357  r=15 az17+90
            1.803115212e-03,  # :6358  r=30 az17+90
            4.026865820e-03,
        ],  # :6359  r=45 az17+90
        [
            6.587520880e-04,  # :6360  r=15 az17+180
            2.230637254e-03,  # :6361  r=30 az17+180
            4.816772104e-03,
        ],  # :6362  r=45 az17+180
        [
            5.882637672e-04,  # :6363  r=15 az17+270
            1.879465733e-03,  # :6364  r=30 az17+270
            4.063626880e-03,
        ],  # :6365  r=45 az17+270
    ]
)
#: 359's spreads (`…step2s-W1.log:6375`) and the brief's stdlib C4 misses.
TWO_PORT_SPREADS = (0.16727, 0.20999, 0.18604)
TWO_PORT_C4_MISSES = {1: 0.17383, 2: 0.00808, 3: 0.16140}
#: ⅓ of 359's spread per radius: the verdict threshold (printed only).
VERDICT_FRACTION = 1.0 / 3.0
# Predictions written first (brief), printed only.
PREDICTED_C4 = (0.3e-2, 3.0e-2)
PREDICTED_SPREAD4 = (0.5e-2, 6.0e-2)
PREDICTED_TWO_PORT_REPRO = 1.0e-9


def _c4_misses(rows):
    """``{m: rel L2 of row m vs row 0}`` for m = 1, 2, 3; rows shape (4, 3)."""
    rows = np.asarray(rows, dtype=float)
    mask = np.ones(rows.shape[1], dtype=bool)
    return {m: _relative_l2(rows[m], rows[0], mask) for m in (1, 2, 3)}


def _spreads(rows):
    """359's (max − min)/mean per column (`…sar10g.py` spread block)."""
    rows = np.asarray(rows, dtype=float)
    return (rows.max(axis=0) - rows.min(axis=0)) / rows.mean(axis=0)


def _zc_rows(table, lattice):
    """The (4, 3) z_c table from a ``{label: value}`` dict, in ``_lattice`` order."""
    zc = [label for label, _r, dz, _d, _c in lattice if dz == 0.0]
    assert len(zc) == 12, len(zc)
    return np.array([table[label] for label in zc], dtype=float).reshape(4, 3)


def _spectrum_share(vec16):
    """Lifted from ``examples/mri/18_f_human_four_port_quadrature.py`` at
    ``63e58597`` (``_spectrum_share``): |X_b|² / Σ over the 16 FFT bins."""
    spec = np.fft.fft(vec16)
    p = np.abs(spec) ** 2
    tot = float(np.sum(p))
    return p / tot if tot else np.full(16, np.nan)


def _ring_vector(currents, order, slots, sign):
    """Lifted from ``examples/mri/18_f_human_four_port_quadrature.py`` at ``63e58597``."""
    vec = np.zeros(16, dtype=np.complex128)
    for idx, pid in enumerate(order):
        k, sg = slots[pid]
        if sg == sign:
            vec[k] = currents[idx]
    return vec


def _four_port_drive(out, comm):
    """The four-port drive on the builder's fixture.

    Provenance: lifted from ``examples/mri/18_f_human_four_port_quadrature.py``
    at commit ``63e58597`` — feeds ``:254–265``, solve / superposition / (X4-P) /
    (X4-L) / (X4-K) ``:477–578`` — with the example's prints and figure removed.
    Returns the superposed field, the partition, kinds, residuals, Y_in, feeds,
    azimuths, the superposed currents and times.  Asserts nothing (the caller
    does).
    """
    import ufl
    from dolfinx import fem

    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.validation.test_birdcage_power_identity import (
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ
    from tests.validation.test_wf10_f_human_copper_esr_root import (
        _source_voltage,
        _thevenin_feed_impedance,
        _z_rc,
    )
    from tests.validation.test_wf10_f_human_tuned_circuit import _respec
    from tests.validation.test_wf12_f_human_copper_matrix import COPPER_SIGMA

    omega, ctx, built = out["omega"], out["ctx"], out["built"]
    msh, cell_tags, sigma = ctx["msh"], ctx["cell_tags"], out["sigma"]
    tags_f = ctx["tags_f"]
    order, slots, feeds2 = out["order"], out["slots"], out["feeds"]

    w4 = quadrature_phase_weights([0, 1, 2, 3], "ccw")
    w4 = w4 / w4[0]

    # ---- the four feeds: slots 0 / 4 / 8 / 12 on the feed ring (example :254–265) --
    ref_k, ref_sign = slots[order[0]]
    feeds4 = []
    for k in (0, 4, 8, 12):
        match = [
            pid for pid, (kk, sg) in slots.items() if sg == ref_sign and kk == (ref_k + k) % 16
        ]
        assert len(match) == 1, (k, match)
        feeds4.append(match[0])
    az_by = {f"P{s['ordinal']}": float(s["azimuth_deg"]) for s in built["sheets"]}
    az4 = [az_by[p] for p in feeds4]

    # ---- the four-port solve (example :477–578) --------------------------------------
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    c_f = C_IN_F_HUMAN_PF * 1e-12
    z_rc = _z_rc(c_f, FREQUENCY_HZ)
    z_th = _thevenin_feed_impedance(z_rc, z0)
    v_s = _source_voltage(z_rc, z0)
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    solver = ctx["solver"]
    base_specs = list(ctx["specs"])
    imp4 = {pid: z_rc for pid in order}
    imp4.update({pid: z_th for pid in feeds4})
    specs4 = _respec(base_specs, None, imp4)
    kappa_off = all(sp.width_correction_kappa is None for sp in specs4)
    specs4_by = {sp.port_id: sp for sp in specs4}
    acct = {"mesh": msh, "facet_tags": tags_f}

    e4, cur4, kinds4, t4 = [], [], [], []
    for j, pid_j in enumerate(feeds4):
        held = j > 0
        sheets = [sp.sheet(driven=(sp.port_id == pid_j)) for sp in specs4]
        dsheet = next(s for s in sheets if s.port_id == pid_j)
        linear = [
            lambda test, _s=dsheet: lumped_port_linear_term(
                msh, tags_f, _s, test, omega_rad_per_s=omega
            )
        ]
        comm.Barrier()
        t0 = time.perf_counter()
        if held:
            fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
        else:
            bil = [wall] + [
                (
                    lambda trial, test, _s=s: lumped_port_bilinear_term(
                        msh, tags_f, _s, trial, test, omega_rad_per_s=omega
                    )
                )
                for s in sheets
            ]
            fields = solver.solve(
                current_density=None,
                project_source=False,
                extra_bilinear_terms=bil,
                extra_linear_terms=linear,
            )
        comm.Barrier()
        t4.append(time.perf_counter() - t0)
        kinds4.append("held" if held else "solve")
        e = fields.e_complex
        currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
        cur4.append(np.array([currents[p] for p in order], dtype=np.complex128))
        cp = fem.Function(e.function_space, name=f"E_{pid_j}_4port")
        cp.x.array[:] = e.x.array
        cp.x.scatter_forward()
        e4.append(cp)
        del fields, e, sheets, dsheet, linear

    t_sup = time.perf_counter()
    ef4 = fem.Function(e4[0].function_space, name="E_4port")
    ef4.x.array[:] = sum(complex(w4[j]) * e4[j].x.array for j in range(4))
    ef4.x.scatter_forward()
    del e4
    cur4_sup = sum(complex(w4[j]) * cur4[j] for j in range(4))
    wf = {pid: complex(w4[j]) for j, pid in enumerate(feeds4)}
    sheets_w = {
        pid: (
            dataclasses.replace(specs4_by[pid].sheet(driven=True), source_voltage_v=wf[pid])
            if pid in wf
            else specs4_by[pid].sheet(driven=False)
        )
        for pid in order
    }
    cur_read = np.array(
        [sheet_terminal_current(msh, tags_f, sheets_w[pid], ef4, comm) for pid in order],
        dtype=np.complex128,
    )
    r_x4l = float(np.max(np.abs(cur_read - cur4_sup)) / np.max(np.abs(cur4_sup)))
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)

    def _vol(e, measure):
        form = fem.form(0.5 * sigma * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    p4_src = float(sum(_source_power_w(acct, sheets_w[pid], ef4, omega) for pid in feeds4))
    p4_sheets = float(
        sum(_sheet_field_dissipation_w(acct, sheets_w[pid], ef4, omega) for pid in order)
    )
    p4_total = _vol(ef4, dx)
    p4_phantom = _vol(ef4, dx(PHANTOM_CELL_TAG))
    p4_surf = _surface_loss_w(msh, built["facets"], z_s, ef4, comm)
    r_x4p = abs(p4_src - p4_sheets - p4_total - p4_surf) / abs(p4_src)
    y4 = [
        complex(
            1.0 / z_rc + cur_read[order.index(pid)] / (wf[pid] - z_th * cur_read[order.index(pid)])
        )
        for pid in feeds4
    ]
    t_sup = comm.allreduce(time.perf_counter() - t_sup, op=MPI.MAX)
    parts = {
        "p_src": p4_src,
        "p_sheets": p4_sheets,
        "p_total": p4_total,
        "p_phantom": p4_phantom,
        "p_rest": p4_total - p4_phantom,
        "p_surf": p4_surf,
    }
    return {
        "ef4": ef4,
        "parts": parts,
        "kinds4": kinds4,
        "r_x4p": r_x4p,
        "r_x4l": r_x4l,
        "y4": y4,
        "feeds4": feeds4,
        "feeds2": list(feeds2),
        "az4": az4,
        "cur4": cur4_sup,
        "t4": [comm.allreduce(t, op=MPI.MAX) for t in t4],
        "t_sup": t_sup,
        "kappa_off": kappa_off,
        "z_rc": z_rc,
        "z_th": z_th,
        "v_s": v_s,
        "ref_sign": ref_sign,
    }


def test_wf10_step5_c4_lattice_coding():
    """(L4): lattice rows are C4 images; four-port weights; planted C4 table;
    359's logged two-port table misses C4 (the backed coding control)."""
    # (a) row m = R(m·90°) row 0.
    for az17 in (11.25, 0.0, 33.7):
        pts = np.array([p for *_x, p in _lattice(az17, 0.0)])[:12].reshape(4, 3, 3)
        d = max(
            float(np.max(np.abs(_rotate_z(pts[0], m * np.pi / 2.0) - pts[m]))) for m in (1, 2, 3)
        )
        print(
            f"[WF-10 step5 L4] (a) az17 {az17}: max |R(m*90) row0 - row m| {d:.3e} m "
            f"(<= 1e-12 ASSERTED)"
        )
        assert d <= 1e-12, d
    # (b) four-port weights.
    w4 = quadrature_phase_weights([0, 1, 2, 3], "ccw")
    w4 = w4 / w4[0]
    x4 = np.zeros(16, dtype=np.complex128)
    for n, k in enumerate((0, 4, 8, 12)):
        x4[k] = w4[n]
    f4 = np.fft.fft(x4)
    print(
        f"[WF-10 step5 L4] (b) four-port FFT bin 0 {abs(f4[0]):.3e}, bin 1 {abs(f4[1]):.3e}, "
        f"|bin 15| {abs(f4[15]):.12f} (ASSERTED 0 / 0 / 4)"
    )
    assert abs(f4[0]) <= 1e-12 and abs(f4[1]) <= 1e-12 and abs(abs(f4[15]) - 4.0) <= 1e-12
    # (c) planted C4-invariant table.
    planted = np.tile(np.array([1.3e-3, 2.7e-3, 5.1e-3]), (4, 1))
    c = _c4_misses(planted)
    print(f"[WF-10 step5 L4] (c) planted C4 table misses {c} (<= 1e-15 ASSERTED)")
    assert max(c.values()) <= 1e-15, c
    # (d) 359's logged two-port table: the coding control.
    c2 = _c4_misses(TWO_PORT_ZC_RECORD)
    s2 = _spreads(TWO_PORT_ZC_RECORD)
    print(
        "[WF-10 step5 L4] (d) 359 two-port C4 misses "
        + " / ".join(f"m={m} {v:.5f}" for m, v in c2.items())
        + " (want "
        + " / ".join(f"{v:.5f}" for v in TWO_PORT_C4_MISSES.values())
        + "); spreads "
        + " / ".join(f"{v:.5f}" for v in s2)
        + f" (want {TWO_PORT_SPREADS}); max_m {max(c2.values()):.5f} > band "
        f"{C4_COVARIANCE_BAND:g} (ASSERTED)"
    )
    for m, want in TWO_PORT_C4_MISSES.items():
        assert abs(c2[m] / want - 1.0) <= 1e-3, (m, c2[m], want)
    for got, want in zip(s2, TWO_PORT_SPREADS):
        assert abs(got / want - 1.0) <= 1e-3, (got, want)
    assert max(c2.values()) > C4_COVARIANCE_BAND, c2


@complex_only
def test_wf10_step5_four_port_sar10g():
    if os.environ.get(STEP5_SAR_ENV, "").strip() != "1":
        pytest.skip(f"{STEP5_SAR_ENV} unset: heavy -n 8 F-human window, environment-selected")

    from fem_em_solver.post.sar import (
        averaging_ball_radius,
        build_density_field,
        mass_averaged_sar,
        mean_sar,
    )
    from tests.mesh.helpers import global_cell_tag_set
    from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_birdcage_f_human_rung import F_HUMAN_RING_RADIUS
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
    from tests.validation.test_mass_averaged_sar_standard_masses import (
        QUADRATURE_DEGREE,
        TEN_GRAM_KG,
    )
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import _residual
    from tests.validation.test_port_drive_superposition import UNIT_WEIGHT_CONTROL_RTOL
    from tests.validation.test_wf7_f_human_sar10g_hotspot import (
        CONDUCTOR_CLEARANCE_M,
        COVERAGE_BALL_PAD_M,
        FIRST_BALL_STOP_S,
        _phantom_geometry,
    )
    from tests.validation.test_wf10_f_human_copper_esr_root import Y0_CODING_BAND
    from tests.validation.test_wf10_f_human_tuned_b1_sar import build_tuned_quadrature
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    # ---- 1. degree 1 only -----------------------------------------------------------
    assert _degree_knob() is None, "FEM_EM_WF10_DEGREE must be unset (degree 1 only)"

    # ---- 2. 358's node record, read BEFORE the builder (which rewrites it) ----------
    assert os.path.exists(NODE_JSON), f"{NODE_JSON} absent: item 422 stops (BLOCKED)"
    rec = None
    if rank0:
        with open(NODE_JSON) as fh:
            rec = json.load(fh)
    rec = comm.bcast(rec, root=0)
    assert abs(rec["c_in_pf"] - C_IN_F_HUMAN_PF) == 0.0, rec["c_in_pf"]
    s_rec = np.array(
        [[complex(a, b) for a, b in row] for row in rec["s_node"]], dtype=np.complex128
    )

    # ---- 3. build + (S_rep) ----------------------------------------------------------
    out = build_tuned_quadrature(comm, with_z50=False)
    say, ctx, times = out["say"], out["ctx"], out["times"]
    msh, cell_tags, sigma = ctx["msh"], ctx["cell_tags"], out["sigma"]
    s_node = np.asarray(out["s_node"], dtype=np.complex128)
    r_srep = _residual(s_node, s_rec)
    say(
        f"(S_rep) s_node vs {os.path.basename(NODE_JSON)} at {C_IN_F_HUMAN_PF} pF: residual "
        f"{r_srep:.6e} (band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED)"
    )

    # ---- 4. geometry -----------------------------------------------------------------
    centre, radius, half_height, lo, hi = _phantom_geometry(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    z_c = float(centre[2])
    a10 = averaging_ball_radius(mass_kg=TEN_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    az17 = float(out["az_feeds"][0])
    say(
        f"phantom bbox lo {lo} hi {hi}: centre xy offset from the coil axis "
        f"({centre[0] * 1e3:.6f}, {centre[1] * 1e3:.6f}) mm, |.| "
        f"{np.hypot(centre[0], centre[1]) * 1e3:.6f} mm (printed; the C4 axis is the coil "
        f"axis); z_c {z_c:.6e} m; R {radius * 1e3:.4f} mm, H/2 {half_height * 1e3:.4f} mm; "
        f"10 g a = {a10 * 1e3:.4f} mm; az17 {az17:.9f} deg"
    )
    all_tags = global_cell_tag_set(msh, cell_tags)
    rho_field = build_density_field(
        msh,
        PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        density_map={int(t): 0.0 for t in sorted(all_tags) if int(t) != PHANTOM_CELL_TAG},
    )

    # ---- 5. memory: keep only ccw ------------------------------------------------------
    e_ccw = out["e"]["ccw"]
    del out["e"]["cw"], out["e17"], out["e21"]

    # ---- 6. four-port drive + (Σ_P)₄ ---------------------------------------------------
    d4 = _four_port_drive(out, comm)
    ef4 = d4["ef4"]
    tagged4 = mean_sar(
        ef4,
        sigma=sigma,
        rho=PHANTOM_RHO_KG_PER_M3,
        cell_tags=cell_tags,
        subdomain_ids=PHANTOM_CELL_TAG,
        comm=comm,
    )
    p4_ph = d4["parts"]["p_phantom"]
    p4_src = d4["parts"]["p_src"]
    r_sigma4 = abs(tagged4["dissipated_power_w"] - p4_ph) / abs(p4_ph)
    az4 = d4["az4"]
    d_az4 = [((az4[(i + 1) % 4] - az4[i]) % 360.0) for i in range(4)]
    az_ok = all(abs(d - 90.0) < AZIMUTH_MATCH_DEG for d in d_az4)
    feeds_ok = d4["feeds4"][0] == d4["feeds2"][0] and d4["feeds4"][1] == d4["feeds2"][1]
    ref_sign = d4["ref_sign"]
    spec4 = {
        s: _spectrum_share(_ring_vector(d4["cur4"], out["order"], out["slots"], s))
        for s in (1.0, -1.0)
    }
    k4 = _signed_k_shares(d4["cur4"], out["order"], out["slots"])

    # ---- 7–8. balls -----------------------------------------------------------------------
    fields = {"4p": ef4, "ccw": e_ccw}
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
            f"  ball {key:<3s} {label:<26s} at ({c[0] * 1e3:8.3f}, {c[1] * 1e3:8.3f}, "
            f"{c[2] * 1e3:8.3f}) mm: SAR10g {r['averaged_sar_w_per_kg']:.9e} W/kg, mass "
            f"{r['mass_kg'] * 1e3:.6f} g, {r['wall_s']:.2f} s"
        )
        return r

    comm.Barrier()
    t_sar = time.perf_counter()
    lattice = _lattice(az17, z_c)
    sar4 = {}
    reduced = False
    for idx, (label, _r, dz, _d, c) in enumerate(lattice):
        if reduced and dz != 0.0:
            continue
        sar4[label] = ball("4p", c, label)["averaged_sar_w_per_kg"]
        if idx == 0 and walls[-1] > FIRST_BALL_STOP_S:
            reduced = True
            say(
                f"REDUCED LATTICE: first ball {walls[-1]:.2f} s > {FIRST_BALL_STOP_S:.0f} s — "
                f"z = z_c rows only (pre-registered, not a red)"
            )
    used4 = [label for label, *_x in lattice if label in sar4]
    key4 = max(used4, key=lambda k: sar4[k])
    coords = {label: np.array(c) for label, *_x, c in lattice}
    sar2 = {}
    for label, _r, dz, _d, c in lattice:
        if dz == 0.0:
            sar2[label] = ball("ccw", c, label)["averaged_sar_w_per_kg"]

    # ---- 9. coverage --------------------------------------------------------------------
    coverage_radius = float(np.hypot(radius, half_height) + COVERAGE_BALL_PAD_M)
    conductor_limit = F_HUMAN_RING_RADIUS - CONDUCTOR_CLEARANCE_M
    t = time.perf_counter()
    whole4 = mass_averaged_sar(
        ef4,
        sigma=sigma,
        rho=rho_field,
        center=tuple(float(v) for v in centre),
        radius=coverage_radius,
        comm=comm,
        quadrature_degree=QUADRATURE_DEGREE,
    )
    t_cov = comm.allreduce(time.perf_counter() - t, op=MPI.MAX)
    t_sar = comm.allreduce(time.perf_counter() - t_sar, op=MPI.MAX)
    power_miss = abs(whole4["dissipated_power_w"] / tagged4["dissipated_power_w"] - 1.0)
    mass_miss = abs(whole4["mass_kg"] / (PHANTOM_RHO_KG_PER_M3 * tagged4["volume_m3"]) - 1.0)

    # ---- readings -------------------------------------------------------------------------
    rows4 = _zc_rows(sar4, lattice)
    rows2 = _zc_rows(sar2, lattice)
    c4_4 = _c4_misses(rows4)
    c4_2 = _c4_misses(rows2)
    c4_max4 = max(c4_4.values())
    c4_max2 = max(c4_2.values())
    sp4 = _spreads(rows4)
    sp2 = _spreads(rows2)
    sp359 = np.array(TWO_PORT_SPREADS)
    repro2 = float(np.max(np.abs(rows2 / TWO_PORT_ZC_RECORD - 1.0)))
    below = [bool(sp4[i] <= VERDICT_FRACTION * _spreads(TWO_PORT_ZC_RECORD)[i]) for i in range(3)]
    verdict = (
        "drive anisotropy"
        if all(below)
        else "lattice resolution next suspect" if not any(below) else "mixed"
    )
    worst_cont = max(x[2] for x in containment)
    mean4 = tagged4["mean_sar_w_per_kg"]
    sar_star4 = sar4[key4]
    _m, _rss_per, rss_sum = _rss_gib(comm)

    if rank0:
        p = d4["parts"]
        print(
            f"[WF-10 step5] four-port feeds {d4['feeds4']} (two-port {d4['feeds2']}), az "
            + " / ".join(f"{a:.9f}" for a in az4)
            + " deg; steps "
            + " / ".join(f"{d:.9f}" for d in d_az4)
            + f" (90 within {AZIMUTH_MATCH_DEG:g} ASSERTED); kappa off {d4['kappa_off']} "
            f"(ASSERTED); Z_RC {d4['z_rc']:.9e}, Z_th {d4['z_th']:.9e} Ohm",
            flush=True,
        )
        print(
            f"[WF-10 step5] (X4-K) kinds {d4['kinds4']} (ASSERTED ['solve','held','held',"
            f"'held']); solves " + ", ".join(f"{t:.2f}" for t in d4["t4"]) + " s",
            flush=True,
        )
        print(
            f"[WF-10 step5] (X4-P) P_src {p['p_src']:.9e} sumP_sheet {p['p_sheets']:.9e} "
            f"P_phantom {p['p_phantom']:.9e} P_rest {p['p_rest']:.9e} P_surf {p['p_surf']:.9e} "
            f"W; rel {d4['r_x4p']:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED; EX-75 "
            f"3.274e-13)",
            flush=True,
        )
        print(
            f"[WF-10 step5] (X4-L) {d4['r_x4l']:.3e} (band {UNIT_WEIGHT_CONTROL_RTOL:g}, "
            f"ASSERTED; EX-75 4.161e-16)",
            flush=True,
        )
        print(
            "[WF-10 step5] Y_in four feeds (printed): "
            + "; ".join(
                f"{pid} {y.real:+.6e} {y.imag:+.6e}j S" for pid, y in zip(d4["feeds4"], d4["y4"])
            ),
            flush=True,
        )
        print(
            f"[WF-10 step5] four-port feed-ring k shares: k = 0 {spec4[ref_sign][0]:.4e} "
            f"(EX-75 4.54e-07), k = -1 (bin 15) {spec4[ref_sign][15]:.6f} (EX-75 0.99945); "
            f"signed lag ring+ {k4[1.0]['lag']:.6f}, ring- {k4[-1.0]['lag']:.6f} (printed)",
            flush=True,
        )
        print(
            f"[WF-10 step5] four-port 10 g table ({len(used4)} points"
            f"{', REDUCED' if reduced else ''}):",
            flush=True,
        )
        for k in used4:
            extra = f"   two-port ccw {sar2[k]:.9e}" if k in sar2 else ""
            print(
                f"[WF-10 step5]   {k:<26s} 4p {sar4[k]:.9e} W/kg{extra}"
                f"{'   <- c*4' if k == key4 else ''}",
                flush=True,
            )
        c = coords[key4]
        print(
            f"[WF-10 step5] c*4 = {key4} ({c[0] * 1e3:.3f}, {c[1] * 1e3:.3f}, {c[2] * 1e3:.3f}) "
            f"mm, SAR10g(c*4) {sar_star4:.9e} W/kg (PREDICTED r = 45 mm, z_c)",
            flush=True,
        )
        print(
            "[WF-10 step5] (C4_4) four-port C4 misses (row m vs row 0, rel L2 over 3 radii): "
            + " / ".join(f"m={m} {v * 100:.4f}%" for m, v in c4_4.items())
            + f"; max {c4_max4 * 100:.4f}% (band {C4_COVARIANCE_BAND * 100:.1f}%, ASSERTED; "
            f"PREDICTED {PREDICTED_C4[0] * 100:g}-{PREDICTED_C4[1] * 100:g}% per m)",
            flush=True,
        )
        print(
            "[WF-10 step5] control: two-port ccw C4 misses "
            + " / ".join(f"m={m} {v * 100:.4f}%" for m, v in c4_2.items())
            + f"; max {c4_max2 * 100:.4f}% (> band {C4_COVARIANCE_BAND * 100:.1f}% ASSERTED; "
            f"359 17.383 / 0.808 / 16.140%)",
            flush=True,
        )
        print(
            f"[WF-10 step5] two-port in-window vs 359's logged values: max rel "
            f"{repro2:.3e} (printed; PREDICTED <= {PREDICTED_TWO_PORT_REPRO:g})",
            flush=True,
        )
        print(
            "[WF-10 step5] spread per radius (the record): "
            + "; ".join(
                f"r={r * 1e3:.0f} mm four-port {sp4[i] * 100:.4f}%, two-port in-window "
                f"{sp2[i] * 100:.4f}%, 359 {sp359[i] * 100:.1f}% (1/3: "
                f"{sp359[i] * 100 / 3:.2f}%)"
                for i, r in enumerate(LATTICE_R_M)
            )
            + f" (PREDICTED four-port {PREDICTED_SPREAD4[0] * 100:g}-"
            f"{PREDICTED_SPREAD4[1] * 100:g}%)",
            flush=True,
        )
        print(
            f"[WF-10 step5] VERDICT (printed only; PREDICTED drive anisotropy): {verdict} "
            f"(per radius <= 1/3 of 359: {below})",
            flush=True,
        )
        print(
            f"[WF-10 step5] (i)4 coverage ball {coverage_radius * 1e3:.4f} mm (< conductor "
            f"limit {conductor_limit * 1e3:.1f} mm ASSERTED): power "
            f"{whole4['dissipated_power_w']:.12e} vs tagged {tagged4['dissipated_power_w']:.12e}"
            f" W, rel {power_miss:.3e}; mass rel {mass_miss:.3e} (band {EXACT_IDENTITY_RTOL:g}, "
            f"ASSERTED); {t_cov:.2f} s",
            flush=True,
        )
        print(
            f"[WF-10 step5] (Sigma_P)4 P_phantom {p4_ph:.12e} vs mean_sar dissipated "
            f"{tagged4['dissipated_power_w']:.12e} W, rel {r_sigma4:.3e} (band "
            f"{EXACT_IDENTITY_RTOL:g}, ASSERTED)",
            flush=True,
        )
        rel_w, _b = out["anchors"]["p_w"]["ccw"]
        r_lin, r_y = out["anchors"]["lin"]["ccw"]
        print(
            f"[WF-10 step5] two-port [ccw] (P_w) {rel_w:.3e}; (Lin) {r_lin:.3e} / Y {r_y:.3e}; "
            f"(P_RC) {out['p_rc_node']:.3e} (ASSERTED)",
            flush=True,
        )
        print(
            f"[WF-10 step5] (iv) worst containment {worst_cont * 1e3:.4f} mm over "
            f"{len(containment)} balls (<= 0 ASSERTED)",
            flush=True,
        )
        print(
            f"[WF-10 step5] four-port per watt: SAR10g(c*4)/P_src,4 {sar_star4 / p4_src:.6e} "
            f"W/kg/W (two-port 5.607e-02); peak-to-mean {sar_star4 / mean4:.4f} (two-port "
            f"1.309); P_phantom/P_src,4 {p4_ph / p4_src:.6e} (two-port 0.0947)",
            flush=True,
        )
        print(
            f"[WF-10 step5] PRICE: {len(walls)} 10 g balls + coverage, per ball mean "
            f"{np.mean(walls):.2f} s max {np.max(walls):.2f} s; SAR phase {t_sar:.1f} s; build "
            f"{times['build']:.1f} s, node pair {times['node_pair']:.1f} s, superposition "
            f"{times['superposition']:.1f} s; four-port solves {sum(d4['t4']):.1f} s, "
            f"partitions {d4['t_sup']:.1f} s; RSS summed {rss_sum:.2f} GiB; window "
            f"{time.perf_counter() - out['t_window']:.1f} s at -n {comm.size}",
            flush=True,
        )

    # ---- asserts ---------------------------------------------------------------------------
    assert r_srep <= REUSE_REPRODUCTION_RTOL, f"(S_rep) {r_srep:.6e}"
    for sense in ("ccw", "cw"):
        rel_w, _b = out["anchors"]["p_w"][sense]
        r_lin, r_y = out["anchors"]["lin"][sense]
        assert rel_w <= DISCRETE_IDENTITY_RTOL, f"(P_w) {sense} {rel_w:.3e}"
        assert r_lin <= UNIT_WEIGHT_CONTROL_RTOL, f"(Lin) {sense} {r_lin:.3e}"
        assert r_y <= Y0_CODING_BAND, f"(Lin) Y_in {sense} {r_y:.3e}"
    assert d4["kappa_off"], "kappa must be off on every sheet"
    assert az_ok, d_az4
    assert feeds_ok, (d4["feeds4"], d4["feeds2"])
    assert d4["kinds4"] == ["solve", "held", "held", "held"], d4["kinds4"]
    assert d4["r_x4p"] <= DISCRETE_IDENTITY_RTOL, f"(X4-P) {d4['r_x4p']:.3e}"
    assert d4["r_x4l"] <= UNIT_WEIGHT_CONTROL_RTOL, f"(X4-L) {d4['r_x4l']:.3e}"
    assert r_sigma4 <= EXACT_IDENTITY_RTOL, f"(Sigma_P)4 {r_sigma4:.3e}"
    assert coverage_radius < conductor_limit, (coverage_radius, conductor_limit)
    assert power_miss <= EXACT_IDENTITY_RTOL, f"(i)4 power {power_miss:.3e}"
    assert mass_miss <= EXACT_IDENTITY_RTOL, f"(i)4 mass {mass_miss:.3e}"
    assert worst_cont <= 0.0, f"(iv) containment {worst_cont:.4e} m"
    assert c4_max2 > C4_COVARIANCE_BAND, f"control: two-port C4 max {c4_max2:.4e} inside band"
    assert c4_max4 <= C4_COVARIANCE_BAND, f"(C4_4) {c4_4} — never re-band"
