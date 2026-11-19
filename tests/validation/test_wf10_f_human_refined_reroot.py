"""`WF-10` step 3a (§10 H18 part 1, §9 item 367a) — the refined-phantom re-root.

H9′'s construction (copper F-human hole mesh, κ off, r = 0.2 Ω on all 32
gaps, the two-port 90° feed) with the phantom meshed at
``PHANTOM_RESOLUTION_M`` (0.0075 m, `WF-7` step 3's constant, imported): the
in-model root on C_in is re-read by three nodes (118.20, 118.2437, 118.30 pF).

Brief: ``docs/planning/chunks/WF-10.md`` "§9 item 367a" and the 2026-10-11
03:00 review rulings beneath it.

Tests:

* ``test_wf10_step3_refined_kwargs_contract`` (K) — no mesh, always runs.
* ``test_wf10_step3_reroot`` (S-root) — selected by ``FEM_EM_WF10_STEP3_ROOT=1``
  (never ``-k``); heavy, ``-n 8``.

``build_refined_fixture`` is module-level so item 367b imports it.

Scope: degree 1, 64 MHz, κ off, one copper hole fixture, one phantom rung at
0.0075 m.  The root is a record at its grid resolution, never asserted; not
f₀ / Q / match / tuning accuracy, not an h-statement, not ``|B₁⁺|`` or SAR
(367b), not degree 2, not AED.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.validation.test_wf7_f_human_sar10g_hotspot import (
    BUILD_STOP_S,
    PHANTOM_RESOLUTION_M,
    RSS_STOP_GIB,
)
from tests.validation.test_wf10_f_human_copper_esr_root import (
    ESR_OHM,
    JSON_DIR,
    NPZ_PATH,
    PREDICTED_SLOPE_FACTOR,
    SENSES,
    _active_y_in,
    _current_route_s,
    _sign_change_roots,
    _source_voltage,
    _thevenin_feed_impedance,
    _z_rc,
)
from tests.validation.test_wf10_f_human_tuned_b1_sar import C_IN_F_HUMAN_PF

STEP3_ROOT_ENV = "FEM_EM_WF10_STEP3_ROOT"
DRIVEN = "P17"
#: The three re-root nodes (brief); the middle one is C_IN_F_HUMAN_PF, imported.
REROOT_NODES_PF = (118.20, C_IN_F_HUMAN_PF, 118.30)
#: Item 354's grid resolution at the default mesh (printed beside the move).
ITEM354_RESOLUTION_PF = 0.10
#: Item 354's four-reading root spread (printed beside this run's).
ITEM354_SPREAD_PF = 0.0378
#: Default-mesh cell / phantom-cell figures, printed beside (A0′), never asserted.
DEFAULT_HOLE_CELLS = 369387  # 20261010T170356Z_WF-10-step2-W1.log (A0)
REFINED_SOLID_PHANTOM_CELLS = 22594  # 20261010T110256Z_GEO-35-s2bi-W1.log:10599
# Predictions written first (brief), printed only.
PREDICTED_MOVE_PF = (-0.03, 0.03)
PREDICTED_IM_Y_MID = 2.0e-4
PREDICTED_CELLS = 402000
PREDICTED_PHANTOM_CELLS = 22600
PREDICTED_BLIND = 7.0e-3
PREDICTED_DZ = 0.23

#: The default (0.015 m) mesh's Y_in at the same C, printed beside, never asserted.
DEFAULT_MESH_Y_IN = {
    # 20261010T093200Z_WF-10-step1e2-W1a.log:6334–6335 (node 118.20 pF)
    118.20: {
        "(1,+1j)": (1.187572591e-02 - 3.017214954e-04j, 1.032548527e-02 - 4.159704807e-04j),
        "(1,-1j)": (1.021210440e-02 - 2.169852168e-04j, 1.200087096e-02 - 5.238960948e-04j),
    },
    # 20261010T170356Z_WF-10-step2-W1.log:6325–6328 (node 118.2437 pF)
    C_IN_F_HUMAN_PF: {
        "(1,+1j)": (1.172678791e-02 + 2.783852939e-06j, 1.021104330e-02 - 1.869819907e-04j),
        "(1,-1j)": (1.010006039e-02 + 1.017991223e-05j, 1.184992763e-02 - 2.164134517e-04j),
    },
    # 20261010T093200Z_WF-10-step1e2-W1a.log:6346–6347 (node 118.30 pF)
    118.30: {
        "(1,+1j)": (1.153823690e-02 + 3.890284628e-04j, 1.006795539e-02 + 1.057884710e-04j),
        "(1,-1j)": (9.959965021e-03 + 3.006225189e-04j, 1.165885509e-02 + 1.735839748e-04j),
    },
}
ROOT_JSON = os.path.join(JSON_DIR, "h18-root.json")


def _node_json(c_pf):
    return os.path.join(JSON_DIR, f"h18-node-{c_pf:.4f}.json")


def test_wf10_step3_refined_kwargs_contract():
    """(K): the refined kwargs are the default hole kwargs plus exactly one key."""
    from tests.validation.test_wf7_f_human_b1_quadrature import _f_human_kwargs

    ref = _f_human_kwargs(PHANTOM_RESOLUTION_M, as_hole=True)
    dflt = _f_human_kwargs(as_hole=True)
    extra = {k: ref[k] for k in ref if k not in dflt}
    same = all(ref[k] == dflt[k] for k in dflt)
    print(
        f"[WF-10 step3a K] PHANTOM_RESOLUTION_M {PHANTOM_RESOLUTION_M}; extra keys {extra}; "
        f"default keys all equal {same}; default has phantom_resolution "
        f"{'phantom_resolution' in dflt} (ASSERTED False)"
    )
    assert "phantom_resolution" not in dflt, "negative control: default dict carries the key"
    assert set(dflt) <= set(ref) and same
    assert extra == {"phantom_resolution": 0.0075}, extra


def build_refined_fixture(comm):
    """Fixture at the refined phantom, (A0′)/(iv), port order, feeds, and a node
    closure.  Never a ``test_*``.

    Lifted from `main` at this item's base (1f96adbf),
    tests/validation/test_wf10_f_human_tuned_b1_sar.py lines 147–489
    (``build_tuned_quadrature``), with the brief's deltas only: (1) the phantom
    at ``PHANTOM_RESOLUTION_M``; (2) (A0′) prints cells and phantom cells, no
    record, no npz cell equality, census (iv) and npz port order asserted;
    (3) no Z50 block; (4) the node is a parameter (``node(c_pf, *,
    keep_fields=False)``), no sense superposition; (5) the returned keys below;
    (6) ``h18-`` JSONs only.  Plus the brief's cost gates after the build.

    Returned keys: ``built``, ``ctx``, ``order``, ``slots``, ``feeds``,
    ``others``, ``az_feeds``, ``az_mirror``, ``wall``, ``z_s``, ``node``,
    ``say``, ``times``, ``cells``, ``phantom_cells``, ``omega``, ``t_window``,
    ``sigma`` (after the first node).
    """
    import ufl
    from dolfinx import fem

    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.mesh.test_birdcage_phantom_resolution import _global_tag_cell_count
    from tests.validation.test_birdcage_f_human_rung import _ring_ports
    from tests.validation.test_birdcage_power_identity import (
        DISCRETE_IDENTITY_RTOL,
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
    from tests.validation.test_port_drive_superposition import _ring_quadrature_slots
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf7_f_human_b1_quadrature import (
        DEGREE,
        FREQUENCY_HZ,
        _build_f_human,
    )
    from tests.validation.test_wf10_f_human_tuned_circuit import (
        QUADRATURE_SLOT_OFFSET,
        _respec,
    )
    from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR, COPPER_SIGMA
    from tests.validation.test_wf12_f_human_hole import _census_checks, _owned_census, _rss_gib

    z0 = float(REFERENCE_IMPEDANCE_OHM)
    t_window = time.perf_counter()
    times = {}

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-10 step3a] {msg}", flush=True)

    assert os.path.exists(
        NPZ_PATH
    ), f"{NPZ_PATH} absent: item 367a stops (BLOCKED, rule (d)); never re-solve"
    with np.load(NPZ_PATH, allow_pickle=False) as npz:
        order_npz = [str(p) for p in npz["order"]]
    say(
        f"start: -n {comm.size}, f = {FREQUENCY_HZ:.3e} Hz, degree {DEGREE}, kappa off, r = "
        f"{ESR_OHM} Ohm on all 32 gaps; phantom_resolution {PHANTOM_RESOLUTION_M} m; nodes "
        f"{REROOT_NODES_PF} pF"
    )

    # ---- fixture, (A0′), (iv) -------------------------------------------------
    t0 = time.perf_counter()
    built = _build_f_human(comm, phantom_resolution=PHANTOM_RESOLUTION_M, as_hole=True)
    times["build"] = time.perf_counter() - t0
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    assert abs(omega / (2.0 * np.pi * float(FREQUENCY_HZ)) - 1.0) < 1e-12, omega
    n_cells_global = int(msh.topology.index_map(msh.topology.dim).size_global)
    phantom_cells = _global_tag_cell_count(msh, cell_tags, PHANTOM_CELL_TAG, comm)
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, built["cells"], _ring_ports())
    _m, rss_per, rss_sum = _rss_gib(comm)
    say(
        f"(A0′) hole cells {built['cells']} (index map {n_cells_global}) beside default "
        f"{DEFAULT_HOLE_CELLS} ({built['cells'] / DEFAULT_HOLE_CELLS - 1.0:+.4%}; PREDICTED ≈ "
        f"{PREDICTED_CELLS}; printed, no record); phantom cells {phantom_cells} beside refined-solid "
        f"{REFINED_SOLID_PHANTOM_CELLS} (PREDICTED ≈ {PREDICTED_PHANTOM_CELLS}; printed); build "
        f"{times['build']:.1f} s (stop {BUILD_STOP_S:g} s); summed RSS {rss_sum:.2f} GiB (stop "
        f"{RSS_STOP_GIB:g} GiB); (iv) failures {census_fails or 'none'} (ASSERTED none)"
    )
    assert not census_fails, census_fails
    assert (
        times["build"] <= BUILD_STOP_S
    ), f"cost gate: build {times['build']:.1f} s > {BUILD_STOP_S} s — park, the review re-prices"
    assert (
        rss_sum <= RSS_STOP_GIB
    ), f"cost gate: summed RSS {rss_sum:.2f} GiB > {RSS_STOP_GIB} GiB — park"

    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}
    solver = ctx["solver"]
    base_specs = list(ctx["specs"])
    keep = {}

    def _vol(sigma_f, e, measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    # `_drive`, verbatim from 358's builder (itself lifted from H17's module).
    def _drive(specs, driven_id, *, held, partition, keep_as=None):
        sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in specs]
        dsheet = next(s for s in sheets if s.port_id == driven_id)
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
        t_solve = time.perf_counter() - t0
        t1 = time.perf_counter()
        e = fields.e_complex
        currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
        v_src = complex(dsheet.source_voltage_v)
        s_col = {}
        for s in sheets:
            v = complex(s.source_voltage_v) - currents[s.port_id] * z_p
            s_col[s.port_id] = (v - z_p * currents[s.port_id]) / v_src
        out = {
            "s": s_col,
            "currents": currents,
            "kind": "held" if held else "solve",
            "t_solve": t_solve,
        }
        if "sigma" not in keep:
            keep["sigma"] = fields.sigma_field
        if partition:
            sigma_f = fields.sigma_field
            by_id = {
                s.port_id: float(_sheet_field_dissipation_w(acct, s, e, omega)) for s in sheets
            }
            out["p_sheet_by_id"] = by_id
            out["p_sheets"] = float(sum(by_id.values()))
            out["p_src"] = float(_source_power_w(acct, dsheet, e, omega))
            p_total = _vol(sigma_f, e, dx)
            out["p_phantom"] = _vol(sigma_f, e, dx(PHANTOM_CELL_TAG))
            out["p_rest"] = p_total - out["p_phantom"]
            out["p_surf"] = _surface_loss_w(msh, built["facets"], z_s, e, comm)
        if keep_as is not None:
            cp = fem.Function(e.function_space, name=keep_as)
            cp.x.array[:] = e.x.array
            cp.x.scatter_forward()
            keep[keep_as] = cp
        del fields, e
        out["t_post"] = time.perf_counter() - t1
        return out

    # ---- order (H17's: [P17, feed2] + rest by number) -------------------------
    az_ref, slots = _ring_quadrature_slots(built["sheets"])
    ref_k, ref_sign = slots[DRIVEN]
    feed2 = next(
        pid
        for pid, (k, sg) in slots.items()
        if sg == ref_sign and k == (ref_k + QUADRATURE_SLOT_OFFSET) % 16
    )
    others = sorted((pid for pid in slots if pid not in (DRIVEN, feed2)), key=lambda p: int(p[1:]))
    order = [DRIVEN, feed2] + others
    feeds = order[:2]
    feeds_idx = [0, 1]
    assert len(order) == 32 and sorted(order) == sorted(sp.port_id for sp in base_specs)
    assert order_npz == order, "npz order differs from this run's order"
    az_by = {f"P{s['ordinal']}": float(s["azimuth_deg"]) for s in built["sheets"]}
    az_feeds = (az_by[feeds[0]], az_by[feeds[1]])
    d_az = (az_feeds[1] - az_feeds[0]) % 360.0
    az_mirror = az_feeds[0] + 45.0
    say(
        f"feeds: {feeds[0]} (slot {ref_k}, ring {ref_sign:+.0f}, az {az_feeds[0]:.9f} deg) and "
        f"{feeds[1]} (slot {slots[feeds[1]][0]}, az {az_feeds[1]:.9f} deg); az21 - az17 = "
        f"{d_az:.9f} deg (90 within {AZIMUTH_MATCH_DEG:g}, ASSERTED); npz port order equal "
        f"(ASSERTED); mirror plane az {az_mirror:.9f} deg; az_ref {az_ref:.9f} deg"
    )
    assert abs(d_az - 90.0) < AZIMUTH_MATCH_DEG, d_az
    assert slots[feeds[0]][0] == 0 and slots[feeds[1]][0] == 4, (slots[feeds[0]], slots[feeds[1]])

    # `_pair`, verbatim from 358's builder.
    def _pair(imp_x, z_feed, v_feed, keep_fields=False):
        specs_x = _respec(base_specs, None, imp_x)
        assert all(sp.width_correction_kappa is None for sp in specs_x)
        cur, knd, tt, parts = [], [], [], []
        for n_drive, pid in enumerate(feeds):
            d = _drive(
                specs_x,
                pid,
                held=(n_drive > 0),
                partition=(n_drive == 0),
                keep_as=(f"E_{pid}" if keep_fields else None),
            )
            knd.append(d["kind"])
            tt.append(d["t_solve"] + d["t_post"])
            cur.append(np.array([d["currents"][p] for p in order], dtype=np.complex128))
            if n_drive == 0:
                parts.append({k: v for k, v in d.items() if k.startswith("p_")})
            del d
        return _current_route_s(cur, feeds_idx, z_feed, v_feed), cur, knd, tt, parts[0], specs_x

    def _partition(p):
        p_vol = p["p_phantom"] + p["p_rest"]
        rel = abs(p["p_src"] - p["p_sheets"] - p_vol - p["p_surf"]) / abs(p["p_src"])
        blind = abs(p["p_src"] - p["p_sheets"] - p_vol) / abs(p["p_src"])
        return rel, blind

    out = {}

    def node(c_pf, *, keep_fields=False):
        """H17's node pair at ``c_pf``: (P_RC) and kinds asserted; ``h18-node``."""
        c_f = float(c_pf) * 1e-12
        z_rc = _z_rc(c_f, FREQUENCY_HZ)
        z_th = _thevenin_feed_impedance(z_rc, z0)
        v_s = _source_voltage(z_rc, z0)
        say(
            f"node C = {c_f:.12e} F ({c_pf:.4f} pF): Z_RC {z_rc:.9e} Ohm, Z_th {z_th:.9e} Ohm, "
            f"V_s {v_s:.9e} V"
        )
        imp = {pid: z_rc for pid in others}
        imp.update({pid: z_th for pid in feeds})
        s_node, cur_node, kinds_n, t_n, pn, specs_n = _pair(imp, z_th, v_s, keep_fields=keep_fields)
        t_pair = float(sum(t_n))
        rel_prc, blind_n = _partition(pn)
        p_rc = sum(v for pid, v in pn["p_sheet_by_id"].items() if pid not in feeds)
        p_vol = pn["p_phantom"] + pn["p_rest"]
        dz = abs(pn["p_src"] - (pn["p_sheets"] - p_rc) - p_vol - pn["p_surf"]) / abs(pn["p_src"])
        say(f"node pair {kinds_n}, {', '.join(f'{t:.2f}' for t in t_n)} s (sum {t_pair:.1f} s)")
        say(
            f"(P_RC) P17 drive: P_src {pn['p_src']:.9e} sumP_sheet {pn['p_sheets']:.9e} (30 RC "
            f"{p_rc:.9e}) P_phantom {pn['p_phantom']:.9e} P_rest {pn['p_rest']:.9e} P_surf "
            f"{pn['p_surf']:.9e} W; rel {rel_prc:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED; "
            f"default mesh 7.236e-13 … 7.946e-13)"
        )
        say(
            f"NEG P_surf-blind miss {blind_n:.4e} = {blind_n / DISCRETE_IDENTITY_RTOL:.1f}x the band "
            f"(bar {BLIND_FACTOR:.0f}x; PREDICTED ≈ {PREDICTED_BLIND:g}, printed); (Dz) "
            f"{dz:.4e} (PREDICTED ≈ {PREDICTED_DZ}, printed)"
        )
        assert kinds_n == ["solve", "held"], kinds_n
        assert (
            rel_prc <= DISCRETE_IDENTITY_RTOL
        ), f"(P_RC) {rel_prc:.3e} at {c_pf} pF — never re-band"
        if "sigma" not in out and "sigma" in keep:
            out["sigma"] = keep["sigma"]
        rec = {
            "c_f": c_f,
            "r_ohm": ESR_OHM,
            "freq_hz": float(FREQUENCY_HZ),
            "order": order,
            "s_node": [[[complex(v).real, complex(v).imag] for v in row] for row in s_node],
            "y_in": {
                lbl: [[y.real, y.imag] for y in _active_y_in(s_node, a, z0)]
                for lbl, a in SENSES.items()
            },
            "partition": {k: v for k, v in pn.items()},
            "rel_p_rc": rel_prc,
            "dz": dz,
            "blind": blind_n,
            "times_s": {"pair": t_n, "build": times["build"]},
            "cells": int(built["cells"]),
            "phantom_cells": int(phantom_cells),
            "phantom_resolution_m": float(PHANTOM_RESOLUTION_M),
        }
        if comm.rank == 0:
            os.makedirs(JSON_DIR, exist_ok=True)
            path = _node_json(c_pf)
            with open(path, "w") as fh:
                json.dump(rec, fh, indent=1)
            print(f"[WF-10 step3a] wrote {path}", flush=True)
        comm.Barrier()
        res = {
            "s_node": s_node,
            "currents": cur_node,
            "partition": pn,
            "rel_p_rc": rel_prc,
            "dz": dz,
            "blind": blind_n,
            "kinds": kinds_n,
            "t_pair": t_pair,
            "specs": specs_n,
        }
        if keep_fields:
            res["e17"] = keep.pop(f"E_{feeds[0]}")
            res["e21"] = keep.pop(f"E_{feeds[1]}")
        return res

    out.update(
        {
            "built": built,
            "ctx": ctx,
            "order": order,
            "slots": slots,
            "feeds": feeds,
            "others": others,
            "az_feeds": az_feeds,
            "az_mirror": az_mirror,
            "wall": wall,
            "z_s": z_s,
            "node": node,
            "say": say,
            "times": times,
            "cells": built["cells"],
            "phantom_cells": phantom_cells,
            "omega": omega,
            "t_window": t_window,
        }
    )
    return out


@complex_only
def test_wf10_step3_reroot():
    if os.environ.get(STEP3_ROOT_ENV, "").strip() != "1":
        pytest.skip(f"{STEP3_ROOT_ENV} unset: heavy -n 8 F-human window, environment-selected")

    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    comm = MPI.COMM_WORLD
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    fx = build_refined_fixture(comm)
    say, feeds, omega = fx["say"], fx["feeds"], fx["omega"]

    y_by_node = {}
    pair_s = []
    for c_pf in REROOT_NODES_PF:
        r = fx["node"](c_pf)
        pair_s.append(r["t_pair"])
        y_by_node[c_pf] = {lbl: _active_y_in(r["s_node"], a, z0) for lbl, a in SENSES.items()}
        for lbl in SENSES:
            for k, name in enumerate(feeds):
                y = y_by_node[c_pf][lbl][k]
                yd = DEFAULT_MESH_Y_IN[c_pf][lbl][k]
                sgn = "+" if y.imag > 0 else "-" if y.imag < 0 else "0"
                say(
                    f"Y_in {c_pf:.4f} pF a = {lbl} {name}: {y.real:+.9e} {y.imag:+.9e}j S (sign "
                    f"Im {sgn}); default mesh {yd.real:+.9e} {yd.imag:+.9e}j S; move Im "
                    f"{y.imag - yd.imag:+.4e} S, Re {y.real / yd.real - 1.0:+.4e} rel (printed)"
                )
        del r
        say(f"t={time.perf_counter() - fx['t_window']:.1f} s")

    c = np.array(REROOT_NODES_PF) * 1e-12
    roots = {}
    for lbl in SENSES:
        for k, name in enumerate(feeds):
            b = np.array([y_by_node[cp][lbl][k].imag for cp in REROOT_NODES_PF])
            slopes = np.diff(b) / np.diff(c)
            say(
                f"{name} a = {lbl}: Im Y_in "
                + ", ".join(f"{x * 1e12:.4f}:{v:+.6e}" for x, v in zip(c, b))
                + " S; interval slopes "
                + ", ".join(f"{s * 1e-12:+.4e}" for s in slopes)
                + " S/pF (= "
                + ", ".join(f"{s / (PREDICTED_SLOPE_FACTOR * omega):.3f}" for s in slopes)
                + " x 16w)"
            )
            found = _sign_change_roots(c, b)
            if not found:
                i_near = int(np.argmin(np.abs(b)))
                say(
                    f"{name} a = {lbl}: NO sign change in the grid; nearest node "
                    f"{c[i_near] * 1e12:.4f} pF (Im {b[i_near]:+.6e} S) — a record"
                )
            for c_a, c_b, rt in found:
                say(
                    f"{name} a = {lbl}: sign change bracket [{c_a * 1e12:.4f}, {c_b * 1e12:.4f}] "
                    f"pF, secant root {rt * 1e12:.6f} pF, resolution {(c_b - c_a) * 1e12:.4f} pF "
                    f"({(c_b - c_a) / rt:.3e} rel)"
                )
            roots[(name, lbl)] = found

    key = (feeds[0], "(1,+1j)")
    all_first = [v[0][2] for v in roots.values() if v]
    spread = (max(all_first) - min(all_first)) if all_first else float("nan")
    rec = {
        "nodes_pf": list(REROOT_NODES_PF),
        "cells": int(fx["cells"]),
        "phantom_cells": int(fx["phantom_cells"]),
        "phantom_resolution_m": float(PHANTOM_RESOLUTION_M),
        "roots_pf": {
            f"{n} {l}": ([[x * 1e12 for x in v[0]]] if v else []) for (n, l), v in roots.items()
        },
        "spread_pf": spread * 1e12,
    }
    if roots[key]:
        c_a, c_b, c_in = roots[key][0]
        move = c_in * 1e12 - C_IN_F_HUMAN_PF
        say(
            f"C_in′ ({key[0]} (1,+1j), item 354 ruling (2)) = {c_in * 1e12:.6f} pF, bracket "
            f"[{c_a * 1e12:.4f}, {c_b * 1e12:.4f}] pF, resolution {(c_b - c_a) * 1e12:.4f} pF "
            f"({(c_b - c_a) / c_in:.3e} rel); four-reading roots "
            + ", ".join(
                f"{n} {l}: " + (f"{v[0][2] * 1e12:.6f}" if v else "none")
                for (n, l), v in roots.items()
            )
            + f" pF; spread {spread * 1e12:.6f} pF beside item 354's {ITEM354_SPREAD_PF} pF; move "
            f"C_in′ - {C_IN_F_HUMAN_PF} = {move:+.6f} pF beside item 354's resolution "
            f"{ITEM354_RESOLUTION_PF} pF (PREDICTED in {PREDICTED_MOVE_PF} pF) — RECORD, never "
            f"asserted"
        )
        rec.update(
            {"c_in_prime_pf": c_in * 1e12, "bracket_pf": [c_a * 1e12, c_b * 1e12], "move_pf": move}
        )
    else:
        say(
            f"C_in′: no sign change on {key} in [{REROOT_NODES_PF[0]}, {REROOT_NODES_PF[-1]}] pF — "
            f"a record (367b BLOCKED until the review re-grids)"
        )
        rec.update({"c_in_prime_pf": None, "bracket_pf": None})
    y_mid = y_by_node[C_IN_F_HUMAN_PF]["(1,+1j)"][0].imag
    say(
        f"Im Y_in(P17, (1,+1j)) at {C_IN_F_HUMAN_PF} pF = {y_mid:+.6e} S (PREDICTED |.| <= "
        f"{PREDICTED_IM_Y_MID:g} S; default +2.78e-6 S; printed)"
    )
    if comm.rank == 0:
        os.makedirs(JSON_DIR, exist_ok=True)
        with open(ROOT_JSON, "w") as fh:
            json.dump(rec, fh, indent=1)
        print(f"[WF-10 step3a] wrote {ROOT_JSON}", flush=True)
    comm.Barrier()

    _m, rss_per, rss_sum = _rss_gib(comm)
    say(
        f"price: build {fx['times']['build']:.1f} s, node pairs "
        f"{', '.join(f'{t:.1f}' for t in pair_s)} s; RSS per rank "
        f"{', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed {rss_sum:.2f} GiB; window "
        f"elapsed {time.perf_counter() - fx['t_window']:.1f} s at -n {comm.size}"
    )
