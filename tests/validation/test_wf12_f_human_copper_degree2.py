"""`WF-12` step 4 (§10 chain H16, §9 item 264) — one degree-2 P17 column on
the 369 387-cell copper hole mesh, priced.

Selected by environment only (``FEM_EM_WF12_STEP4=1``), never by ``-k``.
``FEM_EM_WF12_STEP4_BUILD_ONLY=1`` (§9 item 312, step 4a) returns right after
the (A0) / (Deg) asserts, before any form is compiled or assembled; per-phase
``phase <name>: t … s, RSS …`` prints price the build, solve and post.

Fixture: ``_build_f_human(comm, as_hole=True)`` unchanged (its ``DEGREE = 1``
is not edited or monkeypatched); a fresh ``TimeHarmonicSolver(problem,
degree=2)`` is built on the same problem object (mesh / tags / outer-box PEC
group, cavity not pinned) after the degree-1 solver is deleted.  One
factorising Cu drive of P17 with `TH-14`'s Leontovich copper wall first in the
bilinear list and all 32 ring sheets at 50 Ω, κ off.

Anchors (asserted, every band imported):
  (A0)  hole cells = ``HOLE_CELL_RECORD`` within ``CELL_COUNT_BAND``; (iv) reduced;
  (Deg) solver degree == 2 and global unknowns == 2·E + 2·F exactly (the
        degree-2 Nédélec first-kind count on tets);
  (P)   P_src − ΣP_sheet = P_phantom + P_(Ω∖phantom) + P_surf within
        ``DISCRETE_IDENTITY_RTOL``, and P_surf > 0;
  (Col) column Σ|S_i,17|² ≤ ``COLUMN_PASSIVITY_CEILING``.
Negative control (asserted, rule (e)): the P_surf-blind identity misses by
≥ ``BLIND_FACTOR`` × ``DISCRETE_IDENTITY_RTOL`` (degree 1 on this mesh:
123.3× / 123.6×, `20261002T094751Z_WF-12-step2-W2.log:6353`).
Predicted, printed only: |S17,17(d2) − S17,17(d1)|/|S17,17(d1)| in
[0.5 %, 8 %]; summed RSS 60–90 GiB; build + degree-2 drive 350–550 s; the
coil-loss share beside degree 1's 0.3704; the projected `xl` entry 32 price.

Scope: one column, degree 2, 64 MHz, κ off, this copper fixture.  Not a
degree-2 32×32, not C_mode at degree 2, not an h-statement, not AED.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP4_ENV = "FEM_EM_WF12_STEP4"
# (K) §9 item 312: build-only knob; unset in `xl` entry 32's command (ruling 312-7).
BUILD_ONLY_ENV = "FEM_EM_WF12_STEP4_BUILD_ONLY"
DRIVEN = "P17"
DEGREE_2 = 2
H12C_NPZ = "/workspace/logs/wf12/h12c-cu-s32.npz"

# Degree-1 record on this mesh (`20261002T094751Z_WF-12-step2-W2.log:6411`).
S1717_D1_RECORD = 0.383271 + 0.391623j
D1_COIL_SHARE = 0.3704
# Predictions written first (§9 item 264), printed only.
PREDICTED_D12_MOVE = (0.005, 0.08)
PREDICTED_RSS_GIB = (60.0, 90.0)
PREDICTED_BUILD_DRIVE_S = (350.0, 550.0)
PREDICTED_UNKNOWNS = 2.4e6
HELD_PROBE_CUTOFF_S = 480.0


def _run_step4():
    import ufl
    from dolfinx import fem

    from fem_em_solver.core import TimeHarmonicSolver
    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.lumped import (
        lumped_port_bilinear_term,
        lumped_port_linear_term,
        sheet_terminal_current,
    )
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_birdcage_f_human_rung import _ring_ports
    from tests.validation.test_birdcage_power_identity import (
        DISCRETE_IDENTITY_RTOL,
        _sheet_field_dissipation_w,
        _source_power_w,
    )
    from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
    from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
    from tests.validation.test_port_birdcage_ring_column import COLUMN_PASSIVITY_CEILING
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf7_f_human_b1_quadrature import _build_f_human
    from tests.validation.test_wf12_f_human_copper_matrix import (
        BLIND_FACTOR,
        COPPER_SIGMA,
        HOLE_CELL_RECORD,
    )
    from tests.validation.test_wf12_f_human_hole import (
        _census_checks,
        _owned_census,
        _rss_gib,
    )

    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()

    def say(msg):
        if comm.rank == 0:
            print(f"[WF-12 step4] {msg}", flush=True)

    def _phase(name, t_s):
        """(T) §9 item 312: per-phase price print (collective: `_rss_gib` reduces)."""
        _mine, per, total = _rss_gib(comm)
        per_s = ", ".join(f"{r:.2f}" for r in (per or []))
        say(f"phase {name}: t {t_s:.1f} s, RSS summed {total:.2f} GiB, per rank {per_s}")

    say(f"start: -n {comm.size}, degree {DEGREE_2}, one factorising Cu drive of {DRIVEN}")

    # ---- fixture, (A0), (iv) ---------------------------------------------
    t0 = time.perf_counter()
    built = _build_f_human(comm, as_hole=True)
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, built["cells"], _ring_ports())
    oc = built["outer_census"]
    cell_rel = built["cells"] / HOLE_CELL_RECORD - 1.0
    # Degree-2 solver on the same problem; degree-1 solver deleted unassembled.
    problem = ctx["solver"].problem
    del ctx["solver"]
    solver = TimeHarmonicSolver(problem, degree=DEGREE_2)
    v_space = solver.function_space()
    imap = v_space.dofmap.index_map
    n_dofs = int(imap.size_global * v_space.dofmap.index_map_bs)
    n_edges = int(msh.topology.index_map(1).size_global)
    n_faces = int(msh.topology.index_map(2).size_global)
    t_build = time.perf_counter() - t0
    say(f"(A0) hole cells {built['cells']} vs record {HOLE_CELL_RECORD}, rel {cell_rel:.3e} "
        f"(band {CELL_COUNT_BAND}, ASSERTED); build {t_build:.1f} s; (iv) failures "
        f"{census_fails or 'none'}; outer census {oc}")
    say(f"(Deg) solver.degree {solver.degree}; global unknowns {n_dofs} vs 2E+2F = "
        f"2*{n_edges}+2*{n_faces} = {2 * n_edges + 2 * n_faces} (ASSERTED exact; "
        f"PREDICTED ≈{PREDICTED_UNKNOWNS:.1e})")
    say(f"(Deg) control (PREDICTED, printed; rule (e)): a degree-1 space would read E = "
        f"{n_edges} unknowns, not 2E+2F (the `WF-15` step 1 degree-1 default trap)")
    _phase("build", t_build)
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    assert oc["tag_collision"] == 0 and oc["cavity"] > 0, oc
    assert oc["cavity_exterior"] == oc["cavity"] and oc["outer"] + oc["cavity"] == oc["exterior"], oc
    assert solver.degree == DEGREE_2, solver.degree
    assert n_dofs == 2 * n_edges + 2 * n_faces, (n_dofs, n_edges, n_faces)
    if os.environ.get(BUILD_ONLY_ENV, "").strip() == "1":
        # (K) §9 item 312: build-only knob — no form compiled, nothing assembled.
        say(f"{BUILD_ONLY_ENV}=1: build-only return after (A0)/(Deg); no assembly, "
            f"no factorisation; window {time.perf_counter() - t_window:.1f} s")
        return {"cells": built["cells"], "n_dofs": n_dofs, "t_build": t_build}

    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}

    def _drive(driven_id, *, held, partition):
        """Local copy of H12b's `_drive` (`test_wf12_f_human_copper_matrix.py:
        170–217`) with the wall first in the factorising call; ``partition``
        switches the power accounting off for the held cadence probe."""
        sheets = [sp.sheet(driven=(sp.port_id == driven_id)) for sp in ctx["specs"]]
        dsheet = next(s for s in sheets if s.port_id == driven_id)
        linear = [lambda test, _s=dsheet: lumped_port_linear_term(
            msh, tags_f, _s, test, omega_rad_per_s=omega)]
        comm.Barrier()
        t0 = time.perf_counter()
        if held:
            fields = solver.solve_with_held_factorization(extra_linear_terms=linear)
        else:
            bil = [wall] + [
                (lambda trial, test, _s=s: lumped_port_bilinear_term(
                    msh, tags_f, _s, trial, test, omega_rad_per_s=omega))
                for s in sheets
            ]
            fields = solver.solve(current_density=None, project_source=False,
                                  extra_bilinear_terms=bil, extra_linear_terms=linear)
        comm.Barrier()
        t_solve = time.perf_counter() - t0
        # Flushed before post-processing so a window killed in post keeps the price.
        _phase("held solve" if held else "factorising solve", t_solve)
        t1 = time.perf_counter()
        e = fields.e_complex
        currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
        v_src = complex(dsheet.source_voltage_v)
        s_col = {}
        for s in sheets:
            v = complex(s.source_voltage_v) - currents[s.port_id] * z_p
            s_col[s.port_id] = (v - z_p * currents[s.port_id]) / v_src
        out = {"s": s_col, "t_solve": t_solve}
        if partition:
            sigma_f = fields.sigma_field

            def _vol(measure):
                form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
                return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

            out["p_sheets"] = float(sum(_sheet_field_dissipation_w(acct, s, e, omega)
                                        for s in sheets))
            out["p_src"] = float(_source_power_w(acct, dsheet, e, omega))
            p_total = _vol(dx)
            out["p_phantom"] = _vol(dx(PHANTOM_CELL_TAG))
            out["p_rest"] = p_total - out["p_phantom"]
            out["p_surf"] = _surface_loss_w(msh, built["facets"], z_s, e, comm)
        del fields
        out["t_post"] = time.perf_counter() - t1
        _phase("held post" if held else "post", out["t_post"])
        return out

    # ---- the degree-2 factorising P17 drive ------------------------------
    d = _drive(DRIVEN, held=False, partition=True)
    t_first = d["t_solve"] + d["t_post"]
    _m, rss_per, rss_sum = _rss_gib(comm)
    p_vol = d["p_phantom"] + d["p_rest"]
    rel = abs(d["p_src"] - d["p_sheets"] - p_vol - d["p_surf"]) / abs(d["p_src"])
    blind = abs(d["p_src"] - d["p_sheets"] - p_vol) / abs(d["p_src"])
    share = d["p_surf"] / (d["p_src"] - d["p_sheets"])
    colsum = float(sum(abs(v) ** 2 for v in d["s"].values()))
    s1717 = d["s"][DRIVEN]
    move = abs(s1717 - S1717_D1_RECORD) / abs(S1717_D1_RECORD)
    say(f"Cu d2 drive {DRIVEN}: factorising solve {d['t_solve']:.2f} s + post {d['t_post']:.2f} s "
        f"t={time.perf_counter() - t_window:.1f} s")
    say(f"(P) P_src {d['p_src']:.9e} W sumP_sheet {d['p_sheets']:.9e} P_phantom "
        f"{d['p_phantom']:.9e} P_rest {d['p_rest']:.9e} P_surf {d['p_surf']:.9e} W; rel "
        f"{rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); P_surf > 0 (ASSERTED)")
    say(f"NEG blind {blind:.4e} = {blind / DISCRETE_IDENTITY_RTOL:.1f}x band (bar "
        f"{BLIND_FACTOR:.0f}x, ASSERTED; d1 123.3x / 123.6x)")
    say(f"(Col) column sum|S_i,17|^2 {colsum:.9f} (ceiling {COLUMN_PASSIVITY_CEILING}, ASSERTED)")
    say(f"coil share P_surf/(P_src-sumP_sheet) {share:.4f} (d1 {D1_COIL_SHARE})")
    say(f"S17,17 d2 {s1717:.6f} vs d1 {S1717_D1_RECORD:.6f}: rel move {move * 100:.3f}% "
        f"(PREDICTED in [{PREDICTED_D12_MOVE[0] * 100:.1f}%, {PREDICTED_D12_MOVE[1] * 100:.0f}%], "
        f"printed)")

    # P17-adjacent same-ring entries, by azimuth.
    by_ord = {sh["ordinal"]: sh for sh in built["sheets"]}
    me = by_ord[int(DRIVEN[1:])]
    same = [o for o, sh in by_ord.items() if o != int(DRIVEN[1:]) and sh["z"] * me["z"] > 0.0]

    def _daz(o):
        a = (by_ord[o]["azimuth_deg"] - me["azimuth_deg"]) % 360.0
        return min(a, 360.0 - a)

    adj = sorted(same, key=_daz)[:2]
    d1_col = None
    if comm.rank == 0 and os.path.exists(H12C_NPZ):
        try:
            z = np.load(H12C_NPZ, allow_pickle=True)
            order = [str(x) for x in z["order"]]
            j = order.index(DRIVEN)
            d1_col = {order[i]: complex(z["s32"][i, j]) for i in range(len(order))}
        except Exception as exc:  # noqa: BLE001 — optional read, printed only
            print(f"[WF-12 step4] npz read failed: {exc!r} (not a failure)", flush=True)
    d1_col = comm.bcast(d1_col, root=0)
    if d1_col is not None:
        say(f"d1 column from {H12C_NPZ}: S17,17 {d1_col[DRIVEN]:.9f}")
    for o in adj:
        pid = f"P{o}"
        d2v = d["s"][pid]
        line = f"adjacent {pid} (dAz {_daz(o):.2f} deg): S d2 {d2v:.6f}"
        if d1_col is not None and pid in d1_col:
            d1v = d1_col[pid]
            line += f" vs d1 {d1v:.6f}, rel move {abs(d2v - d1v) / abs(d1v) * 100:.3f}%"
        say(line)
    if d1_col is not None:
        worst_pid = max(d1_col, key=lambda p: abs(d["s"][p] - d1_col[p]) / abs(d1_col[p]))
        say(f"largest relative d1->d2 entry move: {worst_pid} "
            f"{abs(d['s'][worst_pid] - d1_col[worst_pid]) / abs(d1_col[worst_pid]) * 100:.3f}%")

    say(f"PRICE: build {t_build:.1f} s + d2 factorising drive {t_first:.1f} s = "
        f"{t_build + t_first:.1f} s (PREDICTED {PREDICTED_BUILD_DRIVE_S}); unknowns {n_dofs}; "
        f"RSS per rank {', '.join(f'{r:.2f}' for r in (rss_per or []))} GiB, summed "
        f"{rss_sum:.2f} GiB (PREDICTED {PREDICTED_RSS_GIB})")

    # ---- optional held cadence probe -------------------------------------
    elapsed = time.perf_counter() - t_window
    if elapsed <= HELD_PROBE_CUTOFF_S:
        h = _drive(DRIVEN, held=True, partition=False)
        cad = h["t_solve"] + h["t_post"]
        dh = abs(h["s"][DRIVEN] - s1717) / abs(s1717)
        proj = t_build + t_first + 31 * cad
        say(f"held re-solve {DRIVEN}: {h['t_solve']:.2f} s + post {h['t_post']:.2f} s; "
            f"S17,17 held vs factorised rel {dh:.3e} (printed)")
        say(f"PROJECTED xl entry 32 (-n {comm.size}): build {t_build:.1f} + factor {t_first:.1f} "
            f"+ 31 x {cad:.2f} = {proj:.0f} s")
    else:
        say(f"held probe SKIPPED (elapsed {elapsed:.1f} s > {HELD_PROBE_CUTOFF_S:.0f} s); "
            f"PROJECTED xl entry 32 lower bound: build {t_build:.1f} + factor {t_first:.1f} "
            f"+ 31 x (held cadence, unmeasured)")

    _m, rss_per, rss_sum = _rss_gib(comm)
    say(f"window {time.perf_counter() - t_window:.1f} s; RSS summed {rss_sum:.2f} GiB, per rank "
        f"max {max(rss_per) if rss_per else float('nan'):.2f} GiB")

    assert d["p_surf"] > 0.0, d["p_surf"]
    assert rel <= DISCRETE_IDENTITY_RTOL, rel
    assert blind >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, blind
    assert colsum <= COLUMN_PASSIVITY_CEILING, colsum
    return {"rel": rel, "blind": blind, "colsum": colsum, "s1717": s1717, "move": move}


@complex_only
def test_wf12_step4_copper_degree2_column():
    if os.environ.get(STEP4_ENV, "").strip() != "1":
        pytest.skip(f"{STEP4_ENV} unset: WF-12 step 4 is selected by environment only")
    _run_step4()
