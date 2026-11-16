"""`WF-12` step 5a (§10 chain H16b, §9 item 415) — the degree-2 copper 32×32
module for `xl` entry 36: all-columns and save knobs, a build-only path and a
numpy save / load contract.

Selected by environment only, never by ``-k``:
  ``FEM_EM_WF12_STEP5``           unset / empty ⇒ the gated test skips;
                                  ``build`` ⇒ build-only return after (A0) / (Deg) /
                                  (Ord), before any form is compiled or assembled;
                                  ``all`` ⇒ the 32-column degree-2 sweep (1
                                  factorising + 31 held drives, `PORT-19`), saved
                                  first, reloaded, then step 2's matrix gates.
  ``FEM_EM_WF12_STEP5_SAVE_NPZ``  absolute, normalised, under
                                  ``/workspace/logs/wf12/``, ending ``-d2.npz``,
                                  never one of the three degree-1 records
                                  (``REFUSED``); required with ``all``, refused
                                  with ``build``.  Entry 36's path:
                                  ``/workspace/logs/wf12/h16b-cu-s32-d2.npz``.

Provenance: the degree swap is step 4's (`test_wf12_f_human_copper_degree2.py:
113–145`, entry 32's module, not edited); the drive and the matrix gates are
step 2's (`test_wf12_f_human_copper_matrix.py:170–220`, `:279–342`, not
edited), copied here; the save contract is `WF-7` H15a's
(`test_wf7_f_human_ring_matrix.py:300–328`, `:904–991`).  No band, tolerance or
target is set or moved: every gate band is imported.

Anchors (asserted): (A0) 369 387 cells within ``CELL_COUNT_BAND``, (iv) none,
outer census as step 4; (Deg) degree 2 and unknowns == 2E + 2F exactly
(2 570 024 on this mesh at `-n 2` and `-n 16`); (Ord) 32 ring ids, ``P17``
first.  Negative controls (asserted, rule (e)): a degree-1 space has E
unknowns, not 2E + 2F; the knob parsers refuse the listed bad values; the
writer refuses to overwrite; a complex64 round trip is not bitwise equal.

Scope: the module, its knobs and save format, landed and built at `-n 2`.
Not claimed: any degree-2 column, the 32×32 or its gates, ``C_mode`` or the
k-share at degree 2, entry 36's command / price / tier, AED.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP5_ENV = "FEM_EM_WF12_STEP5"
SAVE_ENV = "FEM_EM_WF12_STEP5_SAVE_NPZ"
MODES = ("build", "all")
SAVE_DIR = "/workspace/logs/wf12/"
SAVE_SUFFIX = "-d2.npz"
# The three degree-1 records on disk (item 285 (E5)): never written by this module.
REFUSED = frozenset({"h12c-cu-s32.npz", "h12c-cu-s32-a2.npz", "h12b-prime-held-s32.npz"})
ENTRY36_SAVE = "/workspace/logs/wf12/h16b-cu-s32-d2.npz"
H12C_NPZ = "/workspace/logs/wf12/h12c-cu-s32.npz"
DRIVEN = "P17"
DEGREE_2 = 2
FREQ_HZ = 64e6
# (Deg) records: `20261008T093243Z_WF-12-step4a-W1.log:6221` (-n 2) and
# `20261013T070008Z_WF-12-step4-H16.log:6532` (-n 16); printed beside.
UNKNOWNS_RECORD = 2_570_024
EDGES_RECORD = 500_568
FACES_RECORD = 784_444
# Step 4's degree-2 S17,17 (entry 32), printed beside only.
S1717_D2_STEP4_RECORD = 0.402723 + 0.391780j
PREDICTED_RSS_GIB_BUILD = 1.1
PER_COLUMN_KEYS = ("t_solve", "t_post", "held", "p_src", "p_sheets", "p_phantom",
                   "p_rest", "p_surf")
SCALAR_KEYS = ("degree", "cells", "n_dofs", "n_edges", "n_faces", "freq_hz", "ranks",
               "kappa_applied")


# ---- knobs ---------------------------------------------------------------
def _mode_knob():
    raw = os.environ.get(STEP5_ENV, "")
    val = raw.strip()
    if not val:
        return None
    if val not in MODES:
        raise ValueError(f"{STEP5_ENV} must be one of {MODES} or unset, got {raw!r}")
    return val


def _save_knob():
    raw = os.environ.get(SAVE_ENV, "")
    if not raw.strip():
        return None
    if raw != raw.strip() or not os.path.isabs(raw) or os.path.normpath(raw) != raw:
        raise ValueError(f"{SAVE_ENV} must be an absolute normalised path, got {raw!r}")
    if not raw.startswith(SAVE_DIR) or os.path.dirname(raw) + "/" != SAVE_DIR:
        raise ValueError(f"{SAVE_ENV} must live directly under {SAVE_DIR}, got {raw!r}")
    base = os.path.basename(raw)
    if not base.endswith(SAVE_SUFFIX) or base == SAVE_SUFFIX:
        raise ValueError(f"{SAVE_ENV} must end {SAVE_SUFFIX!r}, got {raw!r}")
    if base in REFUSED:
        raise ValueError(f"{SAVE_ENV} names a degree-1 record ({base}); refused")
    return raw


def _knobs():
    mode, path = _mode_knob(), _save_knob()
    if mode == "all" and path is None:
        raise ValueError(f"{STEP5_ENV}=all needs {SAVE_ENV}: an unsaved degree-2 32x32 is refused")
    if mode == "build" and path is not None:
        raise ValueError(f"{STEP5_ENV}=build writes nothing; unset {SAVE_ENV}")
    return mode, path


# ---- save / load (numpy only; the caller decides which rank writes) -------
def _save_s32(path, s, order, per_col, meta):
    """Write the matrix with mode ``"xb"`` (never overwrites).  ``s`` rows =
    receive, columns = driven, both in ``order``; ``per_col`` maps each of
    ``PER_COLUMN_KEYS`` to a length-32 sequence in ``order``; ``meta`` holds
    ``SCALAR_KEYS``."""
    s = np.asarray(s)
    if s.dtype != np.complex128:
        raise TypeError(f"s32 must be complex128, got {s.dtype}")
    n = len(order)
    if s.shape != (n, n):
        raise ValueError(f"s32 shape {s.shape} != ({n}, {n})")
    arrays = {"s32": s, "order": np.array([str(o) for o in order], dtype=str)}
    for key in PER_COLUMN_KEYS:
        arr = np.asarray(per_col[key], dtype=bool if key == "held" else np.float64)
        if arr.shape != (n,):
            raise ValueError(f"{key} shape {arr.shape} != ({n},)")
        arrays[key] = arr
    for key in SCALAR_KEYS:
        arrays[key] = np.asarray(meta[key])
    with open(path, "xb") as fh:
        np.savez(fh, **arrays)


def _load_s32(path):
    with np.load(path, allow_pickle=False) as z:
        out = {k: z[k] for k in z.files}
    out["order"] = [str(x) for x in out["order"]]
    return out


# ---- the build and the sweep ----------------------------------------------
def _run_step5(mode, save_path):
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
    from tests.validation.test_port_birdcage_ring_column import (
        AZIMUTH_MATCH_DEG,
        OPPOSITE_SPREAD_BAND,
        _reciprocity_ratio,
    )
    from tests.validation.test_port_birdcage_ring_matrix import (
        MATRIX_PASSIVITY_CEILING,
        N_AZIMUTH_CLASSES,
        N_SYMMETRY_CLASSES,
        _azimuth_class,
    )
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
    from tests.validation.test_th14_birdcage_copper import _leontovich_term, _surface_loss_w
    from tests.validation.test_wf7_f_human_b1_quadrature import _build_f_human
    from tests.validation.test_wf12_f_human_copper_matrix import (
        BLIND_FACTOR,
        COPPER_SIGMA,
        HOLE_CELL_RECORD,
        P17_P_SRC_RECORD_W,
        P17_P_SURF_RECORD_W,
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
            print(f"[WF-12 step5] {msg}", flush=True)

    def _phase(name, t_s):
        """Step 4's per-phase price print (collective: `_rss_gib` reduces)."""
        _mine, per, total = _rss_gib(comm)
        per_s = ", ".join(f"{r:.2f}" for r in (per or []))
        say(f"phase {name}: t {t_s:.1f} s, RSS summed {total:.2f} GiB, per rank {per_s}")

    say(f"start: -n {comm.size}, degree {DEGREE_2}, mode {mode!r}, save {save_path!r}")

    # ---- stale-file refusal (all ranks raise; step 2a's pattern) ----------
    if mode == "all":
        exists = comm.bcast(os.path.exists(save_path) if comm.rank == 0 else None, root=0)
        if exists:
            raise FileExistsError(f"{save_path} exists; this module never overwrites")

    # ---- fixture, (A0), (iv), degree-2 solver (step 4 :113–145) ----------
    t0 = time.perf_counter()
    built = _build_f_human(comm, as_hole=True)
    ctx = built["ctx"]
    msh, tags_f, omega, cell_tags = ctx["msh"], ctx["tags_f"], ctx["omega"], ctx["cell_tags"]
    census = _owned_census(msh, cell_tags, comm)
    census_fails = _census_checks(census, built["cells"], _ring_ports())
    oc = built["outer_census"]
    cell_rel = built["cells"] / HOLE_CELL_RECORD - 1.0
    problem = ctx["solver"].problem
    del ctx["solver"]
    solver = TimeHarmonicSolver(problem, degree=DEGREE_2)
    v_space = solver.function_space()
    imap = v_space.dofmap.index_map
    n_dofs = int(imap.size_global * v_space.dofmap.index_map_bs)
    n_edges = int(msh.topology.index_map(1).size_global)
    n_faces = int(msh.topology.index_map(2).size_global)
    t_build = time.perf_counter() - t0
    ring_ports = _ring_ports()
    all_ids = [f"P{i}" for i in ring_ports]
    say(f"(A0) hole cells {built['cells']} vs record {HOLE_CELL_RECORD}, rel {cell_rel:.3e} "
        f"(band {CELL_COUNT_BAND}, ASSERTED); build {t_build:.1f} s; (iv) failures "
        f"{census_fails or 'none'}; outer census {oc}")
    say(f"(Deg) solver.degree {solver.degree}; global unknowns {n_dofs} vs 2E+2F = "
        f"2*{n_edges}+2*{n_faces} = {2 * n_edges + 2 * n_faces} (ASSERTED exact); records "
        f"{UNKNOWNS_RECORD} = 2*{EDGES_RECORD}+2*{FACES_RECORD} (-n 2 step4a W1, -n 16 entry 32; "
        f"printed beside)")
    say(f"(Deg) NEG (ASSERTED, rule (e)): degree-1 count E = {n_edges} != {n_dofs}; "
        f"n(2)/E = {n_dofs / n_edges:.6f} (record ceiling 5.134216)")
    say(f"(Ord) {len(all_ids)} ring ids, first {all_ids[0]} (ASSERTED: 32, {DRIVEN} first)")
    _phase("build", t_build)
    say(f"RSS after degree-2 space: PREDICTED ≈{PREDICTED_RSS_GIB_BUILD} GiB summed (printed)")
    assert abs(cell_rel) <= CELL_COUNT_BAND, cell_rel
    assert not census_fails, census_fails
    assert oc["tag_collision"] == 0 and oc["cavity"] > 0, oc
    assert oc["cavity_exterior"] == oc["cavity"] and oc["outer"] + oc["cavity"] == oc["exterior"], oc
    assert solver.degree == DEGREE_2, solver.degree
    assert n_dofs == 2 * n_edges + 2 * n_faces, (n_dofs, n_edges, n_faces)
    assert n_dofs != n_edges, (n_dofs, n_edges)
    assert n_dofs * n_edges == n_edges * (2 * n_edges + 2 * n_faces), (n_dofs, n_edges)
    assert len(all_ids) == 32 and all_ids.index(DRIVEN) == 0, all_ids[:3]
    base = {"cells": built["cells"], "n_dofs": n_dofs, "n_edges": n_edges,
            "n_faces": n_faces, "t_build": t_build}
    if mode == "build":
        say(f"{STEP5_ENV}=build: return after (A0)/(Deg)/(Ord); no form compiled, nothing "
            f"assembled; window {time.perf_counter() - t_window:.1f} s")
        return base

    # ---- the sweep (step 2 :170–220, wall first, field_hook dropped) ------
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, COPPER_SIGMA)
    wall = _leontovich_term(msh, built["facets"], z_s, omega)
    z_p = complex(TERMINATED_PORT_IMPEDANCE_OHM)
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    acct = {"mesh": msh, "facet_tags": tags_f}

    def _drive(driven_id, *, held):
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
        _phase(f"{'held' if held else 'factorising'} solve {driven_id}", t_solve)
        t1 = time.perf_counter()
        e = fields.e_complex
        currents = {s.port_id: sheet_terminal_current(msh, tags_f, s, e, comm) for s in sheets}
        v_src = complex(dsheet.source_voltage_v)
        s_col = {}
        for s in sheets:
            v = complex(s.source_voltage_v) - currents[s.port_id] * z_p
            s_col[s.port_id] = (v - z_p * currents[s.port_id]) / v_src
        sigma_f = fields.sigma_field

        def _vol(measure):
            form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
            return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

        p_sheets = float(sum(_sheet_field_dissipation_w(acct, s, e, omega) for s in sheets))
        p_src = float(_source_power_w(acct, dsheet, e, omega))
        p_total = _vol(dx)
        p_phantom = _vol(dx(PHANTOM_CELL_TAG))
        p_surf = _surface_loss_w(msh, built["facets"], z_s, e, comm)
        del fields
        t_post = time.perf_counter() - t1
        _phase(f"post {driven_id}", t_post)
        return {"s": s_col, "held": held, "t_solve": t_solve, "t_post": t_post,
                "p_src": p_src, "p_sheets": p_sheets, "p_phantom": p_phantom,
                "p_rest": p_total - p_phantom, "p_surf": p_surf}

    cols = {}
    for k, pid in enumerate(all_ids):
        d = _drive(pid, held=(k > 0))
        p_vol = d["p_phantom"] + d["p_rest"]
        d["rel"] = abs(d["p_src"] - d["p_sheets"] - p_vol - d["p_surf"]) / abs(d["p_src"])
        d["blind"] = abs(d["p_src"] - d["p_sheets"] - p_vol) / abs(d["p_src"])
        cols[pid] = d
        say(f"Cu d2 drive {k + 1:2d}/32 {pid}: {'held' if d['held'] else 'factorising'} "
            f"{d['t_solve']:.2f} s + post {d['t_post']:.2f} s; S_driven {d['s'][pid]:.6f}; "
            f"(P) rel {d['rel']:.3e}; NEG blind {d['blind'] / DISCRETE_IDENTITY_RTOL:.1f}x band "
            f"t={time.perf_counter() - t_window:.1f} s")

    s = np.array([[cols[pj]["s"][pi] for pj in all_ids] for pi in all_ids], dtype=np.complex128)
    per_col = {key: [cols[p][key] for p in all_ids] for key in PER_COLUMN_KEYS}
    meta = {"degree": DEGREE_2, "cells": built["cells"], "n_dofs": n_dofs, "n_edges": n_edges,
            "n_faces": n_faces, "freq_hz": FREQ_HZ, "ranks": comm.size, "kappa_applied": False}

    # ---- save first (H15a): rank 0 writes, every rank raises on failure ---
    err = None
    if comm.rank == 0:
        try:
            _save_s32(save_path, s, all_ids, per_col, meta)
        except Exception as exc:  # noqa: BLE001 — broadcast, then raised on every rank
            err = repr(exc)
    err = comm.bcast(err, root=0)
    if err is not None:
        raise RuntimeError(f"save to {save_path} failed: {err}")
    loaded = comm.bcast(_load_s32(save_path) if comm.rank == 0 else None, root=0)
    say(f"saved {save_path}; (L1) reload dtype {loaded['s32'].dtype}, bitwise "
        f"{np.array_equal(loaded['s32'].view(np.uint64), s.view(np.uint64))}, order equal "
        f"{loaded['order'] == all_ids} (ASSERTED)")
    assert loaded["s32"].dtype == np.complex128
    assert np.array_equal(loaded["s32"].view(np.uint64), s.view(np.uint64))
    assert loaded["order"] == all_ids

    # ---- (P), NEG on every column; (R) printed only ----------------------
    first = cols[DRIVEN]
    r_src = abs(first["p_src"] / P17_P_SRC_RECORD_W - 1.0)
    r_surf = abs(first["p_surf"] / P17_P_SURF_RECORD_W - 1.0)
    worst_rel = max(cols[p]["rel"] for p in all_ids)
    min_blind = min(cols[p]["blind"] for p in all_ids)
    say(f"(R) degree-1 comparand, PRINTED only: P17 P_src rel {r_src:.3e}, P_surf rel {r_surf:.3e}")
    say(f"(P) worst rel {worst_rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); NEG min "
        f"blind {min_blind / DISCRETE_IDENTITY_RTOL:.1f}x band (bar {BLIND_FACTOR:.0f}x, ASSERTED)")
    assert worst_rel <= DISCRETE_IDENTITY_RTOL, worst_rel
    assert min_blind >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL, min_blind

    # ---- matrix gates (step 2 :283–339, verbatim bands) -------------------
    ratio = _reciprocity_ratio(s)
    sigma_max = float(np.linalg.norm(s, 2))
    col_norms = np.sum(np.abs(s) ** 2, axis=0)
    by_ord = {sh["ordinal"]: sh for sh in built["sheets"]}
    az = np.array([by_ord[i]["azimuth_deg"] for i in ring_ports])
    zc = np.array([by_ord[i]["z"] for i in ring_ports])
    n = len(ring_ports)
    classes = {}
    worst_residual_deg = 0.0
    for i in range(n):
        for j in range(n):
            steps, residual_deg = _azimuth_class(az[i], az[j])
            worst_residual_deg = max(worst_residual_deg, residual_deg)
            classes.setdefault((bool(zc[i] * zc[j] > 0.0), steps), []).append((i, j))
    spreads = {}
    for key, ij in classes.items():
        m = np.array([abs(s[i, j]) for i, j in ij])
        spreads[key] = (float((m.max() - m.min()) / m.mean()), len(m), float(m.mean()))
    worst = max(spreads, key=lambda key: spreads[key][0])

    # Degree-1 class means from h12c-cu-s32.npz, read BY NAME (two orders trap).
    d1_means = None
    if comm.rank == 0 and os.path.exists(H12C_NPZ):
        try:
            z1 = _load_s32(H12C_NPZ)
            idx = {p: z1["order"].index(p) for p in all_ids}
            d1_means = {key: float(np.mean([abs(z1["s32"][idx[all_ids[i]], idx[all_ids[j]]])
                                            for i, j in ij])) for key, ij in classes.items()}
        except Exception as exc:  # noqa: BLE001 — optional read, printed only
            print(f"[WF-12 step5] h12c read failed: {exc!r} (not a failure)", flush=True)
    d1_means = comm.bcast(d1_means, root=0)
    say(f"(ii) reciprocity {ratio:.6e} (band {RECIPROCITY_BAND:.0e}, ASSERTED)")
    say(f"(iii) sigma_max {sigma_max:.12f}, column sum|S|^2 max {col_norms.max():.6f} (ceiling "
        f"{MATRIX_PASSIVITY_CEILING}, ASSERTED)")
    say(f"(iv') classes {len(classes)} (expected {N_SYMMETRY_CLASSES}, ASSERTED); azimuth "
        f"residual {worst_residual_deg:.3e} deg (< {AZIMUTH_MATCH_DEG:.0e}, ASSERTED)")
    for same_ring in (True, False):
        for steps in range(N_AZIMUTH_CLASSES):
            key = (same_ring, steps)
            if key not in spreads:
                continue
            sp, cnt, mean = spreads[key]
            d1 = f"  d1 {d1_means[key]:.9e}" if d1_means else ""
            say(f"     {'same ' if same_ring else 'other'} ring, {steps} steps  n = {cnt:4d}  "
                f"mean |S| d2 {mean:.9e}{d1}  spread {sp * 100:8.4f}%")
    say(f"     worst class spread {spreads[worst][0]:.4e} (band {OPPOSITE_SPREAD_BAND}, ASSERTED)")
    say(f"S17,17 {cols[DRIVEN]['s'][DRIVEN]:.6f} beside step 4 {S1717_D2_STEP4_RECORD:.6f} (printed)")
    held_cad = [cols[p]["t_solve"] + cols[p]["t_post"] for p in all_ids[1:]]
    say(f"PRICE: build {t_build:.1f} s + factorising {first['t_solve'] + first['t_post']:.1f} s + "
        f"31 held, mean {np.mean(held_cad):.2f} s, max {max(held_cad):.2f} s; window "
        f"{time.perf_counter() - t_window:.1f} s")
    assert ratio <= RECIPROCITY_BAND, ratio
    assert sigma_max <= MATRIX_PASSIVITY_CEILING, sigma_max
    assert float(col_norms.max()) <= MATRIX_PASSIVITY_CEILING, col_norms.max()
    assert len(classes) == N_SYMMETRY_CLASSES, len(classes)
    assert worst_residual_deg < AZIMUTH_MATCH_DEG, worst_residual_deg
    assert spreads[worst][0] <= OPPOSITE_SPREAD_BAND, spreads[worst]
    return {**base, "s": s, "cols": cols}


# ---- tests -----------------------------------------------------------------
def test_wf12_step5a_knob_coding(monkeypatch):
    """(K0a): the parsers; the rejected values are the negative control."""
    monkeypatch.delenv(STEP5_ENV, raising=False)
    monkeypatch.delenv(SAVE_ENV, raising=False)
    assert _mode_knob() is None and _save_knob() is None and _knobs() == (None, None)
    monkeypatch.setenv(STEP5_ENV, "")
    assert _mode_knob() is None
    for ok, want in (("build", "build"), (" all ", "all")):
        monkeypatch.setenv(STEP5_ENV, ok)
        assert _mode_knob() == want
    for bad in ("1", "P17", "32", "ALL", "x"):
        monkeypatch.setenv(STEP5_ENV, bad)
        with pytest.raises(ValueError):
            _mode_knob()
    monkeypatch.setenv(SAVE_ENV, ENTRY36_SAVE)
    assert _save_knob() == ENTRY36_SAVE
    bad_paths = ["logs/wf12/x-d2.npz", "/workspace/logs/wf7/x-d2.npz",
                 "/workspace/logs/wf12/../wf12/x-d2.npz", "/workspace/logs/wf12/x.npz",
                 "/workspace/logs/wf12/x-d2.npy"] + [SAVE_DIR + b for b in sorted(REFUSED)]
    for bad in bad_paths:
        monkeypatch.setenv(SAVE_ENV, bad)
        with pytest.raises(ValueError):
            _save_knob()
    monkeypatch.delenv(SAVE_ENV)
    monkeypatch.setenv(STEP5_ENV, "all")
    with pytest.raises(ValueError):
        _knobs()
    monkeypatch.setenv(STEP5_ENV, "build")
    assert _knobs() == ("build", None)
    monkeypatch.setenv(SAVE_ENV, ENTRY36_SAVE)
    with pytest.raises(ValueError):
        _knobs()
    monkeypatch.setenv(STEP5_ENV, "all")
    assert _knobs() == ("all", ENTRY36_SAVE)


def test_wf12_step5a_save_round_trip(tmp_path):
    """(K0b): seeded synthetic 32x32, bitwise reload; (n1) no overwrite;
    (n2) a complex64 round trip is not bitwise equal."""
    rng = np.random.default_rng(20261013)
    s = rng.standard_normal((32, 32)) + 1j * rng.standard_normal((32, 32))
    assert s.dtype == np.complex128
    order = [f"P{i}" for i in [17] + [k for k in range(1, 33) if k != 17]]
    per_col = {key: (rng.random(32) > 0.5 if key == "held" else rng.random(32))
               for key in PER_COLUMN_KEYS}
    meta = {"degree": 2, "cells": 369_387, "n_dofs": UNKNOWNS_RECORD, "n_edges": EDGES_RECORD,
            "n_faces": FACES_RECORD, "freq_hz": FREQ_HZ, "ranks": 2, "kappa_applied": False}
    path = tmp_path / "synthetic-d2.npz"
    _save_s32(path, s, order, per_col, meta)
    z = _load_s32(path)
    assert z["s32"].dtype == np.complex128
    assert np.array_equal(z["s32"].view(np.uint64), s.view(np.uint64))
    assert z["order"] == order
    for key in PER_COLUMN_KEYS:
        assert np.array_equal(z[key], np.asarray(per_col[key])), key
    for key in SCALAR_KEYS:
        assert z[key].item() == meta[key], key
    # (n1) the writer never overwrites; the file bytes are unchanged.
    before = path.read_bytes()
    with pytest.raises(FileExistsError):
        _save_s32(path, s, order, per_col, meta)
    assert path.read_bytes() == before
    # (n2) complex64 storage would lose bits: not bitwise equal.
    path64 = tmp_path / "synthetic64-d2.npz"
    np.savez(path64, s32=s.astype(np.complex64))
    with np.load(path64, allow_pickle=False) as z64:
        back = z64["s32"].astype(np.complex128)
    assert not np.array_equal(back.view(np.uint64), s.view(np.uint64))
    # and the writer refuses a non-complex128 matrix outright.
    with pytest.raises(TypeError):
        _save_s32(tmp_path / "c64-d2.npz", s.astype(np.complex64), order, per_col, meta)
    print(f"[WF-12 step5] (K0b) synthetic 32x32 bitwise reload True; overwrite refused; "
          f"complex64 round trip max |d| {np.max(np.abs(back - s)):.3e} (not bitwise)", flush=True)


@complex_only
def test_wf12_step5_copper_s32_degree2():
    mode, path = _knobs()
    if mode is None:
        pytest.skip(f"{STEP5_ENV} unset: WF-12 step 5 is selected by environment only")
    _run_step5(mode, path)
