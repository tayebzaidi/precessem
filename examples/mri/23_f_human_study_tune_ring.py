"""Example (`EX-118`): ``Study.tune`` on the 32-port copper F-human -- one solve, the ring
read returns the H12c ``C_mode``, k = 2 misses, and the field of the k = 1 ring-antisymmetric
drive.

`WF-16` step 12 let ``Study.tune`` read the 32-port ring-gap matrix of the copper
(Leontovich) F-human hole birdcage: ``Z = s_to_z(S)``, the k-block of ``W Z W^H``, the
ring-antisymmetric eigenvalue ``z_k``, ``C = 1/(omega Im z_k)``.  No example called it.  This
one runs the gate module's own Study (degree 1, 64 MHz, kappa off, ``-n 8``) through
``study.solve`` -> ``study.tune``, reproduces the gate's anchors, then re-sweeps the kept
problem once with ``keep_fields=True`` (this Study keeps no field) and superposes the 32 kept
drives into the k = 1 and k = 2 ring-antisymmetric patterns.

Every record, band and helper is imported from the gate module (the ANS-1 rule); record digits
sit only in the example-local ``REC`` dict, each with its ``log:line``.  New values print at
<= 6 significant digits.

Asserted: (K) refusals, (W) the drive pattern, (R0), (I), (a), (S=), (P).
Negative controls (asserted): (c1) k = 2 misses the record; (cW) the k = 2 weights fail the
k = 1 target.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:23 -n 8 -t 590

Not claimed: degree 2; an h-statement; a tuned f0, Q, R_in, match, C_in or resonance; the
field of the coil tuned with ``C_mode`` on its gaps (the output is the matched-port k-pattern);
ESR; the solid sigma = 800 F-human; an external check.
"""

from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.circuit import s_to_z  # noqa: E402
from fem_em_solver.ports.ring_modes import (  # noqa: E402
    _canonical_permutation,
    _dft_matrix,
    antisymmetric_eigenpair,
    k_blocks,
    ring_port_map_from_names,
)
from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep  # noqa: E402
from fem_em_solver.ports.superposition import superpose_drives  # noqa: E402
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.faraday import project_to_cg1  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.instrumentation import report_peak_rss  # noqa: E402

import tests.package.test_wf16_study_tune_f_human_copper as g  # noqa: E402
import tests.validation.test_port_birdcage_ring_matrix as ringm  # noqa: E402
import tests.validation.test_port_lumped_rlc_termination as rlc  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_23_f_human_study_tune_ring"
N_SLOTS = 16

# Gate-log records (docs/testing/logs/20261015T213534Z_WF-16-s12-W1.log), asserted at the
# 2026-10-14 10:30 ruling's form: relative 10^-(d-2), d the printed significant digits.
REC = {
    # :7375  C_mode, 13 digits -> 1e-11 (same value as g.WF12_STEP3_C_MODE_CU_F)
    "c_mode_f": 1.181185799070e-10,
    # :7380  k = 2 C, 13 digits -> 1e-11
    "c2_f": 7.881746394502e-11,
    # :7380  (c1) gap, 7 digits -> 1e-5
    "gap_c1": -3.327259e-01,
    # :7346  (R0) rel and :7375 (a) gap are round-off readings: printed beside, band-only
    "r0_rel": 1.632e-12,
    "gap_a": -3.544942e-13,
    # :7361-7364 z_k (4 decimals) and the off-block (7 digits)
    "off_block": 4.195925e-03,
}


def _rel(a, b):
    return abs(a / b - 1.0)


def _canon_weights(k, sign):
    j = np.arange(N_SLOTS)
    bottom = np.exp(2j * np.pi * j * k / N_SLOTS) / (4.0 * np.sqrt(2.0))
    return np.concatenate([bottom, sign * bottom])


def _to_port_order(a_canonical, perm):
    """Canonical-order weights -> the sweep's port order (``perm[c]`` = port index of c)."""
    out = np.zeros(2 * N_SLOTS, dtype=np.complex128)
    for c in range(2 * N_SLOTS):
        out[int(perm[c])] = a_canonical[c]
    return out


def _centre(fn, comm):
    vals, valid = evaluate_vector_field_parallel(fn, np.zeros((1, 3)))
    vals = np.asarray(vals).reshape(-1, 3)
    assert bool(np.all(np.asarray(valid, dtype=bool))), "origin not found"
    v = vals[0]
    return {"b1p": abs(v[0] + 1j * v[1]) / 2.0, "b1m": abs(v[0] - 1j * v[1]) / 2.0,
            "norm": float(np.linalg.norm(v))}


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    def p0(*args):
        if rank0:
            print(*args, flush=True)

    # ---- structural asserts -------------------------------------------------
    assert comm.size == g.RECORD_RANK_WIDTH, (comm.size, g.RECORD_RANK_WIDTH)
    if not np.issubdtype(default_scalar_type, np.complexfloating):
        raise RuntimeError("complex build required (source dolfinx-complex-mode)")
    # two constants carry the name REUSE_REPRODUCTION_RTOL; the 1e-6 one is imported here
    assert g.REUSE_REPRODUCTION_RTOL is ringm.REUSE_REPRODUCTION_RTOL
    assert g.REDUCTION_BAND is rlc.REDUCTION_BAND
    study = g._study()
    assert study.degree == 1
    assert isinstance(study.conductor, g.study_mod.LeontovichWall)
    assert study.conductor.sigma_s_per_m == g.COPPER_SIGMA
    assert study.geometry.as_hole is True
    assert study.result is None

    # ---- (K) refusals, no mesh, before the solve ----------------------------
    for bad in (0, 3):
        try:
            study.tune(k=bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"(K) tune(k={bad}) did not raise ValueError")
    try:
        study.tune(k=True)
    except TypeError:
        pass
    else:
        raise AssertionError("(K) tune(k=True) did not raise TypeError")
    assert study.result is None
    p0("[(K)] tune(k=0), tune(k=3) raise ValueError; tune(k=True) raises TypeError; "
       "study.result is None")

    npz = None
    if rank0:
        import os
        assert os.path.exists(g.NPZ_PATH), f"{g.NPZ_PATH} absent: stop (BLOCKED, rule (d))"
        with np.load(g.NPZ_PATH, allow_pickle=False) as fh:
            npz = {"s32": np.array(fh["s32"], dtype=np.complex128),
                   "order": [str(x) for x in fh["order"]],
                   "cells": int(fh["cells"]), "c_mode": float(fh["c_mode"])}
    npz = comm.bcast(npz, root=0)

    # ---- (1) anchor: the gate's Study through solve -> tune ------------------
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_solve = time.perf_counter() - t0

    # EX-57 setup figure, right after the mesh is built (no-op unless FEM_EM_SETUP_FIGURES=1).
    local_tags = np.unique(res.cell_tags.values)
    all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist()) for t in part})
    hide = tuple(t for t in all_tags if t not in (PHANTOM_CELL_TAG, CONDUCTOR_CELL_TAG))
    names = {PHANTOM_CELL_TAG: "phantom", CONDUCTOR_CELL_TAG: "conductor (hole route: void)"}
    for t in hide:
        names[t] = "air"
    write_setup_figure(
        res.mesh, res.cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=hide, translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title=("mri:23 -- Study.tune on the 32-port copper F-human hole birdcage "
               "(Leontovich, degree 1), phantom (WF-16 step 12)"),
        comm=comm,
    )

    # every tune read happens before the field-keeping sweep; study.result is never assigned
    t1 = study.tune(k=1, comm=comm)
    t2 = study.tune(k=2, comm=comm)
    t1b = study.tune(k=1, comm=comm)
    assert study.result is res

    names_port = [pd.port_id for pd in res.port_definitions]
    s = np.asarray(res.s_matrix, dtype=np.complex128)
    z0 = float(study.excitation.z0_ohm)
    omega = study.omega
    perm = _canonical_permutation(ring_port_map_from_names(names_port))
    blocks, off_all = k_blocks(s_to_z(s, z0)[np.ix_(perm, perm)])
    anti = {m: antisymmetric_eigenpair(blocks[m])[1] for m in range(9)}

    # (R0)
    idx = {n: i for i, n in enumerate(names_port)}
    assert sorted(names_port) == sorted(npz["order"])
    pidx = np.array([idx[n] for n in npz["order"]], dtype=int)
    r0 = g._residual(s[np.ix_(pidx, pidx)], npz["s32"])
    # (I)
    rows = g._rows(s, names_port, t1.c_mode_farad)
    r_i1 = _rel(t1.c_mode_farad, rows["c_mode_f"])
    r_i2 = _rel(t2.c_mode_farad, rows["c_k2_f"])
    r_io = _rel(t1.off_diagonal_residual, rows["off_block"])
    r_iz = abs(t1.c_mode_farad * 2.0 * np.pi * t1.frequency_hz * t1.z_k[1].imag - 1.0)
    # (a), (c1)
    gap_a = t1.c_mode_farad / g.WF12_STEP3_C_MODE_CU_F - 1.0
    gap_c1 = t2.c_mode_farad / g.WF12_STEP3_C_MODE_CU_F - 1.0

    # ---- (W) the drive is the k-th ring-antisymmetric pattern ----------------
    dft = _dft_matrix()
    target = {}
    for k in (1, 2):
        tgt = np.zeros(2 * N_SLOTS, dtype=np.complex128)
        tgt[k] = 1.0 / np.sqrt(2.0)
        tgt[N_SLOTS + k] = -1.0 / np.sqrt(2.0)
        target[k] = tgt
    a_can = {k: _canon_weights(k, -1.0) for k in (1, 2)}
    w_gap = {k: float(np.max(np.abs(dft @ a_can[k] - target[k]), initial=0.0)) for k in (1, 2)}
    cw_gap = float(np.max(np.abs(dft @ a_can[2] - target[1]), initial=0.0))
    a_port = {k: _to_port_order(a_can[k], perm) for k in (1, 2)}
    # the map back: canonical slot c of the port-order vector is a_port[perm[c]]
    for k in (1, 2):
        assert np.array_equal(a_port[k][perm], a_can[k])

    # ---- the field-keeping re-sweep on the kept problem ----------------------
    gc.collect()
    comm.Barrier()
    t0 = time.perf_counter()
    fres = run_n_port_sparameter_sweep(
        res.problem, res.port_definitions,
        lumped_sheet_ports=res.specs, lumped_sheet_facet_tags=res.facet_tags,
        degree=study.degree, keep_fields=True,
        **({"extra_bilinear_terms": res.conductor_terms} if res.conductor_terms else {}),
    )
    comm.Barrier()
    t_fsweep = time.perf_counter() - t0
    rss_sum = report_peak_rss(comm, label="[EX-118] summed peak ru_maxrss after the field sweep")
    assert list(fres.port_ids) == names_port
    s_tie = g._residual(np.asarray(fres.s_matrix, dtype=np.complex128), s)
    hex_same = float(t1b.c_mode_farad).hex() == float(t1.c_mode_farad).hex()
    assert study.result is res

    # ---- (P) superposed drives and passivity --------------------------------
    t0 = time.perf_counter()
    d1 = superpose_drives(fres, a_port[1])
    d2 = superpose_drives(fres, a_port[2])
    tol_p = 1.0 + g.EXACT_IDENTITY_RTOL
    ok_p = all(0.0 <= d.accepted_power_w <= d.available_power_w * tol_p for d in (d1, d2))
    b1 = project_to_cg1(d1.b_complex, name="B_mode1_cg1")
    b2 = project_to_cg1(d2.b_complex, name="B_mode2_cg1")
    t_proj = time.perf_counter() - t0
    space = b1.function_space
    msh = space.mesh
    c1c, c2c = _centre(b1, comm), _centre(b2, comm)

    space_s = fem.functionspace(msh, ("Lagrange", 1))

    def _scalar(name, values):
        f = fem.Function(space_s, name=name)
        assert values.shape[0] == f.x.array.shape[0], (values.shape, f.x.array.shape)
        f.x.array[:] = values
        f.x.scatter_forward()
        return f

    nb1 = b1.x.array.reshape(-1, 3)
    nb2 = b2.x.array.reshape(-1, 3)
    b1_real = fem.Function(space, name="B_mode1_real")
    b1_real.x.array[:] = np.real(b1.x.array)
    b1_real.x.scatter_forward()
    fields = {
        "B_mode1_real": b1_real,
        "B_mode1_abs": _scalar("B_mode1_abs", np.sqrt(np.sum(np.abs(nb1) ** 2, axis=1))),
        "B1_plus_mode1": _scalar("B1_plus_mode1", np.abs(nb1[:, 0] + 1j * nb1[:, 1]) / 2.0),
        "B1_minus_mode1": _scalar("B1_minus_mode1", np.abs(nb1[:, 0] - 1j * nb1[:, 1]) / 2.0),
        "B_mode2_abs": _scalar("B_mode2_abs", np.sqrt(np.sum(np.abs(nb2) ** 2, axis=1))),
    }
    OUTPUT_DIR.mkdir(exist_ok=True)
    t0 = time.perf_counter()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, res.cell_tags, fields, comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    t_xdmf = time.perf_counter() - t0

    # ---- prints (rank 0; no collective inside) ------------------------------
    if rank0:
        print("=" * 72)
        print("EX-118 -- Study.tune on the 32-port copper F-human: ring read, k = 2 miss, field")
        print("=" * 72)
        print(f"[study] degree 1, -n {comm.size}, f = {study.frequency_hz:.6e} Hz, "
              f"z0 = {z0:g} ohm, ports {names_port[0]}..{names_port[-1]} ({len(names_port)})")
        print(f"[(R0)] cells {res.cells} vs record {g.HOLE_CELL_RECORD} (==) vs npz "
              f"{npz['cells']} (==); ||S_study - S_npz||/||S_npz|| = {r0:.3e} (band "
              f"REUSE_REPRODUCTION_RTOL {g.REUSE_REPRODUCTION_RTOL:.1e}; prior {REC['r0_rel']:.3e}, "
              f"log :7346)")
        print(f"[(I)] rel vs ring_mode_rows: C1 {r_i1:.3e}, C2 {r_i2:.3e}, off {r_io:.3e}; "
              f"C1 w Im z1 - 1 = {r_iz:.3e} (band EXACT_IDENTITY_RTOL {g.EXACT_IDENTITY_RTOL:.1e})")
        print(f"[(a)] C_mode = {t1.c_mode_farad:.12e} F; record {g.WF12_STEP3_C_MODE_CU_F:.12e} F; "
              f"gap {gap_a:+.6e} (band REDUCTION_BAND {g.REDUCTION_BAND:.1e}; prior "
              f"{REC['gap_a']:+.6e}, log :7375); vs C_MODE_CU_PF "
              f"{t1.c_mode_farad / (g.C_MODE_CU_PF * 1e-12) - 1:+.3e}; vs npz c_mode "
              f"{t1.c_mode_farad / npz['c_mode'] - 1:+.3e}")
        print(f"[(a)] vs C_in F-human {t1.c_mode_farad / (g.C_IN_F_HUMAN_PF * 1e-12) - 1:+.4e} "
              f"(a different quantity, the r = 0.2 Ohm in-model root, not a tuning error)")
        print(f"[(c1)] k = 2 C = {t2.c_mode_farad:.12e} F, gap {gap_c1:+.6e} "
              f"({abs(gap_c1) / g.REDUCTION_BAND:.1f}x band; gate {REC['gap_c1']:+.6e}, log "
              f":7380); vs WF12_STEP3_C2_CU_F {_rel(t2.c_mode_farad, g.WF12_STEP3_C2_CU_F):.3e}")
        print(f"[(K)] refusals passed before the solve (see above)")
        print(f"[(W)] max |F a_k - target_k|: k = 1 {w_gap[1]:.3e}, k = 2 {w_gap[2]:.3e} (band "
              f"{g.EXACT_IDENTITY_RTOL:.1e}); (cW) k = 2 weights vs the k = 1 target {cw_gap:.5f} "
              f"(asserted > 0.5)")
        print(f"[(S=)] ||S_fieldsweep - S_study||/||S_study|| = {s_tie:.3e} (band "
              f"REUSE_REPRODUCTION_RTOL {g.REUSE_REPRODUCTION_RTOL:.1e}; predicted <= 1e-14); "
              f"repeat tune(k=1) hex-equal: {hex_same}; fres.port_ids == study port order")
        print(f"[(P)] drive 1 accepted {d1.accepted_power_w:.5e} W of available "
              f"{d1.available_power_w:.5e} W; drive 2 accepted {d2.accepted_power_w:.5e} W of "
              f"available {d2.available_power_w:.5e} W (0 <= acc <= avail (1 + "
              f"{g.EXACT_IDENTITY_RTOL:.0e}))")
        print(f"[off-block] total {t1.off_diagonal_residual:.6e} (gate {REC['off_block']:.6e}, "
              f"log :7364)")
        for m in range(9):
            z = t1.z_k[m]
            cm = f"{1.0 / (omega * z.imag):.5e} F" if z.imag > 0.0 else "n/a (Im <= 0)"
            print(f"[z_k] k = {m}: z = {z.real:+.5g}{z.imag:+.5g}j ohm; C_k = {cm}; "
                  f"anti fraction {anti[m]:.5f}")
        print(f"[z_k] gate z_0 / z_1 / z_2 (log :7361-7363, 4 decimals): "
              + " / ".join(f"{v.real:+.4f}{v.imag:+.4f}j" for v in g.WF12_STEP3_LAMBDA_ANTI))
        for k, c in ((1, c1c), (2, c2c)):
            print(f"[mode {k}] centre |B1+| {c['b1p']:.5e}, |B1-| {c['b1m']:.5e}, "
                  f"|B| {c['norm']:.5e} T per unit incident-wave amplitude weights")
        print(f"[mode] |B_2(0)|/|B_1(0)| {c2c['norm'] / c1c['norm']:.4e} (predicted << 1)")
        print(f"[time] mesh {res.mesh_time:.1f} s, sweep {res.sweep_time:.1f} s, solve total "
              f"{t_solve:.1f} s, field sweep {t_fsweep:.1f} s, superpose+project {t_proj:.1f} s, "
              f"XDMF {t_xdmf:.1f} s")
        print(f"[mem] summed peak {rss_sum / 2**30:.3f} GiB")
        print(f"[paraview] wrote {xdmf_path}")

    # ---- asserts ------------------------------------------------------------
    assert res.cells == g.HOLE_CELL_RECORD
    assert res.cells == npz["cells"]
    assert r0 <= g.REUSE_REPRODUCTION_RTOL, f"(R0) {r0:.3e}"
    assert r_i1 <= g.EXACT_IDENTITY_RTOL and r_i2 <= g.EXACT_IDENTITY_RTOL, (r_i1, r_i2)
    assert r_io <= g.EXACT_IDENTITY_RTOL and r_iz <= g.EXACT_IDENTITY_RTOL, (r_io, r_iz)
    assert t1.z_k[1].imag > 0.0 and t1.k == 1 and t1.frequency_hz == g.FREQUENCY_HZ
    assert abs(gap_a) <= g.REDUCTION_BAND, f"(a) {gap_a:+.6e}"
    # record checks (10^-(d-2)); (R0) rel and (a) gap are band-only
    assert _rel(t1.c_mode_farad, REC["c_mode_f"]) <= 1e-11
    assert _rel(t2.c_mode_farad, REC["c2_f"]) <= 1e-11
    assert _rel(t2.c_mode_farad, g.WF12_STEP3_C2_CU_F) <= 1e-11
    assert _rel(gap_c1, REC["gap_c1"]) <= 1e-5, f"(c1) gap {gap_c1:+.6e}"
    assert _rel(t1.off_diagonal_residual, REC["off_block"]) <= 1e-5
    for m in range(3):
        rec = g.WF12_STEP3_LAMBDA_ANTI[m]
        z = t1.z_k[m]
        assert abs(z.real - rec.real) <= 1e-2 * abs(rec.real), (m, z, rec)
        assert abs(z.imag - rec.imag) <= 1e-4 * abs(rec.imag), (m, z, rec)
    assert abs(gap_c1) > g.REDUCTION_BAND, f"(c1) blind {gap_c1:+.6e}"
    assert w_gap[1] <= g.EXACT_IDENTITY_RTOL and w_gap[2] <= g.EXACT_IDENTITY_RTOL, w_gap
    assert cw_gap > 0.5, f"(cW) {cw_gap}"
    assert s_tie <= g.REUSE_REPRODUCTION_RTOL, f"(S=) {s_tie:.3e}"
    assert hex_same, "repeat tune(k=1) differs"
    assert ok_p, "(P) passivity"

    p0(f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
