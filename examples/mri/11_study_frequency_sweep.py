"""Example (`EX-95`): a frequency sweep from one ``Study`` -- 10 / 64 / 128 MHz on one
mesh, one multi-frequency Touchstone file.

`WF-16` step 2: ``Study.solve_sweep`` solves the F-small four-port birdcage (exact
CAD port rectangles, kappa off, operator ruling 2026-10-02) at several frequencies
on **one mesh** and writes **one** ``.s4p`` with one block per frequency.  This
example builds it at degree 1 on ``-n 2`` -- the width and degree of item 192's
registered records -- and holds every block, in memory and read back through
``scikit-rf``, to the record at that frequency.

Every constant, band and helper is imported from the gate modules (the ANS-1
rule); the example states no sigma, epsilon_r, frequency, geometry dimension or
impedance of its own.  The three frequencies come from ``FREQUENCIES``.

Asserted:

* (W) ``comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH`` before any S comparison;
* (R') per frequency, ``_rel(S_i, P1[k]) <= LEG_D0_REPRODUCTION_BAND``;
* (T') the single ``.s4p`` read back through ``skrf.Network``: three points at
  ``FREQUENCIES`` within ``FREQ_ATOL_HZ``, four ports, max |S_skrf - S| <=
  ``ANCHOR_TOL``, and ``_rel(S_file, P1[k]) <= LEG_D0_REPRODUCTION_BAND``;
* (G') `PORT-11`'s reciprocity / passivity / class-spread / pooled-separation
  gates at each frequency;
* (M) one mesh: every result carries the same mesh object (cells printed);
* (X) a ``keep_fields=True`` sweep on the 64 MHz result's own problem reproduces
  that result's S inside ``LEG_D0_REPRODUCTION_BAND``, so the field written to
  XDMF belongs to the gated S.

Negative controls (asserted, no solve; read files or in-memory data):

* (c4) a frequency mix-up: for every i != j, ``_rel(S_file_i, P1[k_j])`` exceeds
  the band;
* (c5) block order: the file's blocks with 10 and 128 swapped miss their records;
* (c1') printed: each block against the degree-2 record at its frequency.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:11 -n 2 -t 300

Not claimed: resonance, tuning, degree 2, frequency dispersion (the gel is
constant-property across the sweep), an absolute S at any frequency, anything
`WF-16` does not claim.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL, FREQ_ATOL_HZ  # noqa: E402
from tests.package.test_wf16_study_four_port import _rel  # noqa: E402
from tests.package.test_wf16_study_sweep import KEYS, _port_study, _records  # noqa: E402
from tests.validation.test_ans4_exact_strip_records import (  # noqa: E402
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
    FREQUENCIES,
)
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    _circulant_classes,
    _class_spread,
)
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    ADJACENT_SPREAD_BAND,
)
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_11_study_frequency_sweep"
DEGREE = 1  # the records' degree (item 192); the Study default is 2.


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- (W) the records hold at their own width (OPS-41) -------------------
    assert comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 records are asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}; "
        f"got -n {comm.size}"
    )

    # ---- the Study and its sweep: one mesh, three frequencies ---------------
    study, _phantom = _port_study(DEGREE)
    geometry = study.geometry
    freqs = tuple(FREQUENCIES.values())
    t0 = time.perf_counter()
    sw = study.solve_sweep(freqs, comm)
    t_sweep_total = time.perf_counter() - t0
    r0 = sw.results[0]
    msh = r0.mesh

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    local_tags = np.unique(r0.cell_tags.values)
    all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist()) for t in part})
    named = {
        geometry.conductor_tag: f"copper-proxy conductor (sigma {geometry.conductor_sigma_s_per_m:g} S/m)",
        geometry.phantom_tag: "gel phantom",
    }
    air_tags = tuple(t for t in all_tags if t not in named)
    region_names = dict(named)
    for t in air_tags:
        region_names[t] = "air (PEC outer wall)"
    write_setup_figure(
        msh,
        r0.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=air_tags,
        translucent_tags=(geometry.phantom_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:11 -- one mesh, three frequencies: F-small birdcage (WF-16 step 2, degree 1)",
        comm=comm,
    )

    p1 = _records(1)
    p2 = _records(2)

    # ---- (M) one mesh --------------------------------------------------------
    same_mesh = all(r.mesh is r0.mesh for r in sw.results)

    # ---- (R') each S against item 192's record at its frequency -------------
    s_list = [np.asarray(r.s_matrix, dtype=np.complex128) for r in sw.results]
    rel_r = {k: _rel(s, p1[k]) for k, s in zip(KEYS, s_list)}
    c1 = {k: _rel(s, p2[k]) for k, s in zip(KEYS, s_list)}

    # ---- (G') PORT-11's gates per frequency ---------------------------------
    gates = {}
    for k, r in zip(KEYS, sw.results):
        s, z = r.s_matrix, r.z_matrix
        recip = float(np.linalg.norm(s - s.T) / np.linalg.norm(s))
        sigma_sv = np.linalg.svd(s, compute_uv=False)
        power = np.sum(np.abs(s) ** 2, axis=0)
        classes = _circulant_classes(z)
        spreads = {n: _class_spread(v) for n, v in classes.items()}
        pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
        worst = max(spreads.values())
        gates[k] = {
            "recip": recip,
            "sigma_max": float(np.max(sigma_sv)),
            "power": [float(v) for v in power],
            "max_power": float(np.max(power)),
            "spreads": spreads,
            "separation": pooled / worst if worst > 0.0 else np.inf,
        }

    # ---- (T') the one Touchstone file, read back through scikit-rf ----------
    OUTPUT_DIR.mkdir(exist_ok=True)
    s4p_path = OUTPUT_DIR / f"{BASENAME}.s4p"
    t0 = time.perf_counter()
    written = sw.write_touchstone(s4p_path)
    t_ts = time.perf_counter() - t0
    readings = None
    if comm.rank == 0:
        import skrf

        ntw = skrf.Network(str(written))
        sk = np.asarray(ntw.s, dtype=np.complex128)
        readings = {
            "n_f": int(len(ntw.f)),
            "nports": int(ntw.nports),
            "df": [abs(float(ntw.f[i]) - freqs[i]) for i in range(min(len(ntw.f), len(freqs)))],
            "d_skrf": [float(np.max(np.abs(sk[i] - s_list[i]))) for i in range(len(freqs))],
            "rel_file": [_rel(sk[i], p1[KEYS[i]]) for i in range(len(freqs))],
            "c4": {
                (i, j): _rel(sk[i], p1[KEYS[j]])
                for i in range(len(freqs))
                for j in range(len(freqs))
                if i != j
            },
            "c5": [_rel(sk[2], p1[KEYS[0]]), _rel(sk[1], p1[KEYS[1]]), _rel(sk[0], p1[KEYS[2]])],
        }
    readings = comm.bcast(readings, root=0)

    # ---- (X) the 64 MHz field, from that result's own problem and mesh ------
    r64 = sw.results[KEYS.index("64")]
    t0 = time.perf_counter()
    sweep64 = run_n_port_sparameter_sweep(
        r64.problem,
        r64.port_definitions,
        lumped_sheet_ports=r64.specs,
        lumped_sheet_facet_tags=r64.facet_tags,
        degree=r64.degree,
        keep_fields=True,
    )
    t_x = time.perf_counter() - t0
    s_x = np.asarray(sweep64.s_matrix, dtype=np.complex128)
    rel_x = _rel(s_x, np.asarray(r64.s_matrix, dtype=np.complex128))
    fields = sweep64.fields["P1"]
    e = fields.e_complex

    # Every collective on every rank.  The mesh is affine, so the N1curl Function
    # is interpolated, not written: |E| into P1, sigma into DG0.
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    e_abs = fem.Function(p1_space, name="E_abs")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1_space.element.interpolation_points)
    )
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()
    dg0 = fem.functionspace(msh, ("DG", 0))
    sigma_out = fem.Function(dg0, name="sigma")
    sigma_out.interpolate(fields.sigma_field)
    sigma_out.x.array[:] = np.real(sigma_out.x.array)
    sigma_out.x.scatter_forward()

    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        r64.cell_tags,
        {"E_abs": e_abs, "sigma": sigma_out},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- prints (rank 0, no collective inside) ------------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-95 -- a frequency sweep from one Study (F-small birdcage, WF-16 step 2)")
        print("=" * 72)
        print(
            f"[study] degree {study.degree}, -n {comm.size}, frequencies "
            + " / ".join(f"{k} MHz = {f:.6e} Hz" for k, f in zip(KEYS, freqs))
        )
        print(
            f"[(M)] one mesh {same_mesh}; {r0.cells} cells (printed; version-tagged); "
            f"solve_sweep {t_sweep_total:.2f} s (mesh {sw.mesh_time:.2f} s; sweeps "
            + " / ".join(f"{k} MHz {t:.2f} s" for k, t in zip(KEYS, sw.sweep_times))
            + f"); Touchstone {t_ts:.3f} s; keep_fields sweep {t_x:.2f} s"
        )
        print(f"[(W)] comm.size {comm.size} == record width {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}")
        for k in KEYS:
            print(
                f"[(R')] {k} MHz ||S - P1||/||P1|| {rel_r[k]:.3e} (band "
                f"LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:g})"
            )
        for k in KEYS:
            g = gates[k]
            print(
                f"[(G')] {k} MHz reciprocity {g['recip']:.3e} (band {RECIPROCITY_BAND:g}); "
                f"sigma_max {g['sigma_max']:.12f}; max column power {g['max_power']:.12f} "
                f"(<= 1 + {PASSIVITY_SIGMA_TOLERANCE:g}); spreads "
                + ", ".join(f"{n} {v * 100:.4f}%" for n, v in g["spreads"].items())
                + f" (band {ADJACENT_SPREAD_BAND:g}); separation {g['separation']:.2f}x "
                f"(floor {POOLED_SEPARATION_FLOOR:g})"
            )
        min_sep = min(g["separation"] for g in gates.values())
        print(f"[(G')] minimum separation over the sweep {min_sep:.2f}x")
        print(
            f"[(T')] {readings['n_f']} points, nports {readings['nports']}; |df| "
            + " / ".join(f"{d:.3e}" for d in readings["df"])
            + f" Hz (FREQ_ATOL_HZ {FREQ_ATOL_HZ:g}); max|S_skrf - S| "
            + " / ".join(f"{d:.3e}" for d in readings["d_skrf"])
            + f" (ANCHOR_TOL {ANCHOR_TOL:g}); ||S_file - P1||/||P1|| "
            + " / ".join(f"{d:.3e}" for d in readings["rel_file"])
        )
        print(
            "[(c4)] frequency mix-up, file block i vs record j (must exceed "
            f"{LEG_D0_REPRODUCTION_BAND:g}): "
            + ", ".join(
                f"{KEYS[i]}->{KEYS[j]} {v:.4f}" for (i, j), v in readings["c4"].items()
            )
            + f"; minimum {min(readings['c4'].values()):.4f}"
        )
        print(
            "[(c5)] block order, file with 10 and 128 swapped vs records (10 / 64 / 128, "
            "64 unmoved): "
            + " / ".join(f"{v:.4f}" for v in readings["c5"])
        )
        print(
            "[(c1')] each block vs the degree-2 record at its frequency (printed): "
            + " / ".join(f"{k} MHz {c1[k]:.4f}" for k in KEYS)
        )
        print(
            f"[(X)] 64 MHz keep_fields sweep S vs that result's S {rel_x:.3e} (band "
            f"{LEG_D0_REPRODUCTION_BAND:g}; printed, not a record)"
        )
        print(f"[paraview] wrote {xdmf_path}; Touchstone {written}")

    assert same_mesh, "(M) the sweep results do not share one mesh"
    assert readings["n_f"] == len(freqs), f"(T') {readings['n_f']} points"
    assert readings["nports"] == 4, f"(T') nports {readings['nports']}"
    for i, k in enumerate(KEYS):
        assert readings["df"][i] <= FREQ_ATOL_HZ, f"(T') {k} MHz df {readings['df'][i]:.3e}"
        assert readings["d_skrf"][i] <= ANCHOR_TOL, f"(T') {k} |S_skrf - S| {readings['d_skrf'][i]:.3e}"
        assert readings["rel_file"][i] <= LEG_D0_REPRODUCTION_BAND, (
            f"(T') {k} file vs P1 {readings['rel_file'][i]:.3e}"
        )
        assert rel_r[k] <= LEG_D0_REPRODUCTION_BAND, f"(R') {k} rel {rel_r[k]:.3e}"
        g = gates[k]
        assert g["recip"] <= RECIPROCITY_BAND, f"(G') {k} reciprocity {g['recip']:.3e}"
        assert g["sigma_max"] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, f"(G') {k} sigma_max"
        for v in g["power"]:
            assert v <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, f"(G') {k} column power {v!r}"
        for name, v in g["spreads"].items():
            assert v <= ADJACENT_SPREAD_BAND, f"(G') {k} spread {name} {v:.3e}"
        assert g["separation"] >= POOLED_SEPARATION_FLOOR, f"(G') {k} separation"
    for (i, j), v in readings["c4"].items():
        assert v > LEG_D0_REPRODUCTION_BAND, f"(c4) blind {KEYS[i]}->{KEYS[j]}: {v:.3e}"
    # (c5) swaps 10 and 128 only: the two moved blocks must miss; the 64 MHz block
    # stays in place and must still hold.
    for idx in (0, 2):
        assert readings["c5"][idx] > LEG_D0_REPRODUCTION_BAND, (
            f"(c5) blind: {readings['c5'][idx]:.3e}"
        )
    assert readings["c5"][1] <= LEG_D0_REPRODUCTION_BAND, f"(c5) 64 MHz {readings['c5'][1]:.3e}"
    assert rel_x <= LEG_D0_REPRODUCTION_BAND, f"(X) sweep S vs result S {rel_x:.3e}"

    if comm.rank == 0:
        print(
            f"\nAll anchors and controls hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
