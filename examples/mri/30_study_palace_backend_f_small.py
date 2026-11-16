"""Example (`EX-134`): the F-small four-port case through ``Study.solve(backend="palace")``.

A three-stage workflow, run end to end into a fresh case directory:

1. (this script, stage 1) ``Study.solve(comm, backend="palace", case_dir=...)`` writes
   the Palace case and raises ``PalaceRunRequired``;
2. a separate host-side ``fem-em-palace`` harness window runs the printed command
   (the solver container has no docker socket and cannot launch Palace);
3. (this script, stage 3) the same call reads S back from Palace's output.

The stage is chosen by what is on disk: no Palace output in the case directory means
stage 1, an existing output means stage 3.  The writer is serial: ``-n 1`` only.

Every record, band and helper is imported from the gate modules (the ANS-1 rule).  The
Palace figures are a cross-solver identity on our own discretisation at matched order,
not an accuracy level.

Asserted (stage 1): (K) serial; (Arg) the printed argv / cwd; (W-case) the written case
equals item 517's U16 case and `OPS-61`'s inputs byte for byte; (c-stage) no output
exists after the raise.  Asserted (stage 3): (K), (W-case), (R) the result record,
(E2) / (E2-u16) the read equals the recorded Palace reads, (P) `PORT-11`'s identities,
(Tie2) the per-class tie to our degree-2 record, (c-order) the read is not Palace
order 1.  Printed only: (c-freq), the Palace wall time, the sha256 manifest.

Not claimed: order 1, other frequencies or rungs, fields from Palace, a Study that
launches Palace, accuracy, tuning, Q, any AED comparison.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import (  # noqa: E402
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    Study,
)
from fem_em_solver.backends.palace import (  # noqa: E402
    PalaceRunRequired,
    PalaceStudyResult,
    read_port_s,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.materials import GelledSalinePhantomMaterial  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.package.test_arch1_palace_backend import (  # noqa: E402
    C_ORDER_FLOOR,
    D2_CASE,
    OPS61,
    TIE2_BAND,
    TIE2_RECORDED,
    _class_means,
    _flat,
    _fro_rel,
    _p_rec,
)
from tests.validation.test_ans4_exact_strip_records import (  # noqa: E402
    ANS4_EXACT_STRIP_S_P2_X1,
)
from tests.validation.test_lossy_sphere_fullwave import (  # noqa: E402
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_ops61_exact_strip_ports import (  # noqa: E402
    ANS4_EXACT_STRIP_CELL_RECORD,
)
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    LEG_D0_REPRODUCTION_BAND,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_30_study_palace_backend_f_small"
CASE_DIR = OUTPUT_DIR / f"{BASENAME}_case"

# Records printed beside the readings (never asserted by equality), each with its log:line.
REC = {
    "E2": (0.0, "20261018T172228Z_ARCH-1-U16c-W1.log (T4)"),
    "reciprocity": (9.848e-11, "20261018T172228Z_ARCH-1-U16c-W1.log (T4)"),
    "sigma_max": (0.999819887373, "20261018T172228Z_ARCH-1-U16c-W1.log (T4)"),
    "c_order_pct": ((6.6674, 2.5442, 4.5224), "20261018T172228Z_ARCH-1-U16c-W1.log (T4)"),
    "c_freq": (8.515e-01, "20261018T172228Z_ARCH-1-U16c-W1.log (T4)"),
    "palace_total_s": (261.6, "20261018T140040Z_ARCH-1-U16b-W1-palace-d2.log"),
    "write_s": (27.1, "20261018T172228Z_ARCH-1-U16c-W1.log:1833"),
}


def _build_study() -> Study:
    return Study(
        FourPortBirdcageGeometry.f_small(),
        GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ),
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        degree=2,
    )


def _check_w_case(comm) -> dict:
    """(W-case): the written case is item 517's U16 case and OPS-61's inputs, byte for byte."""
    out = {}
    out["manifest_equal"] = (CASE_DIR / "case.sha256").read_text() == (
        D2_CASE / "case.sha256"
    ).read_text()
    out["config_bytes_equal"] = (CASE_DIR / "case.json").read_bytes() == (
        D2_CASE / "case.json"
    ).read_bytes()
    ours = _flat(json.loads((CASE_DIR / "case.json").read_text()))
    ref = _flat(json.loads((OPS61 / "f05x_o2.json").read_text()))
    out["cfg_diff"] = {k for k in set(ours) | set(ref) if ours.get(k) != ref.get(k)}
    out["mesh_bytes_equal"] = (CASE_DIR / "mesh.msh").read_bytes() == (
        OPS61 / "ops61_ans4_x1_strip.msh"
    ).read_bytes()
    assert out["manifest_equal"], "(W-case) manifest differs from the U16 case"
    assert out["config_bytes_equal"], "(W-case) case config bytes differ from the U16 case"
    assert out["cfg_diff"] == {"Problem.Output", "Model.Mesh"}, (
        f"(W-case) config diff vs OPS-61 {sorted(out['cfg_diff'])}"
    )
    assert out["mesh_bytes_equal"], "(W-case) mesh bytes differ from OPS-61's"
    return out


def stage1(comm) -> None:
    started = time.perf_counter()
    assert comm.size == 1, f"(K) the writer is serial: -n 1, got -n {comm.size}"
    study = _build_study()
    t0 = time.perf_counter()
    try:
        study.solve(comm, backend="palace", case_dir=CASE_DIR)
    except PalaceRunRequired as req:
        t_write = time.perf_counter() - t0
        assert req.argv == ["palace", "-np", "8", "case.json"], f"(Arg) argv {req.argv}"
        assert Path(req.cwd).resolve() == CASE_DIR.resolve(), f"(Arg) cwd {req.cwd}"
        assert not (CASE_DIR / "postpro").exists(), "(c-stage) postpro exists after the raise"
        w = _check_w_case(comm)
        adopt_host_ownership(OUTPUT_DIR, comm=comm)
        if comm.rank == 0:
            print("=" * 72)
            print("EX-134 -- stage 1 of 3: write the Palace case through the Study")
            print("=" * 72)
            print(f"[(K)] comm.size {comm.size} (serial writer)")
            print(
                f"[(Arg)] argv {req.argv}; cwd {req.cwd} (resolves to the case directory)"
            )
            print(
                f"[(W-case)] manifest equal {w['manifest_equal']}; config bytes equal "
                f"{w['config_bytes_equal']}; config diff vs OPS-61 {sorted(w['cfg_diff'])}; "
                f"mesh bytes equal {w['mesh_bytes_equal']}"
            )
            print(
                f"[(c-stage)] PalaceRunRequired raised, no postpro directory after it; "
                f"write stage {t_write:.1f} s (record {REC['write_s'][0]} s, "
                f"{REC['write_s'][1]})"
            )
            print(
                "[next] run Palace as a separate host-side window (the solver container "
                "cannot launch it):"
            )
            print(
                "  docker compose -f docker/docker-compose.yml --profile palace run --rm -T "
                "fem-em-palace bash -lc 'cd /workspace/examples/mri/paraview_output/"
                f"{CASE_DIR.name} && timeout -k 30 590 palace -np 8 case.json'"
            )
            print(
                f"\n[stage 1 of 3] All stage-1 assertions hold. Total elapsed "
                f"{time.perf_counter() - started:.1f} s.",
                flush=True,
            )
        return
    raise AssertionError("(c-stage) the call returned without PalaceRunRequired")


def stage3(comm) -> None:
    started = time.perf_counter()
    assert comm.size == 1, f"(K) the reader is serial: -n 1, got -n {comm.size}"
    study = _build_study()
    t0 = time.perf_counter()
    res = study.solve(comm, backend="palace", case_dir=CASE_DIR)
    t_read = time.perf_counter() - t0
    w = _check_w_case(comm)

    # ---- (R) the result record ------------------------------------------------
    assert isinstance(res, PalaceStudyResult), type(res)
    assert res.backend == "palace" and res.order == 2, (res.backend, res.order)
    assert res.cells == ANS4_EXACT_STRIP_CELL_RECORD, res.cells
    assert res.frequency_hz == FREQUENCY_64_HZ, res.frequency_hz
    s = np.asarray(res.s_matrix)

    # ---- (E2) / (E2-u16) the read equals the recorded Palace reads -----------
    s_ops61, _ = read_port_s(OPS61 / "postpro_f05x_o2" / "port-S.csv")
    s_u16, _ = read_port_s(D2_CASE / "postpro" / "port-S.csv")
    e2 = _fro_rel(s, s_ops61)
    e2_max = float(np.max(np.abs(s - s_ops61)))
    e2u = _fro_rel(s, s_u16)
    e2u_max = float(np.max(np.abs(s - s_u16)))
    assert e2 <= LEG_D0_REPRODUCTION_BAND, f"(E2) {e2:.3e}"
    assert e2u <= LEG_D0_REPRODUCTION_BAND, f"(E2-u16) {e2u:.3e}"

    # ---- (P) PORT-11's identities (prints its own line) ----------------------
    _p_rec(s)

    # ---- (Tie2) / (c-order) / (c-freq) ---------------------------------------
    m = _class_means(s)
    ref = _class_means(ANS4_EXACT_STRIP_S_P2_X1["64"])
    tie2 = {k: float(abs(m[k] - ref[k]) / abs(ref[k])) for k in ref}
    for k, v in tie2.items():
        assert v <= TIE2_BAND, f"(Tie2) {k} {v:.3e} > {TIE2_BAND:.1e}"
    o1 = _class_means(read_port_s(OPS61 / "postpro_f05x" / "port-S.csv")[0])
    c_order = {k: float(abs(m[k] - o1[k]) / abs(o1[k])) for k in o1}
    for k, v in c_order.items():
        assert v > C_ORDER_FLOOR, f"(c-order) blind: {k} {v:.3e}"
    s_f10, _ = read_port_s(OPS61 / "postpro_f05x_o2_f10" / "port-S.csv")
    c_freq = _fro_rel(s, s_f10)

    # ---- the DolfinX side for ParaView: geometry and tags only ---------------
    t0 = time.perf_counter()
    geometry = study.geometry
    built = geometry.build(comm)
    t_build = time.perf_counter() - t0
    msh, cell_tags, facet_tags = built["mesh"], built["cell_tags"], built["facet_tags"]

    from dolfinx import fem

    dg0 = fem.functionspace(msh, ("DG", 0))
    sigma_out = fem.Function(dg0, name="sigma_palace")
    sigma_out.x.array[:] = 0.0
    idx = np.asarray(cell_tags.indices)
    val = np.asarray(cell_tags.values)
    sigma_by_tag = {
        geometry.conductor_tag: geometry.conductor_sigma_s_per_m,
        geometry.phantom_tag: SALINE_SIGMA,
    }
    for c, tag in zip(idx, val):
        if int(tag) in sigma_by_tag:
            sigma_out.x.array[dg0.dofmap.cell_dofs(int(c))[0]] = sigma_by_tag[int(tag)]
    sigma_out.x.scatter_forward()

    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cell_tags,
        {"sigma_palace": sigma_out},
        comm=comm,
        facet_tags=facet_tags,
    )

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    local_tags = np.unique(cell_tags.values)
    all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist()) for t in part})
    named = {
        geometry.conductor_tag: (
            f"copper-proxy conductor (sigma {geometry.conductor_sigma_s_per_m:g} S/m)"
        ),
        geometry.phantom_tag: "gel phantom",
    }
    air_tags = tuple(t for t in all_tags if t not in named)
    region_names = dict(named)
    for t in air_tags:
        region_names[t] = "air (PEC outer wall)"
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=air_tags,
        translucent_tags=(geometry.phantom_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:30 -- F-small four-port birdcage through backend='palace' (degree 2)",
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- CSV (rank 0) ----------------------------------------------------------
    csv_path = OUTPUT_DIR / f"{BASENAME}_s_classes.csv"
    if comm.rank == 0:
        lines = [
            "class,palace_abs,palace_arg_rad,p2_record_abs,p2_record_arg_rad,"
            "tie2_gap,tie2_recorded_gap"
        ]
        for (k, v), rec in zip(tie2.items(), TIE2_RECORDED):
            lines.append(
                f"{k},{abs(m[k]):.6g},{np.angle(m[k]):.6g},{abs(ref[k]):.6g},"
                f"{np.angle(ref[k]):.6g},{v:.6g},{rec:.6g}"
            )
        csv_path.write_text("\n".join(lines) + "\n")

    # ---- prints (rank 0, no collective inside) -------------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-134 -- stage 3 of 3: read S back through the Study (backend='palace')")
        print("=" * 72)
        print(f"[(K)] comm.size {comm.size}; read stage {t_read:.1f} s")
        print(
            f"[(W-case)] manifest equal {w['manifest_equal']}; config bytes equal "
            f"{w['config_bytes_equal']}; config diff vs OPS-61 {sorted(w['cfg_diff'])}; "
            f"mesh bytes equal {w['mesh_bytes_equal']}"
        )
        print(
            f"[(R)] {type(res).__name__}; backend {res.backend}; order {res.order}; "
            f"cells {res.cells} (record {ANS4_EXACT_STRIP_CELL_RECORD}); "
            f"f = {res.frequency_hz:.6e} Hz"
        )
        print(
            f"[(E2)] Frobenius-rel vs postpro of the OPS-61 order-2 run {e2:.3e} "
            f"(predicted {REC['E2'][0]:.3e}; band {LEG_D0_REPRODUCTION_BAND:g}); "
            f"entry max {e2_max:.3e}"
        )
        print(
            f"[(E2-u16)] Frobenius-rel vs the U16 case output {e2u:.3e} "
            f"(predicted 0.000e+00; band {LEG_D0_REPRODUCTION_BAND:g}); entry max "
            f"{e2u_max:.3e}"
        )
        print(
            f"[(Tie2)] "
            + ", ".join(f"{k} {v:.3e}" for k, v in tie2.items())
            + f" (band {TIE2_BAND:.1e}; recorded {TIE2_RECORDED})"
        )
        print(
            "[(c-order)] vs Palace order 1: "
            + ", ".join(f"{k} {v * 100:.4f}%" for k, v in c_order.items())
            + f" (floor {C_ORDER_FLOOR:g} relative; recorded {REC['c_order_pct'][0]} %)"
        )
        print(
            f"[(c-freq)] vs the 10 MHz Palace output {c_freq:.4g} "
            f"(predicted ~{REC['c_freq'][0]:.4g}; printed, not asserted)"
        )
        print(
            f"[printed] palace_total_s {res.palace_total_s} (record {REC['palace_total_s'][0]} s); "
            f"sha256 {res.sha256}"
        )
        print(f"[paraview] extra geometry build {t_build:.1f} s; wrote {xdmf_path}; {csv_path}")
        print(
            f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )


def main() -> None:
    comm = MPI.COMM_WORLD
    if (CASE_DIR / "postpro" / "port-S.csv").exists():
        stage3(comm)
    else:
        stage1(comm)


if __name__ == "__main__":
    main()
