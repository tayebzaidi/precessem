"""Example (`EX-115`): ``femem run`` with a top-level ``tune:`` section -- one solve,
``summary["tune"]`` is the same Study's ``tune()``, the written ``.s4p`` alone gives the
same capacitance, and the field of the mode the file tuned.

`WF-1` step 11 (`U49`) let a YAML file carry ``tune: {k: 1}``: ``run_study_file`` then calls
``Study.tune(k)`` on the *solved* Study (no second solve) and writes the result into
``summary.json``.  No example ran such a file.  This one runs the tracked
``tests/package/wf1/f_small_hole_pec_tune.yaml`` in place (the exact-rectangle PEC hole
F-small birdcage, degree 1, 64 MHz, ``-n 2``, kappa off) through ``cli.main(["run", ...])``
in-process, captures the returned ``StudyRun`` with a call-through wrapper on
``fem_em_solver.study_yaml.run_study_file`` (restored in ``finally``), and shows:

* the gate's anchors, re-asserted: (R0) cells + class means, (I=) the summary is the same
  Study's ``tune()`` by ``float.hex`` (in memory and read back from ``summary.json``),
  (a) ``C_mode`` vs ``TH17_STEP3_C_MODE_RECORD_F`` within ``REDUCTION_BAND``;
* (Z) the keys tie together from the summary alone: ``C * 2 pi f * Im z_k = 1``;
* (D) an omitted ``k`` (``tune: {}``) is ``Study.tune``'s default, k = 1;
* (T) the written ``.s4p`` alone, read through scikit-rf and put on a fresh mesh-free
  Study, returns the summary's C;
* (Once) one Study solve;
* the field of the mode the file tuned (``Study.b1_plus()`` is a second, disclosed,
  field-keeping build and sweep, run after every ``tune`` read), written to ParaView.

Negative controls (asserted): (c1) ``tune: {k: 2}`` misses the record; (cX)
``malformed_tune_k_zero.yaml`` refuses at ``tune.k``, rc (2, 2), nothing written.

Every record, band and helper is imported from the gate module (the ANS-1 rule); record
digits sit only in the example-local ``REC`` dict with their ``log:line``.  New values
print at <= 6 significant digits; the summary's S / Z / tune block is never printed whole.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:21 -n 2 -t 420

Not claimed: the installed console script; ``tune`` at degree 2, on copper, on the solid
four-port or on F-human; k = 2 as a gated mode; a tuned f0, Q, match or any absolute
S11 / Z_in / C figure; Touchstone accuracy beyond the round trip; an external check.
"""

from __future__ import annotations

import contextlib
import copy
import io
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

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
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.package.test_wf1_yaml_tune as g  # noqa: E402
from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL, FREQ_ATOL_HZ  # noqa: E402
from tests.validation import test_port_lumped_rlc_termination as rlc  # noqa: E402
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG  # noqa: E402
from tests.validation.test_th17_birdcage_eigenmodes import _four_port_dft  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_21_femem_run_tune_section"
RUN_DIR = OUTPUT_DIR / "mri_21_run"
BAD_DIR = OUTPUT_DIR / "mri_21_bad"

# Gate-log records (tracked), compared at 10^-(d-2) with d the printed digits.
REC = {
    # docs/testing/logs/20261014T110326Z_WF-1-s11-W1.log:1846 (a) C_mode and gap
    "c_mode": (6.626796421e-11, 1.0e-8),
    "gap_a": (-2.976552e-04, 1.0e-5),
    # :1850 (c1) k = 2
    "c_k2": (9.873468900e-11, 1.0e-8),
    "gap_c1": (+4.894874e-01, 1.0e-5),
}
REC_CELLS = 80798  # :1834 (asserted against the gate's ANS6_EXACT_STRIP_HOLE_CELL_RECORD)

# mri:19 mode-field readings (docs/testing/logs/20261014T111653Z_EX-113.log:3627-3630),
# printed beside, never asserted.
MRI19 = {
    1: (3.33666e-10, 3.68784e-08, 5.21563e-08),
    2: (2.67675e-10, 2.63168e-10, 5.36026e-10),
    "ratio": 1.0277e-02,
}


def _within(value, key):
    rec, rtol = REC[key]
    return abs(value / rec - 1.0) <= rtol


def _hex_eq(a, b):
    return float(a).hex() == float(b).hex()


def _tune_equal(a, b):
    """Bitwise equality of two summary ``tune`` blocks (dicts of floats / lists)."""
    if a["k"] != b["k"] or len(a["z_k"]) != len(b["z_k"]):
        return False
    ok = (
        _hex_eq(a["frequency_hz"], b["frequency_hz"])
        and _hex_eq(a["c_mode_farad"], b["c_mode_farad"])
        and _hex_eq(a["off_diagonal_residual"], b["off_diagonal_residual"])
    )
    for (r1, i1), (r2, i2) in zip(a["z_k"], b["z_k"]):
        ok = ok and _hex_eq(r1, r2) and _hex_eq(i1, i2)
    return bool(ok)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- structural asserts -------------------------------------------------
    assert comm.size == g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 record is asserted at -n {g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; "
        f"got -n {comm.size}"
    )
    if not np.issubdtype(default_scalar_type, np.complexfloating):
        raise RuntimeError("complex build required (source dolfinx-complex-mode)")
    assert g.REDUCTION_BAND is rlc.REDUCTION_BAND
    parsed_outputs = g.study_yaml._load(g.TUNE_FILE)[2]
    assert parsed_outputs == {"tune": {"k": 1}}, parsed_outputs
    assert type(parsed_outputs["tune"]["k"]) is int
    probe = g.load_study(g.TUNE_FILE)
    assert probe.degree == 1 and probe.conductor is None
    omega = float(probe.omega)

    # ---- (cX) the malformed file refuses at tune.k; no mesh, nothing written -----
    try:
        g.load_study(g.MALFORMED_TUNE)
    except g.StudyConfigError as exc:
        cx_field = exc.field
        cx_type_ok = type(exc) is g.StudyConfigError
        cx_msg = str(exc)
    else:
        raise AssertionError("(cX) malformed_tune_k_zero.yaml did not refuse")
    OUTPUT_DIR.mkdir(exist_ok=True)
    g._fresh_dir(BAD_DIR, comm)
    rc_bad_run = g.cli.main(["run", str(g.MALFORMED_TUNE), "--output-dir", str(BAD_DIR)])
    rc_bad_check = g.cli.main(["check", str(g.MALFORMED_TUNE)])
    cx_run, cx_check = g._rc_all(rc_bad_run, comm), g._rc_all(rc_bad_check, comm)
    cx_exists = comm.allreduce(bool(BAD_DIR.exists()), op=MPI.LOR)

    # ---- femem check on the tune file (inputs only) --------------------------------
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc_check = g.cli.main(["check", str(g.TUNE_FILE)])
    check_all = g._rc_all(rc_check, comm)
    family_line = next(
        (ln for ln in buf.getvalue().splitlines() if ln.startswith("family:")), "")
    family_line = comm.bcast(family_line, root=0)

    # ---- (1) femem run in-process, captured by a call-through wrapper --------------
    g._fresh_dir(RUN_DIR, comm)
    captured = []
    original = g.study_yaml.run_study_file

    def _call_through(*args, **kwargs):
        out = original(*args, **kwargs)
        captured.append(out)
        return out

    g.study_yaml.run_study_file = _call_through
    t0 = time.perf_counter()
    try:
        rc_run = g.cli.main(["run", str(g.TUNE_FILE), "--output-dir", str(RUN_DIR)])
    finally:
        g.study_yaml.run_study_file = original
    t_cli = time.perf_counter() - t0
    assert g.study_yaml.run_study_file is original
    rc_all = g._rc_all(rc_run, comm)
    n_solves = len(captured)
    assert n_solves == 1, f"(Once) wrapper called {n_solves} times"
    run = captured[0]
    comm.Barrier()

    with open(RUN_DIR / "summary.json") as fh:
        disk = json.load(fh)
    listing = sorted((p.name, p.stat().st_size) for p in RUN_DIR.iterdir())

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    msh0 = run.result.mesh
    local_tags = np.unique(run.result.cell_tags.values)
    all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist()) for t in part})
    air_tags = tuple(t for t in all_tags if t != PHANTOM_CELL_TAG)
    names = {PHANTOM_CELL_TAG: "gel phantom"}
    for t in air_tags:
        names[t] = "air (outer box PEC)"
    FIGURE_DIR.mkdir(exist_ok=True)
    write_setup_figure(
        msh0, run.result.cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=names, hide_tags=air_tags,
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title=("mri:21 -- femem run with a tune: section, exact-rectangle PEC hole "
               "F-small birdcage, four gap ports, phantom (WF-1 step 11)"),
        comm=comm,
    )

    # ---- every tune read first (no solve), the field-keeping build after -----------
    summary = run.summary
    t = summary["tune"]
    kt = int(t["k"])
    res0 = run.study.result
    s_run = np.asarray(run.result.s_matrix, dtype=np.complex128)
    rel_means = g._rel_means(
        g._class_means(run.result.s_matrix),
        g.ANS6_EXACT_STRIP_CLASS_MEANS_P1[g.KEY]["pec"],
    )

    # (I=) summary vs the same Study's tune(), bitwise
    api = run.study.tune(k=1, comm=comm)
    same_after_api = run.study.result is res0
    i_mem = (
        _hex_eq(t["c_mode_farad"], api.c_mode_farad)
        and _hex_eq(t["off_diagonal_residual"], api.off_diagonal_residual)
        and len(t["z_k"]) == len(api.z_k)
        and all(
            _hex_eq(zr, complex(z).real) and _hex_eq(zi, complex(z).imag)
            for (zr, zi), z in zip(t["z_k"], api.z_k)
        )
    )
    i_disk = _tune_equal(disk["tune"], t)

    # (a), (c1)
    c_mode = float(t["c_mode_farad"])
    gap_a = c_mode / g.TH17_STEP3_C_MODE_RECORD_F - 1.0
    m2 = copy.deepcopy(g._mapping(g.TUNE_FILE))
    m2["tune"] = {"k": 2}
    spec2 = g._parse_mapping(m2, "<tune k2>")[2]["tune"]
    t2 = g._tune_summary(run.study, spec2, comm)
    gap_c1 = t2["c_mode_farad"] / g.TH17_STEP3_C_MODE_RECORD_F - 1.0

    # (Z) the keys tie together, from the summary alone (index = the summary's own k)
    z_im_k = float(t["z_k"][kt][1])
    z_tie = abs(c_mode * 2.0 * np.pi * float(t["frequency_hz"]) * z_im_k - 1.0)

    # (D) an omitted k is the default
    m0 = copy.deepcopy(g._mapping(g.TUNE_FILE))
    m0["tune"] = {}
    spec0 = g._parse_mapping(m0, "<tune {}>")[2]["tune"]
    t0s = g._tune_summary(run.study, spec0, comm)
    d_equal = _tune_equal(t0s, t)

    # (T) the written .s4p alone
    s4p = Path(run.paths["touchstone"])
    reading = None
    if comm.rank == 0:
        import skrf

        ntw = skrf.Network(str(s4p))
        reading = {
            "s": np.asarray(ntw.s, dtype=np.complex128),
            "f": np.asarray(ntw.f, dtype=float),
            "nports": int(ntw.nports),
        }
    reading = comm.bcast(reading, root=0)
    s_file = reading["s"]
    t_nf = int(reading["f"].shape[0])
    t_df = abs(float(reading["f"][0]) - g.FREQUENCY_64_HZ) if t_nf >= 1 else float("inf")
    t_ds = (float(np.max(np.abs(s_file[0] - s_run))) if s_file.shape == (1, 4, 4)
            else float("inf"))
    st = g.load_study(g.TUNE_FILE)
    st.result = SimpleNamespace(s_matrix=s_file[0])
    c_file = float(st.tune(k=kt, comm=comm).c_mode_farad)
    t_rel = abs(c_file / c_mode - 1.0)

    same_after_reads = run.study.result is res0

    # ---- prints and asserts for everything above ------------------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-115 -- femem run with a tune: section (mri:21)")
        print("=" * 72)
        print(f"[(W)] comm.size {comm.size} == record width "
              f"{g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; f = {float(t['frequency_hz']):.6e} Hz")
        print(f"[outputs] parsed outputs {parsed_outputs!r}")
        print(f"[check] femem check rc {check_all}; stdout line {family_line!r}")
        print(f"[(cX)] {g.MALFORMED_TUNE.name}: field {cx_field!r}, StudyConfigError exact "
              f"{cx_type_ok}; {cx_msg}")
        print(f"[(cX)] run rc (min, max) {cx_run}, check rc {cx_check}, out dir exists "
              f"{cx_exists} (expected (2, 2), (2, 2), False)")
        print(f"[run] cli.main rc {rc_all}; wall {t_cli:.1f} s; summary elapsed_s "
              f"{float(summary.get('elapsed_s', float('nan'))):.1f}; wrapper calls {n_solves}")
        print("[run] directory: " + ", ".join(f"{n} {sz} B" for n, sz in listing))
        print(f"[(R0)] cells {summary['cells']} (record {g.ANS6_EXACT_STRIP_HOLE_CELL_RECORD}); "
              "class-mean rel vs 'pec' record: "
              + "  ".join(f"{c} {rel_means[c]:.3e}" for c in g.CLASSES)
              + f" (band {g.LEG_D0_REPRODUCTION_BAND:.0e}; gate 5.052e-11 / 3.878e-11 / "
              "5.875e-11)")
        print(f"[(I=)] summary == same Study's tune() bitwise: in memory {i_mem}, "
              f"on disk {i_disk}; study.result unchanged {same_after_api}; k {kt}; "
              f"frequency == FREQUENCY_64_HZ {_hex_eq(t['frequency_hz'], g.FREQUENCY_64_HZ)}")
        print(f"[(a)] C_mode = {c_mode:.9e} F (record {REC['c_mode'][0]:.9e}); gap "
              f"{gap_a:+.6e} (gate {REC['gap_a'][0]:+.6e}; item 429 print "
              f"{g.ITEM429_A_PRINT:+.6e}); band {g.REDUCTION_BAND:.1e}")
        print(f"[(c1)] k = 2: C = {t2['c_mode_farad']:.9e} F (record {REC['c_k2'][0]:.9e}); "
              f"gap {gap_c1:+.6e} (gate {REC['gap_c1'][0]:+.6e}); "
              f"{abs(gap_c1) / g.REDUCTION_BAND:.1f}x band")
        print(f"[(Z)] |C * 2 pi f * Im z_k - 1| = {z_tie:.3e} (k = {kt}; band "
              f"{EXACT_IDENTITY_RTOL:g}; predicted <= 1e-15)")
        print(f"[(D)] tune: {{}} parses to {spec0!r}; block == k = 1 block bitwise {d_equal}; "
              f"k {t0s['k']}")
        print(f"[(T)] {s4p.name}: {t_nf} frequency, {reading['nports']} ports, |f - f_64| "
              f"{t_df:.3e} Hz (atol {FREQ_ATOL_HZ:g}); max|S_file - S_run| {t_ds:.3e} "
              f"(ANCHOR_TOL {ANCHOR_TOL:g}; mri:11 <= 6.855e-14)")
        print(f"[(T)] C from the file's S alone: |C_file / C_summary - 1| = {t_rel:.3e} "
              f"(band {g.LEG_D0_REPRODUCTION_BAND:.0e}; predicted ~1e-12)")
        print(f"[(Once)] wrapper calls {n_solves}; study.result is the first result after "
              f"(I=), (D), (c1), (T): {same_after_reads}")
        print("[tune block] k, f(Hz), residual(4 digits): "
              f"{kt}, {float(t['frequency_hz']):.6g}, {float(t['off_diagonal_residual']):.4g}")
        print("  k    Re z_k (ohm)   Im z_k (ohm)   C_k = 1/(omega Im z_k) (F)")
        for j, (zr, zi) in enumerate(t["z_k"]):
            ck = f"{1.0 / (omega * float(zi)):.6g}" if float(zi) > 0.0 else ""
            print(f"  {j}   {float(zr):+.6g}   {float(zi):+.6g}   {ck}")

    assert cx_type_ok and cx_field == "tune.k", (cx_field, cx_type_ok)
    assert cx_run == (2, 2) and cx_check == (2, 2), (cx_run, cx_check)
    assert not cx_exists, "(cX) the malformed run wrote a directory"
    assert check_all == (0, 0) and family_line == "family: four_port_birdcage", (
        check_all, family_line)
    assert rc_all == (0, 0), rc_all
    assert summary["cells"] == g.ANS6_EXACT_STRIP_HOLE_CELL_RECORD == REC_CELLS
    for c in g.CLASSES:
        assert rel_means[c] <= g.LEG_D0_REPRODUCTION_BAND, f"(R0) {c} {rel_means[c]:.3e}"
    assert same_after_api, "(I=) tune() replaced study.result"
    assert i_mem, "(I=) summary != api tune() bitwise"
    assert i_disk, "(I=) summary.json != in-memory tune block bitwise"
    assert kt == 1 and _hex_eq(t["frequency_hz"], g.FREQUENCY_64_HZ)
    assert abs(gap_a) <= g.REDUCTION_BAND, f"(a) {gap_a:+.6e}"
    assert _within(c_mode, "c_mode") and _within(gap_a, "gap_a"), "(a) record digits"
    assert abs(gap_c1) > g.REDUCTION_BAND, f"(c1) blind {gap_c1:+.6e}"
    assert t2["k"] == 2
    assert _within(t2["c_mode_farad"], "c_k2") and _within(gap_c1, "gap_c1"), "(c1) record"
    assert z_tie <= EXACT_IDENTITY_RTOL, f"(Z) {z_tie:.3e}"
    assert spec0 == {}, spec0
    assert d_equal and t0s["k"] == 1, "(D) omitted k != the k = 1 block"
    assert s4p.name == "f_small_hole_pec_tune.s4p", s4p.name
    assert t_nf == 1 and reading["nports"] == 4, (t_nf, reading["nports"])
    assert t_df <= FREQ_ATOL_HZ, f"(T) frequency {t_df:.3e}"
    assert t_ds <= ANCHOR_TOL, f"(T) max|dS| {t_ds:.3e}"
    assert t_rel <= g.LEG_D0_REPRODUCTION_BAND, f"(T) C rel {t_rel:.3e}"
    assert same_after_reads, "(Once) a read replaced study.result"

    # ---- the mode field (b1_plus is a second, field-keeping build and sweep) --------
    t_f0 = time.perf_counter()
    b1 = run.study.b1_plus(comm=comm)
    t_fields = time.perf_counter() - t_f0
    fres = run.study.fields_result
    port_order = [d.port_id for d in fres.port_definitions]
    assert port_order == [d.port_id for d in run.result.port_definitions]
    assert port_order == sorted(b1, key=lambda p: int(p[1:])), (port_order, sorted(b1))
    s_f = np.asarray(fres.s_matrix, dtype=np.complex128)
    ds_fields = float(np.max(np.abs(s_f - s_run)) / np.max(np.abs(s_run)))

    f_dft = _four_port_dft()
    space = b1[port_order[0]].function_space
    msh = space.mesh

    def _mode(k, name):
        fk = fem.Function(space, name=name)
        acc = np.zeros_like(b1[port_order[0]].x.array)
        for p, pid in enumerate(port_order):
            acc = acc + np.conj(f_dft[k, p]) * b1[pid].x.array
        fk.x.array[:] = acc
        fk.x.scatter_forward()
        return fk

    mode_t = _mode(kt, "B_tuned_mode")
    mode_2 = _mode(2, "B_mode2")

    origin = np.zeros((1, 3))
    centre = {}
    for label, fk in (("tuned", mode_t), ("mode2", mode_2)):
        vals, valid = evaluate_vector_field_parallel(fk, origin)
        vals = np.asarray(vals).reshape(-1, 3)
        assert bool(np.all(np.asarray(valid, dtype=bool))), f"{label}: origin not found"
        v = vals[0]
        centre[label] = {
            "b1p": abs(v[0] + 1j * v[1]) / 2.0,
            "b1m": abs(v[0] - 1j * v[1]) / 2.0,
            "norm": float(np.linalg.norm(v)),
        }

    space_s = fem.functionspace(msh, ("Lagrange", 1))

    def _scalar(name, values):
        f = fem.Function(space_s, name=name)
        assert values.shape[0] == f.x.array.shape[0], (values.shape, f.x.array.shape)
        f.x.array[:] = values
        f.x.scatter_forward()
        return f

    nbt = mode_t.x.array.reshape(-1, 3)
    nb2 = mode_2.x.array.reshape(-1, 3)
    real_t = fem.Function(space, name="B_tuned_mode_real")
    real_t.x.array[:] = np.real(mode_t.x.array)
    real_t.x.scatter_forward()
    fields = {
        "B_tuned_mode_real": real_t,
        "B_tuned_mode_abs": _scalar("B_tuned_mode_abs", np.sqrt(np.sum(np.abs(nbt) ** 2, axis=1))),
        "B1_plus_tuned_mode": _scalar("B1_plus_tuned_mode",
                                      np.abs(nbt[:, 0] + 1j * nbt[:, 1]) / 2.0),
        "B1_minus_tuned_mode": _scalar("B1_minus_tuned_mode",
                                       np.abs(nbt[:, 0] - 1j * nbt[:, 1]) / 2.0),
        "B_mode2_abs": _scalar("B_mode2_abs", np.sqrt(np.sum(np.abs(nb2) ** 2, axis=1))),
    }
    for key, fn in fields.items():
        assert fn.name == key, (key, fn.name)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, fres.cell_tags, fields, comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    if comm.rank == 0:
        print(f"[fields] port order {port_order}; field-keeping sweep {t_fields:.1f} s; "
              f"|S(fields) - S(result)|max/|S|max {ds_fields:.3e} (predicted <= 1e-12; "
              "mri:19 read 1.681e-15)")
        ct, c2 = centre["tuned"], centre["mode2"]
        print(f"[mode k={kt}] centre |B1+| {ct['b1p']:.6g}, |B1-| {ct['b1m']:.6g}, "
              f"|B| {ct['norm']:.6g} T per unit incident-wave mode excitation "
              f"(mri:19 mode 1: {MRI19[1][0]:.6g}, {MRI19[1][1]:.6g}, {MRI19[1][2]:.6g})")
        print(f"[mode k=2] centre |B1+| {c2['b1p']:.6g}, |B1-| {c2['b1m']:.6g}, "
              f"|B| {c2['norm']:.6g} T (mri:19 mode 2: {MRI19[2][0]:.6g}, "
              f"{MRI19[2][1]:.6g}, {MRI19[2][2]:.6g})")
        print(f"[mode] |B_2(0)|/|B_k(0)| {c2['norm'] / ct['norm']:.4e} "
              f"(mri:19: {MRI19['ratio']:.4e})")
        print(f"[paraview] wrote {xdmf_path}")
        print(
            f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
