"""Example (`EX-126`): two ``femem run``s of one study file -- the output bundle and its diff.

`WF-2` (item 346, ``tests/package/test_wf2_bundle_reproducible.py``) gates that the
``femem run`` bundle (the run's JSON summary plus, for a port study, the Touchstone file)
reproduces at fixed width.  This example runs the tracked F-small four-port file
``tests/package/wf1/f_small_four_port.yaml`` in place (exact rectangles, kappa off,
degree 1, 64 MHz, ``-n 2``) twice in one process -- ``run_a`` through ``cli.main``,
``run_b`` through ``run_study_file`` -- and diffs the two bundles key by key at two layers:

* the deterministic layer, bitwise: the summary's keys, the Touchstone header lines and the
  key set (``sha256``, ``cells``, ``degree``, ``ranks`` ...);
* the floating layer, at ``L0_BAND``: S and the parsed Touchstone.

Then it plants a changed degree (c1), a P1<->P2 swap (c2), key edits on the real summary
(c3), a one-ULP sigma config run through a third ``femem run`` (c4) and a 1e-8 resolution
change (c5), and shows which layer catches which.  The diff is written as a CSV beside the
XDMF.

Every record, band and helper is imported from the gate module (the ANS-1 rule); record
digits sit only in the example-local ``REC`` dict with their ``log:line``.  New values print
at <= 6 significant digits; S, Z and the Touchstone data lines are never printed.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:27 -n 2 -t 420

Not claimed: reproduction across widths, machines or images; degree 2; the sphere bundle;
Z reproduction (printed only); the partition in the XDMF as reproducible; the physics or
accuracy of S; any package-side bundle comparator (the comparator stays test-side).
"""

from __future__ import annotations

import copy
import csv
import hashlib
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import load_study  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.instrumentation import report_peak_rss  # noqa: E402

import tests.package.test_wf2_bundle_reproducible as g  # noqa: E402
import tests.package.test_wf16_study_four_port as wf16  # noqa: E402
import tests.validation.test_mat5_em_coupled_sphere as mat5  # noqa: E402
from tests.package.test_wf1_study_yaml import _spec_diff  # noqa: E402
from tests.validation.test_ops61_exact_strip_ports import (  # noqa: E402
    ANS4_EXACT_STRIP_CELL_RECORD,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_27_femem_run_bundle_reproducible"
CONTROL_YAML = Path(__file__).resolve().parent / "27_femem_run_bundle_reproducible_sigma_ulp.yaml"
DIR_A = OUTPUT_DIR / "mri_27_run_a"
DIR_B = OUTPUT_DIR / "mri_27_run_b"
DIR_C = OUTPUT_DIR / "mri_27_run_c"

# Gate-log records (tracked), docs/testing/logs/20261009T053110Z_WF-2-W2-four-port.log.
# Each is (record, relative tolerance 10^-(d-2), d = printed significant digits).
REC = {
    # :3688 sha256 of tests/package/wf1/f_small_four_port.yaml (unchanged since c8d3a532)
    "sha256": "7820a305a7b06c10ba6d5a03df95668bc769001ea631031f846199b2d67fa4f7",
    # :3686 cells (asserted against ANS4_EXACT_STRIP_CELL_RECORD)
    "cells": 117206,
    # :3693 Touchstone header lines (bar generated_utc), one generated_utc line per file
    "header_lines": 5,
    # :3694 (c1) degree plant S miss, d = 3
    "c1": (0.0495, 1.0e-1),
    # :3695 (c2) P1<->P2 swap, d = 4
    "c2": (0.2492, 1.0e-2),
    # :3692 cond(I - S_a), d = 4
    "cond": (19.28, 1.0e-2),
}
# Round-off readings, asserted at their band only, printed beside the gate's lines:
GATE_R, GATE_X, GATE_T, GATE_Z = 5.010e-16, 9.557e-11, 5.521e-14, 5.628e-16  # :3690-3693


def _within(value, key):
    rec, rtol = REC[key]
    return abs(value / rec - 1.0) <= rtol


def _hdr(lines):
    # Copied from tests/package/test_wf2_bundle_reproducible.py:296-297 (a closure there):
    # the header lines, bar the generated_utc line by its exact prefix.
    return [ln for ln in lines if ln[:1] in ("!", "#") and not ln.startswith(g.GENERATED_PREFIX)]


def _n_gen(lines):
    return sum(ln.startswith(g.GENERATED_PREFIX) for ln in lines)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    t_start = time.time()

    # ---- structural asserts (parse only, no mesh) ----------------------------------
    if not np.issubdtype(default_scalar_type, np.complexfloating):
        raise RuntimeError("complex build required (source dolfinx-complex-mode)")
    assert comm.size == g.ANS4_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 record is asserted at -n {g.ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}; "
        f"got -n {comm.size}"
    )
    assert g.L0_BAND is mat5.L0_BAND
    assert g._record is wf16._record
    assert REC["cells"] == ANS4_EXACT_STRIP_CELL_RECORD
    base = load_study(g.BIRDCAGE_FILE)
    assert base.degree == 1 and base.conductor is None
    ctrl = load_study(CONTROL_YAML)
    diff_spec = _spec_diff(ctrl, base)
    assert diff_spec == ["material.sigma"], diff_spec
    assert ctrl.material.sigma == float(np.nextafter(base.material.sigma, 1.0))
    bytes_base, bytes_ctrl = g.BIRDCAGE_FILE.read_bytes(), CONTROL_YAML.read_bytes()
    assert bytes_ctrl != bytes_base
    lines_base = bytes_base.decode().splitlines(keepends=True)
    lines_ctrl = bytes_ctrl.decode().splitlines(keepends=True)
    assert len(lines_base) == len(lines_ctrl)
    differ = [i for i, (x, y) in enumerate(zip(lines_base, lines_ctrl)) if x != y]
    assert len(differ) == 2, differ
    assert lines_ctrl[differ[0]].startswith("# EX-126")
    assert lines_ctrl[differ[1]].strip().startswith("sigma:")
    assert DIR_A.resolve() != DIR_B.resolve() and DIR_B.resolve() != DIR_C.resolve()
    assert DIR_A.resolve() != DIR_C.resolve()
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)

    # ---- (1) run_a through cli.main, run_b through run_study_file -----------------------
    g._fresh_dir(DIR_A, comm)
    g._fresh_dir(DIR_B, comm)
    t0 = time.perf_counter()
    rc = g.cli.main(["run", str(g.BIRDCAGE_FILE), "--output-dir", str(DIR_A)])
    t_a = time.perf_counter() - t0
    rc_a = g._rc_all(rc, comm)
    t0 = time.perf_counter()
    run_b = g.run_study_file(g.BIRDCAGE_FILE, output_dir=DIR_B, comm=comm)
    t_b = time.perf_counter() - t0
    comm.Barrier()
    mt_a = (DIR_A / "summary.json").stat().st_mtime
    mt_b = (DIR_B / "summary.json").stat().st_mtime
    a = g._load(DIR_A / "summary.json")
    b = g._load(DIR_B / "summary.json")

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    geometry = run_b.study.geometry
    msh = run_b.result.mesh
    local_tags = np.unique(run_b.result.cell_tags.values)
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
        run_b.result.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=air_tags,
        translucent_tags=(geometry.phantom_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:27 -- two femem runs of one file: the F-small four-port bundle diffed (WF-2)",
        comm=comm,
    )

    # ---- (D) deterministic layer, bitwise ------------------------------------------------
    digest = hashlib.sha256(g.BIRDCAGE_FILE.read_bytes()).hexdigest()
    d_flags = g._bitwise_flags(a, b)
    # ---- (F) file equals the writer's memory, bitwise (broadcast; every rank compares) ---
    writer = comm.bcast(run_b.summary if comm.rank == 0 else None, root=0)
    f_eq = g._same(b, writer)
    local_off = sorted(
        k for k in set(b) | set(run_b.summary) if not g._same(b.get(k), run_b.summary.get(k))
    )
    offs = comm.gather(local_off, root=0)

    # ---- (R), (X), Z, cond ---------------------------------------------------------------
    s_a, s_b = g._s(a), g._s(b)
    p1 = g._record(1)
    r_ab = g._rel(s_b, s_a)
    x_a, x_b = g._rel(s_a, p1), g._rel(s_b, p1)
    z_r = g._rel(g._z(b), g._z(a))
    cond = float(np.linalg.cond(np.eye(s_a.shape[0]) - s_a))

    # ---- (T) Touchstone ------------------------------------------------------------------
    ta, tb = DIR_A / a["touchstone"], DIR_B / b["touchstone"]
    text_a, text_b = ta.read_text(), tb.read_text()
    lines_a, lines_b = text_a.splitlines(), text_b.splitlines()
    ha, hb = _hdr(lines_a), _hdr(lines_b)
    n_gen = (_n_gen(lines_a), _n_gen(lines_b))
    placeholder = ("PLACEHOLDER" in text_a) or ("PLACEHOLDER" in text_b)
    _, mats_a, _ = g.load_touchstone(ta)
    _, mats_b, _ = g.load_touchstone(tb)
    la = np.asarray(mats_a[0], dtype=np.complex128)
    lb = np.asarray(mats_b[0], dtype=np.complex128)
    t_ab = g._rel(lb, la)
    t_file = g._rel(la, s_a)

    # ---- (c1) degree plant, (c2) swap, (c3) key plants, (c5) resolution -----------------
    plant = copy.deepcopy(b)
    plant["degree"] = 2
    plant["s_matrix"] = [[[float(z.real), float(z.imag)] for z in row] for row in g._record(2)]
    c1_flags = g._bitwise_flags(a, plant)
    c1_miss = g._rel(g._s(plant), s_a)
    perm = [1, 0, 2, 3]
    c2_miss = g._rel(s_b[np.ix_(perm, perm)], s_a)
    swapped = copy.deepcopy(b)
    swapped["s_matrix"] = [
        [[float(z.real), float(z.imag)] for z in row] for row in s_b[np.ix_(perm, perm)]
    ]
    c2_flags = g._bitwise_flags(a, swapped)
    c5_miss = g._rel(s_b * (1.0 + 1e-8), s_a)

    p_cells = copy.deepcopy(b)
    p_cells["cells"] += 1
    p_sha = copy.deepcopy(b)
    p_sha["sha256"] = ("1" if b["sha256"][0] != "1" else "2") + b["sha256"][1:]
    p_ranks = copy.deepcopy(b)
    del p_ranks["ranks"]
    p_elapsed = copy.deepcopy(b)
    p_elapsed["elapsed_s"] = b["elapsed_s"] + 1.0
    c3 = {
        "cells+1": g._bitwise_flags(b, p_cells),
        "sha256 digit": g._bitwise_flags(b, p_sha),
        "ranks deleted": g._bitwise_flags(b, p_ranks),
        "elapsed_s changed": g._bitwise_flags(b, p_elapsed),
    }

    # ---- (c4) one-ULP sigma config: the control file through the CLI --------------------
    g._fresh_dir(DIR_C, comm)
    t0 = time.perf_counter()
    rc = g.cli.main(["run", str(CONTROL_YAML), "--output-dir", str(DIR_C)])
    t_c = time.perf_counter() - t0
    rc_c = g._rc_all(rc, comm)
    mt_c = (DIR_C / "summary.json").stat().st_mtime
    c = g._load(DIR_C / "summary.json")
    c4_flags = g._bitwise_flags(a, c)
    c4_rel = g._rel(g._s(c), s_a)
    tc = DIR_C / c["touchstone"]
    hc = _hdr(tc.read_text().splitlines())

    # ---- prints (new values <= 6 significant digits; no S, Z or data line) --------------
    if comm.rank == 0:
        band = g.L0_BAND
        print("=" * 72)
        print("EX-126 -- two femem runs of one study file: the bundle and its diff (mri:27)")
        print("=" * 72)
        print(f"[(W)] comm.size {comm.size}; parse-only: control differs from base by "
              f"{diff_spec}, file lines {differ}")
        print(f"[runs] run_a (cli.main) {t_a:.1f} s, elapsed_s {a['elapsed_s']:.2f}; run_b "
              f"(run_study_file) {t_b:.1f} s, elapsed_s {b['elapsed_s']:.2f}; cells {a['cells']} "
              f"(W2 :3686: 117206)")
        print(f"[run_c] control {t_c:.1f} s, elapsed_s {c['elapsed_s']:.2f}")
        print(f"  (F) printed: per-rank keys where local StudyRun.summary != file {offs} "
              "(W2 :3687: [[], ['elapsed_s']])")
        print(f"  (D) bitwise flags run_a vs run_b {sorted(d_flags)}; key sets equal "
              f"{set(a) == set(b)}; sha256 {a['sha256']} == recomputed "
              f"{a['sha256'] == digest and b['sha256'] == digest}; record equal "
              f"{a['sha256'] == REC['sha256']} (W2 :3688)")
        print(f"  (F) summary file == writer's summary (float.hex) {f_eq}")
        print(f"  (R) ||S_b - S_a||/||S_a|| {r_ab:.4g} ({r_ab / band:.3g} x L0_BAND; "
              f"W2 :3690 {GATE_R:.3e})")
        print(f"  (X) vs _record(1): run_a {x_a:.4g}, run_b {x_b:.4g} "
              f"({x_a / band:.3g} x, {x_b / band:.3g} x; W2 :3691 {GATE_X:.3e}; "
              f"L0_BAND {band:.0e})")
        print(f"  Z printed: ||Z_b - Z_a||/||Z_a|| {z_r:.4g} ({z_r / band:.3g} x; W2 :3692 "
              f"{GATE_Z:.3e}); cond(I - S_a) {cond:.4g} (W2 :3692 19.28)")
        print(f"  (T) header lines (bar generated_utc) {len(ha)} / {len(hb)}, bitwise equal "
              f"{ha == hb}; generated_utc lines {n_gen}; PLACEHOLDER present {placeholder}; "
              f"parsed rel(b, a) {t_ab:.4g} ({t_ab / band:.3g} x; W2 :3693 {GATE_T:.3e}); "
              f"parsed run_a vs summary S {t_file:.4g} ({t_file / band:.3g} x)")
        for i, ln in enumerate(ha):
            print(f"  (T) run_a header {i}: {ln}")
        for i, ln in enumerate(hb):
            print(f"  (T) run_b header {i}: {ln}")
        print(f"  (c1) degree plant: flags {sorted(c1_flags)}; S miss {c1_miss:.4g} "
              f"({c1_miss / band:.3g} x; W2 :3694 0.0495)")
        print(f"  (c2) P1<->P2 swap: S miss {c2_miss:.4g} ({c2_miss / band:.3g} x; W2 :3695 "
              f"0.2492); bitwise flags on the swapped summary {sorted(c2_flags)} (blind by "
              "construction)")
        print("  (c3) key plants on the real run_b summary: "
              + "; ".join(f"{k} {sorted(v)}" for k, v in c3.items()))
        print(f"  (c4) control sha256 {c['sha256']} differs {c['sha256'] != a['sha256']}; "
              f"flags {sorted(c4_flags)}")
        print(f"  (c4) predicted, printed: S rel(c, a) {c4_rel:.4g} ({c4_rel / band:.3g} x "
              f"L0_BAND; prediction <= 1e-12); control header lines equal to run_a's "
              f"{hc == ha} (predicted True)")
        print(f"  (c5) S_b x (1 + 1e-8): rel {c5_miss:.4g} ({c5_miss / band:.3g} x L0_BAND; "
              "W2 :3696 1.000e-08)")

    # ---- asserts, on every rank, from the same files -------------------------------------
    assert rc_a == (0, 0) and rc_c == (0, 0), (rc_a, rc_c)
    assert mt_a > t_start and mt_b > t_start and mt_c > t_start
    assert d_flags == set(), d_flags
    assert set(a) == set(b)
    assert a["sha256"] == b["sha256"] == digest
    assert a["sha256"] == REC["sha256"], "(D) sha256 != W2 :3688"
    assert a["ranks"] == g.ANS4_EXACT_STRIP_RECORD_RANK_WIDTH and a["degree"] == 1
    assert a["touchstone"] == b["touchstone"] == "f_small_four_port.s4p"
    assert a["cells"] == b["cells"] == ANS4_EXACT_STRIP_CELL_RECORD == REC["cells"]
    assert f_eq, "(F) file != writer's summary"
    assert r_ab <= g.L0_BAND, f"(R) {r_ab:.3e}"
    assert x_a <= g.L0_BAND and x_b <= g.L0_BAND, f"(X) {x_a:.3e} {x_b:.3e}"
    assert ha == hb and len(ha) == len(hb) == REC["header_lines"]
    assert n_gen == (1, 1), n_gen
    assert not placeholder
    assert t_ab <= g.L0_BAND, f"(T) {t_ab:.3e}"
    assert t_file <= g.L0_BAND, f"(T) parsed vs summary {t_file:.3e}"
    assert _within(cond, "cond"), f"cond {cond:.4g}"
    assert c1_flags == {"degree"}, c1_flags
    assert c1_miss > g.L0_BAND and _within(c1_miss, "c1"), f"(c1) {c1_miss:.4g}"
    assert c2_miss > g.L0_BAND and _within(c2_miss, "c2"), f"(c2) {c2_miss:.4g}"
    assert c2_flags == set(), c2_flags
    assert c3["cells+1"] == {"cells"}, c3
    assert c3["sha256 digit"] == {"sha256"}, c3
    assert c3["ranks deleted"] == {g.KEYSET}, c3
    assert c3["elapsed_s changed"] == set(), c3
    assert c["sha256"] != a["sha256"]
    assert c["sha256"] == hashlib.sha256(CONTROL_YAML.read_bytes()).hexdigest()
    assert c4_flags == {"source", "sha256", "touchstone"}, c4_flags
    assert c5_miss > g.L0_BAND, f"(c5) {c5_miss:.3e}"

    # ---- outputs (after every assertion): partition field, combined XDMF, diff CSV ------
    space = fem.functionspace(msh, ("DG", 0))
    owner = fem.Function(space, name="owner_rank")
    n_owned = space.dofmap.index_map.size_local
    owner.x.array[:n_owned] = complex(float(comm.rank), 0.0)
    owner.x.scatter_forward()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, run_b.result.cell_tags,
        {"owner_rank": owner}, comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    peak = report_peak_rss(comm, label="mri:27 peak RSS (summed over ranks)")

    flagged = {}

    def _flag(key, label):
        flagged.setdefault(key, []).append(label)

    _flag("degree", "c1")
    _flag("s_matrix", "c1")
    _flag("s_matrix", "c2")
    _flag("cells", "c3")
    _flag("sha256", "c3")
    _flag("ranks", "c3")
    for k in ("source", "sha256", "touchstone"):
        _flag(k, "c4")
    _flag("s_matrix", "c5")

    def _cmp(key, x, y):
        if key == "elapsed_s":
            return "wall_clock", "n/a"
        if key == "s_matrix":
            return "floating", f"{g._rel(g._s(x), g._s(y)):.4g}"
        if key == "z_matrix":
            return "floating", f"{g._rel(g._z(x), g._z(y)):.4g}"
        return "deterministic", "equal" if g._same(x.get(key), y.get(key)) else "differs"

    if comm.rank == 0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_diff.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["key", "layer", "a_vs_b", "a_vs_c4", "flagged_by"])
            for key in sorted(set(a) | set(c)):
                if key in a and key in b and key in c:
                    layer, ab = _cmp(key, a, b)
                    _, ac = _cmp(key, a, c)
                else:
                    layer, ab, ac = "deterministic", "key missing", "key missing"
                w.writerow([key, layer, ab, ac, "+".join(flagged.get(key, []))])
            hdr_c = set(hc)
            hdr_b = set(hb)
            for ln in ha:
                w.writerow([f"header: {ln}", "deterministic",
                            "equal" if ln in hdr_b else "differs",
                            "equal" if ln in hdr_c else "differs", ""])
        print(f"[paraview] wrote {xdmf_path}; diff table {csv_path}")
        print(f"[memory] summed peak RSS {peak / 2**30:.3g} GiB")
        print(
            f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
