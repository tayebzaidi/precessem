"""WF-2 (§9 item 346): the `femem run` output bundle reproduces at fixed width.

Two runs of the same study file at `-n 2`, complex, degree 1, in one process
(review ruling (3), 2026-10-08 18:00): `run_a` through `cli.main`, `run_b`
through `run_study_file` directly (the function `cli.main` calls) so the
in-memory `StudyRun.summary` is available for (F).  The cross-process half is
carried by the registered records (X).

Labels (docs/planning/chunks/WF-2.md, heading `§9 item 346`):
  (D) deterministic keys bitwise between run_a and run_b; sha256 = recomputed.
  (F) summary.json == StudyRun.summary of the same call, float.hex.
  (R) floating fields run_a vs run_b within L0_BAND.
  (X) each bundle vs MAT4_DEG1_RECORD / _record(1) within L0_BAND.
  (T) Touchstone headers bar `! generated_utc:` bitwise; no PLACEHOLDER;
      parsed matrices within L0_BAND.
  (c1) degree plant, (c2) P1<->P2 swap, (c3) key plants -- asserted.
  (c4) one-ULP sigma config -- sha256/source asserted, float layer printed.
  (c5) 1e-8 resolution plant -- printed only.

No solve output is compared bitwise (-n 2 MUMPS is not bitwise run to run).
Every band and record is imported; none is set here.
"""

from __future__ import annotations

import copy
import hashlib
import json
import time

import numpy as np
from mpi4py import MPI

from fem_em_solver import cli
from fem_em_solver.ports.touchstone import load_touchstone
from fem_em_solver.study_yaml import run_study_file

from tests.complex_mode import complex_only
from tests.package.test_wf15_study_lossy_sphere import MAT4_DEG1_RECORD
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf1_study_yaml import (
    BIRDCAGE_FILE,
    CONTROL_FILE,
    REPO,
    SPHERE_FILE,
    WF15_DEG2_RECORD,
    _fresh_dir,
    _rc_all,
)
from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND

LOGS = REPO / "logs" / "wf2"  # gitignored (logs/)
DETERMINISTIC = ("femem", "source", "sha256", "family", "degree", "frequency_hz", "ranks", "cells")
KEYSET = "<keyset>"
GENERATED_PREFIX = "! generated_utc:"
SAR_KEYS = ("mean_sar_w_per_kg", "dissipated_power_w", "volume_m3")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _same(x, y) -> bool:
    """Bitwise structural equality; floats by float.hex."""
    if type(x) is not type(y):
        return False
    if isinstance(x, float):
        return x.hex() == y.hex()
    if isinstance(x, dict):
        return set(x) == set(y) and all(_same(x[k], y[k]) for k in x)
    if isinstance(x, (list, tuple)):
        return len(x) == len(y) and all(_same(a, b) for a, b in zip(x, y))
    return x == y


def _bitwise_flags(a: dict, b: dict) -> set:
    """The bitwise layer: deterministic keys that differ, plus KEYSET."""
    flags = set()
    if set(a) != set(b):
        flags.add(KEYSET)
    keys = DETERMINISTIC + (("touchstone",) if ("touchstone" in a or "touchstone" in b) else ())
    for k in keys:
        if k in a and k in b and not _same(a[k], b[k]):
            flags.add(k)
    return flags


def _s(summary):
    return np.array([[complex(*p) for p in row] for row in summary["s_matrix"]])


def _z(summary):
    return np.array([[complex(*p) for p in row] for row in summary["z_matrix"]])


def _sar_rel(a, b, key="mean_sar_w_per_kg"):
    return abs(b["mean_sar"][key] - a["mean_sar"][key]) / abs(a["mean_sar"][key])


def _load(path):
    with open(path) as fh:
        return json.load(fh)


def _two_runs(file, comm, t_start):
    """run_a via cli.main, run_b via run_study_file; fresh dirs under logs/wf2."""
    dir_a = LOGS / file.stem / "run_a"
    dir_b = LOGS / file.stem / "run_b"
    assert dir_a.resolve() != dir_b.resolve()
    _fresh_dir(dir_a, comm)
    _fresh_dir(dir_b, comm)
    t0 = time.perf_counter()
    rc = cli.main(["run", str(file), "--output-dir", str(dir_a)])
    t_a = time.perf_counter() - t0
    assert _rc_all(rc, comm) == (0, 0)
    t0 = time.perf_counter()
    run_b = run_study_file(file, output_dir=dir_b, comm=comm)
    t_b = time.perf_counter() - t0
    for d in (dir_a, dir_b):
        mt = (d / "summary.json").stat().st_mtime
        assert mt > t_start, (d, mt, t_start)
    a, b = _load(dir_a / "summary.json"), _load(dir_b / "summary.json")
    _p0(
        f"\n[WF-2] {file.name}: run_a (cli.main) {t_a:.1f} s, elapsed_s {a['elapsed_s']:.2f}; "
        f"run_b (run_study_file) {t_b:.1f} s, elapsed_s {b['elapsed_s']:.2f}; cells {a['cells']}"
    )
    return dir_a, dir_b, a, b, run_b


def _assert_d_and_f(file, a, b, run_b):
    comm = MPI.COMM_WORLD
    digest = hashlib.sha256(file.read_bytes()).hexdigest()
    flags = _bitwise_flags(a, b)
    # (F) The file is rank 0's memory (rank 0 writes).  `elapsed_s` is each
    # rank's own wall clock, so a non-root StudyRun.summary differs from the
    # file there by design -- measured 2026-10-09, WF-2-W1-sphere 052650Z:
    # rank 0 printed (F) True while the other rank's `assert f_eq` failed and
    # the window deadlocked to 124.  Compare the file against the writer's
    # memory, broadcast, on every rank; print the non-root reading.
    writer = comm.bcast(run_b.summary if comm.rank == 0 else None, root=0)
    f_eq = _same(b, writer)
    local_off = sorted(k for k in set(b) | set(run_b.summary)
                       if not _same(b.get(k), run_b.summary.get(k)))
    offs = comm.gather(local_off, root=0)
    _p0(f"  (F) printed: per-rank keys where local StudyRun.summary != file {offs}")
    _p0(
        f"  (D) bitwise flags run_a vs run_b {sorted(flags)}; key sets equal {set(a) == set(b)}; "
        f"sha256 {a['sha256']} == recomputed {a['sha256'] == digest and b['sha256'] == digest}\n"
        f"  (F) summary.json == StudyRun.summary (float.hex) {f_eq}"
    )
    assert flags == set(), flags
    assert set(a) == set(b)
    assert a["sha256"] == digest and b["sha256"] == digest
    assert a["ranks"] == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH and a["degree"] == 1
    assert f_eq


def test_wf2_comparator_flags_planted_changes():
    """(c3): structural key plants on a synthetic summary (solve-free)."""
    base = {
        "femem": 1,
        "source": "/workspace/tests/package/wf1/sphere_mat4.yaml",
        "sha256": "0" * 64,
        "family": "sphere_in_box",
        "degree": 1,
        "frequency_hz": 64.0e6,
        "ranks": 2,
        "cells": 74020,
        "mean_sar": {"mean_sar_w_per_kg": 1.0, "dissipated_power_w": 1.0, "volume_m3": 1.0},
        "elapsed_s": 1.0,
    }
    same = _bitwise_flags(base, copy.deepcopy(base))
    p_cells = copy.deepcopy(base)
    p_cells["cells"] += 1
    p_sha = copy.deepcopy(base)
    p_sha["sha256"] = "1" + base["sha256"][1:]
    p_ranks = copy.deepcopy(base)
    del p_ranks["ranks"]
    p_elapsed = copy.deepcopy(base)
    p_elapsed["elapsed_s"] = 2.0
    f_cells = _bitwise_flags(base, p_cells)
    f_sha = _bitwise_flags(base, p_sha)
    f_ranks = _bitwise_flags(base, p_ranks)
    f_elapsed = _bitwise_flags(base, p_elapsed)
    _p0(
        f"\n[WF-2] (c3) unplanted {sorted(same)}; cells+1 {sorted(f_cells)}; sha256 digit "
        f"{sorted(f_sha)}; ranks deleted {sorted(f_ranks)}; elapsed_s changed {sorted(f_elapsed)}"
    )
    assert same == set()
    assert f_cells == {"cells"}
    assert f_sha == {"sha256"}
    assert f_ranks == {KEYSET}
    assert f_elapsed == set()  # wall clock is excluded by design


@complex_only
def test_wf2_sphere_bundle_reproduces():
    """W1: (D), (F), (R), (X) sphere; (c1) sphere; (c4); (c5) sphere."""
    comm = MPI.COMM_WORLD
    width = ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
    assert comm.size == width, f"record width -n {width}; got -n {comm.size}"
    t_start = time.time()
    _, _, a, b, run_b = _two_runs(SPHERE_FILE, comm, t_start)
    _assert_d_and_f(SPHERE_FILE, a, b, run_b)

    # (R)
    r = {k: _sar_rel(a, b, k) for k in SAR_KEYS}
    vol_bitwise = _same(a["mean_sar"]["volume_m3"], b["mean_sar"]["volume_m3"])
    # (X)
    xa = abs(a["mean_sar"]["mean_sar_w_per_kg"] - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    xb = abs(b["mean_sar"]["mean_sar_w_per_kg"] - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    _p0(
        f"  SAR_a {a['mean_sar']['mean_sar_w_per_kg']!r}, SAR_b {b['mean_sar']['mean_sar_w_per_kg']!r}\n"
        "  (R) rel(b, a): "
        + ", ".join(f"{k} {v:.3e} ({v / L0_BAND:.2e} x L0_BAND)" for k, v in r.items())
        + f"; volume_m3 bitwise equal {vol_bitwise} (predicted True)\n"
        f"  (X) vs MAT4_DEG1_RECORD {MAT4_DEG1_RECORD!r}: run_a {xa:.3e} ({xa / L0_BAND:.2e} x), "
        f"run_b {xb:.3e} ({xb / L0_BAND:.2e} x) (L0_BAND {L0_BAND:.0e})"
    )
    for v in r.values():
        assert v <= L0_BAND
    assert xa <= L0_BAND and xb <= L0_BAND

    # (c1) degree plant
    plant = copy.deepcopy(b)
    plant["degree"] = 2
    plant["mean_sar"]["mean_sar_w_per_kg"] = WF15_DEG2_RECORD
    c1_flags = _bitwise_flags(a, plant)
    c1_miss = _sar_rel(a, plant)
    _p0(f"  (c1) degree plant: flags {sorted(c1_flags)}; SAR miss {c1_miss:.4%} ({c1_miss / L0_BAND:.2e} x)")
    assert c1_flags == {"degree"}
    assert c1_miss > L0_BAND

    # (c5) resolution plant, printed only
    c5 = copy.deepcopy(b)
    c5["mean_sar"]["mean_sar_w_per_kg"] *= 1.0 + 1e-8
    c5_r = _sar_rel(a, c5)
    _p0(f"  (c5) SAR_b x (1 + 1e-8): rel {c5_r:.3e} ({c5_r / L0_BAND:.2f} x L0_BAND) (predicted ~1e-8)")

    # (c4) one-ULP sigma config: the control file through the CLI
    out_c = LOGS / CONTROL_FILE.stem / "run_c"
    _fresh_dir(out_c, comm)
    t0 = time.perf_counter()
    rc = cli.main(["run", str(CONTROL_FILE), "--output-dir", str(out_c)])
    t_c = time.perf_counter() - t0
    assert _rc_all(rc, comm) == (0, 0)
    assert (out_c / "summary.json").stat().st_mtime > t_start
    c = _load(out_c / "summary.json")
    c4_flags = _bitwise_flags(a, c)
    c4_sar = _sar_rel(a, c)
    _p0(
        f"  (c4) control run_c {t_c:.1f} s (elapsed_s {c['elapsed_s']:.2f}): sha256 {c['sha256']} "
        f"differs {c['sha256'] != a['sha256']}; flags {sorted(c4_flags)}\n"
        f"  (c4) predicted, printed: SAR rel(c, a) {c4_sar:.3e} ({c4_sar / L0_BAND:.2e} x L0_BAND; "
        f"prediction <= 1e-10)"
    )
    assert c["sha256"] != a["sha256"]
    assert c4_flags == {"source", "sha256"}


@complex_only
def test_wf2_four_port_bundle_reproduces():
    """W2: (D), (F), (R), (X), (T) four-port; (c1), (c2); Z and (c5) printed."""
    comm = MPI.COMM_WORLD
    width = ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
    assert comm.size == width, f"record width -n {width}; got -n {comm.size}"
    t_start = time.time()
    dir_a, dir_b, a, b, run_b = _two_runs(BIRDCAGE_FILE, comm, t_start)
    _assert_d_and_f(BIRDCAGE_FILE, a, b, run_b)
    assert a["touchstone"] == b["touchstone"] == f"{BIRDCAGE_FILE.stem}.s4p"

    s_a, s_b = _s(a), _s(b)
    p1 = _record(1)
    # (R)
    r = _rel(s_b, s_a)
    # (X)
    xa, xb = _rel(s_a, p1), _rel(s_b, p1)
    # Z printed beside cond(I - S_a) (ruling (4))
    z_r = _rel(_z(b), _z(a))
    cond = float(np.linalg.cond(np.eye(s_a.shape[0]) - s_a))
    _p0(
        f"  (R) ||S_b - S_a||/||S_a|| {r:.3e} ({r / L0_BAND:.2e} x L0_BAND)\n"
        f"  (X) vs _record(1): run_a {xa:.3e} ({xa / L0_BAND:.2e} x), run_b {xb:.3e} "
        f"({xb / L0_BAND:.2e} x) (L0_BAND {L0_BAND:.0e})\n"
        f"  Z printed: ||Z_b - Z_a||/||Z_a|| {z_r:.3e} ({z_r / L0_BAND:.2e} x); cond(I - S_a) {cond:.3e}"
    )
    assert r <= L0_BAND
    assert xa <= L0_BAND and xb <= L0_BAND

    # (T)
    ta, tb = dir_a / a["touchstone"], dir_b / b["touchstone"]
    lines_a, lines_b = ta.read_text().splitlines(), tb.read_text().splitlines()

    def _hdr(lines):
        return [ln for ln in lines if ln[:1] in ("!", "#") and not ln.startswith(GENERATED_PREFIX)]

    ha, hb = _hdr(lines_a), _hdr(lines_b)
    n_gen = (
        sum(ln.startswith(GENERATED_PREFIX) for ln in lines_a),
        sum(ln.startswith(GENERATED_PREFIX) for ln in lines_b),
    )
    placeholder = ("PLACEHOLDER" in ta.read_text()) or ("PLACEHOLDER" in tb.read_text())
    _, mats_a, _ = load_touchstone(ta)
    _, mats_b, _ = load_touchstone(tb)
    la = np.asarray(mats_a[0], dtype=np.complex128)
    lb = np.asarray(mats_b[0], dtype=np.complex128)
    t_r = _rel(lb, la)
    _p0(
        f"  (T) header lines (bar generated_utc) {len(ha)} / {len(hb)}, bitwise equal {ha == hb}; "
        f"generated_utc lines {n_gen}; PLACEHOLDER present {placeholder}; parsed rel(b, a) "
        f"{t_r:.3e} ({t_r / L0_BAND:.2e} x)"
    )
    assert ha == hb and len(ha) > 0
    assert n_gen == (1, 1)
    assert not placeholder
    assert t_r <= L0_BAND

    # (c1) degree plant
    plant = copy.deepcopy(b)
    plant["degree"] = 2
    plant["s_matrix"] = [[[float(z.real), float(z.imag)] for z in row] for row in _record(2)]
    c1_flags = _bitwise_flags(a, plant)
    c1_miss = _rel(_s(plant), s_a)
    # (c2) P1<->P2 swap
    perm = [1, 0, 2, 3]
    c2 = _rel(s_b[np.ix_(perm, perm)], s_a)
    # (c5) printed only
    c5 = _rel(s_b * (1.0 + 1e-8), s_a)
    _p0(
        f"  (c1) degree plant: flags {sorted(c1_flags)}; S miss {c1_miss:.4f} ({c1_miss / L0_BAND:.2e} x)\n"
        f"  (c2) P1<->P2 swap: {c2:.4f} ({c2 / L0_BAND:.2e} x)\n"
        f"  (c5) S_b x (1 + 1e-8): rel {c5:.3e} ({c5 / L0_BAND:.2f} x L0_BAND) (predicted ~1e-8)"
    )
    assert c1_flags == {"degree"}
    assert c1_miss > L0_BAND
    assert c2 > L0_BAND
