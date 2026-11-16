"""`WF-4` step 1 (§9 item 471b, chain U42): three presets under ``examples/presets/``.

Brief: docs/planning/chunks/WF-4.md, "§9 item 471" part 471b.  Each preset is a
YAML file a user passes to ``femem run``; each reproduces one gated record at its
imported band through the CLI, and each is priced by ``femem run --estimate``
inside ``BRACKET`` of its own ``summary["elapsed_s"]`` measured in the same window.

Stage knob ``FEM_EM_WF4_STAGE``: unset => W1 (``debug``, ``-n 2``; (K) x3,
(R-debug), (P-debug), (c1), (c3), (c4)); ``dev`` => W2 (``-n 2``; (R-dev),
(P-dev), (c2)); ``bl`` => W3 (``-n 8``, degree 2; (R-bl), (P-bl), (c5), (c6)).

Disclosed (2026-10-15 03:00 review ruling): the two degree-1 ``PRICE_TABLE`` rows
(P1, P3) are the ``debug`` / ``dev`` fixtures themselves, so their x2 checks are
reproductions; only ``benchmark-lite`` tests the fit.  ``femem`` stdout and
``summary.json`` are rank-0 artefacts: read on rank 0 and broadcast.
"""

from __future__ import annotations

import dataclasses
import json
import os
import shutil
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import cli, load_study
from fem_em_solver import study_yaml

from tests.complex_mode import complex_only
from tests.package.test_wf13_cli_estimate import _diff, _parse, _patch_solves
from tests.package.test_wf13_study_estimate import BRACKET
from tests.package.test_wf15_study_lossy_sphere import MAT4_DEG1_RECORD, _study
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf1_study_yaml import (
    BIRDCAGE_FILE,
    SPHERE_FILE,
    WF15_DEG2_RECORD,
    _py_birdcage,
    _spec_diff,
)
from tests.package.test_wf1_yaml_copper import PEC_CONTROL_FILE, _py_copper
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_CLASS_MEANS_P2,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
    _circulant_classes,
    _class_spread,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

REPO = Path(__file__).resolve().parents[2]
PRESETS = REPO / "examples" / "presets"
DEBUG = PRESETS / "debug.yaml"
DEV = PRESETS / "dev.yaml"
BL = PRESETS / "benchmark-lite.yaml"
LOGS = REPO / "logs" / "wf4"  # gitignored (logs/)
STAGES = {None: ("w1", 2), "dev": ("w2", 2), "bl": ("w3", ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH)}
KEY = "64"


def _stage():
    raw = os.environ.get("FEM_EM_WF4_STAGE")
    if raw not in STAGES:
        raise ValueError(f"FEM_EM_WF4_STAGE must be unset, 'dev' or 'bl'; got {raw!r}")
    return STAGES[raw]


def _need(stage):
    name, width = _stage()
    if name != stage:
        pytest.skip(f"stage {stage} not selected (FEM_EM_WF4_STAGE selects {name})")
    assert MPI.COMM_WORLD.size == width, f"stage {stage} runs at -n {width}; got -n {MPI.COMM_WORLD.size}"


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _rc_all(rc, comm):
    return comm.allreduce(int(rc), op=MPI.MIN), comm.allreduce(int(rc), op=MPI.MAX)


def test_wf4_stage_knob():
    """The knob's _enabled/_need pattern: unknown values raise, never skip silently."""
    old = os.environ.get("FEM_EM_WF4_STAGE")
    try:
        os.environ["FEM_EM_WF4_STAGE"] = "nonsense"
        with pytest.raises(ValueError):
            _stage()
    finally:
        if old is None:
            os.environ.pop("FEM_EM_WF4_STAGE", None)
        else:
            os.environ["FEM_EM_WF4_STAGE"] = old
    assert _stage()[0] in ("w1", "w2", "w3")


def _price_then_run(preset, capsys):
    """(P): ``--estimate`` with solves patched to raise, then the real run.

    Returns (estimate, parsed-print verdict, summary) with summary read on
    rank 0 and broadcast.
    """
    from fem_em_solver.study import ResourceEstimate

    comm = MPI.COMM_WORLD
    out = LOGS / preset.stem
    est_dir = LOGS / f"{preset.stem}-estimate"
    if comm.rank == 0:
        for d in (out, est_dir):
            if d.exists():
                shutil.rmtree(d)
    comm.Barrier()
    api = load_study(preset).estimate(comm)
    names = [f.name for f in dataclasses.fields(ResourceEstimate)]
    with pytest.MonkeyPatch.context() as mp:
        _patch_solves(mp)
        capsys.readouterr()
        rc = cli.main(["run", "--estimate", str(preset), "--output-dir", str(est_dir)])
        printed = capsys.readouterr().out
    verdict = None
    if comm.rank == 0:
        parsed = _parse(printed, names)
        complete = set(parsed) == set(names) and all(len(v) == 1 for v in parsed.values())
        verdict = dict(complete=complete, parsed=parsed, differ=_diff(parsed, api) if complete else None)
    verdict = comm.bcast(verdict, root=0)
    exists = comm.allreduce(bool(est_dir.exists()), op=MPI.LOR)
    assert _rc_all(rc, comm) == (0, 0)
    assert verdict["complete"], verdict
    assert verdict["differ"] == set(), verdict["differ"]
    assert not exists

    rc = cli.main(["run", str(preset), "--output-dir", str(out)])
    assert _rc_all(rc, comm) == (0, 0)
    summary = None
    if comm.rank == 0:
        with open(out / "summary.json") as fh:
            summary = json.load(fh)
    summary = comm.bcast(summary, root=0)
    from fem_em_solver.utils.instrumentation import report_peak_rss

    peak = report_peak_rss(comm, label=f"[WF-4] {preset.stem} summed peak ru_maxrss")
    ratio = api.predicted_seconds / summary["elapsed_s"]
    with capsys.disabled():
        _p0(
            f"\n[WF-4] price {preset.stem}: predicted_seconds {api.predicted_seconds:.2f} / "
            f"elapsed_s {summary['elapsed_s']:.2f} = {ratio:.3f} (BRACKET {BRACKET}); "
            f"predicted_peak_gib {api.predicted_peak_gib:.2f} vs measured "
            f"{(peak or 0) / 2**30:.2f} GiB summed (printed, not asserted); rows_used "
            f"{api.rows_used}; cells {api.cells}, unknowns {api.unknowns}, -n {comm.size}"
        )
    assert 1.0 / BRACKET <= ratio <= BRACKET
    return api, verdict, summary


def _s(summary):
    return np.array([[complex(*p) for p in row] for row in summary["s_matrix"]])


# --------------------------------------------------------------------- W1


def test_wf4_k_parse_identity():
    """(K) x3 and (c3): presets equal their sources / Python Studies; check prints family."""
    if _stage()[0] != "w1":
        pytest.skip("(K) runs in W1 only")
    comm = MPI.COMM_WORLD
    d = {
        "debug~sphere_mat4": _spec_diff(load_study(DEBUG), load_study(SPHERE_FILE)),
        "debug~python": _spec_diff(load_study(DEBUG), _study(degree=1)),
        "dev~f_small_four_port": _spec_diff(load_study(DEV), load_study(BIRDCAGE_FILE)),
        "dev~python": _spec_diff(load_study(DEV), _py_birdcage()),
        "bl~pec_control": _spec_diff(load_study(BL), load_study(PEC_CONTROL_FILE)),
        "bl~python(deg1)": _spec_diff(load_study(BL), _py_copper("pec")),
    }
    _p0(f"\n[WF-4] (K) diffs {d}; benchmark-lite degree {load_study(BL).degree}")
    for k in ("debug~sphere_mat4", "debug~python", "dev~f_small_four_port", "dev~python"):
        assert d[k] == [], k
    assert d["bl~pec_control"] == ["degree"]
    assert d["bl~python(deg1)"] == ["degree"]
    assert load_study(BL).degree == 2

    import io
    import contextlib

    for f in (DEBUG, DEV, BL):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = cli.main(["check", str(f)])
        fam = None
        if comm.rank == 0:
            lines = [ln for ln in buf.getvalue().splitlines() if ln.startswith("family: ")]
            fam = (lines, study_yaml._load(f)[1])
        fam = comm.bcast(fam, root=0)
        _p0(f"  (K) check {f.name}: rc {rc}; family lines {fam[0]} vs _load {fam[1]!r}")
        assert _rc_all(rc, comm) == (0, 0)
        assert fam[0] == [f"family: {fam[1]}"]

    # (c3) one ULP up on dev's sigma -> exactly the material sigma field differs.
    tmp = LOGS / f"c3_rank{comm.rank}.yaml"
    tmp.parent.mkdir(parents=True, exist_ok=True)
    up = float(np.nextafter(0.5, 1.0))
    text = DEV.read_text()
    assert text.count("  sigma: 0.5\n") == 1
    tmp.write_text(text.replace("  sigma: 0.5\n", f"  sigma: {up!r}\n"))
    c3 = _spec_diff(load_study(tmp), load_study(DEV))
    _p0(f"  (c3) dev sigma {up!r}: diff {c3}")
    assert c3 == ["material.sigma"]


@complex_only
def test_wf4_debug_reproduces_and_prices(capsys):
    """(R-debug), (P-debug), (c1), (c4)."""
    _need("w1")
    api, verdict, summary = _price_then_run(DEBUG, capsys)
    sar = summary["mean_sar"]["mean_sar_w_per_kg"]
    r = abs(sar - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    c1 = abs(sar / WF15_DEG2_RECORD - 1.0)
    # (c4) the comparator on a one-field-perturbed estimate names exactly that field.
    bumped = dataclasses.replace(api, predicted_seconds=api.predicted_seconds * 1.1)
    c4 = _diff(verdict["parsed"], bumped)
    with capsys.disabled():
        _p0(
            f"[WF-4] (R-debug) SAR {sar!r} rel MAT4_DEG1_RECORD {r:.3e} (L0_BAND {L0_BAND:.0e}); "
            f"(c1) vs WF15_DEG2_RECORD {c1:.4%} ({c1 / L0_BAND:.2e} x band); (c4) {c4}; "
            f"family {summary['family']}, degree {summary['degree']}, cells {summary['cells']}"
        )
    assert summary["degree"] == 1 and summary["ranks"] == 2
    assert r <= L0_BAND
    assert c1 > L0_BAND
    assert c4 == {"predicted_seconds"}


# --------------------------------------------------------------------- W2


@complex_only
def test_wf4_dev_reproduces_and_prices(capsys):
    """(R-dev) with PORT-11's gates, (P-dev), (c2)."""
    _need("w2")
    api, _, summary = _price_then_run(DEV, capsys)
    s = _s(summary)
    z = np.array([[complex(*p) for p in row] for row in summary["z_matrix"]])
    r = _rel(s, _record(1))
    c2 = _rel(s, _record(2))
    recip = float(np.linalg.norm(s - s.T) / np.linalg.norm(s))
    sigma = np.linalg.svd(s, compute_uv=False)
    power = np.sum(np.abs(s) ** 2, axis=0)
    spreads = {n: _class_spread(v) for n, v in _circulant_classes(z).items()}
    with capsys.disabled():
        _p0(
            f"[WF-4] (R-dev) ||S - P1||/||P1|| {r:.3e} (band {LEG_D0_REPRODUCTION_BAND:.0e}); "
            f"(c2) vs _record(2) {c2:.4f}; reciprocity {recip:.3e} (band {RECIPROCITY_BAND:.0e}); "
            f"sigma_max {float(np.max(sigma)):.12f}; max column power {float(np.max(power)):.12f}; "
            "spreads " + ", ".join(f"{k} {v * 100:.4f}%" for k, v in spreads.items())
        )
    assert summary["degree"] == 1 and summary["family"] == "four_port_birdcage"
    assert r <= LEG_D0_REPRODUCTION_BAND
    assert c2 > LEG_D0_REPRODUCTION_BAND
    assert recip <= RECIPROCITY_BAND
    assert float(np.max(sigma)) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for v in power:
        assert v <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for name, v in spreads.items():
        assert v <= ADJACENT_SPREAD_BAND, name


# --------------------------------------------------------------------- W3


@complex_only
def test_wf4_benchmark_lite_reproduces_and_prices(capsys):
    """(R-bl), (P-bl), (c5), (c6) — the first degree-2 solve through YAML."""
    _need("w3")
    comm = MPI.COMM_WORLD
    api, _, summary = _price_then_run(BL, capsys)
    means = _class_means(_s(summary))
    want = ANS6_EXACT_STRIP_CLASS_MEANS_P2[KEY]["pec"]
    rel = {c: abs(means[c] - want[c]) / abs(want[c]) for c in CLASSES}
    p1 = ANS6_EXACT_STRIP_CLASS_MEANS_P1[KEY]["pec"]
    c6 = {c: abs(means[c] - p1[c]) / abs(p1[c]) for c in CLASSES}
    # (c5) benchmark-lite's price against debug's measured time (W1's summary on disk).
    dbg = None
    if comm.rank == 0:
        p = LOGS / "debug" / "summary.json"
        dbg = json.loads(p.read_text())["elapsed_s"] if p.exists() else None
    dbg = comm.bcast(dbg, root=0)
    c5 = api.predicted_seconds / dbg if dbg else None
    with capsys.disabled():
        _p0(
            f"[WF-4] (R-bl) cells {summary['cells']} (record {ANS6_EXACT_STRIP_HOLE_CELL_RECORD}); "
            "class rel vs P2 pec: " + ", ".join(f"{c} {v:.3e}" for c, v in rel.items())
            + f" (band {LEG_D0_REPRODUCTION_BAND:.0e}); (c6) vs P1 pec: "
            + ", ".join(f"{c} {v:.3e} ({v / LEG_D0_REPRODUCTION_BAND:.2e} x band)" for c, v in c6.items())
            + f"; (c5) predicted {api.predicted_seconds:.2f} / debug elapsed {dbg} = {c5}"
        )
    assert summary["degree"] == 2 and summary["family"] == "four_port_birdcage"
    assert summary["cells"] == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    for c in CLASSES:
        assert rel[c] <= LEG_D0_REPRODUCTION_BAND, c
    assert max(c6.values()) > LEG_D0_REPRODUCTION_BAND
    assert dbg is not None, "W1's logs/wf4/debug/summary.json missing; run W1 first"
    assert not (1.0 / BRACKET <= c5 <= BRACKET)
