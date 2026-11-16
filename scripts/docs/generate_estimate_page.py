"""Generate docs/guide/estimates.md: ``femem run --estimate`` over every shipped study file.

`WF-17` step 30 (§9 item 827, §10 U130).  Brief: docs/planning/chunks/WF-17.md,
"§9 item 827".  Run in-process under ``mpiexec -n 2`` (``--estimate`` prices the
current launch width, ``cli.py``); never shells out to ``femem``, so the build
memo below reaches every call.  Each entry is run twice -- once through
``cli.main(["run", "--estimate", f])`` with stdout / stderr captured, once through
``load_study(f).estimate(comm)`` -- and every geometry's ``build`` is memoised
on ``(repr(self), comm.size)``, so there is one real mesh build per geometry.

Usage::

    mpiexec -n 2 python3 scripts/docs/generate_estimate_page.py --out docs/guide/estimates.md [--check]

Every digit on the page is generated (item 718), never typed.

Step 47 (§9 item 961, §10 U157) appends one section: the measured 3 T box-window
prices of `xl` entries 58 / 60, read from ``docs/testing/xl-ledger.md`` by
string, beside ``Study.estimate``'s FC3 row (``tests/package/wf13/f_clinical_3t.yaml``,
no mesh), with the estimator / measured ratios.  ``render(rows)`` is unchanged.
Step 51 (§9 item 977, §10 U161) adds entry 61 and a per-row ``tag 403`` column.
Step 52 (§9 item 980, §10 U162) inserts a 1.5 T section before the 3 T one: entries
48 / 52 / 59 beside FC1 (``tests/package/wf13/f_clinical_15t.yaml``), through the
same builders with default-preserving keywords; ``render(rows, measured)`` is unchanged.
Step 53 (§9 item 992, §10 U163) appends a degree-2 clinical section after the 3 T one:
`xxl` entry 62's window and its two calls, read from ``docs/testing/xxl-ledger.md`` by
string, beside FC2 (``tests/package/wf13/f_clinical_15t_d2.yaml``), through its own
renderer ``render_measured_d2``; ``render(rows, measured, measured_15t)`` is unchanged.
Step 54 (§9 item 1001, §10 U166) inserts a 3 T degree-2 section between the 3 T and the
degree-2 sections: `xl` entry 64's window from ``docs/testing/xl-ledger.md`` by string,
beside FC4 (``tests/package/wf13/f_clinical_3t_d2.yaml``, not priced: ``nan`` ratios) and
FC3 (degree 1, same mesh and width), through its own renderer ``render_measured_3t_d2``;
every landed constant, default and renderer is unchanged.
"""

from __future__ import annotations

import argparse
import contextlib
import dataclasses
import io
import os
import re
import sys
from pathlib import Path

import yaml
from mpi4py import MPI

REPO = Path(__file__).resolve().parents[2]
DEFAULT_OUT = REPO / "docs" / "guide" / "estimates.md"
SCRATCH = REPO / "logs" / "wf17" / "s30"
LAUNCH_RANKS = 2

# The 14 entries, pinned literally (brief "The file set").
FILES = (
    "examples/mri/09_study_from_yaml_four_port.yaml",
    "examples/mri/09_study_from_yaml_sphere.yaml",
    "examples/mri/13_temperature_from_yaml_sphere.yaml",
    "examples/mri/13_temperature_from_yaml_sphere_control_2g.yaml",
    "examples/mri/27_femem_run_bundle_reproducible_sigma_ulp.yaml",
    "examples/presets/dev.yaml",
    "examples/presets/debug.yaml",
    "examples/presets/benchmark-lite.yaml",
)
FENCES = (
    "capacitor-esr-yaml",
    "implant-rod-file",
    "f-small-hole-pec-tune",
    "f-human-tuned-copper",
    "swept-tune-file",
    "loaded-yaml-file",
)
TREES = ("examples", "docs/guide")
EXCLUDED_PART = "paraview_output"  # gitignored run copies
FENCE_RE = re.compile(r"^```yaml (\S+)$")
TIERS = (("smoke", 30.0), ("standard", 180.0), ("heavy", 1200.0))
MAX_RANKS = 12

# `WF-17` step 47 (§9 item 961, §10 U157): the measured 3 T box-window prices of
# `xl` entries 58 / 60, read from the ledger by string, beside the FC3 row.
LEDGER = REPO / "docs" / "testing" / "xl-ledger.md"
FC3T_FILE = REPO / "tests" / "package" / "wf13" / "f_clinical_3t.yaml"
MEASURED_3T = (("58", "2026-11-09", "GEO-38-C6-3T"), ("60", "2026-11-11", "GEO-38-P33-3T"),
               ("61", "2026-11-12", "GEO-36-step11-I36-3T"))
# Step 51 (§9 item 977, §10 U161): what tag 403 was in each window, restating its
# ledger cell in words (ASCII only: the page has no non-ASCII character).
TAG403 = {"58": "gel", "60": "gel", "61": "sigma_wire 800 S/m (arm M); gel (arm R)"}
MEASURED_HEADING = "## Measured 3 T box-window prices beside the estimator"
# Step 52 (§9 item 980, §10 U162): the green 1.5 T rod-in-box windows beside FC1.
FC15_FILE = REPO / "tests" / "package" / "wf13" / "f_clinical_15t.yaml"
MEASURED_15T = (("48", "2026-11-02", "GEO-36-step9-I36"), ("52", "2026-11-03", "GEO-38-C6"),
                ("59", "2026-11-10", "GEO-38-P33"))
# §9 item 982: the 1.5 T tags in their own dict, so TAG403 (the 3 T module's
# (Tag) premise) stays main's (03:00 11-14 ruling, option (i)).
TAG403_15T = {"48": "sigma_wire 800 S/m (arm M); gel (arm R)", "52": "gel", "59": "gel"}
MEASURED_15T_HEADING = "## Measured 1.5 T box-window prices beside the estimator"
FC15_MESH = "the F-clinical-1.5T r1 mesh (no rod)"
# Step 53 (§9 item 992, §10 U163): `xxl` entry 62 (degree 2, -n 16) beside FC2.
XXL_LEDGER = REPO / "docs" / "testing" / "xxl-ledger.md"
FC2_FILE = REPO / "tests" / "package" / "wf13" / "f_clinical_15t_d2.yaml"
MEASURED_D2 = (("62", "2026-11-14", "GEO-38-C19"),)
MEASURED_D2_HEADING = "## Measured degree-2 clinical prices beside the estimator"
# Its own dict (03:00 11-14 ruling, option (i)); restates log :12412 / :12833 in words.
TAG403_D2 = {"62 call 1": "gel", "62 call 2": "sigma_wire (arm M); gel (arm R)"}
CALL_LABELS_D2 = {"62 call 1": "C6 degree 2 + grid", "62 call 2": "I36 degree 2"}
# Step 54 (§9 item 1001, §10 U166): `xl` entry 64, the 3 T degree-2 root grid, own constants.
FC4_FILE = REPO / "tests" / "package" / "wf13" / "f_clinical_3t_d2.yaml"
MEASURED_3T_D2 = (("64", "2026-11-15", "GEO-38-C5f"),)
MEASURED_3T_D2_HEADING = "## Measured 3 T degree-2 price beside the estimator"
TAG403_3T_D2 = {"64": "absent"}  # `20261115T080008Z_GEO-38-C5f.log:13016`, (G) tag 403 absent
CALL_WINDOW_RE = re.compile(r"window (\d[\d ]*\.\d) s \(\x60:(\d+)\x60\)")

HEADER = (
    "<!-- GENERATED by scripts/docs/generate_estimate_page.py (`WF-17` step 30); "
    "do not edit by hand. -->"
)


class GenerationError(RuntimeError):
    pass


# ------------------------------------------------------------- discovery
def _is_study(text: str) -> bool:
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError:
        return False
    return isinstance(data, dict) and "femem" in data


def fence_blocks(root: Path = REPO) -> dict:
    """tag -> (repo-relative page, body) for every ```yaml <tag> fence under TREES."""
    out = {}
    for tree in TREES:
        for page in sorted((root / tree).rglob("*.md")):
            rel = page.relative_to(root).as_posix()
            if EXCLUDED_PART in rel or page.name.endswith("_private.md"):
                continue
            lines = page.read_text(encoding="utf-8").splitlines()
            for i, line in enumerate(lines):
                m = FENCE_RE.match(line)
                if not m:
                    continue
                body = []
                for ln in lines[i + 1 :]:
                    if ln == "```":
                        break
                    body.append(ln)
                text = "\n".join(body) + "\n"
                if _is_study(text):
                    if m.group(1) in out:
                        raise GenerationError(f"fence tag {m.group(1)!r} appears twice")
                    out[m.group(1)] = (rel, text)
    return out


def discover(root: Path = REPO) -> tuple[list, list]:
    """(study YAML files, study fence tags) found by globbing TREES."""
    files = []
    for tree in TREES:
        for p in sorted((root / tree).rglob("*.yaml")):
            rel = p.relative_to(root).as_posix()
            if EXCLUDED_PART in rel:
                continue
            if _is_study(p.read_text(encoding="utf-8")):
                files.append(rel)
    return sorted(files), sorted(fence_blocks(root))


# ------------------------------------------------------------- tier rule
def tier_for(est) -> str:
    """`xl` when ranks > 12 or seconds > 1200; else the smallest covering tier."""
    if est.ranks > MAX_RANKS or est.predicted_seconds > TIERS[-1][1]:
        return "xl"
    for name, ceiling in TIERS:
        if est.predicted_seconds <= ceiling:
            return name
    raise AssertionError("unreachable")


# ------------------------------------------------------------- build memo
@contextlib.contextmanager
def memoised_builds(comm=MPI.COMM_WORLD):
    """Wrap every geometry ``build`` with a ``(repr(self), comm.size)`` memo."""
    from fem_em_solver import study as st

    classes = [
        c for c in vars(st).values()
        if isinstance(c, type) and c.__module__ == st.__name__ and "build" in vars(c)
    ]
    memo, count = {}, {"builds": 0, "hits": 0, "keys": []}
    saved = {}
    for cls in classes:
        orig = vars(cls)["build"]
        saved[cls] = orig

        def _memo(self, c=MPI.COMM_WORLD, _orig=orig):
            key = (repr(self), c.size)
            if key in memo:
                count["hits"] += 1
                return memo[key]
            count["builds"] += 1
            count["keys"].append(type(self).__name__)
            memo[key] = _orig(self, c)
            return memo[key]

        cls.build = _memo
    try:
        yield count
    finally:
        for cls, orig in saved.items():
            cls.build = orig


# ------------------------------------------------------------- rows
def entries(root: Path = REPO, comm=MPI.COMM_WORLD) -> list:
    """(label, path) for the 14 entries; fences written to SCRATCH by rank 0."""
    blocks = fence_blocks(root)
    missing = [t for t in FENCES if t not in blocks]
    if missing:
        raise GenerationError(f"pinned fence(s) not found: {missing}")
    out = [(f, root / f) for f in FILES]
    if comm.rank == 0:
        SCRATCH.mkdir(parents=True, exist_ok=True)
        for tag in FENCES:
            (SCRATCH / f"{tag}.yaml").write_text(blocks[tag][1], encoding="utf-8")
    comm.Barrier()
    out += [(f"{blocks[t][0]}, fence {t}", SCRATCH / f"{t}.yaml") for t in FENCES]
    return out


def row_for(label: str, path: Path, comm=MPI.COMM_WORLD) -> dict:
    """One entry through the CLI and the API, same process, same comm."""
    from fem_em_solver.cli import main as cli_main
    from fem_em_solver.study import ResourceEstimate
    from fem_em_solver.study_yaml import _load, load_study

    so, se = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(so), contextlib.redirect_stderr(se):
        rc = cli_main(["run", "--estimate", str(path)])
    est, api_error = None, None
    try:
        est = load_study(path).estimate(comm)
    except (ValueError, TypeError) as exc:
        api_error = str(exc)
    kind = _load(path)[1]
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    names = [f.name for f in dataclasses.fields(ResourceEstimate)]
    printed = {}
    for ln in so.getvalue().splitlines():
        key, sep, val = ln.partition(": ")
        if sep and key in names:
            printed.setdefault(key, []).append(val)
    return dict(
        label=label, path=str(path), rc=rc, stdout=so.getvalue(), stderr=se.getvalue(),
        printed=printed, est=est, api_error=api_error, family=kind,
        load_result="load_result" in data,
        thermal_unpriced=("thermal" in data),
    )


def collect(comm=MPI.COMM_WORLD, root: Path = REPO):
    with memoised_builds(comm) as count:
        rows = [row_for(label, path, comm) for label, path in entries(root, comm)]
    return rows, dict(count)


# ------------------------------------------------------------- page
def _cell(text: str) -> str:
    return text.replace("|", "\\|")


COLUMNS = ("degree", "cells", "unknowns", "ranks", "predicted_peak_gib",
           "predicted_seconds", "rows_used", "not_priced")


LOG_CELL_RE = re.compile(r"^`([^`]+\.log)`$")
RANKS_CELL_RE = re.compile(r"^\d+$")
CELLS_CELL_RE = re.compile(r"^([\d ]+) \(")
GIB_CELL_RE = re.compile(r"^\*\*([\d.]+)\*\*")
SECONDS_CELL_RE = re.compile(r"^\d+$")


def ledger_prices(ledger: Path = LEDGER, entries=MEASURED_3T) -> tuple:
    """One dict of strings per ``entries`` item, from its one ledger row (shape-checked)."""
    lines = Path(ledger).read_text(encoding="utf-8").splitlines()
    out = []
    for entry, date, chunk in entries:
        hits = [ln for ln in lines if ln.startswith(f"| {date} | {chunk} | ")]
        if len(hits) != 1:
            raise GenerationError(f"entry {entry}: {len(hits)} ledger rows for {date} {chunk}")
        cells = hits[0].split(" | ")
        if len(cells) < 7:
            raise GenerationError(f"entry {entry}: {len(cells)} cells in the ledger row")
        m_log = LOG_CELL_RE.match(cells[2])
        m_ranks = RANKS_CELL_RE.match(cells[3])
        m_cells = CELLS_CELL_RE.match(cells[4])
        m_gib = GIB_CELL_RE.match(cells[5])
        m_sec = SECONDS_CELL_RE.match(cells[6])
        bad = [name for name, m in (("log", m_log), ("ranks", m_ranks), ("cells", m_cells),
                                    ("peak GiB", m_gib), ("seconds", m_sec)) if m is None]
        if bad:
            raise GenerationError(f"entry {entry}: ledger cell shape failed: {bad}")
        out.append(dict(entry=entry, date=date, chunk=chunk, log=m_log.group(1),
                        ranks=m_ranks.group(0), cells=m_cells.group(1),
                        gib=m_gib.group(1), seconds=m_sec.group(0)))
    return tuple(out)


def measured_3t(comm=MPI.COMM_WORLD, ledger: Path = LEDGER, fc_file: Path = FC3T_FILE,
                entries=MEASURED_3T) -> dict:
    """Ledger strings, the estimator row (CLI / API, no mesh) and the estimator / measured ratios."""
    prices = ledger_prices(ledger, entries)
    row = row_for(Path(fc_file).relative_to(REPO).as_posix(), Path(fc_file), comm)
    # The CLI prints on rank 0 only: every rank renders rank 0's strings.
    printed, stdout = comm.bcast((row["printed"], row["stdout"]) if comm.rank == 0 else None,
                                 root=0)
    row = dict(row, printed=printed, stdout=stdout)
    est = row["est"]
    ratios = []
    for p in prices:
        ratios.append(dict(
            entry=p["entry"],
            gib=f"{est.predicted_peak_gib / float(p['gib']):.4f}",
            seconds=f"{est.predicted_seconds / float(p['seconds']):.4f}",
        ))
    return dict(prices=prices, row=row, ratios=tuple(ratios))


def render_measured(measured, heading: str = MEASURED_HEADING, field: str = "3 T",
                    fc_mesh: str = "the stand-in-cylinder twin", tags=TAG403) -> str:
    """One measured-price section from its heading (every figure from a variable)."""
    prices, row, ratios = measured["prices"], measured["row"], measured["ratios"]
    pr = row["printed"]
    fc_name = ", ".join(row["est"].rows_used)
    ids = [p["entry"] for p in prices]
    entries_txt = ids[0] if len(ids) == 1 else f"{', '.join(ids[:-1])} and {ids[-1]}"
    widths = sorted({p["ranks"] for p in prices})
    lines = [
        heading,
        "",
        f"The estimator's {field} row ({fc_name}, `ranks` {pr['ranks'][0]}) was measured on "
        f"{fc_mesh} and its seconds are one EM solve with the mesh build excluded; "
        f"`xl` entries {entries_txt} ran the rod-in-box {field} mesh (tag 403 as each row states) at "
        f"`-n {', '.join(widths)}`, their seconds are whole-window wall times (build, EM solve "
        "and readouts) and their GiB the cgroup `memory.peak`. The ratios therefore set two "
        "different quantities, on two meshes, at two widths, beside each other. Measured "
        "values are the strings of `docs/testing/xl-ledger.md`; the estimator values are the "
        f"CLI's printed `repr` of `femem run --estimate {row['label']}`.",
        "",
        "| Row | Log | ranks | cells | tag 403 | peak GiB | seconds | ratio GiB | ratio seconds |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for p, r in zip(prices, ratios):
        lines.append(
            f"| `xl` entry {p['entry']} (`{p['chunk']}`) | `{p['log']}` | {p['ranks']} | "
            f"{p['cells']} | {_cell(tags[p['entry']])} | {p['gib']} | {p['seconds']} | "
            f"{r['gib']} | {r['seconds']} |"
        )
    lines.append(
        f"| `Study.estimate` {fc_name} | `{_cell(row['label'])}` | `{_cell(pr['ranks'][0])}` | "
        f"`{_cell(pr['cells'][0])}` | -- | `{_cell(pr['predicted_peak_gib'][0])}` | "
        f"`{_cell(pr['predicted_seconds'][0])}` | -- | -- |"
    )
    lines += [
        "",
        "Ratio = estimator / measured. A record, not a calibration; no row is fitted.",
        "",
    ]
    return "\n".join(line.rstrip() for line in lines)


def ledger_calls(ledger: Path = XXL_LEDGER, entry=MEASURED_D2[0]) -> tuple:
    """Entry 62's two per-call windows (strings) from its one ledger row's Price clause."""
    eid, date, chunk = entry
    lines = Path(ledger).read_text(encoding="utf-8").splitlines()
    hits = [ln for ln in lines if ln.startswith(f"| {date} | {chunk} | ")]
    if len(hits) != 1:
        raise GenerationError(f"entry {eid}: {len(hits)} ledger rows for {date} {chunk}")
    found = CALL_WINDOW_RE.findall(hits[0])
    if len(found) != 2:
        raise GenerationError(f"entry {eid}: {len(found)} per-call windows in the ledger row")
    return tuple(dict(key=f"{eid} call {i}", seconds=s, line=ln)
                 for i, (s, ln) in enumerate(found, start=1))


def measured_d2(comm=MPI.COMM_WORLD, ledger: Path = XXL_LEDGER) -> dict:
    """``measured_3t`` on entry 62 / FC2, plus the per-call windows and their seconds ratios."""
    m = measured_3t(comm, ledger=ledger, fc_file=FC2_FILE, entries=MEASURED_D2)
    calls = ledger_calls(ledger)
    est = m["row"]["est"]
    call_ratios = tuple(
        dict(key=c["key"],
             seconds=f"{est.predicted_seconds / float(c['seconds'].replace(' ', '')):.4f}")
        for c in calls
    )
    return dict(m, calls=calls, call_ratios=call_ratios)


def render_measured_d2(measured) -> str:
    """The degree-2 clinical section (`xxl`, one width); ``render_measured`` is not reused."""
    prices, row, ratios = measured["prices"], measured["row"], measured["ratios"]
    calls, call_ratios = measured["calls"], measured["call_ratios"]
    pr = row["printed"]
    fc_name = ", ".join(row["est"].rows_used)
    p, r = prices[0], ratios[0]
    lines = [
        MEASURED_D2_HEADING,
        "",
        f"The estimator's degree-2 1.5 T row ({fc_name}, `ranks` {pr['ranks'][0]}) was measured "
        f"on {FC15_MESH} and its seconds are one EM solve with the mesh build excluded; "
        f"`xxl` entry {p['entry']} ran the rod-in-box 1.5 T mesh at degree 2 at `-n {p['ranks']}` "
        "in two calls in one window -- call 1 C6 with the box-wide grid, call 2 I36 in two arms "
        "(tag 403 as each row states). Each call's seconds are that call's window (build, EM "
        "solves and readouts), the window's seconds the harness wall time, and the GiB the "
        "cgroup `memory.peak` over the whole window, so the call rows carry no GiB and no GiB "
        "ratio. The ratios therefore set two different quantities, on two meshes, at one "
        "width, beside each other. Measured values are the strings of "
        "`docs/testing/xxl-ledger.md`; the estimator values are the CLI's printed `repr` of "
        f"`femem run --estimate {row['label']}`.",
        "",
        "| Row | Log | ranks | cells | tag 403 | peak GiB | seconds | ratio GiB | ratio seconds |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for c, cr in zip(calls, call_ratios):
        n = c["key"].rsplit(" ", 1)[1]
        lines.append(
            f"| `xxl` entry {p['entry']} call {n} (`{p['chunk']}`, {CALL_LABELS_D2[c['key']]}) | "
            f"`{p['log']}` | {p['ranks']} | {p['cells']} | {_cell(TAG403_D2[c['key']])} | -- | "
            f"{c['seconds']} | -- | {cr['seconds']} |"
        )
    lines.append(
        f"| `xxl` entry {p['entry']} window (`{p['chunk']}`) | `{p['log']}` | {p['ranks']} | "
        f"{p['cells']} | as calls 1 and 2 | {p['gib']} | {p['seconds']} | {r['gib']} | "
        f"{r['seconds']} |"
    )
    lines.append(
        f"| `Study.estimate` {fc_name} | `{_cell(row['label'])}` | `{_cell(pr['ranks'][0])}` | "
        f"`{_cell(pr['cells'][0])}` | -- | `{_cell(pr['predicted_peak_gib'][0])}` | "
        f"`{_cell(pr['predicted_seconds'][0])}` | -- | -- |"
    )
    lines += [
        "",
        "Ratio = estimator / measured. A record, not a calibration; no row is fitted. The "
        "window was lightly contended (an operator `GEO-34` window of 104 s ran inside it, "
        "`xxl-ledger.md:27`).",
        "",
    ]
    return "\n".join(line.rstrip() for line in lines)


def measured_3t_d2(comm=MPI.COMM_WORLD, ledger: Path = LEDGER) -> dict:
    """``measured_3t`` on entry 64 twice: FC4 (degree 2, not priced) and FC3 (degree 1)."""
    return dict(fc4=measured_3t(comm, ledger, FC4_FILE, MEASURED_3T_D2),
                fc3=measured_3t(comm, ledger, FC3T_FILE, MEASURED_3T_D2))


def render_measured_3t_d2(measured) -> str:
    """The 3 T degree-2 section (`xl` entry 64); ``render_measured`` is not reused."""
    fc4, fc3 = measured["fc4"], measured["fc3"]
    p = fc4["prices"][0]
    lines = [
        MEASURED_3T_D2_HEADING,
        "",
        f"`xl` entry {p['entry']} ran the F-clinical-3T r1 coil mesh (no box, tag 403 "
        f"{TAG403_3T_D2[p['entry']]}) at degree 2 at `-n {p['ranks']}`; its seconds are the "
        "whole-window wall time (build, the Z50 pair and three node pairs) and its GiB the "
        "cgroup `memory.peak`. The estimator's 3 T degree-2 row FC4 is not priced "
        "(`Study.estimate` prints its `not_priced` `repr`), so its ratios are `nan`. FC3 "
        "(degree 1) was measured on the same mesh at the same width, so its ratio sets the "
        "degree-1 estimator beside the degree-2 window: one mesh, one width, two degrees, two "
        "quantities. Measured values are the strings of `docs/testing/xl-ledger.md`; the "
        "estimator values are the CLI's printed `repr` of `femem run --estimate <file>` for "
        f"`{fc4['row']['label']}` and `{fc3['row']['label']}`.",
        "",
        "| Row | Log | ranks | cells | tag 403 | peak GiB | seconds | ratio GiB | ratio seconds |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
        f"| `xl` entry {p['entry']} (`{p['chunk']}`) | `{p['log']}` | {p['ranks']} | "
        f"{p['cells']} | {_cell(TAG403_3T_D2[p['entry']])} | {p['gib']} | {p['seconds']} | "
        "-- | -- |",
    ]
    for m, name in ((fc4, "FC4 (degree 2)"), (fc3, "FC3 (degree 1)")):
        row, r = m["row"], m["ratios"][0]
        pr = row["printed"]
        lines.append(
            f"| `Study.estimate` {name} | `{_cell(row['label'])}` | `{_cell(pr['ranks'][0])}` | "
            f"`{_cell(pr['cells'][0])}` | -- | `{_cell(pr['predicted_peak_gib'][0])}` | "
            f"`{_cell(pr['predicted_seconds'][0])}` | {r['gib']} | {r['seconds']} |"
        )
    lines += [
        "",
        f"FC4's `not_priced`: `{_cell(fc4['row']['printed']['not_priced'][0])}`.",
        "",
        "Ratio = estimator / measured. A record, not a calibration; no row is fitted.",
        "",
    ]
    return "\n".join(line.rstrip() for line in lines)


def render(rows, measured=None, measured_15t=None, measured_d2=None,
           measured_3t_d2=None) -> str:
    lines = [
        HEADER,
        "",
        "# Resource estimates for the shipped study files",
        "",
        "`femem run --estimate <file>` on every study file under `examples/` and every "
        "study fence (```` ```yaml <tag> ````) under `docs/guide/`, printed with no solve. "
        f"Launch width: `mpiexec -n {LAUNCH_RANKS}` -- `--estimate` prices the current "
        "launch width, so `ranks` and `predicted_seconds` are this width's. Each value "
        "is the CLI's printed `repr`, asserted equal to `Study.estimate` "
        "(`tests/package/test_wf17_estimate_table_page.py`).",
        "",
        "The tier is the tier **by prediction**: `xl` when ranks > 12 or predicted "
        "seconds > 1200; otherwise the smallest of smoke (30 s), standard (180 s) and "
        "heavy (1200 s) that covers the predicted seconds. It is not a measurement.",
        "",
        "| Entry | Family | " + " | ".join(COLUMNS) + " | Tier | Refusal |",
        "| --- | --- | " + " | ".join("---" for _ in COLUMNS) + " | --- | --- |",
    ]
    for r in rows:
        label = r["label"]
        marks = ""
        if r["load_result"]:
            marks += " [1]"
        if r["thermal_unpriced"] and r["est"] is not None and not r["est"].not_priced:
            marks += " [2]"
        if r["rc"] == 0:
            vals = [f"`{_cell(r['printed'][c][0])}`" for c in COLUMNS]
            tier, refusal = tier_for(r["est"]), "--"
        else:
            vals = ["--"] * len(COLUMNS)
            tier, refusal = "--", f"rc {r['rc']}: {_cell(r['api_error'] or '')}"
        lines.append(
            f"| `{_cell(label)}`{marks} | `{r['family']}` | " + " | ".join(vals)
            + f" | {tier} | {refusal} |"
        )
    lines += [
        "",
        "[1] A `load_result:` file: the row prices a mesh build and solve that "
        "`femem run` never performs (it tunes the saved result). Recorded, not fixed.",
        "",
        "[2] A file with a `thermal:` section whose `not_priced` is empty: the row "
        "prices the EM half only. Recorded, not fixed.",
        "",
    ]
    text = "\n".join(line.rstrip() for line in lines)
    if measured_15t is not None:  # before the 3 T section, whose test reads heading to EOF
        text += "\n" + render_measured(measured_15t, MEASURED_15T_HEADING, "1.5 T", FC15_MESH,
                                       tags=TAG403_15T)
    if measured is not None:
        text += "\n" + render_measured(measured)
    if measured_3t_d2 is not None:  # between the 3 T and d2 sections (U166; d2 stays last)
        text += "\n" + render_measured_3t_d2(measured_3t_d2)
    if measured_d2 is not None:  # after the 3 T section (U163 layout ruling)
        text += "\n" + render_measured_d2(measured_d2)
    return text


def build(comm=MPI.COMM_WORLD, root: Path = REPO):
    rows, count = collect(comm, root)
    measured = measured_3t(comm)  # no mesh on this geometry; outside the build memo
    measured_15t = measured_3t(comm, fc_file=FC15_FILE, entries=MEASURED_15T)
    m_d2 = measured_d2(comm)
    m_3t_d2 = measured_3t_d2(comm)
    text = (render(rows, measured, measured_15t, measured_d2=m_d2, measured_3t_d2=m_3t_d2)
            if comm.rank == 0 else None)
    return comm.bcast(text, root=0), rows, count


def diff_against(out: Path, text: str) -> list:
    """Unified diff of the tracked page ``out`` against a fresh render (empty = identical)."""
    import difflib

    old = out.read_text(encoding="utf-8") if out.is_file() else ""
    return list(difflib.unified_diff(old.splitlines(), text.splitlines(),
                                     fromfile=str(out), tofile=f"{out.name} (fresh)", lineterm=""))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(DEFAULT_OUT))
    parser.add_argument("--check", action="store_true", help="diff against --out, exit 1 on change")
    args = parser.parse_args(argv)
    comm = MPI.COMM_WORLD
    out = Path(args.out)
    if not out.is_absolute():
        out = Path(os.getcwd()) / out
    text, rows, count = build(comm)
    rc = 0
    if comm.rank == 0:
        print(f"Estimate page: builds={count['builds']} hits={count['hits']} "
              f"geometries={count['keys']}", flush=True)
        for r in rows:
            print(f"  {r['label']}: rc={r['rc']} est={r['est']!r} err={r['api_error']!r}", flush=True)
        if args.check:
            same = out.is_file() and out.read_bytes() == text.encode("utf-8")
            print(f"Estimate page check: {'identical' if same else 'DIFFERS'}: {out.name}", flush=True)
            for ln in diff_against(out, text):
                print(ln, flush=True)
            rc = 0 if same else 1
        else:
            out.parent.mkdir(parents=True, exist_ok=True)
            with open(out, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(text)
            print(f"Estimate page: wrote {out.name}", flush=True)
    return comm.bcast(rc, root=0)


if __name__ == "__main__":
    sys.exit(main())
