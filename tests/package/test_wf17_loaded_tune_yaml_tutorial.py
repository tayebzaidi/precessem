"""`WF-17` step 23 (§9 item 681, §10 U101): the YAML ``load_result:`` page.

``docs/guide/tuning_a_saved_result_from_yaml.md`` carries three blocks:
(a) ``loaded-yaml-solve``, the four-port solve that writes the Touchstone file
(U99's route, ``tests/package/test_wf1_yaml_load_result.py:245-249``); (b)
``loaded-yaml-file``, the study file with ``load_result:`` and ``tune: {k: 1}``;
(c) ``loaded-yaml-cli``, ``femem check`` and ``femem run`` on (b).  This module
execs (a) verbatim, writes (b) verbatim and runs (c)'s two ``femem`` lines
through ``cli.main`` (as U99's (CLI = API) does), all in one scratch directory
under the gitignored ``logs/``.

Brief: docs/planning/chunks/WF-17.md, "§9 item 681".

Anchors (every band and record imported, none set here): (X) the CLI
summary's ``C_mode`` vs ``ITEM429_C_MODE_F`` within
``LEG_D0_REPRODUCTION_BAND`` (U99's (CM)), asserted; the ``.9e`` digits
printed.  (S0) block (c) runs under U99's ``_no_mesh`` patches (imported):
0 solve / build / gmsh calls.  (CLI = API) U99's ``_summary_diff`` of the
CLI's ``summary.json`` against ``run_study_file``'s on the same file is
empty (``elapsed_s`` excepted), and U99's ``_tune_block_hex_equal`` against
``load_study`` + ``load_result`` + ``tune(1)`` is empty.  (Y) block (b)
parses to U99's ``_loaded_mapping`` of the tracked tune file.  (Ck) the
page's ``femem check`` line, with the scratch directory substituted, is a
line of ``femem check``'s output.  (D) docrefs 0 dead.  Negative control
(asserted; U90b measured 0.48993 for k = 2 on this fixture): block (b) with
``k: 2`` misses ``ITEM429_C_MODE_F`` by more than the band.

The solve node runs only with ``FEM_EM_WF17_S23=1`` (then it asserts the
complex build and ``-n 2``, never skips).

`WF-17` step 28 (§9 item 744, §10 U114): five stale sentences, four on this page and
one on ``tuning_a_saved_result.md``.  (M) a two-frequency file (written under the
gitignored ``logs/``, so U71's file is not needed) runs with ``tune.at_hz`` and is
refused without it, by a ``ValueError`` naming ``at_hz``; (H) ``load_result:`` on an
``f_human_birdcage`` parses with ``conductor`` beside it, on a stub path, nothing
solved, and ``conductor`` stays refused on a ``four_port_birdcage``.  Negative
control: each old sentence is absent, and a copy of the page with it restored fails.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import os
import shlex
import shutil
import time
from pathlib import Path

import pytest
import yaml
from mpi4py import MPI

from fem_em_solver import Study, StudyConfigError, cli, load_study, run_study_file
from fem_em_solver import study as study_mod
from tests.complex_mode import IS_COMPLEX, complex_only
from tests.package.test_wf1_study_yaml import _rc_all
from tests.package.test_wf1_yaml_load_result import (
    _loaded_mapping,
    _no_mesh,
    _summary_diff,
    _tune_block_hex_equal,
    _write_yaml,
)
from tests.package.test_wf1_yaml_load_result_f_human import _loaded_mapping as _human_mapping
from tests.package.test_wf1_yaml_load_result_sweep import _doc
from tests.package.test_wf16_load_result_sweep import FREQS
from tests.package.test_wf16_study_load_result import _synthetic_s, _write
from tests.package.test_wf16_study_tune_capacitors import ITEM429_C_MODE_F
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.package.test_wf17_tune_tutorial import _assert_gate_inputs
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "tuning_a_saved_result_from_yaml.md"
PAGE_ID = "docs/guide/tuning_a_saved_result_from_yaml.md"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s23"
S23_ENV = "FEM_EM_WF17_S23"
MIN_GUIDE_PAGES = 66

SOLVE_TAG = "loaded-yaml-solve"
FILE_TAG = "loaded-yaml-file"
CLI_TAG = "loaded-yaml-cli"
CHECK_TAG = "loaded-yaml-check"
BANDS_TAG = "loaded-yaml-bands"
TOUCHSTONE = "f_small_pec_64mhz.s4p"
YAML_NAME = "f_small_pec_load.yaml"
WORK_DIR_PLACEHOLDER = "/your/work/dir"
K1_LINE = "  k: 1\n"
K2_LINE = "  k: 2\n"

# `WF-17` step 28 (§9 item 744): the five old sentences, whitespace normalised.
OTHER_PAGE = REPO / "docs" / "guide" / "tuning_a_saved_result.md"
OLD_HERE = (
    "`load_result` checks that the file holds one frequency equal to `study.frequency_hz`",
    "`load_result:` is defined on a `four_port_birdcage` with `excitation.kind: "
    "lumped_port_drive` only",
    "`load_result:` is four-port only",
    "A Touchstone file holding several frequencies.",
)
OLD_OTHER = ("that `load_result` refuses the human-scale Study",)


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _text() -> str:
    return PAGE.read_text()


def _fence_openers(text: str) -> list:
    return [line for line in text.splitlines() if line.startswith("```") and line != "```"]


def _bands_line() -> str:
    return (
        f"Saved-file C_mode from femem run vs the solved-Study record: "
        f"{LEG_D0_REPRODUCTION_BAND:.0e}"
    )


def _femem_commands(block: str) -> list:
    """Block (c)'s ``femem`` argument lists, in page order."""
    lines = [ln for ln in block.splitlines() if ln.strip()]
    assert lines[0] == "source /usr/local/bin/dolfinx-complex-mode", lines
    out = []
    for ln in lines[1:]:
        words = shlex.split(ln)
        assert words[:4] == ["mpiexec", "-n", "2", "femem"], words
        out.append(words[4:])
    return out


def _knob() -> bool:
    raw = os.environ.get(S23_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S23_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _scratch(subdir: str) -> Path:
    comm = MPI.COMM_WORLD
    d = SCRATCH_ROOT / subdir
    if comm.rank == 0:
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
    comm.Barrier()
    return d


def _write_text(path: Path, text: str) -> Path:
    comm = MPI.COMM_WORLD
    if comm.rank == 0:
        path.write_text(text)
    comm.Barrier()
    return path


def _cli(args, capsys):
    """``cli.main(args)`` with its stdout captured (rank 0 prints)."""
    pre = capsys.readouterr().out  # keep this module's own earlier prints in the log
    with capsys.disabled():
        if MPI.COMM_WORLD.rank == 0 and pre:
            print(pre, end="", flush=True)
    rc = cli.main(list(args))
    out = capsys.readouterr().out
    return _rc_all(rc, MPI.COMM_WORLD), out


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s23_blocks_and_check(capsys, monkeypatch):
    """Page shape, (Y), (Ck) and (E); no mesh."""
    comm = MPI.COMM_WORLD
    text = _text()
    openers = _fence_openers(text)
    _p0(f"\n[WF-17 s23] fences {openers}")
    assert openers == [
        f"```python {SOLVE_TAG}",
        f"```yaml {FILE_TAG}",
        f"```bash {CLI_TAG}",
        f"```text {CHECK_TAG}",
        f"```text {BANDS_TAG}",
    ], openers
    solve = extract_tagged_block(text, SOLVE_TAG)
    yml = extract_tagged_block(text, FILE_TAG, lang="yaml")
    cmds = _femem_commands(extract_tagged_block(text, CLI_TAG, lang="bash"))
    check = extract_tagged_block(text, CHECK_TAG, lang="text").strip("\n").splitlines()
    bands = extract_tagged_block(text, BANDS_TAG, lang="text").strip("\n").splitlines()
    for label, b in (("solve", solve), ("file", yml)):
        _p0(f"[WF-17 s23] block {label} sha256 {_sha(b.encode())}")

    # (a) is U99's route: one Study, one solve, one write_touchstone to TOUCHSTONE
    assert solve.count("Study(") == 1 and solve.count(".solve(") == 1
    assert "load_result" not in solve and "tune" not in solve
    assert solve.count(f'write_touchstone("{TOUCHSTONE}")') == 1
    # (Y) (b) is the tracked tune file plus load_result: <TOUCHSTONE>
    parsed = yaml.safe_load(yml)
    want = _loaded_mapping(TOUCHSTONE)
    _p0(
        f"[WF-17 s23] (Y) page file keys {list(parsed)}; equal to U99's _loaded_mapping "
        f"{parsed == want}"
    )
    assert parsed == want
    assert yml.count(K1_LINE) == 1
    # (c): exactly `femem check` then `femem run`, on the page's file name
    _p0(f"[WF-17 s23] (c) femem commands {cmds}")
    assert cmds == [["check", YAML_NAME], ["run", YAML_NAME, "--output-dir", "femem_out"]]
    assert f"Save this as `{YAML_NAME}`" in text
    # (E) the bands block
    _p0(f"[WF-17 s23] (E) bands block {bands} (rendered {[_bands_line()]})")
    assert bands == [_bands_line()]
    for token in ("AED", "HFSS"):
        assert token not in text, token
    for link in ("(tuning_a_saved_result.md)", "(tuning_from_yaml.md)"):
        assert link in text, link

    # (Ck) femem check on the page's file, in a scratch directory (check reads no S file)
    d = _scratch(f"check_n{comm.size}")
    _write_text(d / YAML_NAME, yml)
    calls = {"n": 0}
    monkeypatch.chdir(d)
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        rcs, out = _cli(cmds[0], capsys)
    page_line = check[0].replace(WORK_DIR_PLACEHOLDER, str(d.resolve()), 1)
    _p0(f"[WF-17 s23] (Ck) femem check rc {rcs}:\n{out.rstrip()}\n  page line -> {page_line!r}")
    assert len(check) == 1 and WORK_DIR_PLACEHOLDER in check[0]
    assert rcs == (0, 0)
    if comm.rank == 0:
        assert page_line in out.splitlines()
    assert comm.allreduce(calls["n"], op=MPI.SUM) == 0


def test_wf17_s23_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s23] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
        f"exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


def _stale(text: str, olds) -> list:
    flat = " ".join(text.split())
    return [old for old in olds if old in flat]


@complex_only
def test_wf17_s28_stale_lines_corrected(tmp_path, capsys):
    """The five corrected sentences: their behaviour, then their absence.  No mesh."""
    comm = MPI.COMM_WORLD
    here, other = _text(), OTHER_PAGE.read_text()
    assert _stale(here, OLD_HERE) == [] and _stale(other, OLD_OTHER) == []
    copy = tmp_path / "restored.md"
    missed = []
    for old in OLD_HERE + OLD_OTHER:  # control: a restored sentence fails the check
        copy.write_text(here + "\n" + old + "\n", encoding="utf-8")
        if _stale(copy.read_text(encoding="utf-8"), (old,)) != [old]:
            missed.append(old)
    assert missed == []

    d = _scratch(f"s28_n{comm.size}")
    s = _synthetic_s()
    multi = _write(d / "two.s4p", [(FREQS[1], s), (FREQS[2], 0.9 * s)])
    ok_yml = _write_yaml(d / "at_hz.yaml", _doc(multi, {"k": 1, "at_hz": FREQS[1]}))
    no_yml = _write_yaml(d / "no_at_hz.yaml", _doc(multi, {"k": 1}))
    four_conductor = _doc(multi, {"k": 1, "at_hz": FREQS[1]})
    four_conductor["conductor"] = {}
    four_yml = _write_yaml(d / "four_conductor.yaml", four_conductor)
    human_yml = _write_yaml(d / "human.yaml", _human_mapping("stub_absent.s32p"))
    calls = {"n": 0}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        run_study_file(ok_yml, output_dir=d / "ok", comm=comm)  # (M)
        with pytest.raises(ValueError, match="at_hz") as err:
            run_study_file(no_yml, output_dir=d / "no", comm=comm)
        capsys.readouterr()
        rc_cli = cli.main(["run", str(no_yml), "--output-dir", str(d / "no_cli")])
        err_text = capsys.readouterr().err
        human = load_study(human_yml)  # (H), the stub path is not read at parse
        with pytest.raises(StudyConfigError) as err_four:
            load_study(four_yml)
    n_calls = comm.allreduce(calls["n"], op=MPI.SUM)
    rcs = tuple(comm.allgather(rc_cli))
    named = comm.bcast("tune.at_hz" in err_text, root=0)
    with open(d / "ok" / "summary.json") as fh:
        selected = json.load(fh)["tune"]["frequency_hz"]
    with capsys.disabled():  # printed past capsys
        _p0(
            f"\n[WF-17 s28] old sentences absent {len(OLD_HERE) + len(OLD_OTHER)}, restored "
            f"copies caught {len(OLD_HERE) + len(OLD_OTHER) - len(missed)}"
        )
        _p0(
            f"[WF-17 s28 M] with at_hz: selected {selected:.3e} Hz; without: run_study_file "
            f"{type(err.value).__name__}: {err.value}; cli.main rcs {rcs}; "
            f"StudyConfigError {isinstance(err.value, StudyConfigError)} (printed)"
        )
        _p0(
            f"[WF-17 s28 H] {type(human.geometry).__name__}, conductor "
            f"{type(human.conductor).__name__}, stub on disk {(d / 'stub_absent.s32p').exists()}; "
            f"four_port conductor refused by field {err_four.value.field!r}; "
            f"{n_calls} solve / build / gmsh calls"
        )
    assert selected == FREQS[1]
    assert rcs == (2, 2) and named
    assert not (d / "stub_absent.s32p").exists()
    assert type(human.geometry).__name__ == "FHumanBirdcageGeometry"
    assert human.conductor is not None
    assert err_four.value.field == "conductor"
    assert n_calls == 0


# ------------------------------------------------------------- solve window (-n 2)


def test_wf17_s23_yaml_loaded_tune(capsys, monkeypatch):
    if not _knob():
        pytest.skip(
            f"{S23_ENV} unset: the YAML loaded-tune solve window is environment-selected "
            "(complex build, -n 2)"
        )
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the four-port Study needs the complex build"
    assert comm.size == 2, comm.size
    text = _text()
    solve = extract_tagged_block(text, SOLVE_TAG)
    yml = extract_tagged_block(text, FILE_TAG, lang="yaml")
    cmds = _femem_commands(extract_tagged_block(text, CLI_TAG, lang="bash"))
    assert yml.count(K1_LINE) == 1

    d = _scratch(f"run_n{comm.size}")
    monkeypatch.chdir(d)

    # (a) exec'd verbatim; one solve, on `solved`
    real_solve = Study.solve
    solves = []

    def solve_spy(self, *a, **k):
        solves.append(self)
        return real_solve(self, *a, **k)

    t0 = time.perf_counter()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(Study, "solve", solve_spy)
        ns = {"__name__": "__main__"}
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            exec(compile(solve, f"{PAGE_ID}#{SOLVE_TAG}", "exec"), ns)
    t_solve = time.perf_counter() - t0
    _p0(
        f"\n[WF-17 s23] block (a) stdout (rank 0): {buf.getvalue()!r}; {t_solve:.1f} s; "
        f"solves {len(solves)}"
    )
    assert len(solves) == 1 and solves[0] is ns["solved"]
    _assert_gate_inputs(ns["solved"], "solved")
    s_path = Path(ns["path"]).resolve()
    assert s_path == (d / TOUCHSTONE).resolve() and s_path.exists()

    # (b) written verbatim; its k = 2 twin for the control
    yml_path = _write_text(d / YAML_NAME, yml)
    yml2_path = _write_text(d / "f_small_pec_load_k2.yaml", yml.replace(K1_LINE, K2_LINE, 1))

    # (c) through cli.main, and the API routes, under U99's _no_mesh patches
    calls = {"n": 0}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, calls)
        rc_check, out_check = _cli(cmds[0], capsys)
        t0 = time.perf_counter()
        rc_run, out_run = _cli(cmds[1], capsys)
        t_cli = time.perf_counter() - t0
        run = run_study_file(yml_path, output_dir=d / "api_out", comm=comm)
        st = load_study(yml_path)
        st.load_result(s_path)
        api = st.tune(k=1, comm=comm)
        rc_k2, _ = _cli(["run", yml2_path.name, "--output-dir", "femem_out_k2"], capsys)
    n_calls = comm.allreduce(calls["n"], op=MPI.SUM)

    with open(d / "femem_out" / "summary.json") as fh:
        js_cli = json.load(fh)
    with open(d / "api_out" / "summary.json") as fh:
        js_api = json.load(fh)
    with open(d / "femem_out_k2" / "summary.json") as fh:
        js_k2 = json.load(fh)

    # (S0)
    _p0(
        f"[WF-17 s23] (S0) block (c) under _no_mesh: {n_calls} solve / build / gmsh calls; "
        f"femem check rc {rc_check}, femem run rc {rc_run} ({t_cli:.3f} s), k2 rc {rc_k2}; "
        f"result {type(run.study.result).__name__}; summary keys {sorted(js_cli)}"
    )
    _p0(f"[WF-17 s23] femem check (rank 0):\n{out_check.rstrip()}")

    # (CLI = API)
    diff = _summary_diff(js_api, js_cli)
    bad = _tune_block_hex_equal(js_cli["tune"], api)
    _p0(
        f"[WF-17 s23] (CLI=API) summary fields differing (elapsed_s excepted) {diff}; tune block "
        f"vs load_study + load_result + tune(1) differing {bad}"
    )

    # (X)
    c = js_cli["tune"]["c_mode_farad"]
    rel = abs(c / ITEM429_C_MODE_F - 1.0)
    rel_g = comm.allreduce(rel, op=MPI.MAX)
    _p0(
        f"[WF-17 s23] (X) femem run C_mode {c:.12e} F (hex {float(c).hex()}); '{c:.9e}' vs "
        f"ITEM429_C_MODE_F '{ITEM429_C_MODE_F:.9e}' equal {f'{c:.9e}' == f'{ITEM429_C_MODE_F:.9e}'} "
        f"(printed); rel {rel_g:.3e} (ASSERTED <= LEG_D0_REPRODUCTION_BAND "
        f"{LEG_D0_REPRODUCTION_BAND:.0e})"
    )

    # negative control: the page's file with k: 2
    c2 = js_k2["tune"]["c_mode_farad"]
    miss = abs(c2 / ITEM429_C_MODE_F - 1.0)
    miss_g = comm.allreduce(miss, op=MPI.MIN)
    _p0(
        f"[WF-17 s23] (c-k2) k: 2 file, k {js_k2['tune']['k']}, C {c2:.12e} F: miss {miss_g:.5f} "
        f"(ASSERTED > {LEG_D0_REPRODUCTION_BAND:.0e}; U90b measured 0.48993)"
    )

    assert n_calls == 0
    assert rc_check == (0, 0) and rc_run == (0, 0) and rc_k2 == (0, 0)
    assert isinstance(run.study.result, study_mod.LoadedPortResult)
    assert js_cli["load_result"] == str(s_path)
    assert "cells" not in js_cli and "touchstone" not in js_cli
    assert js_cli["ranks"] == comm.size == 2
    assert diff == []
    assert bad == []
    assert js_cli["tune"]["k"] == 1
    assert rel_g <= LEG_D0_REPRODUCTION_BAND
    assert js_k2["tune"]["k"] == 2
    assert miss_g > LEG_D0_REPRODUCTION_BAND
