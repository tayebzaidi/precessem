"""`WF-17` step 15: the capacitor-sheets-from-YAML tutorial page runs a copper
F-small study file with a 0.2 ohm series-RC sheet on P2 through
``run_study_file`` and reproduces `WF-16` step 10's per-sheet circuit-reduction
identity, with the unterminated 4x4's copper tie to the `ANS-6` degree-1 record.

Brief: docs/planning/chunks/WF-17.md, "§9 items 627, 632, 633, 634" (item 627),
over the template "§9 item 621".  The solve node runs only with
``FEM_EM_WF17_S15=1`` (and then asserts the complex build and ``-n
ANS6_EXACT_STRIP_RECORD_RANK_WIDTH``, never skips).  The other nodes build no
mesh.

Anchors: (X) one block per tag (three tags); the YAML block equals the tracked
file byte for byte; the python block holds one ``run_study_file(`` call whose
argument is the tracked file's repository path and no ``Study(`` call (AST
names); ``femem check`` on the file exits 0 with one ``family:`` line; (In) the
file's parsed Study equals the gate's inputs exactly (no mesh); (E) the bands
block = the three imported constants at ``.0e``; (R0) S0's cells and class
means vs the copper record at ``LEG_D0_REPRODUCTION_BAND``; (I) the run's 3x3
vs the circuit reduction of S0 with Z = r + 1/(jwC) on P2 at
``REDUCTION_BAND``; (D) 0 dead links.  Controls (rule (e)): (N1) a scratch
copy of the file with ``esr_ohm: 0.4`` misses the 0.2 ohm reduction (asserted,
backed by step 14's 10.953x); (c-page) by construction.  Bands only are
asserted; every reading is printed.

Mesh economy (disclosed, as the gate's ``s10`` fixture): inside the solve node
``FourPortBirdcageGeometry.build`` is wrapped so the first call builds the real
mesh and later equal-geometry calls on the same communicator return it.  The
block names the file by its repository path; it is exec'd from the scratch, so
rank 0 copies the tracked file to the same relative path under the scratch
first (disclosed adaptation; the copy's sha256 is printed beside the tracked
one).
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import io
import os
import shutil
import time
from pathlib import Path

import pytest
from mpi4py import MPI

import fem_em_solver
from fem_em_solver import (
    CapacitorSheets,
    FourPortBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    cli,
)
from fem_em_solver.ports.circuit import reduce_terminated_ports
from fem_em_solver.study import COPPER_SIGMA, GATED_WIDTH_FRACTION
from fem_em_solver.study_yaml import load_study
from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf16_study_capacitor_sheets import PORTS, _phantom, _study
from tests.package.test_wf16_study_copper import _geometry
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    extract_tagged_block,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_ans6_exact_strip import (
    ANS6_EXACT_STRIP_CLASS_MEANS_P1,
    ANS6_EXACT_STRIP_HOLE_CELL_RECORD,
    ANS6_EXACT_STRIP_RECORD_RANK_WIDTH,
    CLASSES,
    _class_means,
)
from tests.validation.test_ans7_exact_strip import C_F
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port24_series_rc_sheet import (
    ESR_OHM,
    PLANT_INDEX,
    PLANT_OHM,
    _rc,
    _rel,
)
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "capacitor_sheets_from_yaml.md"
PAGE_ID = "docs/guide/capacitor_sheets_from_yaml.md"
FILE_REL = "tests/package/wf1/f_small_copper_capacitors_esr.yaml"
YAML_FILE = REPO / FILE_REL
BASE_FILE = REPO / "tests" / "package" / "wf1" / "f_small_copper.yaml"
LOSSLESS_REL = "tests/package/wf1/f_small_hole_pec_capacitors.yaml"
SCRATCH_ROOT = REPO / "logs" / "wf17" / "s15"
S15_ENV = "FEM_EM_WF17_S15"
MIN_GUIDE_PAGES = 64

YAML_TAG = "capacitor-esr-yaml"
RUN_TAG = "capacitor-esr-run"
BANDS_TAG = "capacitor-esr-yaml-bands"

# Printed beside the readings (`20261025T004502Z_WF-17-s14-Wesr.log` :941 /
# :942 / :943); never asserted.
S14_R0 = "4.093e-11 / 8.673e-11 / 1.376e-10"
S14_I_P2 = "3.036163e-05"
S14_N1 = "1.095344e-02 = 10.953x"


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _text() -> str:
    return PAGE.read_text()


def _calls_named(block: str, name: str) -> list:
    out = []
    for node in ast.walk(ast.parse(block)):
        if isinstance(node, ast.Call):
            fn = node.func
            fname = fn.id if isinstance(fn, ast.Name) else getattr(fn, "attr", None)
            if fname == name:
                out.append(node)
    return out


def _run_defect(block: str):
    """None when the run block has the (X) shape, else the reason."""
    runs = _calls_named(block, "run_study_file")
    if len(runs) != 1:
        return f"{len(runs)} run_study_file calls; need exactly one"
    arg = ast.literal_eval(runs[0].args[0])
    if arg != FILE_REL:
        return f"run_study_file argument {arg!r} != {FILE_REL!r}"
    for name in ("Study", "solve", "solve_sweep", "tune"):
        if _calls_named(block, name):
            return f"the block calls {name}"
    return None


def _yaml_defect(block: str):
    data = YAML_FILE.read_text()
    if block != data:
        return f"YAML block ({len(block)} chars) != {FILE_REL} ({len(data)} chars)"
    return None


def _bands_lines() -> list:
    return [
        f"Study S vs the circuit reduction with Z = r + 1/(jωC): {REDUCTION_BAND:.0e}",
        f"Unterminated 4×4 vs the copper record: {LEG_D0_REPRODUCTION_BAND:.0e}",
        f"Series resistance r (Ω): {ESR_OHM:.0e}",
    ]


def _bands_defect(block: str):
    got = block.strip("\n").splitlines()
    want = _bands_lines()
    if len(got) != len(want):
        return f"{len(got)} lines != {len(want)}"
    for g, w in zip(got, want):
        if g != w:
            return f"line {g!r} != {w!r}"
    return None


def _knob() -> bool:
    raw = os.environ.get(S15_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S15_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _solve_preamble():
    if not _knob():
        pytest.skip(
            f"{S15_ENV} unset: the tutorial solve window is environment-selected "
            "(complex build, -n 2)"
        )
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the copper four-port capacitor-sheet Study needs the complex build"
    assert comm.size == ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, comm.size
    return comm


def _plant_text() -> str:
    data = YAML_FILE.read_text()
    assert data.count("  esr_ohm: 0.2\n") == 1
    return data.replace("  esr_ohm: 0.2\n", f"  esr_ohm: {PLANT_OHM!r}\n")


def exec_run_block(monkeypatch):
    """Clear the scratch, copy the file to its relative path there, chdir, exec."""
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), RUN_TAG)
    data = YAML_FILE.read_bytes()
    if comm.rank == 0:
        print(
            f"[WF-17 s15] page sha256 {_sha(raw)}; block {RUN_TAG} sha256 "
            f"{_sha(block.encode())}",
            flush=True,
        )
        if SCRATCH_ROOT.exists():
            shutil.rmtree(SCRATCH_ROOT)
        target = SCRATCH_ROOT / FILE_REL
        target.parent.mkdir(parents=True)
        target.write_bytes(data)
        print(
            f"[WF-17 s15] copied {FILE_REL} sha256 {_sha(target.read_bytes())} (tracked "
            f"{_sha(data)})",
            flush=True,
        )
    comm.Barrier()
    monkeypatch.chdir(SCRATCH_ROOT)
    calls = {"n": 0, "paths": []}
    original = fem_em_solver.run_study_file

    def spy(path, *args, **kwargs):
        calls["n"] += 1
        calls["paths"].append(str(path))
        return original(path, *args, **kwargs)

    monkeypatch.setattr(fem_em_solver, "run_study_file", spy)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{RUN_TAG}", "exec")
    with contextlib.redirect_stdout(buf):
        exec(code, namespace)
    out = buf.getvalue()
    _p0(f"[WF-17 s15] block {RUN_TAG} stdout (rank 0): {out!r}")
    monkeypatch.setattr(fem_em_solver, "run_study_file", original)
    return namespace, calls, out


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s15_blocks_and_inputs(tmp_path):
    text = _text()
    yblock = extract_tagged_block(text, YAML_TAG, lang="yaml")
    rblock = extract_tagged_block(text, RUN_TAG)
    extract_tagged_block(text, BANDS_TAG, lang="text")
    y_reason, r_reason = _yaml_defect(yblock), _run_defect(rblock)
    _p0(
        f"[WF-17 s15] (X) {YAML_TAG}: sha256 {_sha(yblock.encode())} vs tracked "
        f"{_sha(YAML_FILE.read_bytes())}; defect {y_reason}"
    )
    _p0(
        f"[WF-17 s15] (X) {RUN_TAG}: run_study_file {len(_calls_named(rblock, 'run_study_file'))}, "
        f"Study {len(_calls_named(rblock, 'Study'))}, solve {len(_calls_named(rblock, 'solve'))}; "
        f"defect {r_reason}"
    )
    assert y_reason is None, y_reason
    assert r_reason is None, r_reason
    assert "Study(" not in rblock and rblock.count("run_study_file(") == 1
    assert f"(../../tests/package/wf1/{YAML_FILE.name})" in text
    assert "(capacitor_sheets_with_esr.md)" in text
    assert LOSSLESS_REL in text and "**lossless**" in text
    assert text.count("```yaml ") == 1

    # the tracked file is f_small_copper.yaml except the header and excitation
    base = BASE_FILE.read_text().splitlines()
    new = YAML_FILE.read_text().splitlines()
    i_b, i_n = base.index("excitation:"), new.index("excitation:")
    _p0(
        f"[WF-17 s15] (X) file vs {BASE_FILE.name}: body before excitation equal "
        f"{base[1:i_b] == new[1:i_n]}; after equal {base[i_b + 2:] == new[i_n + 5:]}; "
        f"excitation {new[i_n:i_n + 5]}"
    )
    assert base[1:i_b] == new[1:i_n]
    assert base[i_b + 1] == "  kind: lumped_port_drive"
    assert base[i_b + 2 :] == new[i_n + 5 :]
    assert len(new[0]) > 0 and new[0].startswith("# ") and not base[0] == new[0]

    # femem check on the file (copper + capacitor_sheets must be accepted)
    comm = MPI.COMM_WORLD
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = cli.main(["check", str(YAML_FILE)])
    lines = buf.getvalue().splitlines()
    rc_all = (comm.allreduce(int(rc), op=MPI.MIN), comm.allreduce(int(rc), op=MPI.MAX))
    _p0(f"[WF-17 s15] (X) femem check rc {rc_all}; lines {lines!r}")
    assert rc_all == (0, 0)
    if comm.rank == 0:
        fam = [ln for ln in lines if ln.startswith("family: ")]
        assert fam == ["family: four_port_birdcage"], lines

    # (c-page), by construction: a removed / doubled tag raises
    for tag, lang in ((YAML_TAG, "yaml"), (RUN_TAG, "python"), (BANDS_TAG, "text")):
        with pytest.raises(LookupError):
            extract_tagged_block(
                text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag, lang=lang
            )
        doubled = (
            text
            + "\n"
            + f"```{lang} {tag}\n"
            + extract_tagged_block(text, tag, lang=lang)
            + "```\n"
        )
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    # a tmp copy of the page with esr_ohm: 0.2 removed from the YAML block fails (X)
    stripped = tmp_path / "page.md"
    stripped.write_text(text.replace("  esr_ohm: 0.2\n", "", 1))
    stripped_reason = _yaml_defect(
        extract_tagged_block(stripped.read_text(), YAML_TAG, lang="yaml")
    )
    other_file = tmp_path / "page2.md"
    other_file.write_text(
        text.replace(f'run_study_file("{FILE_REL}")', f'run_study_file("{LOSSLESS_REL}")', 1)
    )
    other_reason = _run_defect(extract_tagged_block(other_file.read_text(), RUN_TAG))
    _p0(
        f"[WF-17 s15] (c-page) tags removed / doubled raise LookupError for all three blocks; "
        f"esr_ohm removed -> {stripped_reason}; lossless file named -> {other_reason}"
    )
    assert stripped_reason is not None
    assert other_reason is not None and "f_small_hole_pec_capacitors" in other_reason
    for token in ("Study(", ".json", ".npz", ".s4p", ".xdmf", "AED", "HFSS"):
        assert token not in text, token

    # (In): the file's Study is the gate's (parse only, no mesh)
    st = load_study(YAML_FILE)
    want_exc = CapacitorSheets(C_F, ESR_OHM, on=(PORTS[PLANT_INDEX],))
    _p0(
        f"[WF-17 s15] (In) geometry {st.geometry!r}; material {st.material!r}; excitation "
        f"{st.excitation!r}; frequency {float(st.frequency_hz).hex()} (want "
        f"{FREQUENCY_64_HZ.hex()}); degree {st.degree}; conductor {st.conductor!r}; c_f hex "
        f"{float(st.excitation.c_f).hex()} (C_F {C_F.hex()})"
    )
    assert type(st.geometry) is FourPortBirdcageGeometry and st.geometry == _geometry()
    assert st.geometry.strip_fraction == GATED_WIDTH_FRACTION
    assert st.material == _phantom()
    assert float(st.material.sigma).hex() == SALINE_SIGMA.hex()
    assert float(st.material.epsilon_r).hex() == SALINE_EPSILON_R.hex()
    assert type(st.excitation) is CapacitorSheets and st.excitation == want_exc
    assert float(st.excitation.c_f).hex() == C_F.hex()
    assert float(st.excitation.esr_ohm).hex() == ESR_OHM.hex()
    assert tuple(st.excitation.on) == ("P2",) and PORTS[PLANT_INDEX] == "P2"
    assert st.excitation.drive == LumpedPortDrive()
    assert float(st.frequency_hz).hex() == FREQUENCY_64_HZ.hex()
    assert st.degree == 1
    assert type(st.conductor) is LeontovichWall and st.conductor == LeontovichWall(COPPER_SIGMA)
    assert float(st.conductor.sigma_s_per_m).hex() == COPPER_SIGMA.hex()
    # the planted copy parses to the same Study with r = PLANT_OHM only
    plant_path = tmp_path / "plant.yaml"
    plant_path.write_text(_plant_text())
    sp = load_study(plant_path)
    _p0(f"[WF-17 s15] (In) planted copy excitation {sp.excitation!r}")
    assert sp.excitation == CapacitorSheets(C_F, PLANT_OHM, on=(PORTS[PLANT_INDEX],))
    assert sp.geometry == st.geometry and sp.conductor == st.conductor


def test_wf17_s15_bands_on_the_page():
    block = extract_tagged_block(_text(), BANDS_TAG, lang="text")
    _p0(f"[WF-17 s15] (E) bands block {block.strip().splitlines()} (rendered {_bands_lines()})")
    assert _bands_defect(block) is None, _bands_defect(block)
    altered = block.replace(f"{REDUCTION_BAND:.0e}", f"{REDUCTION_BAND * 10:.0e}", 1)
    reason = _bands_defect(altered)
    _p0(f"[WF-17 s15] (c-page) one exponent altered -> {reason}")
    assert reason is not None and "circuit reduction" in reason


def test_wf17_s15_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s15] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
        f"exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve window (-n 2)


def test_wf17_s15_yaml_block_reproduces_wf16_s10(monkeypatch):
    comm = _solve_preamble()
    real_build = FourPortBirdcageGeometry.build
    cache = {"key": None, "built": None, "hits": 0, "real": 0}

    def _build(self, comm_):
        if cache["built"] is not None and cache["key"] == self and comm_ is comm:
            cache["hits"] += 1
            return cache["built"]
        cache["real"] += 1
        out = real_build(self, comm_)
        cache["key"], cache["built"] = self, out
        return out

    t_node = time.perf_counter()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _build)
        t0 = time.perf_counter()
        s0_res = _study(LumpedPortDrive()).solve(comm)
        t_s0 = time.perf_counter() - t0
        t0 = time.perf_counter()
        ns, calls, _ = exec_run_block(monkeypatch)
        t_page = time.perf_counter() - t0
        plant_path = SCRATCH_ROOT / "f_small_copper_capacitors_esr_plant.yaml"
        if comm.rank == 0:
            plant_path.write_text(_plant_text())
        comm.Barrier()
        t0 = time.perf_counter()
        plant = fem_em_solver.run_study_file(
            plant_path, output_dir=SCRATCH_ROOT / "plant", comm=comm
        )
        t_plant = time.perf_counter() - t0
    _p0(
        f"[WF-17 s15] build cache: {cache['real']} real build, {cache['hits']} hits (3 Studies); "
        f"S0 {t_s0:.1f} s, page block {t_page:.1f} s, planted run {t_plant:.1f} s"
    )

    # (R0)
    s0 = s0_res.s_matrix
    means = _class_means(s0)
    want = ANS6_EXACT_STRIP_CLASS_MEANS_P1["64"]["cu"]
    rel_r0 = {c: abs(means[c] - want[c]) / abs(want[c]) for c in CLASSES}
    _p0(
        f"[WF-17 s15] (R0) cells {s0_res.cells} (record {ANS6_EXACT_STRIP_HOLE_CELL_RECORD}); "
        "Cu class means vs record: "
        + " / ".join(f"{c} {rel_r0[c]:.3e}" for c in CLASSES)
        + f" (ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e}; s14-Wesr:941 {S14_R0})"
    )

    # (I) the page's identity
    run = ns["run"]
    result = run.result
    summ = run.summary
    ids = [p.port_id for p in result.port_definitions]
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    pred = reduce_terminated_ports(s0, z0, {PLANT_INDEX: _rc(ESR_OHM)})
    r_i = _rel(result.s_matrix, pred)
    _p0(
        f"[WF-17 s15] (I) block spy calls {calls['n']} {calls['paths']}; family "
        f"{summ['family']!r}; kept {summ['kept_ports']}; capacitors {summ['capacitors']}; "
        f"cells {result.cells}"
    )
    _p0(
        f"[WF-17 s15] (I) page ports {ids}; Z {result.capacitor_impedance_ohm!r} (want "
        f"{_rc(ESR_OHM)!r}); S {result.s_matrix.shape}; rel {r_i:.6e} "
        f"({r_i / REDUCTION_BAND:.4f} x REDUCTION_BAND {REDUCTION_BAND:.0e}; s14-Wesr:942 "
        f"{S14_I_P2})"
    )

    # (N1) planted esr_ohm: 0.4 in a scratch copy of the file
    r_n1 = _rel(plant.result.s_matrix, pred)
    _p0(
        f"[WF-17 s15] (N1) planted esr_ohm {plant.summary['capacitors']['esr_ohm']!r} on "
        f"{plant.summary['capacitors']['sheets']}: {r_n1:.6e} = "
        f"{r_n1 / REDUCTION_BAND:.3f} x REDUCTION_BAND (ASSERTED > {REDUCTION_BAND:.0e}; "
        f"s14-Wesr:943 {S14_N1})"
    )
    _p0(f"[WF-17 s15] W-yaml node {time.perf_counter() - t_node:.1f} s")

    assert cache["real"] == 1 and cache["hits"] == 2
    assert calls["n"] == 1 and calls["paths"] == [FILE_REL]
    assert summ["family"] == "four_port_birdcage"
    assert summ["capacitors"]["sheets"] == ["P2"]
    assert float(summ["capacitors"]["esr_ohm"]).hex() == ESR_OHM.hex()
    assert summ["conductor"] == {"kind": "leontovich", "sigma_s_per_m": COPPER_SIGMA}
    assert plant.summary["capacitors"]["esr_ohm"] == PLANT_OHM
    assert s0_res.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert result.cells == ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    assert max(rel_r0.values()) <= LEG_D0_REPRODUCTION_BAND
    assert ids == ["P1", "P3", "P4"]
    assert result.capacitor_impedance_ohm == _rc(ESR_OHM)
    assert r_i <= REDUCTION_BAND
    assert r_n1 > REDUCTION_BAND
