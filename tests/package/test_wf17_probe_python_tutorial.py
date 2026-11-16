"""`WF-17` step 18: the probes-from-Python page builds `MAT-12`'s sphere
``Study`` with `WF-15`'s probe specification, calls ``Study.probe`` and
reproduces `MAT-19`'s whole-ball probe identity.

Brief: docs/planning/chunks/WF-17.md, "§9 items 627, 632, 633, 634" -> "Item
634" (template: "§9 item 621").  The solve node runs only with
``FEM_EM_WF17_S18=1``; set, it asserts the complex build and ``-n 2`` and never
skips.  The other nodes build no mesh.

Anchors: (X) one block per tag (two tags); the python block holds one
``Study(``, one ``.solve(``, one ``.probe(`` call and no ``run_study_file`` (AST
names); (In) with ``Study.solve`` / ``Study.probe`` patched, the block's Study
inputs, probe list and cooling equal the gate's (`WF-15` step 4), floats by
``float.hex``; (E) the bands block = the two lines rendered from ``B_BAND`` /
``WIDTH_BAND``; on the block's run, as ``test_wf15_study_probe.py`` asserts
them: (Once), (A0) submesh cells == ``N_SUB``, (A) whole-ball value and volume
within ``B_BAND``, (U) values == ``ball_average`` on the Study's own θ; (D) 0
dead links.  Control (rule (e)): (c1) ``C2_OUT`` (c2 moved +0.1 mm) misses by
> ``WIDTH_BAND``, asserted, backed by 7.5333e-03
(`20261018T123436Z_WF-15-s6-RG1-probe-n2.log:143`); (c-page) by construction.
Every band, record, constant and comparator is imported, none restated;
readings are printed.
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
from types import SimpleNamespace

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver import QuasiStaticSphereField, RobinCooling, SphereInBoxGeometry, Study
from fem_em_solver.thermal.probe import ProbeBall, ball_average

from tests.complex_mode import IS_COMPLEX
from tests.package.test_wf15_study_lossy_sphere import GEOMETRY, _material_057
from tests.package.test_wf15_study_probe import C2_OUT, PROBES, WHOLE
from tests.package.test_wf17_first_study_tutorial import LINK_PASS_RE, _checker_stdout, extract_tagged_block
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_lossy_sphere_sar import E0, FREQUENCY_HZ, RHO_KG_M3, SIGMA_HIGH
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, G, H, K, T_INF, _allsum
from tests.validation.test_mat19_virtual_temperature_probe import N_SUB, WIDTH_BAND

REPO = Path(__file__).resolve().parents[2]
PAGE = REPO / "docs" / "guide" / "probes_from_python.md"
PAGE_ID = "docs/guide/probes_from_python.md"
SCRATCH = REPO / "logs" / "wf17" / "s18"
S18_ENV = "FEM_EM_WF17_S18"
MIN_GUIDE_PAGES = 67

PROBE_TAG = "probe-python"
BANDS_TAG = "probe-python-bands"

# Printed beside the readings, never asserted:
# `20261018T123436Z_WF-15-s6-RG1-probe-n2.log` :136 / :143.
RG1_A = "1.963e-15 / 1.221e-15"
RG1_C1 = "7.5333e-03"


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


def _block_defect(block: str):
    """None when the block has the brief's call shape, else the reason."""
    want = {"Study": 1, "solve": 1, "probe": 1, "run_study_file": 0}
    for name, n in want.items():
        got = len(_calls_named(block, name))
        if got != n:
            return f"{got} {name!r} calls; need {n}"
    return None


def _bands_lines() -> list:
    return [
        f"Whole-ball probe vs the volume mean: {B_BAND:.0e}",
        f"Displaced probe must miss by more than: {WIDTH_BAND:.0e}",
    ]


def _lines_defect(block: str, want: list):
    got = block.strip("\n").splitlines()
    if len(got) != len(want):
        return f"{len(got)} lines != {len(want)}"
    for g, w in zip(got, want):
        if g != w:
            return f"line {g!r} != {w!r}"
    return None


def _knob():
    raw = os.environ.get(S18_ENV, "").strip()
    if raw == "":
        return None
    if raw != "1":
        raise ValueError(f"{S18_ENV}={raw!r}: only unset/empty or '1' are accepted")
    return raw


def exec_probe_block(monkeypatch):
    """Clear the scratch on rank 0, chdir there on every rank, exec the block as __main__."""
    comm = MPI.COMM_WORLD
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), PROBE_TAG)
    if comm.rank == 0:
        print(f"[WF-17 s18] page sha256 {_sha(raw)}; block sha256 {_sha(block.encode())}",
              flush=True)
        if SCRATCH.exists():
            shutil.rmtree(SCRATCH)
        SCRATCH.mkdir(parents=True)
    comm.Barrier()
    monkeypatch.chdir(SCRATCH)
    namespace = {"__name__": "__main__"}
    buf = io.StringIO()
    code = compile(block, f"{PAGE_ID}#{PROBE_TAG}", "exec")
    with contextlib.redirect_stdout(buf):
        exec(code, namespace)
    out = buf.getvalue()
    _p0(f"[WF-17 s18] block stdout (rank 0): {out!r}")
    return namespace, out


def _hex_ball(p: ProbeBall):
    return tuple(float(c).hex() for c in p.centre), float(p.radius).hex()


# ------------------------------------------------------------- W0 (no mesh)


def test_wf17_s18_blocks_and_inputs(monkeypatch, tmp_path):
    text = _text()
    block = extract_tagged_block(text, PROBE_TAG)
    reason = _block_defect(block)
    _p0(f"[WF-17 s18] (X) {PROBE_TAG}: Study {len(_calls_named(block, 'Study'))}, solve "
        f"{len(_calls_named(block, 'solve'))}, probe {len(_calls_named(block, 'probe'))}, "
        f"run_study_file {len(_calls_named(block, 'run_study_file'))}; defect {reason}")
    assert reason is None, reason
    assert block.count("Study(") == 1 and "run_study_file" not in block
    extract_tagged_block(text, BANDS_TAG, lang="text")

    # (c-page), by construction: a removed / doubled tag raises
    for tag, lang in ((PROBE_TAG, "python"), (BANDS_TAG, "text")):
        with pytest.raises(LookupError):
            extract_tagged_block(text.replace(f"```{lang} {tag}\n", f"```{lang} other\n"), tag,
                                 lang=lang)
        doubled = (text + "\n" + f"```{lang} {tag}\n"
                   + extract_tagged_block(text, tag, lang=lang) + "```\n")
        with pytest.raises(LookupError):
            extract_tagged_block(doubled, tag, lang=lang)
    # a tmp copy of the page with the probe call removed fails (X), naming it
    assert "readings = study.probe(" in block
    stripped = tmp_path / "page.md"
    stripped.write_text(text.replace("readings = study.probe(", "readings = (", 1))
    stripped_reason = _block_defect(extract_tagged_block(stripped.read_text(), PROBE_TAG))
    _p0(f"[WF-17 s18] (c-page) tags removed / doubled raise LookupError; probe call removed -> "
        f"{stripped_reason}")
    assert stripped_reason is not None and "'probe'" in stripped_reason
    for token in (".json", ".npz", ".s4p", ".xdmf", ".h5", "AED", "HFSS"):
        assert token not in text, token
    assert "(gel_sphere_heating_from_yaml.md)" in text

    # (In): the block's inputs are the gate's; Study.solve / Study.probe patched (no mesh)
    calls = {"solve": [], "probe": []}

    def fake_solve(self, *args, **kwargs):
        calls["solve"].append(self)
        return SimpleNamespace()

    def fake_probe(self, probes, cooling, **kwargs):
        calls["probe"].append((self, list(probes), cooling, kwargs))
        return SimpleNamespace(values=[], volumes=[])

    monkeypatch.setattr(Study, "solve", fake_solve)
    monkeypatch.setattr(Study, "probe", fake_probe)
    t0 = time.perf_counter()
    ns, _ = exec_probe_block(monkeypatch)
    st = ns["study"]
    want_probes = [WHOLE] + PROBES
    want_cooling = RobinCooling(K, H, T_INF, source_scale=G)
    assert len(calls["solve"]) == 1 and calls["solve"][0] is st
    assert len(calls["probe"]) == 1
    p_self, p_list, p_cool, p_kw = calls["probe"][0]
    _p0(f"[WF-17 s18] (In) solve calls {len(calls['solve'])}, probe calls "
        f"{len(calls['probe'])}; geometry {st.geometry!r}; material {st.material!r}; excitation "
        f"{st.excitation!r}; frequency {float(st.frequency_hz).hex()} (want "
        f"{FREQUENCY_HZ.hex()}); density {float(st.density_kg_m3).hex()} (want "
        f"{RHO_KG_M3.hex()}); degree {st.degree}; conductor {st.conductor!r}")
    _p0(f"[WF-17 s18] (In) probes {p_list!r} (want {want_probes!r}); cooling {p_cool!r} (want "
        f"{want_cooling!r}); kwargs {p_kw!r}; exec {time.perf_counter() - t0:.2f} s")
    assert p_self is st and p_kw == {}
    assert type(st.geometry) is SphereInBoxGeometry and st.geometry == GEOMETRY
    assert st.material == _material_057()
    assert float(st.material.sigma).hex() == SIGMA_HIGH.hex()
    assert type(st.excitation) is QuasiStaticSphereField
    assert st.excitation == QuasiStaticSphereField(e0=E0)
    assert float(st.frequency_hz).hex() == FREQUENCY_HZ.hex()
    assert float(st.density_kg_m3).hex() == RHO_KG_M3.hex()
    assert st.degree == 1
    assert st.conductor is None
    assert all(type(p) is ProbeBall for p in p_list)
    assert p_list == want_probes
    assert [_hex_ball(p) for p in p_list] == [_hex_ball(p) for p in want_probes]
    assert type(p_cool) is RobinCooling and p_cool == want_cooling
    assert float(p_cool.source_scale).hex() == float(G).hex()


def test_wf17_s18_bands_on_the_page():
    block = extract_tagged_block(_text(), BANDS_TAG, lang="text")
    _p0(f"[WF-17 s18] (E) bands block {block.strip().splitlines()} (rendered {_bands_lines()})")
    assert _lines_defect(block, _bands_lines()) is None, _lines_defect(block, _bands_lines())
    altered = block.replace(f"{WIDTH_BAND:.0e}", f"{WIDTH_BAND * 10:.0e}", 1)
    reason = _lines_defect(altered, _bands_lines())
    _p0(f"[WF-17 s18] (c-page) one exponent altered -> {reason}")
    assert reason is not None and "Displaced probe" in reason


def test_wf17_s18_docrefs_zero_dead():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(f"[WF-17 s18] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
          f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
          f"exit={counts['exit']} rc={rc}")
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- solve window (-n 2)


def test_wf17_s18_probe_block_reproduces_wf15_s4(monkeypatch):
    if _knob() is None:
        pytest.skip(f"{S18_ENV} unset: the tutorial solve window is environment-selected "
                    "(complex build, -n 2)")
    comm = MPI.COMM_WORLD
    assert IS_COMPLEX, "the gel-sphere EM solve needs the complex build"
    assert comm.size == 2, comm.size

    t0 = time.perf_counter()
    ns, _ = exec_probe_block(monkeypatch)
    t_block = time.perf_counter() - t0
    study, res, pr = ns["study"], ns["result"], ns["readings"]
    once = study.result is res

    # (A0) / (A) / (U) as test_wf15_study_probe.py computes them.
    sub, theta = pr.thermal.submesh, pr.theta
    cells = int(sub.topology.index_map(3).size_global)
    values = np.asarray(pr.values, dtype=np.float64)
    volumes = np.asarray(pr.volumes, dtype=np.float64)
    v_sub = float(np.real(_allsum(comm, fem.Constant(sub, default_scalar_type(1.0))
                                  * ufl.dx(domain=sub))))
    mean = float(np.real(_allsum(comm, theta * ufl.dx(domain=sub)))) / v_sub
    a_val = abs(values[0] - mean) / abs(mean)
    a_vol = abs(volumes[0] / v_sub - 1.0)
    want_probes = [WHOLE] + PROBES
    direct = np.array([ball_average(theta, p)[0] for p in want_probes])
    u_equal = bool(np.array_equal(direct, values))
    probes_equal = list(pr.probes) == want_probes

    # (c1) c2 displaced 0.1 mm outward.
    v_out = ball_average(theta, C2_OUT)[0]
    c1 = abs(v_out - values[3]) / abs(values[3])

    _p0(f"[WF-17 s18] (Once) {once}; probes = [WHOLE] + PROBES {probes_equal}")
    _p0(f"[WF-17 s18] (A0) sub cells {cells} (N_SUB {N_SUB}); (A) |P - int/V|/|int/V| = "
        f"{a_val:.3e}, |V_h/V_sub - 1| = {a_vol:.3e} (band {B_BAND:.0e}; RG1:136 {RG1_A}); "
        f"V_sub {v_sub:.10e}; whole-ball mean {mean:.10e} K")
    _p0(f"[WF-17 s18] (U) Study.probe values {values.tolist()!r}; direct {direct.tolist()!r}; "
        f"equal {u_equal}; volumes {volumes.tolist()!r}")
    _p0(f"[WF-17 s18] (c1) c2 +0.1 mm: |d|/theta = {c1:.4e} (ASSERTED > {WIDTH_BAND:.0e}; "
        f"RG1:143 {RG1_C1}); {c1 / WIDTH_BAND:.3e}x the band")
    _p0(f"[WF-17 s18] W-probe block (solve + probe) {t_block:.1f} s; max_imag_ratio "
        f"{pr.max_imag_ratio:.3e}; width {comm.size}")

    assert once, "(Once) a second EM solve happened"
    assert probes_equal
    assert cells == N_SUB, f"(A0) {cells}"
    assert a_val <= B_BAND, f"(A) value {a_val:.3e}"
    assert a_vol <= B_BAND, f"(A) volume {a_vol:.3e}"
    assert u_equal, "(U) Study.probe != ball_average"
    assert c1 > WIDTH_BAND, f"(c1) {c1:.3e} did not miss {WIDTH_BAND:.0e}"
