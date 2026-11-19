"""`WF-17` step 4 (§9 item 383, chain U39): "the tuned copper birdcage" is a test.

The committed page `docs/guide/tuned_copper_birdcage.md` carries one fenced
block tagged ```` ```python tuned-copper-birdcage ````. This module extracts it
verbatim and ``exec``s it. Asserted (brief: docs/planning/chunks/WF-17.md,
§9 item 383, and the 2026-10-11 18:00 review ruling beneath it):

(X) exactly one block per tag; exactly one ``Study(`` call, carrying ``degree``;
(K) the block configures H9′'s tuned copper F-human Study exactly (no solve);
    the ``Study.sar`` call asks for 10 g, ``("ccw",)`` and the lattice's c*;
(E) the expected-output block states item 382b's record at ``.4e``;
(R) at ``-n 8`` the block's SAR10g(c*) reproduces
    ``H9P_SAR10G_CSTAR_RECORD_W_PER_KG`` within ``REUSE_REPRODUCTION_RTOL``;
(D) ``docrefs --docs-root docs/guide --check-links`` reads dead == 0.

Negative controls (solve-free, asserted): (c1) untagged / doubled page raises;
(c2) ``esr_ohm=0.0`` mismatches exactly ``["excitation"]`` and the mirror-plane
azimuth 56.25° misses the centre check; (c3) step 1's planted dead link.

Disclosed (review ruling 2026-10-11 18:00, (3)): the no-solve stop patches the
module attribute ``fem_em_solver.study._study_solve_tuned_quadrature`` (looked
up at call time by ``_tuned_result``), because ``Study.sar`` on a
``TwoPortQuadratureDrive`` Study reaches neither ``Study.solve`` nor
``TwoPortQuadratureDrive._solve_port_study``. For the same reason this module
does not reuse step 1's ``exec_block`` (its stop patch is ``Study.solve`` and
its spy is ``Study.mean_sar``); it carries the same pattern on ``Study.sar``.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.study as study_module
from fem_em_solver import (
    FHumanBirdcageGeometry,
    GelledSalinePhantomMaterial,
    LeontovichWall,
    Study,
    TwoPortQuadratureDrive,
)
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.package.test_wf16_study_tuned_f_human_b1_sar import (
    H9P_CSTAR_LABEL,
    H9P_SAR10G_CSTAR_RECORD_W_PER_KG,
)
from tests.package.test_wf17_first_study_tutorial import (
    LINK_PASS_RE,
    _checker_stdout,
    _study_calls,
    extract_tagged_block,
)
from tests.unit.test_doc_reference_exit_codes import run_checker
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_wf7_f_human_sar10g_hotspot import _phantom_geometry
from tests.validation.test_wf10_f_human_copper_esr_root import ESR_OHM
from tests.validation.test_wf10_f_human_tuned_b1_sar import C_IN_F_HUMAN_PF
from tests.validation.test_wf10_f_human_tuned_sar10g import AZ17_RECORD_DEG, _lattice

REPO_ROOT = Path(__file__).resolve().parents[2]
PAGE = REPO_ROOT / "docs" / "guide" / "tuned_copper_birdcage.md"
TAG = "tuned-copper-birdcage"
OUTPUT_TAG = "tuned-copper-birdcage-output"
RECORD_RANK_WIDTH = 8  # the width H9′'s records were taken at
#: Pages in docs/guide at step 3 (39) plus this one.
MIN_GUIDE_PAGES = 40
#: z_c used by (K), before any mesh: item 382b's W1 printed the measured
#: phantom centre as "z_c 0.000000e+00 m" beside (S10); (R) re-derives z_c from
#: the live mesh with ``_phantom_geometry`` and re-asserts the centre.
Z_C_NO_MESH = 0.0


class _StopBeforeSolve(Exception):
    pass


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _p0(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(f"[WF-17 s4] {msg}", flush=True)


def _page_and_block():
    raw = PAGE.read_bytes()
    block = extract_tagged_block(raw.decode(), TAG)
    _p0(f"page sha256  = {_sha(raw)}")
    _p0(f"block sha256 = {_sha(block.encode())}")
    return raw.decode(), block


def exec_block(source: str, monkeypatch, *, stop_before_solve: bool):
    """Exec ``source`` as __main__ on every rank; return (recorded, stdout).

    A spy wraps ``Study.sar`` (calls the original; records the instance, the
    arguments and the returned dict). With ``stop_before_solve`` the tuned
    quadrature solve records the study and raises. Neither alters what the
    block computes.
    """
    recorded: dict = {}
    original_sar = Study.sar

    def spy_sar(self, *args, **kwargs):
        recorded["study"] = self
        recorded["args"] = args
        recorded["kwargs"] = dict(kwargs)
        out = original_sar(self, *args, **kwargs)
        recorded["sar"] = out
        return out

    monkeypatch.setattr(Study, "sar", spy_sar)
    if stop_before_solve:

        def stop_solve(study, comm=MPI.COMM_WORLD):
            recorded["solve_study"] = study
            raise _StopBeforeSolve

        monkeypatch.setattr(study_module, "_study_solve_tuned_quadrature", stop_solve)

    buf = io.StringIO()
    code = compile(source, f"docs/guide/tuned_copper_birdcage.md#{TAG}", "exec")
    with contextlib.redirect_stdout(buf):
        try:
            exec(code, {"__name__": "__main__"})
        except _StopBeforeSolve:
            if not stop_before_solve:
                raise
    captured = buf.getvalue()
    _p0(f"block stdout (rank 0): {captured!r}")
    return recorded, captured


def configuration_mismatches(study) -> list[str]:
    """(K)'s comparator: which of H9′'s Study fields differ, exactly."""
    expected = {
        "geometry": FHumanBirdcageGeometry(as_hole=True),
        "material": GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ),
        "excitation": TwoPortQuadratureDrive(C_IN_F_HUMAN_PF * 1e-12, esr_ohm=ESR_OHM),
        "conductor": LeontovichWall(),
        "frequency_hz": FREQUENCY_64_HZ,
        "degree": 1,
    }
    return [name for name, value in expected.items() if getattr(study, name) != value]


def _expected_c_star(z_c):
    rows = {lab: c for lab, *_x, c in _lattice(AZ17_RECORD_DEG, z_c)}
    return np.asarray(rows[H9P_CSTAR_LABEL], dtype=float)


def sar_call_mismatches(recorded, z_c) -> list[str]:
    """(K)'s ``Study.sar`` comparator: mass, drives, the single centre c*."""
    kw = dict(recorded["kwargs"])
    if recorded["args"]:
        kw.setdefault("mass_g", recorded["args"][0])
    out = []
    if kw.get("mass_g") != 10.0:
        out.append("mass_g")
    if tuple(kw.get("drives") or ()) != ("ccw",):
        out.append("drives")
    centres = list(kw.get("centres") or [])
    c_star = _expected_c_star(z_c)
    rel = (
        float(np.linalg.norm(np.asarray(centres[0], dtype=float) - c_star) / np.linalg.norm(c_star))
        if len(centres) == 1
        else float("inf")
    )
    _p0(f"(K) spied centres {centres}; c* {tuple(c_star)}; rel {rel:.3e} (band {IDENTITY_BAND:g})")
    if not rel <= IDENTITY_BAND:
        out.append("centre")
    return out


def _assert_extraction(block):
    assert block.strip(), "(X) the block is empty"
    calls = _study_calls(block)
    assert len(calls) == 1, f"(X) {len(calls)} Study( calls in the block, need exactly 1"
    keywords = [kw.arg for kw in calls[0].keywords]
    _p0(f"(X) Study keywords in the block: {keywords}")
    assert "degree" in keywords, "(X) the block must set degree explicitly"
    return calls[0]


# ---------------------------------------------------------------- W0 (smoke)


def test_wf17_s4_extracts_exactly_one_block():
    text, block = _page_and_block()
    _assert_extraction(block)
    extract_tagged_block(text, OUTPUT_TAG, lang="text")


def test_wf17_s4_block_configures_the_h9p_study(monkeypatch):
    _, block = _page_and_block()
    _assert_extraction(block)
    recorded, _ = exec_block(block, monkeypatch, stop_before_solve=True)
    assert "solve_study" in recorded, "(K) the block never reached the tuned solve"
    assert recorded["solve_study"] is recorded["study"]
    mismatches = configuration_mismatches(recorded["study"])
    _p0(f"(K) Study mismatches: {mismatches}")
    assert mismatches == []
    call = sar_call_mismatches(recorded, Z_C_NO_MESH)
    _p0(f"(K) sar-call mismatches: {call}")
    assert call == []


def test_wf17_s4_planted_page_and_parameter_controls(monkeypatch):
    text, block = _page_and_block()
    # (c1) zero tags and two tags both raise.
    untagged = text.replace(f"```python {TAG}\n", "```python\n")
    with pytest.raises(LookupError):
        extract_tagged_block(untagged, TAG)
    doubled = text + f"\n```python {TAG}\n{block}```\n"
    with pytest.raises(LookupError):
        extract_tagged_block(doubled, TAG)
    # (c2a) ESR 0.2 -> 0.0 is caught before any solve, and only as excitation.
    assert block.count("esr_ohm=0.2") == 1, "(c2) the ESR literal must occur exactly once"
    recorded, _ = exec_block(
        block.replace("esr_ohm=0.2", "esr_ohm=0.0"), monkeypatch, stop_before_solve=True
    )
    mismatches = configuration_mismatches(recorded["study"])
    _p0(f"(c2a) mismatches on esr_ohm=0.0: {mismatches}")
    assert mismatches == ["excitation"]
    assert sar_call_mismatches(recorded, Z_C_NO_MESH) == []
    # (c2b) the feed-bisector mirror azimuth misses the centre check only.
    assert block.count("np.radians(11.25)") == 1, "(c2) the azimuth literal must occur once"
    monkeypatch.undo()
    recorded, _ = exec_block(
        block.replace("np.radians(11.25)", "np.radians(56.25)"), monkeypatch, stop_before_solve=True
    )
    assert configuration_mismatches(recorded["study"]) == []
    call = sar_call_mismatches(recorded, Z_C_NO_MESH)
    _p0(f"(c2b) sar-call mismatches on az 56.25: {call}")
    assert call == ["centre"]


def test_wf17_s4_expected_output_states_the_record():
    text, _ = _page_and_block()
    out = extract_tagged_block(text, OUTPUT_TAG, lang="text")
    expected = f"SAR10g(c*) = {H9P_SAR10G_CSTAR_RECORD_W_PER_KG:.4e} W/kg\n"
    _p0(f"(E) page output {out!r}; expected {expected!r}")
    assert out == expected


def test_wf17_s4_docrefs_guide_has_zero_dead_links():
    if MPI.COMM_WORLD.rank != 0:
        return
    rc, counts = run_checker("--docs-root", "docs/guide", "--check-links")
    link = LINK_PASS_RE.search(_checker_stdout("--docs-root", "docs/guide", "--check-links"))
    assert link is not None, "(D) no 'Link pass:' line"
    print(
        f"[WF-17 s4] (D) dead={counts['dead']} links={link.group(1)} pages={link.group(2)} "
        f"link-dead={link.group(3)} guide={counts['guide']} stale={counts['stale']} "
        f"exit={counts['exit']} rc={rc}"
    )
    assert counts["dead"] == 0
    assert int(link.group(3)) == 0
    assert int(link.group(2)) >= MIN_GUIDE_PAGES


# ------------------------------------------------------------- W1 (-n 8)


@complex_only
def test_wf17_s4_block_reproduces_the_h9p_record(monkeypatch):
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, f"record is -n {RECORD_RANK_WIDTH}; got {comm.size}"
    _, block = _page_and_block()
    _assert_extraction(block)
    comm.Barrier()
    t0 = time.perf_counter()
    recorded, captured = exec_block(block, monkeypatch, stop_before_solve=False)
    comm.Barrier()
    seconds = time.perf_counter() - t0
    study = recorded["study"]
    # (K) on the live instance.
    assert configuration_mismatches(study) == []
    res = study.result
    centre, _radius, _hh, _lo, _hi = _phantom_geometry(
        res.mesh, res.cell_tags, PHANTOM_CELL_TAG, comm
    )
    z_c = float(centre[2])
    _p0(f"(K) live z_c {z_c:.6e} m; cells {res.mesh.topology.index_map(3).size_global}")
    assert sar_call_mismatches(recorded, z_c) == []
    # (R)
    rows = recorded["sar"]["ccw"]
    assert len(rows) == 1
    sar = float(rows[0]["averaged_sar_w_per_kg"])
    rel = abs(sar - H9P_SAR10G_CSTAR_RECORD_W_PER_KG) / H9P_SAR10G_CSTAR_RECORD_W_PER_KG
    _p0(
        f"(R) SAR10g(c*) = {sar!r} ({sar:.9e}) W/kg; record "
        f"{H9P_SAR10G_CSTAR_RECORD_W_PER_KG:.9e}; rel {rel:.3e} (band REUSE_REPRODUCTION_RTOL "
        f"{REUSE_REPRODUCTION_RTOL:g}); mass {rows[0]['mass_kg'] * 1e3:.6f} g"
    )
    _p0(f"block exec incl. build + node pair + one 10 g ball: {seconds:.1f} s at -n {comm.size}")
    assert rel <= REUSE_REPRODUCTION_RTOL
    if comm.rank == 0:
        assert f"SAR10g(c*) = {sar:.4e} W/kg" in captured
