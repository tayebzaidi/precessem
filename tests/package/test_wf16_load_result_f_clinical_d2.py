"""`WF-16` step 20 (§9 item 867, §10 U129): ``Study.load_result`` and ``tune(1)``
on the saved **degree-2** F-clinical 1.5 T 64×64, no solve.

Fixture: ``logs/geo38/c3-cu-s64-d2.npz`` (gitignored; degree 2, 64 MHz), converted
unpermuted as U104 did; comparand C18's ``logs/geo38/c18-c-mode-d2.json``.

Gates (asserted; every band imported, none set here):
  (CM)   the loaded ``C_mode`` equals C18's ``c_mode_f`` within ``REDUCTION_BAND``.
  (CM-r) ``f"{C * 1e12:.4f}"`` equals ``C18_C_MODE_PF_RECORD``
         (`20261107T182319Z_GEO-38-C18-W1.log:74`).
  (S0)   class-level spies on ``Study.solve`` / ``Study._solve_ports``: 0 calls.
  (RT)   max |ΔS| against the npz within ``EXACT_IDENTITY_RTOL``, port order kept.
  (F)    file frequency 64e6 equals the Study's; ``z0_ohm`` equals the JSON's.
Controls:
  (D1)   the degree-1 file loaded into a fresh degree-2 Study misses C18 by more
         than ``REDUCTION_BAND`` -- asserted because C18's W1 printed the d1 -> d2
         move -4.9335 % (`20261107T182319Z_GEO-38-C18-W1.log:78`), above 1e-3; its
         ``.4f`` equals C10's record.  ``load_result`` does not refuse a degree
         mismatch, so this is a numeric miss, weaker than U117's refusal.
  (Sw)   ports 1 <-> 2 swapped misses by more than ``LEG_D0_REPRODUCTION_BAND``.

Opt-in: ``FEM_EM_WF16_S20=1`` (unset => skip; any other value raises; an absent
input under the knob fails, never skips, never re-solves).
Scope: the circuit layer on a saved degree-2 file.  Not claimed: a root, a tuning
statement, any solve.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.study import Study

from tests.complex_mode import complex_only
from tests.package.test_wf13_estimate_f_clinical import _study, no_build  # noqa: F401
from tests.package.test_wf16_load_result_f_clinical import C10_C_MODE_PF_RECORD
from tests.package.test_wf16_load_result_f_human import OUT_DIR, _swap12_n, _write_named
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_geo38_f_clinical_c_mode import WORKSPACE, _load
from tests.validation.test_geo38_f_clinical_matrix import SAVE_PATH, SAVE_PATH_D2
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

S20_ENV = "FEM_EM_WF16_S20"
C18_JSON = WORKSPACE / "logs/geo38/c18-c-mode-d2.json"
#: C18's record at ``.4f`` pF (`20261107T182319Z_GEO-38-C18-W1.log:74`).
C18_C_MODE_PF_RECORD = "185.4048"
#: C18's printed d1 -> d2 relative move of C_mode (`…C18-W1.log:78`), the (D1) licence.
C18_D1_D2_MOVE = -4.9335e-2


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S20_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S20_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _knob_parse(monkeypatch):
    for bad in ("2", "yes"):
        monkeypatch.setenv(S20_ENV, bad)
        with pytest.raises(ValueError, match=S20_ENV):
            _knob()
    monkeypatch.delenv(S20_ENV, raising=False)
    assert _knob() is False


@complex_only
def test_wf16_s20_w0_knob_and_degree(monkeypatch, no_build):  # noqa: F811
    """No mesh: knob parse; ``_study(degree=2)`` carries degree 2, and the
    default ``_study()`` degree 1 (the trap the brief names)."""
    _knob_parse(monkeypatch)
    assert _study(degree=2).degree == 2
    assert _study().degree == 1
    assert no_build == []


def _convert(npz, tag, comm):
    data = _load(npz)
    order = list(data["order"])
    s64 = np.asarray(data["s64"], dtype=np.complex128)
    assert s64.shape == (64, 64) and len(order) == 64
    st = _study(degree=2)
    path = _write_named(OUT_DIR / f"s20_{tag}_n{comm.size}.s64p", st.frequency_hz, s64, order)
    return st, path, s64, order


@complex_only
def test_wf16_s20_tune_on_saved_f_clinical_d2():
    if not _knob():
        pytest.skip(f"{S20_ENV} unset: the saved degree-2 matrix lives in the gitignored logs/")
    comm = MPI.COMM_WORLD
    npz = WORKSPACE / SAVE_PATH_D2
    assert Path(npz).exists(), f"{npz} is absent: restore it, never re-solve"
    assert Path(C18_JSON).exists(), f"{C18_JSON} is absent"
    c18 = json.loads(Path(C18_JSON).read_text())
    c_ref = float(c18["c_mode_f"])

    fresh, path, s64, order = _convert(npz, "f_clinical_cu_d2", comm)
    fhz = fresh.frequency_hz
    _p0(f"\n[WF-16 s20 F] Study degree {fresh.degree}; Study f {fhz:.6e} Hz, C18 f "
        f"{float(c18['freq_hz']):.6e} Hz; z0 {float(fresh.excitation.z0_ohm)} vs C18 "
        f"{float(c18['z0_ohm'])}")
    assert fresh.degree == 2
    assert float(fhz) == 64e6
    assert float(c18["freq_hz"]) == float(fhz)
    assert float(c18["z0_ohm"]) == float(fresh.excitation.z0_ohm)

    calls = {"solve": 0, "_solve_ports": 0}

    def _spy(name):
        def _boom(self, *a, **k):
            calls[name] += 1
            raise AssertionError(f"(S0) the loaded Study called {name}")
        return _boom

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(Study, "solve", _spy("solve"))
        mp.setattr(Study, "_solve_ports", _spy("_solve_ports"))
        loaded = fresh.load_result(path)
        tl = fresh.tune(1, comm=comm)

        # (D1): the degree-1 file into a fresh degree-2 Study.
        d1_npz = WORKSPACE / SAVE_PATH
        assert Path(d1_npz).exists(), f"{d1_npz} is absent: restore it, never re-solve"
        st_d1, p_d1, _, _ = _convert(d1_npz, "f_clinical_cu_d1_into_d2", comm)
        st_d1.load_result(p_d1)
        c_d1 = st_d1.tune(1, comm=comm).c_mode_farad
    _p0(f"[WF-16 s20 S0] class-level spies: solve {calls['solve']}, "
        f"_solve_ports {calls['_solve_ports']}")

    rt = float(np.max(np.abs(loaded.s_matrix - s64)))
    cl = tl.c_mode_farad
    g_cm = abs(cl / c_ref - 1.0)
    pf = f"{cl * 1e12:.4f}"
    _p0(f"[WF-16 s20 RT] max|dS| against the npz {rt:.3e} (EXACT_IDENTITY_RTOL "
        f"{EXACT_IDENTITY_RTOL:.1e}); port order kept {list(loaded.port_ids) == order}")
    _p0(f"[WF-16 s20 CM] loaded C_mode {cl:.12e} F vs C18 {c_ref:.12e} F: rel "
        f"{g_cm:.3e} (REDUCTION_BAND {REDUCTION_BAND:.1e}); off-block {tl.off_diagonal_residual:.6e}")
    _p0(f"[WF-16 s20 CM-r] {pf} pF vs record {C18_C_MODE_PF_RECORD}")

    miss_d1 = abs(c_d1 / c_ref - 1.0)
    pf_d1 = f"{c_d1 * 1e12:.4f}"
    _p0(f"[WF-16 s20 D1] degree-1 file in a degree-{st_d1.degree} Study: C {pf_d1} pF "
        f"(C10 record {C10_C_MODE_PF_RECORD}), miss vs C18 {miss_d1:.6e} vs REDUCTION_BAND "
        f"{REDUCTION_BAND:.1e} (ASSERTED: C18 W1 printed move {C18_D1_D2_MOVE:+.4%})")

    sw = _write_named(OUT_DIR / f"s20_swap12_n{comm.size}.s64p", fhz, _swap12_n(loaded.s_matrix),
                      list(loaded.port_ids))
    st_sw = _study(degree=2)
    st_sw.load_result(sw)
    try:
        cw = st_sw.tune(1, comm=comm).c_mode_farad
        miss = abs(cw / c_ref - 1.0)
        _p0(f"[WF-16 s20 Sw] ports 1<->2 swapped: C {cw:.12e} F, miss {miss:.6e} vs "
            f"LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.1e}")
    except ValueError as exc:  # a non-capacitive read is a miss too; printed
        miss = float("inf")
        _p0(f"[WF-16 s20 Sw] ports 1<->2 swapped: tune raised {exc} (miss = inf)")

    assert calls == {"solve": 0, "_solve_ports": 0}
    assert list(loaded.port_ids) == order
    assert rt <= EXACT_IDENTITY_RTOL
    assert g_cm <= REDUCTION_BAND
    assert pf == C18_C_MODE_PF_RECORD
    assert abs(C18_D1_D2_MOVE) > REDUCTION_BAND
    assert miss_d1 > REDUCTION_BAND
    assert pf_d1 == C10_C_MODE_PF_RECORD
    assert miss > LEG_D0_REPRODUCTION_BAND
