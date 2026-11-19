"""`WF-16` step 19 (§9 item 781, §10 U117): ``Study.load_result`` and ``tune`` on
the F-clinical-3T's saved 64×64 at 128 MHz, with no solve; plus `WF-1`'s YAML half.

Fixture: entry 47's copper matrix ``logs/geo38/c11-cu-s64-d1-128mhz.npz``
(gitignored; degree 1, 128 MHz), read as C15 reads it: ``s64`` in the npz's own
``order`` (P33…P96), no port swap (C15, ``test_geo38_c15_3t_c_mode.py``, uses the
npz unpermuted; ``_swap12_n`` is U104's negative control, not a conversion step).

The 3 T selector: ``FClinicalBirdcageGeometry(field_strength_t=3.0)`` (a frozen
dataclass with no validation, so 3.0 is admitted; printed).  ``load_result`` reads
only the file's frequency and port count, so the selector changes nothing on
this path.

Gates (asserted; every band imported, none set here):
  (CM)   the loaded ``C_mode`` equals C15's ``c_mode_f`` within ``REDUCTION_BAND``,
         and ``f"{C * 1e12:.4f}"`` equals ``C15_C_MODE_PF_RECORD``.
  (S0)   class-level spies on ``Study.solve`` / ``Study._solve_ports`` record 0 calls.
  (F)    the file's frequency is 128e6 Hz.
  (RT)   max |ΔS| against the npz within ``EXACT_IDENTITY_RTOL``.
  (Y=)   the 1.5 T YAML with ``frequency_hz: 128000000.0`` and ``load_result:`` gives
         the same ``C_mode`` by ``float.hex``, no mesh and no solve reached.
Control (asserted, a refusal): the 1.5 T file (C3's 64 MHz npz, converted) loaded
into the 128 MHz Study raises the frequency ``ValueError``; nothing is stored.

Opt-in: ``FEM_EM_WF16_S19=1`` (unset => the saved-matrix node skips; any other
value raises; an absent input under the knob fails, never skips, never re-solves).
The no-mesh node always runs.  No S entry is printed.
Scope: the circuit layer on a saved file.  Not claimed: a 3 T root, any solve.
"""

from __future__ import annotations

import copy
import json
import os
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import LumpedPortDrive, run_study_file
from fem_em_solver import study as study_mod
from fem_em_solver.ports.touchstone import load_touchstone
from fem_em_solver.study import FClinicalBirdcageGeometry, Study
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import _fresh_dir, _mapping
from tests.package.test_wf1_yaml_load_result import _hex, _no_mesh, _write_yaml
from tests.package.test_wf13_estimate_f_clinical import FC, _study
from tests.package.test_wf16_load_result_f_human import OUT_DIR, _write_named
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_geo38_c15_3t_c_mode import C15_JSON
from tests.validation.test_geo38_f_clinical_3t_matrix import C11_SAVE_PATH
from tests.validation.test_geo38_f_clinical_c_mode import WORKSPACE, _load
from tests.validation.test_geo38_f_clinical_matrix import SAVE_PATH
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

S19_ENV = "FEM_EM_WF16_S19"
N_SLOTS = 32
FREQUENCY_3T_HZ = 128e6
#: C15's record at ``.4f`` pF (`20261027T093207Z_GEO-38-C15-W1.log:76`).
C15_C_MODE_PF_RECORD = "43.9719"


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S19_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S19_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _study_3t():
    base = _study()
    st = Study(
        FClinicalBirdcageGeometry(field_strength_t=3.0),
        base.material,
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_3T_HZ,
        degree=1,
    )
    return st


@complex_only
def test_wf16_s19_w0_knob_and_selector(monkeypatch):
    """No mesh: knob parse; the 3 T Study constructs with its selector and frequency."""
    for bad in ("2", "yes"):
        monkeypatch.setenv(S19_ENV, bad)
        with pytest.raises(ValueError, match=S19_ENV):
            _knob()
    monkeypatch.delenv(S19_ENV, raising=False)
    assert _knob() is False
    st = _study_3t()
    _p0(
        f"\n[WF-16 s19 W0] selector: FClinicalBirdcageGeometry(field_strength_t="
        f"{st.geometry.field_strength_t}), frequency_hz {st.frequency_hz!r}, degree {st.degree}"
    )
    assert st.geometry.field_strength_t == 3.0
    assert st.frequency_hz == FREQUENCY_3T_HZ
    assert st.result is None


@complex_only
def test_wf16_s19_tune_on_saved_f_clinical_3t(capsys):
    if not _knob():
        pytest.skip(f"{S19_ENV} unset: the saved C11 matrix lives in the gitignored logs/")
    comm = MPI.COMM_WORLD
    npz = WORKSPACE / C11_SAVE_PATH
    npz15 = WORKSPACE / SAVE_PATH
    for p in (npz, npz15, C15_JSON):
        assert Path(p).exists(), f"{p} is absent: restore it, never re-solve"
    data = _load(npz)
    c15 = json.loads(Path(C15_JSON).read_text())
    c_ref = float(c15["c_mode_f"])
    order = data["order"]
    s64 = np.asarray(data["s64"], dtype=np.complex128)
    fresh = _study_3t()
    fhz = fresh.frequency_hz
    assert s64.shape == (64, 64) and order == [f"P{i}" for i in range(33, 97)]
    assert float(c15["z0_ohm"]) == float(fresh.excitation.z0_ohm), c15["z0_ohm"]
    assert float(c15["freq_hz"]) == float(fhz), c15["freq_hz"]
    assert int(c15["n_slots"]) == N_SLOTS
    path = _write_named(OUT_DIR / f"s19_f_clinical_3t_cu_d1_n{comm.size}.s64p", fhz, s64, order)
    f_file = np.asarray(load_touchstone(path)[0], dtype=np.float64)

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
        tl = fresh.tune(k=1, n_slots=N_SLOTS, comm=comm)

    rt = float(np.max(np.abs(loaded.s_matrix - s64)))
    cl = tl.c_mode_farad
    g_cm = abs(cl / c_ref - 1.0)
    pf = f"{cl * 1e12:.4f}"

    # Control: the 1.5 T file (C3's 64 MHz matrix) into the 128 MHz Study is refused.
    d15 = _load(npz15)
    p15 = _write_named(
        OUT_DIR / f"s19_f_clinical_15t_cu_d1_n{comm.size}.s64p",
        float(_study().frequency_hz),
        np.asarray(d15["s64"], dtype=np.complex128),
        d15["order"],
    )
    st_c = _study_3t()
    with pytest.raises(ValueError, match="file frequency") as ei:
        st_c.load_result(p15)
    ctrl_msg = str(ei.value)
    ctrl_none = st_c.result is None

    # YAML half: the 1.5 T file with frequency_hz 128 MHz and load_result:.
    d = OUT_DIR / f"s19_yaml_n{comm.size}"
    _fresh_dir(d, comm)
    tgt = _write_named(d / "c11_cu_d1_128mhz.s64p", fhz, s64, order)
    m = copy.deepcopy(_mapping(FC))
    m["study"]["frequency_hz"] = 128000000.0
    m["tune"] = {"k": 1}
    m["load_result"] = tgt.name
    yml = _write_yaml(d / "loaded_128mhz.yaml", m)
    ycalls = {"n": 0}
    with pytest.MonkeyPatch.context() as mp:
        _no_mesh(mp, ycalls)
        run = run_study_file(yml, output_dir=d / "out", comm=comm)
    n_y = comm.allreduce(ycalls["n"], op=MPI.SUM)
    cy = run.summary["tune"]["c_mode_farad"]
    with open(d / "out" / "summary.json") as fh:
        cy_js = float(json.load(fh)["tune"]["c_mode_farad"])

    with capsys.disabled():
        _p0(
            f"\n[WF-16 s19 sel] FClinicalBirdcageGeometry(field_strength_t="
            f"{fresh.geometry.field_strength_t}), frequency_hz {fhz!r}"
        )
        _p0(
            f"[WF-16 s19 S0] class-level spies: solve {calls['solve']}, "
            f"_solve_ports {calls['_solve_ports']}"
        )
        _p0(f"[WF-16 s19 F] file frequency {f_file.tolist()} Hz; loaded {loaded.frequency_hz!r}")
        _p0(
            f"[WF-16 s19 RT] max|dS| against the npz {rt:.3e} (EXACT_IDENTITY_RTOL "
            f"{EXACT_IDENTITY_RTOL:.1e}); port order kept {list(loaded.port_ids) == order}; "
            f"port swap applied: none (C15 reads the npz unpermuted)"
        )
        _p0(
            f"[WF-16 s19 CM] loaded C_mode {cl:.12e} F vs C15 {c_ref:.12e} F: rel "
            f"{g_cm:.3e} (REDUCTION_BAND {REDUCTION_BAND:.1e}); off-block "
            f"{tl.off_diagonal_residual:.6e} (C15 {float(c15['off_block']):.6e})"
        )
        _p0(f"[WF-16 s19 CM-r] {pf} pF vs record {C15_C_MODE_PF_RECORD}")
        _p0(
            f"[WF-16 s19 Ctl] 1.5 T file into the 128 MHz Study: ValueError: {ctrl_msg}; "
            f"result None {ctrl_none}"
        )
        _p0(
            f"[WF-16 s19 Y=] YAML C {_hex(cy)} vs api {_hex(cl)}; summary.json {_hex(cy_js)}; "
            f"mesh/solve calls reached {n_y}; result loaded "
            f"{isinstance(run.study.result, study_mod.LoadedPortResult)}"
        )

    assert calls == {"solve": 0, "_solve_ports": 0}
    assert f_file.shape == (1,) and float(f_file[0]) == FREQUENCY_3T_HZ
    assert float(loaded.frequency_hz) == FREQUENCY_3T_HZ
    assert list(loaded.port_ids) == order
    assert rt <= EXACT_IDENTITY_RTOL
    assert g_cm <= REDUCTION_BAND
    assert pf == C15_C_MODE_PF_RECORD
    assert ctrl_none
    assert n_y == 0
    assert isinstance(run.study.result, study_mod.LoadedPortResult)
    assert _hex(cy) == _hex(cl) and _hex(cy_js) == _hex(cl)
