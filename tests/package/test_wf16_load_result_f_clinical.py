"""`WF-16` step 16 (§9 item 696, §10 U104): ``Study.load_result`` and ``tune`` on
the saved F-clinical 64×64, the package's first read of the production fixture.

``FClinicalBirdcageGeometry`` has no ``build`` and a solve stays ``xl`` tier; a
tune on a matrix the ``xl`` tier already produced needs no solve.  Fixture: C3's
copper matrix ``logs/geo38/c3-cu-s64-d1.npz`` (gitignored; degree 1, 64 MHz).

Gates (asserted; every band imported, none set here):
  (CM)   the loaded ``C_mode`` equals C10's ``c_mode_f`` within ``REDUCTION_BAND``.
  (CM-r) ``f"{C * 1e12:.4f}"`` equals ``C10_C_MODE_PF_RECORD``, C10's record
         (`20261024T093056Z_GEO-38-C10-W1.log`), the record a clone can read.
  (S0)   class-level spies on ``Study.solve`` / ``Study._solve_ports`` record 0
         calls on the fresh Study.
  (RT)   max |ΔS| against the npz within ``EXACT_IDENTITY_RTOL``.
Negative control (asserted): ports 1 and 2 swapped in the data, names kept, misses
the record by more than ``LEG_D0_REPRODUCTION_BAND`` (a non-capacitive read is a
miss too); the miss is printed.

Opt-in: ``FEM_EM_WF16_S16=1`` (unset => the saved-matrix node skips; any other
value raises; the npz absent under the knob fails, never skips, never re-solves).
The no-mesh node (W0) always runs.
Scope: a saved 1.5 T 64×64 tunes with no solve.  Not claimed: the 3 T twin, the
YAML route, a solve, any field output, a tuning statement beyond C10's ``C_mode``.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import study as study_mod
from fem_em_solver.ports.circuit import z_to_s
from fem_em_solver.study import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    SphereInBoxGeometry,
    Study,
)
from tests.complex_mode import complex_only
from tests.package.test_wf13_estimate_f_clinical import _study, no_build  # noqa: F401
from tests.package.test_wf16_load_result_f_human import OUT_DIR, _swap12_n, _write_named
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_geo38_f_clinical_c_mode import (
    C10_JSON,
    WORKSPACE,
    _block_circulant,
    _load,
    _sym_blocks,
)
from tests.validation.test_geo38_f_clinical_matrix import SAVE_PATH
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

S16_ENV = "FEM_EM_WF16_S16"
N_SLOTS = 32
#: C10's record at ``.4f`` pF (`20261024T093056Z_GEO-38-C10-W1.log`).
C10_C_MODE_PF_RECORD = "195.0265"


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S16_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S16_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _expect_ref(st, tag):
    """The loaded F-clinical refuses every field method, whichever refusal fires
    first (printed)."""
    for name, call in (
        ("sar", lambda: st.sar(centres=[(0.0, 0.0, 0.0)])),
        ("b1_plus", lambda: st.b1_plus()),
        ("temperature", lambda: st.temperature(None)),
    ):
        with pytest.raises((TypeError, ValueError)) as ei:
            call()
        _p0(
            f"[WF-16 s16 {tag} Ref] {name} on the loaded F-clinical: "
            f"{type(ei.value).__name__}: {ei.value}"
        )


@complex_only
def test_wf16_s16_w0_knob_and_synthetic(monkeypatch, no_build):  # noqa: F811
    """W0, no mesh: knob parse; a synthetic 64-port under P33…P96 loads and tunes to
    its closed form; a 63-port file is refused; the unloaded tune raises the solve's
    refusal."""
    comm = MPI.COMM_WORLD
    for bad in ("2", "yes"):
        monkeypatch.setenv(S16_ENV, bad)
        with pytest.raises(ValueError, match=S16_ENV):
            _knob()
    monkeypatch.delenv(S16_ENV, raising=False)
    assert _knob() is False

    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("W0 must not build a mesh")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry):
        monkeypatch.setattr(cls, "build", stub)

    st = _study()
    c_ref = 195.0e-12
    z = _block_circulant(_sym_blocks(N_SLOTS, st.omega, c_ref), N_SLOTS)
    z0 = float(st.excitation.z0_ohm)
    s = z_to_s(z, z0)
    names = [f"P{i}" for i in range(N_SLOTS + 1, 3 * N_SLOTS + 1)]
    tag = f"r{comm.size}"
    ok = _write_named(OUT_DIR / f"s16_synthetic_{tag}.s64p", st.frequency_hz, s, names)
    res = st.load_result(ok)
    assert isinstance(res, study_mod.LoadedPortResult) and st.result is res
    assert len(res.port_ids) == 64 and res.s_matrix.shape == (64, 64)
    assert list(res.port_ids) == names
    rt = float(np.max(np.abs(res.s_matrix - s)))
    t1 = st.tune(1, n_slots=N_SLOTS, comm=comm)
    e1 = abs(t1.c_mode_farad / c_ref - 1.0)
    _p0(
        f"\n[WF-16 s16 W0] synthetic 64-port: ports {len(res.port_ids)}, max|dS| {rt:.3e}; "
        f"C_1 {t1.c_mode_farad:.12e} F vs closed {c_ref:.12e} F, rel {e1:.3e} "
        f"(EXACT_IDENTITY_RTOL {EXACT_IDENTITY_RTOL:.1e})"
    )
    assert e1 <= EXACT_IDENTITY_RTOL
    _expect_ref(st, "W0")

    st63 = _study()
    p63 = _write_named(
        OUT_DIR / f"s16_synthetic63_{tag}.s63p", st63.frequency_hz, s[:63, :63], names[:63]
    )
    with pytest.raises(ValueError, match="ports") as ei:
        st63.load_result(p63)
    assert st63.result is None
    _p0(f"[WF-16 s16 W0] 63-port file refused: ValueError: {ei.value}")

    unloaded = _study()
    with pytest.raises(ValueError, match="xl tier only") as ei:
        unloaded.tune(1, comm=comm)
    assert unloaded.result is None
    _p0(f"[WF-16 s16 W0] unloaded tune refused: ValueError: {ei.value}")
    assert calls["n"] == 0
    assert no_build == []


@complex_only
def test_wf16_s16_tune_on_saved_f_clinical():
    if not _knob():
        pytest.skip(f"{S16_ENV} unset: the saved C3 matrix lives in the gitignored logs/")
    comm = MPI.COMM_WORLD
    npz = WORKSPACE / SAVE_PATH
    assert Path(npz).exists(), f"{npz} is absent: restore it, never re-solve"
    assert Path(C10_JSON).exists(), f"{C10_JSON} is absent"
    data = _load(npz)
    c10 = json.loads(Path(C10_JSON).read_text())
    order = data["order"]
    s64 = np.asarray(data["s64"], dtype=np.complex128)
    fresh = _study()
    fhz = fresh.frequency_hz
    assert s64.shape == (64, 64) and len(order) == 64
    assert float(c10["z0_ohm"]) == float(fresh.excitation.z0_ohm), c10["z0_ohm"]
    assert float(c10["freq_hz"]) == float(fhz), c10["freq_hz"]
    path = _write_named(OUT_DIR / f"s16_f_clinical_cu_d1_n{comm.size}.s64p", fhz, s64, order)

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
    _p0(
        f"\n[WF-16 s16 S0] class-level spies: solve {calls['solve']}, "
        f"_solve_ports {calls['_solve_ports']}"
    )
    _expect_ref(fresh, "W1")

    rt = float(np.max(np.abs(loaded.s_matrix - s64)))
    cl = tl.c_mode_farad
    g_cm = abs(cl / float(c10["c_mode_f"]) - 1.0)
    pf = f"{cl * 1e12:.4f}"
    _p0(
        f"[WF-16 s16 RT] max|dS| against the npz {rt:.3e} (EXACT_IDENTITY_RTOL "
        f"{EXACT_IDENTITY_RTOL:.1e}); port order kept {list(loaded.port_ids) == order}"
    )
    _p0(
        f"[WF-16 s16 CM] loaded C_mode {cl:.12e} F vs C10 {float(c10['c_mode_f']):.12e} F: rel "
        f"{g_cm:.3e} (REDUCTION_BAND {REDUCTION_BAND:.1e}); off-block {tl.off_diagonal_residual:.6e}"
    )
    _p0(f"[WF-16 s16 CM-r] {pf} pF vs record {C10_C_MODE_PF_RECORD}")

    # Negative control: ports 1 <-> 2 swapped in the data, names kept.
    sw = _write_named(
        OUT_DIR / f"s16_swap12_n{comm.size}.s64p",
        fhz,
        _swap12_n(loaded.s_matrix),
        list(loaded.port_ids),
    )
    st_sw = _study()
    st_sw.load_result(sw)
    try:
        cw = st_sw.tune(1, comm=comm).c_mode_farad
        miss = abs(cw / float(c10["c_mode_f"]) - 1.0)
        _p0(
            f"[WF-16 s16 Sw] ports 1<->2 swapped: C {cw:.12e} F, miss {miss:.6e} vs "
            f"LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.1e}"
        )
    except ValueError as exc:  # a non-capacitive read is a miss too; printed
        miss = float("inf")
        _p0(f"[WF-16 s16 Sw] ports 1<->2 swapped: tune raised {exc} (miss = inf)")

    assert calls == {"solve": 0, "_solve_ports": 0}
    assert list(loaded.port_ids) == order
    assert rt <= EXACT_IDENTITY_RTOL
    assert g_cm <= REDUCTION_BAND
    assert pf == C10_C_MODE_PF_RECORD
    assert miss > LEG_D0_REPRODUCTION_BAND
