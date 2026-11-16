"""`WF-16` step 17 (§9 item 722, §10 U108): ``Study.load_result`` on a multi-frequency file.

A multi-frequency Touchstone file now loads on a four-port ``Study``; ``tune(k, at_hz=f)``
selects the file's nearest frequency and uses that frequency's ω.  The file is U71's sweep
output, ``logs/wf16/wf16_sweep_f_small_p1_n2.s4p`` (gitignored; absent ⇒ fail, never
re-solve), behind ``FEM_EM_WF16_S17=1``; W0 runs the same checks on a synthetic twin.

Anchors, each band imported: (Sel) the selected frequency equals the file's, on an exact
value and on either side of a midpoint; (C=) ``C_mode`` equals the single-frequency tune of
that column, written to a one-frequency file, at ``REDUCTION_BAND``; (ω) at the lowest
frequency C = 1/(2π f Im z_k) on the stored column at ``REDUCTION_BAND``; (S0) 0 solves with
``solve``, ``build`` and ``gmsh.initialize`` patched to raise; (Ref) ``at_hz`` outside the
span, ``tune`` with no ``at_hz``, ``at_hz`` on any other result, a multi-frequency
``.s32p``, and YAML ``load_result:`` on a multi-frequency file.  Negative controls: the
64 MHz column at 128 MHz's ω misses > ``REDUCTION_BAND``; the file with its 64 and 128 MHz
columns exchanged misses > ``LEG_D0_REPRODUCTION_BAND``.

Not claimed: a YAML route, a 32- or 64-port sweep, interpolation between frequencies.
"""

from __future__ import annotations

import os
from pathlib import Path

import gmsh
import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import LumpedPortDrive, Study, run_study_file
from fem_em_solver.ports.circuit import s_to_z
from fem_em_solver.ports.touchstone import load_touchstone
from fem_em_solver.study import FourPortBirdcageGeometry, _four_port_dft

from tests.complex_mode import complex_only
from tests.package.test_wf1_yaml_load_result import _loaded_mapping, _write_yaml
from tests.package.test_wf16_study_load_result import _synthetic_s, _write
from tests.package.test_wf16_study_tune_capacitors import _geometry, _phantom
from tests.package.test_wf16_study_tune_capacitors import _study as _four_port
from tests.package.test_wf16_study_tune_f_human_copper import _study as _f_human
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND

S17_ENV = "FEM_EM_WF16_S17"
SWEEP = Path("/workspace/logs/wf16/wf16_sweep_f_small_p1_n2.s4p")
OUT_DIR = Path("/workspace/logs/wf16_s17")  # gitignored
FREQS = (1.0e7, 6.4e7, 1.28e8)


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _at(f):
    return Study(_geometry(), _phantom(), LumpedPortDrive(), frequency_hz=f, degree=1)


def _zk(s, z0=50.0):
    F = _four_port_dft(4)
    return np.diag(F @ s_to_z(np.asarray(s, np.complex128), z0) @ F.conj().T)


def _boom(*a, **k):
    raise AssertionError("(S0) a multi-frequency tune reached a build or a solve")


def _sweep_checks(path, tmp, tag):
    """(Sel), (C=), (ω), (S0) and both controls on one three-frequency file."""
    comm = MPI.COMM_WORLD
    freqs, stack, _ids = load_touchstone(path)
    freqs = np.asarray(freqs, dtype=np.float64)
    assert tuple(freqs) == FREQS, freqs
    st = _four_port(LumpedPortDrive())
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(Study, "solve", _boom)
        mp.setattr(FourPortBirdcageGeometry, "build", _boom)
        mp.setattr(gmsh, "initialize", _boom)
        res = st.load_result(path)
        assert np.array_equal(res.s_matrix, stack[1]) and res.frequency_hz == FREQS[1]
        assert res.frequencies_hz == FREQS and np.array_equal(res.s_matrices, stack)
        tunes = {}
        for i, f in enumerate(FREQS):
            r = st.tune(k=1, comm=comm, at_hz=f)
            assert r.frequency_hz == f, (f, r.frequency_hz)  # (Sel), exact
            tunes[f] = r
            _p0(f"[WF-16 s17 {tag} Sel] at_hz {f:.3e}: selected {r.frequency_hz:.3e}, C_mode "
                f"{r.c_mode_farad:.9e} F, Im z_1 {r.z_k[1].imag:.6f} ohm")
        for lo, hi in zip(FREQS[:-1], FREQS[1:]):
            below, above = lo + 0.4 * (hi - lo), lo + 0.6 * (hi - lo)
            assert st.tune(k=1, comm=comm, at_hz=below).frequency_hz == lo, below
            assert st.tune(k=1, comm=comm, at_hz=above).frequency_hz == hi, above
    _p0(f"[WF-16 s17 {tag} S0] solve, build and gmsh.initialize patched to raise; every tune "
        "returned")

    worst = 0.0
    for i, f in enumerate(FREQS):  # (C=) one file per rank: the write is COMM_SELF
        one = _write(tmp / f"{tag}_col{i}.s4p", [(f, stack[i])], comm=MPI.COMM_SELF)
        st1 = _at(f)
        st1.load_result(one)
        c1 = st1.tune(k=1, comm=comm).c_mode_farad
        rel = abs(tunes[f].c_mode_farad / c1 - 1.0)
        worst = max(worst, rel)
        _p0(f"[WF-16 s17 {tag} C=] {f:.3e} Hz: sweep {tunes[f].c_mode_farad:.12e} vs one-frequency "
            f"{c1:.12e}, rel {rel:.3e} (REDUCTION_BAND {REDUCTION_BAND:.0e})")
    assert comm.allreduce(worst, op=MPI.MAX) <= REDUCTION_BAND

    f0 = FREQS[0]
    c_w = 1.0 / (2.0 * np.pi * f0 * _zk(stack[0])[1].imag)
    rel_w = abs(tunes[f0].c_mode_farad / c_w - 1.0)
    _p0(f"[WF-16 s17 {tag} omega] {f0:.3e} Hz: C {tunes[f0].c_mode_farad:.12e} vs "
        f"1/(2 pi f Im z_1) {c_w:.12e}, rel {rel_w:.3e} (REDUCTION_BAND {REDUCTION_BAND:.0e})")
    assert rel_w <= REDUCTION_BAND

    # Controls.  The 64 MHz column at 128 MHz's ω: the ceiling, from the file's own z_k.
    c64 = tunes[FREQS[1]].c_mode_farad
    ceil = abs(1.0 / (2.0 * np.pi * FREQS[2] * _zk(stack[1])[1].imag) / c64 - 1.0)
    _p0(f"[WF-16 s17 {tag} ctl] 64 MHz column at 128 MHz's omega: miss {ceil:.6f} (asserted > "
        f"{REDUCTION_BAND:.0e})")
    assert ceil > REDUCTION_BAND
    swapped = _write(tmp / f"{tag}_swapped.s4p",
                     [(FREQS[0], stack[0]), (FREQS[1], stack[2]), (FREQS[2], stack[1])],
                     comm=MPI.COMM_SELF)
    st2 = _four_port(LumpedPortDrive())
    st2.load_result(swapped)
    miss = abs(st2.tune(k=1, comm=comm, at_hz=FREQS[1]).c_mode_farad / c64 - 1.0)
    _p0(f"[WF-16 s17 {tag} ctl] 64 and 128 MHz columns exchanged: miss {miss:.6f} (asserted > "
        f"{LEG_D0_REPRODUCTION_BAND:.0e})")
    assert miss > LEG_D0_REPRODUCTION_BAND


@complex_only
def test_wf16_s17_w0_synthetic_and_refusals(tmp_path):
    comm = MPI.COMM_WORLD
    s = _synthetic_s()
    tag = f"r{comm.size}"
    path = _write(OUT_DIR / f"synthetic_{tag}.s4p",
                  [(FREQS[0], s), (FREQS[1], 0.9 * s), (FREQS[2], 0.8 * s)])
    _sweep_checks(path, tmp_path, "W0")

    st = _four_port(LumpedPortDrive())
    st.load_result(path)
    for bad in (FREQS[0] * 0.5, FREQS[2] * 1.5, float("nan")):
        with pytest.raises(ValueError, match="at_hz"):
            st.tune(k=1, comm=comm, at_hz=bad)
    with pytest.raises(ValueError, match="at_hz"):
        st.tune(k=1, comm=comm)
    _p0("[WF-16 s17 W0 Ref] at_hz outside the span and tune with no at_hz: ValueError naming at_hz")

    # at_hz on any other result, before anything solves.
    fresh = _four_port(LumpedPortDrive())
    with pytest.raises(TypeError, match="at_hz"):
        fresh.tune(k=1, comm=comm, at_hz=FREQS[1])
    assert fresh.result is None
    one = _write(OUT_DIR / f"one_{tag}.s4p", [(FREQS[1], s)])
    fresh.load_result(one)
    assert fresh.result.frequencies_hz == () and fresh.result.s_matrices is None
    with pytest.raises(TypeError, match="at_hz"):
        fresh.tune(k=1, comm=comm, at_hz=FREQS[1])
    ring = _f_human()
    with pytest.raises(TypeError, match="at_hz"):
        ring.tune(k=1, comm=comm, at_hz=FREQS[1])
    _p0("[WF-16 s17 W0 Ref] at_hz on an unloaded, a single-frequency and a ring Study: TypeError")

    # A multi-frequency ring file stays refused, with today's message, nothing stored.
    s32 = 0.1 * np.eye(32, dtype=np.complex128)
    p32 = _write(OUT_DIR / f"ring_{tag}.s32p", [(ring.frequency_hz, s32), (2.0e8, s32)])
    with pytest.raises(ValueError, match="exactly one frequency"):
        ring.load_result(p32)
    assert ring.result is None
    # The Study's own frequency outside the span.
    out = _at(2.0e8)
    with pytest.raises(ValueError, match="span"):
        out.load_result(path)
    assert out.result is None
    _p0("[WF-16 s17 W0 Ref] a multi-frequency .s32p and a Study outside the span: ValueError, "
        "nothing stored")

    # YAML with no `tune.at_hz` (`WF-1` step 27 adds the key): the by-field refusal, a `ValueError` subclass (`OPS-76`).
    yml = _write_yaml(OUT_DIR / f"yaml_{tag}" / "loaded.yaml", _loaded_mapping(str(path)))
    with pytest.raises(ValueError, match="at_hz"):
        run_study_file(yml, output_dir=OUT_DIR / f"yaml_{tag}" / "out", comm=comm)
    _p0("[WF-16 s17 W0 Ref] YAML load_result on a multi-frequency file: ValueError naming at_hz")


@complex_only
@pytest.mark.skipif(os.environ.get(S17_ENV) != "1", reason=f"{S17_ENV}=1 reads U71's sweep file")
def test_wf16_s17_on_the_sweep_file(tmp_path):
    assert SWEEP.is_file(), f"{SWEEP} is absent: this test never re-solves"
    _sweep_checks(SWEEP, tmp_path, "W1")
