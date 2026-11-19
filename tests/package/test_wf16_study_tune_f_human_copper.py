"""`WF-16` step 12 (§9 item 472, chain U53): ``Study.tune(k=1)`` on the copper
F-human hole Study reproduces `WF-12` step 3's H12c mode capacitance.

``Study(FHumanBirdcageGeometry(as_hole=True), gel, LumpedPortDrive(),
conductor=LeontovichWall(COPPER_SIGMA), frequency_hz=64e6, degree=1)`` solves
the 32-port ring-gap matrix (κ off) and ``tune(k)`` reads the ring-antisymmetric
eigenvalue of the k-block of the per-ring DFT of ``Z``; ``C = 1/(omega Im z_k)``.

Anchors (asserted, every band imported, none set here):
  (K0)  no mesh: refusals, and the ring read on a shuffled synthetic 32-port at
        ``EXACT_IDENTITY_RTOL``.
  (R0)  cells == ``HOLE_CELL_RECORD`` == the npz cells; the Study's 32x32 against
        the stored H12c matrix at ``REUSE_REPRODUCTION_RTOL``.
  (I)   the API equals ``ring_mode_rows`` on the same S at ``EXACT_IDENTITY_RTOL``.
  (a)   ``C_mode`` against `WF-12` step 3's printed 1.181185799070e-10 F within
        ``REDUCTION_BAND``.
Negative control (rule (e)), asserted, backed by step 3's own k = 2 control:
  (c1)  k = 2 misses the record by more than ``REDUCTION_BAND``.

Opt-in: ``FEM_EM_WF16_S12=1`` (unset => the solve nodes skip; any other value
raises); width 8 asserted when set.  Scope: degree 1, 64 MHz, κ off, lossless,
one copper hole fixture; not a tuned f0, Q, R_in, C_in or resonance, not ESR,
not degree 2, not AED.
"""

from __future__ import annotations

import os
import time
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.circuit import z_to_s
from fem_em_solver.ports.ring_modes import ring_mode_rows, ring_port_map_from_names
from fem_em_solver.study import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    LeontovichWall,
    LumpedPortDrive,
    SphereInBoxGeometry,
    Study,
    TwoPortQuadratureDrive,
)
from fem_em_solver.utils.instrumentation import report_peak_rss
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_circuit_layer_field import _residual
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE, FREQUENCY_HZ
from tests.validation.test_wf10_f_human_copper_esr_root import C_MODE_CU_PF, ESR_OHM, NPZ_PATH
from tests.validation.test_wf10_f_human_tuned_b1_sar import C_IN_F_HUMAN_PF
from tests.validation.test_wf10_f_human_tuned_circuit import C_MODE_F
from tests.validation.test_wf12_f_human_copper_matrix import COPPER_SIGMA, HOLE_CELL_RECORD

S12_ENV = "FEM_EM_WF16_S12"
RECORD_RANK_WIDTH = 8  # the comparand's width (20261006T170222Z_WF-12-step3-a2-W1.log:219)
# Fixture of every transcribed constant below: degree 1, -n 8, 369 387 cells,
# κ off, Cu sigma = 5.8e7 S/m by Leontovich.
#: H12c k = 1 C_mode (20261006T170222Z_WF-12-step3-a2-W1.log:6366, again :6367).
#: A comparand, never a band; its 13-digit print floor is <= 5e-13 relative.
WF12_STEP3_C_MODE_CU_F = 1.181185799070e-10
#: H12c k = 2 C and |C2/C_mode - 1| (…step3-a2-W1.log:6397): (c1)'s backing, printed.
WF12_STEP3_C2_CU_F = 7.881746394502e-11
WF12_STEP3_K2_SEPARATION = 0.332726
#: off-block norm (:6366), reciprocity and sigma_max (:6359), printed beside.
WF12_STEP3_OFF_BLOCK = 4.195925e-03
WF12_STEP3_RECIPROCITY = 8.137992e-15
WF12_STEP3_SIGMA_MAX = 0.999993339
#: lambda_anti k = 0, 1, 2 (:6361, :6362, :6364), printed beside.
WF12_STEP3_LAMBDA_ANTI = (0.2431 + 14.8162j, 0.3131 + 21.0534j, 0.2555 + 31.5513j)
FEEDS = [(0, "bottom"), (4, "bottom")]
DRIVES = {"(1,+1j)": (1.0, 1.0j)}
OMEGA = 2.0 * np.pi * FREQUENCY_HZ


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S12_ENV, "").strip()
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{S12_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _phantom():
    return GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_HZ)


def _study(geom=None, drive=None):
    return Study(
        FHumanBirdcageGeometry(as_hole=True) if geom is None else geom,
        _phantom(),
        LumpedPortDrive() if drive is None else drive,
        frequency_hz=FREQUENCY_HZ,
        degree=DEGREE,
        conductor=LeontovichWall(study_mod.COPPER_SIGMA),
    )


def _rows(s, names, c_f):
    return ring_mode_rows(
        s,
        port_map=ring_port_map_from_names(names),
        freq_hz=FREQUENCY_HZ,
        c_f=c_f,
        feeds=FEEDS,
        drives=DRIVES,
    )


def _synthetic():
    """A 32-port with a known block-diagonal ring spectrum (no mesh)."""
    rng = np.random.default_rng(472)
    a = np.zeros(16, dtype=np.complex128)
    b = np.zeros(16, dtype=np.complex128)
    for m in range(16):
        b[m] = complex(rng.uniform(0.1, 2.0), rng.uniform(-8.0, 8.0))
        d = complex(rng.uniform(0.1, 1.0), rng.uniform(5.0, 60.0))
        if m == 1:
            d = 0.3 + 21.05j  # ≈ the copper value (…step3-a2-W1.log:6362)
        a[m] = b[m] + d
    zh = np.zeros((32, 32), dtype=np.complex128)
    for m in range(16):
        zh[m, m] = zh[16 + m, 16 + m] = a[m]
        zh[m, 16 + m] = zh[16 + m, m] = b[m]
    j = np.arange(16)
    f = np.exp(-2j * np.pi * np.outer(j, j) / 16) / 4.0
    w = np.kron(np.eye(2), f)
    z = w.conj().T @ zh @ w
    s = z_to_s(z, 50.0)
    names = [f"P{17 + n}" for n in range(32)]
    return s, names, a - b


def test_wf16_s12_k0_refusals_and_synthetic(monkeypatch):
    """(K0): refusals with no mesh, and the ring read on a shuffled synthetic."""
    _knob()  # a bad knob value raises here too
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called by a refused tune")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry):
        monkeypatch.setattr(cls, "build", stub)
    quad = _study(drive=TwoPortQuadratureDrive(118.0e-12, esr_ohm=ESR_OHM))
    with pytest.raises(TypeError):
        quad.tune()
    solid = Study(
        FHumanBirdcageGeometry(),
        _phantom(),
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_HZ,
        degree=DEGREE,
        conductor=LeontovichWall(),
    )
    with pytest.raises(TypeError):
        solid.solve()
    st = _study()
    for bad in (0, 3):
        with pytest.raises(ValueError):
            st.tune(k=bad)
    for bad in (1.0, True):
        with pytest.raises(TypeError):
            st.tune(k=bad)
    assert calls["n"] == 0 and st.result is None

    s, names, d = _synthetic()
    rng = np.random.default_rng(12)
    p = rng.permutation(32)
    s_sh = s[np.ix_(p, p)]
    names_sh = [names[i] for i in p]
    st.result = SimpleNamespace(
        s_matrix=s_sh, port_definitions=[SimpleNamespace(port_id=n) for n in names_sh]
    )
    t1, t2 = st.tune(k=1), st.tune(k=2)
    ref = {1: 1.0 / (OMEGA * d[1].imag), 2: 1.0 / (OMEGA * d[2].imag)}
    rows = _rows(s_sh, names_sh, t1.c_mode_farad)
    e1 = abs(t1.c_mode_farad / ref[1] - 1.0)
    e2 = abs(t2.c_mode_farad / ref[2] - 1.0)
    r1 = abs(t1.c_mode_farad / rows["c_mode_f"] - 1.0)
    r2 = abs(t2.c_mode_farad / rows["c_k2_f"] - 1.0)
    _p0(
        f"\n[WF-16 s12 K0] synthetic C1 {t1.c_mode_farad:.12e} F (rel vs closed {e1:.3e}, vs "
        f"ring_mode_rows {r1:.3e}); C2 {t2.c_mode_farad:.12e} F (rel {e2:.3e}, {r2:.3e}); "
        f"off {t1.off_diagonal_residual:.3e}; band EXACT_IDENTITY_RTOL {EXACT_IDENTITY_RTOL:.1e}"
    )
    assert e1 <= EXACT_IDENTITY_RTOL and e2 <= EXACT_IDENTITY_RTOL
    assert r1 <= EXACT_IDENTITY_RTOL and r2 <= EXACT_IDENTITY_RTOL
    assert t1.off_diagonal_residual <= EXACT_IDENTITY_RTOL
    assert t1.k == 1 and t2.k == 2 and t1.frequency_hz == FREQUENCY_HZ
    plain = _study()
    plain.result = SimpleNamespace(
        s_matrix=s, port_definitions=[SimpleNamespace(port_id=n) for n in names]
    )
    # The brief asked "bit for bit"; measured 20261015T213212Z_WF-16-s12-W0.log:63:
    # 1.1813757652308145e-10 vs …148e-10 (1 ulp, ≈ 2.5e-16 relative) — s_to_z
    # inverts the shuffled and unshuffled S in different row orders, so LAPACK
    # rounds differently.  The bound was wrong, not the map (both reads equal the
    # closed form at 0.000e+00, same log :52); asserted at EXACT_IDENTITY_RTOL.
    for k, t in ((1, t1), (2, t2)):
        assert abs(plain.tune(k=k).c_mode_farad / t.c_mode_farad - 1.0) <= EXACT_IDENTITY_RTOL
    assert calls["n"] == 0


@pytest.fixture(scope="module")
def tuned():
    if not _knob():
        pytest.skip(
            f"{S12_ENV} unset: the copper F-human tune is a heavy -n 8 window, "
            "environment-selected"
        )
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, (comm.size, RECORD_RANK_WIDTH)
    npz = None
    if comm.rank == 0:
        assert os.path.exists(NPZ_PATH), f"{NPZ_PATH} absent: stop (BLOCKED, rule (d))"
        with np.load(NPZ_PATH, allow_pickle=False) as fh:
            npz = {
                "s32": np.array(fh["s32"], dtype=np.complex128),
                "order": [str(x) for x in fh["order"]],
                "cells": int(fh["cells"]),
                "c_mode": float(fh["c_mode"]),
            }
    npz = comm.bcast(npz, root=0)
    study = _study()
    assert study.degree == 1
    assert isinstance(study.conductor, LeontovichWall)
    assert study_mod.COPPER_SIGMA == COPPER_SIGMA
    t0 = time.perf_counter()
    res = study.solve(comm)
    t_solve = time.perf_counter() - t0
    t1 = study.tune(k=1, comm=comm)
    t2 = study.tune(k=2, comm=comm)
    assert study.result is res
    _p0(
        f"\n[WF-16 s12] cells {res.cells}, mesh {res.mesh_time:.1f} s, sweep "
        f"{res.sweep_time:.1f} s, total {t_solve:.1f} s, -n {comm.size}"
    )
    report_peak_rss(comm, label="[WF-16 s12] summed peak ru_maxrss")
    s = np.asarray(res.s_matrix, dtype=np.complex128)
    names = [pd.port_id for pd in res.port_definitions]
    return {"study": study, "res": res, "s": s, "names": names, "t1": t1, "t2": t2, "npz": npz}


@complex_only
def test_wf16_s12_r0_fixture_tie(tuned):
    """(R0): the Study is the H12c fixture and its 32x32 ties the stored matrix."""
    res, npz, s, names = tuned["res"], tuned["npz"], tuned["s"], tuned["names"]
    idx = {n: i for i, n in enumerate(names)}
    p = np.array([idx[n] for n in npz["order"]], dtype=int)
    s_p = s[np.ix_(p, p)]
    rel = _residual(s_p, npz["s32"])
    recip = float(np.linalg.norm(s - s.T) / np.linalg.norm(s))
    smax = float(np.linalg.svd(s, compute_uv=False)[0])
    col = np.linalg.norm(s_p - npz["s32"], axis=0) / np.linalg.norm(npz["s32"], axis=0)
    _p0(
        f"\n[WF-16 s12 R0] cells {res.cells} vs record {HOLE_CELL_RECORD} (==, band "
        f"CELL_COUNT_BAND {CELL_COUNT_BAND} printed) vs npz {npz['cells']} (==); "
        f"||S_study - S_npz||/||S_npz|| = {rel:.3e} (band REUSE_REPRODUCTION_RTOL "
        f"{REUSE_REPRODUCTION_RTOL:.1e}; predicted <= 1e-10); worst column {col.max():.3e}; "
        f"reciprocity {recip:.6e} (step 3 {WF12_STEP3_RECIPROCITY:.6e}), sigma_max "
        f"{smax:.9f} (step 3 {WF12_STEP3_SIGMA_MAX:.9f})"
    )
    assert res.cells == HOLE_CELL_RECORD
    assert res.cells == npz["cells"]
    assert sorted(names) == sorted(npz["order"])
    assert rel <= REUSE_REPRODUCTION_RTOL


@complex_only
def test_wf16_s12_i_api_equals_ring_mode_rows(tuned):
    """(I) and (Im): the API equals ``ring_mode_rows`` on the same S."""
    t1, t2 = tuned["t1"], tuned["t2"]
    rows = _rows(tuned["s"], tuned["names"], t1.c_mode_farad)
    r1 = abs(t1.c_mode_farad / rows["c_mode_f"] - 1.0)
    r2 = abs(t2.c_mode_farad / rows["c_k2_f"] - 1.0)
    ro = abs(t1.off_diagonal_residual / rows["off_block"] - 1.0)
    z1 = t1.z_k[1]
    rz = abs(t1.c_mode_farad * OMEGA * z1.imag - 1.0)
    _p0(
        f"\n[WF-16 s12 I] rel vs ring_mode_rows: C1 {r1:.3e}, C2 {r2:.3e}, off {ro:.3e}; "
        f"C1 w Im z1 - 1 = {rz:.3e} (band EXACT_IDENTITY_RTOL {EXACT_IDENTITY_RTOL:.1e})"
    )
    for k in range(3):
        z = t1.z_k[k]
        _p0(
            f"    z_{k} = {z.real:+.4f}{z.imag:+.4f}j ohm (step 3 "
            f"{WF12_STEP3_LAMBDA_ANTI[k].real:+.4f}{WF12_STEP3_LAMBDA_ANTI[k].imag:+.4f}j)"
        )
    _p0(f"    off-block {t1.off_diagonal_residual:.6e} (step 3 {WF12_STEP3_OFF_BLOCK:.6e})")
    assert r1 <= EXACT_IDENTITY_RTOL
    assert r2 <= EXACT_IDENTITY_RTOL
    assert ro <= EXACT_IDENTITY_RTOL
    assert rz <= EXACT_IDENTITY_RTOL
    assert z1.imag > 0.0
    assert t1.frequency_hz == FREQUENCY_HZ and t1.k == 1


@complex_only
def test_wf16_s12_a_c_mode_equals_wf12_h12c(tuned):
    """(a): C_mode within REDUCTION_BAND of `WF-12` step 3's H12c reading."""
    c = tuned["t1"].c_mode_farad
    gap = c / WF12_STEP3_C_MODE_CU_F - 1.0
    _p0(
        f"\n[WF-16 s12 a] C_mode = {c:.12e} F; record {WF12_STEP3_C_MODE_CU_F:.12e} F; gap "
        f"{gap:+.6e} (band REDUCTION_BAND {REDUCTION_BAND:.1e}; predicted <= 1e-9); "
        f"|gap| vs LEG_D0_REPRODUCTION_BAND {LEG_D0_REPRODUCTION_BAND:.1e} printed; vs "
        f"C_MODE_CU_PF {c / (C_MODE_CU_PF * 1e-12) - 1:+.3e}; vs npz c_mode "
        f"{c / tuned['npz']['c_mode'] - 1:+.3e}; vs sigma = 800 C_MODE_F "
        f"{c / C_MODE_F - 1:+.4e} (different mesh, printed); vs C_in F-human "
        f"{c / (C_IN_F_HUMAN_PF * 1e-12) - 1:+.4e} (different quantity, printed)"
    )
    assert abs(gap) <= REDUCTION_BAND


@complex_only
def test_wf16_s12_c1_mode_index_control(tuned):
    """(c1): k = 2 on the same Study misses the record (predicted -33.27 %)."""
    c2 = tuned["t2"].c_mode_farad
    gap = c2 / WF12_STEP3_C_MODE_CU_F - 1.0
    _p0(
        f"\n[WF-16 s12 c1] k = 2 C = {c2:.12e} F, gap {gap:+.6e} "
        f"({abs(gap) / REDUCTION_BAND:.1f}x band; step 3 separation "
        f"{WF12_STEP3_K2_SEPARATION}); vs step 3 C2 {c2 / WF12_STEP3_C2_CU_F - 1:+.3e}"
    )
    assert abs(gap) > REDUCTION_BAND
