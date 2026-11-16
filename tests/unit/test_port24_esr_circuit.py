"""`PORT-24` step 1 (§9 item 308, chain T25) — capacitor ESR in the circuit layer, no solve.

Operator binding ruling 2026-10-06: every lumped capacitor carries a 0.2 Ω
series resistance.  Step 1 gives ``ring_mode_rows`` / ``ring_mode_rows_v2`` /
``terminated_currents`` an ``esr_ohm`` keyword (default 0) and, per the
(308-scope) review ruling, re-implements the F-small four-port terminations
here (``tuned_input`` is *not* edited) with their r = 0 path asserted bitwise
equal to the existing helper by import.

Anchors (asserted):

* (B) default bitwise — with ``esr_ohm`` unset, every output of the edited
  functions on the stored F-human 32×32 (``logs/ans8/s32.npz``) equals the
  **pre-change** source's output bitwise; the pre-change source is loaded from
  ``git show PRE_CHANGE_SHA:…`` (pinned, never ``HEAD:``).
* (S) at r = 0.2 Ω: ``Re λ_k(r) − Re λ_k(0) = r`` at 1e-12 relative (to r) for
  every k, on F-human (k = 0…8, through ``ring_mode_rows``) and on the F-small
  4×4 at 10 / 64 / 128 MHz (λ_k the N = 4 DFT diagonal, local helper, ruling
  Q2) — and ``Im λ_k``, ``C_mode`` unchanged within 1e-12 relative
  (``SHIFT_RTOL``). The Done-when was re-stated from 'bitwise' by the
  2026-10-10 weekly (T25): ``eig(A + rI)`` has no bitwise contract with
  ``eig(A)``. Item 308 measured 2.889e-15 / 4.571e-16.
* (P) planted non-uniform r (one gap 0.4 Ω, the rest 0.2): the uniform shift
  fails on every mode (predicted ≈ 3 % off on F-human, asserted as the control).

(T) *predicted, printed*: termination-only wiring leaves λ unchanged and misses
(S) by exactly r.
"""

from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import types
from pathlib import Path

import numpy as np
import pytest

from fem_em_solver.ports import ring_modes as rm
from fem_em_solver.ports.circuit import reduce_terminated_ports, s_to_z
from fem_em_solver.ports.lumped import series_rlc_impedance

REPO = Path(__file__).resolve().parents[2]
#: The commit before `PORT-24` step 1 touched ``ring_modes.py`` (pinned).
PRE_CHANGE_SHA = "b2174873"
RING_MODES_PATH = "src/fem_em_solver/ports/ring_modes.py"
ANS8_CASE = REPO / "examples" / "ansys_benchmarks" / "ans8_f_human_quadrature_mode_content_64MHz"
ANS8_SCRIPT = ANS8_CASE / "08_f_human_quadrature_mode_content.py"
S32_NPZ = REPO / "logs" / "ans8" / "s32.npz"
ANS4_METRICS = REPO / "examples" / "ansys_benchmarks" / "ans4_birdcage_four_port_10_64_128MHz" / "metrics.json"
ANS6_METRICS = REPO / "examples" / "ansys_benchmarks" / "ans6_copper_birdcage_four_port_10_64_128MHz" / "metrics.json"
F_SMALL_RUNGS = ("10 MHz", "64 MHz", "128 MHz")
#: `ans:7`'s tuned C (``ans7 metrics.json`` ``c_f``) — the F-small termination C.
F_SMALL_C_F = 6.629e-11
Z0 = 50.0

ESR_OHM = 0.2          # operator ruling 2026-10-06
PLANT_OHM = 0.4        # (P): one gap
SHIFT_RTOL = 1.0e-12   # (S): the Done-when's 1e-12 relative (to r)
HELPER_RTOL = 1.0e-14  # ruling Q2: local DFT helper vs `k_blocks`


# -- fixtures -------------------------------------------------------------------
@pytest.fixture(scope="module")
def pre():
    """The pre-change ``ring_modes`` module, from the pinned commit."""
    src = subprocess.run(
        ["git", "-c", f"safe.directory={REPO}", "-C", str(REPO), "show",
         f"{PRE_CHANGE_SHA}:{RING_MODES_PATH}"],
        check=True, capture_output=True, text=True).stdout
    assert "esr_ohm" not in src, "pinned pre-change source already carries esr_ohm"
    mod = types.ModuleType("fem_em_solver.ports._ring_modes_pre_port24")
    mod.__package__ = "fem_em_solver.ports"
    exec(compile(src, f"{PRE_CHANGE_SHA}:{RING_MODES_PATH}", "exec"), mod.__dict__)
    return mod


@pytest.fixture(scope="module")
def ans8():
    spec = importlib.util.spec_from_file_location("ans8_port24_ut", ANS8_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ans8_port24_ut"] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def f_human(ans8):
    if not S32_NPZ.exists():
        pytest.fail(f"stored F-human record {S32_NPZ} absent")
    saved = np.load(S32_NPZ, allow_pickle=False)
    s32 = np.asarray(saved["s32"], dtype=np.complex128)
    pm = ans8._slot_port_map(saved)
    names = [str(n) for n in saved["order"]]
    freq = float(saved["freq_hz"])
    feeds = [pm[names.index(n)] for n in ans8.SPEC_FEED_NAMES]
    return {"s": s32, "pm": pm, "names": names, "freq": freq, "feeds": feeds}


def _s_from_json(rows):
    def c(x):
        return complex(x["re"], x["im"]) if isinstance(x, dict) else complex(x[0], x[1])
    return np.array([[c(x) for x in row] for row in rows], dtype=np.complex128)


def _f_small_records():
    out = []
    a4 = json.loads(ANS4_METRICS.read_text())["rungs"]
    a6 = json.loads(ANS6_METRICS.read_text())["rungs"]
    for rung in F_SMALL_RUNGS:
        out.append((f"ans4 {rung}", float(a4[rung]["frequency_hz"]), _s_from_json(a4[rung]["s_matrix"])))
        f6 = float(rung.split()[0]) * 1.0e6
        out.append((f"ans6-cu {rung}", f6, _s_from_json(a6[rung]["cu"]["s_matrix"])))
    return out


def _canon(x):
    return json.dumps(x, sort_keys=True, allow_nan=True)


# -- the F-small local helper (ruling Q2) and the re-implemented terminations ----
def dft_k_domain(z, n_slots, n_rings):
    """``W Z Wᴴ`` with ``W = I_rings ⊗ F``, ``F_kj = e^{−2πi jk/N}/√N`` — the
    normalisation ``ring_mode_rows`` uses at N = 16 (`_dft_matrix`)."""
    j = np.arange(n_slots)
    f = np.exp(-2j * np.pi * np.outer(j, j) / n_slots) / np.sqrt(n_slots)
    w = np.kron(np.eye(n_rings), f)
    return w @ np.asarray(z, dtype=np.complex128) @ w.conj().T


def f_small_lambdas(s4, freq_hz, *, esr_ohm=0.0, z0=Z0):
    """``(λ_k for k = 0…3, C_mode)`` of the F-small 4×4: the N = 4 DFT diagonal of
    ``Z + diag(r)`` (ruling Q2), ``C_mode = 1/(ω Im λ₁)``."""
    z = s_to_z(np.asarray(s4, dtype=np.complex128), z0)
    if not (np.ndim(esr_ohm) == 0 and float(esr_ohm) == 0.0):
        z = z + np.diag(np.broadcast_to(np.asarray(esr_ohm, dtype=float), (z.shape[0],)))
    lam = np.diag(dft_k_domain(z, z.shape[0], 1)).copy()
    return lam, 1.0 / (2.0 * np.pi * float(freq_hz) * lam[1].imag)


def capacitor_impedance_esr(frequency_hz, c_f, esr_ohm=0.0):
    """``r + 1/(jωC)``; r = 0 keeps the old expression (a branch, not arithmetic)."""
    if float(esr_ohm) == 0.0:
        return series_rlc_impedance(frequency_hz, c_f=float(c_f))
    return series_rlc_impedance(frequency_hz, r_ohm=float(esr_ohm), c_f=float(c_f))


def tuned_input_esr(s4, frequency_hz, c_f, terminated=(1, 2, 3), esr_ohm=0.0):
    """``tuned_input`` with ``esr_ohm`` (test-local, per the 308-scope ruling)."""
    z_c = capacitor_impedance_esr(frequency_hz, c_f, esr_ohm)
    s_red = reduce_terminated_ports(np.asarray(s4, dtype=np.complex128), Z0,
                                    {k: z_c for k in terminated})
    s11 = complex(s_red[0, 0])
    if abs(1.0 - s11) <= 1.0e-14:
        return complex(np.inf, np.inf), s_red
    return Z0 * (1.0 + s11) / (1.0 - s11), s_red


# -- helper identity (ruling Q2) -------------------------------------------------
def test_local_dft_helper_matches_k_blocks(f_human):
    perm = rm._canonical_permutation(f_human["pm"])
    zc = s_to_z(f_human["s"], Z0)[np.ix_(perm, perm)]
    blocks, _ = rm.k_blocks(zc)
    zh = dft_k_domain(zc, rm.N_SLOTS, 2)
    worst = 0.0
    for k in range(rm.N_SLOTS):
        idx = [k, rm.N_SLOTS + k]
        worst = max(worst, float(np.linalg.norm(zh[np.ix_(idx, idx)] - blocks[k])
                                 / np.linalg.norm(blocks[k])))
    print(f"[PORT-24 Q2] local helper vs k_blocks on F-human 32x32: max rel {worst:.3e} "
          f"(band {HELPER_RTOL:.0e})")
    assert worst <= HELPER_RTOL


# -- (B) default bitwise vs the pinned pre-change source ---------------------------
def test_default_bitwise_f_human(pre, ans8, f_human):
    s, pm, freq, feeds = f_human["s"], f_human["pm"], f_human["freq"], f_human["feeds"]
    controls = ans8.spec_controls(pm, f_human["names"])
    kw = dict(port_map=pm, freq_hz=freq, feeds=feeds, drives=ans8.SPEC_DRIVES)
    for c in (ans8.SPEC_C_F, ans8.C2_RECORD_F):
        old = pre.ring_mode_rows(s, c_f=c, **kw)
        assert _canon(rm.ring_mode_rows(s, c_f=c, **kw)) == _canon(old)
        assert _canon(rm.ring_mode_rows(s, c_f=c, esr_ohm=0.0, **kw)) == _canon(old)
        old2 = pre.ring_mode_rows_v2(s, c_spec_f=c, controls=controls, **kw)
        assert _canon(rm.ring_mode_rows_v2(s, c_spec_f=c, controls=controls, **kw)) == _canon(old2)
    perm = rm._canonical_permutation(pm)
    zc = s_to_z(s, Z0)[np.ix_(perm, perm)]
    tk = dict(freq_hz=freq, c_f=ans8.SPEC_C_F, feed_idx=[0, 4], v_feeds=[1.0, 1.0j])
    v0, i0 = pre.terminated_currents(zc, **tk)
    v1, i1 = rm.terminated_currents(zc, **tk)
    assert np.array_equal(v0, v1) and np.array_equal(i0, i1)
    # The stored record's own path (stage R's rows_v2, as `ans:8` writes it).
    assert _canon(ans8.rows_v2(s, pm, f_human["names"], freq)) == _canon(
        pre.ring_mode_rows_v2(s, c_spec_f=ans8.SPEC_C_F, controls=controls, **kw))
    print("[PORT-24 B] F-human 32x32: ring_mode_rows / _v2 / terminated_currents bitwise "
          f"vs {PRE_CHANGE_SHA} at C = SPEC_C_F and C2_RECORD_F")


def test_default_bitwise_f_small_terminations():
    tif = pytest.importorskip("tests.validation.test_port_circuit_layer_field")
    n = 0
    for label, freq, s4 in _f_small_records():
        assert capacitor_impedance_esr(freq, F_SMALL_C_F) == tif._capacitor_impedance(freq, F_SMALL_C_F)
        for term in (tif.TUNING_TERMINATED_INDICES, tif.TUNING_2X2_TERMINATED_INDICES):
            z_old, s_old = tif.tuned_input(s4, freq, F_SMALL_C_F, terminated=term)
            z_new, s_new = tuned_input_esr(s4, freq, F_SMALL_C_F, terminated=term)
            assert z_new == z_old and np.array_equal(s_new, s_old), label
            n += 1
        # r = 0.2: the termination is exactly r + 1/(jωC).
        dz = capacitor_impedance_esr(freq, F_SMALL_C_F, ESR_OHM) - capacitor_impedance_esr(freq, F_SMALL_C_F)
        assert abs(dz - ESR_OHM) <= SHIFT_RTOL * ESR_OHM, (label, dz)
        z_r, _ = tuned_input_esr(s4, freq, F_SMALL_C_F, esr_ohm=ESR_OHM)
        z_0, _ = tuned_input_esr(s4, freq, F_SMALL_C_F)
        print(f"[PORT-24 B] {label}: Z_in(C) r=0 {z_0:.6g}  r=0.2 {z_r:.6g}")
    print(f"[PORT-24 B] F-small tuned_input r=0 path bitwise on {n} (record, termination) pairs")


# -- (S) uniform shift ---------------------------------------------------------------
def _f_human_modes(f_human, ans8, esr):
    return rm.ring_mode_rows(f_human["s"], port_map=f_human["pm"], freq_hz=f_human["freq"],
                             c_f=ans8.SPEC_C_F, feeds=f_human["feeds"],
                             drives=ans8.SPEC_DRIVES, esr_ohm=esr)


def _shift_table(f_human, ans8):
    """``[(record, k, ΔRe−r rel, Im rel, λ0, λr)]`` and ``[(record, C_mode rel)]``."""
    rows, cmode = [], []
    r0 = _f_human_modes(f_human, ans8, 0.0)
    rr = _f_human_modes(f_human, ans8, ESR_OHM)
    for k in range(rm.N_SLOTS // 2 + 1):
        l0 = complex(*r0["modes"][str(k)]["lambda_anti"])
        lr = complex(*rr["modes"][str(k)]["lambda_anti"])
        rows.append(("F-human", k, abs((lr.real - l0.real) - ESR_OHM) / ESR_OHM,
                     abs(lr.imag - l0.imag) / abs(l0.imag), l0, lr))
    cmode.append(("F-human", abs(rr["c_mode_f"] - r0["c_mode_f"]) / abs(r0["c_mode_f"])))
    for label, freq, s4 in _f_small_records():
        l0, c0 = f_small_lambdas(s4, freq)
        lr, cr = f_small_lambdas(s4, freq, esr_ohm=ESR_OHM)
        for k in range(4):
            rows.append((label, k, abs((lr[k].real - l0[k].real) - ESR_OHM) / ESR_OHM,
                         abs(lr[k].imag - l0[k].imag) / abs(l0[k].imag), l0[k], lr[k]))
        cmode.append((label, abs(cr - c0) / abs(c0)))
    return rows, cmode


def test_uniform_shift_re(f_human, ans8):
    rows, _ = _shift_table(f_human, ans8)
    for rec, k, re_rel, im_rel, l0, lr in rows:
        print(f"[PORT-24 S] {rec:16s} k={k}: lambda0 {l0:.10g}  |dRe-r|/r {re_rel:.3e}")
    worst = max(r[2] for r in rows)
    print(f"[PORT-24 S] max |dRe - r|/r = {worst:.3e} over {len(rows)} (record, k) (band {SHIFT_RTOL:.0e})")
    assert worst <= SHIFT_RTOL
    # (T) predicted, printed: termination-only wiring leaves lambda unchanged,
    # so it misses (S) by |0 - r|/r = 1 exactly on every mode.
    print("[PORT-24 T] termination-only wiring: |dRe - r|/r = 1.000e+00 on every mode (lambda from Z alone)")
    # The v2 composition runs with ESR on the stored record.
    v2 = rm.ring_mode_rows_v2(f_human["s"], port_map=f_human["pm"], freq_hz=f_human["freq"],
                              c_spec_f=ans8.SPEC_C_F, feeds=f_human["feeds"],
                              drives=ans8.SPEC_DRIVES,
                              controls=ans8.spec_controls(f_human["pm"], f_human["names"]),
                              esr_ohm=ESR_OHM)
    print(f"[PORT-24 S] v2 at r=0.2: c_own {v2['c_own_f']:.10g} F  f0 {v2['f0_one_point_hz']:.10g} Hz  "
          f"share_k1(ccw,bottom) {v2['fixed_c']['drives']['ccw']['rings']['bottom']['share_k1']:.6f}")


def test_uniform_shift_im_and_c_mode_unchanged(f_human, ans8):
    rows, cmode = _shift_table(f_human, ans8)
    im_worst = max(r[3] for r in rows)
    c_worst = max(c[1] for c in cmode)
    for rec, val in cmode:
        print(f"[PORT-24 S-Im] {rec:16s}: |dC_mode|/C_mode {val:.3e}")
    print(f"[PORT-24 S-Im] max |dIm|/|Im| = {im_worst:.3e}; max |dC_mode|/C_mode = {c_worst:.3e} "
          "(ASSERTED <= 1e-12 relative, SHIFT_RTOL; item 308 read 2.889e-15 / 4.571e-16)")
    assert im_worst <= SHIFT_RTOL and c_worst <= SHIFT_RTOL, (im_worst, c_worst)


# -- (P) planted non-uniform r ----------------------------------------------------------
def test_planted_nonuniform_breaks_shift(f_human, ans8):
    n = len(f_human["pm"])
    r_plant = np.full(n, ESR_OHM)
    r_plant[f_human["names"].index("P17")] = PLANT_OHM
    r0 = _f_human_modes(f_human, ans8, 0.0)
    rp = _f_human_modes(f_human, ans8, r_plant)
    offs = []
    for k in range(rm.N_SLOTS // 2 + 1):
        l0 = complex(*r0["modes"][str(k)]["lambda_anti"])
        lp = complex(*rp["modes"][str(k)]["lambda_anti"])
        off = abs((lp.real - l0.real) - ESR_OHM) / ESR_OHM
        offs.append(off)
        print(f"[PORT-24 P] F-human k={k}: dRe {lp.real - l0.real:.6f} ohm  |dRe-r|/r {off:.3e}")
    # (P-Im) predicted, printed only: can the re-stated Im clause see a plant?
    im_plant = max(
        abs(complex(*rp["modes"][str(k)]["lambda_anti"]).imag
            - complex(*r0["modes"][str(k)]["lambda_anti"]).imag)
        / abs(complex(*r0["modes"][str(k)]["lambda_anti"]).imag)
        for k in range(rm.N_SLOTS // 2 + 1))
    print(f"[PORT-24 P-Im] F-human plant: max_k |dIm|/|Im| {im_plant:.3e}  "
          f"|C_mode(plant)/C_mode(0) - 1| {abs(rp['c_mode_f'] / r0['c_mode_f'] - 1.0):.3e} (printed)")
    for label, freq, s4 in _f_small_records():
        rpl = np.full(4, ESR_OHM)
        rpl[0] = PLANT_OHM
        l0, _ = f_small_lambdas(s4, freq)
        lp, _ = f_small_lambdas(s4, freq, esr_ohm=rpl)
        for k in range(4):
            off = abs((lp[k].real - l0[k].real) - ESR_OHM) / ESR_OHM
            offs.append(off)
            print(f"[PORT-24 P] {label:16s} k={k}: |dRe-r|/r {off:.3e}")
    print(f"[PORT-24 P] min |dRe - r|/r over every mode = {min(offs):.3e} (must exceed {SHIFT_RTOL:.0e})")
    assert min(offs) > SHIFT_RTOL
