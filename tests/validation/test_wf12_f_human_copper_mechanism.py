"""`WF-12` step 3b (§10 chain H12d, §9 item 313) — the mechanism of H12c's O(1)
circuit ↔ in-model identity at copper's ``C_mode``, on stored data only.

Selected by environment only (``FEM_EM_WF12_STEP3B=1``), never by ``-k``.
Numpy on two stored npz files (opened read-only, ``allow_pickle=False``, never
written, never recomputed): no mesh, no solve.  Degree 1, 64 MHz, one copper
fixture, κ off (operator ruling 2026-10-02; ``kappa`` printed only).

Asserted, in order:
(R0) at δ = 0 the copied circuit current route reproduces the stored ``circ``
     to ``REUSE_REPRODUCTION_RTOL`` and ``_residual(meas, circ)`` re-reads
     1.581167190125247 to 1e-9 relative (a miss is a copy defect: stop);
(b)  item 175's two H8a-configuration pairs (19.7327 pF, 113.1081 pF) re-read
     their stored residuals to ``REUSE_REPRODUCTION_RTOL``;
(a)  the termination fit Z_i → Z_i(1 + δ_i), i = 0…31 (30 gap Z_C, two feed
     Z_th), minimum-norm Levenberg–Marquardt from δ = 0 on a finite-difference
     Jacobian (step 1e-7): fit residual ≤ 1e-2 of the δ = 0 miss, and
     max|δ̂_i| ∈ [1e-3, 3e-2] (the 10-07 weekly's bands, unchanged).
     ``PREDICTED RED: (a) range [3e-2, 3e-1]`` (on max|δ̂|).

Printed only: the k-spectrum of δ̂ per ring; the δ̂ cross-check on the
H8b-configuration pair at 113.108 pF (``meas_e2`` / ``circ_e2``); iteration
count and Jacobian singular values; the uniform-δ control (PREDICTED, ruling
313-4) and the seeded random-δ control (PREDICTED).

Brief: ``docs/planning/chunks/WF-12.md`` "§9 item 313" with the 2026-10-08
03:00 review rulings (313-1 … 313-5).  Not claimed: a corrected circuit model,
a re-banded (2), a moved ``C_mode``, any tuning / f0 / Q figure.
"""

from __future__ import annotations

import os
import time

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

STEP3B_ENV = "FEM_EM_WF12_STEP3B"
A2_NPZ = "/workspace/logs/wf12/h12c-cu-s32-a2.npz"
H12BP_NPZ = "/workspace/logs/wf12/h12b-prime-held-s32.npz"

# The weekly's bands (10-07 H12c ruling), copied unchanged.
DELTA_MAX_BAND = (1.0e-3, 3.0e-2)
FIT_RESIDUAL_FRACTION = 1.0e-2
FD_STEP = 1.0e-7
# Item 285's print (`20261006T170222Z_WF-12-step3-a2-W1.log`).
ITEM285_R_ID = 1.581167190125247
R_ID_REREAD_RTOL = 1.0e-9
ITEM285_R_E2 = 3.518653e-01
# Item 175's stored residuals (`20261002T183901Z_WF-12-step2p-W2.log:6405`, `:6410`).
ITEM175_R_TUNED = 1.726249e-02
ITEM175_R_MODE = 2.747383e-01
PREDICTED_RED_RANGE = (3.0e-2, 3.0e-1)
RANDOM_SEED = 20261008


# ============================================================================
# Copied verbatim from `attempt/WF-12-step3-a2-20261006T1702Z` (tip b3cafbfd),
# tests/validation/test_wf12_f_human_copper_c_mode.py, branch lines 130–177
# (themselves lifted from `attempt/WF-10-step1d-20261001T140319Z`, 11b58277).
# The item-175 precedent: copied, never imported from the branch.
# ============================================================================


def _thevenin_feed_impedance(z_c, z0):
    """Z_th of a z0 source in parallel with the gap's Z_C."""
    return z0 * z_c / (z0 + z_c)


def _source_voltage(z_c, z0):
    """The z0 source's open-circuit V_s for V_th = 1 V across Z_C."""
    return (z0 + z_c) / z_c


def _current_route_s(currents_by_drive, feeds_idx, z_th, v_s):
    """``S`` (len(feeds)²) from port currents: column j from the drive of feed j.

    ``currents_by_drive[j]`` is the full current vector with feed j driven at
    V_th = 1 V through Z_th; V_i = δ_ij − Z_th I_i; S_ij = 2V_i/V_s − δ_ij.
    """
    n = len(feeds_idx)
    s = np.zeros((n, n), dtype=np.complex128)
    for j in range(n):
        cur = currents_by_drive[j]
        for i, idx in enumerate(feeds_idx):
            d = 1.0 if i == j else 0.0
            v = d - z_th * cur[idx]
            s[i, j] = 2.0 * v / v_s - d
    return s


def _parallel_currents(z32, z_load, drive):
    """``(Z + diag Z_load)⁻¹ v`` — circuit-side port currents for Thevenin ``v``."""
    return np.linalg.solve(
        np.asarray(z32) + np.diag(z_load), np.asarray(drive, dtype=np.complex128)
    )


def _parallel_loads(n, feeds_idx, z_c, z_feed):
    z_load = np.full(n, z_c, dtype=np.complex128)
    for i in feeds_idx:
        z_load[i] = z_feed
    return z_load


def _y_route_s(s_red, omega, c_f, z0):
    """Reduced 2×2 S (feeds in z0, rest in C) → Z → Y + jωC·I → S at z0."""
    from fem_em_solver.ports.circuit import s_to_z, z_to_s

    z2 = s_to_z(np.asarray(s_red, dtype=np.complex128), z0)
    n = z2.shape[0]
    y = np.linalg.inv(z2) + 1j * omega * c_f * np.eye(n)
    return z_to_s(np.linalg.inv(y), z0)


# ============================================================================
# End of the copied block.
# ============================================================================


def _perturbed_circ(z32, z_c, z_th, v_s, delta, feeds_idx=(0, 1)):
    """CIRC(δ): the current route with every termination Z_i → Z_i(1 + δ_i).

    At δ = 0 this is the copied ``_parallel_loads`` / ``_parallel_currents`` /
    ``_current_route_s`` chain unchanged.  A feed's perturbed Z_th enters both
    the loaded matrix and its own V_i = δ_ij − Z_th,i I_i; V_s stays the
    reference source voltage.
    """
    n = z32.shape[0]
    z_load = _parallel_loads(n, feeds_idx, z_c, z_th) * (1.0 + np.asarray(delta))
    nf = len(feeds_idx)
    s = np.zeros((nf, nf), dtype=np.complex128)
    for j, jdx in enumerate(feeds_idx):
        e = np.zeros(n, dtype=np.complex128)
        e[jdx] = 1.0
        cur = _parallel_currents(z32, z_load, e)
        for i, idx in enumerate(feeds_idx):
            d = 1.0 if i == j else 0.0
            v = d - z_load[idx] * cur[idx]
            s[i, j] = 2.0 * v / v_s - d
    return s


def _lm_fit(fun, n_unknown, r_target, max_iter=200, tol=1e-13):
    """Minimum-norm Levenberg–Marquardt from 0 on a forward-difference Jacobian.

    ``fun(δ)`` is holomorphic in δ, so the complex FD Jacobian is exact to
    O(step).  Step Δ = −Jᴴ(J Jᴴ + λI)⁻¹ r (the minimum-norm GN step at λ → 0).
    Returns ``(δ̂, iterations, sv_first, sv_last)``.
    """
    delta = np.zeros(n_unknown, dtype=np.complex128)
    r = (fun(delta) - r_target).ravel()
    lam = 1.0e-6
    sv_first = sv_last = None
    it = 0
    for it in range(1, max_iter + 1):
        f0 = fun(delta).ravel()
        jac = np.zeros((r.size, n_unknown), dtype=np.complex128)
        for i in range(n_unknown):
            dd = delta.copy()
            dd[i] += FD_STEP
            jac[:, i] = (fun(dd).ravel() - f0) / FD_STEP
        sv = np.linalg.svd(jac, compute_uv=False)
        if sv_first is None:
            sv_first = sv
        sv_last = sv
        accepted = False
        for _ in range(40):
            a = jac @ jac.conj().T + lam * np.eye(r.size)
            step = -jac.conj().T @ np.linalg.solve(a, r)
            trial = delta + step
            r_new = (fun(trial) - r_target).ravel()
            if np.linalg.norm(r_new) < np.linalg.norm(r):
                delta, r = trial, r_new
                lam = max(lam / 10.0, 1e-15)
                accepted = True
                break
            lam *= 10.0
        if not accepted or np.linalg.norm(r) <= tol * np.linalg.norm(r_target):
            break
    return delta, it, sv_first, sv_last


def _k_spectrum(values, slot_k, slot_sign):
    """``{ring sign: folded |k| = 0…8 power fractions}`` of a per-port vector."""
    out = {}
    for sign in (1.0, -1.0):
        vec = np.zeros(16, dtype=np.complex128)
        for idx in range(len(values)):
            if float(slot_sign[idx]) == sign:
                vec[int(slot_k[idx])] = values[idx]
        p = np.abs(np.fft.fft(vec)) ** 2
        tot = float(np.sum(p))
        p = p / tot if tot else p
        out[sign] = [float(p[0])] + [float(p[k] + p[16 - k]) for k in range(1, 8)] + [float(p[8])]
    return out


@complex_only
def test_wf12_step3b_termination_fit():
    if os.environ.get(STEP3B_ENV, "").strip() != "1":
        pytest.skip(f"{STEP3B_ENV} unset: WF-12 step 3b is selected by environment only")

    from fem_em_solver.ports.circuit import s_to_z
    from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
    from tests.validation.test_port_circuit_layer_field import (
        _capacitor_impedance,
        _residual,
        tuned_input,
    )
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ
    from tests.validation.test_wf10_f_human_tuned_circuit import C_MODE_F, C_TUNED_F

    comm = MPI.COMM_WORLD
    if comm.rank != 0:
        return
    t_window = time.perf_counter()

    def say(msg):
        print(f"[WF-12 step3b] {msg}", flush=True)

    for path in (A2_NPZ, H12BP_NPZ):
        if not os.path.exists(path):
            pytest.fail(f"{path} missing: step 3b never recomputes it — stop")
    a2 = np.load(A2_NPZ, allow_pickle=False)
    hb = np.load(H12BP_NPZ, allow_pickle=False)
    say(f"loaded {A2_NPZ} keys {sorted(a2.files)}")
    say(f"loaded {H12BP_NPZ} keys {sorted(hb.files)}")

    z0 = float(REFERENCE_IMPEDANCE_OHM)
    omega = 2.0 * np.pi * float(FREQUENCY_HZ)
    order = [str(p) for p in a2["order"]]
    assert order[:2] == ["P17", "P21"], order[:2]
    slot_k, slot_sign = a2["slot_k"], a2["slot_sign"]
    s32 = np.asarray(a2["s32"], dtype=np.complex128)
    c_mode = float(a2["c_mode"])
    meas = np.asarray(a2["meas"], dtype=np.complex128)
    circ = np.asarray(a2["circ"], dtype=np.complex128)
    say(
        f"order[:2] {order[:2]}; cells {int(a2['cells'])}; C_mode {c_mode * 1e12:.4f} pF; "
        f"kappa {float(a2['kappa']):.6e} (PRINTED ONLY, not applied)"
    )
    z32 = s_to_z(s32, z0)
    z_c = _capacitor_impedance(FREQUENCY_HZ, c_mode)
    z_th = _thevenin_feed_impedance(z_c, z0)
    v_s = _source_voltage(z_c, z0)
    feeds_idx = (0, 1)

    def circ_at(delta, zc=z_c, zth=z_th, vs=v_s):
        return _perturbed_circ(z32, zc, zth, vs, delta, feeds_idx)

    # ---- (R0) reconstruction -------------------------------------------------
    z_load0 = _parallel_loads(32, feeds_idx, z_c, z_th)
    curs0 = [
        _parallel_currents(z32, z_load0, np.eye(32, dtype=np.complex128)[j]) for j in feeds_idx
    ]
    s_ic = _current_route_s(curs0, feeds_idx, z_th, v_s)
    c0 = circ_at(np.zeros(32))
    r_r0 = max(_residual(s_ic, circ), _residual(c0, circ))
    say(f"(R0) verbatim chain vs perturbed route at delta=0 {_residual(c0, s_ic):.3e} (printed)")
    _, s_red = tuned_input(s32, FREQUENCY_HZ, c_mode, terminated=tuple(range(2, 32)))
    r_r0_y = _residual(_y_route_s(s_red, omega, c_mode, z0), circ)
    r_id = _residual(meas, circ)
    r_id_rel = abs(r_id / ITEM285_R_ID - 1.0)
    say(
        f"(R0) current route at delta=0 vs stored circ {r_r0:.3e} (band "
        f"{REUSE_REPRODUCTION_RTOL:g}, ASSERTED); Y route vs stored circ {r_r0_y:.3e} (printed); "
        f"_residual(meas, circ) {r_id:.15f} vs item 285 {ITEM285_R_ID:.15f}, rel {r_id_rel:.3e} "
        f"(band {R_ID_REREAD_RTOL:g}, ASSERTED)"
    )
    assert r_r0 <= REUSE_REPRODUCTION_RTOL, f"(R0) copy defect: {r_r0:.3e}"
    assert r_id_rel <= R_ID_REREAD_RTOL, f"(R0) re-read {r_id!r}"

    # ---- (b) the detuned re-read (H8a configuration, item 175) ---------------
    hb_order = [str(p) for p in hb["order"]]
    s32_h = np.asarray(hb["s32"], dtype=np.complex128)
    b_rows = []
    for label, c_f, mkey, rkey, lit in (
        ("C_tuned", C_TUNED_F, "meas_tuned", "r_tuned", ITEM175_R_TUNED),
        ("C_mode(sigma=800)", C_MODE_F, "meas_mode", "r_mode", ITEM175_R_MODE),
    ):
        _, pred = tuned_input(s32_h, FREQUENCY_HZ, c_f, terminated=tuple(range(2, 32)))
        r = _residual(np.asarray(hb[mkey], dtype=np.complex128), pred)
        stored = float(hb[rkey])
        rel = abs(r / stored - 1.0)
        b_rows.append((label, rel))
        say(
            f"(b) H8a configuration (feeds at 50 Ohm, NOT (a)'s quantity) {label} "
            f"{c_f * 1e12:.4f} pF: recomputed {r:.9e} vs stored {stored:.9e}, rel {rel:.3e} "
            f"(band {REUSE_REPRODUCTION_RTOL:g}, ASSERTED); brief literal {lit:.6e}"
        )
    say(f"(b) order[:2] {hb_order[:2]}; kappa {float(hb['kappa']):.6e} (printed only)")
    for label, rel in b_rows:
        assert rel <= REUSE_REPRODUCTION_RTOL, f"(b) {label} re-read rel {rel:.3e}"

    # ---- (a) the termination fit ---------------------------------------------
    miss0 = float(np.linalg.norm(c0 - meas))
    t_fit = time.perf_counter()
    d_hat, iters, sv0, sv1 = _lm_fit(circ_at, 32, meas)
    t_fit = time.perf_counter() - t_fit
    fit = float(np.linalg.norm(circ_at(d_hat) - meas))
    ratio = fit / miss0
    dmax = float(np.max(np.abs(d_hat)))
    imax = int(np.argmax(np.abs(d_hat)))
    say(
        f"(a) ||CIRC(0) - MEAS||_F {miss0:.6e}; fit ||CIRC(d) - MEAS||_F {fit:.6e}; ratio "
        f"{ratio:.3e} (band <= {FIT_RESIDUAL_FRACTION:g}, ASSERTED); {iters} iterations, "
        f"{t_fit:.2f} s"
    )
    say(
        f"(a) max|d_i| {dmax:.6e} at index {imax} ({order[imax]}, slot k {int(slot_k[imax])}, "
        f"ring {float(slot_sign[imax]):+.0f}) (band {DELTA_MAX_BAND}, ASSERTED; PREDICTED RED "
        f"range {PREDICTED_RED_RANGE}); ||d||_2 {np.linalg.norm(d_hat):.6e}; feeds d "
        f"{d_hat[0]:.4e}, {d_hat[1]:.4e}"
    )
    for i in range(32):
        say(
            f"  d[{i:2d}] {order[i]:>4s} k={int(slot_k[i]):2d} ring {float(slot_sign[i]):+.0f}: "
            f"{d_hat[i].real:+.6e}{d_hat[i].imag:+.6e}j |d| {abs(d_hat[i]):.4e}"
        )
    say(
        f"(a) Jacobian singular values at d=0: {', '.join(f'{v:.4e}' for v in sv0)}; at last "
        f"iterate: {', '.join(f'{v:.4e}' for v in sv1)} (rank "
        f"{int(np.sum(sv1 > 1e-10 * sv1[0]))} of {sv1.size})"
    )
    spec = _k_spectrum(d_hat, slot_k, slot_sign)
    spec_gap = _k_spectrum(np.where(np.arange(32) < 2, 0.0, d_hat), slot_k, slot_sign)
    feed_sign = float(slot_sign[0])
    for sign in (1.0, -1.0):
        tag = "feed ring" if sign == feed_sign else "other ring"
        say(
            f"(a) k-spectrum of d, ring {sign:+.0f} ({tag}), |k|=0..8: "
            f"{', '.join(f'{p:.4f}' for p in spec[sign])}; gap sheets only: "
            f"{', '.join(f'{p:.4f}' for p in spec_gap[sign])} (PREDICTED: k=0/+-1 if the "
            f"copper-wall mechanism, flat if noise)"
        )

    # ---- cross-check: d applied to the H8b-configuration pair at 113.108 pF ---
    c_e2 = float(a2["c_e2"])
    zc2 = _capacitor_impedance(FREQUENCY_HZ, c_e2)
    zth2 = _thevenin_feed_impedance(zc2, z0)
    vs2 = _source_voltage(zc2, z0)
    meas_e2 = np.asarray(a2["meas_e2"], dtype=np.complex128)
    circ_e2 = np.asarray(a2["circ_e2"], dtype=np.complex128)
    ce2_0 = circ_at(np.zeros(32), zc2, zth2, vs2)
    ce2_d = circ_at(d_hat, zc2, zth2, vs2)
    say(
        f"cross-check (PRINTED ONLY) H8b configuration at {c_e2 * 1e12:.4f} pF: route at d=0 vs "
        f"stored circ_e2 {_residual(ce2_0, circ_e2):.3e}; residual at d=0 "
        f"{_residual(meas_e2, ce2_0):.6e} (item 285 {ITEM285_R_E2:.6e}); after d "
        f"{_residual(meas_e2, ce2_d):.6e} (a physical d shrinks it)"
    )

    # ---- negative controls (PREDICTED, printed — ruling 313-4) ----------------
    def circ_uniform(u):
        d = np.zeros(32, dtype=np.complex128)
        d[2:] = u[0]
        return circ_at(d)

    u_hat, u_it, _, _ = _lm_fit(circ_uniform, 1, meas)
    u_fit = float(np.linalg.norm(circ_uniform(u_hat) - meas))
    say(
        f"uniform-d control (PREDICTED > {FIT_RESIDUAL_FRACTION:g} of the miss, printed, not "
        f"asserted): d_u {u_hat[0]:.6e} ({u_it} it), fit ratio {u_fit / miss0:.4e}"
    )
    rng = np.random.default_rng(RANDOM_SEED)
    d_rand = rng.standard_normal(32) + 1j * rng.standard_normal(32)
    d_rand *= np.linalg.norm(d_hat) / np.linalg.norm(d_rand)
    c_rand = circ_at(d_rand)
    say(
        f"random-d control (PREDICTED: change ~ miss, no reduction; seed {RANDOM_SEED}, norm "
        f"{np.linalg.norm(d_rand):.4e}): ||CIRC(d_r) - CIRC(0)|| / miss "
        f"{np.linalg.norm(c_rand - c0) / miss0:.4e}; ||CIRC(d_r) - MEAS|| / miss "
        f"{np.linalg.norm(c_rand - meas) / miss0:.4e}"
    )
    say(f"window {time.perf_counter() - t_window:.2f} s")

    # ---- (a) asserted last ----------------------------------------------------
    assert ratio <= FIT_RESIDUAL_FRACTION, f"(a) fit residual ratio {ratio:.3e}"
    assert DELTA_MAX_BAND[0] <= dmax <= DELTA_MAX_BAND[1], (
        f"(a) max|d| {dmax:.6e} outside {DELTA_MAX_BAND} (PREDICTED RED range "
        f"{PREDICTED_RED_RANGE})"
    )
