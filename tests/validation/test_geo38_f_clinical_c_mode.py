"""GEO-38 C10a (§9 item 581, §10 chain C10): an N-slot two-ring ``C_mode``
reader, proven on `WF-12`'s degree-1 and H16b's degree-2 copper 32×32 records;
C10's gated test lands skipped.

``ring_mode_rows`` (``src/fem_em_solver/ports/ring_modes.py``) is 16-slot only
(``N_SLOTS = 16``; the port map must hold exactly 32 pairs) and raises on C3's
64 ports.  :func:`ring_c_modes` is the same algebra generalised in N inside
this test — ``Z = s_to_z(S)``, the canonical (bottom 0…N−1, top 0…N−1) order,
``Ẑ = W Z Wᴴ`` with ``W = I₂ ⊗ F_N`` (a copy of ``dft_k_domain``,
``tests/unit/test_port24_esr_circuit.py:125–131``, copied per the `PORT-24`
ruling Q2 precedent), ``antisymmetric_eigenpair`` per k = 0…N/2 and
``C_k = 1/(ω Im λ_anti)`` — and it is proven identical to ``ring_mode_rows``
at N = 16 ((K2), (I16)).  Lossless: no ``esr_ohm`` (a uniform ESR leaves
``C_mode`` unchanged, §9 ESR ruling (i)).  No ``src/`` edit.

Knobs (``FEM_EM_GEO38_C10A``, ``FEM_EM_GEO38_C10``): unset / empty ⇒ the gated
test skips; ``1`` ⇒ it runs and an absent input **fails**; anything else ⇒
``ValueError``.

* ``test_c10a_coding`` (ungated): (K1) closed form at N = 32 and N = 16, ports
  shuffled; (K2) the local reader equals ``ring_mode_rows`` at N = 16; (K3) the
  geometry map on C3a's ideal layout equals the name route P(33+j) / P(65+j),
  64/64, and a 64-port matrix read with N = 16 raises; (K4) the knob parsers.
* ``test_c10a_records`` (``FEM_EM_GEO38_C10A=1``): (R-d1) on
  ``logs/wf12/h12c-cu-s32-a2.npz`` and (R-d2) on ``logs/wf12/h16b-cu-s32-d2.npz``;
  negative controls (rule (e)): k = 2 in place of k = 1 misses the d1 record by
  more than ``REDUCTION_BAND`` (asserted, backed by `WF-12` ``:6397``); ring
  swap, one-slot ring offset and a within-ring scramble printed against their
  algebraic predictions.
* ``test_geo38_c10_c_mode`` (``FEM_EM_GEO38_C10=1``; C10 runs it after entry
  44): the k = 1 ``C_mode`` of C3's saved 64×64.

Bands, all imported: ``EXACT_IDENTITY_RTOL``, ``REUSE_REPRODUCTION_RTOL``,
``REDUCTION_BAND``, ``RECIPROCITY_BAND``, ``MATRIX_PASSIVITY_CEILING``.  No band,
target, node grid or Done-when is set here.  Not claimed: any F-clinical
``C_mode``, root or tuning figure; any AED value.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

import numpy as np
import pytest

import tests.validation.test_geo38_f_clinical_esr_root as c4
import tests.validation.test_geo38_f_clinical_matrix as c3
from fem_em_solver.ports.circuit import s_to_z, z_to_s
from fem_em_solver.ports.ring_modes import (
    antisymmetric_eigenpair,
    ring_mode_rows,
    ring_port_map_from_names,
)
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
from tests.validation.test_wf7_f_human_b1_quadrature import FREQUENCY_HZ

KNOB_C10A = "FEM_EM_GEO38_C10A"
KNOB_C10 = "FEM_EM_GEO38_C10"
RINGS = ("bottom", "top")
WORKSPACE = Path("/workspace")
H12C_A2_NPZ = WORKSPACE / "logs" / "wf12" / "h12c-cu-s32-a2.npz"
H16B_D2_NPZ = WORKSPACE / "logs" / "wf12" / "h16b-cu-s32-d2.npz"
C10_JSON = WORKSPACE / "logs" / "geo38" / "c10-c-mode-d1.json"
C3_KEYS = ("s64", "order", "azimuth_deg", "z_m", "cells", "n_ranks")

# ---- record comparands (not bands), each with its log line --------------------
# H12c attempt 2: docs/testing/logs/20261006T170222Z_WF-12-step3-a2-W1.log.gz.
C_MODE_CU_D1_RECORD_F = 1.181185799070e-10  # :6366 (also test_ans8_copper_degree2_rows.py:90)
LAMBDA_ANTI_CU_D1_RECORD = {  # :6361 (k=0 lambda_1), :6362 (k=1 lambda_0), :6364 (k=2 lambda_0)
    0: 2.431344655e-01 + 1.481624570e01j,
    1: 3.131486739e-01 + 2.105338540e01j,
    2: 2.555062808e-01 + 3.155133217e01j,
}
OFF_BLOCK_CU_D1_RECORD = 4.195925e-03  # :6366
C2_OVER_CMODE_CU_D1_RECORD = 0.332726  # :6397, |C_2/C_mode - 1|
# H25 W1: docs/testing/logs/20261018T095559Z_ANS-8-H25-W1.log:171 (transcribed comparand, never a band).
C_MODE_CU_D2_COMPARAND_F = 1.145221310828e-10

# ---- predictions (printed, never asserted — rule (e)) -------------------------
PRED_SWAP_REL = 1e-15
PRED_SCRAMBLE_OFF_BLOCK_FLOOR = 1e-1
PRED_C3_P33_P41_SLOTS = ((0, "bottom"), (8, "bottom"))


# ---------------------------------------------------------------------------
# Knobs.
# ---------------------------------------------------------------------------
def _knob(name):
    """unset / empty ⇒ False; ``1`` ⇒ True; anything else ⇒ ``ValueError``."""
    raw = os.environ.get(name)
    if raw is None or raw.strip() == "":
        return False
    if raw.strip() != "1":
        raise ValueError(f"{name} must be unset or '1', got {raw!r}")
    return True


def c10a_enabled():
    return _knob(KNOB_C10A)


def c10_enabled():
    return _knob(KNOB_C10)


# ---------------------------------------------------------------------------
# The N-slot reader.
# ---------------------------------------------------------------------------
def dft_k_domain(z, n_slots, n_rings):
    """``W Z Wᴴ`` with ``W = I_rings ⊗ F``, ``F_kj = e^{−2πi jk/N}/√N`` — a copy of
    ``tests/unit/test_port24_esr_circuit.py:125–131`` (the normalisation
    ``ring_mode_rows`` uses at N = 16)."""
    j = np.arange(n_slots)
    f = np.exp(-2j * np.pi * np.outer(j, j) / n_slots) / np.sqrt(n_slots)
    w = np.kron(np.eye(n_rings), f)
    return w @ np.asarray(z, dtype=np.complex128) @ w.conj().T


def canonical_permutation(port_map, n_slots):
    """``perm`` with ``Z[perm][:, perm]`` in (bottom 0…N−1, top 0…N−1) order;
    each (slot, ring) exactly once, else ``ValueError``."""
    pm = [(int(j), str(r)) for j, r in port_map]
    if len(pm) != 2 * n_slots or len(set(pm)) != 2 * n_slots:
        raise ValueError(f"port_map must name each of the {2 * n_slots} (slot, ring) pairs once")
    bad = [p for p in pm if p[1] not in RINGS or not 0 <= p[0] < n_slots]
    if bad:
        raise ValueError(f"port_map entries outside 0..{n_slots - 1} x {RINGS}: {bad[:4]}")
    index = {p: i for i, p in enumerate(pm)}
    return np.array([index[(j, r)] for r in RINGS for j in range(n_slots)], dtype=int)


def ring_c_modes(s, port_map, *, freq_hz, n_slots, z0):
    """k-mode ``C`` of a two-ring ``2N``-port ring-gap S matrix (lossless).

    Returns ``{"modes": {k: {lambda_anti, anti, c_anti_f, eig}}, "c_mode_f",
    "c_k2_f", "off_block", "n_slots", "freq_hz"}`` for k = 0…N/2."""
    s = np.asarray(s, dtype=np.complex128)
    n = int(n_slots)
    if s.shape != (2 * n, 2 * n):
        raise ValueError(f"S is {s.shape}, expected ({2 * n}, {2 * n}) for n_slots = {n}")
    perm = canonical_permutation(port_map, n)
    z = s_to_z(s, z0)
    zh = dft_k_domain(z[np.ix_(perm, perm)], n, 2)
    mask = np.zeros(zh.shape, dtype=bool)
    blocks = {}
    for k in range(n):
        idx = [k, n + k]
        blocks[k] = zh[np.ix_(idx, idx)].copy()
        mask[np.ix_(idx, idx)] = True
    off = float(np.linalg.norm(np.where(mask, 0.0, zh)) / np.linalg.norm(zh))
    omega = 2.0 * np.pi * float(freq_hz)
    modes = {}
    for k in range(n // 2 + 1):
        lam, anti, both = antisymmetric_eigenpair(blocks[k])
        modes[k] = {
            "lambda_anti": lam,
            "anti": anti,
            "c_anti_f": 1.0 / (omega * lam.imag),
            "eig": both,
        }
    return {
        "n_slots": n,
        "freq_hz": float(freq_hz),
        "off_block": off,
        "modes": modes,
        "c_mode_f": modes[1]["c_anti_f"],
        "c_k2_f": modes[2]["c_anti_f"],
    }


# ---------------------------------------------------------------------------
# Port-map routes.
# ---------------------------------------------------------------------------
def map_from_slot_keys(slot_k, slot_sign):
    """``[(slot, ring)]`` from saved slot keys; ring sign −1 = bottom (the
    semantics of the ans8 script's ``_slot_port_map``)."""
    return [(int(k), "top" if float(sg) > 0 else "bottom") for k, sg in zip(slot_k, slot_sign)]


def map_from_geometry(azimuth_deg, z_m):
    """``[(slot, ring)]`` from measured azimuth / z: C4a's ``_slot_of``; bottom when z < 0."""
    return [
        (c4._slot_of(float(a)), "bottom" if float(z) < 0.0 else "top")
        for a, z in zip(azimuth_deg, z_m)
    ]


def map_from_names_64(names):
    """F-clinical name route: P(33+j) → (j, bottom), P(65+j) → (j, top)."""
    out = []
    for name in names:
        n = int(str(name).strip().upper().lstrip("P"))
        if 33 <= n <= 64:
            out.append((n - 33, "bottom"))
        elif 65 <= n <= 96:
            out.append((n - 65, "top"))
        else:
            raise ValueError(f"{name!r} is not an F-clinical ring-gap port name (P33…P96)")
    return out


# ---------------------------------------------------------------------------
# Synthetic closed forms (tests/unit/test_ring_mode_rows.py U2 pattern, N a parameter).
# ---------------------------------------------------------------------------
def _sym_blocks(n, omega, c_mode):
    """Mirror- and ring-symmetric 2×2 blocks; k = 1's antisymmetric λ has
    Im = 1/(ω c_mode); |k| = 2's has Im = 30 + 22."""
    blocks = {}
    for k in range(n):
        kk = min(k, n - k)
        x_anti = 1.0 / (omega * c_mode) if kk == 1 else 30.0 + 11.0 * kk
        r = 1.5 + 0.2 * kk
        a_minus_b = r + 1j * x_anti
        a_plus_b = (r + 0.3) - 1j * (150.0 + 40.0 * kk)
        a, b = (a_plus_b + a_minus_b) / 2, (a_plus_b - a_minus_b) / 2
        blocks[k] = np.array([[a, b], [b, a]])
    return blocks


def _block_circulant(blocks, n):
    """Canonical-order Z = Wᴴ blockdiag(Ẑ_k) W."""
    zh = np.zeros((2 * n, 2 * n), dtype=np.complex128)
    for k, blk in blocks.items():
        idx = [k, n + k]
        zh[np.ix_(idx, idx)] = blk
    j = np.arange(n)
    w = np.kron(np.eye(2), np.exp(-2j * np.pi * np.outer(j, j) / n) / np.sqrt(n))
    return w.conj().T @ zh @ w


def _canonical_map(n):
    return [(j, r) for r in RINGS for j in range(n)]


def _shuffled(z, n, seed):
    perm = np.random.default_rng(seed).permutation(2 * n)
    pm = _canonical_map(n)
    return z[np.ix_(perm, perm)], [pm[p] for p in perm]


def _i16_worst(local, ref):
    """Worst relative difference between the local reader and ``ring_mode_rows``
    over λ_anti k = 0…8, ``c_mode_f``, ``c_k2_f`` and ``off_block``."""
    worst = 0.0
    for k in range(9):
        a = local["modes"][k]["lambda_anti"]
        b = complex(*ref["modes"][str(k)]["lambda_anti"])
        worst = max(worst, abs(a - b) / abs(b))
    for key in ("c_mode_f", "c_k2_f", "off_block"):
        worst = max(worst, abs(local[key] / ref[key] - 1.0))
    return worst


def _ref16(s, port_map, c_f):
    return ring_mode_rows(
        s,
        port_map=port_map,
        freq_hz=FREQUENCY_HZ,
        c_f=c_f,
        feeds=[(0, "bottom"), (4, "bottom")],
        drives={"x": (1.0, 0.0)},
        z0=REFERENCE_IMPEDANCE_OHM,
    )


def _load(path):
    with np.load(path, allow_pickle=False) as z:
        out = {k: z[k] for k in z.files}
    out["order"] = [str(x) for x in out["order"]]
    return out


def _say(msg):
    print(f"[GEO-38 C10a] {msg}", flush=True)


# --------------------------------------------------------------------------- W0
def test_c10a_coding(monkeypatch):
    """Ungated: (K1) closed form N = 32 / 16; (K2) = ring_mode_rows at N = 16;
    (K3) geometry map = name route, 64/64, and N = 16 on 64 ports raises; (K4) knobs."""
    omega = 2.0 * np.pi * FREQUENCY_HZ
    z0 = REFERENCE_IMPEDANCE_OHM

    # (K1) closed form, ports shuffled.
    c0 = 1.0731e-10
    c2_closed = 1.0 / (omega * (30.0 + 22.0))
    k1 = {}
    for n in (32, 16):
        zs, pm = _shuffled(_block_circulant(_sym_blocks(n, omega, c0), n), n, seed=581 + n)
        r = ring_c_modes(z_to_s(zs, z0), pm, freq_hz=FREQUENCY_HZ, n_slots=n, z0=z0)
        k1[n] = (
            abs(r["c_mode_f"] / c0 - 1.0),
            abs(r["c_k2_f"] / c2_closed - 1.0),
            r["modes"][1]["anti"],
            r["off_block"],
            len(r["modes"]),
        )
        _say(
            f"(K1) N={n}: C_mode rel {k1[n][0]:.3e}, C_2 rel {k1[n][1]:.3e} (band "
            f"EXACT_IDENTITY_RTOL {EXACT_IDENTITY_RTOL:g}); anti(k=1) {k1[n][2]:.15f}; "
            f"off-block {k1[n][3]:.3e}; modes k=0..{k1[n][4] - 1}"
        )

    # (K2) the local reader equals ring_mode_rows at N = 16 (a C16-broken Z, off-block O(1e-2)).
    rng = np.random.default_rng(5810)
    z16 = _block_circulant(_sym_blocks(16, omega, c0), 16)
    a = rng.normal(size=(32, 32)) + 1j * rng.normal(size=(32, 32))
    z16 = z16 + 0.5 * (a + a.T)
    zs16, pm16 = _shuffled(z16, 16, seed=17)
    s16 = z_to_s(zs16, z0)
    loc = ring_c_modes(s16, pm16, freq_hz=FREQUENCY_HZ, n_slots=16, z0=z0)
    ref = _ref16(s16, pm16, C_MODE_CU_D1_RECORD_F)
    k2 = _i16_worst(loc, ref)
    _say(
        f"(K2) N=16 local vs ring_mode_rows: worst rel {k2:.3e} over lambda_anti k=0..8, C_mode, "
        f"C_2, off-block ({loc['off_block']:.4e}) (band {EXACT_IDENTITY_RTOL:g}; predicted 0)"
    )

    # (K3) geometry map on C3a's ideal layout = the name route, 64/64; N = 16 on 64 ports raises.
    az, zc = c3._ideal_layout()
    names = [f"P{i}" for i in range(33, 97)]
    geo = map_from_geometry(az, zc)
    nm = map_from_names_64(names)
    agree = sum(g == m for g, m in zip(geo, nm))
    unique = len(set(geo)) == 64 and len(geo) == 64
    _say(
        f"(K3) geometry map vs name route {agree}/64, unique {unique}; P33 -> {geo[0]}, "
        f"P41 -> {geo[8]}, P65 -> {geo[32]}"
    )
    s64 = c3._synthetic_s(az, zc)
    syn32 = ring_c_modes(s64, geo, freq_hz=FREQUENCY_HZ, n_slots=32, z0=z0)
    _say(f"(K3) C3a synthetic 64x64 read at N=32: off-block {syn32['off_block']:.3e} (printed)")
    with pytest.raises(ValueError):
        ring_c_modes(s64, geo, freq_hz=FREQUENCY_HZ, n_slots=16, z0=z0)
    with pytest.raises(ValueError):
        ring_c_modes(s64, geo[:-1] + [geo[0]], freq_hz=FREQUENCY_HZ, n_slots=32, z0=z0)
    with pytest.raises(ValueError):
        map_from_names_64(["P97"])

    # (K4) knob parsers.
    for knob, fn in ((KNOB_C10A, c10a_enabled), (KNOB_C10, c10_enabled)):
        monkeypatch.delenv(knob, raising=False)
        assert fn() is False
        monkeypatch.setenv(knob, "")
        assert fn() is False
        monkeypatch.setenv(knob, "1")
        assert fn() is True
        for bad in ("0", "yes", "2", "build"):
            monkeypatch.setenv(knob, bad)
            with pytest.raises(ValueError):
                fn()
        monkeypatch.delenv(knob, raising=False)
    _say("(K4) knob parsers: unset/empty -> skip, '1' -> run, else ValueError")

    for n in (32, 16):
        assert k1[n][0] <= EXACT_IDENTITY_RTOL and k1[n][1] <= EXACT_IDENTITY_RTOL, k1[n]
        assert k1[n][2] >= 1.0 - EXACT_IDENTITY_RTOL, k1[n]
        assert k1[n][3] <= EXACT_IDENTITY_RTOL, k1[n]
        assert k1[n][4] == n // 2 + 1
    assert k2 <= EXACT_IDENTITY_RTOL, k2
    assert agree == 64 and unique


# --------------------------------------------------------------------------- W1
def test_c10a_records():
    if not c10a_enabled():
        pytest.skip(f"{KNOB_C10A} unset (reads gitignored logs/wf12 npz)")
    t0 = time.perf_counter()
    z0 = REFERENCE_IMPEDANCE_OHM
    for p in (H12C_A2_NPZ, H16B_D2_NPZ):
        assert p.exists(), f"{p} absent: never re-solve"

    # ---- (R-d1) degree-1 copper, slot-key map ----
    a2 = _load(H12C_A2_NPZ)
    s1 = np.array(a2["s32"], dtype=np.complex128)
    pm_slot = map_from_slot_keys(a2["slot_k"], a2["slot_sign"])
    pm_name = ring_port_map_from_names(a2["order"])
    agree1 = sum(a == b for a, b in zip(pm_slot, pm_name))
    _say(
        f"(R-d1) {H12C_A2_NPZ.name}: shape {s1.shape}; slot-key map vs names {agree1}/32, "
        f"unique {len(set(pm_slot))}"
    )
    r1 = ring_c_modes(s1, pm_slot, freq_hz=FREQUENCY_HZ, n_slots=16, z0=z0)
    c1 = r1["c_mode_f"]
    rel_key = abs(c1 / float(a2["c_mode"]) - 1.0)
    rel_rec = abs(c1 / C_MODE_CU_D1_RECORD_F - 1.0)
    _say(
        f"(R-d1) C_mode {c1:.12e} F vs npz c_mode rel {rel_key:.3e}, vs record "
        f"{C_MODE_CU_D1_RECORD_F:.12e} rel {rel_rec:.3e} (band REUSE_REPRODUCTION_RTOL "
        f"{REUSE_REPRODUCTION_RTOL:g}; H25 read 1.776e-15 / 3.515e-13)"
    )
    for k in (0, 1, 2):
        lam = r1["modes"][k]["lambda_anti"]
        rec = LAMBDA_ANTI_CU_D1_RECORD[k]
        _say(
            f"(R-d1)   k={k} lambda_anti {lam.real:+.9e} {lam.imag:+.9e}j Ohm anti "
            f"{r1['modes'][k]['anti']:.6f} C {r1['modes'][k]['c_anti_f'] * 1e12:+.4f} pF beside "
            f"WF-12 {rec.real:+.9e} {rec.imag:+.9e}j (rel {abs(lam - rec) / abs(rec):.2e}, printed)"
        )
    _say(f"(R-d1) off-block {r1['off_block']:.6e} beside WF-12 :6366 {OFF_BLOCK_CU_D1_RECORD:.6e}")
    i16_d1 = _i16_worst(r1, _ref16(s1, pm_name, c1))
    _say(
        f"(R-d1) (I16) local vs ring_mode_rows worst rel {i16_d1:.3e} (band {EXACT_IDENTITY_RTOL:g})"
    )

    # ---- negative controls on the d1 record (rule (e)) ----
    k2_miss = abs(r1["c_k2_f"] / C_MODE_CU_D1_RECORD_F - 1.0)
    _say(
        f"NEG (asserted) k=2 in place of k=1: C_2 {r1['c_k2_f']:.6e} misses the d1 record by "
        f"{k2_miss:.6f} (> REDUCTION_BAND {REDUCTION_BAND:g}; WF-12 :6397 {C2_OVER_CMODE_CU_D1_RECORD})"
    )
    swap = [(j, "top" if r == "bottom" else "bottom") for j, r in pm_slot]
    rs = ring_c_modes(s1, swap, freq_hz=FREQUENCY_HZ, n_slots=16, z0=z0)
    swap_c = abs(rs["c_mode_f"] / c1 - 1.0)
    swap_l = max(
        abs(rs["modes"][k]["lambda_anti"] - r1["modes"][k]["lambda_anti"])
        / abs(r1["modes"][k]["lambda_anti"])
        for k in range(9)
    )
    _say(
        f"NEG (printed) ring swap: C_mode rel {swap_c:.3e}, lambda_anti worst rel {swap_l:.3e} "
        f"(predicted <= {PRED_SWAP_REL:g})"
    )
    for n_slots in (16,):
        off1 = [((j + 1) % n_slots, r) if r == "top" else (j, r) for j, r in pm_slot]
        ro = ring_c_modes(s1, off1, freq_hz=FREQUENCY_HZ, n_slots=n_slots, z0=z0)
        _say(
            f"NEG (printed) top ring offset one slot: C_mode rel {abs(ro['c_mode_f'] / c1 - 1.0):.3e} "
            f"(predicted unchanged); anti(k=1) {ro['modes'][1]['anti']:.6f} vs "
            f"{r1['modes'][1]['anti']:.6f} (predicted cos^2(pi/N) = {np.cos(np.pi / n_slots) ** 2:.4f}; "
            f"N=32 would give {np.cos(np.pi / 32) ** 2:.4f}) — C_mode cannot validate the port map"
        )
    rng = np.random.default_rng(581)
    bottom_slots = rng.permutation(16)
    scr = [(int(bottom_slots[j]), r) if r == "bottom" else (j, r) for j, r in pm_slot]
    rsc = ring_c_modes(s1, scr, freq_hz=FREQUENCY_HZ, n_slots=16, z0=z0)
    _say(
        f"NEG (printed) fixed-seed scramble within the bottom ring: off-block {rsc['off_block']:.4e} "
        f"vs {r1['off_block']:.4e} (predicted >= {PRED_SCRAMBLE_OFF_BLOCK_FLOOR:g})"
    )

    # ---- (R-d2) degree-2 copper, names map ----
    d2 = _load(H16B_D2_NPZ)
    s2 = np.array(d2["s32"], dtype=np.complex128)
    pm2 = ring_port_map_from_names(d2["order"])
    r2 = ring_c_modes(s2, pm2, freq_hz=FREQUENCY_HZ, n_slots=16, z0=z0)
    c2 = r2["c_mode_f"]
    rel_d2 = abs(c2 / C_MODE_CU_D2_COMPARAND_F - 1.0)
    _say(
        f"(R-d2) {H16B_D2_NPZ.name}: shape {s2.shape}; C_mode {c2:.12e} F vs H25 W1 :171 "
        f"{C_MODE_CU_D2_COMPARAND_F:.12e} rel {rel_d2:.3e} (band {REUSE_REPRODUCTION_RTOL:g}); "
        f"anti(k=1) {r2['modes'][1]['anti']:.6f}; off-block {r2['off_block']:.6e}"
    )
    i16_d2 = _i16_worst(r2, _ref16(s2, pm2, c2))
    _say(
        f"(R-d2) (I16) local vs ring_mode_rows worst rel {i16_d2:.3e} (band {EXACT_IDENTITY_RTOL:g})"
    )
    _say(f"test body {time.perf_counter() - t0:.2f} s")

    assert agree1 == 32 and len(set(pm_slot)) == 32, (agree1, pm_slot)
    assert i16_d1 <= EXACT_IDENTITY_RTOL, i16_d1
    assert rel_key <= REUSE_REPRODUCTION_RTOL and rel_rec <= REUSE_REPRODUCTION_RTOL, (
        rel_key,
        rel_rec,
    )
    assert k2_miss > REDUCTION_BAND, k2_miss
    assert i16_d2 <= EXACT_IDENTITY_RTOL, i16_d2
    assert rel_d2 <= REUSE_REPRODUCTION_RTOL, rel_d2


# --------------------------------------------------------------------------- C10
def test_geo38_c10_c_mode():
    """C10 (entry 44's npz): the k = 1 ``C_mode`` of C3's saved 64×64.  Lands skipped."""
    if not c10_enabled():
        pytest.skip(f"{KNOB_C10} unset (C10 runs after xl-pending entry 44)")
    from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
    from tests.validation.test_port_birdcage_ring_matrix import MATRIX_PASSIVITY_CEILING
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    t0 = time.perf_counter()
    path = WORKSPACE / c3.SAVE_PATH
    assert path.exists(), f"{path} absent: entry 44's matrix; never re-solve"
    with np.load(path, allow_pickle=False) as z:
        sv = {k: z[k] for k in z.files}
    missing = [k for k in C3_KEYS if k not in sv]
    assert not missing, f"C3 npz lacks {missing}"
    s = np.array(sv["s64"], dtype=np.complex128)
    order = [str(x) for x in sv["order"]]
    assert s.shape == (64, 64), s.shape
    assert order == [f"P{i}" for i in range(33, 97)], order[:4]
    geo = map_from_geometry(sv["azimuth_deg"], sv["z_m"])
    assert len(geo) == 64 and len(set(geo)) == 64, sorted(geo)
    agree = sum(g == m for g, m in zip(geo, map_from_names_64(order)))
    ratio = _reciprocity_ratio(s)
    sigma_max = float(np.linalg.norm(s, 2))
    r = ring_c_modes(s, geo, freq_hz=FREQUENCY_HZ, n_slots=32, z0=REFERENCE_IMPEDANCE_OHM)
    lam1 = r["modes"][1]["lambda_anti"]
    print(
        f"[GEO-38 C10] {path} ({path.stat().st_size} B); geometry map unique 64/64; name-route "
        f"agreement {agree}/64 (predicted 64/64); P33 -> {geo[0]}, P41 -> {geo[8]} (predicted "
        f"{PRED_C3_P33_P41_SLOTS}); cells {int(sv['cells'])}, n_ranks {int(sv['n_ranks'])}",
        flush=True,
    )
    print(
        f"[GEO-38 C10] reciprocity {ratio:.6e} (band {RECIPROCITY_BAND:.0e}); sigma_max "
        f"{sigma_max:.12f} (ceiling {MATRIX_PASSIVITY_CEILING})",
        flush=True,
    )
    for k in range(17):
        m = r["modes"][k]
        lam = m["lambda_anti"]
        print(
            f"[GEO-38 C10]   k={k:2d} lambda_anti {lam.real:+.9e} {lam.imag:+.9e}j Ohm anti "
            f"{m['anti']:.6f} C {m['c_anti_f'] * 1e12:+.4f} pF",
            flush=True,
        )
    print(
        f"[GEO-38 C10] C_mode {r['c_mode_f']:.12e} F = {r['c_mode_f'] * 1e12:.4f} pF; C_2 "
        f"{r['c_k2_f']:.12e} F; anti(k=1) {r['modes'][1]['anti']:.6f}; off-block "
        f"{r['off_block']:.6e}; body {time.perf_counter() - t0:.2f} s",
        flush=True,
    )
    C10_JSON.parent.mkdir(parents=True, exist_ok=True)
    C10_JSON.write_text(
        json.dumps(
            {
                "source": str(path),
                "freq_hz": FREQUENCY_HZ,
                "z0_ohm": REFERENCE_IMPEDANCE_OHM,
                "n_slots": 32,
                "c_mode_f": r["c_mode_f"],
                "c_k2_f": r["c_k2_f"],
                "off_block": r["off_block"],
                "reciprocity": ratio,
                "sigma_max": sigma_max,
                "name_route_agreement": agree,
                "cells": int(sv["cells"]),
                "n_ranks": int(sv["n_ranks"]),
                "modes": {
                    str(k): {
                        "lambda_anti": [m["lambda_anti"].real, m["lambda_anti"].imag],
                        "anti": m["anti"],
                        "c_anti_f": m["c_anti_f"],
                    }
                    for k, m in r["modes"].items()
                },
            },
            indent=2,
        )
    )
    assert ratio <= RECIPROCITY_BAND, ratio
    assert sigma_max <= MATRIX_PASSIVITY_CEILING, sigma_max
    assert lam1.imag > 0.0, lam1
