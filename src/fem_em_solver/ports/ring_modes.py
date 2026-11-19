"""Ring-mode rows of a 16-leg birdcage's 32-port gap matrix (`ANS-8` step 1).

Pure NumPy (no dolfinx): the one function both halves of `ANS-8` apply to a
32×32 S matrix — ours (``ans:8`` stage R) and AED's exported ``"S"`` (the
private writer) — so the two codes' rows cannot drift apart.  The algebra is
SPEC rows 1–3 (``examples/ansys_benchmarks/ans8_…/SPEC.md``):

* ``Z = z0 (I + S)(I − S)⁻¹`` (`circuit.s_to_z`), ``Y_net = Z⁻¹``;
* the k-blocks ``Ẑ_k`` (2×2, bottom / top ring) of the unitary per-ring DFT
  ``Ẑ = W Z Wᴴ``, ``W = I₂ ⊗ F``, ``F_kj = e^{−2πi jk/16} / 4`` (for an
  N-slot ring, ``n_slots=N``: ``F_kj = e^{−2πi jk/N}/√N``); the off-block
  norm ``‖Ẑ − blockdiag‖_F / ‖Ẑ‖_F`` (C16 breaking); ``C_k = 1/(ω Im λ)`` of the
  eigenpair of ``Ẑ_k`` whose eigenvector is ring-antisymmetric;
* every non-feed gap terminated in ``jωC``, each feed a Norton source
  ``I_s = V_f / z0`` in parallel with ``1/z0 + jωC``; ``(Y_net + Y_load) V =
  I_s``, gap currents ``I = Y_net V``;
* per ring ``ĝ_k = Σ_j I_j e^{−2πi jk/16}`` (NumPy ``fft`` bins): the k = ±1
  share (bins 1 and 15), the k = 0 leak, all 16 fractions and the |k| folds;
* feed rows: ``V_f / I_f`` per feed and the 2×2 feed admittance matrix of the
  network with every other gap in C.

SPEC v2 (§9 item 161): :func:`ring_mode_rows_v2` composes the above into
rows 2b (rows 2–3 at the matrix's own row-1 ``C_mode``), 2c (symmetric-drive
controls at both C) and 4a (``f₀ = f·√(C_mode / C)``), with no new solve.

Ports are identified by ``(slot j, ring)`` with ring ``"bottom"`` or ``"top"``
— :func:`ring_port_map_from_names` maps the SPEC names P(17 + j) → (j,
bottom), P(33 + j) → (j, top).
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

from .circuit import s_to_z

N_SLOTS = 16
RINGS = ("bottom", "top")


def _check_n_slots(n_slots) -> int:
    """``n_slots`` (`WF-16` step 13, §9 item 588): an ``int`` (not ``bool``),
    even, ≥ 4 (k = 2 is read); a ``float`` is refused, never cast."""
    if isinstance(n_slots, bool) or not isinstance(n_slots, (int, np.integer)):
        raise TypeError(f"n_slots: must be an int, got {type(n_slots).__name__} {n_slots!r}")
    if n_slots < 4 or n_slots % 2:
        raise ValueError(f"n_slots: must be an even int >= 4, got {n_slots!r}")
    return int(n_slots)


def ring_port_map_from_names(names: Sequence[str], n_slots: int = N_SLOTS) -> list:
    """``[(slot j, ring)]`` for SPEC port names ``P17…P48`` in ``names`` order.

    ``n_slots=N`` (`WF-16` step 13): P(N+1+j) → (j, bottom), P(2N+1+j) → (j,
    top); the default 16 is P17…P48 as before."""
    n_s = _check_n_slots(n_slots)
    lo, mid, hi = n_s + 1, 2 * n_s + 1, 3 * n_s
    out = []
    for name in names:
        n = int(str(name).strip().upper().lstrip("P"))
        if lo <= n <= mid - 1:
            out.append((n - lo, "bottom"))
        elif mid <= n <= hi:
            out.append((n - mid, "top"))
        else:
            raise ValueError(f"{name!r} is not a ring-gap port name (P{lo}…P{hi}, n_slots={n_s})")
    return out


def _canonical_permutation(port_map, n_slots: int = N_SLOTS) -> np.ndarray:
    """``perm`` with ``Z[perm][:, perm]`` in (bottom 0…N−1, top 0…N−1) order."""
    n_s = _check_n_slots(n_slots)
    pm = [(int(j), str(r)) for j, r in port_map]
    if len(pm) != 2 * n_s or len(set(pm)) != 2 * n_s:
        raise ValueError(
            f"port_map must name each of the {2 * n_s} (slot, ring) pairs once "
            f"(n_slots={n_s}): {pm}"
        )
    index = {p: i for i, p in enumerate(pm)}
    try:
        return np.array([index[(j, r)] for r in RINGS for j in range(n_s)], dtype=int)
    except KeyError as exc:
        raise ValueError(f"port_map lacks {exc.args[0]}") from exc


def _dft_matrix(n_slots: int = N_SLOTS) -> np.ndarray:
    j = np.arange(n_slots)
    f = np.exp(-2j * np.pi * np.outer(j, j) / n_slots) / np.sqrt(n_slots)
    return np.kron(np.eye(2), f)


def k_blocks(z_canonical: np.ndarray, n_slots: int = N_SLOTS):
    """``({k: Ẑ_k (2×2)}, off-block norm)`` of a canonical-order 2N×2N ``Z``."""
    n_s = _check_n_slots(n_slots)
    w = _dft_matrix(n_s)
    zh = w @ np.asarray(z_canonical, dtype=np.complex128) @ w.conj().T
    blocks = {}
    mask = np.zeros(zh.shape, dtype=bool)
    for k in range(n_s):
        idx = [k, n_s + k]
        blocks[k] = zh[np.ix_(idx, idx)].copy()
        mask[np.ix_(idx, idx)] = True
    off = float(np.linalg.norm(np.where(mask, 0.0, zh)) / np.linalg.norm(zh))
    return blocks, off


def antisymmetric_eigenpair(block: np.ndarray):
    """``(λ, anti fraction, both (λ, anti))`` — the ring-antisymmetric eigenpair."""
    lam, vec = np.linalg.eig(np.asarray(block, dtype=np.complex128))
    pairs = []
    for i in range(2):
        v = vec[:, i]
        anti = float(abs(v[0] - v[1]) ** 2 / (2.0 * (abs(v[0]) ** 2 + abs(v[1]) ** 2)))
        pairs.append((complex(lam[i]), anti))
    best = max(pairs, key=lambda p: p[1])
    return best[0], best[1], pairs


def k_spectrum(ring_currents: Sequence[complex], n_slots: int = N_SLOTS) -> dict:
    """NumPy-``fft`` k-fractions of one ring's N (default 16) slot currents."""
    n_s = _check_n_slots(n_slots)
    a = np.asarray(ring_currents, dtype=np.complex128)
    if a.shape != (n_s,):
        raise ValueError(f"n_slots: ring_currents has shape {a.shape}, expected ({n_s},)")
    g = np.fft.fft(a)
    p = np.abs(g) ** 2
    total = float(np.sum(p))
    frac = p / total
    folds = [float(frac[0])] + [float(frac[k] + frac[n_s - k]) for k in range(1, n_s // 2)]
    folds.append(float(frac[n_s // 2]))
    return {
        "share_k1": float(frac[1] + frac[n_s - 1]),
        "leak_k0": float(frac[0]),
        "fractions": [float(x) for x in frac],
        "abs_k_folds": folds,
    }


def _esr_canonical(esr_ohm, perm):
    """``None`` for the lossless default (``esr_ohm`` a scalar 0), else the
    per-gap series resistance in canonical order (`PORT-24` step 1).

    ``esr_ohm`` — a scalar (every gap) or one value per row of ``s`` (the
    caller's port order, permuted here)."""
    if np.ndim(esr_ohm) == 0:
        if float(esr_ohm) == 0.0:
            return None
        return np.full(len(perm), float(esr_ohm))
    r = np.asarray(esr_ohm, dtype=float)
    if r.shape != (len(perm),):
        raise ValueError(f"esr_ohm must be a scalar or {len(perm)} values, got shape {r.shape}")
    return r[perm].copy()


def _gap_admittance(jwc, n, r_canonical):
    """Per-gap termination admittance: ``jωC`` (lossless, the old expression,
    bitwise) or ``1/(r + 1/(jωC))`` (series RC, capacitor ESR)."""
    if r_canonical is None:
        return np.full(n, jwc, dtype=np.complex128)
    return 1.0 / (np.asarray(r_canonical, dtype=float) + 1.0 / jwc)


def terminated_currents(z_canonical, *, freq_hz, c_f, feed_idx, v_feeds, z0=50.0, esr_ohm=None):
    """``(V, I)`` of the Norton-terminated network, canonical order.

    ``esr_ohm`` (`PORT-24` step 1, additive): ``None`` / 0 is the lossless
    ``jωC`` load bitwise; else a scalar or a canonical-order vector of series
    resistances, each gap loaded by ``1/(r + 1/(jωC))``.
    """
    z = np.asarray(z_canonical, dtype=np.complex128)
    y_net = np.linalg.inv(z)
    jwc = 1j * 2.0 * np.pi * float(freq_hz) * float(c_f)
    if esr_ohm is None or (np.ndim(esr_ohm) == 0 and float(esr_ohm) == 0.0):
        y_load = np.full(z.shape[0], jwc, dtype=np.complex128)
    else:
        r = np.broadcast_to(np.asarray(esr_ohm, dtype=float), (z.shape[0],))
        y_load = _gap_admittance(jwc, z.shape[0], r)
    i_s = np.zeros(z.shape[0], dtype=np.complex128)
    for f, v in zip(feed_idx, v_feeds):
        y_load[f] += 1.0 / z0
        i_s[f] = complex(v) / z0
    v = np.linalg.solve(y_net + np.diag(y_load), i_s)
    return v, y_net @ v


def ring_mode_rows(
    s, *, port_map, freq_hz, c_f, feeds, drives, z0=50.0, esr_ohm=0.0, n_slots=N_SLOTS
) -> dict:
    """SPEC rows 1–3 of a 32-port ring-gap S matrix.

    ``s`` — 32×32, in any port order; ``port_map`` — ``(slot j, ring)`` per row
    of ``s``; ``feeds`` — the feeds as ``(slot, ring)`` pairs; ``drives`` —
    ``{label: feed voltages}`` (one complex per feed).  Returns a JSON-ready
    dict (complex numbers as ``[re, im]``).

    ``esr_ohm`` (`PORT-24` step 1, operator ESR ruling 2026-10-06; additive,
    default 0 = the lossless rows bitwise): a series resistance on every gap
    capacitor — a scalar, or one value per row of ``s``.  It enters the loop
    impedance, ``λ_k`` from ``Z + diag(r)`` (so ``C_mode = 1/(ω Im λ₁)`` sees
    the lossy loop), and every termination, ``1/(r + 1/(jωC))``.

    ``n_slots`` (`WF-16` step 13, additive, default 16 = the rows above
    bitwise): an N-slot ring pair, ``s`` 2N×2N, modes k = 0…N/2 with
    ``F_kj = e^{−2πi jk/N}/√N``.
    """
    n_s = _check_n_slots(n_slots)
    perm = _canonical_permutation(port_map, n_s)
    z = s_to_z(np.asarray(s, dtype=np.complex128), z0)
    zc = z[np.ix_(perm, perm)]
    omega = 2.0 * np.pi * float(freq_hz)
    r_c = _esr_canonical(esr_ohm, perm)
    blocks, off = k_blocks(zc if r_c is None else zc + np.diag(r_c), n_s)
    modes = {}
    for k in range(n_s // 2 + 1):
        lam, anti, both = antisymmetric_eigenpair(blocks[k])
        modes[k] = {
            "lambda_anti": [lam.real, lam.imag],
            "anti": anti,
            "c_anti_f": 1.0 / (omega * lam.imag),
            "eig": [[l.real, l.imag, a] for l, a in both],
        }
    canon = {(j, r): RINGS.index(r) * n_s + j for r in RINGS for j in range(n_s)}
    feed_idx = [canon[(int(j), str(r))] for j, r in feeds]
    jwc = 1j * omega * float(c_f)
    out_drives = {}
    for label, v_feeds in drives.items():
        v_feeds = [complex(x) for x in v_feeds]
        if len(v_feeds) != len(feed_idx):
            raise ValueError(
                f"drive {label!r} has {len(v_feeds)} voltages for {len(feed_idx)} feeds"
            )
        v, i = terminated_currents(
            zc, freq_hz=freq_hz, c_f=c_f, feed_idx=feed_idx, v_feeds=v_feeds, z0=z0, esr_ohm=r_c
        )
        rings = {r: k_spectrum(i[n * n_s : (n + 1) * n_s], n_s) for n, r in enumerate(RINGS)}
        feed_rows = []
        for f in feed_idx:
            zf = v[f] / i[f] if i[f] != 0 else complex("nan")
            y_c = jwc if r_c is None else 1.0 / (r_c[f] + 1.0 / jwc)
            yin = (i[f] + y_c * v[f]) / v[f] if v[f] != 0 else complex("nan")
            feed_rows.append(
                {
                    "v": [v[f].real, v[f].imag],
                    "i": [i[f].real, i[f].imag],
                    "z_port_ohm": [zf.real, zf.imag],
                    "y_in_with_c_s": [yin.real, yin.imag],
                }
            )
        out_drives[label] = {
            "v_feeds": [[x.real, x.imag] for x in v_feeds],
            "rings": rings,
            "feeds": feed_rows,
        }
    # Feed 2×2: the network admittance at the feeds, every other gap in C.
    y_net = np.linalg.inv(zc)
    others = [n for n in range(2 * n_s) if n not in feed_idx]
    if r_c is None:
        y_oo = y_net[np.ix_(others, others)] + jwc * np.eye(len(others))
    else:
        y_oo = y_net[np.ix_(others, others)] + np.diag(
            _gap_admittance(jwc, len(others), r_c[others])
        )
    y_ff = y_net[np.ix_(feed_idx, feed_idx)] - y_net[np.ix_(feed_idx, others)] @ np.linalg.solve(
        y_oo, y_net[np.ix_(others, feed_idx)]
    )
    return {
        "freq_hz": float(freq_hz),
        "c_f": float(c_f),
        "z0_ohm": float(z0),
        "feeds": [[int(j), str(r)] for j, r in feeds],
        "off_block": off,
        "modes": {str(k): m for k, m in modes.items()},
        "c_mode_f": modes[1]["c_anti_f"],
        "c_k2_f": modes[2]["c_anti_f"],
        "drives": out_drives,
        "y_feed_2x2_s": [[[complex(x).real, complex(x).imag] for x in row] for row in y_ff],
    }


def ring_mode_rows_v2(
    s, *, port_map, freq_hz, c_spec_f, feeds, drives, controls, z0=50.0, esr_ohm=0.0
) -> dict:
    """SPEC v2 rows 2b / 2c / 4a — composes :func:`ring_mode_rows`, never edits it.

    ``fixed_c`` — rows 1–3 at ``c_spec_f`` (v1's rows, bit-identical);
    ``own_c_mode`` — the same drives with every gap (feeds included) at the
    matrix's **own** row-1 ``C_mode``, wired from ``fixed_c["c_mode_f"]`` and
    never a parameter, so a caller cannot hand one code the other's C_mode
    (row 2b, with row 3 at own C as 3b); ``controls`` — ``{label: (feeds,
    voltages)}`` — each its own ``ring_mode_rows`` call with its own feeds (the
    Norton 1/z0 sits on driven ports only) at both C (row 2c);
    ``f0_one_point_hz = freq_hz·√(C_own / c_spec)`` (row 4a, λ₁ taken ∝ ω).
    ``esr_ohm`` (`PORT-24` step 1) passes through to every inner call; the
    default 0 leaves every row bitwise.
    """
    common = dict(port_map=port_map, freq_hz=freq_hz, z0=z0)
    if not (np.ndim(esr_ohm) == 0 and float(esr_ohm) == 0.0):
        common["esr_ohm"] = esr_ohm
    fixed = ring_mode_rows(s, c_f=c_spec_f, feeds=feeds, drives=drives, **common)
    c_own = fixed["c_mode_f"]
    own = ring_mode_rows(s, c_f=c_own, feeds=feeds, drives=drives, **common)
    ctl = {"fixed_c": {}, "own_c_mode": {}}
    for label, (c_feeds, volts) in controls.items():
        for key, c in (("fixed_c", c_spec_f), ("own_c_mode", c_own)):
            ctl[key][label] = ring_mode_rows(
                s, c_f=c, feeds=c_feeds, drives={label: volts}, **common
            )
    return {
        "fixed_c": fixed,
        "own_c_mode": own,
        "controls": ctl,
        "c_spec_f": float(c_spec_f),
        "c_own_f": float(c_own),
        "anti_k1": fixed["modes"]["1"]["anti"],
        "f0_one_point_hz": float(freq_hz) * float(np.sqrt(c_own / float(c_spec_f))),
    }


__all__ = [
    "N_SLOTS",
    "RINGS",
    "antisymmetric_eigenpair",
    "k_blocks",
    "k_spectrum",
    "ring_mode_rows",
    "ring_mode_rows_v2",
    "ring_port_map_from_names",
    "terminated_currents",
]
