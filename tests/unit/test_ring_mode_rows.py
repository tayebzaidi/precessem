"""`ANS-8` step 1 (§9 item 141): synthetic closed forms for
``fem_em_solver.ports.ring_modes.ring_mode_rows`` — every assertion ≤ 1e-12
unless named otherwise.

(U1) a block-circulant Z under a shuffled port order is recovered, off-block
≈ 0; negative control: a shuffled slot map gives off-block > 0.1.
(U2) ``C_mode = 1/(ω Im λ_anti)`` on a k = 1 block built from a known C.
(U3) D16 network, SPEC loads: four-port k = 0 leak and every even-k bin ≈ 0;
opposite pair even-k ≈ 0; ccw and cw |k|-fold shares equal (mirror symmetry).
(U4) every non-±1 mode open (j·1e9 Ω): four-port share ≥ 1 − 1e-9, leak ≤ 1e-12.
(U5) ideal-source two-port: the dense feed-ring share equals the k-domain
closed form from the per-k 2×2 solves.
(U6) sense: drive (ii)'s weights equal ``quadrature_phase_weights([0, 1],
"ccw")``; ccw sits wholly in NumPy bin 15, cw in bin 1.

SPEC v2 (§9 item 161), ``ring_mode_rows_v2``:
(U7) own C_mode wired from row 1 (bitwise), own-C blocks equal a direct call
and are invariant in ``c_spec_f``; U5's closed form at own C.
(U8) C16 (mirror-broken) selection rules for both controls at both C; open
others ⇒ control share → 1.
(U9) gauge negative control: a single-gap defect gives a control leak > 1e-9.
(U10) λ₁ = R + jωL ⇒ f₀ = 1/(2π√(LC)) and a rebuild at f₀ gives C; control:
a series-C λ₁ misses its true root by > 1e-3.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from fem_em_solver.ports.circuit import z_to_s
from fem_em_solver.ports.ring_modes import (
    N_SLOTS,
    RINGS,
    k_blocks,
    k_spectrum,
    ring_mode_rows,
    ring_port_map_from_names,
    terminated_currents,
)

REPO = Path(__file__).resolve().parents[2]
SCRIPT = (REPO / "examples" / "ansys_benchmarks" / "ans8_f_human_quadrature_mode_content_64MHz"
          / "08_f_human_quadrature_mode_content.py")
TOL = 1e-12
F = 64.0e6
W = 2 * np.pi * F
C = 1.1311e-10
B = lambda j: (j, "bottom")  # noqa: E731


def _w():
    j = np.arange(N_SLOTS)
    f = np.exp(-2j * np.pi * np.outer(j, j) / N_SLOTS) / np.sqrt(N_SLOTS)
    return np.kron(np.eye(2), f)


def _block_circulant(blocks):
    """Canonical-order Z = Wᴴ blockdiag(Ẑ_k) W."""
    zh = np.zeros((32, 32), dtype=np.complex128)
    for k, blk in blocks.items():
        idx = [k, N_SLOTS + k]
        zh[np.ix_(idx, idx)] = blk
    w = _w()
    return w.conj().T @ zh @ w


def _d16_blocks(c_mode=1.13e-10, open_others=None):
    """Mirror-symmetric (Ẑ_k = Ẑ_{16−k}) ring-symmetric 2×2 blocks."""
    blocks = {}
    for k in range(N_SLOTS):
        kk = min(k, N_SLOTS - k)
        if open_others is not None and kk != 1:
            blocks[k] = open_others * np.eye(2)
            continue
        x_anti = 1.0 / (W * c_mode) if kk == 1 else 30.0 + 11.0 * kk
        r = 1.5 + 0.2 * kk
        a_minus_b = r + 1j * x_anti
        a_plus_b = (r + 0.3) - 1j * (150.0 + 40.0 * kk)
        a, b = (a_plus_b + a_minus_b) / 2, (a_plus_b - a_minus_b) / 2
        blocks[k] = np.array([[a, b], [b, a]])
    return blocks


def _canonical_map():
    return [(j, r) for r in RINGS for j in range(N_SLOTS)]


def _shuffled(z, seed=7):
    rng = np.random.default_rng(seed)
    perm = rng.permutation(32)
    pm = _canonical_map()
    return z[np.ix_(perm, perm)], [pm[p] for p in perm]


def test_u1_block_circulant_recovered_under_shuffle():
    blocks = _d16_blocks()
    blocks[3] = blocks[3] + np.array([[0.1j, 0.2], [0.2, -0.05j]])  # non-symmetric ring pair
    blocks[13] = blocks[3]
    z = _block_circulant(blocks)
    zs, pm = _shuffled(z)
    r = ring_mode_rows(z_to_s(zs, 50.0), port_map=pm, freq_hz=F, c_f=C, feeds=[B(0), B(4)],
                       drives={"x": (1, 0)})
    assert r["off_block"] <= TOL
    # The blocks themselves, through the same permutation the function uses.
    from fem_em_solver.ports.ring_modes import _canonical_permutation
    perm = _canonical_permutation(pm)
    got, off = k_blocks(zs[np.ix_(perm, perm)])
    for k in range(N_SLOTS):
        assert np.max(np.abs(got[k] - blocks[k])) / np.max(np.abs(blocks[k])) <= 1e-12
    # Negative control (asserted): a shuffled slot map breaks C16.
    rng = np.random.default_rng(3)
    wrong = [pm[i] for i in rng.permutation(32)]
    bad = ring_mode_rows(z_to_s(zs, 50.0), port_map=wrong, freq_hz=F, c_f=C,
                         feeds=[B(0), B(4)], drives={"x": (1, 0)})
    print(f"U1 off-block {r['off_block']:.3e}; shuffled-map negative control {bad['off_block']:.3f}")
    assert bad["off_block"] > 0.1


def test_u2_c_mode_closed_form():
    c0 = 1.0731e-10
    z = _block_circulant(_d16_blocks(c_mode=c0))
    r = ring_mode_rows(z_to_s(z, 50.0), port_map=_canonical_map(), freq_hz=F, c_f=C,
                       feeds=[B(0), B(4)], drives={"x": (1, 0)})
    rel = abs(r["c_mode_f"] / c0 - 1.0)
    c2_closed = 1.0 / (W * (30.0 + 22.0))
    rel2 = abs(r["c_k2_f"] / c2_closed - 1.0)
    print(f"U2 C_mode relative {rel:.3e}; C_2 relative {rel2:.3e}")
    assert rel <= TOL and rel2 <= TOL
    assert r["modes"]["1"]["anti"] >= 1.0 - TOL


def _rows(z, feeds, drives):
    return ring_mode_rows(z_to_s(z, 50.0), port_map=_canonical_map(), freq_hz=F, c_f=C,
                          feeds=feeds, drives=drives)


def test_u3_d16_parity_and_mirror():
    z = _block_circulant(_d16_blocks())
    four = _rows(z, [B(0), B(4), B(8), B(12)], {"q": (1, 1j, -1, -1j)})["drives"]["q"]
    opp = _rows(z, [B(0), B(8)], {"o": (1, -1)})["drives"]["o"]
    two = _rows(z, [B(0), B(4)], {"ccw": (1, -1j), "cw": (1, 1j)})["drives"]
    worst4 = max(four["rings"][r]["fractions"][k] for r in RINGS for k in range(0, 16, 2))
    worsto = max(opp["rings"][r]["fractions"][k] for r in RINGS for k in range(0, 16, 2))
    mirror = max(abs(a - b) for r in RINGS for a, b in
                 zip(two["ccw"]["rings"][r]["abs_k_folds"], two["cw"]["rings"][r]["abs_k_folds"]))
    print(f"U3 four-port even-k max {worst4:.3e} leak {four['rings']['bottom']['leak_k0']:.3e}; "
          f"opposite even-k max {worsto:.3e}; |ccw - cw| fold max {mirror:.3e}")
    assert four["rings"]["bottom"]["leak_k0"] <= TOL and worst4 <= TOL
    assert worsto <= TOL
    assert mirror <= TOL


def test_u4_only_k1_closed():
    z = _block_circulant(_d16_blocks(open_others=1j * 1e9))
    four = _rows(z, [B(0), B(4), B(8), B(12)], {"q": (1, 1j, -1, -1j)})["drives"]["q"]
    sh = four["rings"]["bottom"]["share_k1"]
    leak = four["rings"]["bottom"]["leak_k0"]
    print(f"U4 four-port share {sh:.15f} leak {leak:.3e}")
    assert sh >= 1.0 - 1e-9
    assert leak <= TOL


def test_u5_ideal_source_k_domain():
    blocks = _d16_blocks()
    z = _block_circulant(blocks)
    jwc = 1j * W * C
    feed_idx, v = [0, 4], (1.0, 1j)
    big = 1e15  # ideal current sources: the feed's 1/z0 vanishes
    _, i_dense = terminated_currents(z, freq_hz=F, c_f=C, feed_idx=feed_idx, v_feeds=v, z0=big)
    dense = k_spectrum(i_dense[:N_SLOTS])["share_k1"]
    # k-domain: Î_s,k on each ring, per-k 2×2 (Ŷ_k + jωC) V̂ = Î_s, Î = Ŷ_k V̂.
    i_s = np.zeros(32, dtype=np.complex128)
    for f, vv in zip(feed_idx, v):
        i_s[f] = vv / big
    w = _w()
    i_hat_s = w @ i_s
    i_hat = np.zeros(32, dtype=np.complex128)
    for k in range(N_SLOTS):
        idx = [k, N_SLOTS + k]
        y = np.linalg.inv(blocks[k])
        vk = np.linalg.solve(y + jwc * np.eye(2), i_hat_s[idx])
        i_hat[idx] = y @ vk
    p = np.abs(i_hat[:N_SLOTS]) ** 2
    closed = float((p[1] + p[15]) / p.sum())
    print(f"U5 dense {dense:.15f} k-domain {closed:.15f}")
    assert abs(dense - closed) <= TOL


def test_u6_sense():
    from tests.validation.test_birdcage_b1_quadrature import quadrature_phase_weights

    spec = importlib.util.spec_from_file_location("ans8_ut", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["ans8_ut"] = mod
    spec.loader.exec_module(mod)
    assert np.max(np.abs(np.array(mod.SPEC_DRIVES["ccw"]) - quadrature_phase_weights([0, 1], "ccw"))) <= TOL
    assert np.max(np.abs(np.array(mod.SPEC_DRIVES["cw"]) - quadrature_phase_weights([0, 1], "cw"))) <= TOL
    j = np.arange(N_SLOTS)
    ccw = k_spectrum(np.exp(-2j * np.pi * j / N_SLOTS))["fractions"]
    cw = k_spectrum(np.exp(+2j * np.pi * j / N_SLOTS))["fractions"]
    assert abs(ccw[15] - 1.0) <= TOL and abs(cw[1] - 1.0) <= TOL
    # The ccw pattern sampled at the feed slots 0 / 4 is drive (ii)'s weights.
    assert abs(np.exp(-2j * np.pi * 4 / N_SLOTS) - mod.SPEC_DRIVES["ccw"][1]) <= TOL
    # Naming map: P17 → (0, bottom), P21 → (4, bottom), P33 → (0, top).
    assert ring_port_map_from_names(["P17", "P21", "P33", "P48"]) == [
        (0, "bottom"), (4, "bottom"), (0, "top"), (15, "top")]
    with pytest.raises(ValueError):
        ring_port_map_from_names(["P16"])


# ---- SPEC v2 (§9 item 161): ring_mode_rows_v2 — U7–U10, each ≤ 1e-12 ----
from fem_em_solver.ports.ring_modes import ring_mode_rows_v2  # noqa: E402

V2_DRIVES = {"p17_alone": (1.0, 0.0), "ccw": (1.0, -1j), "cw": (1.0, 1j)}
V2_CONTROLS = {"opposite": ([B(0), B(8)], (1.0, -1.0)),
               "four_ccw": ([B(0), B(4), B(8), B(12)], (1.0, -1j, -1.0, 1j))}


def _v2(z, c_spec=C, controls=V2_CONTROLS, freq=F):
    return ring_mode_rows_v2(z_to_s(z, 50.0), port_map=_canonical_map(), freq_hz=freq,
                             c_spec_f=c_spec, feeds=[B(0), B(4)], drives=V2_DRIVES,
                             controls=controls)


def _k_domain_share(blocks, c, feed_idx, v, big):
    jwc = 1j * W * c
    i_s = np.zeros(32, dtype=np.complex128)
    for f, vv in zip(feed_idx, v):
        i_s[f] = vv / big
    i_hat_s = _w() @ i_s
    i_hat = np.zeros(32, dtype=np.complex128)
    for k in range(N_SLOTS):
        idx = [k, N_SLOTS + k]
        y = np.linalg.inv(blocks[k])
        i_hat[idx] = y @ np.linalg.solve(y + jwc * np.eye(2), i_hat_s[idx])
    p = np.abs(i_hat[:N_SLOTS]) ** 2
    return float((p[1] + p[15]) / p.sum())


def test_u7_own_c_mode_wiring():
    c0 = 1.0731e-10
    blocks = _d16_blocks(c_mode=c0)
    z = _block_circulant(blocks)
    v2 = _v2(z)
    own, fixed = v2["own_c_mode"], v2["fixed_c"]
    rel = abs(own["c_f"] / c0 - 1.0)
    print(f"U7 own c_f relative {rel:.3e}; bitwise == fixed c_mode_f: {own['c_f'] == fixed['c_mode_f']}; "
          f"anti {v2['anti_k1']:.15f}")
    assert rel <= TOL
    assert own["c_f"] == fixed["c_mode_f"] == v2["c_own_f"]
    direct = ring_mode_rows(z_to_s(z, 50.0), port_map=_canonical_map(), freq_hz=F,
                            c_f=fixed["c_mode_f"], feeds=[B(0), B(4)], drives=V2_DRIVES)
    assert own == direct
    # fixed_c is v1's call, bit-identical.
    assert fixed == ring_mode_rows(z_to_s(z, 50.0), port_map=_canonical_map(), freq_hz=F, c_f=C,
                                   feeds=[B(0), B(4)], drives=V2_DRIVES)
    for c_spec in (1.1311e-10, 0.5e-10):
        other = _v2(z, c_spec=c_spec)
        assert other["own_c_mode"] == own
        assert other["controls"]["own_c_mode"] == v2["controls"]["own_c_mode"]
    # U5's ideal-source closed form, at c0 (R > 0 keeps the k = 1 solve regular).
    big = 1e15
    _, i_dense = terminated_currents(z, freq_hz=F, c_f=own["c_f"], feed_idx=[0, 4],
                                     v_feeds=(1.0, 1j), z0=big)
    dense = k_spectrum(i_dense[:N_SLOTS])["share_k1"]
    closed = _k_domain_share(blocks, own["c_f"], [0, 4], (1.0, 1j), big)
    print(f"U7 ideal-source at own C: dense {dense:.15f} k-domain {closed:.15f}")
    assert abs(dense - closed) <= TOL


def _c16_mirror_broken():
    blocks = _d16_blocks()
    blocks[3] = blocks[3] + np.array([[0.1j, 0.2], [0.2, -0.05j]])
    blocks[13] = blocks[13] + np.array([[-0.07j, 0.05], [0.11, 0.3j]])   # Ẑ_3 ≠ Ẑ_13
    blocks[1] = blocks[1] + np.array([[0.02j, 0.01], [0.03, 0.0]])        # Ẑ_1 ≠ Ẑ_15
    return blocks


def test_u8_c16_selection_rules():
    blocks = _c16_mirror_broken()
    assert np.max(np.abs(blocks[3] - blocks[13])) > 0.1
    z = _block_circulant(blocks)
    ctl = dict(V2_CONTROLS)
    ctl["four_cw"] = ([B(0), B(4), B(8), B(12)], (1.0, 1j, -1.0, -1j))
    v2 = _v2(z, controls=ctl)
    worst = {}
    for key in ("fixed_c", "own_c_mode"):
        c = v2["controls"][key]

        def fr(lbl, r, c=c):
            return c[lbl]["drives"][lbl]["rings"][r]["fractions"]

        worst[(key, "opposite")] = max(fr("opposite", r)[k] for r in RINGS for k in range(0, 16, 2))
        worst[(key, "four_ccw")] = max(fr("four_ccw", r)[k] for r in RINGS for k in range(16)
                                       if k not in (3, 7, 11, 15))
        worst[(key, "four_cw")] = max(fr("four_cw", r)[k] for r in RINGS for k in range(16)
                                      if k not in (1, 5, 9, 13))
    print("U8 forbidden-bin max: " + "; ".join(f"{k[0]}/{k[1]} {v:.3e}" for k, v in worst.items()))
    assert max(worst.values()) <= TOL
    # U4's open-others network: both controls' share → 1.
    v2o = _v2(_block_circulant(_d16_blocks(open_others=1j * 1e9)))
    shares = {(key, lbl): v2o["controls"][key][lbl]["drives"][lbl]["rings"]["bottom"]["share_k1"]
              for key in ("fixed_c", "own_c_mode") for lbl in V2_CONTROLS}
    print("U8 open-others shares: " + "; ".join(f"{k[0]}/{k[1]} {v:.15f}" for k, v in shares.items()))
    assert min(shares.values()) >= 1.0 - 1e-9


def _control_leak(v2, key):
    return max(v2["controls"][key][lbl]["drives"][lbl]["rings"][r]["leak_k0"]
               for lbl in V2_CONTROLS for r in RINGS)


def test_u9_gauge_negative_control():
    z = _block_circulant(_d16_blocks())
    z_def = z.copy()
    z_def[2, 2] *= 1.5          # single-gap 50 % self-impedance defect, non-feed slot 2 (bottom)
    lc = _control_leak(_v2(z), "own_c_mode")
    bad = _v2(z_def)
    lb = _control_leak(bad, "own_c_mode")
    print(f"U9 control leak clean {lc:.3e} vs defect {lb:.3e} (predicted >= 1e-4); "
          f"defect at SPEC C {_control_leak(bad, 'fixed_c'):.3e}")
    assert lc <= TOL
    assert lb > 1e-9


def test_u10_one_point_f0():
    r1 = 1.5 + 0.2  # _d16_blocks' kk = 1 resistance (R > 0)
    c_s = C
    l_h = 1.05 / (W ** 2 * c_s)
    # λ₁ = R + jωL at 64 MHz: _d16_blocks takes c_mode with x_anti = 1/(W c_mode).
    v2 = _v2(_block_circulant(_d16_blocks(c_mode=1.0 / (W * W * l_h))), c_spec=c_s)
    lam = v2["fixed_c"]["modes"]["1"]["lambda_anti"]
    assert abs(lam[0] - r1) <= 1e-9 and abs(lam[1] / (W * l_h) - 1.0) <= TOL
    f0_closed = 1.0 / (2.0 * np.pi * np.sqrt(l_h * c_s))
    rel = abs(v2["f0_one_point_hz"] / f0_closed - 1.0)
    w0 = 2.0 * np.pi * f0_closed
    rebuilt = ring_mode_rows(z_to_s(_block_circulant(_d16_blocks(c_mode=1.0 / (W * w0 * l_h))), 50.0),
                             port_map=_canonical_map(), freq_hz=f0_closed, c_f=c_s,
                             feeds=[B(0), B(4)], drives={"x": (1, 0)})
    rel_c = abs(rebuilt["c_mode_f"] / c_s - 1.0)
    print(f"U10 f0 {v2['f0_one_point_hz']:.6f} vs 1/(2π√(LC)) {f0_closed:.6f}: relative {rel:.3e}; "
          f"rebuilt at f0 C_mode relative {rel_c:.3e}")
    assert rel <= TOL and rel_c <= TOL
    # Control (closed form): λ₁ = R + j(ωL − 1/(ωC_s')), C_s' = 10 C, true root 5 % below 64 MHz.
    c_ser = 10.0 * c_s
    w_t = 0.95 * W
    l2 = (1.0 / c_s + 1.0 / c_ser) / w_t ** 2
    x1 = W * l2 - 1.0 / (W * c_ser)
    v2c = _v2(_block_circulant(_d16_blocks(c_mode=1.0 / (W * x1))), c_spec=c_s)
    own_closed = F * np.sqrt((1.0 / (W * x1)) / c_s)
    rel_own = abs(v2c["f0_one_point_hz"] / own_closed - 1.0)
    miss = abs(v2c["f0_one_point_hz"] / (0.95 * F) - 1.0)
    print(f"U10 control: f0 {v2c['f0_one_point_hz']:.3f} vs own one-point closed form {own_closed:.3f} "
          f"(relative {rel_own:.3e}); vs true root {0.95 * F:.3f}: miss {miss:.3e} (> 1e-3)")
    assert rel_own <= TOL
    assert miss > 1e-3
