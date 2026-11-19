"""`WF-16` step 13 (§9 item 588, §10 U74): ``n_slots`` on ``ring_mode_rows`` and
``Study.tune`` — the N-slot ring read C10a (`GEO-38`) verified test-locally, now in
``src``, with the 16-slot default proven bitwise against the pre-change source
(the `PORT-24` pinned-source pattern).

W0 (``test_wf16_s13_n_slots_coding``, ungated, no mesh): (N16 / N32) the package
read equals C10a's ``ring_c_modes`` exactly on C10a's synthetics and C3a's
ideal-layout 64×64; (U2) the closed form at ``REUSE_REPRODUCTION_RTOL``; (T32)
``Study.tune(n_slots=32)`` on an injected 64-port; (X) refusals naming
``n_slots``; (c-k2) the k = 2 read misses c₀.

W1 (``test_wf16_s13_default_bitwise_records``, ``FEM_EM_WF16_S13=1``): (B-rm)
``ring_mode_rows`` / ``k_blocks`` / the name route bitwise against ``9b90d7b3`` on
the d1 / d2 copper F-human records, with and without ESR; (B-tune) ``Study.tune``
bitwise against the pre-change ``_study_tune_ring`` body, the recorded
``C_mode`` prints reproduced as strings; (c-bit) the comparator sees an ESR change.

Every ``build`` is stubbed to raise; 0 calls asserted in both nodes.
"""

from __future__ import annotations

import json
import os
import subprocess
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

import tests.validation.test_geo38_f_clinical_matrix as c3
from fem_em_solver.ports import ring_modes as rm
from fem_em_solver.ports.circuit import s_to_z, z_to_s
from fem_em_solver.study import (
    FHumanBirdcageGeometry,
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    QuasiStaticSphereField,
    SphereInBoxGeometry,
    Study,
)
from tests.package.test_wf16_study_tune_f_human_copper import DRIVES, FEEDS, _phantom, _study
from tests.post.test_point_owner_rule import _bitwise_equal_across_ranks
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL
from tests.validation.test_geo38_f_clinical_c_mode import (
    C_MODE_CU_D1_RECORD_F,
    C_MODE_CU_D2_COMPARAND_F,
    H12C_A2_NPZ,
    H16B_D2_NPZ,
    _block_circulant,
    _load,
    _shuffled,
    _sym_blocks,
    map_from_geometry,
    map_from_names_64,
    ring_c_modes,
)
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
from tests.validation.test_wf7_f_human_b1_quadrature import DEGREE, FREQUENCY_HZ
from tests.validation.test_wf10_f_human_copper_esr_root import ESR_OHM

S13_ENV = "FEM_EM_WF16_S13"
REPO = Path(__file__).resolve().parents[2]
#: The last commit to touch ``ring_modes.py`` before this step (`PORT-24` step 1), pinned.
PRE_CHANGE_SHA = "9b90d7b3"
RING_MODES_PATH = "src/fem_em_solver/ports/ring_modes.py"
Z0 = 50.0
C0_F = 1.0731e-10  # C10a's (K1) synthetic c₀
OMEGA = 2.0 * np.pi * FREQUENCY_HZ
#: The recorded prints (20261015T213534Z_WF-16-s12-W1.log:7375; C10a W1 :73 / :83).
RECORD_PRINTS = {"d1": "1.181185799070e-10", "d2": "1.145221310828e-10"}
#: Predictions (printed, never asserted — rule (e)).
PRED_U2 = 2e-15
PRED_C_K2_MISS = 0.554


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _knob() -> bool:
    raw = os.environ.get(S13_ENV)
    if raw is None or raw.strip() == "":
        return False
    if raw.strip() == "1":
        return True
    raise ValueError(f"{S13_ENV}={raw!r}: only unset/empty or '1' are accepted")


def _canon(x):
    return json.dumps(x, sort_keys=True, allow_nan=True)


def _stub_builds(monkeypatch):
    calls = {"n": 0}

    def stub(self, comm=MPI.COMM_WORLD):
        calls["n"] += 1
        raise AssertionError("build must not be called (WF-16 step 13: no mesh)")

    for cls in (SphereInBoxGeometry, FHumanBirdcageGeometry, FourPortBirdcageGeometry):
        monkeypatch.setattr(cls, "build", stub)
    return calls


def _worst_vs_reader(local, ref, n):
    """Worst relative difference, package ``ring_mode_rows`` vs C10a's reader,
    over λ_anti k = 0…N/2, ``c_mode_f``, ``c_k2_f``, ``off_block``."""
    worst = 0.0
    for k in range(n // 2 + 1):
        a = local["modes"][k]["lambda_anti"]
        b = complex(*ref["modes"][str(k)]["lambda_anti"])
        worst = max(worst, abs(a - b) / abs(b))
    for key in ("c_mode_f", "c_k2_f", "off_block"):
        worst = max(worst, abs(local[key] - ref[key]) / abs(ref[key]))
    return worst


def _pkg_rows(s, pm, n):
    return rm.ring_mode_rows(
        s,
        port_map=pm,
        n_slots=n,
        freq_hz=FREQUENCY_HZ,
        c_f=C0_F,
        feeds=[(0, "bottom"), (n // 4, "bottom")],
        drives={"x": (1.0, 0.0)},
        z0=Z0,
    )


def _broken(n, seed_noise, seed_shuffle):
    rng = np.random.default_rng(seed_noise)
    z = _block_circulant(_sym_blocks(n, OMEGA, C0_F), n)
    a = rng.normal(size=(2 * n, 2 * n)) + 1j * rng.normal(size=(2 * n, 2 * n))
    zs, pm = _shuffled(z + 0.5 * (a + a.T), n, seed=seed_shuffle)
    return z_to_s(zs, Z0), pm


# --------------------------------------------------------------------------- W0
def test_wf16_s13_n_slots_coding(monkeypatch):
    """(N16 / N32), (U2), (T32), (X), (c-k2); no mesh."""
    _knob()  # a bad knob value raises here too
    comm = MPI.COMM_WORLD
    calls = _stub_builds(monkeypatch)

    # (N16 / N32) on (i) the closed-form synthetics, (ii) C16-broken, (iii) C3a's 64x64.
    cases = {}
    syn = {}
    for n in (16, 32):
        zs, pm = _shuffled(_block_circulant(_sym_blocks(n, OMEGA, C0_F), n), n, seed=581 + n)
        syn[n] = (z_to_s(zs, Z0), pm)
        cases[f"(i) synthetic N={n}"] = (syn[n][0], pm, n)
    cases["(ii) C16-broken N=16"] = (*_broken(16, 5810, 17), 16)
    cases["(ii) C32-broken N=32"] = (*_broken(32, 5820, 33), 32)
    az, zc = c3._ideal_layout()
    geo = map_from_geometry(az, zc)
    s64 = c3._synthetic_s(az, zc)
    cases["(iii) C3a ideal 64x64 N=32"] = (s64, geo, 32)
    worst = {}
    nmodes = {}
    for label, (s, pm, n) in cases.items():
        loc = ring_c_modes(s, pm, freq_hz=FREQUENCY_HZ, n_slots=n, z0=Z0)
        ref = _pkg_rows(s, pm, n)
        worst[label] = _worst_vs_reader(loc, ref, n)
        nmodes[label] = (len(ref["modes"]), n // 2 + 1)
        _p0(
            f"\n[WF-16 s13 N] {label}: worst rel vs C10a ring_c_modes {worst[label]:.3e} "
            f"(asserted == 0; predicted 0); modes {nmodes[label][0]} (N/2+1 = {nmodes[label][1]}); "
            f"off-block {ref['off_block']:.3e}"
        )

    # (U2) the closed form on (i).
    u2 = {}
    for n in (16, 32):
        r = _pkg_rows(*syn[n], n)
        u2[n] = (
            abs(r["c_mode_f"] / C0_F - 1.0),
            abs(r["c_k2_f"] * OMEGA * 52.0 - 1.0),
            r["modes"]["1"]["anti"],
        )
        _p0(
            f"[WF-16 s13 U2] N={n}: |C_mode/c0-1| {u2[n][0]:.3e}, |C2 w 52 - 1| {u2[n][1]:.3e} "
            f"(band REUSE_REPRODUCTION_RTOL {REUSE_REPRODUCTION_RTOL:g}; beside "
            f"EXACT_IDENTITY_RTOL {EXACT_IDENTITY_RTOL:g}; C10a read <= 1.6e-15; predicted "
            f"<= {PRED_U2:g}); anti(k=1) {u2[n][2]:.15f}"
        )

    # (T32) Study.tune(n_slots=32) on the injected shuffled N = 32 synthetic.
    g32 = FHumanBirdcageGeometry(leg_count=32, as_hole=True)
    ids = [f"P{o}" for o in g32.ring_ordinals()]
    assert len(ids) == 64 and ids[0] == "P33" and ids[-1] == "P96"
    z32 = _block_circulant(_sym_blocks(32, OMEGA, C0_F), 32)
    perm = np.random.default_rng(581 + 32).permutation(64)
    s32 = z_to_s(z32[np.ix_(perm, perm)], Z0)
    ids_sh = [ids[p] for p in perm]
    assert map_from_names_64(ids_sh) == syn[32][1]
    assert rm.ring_port_map_from_names(ids_sh, n_slots=32) == map_from_names_64(ids_sh)
    st32 = _study(geom=g32)
    st32.result = SimpleNamespace(
        s_matrix=s32, port_definitions=[SimpleNamespace(port_id=n) for n in ids_sh]
    )
    t1, t2 = st32.tune(k=1, n_slots=32), st32.tune(k=2, n_slots=32)
    t1d = st32.tune(k=1)
    rd = ring_c_modes(
        s32,
        map_from_names_64(ids_sh),
        freq_hz=FREQUENCY_HZ,
        n_slots=32,
        z0=float(st32.excitation.z0_ohm),
    )
    t_e1 = abs(t1.c_mode_farad / C0_F - 1.0)
    t_e2 = abs(t2.c_mode_farad * OMEGA * 52.0 - 1.0)
    bit = _bitwise_equal_across_ranks(comm, np.array([t1.c_mode_farad, t2.c_mode_farad]))
    _p0(
        f"[WF-16 s13 T32] tune(k=1, n_slots=32) {t1.c_mode_farad:.12e} F (rel vs c0 {t_e1:.3e}), "
        f"tune(k=2) {t2.c_mode_farad:.12e} F (|C2 w 52 - 1| {t_e2:.3e}); == reader: "
        f"{t1.c_mode_farad == rd['c_mode_f']} / {t2.c_mode_farad == rd['c_k2_f']}; len(z_k) "
        f"{len(t1.z_k)}; default == explicit {t1d.c_mode_farad == t1.c_mode_farad}; "
        f"bitwise across {comm.size} ranks {bit}"
    )

    # (c-k2) asserted, backed: the k = 2 read misses c0 by far more than REDUCTION_BAND.
    miss = abs(t2.c_mode_farad / C0_F - 1.0)
    _p0(
        f"[WF-16 s13 c-k2] |C(k=2)/c0 - 1| {miss:.4f} (> REDUCTION_BAND {REDUCTION_BAND:g}; "
        f"predicted {PRED_C_K2_MISS})"
    )

    # Printed only: the top ring offset by one slot at N = 32.
    pm_off = [(j if r == "bottom" else (j + 1) % 32, r) for j, r in syn[32][1]]
    r_off = _pkg_rows(syn[32][0], pm_off, 32)
    _p0(
        f"[WF-16 s13 print] top ring offset one slot, N=32: C_mode rel "
        f"{abs(r_off['c_mode_f'] / C0_F - 1.0):.3e}, anti(k=1) {r_off['modes']['1']['anti']:.6f} "
        f"(cos^2(pi/32) = {np.cos(np.pi / 32) ** 2:.6f})"
    )

    # (X) refusals naming n_slots, no build.
    for bad in (2, 3, 15, True, 16.0):
        with pytest.raises((ValueError, TypeError), match="n_slots") as ei:
            (
                _pkg_rows(*syn[16], bad)
                if not isinstance(bad, (bool, float))
                else rm.ring_mode_rows(
                    syn[16][0],
                    port_map=syn[16][1],
                    n_slots=bad,
                    freq_hz=FREQUENCY_HZ,
                    c_f=C0_F,
                    feeds=[(0, "bottom")],
                    drives={"x": (1.0,)},
                    z0=Z0,
                )
            )
        expect = TypeError if isinstance(bad, (bool, float)) else ValueError
        assert ei.type is expect, (bad, ei.type)
    with pytest.raises(ValueError, match="n_slots"):
        _pkg_rows(s64, geo, 16)
    with pytest.raises(ValueError, match="n_slots"):
        rm.ring_port_map_from_names(["P17"], n_slots=32)
    fresh32 = _study(geom=FHumanBirdcageGeometry(leg_count=32, as_hole=True))
    with pytest.raises(ValueError, match="n_slots"):
        fresh32.tune(k=1, n_slots=16)
    assert fresh32.result is None
    held = st32.result
    with pytest.raises(ValueError, match="n_slots"):
        st32.tune(k=1, n_slots=16)
    assert st32.result is held
    four = Study(
        FourPortBirdcageGeometry.f_small(),
        _phantom(),
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_HZ,
        degree=DEGREE,
    )
    sphere = Study(
        SphereInBoxGeometry(0.05, 0.15, 0.01, 0.03),
        _phantom(),
        QuasiStaticSphereField(),
        frequency_hz=FREQUENCY_HZ,
        degree=DEGREE,
    )
    for other in (four, sphere):
        with pytest.raises(TypeError, match="n_slots"):
            other.tune(k=1, n_slots=16)
        assert other.result is None
    _p0(
        f"[WF-16 s13 X] refusals: n_slots in (2, 3, 15, True, 16.0); 64-port at 16; P17 at 32; "
        f"32-leg Study at 16 (result untouched); four-port / sphere TypeError; builds "
        f"{calls['n']}"
    )

    for label in cases:
        assert worst[label] == 0.0, (label, worst[label])
        assert nmodes[label][0] == nmodes[label][1], (label, nmodes[label])
    for n in (16, 32):
        assert u2[n][0] <= REUSE_REPRODUCTION_RTOL and u2[n][1] <= REUSE_REPRODUCTION_RTOL, u2[n]
        assert u2[n][2] >= 1.0 - REUSE_REPRODUCTION_RTOL, u2[n]
    assert t_e1 <= REUSE_REPRODUCTION_RTOL and t_e2 <= REUSE_REPRODUCTION_RTOL, (t_e1, t_e2)
    assert t1.c_mode_farad == rd["c_mode_f"] and t2.c_mode_farad == rd["c_k2_f"]
    assert len(t1.z_k) == 17 and len(t2.z_k) == 17
    assert t1d.c_mode_farad == t1.c_mode_farad and np.array_equal(t1d.z_k, t1.z_k)
    assert bit
    assert miss > REDUCTION_BAND, miss
    assert calls["n"] == 0


# --------------------------------------------------------------------------- W1
@pytest.fixture(scope="module")
def pre():
    """The pre-change ``ring_modes`` module, from the pinned commit."""
    src = subprocess.run(
        [
            "git",
            "-c",
            f"safe.directory={REPO}",
            "-C",
            str(REPO),
            "show",
            f"{PRE_CHANGE_SHA}:{RING_MODES_PATH}",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    assert "n_slots" not in src, "pinned pre-change source already carries n_slots"
    assert "esr_ohm" in src, "pinned pre-change source predates PORT-24 step 1"
    mod = types.ModuleType("fem_em_solver.ports._ring_modes_pre_wf16s13")
    mod.__package__ = "fem_em_solver.ports"
    exec(compile(src, f"{PRE_CHANGE_SHA}:{RING_MODES_PATH}", "exec"), mod.__dict__)
    return mod


def _pre_tune_ring(pre, st, k):
    """The pre-change ``_study_tune_ring`` body (study.py at 9b90d7b3's era,
    :1140–1150), with ``pre``'s functions and the Study's own ``z0`` / ``omega``."""
    names = [pd.port_id for pd in st.result.port_definitions]
    perm = pre._canonical_permutation(pre.ring_port_map_from_names(names))
    s = np.asarray(st.result.s_matrix, dtype=np.complex128)
    z = s_to_z(s, float(st.excitation.z0_ohm))
    blocks, off = pre.k_blocks(z[np.ix_(perm, perm)])
    lam = np.array([pre.antisymmetric_eigenpair(blocks[m])[0] for m in range(9)])
    return float(1.0 / (st.omega * lam[k].imag)), lam, float(off)


def test_wf16_s13_default_bitwise_records(monkeypatch, pre):
    """(B-rm), (B-tune), (c-bit) on the two copper F-human records."""
    if not _knob():
        pytest.skip(f"{S13_ENV} unset: the records node reads the gitignored npz")
    calls = _stub_builds(monkeypatch)
    recs = {}
    for tag, path in (("d1", H12C_A2_NPZ), ("d2", H16B_D2_NPZ)):
        if not path.exists():
            pytest.fail(f"{path} absent: stop (BLOCKED, rule (d)); never re-solve")
        d = _load(path)
        recs[tag] = (np.array(d["s32"], dtype=np.complex128), d["order"])

    rm_ok, kb_ok, map_ok, cbit = [], [], [], []
    for tag, (s, order) in recs.items():
        kw = dict(
            port_map=pre.ring_port_map_from_names(order),
            freq_hz=FREQUENCY_HZ,
            feeds=FEEDS,
            drives=DRIVES,
        )
        map_ok.append(rm.ring_port_map_from_names(order) == pre.ring_port_map_from_names(order))
        for c in (C_MODE_CU_D1_RECORD_F, C_MODE_CU_D2_COMPARAND_F):
            for r in (None, ESR_OHM):
                extra = {} if r is None else {"esr_ohm": r}
                old = _canon(pre.ring_mode_rows(s, c_f=c, **kw, **extra))
                new = _canon(rm.ring_mode_rows(s, c_f=c, **kw, **extra))
                new16 = _canon(rm.ring_mode_rows(s, c_f=c, n_slots=16, **kw, **extra))
                rm_ok.append(old == new == new16)
                if not rm_ok[-1]:
                    _p0(f"[WF-16 s13 B-rm] {tag} c={c!r} r={r}: DIFFER")
            # (c-bit): the comparator sees an ESR change.
            cbit.append(
                _canon(rm.ring_mode_rows(s, c_f=c, esr_ohm=ESR_OHM, **kw))
                != _canon(pre.ring_mode_rows(s, c_f=c, **kw))
            )
        perm = pre._canonical_permutation(pre.ring_port_map_from_names(order))
        zc = s_to_z(s, Z0)[np.ix_(perm, perm)]
        bn, on = rm.k_blocks(zc)
        bo, oo = pre.k_blocks(zc)
        kb_ok.append(all(np.array_equal(bn[k], bo[k]) for k in range(16)) and on == oo)
    _p0(
        f"\n[WF-16 s13 B-rm] ring_mode_rows bitwise vs {PRE_CHANGE_SHA} (default and n_slots=16, "
        f"r in (0, ESR_OHM), two C, two records): {sum(rm_ok)}/{len(rm_ok)}; name route "
        f"{map_ok}; k_blocks + off-block {kb_ok}; (c-bit) ESR-vs-lossless differs {cbit}"
    )

    tune_ok, prints, rel = {}, {}, {}
    for tag, (s, order) in recs.items():
        st = _study()
        st.result = SimpleNamespace(
            s_matrix=s, port_definitions=[SimpleNamespace(port_id=n) for n in order]
        )
        oks = []
        for k in (1, 2):
            t = st.tune(k=k)
            c_old, lam_old, off_old = _pre_tune_ring(pre, st, k)
            oks.append(
                t.c_mode_farad == c_old
                and np.array_equal(t.z_k, lam_old)
                and t.off_diagonal_residual == off_old
            )
            if k == 1:
                prints[tag] = f"{t.c_mode_farad:.12e}"
                const = C_MODE_CU_D1_RECORD_F if tag == "d1" else C_MODE_CU_D2_COMPARAND_F
                rel[tag] = abs(t.c_mode_farad / const - 1.0)
        tune_ok[tag] = oks
        _p0(
            f"[WF-16 s13 B-tune] {tag}: tune(k=1,2) == pre-change body {oks}; C_mode "
            f"{prints[tag]} F (record print {RECORD_PRINTS[tag]}; rel vs constant {rel[tag]:.3e}, "
            f"band REUSE_REPRODUCTION_RTOL {REUSE_REPRODUCTION_RTOL:g}; C10a read 3.515e-13 / "
            f"4.212e-13)"
        )

    assert all(rm_ok) and len(rm_ok) == 8  # 2 records x 2 C x 2 r, each old == new == new16
    assert all(map_ok) and all(kb_ok)
    assert all(cbit)
    for tag in recs:
        assert all(tune_ok[tag]), (tag, tune_ok[tag])
        assert prints[tag] == RECORD_PRINTS[tag], (tag, prints[tag])
        assert rel[tag] <= REUSE_REPRODUCTION_RTOL, (tag, rel[tag])
    assert calls["n"] == 0
