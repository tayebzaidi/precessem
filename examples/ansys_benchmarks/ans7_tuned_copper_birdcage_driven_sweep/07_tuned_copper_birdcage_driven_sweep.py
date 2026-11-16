"""`ANS-7` v2 runnable half (§9 item 63, §10 T10): the tuned copper F-small
birdcage, C on every gap, the drive across C₁ -- Route B.

SPEC v2 (``SPEC.md`` in this directory, operator directive 2026-09-28):
C = 66.29 pF on all four leg gaps, gap 1 driven (C₁ across the port),
parallel resonance read on ``Y_in`` over the FineWindow 62.0 → 67.0 MHz in
0.1 MHz steps; degree 1; Cu and PEC columns.

Route B (§9 item 69): `TH-14`'s ``_build_ladder`` four-port on the F-small hole
mesh (imported, rule (a) -- never re-implemented), κ-corrected sheet widths,
then the circuit layer: ``Z_bare = tuned_input(s4, f, C)[0]`` (`PORT-15`'s
reduction, ports 2–4 in C) and ``Z_in = 1/(1/Z_bare + jωC)``.  Route A
(in-model), the reduction band between routes and the private writer are
stage ``A`` (§9 item 68).

κ is derived once, by `PORT-22`'s ``_kappa_at`` route (`PORT-14` step 3's
pooled ``terminal_form_deficit``) on an **uncorrected Cu P1 solve at 64 MHz**
of this fixture, and the same κ goes on all four sheets, both columns.

Stages (``FEM_EM_ANS7_STAGE``), each saving to the gitignored
``/workspace/logs/ans7/`` (``FEM_EM_ANS7_RESULTS_DIR`` overrides):

* ``B1`` -- anchor 1 (the uncorrected 64 MHz ladder reproduces the tracked
  ``ans:6`` ``metrics.json`` within ``LEG_D0_REPRODUCTION_BAND``), the
  asserted Cu-vs-PEC negative control on that 4×4, κ, then the corrected
  sweep 62.0 → 64.4 MHz (25 points).  Anchor 1 red ⇒ stop before the sweep.
* ``B2`` -- reads κ from ``B1``, sweeps 64.5 → 67.0 MHz (26 points), then
  assembles the 51-point window: anchors (4)–(5), rows (a)–(d), the
  predicted readings and the 2·C control; writes ``metrics.json`` and
  ``COMPARISON.md`` beside this script (AED and Route A cells blank).
* ``A`` -- Route A (in-model C on every sheet) over the same 51 points: κ
  recomputed in-run within ``KAPPA_RECOMPUTE_ATOL`` of ``B1``'s, every
  ``S11_bare`` A − B residual (102, both columns) ≤ ``REDUCTION_BAND``, and
  exactly one − → + zero of ``Im Y_in`` per column; controls (i) (circuit at
  ``CONTROL_I_C_FACTOR``·C vs in-model at C, 64 MHz) and (ii) (uncorrected
  Route A vs the uncorrected circuit) are predicted and printed.  Writes the
  Route A cells into ``metrics.json``, a combined Cu/PEC XDMF to the results
  directory and, unless ``FEM_EM_ANS_SKIP_PRIVATE=1``, the private table.

Asserted at every point, both columns (bands imported, never restated):
(2) reciprocity ≤ ``RECIPROCITY_BAND``, copper identity ≤
``DISCRETE_IDENTITY_RTOL``.  At assembly: (4) exactly one − → + zero of
``Im Y_in`` per column in the window; (5) k = ±1 gap-current share at the
point nearest ``f₀`` > ``K_SHARE_FLOOR``.

*Record, not a gate* (2026-09-29 03:00 review, rule (h)): (3) the S route
equals the SPEC's Z formula (``terminate_z`` on the network Z of the same S)
to ≤ ``S_VS_Z_ROUTE_RTOL`` -- circuit-layer algebra on one S, true by
construction; the script still stops on a miss as a code-consistency check.

*Predicted, printed only:* PEC ``f₀`` in [63.9, 64.4] MHz; k-share ≈ 1; κ
beside the gate mesh's 1.060762e-02; the same 4×4s reduced at 2·C show no
− → + zero of ``Im Y_in`` in the window (k = 1 near 45.3 MHz).  The v1
prediction Q ≈ 68 is withdrawn (mis-imported from v1; 03:00 review); the
v2 readings are Q 918 (Cu) / 1536 (PEC), printed.

Degree-2 exact-rectangle stages (§9 item 230, §10 T20, additive):
``X1``–``X4`` solve the seven ``DEFAULT_NODES_MHZ`` Route B nodes at degree 2 on
the exact-rectangle hole mesh, κ off (two per window; ``X4`` also assembles,
asserts (R2x) against ``ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2``, (N-deg), (K), and
adds the ``degree2_exact`` / ``superseded`` keys to ``metrics.json`` and
regenerates ``COMPARISON.md``); ``XA`` drives Route A at degree 2 (Cu) at
62.6 / 62.7 / 62.8 MHz against the assembled Route B, (AB).  Stage files go to
the gitignored ``/workspace/logs/ans7_p2/``.  Every gate record and band is
imported lazily from the test modules, never restated.  The degree-1 ragged
κ-corrected keys of ``metrics.json`` are superseded and kept in place for their
readers.  Run at ``-n 8``; B1 / B2 / A are not re-run.

Tuned-to-Larmor stages (§9 item 288, §10 T23, additive): ``XT64`` / ``XT128`` re-drive
Route A at degree 2 on the ×1 exact-sheet mesh at C_in (read from the stored
``tune_p2.json`` of the field, never recomputed) at the field's three tune points
(63.9 / 64.0 / 64.1 and 127.64 / 127.74 / 127.84 MHz).  They assert (AB′) / (AB″), (T2′) /
(T2″), the reproduction (R-T) of the registered run, and the controls (N1′) / (N0′) and
(N-C); (N2′) is printed.  They write ``degree2_exact["tuned_larmor"]`` into
``metrics.json`` (with (K)) and a third ``COMPARISON.md`` section.  The field set is
chosen by the stage name alone; every ``FEM_EM_ANS7_EXACT_*`` knob stays unset.

Capacitor ESR stage (§9 item 378, §10 T27, additive): ``XT64R`` re-drives Route A at degree 2
on the ×1 exact-sheet mesh at (r = 0.2 Ω, C_in(r)) (read from the stored ``tune_p2_esr0.2.json``
of item 311, never recomputed) at 63.9 / 64.0 / 64.1 MHz, the series-RC sheet on P2–P4 and the
feed gap's (r + C) branch in parallel at P1 (the circuit layer's new ``esr_ohm`` keyword; 0 keeps
the old expression by a branch, asserted bitwise as (R0)).  It asserts (AB-r), (T2-r) and the
reproduction (R-T-r) of item 311's run; prints (N-L), (N-C-r) (asserted only if item 311's own
separation exceeded 1e-7), Q(r)/Q_0 against the closed form and C_in(r)/C_in(0) − 1; writes
``degree2_exact["tuned_larmor_esr"]["64.0"]`` into ``metrics.json`` (with (K)) and a fourth
``COMPARISON.md`` section.  The lossless ``tuned_larmor`` rows stay, labelled lossless.  Selected
by the stage name alone; every ``FEM_EM_ANS7_EXACT_*`` knob stays unset.

``XT128R`` (§9 item 391, §10 T36) is its 127.74 MHz sibling at 127.64 / 127.74 / 127.84 MHz from
item 387's stored ``logs/ans7_exact/t127.74/`` stage files: (AB-r'), (T2-r') in ``T22P_ZERO_RESOLUTION``
bands, (R0'), (R-T-r'), (N-C-S') (asserted, the stored S11 separation exceeds 1e-7 at every point),
(K) over ``tuned_larmor`` and ``tuned_larmor_esr["64.0"]`` too; (N-C-r') and (N-L') printed
(item 387's reading did not meet their bars).  It writes ``tuned_larmor_esr["127.74"]`` and a
sibling 127.74 MHz ESR section appended to ``COMPARISON.md`` (item 392: its own writer).

Not claimed: any absolute ``Z_in``, a match.  No AED value is read or written.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[3]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from dolfinx import default_scalar_type  # noqa: E402

from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_port_circuit_layer_field import (  # noqa: E402
    TUNING_TERMINATED_INDICES,
    tuned_input,
)
from tests.validation.test_port_package_sparameters import (  # noqa: E402
    REFERENCE_IMPEDANCE_OHM,
)

CASE_DIR = Path(__file__).resolve().parent
FIGURE_DIR = CASE_DIR / "figures"
FIGURE_BASENAME = "ans7_tuned_copper_birdcage_driven_sweep"
ANS6_PATH = (CASE_DIR.parent / "ans6_copper_birdcage_four_port_10_64_128MHz"
             / "06_copper_birdcage_four_port_10_64_128MHz.py")

# SPEC v2: one capacitor value on every gap (``SPEC.md`` "Every gap carries the
# same capacitor, C = 66.29 pF exactly").
C_F = 6.629e-11
STAGE_ENV = "FEM_EM_ANS7_STAGE"
RESULTS_DIR_ENV = "FEM_EM_ANS7_RESULTS_DIR"
DEFAULT_RESULTS_DIR = "/workspace/logs/ans7"
STAGE_B1_MHZ = tuple(round(62.0 + 0.1 * k, 6) for k in range(25))   # 62.0 .. 64.4
STAGE_B2_MHZ = tuple(round(64.5 + 0.1 * k, 6) for k in range(26))   # 64.5 .. 67.0
ANCHOR1_FREQUENCY_HZ = 64.0e6
# Anchor (3): the S route vs the SPEC's Z formula -- an identity of two
# algebraically equal reductions; machine precision with headroom.
S_VS_Z_ROUTE_RTOL = 1.0e-9
# Anchor (5): the k = ±1 share must dominate at f0.
K_SHARE_FLOOR = 0.5
# The asserted negative control's floor (item 50's, same comparison).
NEGATIVE_CONTROL_FLOOR = 1.0e-6
# *Predicted, printed only* (§9 item 63).
PREDICTED_PEC_F0_MHZ = (63.9, 64.4)
PREDICTED_Q = 68.0
GATE_MESH_KAPPA = 1.060762e-02
CONTROL_C_FACTOR = 2.0


def _load_example(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _cplx(z):
    return {"re": float(np.real(z)), "im": float(np.imag(z))}


def _mat(m):
    return [[_cplx(v) for v in row] for row in np.asarray(m)]


def _unmat(rows):
    return np.array([[complex(v["re"], v["im"]) for v in row] for row in rows],
                    dtype=np.complex128)


def _results_dir() -> Path:
    return Path(os.environ.get(RESULTS_DIR_ENV, "") or DEFAULT_RESULTS_DIR)


# -- the SPEC's circuit layer (definitions per SPEC.md; never the AED script) --

def _gap_z(f_hz, c_f=C_F, esr_ohm=0.0):
    """One ring-gap branch (§9 item 378, §10 T27; additive).

    ``esr_ohm == 0`` keeps the old expression ``1/(jωC)`` by a branch, not arithmetic
    (``1/(0 + 1/(jωC))`` is not bitwise ``jωC``; item 308's trap).  ``esr_ohm != 0`` is
    `PORT-24`'s series-RC branch, imported, never ``r + 1/(jωC)`` restated.  Both
    `_route_a_solve` and `route_b_z` take the gap impedance from here.
    """
    if esr_ohm == 0.0:
        return 1.0 / (1j * (2.0 * np.pi * f_hz) * c_f)
    from fem_em_solver.ports.lumped import series_rlc_impedance
    return series_rlc_impedance(f_hz, r_ohm=esr_ohm, l_h=0.0, c_f=c_f)


def z_cap(omega, c_f=C_F, esr_ohm=0.0):
    if esr_ohm == 0.0:
        return 1.0 / (1j * omega * c_f)
    return _gap_z(omega / (2.0 * np.pi), c_f, esr_ohm)


def add_c1(z_bare, omega, c_f=C_F, esr_ohm=0.0):
    if esr_ohm == 0.0:
        return 1.0 / (1.0 / z_bare + 1j * omega * c_f)
    # the feed gap's (r + C) branch in parallel (ESR ruling 311-Q5)
    return 1.0 / (1.0 / z_bare + 1.0 / z_cap(omega, c_f, esr_ohm))


def route_b_s(s4, f_hz, c_f=C_F):
    """Route B, S route: `PORT-15`'s ``tuned_input`` (ports 2–4 in C), then C₁."""
    z_bare = complex(tuned_input(np.asarray(s4), f_hz, c_f)[0])
    return add_c1(z_bare, 2.0 * np.pi * f_hz, c_f), z_bare


def route_b_z(z4, f_hz, c_f=C_F, esr_ohm=0.0):
    """Route B, Z route: (Z_tt + Z_C) I_t = −Z_t1 I₁; returns Z_in, Z_bare, gap currents.

    ``esr_ohm`` (§9 item 378, additive; 0 keeps the old path bitwise): r in series with C
    on every gap, the feed's branch included.
    """
    omega = 2.0 * np.pi * f_hz
    z4 = np.asarray(z4, dtype=np.complex128)
    z_tt = z4[1:, 1:] + _gap_z(f_hz, c_f, esr_ohm) * np.eye(3)
    i_t = -np.linalg.solve(z_tt, z4[1:, 0])
    z_bare = z4[0, 0] + z4[0, 1:] @ i_t
    return add_c1(z_bare, omega, c_f, esr_ohm), z_bare, np.concatenate(([1.0 + 0j], i_t))


def k1_share(gap_currents):
    g = np.fft.fft(np.asarray(gap_currents))
    p = np.abs(g) ** 2
    return float((p[1] + p[3]) / p.sum())


def parallel_rows(freq_hz, z_in):
    """Rows (a)–(c) on Y_in = 1/Z_in: f0 at every − → + zero of Im Y_in."""
    freq_hz = np.asarray(freq_hz, dtype=float)
    y = 1.0 / np.asarray(z_in)
    b = y.imag
    up = [k for k in range(len(b) - 1) if b[k] < 0 <= b[k + 1]]
    down = [k for k in range(len(b) - 1) if b[k] >= 0 > b[k + 1]]

    def at(k):
        t = -b[k] / (b[k + 1] - b[k])
        return t, freq_hz[k] + t * (freq_hz[k + 1] - freq_hz[k])

    rows = {"crossings_up_hz": [float(at(k)[1]) for k in up],
            "crossings_down_hz": [float(at(k)[1]) for k in down]}
    if not up:
        return rows
    k = up[0]
    t, f0 = at(k)
    g0 = y.real[k] + t * (y.real[k + 1] - y.real[k])
    dbdf = (b[k + 1] - b[k]) / (freq_hz[k + 1] - freq_hz[k])
    rows.update(f0_hz=float(f0), r_in_ohm=float(1.0 / g0), dbdf_s_per_hz=float(dbdf),
                q_loaded=float(f0 * dbdf / (2.0 * g0)),
                bracket_hz=[float(freq_hz[k]), float(freq_hz[k + 1])])
    return rows


def z_from_s(s4, z0=None):
    """The 50 Ω network Z of the stored S: Z = z0 (I + S)(I − S)⁻¹.

    W1 (``20260929T021845Z_ANS-7-W1-B1.log:3624–3651``) fed the sweep's own
    ``z_matrix`` to the Z route and missed anchor (3) by factors 6–307: that
    matrix is not the network Z of the returned S at z0 (both columns, every
    point).  The SPEC's Z formula is applied to the S-derived network Z, the
    matrix Ansys Circuit / skrf build from the same S; the sweep's
    ``z_matrix`` is kept in the stage file as ``z_sweep`` for the record.
    """
    z0 = float(REFERENCE_IMPEDANCE_OHM) if z0 is None else float(z0)
    s4 = np.asarray(s4, dtype=np.complex128)
    eye = np.eye(s4.shape[0])
    return z0 * (eye + s4) @ np.linalg.inv(eye - s4)


def _circuit_point(f_hz, s4):
    """Route B at one frequency from the stored 4×4 S: both routes, gap currents."""
    z4 = z_from_s(s4)
    z_in_s, z_bare = route_b_s(s4, f_hz)
    z_in_z, _zb, gaps = route_b_z(z4, f_hz)
    dev = abs(z_in_s - z_in_z) / abs(z_in_z)
    return {"s": _mat(s4), "S11_bare": _cplx(
        (z_bare - REFERENCE_IMPEDANCE_OHM) / (z_bare + REFERENCE_IMPEDANCE_OHM)),
        "Z_bare": _cplx(z_bare), "Z_in": _cplx(z_in_s), "Z_in_zroute": _cplx(z_in_z),
        "s_vs_z_route": float(dev), "gap_currents": [_cplx(g) for g in gaps]}


def _recompute_b1(b1):
    """B1's circuit layer re-derived from its stored S (anchor (3) re-checked)."""
    failures = []
    for p in b1["points"]:
        for name in ("cu", "pec"):
            recip = p[name]["reciprocity"]
            p[name] = {**_circuit_point(p["f_hz"], _unmat(p[name]["s"])), "reciprocity": recip}
            if not p[name]["s_vs_z_route"] <= S_VS_Z_ROUTE_RTOL:
                failures.append(("s-vs-z-route-B1", p["f_hz"], name, p[name]["s_vs_z_route"]))
    if MPI.COMM_WORLD.rank == 0:
        worst = max(p[c]["s_vs_z_route"] for p in b1["points"] for c in ("cu", "pec"))
        print(f"[ans7] B1 re-derived from stored S: worst S-vs-Z route {worst:.3e} "
              f"(RECORD, rule (h); stops on > {S_VS_Z_ROUTE_RTOL:.0e} as a "
              f"code-consistency check)", flush=True)
    return failures


# -- stages --------------------------------------------------------------------

def _point_records(ladder, freqs_mhz, ans6):
    """Per point, both columns: 4×4 S / Z, reciprocity, identity; anchors (2)–(3)."""
    pts, failures = [], []
    for mhz in freqs_mhz:
        key = f"{float(mhz):g}"
        rec = ladder["freqs"][key]
        f_hz = float(rec["f_hz"])
        cols = {"cu": rec["sigma"][ans6.COPPER_SIGMA], "pec": rec["pec"]}
        pt = {"f_hz": f_hz, "identity_residual": float(rec["identity_residual"])}
        if not rec["identity_residual"] <= ans6.DISCRETE_IDENTITY_RTOL:
            failures.append(("identity", key, rec["identity_residual"]))
        for name, sw in cols.items():
            s4 = np.asarray(sw["s"], dtype=np.complex128)
            pt[name] = _circuit_point(f_hz, s4)
            pt[name]["z_sweep"] = _mat(np.asarray(sw["z"], dtype=np.complex128))
            pt[name]["reciprocity"] = float(sw["reciprocity"])
            if not sw["reciprocity"] <= ans6.RECIPROCITY_BAND:
                failures.append(("reciprocity", key, name, sw["reciprocity"]))
            if not pt[name]["s_vs_z_route"] <= S_VS_Z_ROUTE_RTOL:
                failures.append(("s-vs-z-route", key, name, pt[name]["s_vs_z_route"]))
        pts.append(pt)
        if MPI.COMM_WORLD.rank == 0:
            print(f"[ans7] f = {mhz:7.3f} MHz: Z_in Cu {pt['cu']['Z_in']['re']:+.9e} "
                  f"{pt['cu']['Z_in']['im']:+.9e}j  PEC {pt['pec']['Z_in']['re']:+.9e} "
                  f"{pt['pec']['Z_in']['im']:+.9e}j ohm; recip Cu {pt['cu']['reciprocity']:.2e} "
                  f"PEC {pt['pec']['reciprocity']:.2e}; identity {pt['identity_residual']:.2e}; "
                  f"S-vs-Z route Cu {pt['cu']['s_vs_z_route']:.2e} PEC "
                  f"{pt['pec']['s_vs_z_route']:.2e}", flush=True)
    return pts, failures


def _stage_b1(comm, ans6, outdir):
    rank0 = comm.rank == 0
    failures = []
    t0 = time.perf_counter()
    unc = ans6._build_ladder(degree=1, sigmas=(ans6.COPPER_SIGMA,), resolution=None,
                             frequencies_hz=[ANCHOR1_FREQUENCY_HZ], measure_kappa=True)
    comm.Barrier()
    t_unc = time.perf_counter() - t0
    rec = unc["freqs"]["64"]
    cu, pec = rec["sigma"][ans6.COPPER_SIGMA], rec["pec"]
    tracked = json.loads((ans6.CASE_DIR / "metrics.json").read_text())
    view = {"frequency": "64 MHz",
            "cu": {"s_matrix": ans6._matrix_payload(np.asarray(cu["s"]))},
            "pec": {"s_matrix": ans6._matrix_payload(np.asarray(pec["s"]))}}
    worst, where, n = ans6._REPRO.knob_worst_deviation(view, tracked)
    band = ans6._REPRO.LEG_D0_REPRODUCTION_BAND
    cells_ok = int(unc["cells"]) == int(tracked["n_cells"])
    if rank0:
        print(f"[ans7] anchor 1: {n} S leaves at 64 MHz (uncorrected) vs tracked ans:6 "
              f"metrics.json ({unc['cells']} vs {tracked['n_cells']} cells): worst rel "
              f"{worst:.3e} at {where} (ASSERTED <= {band:.0e}); {t_unc:.1f} s", flush=True)
    if not (worst <= band and cells_ok):
        if rank0:
            print(f"[ans7] ANCHOR 1 RED -- stopping before the sweep: {(worst, where, cells_ok)}",
                  flush=True)
        return [("anchor-1", worst, where, cells_ok)]
    ds = cu["class_ds"]
    if rank0:
        print("[ans7] negative control (ASSERTED) |S(Cu)-S(PEC)| at 64 MHz, uncorrected: "
              + ", ".join(f"{c} {ds[c]:.6e}" for c in ans6.CLASSES)
              + f" (> {NEGATIVE_CONTROL_FLOOR:.0e})", flush=True)
    for c in ans6.CLASSES:
        if not ds[c] > NEGATIVE_CONTROL_FLOOR:
            failures.append(("negative-control", c, ds[c]))
    kappa = float(rec["kappa"]["pooled"])
    if rank0:
        print(f"[ans7] kappa (Cu P1, 64 MHz, uncorrected, pooled) = {kappa:.9e}; per sheet "
              f"{rec['kappa']['per_sheet']}; gate mesh's {GATE_MESH_KAPPA:.6e} (PRINTED); "
              f"told width factor 1/(1+kappa) = {1.0 / (1.0 + kappa):.9f}", flush=True)
    anchor = {"kappa": kappa, "kappa_per_sheet": rec["kappa"]["per_sheet"],
              "cells": int(unc["cells"]), "unknowns": int(unc["unknowns"]),
              "anchor1_worst_rel": float(worst), "anchor1_where": str(where),
              "negative_control_class_ds": {c: float(ds[c]) for c in ans6.CLASSES},
              "check_64mhz_uncorrected": {"cu_s": _mat(cu["s"]), "pec_s": _mat(pec["s"])},
              "t_uncorrected_s": t_unc}
    del unc
    t1 = time.perf_counter()
    lad = ans6._build_ladder(degree=1, sigmas=(ans6.COPPER_SIGMA,), resolution=None,
                             frequencies_hz=[f * 1e6 for f in STAGE_B1_MHZ],
                             width_correction_kappa=kappa)
    comm.Barrier()
    t_sweep = time.perf_counter() - t1
    pts, f2 = _point_records(lad, STAGE_B1_MHZ, ans6)
    failures += f2
    if rank0:
        payload = {"stage": "B1", "mpi_ranks": comm.size, "c_f": C_F, "anchor": anchor,
                   "points": pts, "t_sweep_s": t_sweep}
        (outdir / "B1.json").write_text(json.dumps(payload))
        print(f"[ans7] B1: {len(pts)} points, sweep {t_sweep:.1f} s at -n {comm.size}; "
              f"saved {outdir / 'B1.json'}", flush=True)
    return failures


def _stage_b2(comm, ans6, outdir):
    rank0 = comm.rank == 0
    b1 = json.loads((outdir / "B1.json").read_text())
    kappa = float(b1["anchor"]["kappa"])
    if rank0:
        print(f"[ans7] B2: kappa from B1 = {kappa:.9e}", flush=True)
    t1 = time.perf_counter()
    lad = ans6._build_ladder(degree=1, sigmas=(ans6.COPPER_SIGMA,), resolution=None,
                             frequencies_hz=[f * 1e6 for f in STAGE_B2_MHZ],
                             width_correction_kappa=kappa)
    comm.Barrier()
    t_sweep = time.perf_counter() - t1
    pts, failures = _point_records(lad, STAGE_B2_MHZ, ans6)
    if rank0:
        (outdir / "B2.json").write_text(json.dumps(
            {"stage": "B2", "mpi_ranks": comm.size, "kappa": kappa, "points": pts,
             "t_sweep_s": t_sweep}))
        print(f"[ans7] B2: {len(pts)} points, sweep {t_sweep:.1f} s at -n {comm.size}",
              flush=True)
    failures += _recompute_b1(b1)
    failures += _assemble(comm, b1, pts, t_sweep)
    return failures


def _assemble(comm, b1, pts_b2, t_b2):
    """Anchors (4)–(5), rows (a)–(d), predictions, 2·C control; metrics + COMPARISON."""
    failures = []
    pts = list(b1["points"]) + list(pts_b2)
    f_hz = np.array([p["f_hz"] for p in pts])
    if not (len(pts) == 51 and np.all(np.diff(f_hz) > 0)):
        failures.append(("window", len(pts)))
    out = {}
    for name in ("cu", "pec"):
        z_in = np.array([complex(p[name]["Z_in"]["re"], p[name]["Z_in"]["im"]) for p in pts])
        rows = parallel_rows(f_hz, z_in)
        n_up = len(rows["crossings_up_hz"])
        if n_up != 1:
            failures.append(("anchor-4", name, rows["crossings_up_hz"]))
        share = None
        if "f0_hz" in rows:
            j = int(np.argmin(np.abs(f_hz - rows["f0_hz"])))
            gaps = [complex(g["re"], g["im"]) for g in pts[j][name]["gap_currents"]]
            share = k1_share(gaps)
            rows["k_share"] = share
            rows["k_share_at_hz"] = float(f_hz[j])
            if not share > K_SHARE_FLOOR:
                failures.append(("anchor-5", name, share))
        else:
            failures.append(("anchor-5", name, "no f0"))
        # 2·C control (predicted): the same 4×4s reduced with 2·C on every gap.
        z_2c = np.array([route_b_s(_unmat(p[name]["s"]), p["f_hz"],
                                   CONTROL_C_FACTOR * C_F)[0] for p in pts])
        rows_2c = parallel_rows(f_hz, z_2c)
        rows["control_2c_crossings_up_hz"] = rows_2c["crossings_up_hz"]
        out[name] = rows
    rank0 = comm.rank == 0
    if rank0:
        for name, rows in out.items():
            if "f0_hz" in rows:
                print(f"[ans7] {name}: (a) f0 = {rows['f0_hz'] / 1e6:.6f} MHz (bracket "
                      f"{[b / 1e6 for b in rows['bracket_hz']]}), (b) R_in(f0) = "
                      f"{rows['r_in_ohm']:.6e} ohm, (c) Q = {rows['q_loaded']:.4f}; k = +-1 "
                      f"share at {rows['k_share_at_hz'] / 1e6:.1f} MHz = {rows['k_share']:.6f} "
                      f"(ASSERTED > {K_SHARE_FLOOR})", flush=True)
            print(f"[ans7] {name}: - -> + zeros of Im Y_in {[c / 1e6 for c in rows['crossings_up_hz']]}"
                  f" MHz (ASSERTED exactly 1); + -> - {[c / 1e6 for c in rows['crossings_down_hz']]}",
                  flush=True)
            print(f"[ans7] {name}: control 2*C (PREDICTED no - -> + zero in window): "
                  f"{[c / 1e6 for c in rows['control_2c_crossings_up_hz']]} -> "
                  f"{'met' if not rows['control_2c_crossings_up_hz'] else 'NOT met'}", flush=True)
        if "f0_hz" in out["pec"]:
            f0p = out["pec"]["f0_hz"] / 1e6
            lo, hi = PREDICTED_PEC_F0_MHZ
            print(f"[ans7] PREDICTED: PEC f0 {f0p:.6f} MHz in [{lo}, {hi}]: "
                  f"{'met' if lo <= f0p <= hi else 'NOT met'}; Q Cu "
                  f"{out['cu'].get('q_loaded', float('nan')):.3f} / PEC "
                  f"{out['pec'].get('q_loaded', float('nan')):.3f} vs ~{PREDICTED_Q:g}; "
                  f"k-share ~1: Cu {out['cu'].get('k_share')} PEC {out['pec'].get('k_share')}",
                  flush=True)
        if "f0_hz" in out["cu"] and "f0_hz" in out["pec"]:
            dfr = (out["cu"]["f0_hz"] - out["pec"]["f0_hz"]) / out["pec"]["f0_hz"]
            print(f"[ans7] (d) PRINTED, Cu vs PEC (not the AED row): Delta f0/f0 = {dfr:.6e}, "
                  f"Delta C/C equivalent = {-2.0 * dfr:.6e}", flush=True)
        _write_outputs(comm, b1, pts, out, t_b2)
    return failures


def _write_outputs(comm, b1, pts, rows, t_b2):
    metrics = {
        "case": "ans7_tuned_copper_birdcage_driven_sweep",
        "spec": "v2 (2026-09-28)",
        "route": "B (four-port + circuit layer); Route A is §9 item 68",
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "mpi_ranks": comm.size, "degree": 1, "n_cells": b1["anchor"]["cells"],
        "unknowns": b1["anchor"]["unknowns"], "c_f": C_F,
        "kappa": b1["anchor"]["kappa"], "kappa_per_sheet": b1["anchor"]["kappa_per_sheet"],
        "check_64mhz_uncorrected": b1["anchor"]["check_64mhz_uncorrected"],
        "anchor1_worst_rel": b1["anchor"]["anchor1_worst_rel"],
        "negative_control_class_ds": b1["anchor"]["negative_control_class_ds"],
        "terminated_indices": [int(k) for k in TUNING_TERMINATED_INDICES],
        "z0_ohm": float(REFERENCE_IMPEDANCE_OHM),
        "rows": rows,
        "points": [{"f_hz": p["f_hz"], **{c: {"S11_bare": p[c]["S11_bare"],
                                               "Z_in": p[c]["Z_in"]} for c in ("cu", "pec")}}
                   for p in pts],
        "t_sweep_s": {"B1": b1["t_sweep_s"], "B2": t_b2},
    }
    prev_path = CASE_DIR / "metrics.json"
    prev = json.loads(prev_path.read_text()) if prev_path.exists() else {}
    for key in ("degree2_exact", "superseded"):   # §9 item 230: never drop the degree-2 block
        if key in prev:
            metrics[key] = prev[key]
    (CASE_DIR / "metrics.json").write_text(json.dumps(metrics, indent=1) + "\n")
    lines = ["# ANS-7 v2 — comparison (Route B; AED and Route A columns blank)", "",
             "Generated by `07_tuned_copper_birdcage_driven_sweep.py` (stage B2). "
             "No AED value enters this file; the weekly adjudicates at a stated relative "
             "level (T10b).", "",
             f"κ (Cu P1, 64 MHz, uncorrected, pooled) = {b1['anchor']['kappa']:.9e}; "
             f"C = {C_F:.4e} F on every gap; degree 1; {b1['anchor']['cells']} cells.", "",
             "| Row | Ours Cu (B) | Ours PEC (B) | Ours Cu (A) | AED |",
             "|---|---|---|---|---|"]

    def cell(name, key, scale=1.0, fmt="{:.6f}"):
        v = rows[name].get(key)
        return fmt.format(v * scale) if v is not None else "—"

    lines.append(f"| (a) f₀ [MHz] | {cell('cu', 'f0_hz', 1e-6)} | {cell('pec', 'f0_hz', 1e-6)} |  |  |")
    lines.append(f"| (b) R_in(f₀) [Ω] | {cell('cu', 'r_in_ohm', fmt='{:.6e}')} | "
                 f"{cell('pec', 'r_in_ohm', fmt='{:.6e}')} |  |  |")
    lines.append(f"| (c) loaded Q | {cell('cu', 'q_loaded', fmt='{:.4f}')} | "
                 f"{cell('pec', 'q_loaded', fmt='{:.4f}')} |  |  |")
    lines.append("| (d) ΔC/C = −2Δf₀/f₀ (vs AED) |  |  |  |  |")
    lines.append(f"| k = ±1 gap-current share at f₀ | {cell('cu', 'k_share')} | "
                 f"{cell('pec', 'k_share')} | — |  |")
    lines += ["", "Not claimed: any absolute `Z_in`, a match.", ""]
    (CASE_DIR / "COMPARISON.md").write_text("\n".join(lines))
    print(f"[ans7] wrote {CASE_DIR / 'metrics.json'} and {CASE_DIR / 'COMPARISON.md'}",
          flush=True)


# -- stage A: Route A (in-model), the band, the private writer (§9 item 68) --------

#: The operator's AED half, if it has ever run on this box. Gitignored (ANS-1:
#: licence terms restrict disclosure); nothing read from here reaches a
#: tracked file.  Resolved at call time (tests repoint it and ``CASE_DIR``).
AED_RESULTS_PATH = CASE_DIR / "aed_results" / "ans7_aed_results.json"
#: `ANS-4`'s shared knob: ``=1`` skips the private write even if the JSON exists.
SKIP_PRIVATE_ENV = "FEM_EM_ANS_SKIP_PRIVATE"
#: κ read from the tracked metrics.json vs κ recomputed in-run (§9 item 68).
KAPPA_RECOMPUTE_ATOL = 1.0e-9
#: Control (i) (predicted, rule (e)): circuit at this factor · C vs in-model at C.
CONTROL_I_C_FACTOR = 0.5
CONTROL_I_FREQUENCY_HZ = 64.0e6
#: *Predicted*, printed only: control (i) ≥ 10× the band (gate mesh 1116.9×,
#: ``20260920T110503Z_PORT-22.log:2153``); control (ii) ≈ 1e-2
#: (``chunks/PORT-14.md:650``); band residuals 6–8e-5 (``…PORT-22.log:2144–2149``).
PREDICTED_CONTROL_I_RATIO_MIN = 10.0
PREDICTED_CONTROL_II = 1.0e-2
PREDICTED_BAND_RESIDUAL = (6.0e-5, 8.0e-5)
ROWS_PRINTED = ("f0_hz", "r_in_ohm", "q_loaded")


def _s11_from_z(z_bare):
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    return (z_bare - z0) / (z_bare + z0)


def _z_from_s11(s11):
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    return z0 * (1.0 + s11) / (1.0 - s11)


def _band_residual(s_a, s_b):
    """`PORT-15`'s ``_residual`` on the 1×1: |S_A − S_B| / |S_B|."""
    return float(abs(s_a - s_b) / abs(s_b))


def _route_a_solve(th14, geo, specs, f_hz, column, c_f=C_F, return_fields=False,
                   degree=1, esr_ohm=0.0):
    """One in-model drive of P1 with Z_C on P2–P4; ``S11_bare`` from power waves.

    Copied from `PORT-15`'s ``_terminated_kept_network``
    (``test_port_circuit_layer_field.py:522``), 1×1, plus the Cu Leontovich term.
    ``degree`` (§9 item 217, additive; default 1 byte-identical) is forwarded
    to ``run_lumped_sheet_port_case``.
    """
    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.lumped import run_lumped_sheet_port_case
    from fem_em_solver.ports.sparameters import _power_waves

    msh, cell_tags, facet_tags, outer_tags, tags_f, port_defs = geo
    omega = 2.0 * np.pi * f_hz
    z_c = _gap_z(f_hz, c_f, esr_ohm)   # esr_ohm = 0: the old 1/(jωC), bitwise (item 378)
    sp = [dataclasses.replace(s, port_impedance_ohm=z_c) if i in TUNING_TERMINATED_INDICES
          else s for i, s in enumerate(specs)]
    if column == "cu":
        problem = th14._problem(msh, cell_tags, f_hz, outer_tags, (th14.OUTER_BOX_TAG,))
        z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, th14.COPPER_SIGMA)
        extra = [th14._leontovich_term(msh, facet_tags, z_s, omega)]
    else:
        problem = th14._problem(msh, cell_tags, f_hz, facet_tags, None)
        extra = None
    out = run_lumped_sheet_port_case(problem, port_defs, sp, facet_tags=tags_f,
                                     driven_port_id=th14.DRIVEN, verbose=False,
                                     return_fields=return_fields,
                                     extra_bilinear_terms=extra, degree=degree)
    result, fields = out if return_fields else (out, None)
    drive = result.responses[th14.DRIVEN]
    a1, b1 = _power_waves(drive.voltage_v, drive.current_a, float(REFERENCE_IMPEDANCE_OHM))
    assert abs(a1) > 0.0, "incident wave at P1 vanished"
    return complex(b1 / a1), fields


def _kappa_in_run(th14, geo, specs_unc, comm):
    """κ recomputed exactly as `_build_ladder(measure_kappa=True)` does (Cu P1, 64 MHz,
    every sheet at its 50 Ω spec, uncorrected)."""
    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.lumped import run_lumped_sheet_port_case
    from fem_em_solver.ports.shares import terminal_form_deficit

    msh, cell_tags, facet_tags, outer_tags, tags_f, port_defs = geo
    f_hz = ANCHOR1_FREQUENCY_HZ
    omega = 2.0 * np.pi * f_hz
    problem = th14._problem(msh, cell_tags, f_hz, outer_tags, (th14.OUTER_BOX_TAG,))
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, th14.COPPER_SIGMA)
    term = th14._leontovich_term(msh, facet_tags, z_s, omega)
    _p1, fields = run_lumped_sheet_port_case(
        problem, port_defs, specs_unc, facet_tags=tags_f, driven_port_id=th14.DRIVEN,
        verbose=False, return_fields=True, extra_bilinear_terms=[term], degree=1)
    sheets = [sp.sheet(driven=(sp.port_id == th14.DRIVEN)) for sp in specs_unc]
    return float(terminal_form_deficit(msh, tags_f, sheets, fields.e_complex, comm)["pooled"])


def _write_xdmf(comm, msh, fields_by_col, outdir, f_hz):
    """Combined XDMF (Cu and PEC E at the point nearest f0), P1 vector, gitignored.

    P1, not DG1: XDMFFile needs the function's dof layout to match the mesh
    geometry (the first W3, ``20260929T095257Z_ANS-7-item68-W3-A.log``, died here
    on DG1 after all 102 Route A solves).
    """
    import basix.ufl
    from dolfinx import fem, io

    el = basix.ufl.element("Lagrange", msh.basix_cell(), 1, shape=(msh.geometry.dim,))
    v = fem.functionspace(msh, el)
    path = outdir / f"ans7_routeA_E_{f_hz / 1e6:.1f}MHz.xdmf"
    with io.XDMFFile(comm, str(path), "w") as xf:
        xf.write_mesh(msh)
        for name, fields in fields_by_col.items():
            u = fem.Function(v, name=f"E_{name}")
            u.interpolate(fields.e_complex)
            xf.write_function(u)
    return path


def _route_rows(f_hz, s11_bare):
    z_in = np.array([add_c1(_z_from_s11(s), 2.0 * np.pi * f) for f, s in zip(f_hz, s11_bare)])
    return parallel_rows(f_hz, z_in), z_in


def _rel(a, b):
    return float((a - b) / b) if (a is not None and b is not None) else None


def _route_comparison(rows_a, rows_b, z_a, z_b):
    """Route A − B, printed only: df0_rel / dq_rel / dr_in_rel / max_dzin_rel."""
    return {"df0_rel": _rel(rows_a.get("f0_hz"), rows_b.get("f0_hz")),
            "dq_rel": _rel(rows_a.get("q_loaded"), rows_b.get("q_loaded")),
            "dr_in_rel": _rel(rows_a.get("r_in_ohm"), rows_b.get("r_in_ohm")),
            "max_dzin_rel": float(np.max(np.abs(z_a - z_b) / np.abs(z_b)))}


def _load_aed_results(path=None):
    """The AED half if the JSON is on this box, else ``None``.

    Key layout per ``aed/ans7_hfss_pyaedt.py:44–61``: ``"<design>" →
    "<sweep>" → {"freq_hz", "Zin", "rows": {"f0_hz", "r_in_ohm", "q_loaded", ...},
    "k1_share_at_f0", ...}``, plus ``"_route_comparison"`` and
    ``"_column_difference"``; metadata keys start with ``_``.
    """
    path = AED_RESULTS_PATH if path is None else Path(path)
    if not path.exists():
        return None
    return json.loads(path.read_text())


def _private_comparison_text(metrics, aed) -> str:
    """The private table: ours (tracked metrics) beside every AED design/sweep's rows."""
    lines = ["# ANS-7 v2 — PRIVATE comparison (AED values; gitignored, never commit)", "",
             "| Source | (a) f₀ [MHz] | (b) R_in(f₀) [Ω] | (c) loaded Q | k-share |",
             "|---|---|---|---|---|"]

    def fmt(v, scale=1.0, f="{:.6f}"):
        return f.format(v * scale) if isinstance(v, (int, float)) else "—"

    for route in ("B", "A"):
        key = "rows" if route == "B" else "rows_route_a"
        for col in ("cu", "pec"):
            r = metrics.get(key, {}).get(col, {})
            lines.append(f"| ours {col} (route {route}) | {fmt(r.get('f0_hz'), 1e-6)} | "
                         f"{fmt(r.get('r_in_ohm'), f='{:.6e}')} | {fmt(r.get('q_loaded'), f='{:.4f}')} | "
                         f"{fmt(r.get('k_share'))} |")
    for design, block in aed.items():
        if design.startswith("_") or not isinstance(block, dict):
            continue
        for sweep, sw in block.items():
            if not isinstance(sw, dict) or "rows" not in sw:
                continue
            r = sw["rows"]
            lines.append(f"| AED {design} / {sweep} | {fmt(r.get('f0_hz'), 1e-6)} | "
                         f"{fmt(r.get('r_in_ohm'), f='{:.6e}')} | {fmt(r.get('q_loaded'), f='{:.4f}')} | "
                         f"{fmt(sw.get('k1_share_at_f0'))} |")
    lines += ["", "## Relative gaps (ours − AED) / AED", "",
              "| Ours | AED | df0_rel | dr_in_rel | dq_rel |", "|---|---|---|---|---|"]
    for design, block in aed.items():
        if design.startswith("_") or not isinstance(block, dict):
            continue
        for sweep, sw in block.items():
            if not isinstance(sw, dict) or "rows" not in sw:
                continue
            r = sw["rows"]
            for route, key in (("B", "rows"), ("A", "rows_route_a")):
                for col in ("cu", "pec"):
                    o = metrics.get(key, {}).get(col, {})
                    g = [_rel(o.get(k), r.get(k)) for k in ("f0_hz", "r_in_ohm", "q_loaded")]
                    lines.append(f"| {col} route {route} | {design} / {sweep} | "
                                 + " | ".join(fmt(x, f="{:+.3e}") for x in g) + " |")
    return "\n".join(lines) + "\n"


def _write_private_if_available(metrics) -> str:
    """`ANS-6`'s pattern (``06_…py:425–437``): knob first, then the JSON's presence."""
    if os.environ.get(SKIP_PRIVATE_ENV) == "1":
        return f"skipped ({SKIP_PRIVATE_ENV}=1) — private table not written"
    aed = _load_aed_results()
    if aed is None:
        return "no aed_results/ on this box — private table not written"
    path = CASE_DIR / "COMPARISON_private.md"
    path.write_text(_private_comparison_text(metrics, aed))
    return f"{path.name} — AED columns filled (untracked, never commit)"


def _stage_a(comm, outdir):
    """Route A at all 51 points, both columns; the band vs tracked Route B; controls."""
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_port_circuit_layer_field import REDUCTION_BAND

    rank0 = comm.rank == 0
    failures = []
    tracked = json.loads((CASE_DIR / "metrics.json").read_text())
    tracked_text = (CASE_DIR / "metrics.json").read_text()
    kappa = float(tracked["kappa"])
    pts_b = tracked["points"]
    f_hz = np.array([p["f_hz"] for p in pts_b])
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags, _diag, t_mesh = th14._sheets_build(True, as_hole=True)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    ncells = int(msh.topology.index_map(tdim).size_global)
    outer_tags, _census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_u, port_defs, specs_unc = th14._build_ports(msh, cell_tags, comm)
    tags_k, port_defs_k, specs_k = th14._build_ports(msh, cell_tags, comm,
                                                     width_correction_kappa=kappa)
    geo_u = (msh, cell_tags, facet_tags, outer_tags, tags_u, port_defs)
    geo_k = (msh, cell_tags, facet_tags, outer_tags, tags_k, port_defs_k)
    region_names = {2: "air", th14.PHANTOM_CELL_TAG: "phantom"}
    for spec in specs_unc:
        region_names[int(spec.facet_tag)] = f"port {spec.port_id}"
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{FIGURE_BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=(2,),
        translucent_tags=(th14.PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="ans:7 -- tuned copper birdcage F-small hole mesh, four gap ports (tag 401)",
        comm=comm,
    )
    comm.Barrier()
    if rank0:
        print(f"[ans7] A: hole mesh {ncells} cells (tracked {tracked['n_cells']}), mesh "
              f"{t_mesh:.1f} s, setup {time.perf_counter() - t0:.1f} s; kappa (tracked) = "
              f"{kappa:.9e}", flush=True)
    if ncells != int(tracked["n_cells"]):
        failures.append(("cells", ncells, tracked["n_cells"]))
    # κ recomputed in-run (asserted ≤ 1e-9 vs the tracked value).
    t1 = time.perf_counter()
    kappa_run = _kappa_in_run(th14, geo_u, specs_unc, comm)
    dk = abs(kappa_run - kappa)
    if rank0:
        print(f"[ans7] A: kappa in-run = {kappa_run:.12e}; |kappa_run - kappa_tracked| = "
              f"{dk:.3e} (ASSERTED <= {KAPPA_RECOMPUTE_ATOL:.0e}); {time.perf_counter() - t1:.1f} s",
              flush=True)
    if not dk <= KAPPA_RECOMPUTE_ATOL:
        failures.append(("kappa-recompute", kappa_run, kappa))
    # Control (ii) (predicted): one uncorrected Route A Cu solve at 64 MHz vs item
    # 69's uncorrected 4×4 reduced in the circuit layer.
    s11_unc_a, _ = _route_a_solve(th14, geo_u, specs_unc, ANCHOR1_FREQUENCY_HZ, "cu")
    s4_unc = _unmat(tracked["check_64mhz_uncorrected"]["cu_s"])
    _zin, zb_unc = route_b_s(s4_unc, ANCHOR1_FREQUENCY_HZ)
    ctrl_ii = _band_residual(s11_unc_a, _s11_from_z(zb_unc))
    # The 51-point Route A sweep, both columns (κ on all four sheets).
    s11_a = {"cu": [], "pec": []}
    res = {"cu": [], "pec": []}
    t_sweep0 = time.perf_counter()
    for k, p in enumerate(pts_b):
        for col in ("cu", "pec"):
            s_a, _ = _route_a_solve(th14, geo_k, specs_k, p["f_hz"], col)
            s_b = complex(p[col]["S11_bare"]["re"], p[col]["S11_bare"]["im"])
            r = _band_residual(s_a, s_b)
            s11_a[col].append(s_a)
            res[col].append(r)
            if not r <= REDUCTION_BAND:
                failures.append(("band", p["f_hz"], col, r))
        if rank0:
            print(f"[ans7] A f = {p['f_hz'] / 1e6:7.3f} MHz: S11_bare A Cu "
                  f"{s11_a['cu'][-1].real:+.9e}{s11_a['cu'][-1].imag:+.9e}j PEC "
                  f"{s11_a['pec'][-1].real:+.9e}{s11_a['pec'][-1].imag:+.9e}j; A-B residual Cu "
                  f"{res['cu'][-1]:.3e} PEC {res['pec'][-1]:.3e} (ASSERTED <= {REDUCTION_BAND:.0e})",
                  flush=True)
    comm.Barrier()
    t_sweep = time.perf_counter() - t_sweep0
    # Control (i) (predicted): circuit at 0.5·C (item 69's corrected 4×4 at 64.0 MHz,
    # the stage file B1.json) vs in-model at C at 64.0 MHz.
    ctrl_i = None
    b1_path = outdir / "B1.json"
    j64 = int(np.argmin(np.abs(f_hz - CONTROL_I_FREQUENCY_HZ)))
    if b1_path.exists():
        b1 = json.loads(b1_path.read_text())
        p64 = [q for q in b1["points"] if abs(q["f_hz"] - CONTROL_I_FREQUENCY_HZ) < 1.0][0]
        _z, zb_half = route_b_s(_unmat(p64["cu"]["s"]), CONTROL_I_FREQUENCY_HZ,
                                CONTROL_I_C_FACTOR * C_F)
        ctrl_i = _band_residual(s11_a["cu"][j64], _s11_from_z(zb_half))
    # Rows and the − → + zero count per column (asserted: exactly one).
    rows_a, cmp_ab = {}, {}
    for col in ("cu", "pec"):
        ra, za = _route_rows(f_hz, s11_a[col])
        zb = np.array([complex(p[col]["Z_in"]["re"], p[col]["Z_in"]["im"]) for p in pts_b])
        rows_a[col] = ra
        cmp_ab[col] = _route_comparison(ra, tracked["rows"][col], za, zb)
        if len(ra["crossings_up_hz"]) != 1:
            failures.append(("zero-count-A", col, ra["crossings_up_hz"]))
    # Combined XDMF at the point nearest the Route A Cu f0.
    xdmf = None
    if "f0_hz" in rows_a["cu"]:
        jf = int(np.argmin(np.abs(f_hz - rows_a["cu"]["f0_hz"])))
        flds = {col: _route_a_solve(th14, geo_k, specs_k, f_hz[jf], col, return_fields=True)[1]
                for col in ("cu", "pec")}
        xdmf = _write_xdmf(comm, msh, flds, outdir, f_hz[jf])
    comm.Barrier()
    if rank0:
        worst = max((r, c, f_hz[i]) for c in ("cu", "pec") for i, r in enumerate(res[c]))
        print(f"[ans7] A: band worst residual {worst[0]:.6e} ({worst[1]}, "
              f"{worst[2] / 1e6:.1f} MHz), ratio to REDUCTION_BAND {worst[0] / REDUCTION_BAND:.4f} "
              f"(ASSERTED all 102 <= {REDUCTION_BAND:.0e}; PREDICTED {PREDICTED_BAND_RESIDUAL}); "
              f"sweep {t_sweep:.1f} s", flush=True)
        for col in ("cu", "pec"):
            print(f"[ans7] A {col}: min/median/max residual {min(res[col]):.3e} / "
                  f"{float(np.median(res[col])):.3e} / {max(res[col]):.3e}", flush=True)
            ra = rows_a[col]
            print(f"[ans7] A {col}: - -> + zeros of Im Y_in {[c / 1e6 for c in ra['crossings_up_hz']]} "
                  f"MHz (ASSERTED exactly 1)", flush=True)
            print(f"[ans7] A {col} rows (PRINTED): (a) f0 {ra.get('f0_hz', float('nan')) / 1e6:.6f} MHz, "
                  f"(b) R_in {ra.get('r_in_ohm', float('nan')):.6e} ohm, (c) Q "
                  f"{ra.get('q_loaded', float('nan')):.4f}; A - B: " + ", ".join(
                      f"{k} {v:+.3e}" if v is not None else f"{k} n/a"
                      for k, v in cmp_ab[col].items()), flush=True)
        if "f0_hz" in rows_a["cu"] and "f0_hz" in rows_a["pec"]:
            dfr = (rows_a["cu"]["f0_hz"] - rows_a["pec"]["f0_hz"]) / rows_a["pec"]["f0_hz"]
            print(f"[ans7] A (d) PRINTED, Cu vs PEC: Delta f0/f0 = {dfr:.6e}, Delta C/C "
                  f"equivalent = {-2.0 * dfr:.6e}", flush=True)
        if ctrl_i is None:
            print(f"[ans7] control (i) NOT MEASURED: {b1_path} absent", flush=True)
        else:
            print(f"[ans7] control (i) (PREDICTED >= {PREDICTED_CONTROL_I_RATIO_MIN:g}x band): circuit "
                  f"at {CONTROL_I_C_FACTOR}*C vs in-model at C, 64.0 MHz Cu: {ctrl_i:.6e} = "
                  f"{ctrl_i / REDUCTION_BAND:.1f}x band -> "
                  f"{'met' if ctrl_i >= PREDICTED_CONTROL_I_RATIO_MIN * REDUCTION_BAND else 'NOT met'}",
                  flush=True)
        print(f"[ans7] control (ii) (PREDICTED ~{PREDICTED_CONTROL_II:.0e}): uncorrected Route A Cu "
              f"64 MHz vs circuit on item 69's uncorrected 4x4: {ctrl_ii:.6e} = "
              f"{ctrl_ii / REDUCTION_BAND:.1f}x band", flush=True)
        if xdmf is not None:
            print(f"[ans7] A: combined XDMF {xdmf}", flush=True)
    if failures or not rank0:
        return failures
    # Green: write Route A cells into metrics.json; Route B entries unchanged.
    metrics = json.loads(tracked_text)
    metrics["route"] = "A (in-model) and B (four-port + circuit layer)"
    metrics["kappa_in_run_route_a"] = kappa_run
    metrics["rows_route_a"] = rows_a
    metrics["route_comparison_a_minus_b"] = cmp_ab
    metrics["band"] = {"reduction_band": float(REDUCTION_BAND),
                       "residual_cu": res["cu"], "residual_pec": res["pec"],
                       "worst": max(max(res["cu"]), max(res["pec"]))}
    metrics["controls"] = {"i_half_c_64mhz": ctrl_i, "ii_uncorrected_64mhz": ctrl_ii}
    for p, sa, sp in zip(metrics["points"], s11_a["cu"], s11_a["pec"]):
        p["cu"]["S11_bare_route_a"] = _cplx(sa)
        p["pec"]["S11_bare_route_a"] = _cplx(sp)
    metrics["t_route_a_sweep_s"] = t_sweep
    (CASE_DIR / "metrics.json").write_text(json.dumps(metrics, indent=1) + "\n")
    print(f"[ans7] A: wrote Route A cells into {CASE_DIR / 'metrics.json'} (Route B entries "
          f"unchanged); private: {_write_private_if_available(metrics)}", flush=True)
    return failures


# -- §9 item 230 (`ANS-7` step 9, §10 T20): degree 2 on exact rectangles --------

P2_DEFAULT_DIR = "/workspace/logs/ans7_p2"
X_NODE_STAGES = {"X1": (62.0, 62.3), "X2": (63.2, 64.5), "X3": (65.8, 66.7), "X4": (67.0,)}
X_LOSS_NODE_MHZ = 64.5
XA_MHZ = (62.6, 62.7, 62.8)
R2X_RANKS = 8
#: Printed predictions (§9 item 230): bands, never asserted.
PRED_R2X_MAX = 1.0e-11
PRED_AB = (5.0e-6, 3.0e-4)
PRED_DF0_A_B = (-1.0e-5, -1.0e-6)
PRED_L = (1.0e-6, 1.0e-4)
PRED_ROUTE_A_ROWS_REL = 0.02
SUPERSEDED_KEYS = ("rows", "rows_route_a", "points", "kappa", "kappa_per_sheet", "n_cells",
                   "degree", "band", "controls", "route_comparison_a_minus_b")
D2_KEYS = ("degree2_exact", "superseded")


def _p2_dir() -> Path:
    raw = os.environ.get(RESULTS_DIR_ENV, "")
    return Path(raw) if raw else Path(P2_DEFAULT_DIR)


def _say0(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _x_node_path(outdir: Path, mhz: float) -> Path:
    return outdir / f"node_{mhz:.1f}.json"


def _x_guard():
    knobs = sorted(k for k in os.environ
                   if k.startswith("FEM_EM_ANS7_EXACT_") or k.startswith("FEM_EM_ANS7_SPARSE_"))
    if knobs:
        raise RuntimeError(f"an X stage needs every FEM_EM_ANS7_EXACT_* / FEM_EM_ANS7_SPARSE_* "
                           f"knob unset; found {knobs}")


def _x_nodes_stage(comm, stage, outdir):
    """X1–X4: item 196's per-node call at degree 2, exact rectangles, κ off."""
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans7_exact_strip import (
        _assert_tie, _cj, _fixture_tie, _im_y_in, _tie_json)
    from tests.validation.test_ans7_sparse_route_b import (
        DEFAULT_NODES_MHZ, _loss_pair, _p_inc_w, z_from_s)
    from tests.validation.test_port_birdcage_four_port import PASSIVITY_SIGMA_TOLERANCE
    from tests.validation.test_port_circuit_layer_field import REDUCTION_BAND
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    rank0 = comm.rank == 0
    failures = []
    nodes = X_NODE_STAGES[stage]
    if stage == "X4":
        missing = sorted({g for g in ("X1", "X2", "X3") for m in X_NODE_STAGES[g]
                          if not _x_node_path(outdir, m).exists()})
        if missing:
            raise RuntimeError(f"X4 needs the node files of {missing} in {outdir}; run those "
                               f"stages first")
    tie = _fixture_tie() if stage == "X1" else None
    field = X_LOSS_NODE_MHZ in nodes
    t0 = time.perf_counter()
    lad = th14._build_ladder(2, sigmas=(th14.COPPER_SIGMA,),
                             frequencies_hz=[m * 1e6 for m in nodes],
                             width_correction_kappa=None, measure_kappa=False,
                             pec_column=True, field_solve=field,
                             leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
    comm.Barrier()
    t_solve = time.perf_counter() - t0
    recs = sorted(lad["freqs"].values(), key=lambda r: r["f_hz"])
    assert len(recs) == len(nodes)
    cells, unknowns = int(lad["cells"]), int(lad["unknowns"])
    if tie is not None:
        _assert_tie(cells, tie)           # (F)
        _say0(comm, f"[ans7-x] X1 (F) fixture tie green: cells {cells}")
    for mhz, rec in zip(nodes, recs):
        assert abs(rec["f_hz"] - mhz * 1e6) < 1.0
        s, z, recip, smax, im_y = {}, {}, {}, {}, {}
        for c, sw in (("cu", rec["sigma"][th14.COPPER_SIGMA]), ("pec", rec["pec"])):
            s4 = np.asarray(sw["s"], dtype=np.complex128)
            s[c], z[c] = s4, z_from_s(s4)
            recip[c] = float(sw["reciprocity"])
            smax[c] = float(np.max(sw["sigma"]))
            im_y[c], z_in = _im_y_in(z[c], mhz)
            _say0(comm, f"[ans7-x] {stage} node {mhz:.1f} MHz {c}: Im Y_in {im_y[c]:+.9e} S "
                        f"(Z_in {z_in.real:.6e}{z_in.imag:+.6e}j ohm); recip {recip[c]:.3e} "
                        f"(ASSERTED <= {RECIPROCITY_BAND:.0e}); sigma_max {smax[c]:.12f} "
                        f"(ASSERTED <= 1 + {PASSIVITY_SIGMA_TOLERANCE:.0e})")
            if not recip[c] <= RECIPROCITY_BAND:
                failures.append(("reciprocity", mhz, c, recip[c]))
            if not smax[c] <= 1.0 + PASSIVITY_SIGMA_TOLERANCE:
                failures.append(("passivity", mhz, c, smax[c]))
        power = dict(rec["power"]) if field else None
        ident = float(rec["identity_residual"]) if field else None
        l_pair = None
        if field:
            l_s, l_f = _loss_pair(s["cu"], power)
            gap = abs(l_s - l_f)
            l_pair = {"l_s": float(l_s), "l_f": float(l_f), "gap": float(gap),
                      "p_inc_w": float(_p_inc_w())}
            lo, hi = PRED_L
            _say0(comm, f"[ans7-x] {stage} node {mhz:.1f} MHz copper identity residual {ident:.3e} "
                        f"(ASSERTED <= {th14.DISCRETE_IDENTITY_RTOL:g}); (L) L_S {l_s:.9e}, L_F "
                        f"{l_f:.9e}, |L_S - L_F| = {gap:.6e} (ASSERTED <= {REDUCTION_BAND:.0e}; "
                        f"PREDICTED [{lo:g}, {hi:g}]: {'met' if lo <= gap <= hi else 'NOT met'})")
            if not ident <= th14.DISCRETE_IDENTITY_RTOL:
                failures.append(("identity", mhz, ident))
            if not gap <= REDUCTION_BAND:
                failures.append(("L", mhz, gap))
        if rank0:
            _x_node_path(outdir, mhz).write_text(json.dumps({
                "f_hz": rec["f_hz"], "degree": 2, "ranks": comm.size, "cells": cells,
                "unknowns": unknowns, "kappa": None, "stage": stage,
                "tie": None if tie is None else _tie_json(tie), "field_solve": field,
                "s": {c: _cj(s[c]) for c in s}, "z": {c: _cj(z[c]) for c in z},
                "recip": recip, "sigma_max": smax, "im_y_in": im_y, "power": power,
                "identity_residual": ident, "loss_identity": l_pair,
                "t_cu_s": float(rec["sigma"][th14.COPPER_SIGMA]["t"]),
                "t_pec_s": float(rec["t_pec"]), "t_stage_solve_s": t_solve}))
    comm.Barrier()
    _say0(comm, f"[ans7-x] {stage}: nodes {nodes} solved in {t_solve:.1f} s, cells {cells}, "
                f"unknowns {unknowns}, degree 2, -n {comm.size}; files in {outdir}")
    if stage == "X4":
        failures += _x_assemble(comm, outdir, t_solve, DEFAULT_NODES_MHZ)
    return failures


def _snapshot(d, skip=()):
    return {k: json.dumps(v, sort_keys=True) for k, v in d.items() if k not in skip}


def _x_assemble(comm, outdir, t_x4, nodes_mhz):
    """X4's assembly: (R2x), Δ_LOO, crossings, k-share, (N-deg), (K); metrics + COMPARISON."""
    from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
    from tests.validation.test_ans6_exact_strip import ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    from tests.validation.test_ans7_exact_strip import (
        ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2, ROWS, _loo, _uj)
    from tests.validation.test_ans7_sparse_route_b import (
        BANDS, GRID_MHZ, K_SHARE_FLOOR, N_FACTOR, NODE_S_RTOL)

    failures = []
    recs = [json.loads(_x_node_path(outdir, m).read_text()) for m in nodes_mhz]
    cells = {r["cells"] for r in recs}
    if len(cells) != 1 or not all(abs(r["cells"] / ANS6_EXACT_STRIP_HOLE_CELL_RECORD - 1.0)
                                  < CELL_COUNT_BAND for r in recs):
        failures.append(("cells", sorted(cells), ANS6_EXACT_STRIP_HOLE_CELL_RECORD))
    if not all(r["degree"] == 2 and r["kappa"] is None for r in recs):
        failures.append(("degree/kappa", [(r["degree"], r["kappa"]) for r in recs]))
    z_nodes = {c: [_uj(r["z"][c]) for r in recs] for c in ("cu", "pec")}
    full, worst = _loo(list(nodes_mhz), z_nodes)
    rows = {c: dict(full[c]["rows"]) for c in ("cu", "pec")}
    rank0 = comm.rank == 0
    r2x = {c: {} for c in rows}
    for c in rows:
        r = rows[c]
        if len(r.get("crossings_up_hz", [])) != 1:
            failures.append(("crossings", c, r.get("crossings_up_hz")))
        if not r.get("k_share", 0.0) > K_SHARE_FLOOR:
            failures.append(("k-share", c, r.get("k_share")))
        if rank0:
            print(f"[ans7-x] X4 {c} rows: f0 {r.get('f0_hz')!r} Hz, R_in {r.get('r_in_ohm')!r} ohm, "
                  f"Q {r.get('q_loaded')!r}, k-share {r.get('k_share')!r}; up-crossings "
                  f"{r.get('crossings_up_hz')} (ASSERTED exactly 1; k-share > {K_SHARE_FLOOR})",
                  flush=True)
        for k in ROWS:
            want = ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2[c][k]
            r2x[c][k] = float(abs(r[k] / want - 1.0)) if k in r else float("inf")
    for k in ROWS:
        if rank0:
            print(f"[ans7-x] X4 Delta_LOO {k} = {worst[k]:.3e} (ASSERTED <= {BANDS[k]:.0e}; "
                  f"PREDICTED 6.3e-13 / 1.6e-9 / 1.6e-9)", flush=True)
        if not worst[k] <= BANDS[k]:
            failures.append(("delta-loo", k, worst[k]))
    worst_r2x = max(v for c in r2x for v in r2x[c].values())
    if rank0:
        for c in r2x:
            print(f"[ans7-x] X4 (R2x) {c}: " + ", ".join(f"{k} {v:.3e}" for k, v in r2x[c].items())
                  + f" vs ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2 (PREDICTED <= {PRED_R2X_MAX:.0e}); "
                  f"ASSERTED <= NODE_S_RTOL {NODE_S_RTOL:.0e} at -n {R2X_RANKS}", flush=True)
    if comm.size == R2X_RANKS:
        if not worst_r2x <= NODE_S_RTOL:
            failures.append(("R2x", r2x))
    elif rank0:
        print(f"[ans7-x] X4 (R2x) printed, not asserted: -n {comm.size} != {R2X_RANKS} "
              f"(EX-56 rule)", flush=True)
    tracked_text = (CASE_DIR / "metrics.json").read_text()
    tracked = json.loads(tracked_text)
    f0_sup = float(tracked["rows"]["cu"]["f0_hz"])
    n_deg = float(abs(rows["cu"]["f0_hz"] / f0_sup - 1.0)) if "f0_hz" in rows["cu"] else 0.0
    n_deg_floor = N_FACTOR * BANDS["f0_hz"]
    if rank0:
        print(f"[ans7-x] X4 (N-deg) |f0,X(Cu)/f0,superseded(Cu) - 1| = {n_deg:.6e} "
              f"(superseded {f0_sup!r} Hz) = {n_deg / n_deg_floor:.1f} x floor "
              f"(ASSERTED > N_FACTOR {N_FACTOR:g} x {BANDS['f0_hz']:.0e})", flush=True)
    if not n_deg > n_deg_floor:
        failures.append(("N-deg", n_deg))
    rec64 = recs[list(nodes_mhz).index(X_LOSS_NODE_MHZ)]
    loss = rec64["loss_identity"]
    sec = {}
    for r in recs:
        if r["stage"] != "X4":
            sec[r["stage"]] = r["t_stage_solve_s"]
    sec["X4"] = float(t_x4)
    d2 = {"degree": 2, "kappa": None,
          "sheets": "exact CAD rectangles (operator ruling 2026-10-02)",
          "cells": int(recs[0]["cells"]), "unknowns": int(recs[0]["unknowns"]),
          "mpi_ranks": int(comm.size), "c_f": C_F, "nodes_mhz": list(nodes_mhz),
          "grid_mhz": [GRID_MHZ[0], GRID_MHZ[-1], len(GRID_MHZ)],
          "rows": rows, "delta_loo": {k: float(v) for k, v in worst.items()},
          "r2x_rel": r2x, "loss_identity_64p5": loss, "n_deg": n_deg,
          "t_stage_s": sec,
          "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
          "record": "ANS7_EXACT_STRIP_ROUTE_B_ROWS_P2 (item 196 W3) reproduced through the example"}
    sup = {"keys": list(SUPERSEDED_KEYS),
           "statement": "superseded by degree2_exact (§10 T20): degree 1, ragged sheets, "
                        "kappa-corrected (stage B2, 2026-09-29); kept in place, unchanged, for "
                        "their readers",
           "readers": ["test_ans7_sparse_route_b", "test_wf5_q_two_ways", "test_ans6_exact_strip",
                       "test_port23_skrf_route_b"]}
    if failures:
        return failures
    if rank0:
        snap = _snapshot(tracked, skip=D2_KEYS)
        new = dict(tracked)
        new["degree2_exact"], new["superseded"] = d2, sup
        (CASE_DIR / "metrics.json").write_text(json.dumps(new, indent=1) + "\n")
        after = json.loads((CASE_DIR / "metrics.json").read_text())
        bad = [k for k, v in snap.items() if json.dumps(after.get(k), sort_keys=True) != v]
        print(f"[ans7-x] X4 (K): {len(snap)} pre-existing keys re-serialised identically: "
              f"{'yes' if not bad else 'NO ' + str(bad)} (ASSERTED)", flush=True)
        if bad:
            (CASE_DIR / "metrics.json").write_text(tracked_text)
            failures.append(("K", bad))
        else:
            _write_comparison(after)
            print(f"[ans7-x] X4: wrote {CASE_DIR / 'metrics.json'} and "
                  f"{CASE_DIR / 'COMPARISON.md'}", flush=True)
    return failures


# §9 item 514 (§10 T47): the closing line, shared by the writer and both before-tail inserters
# (ruling R1 (A), 2026-10-18 03:00 review). Re-pointed from "no h-statement (T12′); no absolute
# f₀ / Q / R_in", stale since T13 (f₀) and T41 (Q / R_in).
CLOSING_LINE = ("Not claimed: any absolute `Z_in`, a match; an absolute Q or R_in (INCONCLUSIVE, "
                "T48); the tables above are one mesh each — the h-statements are the verdict "
                "section's.")

# §9 item 514 (§10 T47): the public verdict section, rendered above every table. Every figure is
# a literal traced by tests/unit/test_ans7_verdict_section.py. T30 =
# docs/testing/logs/20261011T094910Z_ANS-7-item364-T30.log. No AED value, no Palace sentence.
VERDICT_HEAD = "## Verdict (public; adjudicated by the weekly review)"
VERDICT_SECTION = (
    VERDICT_HEAD,
    "",
    # cells: T30:51 / :52 / :53; h: T30:86
    "- **What was adjudicated:** our degree-2 Route B Cu rows on the exact SPEC rectangles, κ off, "
    "C = 66.29 pF on every gap — the first table below, on the ×1 mesh (h = 0.015 m, 80 798 "
    "cells), and the same rows on two refined rungs (h = 0.01 m, 140 270 cells; h = 0.005 m, "
    "523 880 cells; `tests/validation/test_ans7_t30_rin_q_sequence.py`) — against AED First "
    "Order. The degree-1, C_in, ESR and cross-solver sections below are records, not the "
    "adjudicated pairing.",
    # SPEC.md:50 / :51–52
    "- **Mode and copper contribution** (2026-09-30 weekly, degree 1): AGREE — one k = ±1 "
    "parallel resonance in the window in both codes (share ≥ 0.9999); the conductor's part of "
    "1/Q, `Q_Cu/Q_PEC` and the Cu − PEC shift of f₀ within 5 %.",
    # PROJECT_PLAN.md:316–317 (0.11 % / 0.03 %); T30:93 (values, changes, ratio), :94 (p);
    # Richardson limit at the printed p (recomputed in the test)
    "- **f₀** (T13, 2026-10-10 weekly): AGREE with an h-statement — within 0.11 % of AED First "
    "Order at ×1 and 0.03 % at h = 0.005. Ours: 62.70108 / 62.72685 / 62.75178 MHz on "
    "h = 0.015 / 0.01 / 0.005 m; changes +25.76 / +24.94 kHz, ratio 0.968, printed rate "
    "p = 1.0625 (first order in h); Richardson limit ≈ 62.775 MHz (+0.04 % from h = 0.005).",
    # T30:89 / :90 (Q), :85 / :86 (R_in), :87 / :91 (BRANCH lines, retired rule)
    "- **Q and R_in** (T41, 2026-10-17 weekly): converge in h. Q 977.26 / 970.44 / 965.55 and "
    "R_in 18 441.75 / 18 313.60 / 18 223.23 Ω on h = 0.015 / 0.01 / 0.005 m; successive-change "
    "ratios 0.717 (Q) / 0.705 (R_in) on these steps, where a first-order quantity gives 1.00 and "
    "a second-order one 0.60; printed rates p = 1.6451 (Q) / 1.6786 (R_in); Richardson limits at "
    "those rates Q∞ ≈ 963.25 (−0.24 % from h = 0.005, −1.43 % from ×1), R_in∞ ≈ 18 182 Ω "
    "(−0.23 %). The T30 log's BRANCH lines (`20261011T094910Z_ANS-7-item364-T30.log:87`, `:91`) "
    "print the retired 'second change at most half the first' test; on these steps first order "
    "itself gives 1.00, so that test does not read the sequence and was retired (2026-10-17 "
    "weekly epitaph).",
    # PROJECT_PLAN.md:4682–4687 (Q −3.1 %, R_in −2.7 %, ≈ 3 %, INCONCLUSIVE, T48)
    "- **Against AED:** at h = 0.005 Q sits 3.1 % and R_in 2.7 % below AED First Order (T13), "
    "and the limits sit ≈ 3 % below it. The absolute Q / R_in is **INCONCLUSIVE — not an "
    "agreement and not a disagreement**: AED First Order is one discretisation whose own "
    "pass-to-pass trend is unread (T48, the operator's). `ANS-7` ✅ on its Done-when (Q's "
    "h-statement).",
    "- Not claimed: any AED value; an absolute Q, R_in, `Z_in` or match; tuning accuracy at C_in "
    "or with ESR; human scale (`ANS-8`).",
    "",
)


def _write_comparison(m):
    """COMPARISON.md from metrics.json: degree-2 table first, the superseded table below."""
    d2 = m["degree2_exact"]
    ra_rows = d2.get("route_a", {}).get("rows", {})

    def cell(rows, key, scale=1.0, fmt="{:.6f}"):
        v = rows.get(key)
        return fmt.format(v * scale) if v is not None else "—"

    def table(rows_cu, rows_pec, rows_a):
        return ["| Row | Ours Cu (B) | Ours PEC (B) | Ours Cu (A) | AED |", "|---|---|---|---|---|",
                f"| (a) f₀ [MHz] | {cell(rows_cu, 'f0_hz', 1e-6)} | {cell(rows_pec, 'f0_hz', 1e-6)} | "
                f"{cell(rows_a, 'f0_hz', 1e-6) if rows_a else ''} |  |",
                f"| (b) R_in(f₀) [Ω] | {cell(rows_cu, 'r_in_ohm', fmt='{:.6e}')} | "
                f"{cell(rows_pec, 'r_in_ohm', fmt='{:.6e}')} | "
                f"{cell(rows_a, 'r_in_ohm', fmt='{:.6e}') if rows_a else ''} |  |",
                f"| (c) loaded Q | {cell(rows_cu, 'q_loaded', fmt='{:.4f}')} | "
                f"{cell(rows_pec, 'q_loaded', fmt='{:.4f}')} | "
                f"{cell(rows_a, 'q_loaded', fmt='{:.4f}') if rows_a else ''} |  |",
                "| (d) ΔC/C = −2Δf₀/f₀ (vs AED) |  |  |  |  |",
                f"| k = ±1 gap-current share at f₀ | {cell(rows_cu, 'k_share')} | "
                f"{cell(rows_pec, 'k_share')} | — |  |"]

    lines = ["# ANS-7 v2 — comparison (Route B; AED column blank)", "",
             "Generated by `07_tuned_copper_birdcage_driven_sweep.py` (stages X4 / XA, "
             "§9 item 230). No AED value enters this file; the weekly adjudicates at a "
             "stated relative level (T10b).", "",
             *VERDICT_SECTION,  # §9 item 514: above every table, untouched by the inserters
             "## Degree 2, exact CAD rectangles, κ off (current; §10 T20)", "",
             f"C = {m['c_f']:.4e} F on every gap; degree 2; {d2['cells']} cells; seven Route B "
             f"nodes on the 51-point grid; -n {d2['mpi_ranks']}. Route A (Cu) is the three-point "
             f"in-model drive at 62.6 / 62.7 / 62.8 MHz (blank until stage XA).", ""]
    lines += table(d2["rows"]["cu"], d2["rows"]["pec"], ra_rows)
    lines += ["", "## Superseded — degree 1, ragged sheets, κ-corrected (stage B2, 2026-09-29)", "",
              f"κ (Cu P1, 64 MHz, uncorrected, pooled) = {m['kappa']:.9e}; C = {m['c_f']:.4e} F on "
              f"every gap; degree 1; {m['n_cells']} cells.", ""]
    lines += table(m["rows"]["cu"], m["rows"]["pec"], m.get("rows_route_a", {}).get("cu"))
    tl = d2.get("tuned_larmor")
    if tl:
        def tcell(key, fn):
            e = tl.get(key)
            try:
                return fn(e) if e else "—"
            except (KeyError, TypeError):
                return "—"
        ks = ("64.0", "127.74")
        lines += ["", "## Tuned to the Larmor frequency: C_in, the in-model zero (×1 mesh, degree 2; "
                  "§10 T23) — lossless capacitors (r = 0)", "",
                  "Route A (Cu, in-model) at C_in at three tune points per field; κ off; no AED "
                  "column.", "",
                  "| Row | 64.000 MHz | 127.74 MHz |", "|---|---|---|"]
        for lab, fn in (
                ("C_in [pF]", lambda e: f"{e['c_in_f'] * 1e12:.6f}"),
                ("f₀,A [MHz]", lambda e: f"{e['rows_a']['f0_hz'] * 1e-6:.6f}"),
                ("f₀,A/target − 1", lambda e: f"{e['f0a_rel_target']:+.6e}"),
                ("bands used", lambda e: f"{e['bands_fraction']:.3f}"),
                ("c = C_mode/C_in − 1", lambda e: f"{e['c_rel']:+.6e}"),
                ("r = f₀,A/f₀,B(C_in) − 1", lambda e: f"{e['r']:+.6e}"),
                ("(AB) max", lambda e: f"{max(e['residual_ab']):.3e}"),
                ("R_in [Ω]", lambda e: f"{e['rows_a']['r_in_ohm']:.6e}"),
                ("Q", lambda e: f"{e['rows_a']['q_loaded']:.4f}")):
            lines.append(f"| {lab} | " + " | ".join(tcell(k, fn) for k in ks) + " |")
    te = d2.get("tuned_larmor_esr", {}).get("64.0")
    if te:
        lines += ["", "## Tuned with capacitor ESR (r = 0.2 Ω on every capacitor; §10 T26 / T27)", "",
                  "Route A (Cu, in-model) at (r, C_in(r)) at three tune points, the series-RC sheet "
                  "on P2–P4 and the feed gap's (r + C) branch in parallel; ×1 mesh, degree 2, κ off; "
                  "Q_0 is the lossless-capacitor row above. No AED column: no AED ESR case exists.",
                  "", "| Row | 64.000 MHz |", "|---|---|"]
        for lab, val in (
                ("C_in(r) [pF]", f"{te['c_in_r_f'] * 1e12:.6f}"),
                ("C_in(r)/C_in(0) − 1", f"{te['c_in_rel']:+.6e}"),
                ("f₀,A [MHz]", f"{te['rows_a']['f0_hz'] * 1e-6:.6f}"),
                ("f₀,A/64 MHz − 1", f"{te['f0a_rel_target']:+.6e}"),
                ("bands used", f"{te['bands_fraction']:.3f}"),
                ("(AB) max", f"{max(te['residual_ab']):.3e}"),
                ("R_in [Ω]", f"{te['rows_a']['r_in_ohm']:.6e}"),
                ("Q(r)", f"{te['rows_a']['q_loaded']:.4f}"),
                ("Q(r)/Q_0", f"{te['q_ratio']:.6f}"),
                ("closed-form Q (1/Q = 1/Q_0 + 1/Q_C)", f"{te['q_closed_form']['q']:.4f}")):
            lines.append(f"| {lab} | {val} |")
    lines += ["", CLOSING_LINE, ""]
    (CASE_DIR / "COMPARISON.md").write_text("\n".join(lines))
    if "127.74" in m["degree2_exact"].get("tuned_larmor_esr", {}):
        _write_comparison_esr128(m)  # §9 item 400: keep XT128R's sibling section on a rewrite
    _write_comparison_palace()  # §9 item 460: the public Palace section, before the tail


def _write_comparison_esr128(m):
    """XT128R only (§9 item 392): insert the 127.74 MHz ESR sibling section into COMPARISON.md,
    directly after item 378's 64 MHz ESR section and before the closing 'Not claimed' line.
    Every pre-existing character of the file is kept; an earlier copy of this section is replaced."""
    te = m["degree2_exact"]["tuned_larmor_esr"]["127.74"]
    head = "## Tuned with capacitor ESR at 127.74 MHz (§10 T35 / T36)"
    sec = ["", head, "",
           "Route A (Cu, in-model) at (r, C_in(r)) at three tune points, the series-RC sheet "
           "on P2–P4 and the feed gap's (r + C) branch in parallel; ×1 mesh, degree 2, κ off; "
           "Q_0 is the lossless-capacitor row above. No AED column: no AED ESR case exists.",
           "", "| Row | 127.74 MHz |", "|---|---|"]
    for lab, val in (
            ("C_in(r) [pF]", f"{te['c_in_r_f'] * 1e12:.6f}"),
            ("C_in(r)/C_in(0) − 1", f"{te['c_in_rel']:+.6e}"),
            ("f₀,A [MHz]", f"{te['rows_a']['f0_hz'] * 1e-6:.6f}"),
            ("f₀,A/127.74 MHz − 1", f"{te['f0a_rel_target']:+.6e}"),
            ("bands used", f"{te['bands_fraction']:.3f}"),
            ("(AB) max", f"{max(te['residual_ab']):.3e}"),
            ("R_in [Ω]", f"{te['rows_a']['r_in_ohm']:.6e}"),
            ("Q(r)", f"{te['rows_a']['q_loaded']:.4f}"),
            ("Q(r)/Q_0", f"{te['q_ratio']:.6f}"),
            ("closed-form Q (1/Q = 1/Q_0 + 1/Q_C)", f"{te['q_closed_form']['q']:.4f}")):
        sec.append(f"| {lab} | {val} |")
    path = CASE_DIR / "COMPARISON.md"
    text = path.read_text()
    tail = "\n" + CLOSING_LINE + "\n"
    if head in text:
        i = text.index("\n" + head)
        j = text.index(tail, i)
        text = text[:i] + text[j:]
    assert text.endswith(tail), "COMPARISON.md closing line not found"
    k = len(text) - len(tail)
    path.write_text(text[:k] + "\n".join(sec) + "\n" + text[k:])
    _write_comparison_palace()  # §9 item 460: a stand-alone XT128R rewrite keeps the Palace section


# §9 item 460 (§10 T45): the public "Palace at matched order" section. Every figure is a literal
# copied from a tracked read log (docs/testing/logs/, cited <log>:<line>); traced by
# tests/unit/test_ans7_palace_section.py. No AED column: no AED ESR case exists (§5.4 (iv)).
PALACE_HEAD = "## Palace at matched order (public; §10 T39–T44, `OPS-61`)"
PALACE_SECTION = (
    "",
    PALACE_HEAD,
    "",
    "Palace v0.18.1, order 2 — the same N1curl degree as ours, `ND (p = 2)` 570 332 on the ×1 hole "
    "mesh (80 798 elements) and 968 016 on the h = 0.01 hole mesh (140 270 elements), equal to our "
    "unknown counts — on the meshes our Route A solves, exported and configured by the tracked "
    "`scripts/probes/ops61_step1a_palace_ans4.py` (`export --hole-pec`, `export --hole-h001`; "
    "configs `ans7-configs`, `ans7-sweep-configs`, `ans7-128-configs`, `ans7-h001-configs [--f128]`; "
    "reads `read --ans7-esr`, `--ans7-sweep`, `--ans7-128`, `--ans7-h001 [--f128]`). The capacitor "
    "ESR r = 0.2 Ω enters Palace as its lumped R ‖ C — one exact (R_p, C_p) map per frequency (ESR "
    "ruling (iii); Palace's lumped R, C read as parallel, T39) — and the feed gap's (r + C) branch "
    "is added in post-processing on both sides. f₀, Q and R_in come from `parallel_rows` on the "
    "same three tune points (63.9 / 64.0 / 64.1 and 127.64 / 127.74 / 127.84 MHz) on both sides. "
    "No AED column: no AED ESR case exists.",
    "",
    "Table A — ties, Palace / ours − 1 (signed; S₁₁ columns are |S_Pal − S_ours| / |S_ours|):",
    "",
    "| Rung, frequency | f₀ | Q | R_in | S₁₁ per tune point (max of 3) | lossless-C S₁₁ "
    "| ESR-vs-lossless move |",
    "|---|---|---|---|---|---|---|",
    # 421 W4:48 / :49 / :50 / max(:38–40); 420 W6:46 / :48
    "| ×1, 64.000 MHz | +6.530e-11 | -4.871e-09 | -7.456e-10 | 2.580e-10 | 2.565e-10 | 1.994e-09 |",
    # 436 W6:51 / :52 / :53 / max(:39–41) / :42 / :43
    "| ×1, 127.74 MHz | +6.073e-11 | +5.444e-09 | +2.906e-09 | 2.373e-10 | 2.308e-10 | 5.648e-09 |",
    # 444 W7:51 / :52 / :53 / max(:39–41) / :42 / :43
    "| h = 0.01, 64.000 MHz | +6.437e-11 | +1.829e-09 | +6.507e-10 | 2.591e-10 | 2.499e-10 "
    "| 1.717e-09 |",
    # 448 W6:51 / :52 / :53 / max(:39–41) / :42 / :43
    "| h = 0.01, 127.74 MHz | +5.912e-11 | +1.027e-08 | +4.872e-09 | 2.408e-10 | 2.093e-10 "
    "| 4.356e-09 |",
    "",
    "Table B — the values side by side (Palace / ours):",
    "",
    "| Rung, frequency | Q | f₀/target − 1 | R_in [Ω] |",
    "|---|---|---|---|",
    # 421 W4:49 (.4f) / :48 + this file's ESR section (ours) / :50 (.3f)
    "| ×1, 64.000 MHz | 164.5434 / 164.5434 | -2.583046e-06 / -2.583111e-06 | 3175.979 / 3175.979 |",
    # 436 W6:52 / :51 / :53
    "| ×1, 127.74 MHz | 254.4693 / 254.4693 | -2.301744e-06 / -2.301805e-06 | 9790.704 / 9790.704 |",
    # 444 W7:52 / :51 / :53
    "| h = 0.01, 64.000 MHz | 164.2158 / 164.2158 | -2.588508e-06 / -2.588573e-06 "
    "| 3166.748 / 3166.748 |",
    # 448 W6:52 / :51 / :53
    "| h = 0.01, 127.74 MHz | 252.3012 / 252.3012 | -2.307258e-06 / -2.307317e-06 "
    "| 9698.113 / 9698.113 |",
    "",
    "Table C — each code's own ×1 → 0.01 move, X(h = 0.01) / X(×1) − 1 (Palace / ours):",
    "",
    "| Frequency | f₀ | Q | R_in |",
    "|---|---|---|---|",
    # 444 W7:54
    "| 64.000 MHz | -5.462608e-09 / -5.461672e-09 | -1.991402e-03 / -1.991409e-03 "
    "| -2.906654e-03 / -2.906656e-03 |",
    # 448 W6:54
    "| 127.74 MHz | -5.514397e-09 / -5.512787e-09 | -8.519899e-03 / -8.519904e-03 "
    "| -9.456957e-03 / -9.456959e-03 |",
    "",
    "The Q move is equal to five significant digits at 64 MHz and six at 127.74 MHz "
    "(−0.199 % / −0.852 %). The h-sequence of Q is therefore a property of the discretisation "
    "both codes share; no cross-solver reading shortens Q's h-statement, which is T41's, on ours "
    "alone.",
    "",
    "Not claimed: an h-statement (two rungs, one move each); Q, f₀ or R_in accuracy against any "
    "external value; the h = 0.005 Palace twin (declined, 2026-10-14 weekly); any AED value.",
)


def _write_comparison_palace():
    """§9 item 460: insert ``PALACE_SECTION`` into COMPARISON.md directly before the closing
    'Not claimed' line; an earlier copy of this section is replaced (idempotent)."""
    path = CASE_DIR / "COMPARISON.md"
    text = path.read_text()
    tail = "\n" + CLOSING_LINE + "\n"
    if PALACE_HEAD in text:
        i = text.index("\n" + PALACE_HEAD)
        j = text.index(tail, i)
        text = text[:i] + text[j:]
    assert text.endswith(tail), "COMPARISON.md closing line not found"
    k = len(text) - len(tail)
    path.write_text(text[:k] + "\n".join(PALACE_SECTION) + "\n" + text[k:])


def _x_stage_a(comm, outdir):
    """XA: Route A at degree 2, Cu, at C_F, three GRID_MHZ points bracketing f0."""
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans7_exact_strip import _exact_geo, _uj
    from tests.validation.test_ans7_sparse_route_b import DEFAULT_NODES_MHZ, bary_interp
    from tests.validation.test_port_circuit_layer_field import REDUCTION_BAND
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    rank0 = comm.rank == 0
    failures = []
    metrics_path = CASE_DIR / "metrics.json"
    text = metrics_path.read_text()
    m = json.loads(text)
    if "degree2_exact" not in m:
        raise RuntimeError("XA needs X4's degree2_exact block in metrics.json; run X4 first")
    d2 = m["degree2_exact"]
    recs = [json.loads(_x_node_path(outdir, mm).read_text()) for mm in DEFAULT_NODES_MHZ]
    z_cu = [_uj(r["z"]["cu"]) for r in recs]
    xn = np.asarray(DEFAULT_NODES_MHZ, dtype=float) - 64.5
    t0 = time.perf_counter()
    geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
    t_mesh = time.perf_counter() - t0
    if cells != d2["cells"]:
        failures.append(("cells", cells, d2["cells"]))
    _say0(comm, f"[ans7-x] XA: exact hole mesh {cells} cells (X4 {d2['cells']}), mesh {t_mesh:.1f} s")
    res, s_as, s_bs, z_in_a, t_drive = [], [], [], [], []
    for g in XA_MHZ:
        f = g * 1e6
        omega = 2.0 * np.pi * f
        t1 = time.perf_counter()
        s_a, _ = _route_a_solve(th14, geo, specs, f, "cu", degree=2)
        t_drive.append(time.perf_counter() - t1)
        s_b = _s11_from_z(route_b_z(bary_interp(xn, z_cu, g - 64.5), f, C_F)[1])
        r = _band_residual(s_a, s_b)
        res.append(r)
        s_as.append(complex(s_a))
        s_bs.append(complex(s_b))
        z_in_a.append(add_c1(_z_from_s11(s_a), omega))
        lo, hi = PRED_AB
        _say0(comm, f"[ans7-x] XA (AB) p2 Cu {g:.1f} MHz at C_F: S11_A {complex(s_a):.12e}, S11_B "
                    f"{complex(s_b):.12e}; |S11_A - S11_B|/|S11_B| = {r:.3e} (ASSERTED <= "
                    f"{REDUCTION_BAND:.0e}; PREDICTED [{lo:g}, {hi:g}]: "
                    f"{'met' if lo <= r <= hi else 'NOT met'}); drive {t_drive[-1]:.1f} s")
        if not r <= REDUCTION_BAND:
            failures.append(("AB", g, r))
    f_arr = np.array(XA_MHZ) * 1e6
    rows_a = parallel_rows(f_arr, np.array(z_in_a))
    b = (1.0 / np.asarray(z_in_a)).imag
    ups = [c for c in rows_a["crossings_up_hz"] if 62.6e6 < c < 62.8e6]
    _say0(comm, f"[ans7-x] XA Im Y_in,A at {XA_MHZ} MHz = " + ", ".join(f"{v:+.9e}" for v in b)
          + f" S; rows A: f0 {rows_a.get('f0_hz')!r} Hz, R_in {rows_a.get('r_in_ohm')!r} ohm, Q "
          f"{rows_a.get('q_loaded')!r}; - -> + crossings in (62.6, 62.8) MHz: "
          f"{[c / 1e6 for c in ups]} (ASSERTED exactly 1)")
    if len(ups) != 1 or len(rows_a["crossings_up_hz"]) != 1:
        failures.append(("crossing", rows_a["crossings_up_hz"]))
    rb = d2["rows"]["cu"]
    if "f0_hz" in rows_a:
        dfr = rows_a["f0_hz"] / rb["f0_hz"] - 1.0
        lo, hi = PRED_DF0_A_B
        dr = rows_a["r_in_ohm"] / rb["r_in_ohm"] - 1.0
        dq = rows_a["q_loaded"] / rb["q_loaded"] - 1.0
        ok = abs(dr) <= PRED_ROUTE_A_ROWS_REL and abs(dq) <= PRED_ROUTE_A_ROWS_REL
        _say0(comm, f"[ans7-x] XA f0,A - f0,B = {rows_a['f0_hz'] - rb['f0_hz']:+.3f} Hz, relative "
                    f"{dfr:+.3e} (PREDICTED [{lo:g}, {hi:g}]: {'met' if lo <= dfr <= hi else 'NOT met'}); "
                    f"R_in A/B - 1 = {dr:+.4%}, Q A/B - 1 = {dq:+.4%} (PREDICTED within "
                    f"+-{PRED_ROUTE_A_ROWS_REL:.0%}: {'met' if ok else 'NOT met'})")
    # Control (i), predicted, printed: circuit at 0.5·C_F vs in-model at C_F, 62.7 MHz.
    j = XA_MHZ.index(62.7)
    z627 = bary_interp(xn, z_cu, 62.7 - 64.5)
    s_half = _s11_from_z(route_b_z(z627, 62.7e6, CONTROL_I_C_FACTOR * C_F)[1])
    ctrl_i = _band_residual(s_as[j], s_half)
    _say0(comm, f"[ans7-x] XA control (i) (PREDICTED >= {PREDICTED_CONTROL_I_RATIO_MIN:g} x band): "
                f"circuit at {CONTROL_I_C_FACTOR}*C_F vs in-model at C_F, 62.7 MHz Cu: "
                f"{ctrl_i:.6e} = {ctrl_i / REDUCTION_BAND:.1f} x band -> "
                f"{'met' if ctrl_i >= PREDICTED_CONTROL_I_RATIO_MIN * REDUCTION_BAND else 'NOT met'}")
    if failures or not rank0:
        return failures
    snap = _snapshot(m, skip=("degree2_exact",))
    snap_d2 = _snapshot(d2, skip=("route_a",))
    d2["route_a"] = {"points_mhz": list(XA_MHZ), "degree": 2, "column": "cu", "c_f": C_F,
                     "s11_a": [[s.real, s.imag] for s in s_as],
                     "s11_b": [[s.real, s.imag] for s in s_bs], "residual_ab": res,
                     "reduction_band": float(REDUCTION_BAND), "rows": rows_a,
                     "control_i_half_c": ctrl_i, "t_drive_s": t_drive, "t_mesh_s": t_mesh}
    metrics_path.write_text(json.dumps(m, indent=1) + "\n")
    after = json.loads(metrics_path.read_text())
    bad = [k for k, v in snap.items() if json.dumps(after.get(k), sort_keys=True) != v]
    bad += [("degree2_exact", k) for k, v in snap_d2.items()
            if json.dumps(after["degree2_exact"].get(k), sort_keys=True) != v]
    _say0(comm, f"[ans7-x] XA (K): {len(snap)} top-level + {len(snap_d2)} degree2_exact keys "
                f"re-serialised identically: {'yes' if not bad else 'NO ' + str(bad)} (ASSERTED)")
    if bad:
        metrics_path.write_text(text)
        failures.append(("K", bad))
        return failures
    _write_comparison(after)
    _say0(comm, f"[ans7-x] XA: wrote degree2_exact['route_a'] into {metrics_path} and "
                f"regenerated COMPARISON.md")
    return failures


def _xt_stage(comm, stage):
    """XT64 / XT128 (§9 item 288, §10 T23): Route A at degree 2 at C_in, Cu, three tune
    points; re-asserts T21' / T22' through the example and records the rows.

    The field set is chosen by the stage name alone (no knob helper is called, so
    ``_x_guard`` stays satisfied).  C_in is read from the stored ``tune_p2.json``, never
    recomputed.  Every record, band and helper is imported lazily, never restated.
    """
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans7_exact_strip import (
        OUT_DIR, T21P_ZERO_RESOLUTION, T22P_ZERO_RESOLUTION, TUNE128_GRID_MHZ, TUNE128_MHZ,
        TUNE128_NODES_MHZ, TUNE128_OUT_DIR, TUNE128_TARGET_HZ, TUNE_MHZ, TUNE_TARGET_HZ,
        _assert_tie, _exact_geo, _node128_json, _stored_x1_p2, _tie_from_json, _uj)
    from tests.validation.test_ans7_sparse_route_b import (
        DEFAULT_NODES_MHZ, GRID_MHZ, NODE_S_RTOL, bary_interp, sparse_rows)
    from tests.validation.test_port_circuit_layer_field import REDUCTION_BAND
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    rank0 = comm.rank == 0
    failures = []
    is64 = stage == "XT64"
    key = "64.0" if is64 else "127.74"
    pt, tp, tag = ("T21'", "T2'", "AB'") if is64 else ("T22'", "T2''", "AB''")
    if is64:
        target_hz, tune, nodes, grid, out = TUNE_TARGET_HZ, TUNE_MHZ, DEFAULT_NODES_MHZ, GRID_MHZ, OUT_DIR
        band, band_name = T21P_ZERO_RESOLUTION, "T21P_ZERO_RESOLUTION"
        paths = [OUT_DIR / f"node_{m:.1f}_p2.json" for m in nodes]
        pred = {"c": -1.494697e-5, "rt": -2.535846e-6, "t2": -2.535844e-6,
                "ab": (1.237e-5, 1.258e-5, 1.287e-5), "r_in": 18739.27, "q": 971.07}
    else:
        target_hz, tune, nodes, grid, out = (TUNE128_TARGET_HZ, TUNE128_MHZ, TUNE128_NODES_MHZ,
                                             TUNE128_GRID_MHZ, TUNE128_OUT_DIR)
        band, band_name = T22P_ZERO_RESOLUTION, "T22P_ZERO_RESOLUTION"
        paths = [_node128_json(m, ("cu",)) for m in nodes]
        pred = {"c": -7.754048e-6, "rt": -2.264620e-6, "t2": -2.264621e-6,
                "ab": (1.590e-5, 1.601e-5, 1.616e-5), "r_in": 23793.40, "q": 618.34}
    stored = {n: out / n for n in ("tune_p2.json", "tune_inmodel_p2.json",
                                   "tune_inmodel_cin_p2.json")}
    missing = [str(p) for p in list(stored.values()) + paths if not p.exists()]
    if missing:
        raise RuntimeError(f"{stage}: stored inputs absent (never re-created here): {missing}")
    tj = json.loads(stored["tune_p2.json"].read_text())
    st_cm = json.loads(stored["tune_inmodel_p2.json"].read_text())
    st_ci = json.loads(stored["tune_inmodel_cin_p2.json"].read_text())
    if tj["degree"] != 2 or abs(tj["f_hz"] - target_hz) > 1.0:
        raise RuntimeError(f"{stage}: tune_p2.json degree {tj['degree']} f_hz {tj['f_hz']}")
    c_in, c_mode = float(tj["c_in"]), float(tj["c_mode"]["cu"])
    if is64:
        recs = _stored_x1_p2()
    else:
        recs = [json.loads(p.read_text()) for p in paths]
        bad = [(r["f_hz"], r["degree"], r["cells"], r["kappa"], r.get("resolution"),
                r.get("target_mhz")) for r in recs
               if r["degree"] != 2 or r["kappa"] is not None or r.get("resolution") is not None
               or r.get("target_mhz") != 127.74]
        if bad or len({r["cells"] for r in recs}) != 1:
            raise RuntimeError(f"{stage}: stored node JSONs inconsistent: {bad}")
        for m, r in zip(nodes, recs):
            assert abs(r["f_hz"] - m * 1e6) < 1.0, (m, r["f_hz"])
            _assert_tie(r["cells"], _tie_from_json(r["tie"]))
    z_cu = [_uj(r["z"]["cu"]) for r in recs]
    xn = np.asarray(nodes, dtype=float) - 64.5
    if st_ci["c_in"] != tj["c_in"] or st_ci["degree"] != 2:
        failures.append(("stored-cin-c_in/degree", st_ci["c_in"], tj["c_in"], st_ci["degree"]))
    if st_ci["ranks"] != 8:
        failures.append(("stored-cin-ranks", st_ci["ranks"]))
    _say0(comm, f"[ans7-xt] {stage}: C_in {c_in!r} F and C_mode,Cu {c_mode!r} F read from the stored "
                f"tune_p2.json; stored nodes {nodes}; -n {comm.size}")

    t0 = time.perf_counter()
    geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
    t_mesh = time.perf_counter() - t0
    if cells != st_ci["cells"]:
        failures.append(("cells", cells, st_ci["cells"]))
    _say0(comm, f"[ans7-xt] {stage}: exact hole mesh {cells} cells (stored {st_ci['cells']}), "
                f"mesh {t_mesh:.1f} s")
    res, s_as, s_bs, z_in_a, t_drive = [], [], [], [], []
    for i, g in enumerate(tune):
        f = g * 1e6
        omega = 2.0 * np.pi * f
        t1 = time.perf_counter()
        s_a, _ = _route_a_solve(th14, geo, specs, f, "cu", c_f=c_in, degree=2)
        t_drive.append(time.perf_counter() - t1)
        s_b = _s11_from_z(route_b_z(bary_interp(xn, z_cu, g - 64.5), f, c_in)[1])
        r = _band_residual(s_a, s_b)
        res.append(r)
        s_as.append(complex(s_a))
        s_bs.append(complex(s_b))
        z_in_a.append(add_c1(_z_from_s11(s_a), omega, c_in))
        _say0(comm, f"[ans7-xt] {pt} ({tag}) p2 Cu {g:.2f} MHz at C_in: S11_A {complex(s_a):.12e}, "
                    f"S11_B {complex(s_b):.12e}; |S11_A - S11_B|/|S11_B| = {r:.3e} (ASSERTED <= "
                    f"{REDUCTION_BAND:.0e}; PREDICTED {pred['ab'][i]:.3e}); drive {t_drive[-1]:.1f} s")
        if not r <= REDUCTION_BAND:
            failures.append((tag, g, r))
    rows_a = parallel_rows(np.array(tune) * 1e6, np.array(z_in_a))
    b = (1.0 / np.asarray(z_in_a)).imag
    _say0(comm, f"[ans7-xt] {pt} Im Y_in,A at {list(tune)} MHz = " + ", ".join(f"{v:+.9e}" for v in b)
          + f" S; rows A: f0 {rows_a.get('f0_hz')!r} Hz, R_in {rows_a.get('r_in_ohm')!r} ohm, Q "
          f"{rows_a.get('q_loaded')!r}; up-crossings {rows_a['crossings_up_hz']}")
    f0a = rows_a.get("f0_hz")
    lo_hz, hi_hz = tune[0] * 1e6, tune[-1] * 1e6
    sa = (f0a / target_hz - 1.0) if f0a else float("nan")
    t2_ok = (len(rows_a["crossings_up_hz"]) == 1 and f0a is not None and lo_hz < f0a < hi_hz
             and abs(sa) <= band)
    _say0(comm, f"[ans7-xt] {pt} ({tp}) f0,A(C_in) {f0a!r} Hz; f0,A/target - 1 = {sa:+.6e} (ASSERTED "
                f"exactly one - -> + crossing in ({tune[0]}, {tune[-1]}) MHz and |.| <= {band_name} "
                f"{band:.3e}: {'met' if t2_ok else 'NOT met'}; = {abs(sa) / band:.3f} bands; "
                f"PREDICTED {pred['t2']:+.6e})")
    if not t2_ok:
        failures.append((tp, rows_a["crossings_up_hz"], sa))
    zc = {"cu": z_cu}
    c_rel = c_mode / c_in - 1.0
    rb_ci = sparse_rows(nodes, zc, grid_mhz=grid, columns=("cu",), c_f=c_in)["cu"]["rows"]
    f0b_ci = rb_ci.get("f0_hz")
    r_ab = (f0a / f0b_ci - 1.0) if (f0a and f0b_ci) else float("nan")
    _say0(comm, f"[ans7-xt] {pt} c = C_mode/C_in - 1 = {c_rel:+.6e} (PREDICTED {pred['c']:+.6e}); "
                f"f0,B(C_in) {f0b_ci!r} Hz; r = f0,A/f0,B(C_in) - 1 = {r_ab:+.6e} (PREDICTED "
                f"{pred['rt']:+.6e}); R_in {rows_a.get('r_in_ohm')!r} ohm (PREDICTED {pred['r_in']}), "
                f"Q {rows_a.get('q_loaded')!r} (PREDICTED {pred['q']}); records, not gates")
    # (R-T): reproduction of the registered run, against its stored tune_inmodel_cin_p2.json.
    st_a = st_ci["rows_a"]
    d_f0 = abs(f0a / st_a["f0_hz"] - 1.0) if (f0a and "f0_hz" in st_a) else float("inf")
    d_res = max(abs(a - s) for a, s in zip(res, st_ci["residual_ab"]))
    d_s = max(abs(complex(a) - complex(*s)) / abs(complex(*s)) for a, s in zip(s_as, st_ci["s11_a"]))
    rt = {"f0_rel": d_f0, "residual_abs": d_res, "s11_a_rel": d_s}
    rt_assert = comm.size == 8
    rt_ok = max(d_f0, d_res, d_s) <= NODE_S_RTOL
    _say0(comm, f"[ans7-xt] {pt} (R-T) vs stored tune_inmodel_cin_p2.json: |f0 ratio - 1| {d_f0:.3e}, "
                f"max |res diff| {d_res:.3e}, max S11_A rel {d_s:.3e} ("
                f"{'ASSERTED' if rt_assert else 'printed only, width != 8'} <= {NODE_S_RTOL:.0e}: "
                f"{'met' if rt_ok else 'NOT met'}; PREDICTED <= 1e-12)")
    if rt_assert and not rt_ok:
        failures.append(("R-T", rt))
    # (N-C): the (R-T) comparator separates C_in from C_mode.
    f0_cm = st_cm["rows_a"]["f0_hz"]
    n_c = abs(f0a / f0_cm - 1.0) if f0a else float("nan")
    _say0(comm, f"[ans7-xt] {pt} (N-C) |f0,A(C_in)/f0,A(C_mode, stored) - 1| = {n_c:.6e} = "
                f"{n_c / NODE_S_RTOL:.1f} x NODE_S_RTOL (ASSERTED > 100 x: "
                f"{'met' if n_c > 100.0 * NODE_S_RTOL else 'NOT met'})")
    if not n_c > 100.0 * NODE_S_RTOL:
        failures.append(("N-C", n_c))
    # (N1') at 64 MHz / (N0') at 127.74: the untuned (C_F) reading on the stored nodes.
    n_cf = sparse_rows(nodes, zc, grid_mhz=grid, columns=("cu",), c_f=C_F)["cu"]["rows"]
    n1_off = n0x = None
    if is64:
        n1_off = abs(n_cf["f0_hz"] / target_hz - 1.0) if "f0_hz" in n_cf else float("nan")
        n1_ok = n1_off >= 10.0 * T21P_ZERO_RESOLUTION
        _say0(comm, f"[ans7-xt] {pt} (N1') |f0,B(C_F)/64 MHz - 1| = {n1_off:.6e} = "
                    f"{n1_off / T21P_ZERO_RESOLUTION:.1f} x T21P_ZERO_RESOLUTION (ASSERTED >= 10 x: "
                    f"{'met' if n1_ok else 'NOT met'}; backing 2.029559e-2)")
        if not n1_ok:
            failures.append(("N1'", n1_off))
    else:
        n0x = list(n_cf["crossings_up_hz"]) if "crossings_up_hz" in n_cf else None
        _say0(comm, f"[ans7-xt] {pt} (N0') at C_F on the stored 127.74 MHz Cu grid: up-crossings {n0x} "
                    f"(ASSERTED == [])")
        if n0x != []:
            failures.append(("N0'", n0x))
    # (N2'), predicted, printed: circuit at 0.5*C_in vs in-model at C_in at the target.
    j = [round(g, 6) for g in tune].index(round(target_hz / 1e6, 6))
    z_t = bary_interp(xn, z_cu, target_hz / 1e6 - 64.5)
    s_n2 = _s11_from_z(route_b_z(z_t, target_hz, 0.5 * c_in)[1])
    n2 = float(abs(s_n2 - s_as[j]) / abs(s_as[j]))
    _say0(comm, f"[ans7-xt] {pt} (N2{chr(8242) if is64 else chr(8243)}) circuit at 0.5*C_in vs in-model at "
                f"C_in, target Cu: |dS|/|S| = {n2:.6e} = {n2 / REDUCTION_BAND:.2f} x REDUCTION_BAND "
                f"(PREDICTED >= 10 x: {'met' if n2 >= 10.0 * REDUCTION_BAND else 'NOT met'}; "
                f"printed, never asserted)")
    _say0(comm, f"[ans7-xt] {stage}: total {time.perf_counter() - t0:.1f} s (mesh {t_mesh:.1f} s, drives "
                f"{', '.join(f'{t:.1f}' for t in t_drive)} s)")
    if failures or not rank0:
        return failures
    metrics_path = CASE_DIR / "metrics.json"
    text = metrics_path.read_text()
    m = json.loads(text)
    if "degree2_exact" not in m:
        raise RuntimeError("XT needs X4's degree2_exact block in metrics.json")
    d2 = m["degree2_exact"]
    tl0 = d2.get("tuned_larmor", {})
    snap = _snapshot(m, skip=("degree2_exact",))
    snap_d2 = _snapshot(d2, skip=("tuned_larmor",))
    snap_other = _snapshot(tl0, skip=(key,))
    d2.setdefault("tuned_larmor", {})[key] = {
        "target_hz": target_hz, "c_in_f": c_in, "c_mode_f": c_mode, "c_rel": c_rel,
        "tune_mhz": list(tune), "s11_a": [[s.real, s.imag] for s in s_as],
        "s11_b": [[s.real, s.imag] for s in s_bs], "residual_ab": res,
        "reduction_band": float(REDUCTION_BAND), "rows_a": rows_a, "f0a_rel_target": sa,
        "zero_resolution": {"value": band, "constant": band_name},
        "bands_fraction": abs(sa) / band, "f0b_cin_hz": f0b_ci, "r": r_ab,
        "n1_off": n1_off, "n0_crossings_up_hz": n0x, "n2": n2, "rt": rt,
        "cells": cells, "mpi_ranks": comm.size, "t_mesh_s": t_mesh, "t_drive_s": t_drive,
        "source": ["§9 item 262 (ANS-7-item262 logs)" if is64 else "§9 item 283 (ANS-7-item283 logs)",
                   "stored tune_p2.json / tune_inmodel_p2.json / tune_inmodel_cin_p2.json"]}
    metrics_path.write_text(json.dumps(m, indent=1) + "\n")
    after = json.loads(metrics_path.read_text())
    bad = [k for k, v in snap.items() if json.dumps(after.get(k), sort_keys=True) != v]
    bad += [("degree2_exact", k) for k, v in snap_d2.items()
            if json.dumps(after["degree2_exact"].get(k), sort_keys=True) != v]
    bad += [("tuned_larmor", k) for k, v in snap_other.items()
            if json.dumps(after["degree2_exact"]["tuned_larmor"].get(k), sort_keys=True) != v]
    _say0(comm, f"[ans7-xt] {stage} (K): {len(snap)} top-level + {len(snap_d2)} degree2_exact + "
                f"{len(snap_other)} other-field keys re-serialised identically: "
                f"{'yes' if not bad else 'NO ' + str(bad)} (ASSERTED)")
    if bad:
        metrics_path.write_text(text)
        failures.append(("K", bad))
        return failures
    _write_comparison(after)
    _say0(comm, f"[ans7-xt] {stage}: wrote degree2_exact['tuned_larmor']['{key}'] into {metrics_path} "
                f"and regenerated COMPARISON.md")
    return failures


def _bits(x) -> bytes:
    return np.ascontiguousarray(np.asarray(x, dtype=np.complex128)).tobytes()


def _xt_esr_stage(comm):
    """XT64R (§9 item 378, §10 T27): item 311's T26 (r = 0.2 ohm ESR on every capacitor, both
    sides) re-driven through the example at (r, C_in(r)), 63.9 / 64.0 / 64.1 MHz.

    The field set is chosen by the stage name alone (``_x_guard`` stays satisfied).  C_in(r)
    and r are read from the gitignored ``tune_p2_esr0.2.json`` stage file, never recomputed;
    the ESR terminations are the keyworded circuit layer above.  f0,B is formed as item 311
    forms it (``tune_esr`` in ``test_ans7_exact_strip.py``): ``parallel_rows`` over the three
    tune points of Route B's Z_in.  Every record, band and helper is imported lazily.
    """
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans7_exact_strip import (
        ESR_VALUES, OUT_DIR, T21P_ZERO_RESOLUTION, TUNE_MHZ, TUNE_TARGET_HZ, _exact_geo,
        _stored_x1_p2, _uj)
    from tests.validation.test_ans7_sparse_route_b import (
        DEFAULT_NODES_MHZ, NODE_S_RTOL, bary_interp)
    from tests.validation.test_port_circuit_layer_field import REDUCTION_BAND
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    stage = "XT64R"
    rank0 = comm.rank == 0
    failures = []
    key = "64.0"
    target_hz, tune, nodes = TUNE_TARGET_HZ, TUNE_MHZ, DEFAULT_NODES_MHZ
    band, band_name = T21P_ZERO_RESOLUTION, "T21P_ZERO_RESOLUTION"
    names = {"tp": "tune_p2.json", "te": "tune_p2_esr0.2.json",
             "st": "tune_inmodel_cin_esr0.2_p2.json", "ci0": "tune_inmodel_cin_p2.json"}
    paths = {k: OUT_DIR / v for k, v in names.items()}
    missing = [str(p) for p in paths.values() if not p.exists()]
    if missing:
        raise RuntimeError(f"{stage}: stored inputs absent (never re-created here): {missing}")
    tj = json.loads(paths["tp"].read_text())
    te = json.loads(paths["te"].read_text())
    st = json.loads(paths["st"].read_text())
    ci0 = json.loads(paths["ci0"].read_text())
    r_esr = float(te["esr_ohm"])
    if not (r_esr == float(ESR_VALUES[0]) == 0.2) or abs(te["f_hz"] - target_hz) > 1.0:
        raise RuntimeError(f"{stage}: esr_ohm {te['esr_ohm']!r} / f_hz {te['f_hz']!r} / "
                           f"ESR_VALUES {ESR_VALUES}")
    c_in_r, c_in_0 = float(te["c_in_r"]), float(te["c_in_0"])
    if c_in_0 != float(tj["c_in"]) or tj["degree"] != 2 or st["c_in_r"] != c_in_r:
        failures.append(("stored-inputs", c_in_0, tj["c_in"], tj["degree"], st["c_in_r"]))
    if st["ranks"] != 8 or st["degree"] != 2 or st["esr_ohm"] != r_esr:
        failures.append(("stored-esr-ranks/degree/esr", st["ranks"], st["degree"], st["esr_ohm"]))
    z64 = _uj(tj["z64"]["cu"])
    recs = _stored_x1_p2()
    z_cu = [_uj(r["z"]["cu"]) for r in recs]
    xn = np.asarray(nodes, dtype=float) - 64.5
    _say0(comm, f"[ans7-xtr] {stage}: r = {r_esr!r} ohm (== ESR_VALUES[0] {ESR_VALUES[0]}), C_in(r) "
                f"{c_in_r!r} F, C_in(0) {c_in_0!r} F read from the stored tune_p2_esr0.2.json "
                f"(item 311); -n {comm.size}")

    # (R0): the r = 0 branch is bitwise, on the stored x1 z64 at 64.0 MHz.
    f64 = target_hz
    om64 = 2.0 * np.pi * f64
    old_cap = 1.0 / (1j * om64 * c_in_0)
    zb_u, zb_k = z_cap(om64, c_in_0), z_cap(om64, c_in_0, esr_ohm=0.0)
    ac_u, ac_k = add_c1(z64[0, 0], om64, c_in_0), add_c1(z64[0, 0], om64, c_in_0, esr_ohm=0.0)
    rb_u, rb_k = route_b_z(z64, f64, c_in_0), route_b_z(z64, f64, c_in_0, esr_ohm=0.0)
    r0 = {"z_cap": _bits(zb_u) == _bits(zb_k) == _bits(old_cap),
          "gap_z": _bits(_gap_z(f64, c_in_0)) == _bits(old_cap)
          == _bits(_gap_z(f64, c_in_0, 0.0)),
          "add_c1": _bits(ac_u) == _bits(ac_k),
          "route_b_z": all(_bits(a) == _bits(b) for a, b in zip(rb_u, rb_k))}
    r0_ok = all(r0.values())
    _say0(comm, f"[ans7-xtr] {stage} (R0) the r = 0 branch is bitwise (z_cap, the gap helper, add_c1, "
                f"route_b_z with esr_ohm=0.0 vs the unkeyworded call; z_cap/helper vs the literal old "
                f"1/(jwC)): {r0} (ASSERTED all True: {'met' if r0_ok else 'NOT met'})")
    if not r0_ok:
        failures.append(("R0", r0))
        return failures

    t0 = time.perf_counter()
    geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
    t_mesh = time.perf_counter() - t0
    if cells != st["cells"]:
        failures.append(("cells", cells, st["cells"]))
    _say0(comm, f"[ans7-xtr] {stage}: exact hole mesh {cells} cells (stored {st['cells']}), "
                f"mesh {t_mesh:.1f} s")
    res, s_as, s_bs, z_in_a, z_in_b, t_drive = [], [], [], [], [], []
    for i, g in enumerate(tune):
        f = g * 1e6
        omega = 2.0 * np.pi * f
        t1 = time.perf_counter()
        s_a, _ = _route_a_solve(th14, geo, specs, f, "cu", c_f=c_in_r, degree=2, esr_ohm=r_esr)
        t_drive.append(time.perf_counter() - t1)
        zin_b, zbare_b, _g = route_b_z(bary_interp(xn, z_cu, g - 64.5), f, c_in_r, esr_ohm=r_esr)
        s_b = _s11_from_z(zbare_b)
        r = _band_residual(s_a, s_b)
        res.append(r)
        s_as.append(complex(s_a))
        s_bs.append(complex(s_b))
        z_in_a.append(add_c1(_z_from_s11(s_a), omega, c_in_r, esr_ohm=r_esr))
        z_in_b.append(complex(zin_b))
        _say0(comm, f"[ans7-xtr] T26 (AB-r) p2 Cu {g:.1f} MHz at (r, C_in(r)): S11_A "
                    f"{complex(s_a):.12e}, S11_B {complex(s_b):.12e}; |S11_A - S11_B|/|S11_B| = "
                    f"{r:.3e} (ASSERTED <= REDUCTION_BAND {REDUCTION_BAND:.0e}; PREDICTED "
                    f"{st['residual_ab'][i]:.3e}, item 311); drive {t_drive[-1]:.1f} s")
        if not r <= REDUCTION_BAND:
            failures.append(("AB-r", g, r))
    rows_a = parallel_rows(np.array(tune) * 1e6, np.array(z_in_a))
    rows_b = parallel_rows(np.array(tune) * 1e6, np.array(z_in_b))
    b = (1.0 / np.asarray(z_in_a)).imag
    _say0(comm, f"[ans7-xtr] T26 Im Y_in,A at {list(tune)} MHz = " + ", ".join(f"{v:+.9e}" for v in b)
          + f" S; rows A: f0 {rows_a.get('f0_hz')!r} Hz, R_in {rows_a.get('r_in_ohm')!r} ohm, Q "
          f"{rows_a.get('q_loaded')!r}; up-crossings {rows_a['crossings_up_hz']}")
    f0a, f0b = rows_a.get("f0_hz"), rows_b.get("f0_hz")
    sa = (f0a / target_hz - 1.0) if f0a else float("nan")
    t2_ok = (len(rows_a["crossings_up_hz"]) == 1 and f0a is not None
             and tune[0] * 1e6 < f0a < tune[-1] * 1e6 and abs(sa) <= band)
    _say0(comm, f"[ans7-xtr] T26 (T2-r) f0,A(C_in(r)) {f0a!r} Hz; f0,A/64 MHz - 1 = {sa:+.6e} (ASSERTED "
                f"exactly one - -> + crossing in ({tune[0]}, {tune[-1]}) MHz and |.| <= {band_name} "
                f"{band:.3e}: {'met' if t2_ok else 'NOT met'}; = {abs(sa) / band:.3f} bands; "
                f"PREDICTED -2.583111e-6, item 311)")
    if not t2_ok:
        failures.append(("T2-r", rows_a["crossings_up_hz"], sa))
    sb = (f0b / target_hz - 1.0) if f0b else float("nan")
    r_ab = (f0a / f0b - 1.0) if (f0a and f0b) else float("nan")
    _say0(comm, f"[ans7-xtr] T26 f0,B(C_in(r)) {f0b!r} Hz (Route B Z_in on the three tune points, as "
                f"item 311 forms it); f0,B/64 MHz - 1 = {sb:+.6e} (PREDICTED within +-1e-9); "
                f"r_AB = f0,A/f0,B - 1 = {r_ab:+.6e} (PREDICTED -2.583113e-6, item 311); records, "
                f"not gates")
    # (R-T-r): reproduction of item 311's registered run, against its stored JSON.
    st_a = st["rows_a"]
    d_f0 = abs(f0a / st_a["f0_hz"] - 1.0) if (f0a and "f0_hz" in st_a) else float("inf")
    d_res = max(abs(a - s) for a, s in zip(res, st["residual_ab"]))
    d_s = max(abs(complex(a) - complex(*s)) / abs(complex(*s)) for a, s in zip(s_as, st["s11_a"]))
    rt = {"f0_rel": d_f0, "residual_abs": d_res, "s11_a_rel": d_s}
    rt_assert = comm.size == 8
    rt_ok = max(d_f0, d_res, d_s) <= NODE_S_RTOL
    _say0(comm, f"[ans7-xtr] T26 (R-T-r) vs stored tune_inmodel_cin_esr0.2_p2.json: |f0 ratio - 1| "
                f"{d_f0:.3e}, max |res diff| {d_res:.3e}, max S11_A rel {d_s:.3e} ("
                f"{'ASSERTED' if rt_assert else 'printed only, width != 8'} <= {NODE_S_RTOL:.0e}: "
                f"{'met' if rt_ok else 'NOT met'}; PREDICTED <= 1e-12)")
    if rt_assert and not rt_ok:
        failures.append(("R-T-r", rt))
    # (N-C-r): the (R-T-r) comparator separates the ESR run from the lossless run -- asserted
    # only if item 311's own f0,A(r) and f0,A(0) differ by > 1e-7 (its stored rows).
    f0_ll = ci0["rows_a"]["f0_hz"]
    sep311 = abs(st_a["f0_hz"] / f0_ll - 1.0)
    n_c = abs(f0a / f0_ll - 1.0) if f0a else float("nan")
    n_c_assert = sep311 > 1.0e-7
    n_c_ok = n_c > 100.0 * NODE_S_RTOL
    _say0(comm, f"[ans7-xtr] T26 (N-C-r) |f0,A,ex(r)/f0,A,stored(lossless C_in) - 1| = {n_c:.6e} = "
                f"{n_c / NODE_S_RTOL:.1f} x NODE_S_RTOL (item 311's stored separation {sep311:.6e}; "
                f"{'ASSERTED > 100 x' if n_c_assert else 'PRINTED only, 311 separation <= 1e-7'}: "
                f"{'met' if n_c_ok else 'NOT met'}; PREDICTED O(1/Q_C^2) ~ 3e-5, 311 measured 4.7e-8)")
    if n_c_assert and not n_c_ok:
        failures.append(("N-C-r", n_c))
    # (N-L), printed (review ruling 2026-10-11 18:00): lossless circuit at C_in(0) vs in-model.
    j = [round(g, 6) for g in tune].index(round(target_hz / 1e6, 6))
    s_nl = _s11_from_z(route_b_z(z64, target_hz, c_in_0)[1])
    n_l = float(abs(s_nl - s_as[j]) / abs(s_as[j]))
    _say0(comm, f"[ans7-xtr] T26 (N-L) lossless circuit at C_in(0) vs in-model at (r, C_in(r)), 64.0 MHz "
                f"Cu: |dS11_bare|/|S11| = {n_l:.6e} = {n_l / REDUCTION_BAND:.2f} x REDUCTION_BAND "
                f"(PREDICTED >= 10 x: {'met' if n_l >= 10.0 * REDUCTION_BAND else 'NOT met'}; item "
                f"311 measured {st['n_l'] / REDUCTION_BAND:.2f} x, not met; PRINTED, never asserted)")
    # Q(r)/Q_0 beside the closed form; C_in(r)/C_in(0) - 1 (printed, never asserted).
    d2_in = json.loads((CASE_DIR / "metrics.json").read_text()).get("degree2_exact", {})
    tl64 = d2_in.get("tuned_larmor", {}).get(key)
    if tl64 is None:
        raise RuntimeError(f"{stage}: metrics.json degree2_exact['tuned_larmor']['{key}'] absent")
    q0, r_in0 = tl64["rows_a"]["q_loaded"], tl64["rows_a"]["r_in_ohm"]
    qr = rows_a.get("q_loaded")
    q_c = 1.0 / (2.0 * np.pi * target_hz * c_in_r * r_esr)
    q_cf = 1.0 / (1.0 / q0 + 1.0 / q_c)
    q_ratio = qr / q0 if qr else float("nan")
    c_rel = c_in_r / c_in_0 - 1.0
    _say0(comm, f"[ans7-xtr] T26 Q(r) {qr!r} vs Q_0 = lossless tuned_larmor['64.0'] rows_a.q_loaded "
                f"{q0!r}: Q(r)/Q_0 = {q_ratio:.6f}; closed form Q_C = 1/(omega C_in r) = {q_c:.4f}, "
                f"1/Q = 1/Q_0 + 1/Q_C -> Q {q_cf:.4f}, ratio {q_cf / q0:.6f}; Q(r)/closed form - 1 = "
                f"{(qr / q_cf - 1.0) if qr else float('nan'):+.4%} (PREDICTED ratio 0.169446 vs "
                f"0.167603, item 311; printed, never asserted); R_in(r)/R_in_0 = "
                f"{rows_a.get('r_in_ohm', float('nan')) / r_in0:.6f}")
    _say0(comm, f"[ans7-xtr] T26 C_in(r)/C_in(0) - 1 = {c_rel:+.6e} (PREDICTED +4.602610e-5, item 311; "
                f"printed, never asserted)")
    _say0(comm, f"[ans7-xtr] {stage}: total {time.perf_counter() - t0:.1f} s (mesh {t_mesh:.1f} s, drives "
                f"{', '.join(f'{t:.1f}' for t in t_drive)} s)")
    if failures or not rank0:
        return failures
    metrics_path = CASE_DIR / "metrics.json"
    text = metrics_path.read_text()
    m = json.loads(text)
    d2 = m["degree2_exact"]
    snap = _snapshot(m, skip=("degree2_exact",))
    snap_d2 = _snapshot(d2, skip=("tuned_larmor_esr",))
    snap_tl = _snapshot(d2.get("tuned_larmor", {}))
    d2.setdefault("tuned_larmor_esr", {})[key] = {
        "esr_ohm": r_esr, "target_hz": target_hz, "c_in_r_f": c_in_r, "c_in_0_f": c_in_0,
        "c_in_rel": c_rel, "tune_mhz": list(tune),
        "s11_a": [[s.real, s.imag] for s in s_as], "s11_b": [[s.real, s.imag] for s in s_bs],
        "residual_ab": res, "reduction_band": float(REDUCTION_BAND), "rows_a": rows_a,
        "rows_b": rows_b, "f0a_rel_target": sa,
        "zero_resolution": {"value": band, "constant": band_name},
        "bands_fraction": abs(sa) / band, "r_ab": r_ab, "q_ratio": q_ratio,
        "q_closed_form": {"q_c": q_c, "q": q_cf, "ratio": q_cf / q0, "q_0": q0},
        "lossless_comparand": "tuned_larmor['64.0']", "rt": rt, "n_l": n_l,
        "n_c_r": n_c, "cells": cells, "mpi_ranks": comm.size, "t_mesh_s": t_mesh,
        "t_drive_s": t_drive,
        "source": ["§9 item 311 (ANS-7-item311-W1-cin64-lossless / -W2-cin64-esr logs)",
                   "stored tune_p2.json / tune_p2_esr0.2.json / tune_inmodel_cin_esr0.2_p2.json"]}
    metrics_path.write_text(json.dumps(m, indent=1) + "\n")
    after = json.loads(metrics_path.read_text())
    bad = [k for k, v in snap.items() if json.dumps(after.get(k), sort_keys=True) != v]
    bad += [("degree2_exact", k) for k, v in snap_d2.items()
            if json.dumps(after["degree2_exact"].get(k), sort_keys=True) != v]
    bad += [("tuned_larmor", k) for k, v in snap_tl.items()
            if json.dumps(after["degree2_exact"]["tuned_larmor"].get(k), sort_keys=True) != v]
    _say0(comm, f"[ans7-xtr] {stage} (K): {len(snap)} top-level + {len(snap_d2)} degree2_exact + "
                f"{len(snap_tl)} tuned_larmor keys re-serialised identically: "
                f"{'yes' if not bad else 'NO ' + str(bad)} (ASSERTED)")
    if bad:
        metrics_path.write_text(text)
        failures.append(("K", bad))
        return failures
    _write_comparison(after)
    _say0(comm, f"[ans7-xtr] {stage}: wrote degree2_exact['tuned_larmor_esr']['{key}'] into "
                f"{metrics_path} and regenerated COMPARISON.md")
    return failures


def _xt_esr128_stage(comm):
    """XT128R (§9 item 391, §10 T36): item 387's T35 (r = 0.2 ohm ESR on every capacitor, both
    sides) re-driven through the example at (r, C_in(r)), 127.64 / 127.74 / 127.84 MHz.

    ``XT64R``'s function is left byte-unchanged; this is its 127.74 MHz sibling, reached only by
    the stage name (``_x_guard`` stays satisfied).  r, C_in(r) are read from the gitignored
    ``logs/ans7_exact/t127.74/`` stage files, never recomputed.  f0,B is formed as item 378
    forms it at 64 MHz (item 311's ``parallel_rows`` over the three tune points of Route B's
    Z_in).  ``64.5`` is the interpolation origin at both fields (as in ``XT128``).  Every
    record, band and helper is imported lazily, never restated.
    """
    from tests.validation import test_th14_birdcage_copper as th14
    from tests.validation.test_ans7_exact_strip import (
        ESR_VALUES, T22P_ZERO_RESOLUTION, TUNE128_MHZ, TUNE128_NODES_MHZ, TUNE128_OUT_DIR,
        TUNE128_TARGET_HZ, _assert_tie, _exact_geo, _node128_json, _tie_from_json, _uj)
    from tests.validation.test_ans7_sparse_route_b import NODE_S_RTOL, bary_interp
    from tests.validation.test_port_circuit_layer_field import REDUCTION_BAND
    from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

    stage = "XT128R"
    rank0 = comm.rank == 0
    failures = []
    key = "127.74"
    target_hz, tune, nodes = TUNE128_TARGET_HZ, TUNE128_MHZ, TUNE128_NODES_MHZ
    band, band_name = T22P_ZERO_RESOLUTION, "T22P_ZERO_RESOLUTION"
    out = TUNE128_OUT_DIR
    names = {"tp": "tune_p2.json", "te": "tune_p2_esr0.2.json",
             "st": "tune_inmodel_cin_esr0.2_p2.json", "ci0": "tune_inmodel_cin_p2.json"}
    paths = {k: out / v for k, v in names.items()}
    node_paths = [_node128_json(m, ("cu",)) for m in nodes]
    missing = [str(p) for p in list(paths.values()) + node_paths if not p.exists()]
    if missing:
        raise RuntimeError(f"{stage}: stored inputs absent (never re-created here): {missing}")
    tj = json.loads(paths["tp"].read_text())
    te = json.loads(paths["te"].read_text())
    st = json.loads(paths["st"].read_text())
    ci0 = json.loads(paths["ci0"].read_text())
    r_esr = float(te["esr_ohm"])
    if not (r_esr == float(ESR_VALUES[0]) == 0.2) or abs(te["f_hz"] - target_hz) > 1.0:
        raise RuntimeError(f"{stage}: esr_ohm {te['esr_ohm']!r} / f_hz {te['f_hz']!r} / "
                           f"ESR_VALUES {ESR_VALUES}")
    c_in_r, c_in_0 = float(te["c_in_r"]), float(te["c_in_0"])
    if c_in_0 != float(tj["c_in"]) or tj["degree"] != 2 or st["c_in_r"] != c_in_r:
        failures.append(("stored-inputs", c_in_0, tj["c_in"], tj["degree"], st["c_in_r"]))
    if st["ranks"] != 8 or st["degree"] != 2 or st["esr_ohm"] != r_esr:
        failures.append(("stored-esr-ranks/degree/esr", st["ranks"], st["degree"], st["esr_ohm"]))
    z64 = _uj(tj["z64"]["cu"])
    recs = [json.loads(p.read_text()) for p in node_paths]
    bad = [(r["f_hz"], r["degree"], r["cells"], r["kappa"], r.get("resolution"),
            r.get("target_mhz")) for r in recs
           if r["degree"] != 2 or r["kappa"] is not None or r.get("resolution") is not None
           or r.get("target_mhz") != 127.74]
    if bad or len({r["cells"] for r in recs}) != 1:
        raise RuntimeError(f"{stage}: stored node JSONs inconsistent: {bad}")
    for m, r in zip(nodes, recs):
        assert abs(r["f_hz"] - m * 1e6) < 1.0, (m, r["f_hz"])
        _assert_tie(r["cells"], _tie_from_json(r["tie"]))
    z_cu = [_uj(r["z"]["cu"]) for r in recs]
    xn = np.asarray(nodes, dtype=float) - 64.5
    _say0(comm, f"[ans7-xtr] {stage}: r = {r_esr!r} ohm (== ESR_VALUES[0] {ESR_VALUES[0]}), C_in(r) "
                f"{c_in_r!r} F, C_in(0) {c_in_0!r} F read from the stored tune_p2_esr0.2.json "
                f"(item 387); -n {comm.size}")

    # (R0'): the r = 0 branch is bitwise, on the stored x1 z64 at 127.74 MHz.
    f0hz = target_hz
    om0 = 2.0 * np.pi * f0hz
    old_cap = 1.0 / (1j * om0 * c_in_0)
    zb_u, zb_k = z_cap(om0, c_in_0), z_cap(om0, c_in_0, esr_ohm=0.0)
    ac_u, ac_k = add_c1(z64[0, 0], om0, c_in_0), add_c1(z64[0, 0], om0, c_in_0, esr_ohm=0.0)
    rb_u, rb_k = route_b_z(z64, f0hz, c_in_0), route_b_z(z64, f0hz, c_in_0, esr_ohm=0.0)
    r0 = {"z_cap": _bits(zb_u) == _bits(zb_k) == _bits(old_cap),
          "gap_z": _bits(_gap_z(f0hz, c_in_0)) == _bits(old_cap)
          == _bits(_gap_z(f0hz, c_in_0, 0.0)),
          "add_c1": _bits(ac_u) == _bits(ac_k),
          "route_b_z": all(_bits(a) == _bits(b) for a, b in zip(rb_u, rb_k))}
    r0_ok = all(r0.values())
    _say0(comm, f"[ans7-xtr] {stage} (R0') the r = 0 branch is bitwise at 127.74 MHz (z_cap, the gap "
                f"helper, add_c1, route_b_z with esr_ohm=0.0 vs the unkeyworded call; z_cap/helper vs "
                f"the literal old 1/(jwC)): {r0} (ASSERTED all True: {'met' if r0_ok else 'NOT met'})")
    if not r0_ok:
        failures.append(("R0'", r0))
        return failures

    t0 = time.perf_counter()
    geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
    t_mesh = time.perf_counter() - t0
    if cells != st["cells"]:
        failures.append(("cells", cells, st["cells"]))
    _say0(comm, f"[ans7-xtr] {stage}: exact hole mesh {cells} cells (stored {st['cells']}), "
                f"mesh {t_mesh:.1f} s")
    res, s_as, s_bs, z_in_a, z_in_b, t_drive = [], [], [], [], [], []
    for i, g in enumerate(tune):
        f = g * 1e6
        omega = 2.0 * np.pi * f
        t1 = time.perf_counter()
        s_a, _ = _route_a_solve(th14, geo, specs, f, "cu", c_f=c_in_r, degree=2, esr_ohm=r_esr)
        t_drive.append(time.perf_counter() - t1)
        zin_b, zbare_b, _g = route_b_z(bary_interp(xn, z_cu, g - 64.5), f, c_in_r, esr_ohm=r_esr)
        s_b = _s11_from_z(zbare_b)
        r = _band_residual(s_a, s_b)
        res.append(r)
        s_as.append(complex(s_a))
        s_bs.append(complex(s_b))
        z_in_a.append(add_c1(_z_from_s11(s_a), omega, c_in_r, esr_ohm=r_esr))
        z_in_b.append(complex(zin_b))
        _say0(comm, f"[ans7-xtr] T35 (AB-r') p2 Cu {g:.2f} MHz at (r, C_in(r)): S11_A "
                    f"{complex(s_a):.12e}, S11_B {complex(s_b):.12e}; |S11_A - S11_B|/|S11_B| = "
                    f"{r:.3e} (ASSERTED <= REDUCTION_BAND {REDUCTION_BAND:.0e}; PREDICTED "
                    f"{st['residual_ab'][i]:.3e}, item 387); drive {t_drive[-1]:.1f} s")
        if not r <= REDUCTION_BAND:
            failures.append(("AB-r'", g, r))
    rows_a = parallel_rows(np.array(tune) * 1e6, np.array(z_in_a))
    rows_b = parallel_rows(np.array(tune) * 1e6, np.array(z_in_b))
    b = (1.0 / np.asarray(z_in_a)).imag
    _say0(comm, f"[ans7-xtr] T35 Im Y_in,A at {list(tune)} MHz = " + ", ".join(f"{v:+.9e}" for v in b)
          + f" S; rows A: f0 {rows_a.get('f0_hz')!r} Hz, R_in {rows_a.get('r_in_ohm')!r} ohm, Q "
          f"{rows_a.get('q_loaded')!r}; up-crossings {rows_a['crossings_up_hz']}")
    f0a, f0b = rows_a.get("f0_hz"), rows_b.get("f0_hz")
    sa = (f0a / target_hz - 1.0) if f0a else float("nan")
    t2_ok = (len(rows_a["crossings_up_hz"]) == 1 and f0a is not None
             and tune[0] * 1e6 < f0a < tune[-1] * 1e6 and abs(sa) <= band)
    _say0(comm, f"[ans7-xtr] T35 (T2-r') f0,A(C_in(r)) {f0a!r} Hz; f0,A/127.74 MHz - 1 = {sa:+.6e} "
                f"(ASSERTED exactly one - -> + crossing in ({tune[0]}, {tune[-1]}) MHz and |.| <= "
                f"{band_name} {band:.3e}: {'met' if t2_ok else 'NOT met'}; = {abs(sa) / band:.3f} "
                f"bands; PREDICTED -2.301805e-6, item 387)")
    if not t2_ok:
        failures.append(("T2-r'", rows_a["crossings_up_hz"], sa))
    sb = (f0b / target_hz - 1.0) if f0b else float("nan")
    r_ab = (f0a / f0b - 1.0) if (f0a and f0b) else float("nan")
    _say0(comm, f"[ans7-xtr] T35 f0,B(C_in(r)) {f0b!r} Hz (Route B Z_in on the three tune points, as "
                f"item 378 forms it at 64 MHz); f0,B/127.74 MHz - 1 = {sb:+.6e} (PREDICTED -9.16e-13, "
                f"item 387); r_AB = f0,A/f0,B - 1 = {r_ab:+.6e} (PREDICTED -2.301804e-6, item 387); "
                f"records, not gates")
    # (R-T-r'): reproduction of item 387's registered run, against its stored JSON.
    st_a = st["rows_a"]
    d_f0 = abs(f0a / st_a["f0_hz"] - 1.0) if (f0a and "f0_hz" in st_a) else float("inf")
    d_res = max(abs(a - s) for a, s in zip(res, st["residual_ab"]))
    d_s = max(abs(complex(a) - complex(*s)) / abs(complex(*s)) for a, s in zip(s_as, st["s11_a"]))
    rt = {"f0_rel": d_f0, "residual_abs": d_res, "s11_a_rel": d_s}
    rt_assert = comm.size == 8
    rt_ok = max(d_f0, d_res, d_s) <= NODE_S_RTOL
    _say0(comm, f"[ans7-xtr] T35 (R-T-r') vs stored tune_inmodel_cin_esr0.2_p2.json: |f0 ratio - 1| "
                f"{d_f0:.3e}, max |res diff| {d_res:.3e}, max S11_A rel {d_s:.3e} ("
                f"{'ASSERTED' if rt_assert else 'printed only, width != 8'} <= {NODE_S_RTOL:.0e}: "
                f"{'met' if rt_ok else 'NOT met'}; PREDICTED <= 1e-12)")
    if rt_assert and not rt_ok:
        failures.append(("R-T-r'", rt))
    # (N-C-S'): the S11 comparator separates the ESR run from the lossless run -- asserted only
    # because item 387's stored s11_a and the lossless s11_a differ by > 1e-7 at every point.
    s_ll = ci0["s11_a"]
    sep_pts = [abs(complex(*a) - complex(*b)) / abs(complex(*b)) for a, b in zip(st["s11_a"], s_ll)]
    ncs_assert = all(v > 1.0e-7 for v in sep_pts)
    ncs = max(abs(a - complex(*b)) / abs(complex(*b)) for a, b in zip(s_as, s_ll))
    ncs_ok = ncs > 100.0 * NODE_S_RTOL
    _say0(comm, f"[ans7-xtr] T35 (N-C-S') max |S11_A,ex(r) - S11_A,stored(lossless)|/|.| = {ncs:.6e} = "
                f"{ncs / NODE_S_RTOL:.3e} x NODE_S_RTOL (stored 387-vs-lossless separation per point "
                f"{', '.join(f'{v:.3e}' for v in sep_pts)}; "
                f"{'ASSERTED > 100 x' if ncs_assert else 'PRINTED only, stored separation <= 1e-7'}: "
                f"{'met' if ncs_ok else 'NOT met'}; PREDICTED O(1e-3))")
    if ncs_assert and not ncs_ok:
        failures.append(("N-C-S'", ncs))
    # (N-C-r'): the f0 form -- asserted only if item 387's own f0,A(r) / f0,A(0) differ by > 1e-7.
    f0_ll = ci0["rows_a"]["f0_hz"]
    sep387 = abs(st_a["f0_hz"] / f0_ll - 1.0)
    n_c = abs(f0a / f0_ll - 1.0) if f0a else float("nan")
    n_c_assert = sep387 > 1.0e-7
    n_c_ok = n_c > 100.0 * NODE_S_RTOL
    _say0(comm, f"[ans7-xtr] T35 (N-C-r') |f0,A,ex(r)/f0,A,stored(lossless C_in) - 1| = {n_c:.6e} = "
                f"{n_c / NODE_S_RTOL:.1f} x NODE_S_RTOL (item 387's stored separation {sep387:.6e}; "
                f"{'ASSERTED > 100 x' if n_c_assert else 'PRINTED only, 387 separation <= 1e-7'}: "
                f"{'met' if n_c_ok else 'NOT met'}; PREDICTED 1e-8 to 5e-8, 64 MHz 4.7e-8)")
    if n_c_assert and not n_c_ok:
        failures.append(("N-C-r'", n_c))
    # (N-L'), printed unless item 387's own reading met the >= 10 x bar.
    j = [round(g, 6) for g in tune].index(round(target_hz / 1e6, 6))
    s_nl = _s11_from_z(route_b_z(z64, target_hz, c_in_0)[1])
    n_l = float(abs(s_nl - s_as[j]) / abs(s_as[j]))
    n_l_387_met = st["n_l"] >= 10.0 * REDUCTION_BAND
    n_l_ok = n_l >= 10.0 * REDUCTION_BAND
    _say0(comm, f"[ans7-xtr] T35 (N-L') lossless circuit at C_in(0) vs in-model at (r, C_in(r)), "
                f"127.74 MHz Cu: |dS11_bare|/|S11| = {n_l:.6e} = {n_l / REDUCTION_BAND:.2f} x "
                f"REDUCTION_BAND (>= 10 x: {'met' if n_l_ok else 'NOT met'}; item 387 measured "
                f"{st['n_l'] / REDUCTION_BAND:.2f} x, {'met' if n_l_387_met else 'not met'}; "
                f"{'ASSERTED' if n_l_387_met else 'PRINTED, never asserted'})")
    if n_l_387_met and not n_l_ok:
        failures.append(("N-L'", n_l))
    # Q(r)/Q_0 beside the closed form; C_in(r)/C_in(0) - 1 (printed, never asserted).
    d2_in = json.loads((CASE_DIR / "metrics.json").read_text()).get("degree2_exact", {})
    tl = d2_in.get("tuned_larmor", {}).get(key)
    if tl is None:
        raise RuntimeError(f"{stage}: metrics.json degree2_exact['tuned_larmor']['{key}'] absent")
    q0, r_in0 = tl["rows_a"]["q_loaded"], tl["rows_a"]["r_in_ohm"]
    qr = rows_a.get("q_loaded")
    q_c = 1.0 / (2.0 * np.pi * target_hz * c_in_r * r_esr)
    q_cf = 1.0 / (1.0 / q0 + 1.0 / q_c)
    q_ratio = qr / q0 if qr else float("nan")
    c_rel = c_in_r / c_in_0 - 1.0
    c_rel64 = d2_in.get("tuned_larmor_esr", {}).get("64.0", {}).get("c_in_rel")
    _say0(comm, f"[ans7-xtr] T35 Q(r) {qr!r} vs Q_0 = lossless tuned_larmor['127.74'] rows_a.q_loaded "
                f"{q0!r}: Q(r)/Q_0 = {q_ratio:.6f}; closed form Q_C = 1/(omega C_in r) = {q_c:.4f}, "
                f"1/Q = 1/Q_0 + 1/Q_C -> Q {q_cf:.4f}, ratio {q_cf / q0:.6f}; Q(r)/closed form - 1 = "
                f"{(qr / q_cf - 1.0) if qr else float('nan'):+.4%} (PREDICTED ratio 0.411537 vs "
                f"0.397431, item 387; printed, never asserted); R_in(r)/R_in_0 = "
                f"{rows_a.get('r_in_ohm', float('nan')) / r_in0:.6f}")
    _say0(comm, f"[ans7-xtr] T35 C_in(r)/C_in(0) - 1 = {c_rel:+.6e} beside tuned_larmor_esr['64.0'] "
                f"{c_rel64!r} (PREDICTED +1.344645e-5, item 387; printed, never asserted)")
    _say0(comm, f"[ans7-xtr] {stage}: total {time.perf_counter() - t0:.1f} s (mesh {t_mesh:.1f} s, drives "
                f"{', '.join(f'{t:.1f}' for t in t_drive)} s)")
    if failures or not rank0:
        return failures
    metrics_path = CASE_DIR / "metrics.json"
    text = metrics_path.read_text()
    m = json.loads(text)
    d2 = m["degree2_exact"]
    snap = _snapshot(m, skip=("degree2_exact",))
    snap_d2 = _snapshot(d2, skip=("tuned_larmor_esr",))
    snap_tl = _snapshot(d2.get("tuned_larmor", {}))
    snap_esr = _snapshot(d2.get("tuned_larmor_esr", {}))
    d2.setdefault("tuned_larmor_esr", {})[key] = {
        "esr_ohm": r_esr, "target_hz": target_hz, "c_in_r_f": c_in_r, "c_in_0_f": c_in_0,
        "c_in_rel": c_rel, "tune_mhz": list(tune),
        "s11_a": [[s.real, s.imag] for s in s_as], "s11_b": [[s.real, s.imag] for s in s_bs],
        "residual_ab": res, "reduction_band": float(REDUCTION_BAND), "rows_a": rows_a,
        "rows_b": rows_b, "f0a_rel_target": sa,
        "zero_resolution": {"value": band, "constant": band_name},
        "bands_fraction": abs(sa) / band, "r_ab": r_ab, "q_ratio": q_ratio,
        "q_closed_form": {"q_c": q_c, "q": q_cf, "ratio": q_cf / q0, "q_0": q0},
        "lossless_comparand": "tuned_larmor['127.74']", "rt": rt, "n_l": n_l,
        "n_c_r": n_c, "n_c_s": ncs, "cells": cells, "mpi_ranks": comm.size, "t_mesh_s": t_mesh,
        "t_drive_s": t_drive,
        "source": ["§9 item 387 (ANS-7-item387-W2-esr128 log)",
                   "stored t127.74 tune_p2.json / tune_p2_esr0.2.json / "
                   "tune_inmodel_cin_esr0.2_p2.json"]}
    metrics_path.write_text(json.dumps(m, indent=1) + "\n")
    after = json.loads(metrics_path.read_text())
    bad = [k for k, v in snap.items() if json.dumps(after.get(k), sort_keys=True) != v]
    bad += [("degree2_exact", k) for k, v in snap_d2.items()
            if json.dumps(after["degree2_exact"].get(k), sort_keys=True) != v]
    bad += [("tuned_larmor", k) for k, v in snap_tl.items()
            if json.dumps(after["degree2_exact"]["tuned_larmor"].get(k), sort_keys=True) != v]
    bad += [("tuned_larmor_esr", k) for k, v in snap_esr.items()
            if json.dumps(after["degree2_exact"]["tuned_larmor_esr"].get(k), sort_keys=True) != v]
    _say0(comm, f"[ans7-xtr] {stage} (K): {len(snap)} top-level + {len(snap_d2)} degree2_exact + "
                f"{len(snap_tl)} tuned_larmor + {len(snap_esr)} tuned_larmor_esr keys re-serialised "
                f"identically: {'yes' if not bad else 'NO ' + str(bad)} (ASSERTED)")
    if bad:
        metrics_path.write_text(text)
        failures.append(("K", bad))
        return failures
    _write_comparison_esr128(after)
    _say0(comm, f"[ans7-xtr] {stage}: wrote degree2_exact['tuned_larmor_esr']['{key}'] into "
                f"{metrics_path} and appended its COMPARISON.md section")
    return failures


def main() -> int:
    comm = MPI.COMM_WORLD
    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError("complex DolfinX build required (dolfinx-complex-mode)")
    if os.environ.get("FEM_EM_REQUIRE_COMPLEX") != "1":
        raise RuntimeError("FEM_EM_REQUIRE_COMPLEX=1 required")
    knobs = sorted(k for k in os.environ if k.startswith("FEM_EM_ANS6_"))
    if knobs:
        raise RuntimeError(f"every FEM_EM_ANS6_* knob must be unset; found {knobs}")
    stage = os.environ.get(STAGE_ENV, "").strip()
    if stage not in ("B1", "B2", "A", "X1", "X2", "X3", "X4", "XA", "XT64", "XT128", "XT64R",
                     "XT128R"):
        raise RuntimeError(f"{STAGE_ENV} must be one of B1, B2, A, X1, X2, X3, X4, XA, XT64, "
                           f"XT128, XT64R, XT128R; got {stage!r}")
    is_x = stage.startswith("X")
    if is_x:
        _x_guard()
    outdir = _p2_dir() if is_x else _results_dir()
    if comm.rank == 0:
        outdir.mkdir(parents=True, exist_ok=True)
        print(f"[ans7] stage {stage}; C = {C_F:.4e} F on every gap, gap 1 driven; degree "
              f"{2 if is_x else 1}; -n {comm.size}; results {outdir}", flush=True)
    comm.Barrier()
    ans6 = _load_example(ANS6_PATH, "ans6_for_ans7")
    t0 = time.perf_counter()
    if stage in ("XT64", "XT128"):
        failures = _xt_stage(comm, stage)
    elif stage == "XT64R":
        failures = _xt_esr_stage(comm)
    elif stage == "XT128R":
        failures = _xt_esr128_stage(comm)
    elif stage == "XA":
        failures = _x_stage_a(comm, outdir)
    elif is_x:
        failures = _x_nodes_stage(comm, stage, outdir)
    elif stage == "A":
        failures = _stage_a(comm, outdir)
    elif stage == "B1":
        failures = _stage_b1(comm, ans6, outdir)
    else:
        failures = _stage_b2(comm, ans6, outdir)
    comm.Barrier()
    failures = comm.bcast(failures, root=0)
    if comm.rank == 0:
        print(f"[ans7] stage {stage}: total {time.perf_counter() - t0:.1f} s at -n {comm.size}",
              flush=True)
        if failures:
            print(f"[ans7] ASSERTION FAILURES: {failures}", flush=True)
        else:
            print(f"[ans7] stage {stage}: all asserted checks green", flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
