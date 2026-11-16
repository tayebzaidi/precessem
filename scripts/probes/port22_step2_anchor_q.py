"""§9 item 50 / chain T7 anchor (c), measured: fixed-C Cu and PEC resonance,
loaded Q and half-power width on the F-small hole mesh, by the circuit layer.

A design measurement -- no gate, band or record. Writes nothing (stdout only).

Route (PROJECT_PLAN §9 item 50):

* The 4x4: `ans:6`'s script loaded as a module (``_load_example`` pattern);
  its imported ``_build_ladder`` (``tests/validation/test_th14_birdcage_copper.py``)
  called with ``degree=1, sigmas=(COPPER_SIGMA,), resolution=None`` and the
  additive ``frequencies_hz`` keyword. ``main()`` / ``_knob_main()`` are never
  called, so no tracked or private file is written.
* The circuit layer: ``Z_in(f)`` from `PORT-15`'s ``tuned_input`` with
  ``PORT15_STEP3_C_TUNED_F`` (the SPEC's 15.57 pF, imported) on P2-P4, P1
  driven. Sheet widths uncorrected; kappa ~ 1 % moves f0 by ~kappa/2, which
  does not matter for a band-width reading.
* SPEC rows (a)-(c) (``ans7 .../SPEC.md:89-91``): f0 = zero of Im Z_in by
  linear interpolation between the bracketing grid points; R_in(f0)
  interpolated the same way; Q = f0 * (dIm Z_in/df) / (2 R_in).

Self-checks (ASSERTED; exit status 1 on failure): (1) at 64 MHz, both
columns' S within ``LEG_D0_REPRODUCTION_BAND`` of the tracked ``ans:6``
``metrics.json`` (read only); (2) reciprocity <= ``RECIPROCITY_BAND`` at every
point; (3) ``identity_residual`` <= ``DISCRETE_IDENTITY_RTOL`` at every point.
Negative control (ASSERTED, rule (e)): at 64 MHz, |S(Cu) - S(PEC)| > 1e-6 on
every C4 class.

Usage (complex build, every FEM_EM_ANS6_* knob unset)::

    mpiexec -n 4 python3 scripts/probes/port22_step2_anchor_q.py --window A
    mpiexec -n 4 python3 scripts/probes/port22_step2_anchor_q.py --window B --centre-mhz 63.2
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from dolfinx import default_scalar_type  # noqa: E402

SPEC_GRID_MHZ = np.arange(44.0, 84.0 + 1e-9, 2.0)  # ANS-7 SPEC: 21 points
WINDOW_B_HALF_MHZ = 1.0
WINDOW_B_STEP_MHZ = 0.1
NEGATIVE_CONTROL_FLOOR = 1.0e-6
ANS6_PATH = (REPO / "examples" / "ansys_benchmarks"
             / "ans6_copper_birdcage_four_port_10_64_128MHz"
             / "06_copper_birdcage_four_port_10_64_128MHz.py")


def _load_example(path: Path, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _crossings(f_mhz, x):
    """Every sign change of ``x`` between adjacent grid points.

    Returns (i, kind) with kind ``zero`` (- -> +, Foster reactance rising
    through zero: series resonance) or ``pole`` (+ -> -, the reactance jumping
    through infinity: a zero of 1/Im Z, parallel resonance / another mode).
    """
    out = []
    for i in range(len(f_mhz) - 1):
        if np.sign(x[i]) != np.sign(x[i + 1]) and x[i] != 0.0:
            out.append((i, "zero" if x[i] < 0.0 < x[i + 1] else "pole"))
    return out


def _row_abc(f_mhz, z, i):
    """SPEC rows (a)-(c) on the bracket (i, i+1)."""
    f1, f2 = f_mhz[i] * 1e6, f_mhz[i + 1] * 1e6
    x1, x2 = z[i].imag, z[i + 1].imag
    t = -x1 / (x2 - x1)
    f0 = f1 + t * (f2 - f1)
    r0 = z[i].real + t * (z[i + 1].real - z[i].real)
    slope = (x2 - x1) / (f2 - f1)
    q = f0 * slope / (2.0 * r0)
    return f0, r0, q


def _half_power_points(f_mhz, z, r0):
    """Interpolated f where Im Z_in = -R_in(f0) and +R_in(f0) (|Z| = sqrt2 R)."""
    x = np.array([v.imag for v in z])
    pts = {}
    for target, name in ((-r0, "lo"), (r0, "hi")):
        y = x - target
        for i in range(len(x) - 1):
            if y[i] < 0.0 <= y[i + 1]:
                pts[name] = f_mhz[i] + (-y[i] / (y[i + 1] - y[i])) * (f_mhz[i + 1] - f_mhz[i])
                break
    return pts


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", choices=("A", "B"), required=True)
    ap.add_argument("--centre-mhz", type=float, default=None)
    args = ap.parse_args()
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError("complex DolfinX build required (dolfinx-complex-mode)")
    if os.environ.get("FEM_EM_REQUIRE_COMPLEX") != "1":
        raise RuntimeError("FEM_EM_REQUIRE_COMPLEX=1 required (§9 item 50 trap)")
    knobs = sorted(k for k in os.environ if k.startswith("FEM_EM_ANS6_"))
    if knobs:
        raise RuntimeError(f"every FEM_EM_ANS6_* knob must be unset; found {knobs}")

    if args.window == "A":
        grid_mhz = SPEC_GRID_MHZ.copy()
    else:
        if args.centre_mhz is None:
            raise RuntimeError("--window B needs --centre-mhz (A's bracket)")
        k = np.arange(-10, 11)
        grid_mhz = np.round(args.centre_mhz + WINDOW_B_STEP_MHZ * k, 6)

    ans6 = _load_example(ANS6_PATH, "ans6_for_item50_anchor_q")
    from tests.validation.test_port22_driven_sweep_resonance import PORT15_STEP3_C_TUNED_F
    from tests.validation.test_port_circuit_layer_field import tuned_input

    c_f = float(PORT15_STEP3_C_TUNED_F)
    if rank0:
        print(f"[item50] window {args.window}: {len(grid_mhz)} points "
              f"{grid_mhz[0]:g}..{grid_mhz[-1]:g} MHz; C = {c_f:.4e} F on P2-P4, P1 driven; "
              f"-n {comm.size}; degree 1; hole mesh, pec_facet_tags=None route (TH-15)",
              flush=True)
        print("[item50] sheet widths uncorrected: kappa ~ 1 % moves f0 by ~kappa/2; "
              "immaterial to a band-width reading", flush=True)

    t0 = time.perf_counter()
    ladder = ans6._build_ladder(degree=1, sigmas=(ans6.COPPER_SIGMA,), resolution=None,
                                frequencies_hz=[float(f) * 1e6 for f in grid_mhz])
    comm.Barrier()
    elapsed = time.perf_counter() - t0

    failures = []
    cols = {"cu": [], "pec": []}
    keys = [f"{float(f):g}" for f in grid_mhz]
    for key, f in zip(keys, grid_mhz):
        rec = ladder["freqs"][key]
        cu = rec["sigma"][ans6.COPPER_SIGMA]
        pec = rec["pec"]
        for name, sw in (("cu", cu), ("pec", pec)):
            if not sw["reciprocity"] <= ans6.RECIPROCITY_BAND:
                failures.append(("reciprocity", key, name, sw["reciprocity"]))
            z_in, _ = tuned_input(np.asarray(sw["s"]), float(f) * 1e6, c_f)
            cols[name].append(complex(z_in))
        if not rec["identity_residual"] <= ans6.DISCRETE_IDENTITY_RTOL:
            failures.append(("identity", key, rec["identity_residual"]))
        if rank0:
            print(f"[item50] f = {f:8.3f} MHz: Z_in Cu {cols['cu'][-1].real:+.6e} "
                  f"{cols['cu'][-1].imag:+.6e}j  PEC {cols['pec'][-1].real:+.6e} "
                  f"{cols['pec'][-1].imag:+.6e}j ohm; recip Cu {cu['reciprocity']:.2e} "
                  f"PEC {pec['reciprocity']:.2e}; identity {rec['identity_residual']:.2e}",
                  flush=True)

    # -- self-check (1) and the negative control, at 64 MHz ------------------
    if "64" in ladder["freqs"]:
        import json
        rec = ladder["freqs"]["64"]
        cu, pec = rec["sigma"][ans6.COPPER_SIGMA], rec["pec"]
        tracked = json.loads((ans6.CASE_DIR / "metrics.json").read_text())
        view = {"frequency": "64 MHz",
                "cu": {"s_matrix": ans6._matrix_payload(np.asarray(cu["s"]))},
                "pec": {"s_matrix": ans6._matrix_payload(np.asarray(pec["s"]))}}
        worst, where, n = ans6._REPRO.knob_worst_deviation(view, tracked)
        band = ans6._REPRO.LEG_D0_REPRODUCTION_BAND
        if rank0:
            print(f"[item50] self-check 1: {n} S leaves at 64 MHz vs tracked metrics.json "
                  f"(-n {tracked['mpi_ranks']}; {ladder['cells']} vs {tracked['n_cells']} cells): "
                  f"worst rel {worst:.3e} at {where} (ASSERTED <= {band:.0e})", flush=True)
        if not (worst <= band and int(ladder["cells"]) == int(tracked["n_cells"])):
            failures.append(("reproduction", worst, where))
        ds = cu["class_ds"]
        if rank0:
            print("[item50] negative control |S(Cu)-S(PEC)| at 64 MHz: "
                  + ", ".join(f"{c} {ds[c]:.6e}" for c in ans6.CLASSES)
                  + f" (ASSERTED > {NEGATIVE_CONTROL_FLOOR:.0e})", flush=True)
        for c in ans6.CLASSES:
            if not ds[c] > NEGATIVE_CONTROL_FLOOR:
                failures.append(("negative-control", c, ds[c]))
    elif rank0:
        print("[item50] 64 MHz not on this grid: self-check 1 and the negative control "
              "are window A's", flush=True)

    # -- SPEC rows (a)-(c), per column ---------------------------------------
    if rank0:
        max_rel = (max(r.get("reciprocity", 0.0) for r in
                       [ladder["freqs"][k]["sigma"][ans6.COPPER_SIGMA] for k in keys]
                       + [ladder["freqs"][k]["pec"] for k in keys]))
        max_id = max(ladder["freqs"][k]["identity_residual"] for k in keys)
        print(f"[item50] self-check 2: max reciprocity {max_rel:.3e} (ASSERTED <= "
              f"{ans6.RECIPROCITY_BAND:.0e}); self-check 3: max identity residual "
              f"{max_id:.3e} (ASSERTED <= {ans6.DISCRETE_IDENTITY_RTOL:g})", flush=True)
        for name in ("cu", "pec"):
            z = np.array(cols[name])
            x = z.imag
            cr = _crossings(grid_mhz, x)
            for i, kind in cr:
                print(f"[item50] {name}: Im Z_in sign change {kind} in "
                      f"[{grid_mhz[i]:g}, {grid_mhz[i + 1]:g}] MHz "
                      f"({x[i]:+.4e} -> {x[i + 1]:+.4e} ohm); 1/Im Z_in "
                      f"{1 / x[i]:+.4e} -> {1 / x[i + 1]:+.4e}", flush=True)
            zeros = [i for i, kind in cr if kind == "zero"]
            if not zeros:
                j = int(np.argmin(np.abs(x)))
                print(f"[item50] {name}: NO zero of Im Z_in on this grid; |Im Z_in| "
                      f"smallest at {grid_mhz[j]:g} MHz ({x[j]:+.4e} ohm)", flush=True)
                continue
            for i in zeros:
                f0, r0, q = _row_abc(grid_mhz, z, i)
                width = f0 / q
                lo, hi = f0 - width / 2.0, f0 + width / 2.0
                n_spec = int(np.sum((SPEC_GRID_MHZ * 1e6 >= lo) & (SPEC_GRID_MHZ * 1e6 <= hi)))
                n_here = int(np.sum((grid_mhz * 1e6 >= lo) & (grid_mhz * 1e6 <= hi)))
                hp = _half_power_points(grid_mhz, z, r0)
                print(f"[item50] {name}: window {args.window} f0 = {f0 / 1e6:.6f} MHz, "
                      f"R_in(f0) = {r0:.6e} ohm, Q = {q:.4f}, width f0/Q = "
                      f"{width / 1e6:.6f} MHz; SPEC 2 MHz grid points in f0 +- f0/(2Q): "
                      f"{n_spec}; this window's points in band: {n_here}; "
                      f"half-power (|Im Z|=R_in) points {hp}", flush=True)
        print(f"[item50] window {args.window}: {len(grid_mhz)} frequencies, ladder "
              f"{elapsed:.1f} s at -n {comm.size}", flush=True)

    failures = comm.bcast(failures, root=0)
    if failures:
        if rank0:
            print(f"[item50] SELF-CHECK FAILURES: {failures}", flush=True)
        return 1
    if rank0:
        print(f"[item50] window {args.window}: all self-checks green", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
