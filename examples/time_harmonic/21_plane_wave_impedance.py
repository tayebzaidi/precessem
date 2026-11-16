"""Example (`EX-146`): the wave impedance of a lossy plane wave, from Faraday's law.

`POST-2` (closed 2026-09-30) read the field-consistency diagnostics on `TH-6`'s
lossy plane wave ``E = ẑ exp(−j k x)``. Since ``∇×E = j k exp(−j k x) ŷ``, the
magnetic flux density ``B = ∇×E/(−jω)`` obtained from the *solved* E has
``|E|/|B| = ω/|k|`` at every point. In a conductor ``|k|`` is set by the complex
permittivity, so ω/|k| is **not** the lossless phase speed c/√εᵣ: in this
saline (σ = 0.7 S/m, εᵣ = 78, 127.74 MHz) a σ-blind reading misses by ≈ 21 %.

Everything is imported by name from the gate modules -- the band, the diagnostics
helper and the closed forms -- so the example and the gates cannot drift apart:

* (R)   the ratio of the means |E|/|B| against ω/|k|, inside ``POST2_BAND``;
* (S_E) and (S_B) the E and B span ratios against ``1 − e^{−αΔx_c}``;
* (Bal) the dimensionless balance at v = ω/|k|, inside ``POST2_BAND``;
* (c-σ) the negative control, asserted: against the σ-blind c/√εᵣ the same
  ratio misses by more than ``POST2_BAND``.

Records are printed beside the readings. Nothing asserts digit equality.

Run it::

    ./run_examples.sh -e th:21 -n 2 -t 600

Output lands in ``examples/time_harmonic/paraview_output``: open
``time_harmonic_21_plane_wave_impedance_combined.xdmf`` and colour by
``E_magnitude`` / ``B_magnitude`` (the P1 export of an N1curl or DG field is a
qualitative picture).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post import compute_field_consistency_diagnostics  # noqa: E402
from fem_em_solver.post.faraday import (  # noqa: E402
    magnetic_flux_density_from_e,
    project_to_cg1,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.constants import EPSILON_0, MU_0  # noqa: E402

from tests.validation.test_field_consistency_metrics import (  # noqa: E402
    POST2_BAND,
    _post2_read_diagnostics,
)
from tests.validation.test_lossy_plane_wave import (  # noqa: E402
    EPSILON_R,
    FREQUENCY_HZ,
    SIGMA,
    _analytic_alpha_beta,
    _solve_plane_wave,
    _wavenumber,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "time_harmonic_21_plane_wave_impedance"

#: Records of the printed readings (filled in ``main``).
RECORD: dict = {}


def _paraview_fields(msh, e_complex, b_dg0):
    """CG1 E_real / E_imag / E_magnitude / B_real / B_imag / B_magnitude."""
    v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    s_cg = fem.functionspace(msh, ("Lagrange", 1))

    e_cg = fem.Function(v_cg, name="E_phasor")
    e_cg.interpolate(e_complex)
    e_cg.x.scatter_forward()
    # B is DG0: an L2 projection (mass solve), never an interpolation.
    b_cg = project_to_cg1(b_dg0)

    out = {}
    for tag, f in (("E", e_cg), ("B", b_cg)):
        re = fem.Function(v_cg, name=f"{tag}_real")
        re.x.array[:] = np.real(f.x.array)
        im = fem.Function(v_cg, name=f"{tag}_imag")
        im.x.array[:] = np.imag(f.x.array)
        mag = fem.Function(s_cg, name=f"{tag}_magnitude")
        comps = np.abs(f.x.array.reshape(-1, 3))
        mag.x.array[:] = np.sqrt(np.sum(comps * comps, axis=1))
        for g in (re, im, mag):
            g.x.scatter_forward()
        out[f"{tag}_real"], out[f"{tag}_imag"], out[f"{tag}_magnitude"] = re, im, mag
    return out


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `th:` "
            "group sources it automatically"
        )

    omega = 2.0 * np.pi * FREQUENCY_HZ
    k = _wavenumber(SIGMA)
    alpha_cf, beta_cf = _analytic_alpha_beta(SIGMA)
    alpha = -k.imag
    ratio_cf = omega / abs(k)
    c0 = 1.0 / np.sqrt(MU_0 * EPSILON_0)
    ratio_blind = c0 / np.sqrt(EPSILON_R)

    if rank0:
        print("=" * 72)
        print("EX-146 -- lossy plane wave: |E|/|B| from Faraday's law vs omega/|k|")
        print("=" * 72)
        print(
            f"\n[material] sigma = {SIGMA} S/m, eps_r = {EPSILON_R} at "
            f"f = {FREQUENCY_HZ / 1e6:.2f} MHz"
            f"\n[closed form] omega/|k| = {ratio_cf:.6e} m/s;"
            f" sigma-blind c/sqrt(eps_r) = {ratio_blind:.6e} m/s"
            f"\n[closed form] alpha = {alpha_cf:.6f} Np/m, beta = {beta_cf:.6f} rad/m",
            flush=True,
        )

    # ---- the gate's two rungs (lru-cached helper, positional int) -------------
    rungs = {}
    for n in (24, 12):
        t0 = time.perf_counter()
        diag, e_stats, b_stats, dx_c, a_fit, b_fit, ncells = _post2_read_diagnostics(n)
        seconds = time.perf_counter() - t0
        span_cf = 1.0 - np.exp(-alpha * dx_c)
        rec = {
            "cells": ncells,
            "samples_E": e_stats["count"],
            "samples_B": b_stats["count"],
            "alpha_fit": a_fit,
            "beta_fit": b_fit,
            "mean_E": e_stats["mean"],
            "mean_B": b_stats["mean"],
            "mean_ratio": diag["e_to_b_mean_ratio"],
            "mean_miss": diag["e_to_b_mean_ratio"] / ratio_cf - 1.0,
            "max_miss": diag["e_to_b_max_ratio"] / ratio_cf - 1.0,
            "dx_c": dx_c,
            "span_cf": span_cf,
            "e_span_miss": diag["e_span_ratio"] / span_cf - 1.0,
            "b_span_miss": diag["b_span_ratio"] / span_cf - 1.0,
            "control_miss": diag["e_to_b_mean_ratio"] / ratio_blind - 1.0,
        }
        rungs[n] = rec
        if rank0:
            print(
                f"\n[rung {n}^3] cells = {ncells}, samples E/B = "
                f"{rec['samples_E']}/{rec['samples_B']} ({seconds:.1f} s)"
                f"\n  premise fit: alpha = {a_fit:.5f} (cf {alpha_cf:.5f}),"
                f" beta = {b_fit:.5f} (cf {beta_cf:.5f})"
                f"\n  mean |E| = {rec['mean_E']:.6e} V/m, mean |B| = {rec['mean_B']:.6e} T"
                f"\n  mean ratio = {rec['mean_ratio']:.6e} vs omega/|k| = {ratio_cf:.6e}"
                f"  miss {rec['mean_miss'] * 100:+.4f} %"
                f"\n  max ratio miss {rec['max_miss'] * 100:+.4f} % (printed, not gated)"
                f"\n  dx_c = {dx_c:.6f}, span_cf = {span_cf:.6f}"
                f"\n  E span miss {rec['e_span_miss'] * 100:+.4f} %,"
                f" B span miss {rec['b_span_miss'] * 100:+.4f} %"
                f"\n  control vs c/sqrt(eps_r): miss {rec['control_miss'] * 100:+.4f} %",
                flush=True,
            )

    # ---- (Bal): the dimensionless balance at v = omega/|k| ---------------------
    _d, e_stats, b_stats, *_rest = _post2_read_diagnostics(24)
    at_cf = compute_field_consistency_diagnostics(
        e_stats, b_stats, reference_speed_m_per_s=ratio_cf
    )
    at_blind = compute_field_consistency_diagnostics(
        e_stats, b_stats, reference_speed_m_per_s=ratio_blind
    )
    r = at_cf["e_to_b_mean_ratio"]
    identity = abs(r - ratio_cf) / max(r, ratio_cf)
    balance = at_cf["mean_balance_rel_diff"]
    balance_blind = at_blind["mean_balance_rel_diff"]
    if rank0:
        print(
            f"\n[balance] at v = omega/|k|: {balance:.6e}"
            f" vs |r-v|/max(r,v) = {identity:.6e}"
            f" (gate's 1e-9 identity: ratio - 1 = {balance / identity - 1.0:+.3e}, printed only)"
            f"\n[balance] at v = c/sqrt(eps_r): {balance_blind:.6e}"
            f" ({balance_blind * 100:.4f} %, predicted 21.20 %)",
            flush=True,
        )

    # ---- the picture: one more 24^3 solve, E and B to one XDMF ------------------
    t0 = time.perf_counter()
    _rel, _a, _b, _nc, msh, fields = _solve_plane_wave(24, SIGMA, return_fields=True)
    resolve_seconds = time.perf_counter() - t0
    write_setup_figure(
        msh,
        None,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={0: "box (homogeneous lossy medium, no cell-tag regions)"},
        slice_normal=(0.0, 1.0, 0.0),
        title="th:21 -- plane-wave impedance: 0.1 m cube, homogeneous lossy medium, 24^3",
        comm=comm,
    )
    b_dg0 = magnetic_flux_density_from_e(fields.e_complex, omega)
    pv = _paraview_fields(msh, fields.e_complex, b_dg0)
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, None, pv, comm=comm
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    if rank0:
        print(
            f"\n[paraview] wrote {xdmf_path} (re-solve {resolve_seconds:.1f} s);"
            " the P1 export of an N1curl / DG field is a qualitative picture",
            flush=True,
        )

    m24, m12 = rungs[24], rungs[12]
    RECORD.update(
        omega_over_k=ratio_cf,
        blind=ratio_blind,
        rung24=m24,
        rung12=m12,
        balance=balance,
        balance_blind=balance_blind,
    )

    if rank0:
        print(
            f"\n[record] 24^3 mean {m24['mean_miss'] * 100:+.4f} %,"
            f" E span {m24['e_span_miss'] * 100:+.4f} %,"
            f" B span {m24['b_span_miss'] * 100:+.4f} %,"
            f" control {m24['control_miss'] * 100:+.4f} %"
            f"\n[record] 12^3 mean {m12['mean_miss'] * 100:+.4f} %,"
            f" E span {m12['e_span_miss'] * 100:+.4f} %,"
            f" B span {m12['b_span_miss'] * 100:+.4f} %,"
            f" control {m12['control_miss'] * 100:+.4f} %"
            f"\n\nTotal elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )

    # ---- gates (POST-2's 24^3 asserts) -----------------------------------------
    assert abs(m24["mean_miss"]) < POST2_BAND, (
        f"(R) mean-ratio miss {m24['mean_miss']:.4%} at 24^3"
    )
    assert abs(m24["e_span_miss"]) < POST2_BAND, (
        f"(S_E) E span miss {m24['e_span_miss']:.4%} at 24^3"
    )
    assert abs(m24["b_span_miss"]) < POST2_BAND, (
        f"(S_B) B span miss {m24['b_span_miss']:.4%} at 24^3"
    )
    assert balance < POST2_BAND, f"(Bal) balance {balance:.4%} at v = omega/|k|"
    assert abs(m24["control_miss"]) > POST2_BAND, (
        f"(c-sigma) sigma-blind comparand miss only {m24['control_miss']:.4%}: "
        "sigma did not reach the solve"
    )
    if rank0:
        print("All assertions hold.", flush=True)


if __name__ == "__main__":
    main()
