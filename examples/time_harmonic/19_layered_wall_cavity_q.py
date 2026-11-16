"""Example (`EX-144`): a foil thinner than its skin depth on a cavity wall.

The `TH-14` box (1.0 x 0.8 x 0.6 m, TE101, degree 2, sigma = 1e4 S/m) with the
Leontovich wall (1 + j) R_s replaced, through
``solve_impedance_wall_cavity_mode(surface_impedance_fn=)``, by the layered input
impedance of ``core.layered_impedance`` (Jin 3e 3.4.2 (3.87)-(3.89)): a
conductor foil of fixed thickness t on a PEC backing.  The Q reads the wall's R,
the frequency shift reads its X.  A foil at delta0/2 cuts R to 0.16 R_s while X
stays 0.97 R_s, so Q rises 6.24x and the shift is nearly the Leontovich one.

**It asserts, and it does not re-implement** (the `ANS-1` rule).  Every fixture
constant and band is imported by name from
``tests/validation/test_th18_layered_impedance_cavity.py`` (`TH-18` step 1) and
``Q_TOLERANCE`` from ``tests/validation/test_cavity_leontovich_q.py``.  Asserted
as the gate asserts: (P) Q_leo against its record, (i) five skin depths is
Leontovich, (ii) the thin-foil Q ratio and Q, (iii) the thin-foil frequency
shift ratio on Re omega, and (c-bulk) the negative control: bulk Leontovich
misses the thin-layer prediction.  Records are printed beside the readings and
never asserted for equality.

Run it::

    ./run_examples.sh -e th:19 -n 2 -t 300

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_19_layered_wall_cavity_q_combined.xdmf`` (CG1 ``E_real``,
``E_imag`` of the thin-foil mode).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.core.cavity import (  # noqa: E402
    pozar_te10l_q_conductor,
    solve_impedance_wall_cavity_mode,
)
from fem_em_solver.core.layered_impedance import (  # noqa: E402
    ConductorLayer,
    conductor_skin_depth_m,
    conductor_surface_resistance_ohm,
    layered_surface_impedance,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.constants import ETA_0  # noqa: E402

from tests.validation.test_cavity_leontovich_q import Q_TOLERANCE  # noqa: E402
from tests.validation.test_th18_layered_impedance_cavity import (  # noqa: E402
    DEGREE,
    EDGES,
    FINE,
    NEGATIVE_CONTROL_FLOOR,
    PREMISE_TOLERANCE,
    Q_LEO_RECORD,
    Q_THIN_PREDICTED,
    SHIFT_RATIO,
    SHIFT_TOLERANCE,
    SIGMA,
    T5_RATIO,
    T5_TOLERANCE,
    THIN_RATIO,
    THIN_RATIO_TOLERANCE,
    _omega0,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_19_layered_wall_cavity_q"

#: Gate records -- printed beside this run, a difference reported and never asserted.
#: W1 = 20260929T200249Z_TH-18-step1-W1.log.gz
RECORD = {
    "cells": 2268,                      # W1:132
    "dofs": 15998,                      # W1:132
    "q_pec": 9.900075e18,               # W1:132
    "q_leo": 8.017679e02,               # W1:133
    "q_t5": 8.016672e02,                # W1:134
    "q_thin": 5.005122e03,              # W1:135
    "q_air": 3.992288e02,               # W1:139
    "thin_r_over_rs": 0.1601866860,     # W1:140
    "thin_x_over_rs": 0.9679795966,     # W1:140
    "q_c": 8.016849e02,                 # W1:154
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _cg1(msh, name):
    v = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    return v, fem.Function(v, name=name)


def _export(msh, e_h, comm):
    """CG1 interpolants of Re/Im E (qualitative picture only)."""
    v_vec, e_real = _cg1(msh, "E_real")
    e_real.interpolate(fem.Expression(ufl.real(e_h), v_vec.element.interpolation_points))
    _, e_imag = _cg1(msh, "E_imag")
    e_imag.interpolate(fem.Expression(ufl.imag(e_h), v_vec.element.interpolation_points))
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, None,
        {"E_real": e_real, "E_imag": e_imag}, comm=comm,
    )
    return path


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "This example needs the complex DolfinX build: "
            "source /usr/local/bin/dolfinx-complex-mode (the runner does this "
            "automatically for the `th:` group)."
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-144 -- a foil thinner than its skin depth on a cavity wall (TE101)")
    _say(comm, "=" * 78)
    omega0 = _omega0()
    delta0 = conductor_skin_depth_m(omega0, SIGMA)
    t5 = 5.0 * delta0
    t_half = 0.5 * delta0
    r_s = conductor_surface_resistance_ohm(omega0, SIGMA)
    _say(comm, f"\n[geometry] {EDGES} m box, N1curl degree {DEGREE}, divisions {FINE}, "
               f"sigma {SIGMA:g} S/m (imported from the TH-18 step 1 gate module)")
    _say(comm, f"[foil] omega0 = {omega0:.9e} rad/s, delta0 = {delta0:.9e} m, "
               f"t5 = {t5:.9e} m, delta0/2 = {t_half:.9e} m, R_s = {r_s:.6e} ohm")

    def stack(thickness, backing):
        return lambda w: layered_surface_impedance(w, [ConductorLayer(SIGMA, thickness)], backing)

    cases = (
        ("pec", dict(sigma_s_per_m=None)),
        ("leo", dict(sigma_s_per_m=SIGMA)),
        ("t5", dict(surface_impedance_fn=stack(t5, 0.0))),
        ("thin", dict(surface_impedance_fn=stack(t_half, 0.0), return_mode=True)),
        ("air", dict(surface_impedance_fn=stack(t_half, ETA_0))),
    )
    out = {}
    for key, kw in cases:
        t = time.perf_counter()
        out[key] = solve_impedance_wall_cavity_mode(
            edges=EDGES, divisions=FINE, degree=DEGREE, comm=comm, **kw
        )
        r = out[key]
        _say(comm, f"[{key:4s}] cells={r.n_cells} dofs={r.n_dofs} nconv={r.n_converged} "
                   f"Z_s={r.surface_impedance_ohm} f={r.omega_rad_s.real / (2 * np.pi) / 1e6:.6f} MHz "
                   f"Q={r.q:.6e} ({time.perf_counter() - t:.1f} s)")
    pec, leo, t5r, thin, air = (out[k] for k in ("pec", "leo", "t5", "thin", "air"))

    _say(comm, f"\n[records] cells {pec.n_cells} (record {RECORD['cells']}), dofs {pec.n_dofs} "
               f"(record {RECORD['dofs']})")
    _say(comm, f"[records] Q: pec {pec.q:.6e} ({RECORD['q_pec']:.6e}); leo {leo.q:.6e} "
               f"({RECORD['q_leo']:.6e}); t5 {t5r.q:.6e} ({RECORD['q_t5']:.6e}); thin {thin.q:.6e} "
               f"({RECORD['q_thin']:.6e}); air {air.q:.6e} ({RECORD['q_air']:.6e})")
    for key in ("t5", "thin", "air"):
        z = out[key].surface_impedance_ohm
        extra = (f" (record {RECORD['thin_r_over_rs']:.10f} / {RECORD['thin_x_over_rs']:.10f})"
                 if key == "thin" else "")
        _say(comm, f"[wall] {key:4s}: R/R_s = {z.real / r_s:.10f}, X/R_s = {z.imag / r_s:.10f}{extra}")

    # ---- setup figure (EX-57): the solver owns the mesh, so render on it ----
    write_setup_figure(
        thin.mesh,
        None,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={0: "cavity (vacuum); every wall a delta0/2 foil on PEC"},
        slice_normal=(0.0, 1.0, 0.0),
        title="th:19 -- TH-14 box, delta0/2 copper-like foil on PEC on every face (TE101)",
        comm=comm,
    )

    # ---- gates, as the TH-18 step 1 gate asserts them -------------------------
    q_c = pozar_te10l_q_conductor(EDGES, 1, SIGMA)

    rel_p = abs(leo.q / Q_LEO_RECORD - 1.0)
    _say(comm, f"\n(P)  Q_leo = {leo.q:.7e} vs record {Q_LEO_RECORD:.6e}: rel {rel_p:.3e} "
               f"<= {PREMISE_TOLERANCE:.0e}: ASSERTED")

    ratio_t5 = t5r.q / leo.q
    d_t5 = abs(ratio_t5 - T5_RATIO)
    _say(comm, f"(i)  Q_t5/Q_leo = {ratio_t5:.10f} vs {T5_RATIO}: |diff| {d_t5:.3e} "
               f"<= {T5_TOLERANCE:.0e}: ASSERTED (five skin depths is Leontovich)")

    ratio_thin = thin.q / leo.q
    miss_ratio = ratio_thin / THIN_RATIO - 1.0
    miss_q = thin.q / Q_THIN_PREDICTED - 1.0
    _say(comm, f"(ii) Q_thin/Q_leo = {ratio_thin:.7f} vs {THIN_RATIO}: miss {100 * miss_ratio:+.4f}% "
               f"(band {100 * THIN_RATIO_TOLERANCE:.1f}%): ASSERTED")
    _say(comm, f"(ii) Q_thin = {thin.q:.6e} vs {Q_THIN_PREDICTED}: miss {100 * miss_q:+.4f}% "
               f"(band {100 * Q_TOLERANCE:.0f}%): ASSERTED; Q_c = {q_c:.6e} (record {RECORD['q_c']:.6e})")

    f = {k: out[k].omega_rad_s.real / (2 * np.pi) for k in ("pec", "leo", "thin", "air")}
    shift = (f["thin"] - f["pec"]) / (f["leo"] - f["pec"])
    miss_shift = shift / SHIFT_RATIO - 1.0
    _say(comm, f"(iii) f_pec {f['pec'] / 1e6:.6f}, f_leo {f['leo'] / 1e6:.6f}, f_thin "
               f"{f['thin'] / 1e6:.6f} MHz; (f_thin-f_pec)/(f_leo-f_pec) = {shift:.7f} vs "
               f"{SHIFT_RATIO}: miss {100 * miss_shift:+.4f}% (band {100 * SHIFT_TOLERANCE:.0f}%): ASSERTED")

    gap = 1.0 - leo.q / Q_THIN_PREDICTED
    _say(comm, f"(c-bulk) 1 - Q_leo/Q_thin_predicted = {gap:.4f} > {NEGATIVE_CONTROL_FLOOR}: ASSERTED "
               "(bulk Leontovich misses the thin layer; ceiling 1)")

    _say(comm, f"\n[printed] air (delta0/2 on eta0): Q = {air.q:.6e}, Q_c/Q_air = {q_c / air.q:.4f}, "
               f"f = {f['air'] / 1e6:.6f} MHz; Q_thin/Q_air = {thin.q / air.q:.3f} "
               "(the foil on free space is lossier than bulk metal; on PEC it is less lossy)")

    # ---- ParaView -------------------------------------------------------------
    path = _export(thin.mesh, thin.mode_function, comm)
    _say(comm, f"[output] {path}")

    assert rel_p <= PREMISE_TOLERANCE, leo.q
    assert d_t5 <= T5_TOLERANCE, ratio_t5
    assert abs(miss_ratio) <= THIN_RATIO_TOLERANCE, ratio_thin
    assert abs(miss_q) <= Q_TOLERANCE, thin.q
    assert abs(miss_shift) <= SHIFT_TOLERANCE, shift
    assert gap > NEGATIVE_CONTROL_FLOOR, gap

    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
