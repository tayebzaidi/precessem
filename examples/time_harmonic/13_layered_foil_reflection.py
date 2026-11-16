"""Example (`EX-70`): a foil thinner than its skin depth — Γ of a slab-backed conductor layer on a driven TEM line.

A TEM line (box 1.2 x 0.04 x 0.04 m, E_z, 127.74 MHz) is driven from x = 0 and the
reflected wave is read from the standing pattern on x in [0.15, 0.85].  The line is
air to x = 1.0 and an eps_r = 4 dielectric slab beyond, terminated at x = 1.2 by a
**conductor layer** of conductivity 10 S/m carried as Jin §5.8.3's third-kind surface
term with Z_s from ``core/layered_impedance.py``.  Five terminations:

* ``leo``    bulk Leontovich, Z_s = (1 + j) R_s;
* ``t5``     a foil 5 skin depths thick — reflects as the bulk surface;
* ``thin``   a foil half a skin depth thick — leaks, and Γ moves by about 0.09,
             which the bulk term cannot show;
* ``pec``    a perfect conductor;
* ``copper`` 35 um of 5.8e7 S/m, printed only.

Everything gated is **imported** from ``tests/validation/test_th18_layered_reflection.py``
(``fem_gamma``, ``closed_form_gamma``, ``_deltas``, ``_fixture``, ``_delta``, the bands
``B_MAG`` / ``B_ARG`` / ``PREMISE_CAP`` / ``NEGATIVE_CONTROL_HAND`` and the geometry).
Γ comes only from ``fem_gamma`` (the 41-point least-squares fit); it is never re-fitted
from the exported CG1 interpolants.  Run it through the example runner::

    ./run_examples.sh -e th:13 -n 2 -t 180
"""

from __future__ import annotations

import csv
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

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_th18_layered_reflection import (  # noqa: E402
    AIR_TAG,
    B_ARG,
    B_MAG,
    COPPER_SIGMA,
    COPPER_T,
    LENGTH,
    NEGATIVE_CONTROL_HAND,
    PREMISE_CAP,
    SLAB_TAG,
    WIDTH,
    _delta,
    _fixture,
    _foil_zs,
    _leo_zs,
    closed_form_gamma,
    fem_gamma,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_13_layered_foil_reflection"

CASES = ("leo", "t5", "thin", "pec", "copper")

#: Gate records, 20260930T021508Z_TH-18-step2-W2.log:99-124 — printed beside this
#: run, a difference reported and never asserted.
RECORD = {  # case: (d|G|, d_arg rad)
    "leo": (7.185e-07, 1.354e-05),
    "t5": (7.184e-07, 1.354e-05),
    "thin": (1.176e-06, 1.417e-05),
    "pec": (3.331e-16, 1.423e-05),
}
RECORD_COPPER_DELTA = 5.242e-05
RECORD_CONTROL = 0.092963


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _zs(case):
    delta = _delta()
    return {
        "leo": _leo_zs(),
        "t5": _foil_zs(5.0 * delta),
        "thin": _foil_zs(0.5 * delta),
        "pec": 0.0,
        "copper": _foil_zs(COPPER_T, COPPER_SIGMA),
    }[case]


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
    _say(comm, "EX-70 -- a foil thinner than its skin depth: Gamma on a driven TEM line")
    _say(comm, "=" * 78)
    delta = _delta()
    _say(comm, f"\n[gate] bands imported: B_MAG {B_MAG:.0e}, B_ARG {B_ARG:.0e} rad, "
               f"PREMISE_CAP {PREMISE_CAP:.0e}, NEGATIVE_CONTROL_HAND {NEGATIVE_CONTROL_HAND}")
    _say(comm, f"[gate] skin depth {delta * 1e3:.6f} mm: 5 delta = {5 * delta * 1e3:.3f} mm, "
               f"delta/2 = {0.5 * delta * 1e3:.3f} mm")

    # ---- setup figure (EX-57): the gate's own mesh; the foil has no cells ----
    fx = _fixture()
    write_setup_figure(
        fx["msh"], fx["cell_tags"], FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={AIR_TAG: "air", SLAB_TAG: "dielectric slab eps_r = 4"},
        translucent_tags=(AIR_TAG, SLAB_TAG),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.5 * LENGTH, 0.5 * WIDTH, 0.5 * WIDTH),
        title="th:13 -- TEM line, air to x = 1.0 then eps_r = 4 slab; "
              "conductor layer on the x = 1.2 face; slice z = W/2",
        comm=comm,
    )

    # ---- solves (imported; fields kept for the export) ------------------------
    g_fem, fields = {}, {}
    t = time.perf_counter()
    for case in CASES:
        g_fem[case], fields[case] = fem_gamma(case, return_fields=True)
    _say(comm, f"\n[solves] five terminations ({time.perf_counter() - t:.1f} s)")

    # ---- deltas against the closed form, beside the record --------------------
    d = {}
    rows = []
    _say(comm, "")
    for case in CASES:
        g_cf = closed_form_gamma(case)
        d_mag = abs(abs(g_fem[case]) - abs(g_cf))
        d_arg = abs(float(np.angle(g_fem[case] / g_cf)))
        d[case] = (d_mag, d_arg)
        rec = RECORD.get(case)
        rec_txt = (f"record d|G| {rec[0]:.3e} d_arg {rec[1]:.3e}" if rec else "printed only")
        _say(comm, f"[{case:6s}] FEM |G|={abs(g_fem[case]):.8f} arg={np.angle(g_fem[case]):.8f} | "
                   f"cf |G|={abs(g_cf):.8f} arg={np.angle(g_cf):.8f} | "
                   f"d|G|={d_mag:.3e} d_arg={d_arg:.3e} rad ({rec_txt})")
        rows.append([case, _zs(case), abs(g_fem[case]), float(np.angle(g_fem[case])),
                     abs(g_cf), float(np.angle(g_cf)), d_mag, d_arg])

    # ---- assertions (imported bands, as the gate writes them) ------------------
    assert d["leo"][0] <= PREMISE_CAP and d["leo"][1] <= PREMISE_CAP, d["leo"]
    _say(comm, f"\n(1) leo d|G|, d_arg <= PREMISE_CAP {PREMISE_CAP:.0e}: ASSERTED")
    for n, case in enumerate(("t5", "thin", "pec"), start=2):
        assert d[case][0] <= B_MAG and d[case][1] <= B_ARG, (case, d[case])
        _say(comm, f"({n}) {case} d|G| <= {B_MAG:.0e}, d_arg <= {B_ARG:.0e}: ASSERTED")

    # ---- negative control: the bulk term against the thin-foil closed form ----
    nc = abs(g_fem["leo"] - closed_form_gamma("thin"))
    b = max(B_MAG, B_ARG)
    _say(comm, f"\n[control] |G_leo,FEM - G_thin,cf| = {nc:.6f} "
               f"(gate record {RECORD_CONTROL}; hand {NEGATIVE_CONTROL_HAND}); floor 10*B = {10 * b:.1e}")
    assert nc > 10.0 * b, nc
    _say(comm, "[control] miss > 10*max(B): ASSERTED")
    assert abs(nc - NEGATIVE_CONTROL_HAND) <= b, nc
    _say(comm, "[control] within max(B) of NEGATIVE_CONTROL_HAND: ASSERTED")

    # ---- printed only ------------------------------------------------------------
    cu = abs(g_fem["copper"] - g_fem["pec"])
    cu_cf = abs(closed_form_gamma("copper") - closed_form_gamma("pec"))
    _say(comm, f"\n[printed only] copper 35 um: |G_cu,FEM - G_pec,FEM| = {cu:.3e} "
               f"(record {RECORD_COPPER_DELTA:.3e}); cf |G_cu - G_pec| = {cu_cf:.3e}; never gated")

    # ---- export: gamma table (rank 0) + one combined XDMF ------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    csv_path = OUTPUT_DIR / f"{BASENAME}_gamma_table.csv"
    if comm.rank == 0:
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["case", "Zs_ohm", "fem_abs_gamma", "fem_arg_gamma_rad",
                        "cf_abs_gamma", "cf_arg_gamma_rad", "d_abs_gamma", "d_arg_gamma_rad"])
            for r in rows:
                w.writerow([r[0], f"{complex(r[1]):.8g}"] + [f"{x:.10e}" for x in r[2:]])
    comm.Barrier()
    msh = fx["msh"]
    v_vec = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    funcs = {}
    for case in CASES:
        e = fields[case].e_complex
        e_re = fem.Function(v_vec, name=f"E_real_{case}")
        e_re.interpolate(fem.Expression(ufl.real(e), v_vec.element.interpolation_points))
        e_im = fem.Function(v_vec, name=f"E_imag_{case}")
        e_im.interpolate(fem.Expression(ufl.imag(e), v_vec.element.interpolation_points))
        funcs[f"E_real_{case}"] = e_re
        funcs[f"E_imag_{case}"] = e_im
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, None, funcs, comm=comm,
        facet_tags=fx["facet_tags"],
    )
    _say(comm, f"\n[output] {path}")
    _say(comm, f"[output] {csv_path}")
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
