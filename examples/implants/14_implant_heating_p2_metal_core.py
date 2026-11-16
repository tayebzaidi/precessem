"""Example (`EX-154`): the implant-heated gel ball at P2 -- element order on `imp:11`'s rung.

`imp:11` solves `IMP-7` step 1's metal-core case at P1 on the 1.0 mm rung
(`LADDER[1]`). `IMP-7` step 2 (I40) added a degree-2 route; no thermal example
runs at degree 2. This one solves the same fixture at P1 and then at P2 in the
complex build and narrates the cut in error.

**It asserts, and it does not re-implement** (the `ANS-1` rule). The fixture,
bands, degree route and helpers are imported by name from the gate module
(``test_imp7_implant_heating_fem``) as ``m7``. Two disclosed local pieces: the
gate's ``_metal_rung`` calls ``_errors``, which squares the error as ``d**2`` --
not a norm in the complex build -- so ``_rung`` writes ``_metal_rung``'s body
out with ``rt = DEG_ROUTE[degree]`` and a local ``_errors_c`` using
``ufl.inner(d, d)``, carrying ``rt["qmeta"]`` as the quadrature metadata.

Asserted: (S3) balance <= S3_BAND at P1 and P2; (S2-P1) and (S2-P2) interior
rel L2 <= S2_BAND; (c-beta0) the beta = 0 source at P2 on CONTROL_RES misses
S2_BAND (backed by the gate's own measurement).

Printed only: P2 / P1 error ratios, the l = 2 readout R2_h and the pole /
equator rise at r = 1.1a beside SERIES_READ, cells, build and solve seconds.
No rate is shown (one rung).

Run it::

    ./scripts/run_examples.sh -e imp:14 -n 2 -t 420

Output: ``examples/implants/paraview_output/
implants_14_implant_heating_p2_metal_core_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.validation.test_imp7_implant_heating_fem as m7  # noqa: E402
from tests.validation.test_imp7_implant_heating_closed_form import T_INF  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_14_implant_heating_p2_metal_core"
RES = m7.LADDER[1]
assert m7.CONTROL_RES == RES


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _errors_c(msh, T_h, theta, qmeta):
    """(whole, interior |x|^2 < 0.81 b^2) rel L2 on theta; `inner`, complex-safe."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=qmeta)
    d = (T_h - T_INF) - theta
    chi_in = m7._indicator(msh, lambda m: np.sum(m**2, axis=1) < 0.81 * m7.B**2)
    whole = np.sqrt(max(m7._allsum(comm, ufl.inner(d, d) * dx), 0.0)
                    / m7._allsum(comm, ufl.inner(theta, theta) * dx))
    inner_ = np.sqrt(max(m7._allsum(comm, chi_in * ufl.inner(d, d) * dx), 0.0)
                     / m7._allsum(comm, chi_in * ufl.inner(theta, theta) * dx))
    return float(whole), float(inner_)


def _rung(degree, q_fn=m7._q_shell):
    """`m7._metal_rung`'s body (m7 :335-345), complex-safe."""
    rt = m7.DEG_ROUTE[degree]
    msh, ct, ft, t_build = m7._ball(RES)
    k_f = m7._dg0(msh, ct, m7.K_CORE, m7.K_GEL)
    q_f = m7._source_dg2(msh, ct.find(2), q_fn, dg=rt["q_dg"])
    sol, t_solve = m7._solve(msh, ft, k_f, q_f, degree=degree)
    theta = m7._p2(msh, m7._theta_metal, p=rt["comparand"])
    e_all, e_in = _errors_c(msh, sol.temperature, theta, rt["qmeta"])
    src = float(np.real(sol.source_power))
    bal = abs(float(np.real(sol.convected_power)) - src) / src
    return msh, ct, sol, q_f, theta, t_build, t_solve, e_all, e_in, bal


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-154 -- the implant-heated gel ball at P2 (IMP-7 step 2, one rung)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: S2 interior {m7.S2_BAND:.0e}, S3 balance {m7.S3_BAND:.0e}; "
               f"rung {RES} m; P2 route {m7.DEG_ROUTE[2]}")

    msh, ct, ft, t_build = m7._ball(RES)
    nglob = msh.topology.index_map(3).size_global
    _say(comm, f"[mesh] res {RES}: {nglob} cells, build {t_build:.1f} s")

    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: f"metal core a = 5 mm (k {m7.K_CORE} W/m/K, no source)",
            2: f"gel shell to b = 20 mm (k {m7.K_GEL} W/m/K, dipole SAR)",
        },
        translucent_tags=(2,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:14 -- metal-core sphere in a gel ball, solved at P1 and P2 (IMP-7 fixture)",
        window_size=(1400, 630),
        comm=comm,
    )

    th_pole = float(m7._theta_metal(np.array([[0.0], [0.0], [m7.R0_READ]]))[0])
    th_eq = float(m7._theta_metal(np.array([[m7.R0_READ], [0.0], [0.0]]))[0])
    sr = m7.SERIES_READ["metal"]
    res = {}
    for degree in (1, 2):
        _, _, sol, q_f, theta, tb, ts, e_all, e_in, bal = _rung(degree)
        res[degree] = (sol, q_f, theta, e_all, e_in, bal)
        _say(comm, f"\n[P{degree}] cells {nglob}, solve {ts:.1f} s: whole rel L2 {e_all:.4e}, "
                   f"interior rel L2 {e_in:.4e}  ASSERTED <= {m7.S2_BAND:.0e}; "
                   f"balance {bal:.3e}  ASSERTED <= {m7.S3_BAND:.0e}")
        m7._l2_readout(sol.temperature, f"metal P{degree}", sr)

    p1 = m7.P1_READ["metal"]
    _say(comm, f"[readout] series (R2, pole, equator) = {sr[0]:.4e}, {sr[1]:.6f}, {sr[2]:.6f} "
               f"(closed form pole / equator {th_pole:.6f} / {th_eq:.6f}); gate P1 finest-rung record "
               f"{p1[0]:.4e}, {p1[1]:.6f}, {p1[2]:.6f}. R2_h is a difference of two near-equal "
               f"numbers: a cancellation readout, no band.")
    _say(comm, f"[ratio] P2 / P1 whole L2 {res[2][3] / res[1][3]:.4f}; "
               f"interior {res[2][4] / res[1][4]:.4f} (printed; one rung, no rate shown)")

    # ---- (c-beta0) negative control at P2 -----------------------------------------------------
    _, _, _, _, _, _, ts0, e_all0, e_in0, _ = _rung(2, q_fn=lambda x: m7._q_shell(x, beta=0.0))
    _say(comm, f"[control] (c-beta0) P2, beta = 0 source on {RES} m: interior miss {e_in0:.4e}; "
               f"ASSERTED > {m7.S2_BAND:.0e} ({e_in0 / m7.S2_BAND:.1f}x the band; P1 read "
               f"{m7.P1_BETA0['metal']:.3e}); whole {e_all0:.4e}; solve {ts0:.1f} s")

    # ---- assertions ---------------------------------------------------------------------------
    for d in (1, 2):
        assert res[d][5] <= m7.S3_BAND, f"(S3) P{d} balance {res[d][5]:.3e}"
        assert res[d][4] <= m7.S2_BAND, f"(S2) P{d} interior {res[d][4]:.4e}"
    assert e_in0 > m7.S2_BAND, f"(c-beta0) did not fire: {e_in0:.4e}"

    # ---- outputs ------------------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    V1 = fem.functionspace(msh, ("Lagrange", 1))
    fields = {}
    for d in (1, 2):
        f = fem.Function(V1, name=f"theta_p{d}")
        f.interpolate(res[d][0].temperature)
        f.x.array[:] = np.real(f.x.array) - T_INF
        f.x.scatter_forward()
        fields[f"theta_p{d}"] = f
    th_ex = fem.Function(V1, name="theta_exact")
    th_ex.interpolate(res[2][2])
    th_ex.x.scatter_forward()
    fields["theta_exact"] = th_ex
    xdmf = OUTPUT_DIR / f"{BASENAME}_combined"
    write_xdmf_with_tags(xdmf, msh, ct, fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf}.xdmf -- P1-sampled fields on the {nglob}-cell mesh")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
