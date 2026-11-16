"""Example (`EX-140`): the gel heated by an implant's own SAR.

`imp:1` shows the PEC sphere's SAR field in gel; `mat:12` shows conduction
through a two-layer gel ball. No example shows a temperature rise *driven by an
implant's own SAR*, or the metal core's effect on it. This one runs `IMP-7`
step 1's metal-core case (the middle rung of its ladder, 1.0 mm, `LADDER[1]`) in
the complex build and narrates it: a 5 mm high-k metal core (k = 20 W/m/K) in a
gel shell out to b = 4a with Robin cooling, heated by the beta = 1 quasi-static
dipole SAR of the PEC sphere (shell only), solved at P1 against the closed form.

**It asserts, and it does not re-implement** (the `ANS-1` rule). The fixture,
bands and helpers are imported by name from `IMP-7` step 1
(``test_imp7_implant_heating_fem.py``), its closed-form module
(``test_imp7_implant_heating_closed_form.py``) and `MAT-18`'s FEM module
(``test_mat18_two_layer_sphere_fem.py``). The one local re-definition is
``_errors_c``: the gate's ``_errors`` squares the error as ``d**2``, which in
the complex build is not the norm, so this example computes the identical three
quantities (whole / interior / shell relative L2) with ``ufl.inner(d, d)``.

Asserted: (S2) interior rel L2 <= S2_BAND; (S3) balance <= S3_BAND; (Src')
the unmasked source continuation equals the gated closed-form source to
B1_BAND; (Im) max |Im T_h| = 0; (c-iso) the isotropic source (l = 2 term
dropped) misses the interior band.

Printed only, never gated: the assembled source power beside the closed form
and the gate's record; pole / equator rise at r = 1.1a and the l = 2 readout
(a difference of two near-equal numbers); the core's near-uniform temperature;
the heating factor of the dipole field.

Run it::

    ./run_examples.sh -e imp:11 -n 2 -t 180

Output: ``examples/implants/paraview_output/
implants_11_pec_sphere_implant_heating_in_gel_ball_combined.xdmf``.
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
from fem_em_solver.utils.analytical import (  # noqa: E402
    pec_sphere_implant_radial,
    pec_sphere_implant_source,
    pec_sphere_implant_steady_temperature,
)

from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND  # noqa: E402
from tests.validation.test_imp7_implant_heating_closed_form import (  # noqa: E402
    H,
    K_CORE,
    K_GEL,
    Q_FAR,
    T_INF,
)
from tests.validation.test_imp7_implant_heating_fem import (  # noqa: E402
    A,
    B,
    BETA,
    FIX,
    LADDER,
    R0_READ,
    S2_BAND,
    S3_BAND,
    _ball,
    _coeffs,
    _l2_readout,
    _p2,
    _q_shell,
    _solve,
    _source_dg2,
    _theta_metal,
)
from tests.validation.test_mat18_two_layer_sphere_fem import (  # noqa: E402
    QDEG,
    _allsum,
    _dg0,
    _indicator,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_11_pec_sphere_implant_heating_in_gel_ball"
RES = LADDER[1]

#: Gate records (IMP-7 step 1, real build, 1.0 mm rung), printed beside this run.
RECORD = {
    "e_in": 1.2080e-03,   # 20261022T202402Z_IMP-7-step1-W1.log:67
    "bal": 9.099e-14,     # same log :67
    "src": 6.802147e-02,  # same log :67 (W)
    "c_iso": 2.7326e-02,  # same log :78
    "pole": 0.359626,     # same log :71, finest rung
    "equator": 0.344320,  # same log :71, finest rung
    "heating_factor": 1.1002,  # 20261022T185726Z_IMP-7-step0-W1.log:107
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _errors_c(msh, T_h, theta):
    """(whole, interior |x|^2 < 0.81 b^2, shell) rel L2 on theta; `inner`, complex-safe."""
    comm = msh.comm
    dx = ufl.dx(domain=msh, metadata=QDEG)
    d = (T_h - T_INF) - theta
    chi_in = _indicator(msh, lambda m: np.sum(m**2, axis=1) < 0.81 * B**2)
    out = [np.sqrt(max(_allsum(comm, ufl.inner(d, d) * dx), 0.0)
                   / _allsum(comm, ufl.inner(theta, theta) * dx))]
    out.append(np.sqrt(_allsum(comm, chi_in * ufl.inner(d, d) * dx)
                       / _allsum(comm, chi_in * ufl.inner(theta, theta) * dx)))
    return float(out[0]), float(out[1])


def _max_im(comm, fs):
    m = max(float(np.max(np.abs(np.imag(f.x.array)))) if f.x.array.size else 0.0 for f in fs)
    return comm.allreduce(m, op=MPI.MAX)


def _range(comm, V, T, sel):
    xs = V.tabulate_dof_coordinates()
    n = V.dofmap.index_map.size_local
    r = np.linalg.norm(xs[:n], axis=1)
    v = np.real(T.x.array[:n]) - T_INF
    m = sel(r)
    lo = comm.allreduce(float(v[m].min()) if m.any() else np.inf, op=MPI.MIN)
    hi = comm.allreduce(float(v[m].max()) if m.any() else -np.inf, op=MPI.MAX)
    return lo, hi


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
    _say(comm, "EX-140 -- the gel heated by an implant's own SAR (IMP-7 step 1, one rung)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: S2 interior {S2_BAND:.0e}, S3 balance {S3_BAND:.0e}, "
               f"Src' B1_BAND {B1_BAND:.0e}; core a = {A} m (k {K_CORE}), shell to b = {B} m "
               f"(k {K_GEL}), q_far {Q_FAR}, h {H}, T_inf {T_INF}, beta {BETA}; rung {RES} m")

    msh, ct, ft, t_build = _ball(RES)
    nglob = msh.topology.index_map(3).size_global
    _say(comm, f"[mesh] res {RES}: {nglob} cells, build {t_build:.1f} s")

    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: f"metal core a = 5 mm (k {K_CORE} W/m/K, no source)",
            2: f"gel shell to b = 20 mm (k {K_GEL} W/m/K, dipole SAR)",
        },
        translucent_tags=(2,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:11 -- metal-core sphere in a gel ball, Robin outside (IMP-7 fixture)",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- (Src') the source continuation is the gated closed form ----------------------------
    rng = np.random.default_rng(594)
    r_ = rng.uniform(1.001 * A, B, 512)
    mu = rng.uniform(-1.0, 1.0, 512)
    ph = rng.uniform(0.0, 2.0 * np.pi, 512)
    s_ = np.sqrt(1.0 - mu * mu)
    pts = np.column_stack([r_ * s_ * np.cos(ph), r_ * s_ * np.sin(ph), r_ * mu])
    q_ref = pec_sphere_implant_source(pts, radius=A, q_far=Q_FAR, beta=BETA)
    e_src = float(np.max(np.abs(_q_shell(pts.T) - q_ref)) / np.max(np.abs(q_ref)))
    _say(comm, f"\n[source] (Src') _q_shell vs pec_sphere_implant_source, 512 points: {e_src:.3e}  "
               f"ASSERTED <= B1_BAND {B1_BAND:.0e} (the gate's own check; `_q_shell` is the source used)")

    # ---- the metal-core solve -----------------------------------------------------------------
    k_f = _dg0(msh, ct, K_CORE, K_GEL)
    q_f = _source_dg2(msh, ct.find(2), _q_shell)
    sol, t_solve = _solve(msh, ft, k_f, q_f)
    theta = _p2(msh, _theta_metal)
    e_all, e_in = _errors_c(msh, sol.temperature, theta)
    src = float(np.real(sol.source_power))
    conv = float(np.real(sol.convected_power))
    bal = abs(conv - src) / src
    p_cf = float(_coeffs(K_CORE)["power"])
    _say(comm, f"\n[steady] (S2) interior rel L2 {e_in:.4e} (record {RECORD['e_in']:.4e}); "
               f"whole {e_all:.4e}  ASSERTED interior <= {S2_BAND:.0e}; solve {t_solve:.1f} s")
    _say(comm, f"[steady] (S3) balance {bal:.3e} (record {RECORD['bal']:.3e})  "
               f"ASSERTED <= {S3_BAND:.0e}")
    _say(comm, f"[source] assembled source power {src:.6e} W (gate record {RECORD['src']:.6e} W); "
               f"closed-form shell integral {p_cf:.6e} W; int_q_h/P-1 = {src / p_cf - 1:.4e} "
               f"(printed, mesh volume discretisation)")

    # ---- (c-iso) negative control -------------------------------------------------------------
    q_iso = _source_dg2(msh, ct.find(2), lambda x: _q_shell(x, drop_l2=True))
    sol_iso, _ = _solve(msh, ft, k_f, q_iso)
    _, e_iso = _errors_c(msh, sol_iso.temperature, theta)
    _say(comm, f"[control] (c-iso) isotropic source (l = 2 term dropped): interior miss {e_iso:.4e} "
               f"(record {RECORD['c_iso']:.4e}); ASSERTED > {S2_BAND:.0e}; "
               f"{e_iso / S2_BAND:.1f}x the band")

    # ---- printed only -------------------------------------------------------------------------
    c0 = _coeffs(K_CORE)
    r2_cf = float(pec_sphere_implant_radial(R0_READ, c0, 2))
    th_pole = float(_theta_metal(np.array([[0.0], [0.0], [R0_READ]]))[0])
    th_eq = float(_theta_metal(np.array([[R0_READ], [0.0], [0.0]]))[0])
    _say(comm, "")
    _l2_readout(sol.temperature, "metal", (r2_cf, th_pole, th_eq))
    _say(comm, f"[readout] pole / equator rise at r = 1.1a: closed form {th_pole:.6f} / {th_eq:.6f} K; "
               f"gate finest-rung records {RECORD['pole']:.6f} / {RECORD['equator']:.6f} K. "
               f"R2_h is a difference of two near-equal numbers (pole - equator is ~5 % of the mean), "
               f"so a sub-percent pointwise error reads as a double-digit percent on R2: a "
               f"cancellation readout, not a defect. No band.")
    V1 = sol.temperature.function_space
    cl, ch = _range(comm, V1, sol.temperature, lambda r: r <= A * (1 + 1e-9))
    sl, sh = _range(comm, V1, sol.temperature, lambda r: r >= A * (1 - 1e-9))
    _say(comm, f"[core] core theta range [{cl:.6f}, {ch:.6f}] K, spread {ch - cl:.3e} K; "
               f"shell range [{sl:.6f}, {sh:.6f}] K, spread {sh - sl:.3e} K")
    o = np.zeros((1, 3))
    fx = dict(radius=A, outer_radius=B, k_core=K_CORE, k_gel=K_GEL, q_far=Q_FAR, h=H, t_inf=T_INF)
    c1 = float(pec_sphere_implant_steady_temperature(o, beta=1.0, **fx)[0]) - T_INF
    c0b = float(pec_sphere_implant_steady_temperature(o, beta=0.0, **fx)[0]) - T_INF
    _say(comm, f"[heating] closed-form centre rise beta = 1: {c1:.5f} K; beta = 0: {c0b:.5f} K; "
               f"factor {c1 / c0b:.4f} (IMP-7 step 0 printed {RECORD['heating_factor']})")

    im = _max_im(comm, [sol.temperature, sol_iso.temperature])
    _say(comm, f"[complex] max |Im T_h| over both solves = {im:.3e}  ASSERTED == 0")

    # ---- assertions ---------------------------------------------------------------------------
    assert e_src <= B1_BAND, f"(Src') {e_src:.3e}"
    assert e_in <= S2_BAND, f"(S2) interior {e_in:.4e}"
    assert bal <= S3_BAND, f"(S3) balance {bal:.3e}"
    assert e_iso > S2_BAND, f"(c-iso) did not fire: {e_iso:.4e}"
    assert im == 0.0, f"max |Im T_h| = {im:.3e}"

    # ---- outputs ------------------------------------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    th = fem.Function(V1, name="theta_h")
    th.x.array[:] = np.real(sol.temperature.x.array) - T_INF
    th.x.scatter_forward()
    th_ex = fem.Function(V1, name="theta_exact")
    th_ex.interpolate(theta)
    th_ex.x.scatter_forward()
    err = fem.Function(V1, name="theta_abs_error")
    err.x.array[:] = np.abs(np.real(th.x.array - th_ex.x.array))
    err.x.scatter_forward()
    q_p1 = fem.Function(V1, name="q_source")
    q_p1.interpolate(q_f)
    q_p1.x.scatter_forward()
    xdmf = OUTPUT_DIR / f"{BASENAME}_combined"
    write_xdmf_with_tags(xdmf, msh, ct, {"theta_h": th, "theta_exact": th_ex,
                                          "theta_abs_error": err, "q_source": q_p1}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {xdmf}.xdmf -- P1 fields on the {nglob}-cell mesh")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
