"""Example (`EX-68`): the gel warms and conducts better — the sigma(T)
EM-thermal fixed point, iterate by iterate.

No other example feeds a temperature back into a material. This one runs
`MAT-5` step 2's eight Picard iterates on `mat:2`'s lossy sphere (R = 0.01 m,
eps_r 78, 64 MHz, imposed uniform field): EM solve with a per-cell DG0
sigma on the sphere -> G * 1/2 sigma |E|^2 -> the Robin-cooled ball
(h = 3 W/m^2/K, Bi = 0.05) -> sigma = sigma_0 [1 + alpha theta_K] -> the
next EM solve. On an imposed field the more conductive sphere absorbs *less*
(`mat:3`), so the feedback is negative and the iterates alternate,
contracting by about -0.31 per step. Freezing sigma at sigma_0 overshoots by
29 %.

**It asserts, and it does not re-implement** (the `ANS-1` rule). The helpers
(``_em_solve``, ``_heat_source``, ``_volume_mean``, ``_sub_to_parent``), the
iterate count and the bands are imported from
``tests/validation/test_mat5_em_coupled_sphere.py``; the fixture constants
from ``tests/materials/test_mat5_step2_lumped_root.py`` and the submesh tools
from `MAT-10`'s module. ``src/`` has no coupled-solver API: the loop below is
this script's own, over those helpers and public calls
(``solve_steady_conduction``, ``LinearTemperatureConductivity``,
``robin_sphere_steady_temperature``), and it is held to the gate's bands and
printed record. No gate module is edited.

Asserted at iterate 8: (A') interior theta against the profile at q* <=
``A2_BAND``; increments alternate and |delta_8|/theta_8 <= ``S_BAND``; (B')
<= ``B_BAND``; (C') <= ``C_BAND``; (D') <= ``D_BAND``. Every surface-mean
iterate theta_s must equal the gate's record to its printed six decimals —
a miss means this loop is not the gate's (a divergence finding, not a
tolerance). Negative control, asserted only as "misses ``A2_BAND``": iterate
1, the frozen-sigma_0 solve.

Run it::

    ./run_examples.sh -e mat:5 -n 2 -t 600

Output (on the submesh): ``examples/materials/paraview_output/
materials_05_sigma_of_t_fixed_point_combined.xdmf`` (iterates 1-8 as time
steps: DG0 ``sigma``, DG0 ``q_h``, P1 ``theta``) and
``..._iterates.csv``.
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
    write_xdmf_time_series,
)
from fem_em_solver.materials.temperature import LinearTemperatureConductivity  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal.steady import solve_steady_conduction  # noqa: E402
from fem_em_solver.utils.analytical import robin_sphere_steady_temperature  # noqa: E402

from tests.materials.test_mat5_step2_lumped_root import (  # noqa: E402
    ALPHA,
    G_STEP2,
    H_STEP2,
    Q_STAR,
    SIGMA_0,
    T_REF,
    THETA_STAR,
)
from tests.validation.test_lossy_sphere_sar import (  # noqa: E402
    SIGMA_HIGH,
    SPHERE_RADIUS,
    SPHERE_TAG,
    _solve_lossy_sphere,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (  # noqa: E402
    B_BAND,
    C_BAND,
    K,
    ROBIN_TAG,
    T_INF,
    _allsum,
    _interior_indicator,
    _interior_rel,
    _robin_tags,
    _sphere_submesh,
)
from tests.validation.test_mat5_em_coupled_sphere import (  # noqa: E402
    A2_BAND,
    D_BAND,
    N_PICARD,
    S_BAND,
    _em_solve,
    _heat_source,
    _real,
    _volume_mean,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_05_sigma_of_t_fixed_point"
R = SPHERE_RADIUS

#: Gate record, `20260930T124701Z_MAT-5-step2-W2.log:163-183` — theta_s per
#: iterate (K, six decimals, held to the printed digits), the last ratio and
#: the (A') / (S) / (B') / (C') / (D') / control readings (printed beside).
RECORD_THETA_S = (64.599332, 47.180063, 52.417352, 50.780364,
                  51.286944, 51.129658, 51.178444, 51.163307)
RECORD = {"ratio": -0.31027, "A": 2.3225e-02, "S": 2.959e-04, "B": 2.307e-15,
          "C": 6.344e-15, "D": 0.0, "ctrl": 2.9192e-01}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-68 -- the gel warms and conducts better: the sigma(T) fixed point")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] N = {N_PICARD} undamped iterates; bands A' {A2_BAND:.3f}, S {S_BAND:.0e}, "
               f"B' {B_BAND:.0e}, C' {C_BAND:.0e}, D' {D_BAND:.0e} (imported from the MAT-5 gate)")
    _say(comm, f"[fixture] sigma_0 = {SIGMA_0} S/m at T_ref = {T_REF}, alpha = {ALPHA} /K, "
               f"h = {H_STEP2} W/m^2/K, G = {G_STEP2:.6e}, q* = {Q_STAR:.4f} W/m^3, "
               f"theta* = {THETA_STAR} K")

    # ---- the gate's mesh, through the same builder --------------------------
    t = time.perf_counter()
    run = _solve_lossy_sphere(SIGMA_HIGH, R / 10.0, R / 5.0, return_fields=True)
    msh, cell_tags = run["mesh"], run["cell_tags"]
    t_mesh = time.perf_counter() - t

    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={SPHERE_TAG: "sphere (gel, sigma(T))", 2: "air"},
        hide_tags=(2,),
        translucent_tags=(SPHERE_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mat:5 -- sigma(T) EM-thermal fixed point on the lossy gel ball (MAT-5 step 2 fixture), h = 3",
        comm=comm,
    )

    sub, parent_ids = _sphere_submesh(msh, cell_tags)
    ft = _robin_tags(sub)
    chi = _interior_indicator(sub)
    n_parent = msh.topology.index_map(3).size_global
    n_sub = sub.topology.index_map(3).size_global
    _say(comm, f"\n[mesh] parent {n_parent} cells, sphere submesh {n_sub} cells ({t_mesh:.1f} s incl. "
               "the builder's reference solve)")

    Ws = fem.functionspace(sub, ("DG", 0))
    law = LinearTemperatureConductivity(sigma_ref=SIGMA_0, alpha_per_k=ALPHA, t_ref=T_REF)
    x = ufl.SpatialCoordinate(sub)
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    area = _real(_allsum(comm, fem.Constant(sub, default_scalar_type(1.0)) * ds))

    th_c = float(robin_sphere_steady_temperature(
        0.0, radius=R, k=K, q=Q_STAR, h=H_STEP2, t_inf=T_INF)) - T_INF
    th_R = float(robin_sphere_steady_temperature(
        R, radius=R, k=K, q=Q_STAR, h=H_STEP2, t_inf=T_INF)) - T_INF
    theta_cf = th_c - (th_c - th_R) / R**2 * ufl.inner(x, x)
    _say(comm, f"[profile] at q*: theta(0) = {th_c:.6f} K, theta(R) = {th_R:.6f} K")

    # ---- Picard: theta_0 = 0, N fixed, no relaxation ---------------------------
    T_prev = None
    th_s_hist = [0.0]
    err_hist, sbar_hist, frames = [], [], []
    last = {}
    _say(comm, "\n[picard] n  sigma_bar_in [S/m]  theta_s [K]   record       delta [K]     "
               "interior rel L2 vs profile")
    for n in range(1, N_PICARD + 1):
        sig_sub = fem.Function(Ws, name="sigma")
        if T_prev is None:
            t_cell = np.full(sig_sub.x.array.shape, T_INF, dtype=float)
        else:
            t_dg0 = fem.Function(Ws)
            t_dg0.interpolate(T_prev)
            t_dg0.x.scatter_forward()
            t_cell = np.real(t_dg0.x.array)
        sig_sub.x.array[:] = law.sigma(t_cell)
        sig_sub.x.scatter_forward()

        E_n, sigf_n, sbar_n, p_n = _em_solve(msh, cell_tags, sub, sig_sub, parent_ids, f"ex68it{n}")
        q_n = _heat_source(msh, cell_tags, sub, parent_ids, E_n, sigf_n, G_STEP2)
        q_n.name = "q_h"
        res_n = solve_steady_conduction(sub, ft, k=K, q=q_n, h=H_STEP2, t_inf=T_INF,
                                        robin_tag=ROBIN_TAG, degree=1)
        th_s = _real(_allsum(comm, (res_n.temperature - T_INF) * ds)) / area
        err_n = _interior_rel(sub, chi, (res_n.temperature - T_INF) - theta_cf, theta_cf)
        th_s_hist.append(th_s)
        err_hist.append(err_n)
        sbar_hist.append(sbar_n)
        # Copies before the next solve overwrites anything.
        theta_n = fem.Function(res_n.temperature.function_space, name="theta")
        theta_n.x.array[:] = res_n.temperature.x.array - T_INF
        theta_n.x.scatter_forward()
        frames.append((float(n), {"sigma": sig_sub, "q_h": q_n, "theta": theta_n}))
        _say(comm, f"[picard] {n}  {sbar_n:.9f}       {th_s:.6f}   {RECORD_THETA_S[n - 1]:.6f}   "
                   f"{th_s - th_s_hist[-2]:+.6e}  {err_n:.4e}")
        T_prev = res_n.temperature
        last = dict(sigf=sigf_n, sig_sub=sig_sub, p=p_n, q=q_n, res=res_n)

    deltas = [th_s_hist[i] - th_s_hist[i - 1] for i in range(1, N_PICARD + 1)]
    ratios = [deltas[i] / deltas[i - 1] for i in range(1, N_PICARD)]
    alternate = all(deltas[i] * deltas[i - 1] < 0.0 for i in range(1, N_PICARD))
    s_rel = abs(deltas[-1]) / abs(th_s_hist[-1])
    _say(comm, "[picard] delta_{n+1}/delta_n = " + "  ".join(f"{v:+.5f}" for v in ratios)
               + f"  (record last {RECORD['ratio']:+.5f}; printed)")
    _say(comm, f"(S) alternate = {alternate}, |delta_8|/theta_8 = {s_rel:.3e} (record "
               f"{RECORD['S']:.3e})  ASSERTED <= {S_BAND:.0e}")

    err_A = err_hist[-1]
    _say(comm, f"(A') theta_8 interior rel L2 = {err_A:.4e} (record {RECORD['A']:.4e})  "
               f"ASSERTED <= {A2_BAND:.3f}")

    res8 = last["res"]
    int_q = _real(_allsum(comm, last["q"] * ufl.dx(domain=sub)))
    p_par = G_STEP2 * last["p"]
    err_B = abs(int_q - p_par) / abs(p_par)
    err_C = abs(res8.convected_power - res8.source_power) / abs(res8.source_power)
    xp = ufl.SpatialCoordinate(msh)
    dx1 = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)(SPHERE_TAG)
    m_par = _real(_allsum(comm, last["sigf"] * ufl.inner(xp, xp) * dx1))
    m_sub = _real(_allsum(comm, last["sig_sub"] * ufl.inner(x, x) * ufl.dx(domain=sub)))
    err_D = abs(m_par - m_sub) / abs(m_sub)
    _say(comm, f"(B') transfer rel = {err_B:.3e} (record {RECORD['B']:.3e})  ASSERTED <= {B_BAND:.0e}")
    _say(comm, f"(C') balance rel = {err_C:.3e} (record {RECORD['C']:.3e})  ASSERTED <= {C_BAND:.0e}")
    _say(comm, f"(D') sigma x.x parent vs sub rel = {err_D:.3e} (record {RECORD['D']:.1e})  "
               f"ASSERTED <= {D_BAND:.0e}")
    ctrl = err_hist[0]
    _say(comm, f"[control] iterate 1 (frozen sigma_0): {ctrl:.4e} (record {RECORD['ctrl']:.4e})  "
               f"ASSERTED > {A2_BAND:.3f}")
    q_bar, _vol = _volume_mean(sub, last["q"])
    _say(comm, f"[record] theta_s,8 = {th_s_hist[-1]:.6f} K vs the lumped root theta_L = {THETA_STAR} K "
               f"({th_s_hist[-1] / THETA_STAR - 1.0:+.4%}); mean q_h,8 = {q_bar:.4f} W/m^3 vs q* = "
               f"{Q_STAR:.4f} (printed)")

    # ---- the loop is the gate's: every iterate on the printed record -----------
    diverged = [n for n, (v, rec) in enumerate(zip(th_s_hist[1:], RECORD_THETA_S), start=1)
                if f"{v:.6f}" != f"{rec:.6f}"]
    _say(comm, f"[record] theta_s iterates off the gate's printed digits: {diverged or 'none'}")
    assert not diverged, (
        f"DIVERGENCE: theta_s at iterates {diverged} differs from the gate record beyond its "
        "printed digits -- this loop is not the gate's"
    )

    assert err_A <= A2_BAND, f"(A') {err_A:.4e}"
    assert alternate, f"(S) increments do not alternate: {deltas}"
    assert s_rel <= S_BAND, f"(S) {s_rel:.3e}"
    assert err_B <= B_BAND, f"(B') {err_B:.3e}"
    assert err_C <= C_BAND, f"(C') {err_C:.3e}"
    assert err_D <= D_BAND, f"(D') {err_D:.3e}"
    assert ctrl > A2_BAND, f"control {ctrl:.4e} did not miss {A2_BAND}"

    # ---- outputs (on the submesh) ------------------------------------------------
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_time_series(
        OUTPUT_DIR / f"{BASENAME}_combined", sub, None, frames, comm=comm,
    )
    table = OUTPUT_DIR / f"{BASENAME}_iterates.csv"
    if comm.rank == 0:
        with table.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["iterate", "sigma_bar_in_S_per_m", "theta_s_K", "delta_K",
                        "interior_rel_l2_vs_profile"])
            for n in range(1, N_PICARD + 1):
                w.writerow([n, f"{sbar_hist[n - 1]:.9f}", f"{th_s_hist[n]:.9f}",
                            f"{deltas[n - 1]:+.9e}", f"{err_hist[n - 1]:.6e}"])
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- iterates 1..{N_PICARD} as time steps: DG0 sigma, DG0 q_h, "
               "P1 theta (theta = T - T_inf)")
    _say(comm, f"[output] {table}")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
