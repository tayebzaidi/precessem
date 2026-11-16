"""OPS-83 (§9 item 859): the (S2-tip) comparand's form on a rough DG0 source.

On I32a's 12³ stub (K10es's set-up, `test_geo36_rod_f_clinical_heating.py:656–666`)
one backward-Euler step from T∞ is taken on a uniform DG0 source q_u (G = 1/P_gel)
and on a cube-scale checkerboard q_r = q_u·(1 + ½·s), at Δt₁ and Δt₁/10.  The 1 g
ball means of θ¹ are set against ⟨q⟩_B (the gate as written), ⟨Πq⟩_B (the
consistent-mass L2 projection onto the P1 temperature space) and ⟨Π_L q⟩_B (the
lumped projection).

Asserted (bands imported): (E) on all four solves; the global (S2) at `S2_BAND` on
all four; K10-S's printed tip 8.0400e-06 reproduced on (u) at Δt₁ within 5e-11;
the projection identity ⟨Πq_u⟩_B / ⟨q_u⟩_B = 1 within 1e-10 at both balls (LU
round-off ceiling, not a physics band).

Printed only (rule (e)): R_DG0, R_Π, R_L, their Δt-ratios, and the pre-stated
reading line (form CONFIRMED / REFUTED / INCONCLUSIVE).  Not claimed: a re-pose of
the tip gate, a band, a clinical-scale number.
"""

from __future__ import annotations

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI

from fem_em_solver.post.sar import mean_sar
from tests.complex_mode import complex_only
from tests.validation.test_geo36_rod_f_clinical_drive import _k9_fixture
from tests.validation.test_geo36_rod_f_clinical_heating import (
    DT_1_BOX,
    E_BAND,
    ONE_GRAM_KG,
    PHANTOM_CELL_TAG,
    PHANTOM_RHO_KG_PER_M3,
    QUADRATURE_DEGREE,
    RHO_C,
    ROBIN_TAG,
    S2_BAND,
    STUB_CENTRE,
    STUB_PROBES,
    T_INF,
    _allsum,
    _ball_reader,
    _columns_e_sup,
    _heat_source,
    _merge_rod_into_phantom,
    _outer_robin_tags,
    _say,
    _tag_submesh,
    _transient,
    averaging_ball_radius,
)

H0 = 0.1 / 12  # the stub's create_box cube size (`test_geo36_rod_f_clinical_drive.py:1188`)
K10S_TIP = 8.0400e-06  # `20261101T103929Z_GEO-36-I32a-W0b.log:76`
K10S_TIP_TOL = 5.0e-11  # its last printed digit
PROJECTION_ID_TOL = 1.0e-10  # LU round-off ceiling, not a physics band
CONST_RATIO = (0.5, 2.0)  # "Δt-independent" in the pre-stated reading
SCALING_RATIO = 5.0  # "an O(Δt) residual" in the pre-stated reading


def _checkerboard(q_u):
    """q_u·(1 + ½·s), s = sign(sin(πx/h₀)sin(πy/h₀)sin(πz/h₀)) at each cell midpoint."""
    V0 = q_u.function_space
    sub = V0.mesh
    tdim = sub.topology.dim
    imap = sub.topology.index_map(tdim)
    cells = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
    mids = dmesh.compute_midpoints(sub, tdim, cells)
    s = np.sign(np.sin(np.pi * mids[:, 0] / H0) * np.sin(np.pi * mids[:, 1] / H0)
                * np.sin(np.pi * mids[:, 2] / H0))
    dofs = np.asarray(V0.dofmap.list)[cells, 0]
    q_r = fem.Function(V0)
    q_r.x.array[:] = q_u.x.array
    q_r.x.array[dofs] = q_u.x.array[dofs] * (1.0 + 0.5 * s)
    q_r.x.scatter_forward()
    n_zero = sub.comm.allreduce(int(np.count_nonzero(s[: imap.size_local] == 0.0)), op=MPI.SUM)
    return q_r, n_zero


def _owned_min_max(f):
    comm = f.function_space.mesh.comm
    n = f.function_space.dofmap.index_map.size_local
    v = np.real(f.x.array[:n])
    return (comm.allreduce(float(v.min()) if n else np.inf, op=MPI.MIN),
            comm.allreduce(float(v.max()) if n else -np.inf, op=MPI.MAX))


def _projections(q, V1, tag):
    """Πq (consistent mass, direct LU) and Π_L q = b / (M·1) on V1."""
    u, v = ufl.TrialFunction(V1), ufl.TestFunction(V1)
    dx = ufl.dx(domain=V1.mesh)
    problem = LinearProblem(
        ufl.inner(u, v) * dx, ufl.inner(q, v) * dx, bcs=[],
        petsc_options_prefix=f"fem_em_ops83_proj_{tag}_",
        petsc_options={"ksp_type": "preonly", "pc_type": "lu",
                       "pc_factor_mat_solver_type": "mumps"},
    )
    pq = problem.solve()
    pq.x.scatter_forward()
    one = fem.Constant(V1.mesh, default_scalar_type(1.0))
    b = fem.assemble_vector(fem.form(ufl.inner(q, v) * dx))
    m1 = fem.assemble_vector(fem.form(ufl.inner(one, v) * dx))
    from dolfinx import la

    b.scatter_reverse(la.InsertMode.add)
    m1.scatter_reverse(la.InsertMode.add)
    n = V1.dofmap.index_map.size_local
    plq = fem.Function(V1)
    plq.x.array[:n] = b.array[:n] / m1.array[:n]
    plq.x.scatter_forward()
    return pq, plq


@complex_only
def test_ops83_s2tip_form_on_stub(monkeypatch, tmp_path):
    comm = MPI.COMM_WORLD
    fx = _k9_fixture(comm)
    msh, cell_tags = fx[0], fx[1]
    root = comm.bcast(str(tmp_path), root=0)
    a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    v_ball = 4.0 / 3.0 * np.pi * a_1g ** 3
    balls = {"tip": np.asarray(fx[4]["a"], dtype=float),
             "interior": np.asarray(STUB_PROBES["interior"], dtype=float)}
    dts = {"dt1": DT_1_BOX, "dt1/10": DT_1_BOX / 10.0}
    R, E, S2 = {}, {}, {}
    with monkeypatch.context() as m:
        e_sup, sigma_f, _ = _columns_e_sup(m, fx, root)
        merged = _merge_rod_into_phantom(msh, cell_tags)
        P = float(np.real(mean_sar(e_sup, sigma=sigma_f, rho=PHANTOM_RHO_KG_PER_M3,
                                   cell_tags=merged, comm=comm,
                                   subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"]))
        G = 1.0 / P
        sub, parent_ids = _tag_submesh(msh, merged, PHANTOM_CELL_TAG)
        ft, counts = _outer_robin_tags(sub, STUB_CENTRE)
        one = fem.Constant(sub, default_scalar_type(1.0))
        v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
        q_u = _heat_source(msh, merged, e_sup, sigma_f, G, sub, parent_ids)
        q_r, n_zero = _checkerboard(q_u)
        _say(comm, f"[OPS-83] P_gel {P:.12e} W, G {G:.12e} 1/W, V_sub {v_sub:.12e} m^3, facets "
                   f"{counts} (Robin tag {ROBIN_TAG}); ball quadrature degree {QUADRATURE_DEGREE}; a_1g {a_1g * 1e3:.4f} mm, 4/3 pi a^3 {v_ball:.6e} m^3")
        for name, q in (("u", q_u), ("r", q_r)):
            lo, hi = _owned_min_max(q)
            _say(comm, f"[OPS-83] q_{name} owned min {lo:.12e} max {hi:.12e} W/m^3"
                       + (f"; zero-sign midpoints {n_zero}" if name == "r" else ""))
        V1 = None
        proj = {}
        for name, q in (("u", q_u), ("r", q_r)):
            int_q = float(np.real(_allsum(comm, q * ufl.dx(domain=sub))))
            qbar = int_q / v_sub
            for dname, dt in dts.items():
                res = _transient(sub, ft, q, dt=dt, n_steps=1, theta=1.0)
                if V1 is None:
                    V1 = res.temperature.function_space
                if name not in proj:
                    proj[name] = _projections(q, V1, name)
                pq, plq = proj[name]
                th1 = fem.Function(V1)
                th1.x.array[:] = res.temperature.x.array - T_INF
                th1.x.scatter_forward()
                thbar = float(np.real(_allsum(comm, th1 * ufl.dx(domain=sub)))) / v_sub
                S2[(name, dname)] = abs(RHO_C * thbar / (dt * qbar) - 1.0)
                E[(name, dname)] = float(res.max_balance_residual)
                read = _ball_reader(sub, {"th1": th1, "q": q, "pq": pq, "plq": plq})
                for bname, c in balls.items():
                    vol, mb = read(c, a_1g)
                    lhs = RHO_C * mb["th1"] / dt
                    R[(name, dname, bname)] = {
                        "DG0": abs(lhs / mb["q"] - 1.0), "Pi": abs(lhs / mb["pq"] - 1.0),
                        "L": abs(lhs / mb["plq"] - 1.0), "pq/q": abs(mb["pq"] / mb["q"] - 1.0),
                        "vol": vol,
                    }
                    r = R[(name, dname, bname)]
                    _say(comm, f"[OPS-83] q_{name} {dname} ({dt:.6e} s) ball {bname} {tuple(c)}: vol "
                               f"{vol:.6e} m^3 ({vol / v_ball:.6f} x 4/3 pi a^3); <q> {mb['q']:.12e} "
                               f"<Pq> {mb['pq']:.12e} <PLq> {mb['plq']:.12e}; R_DG0 {r['DG0']:.4e} "
                               f"R_Pi {r['Pi']:.4e} R_L {r['L']:.4e}; |<Pq>/<q> - 1| {r['pq/q']:.3e}")
                _say(comm, f"[OPS-83] q_{name} {dname}: (E) {E[(name, dname)]:.3e} (band {E_BAND:.0e}, "
                           f"ASSERTED); global (S2) {S2[(name, dname)]:.4e} (band {S2_BAND:.0e}, "
                           f"ASSERTED)")
    ratio = {}
    _say(comm, "[OPS-83] table: source ball | R_DG0(dt1) R_DG0(dt1/10) ratio | R_Pi(dt1) "
               "R_Pi(dt1/10) ratio | R_L(dt1) R_L(dt1/10) ratio")
    for name in ("u", "r"):
        for bname in balls:
            a, b = R[(name, "dt1", bname)], R[(name, "dt1/10", bname)]
            row = []
            for k in ("DG0", "Pi", "L"):
                ratio[(name, bname, k)] = a[k] / b[k] if b[k] > 0 else np.inf
                row.append(f"{a[k]:.4e} {b[k]:.4e} {ratio[(name, bname, k)]:.3f}")
            _say(comm, f"[OPS-83] table: q_{name} {bname} | " + " | ".join(row))
    # The pre-stated reading (rule (e): printed, never asserted).
    rough = [R[("r", d, bn)] for d in dts for bn in balls]
    refuted = any(r["Pi"] > S2_BAND for r in rough)
    confirmed = (not refuted
                 and all(R[("r", "dt1", bn)]["DG0"] > S2_BAND for bn in balls)
                 and all(CONST_RATIO[0] <= ratio[("r", bn, "DG0")] <= CONST_RATIO[1] for bn in balls)
                 and all(ratio[("r", bn, "Pi")] >= SCALING_RATIO for bn in balls))
    reading = "form REFUTED" if refuted else ("form CONFIRMED" if confirmed else "INCONCLUSIVE")
    _say(comm, f"[OPS-83] READING (printed, not asserted): {reading}")

    for key in E:
        assert E[key] <= E_BAND, (key, E[key])
        assert S2[key] <= S2_BAND, (key, S2[key])
    tip = R[("u", "dt1", "tip")]["DG0"]
    _say(comm, f"[OPS-83] reproduction control: tip R_DG0(u, dt1) {tip:.6e} vs K10-S {K10S_TIP:.4e}: "
               f"|diff| {abs(tip - K10S_TIP):.3e} (tol {K10S_TIP_TOL:.0e}, ASSERTED); projection "
               f"identity tip {R[('u', 'dt1', 'tip')]['pq/q']:.3e} interior "
               f"{R[('u', 'dt1', 'interior')]['pq/q']:.3e} (tol {PROJECTION_ID_TOL:.0e}, ASSERTED)")
    assert abs(tip - K10S_TIP) <= K10S_TIP_TOL, tip
    for bname in balls:
        assert R[("u", "dt1", bname)]["pq/q"] <= PROJECTION_ID_TOL, (bname, R[("u", "dt1", bname)])
