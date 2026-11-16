"""`MAT-5` step 3 — σ(T) fed back on P1's coil-driven fixture (§9 item 237, chain P6).

Fixture: `MAT-15` step 1's (P1) — `ANS-2`'s 10 MHz loaded F-small birdcage on
the 0.0025 m phantom rung (``build_four_port_sweep(…, build_only=True)``),
drive k = 0 (``P1``) through ``_solve_driven``, the σ = 800 S/m wire, saline
σ₀ = ``SALINE_SIGMA`` 0.5 S/m on tag 3.  Thermal: the tag-3 submesh, Robin on
every exterior facet, gel k 0.6, h 10, T∞ 20 (Bi = 0.5), steady P1 solves.

Map (steady Picard, θ = T − 20): Φ(θ) = S[G·½σ(θ)|E(σ(θ))|²], σ(θ) =
σ₀[1 + αΠ₀θ] on tag 3 only (Π₀ = DG0 interpolant on the sub), written into the
parent DG0 σ through ``TimeHarmonicSolver.solve(sigma_field_hook=…)`` and
``_sub_to_parent`` (cell → dof through both dofmaps).  E is re-solved on every
iterate.  G = 1 W / P₁ is fixed once from iterate 1 (σ ≡ σ₀) and never
renormalised.  θ₀ = 0, N = 8, no relaxation.

Anchors (asserted): (R0a) hook path = P1 path in-run ≤ ``L0_BAND``; (R0b)
α = 0 reproduces P1's records ≤ ``REPRO_BAND`` (mean θ(t*) at ½ ulp 5e-7 K,
V_sub at ``V_BAND``); (R0c) α = 0 map residual ≤ ``L0_BAND`` 1e-9 (§9 item 239,
provisional: item 237 measured 2.747e-12 against the old 1e-12, an in-process
repeat at -n 2 is not bitwise; printed beside it, never asserted: max|σ_it2 −
σ_it1| over all local parent DG0 dofs and ‖E_it2 − E_it1‖/‖E_it1‖ on owned
dofs, the pathologist's separation test); (D0) hook touches only
tag 3 (exact 0); (K) |ρ₄/μ_pred − 1| ≤ 0.10, μ_pred the FD power iteration of
the re-solved-E Jacobian at θ₇ (v₀ ≡ 1, M = 5, εv RMS 1e-2 K); (K-conv)
|μ₅/μ₄ − 1| ≤ 1e-2; (M0) |μ_pred| ≥ 5e-3; (F) ‖δ₈‖/‖θ₈‖ ≤ ``S_BAND``; (E*),
(B′), (C′), (D′) at iterate 8.  Printed: ρ₁…ρ₇, μ₁…μ₅, μ_E (frozen-|E|,
control N2), the α = 0 steady state (θ̄, θ_s, max θ), θ̄₈/θ̄₁ − 1, σ̄₈,
P_diss,8/P₁, EM s/solve, peak RSS.

Scope: one drive, 10 MHz, degree-1 EM, P1 thermal, one rung, α an input,
G = 1 W.  Not: transient coupling, ε(T), perfusion, tissue, h-convergence, an
absolute temperature, the external check (`ANS-9` step 2).
"""

from __future__ import annotations

import resource
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.materials.temperature import LinearTemperatureConductivity
from fem_em_solver.post.sar import mean_sar
from fem_em_solver.thermal.steady import solve_steady_conduction

from tests.complex_mode import complex_only
from tests.materials.test_mat5_step2_lumped_root import ALPHA
from tests.validation.test_birdcage_b1_plus_map import PHANTOM_RHO_KG_PER_M3, _solve_driven
from tests.validation.test_birdcage_sar_1g_rung import PHANTOM_RESOLUTION_1G_RUNG
from tests.validation.test_lossy_sphere_fullwave import SALINE_SIGMA
from tests.validation.test_mat5_em_coupled_sphere import D_BAND, L0_BAND, S_BAND, _sub_to_parent
from tests.validation.test_mat10_sar_to_heat_sphere import (
    B_BAND,
    C_BAND,
    ROBIN_TAG,
    _allsum,
    _robin_tags,
    _say,
    _to_sub,
)
from tests.validation.test_mat12_transient_sar_heating_sphere import E_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    H_GEL,
    K_GEL,
    P1_RECORD_W,
    T_INF,
    V_BAND,
    V_SUB_RECORD_M3,
    _cell_average_dg0_tag,
    _cn,
    _e_prime,
    _tag_submesh,
)
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

# P1's (E′) split on drive k = 0, CN N = 40 (`20261001T123309Z_MAT-15-step1-W1.log:1930`;
# two other runs read 7.657997570404 / 405e-01, ≤ 4e-13).
P1_E_STORED_W = 7.657997570402e-01
P1_E_SURFACE_W = 2.342002429598e-01
P1_MEAN_THETA_K = 1.153314  # `:1929`, printed to 6 decimals
MEAN_THETA_HALF_ULP_K = 5.0e-7  # ½ ulp of the print, never tighter
# Relative.  ≥ ½ ulp of every 13-digit print (≤ 1e-13) and ≈ 190× the measured
# cross-run spread of P₁ over six runs (≤ 5.3e-12); complex MUMPS at -n 2 is
# not bit-reproducible across processes (OPS-43), so no cross-run compare is bitwise.
REPRO_BAND = 1.0e-9
# §9 item 239 (provisional, the weekly ratifies or reverts): step 2's in-run
# α = 0 band.  Item 237 measured 2.747e-12 against the old 1e-12, which assumed a
# bitwise in-process repeat that -n 2 does not give ((R0a): P 1.7e-12 apart;
# `20261005T050759Z_MAT-5-step3-W1-alpha0.log`).
R0C_BAND = L0_BAND
K_BAND = 0.10  # §10: within 10 %
KCONV_BAND = 1.0e-2
M0_FLOOR = 5.0e-3
N_PICARD = 8
M_POWER = 5
FD_RMS_K = 1.0e-2  # εv RMS; 1e-6 K is solver noise, 1 K is nonlinear

_ = E_BAND  # imported per the brief; the steady runs report no per-step balance


def _real(v):
    return float(np.real(v))


def _rss_summed_gib(comm):
    """Peak RSS summed over ranks; every rank calls (a comm call under rank 0 deadlocks)."""
    kib = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(float(kib), op=MPI.SUM) / 1024.0**2


class _Fixture:
    """Mesh, submesh and the Φ-map pieces shared by both windows."""

    def __init__(self, alpha):
        self.comm = MPI.COMM_WORLD
        t0 = time.perf_counter()
        self.sweep = build_four_port_sweep(phantom_resolution=PHANTOM_RESOLUTION_1G_RUNG,
                                           build_only=True)
        self.t_mesh = time.perf_counter() - t0
        self.msh, self.cell_tags = self.sweep["mesh"], self.sweep["cell_tags"]
        self.sub, self.parent_ids = _tag_submesh(self.msh, self.cell_tags, PHANTOM_CELL_TAG)
        self.ft = _robin_tags(self.sub)
        self.ds = ufl.Measure("ds", domain=self.sub, subdomain_data=self.ft)(ROBIN_TAG)
        self.dx = ufl.dx(domain=self.sub)
        one = fem.Constant(self.sub, default_scalar_type(1.0))
        self.v_sub = _real(_allsum(self.comm, one * self.dx))
        self.area = _real(_allsum(self.comm, one * self.ds))
        self.Ws = fem.functionspace(self.sub, ("DG", 0))
        self.V1 = fem.functionspace(self.sub, ("Lagrange", 1))
        self.law = LinearTemperatureConductivity(sigma_ref=SALINE_SIGMA, alpha_per_k=alpha,
                                                 t_ref=T_INF)
        self.alpha = alpha
        self.G = None
        self.em_secs = []
        _say(f"\nMAT5s3 alpha = {alpha}  mesh {int(self.sweep['cells'])} cells, sub "
             f"{self.sub.topology.index_map(3).size_global}; mesh+rung {self.t_mesh:.1f} s; "
             f"V_sub = {self.v_sub:.12e} m^3  oint ds = {self.area:.9e} m^2")

    # -- sub-mesh algebra (all reductions real, allreduced) --------------------
    def ip(self, a, b):
        return _real(_allsum(self.comm, ufl.inner(a, b) * self.dx))

    def norm(self, a):
        return float(np.sqrt(self.ip(a, a)))

    def p1(self, arr):
        f = fem.Function(self.V1)
        f.x.array[:] = arr
        f.x.scatter_forward()
        return f

    def sigma_sub(self, T):
        """σ(Π₀T) on the sub (DG0); ``T is None`` is θ = 0."""
        s = fem.Function(self.Ws)
        if T is None:
            t_cell = np.full(s.x.array.shape, T_INF, dtype=float)
        else:
            t_dg0 = fem.Function(self.Ws)
            t_dg0.interpolate(T)
            t_dg0.x.scatter_forward()
            t_cell = np.real(t_dg0.x.array)
        s.x.array[:] = self.law.sigma(t_cell)
        s.x.scatter_forward()
        return s

    # -- EM ---------------------------------------------------------------------
    def em(self, sig_sub=None):
        hook = None
        if sig_sub is not None:
            def hook(sf, _s=sig_sub):
                _sub_to_parent(sf, _s, self.parent_ids)
        t1 = time.perf_counter()
        solved = _solve_driven(self.sweep, "P1", sigma_field_hook=hook)
        self.em_secs.append(time.perf_counter() - t1)
        f = solved["fields"]
        P = float(mean_sar(f.e_complex, sigma=f.sigma_field, rho=PHANTOM_RHO_KG_PER_M3,
                           cell_tags=self.cell_tags, comm=self.comm,
                           subdomain_ids=PHANTOM_CELL_TAG)["dissipated_power_w"])
        return f, P

    def q_h(self, f):
        Wp, q_par = _cell_average_dg0_tag(
            self.msh, self.cell_tags, PHANTOM_CELL_TAG,
            self.G * 0.5 * f.sigma_field * ufl.inner(f.e_complex, f.e_complex))
        return _to_sub(Wp, q_par, self.sub, self.parent_ids)

    def steady(self, q):
        return solve_steady_conduction(self.sub, self.ft, k=K_GEL, q=q, h=H_GEL, t_inf=T_INF,
                                       robin_tag=ROBIN_TAG, degree=1)

    def phi(self, T):
        """One Picard map: σ(T) → E → q → steady T.  Sets G on first use."""
        sig = self.sigma_sub(T)
        f, P = self.em(sig)
        if self.G is None:
            self.G = 1.0 / P
        q = self.q_h(f)
        res = self.steady(q)
        Tn = self.p1(res.temperature.x.array)
        return dict(T=Tn, f=f, P=P, q=q, res=res, sig=sig)

    def theta(self, T):
        return self.p1(T.x.array - T_INF)

    def records(self, T, tag):
        th = self.theta(T)
        n_own = self.V1.dofmap.index_map.size_local
        loc = float(np.max(np.real(th.x.array[:n_own]))) if n_own else -np.inf
        tmax = self.comm.allreduce(loc, op=MPI.MAX)
        tbar = _real(_allsum(self.comm, th * self.dx)) / self.v_sub
        ts = _real(_allsum(self.comm, th * self.ds)) / self.area
        _say(f"MAT5s3 record {tag}: theta_bar = {tbar:.9f} K  theta_s = {ts:.9f} K  "
             f"max theta = {tmax:.9f} K")
        return tbar, ts, tmax


def _owned_cell_dofs(fx):
    """Parent DG0 owned-cell dofs split into tag 3 / not tag 3."""
    msh, cell_tags = fx.msh, fx.cell_tags
    n_own = msh.topology.index_map(3).size_local
    tags = np.full(n_own, -1, dtype=np.int32)
    keep = cell_tags.indices < n_own
    tags[cell_tags.indices[keep]] = cell_tags.values[keep]
    cells = np.arange(n_own, dtype=np.int32)
    return cells[tags == PHANTOM_CELL_TAG], cells[tags != PHANTOM_CELL_TAG]


@pytest.fixture(scope="module")
def alpha0():
    fx = _Fixture(alpha=0.0)
    comm = fx.comm
    out = {}

    # (R0a) P1's path (no hook) against the hook writing σ ≡ σ₀.
    f_p1, P_p1 = fx.em(None)
    fx.G = 1.0 / P_p1
    q_p1 = fx.q_h(f_p1)
    it1 = fx.phi(None)  # θ₀ = 0, hooked σ ≡ σ₀, G already fixed from P1's path
    out["R0a_P"] = abs(it1["P"] / P_p1 - 1.0)
    diff = fem.Function(fx.Ws)
    diff.x.array[:] = it1["q"].x.array - q_p1.x.array
    out["R0a_q"] = fx.norm(diff) / fx.norm(q_p1)
    _say(f"MAT5s3 (R0a) P_hook = {it1['P']:.15e}  P_P1path = {P_p1:.15e}  "
         f"|ratio-1| = {out['R0a_P']:.3e};  ||q_hook - q||/||q|| = {out['R0a_q']:.3e}  "
         f"(bound {L0_BAND:.0e}, predicted 0.0)")

    # (D0) non-tag-3 owned cells unchanged; tag-3 cells printed.
    c3, cn3 = _owned_cell_dofs(fx)
    dm = f_p1.sigma_field.function_space.dofmap.list
    d_n3 = dm[cn3, 0]
    d_3 = dm[c3, 0]
    a_b, a_h = f_p1.sigma_field.x.array, it1["f"].sigma_field.x.array
    loc = float(np.max(np.abs(a_h[d_n3] - a_b[d_n3]))) if d_n3.size else 0.0
    out["D0"] = comm.allreduce(loc, op=MPI.MAX)
    loc3 = float(np.max(np.abs(a_h[d_3] - SALINE_SIGMA))) if d_3.size else 0.0
    d0_3 = comm.allreduce(loc3, op=MPI.MAX)
    loc_hi = float(np.max(np.real(a_b[d_n3]))) if d_n3.size else 0.0
    hi = comm.allreduce(loc_hi, op=MPI.MAX)
    _say(f"MAT5s3 (D0) max|sigma_hooked - sigma_built| over non-tag-3 owned cells = "
         f"{out['D0']:.3e} (exact 0 required; max non-tag-3 sigma {hi:.1f} S/m);  "
         f"tag 3 max|sigma_hooked - sigma0| = {d0_3:.3e} (printed)")

    # (R0b) records.
    out["R0b_P"] = abs(P_p1 / P1_RECORD_W - 1.0)
    out["R0b_V"] = abs(fx.v_sub / V_SUB_RECORD_M3 - 1.0)
    t1 = time.perf_counter()
    res_cn = _cn(fx.sub, fx.ft, it1["q"])
    t_cn = time.perf_counter() - t1
    ep, stor, conv = _e_prime(fx.sub, fx.ds, res_cn, fx.G * it1["P"])
    out["R0b_stor"] = abs(stor / P1_E_STORED_W - 1.0)
    out["R0b_surf"] = abs(conv / P1_E_SURFACE_W - 1.0)
    mean_th = _real(_allsum(comm, (res_cn.temperature - T_INF) * fx.dx)) / fx.v_sub
    out["R0b_theta"] = abs(mean_th - P1_MEAN_THETA_K)
    _say(f"MAT5s3 (R0b) P_1 = {P_p1:.12e} vs record {P1_RECORD_W:.12e}: rel {out['R0b_P']:.3e}  "
         f"(bound {REPRO_BAND:.0e})")
    _say(f"MAT5s3 (R0b) CN N=40 stored {stor:.12e} W vs {P1_E_STORED_W:.12e}: rel "
         f"{out['R0b_stor']:.3e};  surface {conv:.12e} W vs {P1_E_SURFACE_W:.12e}: rel "
         f"{out['R0b_surf']:.3e}  (bound {REPRO_BAND:.0e}); (E') {ep:.3e}; CN {t_cn:.1f} s")
    _say(f"MAT5s3 (R0b) mean theta(t*) = {mean_th:.9f} K vs {P1_MEAN_THETA_K:.6f}: "
         f"|diff| = {out['R0b_theta']:.3e} K (bound {MEAN_THETA_HALF_ULP_K:.0e} K);  "
         f"V_sub rel {out['R0b_V']:.3e} (bound {V_BAND:.0e})")

    # (R0c) two α = 0 iterates.
    it2 = fx.phi(it1["T"])
    th1 = fx.theta(it1["T"])
    d = fx.p1(it2["T"].x.array - it1["T"].x.array)
    out["R0c"] = fx.norm(d) / fx.norm(th1)
    _say(f"MAT5s3 (R0c) ||theta_2 - theta_1||/||theta_1|| = {out['R0c']:.3e}  "
         f"(bound {R0C_BAND:.0e}, predicted 0.0)")
    # Separation diagnostics (§9 item 239), printed, never asserted.
    s1, s2 = it1["f"].sigma_field.x.array, it2["f"].sigma_field.x.array
    loc = float(np.max(np.abs(s2 - s1))) if s1.size else 0.0
    out["R0c_dsigma"] = comm.allreduce(loc, op=MPI.MAX)
    e1f, e2f = it1["f"].e_complex, it2["f"].e_complex
    n_e = e1f.function_space.dofmap.index_map.size_local * e1f.function_space.dofmap.index_map_bs
    e1, e2 = e1f.x.array[:n_e], e2f.x.array[:n_e]
    de2 = comm.allreduce(float(np.sum(np.abs(e2 - e1) ** 2)), op=MPI.SUM)
    ee2 = comm.allreduce(float(np.sum(np.abs(e1) ** 2)), op=MPI.SUM)
    out["R0c_dE"] = float(np.sqrt(de2 / ee2))
    _say(f"MAT5s3 (R0c diag) max|sigma_it2 - sigma_it1| over all local parent DG0 dofs = "
         f"{out['R0c_dsigma']:.3e} (predicted 0.0 exactly);  ||E_it2 - E_it1||/||E_it1|| owned = "
         f"{out['R0c_dE']:.3e} (predicted 1e-12 ... 1e-9; printed)")
    fx.records(it1["T"], "alpha=0 steady state (1 W)")
    rss = _rss_summed_gib(comm)
    _say(f"MAT5s3 W1 clocks: mesh+rung {fx.t_mesh:.1f} s;  EM s/solve "
         + "  ".join(f"{s:.1f}" for s in fx.em_secs)
         + f";  peak RSS summed over ranks {rss:.2f} GiB")
    out["t_mesh"] = fx.t_mesh
    out["em_max"] = max(fx.em_secs)
    return out


@pytest.fixture(scope="module")
def coupled():
    fx = _Fixture(alpha=ALPHA)
    comm = fx.comm
    out = {}
    t0 = time.perf_counter()

    # ---- Picard, θ₀ = 0, N = 8, G fixed from iterate 1.
    T_hist, its = [None], []
    for n in range(1, N_PICARD + 1):
        it = fx.phi(T_hist[-1])
        its.append(it)
        T_hist.append(it["T"])
        ts = _real(_allsum(comm, fx.theta(it["T"]) * fx.ds)) / fx.area
        sbar = _real(_allsum(comm, it["sig"] * fx.dx)) / fx.v_sub
        _say(f"MAT5s3 iterate {n}: sigma_bar_in = {sbar:.9f}  P_diss = {it['P']:.12e} W  "
             f"theta_s = {ts:.9f} K  EM {fx.em_secs[-1]:.1f} s")
    P1 = its[0]["P"]
    _say(f"MAT5s3 G = 1 W / P_1 = {fx.G:.12e} 1/W (fixed from iterate 1)")

    thetas = [fx.theta(T) for T in T_hist[1:]]  # θ₁ … θ₈
    deltas = [thetas[0]] + [fx.p1(thetas[k].x.array - thetas[k - 1].x.array)
                            for k in range(1, N_PICARD)]  # δ₁ … δ₈
    rho = [fx.ip(deltas[k + 1], deltas[k]) / fx.ip(deltas[k], deltas[k])
           for k in range(N_PICARD - 1)]  # ρ₁ … ρ₇
    ts_hist = [_real(_allsum(comm, th * fx.ds)) / fx.area for th in thetas]
    d_s = [ts_hist[0]] + [ts_hist[k] - ts_hist[k - 1] for k in range(1, N_PICARD)]
    s_ratios = [d_s[k + 1] / d_s[k] for k in range(N_PICARD - 1)]
    out["rho"] = rho
    _say("MAT5s3 rho_k (L2(sub)) k=1..7 = " + "  ".join(f"{r:+.6e}" for r in rho))
    _say("MAT5s3 surface-mean increment ratios k=1..7 = " + "  ".join(f"{r:+.6e}" for r in s_ratios))
    _say("MAT5s3 ||delta_k|| k=1..8 = " + "  ".join(f"{fx.norm(d):.4e}" for d in deltas))
    out["F"] = fx.norm(deltas[-1]) / fx.norm(thetas[-1])
    _say(f"MAT5s3 (F) ||delta_8||/||theta_8|| = {out['F']:.3e}  (bound {S_BAND:.0e}, "
         f"predicted < 1e-6)")

    # ---- Energy identity at iterate 8.
    it8 = its[-1]
    dep = fx.G * it8["P"]
    conv_int = _real(_allsum(comm, H_GEL * thetas[-1] * fx.ds))
    out["Estar"] = abs(conv_int - dep) / abs(dep)
    int_q = _real(_allsum(comm, it8["q"] * fx.dx))
    out["Bp"] = abs(int_q - dep) / abs(dep)
    r8 = it8["res"]
    out["Cp"] = abs(r8.convected_power - r8.source_power) / abs(r8.source_power)
    xp = ufl.SpatialCoordinate(fx.msh)
    xs = ufl.SpatialCoordinate(fx.sub)
    dx3 = ufl.Measure("dx", domain=fx.msh, subdomain_data=fx.cell_tags)(PHANTOM_CELL_TAG)
    m_par = _real(_allsum(comm, it8["f"].sigma_field * ufl.inner(xp, xp) * dx3))
    m_sub = _real(_allsum(comm, it8["sig"] * ufl.inner(xs, xs) * fx.dx))
    out["Dp"] = abs(m_par - m_sub) / abs(m_sub)
    _say(f"MAT5s3 (E*) oint h theta_8 = {conv_int:.15e}  G*P_diss,8 = {dep:.15e}  rel "
         f"{out['Estar']:.3e} (bound {C_BAND:.0e})")
    _say(f"MAT5s3 (B') int_sub q_8 = {int_q:.15e}  rel {out['Bp']:.3e} (bound {B_BAND:.0e})")
    _say(f"MAT5s3 (C') conv = {r8.convected_power:.15e}  src = {r8.source_power:.15e}  rel "
         f"{out['Cp']:.3e} (bound {C_BAND:.0e})")
    _say(f"MAT5s3 (D') int sigma x.x: parent = {m_par:.15e}  sub = {m_sub:.15e}  rel "
         f"{out['Dp']:.3e} (bound {D_BAND:.0e}; ghost propagation not gated)")

    # ---- μ_pred: FD power iteration of the re-solved-E Jacobian at θ_b = θ₇.
    T7, T8 = T_hist[N_PICARD - 1], T_hist[N_PICARD]
    eps = FD_RMS_K * np.sqrt(fx.v_sub)  # ‖v‖ = 1 ⇒ RMS(εv) = 1e-2 K
    v = fx.p1(np.ones_like(np.real(T7.x.array)))
    v = fx.p1(v.x.array / fx.norm(v))
    mu = []
    for j in range(M_POWER):
        Tp = fx.p1(T7.x.array + eps * v.x.array)
        itp = fx.phi(Tp)
        Jv = fx.p1((itp["T"].x.array - T8.x.array) / eps)
        mu.append(fx.ip(Jv, v) / fx.ip(v, v))
        _say(f"MAT5s3 power iteration j={j + 1}: mu = {mu[-1]:+.9e}  ||Jv|| = {fx.norm(Jv):.6e}  "
             f"EM {fx.em_secs[-1]:.1f} s")
        v = fx.p1(Jv.x.array / fx.norm(Jv))
    out["mu"] = mu
    mu_pred = mu[-1]

    # ---- μ_E (N2): frozen-|E| thermal-only power iteration, E₈.
    f8 = it8["f"]
    Wp, e2_par = _cell_average_dg0_tag(fx.msh, fx.cell_tags, PHANTOM_CELL_TAG,
                                       fx.G * 0.5 * ufl.inner(f8.e_complex, f8.e_complex))
    e2_sub = _to_sub(Wp, e2_par, fx.sub, fx.parent_ids)
    v = fx.p1(np.ones_like(np.real(T7.x.array)))
    v = fx.p1(v.x.array / fx.norm(v))
    mu_E = []
    for j in range(M_POWER):
        v_dg0 = fem.Function(fx.Ws)
        v_dg0.interpolate(v)
        v_dg0.x.scatter_forward()
        q = fem.Function(fx.Ws)
        q.x.array[:] = SALINE_SIGMA * ALPHA * np.real(v_dg0.x.array) * e2_sub.x.array
        q.x.scatter_forward()
        Jv = fx.theta(fx.steady(q).temperature)
        mu_E.append(fx.ip(Jv, v) / fx.ip(v, v))
        v = fx.p1(Jv.x.array / fx.norm(Jv))
    out["mu_E"] = mu_E[-1]

    _say("MAT5s3 mu_j j=1..5 = " + "  ".join(f"{m:+.9e}" for m in mu))
    _say("MAT5s3 mu_E,j j=1..5 = " + "  ".join(f"{m:+.9e}" for m in mu_E))
    out["K"] = abs(rho[3] / mu_pred - 1.0)
    out["Kconv"] = abs(mu[-1] / mu[-2] - 1.0)
    _say(f"MAT5s3 (K) rho_4 = {rho[3]:+.9e}  mu_pred = {mu_pred:+.9e}  |rho_4/mu_pred - 1| = "
         f"{out['K']:.4e}  (bound {K_BAND}; predicted mu_pred in [+0.03, +0.14])")
    _say(f"MAT5s3 (K-conv) |mu_5/mu_4 - 1| = {out['Kconv']:.3e}  (bound {KCONV_BAND:.0e}, "
         f"predicted <= 1e-3)")
    _say(f"MAT5s3 (M0) |mu_pred| = {abs(mu_pred):.4e}  (floor {M0_FLOOR:.0e})")
    _say(f"MAT5s3 (N1) frozen-sigma0 miss of (K) = |0 - mu_pred|/|mu_pred| = 1.000 (by (R0c), "
         f"its residual <= {R0C_BAND:.0e}, << M0 floor {M0_FLOOR:.0e})")
    _say(f"MAT5s3 (N2) mu_E = {mu_E[-1]:+.9e}  |mu_E/mu_pred - 1| = "
         f"{abs(mu_E[-1] / mu_pred - 1.0):.4e}  (predicted 0-50 %, printed)")

    # ---- Records.
    tb1, _, _ = fx.records(T_hist[1], "iterate 1 (sigma0)")
    tb8, _, _ = fx.records(T8, "iterate 8 (coupled)")
    sig8 = fx.sigma_sub(T8)
    sbar8 = _real(_allsum(comm, sig8 * fx.dx)) / fx.v_sub
    sbar8_in = _real(_allsum(comm, it8["sig"] * fx.dx)) / fx.v_sub
    _say(f"MAT5s3 record theta_bar_8/theta_bar_1 - 1 = {tb8 / tb1 - 1.0:+.6e} "
         f"(predicted mu/(1-mu) ≈ 3-16 %, mu_pred/(1-mu_pred) = {mu_pred / (1 - mu_pred):+.4e})")
    _say(f"MAT5s3 record sigma_bar(theta_8) = {sbar8:.9f} S/m  (iterate-8 input "
         f"{sbar8_in:.9f});  P_diss,8/P_1 = {it8['P'] / P1:.9f}")
    rss = _rss_summed_gib(comm)
    _say(f"MAT5s3 W2 clocks: mesh+rung {fx.t_mesh:.1f} s;  EM s/solve mean "
         f"{np.mean(fx.em_secs):.1f} (max {max(fx.em_secs):.1f}, n = {len(fx.em_secs)});  "
         f"loop {time.perf_counter() - t0:.1f} s;  peak RSS summed over ranks {rss:.2f} GiB")
    return out


# ---- W1 -------------------------------------------------------------------------
@complex_only
def test_r0a_hook_path_equals_p1_path(alpha0):
    assert alpha0["R0a_P"] <= L0_BAND, f"(R0a) P {alpha0['R0a_P']:.3e}"
    assert alpha0["R0a_q"] <= L0_BAND, f"(R0a) q {alpha0['R0a_q']:.3e}"


@complex_only
def test_r0b_alpha0_reproduces_p1_records(alpha0):
    assert alpha0["R0b_P"] <= REPRO_BAND, f"(R0b) P_1 {alpha0['R0b_P']:.3e}"
    assert alpha0["R0b_stor"] <= REPRO_BAND, f"(R0b) stored {alpha0['R0b_stor']:.3e}"
    assert alpha0["R0b_surf"] <= REPRO_BAND, f"(R0b) surface {alpha0['R0b_surf']:.3e}"
    assert alpha0["R0b_theta"] <= MEAN_THETA_HALF_ULP_K, f"(R0b) theta {alpha0['R0b_theta']:.3e}"
    assert alpha0["R0b_V"] <= V_BAND, f"(R0b) V_sub {alpha0['R0b_V']:.3e}"


@complex_only
def test_r0c_alpha0_map_residual_vanishes(alpha0):
    assert alpha0["R0c"] <= R0C_BAND, f"(R0c) {alpha0['R0c']:.3e}"


@complex_only
def test_d0_hook_touches_only_phantom_cells(alpha0):
    assert alpha0["D0"] == 0.0, f"(D0) {alpha0['D0']:.3e}"


# ---- W2 -------------------------------------------------------------------------
@complex_only
def test_k_contraction_matches_linearised_prediction(coupled):
    assert coupled["K"] <= K_BAND, f"(K) |rho_4/mu_pred - 1| = {coupled['K']:.4e}"


@complex_only
def test_k_power_iteration_converged(coupled):
    assert coupled["Kconv"] <= KCONV_BAND, f"(K-conv) {coupled['Kconv']:.3e}"


@complex_only
def test_m0_coupling_measurable(coupled):
    assert abs(coupled["mu"][-1]) >= M0_FLOOR, f"(M0) |mu_pred| = {abs(coupled['mu'][-1]):.3e}"


@complex_only
def test_f_fixed_point_converged(coupled):
    assert coupled["F"] <= S_BAND, f"(F) {coupled['F']:.3e}"


@complex_only
def test_energy_identity_at_fixed_point(coupled):
    assert coupled["Estar"] <= C_BAND, f"(E*) {coupled['Estar']:.3e}"
    assert coupled["Bp"] <= B_BAND, f"(B') {coupled['Bp']:.3e}"
    assert coupled["Cp"] <= C_BAND, f"(C') {coupled['Cp']:.3e}"


@complex_only
def test_dprime_sigma_transfer(coupled):
    assert coupled["Dp"] <= D_BAND, f"(D') {coupled['Dp']:.3e}"
