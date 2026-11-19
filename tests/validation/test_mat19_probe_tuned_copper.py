"""`MAT-19` step 3 -- the probe on the tuned copper F-small birdcage at 64 MHz (§9 item 265).

Brief: ``docs/planning/chunks/MAT-19.md`` (§9 item 265, chain P10).  Item 258's
probe series (``thermal.probe``, unchanged; ``_probe_set`` / ``A_PROBE`` /
``OMEGA_RADIUS`` / ``_lattice_q_mean`` imported) on `MAT-15` step 3's source
(item 236): the ×1 exact-strip hole mesh (80 798 cells), Leontovich copper,
Z_C(C_mode,Cu) on P2–P4 and C₁ across P1, one P1 drive at 64.0 MHz through
``_route_a_solve(..., degree=2)``, G = 1 W / P_acc, P1 thermal on the tag-3
submesh.  κ off, exact rectangles.

Item 236's ``_run()`` returns only scalars (no submesh / q_h), so its steps are
rebuilt here in the same order from the imported helpers (disclosed mapping):
``_read_t21`` → ``_fixture_tie`` → ``_exact_geo`` → (Φ) ``_sheet_azimuth_deg``
→ ``_route_a_solve(degree=2)`` → (W) → TH-14 (b) power block on the Z_C-replaced
specs → ``_tag_submesh`` → ``_robin_tags`` → ``_cell_average_dg0_tag`` →
``_to_sub``.  Item 236's module is not edited.

Asserted: (A0) tie, (Φ), 537 tag-3 cells, snapshot 0 = 0.0, (M); (Src) (W) and
(P); (A) P_Ω = ∫θ/V_sub ≤ ``B_BAND``; (U4) bitwise at width > 1; (W) width
invariance ≤ ``WIDTH_BAND`` at -n 4 / -n 2 against the -n 8 npz; (Σ) C / C′
≤ ``C4_COVARIANCE_BAND``.  Negative controls (rule (e), all *predicted*,
miss-only asserted): (U3c) iff G_w > 0, (C1″), (f″).  Printed: ball/point,
64 / 10 MHz ratios (i) and (ii), series table, ⟨θ⟩ vs item 236, (S), clocks.

Opt-in knob ``FEM_EM_MAT19_STEP3=1``; unset ⇒ skip; any other value raises.
"""

from __future__ import annotations

import dataclasses
import hashlib
import os
import resource
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.core.cavity import surface_resistance_ohm
from fem_em_solver.thermal import probe as _probe_mod
from fem_em_solver.thermal import probe_series
from fem_em_solver.thermal.probe import _ball_average_on_cells, _lattice
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import _sheet_azimuth_deg
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT
from tests.validation import test_mat15_tuned_copper_heating as _item236
from tests.validation import test_th14_birdcage_copper as th14
from tests.validation.test_ans7_exact_strip import (
    GATED_WIDTH_FRACTION,
    _assert_tie,
    _exact_geo,
    _fixture_tie,
    _route_a_solve,
)
from tests.validation.test_ans7_sparse_route_b import _ANS7, NODE_S_RTOL
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_mat10_sar_to_heat_sphere import _robin_tags, _to_sub
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    B_BAND,
    C4_COVARIANCE_BAND,
    DT_1,
    E_BAND,
    F_CONTROL_FLOOR,
    N_MAIN,
    RHO_C,
    T_INF,
    T_STAR,
    A,
    _allsum,
    _be_one,
    _cell_average_dg0_tag,
    _say,
    _tag_submesh,
    _transient,
)
from tests.validation.test_mat15_tuned_copper_heating import (
    DRIVE_HZ,
    PHI_BAND,
    _read_t21,
)
from tests.validation.test_mat19_probe_coil_phantom import (  # noqa: F401  (Ω's radius, used through _probe_set)
    A_PROBE,
    OMEGA_RADIUS,
    _lattice_q_mean,
    _probe_set,
)
from tests.validation.test_mat19_virtual_temperature_probe import WIDTH_BAND

KNOB = "FEM_EM_MAT19_STEP3"
N_SUB = 537  # 20261005T040135Z_MAT-15-step3-W1.log:1926 (item 236's printed record)
DT = T_STAR / N_MAIN
NPZ = "/workspace/logs/mat19/step3-n8.npz"
NPZ_WIDTH = 8
# Item 236's printed record, 20261005T040135Z_MAT-15-step3-W1.log:1930 (printed only).
MEAN_THETA_236_K = 0.687806
# Step 2's 10 MHz rises at t* per DEPOSITED watt (G = 1 W / P_1), printed records only,
# 20261006T005727Z_MAT-19-step2-a2-W1.log:1970 (the `MAT19-2   n=40` row).  Never comparands.
STEP2_N40_PER_DEP_W = {
    "A": 9.123011e-01,
    "B": 3.520711e00,
    "C": 6.976862e-01,
    "C'": 6.980991e-01,
    "D": 6.944133e-01,
    "E": 9.613614e-01,
    "Omega": 1.153314e00,
    "B+": 3.536194e00,
}


def _enabled():
    raw = os.environ.get(KNOB, "").strip()
    if not raw:
        return False
    if raw != "1":
        raise ValueError(f"{KNOB} must be unset or '1', got {raw!r}")
    return True


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _point_values(theta, centres, keys):
    """P1 θ at each centre: mean over owned containing cells (λ ≥ −1e-12); one allreduce."""
    V = theta.function_space
    comm = V.mesh.comm
    n_own = V.mesh.topology.index_map(3).size_local
    xdof = np.real(V.tabulate_dof_coordinates()).astype(np.float64)
    dm = V.dofmap.list[:n_own]
    th = np.real(theta.x.array)
    loc = np.zeros(2 * len(keys))
    if n_own:
        X = xdof[dm]  # (cells, 4, 3)
        Tm = np.transpose(X[:, 1:, :] - X[:, :1, :], (0, 2, 1))  # (cells, 3, 3)
        Tinv = np.linalg.inv(Tm)
        for i, k in enumerate(keys):
            p = np.asarray(centres[k], dtype=np.float64)
            l123 = np.einsum("cij,cj->ci", Tinv, p - X[:, 0, :])
            lam = np.concatenate([1.0 - l123.sum(axis=1)[:, None], l123], axis=1)
            hit = np.all(lam >= -1.0e-12, axis=1)
            if np.any(hit):
                loc[2 * i] = float(np.sum(lam[hit] * th[dm[hit]]))
                loc[2 * i + 1] = float(np.count_nonzero(hit))
    tot = comm.allreduce(loc, op=MPI.SUM)
    return {
        k: (
            tot[2 * i] / tot[2 * i + 1] if tot[2 * i + 1] > 0 else float("nan"),
            int(tot[2 * i + 1]),
        )
        for i, k in enumerate(keys)
    }


@pytest.fixture(scope="module")
def run():
    if not _enabled():
        pytest.skip(f"{KNOB} unset (heavy window, gitignored T21 inputs)")
    comm = MPI.COMM_WORLD
    clock = {}
    t0 = time.perf_counter()
    c_mode, s11_ref, _inm = _read_t21()
    tie = _fixture_tie()
    geo, specs, cells = _exact_geo(th14, comm, GATED_WIDTH_FRACTION)
    msh, cell_tags, facet_tags, _outer, tags_f, _pdefs = geo
    phi1 = _sheet_azimuth_deg(msh, tags_f, SHEET_IFACE + 1, comm)
    sin_phi = abs(np.sin(np.radians(phi1)))
    clock["mesh+tie"] = time.perf_counter() - t0
    _say(
        f"\nMAT19-3 C_mode,Cu = {c_mode:.15e} F (tune_p2.json); (F) cells {cells}; (Phi) phi_1 = "
        f"{phi1:.12f} deg |sin| = {sin_phi:.3e} (bound {PHI_BAND:.0e})"
    )

    t1 = time.perf_counter()
    s11, fields = _route_a_solve(
        th14, geo, specs, DRIVE_HZ, "cu", c_f=c_mode, return_fields=True, degree=2
    )
    comm.Barrier()
    clock["EM drive"] = time.perf_counter() - t1
    w_src = abs(s11 - s11_ref) / abs(s11_ref)
    _say(
        f"MAT19-3 (Src/W) S11 = {s11.real:+.15e}{s11.imag:+.15e}j rel to T21 {w_src:.3e} "
        f"(bound {NODE_S_RTOL:.0e}); EM drive {clock['EM drive']:.1f} s"
    )

    # (P) -- TH-14's (b) block on the Z_C-replaced specs (item 236's, same order).
    t2 = time.perf_counter()
    omega = 2.0 * np.pi * DRIVE_HZ
    z_c = _ANS7.z_cap(omega, c_mode)
    sp_list = [
        (
            dataclasses.replace(s, port_impedance_ohm=z_c)
            if i in _ANS7.TUNING_TERMINATED_INDICES
            else s
        )
        for i, s in enumerate(specs)
    ]
    e = fields.e_complex
    sigma_f = fields.sigma_field
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)

    def _vol(measure):
        form = fem.form(0.5 * sigma_f * ufl.inner(e, e) * measure)
        return float(np.real(comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)))

    acct = {"mesh": msh, "facet_tags": tags_f}
    sheet_objs = {sp.port_id: sp.sheet(driven=(sp.port_id == th14.DRIVEN)) for sp in sp_list}
    p_sheets = sum(
        float(_sheet_field_dissipation_w(acct, sh, e, omega)) for sh in sheet_objs.values()
    )
    p_src = float(_source_power_w(acct, sheet_objs[th14.DRIVEN], e, omega))
    p_total = _vol(dx)
    p_phantom = _vol(dx(th14.PHANTOM_CELL_TAG))
    z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, th14.COPPER_SIGMA)
    p_surf = th14._surface_loss_w(msh, facet_tags, z_s, e, comm)
    p_non = p_total - p_phantom
    p_rel = abs(p_src - p_sheets - p_phantom - p_non - p_surf) / abs(p_src)
    p_acc = p_src - p_sheets
    G = 1.0 / p_acc
    eta = p_phantom / p_acc
    _say(
        f"MAT19-3 (Src/P) identity rel {p_rel:.3e} (bound {DISCRETE_IDENTITY_RTOL:.0e}); P_acc = "
        f"{p_acc:.15e} W  G = {G:.12e} 1/W  eta = {eta:.9f}"
    )

    sub, parent_ids = _tag_submesh(msh, cell_tags, th14.PHANTOM_CELL_TAG)
    ft = _robin_tags(sub)
    n_sub = sub.topology.index_map(3).size_global
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    Wp, q_par = _cell_average_dg0_tag(
        msh, cell_tags, th14.PHANTOM_CELL_TAG, G * 0.5 * sigma_f * ufl.inner(e, e)
    )
    q_h = _to_sub(Wp, q_par, sub, parent_ids)
    dep = G * p_phantom
    int_q = float(np.real(_allsum(comm, q_h * ufl.dx(domain=sub))))
    _say(
        f"MAT19-3 sub cells {n_sub} (record {N_SUB}); V_sub = {v_sub:.12e} m^3; (B) int q_h vs "
        f"G P_phantom rel {abs(int_q - dep) / abs(dep):.3e} (printed)"
    )
    clock["power+source transfer"] = time.perf_counter() - t2

    t3 = time.perf_counter()
    res = _transient(
        sub,
        ft,
        q_h,
        dt=DT,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        snapshot_steps=range(N_MAIN + 1),
    )
    res_be = _be_one(sub, ft, q_h)
    clock["transient + BE"] = time.perf_counter() - t3

    V = res.snapshots[N_MAIN].function_space
    th = {}
    for n in range(N_MAIN + 1):
        f = fem.Function(V)
        f.x.array[:] = res.snapshots[n].x.array - T_INF
        f.x.scatter_forward()
        th[n] = f
    th_be = fem.Function(V)
    th_be.x.array[:] = res_be.temperature.x.array - T_INF
    th_be.x.scatter_forward()

    names, order, probes, centres = _probe_set()
    t4 = time.perf_counter()
    lam = _lattice(3)
    s = probe_series(th, DT, probes)
    s_be = probe_series({1: th_be}, DT_1, probes)
    clock["probe series"] = time.perf_counter() - t4

    t5 = time.perf_counter()
    qbar = {}
    for k, p in zip(order, probes):
        qbar[k], _v = _lattice_q_mean(V, q_h, p, lam)
    point = _point_values(th[N_MAIN], centres, names)
    clock["(S) + point"] = time.perf_counter() - t5

    means = np.array(
        [
            float(np.real(_allsum(comm, th[n] * ufl.dx(domain=sub)))) / v_sub
            for n in range(N_MAIN + 1)
        ]
    )
    mean_be = float(np.real(_allsum(comm, th_be * ufl.dx(domain=sub)))) / v_sub

    imap = sub.topology.index_map(3)
    part = comm.allgather((int(imap.size_local), int(imap.num_ghosts)))
    o_w = sum(1 for so, _g in part if so > 0)
    g_w = sum(g for _so, g in part)
    _say(
        f"MAT19-3 (part) -n {comm.size}: sub owned {[p[0] for p in part]} ghosts "
        f"{[p[1] for p in part]} O_w {o_w} G_w {g_w}; partition exercised: {o_w >= 2}"
    )

    ghost = None
    if comm.size > 1:
        allc = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
        gval, gvol = _ball_average_on_cells(th[N_MAIN], probes[order.index("Omega")], allc)
        ghost = (abs(gvol / v_sub - 1.0), abs(gval / means[N_MAIN] - 1.0))

    allrows = comm.allgather(s.values[N_MAIN].copy())
    u4 = all(np.array_equal(allrows[0], r) for r in allrows)

    rss = comm.allreduce(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, op=MPI.SUM)
    hashes = (_sha(Path(__file__).resolve()), _sha(_probe_mod.__file__), _sha(_item236.__file__))
    clock["fixture total"] = time.perf_counter() - t0
    return dict(
        cells=cells,
        tie=tie,
        phi1=phi1,
        sin_phi=sin_phi,
        w_src=w_src,
        p_rel=p_rel,
        eta=eta,
        G=G,
        n_sub=n_sub,
        v_sub=v_sub,
        s=s,
        s_be=s_be,
        qbar=qbar,
        point=point,
        means=means,
        mean_be=mean_be,
        ghost=ghost,
        u4=u4,
        names=names,
        order=order,
        centres=centres,
        hashes=hashes,
        clock=clock,
        comm=comm,
        o_w=o_w,
        g_w=g_w,
        rss_kib=rss,
    )


@complex_only
def test_mat19_step3_probe_tuned_copper(run):  # noqa: C901
    comm = run["comm"]
    s, s_be = run["s"], run["s_be"]
    order, names, centres = run["order"], run["names"], run["centres"]
    v = s.values
    j = {k: order.index(k) for k in order}
    om = "Omega"

    _say(f"\nMAT19-3 width -n {comm.size}; dt = {DT:.4f} s; t* = {T_STAR:.3f} s")
    _say(
        "MAT19-3 centres [m]: "
        + "  ".join(
            f"{k} ({centres[k][0]:+.4f},{centres[k][1]:+.4f},{centres[k][2]:+.4f})" for k in order
        )
    )
    _say(
        "MAT19-3 ball V_h/(4 pi a^3/3) - 1: "
        + "  ".join(
            f"{k} {s.volumes[j[k]] / s.nominal_volumes[j[k]] - 1.0:+.3e}" for k in order if k != om
        )
        + " (predicted |.| <= 3e-2, printed)"
    )
    m_ok = all(
        np.hypot(centres[k][0], centres[k][1]) + A_PROBE <= A
        and abs(centres[k][2]) + A_PROBE <= PHANTOM_HEIGHT / 2.0
        for k in order
        if k != om
    )
    _say(
        f"MAT19-3 (A0) cells {run['cells']}; sub cells {run['n_sub']} (record {N_SUB}); "
        f"(M) every probe ball inside the phantom: {m_ok}"
    )

    PO = v[:, j[om]]
    relA = np.abs(PO[1:] - run["means"][1:]) / np.abs(run["means"][1:])
    relA_be = abs(s_be.values[0, j[om]] - run["mean_be"]) / abs(run["mean_be"])
    volA = abs(s.volumes[j[om]] / run["v_sub"] - 1.0)
    _say(
        f"MAT19-3 (A) P_Omega vs int theta/V_sub: max n=1..40 {np.max(relA):.3e}, BE theta^1 "
        f"{relA_be:.3e}; V_Omega/V_sub - 1 = {volA:.3e} (bound {B_BAND:.0e}); n=0: P_Omega "
        f"{PO[0]!r} mean {run['means'][0]!r}"
    )
    r0 = abs(run["means"][N_MAIN] / MEAN_THETA_236_K - 1.0)
    _say(
        f"MAT19-3 <theta>(t*) per accepted watt = {run['means'][N_MAIN]:.6f} K vs item 236's "
        f"{MEAN_THETA_236_K} rel {r0:.3e} (expected <= {E_BAND:.0e} at 6-digit record, printed)"
    )

    _say("MAT19-3 series theta_bar [K per accepted W] (n: " + " ".join(order) + "):")
    for n in (1, 10, 20, 30, 40):
        _say(f"MAT19-3   n={n:2d} t={n * DT:7.2f}  " + " ".join(f"{v[n, j[k]]:.6e}" for k in order))
    _say(f"MAT19-3 max_imag_ratio {s.max_imag_ratio:.3e}; BE {s_be.max_imag_ratio:.3e}")

    # (P4′) ball/point − 1 at t*
    _say(
        "MAT19-3 (P4') ball(t*)/point(t*) - 1: "
        + "  ".join(
            f"{k} {v[N_MAIN, j[k]] / run['point'][k][0] - 1.0:+.4e} (point {run['point'][k][0]:.6e}, "
            f"{run['point'][k][1]} cells)"
            for k in names
        )
        + " (predicted |.| <= 2 %, printed)"
    )
    # 64 / 10 MHz ratios
    keys_r = names + ["C'", om]
    _say(
        "MAT19-3 ratio (i) [MIXED normalisation: 64 MHz per ACCEPTED W / 10 MHz per DEPOSITED W, "
        "step 2 log 20261006T005727Z_MAT-19-step2-a2-W1.log:1970]: "
        + "  ".join(f"{k} {v[N_MAIN, j[k]] / STEP2_N40_PER_DEP_W[k]:.6f}" for k in keys_r)
    )
    _say(
        f"MAT19-3 ratio (ii) [per DEPOSITED W at both, 64 MHz rise / eta = {run['eta']:.9f}; "
        "shape-only]: "
        + "  ".join(
            f"{k} {v[N_MAIN, j[k]] / run['eta'] / STEP2_N40_PER_DEP_W[k]:.6f}" for k in keys_r
        )
    )

    c1 = abs(v[N_MAIN, j["B+"]] - v[N_MAIN, j["B"]]) / v[N_MAIN, j["B"]]
    _say(
        f"MAT19-3 (C1'') |B+ - B|/B at t* = {c1:.4e} (predicted [1e-4, 1e-2]; asserted > "
        f"{WIDTH_BAND:.0e})"
    )
    sig = abs(v[N_MAIN, j["C"]] - v[N_MAIN, j["C'"]]) / v[N_MAIN, j["C"]]
    qc, qcp = run["qbar"]["C"], run["qbar"]["C'"]
    vc, vcp = s.volumes[j["C"]], s.volumes[j["C'"]]
    _say(
        f"MAT19-3 (Sigma) |C - C'|/C at t* = {sig:.4e} (bound {C4_COVARIANCE_BAND:.0e}); q_h ball "
        f"pair C {qc:.6e} C' {qcp:.6e} rel {abs(qc - qcp) / abs(qc):.4e}; ball volumes C "
        f"{vc:.6e} C' {vcp:.6e} (printed)"
    )
    ff = abs(v[N_MAIN, j["C"]] - v[N_MAIN, j["B"]]) / v[N_MAIN, j["B"]]
    _say(
        f"MAT19-3 (f'') |C - B|/B at t* = {ff:.4f} (predicted [0.1, 1.0], asserted > "
        f"{F_CONTROL_FLOOR:.0e})"
    )
    _say(
        "MAT19-3 (S) rho c theta_bar(theta^1)/(dt_1 q_bar) - 1: "
        + "  ".join(
            f"{k} {RHO_C * s_be.values[0, j[k]] / (DT_1 * run['qbar'][k]) - 1.0:+.4e}"
            for k in order
        )
        + " (printed)"
    )
    if run["ghost"] is not None:
        _say(
            f"MAT19-3 (U3c) ghost-included Omega: volume miss {run['ghost'][0]:.4e}, value miss "
            f"{run['ghost'][1]:.4e} (asserted > {B_BAND:.0e} iff G_w > 0; predicted G_w/{N_SUB} = "
            f"{run['g_w'] / N_SUB:.4e} within x[0.5, 2], printed)"
        )
        if run["g_w"] == 0:
            _say(f"MAT19-3 (U3c) inapplicable at -n {comm.size}: G_w = 0 (a record, not a red)")
    if run["o_w"] < 2:
        _say(
            f"MAT19-3 (A)/(U4) on a submesh partition: not exercised at -n {comm.size} "
            f"(O_w = {run['o_w']}; a record, not a red)"
        )
    _say(f"MAT19-3 (U4) bitwise rank agreement at n=40: {run['u4']} (-n {comm.size})")

    w = None
    if comm.size == NPZ_WIDTH:
        if comm.rank == 0:
            os.makedirs(os.path.dirname(NPZ), exist_ok=True)
            np.savez(
                NPZ,
                series=v,
                volumes=s.volumes,
                cells=np.array([run["n_sub"]]),
                sha_module=run["hashes"][0],
                sha_probe=run["hashes"][1],
                sha_item236=run["hashes"][2],
            )
            print(f"MAT19-3 (W) wrote {NPZ} (sha module {run['hashes'][0][:12]})", flush=True)
        comm.Barrier()
    else:
        if not os.path.exists(NPZ):
            pytest.fail(f"(W) {NPZ} missing: run W1 (-n 8) first; never recompute")
        z = np.load(NPZ)
        if (
            str(z["sha_module"]) != run["hashes"][0]
            or str(z["sha_probe"]) != run["hashes"][1]
            or str(z["sha_item236"]) != run["hashes"][2]
        ):
            pytest.fail("(W) npz hash differs from this module / probe.py / item 236's module")
        cells_eq = int(z["cells"][0]) == run["n_sub"]
        rel = np.abs(v[1:] / z["series"][1:] - 1.0)
        w = float(np.max(rel))
        _say(
            f"MAT19-3 (W) -n {comm.size} vs -n {NPZ_WIDTH}: cells equal {cells_eq}; series max rel "
            f"{w:.3e} (band {WIDTH_BAND:.0e}); per probe max: "
            + "  ".join(f"{k} {np.max(rel[:, j[k]]):.2e}" for k in order)
        )
    _say(
        "MAT19-3 phase clocks: "
        + "  ".join(f"{k} {t:.1f} s" for k, t in run["clock"].items())
        + f"; summed ru_maxrss {run['rss_kib'] / 2**20:.2f} GiB"
    )

    assert run["sin_phi"] <= PHI_BAND, f"(Phi) {run['sin_phi']:.3e}"
    _assert_tie(run["cells"], run["tie"])
    assert run["n_sub"] == N_SUB, f"(A0) sub cells {run['n_sub']} != {N_SUB}"
    assert m_ok, "(M) a probe ball leaves the phantom"
    assert np.all(v[0] == 0.0), "(A0) snapshot 0 does not read 0.0"
    assert run["w_src"] <= NODE_S_RTOL, f"(Src/W) {run['w_src']:.3e}"
    assert run["p_rel"] <= DISCRETE_IDENTITY_RTOL, f"(Src/P) {run['p_rel']:.3e}"
    assert PO[0] == 0.0 and run["means"][0] == 0.0
    assert np.max(relA) <= B_BAND, f"(A) {np.max(relA):.3e}"
    assert relA_be <= B_BAND, f"(A) BE {relA_be:.3e}"
    assert volA <= B_BAND, f"(A) volume {volA:.3e}"
    if comm.size > 1:
        assert run["u4"], "(U4) ranks disagree bitwise"
    if run["ghost"] is not None and run["g_w"] > 0:
        assert run["ghost"][0] > B_BAND, "(U3c) ghost-included volume did not miss"
    assert sig <= C4_COVARIANCE_BAND, f"(Sigma) {sig:.4e}"
    assert c1 > WIDTH_BAND, f"(C1'') {c1:.3e}"
    assert ff > F_CONTROL_FLOOR, f"(f'') {ff:.4f}"
    if w is not None:
        assert cells_eq, "(W) cells differ"
        assert w <= WIDTH_BAND, f"(W) {w:.3e} > {WIDTH_BAND:.0e}"
