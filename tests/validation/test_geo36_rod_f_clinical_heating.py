"""GEO-36 I32a (§9 item 707, §10 chain I32a): the clinical implant-heating body —
arm rod's and arm gel twin's (1, +1j) SAR heating the ASTM gel box around the rod
inside F-clinical-1.5T, composed and stub-run before ``xl`` entry 55.

Brief: ``docs/planning/chunks/GEO-36.md`` "§9 item 707".  Composition: I36b's
drive pieces by name (`test_geo36_rod_f_clinical_drive.py` stays byte-identical;
`_arm` frees its fields, so `_heat_columns` holds the two `_drive` outputs and
forms ``e_sup`` as `_arm` does), I39's heating arm re-centred on the clinical rod,
and I41's matched thermal domain.

The rod is kept cells at σ_wire (tag 403), not a PEC hole; it carries no source,
it is outside the thermal domain, and its wall is adiabatic.  The gel-twin arm
heats tags 3 + 403 as one gel (`_merge_rod_into_phantom`), so it has no wall.
One G = 1 W / P_gel(gel twin) on both arms.  I36 calls these arms M / R and I41
calls them R / M, so this module says rod / gel twin.

MAT-21's time base, not MAT-15's: t* = 0.2 T_C, Δt = t*/N, Δt₁ from
`test_mat21_box_fem.py`.  `_cn` / `_be_one` / `_e_prime` close over MAT-15's t*,
so `_transient` is called directly.

Selected by ``FEM_EM_GEO36_S8_HEAT`` only (never by ``-k``): unset ⇒ the gated
test skips and the coding tests run (no gmsh); ``1`` ⇒ `_run_heat`, ``xl`` entry
55 only, with the node ``FEM_EM_GEO36_S8_DRIVE_C_PF`` read through `_node_knob`;
anything else ⇒ ``ValueError``.  Degree 1 only: a set
``FEM_EM_GEO36_S8_DRIVE_DEGREE`` is refused.

Asserted on each arm (bands imported): (B), (E), (E′), (S2); on the rod arm
(S2-tip); on the mesh (S0), the merged arm's tag 3 = the rod arm's 3 + 403; I36b's
(M2), (Z), (R1), (P) and NEG on both arms.  Every gate is collected, then asserted.
Printed: θ(t*) at both tip balls and P37's six points on both arms, the tip θ
ratio rod / gel twin beside I41's 1.10, the rod's power share, per-phase time.
θ here is on the (1, +1j) sense, entry 48's; entry 52's lattice reads ccw, so a θ
figure from this module does not compare with C6's.

Scope: the body and its stub proof.  Not claimed: any clinical window, a price,
a θ figure, a tissue or published comparand.
"""

from __future__ import annotations

import dataclasses
import gc
import os
import sys
import time

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.post.sar import averaging_ball_radius, mean_sar
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import ROD_LENGTH
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_geo36_rod_f_clinical_drive import (
    S8_DRIVE_DEGREE_ENV,
    _a0m_failures,
    _drive,
    _drive_env,
    _gel_twin_map,
    _guard_foreign,
    _k9_fixture,
    _k9_stubs,
    _m2_failures,
    _node_knob,
    _partition,
    _probe_points,
    _route_drive,
    _vol_powers,
)
from tests.validation.test_imp3_rod_heating_matched import (
    _merge_rod_into_phantom,
    _owned_tag_count,
)
from tests.validation.test_imp3_rod_heating_phantom import WALL_SELECT_M
from tests.validation.test_mass_averaged_sar_standard_masses import (
    ONE_GRAM_KG,
    QUADRATURE_DEGREE,
)
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, ROBIN_TAG, _allsum, _to_sub
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    E_BAND,
    H_GEL,
    K_GEL,
    N_MAIN,
    RHO_C,
    S2_BAND,
    T_INF,
    _cell_average_dg0_tag,
    _max_imag,
    _tag_submesh,
    _transient,
)
from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX
from tests.validation.test_mat21_box_fem import DT_1 as DT_1_BOX
from tests.validation.test_mat21_box_fem import T_STAR as T_STAR_BOX
from tests.validation.test_mat21_box_probe_points_fem import PROBE_LABELS, PROBE_PTS
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

HEAT_ENV = "FEM_EM_GEO36_S8_HEAT"
SENSE = "(1,+1j)"
#: I41's matched tip θ ratio, the brief's literal; printed beside ours, never asserted.
I41_TIP_RATIO = 1.10
#: Named controls (brief): the perturbed source's size and the e₂ = 0 miss.
PERTURBED_SOURCE = 1.0e-3
PERTURBED_B_RATIO = PERTURBED_SOURCE / B_BAND
SINGLE_COLUMN_MISS = 0.5
SINGLE_COLUMN_CEILING = 1.0e-10  # `_k8_arm`'s, `test_geo36_rod_f_clinical_drive.py:1090`


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-36 I32] {msg}", flush=True)


def _heat_mode():
    raw = os.environ.get(HEAT_ENV, "").strip()
    if raw == "":
        return None
    if raw == "1":
        return raw
    raise ValueError(f"{HEAT_ENV} must be unset or '1', got {raw!r}")


def _refuse_degree():
    raw = os.environ.get(S8_DRIVE_DEGREE_ENV, "").strip()
    if raw != "":
        raise ValueError(f"{S8_DRIVE_DEGREE_ENV}={raw!r} is refused: I32 is degree 1 only")


def _box_preconditions():
    """``H_GEL`` / ``K_GEL`` / ρc are the box's: `_transient` carries MAT-15's."""
    from tests.validation import test_mat21_box_closed_form as box

    assert H_GEL == box.H, (H_GEL, box.H)
    assert K_GEL == box.K, (K_GEL, box.K)
    assert RHO_C == box.RHO_C, (RHO_C, box.RHO_C)


def _map_box_point(point, off):
    """A point of `MAT-21`'s frame (A, B, C on x, y, z) in the fixture's (Lz = 2A,
    Lx = 2B, Ly = 2C, centred at (0, off, 0)): (a, b, c) → (b, off + c, a)."""
    a, b, c = (float(v) for v in point)
    mapped = (b, off + c, a)
    inside = abs(mapped[0]) <= B_BOX and abs(mapped[1] - off) <= C_BOX and abs(mapped[2]) <= A_BOX
    assert inside, (point, mapped)
    return mapped


def _p37_fixture(off):
    """P37's six points in the fixture frame, keyed by label."""
    return {lab: _map_box_point(p, off) for lab, p in zip(PROBE_LABELS, PROBE_PTS)}


def _dist_to_rod_axis(x, centre):
    """Distance of points x (3, n) from the clinical rod's axis segment (on ẑ)."""
    half = 0.5 * ROD_LENGTH
    zc = np.clip(x[2] - centre[2], -half, half)
    return np.sqrt((x[0] - centre[0]) ** 2 + (x[1] - centre[1]) ** 2 + (x[2] - centre[2] - zc) ** 2)


def _outer_robin_tags(sub, centre):
    """I39's `_outer_robin_tags` on ``centre``: Robin away from the rod wall."""
    tdim = sub.topology.dim
    fdim = tdim - 1
    sub.topology.create_connectivity(fdim, tdim)
    ext = np.sort(np.asarray(dmesh.exterior_facet_indices(sub.topology), dtype=np.int32))
    mids = dmesh.compute_midpoints(sub, fdim, ext) if ext.size else np.zeros((0, 3))
    keep = _dist_to_rod_axis(mids.T, centre) > WALL_SELECT_M if ext.size else np.zeros(0, bool)
    outer = ext[keep].astype(np.int32)
    wall = ext[~keep]
    n_owned = sub.topology.index_map(fdim).size_local
    comm = sub.comm
    counts = {
        "exterior": comm.allreduce(int(np.count_nonzero(ext < n_owned)), op=MPI.SUM),
        "outer": comm.allreduce(int(np.count_nonzero(outer < n_owned)), op=MPI.SUM),
        "wall": comm.allreduce(int(np.count_nonzero(wall < n_owned)), op=MPI.SUM),
    }
    tags = dmesh.meshtags(sub, fdim, outer, np.full(outer.size, ROBIN_TAG, dtype=np.int32))
    return tags, counts


def _heat_columns(env, solver, specs, feeds, order, z_th, v_s):
    """`_arm`'s column loop (`test_geo36_rod_f_clinical_drive.py:432–447`), with the
    fields kept: S, kinds, times, the partition parts, ``e_sup`` and σ."""
    from tests.validation import test_wf10_f_human_copper_esr_root as h17

    cur, kinds, times, parts, es = [], [], [], [], []
    sigma_f = None
    for n_drive, pid in enumerate(feeds):
        d = _drive(env, solver, specs, pid, held=(n_drive > 0))
        kinds.append(d["kind"])
        times.append(d["t_solve"] + d["t_post"])
        cur.append(np.array([d["currents"][p] for p in order], dtype=np.complex128))
        parts.append({k: v for k, v in d.items() if k.startswith("p_")})
        es.append(d["e"])
        sigma_f = d["sigma"]
        del d
    s2 = h17._current_route_s(cur, [0, 1], z_th, v_s)
    w = h17.SENSES[SENSE]
    e_sup = es[0].copy()
    e_sup.x.array[:] = w[0] * es[0].x.array + w[1] * es[1].x.array
    e_sup.x.scatter_forward()
    return {
        "s2": s2,
        "kinds": kinds,
        "times": times,
        "parts": parts,
        "e_sup": e_sup,
        "sigma": sigma_f,
    }


def _heat_source(msh, cell_tags, e_sup, sigma_f, G, sub, parent_ids):
    """q = G · (DG0 cell average of ½σ|E|² on the gel tag), moved to the submesh."""
    Wp, q_par = _cell_average_dg0_tag(
        msh, cell_tags, PHANTOM_CELL_TAG, G * 0.5 * sigma_f * ufl.inner(e_sup, e_sup)
    )
    return _to_sub(Wp, q_par, sub, parent_ids)


def _e_prime_box(sub, ds, res, deposited, dt):
    """`_e_prime` on this module's Δt: ρc∫(θᴺ−θᴺ⁻¹)/Δt + ∮hθ^(n+½) against G·P."""
    comm = sub.comm
    tn, tm = res.snapshots[N_MAIN], res.snapshots[N_MAIN - 1]
    stor = float(np.real(_allsum(comm, RHO_C * (tn - tm) / dt * ufl.dx(domain=sub))))
    conv = float(np.real(_allsum(comm, H_GEL * (0.5 * (tn + tm) - T_INF) * ds)))
    return abs(stor + conv - deposited) / abs(deposited), stor, conv


def _ball_reader(sub, fields):
    """``(centre, radius) -> (volume, {name: ball mean})`` over ``fields`` on ``sub``."""
    comm = sub.comm
    centre_c = fem.Constant(sub, np.zeros(3, dtype=default_scalar_type))
    rad_c = fem.Constant(sub, default_scalar_type(0.0))
    off = ufl.SpatialCoordinate(sub) - centre_c
    inside = ufl.conditional(ufl.lt(ufl.real(ufl.dot(off, off)), ufl.real(rad_c) ** 2), 1.0, 0.0)
    dxq = ufl.dx(domain=sub, metadata={"quadrature_degree": int(QUADRATURE_DEGREE)})
    f_vol = fem.form(inside * dxq)
    forms = {name: fem.form(f * inside * dxq) for name, f in fields.items()}

    def _read(centre, radius):
        centre_c.value[:] = centre
        rad_c.value = radius
        vol = float(np.real(comm.allreduce(fem.assemble_scalar(f_vol), op=MPI.SUM)))
        return vol, {
            name: float(np.real(comm.allreduce(fem.assemble_scalar(f), op=MPI.SUM))) / vol
            for name, f in forms.items()
        }

    return _read


def _probe_theta(T, pts, parent):
    """θ = Re(T) − T_INF at ``pts`` on a submesh (item 700: ``parent=``), rank-safe."""
    comm = T.function_space.mesh.comm
    pts = np.ascontiguousarray(np.asarray(pts, dtype=np.float64))
    vals, ok = evaluate_vector_field_parallel(T, pts, owner="lowest_original_cell", parent=parent)
    vals = comm.bcast(vals if comm.rank == 0 else None, root=0)
    ok = np.asarray(comm.bcast(ok if comm.rank == 0 else None, root=0), dtype=bool)
    th = np.array(
        [float(np.real(np.ravel(vals[i])[0])) - T_INF if ok[i] else np.nan for i in range(len(pts))]
    )
    return th, ok


def _projected_ball_comparands(sub, q, V1, centres, radius, tag):
    """I32b's shared tip comparand (§9 item 870; U127′ imports it by name):
    ``{name: {"vol", "q", "pq"}}`` — the ball means of the DG0 source q and of its
    consistent-mass L2 projection Πq onto V1, through `OPS-83`'s `_projections`
    (imported, never copied).  ``tag`` must be unique per call and space-free:
    `_projections` builds its PETSc options prefix from it."""
    # Function-level import: `test_ops83_s2tip_form_on_stub` imports this module at its top.
    from tests.validation.test_ops83_s2tip_form_on_stub import _projections

    pq, _ = _projections(q, V1, tag)
    read = _ball_reader(sub, {"q": q, "pq": pq})
    out = {}
    for name, c in centres.items():
        vol, m = read(np.asarray(c, dtype=float), radius)
        out[name] = {"vol": vol, "q": m["q"], "pq": m["pq"]}
    return out


def _tip_slope_rows(theta1_means, dt, comparands, rho_c):
    """Per ball: ``R_Pi`` = |ρc⟨θ¹⟩/(Δt⟨Πq⟩) − 1|, ``R_DG0`` = |ρc⟨θ¹⟩/(Δt⟨q⟩) − 1|
    (the retired form), ``pq_q`` = ⟨Πq⟩/⟨q⟩ − 1 (signed).  Plain arithmetic, no mesh."""
    rows = {}
    for name, th in theta1_means.items():
        c = comparands[name]
        lhs = rho_c * th / dt
        rows[name] = {
            "R_Pi": abs(lhs / c["pq"] - 1.0),
            "R_DG0": abs(lhs / c["q"] - 1.0),
            "pq_q": c["pq"] / c["q"] - 1.0,
        }
    return rows


#: The F2182-shaped read (printed, never asserted): quoted from memory, unconfirmed (W3).
F2182_T_S = 900.0  # quoted from memory, unconfirmed (W3)
F2182_OFFSET_M = 1.0e-3  # quoted from memory, unconfirmed (W3)


def _thermal_arm(
    label,
    msh,
    cell_tags,
    e_sup,
    sigma_f,
    G,
    centre,
    tips,
    probes,
    *,
    merged,
    tip_projection=False,
    f2182=None,
):
    """I39's `_heat_arm` on the box's time base; ``G=None`` sets G = 1 / P_gel here.

    I32b keywords, both default off (P61 / C31 / K10es reach the unedited path):
    ``tip_projection=True`` adds one BE step at Δt₁/10 and the projected tip
    comparand (``tip_s2_pi`` = R_Π at Δt₁/10, the gate); ``f2182`` (name → point)
    adds one CN run to `F2182_T_S`, before the t* CN so that stays the last θ = ½ call."""
    comm = msh.comm
    t0 = time.perf_counter()
    P = float(
        np.real(
            mean_sar(
                e_sup,
                sigma=sigma_f,
                rho=PHANTOM_RHO_KG_PER_M3,
                cell_tags=cell_tags,
                comm=comm,
                subdomain_ids=PHANTOM_CELL_TAG,
            )["dissipated_power_w"]
        )
    )
    if G is None:
        G = 1.0 / P
    dep = G * P
    sub, parent_ids = _tag_submesh(msh, cell_tags, PHANTOM_CELL_TAG)
    ft, counts = _outer_robin_tags(sub, centre)
    one = fem.Constant(sub, default_scalar_type(1.0))
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    area = float(np.real(_allsum(comm, one * ds)))
    q = _heat_source(msh, cell_tags, e_sup, sigma_f, G, sub, parent_ids)
    int_q = float(np.real(_allsum(comm, q * ufl.dx(domain=sub))))
    err_B = abs(int_q - dep) / abs(dep)
    t_source = time.perf_counter() - t0
    _say(
        comm,
        f"[{label}] sense {SENSE}: P_gel {P:.12e} W; G {G:.12e} 1/W; G*P {dep:.12e} W; V_sub "
        f"{v_sub:.12e} m^3; Robin area {area:.9e} m^2; facets exterior {counts['exterior']} = "
        f"outer {counts['outer']} + rod wall {counts['wall']} (adiabatic)",
    )
    _say(
        comm,
        f"[{label}] (B) int_sub q {int_q:.15e} vs G*P {dep:.15e}: rel {err_B:.3e} "
        f"(band {B_BAND:.0e}, ASSERTED)",
    )

    f2182_th = None
    if f2182:
        res900 = _transient(
            sub, ft, q, dt=F2182_T_S / N_MAIN, n_steps=N_MAIN, theta=0.5, startup_half_steps=2
        )
        f_names = list(f2182)
        f_th, f_ok = _probe_theta(
            res900.temperature, [f2182[n] for n in f_names], (msh, parent_ids)
        )
        f2182_th = dict(zip(f_names, f_th))
        for n, v, valid in zip(f_names, f_th, f_ok):
            _say(
                comm,
                f"[{label}] F2182-shaped probe {n} {tuple(round(float(x), 6) for x in f2182[n])} "
                f"({F2182_OFFSET_M * 1e3:.1f} mm beyond the apex), sense {SENSE}: theta"
                f"({F2182_T_S:.0f} s) {v:+.9e} K (valid {bool(valid)}; placement and time "
                f"quoted from memory, unconfirmed; W3; printed)",
            )
        del res900

    dt = T_STAR_BOX / N_MAIN
    t1 = time.perf_counter()
    res = _transient(
        sub,
        ft,
        q,
        dt=dt,
        n_steps=N_MAIN,
        theta=0.5,
        startup_half_steps=2,
        snapshot_steps=(N_MAIN - 1, N_MAIN),
    )
    t_cn = time.perf_counter() - t1
    ep, stor, conv = _e_prime_box(sub, ds, res, dep, dt)
    res1 = _transient(sub, ft, q, dt=DT_1_BOX, n_steps=1, theta=1.0)
    thbar1 = float(np.real(_allsum(comm, (res1.temperature - T_INF) * ufl.dx(domain=sub)))) / v_sub
    err_S2 = abs(RHO_C * thbar1 / (DT_1_BOX * dep / v_sub) - 1.0)
    imag = max(_max_imag(res.temperature), _max_imag(res1.temperature))
    res10 = None
    if tip_projection:
        res10 = _transient(sub, ft, q, dt=DT_1_BOX / 10.0, n_steps=1, theta=1.0)
        imag = max(imag, _max_imag(res10.temperature))

    V1 = res.temperature.function_space
    th, th1 = fem.Function(V1), fem.Function(V1)
    th.x.array[:] = res.temperature.x.array - T_INF
    th1.x.array[:] = res1.temperature.x.array - T_INF
    th.x.scatter_forward()
    th1.x.scatter_forward()
    fields = {"th": th, "th1": th1, "q": q}
    if tip_projection:
        th10 = fem.Function(V1)
        th10.x.array[:] = res10.temperature.x.array - T_INF
        th10.x.scatter_forward()
        fields["th10"] = th10
    read = _ball_reader(sub, fields)
    a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    tip_th, tip_s2 = {}, {}
    tip_m = {}
    for name, c in tips.items():
        vol, m = read(np.asarray(c, dtype=float), a_1g)
        tip_m[name] = m
        tip_th[name] = m["th"]
        tip_s2[name] = abs(RHO_C * m["th1"] / (DT_1_BOX * m["q"]) - 1.0)
        if tip_projection:
            _say(
                comm,
                f"[{label}] 1 g ball (a {a_1g * 1e3:.4f} mm) at {name}, sense {SENSE}: theta(t*) "
                f"{m['th']:.9e} K, ball volume in the thermal domain {vol:.6e} m^3; retired form "
                f"(S2-tip, DG0) {tip_s2[name]:.4e} (band {S2_BAND:.0e}, printed)",
            )
        else:
            _say(
                comm,
                f"[{label}] 1 g ball (a {a_1g * 1e3:.4f} mm) at {name}, sense {SENSE}: theta(t*) "
                f"{m['th']:.9e} K, ball volume in the thermal domain {vol:.6e} m^3; (S2-tip) "
                f"{tip_s2[name]:.4e} (band {S2_BAND:.0e})",
            )
    tip_rows = tip_s2_pi = None
    if tip_projection:
        comps = _projected_ball_comparands(
            sub, q, V1, tips, a_1g, "i32b_" + label.replace(" ", "_")
        )
        tip_rows = {
            "dt1": _tip_slope_rows({n: tip_m[n]["th1"] for n in tips}, DT_1_BOX, comps, RHO_C),
            "dt1/10": _tip_slope_rows(
                {n: tip_m[n]["th10"] for n in tips}, DT_1_BOX / 10.0, comps, RHO_C
            ),
        }
        tip_s2_pi = {n: tip_rows["dt1/10"][n]["R_Pi"] for n in tips}
        for n in tips:
            a, b = tip_rows["dt1"][n], tip_rows["dt1/10"][n]
            _say(
                comm,
                f"[{label}] tip ball {n}, sense {SENSE}: <q> {comps[n]['q']:.7e} <Pq> "
                f"{comps[n]['pq']:.7e} W/m^3, <Pq>/<q> - 1 {a['pq_q']:+.4e}; R_Pi dt1 "
                f"{a['R_Pi']:.4e}, R_Pi dt1/10 {b['R_Pi']:.4e} ((S2-tip) gate, band "
                f"{S2_BAND:.0e}); R_DG0 dt1 {a['R_DG0']:.4e} (retired form, printed)",
            )
    labels = list(probes)
    probe_th, probe_ok = _probe_theta(
        res.temperature, [probes[n] for n in labels], (msh, parent_ids)
    )
    for n, v, valid in zip(labels, probe_th, probe_ok):
        _say(
            comm,
            f"[{label}] probe {n} {tuple(round(x, 6) for x in probes[n])}, sense {SENSE}: "
            f"theta(t*) {v:+.9e} K (valid {bool(valid)}; printed)",
        )
    t_total = time.perf_counter() - t0
    _say(
        comm,
        f"[{label}] CN N={N_MAIN} to t* {T_STAR_BOX:.3f} s, dt {dt:.6f} s (T_STAR_BOX / N_MAIN; "
        f"BE dt_1 {DT_1_BOX:.6f} s): (E) balance CN {res.max_balance_residual:.3e} BE "
        f"{res1.max_balance_residual:.3e} (band {E_BAND:.0e}, ASSERTED); (E') stored {stor:.12e} "
        f"+ surface {conv:.12e} W: rel {ep:.3e} (band {E_BAND:.0e}, ASSERTED); (S2) {err_S2:.4e} "
        f"(band {S2_BAND:.0e}, ASSERTED); max|Im theta| {imag:.3e}; source {t_source:.1f} s, CN "
        f"{t_cn:.1f} s, arm {t_total:.1f} s",
    )
    out = {
        "label": label,
        "G": G,
        "P": P,
        "B": err_B,
        "Ep": ep,
        "S2": err_S2,
        "dt": dt,
        "E": max(res.max_balance_residual, res1.max_balance_residual),
        "counts": counts,
        "merged": merged,
        "tip_th": tip_th,
        "tip_s2": tip_s2,
        "probe_ok": probe_ok,
        "t": t_total,
        "probe_th": dict(zip(labels, probe_th)),
    }
    if tip_projection:
        out["E"] = max(out["E"], res10.max_balance_residual)
        out["E10"] = res10.max_balance_residual
        out["tip_rows"] = tip_rows
        out["tip_s2_pi"] = tip_s2_pi
        out["comparands"] = comps
    if f2182_th is not None:
        out["f2182_th"] = f2182_th
    return out


def _thermal_failures(arm):
    """Every thermal gate of one arm, as names (collected, item 668)."""
    fails = []
    tag = arm["label"]
    if not arm["B"] <= B_BAND:
        fails.append(f"{tag} (B) {arm['B']:.3e}")
    if not arm["E"] <= E_BAND:
        fails.append(f"{tag} (E) {arm['E']:.3e}")
    if not arm["Ep"] <= E_BAND:
        fails.append(f"{tag} (E') {arm['Ep']:.3e}")
    if not arm["S2"] <= S2_BAND:
        fails.append(f"{tag} (S2) {arm['S2']:.3e}")
    if arm["merged"]:
        if arm["counts"]["wall"] != 0:
            fails.append(f"{tag} wall facets {arm['counts']['wall']} on the merged arm")
    else:
        if not arm["counts"]["wall"] > 0:
            fails.append(f"{tag} no adiabatic wall on the rod arm")
        # I32b: the projected comparand (R_Π at Δt₁/10) when the arm carries it.
        worst = max(arm.get("tip_s2_pi", arm["tip_s2"]).values())
        if not worst <= S2_BAND:
            fails.append(f"{tag} (S2-tip) {worst:.3e}")
    return fails


def _em_failures(comm, arm, feeds):
    """(R1), (P), NEG and the kinds on one arm, as `_run_drive` asserts them
    (`test_geo36_rod_f_clinical_drive.py:561–583`)."""
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_geo38_f_clinical_matrix import BLIND_FACTOR
    from tests.validation.test_port_birdcage_ring_column import _reciprocity_ratio
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    fails = []
    tag = arm["label"]
    recip = _reciprocity_ratio(arm["s2"])
    _say(
        comm,
        f"arm {tag}: kinds {arm['kinds']}, (R1) {recip:.3e} (band {RECIPROCITY_BAND:.0e}, "
        f"ASSERTED); column times {', '.join(f'{t:.2f}' for t in arm['times'])} s",
    )
    if arm["kinds"] != ["solve", "held"]:
        fails.append(f"{tag} kinds {arm['kinds']}")
    if not recip <= RECIPROCITY_BAND:
        fails.append(f"{tag} (R1) {recip:.3e}")
    for j, p in enumerate(arm["parts"]):
        rel, blind = _partition(p)
        _say(
            comm,
            f"  (P) column {feeds[j]}: rel {rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); "
            f"NEG blind {blind:.4e} (>= {BLIND_FACTOR:g}x band, ASSERTED)",
        )
        if not rel <= DISCRETE_IDENTITY_RTOL:
            fails.append(f"{tag} (P) {feeds[j]} {rel:.3e}")
        if not blind >= BLIND_FACTOR * DISCRETE_IDENTITY_RTOL:
            fails.append(f"{tag} NEG {feeds[j]} {blind:.3e}")
        if not p["p_surf"] > 0.0:
            fails.append(f"{tag} P_surf {p['p_surf']!r}")
    return fails


def _s0_failures(census_rod, census_twin, merged_gel, merged_rod):
    """(S0): `_a0m_failures`, plus the merged arm's tag 3 = the rod arm's 3 + 403."""
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG

    fails = _a0m_failures(census_rod, census_twin)
    expected = census_rod.get(PHANTOM_CELL_TAG, 0) + census_rod.get(BIRDCAGE_IMPLANT_CELL_TAG, 0)
    if merged_gel != expected:
        fails.append(f"merged tag 3 {merged_gel} vs rod arm's 3 + 403 {expected}")
    if merged_rod != 0:
        fails.append(f"merged tag 403 {merged_rod}")
    return fails


# ---------------------------------------------------------------------------
# I87' (§9 item 968): the ``ROD_THERMAL`` knob -- rod cells inside the thermal
# domain with I66's coefficients, printed after the exclusion arm.  Unset (the
# default) leaves `_run_heat` exactly as before and prints nothing.
# ---------------------------------------------------------------------------
ROD_THERMAL_ENV = "FEM_EM_GEO36_S8_ROD_THERMAL"
ROD_THERMAL_ARMS = ("gel", "ti", "hole")


def _rod_thermal():
    """Unset or empty => ``None``; one of `ROD_THERMAL_ARMS` => itself; else ``ValueError``."""
    raw = os.environ.get(ROD_THERMAL_ENV, "").strip()
    if raw == "":
        return None
    if raw in ROD_THERMAL_ARMS:
        return raw
    raise ValueError(f"{ROD_THERMAL_ENV} must be unset or one of {ROD_THERMAL_ARMS}, got {raw!r}")


def _rod_coefficients(sub, parent_ids, cell_tags, arm):
    """DG0 ``k``, ``rho_c`` and the rod indicator on the thermal submesh, plus the
    owned rod-cell count (allreduced).  Rod cells: parent tag 403 on the original
    tags.  I66's rule (`test_imp8_rod_conducting_body_fem.py:_coefficients`)."""
    # Function-level imports: the IMP-8 modules import this one (the cycle).
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from tests.validation.test_imp8_composite_cylinder_closed_form import ROD_K, ROD_RHO_C
    from tests.validation.test_imp8_rod_conducting_body_fem import EPS_HOLE

    if arm not in ROD_THERMAL_ARMS:
        raise ValueError(f"rod thermal arm {arm!r} not in {ROD_THERMAL_ARMS}")
    pmsh_tdim = sub.topology.dim
    pimap = cell_tags.topology.index_map(pmsh_tdim)
    parent_values = np.full(pimap.size_local + pimap.num_ghosts, -1, dtype=np.int32)
    parent_values[cell_tags.indices] = cell_tags.values
    rod = parent_values[parent_ids] == BIRDCAGE_IMPLANT_CELL_TAG
    simap = sub.topology.index_map(pmsh_tdim)
    n_all = simap.size_local + simap.num_ghosts
    W = fem.functionspace(sub, ("DG", 0))
    dofs = np.asarray(W.dofmap.list).reshape(-1)[:n_all]
    k, rho_c, ind = fem.Function(W), fem.Function(W), fem.Function(W)
    k.x.array[:] = K_GEL
    rho_c.x.array[:] = RHO_C
    ind.x.array[:] = 0.0
    idx = np.flatnonzero(rod[:n_all])
    if arm == "ti":
        k.x.array[dofs[idx]] = ROD_K
        rho_c.x.array[dofs[idx]] = ROD_RHO_C
    elif arm == "hole":
        k.x.array[dofs[idx]] = EPS_HOLE * K_GEL
        rho_c.x.array[dofs[idx]] = EPS_HOLE * RHO_C
    ind.x.array[dofs[idx]] = 1.0
    for f in (k, rho_c, ind):
        f.x.scatter_forward()
    n_rod = sub.comm.allreduce(int(np.count_nonzero(rod[: simap.size_local])), op=MPI.SUM)
    return k, rho_c, ind, n_rod


def _thermal_arm_rod_cells(
    label, msh, cell_tags, e_sup, sigma_f, G, centre, tips, probes, *, rod_thermal, f2182=None
):
    """`_thermal_arm` (tip_projection on, merged domain) with the rod cells inside the
    thermal domain at I66's DG0 coefficients.  P and q on the original tags (q = 0 on
    rod cells), so G*P is the exclusion arm's.  Never calls `_transient`."""
    from fem_em_solver.thermal import solve_transient_conduction

    comm = msh.comm
    t0 = time.perf_counter()
    P = float(
        np.real(
            mean_sar(
                e_sup,
                sigma=sigma_f,
                rho=PHANTOM_RHO_KG_PER_M3,
                cell_tags=cell_tags,
                comm=comm,
                subdomain_ids=PHANTOM_CELL_TAG,
            )["dissipated_power_w"]
        )
    )
    if G is None:
        G = 1.0 / P
    dep = G * P
    merged = _merge_rod_into_phantom(msh, cell_tags)
    sub, parent_ids = _tag_submesh(msh, merged, PHANTOM_CELL_TAG)
    ft, counts = _outer_robin_tags(sub, centre)
    one = fem.Constant(sub, default_scalar_type(1.0))
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)
    v_sub = float(np.real(_allsum(comm, one * ufl.dx(domain=sub))))
    area = float(np.real(_allsum(comm, one * ds)))
    k, rho_c, rod_ind, n_rod = _rod_coefficients(sub, parent_ids, cell_tags, rod_thermal)
    q = _heat_source(msh, cell_tags, e_sup, sigma_f, G, sub, parent_ids)
    int_q = float(np.real(_allsum(comm, q * ufl.dx(domain=sub))))
    err_B = abs(int_q - dep) / abs(dep)
    q_m = _heat_source(msh, merged, e_sup, sigma_f, G, sub, parent_ids)
    int_qm = float(np.real(_allsum(comm, q_m * ufl.dx(domain=sub))))
    n_src = abs(int_qm - dep) / abs(dep)
    del q_m
    t_source = time.perf_counter() - t0
    _say(
        comm,
        f"[{label}] sense {SENSE}: P_gel {P:.12e} W; G {G:.12e} 1/W; G*P {dep:.12e} W; V_sub "
        f"{v_sub:.12e} m^3; Robin area {area:.9e} m^2; facets exterior {counts['exterior']} = "
        f"outer {counts['outer']} + rod wall {counts['wall']} (rod cells in the domain, "
        f"{rod_thermal}; n_rod {n_rod})",
    )
    _say(
        comm,
        f"[{label}] (B) int_sub q {int_q:.15e} vs G*P {dep:.15e}: rel {err_B:.3e} "
        f"(band {B_BAND:.0e}, ASSERTED); (N-src) source on the merged tags: int q_m "
        f"{int_qm:.15e}, rel {n_src:.6e} (> band, ASSERTED; predicted 108/10260 on the stub)",
    )

    def _solve(dt, n_steps, theta, startup_half_steps=0, snapshot_steps=()):
        return solve_transient_conduction(
            sub,
            ft,
            k=k,
            rho_c=rho_c,
            q=q,
            h=H_GEL,
            t_inf=T_INF,
            t_init=T_INF,
            robin_tag=ROBIN_TAG,
            dt=dt,
            n_steps=n_steps,
            theta=theta,
            startup_half_steps=startup_half_steps,
            degree=1,
            snapshot_steps=snapshot_steps,
        )

    f2182_th = None
    if f2182:
        res900 = _solve(F2182_T_S / N_MAIN, N_MAIN, 0.5, startup_half_steps=2)
        f_names = list(f2182)
        f_th, f_ok = _probe_theta(
            res900.temperature, [f2182[n] for n in f_names], (msh, parent_ids)
        )
        f2182_th = dict(zip(f_names, f_th))
        for n, v, valid in zip(f_names, f_th, f_ok):
            _say(
                comm,
                f"[{label}] F2182-shaped probe {n} {tuple(round(float(x), 6) for x in f2182[n])} "
                f"({F2182_OFFSET_M * 1e3:.1f} mm beyond the apex), sense {SENSE}: theta"
                f"({F2182_T_S:.0f} s) {v:+.9e} K (valid {bool(valid)}; placement and time "
                f"quoted from memory, unconfirmed; W3; printed)",
            )
        del res900

    dt = T_STAR_BOX / N_MAIN
    t1 = time.perf_counter()
    res = _solve(dt, N_MAIN, 0.5, startup_half_steps=2, snapshot_steps=(N_MAIN - 1, N_MAIN))
    t_cn = time.perf_counter() - t1
    tn, tm = res.snapshots[N_MAIN], res.snapshots[N_MAIN - 1]
    stor_rc = float(np.real(_allsum(comm, rho_c * (tn - tm) / dt * ufl.dx(domain=sub))))
    conv = float(np.real(_allsum(comm, H_GEL * (0.5 * (tn + tm) - T_INF) * ds)))
    ep_rc = abs(stor_rc + conv - dep) / abs(dep)
    ep_sc, stor_sc, _conv_sc = _e_prime_box(sub, ds, res, dep, dt)
    res1 = _solve(DT_1_BOX, 1, 1.0)
    int_th1 = float(np.real(_allsum(comm, (res1.temperature - T_INF) * ufl.dx(domain=sub))))
    int_rc_th1 = float(
        np.real(_allsum(comm, rho_c * (res1.temperature - T_INF) * ufl.dx(domain=sub)))
    )
    s2_rc = abs(int_rc_th1 / (DT_1_BOX * dep) - 1.0)
    s2_sc = abs(RHO_C * int_th1 / (DT_1_BOX * dep) - 1.0)
    res10 = _solve(DT_1_BOX / 10.0, 1, 1.0)
    imag = max(
        _max_imag(res.temperature), _max_imag(res1.temperature), _max_imag(res10.temperature)
    )

    V1 = res.temperature.function_space
    th, th1, th10 = fem.Function(V1), fem.Function(V1), fem.Function(V1)
    th.x.array[:] = res.temperature.x.array - T_INF
    th1.x.array[:] = res1.temperature.x.array - T_INF
    th10.x.array[:] = res10.temperature.x.array - T_INF
    for f in (th, th1, th10):
        f.x.scatter_forward()
    fields = {
        "th": th,
        "th1": th1,
        "q": q,
        "th10": th10,
        "rc_th1": rho_c * th1,
        "rc_th10": rho_c * th10,
        "rod": rod_ind,
    }
    read = _ball_reader(sub, fields)
    a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    tip_th, tip_s2, tip_m, v_rod_ball = {}, {}, {}, {}
    for name, c in tips.items():
        vol, m = read(np.asarray(c, dtype=float), a_1g)
        tip_m[name] = m
        tip_th[name] = m["th"]
        tip_s2[name] = abs(RHO_C * m["th1"] / (DT_1_BOX * m["q"]) - 1.0)
        v_rod_ball[name] = vol * m["rod"]
        _say(
            comm,
            f"[{label}] 1 g ball (a {a_1g * 1e3:.4f} mm) at {name}, sense {SENSE}: theta(t*) "
            f"{m['th']:.9e} K, ball volume in the thermal domain {vol:.6e} m^3; V(B cap rod) "
            f"{v_rod_ball[name]:.6e} m^3; retired form (S2-tip, DG0) {tip_s2[name]:.4e} (band "
            f"{S2_BAND:.0e}, printed)",
        )
    comps = _projected_ball_comparands(sub, q, V1, tips, a_1g, "i32b_" + label.replace(" ", "_"))
    rows_rc = _tip_slope_rows({n: tip_m[n]["rc_th10"] for n in tips}, DT_1_BOX / 10.0, comps, 1.0)
    rows_sc = _tip_slope_rows({n: tip_m[n]["th10"] for n in tips}, DT_1_BOX / 10.0, comps, RHO_C)
    tip_s2_pi = {n: rows_rc[n]["R_Pi"] for n in tips}
    tip_s2_pi_sc = {n: rows_sc[n]["R_Pi"] for n in tips}
    for n in tips:
        _say(
            comm,
            f"[{label}] tip ball {n}, sense {SENSE}: <q> {comps[n]['q']:.7e} <Pq> "
            f"{comps[n]['pq']:.7e} W/m^3, <Pq>/<q> - 1 {rows_rc[n]['pq_q']:+.4e}; R_Pi dt1/10 "
            f"rho_c form {tip_s2_pi[n]:.4e} ((S2-tip) gate, band {S2_BAND:.0e}); scalar form "
            f"{tip_s2_pi_sc[n]:.4e} (printed)",
        )
    labels = list(probes)
    probe_th, probe_ok = _probe_theta(
        res.temperature, [probes[n] for n in labels], (msh, parent_ids)
    )
    for n, v, valid in zip(labels, probe_th, probe_ok):
        _say(
            comm,
            f"[{label}] probe {n} {tuple(round(x, 6) for x in probes[n])}, sense {SENSE}: "
            f"theta(t*) {v:+.9e} K (valid {bool(valid)}; printed)",
        )
    t_total = time.perf_counter() - t0
    e_all = max(res.max_balance_residual, res1.max_balance_residual, res10.max_balance_residual)
    _say(
        comm,
        f"[{label}] CN N={N_MAIN} to t* {T_STAR_BOX:.3f} s, dt {dt:.6f} s (T_STAR_BOX / N_MAIN; "
        f"BE dt_1 {DT_1_BOX:.6f} s): (E) balance CN {res.max_balance_residual:.3e} BE "
        f"{res1.max_balance_residual:.3e} BE/10 {res10.max_balance_residual:.3e} (band "
        f"{E_BAND:.0e}, ASSERTED); (E') rho_c form stored {stor_rc:.12e} + surface {conv:.12e} W: "
        f"rel {ep_rc:.3e} (band {E_BAND:.0e}, ASSERTED), scalar form stored {stor_sc:.12e}: rel "
        f"{ep_sc:.3e} (printed); (S2) rho_c form {s2_rc:.4e} (band {S2_BAND:.0e}, ASSERTED), "
        f"scalar form {s2_sc:.4e} (printed); max|Im theta| {imag:.3e}; source {t_source:.1f} s, "
        f"CN {t_cn:.1f} s, arm {t_total:.1f} s",
    )
    out = {
        "label": label,
        "G": G,
        "P": P,
        "B": err_B,
        "Ep": ep_rc,
        "S2": s2_rc,
        "dt": dt,
        "E": e_all,
        "counts": counts,
        "merged": True,
        "tip_th": tip_th,
        "tip_s2": tip_s2,
        "probe_ok": probe_ok,
        "t": t_total,
        "probe_th": dict(zip(labels, probe_th)),
        "E10": res10.max_balance_residual,
        "tip_rows": {"dt1/10": rows_rc, "dt1/10 scalar": rows_sc},
        "comparands": comps,
        "Ep_rc": ep_rc,
        "Ep_sc": ep_sc,
        "S2_rc": s2_rc,
        "S2_sc": s2_sc,
        "tip_s2_pi": tip_s2_pi,
        "tip_s2_pi_sc": tip_s2_pi_sc,
        "v_rod_ball": v_rod_ball,
        "n_rod": n_rod,
        "n_src": n_src,
        "v_sub": v_sub,
        "rod_thermal": rod_thermal,
    }
    if f2182_th is not None:
        out["f2182_th"] = f2182_th
    return out


def _rod_cells_failures(arm):
    """The in-domain rod arm's gates, as names: (B), (E), (E')rho_c, (S2)rho_c, no
    wall, and the rho_c-form (S2-tip) at Δt₁/10."""
    fails = []
    tag = arm["label"]
    if not arm["B"] <= B_BAND:
        fails.append(f"{tag} (B) {arm['B']:.3e}")
    if not arm["E"] <= E_BAND:
        fails.append(f"{tag} (E) {arm['E']:.3e}")
    if not arm["Ep_rc"] <= E_BAND:
        fails.append(f"{tag} (E')rho_c {arm['Ep_rc']:.3e}")
    if not arm["S2_rc"] <= S2_BAND:
        fails.append(f"{tag} (S2)rho_c {arm['S2_rc']:.3e}")
    if arm["counts"]["wall"] != 0:
        fails.append(f"{tag} wall facets {arm['counts']['wall']}")
    worst = max(arm["tip_s2_pi"].values())
    if not worst <= S2_BAND:
        fails.append(f"{tag} (S2-tip)rho_c {worst:.3e}")
    return fails


def _run_heat(comm, monkeypatch, c_pf):
    """I32's body: both arms' columns on one mesh, both solvers freed, then the
    gel twin's thermal arm (it sets G) and the rod's.  Gates asserted at the end."""
    from fem_em_solver.core import TimeHarmonicSolver
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo38_f_clinical_column as col
    from tests.validation import test_geo38_f_clinical_esr_root as c4
    from tests.validation import test_wf10_f_human_copper_esr_root as h17
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM
    from tests.validation.test_wf12_f_human_hole import _rss_gib

    t_window = time.perf_counter()
    _box_preconditions()
    rt = _rod_thermal()
    if rt is not None:
        _say(
            comm,
            f"rod thermal: {rt} (rod cells in the thermal domain, I66's coefficients; printed "
            f"after the exclusion arm)",
        )
    z0 = float(REFERENCE_IMPEDANCE_OHM)
    _say(
        comm,
        f"\nstart: -n {comm.size}, I32 body; node C {c_pf:.4f} pF, degree 1, sense {SENSE}; rod "
        f"(tag 403 sigma_wire, kept cells) then gel twin (403 -> gel, same mesh)",
    )
    records = _route_drive(monkeypatch)
    monkeypatch.delenv(c4.C4_DEGREE_ENV, raising=False)
    fixture = c4._build_c4_fixture(comm)  # (A0) / census / (U1) / (D-rod) / (F) asserted inside
    assert fixture["degree"] == 1, fixture["degree"]
    assert len(records) == 1 and records[0][0] == 1, [r[0] for r in records]
    fails = list(_m2_failures([(1, records[0][1]), (2, records[0][1])]))
    feeds, ids = fixture["feeds"], fixture["ids"]
    ctx = fixture["built"]["ctx"]
    env = _drive_env(fixture, comm)
    frequency_hz = env["omega"] / (2.0 * np.pi)
    others = sorted((pid for pid in ids if pid not in feeds), key=lambda p: int(p[1:]))
    order = list(feeds) + others

    c_f = c_pf * 1e-12
    base_specs = list(ctx["specs"])
    base_ids = [sp.port_id for sp in base_specs]
    specs, z_rc, z_th = c4._esr_specs(base_specs, feeds, c_f, frequency_hz, z0)
    z_fails, im_worst, th_worst = c4._z_failures(specs, base_ids, feeds, c_f, frequency_hz, z0)
    fails += list(z_fails)
    v_s = h17._source_voltage(z_rc, z0)
    _say(
        comm,
        f"(M2) / (Z) node {c_pf:.4f} pF: Z_RC {z_rc:.9e} Ohm, Z_th {z_th:.9e} Ohm, V_s {v_s:.9e} "
        f"V; Im rel {im_worst:.1e}, Z_th rel {th_worst:.1e}; failures {fails or 'none'} (ASSERTED)",
    )

    msh, cell_tags = env["msh"], env["cell_tags"]
    centre = i35._rod_centre()
    tips_all, _t1 = _probe_points(centre)
    tips = {n: c for n, c in tips_all.items() if not n.endswith("_apex")}
    probes = _p37_fixture(centre[1])

    t_em = time.perf_counter()
    solver_rod = ctx["solver"]
    problem_rod = solver_rod.problem
    census_rod = col._owned_tag_counts(msh, problem_rod.cell_tags, comm)
    arm_rod = _heat_columns(env, solver_rod, specs, feeds, order, z_th, v_s)
    arm_rod["label"] = "rod"
    del solver_rod
    del ctx["solver"]
    gc.collect()

    problem_twin = dataclasses.replace(
        problem_rod, material_map=_gel_twin_map(dict(problem_rod.material_map))
    )
    census_twin = col._owned_tag_counts(msh, problem_twin.cell_tags, comm)
    solver_twin = TimeHarmonicSolver(problem_twin, degree=1)
    arm_twin = _heat_columns(env, solver_twin, specs, feeds, order, z_th, v_s)
    arm_twin["label"] = "gel twin"
    del solver_twin
    gc.collect()
    t_em = time.perf_counter() - t_em

    fails += _em_failures(comm, arm_rod, feeds)
    fails += _em_failures(comm, arm_twin, feeds)
    p_tot, p_gel, p_rod = _vol_powers(env, arm_rod["sigma"], arm_rod["e_sup"])
    _say(
        comm,
        f"rod share of P_vol, sense {SENSE}: {p_rod / p_tot:.6e} (P_vol {p_tot:.9e} W, gel "
        f"{p_gel:.9e}, rod {p_rod:.9e}; printed)",
    )

    merged = _merge_rod_into_phantom(msh, cell_tags)
    merged_gel = _owned_tag_count(msh, merged, PHANTOM_CELL_TAG)
    merged_rod = _owned_tag_count(msh, merged, BIRDCAGE_IMPLANT_CELL_TAG)
    s0_fails = _s0_failures(census_rod, census_twin, merged_gel, merged_rod)
    fails += s0_fails
    _say(
        comm,
        f"(S0) census rod {census_rod} vs twin {census_twin}; merged tag 3 {merged_gel}, tag 403 "
        f"{merged_rod}: failures {s0_fails or 'none'} (ASSERTED exact)",
    )

    f2182 = {}
    for n in tips:
        if f"{n}_apex" not in tips_all:
            _say(comm, f"F2182-shaped probe at tip {n}: no apex on this fixture (skipped)")
            continue
        apex = np.asarray(tips_all[f"{n}_apex"], dtype=float)
        axis = apex - np.asarray(tips_all[n], dtype=float)
        f2182[n] = tuple(apex + F2182_OFFSET_M * axis / np.linalg.norm(axis))
    _say(
        comm,
        f"F2182-shaped read: t {F2182_T_S:.0f} s, offset {F2182_OFFSET_M * 1e3:.1f} mm beyond each "
        f"apex (placement and time quoted from memory, unconfirmed; W3; printed, never asserted)",
    )

    heat_twin = _thermal_arm(
        "gel twin",
        msh,
        merged,
        arm_twin.pop("e_sup"),
        arm_twin.pop("sigma"),
        None,
        centre,
        tips,
        probes,
        merged=True,
        tip_projection=True,
        f2182=f2182,
    )
    rod_fields = (arm_rod["e_sup"], arm_rod["sigma"]) if rt is not None else None
    heat_rod = _thermal_arm(
        "rod",
        msh,
        cell_tags,
        arm_rod.pop("e_sup"),
        arm_rod.pop("sigma"),
        heat_twin["G"],
        centre,
        tips,
        probes,
        merged=False,
        tip_projection=True,
        f2182=f2182,
    )
    fails += _thermal_failures(heat_twin) + _thermal_failures(heat_rod)
    if rt is not None:
        heat_rod_cells = _thermal_arm_rod_cells(
            f"rod {rt}",
            msh,
            cell_tags,
            *rod_fields,
            heat_twin["G"],
            centre,
            tips,
            probes,
            rod_thermal=rt,
            f2182=f2182,
        )
        fails += _rod_cells_failures(heat_rod_cells)
        _say(
            comm,
            f"rod thermal {rt} over the exclusion arm - 1, sense {SENSE}: "
            + ", ".join(
                f"{n} theta(t*) {heat_rod_cells['tip_th'][n] / heat_rod['tip_th'][n] - 1.0:+.6e}"
                for n in tips
            )
            + "; "
            + ", ".join(
                f"{n} F2182 theta({F2182_T_S:.0f} s) "
                f"{heat_rod_cells['f2182_th'][n] / heat_rod['f2182_th'][n] - 1.0:+.6e}"
                for n in heat_rod_cells.get("f2182_th", {})
            )
            + " (printed)",
        )

    ratios = {n: heat_rod["tip_th"][n] / heat_twin["tip_th"][n] for n in tips}
    _say(
        comm,
        f"tip theta(t*) ratio rod / gel twin, sense {SENSE}: "
        + ", ".join(f"{n} {v:.6f}" for n, v in ratios.items())
        + f" beside I41's {I41_TIP_RATIO} (F-small, ccw; printed)",
    )
    comm.Barrier()
    _m, per, tot = _rss_gib(comm)
    _say(
        comm,
        f"price: build {fixture['t_build']:.1f} s; EM {t_em:.1f} s; thermal twin "
        f"{heat_twin['t']:.1f} s, rod {heat_rod['t']:.1f} s; RSS per rank "
        f"{', '.join(f'{r:.2f}' for r in (per or []))} GiB, summed {tot:.3f} GiB; window "
        f"{time.perf_counter() - t_window:.1f} s",
    )
    _say(comm, f"gates: failures {fails or 'none'} (ASSERTED)")
    assert not fails, fails
    if rt is not None:
        return {"rod": heat_rod, "twin": heat_twin, "ratios": ratios, "rod_cells": heat_rod_cells}
    return {"rod": heat_rod, "twin": heat_twin, "ratios": ratios}


# ---------------------------------------------------------------------------
# The gated test.
# ---------------------------------------------------------------------------
@complex_only
def test_geo36_i32_rod_heating(monkeypatch):
    mode = _heat_mode()
    c_pf = _node_knob(mode)  # missing under 1 => ValueError, before gmsh
    if mode is None:
        pytest.skip(f"{HEAT_ENV} unset: GEO-36 I32 is selected by environment only")
    _guard_foreign()
    _refuse_degree()
    _run_heat(MPI.COMM_WORLD, monkeypatch, c_pf)


# ---------------------------------------------------------------------------
# Coding tests (no gmsh, stubbed `_drive`); skip whenever HEAT_ENV is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(
    os.environ.get(HEAT_ENV, "").strip() != "",
    reason=f"{HEAT_ENV} set: coding tests skip in windows",
)
STUB_E0 = 3.0
STUB_SIGMA = 0.5
STUB_WALL_SELECT_M = 0.03  # the stub's jagged 403 block reaches ≈ 24 mm; the box face is at 50 mm
STUB_CENTRE = (0.0, 0.0, 0.0)
STUB_PROBES = {"interior": (0.03, 0.03, 0.0), "face": (0.05, 0.02, 0.01)}


def _stubs(m, fx, json_dir):
    """`_k9_stubs` (`_run_drive`'s surroundings on `_k8_arm`'s mesh), with the
    stubbed names rebound here: this module imports them by name.  Rod columns
    are E_u, gel-twin columns 0.5 E_u.  Returns the recorded prints."""
    import tests.validation.test_geo36_rod_f_clinical as i35

    this = sys.modules[_heat_columns.__module__]
    drive_module = sys.modules[_drive.__module__]
    _k9_stubs(m, fx, json_dir)
    for name in ("_route_drive", "_drive_env", "_drive", "_probe_points"):
        m.setattr(this, name, getattr(drive_module, name))
    stub_points = drive_module._probe_points

    def _points_with_apex(centre):
        # I32b: planted a_apex / b_apex = tip + (0, 0, 1 mm), so K10-R′ runs the F2182 branch.
        tips, t1 = stub_points(centre)
        tips = dict(tips)
        for n in ("a", "b"):
            if n in tips:
                tips[f"{n}_apex"] = tuple(np.asarray(tips[n], dtype=float) + (0.0, 0.0, 1.0e-3))
        return tips, t1

    m.setattr(this, "_probe_points", _points_with_apex)
    m.setattr(this, "WALL_SELECT_M", STUB_WALL_SELECT_M)
    m.setattr(this, "_p37_fixture", lambda _off: STUB_PROBES)
    m.setattr(i35, "_rod_centre", lambda: STUB_CENTRE)
    lines = []
    real_say = this._say

    def _say_recorded(comm_, msg):
        lines.append(msg)
        real_say(comm_, msg)

    m.setattr(this, "_say", _say_recorded)
    return lines


def _columns_e_sup(m, fx, json_dir, *, second_zero=False):
    """``e_sup`` and σ from `_heat_columns` on the stub, optionally with e₂ = 0."""
    from tests.validation import test_wf10_f_human_copper_esr_root as h17

    this = sys.modules[_heat_columns.__module__]
    _stubs(m, fx, json_dir)
    stubbed = this._drive

    def _column(env, solver, specs, pid, *, held):
        d = stubbed(env, solver, specs, pid, held=held)
        if held and second_zero:
            d["e"].x.array[:] = 0.0
        return d

    m.setattr(this, "_drive", _column)
    z_rc = h17._z_rc(193.5996e-12, 64.0e6)
    z_th = h17._thevenin_feed_impedance(z_rc, 50.0)
    v_s = h17._source_voltage(z_rc, 50.0)
    env = this._drive_env(None, fx[0].comm)
    out = _heat_columns(env, object(), [], ["P1", "P2"], ["P1", "P2", "P3", "P4"], z_th, v_s)
    return out["e_sup"], out["sigma"], out


def _merged_q(msh, cell_tags, e_sup, sigma_f, G=1.0):
    merged = _merge_rod_into_phantom(msh, cell_tags)
    sub, parent_ids = _tag_submesh(msh, merged, PHANTOM_CELL_TAG)
    return merged, sub, parent_ids, _heat_source(msh, merged, e_sup, sigma_f, G, sub, parent_ids)


def _owned_extrema(q):
    comm = q.function_space.mesh.comm
    n = q.function_space.dofmap.index_map.size_local
    vals = np.real(q.x.array[:n])
    return (
        comm.allreduce(float(vals.min()) if n else np.inf, op=MPI.MIN),
        comm.allreduce(float(vals.max()) if n else -np.inf, op=MPI.MAX),
    )


@_coding
def test_i32a_k0_knobs_and_box_preconditions(monkeypatch):
    monkeypatch.delenv(HEAT_ENV, raising=False)
    assert _heat_mode() is None
    monkeypatch.setenv(HEAT_ENV, "1")
    assert _heat_mode() == "1"
    for bad in ("0", "build", "x"):
        monkeypatch.setenv(HEAT_ENV, bad)
        with pytest.raises(ValueError):
            _heat_mode()
    monkeypatch.delenv(S8_DRIVE_DEGREE_ENV, raising=False)
    _refuse_degree()
    monkeypatch.setenv(S8_DRIVE_DEGREE_ENV, "2")
    with pytest.raises(ValueError):
        _refuse_degree()
    _box_preconditions()
    print(
        f"\n[GEO-36 I32a K0] box time base: T_STAR {T_STAR_BOX:.3f} s, dt "
        f"{T_STAR_BOX / N_MAIN:.6f} s, dt_1 {DT_1_BOX:.6f} s",
        flush=True,
    )
    assert (
        abs(T_STAR_BOX - 2821.5) <= 1.0
    ), T_STAR_BOX  # the brief's 2 821.5 s, MAT-15's would be 1 254 s


@complex_only
@_coding
def test_i32a_k10q_source_closed_form(monkeypatch, tmp_path):
    """(K10-q): e_sup = (1 + j) E₀ x̂, so ½σ|e_sup|² = σE₀² on every cell, and q/ρ is
    the SAR `mass_averaged_sar` reads.  Control: e₂ = 0 misses the closed form by 0.5."""
    from fem_em_solver.post.sar import build_density_field, mass_averaged_sar
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL

    comm = MPI.COMM_WORLD
    fx = _k9_fixture(comm)
    msh, cell_tags = fx[0], fx[1]
    root = comm.bcast(str(tmp_path), root=0)
    q_exact = STUB_SIGMA * STUB_E0**2
    with monkeypatch.context() as m:
        e_sup, sigma_f, out = _columns_e_sup(m, fx, root)
        _, _sub, _ids, q = _merged_q(msh, cell_tags, e_sup, sigma_f)
        lo, hi = _owned_extrema(q)
        a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
        rho_field = build_density_field(msh, PHANTOM_RHO_KG_PER_M3)
        ball = mass_averaged_sar(
            e_sup,
            sigma=sigma_f,
            rho=rho_field,
            center=fx[4]["a"],
            radius=a_1g,
            quadrature_degree=QUADRATURE_DEGREE,
        )
    miss = max(abs(lo / q_exact - 1.0), abs(hi / q_exact - 1.0))
    r_sar = abs(float(ball["averaged_sar_w_per_kg"]) * PHANTOM_RHO_KG_PER_M3 / q_exact - 1.0)
    with monkeypatch.context() as m:
        e_one, sigma_one, _ = _columns_e_sup(m, fx, root, second_zero=True)
        _, _sub, _ids, q_one = _merged_q(msh, cell_tags, e_one, sigma_one)
        lo_c, hi_c = _owned_extrema(q_one)
    miss_c = max(abs(lo_c / q_exact - 1.0), abs(hi_c / q_exact - 1.0))
    if comm.rank == 0:
        print(
            f"\n[GEO-36 I32a K10-q] q on [{lo:.12e}, {hi:.12e}] vs sigma E0^2 {q_exact:.12e}: miss "
            f"{miss:.3e}; SAR*rho miss {r_sar:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); kinds "
            f"{out['kinds']}; control e2 = 0 miss {miss_c:.12f} (predicted {SINGLE_COLUMN_MISS}, "
            f"within {SINGLE_COLUMN_CEILING:g}, ASSERTED)",
            flush=True,
        )
    assert miss <= DISCRETE_IDENTITY_RTOL, miss
    assert r_sar <= DISCRETE_IDENTITY_RTOL, r_sar
    assert out["kinds"] == ["solve", "held"], out["kinds"]
    assert abs(miss_c - SINGLE_COLUMN_MISS) <= SINGLE_COLUMN_CEILING, miss_c


@complex_only
@_coding
def test_i32a_k10es_energy_and_slope_on_uniform_source(monkeypatch, tmp_path):
    """(K10-E) / (K10-S): `_transient` on the uniform source with Robin on the outer
    faces, the merged arm's route: (B), balance, (E′) ≤ their bands, the slope far
    inside ``S2_BAND``.  Control: q(1 + 1e-3) misses (B) by 1e-3, 1e7 × ``B_BAND``."""
    comm = MPI.COMM_WORLD
    fx = _k9_fixture(comm)
    msh, cell_tags = fx[0], fx[1]
    root = comm.bcast(str(tmp_path), root=0)
    this = sys.modules[_heat_columns.__module__]
    tips = {"a": fx[4]["a"]}
    with monkeypatch.context() as m:
        e_sup, sigma_f, _ = _columns_e_sup(m, fx, root)
        merged = _merge_rod_into_phantom(msh, cell_tags)
        arm = _thermal_arm(
            "K10", msh, merged, e_sup, sigma_f, None, STUB_CENTRE, tips, STUB_PROBES, merged=True
        )
        real_source = this._heat_source

        def _perturbed(*a):
            q = real_source(*a)
            q.x.array[:] *= 1.0 + PERTURBED_SOURCE
            return q

        m.setattr(this, "_heat_source", _perturbed)
        neg = _thermal_arm(
            "K10 NEG",
            msh,
            merged,
            e_sup,
            sigma_f,
            None,
            STUB_CENTRE,
            tips,
            STUB_PROBES,
            merged=True,
        )
    fails, neg_fails = _thermal_failures(arm), _thermal_failures(neg)
    if comm.rank == 0:
        print(
            f"\n[GEO-36 I32a K10-E] (B) {arm['B']:.3e}, (E) {arm['E']:.3e}, (E') {arm['Ep']:.3e}, "
            f"dt {arm['dt']:.6f} s; (K10-S) slope {arm['S2']:.4e}, tip {arm['tip_s2']['a']:.4e} (band "
            f"{S2_BAND:.0e}), failures {fails or 'none'} (ASSERTED)"
        )
        print(
            f"[GEO-36 I32a K10-E] control q(1 + 1e-3): (B) {neg['B']:.9e} = "
            f"{neg['B'] / B_BAND:.4e} x B_BAND (named ratio {PERTURBED_B_RATIO:.4e}, ASSERTED); failures "
            f"{neg_fails}",
            flush=True,
        )
    assert not fails, fails
    assert arm["dt"] == T_STAR_BOX / N_MAIN, arm["dt"]
    assert arm["counts"]["wall"] == 0 and arm["counts"]["outer"] == arm["counts"]["exterior"]
    assert arm["S2"] <= 0.5 * S2_BAND, arm["S2"]
    assert all(arm["probe_ok"]), arm["probe_ok"]
    assert abs(neg["B"] / PERTURBED_SOURCE - 1.0) <= 1.0e-9, neg["B"]
    assert neg["B"] / B_BAND >= 0.999 * PERTURBED_B_RATIO, neg["B"]
    assert any("(B)" in f for f in neg_fails), neg_fails


@complex_only
@_coding
def test_i32a_k10c_constant_curl_closed_form():
    """(K10-C): E₀(−y, x, 0)/L, a constant-curl field N1curl-1 holds exactly, on
    (1 + j): (B) against ∫½σ|E|² = σE₀²(16/3)a⁵/L² over the cube of half-width a."""
    comm = MPI.COMM_WORLD
    msh, cell_tags, e_u, sigma_f, _tips, _t1 = _k9_fixture(comm)
    a, length = 0.05, 0.05
    e_sup = fem.Function(e_u.function_space)
    e_sup.interpolate(
        lambda x: (1.0 + 1.0j) * STUB_E0 / length * np.vstack([-x[1], x[0], 0.0 * x[0]])
    )
    e_sup.x.scatter_forward()
    _merged, sub, _ids, q = _merged_q(msh, cell_tags, e_sup, sigma_f)
    int_q = float(np.real(_allsum(comm, q * ufl.dx(domain=sub))))
    exact = STUB_SIGMA * STUB_E0**2 * (16.0 / 3.0) * a**5 / length**2
    rel = abs(int_q / exact - 1.0)
    if comm.rank == 0:
        print(
            f"\n[GEO-36 I32a K10-C] int q {int_q:.15e} vs closed form {exact:.15e}: rel {rel:.3e} (band "
            f"{B_BAND:.0e}, ASSERTED)",
            flush=True,
        )
    assert rel <= B_BAND, rel


@_coding
def test_i32a_k10p_fixture_map():
    """(K10-p): P0 … P5 and three planted points, exactly; a point off the box refuses."""
    off = 0.0123
    pts = _p37_fixture(off)
    expected = {
        "P0 centre": (0.0, off, 0.0),
        "P1 +y wall": (B_BOX, off, 0.0),
        "P2 -y wall": (-B_BOX, off, 0.0),
        "P3 +x wall": (0.0, off, A_BOX),
        "P4 -x wall": (0.0, off, -A_BOX),
        "P5 top face": (0.0, off + C_BOX, 0.0),
    }
    assert pts == expected, pts
    assert _map_box_point((A_BOX, B_BOX, C_BOX), off) == (B_BOX, off + C_BOX, A_BOX)
    assert _map_box_point((-A_BOX, 0.0, C_BOX), off) == (0.0, off + C_BOX, -A_BOX)
    assert _map_box_point((0.1, 0.2, 0.03), off) == (0.2, off + 0.03, 0.1)
    with pytest.raises(AssertionError):
        _map_box_point((0.0, 0.0, 2.0 * C_BOX), off)
    print(f"\n[GEO-36 I32a K10-p] mapped {pts}", flush=True)


@complex_only
@_coding
def test_i32a_k10m_merged_census():
    """(K10-M): on the small mesh the merged tag 3 is the rod arm's 3 + 403, and a
    census one cell short fails it."""
    from tests.validation import test_geo38_f_clinical_column as col

    comm = MPI.COMM_WORLD
    msh, cell_tags = _k9_fixture(comm)[:2]
    census = col._owned_tag_counts(msh, cell_tags, comm)
    merged = _merge_rod_into_phantom(msh, cell_tags)
    gel, rod = _owned_tag_count(msh, merged, PHANTOM_CELL_TAG), _owned_tag_count(msh, merged, 403)
    if comm.rank == 0:
        print(
            f"\n[GEO-36 I32a K10-M] census {census}, merged tag 3 {gel}, tag 403 {rod}", flush=True
        )
    assert census[403] > 0
    assert gel == msh.topology.index_map(3).size_global
    assert _s0_failures(census, dict(census), gel, rod) == []
    assert _s0_failures(census, dict(census), gel - 1, rod)
    assert _s0_failures(census, dict(census), gel, 1)


@complex_only
@_coding
def test_i32a_k10r_run_heat_end_to_end_on_stubs(monkeypatch, tmp_path):
    """`_run_heat` end to end, fixture and `_drive` stubbed (items 701 / 702's
    pattern).  One G on both arms, and with the twin's field half the rod's the
    deposited powers read G·P(rod) = 4 V₃ / V, the mesh's own volumes."""
    from tests.validation.test_geo36_rod_f_clinical_drive import K9_NODE_PF

    comm = MPI.COMM_WORLD
    fx = _k9_fixture(comm)
    msh, cell_tags = fx[0], fx[1]
    root = comm.bcast(str(tmp_path), root=0)
    with monkeypatch.context() as m:
        lines = _stubs(m, fx, root)
        out = _run_heat(comm, m, K9_NODE_PF)
    rod, twin = out["rod"], out["twin"]
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    one = fem.Constant(msh, default_scalar_type(1.0))
    v3 = float(np.real(_allsum(comm, one * dx(PHANTOM_CELL_TAG))))
    v_all = float(np.real(_allsum(comm, one * dx)))
    dep_rod = rod["G"] * rod["P"]
    r_dep = abs(dep_rod / (4.0 * v3 / v_all) - 1.0)
    heads = ("(M2)", "(S0)", "rod share", "tip theta(t*) ratio", "price:", "gates: failures none")
    missing = [h for h in heads if not any(ln.startswith(h) for ln in lines)]
    sense_free = [ln for ln in lines if "theta(t*)" in ln and SENSE not in ln]
    if comm.rank == 0:
        print(
            f"\n[GEO-36 I32a K10-R] G rod {rod['G']:.12e} = twin {twin['G']:.12e}; G*P twin "
            f"{twin['G'] * twin['P']:.12f}, rod {dep_rod:.12f} vs 4 V3/V: rel {r_dep:.3e}; walls rod "
            f"{rod['counts']['wall']}, twin {twin['counts']['wall']}; ratios {out['ratios']}; probes "
            f"valid {list(rod['probe_ok'])} / {list(twin['probe_ok'])}; heads missing {missing or 'none'}",
            flush=True,
        )
    assert rod["G"] == twin["G"]
    assert abs(twin["G"] * twin["P"] - 1.0) <= B_BAND
    assert r_dep <= B_BAND, r_dep
    assert rod["counts"]["wall"] > 0 and twin["counts"]["wall"] == 0
    assert all(rod["probe_ok"]) and all(twin["probe_ok"])
    assert set(out["ratios"]) == set(fx[4])
    assert not missing, missing
    assert not sense_free, sense_free


# ---------------------------------------------------------------------------
# I32b (§9 item 870): the (S2-tip) comparand re-posed on Πq, stub-proved.
# ---------------------------------------------------------------------------
#: `OPS-83`'s printed R_Π at Δt₁/10 and |⟨Πq⟩/⟨q⟩ − 1| (`20261106T120232Z_OPS-83-W1.log:69`, `:75–76`).
OPS83_R_PI_R = {"a": 6.8552e-07, "interior": 8.5219e-07}
OPS83_PQ_Q_R = {"a": 1.573e-03, "interior": 2.115e-03}
OPS83_R_PI_U_TIP = 8.0980e-07
REP_TOL_4E = 5.0e-12  # half a unit of the `.4e` print's last digit (reproduction control)
REP_TOL_3E = 5.0e-07  # half a unit of the `.3e` print's last digit
#: The rank width of `OPS-83`'s record (`OPS-83-W1.log:12`, ``mpiexec -n 2``).  The
#: literal R_Π reproductions (K11-rep's R_Π, K11-u's 8.0980e-07) are asserted at this
#: width only.  Measured (I32b W0, `20261108T103328Z_GEO-36-I32b-W0.log:141`): at
#: ``-n 1`` they read |diff| 8.598e-10 (tip) / 5.407e-11 (interior) / 2.641e-10 (u tip),
#: while at ``-n 2`` (W0b, `20261108T103407Z_GEO-36-I32b-W0b.log:161`) 4.402e-12 /
#: 1.436e-12 / 3.954e-12.  The solve carries T absolute (T∞ ≈ 293 K against a θ¹ of a
#: few µK at Δt₁/10), so R at 1e-9 sits on the partition-dependent round-off floor the
#: (E) residual also shows (≈ 1e-10 … 5e-10); K10-S's own tip reads 8.0402e-06 at
#: ``-n 1`` against 8.0400e-06 at ``-n 2`` (`20261101T103903Z_GEO-36-I32a-W0.log:63`).
#: 5e-12 is a fixed-partition reproduction, not a physics bound; at other widths the
#: values are printed, and K11-Π, the |pq_q| reproduction, the K11-u identity and
#: K11-c stay asserted.
#: Width-scoped by the 2026-11-08 10:30 review ruling (`GEO-36.md`, provisional; the 11-09 weekly ratifies).
OPS83_RECORD_RANKS = 2


@complex_only
@_coding
def test_i32b_k11_projected_tip_gate_on_checkerboard(monkeypatch, tmp_path):
    """(K11): `_thermal_arm(..., tip_projection=True)` on K10es's set-up, on q_u and on
    `OPS-83`'s checkerboard q_r, reproduces `OPS-83`'s R_Π; the retired DG0 form misses
    on q_r (K11-c, the negative control)."""
    from tests.validation.test_ops83_s2tip_form_on_stub import PROJECTION_ID_TOL, _checkerboard

    comm = MPI.COMM_WORLD
    fx = _k9_fixture(comm)
    msh, cell_tags = fx[0], fx[1]
    root = comm.bcast(str(tmp_path), root=0)
    this = sys.modules[_heat_columns.__module__]
    tips = {"a": fx[4]["a"], "interior": STUB_PROBES["interior"]}
    with monkeypatch.context() as m:
        e_sup, sigma_f, _ = _columns_e_sup(m, fx, root)
        merged = _merge_rod_into_phantom(msh, cell_tags)
        arm_u = _thermal_arm(
            "K11 u",
            msh,
            merged,
            e_sup,
            sigma_f,
            None,
            STUB_CENTRE,
            tips,
            STUB_PROBES,
            merged=True,
            tip_projection=True,
        )
        real_source = this._heat_source

        def _rough(*a):
            return _checkerboard(real_source(*a))[0]

        m.setattr(this, "_heat_source", _rough)
        arm_r = _thermal_arm(
            "K11 r",
            msh,
            merged,
            e_sup,
            sigma_f,
            None,
            STUB_CENTRE,
            tips,
            STUB_PROBES,
            merged=True,
            tip_projection=True,
        )
    fails = []
    for arm in (arm_u, arm_r):
        if not arm["E"] <= E_BAND:
            fails.append(f"{arm['label']} (E) {arm['E']:.3e}")
        if not arm["S2"] <= S2_BAND:
            fails.append(f"{arm['label']} (S2) {arm['S2']:.3e}")
    if not arm_u["B"] <= B_BAND:
        fails.append(f"K11 u (B) {arm_u['B']:.3e}")
    ru, rr = arm_u["tip_rows"]["dt1/10"], arm_r["tip_rows"]["dt1/10"]
    literal = comm.size == OPS83_RECORD_RANKS
    for n in tips:
        if not arm_r["tip_s2_pi"][n] <= S2_BAND:
            fails.append(f"(K11-Pi) {n} {arm_r['tip_s2_pi'][n]:.4e}")
        if literal and not abs(rr[n]["R_Pi"] - OPS83_R_PI_R[n]) <= REP_TOL_4E:
            fails.append(f"(K11-rep) R_Pi {n} {rr[n]['R_Pi']:.6e} vs {OPS83_R_PI_R[n]:.4e}")
        if not abs(abs(rr[n]["pq_q"]) - OPS83_PQ_Q_R[n]) <= REP_TOL_3E:
            fails.append(f"(K11-rep) pq_q {n} {rr[n]['pq_q']:.6e} vs {OPS83_PQ_Q_R[n]:.3e}")
        if not abs(ru[n]["R_Pi"] - ru[n]["R_DG0"]) <= PROJECTION_ID_TOL:
            fails.append(f"(K11-u) R_Pi - R_DG0 {n} {ru[n]['R_Pi'] - ru[n]['R_DG0']:.3e}")
    if literal and not abs(ru["a"]["R_Pi"] - OPS83_R_PI_U_TIP) <= REP_TOL_4E:
        fails.append(f"(K11-u) R_Pi tip {ru['a']['R_Pi']:.6e} vs {OPS83_R_PI_U_TIP:.4e}")
    neg = arm_r["tip_s2"]["a"]
    if comm.rank == 0:
        for arm in (arm_u, arm_r):
            print(
                f"\n[GEO-36 I32b K11] [{arm['label']}] (B) {arm['B']:.3e} (band {B_BAND:.0e}, "
                f"{'ASSERTED' if arm is arm_u else 'printed'}); (E) {arm['E']:.3e} (E dt1/10 "
                f"{arm['E10']:.3e}; band {E_BAND:.0e}, ASSERTED); global (S2) {arm['S2']:.4e} (band "
                f"{S2_BAND:.0e}, ASSERTED)",
                flush=True,
            )
            for n in tips:
                a, b = arm["tip_rows"]["dt1"][n], arm["tip_rows"]["dt1/10"][n]
                print(
                    f"[GEO-36 I32b K11] [{arm['label']}] ball {n}: R_Pi dt1 {a['R_Pi']:.6e} dt1/10 "
                    f"{b['R_Pi']:.6e}; R_DG0 dt1 {a['R_DG0']:.6e} (retired form) dt1/10 "
                    f"{b['R_DG0']:.6e}; <Pq>/<q> - 1 {b['pq_q']:+.6e}",
                    flush=True,
                )
        rep = ", ".join(f"{n} |diff| {abs(rr[n]['R_Pi'] - OPS83_R_PI_R[n]):.3e}" for n in tips)
        rep_pq = ", ".join(f"{n} {abs(abs(rr[n]['pq_q']) - OPS83_PQ_Q_R[n]):.3e}" for n in tips)
        lit_word = (
            "ASSERTED"
            if literal
            else "printed at -n %d (record is -n %d)" % (comm.size, OPS83_RECORD_RANKS)
        )
        ident = ", ".join(f"{n} {abs(ru[n]['R_Pi'] - ru[n]['R_DG0']):.3e}" for n in tips)
        print(
            f"[GEO-36 I32b K11] (K11-rep) vs OPS-83 R_Pi(r, dt1/10): {rep} (tol {REP_TOL_4E:.0e}); "
            f"|pq_q| |diff| {rep_pq} (tol {REP_TOL_3E:.0e}); (K11-u) tip R_Pi |diff| "
            f"{abs(ru['a']['R_Pi'] - OPS83_R_PI_U_TIP):.3e}, R_Pi - R_DG0 {ident} (tol "
            f"{PROJECTION_ID_TOL:.0e}); literal R_Pi reproductions {lit_word}, |pq_q| and identity "
            f"ASSERTED",
            flush=True,
        )
        print(
            f"[GEO-36 I32b K11] (K11-c) control: retired DG0 tip_s2[a] on q_r {neg:.4e} = "
            f"{neg / S2_BAND:.2f} x S2_BAND (> band, ASSERTED; OPS-83 read 1.5793e-03); failures "
            f"{fails or 'none'}",
            flush=True,
        )
    assert not fails, fails
    assert neg > S2_BAND, neg


@_coding
def test_i32b_k11f_thermal_failures_gate_projected_tip():
    """(K11-f): on a non-merged arm `_thermal_failures` gates ``tip_s2_pi`` when present."""
    base = {
        "label": "rod",
        "B": 0.0,
        "E": 0.0,
        "Ep": 0.0,
        "S2": 0.0,
        "merged": False,
        "counts": {"wall": 4},
        "tip_s2": {"a": 1.6e-03},
    }
    green = _thermal_failures(dict(base, tip_s2_pi={"a": 7.0e-07}))
    red = _thermal_failures(dict(base, tip_s2_pi={"a": 2.0e-04}))
    retired = _thermal_failures(dict(base))
    print(
        f"\n[GEO-36 I32b K11-f] green {green or 'none'}; red {red}; no tip_s2_pi key (old form) "
        f"{retired}",
        flush=True,
    )
    assert not any("(S2-tip)" in f for f in green), green
    assert len(red) == 1 and "(S2-tip)" in red[0], red
    assert len(retired) == 1 and "(S2-tip)" in retired[0], retired


# ---------------------------------------------------------------------------
# I87' (§9 item 968): the ROD_THERMAL knob, stub-proved.
# ---------------------------------------------------------------------------
#: L-W0b's literals (`20261108T181211Z_GEO-36-I32b-L-W0b.log:104`, `:106`, `:117`).
I87P_N_ROD = 108
I87P_GP_STR = "3.958333333333e+00"
I87P_VSUB_STR = "9.999999999999e-04"
#: (N-src)'s ceiling, the exact separation V403 / V3 on the stub (`L-W0b.log:103–104`); predicted.
I87P_NSRC_PRED = 108.0 / 10260.0
#: (N-ρc) predicted factors (item 933's unlogged arithmetic); printed beside, never a band.
I87P_NRC_PRED = {"ti": 3.3e-3, "hole": 7.4e-3}
#: The weekly's predicted ranges for the former gates (`PROJECT_PLAN.md:4935–4938`), printed only.
I87P_FORMER = (
    "(K-gel) |gel / default - 1| 4e-3 ... 7e-3",
    "(K-Ti) Ti / gel - 1 +1.5e-3 ... +3.5e-3",
    "(K-hole) hole / gel - 1 +4e-3 ... +7e-3",
    "hole / default - 1 <= 1e-4",
)


def _line_token(lines, head, key):
    """The whitespace token after ``key`` on the first recorded line starting with ``head``."""
    for ln in lines:
        s = ln.lstrip("\n")
        if s.startswith(head) and key in s:
            return s.split(key, 1)[1].split()[0].rstrip(";,")
    return None


@complex_only
@_coding
def test_i87p_k0_knob_parsing(monkeypatch):
    monkeypatch.delenv(ROD_THERMAL_ENV, raising=False)
    assert _rod_thermal() is None
    monkeypatch.setenv(ROD_THERMAL_ENV, "")
    assert _rod_thermal() is None
    for arm in ROD_THERMAL_ARMS:
        monkeypatch.setenv(ROD_THERMAL_ENV, arm)
        assert _rod_thermal() == arm
    refused = []
    for bad in ("TI", "1", "x", "Gel", "rod"):
        monkeypatch.setenv(ROD_THERMAL_ENV, bad)
        with pytest.raises(ValueError, match=ROD_THERMAL_ENV):
            _rod_thermal()
        refused.append(bad)
    print(
        f"\n[GEO-36 I87p K0] unset / '' -> None; arms {ROD_THERMAL_ARMS} accepted; refused {refused} "
        f"(ValueError naming {ROD_THERMAL_ENV})",
        flush=True,
    )


@complex_only
@_coding
def test_i87p_rod_thermal_arms_on_stubs(monkeypatch, tmp_path):
    """(K-def)(b), (K-arm), (K-gel=), (N-src), (N-ρc), (K-c) on the stub; the three
    former (K-…) gates printed with the weekly's predicted ranges, never asserted."""
    from tests.validation.test_geo36_rod_f_clinical_3t_heating import T_REP, _find
    from tests.validation.test_geo36_rod_f_clinical_drive import K9_NODE_PF

    comm = MPI.COMM_WORLD
    fx = _k9_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    literal = comm.size == OPS83_RECORD_RANKS
    runs = {}
    with monkeypatch.context() as m:
        m.delenv(ROD_THERMAL_ENV, raising=False)
        lines = _stubs(m, fx, root)
        runs[None] = (lines, _run_heat(comm, m, K9_NODE_PF))
    for arm in ROD_THERMAL_ARMS:
        with monkeypatch.context() as m:
            lines = _stubs(m, fx, root)
            m.setenv(ROD_THERMAL_ENV, arm)
            runs[arm] = (lines, _run_heat(comm, m, K9_NODE_PF))

    fails, notes = [], []
    # (K-def)(b): the unset run reproduces T_REP and carries no knob trace.
    lines0, out0 = runs[None]
    missing = {h: _find(lines0, h, toks) for h, toks in T_REP}
    missing = {h: t for h, t in missing.items() if t}
    if missing:
        (fails if literal else notes).append(f"(K-def)(b) T_REP missing {missing}")
    rt_lines = [ln for ln in lines0 if ln.lstrip("\n").startswith("rod thermal")]
    if rt_lines:
        fails.append(f"(K-def)(b) unset run printed {rt_lines}")
    if "rod_cells" in out0:
        fails.append("(K-def)(b) rod_cells in the unset run's out")
    tips = sorted(out0["rod"]["tip_th"])
    default = {
        "t*": {n: out0["rod"]["tip_th"][n] for n in tips},
        "F2182": {n: out0["rod"]["f2182_th"][n] for n in tips},
    }
    for arm in ROD_THERMAL_ARMS:
        exc = runs[arm][1]["rod"]
        for kind, key in (("t*", "tip_th"), ("F2182", "f2182_th")):
            for n in tips:
                s0, s1 = f"{default[kind][n]:.9e}", f"{exc[key][n]:.9e}"
                if s0 != s1:
                    (fails if literal else notes).append(
                        f"(K-def)(b) {arm} run's exclusion {kind} {n} {s1} vs unset {s0}"
                    )

    # (K-arm), (K-gel=), (N-src), (N-rho_c), per arm.
    per = {}
    for arm in ROD_THERMAL_ARMS:
        lines, out = runs[arm]
        rc, exc = out.get("rod_cells"), out["rod"]
        if rc is None:
            fails.append(f"(K-arm) {arm}: no rod_cells in out")
            continue
        gp = rc["G"] * rc["P"]
        gp_exc = exc["G"] * exc["P"]
        vsub_twin = _line_token(lines, "[gel twin] sense", "V_sub ")
        vsub_rc = f"{rc['v_sub']:.12e}"
        per[arm] = rc
        if gp != gp_exc:
            fails.append(f"(K-arm) {arm} G*P {gp!r} != exclusion {gp_exc!r}")
        if literal and f"{gp:.12e}" != I87P_GP_STR:
            fails.append(f"(K-arm) {arm} G*P {gp:.12e} vs {I87P_GP_STR}")
        if rc["counts"]["wall"] != 0:
            fails.append(f"(K-arm) {arm} wall {rc['counts']['wall']}")
        if vsub_rc != vsub_twin:
            fails.append(f"(K-arm) {arm} V_sub {vsub_rc} vs twin {vsub_twin}")
        if literal and vsub_rc != I87P_VSUB_STR:
            fails.append(f"(K-arm) {arm} V_sub {vsub_rc} vs {I87P_VSUB_STR}")
        if rc["n_rod"] != I87P_N_ROD:
            fails.append(f"(K-arm) {arm} n_rod {rc['n_rod']}")
        for n in tips:
            if rc["v_rod_ball"][n] != 0.0:
                fails.append(f"(K-arm) {arm} V(B cap rod) {n} {rc['v_rod_ball'][n]!r}")
            d_tip = abs(rc["tip_s2_pi"][n] - rc["tip_s2_pi_sc"][n])
            if not d_tip <= S2_BAND:
                fails.append(f"(K-gel=) {arm} tip {n} rho_c vs scalar form {d_tip:.3e}")
        if not rc["n_src"] > B_BAND:
            fails.append(f"(N-src) {arm} {rc['n_src']:.3e}")
        if arm == "gel":
            if not abs(rc["Ep_rc"] - rc["Ep_sc"]) <= E_BAND:
                fails.append(f"(K-gel=) Ep_rc - Ep_sc {rc['Ep_rc'] - rc['Ep_sc']:.3e}")
            if not abs(rc["S2_rc"] - rc["S2_sc"]) <= S2_BAND:
                fails.append(f"(K-gel=) S2_rc - S2_sc {rc['S2_rc'] - rc['S2_sc']:.3e}")
        else:
            if not rc["Ep_sc"] > E_BAND:
                fails.append(f"(N-rho_c) {arm} Ep_sc {rc['Ep_sc']:.3e}")
        if comm.rank == 0:
            print(
                f"\n[GEO-36 I87p] [{arm}] (K-arm) G*P {gp:.12e} (exclusion {gp_exc:.12e}, float equal "
                f"{gp == gp_exc}); V_sub {vsub_rc} (twin {vsub_twin}); wall {rc['counts']['wall']}; n_rod "
                f"{rc['n_rod']}; V(B cap rod) "
                + ", ".join(f"{n} {rc['v_rod_ball'][n]!r}" for n in tips),
                flush=True,
            )
            print(
                f"[GEO-36 I87p] [{arm}] (B) {rc['B']:.3e}, (E) {rc['E']:.3e}, (E')rho_c {rc['Ep_rc']:.3e}, "
                f"(E') scalar {rc['Ep_sc']:.6e}, (S2)rho_c {rc['S2_rc']:.6e}, (S2) scalar "
                f"{rc['S2_sc']:.6e}; (S2-tip) rho_c form "
                + ", ".join(f"{n} {rc['tip_s2_pi'][n]:.4e}" for n in tips)
                + "; scalar form "
                + ", ".join(f"{n} {rc['tip_s2_pi_sc'][n]:.4e}" for n in tips),
                flush=True,
            )
            print(
                f"[GEO-36 I87p] [{arm}] (N-src) measured {rc['n_src']:.6e} vs predicted {I87P_NSRC_PRED:.6e} "
                f"(> B_BAND {B_BAND:.0e}, ASSERTED; predicted, the form is new)"
                + (
                    f"; (N-rho_c) |Ep_sc| {rc['Ep_sc']:.6e} vs predicted ~{I87P_NRC_PRED[arm]:.1e} "
                    f"(> E_BAND {E_BAND:.0e}, ASSERTED)"
                    if arm in I87P_NRC_PRED
                    else f"; (K-gel=) |Ep_rc - Ep_sc| {abs(rc['Ep_rc'] - rc['Ep_sc']):.3e} (<= E_BAND), "
                    f"|S2_rc - S2_sc| {abs(rc['S2_rc'] - rc['S2_sc']):.3e} (<= S2_BAND) ASSERTED"
                ),
                flush=True,
            )

    # (K-c): one θ digit of the `[rod] 1 g ball … at a,` token, in memory.
    head_c = "[rod] 1 g ball (a 6.2035 mm) at a,"
    tok_c = dict(T_REP)[head_c][0]
    bad_tok = tok_c.replace("2.424561997e+00", "2.424561998e+00")
    kc = _find(lines0, head_c, (bad_tok,))
    if kc != [bad_tok] or bad_tok == tok_c:
        fails.append(f"(K-c) _find did not miss {bad_tok!r}: {kc}")

    # The three former gates: printed with the weekly's ranges, never asserted.
    if comm.rank == 0 and set(per) == set(ROD_THERMAL_ARMS):
        print(
            "[GEO-36 I87p] former gates (printed, never asserted): " + "; ".join(I87P_FORMER),
            flush=True,
        )
        for kind, key in (("theta(t*)", "tip_th"), ("F2182 theta(900 s)", "f2182_th")):
            dk = "t*" if key == "tip_th" else "F2182"
            for n in tips:
                d = default[dk][n]
                g, ti, ho = (per[a][key][n] for a in ROD_THERMAL_ARMS)
                print(
                    f"[GEO-36 I87p] former gates at {n}, {kind}: default {d:.9e} gel {g:.9e} ti {ti:.9e} "
                    f"hole {ho:.9e} K; (K-gel) |gel/default - 1| {abs(g / d - 1.0):.4e}; (K-Ti) ti/gel - 1 "
                    f"{ti / g - 1.0:+.4e}; (K-hole) hole/gel - 1 {ho / g - 1.0:+.4e}; hole/default - 1 "
                    f"{ho / d - 1.0:+.4e} (printed)",
                    flush=True,
                )
                vals = {"default": d, "gel": g, "ti": ti, "hole": ho}
                order = sorted(vals, key=lambda a: vals[a])
                print(
                    f"[GEO-36 I87p] ordering at {n} ({kind}): " + " < ".join(order) + " (printed)",
                    flush=True,
                )
        print(
            f"[GEO-36 I87p] (K-def)(b) unset run: T_REP missing {missing or 'none'} (asserted at -n "
            f"{OPS83_RECORD_RANKS}, -n {comm.size} here); rod thermal lines {len(rt_lines)}; rod_cells in "
            f"out {'rod_cells' in out0}; exclusion arms agree by string with the unset [rod]: "
            f"{not any('exclusion' in f for f in fails + notes)}; (K-c) {bad_tok!r} reported missing "
            f"{kc == [bad_tok]}; notes {notes or 'none'}; failures {fails or 'none'}",
            flush=True,
        )
    assert not fails, fails
