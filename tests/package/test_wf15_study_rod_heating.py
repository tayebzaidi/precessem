"""`WF-15` step 14 (§9 item 831, §10 U127): the rod's heating through the ``Study``.

``Study.temperature(cooling, transient=..., balls=..., drive="P1")`` on U97′'s
PEC implant-rod four-port (``FourPortBirdcageGeometry(implant_rod=ROD)``,
degree 1, 64 MHz, exact-rectangle lumped ports, κ off), on `IMP-3` I41's
route under the 2026-11-05 03:00 review ruling: complex build, one port
(``P1`` at 1 V) — the four-port ``Study`` has no ccw superposition — arm R
through the ``Study``, arm M (the kept rod with tag 403 as gel) built test-side
and heated with the same private core.  Every band is imported.

Knob ``FEM_EM_WF15_S14`` = ``1``; unset, the coding tests run (no mesh) and the
solve node skips.  Not claimed: a ccw drive, a converged ratio, a tissue model,
compliance, the copper rod, or a sweep.
"""

from __future__ import annotations

import gc
import os
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from mpi4py import MPI

import fem_em_solver.study as study_module
from fem_em_solver import (
    CapacitorSheets,
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    RobinCooling,
    Study,
    TransientSpec,
)
from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
from fem_em_solver.post.sar import mean_sar
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import ROD, _build_kept_rod
from tests.mesh.test_birdcage_leg_offset import IDENTITY_BAND
from tests.package.test_wf15_study_implant_rod import ROD_EXACT_CELL_RECORD, _phantom
from tests.validation.test_birdcage_power_identity import STEP2D_FREQUENCY_HZ
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
from tests.validation.test_geo36_rod_f_clinical_heating import (
    _projected_ball_comparands,
    _tip_slope_rows,
)
from tests.validation.test_imp3_rod_birdcage import TIPS
from tests.validation.test_imp3_rod_heating_matched import (
    _arm_with_403,
    _merge_rod_into_phantom,
    _owned_tag_count,
)
from tests.validation.test_imp3_rod_heating_phantom import _parent_wall_count
from tests.validation.test_mass_averaged_sar_standard_masses import QUADRATURE_DEGREE
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    DT_1,
    E_BAND,
    H_GEL,
    K_GEL,
    N_MAIN,
    RHO_C,
    S2_BAND,
    T_INF,
    T_STAR,
)
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION

KNOB = "FEM_EM_WF15_S14"
DRIVE = "P1"
MASS_G = 1.0
#: I41's tip θ̄(t*) R / M (`20261025T201613Z_IMP-3-I41-W1b.log:3776–3777`), ccw, ragged mesh.
I41_TIP_RATIO = {"top": 1.104515, "bottom": 1.108684}

needs_knob = pytest.mark.skipif(os.environ.get(KNOB) != "1", reason=f"{KNOB} unset")


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _study(geometry, excitation=None, degree=1):
    return Study(
        geometry,
        _phantom(),
        excitation or LumpedPortDrive(),
        frequency_hz=STEP2D_FREQUENCY_HZ,
        degree=degree,
    )


def _cooling(scale=1.0):
    return RobinCooling(k_w_per_m_k=K_GEL, h_w_per_m2_k=H_GEL, t_inf_c=T_INF, source_scale=scale)


def _cn():
    return TransientSpec(
        dt_s=T_STAR / N_MAIN,
        n_steps=N_MAIN,
        rho_c_j_per_m3_k=RHO_C,
        theta=0.5,
        startup_half_steps=2,
    )


def _be(dt=DT_1):
    return TransientSpec(
        dt_s=dt, n_steps=1, rho_c_j_per_m3_k=RHO_C, theta=1.0, startup_half_steps=0
    )


def _balls():
    return {"mass_g": MASS_G, "centres": [list(c) for c in TIPS.values()]}


def _no_mesh(monkeypatch):
    def _boom(*a, **k):
        raise AssertionError("a mesh was built")

    monkeypatch.setattr(FourPortBirdcageGeometry, "build", _boom)


# ---------------------------------------------------------------------------
# Coding tests (W0): refusals before any mesh.
# ---------------------------------------------------------------------------


@complex_only
def test_s14_refusals_before_any_mesh(monkeypatch):
    _no_mesh(monkeypatch)
    rod = FourPortBirdcageGeometry(implant_rod=ROD)
    kw = {"transient": _cn(), "balls": _balls(), "drive": DRIVE}
    with pytest.raises(ValueError, match="degree 2"):
        _study(rod, degree=2).temperature(_cooling(), **kw)
    with pytest.raises(TypeError, match="CapacitorSheets"):
        _study(rod, excitation=CapacitorSheets(c_f=1.0e-11, esr_ohm=0.2, on=("P2",))).temperature(
            _cooling(), **kw
        )
    # (c-free) the rod-free four-port keeps its TypeError (transient omitted so the
    # unchanged port refusal is the one reached).
    with pytest.raises(TypeError, match="port Study"):
        _study(FourPortBirdcageGeometry()).temperature(_cooling())
    # drive off the rod is refused.
    with pytest.raises(ValueError, match="drive is defined"):
        _study(FourPortBirdcageGeometry()).temperature(_cooling(), drive=DRIVE)
    st = _study(rod)
    with pytest.raises(ValueError, match="drive must name"):
        st.temperature(_cooling(), transient=_cn(), balls=_balls())
    with pytest.raises(ValueError, match="initial_gradient"):
        st.temperature(
            _cooling(),
            transient=TransientSpec(
                dt_s=1.0, n_steps=1, rho_c_j_per_m3_k=RHO_C, initial_gradient_k_per_m=1.0
            ),
            balls=_balls(),
            drive=DRIVE,
        )
    with pytest.raises(TypeError, match="TransientSpec"):
        st.temperature(_cooling(), balls=_balls(), drive=DRIVE)
    with pytest.raises(ValueError, match="balls"):
        st.temperature(_cooling(), transient=_cn(), balls={"mass_g": 1.0}, drive=DRIVE)
    with pytest.raises(ValueError, match="sense"):
        st.temperature(_cooling(), transient=_cn(), balls=_balls(), drive=DRIVE, sense="ccw")
    _p0(
        "\n[WF-15 s14] refusals: degree 2, caps, rod-free, drive off the rod, gradient, "
        "transient, balls, sense -- all before any mesh"
    )


YAML = Path(__file__).parent / "wf15" / "rod_heating.yaml"
ROD_FREE_THERMAL_MSG = "temperature is not defined on a port Study (four_port_birdcage)"


def _yaml_mapping():
    import yaml

    with open(YAML) as fh:
        return yaml.safe_load(fh)


@complex_only
def test_s14_yaml_parse_and_refusals(monkeypatch):
    """§9 item 884 (U127′): the rod four-port ``thermal:`` parses to the test's own
    Study / cooling / transient / balls / drive; the rod-free four-port's message is
    unchanged; ``thermal:`` beside ``study.sweep_hz`` is refused.  No mesh."""
    from fem_em_solver import study_yaml as sy

    _no_mesh(monkeypatch)
    st, kind, _outputs, th = sy._parse_mapping(_yaml_mapping(), str(YAML))
    ref = _study(FourPortBirdcageGeometry(implant_rod=ROD))
    assert kind == "four_port_birdcage"
    assert isinstance(th, sy._RodThermal)
    assert st.geometry == ref.geometry and dict(st.geometry.implant_rod) == dict(ROD)
    assert st.degree == 1 and st.frequency_hz == STEP2D_FREQUENCY_HZ
    assert st.material == ref.material
    assert th.drive == DRIVE
    c = th.cooling
    assert (c.k_w_per_m_k, c.h_w_per_m2_k, c.t_inf_c, c.degree, c.source_scale) == (
        K_GEL,
        H_GEL,
        T_INF,
        1,
        1.0,
    )
    cn = _cn()
    t = th.transient
    assert abs(t.dt_s / cn.dt_s - 1.0) <= IDENTITY_BAND, (t.dt_s, cn.dt_s)
    assert (t.n_steps, t.rho_c_j_per_m3_k, t.theta, t.startup_half_steps) == (
        cn.n_steps,
        cn.rho_c_j_per_m3_k,
        cn.theta,
        cn.startup_half_steps,
    )
    assert th.balls["mass_g"] == MASS_G
    assert [tuple(x) for x in th.balls["centres"]] == [tuple(x) for x in TIPS.values()]
    # (c-free) the rod-free four-port's message, byte-identical.
    m = _yaml_mapping()
    del m["implant"]
    with pytest.raises(sy.StudyConfigError) as ei:
        sy._parse_mapping(m, "<t>")
    assert ei.value.field == "thermal" and ei.value.reason == ROD_FREE_THERMAL_MSG
    # thermal beside study.sweep_hz: refused (whichever field fires first is asserted).
    m = _yaml_mapping()
    m["study"]["sweep_hz"] = [10.0e6, 64.0e6]
    with pytest.raises(sy.StudyConfigError) as ei:
        sy._parse_mapping(m, "<t>")
    assert ei.value.field == "thermal" and "sweep_hz" in ei.value.reason, ei.value
    # drive missing / window_s (an F-human key) refused at parse.
    m = _yaml_mapping()
    del m["thermal"]["drive"]
    with pytest.raises(sy.StudyConfigError, match="thermal.drive"):
        sy._parse_mapping(m, "<t>")
    m = _yaml_mapping()
    m["thermal"]["window_s"] = 360.0
    with pytest.raises(sy.StudyConfigError) as ei:
        sy._parse_mapping(m, "<t>")
    assert ei.value.field == "thermal.window_s" and "unknown key" in ei.value.reason, ei.value
    _p0(
        f"[WF-15 s14p] YAML parse: {YAML.name} = the test's Study / cooling / CN / balls / "
        f"drive {DRIVE}; rod-free message unchanged; sweep, missing drive, window_s refused"
    )


# ---------------------------------------------------------------------------
# The solve (W1).
# ---------------------------------------------------------------------------


def _theta_tip(tr):
    names = list(TIPS)
    return {names[i]: float(b["theta_k"][-1]) for i, b in enumerate(tr.balls)}


def _run():
    comm = MPI.COMM_WORLD
    t0 = time.perf_counter()
    captured = {}
    real_build = FourPortBirdcageGeometry.build
    solves = {"n": 0}
    real_solve = Study._solve_port_fields

    def _capturing_build(self, comm_=MPI.COMM_WORLD):
        out = real_build(self, comm_)
        captured["rod"] = out
        return out

    def _counting_solve(self, comm_=MPI.COMM_WORLD):
        if getattr(self, "fields_result", None) is None:
            solves["n"] += 1
        return real_solve(self, comm_)

    out = {}
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(FourPortBirdcageGeometry, "build", _capturing_build)
        mp.setattr(Study, "_solve_port_fields", _counting_solve)  # (S0) class spy
        st = _study(FourPortBirdcageGeometry(implant_rod=ROD))
        sar = st.sar(MASS_G, centres=list(TIPS.values()), drives=[DRIVE], comm=comm)
        res = st.fields_result
        t_em = time.perf_counter() - t0
        f = res.fields[DRIVE]
        P = float(
            np.real(
                mean_sar(
                    f.e_complex,
                    sigma=f.sigma_field,
                    rho=PHANTOM_RHO_KG_PER_M3,
                    cell_tags=res.cell_tags,
                    comm=comm,
                    subdomain_ids=PHANTOM_CELL_TAG,
                )["dissipated_power_w"]
            )
        )
        G = 1.0 / P
        t1 = time.perf_counter()
        main = st.temperature(_cooling(G), transient=_cn(), balls=_balls(), drive=DRIVE, comm=comm)
        be = st.temperature(_cooling(G), transient=_be(), balls=_balls(), drive=DRIVE, comm=comm)
        be10 = st.temperature(
            _cooling(G), transient=_be(DT_1 / 10.0), balls=_balls(), drive=DRIVE, comm=comm
        )
        t_th = time.perf_counter() - t1
    out.update(
        st=st,
        res=res,
        sar=sar,
        P=P,
        G=G,
        main=main,
        be=be,
        be10=be10,
        solves=solves["n"],
        t_em=t_em,
        t_th=t_th,
    )

    # (W0) the Robin census on arm R, and the parent wall count.
    sub, _pid = study_module._tag_submesh(res.mesh, res.cell_tags, PHANTOM_CELL_TAG)
    _outer, counts = study_module._rod_outer_robin_counts(sub, dict(ROD))
    out["counts"] = counts
    out["parent_wall"] = _parent_wall_count(res.mesh, captured["rod"]["mesh_facet_tags"])
    # (c-Robin), printed: Robin on every exterior facet.
    allext = study_module._phantom_transient_core(
        st,
        res.mesh,
        res.cell_tags,
        res.comm,
        f,
        _cooling(G),
        _cn(),
        MASS_G,
        [tuple(c) for c in TIPS.values()],
        None,
        study_module._exterior_robin_tags,
        int(res.cells),
        DRIVE,
    )
    out["allext"] = allext
    out["t3_R"] = _owned_tag_count(res.mesh, res.cell_tags, PHANTOM_CELL_TAG)
    gc.collect()

    # (Y) the tracked YAML through `run_study_file`, with thermal.source_scale = G
    # injected (the one key the file cannot know): a second EM solve plus a CN.
    import yaml

    from fem_em_solver.study_yaml import run_study_file

    t_y = time.perf_counter()
    ydir = Path("logs/wf15/s14p")
    ypath = ydir / "rod_heating_G.yaml"
    if comm.rank == 0:
        ydir.mkdir(parents=True, exist_ok=True)
        m = _yaml_mapping()
        m["thermal"]["source_scale"] = float(G)
        with open(ypath, "w") as fh:
            yaml.safe_dump(m, fh)
    comm.Barrier()
    run = run_study_file(ypath, output_dir=ydir / "out", comm=comm)
    out["yaml"] = run.summary
    out["t_yaml"] = time.perf_counter() - t_y
    del run
    gc.collect()

    # Arm M (test side): the kept rod, 403 -> gel, field P1, the same private core.
    t2 = time.perf_counter()
    built_m = _build_kept_rod(leg_sheet_strip_fraction=GATED_WIDTH_FRACTION)
    msh_m, ct_m = built_m[0], built_m[1]
    c_m = {
        "cells": int(msh_m.topology.index_map(msh_m.topology.dim).size_global),
        "t3": _owned_tag_count(msh_m, ct_m, PHANTOM_CELL_TAG),
        "t403": _owned_tag_count(msh_m, ct_m, BIRDCAGE_IMPLANT_CELL_TAG),
    }
    out["c_m"] = c_m
    _p0(
        f"[WF-15 s14] (A0m) arm M cells {c_m['cells']}, tag-3 {c_m['t3']} (arm R {out['t3_R']}), "
        f"tag-403 {c_m['t403']}; expected global {ROD_EXACT_CELL_RECORD} + {c_m['t403']}"
    )
    # (A0m) asserted before any thermal call on M (collective values: same on every rank).
    assert c_m["t3"] == out["t3_R"], (c_m, out["t3_R"])
    assert c_m["cells"] == ROD_EXACT_CELL_RECORD + c_m["t403"], c_m
    heat = {}

    def hook_m(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        merged = _merge_rod_into_phantom(msh, cell_tags)
        fm = SimpleNamespace(e_complex=e_of[DRIVE], sigma_field=sigma_field)
        pm = float(
            np.real(
                mean_sar(
                    e_of[DRIVE],
                    sigma=sigma_field,
                    rho=PHANTOM_RHO_KG_PER_M3,
                    cell_tags=merged,
                    comm=comm,
                    subdomain_ids=PHANTOM_CELL_TAG,
                )["dissipated_power_w"]
            )
        )
        heat["P"] = pm
        heat["main"] = study_module._phantom_transient_core(
            st,
            msh,
            merged,
            comm,
            fm,
            _cooling(G),
            _cn(),
            MASS_G,
            [tuple(c) for c in TIPS.values()],
            None,
            lambda s: study_module._rod_outer_robin_tags(s, dict(ROD)),
            c_m["cells"],
            DRIVE,
        )

    _arm_with_403(built_m, "WF-15 s14 arm M (rod as gel)", hook_m, map_gel=True)
    out["M"] = heat
    out["t_m"] = time.perf_counter() - t2
    out["t_all"] = time.perf_counter() - t0
    return out


@pytest.fixture(scope="module")
def s14():
    if os.environ.get(KNOB) != "1":
        pytest.skip(f"{KNOB} unset")
    return _run()


@needs_knob
@complex_only
def test_s14_all_gates(s14):
    """Collect every gate before asserting (item 668), so one red hides no other."""
    comm = MPI.COMM_WORLD
    res, main, be, be10 = s14["res"], s14["main"], s14["be"], s14["be10"]
    P, G = s14["P"], s14["G"]
    counts = s14["counts"]
    gates = {}
    gates["A0"] = res.cells == ROD_EXACT_CELL_RECORD
    gates["W0 ext=outer+wall"] = counts["exterior"] == counts["outer"] + counts["wall"]
    gates["W0 wall=parent402"] = counts["wall"] == s14["parent_wall"]
    gates["W0 wall>0"] = counts["wall"] > 0
    err_b = abs(main.source_power_w / (G * P) - 1.0)
    gates["B"] = err_b <= B_BAND
    gates["E cn"] = main.max_balance_residual <= E_BAND
    gates["E be"] = be.max_balance_residual <= E_BAND
    v_sub = None
    # (S2) volume slope: rho c theta-bar^1 / (dt1 G P / V).
    import ufl
    from dolfinx import default_scalar_type, fem

    sub, _pid = study_module._tag_submesh(res.mesh, res.cell_tags, PHANTOM_CELL_TAG)
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = float(
        np.real(comm.allreduce(fem.assemble_scalar(fem.form(one * ufl.dx(domain=sub))), op=MPI.SUM))
    )
    s2 = abs(RHO_C * be.volume_mean_k[1] / (DT_1 * G * P / v_sub) - 1.0)
    s2_10 = abs(RHO_C * be10.volume_mean_k[1] / (DT_1 / 10.0 * G * P / v_sub) - 1.0)
    gates["S2"] = s2 <= S2_BAND
    # (S2-tip), re-posed (§9 item 884, U127′): the consistent-mass form.  R_Pi =
    # |rho c <theta^1>_B / (dt <Pi q>_B) - 1| at dt1/10 is the gate; R_Pi at dt1,
    # R_DG0 (the retired form) and <Pi q>/<q> - 1 are printed.  Helpers are item
    # 870's, imported by name; the radius is the Study's own.
    comps = _ball_comparands(s14, sub)
    names = list(TIPS)
    rows10 = _tip_slope_rows(
        {n: be10.balls[i]["theta_k"][1] for i, n in enumerate(names)}, DT_1 / 10.0, comps, RHO_C
    )
    rows1 = _tip_slope_rows(
        {n: be.balls[i]["theta_k"][1] for i, n in enumerate(names)}, DT_1, comps, RHO_C
    )
    # (B-vol), a precondition: both ball means integrate over one set.
    bvol = {}
    for i, n in enumerate(names):
        bvol[n] = abs(comps[n]["vol"] / float(be10.balls[i]["volume_m3"]) - 1.0)
        gates[f"B-vol {n}"] = bvol[n] <= 1.0e-12
    gates["B-vol quad 16=16"] = int(study_module.SAR_QUADRATURE_DEGREE) == int(QUADRATURE_DEGREE)
    for n in names:
        gates[f"S2-tip {n}"] = rows10[n]["R_Pi"] <= S2_BAND
    gates["S0"] = s14["solves"] == 1
    th_r = _theta_tip(main)
    th_m = _theta_tip(s14["M"]["main"])
    th_all = _theta_tip(s14["allext"])
    for name in TIPS:
        gates[f"c-M {name}"] = abs(th_r[name] / th_m[name] - 1.0) > IDENTITY_BAND
    # (Y) the YAML's tip theta-bar(t*) equals the Study's.
    ys = s14["yaml"]["temperature"]
    th_y = {n: float(ys["balls"][i]["theta_k"][-1]) for i, n in enumerate(names)}
    err_y = {n: abs(th_y[n] / th_r[n] - 1.0) for n in names}
    for n in names:
        gates[f"Y {n}"] = err_y[n] <= IDENTITY_BAND
    gates["Y drive/cells"] = ys["drive"] == DRIVE and s14["yaml"]["cells"] == ROD_EXACT_CELL_RECORD

    _p0(
        f"\n[WF-15 s14] -n {comm.size}: EM (solve + sar) {s14['t_em']:.1f} s, three temperature "
        f"calls {s14['t_th']:.1f} s, arm M {s14['t_m']:.1f} s, total {s14['t_all']:.1f} s"
    )
    _p0(
        f"[WF-15 s14] (A0) cells {res.cells} (record {ROD_EXACT_CELL_RECORD}); (W0) exterior "
        f"{counts['exterior']} = outer {counts['outer']} + wall {counts['wall']}; parent tag-402 "
        f"owned {s14['parent_wall']}"
    )
    _p0(
        f"[WF-15 s14] drive {DRIVE} at 1 V: P_phantom {P:.12e} W, G {G:.12e} 1/W; (B) int q "
        f"{main.source_power_w:.15e} vs G*P {G * P:.15e}: rel {err_b:.3e} (B_BAND {B_BAND:.0e})"
    )
    _p0(
        f"[WF-15 s14] (E) CN {main.max_balance_residual:.3e}, BE {be.max_balance_residual:.3e} "
        f"(E_BAND {E_BAND:.0e}); (S0) port solves {s14['solves']}"
    )
    _p0(
        f"[WF-15 s14] (S2) volume slope |r-1| {s2:.4e} (S2_BAND {S2_BAND:.0e}); at dt1/10 "
        f"{s2_10:.4e} (printed); V_sub {v_sub:.12e}"
    )
    for name in TIPS:
        _p0(
            f"[WF-15 s14p] (S2-tip) {name}: R_Pi at dt1/10 {rows10[name]['R_Pi']:.4e} "
            f"(S2_BAND {S2_BAND:.0e}, asserted); R_Pi at dt1 {rows1[name]['R_Pi']:.4e}, ratio "
            f"dt1/(dt1/10) {rows1[name]['R_Pi'] / rows10[name]['R_Pi']:.3f} (printed); retired form "
            f"R_DG0 at dt1/10 {rows10[name]['R_DG0']:.4e}, at dt1 {rows1[name]['R_DG0']:.4e} "
            f"(printed); <q>_B {comps[name]['q']:.12e} W/m^3, <Pi q>_B/<q>_B - 1 "
            f"{rows10[name]['pq_q']:+.4e}; (B-vol) helper vol {comps[name]['vol']:.12e} m^3 vs "
            f"Study {s14['be10'].balls[names.index(name)]['volume_m3']:.12e}: rel {bvol[name]:.3e}"
        )
    _p0(
        f"[WF-15 s14p] (B-vol) quadrature Study {int(study_module.SAR_QUADRATURE_DEGREE)} = "
        f"helper {int(QUADRATURE_DEGREE)}; predicted (not asserted) R_Pi at dt1/10 at the 1e-5 level"
    )
    for name in TIPS:
        _p0(
            f"[WF-15 s14] record {name} tip 1 g theta(t*) R {th_r[name]:.12e} K, M "
            f"{th_m[name]:.12e} K: R/M = {th_r[name] / th_m[name]:.6f} (I41 {I41_TIP_RATIO[name]} "
            f"under ccw on the ragged mesh; here exact rectangles, drive {DRIVE}); (c-M) "
            f"|R/M-1| {abs(th_r[name] / th_m[name] - 1.0):.3e} > IDENTITY_BAND {IDENTITY_BAND:.0e}, "
            f"predicted size printed; (c-Robin) all-exterior Robin {th_all[name]:.12e} K, "
            f"change {th_all[name] / th_r[name] - 1.0:+.4e}"
        )
    _p0(f"[WF-15 s14] arm M P_phantom {s14['M']['P']:.12e} W (G from arm R)")
    for n in names:
        _p0(
            f"[WF-15 s14p] (Y) {n}: YAML theta(t*) {th_y[n]:.15e} K vs Study {th_r[n]:.15e} K: "
            f"rel {err_y[n]:.3e} (IDENTITY_BAND {IDENTITY_BAND:.0e}); YAML run {s14['t_yaml']:.1f} s"
        )
    red = {k: v for k, v in gates.items() if not v}
    _p0(f"[WF-15 s14] gates: {len(gates) - len(red)}/{len(gates)} green; red: {sorted(red)}")
    assert not red, red


def _ball_comparands(s14, sub):
    """U127′'s tip comparands (§9 item 884): q = G*1/2 sigma|E|^2 built on ``sub2``
    exactly as the retired ``_ball_qbar`` did, then item 870's
    ``_projected_ball_comparands`` with V1 = P1 (the Study's thermal degree) and the
    Study's own ball radius; tag ``"u127p_R"`` (a reused PETSc prefix collides)."""
    import ufl
    from dolfinx import fem

    res, G = s14["res"], s14["G"]
    f = res.fields[DRIVE]
    Wp, avg = study_module._dg0_cell_average(
        res.mesh,
        res.cell_tags,
        0.5 * f.sigma_field * ufl.inner(f.e_complex, f.e_complex),
        PHANTOM_CELL_TAG,
    )
    avg *= G
    sub2, pid = study_module._tag_submesh(res.mesh, res.cell_tags, PHANTOM_CELL_TAG)
    q = study_module._dg0_to_submesh(Wp, avg, sub2, pid)
    V1 = fem.functionspace(sub2, ("Lagrange", 1))
    radius = float(s14["be10"].balls[0]["radius_m"])
    return _projected_ball_comparands(sub2, q, V1, TIPS, radius, "u127p_R")
