"""GEO-38 P61 (§9 item 777, §10 chain P): P33's body — the ASTM gel box inside
F-clinical-1.5T heated by the coil, composed and stub-run before ``xl`` entry 59.

Brief: ``docs/planning/chunks/GEO-38.md`` "§9 item 777".  Composition, no
imported module edited:

* **Field half:** C6a's `_run_box` (arm R, tag 403 mapped to gel, no rod),
  unchanged.  Its fields are taken by wrapping two of C6a's module globals for
  the call: `_box_readouts` (hands over ``e_ccw``, ``e_cw``, σ and C6's c*) and
  `_p_acc_spec` (hands over S, V_s, z0 for the lumped line cut).
* **Thermal half:** I32a's `_thermal_arm` (its `_heat_source` q = ½σ|E|², the
  ρ·SAR of peak phasors; MAT-15's `_cell_average_dg0_tag` / `_tag_submesh` /
  `_transient`, CN θ = 0.5 with two Rannacher half-steps then the BE step at
  ``DT_1_BOX``; Robin ``H_GEL`` on the box faces; t* = ``T_STAR_BOX``, Δt = t*/40;
  P37's six points) on tags (3, 403) merged, with `_transient` wrapped to keep
  every CN snapshot so the history at the hotspot exists.
* **Dose:** `cem43` and `running_average` (W = 360 s) from
  `fem_em_solver.thermal.dose` on the point history at C6's c* (the ccw SAR10g
  argmax of the 18-point lattice — which hotspot entry 59 reads is the weekly's).

**Sense.** The thermal arm runs the weekly's (1, +1j) — C6's *cw*, I32a's
``SENSE`` — so entries 55 and 59 share a drive; the ccw (1, −1j) source's P_gel
is printed beside it.  **Normaliser.** The arm runs at G = 1 / P_gel(cw)
(I32a's (B) closes on G·P = 1); every θ, the dose and the running mean are then
reported per accepted watt, θ_acc = θ · P_gel / P_acc(cw), with P_acc the lumped
line cut (ruling 2026-10-12) and the factor printed.  T = ``T_INF`` + θ_acc.

**Sense knob (P61b, §9 item 818).** ``FEM_EM_GEO38_C6_HEAT_SENSE``: unset or empty
⇒ ``("cw",)`` (the landed behaviour above); ``both`` ⇒ ``("ccw", "cw")``, ccw first;
``ccw`` / ``cw`` ⇒ that sense alone; anything else ⇒ ``ValueError``.  Each sense runs
its own thermal arm at G = 1 / P_gel(s) at the same c*, with its own per-accepted-
watt factor P_gel(s) / P_acc(s) and dose record; whenever cw runs, the top-level
``arm`` / ``dose`` / ``scale`` are cw's.

Selected by ``FEM_EM_GEO38_C6_HEAT`` only: unset ⇒ the gated test skips and the
coding tests run (no gmsh); ``1`` ⇒ the body, with C6's node knob
``FEM_EM_GEO36_S8_DRIVE_C_PF`` required; anything else ⇒ ``ValueError``.  Degree
1 by default (I32a's `_refuse_degree` still refuses the rod's S8 knob); P87 (§9
item 969): C34's ``FEM_EM_GEO38_C6_DEGREE=2`` routes degree 2 through
`_p61_degree` to C6a's `_run_box(degree=2)`, and C52's 3 T knob
``FEM_EM_GEO38_C21_DEGREE`` is refused here.  This is the 1.5 T fixture: C6a's
64 MHz ω is right here.

Scope: the body and its stub proof.  Not claimed: any temperature on the
clinical fixture (that is entry 59), any solve.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.thermal.dose import cem43, running_average
from tests.complex_mode import complex_only
from tests.validation import test_geo36_rod_f_clinical_heating as i32
from tests.validation import test_geo38_f_clinical_box_loaded as c6
from tests.validation.test_imp3_rod_heating_matched import (
    _merge_rod_into_phantom,
    _owned_tag_count,
)
from tests.validation.test_mat10_sar_to_heat_sphere import B_BAND, _allsum
from tests.validation.test_mat15_coil_driven_heating_phantom import (
    E_BAND,
    N_MAIN,
    S2_BAND,
    T_INF,
)
from tests.validation.test_mat20_cem43_closed_form import BAR as DOSE_BAR
from tests.validation.test_mat21_box_fem import DT_1 as DT_1_BOX
from tests.validation.test_mat21_box_fem import T_STAR as T_STAR_BOX
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG

HEAT_ENV = "FEM_EM_GEO38_C6_HEAT"
HEAT_SENSE_ENV = "FEM_EM_GEO38_C6_HEAT_SENSE"
SENSE = "(1,+1j)"
#: The drive tuple each sense key prints (C6's K5: ccw (1, -1j), cw (1, +1j)).
SENSE_TUPLE = {"ccw": "(1,-1j)", "cw": "(1,+1j)"}
#: (R2-c), item 818: |P_acc(ccw) / P_acc(cw) - 1| on C6a's `_k8_stubs`, from
#: 3.425774915606e-02 / 3.416698509346e-02 - 1 (`20261103T104609Z_GEO-38-P61-W0c.log:129`).
SENSE_NORMALISER_MISS = 2.6565e-03
#: MAT-20's running-mean window (s).
RUNNING_WINDOW_S = 360.0
#: Variables a P61 window never reads: C6a's own knob, its foreign list, I32a's knob.
FOREIGN_ENVS = (c6.C6_ENV, i32.HEAT_ENV) + tuple(c6.FOREIGN_ENVS)
#: Named controls (brief), I32a's.
PERTURBED_SOURCE = i32.PERTURBED_SOURCE
SINGLE_COLUMN_MISS = i32.SINGLE_COLUMN_MISS
SINGLE_COLUMN_CEILING = i32.SINGLE_COLUMN_CEILING


def _say(comm, msg):
    if comm.rank == 0:
        print(f"[GEO-38 P61] {msg}", flush=True)


def _heat_mode():
    raw = os.environ.get(HEAT_ENV, "").strip()
    if raw == "":
        return None
    if raw == "1":
        return raw
    raise ValueError(f"{HEAT_ENV} must be unset or '1', got {raw!r}")


def _sense_knob():
    """The thermal senses, in run order (ccw first): unset/empty => ``("cw",)``."""
    raw = os.environ.get(HEAT_SENSE_ENV, "").strip()
    table = {"": ("cw",), "both": ("ccw", "cw"), "ccw": ("ccw",), "cw": ("cw",)}
    if raw not in table:
        raise ValueError(f"{HEAT_SENSE_ENV} must be unset, 'both', 'ccw' or 'cw', got {raw!r}")
    return table[raw]


def _knobs():
    """Foreign guard, degree refusal, then ``(mode, c_pf)`` (C6's node knob)."""
    from tests.validation import test_geo36_rod_f_clinical_drive as i36

    bad = [k for k in FOREIGN_ENVS if os.environ.get(k) is not None]
    if bad:
        raise ValueError(f"P61 never reads these variables; unset them: {bad}")
    mode = _heat_mode()
    i32._refuse_degree()
    return mode, i36._node_knob(mode)


#: P87 (§9 item 969): C52's 3 T degree knob, refused by a 1.5 T P61 window.
C21_DEGREE_ENV_3T = "FEM_EM_GEO38_C21_DEGREE"


def _p61_degree(mode):
    """P87: C34's ``FEM_EM_GEO38_C6_DEGREE`` through `c6._c6_degree` (unset / empty ⇒
    ``None``, ``"2"`` ⇒ 2, anything else ⇒ ValueError); C52's 3 T knob refused."""
    if os.environ.get(C21_DEGREE_ENV_3T) is not None:
        raise ValueError(f"P61 never routes a 3 T degree; unset {C21_DEGREE_ENV_3T}")
    return c6._c6_degree(mode)


def _probes():
    """P37's six points in the fixture frame, on the box's own y offset."""
    return i32._p37_fixture(float(c6._box()[0][1]))


def _thermal_half(label, msh, merged, e_sup, sigma_f, centre, tips, probes, hot, *, sense=SENSE):
    """I32a's `_thermal_arm` (G = 1 / P_gel), every CN snapshot kept, then the point
    history at ``hot``.  Returns the arm dict plus ``times``, ``hist`` (θ per gel
    watt) and ``probe_last`` (θ at ``hot`` read on the final field).  ``sense`` is
    the drive tuple I32a's print lines carry (its ``SENSE``; a label only)."""
    keep = {}
    real_transient, real_probe = i32._transient, i32._probe_theta

    def _transient_all(sub, ft, q, **kw):
        if kw.get("theta") == 0.5:
            kw["snapshot_steps"] = tuple(range(N_MAIN + 1))
            res = real_transient(sub, ft, q, **kw)
            keep["res"] = res
            return res
        return real_transient(sub, ft, q, **kw)

    def _probe_kept(T, pts, parent):
        keep["parent"] = parent
        return real_probe(T, pts, parent)

    with pytest.MonkeyPatch.context() as m:
        m.setattr(i32, "_transient", _transient_all)
        m.setattr(i32, "_probe_theta", _probe_kept)
        m.setattr(i32, "SENSE", sense)
        arm = i32._thermal_arm(
            label, msh, merged, e_sup, sigma_f, None, centre, tips, probes, merged=True
        )
    res, parent = keep["res"], keep["parent"]
    # P87: the volume-mean theta(t*) on the thermal submesh (per deposited watt, G·P = 1);
    # real part, allreduced through `_allsum`.
    sub = res.temperature.function_space.mesh
    one = fem.Constant(sub, default_scalar_type(1.0))
    v_sub = float(np.real(_allsum(sub.comm, one * ufl.dx(domain=sub))))
    vol_mean = (
        float(np.real(_allsum(sub.comm, (res.temperature - T_INF) * ufl.dx(domain=sub)))) / v_sub
    )
    arm.update({"vol_mean": vol_mean, "v_sub": v_sub})
    hot = np.asarray(hot, dtype=np.float64).reshape(1, 3)
    hist, ok = [], []
    for n in range(N_MAIN + 1):
        th, valid = real_probe(res.snapshots[n], hot, parent)
        hist.append(float(th[0]))
        ok.append(bool(valid[0]))
    last, _ok = real_probe(res.temperature, hot, parent)
    res.snapshots.clear()
    arm.update(
        {
            "times": np.arange(N_MAIN + 1) * arm["dt"],
            "hist": np.asarray(hist),
            "hist_ok": ok,
            "probe_last": float(last[0]),
        }
    )
    return arm


def _dose_record(arm, scale):
    """θ_acc = scale · θ; cem43 on T_INF + θ_acc and the 360 s running mean of θ_acc."""
    th_acc = scale * arm["hist"]
    dose = cem43(arm["times"], T_INF + th_acc)
    _means, ra_max, ra_tmax = running_average(arm["times"], th_acc, window=RUNNING_WINDOW_S)
    return {
        "th_acc": th_acc,
        "cem43": float(dose),
        "ra_max": float(ra_max),
        "ra_tmax": float(ra_tmax),
    }


def _p_gel(e, sigma_f, merged, comm):
    from fem_em_solver.post.sar import mean_sar
    from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3

    return float(
        np.real(
            mean_sar(
                e,
                sigma=sigma_f,
                rho=PHANTOM_RHO_KG_PER_M3,
                cell_tags=merged,
                comm=comm,
                subdomain_ids=PHANTOM_CELL_TAG,
            )["dissipated_power_w"]
        )
    )


def _run_heat_box(comm, monkeypatch, c_pf, *, senses=("cw",), degree=None):
    """Entry 59's body: C6a's `_run_box` (gates asserted inside), then the gel twin's
    thermal arm on each of ``senses`` (ccw first), the dose at c*, per accepted watt.
    Gates collected per sense, then asserted."""
    from fem_em_solver.io.mesh import BIRDCAGE_IMPLANT_CELL_TAG
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation.test_wf10_f_human_tuned_b1_sar import _sense_weights

    senses = tuple(senses)
    if not senses or len(set(senses)) != len(senses) or any(s not in SENSE_TUPLE for s in senses):
        raise ValueError(f"senses must be distinct keys of {tuple(SENSE_TUPLE)}, got {senses!r}")
    i32._box_preconditions()
    cap = {}
    real_readouts, real_p_acc = c6._box_readouts, c6._p_acc_spec

    def _readouts(msh, cell_tags, e_ccw, e_cw, sigma_f, omega, *a, **k):
        out = real_readouts(msh, cell_tags, e_ccw, e_cw, sigma_f, omega, *a, **k)
        cap.update(
            msh=msh,
            cell_tags=cell_tags,
            e_ccw=e_ccw,
            e_cw=e_cw,
            sigma=sigma_f,
            omega=omega,
            c_star=out["c_star"],
            star_label=out["star_label"],
        )
        return out

    def _p_acc(s2, w, v_s, z0):
        cap.update(s2=np.asarray(s2), v_s=v_s, z0=z0)
        return real_p_acc(s2, w, v_s, z0)

    _say(
        comm,
        f"\nstart: -n {comm.size}, P61 body; node C {c_pf:.4f} pF, degree {degree or 1}; field half C6a "
        f"_run_box (arm R, gel twin); thermal senses {list(senses)} (ccw (1,-1j), cw {SENSE}), "
        f"ccw P_gel printed beside",
    )
    monkeypatch.setattr(c6, "_box_readouts", _readouts)
    monkeypatch.setattr(c6, "_p_acc_spec", _p_acc)
    c6._run_box(comm, monkeypatch, c_pf, degree=degree)
    monkeypatch.setattr(c6, "_box_readouts", real_readouts)
    monkeypatch.setattr(c6, "_p_acc_spec", real_p_acc)

    msh, cell_tags = cap["msh"], cap["cell_tags"]
    w = _sense_weights()
    p_acc = {s: real_p_acc(cap["s2"], w[s], cap["v_s"], cap["z0"]) for s in ("ccw", "cw")}
    merged = _merge_rod_into_phantom(msh, cell_tags)
    merged_rod = _owned_tag_count(msh, merged, BIRDCAGE_IMPLANT_CELL_TAG)
    p_gel_ccw = _p_gel(cap["e_ccw"], cap["sigma"], merged, comm)
    hot = cap["c_star"]  # C6's c*, the ccw lattice argmax, for both senses
    tips = {f"c* {cap['star_label']}": hot}
    by_sense, fails = {}, []
    for s in senses:
        st = SENSE_TUPLE[s]
        arm = _thermal_half(
            f"gel twin {s}",
            msh,
            merged,
            cap[f"e_{s}"],
            cap["sigma"],
            i35._rod_centre(),
            tips,
            _probes(),
            hot,
            sense=st,
        )
        scale = arm["P"] / p_acc[s]
        dose = _dose_record(arm, scale)
        by_sense[s] = {"arm": arm, "dose": dose, "scale": scale}
        _say(
            comm,
            f"sources, sense {st} ({s}, heated): P_gel {arm['P']:.12e} W, P_acc {p_acc[s]:.12e} W "
            f"(lumped line cut, ruling 2026-10-12); per-accepted-watt factor P_gel/P_acc ({s}) "
            f"{scale:.9e}",
        )
        _say(
            comm,
            f"hotspot {list(tips)[0]} {tuple(round(v, 6) for v in hot)}, sense {st}: theta(t*) per "
            f"P_acc watt {dose['th_acc'][-1]:.9e} K (point; final-field read "
            f"{scale * arm['probe_last']:.9e} K); cem43 per P_acc watt {dose['cem43']:.9e} s "
            f"equivalent at 43 C (the function's units: times in s); running mean (W "
            f"{RUNNING_WINDOW_S:.0f} s) max {dose['ra_max']:.9e} K at t {dose['ra_tmax']:.1f} s "
            f"(printed)",
        )
        for n, v in arm["probe_th"].items():
            _say(
                comm,
                f"probe {n}, sense {st}: theta(t*) per P_acc watt {scale * v:+.9e} K (printed)",
            )
        fails += [f"{s}: {f}" for f in i32._thermal_failures(arm)]
        if not all(arm["hist_ok"]):
            fails.append(f"{s}: hotspot history not located")
    # P87: one volume-mean line per sense, printed after every sense's thermal block so
    # the pre-existing lines keep their positions (K10-R2's (R2-o) indices are (L) tokens).
    for s in senses:
        arm, scale = by_sense[s]["arm"], by_sense[s]["scale"]
        vm = arm["vol_mean"]
        frac = i32.RHO_C * arm["v_sub"] * vm / (arm["G"] * arm["P"] * T_STAR_BOX)
        arm["stored_frac"] = frac
        _say(
            comm,
            f"volume-mean theta(t*), sense {SENSE_TUPLE[s]}: {vm:.9e} K per deposited watt "
            f"(ANS-13 row 5); per P_acc watt {scale * vm:.9e} K; stored fraction rho c V <theta>/"
            f"(G P t*) {frac:.6e} (printed)",
        )
    _say(
        comm,
        f"sources: P_gel ccw (1,-1j) {p_gel_ccw:.12e} W (printed beside); P_acc cw "
        f"{p_acc['cw']:.12e} W, ccw {p_acc['ccw']:.12e} W (lumped line cut, ruling 2026-10-12)",
    )
    if merged_rod != 0:
        fails.append(f"merged tag 403 {merged_rod}")
    _say(comm, f"gates: failures {fails or 'none'} (ASSERTED)")
    assert not fails, fails
    top = by_sense["cw"] if "cw" in by_sense else by_sense[senses[-1]]
    return {
        "arm": top["arm"],
        "dose": top["dose"],
        "p_acc": p_acc,
        "p_gel_ccw": p_gel_ccw,
        "scale": top["scale"],
        "merged_rod": merged_rod,
        "by_sense": by_sense,
        "degree": degree or 1,
    }


@complex_only
def test_geo38_p61_box_heating(monkeypatch):
    mode, c_pf = _knobs()  # node knob missing under 1 => ValueError, before gmsh
    senses = _sense_knob()  # a bad sense raises before the mode skip
    deg = _p61_degree(mode)  # P87: a bad degree / the 3 T knob raises before the mode skip
    if mode is None:
        pytest.skip(f"{HEAT_ENV} unset: GEO-38 P61 is selected by environment only")
    _run_heat_box(MPI.COMM_WORLD, monkeypatch, c_pf, senses=senses, degree=deg)


# ---------------------------------------------------------------------------
# Coding tests (no gmsh); skip whenever HEAT_ENV is set.
# ---------------------------------------------------------------------------
_coding = pytest.mark.skipif(
    os.environ.get(HEAT_ENV, "").strip() != "",
    reason=f"{HEAT_ENV} set: coding tests skip in windows",
)


@_coding
def test_p61_k0_knobs_and_box_time_base(monkeypatch):
    from tests.validation.test_geo36_rod_f_clinical_drive import (
        S8_DRIVE_C_PF_ENV,
        S8_DRIVE_DEGREE_ENV,
    )

    for k in (HEAT_ENV, S8_DRIVE_C_PF_ENV, S8_DRIVE_DEGREE_ENV) + FOREIGN_ENVS:
        monkeypatch.delenv(k, raising=False)
    assert _knobs() == (None, None)
    monkeypatch.setenv(HEAT_ENV, "1")
    with pytest.raises(ValueError, match=S8_DRIVE_C_PF_ENV):
        _knobs()
    monkeypatch.setenv(S8_DRIVE_C_PF_ENV, "193.5996")
    assert _knobs() == ("1", 193.5996)
    for bad in ("0", "build", "x"):
        monkeypatch.setenv(HEAT_ENV, bad)
        with pytest.raises(ValueError):
            _heat_mode()
    monkeypatch.setenv(HEAT_ENV, "1")
    monkeypatch.setenv(S8_DRIVE_DEGREE_ENV, "2")
    with pytest.raises(ValueError):
        _knobs()
    with pytest.raises(ValueError, match="refused"):
        i32._refuse_degree()
    monkeypatch.delenv(S8_DRIVE_DEGREE_ENV)
    for k in FOREIGN_ENVS:
        with monkeypatch.context() as m:
            m.setenv(k, "1")
            with pytest.raises(ValueError):
                _knobs()
    i32._box_preconditions()
    print(
        f"\n[GEO-38 P61 K0] modes unset/1 accepted, 0/build/x refused; degree 2 refused; node "
        f"required under 1; {len(FOREIGN_ENVS)} foreign variables refused; box time base T_STAR "
        f"{T_STAR_BOX:.3f} s, dt {T_STAR_BOX / N_MAIN:.6f} s, dt_1 {DT_1_BOX:.6f} s (ASSERTED)",
        flush=True,
    )
    assert abs(T_STAR_BOX - 2821.5) <= 1.0, T_STAR_BOX


def _k10_field(comm, factor):
    """`_k9_fixture`'s mesh, e = factor · E₀ x̂ on its N1curl space."""
    from tests.validation.test_geo36_rod_f_clinical_drive import _k9_fixture

    msh, cell_tags, e_u, sigma_f, tips, _t1 = _k9_fixture(comm)
    e = fem.Function(e_u.function_space)
    e.interpolate(lambda x: np.vstack([factor * i32.STUB_E0 + 0.0 * x[0], 0.0 * x[0], 0.0 * x[0]]))
    e.x.scatter_forward()
    return msh, cell_tags, e, sigma_f, tips


@complex_only
@_coding
def test_p61_k10q_source_closed_form():
    """(K10-q): e = (1 + j)E₀x̂ ⇒ q = σE₀² on every cell, and SAR·ρ = q.  Control:
    e₂ = 0 (e = E₀x̂) misses by 0.5."""
    from fem_em_solver.post.sar import averaging_ball_radius, build_density_field, mass_averaged_sar
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL
    from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3
    from tests.validation.test_mass_averaged_sar_standard_masses import (
        ONE_GRAM_KG,
        QUADRATURE_DEGREE,
    )

    comm = MPI.COMM_WORLD
    msh, cell_tags, e, sigma_f, tips = _k10_field(comm, 1.0 + 1.0j)
    q_exact = i32.STUB_SIGMA * i32.STUB_E0**2
    _m, _s, _i, q = i32._merged_q(msh, cell_tags, e, sigma_f)
    lo, hi = i32._owned_extrema(q)
    miss = max(abs(lo / q_exact - 1.0), abs(hi / q_exact - 1.0))
    a_1g = averaging_ball_radius(mass_kg=ONE_GRAM_KG, rho=PHANTOM_RHO_KG_PER_M3)
    ball = mass_averaged_sar(
        e,
        sigma=sigma_f,
        rho=build_density_field(msh, PHANTOM_RHO_KG_PER_M3),
        center=tips["a"],
        radius=a_1g,
        quadrature_degree=QUADRATURE_DEGREE,
    )
    r_sar = abs(float(ball["averaged_sar_w_per_kg"]) * PHANTOM_RHO_KG_PER_M3 / q_exact - 1.0)
    e_one = fem.Function(e.function_space)  # e₂ = 0: the (1 + j) field's first column alone, E₀x̂
    e_one.x.array[:] = e.x.array / (1.0 + 1.0j)
    e_one.x.scatter_forward()
    _m, _s, _i, q1 = i32._merged_q(msh, cell_tags, e_one, sigma_f)
    lo_c, hi_c = i32._owned_extrema(q1)
    miss_c = max(abs(lo_c / q_exact - 1.0), abs(hi_c / q_exact - 1.0))
    if comm.rank == 0:
        print(
            f"\n[GEO-38 P61 K10-q] q on [{lo:.12e}, {hi:.12e}] vs sigma E0^2 {q_exact:.12e}: miss "
            f"{miss:.3e}; SAR*rho miss {r_sar:.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); control "
            f"e2 = 0 miss {miss_c:.12f} (named {SINGLE_COLUMN_MISS}, within {SINGLE_COLUMN_CEILING:g}, "
            f"ASSERTED)",
            flush=True,
        )
    assert miss <= DISCRETE_IDENTITY_RTOL, miss
    assert r_sar <= DISCRETE_IDENTITY_RTOL, r_sar
    assert abs(miss_c - SINGLE_COLUMN_MISS) <= SINGLE_COLUMN_CEILING, miss_c


@complex_only
@_coding
def test_p61_k10es_energy_slope_and_history_on_uniform_source(monkeypatch):
    """(K10-es): `_thermal_half` on the uniform source: (B), (E), (E′) in band, the
    slope within 0.5·``S2_BAND``; the history starts at 0, is monotone, and its last
    point equals the final-field read.  Control: q(1 + 1e-3) misses (B) by 1e-3."""
    comm = MPI.COMM_WORLD
    msh, cell_tags, e, sigma_f, tips = _k10_field(comm, 1.0 + 1.0j)
    merged = _merge_rod_into_phantom(msh, cell_tags)
    hot = tips["a"]
    with monkeypatch.context() as m:
        m.setattr(i32, "WALL_SELECT_M", i32.STUB_WALL_SELECT_M)
        arm = _thermal_half(
            "K10", msh, merged, e, sigma_f, i32.STUB_CENTRE, {"a": hot}, i32.STUB_PROBES, hot
        )
        real_source = i32._heat_source

        def _perturbed(*a):
            q = real_source(*a)
            q.x.array[:] *= 1.0 + PERTURBED_SOURCE
            return q

        m.setattr(i32, "_heat_source", _perturbed)
        neg = _thermal_half(
            "K10 NEG", msh, merged, e, sigma_f, i32.STUB_CENTRE, {"a": hot}, i32.STUB_PROBES, hot
        )
    fails, neg_fails = i32._thermal_failures(arm), i32._thermal_failures(neg)
    h = arm["hist"]
    r_last = abs(h[-1] / arm["probe_last"] - 1.0)
    if comm.rank == 0:
        print(
            f"\n[GEO-38 P61 K10-es] (B) {arm['B']:.3e}, (E) {arm['E']:.3e}, (E') {arm['Ep']:.3e} (bands "
            f"{B_BAND:.0e} / {E_BAND:.0e}); slope (S2) {arm['S2']:.4e} (<= 0.5 x {S2_BAND:.0e}); dt "
            f"{arm['dt']:.6f} s; history n {len(h)}, theta0 {h[0]:.3e}, last {h[-1]:.9e} vs final read "
            f"{arm['probe_last']:.9e} rel {r_last:.3e}; failures {fails or 'none'} (ASSERTED)"
        )
        print(
            f"[GEO-38 P61 K10-es] control q(1 + 1e-3): (B) {neg['B']:.9e} = {neg['B'] / B_BAND:.4e} x "
            f"B_BAND; failures {neg_fails} (ASSERTED)",
            flush=True,
        )
    assert not fails, fails
    assert arm["dt"] == T_STAR_BOX / N_MAIN, arm["dt"]
    assert arm["S2"] <= 0.5 * S2_BAND, arm["S2"]
    assert arm["counts"]["wall"] == 0
    assert len(h) == N_MAIN + 1 and all(arm["hist_ok"])
    assert abs(h[0]) <= 1e-12, h[0]
    assert np.all(np.diff(h) > 0.0), h
    assert r_last <= 1e-12, r_last
    assert abs(neg["B"] / PERTURBED_SOURCE - 1.0) <= 1.0e-9, neg["B"]
    assert any("(B)" in f for f in neg_fails), neg_fails


@complex_only
@_coding
def test_p61_k10c_constant_curl_closed_form():
    """(K10-C): E₀(−y, x, 0)/L on (1 + j): ∫q = σE₀²(16/3)a⁵/L² over the cube."""
    from tests.validation.test_geo36_rod_f_clinical_drive import _k9_fixture

    comm = MPI.COMM_WORLD
    msh, cell_tags, e_u, sigma_f, _tips, _t1 = _k9_fixture(comm)
    a, length = 0.05, 0.05
    e = fem.Function(e_u.function_space)
    e.interpolate(
        lambda x: (1.0 + 1.0j) * i32.STUB_E0 / length * np.vstack([-x[1], x[0], 0.0 * x[0]])
    )
    e.x.scatter_forward()
    _m, sub, _i, q = i32._merged_q(msh, cell_tags, e, sigma_f)
    int_q = float(np.real(_allsum(comm, q * ufl.dx(domain=sub))))
    exact = i32.STUB_SIGMA * i32.STUB_E0**2 * (16.0 / 3.0) * a**5 / length**2
    rel = abs(int_q / exact - 1.0)
    if comm.rank == 0:
        print(
            f"\n[GEO-38 P61 K10-C] int q {int_q:.15e} vs closed form {exact:.15e}: rel {rel:.3e} (band "
            f"{B_BAND:.0e}, ASSERTED)",
            flush=True,
        )
    assert rel <= B_BAND, rel


@_coding
def test_p61_k10p_point_map_on_the_box():
    """(K10-p): P37's six points on the box's own offset, exactly."""
    from tests.validation.test_mat21_box_closed_form import A_BOX, B_BOX, C_BOX

    off = float(c6._box()[0][1])
    pts = _probes()
    expected = {
        "P0 centre": (0.0, off, 0.0),
        "P1 +y wall": (B_BOX, off, 0.0),
        "P2 -y wall": (-B_BOX, off, 0.0),
        "P3 +x wall": (0.0, off, A_BOX),
        "P4 -x wall": (0.0, off, -A_BOX),
        "P5 top face": (0.0, off + C_BOX, 0.0),
    }
    print(f"\n[GEO-38 P61 K10-p] box y offset {off:.6f} m; mapped {pts} (ASSERTED)", flush=True)
    assert pts == expected, pts


@_coding
def test_p61_k10d_cem43_constant_histories():
    """(K10-D): `cem43` on a constant 44 °C history is 2 × the exposure, on 42 °C
    0.25 ×, in the function's own units; the 360 s running mean of a constant is it."""
    t = np.linspace(0.0, 600.0, 11)
    exposure = t[-1] - t[0]
    r44 = abs(cem43(t, np.full_like(t, 44.0)) / (2.0 * exposure) - 1.0)
    r42 = abs(cem43(t, np.full_like(t, 42.0)) / (0.25 * exposure) - 1.0)
    _means, mx, _tm = running_average(t, np.full_like(t, 3.0), window=RUNNING_WINDOW_S)
    r_ra = abs(mx / 3.0 - 1.0)
    print(
        f"\n[GEO-38 P61 K10-D] cem43 44 C rel {r44:.3e}, 42 C rel {r42:.3e}, running mean of 3 rel "
        f"{r_ra:.3e} (bar {DOSE_BAR:g}, ASSERTED)",
        flush=True,
    )
    assert r44 <= DOSE_BAR, r44
    assert r42 <= DOSE_BAR, r42
    assert r_ra <= DOSE_BAR, r_ra


@complex_only
@_coding
def test_p61_k10r_run_heat_box_end_to_end_on_stubs(monkeypatch, tmp_path):
    """(K10-R): `_run_heat_box` end to end on C6a's `_k8_stubs` (column 2 = 0.5j ×
    column 1, so ccw = 1.5 E and cw = 0.5 E): P_gel ccw / cw = 9; G·P_gel = 1; the
    thermal gates; the per-accepted-watt factor applied to the history exactly."""
    from tests.validation.test_wf10_f_human_tuned_b1_sar import _sense_weights

    comm = MPI.COMM_WORLD
    fx = c6._k6_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    this = sys.modules[_run_heat_box.__module__]
    with monkeypatch.context() as m:
        lines, _seen, _z = c6._k8_stubs(m, fx, root)
        m.setattr(i32, "WALL_SELECT_M", -1.0)  # no rod wall on K6's box: every face is Robin
        m.setattr(this, "_probes", lambda: i32.STUB_PROBES)
        out = _run_heat_box(comm, m, c6.K8_NODE_PF)
    arm, dose = out["arm"], out["dose"]
    r9 = abs(out["p_gel_ccw"] / arm["P"] / 9.0 - 1.0)
    r_gp = abs(arm["G"] * arm["P"] - 1.0)
    r_scale = float(np.max(np.abs(dose["th_acc"] - out["scale"] * arm["hist"])))
    w = _sense_weights()
    rec = c6._k8_json(comm, root)
    r_pacc = abs(out["p_acc"]["ccw"] / rec["p_acc"] - 1.0)
    heads = ("(i) ball radius", "normalisers:", "price:")
    missing = [h for h in heads if not any(ln.startswith(h) for ln in lines)]
    if comm.rank == 0:
        print(
            f"\n[GEO-38 P61 K10-R] P_gel ccw/cw miss from 9 {r9:.3e}; G*P_gel miss {r_gp:.3e}; P_acc "
            f"ccw vs C6a's JSON {r_pacc:.3e}; P_acc cw {out['p_acc']['cw']:.9e} W (weights {w['cw']}); "
            f"scale {out['scale']:.9e}, applied miss {r_scale:.3e}; cem43 {dose['cem43']:.6e}, running "
            f"max {dose['ra_max']:.6e} K; wall {arm['counts']['wall']}, probes valid "
            f"{list(arm['probe_ok'])}; C6a heads missing {missing or 'none'} (ASSERTED)",
            flush=True,
        )
    assert r9 <= B_BAND, r9
    assert r_gp <= B_BAND, r_gp
    assert r_pacc <= 1e-12, r_pacc
    assert r_scale == 0.0, r_scale
    assert out["merged_rod"] == 0
    assert arm["counts"]["wall"] == 0 and all(arm["probe_ok"]) and all(arm["hist_ok"])
    assert not i32._thermal_failures(arm)
    assert np.isfinite(dose["cem43"]) and dose["cem43"] > 0.0
    assert not missing, missing


@_coding
def test_p61b_k0s_sense_knob(monkeypatch):
    """(K0s), item 818: the four accepted values and three refusals."""
    monkeypatch.delenv(HEAT_SENSE_ENV, raising=False)
    got = {"(unset)": _sense_knob()}
    for raw in ("", "both", "ccw", "cw"):
        monkeypatch.setenv(HEAT_SENSE_ENV, raw)
        got[repr(raw)] = _sense_knob()
    refused = []
    for bad in ("x", "CW", "ccw,cw"):
        monkeypatch.setenv(HEAT_SENSE_ENV, bad)
        with pytest.raises(ValueError, match=HEAT_SENSE_ENV):
            _sense_knob()
        refused.append(bad)
    print(f"\n[GEO-38 P61b K0s] {HEAT_SENSE_ENV}: {got}; refused {refused} (ASSERTED)", flush=True)
    assert got == {
        "(unset)": ("cw",),
        "''": ("cw",),
        "'both'": ("ccw", "cw"),
        "'ccw'": ("ccw",),
        "'cw'": ("cw",),
    }, got


@complex_only
@_coding
def test_p61b_k10r2_both_senses_on_stubs(monkeypatch, tmp_path):
    """(K10-R2), item 818: `_run_heat_box(..., senses=("ccw", "cw"))` on C6a's
    `_k8_stubs` (ccw = 1.5 E, cw = 0.5 E, so q_ccw = 9 q_cw).  (R2-p) P_gel ratio 9;
    (R2-h) identical θ-per-gel-watt histories; (R2-s) each scale applied exactly and
    θ_acc ratio 9·P_acc(cw)/P_acc(ccw); (R2-o) the ccw thermal lines print first;
    (R2-g) both arms' thermal gates empty.  Control (R2-c): the ccw arm normalised by
    P_acc(cw) misses (R2-s) by |P_acc(ccw)/P_acc(cw) - 1| > ``B_BAND``."""
    comm = MPI.COMM_WORLD
    fx = c6._k6_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    this = sys.modules[_run_heat_box.__module__]
    said = []
    real_say = this._say

    def _rec(c, msg):
        said.append(msg)
        real_say(c, msg)

    with monkeypatch.context() as m:
        _lines, _seen, _z = c6._k8_stubs(m, fx, root)
        m.setattr(i32, "WALL_SELECT_M", -1.0)  # no rod wall on K6's box: every face is Robin
        m.setattr(this, "_probes", lambda: i32.STUB_PROBES)
        m.setattr(this, "_say", _rec)
        out = _run_heat_box(comm, m, c6.K8_NODE_PF, senses=("ccw", "cw"))
    bs, p_acc = out["by_sense"], out["p_acc"]
    a_cc, a_cw = bs["ccw"]["arm"], bs["cw"]["arm"]
    d_cc, d_cw = bs["ccw"]["dose"], bs["cw"]["dose"]
    r_p = abs(a_cc["P"] / a_cw["P"] / 9.0 - 1.0)
    r_h = float(np.max(np.abs(a_cc["hist"] - a_cw["hist"])) / np.max(np.abs(a_cw["hist"])))
    applied = {
        s: float(np.max(np.abs(bs[s]["dose"]["th_acc"] - bs[s]["scale"] * bs[s]["arm"]["hist"])))
        for s in ("ccw", "cw")
    }
    scale_def = {s: abs(bs[s]["scale"] - bs[s]["arm"]["P"] / p_acc[s]) for s in ("ccw", "cw")}
    predicted = 9.0 * p_acc["cw"] / p_acc["ccw"]
    r_s = abs(d_cc["th_acc"][-1] / d_cw["th_acc"][-1] / predicted - 1.0)
    wrong = a_cc["hist"][-1] * a_cc["P"] / p_acc["cw"]  # control: ccw arm on cw's normaliser
    miss_c = abs(wrong / d_cw["th_acc"][-1] / predicted - 1.0)
    pred_c = abs(p_acc["ccw"] / p_acc["cw"] - 1.0)
    thermal = [i for i, ln in enumerate(said) if ln.startswith(("hotspot ", "probe "))]
    first = {t: next((i for i in thermal if t in said[i]), None) for t in ("(1,-1j)", "(1,+1j)")}
    order_ok = (
        first["(1,-1j)"] is not None
        and first["(1,+1j)"] is not None
        and first["(1,-1j)"] < first["(1,+1j)"]
    )
    g = {s: i32._thermal_failures(bs[s]["arm"]) for s in ("ccw", "cw")}
    top_ok = out["arm"] is a_cw and out["dose"] is d_cw and out["scale"] == bs["cw"]["scale"]
    if comm.rank == 0:
        print(
            f"\n[GEO-38 P61b K10-R2] (R2-p) P_gel ccw {a_cc['P']:.12e} / cw {a_cw['P']:.12e} W miss from 9 "
            f"{r_p:.3e}; (R2-h) max|h_ccw - h_cw|/max|h_cw| {r_h:.3e}; (R2-s) applied {applied}, scale "
            f"defs {scale_def}, theta_acc(t*) ccw {d_cc['th_acc'][-1]:.9e} / cw {d_cw['th_acc'][-1]:.9e} K "
            f"vs 9 P_acc(cw)/P_acc(ccw) {predicted:.12e}: miss {r_s:.3e} (band {B_BAND:.0e}); (R2-o) first "
            f"thermal line (1,-1j) {first['(1,-1j)']} < (1,+1j) {first['(1,+1j)']}: {order_ok}; (R2-g) "
            f"failures {g}; top-level = cw {top_ok} (ASSERTED)"
        )
        print(
            f"[GEO-38 P61b K10-R2] control (R2-c) ccw on P_acc(cw): miss {miss_c:.6e} vs |P_acc(ccw)/"
            f"P_acc(cw) - 1| {pred_c:.6e} (brief {SENSE_NORMALISER_MISS:.4e}); {miss_c / B_BAND:.4e} x "
            f"B_BAND (> B_BAND, ASSERTED)",
            flush=True,
        )
    assert r_p <= B_BAND, r_p
    assert r_h <= B_BAND, r_h
    assert all(v == 0.0 for v in applied.values()), applied
    assert all(v == 0.0 for v in scale_def.values()), scale_def
    assert r_s <= B_BAND, r_s
    assert order_ok, first
    assert not g["ccw"] and not g["cw"], g
    assert all(a_cc["hist_ok"]) and all(a_cw["hist_ok"])
    assert top_ok
    assert miss_c > B_BAND, miss_c


# ---------------------------------------------------------------------------
# P87 (§9 item 969): the degree knob routed through P61; the volume-mean print.
# ---------------------------------------------------------------------------
#: (Deg-arith): I36b's build, E and F (`20261027T170559Z_GEO-36-I36b-W1-build.log.gz:11989`).
ROD_E_15T, ROD_F_15T = 1_606_744, 2_556_953
#: (N-B-d2): |B / PERTURBED_SOURCE - 1|, K10-es's identity (`:409`).
PERTURBED_B_RTOL = 1.0e-9
_VM_HEAD = "volume-mean theta(t*), sense "


def _k6_fixture_d2(comm):
    """C6a's `_k6_fixture` with both fields re-interpolated into ``("N1curl", 2)`` on
    the same mesh (K6's E = factor · ½c x̂ × r lies in both spaces).  Same 7-tuple."""
    msh, cell_tags, _e1, _e2, sigma_f, c, half = c6._k6_fixture(comm)
    v2 = fem.functionspace(msh, ("N1curl", 2))

    def _field(factor):
        e = fem.Function(v2)
        e.interpolate(lambda x: factor * 0.5 * c * np.vstack([0.0 * x[0], -x[2], x[1]]))
        e.x.scatter_forward()
        return e

    return msh, cell_tags, _field(1.0), _field(2.0), sigma_f, c, half


def _element_degree(e):
    """The installed DOLFINx's accessor: ``V.ufl_element().degree``."""
    return int(e.function_space.ufl_element().degree)


def _vm_lines(said):
    return [ln for ln in said if ln.startswith(_VM_HEAD)]


def _vm_acc_printed(line):
    """The per-P_acc-watt token of a volume-mean line."""
    return float(line.split("per P_acc watt ")[1].split(" K")[0])


def _stub_heat_run(comm, m, fx, json_dir, *, degree, senses, recorders):
    """K10-R2's stubs on ``fx`` (C6a's `_k8_stubs`, every face Robin, the stub probes),
    optionally C34's recorders at ``degree``, `_say` recorded; then `_run_heat_box`."""
    this = sys.modules[_run_heat_box.__module__]
    said = []
    real_say = this._say

    def _rec(c, msg):
        said.append(msg)
        real_say(c, msg)

    lines, _seen, _z = c6._k8_stubs(m, fx, json_dir)
    m.setattr(i32, "WALL_SELECT_M", -1.0)  # no rod wall on K6's box: every face is Robin
    m.setattr(this, "_probes", lambda: i32.STUB_PROBES)
    m.setattr(this, "_say", _rec)
    env_seen, degrees = c6._c34_recorders(m, degree=degree) if recorders else (None, None)
    out = this._run_heat_box(comm, m, c6.K8_NODE_PF, senses=senses, degree=degree)
    return out, lines, said, env_seen, degrees


def _r2_checks(out, said):
    """K10-R2's (R2-p) / (R2-h) / (R2-s) / (R2-o) / (R2-g) readings on a both-sense run,
    plus (V-acc) and (V-s) on the volume mean."""
    bs, p_acc = out["by_sense"], out["p_acc"]
    a_cc, a_cw = bs["ccw"]["arm"], bs["cw"]["arm"]
    d_cc, d_cw = bs["ccw"]["dose"], bs["cw"]["dose"]
    r_p = abs(a_cc["P"] / a_cw["P"] / 9.0 - 1.0)
    r_h = float(np.max(np.abs(a_cc["hist"] - a_cw["hist"])) / np.max(np.abs(a_cw["hist"])))
    applied = max(
        float(np.max(np.abs(bs[s]["dose"]["th_acc"] - bs[s]["scale"] * bs[s]["arm"]["hist"])))
        for s in ("ccw", "cw")
    )
    predicted = 9.0 * p_acc["cw"] / p_acc["ccw"]
    r_s = abs(d_cc["th_acc"][-1] / d_cw["th_acc"][-1] / predicted - 1.0)
    thermal = [i for i, ln in enumerate(said) if ln.startswith(("hotspot ", "probe "))]
    first = {t: next((i for i in thermal if t in said[i]), None) for t in ("(1,-1j)", "(1,+1j)")}
    order_ok = (
        first["(1,-1j)"] is not None
        and first["(1,+1j)"] is not None
        and first["(1,-1j)"] < first["(1,+1j)"]
    )
    g = {s: i32._thermal_failures(bs[s]["arm"]) for s in ("ccw", "cw")}
    vm = {s: bs[s]["arm"]["vol_mean"] for s in ("ccw", "cw")}
    r_vs = abs(vm["ccw"] / vm["cw"] - 1.0)
    vlines = _vm_lines(said)
    v_acc = {}
    for s, st in SENSE_TUPLE.items():
        ln = next((x for x in vlines if x.startswith(f"{_VM_HEAD}{st}:")), None)
        want = float(f"{bs[s]['scale'] * vm[s]:.9e}")
        v_acc[s] = None if ln is None else abs(_vm_acc_printed(ln) - want)
    return {
        "r_p": r_p,
        "r_h": r_h,
        "applied": applied,
        "r_s": r_s,
        "order_ok": order_ok,
        "first": first,
        "g": g,
        "vm": vm,
        "r_vs": r_vs,
        "v_acc": v_acc,
        "vlines": vlines,
        "B": {s: bs[s]["arm"]["B"] for s in ("ccw", "cw")},
    }


def _r2_asserts(r):
    assert r["r_p"] <= B_BAND, r["r_p"]
    assert r["r_h"] <= B_BAND, r["r_h"]
    assert r["applied"] == 0.0, r["applied"]
    assert r["r_s"] <= B_BAND, r["r_s"]
    assert r["order_ok"], r["first"]
    assert not r["g"]["ccw"] and not r["g"]["cw"], r["g"]  # (B), (E), (E'), (S2) in band
    assert all(v == 0.0 for v in r["v_acc"].values()), r["v_acc"]  # (V-acc)
    assert r["r_vs"] <= B_BAND, r["r_vs"]  # (V-s)


@_coding
def test_p87_kdeg_p61_knob_parsing(monkeypatch):
    """(K-deg-P61): unset / '' / 2 accepted; 1 / x / 3 refused; C52's 3 T knob refused."""
    monkeypatch.delenv(c6.C6_DEGREE_ENV, raising=False)
    monkeypatch.delenv(C21_DEGREE_ENV_3T, raising=False)
    got = {"(unset)": _p61_degree("1")}
    for raw in ("", "2"):
        monkeypatch.setenv(c6.C6_DEGREE_ENV, raw)
        got[repr(raw)] = _p61_degree("1")
    refused = []
    for bad in ("1", "x", "3"):
        monkeypatch.setenv(c6.C6_DEGREE_ENV, bad)
        with pytest.raises(ValueError, match="only unset"):
            _p61_degree("1")
        refused.append(bad)
    monkeypatch.delenv(c6.C6_DEGREE_ENV)
    monkeypatch.setenv(C21_DEGREE_ENV_3T, "2")
    with pytest.raises(ValueError, match=C21_DEGREE_ENV_3T):
        _p61_degree("1")
    print(
        f"\n[GEO-38 P87 K-deg-P61] {c6.C6_DEGREE_ENV}: {got}; refused {refused}; "
        f"{C21_DEGREE_ENV_3T}=2 refused (ASSERTED)",
        flush=True,
    )
    assert got == {"(unset)": None, "''": None, "'2'": 2}, got


@_coding
def test_p87_k69_entry69_environment_shape(monkeypatch):
    """(K-69): entry 69's environment accepted alone; the S8 degree knob added refused."""
    from tests.validation.test_geo36_rod_f_clinical_drive import (
        S8_DRIVE_C_PF_ENV,
        S8_DRIVE_DEGREE_ENV,
    )

    for k in (
        HEAT_ENV,
        HEAT_SENSE_ENV,
        S8_DRIVE_C_PF_ENV,
        S8_DRIVE_DEGREE_ENV,
        c6.C6_DEGREE_ENV,
        C21_DEGREE_ENV_3T,
    ) + FOREIGN_ENVS:
        monkeypatch.delenv(k, raising=False)
    monkeypatch.setenv(HEAT_ENV, "1")
    monkeypatch.setenv(HEAT_SENSE_ENV, "both")
    monkeypatch.setenv(c6.C6_DEGREE_ENV, "2")
    monkeypatch.setenv(S8_DRIVE_C_PF_ENV, "185.7629")
    got = (_knobs(), _sense_knob(), _p61_degree("1"))
    monkeypatch.setenv(S8_DRIVE_DEGREE_ENV, "2")
    with pytest.raises(ValueError, match="never reads"):
        _knobs()
    with pytest.raises(ValueError, match="refused"):
        i32._refuse_degree()
    print(
        f"\n[GEO-38 P87 K-69] {HEAT_ENV}=1 {HEAT_SENSE_ENV}=both {c6.C6_DEGREE_ENV}=2 "
        f"{S8_DRIVE_C_PF_ENV}=185.7629: {got}; + {S8_DRIVE_DEGREE_ENV}=2 refused (foreign guard and "
        f"_refuse_degree) (ASSERTED)",
        flush=True,
    )
    assert got == (("1", 185.7629), ("ccw", "cw"), 2), got


@complex_only
@_coding
def test_p87_k10r_d2_run_heat_box_on_degree2_stubs(monkeypatch, tmp_path):
    """(K10-R-d2) + (V12): `_run_heat_box(…, degree=2)` on `_k6_fixture_d2` with C34's
    recorders — (Deg) exact, (Deg-arith), the ``-d2`` JSON, (B) and (R2-*) at
    ``B_BAND``, (V-acc), (V-s); then the same call at degree 1 on `c6._k6_fixture`
    and the degree-1 → degree-2 identity at ``B_BAND`` (P_gel, history, volume mean)."""
    import json

    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo38_f_clinical_esr_root as c4

    comm = MPI.COMM_WORLD
    fx2 = _k6_fixture_d2(comm)
    fx1 = c6._k6_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    d2_dir, d1_dir = os.path.join(root, "d2"), os.path.join(root, "d1")
    el_deg = (_element_degree(fx2[2]), _element_degree(fx2[3]), _element_degree(fx1[2]))
    with monkeypatch.context() as m:
        m.delenv(c6.C6_DEGREE_ENV, raising=False)
        m.delenv(c4.C4_DEGREE_ENV, raising=False)
        out2, lines, said2, env_seen, degrees = _stub_heat_run(
            comm, m, fx2, d2_dir, degree=2, senses=("ccw", "cw"), recorders=True
        )
        env_after = os.environ.get(c4.C4_DEGREE_ENV)
    with monkeypatch.context() as m:
        m.delenv(c6.C6_DEGREE_ENV, raising=False)
        out1, _l1, said1, _e, _d = _stub_heat_run(
            comm, m, fx1, d1_dir, degree=None, senses=("ccw", "cw"), recorders=False
        )
    rec, d1_name = None, None
    if comm.rank == 0:
        with open(os.path.join(d2_dir, f"c6-node-{c6.K8_NODE_PF:.4f}-d2.json")) as fh:
            rec = json.load(fh)
        d1_name = os.path.exists(os.path.join(d2_dir, f"c6-node-{c6.K8_NODE_PF:.4f}.json"))
    rec, d1_name = comm.bcast((rec, d1_name), root=0)
    u2 = i35.F_CLINICAL_ROD_UNKNOWNS_RECORD[1]
    deg_want = f"(Deg) degree 2 global unknowns {u2} vs record {u2} (EXACT, ASSERTED)"
    deg_lines = [ln for ln in lines if ln.startswith("(Deg)")]
    arith = 2 * ROD_E_15T + 2 * ROD_F_15T
    r2 = _r2_checks(out2, said2)
    # (V12): the degree-1 → degree-2 identity on a field both spaces hold.
    bs1, bs2 = out1["by_sense"], out2["by_sense"]
    v12 = {}
    for s in ("ccw", "cw"):
        a1, a2 = bs1[s]["arm"], bs2[s]["arm"]
        v12[f"P_gel {s}"] = abs(a2["P"] / a1["P"] - 1.0)
        v12[f"hist {s}"] = float(
            np.max(np.abs(a2["hist"] - a1["hist"])) / np.max(np.abs(a1["hist"]))
        )
        v12[f"vol_mean {s}"] = abs(a2["vol_mean"] / a1["vol_mean"] - 1.0)
    if comm.rank == 0:
        print(
            f"\n[GEO-38 P87 K10-R-d2] element degree V.ufl_element().degree: d2 fields {el_deg[:2]}, d1 "
            f"field {el_deg[2]}; C4 variable at the fixture call {env_seen}, after {env_after!r}; solver "
            f"degrees {degrees}; returned degree {out2['degree']}"
        )
        print(
            f"[GEO-38 P87 K10-R-d2] {deg_lines[0] if deg_lines else '(Deg) line MISSING'}; (Deg-arith) "
            f"2*{ROD_E_15T} + 2*{ROD_F_15T} = {arith} vs record {u2} (ASSERTED exact)"
        )
        print(
            f"[GEO-38 P87 K10-R-d2] JSON -d2 keys == K8_JSON_KEYS {set(rec) == c6.K8_JSON_KEYS}; degree-1 "
            f"name written {d1_name}; (B) ccw {r2['B']['ccw']:.3e} cw {r2['B']['cw']:.3e} (band "
            f"{B_BAND:.0e}); failures {r2['g']}; (R2-p) {r2['r_p']:.3e}, (R2-h) {r2['r_h']:.3e}, (R2-s) "
            f"applied {r2['applied']} miss {r2['r_s']:.3e}, (R2-o) {r2['order_ok']}"
        )
        print(
            f"[GEO-38 P87 K10-R-d2] (V-acc) |printed - scale*vol_mean| {r2['v_acc']}; (V-s) vol_mean ccw "
            f"{r2['vm']['ccw']:.12e} / cw {r2['vm']['cw']:.12e} K per W: miss {r2['r_vs']:.3e} (band "
            f"{B_BAND:.0e}) (ASSERTED)"
        )
        print(
            f"[GEO-38 P87 V12] degree 2 vs degree 1 on K6's field: {v12} (band {B_BAND:.0e}, ASSERTED)",
            flush=True,
        )
    assert el_deg == (2, 2, 1), el_deg
    assert env_seen == ["2"], env_seen
    assert env_after is None, env_after
    assert degrees == [2], degrees
    assert out2["degree"] == 2 and out1["degree"] == 1
    assert deg_lines == [deg_want], deg_lines
    assert arith == 8_327_394 == u2, (arith, u2)
    assert set(rec) == c6.K8_JSON_KEYS, sorted(rec)
    assert not d1_name
    _r2_asserts(r2)
    assert len(r2["vlines"]) == 2, r2["vlines"]
    assert max(v12.values()) <= B_BAND, v12
    assert len(_vm_lines(said1)) == 2


@complex_only
@_coding
def test_p87_k10vm_d1_volume_mean_line(monkeypatch, tmp_path):
    """(K10-vm-d1): the degree-1 K10-R path prints the volume-mean line: present,
    finite, and the stored fraction in (0, 1)."""
    comm = MPI.COMM_WORLD
    fx = c6._k6_fixture(comm)
    root = comm.bcast(str(tmp_path), root=0)
    with monkeypatch.context() as m:
        m.delenv(c6.C6_DEGREE_ENV, raising=False)
        out, _l, said, _e, _d = _stub_heat_run(
            comm, m, fx, root, degree=None, senses=("cw",), recorders=False
        )
    arm = out["arm"]
    vl = _vm_lines(said)
    frac = arm["stored_frac"]
    if comm.rank == 0:
        print(
            f"\n[GEO-38 P87 K10-vm-d1] {vl}; vol_mean {arm['vol_mean']:.12e} K per W, V_sub "
            f"{arm['v_sub']:.12e} m^3, stored fraction {frac:.9e} (in (0, 1), ASSERTED)",
            flush=True,
        )
    assert len(vl) == 1 and vl[0].startswith(f"{_VM_HEAD}(1,+1j):"), vl
    assert np.isfinite(arm["vol_mean"]) and arm["vol_mean"] > 0.0
    assert np.isfinite(_vm_acc_printed(vl[0]))
    assert 0.0 < frac < 1.0, frac


@complex_only
@_coding
def test_p87_negative_controls_d2c_and_perturbed_source(monkeypatch, tmp_path):
    """(d2-c-P61): a recorder reporting degree 1 under ``degree=2`` raises at C6a's
    degree assert before any thermal arm (no ``[GEO-36 I32]`` line).  (N-B-d2):
    `_thermal_half` on `_k6_fixture_d2`'s field with q × (1 + ``PERTURBED_SOURCE``)
    gives |B / PERTURBED_SOURCE − 1| ≤ 1e-9 with "(B)" among the failures."""
    from tests.validation import test_geo36_rod_f_clinical as i35
    from tests.validation import test_geo38_f_clinical_esr_root as c4

    comm = MPI.COMM_WORLD
    fx2 = _k6_fixture_d2(comm)
    root = comm.bcast(str(tmp_path), root=0)
    i32_said = []
    with monkeypatch.context() as m:
        m.delenv(c6.C6_DEGREE_ENV, raising=False)
        m.delenv(c4.C4_DEGREE_ENV, raising=False)
        this = sys.modules[_run_heat_box.__module__]
        c6._k8_stubs(m, fx2, root)
        m.setattr(i32, "WALL_SELECT_M", -1.0)
        m.setattr(this, "_probes", lambda: i32.STUB_PROBES)
        _env_c, degrees_c = c6._c34_recorders(m, degree=1)
        real_i32_say = i32._say

        def _i32_rec(c, msg):
            i32_said.append(msg)
            real_i32_say(c, msg)

        m.setattr(i32, "_say", _i32_rec)
        with pytest.raises(AssertionError) as d2c:
            this._run_heat_box(comm, m, c6.K8_NODE_PF, senses=("ccw", "cw"), degree=2)

    msh, cell_tags, e2, _e2b, sigma_f, _c, _half = fx2
    merged = _merge_rod_into_phantom(msh, cell_tags)
    hot = (0.0, 0.0, 0.1)
    with monkeypatch.context() as m:
        m.setattr(i32, "WALL_SELECT_M", -1.0)
        real_source = i32._heat_source

        def _perturbed(*a):
            q = real_source(*a)
            q.x.array[:] *= 1.0 + PERTURBED_SOURCE
            return q

        m.setattr(i32, "_heat_source", _perturbed)
        neg = _thermal_half(
            "P87 NEG d2",
            msh,
            merged,
            e2,
            sigma_f,
            i35._rod_centre(),
            {"a": hot},
            i32.STUB_PROBES,
            hot,
        )
    neg_fails = i32._thermal_failures(neg)
    r_neg = abs(neg["B"] / PERTURBED_SOURCE - 1.0)
    if comm.rank == 0:
        print(
            f"\n[GEO-38 P87 d2-c] degree 1 under degree=2 raised {str(d2c.value).splitlines()[0]}; "
            f"solvers {degrees_c}; [GEO-36 I32] lines recorded {len(i32_said)} (ASSERTED 0)"
        )
        print(
            f"[GEO-38 P87 N-B-d2] q(1 + {PERTURBED_SOURCE:g}) on the degree-2 field: (B) {neg['B']:.9e}, "
            f"|B/{PERTURBED_SOURCE:g} - 1| {r_neg:.3e} (<= {PERTURBED_B_RTOL:g}); {neg['B'] / B_BAND:.4e} x "
            f"B_BAND; failures {neg_fails} (ASSERTED)",
            flush=True,
        )
    assert str(d2c.value).startswith("(1, 2)"), str(d2c.value)
    assert degrees_c == [], degrees_c
    assert i32_said == [], i32_said
    assert r_neg <= PERTURBED_B_RTOL, neg["B"]
    assert any("(B)" in f for f in neg_fails), neg_fails
