"""`IMP-3` step 3 (§9 item 235, chain I9): ours at degree 1 against Palace order 1.

Same two exact-rectangle builds the probe ``scripts/probes/imp3_step3_palace_rod.py``
exported (``imp3x_off`` / ``imp3x_rod``), swept through
``build_four_port_sweep(..., built=..., build_only=True)`` +
``run_n_port_sparameter_sweep(..., degree=1)``; a third sweep leaves the rod
wall natural (``replace(problem, facet_tags=outer, pec_facet_tags=(499,))``,
`IMP-3` C1's construction, all four drives).  κ off.

Asserted: (F) N_i areas and cells equal ``imp3x_meta.json``; (i)/(ii)
reciprocity / passivity on all three sweeps; (E) rod PEC wall max |E| <=
``PEC_WALL_MAX_E``; (Kw) natural wall P1 max |E| > ``PEC_WALL_MAX_E``; (M)
σ_y pairs on the rod Z; (S) without-rod separation.  Printed, never gated:
(T) the matched-order tie per arm, ΔS(rod) and ΔS_wall ours vs Palace, the
negative controls (C1) / (C2), and the ANS-4 record free tie.
"""

from __future__ import annotations

import gc
import importlib.util
import json
from dataclasses import replace
from pathlib import Path

import dolfinx
import numpy as np
from mpi4py import MPI

from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheets import SHEET_IFACE
from tests.validation.test_birdcage_power_identity import STEP2D_FREQUENCY_HZ
from tests.validation.test_imp3_rod_birdcage import (
    PEC_WALL_MAX_E,
    PORT_IDS,
    _matrix_quantities,
    _outer_box_tags,
    _wall_dofs,
    _wall_max_e,
)
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
from tests.validation.test_port_birdcage_four_port import (
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    build_four_port_sweep,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND
from tests.validation.test_th14_birdcage_copper import OUTER_BOX_TAG

REPO = Path(__file__).resolve().parents[2]
NOMINAL_STRIP_AREA = 5.6e-5


def _probe():
    path = REPO / "scripts" / "probes" / "imp3_step3_palace_rod.py"
    spec = importlib.util.spec_from_file_location("imp3_step3_palace_rod", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _sheet_areas(msh, facet_tags, comm):
    import ufl
    from dolfinx import fem

    ds = ufl.Measure("dS", domain=msh, subdomain_data=facet_tags)
    out = {}
    for i in range(1, 5):
        form = fem.form(fem.Constant(msh, dolfinx.default_scalar_type(1.0)) * ds(SHEET_IFACE + i))
        out[i] = comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM).real
    return out


def _sweep(case, problem=None, keep=False):
    return run_n_port_sparameter_sweep(
        problem if problem is not None else case["problem"],
        case["port_defs"],
        lumped_sheet_ports=case["specs"],
        lumped_sheet_facet_tags=case["facet_tags"],
        keep_fields=keep,
        degree=1,
    )


def _ids(label, mq):
    _say(
        f"[IMP-3 3 W6] {label}: recip {mq['recip']:.3e} sigma_max {float(mq['sigma'][0]):.12f} "
        f"colpow {np.array2string(mq['colpow'], precision=12)} Z spreads "
        + "  ".join(f"{k} {v:.3e}" for k, v in mq["spreads"].items())
        + f" separation {mq['separation']:.2f}x; σ_y Z pairs "
        + ", ".join(f"{v:.3e}" for v in mq["mirror"])
    )
    assert mq["recip"] <= RECIPROCITY_BAND, (label, mq["recip"])
    assert float(mq["sigma"][0]) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (label, mq["sigma"])
    assert np.all(mq["colpow"] <= 1.0), (label, mq["colpow"])


@complex_only
def test_imp3_palace_twin_degree1():
    comm = MPI.COMM_WORLD
    probe = _probe()
    meta = json.loads((probe.OUT_DIR / "imp3x_meta.json").read_text())
    ours = {}
    wall = {}
    for arm in probe.ARMS:
        built = probe.build_arm(arm)
        case = build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True)
        msh = case["mesh"]
        ncells = int(case["cells"])
        areas = _sheet_areas(msh, case["facet_tags"], comm)
        _say(
            f"\n[IMP-3 3 W6] arm {arm}: {ncells} cells (meta {meta[arm]['cells']}); N_i areas "
            + ", ".join(f"P{i} {a:.12e}" for i, a in areas.items())
        )
        # (F)
        assert ncells == meta[arm]["cells"], (arm, ncells, meta[arm]["cells"])
        for i, a in areas.items():
            assert abs(a - NOMINAL_STRIP_AREA) <= STRIP_BAND, (arm, i, a)
            assert abs(a / meta[arm]["ports"][f"P{i}"]["h"] - 7.0e-3) <= 1e-9, (arm, i)
        if arm == "off":
            res = _sweep(case)
            assert tuple(res.port_ids) == PORT_IDS
            mq = _matrix_quantities(res)
            _ids("off", mq)
            assert mq["separation"] >= POOLED_SEPARATION_FLOOR, mq["separation"]  # (S)
            ours["imp3x_off"] = mq["s"]
            del res, case, built, msh
            gc.collect()
            continue
        res = _sweep(case, keep=True)
        mq = _matrix_quantities(res)
        _ids("rod PEC", mq)
        for v in mq["mirror"]:  # (M)
            assert v <= ADJACENT_SPREAD_BAND, mq["mirror"]
        v_space = res.fields["P1"].e_complex.function_space
        owned, n_wall = _wall_dofs(msh, built[2], v_space)
        assert n_wall > 0
        wall["pec"] = {p: _wall_max_e(res.fields[p].e_complex, owned, comm) for p in PORT_IDS}
        ours["imp3x_rod"] = mq["s"]
        del res
        gc.collect()
        outer, census = _outer_box_tags(msh, built[2], comm)
        _say(f"[IMP-3 3 W6] natural-wall census {census}; wall dofs {n_wall}")
        assert census["tag_collision"] == 0, census
        assert census["exterior"] == census["outer"] + census["wall_exterior"], census
        problem = replace(case["problem"], facet_tags=outer, pec_facet_tags=(OUTER_BOX_TAG,))
        res = _sweep(case, problem=problem, keep=True)
        mq = _matrix_quantities(res)
        _ids("rod natural", mq)
        v_space = res.fields["P1"].e_complex.function_space
        owned, _ = _wall_dofs(msh, built[2], v_space)
        wall["nat"] = {p: _wall_max_e(res.fields[p].e_complex, owned, comm) for p in PORT_IDS}
        ours["imp3x_rodnat"] = mq["s"]
        del res
    _say(f"[IMP-3 3 W6] (E) rod PEC wall max |E| {wall['pec']} <= {PEC_WALL_MAX_E}")
    _say(f"[IMP-3 3 W6] (Kw) natural wall max |E| {wall['nat']} (P1 must exceed {PEC_WALL_MAX_E})")
    for p, v in wall["pec"].items():
        assert v <= PEC_WALL_MAX_E, (p, v)
    assert wall["nat"]["P1"] > PEC_WALL_MAX_E, wall["nat"]

    if comm.rank == 0:
        with open(probe.OUT_DIR / "ours_imp3x_p1.json", "w") as fh:
            json.dump(
                {
                    k: [[[c.real, c.imag] for c in row] for row in v.tolist()]
                    for k, v in ours.items()
                },
                fh,
            )
    pal = probe.read_palace()
    line = probe.line
    for name in probe.RUNS:
        _say(
            f"[IMP-3 3 W6] (T) {name} ours vs Palace per class {line(probe.rel_class_gap(ours[name], pal[name]))} (predicted <= 1e-6)"
        )
    d_o = probe.class_delta(ours["imp3x_rod"], ours["imp3x_off"])
    d_p = probe.class_delta(pal["imp3x_rod"], pal["imp3x_off"])
    w_o = probe.class_delta(ours["imp3x_rod"], ours["imp3x_rodnat"])
    w_p = probe.class_delta(pal["imp3x_rod"], pal["imp3x_rodnat"])
    mo = probe.class_means(ours["imp3x_off"])
    _say(
        f"[IMP-3 3 W6] ours |Δ_k(rod)|/|mean_k| {line({k: abs(v) / abs(mo[k]) for k, v in d_o.items()})} (predicted [1e-5, 1e-3])"
    )
    _say(
        f"[IMP-3 3 W6] ours |Δ_k,wall|/|mean_k| {line({k: abs(v) / abs(mo[k]) for k, v in w_o.items()})} (predicted [1e-6, 1e-3])"
    )
    _say(
        f"[IMP-3 3 W6] ΔS(rod) ours vs Palace per class {line(probe.delta_vs(d_o, d_p))} (predicted <= 1e-4)"
    )
    ds_o = ours["imp3x_rod"] - ours["imp3x_off"]
    ds_p = pal["imp3x_rod"] - pal["imp3x_off"]
    _say(
        f"[IMP-3 3 W6] ΔS(rod) entrywise max|ΔS_ours - ΔS_Pal|/max|ΔS_Pal| "
        f"{np.max(np.abs(ds_o - ds_p)) / np.max(np.abs(ds_p)):.3e}"
    )
    dw_o = ours["imp3x_rod"] - ours["imp3x_rodnat"]
    dw_p = pal["imp3x_rod"] - pal["imp3x_rodnat"]
    _say(
        f"[IMP-3 3 W6] ΔS_wall ours vs Palace per class {line(probe.delta_vs(w_o, w_p))}; entrywise "
        f"{np.max(np.abs(dw_o - dw_p)) / np.max(np.abs(dw_p)):.3e} (predicted <= 1e-3)"
    )
    _say(
        f"[IMP-3 3 W6] (C1) ours rod PEC vs Palace imp3x_rodnat {line(probe.rel_class_gap(ours['imp3x_rod'], pal['imp3x_rodnat']))} (predicted [1e-6, 1e-3])"
    )
    _say(
        f"[IMP-3 3 W6] (C2) ours off vs Palace imp3x_rod {line(probe.rel_class_gap(ours['imp3x_off'], pal['imp3x_rod']))} (predicted [1e-5, 1e-3])"
    )
    try:
        from tests.validation.test_ans4_exact_strip_records import (
            ANS4_EXACT_STRIP_S_P1_RECORDS,
        )

        rec = np.asarray(ANS4_EXACT_STRIP_S_P1_RECORDS["64"], dtype=np.complex128)
        _say(
            f"[IMP-3 3 W6] free tie ours off vs ANS4_EXACT_STRIP_S_P1_RECORDS['64'] "
            f"{line(probe.rel_class_gap(ours['imp3x_off'], rec))} (predicted <= 1e-6)"
        )
    except Exception as exc:  # a reading about the builders, never a stop
        _say(f"[IMP-3 3 W6] free tie not measured: {exc!r}")
