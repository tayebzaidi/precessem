"""Example (`EX-89`): why a lead's tip response cannot be read off the lead-free field.

One mesh (the `WF-9` step 1 fixture: the F-small untuned 4-leg birdcage with 50 Ohm
lumped sheets, a gel phantom and a sigma = 800 S/m capped rod, r = 1 mm, tip to tip
32 mm, centre (15, 0, 0) mm, kept as cell tag 403) at 64 MHz, degree 1.

* **One adjoint solve gives the TF.** A unit tip load (mean E_z in a 3 mm ball at the
  top tip) solved *with the lead present* gives w; TF(z) = -j w mu0 dsigma w_z along
  the lead. Folded with each port drive's lead-free field it predicts the with-lead
  tip response (T3), and as a 1-D path integral within the weekly's 5 % band (PF).
* **A drive the adjoint never saw.** The clinical quadrature drive (all four ports at
  once, phases from the fixture's own azimuths) is solved fresh with the lead; the
  same w folded with the superposed lead-free field predicts its tip response (X3).
* **Why the lead-free field is not enough.** The Born fold (adjoint solved without the
  lead) overstates the response by about 14x; the fold of w with the with-lead field
  (C2) captures only about 6.5 %.

**Import, never restate.** ``_run_step1_identities`` is the body of the gate's own
``test_wf9_step1_identities`` lifted to module level (rule (a), additive: one trailing
``return``); every band and record comes from the gate module and what it imports.
The example re-asserts on the returned case and adds (X0)-(X3) and the C4 control.
Round-off readings ((A1), (i'), (H), (T3), (X1)-(X3)) vary run to run and are asserted
at their bands only.

**The insulated lead (`WF-9` step 3, §9 item 332).** After the step 1 section the same
adjoint-fold identities run on the insulated lead (rod tag 403 inside an eps_r 3 shaft
sleeve, tag 404, caps bare) on `GEO-36` step 5's coated build (190 747 cells), with the
bare arm re-solved on that mesh for comparison: ``_run_step3_identities`` is the body of
``test_wf9_step3_identities`` lifted to module level (rule (a), additive: a ``shaft_tip``
dict and one trailing ``return``), called once. The example adds no form, solve or field
write on that mesh; every coated number comes from the returned dict and is re-asserted
against ``20261008T142909Z_WF-9-step3-W1.log``.

Scope: degree 1, 64 MHz, untuned solid coil, 50 Ohm sheets, one sigma = 800 S/m rod,
a linear tip functional (mean E_z), four port drives and one quadrature drive. Not
claimed: a measured / ISO 10974 TF, heating or SAR, a PEC or real lead (the insulated
lead only as `WF-9` step 3's record on one mesh), F-human, a tuned coil, a B1+
polarisation statement, any S / Z figure, the path-fold band (the weekly's),
h-convergence.

Run it through the example runner (the ``imp:`` group sources the complex build)::

    ./run_examples.sh -e imp:7 -n 2 -t 540

Output: ``examples/implants/paraview_output/implants_07_lead_tf_in_the_coil_combined.xdmf``.
"""

from __future__ import annotations

import dataclasses
import gc
import sys
import time
from pathlib import Path

import dolfinx
import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem, la

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.core import TimeHarmonicSolver  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.lumped import (  # noqa: E402
    lumped_port_bilinear_term,
    lumped_port_linear_term,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_port_sheets import (  # noqa: E402
    PORT_LOWER,
    PORT_UPPER,
    SHEET_IFACE,
)
from tests.validation.test_birdcage_b1_quadrature import (  # noqa: E402
    QUADRATURE_STEP_DEG,
    _port_index,
    quadrature_phase_weights,
)
from tests.validation.test_wf9_lead_in_coil_phantom import (  # noqa: E402
    BALL_RADIUS,
    CELL_COUNT_BAND,
    DISCRETE_IDENTITY_RTOL,
    FOLD_NS,
    KEPT_ROD_CELL_RECORD,
    LOAD,
    N_SLABS,
    OMEGA,
    PATH_FOLD_BAND,
    PORT_IDS,
    RECIPROCITY_BAND,
    SLAB_DZ,
    SLAB_Z0,
    SYMMETRY_RTOL,
    TIPS,
    _build_case,
    _cell_values_with_ghosts,
    _free,
    _run_step1_identities,
    _scalar,
    _sheets,
)
from tests.validation import test_wf9_insulated_lead as _step3  # noqa: E402
from tests.validation.test_wf9_insulated_lead import (  # noqa: E402
    COATED_CELL_RECORD,
    R1_RTOL,
    _build_case as _build_coated_case,
    _run_step3_identities,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "implants_07_lead_tf_in_the_coil"

# Records from 20261004T201655Z_WF-9-step1-W1.log (printed comparators; digit checks).
W1_LOG = "20261004T201655Z_WF-9-step1-W1.log"
REC_CELLS = (15000, 7605)  # lead / tip-ball owned counts (W1:1889)
REC = {  # q_without, q_with (W1:1900, 1913, 1926, 1939), PF (:1903...), BF % (:1904...)
    "P1": (-1.2753864198e-01 - 1.8713433801e-01j, -3.1505058232e-01 - 4.2222761708e-01j, 1.441865e-04, 1399.153245),
    "P2": (-1.9768599288e-02 - 2.2756234126e-02j, -4.8687195593e-02 - 5.2392004656e-02j, 1.473321e-03, 1434.958086),
    "P3": (4.1543136576e-02 + 4.8934145739e-02j, 1.0469981634e-01 + 1.1224765880e-01j, 1.509399e-04, 1388.529790),
    "P4": (-2.1013373374e-02 - 2.4680695145e-02j, -5.0248941990e-02 - 5.4629859795e-02j, 1.022076e-03, 1344.164395),
}
REC_DIGIT_RTOL = 1.0e-8  # q records are printed to 11 significant digits
REC_PRINT_RTOL = 1.0e-5  # PF / BF records are printed to 7 digits
# Records from 20261008T142909Z_WF-9-step3-W1.log (`-n 2`; the coated-lead rows, item 332).
# Compared against W1 only: W2 (`-n 4`) differs in the last printed digit of coated P1 q_with.
W3_LOG = "20261008T142909Z_WF-9-step3-W1.log"
REC_COATED_OWNED = (15035, 18495)  # tags 403 / 404 owned counts (W1:1948)
REC_SHAFT_TIP = {"bare": 8.670205, "coated": 3.043634}  # W1:2098, W1:2099 (.6f)
# q at .10e, ratio at .6f, PF at .6e, BF and C3c at .6f %, coat share at .6f (one source: W1).
REC_COATED = {
    "bare": {
        # W1:1966 q, :1967 ratio, :1970 PF, :1971 BF, :1972 C3c + share
        "P1": dict(q_without=(-0.12507997183-0.18319294041j), q_with=(-0.33454696024-0.44669665165j), dq=(-0.20946698841-0.26350371124j),
                   ratio=2.515925, pf=3.089708e-04, bf=1574.706958, c3=0.000000, share=0.000000),
        # W1:1979 q, :1980 ratio, :1983 PF, :1984 BF, :1985 C3c + share
        "P2": dict(q_without=(-0.020336247385-0.024213736482j), q_with=(-0.052514024761-0.057222548774j), dq=(-0.032177777376-0.033008812292j),
                   ratio=2.456204, pf=1.536992e-03, bf=1666.679940, c3=0.000000, share=0.000000),
        # W1:1992 q, :1993 ratio, :1996 PF, :1997 BF, :1998 C3c + share
        "P3": dict(q_without=(0.040134706451+0.0470540373j), q_with=(0.11056311762+0.11772411862j), dq=(0.070428411172+0.070670081318j),
                   ratio=2.611388, pf=8.678799e-04, bf=1578.033775, c3=0.000000, share=0.000000),
        # W1:2005 q, :2006 ratio, :2009 PF, :2010 BF, :2011 C3c + share
        "P4": dict(q_without=(-0.020981897339-0.024487533345j), q_with=(-0.054449191033-0.059198460039j), dq=(-0.033467293694-0.034710926695j),
                   ratio=2.494209, pf=1.718626e-03, bf=1493.126591, c3=0.000000, share=0.000000),
    },
    "coated": {
        # W1:2018 q, :2019 ratio, :2021 PF, :2022 BF, :2023 C3c + share
        "P1": dict(q_without=(-0.12507997183-0.18319294041j), q_with=(-0.23689484167-0.34898329039j), dq=(-0.11181486984-0.16579034998j),
                   ratio=1.901495, pf=1.224651e-03, bf=2712.936245, c3=1.215826, share=0.012158),
        # W1:2030 q, :2031 ratio, :2033 PF, :2034 BF, :2035 C3c + share
        "P2": dict(q_without=(-0.020336247385-0.024213736482j), q_with=(-0.037882011976-0.045203845295j), dq=(-0.01754576459-0.020990108813j),
                   ratio=1.865179, pf=2.544819e-03, bf=2870.440874, c3=1.437765, share=0.014378),
        # W1:2042 q, :2043 ratio, :2045 PF, :2046 BF, :2047 C3c + share
        "P3": dict(q_without=(0.040134706451+0.0470540373j), q_with=(0.078390444443+0.092040439773j), dq=(0.038255737992+0.044986402473j),
                   ratio=1.954848, pf=2.867874e-03, bf=2728.955682, c3=1.368714, share=0.013687),
        # W1:2054 q, :2055 ratio, :2057 PF, :2058 BF, :2059 C3c + share
        "P4": dict(q_without=(-0.020981897339-0.024487533345j), q_with=(-0.039198708404-0.046563466446j), dq=(-0.018216811065-0.022075933101j),
                   ratio=1.887490, pf=1.875099e-03, bf=2577.939290, c3=0.864207, share=0.008642),
    },
}

PRED_DQ_CCW = 0.390
PRED_RATIO_CCW = (2.0, 2.7)
PRED_PF_CCW = (5.0e-5, 1.0e-2)
PRED_BF_CCW = (100.0, 3000.0)  # percent
PRED_CW_CCW = (1.0e-4, 1.0e-2)
PRED_INLEAD = (0.03, 0.2)
PRED_C4 = (10.0, 40.0)  # percent
PRED_PEAK_K = 27
PRED_TF_ARG = (-7.0, -1.5)


def _band(x, lo, hi):
    return "IN" if lo <= x <= hi else "OUT"


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )

    def say(msg):
        if rank0:
            print(msg, flush=True)

    say("=" * 72)
    say("EX-89 -- why a lead's tip response cannot be read off the lead-free field")
    say("=" * 72)
    say(
        "\n[scope] degree 1, 64 MHz, untuned solid F-small 4-leg birdcage, 50 Ohm sheets, "
        "one sigma = 800 S/m capped rod, linear tip functional (mean E_z), four port "
        "drives + one quadrature drive. No measured / ISO 10974 TF, no heating or SAR, "
        "no tuned coil, no B1+ polarisation statement, no path-fold band of our own. "
        "The insulated lead appears only as WF-9 step 3's record on one mesh (no real lead)."
    )

    # ---- the fixture -----------------------------------------------------------
    t0 = time.perf_counter()
    fx = _build_case()
    msh = fx["mesh"]
    tdim = msh.topology.dim
    n_owned = fx["n_owned"]
    cells = fx["cells"]
    n_local = len(cells)
    flags = fx["flags"]
    tol = DISCRETE_IDENTITY_RTOL
    say(f"\n[fixture reproduced] _build_case in {time.perf_counter() - t0:.1f} s")

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) --------------------
    # Figure-only tags, never passed to a solve: tip-ball cells set to 404.
    tag_vals = np.asarray(_cell_values_with_ghosts(msh, fx["built"][1]), dtype=np.int32).copy()
    tag_vals[flags["ball"]] = 404
    fig_tags = dolfinx.mesh.meshtags(msh, tdim, cells, tag_vals)
    port_names = {PORT_LOWER + i: f"port P{i} box (lower)" for i in range(1, 5)}
    port_names.update({PORT_UPPER + i: f"port P{i} box (upper)" for i in range(1, 5)})
    write_setup_figure(
        msh,
        fig_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: "conductor (birdcage)",
            2: "air",
            3: "gel phantom (sigma 0.5 S/m)",
            403: "lead (kept rod cells, sigma 800 S/m)",
            404: "tip ball r = 3 mm (q readout)",
            **port_names,
        },
        hide_tags=(2,),
        translucent_tags=(3,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:7 -- lead transfer function in the coil-driven gel phantom (WF-9 step 1)",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- the gate's own body, verbatim ----------------------------------------
    t0 = time.perf_counter()
    case = _run_step1_identities(fx)
    say(f"\n[gate reproduced] _run_step1_identities in {time.perf_counter() - t0:.1f} s")
    q = case["q"]
    w = case["w"]
    w0 = case["w0"]
    e_o = case["e_without"]
    e_w = case["e_with"]
    res = case["res"]
    tf_k = case["tf_k"]
    v_k = case["v_k"]

    dx = ufl.dx(domain=msh)
    chi_lead = fx["chi_lead"]
    ds_ = fx["delta_sigma"]

    # ---- anchors ---------------------------------------------------------------
    say(f"\n[anchors] asserted at the imported bands, on the returned case (records: {W1_LOG}):")
    rel_cells = abs(fx["n_cells"] / KEPT_ROD_CELL_RECORD - 1.0)
    n_lead = fx["owned_count"](flags["lead"])
    n_ball = fx["owned_count"](flags["ball"])
    say(
        f"  (A0) cells={fx['n_cells']} vs record {KEPT_ROD_CELL_RECORD} (rel {rel_cells:.3e}, band "
        f"{CELL_COUNT_BAND:g}); lead / ball {n_lead} / {n_ball} (record {REC_CELLS[0]} / {REC_CELLS[1]})"
    )
    assert rel_cells <= CELL_COUNT_BAND, rel_cells
    assert (n_lead, n_ball) == REC_CELLS, (n_lead, n_ball)
    say(f"  (A1) ||A-A^T||/||A|| = {case['sym_ratio']:.3e} (band {SYMMETRY_RTOL:g}; record 9.757e-17)")
    assert case["sym_ratio"] <= SYMMETRY_RTOL
    rc = case["recip"]
    say(
        f"  (i') with {rc['with']:.3e}, without {rc['without']:.3e} (band {RECIPROCITY_BAND:g}; "
        "record 1.505e-15 / 1.796e-15)"
    )
    assert rc["with"] <= RECIPROCITY_BAND and rc["without"] <= RECIPROCITY_BAND
    say(f"  (H)  {case['h_dev']:.3e} (band {tol:g}; record 0.000e+00)")
    assert case["h_dev"] <= tol
    for d in PORT_IDS:
        r = res[d]
        dq = r["dq"]
        qo_rec, qw_rec, pf_rec, bf_rec = REC[d]
        t3 = abs(r["t3"] - dq) / abs(dq)
        pf = abs(r["folds"][32][0] - dq) / abs(dq)
        bf = abs(r["born"] - dq) / abs(dq)
        say(
            f"  {d}: q_without={r['q_without']:.10e} (record {qo_rec:.10e}); "
            f"q_with={r['q_with']:.10e} (record {qw_rec:.10e})\n"
            f"      (T3) {t3:.3e} (band {tol:g}); (PF) {pf:.6e} (record {pf_rec:.6e}, band "
            f"{PATH_FOLD_BAND:g}); (BF) {bf:.6%} (record {bf_rec:.6f} %, must exceed "
            f"{PATH_FOLD_BAND:.0%})"
        )
        assert abs(r["q_without"] - qo_rec) <= REC_DIGIT_RTOL * abs(qo_rec), (d, "q_without")
        assert abs(r["q_with"] - qw_rec) <= REC_DIGIT_RTOL * abs(qw_rec), (d, "q_with")
        assert abs(pf / pf_rec - 1.0) <= REC_PRINT_RTOL, (d, "PF record")
        assert abs(100.0 * bf / bf_rec - 1.0) <= REC_PRINT_RTOL, (d, "BF record")
        assert abs(r["t3"] - dq) <= tol * abs(dq), (d, "T3")
        assert abs(r["folds"][32][0] - dq) <= PATH_FOLD_BAND * abs(dq), (d, "PF")
        assert abs(r["born"] - dq) > PATH_FOLD_BAND * abs(dq), (d, "BF does not miss")

    # ---- (X0) the phase pattern ------------------------------------------------
    say("\n[example-side, asserted]")
    az = {f"P{s['tag'] - SHEET_IFACE}": float(s["azimuth_deg"]) for s in fx["sweep"]["sheets"]}
    ks = [_port_index(az[p], az["P1"], QUADRATURE_STEP_DEG) for p in PORT_IDS]
    a_ccw = quadrature_phase_weights(ks, "ccw")
    a_cw = quadrature_phase_weights(ks, "cw")
    say(
        "  (X0) azimuths [deg] " + ", ".join(f"{p} {az[p]:.3f}" for p in PORT_IDS)
        + f"; ks = {ks}; a_ccw = " + ", ".join(f"{x:.3f}" for x in a_ccw)
    )
    assert sorted(ks) == [0, 1, 2, 3], ks

    # ---- DG0 helpers -----------------------------------------------------------
    dg0 = fem.functionspace(msh, ("DG", 0))
    dof_of_cell = dg0.dofmap.list[:, 0]

    def field_fn(per_cell, name):
        f = fem.Function(dg0, name=name)
        f.x.array[:] = 0.0
        f.x.array[dof_of_cell[cells]] = np.asarray(per_cell, dtype=float)
        f.x.scatter_forward()
        return f

    def cellwise(form_expr):
        """Per-cell integral of ``form_expr(v)`` against a DG0 test fn (owned+ghost, complex)."""
        v = ufl.TestFunction(dg0)
        b = fem.assemble_vector(fem.form(form_expr(v)))
        b.scatter_reverse(la.InsertMode.add)
        b.scatter_forward()
        return np.array(b.array)[dof_of_cell[cells]]

    vol = np.real(cellwise(lambda v: ufl.conj(v) * dx))

    # ---- (X1) the TF-density map -----------------------------------------------
    tfd = cellwise(lambda v: LOAD * ds_ * chi_lead * w[2] * ufl.conj(v) * dx)
    x1 = 0.0
    for k in range(N_SLABS):
        s = comm.allreduce(complex(np.sum(tfd[:n_owned][flags["slabs"][k][:n_owned]])), op=MPI.SUM)
        x1 = max(x1, abs(s - tf_k[k]) / abs(tf_k[k]))
    say(f"  (X1) TF-density map: max_k |sum_cells - tf_k|/|tf_k| = {x1:.3e} over {N_SLABS} slabs (band {tol:g})")
    assert x1 <= tol, x1

    # ---- (X2) the fold-integrand maps ------------------------------------------
    fold_maps = {}
    for d in PORT_IDS:
        eo = e_o[d]
        vec = cellwise(lambda v, e=eo: LOAD * ds_ * chi_lead * ufl.dot(w, e) * ufl.conj(v) * dx)
        fold_maps[d] = vec
        s = comm.allreduce(complex(np.sum(vec[:n_owned])), op=MPI.SUM)
        x2 = abs(s - res[d]["t3"]) / abs(res[d]["t3"])
        say(f"  (X2) {d} fold-integrand map: |sum - t3|/|t3| = {x2:.3e} (band {tol:g})")
        assert x2 <= tol, (d, x2)

    # ---- (X3) the quadrature drive, fresh, with the lead -----------------------
    sweep = fx["sweep"]
    tags_f = sweep["facet_tags"]
    specs = sweep["specs"]
    for spec in specs:
        assert complex(spec.drive_voltage_v) == 1.0 + 0.0j, (spec.port_id, spec.drive_voltage_v)
    held_sheets = _sheets(specs)
    bil = [
        (lambda trial, test, _s=s: lumped_port_bilinear_term(
            msh, tags_f, _s, trial, test, omega_rad_per_s=OMEGA))
        for s in held_sheets
    ]
    a_of = dict(zip(PORT_IDS, a_ccw))
    lin = []
    for spec in specs:
        drv = dataclasses.replace(spec.sheet(driven=True), source_voltage_v=complex(a_of[spec.port_id]))
        lin.append(lambda v, _s=drv: lumped_port_linear_term(
            msh, tags_f, _s, v, omega_rad_per_s=OMEGA))
    t0 = time.perf_counter()
    solver = TimeHarmonicSolver(fx["p_with"], degree=1)
    e_with_ccw = solver.solve(
        current_density=None, project_source=False,
        extra_bilinear_terms=bil, extra_linear_terms=lin,
    ).e_complex
    t_x3 = time.perf_counter() - t0
    _free(solver)

    space = e_o["P1"].function_space
    e_wo_ccw = fem.Function(space, name="E_without_ccw")
    e_wo_ccw.x.array[:] = sum(a_of[d] * e_o[d].x.array for d in PORT_IDS)
    e_wo_ccw.x.scatter_forward()
    e_wo_sum = fem.Function(space, name="E_without_unphased_sum")
    e_wo_sum.x.array[:] = sum(e_o[d].x.array for d in PORT_IDS)
    e_wo_sum.x.scatter_forward()

    q_wo_ccw = sum(a_of[d] * res[d]["q_without"] for d in PORT_IDS)
    q_wo_ccw_direct = q(e_wo_ccw)
    q_with_ccw = q(e_with_ccw)
    dq_ccw = q_with_ccw - q_wo_ccw
    fold_ccw = LOAD * ds_ * _scalar(comm, chi_lead * ufl.dot(w, e_wo_ccw) * dx)
    fold_sum_t3 = sum(a_of[d] * res[d]["t3"] for d in PORT_IDS)
    x3 = abs(fold_ccw - dq_ccw) / abs(dq_ccw)
    x3b = abs(fold_ccw - fold_sum_t3) / abs(fold_ccw)
    say(
        f"  (X3) quadrature ccw, one fresh factorisation with the lead ({t_x3:.1f} s):\n"
        f"       q_without,ccw = {q_wo_ccw:.10e} (superposed; direct read {q_wo_ccw_direct:.10e}); "
        f"q_with,ccw = {q_with_ccw:.10e}\n"
        f"       dq_ccw = {dq_ccw:.10e}; fold(w, E_without,ccw) = {fold_ccw:.10e}; "
        f"|fold - dq|/|dq| = {x3:.3e} (band {tol:g}, predicted <= 1e-11)\n"
        f"       |fold - sum_d a_d t3_d|/|fold| = {x3b:.3e} (band {tol:g}; linearity wiring check)"
    )
    assert abs(q_wo_ccw_direct - q_wo_ccw) <= tol * abs(q_wo_ccw), "superposition wiring"
    assert abs(fold_ccw - dq_ccw) <= tol * abs(dq_ccw), (fold_ccw, dq_ccw, x3)
    assert abs(fold_ccw - fold_sum_t3) <= tol * abs(fold_ccw), x3b
    say("  all of (X0)-(X3) hold")

    # ---- controls --------------------------------------------------------------
    say("\n[controls] each ASSERTED wired:")
    born_ccw = sum(a_of[d] * res[d]["born"] for d in PORT_IDS)
    bf_ccw = abs(born_ccw - dq_ccw) / abs(dq_ccw)
    say(
        f"  (BF) ccw Born fold miss |B - dq|/|dq| = {bf_ccw:.6%} (must exceed {PATH_FOLD_BAND:.0%}; "
        f"predicted [{PRED_BF_CCW[0]:.0f}%, {PRED_BF_CCW[1]:.0f}%], {_band(100 * bf_ccw, *PRED_BF_CCW)}); "
        f"|B|/|dq| = {abs(born_ccw) / abs(dq_ccw):.4g}"
    )
    assert abs(born_ccw - dq_ccw) > PATH_FOLD_BAND * abs(dq_ccw), "BF ccw does not miss"
    for d in PORT_IDS:
        r = res[d]
        dq = r["dq"]
        c1v = abs(r["c1"] - dq) / abs(dq)
        c2v = abs(r["c2"] - dq) / abs(dq)
        say(
            f"  {d}: C1 Hermitian read-back {c1v:.6%} (ceiling 200%); "
            f"C2 total-field fold {c2v:.6%} (record approx. 93.5%)"
        )
        assert abs(r["c1"] - dq) > tol * abs(dq), (d, "C1 not wired")
        assert abs(r["c2"] - dq) > tol * abs(dq), (d, "C2 not wired")
    fold_unphased = sum(res[d]["t3"] for d in PORT_IDS)
    c4 = abs(fold_unphased - dq_ccw) / abs(dq_ccw)
    say(
        f"  C4 unphased drive (fold of w with sum_d E_without,d) vs the fresh ccw dq: {c4:.6%} "
        f"(predicted [{PRED_C4[0]:.0f}%, {PRED_C4[1]:.0f}%], {_band(100 * c4, *PRED_C4)}; no factor asserted)"
    )
    assert abs(fold_unphased - dq_ccw) > tol * abs(dq_ccw), "C4 not wired: the phases did not reach the solve"

    # ---- printed ---------------------------------------------------------------
    say("\n[printed, no band] ccw predictions from the per-port records:")
    ratio = abs(q_with_ccw) / abs(q_wo_ccw)
    say(
        f"  |dq_ccw| = {abs(dq_ccw):.4f} V/m per V (predicted {PRED_DQ_CCW}); "
        f"|q_with/q_without| = {ratio:.4f} (predicted [{PRED_RATIO_CCW[0]:g}, {PRED_RATIO_CCW[1]:g}], "
        f"{_band(ratio, *PRED_RATIO_CCW)})"
    )
    q32_ccw = sum(a_of[d] * res[d]["folds"][32][0] for d in PORT_IDS)
    pf_ccw = abs(q32_ccw - dq_ccw) / abs(dq_ccw)
    say(
        f"  path fold N=32, ccw: |Q32 - dq|/|dq| = {pf_ccw:.6e} (predicted [{PRED_PF_CCW[0]:g}, "
        f"{PRED_PF_CCW[1]:g}], {_band(pf_ccw, *PRED_PF_CCW)}; the 5 % band is the weekly's)"
    )
    for n_seg in FOLD_NS:
        qn_ccw = sum(a_of[d] * res[d]["folds"][n_seg][0] for d in PORT_IDS)
        per_port = "  ".join(
            f"{d} {abs(res[d]['folds'][n_seg][0] - res[d]['dq']) / abs(res[d]['dq']):.3e}" for d in PORT_IDS
        )
        say(f"  N={n_seg:2d}: ccw {abs(qn_ccw - dq_ccw) / abs(dq_ccw):.3e}  |  {per_port}")
    dq_cw = sum(a_cw[i] * res[d]["dq"] for i, d in enumerate(PORT_IDS))
    cwccw = abs(dq_cw / dq_ccw - 1.0)
    say(
        f"  |dq_cw/dq_ccw - 1| = {cwccw:.3e} (predicted [{PRED_CW_CCW[0]:g}, {PRED_CW_CCW[1]:g}], "
        f"{_band(cwccw, *PRED_CW_CCW)}): the rod lies on leg 1's azimuth in the y = 0 mirror "
        "plane (dq_P2 approx. dq_P4), so the lead cannot tell the two rotation senses apart"
    )
    ez_with = _scalar(comm, chi_lead * e_with_ccw[2] * dx)
    ez_wo = _scalar(comm, chi_lead * e_wo_ccw[2] * dx)
    inlead = abs(ez_with) / abs(ez_wo)
    say(
        f"  in-lead field |int chi_lead E_with,ccw,z| / |int chi_lead E_without,ccw,z| = {inlead:.4f} "
        f"(predicted [{PRED_INLEAD[0]:g}, {PRED_INLEAD[1]:g}], {_band(inlead, *PRED_INLEAD)}; "
        "quasi-static depolarisation 0.07 at aspect 16)"
    )
    say(f"\n  TF(z) from one adjoint solve (k, z_mid mm, |TF_k|, |TF_k|/max, arg deg; W1:1952-1983):")
    tmax = float(np.max(np.abs(tf_k)))
    args = np.degrees(np.angle(tf_k))
    for k in range(N_SLABS):
        say(
            f"  {k:2d} {(SLAB_Z0 + (k + 0.5) * SLAB_DZ) * 1e3:6.1f} {abs(tf_k[k]):.6e} "
            f"{abs(tf_k[k]) / tmax:.4f} {args[k]:9.4f}"
        )
    kpk = int(np.argmax(np.abs(tf_k)))
    say(
        f"  peak at k = {kpk} (predicted {PRED_PEAK_K}, {'IN' if kpk == PRED_PEAK_K else 'OUT'}; "
        f"|TF| {abs(tf_k[kpk]):.6e}); arg range [{args.min():.2f}, {args.max():.2f}] deg "
        f"(predicted within [{PRED_TF_ARG[0]:g}, {PRED_TF_ARG[1]:g}], "
        f"{_band(args.min(), *PRED_TF_ARG)}/{_band(args.max(), *PRED_TF_ARG)})"
    )
    ta = case["tf_agreement"]
    if ta is None:
        say("  W1/W2 TF agreement: no other-width TF file yet")
    else:
        say(f"  W1/W2 TF agreement vs n={ta[0]}: max rel = {ta[1]:.3e} (predicted <= 1e-10)")
    say(f"  gate body: dofs={case['ndofs']}, elapsed {case['elapsed']:.1f} s, summed max RSS {case['rss']:.3f} GiB")

    # ---- ParaView --------------------------------------------------------------
    mag = np.abs(tfd)
    phase = np.where(mag > 0.0, np.degrees(np.angle(tfd)), 0.0)
    safe_vol = np.where(vol > 0.0, vol, 1.0)

    def rms(fn):
        e2 = np.real(cellwise(lambda v: ufl.inner(fn, fn) * ufl.conj(v) * dx))
        return np.sqrt(np.maximum(e2 / safe_vol, 0.0))

    fold_ccw_map = cellwise(lambda v: LOAD * ds_ * chi_lead * ufl.dot(w, e_wo_ccw) * ufl.conj(v) * dx)
    born_ccw_map = cellwise(lambda v: LOAD * ds_ * chi_lead * ufl.dot(w0, e_wo_ccw) * ufl.conj(v) * dx)
    slab_idx = np.zeros(n_local)
    for k in range(N_SLABS):
        slab_idx[flags["slabs"][k]] = k + 1
    fields = {
        "tf_density_magnitude": field_fn(mag, "tf_density_magnitude"),
        "tf_density_phase_deg": field_fn(phase, "tf_density_phase_deg"),
        "adjoint_w_magnitude": field_fn(rms(w), "adjoint_w_magnitude"),
        "adjoint_w0_magnitude": field_fn(rms(w0), "adjoint_w0_magnitude"),
        "E_without_magnitude_ccw": field_fn(rms(e_wo_ccw), "E_without_magnitude_ccw"),
        "E_with_magnitude_ccw": field_fn(rms(e_with_ccw), "E_with_magnitude_ccw"),
        "born_integrand_magnitude_ccw": field_fn(np.abs(born_ccw_map), "born_integrand_magnitude_ccw"),
        "fold_integrand_magnitude_ccw": field_fn(np.abs(fold_ccw_map), "fold_integrand_magnitude_ccw"),
        "lead_slab_index": field_fn(slab_idx, "lead_slab_index"),
        "tip_ball_mask": field_fn(flags["ball"].astype(float), "tip_ball_mask"),
    }
    for d in PORT_IDS:
        fields[f"fold_integrand_magnitude_{d}"] = field_fn(
            np.abs(fold_maps[d]), f"fold_integrand_magnitude_{d}"
        )
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        fx["built"][1],
        fields,
        comm=comm,
        facet_tags=fx["built"][2],
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    say(f"\n[paraview] wrote {xdmf_path}")

    # ---- the insulated lead (WF-9 step 3), second mesh -------------------------
    # Hygiene, not a gate: free the step 1 objects before the second build.
    del (fx, msh, case, q, w, w0, e_o, e_w, eo, dx, chi_lead, sweep, tags_f, specs, held_sheets,
         bil, lin, drv, e_with_ccw, space, e_wo_ccw, e_wo_sum, dg0, field_fn, cellwise, rms,
         fields, fig_tags, solver, ds_, flags, cells, dof_of_cell)
    gc.collect()
    assert _step3.DISCRETE_IDENTITY_RTOL is DISCRETE_IDENTITY_RTOL, "one band object, two importers"
    assert _step3.PATH_FOLD_BAND is PATH_FOLD_BAND, "one band object, two importers"
    say(
        f"\n[WF-9 step3 reproduced] the gate's own body on the coated build (records: {W3_LOG}); "
        "the block below is the gate's print, verbatim. The coated C1 / C2 were OUT of step 1's "
        "predicted bands: records, printed only, no band."
    )
    t0 = time.perf_counter()
    fx_c = _build_coated_case()
    case_c = _run_step3_identities(fx_c)
    t_c = time.perf_counter() - t0
    res_c = case_c["res"]

    say(f"\n[coated lead] asserted on the returned case_c (records: {W3_LOG}):")
    rel_c = abs(fx_c["n_cells"] / COATED_CELL_RECORD - 1.0)
    owned = (
        fx_c["owned_count"](fx_c["flags"]["rod"]),
        fx_c["owned_count"](fx_c["flags"]["coat"]),
    )
    say(
        f"  (A0c) cells={fx_c['n_cells']} vs COATED_CELL_RECORD {COATED_CELL_RECORD} (rel {rel_c:.3e}, "
        f"band {CELL_COUNT_BAND:g}); owned 403 / 404 = {owned[0]} / {owned[1]} "
        f"(record {REC_COATED_OWNED[0]} / {REC_COATED_OWNED[1]})"
    )
    assert rel_c <= CELL_COUNT_BAND, rel_c
    assert owned == REC_COATED_OWNED, owned
    for arm in ("bare", "coated"):
        assert case_c["r1"][arm][0] <= R1_RTOL and case_c["r1"][arm][1] == 0.0, (arm, case_c["r1"][arm])
    assert all(v <= SYMMETRY_RTOL for v in case_c["sym"].values()), case_c["sym"]
    assert case_c["h_dev"] <= tol, case_c["h_dev"]
    say(
        f"  (R1) {R1_RTOL:g}: bare {case_c['r1']['bare'][0]:.3e}, coated {case_c['r1']['coated'][0]:.3e}; "
        "(A1) " + ", ".join(f"{k} {v:.3e}" for k, v in case_c["sym"].items())
        + f" (band {SYMMETRY_RTOL:g}); (H) {case_c['h_dev']:.3e} (band {tol:g})"
    )
    for arm in ("bare", "coated"):
        for d in PORT_IDS:
            r = res_c[arm][d]
            rec = REC_COATED[arm][d]
            dq = r["dq"]
            t3 = abs(r["t3"] - dq) / abs(dq)
            pf = abs(r["folds"][32] - dq) / abs(dq)
            bf = abs(r["born"] - dq) / abs(dq)
            c3 = abs(r["rod_only"] - dq) / abs(dq)
            share = abs(r["coat_part"]) / abs(r["t3"])
            ratio = abs(r["q_with"]) / abs(r["q_without"])
            say(
                f"  {arm} {d}: q_with={r['q_with']:.10e} (record {rec['q_with']:.10e}); "
                f"dq={dq:.10e} (record {rec['dq']:.10e})\n"
                f"      (T3c) {t3:.3e} (band {tol:g}); (PF) {pf:.6e} (record {rec['pf']:.6e}, band "
                f"{PATH_FOLD_BAND:g}); (BF) {bf:.6%} (record {rec['bf']:.6f} %, must exceed "
                f"{PATH_FOLD_BAND:.0%}); (C3c) {c3:.6%} (record {rec['c3']:.6f} %), coat share {share:.6f} "
                f"(record {rec['share']:.6f})"
            )
            for key in ("q_without", "q_with", "dq"):
                assert abs(r[key] - rec[key]) <= REC_DIGIT_RTOL * abs(rec[key]), (arm, d, key)
            assert abs(ratio - rec["ratio"]) <= 1.0e-6, (arm, d, "ratio")
            assert abs(pf / rec["pf"] - 1.0) <= REC_PRINT_RTOL, (arm, d, "PF record")
            assert abs(100.0 * bf / rec["bf"] - 1.0) <= REC_PRINT_RTOL, (arm, d, "BF record")
            assert abs(r["t3"] - dq) <= tol * abs(dq), (arm, d, "T3c")
            assert abs(r["folds"][32] - dq) <= PATH_FOLD_BAND * abs(dq), (arm, d, "PF")
            assert abs(r["born"] - dq) > PATH_FOLD_BAND * abs(dq), (arm, d, "BF does not miss")
            assert abs(r["c1"] - dq) > tol * abs(dq), (arm, d, "C1 not wired")
            assert abs(r["c2"] - dq) > tol * abs(dq), (arm, d, "C2 not wired")
            if arm == "coated":
                assert abs(r["rod_only"] - dq) > tol * abs(dq), (d, "C3c not wired")
                assert abs(100.0 * c3 / rec["c3"] - 1.0) <= REC_PRINT_RTOL, (d, "C3c record")
                assert abs(share - rec["share"]) <= 1.0e-6, (d, "coat share record")
    for arm in ("bare", "coated"):
        st = case_c["shaft_tip"][arm]
        say(f"  shaft/tip |TF| {arm} = {st:.6f} (record {REC_SHAFT_TIP[arm]:.6f})")
        assert abs(st - REC_SHAFT_TIP[arm]) <= 1.0e-6, (arm, "shaft/tip record")
    say(
        f"  coated section {t_c:.1f} s; gate body dofs={case_c['ndofs']}, elapsed {case_c['elapsed']:.1f} s, "
        f"summed max RSS {case_c['rss']:.3f} GiB (a process peak, includes the step 1 section)"
    )
    ta_c = case_c["tf_agreement"]
    if ta_c is None:
        say("  W1/W2 TF agreement (step 3): no other-width TF file yet")
    else:
        say(f"  W1/W2 TF agreement (step 3) vs n={ta_c[0]}: max rel = {ta_c[1]:.3e} (printed)")
    say("  all (A0c), (R1), (A1), (H), (T3c), (PF) hold; (BF) misses on both arms and every drive")

    say(f"\nAll anchors and controls hold. Total elapsed {time.perf_counter() - started:.1f} s.")


if __name__ == "__main__":
    main()
