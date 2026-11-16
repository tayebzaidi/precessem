"""Example (`EX-83`): a PEC rod in the coil-driven phantom -- the 1 g tip SAR with and without it.

`imp:1` / `imp:2` impose a uniform field on a metal body in gel; `mri:3` drives
the F-small birdcage into a saline phantom with no implant in it. This example
joins them: the gapped 4-port birdcage at 64 MHz (degree 1, the 1 g phantom
rung) is solved twice -- once as built, once with a 2 mm diameter, 32 mm
tip-to-tip PEC rod standing axially in the phantom on leg 1's azimuth -- and the
1 g mass-averaged SAR is read at the two rod tips on the same six drives (four
single-port, plus the quadrature ccw / cw). The ratio Q (with rod / without) is
the headline, set beside `IMP-2`'s closed-form estimate.

**Import, never restate.** The two arms, every band and every record come from
``tests/validation/test_imp3_rod_birdcage.py`` (`IMP-3`) and the modules it
imports them from. The construction is *called*, not re-implemented:
``_run_one_gram_rung`` is the body of the gate's own ``test_imp3_one_gram_rung``
lifted to module level (rule (a), additive; ``field_hook=None`` leaves every
print of the gate unchanged -- the gate was re-run from this tree in the same
slot, digit-identical). The only new code is a hook, called inside each arm
just before its solved fields are freed, that keeps DG0 maps and profile
samples (never the fields themselves).

**Anchors (asserted at the imported bands, on the returned case).**
(A0) the without-rod cell count within ``CELL_COUNT_BAND`` of
``ONE_GRAM_RUNG_CELL_RECORD`` and the rod arm under ``ROD_CELL_GUARD_W2``;
(i) reciprocity and (ii) passivity of the 4x4 on both arms; (P) the power
identity on every port drive of both arms at ``DISCRETE_IDENTITY_RTOL``;
(M) the rod-arm mirror pairs at ``ADJACENT_SPREAD_BAND``; (E) the rod wall's
tangential E, count > 0 and max <= ``PEC_WALL_MAX_E``; (S) the without-rod
separation >= ``POOLED_SEPARATION_FLOOR``. Round-off readings ((i), (P), (E))
vary run to run and are asserted at their bands only.

**Example-side checks (asserted, true by construction).** On both arms the
owned sum of the P1 DG0 phantom-loss vector, assembled from
``0.5 sigma (E . conj E) conj(v) dx(PHANTOM_CELL_TAG)``, equals the gate's
``shares["P1"]["phantom"]`` within ``DISCRETE_IDENTITY_RTOL``; every profile
point's ``ok`` mask is True.

**Negative control (asserted; backed on this fixture).** The rod-blind arm
under-reads every tip, so Q > 1 on all 6 drives x 2 tips. The size is printed
(min Q - 1 and 1 - 1/Q; the latter has a ceiling of 100 %, Q has none).

**Printed, no band.** The gate's own table (Q, eta_z, M_cf, the estimate and
Q / estimate; P2 / P4 Q read OUT of the gate's ``PRED_Q`` because their field on
the axis is transverse); |E|_rod / |E|_without for P1 and ccw along the axial
line beyond the top tip and along the midplane line.

**Scope.** The 1 g rung, degree 1, 64 MHz, the untuned solid coil with 50 Ohm
sheets, one axial rod on leg 1's azimuth. Not claimed: the accuracy of Q
(`ANS-10`), absolute or pointwise SAR, heating or compliance, C1 on this rung,
a without-rod map on the rod mesh, the tuned or copper coil, other placements.

Run it through the example runner (the ``imp:`` group sources the complex
build)::

    ./run_examples.sh -e imp:3 -n 4 -t 480

Output lands in ``examples/implants/paraview_output/``:
``implants_03_rod_in_birdcage_tip_sar_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem, la

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE, ROD_RADIUS  # noqa: E402
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND  # noqa: E402
from tests.mesh.test_birdcage_port_tags import PHANTOM_HEIGHT, PHANTOM_RADIUS  # noqa: E402
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_birdcage_sar_1g_rung import ONE_GRAM_RUNG_CELL_RECORD  # noqa: E402
from tests.validation.test_birdcage_sar_mass_averaged import PHANTOM_RHO_KG_PER_M3  # noqa: E402
from tests.validation.test_imp3_rod_birdcage import (  # noqa: E402
    DRIVES,
    PEC_WALL_MAX_E,
    PORT_IDS,
    PRED_Q,
    PRED_Q_OVER_ESTIMATE,
    ROD_CELL_GUARD_W2,
    TIPS,
    _run_one_gram_rung,
)
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
)
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    ADJACENT_SPREAD_BAND,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "implants_03_rod_in_birdcage_tip_sar"

N_PROFILE = 12
PROFILE_DRIVES = ("P1", "ccw")
OFF_WALL_M = 0.5e-3
PRED_AXIAL_FIRST = (2.0, 20.0)
PRED_MIDPLANE_FIRST = (0.05, 0.9)
MAP_DRIVES = ("P1", "ccw")

# The hook's yield, keyed by arm ("off" / "rod").
STORE: dict = {}


def _profile_points():
    tip_z = TIPS["top"][2]
    z = np.linspace(tip_z + OFF_WALL_M, 0.5 * PHANTOM_HEIGHT - 2.0e-3, N_PROFILE)
    axial = np.column_stack([np.full(N_PROFILE, ROD_CENTRE[0]), np.zeros(N_PROFILE), z])
    x = np.linspace(ROD_CENTRE[0] + ROD_RADIUS + OFF_WALL_M, PHANTOM_RADIUS - 2.0e-3, N_PROFILE)
    mid = np.column_stack([x, np.zeros(N_PROFILE), np.zeros(N_PROFILE)])
    return {"axial": axial, "midplane": mid}


def _cellwise(msh, form_expr):
    """Per-cell integral of ``form_expr`` against a DG0 test function (owned+ghost)."""
    dg0 = fem.functionspace(msh, ("DG", 0))
    v = ufl.TestFunction(dg0)
    b = fem.assemble_vector(fem.form(form_expr(v)))
    b.scatter_reverse(la.InsertMode.add)
    b.scatter_forward()
    return dg0, np.real(np.array(b.array))


def _field_hook(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
    """Keep profiles and DG0 maps from one arm; never keeps ``e`` (it is freed next)."""
    comm = msh.comm
    arm = "rod" if rod else "off"
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags)
    sig = ufl.real(sigma_field)
    n_loc = msh.topology.index_map(msh.topology.dim).size_local

    prof = {}
    for name, pts in _profile_points().items():
        for d in PROFILE_DRIVES:
            vals, ok = evaluate_vector_field_parallel(e_of[d], pts)
            prof[(name, d)] = (np.sqrt(np.sum(np.abs(vals) ** 2, axis=1)), np.asarray(ok))

    entry = {"label": label, "profiles": prof, "loss_sum": {}, "fields": {}}
    dxa = ufl.Measure("dx", domain=msh)
    dg0, vol = _cellwise(msh, lambda v: ufl.conj(v) * dxa)
    for d in MAP_DRIVES:
        e = e_of[d]
        _, loss = _cellwise(
            msh, lambda v, e=e: 0.5 * sig * ufl.inner(e, e) * ufl.conj(v) * dx(PHANTOM_CELL_TAG)
        )
        entry["loss_sum"][d] = float(comm.allreduce(float(np.sum(loss[:n_loc])), op=MPI.SUM))
        if rod:
            _, e2 = _cellwise(msh, lambda v, e=e: ufl.inner(e, e) * ufl.conj(v) * dxa)
            emag = fem.Function(dg0, name=f"E_magnitude_rod_{d}")
            emag.x.array[:] = np.sqrt(np.maximum(e2 / vol, 0.0))
            emag.x.scatter_forward()
            lsar = fem.Function(dg0, name=f"local_sar_rod_{d}")
            lsar.x.array[:] = loss / (PHANTOM_RHO_KG_PER_M3 * vol)
            lsar.x.scatter_forward()
            entry["fields"][d] = (emag, lsar)
    STORE[arm] = entry


def _ball_mask(msh, dg0, centre, radius):
    """DG0 0/1 field: cell midpoint inside the ball (numpy, no UFL comparison)."""
    gdm = msh.geometry.dofmap
    xg = msh.geometry.x
    mid = xg[gdm].mean(axis=1)
    inside = np.sum((mid - np.asarray(centre)) ** 2, axis=1) < radius**2
    f = fem.Function(dg0)
    f.x.array[:] = 0.0
    f.x.array[dg0.dofmap.list[:, 0]] = inside.astype(float)
    f.x.scatter_forward()
    return f


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
    say("EX-83 -- a PEC rod in the coil-driven phantom: 1 g tip SAR with / without")
    say("=" * 72)
    say(
        "\n[scope] 1 g rung, degree 1, 64 MHz, untuned solid coil, 50 Ohm sheets, one "
        "axial rod on leg 1's azimuth. No accuracy of Q (ANS-10), no absolute SAR, no "
        "heating, no compliance."
    )

    # ---- the gate's own body, verbatim, with the hook --------------------------
    t0 = time.perf_counter()
    case = _run_one_gram_rung(field_hook=_field_hook)
    say(f"\n[gate reproduced] _run_one_gram_rung in {time.perf_counter() - t0:.1f} s")
    off, rod = case["off"], case["rod"]
    rod_case = rod["case"]
    msh = rod_case["mesh"]
    cell_tags = rod_case["cell_tags"]

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) --------------------
    port_names = {100 + i: f"port P{i} box (lower)" for i in range(1, 5)}
    port_names.update({200 + i: f"port P{i} box (upper)" for i in range(1, 5)})
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "conductor", 2: "air", 3: "phantom (saline)", **port_names},
        hide_tags=(2,),
        translucent_tags=(3,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:3 -- PEC rod in the driven birdcage phantom (IMP-3, 1 g rung)",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- anchors ---------------------------------------------------------------
    say("\n[anchors] asserted at the imported bands (compare 20261002T170843Z_IMP-3-step1-W2.log):")
    a0 = abs(off["cells"] / ONE_GRAM_RUNG_CELL_RECORD - 1.0)
    say(
        f"  (A0) without {off['cells']} vs {ONE_GRAM_RUNG_CELL_RECORD} ({a0:.3e}, band "
        f"{CELL_COUNT_BAND:.0%}); with {rod['cells']} (guard {ROD_CELL_GUARD_W2})"
    )
    assert a0 <= CELL_COUNT_BAND, a0
    assert rod["cells"] < ROD_CELL_GUARD_W2, rod["cells"]
    for arm in (off, rod):
        mq = arm["mq"]
        say(
            f"  (i)  {arm['label']}: ||S-S^T||/||S|| = {mq['recip']:.3e} (band {RECIPROCITY_BAND:g}); "
            f"(ii) sigma_max = {float(mq['sigma'][0]):.9f} (<= 1 + {PASSIVITY_SIGMA_TOLERANCE:g}), "
            f"max column power {float(np.max(mq['colpow'])):.9f} (<= 1)"
        )
        assert mq["recip"] <= RECIPROCITY_BAND, mq["recip"]
        assert float(mq["sigma"][0]) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, mq["sigma"]
        assert np.all(mq["colpow"] <= 1.0), mq["colpow"]
        for pid in PORT_IDS:
            ex = arm["shares"][pid]
            resid = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
            say(f"  (P)  {arm['label']} {pid}: rel {resid:.3e} (band {DISCRETE_IDENTITY_RTOL:g})")
            assert resid <= DISCRETE_IDENTITY_RTOL, (arm["label"], pid, resid)
    say(
        f"  (M)  rod-arm sigma_y pairs max {max(rod['mq']['mirror']):.4%} (band "
        f"{ADJACENT_SPREAD_BAND:.2%})"
    )
    for v in rod["mq"]["mirror"]:
        assert v <= ADJACENT_SPREAD_BAND, rod["mq"]["mirror"]
    say(
        f"  (E)  tag-402 dofs {rod['wall_dofs']}; max |E| per drive "
        + ", ".join(f"{p} {v:.3e}" for p, v in rod["wall_max"].items())
        + f" (bound {PEC_WALL_MAX_E:g})"
    )
    assert rod["wall_dofs"] > 0
    for pid, v in rod["wall_max"].items():
        assert v <= PEC_WALL_MAX_E, (pid, v)
    say(
        f"  (S)  without-rod separation {off['mq']['separation']:.1f}x (floor "
        f"{POOLED_SEPARATION_FLOOR:g})"
    )
    assert off["mq"]["separation"] >= POOLED_SEPARATION_FLOOR, off["mq"]["separation"]

    # ---- example-side checks ---------------------------------------------------
    say("\n[example-side, asserted] DG0 phantom-loss map sum vs the gate's P1 phantom power:")
    for key, arm in (("off", off), ("rod", rod)):
        ref = arm["shares"]["P1"]["phantom"]
        got = STORE[key]["loss_sum"]["P1"]
        rel = abs(got - ref) / abs(ref)
        say(f"  {arm['label']}: map {got:.9e} W vs shares {ref:.9e} W, rel {rel:.3e} (band {DISCRETE_IDENTITY_RTOL:g})")
        assert rel <= DISCRETE_IDENTITY_RTOL, (key, rel)
        for (name, d), (_, ok) in STORE[key]["profiles"].items():
            assert bool(np.all(ok)), (key, name, d, ok)
    say("  every profile point's ok mask is True on both arms")

    # ---- negative control: the rod-blind arm under-reads every tip --------------
    say("\n[control] Q = SAR_1g(rod)/SAR_1g(without): the rod-blind arm under-reads, Q > 1 ASSERTED:")
    qs = {(tip, d): rod["sar"][(tip, d)] / off["sar"][(tip, d)] for tip in TIPS for d in DRIVES}
    kmin = min(qs, key=qs.get)
    qmin = qs[kmin]
    for d in DRIVES:
        say(
            f"  {d:>3s}: top Q {qs[('top', d)]:.6f}  bottom Q {qs[('bottom', d)]:.6f}  "
            f"(gate PRED_Q [{PRED_Q[0]:g}, {PRED_Q[1]:g}]; Q/estimate band "
            f"[{PRED_Q_OVER_ESTIMATE[0]:g}, {PRED_Q_OVER_ESTIMATE[1]:g}] is the gate's)"
        )
    say(
        f"  min Q - 1 = {qmin - 1.0:.6f} ({kmin[1]} {kmin[0]}), 1 - 1/Q = {1.0 - 1.0 / qmin:.4%} "
        "(ceiling 100 %); W2 record minimum: P2 bottom Q = 1.576124 "
        "(20261002T170843Z_IMP-3-step1-W2.log:3867) -- printed comparator, not asserted"
    )
    for k, q in qs.items():
        assert q > 1.0, (k, q)
    say("  all 12 Q > 1 -- ASSERTED")

    # ---- printed: |E| ratios on two lines ---------------------------------------
    say("\n[printed, no band] |E|_rod / |E|_without at the same points:")
    firsts = {}
    for name in ("axial", "midplane"):
        pts = _profile_points()[name]
        for d in PROFILE_DRIVES:
            r = STORE["rod"]["profiles"][(name, d)][0] / STORE["off"]["profiles"][(name, d)][0]
            firsts[(name, d)] = float(r[0])
            say(f"  {name:>8s} {d:>3s}: " + " ".join(f"{x:.3f}" for x in r))
        if rank0:
            if name == "axial":
                say(f"           points z (mm): " + " ".join(f"{1e3 * z:.1f}" for z in pts[:, 2]))
            else:
                say(f"           points x (mm): " + " ".join(f"{1e3 * x:.1f}" for x in pts[:, 0]))
    ra = firsts[("axial", "P1")]
    rm = firsts[("midplane", "P1")]
    say(
        f"  first P1 point, axial: {ra:.3f} (predicted [{PRED_AXIAL_FIRST[0]:g}, "
        f"{PRED_AXIAL_FIRST[1]:g}], {'IN' if PRED_AXIAL_FIRST[0] <= ra <= PRED_AXIAL_FIRST[1] else 'OUT'}); "
        f"midplane: {rm:.3f} (predicted [{PRED_MIDPLANE_FIRST[0]:g}, {PRED_MIDPLANE_FIRST[1]:g}], "
        f"{'IN' if PRED_MIDPLANE_FIRST[0] <= rm <= PRED_MIDPLANE_FIRST[1] else 'OUT'})"
    )

    # ---- ParaView (rod mesh only) -----------------------------------------------
    dg0 = fem.functionspace(msh, ("DG", 0))
    ball_top = _ball_mask(msh, dg0, TIPS["top"], off["a_1g"])
    ball_top.name = "one_gram_ball_top"
    ball_bot = _ball_mask(msh, dg0, TIPS["bottom"], off["a_1g"])
    ball_bot.name = "one_gram_ball_bottom"
    fields = {}
    for d, (emag, lsar) in STORE["rod"]["fields"].items():
        suffix = "P1" if d == "P1" else "ccw"
        fields[f"E_magnitude_rod_{suffix}"] = emag
        fields[f"local_sar_rod_{suffix}"] = lsar
    fields["one_gram_ball_top"] = ball_top
    fields["one_gram_ball_bottom"] = ball_bot
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cell_tags,
        fields,
        comm=comm,
        facet_tags=rod["built"][2],
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    say(
        f"\n[paraview] wrote {xdmf_path} (rod mesh only; the without-rod arm is on a "
        "different mesh and is shown as profiles)"
        f"\n\nAll anchors and controls hold. Total elapsed {time.perf_counter() - started:.1f} s."
    )


if __name__ == "__main__":
    main()
