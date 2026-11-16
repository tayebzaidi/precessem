"""Example (`EX-84`): one adjoint solve gives an implant lead's transfer function.

`imp:1` / `imp:2` / `imp:3` read the field near a PEC body. This example does the
other thing an implant study needs: it computes the *transfer function* TF(z) of a
sigma = 800 S/m lead (2 x 2 x 24 mm) in `IMP-1`'s lossy host at 64 MHz (degree 1,
complex, one structured box) with **one adjoint solve** read along the lead, and
then folds that TF against the lead-free field of three drives -- the gate's D1
(uniform Dirichlet E0) and D2 (an impressed z-block) plus a new drive D3 (D2's block
mirrored through x = 0 and z = 0) that the adjoint never saw -- and checks each fold
against the direct with-lead tip response q_with - q_without.

**Import, never restate.** The fixture, the adjoint solve, every band and every
record come from ``tests/validation/test_wf9_transfer_function_identity.py``
(`WF-9` step 0, §9 item 183).  ``_run_identities`` is the body of the gate's own
``test_wf9_identities`` lifted to module level (rule (a), additive: a trailing
``return`` and nothing else; the gate was re-run from this tree at ``-n 2`` and
``-n 4`` in the same slot, digit-identical).  ``DISCRETE_IDENTITY_RTOL`` comes from
the module the gate imports it from.

**Anchors (asserted on the returned case, bands imported).**
(A1) ||A - A^T||_F/||A||_F <= ``SYMMETRY_RTOL``; (T2) TF_step vs TF_adj and (H)
held vs fresh <= ``DISCRETE_IDENTITY_RTOL``; (T1) superposition and (T3) the row's
identity -j w mu0 dsigma int chi_lead w . E_inc = q_with - q_without on D1 and D2 at
``DISCRETE_IDENTITY_RTOL`` x |dq|.  Round-off readings vary run to run and are
asserted at their bands only.

**Example-side checks (asserted, true by construction).**
(X0) the D3 block census (owned count = ``BLOCK_CELL_COUNT``, disjoint from the lead,
the tip ball and D2's block); (X1) the per-cell TF-density vector sums per layer to
``tf_adj``; (X2) the per-cell fold-integrand vector sums to the (T3) fold on D1 and
D2; (X3) the new drive: a lead-out held back-substitution gives E_inc,D3, a fresh
lead-in solve gives E_with,D3, and the fold of the gate's adjoint field against
E_inc,D3 reproduces q_with - q_without on D3.

**Negative controls (asserted, each wired).** C1: the Hermitian read-back differs
from the adjoint TF (A = A^T, not A = A^H; ceiling 200 %, no factor claimed).  C2: the
total-field fold misses on D1, D2 and D3.  C3 (the Born control): the lead-free
adjoint w0 folded against the lead-free field misses dq on D1, D2 and D3 -- the TF
must come from an adjoint solve with the lead present.  Factors are printed against
predictions, never asserted.

**Printed, no band.** The TF(z) table; D3's |q_with/q_without| and its ISO path fold
at N = 24 with the transverse term beside it (the path fold's band is the weekly's);
the gate's own W1/W2 TF agreement.

**Scope.** Degree 1, one structured box, one square sigma = 800 S/m lead in `IMP-1`'s
host (sigma 0.57 S/m, eps_r 78), 64 MHz, three drives, a linear tip functional.  Not
claimed: the path fold's band, ISO 10974's measured TF or its calibration, a PEC or
insulated lead, a coil or F-small drive, a gel-phantom fixture, heating, SAR or
compliance, h-convergence, `WF-9`'s check mark.

Run it through the example runner (the ``imp:`` group sources the complex build)::

    ./run_examples.sh -e imp:4 -n 2 -t 420

Output lands in ``examples/implants/paraview_output/``:
``implants_04_lead_transfer_function_adjoint_combined.xdmf``.
"""

from __future__ import annotations

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
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_wf9_transfer_function_identity import (  # noqa: E402
    BLOCK_CELL_COUNT,
    DELTA_SIGMA,
    HOST_TAG,
    J_BLOCK,
    LEAD_HALF_Z,
    LEAD_TAG,
    LOAD,
    N_LAYERS,
    SYMMETRY_RTOL,
    _free,
    _problem,
    _run_identities,
    _scalar,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "implants_04_lead_transfer_function_adjoint"

# D3: D2's block (4 < x < 6, |y| < 1, 2 < z < 10 mm) mirrored through x = 0 and z = 0.
D3_X = (-0.006, -0.004)
D3_Y = 0.001
D3_Z = (-0.010, -0.002)
PRED_PEAK_LAYER = 19
PRED_TAIL_PHASE_WINDOW_DEG = 20.0
PRED_D3_RATIO = (0.05, 50.0)
PRED_PATH_FOLD = (1.0e-4, 1.0e-1)
PRED_C1 = (1.90, 2.00)
PRED_C2 = (0.80, 1.05)
PRED_C3_D1 = (3.0, 300.0)
PRED_C3_D23 = (1.0, 1000.0)


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

    def band(x, lo, hi):
        return "IN" if lo <= x <= hi else "OUT"

    say("=" * 72)
    say("EX-84 -- one adjoint solve gives an implant lead's transfer function")
    say("=" * 72)
    say(
        "\n[scope] degree 1, one box, one sigma = 800 S/m lead in IMP-1's host, 64 MHz, "
        "three drives, a linear tip functional. No path-fold band, no ISO 10974 TF, no "
        "coil, no heating, no compliance."
    )

    # ---- the gate's own body, verbatim ----------------------------------------
    t0 = time.perf_counter()
    case = _run_identities()
    say(f"\n[gate reproduced] _run_identities in {time.perf_counter() - t0:.1f} s")
    fx = case["fx"]
    msh = fx["mesh"]
    flags = fx["flags"]
    n_owned = fx["n_owned"]
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_local = n_owned + imap.num_ghosts
    cells = np.arange(n_local, dtype=np.int32)
    q = case["q"]
    w = case["w"]
    results = case["results"]
    tf_adj = case["tf_adj"]
    tf_step = case["tf_step"]
    dx = ufl.dx(domain=msh)
    zhat = ufl.as_vector((0.0, 0.0, 1.0))
    chi_lead = fx["chi_lead"]
    tol = DISCRETE_IDENTITY_RTOL

    dg0 = fem.functionspace(msh, ("DG", 0))
    dof_of_cell = dg0.dofmap.list[:, 0]

    def mask_fn(flag_arr, name):
        f = fem.Function(dg0, name=name)
        f.x.array[:] = 0.0
        f.x.array[dof_of_cell[cells]] = flag_arr.astype(float)
        f.x.scatter_forward()
        return f

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

    # ---- anchors ---------------------------------------------------------------
    rec = "20261003T003327Z_WF-9-step0-W1.log:96-124"
    say(f"\n[anchors] asserted at the imported bands, on the returned case (compare {rec}):")
    say(f"  (A1) ||A-A^T||_F/||A||_F = {case['sym_ratio']:.3e} (band {SYMMETRY_RTOL:g}; record 1.102e-16)")
    say(f"  (T2) max|TF_step-TF_adj|/max|TF_adj| = {case['t2_dev']:.3e} (band {tol:g}; record 1.478e-12)")
    say(f"  (H)  |TF23_fresh-TF23_held|/|TF23| = {case['h_dev']:.3e} (band {tol:g}; record 1.234e-14)")
    assert case["sym_ratio"] <= SYMMETRY_RTOL, case["sym_ratio"]
    assert case["t2_dev"] <= tol, case["t2_dev"]
    assert case["h_dev"] <= tol, case["h_dev"]
    for name in ("D1", "D2"):
        r = results[name]
        dq = r["dq"]
        t1 = abs(r["t1_sum"] - dq) / abs(dq)
        t3 = abs(r["t3_fold"] - dq) / abs(dq)
        say(
            f"  {name}: q_without={r['q_without']:.10e} q_with={r['q_with']:.10e}\n"
            f"      (T1) {t1:.3e}  (T3) {t3:.3e}  (band {tol:g} x |dq|)"
        )
        assert abs(r["t1_sum"] - dq) <= tol * abs(dq), (name, "T1", t1)
        assert abs(r["t3_fold"] - dq) <= tol * abs(dq), (name, "T3", t3)
    say(
        "  records: D1 (T1) 1.734e-13 (T3) 1.447e-12; D2 (T1) 4.770e-14 (T3) 1.096e-12; "
        "D1 q_without 1.0024793765e+00-5.1409041670e-03j, D2 q_without 4.7106960161e-03-2.3488734907e-03j"
    )

    # ---- (X0) the D3 census ----------------------------------------------------
    say("\n[example-side, asserted]")
    mid = dolfinx.mesh.compute_midpoints(msh, tdim, cells)
    x, y, z = mid[:, 0], mid[:, 1], mid[:, 2]
    d3 = (x > D3_X[0]) & (x < D3_X[1]) & (np.abs(y) < D3_Y) & (z > D3_Z[0]) & (z < D3_Z[1])
    n_d3 = comm.allreduce(int(np.count_nonzero(d3[:n_owned])), op=MPI.SUM)
    disjoint = comm.allreduce(
        bool(
            not (
                np.any(d3 & flags["lead"])
                or np.any(d3 & flags["ball"])
                or np.any(d3 & flags["block"])
            )
        ),
        op=MPI.LAND,
    )
    say(
        f"  (X0) D3 block: {n_d3} owned cells (= BLOCK_CELL_COUNT {BLOCK_CELL_COUNT}); "
        f"disjoint from lead, tip ball and D2's block: {disjoint}"
    )
    assert n_d3 == BLOCK_CELL_COUNT, n_d3
    assert disjoint
    chi_d3 = mask_fn(d3, "chi_d3_mask")

    # ---- (X1) the TF-density map -----------------------------------------------
    tfd = cellwise(lambda v: LOAD * DELTA_SIGMA * chi_lead * w[2] * ufl.conj(v) * dx)
    x1 = 0.0
    for j in range(N_LAYERS):
        s = comm.allreduce(complex(np.sum(tfd[:n_owned][flags["layers"][j][:n_owned]])), op=MPI.SUM)
        x1 = max(x1, abs(s - tf_adj[j]) / abs(tf_adj[j]))
    say(f"  (X1) TF-density map: max_j |sum_cells - tf_adj[j]|/|tf_adj[j]| = {x1:.3e} over {N_LAYERS} layers (band {tol:g})")
    assert x1 <= tol, x1

    # ---- (X2) the fold-integrand maps ------------------------------------------
    fold_maps = {}
    for name in ("D1", "D2"):
        e_inc = case["drives"][name]["inc"]
        vec = cellwise(lambda v, e=e_inc: LOAD * DELTA_SIGMA * chi_lead * ufl.dot(w, e) * ufl.conj(v) * dx)
        fold_maps[name] = vec
        s = comm.allreduce(complex(np.sum(vec[:n_owned])), op=MPI.SUM)
        x2 = abs(s - results[name]["t3_fold"]) / abs(results[name]["t3_fold"])
        say(f"  (X2) {name} fold-integrand map: |sum - t3_fold|/|t3_fold| = {x2:.3e} (band {tol:g})")
        assert x2 <= tol, (name, x2)

    # ---- (X3) the new drive D3 -------------------------------------------------
    t0 = time.perf_counter()
    s_out = TimeHarmonicSolver(_problem(fx, with_lead=False, dirichlet=False), degree=1)
    w0 = s_out.solve(extra_linear_terms=[case["ell"]]).e_complex
    e_inc_d3 = s_out.solve_with_held_factorization(
        extra_linear_terms=[lambda v: LOAD * ufl.inner(J_BLOCK * chi_d3 * zhat, v) * ufl.dx]
    ).e_complex
    _free(s_out)
    s_in = TimeHarmonicSolver(_problem(fx, with_lead=True, dirichlet=False), degree=1)
    e_with_d3 = s_in.solve(
        current_density=lambda xx: J_BLOCK * chi_d3 * zhat, project_source=False
    ).e_complex
    _free(s_in)
    t_x3 = time.perf_counter() - t0

    q_without3 = q(e_inc_d3)
    q_with3 = q(e_with_d3)
    dq3 = q_with3 - q_without3
    fold3 = LOAD * DELTA_SIGMA * _scalar(comm, chi_lead * ufl.dot(w, e_inc_d3) * dx)
    x3 = abs(fold3 - dq3) / abs(dq3)
    say(
        f"  (X3) D3 (two solves, {t_x3:.1f} s): q_without={q_without3:.10e} q_with={q_with3:.10e}\n"
        f"       dq = {dq3:.10e}; fold(w, E_inc,D3) = {fold3:.10e}; |fold - dq|/|dq| = {x3:.3e} (band {tol:g})"
    )
    assert abs(fold3 - dq3) <= tol * abs(dq3), (fold3, dq3, x3)
    say("  all of (X0)-(X3) hold")

    # ---- negative controls -----------------------------------------------------
    say("\n[controls] each ASSERTED wired:")
    c1 = case["c1_dev"]
    say(
        f"  C1 Hermitian read-back deviation = {c1:.6%} (predicted [{PRED_C1[0]:.0%}, {PRED_C1[1]:.0%}], "
        f"{band(c1, *PRED_C1)}; record 199.871741%; ceiling 200%, no factor claimed)"
    )
    assert c1 > tol
    dq_all = {"D1": results["D1"]["dq"], "D2": results["D2"]["dq"], "D3": dq3}
    e_inc_all = {"D1": case["drives"]["D1"]["inc"], "D2": case["drives"]["D2"]["inc"], "D3": e_inc_d3}
    e_with_all = {"D1": case["drives"]["D1"]["with"], "D2": case["drives"]["D2"]["with"], "D3": e_with_d3}
    c3_ratio = {}
    for name in ("D1", "D2", "D3"):
        dq = dq_all[name]
        c2_fold = LOAD * DELTA_SIGMA * _scalar(comm, chi_lead * ufl.dot(w, e_with_all[name]) * dx)
        c2 = abs(c2_fold - dq) / abs(dq)
        qb = LOAD * DELTA_SIGMA * _scalar(comm, chi_lead * ufl.dot(w0, e_inc_all[name]) * dx)
        c3 = abs(qb - dq) / abs(dq)
        c3_ratio[name] = abs(qb) / abs(dq)
        pr = PRED_C3_D1 if name == "D1" else PRED_C3_D23
        say(
            f"  {name}: C2 total-field fold miss = {c2:.6%} (predicted [{PRED_C2[0]:.0%}, {PRED_C2[1]:.0%}], "
            f"{band(c2, *PRED_C2)}); C3 Born miss |Q_B-dq|/|dq| = {c3:.4g}, |Q_B|/|dq| = {c3_ratio[name]:.4g} "
            f"(predicted [{pr[0]:g}, {pr[1]:g}], {band(c3_ratio[name], *pr)})"
        )
        assert abs(c2_fold - dq) > tol * abs(dq), (name, "C2 not wired")
        assert abs(qb - dq) > tol * abs(dq), (name, "C3 not wired")

    # ---- printed: the TF(z) table ----------------------------------------------
    say("\n[printed, no band] TF(z) from one adjoint solve (j, z_mid mm, |TF_adj|, |TF|/max, arg deg):")
    tmax = float(np.max(np.abs(tf_adj)))
    for j in range(N_LAYERS):
        say(
            f"  {j:2d} {-1e3 * LEAD_HALF_Z + 0.5 + j:6.1f} {abs(tf_adj[j]):.6e} "
            f"{abs(tf_adj[j]) / tmax:.4f} {np.degrees(np.angle(tf_adj[j])):9.4f}"
        )
    jpk = int(np.argmax(np.abs(tf_adj)))
    ph = {j: float(np.degrees(np.angle(tf_adj[j]))) for j in (22, 23)}
    near = {j: min(abs(abs(p) - 180.0), 360.0) for j, p in ph.items()}
    say(
        f"  peak at j = {jpk} (predicted {PRED_PEAK_LAYER}, {'IN' if jpk == PRED_PEAK_LAYER else 'OUT'}; W1:145); "
        f"j=22 / 23 phase {ph[22]:.2f} / {ph[23]:.2f} deg, distance from +-180: {near[22]:.2f} / {near[23]:.2f} "
        f"(predicted <= {PRED_TAIL_PHASE_WINDOW_DEG:g}, "
        f"{'IN' if max(near.values()) <= PRED_TAIL_PHASE_WINDOW_DEG else 'OUT'}; W1:148-149)"
    )

    r3 = abs(q_with3) / abs(q_without3)
    ez = np.array([_scalar(comm, fx["chi_layers"][j] * e_inc_d3[2] * dx) for j in range(N_LAYERS)])
    lv = np.array([_scalar(comm, fx["chi_layers"][j] * dx).real for j in range(N_LAYERS)])
    q24 = 0.0 + 0.0j
    for k in range(N_LAYERS):  # N = 24: one segment per layer
        q24 += tf_step[k] * (ez[k] / lv[k])
    transverse3 = LOAD * DELTA_SIGMA * _scalar(
        comm, chi_lead * (w[0] * e_inc_d3[0] + w[1] * e_inc_d3[1]) * dx
    )
    pf = abs(q24 - dq3) / abs(dq3)
    say(
        f"\n  D3: |q_with/q_without| = {r3:.6f} (predicted [{PRED_D3_RATIO[0]:g}, {PRED_D3_RATIO[1]:g}], "
        f"{band(r3, *PRED_D3_RATIO)})"
    )
    say(
        f"  D3 ISO path fold N=24: Q_N = {q24:.10e}, |Q_N-dq|/|dq| = {pf:.6e} "
        f"(predicted [{PRED_PATH_FOLD[0]:g}, {PRED_PATH_FOLD[1]:g}], {band(pf, *PRED_PATH_FOLD)}); "
        f"transverse term / |dq| = {abs(transverse3) / abs(dq3):.6e}  (the band is the weekly's)"
    )
    say("  W1/W2 TF agreement: the gate body prints it above, against the other width's exchange file.")
    say(f"  gate body: dofs={case['ndofs']}, elapsed {case['elapsed']:.1f} s, summed max RSS {case['rss']:.3f} GiB")

    # ---- setup figure (no-op unless FEM_EM_SETUP_FIGURES=1) --------------------
    # Figure-only tags, never passed to a solve: 1 host, 2 lead, 3 tip ball, 4 D2, 5 D3.
    tag_vals = np.full(n_local, HOST_TAG, dtype=np.int32)
    tag_vals[flags["block"]] = 4
    tag_vals[d3] = 5
    tag_vals[flags["ball"]] = 3
    tag_vals[flags["lead"]] = LEAD_TAG
    fig_tags = dolfinx.mesh.meshtags(msh, tdim, cells, tag_vals)
    write_setup_figure(
        msh,
        fig_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: "host phantom (sigma 0.57 S/m)",
            2: "lead wire (sigma 800 S/m)",
            3: "tip ball (q readout)",
            4: "drive block D2",
            5: "drive block D3",
        },
        translucent_tags=(1,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:4 -- lead transfer function from one adjoint solve (WF-9 fixture)",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- ParaView ----------------------------------------------------------------
    mag = np.abs(tfd)
    phase = np.where(mag > 0.0, np.degrees(np.angle(tfd)), 0.0)
    safe_vol = np.where(vol > 0.0, vol, 1.0)
    w2 = np.real(cellwise(lambda v: ufl.inner(w, w) * ufl.conj(v) * dx))
    ei3 = np.real(cellwise(lambda v: ufl.inner(e_inc_d3, e_inc_d3) * ufl.conj(v) * dx))
    fold_d3 = cellwise(lambda v: LOAD * DELTA_SIGMA * chi_lead * ufl.dot(w, e_inc_d3) * ufl.conj(v) * dx)
    layer_idx = np.zeros(n_local)
    for j in range(N_LAYERS):
        layer_idx[flags["layers"][j]] = j + 1
    fields = {
        "tf_density_magnitude": field_fn(mag, "tf_density_magnitude"),
        "tf_density_phase_deg": field_fn(phase, "tf_density_phase_deg"),
        "adjoint_w_magnitude": field_fn(np.sqrt(np.maximum(w2 / safe_vol, 0.0)), "adjoint_w_magnitude"),
        "fold_integrand_magnitude_D2": field_fn(np.abs(fold_maps["D2"]), "fold_integrand_magnitude_D2"),
        "fold_integrand_magnitude_D3": field_fn(np.abs(fold_d3), "fold_integrand_magnitude_D3"),
        "E_inc_magnitude_D3": field_fn(np.sqrt(np.maximum(ei3 / safe_vol, 0.0)), "E_inc_magnitude_D3"),
        "lead_layer_index": field_fn(layer_idx, "lead_layer_index"),
        "tip_ball_mask": mask_fn(flags["ball"], "tip_ball_mask"),
        "drive_block_D2": mask_fn(flags["block"], "drive_block_D2"),
        "drive_block_D3": mask_fn(d3, "drive_block_D3"),
    }
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        fx["cell_tags"],
        fields,
        comm=comm,
        facet_tags=None,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    say(
        f"\n[paraview] wrote {xdmf_path}"
        f"\n\nAll anchors, checks and controls hold. Total elapsed {time.perf_counter() - started:.1f} s."
    )


if __name__ == "__main__":
    main()
