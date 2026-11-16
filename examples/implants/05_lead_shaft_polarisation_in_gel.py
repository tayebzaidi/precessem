"""Example (`EX-87`): a lead shaft in the gel -- the wire shorts the field that runs
along it and doubles the field that crosses it.

`IMP-4` step 1 gates an infinite, 2D-invariant PEC circular cylinder (a = 1 mm) in
`IMP-1`'s lossy host (sigma 0.57 S/m, eps_r 78) at 64 MHz and *finite k*, cut out of a
coaxial c = 10 mm, L = 4 mm slab as a hole (degree 2, 25 568 cells), against the
closed-form field for both drives: TM_z n = 0 (E along the shaft) and TE_z n = 1
(E across it).  This example runs those two gate solves through the gate's own
``_solve_cylinder``, asserts the gate's three conditions at the imported bands, and
shows what they mean: SAR in the annulus a < rho < 1.5a is ~0.4 % of open gel for the
axial drive and ~1.44x for the transverse one, a ratio of ~385.

**Import, never restate.** The fixture, solve, helpers and bands come from
``tests/validation/test_imp4_pec_cylinder_fem.py``.  Rule (a), disclosed: ``_solve_cylinder``
gained keyword-only ``return_fields`` / ``trace`` (defaults change nothing; the module
was re-gated green at -n 2, digit-identical to its W1 / W2 records).

**Anchors (asserted, imported bands).**  Each drive: |beta - 1| <= ``TH15_BETA_BAND``;
|Phi_h / Phi_cf,h - 1| <= ``MAT4_BAND``; cavity dofs > 0 with max |E| < 1e-12.  Each
value is printed beside the gate's W1 record.

**Example-side checks (asserted, true by construction).** (X1) the DG0 FEM SAR-ratio
map sums over the annulus to Phi_h; (X2) the DG0 closed-form map to Phi_cf,h (both within
1e-12); (X3) all solves report one cell count.

**Negative controls (asserted, wired).** The swapped end condition: C1 TM with natural
ends misses the Phi band; C2 TE with zeroed (PEC) ends misses the beta and Phi bands.
Their factors are printed beside predictions, never asserted.

**Printed, no band.** The polarisation contrast; pointwise |E|^2/E0^2 at z = 0; finite-k
against quasi-static for TE; a closed-form radius sweep for TM (numpy only); z-invariance.

Scope: one mesh, degree 2, a = 1 mm, an exact imposed trace, single-mode drives, ``-n 2``.
Not claimed: a finite, insulated or coated lead, a tip, a coil, heating, convergence, or
FEM at any other radius.

Run it through the example runner (the ``imp:`` group sources the complex build)::

    ./scripts/run_examples.sh -e imp:5 -n 2 -t 420

Output lands in ``examples/implants/paraview_output/``:
``implants_05_lead_shaft_polarisation_in_gel_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import basix
import dolfinx
import numpy as np
import ufl
from mpi4py import MPI
from petsc4py import PETSc

from dolfinx import default_scalar_type, fem, la

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.analytical import (  # noqa: E402
    pec_cylinder_te1_field,
    pec_cylinder_tm0_current,
    pec_cylinder_tm0_field,
)

from tests.validation.test_imp4_pec_cylinder_fem import (  # noqa: E402
    A,
    C,
    CAVITY_TAG,
    E0,
    K,
    LENGTH,
    MASK_FACTOR,
    MAT4_BAND,
    OMEGA,
    OUTER_BOUNDARY_TAG,
    PEC_TAGS,
    QDEG,
    TH15_BETA_BAND,
    _annulus_mask,
    _cf,
    _continuous_annulus_mean,
    _dirichlet,
    _eval_complex,
    _mesh,
    _owned_masked_cells,
    _rss_gib,
    _solve_cylinder,
    _static_te_numpy,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "implants_05_lead_shaft_polarisation_in_gel"
W1 = "20261003T200454Z_IMP-4-step1-W1.log"
END_TOL = 1.0e-6 * LENGTH
MAP_BAND = 1.0e-12
SWEEP_A_MM = (0.1, 0.25, 0.5, 1.0, 2.0)
SWEEP_A_NM = 1.0e-9
PRED_SWEEP_EZ = (0.0622, 0.0722, 0.0821, 0.0952, 0.1129, 0.0226)
PRED_SWEEP_I = (1.907e-3, 2.214e-3, 2.519e-3, 2.918e-3, 3.462e-3, 6.92e-4)

# (W1 record) beta, Phi_h, Phi_cf,h, n_cavity, cavity_max, z_inv
RECORD = {
    "tm": dict(beta=0.999655 - 0.000093j, phi_h=3.742390e-03, phi_cf_h=3.704934e-03,
               n_cavity=10704, cavity_max=3.726e-20, z_inv=7.335e-03, lines="133-137"),
    "te": dict(beta=0.992258 - 0.000001j, phi_h=1.440373, phi_cf_h=1.443759,
               n_cavity=10704, cavity_max=2.095e-34, z_inv=1.333e-02, lines="139-142"),
}


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
    say("EX-87 -- a lead shaft in the gel: shorted along it, doubled across it")
    say("=" * 72)
    say(
        f"\n[scope] PEC cylinder a = {A * 1e3:.1f} mm in c = {C * 1e3:.1f} mm, L = "
        f"{LENGTH * 1e3:.1f} mm slab as a hole, IMP-1 host at 64 MHz, k = {K:.6g} 1/m "
        f"(|k|a = {abs(K) * A:.4f}), degree 2, imposed exact trace, -n {comm.size}. "
        "No finite/insulated lead, no tip, no coil, no heating, no convergence."
    )

    msh, cell_tags, facet_tags = _mesh()
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_owned = imap.size_local
    n_local = n_owned + imap.num_ghosts
    all_cells = np.arange(n_local, dtype=np.int32)

    # ---- figure-only tags (never passed to a solve) ---------------------------
    mid = dolfinx.mesh.compute_midpoints(msh, tdim, all_cells)
    in_annulus = (mid[:, 0] ** 2 + mid[:, 1] ** 2) < (MASK_FACTOR * A) ** 2
    fig_vals = np.where(in_annulus, 3, 2).astype(np.int32)
    fig_tags = dolfinx.mesh.meshtags(msh, tdim, all_cells, fig_vals)
    write_setup_figure(
        msh,
        fig_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            2: "gel host sigma 0.57 S/m (1 mm PEC lead shaft cut out as a hole)",
            3: "annulus a < rho < 1.5a (Phi readout)",
        },
        translucent_tags=(2,),
        clip_normal=(0.0, -1.0, 0.0),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:5 -- 1 mm PEC lead shaft cut out of a gel host: z = 0 slice, graded ring",
        window_size=(1400, 630),
        comm=comm,
    )

    # ---- DG0 helpers -----------------------------------------------------------
    dg0 = fem.functionspace(msh, ("DG", 0))
    dof_of_cell = dg0.dofmap.list[:, 0]
    dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": QDEG})

    def cellwise(form_expr):
        v = ufl.TestFunction(dg0)
        b = fem.assemble_vector(fem.form(form_expr(v)))
        b.scatter_reverse(la.InsertMode.add)
        b.scatter_forward()
        return np.array(b.array)[dof_of_cell[all_cells]]

    def dg0_fn(per_cell, name):
        f = fem.Function(dg0, name=name)
        f.x.array[:] = 0.0
        f.x.array[dof_of_cell[all_cells]] = np.asarray(per_cell, dtype=float)
        f.x.scatter_forward()
        return f

    vol = np.real(cellwise(lambda v: ufl.conj(v) * dx))
    safe_vol = np.where(vol > 0.0, vol, 1.0)

    chi = _annulus_mask(msh)
    ann_cells = _owned_masked_cells(msh, chi)

    pts_ref, wts = basix.make_quadrature(basix.CellType.tetrahedron, QDEG)
    xg = msh.geometry.x
    gdm = msh.geometry.dofmap

    def cf_cell_map(field_fn):
        """Per-cell mean |E|^2/E0^2 of a closed form: the gate's basix rule, per cell."""
        v = xg[gdm[all_cells]]
        jac = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
        det = np.abs(np.linalg.det(jac))
        phys = v[:, 0][:, None, :] + np.einsum("cij,qj->cqi", jac, pts_ref)
        e = field_fn(phys.reshape(-1, 3)).reshape(len(all_cells), len(wts), 3)
        e2 = np.sum(np.abs(e) ** 2, axis=2)
        return np.sum(det[:, None] * wts[None, :] * e2, axis=1) / (
            det * np.sum(wts)
        ) / E0**2

    def masked_mean(per_cell):
        num = float(np.sum((vol * per_cell)[ann_cells]))
        den = float(np.sum(vol[ann_cells]))
        return comm.allreduce(num, op=MPI.SUM) / comm.allreduce(den, op=MPI.SUM)

    ann_flag = np.zeros(n_local)
    ann_flag[ann_cells] = 1.0
    annulus_mask_out = dg0_fn(ann_flag, "annulus_mask")

    probe_specs = [(1.2, 0.0), (1.2, np.pi / 2), (2.0, 0.0), (2.0, np.pi / 2)]
    PROBES = np.array([[s * A * np.cos(p), s * A * np.sin(p), 0.0] for s, p in probe_specs])

    outputs = {"annulus_mask": annulus_mask_out}
    results = {}
    cf_maps = {}

    def run(label, drive, tags, out_name, *, trace=None, with_e_real=False):
        t0 = time.perf_counter()
        d = _solve_cylinder(drive, tags, return_fields=True, trace=trace)
        fields = d.pop("fields")
        e = fields.e_complex
        num = cellwise(lambda v: ufl.inner(e, e) * ufl.conj(v) * dx)
        ratio = np.real(num) / (E0**2 * safe_vol)
        d["ratio_cells"] = ratio
        d["x1"] = abs(masked_mean(ratio) / d["phi_h"] - 1.0)
        outputs[f"sar_ratio_{out_name}"] = dg0_fn(ratio, f"sar_ratio_{out_name}")
        if with_e_real:
            v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
            e_cg = fem.Function(v_cg)
            e_cg.interpolate(e)
            e_re = fem.Function(v_cg, name=f"E_real_{out_name}")
            e_re.x.array[:] = np.real(e_cg.x.array)
            e_re.x.scatter_forward()
            outputs[f"E_real_{out_name}"] = e_re
        d["pt_e"] = _eval_complex(fields, PROBES, comm)
        d["wall"] = time.perf_counter() - t0
        del fields, e
        return d

    # ---- the two gate solves --------------------------------------------------
    for drive in ("tm", "te"):
        d = run("gate", drive, PEC_TAGS[drive], drive, with_e_real=True)
        cf = cf_cell_map(_cf(drive, "continue"))
        cf_maps[drive] = cf
        outputs[f"sar_ratio_{drive}_closed_form"] = dg0_fn(cf, f"sar_ratio_{drive}_closed_form")
        d["x2"] = abs(masked_mean(cf) / d["phi_cf_h"] - 1.0)
        results[drive] = d
    # the continuous Phi comes with the gate dict
    for drive, d in results.items():
        rec = RECORD[drive]
        say(
            f"\n[gate {drive}] pec tags {PEC_TAGS[drive]}: cells {d['ncells']} (W1 25568), "
            f"dofs {d['ndofs']} (174112), annulus cells {d['n_mask']} (9735), "
            f"h_wall/a {d['h_wall']:.3f} (0.396), {d['seconds']:.1f} s, "
            f"RSS(sum) {d['rss']:.2f} GiB"
            f"\n  (i) beta = {d['beta'].real:.6f}{d['beta'].imag:+.6f}j, |beta - 1| = "
            f"{abs(d['beta'] - 1):.4%} vs band {TH15_BETA_BAND:.4%} "
            f"(W1 record {rec['beta'].real:.6f}{rec['beta'].imag:+.6f}j, {W1}:{rec['lines']})"
            f"\n  (ii) Phi_h = {d['phi_h']:.6e}, Phi_cf,h = {d['phi_cf_h']:.6e}, "
            f"|Phi_h/Phi_cf,h - 1| = {abs(d['phi_h'] / d['phi_cf_h'] - 1):.4%} vs band "
            f"{MAT4_BAND:.0%} (record {rec['phi_h']:.6e} / {rec['phi_cf_h']:.6e})"
            f"\n  (iii) cavity dofs {d['n_cavity']} (record {rec['n_cavity']}), max |E| "
            f"{d['cavity_max']:.3e} (record {rec['cavity_max']:.3e})"
            f"\n  z-invariance {d['z_inv']:.3e} (record {rec['z_inv']:.3e}; cause untested)"
        )
        assert d["n_cavity"] > 0, drive
        assert d["cavity_max"] < 1e-12, (drive, d["cavity_max"])
        assert abs(d["beta"] - 1.0) <= TH15_BETA_BAND, (
            drive, d["beta"], "an example/test divergence, not a tolerance to move")
        assert abs(d["phi_h"] / d["phi_cf_h"] - 1.0) <= MAT4_BAND, (
            drive, d["phi_h"], d["phi_cf_h"], "an example/test divergence")
    i_cf = pec_cylinder_tm0_current(k=K, radius=A, omega=OMEGA, e0=E0)
    i_fit = results["tm"]["beta"] * i_cf
    say(
        f"  TM current beta*I = {abs(i_fit):.6g} A per V/m, |I| (record 2.91733e-03, "
        f"{W1}:137), closed form |I_cf| {abs(i_cf):.6g}"
    )

    say("\n[example-side, asserted]")
    ncells = {d["ncells"] for d in results.values()}
    for drive, d in results.items():
        say(
            f"  (X1) {drive}: DG0 FEM map annulus mean vs Phi_h: {d['x1']:.3e} (band "
            f"{MAP_BAND:g}); (X2) closed-form map vs Phi_cf,h: {d['x2']:.3e} (band {MAP_BAND:g})"
        )
        assert d["x1"] <= MAP_BAND, (drive, "X1", d["x1"])
        assert d["x2"] <= MAP_BAND, (drive, "X2", d["x2"])

    # ---- negative controls ----------------------------------------------------
    c1 = run("C1", "tm", (OUTER_BOUNDARY_TAG, CAVITY_TAG), "tm_natural_ends")
    ncells.add(c1["ncells"])

    def t0_trace(x):
        return _dirichlet("te")(x) * (np.abs(x[2]) < LENGTH / 2 - END_TOL)

    def trace_zeroed(x):
        out = np.asarray(t0_trace(x))
        return np.ascontiguousarray(out).astype(PETSc.ScalarType)

    c2 = run("C2", "te", PEC_TAGS["tm"], "te_pec_ends", trace=trace_zeroed)
    ncells.add(c2["ncells"])
    say(f"  (X3) all four solves report cells {sorted(ncells)} (one cached mesh)")
    assert len(ncells) == 1, ncells

    miss1 = c1["phi_h"] / c1["phi_cf_h"] - 1.0
    say(
        "\n[control C1] TM, natural ends (pec tags (1, 2)): the wrong end condition for E_z"
        f"\n  Phi_h = {c1['phi_h']:.6e} vs Phi_cf,h {c1['phi_cf_h']:.6e}: miss {miss1:+.4%} "
        f"(band {MAT4_BAND:.0%}; predicted [-100 %, -90 %], "
        f"{'IN' if -1.0 <= miss1 <= -0.9 else 'OUT'}; asserted only as > band)"
        f"\n  beta = {c1['beta'].real:.6f}{c1['beta'].imag:+.6f}j, |beta - 1| = "
        f"{abs(c1['beta'] - 1):.4f} (printed, predicted ~0.11 [0.05, 0.2], "
        f"{'IN' if 0.05 <= abs(c1['beta'] - 1) <= 0.2 else 'OUT'}); z-invariance "
        f"{c1['z_inv']:.3e} (predicted [0.1, 1.0], {'IN' if 0.1 <= c1['z_inv'] <= 1.0 else 'OUT'}); cavity max |E| {c1['cavity_max']:.3e}"
    )
    assert abs(miss1) > MAT4_BAND, "C1 is inside the Phi band: the swapped end condition is not reaching the probes"

    miss2 = c2["phi_h"] / c2["phi_cf_h"] - 1.0
    say(
        "\n[control C2] TE, PEC ends with the trace zeroed on the end planes (pec tags (1, 2, 3))"
        f"\n  beta = {c2['beta'].real:.6f}{c2['beta'].imag:+.6f}j, |beta - 1| = "
        f"{abs(c2['beta'] - 1):.4f} vs band {TH15_BETA_BAND:.4%} (predicted [0.9, 1.1], "
        f"{'IN' if 0.9 <= abs(c2['beta'] - 1) <= 1.1 else 'OUT'}; asserted > band)"
        f"\n  Phi_h = {c2['phi_h']:.6e} vs Phi_cf,h {c2['phi_cf_h']:.6e}: miss {miss2:+.4%} "
        f"(band {MAT4_BAND:.0%}; predicted [-100 %, -90 %], "
        f"{'IN' if -1.0 <= miss2 <= -0.9 else 'OUT'}; asserted > band); z-invariance "
        f"{c2['z_inv']:.3e} (predicted [0.1, 1.0], {'IN' if 0.1 <= c2['z_inv'] <= 1.0 else 'OUT'}); cavity max |E| {c2['cavity_max']:.3e}"
        "\n  (the gate's own TE natural-cavity control can miss only on beta; zeroing the "
        "trace removes the field outright, so Phi discriminates here)"
    )
    assert abs(c2["beta"] - 1.0) > TH15_BETA_BAND, "C2 beta inside band: zeroing not applied"
    assert abs(miss2) > MAT4_BAND, "C2 Phi inside band: zeroing not applied"

    # ---- printed: polarisation contrast ----------------------------------------
    tm, te = results["tm"], results["te"]
    say(
        "\n[printed] polarisation contrast Phi_TE / Phi_TM: h "
        f"{te['phi_h'] / tm['phi_h']:.1f} (predicted 384.9), closed form on the cells "
        f"{te['phi_cf_h'] / tm['phi_cf_h']:.1f} (389.7), continuous "
        f"{te['phi_cont'] / tm['phi_cont']:.1f} (396.9)"
    )

    # ---- printed: pointwise |E|^2/E0^2 -----------------------------------------
    cf_pts = {d_: np.sum(np.abs(_cf(d_)(PROBES)) ** 2, axis=1) / E0**2 for d_ in ("tm", "te")}
    fem_pts = {d_: np.sum(np.abs(results[d_]["pt_e"]) ** 2, axis=1) / E0**2 for d_ in ("tm", "te")}
    pred_pts = {
        "te": (0.0934, 2.871, 0.5625, 1.5625),
        "tm": (1.83e-3, 1.83e-3, 2.65e-2, 2.65e-2),
    }
    say("\n[printed] pointwise |E|^2/E0^2 at z = 0 (phi = pi/2 faces the field, phi = 0 is the flank):")
    for drive in ("te", "tm"):
        for i, (s, p) in enumerate(probe_specs):
            rel = fem_pts[drive][i] / cf_pts[drive][i] - 1.0
            tol = 0.10 if (drive == "te" and p == 0.0) else 0.05
            say(
                f"  {drive} rho {s}a phi {np.degrees(p):4.0f} deg: FEM {fem_pts[drive][i]:.5g}, "
                f"closed form {cf_pts[drive][i]:.5g} (predicted {pred_pts[drive][i]:.4g}), "
                f"FEM/cf - 1 = {rel:+.3%} (predicted within +-{tol:.0%}, "
                f"{'IN' if abs(rel) <= tol else 'OUT'})"
            )
    surf = np.array([[0.0, A, 0.0]])
    te_surf = float(np.sum(np.abs(_cf("te", "continue")(surf)) ** 2) / E0**2)
    tm_surf = float(np.sum(np.abs(_cf("tm", "continue")(surf)) ** 2) / E0**2)
    say(f"  surface rho = a, phi = pi/2: TE {te_surf:.4f} (predicted 4.0009), TM {tm_surf:.3e} (predicted 0)")

    # ---- printed: finite k vs quasi-static, TE ---------------------------------
    te_fk = _continuous_annulus_mean(_cf("te", "continue"))
    te_st = _continuous_annulus_mean(_static_te_numpy)
    rel = te_fk / te_st - 1.0
    say(
        f"\n[printed] TE annulus mean, continuous: finite k {te_fk:.6f}, quasi-static "
        f"{te_st:.6f} (predicted 1.444444); finite-k/static - 1 = {rel:+.3e} (predicted "
        f"+1.9e-4 [+1.5e-4, +2.3e-4], {'IN' if 1.5e-4 <= rel <= 2.3e-4 else 'OUT'})"
    )

    # ---- printed: closed-form radius sweep (numpy only) ------------------------
    say(
        "\n[printed] closed-form radius sweep at fixed k (numpy, no FEM); a = 1 nm is "
        "k -> 1e-6 k at 1 mm:"
    )
    radii = [a_mm * 1e-3 for a_mm in SWEEP_A_MM] + [SWEEP_A_NM]
    labels = [f"{a_mm:g} mm" for a_mm in SWEEP_A_MM] + ["1 nm"]
    for r, lab, p_ez, p_i in zip(radii, labels, PRED_SWEEP_EZ, PRED_SWEEP_I):
        p = np.array([[1.5 * r, 0.0, 0.0]])
        ez = float(abs(pec_cylinder_tm0_field(p, k=K, radius=r, e0=E0, interior="continue")[0, 2]) / E0)
        cur = float(abs(pec_cylinder_tm0_current(k=K, radius=r, omega=OMEGA, e0=E0)))
        pe = np.array([[0.0, 1.5 * r, 0.0]])
        te_v = float(np.linalg.norm(
            pec_cylinder_te1_field(pe, k=K, radius=r, omega=OMEGA, e0=E0, interior="continue")[0]
        ) / E0)
        tol = 0.05 if lab == "1 nm" else 0.03
        say(
            f"  a = {lab:>7}: TM |E_z(1.5a)|/E0 = {ez:.4f} (predicted {p_ez}, "
            f"{'IN' if abs(ez / p_ez - 1) <= tol else 'OUT'}); |I| = {cur:.4g} A per V/m "
            f"(predicted {p_i:.4g}, {'IN' if abs(cur / p_i - 1) <= tol else 'OUT'}); "
            f"TE |E(1.5a, pi/2)|/E0 = {te_v:.5f} (predicted 1.4444, "
            f"{'IN' if abs(te_v - 1.4444) <= 1e-3 * 1.4444 else 'OUT'})"
        )

    rss_all = _rss_gib(comm)  # collective

    # ---- ParaView --------------------------------------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cell_tags,
        outputs,
        comm=comm,
        facet_tags=facet_tags,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    say(
        f"\n[paraview] wrote {xdmf_path}"
        "\n[paraview] Slice z = 0; colour `sar_ratio_te` on a log scale, then `sar_ratio_tm`."
        f"\n\nAll assertions hold. RSS(sum) {rss_all:.2f} GiB. Total elapsed "
        f"{time.perf_counter() - started:.1f} s."
    )


if __name__ == "__main__":
    main()
