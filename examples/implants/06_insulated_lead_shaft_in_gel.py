"""Example (`EX-88`): an insulated lead shaft in the gel -- the sleeve takes the field that
crosses the wire and moves its hot spot 90 degrees; the field that runs along it does not
notice the sleeve.

`IMP-5` step 1 gates an infinite PEC circular cylinder (a = 1 mm) wearing a lossless
eps_r 3 sleeve (a < rho < b = 1.5 mm, its own conforming cell tag 4) in `IMP-1`'s lossy
host (sigma 0.57 S/m, eps_r 78) at 64 MHz and finite k, cut out of a coaxial c = 10 mm,
L = 4 mm slab as a hole (degree 2, 25 878 cells), against step 0's exact coated
Bessel/Hankel series for both drives: TM_z n = 0 (E along the shaft) and TE_z n = 1
(E across it).  This example runs those gate solves through the gate's own
``_solve_coated``, asserts the gate's conditions at the imported bands, and shows what the
sleeve does: normal D is continuous at rho = b, so the normal E jumps by the ratio of the
complex permittivities (about 59) and the field crowds into the insulation.

**Import, never restate.** Fixture, solve, helpers and bands come from
``tests/validation/test_imp5_coated_cylinder_fem.py`` (and the bare-wire closed form and
trace from ``tests/validation/test_imp4_pec_cylinder_fem.py``).  Rule (a), disclosed:
``_solve_coated`` gained keyword-only ``return_fields`` / ``trace`` (defaults change
nothing; the module was re-gated green at -n 2, digit-identical to its records).

**Anchors (asserted, imported bands).** Each drive: |beta - 1| <= ``TH15_BETA_BAND``;
|Phi_h / Phi_cf,h - 1| <= ``MAT4_BAND`` on the sleeve (tag 4) and the host annulus
(tag 2, rho < 1.5b); cavity dofs > 0 with max |E| < 1e-12.

**Example-side checks (asserted, true by construction).** (X1) DG0 FEM |E|^2 map means
equal the returned Phi_h; (X2) the DG0 closed-form map (region forced by cell tag) means
equal Phi_cf,h (both within 1e-12); (X3) all solves report one cell count.

**Negative controls (asserted).** C1 TE with the coating omitted (the gate's own control);
B0 TE bare wire on the same mesh (bare trace) -- each must miss the coated gate.

Scope: one mesh, degree 2, a = 1 mm, b = 1.5 mm, lossless eps_r 3, an exact imposed trace,
single-mode k_z = 0 drives, ``-n 2``.  Not claimed: k_z != 0 guided waves along a real
insulated lead, lossy / thin insulation, a finite lead, a tip, a coil, heating, ISO 10974,
convergence.

Run it through the example runner (the ``imp:`` group sources the complex build)::

    ./scripts/run_examples.sh -e imp:6 -n 2 -t 600

Output lands in ``examples/implants/paraview_output/``:
``implants_06_insulated_lead_shaft_in_gel_combined.xdmf``.
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
    coated_pec_cylinder_tm0_current,
    pec_cylinder_tm0_current,
)

from tests.validation.test_imp4_pec_cylinder_fem import (  # noqa: E402
    _cf as _cf_bare,
    _dirichlet as _dirichlet_bare,
)
from tests.validation.test_imp5_coated_cylinder_fem import (  # noqa: E402
    A,
    B,
    C,
    COAT_TAG,
    E0,
    EPS_R_COAT,
    HOST_MASK_FACTOR,
    HOST_TAG,
    K,
    K_COAT,
    LENGTH,
    MASK_RANGE,
    MAT4_BAND,
    OMEGA,
    PEC_TAGS,
    QDEG,
    TH15_BETA_BAND,
    _cf,
    _continuous_mean,
    _eval_complex,
    _incident,
    _masks,
    _mesh,
    _owned_masked_cells,
    _probe_points,
    _rss_gib,
    _solve_coated,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "implants_06_insulated_lead_shaft_in_gel"
W1A = "20261004T035027Z_IMP-5-step1-W1a.log"
W1B = "20261004T035247Z_IMP-5-step1-W1b.log"
W2A = "20261004T035543Z_IMP-5-step1-W2a.log"
STEP0 = "20261004T003556Z_IMP-5-step0-W1b.log"
MAP_BAND = 1.0e-12

# (records) beta, Phi_h, Phi_cf,h per mask, cavity, z_inv
RECORD = {
    "tm": dict(beta=0.999703 - 0.000089j, coat=(3.634296e-03, 3.596754e-03),
               host=(2.315906e-02, 2.307528e-02), n_cavity=10704, cavity_max=0.0,
               z_inv=4.868e-03, src=f"{W1A}:145-150"),
    "te": dict(beta=0.994548 - 0.000062j, coat=(1.799486e01, 1.802471e01),
               host=(1.410011, 1.412624), n_cavity=10704, cavity_max=7.272e-34,
               z_inv=2.364e-02, src=f"{W1B}:145-149"),
}
RECORD_C1 = dict(beta_miss=1.4735, coat=-0.914590, host=-0.180148, src=f"{W2A}:146-148")

# pointwise predictions: |E|^2/E0^2 at z = 0, (rho label, phi) -> (coated, bare)
PT_SPECS = [("1.1a", 1.1 * A), ("0.97b", 0.97 * B), ("1.03b", 1.03 * B), ("2b", 2.0 * B)]
PT_PHIS = [("facing", np.pi / 2), ("flank", 0.0)]
PT_PRED = {
    ("1.1a", "facing"): (41.56, 3.336), ("1.1a", "flank"): (0.375, 0.0301),
    ("0.97b", "facing"): (27.01, 2.168), ("0.97b", "flank"): (3.468, 0.278),
    ("1.03b", "facing"): (0.0143, 2.013), ("1.03b", "flank"): (3.632, 0.338),
    ("2b", "facing"): (0.578, 1.235), ("2b", "flank"): (1.538, 0.790),
}
PT_SMALL = {("1.03b", "facing", "coated"), ("1.1a", "flank", "coated"),
            ("1.1a", "flank", "bare")}


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
    say("EX-88 -- an insulated lead shaft in the gel: the sleeve rotates the hot spot")
    say("=" * 72)
    say(
        f"\n[scope] PEC cylinder a = {A * 1e3:.1f} mm, eps_r {EPS_R_COAT:g} lossless sleeve to "
        f"b = {B * 1e3:.2f} mm, in c = {C * 1e3:.1f} mm, L = {LENGTH * 1e3:.1f} mm slab as a "
        f"hole, IMP-1 host at 64 MHz, k = {K:.6g} 1/m, degree 2, imposed exact trace, "
        f"-n {comm.size}. k_z = 0 only: no k_z != 0 guided wave, no tip, no coil, no heating."
    )

    msh, cell_tags, facet_tags = _mesh()
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_local = imap.size_local + imap.num_ghosts
    all_cells = np.arange(n_local, dtype=np.int32)
    tag_arr = np.zeros(n_local, dtype=np.int32)
    tag_arr[cell_tags.indices] = cell_tags.values
    coat_cells = np.flatnonzero(tag_arr == COAT_TAG)
    host_cells = np.flatnonzero(tag_arr == HOST_TAG)

    # ---- figure-only tags (never passed to a solve) ---------------------------
    mid = dolfinx.mesh.compute_midpoints(msh, tdim, all_cells)
    rho2 = mid[:, 0] ** 2 + mid[:, 1] ** 2
    fig_vals = np.full(n_local, 2, dtype=np.int32)
    fig_vals[host_cells[rho2[host_cells] < (HOST_MASK_FACTOR * B) ** 2]] = 3
    fig_vals[coat_cells] = 4
    fig_tags = dolfinx.mesh.meshtags(msh, tdim, all_cells, fig_vals)
    write_setup_figure(
        msh,
        fig_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            2: "gel host sigma 0.57 S/m (1 mm PEC lead cut out as a hole)",
            3: "host annulus b < rho < 1.5b (Phi readout)",
            4: "insulation sleeve eps_r 3, sigma 0, a < rho < b = 1.5 mm",
        },
        translucent_tags=(2,),
        clip_normal=(0.0, -1.0, 0.0),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:6 -- 1 mm PEC lead with an eps_r 3 sleeve in gel: z = 0 slice",
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

    masks = _masks(msh, cell_tags)
    mask_cells = {name: _owned_masked_cells(msh, chi) for name, chi in masks.items()}

    pts_ref, wts = basix.make_quadrature(basix.CellType.tetrahedron, QDEG)
    xg = msh.geometry.x
    gdm = msh.geometry.dofmap

    def cf_cell_map(field_fn, cells):
        """Per-cell mean |E|^2/E0^2 of a closed form on ``cells`` (the gate's basix rule)."""
        out = np.zeros(n_local)
        if cells.size == 0:
            return out
        v = xg[gdm[cells]]
        jac = np.stack([v[:, 1] - v[:, 0], v[:, 2] - v[:, 0], v[:, 3] - v[:, 0]], axis=2)
        det = np.abs(np.linalg.det(jac))
        phys = v[:, 0][:, None, :] + np.einsum("cij,qj->cqi", jac, pts_ref)
        e = field_fn(phys.reshape(-1, 3)).reshape(len(cells), len(wts), 3)
        e2 = np.sum(np.abs(e) ** 2, axis=2)
        out[cells] = np.sum(det[:, None] * wts[None, :] * e2, axis=1) / (
            det * np.sum(wts)) / E0**2
        return out

    def masked_mean(per_cell, name):
        cells = mask_cells[name]
        num = float(np.sum((vol * per_cell)[cells]))
        den = float(np.sum(vol[cells]))
        return comm.allreduce(num, op=MPI.SUM) / comm.allreduce(den, op=MPI.SUM)

    def cf_by_tag(drive):
        """Closed-form map, the region forced by cell tag (never by radius)."""
        return (cf_cell_map(_cf(drive, region="coat"), coat_cells)
                + cf_cell_map(_cf(drive, region="host"), host_cells))

    ann_flag = np.zeros(n_local)
    ann_flag[host_cells[rho2[host_cells] < (HOST_MASK_FACTOR * B) ** 2]] = 1.0
    outputs = {"host_annulus_mask": dg0_fn(ann_flag, "host_annulus_mask")}

    # probe points (z = 0)
    probe_list = [(lab, phi_lab, np.array([[r * np.cos(p), r * np.sin(p), 0.0]]))
                  for lab, r in PT_SPECS for phi_lab, p in PT_PHIS]
    PROBES = np.vstack([p for _, _, p in probe_list])
    JUMP_PTS = np.array([[0.0, 0.97 * B, 0.0], [0.0, 1.03 * B, 0.0]])
    gate_pts = _probe_points()

    def run(drive, tags, out_name, *, coated=True, trace=None, with_e_real=False,
            store_ratio=True):
        t0 = time.perf_counter()
        d = _solve_coated(drive, tags, coated=coated, return_fields=True, trace=trace)
        fields = d.pop("fields")
        e = fields.e_complex
        num = cellwise(lambda v: ufl.inner(e, e) * ufl.conj(v) * dx)
        ratio = np.real(num) / (E0**2 * safe_vol)
        d["ratio_cells"] = ratio
        d["x1"] = max(abs(masked_mean(ratio, m) / d["phi"][m]["h"] - 1.0) for m in masks)
        if store_ratio:
            outputs[f"e2_ratio_{out_name}"] = dg0_fn(ratio, f"e2_ratio_{out_name}")
        if with_e_real:
            v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
            e_cg = fem.Function(v_cg)
            e_cg.interpolate(e)
            e_re = fem.Function(v_cg, name=f"E_real_{out_name}")
            e_re.x.array[:] = np.real(e_cg.x.array)
            e_re.x.scatter_forward()
            outputs[f"E_real_{out_name}"] = e_re
        d["pt_e"] = _eval_complex(fields, PROBES, comm)
        d["jump_e"] = _eval_complex(fields, JUMP_PTS, comm)
        d["gate_e"] = _eval_complex(fields, gate_pts, comm)
        d["wall"] = time.perf_counter() - t0
        del fields, e
        return d

    # ---- the two gate solves (coated) ------------------------------------------
    results = {}
    for drive in ("tm", "te"):
        d = run(drive, PEC_TAGS[drive], drive, with_e_real=True)
        cf = cf_by_tag(drive)
        d["x2"] = max(abs(masked_mean(cf, m) / d["phi"][m]["cf"] - 1.0) for m in masks)
        if drive == "te":
            outputs["e2_ratio_te_closed_form"] = dg0_fn(cf, "e2_ratio_te_closed_form")
            sar = np.where(tag_arr == HOST_TAG, d["ratio_cells"], 0.0)
            outputs["sar_ratio_te"] = dg0_fn(sar, "sar_ratio_te")
        results[drive] = d

    for drive, d in results.items():
        rec = RECORD[drive]
        say(
            f"\n[gate {drive}] pec tags {PEC_TAGS[drive]}: cells {d['ncells']} (W 25878), per tag "
            f"{d['ntag']}, dofs {d['ndofs']} (176216), h_wall/a {d['h_wall']:.3f} (0.383), "
            f"sleeve max edge/(b-a) {d['sleeve_edge']:.3f} (0.858), {d['seconds']:.1f} s, "
            f"RSS(sum) {d['rss']:.2f} GiB"
            f"\n  (i) beta = {d['beta'].real:.6f}{d['beta'].imag:+.6f}j, |beta - 1| = "
            f"{abs(d['beta'] - 1):.4%} vs band {TH15_BETA_BAND:.4%} (record "
            f"{rec['beta'].real:.6f}{rec['beta'].imag:+.6f}j, {rec['src']})"
        )
        for name, label in (("coat", "sleeve"), ("host", "host annulus")):
            p = d["phi"][name]
            say(
                f"  (ii/iii) {label}: Phi_h = {p['h']:.6e}, Phi_cf,h = {p['cf']:.6e}, "
                f"Phi_h/Phi_cf,h - 1 = {p['h'] / p['cf'] - 1:+.4%} vs band {MAT4_BAND:.0%} "
                f"(record {rec[name][0]:.6e} / {rec[name][1]:.6e}); continuous Phi "
                f"{p['cont']:.6e}"
            )
        say(
            f"  (iv) cavity dofs {d['n_cavity']} (record {rec['n_cavity']}), max |E| "
            f"{d['cavity_max']:.3e} (record {rec['cavity_max']:.3e}); z-invariance "
            f"{d['z_inv']:.3e} (record {rec['z_inv']:.3e}; cause untested)"
        )
        assert d["n_cavity"] > 0, drive
        assert d["cavity_max"] < 1e-12, (drive, d["cavity_max"])
        assert abs(d["beta"] - 1.0) <= TH15_BETA_BAND, (
            drive, d["beta"], "an example/test divergence, not a tolerance to move")
        for name in ("coat", "host"):
            p = d["phi"][name]
            assert abs(p["h"] / p["cf"] - 1.0) <= MAT4_BAND, (
                drive, name, p, "an example/test divergence")

    say("\n[example-side, asserted]")
    for drive, d in results.items():
        say(
            f"  (X1) {drive}: DG0 FEM |E|^2 map means vs Phi_h (both masks): {d['x1']:.3e}; "
            f"(X2) closed-form map (region by tag) vs Phi_cf,h: {d['x2']:.3e} "
            f"(band {MAP_BAND:g})"
        )
        assert d["x1"] <= MAP_BAND, (drive, "X1", d["x1"])
        assert d["x2"] <= MAP_BAND, (drive, "X2", d["x2"])

    te = results["te"]
    tm = results["tm"]

    # ---- C1: coating omitted ------------------------------------------------------
    c1 = run("te", PEC_TAGS["te"], "te_coating_omitted", coated=False)
    c1_miss = {m: c1["phi"][m]["h"] / c1["phi"][m]["cf"] - 1.0 for m in masks}
    say(
        "\n[control C1] TE, coating omitted (host material in tag 4, the coated trace on the wall)"
        f"\n  |beta - 1| = {abs(c1['beta'] - 1):.4f} (record {RECORD_C1['beta_miss']}, "
        f"{RECORD_C1['src']}); band {TH15_BETA_BAND:.4%}; sleeve miss {c1_miss['coat']:+.4%} "
        f"(record {RECORD_C1['coat']:+.4%}); host miss {c1_miss['host']:+.4%} (record "
        f"{RECORD_C1['host']:+.4%}); asserted only as > band"
    )
    assert abs(c1["beta"] - 1.0) > TH15_BETA_BAND, "C1 beta inside band: material_map not reaching tag 4"
    for m in masks:
        assert abs(c1_miss[m]) > MAT4_BAND, ("C1", m, "inside the band")
    c1_ncells = c1["ncells"]
    del c1

    # ---- B0: the bare wire on the same mesh --------------------------------------
    b0 = run("te", PEC_TAGS["te"], "te_bare", coated=False, trace=_dirichlet_bare("te"),
             with_e_real=True)
    b0_miss = {m: b0["phi"][m]["h"] / b0["phi"][m]["cf"] - 1.0 for m in masks}
    cf_b0_fn = _cf_bare("te", "continue")
    cf_b0 = cf_cell_map(cf_b0_fn, all_cells)
    outputs["e2_ratio_te_bare_closed_form"] = dg0_fn(cf_b0, "e2_ratio_te_bare_closed_form")
    outputs["sar_ratio_te_bare"] = dg0_fn(b0["ratio_cells"], "sar_ratio_te_bare")
    b_lo, b_hi = 1.3, 1.65
    say(
        "\n[control B0] TE, the bare wire on the same mesh (sleeve cells = gel, bare trace)"
        f"\n  against the COATED closed form: beta = {b0['beta'].real:.6f}{b0['beta'].imag:+.6f}j "
        f"(predicted -0.4606-0.0354j), |beta - 1| = {abs(b0['beta'] - 1):.4f} (predicted "
        f"[{b_lo}, {b_hi}], {'IN' if b_lo <= abs(b0['beta'] - 1) <= b_hi else 'OUT'}; "
        f"asserted > band {TH15_BETA_BAND:.4%})"
        f"\n  sleeve miss {b0_miss['coat']:+.4%} (predicted [-95, -88] %, "
        f"{'IN' if -0.95 <= b0_miss['coat'] <= -0.88 else 'OUT'}); host miss "
        f"{b0_miss['host']:+.4%} (predicted [-27, -19] %, "
        f"{'IN' if -0.27 <= b0_miss['host'] <= -0.19 else 'OUT'}); both asserted > band "
        f"{MAT4_BAND:.0%}; z-invariance {b0['z_inv']:.3e}; cavity max |E| {b0['cavity_max']:.3e}"
    )
    assert abs(b0["beta"] - 1.0) > TH15_BETA_BAND, (
        "B0 beta inside band: the bare trace is not applied", b0["beta"], b0["z_inv"])
    for m in masks:
        assert abs(b0_miss[m]) > MAT4_BAND, ("B0", m, b0["phi"][m], b0["z_inv"])

    # B0 against its own closed form (printed only)
    e_s_b = cf_b0_fn(gate_pts) - _incident("te")(gate_pts)
    beta_b = complex(np.sum(np.conj(e_s_b) * (b0["gate_e"] - _incident("te")(gate_pts)))
                     / np.sum(np.abs(e_s_b) ** 2))
    say(
        f"  against its OWN bare closed form (printed, no band): beta_bare = "
        f"{beta_b.real:.6f}{beta_b.imag:+.6f}j, |beta_bare - 1| = {abs(beta_b - 1):.4%} "
        f"(predicted <~ 2 %, {'IN' if abs(beta_b - 1) <= 0.02 else 'OUT'}; IMP-4 TE 0.774 %)"
    )
    for m, pred_cont in (("coat", 1.4447), ("host", 1.0878)):
        cf_mean = masked_mean(cf_b0, m)
        cont = _continuous_mean(cf_b0_fn, *MASK_RANGE[m])
        mm = b0["phi"][m]["h"] / cf_mean - 1.0
        say(
            f"  bare {m}: Phi_h {b0['phi'][m]['h']:.6f}, own-cf cell mean {cf_mean:.6f}, "
            f"Phi_h/Phi_cf - 1 = {mm:+.3%} (predicted within +-2 %, "
            f"{'IN' if abs(mm) <= 0.02 else 'OUT'}); continuous Phi {cont:.4f} "
            f"(predicted {pred_cont})"
        )

    ncells = {tm["ncells"], te["ncells"], c1_ncells, b0["ncells"]}
    say(f"\n  (X3) all four solves report cells {sorted(ncells)} (one cached mesh)")
    assert len(ncells) == 1, ncells

    # ---- printed: the interface jump ------------------------------------------------
    pb = np.array([[0.0, B, 0.0]])
    ey_in = abs(_cf("te", region="coat")(pb)[0, 1])
    ey_out = abs(_cf("te", region="host")(pb)[0, 1])
    jump_cf = ey_in / ey_out
    cf_in = abs(_cf("te", region="coat")(JUMP_PTS[:1])[0, 1])
    cf_out = abs(_cf("te", region="host")(JUMP_PTS[1:])[0, 1])
    jump_cf2 = cf_in / cf_out
    jump_fem = abs(te["jump_e"][0, 1]) / abs(te["jump_e"][1, 1])
    bare_fn = cf_b0_fn
    bare_ratio = abs(bare_fn(JUMP_PTS[:1])[0, 1]) / abs(bare_fn(JUMP_PTS[1:])[0, 1])
    say(
        "\n[printed] TE normal-E jump at phi = pi/2 (E_y), coated:"
        f"\n  closed form at exactly rho = b, coat/host: {jump_cf:.2f} (predicted 59.36, {STEP0}:115, "
        f"{'IN' if abs(jump_cf / 59.36 - 1) < 0.01 else 'OUT'})"
        f"\n  closed form 0.97b / 1.03b: {jump_cf2:.2f} (predicted 43.5 [42, 45], "
        f"{'IN' if 42 <= jump_cf2 <= 45 else 'OUT'}); FEM {jump_fem:.2f}, FEM/cf - 1 = "
        f"{jump_fem / jump_cf2 - 1:+.3%} (predicted within +-15 %, "
        f"{'IN' if abs(jump_fem / jump_cf2 - 1) <= 0.15 else 'OUT'})"
        f"\n  bare wire (B0 closed form) at the same two points: {bare_ratio:.3f} "
        f"(predicted 1.038, {'IN' if abs(bare_ratio - 1.038) < 0.01 else 'OUT'}; no jump)"
    )

    # ---- printed: pointwise |E|^2/E0^2 ----------------------------------------------
    say("\n[printed] pointwise |E|^2/E0^2 at z = 0, FEM beside closed form (coated TE | bare B0):")
    for i, (lab, phi_lab, pt) in enumerate(probe_list):
        region = "coat" if lab in ("1.1a", "0.97b") else "host"
        cf_c = float(np.sum(np.abs(_cf("te", region=region)(pt)) ** 2) / E0**2)
        cf_bb = float(np.sum(np.abs(bare_fn(pt)) ** 2) / E0**2)
        fem_c = float(np.sum(np.abs(te["pt_e"][i]) ** 2) / E0**2)
        fem_b = float(np.sum(np.abs(b0["pt_e"][i]) ** 2) / E0**2)
        pc, pbare = PT_PRED[(lab, phi_lab)]
        for kind, f_v, c_v, p_v in (("coated", fem_c, cf_c, pc), ("bare", fem_b, cf_bb, pbare)):
            tol = 0.15 if (lab, phi_lab, kind) in PT_SMALL else 0.05
            ctol = 0.05 if (lab, phi_lab, kind) == ("1.03b", "facing", "coated") else 0.01
            rel = f_v / c_v - 1.0
            say(
                f"  {kind:6s} rho {lab:>5s} {phi_lab:6s}: FEM {f_v:.5g}, closed form {c_v:.5g} "
                f"(predicted {p_v}, cf {'IN' if abs(c_v / p_v - 1) <= ctol else 'OUT'}), "
                f"FEM/cf - 1 = {rel:+.3%} (+-{tol:.0%}, {'IN' if abs(rel) <= tol else 'OUT'})"
            )

    # ---- printed: coated against bare on one mesh -----------------------------------
    r_host = te["phi"]["host"]["h"] / b0["phi"]["host"]["h"]
    r_coat = te["phi"]["coat"]["h"] / b0["phi"]["coat"]["h"]
    say(
        "\n[printed] coated against bare on one mesh (TE):"
        f"\n  host-annulus Phi_h coated/bare = {r_host:.3f} (predicted 1.297 [1.26, 1.34], "
        f"{'IN' if 1.26 <= r_host <= 1.34 else 'OUT'}; {STEP0}:116)"
        f"\n  sleeve-region ratio = {r_coat:.2f} (predicted 12.45 [12.0, 12.9], "
        f"{'IN' if 12.0 <= r_coat <= 12.9 else 'OUT'}); this is |E|^2 in a sigma = 0 "
        "dielectric (dielectric stress), not SAR"
    )

    # ---- printed: TM, said honestly ---------------------------------------------------
    i_c = coated_pec_cylinder_tm0_current(k_host=K, k_coat=K_COAT, radius=A, coat_radius=B,
                                          omega=OMEGA, e0=E0)
    i_b = pec_cylinder_tm0_current(k=K, radius=A, omega=OMEGA, e0=E0)
    di = abs(i_c / i_b - 1.0)
    fem_tm = abs(tm["beta"] - 1.0)
    say(
        "\n[printed] TM (E along the shaft): a non-magnetic sleeve is nearly invisible at k_z = 0"
        f"\n  closed form |I_coated/I_bare - 1| = {di:.3e} (predicted 4.27e-5, {STEP0}:113); "
        f"s shift 3.084e-6 (cited, same line)"
        f"\n  FEM gate |beta - 1| = {fem_tm:.3e} = {fem_tm / di:.1f}x the current shift and "
        f"{fem_tm / 3.084e-6:.0f}x the s shift (predicted ~7x, ~100x): the FEM cannot see the "
        "coating in TM; the TM gate is not a coating test"
        f"\n  TM sleeve continuous Phi = {tm['phi']['coat']['cont']:.4e} (predicted 3.640e-3)"
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
        "\n[paraview] Slice z = 0; colour `e2_ratio_te` on a log scale (the jump at rho = b); "
        "then `sar_ratio_te` and `sar_ratio_te_bare` on one linear scale. `E_real_*` are CG1 "
        "and smear the jump: read it from the DG0 maps."
        f"\n\nAll assertions hold. RSS(sum) {rss_all:.2f} GiB. Total elapsed "
        f"{time.perf_counter() - started:.1f} s."
    )


if __name__ == "__main__":
    main()
