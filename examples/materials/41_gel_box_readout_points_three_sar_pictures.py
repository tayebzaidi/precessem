"""Example (`EX-192`): entry 59's six readout points and C6's 18-point lattice under three
closed-form SAR pictures, per deposited watt (`MAT-21` steps 42, 44, 45 / P75, P80, P81).

`mat:23` shows the eddy source's FEM and `mat:25` its duty-cycled FEM at the six probes.  Neither
shows C6's lattice, the end-weighted source at the probes, or entry 52's SAR10g shape beside
the closed-form ball shapes.  This page does.  There is no solve: the box is a tetrahedral
`create_box` mesh (m = 2, (26, 16, 4) cells) that carries three P1 fields -- the unwalled
six-face box's theta(t*) per deposited watt at t* = 2 821.5 s for the uniform source (P75's U3),
P36's eddy source (P80) and step 6's end-weighted source (P81).  All three deposit the same
watt, so each column is theta per deposited watt (the gel watt; entry 59 is per accepted watt).

Imports, by module object: P75 (``test_mat21_wall_uniform_source_probe_points_closed_form``),
P80 (``test_mat21_eddy_source_probe_points_closed_form``), P81
(``test_mat21_end_source_probe_points_closed_form``).  Every constant, band and helper is read
from them; no ``test_*`` function and no fixture is called.  No gate edit, no ``src/`` edit.

Asserted (all collected first, then asserted): P80's (P-eq), (Q2-p) bitwise on P75's batch,
(Rep-s), (Lim), (Ball), (K-lat), (c-u), (c-ball); P81's (P-eq), (Q2-s) bitwise, (Rep-s) on step
6's own scalar calls, (Lim), (Lin), (Ball), (c-u), (c-e), (c-ball); (V-pt) at the 24 readout
vertices each P1 field read through ``evaluate_vector_field_parallel(owner=
"lowest_original_cell")`` equals its series within ``B1_BAND`` relative; (cV) the control that
at six z-edge midpoints the P1 field misses its series by more than ``B1_BAND`` for all three
fields; max|Im theta| = 0.  Printed only: the six-probe table, the lattice table with entry 52's
SAR10g shape, the x+0mm column profile, tie sets within 1e-4 relative.

Not claimed: which picture entry 59 follows, a fitted or mixed source, any entry-59 /
entry-60 figure (P76 reads them), a dose band, a compliance figure.

Run it::

    ./run_examples.sh -e mat:41 -n 2 -t 480

Output: ``examples/materials/paraview_output/
materials_41_gel_box_readout_points_three_sar_pictures_combined.xdmf`` (the three per-watt
theta(t*) fields on P1).
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.utils.analytical import (  # noqa: E402
    robin_box_quadratic_steady_temperature,
    robin_box_separable_steady_temperature,
    robin_box_transient_temperature,
)

from tests.validation import test_geo38_f_clinical_box_loaded as c6  # noqa: E402
from tests.validation import test_mat21_end_source_probe_points_closed_form as p81  # noqa: E402
from tests.validation import test_mat21_eddy_source_probe_points_closed_form as p80  # noqa: E402
from tests.validation import (  # noqa: E402
    test_mat21_wall_uniform_source_probe_points_closed_form as p75,
)
from tests.validation.test_imp4_pec_cylinder_closed_form import B1_BAND  # noqa: E402
from tests.validation.test_mat13_transient_perfusion_sphere import INTERIOR_BAND  # noqa: E402
from tests.validation.test_mat21_box_closed_form import (  # noqa: E402
    A_BOX,
    B_BOX,
    C_BOX,
    HALF,
    SERIES_TERMS,
    T_C,
    T_INF,
    H,
    K,
    Q,
)
from tests.validation.test_mat21_box_closed_form import RHO_C as RHO_C_0  # noqa: E402
from tests.validation.test_mat21_box_quadratic_closed_form import (  # noqa: E402
    QUAD_TERMS,
    UNIT,
    WEIGHTS,
)
from tests.validation.test_mat21_box_separable_closed_form import (  # noqa: E402
    FLIPPED,
    PROFILES,
    SEP_TERMS,
    UNIFORM,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_41_gel_box_readout_points_three_sar_pictures"
EPS = float(np.finfo(float).eps)
CELLS = (26, 16, 4)  # step 3's m = 2: BASE (13, 8, 2) x 2
CHUNK = 2000  # points per series call, the pattern of test_mat21_box_fem.py:94-108
PICS = ("U3", "eddy", "end")


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _u3_at(P, t):
    P = np.asarray(P, dtype=float)
    return np.asarray(robin_box_transient_temperature(
        P[:, 0], P[:, 1], P[:, 2], t, half_widths=HALF, k=K, q=Q, h=H, rho_c=RHO_C_0,
        t_inf=0.0, n_terms=SERIES_TERMS), dtype=float)


def _series(pic, P):
    """The series at points P (n, 3), per deposited watt, at t*."""
    t = p75.T_STAR
    if pic == "U3":
        v = _u3_at(P, t)
    elif pic == "eddy":
        v = p80._eddy_at(P, t)
    else:
        v = p81._end_at(P, t)
    return v / p75.P_DEP


def _interp(msh, pic):
    """Chunked P1 interpolation of the series (disclosed local copy of the _exact pattern)."""
    V = fem.functionspace(msh, ("Lagrange", 1))
    g = fem.Function(V, name=f"theta_{pic}_per_W_t_star")

    def ev(x):
        n = x.shape[1]
        out = np.empty(n)
        for s in range(0, n, CHUNK):
            out[s:s + CHUNK] = _series(pic, x[:, s:s + CHUNK].T)
        return out

    g.interpolate(ev)
    g.x.scatter_forward()
    return g


def _read(f, P):
    vals, ok = evaluate_vector_field_parallel(f, np.ascontiguousarray(P),
                                              owner="lowest_original_cell")
    return np.real(vals[:, 0]), bool(np.all(ok))


def closed_form(comm, gates):
    """P80's and P81's gates re-asserted on the page, through their imported helpers."""
    d80 = p80._data()
    d81 = p81._data()
    P = p80._probe_arr()
    ts = p75.T_STAR
    labs6 = p75.PROBE_LABELS
    r = float(c6._a10())
    lat_pts = d80["lat_pts"]
    lat_labs = [lab for lab, _p, _w in d80["lat"]]

    # (P-eq): three pictures, one watt
    p_eddy = Q * 8.0 * A_BOX * B_BOX * C_BOX * (WEIGHTS[0] + (WEIGHTS[1] + WEIGHTS[2]
                                                                + WEIGHTS[3]) / 3.0)
    m = 1.0
    for a, b, j in PROFILES:
        m *= a if b == 0.0 else a + b * np.sin(j * np.pi) / (j * np.pi)
    p_end = Q * 8.0 * A_BOX * B_BOX * C_BOX * m
    pe_e, pe_n = abs(p_eddy / p75.P_DEP - 1.0), abs(p_end / p75.P_DEP - 1.0)
    gates["P-eq eddy"] = pe_e <= 4 * EPS
    gates["P-eq end"] = pe_n <= 4 * EPS
    _say(comm, f"[P-eq] P_dep {p75.P_DEP:.6f} W (uniform); eddy |ratio-1| {pe_e:.3e}, end "
               f"|ratio-1| {pe_n:.3e} (bound 4 eps {4 * EPS:.3e}, asserted)")
    _say(comm, "[norm] every column is per deposited gel watt; entry 59 is per accepted watt "
               "(theta_acc = theta * P_gel / P_acc)")

    # (Q2-p) / (Q2-s): the unit-weight / uniform member is P75's U3, bitwise, on P75's batch
    u3 = p75._u3(ts)
    q2p = np.array_equal(p80._eddy_at(P, ts, weights=UNIT, n_terms=SERIES_TERMS), u3)
    q2s = np.array_equal(p81._end_at(P, ts, profiles=UNIFORM, n_terms=SERIES_TERMS), u3)
    gates["Q2-p"], gates["Q2-s"] = bool(q2p), bool(q2s)
    _say(comm, f"[Q2] (Q2-p) eddy(UNIT) == p75._u3 bitwise {q2p}; (Q2-s) end(UNIFORM) == p75._u3 "
               f"bitwise {q2s} (both asserted)")

    # (Rep-s): P80's six theta_s(t*) strings; P81's on step 6's own scalar calls
    th = tuple(f"{v:.6f}" for v in d80["E_abs"][-1])
    gates["Rep-s P80"] = th == p80.THETA_TSTAR
    pts = ((0.0, 0.0, 0.0), (0.9 * A_BOX, 0.0, 0.0), (0.9 * A_BOX, 0.9 * B_BOX, 0.9 * C_BOX))
    st = tuple(f"{float(p81._sep(*p, n_terms=SEP_TERMS) - T_INF):.6f}" for p in pts)
    tt = tuple(f"{float(p81._sep_t(x, 0.0, 0.0, 0.2 * T_C, n_terms=SEP_TERMS) - T_INF):.6f}"
               for x in (0.0, 0.9 * A_BOX))
    gates["Rep-s P81 steady"] = st == p81.STEP6_STEADY
    gates["Rep-s P81 t*"] = tt == p81.STEP6_TSTAR
    _say(comm, f"[Rep-s] P80 theta_s(t*) {th}: {th == p80.THETA_TSTAR}; P81 step 6 steady {st}: "
               f"{st == p81.STEP6_STEADY}, t* {tt}: {tt == p81.STEP6_TSTAR} (strings, asserted)")

    # (Lim)
    ss_e = np.asarray(robin_box_quadratic_steady_temperature(
        P[:, 0], P[:, 1], P[:, 2], half_widths=HALF, k=K, q=Q, h=H, weights=WEIGHTS, t_inf=0.0,
        n_terms=QUAD_TERMS), dtype=float)
    lim_e = float(np.max(np.abs(p80._eddy_at(P, 1000.0 * T_C) - ss_e) / np.abs(ss_e)))
    ss_n = np.asarray(robin_box_separable_steady_temperature(
        P[:, 0], P[:, 1], P[:, 2], half_widths=HALF, k=K, q=Q, h=H, profiles=PROFILES,
        t_inf=0.0, n_terms=SEP_TERMS), dtype=float)
    lim_n = float(np.max(np.abs(p81._end_at(P, 1000.0 * T_C) - ss_n) / np.abs(ss_n)))
    gates["Lim eddy"], gates["Lim end"] = lim_e <= B1_BAND, lim_n <= B1_BAND
    _say(comm, f"[Lim] transient(1000 T_C) vs steady, six probes: eddy {lim_e:.3e}, end "
               f"{lim_n:.3e} (band {B1_BAND:g}, asserted)")

    # (Lin)
    tp, tf, tu = (p81._end_at(P, ts, profiles=x) for x in (PROFILES, FLIPPED, UNIFORM))
    lin = float(np.max(np.abs(tp + tf - 2.0 * tu)) / np.max(np.abs(tu)))
    gates["Lin"] = lin <= B1_BAND
    _say(comm, f"[Lin] max|theta(PROFILES) + theta(FLIPPED) - 2 theta(UNIFORM)| / max theta"
               f"(UNIFORM) at t*: {lin:.3e} (band {B1_BAND:g}, asserted)")

    # (Ball)
    w_e = max(abs(p80._ball_cub(p, r) - p80._ball_mean(p, r)) / p80._ball_mean(p, r)
              for p in lat_pts)
    w_n = max(abs(p81._ball_cub_end(p, r) - p81._ball_mean_end(p, r)) / p81._ball_mean_end(p, r)
              for p in lat_pts)
    gates["Ball eddy"], gates["Ball end"] = w_e <= B1_BAND, w_n <= B1_BAND
    _say(comm, f"[Ball] |cub - closed| / closed over 18 lattice balls, r {r:.6f} m: eddy "
               f"{w_e:.3e}, end {w_n:.3e} (band {B1_BAND:g}, asserted)")

    # (K-lat): custody of entry 52's lines
    lat52, star = p80._parse_c52()
    labs_c6 = {lab for lab, _X in c6._lattice()}
    fmt_bad = [x for x in lat52 if f"{float(x[2]):.9e}" != x[2]]
    tok = {x[0]: x[2] for x in lat52}
    star_ok = (len(star) == 1 and star[0][0] == p80.C52_STAR
               and star[0][1] == tok.get(p80.C52_STAR))
    gates["K-lat n"] = len(lat52) == 18
    gates["K-lat labels"] = {x[0] for x in lat52} == labs_c6 and len({x[0] for x in lat52}) == 18
    gates["K-lat fmt"] = not fmt_bad
    gates["K-lat c*"] = bool(star_ok)
    _say(comm, f"[K-lat] {p80.C52_LOG.name}: {len(lat52)} lattice matches (18 asserted), labels == "
               f"c6._lattice() {gates['K-lat labels']}, .9e mismatches {len(fmt_bad)}, c* "
               f"{p80.C52_STAR!r} token custody {star_ok} (asserted)")

    # controls
    cu_e = [abs(float(d80["E"][-1, i] / d80["U3"][i]) - 1.0) for i in range(6)]
    cu_n = [abs(float(d81["S"][-1, i] / d80["U3"][i]) - 1.0) for i in range(6)]
    ce_n = [abs(float(d81["S"][-1, i] / d80["E"][-1, i]) - 1.0) for i in range(6)]
    gates["c-u eddy"] = all(v > INTERIOR_BAND for v in cu_e)
    gates["c-u end"] = all(v > INTERIOR_BAND for v in cu_n)
    gates["c-e end"] = all(v > INTERIOR_BAND for v in ce_n)
    _say(comm, f"[c-u] |eddy/U3 - 1| least {min(cu_e):.4f}, |end/U3 - 1| least {min(cu_n):.4f}; "
               f"[c-e] |end/eddy - 1| least {min(ce_n):.4f} (all > {INTERIOR_BAND:g}, asserted)")
    cb_e = [abs(p80._ball_mean(p, r, k5=3.0) / p80._ball_cub(p, r) - 1.0) for p in lat_pts]
    cb_n = [abs(p81._ball_mean_end(p, r, shell=True) / p81._ball_cub_end(p, r) - 1.0)
            for p in lat_pts]
    gates["c-ball eddy"], gates["c-ball end"] = min(cb_e) > B1_BAND, min(cb_n) > B1_BAND
    _say(comm, f"[c-ball] shell-for-ball miss, 18 points: eddy least {min(cb_e):.4e}, end least "
               f"{min(cb_n):.4e} (all > {B1_BAND:g}, asserted)")

    # printed: the six-probe table
    _say(comm, "[TABLE-6] probe | U3/P t* K/W | eddy/P t* | end/P t* | eddy/U3 | end/U3")
    for i, lab in enumerate(labs6):
        u, e, n = d80["U3"][i], d80["E"][-1, i], d81["S"][-1, i]
        _say(comm, f"[TABLE-6] {lab} | {u:.7g} | {e:.7g} | {n:.7g} | {e / u:.7g} | {n / u:.7g}")

    # printed: the lattice table with entry 52
    sar = {x[0]: float(x[2]) for x in lat52}
    ci = p80._cstar_idx()
    bm_e = np.array([p80._ball_mean(p, r) for p in lat_pts])
    bm_n = np.array([p81._ball_mean_end(p, r) for p in lat_pts])
    u_l, e_l, n_l = d80["U3_lat"], d80["L"][-1], d81["SL"][-1]
    _say(comm, "[TABLE-18] label | U3/P t* | eddy/P t* | end/P t* | eddy ball/c* | end ball/c* | "
               "entry 52 SAR10g/c*")
    for j, lab in enumerate(lat_labs):
        _say(comm, f"[TABLE-18] {lab} | {u_l[j]:.7g} | {e_l[j]:.7g} | {n_l[j]:.7g} | "
                   f"{bm_e[j] / bm_e[ci]:.4g} | {bm_n[j] / bm_n[ci]:.4g} | "
                   f"{sar[lab] / sar[p80.C52_STAR]:.4g}")
    col = {lab: j for j, lab in enumerate(lat_labs) if lab.startswith("x+0mm")}
    for tag, bm in (("eddy ball", bm_e), ("end ball", bm_n)):
        c50 = bm[col["x+0mm dz-50mm"]] / bm[col["x+0mm dz-250mm"]]
        c150 = bm[col["x+0mm dz-150mm"]] / bm[col["x+0mm dz-250mm"]]
        _say(comm, f"[COL] x+0mm {tag} mean, dz-50/dz-250 {c50:.4f}, dz-150/dz-250 {c150:.4f}")
    s_m = (sar["x+0mm dz-50mm"] / sar["x+0mm dz-250mm"],
           sar["x+0mm dz-150mm"] / sar["x+0mm dz-250mm"])
    s_p = (sar["x+0mm dz+50mm"] / sar["x+0mm dz+250mm"],
           sar["x+0mm dz+150mm"] / sar["x+0mm dz+250mm"])
    _say(comm, f"[COL] x+0mm entry 52 SAR10g: - side {s_m[0]:.4f} / {s_m[1]:.4f}, + side "
               f"{s_p[0]:.4f} / {s_p[1]:.4f} (printed only; none fitted)")
    for tag, v in (("U3", u_l), ("eddy", e_l), ("end", n_l)):
        tie = [lat_labs[j] for j in range(len(lat_labs)) if v[j] >= np.max(v) * (1.0 - 1e-4)]
        _say(comm, f"[TIE] {tag}/P t* tie set within 1e-4 rel of max: {tie} (no exact argmax)")


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example expects the complex DolfinX build the runner's `mat:` group "
            "sources (source /usr/local/bin/dolfinx-complex-mode)"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-192 -- entry 59's readout points under three closed-form SAR pictures "
               "(MAT-21 P75 / P80 / P81)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] imported: B1_BAND {B1_BAND:g}, INTERIOR_BAND {INTERIOR_BAND:.0e}; "
               f"T_STAR {p75.T_STAR:.4f} s; SERIES_TERMS {SERIES_TERMS}, QUAD_TERMS {QUAD_TERMS}, "
               f"SEP_TERMS {SEP_TERMS}")

    gates: dict = {}
    closed_form(comm, gates)

    msh = dmesh.create_box(
        comm, [np.array([-A_BOX, -B_BOX, -C_BOX]), np.array([A_BOX, B_BOX, C_BOX])],
        list(CELLS), cell_type=dmesh.CellType.tetrahedron)
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n_loc = imap.size_local
    ct = dmesh.meshtags(msh, tdim, np.arange(n_loc, dtype=np.int32),
                        np.full(n_loc, 1, dtype=np.int32))
    n_cells = imap.size_global
    write_setup_figure(
        msh,
        ct,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel, three source pictures"},
        title="mat:41 -- unwalled gel box, P1 fields at t*: y = 0 slice",
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        comm=comm,
    )
    _say(comm, f"[mesh] create_box {CELLS} tetrahedra: {n_cells} cells; ranks {comm.size}")

    tq = time.perf_counter()
    fields = {pic: _interp(msh, pic) for pic in PICS}
    _say(comm, f"[fields] three P1 interpolations: {time.perf_counter() - tq:.1f} s")

    P6 = p80._probe_arr()
    lat_pts = p80._data()["lat_pts"]
    P24 = np.vstack([P6, np.asarray(lat_pts, dtype=float)])
    # (V-pt): the 24 readout points are mesh vertices; one batch of 6, one batch of 18
    worst_v, ok_all = {}, True
    for pic, f in fields.items():
        got6, ok6 = _read(f, P6)
        got18, ok18 = _read(f, np.asarray(lat_pts, dtype=float))
        got = np.concatenate([got6, got18])
        ref = _series(pic, P24)
        rel = np.abs(got - ref) / np.abs(ref)
        worst_v[pic] = float(np.max(rel))
        ok_all = ok_all and ok6 and ok18
        gates[f"V-pt {pic}"] = worst_v[pic] <= B1_BAND
    gates["valid points"] = ok_all
    _say(comm, "[V-pt] 24 vertices, P1 field via owner route vs series, max rel: "
               + ", ".join(f"{p} {worst_v[p]:.3e}" for p in PICS)
               + f" (band {B1_BAND:g}, asserted)")

    # (cV): z-edge midpoints, the six probe points shifted by half the z spacing
    dz = C_BOX / CELLS[2] / 2.0
    Pm = P6.copy()
    Pm[:5, 2] = dz
    Pm[5, 2] = C_BOX - dz
    miss = {}
    for pic, f in fields.items():
        got, ok = _read(f, Pm)
        ref = _series(pic, Pm)
        rel = np.abs(got - ref) / np.abs(ref)
        miss[pic] = rel
        gates[f"cV {pic}"] = bool(np.all(rel > B1_BAND)) and ok
        _say(comm, f"[cV] {pic}: P1 vs series at six z-edge midpoints, rel miss "
                   + ", ".join(f"{v:.3e}" for v in rel)
                   + f" (all > {B1_BAND:g} asserted; no ceiling, P1 is not exact on a curved theta)")

    imag = 0.0
    for f in fields.values():
        n_own = f.function_space.dofmap.index_map.size_local
        loc = float(np.max(np.abs(np.imag(f.x.array[:n_own])))) if n_own else 0.0
        imag = max(imag, comm.allreduce(loc, op=MPI.MAX))
    gates["max_imag_zero"] = imag == 0.0
    _say(comm, f"max|Im theta| = {imag:.3e}")

    outs = {f.name: f for f in fields.values()}
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, ct, outs, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[output] {path} -- theta(t*) per deposited watt: uniform, eddy, end-weighted")

    bad = {k: v for k, v in gates.items() if not v}
    _say(comm, f"[gates] {len(gates)} collected, {len(bad)} red: {bad}")
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    assert not bad, bad
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
