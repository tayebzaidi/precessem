"""Example (`EX-125`): the Pennes sink w(T - T_inf) on the coil-driven gel phantom.

`mat:6` / `mat:7` put the sink under a *uniform* source on a ball; `mat:9`, `mat:15`
and `mat:18` heat this phantom from the coil with no sink.  This example puts the
sink under the **coil-driven** source: `MAT-15` step 4's own body
(``_run_perfused``), on `ANS-2`'s F-small at 10 MHz, sigma = 800 S/m, four single-leg
drives (G = 1 W / P_1), P1 thermal, gel constants, at w = 24 / 240 / 2400
W/m^3/K and the limit value w_eps = 1e-8.  **Gel has no perfusion**: every w is a
numerical sink term, not a gel or tissue property.

Nothing restated.  ``_run_perfused(return_fields=True)`` is the gate's body with
one additive keyword (rule (a) lift, EX-125: it also returns ``base`` and the CN
``results`` keyed (k, w)); every band and helper is imported by name from the gate
module.  The example adds the heat budget at t*, the solver's ``perfused_power``
against the gate's integral, and ParaView fields.

Asserted, each beside ``20261008T140322Z_MAT-15-step4-W1.log``: (F) fixture, (E')
with the sink, (E), (L0), (Mo), (C4) at the gate's bands; (R) the returned runs are
the gate's runs; (Q-w) the solver's ``perfused_power`` equals the gate's
integral; the records; negative controls (d), (P-ctl), (L0-ctl).

Scope: one fixture (untuned F-small, 10 MHz, 1 W normalisation); no gel or tissue
perfusion value, no Pennes closed form, no tuned / 64 MHz, no compliance, no
external check.

Run through the harness (the ``mat:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mat:22 -n 2 -t 590

Output: ``examples/materials/paraview_output/materials_22_pennes_sink_coil_driven_phantom_combined.xdmf``
and ``materials_22_pennes_sink_coil_driven_phantom_sink_table.csv``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

from fem_em_solver.io.paraview_utils import adopt_host_ownership, write_xdmf_with_tags
from fem_em_solver.post.setup_figure import write_setup_figure

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.mesh.test_birdcage_port_tags import LEG_COUNT  # noqa: E402
from tests.validation import test_birdcage_b1_plus_map as b1map  # noqa: E402
from tests.validation import test_mat10_sar_to_heat_sphere as m10  # noqa: E402
from tests.validation import test_mat11_perfusion_sphere as m11  # noqa: E402
from tests.validation import test_mat12_transient_sar_heating_sphere as m12  # noqa: E402
from tests.validation import test_mat15_coil_driven_heating_phantom as m15  # noqa: E402
from tests.validation import test_mat15_pennes_coil_driven as m15s4  # noqa: E402
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    PHANTOM_CELL_TAG,
    V_BAND,
)
from tests.validation.test_mat15_pennes_coil_driven import (  # noqa: E402
    B_BAND,
    C4_COVARIANCE_BAND,
    D_CONTROL_FLOOR,
    E_BAND,
    MO_GAP_FLOOR,
    N_MAIN,
    ROBIN_TAG,
    T_INF,
    T_STAR,
    W_EPS,
    W_VALUES,
    _allsum,
    _e_prime_w,
    _run_perfused,
)
from tests.validation.test_port_birdcage_lumped_column import CONDUCTOR_CELL_TAG  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "materials_22_pennes_sink_coil_driven_phantom"
W1 = "20261008T140322Z_MAT-15-step4-W1.log"
SUB_CELLS = 58866


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `mat:` "
            "group sources it automatically"
        )

    structural = {
        "complex build": True,
        "comm.size == 2": comm.size == 2,
        "m15s4.E_BAND is m12.E_BAND": m15s4.E_BAND is m12.E_BAND,
        "m15s4.B_BAND is m10.B_BAND": m15s4.B_BAND is m10.B_BAND,
        "m15s4.D_CONTROL_FLOOR is m15.D_CONTROL_FLOOR":
            m15s4.D_CONTROL_FLOOR is m15.D_CONTROL_FLOOR,
        "m15s4._run_heated is m15._run_heated": m15s4._run_heated is m15._run_heated,
        "m15s4.T_STAR is m15.T_STAR": m15s4.T_STAR is m15.T_STAR,
        "m15s4.W is m11.W": m15s4.W is m11.W,
        "m15s4.W_STRONG is m11.W_STRONG": m15s4.W_STRONG is m11.W_STRONG,
        "m15s4.C4_COVARIANCE_BAND is b1map.C4_COVARIANCE_BAND":
            m15s4.C4_COVARIANCE_BAND is b1map.C4_COVARIANCE_BAND,
        "W_VALUES == (W/10, W, W_STRONG)":
            m15s4.W_VALUES == (m15s4.W / 10.0, m15s4.W, m15s4.W_STRONG),
        "W_EPS > 0": m15s4.W_EPS > 0.0,
    }
    for name, ok in structural.items():
        assert ok, f"structural import check failed: {name}"
    if rank0:
        print("=" * 72)
        print("EX-125 -- the Pennes sink on the coil-driven gel phantom (MAT-15 step 4, re-executed)")
        print("=" * 72)
        print("[structural] " + "; ".join(f"{k}: {v}" for k, v in structural.items()))
        print(f"[records] every anchor below is beside {W1}", flush=True)

    # ---- the gate's own body, with its fields (called once) -----------------
    t_all = time.perf_counter()
    r = _run_perfused(return_fields=True)
    t_gate = time.perf_counter() - t_all
    base, results = r["base"], r["results"]
    sub, ft, G, P = base["sub"], base["ft"], base["G"], base["P"]
    v_sub = base["v_sub"]
    msh, cell_tags = base["msh"], base["cell_tags"]
    ds = ufl.Measure("ds", domain=sub, subdomain_data=ft)(ROBIN_TAG)  # gate :120

    # EX-57 setup figure (after the fixture, which builds the mesh); no-op unless
    # FEM_EM_SETUP_FIGURES=1.
    t1 = time.perf_counter()
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            CONDUCTOR_CELL_TAG: "conductor",
            2: "air",
            PHANTOM_CELL_TAG: "phantom (gel)",
            **{PORT_LOWER + i: f"leg L{i} port lower" for i in range(1, LEG_COUNT + 1)},
            **{PORT_UPPER + i: f"leg L{i} port upper" for i in range(1, LEG_COUNT + 1)},
        },
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="mat:22 -- the Pennes sink w(T - T_inf) on the coil-driven gel phantom "
              "(MAT-15 step 4; w a numerical sink)",
        comm=comm,
    )
    t_fig = time.perf_counter() - t1

    def _mean_theta(T):  # gate :150-151
        return float(np.real(_allsum(comm, (T - T_INF) * ufl.dx(domain=sub)))) / v_sub

    ks = range(4)
    seq = (0.0,) + tuple(W_VALUES)

    # ---- (R) and (Q-w): one pass over the 12 (k, w) -------------------------
    t1 = time.perf_counter()
    budget = {}   # (k, w) -> (rel, stor, conv, perf)
    r_mean = {}   # (k, w) -> |<theta>_returned / <theta>_gate - 1|
    r_ep = {}
    q_rel = {}
    for k in ks:
        dep = G * P[k]
        for w in W_VALUES:
            res = results[(k, w)]
            r_mean[(k, w)] = abs(_mean_theta(res.temperature) / r["mean"][(k, w)] - 1.0)
            ep, stor, conv, perf = _e_prime_w(sub, ds, res, dep, w)
            r_ep[(k, w)] = ep
            budget[(k, w)] = (ep, stor, conv, perf)
            q_rel[(k, w)] = abs(perf - res.perfused_power) / abs(res.perfused_power)
        r_mean[(k, W_EPS)] = abs(
            _mean_theta(results[(k, W_EPS)].temperature) / r["mean"][(k, 0.0)] - 1.0)
    t_rq = time.perf_counter() - t1

    # ---- negative control (L0-ctl) ------------------------------------------
    l0ctl = {k: abs(r["mean"][(k, W_VALUES[0])] / r["mean"][(k, 0.0)] - 1.0) for k in ks}

    # ---- nodal peak theta, k = 0 (owned dofs, reduced) ----------------------
    def _peak(T):
        n_own = T.function_space.dofmap.index_map.size_local * T.function_space.dofmap.index_map_bs
        loc = float(np.max(np.real(T.x.array[:n_own]) - T_INF)) if n_own else -np.inf
        return float(comm.allreduce(loc, op=MPI.MAX))

    peak = {0.0: _peak(base["T_star"][0])}
    for w in W_VALUES:
        peak[w] = _peak(results[(0, w)].temperature)
    peak_ratio = peak[W_VALUES[-1]] / peak[0.0]
    removed_mean = r["mean"][(0, 0.0)] - r["mean"][(0, W_VALUES[-1])]

    # ---- records (asserted at the gate's printed digits) --------------------
    REC = []  # (name, got, want, tol_rel, W1 line)

    def rec(name, got, want, tol, line):
        REC.append((name, float(got), float(want), tol, line))

    rec("V_sub [m^3]", v_sub, 2.259978941111e-04, 1e-11, ":1926")
    for k, want in zip(ks, (5.503204728587e-08, 5.498192946991e-08,
                            5.500091017251e-08, 5.496862018702e-08)):
        rec(f"P_{k + 1} [W]", P[k], want, 1e-11, ":1925")
    for k, want, line in zip(ks, (1.153313650151, 1.152271458961, 1.152660000060,
                                  1.151992598349), (":1957", ":1965", ":1973", ":1981")):
        rec(f"<theta>_P1(t*) k={k} [K]", r["mean"][(k, 0.0)], want, 1e-11, line)
    for w, want, tol, line in zip(W_VALUES, (1.149352216, 1.114522911, 0.835062640),
                                  (1e-8, 1e-8, 1e-7), (":1960", ":1962", ":1964")):
        rec(f"<theta>(t*) k=0 w={w:g} [K]", r["mean"][(0, w)], want, tol, line)
    for w, want, line in zip(seq, (2.587441071, 2.579306150, 2.507705061, 1.927156653),
                             (":1957", ":1960", ":1962", ":1964")):
        rec(f"10 g theta k=0 w={w:g} [K]", r["th"][(0, w)], want, 1e-8, line)
    for w, (s_w, c_w, p_w), line in zip(
            W_VALUES,
            ((7.603747991358e-01, 2.334596062651e-01, 6.165594599057e-03),
             (7.132479311534e-01, 2.269428325893e-01, 5.980923625733e-02),
             (3.761836210292e-01, 1.742682763548e-01, 4.495481026160e-01)),
            (":1959", ":1961", ":1963")):
        _, stor, conv, _ = budget[(0, w)]
        rec(f"stor k=0 w={w:g} [W]", stor, s_w, 1e-11, line)
        rec(f"conv k=0 w={w:g} [W]", conv, c_w, 1e-11, line)
        rec(f"perfused_power k=0 w={w:g} [W]", results[(0, w)].perfused_power, p_w, 1e-11, line)
    for w, want, line in zip(W_VALUES, (0.996565, 0.966366, 0.724055),
                             (":2002", ":2003", ":2004")):
        ratio = float(np.mean([r["mean"][(k, w)] / r["mean"][(k, 0.0)] for k in ks]))
        rec(f"drive-mean <theta>/P1 w={w:g}", ratio, want, 1e-4, line)
    rec("(P-ctl) k=0", r["pctl"][0], 0.059809, 1e-3, ":1962")
    rec("(d)", r["ctl_d"], 1.000000, 1e-5, ":2005")

    n_sub = sub.topology.index_map(3).size_global

    if rank0:
        print(f"\n[anchor] (F) fixture: sub cells {n_sub} (W1 {SUB_CELLS}, :1926); "
              f"(B) {max(base['B'][k] for k in ks):.3g} <= {B_BAND:.0e}; "
              f"(E') step 1 {max(base['Ep'][k] for k in ks):.3g} <= {E_BAND:.0e}; "
              f"V_sub rel {base['err_V']:.3g} <= {V_BAND:.0e}  (W1 :1926-1945)")
        print(f"[anchor] (E') with the sink: worst over 12 {max(r['Ep'].values()):.4g} "
              f"<= {E_BAND:.0e}  (W1 <= 1.623e-14, :1985)")
        wk = max(r["E"], key=r["E"].get)
        print(f"[anchor] (E) worst over {len(r['E'])} runs {r['E'][wk]:.4g} on {wk} "
              f"<= {E_BAND:.0e}  (W1 8.659e-14, :2006)")
        print("[anchor] (L0) w_eps vs P1: " + "  ".join(f"k{k} {r['L0'][k]:.4g}" for k in ks)
              + f" <= {B_BAND:.0e}  (W1 1.506e-12 ... 1.526e-12, :1958-1982)")
        print(f"[anchor] (Mo) smallest gap {min(min(g) for g in r['mo_gaps'].values()):.4g} "
              f"> {MO_GAP_FLOOR:.0e}  (W1 :1989-1996)")
        print(f"[anchor] (C4) worst pair {max(r['pairs'].values()) * 100:.4g} % "
              f"<= {C4_COVARIANCE_BAND * 100:.1f} %  (W1 <= 0.0489 %, :1998-2000)")
        print(f"[anchor] (R) returned runs are the gate's: <theta> rel worst "
              f"{max(r_mean.values()):.4g} <= {B_BAND:.0e} (predicted 0); (E') recomputed worst "
              f"{max(r_ep.values()):.4g} <= {E_BAND:.0e}, gate's worst "
              f"{max(r['Ep'].values()):.4g}")
        print(f"[anchor] (Q-w) perfused_power vs the gate's integral, worst of 12: "
              f"{max(q_rel.values()):.4g} <= {E_BAND:.0e}  (predicted <= 1e-14; W1 equal at 13 "
              f"digits, :1959-1963)")
        print(f"[control] (d) {r['ctl_d']:.6f} > {D_CONTROL_FLOOR:.0e}  (W1 1.000000, :2005)")
        print("[control] (P-ctl) " + "  ".join(f"k{k} {r['pctl'][k]:.6f}" for k in ks)
              + f" > {D_CONTROL_FLOOR:.0e}  (W1 0.059809-0.059810, :1962-1986; predicted "
                "[0.04, 0.08])")
        print("[control] (L0-ctl) w=24 vs P1 "
              + "  ".join(f"k{k} {l0ctl[k]:.6e}" for k in ks)
              + f" > {B_BAND:.0e}; factor over band "
              f"{min(l0ctl.values()) / B_BAND:.3g}x  (W1 3.434828e-03 ... 3.434838e-03, :1989-1995)")
        print("\n[records] name | got | W1 | rel | tol | W1 line")
        for name, got, want, tol, line in REC:
            print(f"[record] {name} | {got:.6g} | {want:.6g} | {abs(got / want - 1.0):.3g} | "
                  f"{tol:.0e} | {W1}{line}")
        print("\n[budget] share of G*P_k at t*: storage / Robin wall / perfusion (%), "
              "w = 0 row is step 1's own split 76.5800 / 23.4200 / 0 (W1 :1930)")
        for k in ks:
            dep = G * P[k]
            for w in W_VALUES:
                _, stor, conv, perf = budget[(k, w)]
                print(f"[budget] k={k} w={w:g}: {100 * stor / dep:.6g} / {100 * conv / dep:.6g} / "
                      f"{100 * perf / dep:.6g}")
        print("[printed] k=0 predicted perfusion 0.616559 / 5.98092 / 44.9548 %, storage "
              "76.0375 / 71.3248 / 37.6184 % (records)")
        print("[printed] <theta>(t*)/P1 per w, drive mean: "
              + " / ".join(f"{float(np.mean([r['mean'][(k, w)] / r['mean'][(k, 0.0)] for k in ks])):.6f}"
                           for w in W_VALUES)
              + "  (MAT-15.md:279 lumped predictions 0.996-0.998 / 0.955-0.975 / 0.68-0.77)")
        print("[printed] nodal peak theta k=0 (K): "
              + "  ".join(f"w={w:g} {peak[w]:.6g}" for w in seq)
              + f"; theta_max({W_VALUES[-1]:g})/theta_max(0) = {peak_ratio:.6g} "
                "(predicted [0.70, 0.80]; unlogged drafter estimate)")
        print(f"[printed] heat removed <theta_P1> - <theta_{W_VALUES[-1]:g}> at k=0 = "
              f"{removed_mean:.6g} K (predicted 0.318251)")
        print(f"[printed] max|Im theta| = {r['max_imag']:.3e} (predicted 0, W1 :2006)")
        print(f"[printed] wall: _run_perfused {t_gate:.1f} s; (R)/(Q-w) {t_rq:.1f} s; "
              f"setup figure {t_fig:.1f} s", flush=True)

    # ---- asserts (every rank) -------------------------------------------------
    assert comm.size == 2
    assert n_sub == SUB_CELLS, f"(F) sub cells {n_sub}"
    for k in ks:
        assert base["B"][k] <= B_BAND, f"(F) B drive {k}"
        assert base["Ep"][k] <= E_BAND, f"(F) E' step 1 drive {k}"
    assert base["err_V"] <= V_BAND, "(F) V_sub"
    for k in ks:
        for w in W_VALUES:
            assert r["Ep"][(k, w)] <= E_BAND, f"(E') k={k} w={w:g}"
    assert len(r["E"]) == 16
    for name, v in r["E"].items():
        assert v <= E_BAND, f"(E) {name}: {v:.3e}"
    for k in ks:
        assert r["L0"][k] <= B_BAND, f"(L0) k={k}: {r['L0'][k]:.3e}"
        for key in ("mean", "th"):
            for g in r["mo_gaps"][(k, key)]:
                assert g > MO_GAP_FLOOR, f"(Mo) k={k} {key}: {g:.3e}"
    for w in W_VALUES:
        for k in ks:
            assert r["pairs"][(w, k)] <= C4_COVARIANCE_BAND, f"(C4) w={w:g} k={k}"
    for key, v in r_mean.items():
        assert v <= B_BAND, f"(R) <theta> {key}: {v:.3e}"
    for key, v in r_ep.items():
        assert v <= E_BAND, f"(R) E' {key}: {v:.3e}"
    for key, v in q_rel.items():
        assert v <= E_BAND, f"(Q-w) {key}: {v:.3e}"
    assert r["ctl_d"] > D_CONTROL_FLOOR, f"(d) did not miss: {r['ctl_d']:.3e}"
    for k in ks:
        assert r["pctl"][k] > D_CONTROL_FLOOR, f"(P-ctl) k={k} did not miss"
        assert l0ctl[k] > B_BAND, f"(L0-ctl) k={k} did not separate"
    for name, got, want, tol, line in REC:
        assert abs(got / want - 1.0) <= tol, f"record {name}: {got!r} vs {want!r} ({W1}{line})"

    # ---- outputs -------------------------------------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    V1 = results[(0, W_VALUES[0])].temperature.function_space
    w_hi = W_VALUES[-1]

    def _field(name, arr):
        f = fem.Function(V1, name=name)
        f.x.array[:] = np.real(arr)
        f.x.scatter_forward()
        return f

    th0 = np.real(base["T_star"][0].x.array) - T_INF
    th_by_w = {w: np.real(results[(0, w)].temperature.x.array) - T_INF for w in W_VALUES}
    fields = {"theta_w0": _field("theta_w0", th0)}
    for w in W_VALUES:
        fields[f"theta_w{w:g}"] = _field(f"theta_w{w:g}", th_by_w[w])
    fields[f"theta_removed_w{w_hi:g}"] = _field(f"theta_removed_w{w_hi:g}", th0 - th_by_w[w_hi])
    fields[f"sink_density_w{w_hi:g}"] = _field(f"sink_density_w{w_hi:g}", w_hi * th_by_w[w_hi])
    q_k0 = base["q_h"][0]
    q_k0.name = "q_h_k0"  # the XDMF dataset name is the Function's, not the dict key
    fields["q_h_k0"] = q_k0
    t1 = time.perf_counter()
    write_xdmf_with_tags(OUTPUT_DIR / f"{BASENAME}_combined", sub, None, fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    t_xdmf = time.perf_counter() - t1

    if rank0:
        csv_name = f"{BASENAME}_sink_table.csv"
        with open(OUTPUT_DIR / csv_name, "w") as fh:
            fh.write("k,w,mean_theta_K,ratio_to_P1,theta_10g_K,ratio_10g,storage_share,"
                     "convection_share,perfusion_share,perfused_power_W\n")
            for k in ks:
                dep = G * P[k]
                for w in seq:
                    mean_w, th_w = r["mean"][(k, w)], r["th"][(k, w)]
                    head = (f"{k},{w:g},{mean_w!r},{mean_w / r['mean'][(k, 0.0)]!r},"
                            f"{th_w!r},{th_w / r['th'][(k, 0.0)]!r}")
                    if w == 0.0:
                        fh.write(head + ",,,,\n")
                    else:
                        _, stor, conv, perf = budget[(k, w)]
                        fh.write(head + f",{stor / dep!r},{conv / dep!r},{perf / dep!r},"
                                 f"{results[(k, w)].perfused_power!r}\n")
        adopt_host_ownership(OUTPUT_DIR, comm=MPI.COMM_SELF)
        print(f"\n[output] {BASENAME}_combined.xdmf: P1 " + ", ".join(
            n for n in fields if n != "q_h_k0") + f" (K, or W/m^3 per watt, drive k=0, t* = "
              f"{T_STAR:.1f} s) and DG0 q_h_k0 (W/m^3 at G = 1 W/P_1); {t_xdmf:.1f} s")
        print(f"[output] {csv_name}: {4 * len(seq)} rows")
        print("EX-125 anchors, records and controls: All assertions hold")


if __name__ == "__main__":
    main()
