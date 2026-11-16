"""Example (`EX-116`): the three shipped presets through ``femem`` -- all three checked
and priced by ``femem run --estimate`` at one width, the two degree-1 presets run, the
degree-2 one priced but never solved.

`WF-4` shipped ``examples/presets/debug.yaml``, ``dev.yaml`` and ``benchmark-lite.yaml``
and gated that each reproduces its record through ``femem run`` and that each price lands
inside ``BRACKET`` of its own run, one preset per window, ``benchmark-lite`` at ``-n 8``.
No example read them.  This one reads them in place (no copy) at ``-n 2`` and shows the
user's workflow: ``femem check`` and ``femem run --estimate`` on all three, the price
ladder, then ``femem run`` on the two that are cheap at this width.  ``benchmark-lite``
stays priced and unsolved: its price at ``-n 2`` is the width law applied to the same
fitted row set (``Study.estimate`` divides the fitted ``time * ranks`` by the width).

Asserted, every band imported from the gate module (the ANS-1 rule):

* (Pd) per preset: the printed estimate equals ``load_study(p).estimate(comm)`` field by
  field (``g._diff``), rc (0, 0), no estimate directory, and the solve wrapper is never
  called by any ``--estimate``;
* (R-debug), (R-dev): the gate's reproductions (mean SAR vs ``MAT4_DEG1_RECORD`` within
  ``L0_BAND``; S vs ``_record(1)`` within ``LEG_D0_REPRODUCTION_BAND``, reciprocity,
  passivity, circulant spreads);
* (P-debug), (P-dev): predicted / measured seconds inside ``[1/BRACKET, BRACKET]``;
* (E) the estimate builds the mesh the run solves (``est.cells == summary["cells"]``;
  ``benchmark-lite``: the gate's hole-cell record and 570 332 unknowns);
* (W) a price is per width; (L) the ladder is strictly ordered.

Negative controls (asserted): (c1) ``debug`` vs the degree-2 SAR record; (c2) ``dev`` vs
``_record(2)``; (cL) ``benchmark-lite``'s price against ``debug``'s run is outside the
bracket (an estimator returning one constant would pass it); (cW) the ``-n 8`` price
applied at ``-n 2`` unscaled misses the width law.

Per the 2026-10-15 18:00 ruling, a live price is never pinned: ``PRICE_TABLE`` has moved
since `WF-4`'s logs (a third degree-1 row, `U55a`).  The `WF-4` prices are printed beside
the live ones, not asserted.  Nothing prints the S or Z matrices.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:22 -n 2 -t 590

Not claimed: a ``--preset`` flag or presets in the wheel; the price at any width other
than the one run; a memory claim; ``benchmark-lite``'s S at ``-n 2``; an external check.
"""

from __future__ import annotations

import contextlib
import dataclasses
import io
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.study import ResourceEstimate, _predict  # noqa: E402
from fem_em_solver.utils.instrumentation import report_peak_rss  # noqa: E402

import tests.package.test_wf4_presets as g  # noqa: E402
from tests.validation.test_birdcage_sar_mass_averaged import EXACT_IDENTITY_RTOL  # noqa: E402,F401

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "paraview_output"
FIGURE_DIR = HERE / "figures"
BASENAME = "mri_22_presets_priced_then_run"
AIR_TAG = 2

# WF-4's printed prices at -n 2 (debug, dev) and -n 8 (benchmark-lite), printed beside the
# live prices, never asserted (the 2026-10-15 ruling: PRICE_TABLE moved with U55a).
REC = {
    # docs/testing/logs/20261015T202725Z_WF-4-s1b-W1-debug.log:284
    "price_debug_s": 13.70,
    # docs/testing/logs/20261015T202754Z_WF-4-s1b-W2-dev.log:5404
    "price_dev_s": 29.77,
    # docs/testing/logs/20261015T202931Z_WF-4-s1b-W3-benchmark-lite.log:2691
    "price_bl_n8_s": 76.65,
    "bl_unknowns": 570332,
}
PRESETS = (("debug", g.DEBUG), ("dev", g.DEV), ("benchmark-lite", g.BL))


def _fresh(path: Path, comm) -> None:
    if comm.rank == 0 and path.exists():
        shutil.rmtree(path)
    comm.Barrier()


def _listing(path: Path):
    return sorted((p.name, p.stat().st_size) for p in path.iterdir())


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- structural asserts -------------------------------------------------
    assert comm.size == 2, f"mri:22 runs at -n 2; got -n {comm.size}"
    if not np.issubdtype(default_scalar_type, np.complexfloating):
        raise RuntimeError("complex build required (source dolfinx-complex-mode)")
    assert g.load_study(g.BL).degree == 2
    assert g.load_study(g.DEBUG).degree == 1 and g.load_study(g.DEV).degree == 1
    assert g.ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH == 8
    names = [f.name for f in dataclasses.fields(ResourceEstimate)]
    OUTPUT_DIR.mkdir(exist_ok=True)

    # ---- (1) check + price all three; the solve wrapper must stay at zero ----------
    captured = []
    original = g.study_yaml.run_study_file

    def _call_through(*args, **kwargs):
        out = original(*args, **kwargs)
        captured.append(out)
        return out

    rows = {}
    g.study_yaml.run_study_file = _call_through
    try:
        for label, path in PRESETS:
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                rc_check = g.cli.main(["check", str(path)])
            check_all = g._rc_all(rc_check, comm)
            family = next(
                (ln for ln in buf.getvalue().splitlines() if ln.startswith("family:")), "")
            family = comm.bcast(family, root=0)

            est_dir = OUTPUT_DIR / f"mri_22_{label}_estimate"
            _fresh(est_dir, comm)
            api = g.load_study(path).estimate(comm)
            buf = io.StringIO()
            t0 = time.perf_counter()
            with contextlib.redirect_stdout(buf):
                rc_est = g.cli.main(["run", "--estimate", str(path), "--output-dir", str(est_dir)])
            t_est = time.perf_counter() - t0
            est_all = g._rc_all(rc_est, comm)
            verdict = None
            if comm.rank == 0:
                parsed = g._parse(buf.getvalue(), names)
                complete = set(parsed) == set(names) and all(len(v) == 1 for v in parsed.values())
                verdict = dict(complete=complete,
                               differ=sorted(g._diff(parsed, api)) if complete else None)
            verdict = comm.bcast(verdict, root=0)
            exists = comm.allreduce(bool(est_dir.exists()), op=MPI.LOR)
            assert check_all == (0, 0), (label, check_all)
            assert est_all == (0, 0), (label, est_all)
            assert verdict["complete"] and verdict["differ"] == [], (label, verdict)
            assert not exists, f"(Pd) {label}: --estimate wrote a directory"
            rows[label] = dict(est=api, family=family, t_est=t_est, check=check_all)
        n_solves_estimating = len(captured)
        assert n_solves_estimating == 0, f"(Pd) --estimate reached a solve ({n_solves_estimating})"

        # ---- (2) run the two degree-1 presets -------------------------------------
        runs = {}
        for label, path in PRESETS[:2]:
            out = OUTPUT_DIR / f"mri_22_{label}_run"
            _fresh(out, comm)
            before = len(captured)
            t0 = time.perf_counter()
            rc_run = g.cli.main(["run", str(path), "--output-dir", str(out)])
            t_cli = time.perf_counter() - t0
            assert g._rc_all(rc_run, comm) == (0, 0), (label, rc_run)
            assert len(captured) == before + 1, (label, len(captured), before)
            run = captured[-1]
            comm.Barrier()
            summary = None
            if comm.rank == 0:
                with open(out / "summary.json") as fh:
                    summary = json.load(fh)
            summary = comm.bcast(summary, root=0)
            listing = _listing(out) if comm.rank == 0 else None
            peak = report_peak_rss(comm, label=f"[mri:22] {label} summed peak ru_maxrss")
            runs[label] = dict(run=run, summary=summary, t_cli=t_cli, peak=peak, listing=listing)
    finally:
        g.study_yaml.run_study_file = original
    assert g.study_yaml.run_study_file is original

    est_d, est_v, est_b = (rows[k]["est"] for k in ("debug", "dev", "benchmark-lite"))
    sum_d, sum_v = runs["debug"]["summary"], runs["dev"]["summary"]
    run_d = runs["debug"]["run"]

    # ---- EX-57 setup figure on the debug run's mesh (no-op unless flagged) ---------
    msh0 = run_d.result.mesh
    sphere_tag = run_d.study.geometry.sphere_tag
    FIGURE_DIR.mkdir(exist_ok=True)
    write_setup_figure(
        msh0, run_d.result.cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={sphere_tag: "gel sphere", AIR_TAG: "air box (imposed-field walls)"},
        hide_tags=(AIR_TAG,), translucent_tags=(sphere_tag,),
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        title=("mri:22 -- the debug preset's lossy sphere, first of three presets "
               "priced by femem run --estimate (WF-4)"),
        comm=comm,
    )

    # ---- anchors --------------------------------------------------------------------
    sar = sum_d["mean_sar"]["mean_sar_w_per_kg"]
    r_debug = abs(sar - g.MAT4_DEG1_RECORD) / g.MAT4_DEG1_RECORD
    c1 = abs(sar / g.WF15_DEG2_RECORD - 1.0)

    s = g._s(sum_v)
    z = np.array([[complex(*p) for p in row] for row in sum_v["z_matrix"]])
    r_dev = g._rel(s, g._record(1))
    c2 = g._rel(s, g._record(2))
    recip = float(np.linalg.norm(s - s.T) / np.linalg.norm(s))
    sigma = np.linalg.svd(s, compute_uv=False)
    power = np.sum(np.abs(s) ** 2, axis=0)
    spreads = {n: g._class_spread(v) for n, v in g._circulant_classes(z).items()}

    p_debug = est_d.predicted_seconds / sum_d["elapsed_s"]
    p_dev = est_v.predicted_seconds / sum_v["elapsed_s"]
    cl = est_b.predicted_seconds / sum_d["elapsed_s"]

    # (W) per width: the -n 8 price times 8 equals the -n 2 price times 2 (same fit)
    price8 = _predict(est_b.unknowns, 2, int(g.ANS6_EXACT_STRIP_P2_RECORD_RANK_WIDTH))[1]
    w_rel = abs(price8 * 8 / (est_b.predicted_seconds * est_b.ranks) - 1.0)
    cw_rel = abs(price8 / est_b.predicted_seconds - 1.0)
    ladder = (est_d.predicted_seconds, est_v.predicted_seconds, est_b.predicted_seconds)

    if comm.rank == 0:
        print("=" * 72)
        print("EX-116 -- the three shipped presets priced, then two run (mri:22)")
        print("=" * 72)
        for label, _ in PRESETS:
            r = rows[label]
            e = r["est"]
            print(f"[check] {label}: rc {r['check']}, {r['family']!r}; --estimate rc (0, 0), "
                  f"wall {r['t_est']:.1f} s, no output directory (Pd)")
        print(f"[(Pd)] wrapper calls across the three --estimate: {n_solves_estimating}")
        print("[price table] preset degree cells unknowns rows_used predicted_s peak_GiB "
              "| measured_s ratio peak_RSS_GiB")
        for label, _ in PRESETS:
            e = rows[label]["est"]
            meas = ""
            if label in runs:
                sm = runs[label]["summary"]
                pk = (runs[label]["peak"] or 0) / 2**30
                ratio = e.predicted_seconds / sm["elapsed_s"]
                meas = f"| {sm['elapsed_s']:.4g} {ratio:.4g} {pk:.4g}"
            print(f"  {label} {e.degree} {e.cells} {e.unknowns} {e.rows_used} "
                  f"{e.predicted_seconds:.4g} {e.predicted_peak_gib:.4g} {meas}")
        print(f"[WF-4 prices, printed beside] debug {REC['price_debug_s']:.2f} s (-n 2), "
              f"dev {REC['price_dev_s']:.2f} s (-n 2), benchmark-lite {REC['price_bl_n8_s']:.2f} s "
              f"(-n 8); live table has moved (U55a), not asserted")
        print(f"[(R-debug)] SAR {sar!r} rel MAT4_DEG1_RECORD {r_debug:.3e} (L0_BAND "
              f"{g.L0_BAND:.0e}); degree {sum_d['degree']}, ranks {sum_d['ranks']}, "
              f"cells {sum_d['cells']}; wall {runs['debug']['t_cli']:.1f} s")
        print(f"[(R-dev)] ||S - P1||/||P1|| {r_dev:.3e} (band "
              f"{g.LEG_D0_REPRODUCTION_BAND:.0e}); reciprocity {recip:.3e} (band "
              f"{g.RECIPROCITY_BAND:.0e}); sigma_max {float(np.max(sigma)):.12f}; max column "
              f"power {float(np.max(power)):.12f}; spreads "
              + ", ".join(f"{k} {v * 100:.4f}%" for k, v in spreads.items())
              + f"; cells {sum_v['cells']}; wall {runs['dev']['t_cli']:.1f} s")
        print(f"[(P)] debug {p_debug:.3f}, dev {p_dev:.3f} (BRACKET {g.BRACKET})")
        print(f"[(E)] debug cells est {est_d.cells} run {sum_d['cells']}; dev est {est_v.cells} "
              f"run {sum_v['cells']}; benchmark-lite est cells {est_b.cells} (record "
              f"{g.ANS6_EXACT_STRIP_HOLE_CELL_RECORD}), unknowns {est_b.unknowns} "
              f"(record {REC['bl_unknowns']})")
        print(f"[(W)] benchmark-lite price at -n {est_b.ranks} {est_b.predicted_seconds:.4g} s; "
              f"same fit at -n 8 {price8:.4g} s; |t8*8 / (t2*2) - 1| = {w_rel:.3e}")
        print(f"[(L)] ladder predicted_s {ladder[0]:.4g} < {ladder[1]:.4g} < {ladder[2]:.4g}; "
              f"ratios dev/debug {ladder[1] / ladder[0]:.4g}, bl/dev {ladder[2] / ladder[1]:.4g}")
        print(f"[(c1)] debug SAR vs WF15_DEG2_RECORD {c1:.4%} ({c1 / g.L0_BAND:.2e} x band)")
        print(f"[(c2)] dev S vs _record(2) {c2:.4f} (band {g.LEG_D0_REPRODUCTION_BAND:.0e})")
        print(f"[(cL)] benchmark-lite price / debug elapsed = {cl:.3f}, outside "
              f"[{1 / g.BRACKET:.3g}, {g.BRACKET:.3g}]")
        print(f"[(cW)] -n 8 price applied at -n 2 unscaled misses by {cw_rel:.4f} "
              f"(record check 1e-2)")
        for label in ("debug", "dev"):
            print(f"[{label} run dir] " + ", ".join(
                f"{n} {sz} B" for n, sz in runs[label]["listing"]))

    assert est_b.ranks == 2 and est_b.degree == 2
    assert est_d.cells == sum_d["cells"] and est_v.cells == sum_v["cells"], "(E) cells"
    assert est_b.cells == g.ANS6_EXACT_STRIP_HOLE_CELL_RECORD, "(E) bl cells"
    assert est_b.unknowns == REC["bl_unknowns"], "(E) bl unknowns"
    assert sum_d["degree"] == 1 and sum_d["ranks"] == 2
    assert r_debug <= g.L0_BAND, f"(R-debug) {r_debug:.3e}"
    assert sum_v["degree"] == 1 and sum_v["family"] == "four_port_birdcage"
    assert r_dev <= g.LEG_D0_REPRODUCTION_BAND, f"(R-dev) {r_dev:.3e}"
    assert recip <= g.RECIPROCITY_BAND
    assert float(np.max(sigma)) <= 1.0 + g.PASSIVITY_SIGMA_TOLERANCE
    assert all(v <= 1.0 + g.PASSIVITY_SIGMA_TOLERANCE for v in power)
    for name, v in spreads.items():
        assert v <= g.ADJACENT_SPREAD_BAND, name
    assert 1.0 / g.BRACKET <= p_debug <= g.BRACKET, f"(P-debug) {p_debug:.3f}"
    assert 1.0 / g.BRACKET <= p_dev <= g.BRACKET, f"(P-dev) {p_dev:.3f}"
    assert w_rel <= 1e-9, f"(W) {w_rel:.3e}"
    assert ladder[0] < ladder[1] < ladder[2], f"(L) {ladder}"
    assert c1 > g.L0_BAND, "(c1) blind"
    assert c2 > g.LEG_D0_REPRODUCTION_BAND, "(c2) blind"
    assert not (1.0 / g.BRACKET <= cl <= g.BRACKET), f"(cL) blind {cl:.3f}"
    assert cw_rel > 1e-2, f"(cW) blind {cw_rel:.3e}"

    # ---- fields of the debug run ----------------------------------------------------
    res = run_d.result
    msh = res.mesh
    e = res.fields.e_complex
    p1_space = fem.functionspace(msh, ("Lagrange", 1))
    e_abs = fem.Function(p1_space, name="E_abs")
    e_abs.interpolate(
        fem.Expression(ufl.sqrt(ufl.real(ufl.inner(e, e))), p1_space.element.interpolation_points)
    )
    e_abs.x.array[:] = np.real(e_abs.x.array)
    e_abs.x.scatter_forward()
    dg0 = fem.functionspace(msh, ("DG", 0))
    sar_field = fem.Function(dg0, name="SAR")
    sar_field.interpolate(
        fem.Expression(
            res.fields.sigma_field * ufl.real(ufl.inner(e, e)) / (2.0 * run_d.study.density_kg_m3),
            dg0.element.interpolation_points,
        )
    )
    sar_field.x.array[:] = np.real(sar_field.x.array)
    sar_field.x.scatter_forward()
    fields = {"E_abs": e_abs, "SAR": sar_field}
    for key, fn in fields.items():
        assert fn.name == key, (key, fn.name)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, res.cell_tags, fields, comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    if comm.rank == 0:
        print(f"[paraview] wrote {xdmf_path}")
        print(f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.",
              flush=True)


if __name__ == "__main__":
    main()
