"""Example (`EX-132`, `mri:29`): ``femem check`` / ``femem run`` on the three tracked F2182 gel-box
``thermal_box`` files -- steady m = 4, the h-doubled wrong-wall file, and the duty-cycled m = 6 file --
with the temperature result read back from ``summary.json`` alone.

`WF-1` step 15 (`U59`, item 521b) let a YAML file carry a thermal-only ``thermal_box`` study: no EM
solve, no ``frequency_hz``, no ``material:``.  No example ran such a file; ``mat:24`` / ``mat:25``
build Python objects, ``mri:20`` / ``21`` / ``26`` / ``27`` run EM birdcage files.  This one runs
``tests/package/wf1/thermal_box_steady_m4.yaml``, ``thermal_box_wrong_wall_m4.yaml`` and
``thermal_box_duty_m6.yaml`` in place through ``cli.main(["run", ...])`` in-process (the ``mri:21``
call-through wrapper on ``fem_em_solver.study_yaml.run_study_file``, restored in ``finally``), and
shows:

* the gate's anchors, re-asserted: (K) ``femem check`` rc and its two stdout lines, (UT) cells and
  geometry against the Python Study (and the 41 duty multipliers), (I-Py) theta_file = theta_Python at
  ``L0_BAND``, (R) the MAT-21 records at ``.4e`` with ``INTERIOR_BAND`` / ``E_BAND`` /
  ``S2_GAUSS_BAND``, (I=) the summary = its own thermal result and ``summary.json`` = the summary;
* new, read from the ``summary.json`` loaded from disk: (P-json) ``source_power_w`` = Q * V,
  (S0-json) the duty multipliers, (F-diff) the files differ only where the physics says, and
  (G-json) the steady balance ``convected_power_w`` = ``source_power_w`` (printed, not asserted).

Negative controls: (c-wall) the h-doubled file misses (asserted); (c3) the duty file against the
continuous-source step response misses (asserted); (cX) the malformed duty file refuses at
``thermal.duty``, rc (2, 2), nothing written.

Every record, band and helper is imported from the gate module (the ANS-1 rule); record digits sit
only in the example-local ``REC`` dict with their ``log:line``.  Temperatures print at <= 6
significant digits; the summary dict and the JSON are never printed whole.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:29 -n 2 -t 590

Not claimed: the installed console script or ``--ranks``; ``--estimate``; the transverse source from
YAML; P2; probes or point values; a container wall; the real build; perfusion, EM or a coil; bundle
reproducibility (``mri:27``); F2182 compliance; an external check.
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

import fem_em_solver.study_yaml as study_yaml  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.package.test_wf1_yaml_thermal_box as g  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_29_femem_run_f2182_gel_box_files"
RUN_STEADY = OUTPUT_DIR / "mri_29_run_steady"
RUN_WRONG = OUTPUT_DIR / "mri_29_run_wrong_wall"
RUN_DUTY = OUTPUT_DIR / "mri_29_run_duty"
BAD_DIR = OUTPUT_DIR / "mri_29_bad"

W1A = "20261018T185143Z_WF-1-s15-W1a-steady.log"
W1B = "20261018T185201Z_WF-1-s15-W1b-duty.log"

# Records, printed beside readings only (the gate asserts the .4e strings).
REC = {
    "steady_file_s": 0.8,                                  # W1a:63
    "i_py_steady": 1.479e-15,                              # W1a:64
    "whole": 3.4599e-03, "interior": 3.223613e-03,         # W1a:65
    "gauss": 2.708e-13,                                    # W1a:65
    "c_wall": 3.8048e-01,                                  # W1a:66
    "theta_max_steady": 12.0074, "power_steady": 49.1400,  # W1a:67
    "duty_file_s": 9.4,                                    # W1b:63
    "i_py_duty": 4.792e-14,                                # W1b:64
    "t_duty": 9.4014e-04, "balance": 2.809e-12,            # W1b:65
    "c3": 5.0282e-01,                                      # W1b:66
    "theta_max_duty": 0.667259,                            # W1b:67
}


def fmt(v: float) -> str:
    """theta at <= 6 significant digits (the leak guard matches at 8)."""
    v = float(v)
    return f"{v:.4f}" if abs(v) >= 10.0 else f"{v:.5g}"


def hexeq(a, b) -> bool:
    return float(a).hex() == float(b).hex()


def flat_diff(a, b, prefix=""):
    """Dotted keys at which two YAML mappings differ (a key on one side only is a diff at that key)."""
    out = set()
    for k in set(a) | set(b):
        path = f"{prefix}{k}"
        if k not in a or k not in b:
            out.add(path)
        elif isinstance(a[k], dict) and isinstance(b[k], dict):
            out |= flat_diff(a[k], b[k], path + ".")
        elif a[k] != b[k]:
            out.add(path)
    return out


def cli_run(path, out_dir):
    """``cli.main(["run", ...])`` with the StudyRun captured by a call-through wrapper."""
    captured = []
    original = study_yaml.run_study_file

    def _call_through(*args, **kwargs):
        res = original(*args, **kwargs)
        captured.append(res)
        return res

    study_yaml.run_study_file = _call_through
    try:
        rc = g.cli.main(["run", str(path), "--output-dir", str(out_dir)])
    finally:
        study_yaml.run_study_file = original
    assert study_yaml.run_study_file is original
    return rc, captured


def cli_check(path):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = g.cli.main(["check", str(path)])
    return rc, buf.getvalue().splitlines()


def read_summary(run, comm):
    comm.Barrier()
    with open(run.paths["summary"]) as fh:
        return json.load(fh)


def ones_tags(msh):
    tdim = msh.topology.dim
    imap = msh.topology.index_map(tdim)
    n = imap.size_local + imap.num_ghosts
    return dmesh.meshtags(msh, tdim, np.arange(n, dtype=np.int32), np.full(n, 1, dtype=np.int32))


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    def say(msg=""):
        if rank0:
            print(msg, flush=True)

    # ---- (K) structure and femem check --------------------------------------------------------
    assert g.IS_COMPLEX, "complex DolfinX build required (the runner's `mri:` group sources it)"
    assert comm.size == 2, f"comm.size {comm.size}"
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    say("=" * 72)
    say("EX-132 -- femem check / run on the three F2182 gel-box thermal_box files (mri:29)")
    say("=" * 72)
    say(f"[records] beside {W1A} and {W1B}; ranks {comm.size}")
    say(f"[bands] imported: L0_BAND {g.L0_BAND:.0e}, INTERIOR_BAND {g.INTERIOR_BAND:.0e}, "
        f"E_BAND {g.E_BAND:.0e}, S2_GAUSS_BAND {g.S2_GAUSS_BAND:.0e}")

    check_rows = {}
    for f in (g.STEADY, g.WRONG_WALL, g.DUTY):
        rc, lines = cli_check(f)
        rc_all = g._rc_all(rc, comm)
        ok_lines = (comm.rank != 0) or (
            [ln for ln in lines if ln.startswith("family: ")] == ["family: thermal_box"]
            and "excitation: UniformHeatSource" in lines
        )
        ok_lines = g._all(comm, ok_lines)
        check_rows[f.name] = (rc_all, ok_lines)
        say(f"[(K)] femem check {f.name}: rc {rc_all}; one 'family: thermal_box' and "
            f"'excitation: UniformHeatSource' on rank 0: {ok_lines}")
    for name, (rc_all, ok_lines) in check_rows.items():
        assert rc_all == (0, 0), (name, rc_all)
        assert ok_lines, name

    # (F-diff) the files differ exactly where the physics says
    m_st, m_du, m_ww = (g._mapping(f) for f in (g.STEADY, g.DUTY, g.WRONG_WALL))
    d_sd = flat_diff(m_st, m_du)
    d_sw = flat_diff(m_st, m_ww)
    say(f"[(F-diff)] steady vs duty: {sorted(d_sd)}; steady vs wrong wall: {sorted(d_sw)}")
    assert d_sd == {"geometry.divisions", "thermal.transient", "thermal.duty"}, d_sd
    assert d_sw == {"thermal.h_w_per_m2_k"}, d_sw

    # (cX) the malformed file: refusal, rc (2, 2), nothing written
    cx_field = g._field(lambda: g.load_study(g.MALFORMED_DUTY))
    g._fresh_dir(BAD_DIR, comm)
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        rc_bad_run = g.cli.main(["run", str(g.MALFORMED_DUTY), "--output-dir", str(BAD_DIR)])
        rc_bad_check = g.cli.main(["check", str(g.MALFORMED_DUTY)])
    cx_run, cx_check = g._rc_all(rc_bad_run, comm), g._rc_all(rc_bad_check, comm)
    cx_exists = comm.allreduce(bool(BAD_DIR.exists()), op=MPI.LOR)
    named = comm.allreduce(int("thermal.duty" in err.getvalue()), op=MPI.MAX)
    say(f"[(cX)] {g.MALFORMED_DUTY.name}: field {cx_field!r}; run rc {cx_run}, check rc {cx_check}; "
        f"out dir exists {cx_exists}; 'thermal.duty' on stderr of some rank {bool(named)} "
        "(expected 'thermal.duty', (2, 2), (2, 2), False)")
    assert cx_field == "thermal.duty"
    assert cx_run == (2, 2) and cx_check == (2, 2)
    assert not cx_exists

    q_v = g.Q * 8.0 * float(np.prod(g.HALF))  # Q * V = 49.14 W

    # ==== steady leg: the gate's W1a body with cli.main and RUN_DIR substituted ==============
    t00 = time.perf_counter()
    g._fresh_dir(RUN_STEADY, comm)
    rc, cap = cli_run(g.STEADY, RUN_STEADY)
    t_file = time.perf_counter() - t00
    assert g._rc_all(rc, comm) == (0, 0) and len(cap) == 1, (rc, len(cap))
    run = cap[0]
    th = run.thermal
    py = g._py_study(g.DIV4).temperature(g.RobinCooling(g.K_GEL, g.H_GEL, g.T_INF))
    msh = th.submesh

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    cell_tags = ones_tags(msh)
    write_setup_figure(
        msh, cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel (650 x 420 x 90 mm)"},
        title="mri:29 -- F2182-sized gel box from femem run (thermal_box file), uniform q, Robin on "
              "all six faces: z = 0 slice",
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        clip_normal=(0.0, 0.0, -1.0), comm=comm,
    )

    ca, cb, geo = g._same_mesh(msh, py.submesh, comm)
    T_py = g._on(msh, py.temperature)
    ip = g._theta_owned_rel(T_py, th.temperature, comm)
    n_own = th.temperature.function_space.dofmap.index_map.size_local
    hx = g._all(comm, np.array_equal(np.asarray(T_py.x.array[:n_own]),
                                     np.asarray(th.temperature.x.array[:n_own])))
    th_h = th.temperature - g.T_INF
    ex, npts = g._exact(msh, g._steady_ex())
    whole = g._rel(msh, th_h, ex)
    chi = g._interior_chi(msh)
    interior = g._rel(msh, th_h, ex, chi)
    ds = ufl.Measure("ds", domain=msh, subdomain_data=th.facet_tags)(g.ROBIN_TAG)
    conv = float(np.real(g._allsum(comm, g.H_GEL * th_h * ds)))
    one = fem.Constant(msh, default_scalar_type(1.0))
    vh = float(np.real(g._allsum(comm, one * ufl.dx(domain=msh))))
    gauss = abs(conv - g.Q * vh) / (g.Q * vh)
    # (c-wall) through the tracked wrong-wall file
    g._fresh_dir(RUN_WRONG, comm)
    rc_w, cap_w = cli_run(g.WRONG_WALL, RUN_WRONG)
    assert g._rc_all(rc_w, comm) == (0, 0) and len(cap_w) == 1
    run_w = cap_w[0]
    ca_w, _cb_w, geo_w = g._same_mesh(run_w.thermal.submesh, msh, comm)
    T_w = g._on(msh, run_w.thermal.temperature)
    c_wall = g._rel(msh, T_w - g.T_INF, ex, chi)

    # (I=) collective; its own prints are held back (they carry the unrounded theta_max)
    with contextlib.redirect_stdout(io.StringIO()):
        g._check_summary(run, comm, "steady", transient=False)
    disk_s = read_summary(run, comm)
    ts = disk_s["temperature"]
    p_json_s = abs(ts["source_power_w"] / q_v - 1.0)
    g_json = abs(ts["convected_power_w"] - ts["source_power_w"]) / ts["source_power_w"]

    say()
    say(f"[(UT)] steady: summary cells {run.summary['cells']}, cells {ca} / {cb} (record "
        f"{g.STEP1_M4_CELLS}; {W1A}:63 79 872); geometry equal {geo}; file run {t_file:.1f} s "
        f"(record {REC['steady_file_s']} s)")
    say(f"[(I-Py)] steady: owned max|dtheta|/max|theta| = {ip:.3e} (record "
        f"{REC['i_py_steady']:.3e}, {W1A}:64; band {g.L0_BAND:.0e}); owned arrays bitwise equal {hx} "
        "(printed, not asserted)")
    say(f"[(R)] steady: whole L2 rel {whole:.6e} -> {whole:.4e} (record {g.STEP1_M4_WHOLE}); "
        f"interior {interior:.6e} (record {REC['interior']:.6e}; band {g.INTERIOR_BAND:.0e}); "
        f"Gauss {gauss:.3e} (record {REC['gauss']:.3e}; band {g.S2_GAUSS_BAND:.0e}); {npts} P2 points "
        f"({W1A}:65)")
    say(f"[(c-wall)] wrong-wall file: interior miss {c_wall:.4e} (record {REC['c_wall']:.4e}, "
        f"{W1A}:66; 521a API {g.WF15_S5_C_WALL:.4e}; asserted > {g.INTERIOR_BAND:.0e}, < 1); "
        f"mesh equal {geo_w}, cells {ca_w}")
    say(f"[(I=)] steady: summary == own thermal result by hex and summary.json == summary "
        f"(asserted inside the gate helper); summary keys {sorted(run.summary)}; temperature keys "
        f"{sorted(ts)}")
    say(f"[(P-json)] steady: source_power_w {ts['source_power_w']:.4f} W vs Q*V {q_v:.4f} W: rel "
        f"{p_json_s:.3e} (band {g.L0_BAND:.0e}; record {REC['power_steady']:.4f} W, {W1A}:67; "
        "predicted ~2e-13)")
    say(f"[(G-json)] steady: |convected_power_w - source_power_w| / source_power_w = {g_json:.3e} "
        f"(printed only; predicted ~3e-13, range [0, 1e-11]; within S2_GAUSS_BAND "
        f"{g_json <= g.S2_GAUSS_BAND})")
    say(f"[theta] steady theta_max from the JSON {fmt(ts['theta_max_k'])} K (record "
        f"{fmt(REC['theta_max_steady'])} K, {W1A}:67)")

    assert run.summary["cells"] == ca == cb == g.STEP1_M4_CELLS
    assert geo and geo_w
    assert ip <= g.L0_BAND, f"(I-Py) steady {ip:.3e}"
    assert f"{whole:.4e}" == g.STEP1_M4_WHOLE, f"(R) whole {whole:.4e}"
    assert interior <= g.INTERIOR_BAND
    assert gauss <= g.S2_GAUSS_BAND
    assert run_w.summary["temperature"]["transient"] is None
    assert c_wall > g.INTERIOR_BAND, f"(c-wall) did not miss: {c_wall:.4e}"
    assert c_wall < 1.0, f"(c-wall) above ceiling: {c_wall:.4e}"
    assert p_json_s <= g.L0_BAND, f"(P-json) steady {p_json_s:.3e}"
    t_steady_leg = time.perf_counter() - t00

    # ---- steady fields to ParaView (P1, real parts) ------------------------------------------
    V1 = fem.functionspace(msh, ("Lagrange", 1))

    def p1(name, src, offset=0.0):
        f = fem.Function(V1, name=name)
        f.interpolate(src)
        f.x.array[:] = np.real(f.x.array) - offset
        f.x.scatter_forward()
        return f

    f_file = p1("theta_file_steady", th.temperature, g.T_INF)
    f_ser = p1("theta_series_steady", ex)
    f_dif = fem.Function(V1, name="theta_file_minus_series")
    f_dif.x.array[:] = np.real(f_file.x.array) - np.real(f_ser.x.array)
    f_dif.x.scatter_forward()
    f_ww = p1("theta_wrong_wall", T_w, g.T_INF)
    ft = th.facet_tags
    ft.name = "robin_facets"
    steady_fields = {f.name: f for f in (f_file, f_ser, f_dif, f_ww)}
    xdmf_s, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_steady_combined", msh, cell_tags, steady_fields, comm=comm,
        facet_tags=ft)

    # ==== duty leg: the gate's W1b body with cli.main and RUN_DIR substituted ================
    t01 = time.perf_counter()
    g._fresh_dir(RUN_DUTY, comm)
    rc_d, cap_d = cli_run(g.DUTY, RUN_DUTY)
    t_file_d = time.perf_counter() - t01
    assert g._rc_all(rc_d, comm) == (0, 0) and len(cap_d) == 1, (rc_d, len(cap_d))
    run_d = cap_d[0]
    thd = run_d.thermal
    pyd = g._py_study(g.DIV6).temperature(
        g.RobinCooling(g.K_GEL, g.H_GEL, g.T_INF),
        transient=g.TransientSpec(g.DT_BOX, g.N_MAIN, g.RHO_C),
        duty=g.DutyCycle(period_s=2.0 * g.HALF_T, duty=0.5),
    )
    mshd = thd.mesh
    cad, cbd, geod = g._same_mesh(mshd, pyd.mesh, comm)
    mult = np.array_equal(np.asarray(thd.source_multipliers), np.asarray(pyd.source_multipliers))
    T_pyd = g._on(mshd, pyd.temperature)
    ipd = g._theta_owned_rel(T_pyd, thd.temperature, comm)
    th_hd = thd.temperature - g.T_INF
    chid = g._interior_chi(mshd)
    exd, nptsd = g._exact(mshd, lambda x, y, z: g._duty_cf(x, y, z, g.T_STAR))
    t_duty = g._rel(mshd, th_hd, exd, chid)
    ex_s, _ = g._exact(mshd, g._trans_ex(g.T_STAR))
    c3 = g._rel(mshd, th_hd, ex_s, chid)
    bal = run_d.summary["temperature"]["max_balance_residual"]

    with contextlib.redirect_stdout(io.StringIO()):
        g._check_summary(run_d, comm, "duty", transient=True)
    disk_d = read_summary(run_d, comm)
    td = disk_d["temperature"]
    p_json_d = abs(td["source_power_w"] / q_v - 1.0)
    mults_json = td["source_multipliers"]
    mults_mem = [float(m) for m in thd.source_multipliers]
    mean_mult = float(np.mean(mults_json))

    say()
    say(f"[(UT)] duty: summary cells {run_d.summary['cells']}, cells {cad} / {cbd} (record "
        f"{g.STEP1_M6_CELLS}; {W1B}:63 269 568); geometry equal {geod}; multipliers equal {mult} "
        f"({len(thd.source_multipliers)} entries; record 41); file run {t_file_d:.1f} s (record "
        f"{REC['duty_file_s']} s)")
    say(f"[(I-Py)] duty: owned max|dtheta|/max|theta| = {ipd:.3e} (record {REC['i_py_duty']:.3e}, "
        f"{W1B}:64; band {g.L0_BAND:.0e})")
    say(f"[(R)] duty: (T-duty) interior rel {t_duty:.6e} -> {t_duty:.4e} (record {g.P26_T_DUTY}; "
        f"{REC['t_duty']:.4e}; band {g.INTERIOR_BAND:.0e}); balance {bal:.3e} (record "
        f"{REC['balance']:.3e}; band {g.E_BAND:.0e}); {nptsd} P2 points ({W1B}:65)")
    say(f"[(c3)] step response vs duty file: miss {c3:.4e} (record {REC['c3']:.4e}, {W1B}:66; "
        f"P26 {g.P26_C3:.4e}; asserted > {g.INTERIOR_BAND:.0e})")
    say(f"[(I=)] duty: summary == own thermal result by hex and summary.json == summary "
        f"(asserted inside the gate helper); temperature keys {sorted(td)}")
    say(f"[(P-json)] duty: source_power_w {td['source_power_w']:.4f} W vs Q*V {q_v:.4f} W: rel "
        f"{p_json_d:.3e} (band {g.L0_BAND:.0e}; predicted ~2.6e-12, {W1B}:67)")
    say(f"[(S0-json)] duty: source_multipliers in the JSON {len(mults_json)} entries == in-memory "
        f"{mults_json == mults_mem}; mean multiplier {mean_mult:.4f} (predicted ~0.5; printed only)")
    say(f"[theta] duty theta_max from the JSON {fmt(td['theta_max_k'])} K (record "
        f"{fmt(REC['theta_max_duty'])} K, {W1B}:67); steady-to-duty theta_max ratio "
        f"{ts['theta_max_k'] / td['theta_max_k']:.4g} (no band)")

    assert run_d.summary["cells"] == cad == cbd == g.STEP1_M6_CELLS
    assert geod
    assert mult
    assert run_d.summary["temperature"]["duty"] == {"period_s": 2.0 * g.HALF_T, "duty": 0.5}
    assert ipd <= g.L0_BAND, f"(I-Py) duty {ipd:.3e}"
    assert t_duty <= g.INTERIOR_BAND
    assert f"{t_duty:.4e}" == g.P26_T_DUTY, f"(T-duty) {t_duty:.4e}"
    assert bal <= g.E_BAND
    assert c3 > g.INTERIOR_BAND, f"(c3) did not miss: {c3:.4e}"
    assert len(mults_json) == 41 and mults_json == mults_mem
    assert p_json_d <= g.L0_BAND, f"(P-json) duty {p_json_d:.3e}"

    # ---- duty fields to ParaView ---------------------------------------------------------------
    V1d = fem.functionspace(mshd, ("Lagrange", 1))

    def p1d(name, src, offset=0.0):
        f = fem.Function(V1d, name=name)
        f.interpolate(src)
        f.x.array[:] = np.real(f.x.array) - offset
        f.x.scatter_forward()
        return f

    d_file = p1d("theta_file_duty_tstar", thd.temperature, g.T_INF)
    d_ser = p1d("theta_duty_series_tstar", exd)
    d_dif = fem.Function(V1d, name="theta_duty_minus_series")
    d_dif.x.array[:] = np.real(d_file.x.array) - np.real(d_ser.x.array)
    d_dif.x.scatter_forward()
    d_step = p1d("theta_step_series_tstar", ex_s)
    duty_fields = {f.name: f for f in (d_file, d_ser, d_dif, d_step)}
    xdmf_d, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_duty_combined", mshd, ones_tags(mshd), duty_fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- CSV: the summary read-back, rank 0 ------------------------------------------------------
    csv_path = OUTPUT_DIR / f"{BASENAME}_summary_readback.csv"
    rows = []
    for label, disk, mem in (("steady", disk_s, run.summary), ("duty", disk_d, run_d.summary)):
        pairs = [("cells", disk["cells"], mem["cells"])]
        for key, v in disk["temperature"].items():
            if isinstance(v, float):
                pairs.append((key, v, mem["temperature"][key]))
        for key, v_json, v_mem in pairs:
            is_t = key == "theta_max_k"
            sj = fmt(v_json) if is_t else f"{v_json:.6g}"
            sm = fmt(v_mem) if is_t else f"{v_mem:.6g}"
            rows.append((f"{label}", key, sj, sm, hexeq(v_json, v_mem)))
    if rank0:
        with open(csv_path, "w") as fh:
            fh.write("file,key,json_value,in_memory_value,hex_equal\n")
            for r in rows:
                fh.write(",".join(str(x) for x in r) + "\n")
    comm.Barrier()
    assert all(r[4] for r in rows), [r for r in rows if not r[4]]

    say()
    say(f"[output] {xdmf_s.name if rank0 else ''}: P1 theta_file_steady, theta_series_steady, "
        "theta_file_minus_series, theta_wrong_wall, Robin facets")
    say(f"[output] {xdmf_d.name if rank0 else ''}: P1 theta_file_duty_tstar, "
        "theta_duty_series_tstar, theta_duty_minus_series, theta_step_series_tstar")
    say(f"[output] {csv_path.name} ({len(rows)} rows, every JSON value hex-equal to the in-memory "
        "summary)")
    say(f"[time] steady leg {t_steady_leg:.1f} s; duty leg {time.perf_counter() - t01:.1f} s; "
        f"total {time.perf_counter() - started:.1f} s")
    say("[scope] three thermal_box files through cli.main, steady m = 4 and duty-cycled CN m = 6, "
        "uniform source, Robin on six faces, P1, complex build, -n 2; summary read back from disk; "
        "no EM, no coil, no probes, no P2, no compliance, no external check claimed")
    say("All assertions hold")


if __name__ == "__main__":
    main()
