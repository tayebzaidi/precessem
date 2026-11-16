"""Example (`EX-136`, `mri:32`): ``femem run`` on the tracked **duty-cycled transverse-quadratic**
``thermal_box`` file (`WF-1` step 18, §9 item 540) and the tracked **ball-probe** ``thermal_ball`` files
(`WF-1` step 19, §9 item 545), with the temperature, the schedule and the ball probe table read back
from the run's JSON summary on disk.

No example before this one ran a ``thermal_ball`` file or a duty-cycled non-uniform source from a file.
``mri:31`` runs the quadratic source without ``thermal.duty`` and probes on the box; ``mri:25`` builds
Python objects; ``mat:13`` reads MAT-19's probes from no file.  This one runs
``tests/package/wf1/thermal_box_quadratic_duty_m4.yaml``, ``thermal_ball_m19_probes.yaml`` and
``thermal_ball_duty_a_probes.yaml`` in place through ``cli.main(["run", ...])`` in-process (the
``mri:21`` / ``mri:29`` / ``mri:31`` call-through wrapper on ``fem_em_solver.study_yaml.run_study_file``,
restored in ``finally``), plus ``femem check`` on those three and on the malformed ball file.  **No
tracked file carries both** ``excitation.source`` and ``outputs.probes``, and none combines a wall with
probes; those combinations are not gated and are not run here.

The three solve-node bodies of the two gate modules
(``test_wf1_yaml_thermal_box_quadratic_duty``: quadratic duty; ``test_wf1_yaml_thermal_ball_probes``:
m19 and duty) are copied with six substitutions (``run_study_file`` -> the captured ``cli_run``;
``LOGS`` -> this example's run directories; the real-build digit equality of (T*-ball) -> a printed
equality; full-``repr`` probe-value prints -> ``.6g``; the box ``_check_summary`` run with its stdout
held back; every bare name qualified by its module alias).  Every record, band and helper is imported
from the gate modules (the ANS-1 rule); record digits sit only in the example-local ``REC`` dict with
their ``log:line``.

The quadratic-duty leg is a complex-build gate and runs in the build of record.  The ball-probe gate ran
in the real build only; here it runs in the complex build, so its **bands** are asserted and its
**digits** are printed beside U60's (the item 546 / `EX-135` pattern).  These are example records, not a
complex gate of `WF-1` step 19.

Run it through the example runner::

    ./scripts/run_examples.sh -e mri:32 -n 2 -t 590

Not claimed: a source with probes in one file (ungated); a walled ball with probes (refused); any
comparand on the duty-ball probes; probes at P2; scalar point values; a complex gate of the ball-probe
file route; the real build; the installed console script, ``--ranks`` or ``--estimate``; EM, a coil or
SAR; perfusion; bundle reproducibility (``mri:27``); F2182 compliance; an external check.
"""

from __future__ import annotations

import contextlib
import io
import json
import sys
import time
from pathlib import Path

import numpy as np
from dolfinx import fem
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

import tests.package.test_wf1_yaml_thermal_ball_probes as b  # noqa: E402
import tests.package.test_wf1_yaml_thermal_box_quadratic as q  # noqa: E402
import tests.package.test_wf1_yaml_thermal_box_quadratic_duty as d  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_32_femem_run_ball_probes_and_duty_quadratic_files"
RUN_QD = OUTPUT_DIR / "mri_32_run_quadratic_duty"
RUN_M19 = OUTPUT_DIR / "mri_32_run_ball_m19"
RUN_BD = OUTPUT_DIR / "mri_32_run_ball_duty"
BAD_DIR = OUTPUT_DIR / "mri_32_bad"

D18 = "20261019T094516Z_WF-1-s18-W1-qduty.log"
D0 = "20261019T094458Z_WF-1-s18-W0.log"
B19 = "20261019T110412Z_WF-1-s19-W1-ball-probes.log"
B0 = "20261019T110312Z_WF-1-s19-W0.log"
E123 = "20261017T005319Z_EX-123.log"
E135 = "20261019T123750Z_EX-135.log"

# Records, printed beside readings only (the copied quadratic-duty body asserts the .6e P38 string).
REC = {
    "qd_file_s": 33.9, "qd_py_s": 32.8, "qd_leg_s": 117.5,   # D18:63, :71
    "i_py_qd": 2.221e-14, "e_qd": 2.288e-13, "c_step": 5.0333e-01,  # D18:66, :68, :69
    "src_power_rel": 1.157e-13, "theta_max_qd": 1.43990,       # D18:70
    "theta_max_cn_unsched": 2.9282,                              # E135:86
    "m19_file_s": 4.1, "m19_leg_s": 6.7,                         # B19:80, :90
    "i_py_m19": 2.973e-13,                                       # B19:81
    "a_ball": (4.511e-16, 1.332e-15),                            # B19:84
    "c_disp": 5.4208e-02, "sep": 5.4663e-02,                     # B19:88
    "i_du": 8.146e-15, "a_du": (8.476e-15, 6.994e-15),           # B19:93, :94
    "duty_leg_s": 51.2, "bal_du": 2.632e-14,                     # B19:98, E123:45
}


def fmt(v: float) -> str:
    """theta at <= 6 significant digits (the leak guard matches at 8)."""
    v = float(v)
    return f"{v:.4f}" if abs(v) >= 10.0 else f"{v:.5g}"


def hexeq(a, b_) -> bool:
    return float(a).hex() == float(b_).hex()


def flat_diff(a, b_, prefix=""):
    """Dotted keys at which two YAML mappings differ (a key on one side only is a diff at that key)."""
    out = set()
    for k in set(a) | set(b_):
        path = f"{prefix}{k}"
        if k not in a or k not in b_:
            out.add(path)
        elif isinstance(a[k], dict) and isinstance(b_[k], dict):
            out |= flat_diff(a[k], b_[k], path + ".")
        elif a[k] != b_[k]:
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
        rc = d.cli.main(["run", str(path), "--output-dir", str(out_dir)])
    finally:
        study_yaml.run_study_file = original
    assert study_yaml.run_study_file is original
    return rc, captured


def cli_check(path):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = d.cli.main(["check", str(path)])
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


def p1_field(msh, name, src, offset=0.0):
    """A P1 Function on ``msh`` interpolated from ``src`` (a Function or callable), real part."""
    V1 = fem.functionspace(msh, ("Lagrange", 1))
    f = fem.Function(V1, name=name)
    f.interpolate(src)
    f.x.array[:] = np.real(f.x.array) - offset
    f.x.scatter_forward()
    return f


def diff_field(msh, name, a, b_):
    V1 = fem.functionspace(msh, ("Lagrange", 1))
    f = fem.Function(V1, name=name)
    f.x.array[:] = np.real(a.x.array) - np.real(b_.x.array)
    f.x.scatter_forward()
    return f


def json_hex(js):
    return [([float(x).hex() for x in c], float(r).hex())
            for c, r in zip(js["centres"], js["radii"])]


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    def say(msg=""):
        if rank0:
            print(msg, flush=True)

    _p0 = say

    L0_BAND, INTERIOR_BAND, E_BAND, B_BAND, P_BAND = (d.L0_BAND, d.INTERIOR_BAND, d.E_BAND,
                                                      b.B_BAND, b.P_BAND)
    assert d.L0_BAND == b.L0_BAND

    # ---- (K) structure and femem check --------------------------------------------------------
    assert d.IS_COMPLEX, "complex DolfinX build required (the runner's `mri:` group sources it)"
    assert comm.size == 2, f"comm.size {comm.size}"
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    say("=" * 72)
    say("EX-136 -- femem check / run on the duty-cycled quadratic thermal_box file and the two "
        "ball-probe thermal_ball files (mri:32)")
    say("=" * 72)
    say(f"[records] beside {D18} (complex), {B19} (real build), {E123}, {E135}; ranks {comm.size}")
    say(f"[bands] imported: L0_BAND {L0_BAND:.0e}, INTERIOR_BAND {INTERIOR_BAND:.0e}, "
        f"E_BAND {E_BAND:.0e}, B_BAND {B_BAND:.0e}, P_BAND {P_BAND:.0e}")

    check_rows = {}
    for f, fam, excit in ((d.QUAD_DUTY, "thermal_box", "TransverseQuadraticHeatSource"),
                          (b.M19, "thermal_ball", None), (b.DUTY, "thermal_ball", None)):
        rc, lines = cli_check(f)
        rc_all = d._rc_all(rc, comm)
        ok = (comm.rank != 0) or (
            [ln for ln in lines if ln.startswith("family: ")] == [f"family: {fam}"]
            and (excit is None or f"excitation: {excit}" in lines)
        )
        ok = b._all(comm, ok)
        exline = [ln for ln in lines if ln.startswith("excitation: ")] if rank0 else []
        check_rows[f.name] = (rc_all, ok)
        say(f"[(K)] femem check {f.name}: rc {rc_all}; one 'family: {fam}'"
            f"{'' if excit is None else ' and ' + repr('excitation: ' + excit)} on rank 0: {ok}; "
            f"rank-0 excitation line {exline}")
    for name, (rc_all, ok) in check_rows.items():
        assert rc_all == (0, 0), (name, rc_all)
        assert ok, name

    # (F-diff) the files differ exactly where the physics says
    m_qd, m_qc = d._mapping(d.QUAD_DUTY), d._mapping(q.QUAD_CN)
    m_m19, m_du, m_dp, m_bad = (d._mapping(f) for f in (b.M19, b.DUTY, b.DUTY_PLAIN, b.MALFORMED))
    diffs = {
        "quad-duty vs quad-cn": (flat_diff(m_qd, m_qc), {"thermal.duty"}),
        "ball-duty-probes vs ball-duty": (flat_diff(m_du, m_dp), {"outputs"}),
        "ball-m19 vs malformed": (flat_diff(m_m19, m_bad), {"outputs.probes"}),
    }
    for label, (got, want) in diffs.items():
        say(f"[(F-diff)] {label}: {sorted(got)}")
    for label, (got, want) in diffs.items():
        assert got == want, (label, sorted(got), sorted(want))
    # (F-none) no file the example runs carries both capabilities, and no ball file a wall
    all_maps = {f.name: d._mapping(f) for f in (d.QUAD_DUTY, q.QUAD_CN, b.M19, b.DUTY, b.DUTY_PLAIN,
                                                b.MALFORMED)}
    both = [n for n, m in all_maps.items() if "source" in m.get("excitation", {}) and "outputs" in m]
    walled = [n for n, m in all_maps.items()
              if n.startswith("thermal_ball") and "wall" in m.get("geometry", {})]
    say(f"[(F-none)] files carrying both excitation.source and outputs: {both}; ball files carrying "
        f"geometry.wall: {walled} (over {len(all_maps)} mappings; expected [] and [])")
    assert both == [] and walled == []

    # (cX) the malformed ball probe file: refusal, rc (2, 2), nothing written (copied b:205-218)
    assert b._field(lambda: b.load_study(b.MALFORMED)) == "outputs.probes[1].radius"
    cx_field = b._field(lambda: b.load_study(b.MALFORMED))
    d._fresh_dir(BAD_DIR, comm)
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        rc_run = d.cli.main(["run", str(b.MALFORMED), "--output-dir", str(BAD_DIR)])
        rc_check = d.cli.main(["check", str(b.MALFORMED)])
    run_all, check_all = d._rc_all(rc_run, comm), d._rc_all(rc_check, comm)
    exists = comm.allreduce(bool(BAD_DIR.exists()), op=MPI.LOR)
    named = comm.allreduce(int("outputs.probes" in err.getvalue()), op=MPI.MAX)
    say(f"[(cX)] {b.MALFORMED.name}: field {cx_field!r}; run rc {run_all}, check rc {check_all}; out dir "
        f"exists {exists}; 'outputs.probes' on stderr of some rank {bool(named)} (printed, not "
        f"asserted; expected 'outputs.probes[1].radius', (2, 2), (2, 2), False; {B0}:69)")
    assert cx_field == "outputs.probes[1].radius"
    assert run_all == (2, 2) and check_all == (2, 2)
    assert not exists

    # (cX-duty) duty without transient raises at thermal.duty (copied d:180-182); no mesh built
    got = d._field(lambda: d._study_from_mapping(d._without(d._mapping(d.QUAD_DUTY), "thermal",
                                                            "transient")))
    say(f"[(cX-duty)] quadratic duty file without transient: want 'thermal.duty' got {got!r} ({D0}:68)")
    assert got == "thermal.duty"

    q_v = d.Q * 8.0 * float(np.prod(d.HALF))  # Q * V

    # ==== quadratic duty leg: the gate's W1 body (d:205-254) ===================================
    t00 = time.perf_counter()
    d._fresh_dir(RUN_QD, comm)
    rc, cap = cli_run(d.QUAD_DUTY, RUN_QD)
    t_file = time.perf_counter() - t00
    assert d._rc_all(rc, comm) == (0, 0) and len(cap) == 1, (rc, len(cap))
    run = cap[0]
    th = run.thermal
    msh = th.mesh
    py = d._quad_study().temperature(d._cooling(), transient=d._transient(), duty=d._duty())
    t_py = time.perf_counter() - t00 - t_file
    ca, cb, geo = d._same_mesh(msh, py.mesh, comm)
    T_py = d._on_p(msh, py.temperature)
    ip = d._owned_rel(T_py, th.temperature, comm)
    dofs = [d._ndofs(th.temperature), d._ndofs(py.temperature)]
    st = run.summary["temperature"]
    s0 = "".join(str(int(m)) for m in st["source_multipliers"])
    mult = np.array_equal(np.asarray(th.source_multipliers), np.asarray(py.source_multipliers))

    th_t = th.temperature - d.T_INF
    chi = d._interior_chi(msh)
    ex_d, npts = d._exact_p(msh, lambda x, y, z: d._quad_duty(x, y, z, d.T_STAR))
    t_qduty = d._rel_p(msh, th_t, ex_d, chi)
    ex_s, _ = d._exact_p(msh, d._quad_trans_ex(d.T_STAR))
    c_step = d._rel_p(msh, th_t, ex_s, chi)
    bal = st["max_balance_residual"]

    say()
    _p0(f"WF1S18 (Pdeg) P2 dofs file / Python {dofs} (P38 {d.P38_DOFS}; {D18}:63); file run "
        f"{t_file:.1f}s (record {REC['qd_file_s']} s), Python run {t_py:.1f}s (record "
        f"{REC['qd_py_s']} s); summary source {st['source']!r}; duty {st['duty']!r}")
    _p0(f"WF1S18 (UT) cells file / Python {ca} / {cb} (record {d.CELLS_M4}); geometry bitwise {geo}")
    _p0(f"WF1S18 (S0) schedule {s0} (SCHEDULE {d.SCHEDULE}); multipliers equal Python {mult}")
    _p0(f"WF1S18 (I-Py) owned max|dtheta|/max|theta| = {ip:.3e} (hex {float(ip).hex()}; band "
        f"{L0_BAND:.0e}; record {REC['i_py_qd']:.3e}, {D18}:66)")
    _p0(f"WF1S18 (T-qduty) interior vs duty quadratic series {t_qduty:.6e} (P38 {d.P38_T_QDUTY}; "
        f"band {INTERIOR_BAND:.0e}); {npts} P3 points")
    _p0(f"WF1S18 (E-qduty) CN max balance residual {bal:.3e} (P38 {d.P38_E:.3e}; band {E_BAND:.0e})")
    _p0(f"WF1S18 control (c-step) duty file vs unscheduled transient series {c_step:.4e} "
        f"(P38 {d.P38_C_STEP:.4e}; asserted > {INTERIOR_BAND:.0e})")
    with contextlib.redirect_stdout(io.StringIO()):
        d._check_summary(run, comm, "quadratic duty", transient=True)
    say(f"  (I=) quadratic duty: JSON == summary asserted inside the gate helper ({D18}:70); its own "
        "line held back")

    assert dofs == [d.P38_DOFS] * 2
    assert ca == cb == d.CELLS_M4 and geo
    assert s0 == d.SCHEDULE
    assert mult
    assert st["source"] == "transverse_quadratic"
    assert d._hex_dict_equal(st["duty"], {"period_s": 2.0 * d.HALF_T, "duty": d.DUTY})
    assert ip <= L0_BAND, f"(I-Py) {ip:.3e}"
    assert t_qduty <= INTERIOR_BAND
    assert f"{t_qduty:.6e}" == d.P38_T_QDUTY
    assert bal <= E_BAND
    assert c_step > INTERIOR_BAND
    t_qd = time.perf_counter() - t00
    _p0(f"WF1S18 W1 leg {t_qd:.1f}s (record {REC['qd_leg_s']} s, {D18}:71)")

    # ---- the summary read back from disk -------------------------------------------------------
    disk_qd = read_summary(run, comm)
    tq = disk_qd["temperature"]
    s_json = "".join(str(int(m)) for m in tq["source_multipliers"])
    p_json = abs(tq["source_power_w"] / q_v - 1.0)
    ratio_unsched = tq["theta_max_k"] / REC["theta_max_cn_unsched"]
    say()
    say(f"[(S-json)] quadratic duty: source {tq['source']!r}, thermal_degree {tq['thermal_degree']}, "
        f"duty hex-equal {{2 HALF_T, DUTY}} {d._hex_dict_equal(tq['duty'], {'period_s': 2.0 * d.HALF_T, 'duty': d.DUTY})} "
        "(from disk)")
    say(f"[(Sch-json)] schedule from {len(tq['source_multipliers'])} source_multipliers on disk == "
        f"SCHEDULE {s_json == d.SCHEDULE}")
    say(f"[(P-json)] source_power_w {tq['source_power_w']:.4f} W (the full-on power, not the duty mean) "
        f"vs Q*V {q_v:.4f} W: rel {p_json:.3e} (band {L0_BAND:.0e}; record +{REC['src_power_rel']:.3e}, "
        f"{D18}:70)")
    say(f"[(E-json)] max_balance_residual {tq['max_balance_residual']:.3e} (band {E_BAND:.0e}; record "
        f"{REC['e_qd']:.3e}, {D18}:68)")
    say(f"[theta] quadratic duty theta_max at t* {fmt(tq['theta_max_k'])} K (record "
        f"{REC['theta_max_qd']:.5f}); ratio to mri:31's unscheduled CN theta_max "
        f"{REC['theta_max_cn_unsched']:.4f} K: {ratio_unsched:.4f} (printed only; predicted 0.4917, "
        f"range [0.45, 0.55]; in range {0.45 <= ratio_unsched <= 0.55})")
    assert tq["source"] == "transverse_quadratic" and tq["thermal_degree"] == d.DEGREE
    assert d._hex_dict_equal(tq["duty"], {"period_s": 2.0 * d.HALF_T, "duty": d.DUTY})
    assert len(tq["source_multipliers"]) == 41 and s_json == d.SCHEDULE
    assert p_json <= L0_BAND, f"(P-json) {p_json:.3e}"
    assert tq["max_balance_residual"] <= E_BAND

    # ---- quadratic duty fields to ParaView (P1, real parts) ------------------------------------
    f_qfile = p1_field(msh, "theta_quadratic_duty_tstar", th.temperature, d.T_INF)
    f_qser = p1_field(msh, "theta_quadratic_duty_series_tstar", ex_d)
    f_qdif = diff_field(msh, "theta_quadratic_duty_minus_series", f_qfile, f_qser)
    f_unsched = p1_field(msh, "theta_unscheduled_series_tstar", ex_s)
    b2c2 = d.HALF[1] ** 2 + d.HALF[2] ** 2
    f_qsrc = p1_field(msh, "q_transverse_quadratic",
                      lambda x: 3.0 * d.Q * (x[1] ** 2 + x[2] ** 2) / b2c2)
    ft = th.facet_tags
    ft.name = "robin_facets"
    qd_fields = {f.name: f for f in (f_qfile, f_qser, f_qdif, f_unsched, f_qsrc)}
    xdmf_qd, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_quadratic_duty_combined", msh, ones_tags(msh), qd_fields, comm=comm,
        facet_tags=ft)

    # ==== ball m19 leg: the gate's W1 body (b:278-329) =========================================
    t01 = time.perf_counter()
    d._fresh_dir(RUN_M19, comm)
    rc_m, cap_m = cli_run(b.M19, RUN_M19)
    t_file_m = time.perf_counter() - t01
    assert d._rc_all(rc_m, comm) == (0, 0) and len(cap_m) == 1, (rc_m, len(cap_m))
    run_m = cap_m[0]
    pr = run_m.probes
    qcf = b._q_cf()
    pym = b._py_m19().probe([b.WHOLE_BALL, *b.BALL_PROBES, b.DISP_PROBE],
                            b.RobinCooling(b.K, b.H, b.T_INF),
                            transient=b.TransientSpec(b.T_STAR / b.N_MAIN, b.N_MAIN, b.RHO_C))
    m_f = pr.theta.function_space.mesh
    m_p = pym.theta.function_space.mesh
    cam, cbm, geom = b._same_mesh(m_f, m_p, comm)

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    write_setup_figure(
        m_f, ones_tags(m_f), FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel ball (R = 10 mm)"},
        title="mri:32 -- gel ball from femem run: probe balls read back from the JSON summary; "
              "duty-cycled quadratic box beside, z = 0 slice",
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        clip_normal=(0.0, 0.0, -1.0), comm=comm,
    )

    ipm = b._theta_owned_rel(b._on(m_p, pr.thermal.temperature), pym.thermal.temperature, comm)
    two = b._theta_owned_rel(b._on(run_m.thermal.temperature.function_space.mesh,
                                   pr.thermal.temperature), run_m.thermal.temperature, comm)
    vals_f = [float(v) for v in pr.values]
    vals_p = [float(v) for v in pym.values]
    say()
    _p0(f"WF1S19 (UT) m19: cells {cam} / {cbm} (U60 {b.U60_BALL_CELLS}; equal "
        f"{cam == b.U60_BALL_CELLS}); geometry equal {geom}; file run {t_file_m:.1f}s (record "
        f"{REC['m19_file_s']} s, real build, {B19}:80)")
    _p0(f"WF1S19 (I-Py) m19: owned max|dT|/max|T-Tinf| = {ipm:.3e} (band {L0_BAND:.0e}; record "
        f"{REC['i_py_m19']:.3e}, {B19}:81 real); temperature block vs probe re-solve {two:.3e}")
    _p0("WF1S19 probe values file " + ", ".join(f"{v:.6g}" for v in vals_f))
    _p0("WF1S19 probe values py   " + ", ".join(f"{v:.6g}" for v in vals_p) + " (printed only)")
    mean_h, v_h = b._domain_mean(pr.theta, comm)
    a_val = abs(vals_f[0] - mean_h) / abs(mean_h)
    a_vol = abs(float(pr.volumes[0]) / v_h - 1.0)
    _p0(f"WF1S19 (A-ball) |P0 - mean| / |mean| = {a_val:.3e}, |V0/V_h - 1| = {a_vol:.3e} (band "
        f"{B_BAND:.0e}; U60 {b.U60_A_BALL:.3e}; record {REC['a_ball'][0]:.3e} / "
        f"{REC['a_ball'][1]:.3e}, {B19}:84)")
    ratios, digits, comps = [], [], []
    for j, c in enumerate(b.CENTRES):
        cf = b._ball_mean(float(np.linalg.norm(c)), b.A_PROBE, b.T_STAR, qcf)
        comps.append(cf)
        r = abs(vals_f[1 + j] / cf - 1.0)
        ratios.append(r)
        digits.append(f"{r:.4e}")
        _p0(f"WF1S19 (T*-ball) c{j}: P = {vals_f[1 + j]:.6g}, ball mean = {cf:.6g}, rel = "
            f"{r:.4e} (U60 {b.U60_T_BALL[j]}; band {P_BAND:.0e})")
    cf2 = comps[2]
    c_disp = abs(vals_f[-1] / cf2 - 1.0)
    sep = abs(b._ball_mean(float(np.linalg.norm(b.DISPLACED)), b.A_PROBE, b.T_STAR, qcf) / cf2 - 1.0)
    _p0(f"WF1S19 control (c-disp) displaced vs c2's comparand: miss {c_disp:.4e} (U60 "
        f"{b.U60_C_DISP:.4e}; equal digits {f'{c_disp:.4e}' == f'{b.U60_C_DISP:.4e}'}; comparand "
        f"separation {sep:.4e}, record {REC['sep']:.4e}; asserted > {P_BAND:.0e})")
    b._check_summary(run_m, comm, "m19")
    assert cam == cbm and geom
    assert ipm <= L0_BAND, f"(I-Py) {ipm:.3e}"
    assert two <= L0_BAND, f"temperature block vs probe {two:.3e}"
    assert a_val <= B_BAND and a_vol <= B_BAND
    for r in ratios:
        assert r <= P_BAND
    # substitution 3: the real-build digit assertion is a printed equality in the complex build
    say(f"[(T*-ball) digits] complex-build digits {tuple(digits)}; equal U60's "
        f"{tuple(digits) == b.U60_T_BALL} (printed, not asserted)")
    assert c_disp > P_BAND
    t_m19 = time.perf_counter() - t01
    _p0(f"WF1S19 W1 m19 leg {t_m19:.1f}s (record {REC['m19_leg_s']} s, real build, {B19}:90)")

    # ---- the ball probe table, read from the JSON on disk --------------------------------------
    disk_m = read_summary(run_m, comm)
    js_pr = disk_m["probes"]
    jv = js_pr["values_k"]
    rj, cj = [], []
    for j, c in enumerate(b.CENTRES):
        cf = b._ball_mean(float(np.linalg.norm(c)), b.A_PROBE, b.T_STAR, qcf)
        cj.append(cf)
        rj.append(abs(jv[1 + j] / cf - 1.0))
    a_json = abs(jv[0] - mean_h) / abs(mean_h)
    a_json_vol = abs(js_pr["volumes_m3"][0] / v_h - 1.0)
    c_disp_json = abs(jv[4] / cf2 - 1.0)
    rj_digits = tuple(f"{r:.4e}" for r in rj)
    mults_m_json = disk_m["temperature"]["source_multipliers"]
    mults_m_mem = [float(m) for m in run_m.thermal.source_multipliers]
    say()
    say(f"[(S-json)] ball m19: family {disk_m['family']!r}, 'probes' present {'probes' in disk_m}; "
        f"temperature keys {sorted(disk_m['temperature'])} (no 'source' key, printed); probes keys "
        f"{sorted(js_pr)}")
    say(f"[(R-json)] (A-ball) from the JSON values_k[0]: {a_json:.3e}, volumes_m3[0]/V_h - 1 = "
        f"{a_json_vol:.3e} (band {B_BAND:.0e})")
    for j, c in enumerate(b.CENTRES):
        say(f"[(R-json)] c{j} {c}: JSON ball value {jv[1 + j]:.5g} K vs ball mean ratio {rj[j]:.4e} "
            f"(band {P_BAND:.0e}; U60 {b.U60_T_BALL[j]})")
    say(f"[(R-json)] (c-disp) from the JSON values_k[4] vs c2's comparand: {c_disp_json:.4e} (asserted > "
        f"{P_BAND:.0e}; U60 {b.U60_C_DISP:.4e})")
    say(f"[(R-json) digits] from disk {rj_digits}; equal U60's {rj_digits == b.U60_T_BALL} "
        "(printed, not asserted)")
    say(f"[(M-json)] m19 source_multipliers on disk == in-memory {mults_m_json == mults_m_mem}; "
        f"length {len(mults_m_json)} (printed, not asserted)")
    c_m19 = json_hex(js_pr) == b._probe_hex([b.WHOLE_BALL, *b.BALL_PROBES, b.DISP_PROBE])
    say(f"[(C-json)] probes.centres / probes.radii by float.hex == [WHOLE_BALL, *BALL_PROBES, "
        f"DISP_PROBE] {c_m19}")
    assert disk_m["family"] == "thermal_ball" and "probes" in disk_m
    assert a_json <= B_BAND and a_json_vol <= B_BAND
    for r in rj:
        assert r <= P_BAND
    assert c_disp_json > P_BAND
    assert mults_m_json == mults_m_mem
    assert c_m19

    # ---- m19 fields to ParaView (P1 on the probe mesh) -----------------------------------------
    f_mfile = p1_field(m_f, "theta_ball_probe_file", pr.theta)
    f_mpy = p1_field(m_f, "theta_ball_python", b._on(m_f, pym.theta))
    cen = np.array([*b.CENTRES, b.DISPLACED], dtype=float)

    def balls(x):
        dist = np.sqrt((x[0][:, None] - cen[None, :, 0]) ** 2 + (x[1][:, None] - cen[None, :, 1]) ** 2
                       + (x[2][:, None] - cen[None, :, 2]) ** 2)
        return (dist.min(axis=1) <= b.A_PROBE).astype(float)

    f_balls = p1_field(m_f, "probe_balls", balls)
    m19_fields = {f.name: f for f in (f_mfile, f_mpy, f_balls)}
    xdmf_m, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_ball_m19_combined", m_f, ones_tags(m_f), m19_fields, comm=comm)

    # ==== ball duty leg: the gate's duty body (b:334-370) =======================================
    t02 = time.perf_counter()
    d._fresh_dir(RUN_BD, comm)
    rc_b, cap_b = cli_run(b.DUTY, RUN_BD)
    assert d._rc_all(rc_b, comm) == (0, 0) and len(cap_b) == 1, (rc_b, len(cap_b))
    run_b = cap_b[0]
    prd = run_b.probes
    pyd = b._py_duty().probe(
        b.DUTY_PROBES, b.RobinCooling(b.K9, b.H9, b.T_INF9),
        transient=b.TransientSpec(150.0, b._n(b.T_END_A, 150.0), b.RHO_C9),
        duty=b.DutyCycle(*b.CASES["A"]),
    )
    m_fd = prd.theta.function_space.mesh
    m_pd = pyd.theta.function_space.mesh
    cad, cbd, geod = b._same_mesh(m_fd, m_pd, comm)
    ipd = b._theta_owned_rel(b._on(m_pd, prd.thermal.temperature), pyd.thermal.temperature, comm)
    multd = np.array_equal(np.asarray(prd.thermal.source_multipliers),
                           np.asarray(pyd.thermal.source_multipliers))
    mult_sd = np.array_equal(np.asarray(run_b.summary["temperature"]["source_multipliers"]),
                             np.asarray(pyd.thermal.source_multipliers))
    mean_d, v_d = b._domain_mean(prd.theta, comm)
    vals_d = [float(v) for v in prd.values]
    a_vald = abs(vals_d[0] - mean_d) / abs(mean_d)
    a_vold = abs(float(prd.volumes[0]) / v_d - 1.0)
    say()
    _p0(f"WF1S19 (I-du) duty: cells {cad} / {cbd}; geometry equal {geod}; owned max|dT|/max|T-Tinf| "
        f"= {ipd:.3e} (band {L0_BAND:.0e}; record {REC['i_du']:.3e}, {B19}:93); multipliers equal "
        f"{multd} / summary {mult_sd} ({len(pyd.thermal.source_multipliers)} entries)")
    _p0(f"WF1S19 (I-du) whole ball |P0 - mean| / |mean| = {a_vald:.3e}, |V0/V_h - 1| = "
        f"{a_vold:.3e} (band {B_BAND:.0e}; record {REC['a_du'][0]:.3e} / {REC['a_du'][1]:.3e}, "
        f"{B19}:94)")
    _p0("WF1S19 (I-du) probe values file " + ", ".join(f"{v:.6g}" for v in vals_d))
    _p0("WF1S19 (I-du) probe values py   " + ", ".join(f"{float(v):.6g}" for v in pyd.values)
        + " (printed; no comparand is posed)")
    b._check_summary(run_b, comm, "duty")
    assert cad == cbd and geod
    assert ipd <= L0_BAND, f"(I-du) {ipd:.3e}"
    assert multd and mult_sd
    assert a_vald <= B_BAND and a_vold <= B_BAND
    assert run_b.summary["temperature"]["duty"] == {"period_s": b.CASES["A"][0],
                                                    "duty": b.CASES["A"][1]}
    t_bd = time.perf_counter() - t02
    _p0(f"WF1S19 W1 duty leg {t_bd:.1f}s (record {REC['duty_leg_s']} s, real build, {B19}:98)")

    disk_b = read_summary(run_b, comm)
    tb = disk_b["temperature"]
    mults_b_json = tb["source_multipliers"]
    mults_b_mem = [float(m) for m in run_b.thermal.source_multipliers]
    c_du = json_hex(disk_b["probes"]) == b._probe_hex(b.DUTY_PROBES)
    say(f"[(S-json)] ball duty: family {disk_b['family']!r}, 'probes' present {'probes' in disk_b}; "
        f"duty {tb['duty']!r}")
    say(f"[(M-json)] ball duty source_multipliers on disk: {len(mults_b_json)} entries == in-memory "
        f"{mults_b_json == mults_b_mem}")
    say(f"[(C-json)] ball duty probes.centres / probes.radii by float.hex == DUTY_PROBES {c_du}")
    say(f"[balance] ball duty max_balance_residual {tb['max_balance_residual']:.3e} (printed only; "
        f"predicted ~3e-14 from {E123}:45 / record {REC['bal_du']:.3e}, range [0, 1e-11]; band "
        f"{E_BAND:.0e})")
    say("[duty values_k from disk] " + ", ".join(f"{v:.6g}" for v in disk_b["probes"]["values_k"])
        + f" K (printed only; no comparand posed; record row in {B19}:95)")
    assert disk_b["family"] == "thermal_ball" and "probes" in disk_b
    assert tb["duty"] == {"period_s": b.CASES["A"][0], "duty": b.CASES["A"][1]}
    assert len(mults_b_json) == 95 and mults_b_json == mults_b_mem
    assert c_du

    f_dfile = p1_field(m_fd, "theta_ball_duty_file_tend", prd.theta)
    xdmf_b, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_ball_duty_combined", m_fd, ones_tags(m_fd),
        {f_dfile.name: f_dfile}, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- CSV: the summary read-back, rank 0 ------------------------------------------------------
    csv_path = OUTPUT_DIR / f"{BASENAME}_summary_readback.csv"
    rows = []
    for label, disk, mem in (("quadratic_duty", disk_qd, run.summary),
                             ("ball_m19", disk_m, run_m.summary),
                             ("ball_duty", disk_b, run_b.summary)):
        pairs = [("cells", disk["cells"], mem["cells"])]
        for key, v in disk["temperature"].items():
            if isinstance(v, float):
                pairs.append((key, v, mem["temperature"][key]))
        for key, v_json, v_mem in pairs:
            is_t = key == "theta_max_k"
            sj = fmt(v_json) if is_t else f"{v_json:.6g}"
            sm = fmt(v_mem) if is_t else f"{v_mem:.6g}"
            rows.append((label, key, sj, sm, hexeq(v_json, v_mem)))
    if rank0:
        with open(csv_path, "w") as fh:
            fh.write("file,key,json_value,in_memory_value,hex_equal\n")
            for r in rows:
                fh.write(",".join(str(x) for x in r) + "\n")

    # ---- CSV: the probe table, rank 0 ------------------------------------------------------------
    pt_path = OUTPUT_DIR / f"{BASENAME}_probe_table.csv"
    prow = []
    for label, js in (("ball_m19", js_pr), ("ball_duty", disk_b["probes"])):
        for i, (c, rad, val, vol) in enumerate(zip(js["centres"], js["radii"], js["values_k"],
                                                   js["volumes_m3"])):
            comp_s = ratio_s = u60 = eq = ""
            if label == "ball_m19":
                if i == 0:
                    comp_s, ratio_s = f"{mean_h:.6g}", f"{a_json:.4e}"
                elif i <= 3:
                    comp_s, ratio_s = f"{cj[i - 1]:.6g}", f"{rj[i - 1]:.4e}"
                    u60 = b.U60_T_BALL[i - 1]
                    eq = ratio_s == u60
                else:
                    comp_s, ratio_s = f"{cf2:.6g}", f"{c_disp_json:.4e}"
                    u60 = f"{b.U60_C_DISP:.4e}"
                    eq = ratio_s == u60
            prow.append((label, i, "(" + " ".join(f"{x:.6g}" for x in c) + ")", f"{rad:.6g}",
                         f"{val:.6g}", f"{vol:.6g}", comp_s, ratio_s, u60, eq))
    if rank0:
        with open(pt_path, "w") as fh:
            fh.write("file,probe,centre_m,radius_m,value_k,volume_m3,comparand_k,ratio,u60_digit,"
                     "digits_equal\n")
            for r in prow:
                fh.write(",".join(str(x) for x in r) + "\n")
    comm.Barrier()
    assert all(r[4] for r in rows), [r for r in rows if not r[4]]

    say()
    say(f"[output] {xdmf_qd.name if rank0 else ''}: P1 theta_quadratic_duty_tstar, "
        "theta_quadratic_duty_series_tstar, theta_quadratic_duty_minus_series, "
        "theta_unscheduled_series_tstar, q_transverse_quadratic, Robin facets")
    say(f"[output] {xdmf_m.name if rank0 else ''}: P1 theta_ball_probe_file, theta_ball_python, "
        "probe_balls")
    say(f"[output] {xdmf_b.name if rank0 else ''}: P1 theta_ball_duty_file_tend")
    say(f"[output] {csv_path.name} ({len(rows)} rows, every JSON value hex-equal to the in-memory "
        "summary)")
    say(f"[output] {pt_path.name} ({len(prow)} rows)")
    say(f"[time] quadratic duty leg {t_qd:.1f} s; ball m19 leg {t_m19:.1f} s; ball duty leg "
        f"{t_bd:.1f} s; total {time.perf_counter() - started:.1f} s")
    say("[scope] one duty-cycled transverse-quadratic thermal_box file (m = 4, P2) and two ball-probe "
        "thermal_ball files through cli.main plus one malformed: complex build, -n 2; summaries and the "
        "ball probe table read back from disk; no file carries both excitation.source and "
        "outputs (not gated) and no ball file a wall; the ball-probe leg is an example record in the "
        "complex build, not a gate; no EM, no coil, no compliance, no external check claimed")
    say("All assertions hold")


if __name__ == "__main__":
    main()
