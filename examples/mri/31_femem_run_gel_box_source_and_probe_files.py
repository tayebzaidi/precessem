"""Example (`EX-135`, `mri:31`): ``femem run`` on the tracked F2182-sized gel-box ``thermal_box`` files
that carry a **prescribed non-uniform heat source** (``excitation.source: {kind: transverse_quadratic}``,
`WF-1` step 16, §9 item 528) and the files that carry a **probe table** (``outputs.probes``, `WF-1`
step 17, §9 item 539), with the quadratic-source temperature and the five ball-mean readings read back
from the run's JSON summary on disk.

No example before this one set a non-uniform source in a file or read a ``probes`` block back from a
run.  ``mri:29`` runs uniform-source files with no ``outputs``; ``mat:24`` / ``mat:25`` build Python
objects.  This one runs ``tests/package/wf1/thermal_box_quadratic_steady_m4.yaml``,
``thermal_box_quadratic_cn_m4.yaml``, ``thermal_box_uniform_steady_m4_p2.yaml`` (the P2 control),
``thermal_box_steady_m4_probes.yaml``, ``thermal_box_wrong_wall_m4_probes.yaml`` and
``thermal_box_duty_m4_probes.yaml`` in place through ``cli.main(["run", ...])`` in-process (the
``mri:21`` / ``mri:29`` call-through wrapper on ``fem_em_solver.study_yaml.run_study_file``, restored in
``finally``).  **No tracked file carries both** ``excitation.source`` and ``outputs``; that combination
is not gated and is not run here.

The four solve-node bodies of the two gate modules (``test_wf1_yaml_thermal_box_quadratic``: steady and
CN; ``test_wf1_yaml_thermal_box_probes``: steady and duty) are copied with four substitutions
(``run_study_file`` -> the captured ``cli_run``; ``LOGS`` -> this example's run directories; the
real-build digit equality of (R-box) -> a printed equality; full-``repr`` probe-value prints -> ``.6g``).
Every record, band and helper is imported from the gate modules (the ANS-1 rule); record digits sit only
in the example-local ``REC`` dict with their ``log:line``.

The quadratic leg is a complex-build gate and runs in the build of record.  The probe leg's gate ran in
the real build; here it runs in the complex build, so its **bands** are asserted and its **digits** are
printed beside U60's and ``EX-130``'s complex readings (the ``EX-130`` precedent).  These are example
records, not a complex gate of the `WF-1` step 17 route.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:31 -n 2 -t 590

Not claimed: the quadratic source together with probes in one file (ungated); a duty-cycled transverse
source; any comparand on a transient or duty probe; probes at P2; scalar point values; a complex gate of
the probe-file route; the real build; the installed console script, ``--ranks`` or ``--estimate``; EM,
a coil or a box EM SAR; perfusion; bundle reproducibility (``mri:27``); F2182 compliance; an external
check.
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
import tests.package.test_wf1_yaml_thermal_box_probes as p  # noqa: E402
import tests.package.test_wf1_yaml_thermal_box_quadratic as q  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_31_femem_run_gel_box_source_and_probe_files"
RUN_QS = OUTPUT_DIR / "mri_31_run_quadratic_steady"
RUN_UP2 = OUTPUT_DIR / "mri_31_run_uniform_p2"
RUN_QCN = OUTPUT_DIR / "mri_31_run_quadratic_cn"
RUN_PS = OUTPUT_DIR / "mri_31_run_probes_steady"
RUN_PW = OUTPUT_DIR / "mri_31_run_probes_wrong_wall"
RUN_PD = OUTPUT_DIR / "mri_31_run_probes_duty"
BAD_DIR = OUTPUT_DIR / "mri_31_bad"

Q1A = "20261019T003608Z_WF-1-s16-W1a-steady.log"
Q1B = "20261019T003738Z_WF-1-s16-W1b-cn.log"
Q0 = "20261019T003542Z_WF-1-s16-W0.log"
P1L = "20261019T033412Z_WF-1-s17-W1-probes.log"
P0L = "20261019T033322Z_WF-1-s17-W0.log"
E130 = "20261018T214035Z_EX-130.log"

# Records, printed beside readings only (the copied bodies assert the .6e strings of P36 / step 3).
REC = {
    "q_file_s": 11.5, "q_leg_s": 74.9,                      # Q1a:63, :71
    "i_py_q": 1.032e-15,                                     # Q1a:64
    "gauss": 9.775e-14, "src_power_rel": 1.157e-13,          # Q1a:66
    "c_srcp": 5.1471e-01, "c_src": 4.4856e-01, "c_z": 4.7893e-02,  # Q1a:68-70
    "cn_file_s": 32.5, "cn_leg_s": 95.6,                     # Q1b:63, :68
    "i_py_cn": 1.213e-14,                                    # Q1b:64
    "e_quad": 2.288e-13, "c_src_t": 6.9562e-01,              # Q1b:66-67
    "pr_file_s": 6.6, "pr_leg_s": 13.9,                      # P1:61, :74
    "i_py_pr": 1.184e-15,                                    # P1:62 (E130:63 complex)
    "a_box": (2.495e-13, 2.509e-13),                         # P1:65 (E130:64)
    "c_wall": 3.5712e-01,                                    # P1:71 (E130:82 complex)
    "i_tr": 6.912e-14, "duty_leg_s": 8.7,                    # P1:77, :81
}
# E130:65-69, complex Study.probe digits (equal U60's p.U60_R_BOX).
E130_R_BOX = ("2.6992e-03", "2.6936e-03", "2.7066e-03", "2.6468e-03", "2.6443e-03")


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
        rc = q.cli.main(["run", str(path), "--output-dir", str(out_dir)])
    finally:
        study_yaml.run_study_file = original
    assert study_yaml.run_study_file is original
    return rc, captured


def cli_check(path):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = q.cli.main(["check", str(path)])
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


def diff_field(msh, name, a, b):
    V1 = fem.functionspace(msh, ("Lagrange", 1))
    f = fem.Function(V1, name=name)
    f.x.array[:] = np.real(a.x.array) - np.real(b.x.array)
    f.x.scatter_forward()
    return f


def main() -> None:
    comm = MPI.COMM_WORLD
    rank0 = comm.rank == 0
    started = time.perf_counter()

    def say(msg=""):
        if rank0:
            print(msg, flush=True)

    _p0 = say  # the copied bodies print through `_p0`

    # names the copied bodies use unqualified (all imported through the gate modules)
    Q, HALF, T_INF, K_GEL, H_GEL = q.Q, q.HALF, q.T_INF, q.K_GEL, q.H_GEL
    L0_BAND, INTERIOR_BAND, E_BAND, S2_GAUSS_BAND = q.L0_BAND, q.INTERIOR_BAND, q.E_BAND, q.S2_GAUSS_BAND
    B_BAND = p.B_BAND
    P36_S1, P36_T, P36_E, P36_S2, P36_DOFS = q.P36_S1, q.P36_T, q.P36_E, q.P36_S2, q.P36_DOFS
    P36_C_SRCP, P36_C_SRC, P36_C_Z = q.P36_C_SRCP, q.P36_C_SRC, q.P36_C_Z
    STEP3_M4, DEGREE, N_MAIN, T_STAR, RHO_C = q.STEP3_M4, q.DEGREE, q.N_MAIN, q.T_STAR, q.RHO_C
    ROBIN_TAG, SERIES_TERMS, WEIGHTS_NOZ = q.ROBIN_TAG, q.SERIES_TERMS, q.WEIGHTS_NOZ
    _owned_rel, _ndofs, _on_p, _all = q._owned_rel, q._ndofs, q._on_p, q._all
    _exact_p, _rel_p, _interior_chi, _allsum = q._exact_p, q._rel_p, q._interior_chi, q._allsum
    _fresh_dir, _mapping, _rc_all = q._fresh_dir, q._mapping, q._rc_all
    _cooling, _quad_study = q._cooling, q._quad_study
    _steady_ex, _trans_ex, _quad_ex, _quad_trans_ex = q._steady_ex, q._trans_ex, q._quad_ex, q._quad_trans_ex
    TransientSpec, RobinCooling, DutyCycle = q.TransientSpec, q.RobinCooling, p.DutyCycle
    U60_R_BOX, U60_C_WALL, U60_A_BOX = p.U60_R_BOX, p.U60_C_WALL, p.U60_A_BOX
    BOX_CENTRES, BOX_PROBES, WHOLE_BOX, A_BOXPROBE = p.BOX_CENTRES, p.BOX_PROBES, p.WHOLE_BOX, p.A_BOXPROBE
    _domain_mean, _series_at, _same_mesh = p._domain_mean, p._series_at, p._same_mesh
    _on, _theta_owned_rel, _py_study = p._on, p._theta_owned_rel, p._py_study
    DT_BOX, HALF_T, N_STEPS = p.DT_BOX, p.HALF_T, p.N_STEPS

    # ---- (K) structure and femem check --------------------------------------------------------
    assert q.IS_COMPLEX, "complex DolfinX build required (the runner's `mri:` group sources it)"
    assert comm.size == 2, f"comm.size {comm.size}"
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)
    say("=" * 72)
    say("EX-135 -- femem check / run on the transverse-quadratic and ball-probe thermal_box files "
        "(mri:31)")
    say("=" * 72)
    say(f"[records] beside {Q1A}, {Q1B}, {P1L} (real build), {E130} (complex); ranks {comm.size}")
    say(f"[bands] imported: L0_BAND {L0_BAND:.0e}, INTERIOR_BAND {INTERIOR_BAND:.0e}, "
        f"E_BAND {E_BAND:.0e}, S2_GAUSS_BAND {S2_GAUSS_BAND:.0e}, B_BAND {B_BAND:.0e}")

    check_rows = {}
    for f, excit in ((q.QUAD_STEADY, "TransverseQuadraticHeatSource"),
                     (q.QUAD_CN, "TransverseQuadraticHeatSource"),
                     (q.UNIFORM_P2, "UniformHeatSource"),
                     (p.STEADY, None), (p.WRONG_WALL, None), (p.DUTY, None)):
        rc, lines = cli_check(f)
        rc_all = _rc_all(rc, comm)
        ok = (comm.rank != 0) or (
            [ln for ln in lines if ln.startswith("family: ")] == ["family: thermal_box"]
            and (excit is None or f"excitation: {excit}" in lines)
        )
        ok = _all(comm, ok)
        exline = [ln for ln in lines if ln.startswith("excitation: ")] if rank0 else []
        check_rows[f.name] = (rc_all, ok)
        say(f"[(K)] femem check {f.name}: rc {rc_all}; one 'family: thermal_box'"
            f"{'' if excit is None else ' and ' + repr('excitation: ' + excit)} on rank 0: {ok}; "
            f"rank-0 excitation line {exline}")
    for name, (rc_all, ok) in check_rows.items():
        assert rc_all == (0, 0), (name, rc_all)
        assert ok, name

    # (F-diff) the files differ exactly where the physics says
    m_qs, m_qc, m_up, m_gs = (_mapping(f) for f in (q.QUAD_STEADY, q.QUAD_CN, q.UNIFORM_P2, g.STEADY))
    m_ps, m_pw, m_pd = (_mapping(f) for f in (p.STEADY, p.WRONG_WALL, p.DUTY))
    diffs = {
        "quad-steady vs uniform-p2": (flat_diff(m_qs, m_up), {"excitation.source"}),
        "quad-steady vs quad-cn": (flat_diff(m_qs, m_qc), {"thermal.transient"}),
        "uniform-p2 vs uniform-m4": (flat_diff(m_up, m_gs), {"thermal.degree"}),
        "probes-steady vs uniform-m4": (flat_diff(m_ps, m_gs), {"outputs"}),
        "probes-steady vs wrong-wall": (flat_diff(m_ps, m_pw),
                                        {"thermal.h_w_per_m2_k", "outputs.probes"}),
        "probes-steady vs duty": (flat_diff(m_ps, m_pd),
                                  {"thermal.transient", "thermal.duty", "outputs.probes"}),
    }
    for label, (got, want) in diffs.items():
        say(f"[(F-diff)] {label}: {sorted(got)}")
    for label, (got, want) in diffs.items():
        assert got == want, (label, sorted(got), sorted(want))
    # (F-none) no file the example runs carries both capabilities
    all_maps = {f.name: _mapping(f) for f in (q.QUAD_STEADY, q.QUAD_CN, q.UNIFORM_P2, p.STEADY,
                                              p.WRONG_WALL, p.DUTY, p.MALFORMED)}
    both = [n for n, m in all_maps.items()
            if "source" in m.get("excitation", {}) and "outputs" in m]
    say(f"[(F-none)] files carrying both excitation.source and outputs: {both} "
        f"(over {len(all_maps)} mappings; expected [])")
    assert both == []
    src_files = sorted(n for n, m in all_maps.items() if "source" in m.get("excitation", {}))
    out_files = sorted(n for n, m in all_maps.items() if "outputs" in m)
    say(f"[(F-none)] excitation.source in {src_files}; outputs in {out_files}")

    # (cX) the malformed probe file: refusal, rc (2, 2), nothing written (copied p:183-196)
    assert p._field(lambda: p.load_study(p.MALFORMED)) == "outputs.probes[1].radius"
    cx_field = p._field(lambda: p.load_study(p.MALFORMED))
    _fresh_dir(BAD_DIR, comm)
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        rc_run = q.cli.main(["run", str(p.MALFORMED), "--output-dir", str(BAD_DIR)])
        rc_check = q.cli.main(["check", str(p.MALFORMED)])
    run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
    exists = comm.allreduce(bool(BAD_DIR.exists()), op=MPI.LOR)
    named = comm.allreduce(int("outputs.probes" in err.getvalue()), op=MPI.MAX)
    say(f"[(cX)] {p.MALFORMED.name}: field {cx_field!r}; run rc {run_all}, check rc {check_all}; out dir "
        f"exists {exists}; 'outputs.probes' on stderr of some rank {bool(named)} (printed, not "
        f"asserted; expected 'outputs.probes[1].radius', (2, 2), (2, 2), False; {P0L}:71)")
    assert cx_field == "outputs.probes[1].radius"
    assert run_all == (2, 2) and check_all == (2, 2)
    assert not exists

    q_v = Q * 8.0 * float(np.prod(HALF))  # Q * V = 49.14 W

    # ==== quadratic steady leg: the gate's W1a body (q:270-338) ================================
    t00 = time.perf_counter()
    out = RUN_QS
    _fresh_dir(out, comm)
    rc, cap = cli_run(q.QUAD_STEADY, out)
    t_file = time.perf_counter() - t00
    assert _rc_all(rc, comm) == (0, 0) and len(cap) == 1, (rc, len(cap))
    run = cap[0]
    th = run.thermal
    msh = th.submesh

    # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
    cell_tags_q = ones_tags(msh)
    write_setup_figure(
        msh, cell_tags_q, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "gel (650 x 420 x 90 mm)"},
        title="mri:31 -- F2182-sized gel box from femem run: transverse-quadratic q (P2) file and "
              "ball-probe file, z = 0 slice",
        slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
        clip_normal=(0.0, 0.0, -1.0), comm=comm,
    )

    py = _quad_study().temperature(_cooling())
    T_py = _on_p(msh, py.temperature)
    ip = _owned_rel(T_py, th.temperature, comm)
    out_u = RUN_UP2
    _fresh_dir(out_u, comm)
    rc_u, cap_u = cli_run(q.UNIFORM_P2, out_u)
    assert _rc_all(rc_u, comm) == (0, 0) and len(cap_u) == 1, (rc_u, len(cap_u))
    run_u = cap_u[0]
    T_u = _on_p(msh, run_u.thermal.temperature)
    t_solves = time.perf_counter() - t00
    dofs = [_ndofs(th.temperature), _ndofs(py.temperature), _ndofs(run_u.thermal.temperature)]

    th_q = th.temperature - T_INF
    th_u = T_u - T_INF
    chi = _interior_chi(msh)
    # (S2-quad) P36's Gauss reading
    ds = ufl.Measure("ds", domain=msh, subdomain_data=th.facet_tags)(ROBIN_TAG)
    conv = float(np.real(_allsum(comm, H_GEL * th_q * ds)))
    sp = float(np.real(th.source_power))
    gauss = abs(conv - sp) / sp
    # (U) step 3's digits through the uniform P2 file
    ex_u, npts = _exact_p(msh, _steady_ex())
    u_whole = _rel_p(msh, th_u, ex_u)
    u_int = _rel_p(msh, th_u, ex_u, chi)
    # (S1-quad)
    ex_q, _ = _exact_p(msh, _quad_ex())
    s1 = _rel_p(msh, th_q, ex_q, chi)
    # controls
    c_srcp = _rel_p(msh, th_u, ex_q, chi)
    c_src = _rel_p(msh, th_q, ex_u, chi)
    ex_z, _ = _exact_p(msh, _quad_ex(WEIGHTS_NOZ, SERIES_TERMS))
    c_z = _rel_p(msh, th_q, ex_z, chi)

    say()
    _p0(f"WF1S16 (Pdeg) steady P2 dofs file / Python / uniform file {dofs} (P36 {P36_DOFS}; "
        f"{Q1A}:63); file run {t_file:.1f}s (record {REC['q_file_s']} s), three solves "
        f"{t_solves:.1f}s; summary source "
        f"{run.summary['temperature']['source']!r} / {run_u.summary['temperature']['source']!r}")
    _p0(f"WF1S16 (I-Py) steady: owned max|dtheta|/max|theta| = {ip:.3e} (hex {float(ip).hex()}; band "
        f"{L0_BAND:.0e}; record {REC['i_py_q']:.3e}, {Q1A}:64)")
    _p0(f"WF1S16 (S1-quad) steady interior vs quadratic series {s1:.6e} (P36 {P36_S1}; band "
        f"{INTERIOR_BAND:.0e}); {npts} P3 points")
    _p0(f"WF1S16 (S2-quad) Gauss |h oint theta - source_power|/source_power = {gauss:.3e} (P36 "
        f"{P36_S2:.3e}; band {S2_GAUSS_BAND:.0e}); source_power/(Q V) - 1 = "
        f"{sp / (Q * 8.0 * HALF[0] * HALF[1] * HALF[2]) - 1.0:+.3e} (record "
        f"+{REC['src_power_rel']:.3e}, {Q1A}:66)")
    _p0(f"WF1S16 (U) uniform P2 file: whole {u_whole:.6e} (step 3 {STEP3_M4['whole']}); interior "
        f"{u_int:.6e} (step 3 {STEP3_M4['interior']})")
    _p0(f"WF1S16 control (c-src') uniform file vs quadratic series {c_srcp:.4e} (P36 "
        f"{P36_C_SRCP:.4e}; asserted > {INTERIOR_BAND:.0e})")
    _p0(f"WF1S16 control (c-src) quadratic file vs uniform series {c_src:.4e} (P36 "
        f"{P36_C_SRC:.4e}; asserted > {INTERIOR_BAND:.0e})")
    _p0(f"WF1S16 control (c-z) quadratic file vs no-z series {c_z:.4e} (P36 {P36_C_Z:.4e}; "
        f"printed only, ratio to band {c_z / INTERIOR_BAND:.2f})")
    assert dofs == [P36_DOFS] * 3
    assert run.summary["temperature"]["source"] == "transverse_quadratic"
    assert run_u.summary["temperature"]["source"] == "uniform"
    assert run.summary["temperature"]["thermal_degree"] == DEGREE
    assert ip <= L0_BAND, f"(I-Py) steady {ip:.3e}"
    assert s1 <= INTERIOR_BAND
    assert f"{s1:.6e}" == P36_S1
    assert gauss <= S2_GAUSS_BAND
    assert f"{u_whole:.6e}" == STEP3_M4["whole"]
    assert f"{u_int:.6e}" == STEP3_M4["interior"]
    assert c_srcp > INTERIOR_BAND
    assert c_src > INTERIOR_BAND
    t_qsteady = time.perf_counter() - t00
    _p0(f"WF1S16 W1a leg {t_qsteady:.1f}s (record {REC['q_leg_s']} s, {Q1A}:71)")

    # (I=) collective; its own prints are held back (they carry the unrounded theta_max)
    with contextlib.redirect_stdout(io.StringIO()):
        g._check_summary(run, comm, "quad steady", transient=False)
        g._check_summary(run_u, comm, "uniform p2", transient=False)
    disk_qs = read_summary(run, comm)
    disk_up = read_summary(run_u, comm)
    tqs, tup = disk_qs["temperature"], disk_up["temperature"]
    p_json_qs = abs(tqs["source_power_w"] / q_v - 1.0)
    g_json_qs = abs(tqs["convected_power_w"] - tqs["source_power_w"]) / tqs["source_power_w"]
    say(f"[(S-json)] quadratic steady: source {tqs['source']!r}, thermal_degree {tqs['thermal_degree']}; "
        f"uniform P2: source {tup['source']!r}, thermal_degree {tup['thermal_degree']} (from disk)")
    say(f"[(P-json)] quadratic steady: source_power_w {tqs['source_power_w']:.4f} W vs Q*V {q_v:.4f} W: "
        f"rel {p_json_qs:.3e} (band {L0_BAND:.0e}; record +{REC['src_power_rel']:.3e}, {Q1A}:66; "
        "predicted ~1e-13)")
    say(f"[(G-json)] quadratic steady: |convected_power_w - source_power_w| / source_power_w = "
        f"{g_json_qs:.3e} (printed only; predicted ~1e-13, range [0, 1e-11]; within S2_GAUSS_BAND "
        f"{g_json_qs <= S2_GAUSS_BAND})")
    say(f"[theta] quadratic steady theta_max {fmt(tqs['theta_max_k'])} K, uniform P2 "
        f"{fmt(tup['theta_max_k'])} K; quadratic/uniform ratio {tqs['theta_max_k'] / tup['theta_max_k']:.4g} "
        "(no band)")
    assert tqs["source"] == "transverse_quadratic" and tqs["thermal_degree"] == DEGREE
    assert tup["source"] == "uniform" and tup["thermal_degree"] == DEGREE
    assert p_json_qs <= L0_BAND, f"(P-json) quadratic steady {p_json_qs:.3e}"

    # ---- quadratic steady fields to ParaView (P1, real parts) -------------------------------------
    f_qfile = p1_field(msh, "theta_quadratic_file", th.temperature, T_INF)
    f_qser = p1_field(msh, "theta_quadratic_series", ex_q)
    f_qdif = diff_field(msh, "theta_quadratic_minus_series", f_qfile, f_qser)
    f_ufile = p1_field(msh, "theta_uniform_p2_file", T_u, T_INF)
    b2c2 = HALF[1] ** 2 + HALF[2] ** 2
    f_qsrc = p1_field(msh, "q_transverse_quadratic",
                      lambda x: 3.0 * Q * (x[1] ** 2 + x[2] ** 2) / b2c2)
    ft = th.facet_tags
    ft.name = "robin_facets"
    qs_fields = {f.name: f for f in (f_qfile, f_qser, f_qdif, f_ufile, f_qsrc)}
    xdmf_qs, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_quadratic_steady_combined", msh, cell_tags_q, qs_fields, comm=comm,
        facet_tags=ft)

    # ==== quadratic CN leg: the gate's W1b body (q:343-380) ====================================
    t01 = time.perf_counter()
    out = RUN_QCN
    _fresh_dir(out, comm)
    rc_c, cap_c = cli_run(q.QUAD_CN, out)
    t_file = time.perf_counter() - t01
    assert _rc_all(rc_c, comm) == (0, 0) and len(cap_c) == 1, (rc_c, len(cap_c))
    run_c = cap_c[0]
    th = run_c.thermal
    msh = th.mesh
    py = _quad_study().temperature(_cooling(),
                                   transient=TransientSpec(T_STAR / N_MAIN, N_MAIN, RHO_C))
    T_py = _on_p(msh, py.temperature)
    ip = _owned_rel(T_py, th.temperature, comm)
    dofs = [_ndofs(th.temperature), _ndofs(py.temperature)]
    mult = np.array_equal(np.asarray(th.source_multipliers), np.asarray(py.source_multipliers))
    th_t = th.temperature - T_INF
    chi = _interior_chi(msh)
    ex_t, npts = _exact_p(msh, _quad_trans_ex(T_STAR))
    t_int = _rel_p(msh, th_t, ex_t, chi)
    ex_ut, _ = _exact_p(msh, _trans_ex(T_STAR))
    c_src_t = _rel_p(msh, th_t, ex_ut, chi)
    bal = run_c.summary["temperature"]["max_balance_residual"]
    say()
    _p0(f"WF1S16 (Pdeg) CN P2 dofs file / Python {dofs} (P36 {P36_DOFS}; {Q1B}:63); multipliers equal "
        f"{mult}; file run {t_file:.1f}s (record {REC['cn_file_s']} s); summary source "
        f"{run_c.summary['temperature']['source']!r}")
    _p0(f"WF1S16 (I-Py) CN: owned max|dtheta|/max|theta| = {ip:.3e} (hex {float(ip).hex()}; band "
        f"{L0_BAND:.0e}; record {REC['i_py_cn']:.3e}, {Q1B}:64)")
    _p0(f"WF1S16 (T-quad) CN N={N_MAIN} at t*: interior vs quadratic series {t_int:.6e} (P36 "
        f"{P36_T}; band {INTERIOR_BAND:.0e}); {npts} P3 points")
    _p0(f"WF1S16 (E-quad) CN max balance residual {bal:.3e} (P36 {P36_E:.3e}; record "
        f"{REC['e_quad']:.3e}, {Q1B}:66; band {E_BAND:.0e})")
    _p0(f"WF1S16 control (c-src-T) CN quadratic vs uniform transient series {c_src_t:.4e} "
        f"(P36 {REC['c_src_t']:.4e}; printed only)")
    assert dofs == [P36_DOFS] * 2
    assert mult
    assert run_c.summary["temperature"]["source"] == "transverse_quadratic"
    assert run_c.summary["temperature"]["duty"] is None
    assert ip <= L0_BAND, f"(I-Py) CN {ip:.3e}"
    assert t_int <= INTERIOR_BAND
    assert f"{t_int:.6e}" == P36_T
    assert bal <= E_BAND
    t_qcn = time.perf_counter() - t01
    _p0(f"WF1S16 W1b leg {t_qcn:.1f}s (record {REC['cn_leg_s']} s, {Q1B}:68)")

    with contextlib.redirect_stdout(io.StringIO()):
        g._check_summary(run_c, comm, "quad cn", transient=True)
    disk_qc = read_summary(run_c, comm)
    tqc = disk_qc["temperature"]
    p_json_qc = abs(tqc["source_power_w"] / q_v - 1.0)
    mults_json = tqc["source_multipliers"]
    mults_mem = [float(m) for m in th.source_multipliers]
    say(f"[(S-json)] quadratic CN: source {tqc['source']!r}, thermal_degree {tqc['thermal_degree']} "
        "(from disk)")
    say(f"[(P-json)] quadratic CN: source_power_w {tqc['source_power_w']:.4f} W vs Q*V {q_v:.4f} W: rel "
        f"{p_json_qc:.3e} (printed only; predicted ~1e-13, range [0, 1e-11])")
    say(f"[(E-json)] quadratic CN: max_balance_residual {tqc['max_balance_residual']:.3e} (band "
        f"{E_BAND:.0e}); source_multipliers {len(mults_json)} entries (a constant-source CN records "
        f"none; the 41 are the duty file's) == in-memory {mults_json == mults_mem}")
    say(f"[theta] quadratic CN theta_max at t* {fmt(tqc['theta_max_k'])} K (no band)")
    assert tqc["source"] == "transverse_quadratic" and tqc["thermal_degree"] == DEGREE
    assert tqc["max_balance_residual"] <= E_BAND
    assert mults_json == mults_mem

    f_cn = p1_field(msh, "theta_quadratic_cn_tstar", th.temperature, T_INF)
    f_cns = p1_field(msh, "theta_quadratic_cn_series_tstar", ex_t)
    f_cnd = diff_field(msh, "theta_quadratic_cn_minus_series", f_cn, f_cns)
    f_cnu = p1_field(msh, "theta_uniform_cn_series_tstar", ex_ut)
    cn_fields = {f.name: f for f in (f_cn, f_cns, f_cnd, f_cnu)}
    xdmf_qc, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_quadratic_cn_combined", msh, ones_tags(msh), cn_fields, comm=comm)

    # ==== probe steady leg: the gate's W1 body (p:258-312) =====================================
    t02 = time.perf_counter()
    out = RUN_PS
    _fresh_dir(out, comm)
    rc_p, cap_p = cli_run(p.STEADY, out)
    t_file = time.perf_counter() - t02
    assert _rc_all(rc_p, comm) == (0, 0) and len(cap_p) == 1, (rc_p, len(cap_p))
    run_ps = cap_p[0]
    pr = run_ps.probes
    py = _py_study().probe([WHOLE_BOX, *BOX_PROBES], RobinCooling(K_GEL, H_GEL, T_INF))
    m_f = pr.theta.function_space.mesh
    m_p = py.theta.function_space.mesh
    ca, cb, geo = _same_mesh(m_f, m_p, comm)
    # (I-Py) on the temperatures (T - T_INF normalises; θ itself is the rise).
    ip = _theta_owned_rel(_on(m_p, pr.thermal.temperature), py.thermal.temperature, comm)
    # run.thermal (the temperature block) vs the probe's re-solve
    two = _theta_owned_rel(_on(run_ps.thermal.temperature.function_space.mesh,
                               pr.thermal.temperature), run_ps.thermal.temperature, comm)
    vals_f = [float(v) for v in pr.values]
    vals_p = [float(v) for v in py.values]
    say()
    _p0(f"WF1S17 (UT) steady: cells {ca} / {cb}; geometry equal {geo}; file run {t_file:.1f}s (record "
        f"{REC['pr_file_s']} s, {P1L}:61, real build; cells record 79 872)")
    _p0(f"WF1S17 (I-Py) steady: owned max|dT|/max|T-Tinf| = {ip:.3e} (band {L0_BAND:.0e}; record "
        f"{REC['i_py_pr']:.3e}, {P1L}:62 real, {E130}:63 complex); temperature block vs probe re-solve "
        f"{two:.3e}")
    _p0("WF1S17 probe values file " + ", ".join(f"{v:.6g}" for v in vals_f))
    _p0("WF1S17 probe values py   " + ", ".join(f"{v:.6g}" for v in vals_p) + " (printed only)")
    mean_h, v_h = _domain_mean(pr.theta, comm)
    a_val = abs(vals_f[0] - mean_h) / abs(mean_h)
    a_vol = abs(float(pr.volumes[0]) / v_h - 1.0)
    _p0(f"WF1S17 (A-box) |P0 - mean| / |mean| = {a_val:.3e}, |V0/V_h - 1| = {a_vol:.3e} (band "
        f"{B_BAND:.0e}; U60 {U60_A_BOX:.3e}; record {REC['a_box'][0]:.3e} / {REC['a_box'][1]:.3e}, "
        f"{P1L}:65, {E130}:64)")
    mv = Q * A_BOXPROBE**2 / (10.0 * K_GEL)
    ratios, digits = [], []
    for j, c in enumerate(BOX_CENTRES):
        cf = _series_at(c) - mv
        r = abs(vals_f[1 + j] / cf - 1.0)
        ratios.append(r)
        digits.append(f"{r:.4e}")
        _p0(f"WF1S17 (R-box) c{j} = {c}: P = {vals_f[1 + j]:.5e}, ball mean = {cf:.5e}, rel = "
            f"{r:.4e} (U60 {U60_R_BOX[j]}; E130 {E130_R_BOX[j]}; band {INTERIOR_BAND:.0e})")
    # (c-wall) through the tracked file
    out_w = RUN_PW
    _fresh_dir(out_w, comm)
    rc_pw, cap_pw = cli_run(p.WRONG_WALL, out_w)
    assert _rc_all(rc_pw, comm) == (0, 0) and len(cap_pw) == 1, (rc_pw, len(cap_pw))
    run_w = cap_pw[0]
    c_wall = abs(float(run_w.probes.values[0]) / (_series_at(BOX_CENTRES[0]) - mv) - 1.0)
    _p0(f"WF1S17 control (c-wall) wrong-wall file at c0: miss {c_wall:.4e} (U60 "
        f"{U60_C_WALL:.4e}; equal digits {f'{c_wall:.4e}' == f'{U60_C_WALL:.4e}'}; asserted > "
        f"{INTERIOR_BAND:.0e})")
    with contextlib.redirect_stdout(io.StringIO()):
        p._check_summary(run_ps, comm, "steady")
        p._check_summary(run_w, comm, "wrong wall")
    say(f"  (U') steady and wrong wall: values == run.probes.values and JSON == summary, asserted inside "
        f"the gate helper ({P1L}:72-73); its own line held back")
    assert ca == cb and geo
    assert ip <= L0_BAND, f"(I-Py) {ip:.3e}"
    assert two <= L0_BAND, f"temperature block vs probe {two:.3e}"
    assert a_val <= B_BAND and a_vol <= B_BAND
    for r in ratios:
        assert r <= INTERIOR_BAND
    # substitution 3: the real-build digit assertion is a printed equality in the complex build
    say(f"[(R-box) digits] complex-build digits {tuple(digits)}; equal U60's {tuple(digits) == U60_R_BOX}; "
        f"equal {E130}:65-69's {tuple(digits) == E130_R_BOX} (printed, not asserted)")
    assert c_wall > INTERIOR_BAND
    t_psteady = time.perf_counter() - t02
    _p0(f"WF1S17 W1 steady leg {t_psteady:.1f}s (record {REC['pr_leg_s']} s, real build, {P1L}:74)")

    # ---- the probe table, read from the JSON on disk ------------------------------------------------
    disk_ps = read_summary(run_ps, comm)
    disk_pw = read_summary(run_w, comm)
    js_pr = disk_ps["probes"]
    tps = disk_ps["temperature"]
    jv = js_pr["values_k"]
    # (R-json)
    rj = []
    for j, c in enumerate(BOX_CENTRES):
        cf = _series_at(c) - mv
        rj.append(abs(jv[1 + j] / cf - 1.0))
    a_json = abs(jv[0] - mean_h) / abs(mean_h)
    a_json_vol = abs(js_pr["volumes_m3"][0] / v_h - 1.0)
    rj_digits = [f"{r:.4e}" for r in rj]
    say()
    say(f"[(S-json)] probes steady: source {tps['source']!r}, 'probes' present {'probes' in disk_ps}; "
        f"wrong wall: source {disk_pw['temperature']['source']!r}, 'probes' present "
        f"{'probes' in disk_pw}")
    say(f"[(R-json)] (A-box) from the JSON values_k[0]: {a_json:.3e}, volumes_m3[0]/V_h - 1 = "
        f"{a_json_vol:.3e} (band {B_BAND:.0e})")
    for j, c in enumerate(BOX_CENTRES):
        say(f"[(R-json)] c{j} {c}: JSON ball value {jv[1 + j]:.5g} K vs series-minus-correction ratio "
            f"{rj[j]:.4e} (band {INTERIOR_BAND:.0e}; U60 {U60_R_BOX[j]}; E130 {E130_R_BOX[j]})")
    say(f"[(R-json) digits] from disk {tuple(rj_digits)}; equal U60's {tuple(rj_digits) == U60_R_BOX} "
        "(printed, not asserted)")
    assert tps["source"] == "uniform" and "probes" in disk_ps
    assert disk_pw["temperature"]["source"] == "uniform" and "probes" in disk_pw
    assert a_json <= B_BAND and a_json_vol <= B_BAND
    for r in rj:
        assert r <= INTERIOR_BAND

    # (C-json) the centres / radii written back equal the parse node's hex lists
    def json_hex(js):
        return [([float(x).hex() for x in c], float(r).hex())
                for c, r in zip(js["centres"], js["radii"])]

    c_steady = json_hex(js_pr) == p._probe_hex([WHOLE_BOX, *BOX_PROBES])
    c_wall_j = json_hex(disk_pw["probes"]) == p._probe_hex(BOX_PROBES[:1])
    say(f"[(C-json)] probes.centres / probes.radii by float.hex: steady == [WHOLE_BOX, *BOX_PROBES] "
        f"{c_steady}; wrong wall == BOX_PROBES[:1] {c_wall_j}")
    assert c_steady and c_wall_j

    # ==== probe duty leg: the gate's duty body (p:317-345) ======================================
    t03 = time.perf_counter()
    out = RUN_PD
    _fresh_dir(out, comm)
    rc_d, cap_d = cli_run(p.DUTY, out)
    assert _rc_all(rc_d, comm) == (0, 0) and len(cap_d) == 1, (rc_d, len(cap_d))
    run_d = cap_d[0]
    prd = run_d.probes
    pyd = _py_study().probe(
        BOX_PROBES, RobinCooling(K_GEL, H_GEL, T_INF),
        transient=TransientSpec(DT_BOX, N_STEPS, RHO_C), duty=DutyCycle(2.0 * HALF_T, 0.5),
    )
    m_fd = prd.theta.function_space.mesh
    m_pd2 = pyd.theta.function_space.mesh
    cad, cbd, geod = _same_mesh(m_fd, m_pd2, comm)
    ipd = _theta_owned_rel(_on(m_pd2, prd.thermal.temperature), pyd.thermal.temperature, comm)
    multd = np.array_equal(np.asarray(prd.thermal.source_multipliers),
                           np.asarray(pyd.thermal.source_multipliers))
    mult_s = np.array_equal(np.asarray(run_d.summary["temperature"]["source_multipliers"]),
                            np.asarray(pyd.thermal.source_multipliers))
    say()
    _p0(f"WF1S17 (I-tr) duty: cells {cad} / {cbd}; geometry equal {geod}; owned max|dT|/max|T-Tinf| "
        f"= {ipd:.3e} (band {L0_BAND:.0e}; record {REC['i_tr']:.3e}, {P1L}:77); multipliers equal "
        f"{multd} / summary {mult_s} ({len(pyd.thermal.source_multipliers)} entries)")
    _p0("WF1S17 (I-tr) probe values file " + ", ".join(f"{float(v):.6g}" for v in prd.values))
    _p0("WF1S17 (I-tr) probe values py   " + ", ".join(f"{float(v):.6g}" for v in pyd.values)
        + " (printed; no comparand is posed)")
    with contextlib.redirect_stdout(io.StringIO()):
        p._check_summary(run_d, comm, "duty")
    say(f"  (U') duty: values == run.probes.values and JSON == summary, asserted inside the gate helper "
        f"({P1L}:80); its own line held back")
    assert cad == cbd and geod
    assert ipd <= L0_BAND, f"(I-tr) {ipd:.3e}"
    assert multd and mult_s
    assert run_d.summary["temperature"]["duty"] == {"period_s": 2.0 * HALF_T, "duty": 0.5}
    t_pduty = time.perf_counter() - t03
    _p0(f"WF1S17 W1 duty leg {t_pduty:.1f}s (record {REC['duty_leg_s']} s, real build, {P1L}:81)")

    disk_pd = read_summary(run_d, comm)
    c_duty = json_hex(disk_pd["probes"]) == p._probe_hex(BOX_PROBES)
    say(f"[(C-json)] duty: probes.centres / probes.radii by float.hex == BOX_PROBES {c_duty}; "
        f"duty source {disk_pd['temperature']['source']!r}")
    say("[duty values_k from disk] " + ", ".join(f"{v:.6g}" for v in disk_pd["probes"]["values_k"])
        + " K (printed only; no comparand posed; record row in " + f"{P1L}:78)")
    assert c_duty
    assert disk_pd["temperature"]["source"] == "uniform"

    # ---- probe fields to ParaView (P1 on the probe mesh) ------------------------------------------
    f_pfile = p1_field(m_f, "theta_probe_file_steady", pr.theta)
    f_pw = p._on(m_f, run_w.probes.theta)
    f_pw.name = "theta_probe_wrong_wall"
    cen = np.array(BOX_CENTRES, dtype=float)

    def balls(x):
        d = np.sqrt((x[0][:, None] - cen[None, :, 0]) ** 2 + (x[1][:, None] - cen[None, :, 1]) ** 2
                    + (x[2][:, None] - cen[None, :, 2]) ** 2)
        return (d.min(axis=1) <= A_BOXPROBE).astype(float)

    f_balls = p1_field(m_f, "probe_balls", balls)
    pr_fields = {f.name: f for f in (f_pfile, f_pw, f_balls)}
    xdmf_pr, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_probes_combined", m_f, ones_tags(m_f), pr_fields, comm=comm)
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- CSV: the summary read-back, rank 0 ------------------------------------------------------
    csv_path = OUTPUT_DIR / f"{BASENAME}_summary_readback.csv"
    rows = []
    for label, disk, mem in (("quadratic_steady", disk_qs, run.summary),
                             ("uniform_p2", disk_up, run_u.summary),
                             ("quadratic_cn", disk_qc, run_c.summary),
                             ("probes_steady", disk_ps, run_ps.summary),
                             ("probes_wrong_wall", disk_pw, run_w.summary),
                             ("probes_duty", disk_pd, run_d.summary)):
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

    # ---- CSV: the probe table, rank 0 -----------------------------------------------------------
    pt_path = OUTPUT_DIR / f"{BASENAME}_probe_table.csv"
    prow = []
    for i, (c, rad, val, vol) in enumerate(zip(js_pr["centres"], js_pr["radii"], jv,
                                               js_pr["volumes_m3"])):
        if i == 0:
            comp, ratio, u60, e130, eq = mean_h, a_json, "", "", ""
            comp_s = f"{comp:.6g}"
        else:
            comp = _series_at(BOX_CENTRES[i - 1]) - mv
            ratio = rj[i - 1]
            u60, e130 = U60_R_BOX[i - 1], E130_R_BOX[i - 1]
            eq = f"{ratio:.4e}" == u60
            comp_s = f"{comp:.6g}"
        prow.append((i, "(" + " ".join(f"{x:.6g}" for x in c) + ")", f"{rad:.6g}", f"{val:.6g}",
                     f"{vol:.6g}", comp_s, f"{ratio:.4e}", u60, e130, eq))
    if rank0:
        with open(pt_path, "w") as fh:
            fh.write("probe,centre_m,radius_m,value_k,volume_m3,comparand_k,ratio,u60_digit,"
                     "ex130_digit,digits_equal\n")
            for r in prow:
                fh.write(",".join(str(x) for x in r) + "\n")
    comm.Barrier()
    assert all(r[4] for r in rows), [r for r in rows if not r[4]]

    say()
    say(f"[output] {xdmf_qs.name if rank0 else ''}: P1 theta_quadratic_file, theta_quadratic_series, "
        "theta_quadratic_minus_series, theta_uniform_p2_file, q_transverse_quadratic, Robin facets")
    say(f"[output] {xdmf_qc.name if rank0 else ''}: P1 theta_quadratic_cn_tstar, "
        "theta_quadratic_cn_series_tstar, theta_quadratic_cn_minus_series, "
        "theta_uniform_cn_series_tstar")
    say(f"[output] {xdmf_pr.name if rank0 else ''}: P1 theta_probe_file_steady, theta_probe_wrong_wall, "
        "probe_balls")
    say(f"[output] {csv_path.name} ({len(rows)} rows, every JSON value hex-equal to the in-memory "
        "summary)")
    say(f"[output] {pt_path.name} ({len(prow)} rows)")
    say(f"[time] quadratic steady leg {t_qsteady:.1f} s; CN leg {t_qcn:.1f} s; probe steady leg "
        f"{t_psteady:.1f} s; probe duty leg {t_pduty:.1f} s; total {time.perf_counter() - started:.1f} s")
    say("[scope] six thermal_box files through cli.main plus one malformed: transverse-quadratic and "
        "uniform P2 steady, quadratic CN, ball probes steady / wrong wall / duty; m = 4, complex "
        "build, -n 2; summaries and the probe table read back from disk; no file carries both "
        "excitation.source and outputs (not gated); the probe leg is an example record in the complex "
        "build, not a gate; no EM, no coil, no compliance, no external check claimed")
    say("All assertions hold")


if __name__ == "__main__":
    main()
