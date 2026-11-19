"""`WF-15` step 4 (§9 item 326, chain U30): ``Study.probe(probes, cooling)``.

`MAT-19`'s ball-average virtual temperature probe through the API, on
`MAT-12`'s fixture (the same EM mesh, the 19 496-cell tag-1 submesh) under the
steady Robin temperature the Study computes (item 260), EM degree 1, P1
thermal.  Anchors (brief: docs/planning/chunks/WF-15.md, "§9 item 326"):

(A)   P8's quadrature identity through the API: the whole-sphere probe reads
      ∫θ dx / V_sub and V_sub within ``B_BAND`` (both widths);
(A0)  submesh global cells == 19 496 (both widths);
(U)   ``Study.probe(...).values`` == ``ball_average`` called directly on the
      Study's own θ (same process), and (Once): no second EM solve;
(W)   ``-n 2`` writes an npz, ``-n 1`` loads it: every value within
      ``WIDTH_BAND`` relative.
Controls: (c1) c₂ displaced 0.1 mm outward misses ``WIDTH_BAND`` (asserted);
(c2) a ghost-included sum misses (A)'s volume identity iff G_w > 0
(asserted, else "inapplicable"); (c3) T instead of θ (printed).

Every band, constant, helper and comparator is imported, none restated.
"""

from __future__ import annotations

import hashlib
import os
import time
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, RobinCooling, Study
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form
from fem_em_solver.thermal.probe import ProbeBall, _ball_average_on_cells, ball_average
from fem_em_solver.utils.constants import EPSILON_0
from tests.complex_mode import complex_only
from tests.package.test_wf15_study_lossy_sphere import _material_057, _study, _summed_rss_gib
from tests.package.test_wf15_study_temperature import MAT10_A_PRINT
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (
    B_BAND,
    T_INF,
    G,
    H,
    K,
    R,
    _allsum,
    _interior_indicator,
    _interior_rel,
    _say,
    _theta_cf,
)
from tests.validation.test_mat19_virtual_temperature_probe import (
    A_PROBE,
    CENTRES,
    N_SUB,
    WIDTH_BAND,
)

NPZ = "/workspace/logs/wf15/probe-n2.npz"
WHOLE = ProbeBall((0.0, 0.0, 0.0), 2.0 * R)  # every submesh vertex lies at r <= R
PROBES = [ProbeBall(c, A_PROBE) for c in CENTRES]
C2_OUT = ProbeBall(tuple((6.6 / 6.5) * np.asarray(CENTRES[2])), A_PROBE)  # (c1): 0.1 mm out


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _hashes():
    return _sha(study_mod.__file__), _sha(__file__)


def _q_cf():
    cf = uniform_sphere_sar_closed_form(
        e0=E0,
        epsilon_r=EPSILON_R_SPHERE,
        sigma=SIGMA_HIGH,
        omega=OMEGA,
        rho=RHO_KG_M3,
        epsilon_0=EPSILON_0,
    )
    return G * RHO_KG_M3 * cf["sar_w_per_kg"]


def _ss_ball_mean(d, a, q):
    """Steady closed-form ball mean, `MAT-19.md:69` at t → ∞ (printed only)."""
    return q * R / (3.0 * H) + q * (R * R - d * d - 3.0 * a * a / 5.0) / (6.0 * K)


def _probe_leg(comm):
    rd = {}
    study = _study(degree=1)
    t0 = time.perf_counter()
    res = study.solve(comm)
    rd["t_em"] = time.perf_counter() - t0
    spec = RobinCooling(K, H, T_INF, source_scale=G)
    t0 = time.perf_counter()
    pr = study.probe([WHOLE] + PROBES, spec)
    rd["t_probe"] = time.perf_counter() - t0
    rd["once"] = study.result is res
    sub = pr.thermal.submesh
    theta = pr.theta
    rd["cells"] = int(sub.topology.index_map(3).size_global)
    rd["values"] = np.asarray(pr.values, dtype=np.float64)
    rd["volumes"] = np.asarray(pr.volumes, dtype=np.float64)

    # (A) P8's identity through the API.
    v_sub = float(
        np.real(_allsum(comm, fem.Constant(sub, default_scalar_type(1.0)) * ufl.dx(domain=sub)))
    )
    int_th = float(np.real(_allsum(comm, theta * ufl.dx(domain=sub))))
    mean = int_th / v_sub
    rd["A_val"] = abs(pr.values[0] - mean) / abs(mean)
    rd["A_vol"] = abs(pr.volumes[0] / v_sub - 1.0)
    _say(
        f"WF15P (A0) sub cells {rd['cells']} (record {N_SUB}); (A) |P - int/V|/|int/V| = "
        f"{rd['A_val']:.3e}, |V_h/V_sub - 1| = {rd['A_vol']:.3e} (band {B_BAND:.0e}); "
        f"V_sub {v_sub:.10e}"
    )

    # (U) the API is MAT-19's function, same process.
    direct = np.array([ball_average(theta, p)[0] for p in [WHOLE] + PROBES])
    rd["U_equal"] = bool(np.array_equal(direct, rd["values"]))
    _say(
        f"WF15P (U) Study.probe values {rd['values'].tolist()!r}; direct {direct.tolist()!r}; "
        f"equal {rd['U_equal']}; (Once) {rd['once']}"
    )

    # Item 260's (A) on the same thermal field, imported comparator.
    q_cf = _q_cf()
    th_cf = _theta_cf(sub, q_cf)
    rd["err_A"] = _interior_rel(
        sub, _interior_indicator(sub), (pr.thermal.temperature - T_INF) - th_cf, th_cf
    )
    _say(f"WF15P item-260 (A) = {rd['err_A']:.10e} [{rd['err_A']:.4e}] vs {MAT10_A_PRINT}")

    # Printed: closed-form steady ball means, max_imag_ratio, (c3).
    for j, p in enumerate(PROBES):
        d = float(np.linalg.norm(p.centre))
        cf = _ss_ball_mean(d, A_PROBE, q_cf)
        _say(
            f"WF15P probe c{j} d={d * 1e3:.2f} mm: FEM {rd['values'][j + 1]:.6f} K, closed form "
            f"{cf:.6f} K, offset {rd['values'][j + 1] / cf - 1.0:+.4%}; volume "
            f"{rd['volumes'][j + 1]:.6e} / nominal {pr.nominal_volumes[j + 1]:.6e}"
        )
    T_whole = ball_average(pr.thermal.temperature, WHOLE)[0]
    _say(
        f"WF15P (c3) whole-ball T / theta = {T_whole / rd['values'][0]:.1f} (printed); "
        f"max_imag_ratio {pr.max_imag_ratio:.3e}"
    )
    rd["study"] = study
    rd["pr"] = pr
    return rd


def test_wf15_probe_spec_and_refusals():
    """Port Study ⇒ TypeError with no build; empty / non-ProbeBall ⇒ refusal."""
    spec = RobinCooling(K, H, T_INF)
    port = Study(
        FourPortBirdcageGeometry.f_small(), _material_057(), LumpedPortDrive(), frequency_hz=64e6
    )
    with pytest.raises(TypeError):
        port.probe(PROBES, spec)
    assert port.result is None
    sphere = _study(degree=1)
    with pytest.raises(ValueError):
        sphere.probe([], spec)
    with pytest.raises(TypeError):
        sphere.probe([((0.0, 0.0, 0.0), A_PROBE)], spec)
    assert sphere.result is None  # refusals precede any solve
    _say("WF15P (K) port Study TypeError, result None; empty ValueError; tuple TypeError")


@complex_only
def test_wf15_probe_identity_n2():
    """W1 (-n 2): (A), (A0), (U), (Once), item-260 (A), (c1), (c2); writes (W)'s npz."""
    comm = MPI.COMM_WORLD
    rd = _probe_leg(comm)
    pr = rd["pr"]
    theta = pr.theta

    # (c1) displaced probe.
    v_out = ball_average(theta, C2_OUT)[0]
    rd["c1"] = abs(v_out - rd["values"][3]) / abs(rd["values"][3])
    _say(
        f"WF15P (c1) c2 +0.1 mm: |d|/theta = {rd['c1']:.4e} (predicted ≈ 7.4e-3; band "
        f"{WIDTH_BAND:.0e})"
    )

    # (c2) ghost-included sum on the whole ball.
    sub = pr.thermal.submesh
    imap = sub.topology.index_map(3)
    g_w = int(comm.allreduce(int(imap.num_ghosts), op=MPI.SUM))
    cells_all = np.arange(imap.size_local + imap.num_ghosts, dtype=np.int32)
    _, vol_g = _ball_average_on_cells(theta, WHOLE, cells_all)
    c2 = abs(vol_g / rd["volumes"][0] - 1.0)
    if g_w > 0:
        _say(
            f"WF15P (c2) G_w = {g_w}; ghost-included volume miss {c2:.3e} (predicted ≈ "
            f"{g_w / N_SUB:.1e}; band {B_BAND:.0e})"
        )
    else:
        _say(f"WF15P (c2) G_w = 0: inapplicable (record); miss {c2:.3e}")
    _say(
        f"WF15P EM solve + mesh {rd['t_em']:.1f} s; temperature + probe {rd['t_probe']:.1f} s; "
        f"summed ru_maxrss {_summed_rss_gib(comm):.2f} GiB; width {comm.size}"
    )

    if comm.rank == 0:
        os.makedirs(os.path.dirname(NPZ), exist_ok=True)
        h_src, h_test = _hashes()
        np.savez(
            NPZ,
            cells=rd["cells"],
            values=rd["values"],
            volumes=rd["volumes"],
            sha_study=h_src,
            sha_test=h_test,
        )
    comm.Barrier()

    assert rd["once"], "(Once) a second EM solve happened"
    assert rd["cells"] == N_SUB, f"(A0) {rd['cells']}"
    assert rd["A_val"] <= B_BAND, f"(A) value {rd['A_val']:.3e}"
    assert rd["A_vol"] <= B_BAND, f"(A) volume {rd['A_vol']:.3e}"
    assert rd["U_equal"], "(U) Study.probe != ball_average"
    assert f"{rd['err_A']:.4e}" == MAT10_A_PRINT, f"item-260 (A) {rd['err_A']:.4e}"
    assert rd["c1"] > WIDTH_BAND, f"(c1) {rd['c1']:.3e} did not miss {WIDTH_BAND:.0e}"
    if g_w > 0:
        assert c2 > B_BAND, f"(c2) ghost-included sum {c2:.3e} did not miss {B_BAND:.0e}"


@complex_only
def test_wf15_probe_width_n1():
    """W2 (-n 1): (A), (A0), (U), (Once) and (W) against W1's npz."""
    comm = MPI.COMM_WORLD
    if not os.path.exists(NPZ):
        pytest.fail(f"(W) {NPZ} missing: run test_wf15_probe_identity_n2 at -n 2 first")
    ref = np.load(NPZ)
    h_src, h_test = _hashes()
    if str(ref["sha_study"]) != h_src or str(ref["sha_test"]) != h_test:
        pytest.fail("(W) study.py or the test module changed since the -n 2 npz was written")
    rd = _probe_leg(comm)
    rel = np.abs(rd["values"] - ref["values"]) / np.abs(ref["values"])
    _say(
        f"WF15P (W) width {comm.size} vs 2: cells {rd['cells']} / {int(ref['cells'])}; "
        f"rel |d| per probe {rel.tolist()!r}; max {float(np.max(rel)):.3e} (band "
        f"{WIDTH_BAND:.0e})"
    )
    _say(
        f"WF15P EM solve + mesh {rd['t_em']:.1f} s; temperature + probe {rd['t_probe']:.1f} s; "
        f"summed ru_maxrss {_summed_rss_gib(comm):.2f} GiB; width {comm.size}"
    )
    assert rd["once"], "(Once) a second EM solve happened"
    assert rd["cells"] == N_SUB == int(ref["cells"]), f"(A0) {rd['cells']}"
    assert rd["A_val"] <= B_BAND, f"(A) value {rd['A_val']:.3e}"
    assert rd["A_vol"] <= B_BAND, f"(A) volume {rd['A_vol']:.3e}"
    assert rd["U_equal"], "(U) Study.probe != ball_average"
    assert f"{rd['err_A']:.4e}" == MAT10_A_PRINT, f"item-260 (A) {rd['err_A']:.4e}"
    assert float(np.max(rel)) <= WIDTH_BAND, f"(W) {float(np.max(rel)):.3e}"
