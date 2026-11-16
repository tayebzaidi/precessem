"""Example (`EX-97`): a temperature from a config file -- ``femem run`` on a YAML with
a ``thermal:`` section reproduces `MAT-10`'s records (`WF-1` step 3).

`WF-1`'s ``fem_em_solver.cli`` reads a study from a YAML file and runs it; since step 3
the file may carry a ``thermal:`` section (k, h, T_inf, optional source_scale / degree)
that drives ``Study.temperature(RobinCooling(...))``.  This example drives
``cli.main(["run", <yaml>, "--output-dir", <dir>])`` in-process on the sphere file
beside it (`MAT-4` step 1's lossy sphere plus `MAT-10`'s Robin cooling, degree 1) and
reads the ``temperature`` block of the JSON summary.  It then calls the exported
``run_study_file`` on the same file (the function ``femem run`` calls), which returns
the temperature field ``cli.main`` discards, and reads `MAT-10`'s (A) / (P) / (B) / (C)
off that field term for term as ``test_wf1_thermal_yaml_reproduces_mat10_record``.

Every constant, band and helper is imported from the gate modules (the ANS-1 rule);
k, h, T_inf and the source scale come from the YAML or the imports.  The script has
no numeric literal of its own.

Asserted:

* (K)      each YAML copy is spec- and cooling-equal (``float.hex``) to its
           ``tests/package/wf1/`` original, and the sphere copy to the Python-built
           ``Study`` and ``RobinCooling(K, H, T_INF, source_scale=G)``;
* (R-cli)  ``run`` returns 0 on every rank; the JSON summary names the family, degree
           and ranks, and its three ``temperature`` values are within ``L0_BAND`` of
           ``run_study_file``'s (two EM solves of one file);
* (R)      one EM solve; the sphere submesh has `MAT-10`'s cell count; (A) and (P) at
           ``.4e`` equal `MAT-10`'s prints; (A) / (P) / (B) / (C) inside their bands;
* (IT)     the YAML temperature equals an independently solved Python ``Study``'s
           within ``L0_BAND`` (owned dofs).

Negative controls:

* (c1)  asserted: the control YAML (``source_scale`` 2G, the missing 1/2) applied to
        the already-solved Study -- the (A) comparator exceeds ``A_BAND`` (its ``.4e``
        string is printed beside `MAT-10`'s record, not asserted);
* (c1') printed: the source-power ratio of the control to the file;
* (X)   asserted: ``thermal:`` on a port geometry raises ``StudyConfigError`` naming
        ``thermal``; ``run`` and ``check`` return 2 on every rank; no output dir.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:13 -n 2 -t 180

Not claimed: degree 2 through YAML; transient heating, Pennes or duty cycle; thermal on
a port Study; tissue or a coil; a stable schema; the installed-wheel ``femem`` with a
solve or ``femem run --ranks``; anything `WF-1` / `WF-15` / `MAT-10` do not claim.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver import (  # noqa: E402
    RobinCooling,
    StudyConfigError,
    cli,
    load_study,
    run_study_file,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.sar import uniform_sphere_sar_closed_form  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402
from fem_em_solver.thermal.steady import solve_steady_conduction  # noqa: E402
from fem_em_solver.utils.constants import EPSILON_0  # noqa: E402

from tests.package.test_wf15_study_lossy_sphere import _study  # noqa: E402
from tests.package.test_wf15_study_temperature import (  # noqa: E402
    MAT10_A_PRINT,
    MAT10_B_MISSING_HALF_PRINT,
    MAT10_P_PRINT,
    MAT10_SUB_CELLS,
    _theta_owned_rel,
)
from tests.package.test_wf1_study_yaml import (  # noqa: E402
    WF1,
    _fresh_dir,
    _rc_all,
    _spec_diff,
)
from tests.package.test_wf1_yaml_thermal import (  # noqa: E402
    CONTROL_2G_FILE,
    MALFORMED_THERMAL,
    THERMAL_FILE,
    _cooling,
    _hex_diff,
)
from tests.validation.test_lossy_sphere_sar import (  # noqa: E402
    E0,
    EPSILON_R_SPHERE,
    OMEGA,
    RHO_KG_M3,
)
from tests.validation.test_mat10_sar_to_heat_sphere import (  # noqa: E402
    A_BAND,
    B_BAND,
    C_BAND,
    G,
    H,
    K,
    P_BAND,
    ROBIN_TAG,
    T_INF,
    _allsum,
    _interior_indicator,
    _interior_rel,
    _theta_cf,
)
from tests.validation.test_lossy_sphere_sar import SIGMA_HIGH  # noqa: E402
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND  # noqa: E402

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "paraview_output"
FIGURE_DIR = HERE / "figures"
BASENAME = "mri_13_temperature_from_yaml"
EX_YAML = HERE / "13_temperature_from_yaml_sphere.yaml"
EX_CONTROL = HERE / "13_temperature_from_yaml_sphere_control_2g.yaml"
MALFORMED_NAME = "malformed_thermal_on_port.yaml"
TEMPERATURE_KEYS = {"source_power_w", "convected_power_w", "theta_max_k"}
AIR_TAG = 2


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    OUTPUT_DIR.mkdir(exist_ok=True)
    comm.Barrier()

    # ---- (K): the copies equal the originals, no solve ----------------------
    ex_study = load_study(EX_YAML)
    ref_cooling = RobinCooling(K, H, T_INF, source_scale=G)
    k_spec_orig = _spec_diff(ex_study, load_study(THERMAL_FILE))
    k_cool_orig = _hex_diff(_cooling(EX_YAML), _cooling(THERMAL_FILE))
    k_spec_ctl = _spec_diff(load_study(EX_CONTROL), load_study(CONTROL_2G_FILE))
    k_cool_ctl = _hex_diff(_cooling(EX_CONTROL), _cooling(CONTROL_2G_FILE))
    k_spec_py = _spec_diff(ex_study, _study(degree=1))
    k_cool_py = _hex_diff(_cooling(EX_YAML), ref_cooling)

    # ---- (X): thermal on a port geometry, no solve --------------------------
    bad = WF1 / MALFORMED_NAME
    expected_field = MALFORMED_THERMAL[MALFORMED_NAME]
    try:
        load_study(bad)
        x_message, x_field = None, None
    except StudyConfigError as exc:
        x_message, x_field = str(exc), exc.field
    bad_out = OUTPUT_DIR / f"{BASENAME}_malformed"
    _fresh_dir(bad_out, comm)
    rc_run = _rc_all(cli.main(["run", str(bad), "--output-dir", str(bad_out)]), comm)
    rc_check = _rc_all(cli.main(["check", str(bad)]), comm)
    bad_exists = comm.allreduce(bool(bad_out.exists()), op=MPI.LOR)

    # ---- (R-cli): femem run on the sphere YAML ------------------------------
    cli_out = OUTPUT_DIR / f"{BASENAME}_cli"
    _fresh_dir(cli_out, comm)
    t0 = time.perf_counter()
    rc_cli = _rc_all(cli.main(["run", str(EX_YAML), "--output-dir", str(cli_out)]), comm)
    t_cli = time.perf_counter() - t0
    comm.Barrier()
    with open(cli_out / "summary.json") as fh:
        cli_summary = json.load(fh)

    # ---- the API route: run_study_file returns the field cli.main discards --
    api_out = OUTPUT_DIR / f"{BASENAME}_api"
    _fresh_dir(api_out, comm)
    t0 = time.perf_counter()
    run = run_study_file(EX_YAML, output_dir=api_out, comm=comm)
    t_api = time.perf_counter() - t0
    study, res, th = run.study, run.result, run.thermal
    assert th is not None
    once_api = study.result is res
    api_temp = run.summary["temperature"]
    cli_temp = cli_summary["temperature"]
    cli_keys = set(cli_temp)
    cli_rel = {
        key: abs(cli_temp[key] - api_temp[key]) / abs(api_temp[key]) for key in sorted(api_temp)
    }

    cooling = th.cooling
    write_setup_figure(
        res.mesh,
        res.cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            study.geometry.sphere_tag: (
                f"gel sphere (Robin-cooled, h {cooling.h_w_per_m2_k:g} W/m²K)"
            ),
            AIR_TAG: "air box (imposed-field walls)",
        },
        hide_tags=(AIR_TAG,),
        translucent_tags=(study.geometry.sphere_tag,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="mri:13 -- the sphere Study with thermal: loaded from YAML (WF-1 step 3)",
        comm=comm,
    )

    # ---- (R): MAT-10's records on the YAML theta, term for term -------------
    cf = uniform_sphere_sar_closed_form(
        e0=E0, epsilon_r=EPSILON_R_SPHERE, sigma=SIGMA_HIGH, omega=OMEGA,
        rho=RHO_KG_M3, epsilon_0=EPSILON_0,
    )
    q_cf = G * RHO_KG_M3 * cf["sar_w_per_kg"]
    sub, ft = th.submesh, th.facet_tags
    n_sub = int(sub.topology.index_map(sub.topology.dim).size_global)
    chi = _interior_indicator(sub)
    theta_cf = _theta_cf(sub, q_cf)
    err_a = _interior_rel(sub, chi, (th.temperature - T_INF) - theta_cf, theta_cf)
    res_p = solve_steady_conduction(
        sub, ft, k=K, q=q_cf, h=H, t_inf=T_INF, robin_tag=ROBIN_TAG, degree=1
    )
    err_p = _interior_rel(sub, chi, (res_p.temperature - T_INF) - theta_cf, theta_cf)
    int_q = float(np.real(_allsum(comm, th.heat_source * ufl.dx(domain=sub))))
    p_parent = cooling.source_scale * study.mean_sar()["dissipated_power_w"]
    err_b = abs(int_q - p_parent) / abs(p_parent)
    err_c = abs(th.convected_power - th.source_power) / abs(th.source_power)

    # ---- (c1): the control file's cooling on the already-solved Study -------
    th2 = study.temperature(_cooling(EX_CONTROL))
    once_c1 = study.result is res
    sub2 = th2.submesh
    theta2 = _theta_cf(sub2, q_cf)
    c1 = _interior_rel(sub2, _interior_indicator(sub2), (th2.temperature - T_INF) - theta2, theta2)
    c1_ratio = float(np.real(th2.source_power)) / float(np.real(th.source_power))

    # ---- (IT): an independently solved Python Study -------------------------
    t0 = time.perf_counter()
    study_py = _study(degree=1)
    study_py.solve(comm)
    th_py = study_py.temperature(ref_cooling)
    t_py = time.perf_counter() - t0
    it, th_max = _theta_owned_rel(th.temperature, th_py.temperature, comm)

    # ---- output fields on the sphere submesh (collective on every rank) -----
    theta_out = fem.Function(th.temperature.function_space, name="theta")
    theta_out.x.array[:] = np.real(th.temperature.x.array) - cooling.t_inf_c
    theta_out.x.scatter_forward()
    q_out = fem.Function(th.heat_source.function_space, name="q")
    q_out.x.array[:] = np.real(th.heat_source.x.array)
    q_out.x.scatter_forward()
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        sub,
        None,
        {"theta": theta_out, "q": q_out},
        comm=comm,
        facet_tags=ft,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- prints (rank 0, no collective inside) ------------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-97 -- a temperature from a config file (WF-1 step 3)")
        print("=" * 72)
        print(
            f"[(K)] sphere vs original spec {k_spec_orig} cooling {k_cool_orig}; control vs "
            f"original spec {k_spec_ctl} cooling {k_cool_ctl}; sphere vs Python spec "
            f"{k_spec_py} cooling {k_cool_py}"
        )
        print(
            f"[(R-cli)] run rc (min, max) {rc_cli}, {t_cli:.1f} s; family "
            f"{cli_summary['family']!r}, degree {cli_summary['degree']}, ranks "
            f"{cli_summary['ranks']}; keys {sorted(cli_keys)}"
        )
        for key in sorted(api_temp):
            print(
                f"[(R-cli)] {key}: cli {cli_temp[key]!r}; api {api_temp[key]!r}; "
                f"rel {cli_rel[key]:.3e} (band L0_BAND {L0_BAND:g})"
            )
        print(
            f"[(R)] run_study_file {t_api:.1f} s; one EM solve {once_api}; sub cells {n_sub} "
            f"(record {MAT10_SUB_CELLS}); (A) {err_a:.4e} vs {MAT10_A_PRINT}; (P) {err_p:.4e} "
            f"vs {MAT10_P_PRINT}"
        )
        print(
            f"[(R)] (A) = {err_a:.10e} (band {A_BAND}); (P) = {err_p:.10e} (band {P_BAND:.0e}); "
            f"(B) = {err_b:.3e} (band {B_BAND:.0e}); (C) = {err_c:.3e} (band {C_BAND:.0e})"
        )
        print(
            f"[(IT)] owned max|dT| / max|T - Tinf| = {it:.3e} (band L0_BAND {L0_BAND:g}); "
            f"theta_py max {th_max:.6f} K; Python Study + thermal {t_py:.1f} s"
        )
        print(
            f"[(c1)] control file source_scale {th2.cooling.source_scale!r}: (A) comparator "
            f"{c1:.4e} = {c1:.10e} (MAT-10 record {MAT10_B_MISSING_HALF_PRINT}; string equal "
            f"{f'{c1:.4e}' == MAT10_B_MISSING_HALF_PRINT}, printed not asserted; must exceed "
            f"A_BAND {A_BAND}); no re-solve {once_c1}"
        )
        print(f"[(c1')] source power control / file = {c1_ratio!r} (printed, not asserted)")
        print(f"[(X)] field {x_field!r} (expected {expected_field!r}); message: {x_message}")
        print(f"[(X)] run rc (min, max) {rc_run}, check rc {rc_check}; output dir exists {bad_exists}")
        print(f"[paraview] wrote {xdmf_path}")

    assert k_spec_orig == [] and k_cool_orig == [], (k_spec_orig, k_cool_orig)
    assert k_spec_ctl == [] and k_cool_ctl == [], (k_spec_ctl, k_cool_ctl)
    assert k_spec_py == [] and k_cool_py == [], (k_spec_py, k_cool_py)
    assert x_message is not None and x_field == expected_field, (x_field, expected_field)
    assert rc_run == (2, 2) and rc_check == (2, 2), (rc_run, rc_check)  # exit-code contract
    assert not bad_exists, "(X) output dir exists"
    assert rc_cli == (0, 0), rc_cli
    assert cli_summary["family"] == "sphere_in_box" and cli_summary["degree"] == study.degree
    assert cli_summary["ranks"] == comm.size
    assert cli_keys == TEMPERATURE_KEYS, cli_keys
    for key, rel in cli_rel.items():
        assert rel <= L0_BAND, f"(R-cli) {key} rel {rel:.3e} > {L0_BAND:g}"
    assert once_api, "(R) a second EM solve happened"
    assert n_sub == MAT10_SUB_CELLS, f"(R) sub cells {n_sub}"
    assert f"{err_a:.4e}" == MAT10_A_PRINT, f"(R) (A) {err_a:.4e}"
    assert f"{err_p:.4e}" == MAT10_P_PRINT, f"(R) (P) {err_p:.4e}"
    assert err_a <= A_BAND, f"(R) (A) {err_a:.4e} > {A_BAND}"
    assert err_p <= P_BAND, f"(R) (P) {err_p:.4e} > {P_BAND}"
    assert err_b <= B_BAND, f"(R) (B) {err_b:.3e}"
    assert err_c <= C_BAND, f"(R) (C) {err_c:.3e}"
    assert it <= L0_BAND, f"(IT) {it:.3e} > {L0_BAND:g}"
    assert once_c1, "(c1) re-solved"
    assert c1 > A_BAND, f"(c1) blind: {c1:.4e} <= {A_BAND}"

    if comm.rank == 0:
        print(
            f"\nAll anchors and controls hold. Total elapsed "
            f"{time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
