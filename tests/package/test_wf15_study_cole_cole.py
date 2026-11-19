"""`WF-15` step 8 (§9 item 653, U92): a Cole–Cole medium through `Study`.

A ``ColeColeMaterial`` handed to a ``Study`` on the imposed-field sphere is
resolved at the Study's frequency to the real (σ_eff, εᵣ) pair and gated
against `MAT-4`'s quasi-static closed form at that pair. Brief:
docs/planning/chunks/WF-15.md "§9 item 653".

Asserted:
(Ev-S) ``Study.epsilon_c()`` equals ``COLE_COLE.complex_permittivity(f)`` to
       ≤ 1e-12 relative (`MAT-3`'s (Ev) bound,
       test_mat3_cole_cole_plane_wave.py:109);
(Id0)  the zero-pole medium resolves to WF-15's gel (σ, εᵣ) to ≤ 1e-12;
(Ref)  a ColeColeMaterial on a port geometry raises TypeError before any
       mesh; on a ThermalBallGeometry the existing TypeError holds;
(Ceil) the two control ceilings exceed 2 × ``MAT4_BAND`` from closed forms;
(G-cc) one degree-2 Study solve within ``MAT4_BAND`` of the closed form at
       COLE_COLE's resolved pair (mean SAR, Im/Re ⟨E_z⟩), spread < SPREAD_BAND;
(c-D)  the same FEM SAR misses the Debye-reduced closed form by > MAT4_BAND;
(c-S)  the same FEM SAR misses the poles-dropped closed form by > MAT4_BAND.
Printed, not a control: ε̂ frozen at 2f, solved at f (closed forms only).
"""

from __future__ import annotations

import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    FourPortBirdcageGeometry,
    QuasiStaticSphereField,
    Study,
    ThermalBallGeometry,
    UniformHeatSource,
)
from fem_em_solver.materials import ColeColeMaterial
from tests.complex_mode import complex_only
from tests.materials.test_phantom_material_model import MAT4_BAND
from tests.package.test_wf15_study_lossy_sphere import (
    GEOMETRY,
    SPREAD_BAND,
    _closed,
    _ez_probe,
)
from tests.validation.test_lossy_sphere_sar import (
    E0,
    EPSILON_R_SPHERE,
    FREQUENCY_HZ,
    OMEGA,
    RHO_KG_M3,
    SIGMA_HIGH,
)
from tests.validation.test_mat3_cole_cole_plane_wave import (
    COLE_COLE,
    DEBYE_REDUCED,
    EPS_INF,
    SIGMA_IONIC,
)

# `MAT-3`'s (Ev) bound, by value with citation: test_mat3_cole_cole_plane_wave.py:109.
EV_BOUND = 1e-12
POLES_DROPPED = ColeColeMaterial(epsilon_inf=EPS_INF, sigma_ionic=SIGMA_IONIC)


def _study(material, geometry=GEOMETRY, **kwargs):
    return Study(
        geometry,
        material,
        QuasiStaticSphereField(e0=E0),
        frequency_hz=FREQUENCY_HZ,
        density_kg_m3=RHO_KG_M3,
        **kwargs,
    )


def _closed_at(material, frequency_hz=FREQUENCY_HZ):
    """`_closed` at ``material``'s (σ_eff, εᵣ) resolved at ``frequency_hz``."""
    sigma_eff, eps_r = material.sigma_epsilon_at(frequency_hz)
    return _closed(sigma_eff, eps_r), sigma_eff, eps_r


def _ceilings():
    c_cc = _closed_at(COLE_COLE)[0]["sar_w_per_kg"]
    c_d = _closed_at(DEBYE_REDUCED)[0]["sar_w_per_kg"]
    c_s = _closed_at(POLES_DROPPED)[0]["sar_w_per_kg"]
    c_2f = _closed_at(COLE_COLE, 2.0 * FREQUENCY_HZ)[0]["sar_w_per_kg"]
    return {
        "cc": c_cc,
        "D": c_d,
        "S": c_s,
        "2f": c_2f,
        "ceil_D": abs(c_cc - c_d) / c_d,
        "ceil_S": abs(c_cc - c_s) / c_s,
        "miss_2f": abs(c_cc - c_2f) / c_2f,
    }


@complex_only
def test_wf15_s8_resolution_refusals_and_ceilings():
    """W0: (Ev-S), (Id0), (Ref), (Ceil); no solve, no mesh."""
    comm = MPI.COMM_WORLD
    out = ["\n[WF-15 s8] Cole-Cole through Study, resolution and refusals (no solve)"]

    # (Ev-S)
    study = _study(COLE_COLE)
    eps_s = study.epsilon_c()
    eps_m = COLE_COLE.complex_permittivity(FREQUENCY_HZ)
    ev = abs(eps_s - eps_m) / abs(eps_m)
    assert study.dispersive_material is COLE_COLE
    assert _study("mid").dispersive_material is None
    out.append(
        f"  (Ev-S) Study.epsilon_c() = {eps_s!r}, COLE_COLE.complex_permittivity = {eps_m!r}, "
        f"rel {ev:.3e} (bound {EV_BOUND:.0e}); resolved sigma_eff = {study.material.sigma:.6g} "
        f"S/m, eps_r = {study.material.epsilon_r:.6g} (predicted 1.00980 / 181.279); "
        f"omega check {study.omega == OMEGA}"
    )
    assert ev <= EV_BOUND, f"(Ev-S) {ev:.3e}"

    # (Id0)
    zero = _study(ColeColeMaterial(epsilon_inf=EPSILON_R_SPHERE, sigma_ionic=SIGMA_HIGH))
    d_sig = abs(zero.material.sigma - SIGMA_HIGH) / SIGMA_HIGH
    d_eps = abs(zero.material.epsilon_r - EPSILON_R_SPHERE) / EPSILON_R_SPHERE
    out.append(
        f"  (Id0) zero-pole medium: sigma {zero.material.sigma!r} (rel {d_sig:.3e}), "
        f"eps_r {zero.material.epsilon_r!r} (rel {d_eps:.3e})"
    )
    assert d_sig <= EV_BOUND and d_eps <= EV_BOUND, f"(Id0) {d_sig:.3e} / {d_eps:.3e}"

    # (Ref)
    with pytest.raises(TypeError) as port_err:
        Study(FourPortBirdcageGeometry.f_small(), COLE_COLE, None, frequency_hz=FREQUENCY_HZ)
    assert "material" in str(port_err.value)
    with pytest.raises(TypeError) as ball_err:
        Study(ThermalBallGeometry(0.1, 0.005), COLE_COLE, UniformHeatSource(2000.0))
    assert "material" in str(ball_err.value)
    out.append(f"  (Ref) port: TypeError({port_err.value}); ball: TypeError({ball_err.value})")

    # (Ceil)
    c = _ceilings()
    out.append(
        f"  (Ceil) closed SAR: Cole-Cole {c['cc']:.6e} (predicted 3.98508e-08), Debye-reduced "
        f"{c['D']:.6e}, poles-dropped {c['S']:.6e} W/kg; ceilings (c-D) {c['ceil_D']:.4%} "
        f"(predicted 160.9 %), (c-S) {c['ceil_S']:.4%} (predicted 607.7 %), bound "
        f"2 x MAT4_BAND = {2.0 * MAT4_BAND:.0%}"
    )
    out.append(
        f"  printed, not a control: eps_hat frozen at 2f, solved at f, closed-form miss "
        f"{c['miss_2f']:.4%} (predicted 0.6 %)"
    )
    if comm.rank == 0:
        print("\n".join(out))
    assert c["ceil_D"] > 2.0 * MAT4_BAND, f"(Ceil) Debye-reduced {c['ceil_D']:.3%}"
    assert c["ceil_S"] > 2.0 * MAT4_BAND, f"(Ceil) poles-dropped {c['ceil_S']:.3%}"


@complex_only
def test_wf15_s8_cole_cole_sphere_closed_form():
    """W1: (G-cc), (c-D), (c-S) on one degree-2 Study solve."""
    comm = MPI.COMM_WORLD
    study = _study(COLE_COLE)
    assert study.degree == 2
    t0 = time.perf_counter()
    res = study.solve(comm)
    seconds = time.perf_counter() - t0
    sar_dict = study.mean_sar()
    sar = float(sar_dict["mean_sar_w_per_kg"])
    ez = _ez_probe(res.fields, comm)
    ez_mean = complex(np.mean(ez))
    n_cells = int(res.mesh.topology.index_map(res.mesh.topology.dim).size_global)

    c, sigma_eff, eps_r = _closed_at(COLE_COLE)
    err = abs(sar - c["sar_w_per_kg"]) / c["sar_w_per_kg"]
    im_re = ez_mean.imag / ez_mean.real
    im_re_err = abs(im_re - c["im_re"]) / abs(c["im_re"])
    spread = float(np.max(np.abs(ez - ez_mean)) / abs(ez_mean))
    c_d = _closed_at(DEBYE_REDUCED)[0]["sar_w_per_kg"]
    c_s = _closed_at(POLES_DROPPED)[0]["sar_w_per_kg"]
    miss_d = abs(sar - c_d) / c_d
    miss_s = abs(sar - c_s) / c_s
    if comm.rank == 0:
        print(
            f"\n[WF-15 s8] degree-{res.degree} Cole-Cole Study: {n_cells} cells, solve+mesh "
            f"{seconds:.1f} s; sigma_eff = {sigma_eff:.6g} S/m, eps_r = {eps_r:.6g}\n"
            f"  (G-cc) closed SAR = {c['sar_w_per_kg']:.6e}, FEM = {sar:.6e} err {err:.4%} "
            f"(predicted <= 0.2 %); Im/Re closed {c['im_re']:.5f}, FEM {im_re:.5f} err "
            f"{im_re_err:.4%}; spread {spread:.4%} (predicted <= 1 %) (bands {MAT4_BAND:.0%} / "
            f"{MAT4_BAND:.0%} / {SPREAD_BAND:.0%})\n"
            f"  (c-D) FEM vs Debye-reduced closed {c_d:.6e}: miss {miss_d:.4%} (> {MAT4_BAND:.0%})\n"
            f"  (c-S) FEM vs poles-dropped closed {c_s:.6e}: miss {miss_s:.4%} (> {MAT4_BAND:.0%})"
        )
    assert res.degree == 2
    assert err <= MAT4_BAND, f"(G-cc) mean SAR {err:.2%} off the closed form"
    assert im_re_err <= MAT4_BAND, f"(G-cc) Im/Re {im_re_err:.2%} off the closed form"
    assert spread < SPREAD_BAND, f"(G-cc) interior spread {spread:.2%}"
    assert miss_d > MAT4_BAND, f"(c-D) blind: FEM within {miss_d:.2%} of Debye-reduced"
    assert miss_s > MAT4_BAND, f"(c-S) blind: FEM within {miss_s:.2%} of poles-dropped"
