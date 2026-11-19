"""`WF-16` step 2 (§9 item 259, chain U9): the frequency sweep through the ``Study``.

One ``Study`` on item 192's F-small ×1 fixture (exact port rectangles, κ off)
is solved at 10 / 64 / 128 MHz on **one mesh** by ``Study.solve_sweep`` and
writes **one** multi-frequency ``.s4p``.  Each frequency's S, in memory and
parsed back from the file, is held to `ANS-4` step 5's registered
exact-rectangle record at that frequency within ``LEG_D0_REPRODUCTION_BAND``,
at the record's own degree and width (degree 1 at ``-n 2``: ``P1``; degree 2,
the Study's default, at ``-n 8``: ``P2_X1``).  `PORT-11`'s gates are imported
per frequency.  Every band is imported; none is set here.

Reuses step 1's knob ``FEM_EM_WF16_DEGREE`` (helpers only are imported from
step 1's module; never its module-scoped fixture, which would re-solve).
"""

from __future__ import annotations

import dataclasses
import inspect
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    FourPortBirdcageGeometry,
    LumpedPortDrive,
    QuasiStaticSphereField,
    SphereInBoxGeometry,
    Study,
)
from fem_em_solver.core import build_material_fields, build_mu_r_field
from fem_em_solver.materials import GelledSalinePhantomMaterial
from fem_em_solver.ports.touchstone import load_touchstone
from fem_em_solver.study import PortStudySweepResult
from fem_em_solver.utils.instrumentation import report_peak_rss
from tests.complex_mode import complex_only
from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL, FREQ_ATOL_HZ
from tests.materials.test_phantom_material_model import _dg0_owned_diff
from tests.mesh.test_birdcage_leg_offset import SHEET_AREA_BAND
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_tags import RING_RADIUS
from tests.package.test_wf16_study_four_port import (
    DEGREE2_RECORD_WIDTH,
    OUT_DIR,
    _degree,
    _p0,
    _rel,
    needs_degree,
)
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
    ANS4_EXACT_STRIP_S_P1_RECORDS,
    ANS4_EXACT_STRIP_S_P2_X1,
    FREQUENCIES,
    SPEC_AREA_M2,
    SPEC_HEIGHT_M,
    SPEC_WIDTH_M,
)
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    FREQUENCY_128_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_ops61_exact_strip_ports import (
    ANS4_EXACT_STRIP_CELL_RECORD,
    STRIP_BAND,
)
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    _circulant_classes,
    _class_spread,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_gap_voltage_impedance import FREQUENCY_HZ
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

KEYS = tuple(FREQUENCIES)

# (c4) hand arithmetic from the brief (item 259), ||rec_i - rec_j|| / ||rec_j||
# on the registered records; the in-test matrix must reproduce it to two
# digits, otherwise the comparator is mis-wired.  Not a band.
C4_HAND_TABLE = {
    1: {
        ("10", "64"): 0.8103,
        ("10", "128"): 1.2367,
        ("64", "10"): 0.8120,
        ("64", "128"): 0.5524,
        ("128", "10"): 1.2763,
        ("128", "64"): 0.5689,
    },
    2: {
        ("10", "64"): 0.8395,
        ("10", "128"): 1.2592,
        ("64", "10"): 0.8515,
        ("64", "128"): 0.5594,
        ("128", "10"): 1.3138,
        ("128", "64"): 0.5755,
    },
}
C4_TWO_DIGITS = 0.005


def _records(degree):
    table = ANS4_EXACT_STRIP_S_P1_RECORDS if degree == 1 else ANS4_EXACT_STRIP_S_P2_X1
    return {k: np.asarray(table[k], dtype=np.complex128) for k in KEYS}


def _phantom():
    return GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)


def _port_study(degree=None):
    kwargs = {"frequency_hz": FREQUENCY_64_HZ}
    if degree == 1:
        kwargs["degree"] = 1  # degree 2 is the Study's default: no argument.
    phantom = _phantom()
    return (
        Study(FourPortBirdcageGeometry.f_small(), phantom, LumpedPortDrive(), **kwargs),
        phantom,
    )


@pytest.fixture(scope="module")
def sweep_leg():
    degree = _degree()
    comm = MPI.COMM_WORLD
    width = ANS4_EXACT_STRIP_RECORD_RANK_WIDTH if degree == 1 else DEGREE2_RECORD_WIDTH
    # Width guard, asserted (a record holds at its own width, `OPS-41`).
    assert (
        comm.size == width
    ), f"degree {degree} record is asserted at -n {width}; got -n {comm.size}"
    study, phantom = _port_study(degree)
    t0 = time.perf_counter()
    sw = study.solve_sweep(tuple(FREQUENCIES.values()), comm)
    t_total = time.perf_counter() - t0
    path = OUT_DIR / f"wf16_sweep_f_small_p{degree}_n{comm.size}.s4p"
    t0 = time.perf_counter()
    written = sw.write_touchstone(path)
    t_ts = time.perf_counter() - t0
    _p0(
        f"\n[WF-16 s2] degree {degree} -n {comm.size}: solve_sweep {t_total:.2f} s "
        f"(mesh {sw.mesh_time:.2f} s; sweeps "
        + " / ".join(f"{k} MHz {t:.2f} s" for k, t in zip(KEYS, sw.sweep_times))
        + f"); Touchstone {t_ts:.3f} s"
    )
    report_peak_rss(comm, label="[WF-16 s2] summed peak ru_maxrss")
    return {
        "degree": degree,
        "study": study,
        "phantom": phantom,
        "sw": sw,
        "path": written,
        "rec": _records(degree),
        "other": _records(3 - degree),
    }


def test_wf16_sweep_api_and_refusals():
    """(K′) scalar constructor pinned (`conductor` added by `df25001f`, ruled
    2026-10-09); refusals before any mesh; retag exact."""
    params = inspect.signature(Study.__init__).parameters
    assert tuple(params) == (
        "self",
        "geometry",
        "material",
        "excitation",
        "frequency_hz",
        "degree",
        "density_kg_m3",
        "conductor",
    )
    assert params["degree"].default == 2
    assert params["density_kg_m3"].default == 1000.0
    assert params["conductor"].default is None
    assert params["conductor"].kind is inspect.Parameter.KEYWORD_ONLY
    study, phantom = _port_study()
    bad = (
        (),
        (FREQUENCY_128_HZ, FREQUENCY_HZ),
        (FREQUENCY_64_HZ, FREQUENCY_64_HZ),
        (0.0,),
        (float("nan"),),
    )
    for freqs in bad:
        with pytest.raises(ValueError):
            study.solve_sweep(freqs)
    assert not hasattr(study, "sweep_result")
    assert study.result is None
    sphere = Study(
        SphereInBoxGeometry(0.05, 0.1, 0.01, 0.02),
        _phantom(),
        QuasiStaticSphereField(),
        frequency_hz=FREQUENCY_64_HZ,
    )
    with pytest.raises(TypeError):
        sphere.solve_sweep(tuple(FREQUENCIES.values()))
    for f in FREQUENCIES.values():
        retag = dataclasses.replace(phantom, frequency_hz=f)
        assert retag.sigma == phantom.sigma
        assert retag.epsilon_r == phantom.epsilon_r
        assert retag.mu_r == phantom.mu_r
        assert retag.frequency_hz == f
    assert hasattr(Study, "solve_sweep")
    assert dataclasses.is_dataclass(PortStudySweepResult)


@complex_only
@needs_degree
def test_wf16_sweep_one_mesh_per_frequency_problems(sweep_leg):
    """(F′) one mesh, SPEC rectangles; (U′) per-frequency problems; state unmutated."""
    comm = MPI.COMM_WORLD
    sw, study, phantom = sweep_leg["sw"], sweep_leg["study"], sweep_leg["phantom"]
    r0 = sw.results[0]
    assert len(sw.results) == len(FREQUENCIES)
    assert all(r.mesh is r0.mesh for r in sw.results)
    assert isinstance(sw.mesh_time, float)
    assert len(sw.sweep_times) == len(FREQUENCIES)
    # (F′)
    ratio = r0.cells / ANS4_EXACT_STRIP_CELL_RECORD
    _p0(f"[WF-16 s2] (F') cells {r0.cells} / record {ANS4_EXACT_STRIP_CELL_RECORD}")
    assert abs(ratio - 1.0) <= CELL_COUNT_BAND
    for s in r0.sheets:
        _p0(
            f"    sheet {s['tag']}: {s['facets']} facets area {s['area']:.12e} "
            f"w {s['w']:.12e} h {s['h']:.12e} centre {s['centre_radial']:.12e}"
        )
        assert abs(s["area"] / SPEC_AREA_M2 - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["w"] / SPEC_WIDTH_M - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["w_bbox"] / SPEC_WIDTH_M - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["h"] / SPEC_HEIGHT_M - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["centre_radial"] / RING_RADIUS - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["area_ratio_full"] - 1.0) < SHEET_AREA_BAND, s["tag"]
    # (U′)
    fields_by_key = {}
    for i, (k, f) in enumerate(FREQUENCIES.items()):
        r = sw.results[i]
        assert sw.frequencies_hz[i] == f
        assert r.problem.frequency_hz == f
        assert r.sweep.frequency_hz == f
        assert r.problem.phantom_material == dataclasses.replace(phantom, frequency_hz=f)
        assert r.sweep.is_placeholder is False
        assert r.degree == sweep_leg["degree"]
        assert len(r.specs) == len(r0.specs) == 4
        for a, b in zip(r.specs, r0.specs):
            assert a == b
        for a, b in zip(r.port_definitions, r0.port_definitions):
            assert a == b
        p = r.problem
        sig, eps = build_material_fields(
            r.mesh,
            p.material,
            cell_tags=r.cell_tags,
            material_map=p.material_map,
            phantom_material=p.phantom_material,
            phantom_tag=p.phantom_tag,
        )
        mu = build_mu_r_field(
            r.mesh, p.material, cell_tags=r.cell_tags, material_map=p.material_map
        )
        fields_by_key[k] = (sig, eps, mu)
    ref = fields_by_key["64"]
    for k, (sig, eps, mu) in fields_by_key.items():
        ds = _dg0_owned_diff(sig, ref[0], comm)
        de = _dg0_owned_diff(eps, ref[1], comm)
        dm = _dg0_owned_diff(mu, ref[2], comm)
        _p0(f"[WF-16 s2] (U') {k} MHz vs 64 MHz DG0 max|d| sigma {ds}, eps_r {de}, mu_r {dm}")
        assert ds[0] == 0.0 and de[0] == 0.0 and dm[0] == 0.0
    # Scalar state not mutated.
    assert study.frequency_hz == FREQUENCY_64_HZ
    assert study.material is phantom
    assert study.result is None
    assert study.sweep_result is sw


@complex_only
@needs_degree
def test_wf16_sweep_s_equals_records(sweep_leg):
    """(R′) S_k == rec[k] per frequency; (c1′) cross-degree record per frequency."""
    d = sweep_leg["degree"]
    s_all = sweep_leg["sw"].s_matrices
    assert s_all.shape == (len(KEYS), 4, 4) and s_all.dtype == np.complex128
    for i, k in enumerate(KEYS):
        rel = _rel(s_all[i], sweep_leg["rec"][k])
        c1 = _rel(s_all[i], sweep_leg["other"][k])
        _p0(
            f"[WF-16 s2] (R'{d}) {k} MHz ||S - S_rec||/||S_rec|| {rel:.3e} (band "
            f"{LEG_D0_REPRODUCTION_BAND:.0e}); (c1') vs degree {3 - d} record {c1:.4f}"
        )
        assert rel <= LEG_D0_REPRODUCTION_BAND, k
        assert c1 > LEG_D0_REPRODUCTION_BAND, k


@complex_only
@needs_degree
def test_wf16_sweep_touchstone_equals_records(sweep_leg):
    """(T1′) one file == in-memory (both readers); (T2′) vs records; (c4)–(c6)."""
    import skrf

    d = sweep_leg["degree"]
    sw, rec = sweep_leg["sw"], sweep_leg["rec"]
    freqs, mats, ports = load_touchstone(sweep_leg["path"])
    loaded = np.asarray(mats, dtype=np.complex128)
    ntw = skrf.Network(str(sweep_leg["path"]))
    sk = np.asarray(ntw.s, dtype=np.complex128)
    # (T1′)
    assert len(freqs) == len(KEYS) and len(ntw.f) == len(KEYS)
    assert loaded.shape == (len(KEYS), 4, 4) and sk.shape == (len(KEYS), 4, 4)
    for i, f in enumerate(FREQUENCIES.values()):
        assert abs(float(freqs[i]) - f) <= FREQ_ATOL_HZ
        assert abs(float(ntw.f[i]) - f) <= FREQ_ATOL_HZ
    assert ntw.nports == 4 and tuple(ports) == ("P1", "P2", "P3", "P4")
    assert np.all(np.asarray(ntw.z0) == 50.0)
    for i, k in enumerate(KEYS):
        d_load = float(np.max(np.abs(loaded[i] - sw.s_matrices[i])))
        d_skrf = float(np.max(np.abs(sk[i] - sw.s_matrices[i])))
        _p0(f"[WF-16 s2] (T1') {k} MHz max|dS| load_touchstone {d_load:.3e}, skrf {d_skrf:.3e}")
        assert d_load <= ANCHOR_TOL and d_skrf <= ANCHOR_TOL, k
    # (T2′) and (c4): the confusion matrix, printed before it is asserted.
    conf = {(ki, kj): _rel(loaded[i], rec[kj]) for i, ki in enumerate(KEYS) for kj in KEYS}
    _p0(f"[WF-16 s2] (c4) confusion ||S_file,i - rec_j||/||rec_j|| (degree {d}):")
    for ki in KEYS:
        _p0("    " + ki + " -> " + "  ".join(f"{kj}: {conf[(ki, kj)]:.4e}" for kj in KEYS))
    for ki in KEYS:
        assert conf[(ki, ki)] <= LEG_D0_REPRODUCTION_BAND, ("T2'", ki)
    hand = C4_HAND_TABLE[d]
    for (ki, kj), v in conf.items():
        if ki == kj:
            continue
        assert v > LEG_D0_REPRODUCTION_BAND, ("c4", ki, kj)
        assert abs(v - hand[(ki, kj)]) <= C4_TWO_DIGITS, ("c4 mis-wired", ki, kj, v)
    # (c5) planted 10 <-> 128 MHz block swap, frequency labels kept.
    swapped = loaded.copy()
    i10, i128 = KEYS.index("10"), KEYS.index("128")
    swapped[[i10, i128]] = swapped[[i128, i10]]
    c5 = {k: _rel(swapped[i], rec[k]) for i, k in enumerate(KEYS)}
    _p0(
        "[WF-16 s2] (c5) swapped blocks vs records: "
        + ", ".join(f"{k} {v:.4e}" for k, v in c5.items())
    )
    assert c5["10"] > LEG_D0_REPRODUCTION_BAND
    assert c5["128"] > LEG_D0_REPRODUCTION_BAND
    assert c5["64"] <= LEG_D0_REPRODUCTION_BAND
    # (c6) band resolution: delta = 3 * band * ||rec_64|| on entry (0, 1).
    i64 = KEYS.index("64")
    planted = loaded[i64].copy()
    planted[0, 1] += 3.0 * LEG_D0_REPRODUCTION_BAND * np.linalg.norm(rec["64"])
    c6 = _rel(planted, rec["64"])
    _p0(f"[WF-16 s2] (c6) planted 3x-band delta: {c6:.4e} (band {LEG_D0_REPRODUCTION_BAND:.0e})")
    assert c6 > LEG_D0_REPRODUCTION_BAND


@complex_only
@needs_degree
def test_wf16_sweep_port11_gates(sweep_leg):
    """(G′) `PORT-11` (i), (ii), (iii′) per frequency on each result's own S and Z."""
    for k, r in zip(KEYS, sweep_leg["sw"].results):
        s, z = r.s_matrix, r.z_matrix
        recip = float(np.linalg.norm(s - s.T) / np.linalg.norm(s))
        sigma = np.linalg.svd(s, compute_uv=False)
        power = np.sum(np.abs(s) ** 2, axis=0)
        classes = _circulant_classes(z)
        spreads = {n: _class_spread(v) for n, v in classes.items()}
        pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
        worst = max(spreads.values())
        separation = pooled / worst if worst > 0.0 else np.inf
        _p0(
            f"[WF-16 s2] (G') {k} MHz reciprocity {recip:.3e}; sigma_max "
            f"{float(np.max(sigma)):.12f}; max column power {float(np.max(power)):.12f}; "
            "spreads "
            + ", ".join(f"{n} {v * 100:.4f}%" for n, v in spreads.items())
            + f"; separation {separation:.2f}x"
        )
        assert recip <= RECIPROCITY_BAND, k
        assert float(np.max(sigma)) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, k
        for v in power:
            assert v <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, k
        for name, v in spreads.items():
            assert v <= ADJACENT_SPREAD_BAND, (k, name)
        assert separation >= POOLED_SEPARATION_FLOOR, k
