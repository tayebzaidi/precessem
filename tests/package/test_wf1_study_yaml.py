"""`WF-1` step 1 (§9 item 256, chain U5): ``femem run <file>.yaml``.

A schema-v0 YAML file parses into the same ``Study`` `WF-15` / `WF-16` build in
Python (brief: docs/planning/chunks/WF-1.md, "§9 item 256").  Asserted:

(K)   the YAML-built Study equals the Python-built one field by field, floats
      by ``float.hex`` (sphere, F-small, and the "mid" preset route);
(Deg) degree 2 when ``study.degree`` is omitted;
(X)   three malformed files raise ``StudyConfigError`` naming the field, and
      ``cli.main(["run", ...])`` returns 2 on every rank with nothing written;
      plus the duplicate-key guard;
(U)   DG0 sigma / eps_r / mu_r from the YAML material == the Python one;
(S)   the CLI's summary.json is the run's (rc, family, degree, ranks, cells);
(R)   `MAT-4`'s / item 192's records through the CLI;
(T)   the CLI's Touchstone file (both readers; record);
(I)   the CLI run vs the Python Study in the same process;
(G)   `MAT-4`'s band / `PORT-11`'s gates on the CLI's outputs.

Negative controls (asserted): (c1) one-ULP sigma through the parser and on the
DG0 field; (c2) sphere vs the degree-2 record; (c3) birdcage vs the other
degree's record; (c4) Touchstone P1<->P2 swap; (c5) good files check 0.

No solve output is compared bitwise (known-issues `MAT-1` entry,
``SOLVE_FLOOR_RECORD``); bitwise lives upstream of the solve only.
"""

from __future__ import annotations

import copy
import dataclasses
import hashlib
import json
import shutil
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import (
    FourPortBirdcageGeometry,
    GelledSalinePhantomMaterial,
    LumpedPortDrive,
    Study,
    StudyConfigError,
    cli,
    load_study,
)
from fem_em_solver.core import HomogeneousMaterial, build_material_fields, build_mu_r_field
from fem_em_solver.ports.touchstone import load_touchstone
from fem_em_solver.study_yaml import _study_from_mapping
from fem_em_solver.utils.instrumentation import report_peak_rss
from tests.complex_mode import complex_only
from tests.io.test_touchstone_skrf_crosscheck import ANCHOR_TOL, FREQ_ATOL_HZ
from tests.materials.test_phantom_material_model import MAT4_BAND, _dg0_owned_diff
from tests.package.test_wf15_study_lossy_sphere import (
    GEOMETRY,
    MAT4_DEG1_PRINT,
    MAT4_DEG1_RECORD,
    _closed,
    _study,
)
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    SALINE_EPSILON_R,
    SALINE_SIGMA,
)
from tests.validation.test_lossy_sphere_sar import EPSILON_R_SPHERE, SIGMA_HIGH, SPHERE_TAG
from tests.validation.test_mat5_em_coupled_sphere import L0_BAND
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_REPRODUCTION_BAND,
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    _circulant_classes,
    _class_spread,
)
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

REPO = Path(__file__).resolve().parents[2]
WF1 = Path(__file__).resolve().parent / "wf1"
LOGS = REPO / "logs" / "wf1"  # gitignored (.gitignore:82)
SPHERE_FILE = WF1 / "sphere_mat4.yaml"
BIRDCAGE_FILE = WF1 / "f_small_four_port.yaml"
CONTROL_FILE = WF1 / "control_sigma_ulp.yaml"
DUPLICATE_FILE = WF1 / "malformed_duplicate_key.yaml"
MALFORMED = {
    "malformed_unknown_key.yaml": "geometry.sphere_raduis",
    "malformed_wrong_type.yaml": "study.frequency_hz",
    "malformed_missing_section.yaml": "excitation",
}
# `WF-15`'s degree-2 Study record on this fixture
# (`20261005T170537Z_WF-15-r2-W2a-deg2-study.log`, PROJECT_PLAN.md `WF-15` row),
# restated with its citation as `WF-15` restates MAT4_DEG1_RECORD; (c2) only.
WF15_DEG2_RECORD = 8.0116577245e-08


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _py_birdcage(**kwargs):
    """`test_wf16_study_four_port.py:139–145` at degree 1."""
    phantom = GelledSalinePhantomMaterial(SALINE_SIGMA, SALINE_EPSILON_R, FREQUENCY_64_HZ)
    return Study(
        FourPortBirdcageGeometry.f_small(),
        phantom,
        LumpedPortDrive(),
        frequency_hz=FREQUENCY_64_HZ,
        degree=1,
        **kwargs,
    )


def _cmp(x, y, name, out):
    if type(x) is not type(y):
        out.append(name)
    elif dataclasses.is_dataclass(x):
        for f in dataclasses.fields(x):
            _cmp(getattr(x, f.name), getattr(y, f.name), f"{name}.{f.name}", out)
    elif isinstance(x, float):
        if x.hex() != y.hex():
            out.append(name)
    elif isinstance(x, complex):
        if (x.real.hex(), x.imag.hex()) != (y.real.hex(), y.imag.hex()):
            out.append(name)
    elif isinstance(x, tuple):
        sub = []
        if len(x) != len(y):
            out.append(name)
            return
        for i, (a, b) in enumerate(zip(x, y)):
            _cmp(a, b, f"{name}[{i}]", sub)
        if sub:
            out.append(name)
    elif x != y:
        out.append(name)


def _spec_diff(a: Study, b: Study) -> list:
    """Dotted names that differ between two Study specs (floats by float.hex)."""
    out = []
    for attr in ("geometry", "material", "excitation"):
        _cmp(getattr(a, attr), getattr(b, attr), attr, out)
    for attr in ("frequency_hz", "degree", "density_kg_m3"):
        _cmp(getattr(a, attr), getattr(b, attr), attr, out)
    # `WF-1` step 4 (§9 item 343): None vs a LeontovichWall differs by type.
    _cmp(a.conductor, b.conductor, "conductor", out)
    return out


def _mapping(path):
    import yaml

    return yaml.safe_load(Path(path).read_text())


def _fresh_dir(path, comm):
    if comm.rank == 0 and path.exists():
        shutil.rmtree(path)
    comm.Barrier()


def _rc_all(rc, comm):
    return comm.allreduce(int(rc), op=MPI.MIN), comm.allreduce(int(rc), op=MPI.MAX)


def test_wf1_good_files_build_the_python_studies():
    """(K), (Deg), (c1) parse half, (c5) good half, duplicate-key guard."""
    import yaml

    comm = MPI.COMM_WORLD
    _p0(f"\n[WF-1] yaml.__version__ = {yaml.__version__}")
    for f in sorted(WF1.glob("*.yaml")):
        _p0(f"  sha256 {f.name} {hashlib.sha256(f.read_bytes()).hexdigest()}")

    # (K) sphere
    s_yaml, s_py = load_study(SPHERE_FILE), _study(degree=1)
    d_sphere = _spec_diff(s_yaml, s_py)
    # (K) F-small
    b_yaml, b_py = load_study(BIRDCAGE_FILE), _py_birdcage()
    d_cage = _spec_diff(b_yaml, b_py)
    # (K) preset route through the private mapping entry point (no extra file)
    m = _mapping(SPHERE_FILE)
    m["material"] = {"preset": "mid"}
    p_yaml, p_py = _study_from_mapping(m), _study(material="mid", degree=1)
    d_preset = _spec_diff(p_yaml, p_py)
    _p0(f"  (K) sphere diff {d_sphere}; F-small diff {d_cage}; preset 'mid' diff {d_preset}")
    for label, a, b in (("yaml", s_yaml, s_py), ("python", s_py, s_yaml)):
        _p0(
            f"  (K) {label}: resolution_sphere {a.geometry.resolution_sphere.hex()}, "
            f"resolution_far {a.geometry.resolution_far.hex()}, frequency_hz "
            f"{a.frequency_hz.hex()}, sigma {a.material.sigma.hex()}"
        )
    assert d_sphere == [], d_sphere
    assert d_cage == [], d_cage
    assert d_preset == [], d_preset
    assert GEOMETRY == s_yaml.geometry

    # (Deg) degree 2 when omitted, both good mappings.
    degs = []
    for f in (SPHERE_FILE, BIRDCAGE_FILE):
        mm = copy.deepcopy(_mapping(f))
        del mm["study"]["degree"]
        degs.append(_study_from_mapping(mm).degree)
    _p0(f"  (Deg) degree with study.degree omitted: {degs}")
    assert degs == [2, 2]

    # (c1) parse half: one ULP through the parser, nothing else moves.
    c_yaml = load_study(CONTROL_FILE)
    d_c1 = _spec_diff(c_yaml, s_py)
    _p0(
        f"  (c1) control diff {d_c1}; sigma {c_yaml.material.sigma!r} == nextafter(0.57, 1): "
        f"{c_yaml.material.sigma == float(np.nextafter(0.57, 1.0))}"
    )
    assert float("0.5700000000000001") == float(np.nextafter(0.57, 1.0))
    assert d_c1 == ["material.sigma"], d_c1

    # Duplicate-key guard (printed and asserted beside the three).
    with pytest.raises(StudyConfigError) as ei:
        load_study(DUPLICATE_FILE)
    _p0(f"  duplicate key: field {ei.value.field!r}; {ei.value}")
    assert type(ei.value) is StudyConfigError
    assert ei.value.field == "material.sigma"

    # (c5) good half: check returns 0 on the good files and the control.
    rcs = {f.name: cli.main(["check", str(f)]) for f in (SPHERE_FILE, BIRDCAGE_FILE, CONTROL_FILE)}
    _p0(f"  (c5) cli check rc on good files: {rcs}")
    for name, rc in rcs.items():
        assert _rc_all(rc, comm) == (0, 0), name


def test_wf1_malformed_files_name_their_field():
    """(X): three malformed files, each with its named field; (c5) bad half."""
    comm = MPI.COMM_WORLD
    out = LOGS / "malformed"
    fields = []
    for name, expected in MALFORMED.items():
        path = WF1 / name
        with pytest.raises(StudyConfigError) as ei:
            load_study(path)
        exc = ei.value
        _p0(f"\n[WF-1] (X) {name}: field {exc.field!r} (expected {expected!r}); {exc}")
        assert type(exc) is StudyConfigError
        assert exc.field == expected
        assert expected in str(exc)
        assert exc.source.endswith(name)
        fields.append(exc.field)

        _fresh_dir(out, comm)
        rc_run = cli.main(["run", str(path), "--output-dir", str(out)])
        rc_check = cli.main(["check", str(path)])
        run_all, check_all = _rc_all(rc_run, comm), _rc_all(rc_check, comm)
        exists = comm.allreduce(bool(out.exists()), op=MPI.LOR)
        _p0(f"  (X) {name}: run rc (min, max) {run_all}, check rc {check_all}, out exists {exists}")
        assert run_all == (2, 2)
        assert check_all == (2, 2)
        assert not exists
    assert len(set(fields)) == 3, fields


@complex_only
def test_wf1_cli_sphere_reproduces_wf15():
    """W1: (S), (U), (R), (I), (G), (c1) downstream, (c2)."""
    comm = MPI.COMM_WORLD
    width = ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
    assert comm.size == width, f"record width -n {width}; got -n {comm.size}"
    out = LOGS / SPHERE_FILE.stem
    _fresh_dir(out, comm)
    t0 = time.perf_counter()
    rc = cli.main(["run", str(SPHERE_FILE), "--output-dir", str(out)])
    t_cli = time.perf_counter() - t0
    assert _rc_all(rc, comm) == (0, 0)
    with open(out / "summary.json") as fh:
        summary = json.load(fh)

    t0 = time.perf_counter()
    study_py = _study(degree=1)
    res = study_py.solve(comm)
    sar_py = float(study_py.mean_sar()["mean_sar_w_per_kg"])
    t_py = time.perf_counter() - t0
    n_cells = int(res.mesh.topology.index_map(res.mesh.topology.dim).size_global)
    _p0(
        f"\n[WF-1] sphere: CLI run {t_cli:.1f} s (summary elapsed {summary['elapsed_s']:.1f} s), "
        f"Python Study {t_py:.1f} s; cells CLI {summary['cells']} / Python {n_cells} "
        f"(version-tagged record 74020)"
    )
    report_peak_rss(comm, label="[WF-1] summed peak ru_maxrss")

    # (S)
    assert summary["femem"] == 1 and summary["degree"] == 1 and summary["ranks"] == 2
    assert summary["family"] == "sphere_in_box"
    assert summary["cells"] == n_cells

    # (U) DG0 from the CLI's own parse vs the Python material, on the solved mesh.
    background = HomogeneousMaterial(sigma=0.0, epsilon_r=1.0)

    def _dg0(material):
        sig, eps = build_material_fields(
            res.mesh,
            background,
            cell_tags=res.cell_tags,
            phantom_material=material,
            phantom_tag=SPHERE_TAG,
        )
        mu = build_mu_r_field(res.mesh, background, cell_tags=res.cell_tags)
        return sig, eps, mu

    yaml_f = _dg0(load_study(SPHERE_FILE).material)
    py_f = _dg0(_study(degree=1).material)
    u = [_dg0_owned_diff(a, b, comm) for a, b in zip(yaml_f, py_f)]
    _p0(f"  (U) DG0 sigma {u[0]}, eps_r {u[1]}, mu_r {u[2]} (max|d|, differing, owned)")
    for dmax, ndiff, nown in u:
        assert dmax == 0.0 and ndiff == 0 and nown == n_cells
    # (c1) downstream: the control file's sigma reads exactly 2^-53.
    c_sig = _dg0(load_study(CONTROL_FILE).material)[0]
    c1 = _dg0_owned_diff(c_sig, py_f[0], comm)
    _p0(f"  (c1) control sigma: max|d| {c1[0]!r} (2^-53 = {2.0**-53!r}) on {c1[1]} owned dofs")
    assert c1[0] == 2.0**-53 and c1[1] > 0

    sar_yaml = summary["mean_sar"]["mean_sar_w_per_kg"]
    # (R)
    r1 = abs(sar_yaml - MAT4_DEG1_RECORD) / MAT4_DEG1_RECORD
    r2 = f"{sar_yaml:.6e}"
    # (I)
    i_rel = abs(sar_yaml - sar_py) / sar_py
    # (G)
    c = _closed(SIGMA_HIGH, EPSILON_R_SPHERE)
    g_err = abs(sar_yaml - c["sar_w_per_kg"]) / c["sar_w_per_kg"]
    # (c2)
    c2 = abs(sar_yaml / WF15_DEG2_RECORD - 1.0)
    _p0(
        f"  (R1) SAR_yaml {sar_yaml!r} rel to MAT4_DEG1_RECORD {r1:.3e} (L0_BAND {L0_BAND:.0e}); "
        f"(R2) '{r2}' vs '{MAT4_DEG1_PRINT}'\n"
        f"  (I) SAR_yaml vs SAR_py {sar_py!r}: rel {i_rel:.3e}\n"
        f"  (G) closed {c['sar_w_per_kg']:.6e}: err {g_err:.3%} (MAT4_BAND {MAT4_BAND:.0%})\n"
        f"  (c2) |SAR_yaml / WF15_DEG2_RECORD - 1| = {c2:.4%} ({c2 / L0_BAND:.2e} x L0_BAND)"
    )
    assert r1 <= L0_BAND
    assert r2 == MAT4_DEG1_PRINT
    assert i_rel <= L0_BAND
    assert g_err <= MAT4_BAND
    assert c2 > L0_BAND


@complex_only
def test_wf1_cli_f_small_reproduces_wf16():
    """W2: (S), (R), (T1), (T2), (I), (G), (c3), (c4)."""
    import skrf

    comm = MPI.COMM_WORLD
    width = ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
    assert comm.size == width, f"record width -n {width}; got -n {comm.size}"
    out = LOGS / BIRDCAGE_FILE.stem
    _fresh_dir(out, comm)
    t0 = time.perf_counter()
    rc = cli.main(["run", str(BIRDCAGE_FILE), "--output-dir", str(out)])
    t_cli = time.perf_counter() - t0
    assert _rc_all(rc, comm) == (0, 0)
    with open(out / "summary.json") as fh:
        summary = json.load(fh)
    s_yaml = np.array([[complex(*p) for p in row] for row in summary["s_matrix"]])
    z_yaml = np.array([[complex(*p) for p in row] for row in summary["z_matrix"]])

    t0 = time.perf_counter()
    res = _py_birdcage().solve(comm)
    t_py = time.perf_counter() - t0
    _p0(
        f"\n[WF-1] F-small: CLI run {t_cli:.1f} s (summary elapsed {summary['elapsed_s']:.1f} s), "
        f"Python Study {t_py:.1f} s (mesh {res.mesh_time:.1f} s, sweep {res.sweep_time:.1f} s); "
        f"cells CLI {summary['cells']} / Python {res.cells} (version-tagged record 117206)"
    )
    report_peak_rss(comm, label="[WF-1] summed peak ru_maxrss")

    # (S)
    assert summary["femem"] == 1 and summary["degree"] == 1 and summary["ranks"] == 2
    assert summary["family"] == "four_port_birdcage"
    assert summary["cells"] == res.cells
    assert summary["touchstone"] == "f_small_four_port.s4p"
    assert s_yaml.shape == (4, 4)

    p1 = _record(1)
    # (R1), (c3)
    r1 = _rel(s_yaml, p1)
    c3 = _rel(s_yaml, _record(2))
    # (T1)
    ts = out / summary["touchstone"]
    freqs, mats, ports = load_touchstone(ts)
    loaded = np.asarray(mats[0], dtype=np.complex128)
    ntw = skrf.Network(str(ts))
    sk = np.asarray(ntw.s[0], dtype=np.complex128)
    d_load = float(np.max(np.abs(loaded - s_yaml)))
    d_skrf = float(np.max(np.abs(sk - s_yaml)))
    # (T2), (c4)
    t2 = _rel(loaded, p1)
    perm = [1, 0, 2, 3]
    c4 = _rel(loaded[np.ix_(perm, perm)], p1)
    # (I)
    i_rel = _rel(s_yaml, res.s_matrix)
    # (G)
    recip = float(np.linalg.norm(s_yaml - s_yaml.T) / np.linalg.norm(s_yaml))
    sigma = np.linalg.svd(s_yaml, compute_uv=False)
    power = np.sum(np.abs(s_yaml) ** 2, axis=0)
    classes = _circulant_classes(z_yaml)
    spreads = {n: _class_spread(v) for n, v in classes.items()}
    pooled = _class_spread(np.concatenate([classes["adjacent"], classes["opposite"]]))
    worst = max(spreads.values())
    separation = pooled / worst if worst > 0.0 else np.inf
    _p0(
        f"  (R1) ||S_yaml - P1||/||P1|| {r1:.3e} (band {LEG_D0_REPRODUCTION_BAND:.0e}); "
        f"(c3) vs P2_X1 {c3:.4f}\n"
        f"  (T1) max|dS| load_touchstone {d_load:.3e}, skrf {d_skrf:.3e} (ANCHOR_TOL "
        f"{ANCHOR_TOL:.0e}); f {float(freqs[0])!r}; ports {tuple(ports)}; nports {ntw.nports}\n"
        f"  (T2) ||S_file - P1||/||P1|| {t2:.3e}; (c4) P1<->P2 swap {c4:.4f}\n"
        f"  (I) ||S_yaml - S_py||/||S_py|| {i_rel:.3e}\n"
        f"  (G) reciprocity {recip:.3e}; sigma_max {float(np.max(sigma)):.12f}; max column "
        f"power {float(np.max(power)):.12f}; spreads "
        + ", ".join(f"{k} {v * 100:.4f}%" for k, v in spreads.items())
        + f"; separation {separation:.2f}x"
    )
    assert r1 <= LEG_D0_REPRODUCTION_BAND
    assert c3 > LEG_D0_REPRODUCTION_BAND
    assert d_load <= ANCHOR_TOL and d_skrf <= ANCHOR_TOL
    assert abs(float(freqs[0]) - FREQUENCY_64_HZ) <= FREQ_ATOL_HZ
    assert abs(float(ntw.f[0]) - FREQUENCY_64_HZ) <= FREQ_ATOL_HZ
    assert ntw.nports == 4 and tuple(ports) == ("P1", "P2", "P3", "P4")
    assert np.all(np.asarray(ntw.z0) == 50.0)
    assert t2 <= LEG_D0_REPRODUCTION_BAND
    assert c4 > LEG_D0_REPRODUCTION_BAND
    assert i_rel <= LEG_D0_REPRODUCTION_BAND
    assert recip <= RECIPROCITY_BAND
    assert float(np.max(sigma)) <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for v in power:
        assert v <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for name, v in spreads.items():
        assert v <= ADJACENT_SPREAD_BAND, name
    assert separation >= POOLED_SEPARATION_FLOOR
