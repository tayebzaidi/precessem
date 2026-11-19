"""`WF-1` step 5 (§9 item 344, chain U33): ``outputs.sar`` / ``outputs.b1_plus``.

Schema v0's ``outputs`` gains ``sar: {mass_g, centres}`` and ``b1_plus:
{points}`` on a port Study; ``run_study_file`` drives ``Study.sar`` (every drive
at every centre) and ``Study.b1_plus_at`` after the S block and writes both into
the JSON summary (brief: docs/planning/chunks/WF-1.md, "§9 item 344").  The
``{ring_radius_m: r}`` centre form is an API convenience (ruling 2026-10-08
(3)): one centre per port at the solved sheets' azimuths.  Two solves
(``Study.solve`` + the cached field-keeping solve) are disclosed (ruling (4)).

Asserted, on item 192's fixture at degree 1, 64 MHz, ``-n 2``:
(K)   spec diff vs ``_port_study(degree=1)`` is ``[]``; ``mass_g`` and
      ``ring_radius_m`` by ``float.hex``; the summary centres equal the
      sheet-derived centres (`test_wf16_study_sar.py`'s expression) by ``float.hex``;
(R)   ``summary["s_matrix"]`` within ``LEG_D0_REPRODUCTION_BAND`` of ``_record(1)``;
(ii)  the 10 g C4 identity read from the summary, each pair <= ``C4_COVARIANCE_BAND``;
(iii) the summary SAR == ``Study.sar`` per drive at ``EXACT_IDENTITY_RTOL``; JSON == memory;
(B1)  ``|B₁⁺|`` C4 covariance on the YAML-solved Study; (B3) every point valid;
(B2′) the summary ``b1_plus`` == ``Study.b1_plus_at`` at ``EXACT_IDENTITY_RTOL``;
(X)   two malformed files raise ``StudyConfigError`` naming ``outputs.sar`` /
      ``outputs.sar.mass_g``; ``cli.main`` run / check rc 2 on every rank.
Negative controls (asserted): (c1) the mis-paired far-side ball and (c2) the
mis-rotated ``|B₁⁺|`` image miss ``C4_COVARIANCE_BAND``.

Every band, record and helper is imported, none restated.
"""

from __future__ import annotations

import json
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver import StudyConfigError, cli, load_study, run_study_file
from tests.complex_mode import complex_only
from tests.package.test_wf1_study_yaml import LOGS, WF1, _fresh_dir, _rc_all, _spec_diff
from tests.package.test_wf16_study_four_port import _record, _rel
from tests.package.test_wf16_study_sar import _port_study
from tests.validation.test_ans4_exact_strip_records import (
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
)
from tests.validation.test_birdcage_b1_plus_map import (
    C4_COVARIANCE_BAND,
    _relative_l2,
    _rotate_z,
    _sample_points,
)
from tests.validation.test_birdcage_b1_quadrature import QUADRATURE_STEP_DEG, _port_index
from tests.validation.test_birdcage_sar_mass_averaged import (
    CENTRE_RADIUS_M,
    EXACT_IDENTITY_RTOL,
)
from tests.validation.test_mass_averaged_sar_standard_masses import TEN_GRAM_KG
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

SAR = "averaged_sar_w_per_kg"
SAR_B1_FILE = WF1 / "f_small_sar_b1.yaml"
MALFORMED_SAR = {
    "malformed_sar_on_sphere.yaml": "outputs.sar",
    "malformed_sar_mass_zero.yaml": "outputs.sar.mass_g",
}


def _p0(*args):
    if MPI.COMM_WORLD.rank == 0:
        print(*args, flush=True)


def _hex3(c):
    return [float(x).hex() for x in c]


def test_wf1_sar_b1_yaml_parse_and_refusals():
    """(K) parse half and (X): no mesh, no solve."""
    comm = MPI.COMM_WORLD
    from fem_em_solver.study_yaml import _load

    study, kind, outputs, _ = _load(SAR_B1_FILE)
    diff = _spec_diff(study, _port_study(degree=1))
    sar = outputs["sar"]
    _p0(
        f"\n[WF-1 s5] (K) spec diff {diff}; mass_g hex {sar['mass_g'].hex()} vs "
        f"1000*TEN_GRAM_KG {(1000.0 * TEN_GRAM_KG).hex()}; ring_radius_m hex "
        f"{sar['centres']['ring_radius_m'].hex()} vs CENTRE_RADIUS_M {CENTRE_RADIUS_M.hex()}; "
        f"b1 points {outputs['b1_plus']['points']}"
    )
    assert kind == "four_port_birdcage"
    assert diff == [], diff
    assert sar["mass_g"].hex() == (10.0).hex() == (1000.0 * TEN_GRAM_KG).hex()
    assert sar["centres"]["ring_radius_m"].hex() == float(CENTRE_RADIUS_M).hex()
    assert outputs["b1_plus"]["points"] == [(0.0, 0.0, 0.0)]

    out = LOGS / "malformed_sar"
    fields = []
    for name, expected in MALFORMED_SAR.items():
        path = WF1 / name
        with pytest.raises(StudyConfigError) as ei:
            load_study(path)
        exc = ei.value
        _p0(f"  (X) {name}: field {exc.field!r} (expected {expected!r}); {exc}")
        assert type(exc) is StudyConfigError
        assert exc.field == expected
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
    assert fields == list(MALFORMED_SAR.values()), fields


@pytest.fixture(scope="module")
def case():
    comm = MPI.COMM_WORLD
    assert (
        comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
    ), f"degree 1 record is asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}; got -n {comm.size}"
    d = LOGS / SAR_B1_FILE.stem
    _fresh_dir(d, comm)
    t0 = time.perf_counter()
    run = run_study_file(SAR_B1_FILE, output_dir=d, comm=comm)
    t_run = time.perf_counter() - t0
    with open(d / "summary.json") as fh:
        summary = json.load(fh)
    study = run.study
    res = study.fields_result
    iface = study.geometry.sheet_iface
    order = [f"P{s['tag'] - iface}" for s in res.sheets]
    az = {f"P{s['tag'] - iface}": float(s["azimuth_deg"]) for s in res.sheets}
    idx = {pid: _port_index(az[pid], az["P1"], QUADRATURE_STEP_DEG) for pid in az}
    by_k = {k: pid for pid, k in idx.items()}
    assert sorted(by_k) == [0, 1, 2, 3]
    # The test's own centres: test_wf16_study_sar.py's expression, per sheet.
    own = {
        pid: (
            CENTRE_RADIUS_M * float(np.cos(np.radians(az[pid]))),
            CENTRE_RADIUS_M * float(np.sin(np.radians(az[pid]))),
            0.0,
        )
        for pid in order
    }

    t0 = time.perf_counter()
    api_sar = {
        pid: study.sar(mass_g=10.0, centres=[own[pid]], drives=[pid], comm=comm)[pid][0]
        for pid in order
    }
    t_sar = time.perf_counter() - t0

    delta = np.radians((az["P2"] - az["P1"]) % 360.0)
    points = _sample_points({"mesh": res.mesh, "cell_tags": res.cell_tags})
    images = {
        "P1@0": ("P1", points),
        "P2@+D": ("P2", _rotate_z(points, delta)),
        "P4@-D": ("P4", _rotate_z(points, -delta)),
        "P3@2D": ("P3", _rotate_z(points, 2.0 * delta)),
        "P2@-D (control)": ("P2", _rotate_z(points, -delta)),
    }
    t0 = time.perf_counter()
    read = {
        label: study.b1_plus_at(pts, drives=[pid], comm=comm)[pid]
        for label, (pid, pts) in images.items()
    }
    api_b1 = study.b1_plus_at(summary["b1_plus"]["points"], comm=comm)
    t_b1 = time.perf_counter() - t0
    _p0(
        f"\n[WF-1 s5] run_study_file {t_run:.1f} s (summary elapsed "
        f"{summary['elapsed_s']:.1f} s); {res.cells} cells; API sar 4 calls {t_sar:.2f} s; "
        f"b1 evaluations {t_b1:.2f} s; sheet order {order}; Δ {np.degrees(delta):.6f} deg; "
        f"{points.shape[0]} sample points; ball radius {summary['sar']['radius_m']:.9e} m"
    )
    return {
        "run": run,
        "summary": summary,
        "order": order,
        "by_k": by_k,
        "own": own,
        "api_sar": api_sar,
        "read": read,
        "api_b1": api_b1,
        "points": points,
    }


def _diag(case, k, centre_of=None):
    """SAR(c_{centre_of(k)}; P_k) from the summary (default: the diagonal)."""
    by_k, order, vals = case["by_k"], case["order"], case["summary"]["sar"]["values"]
    pid = by_k[k]
    cpid = by_k[k if centre_of is None else centre_of(k)]
    return vals[pid][order.index(cpid)][SAR]


@complex_only
def test_wf1_sar_b1_yaml_k_centres_bitwise(case):
    """(K) after the run: summary centres == the sheet-derived centres by float.hex."""
    got = case["summary"]["sar"]["centres"]
    want = [case["own"][pid] for pid in case["order"]]
    _p0(f"\n[WF-1 s5] (K) summary centres {got}")
    assert [_hex3(c) for c in got] == [_hex3(c) for c in want]
    assert case["summary"]["sar"]["mass_g"] == 10.0
    assert set(case["summary"]["sar"]["values"]) == set(case["order"])
    for pid, per in case["summary"]["sar"]["values"].items():
        assert len(per) == 4, (pid, len(per))


@complex_only
def test_wf1_sar_b1_yaml_r_same_fixture(case):
    """(R) the YAML S matrix vs item 192's degree-1 record."""
    s = np.array([[complex(re, im) for re, im in row] for row in case["summary"]["s_matrix"]])
    rel = _rel(s, _record(1))
    _p0(
        f"\n[WF-1 s5] (R) ||S_yaml - P1|| / ||P1|| = {rel:.3e} "
        f"(ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e})"
    )
    assert rel <= LEG_D0_REPRODUCTION_BAND


@complex_only
def test_wf1_sar_b1_yaml_ii_ten_gram_c4(case):
    """(ii) the 10 g C4 identity, read from the summary."""
    d = {k: _diag(case, k) for k in range(4)}
    pairs = {k: abs(d[(k + 1) % 4] - d[k]) / d[k] for k in range(4)}
    _p0(
        "\n[WF-1 s5] (ii) diagonal (printed, per 1 V): "
        + ", ".join(f"{case['by_k'][k]} {d[k]:.9e} W/kg" for k in range(4))
        + "\n    cyclic pairs: "
        + ", ".join(f"{k}->{(k + 1) % 4} {pairs[k] * 100:.4f}%" for k in range(4))
        + f" (ASSERTED <= {C4_COVARIANCE_BAND * 100:.1f}%)"
    )
    for k in range(4):
        assert pairs[k] <= C4_COVARIANCE_BAND, (k, pairs[k])


@complex_only
def test_wf1_sar_b1_yaml_iii_yaml_equals_api(case):
    """(iii) summary SAR == Study.sar per drive; JSON == in-memory."""
    vals = case["summary"]["sar"]["values"]
    worst = 0.0
    for pid in case["order"]:
        got = vals[pid][case["order"].index(pid)]
        api = case["api_sar"][pid]
        assert set(got) == set(api)
        for key, value in api.items():
            if value != 0.0:
                worst = max(worst, abs(got[key] - value) / abs(value))
    _p0(
        f"\n[WF-1 s5] (iii) summary vs Study.sar: max rel {worst:.3e} "
        f"(ASSERTED <= {EXACT_IDENTITY_RTOL:.0e})"
    )
    assert worst <= EXACT_IDENTITY_RTOL
    mem = case["run"].summary
    assert case["summary"]["sar"] == mem["sar"]
    assert case["summary"]["b1_plus"] == mem["b1_plus"]


@complex_only
def test_wf1_sar_b1_yaml_b1_covariance(case):
    """(B1) |B1+| C4 covariance on the YAML-solved Study; (B3) every point valid."""
    read, n = case["read"], case["points"].shape[0]
    counts = {label: int(v.sum()) for label, (_, v) in read.items()}
    mask = np.logical_and.reduce([v for _, v in read.values()])
    ref = read["P1@0"][0]
    table = {
        label: _relative_l2(vals, ref, mask)
        for label, (vals, _) in read.items()
        if label not in ("P1@0", "P2@-D (control)")
    }
    _p0(
        f"\n[WF-1 s5] (B3) valid per image (of {n}): {counts}\n[WF-1 s5] (B1) "
        + ", ".join(f"{k} {v * 100:.4f}%" for k, v in table.items())
        + f" (ASSERTED <= {C4_COVARIANCE_BAND * 100:.1f}%)"
    )
    assert n > 0
    for label, c in counts.items():
        assert c == n, (label, c, n)
    for label, value in table.items():
        assert value <= C4_COVARIANCE_BAND, (label, value)


@complex_only
def test_wf1_sar_b1_yaml_b2_prime(case):
    """(B2′) summary b1_plus == Study.b1_plus_at on the written points; all valid."""
    b1 = case["summary"]["b1_plus"]
    worst = 0.0
    for pid, (vals, ok) in case["api_b1"].items():
        got = np.asarray(b1["values"][pid])
        assert all(b1["valid"][pid]) and bool(np.all(ok))
        worst = max(worst, float(np.max(np.abs(got - vals)) / np.max(np.abs(vals))))
    _p0(
        f"\n[WF-1 s5] (B2') summary |B1+| at {b1['points']} (printed, per 1 V): "
        + ", ".join(f"{p} {v[0]:.9e} T" for p, v in b1["values"].items())
        + f"; max|Δ| / max|B1+| = {worst:.3e} (ASSERTED <= {EXACT_IDENTITY_RTOL:.0e})"
    )
    assert set(b1["values"]) == set(case["order"])
    assert worst <= EXACT_IDENTITY_RTOL


@complex_only
def test_wf1_sar_b1_yaml_c1_mispaired(case):
    """(c1) the far-side ball from the summary misses the band at every k."""
    miss = {}
    for k in range(4):
        diag = _diag(case, k)
        far = _diag(case, k, centre_of=lambda j: (j + 2) % 4)
        miss[k] = abs(far - diag) / diag
    _p0(
        "\n[WF-1 s5] (c1) |SAR(c_k+2;k) - SAR(c_k;k)| / SAR(c_k;k): "
        + ", ".join(f"k={k} {miss[k] * 100:.4f}%" for k in range(4))
        + f" (sign ASSERTED > {C4_COVARIANCE_BAND * 100:.1f}%; item 321 read 88.12 / 87.94 / "
        "88.66 / 88.93%)"
    )
    for k in range(4):
        assert miss[k] > C4_COVARIANCE_BAND, (k, miss[k])


@complex_only
def test_wf1_sar_b1_yaml_c2_misrotation(case):
    """(c2) P2 at the −Δ points vs P1 misses the band."""
    read = case["read"]
    mask = np.logical_and.reduce([v for _, v in read.values()])
    miss = _relative_l2(read["P2@-D (control)"][0], read["P1@0"][0], mask)
    _p0(
        f"\n[WF-1 s5] (c2) P2@-D vs P1: {miss * 100:.4f}% (sign ASSERTED > "
        f"{C4_COVARIANCE_BAND * 100:.1f}%; item 322 read 35.0040%)"
    )
    assert miss > C4_COVARIANCE_BAND
