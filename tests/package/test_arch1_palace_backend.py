"""`ARCH-1` step 1a (§9 item 517, U16a): the Palace backend's write and read halves.

``Study.solve(backend="palace")`` / ``femem run --backend palace`` write a
Palace case that equals `OPS-61`'s ×1 exact-strip ``f05x_o2`` / ``f05x`` inputs
(config keys; mesh byte for byte), and read S back through the lifted reader on
`OPS-61`'s recorded output ((Rd) ``np.array_equal`` against the probe's
``_read_port_s``). T1 and the CLI write stage write to ``tmp_path`` (item
533), so they stay repeatable once the on-disk cases carry Palace output.
T4 / T5 (item 533, U16c; skipped unless ``FEM_EM_U16_PALACE_READ=1``) read S
back through ``Study.solve`` / ``femem run`` on the Palace output of item
517b's runs and compare it with `OPS-61`'s recorded read on the same input
((E2) / (E1) at ``LEG_D0_REPRODUCTION_BAND``), (P), (Tie2) and (c-order).
Brief: docs/planning/chunks/ARCH-1.md.
"""

import importlib.util
import os
from pathlib import Path

import numpy as np
import pytest

from fem_em_solver.backends.palace import read_port_s

ROOT = Path(__file__).resolve().parents[2]
OPS61 = ROOT / "logs" / "ops61"
PROBE = ROOT / "scripts" / "probes" / "ops61_step1a_palace_ans4.py"


def _load_probe():
    spec = importlib.util.spec_from_file_location("ops61_probe", PROBE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


U16_DIR = ROOT / "logs" / "u16"
D2_CASE = U16_DIR / "f_small_d2_64"


def _flat(d, prefix=""):
    out = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flat(v, key + "."))
        else:
            out[key] = v
    return out


def test_u16_case_equals_ops61(tmp_path):
    import hashlib
    import json
    import time

    from mpi4py import MPI

    from fem_em_solver.backends.palace import PalaceRunRequired
    from fem_em_solver.materials import GelledSalinePhantomMaterial
    from fem_em_solver.study import FourPortBirdcageGeometry, LumpedPortDrive, Study

    ref_cfg = OPS61 / "f05x_o2.json"
    ref_msh = OPS61 / "ops61_ans4_x1_strip.msh"
    for p in (ref_cfg, ref_msh):
        if not p.exists():
            pytest.skip(f"BLOCKED: {p} absent")
    study = Study(
        FourPortBirdcageGeometry.f_small(),
        GelledSalinePhantomMaterial(sigma=0.5, epsilon_r=78.0, frequency_hz=64e6),
        LumpedPortDrive(),
        frequency_hz=64e6,
        degree=2,
    )
    case = tmp_path / "d2"
    t0 = time.perf_counter()
    with pytest.raises(PalaceRunRequired) as exc:
        study.solve(MPI.COMM_WORLD, backend="palace", case_dir=case)
    print(f"[U16a T1] write stage {time.perf_counter() - t0:.1f} s; {exc.value}", flush=True)
    assert exc.value.argv == ["palace", "-np", "8", "case.json"]
    assert not (case / "postpro").exists()
    ours = _flat(json.loads((case / "case.json").read_text()))
    ref = _flat(json.loads(ref_cfg.read_text()))
    diff = {k for k in set(ours) | set(ref) if ours.get(k) != ref.get(k)}
    print(f"[U16a W-cfg] differing keys: {sorted(diff)}", flush=True)
    assert diff == {"Problem.Output", "Model.Mesh"}
    a = (case / "mesh.msh").read_bytes()
    b = ref_msh.read_bytes()
    print(
        f"[U16a W-msh] sha256 ours {hashlib.sha256(a).hexdigest()} "
        f"ops61 {hashlib.sha256(b).hexdigest()}",
        flush=True,
    )
    ntet = sum(1 for ln in a.decode().splitlines() if ln.split()[1:2] == ["4"])
    print(f"[U16a W-msh] tets {ntet}", flush=True)
    assert ntet == 117206
    assert a == b


CLI_OUT = U16_DIR / "f_small_d1_64_cli"
CLI_CASE = CLI_OUT / "palace_case"
YAML = ROOT / "tests" / "package" / "wf1" / "f_small_four_port.yaml"
READ_ENV = "FEM_EM_U16_PALACE_READ"


def _d2_study():
    from fem_em_solver.materials import GelledSalinePhantomMaterial
    from fem_em_solver.study import FourPortBirdcageGeometry, LumpedPortDrive, Study

    return Study(
        FourPortBirdcageGeometry.f_small(),
        GelledSalinePhantomMaterial(sigma=0.5, epsilon_r=78.0, frequency_hz=64e6),
        LumpedPortDrive(),
        frequency_hz=64e6,
        degree=2,
    )


def test_u16_backend_refusals(tmp_path):
    import time

    from mpi4py import MPI

    from fem_em_solver.study import FourPortBirdcageGeometry, LumpedPortDrive

    cases = []
    s = _d2_study()
    s.geometry = FourPortBirdcageGeometry(as_hole=True)
    cases.append(("as_hole", s, ValueError))
    s = _d2_study()
    s.conductor = object()
    cases.append(("conductor", s, ValueError))
    s = _d2_study()
    s.excitation = LumpedPortDrive(z0_ohm=75.0)
    cases.append(("z0 != R", s, ValueError))
    s = _d2_study()
    s.excitation = object()
    cases.append(("excitation", s, TypeError))
    s = _d2_study()
    s.geometry = object()
    cases.append(("geometry", s, TypeError))
    for name, study, err in cases:
        t0 = time.perf_counter()
        with pytest.raises(err):
            study.solve(MPI.COMM_WORLD, backend="palace", case_dir=tmp_path / "c")
        dt = time.perf_counter() - t0
        print(f"[U16a T0] {name}: {err.__name__} in {dt:.3f} s", flush=True)
        assert dt < 2.0
    with pytest.raises(ValueError):
        _d2_study().solve(MPI.COMM_WORLD, backend="hfss")

    class _Two:
        size = 2
        rank = 0

    with pytest.raises(ValueError, match="serial"):
        _d2_study().solve(_Two(), backend="palace", case_dir=tmp_path / "c")
    assert not any(tmp_path.iterdir())


def _p_rec(s):
    from tests.validation.test_port_birdcage_four_port import (
        PASSIVITY_SIGMA_TOLERANCE,
        POOLED_SEPARATION_FLOOR,
        _circulant_classes,
        _class_spread,
    )
    from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
    from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

    recip = float(np.max(np.abs(s - s.T)) / np.max(np.abs(s)))
    sig = float(np.linalg.svd(s, compute_uv=False)[0])
    colp = np.sum(np.abs(s) ** 2, axis=0)
    cls = _circulant_classes(s)
    spreads = {k: _class_spread(v) for k, v in cls.items()}
    pooled = _class_spread(np.concatenate([cls["adjacent"], cls["opposite"]]))
    worst = max(spreads.values())
    print(
        f"[U16a P-rec] reciprocity {recip:.3e} (band {RECIPROCITY_BAND}); "
        f"sigma_max {sig:.12f} (item 194: 0.999819887373 on f05x_o2); "
        f"col powers {np.round(colp, 9)}; spreads "
        f"{ {k: f'{v * 100:.4f}%' for k, v in spreads.items()} }; pooled "
        f"{pooled * 100:.4f}% = {pooled / worst:.1f}x worst",
        flush=True,
    )
    assert recip <= RECIPROCITY_BAND
    assert sig <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    assert np.all(colp <= 1.0 + PASSIVITY_SIGMA_TOLERANCE)
    assert all(v <= ADJACENT_SPREAD_BAND for v in spreads.values())
    assert pooled >= POOLED_SEPARATION_FLOOR * worst


def _plant(case_src, postpro_src, dst):
    import shutil

    for name in ("case.json", "mesh.msh", "case.sha256"):
        shutil.copy2(case_src / name, dst / name)
    shutil.copytree(postpro_src, dst / "postpro")


def test_u16_cli_write_stage_and_controls(tmp_path):
    import json

    from fem_em_solver import cli

    for p in (
        YAML,
        OPS61 / "f05x.json",
        OPS61 / "ops61_ans4_x1.msh",
        OPS61 / "ops61_ans4_x1_strip.msh",
    ):
        if not p.exists():
            pytest.skip(f"BLOCKED: {p} absent")
    out = tmp_path / "cli"
    case = out / "palace_case"
    rc = cli.main(["run", str(YAML), "--output-dir", str(out), "--backend", "palace"])
    assert rc == 3
    assert not (case / "postpro").exists()
    assert not (out / "summary.json").exists()
    ours = _flat(json.loads((case / "case.json").read_text()))
    ref = _flat(json.loads((OPS61 / "f05x.json").read_text()))
    diff = {k for k in set(ours) | set(ref) if ours.get(k) != ref.get(k)}
    print(f"[U16a W-cfg cli] vs f05x.json differing keys: {sorted(diff)}", flush=True)
    assert diff == {"Problem.Output", "Model.Mesh"}
    # (c1) the order mapping reaches the config
    ref2 = _flat(json.loads((OPS61 / "f05x_o2.json").read_text()))
    diff2 = {k for k in set(ours) | set(ref2) if ours.get(k) != ref2.get(k)}
    print(f"[U16a c1] vs f05x_o2.json differing keys: {sorted(diff2)}", flush=True)
    assert diff2 == {"Problem.Output", "Model.Mesh", "Solver.Order"}
    a = (case / "mesh.msh").read_bytes()
    assert a == (OPS61 / "ops61_ans4_x1_strip.msh").read_bytes()
    # (c2) the ragged control differs
    rag = (OPS61 / "ops61_ans4_x1.msh").read_text().splitlines()
    ntet_rag = sum(1 for ln in rag if ln.split()[1:2] == ["4"])
    print(f"[U16a c2] ragged control tets {ntet_rag} vs 117206", flush=True)
    assert ntet_rag == 116085 and a != (OPS61 / "ops61_ans4_x1.msh").read_bytes()


def test_u16_read_path_on_recorded_output(tmp_path):
    from mpi4py import MPI

    from fem_em_solver.backends.palace import PalaceStudyResult

    src = OPS61 / "postpro_f05x_o2"
    if not (D2_CASE / "case.json").exists() or not src.exists():
        pytest.skip("BLOCKED: D2_CASE or postpro_f05x_o2 absent")
    case = tmp_path / "case"
    case.mkdir()
    _plant(D2_CASE, src, case)
    res = _d2_study().solve(MPI.COMM_WORLD, backend="palace", case_dir=case)
    assert isinstance(res, PalaceStudyResult)
    s_ref, _ = read_port_s(src / "port-S.csv")
    assert np.array_equal(res.s_matrix, s_ref)
    print(
        f"[U16a T2] read back; palace_total_s {res.palace_total_s}; cells " f"{res.cells}",
        flush=True,
    )
    assert res.cells == 117206 and res.order == 2
    _p_rec(res.s_matrix)
    # (c3) a 10 MHz output is refused by the frequency check
    case3 = tmp_path / "case3"
    case3.mkdir()
    _plant(D2_CASE, OPS61 / "postpro_f05x_o2_f10", case3)
    with pytest.raises(ValueError, match="0.01") as e3:
        _d2_study().solve(MPI.COMM_WORLD, backend="palace", case_dir=case3)
    assert "0.064" in str(e3.value)
    print(f"[U16a c3] {e3.value}", flush=True)
    # (c4) a one-byte edit to case.json is caught by the byte check
    cj = case / "case.json"
    raw = bytearray(cj.read_bytes())
    raw[-2] = ord(" ") if raw[-2] != ord(" ") else ord("\t")
    cj.write_bytes(bytes(raw))
    with pytest.raises(FileExistsError, match="case.json") as e4:
        _d2_study().solve(MPI.COMM_WORLD, backend="palace", case_dir=case)
    print(f"[U16a c4] {e4.value}", flush=True)


def test_u16_cli_read_path_on_recorded_output(tmp_path):
    import json
    import shutil

    from fem_em_solver import cli

    src = OPS61 / "postpro_f05x"
    if not (CLI_CASE / "case.json").exists() or not src.exists():
        pytest.skip("BLOCKED: CLI_CASE or postpro_f05x absent")
    out = tmp_path / "cli"
    (out / "palace_case").mkdir(parents=True)
    _plant(CLI_CASE, src, out / "palace_case")
    rc = cli.main(["run", str(YAML), "--output-dir", str(out), "--backend", "palace"])
    assert rc == 0
    summary = json.loads((out / "summary.json").read_text())
    assert summary["backend"] == "palace"
    assert "z_matrix" not in summary and "touchstone" not in summary
    s_ref, _ = read_port_s(src / "port-S.csv")
    got = np.array([[complex(*z) for z in row] for row in summary["s_matrix"]])
    assert np.array_equal(got, s_ref)
    print(f"[U16a T3] summary palace block {summary['palace']}", flush=True)
    shutil.rmtree(out)


# (Tie2) ceiling: item 476's ratified public ceiling, asserted per the
# 2026-10-18 03:00 review's R3 (a) ruling (ARCH-1.md), never the weekly's 1.9e-10.
TIE2_BAND = 2.6e-10
# Recorded per-class (Tie2) values, printed beside
# (`20261004T214342Z_OPS-61-item229-W4-read.log:41`).
TIE2_RECORDED = (2.464e-10, 9.668e-11, 9.883e-11)
C_ORDER_FLOOR = 1.0e-6


def _class_means(s):
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    return {k: np.mean(v) for k, v in _circulant_classes(s).items()}


def _fro_rel(s, ref):
    return float(np.linalg.norm(s - ref) / np.linalg.norm(ref))


@pytest.mark.skipif(os.environ.get(READ_ENV) != "1", reason=f"item 517b: {READ_ENV}=1")
def test_u16_palace_read_back():
    from mpi4py import MPI

    from fem_em_solver.backends.palace import PalaceStudyResult
    from tests.validation.test_ans4_exact_strip_records import ANS4_EXACT_STRIP_S_P2_X1
    from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

    for p in (
        D2_CASE / "case.json",
        D2_CASE / "postpro" / "port-S.csv",
        OPS61 / "postpro_f05x_o2" / "port-S.csv",
        OPS61 / "postpro_f05x" / "port-S.csv",
        OPS61 / "postpro_f05x_o2_f10" / "port-S.csv",
    ):
        if not p.exists():
            pytest.skip(f"BLOCKED: {p} absent (comparand absent)")
    res = _d2_study().solve(MPI.COMM_WORLD, backend="palace", case_dir=D2_CASE)
    assert isinstance(res, PalaceStudyResult)
    assert res.cells == 117206 and res.order == 2, (res.cells, res.order)
    s = np.asarray(res.s_matrix)
    s_ops61, _ = read_port_s(OPS61 / "postpro_f05x_o2" / "port-S.csv")
    e2 = _fro_rel(s, s_ops61)
    print(
        f"[U16c T4 E2] Study read vs postpro_f05x_o2: Frobenius-rel {e2:.3e} (band "
        f"{LEG_D0_REPRODUCTION_BAND:.0e}); entry max {np.max(np.abs(s - s_ops61)):.3e}; "
        f"palace_total_s {res.palace_total_s}",
        flush=True,
    )
    _p_rec(s)
    m = _class_means(s)
    ref = _class_means(ANS4_EXACT_STRIP_S_P2_X1["64"])
    tie2 = {k: float(abs(m[k] - ref[k]) / abs(ref[k])) for k in ref}
    print(
        "[U16c T4 Tie2] "
        + ", ".join(f"{k} {v:.3e}" for k, v in tie2.items())
        + f" (band {TIE2_BAND:.1e}; recorded {TIE2_RECORDED})",
        flush=True,
    )
    o1 = _class_means(read_port_s(OPS61 / "postpro_f05x" / "port-S.csv")[0])
    c_order = {k: float(abs(m[k] - o1[k]) / abs(o1[k])) for k in o1}
    print(
        "[U16c T4 c-order] "
        + ", ".join(f"{k} {v * 100:.4f}%" for k, v in c_order.items())
        + f" (> {C_ORDER_FLOOR:.0e} each)",
        flush=True,
    )
    s10 = read_port_s(OPS61 / "postpro_f05x_o2_f10" / "port-S.csv")[0]
    print(
        f"[U16c T4 c-freq] Frobenius-rel vs 10 MHz o2 {_fro_rel(s, s10):.3e} (predicted "
        f"O(1e-2 .. 1e-1); printed only)",
        flush=True,
    )
    assert e2 <= LEG_D0_REPRODUCTION_BAND, e2
    for k, v in tie2.items():
        assert v <= TIE2_BAND, (k, v)
    for k, v in c_order.items():
        assert v > C_ORDER_FLOOR, (k, v)


@pytest.mark.skipif(os.environ.get(READ_ENV) != "1", reason=f"item 517b: {READ_ENV}=1")
def test_u16_cli_palace_read_back():
    import json

    from fem_em_solver import cli
    from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

    ref_csv = OPS61 / "postpro_f05x" / "port-S.csv"
    for p in (YAML, CLI_CASE / "case.json", CLI_CASE / "postpro" / "port-S.csv", ref_csv):
        if not p.exists():
            pytest.skip(f"BLOCKED: {p} absent (comparand absent)")
    stale = CLI_OUT / "summary.json"
    if stale.exists():
        stale.unlink()
    rc = cli.main(["run", str(YAML), "--output-dir", str(CLI_OUT), "--backend", "palace"])
    assert rc == 0
    summary = json.loads((CLI_OUT / "summary.json").read_text())
    assert summary["backend"] == "palace"
    got = np.array([[complex(*z) for z in row] for row in summary["s_matrix"]])
    s_ref, _ = read_port_s(ref_csv)
    e1 = _fro_rel(got, s_ref)
    print(
        f"[U16c T5 E1] CLI summary vs postpro_f05x: Frobenius-rel {e1:.3e} (band "
        f"{LEG_D0_REPRODUCTION_BAND:.0e}); entry max {np.max(np.abs(got - s_ref)):.3e}; "
        f"palace block {summary.get('palace')}",
        flush=True,
    )
    assert e1 <= LEG_D0_REPRODUCTION_BAND, e1
    _p_rec(got)


@pytest.mark.parametrize(
    "run, freq_ghz",
    [("postpro_f05x_o2", 0.064), ("postpro_f05x", 0.064), ("postpro_f05x_o2_f10", 0.01)],
)
def test_u16_reader_equals_probe(run, freq_ghz):
    path = OPS61 / run / "port-S.csv"
    if not path.exists():
        pytest.skip(f"BLOCKED: {path} absent")
    probe_s, _header = _load_probe()._read_port_s(path)
    s, f = read_port_s(path)
    print(
        f"[U16a Rd] {run}: shape {s.shape}, f = {f} GHz, " f"max|S| {np.max(np.abs(s)):.12f}",
        flush=True,
    )
    assert s.shape == (4, 4)
    assert np.array_equal(s, probe_s)
    assert f == freq_ghz
