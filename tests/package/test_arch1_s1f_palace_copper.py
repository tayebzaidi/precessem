"""`ARCH-1` step 1f (§9 item 631, §10 U91): the copper F-small x1 hole on the Palace backend.

``Study.solve(backend="palace")`` accepts ``FourPortBirdcageGeometry(as_hole=True)``
together with ``conductor=LeontovichWall(sigma)`` -- only that combination -- and
writes the hole-wall facet attribute plus a Palace ``Boundaries.Conductivity``
block, ported from ``scripts/probes/ops61_step1a_palace_ans4.py``
(``export_hole_pec`` :1106, ``_hole_wall_attr`` :785, the Conductivity block :1533).

Anchors: (E2-Cu) the Study's read of its own Palace run against the probe's
``logs/ops61/postpro_h05x_cu_o2/port-S.csv`` (`OPS-61` step 2b, §9 item 241),
entry by entry at ``LEG_D0_REPRODUCTION_BAND``; (Tie-Cu) its class means against
``ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["cu"]``; the PEC-hole control (the probe's
``postpro_h05x_pec_o2``, same mesh and order, wall PEC instead of Conductivity)
misses by more than ``NEGATIVE_CONTROL_FLOOR`` (asserted). Palace v0.18.1.
"""

import json
import os
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
OPS61 = ROOT / "logs" / "ops61"
CU_CASE = ROOT / "logs" / "arch1" / "s1f" / "cu_case"
PROBE_CFG = OPS61 / "h05x_cu_o2.json"
PROBE_MESH = OPS61 / "ops61_ans6_x1_hole_strip.msh"
PROBE_CU = OPS61 / "postpro_h05x_cu_o2" / "port-S.csv"
PROBE_PEC = OPS61 / "postpro_h05x_pec_o2" / "port-S.csv"
READ_ENV = "FEM_EM_S1F_PALACE_READ"
# WF-1's degree-1 PEC-vs-copper class miss, printed beside the degree-2 one
# (`20261020T202415Z_WF-1-s21-RGCu-copper.log:2716`).
WF1_D1_MISS = (6.422e-04, 2.812e-04, 4.863e-04)
# Step 2b's (T4) read, printed beside (`20261005T093623Z_OPS-61-item241-W6-read.log.gz:56`).
T4_RECORDED = (1.612e-10, 1.586e-10, 2.033e-10)


def _cu_study(geometry=None, conductor="cu"):
    from fem_em_solver.materials import GelledSalinePhantomMaterial
    from fem_em_solver.study import (
        FourPortBirdcageGeometry,
        LeontovichWall,
        LumpedPortDrive,
        Study,
    )

    s = Study(
        FourPortBirdcageGeometry(as_hole=True) if geometry is None else geometry,
        GelledSalinePhantomMaterial(sigma=0.5, epsilon_r=78.0, frequency_hz=64e6),
        LumpedPortDrive(),
        frequency_hz=64e6,
        degree=2,
    )
    s.conductor = LeontovichWall(5.8e7) if conductor == "cu" else conductor
    return s


def _census(msh_path):
    """(nodes, Counter{(gmsh type, attr): count}) of a Gmsh 2.2 ASCII mesh."""
    counts = Counter()
    nodes = None
    with open(msh_path) as fh:
        mode = None
        for line in fh:
            t = line.strip()
            if t in ("$Nodes", "$Elements"):
                mode, first = t, True
                continue
            if t.startswith("$End"):
                mode = None
                continue
            if mode == "$Nodes" and first:
                nodes, first = int(t), False
            elif mode == "$Elements":
                if first:
                    first = False
                    continue
                p = t.split()
                counts[(int(p[1]), int(p[3]))] += 1
    return nodes, counts


def _class_means(s):
    from tests.validation.test_port_birdcage_four_port import _circulant_classes

    return {k: np.mean(v) for k, v in _circulant_classes(s).items()}


def test_s1f_refusals_and_acceptance(tmp_path, monkeypatch):
    """W0: every other combination still refused; as_hole + Leontovich reaches the build."""
    from mpi4py import MPI

    from fem_em_solver.backends.palace import palace_config
    from fem_em_solver.study import FourPortBirdcageGeometry

    cases = [
        ("as_hole, no conductor", _cu_study(conductor=None)),
        ("Leontovich on a solid geometry", _cu_study(geometry=FourPortBirdcageGeometry())),
        ("non-Leontovich conductor on the hole", _cu_study(conductor=object())),
    ]
    for name, study in cases:
        with pytest.raises(ValueError) as ei:
            study.solve(MPI.COMM_WORLD, backend="palace", case_dir=tmp_path / "c")
        print(f"[s1f W0] {name}: ValueError({ei.value})", flush=True)
    assert not any(tmp_path.iterdir())

    class _Reached(Exception):
        pass

    def _no_build(self, comm=None):
        raise _Reached

    monkeypatch.setattr(FourPortBirdcageGeometry, "build", _no_build)
    with pytest.raises(_Reached):
        _cu_study().solve(MPI.COMM_WORLD, backend="palace", case_dir=tmp_path / "c")
    print("[s1f W0] as_hole + LeontovichWall(5.8e7): passes every refusal, reaches build",
          flush=True)

    # The Boundaries block equals the probe's own h05x_cu_o2.json (when on disk).
    cfg = palace_config(
        mesh_name="mesh.msh", output="postpro", materials=[], ports=[1, 2, 3, 4],
        order=2, freq_ghz=0.064, r_ohm=50.0,
        conductivity=[{"Attributes": [401], "Conductivity": 5.8e7, "Permeability": 1.0}],
    )
    if PROBE_CFG.exists():
        probe = json.loads(PROBE_CFG.read_text())
        assert cfg["Boundaries"] == probe["Boundaries"], cfg["Boundaries"]
        assert cfg["Solver"] == probe["Solver"]
        print("[s1f W0] Boundaries / Solver == h05x_cu_o2.json", flush=True)
    assert "Conductivity" not in palace_config(
        mesh_name="m", output="o", materials=[], ports=[1], order=2, freq_ghz=0.064,
        r_ohm=50.0)["Boundaries"]


def test_s1f_write_case():
    """Wwrite (-n 1): the Study writes the copper hole case; it matches the probe's."""
    from mpi4py import MPI

    from fem_em_solver.backends.palace import PalaceRunRequired

    if MPI.COMM_WORLD.size != 1:
        pytest.skip("the Palace case writer is serial")
    for p in (PROBE_CFG, PROBE_MESH, PROBE_CU):
        if not p.exists():
            pytest.skip(f"BLOCKED: {p} absent (comparand absent)")
    if (CU_CASE / "postpro" / "port-S.csv").exists():
        pytest.skip("case already run: Wread's")
    with pytest.raises(PalaceRunRequired) as ei:
        _cu_study().solve(MPI.COMM_WORLD, backend="palace", case_dir=CU_CASE)
    print(f"[s1f Wwrite] PalaceRunRequired: {ei.value}", flush=True)
    ours = json.loads((CU_CASE / "case.json").read_text())
    probe = json.loads(PROBE_CFG.read_text())
    for key in ("Domains", "Boundaries", "Solver"):
        print(f"[s1f Wwrite] case.json {key} == probe: {ours[key] == probe[key]}", flush=True)
    n_o, c_o = _census(CU_CASE / "mesh.msh")
    n_p, c_p = _census(PROBE_MESH)
    print(f"[s1f Wwrite] nodes ours {n_o} probe {n_p}", flush=True)
    print(f"[s1f Wwrite] census ours  {sorted(c_o.items())}", flush=True)
    print(f"[s1f Wwrite] census probe {sorted(c_p.items())}", flush=True)
    assert ours["Boundaries"] == probe["Boundaries"]
    assert ours["Solver"] == probe["Solver"]
    assert ours["Domains"] == probe["Domains"]
    assert n_o == n_p
    used = {(2, 1), (2, 401), (2, 21), (2, 22), (2, 23), (2, 24)}
    for key in used:
        assert c_o[key] == c_p[key], (key, c_o[key], c_p[key])
    tets_o = {k: v for k, v in c_o.items() if k[0] == 4}
    tets_p = {k: v for k, v in c_p.items() if k[0] == 4}
    assert tets_o == tets_p


@pytest.mark.skipif(os.environ.get(READ_ENV) != "1", reason=f"item 631: {READ_ENV}=1")
def test_s1f_read_back():
    """Wread: (E2-Cu), (Tie-Cu) and the PEC-hole control."""
    from mpi4py import MPI

    from fem_em_solver.backends.palace import PalaceStudyResult, read_port_s
    from tests.package.test_arch1_palace_backend import TIE2_BAND
    from tests.package.test_wf16_study_copper import NEGATIVE_CONTROL_FLOOR
    from tests.validation.test_ans6_exact_strip import ANS6_EXACT_STRIP_CLASS_MEANS_P2
    from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND

    for p in (CU_CASE / "postpro" / "port-S.csv", PROBE_CU, PROBE_PEC):
        if not p.exists():
            pytest.skip(f"BLOCKED: {p} absent (comparand absent)")
    res = _cu_study().solve(MPI.COMM_WORLD, backend="palace", case_dir=CU_CASE)
    assert isinstance(res, PalaceStudyResult) and res.order == 2
    s = np.asarray(res.s_matrix)
    s_cu, _ = read_port_s(PROBE_CU)
    s_pec, _ = read_port_s(PROBE_PEC)
    ent = np.abs(s - s_cu) / np.abs(s_cu)
    print(f"[s1f Wread E2-Cu] entry-wise |S - S_probe|/|S_probe| max {ent.max():.3e} "
          f"(band {LEG_D0_REPRODUCTION_BAND:.0e}); Frobenius-rel "
          f"{np.linalg.norm(s - s_cu) / np.linalg.norm(s_cu):.3e}; cells {res.cells}; "
          f"palace_total_s {res.palace_total_s}", flush=True)
    m = _class_means(s)
    ref = ANS6_EXACT_STRIP_CLASS_MEANS_P2["64"]["cu"]
    tie = {k: float(abs(m[k] - ref[k]) / abs(ref[k])) for k in ref}
    print("[s1f Wread Tie-Cu] " + ", ".join(f"{k} {v:.3e}" for k, v in tie.items())
          + f" (band TIE2_BAND {TIE2_BAND:.1e} from tests/package/test_arch1_palace_backend.py;"
          f" step 2b (T4) printed {T4_RECORDED}, predicted <= 1e-6)", flush=True)
    mp = _class_means(s_pec)
    miss = {k: float(abs(mp[k] - m[k]) / abs(m[k])) for k in m}
    print("[s1f Wread control] PEC hole (postpro_h05x_pec_o2) vs Study Cu, degree 2: "
          + ", ".join(f"{k} {v:.3e}" for k, v in miss.items())
          + f" (> NEGATIVE_CONTROL_FLOOR {NEGATIVE_CONTROL_FLOOR:.0e}; WF-1 degree-1 miss "
          f"{WF1_D1_MISS})", flush=True)
    assert ent.max() <= LEG_D0_REPRODUCTION_BAND, ent.max()
    for k, v in tie.items():
        assert v <= TIE2_BAND, (k, v)
    for k, v in miss.items():
        assert v > NEGATIVE_CONTROL_FLOOR, (k, v)
