"""`ANS-7` step 3 (i), degree-1 leg -- the 7-node sparse Route B (§9 item 127).

The S-derived network 4×4 Z (``z_from_s``) is solved at seven Chebyshev-extrema
nodes of the 62–67 MHz window (rounded to the 0.1 MHz grid), interpolated
entrywise in f by a barycentric degree-6 polynomial, terminated with
``route_b_z`` at C = ``C_F`` on the same 51 grid frequencies the dense records
use, and reduced by the same ``parallel_rows``.  The rows must reproduce the
tracked dense Route B rows (``metrics.json`` "rows"; our numbers) inside bands
frozen from W1 by `TH-18` step 2's rule: ``B = ceil_1sf(max(3·Δ_LOO, 1e-6))``,
Δ_LOO the largest relative row change, over the five interior nodes and both
columns, when that node is dropped (end nodes are never dropped).

Pure functions only are imported from the ``ans:7`` script; its stage /
assembly / output functions (which rewrite tracked files) are never touched.
Per-node records go to ``logs/ans7_sparse/`` -- never ``logs/ans7/B*.json``.

Knobs: ``FEM_EM_ANS7_SPARSE_MHZ`` (comma list), ``FEM_EM_ANS7_SPARSE_DEGREE``,
``FEM_EM_ANS7_SPARSE_RESOLUTION``; and (§9 item 132, degree-2 leg)
``FEM_EM_ANS7_SPARSE_FIELD_SOLVE`` (``0`` skips TH-14's P1 field solve, default
``1``) and ``FEM_EM_ANS7_SPARSE_COLUMNS`` (``cu`` skips the PEC column, default
``cu,pec``). The degree-1 tests below run only on the full 7-node, two-column,
degree-1 configuration; any other configuration solves its nodes, asserts the
per-node anchors and writes per-node JSON, and the run whose node list contains
the last node (67.0 MHz) assembles all seven JSONs at its degree
(``test_assemble_from_node_json``).

§9 item 157 (step 3 (iii), chain T11b): ``FEM_EM_ANS7_SPARSE_KAPPA`` -- unset: the
held degree-1 κ from ``metrics.json`` (behaviour, JSON names and print text
unchanged); ``in_run``: before the grid, an uncorrected Cu P1 pass at 64.0 MHz
through ``th14._build_ladder(measure_kappa=True)`` (PEC column and field solve at
their defaults) supplies the pooled κ at this degree, written by rank 0 to
``logs/ans7_sparse/kappa_p<deg>_h<tag>.json``; ``stored``: read that JSON (degree
and resolution must match; missing => fail, never the held κ). Any other value
raises ``ValueError``. With the knob set, per-node JSON carry a ``_kin`` suffix.
``test_kappa_in_run_discriminator`` asserts (K1) / (R2), (L) and (N).

§9 item 173 (`ANS-7` step 4 prerequisite, `xl` entry 24's assembly path):
``FEM_EM_ANS7_SPARSE_ASSEMBLY`` -- unset: the behaviour above; ``solved``:
``test_assemble_solved`` assembles this run's own nodes and columns in memory
(62.0 and 67.0 MHz required), Cu-only per-node JSON are written with a ``_cu``
suffix, a Cu-only run may take the kappa knob provided the field solve is off,
and the kappa pass's ladder is released before the grid factorises;
``test_assemble_from_node_json`` skips. Any other value raises ``ValueError``.
``test_assembly_subset_anchor`` (pure) checks the five-node subset
``ASSEMBLY5_NODES_MHZ`` against the seven-node rows on the stored JSONs.
"""

from __future__ import annotations

import gc
import importlib.util
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only

REPO = Path(__file__).resolve().parents[2]
CASE_DIR = REPO / "examples" / "ansys_benchmarks" / "ans7_tuned_copper_birdcage_driven_sweep"
SCRIPT = CASE_DIR / "07_tuned_copper_birdcage_driven_sweep.py"
STAGE_DIR = REPO / "logs" / "ans7"
OUT_DIR = REPO / "logs" / "ans7_sparse"

DEFAULT_NODES_MHZ = (62.0, 62.3, 63.2, 64.5, 65.8, 66.7, 67.0)
CONTROL_NODES_MHZ = (62.0, 64.5, 67.0)  # negative control (predicted)
SUBSET5_NODES_MHZ = (62.0, 62.7, 64.5, 66.3, 67.0)  # printed only (xl entry 24's count)
# §9 item 173: five of the seven defaults (both window ends and 64.5 MHz kept, so
# CONTROL_NODES_MHZ is a subset); xl entry 24's Cu-only node set.
ASSEMBLY5_NODES_MHZ = (62.0, 63.2, 64.5, 65.8, 67.0)
GRID_MHZ = tuple(round(62.0 + 0.1 * k, 6) for k in range(51))
ROWS = ("f0_hz", "r_in_ohm", "q_loaded")

# Tracked dense Route B rows (metrics.json:186–206, our numbers, public).
DENSE_ROWS = {
    "cu": {"f0_hz": 63983031.49581605, "r_in_ohm": 16998.7066800671, "q_loaded": 918.440964550853},
    "pec": {
        "f0_hz": 63997103.30801084,
        "r_in_ohm": 28425.152089802326,
        "q_loaded": 1535.9515756955145,
    },
}
NODE_S_RTOL = 1.0e-9  # backing 5.042e-14, 20260929T093403Z_ANS-7-item69-W1-B1.log:941

# Bands frozen from W1 (rule: B = ceil_1sf(max(3·Δ_LOO, 1e-6)), TH-18 step 2's).
# Frozen from W1, 20261001T095511Z_ANS-7-step3-W1.log:1514–1516 (7 nodes, -n 4,
# degree 1): Delta_LOO f0 2.503043e-14, R_in 1.943342e-11, Q 1.975765e-11 -> every
# 3·Delta_LOO is below the 1e-6 floor, so B = 1e-6 on all three rows. Never moved.
BANDS = {"f0_hz": 1e-06, "r_in_ohm": 1e-06, "q_loaded": 1e-06}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _ceil_1sf(x):
    e = np.floor(np.log10(x))
    return float(np.ceil(x / 10.0**e - 1e-9) * 10.0**e)


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


_ANS7 = _load(SCRIPT, "ans7_for_sparse")
z_from_s = _ANS7.z_from_s
route_b_z = _ANS7.route_b_z
parallel_rows = _ANS7.parallel_rows
k1_share = _ANS7.k1_share
C_F = _ANS7.C_F
K_SHARE_FLOOR = _ANS7.K_SHARE_FLOOR
KAPPA_FREQUENCY_HZ = _ANS7.ANCHOR1_FREQUENCY_HZ  # 64.0e6 (07_...py:104)
KAPPA_RECOMPUTE_ATOL = _ANS7.KAPPA_RECOMPUTE_ATOL  # 1e-9 (07_...py:487)

# §9 item 157 constants.
KAPPA_ENV = "FEM_EM_ANS7_SPARSE_KAPPA"
HELD_KAPPA = 1.0599570850048146e-02  # metrics.json "kappa" (degree 1)
# Item 132's held-kappa degree-2 rows (20261001T131348Z_ANS-7-step3d2-W3.log:1186-1192).
HELD_KAPPA_P2_ROWS = {
    "cu": {
        "f0_hz": 62657821.491631046,
        "r_in_ohm": 355.1839277557437,
        "q_loaded": 18.136752918228133,
    },
    "pec": {
        "f0_hz": 62672359.551981606,
        "r_in_ohm": 357.8640796561647,
        "q_loaded": 18.286075024732494,
    },
}
# (R2): ANS-6 step 2's degree-2 64 MHz line, uncorrected Cu P1 solve
# (20260921T100022Z_ANS-6-step2-degree2-64MHz.log:937; 10-digit print => ≈ 3.5e-10).
R2_REFERENCE_W = {"p_phantom": 1.574450007e-06, "p_surf": 1.419459787e-06}
V_SRC = 1.0  # drive_voltage_v of every spec (test_th14_birdcage_copper.py:224)
N_FACTOR = 10.0  # (N): uncorrected gap > N_FACTOR * REDUCTION_BAND
CORRECTED_NODE_MHZ = 64.5

# §9 item 173 constants.
ASSEMBLY_ENV = "FEM_EM_ANS7_SPARSE_ASSEMBLY"
# (R2) references keyed by resolution. 0.005: entry 14's uncorrected line, same
# ladder, 64 MHz, -n 16 (20260929T070008Z_ANS-6-step2c.log:1034).
R2_REFERENCE_BY_RES = {
    None: R2_REFERENCE_W,
    0.005: {"p_phantom": 1.696818288e-06, "p_surf": 1.431230505e-06},
}
# Item 157's x1 in-run p2 Cu rows (20261002T005913Z_ANS-7-step3iii-W2.log:1325-1327),
# the h-move reference for a resolution-0.005 run.
ITEM157_P2_KIN_CU_ROWS = {
    "f0_hz": 62704215.74044876,
    "r_in_ohm": 18569.60744330647,
    "q_loaded": 984.2226364364154,
}


def _nodes_mhz():
    raw = os.environ.get("FEM_EM_ANS7_SPARSE_MHZ", "").strip()
    if not raw:
        return DEFAULT_NODES_MHZ
    return tuple(float(v) for v in raw.split(",") if v.strip())


def _degree():
    return int(os.environ.get("FEM_EM_ANS7_SPARSE_DEGREE", "1") or 1)


def _resolution():
    raw = os.environ.get("FEM_EM_ANS7_SPARSE_RESOLUTION", "").strip()
    return float(raw) if raw else None


def _field_solve():
    return os.environ.get("FEM_EM_ANS7_SPARSE_FIELD_SOLVE", "1").strip() != "0"


def _columns():
    raw = os.environ.get("FEM_EM_ANS7_SPARSE_COLUMNS", "").strip()
    cols = ("cu", "pec") if not raw else tuple(c.strip() for c in raw.split(",") if c.strip())
    assert set(cols) <= {"cu", "pec"} and "cu" in cols, cols
    return cols


def _full_degree1(data):
    """127's configuration: degree 1, both columns, the seven default nodes."""
    return (
        data["degree"] == 1
        and tuple(data["columns"]) == ("cu", "pec")
        and data["field_solve"]
        and tuple(np.round(data["nodes_mhz"], 6)) == DEFAULT_NODES_MHZ
    )


def _need_full_degree1(data):
    if not _full_degree1(data):
        pytest.skip(
            "degree-1 leg test: runs only on 127's full 7-node, two-column, "
            "field-solve degree-1 configuration"
        )


def _node_json(mhz, deg, res, kin=False, cu_only=False):
    tag = "def" if res is None else f"{res:g}"
    return (
        OUT_DIR / f"f_{mhz:.1f}_p{deg}_h{tag}{'_kin' if kin else ''}"
        f"{'_cu' if cu_only else ''}.json"
    )


def _kappa_json(deg, res):
    tag = "def" if res is None else f"{res:g}"
    return OUT_DIR / f"kappa_p{deg}_h{tag}.json"


def _kappa():
    return float(json.loads((CASE_DIR / "metrics.json").read_text())["kappa"])


def _item157_constants():
    """Imported lazily (only with the knob set) so the unset path's import graph is
    unchanged: REDUCTION_BAND (test_port_lumped_rlc_termination.py:88),
    REFERENCE_IMPEDANCE_OHM (test_port_package_sparameters.py:140),
    PASSIVITY_SIGMA_TOLERANCE (test_port_birdcage_four_port.py:129)."""
    from tests.validation.test_port_birdcage_four_port import PASSIVITY_SIGMA_TOLERANCE
    from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND
    from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

    return REDUCTION_BAND, REFERENCE_IMPEDANCE_OHM, PASSIVITY_SIGMA_TOLERANCE


def _p_inc_w():
    """|V_src|^2 / (8 Z0) = 2.5e-3 W."""
    return abs(V_SRC) ** 2 / (8.0 * _item157_constants()[1])


def _kappa_mode():
    """``None`` (unset), ``"in_run"`` or ``"stored"``; anything else raises."""
    raw = os.environ.get(KAPPA_ENV, "").strip()
    if not raw:
        return None
    if raw not in ("in_run", "stored"):
        raise ValueError(f"{KAPPA_ENV}={raw!r}: expected unset, 'in_run' or 'stored'")
    return raw


def _assembly_mode():
    """``None`` (unset) or ``"solved"``; anything else raises."""
    raw = os.environ.get(ASSEMBLY_ENV, "").strip()
    if not raw:
        return None
    if raw != "solved":
        raise ValueError(f"{ASSEMBLY_ENV}={raw!r}: expected unset or 'solved'")
    return raw


def _column_power(s4):
    """Sum_i |S_i1|^2 -- the driven (P1) column's returned power fraction."""
    s4 = np.asarray(s4, dtype=np.complex128)
    return float(np.sum(np.abs(s4[:, 0]) ** 2))


def _loss_pair(s4, power):
    """(L_S, L_F): 1 - Sum|S_i1|^2 and the field-balance loss over P_inc."""
    l_s = 1.0 - _column_power(s4)
    l_f = (power["p_phantom"] + power["p_non_phantom"] + power["p_surf"]) / _p_inc_w()
    return l_s, l_f


def _kappa_in_run(th14, deg, res, comm, release=False):
    """Uncorrected Cu P1 pass at 64.0 MHz; pooled terminal-form deficit = kappa."""
    t0 = time.perf_counter()
    klad = th14._build_ladder(
        deg,
        sigmas=(th14.COPPER_SIGMA,),
        resolution=res,
        frequencies_hz=[KAPPA_FREQUENCY_HZ],
        width_correction_kappa=None,
        measure_kappa=True,
    )
    comm.Barrier()
    krec = list(klad["freqs"].values())[0]
    kappa = comm.bcast(float(krec["kappa"]["pooled"]), root=0)
    s_cu = np.asarray(krec["sigma"][th14.COPPER_SIGMA]["s"], dtype=np.complex128)
    s_pec = np.asarray(krec["pec"]["s"], dtype=np.complex128)
    kp = {
        "kappa": kappa,
        "per_sheet": dict(krec["kappa"]["per_sheet"]),
        "s_cu": s_cu,
        "s_pec": s_pec,
        "power": dict(krec["power"]),
        "identity_residual": float(krec["identity_residual"]),
        "cells": int(klad["cells"]),
        "t_s": time.perf_counter() - t0,
    }
    if release:
        # §9 item 173: drop the kappa pass's ladder (and its factor) before the grid
        # factorises -- held + grid factors would otherwise peak above the xl limit.
        del klad, krec
        gc.collect()
    if comm.rank == 0:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        cplx = lambda m: [[[complex(v).real, complex(v).imag] for v in row] for row in m]
        _kappa_json(deg, res).write_text(
            json.dumps(
                {
                    "kappa": kappa,
                    "kappa_repr": repr(kappa),
                    "per_sheet": kp["per_sheet"],
                    "degree": deg,
                    "resolution": res,
                    "mpi_ranks": comm.size,
                    "cells": kp["cells"],
                    "f_hz": KAPPA_FREQUENCY_HZ,
                    "s_cu_uncorrected": cplx(s_cu),
                    "s_pec_uncorrected": cplx(s_pec),
                    "power": kp["power"],
                    "identity_residual": kp["identity_residual"],
                }
            )
        )
    comm.Barrier()
    _say(
        f"[ans7-sparse] kappa pass (uncorrected, 64.0 MHz, p{deg}): kappa_in_run = "
        f"{kappa!r}; per sheet {kp['per_sheet']}; {kp['t_s']:.1f} s; written "
        f"{_kappa_json(deg, res)}"
    )
    return kappa, kp


def _kappa_stored(deg, res):
    p = _kappa_json(deg, res)
    if not p.exists():
        pytest.fail(f"{KAPPA_ENV}=stored: {p} absent (run in_run first); no held-kappa fallback")
    rec = json.loads(p.read_text())
    if rec["degree"] != deg or rec["resolution"] != res:
        pytest.fail(f"{p}: degree/resolution {rec['degree']}/{rec['resolution']} != {deg}/{res}")
    return float(rec["kappa"])


# -- interpolation + reduction (pure) -----------------------------------------


def bary_interp(x_nodes, y_nodes, x):
    """Barycentric Lagrange interpolation of arrays ``y_nodes[j]`` at scalar ``x``."""
    xn = np.asarray(x_nodes, dtype=float)
    w = np.array(
        [1.0 / np.prod([xj - xk for k, xk in enumerate(xn) if k != j]) for j, xj in enumerate(xn)]
    )
    d = x - xn
    hit = np.where(np.abs(d) < 1e-12)[0]
    if hit.size:
        return np.asarray(y_nodes[int(hit[0])])
    t = w / d
    num = sum(t[j] * np.asarray(y_nodes[j]) for j in range(len(xn)))
    return num / t.sum()


def sparse_rows(nodes_mhz, z_nodes, grid_mhz=GRID_MHZ, columns=("cu", "pec"), c_f=C_F):
    """Interpolated Z on the grid -> route_b_z -> parallel_rows, per column.

    ``c_f`` (§9 item 217, additive; default ``C_F`` byte-identical): the tuning
    capacitor forwarded to ``route_b_z``."""
    xn = np.asarray(nodes_mhz, dtype=float) - 64.5
    f_hz = np.array([g * 1e6 for g in grid_mhz])
    out = {}
    for name in columns:
        zs = [bary_interp(xn, z_nodes[name], g - 64.5) for g in grid_mhz]
        red = [route_b_z(z, f, c_f) for z, f in zip(zs, f_hz)]
        z_in = np.array([r[0] for r in red])
        rows = parallel_rows(f_hz, z_in)
        if "f0_hz" in rows:
            j = int(np.argmin(np.abs(f_hz - rows["f0_hz"])))
            rows["k_share"] = k1_share(red[j][2])
        out[name] = {"rows": rows, "z": zs, "z_in": z_in}
    return out


def _rel(rows, ref, key):
    if key not in rows:
        return math.inf
    return abs(rows[key] - ref[key]) / abs(ref[key])


# -- the solve (module-scoped) --------------------------------------------------


def _stored_stage_points():
    pts = {}
    for name in ("B1.json", "B2.json"):
        p = STAGE_DIR / name
        if not p.exists():
            return None
        for pt in json.loads(p.read_text())["points"]:
            pts[round(pt["f_hz"] / 1e6, 6)] = pt
    return pts


@pytest.fixture(scope="module")
def sparse():
    from tests.validation import test_th14_birdcage_copper as th14

    comm = MPI.COMM_WORLD
    nodes = _nodes_mhz()
    mode = _kappa_mode()
    deg, res = _degree(), _resolution()
    columns, field_solve = _columns(), _field_solve()
    assembly = _assembly_mode()
    kpass = None
    if mode is None:
        kappa = _kappa()
        _say(
            f"[ans7-sparse] nodes {nodes} MHz, degree {deg}, resolution {res}, "
            f"kappa {kappa:.16e} (degree-1 value, held), columns {columns}, "
            f"field_solve {field_solve}, -n {comm.size}"
        )
    else:
        if assembly == "solved":
            if "pec" not in columns and field_solve:
                # Trap (item 157): pec_column=False with field_solve=True raises on rank 0 only.
                raise ValueError(
                    f"{KAPPA_ENV}={mode} with columns {columns} needs "
                    f"FEM_EM_ANS7_SPARSE_FIELD_SOLVE=0"
                )
        elif "pec" not in columns:
            # Trap (item 157): pec_column=False with field_solve=True raises on rank 0 only.
            raise ValueError(f"{KAPPA_ENV}={mode} needs both columns (got {columns})")
        if mode == "in_run":
            kappa, kpass = _kappa_in_run(th14, deg, res, comm, release=(assembly == "solved"))
        else:
            kappa = _kappa_stored(deg, res)
        _say(
            f"[ans7-sparse] nodes {nodes} MHz, degree {deg}, resolution {res}, "
            f"kappa {kappa!r} ({KAPPA_ENV}={mode}, degree {deg}), columns {columns}, "
            f"field_solve {field_solve}, -n {comm.size}"
        )
    t0 = time.perf_counter()
    lad = th14._build_ladder(
        deg,
        sigmas=(th14.COPPER_SIGMA,),
        resolution=res,
        frequencies_hz=[m * 1e6 for m in nodes],
        width_correction_kappa=kappa,
        pec_column=("pec" in columns),
        field_solve=field_solve,
    )
    comm.Barrier()
    t_solve = time.perf_counter() - t0
    recs = sorted(lad["freqs"].values(), key=lambda r: r["f_hz"])
    assert len(recs) == len(nodes)
    data = {
        "nodes_mhz": nodes,
        "s": {c: [] for c in columns},
        "z": {c: [] for c in columns},
        "recip": {c: [] for c in columns},
        "identity": [],
        "cells": int(lad["cells"]),
        "t_solve_s": t_solve,
        "degree": deg,
        "resolution": res,
        "kappa": kappa,
        "columns": columns,
        "field_solve": field_solve,
        "kappa_mode": mode,
        "kin": mode is not None,
        "kappa_pass": kpass,
        "assembly": assembly,
        "sigma_max": {c: [] for c in columns},
        "power": [],
    }
    if comm.rank == 0:
        OUT_DIR.mkdir(parents=True, exist_ok=True)
    for mhz, rec in zip(nodes, recs):
        assert abs(rec["f_hz"] - mhz * 1e6) < 1.0
        cols = {"cu": rec["sigma"][th14.COPPER_SIGMA]}
        if "pec" in columns:
            cols["pec"] = rec["pec"]
        ident = float(rec["identity_residual"]) if field_solve else None
        payload = {
            "f_hz": rec["f_hz"],
            "degree": deg,
            "resolution": res,
            "kappa": kappa,
            "mpi_ranks": comm.size,
            "cells": int(lad["cells"]),
            "field_solve": field_solve,
            "identity_residual": ident,
        }
        if ident is not None:
            data["identity"].append(ident)
        if field_solve:
            data["power"].append(dict(rec["power"]))
        if mode is not None:
            payload["kappa_mode"] = mode
            payload["power"] = dict(rec["power"]) if field_solve else None
        for name, sw in cols.items():
            s4 = np.asarray(sw["s"], dtype=np.complex128)
            data["s"][name].append(s4)
            data["z"][name].append(z_from_s(s4))
            data["recip"][name].append(float(sw["reciprocity"]))
            data["sigma_max"][name].append(float(np.max(sw["sigma"])))
            payload[name] = {
                "s": [[[complex(v).real, complex(v).imag] for v in row] for row in s4],
                "reciprocity": float(sw["reciprocity"]),
            }
            if mode is not None:
                payload[name]["sigma_max"] = float(np.max(sw["sigma"]))
                z_in = route_b_z(data["z"][name][-1], mhz * 1e6, C_F)[0]
                _say(
                    f"[ans7-sparse] node {mhz:.1f} MHz p{deg} {name}: Im Y_in = "
                    f"{(1.0 / z_in).imag:+.9e} S (Z_in {z_in.real:.6e}{z_in.imag:+.6e}j ohm); "
                    f"column power sum|S_i1|^2 = {_column_power(s4):.9e}"
                )
        if comm.rank == 0 and tuple(columns) == ("cu", "pec"):
            _node_json(mhz, deg, res, kin=data["kin"]).write_text(json.dumps(payload))
        elif comm.rank == 0 and assembly == "solved" and tuple(columns) == ("cu",):
            _node_json(mhz, deg, res, kin=data["kin"], cu_only=True).write_text(json.dumps(payload))
    comm.Barrier()
    _say(
        f"[ans7-sparse] {len(nodes)} nodes solved in {t_solve:.1f} s "
        f"({data['cells']} cells); per-node JSON in {OUT_DIR}"
    )
    return data


def _subset(data, sub, columns=("cu", "pec")):
    idx = [list(np.round(data["nodes_mhz"], 6)).index(round(m, 6)) for m in sub]
    return {name: [data["z"][name][i] for i in idx] for name in columns}


def _loo(data, columns=("cu", "pec")):
    nodes = list(data["nodes_mhz"])
    full = sparse_rows(nodes, data["z"], columns=columns)
    worst = {k: 0.0 for k in ROWS}
    for i in range(1, len(nodes) - 1):  # interior nodes only
        keep = [m for j, m in enumerate(nodes) if j != i]
        loo = sparse_rows(keep, _subset(data, keep, columns), columns=columns)
        for name in columns:
            for k in ROWS:
                d = _rel(loo[name]["rows"], full[name]["rows"], k)
                _say(f"[ans7-sparse] LOO drop {nodes[i]:.1f} MHz {name} {k}: rel {d:.3e}")
                worst[k] = max(worst[k], d)
    return full, worst


# -- tests ---------------------------------------------------------------------


@complex_only
def test_band_candidate(sparse):
    _need_full_degree1(sparse)
    full, worst = _loo(sparse)
    for name in ("cu", "pec"):
        r = full[name]["rows"]
        _say(
            f"[ans7-sparse] {name} sparse rows: f0 {r.get('f0_hz')!r} Hz, R_in "
            f"{r.get('r_in_ohm')!r} ohm, Q {r.get('q_loaded')!r}"
        )
    cand = {k: _ceil_1sf(max(3.0 * worst[k], 1e-6)) for k in ROWS}
    for k in ROWS:
        _say(f"[ans7-sparse] Delta_LOO {k} = {worst[k]:.6e} -> band candidate B = {cand[k]:.0e}")
    assert all(np.isfinite(worst[k]) for k in ROWS), worst


@complex_only
def test_sparse_reproduces_dense(sparse):
    _need_full_degree1(sparse)
    assert all(BANDS[k] is not None for k in ROWS), "bands not frozen (run W1 first)"
    nodes = list(sparse["nodes_mhz"])
    # Per-node anchors: reciprocity, identity, node S vs stored dense S.
    for i, mhz in enumerate(nodes):
        for name in ("cu", "pec"):
            assert sparse["recip"][name][i] <= _th14().RECIPROCITY_BAND, (mhz, name)
        assert sparse["identity"][i] <= _th14().DISCRETE_IDENTITY_RTOL, (mhz, sparse["identity"][i])
    stored = _stored_stage_points()
    dense_pts = json.loads((CASE_DIR / "metrics.json").read_text())["points"]
    if stored is None:
        _say(
            "[ans7-sparse] logs/ans7/B1.json or B2.json absent: node-S leg skipped; "
            "Z_in at nodes vs metrics.json points instead"
        )
        by_f = {round(p["f_hz"] / 1e6, 6): p for p in dense_pts}
        for i, mhz in enumerate(nodes):
            for name in ("cu", "pec"):
                z_in = route_b_z(sparse["z"][name][i], mhz * 1e6, C_F)[0]
                zd = by_f[round(mhz, 6)][name]["Z_in"]
                zd = complex(zd["re"], zd["im"])
                dev = abs(z_in - zd) / abs(zd)
                _say(f"[ans7-sparse] node {mhz:.1f} {name} |dZ_in|/|Z_in| = {dev:.3e}")
                assert dev <= NODE_S_RTOL, (mhz, name, dev)
    else:
        worst_s = 0.0
        for i, mhz in enumerate(nodes):
            pt = stored[round(mhz, 6)]
            for name in ("cu", "pec"):
                s_old = np.array(
                    [[complex(v["re"], v["im"]) for v in row] for row in pt[name]["s"]]
                )
                dev = float(np.max(np.abs(sparse["s"][name][i] - s_old)) / np.max(np.abs(s_old)))
                worst_s = max(worst_s, dev)
                _say(f"[ans7-sparse] node {mhz:.1f} MHz {name}: max|S-S_stored|/max|S| = {dev:.3e}")
        _say(
            f"[ans7-sparse] node-S reproduction worst {worst_s:.3e} (ASSERTED <= {NODE_S_RTOL:.0e})"
        )
        assert worst_s <= NODE_S_RTOL, worst_s

    full = sparse_rows(nodes, sparse["z"])
    fails = []
    for name in ("cu", "pec"):
        rows = full[name]["rows"]
        n_up = len(rows["crossings_up_hz"])
        _say(
            f"[ans7-sparse] {name}: - -> + zeros {rows['crossings_up_hz']} (ASSERTED exactly 1); "
            f"k = +-1 share {rows.get('k_share')} (ASSERTED > {K_SHARE_FLOOR})"
        )
        if n_up != 1:
            fails.append(("crossings", name, n_up))
        if not rows.get("k_share", 0.0) > K_SHARE_FLOOR:
            fails.append(("k-share", name, rows.get("k_share")))
        for k in ROWS:
            rel = _rel(rows, DENSE_ROWS[name], k)
            _say(
                f"[ans7-sparse] {name} {k}: sparse {rows.get(k)!r} dense {DENSE_ROWS[name][k]!r} "
                f"rel {rel:.3e} (ASSERTED <= B {BANDS[k]:.0e})"
            )
            if not rel <= BANDS[k]:
                fails.append((name, k, rel))
        if not all(_rel(rows, DENSE_ROWS[name], k) <= BANDS[k] for k in ROWS):
            for i, mhz in enumerate(nodes):
                _say(f"[ans7-sparse] z_k {name} {mhz:.1f}: {sparse['z'][name][i].tolist()}")

    # Printed only: 5-node subset, off-node Z error, 51-point Z_in error.
    if all(any(abs(m - n) < 1e-6 for n in nodes) for m in SUBSET5_NODES_MHZ):
        sub = sparse_rows(SUBSET5_NODES_MHZ, _subset(sparse, SUBSET5_NODES_MHZ))
        for name in ("cu", "pec"):
            _say(
                f"[ans7-sparse] PRINTED 5-node {name}: "
                + ", ".join(
                    f"{k} rel {_rel(sub[name]['rows'], DENSE_ROWS[name], k):.3e}" for k in ROWS
                )
            )
    else:
        _say(
            "[ans7-sparse] PRINTED 5-node subset: 62.7 / 66.3 MHz not among the solved nodes "
            "-- not computable from this run"
        )
    if stored is not None:
        off_worst, zin_worst = 0.0, 0.0
        for g_i, g in enumerate(GRID_MHZ):
            pt = stored[round(g, 6)]
            for name in ("cu", "pec"):
                s_d = np.array([[complex(v["re"], v["im"]) for v in row] for row in pt[name]["s"]])
                z_d = z_from_s(s_d)
                if not any(abs(g - n) < 1e-6 for n in nodes):
                    off_worst = max(
                        off_worst,
                        float(np.max(np.abs(full[name]["z"][g_i] - z_d)) / np.max(np.abs(z_d))),
                    )
                zd = complex(pt[name]["Z_in"]["re"], pt[name]["Z_in"]["im"])
                zin_worst = max(zin_worst, abs(full[name]["z_in"][g_i] - zd) / abs(zd))
        _say(
            f"[ans7-sparse] PRINTED off-node (44 pts) max entrywise |dZ|/max|Z| = {off_worst:.3e}; "
            f"51-point max |dZ_in|/|Z_in| = {zin_worst:.3e}"
        )
    assert not fails, fails


@complex_only
def test_three_points_miss(sparse):
    """Negative control (predicted, rule (e)): the 3-node quadratic misses a band."""
    _need_full_degree1(sparse)
    assert all(BANDS[k] is not None for k in ROWS), "bands not frozen (run W1 first)"
    ctl = sparse_rows(CONTROL_NODES_MHZ, _subset(sparse, CONTROL_NODES_MHZ))
    for name in ("cu", "pec"):
        rows = ctl[name]["rows"]
        rels = {k: _rel(rows, DENSE_ROWS[name], k) for k in ROWS}
        _say(
            f"[ans7-sparse] control 3-node {name}: "
            + ", ".join(f"{k} rel {rels[k]:.3e} (B {BANDS[k]:.0e})" for k in ROWS)
        )
        missed = rels["r_in_ohm"] > BANDS["r_in_ohm"] or rels["q_loaded"] > BANDS["q_loaded"]
        assert missed, (name, rels)


@complex_only
def test_node_anchors(sparse):
    """Per-node reciprocity (every column solved) and, if the P1 solve ran, identity."""
    th14 = _th14()
    for i, mhz in enumerate(sparse["nodes_mhz"]):
        for name in sparse["columns"]:
            r = sparse["recip"][name][i]
            _say(
                f"[ans7-sparse] node {mhz:.1f} MHz p{sparse['degree']} {name}: "
                f"||S-S^T||/||S|| = {r:.3e} (ASSERTED <= {th14.RECIPROCITY_BAND:.0e})"
            )
            assert r <= th14.RECIPROCITY_BAND, (mhz, name, r)
            if sparse["kappa_mode"] is not None:
                PASSIVITY_SIGMA_TOLERANCE = _item157_constants()[2]
                smax = sparse["sigma_max"][name][i]
                _say(
                    f"[ans7-sparse] node {mhz:.1f} MHz p{sparse['degree']} {name}: sigma_max(S) "
                    f"= {smax:.12f} (ASSERTED <= 1 + {PASSIVITY_SIGMA_TOLERANCE:.0e})"
                )
                assert smax <= 1.0 + PASSIVITY_SIGMA_TOLERANCE, (mhz, name, smax)
    for i, ident in enumerate(sparse["identity"]):
        if sparse["kappa_mode"] is not None:
            _say(
                f"[ans7-sparse] node {sparse['nodes_mhz'][i]:.1f} MHz identity residual "
                f"{ident:.3e} (ASSERTED <= {th14.DISCRETE_IDENTITY_RTOL:g})"
            )
        assert ident <= th14.DISCRETE_IDENTITY_RTOL, (i, ident)


@complex_only
def test_wiring_control(sparse):
    """§9 item 132 wiring control: degree 1, Cu, ``pec_column=False, field_solve=False``
    reproduces 127's stored Cu S at 64.5 MHz within NODE_S_RTOL (backing 5.042e-14,
    20260929T093403Z_ANS-7-item69-W1-B1.log:941)."""
    if sparse["degree"] != 1 or not any(abs(m - 64.5) < 1e-6 for m in sparse["nodes_mhz"]):
        pytest.skip("wiring control: degree 1 with 64.5 MHz among the nodes only")
    i = [round(m, 6) for m in sparse["nodes_mhz"]].index(64.5)
    _say(
        f"[ans7-sparse] wiring control: columns {sparse['columns']}, field_solve "
        f"{sparse['field_solve']}"
    )
    stored = _stored_stage_points()
    if stored is None:
        pt = {
            round(p["f_hz"] / 1e6, 6): p
            for p in json.loads((CASE_DIR / "metrics.json").read_text())["points"]
        }[64.5]
        z_in = route_b_z(sparse["z"]["cu"][i], 64.5e6, C_F)[0]
        zd = complex(pt["cu"]["Z_in"]["re"], pt["cu"]["Z_in"]["im"])
        dev = abs(z_in - zd) / abs(zd)
        _say(
            f"[ans7-sparse] wiring control (stage JSON absent): |dZ_in|/|Z_in| = {dev:.3e} "
            f"(ASSERTED <= {NODE_S_RTOL:.0e})"
        )
    else:
        s_old = np.array(
            [[complex(v["re"], v["im"]) for v in row] for row in stored[64.5]["cu"]["s"]]
        )
        dev = float(np.max(np.abs(sparse["s"]["cu"][i] - s_old)) / np.max(np.abs(s_old)))
        _say(
            f"[ans7-sparse] wiring control: Cu 64.5 MHz max|S-S_stored|/max|S| = {dev:.3e} "
            f"(ASSERTED <= {NODE_S_RTOL:.0e})"
        )
    assert dev <= NODE_S_RTOL, dev


@complex_only
def test_assemble_from_node_json(sparse):
    """§9 item 132: the window whose node list holds the last node (67.0 MHz) loads the
    seven per-node JSONs at its degree and asserts crossings, k-share and
    Delta_LOO <= the frozen degree-1 BANDS (the grid's licence at this degree);
    prints the rows and the order move (p - p1)/p1 against the degree-1 rows."""
    if sparse["assembly"] == "solved":
        pytest.skip(f"{ASSEMBLY_ENV}=solved: test_assemble_solved assembles in memory")
    if not any(abs(m - DEFAULT_NODES_MHZ[-1]) < 1e-6 for m in sparse["nodes_mhz"]):
        pytest.skip("assembly runs in the window that solves 67.0 MHz")
    if tuple(sparse["columns"]) != ("cu", "pec"):
        pytest.skip("assembly needs both columns")
    deg, res, kappa = sparse["degree"], sparse["resolution"], sparse["kappa"]
    kin = sparse["kin"]
    missing = [m for m in DEFAULT_NODES_MHZ if not _node_json(m, deg, res, kin=kin).exists()]
    assert not missing, f"per-node JSON absent at p{deg}: {missing}"
    z = {"cu": [], "pec": []}
    data = {"nodes_mhz": DEFAULT_NODES_MHZ, "z": z}
    th14 = _th14()
    for m in DEFAULT_NODES_MHZ:
        rec = json.loads(_node_json(m, deg, res, kin=kin).read_text())
        assert rec["degree"] == deg and abs(rec["kappa"] - kappa) <= 1e-15 * abs(kappa), (m, rec)
        for name in ("cu", "pec"):
            s4 = np.array([[complex(a, b) for a, b in row] for row in rec[name]["s"]])
            assert rec[name]["reciprocity"] <= th14.RECIPROCITY_BAND, (m, name)
            z[name].append(z_from_s(s4))
    full, worst = _loo(data)
    fails = []
    for name in ("cu", "pec"):
        rows = full[name]["rows"]
        _say(
            f"[ans7-sparse] p{deg} {name}: - -> + zeros {rows['crossings_up_hz']} (ASSERTED "
            f"exactly 1); k-share {rows.get('k_share')} (ASSERTED > {K_SHARE_FLOOR})"
        )
        if len(rows["crossings_up_hz"]) != 1:
            fails.append(("crossings", name))
        if not rows.get("k_share", 0.0) > K_SHARE_FLOOR:
            fails.append(("k-share", name, rows.get("k_share")))
        for k in ROWS:
            _say(
                f"[ans7-sparse] PRINTED p{deg} {name} {k} = {rows.get(k)!r}; order move "
                f"(p{deg} - p1)/p1 = {(rows.get(k, math.nan) - DENSE_ROWS[name][k]) / DENSE_ROWS[name][k]:+.6e}"
            )
    for k in ROWS:
        _say(f"[ans7-sparse] p{deg} Delta_LOO {k} = {worst[k]:.6e} (ASSERTED <= B {BANDS[k]:.0e})")
        if not worst[k] <= BANDS[k]:
            fails.append(("Delta_LOO", k, worst[k]))
    if kin:
        _kin_rows_report(full, deg, kappa)
    assert not fails, fails


def _kin_rows_report(full, deg, kappa, columns=("cu", "pec")):
    """Item 157, printed only: rows beside the degree-1 and held-kappa values, both
    moves, the 1/Q split (PEC column present only), and the pre-registered
    reading-rule category."""
    held = HELD_KAPPA_P2_ROWS if deg == 2 else None
    _say(f"[ans7-sparse] PRINTED item-157 rows at p{deg}, kappa {kappa!r} (in-run)")
    cats = {}
    for name in columns:
        rows = full[name]["rows"]
        for k in ROWS:
            v = rows.get(k, math.nan)
            m1 = (v - DENSE_ROWS[name][k]) / DENSE_ROWS[name][k]
            line = (
                f"[ans7-sparse] PRINTED {name} {k}: p{deg} in-run {v!r} | p1 "
                f"{DENSE_ROWS[name][k]!r} (m = {m1:+.6e})"
            )
            if held is not None:
                mh = (v - held[name][k]) / held[name][k]
                line += f" | p2 held-kappa {held[name][k]!r} (move vs held {mh:+.6e})"
            _say(line)
            if k in ("r_in_ohm", "q_loaded"):
                if abs(m1) <= 0.30:
                    cats[(name, k)] = "kappa confound"
                elif m1 <= -0.90:
                    cats[(name, k)] = "genuine order move"
                else:
                    cats[(name, k)] = "partial"
    if "pec" in columns:
        _q_split_report(full, deg, held)
    for key, c in cats.items():
        _say(f"[ans7-sparse] reading per row: {key[0]} {key[1]} -> {c}")
    uniq = set(cats.values())
    verdict = uniq.pop() if len(uniq) == 1 else "mixed"
    _say(f"[ans7-sparse] PRINTED reading-rule category (given (L) green at W0): {verdict}")
    _say(
        "[ans7-sparse] PREDICTED (printed, never asserted): Q_Cu in [700, 1250], Q_PEC in "
        "[1100, 2200]; R_in,Cu in [13000, 24000] ohm, R_in,PEC in [21000, 43000] ohm; "
        "f0 in [61.1, 64.2] MHz; k-share >= 0.999; Delta_LOO <= 1e-9"
    )
    for name in columns:
        r = full[name]["rows"]
        _say(
            f"[ans7-sparse] MEASURED {name}: Q {r.get('q_loaded')!r}, R_in {r.get('r_in_ohm')!r}, "
            f"f0 {r.get('f0_hz')!r}, k-share {r.get('k_share')!r}"
        )


def _q_split_report(full, deg, held):
    q = {n: full[n]["rows"].get("q_loaded", math.nan) for n in ("cu", "pec")}
    q1 = {n: DENSE_ROWS[n]["q_loaded"] for n in ("cu", "pec")}
    _say(
        f"[ans7-sparse] PRINTED 1/Q split p{deg}: 1/Q_Cu {1 / q['cu']:.6e} = 1/Q_PEC "
        f"{1 / q['pec']:.6e} + copper {1 / q['cu'] - 1 / q['pec']:.6e}; p1: 1/Q_Cu "
        f"{1 / q1['cu']:.6e} = 1/Q_PEC {1 / q1['pec']:.6e} + copper "
        f"{1 / q1['cu'] - 1 / q1['pec']:.6e}"
    )
    if held is not None:
        qh = {n: held[n]["q_loaded"] for n in ("cu", "pec")}
        _say(
            f"[ans7-sparse] PRINTED 1/Q split p2 held-kappa: 1/Q_Cu {1 / qh['cu']:.6e} = "
            f"1/Q_PEC {1 / qh['pec']:.6e} + copper {1 / qh['cu'] - 1 / qh['pec']:.6e}"
        )


@complex_only
def test_kappa_in_run_discriminator(sparse):
    """§9 item 157: (K1) at degree 1 / (R2) at degree 2 on the uncorrected kappa pass;
    (L) |L_S - L_F| <= REDUCTION_BAND on the corrected 64.5 MHz node; (N) the
    uncorrected pass's gap > 10 x REDUCTION_BAND (negative control, asserted)."""
    if sparse["kappa_mode"] != "in_run":
        pytest.skip(f"discriminator: {KAPPA_ENV}=in_run only")
    kp, deg = sparse["kappa_pass"], sparse["degree"]
    kappa = kp["kappa"]
    REDUCTION_BAND = _item157_constants()[0]
    P_INC_W = _p_inc_w()
    fails = []
    _say(
        f"[ans7-sparse] kappa_in_run p{deg} = {kappa!r} ({kp['cells']} cells); per sheet "
        f"{kp['per_sheet']}; uncorrected pass identity residual {kp['identity_residual']:.3e}"
    )
    if deg == 1:
        dk = abs(kappa - HELD_KAPPA)
        _say(
            f"[ans7-sparse] (K1) |kappa_in_run - {HELD_KAPPA!r}| = {dk:.3e} "
            f"(ASSERTED <= {KAPPA_RECOMPUTE_ATOL:.0e})"
        )
        if not dk <= KAPPA_RECOMPUTE_ATOL:
            fails.append(("K1", dk))
    r2_ref = R2_REFERENCE_BY_RES.get(sparse["resolution"])
    if deg == 2 and r2_ref is None:
        _say(f"[ans7-sparse] (R2) resolution {sparse['resolution']}: no reference")
    if deg == 2 and r2_ref is not None:
        _say(
            f"[ans7-sparse] PREDICTED kappa_deg2 in [4.5e-2, 6.0e-2], point 5.2e-2; "
            f"MEASURED {kappa:.9e} ({'inside' if 4.5e-2 <= kappa <= 6.0e-2 else 'OUTSIDE'})"
        )
        for key, ref in r2_ref.items():
            got = kp["power"][key]
            rel = abs(got - ref) / abs(ref)
            _say(
                f"[ans7-sparse] (R2) uncorrected p2 64 MHz {key} = {got:.12e} W vs ANS-6 "
                f"step 2 {ref:.9e} W: rel {rel:.3e} (ASSERTED <= {NODE_S_RTOL:.0e})"
            )
            if not rel <= NODE_S_RTOL:
                fails.append(("R2", key, rel))
    # (N): the uncorrected pass.
    l_s, l_f = _loss_pair(kp["s_cu"], kp["power"])
    gap_n = abs(l_s - l_f)
    _say(
        f"[ans7-sparse] (N) uncorrected 64.0 MHz p{deg}: column power Cu "
        f"{_column_power(kp['s_cu']):.9e}, PEC {_column_power(kp['s_pec']):.9e}; L_S {l_s:.9e}, "
        f"L_F {l_f:.9e} (P_inc {P_INC_W:.3e} W; power {kp['power']}); |L_S - L_F| = "
        f"{gap_n:.6e} = {gap_n / REDUCTION_BAND:.1f} x band (ASSERTED > "
        f"{N_FACTOR:g} x {REDUCTION_BAND:.0e})"
    )
    if not gap_n > N_FACTOR * REDUCTION_BAND:
        fails.append(("N", gap_n))
    # (L): the corrected node (needs the grid's field balance).
    if not sparse["field_solve"]:
        _say(
            "[ans7-sparse] (L) SKIPPED: FEM_EM_ANS7_SPARSE_FIELD_SOLVE=0, no field balance "
            "on the grid nodes"
        )
        assert not fails, fails
        return
    nodes = [round(m, 6) for m in sparse["nodes_mhz"]]
    assert CORRECTED_NODE_MHZ in nodes, f"(L) needs the {CORRECTED_NODE_MHZ} MHz node: {nodes}"
    i = nodes.index(CORRECTED_NODE_MHZ)
    pw = sparse["power"][i]
    l_s, l_f = _loss_pair(sparse["s"]["cu"][i], pw)
    gap = abs(l_s - l_f)
    cp_cu = _column_power(sparse["s"]["cu"][i])
    cp_pec = _column_power(sparse["s"]["pec"][i])
    _say(
        f"[ans7-sparse] (L) corrected {CORRECTED_NODE_MHZ} MHz p{deg}, kappa {kappa!r}: column "
        f"power Cu {cp_cu:.9e}, PEC {cp_pec:.9e}; L_S {l_s:.9e}, L_F {l_f:.9e} (power {pw}); "
        f"|L_S - L_F| = {gap:.6e} (ASSERTED <= {REDUCTION_BAND:.0e})"
    )
    if deg == 2:
        _say(
            f"[ans7-sparse] PREDICTED (printed): column power Cu in [0.995, 0.9995] -> "
            f"{cp_cu:.6f}; PEC in [0.996, 1.000] -> {cp_pec:.6f}; gap <= 1e-4 -> {gap:.3e}"
        )
    if not gap <= REDUCTION_BAND:
        fails.append(("L", gap, kappa, cp_cu, cp_pec))
    assert not fails, fails


def _cu_z_from_json(paths):
    recs = [json.loads(p.read_text()) for p in paths]
    zs = [
        z_from_s(np.array([[complex(a, b) for a, b in row] for row in r["cu"]["s"]])) for r in recs
    ]
    return recs, zs


@complex_only
def test_assemble_solved(sparse):
    """§9 item 173: assemble this run's own nodes and columns in memory. Asserts one
    - -> + zero and k-share > K_SHARE_FLOOR per column, and (where a non-``_cu`` JSON
    set at the same degree / resolution / kappa exists) solved rows = JSON-assembled
    rows within NODE_S_RTOL. Everything else is printed."""
    if sparse["assembly"] != "solved":
        pytest.skip(f"{ASSEMBLY_ENV}=solved only")
    nodes = [round(m, 6) for m in sparse["nodes_mhz"]]
    for need in (DEFAULT_NODES_MHZ[0], DEFAULT_NODES_MHZ[-1]):
        if need not in nodes:
            raise ValueError(f"{ASSEMBLY_ENV}=solved needs {need} MHz among the nodes: {nodes}")
    columns = tuple(sparse["columns"])
    deg, res, kappa = sparse["degree"], sparse["resolution"], sparse["kappa"]
    full, worst = _loo(sparse, columns)
    fails = []
    for name in columns:
        rows = full[name]["rows"]
        _say(
            f"[ans7-sparse] solved p{deg} {name}: - -> + zeros {rows['crossings_up_hz']} "
            f"(ASSERTED exactly 1); k-share {rows.get('k_share')} (ASSERTED > {K_SHARE_FLOOR})"
        )
        if len(rows["crossings_up_hz"]) != 1:
            fails.append(("crossings", name, rows["crossings_up_hz"]))
        if not rows.get("k_share", 0.0) > K_SHARE_FLOOR:
            fails.append(("k-share", name, rows.get("k_share")))
        for k in ROWS:
            _say(
                f"[ans7-sparse] PRINTED solved p{deg} {name} {k} = {rows.get(k)!r}; order move "
                f"(p{deg} - p1)/p1 = "
                f"{(rows.get(k, math.nan) - DENSE_ROWS[name][k]) / DENSE_ROWS[name][k]:+.6e}"
            )
    if sparse["kin"]:
        _kin_rows_report(full, deg, kappa, columns)
    i62 = nodes.index(DEFAULT_NODES_MHZ[0])
    for name in columns:
        z_in = route_b_z(sparse["z"][name][i62], DEFAULT_NODES_MHZ[0] * 1e6, C_F)[0]
        _say(f"[ans7-sparse] PRINTED {name} Im Y_in(62.0 MHz) = {(1.0 / z_in).imag:+.9e} S")
    for k in ROWS:
        b5 = _ceil_1sf(max(3.0 * worst[k], 1e-6))
        _say(
            f"[ans7-sparse] PRINTED {len(nodes)}-node Delta_LOO {k} = {worst[k]:.6e} "
            f"(interior drops, {len(nodes) - 1}-node interpolants) -> B{len(nodes)} = {b5:.0e}"
        )
    if all(m in nodes for m in CONTROL_NODES_MHZ):
        ctl = sparse_rows(
            CONTROL_NODES_MHZ, _subset(sparse, CONTROL_NODES_MHZ, columns), columns=columns
        )
        for name in columns:
            _say(
                f"[ans7-sparse] PRINTED control 3-node {name}: "
                + ", ".join(
                    f"{k} {ctl[name]['rows'].get(k)!r} (rel vs solved "
                    f"{_rel(ctl[name]['rows'], full[name]['rows'], k):.3e})"
                    for k in ROWS
                )
            )
    if res is not None and abs(res - 0.005) < 1e-12:
        r = full["cu"]["rows"]
        _say(
            "[ans7-sparse] PRINTED h-move vs item 157 x1 in-run p2 Cu rows: "
            + ", ".join(
                f"{k} {(r.get(k, math.nan) - ITEM157_P2_KIN_CU_ROWS[k]) / ITEM157_P2_KIN_CU_ROWS[k]:+.6e}"
                for k in ROWS
            )
        )
    # Solved rows vs JSON-assembled rows (non-_cu set, same degree/resolution/kappa).
    matched = None
    for kin in (True, False):
        paths = [_node_json(m, deg, res, kin=kin) for m in nodes]
        if not all(p.exists() for p in paths):
            continue
        recs, zs = _cu_z_from_json(paths)
        if all(
            r["degree"] == deg and r["resolution"] == res and abs(r["kappa"] - kappa) <= 1e-15
            for r in recs
        ):
            matched = (paths[0].name, zs)
            break
    if matched is None:
        _say(
            f"[ans7-sparse] no non-_cu JSON set at p{deg}, resolution {res}, kappa {kappa!r}: "
            "solved-vs-JSON check not run"
        )
    else:
        rj = sparse_rows(nodes, {"cu": matched[1]}, columns=("cu",))["cu"]["rows"]
        for k in ROWS:
            rel = _rel(full["cu"]["rows"], rj, k)
            _say(
                f"[ans7-sparse] solved vs JSON ({matched[0]} set) cu {k}: solved "
                f"{full['cu']['rows'].get(k)!r} json {rj.get(k)!r} rel {rel:.3e} "
                f"(ASSERTED <= {NODE_S_RTOL:.0e})"
            )
            if not rel <= NODE_S_RTOL:
                fails.append(("solved-vs-json", k, rel))
    assert not fails, fails


def test_assembly_subset_anchor():
    """§9 item 173 anchor (pure): on the stored seven-node JSONs (p1 held kappa, p2
    in-run kappa), the five-node Cu rows reproduce the seven-node Cu rows within
    BANDS; negative control: the 3-node Cu rows miss B on R_in or Q (asserted at p1,
    printed at p2)."""
    if _assembly_mode() != "solved":
        pytest.skip(f"{ASSEMBLY_ENV}=solved only")
    fails = []
    for label, deg, kin in (("p1", 1, False), ("p2_kin", 2, True)):
        paths = [_node_json(m, deg, None, kin=kin) for m in DEFAULT_NODES_MHZ]
        missing = [p.name for p in paths if not p.exists()]
        if missing:
            pytest.fail(f"{label}: stored seven-node JSON absent {missing} -- stop")
        recs, zs = _cu_z_from_json(paths)
        kappas = {r["kappa"] for r in recs}
        assert all(r["degree"] == deg for r in recs) and len(kappas) == 1, (label, kappas)
        data = {"nodes_mhz": DEFAULT_NODES_MHZ, "z": {"cu": zs}}
        r7 = sparse_rows(DEFAULT_NODES_MHZ, data["z"], columns=("cu",))["cu"]["rows"]
        r5 = sparse_rows(
            ASSEMBLY5_NODES_MHZ, _subset(data, ASSEMBLY5_NODES_MHZ, ("cu",)), columns=("cu",)
        )["cu"]["rows"]
        _say(
            f"[ans7-sparse] anchor {label} (kappa {kappas.pop()!r}): 7-node cu rows "
            f"{ {k: r7.get(k) for k in ROWS} }"
        )
        for k in ROWS:
            rel = _rel(r5, r7, k)
            _say(
                f"[ans7-sparse] anchor {label} cu {k}: 5-node {r5.get(k)!r} 7-node {r7.get(k)!r} "
                f"|row5 - row7|/|row7| = {rel:.3e} (ASSERTED <= B {BANDS[k]:.0e})"
            )
            if not rel <= BANDS[k]:
                fails.append(("anchor", label, k, rel))
        ctl = sparse_rows(
            CONTROL_NODES_MHZ, _subset(data, CONTROL_NODES_MHZ, ("cu",)), columns=("cu",)
        )["cu"]["rows"]
        rels = {k: _rel(ctl, r7, k) for k in ROWS}
        missed = rels["r_in_ohm"] > BANDS["r_in_ohm"] or rels["q_loaded"] > BANDS["q_loaded"]
        tag = "ASSERTED misses" if label == "p1" else "PREDICTED misses, printed"
        _say(
            f"[ans7-sparse] control 3-node {label} cu vs 7-node: "
            + ", ".join(f"{k} rel {rels[k]:.3e}" for k in ROWS)
            + f" (B 1e-06; {tag}) -> "
            f"{'misses' if missed else 'INSIDE B'}"
        )
        if label == "p1" and not missed:
            fails.append(("control inside B", label, rels))
    assert not fails, fails


def _th14():
    from tests.validation import test_th14_birdcage_copper as th14

    return th14
