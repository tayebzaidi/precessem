"""`GEO-34` step 1b — `PORT-9`'s four-port gates at 10 MHz on the STEP import.

Step 1a (`tests/mesh/test_step_import.py`) exports the F-small gapped, sheeted
four-leg birdcage as STEP through gmsh's OCC writer and re-imports it through
`io/step_import.import_step` with the JSON group map
`tests/mesh/data/geo34_fsmall_birdcage_step.json`.  This module feeds that
imported mesh to `PORT-9`'s own gated construction —
`tests/validation/test_port_birdcage_four_port.build_four_port_sweep` through
its additive ``built=`` keyword, which replaces only the `_build` call, so the
sheet rebuild, the ``f = 0.5`` narrowing, the specs and the four-drive sweep
are the gate's own — and asserts, every band imported and none restated:

* (iv) cyclic-order **pre-gate**: sheet ``i``'s measured azimuth equals
  ``90°·(i − 1) mod 360`` within ``AZIMUTH_MATCH_DEG`` (`_circulant_classes` is
  pure index arithmetic, so a non-cyclic config would silently mix classes);
* (i) reciprocity ``‖S − Sᵀ‖/‖S‖ ≤ RECIPROCITY_BAND``;
* (ii) ``σ_max(S)`` and every column power sum ``≤ 1 + PASSIVITY_SIGMA_TOLERANCE``;
* (iii′) Z self / adjacent / opposite class spreads ``≤ ADJACENT_SPREAD_BAND``
  and control (2)'s pooled / worst separation ``≥ POOLED_SEPARATION_FLOOR``.

Negative control (asserted, binary): the config with the P2 and P3 port
groups swapped fails pre-gate (iv), both offending azimuths named.  Printed
beside it (predicted, never asserted — §9 rule (e)): the same solve's Z with
indices 2 ↔ 3 permuted and its class spreads recomputed (≈ 5.9 % predicted
from the primitive record's adjacent 17.0085 Ω vs opposite 16.0605 Ω).
Control (1), the leg-(d0) column at 1e-9, is a digit record of the primitive
mesh and is **printed only** ("printed, mesh differs").  Cell count, timings
and all values are printed beside the primitive's, never gated.

Step 1b′ (``FEM_EM_GEO34_CONGRUENT=1`` in addition): the same tests on the
import through ``geo34_fsmall_birdcage_step_congruent.json`` — 1a's config
plus ``congruent_sheets``, which binds the four gap sheets to P1's under the
snapped 90° z-rotation (`GEO-32`'s cut) — plus a structural assert that the
binding held (equal facet counts, centroid Hausdorff ≤ 1e-12 m, `GEO-31`'s
metric) and an asserted negative control: 1a's config on the same STEP file
misses the self band (1b measured 3.3151 %).

Gated by ``FEM_EM_GEO34_STEP1B=1`` (standard tier, ≈ 2–3 min at ``-n 2``;
≈ 3–4 min with the step 1b′ knob).
Complex build required::

    scripts/testing/run_and_log.sh GEO-34-step1b "docker compose exec -T fem-em-solver \\
      bash -lc 'cd /workspace && source /usr/local/bin/dolfinx-complex-mode && \\
       PYTHONPATH=/workspace/src FEM_EM_REQUIRE_COMPLEX=1 FEM_EM_GEO34_STEP1B=1 \\
       timeout -k 30 300 mpiexec -n 2 python3 -m pytest \\
       tests/mesh/test_step_import_four_port.py -v -s --tb=short'"
"""

from __future__ import annotations

import os
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.io.mesh import _interface_facet_tags
from fem_em_solver.io.step_import import import_step, load_step_config
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheets import (
    PORT_LOWER,
    PORT_UPPER,
    SHEET_IFACE,
    _build,
)
from tests.mesh.test_birdcage_port_tags import LEG_COUNT
from tests.mesh.test_step_import import CONFIG_PATH
from tests.validation.test_port_birdcage_four_port import (
    LEG_D0_Z_COLUMN,
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
    _circulant_classes,
    _class_spread,
    build_four_port_sweep,
)
from tests.validation.test_port_birdcage_lumped_column import (
    ADJACENT_SPREAD_BAND,
    STEP2_CELL_COUNT,
    _sheet_bbox_centre,
)
from tests.validation.test_port_birdcage_ring_column import AZIMUTH_MATCH_DEG
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

_PROBES = Path(__file__).resolve().parents[2] / "scripts" / "probes"
if str(_PROBES) not in sys.path:
    sys.path.insert(0, str(_PROBES))
# `GEO-31`'s sheet-congruence metric (facet counts; centroid Hausdorff to
# the rotated sheet 1), imported, never restated.
import geo31_birdcage_port_sheet_congruence as geo31  # noqa: E402

pytestmark = pytest.mark.skipif(
    os.environ.get("FEM_EM_GEO34_STEP1B") != "1",
    reason="GEO-34 step 1b: set FEM_EM_GEO34_STEP1B=1 (standard tier, ~2-3 min)",
)

# `GEO-34` step 1b′: ``FEM_EM_GEO34_CONGRUENT=1`` runs every 1b test on the
# import whose config binds the four gap sheets C4-congruent
# (``congruent_sheets``, `GEO-32`'s cut), and adds the binding's structural
# assert and the unbound-config negative control.  Unset, the module is 1b
# as parked.
CONGRUENT = os.environ.get("FEM_EM_GEO34_CONGRUENT") == "1"
CONGRUENT_CONFIG_PATH = CONFIG_PATH.parent / "geo34_fsmall_birdcage_step_congruent.json"
ACTIVE_CONFIG = CONGRUENT_CONFIG_PATH if CONGRUENT else CONFIG_PATH
congruent_only = pytest.mark.skipif(
    not CONGRUENT, reason="GEO-34 step 1b': set FEM_EM_GEO34_CONGRUENT=1"
)
# The binding's structural bar: `GEO-31`'s round-off threshold on the
# centroid Hausdorff distance (the probes' ``D_ROUNDOFF``), in value.
CENTROID_HAUSDORFF_BAR_M = 1.0e-12

# The primitive mesh's record (`PORT-9` module docstring, `GEO-19` step B),
# printed beside the imported values — never gated here.
PRIMITIVE_RECORD = {
    "reciprocity": 1.044e-14,
    "sigma_max": 0.999992805,
    "column_power_max": 0.7938,
    "spreads_pct": (0.0553, 0.0353, 0.0214),
}


def _azimuth_error_deg(measured, target):
    """Wrap-aware |measured − target| in degrees (`ring_column`'s formula)."""
    return abs((measured - target + 180.0) % 360.0 - 180.0)


def _sheet_azimuths(msh, cell_tags):
    """Measured azimuth [deg] of each port sheet ``SHEET_IFACE + i``."""
    comm = msh.comm
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    ports = list(range(1, LEG_COUNT + 1))
    tags_f = _interface_facet_tags(
        msh, cell_tags,
        {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ports},
    )
    out = {}
    for i in ports:
        c = _sheet_bbox_centre(msh, tags_f, SHEET_IFACE + i, comm)
        out[i] = float(np.degrees(np.arctan2(c[1], c[0])) % 360.0)
    return out


def _cyclic_misses(azimuths):
    """Ports whose azimuth is not ``90°·(i − 1)`` within ``AZIMUTH_MATCH_DEG``."""
    return {
        i: az
        for i, az in azimuths.items()
        if _azimuth_error_deg(az, 90.0 * (i - 1)) >= AZIMUTH_MATCH_DEG
    }


def _swapped_config():
    """Step 1a's config with the P2 and P3 port groups' selectors exchanged."""
    config = load_step_config(CONFIG_PATH)
    groups = {g["name"]: g for g in config["volume_groups"]}
    for half in ("lower", "upper"):
        a, b = groups[f"port_P2_{half}"], groups[f"port_P3_{half}"]
        a["selectors"], b["selectors"] = b["selectors"], a["selectors"]
    return config


def _congruence(msh, cell_tags, comm):
    """`GEO-31`'s sheet metrics on one mesh, through `GEO-31`'s own helpers:
    global facet counts per sheet, and the centroid Hausdorff distance [m] of
    sheet i to sheet 1 rotated by 90°·(i − 1)."""
    rung = geo31._measure_rung(msh, cell_tags, comm)
    d_centroid = None
    if comm.rank == 0:
        ref = rung["sheets"][0]["centroids"]
        d_centroid = [
            geo31._hausdorff(sheet["centroids"], geo31._rotate_z(ref, turns))
            for turns, sheet in enumerate(rung["sheets"])
        ]
    return {
        "counts": [int(c) for c in rung["counts"]],
        "d_centroid": comm.bcast(d_centroid, root=0),
    }


@pytest.fixture(scope="module")
def step_file():
    comm = MPI.COMM_WORLD
    step_dir = tempfile.mkdtemp(prefix="geo34b_") if comm.rank == 0 else None
    step_dir = comm.bcast(step_dir, root=0)
    step_path = os.path.join(step_dir, "fsmall_birdcage.step")

    t0 = time.perf_counter()
    gen = _build(True, step_export_path=step_path)
    t_gen = time.perf_counter() - t0
    del gen
    return step_path, t_gen


@pytest.fixture(scope="module")
def imported_sweep(step_file):
    comm = MPI.COMM_WORLD
    step_path, t_gen = step_file

    t0 = time.perf_counter()
    msh, cell_tags, facet_tags, diag = import_step(step_path, ACTIVE_CONFIG, comm=comm)
    t_import = time.perf_counter() - t0
    azimuths = _sheet_azimuths(msh, cell_tags)
    congruence = _congruence(msh, cell_tags, comm)

    sweep = build_four_port_sweep(
        built=(msh, cell_tags, facet_tags, diag, t_import)
    )

    t0 = time.perf_counter()
    s_msh, s_cells, _s_facets, _s_diag = import_step(
        step_path, _swapped_config(), comm=comm
    )
    swapped_azimuths = _sheet_azimuths(s_msh, s_cells)
    t_swapped = time.perf_counter() - t0

    if comm.rank == 0:
        print(
            f"\n[GEO-34 1b] config {ACTIVE_CONFIG.name}; generator rung + STEP "
            f"export {t_gen:.2f} s; import "
            f"{diag['import_wall_time_s']:.2f} s + mesh {diag['mesh_wall_time_s']:.2f} s "
            f"(total {t_import:.2f} s); sweep {sweep['sweep_time']:.2f} s at "
            f"-n {comm.size}; swapped-config import {t_swapped:.2f} s"
            f"\n[GEO-34 1b] cells {sweep['cells']} vs primitive {STEP2_CELL_COUNT} "
            f"(ratio {sweep['cells'] / STEP2_CELL_COUNT:.6f}; printed, not gated)"
            f"\n[GEO-34 1b] sheet facet counts {congruence['counts']}; centroid "
            f"Hausdorff to rotated P1 "
            + ", ".join(f"{d:.3e}" for d in congruence["d_centroid"])
            + " m",
            flush=True,
        )
    return {
        "sweep": sweep,
        "azimuths": azimuths,
        "swapped_azimuths": swapped_azimuths,
        "congruence": congruence,
        "diag": diag,
    }


@pytest.fixture(scope="module")
def unbound_sweep(step_file):
    """Step 1a's config — no ``congruent_sheets`` — on the same STEP file."""
    comm = MPI.COMM_WORLD
    step_path, _t_gen = step_file
    t0 = time.perf_counter()
    msh, cell_tags, facet_tags, diag = import_step(step_path, CONFIG_PATH, comm=comm)
    t_import = time.perf_counter() - t0
    congruence = _congruence(msh, cell_tags, comm)
    sweep = build_four_port_sweep(built=(msh, cell_tags, facet_tags, diag, t_import))
    if comm.rank == 0:
        print(
            f"\n[GEO-34 1b'] unbound (1a config) import: cells {sweep['cells']}; "
            f"sheet facet counts {congruence['counts']}; centroid Hausdorff "
            + ", ".join(f"{d:.3e}" for d in congruence["d_centroid"])
            + " m",
            flush=True,
        )
    return {"sweep": sweep, "congruence": congruence}


@complex_only
def test_iv_sheets_are_in_cyclic_order(imported_sweep):
    """Pre-gate (iv): sheet i sits at 90°·(i − 1) on the imported mesh."""
    az = imported_sweep["azimuths"]
    sweep_az = {s["tag"] - SHEET_IFACE: s["azimuth_deg"] for s in imported_sweep["sweep"]["sheets"]}
    if MPI.COMM_WORLD.rank == 0:
        for i in sorted(az):
            print(
                f"[GEO-34 1b] (iv) P{i} azimuth {az[i]:.12f} deg (narrowed sheet "
                f"{sweep_az[i]:.12f}) target {90.0 * (i - 1):.1f} err "
                f"{_azimuth_error_deg(az[i], 90.0 * (i - 1)):.3e} "
                f"(band {AZIMUTH_MATCH_DEG:.0e})",
                flush=True,
            )
    misses = _cyclic_misses(az)
    assert not misses, f"imported sheets not in cyclic order: {misses}"
    assert not _cyclic_misses(sweep_az), f"narrowed sheets out of order: {sweep_az}"


@complex_only
def test_negative_control_swapped_p2_p3_fails_the_pre_gate(imported_sweep):
    """Binary control: P2 ↔ P3 swapped in the config ⇒ (iv) fails, both named."""
    az = imported_sweep["swapped_azimuths"]
    misses = _cyclic_misses(az)
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[GEO-34 1b] control (P2<->P3 swapped): azimuths "
            + ", ".join(f"P{i}={v:.6f}" for i, v in sorted(az.items()))
            + f"; pre-gate misses {sorted(misses)}",
            flush=True,
        )
    assert set(misses) == {2, 3}, (
        f"the swapped config should fail (iv) at P2 and P3 exactly; misses {misses}"
    )


@complex_only
def test_i_reciprocity(imported_sweep):
    sweep = imported_sweep["sweep"]
    ratio = sweep["reciprocity"]
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[GEO-34 1b] (i) ||S-S^T||/||S|| = {ratio:.9e} (band "
            f"{RECIPROCITY_BAND:.0e}; primitive {PRIMITIVE_RECORD['reciprocity']:.3e})",
            flush=True,
        )
    assert ratio <= RECIPROCITY_BAND


@complex_only
def test_ii_passivity(imported_sweep):
    sweep = imported_sweep["sweep"]
    sigma_max = float(np.max(sweep["sigma"]))
    col = sweep["column_power"]
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[GEO-34 1b] (ii) sigma_max(S) = {sigma_max:.9f} (primitive "
            f"{PRIMITIVE_RECORD['sigma_max']:.9f}); column sums "
            + ", ".join(f"{v:.9f}" for v in col)
            + f" (primitive max {PRIMITIVE_RECORD['column_power_max']}); "
            f"band 1 + {PASSIVITY_SIGMA_TOLERANCE:.0e}",
            flush=True,
        )
    assert sigma_max <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for v in col:
        assert v <= 1.0 + PASSIVITY_SIGMA_TOLERANCE


@complex_only
def test_iii_c4_classes_and_separation(imported_sweep):
    sweep = imported_sweep["sweep"]
    spreads = sweep["spreads"]
    worst = max(spreads.values())
    separation = sweep["pooled"] / worst if worst > 0.0 else np.inf

    # Printed controls (never asserted): the 2 <-> 3 permuted Z, and control (1).
    z = sweep["z"]
    perm = [0, 2, 1, 3]
    z_perm = z[np.ix_(perm, perm)]
    perm_spreads = {k: _class_spread(v) for k, v in _circulant_classes(z_perm).items()}
    column = z[:, 0]
    if MPI.COMM_WORLD.rank == 0:
        p = PRIMITIVE_RECORD["spreads_pct"]
        print(
            f"[GEO-34 1b] (iii') spreads self {spreads['self'] * 100:.4f}% "
            f"adjacent {spreads['adjacent'] * 100:.4f}% opposite "
            f"{spreads['opposite'] * 100:.4f}% (band {ADJACENT_SPREAD_BAND * 100:.1f}%; "
            f"primitive {p[0]} / {p[1]} / {p[2]} %); pooled "
            f"{sweep['pooled'] * 100:.4f}% => separation {separation:.4f}x "
            f"(floor {POOLED_SEPARATION_FLOOR:.0f}x)",
            flush=True,
        )
        for name in ("self", "adjacent", "opposite"):
            print(
                f"[GEO-34 1b]     {name} mean |Z| = "
                f"{np.mean(np.abs(sweep['classes'][name])):.9e} Ohm",
                flush=True,
            )
        print(
            f"[GEO-34 1b] printed control (predicted ~5.9%, never asserted): Z with "
            f"2<->3 permuted — self {perm_spreads['self'] * 100:.4f}% adjacent "
            f"{perm_spreads['adjacent'] * 100:.4f}% opposite "
            f"{perm_spreads['opposite'] * 100:.4f}%",
            flush=True,
        )
        for k, (zk, zr) in enumerate(zip(column, LEG_D0_Z_COLUMN), start=1):
            print(
                f"[GEO-34 1b] control (1) printed, mesh differs: Z_{k}1 {zk:+.9e} vs "
                f"primitive record {zr:+.9e} rel. dev {abs(zk - zr) / abs(zr):.3e}",
                flush=True,
            )
    for name, value in spreads.items():
        assert value <= ADJACENT_SPREAD_BAND, f"{name} class spreads {value * 100:.4f}%"
    assert separation >= POOLED_SEPARATION_FLOOR


@congruent_only
@complex_only
def test_congruent_binding_held(imported_sweep):
    """Step 1b′ structural assert: the ``congruent_sheets`` binding held —
    four equal sheet facet counts, and every sheet's centroid Hausdorff to
    sheet 1 rotated by 90°·(i − 1) ≤ 1e-12 m (`GEO-31`'s metric)."""
    congruence = imported_sweep["congruence"]
    counts, d = congruence["counts"], congruence["d_centroid"]
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[GEO-34 1b'] binding: surfaces "
            f"{imported_sweep['diag'].get('congruent_sheet_surfaces')}; facet counts "
            f"{counts}; centroid Hausdorff "
            + ", ".join(f"{v:.3e}" for v in d)
            + f" m (bar {CENTROID_HAUSDORFF_BAR_M:.0e})",
            flush=True,
        )
    assert imported_sweep["diag"].get("congruent_sheet_surfaces") is not None
    assert len(set(counts)) == 1, f"sheet facet counts differ: {counts}"
    assert max(d) <= CENTROID_HAUSDORFF_BAR_M, f"centroid Hausdorff {d}"


@congruent_only
@complex_only
def test_negative_control_unbound_config_misses_self_band(unbound_sweep):
    """Asserted control (rule (e): 1b measured it, self 3.3151 %,
    `20260925T094541Z_GEO-34-step1b.log:5543`): step 1a's config, key absent,
    on the same STEP file still misses (iii′)'s self band."""
    spreads = unbound_sweep["sweep"]["spreads"]
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"[GEO-34 1b'] control (unbound 1a config): spreads self "
            f"{spreads['self'] * 100:.4f}% adjacent {spreads['adjacent'] * 100:.4f}% "
            f"opposite {spreads['opposite'] * 100:.4f}% (1b record 3.3151 / "
            f"1.1102 / 2.2713 %; band {ADJACENT_SPREAD_BAND * 100:.1f}%)",
            flush=True,
        )
    assert spreads["self"] > ADJACENT_SPREAD_BAND
