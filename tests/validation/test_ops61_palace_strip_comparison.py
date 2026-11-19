"""`OPS-61` step 1c, first half (§9 item 159, chain F10b-i) — ours on the strip mesh.

F10a (§9 item 148) made each narrowed port sheet N_i an exact CAD rectangle
(``leg_sheet_strip_fraction`` = ``GATED_WIDTH_FRACTION``; 117 206 tets, 40
triangles, 5.6e-5 m^2 on every port) and exported that mesh to Palace as
``f05x``. This module drives **our** solver on the same strip mesh at 64 MHz,
κ off and κ in-run, so the ours-vs-Palace class gaps are no longer confounded
by port width.

Controlled by ``FEM_EM_OPS61_DEGREE`` ∈ {1, 2}; unset, every test skips (the
`PORT-22` precedent). Palace's S is read through the probe
(``scripts/probes/ops61_step1a_palace_ans4.py``), loaded by ``importlib``.

Anchors (asserted, bands imported, never restated):

* (F) the solved mesh is Palace's mesh: cell count equals
  ``ANS4_EXACT_STRIP_CELL_RECORD`` within ``CELL_COUNT_BAND``; every narrowed
  sheet area is 5.6e-5 and ``w = A/h`` is 7.0e-3, both within ``STRIP_BAND``.
* (G) `PORT-11`'s gates (i), (ii), (iii′) on the κ-off and κ-on rungs.
* (K) the κ-on vs κ-off self-class S move exceeds ``KNOB_REACHES_THE_SOLVE``.

Printed, never gated: ours vs Palace ``f05x`` per class (κ off and on); our
spreads beside Palace's; ours-strip vs ours-ragged ``PORT11_S_64MHZ_RECORD``;
the in-run κ pooled and per sheet; C2 (the κ move against its prediction).
Rank 0 writes the matrices and readings to ``logs/ops61/ours_f05x_p{d}.json``
(gitignored). Never reads ``z_matrix_ohm`` from Palace, nor any ``aed`` key.

`OPS-61` step 1d (§9 item 194): ``FEM_EM_OPS61_FREQUENCY_MHZ`` ∈ {10, 64,
128} (anything else raises). Unset, the behaviour above is byte-identical.
Set, the rung runs at that frequency, compares with the order-matched Palace
run (``PALACE_RUN_BY_MHZ_DEGREE``; no match => Palace prints skipped), writes
``ours_f05x_p{d}_f{MHz}.json``, compares with ``LEG_D_S_MATRIX_10MHZ`` at
10 MHz / ``PORT11_S_64MHZ_RECORD`` at 64 MHz / nothing at 128 MHz, and prints
(K) without asserting it at 10 and 128 MHz (unmeasured, rule (e)).
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.validation.test_ans4_resolution_ladder import PORT11_S_64MHZ_RECORD
from tests.validation.test_lossy_sphere_fullwave import (
    FREQUENCY_64_HZ,
    FREQUENCY_128_HZ,
)
from tests.validation.test_ops61_exact_strip_ports import (
    ANS4_EXACT_STRIP_CELL_RECORD,
    STRIP_BAND,
)
from tests.validation.test_port21_feed_sensitivity import (
    KNOB_REACHES_THE_SOLVE,
    _class_means,
    _relative_move,
)
from tests.validation.test_port_birdcage_four_port import (
    PASSIVITY_SIGMA_TOLERANCE,
    POOLED_SEPARATION_FLOOR,
)
from tests.validation.test_port_birdcage_larmor_gate import (
    LEG_D_S_MATRIX_10MHZ,
    _reuse_of,
)
from tests.validation.test_port_birdcage_leg_offset_sweep import _four_port_rung
from tests.validation.test_port_birdcage_lumped_column import ADJACENT_SPREAD_BAND
from tests.validation.test_port_gap_voltage_impedance import FREQUENCY_HZ
from tests.validation.test_port_lumped_narrowed_sheet import GATED_WIDTH_FRACTION
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

REPO = Path(__file__).resolve().parents[2]
PROBE_PATH = REPO / "scripts" / "probes" / "ops61_step1a_palace_ans4.py"
OUT_DIR = REPO / "logs" / "ops61"
DEGREE_ENV = "FEM_EM_OPS61_DEGREE"

# Strip geometry from F10a (item 148): N_i = 0.5 * dx * dz on the 0.014 x
# 0.008 sheet; w = A/h = 0.007.
STRIP_AREA_M2 = 5.6e-5
STRIP_WIDTH_M = 7.0e-3

# C2 (predicted, printed): the ragged degree-1 κ move,
# `20261001T171613Z_PORT-21-step5-W2.log:5654`.
C2_PREDICTED_KAPPA_MOVE = 1.552110e-02


# `OPS-61` step 1d (§9 item 194): the frequency knob. Unset => the item-159
# behaviour byte-identical (64 MHz, ``postpro_f05x``, ``ours_f05x_p{d}.json``).
FREQUENCY_ENV = "FEM_EM_OPS61_FREQUENCY_MHZ"
FREQUENCY_BY_MHZ = {10: FREQUENCY_HZ, 64: FREQUENCY_64_HZ, 128: FREQUENCY_128_HZ}
# Palace runs by (MHz, degree) — order-matched only.
PALACE_RUN_BY_MHZ_DEGREE = {
    (10, 1): "postpro_f05x_f10",
    (64, 1): "postpro_f05x",
    (64, 2): "postpro_f05x_o2",
    (128, 1): "postpro_f05x_f128",
}


def _frequency_mhz():
    raw = os.environ.get(FREQUENCY_ENV, "").strip()
    if not raw:
        return None
    mhz = int(raw)
    if mhz not in FREQUENCY_BY_MHZ:
        raise ValueError(f"{FREQUENCY_ENV}={raw!r}: only 10, 64 or 128")
    return mhz


def _degree():
    raw = os.environ.get(DEGREE_ENV, "").strip()
    if not raw:
        return None
    d = int(raw)
    if d not in (1, 2):
        raise ValueError(f"{DEGREE_ENV}={raw!r}: only 1 or 2")
    return d


needs_degree = pytest.mark.skipif(_degree() is None, reason=f"{DEGREE_ENV} unset (opt-in window)")


def _load_probe():
    name = "ops61_step1a_palace_ans4"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, PROBE_PATH)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _cjson(m):
    return [[[float(z.real), float(z.imag)] for z in row] for row in np.asarray(m)]


@pytest.fixture(scope="module")
def strip_rungs():
    d = _degree()
    mhz = _frequency_mhz()
    if mhz is None:  # item 159 default, byte-identical
        f_mhz, palace_run, out_name = 64, "postpro_f05x", f"ours_f05x_p{d}.json"
        ragged, ragged_label = PORT11_S_64MHZ_RECORD, "PORT11"
    else:
        f_mhz = mhz
        palace_run = PALACE_RUN_BY_MHZ_DEGREE.get((mhz, d))
        out_name = f"ours_f05x_p{d}_f{mhz}.json"
        ragged, ragged_label = {
            10: (LEG_D_S_MATRIX_10MHZ, "LEG_D_S_MATRIX_10MHZ"),
            64: (PORT11_S_64MHZ_RECORD, "PORT11"),
            128: (None, "n/a"),
        }[mhz]
    frequency_hz = FREQUENCY_BY_MHZ[f_mhz]
    off = _four_port_rung(
        f"strip {f_mhz} MHz",
        np.zeros(4),
        frequency_hz=frequency_hz,
        degree=d,
        leg_sheet_strip_fraction=GATED_WIDTH_FRACTION,
    )
    on = _four_port_rung(
        f"strip {f_mhz} MHz kappa in-run",
        np.zeros(4),
        frequency_hz=frequency_hz,
        degree=d,
        reuse=_reuse_of(off),
        width_correction_kappa="in_run",
    )

    probe = _load_probe()
    off_means = _class_means(off["s"])
    on_means = _class_means(on["s"])
    readings = {}
    ident_palace = None
    if palace_run is not None:
        palace_csv = OUT_DIR / palace_run / "port-S.csv"
        s_palace, _header = probe._read_port_s(palace_csv)
        ident_palace = probe._identities(s_palace)
        palace_means = _class_means(s_palace)
        readings["off_vs_palace"] = _relative_move(off_means, palace_means)
        readings["on_vs_palace"] = _relative_move(on_means, palace_means)
    if ragged is not None:
        readings["strip_vs_ragged"] = _relative_move(off_means, _class_means(ragged))
    readings["on_vs_off"] = _relative_move(on_means, off_means)

    comm = MPI.COMM_WORLD
    if comm.rank == 0:
        if mhz is None:
            print(f"\n[OPS-61 1c F10b] degree {d}, -n {comm.size}", flush=True)
        else:
            print(
                f"\n[OPS-61 1d item194] degree {d}, -n {comm.size}, {mhz} MHz; "
                f"Palace run {palace_run or 'no order-matched Palace run'}; "
                f"ragged comparand {ragged_label}",
                flush=True,
            )
        print(
            f"    in-run kappa (pooled) {on['kappa']:.10e}; per sheet "
            + ", ".join(f"{k} {v:.6e}" for k, v in on["kappa_per_sheet"].items())
            + " (predicted same order as ragged 1.0607621546e-02 at degree 1)",
            flush=True,
        )
        for label, rung in (("kappa-off", off), ("kappa-on", on)):
            print(
                f"    ours {label} Z-class spreads "
                + "  ".join(f"{k} {v * 100:.4f}%" for k, v in rung["spreads"].items())
                + f"; sigma_max {float(np.max(rung['sigma'])):.9f}",
                flush=True,
            )
        pal = "f05x" if mhz is None else str(palace_run).replace("postpro_", "")
        if ident_palace is not None:
            print(
                f"    Palace {pal} S-class spreads "
                + "  ".join(f"{k} {v * 100:.4f}%" for k, v in ident_palace["spreads"].items()),
                flush=True,
            )
        else:
            print("    no order-matched Palace run: Palace prints skipped", flush=True)
        if ragged is None:
            print("    ours-strip vs ours-ragged: n/a (no ragged record)", flush=True)
        for key, text in (
            ("off_vs_palace", f"ours-strip p{{d}} kappa-off vs Palace {pal}"),
            ("on_vs_palace", f"ours-strip p{{d}} kappa-on vs Palace {pal}"),
            (
                "strip_vs_ragged",
                f"ours-strip p{{d}} kappa-off vs ours-ragged {ragged_label}",
            ),
            ("on_vs_off", "ours-strip p{d} kappa-on vs kappa-off"),
        ):
            if key not in readings:
                continue
            print(
                f"    {text.format(d=d)}: |d mean|/|mean| "
                + "  ".join(f"{k} {v * 100:.4f}%" for k, v in readings[key].items())
                + " (printed, never gated)",
                flush=True,
            )
        print(
            f"    C2 (predicted, not asserted): kappa self-class S move "
            f"{readings['on_vs_off']['self']:.6e} vs predicted "
            f"~{C2_PREDICTED_KAPPA_MOVE:.6e} (ragged, degree 1)",
            flush=True,
        )
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        out_path = OUT_DIR / out_name
        out_path.write_text(
            json.dumps(
                {
                    "degree": d,
                    "ranks": comm.size,
                    "cells": int(off["cells"]),
                    "s_kappa_off": _cjson(off["s"]),
                    "s_kappa_on": _cjson(on["s"]),
                    "kappa": on["kappa"],
                    "kappa_per_sheet": on["kappa_per_sheet"],
                    "spreads_kappa_off": {k: float(v) for k, v in off["spreads"].items()},
                    "spreads_kappa_on": {k: float(v) for k, v in on["spreads"].items()},
                    "palace_spreads": (
                        ident_palace["spreads"] if ident_palace is not None else None
                    ),
                    "readings": readings,
                },
                indent=2,
            )
            + "\n"
        )
        print(f"    wrote {out_path.relative_to(REPO)}", flush=True)
    return {"off": off, "on": on, "readings": readings}


@complex_only
@needs_degree
def test_f_the_solved_mesh_is_palaces_mesh(strip_rungs):
    """(F) cell count = ``ANS4_EXACT_STRIP_CELL_RECORD``; every N_i exact."""
    off = strip_rungs["off"]
    ratio = off["cells"] / ANS4_EXACT_STRIP_CELL_RECORD
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[OPS-61 1c F10b] (F) cells {off['cells']} vs record "
            f"{ANS4_EXACT_STRIP_CELL_RECORD} (ratio {ratio:.6f}, band "
            f"{CELL_COUNT_BAND})",
            flush=True,
        )
        for s in off["sheets"]:
            print(
                f"    sheet {s['tag']}: {s['facets']} facets, area "
                f"{s['area']:.12e} (rel {s['area'] / STRIP_AREA_M2 - 1:+.3e}), "
                f"w = A/h {s['w']:.12e} (rel {s['w'] / STRIP_WIDTH_M - 1:+.3e})",
                flush=True,
            )
    assert abs(ratio - 1.0) <= CELL_COUNT_BAND, f"cells {off['cells']}"
    for s in off["sheets"]:
        assert abs(s["area"] / STRIP_AREA_M2 - 1.0) <= STRIP_BAND, s
        assert abs(s["w"] / STRIP_WIDTH_M - 1.0) <= STRIP_BAND, s


@complex_only
@needs_degree
@pytest.mark.parametrize("which", ["off", "on"])
def test_g_port11_gates_on_the_strip_mesh(strip_rungs, which):
    """(G) gates (i), (ii), (iii′) on the κ-off and κ-on rungs."""
    rung = strip_rungs[which]
    spreads = rung["spreads"]
    worst = max(spreads.values())
    separation = rung["pooled"] / worst if worst > 0.0 else np.inf
    sigma_max = float(np.max(rung["sigma"]))
    power_max = float(np.max(rung["column_power"]))
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[OPS-61 1c F10b] (G) kappa-{which}: reciprocity "
            f"{rung['reciprocity']:.3e} (band {RECIPROCITY_BAND:.0e}); sigma_max "
            f"{sigma_max:.12f}, max column power {power_max:.12f} (1 + "
            f"{PASSIVITY_SIGMA_TOLERANCE:.0e}); Z spreads "
            + "  ".join(f"{k} {v * 100:.4f}%" for k, v in spreads.items())
            + f" (band {ADJACENT_SPREAD_BAND * 100:.1f}%); pooled "
            f"{rung['pooled'] * 100:.4f}% => separation {separation:.2f}x (floor "
            f"{POOLED_SEPARATION_FLOOR:.0f}x)",
            flush=True,
        )
    assert rung["reciprocity"] <= RECIPROCITY_BAND
    assert sigma_max <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for value in rung["column_power"]:
        assert value <= 1.0 + PASSIVITY_SIGMA_TOLERANCE
    for name, value in spreads.items():
        assert value <= ADJACENT_SPREAD_BAND, f"{name} spread {value:.3e}"
    assert separation >= POOLED_SEPARATION_FLOOR


@complex_only
@needs_degree
def test_k_kappa_reaches_the_solve(strip_rungs):
    """(K) the κ-on vs κ-off self-class S move > ``KNOB_REACHES_THE_SOLVE``."""
    move = strip_rungs["readings"]["on_vs_off"]["self"]
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[OPS-61 1c F10b] (K) kappa self-class S move {move:.6e} vs floor "
            f"{KNOB_REACHES_THE_SOLVE:.0e}",
            flush=True,
        )
    mhz = _frequency_mhz()
    if mhz in (10, 128):
        # Item 194: (K)'s size at 10 / 128 MHz is unmeasured (rule (e)) —
        # printed above, not asserted. At 64 MHz (and unset) it stays asserted.
        if MPI.COMM_WORLD.rank == 0:
            print(f"    (K) at {mhz} MHz: printed, not asserted (item 194)", flush=True)
        return
    assert move > KNOB_REACHES_THE_SOLVE
