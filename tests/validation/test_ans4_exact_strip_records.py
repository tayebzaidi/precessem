"""`ANS-4` step 5 (§9 item 192, chain F15) — the x1 four-port on exact SPEC
port rectangles at 10 / 64 / 128 MHz, re-recorded.

The SPEC (`examples/ansys_benchmarks/ans4_birdcage_four_port_10_64_128MHz/
SPEC.md:85-96`) defines each port as the interior half (f = 0.5) of the 0.014 x
0.008 box mid-plane: 0.007 wide centred on the leg axis, 0.008 tall, area
5.6e-5 m^2. Every earlier `ANS-4` record sits on the *ragged* midpoint-filtered
set (+4.20e-2 off). ``leg_sheet_strip_fraction=GATED_WIDTH_FRACTION`` (item
148) makes the strip edges CAD edges, so the filter selects the exact
rectangle. Operator ruling 2026-10-02: exact rectangles, κ off
(``width_correction_kappa=None``) on every gated rung; κ is printed only.

Controlled by ``FEM_EM_ANS4_EXACT_DEGREE`` ∈ {1, 2}; unset, every test skips
(the `PORT-22` / `OPS-61` precedent); any other value raises ``ValueError``.

Anchors (asserted, bands imported, never restated):

* (F) the SPEC rectangle by CAD identity — cells vs
  ``ANS4_EXACT_STRIP_CELL_RECORD`` within ``CELL_COUNT_BAND``; per port the
  narrowed area, w = A/h, radial span, axial span h and centre radius against
  the SPEC row within ``STRIP_BAND``; every F_i = dx*dz within
  ``SHEET_AREA_BAND``.
* (G) per frequency on the κ-off rung, `PORT-11`'s gates (i), (ii), (iii′).
* (R) degree 1 at -n 2: ‖S − S_rec‖/‖S_rec‖ ≤ ``STEP3A_RECORD_RTOL`` per
  frequency (skips while the record is ``None``).
* (X) the 64 MHz κ-off S reproduces item 159 / 166's ``s_kappa_off``
  (``logs/ops61/ours_f05x_p{d}.json``) within ``STEP3A_RECORD_RTOL``.
* (N) negative control, degree 1, asserted: exact vs ragged
  ``PORT11_S_64MHZ_RECORD`` S class means move > ``KNOB_REACHES_THE_SOLVE`` in
  every class (backed: 1.539 / 0.858 / 0.713 %,
  `20261002T022338Z_OPS-61-F10b-W6-ours-p1.log`).

Printed, never gated: the 10 / 128 MHz moves against the ragged records, κ
pooled and per sheet, κ-on vs κ-off moves, p2-vs-p1 moves. Rank 0 writes
``logs/ans4_exact/ans4_exact_p{d}.json`` (gitignored). Nothing here writes to
``logs/ops61/``. Gates read Z-class spreads; moves read S-class means.

`ANS-4` step 6a (§9 item 212, chain F18a) — three additive knobs; all unset
reproduces item 192's call sequence exactly:

* ``FEM_EM_ANS4_EXACT_RESOLUTION`` (m, > 0): the base rung passes
  ``resolution=<value>``.
* ``FEM_EM_ANS4_EXACT_CONDUCTOR_FACTOR`` (in (0, 1]): the base rung passes
  ``conductor_resolution = factor x CONDUCTOR_RESOLUTION`` (the ladder
  module's convention).
* ``FEM_EM_ANS4_EXACT_MHZ`` (subset of ``10 64 128``, space/comma separated):
  frequency selector; the first selected key is the base rung, the rest reuse
  its mesh. Not a mesh knob.

§9 item 225 adds a fourth, also not a mesh knob (``_knobbed()`` ignores it):

* ``FEM_EM_ANS4_EXACT_KAPPA`` ∈ {unset, ``off``}: ``off`` skips the in-run κ
  diagnostic sweep(s) (``kappa_keys = ()``) and appends ``_koff`` to the
  output name; unset keeps today's κ rungs. Any other value raises.

"Knobbed" = either mesh knob set: (F) asserts the knob reached the mesh
(cells > 1.02 x record) in place of the record band, (R) / (X) skip, the
class moves against item 192's x1 at the same degree are printed, and
``p2_vs_p1`` is suppressed. Knobbed or selected runs write
``ans4_exact_p{d}_r<res|x>_c<cf|x>_f<keys>.json`` so item 192's files are
never overwritten.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.utils.instrumentation import report_peak_rss
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import CONDUCTOR_RESOLUTION, SHEET_AREA_BAND
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_tags import RING_RADIUS
from tests.validation.test_ans4_resolution_ladder import (
    PORT11_S_64MHZ_RECORD,
    STEP3A_RECORD_RTOL,
)
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
OUT_DIR = REPO / "logs" / "ans4_exact"
OPS61_DIR = REPO / "logs" / "ops61"
DEGREE_ENV = "FEM_EM_ANS4_EXACT_DEGREE"
RESOLUTION_ENV = "FEM_EM_ANS4_EXACT_RESOLUTION"
CONDUCTOR_FACTOR_ENV = "FEM_EM_ANS4_EXACT_CONDUCTOR_FACTOR"
MHZ_ENV = "FEM_EM_ANS4_EXACT_MHZ"
KAPPA_ENV = "FEM_EM_ANS4_EXACT_KAPPA"

# (F′) under a mesh knob: the knob reached the mesh (step 4's guard).
KNOB_CELL_FACTOR = 1.02
# Predicted cell brackets (printed, never asserted; §9 item 212):
# x0.5 conductor — congruent 259 543 plus the strip's ≈ +1 100;
# h = 0.0075 global — 281 728 / 281 731 plus the strip.
PREDICTED_CELL_BRACKETS = {
    (None, 0.5): (255_000, 270_000),
    (0.0075, None): (281_000, 287_000),
}

# The SPEC table row (SPEC.md:85-96): 0.007 x 0.008 rectangle, centre on the
# leg axis at RING_RADIUS.
SPEC_WIDTH_M = 7.0e-3
SPEC_HEIGHT_M = 8.0e-3
SPEC_AREA_M2 = SPEC_WIDTH_M * SPEC_HEIGHT_M

FREQUENCIES = {"10": FREQUENCY_HZ, "64": FREQUENCY_64_HZ, "128": FREQUENCY_128_HZ}

# (R) — degree-1 κ-off 4x4 records on the exact rectangles (117 206 cells),
# registered from W1's 10-digit prints. Version tag: complex 0.11 image,
# `20261003T093408Z_ANS-4-item192-W1-p1.log:2117-2132`; (1*) licence —
# re-record on an image change, never re-band.
ANS4_EXACT_STRIP_S_P1_RECORDS = {
    "10": np.array(
        [
            [
                -3.842800723e-01 + 1.446742900e-01j,
                +4.716956786e-01 - 4.260817999e-02j,
                +4.408055448e-01 - 7.282964157e-02j,
                +4.716830532e-01 - 4.254486412e-02j,
            ],
            [
                +4.716956786e-01 - 4.260817999e-02j,
                -3.845018140e-01 + 1.447068947e-01j,
                +4.716823857e-01 - 4.256336451e-02j,
                +4.410278495e-01 - 7.284668757e-02j,
            ],
            [
                +4.408055448e-01 - 7.282964157e-02j,
                +4.716823857e-01 - 4.256336451e-02j,
                -3.841593347e-01 + 1.445957964e-01j,
                +4.715764299e-01 - 4.250994138e-02j,
            ],
            [
                +4.716830532e-01 - 4.254486412e-02j,
                +4.410278495e-01 - 7.284668757e-02j,
                +4.715764299e-01 - 4.250994138e-02j,
                -3.843800239e-01 + 1.445888499e-01j,
            ],
        ]
    ),
    "64": np.array(
        [
            [
                +4.157171200e-02 + 5.886182510e-01j,
                +3.684563461e-01 - 2.077504868e-01j,
                +2.173343512e-01 - 2.589398483e-01j,
                +3.686888008e-01 - 2.076586280e-01j,
            ],
            [
                +3.684563461e-01 - 2.077504868e-01j,
                +4.144353278e-02 + 5.886524504e-01j,
                +3.685174622e-01 - 2.076480570e-01j,
                +2.176401263e-01 - 2.589990095e-01j,
            ],
            [
                +2.173343512e-01 - 2.589398483e-01j,
                +3.685174622e-01 - 2.076480570e-01j,
                +4.157082145e-02 + 5.883320455e-01j,
                +3.686345232e-01 - 2.074715024e-01j,
            ],
            [
                +3.686888008e-01 - 2.076586280e-01j,
                +2.176401263e-01 - 2.589990095e-01j,
                +3.686345232e-01 - 2.074715024e-01j,
                +4.117256303e-02 + 5.883124674e-01j,
            ],
        ]
    ),
    "128": np.array(
        [
            [
                +4.767733777e-01 + 5.851948383e-01j,
                +2.273832946e-01 - 2.693173629e-01j,
                +5.172992006e-02 - 2.217309528e-01j,
                +2.276764567e-01 - 2.694321035e-01j,
            ],
            [
                +2.273832946e-01 - 2.693173629e-01j,
                +4.766576688e-01 + 5.852338181e-01j,
                +2.276330310e-01 - 2.692560161e-01j,
                +5.191824593e-02 - 2.219940352e-01j,
            ],
            [
                +5.172992006e-02 - 2.217309528e-01j,
                +2.276330310e-01 - 2.692560161e-01j,
                +4.763785516e-01 + 5.850151954e-01j,
                +2.278661533e-01 - 2.693169518e-01j,
            ],
            [
                +2.276764567e-01 - 2.694321035e-01j,
                +5.191824593e-02 - 2.219940352e-01j,
                +2.278661533e-01 - 2.693169518e-01j,
                +4.762691255e-01 + 5.851283409e-01j,
            ],
        ]
    ),
}
ANS4_EXACT_STRIP_RECORD_RANK_WIDTH = 2

# Item 192's degree-2 κ-off 4x4 on the exact rectangles (x1, 117 206 cells).
# Version tag: complex 0.11 image, `-n 8`,
# `20261003T094036Z_ANS-4-item192-W3-p2.log:2199-2216`. Printed comparand only
# (step 6a's class moves under knobs), never asserted.
ANS4_EXACT_STRIP_S_P2_X1 = {
    "10": np.array(
        [
            [
                -3.794881923e-01 + 1.512025203e-01j,
                +4.702541630e-01 - 4.401493077e-02j,
                +4.388917071e-01 - 7.541584511e-02j,
                +4.702628332e-01 - 4.401806430e-02j,
            ],
            [
                +4.702541630e-01 - 4.401493077e-02j,
                -3.797009221e-01 + 1.512484230e-01j,
                +4.702080184e-01 - 4.399813995e-02j,
                +4.391589827e-01 - 7.548438455e-02j,
            ],
            [
                +4.388917071e-01 - 7.541584511e-02j,
                +4.702080184e-01 - 4.399813995e-02j,
                -3.793154173e-01 + 1.511704474e-01j,
                +4.701364391e-01 - 4.400081042e-02j,
            ],
            [
                +4.702628332e-01 - 4.401806430e-02j,
                +4.391589827e-01 - 7.548438455e-02j,
                +4.701364391e-01 - 4.400081042e-02j,
                -3.796381095e-01 + 1.512541778e-01j,
            ],
        ]
    ),
    "64": np.array(
        [
            [
                +7.517414498e-02 + 6.086614263e-01j,
                +3.596443853e-01 - 2.136356899e-01j,
                +2.022430908e-01 - 2.603315610e-01j,
                +3.596403574e-01 - 2.136372440e-01j,
            ],
            [
                +3.596443853e-01 - 2.136356899e-01j,
                +7.510059668e-02 + 6.087774258e-01j,
                +3.596674948e-01 - 2.136061775e-01j,
                +2.022799188e-01 - 2.604922388e-01j,
            ],
            [
                +2.022430908e-01 - 2.603315610e-01j,
                +3.596674948e-01 - 2.136061775e-01j,
                +7.523181989e-02 + 6.085913708e-01j,
                +3.595648287e-01 - 2.135901912e-01j,
            ],
            [
                +3.596403574e-01 - 2.136372440e-01j,
                +2.022799188e-01 - 2.604922388e-01j,
                +3.595648287e-01 - 2.135901912e-01j,
                +7.520562427e-02 + 6.087653218e-01j,
            ],
        ]
    ),
    "128": np.array(
        [
            [
                +5.202493801e-01 + 5.794475976e-01j,
                +2.121107652e-01 - 2.661131588e-01j,
                +4.168940818e-02 - 2.088773119e-01j,
                +2.121091692e-01 - 2.661069470e-01j,
            ],
            [
                +2.121107652e-01 - 2.661131588e-01j,
                +5.202339999e-01 + 5.795210963e-01j,
                +2.121228638e-01 - 2.661251292e-01j,
                +4.166524280e-02 - 2.089430895e-01j,
            ],
            [
                +4.168940818e-02 - 2.088773119e-01j,
                +2.121228638e-01 - 2.661251292e-01j,
                +5.203041173e-01 + 5.794115210e-01j,
                +2.120444600e-01 - 2.660502814e-01j,
            ],
            [
                +2.121091692e-01 - 2.661069470e-01j,
                +4.166524280e-02 - 2.089430895e-01j,
                +2.120444600e-01 - 2.660502814e-01j,
                +5.203143604e-01 + 5.794476010e-01j,
            ],
        ]
    ),
}

# Printed comparands only (no full ragged 128 MHz 4x4 on the default cut):
# the x1 driven column (S11, S21, S31) of the ragged rung.
# Degree 1: `20260909T093534Z_ANS-4-step2a.log` (step-2a outcome, x1 entries).
RAGGED_128MHZ_DRIVEN_COLUMN_P1 = {
    "self": +4.753451182e-01 + 5.808054918e-01j,
    "adjacent": +2.276311790e-01 - 2.671367065e-01j,
    "opposite": +5.270667676e-02 - 2.216416825e-01j,
}
# Degree 2: `20260909T143756Z_ANS-4-step2b.log`.
RAGGED_128MHZ_DRIVEN_COLUMN_P2 = {
    "self": +5.154873290e-01 + 5.589486684e-01j,
    "adjacent": +2.122248726e-01 - 2.561950360e-01j,
    "opposite": +4.558023712e-02 - 2.081539784e-01j,
}


def _degree():
    raw = os.environ.get(DEGREE_ENV, "").strip()
    if not raw:
        return None
    d = int(raw)
    if d not in (1, 2):
        raise ValueError(f"{DEGREE_ENV}={raw!r}: only 1 or 2")
    return d


def _resolution():
    raw = os.environ.get(RESOLUTION_ENV, "").strip()
    if not raw:
        return None
    v = float(raw)
    if not (np.isfinite(v) and v > 0.0):
        raise ValueError(f"{RESOLUTION_ENV}={raw!r}: must be a finite float > 0")
    return v


def _conductor_factor():
    raw = os.environ.get(CONDUCTOR_FACTOR_ENV, "").strip()
    if not raw:
        return None
    v = float(raw)
    if not (np.isfinite(v) and 0.0 < v <= 1.0):
        raise ValueError(f"{CONDUCTOR_FACTOR_ENV}={raw!r}: must lie in (0, 1]")
    return v


def _selected_keys():
    """Frequency keys in canonical order; unset ⇒ all three."""
    raw = os.environ.get(MHZ_ENV, "").strip()
    if not raw:
        return tuple(FREQUENCIES)
    tokens = raw.replace(",", " ").split()
    bad = [t for t in tokens if t not in FREQUENCIES]
    if bad or not tokens:
        raise ValueError(f"{MHZ_ENV}={raw!r}: tokens must be a subset of 10 64 128")
    return tuple(k for k in FREQUENCIES if k in tokens)


def _kappa_off():
    """``FEM_EM_ANS4_EXACT_KAPPA``: unset ⇒ False (κ diagnostic runs); ``off`` ⇒ True."""
    raw = os.environ.get(KAPPA_ENV, "").strip()
    if not raw:
        return False
    if raw != "off":
        raise ValueError(f"{KAPPA_ENV}={raw!r}: only unset or 'off'")
    return True


def _knobbed():
    return _resolution() is not None or _conductor_factor() is not None


needs_degree = pytest.mark.skipif(
    _degree() is None, reason=f"{DEGREE_ENV} unset (opt-in window)"
)


def _cjson(m):
    return [[[float(z.real), float(z.imag)] for z in row] for row in np.asarray(m)]


def _from_cjson(m):
    return np.array([[complex(re, im) for re, im in row] for row in m])


def _fmt_moves(moves):
    return "  ".join(f"{k} {v * 100:.4f}%" for k, v in moves.items())


@pytest.fixture(scope="module")
def exact_rungs():
    d = _degree()
    comm = MPI.COMM_WORLD
    # Step 6a knobs: parsed (and any bad value raised) on every rank before any mesh.
    keys = _selected_keys()
    res = _resolution()
    cf = _conductor_factor()
    kappa_off = _kappa_off()
    knobbed = res is not None or cf is not None
    selected = bool(os.environ.get(MHZ_ENV, "").strip())
    h_c = None if cf is None else cf * CONDUCTOR_RESOLUTION
    if comm.rank == 0:
        print(
            f"\n[ANS-4 step6a] knobs: resolution {res if res is not None else 'unset'}, "
            f"conductor factor {cf if cf is not None else 'unset'}, h_c "
            + (f"{h_c:.6e}" if h_c is not None else f"unset ({CONDUCTOR_RESOLUTION:.6e})")
            + f", frequencies {' '.join(keys)}"
            + f", kappa diagnostic {'off' if kappa_off else 'on'}",
            flush=True,
        )
        if kappa_off:
            print(
                f"[ANS-4 step6a] kappa diagnostic skipped ({KAPPA_ENV}=off)",
                flush=True,
            )
    build = {"leg_sheet_strip_fraction": GATED_WIDTH_FRACTION}
    if res is not None:
        build["resolution"] = res
    if cf is not None:
        build["conductor_resolution"] = h_c
    base_key = keys[0]
    off = {}
    off[base_key] = _four_port_rung(
        f"ANS-4 exact {base_key} MHz",
        np.zeros(4),
        frequency_hz=FREQUENCIES[base_key],
        degree=d,
        **build,
    )
    report_peak_rss(comm, label=f"ANS-4 step6a after base rung {base_key} MHz")
    for key in keys[1:]:
        off[key] = _four_port_rung(
            f"ANS-4 exact {key} MHz",
            np.zeros(4),
            frequency_hz=FREQUENCIES[key],
            degree=d,
            reuse=_reuse_of(off[base_key]),
        )
    kappa_keys = (
        ()
        if kappa_off
        else tuple(k for k in (("10", "128") if d == 1 else ("128",)) if k in keys)
    )
    on = {}
    for key in kappa_keys:
        on[key] = _four_port_rung(
            f"ANS-4 exact {key} MHz kappa in-run (diagnostic)",
            np.zeros(4),
            frequency_hz=FREQUENCIES[key],
            degree=d,
            reuse=_reuse_of(off[base_key]),
            width_correction_kappa="in_run",
        )

    means = {k: _class_means(r["s"]) for k, r in off.items()}
    kappa_moves = {k: _relative_move(_class_means(on[k]["s"]), means[k]) for k in on}
    readings = {}
    if "64" in means:
        readings["exact_vs_ragged_64"] = _relative_move(
            means["64"], _class_means(PORT11_S_64MHZ_RECORD)
        )
    if "10" in means:
        readings["exact_vs_ragged_10"] = _relative_move(
            means["10"], _class_means(LEG_D_S_MATRIX_10MHZ)
        )
    if "128" in means:
        readings["exact_vs_ragged_128_column"] = _relative_move(
            means["128"],
            RAGGED_128MHZ_DRIVEN_COLUMN_P1 if d == 1 else RAGGED_128MHZ_DRIVEN_COLUMN_P2,
        )
    if knobbed:
        x1 = ANS4_EXACT_STRIP_S_P1_RECORDS if d == 1 else ANS4_EXACT_STRIP_S_P2_X1
        readings["vs_item192_x1"] = {
            k: _relative_move(means[k], _class_means(x1[k])) for k in means
        }
    p1_path = OUT_DIR / "ans4_exact_p1.json"
    if d == 2 and not knobbed and p1_path.exists():
        p1 = json.loads(p1_path.read_text())
        readings["p2_vs_p1"] = {
            k: _relative_move(means[k], _class_means(_from_cjson(p1["s_kappa_off"][k])))
            for k in means
        }

    if comm.rank == 0:
        print(
            f"\n[ANS-4 step5] degree {d}, -n {comm.size}, cells {off[base_key]['cells']}"
        )
        for key, rung in off.items():
            print(f"[ANS-4 step5] {key} MHz kappa-off S (10 digits):", flush=True)
            for row in range(4):
                print(
                    f"    S_{row + 1}k = "
                    + "  ".join(f"{z:+.9e}" for z in rung["s"][row]),
                    flush=True,
                )
            print(
                f"    Z-class spreads {_fmt_moves(rung['spreads'])}; sigma_max "
                f"{float(np.max(rung['sigma'])):.12f}; reciprocity "
                f"{rung['reciprocity']:.3e}; sweep {rung['sweep_time']:.2f} s",
                flush=True,
            )
        for key, rung in on.items():
            print(
                f"[ANS-4 step5] {key} MHz in-run kappa (diagnostic, never gated) "
                f"pooled {rung['kappa']:.10e}; per sheet "
                + ", ".join(f"{k} {v:.6e}" for k, v in rung["kappa_per_sheet"].items())
                + " (predicted [5e-6, 5e-5])",
                flush=True,
            )
            print(
                f"    kappa-on vs kappa-off S class move {_fmt_moves(kappa_moves[key])}"
                + (" (predicted self in [5e-6, 1e-4])" if d == 1 else ""),
                flush=True,
            )
        if "exact_vs_ragged_64" in readings:
            print(
                "[ANS-4 step5] exact vs ragged PORT11_S_64MHZ_RECORD "
                f"{_fmt_moves(readings['exact_vs_ragged_64'])}"
                + (" (asserted > floor at degree 1)" if d == 1 else " (printed)"),
                flush=True,
            )
        if "exact_vs_ragged_10" in readings:
            print(
                "[ANS-4 step5] exact vs ragged LEG_D_S_MATRIX_10MHZ "
                f"{_fmt_moves(readings['exact_vs_ragged_10'])} "
                "(predicted 0.5-6 %, printed)",
                flush=True,
            )
        if "exact_vs_ragged_128_column" in readings:
            print(
                f"[ANS-4 step5] exact vs RAGGED_128MHZ_DRIVEN_COLUMN_P{d} "
                f"{_fmt_moves(readings['exact_vs_ragged_128_column'])} "
                "(predicted 0.5-6 %, printed)",
                flush=True,
            )
        if "vs_item192_x1" in readings:
            comparand = (
                "ANS4_EXACT_STRIP_S_P1_RECORDS" if d == 1 else "ANS4_EXACT_STRIP_S_P2_X1"
            )
            for key, moves in readings["vs_item192_x1"].items():
                predicted = (
                    "predicted self 0.5-3.5 %, couplings 0.2-1.5 %"
                    if key == "10"
                    else "unbracketed"
                )
                print(
                    f"[ANS-4 step6a] knobbed vs item 192 x1 ({comparand}) {key} MHz "
                    f"S class move {_fmt_moves(moves)} ({predicted}; printed, never "
                    "asserted)",
                    flush=True,
                )
        if knobbed:
            print("[ANS-4 step6a] p2 vs p1: suppressed under mesh knobs (mixes mesh and order)")
        elif "p2_vs_p1" in readings:
            for key, moves in readings["p2_vs_p1"].items():
                print(
                    f"[ANS-4 step5] p2 vs p1 exact {key} MHz {_fmt_moves(moves)} "
                    "(predicted 2-8 % at 10/128; printed)",
                    flush=True,
                )
        elif d == 2:
            print("[ANS-4 step5] p2 vs p1: ans4_exact_p1.json absent, not measured")
        OUT_DIR.mkdir(parents=True, exist_ok=True)
        if knobbed or selected or kappa_off:
            out_path = OUT_DIR / (
                f"ans4_exact_p{d}_r{res if res is not None else 'x'}"
                f"_c{cf if cf is not None else 'x'}_f{'-'.join(keys)}"
                f"{'_koff' if kappa_off else ''}.json"
            )
        else:
            out_path = OUT_DIR / f"ans4_exact_p{d}.json"
        out_path.write_text(
            json.dumps(
                {
                    "degree": d,
                    "ranks": comm.size,
                    "resolution": res,
                    "conductor_factor": cf,
                    "frequencies": list(keys),
                    "cells": int(off[base_key]["cells"]),
                    "s_kappa_off": {k: _cjson(r["s"]) for k, r in off.items()},
                    "spreads": {
                        k: {n: float(v) for n, v in r["spreads"].items()}
                        for k, r in off.items()
                    },
                    "pooled": {k: float(r["pooled"]) for k, r in off.items()},
                    "sigma": {k: [float(v) for v in r["sigma"]] for k, r in off.items()},
                    "column_power": {
                        k: [float(v) for v in r["column_power"]] for k, r in off.items()
                    },
                    "kappa": {k: r["kappa"] for k, r in on.items()},
                    "kappa_per_sheet": {k: r["kappa_per_sheet"] for k, r in on.items()},
                    "kappa_on_vs_off": kappa_moves,
                    "readings": readings,
                },
                indent=2,
            )
            + "\n"
        )
        print(f"    wrote {out_path.relative_to(REPO)}", flush=True)
    return {
        "off": off,
        "on": on,
        "readings": readings,
        "degree": d,
        "base": base_key,
        "knobbed": knobbed,
        "resolution": res,
        "conductor_factor": cf,
    }


@complex_only
@needs_degree
def test_f_the_spec_rectangle_by_cad_identity(exact_rungs):
    """(F) cells, and every port's rectangle against the SPEC row."""
    off = exact_rungs["off"][exact_rungs["base"]]
    ratio = off["cells"] / ANS4_EXACT_STRIP_CELL_RECORD
    knobbed = exact_rungs["knobbed"]
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[ANS-4 step5] (F) cells {off['cells']} vs record "
            f"{ANS4_EXACT_STRIP_CELL_RECORD} (ratio {ratio:.6f}, band {CELL_COUNT_BAND})"
        )
        if knobbed:
            bracket = PREDICTED_CELL_BRACKETS.get(
                (exact_rungs["resolution"], exact_rungs["conductor_factor"])
            )
            print(
                f"[ANS-4 step6a] (F′) knobbed cells {off['cells']} > "
                f"{KNOB_CELL_FACTOR} x {ANS4_EXACT_STRIP_CELL_RECORD} (asserted); "
                + (
                    f"predicted [{bracket[0]}, {bracket[1]}] "
                    f"({'inside' if bracket[0] <= off['cells'] <= bracket[1] else 'OUTSIDE'}"
                    ", printed)"
                    if bracket
                    else "unbracketed"
                ),
                flush=True,
            )
        for s in off["sheets"]:
            print(
                f"    sheet {s['tag']}: {s['facets']} facets; area {s['area']:.12e} "
                f"(rel {s['area'] / SPEC_AREA_M2 - 1:+.3e}); w=A/h {s['w']:.12e}; "
                f"radial span {s['w_bbox']:.12e}; h {s['h']:.12e}; centre "
                f"{s['centre_radial']:.12e}; F/dxdz {s['area_ratio_full']:.12f}",
                flush=True,
            )
    if knobbed:
        # (F′) negative control: a misspelt/ignored knob reproduces 117 206.
        assert off["cells"] > KNOB_CELL_FACTOR * ANS4_EXACT_STRIP_CELL_RECORD, (
            f"cells {off['cells']}: mesh knob did not reach the mesh"
        )
    else:
        assert abs(ratio - 1.0) <= CELL_COUNT_BAND, f"cells {off['cells']}"
    for s in off["sheets"]:
        assert abs(s["area"] / SPEC_AREA_M2 - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["w"] / SPEC_WIDTH_M - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["w_bbox"] / SPEC_WIDTH_M - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["h"] / SPEC_HEIGHT_M - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["centre_radial"] / RING_RADIUS - 1.0) <= STRIP_BAND, s["tag"]
        assert abs(s["area_ratio_full"] - 1.0) < SHEET_AREA_BAND, s["tag"]


@complex_only
@needs_degree
@pytest.mark.parametrize("key", list(FREQUENCIES))
def test_g_port11_gates_on_the_exact_rectangles(exact_rungs, key):
    """(G) `PORT-11`'s gates (i), (ii), (iii′) on the κ-off rung."""
    if key not in exact_rungs["off"]:
        pytest.skip(f"{key} MHz not selected ({MHZ_ENV})")
    rung = exact_rungs["off"][key]
    spreads = rung["spreads"]
    worst = max(spreads.values())
    separation = rung["pooled"] / worst if worst > 0.0 else np.inf
    sigma_max = float(np.max(rung["sigma"]))
    power_max = float(np.max(rung["column_power"]))
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[ANS-4 step5] (G) {key} MHz kappa-off: reciprocity "
            f"{rung['reciprocity']:.3e} (band {RECIPROCITY_BAND:.0e}); sigma_max "
            f"{sigma_max:.12f}, max column power {power_max:.12f} (1 + "
            f"{PASSIVITY_SIGMA_TOLERANCE:.0e}); Z spreads {_fmt_moves(spreads)} "
            f"(band {ADJACENT_SPREAD_BAND * 100:.1f}%); pooled "
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
@pytest.mark.parametrize("key", list(FREQUENCIES))
def test_r_degree1_records(exact_rungs, key):
    """(R) degree 1 at the record width: S within ``STEP3A_RECORD_RTOL``."""
    if exact_rungs["degree"] != 1:
        pytest.skip("records are degree 1 only")
    if exact_rungs["knobbed"]:
        pytest.skip("records are the x1 mesh; skipped under mesh knobs")
    if key not in exact_rungs["off"]:
        pytest.skip(f"{key} MHz not selected ({MHZ_ENV})")
    if MPI.COMM_WORLD.size != ANS4_EXACT_STRIP_RECORD_RANK_WIDTH:
        pytest.skip(f"records asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}")
    if ANS4_EXACT_STRIP_S_P1_RECORDS is None:
        pytest.skip("record not yet registered (W1 measures it)")
    s = exact_rungs["off"][key]["s"]
    rec = np.asarray(ANS4_EXACT_STRIP_S_P1_RECORDS[key], dtype=np.complex128)
    rel = float(np.linalg.norm(s - rec) / np.linalg.norm(rec))
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[ANS-4 step5] (R) {key} MHz ||S - S_rec||/||S_rec|| {rel:.3e} "
            f"(band {STEP3A_RECORD_RTOL:.0e})",
            flush=True,
        )
    assert rel <= STEP3A_RECORD_RTOL


@complex_only
@needs_degree
def test_x_reproduces_items_159_166(exact_rungs):
    """(X) 64 MHz κ-off S equals ``ours_f05x_p{d}.json``'s ``s_kappa_off``."""
    d = exact_rungs["degree"]
    if exact_rungs["knobbed"]:
        pytest.skip("item 159 / 166 is the x1 mesh; skipped under mesh knobs")
    if "64" not in exact_rungs["off"]:
        pytest.skip(f"64 MHz not selected ({MHZ_ENV})")
    path = OPS61_DIR / f"ours_f05x_p{d}.json"
    if not path.exists():
        if MPI.COMM_WORLD.rank == 0:
            print(f"\n[ANS-4 step5] (X) {path.name}: not measured")
        pytest.skip(f"{path.name} absent")
    data = json.loads(path.read_text())
    if int(data["ranks"]) != MPI.COMM_WORLD.size:
        pytest.skip(f"{path.name} was measured at -n {data['ranks']}")
    ref = _from_cjson(data["s_kappa_off"])
    s = exact_rungs["off"]["64"]["s"]
    rel = float(np.linalg.norm(s - ref) / np.linalg.norm(ref))
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[ANS-4 step5] (X) 64 MHz vs {path.name} s_kappa_off: rel {rel:.3e} "
            f"(band {STEP3A_RECORD_RTOL:.0e})",
            flush=True,
        )
        if rel > STEP3A_RECORD_RTOL:
            for row in range(4):
                print("    ours " + "  ".join(f"{z:+.9e}" for z in s[row]))
                print("    ref  " + "  ".join(f"{z:+.9e}" for z in ref[row]))
    assert rel <= STEP3A_RECORD_RTOL


@complex_only
@needs_degree
def test_n_exact_vs_ragged_64mhz_moves(exact_rungs):
    """(N) degree 1: every S class mean moves > ``KNOB_REACHES_THE_SOLVE``."""
    if exact_rungs["degree"] != 1:
        pytest.skip("(N) is asserted at degree 1 only")
    if "exact_vs_ragged_64" not in exact_rungs["readings"]:
        pytest.skip(f"64 MHz not selected ({MHZ_ENV})")
    moves = exact_rungs["readings"]["exact_vs_ragged_64"]
    if MPI.COMM_WORLD.rank == 0:
        print(
            f"\n[ANS-4 step5] (N) exact vs ragged 64 MHz {_fmt_moves(moves)} vs floor "
            f"{KNOB_REACHES_THE_SOLVE:.0e}",
            flush=True,
        )
    for name, value in moves.items():
        assert value > KNOB_REACHES_THE_SOLVE, f"{name} move {value:.3e}"
