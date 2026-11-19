"""`PORT-21` step 1 — the public sensitivity table of the birdcage 4x4's class
entries to three port-model conventions.

**What this is.**  `ANS-4`'s order-matched comparison disagrees with the
independent code on the *self* class at 10 and 64 MHz and agrees on the
couplings (known-issues 2026-09-19; verdict qualitative, numbers private).
This module measures how much **our own** 4x4 moves under each port-model
convention an independent code could legitimately differ in, so a later
(private) reading can be taken against a public table.  **Nothing absolute is
claimed here** — no `S11`, no `Z_in`, no resonance, no tuning.  The identities
below are self-consistency statements on one fixture at degree 1.

**The control comes first.**  If the control solve does not reproduce the gate
record at rtol 1e-9, no variant is read (the step's own negative result).

**Variants, one knob each, all of them knobs
`test_port_birdcage_four_port.build_four_port_sweep` already exposes:**

  (a) the sheet width told as ``w/(1 + kappa)`` with the **in-run** kappa —
      `PORT-14` step 3's default-off opt-in (``width_correction_kappa``);
  (b) the full-width sheet ``f = 1.0`` against the gated interior half
      ``f = 0.5`` — **SKIPPED AND NAMED**: `build_four_port_sweep` narrows every
      sheet with the module-level ``GATED_WIDTH_FRACTION`` and exposes **no**
      width-fraction parameter (`:351-361` of that module, the ninth additive
      parameter is ``width_correction_kappa``, not the fraction).  The step
      forbids building the knob ("a variant whose knob the generator does not
      already expose is skipped and named, not built — this step adds no `src/`
      surface"), so this row of the table is empty by construction, not by
      failure;
  (c) the port-box conductor-side grading halved — ``conductor_resolution =
      0.5 * CONDUCTOR_RESOLUTION`` (`PORT-14` step 1b's fourth additive
      parameter).  This **re-meshes**, so the cell count is *not* the gate
      mesh's and is printed rather than pinned.

      **Step 1b (2026-09-21 review ruling): (c) is REPLACED by (c').**  As
      first built — halved grading on the *uncongruent* cut, 259 509 cells —
      its C4 class spread (self 2.3607 % at 10 MHz) was up to 4x the move it
      was to measure: mesh asymmetry (`GEO-30`..`-32`), not a feed
      sensitivity.  Those numbers stay in the table as a PRINTED RECORD with
      the three parked logs cited (``UNCONGRUENT_C_RECORD``); their gate
      assert is removed *because the variant is replaced*, not because the
      band moved (it did not — 0.5 % stays asserted on every run).  (c') is
      ``c4_congruent_sheets=True`` plus the halved grading, and its
      ``|dS|/|S|`` is read against its own **control'**
      (``c4_congruent_sheets=True`` alone), never against the gate-record
      control.  Both control' and (c') carry every imported identity gate.

**Asserted (identity anchors, every band imported, none restated):**
  * the control reproduces leg (d0)'s registered ``Z`` column at
    ``LEG_D0_REPRODUCTION_BAND`` (1e-9) at 10 MHz, and meshes
    ``STEP1_CELL_RECORD`` cells at every frequency;
  * **every** variant run still passes the imported reciprocity
    (``RECIPROCITY_BAND``), passivity (``PASSIVITY_SIGMA_TOLERANCE``) and C4
    class-spread (``ADJACENT_SPREAD_BAND``) gates;
  * variant (a) reproduces `PORT-14` step 3's registered residual — the
    terminated ``C = 100 pF`` 3x3 against the circuit reduction at
    ``REDUCTION_BAND`` — at its registered frequency
    ``STEP3_REGISTERED_FREQUENCY_HZ`` (64 MHz).

**Negative control (asserted).**  Each variant must move the **self** class by
more than ``KNOB_REACHES_THE_SOLVE`` (1e-6) relative.  A variant that moves
nothing did not reach the solve (the `OPS-41` knob-reaches-the-solve pattern).
The ceiling is O(1e-2) from `PORT-14` step 3's kappa ~ 1 %, so 1e-6 carries
>= 1e4 of headroom; this floor is a wiring check, never a physics band.

**Printed, never asserted:** ``|dS|/|S|`` per class per variant per frequency —
the sensitivity table — and whether each variant's self-class move *falls with
frequency*.

**Trap, from the item:** (b) at ``f = 1.0`` was the pre-`PORT-9` configuration;
had it been runnable, a red gate there would be the row's printed finding, not
a defect to fix.  It is not runnable, so it is named and skipped.

Gated on ``FEM_EM_PORT21_FREQUENCIES_MHZ`` (unset: every test skips, so
`main`'s default collection time and red set are unchanged — the `PORT-22`
step 1 precedent, 2026-09-20).  Per-frequency readings are persisted as JSON
under the gitignored ``logs/port21/`` so a second window can print the table
across frequencies the first window produced.
"""

from __future__ import annotations

import json
import os
import resource
import time
from pathlib import Path

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.ports.shares import terminal_form_deficit
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_port_sheet_prerequisite import CONDUCTOR_RESOLUTION
from tests.validation.test_port_birdcage_four_port import (
    ADJACENT_SPREAD_BAND,
    LEG_D0_REPRODUCTION_BAND,
    LEG_D0_Z_COLUMN,
    PASSIVITY_SIGMA_TOLERANCE,
    RECIPROCITY_BAND,
    _circulant_classes,
    build_four_port_sweep,
)
from tests.validation.test_port_lumped_rlc_termination import (
    REDUCTION_BAND,
    STEP1_CELL_RECORD,
    STEP1B_TERMINATIONS,
    STEP3_REGISTERED_FREQUENCY_HZ,
    TERMINATED_PORT_INDEX,
    _terminated_three_port,
    reduce_terminated_ports,
    series_rlc_impedance,
)
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM

FREQUENCIES_ENV = "FEM_EM_PORT21_FREQUENCIES_MHZ"

# `PORT-21` step 3 (§9 item 53, chain F6): the degree and resolution knobs.
# Unset means step 1 exactly (degree 1, the module RESOLUTION, all four runs).
# These are the names `xl` entry 19's command already spells.
DEGREE_ENV = "FEM_EM_PORT21_DEGREE"
RESOLUTION_ENV = "FEM_EM_PORT21_RESOLUTION"
VARIANTS_ENV = "FEM_EM_PORT21_VARIANTS"
ALL_VARIANTS = ("control", "a", "control_prime", "c_prime")

# Step 1's printed self-class moves of (a) against the gate-record control
# (`20260921T094335Z_PORT-21.log:11532`, `:11534`); step 3's W1 must reproduce
# them to the printed digit with the knobs unset (the plumbing moved nothing).
STEP1_A_SELF_MOVE_PRINTED = {10.0e6: "3.411387e-02", 64.0e6: "1.552110e-02"}

# `PORT-21` step 4 (§9 item 126, chain F11): F6's W3 degree-2 10 MHz records
# on the x1 congruent cut, `-n 8`, copied from the log (ten printed digits,
# <= 2e-10 relative rounding, 5x inside LEG_D0_REPRODUCTION_BAND).  Pinned
# here because W_R overwrites `logs/port21/f_10.000MHz_p2_hdefault.json`.
W3_LOG = "20260927T110931Z_PORT-21-step3-W3.log"
W3_FREQUENCY_HZ = 10.0e6
W3_CONTROL_PRIME_S_DIAGONAL = (  # `:2414-2417`
    -3.175079770e-01 + 1.375124242e-01j,
    -3.175270546e-01 + 1.375125017e-01j,
    -3.173720031e-01 + 1.374924189e-01j,
    -3.174169301e-01 + 1.375149840e-01j,
)
W3_A_S_DIAGONAL = (  # `:2932-2935`
    -3.798849004e-01 + 1.512793449e-01j,
    -3.799058269e-01 + 1.512796736e-01j,
    -3.797358108e-01 + 1.512558674e-01j,
    -3.797848018e-01 + 1.512811964e-01j,
)
W3_KAPPA = 5.209036705e-02  # `:2942`

# The `OPS-41` knob-reaches-the-solve floor.  Not a physics band: the moves are
# O(1e-2) by `PORT-14` step 3's kappa ~ 1 %, so this sits >= 1e4 below them and
# only catches a knob that was silently dropped on the way to the solve.
KNOB_REACHES_THE_SOLVE = 1.0e-6

# The frequency at which leg (d0)'s `Z` column was registered (the gate record
# the control must reproduce); at the two Larmor frequencies no `Z` record of
# this column exists, so the control there is the baseline of the table only.
GATE_RECORD_FREQUENCY_HZ = 10.0e6

# Variant (b)'s knob, named for the report and for any later slot that adds it.
SKIPPED_VARIANT_B = (
    "(b) full-width sheet f = 1.0: `build_four_port_sweep` narrows with the "
    "module-level GATED_WIDTH_FRACTION and exposes no width-fraction "
    "parameter; the step forbids building one (no `src/` surface here)"
)

_READINGS_DIR = Path("logs/port21")

# control' (the C4-congruent cut, `c4_congruent_sheets=True`): its cell count
# and class spreads on record, PRINTED BESIDE this run's, never asserted (the
# records were taken at a different rank width).  No congruent-cut record
# exists at 10 MHz; that reading is new here and is printed.
CONTROL_PRIME_CELL_RECORD = 116118
CONTROL_PRIME_SPREAD_RECORDS = {
    64.0e6: ((0.0475, 0.0529, 0.0366), "20260920T070008Z_ANS-4-step3d.log:2548"),
    128.0e6: ((0.0739, 0.0822, 0.0506), "20260917T070008Z_ANS-4-step3b.log:2548"),
}

# Variant (c) AS FIRST BUILT (halved conductor grading on the *uncongruent*
# cut, 259 509 cells), kept as a PRINTED RECORD only.  Its gate assert is
# removed because the variant is REPLACED by (c') - the 2026-09-21 review
# ruled its class spread (self 2.3607 % at 10 MHz, up to 4x the move it was to
# measure) mesh asymmetry (`GEO-30`..`-32`'s mechanism), not a feed
# sensitivity; the 0.5 % band was not moved.  Values copied from the three
# parked-attempt logs, per frequency: (self, adjacent, opposite) moves against
# the gate-record control, (self, adjacent, opposite) spreads in %, log.
UNCONGRUENT_C_RECORD = {
    10.0e6: (
        (5.982326e-03, 2.205190e-03, 7.502789e-04),
        (2.3607, 1.2387, 1.0399),
        "20260920T112548Z_PORT-21.log (also 20260920T111725Z_PORT-21.log)",
    ),
    64.0e6: (
        (1.607415e-02, 6.406405e-03, 1.110089e-02),
        (1.0508, 0.9719, 0.9038),
        "20260920T112548Z_PORT-21.log",
    ),
    128.0e6: (
        (1.668383e-02, 1.296737e-02, 1.846200e-02),
        (0.5553, 0.6751, 0.6914),
        "20260920T112242Z_PORT-21.log",
    ),
}
UNCONGRUENT_C_CELLS = 259509


def _frequencies_hz():
    raw = os.environ.get(FREQUENCIES_ENV, "").strip()
    if raw in ("", "0"):
        return ()
    # Commas or whitespace (step 3: `xl` entry 19 passes "10 64").
    return tuple(float(tok) * 1.0e6 for tok in raw.replace(",", " ").split())


def _degree():
    raw = os.environ.get(DEGREE_ENV, "").strip()
    return 1 if raw == "" else int(raw)


def _resolution():
    raw = os.environ.get(RESOLUTION_ENV, "").strip()
    return None if raw == "" else float(raw)


def _variants():
    raw = os.environ.get(VARIANTS_ENV, "").strip()
    if raw == "":
        return ALL_VARIANTS
    tokens = tuple(raw.replace(",", " ").split())
    unknown = [t for t in tokens if t not in ALL_VARIANTS]
    if unknown:
        raise ValueError(f"{VARIANTS_ENV}: unknown token(s) {unknown}; allowed {ALL_VARIANTS}")
    return tokens


def _knobs_set():
    """True when a step-3 knob moves the fixture off step 1's mesh or order."""
    return _degree() != 1 or _resolution() is not None


def _readings_path(frequency_hz, degree, resolution):
    """Per-frequency JSON, tagged with degree and resolution (step 3: a stale
    128 MHz row was once read from an older window's untagged JSON,
    `20260921T094335Z_PORT-21.log:11536`)."""
    res = "default" if resolution is None else f"{resolution:.4g}"
    return _READINGS_DIR / f"f_{frequency_hz / 1.0e6:.3f}MHz_p{degree}_h{res}.json"


def _rss_report(comm, phase):
    """Peak and live RSS after a phase, reduced over ranks (rank 0 prints).

    ``ru_maxrss`` (KiB on Linux) is the process high-water mark, so a phase's
    peak includes every earlier phase; ``VmRSS`` is the live set, which says
    whether an earlier sweep's factor is still held while the next one runs.
    """
    peak_kib = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    live_kib = 0.0
    try:
        with open("/proc/self/status") as fh:
            for line in fh:
                if line.startswith("VmRSS:"):
                    live_kib = float(line.split()[1])
                    break
    except OSError:
        pass
    out = {
        "peak_sum_gib": comm.allreduce(peak_kib, op=MPI.SUM) / 2**20,
        "peak_max_rank_gib": comm.allreduce(peak_kib, op=MPI.MAX) / 2**20,
        "live_sum_gib": comm.allreduce(live_kib, op=MPI.SUM) / 2**20,
    }
    if comm.rank == 0:
        print(
            f"[PORT-21 step3] RSS after {phase}: peak sum over ranks "
            f"{out['peak_sum_gib']:.3f} GiB (max rank "
            f"{out['peak_max_rank_gib']:.3f} GiB), live sum "
            f"{out['live_sum_gib']:.3f} GiB  (-n {comm.size})",
            flush=True,
        )
    return out


def _class_means(s_matrix):
    """Mean complex entry of each C4 class of ``S`` (self / adjacent / opposite)."""
    classes = _circulant_classes(np.asarray(s_matrix, dtype=np.complex128))
    return {name: complex(np.mean(v)) for name, v in classes.items()}


def _s_diagonal(s_matrix):
    """The 4x4 ``S`` diagonal as ``[re, im]`` pairs (`PORT-21` step 4: the
    weekly's gap — the JSON kept (a)'s move, not its entries)."""
    diag = np.diag(np.asarray(s_matrix, dtype=np.complex128))
    return [[float(z.real), float(z.imag)] for z in diag]


def _means_json(means):
    return {k: [v.real, v.imag] for k, v in means.items()}


def _relative_move(variant_means, control_means):
    return {
        name: float(abs(variant_means[name] - control_means[name]) / abs(control_means[name]))
        for name in ("self", "adjacent", "opposite")
    }


def _gate_readings(sweep):
    return {
        "reciprocity": float(sweep["reciprocity"]),
        "sigma_max": float(np.max(sweep["sigma"])),
        "column_power_max": float(np.max(sweep["column_power"])),
        "spreads": {k: float(v) for k, v in sweep["spreads"].items()},
        "cells": int(sweep["cells"]),
    }


def _run_one_frequency(frequency_hz):
    """The selected runs at one frequency.  (b) is skipped, named.

    Step 3: ``FEM_EM_PORT21_DEGREE`` / ``FEM_EM_PORT21_RESOLUTION`` are
    forwarded to every build and every re-solve (the kappa solve and the
    terminated 3x3 take the sweep's own ``"degree"``), and
    ``FEM_EM_PORT21_VARIANTS`` selects the runs.  Unset, the runs, their order
    and every keyword reaching the solve are step 1's.  With ``control`` not
    selected, (a) is built on the C4-congruent control-prime mesh and read against
    control' (the item's trap: keep (a) on the congruent cut).
    """
    comm = MPI.COMM_WORLD
    degree = _degree()
    resolution = _resolution()
    variants = _variants()
    record = {
        "frequency_hz": float(frequency_hz),
        "skipped": [SKIPPED_VARIANT_B],
        "degree": int(degree),
        "resolution": resolution,
        "variants": list(variants),
        "rss": {},
    }
    t_control = t_a = t_cp0 = t_c = float("nan")
    kappa = float("nan")
    control = control_prime = None
    control_means = control_prime_means = None
    record["rss"]["start"] = _rss_report(comm, "start")

    def _control_prime():
        comm.Barrier()
        t = time.perf_counter()
        sweep = build_four_port_sweep(
            frequency_hz=frequency_hz,
            c4_congruent_sheets=True,
            resolution=resolution,
            degree=degree,
        )
        comm.Barrier()
        t = time.perf_counter() - t
        means = _class_means(sweep["s"])
        record["control_prime"] = _gate_readings(sweep)
        record["control_prime"]["seconds"] = float(t)
        record["control_prime"]["class_means"] = {k: [v.real, v.imag] for k, v in means.items()}
        record["control_prime"]["s_diagonal"] = _s_diagonal(sweep["s"])
        record["rss"]["control_prime"] = _rss_report(comm, "control' sweep")
        return sweep, means, t

    if "control" in variants:
        comm.Barrier()
        t0 = time.perf_counter()
        control = build_four_port_sweep(
            frequency_hz=frequency_hz, resolution=resolution, degree=degree
        )
        comm.Barrier()
        t_control = time.perf_counter() - t0

        control_means = _class_means(control["s"])
        record["control"] = _gate_readings(control)
        record["control"]["seconds"] = float(t_control)
        record["control"]["z_column1"] = [
            [float(z.real), float(z.imag)] for z in control["z"][:, 0]
        ]
        record["control"]["class_means"] = {k: [v.real, v.imag] for k, v in control_means.items()}
        record["control"]["s_diagonal"] = _s_diagonal(control["s"])
        record["rss"]["control"] = _rss_report(comm, "control sweep")

    need_prime = (
        "control_prime" in variants
        or "c_prime" in variants
        or ("a" in variants and control is None)
    )
    if need_prime and control is None:
        # (a) needs its base first when the gate-record control is not run.
        control_prime, control_prime_means, t_cp0 = _control_prime()

    # --- variant (a): the width told as w/(1 + kappa), kappa computed in-run ---
    if "a" in variants:
        from tests.validation.test_birdcage_b1_plus_map import _solve_driven

        base, base_means, base_key = (
            (control, control_means, "control")
            if control is not None
            else (control_prime, control_prime_means, "control_prime")
        )
        comm.Barrier()
        t1 = time.perf_counter()
        solved = _solve_driven(base, "P1")
        sheets = [spec.sheet(driven=(spec.port_id == "P1")) for spec in base["specs"]]
        kappa = float(
            terminal_form_deficit(
                base["mesh"],
                base["facet_tags"],
                sheets,
                solved["fields"].e_complex,
                comm,
            )["pooled"]
        )
        del solved
        record["rss"]["kappa_solve"] = _rss_report(comm, "(a) kappa solve")
        reuse = {
            "mesh": base["mesh"],
            "cell_tags": base["cell_tags"],
            "facet_tags": base["facet_tags"],
            "sheets": base["sheets"],
            "halves": base["halves"],
            "cells": base["cells"],
        }
        variant_a = build_four_port_sweep(
            frequency_hz=frequency_hz,
            reuse=reuse,
            width_correction_kappa=kappa,
            degree=degree,
        )
        comm.Barrier()
        t_a = time.perf_counter() - t1

        record["a"] = _gate_readings(variant_a)
        record["a"]["seconds"] = float(t_a)
        record["a"]["kappa"] = kappa
        record["a"]["base"] = base_key
        a_means = _class_means(variant_a["s"])
        record["a"]["class_means"] = _means_json(a_means)
        record["a"]["s_diagonal"] = _s_diagonal(variant_a["s"])
        record["a"]["move"] = _relative_move(a_means, base_means)
        record["rss"]["a"] = _rss_report(comm, "(a) sweep")

        # `PORT-14` step 3's registered residual, at its registered frequency
        # only, and only on its own fixture (the gate-record control's mesh at
        # degree 1) — elsewhere it is not that registered reading.
        if (
            abs(frequency_hz - STEP3_REGISTERED_FREQUENCY_HZ) <= 1.0e-9 * frequency_hz
            and base_key == "control"
            and not _knobs_set()
        ):
            label, element = next(
                (lab, el) for lab, el in STEP1B_TERMINATIONS if lab.startswith("C")
            )
            z_term = series_rlc_impedance(frequency_hz, **element)
            measured, _kept, _res = _terminated_three_port(variant_a, z_term)
            predicted = reduce_terminated_ports(
                variant_a["s"],
                float(REFERENCE_IMPEDANCE_OHM),
                {TERMINATED_PORT_INDEX: z_term},
            )
            record["a"]["registered_residual"] = {
                "label": label,
                "value": float(np.linalg.norm(measured - predicted) / np.linalg.norm(predicted)),
            }
        del variant_a

    # --- variant (c'): the port-box conductor-side grading halved, measured on
    # the C4-congruent cut against its own congruent control' (step 1b, the
    # 2026-09-21 review ruling).  Both sides of the difference re-mesh together.
    if need_prime and control_prime is None:
        control_prime, control_prime_means, t_cp0 = _control_prime()
    control_prime = None

    if "c_prime" in variants:
        comm.Barrier()
        t2 = time.perf_counter()
        variant_c = build_four_port_sweep(
            frequency_hz=frequency_hz,
            c4_congruent_sheets=True,
            conductor_resolution=0.5 * CONDUCTOR_RESOLUTION,
            resolution=resolution,
            degree=degree,
        )
        comm.Barrier()
        t_c = time.perf_counter() - t2

        record["c_prime"] = _gate_readings(variant_c)
        record["c_prime"]["seconds"] = float(t_c)
        record["c_prime"]["conductor_resolution"] = float(0.5 * CONDUCTOR_RESOLUTION)
        c_means = _class_means(variant_c["s"])
        record["c_prime"]["class_means"] = _means_json(c_means)
        record["c_prime"]["s_diagonal"] = _s_diagonal(variant_c["s"])
        record["c_prime"]["move"] = _relative_move(c_means, control_prime_means)
        del variant_c
        record["rss"]["c_prime"] = _rss_report(comm, "(c') sweep")

    if comm.rank == 0:
        cells = {k: record[k]["cells"] for k in ALL_VARIANTS if k in record}
        print(
            f"\n[PORT-21 step1] f = {frequency_hz:.6e} Hz  -n {comm.size}: "
            f"degree {degree}, resolution "
            f"{'module default' if resolution is None else f'{resolution:.4g} m'}, "
            f"runs {list(variants)}; "
            f"control {t_control:.2f} s, "
            f"(a) kappa route {t_a:.2f} s (kappa = {kappa:.9e}), "
            f"control' (C4-congruent cut) {t_cp0:.2f} s (cell record "
            f"{CONTROL_PRIME_CELL_RECORD}, printed not asserted), "
            f"(c') halved conductor grading on the congruent cut {t_c:.2f} s "
            f"(h_c = {0.5 * CONDUCTOR_RESOLUTION:.4e} m vs "
            f"{CONDUCTOR_RESOLUTION:.4e} m); cells {cells}\n"
            f"    SKIPPED, not built: {SKIPPED_VARIANT_B}",
            flush=True,
        )
        _READINGS_DIR.mkdir(parents=True, exist_ok=True)
        path = _readings_path(frequency_hz, degree, resolution)
        path.write_text(json.dumps(record, indent=2))
        print(f"[PORT-21 step1] readings persisted to {path}", flush=True)
    return record


@pytest.fixture(scope="module")
def port21_runs():
    frequencies = _frequencies_hz()
    if not frequencies:
        pytest.skip(f"{FREQUENCIES_ENV} unset — `PORT-21` step 1 runs only when a slot asks")
    return [_run_one_frequency(f) for f in frequencies]


@complex_only
def test_the_control_reproduces_the_gate_record(port21_runs):
    """**The control.**  Nothing else is read if this fails (the step's own rule).

    At 10 MHz the control column must land on leg (d0)'s registered ``Z``
    column at the imported ``LEG_D0_REPRODUCTION_BAND`` (1e-9); at every
    frequency the control must mesh the imported ``STEP1_CELL_RECORD`` cells,
    or the variants below are not perturbations of the gate fixture at all.
    """
    comm = MPI.COMM_WORLD
    if _knobs_set() or "control" not in _variants():
        pytest.skip(
            "step 3 knob set (degree/resolution) or `control` not selected: "
            "leg (d0)'s Z record and STEP1_CELL_RECORD belong to the degree-1 "
            "gate mesh only, so this control is not the reference here"
        )
    for record in port21_runs:
        cells = record["control"]["cells"]
        assert cells == STEP1_CELL_RECORD, (
            f"the control at f = {record['frequency_hz']:.3e} Hz meshed {cells} "
            f"cells against the gate record {STEP1_CELL_RECORD} — this is not "
            "the gate fixture and no variant below is readable"
        )

    at_record = [
        r
        for r in port21_runs
        if abs(r["frequency_hz"] - GATE_RECORD_FREQUENCY_HZ) <= 1.0e-9 * GATE_RECORD_FREQUENCY_HZ
    ]
    if not at_record:
        pytest.skip(
            "this window carries no 10 MHz rung, so leg (d0)'s registered Z "
            "column is not reproducible here; the cell-count control above ran"
        )
    column = np.array(
        [complex(re, im) for re, im in at_record[0]["control"]["z_column1"]],
        dtype=np.complex128,
    )
    if comm.rank == 0:
        print(
            f"\n[PORT-21 step1] CONTROL: the P1-driven column vs leg (d0)'s "
            f"record, imported band {LEG_D0_REPRODUCTION_BAND:.0e}:",
            flush=True,
        )
        for k, (z, z_rec) in enumerate(zip(column, LEG_D0_Z_COLUMN), start=1):
            print(
                f"    Z_{k}1 {z:+.9e} vs record {z_rec:+.9e}  rel "
                f"{abs(z - z_rec) / abs(z_rec):.3e}",
                flush=True,
            )
    for k, (z, z_rec) in enumerate(zip(column, LEG_D0_Z_COLUMN), start=1):
        err = abs(z - z_rec) / abs(z_rec)
        assert err < LEG_D0_REPRODUCTION_BAND, (
            f"Z_{k}1 deviates {err:.3e} from leg (d0)'s record against the "
            f"imported {LEG_D0_REPRODUCTION_BAND:.0e} band — record drift on "
            "`main`; no variant in this table is read (known-issues, stop)"
        )


@complex_only
def test_every_variant_still_passes_the_imported_identity_gates(port21_runs):
    """**Anchor (asserted).**  Reciprocity, passivity and C4 hold under each knob.

    Every band is imported from `test_port_birdcage_four_port`; none is
    restated and none is moved.  A red row here would be this row's printed
    finding about that convention, not a defect to fix — but it is still
    asserted, because the step's anchor says every variant run passes.
    """
    comm = MPI.COMM_WORLD
    # Collected, not short-circuited: the step wants the *whole* table printed
    # even when a row is red, because a red row is a printed finding about that
    # convention.  Every band below is still the imported one and every failure
    # below is still a failure — only the reporting order changes.
    misses = []
    for record in port21_runs:
        for key in [k for k in ALL_VARIANTS if k in record]:
            g = record[key]
            if comm.rank == 0:
                print(
                    f"\n[PORT-21 step1] gates, f = {record['frequency_hz']:.3e} Hz, "
                    f"variant {key} ({g['cells']} cells):\n"
                    f"    ||S - S^T||/||S|| = {g['reciprocity']:.9e} "
                    f"(band {RECIPROCITY_BAND:.0e})\n"
                    f"    sigma_max = {g['sigma_max']:.9f}, max column power "
                    f"{g['column_power_max']:.9f} "
                    f"(tol {PASSIVITY_SIGMA_TOLERANCE:.0e})\n"
                    f"    class spreads: self {g['spreads']['self'] * 100:.4f}%  "
                    f"adjacent {g['spreads']['adjacent'] * 100:.4f}%  "
                    f"opposite {g['spreads']['opposite'] * 100:.4f}%  "
                    f"(band {ADJACENT_SPREAD_BAND * 100:.1f}%)",
                    flush=True,
                )
                if key == "control_prime":
                    rec = CONTROL_PRIME_SPREAD_RECORDS.get(record["frequency_hz"])
                    if rec:
                        beside = (
                            f"{rec[0][0]:.4f} / {rec[0][1]:.4f} / " f"{rec[0][2]:.4f}% ({rec[1]})"
                        )
                    else:
                        beside = (
                            "NONE - no congruent-cut record exists at this "
                            "frequency; this reading is new"
                        )
                    print(
                        "    control' spreads on record (printed beside, not "
                        f"asserted; different rank width): {beside}; cells on "
                        f"record {CONTROL_PRIME_CELL_RECORD}",
                        flush=True,
                    )
            tag = f"variant {key} at {record['frequency_hz']:.3e} Hz"
            if g["reciprocity"] > RECIPROCITY_BAND:
                misses.append(
                    f"{tag}: ||S - S^T||/||S|| = {g['reciprocity']:.9e} against "
                    f"the imported {RECIPROCITY_BAND:.0e}"
                )
            if g["sigma_max"] > 1.0 + PASSIVITY_SIGMA_TOLERANCE:
                misses.append(f"{tag}: sigma_max = {g['sigma_max']:.9f} — the network is active")
            if g["column_power_max"] > 1.0 + PASSIVITY_SIGMA_TOLERANCE:
                misses.append(
                    f"{tag}: column power {g['column_power_max']:.9f} > 1 — it "
                    "returns more power than it is fed"
                )
            for name, spread in g["spreads"].items():
                if spread > ADJACENT_SPREAD_BAND:
                    misses.append(
                        f"{tag}: the {name} class spreads {spread * 100:.4f}% "
                        f"against the imported {ADJACENT_SPREAD_BAND * 100:.1f}%"
                    )
    assert not misses, (
        "imported identity gates missed under a port-model convention (bands "
        "imported and unmoved; a red row is a finding about that convention, "
        "never a licence to widen):\n  " + "\n  ".join(misses)
    )


@complex_only
def test_variant_a_reproduces_the_port14_step3_registered_residual(port21_runs):
    """**Anchor (asserted).**  (a)'s terminated ``C`` residual <= ``REDUCTION_BAND``.

    `PORT-14` step 3's registered reading, at its registered frequency, on the
    identical corrected route: the terminated 3x3 measured in the field against
    the circuit reduction of the corrected 4x4.  Band imported, not restated.
    """
    comm = MPI.COMM_WORLD
    rows = [r for r in port21_runs if "registered_residual" in r.get("a", {})]
    if not rows:
        pytest.skip(
            "this window carries no "
            f"{STEP3_REGISTERED_FREQUENCY_HZ:.3e} Hz rung, so `PORT-14` step 3's "
            "registered residual is not reproducible here"
        )
    for record in rows:
        entry = record["a"]["registered_residual"]
        if comm.rank == 0:
            print(
                f"\n[PORT-21 step1] anchor (a), f = {record['frequency_hz']:.3e} Hz: "
                f"`PORT-14` step 3's registered residual, termination "
                f"{entry['label']}: {entry['value']:.6e} "
                f"({entry['value'] / REDUCTION_BAND:.3f}x the imported band "
                f"{REDUCTION_BAND:.0e})   kappa in-run = {record['a']['kappa']:.9e}",
                flush=True,
            )
        assert entry["value"] <= REDUCTION_BAND, (
            f"variant (a) reads residual {entry['value']:.6e} against the "
            f"imported REDUCTION_BAND {REDUCTION_BAND:.0e} — the corrected "
            "route here is not `PORT-14` step 3's"
        )


@complex_only
def test_each_variant_reached_the_solve(port21_runs):
    """**Negative control (asserted).**  Each knob moves the self class > 1e-6.

    The `OPS-41` knob-reaches-the-solve pattern: a variant that moves nothing
    was dropped somewhere between the keyword and the assembled form, and its
    table row would read as "this convention does not matter" when in fact it
    was never applied.  The floor is a wiring check with >= 1e4 of headroom
    against the O(1e-2) moves `PORT-14` step 3's kappa ~ 1 % implies.
    """
    comm = MPI.COMM_WORLD
    for record in port21_runs:
        for key in [k for k in ("a", "c_prime") if k in record]:
            move = record[key]["move"]["self"]
            if comm.rank == 0:
                print(
                    f"[PORT-21 step1] knob-reaches-the-solve, f = "
                    f"{record['frequency_hz']:.3e} Hz, variant {key}: self-class "
                    f"move {move:.6e} against the {KNOB_REACHES_THE_SOLVE:.0e} "
                    f"floor ({move / KNOB_REACHES_THE_SOLVE:.3e}x)",
                    flush=True,
                )
            assert move > KNOB_REACHES_THE_SOLVE, (
                f"variant {key} at {record['frequency_hz']:.3e} Hz moved the "
                f"self class by {move:.6e} <= {KNOB_REACHES_THE_SOLVE:.0e} — the "
                "knob did not reach the solve, so its table row is not a "
                "measurement of that convention"
            )


@complex_only
def test_the_sensitivity_table_is_printed(port21_runs):
    """**Printed, never asserted.**  ``|dS|/|S|`` per class per variant per frequency.

    Reads every per-frequency JSON this slot's windows persisted, so the last
    window prints the whole table even though each window solved part of it.
    The frequency trend of each variant's self-class move — whether it *falls*
    with frequency — is printed as a word, and asserted nowhere: it is the
    reading a later private step takes against the residual, not a gate.
    """
    comm = MPI.COMM_WORLD
    if comm.rank != 0:
        return
    rows = []
    tag = _readings_path(0.0, _degree(), _resolution()).name.split("MHz_", 1)[1]
    for path in sorted(_READINGS_DIR.glob(f"f_*MHz_{tag}")):
        rows.append(json.loads(path.read_text()))
    rows.sort(key=lambda r: r["frequency_hz"])
    print(
        "\n[PORT-21 step1] SENSITIVITY TABLE — |dS|/|S| per C4 class, variant "
        "against the control on the same fixture (PRINTED, NEVER ASSERTED; no "
        "absolute S11 claim is licensed here — known-issues 2026-09-19):",
        flush=True,
    )
    print(
        f"    {'f (MHz)':>9s}  {'variant':<8s}  {'self':>13s}  {'adjacent':>13s}  "
        f"{'opposite':>13s}   cells",
        flush=True,
    )
    for record in rows:
        for key in ("a", "c_prime"):
            if key not in record:
                continue
            m = record[key]["move"]
            print(
                f"    {record['frequency_hz'] / 1e6:9.3f}  {key:<8s}  "
                f"{m['self']:13.6e}  {m['adjacent']:13.6e}  "
                f"{m['opposite']:13.6e}   {record[key]['cells']}",
                flush=True,
            )
    for key in ("a", "c_prime"):
        series = [(r["frequency_hz"] / 1e6, r[key]["move"]["self"]) for r in rows if key in r]
        if len(series) < 2:
            print(
                f"    variant {key}: frequency trend needs >= 2 rungs; this run "
                f"has {len(series)}",
                flush=True,
            )
            continue
        falls = all(b[1] < a[1] for a, b in zip(series, series[1:]))
        rises = all(b[1] > a[1] for a, b in zip(series, series[1:]))
        trend = (
            "FALLS with frequency"
            if falls
            else ("RISES with frequency" if rises else "is non-monotonic in frequency")
        )
        print(
            f"    variant {key}: the self-class move {trend} — "
            + ", ".join(f"{f:.0f} MHz {v:.3e}" for f, v in series),
            flush=True,
        )
    print(
        "    (c') is read against control' (the C4-congruent cut), never against "
        "the gate-record control; (a) against the gate-record control.",
        flush=True,
    )
    print(
        "    RECORD, not re-run - variant (c) as first built on the UNCONGRUENT "
        f"cut ({UNCONGRUENT_C_CELLS} cells), replaced by (c'); its class spread "
        "exceeds the move it was to measure (mesh asymmetry, 2026-09-21 ruling):",
        flush=True,
    )
    for f_hz, (mv, sp, log) in sorted(UNCONGRUENT_C_RECORD.items()):
        print(
            f"    {f_hz / 1e6:9.3f}  c(rec)    {mv[0]:13.6e}  {mv[1]:13.6e}  "
            f"{mv[2]:13.6e}   spreads {sp[0]:.4f}/{sp[1]:.4f}/{sp[2]:.4f}%  "
            f"[{log}]",
            flush=True,
        )
    print(f"    SKIPPED, not built: {SKIPPED_VARIANT_B}", flush=True)


@complex_only
def test_step1_is_reproduced_with_the_knobs_unset(port21_runs):
    """**Step 3 anchor (asserted).**  Unset knobs reproduce step 1's (a) moves.

    (a)'s self-class move against the gate-record control, to the digit step 1
    printed (`20260921T094335Z_PORT-21.log:11532`, `:11534`).  A mismatch means
    the step-3 plumbing moved step 1 (the item's negative result: revert).
    """
    comm = MPI.COMM_WORLD
    if _knobs_set() or not {"control", "a"} <= set(_variants()):
        pytest.skip("a step-3 knob is set or (a)/control not selected")
    rows = [r for r in port21_runs if r["frequency_hz"] in STEP1_A_SELF_MOVE_PRINTED]
    if not rows:
        pytest.skip("no 10 or 64 MHz rung in this window")
    for record in rows:
        printed = f"{record['a']['move']['self']:.6e}"
        expected = STEP1_A_SELF_MOVE_PRINTED[record["frequency_hz"]]
        if comm.rank == 0:
            print(
                f"[PORT-21 step3] W1 anchor, f = {record['frequency_hz']:.3e} Hz: "
                f"(a) self-class move {printed} vs step 1's printed {expected}",
                flush=True,
            )
        assert printed == expected, (
            f"(a)'s self move {printed} != step 1's {expected}: the step-3 " "plumbing moved step 1"
        )


@complex_only
def test_the_step3_knobs_reach_the_solve(port21_runs):
    """**Step 3 negative control (asserted; the `ANS-6` precedent).**

    Degree knob: the control-prime class means at this degree move by more than
    ``KNOB_REACHES_THE_SOLVE`` from the degree-1 reading on the same
    resolution, read from that run's tagged JSON (W1 writes it at 10 MHz).
    Resolution knob: control' (and the control, if run) meshes a cell count
    different from its degree-1 record.  Every move and kappa is printed.
    """
    comm = MPI.COMM_WORLD
    degree, resolution = _degree(), _resolution()
    if not _knobs_set():
        pytest.skip("no step-3 knob is set; W1's anchors carry this window")
    checked = 0
    for record in port21_runs:
        f_hz = record["frequency_hz"]
        if resolution is not None:
            for key, rec in (
                ("control_prime", CONTROL_PRIME_CELL_RECORD),
                ("control", STEP1_CELL_RECORD),
            ):
                if key not in record:
                    continue
                cells = record[key]["cells"]
                if comm.rank == 0:
                    print(
                        f"[PORT-21 step3] resolution knob {resolution:.4g} m, "
                        f"f = {f_hz:.3e} Hz, {key}: {cells} cells vs the "
                        f"default-resolution record {rec}",
                        flush=True,
                    )
                assert cells != rec, (
                    f"{key} meshed the default-resolution record {rec} cells "
                    f"under {RESOLUTION_ENV}={resolution}: the knob did not "
                    "reach the mesher"
                )
                checked += 1
        if degree != 1 and "control_prime" in record:
            ref_path = _readings_path(f_hz, 1, resolution)
            if not ref_path.exists():
                pytest.fail(
                    f"no degree-1 reading at {ref_path} to compare control' "
                    "against: run W1 (knobs unset) at this frequency first"
                )
            ref = json.loads(ref_path.read_text())
            ref_means = {k: complex(*v) for k, v in ref["control_prime"]["class_means"].items()}
            means = {k: complex(*v) for k, v in record["control_prime"]["class_means"].items()}
            move = _relative_move(means, ref_means)
            if comm.rank == 0:
                print(
                    f"[PORT-21 step3] degree knob, f = {f_hz:.3e} Hz: control' "
                    f"degree {degree} vs degree 1 ({ref_path.name}, "
                    f"{ref['control_prime']['cells']} cells): self "
                    f"{move['self']:.6e}  adjacent {move['adjacent']:.6e}  "
                    f"opposite {move['opposite']:.6e}; kappa at degree "
                    f"{degree} = {record.get('a', {}).get('kappa', float('nan')):.9e}"
                    f" (degree 1: {ref.get('a', {}).get('kappa', float('nan')):.9e})",
                    flush=True,
                )
                if "a" in record:
                    m = record["a"]["move"]
                    print(
                        f"[PORT-21 step3] degree {degree}, f = {f_hz:.3e} Hz: (a) "
                        f"against {record['a']['base']}: self {m['self']:.6e}  "
                        f"adjacent {m['adjacent']:.6e}  opposite "
                        f"{m['opposite']:.6e} (PRINTED; no degree-2 PORT-21 "
                        "measurement existed before this)",
                        flush=True,
                    )
            assert move["self"] > KNOB_REACHES_THE_SOLVE, (
                f"control' self class moved {move['self']:.6e} <= "
                f"{KNOB_REACHES_THE_SOLVE:.0e} from degree 1 to degree {degree}: "
                "the degree knob did not reach the solve"
            )
            checked += 1
    assert checked > 0, "a knob is set but no run carried a checkable reading"


@complex_only
def test_degree2_10mhz_reproduces_w3(port21_runs):
    """**`PORT-21` step 4 anchor (asserted).**  F6's W3, re-run, to 1e-9.

    Active only at degree 2, resolution unset, a 10 MHz rung, and both
    control' and (a) selected.  The in-run ``S`` diagonals of control' and (a)
    and the in-run kappa must land on W3's printed records within the imported
    ``LEG_D0_REPRODUCTION_BAND`` (1e-9, relative per entry).  Compared against
    the pinned constants, never the JSON (this window overwrites W3's file).

    **Negative control (asserted; backed by W3, `:2414` vs `:2932`).**  The
    in-run (a) diagonal read against W3's *control'* diagonal must MISS the
    band (measured gap ~ 0.185 relative): the comparison can fail.
    """
    comm = MPI.COMM_WORLD
    if _degree() != 2 or _resolution() is not None:
        pytest.skip("W3 reproduction is degree 2 at the module resolution only")
    if not {"control_prime", "a"} <= set(_variants()):
        pytest.skip("W3 reproduction needs control' and (a) both selected")
    rows = [
        r
        for r in port21_runs
        if abs(r["frequency_hz"] - W3_FREQUENCY_HZ) <= 1.0e-9 * W3_FREQUENCY_HZ
    ]
    if not rows:
        pytest.skip("no 10 MHz rung in this window; W3 is not reproducible here")
    record = rows[0]

    def _rel(measured, pinned):
        return [abs(m - p) / abs(p) for m, p in zip(measured, pinned)]

    cp = [complex(*v) for v in record["control_prime"]["s_diagonal"]]
    a = [complex(*v) for v in record["a"]["s_diagonal"]]
    kappa = float(record["a"]["kappa"])
    err_cp = _rel(cp, W3_CONTROL_PRIME_S_DIAGONAL)
    err_a = _rel(a, W3_A_S_DIAGONAL)
    err_k = abs(kappa - W3_KAPPA) / abs(W3_KAPPA)
    neg = _rel(a, W3_CONTROL_PRIME_S_DIAGONAL)
    if comm.rank == 0:
        print(
            f"\n[PORT-21 step4] W_R: degree-2 10 MHz vs {W3_LOG}, imported band "
            f"{LEG_D0_REPRODUCTION_BAND:.0e} (-n {comm.size}):",
            flush=True,
        )
        for k, (z, p, e) in enumerate(zip(cp, W3_CONTROL_PRIME_S_DIAGONAL, err_cp), 1):
            print(
                f"    control' S_{k}{k} {z:+.9e} vs W3 {p:+.9e}  rel {e:.3e}",
                flush=True,
            )
        for k, (z, p, e) in enumerate(zip(a, W3_A_S_DIAGONAL, err_a), 1):
            print(
                f"    (a)      S_{k}{k} {z:+.9e} vs W3 {p:+.9e}  rel {e:.3e}",
                flush=True,
            )
        print(
            f"    kappa {kappa:.9e} vs W3 {W3_KAPPA:.9e}  rel {err_k:.3e}\n"
            f"    max over the 8 diagonals {max(err_cp + err_a):.3e}\n"
            f"    NEGATIVE CONTROL: (a) vs W3's control' diagonal, min rel "
            f"{min(neg):.6e} / max rel {max(neg):.6e} (must miss "
            f"{LEG_D0_REPRODUCTION_BAND:.0e})",
            flush=True,
        )
    assert min(neg) > LEG_D0_REPRODUCTION_BAND, (
        f"negative control: (a) lands on W3's control' diagonal ({min(neg):.3e}); "
        "the comparison cannot discriminate"
    )
    for name, errs in (("control'", err_cp), ("(a)", err_a)):
        assert max(errs) < LEG_D0_REPRODUCTION_BAND, (
            f"{name} S diagonal misses W3 by {max(errs):.3e} against "
            f"{LEG_D0_REPRODUCTION_BAND:.0e}: the degree-2 solve moved"
        )
    assert (
        err_k < LEG_D0_REPRODUCTION_BAND
    ), f"kappa {kappa:.9e} misses W3's {W3_KAPPA:.9e} by {err_k:.3e}"
