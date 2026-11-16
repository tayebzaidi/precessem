"""Example (`EX-152`): the lumped-sheet power identity per sheet, and the sheet
excess against the sheet width, on the F-small exact-rectangle four-port.

`ANS-4` H33a (gated, 2026-10-25) splits each lumped sheet's field dissipation
into the part its terminal accounts for, ``1/2 Re(1/Z_p)|V_s|^2``, and an
**excess** ``X_s = P_sheet,field - 1/2 Re(1/Z_p)|V_s|^2``. ``CS_s`` is the same
excess built from the sheet's facet integrals (area, ``int E.h``, ``int |n x E|^2``,
``int |E.h|^2``) by a second code path, so ``CS_s / X_s = 1`` compares two
independent computations. This example re-computes H33a's gates at two of its
three widths, ``f = 0.25`` (w0/2) and ``f = 0.5`` (w0, the record width), at
degree 1, ``-n 2``, 64 MHz, P1 driven alone.

**Imported, never restated** (the `ANS-1` rule), from
``tests.validation.test_ans4_h33_sheet_width_excess`` as ``h33``: ``_width_case``,
``excess_w``, ``gap_voltage``, ``p_acc_wave``, ``width_exponent``, ``PORT_IDS``,
``DRIVEN``, ``PLANTED_V_FACTOR``, ``RECORD_FRACTION``; the bands from the modules
``h33`` takes them from: ``DISCRETE_IDENTITY_RTOL``, ``STRIP_BAND``,
``BLIND_FACTOR``, ``LEG_D0_REPRODUCTION_BAND``, ``ANS4_EXACT_STRIP_S_P1_RECORDS``.
No ``test_*`` name is imported; no gate or ``src/`` edit.

**Asserted** (both widths unless noted): (Wd) sheet area = f * port_dx * gap length
within ``STRIP_BAND``; (CS) ``|CS_s/X_s - 1| <= DISCRETE_IDENTITY_RTOL`` on all four
sheets, with (c-V) the planted x1.01 on V missing by >= ``BLIND_FACTOR`` *
``DISCRETE_IDENTITY_RTOL``; (P_w) ``P_src = P_vol + sum P_sheet,field`` and
(P_acc-1) line cut = wave form, each within ``DISCRETE_IDENTITY_RTOL``, with
(c-blind) the conductor-blind sum missing by >= the same floor; (R0) at
f = ``RECORD_FRACTION`` the rung's S equals ``ANS4_EXACT_STRIP_S_P1_RECORDS["64"]``
within ``LEG_D0_REPRODUCTION_BAND`` (``-n 2`` only).

**Printed only, never gated:** X_s / P_acc per sheet, sum X_s / P_acc per width
beside H33a's records, the two-point width exponent (a two-point slope), and
sigma_max(S) per width. Record digits are printed beside readings, never asserted.

**Not shown:** degree 2 (H33b), the w0/4 width, Palace's excess, any cause of the
excess.

Negative result protocol: a gate red here while H33a's W1 log was green is an
image or path divergence -- known-issues entry, stop, never a re-record.

Needs the complex DolfinX build; the runner sources it for the ``ports:`` group::

    ./run_examples.sh -e ports:23 -n 2 -t 600

Output (``paraview_output/``, gitignored): the record-width P1 drive's ``|E|`` as
``ports_23_lumped_sheet_excess_vs_width_combined.xdmf`` with the four sheet facet
tags, and the per-sheet table ``ports_23_lumped_sheet_excess_vs_width_table.csv``.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

import ufl
from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.validation.test_ans4_h33_sheet_width_excess as h33  # noqa: E402
from tests.mesh.test_birdcage_leg_gaps import LEG_GAP_LENGTH  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.mesh.test_birdcage_port_tags import LEG_COUNT  # noqa: E402
from tests.validation.test_ans4_exact_strip_records import (  # noqa: E402
    ANS4_EXACT_STRIP_RECORD_RANK_WIDTH,
    ANS4_EXACT_STRIP_S_P1_RECORDS,
)
from tests.validation.test_birdcage_b1_plus_map import _solve_driven  # noqa: E402
from tests.validation.test_birdcage_power_identity import (  # noqa: E402
    DISCRETE_IDENTITY_RTOL,
)
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND  # noqa: E402
from tests.validation.test_port_birdcage_four_port import (  # noqa: E402
    LEG_D0_REPRODUCTION_BAND,
)
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_wf12_f_human_copper_matrix import BLIND_FACTOR  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
BASENAME = "ports_23_lumped_sheet_excess_vs_width"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
SETUP_BASENAME = "ports_23_lumped_sheet_excess_vs_width_setup"

FRACTIONS = (0.25, 0.5)
DEGREE = 1
FREQ_KEY = "64"

# H33a records (`20261025T095010Z_ANS-4-H33a-W1.log:5660-5665`), printed beside
# the readings -- never asserted for equality.
RECORD_SUM_X = {0.25: 9.957736e-05, 0.5: 4.135372e-05}
RECORD_EXPONENT_3PT = -1.6848
RECORD_R0 = 9.557e-11


def _e_mag_dg0(e_complex, msh, name):
    """``|E|`` as a real DG0 field."""
    dg0 = fem.functionspace(msh, ("DG", 0))
    tmp = fem.Function(dg0)
    tmp.interpolate(
        fem.Expression(
            ufl.inner(e_complex, ufl.conj(e_complex)), dg0.element.interpolation_points
        )
    )
    out = fem.Function(dg0, name=name)
    out.x.array[:] = np.sqrt(np.maximum(np.real(tmp.x.array), 0.0))
    out.x.scatter_forward()
    return out


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()
    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "This example needs the complex DolfinX build: "
            "source /usr/local/bin/dolfinx-complex-mode (the runner does this "
            "automatically for the `ports:` group)."
        )
    if comm.rank == 0:
        print("=" * 78, flush=True)
        print(
            "EX-152 -- lumped-sheet identity per sheet and excess vs width, "
            "F-small exact-rectangle four-port, 64 MHz, degree 1",
            flush=True,
        )
        print("=" * 78, flush=True)

    cases = {}
    for f in FRACTIONS:
        cases[f] = h33._width_case(f, DEGREE, comm)
        if comm.rank == 0:
            c = cases[f]
            print(
                f"\n[fixture] f = {f}: {c['rung']['cells']} cells, {c['time_s']:.1f} s "
                f"(rung + P1 re-solve + diagnostic) at -n {comm.size}",
                flush=True,
            )

    rec_case = cases[h33.RECORD_FRACTION]
    write_setup_figure(
        rec_case["rung"]["mesh"],
        rec_case["rung"]["cell_tags"],
        FIGURE_DIR / f"{SETUP_BASENAME}.png",
        region_names={
            CONDUCTOR_CELL_TAG: "conductor",
            2: "air",
            PHANTOM_CELL_TAG: "phantom",
            **{PORT_LOWER + i: f"leg L{i} port lower" for i in range(1, LEG_COUNT + 1)},
            **{PORT_UPPER + i: f"leg L{i} port upper" for i in range(1, LEG_COUNT + 1)},
        },
        hide_tags=(2,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="ports:23 -- F-small four-port, exact-rectangle sheets (ANS-4 H33a)",
        comm=comm,
    )

    # ------------------------------------------------------------- the gates
    readings = {}
    failures = []
    floor = BLIND_FACTOR * DISCRETE_IDENTITY_RTOL
    for f, c in cases.items():
        dx = c["rung"]["sheet_analytic"] / LEG_GAP_LENGTH
        target = f * dx * LEG_GAP_LENGTH
        wd = [abs(s["area"] / target - 1.0) for s in c["rung"]["sheets"]]
        cs = {pid: abs(r["cs"] / r["x"] - 1.0) for pid, r in c["per"].items()}
        cv = {pid: abs(r["cs"] / r["x_bad"] - 1.0) for pid, r in c["per"].items()}
        ex = c["exact"]
        rel_w = abs(ex["p_src"] - ex["p_vol"] - ex["sheet_field_total"]) / abs(ex["p_src"])
        blind = abs(ex["p_src"] - ex["phantom"] - ex["sheet_field_total"]) / abs(ex["p_src"])
        rel_acc = abs(c["p_acc_line"] - c["p_acc_wave"]) / abs(c["p_acc_wave"])
        xs = {pid: r["x"] / c["p_acc_wave"] for pid, r in c["per"].items()}
        readings[f] = {
            "wd": wd, "cs": cs, "cv": cv, "rel_w": rel_w, "blind": blind,
            "rel_acc": rel_acc, "xs": xs, "sum_x": sum(xs.values()),
        }
        if max(wd) > STRIP_BAND:
            failures.append(f"(Wd) f={f}: {max(wd):.3e} > {STRIP_BAND:g}")
        for pid in h33.PORT_IDS:
            if cs[pid] > DISCRETE_IDENTITY_RTOL:
                failures.append(f"(CS) f={f} {pid}: {cs[pid]:.3e}")
            if cv[pid] < floor:
                failures.append(f"(c-V) f={f} {pid}: miss {cv[pid]:.3e} < {floor:g}")
        if rel_w > DISCRETE_IDENTITY_RTOL:
            failures.append(f"(P_w) f={f}: {rel_w:.3e}")
        if blind < floor:
            failures.append(f"(c-blind) f={f}: miss {blind:.3e} < {floor:g}")
        if rel_acc > DISCRETE_IDENTITY_RTOL:
            failures.append(f"(P_acc-1) f={f}: {rel_acc:.3e}")

    s_rec = np.asarray(ANS4_EXACT_STRIP_S_P1_RECORDS[FREQ_KEY], dtype=np.complex128)
    s_now = rec_case["rung"]["s"]
    r0 = float(np.linalg.norm(s_now - s_rec) / np.linalg.norm(s_rec))
    r0_applies = comm.size == ANS4_EXACT_STRIP_RECORD_RANK_WIDTH
    if r0_applies and r0 > LEG_D0_REPRODUCTION_BAND:
        failures.append(f"(R0) f={h33.RECORD_FRACTION}: {r0:.3e} > {LEG_D0_REPRODUCTION_BAND:g}")

    sums = [readings[f]["sum_x"] for f in FRACTIONS]
    expo = h33.width_exponent(list(FRACTIONS), sums)

    if comm.rank == 0:
        for f in FRACTIONS:
            r = readings[f]
            print(f"\n[f = {f}]", flush=True)
            print(
                f"  (Wd) max sheet-area rel {max(r['wd']):.3e} (band {STRIP_BAND:g})",
                flush=True,
            )
            for pid in h33.PORT_IDS:
                print(
                    f"  (CS) {pid}: |CS/X - 1| {r['cs'][pid]:.3e} (band "
                    f"{DISCRETE_IDENTITY_RTOL:g}); (c-V) x{h33.PLANTED_V_FACTOR} on V "
                    f"misses by {r['cv'][pid]:.4e} (floor {floor:g})",
                    flush=True,
                )
            print(
                f"  (P_w) rel {r['rel_w']:.3e}; (c-blind) conductor-blind miss "
                f"{r['blind']:.4e}; (P_acc-1) rel {r['rel_acc']:.3e}",
                flush=True,
            )
            print(
                "  PRINTED X_s/P_acc: "
                + " ".join(f"{pid} {v:+.6e}" for pid, v in r["xs"].items())
                + f"; sum {r['sum_x']:+.6e} (H33a record {RECORD_SUM_X[f]:.6e}); "
                f"sigma_max(S) {float(cases[f]['rung']['sigma'][0]):.9f}",
                flush=True,
            )
        print(
            f"\n(R0) f={h33.RECORD_FRACTION} ||S - S_rec||/||S_rec|| {r0:.3e} "
            f"(band {LEG_D0_REPRODUCTION_BAND:g}; H33a record {RECORD_R0:.3e}; "
            f"asserted at -n {ANS4_EXACT_STRIP_RECORD_RANK_WIDTH}: {r0_applies})",
            flush=True,
        )
        print(
            f"PRINTED two-point width exponent of sum X_s/P_acc (f = {FRACTIONS[0]}, "
            f"{FRACTIONS[1]}; a two-point slope): {expo:+.4f} (H33a three-point "
            f"{RECORD_EXPONENT_3PT:+.4f}; H33b degree-2 approx -2.09)",
            flush=True,
        )

    assert not failures, "example/test divergence -- known-issues entry, stop: " + "; ".join(
        failures
    )

    # ------------------------------------------------- the picture and table
    sweep = {
        "mesh": rec_case["rung"]["mesh"],
        "cell_tags": rec_case["rung"]["cell_tags"],
        "problem": rec_case["rung"]["problem"],
        "facet_tags": rec_case["rung"]["facet_tags"],
        "specs": rec_case["rung"]["specs"],
        "degree": DEGREE,
    }
    solved = _solve_driven(sweep, h33.DRIVEN)
    fields = {
        "E_mag_P1": _e_mag_dg0(solved["fields"].e_complex, sweep["mesh"], "E_mag_P1")
    }
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        sweep["mesh"],
        sweep["cell_tags"],
        fields,
        comm=comm,
        facet_tags=sweep["facet_tags"],
    )
    if comm.rank == 0:
        csv_path = OUTPUT_DIR / f"{BASENAME}_table.csv"
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["f", "port", "area_m2", "CS_W", "X_W", "CS_over_X",
                        "X_over_P_acc", "P_acc_W"])
            for f in FRACTIONS:
                c = cases[f]
                for pid in h33.PORT_IDS:
                    r = c["per"][pid]
                    w.writerow([f, pid, f"{r['area']:.9e}", f"{r['cs']:.9e}",
                                f"{r['x']:.9e}", f"{r['cs'] / r['x']:.12f}",
                                f"{r['x'] / c['p_acc_wave']:.9e}",
                                f"{c['p_acc_wave']:.9e}"])
    comm.Barrier()
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    if comm.rank == 0:
        print(
            f"\n[paraview] wrote {path} (|E| of the f = {h33.RECORD_FRACTION} P1 drive, "
            f"DG0, with the sheet facet tags) and {csv_path.name}",
            flush=True,
        )
        print(
            f"\nAll gates hold. Total elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
