"""Example (`EX-150`): the power-wave identity per drive on the 4-leg birdcage.

`POST-6` step 1b (gated, 2026-09-23) closes a question the drive-level print of
`ports:13` leaves open. For a single drive ``k`` of the assembled 4x4 S-matrix,
the power the package says was *accepted*,

    P_acc,k = 1/2 |a_k|^2 (1 - sum_i |S_ik|^2),

equals the power the terminal accounting says was *delivered*,

    P_net,k = supplied_k - sum_i sheets_ik,

with ``supplied_k = 1/2 Re(V_src I_k*)`` and ``sheets_ik = 1/2 |I_i|^2 Re Z_p``
on all four sheets (the driven one included). At ``z0 = Re Z_p`` both sides are
built from the *same* stored terminal readings, so the identity is algebraic and
1e-6 tests the assembly, not the discretisation. This example writes that
arithmetic out per drive, and sets beside it the **S-diagonal control**: delete
``S_kk`` from the column and ``P_acc,k`` grows by exactly ``1/2 |a_k|^2 |S_kk|^2``,
a miss of ``|S_kk|^2 / (1 - sum_i |S_ik|^2)`` against ``P_net,k`` -- the
closed-form ceiling, computed from S alone.

**Imported, never restated** (the `ANS-1` rule): ``STEP1B_POWER_WAVE_RTOL`` and
``STEP1B_CONTROL_MIN_MISS`` from ``tests/validation/test_port_drive_superposition.py``;
``build_four_port_sweep`` from ``tests/validation/test_port_birdcage_four_port.py``;
``_power_shares`` and ``POWER_BALANCE_BAND`` from
``tests/validation/test_birdcage_b1_plus_map.py``;
``run_n_port_sparameter_sweep`` / ``superpose_drives`` from ``fem_em_solver.ports``.
No ``test_*`` name and no pytest fixture is imported.

**Asserted per drive P1..P4:**

* (K) the package's ``accepted_power_w`` equals the closed form at rel 1e-12
  (arithmetic check, as the gate states it -- not a physics band);
* (N) ``P_net > 0``;
* (1b-i) ``|P_acc - P_net| / |P_net| <= STEP1B_POWER_WAVE_RTOL``;
* (c-Skk) the control: ceiling >= ``STEP1B_CONTROL_MIN_MISS``; control miss >=
  ``STEP1B_CONTROL_MIN_MISS``; control miss equals the ceiling at rel
  ``2 * STEP1B_POWER_WAVE_RTOL``;
* cross-check: the single-drive power residual
  ``|supplied - (phantom + conductor + sheets)| / supplied <= POWER_BALANCE_BAND``.

**Printed only, never gated:** ``|a_k|``, ``sum_i|S_ik|^2``, ``|S_kk|``,
supplied, sheets, ``P_vol,k = 1/2 int sigma|E_k|^2`` and
``|P_acc - P_vol|/P_acc``; the records from the gate log beside the readings
(no digit equality is asserted).

**Not shown:** no absolute S accuracy, no resonance or tuning claim, no AED
comparison, no 64 / 128 MHz; `POST-6` step 1b-ii's regression pin is not
asserted. Superposition twin: `ports:13`.

Negative result protocol: a (1b-i) or (c-Skk) miss through this path while the
gate module holds is an example/test divergence -- known-issues entry, stop,
never a re-record from this side.

Needs the complex DolfinX build; the runner sources it for the ``ports:`` group::

    ./run_examples.sh -e ports:22 -n 2 -t 600

Output: ``paraview_output/ports_22_power_wave_identity_per_drive_combined.xdmf``
-- mesh, ``CellTags``, the four sheet facet tags, and ``E_sq_P1`` .. ``E_sq_P4``
(``|E_k|^2`` as DG0, one distinctly named array per drive).
"""

from __future__ import annotations

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
from fem_em_solver.ports import run_n_port_sparameter_sweep, superpose_drives  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER  # noqa: E402
from tests.mesh.test_birdcage_port_tags import LEG_COUNT  # noqa: E402
from tests.validation.test_birdcage_b1_plus_map import (  # noqa: E402
    POWER_BALANCE_BAND,
    _power_shares,
)
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_drive_superposition import (  # noqa: E402
    STEP1B_CONTROL_MIN_MISS,
    STEP1B_POWER_WAVE_RTOL,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
BASENAME = "ports_22_power_wave_identity_per_drive"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
SETUP_BASENAME = "ports_22_power_wave_identity_per_drive_setup"

# The gate's own literal arithmetic check (package accepted_power_w vs the
# closed form); not a physics band.
CLOSED_FORM_RTOL = 1.0e-12

# Gate records (POST-6 step 1b, `20260923T095500Z_POST-6-step3-rerun-WF-7-step2`
# and `EX-49`), printed beside the readings -- never asserted for equality.
RECORD_IDENTITY_DEV = {"P1": 8.414e-16, "P2": 8.412e-16, "P3": 1.679e-15, "P4": 8.409e-16}
RECORD_CONTROL_MISS = {
    "P1": 7.659732e-01,
    "P2": 7.656687e-01,
    "P3": 7.634896e-01,
    "P4": 7.652279e-01,
}
RECORD_RESIDUAL = {
    "P1": 9.795751e-03,
    "P2": 9.796209e-03,
    "P3": 9.794985e-03,
    "P4": 9.795283e-03,
}


def _e_sq_dg0(e_complex, msh, name):
    """``|E|^2`` as a real DG0 field (interpolated polynomial expression)."""
    dg0 = fem.functionspace(msh, ("DG", 0))
    tmp = fem.Function(dg0)
    tmp.interpolate(
        fem.Expression(
            ufl.inner(e_complex, ufl.conj(e_complex)), dg0.element.interpolation_points
        )
    )
    out = fem.Function(dg0, name=name)
    out.x.array[:] = np.maximum(np.real(tmp.x.array), 0.0)
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
            "EX-150 -- the power-wave identity per drive, 4-leg birdcage, 10 MHz",
            flush=True,
        )
        print("=" * 78, flush=True)

    sweep = build_four_port_sweep()

    write_setup_figure(
        sweep["mesh"],
        sweep["cell_tags"],
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
        title="ports:22 -- the power-wave identity per drive (POST-6 step 1b gate)",
        comm=comm,
    )

    comm.Barrier()
    t0 = time.perf_counter()
    result = run_n_port_sparameter_sweep(
        sweep["problem"],
        sweep["port_defs"],
        lumped_sheet_ports=sweep["specs"],
        lumped_sheet_facet_tags=sweep["facet_tags"],
        keep_fields=True,
    )
    comm.Barrier()
    t_fields = time.perf_counter() - t0

    assert not result.is_placeholder, (
        "the field-keeping sweep fell back to the PORT-0 coupling heuristic; "
        "no solved field to read"
    )
    port_ids = list(result.port_ids)
    assert port_ids == ["P1", "P2", "P3", "P4"], f"port order is {port_ids}"

    if comm.rank == 0:
        print(
            f"\n[fixture] {sweep['cells']} cells (gate record 116085), "
            f"f = {sweep['problem'].frequency_hz:.3e} Hz, degree 1; "
            f"field-keeping sweep in {t_fields:.1f} s at -n {comm.size}",
            flush=True,
        )

    v_src = {spec.port_id: complex(spec.drive_voltage_v) for spec in sweep["specs"]}
    s_matrix = np.asarray(result.s_matrix, dtype=np.complex128)

    rows = {}
    for index, pid in enumerate(port_ids):
        responses = result.excitation_results[pid].responses
        solved = {
            "driven": pid,
            "fields": result.fields[pid],
            "source_voltage_v": v_src[pid],
            "currents": {p: complex(responses[p].current_a) for p in port_ids},
        }
        sh = _power_shares(sweep, solved)
        net = float(sh["supplied"]) - float(sh["sheet_total"])
        total = sh["phantom"] + sh["conductor"] + sh["sheet_total"]
        residual = abs(sh["supplied"] - total) / abs(sh["supplied"])
        vol = float(sh["phantom"] + sh["conductor"])

        unit = np.zeros(len(port_ids), dtype=np.complex128)
        unit[index] = 1.0
        drive_k = superpose_drives(result, unit)
        a_k = complex(drive_k.incident_waves[index])
        column = s_matrix[:, index]
        reflected = float(np.sum(np.abs(column) ** 2))
        closed_form = 0.5 * abs(a_k) ** 2 * (1.0 - reflected)
        control_column = column.copy()
        control_column[index] = 0.0
        control_acc = 0.5 * abs(a_k) ** 2 * (
            1.0 - float(np.sum(np.abs(control_column) ** 2))
        )
        acc = float(drive_k.accepted_power_w)
        rows[pid] = {
            "a": a_k,
            "acc": acc,
            "closed": closed_form,
            "net": net,
            "supplied": float(sh["supplied"]),
            "sheets": float(sh["sheet_total"]),
            "vol": vol,
            "reflected": reflected,
            "skk": abs(s_matrix[index, index]),
            "identity_dev": abs(acc - net) / abs(net),
            "vol_residual": abs(acc - vol) / abs(acc),
            "control_acc": control_acc,
            "control_miss": abs(control_acc - net) / abs(net),
            "ceiling": float(abs(s_matrix[index, index]) ** 2) / (1.0 - reflected),
            "residual": residual,
        }

    if comm.rank == 0:
        print(
            "\n[identity] P_acc,k = 1/2|a_k|^2(1 - sum_i|S_ik|^2) vs "
            "P_net,k = supplied_k - sum_i sheets_ik "
            f"(asserted rel <= {STEP1B_POWER_WAVE_RTOL:.0e}); records beside "
            "readings, no digit equality asserted",
            flush=True,
        )
        for pid in port_ids:
            q = rows[pid]
            print(
                f"  {pid}  P_acc {q['acc']:.9e} W   P_net {q['net']:.9e} W   "
                f"rel dev {q['identity_dev']:.3e} (record {RECORD_IDENTITY_DEV[pid]:.3e})\n"
                f"      |a_k| {abs(q['a']):.9e}, sum_i|S_ik|^2 {q['reflected']:.9f}, "
                f"|S_kk| {q['skk']:.9f}, supplied {q['supplied']:.9e} W, "
                f"sheets {q['sheets']:.9e} W\n"
                f"      P_vol,k = 1/2 int sigma|E_k|^2 = {q['vol']:.9e} W   "
                f"|P_acc-P_vol|/P_acc = {q['vol_residual']:.6e} (PRINTED, not asserted)\n"
                f"      S_kk deleted: P_acc,k -> {q['control_acc']:.9e} W, miss "
                f"{q['control_miss']:.6e} (record {RECORD_CONTROL_MISS[pid]:.6e}); "
                f"ceiling |S_kk|^2/(1-sum_i|S_ik|^2) = {q['ceiling']:.6e}; "
                f"floor {STEP1B_CONTROL_MIN_MISS:.0e} "
                f"({q['control_miss'] / STEP1B_CONTROL_MIN_MISS:.0f}x)\n"
                f"      single-drive power residual {q['residual']:.6e} "
                f"(record {RECORD_RESIDUAL[pid]:.6e}, band {POWER_BALANCE_BAND:.0e})",
                flush=True,
            )

    for pid in port_ids:
        q = rows[pid]
        assert abs(q["acc"] - q["closed"]) <= CLOSED_FORM_RTOL * abs(q["closed"]), (
            f"[{pid}] (K) package accepted_power_w {q['acc']:.9e} W vs closed "
            f"form {q['closed']:.9e} W: the example's arithmetic is not the package's"
        )
        assert q["net"] > 0.0, f"[{pid}] (N) supplied - sheets = {q['net']:.9e} W"
        assert q["identity_dev"] <= STEP1B_POWER_WAVE_RTOL, (
            f"[{pid}] (1b-i) the power-wave identity misses by "
            f"{q['identity_dev']:.6e} (band {STEP1B_POWER_WAVE_RTOL:.0e}); "
            "example/test divergence -- known-issues entry, stop"
        )
        assert q["ceiling"] >= STEP1B_CONTROL_MIN_MISS, (
            f"[{pid}] (c-Skk) ceiling {q['ceiling']:.6e} below the floor"
        )
        assert q["control_miss"] >= STEP1B_CONTROL_MIN_MISS, (
            f"[{pid}] (c-Skk) control miss {q['control_miss']:.6e} below the floor"
        )
        assert abs(q["control_miss"] - q["ceiling"]) <= (
            2.0 * STEP1B_POWER_WAVE_RTOL * abs(q["ceiling"])
        ), (
            f"[{pid}] (c-Skk) control miss {q['control_miss']:.9e} != ceiling "
            f"{q['ceiling']:.9e}"
        )
        assert q["residual"] <= POWER_BALANCE_BAND, (
            f"[{pid}] single-drive power residual {q['residual']:.6e} outside "
            f"{POWER_BALANCE_BAND:.0e}"
        )

    # ---- the picture: |E_k|^2 per drive, the four sheet facet tags ----------
    fields = {
        f"E_sq_{pid}": _e_sq_dg0(result.fields[pid].e_complex, sweep["mesh"], f"E_sq_{pid}")
        for pid in port_ids
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
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    if comm.rank == 0:
        print(
            f"\n[paraview] wrote {path}: mesh, CellTags, the four sheet facet "
            "tags, and E_sq_P1..E_sq_P4 (|E_k|^2, DG0).",
            flush=True,
        )
        print(
            f"\nAll gates hold. Total elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
