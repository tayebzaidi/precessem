"""`PORT-8` — reversing a lumped-sheet port's drive conjugates the 4x4 by D.

Replaces the two placeholder tests of chunk D6 (a monkeypatched solver and the
`PORT-0` string-label sign heuristic; the second of them was known-issues
entry 3).  On the lumped-sheet route, a port's orientation is
``LumpedSheetPortSpec.drive_direction``; ``PortDefinition.orientation`` is left
alone.

**Anchor (exact identity).**  Reversing the drive direction on P2 and P4 of the
gapped 4-leg birdcage gives ``S' = D S D`` and ``Z' = D Z D`` with
``D = diag(1, -1, 1, -1)``:

* the sheet term in the bilinear form does not depend on direction
  (``ports/lumped.py`` sheet assembly);
* the impressed source is proportional to the drive direction;
* the port current is read back along the drive direction;
* ``V = V_src - I Z_p``, and the power waves follow from V and I.

So each flipped port's V and I both change sign, which negates every entry
coupling a flipped port to an unflipped one (the adjacent class) and leaves the
self and opposite classes unchanged.

**Premise (asserted first).**  The cell count equals ``STEP1_CELL_RECORD``
(116 085) and the unflipped S reproduces ``S_10MHZ_EPS0_RECORD`` at
``RECORD_RTOL`` (both `-n 2` records, ``RECORD_RANK_WIDTH``).

**Band.**  ``‖S' − DSD‖_F/‖S‖_F ≤ 1e-12`` and the same on Z.  MUMPS drifts by
1 ULP across processes (known-issues), so the gate is a band, never bit
equality; the predicted residual is ≈ 1e-14 (the fixture's reciprocity floor is
2.152e-14).

**Negative control (asserted).**  The orientation-blind residual
``‖DSD − S‖_F/‖S‖_F`` on the in-run S is asserted ``≥ 0.5``: a gate that could
not tell S from DSD would pass a drive direction that is ignored.  Predicted
from the registered record: ``2‖S_adj‖/‖S‖ = 1.489`` (ceiling 2).

Scope: one fixture, 10 MHz, degree 1, the lumped route only.  No gap-voltage
route, no ring ports, no ``src/`` change.

Run (complex build, ``-n 2`` only)::

    scripts/testing/run_and_log.sh PORT-8-W1 "docker compose exec -T fem-em-solver \\
      bash -lc 'source /usr/local/bin/dolfinx-complex-mode && cd /workspace && \\
      PYTHONPATH=/workspace/src:/workspace FEM_EM_REQUIRE_COMPLEX=1 timeout -k 30 400 \\
      mpiexec -n 2 python3 -m pytest tests/environment \\
      tests/ports/test_port_orientation_sensitivity.py -v -s --tb=short'"
"""

from __future__ import annotations

import dataclasses
import time

import numpy as np
from mpi4py import MPI

from fem_em_solver.ports import run_n_port_sparameter_sweep

from tests.complex_mode import complex_only
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_circuit_layer_field import (
    RECORD_RANK_WIDTH,
    RECORD_RTOL,
    S_10MHZ_EPS0_RECORD,
    _residual,
)
from tests.validation.test_port_lumped_rlc_termination import STEP1_CELL_RECORD

FLIPPED_PORT_IDS = ("P2", "P4")
D = np.diag([1.0, -1.0, 1.0, -1.0]).astype(np.complex128)

ANCHOR_BAND = 1.0e-12
PREDICTED_ANCHOR_RESIDUAL = 1.0e-14  # the fixture's reciprocity floor is 2.152e-14

NEGATIVE_CONTROL_FLOOR = 0.5  # asserted
NEGATIVE_CONTROL_PREDICTED = 1.489  # 2‖S_adj‖/‖S‖ on S_10MHZ_EPS0_RECORD; printed

_ADJACENT = [(i, j) for i in range(4) for j in range(4) if (i - j) % 4 in (1, 3)]
_OPPOSITE = [(i, j) for i in range(4) for j in range(4) if (i - j) % 4 == 2]
_SELF = [(i, i) for i in range(4)]


def _rel(a, b, ref):
    return float(np.linalg.norm(a - b) / np.linalg.norm(ref))


@complex_only
def test_reversing_p2_p4_drive_conjugates_s_and_z_by_d():
    comm = MPI.COMM_WORLD
    assert comm.size == RECORD_RANK_WIDTH, (
        f"records are -n {RECORD_RANK_WIDTH} records; run at that width"
    )

    b = build_four_port_sweep(build_only=True)
    ncells = int(b["cells"])
    specs = list(b["specs"])
    ids = [s.port_id for s in specs]
    assert ids == ["P1", "P2", "P3", "P4"], ids

    flipped_specs = [
        dataclasses.replace(s, drive_direction=(0.0, 0.0, -1.0))
        if s.port_id in FLIPPED_PORT_IDS
        else s
        for s in specs
    ]
    for s in specs:
        assert tuple(s.drive_direction) == (0.0, 0.0, 1.0)

    def _sweep(port_specs):
        comm.Barrier()
        t0 = time.perf_counter()
        r = run_n_port_sparameter_sweep(
            b["problem"],
            b["port_defs"],
            lumped_sheet_ports=port_specs,
            lumped_sheet_facet_tags=b["facet_tags"],
            degree=b["degree"],
        )
        comm.Barrier()
        assert r.z_matrix is not None
        return (
            np.asarray(r.s_matrix, dtype=np.complex128),
            np.asarray(r.z_matrix, dtype=np.complex128),
            time.perf_counter() - t0,
        )

    s, z, t_s = _sweep(specs)

    # --- premise (asserted first) -----------------------------------------
    premise = _residual(s, np.asarray(S_10MHZ_EPS0_RECORD))
    if comm.rank == 0:
        print(
            f"\n[PORT-8] premise: {ncells} cells (record {STEP1_CELL_RECORD}), "
            f"mesh {b['mesh_time']:.2f} s; unflipped sweep {t_s:.2f} s at "
            f"-n {comm.size}; ‖S − S_10MHZ_EPS0_RECORD‖_F/‖rec‖_F = "
            f"{premise:.3e} (ASSERTED ≤ {RECORD_RTOL:.0e})",
            flush=True,
        )
    assert ncells == STEP1_CELL_RECORD
    assert premise <= RECORD_RTOL

    # --- the flipped sweep on the same mesh --------------------------------
    s_f, z_f, t_f = _sweep(flipped_specs)

    dsd = D @ s @ D
    dzd = D @ z @ D
    anchor_s = _rel(s_f, dsd, s)
    anchor_z = _rel(z_f, dzd, z)
    control = _rel(dsd, s, s)

    def _class_rows(m_new, m_old, entries):
        ratios = np.array([m_new[i, j] / m_old[i, j] for i, j in entries])
        return (
            float(np.max(np.abs(ratios - 1.0))),
            float(np.max(np.abs(ratios + 1.0))),
        )

    if comm.rank == 0:
        print(f"[PORT-8] flipped (P2, P4) sweep {t_f:.2f} s", flush=True)
        for name, m_new, m_old in (("S", s_f, s), ("Z", z_f, z)):
            for cls, entries in (
                ("self", _SELF),
                ("opposite", _OPPOSITE),
                ("adjacent", _ADJACENT),
            ):
                to_same, to_neg = _class_rows(m_new, m_old, entries)
                print(
                    f"    {name} {cls:8s}: max|X'/X − 1| = {to_same:.3e}   "
                    f"max|X'/X + 1| = {to_neg:.3e}",
                    flush=True,
                )
        print(
            f"[PORT-8] anchor ‖S' − DSD‖_F/‖S‖_F = {anchor_s:.3e}, "
            f"‖Z' − DZD‖_F/‖Z‖_F = {anchor_z:.3e} (ASSERTED ≤ {ANCHOR_BAND:.0e}; "
            f"predicted ≈ {PREDICTED_ANCHOR_RESIDUAL:.0e})\n"
            f"[PORT-8] negative control ‖DSD − S‖_F/‖S‖_F = {control:.4f} "
            f"(ASSERTED ≥ {NEGATIVE_CONTROL_FLOOR}; predicted "
            f"{NEGATIVE_CONTROL_PREDICTED} from the record, ceiling 2)",
            flush=True,
        )

    assert control >= NEGATIVE_CONTROL_FLOOR
    assert anchor_s <= ANCHOR_BAND
    assert anchor_z <= ANCHOR_BAND
