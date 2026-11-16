#!/usr/bin/env python3
"""Run a STEP case: ``run_step_case.py STEP CONFIG --preflight --out DIR``.

`GEO-34` step 2: only the preflight exists — import the STEP through the
config's group map, tag it, write the combined XDMF / HDF5 (cell tags and the
wall + port-sheet facet tags) and the setup figure, and stop.  **No solve.**
Without ``--preflight`` the runner exits 2: its solve half is not written.

Collective under ``mpiexec``: every rank calls :func:`main`.
"""

from __future__ import annotations

import argparse
import sys

from mpi4py import MPI


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("step", help="STEP file")
    parser.add_argument("config", help="JSON group map (import_step + port_sheets)")
    parser.add_argument("--preflight", action="store_true",
                        help="mesh, tag, write XDMF + setup PNG; no solve")
    parser.add_argument("--out", required=True, help="output directory")
    args = parser.parse_args(argv)

    comm = MPI.COMM_WORLD
    if not args.preflight:
        if comm.rank == 0:
            print("run_step_case: only --preflight is implemented (no solve mode)",
                  file=sys.stderr, flush=True)
        return 2

    from fem_em_solver.io.step_preflight import write_step_preflight

    write_step_preflight(args.step, args.config, args.out, comm)
    return 0


if __name__ == "__main__":
    sys.exit(main())
