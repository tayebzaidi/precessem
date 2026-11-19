"""`femem` command line (`WF-1` step 1, §9 item 256).

``femem run <file> [--output-dir DIR]`` solves a schema-v0 study file and
writes its outputs; ``femem check <file>`` parses only (on a sweep file it also
prints ``sweep_hz``, `WF-1` step 23).  ``femem run --ranks N <file>``
re-launches the run under ``mpiexec -n N`` (``FEMEM_MPIEXEC`` overrides the
launcher) through ``python -m fem_em_solver`` (`WF-13` step 3).
``femem figure {sparams,slice,probe} <source> <out>`` draws a record
through :mod:`fem_em_solver.post.figures` (`WF-8` step 1).
``femem run --estimate <file>`` prints ``Study.estimate()`` (one ``name: repr``
line per ``ResourceEstimate`` field, rank 0) and returns 0 without solving
(`WF-13` step 4); a file the estimator cannot price returns 2 with the API's
message unchanged.  A malformed file
prints ``femem: error: <message>`` on rank 0 and returns 2 on every rank.

``main`` returns an exit code and never calls ``sys.exit`` (the console-script
wrapper does), so it can be called in-process under ``mpiexec``.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from typing import Optional, Sequence

__all__ = ["main"]


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="femem", description="Run a fem-em-solver study file.")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="solve a YAML study file and write its outputs")
    run.add_argument("file")
    run.add_argument("--output-dir", default=None)
    run.add_argument(
        "--ranks",
        type=int,
        default=None,
        help="re-launch this run under mpiexec with N ranks (N >= 1)",
    )
    run.add_argument(
        "--estimate",
        action="store_true",
        help="print the resource estimate (cells, unknowns, predicted peak/time) and exit; no solve",
    )
    run.add_argument("--backend", choices=["dolfinx", "palace"], default="dolfinx")
    check = sub.add_parser("check", help="parse a YAML study file only (no mesh, no solve)")
    check.add_argument("file")
    fig = sub.add_parser(
        "figure", help="draw a figure from a written record (WF-8): sparams | slice | probe"
    )
    fig.add_argument("kind", choices=["sparams", "slice", "probe"])
    fig.add_argument("source", help="Touchstone file, XDMF bundle, or probe-series CSV")
    fig.add_argument("out", help="output image path")
    fig.add_argument("--name", default="real_SAR", help="slice: cell-centred field name")
    fig.add_argument("--axis", default="z", choices=["x", "y", "z"], help="slice: slab normal")
    fig.add_argument("--column", default="probe_c0_K", help="probe: CSV column")
    return parser


_LAUNCHER_ENV = ("PMI_SIZE", "PMI_RANK", "OMPI_COMM_WORLD_SIZE", "PMIX_RANK")


def _launch(args) -> int:
    """``femem run --ranks N`` (`WF-13` step 3, §9 item 294).

    Re-executes the same run as ``<launcher> -n N <python> -m fem_em_solver run
    <file> [--output-dir DIR]`` and returns the launcher's exit code.  Decided
    before ``mpi4py`` is imported, so no singleton PMI state reaches the child.
    The child stays in this process group (a group kill reaches it).
    """
    present = [k for k in _LAUNCHER_ENV if k in os.environ]
    if present:
        print(
            f"femem: error: --ranks inside an MPI launch ({', '.join(present)} set)",
            file=sys.stderr,
            flush=True,
        )
        return 2
    if args.ranks < 1:
        print(f"femem: error: --ranks must be >= 1, got {args.ranks}", file=sys.stderr, flush=True)
        return 2
    requested = os.environ.get("FEMEM_MPIEXEC") or "mpiexec"
    launcher = shutil.which(requested)
    if launcher is None:
        print(f"femem: error: MPI launcher not found: {requested!r}", file=sys.stderr, flush=True)
        return 2
    cmd = [launcher, "-n", str(args.ranks), sys.executable, "-m", "fem_em_solver", "run", args.file]
    if args.output_dir is not None:
        cmd += ["--output-dir", args.output_dir]
    print("femem: launching: " + " ".join(shlex.quote(c) for c in cmd), flush=True)
    sys.stderr.flush()
    return subprocess.run(cmd).returncode


def _estimate(path, comm) -> int:
    """``femem run --estimate`` (`WF-13` step 4, §9 item 325): prints
    ``Study.estimate(comm)`` on rank 0, no solve, no output directory.  A
    ``ValueError`` from the estimator (e.g. fewer than two priced rows at the
    Study's degree) is surfaced unchanged and returns 2 on every rank (review
    ruling 2026-10-08).  F-human Studies are priced since `WF-13` step 10 (FH1,
    degree 1) and step 11 (FH2, degree 2)."""
    import dataclasses

    from .study_yaml import load_study

    study = load_study(path)
    try:
        est = study.estimate(comm)
    except (ValueError, TypeError) as exc:
        if comm.rank == 0:
            print(f"femem: error: {exc}", file=sys.stderr, flush=True)
        return 2
    if comm.rank == 0:
        for f in dataclasses.fields(est):
            print(f"{f.name}: {getattr(est, f.name)!r}")
        sys.stdout.flush()
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)

    if (
        args.command == "run"
        and args.backend == "palace"
        and (args.ranks is not None or args.estimate)
    ):
        print(
            "femem: error: --backend palace cannot be combined with --ranks or " "--estimate",
            file=sys.stderr,
            flush=True,
        )
        return 2
    if args.command == "run" and args.ranks is not None:
        if args.estimate:
            # An estimate is priced at this process's own width; refuse before
            # any launch (`WF-13` step 4, §9 item 325).
            print(
                "femem: error: --ranks and --estimate cannot be combined "
                "(--estimate prices the current launch width)",
                file=sys.stderr,
                flush=True,
            )
            return 2
        return _launch(args)

    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    if args.command == "figure":
        from .post import figures

        if args.kind == "sparams":
            res = figures.plot_sparameters(args.source, args.out, comm=comm)
        elif args.kind == "slice":
            res = figures.plot_xdmf_slice(
                args.source, args.name, args.out, axis=args.axis, comm=comm
            )
        else:
            res = figures.plot_probe_series(args.source, args.column, args.out, comm=comm)
        if res is not None:
            print(f"femem: wrote {res.path}", flush=True)
        return 0

    from .study_yaml import StudyConfigError, _load, run_study_file

    try:
        if args.command == "check":
            # the parsed geometry kind, the family run_study_file records (§9 item 401)
            study, kind, outputs = _load(args.file)[:3]
            if kind == "thermal_box":  # `WF-1` step 15 (§9 item 521b): no EM material
                if comm.rank == 0:
                    print(f"family: {kind}")
                    print(f"excitation: {type(study.excitation).__name__}")
                    sys.stdout.flush()
                return 0
            if kind == "thermal_ball":  # `WF-1` step 14 (§9 item 500): no EM material
                if comm.rank == 0:
                    print(f"family: {kind}")
                    print(f"excitation: {type(study.excitation).__name__}")
                    sys.stdout.flush()
                return 0
            if comm.rank == 0:
                m = study.material
                print(f"family: {kind}")
                print(f"degree: {study.degree}")
                print(f"frequency_hz: {study.frequency_hz!r}")
                if "sweep_hz" in outputs:  # WF-1 step 23 (§9 item 589)
                    print(f"sweep_hz: {list(outputs['sweep_hz'])!r}")
                if "load_result" in outputs:  # WF-1 step 24 (§9 item 674)
                    print(f"load_result: {str(outputs['load_result'])!r}")
                print(f"material.sigma: {m.sigma!r}")
                print(f"material.epsilon_r: {m.epsilon_r!r}")
                print(f"excitation: {type(study.excitation).__name__}")
                sys.stdout.flush()
            return 0
        if args.estimate:
            return _estimate(args.file, comm)
        if args.backend == "palace":  # `ARCH-1` step 1a (§9 item 517)
            from .backends.palace import PalaceRunRequired
            from .study_yaml import _run_study_file_palace

            try:
                _run_study_file_palace(args.file, output_dir=args.output_dir, comm=comm)
            except PalaceRunRequired as exc:
                if comm.rank == 0:
                    print(
                        f"femem: palace run required: (cd {exc.cwd} && " f"{' '.join(exc.argv)})",
                        flush=True,
                    )
                return 3
            return 0
        run_study_file(args.file, output_dir=args.output_dir, comm=comm)
        return 0
    except StudyConfigError as exc:
        if comm.rank == 0:
            print(f"femem: error: {exc}", file=sys.stderr, flush=True)
        return 2
