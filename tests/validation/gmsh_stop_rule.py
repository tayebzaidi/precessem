"""A wall-clock stop rule around gmsh's ``generate`` / ``optimize`` calls.

`IMP-6` step 4 attempt 2 (§9 item 522, §10 I31). The 10-17 weekly priced the
disc-in-box build with a 450 s gmsh stop rule (`PROJECT_PLAN.md:4679–4683`):
a heavy window that crosses it is stopped and recorded as the build's price,
not killed at the container ceiling. Off unless ``FEM_EM_GMSH_STOP_RULE=1``;
never set in a re-gate or an ``xl`` command.

Stdlib only, and neither ``mpi4py`` nor ``gmsh`` is imported at module level,
so the (N5) controls can import this module in a bare subprocess.
"""

from __future__ import annotations

import faulthandler
import os
import sys
import threading
import time

GMSH_STOP_ENV = "FEM_EM_GMSH_STOP_RULE"
#: The 10-17 weekly's rule (`PROJECT_PLAN.md:4679–4683`).
GMSH_STOP_RULE_S = 450.0
GMSH_STOP_RC = 86
GMSH_STOP_MARKER = "[gmsh-stop-rule] STOPPED"


def gmsh_stop_armed() -> bool:
    """Unset or empty ⇒ False; ``"1"`` ⇒ True; anything else ⇒ ValueError."""
    raw = os.environ.get(GMSH_STOP_ENV, "")
    if raw == "":
        return False
    if raw == "1":
        return True
    raise ValueError(f"{GMSH_STOP_ENV}={raw!r}: only unset, empty or '1' are accepted")


class _StopWatch:
    """A daemon Timer that prints the marker and ``os._exit``s, plus a
    faulthandler fallback 30 s later that does not need the GIL."""

    def __init__(self, seconds: float, label: str):
        self.seconds = float(seconds)
        self.label = str(label)
        self._timer = None
        self.t0 = None

    def _fire(self):
        sys.stdout.write(
            f"\n{GMSH_STOP_MARKER} {self.label}: rule {self.seconds:g} s crossed "
            f"(price > {self.seconds:g} s)\n"
        )
        sys.stdout.flush()
        os._exit(GMSH_STOP_RC)

    def arm(self):
        self.t0 = time.perf_counter()
        self._timer = threading.Timer(self.seconds, self._fire)
        self._timer.daemon = True
        self._timer.start()
        faulthandler.dump_traceback_later(self.seconds + 30.0, exit=True)

    def cancel(self):
        if self._timer is not None:
            self._timer.cancel()
            self._timer = None
        faulthandler.cancel_dump_traceback_later()

    def elapsed(self) -> float:
        return 0.0 if self.t0 is None else time.perf_counter() - self.t0


def arm_gmsh_stop_rule(monkeypatch, comm, label, seconds: float = GMSH_STOP_RULE_S):
    """Wrap ``gmsh.model.mesh.generate`` / ``optimize`` on rank 0 only.

    ``generate`` arms the watch; ``optimize`` calls through, then cancels.
    Either wrapper cancels on an exception and re-raises, so a gmsh failure is
    never mislabelled as a stop.
    """
    if comm.rank != 0:
        return None
    import gmsh

    watch = _StopWatch(seconds, label)
    orig_generate = gmsh.model.mesh.generate
    orig_optimize = gmsh.model.mesh.optimize

    def _generate(*args, **kwargs):
        watch.arm()
        print(
            f"[gmsh-stop-rule] armed {seconds:g} s at generate" f"({', '.join(map(repr, args))})",
            flush=True,
        )
        try:
            return orig_generate(*args, **kwargs)
        except BaseException:
            watch.cancel()
            raise

    def _optimize(*args, **kwargs):
        try:
            out = orig_optimize(*args, **kwargs)
        except BaseException:
            watch.cancel()
            raise
        watch.cancel()
        print(
            f"[gmsh-stop-rule] under the rule: {watch.elapsed():.2f} s of {seconds:g}", flush=True
        )
        return out

    monkeypatch.setattr(gmsh.model.mesh, "generate", staticmethod(_generate))
    monkeypatch.setattr(gmsh.model.mesh, "optimize", staticmethod(_optimize))
    return watch
