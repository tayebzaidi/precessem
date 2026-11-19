"""`WF-15` step 19 (§9 item 937, §10 U147): the tuned two-port ``duty`` refusal.

Before this step ``Study.temperature(..., duty=...)`` on a
``TwoPortQuadratureDrive`` Study dropped ``duty``: ``_study_temperature``
called the excitation's ``_temperature_port_study`` hook without it, and
neither the hook nor ``_tuned_quadrature_temperature`` takes one, so the run
was the continuous one.  The route now raises a named ``TypeError`` from
``_study_temperature`` before the hook is called.

`OPS-85`'s pattern (``tests/package/test_ops85_tuned_gradient_refusal.py``):
no mesh and no solve.  ``study._tuned_result`` is replaced by a recorder that
raises ``_Reached``, ``MeshGenerator.birdcage_port_domain`` by a counter that
raises ``_Tripwire``, and ``TwoPortQuadratureDrive._temperature_port_study`` is
wrapped in a counter that delegates to the original.  No form is built, so no
rank-0-only ``ArityMismatch`` can hang the peer.
"""

from __future__ import annotations

import inspect

import pytest
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.study import DutyCycle, RobinCooling, TransientSpec, TwoPortQuadratureDrive
from tests.package.test_wf16_study_tuned_f_human import _study
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, T_INF

MESSAGE = "duty is not defined on a tuned two-port Study (WF-15 step 19)"
FRAME = "_study_temperature"


class _Reached(Exception):
    pass


class _Tripwire(Exception):
    pass


def _say(msg):
    print(f"[WF-15-s19 rank {MPI.COMM_WORLD.rank}] {msg}", flush=True)


@pytest.fixture
def counters(monkeypatch):
    calls = {"hook": 0, "tuned_result": 0, "builder": 0}
    original_hook = TwoPortQuadratureDrive._temperature_port_study

    def hook_counter(self, *a, **k):
        calls["hook"] += 1
        return original_hook(self, *a, **k)

    def tuned_recorder(study, comm):
        calls["tuned_result"] += 1
        raise _Reached

    def builder(*a, **k):
        calls["builder"] += 1
        raise _Tripwire

    # The module object whose globals `_tuned_quadrature_temperature` reads.
    monkeypatch.setattr(study_mod, "_tuned_result", tuned_recorder)
    monkeypatch.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(builder))
    monkeypatch.setattr(TwoPortQuadratureDrive, "_temperature_port_study", hook_counter)
    return calls


_SENTINEL = object()


def _call(st, d=_SENTINEL):
    kw = {} if d is _SENTINEL else {"duty": d}
    return st.temperature(
        RobinCooling(K_GEL, H_GEL, T_INF),
        transient=TransientSpec(dt_s=1.0, n_steps=1, rho_c_j_per_m3_k=RHO_C),
        balls={"mass_g": 10.0, "centres": [[0.0, 0.0, 0.0]]},
        **kw,
    )


def _assert_refused(counters, d, arm):
    st = _study()
    with pytest.raises(TypeError) as excinfo:
        _call(st, d)
    entry = excinfo.traceback[-1]
    _say(
        f"{arm}: duty={d!r}: {str(excinfo.value)!r} from {entry.name} line {entry.lineno + 1}; "
        f"hook calls {counters['hook']}, _tuned_result calls {counters['tuned_result']}, "
        f"builder calls {counters['builder']}"
    )
    assert str(excinfo.value) == MESSAGE
    assert entry.name == FRAME
    assert counters["hook"] == 0
    assert counters["tuned_result"] == 0
    assert counters["builder"] == 0
    _say(f"VERDICT {arm} PASS")


def test_wf15_s19_tuned_duty_refused_at_its_frame(counters):
    """(Ref): ``DutyCycle(8.0, 0.5)`` refused from `_study_temperature`, before the hook."""
    _assert_refused(counters, DutyCycle(8.0, 0.5), "Ref")


@pytest.mark.parametrize("d", [DutyCycle(8.0, 1.0), 0.5], ids=["Ref-D1", "Ref-type"])
def test_wf15_s19_tuned_duty_refused_for_any_value(counters, d):
    """(Ref-D1) / (Ref-type): the keyword is refused, not the schedule, before any type check."""
    _assert_refused(counters, d, "Ref-D1" if isinstance(d, DutyCycle) else "Ref-type")


def test_wf15_s19_hook_signature_has_no_duty():
    """Record: the hook (unwrapped) takes no ``duty``."""
    params = list(inspect.signature(TwoPortQuadratureDrive._temperature_port_study).parameters)
    _say(f"hook signature parameters: {params}")
    assert "duty" not in params


def test_wf15_s19_duty_free_call_reaches_tuned_result(counters):
    """(Ctl): ``duty`` omitted, the same call reaches the hook and `_tuned_result`."""
    st = _study()
    with pytest.raises(_Reached):
        _call(st)
    _say(
        f"Ctl: duty omitted: reached _tuned_result; hook calls {counters['hook']}, "
        f"_tuned_result calls {counters['tuned_result']}, builder calls {counters['builder']}"
    )
    assert counters["hook"] == 1
    assert counters["tuned_result"] == 1
    assert counters["builder"] == 0
    _say("VERDICT Ctl PASS")
