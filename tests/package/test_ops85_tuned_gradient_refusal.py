"""`OPS-85` (§9 item 879): an exercising node for the tuned two-port gradient refusal.

U125 (`WF-1` step 28, §9 item 830) added ``initial_gradient_k_per_m`` to
``TransientSpec`` and refused it in ``_tuned_quadrature_temperature``
(``src/fem_em_solver/study.py``), but no node reached that site: its Outcome
(``docs/planning/chunks/WF-1.md:3590``), disclosure (3).  The same message
also stands in ``_study_temperature`` (the EM sphere) and
``_study_transient_temperature``, so (Ref) asserts the raising *frame*, not
just the message.

No mesh and no solve: ``study._tuned_result`` is replaced by a recorder that
raises ``_Reached``, and ``MeshGenerator.birdcage_port_domain`` by a counter
that raises ``_Tripwire`` (``tests/validation/test_geo38_f_clinical_3t_matrix.py``'s
form).  The tuned Study is ``test_wf16_study_tuned_f_human._study()``.  No
form is built, so no rank-0-only ``ArityMismatch`` can hang the peer.
"""

from __future__ import annotations

import pytest
from mpi4py import MPI

import fem_em_solver.study as study_mod
from fem_em_solver.io.mesh import MeshGenerator
from fem_em_solver.study import RobinCooling
from tests.package.test_wf16_study_tuned_f_human import _study
from tests.validation.test_mat12_transient_sar_heating_sphere import RHO_C
from tests.validation.test_mat15_coil_driven_heating_phantom import H_GEL, K_GEL, T_INF

MESSAGE = "initial_gradient_k_per_m is defined on a ThermalBoxGeometry only"
FRAME = "_tuned_quadrature_temperature"


class _Reached(Exception):
    pass


class _Tripwire(Exception):
    pass


def _say(msg):
    print(f"[OPS-85 rank {MPI.COMM_WORLD.rank}] {msg}", flush=True)


@pytest.fixture
def counters(monkeypatch):
    calls = {"tuned_result": 0, "builder": 0}

    def tuned_recorder(study, comm):
        calls["tuned_result"] += 1
        raise _Reached

    def builder(*a, **k):
        calls["builder"] += 1
        raise _Tripwire

    # The module object `_tuned_quadrature_temperature`'s globals are.
    monkeypatch.setattr(study_mod, "_tuned_result", tuned_recorder)
    monkeypatch.setattr(MeshGenerator, "birdcage_port_domain", staticmethod(builder))
    return calls


def _call(st, g):
    return st.temperature(
        RobinCooling(K_GEL, H_GEL, T_INF),
        transient=study_mod.TransientSpec(
            dt_s=1.0, n_steps=1, rho_c_j_per_m3_k=RHO_C, initial_gradient_k_per_m=g),
        balls={"mass_g": 10.0, "centres": [[0.0, 0.0, 0.0]]},
    )


@pytest.mark.parametrize("g", [10.0, 0.0], ids=["Ref", "Ctl-g0"])
def test_ops85_tuned_gradient_refused_at_its_frame(counters, g):
    """(Ref) / (Ctl-g0): exact message, raised from `_tuned_quadrature_temperature`."""
    st = _study()
    with pytest.raises(ValueError) as excinfo:
        _call(st, g)
    entry = excinfo.traceback[-1]
    _say(f"g={g!r}: {str(excinfo.value)!r} from {entry.name} line {entry.lineno + 1}; "
         f"_tuned_result calls {counters['tuned_result']}, builder calls {counters['builder']}")
    assert str(excinfo.value) == MESSAGE
    assert entry.name == FRAME
    assert counters["tuned_result"] == 0
    assert counters["builder"] == 0


def test_ops85_gradient_free_spec_reaches_tuned_result(counters):
    """(Ctl): the same call without a gradient passes every check and reaches `_tuned_result`."""
    st = _study()
    with pytest.raises(_Reached):
        _call(st, None)
    _say(f"g=None: reached _tuned_result; _tuned_result calls {counters['tuned_result']}, "
         f"builder calls {counters['builder']}")
    assert counters["tuned_result"] == 1
    assert counters["builder"] == 0
