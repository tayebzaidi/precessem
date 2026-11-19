"""`OPS-81` (§9 item 828): the two-tag loss helper's scope guard, on synthetic cases.

``_loss_power_w`` (``test_port_drive_superposition.py``) sums ``½∫σ|E|²`` over
tags 3 (phantom) and 1 (conductor) only.  Its first statement,
``_assert_two_tag_scope``, must raise on a lossy cell outside (1, 3) — naming
the tag and the all-cells form — unless the case dict carries
``TWO_TAG_ACK_KEY``.  No mesh is built: the guard runs before any assembly, so
calling ``_loss_power_w(case, None, None)`` on a refused case proves the guard
is wired first (a passed-through case would fail on ``sweep["cell_tags"]``
with ``KeyError`` instead).

Positive controls (raise): 403 lossy; default material σ = 0.5 with no 403; a
material without ``sigma``.  Negative controls (do not raise): 403 at σ = 0;
403 lossy with the ack key; a sweep without ``problem``; the stock tags-1/3 map.
"""

from __future__ import annotations

import pytest

from fem_em_solver.core import HomogeneousMaterial
from fem_em_solver.core.time_harmonic import TimeHarmonicProblem
from tests.complex_mode import complex_only
from tests.validation.test_port_drive_superposition import (
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
    TWO_TAG_ACK_KEY,
    _assert_two_tag_scope,
    _loss_power_w,
)

pytestmark = complex_only

ROD_TAG = 403
VACUUM = HomogeneousMaterial(0.0, 1.0, 1.0)
SIGMA_WIRE = HomogeneousMaterial(5.8e7, 1.0, 1.0)
GEL = HomogeneousMaterial(0.5, 78.0, 1.0)


class _NoSigma:
    """A material-like object with no ``sigma`` attribute."""

    epsilon_r = 1.0
    mu_r = 1.0


def _case(extra=None, default=VACUUM):
    mmap = {CONDUCTOR_CELL_TAG: SIGMA_WIRE, PHANTOM_CELL_TAG: GEL}
    if extra is not None:
        mmap[ROD_TAG] = extra
    return {
        "problem": TimeHarmonicProblem(
            mesh=None, frequency_hz=64e6, material=default, material_map=mmap
        )
    }


def test_tags_are_the_stock_pair():
    assert {CONDUCTOR_CELL_TAG, PHANTOM_CELL_TAG} == {1, 3}


def test_positive_lossy_403_raises_through_loss_power_w():
    with pytest.raises(ValueError) as exc:
        _loss_power_w(_case(GEL), None, None)
    msg = str(exc.value)
    print(f"[OPS-81 W0] positive control message: {msg}", flush=True)
    assert "403" in msg
    assert "all-cells" in msg
    assert TWO_TAG_ACK_KEY in msg


def test_positive_lossy_default_material_raises():
    with pytest.raises(ValueError, match="default material"):
        _loss_power_w(_case(None, default=HomogeneousMaterial(0.5, 1.0, 1.0)), None, None)


def test_positive_material_without_sigma_raises():
    with pytest.raises(ValueError, match="403"):
        _loss_power_w(_case(_NoSigma()), None, None)


@pytest.mark.parametrize(
    "label",
    ["403_sigma_zero", "403_lossy_acknowledged", "no_problem", "stock_1_3"],
)
def test_negative_controls_pass_the_guard(label):
    if label == "403_sigma_zero":
        case = _case(VACUUM)
    elif label == "403_lossy_acknowledged":
        case = _case(SIGMA_WIRE)
        case[TWO_TAG_ACK_KEY] = "synthetic: (P) gated over all cells elsewhere"
    elif label == "no_problem":
        case = {}
    else:
        case = _case(None)
    assert _assert_two_tag_scope(case) is None
    # Through the helper: the guard passes, so the failure (if any) is the
    # assembly's missing inputs — a KeyError on the absent mesh/tags, never the
    # guard's ValueError.
    with pytest.raises((KeyError, AttributeError, TypeError)):
        _loss_power_w(case, None, None)
