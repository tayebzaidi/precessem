"""`WF-7` step 1 — the F-human 32x32 at 64 MHz, degree 1, gated on `PORT-13`'s
identities and `PORT-16`'s exact discrete power identity.

§9 item 21 (2026-09-21 review, chain step H3).  One window, heavy tier,
``-n 8``, complex build.  Self-consistency identities on one fixture at degree 1
— **no field, no B1+, no SAR, no absolute S**; degree 1 at human scale carries
the 2026-09-19 weekly's 5.50 % order caveat (`WF-7` step 0b: ``|S_2 - S_1|/|S_1|``
on ``S_driven(P17)``).

Fixture: built exactly as ``scripts/probes/wf7_step0_f_human_cost.py:229–293``
builds it — `GEO-25`'s F-human rung (``_params(F_HUMAN_RING_RADIUS,
scale_sizing=False)``, branch B) with ``ring_sheet_orientation="longitudinal"``,
all 32 ring sheets at the terminated ``z0``, `PORT-13`'s materials, PEC outer
boundary, 64 MHz.  All 32 ring ports are driven in turn through
`test_port_birdcage_ring_column._solve_one_drive`, the first factorising and the
other 31 back-substituting the held MUMPS factor (`PORT-19` step 3).  Per drive,
**before the fields are dropped**, the exact identity
``P_src,exact = P_vol + P_sheet,exact`` is evaluated with `PORT-16`'s own
helpers (`test_birdcage_power_identity`, imported) and the package's
`terminal_form_deficit` is printed beside it.

Anchors (asserted, every band imported, none restated):
  (i)   the cell record inside ``CELL_COUNT_BAND``;
  (ii)  32x32 reciprocity ``<= RECIPROCITY_BAND``;
  (iii) ``sigma_max <= MATRIX_PASSIVITY_CEILING`` (exactly 1.0, `PORT-13`);
  (iv)  exactly ``N_SYMMETRY_CLASSES`` (18) C16 x mirror classes from the
        **measured** sheet azimuths and ring membership, every pair within
        ``AZIMUTH_MATCH_DEG`` of a class, worst spread ``<= OPPOSITE_SPREAD_BAND``;
  (v)   the exact discrete power identity ``<= DISCRETE_IDENTITY_RTOL`` on the
        driven columns (all 32, unless the first identity evaluation costs more
        than ``IDENTITY_PRICE_STOP_S``, or the projected window exceeds
        ``WINDOW_BUDGET_S`` — then drives 1, 9, 17, 25 only, and the window
        says so).

Negative controls:
  * *asserted* (mesh level, no solve on the degenerate build): the same
    fixture with ``ring_sheet_orientation`` at its default has every ring
    sheet's extent along its drive direction ``< DEGENERATE_EXTENT_M``, while
    the longitudinal build's gap chord and every sheet's extent along its drive
    direction are ``> DEGENERATE_EXTENT_M``;
  * *asserted* (arithmetic): one column x ``CONTROL_COLUMN_SCALE`` raises the
    reciprocity ratio by more than ``CONTROL_RATIO_LIFT`` x the measured one;
  * *predicted, printed, never asserted* (rule (e)): the same control against
    ``MATRIX_CONTROL_MARGIN x RECIPROCITY_BAND`` — the review's arithmetic puts
    this fixture's ceiling at ~1.5–2.1e-3, at or past the 2e-3 bar.

`POST-7` step 1 (§9 item 58, chain H7) — two knobs, read inside the test:
``FEM_EM_WF7_DEGREE`` in {1, 2} (the solver's degree) and
``FEM_EM_WF7_IDENTITY_COLUMNS`` = N (drive P17.. the first N ports only, assert
(v) on each, bypass the stop rule; (ii)–(iv) and the matrix controls print
SKIPPED).  Unset = the landed run.  At degree 1, P17's ``P_src,exact`` and
pooled kappa reproduce the landed window within ``REUSE_REPRODUCTION_RTOL``
(asserted; printed only at degree != 1).  Negative control (rule (e)): (v)
with the undriven sheet of smallest ``P_sheet,exact`` omitted — printed on the
landed route (W1 2026-09-28: 479x the band), asserted >= 10x the band on the
identity-columns route (`xl-pending.md` entry 21's degree-2 command).

Printed, never gated: per-drive terminal residual beside the in-run pooled
``C/terminal - 1``; full-column ``sum_i |S_ij|^2``; ``S_driven(P17)`` beside the
step-0 record.
"""

from __future__ import annotations

import gzip
import os
import re
import resource
import time

import numpy as np
import pytest
from mpi4py import MPI

from fem_em_solver.core import (
    HomogeneousMaterial,
    TimeHarmonicProblem,
    TimeHarmonicSolver,
)
from fem_em_solver.io.mesh import MeshGenerator, _interface_facet_tags
from fem_em_solver.ports.lumped import LumpedSheetPortSpec
from fem_em_solver.ports.shares import terminal_form_deficit

from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_leg_offset import _sheet_azimuth_deg, _sheet_nodes
from tests.mesh.test_birdcage_port_sheet_prerequisite import CELL_COUNT_BAND
from tests.mesh.test_birdcage_port_sheets import PORT_LOWER, PORT_UPPER, SHEET_IFACE
from tests.mesh.test_birdcage_port_terminals import _interface_area_or_zero
from tests.mesh.test_birdcage_ring_gaps_scaleup import _ring_gap_frame
from tests.mesh.test_birdcage_ring_sheet_orientation import (
    DEGENERATE_EXTENT_M,
    _sheet_extent_along,
)
from tests.validation.test_birdcage_f_human_rung import (
    F_HUMAN_BRANCH_B_CELL_RECORD,
    F_HUMAN_RING_RADIUS,
    LEG_COUNT,
    _params,
    _ring_ports,
)
from tests.validation.test_birdcage_power_identity import (
    DISCRETE_IDENTITY_RTOL,
    _sheet_field_dissipation_w,
    _source_power_w,
)
from tests.validation.test_lossy_sphere_fullwave import SALINE_EPSILON_R, SALINE_SIGMA
from tests.validation.test_port_birdcage_four_port import TERMINATED_PORT_IMPEDANCE_OHM
from tests.validation.test_port_birdcage_lumped_column import (
    CONDUCTOR_CELL_TAG,
    PHANTOM_CELL_TAG,
)
from tests.validation.test_port_birdcage_ring_column import (
    AZIMUTH_MATCH_DEG,
    CONTROL_COLUMN_SCALE,
    OPPOSITE_SPREAD_BAND,
    _reciprocity_ratio,
    _solve_one_drive,
)
from tests.validation.test_port_birdcage_ring_matrix import (
    AZIMUTH_STEP_DEG,
    MATRIX_CONTROL_MARGIN,
    MATRIX_PASSIVITY_CEILING,
    N_AZIMUTH_CLASSES,
    N_SYMMETRY_CLASSES,
    REUSE_REPRODUCTION_RTOL,
    _azimuth_class,
)
from tests.validation.test_port_drive_superposition import _ring_quadrature_slots
from tests.validation.test_port_gap_voltage_impedance import SIGMA_WIRE_S_PER_M
from tests.validation.test_port_lumped_sheet_sweep import RECIPROCITY_BAND

FREQUENCY_HZ = 64.0e6
DEGREE = 1

# §9 item 21's cost rule: the identity's facet assembly per drive is
# unmeasured; if the first one exceeds this, the identity runs on drives
# 1, 9, 17, 25 only and the window prints that it did.
IDENTITY_PRICE_STOP_S = 15.0
IDENTITY_FALLBACK_DRIVES = (0, 8, 16, 24)
# The slot's wall budget for this one window (container ``timeout`` 590 s,
# the scheduled slot's cap on a single foreground command), and the two
# measured prices the projection uses: a held back-substitution 0.59 s and the
# F-human build 120.98 s, both at -n 8
# (`20260916T093548Z_WF-7-step0c.log:10419–10431`), rounded up.
WINDOW_BUDGET_S = 540.0
HELD_SOLVE_ESTIMATE_S = 1.0
CONTROL_BUILD_ESTIMATE_S = 135.0

# The asserted arithmetic control: a 1 % scale on one column must lift the
# reciprocity ratio by more than this factor over the measured ratio (§9 item
# 21; step 0c measured 1.084e-14 at -n 16, so the lift available is ~1e11).
CONTROL_RATIO_LIFT = 1.0e6

# The step-0 record, printed beside S_driven(P17) — never compared here
# (version- and width-tagged, `OPS-41`; `WF-7` step 0c anchor (a) owns it).
#
# `GEO-35` step 2a: leg-box-free fixture, `20261007T005016Z_GEO-35-s1-W2-wf7.log:9757`;
# boxed fixture, pre-`GEO-35`: "0.407423+0.344417j" (`WF-7` step 0c).
S_DRIVEN_STEP0_RECORD = "0.407315+0.344327j"
S_DRIVEN_STEP0_RECORD_BOXED = "0.407423+0.344417j"

# `POST-7` step 1 (§9 item 58, chain H7): the degree-1 reproduction anchors —
# P17's exact source power and in-run pooled kappa from the landed window
# (`20260921T093429Z_WF-7-step1.log:10484`), asserted at degree 1 within the
# imported `REUSE_REPRODUCTION_RTOL`; printed only at degree != 1.
#
# `GEO-35` step 2a (§9 item 351, (1*) licence per the 2026-10-09 18:00 review
# Q1): re-recorded on the leg-box-free fixture from step 1's knob-on reading
# `20261007T005016Z_GEO-35-s1-W2-wf7.log:9733` (2026-10-07).  The old values
# are the boxed fixture, pre-`GEO-35` (`20260921T093429Z_WF-7-step1.log:10484`),
# kept executable as `*_BOXED` for control (c2) and FEM_EM_WF7_LEG_PORTS=boxed.
P17_P_SRC_EXACT_RECORD_W = 7.036574789e-03
P17_POOLED_KAPPA_RECORD = 4.805416e-03
P17_P_SRC_EXACT_RECORD_W_BOXED = 7.037117061e-03
P17_POOLED_KAPPA_RECORD_BOXED = 4.590991e-03
# The negative control (rule (e)): the identity re-evaluated with the undriven
# sheet of smallest P_sheet,exact omitted must miss by at least this multiple
# of `DISCRETE_IDENTITY_RTOL`.  Predicted >= 460x (`PORT-16` (ii): exact >=
# terminal; 0.03597^2 * 2.5e-3 / 7.037e-3 ~ 4.6e-4, log :10532, :10484).
# Asserted on the FEM_EM_WF7_IDENTITY_COLUMNS route, printed on the landed one.
IDENTITY_NEGATIVE_CONTROL_FACTOR = 10.0

# `WF-7` step 0e, save-knob half (§9 item 213, chain H15a): the degree-1
# 32x32 reproduced **through the saved file** — sigma_max and the 18 C16 x
# mirror class means |S| from `POST-7` step 1 W1
# (`20260928T093759Z_POST-7-step1-W1.log`), asserted at degree 1 within the
# imported `REUSE_REPRODUCTION_RTOL`, printed only at degree != 1.
#
# `GEO-35` step 2a: re-recorded on the leg-box-free fixture from
# `20261007T005016Z_GEO-35-s1-W2-wf7.log:9736` (sigma_max) and `:9762–9779`
# (class means), 2026-10-07; the `_BOXED` dicts below are the boxed fixture,
# pre-`GEO-35`, from the POST-7 step 1 W1 lines cited there.
SIGMA_MAX_D1_RECORD = 0.999914557
CLASS_MEAN_D1_RECORD = {
    (True, 0): 5.330972338e-01,   # :9762
    (True, 1): 1.969175237e-01,   # :9763
    (True, 2): 1.018249854e-01,   # :9764
    (True, 3): 6.952878121e-02,   # :9765
    (True, 4): 6.068218153e-02,   # :9766
    (True, 5): 5.985549118e-02,   # :9767
    (True, 6): 6.098877066e-02,   # :9768
    (True, 7): 6.208144512e-02,   # :9769
    (True, 8): 6.249638442e-02,   # :9770
    (False, 0): 6.471881624e-01,  # :9771
    (False, 1): 6.024754705e-02,  # :9772
    (False, 2): 3.736485525e-02,  # :9773
    (False, 3): 5.568838659e-02,  # :9774
    (False, 4): 6.024391493e-02,  # :9775
    (False, 5): 6.171649476e-02,  # :9776
    (False, 6): 6.267984815e-02,  # :9777
    (False, 7): 6.339491565e-02,  # :9778
    (False, 8): 6.366749492e-02,  # :9779
}
SIGMA_MAX_D1_RECORD_BOXED = 0.999914648  # :10641 (= WF-7 step 1 :10519)
CLASS_MEAN_D1_RECORD_BOXED = {
    (True, 0): 5.332509183e-01,   # :10643
    (True, 1): 1.969479947e-01,   # :10644
    (True, 2): 1.018378835e-01,   # :10645
    (True, 3): 6.953095859e-02,   # :10646
    (True, 4): 6.068100879e-02,   # :10647
    (True, 5): 5.985427050e-02,   # :10648
    (True, 6): 6.098793026e-02,   # :10649
    (True, 7): 6.208078058e-02,   # :10650
    (True, 8): 6.249575359e-02,   # :10651
    (False, 0): 6.470569506e-01,  # :10652
    (False, 1): 6.030242376e-02,  # :10653
    (False, 2): 3.735191926e-02,  # :10654
    (False, 3): 5.568363668e-02,  # :10655
    (False, 4): 6.024289250e-02,  # :10656
    (False, 5): 6.171628899e-02,  # :10657
    (False, 6): 6.267954521e-02,  # :10658
    (False, 7): 6.339437979e-02,  # :10659
    (False, 8): 6.366685653e-02,  # :10660
}
# The asserted column-scale control on the reloaded array: P17 x
# CONTROL_COLUMN_SCALE must move the worst class-mean reproduction by at least
# this multiple of `REUSE_REPRODUCTION_RTOL` (ceiling ~0.01/32 ~ 3.1e-4, ~310x;
# predicted [2e-4, 4e-4]).
SAVE_CONTROL_FACTOR = 10.0
SAVE_ROOT = "/workspace/logs/"

# `POST-7` step 3 (§9 item 269, chain H7'): the per-column terminal accounting
# residual of the F-human 32x32 at 64 MHz registered as RECORDS at both element
# orders — never bands: nothing in this module compares a fresh residual
# against them.  (min, max) over the 32 `phase 3 drive N accounting` columns.
TERMINAL_RESIDUAL_RECORD_DEG1 = (
    4.061e-03,  # drive 15, 20260921T070008Z_WF-7-step0c.log:11187
    4.624e-03,  # drive 20, 20260921T070008Z_WF-7-step0c.log:11397
)
TERMINAL_RESIDUAL_RECORD_DEG2 = (
    1.081e-02,  # drive 15, 20260926T070009Z_WF-7-step0d.log:11187
    1.179e-02,  # drive 20, 20260926T070009Z_WF-7-step0d.log:11397
)
TERMINAL_RESIDUAL_P17_RECORD = {
    1: 4.187e-03,  # drive 1, 20260921T070008Z_WF-7-step0c.log:10599
    2: 1.104e-02,  # 20261004T070008Z_POST-7-step2.log:10955 (:10960); = step0d :10599
}
POST7_STEP0C_LOG = "20260921T070008Z_WF-7-step0c.log"
POST7_STEP0D_LOG = "20260926T070009Z_WF-7-step0d.log"
POST7_STEP2_LOG = "20261004T070008Z_POST-7-step2.log"


# `GEO-35` step 1 (§9 item 275): the F-human cell record this module's own
# landed degree-1 window measured (`20261005T003402Z_WF-7-H15a-r3-W2.log:10583`),
# printed beside the knob-on count only.
# `GEO-35` step 2a: the module's own default build is leg-box-free (502 704,
# `20261007T005016Z_GEO-35-s1-W2-wf7.log:9694`); the boxed fixture, pre-`GEO-35`,
# was 507 266 (`…GEO-35-s1-W1-mesh.log:19883`).  Both printed only.
WF7_CELL_RECORD = 502704
WF7_BOXED_CELL_RECORD = 507266

# The two saved degree-1 / degree-2 matrices on the boxed mesh (H15a, H15b):
# never overwritten (`GEO-35` step 2a, consumer rule (i)).
BOXED_NPZ_BASENAMES = ("wf7-s32-degree1-n8.npz", "wf7-s32-degree2-n16.npz")


def _boxed_knob() -> bool:
    """``FEM_EM_WF7_LEG_PORTS`` (`GEO-35` step 2a), read at call time.

    Unset or empty ⇒ ``False``: the default leg-box-free build
    (``emit_leg_ports=False``), every reproduction asserted against the
    re-recorded constants.  ``boxed`` ⇒ ``True``: no key (the pre-switch
    build); the reproductions are printed against the ``*_BOXED`` constants,
    while reciprocity, passivity, C16 x mirror and the exact identity stay
    asserted.  Any other value raises.
    """
    raw = os.environ.get("FEM_EM_WF7_LEG_PORTS", "").strip()
    if raw == "":
        return False
    if raw != "boxed":
        raise ValueError(f"FEM_EM_WF7_LEG_PORTS must be 'boxed' or unset, got {raw!r}")
    if os.environ.get("FEM_EM_WF7_SAVE_NPZ", "").strip():
        raise ValueError(
            "FEM_EM_WF7_LEG_PORTS=boxed and FEM_EM_WF7_SAVE_NPZ are exclusive: the "
            "boxed matrices are already saved and never overwritten"
        )
    return True


def _save_knob() -> str | None:
    """``FEM_EM_WF7_SAVE_NPZ``, read at call time (`xl-pending.md` entry 29).

    Unset or empty ⇒ ``None`` (no file).  Otherwise an absolute path under
    ``/workspace/logs/`` ending ``.npz``; anything else, or the knob set
    together with ``FEM_EM_WF7_IDENTITY_COLUMNS``, raises ``ValueError`` (the
    file holds only the full 32x32).
    """
    raw = os.environ.get("FEM_EM_WF7_SAVE_NPZ", "").strip()
    if not raw:
        return None
    if not os.path.isabs(raw):
        raise ValueError(f"FEM_EM_WF7_SAVE_NPZ must be absolute, got {raw!r}")
    norm = os.path.normpath(raw)
    if not norm.startswith(SAVE_ROOT) or norm != raw:
        raise ValueError(f"FEM_EM_WF7_SAVE_NPZ must lie under {SAVE_ROOT}, got {raw!r}")
    if not norm.endswith(".npz"):
        raise ValueError(f"FEM_EM_WF7_SAVE_NPZ must end .npz, got {raw!r}")
    if os.path.basename(norm) in BOXED_NPZ_BASENAMES:
        raise ValueError(
            f"FEM_EM_WF7_SAVE_NPZ refuses {os.path.basename(norm)!r}: a boxed-mesh "
            "record file (GEO-35 step 2a), never overwritten"
        )
    if os.environ.get("FEM_EM_WF7_IDENTITY_COLUMNS", "").strip():
        raise ValueError(
            "FEM_EM_WF7_SAVE_NPZ and FEM_EM_WF7_IDENTITY_COLUMNS are exclusive: the "
            "file holds only the full 32x32"
        )
    return norm


def _class_means(s, az, zc):
    """``(classes, spreads, worst_residual_deg)`` — the (iv) class table, lifted
    from the inline block; ``spreads[key] = (spread, n, mean |S|)``."""
    n = len(az)
    classes: dict = {}
    worst_residual_deg = 0.0
    for i in range(n):
        for j in range(n):
            steps, residual_deg = _azimuth_class(az[i], az[j])
            worst_residual_deg = max(worst_residual_deg, residual_deg)
            classes.setdefault((bool(zc[i] * zc[j] > 0.0), steps), []).append(
                float(abs(s[i, j]))
            )
    spreads = {}
    for key, mags in classes.items():
        m = np.array(mags)
        spreads[key] = (float((m.max() - m.min()) / m.mean()), len(m), float(m.mean()))
    return classes, spreads, worst_residual_deg


def _class_mean_worst_rel(spreads) -> float:
    return max(abs(spreads[k][2] / CLASS_MEAN_D1_RECORD[k] - 1.0)
               for k in CLASS_MEAN_D1_RECORD)


ORIENTATION_MODES = ("inline", "split", "only")


def _orientation_knob() -> str:
    """``FEM_EM_WF7_ORIENTATION_CONTROL``, read at call time (§9 item 227).

    Unset or empty ⇒ ``"inline"`` (the landed path: the default-orientation
    build runs inside the matrix test, after the matrix).  ``split`` ⇒ the
    matrix test skips that build and its one assert; ``only`` ⇒ the matrix test
    skips before any mesh build and ``test_wf7_default_orientation_control``
    runs.  Anything else, or ``only`` with ``FEM_EM_WF7_SAVE_NPZ`` set (no
    matrix to save), raises ``ValueError``.
    """
    raw = os.environ.get("FEM_EM_WF7_ORIENTATION_CONTROL", "").strip()
    mode = raw if raw else "inline"
    if mode not in ("split", "only") and raw:
        raise ValueError(
            f"FEM_EM_WF7_ORIENTATION_CONTROL must be unset, 'split' or 'only', got {raw!r}"
        )
    if mode == "only" and os.environ.get("FEM_EM_WF7_SAVE_NPZ", "").strip():
        raise ValueError(
            "FEM_EM_WF7_ORIENTATION_CONTROL=only and FEM_EM_WF7_SAVE_NPZ are exclusive: "
            "there is no matrix to save"
        )
    return mode


def _default_orientation_extent(comm, ring_ports):
    """The mesh-level orientation control's build (no solve), lifted verbatim
    from the inline path: ``(default_extent, default_cells, build_s)``."""
    t0 = time.perf_counter()
    kwargs_default = _params(F_HUMAN_RING_RADIUS, scale_sizing=False)
    assert "ring_sheet_orientation" not in kwargs_default
    msh_d, cell_tags_d, _f_d, _diag_d = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kwargs_default
    )
    tags_d = _interface_facet_tags(
        msh_d, cell_tags_d,
        {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports},
    )
    default_extent = {
        i: _sheet_extent_along(msh_d, tags_d, SHEET_IFACE + i,
                               np.asarray(_ring_gap_frame(i, LEG_COUNT)[0], dtype=float),
                               comm)
        for i in ring_ports
    }
    default_cells = msh_d.topology.index_map(3).size_global
    return default_extent, default_cells, time.perf_counter() - t0


def test_wf7_orientation_knob_contract(monkeypatch):
    """W0: ``FEM_EM_WF7_ORIENTATION_CONTROL``'s value convention; no mesh."""
    monkeypatch.delenv("FEM_EM_WF7_ORIENTATION_CONTROL", raising=False)
    monkeypatch.delenv("FEM_EM_WF7_SAVE_NPZ", raising=False)
    monkeypatch.delenv("FEM_EM_WF7_IDENTITY_COLUMNS", raising=False)
    assert _orientation_knob() == "inline"
    monkeypatch.setenv("FEM_EM_WF7_ORIENTATION_CONTROL", "")
    assert _orientation_knob() == "inline"
    for ok in ("split", "only"):
        monkeypatch.setenv("FEM_EM_WF7_ORIENTATION_CONTROL", ok)
        assert _orientation_knob() == ok
    for bad in ("yes", "0", "inline"):
        monkeypatch.setenv("FEM_EM_WF7_ORIENTATION_CONTROL", bad)
        try:
            _orientation_knob()
        except ValueError:
            pass
        else:
            raise AssertionError(f"_orientation_knob accepted {bad!r}")
    monkeypatch.setenv("FEM_EM_WF7_ORIENTATION_CONTROL", "split")
    monkeypatch.setenv("FEM_EM_WF7_SAVE_NPZ", "/workspace/logs/wf7/wf7-s32-degree1-n8-noleg.npz")
    assert _orientation_knob() == "split"
    monkeypatch.setenv("FEM_EM_WF7_ORIENTATION_CONTROL", "only")
    try:
        _orientation_knob()
    except ValueError:
        pass
    else:
        raise AssertionError("_orientation_knob accepted only + FEM_EM_WF7_SAVE_NPZ")
    print("[WF-7 H15a] W0 orientation-knob contract: unset/empty -> inline, split/only "
          "accepted, 'yes'/'0'/'inline' and only + save path raise ValueError", flush=True)


@complex_only
def test_wf7_default_orientation_control():
    """The mesh-level orientation control in its own window
    (``FEM_EM_WF7_ORIENTATION_CONTROL=only``; skipped otherwise)."""
    if _orientation_knob() != "only":
        pytest.skip("FEM_EM_WF7_ORIENTATION_CONTROL != only: the control runs inline "
                    "or is split out of this window")
    comm = MPI.COMM_WORLD
    ring_ports = _ring_ports()
    assert len(ring_ports) == 2 * LEG_COUNT
    default_extent, default_cells, build_s = _default_orientation_extent(comm, ring_ports)
    rss = _rss_summed_gib(comm)  # collective: every rank, never under the rank-0 branch
    if comm.rank == 0:
        print(f"[WF-7 step1] orientation control: default build {default_cells} cells "
              f"({build_s:.2f} s); max extent along phi_hat "
              f"{max(default_extent.values()):.3e} m (< {DEGENERATE_EXTENT_M:.0e}, "
              f"ASSERTED); build_s {build_s:.2f}  -n {comm.size}, summed ru_maxrss "
              f"{rss:.3f} GiB", flush=True)
    assert max(default_extent.values()) < DEGENERATE_EXTENT_M, (
        "the default-orientation ring sheets have extent along phi_hat — the "
        "mis-paired orientation control does not hold"
    )


def test_wf7_save_knob_contract(monkeypatch):
    """W0: the knob's value convention; no mesh, no file written."""
    monkeypatch.delenv("FEM_EM_WF7_SAVE_NPZ", raising=False)
    monkeypatch.delenv("FEM_EM_WF7_IDENTITY_COLUMNS", raising=False)
    assert _save_knob() is None
    monkeypatch.setenv("FEM_EM_WF7_SAVE_NPZ", "")
    assert _save_knob() is None
    # `GEO-35` step 2a: the good path moved to the leg-box-free file (a contract
    # literal, not a band); the two boxed record files are refused.
    good = "/workspace/logs/wf7/wf7-s32-degree1-n8-noleg.npz"
    monkeypatch.setenv("FEM_EM_WF7_SAVE_NPZ", good)
    assert _save_knob() == good
    for bad in ("logs/wf7/x.npz", "/workspace/docs/testing/logs/x.npz",
                "/workspace/logs/../docs/x.npz", "/workspace/logs/wf7/x.npy",
                "/tmp/x.npz", "/workspace/logs/wf7/wf7-s32-degree1-n8.npz",
                "/workspace/logs/wf7/wf7-s32-degree2-n16.npz"):
        monkeypatch.setenv("FEM_EM_WF7_SAVE_NPZ", bad)
        try:
            _save_knob()
        except ValueError:
            pass
        else:
            raise AssertionError(f"_save_knob accepted {bad!r}")
    monkeypatch.setenv("FEM_EM_WF7_SAVE_NPZ", good)
    monkeypatch.setenv("FEM_EM_WF7_IDENTITY_COLUMNS", "1")
    try:
        _save_knob()
    except ValueError:
        pass
    else:
        raise AssertionError("_save_knob accepted the identity-columns route")
    print(f"[WF-7 H15a] W0 save-knob contract: unset/empty -> None, {good} accepted, "
          "7 bad paths (incl. the 2 boxed record files) + identity-columns raise "
          "ValueError; no file written", flush=True)


def test_wf7_leg_port_knob_contract(monkeypatch):
    """W0 (`GEO-35` step 2a): ``FEM_EM_WF7_LEG_PORTS``'s convention; no mesh."""
    monkeypatch.delenv("FEM_EM_WF7_LEG_PORTS", raising=False)
    monkeypatch.delenv("FEM_EM_WF7_SAVE_NPZ", raising=False)
    assert _boxed_knob() is False
    monkeypatch.setenv("FEM_EM_WF7_LEG_PORTS", "")
    assert _boxed_knob() is False
    monkeypatch.setenv("FEM_EM_WF7_LEG_PORTS", "boxed")
    assert _boxed_knob() is True
    for bad in ("1", "noleg", "Boxed", "true"):
        monkeypatch.setenv("FEM_EM_WF7_LEG_PORTS", bad)
        with pytest.raises(ValueError):
            _boxed_knob()
    monkeypatch.setenv("FEM_EM_WF7_LEG_PORTS", "boxed")
    monkeypatch.setenv("FEM_EM_WF7_SAVE_NPZ", "/workspace/logs/wf7/wf7-s32-degree1-n8-noleg.npz")
    with pytest.raises(ValueError):
        _boxed_knob()
    # The retired step-1 knob is gone (no silent half-honoured setting).
    assert "_geo35_knob" not in globals()
    # (c2)'s comparator inputs exist and differ from the re-records.
    assert P17_P_SRC_EXACT_RECORD_W != P17_P_SRC_EXACT_RECORD_W_BOXED
    assert set(CLASS_MEAN_D1_RECORD) == set(CLASS_MEAN_D1_RECORD_BOXED)
    print("[GEO-35 s2a] W0 leg-port knob contract: unset/empty -> dropped, 'boxed' "
          "-> boxed, 4 bad values + boxed with a save path raise ValueError", flush=True)


_POST7_ACCOUNTING_RE = re.compile(
    r"phase 3 drive (\d+) accounting \(printed.*?\bresidual (\d\.\d{3}e[-+]\d{2})\b")
_POST7_P17_RE = re.compile(
    r"\[WF-7 step1\] drive\s+1/1 P17:.*?terminal residual (\d\.\d{3}e[-+]\d{2})\b")


def _post7_open_log(name: str):
    """Open a tracked harness log, gzipped (``.log.gz``) or plain (``.log``)."""
    logs = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__)))), "docs", "testing", "logs")
    gz = os.path.join(logs, name + ".gz")
    if os.path.exists(gz):
        return gzip.open(gz, "rt")
    return open(os.path.join(logs, name), "rt")


def _post7_parse_columns(name: str) -> dict[int, str]:
    """{drive: residual token} from a step-0 log; tokens kept as printed strings."""
    out: dict[int, str] = {}
    with _post7_open_log(name) as fh:
        for line in fh:
            m = _POST7_ACCOUNTING_RE.search(line)
            if m:
                drive = int(m.group(1))
                assert drive not in out, f"{name}: drive {drive} printed twice"
                out[drive] = m.group(2)
    return out


def _post7_parse_p17(name: str) -> list[str]:
    with _post7_open_log(name) as fh:
        return [m.group(1) for line in fh if (m := _POST7_P17_RE.search(line))]


def _post7_minmax(cols: dict[int, str]) -> tuple[str, str]:
    lo = min(cols.values(), key=float)
    hi = max(cols.values(), key=float)
    return lo, hi


def test_post7_terminal_residual_records():
    """`POST-7` step 3: the terminal-residual records equal their log prints.

    Anchor at print precision (string equality of ``f"{v:.3e}"``), the P17
    in-range identities, and the asserted swapped-order negative control.
    """
    comm = MPI.COMM_WORLD
    parsed = None
    if comm.rank == 0:
        parsed = (
            _post7_parse_columns(POST7_STEP0C_LOG),
            _post7_parse_columns(POST7_STEP0D_LOG),
            _post7_parse_p17(POST7_STEP2_LOG),
        )
    d1, d2, p17_step2 = comm.bcast(parsed, root=0)

    # Count: exactly 32 columns per step-0 log, drives 1..32.
    assert sorted(d1) == list(range(1, 33)), f"step0c columns {sorted(d1)}"
    assert sorted(d2) == list(range(1, 33)), f"step0d columns {sorted(d2)}"
    assert len(p17_step2) == 1, f"step2 P17 tokens {p17_step2}"

    d1_lo, d1_hi = _post7_minmax(d1)
    d2_lo, d2_hi = _post7_minmax(d2)
    fmt = "{:.3e}".format
    print(f"[POST-7 step3] degree 1 (step0c): 32 columns, min {d1_lo} max {d1_hi} "
          f"P17 {d1[1]} | record {fmt(TERMINAL_RESIDUAL_RECORD_DEG1[0])}.."
          f"{fmt(TERMINAL_RESIDUAL_RECORD_DEG1[1])} P17 {fmt(TERMINAL_RESIDUAL_P17_RECORD[1])}",
          flush=True)
    print(f"[POST-7 step3] degree 2 (step0d): 32 columns, min {d2_lo} max {d2_hi} "
          f"P17 {d2[1]}; step2 P17 {p17_step2[0]} | record "
          f"{fmt(TERMINAL_RESIDUAL_RECORD_DEG2[0])}..{fmt(TERMINAL_RESIDUAL_RECORD_DEG2[1])} "
          f"P17 {fmt(TERMINAL_RESIDUAL_P17_RECORD[2])}", flush=True)

    # Anchor: string equality at print precision.
    assert (fmt(TERMINAL_RESIDUAL_RECORD_DEG1[0]), fmt(TERMINAL_RESIDUAL_RECORD_DEG1[1])) \
        == (d1_lo, d1_hi)
    assert (fmt(TERMINAL_RESIDUAL_RECORD_DEG2[0]), fmt(TERMINAL_RESIDUAL_RECORD_DEG2[1])) \
        == (d2_lo, d2_hi)
    assert fmt(TERMINAL_RESIDUAL_P17_RECORD[1]) == d1[1]
    assert fmt(TERMINAL_RESIDUAL_P17_RECORD[2]) == p17_step2[0]

    # Record identities (exact): P17 inside its order's range; step 2 reprints step 0d drive 1.
    for order, (lo, hi) in ((1, TERMINAL_RESIDUAL_RECORD_DEG1),
                            (2, TERMINAL_RESIDUAL_RECORD_DEG2)):
        assert lo <= TERMINAL_RESIDUAL_P17_RECORD[order] <= hi, order
    assert p17_step2[0] == d2[1], (p17_step2[0], d2[1])

    # Negative control (asserted): degree-1 constants against the degree-2 log miss.
    miss_lo = fmt(TERMINAL_RESIDUAL_RECORD_DEG1[0]) != d2_lo
    miss_hi = fmt(TERMINAL_RESIDUAL_RECORD_DEG1[1]) != d2_hi
    disjoint = TERMINAL_RESIDUAL_RECORD_DEG1[1] < TERMINAL_RESIDUAL_RECORD_DEG2[0]
    ratio = float(d2_lo) / TERMINAL_RESIDUAL_RECORD_DEG1[0]
    print(f"[POST-7 step3] NEG swapped (deg-1 constants vs step0d): min "
          f"{fmt(TERMINAL_RESIDUAL_RECORD_DEG1[0])} vs {d2_lo} unequal={miss_lo} "
          f"({ratio:.2f}x), max {fmt(TERMINAL_RESIDUAL_RECORD_DEG1[1])} vs {d2_hi} "
          f"unequal={miss_hi}; ranges disjoint={disjoint} (ASSERTED)", flush=True)
    assert miss_lo and miss_hi and disjoint


def _post7_knobs() -> tuple[int, int | None]:
    """``(degree, identity_columns)`` from the environment, read at call time.

    ``FEM_EM_WF7_DEGREE`` in {1, 2} (unset = 1, the landed run);
    ``FEM_EM_WF7_IDENTITY_COLUMNS`` = N drives ports 1..N only, asserts the
    identity on every one and bypasses the stop rule; the matrix gates
    (ii)-(iv) and their controls print SKIPPED (unset = all 32 drives).  The
    names are `xl-pending.md` entry 21's contract.
    """
    raw_deg = os.environ.get("FEM_EM_WF7_DEGREE", "").strip()
    degree = int(raw_deg) if raw_deg else DEGREE
    if degree not in (1, 2):
        raise ValueError(f"FEM_EM_WF7_DEGREE must be 1 or 2, got {raw_deg!r}")
    raw_cols = os.environ.get("FEM_EM_WF7_IDENTITY_COLUMNS", "").strip()
    cols = int(raw_cols) if raw_cols else None
    if cols is not None and not 1 <= cols <= 2 * LEG_COUNT:
        raise ValueError(
            f"FEM_EM_WF7_IDENTITY_COLUMNS must be in 1..{2 * LEG_COUNT}, got {raw_cols!r}"
        )
    return degree, cols


def _sheet_centre_z(msh, facet_tags, tag, comm) -> float:
    """``z`` of one sheet's bounding-box centre [m], reduced over ranks."""
    pts = _sheet_nodes(msh, facet_tags, tag)
    lo = float(pts[:, 2].min()) if pts.shape[0] else np.inf
    hi = float(pts[:, 2].max()) if pts.shape[0] else -np.inf
    lo = comm.allreduce(lo, op=MPI.MIN)
    hi = comm.allreduce(hi, op=MPI.MAX)
    return 0.5 * (lo + hi)


def _rss_summed_gib(comm) -> float:
    local = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024.0 / 2**30
    return float(comm.allreduce(local, op=MPI.SUM))


@complex_only
def test_the_f_human_ring_matrix_carries_port13_and_port16_identities():
    comm = MPI.COMM_WORLD
    t_window = time.perf_counter()
    degree, identity_columns = _post7_knobs()
    save_path = _save_knob()
    orientation = _orientation_knob()
    if orientation == "only":
        pytest.skip("FEM_EM_WF7_ORIENTATION_CONTROL=only: the matrix runs in its own "
                    "window; test_wf7_default_orientation_control runs here")

    def say(msg: str) -> None:
        if comm.rank == 0:
            print(f"[WF-7 step1] {msg}", flush=True)

    say(f"start: -n {comm.size}, f = {FREQUENCY_HZ:.3e} Hz, degree {degree}, "
        f"ring_radius {F_HUMAN_RING_RADIUS} m, branch B, longitudinal ring sheets; "
        f"FEM_EM_WF7_DEGREE={degree} FEM_EM_WF7_IDENTITY_COLUMNS={identity_columns} "
        f"FEM_EM_WF7_SAVE_NPZ={save_path} FEM_EM_WF7_ORIENTATION_CONTROL={orientation}")

    # ---- the fixture, as the step-0 probe builds it (:229–293) -----------
    kwargs = _params(F_HUMAN_RING_RADIUS, scale_sizing=False)
    kwargs["ring_sheet_orientation"] = "longitudinal"
    # `GEO-35` step 2a: the default build drops the leg port boxes; the
    # FEM_EM_WF7_LEG_PORTS=boxed knob rebuilds the pre-switch fixture (no key)
    # with the reproductions printed against the `*_BOXED` constants.
    geo35 = _boxed_knob()
    if not geo35:
        kwargs["emit_leg_ports"] = False
    t0 = time.perf_counter()
    msh, cell_tags, _facets, diag = MeshGenerator.birdcage_port_domain(
        comm=comm, return_diagnostics=True, **kwargs
    )
    build_s = time.perf_counter() - t0
    n_cells = msh.topology.index_map(3).size_global
    rel = n_cells / F_HUMAN_BRANCH_B_CELL_RECORD - 1.0
    say(f"mesh: cells {n_cells} vs record {F_HUMAN_BRANCH_B_CELL_RECORD}, relative "
        f"{rel:.3e} (band {CELL_COUNT_BAND}, imported)  build {build_s:.2f} s")
    say(f"GEO-35 s2a: {'BOXED knob (pre-switch build)' if geo35 else 'leg-box-free default'}"
        f": cells {n_cells} vs WF7_CELL_RECORD {WF7_CELL_RECORD} rel "
        f"{n_cells / WF7_CELL_RECORD - 1.0:+.4e}, vs WF7_BOXED_CELL_RECORD "
        f"{WF7_BOXED_CELL_RECORD} rel {n_cells / WF7_BOXED_CELL_RECORD - 1.0:+.4e} (printed)")
    from tests.mesh.test_geo35_leg_port_boxes import _owned_census

    leg_tags = ([PORT_LOWER + i for i in range(1, LEG_COUNT + 1)]
                + [PORT_UPPER + i for i in range(1, LEG_COUNT + 1)])
    leg_census_total = sum(_owned_census(msh, cell_tags, leg_tags, comm).values())
    say(f"(N') leg port-box census (allreduced, owned cells, tags PORT_LOWER/UPPER + "
        f"1..{LEG_COUNT}): {leg_census_total} "
        + ("(printed: boxed knob)" if geo35 else "(== 0 ASSERTED)"))

    ring_ports = _ring_ports()
    assert len(ring_ports) == 2 * LEG_COUNT
    tags_f = _interface_facet_tags(
        msh, cell_tags,
        {SHEET_IFACE + i: (PORT_LOWER + i, PORT_UPPER + i) for i in ring_ports},
    )
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    chord = float(diag["ring_port_layout"]["ring_port_gap_chord_m"])
    areas = {i: _interface_area_or_zero(msh, tags_f, SHEET_IFACE + i, comm)
             for i in ring_ports}
    specs = []
    phi_hats = {}
    for i in ring_ports:
        phi_hat, _centre = _ring_gap_frame(i, LEG_COUNT)
        phi_hats[i] = np.asarray(phi_hat, dtype=float)
        specs.append(
            LumpedSheetPortSpec(
                port_id=f"P{i}",
                facet_tag=int(SHEET_IFACE + i),
                port_impedance_ohm=TERMINATED_PORT_IMPEDANCE_OHM,
                gap_height_m=chord,
                sheet_width_m=areas[i] / chord,
                drive_direction=tuple(float(c) for c in phi_hat),
                drive_voltage_v=1.0 + 0.0j,
                interior=True,
            )
        )
    # Measured geometry for the class table and the orientation control.
    azimuth = {i: _sheet_azimuth_deg(msh, tags_f, SHEET_IFACE + i, comm)
               for i in ring_ports}
    centre_z = {i: _sheet_centre_z(msh, tags_f, SHEET_IFACE + i, comm)
                for i in ring_ports}
    long_extent = {
        i: _sheet_extent_along(msh, tags_f, SHEET_IFACE + i, phi_hats[i], comm)
        for i in ring_ports
    }

    problem = TimeHarmonicProblem(
        mesh=msh,
        frequency_hz=FREQUENCY_HZ,
        material=HomogeneousMaterial(sigma=0.0, epsilon_r=1.0, mu_r=1.0),
        cell_tags=cell_tags,
        material_map={
            CONDUCTOR_CELL_TAG: HomogeneousMaterial(
                sigma=SIGMA_WIRE_S_PER_M, epsilon_r=1.0, mu_r=1.0
            ),
            PHANTOM_CELL_TAG: HomogeneousMaterial(
                sigma=SALINE_SIGMA, epsilon_r=SALINE_EPSILON_R, mu_r=1.0
            ),
        },
        boundary_condition="pec_zero_tangential_a",
    )
    omega = 2.0 * np.pi * FREQUENCY_HZ
    ctx = {
        "comm": comm,
        "msh": msh,
        "tags_f": tags_f,
        "cell_tags": cell_tags,
        "omega": omega,
        "specs": specs,
        "solver": TimeHarmonicSolver(problem, degree=degree),
    }
    # `PORT-16`'s helpers read the fixture through this dict shape.
    sweep = {"mesh": msh, "facet_tags": tags_f, "cell_tags": cell_tags, "specs": specs}

    # ---- the 32 drives ----------------------------------------------------
    all_ids = [f"P{i}" for i in ring_ports]
    driven_ids = all_ids if identity_columns is None else all_ids[:identity_columns]
    columns: dict[str, dict] = {}
    identity: dict[str, dict] = {}
    identity_drives = set(range(len(driven_ids)))
    identity_mode = ("all 32 drives" if identity_columns is None else
                     f"drives 1..{identity_columns} ({', '.join(driven_ids)}), "
                     "FEM_EM_WF7_IDENTITY_COLUMNS, stop rule bypassed")
    for k, driven in enumerate(driven_ids):
        col = _solve_one_drive(ctx, driven, reuse_factorization=(k > 0))
        fields = col.pop("fields")
        if k in identity_drives:
            t_id = time.perf_counter()
            sheets = [spec.sheet(driven=(spec.port_id == driven)) for spec in specs]
            driven_sheet = next(s for s in sheets if s.port_id == driven)
            e = fields.e_complex
            p_src = _source_power_w(sweep, driven_sheet, e, omega)
            per_sheet = {
                s.port_id: float(_sheet_field_dissipation_w(sweep, s, e, omega))
                for s in sheets
            }
            p_sheet = float(sum(per_sheet.values()))
            # Negative control: the undriven sheet of smallest dissipation dropped.
            omit_id = min((pid for pid in per_sheet if pid != driven),
                          key=lambda pid: per_sheet[pid])
            # P_vol is `_loss_power_w`'s ½∫σ|E|² over phantom + conductor, which
            # `_solve_one_drive` already evaluated through the same `mean_sar`.
            p_vol = float(col["phantom"] + col["conductor"])
            deficit = terminal_form_deficit(msh, tags_f, sheets, e, comm)
            identity[driven] = {
                "p_src": p_src,
                "p_vol": p_vol,
                "p_sheet": p_sheet,
                "rel": abs(p_src - p_vol - p_sheet) / abs(p_src),
                "pooled": float(deficit["pooled"]),
                "omit_id": omit_id,
                "omit_w": per_sheet[omit_id],
                "neg_rel": abs(p_src - p_vol - (p_sheet - per_sheet[omit_id]))
                / abs(p_src),
            }
            id_s = time.perf_counter() - t_id
            if k == 0 and identity_columns is not None:
                say(f"identity price, drive 1: {id_s:.2f} s (stop rule bypassed, "
                    "FEM_EM_WF7_IDENTITY_COLUMNS)")
            if k == 0 and identity_columns is None:
                say(f"identity price, drive 1: {id_s:.2f} s (stop rule "
                    f"{IDENTITY_PRICE_STOP_S:.0f} s)")
                projected = (time.perf_counter() - t_window
                             + 31 * (HELD_SOLVE_ESTIMATE_S + id_s)
                             + CONTROL_BUILD_ESTIMATE_S)
                say(f"projected window with the identity on all 32: {projected:.0f} s "
                    f"(budget {WINDOW_BUDGET_S:.0f} s)")
                if id_s > IDENTITY_PRICE_STOP_S or projected > WINDOW_BUDGET_S:
                    identity_drives = set(IDENTITY_FALLBACK_DRIVES)
                    identity_mode = (
                        "drives 1, 9, 17, 25 only (first identity cost "
                        f"{id_s:.2f} s vs stop rule {IDENTITY_PRICE_STOP_S:.0f} s; "
                        f"projected window {projected:.0f} s vs budget "
                        f"{WINDOW_BUDGET_S:.0f} s)"
                    )
                    say(f"identity restricted to {identity_mode}")
        del fields  # 32 columns of F-human fields is the memory trap
        columns[driven] = col
        line = (f"drive {k + 1:2d}/{len(driven_ids)} {driven}: {col['solve_kind']} "
                f"{col['solve_time']:.2f} s  S_driven {col['s_column'][driven]:.6f}  "
                f"terminal residual {col['residual']:.3e}  "
                f"sum|S_ij|^2 {sum(abs(s) ** 2 for s in col['s_column'].values()):.6f}")
        if driven in identity:
            d = identity[driven]
            line += (f"  | P_src,exact {d['p_src']:.9e} W  P_vol {d['p_vol']:.9e} W  "
                     f"P_sheet,exact {d['p_sheet']:.9e} W  identity rel {d['rel']:.3e}  "
                     f"pooled C/terminal-1 {d['pooled']:.6e}  | NEG omit "
                     f"{d['omit_id']} ({d['omit_w']:.3e} W): rel {d['neg_rel']:.3e} = "
                     f"{d['neg_rel'] / DISCRETE_IDENTITY_RTOL:.3e}x the band")
        say(line + f"  t={time.perf_counter() - t_window:.1f} s")

    kinds = [columns[p]["solve_kind"] for p in driven_ids]
    say(f"solve kinds: {kinds.count('solve')} factorising, {kinds.count('held')} held; "
        f"summed ru_maxrss {_rss_summed_gib(comm):.3f} GiB")
    del ctx, problem

    # ---- (v), the reproduction and the negative control (POST-7 step 1) ----
    id_rels = [identity[p]["rel"] for p in identity]
    pooled = [identity[p]["pooled"] for p in identity]
    resid = [columns[p]["residual"] for p in driven_ids]
    neg_min = min(identity[p]["neg_rel"] for p in identity)
    say(f"(v) exact discrete power identity on {identity_mode}: max rel "
        f"{max(id_rels):.3e} (band {DISCRETE_IDENTITY_RTOL:g}, ASSERTED); "
        f"pooled C/terminal-1 {min(pooled):.6e}..{max(pooled):.6e} beside the "
        f"terminal residual {min(resid):.3e}..{max(resid):.3e} (printed), degree {degree}")
    say(f"NEG (rule (e)): one undriven sheet's dissipation omitted, min rel "
        f"{neg_min:.3e} = {neg_min / DISCRETE_IDENTITY_RTOL:.3e}x the band (bar "
        f"{IDENTITY_NEGATIVE_CONTROL_FACTOR:.0f}x, predicted >= 460x; "
        + ("ASSERTED" if identity_columns is not None else "printed") + ")")
    repro = None
    if "P17" in identity:
        d17 = identity["P17"]
        repro = (abs(d17["p_src"] / P17_P_SRC_EXACT_RECORD_W - 1.0),
                 abs(d17["pooled"] / P17_POOLED_KAPPA_RECORD - 1.0))
        say(f"reproduction P17: P_src,exact {d17['p_src']:.9e} W vs record "
            f"{P17_P_SRC_EXACT_RECORD_W:.9e} W, rel {repro[0]:.3e}; pooled kappa "
            f"{d17['pooled']:.6e} vs record {P17_POOLED_KAPPA_RECORD:.6e}, rel "
            f"{repro[1]:.3e} (band {REUSE_REPRODUCTION_RTOL:g}, "
            + ("printed only: boxed knob, records are leg-box-free" if geo35 else
               "ASSERTED, degree 1" if degree == 1 else
               f"printed only at degree {degree}: the record is degree-1") + ")")
        c2_p17 = abs(d17["p_src"] / P17_P_SRC_EXACT_RECORD_W_BOXED - 1.0)
        say(f"(c2) P17 P_src,exact vs the BOXED record {P17_P_SRC_EXACT_RECORD_W_BOXED:.9e} W: "
            f"rel {c2_p17:.3e} = {c2_p17 / REUSE_REPRODUCTION_RTOL:.1f}x the band; pooled "
            f"kappa vs {P17_POOLED_KAPPA_RECORD_BOXED:.6e} rel "
            f"{abs(d17['pooled'] / P17_POOLED_KAPPA_RECORD_BOXED - 1.0):.3e} "
            + ("(printed: boxed knob)" if geo35 else
               "(> 1x ASSERTED at degree 1; predicted ~77x)"))
    if not (identity_columns is None):
        say(f"(i) cells {n_cells}, relative {rel:.3e} to {F_HUMAN_BRANCH_B_CELL_RECORD} "
            f"(band {CELL_COUNT_BAND}, ASSERTED)")

    matrix = identity_columns is None
    if matrix:
        # ---- the 32x32 and its identities -------------------------------------
        s = np.array(
            [[columns[pj]["s_column"][pi] for pj in driven_ids] for pi in driven_ids],
            dtype=complex,
        )
        saved = None
        if save_path is not None:
            # ---- H15a: save first, so a later red or kill leaves the matrix ----
            t_save = time.perf_counter()
            sheets_geo = [{"ordinal": i, "tag": SHEET_IFACE + i,
                           "azimuth_deg": azimuth[i], "z": centre_z[i]}
                          for i in ring_ports]
            _az_ref, slots = _ring_quadrature_slots(sheets_geo)
            sigma_pre = float(np.linalg.norm(s, 2))
            err = None
            if comm.rank == 0:
                try:
                    os.makedirs(os.path.dirname(save_path), exist_ok=True)
                    np.savez(
                        save_path,
                        s32=np.asarray(s, dtype=np.complex128),
                        order=np.array(driven_ids),
                        slot_k=np.array([slots[p][0] for p in driven_ids]),
                        slot_sign=np.array([slots[p][1] for p in driven_ids]),
                        kappa=float(identity["P17"]["pooled"]) if "P17" in identity
                        else float("nan"),
                        cells=n_cells, freq_hz=FREQUENCY_HZ, degree=degree,
                        reciprocity=_reciprocity_ratio(s), sigma_max=sigma_pre,
                        kappa_applied=False,
                        pooled_kappa=np.array([identity[p]["pooled"] if p in identity
                                               else np.nan for p in driven_ids]),
                        ranks=comm.size,
                    )
                except Exception as exc:  # noqa: BLE001 — every rank must raise
                    err = f"{type(exc).__name__}: {exc}"
            err = comm.bcast(err, root=0)
            comm.Barrier()
            if err is not None:
                raise RuntimeError(f"FEM_EM_WF7_SAVE_NPZ write failed on rank 0: {err}")
            with np.load(save_path, allow_pickle=False) as z:
                loaded = {k: z[k] for k in z.files}
            save_s = time.perf_counter() - t_save
            s_l = loaded["s32"]
            l1 = {
                "dtype": s_l.dtype == np.complex128,
                "equal": bool(np.array_equal(s_l, s)),
                "order": [str(x) for x in loaded["order"]] == list(driven_ids),
                "slots": all(
                    sorted(loaded["slot_k"][loaded["slot_sign"] == sg].tolist())
                    == list(range(LEG_COUNT)) for sg in (1.0, -1.0)
                ),
            }
            az_l = np.array([azimuth[i] for i in ring_ports])
            zc_l = np.array([centre_z[i] for i in ring_ports])
            _c, spreads_l, _r = _class_means(s_l, az_l, zc_l)
            sigma_l = float(np.linalg.norm(s_l, 2))
            sigma_rel = abs(sigma_l / SIGMA_MAX_D1_RECORD - 1.0)
            class_rels = {k: abs(spreads_l[k][2] / CLASS_MEAN_D1_RECORD[k] - 1.0)
                          for k in CLASS_MEAN_D1_RECORD}
            worst_class_rel = max(class_rels.values())
            # Asserted control: P17 x CONTROL_COLUMN_SCALE on the reloaded array.
            ctl = s_l.copy()
            ctl[:, driven_ids.index("P17")] *= CONTROL_COLUMN_SCALE
            _c, spreads_c, _r = _class_means(ctl, az_l, zc_l)
            ctl_shift = _class_mean_worst_rel(spreads_c)
            # Predicted, printed (rule (e)): a single-precision save.
            _c, spreads_32, _r = _class_means(s_l.astype(np.complex64), az_l, zc_l)
            f32_shift = _class_mean_worst_rel(spreads_32)
            saved = {"l1": l1, "sigma_rel": sigma_rel, "class_rels": class_rels,
                     "worst_class_rel": worst_class_rel, "ctl_shift": ctl_shift}
            size_b = comm.bcast(os.path.getsize(save_path) if comm.rank == 0 else None,
                                root=0)
            say(f"H15a save: {save_path} ({size_b} B = {size_b / 1024:.1f} KiB), save + "
                f"reload {save_s:.3f} s, keys {sorted(loaded)}")
            say(f"(L1) round trip: s32 dtype {s_l.dtype} complex128 {l1['dtype']}, "
                f"array_equal {l1['equal']}, order equal {l1['order']}, 16 slots per "
                f"ring {l1['slots']} (ASSERTED)")
            tag = ("ASSERTED, degree 1" if degree == 1 else
                   f"printed only at degree {degree}: the record is degree-1")
            say(f"(L2) through the file: sigma_max {sigma_l:.9f} vs record "
                f"{SIGMA_MAX_D1_RECORD:.9f}, rel {sigma_rel:.3e}; worst class-mean rel "
                f"{worst_class_rel:.3e} (band {REUSE_REPRODUCTION_RTOL:g}, {tag})")
            for k in CLASS_MEAN_D1_RECORD:
                say(f"     {'same ' if k[0] else 'other'} ring, {k[1]} steps  mean |S| "
                    f"{spreads_l[k][2]:.9e} vs record {CLASS_MEAN_D1_RECORD[k]:.9e}, "
                    f"rel {class_rels[k]:.3e}")
            say(f"save control (asserted): P17 x {CONTROL_COLUMN_SCALE} on the reloaded "
                f"array moves the worst class-mean rel to {ctl_shift:.3e} = "
                f"{ctl_shift / REUSE_REPRODUCTION_RTOL:.1f}x the band (bar "
                f"{SAVE_CONTROL_FACTOR:.0f}x; predicted [2e-4, 4e-4])")
            say(f"save control (predicted, printed, never asserted — rule (e)): complex64 "
                f"copy, worst class-mean rel {f32_shift:.3e} (predicted ~1e-8..1e-7, "
                f"below the {REUSE_REPRODUCTION_RTOL:g} band)")
        ratio = _reciprocity_ratio(s)
        control = s.copy()
        control[:, 0] *= CONTROL_COLUMN_SCALE
        control_ratio = _reciprocity_ratio(control)
        sigma_max = float(np.linalg.norm(s, 2))
        col_norms = np.sum(np.abs(s) ** 2, axis=0)

        az = np.array([azimuth[i] for i in ring_ports])
        zc = np.array([centre_z[i] for i in ring_ports])
        n = len(ring_ports)
        classes: dict = {}
        worst_residual_deg = 0.0
        for i in range(n):
            for j in range(n):
                steps, residual_deg = _azimuth_class(az[i], az[j])
                worst_residual_deg = max(worst_residual_deg, residual_deg)
                classes.setdefault((bool(zc[i] * zc[j] > 0.0), steps), []).append(
                    float(abs(s[i, j]))
                )
        spreads = {}
        for key, mags in classes.items():
            m = np.array(mags)
            spreads[key] = (float((m.max() - m.min()) / m.mean()), len(m), float(m.mean()))
        worst = max(spreads, key=lambda key: spreads[key][0])

        say(f"(i) cells {n_cells}, relative {rel:.3e} to {F_HUMAN_BRANCH_B_CELL_RECORD} "
            f"(band {CELL_COUNT_BAND}, "
            + ("printed only: GEO-35 knob on" if geo35 else "ASSERTED") + ")")
        say(f"(ii) reciprocity ||S-S^T||_F/||S||_F = {ratio:.6e} (band "
            f"{RECIPROCITY_BAND:.0e}, ASSERTED); ||S||_F = {np.linalg.norm(s, 'fro'):.6f}")
        say(f"(iii) sigma_max = {sigma_max:.9f} (ceiling {MATRIX_PASSIVITY_CEILING}, "
            f"ASSERTED), margin {MATRIX_PASSIVITY_CEILING - sigma_max:+.3e}; full-column "
            f"sum|S_ij|^2 min {col_norms.min():.6f} max {col_norms.max():.6f} (printed)")
        say(f"(iv) C16 x mirror classes: {len(classes)} (expected {N_SYMMETRY_CLASSES}, "
            f"ASSERTED); worst azimuth residual {worst_residual_deg:.3e} deg "
            f"(< {AZIMUTH_MATCH_DEG:.0e}, ASSERTED); step {AZIMUTH_STEP_DEG:.3f} deg")
        for same_ring in (True, False):
            for steps in range(N_AZIMUTH_CLASSES):
                if (same_ring, steps) not in spreads:
                    say(f"     {'same ' if same_ring else 'other'} ring, {steps} steps: EMPTY")
                    continue
                sp, cnt, mean = spreads[(same_ring, steps)]
                say(f"     {'same ' if same_ring else 'other'} ring, {steps} steps  "
                    f"n = {cnt:4d}  mean |S| = {mean:.9e}  spread {sp * 100:8.4f}%"
                    + ("   <-- worst" if (same_ring, steps) == worst else ""))
        say(f"     worst class spread {spreads[worst][0]:.4e} (band "
            f"{OPPOSITE_SPREAD_BAND}, ASSERTED)")
        say(f"S_driven(P17) {columns['P17']['s_column']['P17']:.6f} beside the step-0 "
            f"record {S_DRIVEN_STEP0_RECORD} (printed, -n {comm.size})")
        say(f"control (asserted): column {driven_ids[0]} x {CONTROL_COLUMN_SCALE}: "
            f"ratio {control_ratio:.6e} = {control_ratio / max(ratio, 1e-300):.3e}x the "
            f"measured (bar {CONTROL_RATIO_LIFT:.0e}x)")
        bar = MATRIX_CONTROL_MARGIN * RECIPROCITY_BAND
        say(f"control (predicted, printed, never asserted — rule (e)): {control_ratio:.6e} "
            f"vs MATRIX_CONTROL_MARGIN x RECIPROCITY_BAND = {bar:.1e}: "
            f"{'at/above' if control_ratio >= bar else 'BELOW'} the bar "
            f"(predicted ~1.5-2.1e-3, at or past the fixture's ceiling)")
        # Printed on both routes: this run against the *other* fixture's records
        # (on the default route this is control (c2), asserted below).
        s_ref = S_DRIVEN_STEP0_RECORD if geo35 else S_DRIVEN_STEP0_RECORD_BOXED
        sig_ref = SIGMA_MAX_D1_RECORD if geo35 else SIGMA_MAX_D1_RECORD_BOXED
        cls_ref = CLASS_MEAN_D1_RECORD if geo35 else CLASS_MEAN_D1_RECORD_BOXED
        s17 = complex(columns["P17"]["s_column"]["P17"])
        s17_rec = complex(s_ref)
        say(f"GEO-35 move: S_driven(P17) {s17:.6f} vs {s_ref}, "
            f"rel |dS| {abs(s17 - s17_rec) / abs(s17_rec):.3e}; sigma_max "
            f"{sigma_max:.9f} vs {sig_ref:.9f}, rel "
            f"{abs(sigma_max / sig_ref - 1.0):.3e}")
        for k in cls_ref:
            say(f"GEO-35 move: {'same ' if k[0] else 'other'} ring, {k[1]} steps "
                f"mean |S| {spreads[k][2]:.9e} vs {cls_ref[k]:.9e}, "
                f"rel {abs(spreads[k][2] / cls_ref[k] - 1.0):.3e}")
        c2_class = max(abs(spreads[k][2] / cls_ref[k] - 1.0) for k in cls_ref)
        say(f"GEO-35 move: worst class-mean rel vs the "
            f"{'leg-box-free' if geo35 else 'BOXED'} records {c2_class:.3e} = "
            f"{c2_class / REUSE_REPRODUCTION_RTOL:.1f}x REUSE_REPRODUCTION_RTOL "
            + ("(printed)" if geo35 else "((c2) > 1x ASSERTED; predicted ~910x)"))
    else:
        say("(ii)-(iv), the arithmetic column-scale control and the rule-(e) matrix "
            f"control: SKIPPED (FEM_EM_WF7_IDENTITY_COLUMNS={identity_columns})")

    # ---- the mesh-level orientation control (no solve) --------------------
    default_extent = None
    if orientation == "inline":
        default_extent, default_cells, d_build_s = _default_orientation_extent(
            comm, ring_ports)
        say(f"orientation control: default build {default_cells} cells "
            f"({d_build_s:.2f} s); max extent along phi_hat "
            f"{max(default_extent.values()):.3e} m (< {DEGENERATE_EXTENT_M:.0e}, "
            f"ASSERTED); longitudinal chord {chord:.6e} m, min extent along phi_hat "
            f"{min(long_extent.values()):.6e} m (> {DEGENERATE_EXTENT_M:.0e}, ASSERTED)")
    else:
        say("orientation control: SPLIT — asserted by "
            "test_wf7_default_orientation_control in its own window; longitudinal "
            f"chord {chord:.6e} m, min extent along phi_hat "
            f"{min(long_extent.values()):.6e} m (> {DEGENERATE_EXTENT_M:.0e}, ASSERTED)")
    say(f"window elapsed {time.perf_counter() - t_window:.1f} s, summed ru_maxrss "
        f"{_rss_summed_gib(comm):.3f} GiB")

    # ---- the asserts, in the item's order -----------------------------------
    assert abs(rel) < CELL_COUNT_BAND, (
        f"(i) {n_cells} cells, relative {rel:.3e} to {F_HUMAN_BRANCH_B_CELL_RECORD}"
    )
    if not geo35:
        assert leg_census_total == 0, (
            f"(N') {leg_census_total} owned cells still carry a leg port-box tag"
        )
    if matrix:
        assert ratio <= RECIPROCITY_BAND, (
            f"(ii) reciprocity {ratio:.6e} > {RECIPROCITY_BAND}"
        )
        assert float(col_norms.max()) <= MATRIX_PASSIVITY_CEILING
        assert sigma_max <= MATRIX_PASSIVITY_CEILING, (
            f"(iii) sigma_max {sigma_max:.9f} > {MATRIX_PASSIVITY_CEILING}"
        )
        assert worst_residual_deg < AZIMUTH_MATCH_DEG, (
            f"(iv) a pair sits {worst_residual_deg:.3e} deg off every C16 class"
        )
        assert len(classes) == N_SYMMETRY_CLASSES, (
            f"(iv) {len(classes)} classes, not {N_SYMMETRY_CLASSES}"
        )
        assert sum(len(v) for v in classes.values()) == n * n
        assert spreads[worst][0] <= OPPOSITE_SPREAD_BAND, (
            f"(iv) class {worst} spreads {spreads[worst][0]:.4e} > {OPPOSITE_SPREAD_BAND}"
        )
    for pid, d in identity.items():
        assert d["p_src"] > 0.0, f"(v) {pid}: P_src,exact {d['p_src']:.6e} <= 0"
        assert d["rel"] <= DISCRETE_IDENTITY_RTOL, (
            f"(v) {pid}: identity misses by {d['rel']:.6e} — P_src {d['p_src']:.9e}, "
            f"P_vol {d['p_vol']:.9e}, P_sheet {d['p_sheet']:.9e}"
        )
    if degree == 1 and not geo35:
        assert repro is not None, "reproduction: P17 is not among the driven columns"
        assert repro[0] <= REUSE_REPRODUCTION_RTOL, (
            f"reproduction: P17 P_src,exact rel {repro[0]:.3e} > {REUSE_REPRODUCTION_RTOL}"
        )
        assert repro[1] <= REUSE_REPRODUCTION_RTOL, (
            f"reproduction: P17 pooled kappa rel {repro[1]:.3e} > {REUSE_REPRODUCTION_RTOL}"
        )
        # (c2): the 1e-6 comparator tells the meshes apart (rule (e), asserted).
        assert c2_p17 > REUSE_REPRODUCTION_RTOL, (
            f"(c2) P17 P_src,exact within {c2_p17:.3e} of the BOXED record: the "
            "switch did not take effect"
        )
        if matrix:
            assert c2_class > REUSE_REPRODUCTION_RTOL, (
                f"(c2) class means within {c2_class:.3e} of the BOXED records"
            )
    if identity_columns is not None:
        assert neg_min >= IDENTITY_NEGATIVE_CONTROL_FACTOR * DISCRETE_IDENTITY_RTOL, (
            f"NEG: one sheet's dissipation omitted misses only by {neg_min:.3e}"
        )
    assert chord > DEGENERATE_EXTENT_M
    assert min(long_extent.values()) > DEGENERATE_EXTENT_M
    if default_extent is not None:
        assert max(default_extent.values()) < DEGENERATE_EXTENT_M, (
            "the default-orientation ring sheets have extent along phi_hat — the "
            "mis-paired orientation control does not hold"
        )
    if matrix:
        assert control_ratio > CONTROL_RATIO_LIFT * ratio, (
            f"a 1% column scale lifts the reciprocity ratio only to {control_ratio:.6e} "
            f"from {ratio:.6e}"
        )
    if matrix and saved is not None:
        assert all(saved["l1"].values()), f"(L1) round trip failed: {saved['l1']}"
        if degree == 1:
            assert saved["sigma_rel"] <= REUSE_REPRODUCTION_RTOL, (
                f"(L2) sigma_max through the file rel {saved['sigma_rel']:.3e}"
            )
            assert saved["worst_class_rel"] <= REUSE_REPRODUCTION_RTOL, (
                f"(L2) class means through the file: {saved['class_rels']}"
            )
        assert saved["ctl_shift"] >= SAVE_CONTROL_FACTOR * REUSE_REPRODUCTION_RTOL, (
            f"save control: P17 x {CONTROL_COLUMN_SCALE} moves the class means only "
            f"{saved['ctl_shift']:.3e}"
        )
