"""ports:21 -- the capacitor ESR sheet on the copper birdcage (EX-105).

What a 0.2 ohm ESR does to the copper birdcage's loss partition.  On the
F-small exact-rectangle hole mesh at 64 MHz (degree 1, kappa off, copper
Leontovich), each lumped capacitor sheet is a series-RC impedance,
``LumpedSheetPortSpec(port_impedance_ohm=series_rlc_impedance(f, r_ohm=r,
l_h=0, c_f=C_F))`` -- the operator's ESR ruling of 2026-10-06 made runnable.

Nothing here is new physics and nothing is restated: the three arms are the
``PORT-24`` step 2 gate's own, called through the lifted module-level
``_run_port24_arm`` of ``tests/validation/test_port24_series_rc_sheet.py``,
and every band and constant is imported from the modules the gate imports it
from.  The script re-makes the gate's assertions on each returned dict:

* (I)  every sheet's field-solved 3x3 equals the circuit reduction with that
  sheet in its series-RC, at ``REDUCTION_BAND`` (arms r02 and r0);
* (Z)  at r = 0 the sheet reproduces the capacitor-only sheet, at
  ``LEG_D0_REPRODUCTION_BAND`` (arm r0);
* (P)  the exact power identity, the RC sheets' own dissipation counted, at
  ``DISCRETE_IDENTITY_RTOL`` (every arm);
* (R0) the r = 0 sheets dissipate exactly 0 W (arm r0);
* (Pl) negative control: a 0.4 ohm plant on P2 breaks the 0.2 ohm identity
  there and only there (arm plant) -- direction only, no ceiling.

Printed, never asserted: the ESR loss share, where it sits, the copper /
phantom split, and the loss-partition reading Q(r)/Q(0) beside
1/Q = 1/Q_rest + 1/Q_C.  Q_rest is copper + phantom (not Q_Cu alone) and is
inferred from the same partition, so that agreement is true by construction.

Not claimed: an absolute Q, f0, S11, Z_in or tuning figure; a measured ESR.

Needs the complex DolfinX build; the runner sources it for ``ports:``::

    ./run_examples.sh -e ports:21 -n 2 -t 590
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.lumped import series_rlc_impedance  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.validation.test_port24_series_rc_sheet as gate  # noqa: E402
from tests.validation.test_port24_series_rc_sheet import (  # noqa: E402
    DEGREE,
    ESR_OHM,
    FREQUENCY_HZ,
    PLANT_INDEX,
    PLANT_OHM,
    POWER_DRIVEN,
    _rc,
    _run_port24_arm,
)
from tests.validation.test_ans7_exact_strip import C_F  # noqa: E402
from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL  # noqa: E402
from tests.validation.test_port_birdcage_four_port import LEG_D0_REPRODUCTION_BAND  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG  # noqa: E402
from tests.validation.test_port_lumped_rlc_termination import REDUCTION_BAND  # noqa: E402
from tests.validation.test_port_package_sparameters import REFERENCE_IMPEDANCE_OHM  # noqa: E402

# The bands are the gate's own objects, not copies.
assert gate.REDUCTION_BAND is REDUCTION_BAND
assert gate.LEG_D0_REPRODUCTION_BAND is LEG_D0_REPRODUCTION_BAND
assert gate.DISCRETE_IDENTITY_RTOL is DISCRETE_IDENTITY_RTOL
assert gate.C_F is C_F
assert gate.REFERENCE_IMPEDANCE_OHM is REFERENCE_IMPEDANCE_OHM

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "ports_21_capacitor_esr_sheet"
AIR_CELL_TAG = 2   # display input only (ports:18 :222, :230), not a record

FIELD_KEYS = ("msh", "cell_tags", "facet_tags", "specs", "e")


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _e_magnitude_dg0(e_complex, msh):
    """``|E|`` as a real DG0 field (the ``ports:18`` construction)."""
    dg0 = fem.functionspace(msh, ("DG", 0))
    e_sq = fem.Function(dg0)
    e_sq.interpolate(
        fem.Expression(ufl.inner(e_complex, e_complex), dg0.element.interpolation_points)
    )
    e_mag = fem.Function(dg0, name="E_magnitude")
    e_mag.x.array[:] = np.sqrt(np.maximum(np.real(e_sq.x.array), 0.0))
    e_mag.x.scatter_forward()
    return e_mag


def _assert_sheet_identity(res):
    """(I), plus (Z) and (R0) at r = 0, on a returned arm dict."""
    arm = res["arm"]
    rows = res["ident"]
    assert len(rows) == 4
    for row in rows:
        # printed Z column must be the gate's own circuit element
        assert row["z"] == _rc(row["r"])
        assert row["z"] == series_rlc_impedance(FREQUENCY_HZ, r_ohm=row["r"], l_h=0.0, c_f=C_F)
    if arm == "plant":
        return
    for row in rows:
        assert row["res"] <= REDUCTION_BAND, (arm, row["i"], row["res"])
    if arm == "r0":
        for row in rows:
            assert row["cap_res"] <= LEG_D0_REPRODUCTION_BAND, (row["i"], row["cap_res"])
        assert res["p_rc"] == 0.0, res["p_rc"]


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "This example needs the complex DolfinX build: "
            "source /usr/local/bin/dolfinx-complex-mode (the runner does this "
            "automatically for the `ports:` group)."
        )
    _say("=" * 78)
    _say("EX-105 -- the capacitor ESR sheet on the copper birdcage, F-small, 64 MHz")
    _say("=" * 78)

    # ---- arm r02: every capacitor sheet at the ESR; fields kept for the outputs
    r02 = _run_port24_arm("r02", return_fields=True)
    _assert_sheet_identity(r02)
    assert r02["resid"] <= DISCRETE_IDENTITY_RTOL, r02["resid"]

    msh = r02["msh"]
    cell_tags = r02["cell_tags"]
    facet_tags = r02["facet_tags"]
    OUTPUT_DIR.mkdir(exist_ok=True)
    FIGURE_DIR.mkdir(exist_ok=True)

    region_names = {AIR_CELL_TAG: "air", PHANTOM_CELL_TAG: "phantom"}
    for spec in r02["specs"]:
        if spec.port_id == POWER_DRIVEN:
            region_names[int(spec.facet_tag)] = f"{spec.port_id} 50 ohm feed"
        else:
            region_names[int(spec.facet_tag)] = f"{spec.port_id} series-RC r={ESR_OHM:g} ohm"
    write_setup_figure(
        msh,
        cell_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names=region_names,
        hide_tags=(AIR_CELL_TAG,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title=(f"ports:21 -- F-small copper birdcage, exact-rectangle sheets: P1 50 ohm feed, "
               f"P2-P4 series-RC (r = {ESR_OHM:g} ohm, C_F)"),
        comm=comm,
    )

    e_mag = _e_magnitude_dg0(r02["e"], msh)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        cell_tags,
        {"E_magnitude": e_mag},
        comm=comm,
        facet_tags=facet_tags,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(f"\n[output] {path} -- CellTags, E_magnitude (DG0), facet_tags grid "
         f"(sheet tags and tag 401); r02 arm, P1 drive")

    for k in FIELD_KEYS:
        r02.pop(k)
    del msh, cell_tags, facet_tags, e_mag

    # ---- arm r0: the lossless capacitor-only control for the same sheets
    r0 = _run_port24_arm("r0", return_fields=True)
    _assert_sheet_identity(r0)
    assert r0["resid"] <= DISCRETE_IDENTITY_RTOL, r0["resid"]
    for k in FIELD_KEYS:
        r0.pop(k)

    # ---- arm plant: negative control, a 0.4 ohm P2 held against the 0.2 ohm identity
    plant = _run_port24_arm("plant")
    _assert_sheet_identity(plant)
    planted = plant["ident"][PLANT_INDEX]
    others = [r for r in plant["ident"] if r["i"] != PLANT_INDEX]
    unplanted_max = max(r["res"] for r in others)
    assert planted["res"] > REDUCTION_BAND, planted["res"]
    for r in others:
        assert r["res"] <= REDUCTION_BAND, (r["i"], r["res"])
    assert plant["resid"] <= DISCRETE_IDENTITY_RTOL, plant["resid"]

    # ---- the summary: arithmetic on returned values, printed only
    ps = r02["p_sheet"]
    p_rc_sheets = {k: v for k, v in ps.items() if k != POWER_DRIVEN}
    q_ratio, q_rest, q_c = r02["q_ratio"], r02["q_rest"], r02["q_c"]
    q_r = 1.0 / (1.0 / q_rest + 1.0 / q_c)
    p_rest_02, p_rest_0 = r02["p_rest"], r0["p_rest"]

    _say("")
    _say("[ESR on the copper birdcage] (printed; Q_rest = copper + phantom, not Q_Cu)")
    _say(f"  f = {FREQUENCY_HZ:.6e} Hz, degree {DEGREE}, kappa off, C_F = {C_F!r} F, "
         f"Z_0 = {float(REFERENCE_IMPEDANCE_OHM):g} ohm, mesh {r02['cells']} cells, "
         f"-n {comm.size}")
    for row in r02["ident"]:
        z = row["z"]
        _say(f"  (I) r02 sheet P{row['i'] + 1}: Z = {z.real:+.9e}{z.imag:+.9e}j ohm, "
             f"residual {row['res']:.6e} = {row['res'] / REDUCTION_BAND:.4f} x REDUCTION_BAND "
             f"(ASSERTED <= band)")
    for row in r0["ident"]:
        _say(f"  (I) r0  sheet P{row['i'] + 1}: residual {row['res']:.6e} "
             f"(ASSERTED <= band); (Z) vs capacitor-only {row['cap_res']:.6e} "
             f"(ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e}, round-off, compared at the band only)")
    _say(f"  (P) power identity residual/P_src: r02 {r02['resid']:.3e}, r0 {r0['resid']:.3e}, "
         f"plant {plant['resid']:.3e} (ASSERTED <= {DISCRETE_IDENTITY_RTOL:g} each)")
    _say(f"  (R0) r = 0 sheets dissipate {r0['p_rc']:.9e} W (ASSERTED == 0)")
    _say(f"  (Pl) planted P{PLANT_INDEX + 1} at {PLANT_OHM:g} ohm held against the "
         f"{ESR_OHM:g} ohm identity: residual {planted['res']:.6e} = "
         f"{planted['res'] / REDUCTION_BAND:.2f} x REDUCTION_BAND (ASSERTED > band); "
         f"unplanted max {unplanted_max:.6e} (ASSERTED <= band); no ceiling asserted")
    _say(f"  ESR loss share P_RC/(P_RC + P_vol + P_surf) = {r02['share']:.6e} "
         f"({r02['share'] * 100:.1f} %)")
    _say("  where the ESR loss sits (P1 drive): "
         + ", ".join(f"{k} {v:.9e} W ({v / r02['p_rc'] * 100:.4f} % of P_RC)"
                     for k, v in p_rc_sheets.items()))
    _say(f"  the rest: P_vol (phantom) {r02['p_vol']:.9e} W, P_surf (copper, tag 401) "
         f"{r02['p_surf']:.9e} W, copper share of P_rest = "
         f"{r02['p_surf'] / p_rest_02 * 100:.1f} %")
    _say(f"  Q reading: Q(r)/Q(0) = {q_ratio:.6e}; Q_C = 1/(omega C r) = {q_c:.6e}; "
         f"Q_rest = {q_rest:.6e}")
    _say(f"  closed form 1/Q = 1/Q_rest + 1/Q_C => Q(r) = {q_r:.4f}; "
         f"Q(r)/Q_rest = {q_r / q_rest:.6e} (true by construction, printed only)")
    _say(f"  fixed-stored-energy size: P_rest(r={ESR_OHM:g})/P_rest(r=0) = "
         f"{p_rest_02:.9e}/{p_rest_0:.9e} = {p_rest_02 / p_rest_0:.6f}")
    _say(f"  (Dz) RC dissipation omitted: residual/P_src {r02['resid_dz']:.6e} at r = {ESR_OHM:g} "
         f"vs {r0['resid_dz']:.3e} at r = 0 (printed; P_src is dominated by P1's own "
         f"{float(REFERENCE_IMPEDANCE_OHM):g} ohm sheet, {ps[POWER_DRIVEN]:.9e} of "
         f"{r02['p_src']:.9e} W)")

    elapsed = time.perf_counter() - started
    _say(f"\n[timing] total {elapsed:.1f} s at -n {comm.size}")
    _say("\nAll imported gates hold.")


if __name__ == "__main__":
    main()
