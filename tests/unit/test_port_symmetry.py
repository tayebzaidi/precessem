"""Unit: symmetry-plane port rescaling (`TH-16` step 2).

R = Z_p·w/h and E_src = V/h are invariant under either cut (to 1e-15); a cut
whose kind disagrees with ĥ raises; ``pmc`` then ``pec`` is the quarter.
"""

import pytest

from fem_em_solver.ports.lumped import LumpedSheetPortSpec, sheet_resistivity_ohm_per_square
from fem_em_solver.ports.symmetry import full_impedance, halve_sheet_spec

X_NORMAL = (1.0, 0.0, 0.0)
Y_NORMAL = (0.0, 1.0, 0.0)
INVARIANCE_RTOL = 1e-15


def _spec():
    return LumpedSheetPortSpec(
        port_id="P1",
        facet_tag=10,
        port_impedance_ohm=50.0 + 0.0j,
        gap_height_m=0.16,
        sheet_width_m=0.16,
        drive_direction=(0.0, 1.0, 0.0),
        drive_voltage_v=1.0 + 0.0j,
    )


def _r(s):
    return sheet_resistivity_ohm_per_square(
        s.port_impedance_ohm, gap_height_m=s.gap_height_m, sheet_width_m=s.sheet_width_m
    )


def _esrc(s):
    return complex(s.drive_voltage_v) / s.gap_height_m


@pytest.mark.parametrize("cut, normal", [("pmc", X_NORMAL), ("pec", Y_NORMAL)])
def test_sheet_law_invariant(cut, normal):
    full = _spec()
    half = halve_sheet_spec(full, cut, plane_normal=normal)
    assert abs(_r(half) / _r(full) - 1.0) <= INVARIANCE_RTOL
    assert abs(_esrc(half) / _esrc(full) - 1.0) <= INVARIANCE_RTOL


def test_wrong_orientation_raises():
    with pytest.raises(ValueError):
        halve_sheet_spec(_spec(), "pmc", plane_normal=Y_NORMAL)
    with pytest.raises(ValueError):
        halve_sheet_spec(_spec(), "pec", plane_normal=X_NORMAL)
    with pytest.raises(ValueError):
        halve_sheet_spec(_spec(), "pec", plane_normal=(1.0, 1.0, 0.0))


def test_pmc_then_pec_is_quarter():
    full = _spec()
    quarter = halve_sheet_spec(
        halve_sheet_spec(full, "pmc", plane_normal=X_NORMAL), "pec", plane_normal=Y_NORMAL
    )
    assert quarter.sheet_width_m == full.sheet_width_m / 2.0
    assert quarter.gap_height_m == full.gap_height_m / 2.0
    assert quarter.port_impedance_ohm == full.port_impedance_ohm
    assert quarter.drive_voltage_v == full.drive_voltage_v / 2.0
    assert abs(_r(quarter) / _r(full) - 1.0) <= INVARIANCE_RTOL
    z = 12.5 - 3.0j
    assert full_impedance(full_impedance(z, "pec"), "pmc") == z
