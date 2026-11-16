"""Example (`EX-103`): the absorbing boundary on a small box -- a reflector at 10 MHz.

`PORT-1`'s two tori (tag 1 driven, 10 MHz, degree 1) sit in an air box whose
outer wall (facet tag 1) is solved three ways: PEC, natural, and the
first-order ABC (Jin 3e 9.2 (9.44), s = 1, `first_order_abc_term` through
`TimeHarmonicSolver.solve(extra_bilinear_terms=...)`).  At k0 * half-width
of 0.02-0.04 the ABC is not an open boundary: it acts as the natural
(magnetic) wall, so its box term B = d - D_inf^PEC is positive where PEC's is
negative.  The 1/3 PEC + 2/3 natural post-solve mix cancels that term.

**It asserts, and it does not re-implement** (the `ANS-1` rule).  The solver
call, the paddings, the weights, the divisor and every band are imported from
`tests/validation/test_th5_abc_two_torus.py` (`TH-5` step 2a) and
`tests/validation/test_th5_mixed_wall_two_torus.py` (`TH-5` step 2b).  The 2a
gate gained, additively, ``_solve_z_column(..., return_fields=False)``; the
gate was re-run green in the same slot.  2b's own module was NOT re-run: its
(M) / (c2) clauses are evaluated here from this run's PEC and natural 0.08 m
solves, after asserting both have the 119 533 cells the 2b window recorded.

Asserted: (P0), (E), (S), (c1), (M), (c2).  Printed beside the gate record
and never asserted: r(W) = |B^ABC| / |B^PEC| ((R) was retired by the
2026-10-07 weekly), the power per wall, (c3), the ABC free-exponent fit.

Run it::

    ./scripts/run_examples.sh -e th:16

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_16_abc_two_torus_{pec,abc,natural}_combined.xdmf``
(CG1 ``E_real``, ``E_imag``; facet tags in the same file).
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

from fem_em_solver.io.mesh import (  # noqa: E402
    PORT1_PADDING_RUNGS_M,
    _free_exponent_fit,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation import test_th5_abc_two_torus as gate_a  # noqa: E402
from tests.validation import test_th5_mixed_wall_two_torus as gate_b  # noqa: E402
from tests.validation.test_geo4_box_padding_heuristic import (  # noqa: E402
    P0_BAND_PP,
    P0_RUNG_PP,
)
from tests.validation.test_port_reaction_impedance import (  # noqa: E402
    FREQUENCY_HZ,
    IDENTITY_TOLERANCE,
    MAJOR_RADIUS,
    SEPARATION,
    TORUS_TAGS,
)
from tests.validation.test_th5_abc_two_torus import (  # noqa: E402
    ABC_PADDINGS,
    ORIGINAL_FALL_FACTOR,
    OUTER_BOUNDARY_TAG,
    PREMISE_PADDING,
    _b_pec,
    _d_inf_pec,
    _solve_z_column,
)
from tests.validation.test_th5_mixed_wall_two_torus import (  # noqa: E402
    MIX_BAR_DIVISOR,
    W_NAT,
    W_PEC,
    _mix_deficit,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_16_abc_two_torus"

#: Cell count both 2b-window solves recorded (20261008T183455Z_TH-5-step2c-W2.log:238-239);
#: asserted on this run's PEC and natural meshes before (M) / (c2) use them.
CELLS_008 = 119533

#: Gate records -- printed beside this run, a difference reported and never asserted.
#: (log lines are 20261008T183105Z_TH-5-step2c-W1.log unless marked W2 =
#: 20261008T183455Z_TH-5-step2c-W2.log)
RECORD = {
    "p0_delta_pp": 0.0080,                              # W1:928
    "re_z11": (3.784205e-03, 2.620085e-03, 1.923719e-03),   # W1:932-934
    "d_abc": (6.4067, 4.0899, 2.8190),                  # W1:938-940
    "b_abc": (4.7133, 2.3965, 1.1256),                  # W1:938-940
    "b_pec_fit": (-9.7258, -6.7190, -4.9667),           # W1:938-940
    "r": (0.4846, 0.3567, 0.2266),                      # W1:938-940
    "cells": (119533, 135685, 153952),                  # W1:232,574,757
    "power_abc": 3.553290e-03,                          # W1:932
    "power_pec": 4.171602e-37,                          # W1:232
    "power_nat": 2.581840e-01,                          # W1:927
    "c1_sep_pp": 14.4391,                               # W1:941
    "c3_pp": 0.224456,                                  # W1:942
    "dnat_minus_dpec_pp": 14.2227,                      # W2:605
    "b_pec_solve": -9.7338,                             # W2:604
    "b_nat": 4.4889,                                    # W2:604
    "r_nat": 0.4612,                                    # W2:604
    "b_mix": -0.2520,                                   # W2:606
    "bar": 1.9468,                                      # W2:606
    "b_swap": -4.9929,                                  # W2:607
    "swap_ratio": 2.5647,                               # W2:607
    "b_half": -2.6225,                                  # W2:607
    "w_star": 0.3156,                                   # W2:608
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _export(name, run, comm):
    """CG1 interpolants of Re E and Im E (interpolants: nothing asserted on them)."""
    msh, e = run["msh"], run["e"]
    v_vec = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    e_real = fem.Function(v_vec, name="E_real")
    e_real.interpolate(fem.Expression(ufl.real(e), v_vec.element.interpolation_points))
    e_imag = fem.Function(v_vec, name="E_imag")
    e_imag.interpolate(fem.Expression(ufl.imag(e), v_vec.element.interpolation_points))
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_{name}_combined", msh, run["cell_tags"],
        {"E_real": e_real, "E_imag": e_imag}, comm=comm, facet_tags=run["facet_tags"],
    )
    return path


def _release(run):
    for k in ("e", "msh", "cell_tags", "facet_tags"):
        run.pop(k, None)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "This example needs the complex DolfinX build: "
            "source /usr/local/bin/dolfinx-complex-mode (the runner does this "
            "automatically for the `th:` group)."
        )
    # One band object, two importers.
    assert gate_a.P0_BAND_PP is gate_b.P0_BAND_PP is P0_BAND_PP

    _say(comm, "=" * 78)
    _say(comm, "EX-103 -- the absorbing boundary on PORT-1's two-torus box (10 MHz)")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] P0_BAND_PP {P0_BAND_PP}, IDENTITY_TOLERANCE {IDENTITY_TOLERANCE}, "
               f"MIX_BAR_DIVISOR {MIX_BAR_DIVISOR:g}, weights PEC {W_PEC:.4f} / NAT {W_NAT:.4f}, "
               f"paddings {ABC_PADDINGS} m, f = {FREQUENCY_HZ:g} Hz, d = {SEPARATION} m, "
               f"R = {MAJOR_RADIUS} m (all imported from the TH-5 gate modules)")
    assert tuple(ABC_PADDINGS) == tuple(PORT1_PADDING_RUNGS_M)
    assert len(TORUS_TAGS) == 2

    # ---- the five solves, the gate fixture's own sequence --------------------
    outputs = []
    pec = _solve_z_column(PREMISE_PADDING, "pec", return_fields=True)
    outputs.append(_export("pec", pec, comm))
    _release(pec)  # keeps scalars (and omega_m12) only

    abc = {}
    for w in ABC_PADDINGS:
        first = w == PREMISE_PADDING
        abc[w] = _solve_z_column(w, "abc", return_fields=first)
        if first:
            # ---- setup figure (EX-57): the ABC 0.08 m mesh ------------------
            write_setup_figure(
                abc[w]["msh"], abc[w]["cell_tags"], FIGURE_DIR / f"{BASENAME}_setup.png",
                region_names={1: "torus 1 (driven)", 2: "torus 2 (receive)",
                              3: "air box (ABC on the outer wall)"},
                translucent_tags=(3,),
                slice_normal=(0.0, 1.0, 0.0),
                slice_origin=(0.0, 0.0, 0.0),
                title="th:16 -- two tori in a 0.08 m-padded air box, 10 MHz, "
                      "first-order ABC on the outer wall; slice y = 0",
                comm=comm,
            )
            outputs.append(_export("abc", abc[w], comm))
            _release(abc[w])

    nat = _solve_z_column(PREMISE_PADDING, "natural", return_fields=True)
    outputs.append(_export("natural", nat, comm))
    _release(nat)

    abc8 = abc[PREMISE_PADDING]
    d_inf = _d_inf_pec()

    # ---- (P0) premise ---------------------------------------------------------
    delta = abs(pec["deficit_pp"] - P0_RUNG_PP)
    _say(comm, f"\n(P0) PEC 0.08 m deficit {pec['deficit_pp']:+.4f} pp vs {P0_RUNG_PP:+.4f} pp: "
               f"|delta| {delta:.4f} pp (record {RECORD['p0_delta_pp']:.4f}; "
               f"ASSERTED <= {P0_BAND_PP})")
    assert delta <= P0_BAND_PP

    # ---- (E) the term is in the matrix ----------------------------------------
    rel = abs(pec["z21"].real) / abs(pec["z21"])
    _say(comm, f"(E) PEC |Re Z21|/|Z21| = {rel:.3e} (record 0.000e+00; "
               f"ASSERTED <= {IDENTITY_TOLERANCE})")
    assert rel <= IDENTITY_TOLERANCE
    for i, w in enumerate(ABC_PADDINGS):
        z = abc[w]["z11"].real
        _say(comm, f"(E) ABC {w} m: Re Z11 = {z:+.6e} Ohm (record "
                   f"{RECORD['re_z11'][i]:+.6e}; ASSERTED > 0)")
        assert z > 0.0

    # ---- (S) sign of the box term, and the headline table ---------------------
    _say(comm, "\n(S) headline table  [r(W) is a record; (R) retired, weekly 2026-10-07; "
               "sphere model 0.5; never asserted]")
    _say(comm, f"  D_inf^PEC {d_inf:+.6f} pp (fit over the PORT-1 PEC rungs)")
    b_abc_all = []
    for i, w in enumerate(ABC_PADDINGS):
        d = abc[w]["deficit_pp"]
        b_pec = _b_pec(w)
        b_abc = d - d_inf
        b_abc_all.append(b_abc)
        r = abs(b_abc) / abs(b_pec)
        orig = abs(b_abc) <= abs(b_pec) / ORIGINAL_FALL_FACTOR
        _say(comm, f"  W {w} m: cells {abc[w]['ncells']} (record {RECORD['cells'][i]}), "
                   f"d^ABC {d:+.4f} pp (record {RECORD['d_abc'][i]:+.4f}), "
                   f"B^ABC {b_abc:+.4f} pp (record {RECORD['b_abc'][i]:+.4f}), "
                   f"B^PEC {b_pec:+.4f} pp (record {RECORD['b_pec_fit'][i]:+.4f}), "
                   f"r(W) = {r:.4f} (record {RECORD['r'][i]:.4f}); "
                   f"original bar |B^ABC| <= |B^PEC|/{ORIGINAL_FALL_FACTOR:g}: {orig} (printed only)")
    for w in ABC_PADDINGS:
        assert _b_pec(w) < 0.0
    for b in b_abc_all:
        assert b > 0.0
    _say(comm, "(S) ASSERTED: B^PEC < 0 < B^ABC at every padding")

    # ---- power per wall (printed) --------------------------------------------
    lhs = abc8["z11"].real * abc8["i1"] ** 2
    resid = abs(lhs - abc8["boundary_power"]) / abs(abc8["boundary_power"])
    _say(comm, "\n[power at 0.08 m] absorbing / pinned / reflecting (printed only)")
    _say(comm, f"  ABC: Re Z11*I1^2 = {lhs:.6e} W vs oint|n x E|^2/eta0 = "
               f"{abc8['boundary_power']:.6e} W, rel residual {resid:.3e} "
               f"(record {RECORD['power_abc']:.6e} W, residual 1.216e-11)")
    _say(comm, f"  PEC: oint|n x E|^2/eta0 = {pec['boundary_power']:.6e} W "
               f"(record {RECORD['power_pec']:.6e}), Re Z11 = {pec['z11'].real:+.6e} Ohm")
    _say(comm, f"  natural: oint|n x E|^2/eta0 = {nat['boundary_power']:.6e} W "
               f"(record {RECORD['power_nat']:.6e}), Re Z11 = {nat['z11'].real:+.6e} Ohm")

    # ---- negative control (c1): PEC / natural separation ----------------------
    sep = abs(b_abc_all[0] - _b_pec(PREMISE_PADDING))
    _say(comm, f"\n(c1) |B^ABC - B^PEC| at 0.08 m = {sep:.4f} pp (record {RECORD['c1_sep_pp']:.4f}; printed)")
    dn = abs(nat["deficit_pp"] - pec["deficit_pp"])
    _say(comm, f"(c1) |d^NAT - d^PEC| = {dn:.4f} pp (record {RECORD['dnat_minus_dpec_pp']:.4f}; "
               f"ASSERTED > P0_BAND_PP {P0_BAND_PP})")
    assert dn > P0_BAND_PP

    c3 = abs(abc8["deficit_pp"] - nat["deficit_pp"])
    _say(comm, f"(c3) |d^ABC - d^NAT| at 0.08 m = {c3:.6f} pp (record {RECORD['c3_pp']:.6f}; printed only)")

    # ---- ABC free-exponent fit (printed) --------------------------------------
    abc_def = tuple(abc[w]["deficit_pp"] for w in ABC_PADDINGS)
    d_abc, p_abc = _free_exponent_fit(ABC_PADDINGS, abc_def)
    if d_abc is None:
        _say(comm, f"[fit] free-exponent fit failed on the ABC ladder {abc_def} "
                   "(record: failed; printed only)")
    else:
        _say(comm, f"[fit] D_inf^ABC = {d_abc:+.6f} pp, p^ABC = {p_abc:.4f} "
                   "(record: failed; a success is printed and reported, never asserted)")

    # ---- (M) / (c2): the post-solve mix, from this run's PEC and natural solves
    assert pec["ncells"] == CELLS_008 and nat["ncells"] == CELLS_008, (
        pec["ncells"], nat["ncells"])
    b_pec_s = pec["deficit_pp"] - d_inf          # solve-referenced (2b's form)
    b_nat = nat["deficit_pp"] - d_inf
    r_nat = abs(b_nat) / abs(b_pec_s)
    b_mix = _mix_deficit(pec, nat, W_PEC, W_NAT) - d_inf
    b_swap = _mix_deficit(pec, nat, W_NAT, W_PEC) - d_inf
    b_half = _mix_deficit(pec, nat, 0.5, 0.5) - d_inf
    bar = abs(b_pec_s) / MIX_BAR_DIVISOR
    w_star = r_nat / (1.0 + r_nat)
    _say(comm, f"\n[M] cells on the PEC and natural 0.08 m meshes == {CELLS_008} "
               "(ASSERTED; 2b's own module was not re-run)")
    _say(comm, f"[M] B^PEC (solve-referenced) {b_pec_s:+.4f} pp (record {RECORD['b_pec_solve']:+.4f}); "
               f"fit-referenced {_b_pec(PREMISE_PADDING):+.4f} pp; B^NAT {b_nat:+.4f} pp "
               f"(record {RECORD['b_nat']:+.4f}); r_NAT {r_nat:.4f} (record {RECORD['r_nat']:.4f})")
    _say(comm, f"(M) B^mix {b_mix:+.4f} pp (record {RECORD['b_mix']:+.4f}) against bar "
               f"|B^PEC|/{MIX_BAR_DIVISOR:g} = {bar:.4f} pp (record {RECORD['bar']:.4f}); "
               f"ASSERTED |B^mix| <= bar")
    assert abs(b_mix) <= bar
    _say(comm, f"(c2) swapped weights: B^swap {b_swap:+.4f} pp (record {RECORD['b_swap']:+.4f}), "
               f"|B^swap|/bar = {abs(b_swap) / bar:.4f} (record {RECORD['swap_ratio']:.4f}); "
               f"ASSERTED |B^swap| > bar")
    assert abs(b_swap) > bar
    _say(comm, f"[M] equal weights: B^half {b_half:+.4f} pp (record {RECORD['b_half']:+.4f}); "
               f"w_PEC* = {w_star:.4f} (record {RECORD['w_star']:.4f}) (printed only)")

    # ---- outputs ---------------------------------------------------------------
    for path in outputs:
        _say(comm, f"[output] {path}")
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
