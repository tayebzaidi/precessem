"""`IMP-1` step 2b: the β-removed shell residual of the PEC-sphere implant at
three rungs in one process, degree 2 (PROJECT_PLAN §9 item 106; step 2's
first attempt, item 102, is banked below).

Step 1's fixture and anchors (`test_imp1_pec_sphere_lossy_host.py`) at
(a/6, a/3), (a/10, a/5) and the out-of-sample (a/12, a/6).  On the discrete
mask the closed form is ``Φ_cf,h(β) = 1 + 2βX_h + β²Y_h`` — the cross term
``X_h`` vanishes only on the exact shell — so the residual with the rung's own
fitted β removed is ``ρ_h = |Φ_h/Φ_cf,h(β_h) − 1|``.

At all three rungs: (i) ``|β − 1|`` within `TH-15`'s band (imported);
(ii) ``|Φ_h/Φ_cf,h(1) − 1|`` within `MAT-4`'s 10 %; (iii) the cavity dofs are
zero; (M) this module's ``Φ_cf,h(1)`` equals `_solve`'s to 1e-12.  Between
rungs: (iv′) ``|β − 1|`` at a/12 ≤ at a/6.

Banked 2026-09-30 (item 102, `20260930T111403Z_IMP-1-step2-W1.log:449,456`):
the old (v), total Φ miss at a/10 ≤ at a/6, read 0.0743 % → 0.6915 % — the
a/6 figure is a cancellation of the β term against the residual.  It is
printed, not asserted.

Banked 2026-09-30 18:00 review (item 106,
`20260930T170159Z_IMP-1-step2b-W1.log:539`): (v′), ``ρ`` at a/12 ≤ at a/10,
read 1.1650 % → 0.4473 % → 0.5694 % — the fixture's measured floor near
0.5 %, withdrawn as a gate (bands-only closure ratified by the 2026-09-30
weekly).  ``ρ`` is printed beside ``BANKED_RHO`` at all three rungs, not
asserted; no convergence claim is made.

Negative control (`test_control`): the natural cavity at a/10,
``pec_facet_tags=(1,)`` — asserted only to miss both bands.
"""

from __future__ import annotations

import math
import resource
import time

import pytest
from dolfinx import fem
from mpi4py import MPI

from tests.complex_mode import complex_only
from tests.validation.test_imp1_pec_sphere_lossy_host import (
    BETA_PEC,
    BETA_VOID,
    CAVITY_TAG,
    DEGREE,
    MAT4_BAND,
    OUTER_BOUNDARY_TAG,
    PHI_CLOSED_FORM_PEC,
    PHI_CLOSED_FORM_VOID,
    RESOLUTION_FAR,
    RESOLUTION_SPHERE,
    SPHERE_RADIUS,
    TH15_BETA_BAND,
    _hole_mesh,
    _phi_closed_form_on_mask,
    _sar_ratio_records,
    _shell_mask,
    _solve,
)

RUNG_COARSE = (RESOLUTION_SPHERE, RESOLUTION_FAR)  # a/6, a/3 — step 1's
RUNG_FINE = (SPHERE_RADIUS / 10.0, SPHERE_RADIUS / 5.0)
RUNG_OUT = (SPHERE_RADIUS / 12.0, SPHERE_RADIUS / 6.0)  # out of sample (item 106)
CELL_GUARD = 105_000  # item 106: a cost guard, not a gate — over this, do not solve

# Step 1's record (20260930T033210Z_IMP-1-step1-W1.log:372–373), printed beside
# the in-window a/6 values as a reference only — never asserted (cross-run drift).
STEP1_BETA_MISS = 1.4404e-2
STEP1_PHI_MISS = 0.0743e-2

# Predicted before the window (item 106), printed beside the measurement, never
# asserted: ρ at a/6 and a/10 from `20260930T111403Z_IMP-1-step2-W1.log:372–373,
# 443–444` with the mask's cross term taken as 0; a/12 from their two-point rate.
RHO_PREDICTED = {"a/6": 1.17e-2, "a/10": 0.448e-2, "a/12": 0.32e-2}
# Item 102's banked (v): total Φ miss a/6 -> a/10 (…111403Z…:449,456).
BANKED_PHI_MISS = (0.0743e-2, 0.6915e-2)
# Item 106's (v'), banked 2026-09-30 18:00 review as the fixture's measured floor
# record (`20260930T170159Z_IMP-1-step2b-W1.log:539`): ρ grew a/10 -> a/12; the
# model error ((|k|a)² ≈ 0.8 %) and rung-to-rung scatter are not separable at
# three rungs.  Printed beside the in-window ρ, never asserted (item 123).
BANKED_RHO = {"a/6": 1.1650e-2, "a/10": 0.4473e-2, "a/12": 0.5694e-2}


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _peak_rss_gib():
    """Summed process peak RSS over ranks [GiB] (``ru_maxrss`` is kB on Linux)."""
    local = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return MPI.COMM_WORLD.allreduce(local, op=MPI.SUM) / 1024.0**2


def _announce(label, rung):
    """Cells and degree-2 dofs of a rung, printed before its solve."""
    comm = MPI.COMM_WORLD
    msh, _, _ = _hole_mesh(*rung)
    ncells = msh.topology.index_map(msh.topology.dim).size_global
    v = fem.functionspace(msh, ("N1curl", DEGREE))
    ndofs = v.dofmap.index_map.size_global * v.dofmap.index_map_bs
    _say(
        f"\n[IMP-1 step 2] {label}: resolution_sphere = {rung[0]:.6e} m, "
        f"resolution_far = {rung[1]:.6e} m, {ncells} cells, {ndofs} degree-{DEGREE} "
        f"dofs, -n {comm.size} (before the solve)"
    )
    return ncells, ndofs


def _timed_solve(label, pec_facet_tags, rung):
    t0 = time.perf_counter()
    out = _solve(pec_facet_tags, DEGREE, resolution_sphere=rung[0], resolution_far=rung[1])
    out["seconds"] = time.perf_counter() - t0
    out["rss_gib"] = _peak_rss_gib()
    _say(
        f"[IMP-1 step 2] {label}: solve + post {out['seconds']:.1f} s, summed peak RSS "
        f"{out['rss_gib']:.2f} GiB (process peak so far)"
    )
    return out


def _print_rung(label, d):
    beta_miss = abs(d["beta"] - BETA_PEC)
    phi_miss = abs(d["phi_h"] / d["phi_cf_h"] - 1.0)
    _say(
        f"  {label}: beta = {d['beta']:.6f}  |beta - 1| = {beta_miss:.4%} "
        f"(band {TH15_BETA_BAND:.4%})"
    )
    _say(
        f"  {label}: Phi_h = {d['phi_h']:.6f}, Phi_cf,h = {d['phi_cf_h']:.6f}, "
        f"Phi_h/Phi_cf,h - 1 = {d['phi_h'] / d['phi_cf_h'] - 1.0:+.4%} (band {MAT4_BAND:.0%})"
    )
    _say(
        f"  {label} [record]: Phi_cf,h vs continuous closed form {PHI_CLOSED_FORM_PEC:.6f}: "
        f"{d['phi_cf_h'] / PHI_CLOSED_FORM_PEC - 1:+.4%}; |Im E|/|Re E| = "
        f"{d['imag_ratio']:.4e}; cavity dofs = {d['n_cavity']}, max |E| = "
        f"{d['cavity_max']:.3e}; {d['ncells']} cells"
    )
    for theta_deg, ratio, closed in _sar_ratio_records(d["fields"]):
        _say(
            f"  {label} [record]: SAR/SAR0 at r = 1.2a, polar {theta_deg:.0f} deg: "
            f"{ratio:.4f} (quasi-static closed form {closed:.4f})"
        )
    return beta_miss, phi_miss


def _residual(label, rung, d):
    """ρ_h and its companions for one rung: the closed form on the rung's own
    mask at β = 1 and at the rung's fitted β (one JIT compile per rung)."""
    msh = _hole_mesh(*rung)[0]  # cache hit
    chi = _shell_mask(msh)
    phi_cf_1 = _phi_closed_form_on_mask(msh, chi, BETA_PEC)
    phi_cf_beta = _phi_closed_form_on_mask(msh, chi, d["beta"])
    no_cross = 1.0 + d["beta"] ** 2 * (phi_cf_1 - 1.0)
    signed_total = d["phi_h"] / phi_cf_1 - 1.0
    signed_rho = d["phi_h"] / phi_cf_beta - 1.0
    m_defect = abs(phi_cf_1 / d["phi_cf_h"] - 1.0)
    _say(
        f"  {label}: Phi_h/Phi_cf,h(1) - 1 = {signed_total:+.4%}; "
        f"Phi_h/Phi_cf,h(beta_h) - 1 = {signed_rho:+.4%} "
        f"(rho predicted ≈ {RHO_PREDICTED[label]:.3%}, not asserted)"
    )
    _say(
        f"  {label} [record]: Phi_cf,h(beta_h) = {phi_cf_beta:.6f} beside "
        f"1 + beta_h^2 (Phi_cf,h(1) - 1) = {no_cross:.6f}: the mask's cross term "
        f"{phi_cf_beta - no_cross:+.3e}"
    )
    _say(
        f"  {label} [record]: (M) Phi_cf,h(1) here {phi_cf_1:.12f} vs _solve's "
        f"{d['phi_cf_h']:.12f}: {m_defect:.2e} relative"
    )
    return dict(
        phi_cf_1=phi_cf_1,
        signed_total=signed_total,
        signed_rho=signed_rho,
        rho=abs(signed_rho),
        m_defect=m_defect,
    )


def _rate(rho_a, n_a, rho_b, n_b):
    """Two-point rate of ρ between a/n_a and a/n_b (nan if either is zero)."""
    if rho_a <= 0.0 or rho_b <= 0.0:
        return float("nan")
    return math.log(rho_a / rho_b) / math.log(n_b / n_a)


@complex_only
@pytest.mark.integration
def test_three_rung_residual():
    pec = (OUTER_BOUNDARY_TAG, CAVITY_TAG)
    rungs = (("a/6", RUNG_COARSE, 6.0), ("a/10", RUNG_FINE, 10.0), ("a/12", RUNG_OUT, 12.0))
    out = {}

    for label, rung, _ in rungs:
        ncells, _ = _announce(label, rung)
        if label == "a/12":
            # The only assert ahead of the a/12 solve: the cost guard.
            assert (
                ncells <= CELL_GUARD
            ), f"a/12 mesh has {ncells} cells > {CELL_GUARD}: not solved (item 106)"
        d = _timed_solve(label, pec, rung)
        beta_miss, phi_miss = _print_rung(label, d)
        out[label] = dict(d=d, beta_miss=beta_miss, phi_miss=phi_miss, **_residual(label, rung, d))
        if label == "a/6":
            _say(
                f"  a/6 [reference, step 1's log, not asserted]: |beta - 1| = "
                f"{STEP1_BETA_MISS:.4%}, Phi miss = {STEP1_PHI_MISS:.4%}"
            )

    c, f, o = out["a/6"], out["a/10"], out["a/12"]
    _say(
        f"\n[IMP-1 step 2b] rho a/6 -> a/10 -> a/12 = {c['rho']:.4%} -> {f['rho']:.4%} "
        f"-> {o['rho']:.4%} (signed {c['signed_rho']:+.4%}, {f['signed_rho']:+.4%}, "
        f"{o['signed_rho']:+.4%}; predicted ≈ 1.17 % -> 0.448 % -> 0.32 %)"
    )
    for (la, ra, na), (lb, rb, nb) in (
        (("a/6", c, 6.0), ("a/10", f, 10.0)),
        (("a/10", f, 10.0), ("a/12", o, 12.0)),
        (("a/6", c, 6.0), ("a/12", o, 12.0)),
    ):
        _say(
            f"  [record] pairwise rate of rho, {la} -> {lb}: "
            f"{_rate(ra['rho'], na, rb['rho'], nb):.3f} "
            f"(the a/6 -> a/10 prediction was 1.87; not asserted)"
        )
    _say(
        f"  [record] |beta - 1| a/6 -> a/10 -> a/12 = {c['beta_miss']:.4%} -> "
        f"{f['beta_miss']:.4%} -> {o['beta_miss']:.4%} (signed beta - 1: "
        f"{c['d']['beta'] - 1:+.4%}, {f['d']['beta'] - 1:+.4%}, "
        f"{o['d']['beta'] - 1:+.4%}); a/10 -> a/12 is not asserted"
    )
    _say(
        f"  [banked 2026-09-30, not asserted] old (v), total Phi miss a/6 -> a/10: "
        f"{c['phi_miss']:.4%} -> {f['phi_miss']:.4%} (record "
        f"{BANKED_PHI_MISS[0]:.4%} -> {BANKED_PHI_MISS[1]:.4%}); a/12: "
        f"{o['phi_miss']:.4%}"
    )

    for label, _, _ in rungs:
        r = out[label]
        d = r["d"]
        assert d["n_cavity"] > 0, f"(iii) {label}: no dofs located on cavity facet tag 2"
        assert (
            d["cavity_max"] < 1e-12
        ), f"(iii) {label}: cavity dofs carry max |E| = {d['cavity_max']:.3e}"
        assert r["m_defect"] <= 1e-12, (
            f"(M) {label}: Phi_cf,h(1) = {r['phi_cf_1']:.12f} differs from _solve's "
            f"{d['phi_cf_h']:.12f} by {r['m_defect']:.2e} relative"
        )
        assert r["beta_miss"] <= TH15_BETA_BAND, (
            f"(i) {label}: beta = {d['beta']:.6f} misses 1 by {r['beta_miss']:.4%} > "
            f"{TH15_BETA_BAND:.4%}"
        )
        assert (
            r["phi_miss"] <= MAT4_BAND
        ), f"(ii) {label}: Phi miss {r['phi_miss']:.4%} > {MAT4_BAND:.0%}"
    assert o["beta_miss"] <= c["beta_miss"], (
        f"(iv') |beta - 1| at a/12 above a/6: {c['beta_miss']:.4%} (a/6) -> "
        f"{o['beta_miss']:.4%} (a/12)"
    )
    _say(
        f"  [banked 2026-09-30, not asserted] (v'), rho a/6 -> a/10 -> a/12: "
        f"{c['rho']:.4%} -> {f['rho']:.4%} -> {o['rho']:.4%} (record "
        f"{BANKED_RHO['a/6']:.4%} -> {BANKED_RHO['a/10']:.4%} -> "
        f"{BANKED_RHO['a/12']:.4%}; the fixture's floor near 0.5 %)"
    )


@complex_only
@pytest.mark.integration
def test_control():
    """Natural cavity at a/10: drop tag 2 from the Dirichlet set."""
    ncells, _ = _announce("a/10 control", RUNG_FINE)
    assert ncells <= CELL_GUARD, f"a/10 mesh has {ncells} cells > {CELL_GUARD}"
    ctl = _timed_solve("a/10 control", (OUTER_BOUNDARY_TAG,), RUNG_FINE)
    delta_beta = abs(ctl["beta"] - BETA_PEC)
    phi_miss = abs(ctl["phi_h"] / ctl["phi_cf_h"] - 1.0)
    _say(
        f"[IMP-1 step 2] negative control, pec_facet_tags=(1,), degree {DEGREE}, a/10: "
        f"beta = {ctl['beta']:.6f} (void {BETA_VOID:+.2f}; predicted ≈ -0.43), |dbeta| = "
        f"{delta_beta:.4f}"
    )
    _say(
        f"  Phi_h = {ctl['phi_h']:.6f} (void closed form {PHI_CLOSED_FORM_VOID:.4f}) vs "
        f"Phi_cf,h(PEC) = {ctl['phi_cf_h']:.6f}: miss {phi_miss:.4%} (predicted ≈ 44 %)"
    )
    assert (
        delta_beta > TH15_BETA_BAND
    ), f"natural cavity beta = {ctl['beta']:.6f} within the band of 1"
    assert phi_miss > MAT4_BAND, f"natural cavity Phi miss {phi_miss:.4%} not over {MAT4_BAND:.0%}"
