"""`IMP-1` step 3(b): the three rungs with the exact finite-k trace, degree 2
(PROJECT_PLAN §9 item 146; brief `docs/planning/chunks/IMP-1.md`).

Step 2b's fixture (`test_imp1_step2_refinement.py`) at a/6, a/10 and a/12, with
the box and cavity Dirichlet trace switched from the quasi-static field
(`_pec_exterior_numpy`) to the exact l = 1 TM field of step 3(a)
(`pec_sphere_l1_tm_field`, k = `K_FIXTURE`), and the comparand switched to the
same exact field integrated over the rung's discrete shell mask:

    Φ_ex,h = Σ_owned masked cells ∫ |E_reg + β_k E_out|² / Σ vol   (numpy, qdeg 8)

on the *unmasked* `pec_sphere_l1_tm_parts` — the same analytic continuation
`_phi_closed_form_on_mask` uses on chord points.  ρ^ex = Φ_h/Φ_ex,h − 1 (no β
removed: the exact field carries its own β_k).

Anchors (asserted): (L) the trace at s·k (s = 1e-6) equals
`_pec_exterior_numpy()` to 1e-9 on 1 000 points in the (3, N) layout;
(Mesh) the banked cell counts; (M′) the numpy quadrature of the quasi-static
field at β = 1 equals `_phi_closed_form_on_mask(…, 1.0)` to 1e-12; per rung
(i) |β_h − 1| ≤ `TH15_BETA_BAND`, (ii) |Φ_h/Φ_ex,h − 1| ≤ 10 %, (iii) the
cavity holds.  Item 146's two decision tests are replaced by the record
``test_rho_exact_recorded`` (step 3(c), §9 item 165).  Step 3(c) adds
``test_a6_degree3_exact_trace`` (a/6 at degree 3, both traces; reading printed,
no decision assert).  Negative control: the natural cavity at a/10 with the
exact trace, asserted only above the bands.
"""

from __future__ import annotations

import resource
import time

import basix
import numpy as np
import pytest
from mpi4py import MPI
from petsc4py import PETSc

from fem_em_solver.utils.analytical import (
    pec_sphere_l1_tm_beta,
    pec_sphere_l1_tm_field,
    pec_sphere_l1_tm_parts,
)
from tests.complex_mode import complex_only
from tests.validation.test_imp1_pec_sphere_lossy_host import (
    BETA_PEC,
    BOX_HALF_WIDTH,
    CAVITY_TAG,
    DEGREE,
    E0,
    MAT4_BAND,
    OUTER_BOUNDARY_TAG,
    SPHERE_RADIUS,
    TH15_BETA_BAND,
    _dipole_basis,
    _hole_mesh,
    _pec_exterior_numpy,
    _phi_closed_form_on_mask,
    _shell_mask,
    _solve,
)
from tests.validation.test_imp1_step2_refinement import (
    BANKED_RHO,
    RUNG_COARSE,
    RUNG_FINE,
    RUNG_OUT,
)
from tests.validation.test_imp1_step3_exact_comparand import K_FIXTURE

A = SPHERE_RADIUS
QDEG = 8
RHO_EX_BAND = 0.1e-2  # the chain's "under 0.1 %" (item 146 decision test, recorded red)
# Step 3(b)'s record (`20261001T213412Z_IMP-1-step3b-W1.log:528`), signed, %.
BANKED_RHO_EX = {"a/6": -0.9522, "a/10": -0.3400, "a/12": -0.2370}
# Step 3(b)'s exact-trace Φ_h as printed (`…213412Z_IMP-1-step3b-W1.log:374,443,522`).
W1_PHI_H = {"a/6": 1.592857, "a/10": 1.591119, "a/12": 1.592182}
PHI_REPRO_BAND = 1e-7  # record-reproduction band (item 165), not a gate

# Banked cell counts (`20261001T033100Z_IMP-1-land-W1.log:371,445,528`).
BANKED_CELLS = {"a/6": 13467, "a/10": 54524, "a/12": 90713}
# Banked quasi-static-trace Φ_h (`20261001T033100Z_IMP-1-land-W1.log:374,448,531`).
BANKED_PHI_H = {"a/6": 1.606456, "a/10": 1.604703, "a/12": 1.605776}
# Step 3(a)'s Re-E LSQ β on the exact field (item 146 brief), printed only.
BETA_LSQ_EXACT = 1.001933
# Item 146's predictions, printed beside the measurement, never asserted.
RHO_EX_PREDICTED = {"a/6": 0.35e-2, "a/10": 0.27e-2, "a/12": 0.38e-2}
RHO_QS_PREDICTED = {"a/6": 1.01e-2, "a/10": 0.30e-2, "a/12": 0.42e-2}

RUNGS = (("a/6", RUNG_COARSE), ("a/10", RUNG_FINE), ("a/12", RUNG_OUT))


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _exact_trace(k=K_FIXTURE):
    """The exact field as a dolfinx interpolation callable, (3, N) in and out;
    zero for r < a(1 − 1e-12) via `pec_sphere_l1_tm_field`'s mask."""

    def field(x: np.ndarray) -> np.ndarray:
        n = x.shape[1]
        if n == 0:
            return np.zeros((3, 0), dtype=PETSc.ScalarType)
        values = pec_sphere_l1_tm_field(np.ascontiguousarray(x.T), k=k, radius=A, e0=E0).T
        return np.ascontiguousarray(values, dtype=PETSc.ScalarType)

    return field


def _masked_quadrature(msh, chi):
    """Physical quadrature points and weights (degree 8, basix default scheme —
    UFL's for ``quadrature_degree: 8``) over the owned cells with χ = 1."""
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    geo_dofs = np.asarray(msh.geometry.dofmap)[:n_owned]
    if geo_dofs.shape[1] != 4:
        raise RuntimeError(f"non-affine geometry ({geo_dofs.shape[1]} nodes per cell)")
    chi_dofs = np.asarray(chi.function_space.dofmap.list)[:n_owned, 0]
    masked = np.real(chi.x.array[chi_dofs]) > 0.5
    coords = msh.geometry.x[geo_dofs[masked]]  # (m, 4, 3)
    ref_pts, ref_w = basix.make_quadrature(basix.CellType.tetrahedron, QDEG)
    v0 = coords[:, 0, :]
    jac = np.stack(
        [coords[:, 1, :] - v0, coords[:, 2, :] - v0, coords[:, 3, :] - v0], axis=2
    )  # (m, 3, 3), columns are the edge vectors
    det = np.abs(np.linalg.det(jac)) if len(coords) else np.zeros(0)
    pts = v0[:, None, :] + np.einsum("mij,qj->mqi", jac, ref_pts)
    weights = det[:, None] * ref_w[None, :]
    return pts.reshape(-1, 3), weights.reshape(-1)


def _phi_numpy(comm, pts, weights, field_values):
    num_local = float(np.sum(weights * np.sum(np.abs(field_values) ** 2, axis=1)))
    vol_local = float(np.sum(weights))
    num = comm.allreduce(num_local, op=MPI.SUM)
    vol = comm.allreduce(vol_local, op=MPI.SUM)
    return num / (E0**2 * vol)


def _exact_values(pts):
    if pts.shape[0] == 0:
        return np.zeros((0, 3), dtype=complex)
    beta_k = pec_sphere_l1_tm_beta(K_FIXTURE, A)
    e_reg, e_out = pec_sphere_l1_tm_parts(pts, k=K_FIXTURE, radius=A)
    return E0 * (e_reg + beta_k * e_out)


def _qs_values(pts):
    if pts.shape[0] == 0:
        return np.zeros((0, 3), dtype=complex)
    e = BETA_PEC * _dipole_basis(pts).astype(complex)
    e[:, 2] += E0
    return e


def _comparands(rung):
    """Φ_ex,h (exact, unmasked continuation) and the (M′) pair on the rung's mask."""
    msh = _hole_mesh(*rung)[0]
    comm = msh.comm
    chi = _shell_mask(msh)
    pts, w = _masked_quadrature(msh, chi)
    phi_ex = _phi_numpy(comm, pts, w, _exact_values(pts))
    phi_qs_numpy = _phi_numpy(comm, pts, w, _qs_values(pts))
    phi_qs_ufl = _phi_closed_form_on_mask(msh, chi, 1.0)
    return dict(msh=msh, chi=chi, phi_ex=phi_ex, phi_qs_numpy=phi_qs_numpy, phi_qs_ufl=phi_qs_ufl)


@pytest.fixture(scope="module")
def rungs():
    comm = MPI.COMM_WORLD
    pec = (OUTER_BOUNDARY_TAG, CAVITY_TAG)
    out = {}
    for label, rung in RUNGS:
        t0 = time.perf_counter()
        d = _solve(
            pec,
            DEGREE,
            resolution_sphere=rung[0],
            resolution_far=rung[1],
            dirichlet_e_field=_exact_trace(),
        )
        cmp_ = _comparands(rung)
        phi_cf_beta = _phi_closed_form_on_mask(cmp_["msh"], cmp_["chi"], d["beta"])
        seconds = time.perf_counter() - t0
        signed_rho_ex = d["phi_h"] / cmp_["phi_ex"] - 1.0
        signed_rho_qs = d["phi_h"] / phi_cf_beta - 1.0
        m_prime = abs(cmp_["phi_qs_numpy"] / cmp_["phi_qs_ufl"] - 1.0)
        out[label] = dict(
            d=d,
            signed_rho_ex=signed_rho_ex,
            rho_ex=abs(signed_rho_ex),
            signed_rho_qs=signed_rho_qs,
            m_prime=m_prime,
            **cmp_,
        )
        _say(
            f"\n[IMP-1 step 3b] {label}: {d['ncells']} cells (banked "
            f"{BANKED_CELLS[label]}), -n {comm.size}, solve + post {seconds:.1f} s"
        )
        _say(
            f"  {label}: rho^ex = Phi_h/Phi_ex,h - 1 = {signed_rho_ex:+.4%} "
            f"(Phi_h = {d['phi_h']:.6f}, Phi_ex,h = {cmp_['phi_ex']:.6f}; predicted "
            f"|rho^ex| ≈ {RHO_EX_PREDICTED[label]:.2%}, not asserted)"
        )
        _say(
            f"  {label}: Phi_h = {d['phi_h']:.6f} beside banked quasi-static-trace "
            f"Phi_h {BANKED_PHI_H[label]:.6f}: Delta_trace = "
            f"{d['phi_h'] / BANKED_PHI_H[label] - 1:+.4%}"
        )
        _say(
            f"  {label}: beta_h = {d['beta']:.6f} beside the exact field's LSQ beta "
            f"{BETA_LSQ_EXACT:.6f} ({d['beta'] / BETA_LSQ_EXACT - 1:+.4%}); "
            f"|beta_h - 1| = {abs(d['beta'] - 1):.4%} (band {TH15_BETA_BAND:.4%})"
        )
        _say(
            f"  {label} [record]: beta-removed quasi-static rho_h = Phi_h/Phi_cf,h(beta_h) "
            f"- 1 = {signed_rho_qs:+.4%} beside BANKED_RHO {BANKED_RHO[label]:.4%} "
            f"(predicted ≈ {RHO_QS_PREDICTED[label]:.2%}, not asserted)"
        )
        _say(
            f"  {label} [record]: (M') numpy Phi_qs(1) = {cmp_['phi_qs_numpy']:.12f} vs "
            f"UFL {cmp_['phi_qs_ufl']:.12f}: {m_prime:.2e} relative; Phi_ex,h/Phi_qs,h(1) "
            f"- 1 = {cmp_['phi_ex'] / cmp_['phi_qs_ufl'] - 1:+.4%}; cavity dofs = "
            f"{d['n_cavity']}, max |E| = {d['cavity_max']:.3e}; |Im E|/|Re E| = "
            f"{d['imag_ratio']:.4e}"
        )
    c, f, o = out["a/6"], out["a/10"], out["a/12"]
    _say(
        f"\n[IMP-1 step 3b] signed rho^ex a/6 -> a/10 -> a/12 = {c['signed_rho_ex']:+.4%} "
        f"-> {f['signed_rho_ex']:+.4%} -> {o['signed_rho_ex']:+.4%}; "
        f"rho^ex(a/12) - rho^ex(a/10) = {o['rho_ex'] - f['rho_ex']:+.4%}"
    )
    return out


@complex_only
@pytest.mark.integration
def test_three_rung_exact_trace(request):
    # (L): the trace at s·k equals the quasi-static trace, (3, N) layout.
    rng = np.random.default_rng(146)
    x = rng.uniform(-BOX_HALF_WIDTH, BOX_HALF_WIDTH, size=(3, 1000))
    lim = _exact_trace(1e-6 * K_FIXTURE)(x)
    qs = _pec_exterior_numpy()(x)
    l_defect = float(np.max(np.abs(lim - qs)))
    assert lim.shape == (3, 1000)
    empty = _exact_trace()(np.zeros((3, 0)))
    assert empty.shape == (3, 0)
    _say(
        f"\n[IMP-1 step 3b] (L) max |trace(1e-6 k) - quasi-static| on 1000 points = "
        f"{l_defect:.3e} (bound 1e-9)"
    )
    assert l_defect <= 1e-9, f"(L) trace(1e-6 k) differs from quasi-static by {l_defect:.3e}"

    out = request.getfixturevalue("rungs")
    for label, _ in RUNGS:
        r = out[label]
        d = r["d"]
        assert (
            d["ncells"] == BANKED_CELLS[label]
        ), f"(Mesh) {label}: {d['ncells']} cells != banked {BANKED_CELLS[label]}"
        assert r["m_prime"] <= 1e-12, (
            f"(M') {label}: numpy {r['phi_qs_numpy']:.12f} vs UFL {r['phi_qs_ufl']:.12f}: "
            f"{r['m_prime']:.2e}"
        )
        assert (
            abs(d["beta"] - BETA_PEC) <= TH15_BETA_BAND
        ), f"(i) {label}: beta = {d['beta']:.6f} outside {TH15_BETA_BAND:.4%}"
        assert (
            r["rho_ex"] <= MAT4_BAND
        ), f"(ii) {label}: |Phi_h/Phi_ex,h - 1| = {r['rho_ex']:.4%} > {MAT4_BAND:.0%}"
        assert d["n_cavity"] > 0, f"(iii) {label}: no cavity dofs"
        assert (
            d["cavity_max"] < 1e-12
        ), f"(iii) {label}: cavity dofs carry max |E| = {d['cavity_max']:.3e}"


@complex_only
@pytest.mark.integration
def test_rho_exact_recorded(rungs):
    """Step 3(b)'s record, re-read (§9 item 165, step 3(c)).

    Item 146 pre-registered two decision tests.  The level test
    (ρ^ex(a/12) ≤ 0.1 %) read 0.2370 % against 0.1 % — red, inside its
    pre-registered range [0.15 %, 0.65 %] — and was recorded as that rule's
    disposition required (`20261001T213412Z_IMP-1-step3b-W1.log:528,543`); the
    fall test read green on |ρ| only (signed ρ^ex rises a/10 → a/12).  Both are
    replaced by this record: signed ρ^ex is printed beside `BANKED_RHO_EX`, and the
    only assertion is that each rung's Φ_h reproduces W1's printed Φ_h — a
    record-reproduction band, not a gate.  W1 printed Φ_h to six decimals, so the
    1e-7 relative band is applied to Φ_h rounded to those six decimals (the band
    cannot be met by an unrounded value against a six-decimal print).
    """
    for label, _ in RUNGS:
        r = rungs[label]
        phi_h = r["d"]["phi_h"]
        _say(
            f"\n[IMP-1 step 3c] {label}: signed rho^ex = {r['signed_rho_ex']:+.4%} beside "
            f"BANKED_RHO_EX {BANKED_RHO_EX[label]:+.4f}%; Phi_h = {phi_h:.10f} beside "
            f"W1 Phi_h {W1_PHI_H[label]:.6f} (rel. {phi_h / W1_PHI_H[label] - 1:+.3e} "
            f"unrounded; {round(phi_h, 6) / W1_PHI_H[label] - 1:+.3e} at six decimals)"
        )
    for label, _ in RUNGS:
        phi_h = rungs[label]["d"]["phi_h"]
        rel = abs(round(phi_h, 6) / W1_PHI_H[label] - 1.0)
        assert rel <= PHI_REPRO_BAND, (
            f"{label}: Phi_h {phi_h:.10f} does not reproduce W1's {W1_PHI_H[label]:.6f} "
            f"({rel:.3e} > {PHI_REPRO_BAND:.0e})"
        )


def _rss_gib_summed(comm):
    """Process peak RSS (ru_maxrss, KiB on Linux) summed over ranks, GiB."""
    local = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return comm.allreduce(local, op=MPI.SUM) / 1024.0**2


@complex_only
@pytest.mark.integration
def test_a6_degree3_exact_trace():
    """a/6 at degree 3, exact and quasi-static traces (§9 item 165, step 3(c)).

    Discriminates geometric h² (faceting / chord sliver: ρ³ stays near the
    degree-2 −0.9522 %) from FE approximation error (ρ³ falls).  The reading
    rule is printed; there is no decision assert.  Anchors: (Mesh), (M′), and per
    trace (i)–(iii).
    """
    comm = MPI.COMM_WORLD
    pec = (OUTER_BOUNDARY_TAG, CAVITY_TAG)
    runs = {}
    for name, trace in (("exact", _exact_trace()), ("quasi-static", None)):
        t0 = time.perf_counter()
        d = _solve(
            pec,
            3,
            resolution_sphere=RUNG_COARSE[0],
            resolution_far=RUNG_COARSE[1],
            dirichlet_e_field=trace,
        )
        seconds = time.perf_counter() - t0
        v = d["fields"].e_complex.function_space
        ndofs = v.dofmap.index_map.size_global * v.dofmap.index_map_bs
        rss = _rss_gib_summed(comm)
        runs[name] = d
        _say(
            f"\n[IMP-1 step 3c] a/6 degree 3, {name} trace: {d['ncells']} cells, "
            f"{ndofs} dofs, solve + post {seconds:.1f} s, peak RSS summed over "
            f"{comm.size} ranks {rss:.2f} GiB (process peak so far)"
        )

    cmp_ = _comparands(RUNG_COARSE)
    m_prime = abs(cmp_["phi_qs_numpy"] / cmp_["phi_qs_ufl"] - 1.0)
    ex, qs = runs["exact"], runs["quasi-static"]
    rho3 = ex["phi_h"] / cmp_["phi_ex"] - 1.0
    delta_trace = ex["phi_h"] / qs["phi_h"] - 1.0
    if abs(rho3) >= 0.80e-2:
        branch = "GEOMETRIC (|rho3| >= 0.80 %)"
    elif abs(rho3) <= 0.30e-2:
        branch = "FE APPROXIMATION (|rho3| <= 0.30 %)"
    else:
        branch = "MIXED (0.30 % < |rho3| < 0.80 %)"
    _say(
        f"\n[IMP-1 step 3c] rho3 = Phi_h/Phi_ex,h - 1 (a/6, degree 3, exact trace) = "
        f"{rho3:+.4%} (Phi_h = {ex['phi_h']:.6f}, Phi_ex,h = {cmp_['phi_ex']:.6f}; "
        f"degree 2 {BANKED_RHO_EX['a/6']:+.4f}%)"
    )
    _say(f"[IMP-1 step 3c] pre-registered reading: {branch}")
    _say(
        f"  Delta_trace (degree 3) = Phi_h(exact)/Phi_h(quasi-static) - 1 = "
        f"{delta_trace:+.4%} (Phi_h qs = {qs['phi_h']:.6f}; predicted -0.8465% ± 0.01%)"
    )
    for name, d in runs.items():
        _say(
            f"  {name}: beta_h = {d['beta']:.6f} (degree 2 exact-trace 0.981501, exact "
            f"LSQ {BETA_LSQ_EXACT:.6f}); |beta_h - 1| = {abs(d['beta'] - 1):.4%} "
            f"(band {TH15_BETA_BAND:.4%}); |Phi_h/Phi_ex,h - 1| = "
            f"{abs(d['phi_h'] / cmp_['phi_ex'] - 1):.4%} (band {MAT4_BAND:.0%}); "
            f"cavity dofs {d['n_cavity']}, max |E| = {d['cavity_max']:.3e}; "
            f"|Im E|/|Re E| = {d['imag_ratio']:.4e}"
        )
    _say(
        f"  quasi-static record: Phi_h/Phi_cf,h(1) - 1 = "
        f"{qs['phi_h'] / qs['phi_cf_h'] - 1:+.4%}"
    )
    _say(
        f"  (M') numpy Phi_qs(1) = {cmp_['phi_qs_numpy']:.12f} vs UFL "
        f"{cmp_['phi_qs_ufl']:.12f}: {m_prime:.2e} relative"
    )

    assert ex["ncells"] == BANKED_CELLS["a/6"], f"(Mesh) {ex['ncells']} cells"
    assert m_prime <= 1e-12, f"(M') {m_prime:.2e}"
    for name, d in runs.items():
        assert (
            abs(d["beta"] - BETA_PEC) <= TH15_BETA_BAND
        ), f"(i) {name}: beta = {d['beta']:.6f} outside {TH15_BETA_BAND:.4%}"
        miss = abs(d["phi_h"] / cmp_["phi_ex"] - 1.0)
        assert miss <= MAT4_BAND, f"(ii) {name}: {miss:.4%} > {MAT4_BAND:.0%}"
        assert d["n_cavity"] > 0, f"(iii) {name}: no degree-3 cavity dofs located"
        assert d["cavity_max"] < 1e-12, f"(iii) {name}: cavity max |E| = {d['cavity_max']:.3e}"


@complex_only
@pytest.mark.integration
def test_control_exact_trace():
    """Natural cavity at a/10 with the exact trace (rule (e) negative control)."""
    t0 = time.perf_counter()
    ctl = _solve(
        (OUTER_BOUNDARY_TAG,),
        DEGREE,
        resolution_sphere=RUNG_FINE[0],
        resolution_far=RUNG_FINE[1],
        dirichlet_e_field=_exact_trace(),
    )
    cmp_ = _comparands(RUNG_FINE)
    seconds = time.perf_counter() - t0
    delta_beta = abs(ctl["beta"] - BETA_PEC)
    miss = abs(ctl["phi_h"] / cmp_["phi_ex"] - 1.0)
    _say(
        f"\n[IMP-1 step 3b] control, pec_facet_tags=(1,), exact trace, a/10 "
        f"({ctl['ncells']} cells, {seconds:.1f} s): beta = {ctl['beta']:.6f}, "
        f"|dbeta| = {delta_beta:.4f} (predicted ≈ 1.44; quasi-static record -0.438932)"
    )
    _say(
        f"  Phi_h = {ctl['phi_h']:.6f} vs Phi_ex,h = {cmp_['phi_ex']:.6f}: miss "
        f"{miss:.4%} (predicted 40-48 %; quasi-static record 43.94 %)"
    )
    assert (
        delta_beta > TH15_BETA_BAND
    ), f"natural cavity beta = {ctl['beta']:.6f} within the band of 1"
    assert miss > MAT4_BAND, f"natural cavity Phi miss {miss:.4%} not over {MAT4_BAND:.0%}"
