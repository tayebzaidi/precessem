"""Example (`EX-113`): ``Study.tune`` -- the mode capacitance is read from Im z_k only,
cancels mode 1 exactly, and the mode fields show which mode it tuned.

`WF-16` step 6(a) added ``Study.tune`` (``C = 1/(omega Im z_k)``, ``z_k = diag(F Z F^H)``).
No example called it.  This one runs the gate's own two Studies through the API (the
exact-rectangle PEC hole F-small birdcage and the solid sigma = 800 S/m control, degree 1,
64 MHz, ``-n 2``, kappa off), reproduces the gate's anchors, and adds three readings:

* a uniform series resistance r on all four ports shifts every ``z_k`` by exactly r and
  leaves ``C_1`` unchanged; a uniform reactance jX moves it (circuit identities, asserted);
* the series capacitor ``1/(j omega C_1)`` on every port cancels ``Im z_1`` and ``Im z_3``
  exactly and leaves mode 2 with ``Im z_2 - 1/(omega C_1) < 0`` (asserted / printed; a
  circuit-layer statement, **not** a resonance, f0, Q or tuning-accuracy claim);
* the mode fields: the incident-wave superposition of ``Study.b1_plus()``'s per-port CG1
  ``B`` with the DFT weights, printed at the coil centre and written to ParaView.

Every record, band and helper is imported from the gate module (the ANS-1 rule); record
digits sit only in the example-local ``REC`` dict, each with its ``log:line``, and are
never asserted by equality.  New values print at <= 6 significant digits.

Asserted: (W) rank width; (R0) cells + class means; (I) API = the `TH-17` modal route;
(a) ``C_1`` vs ``TH17_STEP3_C_MODE_RECORD_F`` within ``REDUCTION_BAND``; (ESR); (Tuned);
the mode-field port order.
Negative controls (asserted): (c1) k = 2 misses the record; (c2) the solid Study ties to
its own record and misses `TH-17`'s; (ESR-c) a uniform jX moves ``C_1``; (Tuned-c) the
``C_2`` capacitor does not cancel ``Im z_1``.

Run it through the example runner (the ``mri:`` group sources the complex build)::

    ./scripts/run_examples.sh -e mri:19 -n 2 -t 420

Not claimed: a resonance, f0, Q, S11 / Z_in or tuning-accuracy figure; copper tuning;
degree 2; 10 / 128 MHz; F-human; an external check.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.circuit import z_to_s  # noqa: E402
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.package.test_wf16_study_tune as g  # noqa: E402
from tests.validation import test_port_lumped_rlc_termination as rlc  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG  # noqa: E402
from tests.validation.test_th17_birdcage_eigenmodes import _four_port_dft  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "mri_19_study_tune_mode_capacitance"

# Example-local inputs (cited, not records): 0.2 ohm is the F-human ESR the gate's
# refusal test uses (tests/package/test_wf16_study_tune.py:108); 1.0 is a round number.
R_SERIES_OHM = (0.2, 1.0)
X_SERIES_OHM = 1.0

# Gate-log records, for printed comparison only (never asserted by equality).
REC = {
    # docs/testing/logs/20261014T004318Z_WF-16-s6-W1.log:944-947 (pec), :2760-2763 (solid)
    "pec_c1": 6.626796421e-11,
    "pec_c2": 9.873468900e-11,
    "solid_c1": 6.166035422e-11,
    # :2773, :2777, :2781
    "gap_a": -2.976552e-04,
    "gap_c1": +4.894874e-01,
    "gap_c2": -6.980694e-02,
}


def _modes(z4):
    zk, _ = g._modal_z(z4)
    return zk


def _c_of(zk_k, omega):
    return 1.0 / (omega * zk_k.imag) if zk_k.imag > 0.0 else float("nan")


def _fresh_tune(s_prime, k):
    """``tune`` on a fresh PEC Study carrying ``s_prime`` (the gate's W0 pattern: no mesh)."""
    st = g._py_copper("pec")
    st.result = SimpleNamespace(s_matrix=s_prime)
    return st.tune(k=k)


def _zline(zk):
    return ", ".join(f"{z.real:+.5e}{z.imag:+.5e}j" for z in zk)


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    # ---- structural asserts -------------------------------------------------
    assert comm.size == g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH, (
        f"(W) the degree-1 record is asserted at -n {g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}; "
        f"got -n {comm.size}"
    )
    if not np.issubdtype(default_scalar_type, np.complexfloating):
        raise RuntimeError("complex build required (source dolfinx-complex-mode)")
    assert g.REDUCTION_BAND is rlc.REDUCTION_BAND
    probe = g._py_copper("pec")
    assert probe.degree == 1 and probe.conductor is None
    z0 = g.REFERENCE_IMPEDANCE_OHM
    assert float(probe.excitation.z0_ohm) == z0
    omega = probe.omega
    ident = np.eye(4, dtype=np.complex128)

    # ---- (1) anchor: the gate's two Studies through the API -----------------
    out = {}
    for name, make in (("pec", lambda: g._py_copper("pec")), ("solid", g._solid_study)):
        study = make()
        assert study.degree == 1
        t0 = time.perf_counter()
        res = study.solve(comm)
        t_solve = time.perf_counter() - t0
        s = np.asarray(res.s_matrix, dtype=np.complex128)
        out[name] = {
            "study": study, "res": res, "s": s, "means": g._class_means(s),
            "t1": study.tune(k=1, comm=comm), "t2": study.tune(k=2, comm=comm),
            "t": t_solve,
        }
        if name == "pec":
            # EX-57 setup figure; no-op unless FEM_EM_SETUP_FIGURES=1.
            msh0 = res.mesh
            local_tags = np.unique(res.cell_tags.values)
            all_tags = sorted({int(t) for part in comm.allgather(local_tags.tolist())
                               for t in part})
            air_tags = tuple(t for t in all_tags if t != PHANTOM_CELL_TAG)
            names = {PHANTOM_CELL_TAG: "gel phantom"}
            for t in air_tags:
                names[t] = "air (outer box PEC)"
            write_setup_figure(
                msh0, res.cell_tags, FIGURE_DIR / f"{BASENAME}_setup.png",
                region_names=names, hide_tags=air_tags,
                translucent_tags=(PHANTOM_CELL_TAG,),
                slice_normal=(0.0, 0.0, 1.0), slice_origin=(0.0, 0.0, 0.0),
                title=("mri:19 -- Study.tune on the exact-rectangle PEC hole F-small "
                       "birdcage, four gap ports, phantom (WF-16 step 6a)"),
                comm=comm,
            )

    pec, solid = out["pec"], out["solid"]
    s_pec = pec["s"]
    t1, t2 = pec["t1"], pec["t2"]
    rel_means = g._rel_means(pec["means"], g.ANS6_EXACT_STRIP_CLASS_MEANS_P1[g.KEY]["pec"])
    zk_api = np.asarray(t1.z_k)
    z_pec = g.s_to_z(s_pec, z0)
    zk_ref, resid_ref = g._modal_z(z_pec)
    dz = float(np.max(np.abs(zk_api - zk_ref) / np.abs(zk_ref)))
    dr = abs(t1.off_diagonal_residual - resid_ref) / max(abs(resid_ref), np.finfo(float).tiny)
    c_closure = abs(t1.c_mode_farad * omega * zk_ref[1].imag - 1.0)
    gap_a = t1.c_mode_farad / g.TH17_STEP3_C_MODE_RECORD_F - 1.0
    gap_c1 = t2.c_mode_farad / g.TH17_STEP3_C_MODE_RECORD_F - 1.0
    tie_solid = g._rel(solid["s"], g._record(1))
    gap_c2 = solid["t1"].c_mode_farad / g.TH17_STEP3_C_MODE_RECORD_F - 1.0
    zk_solid = np.asarray(solid["t1"].z_k)

    # ---- example-side circuit anchors (numpy, no solve) ---------------------
    esr = {}
    for r in R_SERIES_OHM:
        s_r = z_to_s(z_pec + r * ident, z0)
        tr = _fresh_tune(s_r, 1)
        shift = float(np.max(np.abs(np.asarray(tr.z_k) - zk_api - r) / np.abs(zk_api)))
        esr[r] = (abs(tr.c_mode_farad / t1.c_mode_farad - 1.0), shift)
    s_x = z_to_s(z_pec + 1j * X_SERIES_OHM * ident, z0)
    c_x = _fresh_tune(s_x, 1).c_mode_farad
    esr_c = abs(c_x / t1.c_mode_farad - 1.0)

    z_t = z_pec + ident / (1j * omega * t1.c_mode_farad)
    zt = _modes(z_t)
    tuned_k1 = abs(zt[1].imag) / zk_ref[1].imag
    tuned_k3 = abs(zt[3].imag) / zk_ref[3].imag
    z_t2 = z_pec + ident / (1j * omega * t2.c_mode_farad)
    tuned_c = abs(zk_ref[1].imag - 1.0 / (omega * t2.c_mode_farad)) / zk_ref[1].imag

    # ---- the mode fields (b1_plus triggers a second, field-keeping sweep) ---
    t0 = time.perf_counter()
    b1 = pec["study"].b1_plus(comm=comm)
    t_fields = time.perf_counter() - t0
    fres = pec["study"].fields_result
    port_order = [d.port_id for d in fres.port_definitions]
    assert port_order == [d.port_id for d in pec["res"].port_definitions]
    assert port_order == sorted(b1, key=lambda p: int(p[1:])), (port_order, sorted(b1))
    s_f = np.asarray(fres.s_matrix, dtype=np.complex128)
    ds_fields = float(np.max(np.abs(s_f - s_pec)) / np.max(np.abs(s_pec)))

    f_dft = _four_port_dft()
    space = b1[port_order[0]].function_space
    msh = space.mesh
    mode_b = {}
    for k in (1, 2, 3):
        fk = fem.Function(space, name=f"B_mode{k}")
        acc = np.zeros_like(b1[port_order[0]].x.array)
        for p, pid in enumerate(port_order):
            acc = acc + np.conj(f_dft[k, p]) * b1[pid].x.array
        fk.x.array[:] = acc
        fk.x.scatter_forward()
        mode_b[k] = fk

    origin = np.zeros((1, 3))
    centre = {}
    for k in (1, 2, 3):
        vals, valid = evaluate_vector_field_parallel(mode_b[k], origin)
        vals = np.asarray(vals).reshape(-1, 3)
        assert bool(np.all(np.asarray(valid, dtype=bool))), f"mode {k}: origin not found"
        v = vals[0]
        centre[k] = {
            "b1p": abs(v[0] + 1j * v[1]) / 2.0,
            "b1m": abs(v[0] - 1j * v[1]) / 2.0,
            "norm": float(np.linalg.norm(v)),
        }

    # ---- outputs -------------------------------------------------------------
    space_s = fem.functionspace(msh, ("Lagrange", 1))

    def _vec_real(k):
        f = fem.Function(space, name=f"B_mode{k}_real")
        f.x.array[:] = np.real(mode_b[k].x.array)
        f.x.scatter_forward()
        return f

    def _scalar(name, values):
        f = fem.Function(space_s, name=name)
        assert values.shape[0] == f.x.array.shape[0], (values.shape, f.x.array.shape)
        f.x.array[:] = values
        f.x.scatter_forward()
        return f

    def _nodes(k):
        return mode_b[k].x.array.reshape(-1, 3)

    nb1, nb2, nb3 = _nodes(1), _nodes(2), _nodes(3)
    fields = {
        "B_mode1_real": _vec_real(1),
        "B_mode2_real": _vec_real(2),
        "B_mode1_abs": _scalar("B_mode1_abs", np.sqrt(np.sum(np.abs(nb1) ** 2, axis=1))),
        "B_mode2_abs": _scalar("B_mode2_abs", np.sqrt(np.sum(np.abs(nb2) ** 2, axis=1))),
        "B1_plus_mode1": _scalar("B1_plus_mode1", np.abs(nb1[:, 0] + 1j * nb1[:, 1]) / 2.0),
        "B1_minus_mode1": _scalar("B1_minus_mode1", np.abs(nb1[:, 0] - 1j * nb1[:, 1]) / 2.0),
        "B1_plus_mode3": _scalar("B1_plus_mode3", np.abs(nb3[:, 0] + 1j * nb3[:, 1]) / 2.0),
    }
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined", msh, fres.cell_tags, fields, comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    csv_path = OUTPUT_DIR / f"{BASENAME}_modes.csv"
    if comm.rank == 0:
        with open(csv_path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["study", "k", "variant", "Re_zk_ohm", "Im_zk_ohm", "C_k_F"])
            for sname in ("pec", "solid"):
                zs = g.s_to_z(out[sname]["s"], z0)
                c1_s, c2_s = out[sname]["t1"].c_mode_farad, out[sname]["t2"].c_mode_farad
                variants = {"base": zs}
                for r in R_SERIES_OHM:
                    variants[f"r={r:g}"] = zs + r * ident
                variants[f"X={X_SERIES_OHM:g}"] = zs + 1j * X_SERIES_OHM * ident
                variants["tuned_C1"] = zs + ident / (1j * omega * c1_s)
                variants["tuned_C2"] = zs + ident / (1j * omega * c2_s)
                for vname, zv in variants.items():
                    zkv = _modes(zv)
                    for k in range(4):
                        c = _c_of(zkv[k], omega)
                        w.writerow([sname, k, vname, f"{zkv[k].real:.6e}", f"{zkv[k].imag:.6e}",
                                    "" if np.isnan(c) else f"{c:.6e}"])

    # ---- prints (rank 0; no collective inside) --------------------------------
    if comm.rank == 0:
        print("=" * 72)
        print("EX-113 -- Study.tune: C from Im z_k, the C_1 closure, the mode fields")
        print("=" * 72)
        print(f"[study] degree 1, -n {comm.size}, f = {pec['study'].frequency_hz:.6e} Hz, "
              f"z0 = {z0:g} ohm, omega = {omega:.6e} rad/s")
        print(f"[(W)] comm.size {comm.size} == record width "
              f"{g.ANS6_EXACT_STRIP_RECORD_RANK_WIDTH}")
        print(f"[pec] cells {pec['res'].cells} (record {g.ANS6_EXACT_STRIP_HOLE_CELL_RECORD}), "
              f"mesh {pec['res'].mesh_time:.1f} s, sweep {pec['res'].sweep_time:.1f} s, "
              f"total {pec['t']:.1f} s")
        print(f"[solid] cells {solid['res'].cells}, mesh {solid['res'].mesh_time:.1f} s, "
              f"sweep {solid['res'].sweep_time:.1f} s, total {solid['t']:.1f} s")
        print("[(R0)] pec class-mean rel vs 'pec' record: "
              + "  ".join(f"{c} {rel_means[c]:.3e}" for c in g.CLASSES)
              + f" (band {g.LEG_D0_REPRODUCTION_BAND:.0e})")
        print(f"[(I)] max rel |z_k(api) - z_k(TH-17 route)| {dz:.3e}; residual rel {dr:.3e} "
              f"(api residual {t1.off_diagonal_residual:.3e}); |C_1 omega Im z_1 - 1| "
              f"{c_closure:.3e} (band {g.EXACT_IDENTITY_RTOL:g})")
        print("[z_k] pec   (Re, Im, ohm): " + _zline(zk_api))
        print("[z_k] solid (Re, Im, ohm): " + _zline(zk_solid))
        print(f"[z_k] Re z_1 solid/pec {zk_solid[1].real / zk_api[1].real:.5g} "
              f"(predicted 239); C_1 solid/pec "
              f"{solid['t1'].c_mode_farad / t1.c_mode_farad:.5g} (predicted 0.9305)")
        print(f"[C] pec C_1 = {t1.c_mode_farad:.9e} F, C_2 = {t2.c_mode_farad:.9e} F; "
              f"solid C_1 = {solid['t1'].c_mode_farad:.9e} F, C_2 = "
              f"{solid['t2'].c_mode_farad:.5e} F")
        print(f"[C] vs REC (printed): pec C_1 {t1.c_mode_farad / REC['pec_c1'] - 1:+.2e}, "
              f"pec C_2 {t2.c_mode_farad / REC['pec_c2'] - 1:+.2e}, "
              f"solid C_1 {solid['t1'].c_mode_farad / REC['solid_c1'] - 1:+.2e}")
        print(f"[(a)] C_1 vs TH17_STEP3_C_MODE_RECORD_F gap {gap_a:+.6e} (gate record "
              f"{REC['gap_a']:+.6e}); band {g.REDUCTION_BAND:.1e}")
        print(f"[(c1)] k = 2 gap {gap_c1:+.6e} (gate record {REC['gap_c1']:+.6e}); "
              f"{abs(gap_c1) / g.REDUCTION_BAND:.3g}x band")
        print(f"[(c2)] solid tie to its record {tie_solid:.3e} (band "
              f"{g.LEG_D0_REPRODUCTION_BAND:.0e}); gap {gap_c2:+.6e} (gate record "
              f"{REC['gap_c2']:+.6e}); {abs(gap_c2) / g.REDUCTION_BAND:.3g}x band")
        for r, (dc, shift) in esr.items():
            print(f"[(ESR)] r = {r:g} ohm: |C_1'/C_1 - 1| {dc:.3e}; max |z_k' - z_k - r|/|z_k| "
                  f"{shift:.3e} (band {g.EXACT_IDENTITY_RTOL:g})")
        print(f"[(ESR-c)] jX = {X_SERIES_OHM:g} ohm: |C_1'/C_1 - 1| {esr_c:.5g} "
              f"(predicted 0.02596 in [0.025, 0.027]; asserted only > band)")
        print(f"[(Tuned)] C_1 series capacitor: |Im zt_1|/Im z_1 {tuned_k1:.3e}, "
              f"|Im zt_3|/Im z_3 {tuned_k3:.3e} (band {g.EXACT_IDENTITY_RTOL:g}); "
              f"Re zt_1 - Re z_1 {abs(zt[1].real - zk_ref[1].real):.3e}")
        print(f"[(Tuned)] Im zt_2 = {zt[2].imag:+.5g} ohm (predicted -12.3); "
              f"Im zt_0 = {zt[0].imag:+.5g}")
        print(f"[(Tuned-c)] C_2 capacitor: |Im z_1 - 1/(omega C_2)|/Im z_1 {tuned_c:.5g} "
              f"(predicted 0.3288 in [0.32, 0.34]; asserted only > band); "
              f"Im zt'_1 = {_modes(z_t2)[1].imag:+.5g} ohm")
        print(f"[fields] port order {port_order}; field-keeping sweep {t_fields:.1f} s; "
              f"|S(fields) - S(plain)|max/|S|max {ds_fields:.3e} (predicted <= 1e-12)")
        for k in (1, 2, 3):
            c = centre[k]
            print(f"[mode {k}] centre |B1+| {c['b1p']:.5e}, |B1-| {c['b1m']:.5e}, "
                  f"|B| {c['norm']:.5e} T per unit incident-wave mode excitation")
        print(f"[mode] |B_2(0)|/|B_1(0)| {centre[2]['norm'] / centre[1]['norm']:.4e} "
              f"(predicted <~ 1e-2)")
        print(f"[paraview] wrote {xdmf_path}; {csv_path}")

    # ---- asserts -------------------------------------------------------------
    assert pec["res"].cells == g.ANS6_EXACT_STRIP_HOLE_CELL_RECORD
    for c in g.CLASSES:
        assert rel_means[c] <= g.LEG_D0_REPRODUCTION_BAND, f"(R0) {c} {rel_means[c]:.3e}"
    assert dz <= g.EXACT_IDENTITY_RTOL, f"(I) z_k {dz:.3e}"
    assert dr <= g.EXACT_IDENTITY_RTOL, f"(I) residual {dr:.3e}"
    assert t1.k == 1 and t1.frequency_hz == g.FREQUENCY_64_HZ
    assert zk_api[1].imag > 0.0
    assert c_closure <= g.EXACT_IDENTITY_RTOL, f"(I) C closure {c_closure:.3e}"
    assert abs(gap_a) <= g.REDUCTION_BAND, f"(a) {gap_a:+.6e}"
    assert abs(gap_c1) > g.REDUCTION_BAND, f"(c1) blind {gap_c1:+.6e}"
    assert tie_solid <= g.LEG_D0_REPRODUCTION_BAND, f"(c2) solid tie {tie_solid:.3e}"
    assert abs(gap_c2) > g.REDUCTION_BAND, f"(c2) blind {gap_c2:+.6e}"
    for r, (dc, shift) in esr.items():
        assert dc <= g.EXACT_IDENTITY_RTOL, f"(ESR) r={r} C shift {dc:.3e}"
        assert shift <= g.EXACT_IDENTITY_RTOL, f"(ESR) r={r} z_k shift {shift:.3e}"
    assert esr_c > g.EXACT_IDENTITY_RTOL, f"(ESR-c) blind {esr_c:.3e}"
    assert tuned_k1 <= g.EXACT_IDENTITY_RTOL, f"(Tuned) k=1 {tuned_k1:.3e}"
    assert tuned_k3 <= g.EXACT_IDENTITY_RTOL, f"(Tuned) k=3 {tuned_k3:.3e}"
    assert tuned_c > g.EXACT_IDENTITY_RTOL, f"(Tuned-c) blind {tuned_c:.3e}"

    if comm.rank == 0:
        print(
            f"\nAll assertions hold. Total elapsed {time.perf_counter() - started:.1f} s.",
            flush=True,
        )


if __name__ == "__main__":
    main()
