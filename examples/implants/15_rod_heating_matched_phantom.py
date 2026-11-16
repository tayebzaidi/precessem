"""Example (`EX-155`): the rod heats the coil-driven phantom -- like for like on a matched mesh pair.

`imp:3` shows the rod's SAR in the driven birdcage phantom; no example shows it
change the *temperature*. `IMP-3` step 11 (I39) did that with two separately
meshed phantoms (537 vs 51 693 cells); I41 repeated it on a **matched pair**:
arm M is `GEO-36`'s kept-cells twin of the rod build (rod cells tagged 403,
mapped to the gel, so the rod is just more phantom) and arm R is the PEC-hole
build. The two share every cell outside the rod. Both take the ccw quadrature
superposition at 64 MHz, one normalisation G = 1 W / P_phantom(arm M), and
`MAT-15`'s thermal route (Crank-Nicolson to t*, P1, Robin outside).

**It asserts, and it does not re-implement** (the `ANS-1` rule). Every helper,
record and band is imported by name from the gate modules (I41 as ``i41``).
I41's ``_run`` is not called: it always solves a third, vacuum-control arm
(403 unmapped, ~61 s) that this example leaves out, so the arm sequence is
written out here with the module's own helpers, unchanged.

Asserted: (A0m) arm M cells = ``KEPT_ROD_CELL_RECORD``, tag-403 cells = 15 000,
tag-3 cells equal to arm R's, all before any solve; (W0) arm M excludes 0
rod-wall facets, arm R exactly its parent tag-402 count; (B), (E), (E'), (S2)
on both arms at the imported bands; (C4) on arm M at ``C4_COVARIANCE_BAND``.

Printed only (records, never digit-asserted, item 537): the tip 1 g theta ratio
R / M, the 1 g SAR ratio Q, the max / mean theta ratios, P(R)/P(M) - 1, arm R's
(C4) pairs, cells and seconds.

Run it::

    ./scripts/run_examples.sh -e imp:15 -n 2 -t 420

Output: ``examples/implants/paraview_output/implants_15_rod_heating_matched_phantom_armM_combined.xdmf``
and ``..._armR_combined.xdmf`` (cell tags + the DG0 SAR source per arm; the
thermal route returns scalars, so theta is not exported).
"""

from __future__ import annotations

import gc
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from dolfinx import default_scalar_type, fem
from mpi4py import MPI

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.validation.test_imp3_rod_birdcage as _rb  # noqa: E402
import tests.validation.test_imp3_rod_heating_matched as i41  # noqa: E402
from tests.mesh.test_birdcage_implant_rod import (  # noqa: E402
    KEPT_ROD_CELL_RECORD,
    _build_kept_rod,
    _build_rod,
)
from tests.validation.test_imp3_rod_birdcage import TIPS  # noqa: E402
from tests.validation.test_imp3_rod_heating_phantom import (  # noqa: E402
    SOURCE,
    _heat_arm,
    _parent_wall_count,
)
from tests.validation.test_mat15_coil_driven_heating_phantom import (  # noqa: E402
    _cell_average_dg0_tag,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "implants_15_rod_heating_matched_phantom"
PHANTOM_TAG = i41.PHANTOM_CELL_TAG
ROD_TAG = i41.BIRDCAGE_IMPLANT_CELL_TAG
# I41 W1b record (20261025T201613Z_IMP-3-I41-W1b.log:3776-3777): printed beside, never asserted.
I41_RECORD_TH = {"top": 1.104515, "bottom": 1.108684}
I41_RECORD_Q = {"top": 4.090342, "bottom": 4.113695}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _export_arm(tag, msh, cell_tags, e_ccw, sigma_field, G):
    """Cell tags + the DG0 SAR source (W/m^3, as the heating route builds it, with G)."""
    comm = MPI.COMM_WORLD
    W, q_par = _cell_average_dg0_tag(
        msh, cell_tags, PHANTOM_TAG, G * 0.5 * sigma_field * ufl.inner(e_ccw, e_ccw))
    q = fem.Function(W, name="sar_source_w_per_m3")
    q.x.array[:] = np.real(q_par)
    q.x.scatter_forward()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = OUTPUT_DIR / f"{BASENAME}_arm{tag}_combined"
    write_xdmf_with_tags(path, msh, cell_tags, {"sar_source_w_per_m3": q}, comm=comm)
    return path


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "this example needs the complex DolfinX build "
            "(source /usr/local/bin/dolfinx-complex-mode); the runner's `imp:` "
            "group sources it automatically"
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-155 -- the rod heats the coil-driven phantom (IMP-3 step 11 / I41, matched pair)")
    _say(comm, "=" * 78)
    _say(comm, f"[gate] bands imported: B {i41.B_BAND:.0e}, E {i41.E_BAND:.0e}, S2 {i41.S2_BAND:.0e}, "
               f"C4 {i41.C4_COVARIANCE_BAND:.1%}; source {SOURCE}; tips {list(TIPS)}")

    heat, files = {}, {}
    built_rod = _build_rod()
    c_rod = i41._census(built_rod)
    rod_wall = _parent_wall_count(built_rod[0], built_rod[2])
    built_m = _build_kept_rod()
    c_m = i41._census(built_m)
    m_wall = _parent_wall_count(built_m[0], built_m[2])
    _say(comm, f"\n(A0m) census: arm M {c_m['cells']} cells (record {KEPT_ROD_CELL_RECORD}), "
               f"tag-3 {c_m['t3']}, tag-403 {c_m['t403']} (expected {i41.ROD_CELLS}); arm R "
               f"{c_rod['cells']} cells (record {_rb.ROD_CELL_RECORD_W1}), tag-3 {c_rod['t3']}; "
               f"parent tag-402 owned facets M {m_wall} (interior) / R {rod_wall} (exterior)")
    assert c_m["cells"] == KEPT_ROD_CELL_RECORD, c_m
    assert c_m["t403"] == i41.ROD_CELLS, c_m
    assert c_m["t3"] == c_rod["t3"], (c_m, c_rod)

    msh_m, ct_m = built_m[0], built_m[1]
    port_names = {100 + i: f"port P{i} box (lower)" for i in range(1, 5)}
    port_names.update({200 + i: f"port P{i} box (upper)" for i in range(1, 5)})
    write_setup_figure(
        msh_m,
        ct_m,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "conductor", 2: "air", PHANTOM_TAG: "phantom (saline)",
                      ROD_TAG: "implant rod cells (kept, mapped to gel in arm M)", **port_names},
        hide_tags=(2,),
        translucent_tags=(PHANTOM_TAG,),
        slice_normal=(0.0, 1.0, 0.0),
        slice_origin=(0.0, 0.0, 0.0),
        title="imp:15 -- matched pair: rod cells kept (arm M) vs PEC hole (arm R)",
        window_size=(1400, 630),
        comm=comm,
    )

    def hook_m(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        merged = i41._merge_rod_into_phantom(msh, cell_tags)
        heat["M"] = _heat_arm(label, rod, msh, merged, e_of[SOURCE], sigma_field, None, m_wall)
        files["M"] = _export_arm("M", msh, cell_tags, e_of[SOURCE], sigma_field, heat["M"]["G"])

    def hook_r(label, rod, msh, cell_tags, e_of, sigma_field, rho_field):
        heat["R"] = _heat_arm(label, rod, msh, cell_tags, e_of[SOURCE], sigma_field,
                              heat["M"]["G"], rod_wall)
        files["R"] = _export_arm("R", msh, cell_tags, e_of[SOURCE], sigma_field, heat["R"]["G"])

    m = i41._arm_with_403(built_m, "EX-155 arm M (rod as gel)", hook_m, map_gel=True)
    del built_m, msh_m, ct_m
    gc.collect()
    r = _rb._arm(built_rod, "EX-155 arm R (with rod)", rod=True, field_hook=hook_r)
    del built_rod
    gc.collect()

    hm, hr = heat["M"], heat["R"]
    q_sar = {t: r["sar"][(t, SOURCE)] / m["sar"][(t, SOURCE)] for t in TIPS}
    q_th = {t: hr["tip_th"][t] / hm["tip_th"][t] for t in TIPS}
    _say(comm, f"\n[material] tag 403 in arm M: had {m['material_403']['had_403']} -> "
               f"mapped {m['material_403']['mapped']}")
    for t in TIPS:
        _say(comm, f"[record] {t} tip, {SOURCE}: 1 g theta(t*) R/M = {q_th[t]:.6f} (I41 record "
                   f"{I41_RECORD_TH[t]:.6f}; I39 mesh pair {i41.I39_ROD_TH_RATIO[t]:.4f}); 1 g SAR "
                   f"ratio Q R/M = {q_sar[t]:.6f} (I41 record {I41_RECORD_Q[t]:.6f}); ball masses "
                   f"M {m['mass'][t] * 1e3:.6f} g / R {r['mass'][t] * 1e3:.6f} g")
    _say(comm, f"[record] max theta R/M = {hr['th_max'] / hm['th_max']:.6f} ({hr['th_max']:.6e} / "
               f"{hm['th_max']:.6e} K); mean theta R/M {hr['mean_th'] / hm['mean_th']:.6f}; "
               f"P_phantom M {hm['P']:.12e} W, R {hr['P']:.12e} W: P(R)/P(M) - 1 = "
               f"{hr['P'] / hm['P'] - 1.0:+.4%}")
    _say(comm, "[record] arm R (C4) pairs % (the rod breaks the C4 symmetry; printed, not asserted): "
               + ", ".join(f"{hr['pairs'][k] * 100:.4f}" for k in range(4)))

    # ---- assertions ---------------------------------------------------------------------------
    _say(comm, "\n[asserted]")
    _say(comm, f"  (W0) rod-wall facets excluded: arm M {hm['counts']['wall']} (must be 0), arm R "
               f"{hr['counts']['wall']} (must equal parent tag-402 count {hr['parent_wall']})")
    assert m["material_403"]["mapped"]
    assert hm["counts"]["wall"] == 0, hm["counts"]
    assert hr["parent_wall"] > 0 and hr["counts"]["wall"] == hr["parent_wall"], hr["counts"]
    for arm, h in (("M", hm), ("R", hr)):
        _say(comm, f"  arm {arm}: (B) {h['B']:.3e} <= {i41.B_BAND:.0e}; (E) {h['E']:.3e} and "
                   f"(E') {h['Ep']:.3e} <= {i41.E_BAND:.0e}; (S2) {h['S2']:.4e} <= {i41.S2_BAND:.0e}")
        assert h["B"] <= i41.B_BAND, f"(B) {arm}: {h['B']:.3e}"
        assert h["E"] <= i41.E_BAND, f"(E) {arm}: {h['E']:.3e}"
        assert h["Ep"] <= i41.E_BAND, f"(E') {arm}: {h['Ep']:.3e}"
        assert h["S2"] <= i41.S2_BAND, f"(S2) {arm}: {h['S2']:.4e}"
    worst = max(hm["pairs"][k] for k in range(4))
    _say(comm, f"  (C4) arm M worst pair {worst * 100:.4f}% <= {i41.C4_COVARIANCE_BAND:.1%}")
    for k in range(4):
        assert hm["pairs"][k] <= i41.C4_COVARIANCE_BAND, f"(C4) arm M pair {k}: {hm['pairs'][k] * 100:.4f}%"

    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    for arm in ("M", "R"):
        _say(comm, f"[output] {files[arm]}.xdmf -- cell tags + DG0 SAR source (theta not exported)")
    _say(comm, f"[timing] arms M {m['t_arm']:.1f} / R {r['t_arm']:.1f} s; cells M {m['cells']} / "
               f"R {r['cells']}; total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
