"""Unit tests for coil+phantom region mesh resolution policy helpers."""

import pytest

from fem_em_solver.io.mesh import MeshGenerator


def test_region_resolution_policy_defaults_to_global_resolution():
    policy = MeshGenerator.coil_phantom_region_resolution_policy(resolution=0.015)

    assert policy["coil_resolution_m"] == pytest.approx(0.015)
    assert policy["phantom_resolution_m"] == pytest.approx(0.015)
    assert policy["air_resolution_m"] == pytest.approx(0.015)
    assert policy["min_resolution_m"] == pytest.approx(0.015)
    assert policy["max_resolution_m"] == pytest.approx(0.015)


def test_region_resolution_policy_prefers_explicit_kwargs_over_mapping():
    policy = MeshGenerator.coil_phantom_region_resolution_policy(
        resolution=0.02,
        coil_resolution=0.01,
        region_resolutions={"coil": 0.03, "phantom": 0.015, "air": 0.025},
    )

    assert policy["coil_resolution_m"] == pytest.approx(0.01)
    assert policy["phantom_resolution_m"] == pytest.approx(0.015)
    assert policy["air_resolution_m"] == pytest.approx(0.025)


def test_region_resolution_policy_rejects_unknown_mapping_keys():
    with pytest.raises(ValueError, match="unsupported keys"):
        MeshGenerator.coil_phantom_region_resolution_policy(
            resolution=0.015,
            region_resolutions={"coil": 0.01, "vacuum": 0.02},
        )


# ---------------------------------------------------------------------------
# `GEO-5` (2026-09-30, §9 item 103): the policy on a real mesh.
#
# The three tests above check the resolver's dict arithmetic and build no mesh.
# The test below builds `coil_phantom_domain` on a phantom-size ladder and gates
# what the phantom request is *for*. Imports sit here, below the dict tests, so
# this section is purely additive to the file.
# ---------------------------------------------------------------------------
import math  # noqa: E402

import numpy as np  # noqa: E402
from mpi4py import MPI  # noqa: E402

from tests.mesh.helpers import (  # noqa: E402
    REQUIRED_COIL_PHANTOM_TAGS,
    assert_tag_volumes_partition_domain,
)
from tests.mesh.test_mesh_tag_integrity import (  # noqa: E402
    CAD_VOLUMES,
    GEOMETRY,
    POLICY_RESOLUTIONS,
)

# Requested phantom sizes [m]; coil 0.012 and air 0.020 are held at the
# module policy (`POLICY_RESOLUTIONS`), so rung 2 (0.010) IS the policy mesh.
GEO5_PHANTOM_LADDER = (0.015, 0.010, 0.0075, 0.005)

# Pre-registered by the 2026-09-30 03:00 review (§9 item 103), a priori:
# (iii) the inscribed-volume deficit of a cylinder with n lateral chords is
# 1 - sin(2 pi/n)/(2 pi/n) ~ (2 pi/n)^2/6, second order in the requested size;
# integer-n rounding moves each d by <= +-5 %, <= +-0.1 on a slope over ln 3.
GEO5_RATE_BAND = (1.7, 2.3)
# (iv) equivalent regular-tet edge over the requested size; priors 1.081 at
# 0.015 and 1.165 at 0.010 (`20260826T033758Z_GEO-17-mesh5-control-rechoice`).
GEO5_EDGE_RATIO_BAND = (0.8, 1.5)

# Predicted, printed, never asserted: the `GEO-17` record of the policy mesh's
# phantom volume, and the item's predicted deficits per rung.
GEO5_RUNG2_V3_RECORD = 4.990112950e-04
GEO5_PREDICTED_DEFICIT = (0.0165, 0.0072, 0.0041, 0.0018)
# The finer requests rung 1's cells are scored against in control (n2).
GEO5_N2_REQUESTS = (0.0075, 0.005)


def _owned_cell_counts(mesh, cell_tags, tags, comm):
    """Owned cells per tag, reduced.

    ``cell_tags`` is rank-local and may carry ghost cells: only indices below
    the cell index map's ``size_local`` are owned, and only those are summed.
    """
    tdim = mesh.topology.dim
    size_local = mesh.topology.index_map(tdim).size_local
    values = cell_tags.values[cell_tags.indices < size_local]
    counts = {}
    for tag in tags:
        local = int(np.count_nonzero(values == int(tag)))
        counts[int(tag)] = int(comm.allreduce(local, op=MPI.SUM))
    return counts


def _lsq_slope(x, y):
    """Least-squares slope of ``ln y`` on ``ln x``."""
    return float(np.polyfit(np.log(np.asarray(x)), np.log(np.asarray(y)), 1)[0])


def _in_band(value, band):
    return band[0] <= value <= band[1]


def test_phantom_request_is_realised_and_inscribed_deficit_falls_second_order():
    """The phantom request reaches the mesh, and its volume converges at O(h^2).

    `GEO-5` (§9 item 103). `GEO-17` gates the policy mesh by sign and a
    one-sided recovery floor; this test gates the phantom leg quantitatively
    against the cylinder's closed form. A linear-tet mesh inscribes the phantom
    (r = 0.04, H = 0.10), so ``d(h) = 1 - V_3/(pi r^2 H) > 0`` and falls at
    second order in the *requested* size ``h_p``.

    Asserted: (i) the tags partition every rung's mesh to 1e-9; (ii) d > 0 and
    strictly decreasing; (iii) LSQ slope of ln d on ln h_p in [1.7, 2.3];
    (iv) a_eq/h_p in [0.8, 1.5] at every rung, a_eq = (6 sqrt(2) V_3/N_3)^(1/3)
    the edge of the regular tet with the phantom's mean cell volume.

    Negative controls, asserted only as "misses the band": (n1) the coil's
    deficit (tag 1, whose request is held) regressed on the same ladder;
    (n2) rung 1's a_eq scored against 0.0075 and 0.005 — the reading a
    request that was never applied (`GEO-17`'s before-state) would give.

    Phantom leg only: the coil torus is pre-asymptotic at every priced size and
    the air has no closed form. Meshing only — no solve.
    """
    comm = MPI.COMM_WORLD
    tags = REQUIRED_COIL_PHANTOM_TAGS
    v_cyl = CAD_VOLUMES[3]

    rungs = []
    for h_p in GEO5_PHANTOM_LADDER:
        resolutions = dict(POLICY_RESOLUTIONS)
        resolutions["phantom_resolution"] = h_p
        mesh, cell_tags, _ = MeshGenerator.coil_phantom_domain(comm=comm, **resolutions, **GEOMETRY)
        # (i) — asserted inside the helper, at every rung.
        volumes = assert_tag_volumes_partition_domain(
            mesh, cell_tags, tags, comm=comm, label=f"GEO-5 h_p={h_p:g}"
        )
        counts = _owned_cell_counts(mesh, cell_tags, tags, comm)
        a_eq = (6.0 * math.sqrt(2.0) * volumes[3] / counts[3]) ** (1.0 / 3.0)
        rungs.append(
            dict(
                h_p=h_p,
                volumes=volumes,
                counts=counts,
                d=1.0 - volumes[3] / v_cyl,
                d_coil=1.0 - volumes[1] / CAD_VOLUMES[1],
                a_eq=a_eq,
            )
        )

    h = [r["h_p"] for r in rungs]
    d = [r["d"] for r in rungs]
    nan = float("nan")
    rate = _lsq_slope(h, d) if min(d) > 0.0 else nan
    pairwise = [
        (
            math.log(d[i] / d[i + 1]) / math.log(h[i] / h[i + 1])
            if d[i] > 0.0 and d[i + 1] > 0.0
            else nan
        )
        for i in range(len(rungs) - 1)
    ]
    d_coil = [r["d_coil"] for r in rungs]
    coil_rate = _lsq_slope(h, d_coil) if min(d_coil) > 0.0 else nan
    n2_ratios = [rungs[0]["a_eq"] / h_other for h_other in GEO5_N2_REQUESTS]

    # Every record is printed before anything beyond (i) is asserted, so a red
    # band still leaves V_3, N_3 and the pairwise rates in the log.
    if comm.rank == 0:
        print(f"\n[GEO-5] pi r^2 H = {v_cyl:.9e} m^3", flush=True)
        for r, d_pred in zip(rungs, GEO5_PREDICTED_DEFICIT):
            vol, cnt = r["volumes"], r["counts"]
            print(
                f"[GEO-5] h_p {r['h_p']:.4f}: N_3 {cnt[3]:6d}, "
                f"V_3 {vol[3]:.9e}, d {r['d']:.6%} (pred. {d_pred:.2%}), "
                f"a_eq {r['a_eq']:.6f}, a_eq/h_p {r['a_eq'] / r['h_p']:.4f} | "
                f"V_1 {vol[1]:.9e}, V_2 {vol[2]:.9e}, V_4 {vol[4]:.9e} | "
                f"N_1 {cnt[1]}, N_2 {cnt[2]}, N_4 {cnt[4]}, "
                f"N_total {sum(cnt.values())}",
                flush=True,
            )
        print(
            f"[GEO-5] LSQ rate of ln d on ln h_p: {rate:.4f} (pred. 2.0, band "
            f"{GEO5_RATE_BAND}); pairwise " + " / ".join(f"{p:.4f}" for p in pairwise),
            flush=True,
        )
        v3_policy = rungs[1]["volumes"][3]
        print(
            f"[GEO-5] rung 2 (policy mesh) V_3 {v3_policy:.9e} vs GEO-17 record "
            f"{GEO5_RUNG2_V3_RECORD:.9e} "
            f"(rel {v3_policy / GEO5_RUNG2_V3_RECORD - 1.0:+.3e})",
            flush=True,
        )
        print(
            "[GEO-5] control n1 (coil deficit, tag 1, on the phantom ladder): "
            + " / ".join(f"{x:.6%}" for x in d_coil)
            + f"; slope {coil_rate:+.4f} (pred. 0)",
            flush=True,
        )
        print(
            f"[GEO-5] control n2 (rung 1 a_eq {rungs[0]['a_eq']:.6f} scored "
            f"against 0.0075 / 0.005): {n2_ratios[0]:.4f} / {n2_ratios[1]:.4f} "
            "(pred. 2.16 / 3.24)",
            flush=True,
        )

    # (ii) inscribed: never past the CAD volume, and strictly closer with h.
    for r in rungs:
        assert r["d"] > 0.0, (
            f"h_p = {r['h_p']}: phantom volume {r['volumes'][3]:.9e} m^3 is not "
            f"below pi r^2 H = {v_cyl:.9e} m^3 — an inscribing linear-tet mesh "
            "cannot reach the CAD volume"
        )
    for coarse, fine in zip(rungs[:-1], rungs[1:]):
        assert fine["d"] < coarse["d"], (
            f"phantom deficit did not fall from h_p = {coarse['h_p']} "
            f"({coarse['d']:.6%}) to h_p = {fine['h_p']} ({fine['d']:.6%})"
        )

    # (iii) second order in the requested size.
    assert _in_band(rate, GEO5_RATE_BAND), (
        f"LSQ rate {rate:.4f} of the phantom deficit on the requested size is "
        f"outside {GEO5_RATE_BAND}; pairwise {pairwise}"
    )

    # (iv) the request is realised: the phantom's cells are the size asked for.
    for r in rungs:
        ratio = r["a_eq"] / r["h_p"]
        assert _in_band(ratio, GEO5_EDGE_RATIO_BAND), (
            f"h_p = {r['h_p']}: equivalent edge {r['a_eq']:.6f} m is "
            f"{ratio:.4f} x the request, outside {GEO5_EDGE_RATIO_BAND} "
            f"(N_3 = {r['counts'][3]}, V_3 = {r['volumes'][3]:.9e})"
        )

    # (n1) region-blind comparand: the coil's deficit does not follow the
    # phantom ladder, so its slope must MISS the rate band.
    assert not _in_band(coil_rate, GEO5_RATE_BAND), (
        f"control n1: the coil's deficit regressed on the phantom ladder has "
        f"slope {coil_rate:.4f}, inside {GEO5_RATE_BAND} — the rate gate does "
        "not distinguish the phantom from a region whose size was held"
    )

    # (n2) request-not-applied: rung 1's cells scored against the finer
    # requests must MISS the edge band.
    for h_other, ratio in zip(GEO5_N2_REQUESTS, n2_ratios):
        assert not _in_band(ratio, GEO5_EDGE_RATIO_BAND), (
            f"control n2: the h_p = 0.015 mesh scored against a request of "
            f"{h_other} reads {ratio:.4f}, inside {GEO5_EDGE_RATIO_BAND} — the "
            "edge gate would pass a request that was never applied"
        )
