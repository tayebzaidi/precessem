"""`IMP-3` step 5 (§9 item 252, chain I13): Palace order 2 against ours at degree 2.

Item 235's two exact-rectangle builds (off, rod PEC), degree 2, κ off,
``keep_fields=True``.  Per arm and drive: E at item 242's six points through
``evaluate_vector_field_parallel(..., owner="lowest_original_cell")`` and the
phantom ∫ e·ē dx (quadrature degree 4: exact for |e|² of a degree-2 field on
affine cells), allreduced; the global dof count.

Asserted: (F) cells = ``imp3x_meta.json`` and N_i areas; (R2) the off arm's S
reproduces ``logs/ops61/ours_f05x_p2.json`` ``s_kappa_off`` within
``REUSE_REPRODUCTION_RTOL`` (file ranks == comm.size, degree == 2); (i)/(ii)
on both sweeps; (E) rod PEC wall on the degree-2 wall dofs; (V) all six points
valid.  Printed, never gated: (D), (B), (T_E2), (T_W2), (T_S2), the order 1 -> 2
shifts per code, the closed form beside η_rod, (C1')-(C3').
"""

from __future__ import annotations

import gc
import json

import numpy as np
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.ports.sparameters import run_n_port_sparameter_sweep
from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.utils.analytical import pec_prolate_spheroid_field
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE
from tests.validation.test_birdcage_power_identity import STEP2D_FREQUENCY_HZ
from tests.validation.test_imp3_palace_fields import FP, PNAMES, PTS, _owned_hits
from tests.validation.test_imp3_palace_twin import _ids, _sheet_areas
from tests.validation.test_imp3_rod_birdcage import (
    PEC_WALL_MAX_E,
    PORT_IDS,
    _matrix_quantities,
    _wall_dofs,
    _wall_max_e,
)
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL

NOMINAL_STRIP_AREA = 5.6e-5
ARMS = ("off", "rod")
TAG = "[IMP-3 5 W5]"


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _sweep2(case):
    return run_n_port_sparameter_sweep(
        case["problem"],
        case["port_defs"],
        lumped_sheet_ports=case["specs"],
        lumped_sheet_facet_tags=case["facet_tags"],
        keep_fields=True,
        degree=2,
    )


def _readings2(res, msh, cell_tags, comm):
    """E (6, 3, 4) and W (4,) for every drive; quadrature 4 for degree 2."""
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags, metadata={"quadrature_degree": 4})
    e = np.zeros((len(PTS), 3, len(PORT_IDS)), dtype=complex)
    w = np.zeros(len(PORT_IDS))
    valid = True
    for j, p in enumerate(PORT_IDS):
        f = res.fields[p].e_complex
        vals, mask = evaluate_vector_field_parallel(f, PTS, owner="lowest_original_cell")
        vals = comm.bcast(vals if comm.rank == 0 else None, root=0)
        mask = comm.bcast(mask if comm.rank == 0 else None, root=0)
        valid = valid and bool(np.all(mask))
        e[:, :, j] = np.asarray(vals).reshape(len(PTS), 3)
        form = fem.form(ufl.inner(f, f) * dx(PHANTOM_CELL_TAG))
        w[j] = comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM).real
    return e, w, valid


def _cplx(arr):
    a = np.array(arr)
    return a[..., 0] + 1j * a[..., 1]


@complex_only
def test_imp3_palace_fields_degree2():
    comm = MPI.COMM_WORLD
    s3 = FP._step3()
    out = s3.OUT_DIR
    meta = json.loads((out / "imp3x_meta.json").read_text())
    rec = json.loads((FP.REPO / "logs" / "ops61" / "ours_f05x_p2.json").read_text())
    assert rec["ranks"] == comm.size, ("(R2) ranks", rec["ranks"], comm.size)
    assert rec["degree"] == 2, ("(R2) degree", rec["degree"])
    E, W, B, S, ND = {}, {}, {}, {}, {}
    wall = {}
    fails = []
    for arm in ARMS:
        built = s3.build_arm(arm)
        case = build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True)
        msh = case["mesh"]
        ncells = int(case["cells"])
        areas = _sheet_areas(msh, case["facet_tags"], comm)
        _say(
            f"\n{TAG} arm {arm}: {ncells} cells (meta {meta[arm]['cells']}); N_i areas "
            + ", ".join(f"P{i} {a:.12e}" for i, a in areas.items())
        )
        assert ncells == meta[arm]["cells"], (arm, ncells)  # (F)
        for i, a in areas.items():
            assert abs(a - NOMINAL_STRIP_AREA) <= STRIP_BAND, (arm, i, a)
        B[arm] = _owned_hits(msh, comm)
        _say(
            f"{TAG} (B) {arm} owned-cell hits per point "
            + "  ".join(f"{n} {int(h)}" for n, h in zip(PNAMES, B[arm]))
            + " (predicted 1)"
        )
        res = _sweep2(case)
        assert tuple(res.port_ids) == PORT_IDS
        mq = _matrix_quantities(res)
        _ids(arm, mq)
        S[arm] = mq["s"]
        if arm == "off":
            ref = np.array([[complex(*z) for z in row] for row in rec["s_kappa_off"]])
            r = float(np.max(np.abs(S[arm] - ref)) / np.max(np.abs(ref)))
            _say(
                f"{TAG} (R2) off S vs ours_f05x_p2.json s_kappa_off entrywise {r:.3e} "
                f"<= {REUSE_REPRODUCTION_RTOL}"
            )
            assert r <= REUSE_REPRODUCTION_RTOL, ("R2", r)
        v_space = res.fields["P1"].e_complex.function_space
        ND[arm] = int(v_space.dofmap.index_map.size_global * v_space.dofmap.index_map_bs)
        _say(
            f"{TAG} (D) {arm} global dofs {ND[arm]} (Palace ND predicted "
            f"{753062 if arm == 'off' else 1090614})"
        )
        E[arm], W[arm], valid = _readings2(res, msh, built[1], comm)
        _say(f"{TAG} (V) {arm} valid_mask all True: {valid}")
        if not valid:
            fails.append(f"V {arm}")
        if arm == "rod":
            owned, n_wall = _wall_dofs(msh, built[2], v_space)
            assert n_wall > 0
            wall = {p: _wall_max_e(res.fields[p].e_complex, owned, comm) for p in PORT_IDS}
        del res, v_space, case, built, msh
        gc.collect()
    _say(f"{TAG} (E) rod PEC wall max |E| {wall} <= {PEC_WALL_MAX_E}")
    for p, v in wall.items():
        assert v <= PEC_WALL_MAX_E, (p, v)
    assert not fails, fails

    if comm.rank == 0:
        with open(out / "ours_imp3q_fields.json", "w") as fh:
            json.dump(
                {
                    "E": {
                        a: [
                            [[[z.real, z.imag] for z in E[a][k, c, :]] for c in range(3)]
                            for k in range(len(PTS))
                        ]
                        for a in E
                    },
                    "W": {a: W[a].tolist() for a in W},
                    "ND": ND,
                },
                fh,
            )

    # ---------------------------------------------------------------- printed
    pq = json.loads((out / "palace_imp3q_fields.json").read_text())
    p1 = json.loads((out / "palace_imp3p_fields.json").read_text())
    o1 = json.loads((out / "ours_imp3p_fields.json").read_text())
    assert pq["probe_names"] == list(PNAMES)
    assert np.array_equal(np.array(pq["points"]), PTS)
    EP = {a: _cplx(pq["E"][a]) for a in ARMS}
    WP = {a: np.array(pq["W_elec"][a]) for a in ARMS}
    EP1 = {a: _cplx(p1["E"][a]) for a in ARMS}
    WP1 = {a: np.array(p1["W_elec"][a]) for a in ARMS}
    EO1 = {a: _cplx(o1["E"][a]) for a in ARMS}
    WO1 = {a: np.array(o1["W"][a]) for a in ARMS}
    j = 0
    iF1, iT1 = PNAMES.index("F1"), PNAMES.index("T1")
    cstar = int(np.argmax(np.abs(EP["off"][iF1, :, j])))
    g2 = E["off"][iF1, cstar, j] / EP["off"][iF1, cstar, j]
    g2c = E["off"][iF1, cstar, j] / np.conj(EP["off"][iF1, cstar, j])
    g1 = EO1["off"][iF1, cstar, j] / EP1["off"][iF1, cstar, j]
    _say(
        f"{TAG} g2 = {g2.real:+.9e}{g2.imag:+.9e}j (c* {'xyz'[cstar]}); order-1 g "
        f"{g1.real:+.9e}{g1.imag:+.9e}j; |g2/g - 1| {abs(g2 / g1 - 1):.3e} "
        f"(predicted <= 1e-6)"
    )
    amb = {a: [k for k in range(len(PTS)) if int(B[a][k]) != 1] for a in ARMS}

    def tie(a, ours, pal, gg, jj=j, conj=False):
        worst, where = 0.0, None
        for k in range(len(PTS)):
            if k in amb[a]:
                continue
            ref = gg * (np.conj(pal[k, :, jj]) if conj else pal[k, :, jj])
            nrm = np.linalg.norm(ref)
            for c in range(3):
                if a == "off" and k == iF1 and c == cstar:
                    continue
                v = abs(ours[k, c, jj] - ref[c]) / nrm
                if v > worst:
                    worst, where = v, f"{PNAMES[k]}.{'xyz'[c]}"
        return worst, where

    te = {}
    for a in ARMS:
        te[a], where = tie(a, E[a], EP[a], g2)
        tc, wc = tie(a, E[a], EP[a], g2c, conj=True)
        per = [tie(a, E[a], EP[a], g2, jj=jj)[0] for jj in range(len(PORT_IDS))]
        _say(
            f"{TAG} (T_E2) {a} P1 max rel {te[a]:.3e} at {where} (predicted <= 1e-6); conj "
            f"{tc:.3e} at {wc} (predicted O(1)); P1-P4 " + ", ".join(f"{v:.3e}" for v in per)
        )
    q = {a: W[a][j] / WP[a][j] for a in ARMS}
    _say(
        f"{TAG} (T_W2) q_rod/q_off - 1 {q['rod'] / q['off'] - 1:+.3e} (predicted <= 1e-6); "
        f"q_off/|g2|^2 {q['off'] / abs(g2) ** 2:.9e} (item 242: 2.896e9)"
    )
    for a in ARMS:
        _say(
            f"{TAG} (T_S2) {a} ours S P1 column "
            + ", ".join(f"{z.real:+.6e}{z.imag:+.6e}j" for z in S[a][:, 0])
        )

    mo2 = {a: np.linalg.norm(E[a], axis=1) for a in ARMS}
    mp2 = {a: np.linalg.norm(EP[a], axis=1) for a in ARMS}
    mo1 = {a: np.linalg.norm(EO1[a], axis=1) for a in ARMS}
    mp1 = {a: np.linalg.norm(EP1[a], axis=1) for a in ARMS}
    for a in ARMS:
        do = mo2[a][:, j] / mo1[a][:, j] - 1
        dp = mp2[a][:, j] / mp1[a][:, j] - 1
        _say(
            f"{TAG} δ_p {a} P1 ours / Palace / gap: "
            + "  ".join(
                f"{n} {do[k]:+.6e}/{dp[k]:+.6e}/{abs(do[k] - dp[k]):.2e}"
                for k, n in enumerate(PNAMES)
            )
        )
    eo2, ep2 = mo2["rod"] / mo2["off"], mp2["rod"] / mp2["off"]
    eo1, ep1 = mo1["rod"] / mo1["off"], mp1["rod"] / mp1["off"]
    for jj in range(len(PORT_IDS)):
        _say(
            f"{TAG} η_rod^o2 P{jj + 1} ours / Palace / d_η/η: "
            + "  ".join(
                f"{n} {eo2[k, jj]:.6f}/{ep2[k, jj]:.6f}/{abs(eo2[k, jj] - ep2[k, jj]) / ep2[k, jj]:.2e}"
                for k, n in enumerate(PNAMES)
            )
        )
        _say(
            f"{TAG} η_rod^o2/η_rod^o1 - 1 P{jj + 1} ours / Palace: "
            + "  ".join(
                f"{n} {eo2[k, jj] / eo1[k, jj] - 1:+.6e}/{ep2[k, jj] / ep1[k, jj] - 1:+.6e}"
                for k, n in enumerate(PNAMES)
            )
        )
        dwo2 = W["rod"][jj] / W["off"][jj] - 1
        dwp2 = WP["rod"][jj] / WP["off"][jj] - 1
        dwo1 = WO1["rod"][jj] / WO1["off"][jj] - 1
        dwp1 = WP1["rod"][jj] / WP1["off"][jj] - 1
        _say(
            f"{TAG} ΔW(rod/off) P{jj + 1}: o2 ours {dwo2:+.6e} Palace {dwp2:+.6e} gap "
            f"{abs(dwo2 - dwp2):.2e}; o1 ours {dwo1:+.6e} Palace {dwp1:+.6e}"
        )
    cf = np.linalg.norm(
        pec_prolate_spheroid_field(
            PTS - np.array(ROD_CENTRE), semi_major=0.016, semi_minor=0.001, e0=1.0
        ),
        axis=1,
    )
    _say(
        f"{TAG} closed form |E_cf|/E0 vs η_rod^o2(p)/η_rod^o2(F1) (ours, P1), ratio "
        "(predicted [0.3, 3]): "
        + "  ".join(
            f"{PNAMES[k]} {cf[k]:.4f}/{eo2[k, j] / eo2[iF1, j]:.4f}="
            f"{cf[k] / (eo2[k, j] / eo2[iF1, j]):.3f}"
            for k in (0, 1, 2, 3)
        )
    )
    shift = abs(ep2[iT1, j] / ep1[iT1, j] - 1)
    d1 = abs(eo1[iT1, j] - ep1[iT1, j]) / ep1[iT1, j]
    d2 = abs(eo2[iT1, j] - ep2[iT1, j]) / ep2[iT1, j]
    dmax = max(d1, d2)
    _say(
        f"{TAG} (C1') |η^o2_rod(T1)/η^o1_rod(T1) - 1| / max(d_η^o1, d_η^o2)/η = "
        f"{shift:.6e} / {dmax:.3e} = {shift / dmax if dmax > 0 else float('inf'):.3e}x "
        f"(predicted >= 1e3, ceiling ≈ 1e8)"
    )
    ref = g2 * EP["rod"][iT1, :, j]
    c2 = np.linalg.norm(E["off"][iT1, :, j] - ref) / np.linalg.norm(ref)
    _say(
        f"{TAG} (C2') |E^ours,o2_off(T1) - g2·E^Pal,o2_rod(T1)| / ‖·‖ = {c2:.3e} "
        f"(predicted >= 0.3)"
    )
    ref = g1 * EP1["rod"][iT1, :, j]
    c3 = np.linalg.norm(E["rod"][iT1, :, j] - ref) / np.linalg.norm(ref)
    _say(
        f"{TAG} (C3') |E^ours,o2_rod(T1) - g·E^Pal,o1_rod(T1)| / ‖·‖ = {c3:.3e} "
        f"(predicted >= 1e-3); over (T_E2 rod) "
        f"{c3 / te['rod'] if te['rod'] > 0 else float('inf'):.3e}x"
    )
