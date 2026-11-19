"""`IMP-3` step 4 (§9 item 242, chain I12): the near-tip field, ours against Palace.

Same two exact-rectangle builds and same three sweeps as
``test_imp3_palace_twin.py`` (item 235), degree 1, κ off, ``keep_fields=True``
on all three.  Per arm and drive: E at six fixed phantom points through
``evaluate_vector_field_parallel(..., owner="lowest_original_cell")`` and the
phantom ∫ e·ē dx (quadrature degree pinned, allreduced).

Asserted: (F) cells = ``imp3x_meta.json`` exactly and N_i areas; (R) S per arm
reproduces item 235's ``ours_imp3x_p1.json`` within ``REUSE_REPRODUCTION_RTOL``;
(i)/(ii) on all three sweeps; (E) rod PEC wall; (Kw) natural wall P1; (V) all
six points valid on every arm.  Printed, never gated: (B) owned-cell counts,
(T_E) field tie, (T_W) energy tie, η_wall / η_rod per code, ΔW, the closed form
beside η, and the negative controls (C1)-(C3).  Palace values are read from
``logs/imp3/palace_imp3p_fields.json`` (written by the probe's W5).
"""

from __future__ import annotations

import gc
import importlib.util
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import ufl
from dolfinx import fem, geometry
from mpi4py import MPI

from fem_em_solver.post.evaluation import evaluate_vector_field_parallel
from fem_em_solver.utils.analytical import pec_prolate_spheroid_field
from tests.complex_mode import complex_only
from tests.mesh.test_birdcage_implant_rod import ROD_CENTRE
from tests.validation.test_birdcage_power_identity import STEP2D_FREQUENCY_HZ
from tests.validation.test_imp3_palace_twin import _ids, _sheet_areas, _sweep
from tests.validation.test_imp3_rod_birdcage import (
    PEC_WALL_MAX_E,
    PORT_IDS,
    _matrix_quantities,
    _outer_box_tags,
    _wall_dofs,
    _wall_max_e,
)
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND
from tests.validation.test_port_birdcage_four_port import build_four_port_sweep
from tests.validation.test_port_birdcage_lumped_column import PHANTOM_CELL_TAG
from tests.validation.test_port_birdcage_ring_matrix import REUSE_REPRODUCTION_RTOL
from tests.validation.test_th14_birdcage_copper import OUTER_BOX_TAG

REPO = Path(__file__).resolve().parents[2]
NOMINAL_STRIP_AREA = 5.6e-5
ARMS = ("off", "rod", "rodnat")
RUN = {"off": "imp3x_off", "rod": "imp3x_rod", "rodnat": "imp3x_rodnat"}


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, REPO / rel)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


FP = _load("imp3_step4_palace_fields", "scripts/probes/imp3_step4_palace_fields.py")
PNAMES = FP.PROBE_NAMES
PTS = np.ascontiguousarray(np.array(FP.PROBE_POINTS, dtype=np.float64))


def _say(msg):
    if MPI.COMM_WORLD.rank == 0:
        print(msg, flush=True)


def _owned_hits(msh, comm):
    """(B) per point: number of owned cells colliding, summed over ranks."""
    tree = geometry.bb_tree(msh, msh.topology.dim)
    cand = geometry.compute_collisions_points(tree, PTS)
    coll = geometry.compute_colliding_cells(msh, cand, PTS)
    nloc = msh.topology.index_map(msh.topology.dim).size_local
    local = np.array([int(np.sum(coll.links(i) < nloc)) for i in range(len(PTS))], dtype=np.int64)
    return comm.allreduce(local, op=MPI.SUM)


def _readings(res, msh, cell_tags, comm):
    """E (6, 3, 4) and W (4,) for every drive of one sweep."""
    dx = ufl.Measure("dx", domain=msh, subdomain_data=cell_tags, metadata={"quadrature_degree": 2})
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


def _pal_fields():
    d = json.loads((FP._step3().OUT_DIR / FP.FIELDS_JSON).read_text())
    e = {
        a: np.array([[[complex(*z) for z in comp] for comp in pt] for pt in d["E"][a]])
        for a in ARMS
    }
    w = {a: np.array(d["W_elec"][a]) for a in ARMS}
    assert d["probe_names"] == list(PNAMES)
    assert np.array_equal(np.array(d["points"]), PTS)
    return e, w


@complex_only
def test_imp3_palace_fields_degree1():
    comm = MPI.COMM_WORLD
    s3 = FP._step3()
    meta = json.loads((s3.OUT_DIR / "imp3x_meta.json").read_text())
    prior = json.loads((s3.OUT_DIR / "ours_imp3x_p1.json").read_text())
    E, W, B, S = {}, {}, {}, {}
    wall = {}
    fails = []

    def check_r(name, s):
        ref = np.array([[complex(*z) for z in row] for row in prior[name]])
        r = float(np.max(np.abs(s - ref)) / np.max(np.abs(ref)))
        _say(
            f"[IMP-3 4 W6] (R) {name} S vs ours_imp3x_p1.json entrywise {r:.3e} "
            f"<= {REUSE_REPRODUCTION_RTOL}"
        )
        assert r <= REUSE_REPRODUCTION_RTOL, (name, r)

    for arm in ("off", "rod"):
        built = s3.build_arm(arm)
        case = build_four_port_sweep(frequency_hz=STEP2D_FREQUENCY_HZ, built=built, build_only=True)
        msh = case["mesh"]
        ncells = int(case["cells"])
        areas = _sheet_areas(msh, case["facet_tags"], comm)
        _say(
            f"\n[IMP-3 4 W6] arm {arm}: {ncells} cells (meta {meta[arm]['cells']}); N_i areas "
            + ", ".join(f"P{i} {a:.12e}" for i, a in areas.items())
        )
        assert ncells == meta[arm]["cells"], (arm, ncells)  # (F)
        for i, a in areas.items():
            assert abs(a - NOMINAL_STRIP_AREA) <= STRIP_BAND, (arm, i, a)
        hits = _owned_hits(msh, comm)
        _say(
            f"[IMP-3 4 W6] (B) {arm} mesh owned-cell hits per point "
            + "  ".join(f"{n} {int(h)}" for n, h in zip(PNAMES, hits))
            + " (predicted 1)"
        )
        labels = ("off",) if arm == "off" else ("rod", "rodnat")
        for lab in labels:
            B[lab] = hits
            if lab == "rodnat":
                outer, census = _outer_box_tags(msh, built[2], comm)
                _say(f"[IMP-3 4 W6] natural-wall census {census}")
                assert census["tag_collision"] == 0, census
                assert census["exterior"] == census["outer"] + census["wall_exterior"], census
                problem = replace(
                    case["problem"], facet_tags=outer, pec_facet_tags=(OUTER_BOX_TAG,)
                )
                res = _sweep(case, problem=problem, keep=True)
            else:
                res = _sweep(case, keep=True)
            assert tuple(res.port_ids) == PORT_IDS
            mq = _matrix_quantities(res)
            _ids(lab, mq)
            check_r(RUN[lab], mq["s"])
            S[lab] = mq["s"]
            E[lab], W[lab], valid = _readings(res, msh, built[1], comm)
            _say(f"[IMP-3 4 W6] (V) {lab} valid_mask all True: {valid}")
            if not valid:
                fails.append(f"V {lab}")
            if arm == "rod":
                v_space = res.fields["P1"].e_complex.function_space
                owned, n_wall = _wall_dofs(msh, built[2], v_space)
                assert n_wall > 0
                wall[lab] = {p: _wall_max_e(res.fields[p].e_complex, owned, comm) for p in PORT_IDS}
            del res
            gc.collect()
        del case, built, msh
        gc.collect()
    _say(f"[IMP-3 4 W6] (E) rod PEC wall max |E| {wall['rod']} <= {PEC_WALL_MAX_E}")
    _say(f"[IMP-3 4 W6] (Kw) natural wall max |E| {wall['rodnat']} (P1 > {PEC_WALL_MAX_E})")
    for p, v in wall["rod"].items():
        assert v <= PEC_WALL_MAX_E, (p, v)
    assert wall["rodnat"]["P1"] > PEC_WALL_MAX_E, wall["rodnat"]
    assert not fails, fails

    if comm.rank == 0:
        with open(s3.OUT_DIR / "ours_imp3p_fields.json", "w") as fh:
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
                },
                fh,
            )

    # ---------------------------------------------------------------- printed
    EP, WP = _pal_fields()
    j = 0  # P1, the reading of record
    iF1 = PNAMES.index("F1")
    cstar = int(np.argmax(np.abs(EP["off"][iF1, :, j])))
    g = E["off"][iF1, cstar, j] / EP["off"][iF1, cstar, j]
    gc_ = E["off"][iF1, cstar, j] / np.conj(EP["off"][iF1, cstar, j])
    _say(
        f"[IMP-3 4 W6] g = E^ours_{'xyz'[cstar]}(F1)/E^Pal_{'xyz'[cstar]}(F1) (off, P1) = "
        f"{g.real:+.9e}{g.imag:+.9e}j  |g| {abs(g):.9e}"
    )
    ambiguous = {a: [k for k in range(len(PTS)) if int(B[a][k]) != 1] for a in ARMS}
    _say(
        "[IMP-3 4 W6] (B) face-ambiguous points per arm "
        + str({a: [PNAMES[k] for k in v] for a, v in ambiguous.items()})
    )

    def tie(a, ours, pal, jj=j, conj=False):
        worst, where = 0.0, None
        gg = gc_ if conj else g
        for k in range(len(PTS)):
            if k in ambiguous[a]:
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
        te[a], where = tie(a, E[a], EP[a])
        tc, wc = tie(a, E[a], EP[a], conj=True)
        per_drive = [tie(a, E[a], EP[a], jj=jj)[0] for jj in range(len(PORT_IDS))]
        _say(
            f"[IMP-3 4 W6] (T_E) {a} P1 max rel {te[a]:.3e} at {where} (predicted <= 1e-6); "
            f"conj form {tc:.3e} at {wc} (predicted O(1)); P1-P4 "
            + ", ".join(f"{v:.3e}" for v in per_drive)
        )
    q = {a: W[a][j] / WP[a][j] for a in ARMS}
    tw = max(abs(q[a] / q["off"] - 1.0) for a in ARMS)
    _say(
        "[IMP-3 4 W6] (T_W) q_a/q_off - 1 "
        + "  ".join(f"{a} {q[a] / q['off'] - 1:+.3e}" for a in ARMS)
        + f"; max {tw:.3e} (predicted <= 1e-6); q_off/|g|^2 {q['off'] / abs(g) ** 2:.9e}"
    )

    mo = {a: np.linalg.norm(E[a], axis=1) for a in ARMS}
    mp = {a: np.linalg.norm(EP[a], axis=1) for a in ARMS}
    eta = {}
    for lab, num, den in (("wall", "rod", "rodnat"), ("rod", "rod", "off")):
        eo, ep = mo[num] / mo[den], mp[num] / mp[den]
        eta[lab] = (eo, ep)
        for jj in range(len(PORT_IDS)):
            _say(
                f"[IMP-3 4 W6] η_{lab} P{jj + 1} ours / Palace / d_η/η: "
                + "  ".join(
                    f"{n} {eo[k, jj]:.6f}/{ep[k, jj]:.6f}/"
                    f"{abs(eo[k, jj] - ep[k, jj]) / ep[k, jj]:.2e}"
                    for k, n in enumerate(PNAMES)
                )
            )
    _say(
        "[IMP-3 4 W6] η_wall predicted (P1): T1 [1.5, 20], T2 [1.2, 8], T3 [1.05, 4], "
        "S1 [0.2, 1.5], F1/F2 |η-1| in [3e-3, 1e-1]; η_rod F1/F2 |η-1| in [1e-3, 0.3]; "
        "d_η <= 2e-6·η"
    )
    for jj in range(len(PORT_IDS)):
        parts = []
        for lab, num, den in (("rod/rodnat", "rod", "rodnat"), ("rod/off", "rod", "off")):
            do = W[num][jj] / W[den][jj] - 1.0
            dp = WP[num][jj] / WP[den][jj] - 1.0
            parts.append(f"{lab} - 1 ours {do:+.6e} Palace {dp:+.6e} gap {abs(do - dp):.2e}")
        _say(
            f"[IMP-3 4 W6] ΔW P{jj + 1}: "
            + "; ".join(parts)
            + " (predicted |Δ| in [1e-3, 0.1], gap <= 1e-6)"
        )
    cf = np.linalg.norm(
        pec_prolate_spheroid_field(
            PTS - np.array(ROD_CENTRE), semi_major=0.016, semi_minor=0.001, e0=1.0
        ),
        axis=1,
    )
    eo_rod = eta["rod"][0][:, j]
    _say(
        "[IMP-3 4 W6] closed form |E_cf|/E0 vs η_rod(p)/η_rod(F1) (ours, P1), ratio "
        "(predicted [0.3, 3]): "
        + "  ".join(
            f"{PNAMES[k]} {cf[k]:.4f}/{eo_rod[k] / eo_rod[iF1]:.4f}="
            f"{cf[k] / (eo_rod[k] / eo_rod[iF1]):.3f}"
            for k in (0, 1, 2, 3)
        )
    )

    iT1 = PNAMES.index("T1")
    eo_w, ep_w = eta["wall"]
    d_eta = abs(eo_w[iT1, j] - ep_w[iT1, j])
    c1 = abs(ep_w[iT1, j] - 1.0) / d_eta if d_eta > 0 else float("inf")
    _say(
        f"[IMP-3 4 W6] (C1) |η^Pal_wall(T1)-1| / d_η(T1) = {abs(ep_w[iT1, j] - 1):.6f} / "
        f"{d_eta:.3e} = {c1:.3e}x (predicted >= 1e5, ceiling ≈ 2.4e9); tie_off floor "
        f"(T_E off) {te['off']:.3e}"
    )
    far = [PNAMES.index("F1"), PNAMES.index("F2")]
    for code, (eo, ep) in (("ours", eta["wall"]), ("Palace", eta["wall"])):
        e_ = eo if code == "ours" else ep
        c2 = abs(e_[iT1, j] - 1) / max(abs(e_[k, j] - 1) for k in far)
        _say(
            f"[IMP-3 4 W6] (C2) {code} |η_wall(T1)-1| / max_F|η_wall(F)-1| = {c2:.3f}x "
            f"(predicted >= 5, ceiling ≈ 6e3)"
        )
    ref = g * EP["rod"][iT1, :, j]
    c3 = np.linalg.norm(E["off"][iT1, :, j] - ref) / np.linalg.norm(ref)
    _say(
        f"[IMP-3 4 W6] (C3) |E^ours_off(T1) - g·E^Pal_rod(T1)| / ‖g·E^Pal_rod(T1)‖ = {c3:.3e} "
        f"(predicted >= 0.3); separation over (T_E off) "
        f"{c3 / te['off'] if te['off'] > 0 else float('inf'):.3e}x"
    )
