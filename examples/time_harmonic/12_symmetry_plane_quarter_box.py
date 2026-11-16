"""Example (`EX-66`): quarter the box — a PMC plane and a PEC plane, the port rescaled.

No example cuts a model on a symmetry plane. This one solves one saline box
with one interior lumped sheet (ĥ = ŷ, 50 Ω, 1 V, 64 MHz, degree 1) four ways
on mirror-exact meshes — full, PMC half, PEC half, quarter — and shows the
rule (Jin §8.7.5):

* a cut plane **containing** ĥ (x = 0.24) is PMC: the natural condition, the
  half-sheet sees Z doubled and I halved (``full_impedance(·, "pmc")``);
* a cut plane **normal** to ĥ (y = 0.20) is PEC: pinned, Z and V halve
  (``full_impedance(·, "pec")``);
* the quarter (both cuts) returns Z_f on a quarter of the cells.

**It asserts, and it does not re-implement** (the `ANS-1` rule). The meshes,
specs, solve, points, mirror maps, counts and bands are imported from
``tests/validation/test_th16_symmetry_plane_driven.py`` (`TH-16` step 2). The
gate gains, additively, ``_solve(..., return_fields=False)`` and a
module-level ``build_cases()`` (the fixture's case table, moved verbatim).

* *Structure:* vertices / cells / owned sheet facets ``== EXPECTED``, PMC-cut
  facets ``== EXPECTED_PMC_CUT_FACETS``.
* *Identities, each* ``<= LEG_D0_REPRODUCTION_BAND``: (1) the full field is
  mirror-symmetric at ``POINTS`` and their three images; (2) PMC half Z, V,
  I/2, P/2, E; (3) PEC half Z, V/2, I, P/2; (4) quarter Z, P/4.
* *Negative control:* the PMC half with its cut pinned misses 2Z_f by more
  than ``NC1_FLOOR``.

Run it::

    ./run_examples.sh -e th:12

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_12_symmetry_plane_quarter_box_{full,pmc_half,pec_half,quarter}_combined.xdmf``
(CG1 ``E_real``, ``E_imag``; facet tags in the same file) and
``..._port_table.csv``.
"""

from __future__ import annotations

import csv
import sys
import time
from pathlib import Path

import numpy as np
import ufl
from mpi4py import MPI

from dolfinx import default_scalar_type, fem
from dolfinx import mesh as dmesh

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.symmetry import full_impedance  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_th16_symmetry_plane_driven import (  # noqa: E402
    ALL_FACES,
    EXPECTED,
    EXPECTED_PMC_CUT_FACETS,
    LEG_D0_REPRODUCTION_BAND,
    NC1_FLOOR,
    POINTS,
    S_X,
    S_Y,
    SHEET_TAG,
    SHEET_Z,
    X1,
    X_CUT,
    Y_CUT,
    _distribute,
    _mx,
    _my,
    _owned_count,
    _solve,
    _tag,
    build_cases,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_12_symmetry_plane_quarter_box"

#: Gate record, 20260929T111548Z_TH-16-step2-W2.log:96-132 — printed beside
#: this run, a difference reported and never asserted.
RECORD = {
    "Z_f": 1.082882778204e01 + 3.482120962988e00j,
    "P_f": 1.458517771552e-03,
    "mirror": 1.038e-14,
    "nc1": 5.786521e-01,
}

#: Solve order: the full model first, so every identity has its reference.
ORDER = ("full", "pmc_half", "pec_half", "quarter")


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _rel(a, b):
    return abs(a / b - 1.0)


def _quarter_cell_tags(msh):
    """Figure-only cell tags on the full box: 1 = kept quarter, 2 = mirrored away."""
    n = msh.topology.index_map(3).size_local + msh.topology.index_map(3).num_ghosts
    cells = np.arange(n, dtype=np.int32)
    mid = dmesh.compute_midpoints(msh, 3, cells)
    kept = (mid[:, 0] < X_CUT) & (mid[:, 1] < Y_CUT)
    return dmesh.meshtags(msh, 3, cells, np.where(kept, 1, 2).astype(np.int32))


def _export(name, msh, tags, fields, comm):
    """CG1 interpolants of Re E and Im E (interpolants: nothing asserted on them)."""
    e = fields.e_complex
    v_vec = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    e_real = fem.Function(v_vec, name="E_real")
    e_real.interpolate(fem.Expression(ufl.real(e), v_vec.element.interpolation_points))
    e_imag = fem.Function(v_vec, name="E_imag")
    e_imag.interpolate(fem.Expression(ufl.imag(e), v_vec.element.interpolation_points))
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_{name}_combined", msh, None,
        {"E_real": e_real, "E_imag": e_imag}, comm=comm, facet_tags=tags,
    )
    return path


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "This example needs the complex DolfinX build: "
            "source /usr/local/bin/dolfinx-complex-mode (the runner does this "
            "automatically for the `th:` group)."
        )

    _say(comm, "=" * 78)
    _say(comm, "EX-66 -- quarter the box: a PMC plane and a PEC plane, the port rescaled")
    _say(comm, "=" * 78)
    _say(comm, f"\n[gate] band {LEG_D0_REPRODUCTION_BAND:.0e} on every identity, control floor "
               f"{NC1_FLOOR:.0%} (imported from the TH-16 step 2 gate module)")

    cases, _ = build_cases()
    full_points = np.ascontiguousarray(
        np.vstack([POINTS, _mx(POINTS), _my(POINTS), _mx(_my(POINTS))])
    )
    res, outputs = {}, []
    nc1 = None
    for name in ORDER:
        x, c, extent, spec, pec_tags = cases[name]
        msh = _distribute(x, c, comm)
        tags = _tag(msh, extent)
        nv = msh.topology.index_map(0).size_global
        nc = msh.topology.index_map(3).size_global
        ns = _owned_count(msh, tags, SHEET_TAG)
        _say(comm, f"\n[{name:8s}] vertices {nv}, cells {nc}, owned sheet facets {ns} "
                   f"(ASSERTED == {EXPECTED[name]}); pinned faces {pec_tags}")
        assert (nv, nc, ns) == EXPECTED[name], (name, nv, nc, ns)
        if name in EXPECTED_PMC_CUT_FACETS:
            npmc = _owned_count(msh, tags, X1)
            _say(comm, f"[{name:8s}] natural (PMC) cut facets on x = {X_CUT} : {npmc} "
                       f"(ASSERTED == {EXPECTED_PMC_CUT_FACETS[name]})")
            assert npmc == EXPECTED_PMC_CUT_FACETS[name], (name, npmc)

        if name == "full":
            # ---- setup figure (EX-57): the full box, the kept quarter tagged ----
            write_setup_figure(
                msh, _quarter_cell_tags(msh), FIGURE_DIR / f"{BASENAME}_setup.png",
                region_names={1: "kept quarter (saline)", 2: "mirrored away (saline)"},
                translucent_tags=(1, 2),
                slice_normal=(0.0, 0.0, 1.0),
                slice_origin=(X_CUT, Y_CUT, SHEET_Z),
                title="th:12 -- saline box 0.48x0.40x0.48 m, sheet at z = 0.24, "
                      "PMC cut x = 0.24, PEC cut y = 0.20; slice z = 0.24",
                comm=comm,
            )

        t = time.perf_counter()
        pts = full_points if name == "full" else POINTS
        r = _solve(msh, tags, spec, pec_tags, pts, return_fields=True)
        res[name] = r
        _say(comm, f"[{name:8s}] V={r['V']:.12e} I={r['I']:.12e} Z={r['Z']:.12e} "
                   f"P={r['P']:.12e} valid={int(np.sum(r['valid']))}/{len(pts)} "
                   f"({time.perf_counter() - t:.1f}s)")
        assert r["valid"].all(), name
        outputs.append(_export(name, msh, tags, r["fields"], comm))
        if name == "pmc_half":
            nc1 = _solve(msh, tags, spec, ALL_FACES, POINTS)

    f = res["full"]
    _say(comm, f"\n[record] Z_f {RECORD['Z_f']:.12e} Ohm, P_f {RECORD['P_f']:.12e} W "
               f"(20260929T111548Z_TH-16-step2-W2.log:110); this run minus record: "
               f"Z {f['Z'] - RECORD['Z_f']:+.2e}, P {f['P'] - RECORD['P_f']:+.2e} "
               "(reported, not asserted)")

    # (1) mirror symmetry of the full field
    n = POINTS.shape[0]
    e = f["E"]
    e0, ex, ey, exy = e[:n], e[n:2 * n], e[2 * n:3 * n], e[3 * n:]
    scale = np.max(np.linalg.norm(e0, axis=1))
    mirror = max(
        np.max(np.linalg.norm(ex - S_X * e0, axis=1)),
        np.max(np.linalg.norm(ey - S_Y * e0, axis=1)),
        np.max(np.linalg.norm(exy - S_X * S_Y * e0, axis=1)),
    ) / scale
    _say(comm, f"\n(1) full-field mirror residual {mirror:.3e} (record {RECORD['mirror']:.3e}; "
               f"ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e})")
    assert mirror <= LEG_D0_REPRODUCTION_BAND

    # (2) PMC half, (3) PEC half, (4) quarter
    h, p, q = res["pmc_half"], res["pec_half"], res["quarter"]
    e_rel = np.max(np.linalg.norm(h["E"] - e0, axis=1)) / scale
    identities = {
        "(2) PMC half": {
            "Z=full_impedance(Z_h,'pmc')": _rel(full_impedance(h["Z"], "pmc"), f["Z"]),
            "V": _rel(h["V"], f["V"]), "I/2": _rel(h["I"], f["I"] / 2.0),
            "P/2": _rel(h["P"], f["P"] / 2.0), "E": e_rel,
        },
        "(3) PEC half": {
            "Z=full_impedance(Z_h,'pec')": _rel(full_impedance(p["Z"], "pec"), f["Z"]),
            "V/2": _rel(p["V"], f["V"] / 2.0), "I": _rel(p["I"], f["I"]),
            "P/2": _rel(p["P"], f["P"] / 2.0),
        },
        "(4) quarter": {"Z": _rel(q["Z"], f["Z"]), "P/4": _rel(q["P"], f["P"] / 4.0)},
    }
    for label, errs in identities.items():
        worst = max(errs.values())
        body = ", ".join(f"{k} {v:.2e}" for k, v in errs.items())
        _say(comm, f"{label}: rel errors {body}; worst {worst:.2e} "
                   f"(ASSERTED <= {LEG_D0_REPRODUCTION_BAND:.0e})")
        for k, v in errs.items():
            assert v <= LEG_D0_REPRODUCTION_BAND, (label, k, v)

    # ---- negative control: the PMC cut pinned (the wrong condition) ----------
    miss = abs(nc1["Z"] / (2.0 * f["Z"]) - 1.0)
    _say(comm, f"\n[control] PMC half, cut pinned: |Z/(2Z_f) - 1| = {miss:.6e} "
               f"(ASSERTED > {NC1_FLOOR}); gate record {RECORD['nc1']:.6e}")
    assert miss > NC1_FLOOR, miss

    # ---- port table -----------------------------------------------------------
    if comm.rank == 0:
        OUTPUT_DIR.mkdir(exist_ok=True)
        table = OUTPUT_DIR / f"{BASENAME}_port_table.csv"
        with table.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["model", "V_re", "V_im", "I_re", "I_im", "Z_re", "Z_im", "P_W",
                        "Z_full_re", "Z_full_im"])
            restore = {"full": lambda z: z, "pmc_half": lambda z: full_impedance(z, "pmc"),
                       "pec_half": lambda z: full_impedance(z, "pec"), "quarter": lambda z: z}
            for name in ORDER:
                r = res[name]
                zf = restore[name](r["Z"])
                w.writerow([name, f"{r['V'].real:.12e}", f"{r['V'].imag:.12e}",
                            f"{r['I'].real:.12e}", f"{r['I'].imag:.12e}",
                            f"{r['Z'].real:.12e}", f"{r['Z'].imag:.12e}", f"{r['P']:.12e}",
                            f"{zf.real:.12e}", f"{zf.imag:.12e}"])
        outputs.append(table)

    for path in outputs:
        _say(comm, f"[output] {path}")
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
