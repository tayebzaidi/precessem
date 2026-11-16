"""Example (`EX-153`): a graded, PEC-backed PML at 60 degrees -- TE10 guide, set against the ABC.

The `th:18` guide (a = lambda0 / (2 sin 60 deg), TE10 drive at x = 0, PEC y and z
walls, physical length L_O = 1.5 m) is extended by a PML slab L_O < x < L_O + PML_D
(Jin 3e 9.6.2 (9.152)-(9.154), FE form (9.157), profile (9.165), m = 2,
s_max = S_ABSORB, PEC at the far end).  The metal-backed PML reflects
|Gamma_c(theta)| = R(0)^cos(theta) (Jin 3e 9.6.1 (9.133) integrated over the
profile), so at 60 degrees, R(0) = 1e-6: Gamma_c = +1e-3.  The first-order ABC at
the same angle reflects -1/3.

**It asserts, and it does not re-implement** (the `ANS-1` rule).  Every constant,
band and helper comes by name from `tests/validation/test_th5_pml_oblique.py`
(as `t3b`), `test_th5_pml_reflection.py` and `test_th5_abc_reflection.py`; the
terms from `fem_em_solver.core.boundary_terms`.  No gate-module or `src/` edit.
Asserted (the gate's asserts): Gamma_c = 1e-3 within 1e-12; the angular identity
<= IDENTITY_TOLERANCE; k_fit scan minimum interior; (B) degree 2
|Gamma - Gamma_c| <= PML_REFLECTION_BAR; (Cabc) |Gamma_ABC|/|Gamma_PML| >
CABC_RATIO_BAR; (NC-0) s_max = 0 reads |Gamma| > NC0_BAR.  Records are printed
beside the readings and never asserted for equality.

Run it::

    ./run_examples.sh -e th:22 -n 2 -t 420

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_22_pml_oblique_incidence_pml_combined.xdmf`` and
``..._abc_combined.xdmf`` (CG1 ``E_real``, ``E_imag``; PML slab cell tag 2 and
facet tags in the PML file).
"""

from __future__ import annotations

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

from fem_em_solver.core.boundary_terms import (  # noqa: E402
    pml_reflection_pec_backed,
    pml_slab_term,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

import tests.validation.test_th5_pml_oblique as t3b  # noqa: E402
from tests.validation.test_port_reaction_impedance import IDENTITY_TOLERANCE  # noqa: E402
from tests.validation.test_th5_abc_reflection import (  # noqa: E402
    A_O,
    BETA_O,
    FIT_O,
    K0,
    L_O,
    LAMBDA0,
    OMEGA,
    TAG_Y,
    TAG_Z,
    _box,
    _Projection,
    _solve,
)
from tests.validation.test_th5_pml_reflection import (  # noqa: E402
    PML_D,
    PML_M,
    PML_REFLECTION_BAR,
    R_ABSORB,
    S_ABSORB,
    _solve_terms,
)

COS_O = t3b.COS_O
N_GATE = t3b.N_GATE
CABC_RATIO_BAR = t3b.CABC_RATIO_BAR
NC0_BAR = t3b.NC0_BAR
_drive = t3b._drive

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_22_pml_oblique_incidence"

#: Gate records -- printed beside this run, a difference reported and never asserted.
#: W1 = 20261025T111709Z_TH-5-step3b-W1.log
RECORD = {
    "ident": 1.110e-15,         # W1:52
    "err_deg2": 1.1031e-07,     # W1:53
    "err_k0cos_deg2": 2.8512e-07,  # W1:53
    "ratio": 333.3,             # W1:59
    "nc0": 1.000000,            # W1:60
    "gamma_abc_th18": -0.33325,  # th:18 guide
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _read(e_h, msh, t0, ncells):
    x = ufl.SpatialCoordinate(msh)
    proj = _Projection(e_h[2], msh, FIT_O[0], FIT_O[1], BETA_O,
                       w=ufl.sin(ufl.pi * x[1] / A_O))
    out = {"cells": ncells, "proj": proj.gamma(BETA_O, L_O)}
    out["k_fit"], out["gamma_k_fit"], out["scan_interior"] = proj.k_fit(BETA_O, L_O)
    out["seconds"] = time.perf_counter() - t0
    return out


def _guide_fields(comm, s_max):
    """`t3b._guide` written out (:76-96) keeping the field.  s_max None: the ABC guide."""
    t0 = time.perf_counter()
    h = LAMBDA0 / N_GATE
    ny = int(A_O / h)
    wz = 2.0 * h
    if s_max is None:
        msh, tags, ncells = _box(comm, L_O, A_O, wz, int(round(L_O / h)), ny, 2)
        e_h = _solve(msh, tags, 2, _drive, [TAG_Y, TAG_Z], "abc")
    else:
        length = L_O + PML_D
        msh, tags, ncells = _box(comm, length, A_O, wz, int(round(length / h)), ny, 2)
        term = pml_slab_term(msh, L_O, PML_D, s_max, PML_M, OMEGA)
        e_h = _solve_terms(msh, tags, 2, _drive, [TAG_Y, TAG_Z], [term])
    out = _read(e_h, msh, t0, ncells)
    return out, e_h, msh, tags


def _region_tags(msh):
    """Cell tags: 1 guide (physical), 2 the PML slab x > L_O."""
    tdim = msh.topology.dim
    n_owned = msh.topology.index_map(tdim).size_local
    idx = np.arange(n_owned, dtype=np.int32)
    mids = dmesh.compute_midpoints(msh, tdim, idx)
    vals = np.ones(n_owned, dtype=np.int32)
    vals[mids[:, 0] > L_O] = 2
    return dmesh.meshtags(msh, tdim, idx, vals)


def _export(msh, ct, ft, e_h, comm, suffix):
    v_vec = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    e_real = fem.Function(v_vec, name="E_real")
    e_real.interpolate(fem.Expression(ufl.real(e_h), v_vec.element.interpolation_points))
    e_imag = fem.Function(v_vec, name="E_imag")
    e_imag.interpolate(fem.Expression(ufl.imag(e_h), v_vec.element.interpolation_points))
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_{suffix}_combined", msh, ct,
        {"E_real": e_real, "E_imag": e_imag}, comm=comm, facet_tags=ft,
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
    _say(comm, "EX-153 -- a graded PEC-backed PML at 60 degrees, against the ABC")
    _say(comm, "=" * 78)
    gc = pml_reflection_pec_backed(S_ABSORB, PML_D, PML_M, K0, cos_theta=COS_O)
    gc0 = pml_reflection_pec_backed(S_ABSORB, PML_D, PML_M, K0)
    ident = abs(abs(gc) / abs(gc0) ** COS_O - 1.0)
    _say(comm, f"\n[gate] PML_D {PML_D}, PML_M {PML_M}, S_ABSORB {S_ABSORB:.6f}, "
               f"R(0) {R_ABSORB}, bar {PML_REFLECTION_BAR}, lambda0/{N_GATE}, degree 2")
    _say(comm, f"[closed form] cos = {COS_O:.15f}; Gamma_c(60) = {gc:.6e}, |Gamma_c| = "
               f"{abs(gc):.6e} (R(0)^cos = {R_ABSORB ** COS_O:.6e}); "
               f"| |Gc|/|Gc0|^cos - 1 | = {ident:.3e} (record {RECORD['ident']:.3e})")

    pml, e_pml, msh, ft = _guide_fields(comm, S_ABSORB)
    _say(comm, f"[PML] degree 2: cells {pml['cells']}, {pml['seconds']:.1f} s; k_fit - beta = "
               f"{pml['k_fit'] - BETA_O:+.4e} rad/m (interior {pml['scan_interior']}); "
               f"k_fit read G = {pml['gamma_k_fit']:.6e}, |G - Gc| = "
               f"{abs(pml['gamma_k_fit'] - gc):.4e} (record {RECORD['err_deg2']:.4e}); "
               f"k0cos read |G - Gc| = {abs(pml['proj'] - gc):.4e} "
               f"(record {RECORD['err_k0cos_deg2']:.4e})")

    ct = _region_tags(msh)
    write_setup_figure(
        msh, ct, FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={1: "guide interior (PEC y, z walls; drive at x = 0)",
                      2: "PML slab 1.5 < x < 3.0 m (PEC at x = 3.0 m)"},
        translucent_tags=(1,),
        slice_normal=(0.0, 0.0, 1.0),
        slice_origin=(0.0, 0.0, 0.5 * float(2.0 * LAMBDA0 / N_GATE)),
        title="th:22 -- TE10 guide a = 0.577 m (60 deg), graded PML 1.5 lambda0 "
              "behind x = 1.5 m; slice z = mid",
        comm=comm,
    )
    path_pml = _export(msh, ct, ft, e_pml, comm, "pml")
    _say(comm, f"[output] {path_pml}")
    del e_pml, msh, ft, ct

    abc, e_abc, msh_a, ft_a = _guide_fields(comm, None)
    path_abc = _export(msh_a, None, ft_a, e_abc, comm, "abc")
    _say(comm, f"[ABC] degree 2: cells {abc['cells']}, {abc['seconds']:.1f} s; k_fit read G = "
               f"{abc['gamma_k_fit']:.6e} (interior {abc['scan_interior']}); "
               f"k0cos read G = {abc['proj']:.6e} (th:18: {RECORD['gamma_abc_th18']}, "
               f"closed form -1/3)")
    _say(comm, f"[output] {path_abc}")
    del e_abc, msh_a, ft_a

    nc0 = t3b._guide(comm, N_GATE, 2, 0.0)
    _say(comm, f"[NC-0] s_max = 0: {nc0['seconds']:.1f} s; k_fit read |G| = "
               f"{abs(nc0['gamma_k_fit']):.6f}, k0cos read |G| = {abs(nc0['proj']):.6f} "
               f"(record {RECORD['nc0']:.6f})")

    err = abs(pml["gamma_k_fit"] - gc)
    ratio = abs(abc["gamma_k_fit"]) / abs(pml["gamma_k_fit"])
    _say(comm, "")
    _say(comm, f"(Gc) |Gc - 1e-3| = {abs(gc - 1e-3):.3e} <= 1e-12: ASSERTED")
    _say(comm, f"(ident) {ident:.3e} <= IDENTITY_TOLERANCE {IDENTITY_TOLERANCE}: ASSERTED")
    _say(comm, f"(interior) k_fit scan minimum interior = {pml['scan_interior']}: ASSERTED")
    _say(comm, f"(B) |G - Gc| = {err:.4e} <= PML_REFLECTION_BAR {PML_REFLECTION_BAR:.1e}: ASSERTED")
    _say(comm, f"(Cabc) |G_ABC|/|G_PML| = {ratio:.1f} (record {RECORD['ratio']}) > "
               f"{CABC_RATIO_BAR}: ASSERTED")
    _say(comm, f"(NC-0) |G| = {abs(nc0['gamma_k_fit']):.6f} > {NC0_BAR}: ASSERTED")

    assert abs(gc - 1e-3) <= 1e-3 * 1e-9, gc
    assert ident <= IDENTITY_TOLERANCE, ident
    assert pml["scan_interior"], pml["k_fit"] - BETA_O
    assert ratio > CABC_RATIO_BAR, ratio
    assert abs(nc0["gamma_k_fit"]) > NC0_BAR, nc0["gamma_k_fit"]
    assert err <= PML_REFLECTION_BAR, err

    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
