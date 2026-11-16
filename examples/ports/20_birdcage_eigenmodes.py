"""Example (`EX-76`): the tuned birdcage's own modes in ParaView.

`TH-17` (closed) solved the generalised eigenproblem (G1) of the gapped 4-leg
birdcage on `TH-15` step 3's **PEC-hole** mesh, with the four lumped-sheet
gaps carrying a series capacitor ``C_mode`` chosen so the k = 1 circuit mode
sits on 64 MHz. The gate module only ever reduced the answer to frequencies;
this example keeps the eigenvectors and shows ``|E|`` of both k = 1 partners
and the k = 2 singlet above them.

**It asserts, and it does not re-implement** (the `ANS-1` rule). Every
helper, band and record is imported from
``tests/validation/test_th17_birdcage_eigenmodes.py`` (the fixture, the (G1)
solve, the driven sweep that sets ``C_mode``, the frequency map) and from the
modules that gate imports its bands from: ``KAPPA``
(``STEP2D_C_OVER_TERMINAL_64MHZ_P1``) and ``ADJACENT_SPREAD_BAND``.

**Rule (c), `src/`, disclosed.** ``solve_general_mesh_modes`` returned no
eigenvectors. An additive ``return_vectors`` keyword (default ``False``, the
unset path byte-identical; `EX-60`'s precedent) passes through to
``_solve_pencil_nonhermitian``. The test helper ``solve_modes`` got the same
keyword; it appends one dict ``{"vectors", "V", "eps", "forms", "k0_sq"}``.
The gate's step-3 test was re-run green in the same slot.

**Asserted (imported bands).** (a) ``|Re f - 64 MHz| / 64 MHz <= KAPPA`` for
both k = 1 members; (b) their split ``<= ADJACENT_SPREAD_BAND`` and
``f_pair < f2`` (the k = 2 mode); (d') ``|(f_pair(2C)/f_pair(C)) sqrt2 - 1|
<= KAPPA``; the negative control — at ``2 C_mode`` the pair misses 64 MHz by
more than ``KAPPA``; and the mode identity ``||A x - lambda B x|| /
||lambda B x|| <= 1e-6`` for each exported vector, with ``(A, B)``
re-assembled by ``general_mesh_pencil`` at ``bc_diagonal = 1e4 k0^2``.

**Printed, not asserted.** ``C_mode`` beside the `TH-17` W3 record and the
pair / split / ``f2`` beside ``GATE_E_STEP3_PAIR_HZ`` and the W3 log.

**Scope.** ``C_mode``, degree 1. Not gate (c)'s circuit scan, not eigenvector
class labels (the k = 1 partners are an arbitrary basis of a near-degenerate
pair, their polarisation is not labelled), no AED comparison, no B1, no Q
claim. The hole mesh has no conductor cell tag.

**Negative result protocol.** A red imported band here while the gate module
holds it is an example/test divergence: a known-issues entry and a stop,
never a re-record from the example side.

Needs the complex DolfinX build; the runner sources it for ``ports:``::

    ./run_examples.sh -e ports:20 -n 4 -t 300
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
from mpi4py import MPI

from dolfinx import default_scalar_type, fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.core.cavity import general_mesh_pencil  # noqa: E402
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.test_birdcage_conductor_hole import AIR_CELL_TAG  # noqa: E402
from tests.mesh.test_birdcage_port_sheets import (  # noqa: E402
    PORT_LOWER,
    PORT_UPPER,
)
from tests.mesh.test_birdcage_port_tags import LEG_COUNT  # noqa: E402
from tests.validation.test_port_birdcage_lumped_column import (  # noqa: E402
    ADJACENT_SPREAD_BAND,
    PHANTOM_CELL_TAG,
    SALINE_SIGMA,
)
from tests.validation.test_port_lumped_rlc_termination import (  # noqa: E402
    STEP2D_C_OVER_TERMINAL_64MHZ_P1 as KAPPA,
)
from tests.validation.test_th17_birdcage_eigenmodes import (  # noqa: E402
    GATE_E_STEP3_PAIR_HZ,
    TARGET_FREQUENCY_HZ,
    _driven_hole_z,
    _modal_z,
    _physical_freqs,
    build_hole_fixture,
    frequency_of,
    solve_modes,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
BASENAME = "ports_20_birdcage_eigenmodes"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
SETUP_BASENAME = "ports_20_birdcage_eigenmodes_setup"

#: Mode-identity band (the brief's, a residual of the solver's own pencil).
RESIDUAL_BAND = 1.0e-6
#: Printed-only W3 records (`20261001T124845Z_TH-17-gate-c-W3.log`).
W3_C_MODE_F = 6.628769508e-11
W3_SPLIT = 5.573e-05
W3_F2_HZ = 7.838759746e+07
W3_NEG_MISS = 2.8653e-01
BC_DIAGONAL_FACTOR = 1.0e4  # `solve_general_mesh_modes` default, as the gate runs it


def _residual(A, B, lam, vec):
    """``||A x - lam B x|| / ||lam B x||`` — the residual form, not ``x.dot(y)``."""
    ax = A.createVecLeft()
    bx = A.createVecLeft()
    A.mult(vec, ax)
    B.mult(vec, bx)
    ax.axpy(-lam, bx)
    return float(ax.norm() / (abs(lam) * bx.norm()))


def _abs_e_dg0(V, vec, name, comm):
    """Unit-peak ``|E|`` as a real DG0 field from an eigenvector.

    Interpolated into ``("DG", 0, (3,))``; the modulus is taken in numpy and
    normalised by the global maximum (eigenvector phase and scale are arbitrary).
    """
    msh = V.mesh
    u = fem.Function(V)
    n_loc = V.dofmap.index_map.size_local * V.dofmap.index_map_bs
    u.x.array[:n_loc] = vec.array_r
    u.x.scatter_forward()
    W = fem.functionspace(msh, ("DG", 0, (3,)))
    w = fem.Function(W)
    w.interpolate(fem.Expression(u, W.element.interpolation_points))
    arr = w.x.array.reshape(-1, 3)
    mag = np.sqrt(np.sum(np.abs(arr) ** 2, axis=1))
    n_cells = msh.topology.index_map(msh.topology.dim).size_local
    peak = comm.allreduce(float(mag[:n_cells].max()) if n_cells else 0.0, op=MPI.MAX)
    Q = fem.functionspace(msh, ("DG", 0))
    out = fem.Function(Q, name=name)
    out.x.array[:] = mag / peak
    out.x.scatter_forward()
    return out


def main() -> None:
    comm = MPI.COMM_WORLD
    started = time.perf_counter()

    if not np.issubdtype(np.dtype(default_scalar_type), np.complexfloating):
        raise RuntimeError(
            "This example needs the complex DolfinX build: "
            "source /usr/local/bin/dolfinx-complex-mode (the runner does this "
            "automatically for the `ports:` group)."
        )

    def say(msg):
        if comm.rank == 0:
            print(msg, flush=True)

    say("=" * 78)
    say("EX-76 -- the tuned birdcage's own modes (TH-17 hole mesh, C_mode)")
    say("=" * 78)

    f64 = TARGET_FREQUENCY_HZ
    omega_64 = 2.0 * np.pi * f64
    fx = build_hole_fixture()
    say(f"\n[fixture] {fx['n_cells']} cells, mesh {fx['t_mesh_s']:.2f} s, "
        f"kappa = {KAPPA:.6e}, spread band = {ADJACENT_SPREAD_BAND}")

    # EX-57 setup figure, right after the mesh is built: air hidden, phantom
    # translucent, slice z = 0 (the port-gap mid-plane), the `ports:14`
    # pattern on the same hole mesh (no conductor cell tag exists on it).
    write_setup_figure(
        fx["mesh"],
        fx["cell_tags"],
        FIGURE_DIR / f"{SETUP_BASENAME}.png",
        region_names={
            AIR_CELL_TAG: "air",
            PHANTOM_CELL_TAG: "phantom",
            **{PORT_LOWER + i: f"leg L{i} port lower" for i in range(1, LEG_COUNT + 1)},
            **{PORT_UPPER + i: f"leg L{i} port upper" for i in range(1, LEG_COUNT + 1)},
        },
        hide_tags=(AIR_CELL_TAG,),
        translucent_tags=(PHANTOM_CELL_TAG,),
        slice_normal=(0.0, 0.0, 1.0),
        title="ports:20 -- the TH-17 hole mesh, four sheets at C_mode (no conductor cell tag)",
        comm=comm,
    )

    # --- C_mode from one 64 MHz sweep (gate c's lines, composed from imports) ---
    t0 = time.perf_counter()
    z = _driven_hole_z(fx, f64, sigma_s_per_m=SALINE_SIGMA)
    zk64, off64 = _modal_z(z)
    z1 = complex(0.5 * (zk64[1] + zk64[3]))
    assert z1.imag > 0.0, f"Im z1(64 MHz) = {z1.imag:.6e} <= 0: no C_mode"
    c_mode = 1.0 / (omega_64 * z1.imag)
    say(f"[C_mode] 64 MHz sweep ({time.perf_counter() - t0:.2f} s): off-diag {off64:.2e}; "
        f"z1 = {z1:.9e}; C_mode = {c_mode:.9e} F (W3 record {W3_C_MODE_F:.9e} F, "
        f"rel diff {abs(c_mode - W3_C_MODE_F) / W3_C_MODE_F:.3e}; printed)")

    # --- eigen solve at C_mode, eigenvectors kept ------------------------------
    values, n_conv, n_dofs, wall_e, _nbc, ex = solve_modes(
        fx, c_mode, omega_lin_rad_s=omega_64, target_hz=f64, return_vectors=True)
    V, eps, forms, k0_sq, vectors = ex["V"], ex["eps"], ex["forms"], ex["k0_sq"], ex["vectors"]
    assert len(vectors) == len(values)
    # physical eigenvalues (Re f above the gradient cluster), sorted by Re f,
    # keeping the index into the solver's order
    phys_idx = sorted(
        (i for i, lam in enumerate(values) if frequency_of(lam).real > 1.0e6),
        key=lambda i: frequency_of(values[i]).real,
    )
    phys = _physical_freqs(values)
    assert len(phys_idx) == len(phys) and len(phys) >= 3, (
        "need the k = 1 pair and the k = 2 singlet")
    say(f"[eigen] at C_mode (target 64 MHz): dofs {n_dofs}, converged {n_conv}, {wall_e:.2f} s")
    for f in phys:
        say(f"    f = {f.real:.9e} {f.imag:+.9e}j Hz")
    pair = [phys[0].real, phys[1].real]
    f2 = phys[2].real
    f_pair = 0.5 * (pair[0] + pair[1])
    errs = [abs(p - f64) / f64 for p in pair]
    split = abs(pair[1] - pair[0]) / f_pair
    say(f"[(a)] |Re f - 64|/64 = {errs[0]:.4e} / {errs[1]:.4e} "
        f"({errs[0] / KAPPA:.3f} / {errs[1] / KAPPA:.3f} kappa); pair {pair[0]:.9e} / "
        f"{pair[1]:.9e} Hz beside GATE_E_STEP3_PAIR_HZ {GATE_E_STEP3_PAIR_HZ[0]:.9e} / "
        f"{GATE_E_STEP3_PAIR_HZ[1]:.9e} (rel diff "
        f"{abs(pair[0] - GATE_E_STEP3_PAIR_HZ[0]) / GATE_E_STEP3_PAIR_HZ[0]:.2e} / "
        f"{abs(pair[1] - GATE_E_STEP3_PAIR_HZ[1]) / GATE_E_STEP3_PAIR_HZ[1]:.2e}; printed)")
    say(f"[(b)] split {split:.3e} (W3 {W3_SPLIT:.3e}; band {ADJACENT_SPREAD_BAND}); "
        f"f2 = {f2:.9e} Hz (W3 {W3_F2_HZ:.9e}, rel diff "
        f"{abs(f2 - W3_F2_HZ) / W3_F2_HZ:.2e}; printed)")

    # --- mode identity: residuals on the re-assembled pencil -------------------
    A, B, _n = general_mesh_pencil(
        V, eps, bc_diagonal=BC_DIAGONAL_FACTOR * k0_sq, sheet_mass_forms=forms)
    names = ("E_mode1a", "E_mode1b", "E_mode2")
    fields = {}
    residuals = {}
    for name, i in zip(names, phys_idx[:3]):
        residuals[name] = _residual(A, B, complex(values[i]), vectors[i])
        fields[name] = _abs_e_dg0(V, vectors[i], name, comm)
    A.destroy()
    B.destroy()
    say("[identity] ||Ax - lam Bx|| / ||lam Bx|| = "
        + ", ".join(f"{n} {r:.3e}" for n, r in residuals.items())
        + f" (band {RESIDUAL_BAND:.0e})")

    # --- negative control + (d'): the pencil at 2 C_mode -----------------------
    values_n, _nn, _dn, wall_n, _b = solve_modes(
        fx, 2.0 * c_mode, omega_lin_rad_s=omega_64, target_hz=f64 / np.sqrt(2.0))
    phys_n = _physical_freqs(values_n)
    assert len(phys_n) >= 2, "2*C_mode: fewer than two physical eigenvalues"
    for f in phys_n:
        say(f"    2C f = {f.real:.9e} {f.imag:+.9e}j Hz")
    pair_n = [phys_n[0].real, phys_n[1].real]
    f_pair_n = 0.5 * (pair_n[0] + pair_n[1])
    miss = min(abs(p - f64) / f64 for p in pair_n)
    d_prime = abs((f_pair_n / f_pair) * np.sqrt(2.0) - 1.0)
    say(f"[(d')] {d_prime:.4e} ({d_prime / KAPPA:.3f} kappa; pair(2C) {f_pair_n:.9e}, "
        f"pair(C) {f_pair:.9e}; {wall_n:.2f} s)")
    say(f"[NEG] 2*C_mode pair misses 64 MHz by {miss:.4e} (= {miss / KAPPA:.2f} kappa; "
        f"W3 record {W3_NEG_MISS:.4e} = 27.01 kappa)")

    # --- the asserted bands (imported) ----------------------------------------
    for p, e in zip(pair, errs):
        assert e <= KAPPA, f"(a) {p:.9e} Hz: |Re f - 64|/64 = {e:.4e} > kappa {KAPPA:.4e}"
    assert split <= ADJACENT_SPREAD_BAND, f"(b) pair split {split:.3e}"
    assert f_pair < f2, f"(b) f_pair {f_pair:.9e} >= f2 {f2:.9e}"
    assert d_prime <= KAPPA, f"(d') {d_prime:.4e} > kappa {KAPPA:.4e}"
    assert miss > KAPPA, f"negative control: 2*C_mode pair within kappa ({miss:.3e})"
    for n, r in residuals.items():
        assert r <= RESIDUAL_BAND, f"mode identity {n}: {r:.3e} > {RESIDUAL_BAND:.0e}"

    # --- combined XDMF ---------------------------------------------------------
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        fx["mesh"],
        fx["cell_tags"],
        fields,
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    say(f"\n[output] {path} -- CellTags, DG0 E_mode1a / E_mode1b / E_mode2 (unit peak |E|)")
    say(f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    say("\nAll imported bands hold.")


if __name__ == "__main__":
    main()
