"""Example (`EX-101`): a sheet drawn inside a box becomes an integration surface.

For the user who will draw a port as a literal 2D sheet in an external CAD tool
(Phase 9 (9b), "ports as literal 2D sheets"): the fixture is `GEO-37`'s unit box
with one 0.8 x 0.6 rectangle embedded by ``occ.fragment`` and tagged as a dim-2
physical group (201); the box faces are tagged 101. The model is read through
the project shim ``_model_to_mesh`` with an explicit ``shared_facet`` ghosted
cell partitioner, as `GEO-37` step 2's probe does, and the sheet is then
*integrated over* with ``dS(201)``.

Everything gated is **imported** from the `GEO-37` modules (the `ANS-1` rule):
``STRIP_BAND``, the closed forms ``I1_CF`` .. ``I4_CF``, the fixture builder,
the quadrature pin and the assemble helper. The only literals here are the
planted fields' coefficients (g = 2 + 3x + 5y, DG0 side values 1 / 3), which
are inputs, written exactly as the probe writes them
(``probe_geo37_interior_dS.py:114, 120``) and encoded in ``I2_CF`` / ``I3_CF``.

Asserted: (G) the mesh is ghosted before any ``dS`` is assembled; (T) the
facet-tag set is ``[101, 201]`` and the owned tag-201 count equals gmsh's own
triangle count; (A) area, (F) the linear and two-sided piecewise-constant
integrals, (I) the exterior ``ds(101)`` instrument, all within ``STRIP_BAND``;
(V) the tag validator raises on an absent tag. Controls: (c1) ``dS(999)``
assembles to exactly 0.0 (asserted), (c2) one planted dropped facet leaves the
band (asserted), (c3) a one-sided ``u('+')`` integral is printed, never
asserted.

**Not a port, not a STEP import, no unghosted run.** The face is drawn in gmsh
in-process; ``PortDefinition`` is used only as a tag-presence check.

Run it through the example runner::

    ./run_examples.sh -e mesh:13

Output: ``examples/meshing/paraview_output/meshing_13_interior_sheet_surface_integral_combined.xdmf``.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import gmsh
import numpy as np
import ufl
from mpi4py import MPI

import dolfinx
from dolfinx import fem

_REPO_ROOT = Path(__file__).resolve().parents[2]
# The runner puts only ``src`` on PYTHONPATH; the repo root goes on sys.path so
# the gate's constants and helpers can be imported rather than restated.
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from fem_em_solver.io.mesh import (  # noqa: E402
    _model_to_mesh,
    _raise_geometry_failure_on_every_rank,
)
from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.ports.definitions import (  # noqa: E402
    PortDefinition,
    validate_required_port_tags_exist,
)
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.mesh.probe_geo37_interior_dS import QD, TAG_ABSENT, _assemble  # noqa: E402
from tests.mesh.probe_geo37_interior_facet import (  # noqa: E402
    H,
    RECT,
    TAG_BOX,
    TAG_RECT,
    _build_model,
)
from tests.mesh.test_geo37_interior_dS import I1_CF, I2_CF, I3_CF, I4_CF  # noqa: E402
from tests.validation.test_ops61_exact_strip_ports import STRIP_BAND  # noqa: E402

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"  # committed, EX-57
BASENAME = "meshing_13_interior_sheet_surface_integral"


def _rel(value, closed_form):
    return value / closed_form - 1.0


def _mask_ok(comm, ok, what):
    """Collective check: raise on every rank if any rank's ``ok`` is False."""
    if not comm.allreduce(bool(ok), op=MPI.MIN):
        raise AssertionError(what)


def main() -> None:
    comm = MPI.COMM_WORLD
    rank = comm.rank
    started = time.perf_counter()

    if rank == 0:
        print("=" * 72)
        print("EX-101 — an interior sheet becomes an integration surface (dS)")
        print("=" * 72)
        print(
            f"\n[fixture] GEO-37 unit box, embedded rectangle RECT={RECT}, "
            f"h={H}; facet tags: sheet={TAG_RECT}, box faces={TAG_BOX}"
            f"\n[gate]    STRIP_BAND={STRIP_BAND:.3e}; ghost mode shared_facet; "
            f"ranks={comm.size}",
            flush=True,
        )

    # ---- build (rank 0) and import (collective), ghosted --------------------
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    t0 = time.perf_counter()
    err, n_gmsh = None, None
    if rank == 0:
        try:
            n_gmsh = _build_model()
        except BaseException as e:  # noqa: BLE001 -- GEO-23 pattern
            err = e
            gmsh.finalize()
    _raise_geometry_failure_on_every_rank(comm, 0, "EX-101 sheet fixture", err, H)
    n_gmsh = comm.bcast(n_gmsh, root=0)
    t1 = time.perf_counter()

    # 0.11: max_facet_to_cell_links is required; 2 for a tet manifold mesh.
    partitioner = dolfinx.mesh.create_cell_partitioner(
        dolfinx.mesh.GhostMode.shared_facet, 2
    )
    mesh, cell_tags, facet_tags = _model_to_mesh(
        gmsh.model, comm, 0, gdim=3, partitioner=partitioner
    )
    gmsh.finalize()
    comm.Barrier()
    t2 = time.perf_counter()

    tdim = mesh.topology.dim
    fdim = tdim - 1
    imap_c = mesh.topology.index_map(tdim)
    n_owned_c = imap_c.size_local
    n_cells = comm.allreduce(n_owned_c, op=MPI.SUM)
    n_ghost_cells = comm.allreduce(imap_c.num_ghosts, op=MPI.SUM)
    mesh.topology.create_connectivity(fdim, tdim)
    n_owned_f = mesh.topology.index_map(fdim).size_local

    # ---- (G) ghosted before any dS is assembled -----------------------------
    ghost_total = comm.allreduce(imap_c.num_ghosts, op=MPI.SUM)
    ghosted = comm.size == 1 or ghost_total > 0
    _mask_ok(comm, ghosted, "(G) the mesh carries no ghost cells at >1 rank: "
             "dS would be assembled unghosted (the PORT-1 3b-iv hang)")

    # ---- (T) the sheet arrives as a tag, complete and unduplicated ----------
    tag_set = sorted({int(t) for part in comm.allgather(
        sorted({int(v) for v in facet_tags.values})) for t in part})
    local201 = facet_tags.find(TAG_RECT)
    owned201 = local201[local201 < n_owned_f]
    count_201 = int(comm.allreduce(len(owned201), op=MPI.SUM))
    unfiltered_201 = int(comm.allreduce(len(local201), op=MPI.SUM))
    per_rank_owned = comm.gather(int(len(owned201)), root=0)
    per_rank_local = comm.gather(int(len(local201)), root=0)
    f2c = mesh.topology.connectivity(fdim, tdim)
    one_sided = sum(1 for f in local201 if len(f2c.links(int(f))) < 2)
    per_rank_one_sided = comm.gather(int(one_sided), root=0)

    # ---- planted fields and measures ----------------------------------------
    dS = ufl.Measure("dS", domain=mesh, subdomain_data=facet_tags, metadata=QD)
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags, metadata=QD)
    x = ufl.SpatialCoordinate(mesh)
    g = 2 + 3 * x[0] + 5 * x[1]  # probe_geo37_interior_dS.py:114

    V0 = fem.functionspace(mesh, ("DG", 0))
    u = fem.Function(V0, name="u_side")
    n_all_cells = imap_c.size_local + imap_c.num_ghosts
    mids = dolfinx.mesh.compute_midpoints(
        mesh, tdim, np.arange(n_all_cells, dtype=np.int32)
    )
    # probe_geo37_interior_dS.py:120 -- DG0 u = 1 below / 3 above the sheet plane
    u.x.array[:n_all_cells] = np.where(mids[:, 2] > RECT[2], 3.0, 1.0)
    u.x.scatter_forward()
    one = fem.Constant(mesh, dolfinx.default_scalar_type(1.0))

    # ---- the integrals (every one allreduced inside _assemble) --------------
    t3 = time.perf_counter()
    I1, e1 = _assemble(comm, one * dS(TAG_RECT))
    I2, e2 = _assemble(comm, g("+") * dS(TAG_RECT))
    I3, e3 = _assemble(comm, (u("+") + u("-")) * dS(TAG_RECT))
    I3p, e3p = _assemble(comm, u("+") * dS(TAG_RECT))
    I4, e4 = _assemble(comm, one * ds(TAG_BOX))
    c1, ec1 = _assemble(comm, one * dS(TAG_ABSENT))

    # (c2) one planted dropped facet: lowest rank owning any 201 facet drops its
    # first owned one; the retag is built on every rank (sorted, unique int32).
    owners = comm.allgather(len(owned201) > 0)
    drop_rank = owners.index(True) if True in owners else -1
    keep = np.ones(len(facet_tags.indices), dtype=bool)
    if rank == drop_rank:
        first = int(owned201[0])
        keep &= ~((facet_tags.indices == first) & (facet_tags.values == TAG_RECT))
    ft_drop = dolfinx.mesh.meshtags(
        mesh,
        fdim,
        np.asarray(facet_tags.indices[keep], dtype=np.int32),
        np.asarray(facet_tags.values[keep], dtype=np.int32),
    )
    dS_drop = ufl.Measure("dS", domain=mesh, subdomain_data=ft_drop, metadata=QD)
    I1_drop, ed = _assemble(comm, one * dS_drop(TAG_RECT))
    t4 = time.perf_counter()

    # ---- (V) the validator, caught on every rank ----------------------------
    validate_required_port_tags_exist(
        [PortDefinition("P1", TAG_RECT, TAG_BOX, "native_interior_sheet")], tag_set
    )
    raised_text = None
    try:
        validate_required_port_tags_exist(
            [PortDefinition("P1", TAG_ABSENT, TAG_BOX, "absent_tag_control")], tag_set
        )
    except ValueError as e:
        raised_text = str(e)

    # ---- XDMF: planted fields + rank + facet tags ----------------------------
    V1 = fem.functionspace(mesh, ("Lagrange", 1))
    g_planted = fem.Function(V1, name="g_planted")
    g_planted.interpolate(lambda p: 2 + 3 * p[0] + 5 * p[1])
    rank_f = fem.Function(V0, name="rank")
    rank_f.x.array[:] = 0.0
    rank_f.x.array[:n_owned_c] = float(rank)
    rank_f.x.scatter_forward()
    OUTPUT_DIR.mkdir(exist_ok=True)
    xdmf_path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        mesh,
        cell_tags,
        {"g_planted": g_planted, "u_side": u, "rank": rank_f},
        comm=comm,
        facet_tags=facet_tags,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)

    # ---- setup figure (EX-57; display-only cell tags, never in an integral) --
    c_idx = mesh.topology.connectivity(fdim, tdim)
    touching = set()
    for f in local201:
        for c in c_idx.links(int(f)):
            if int(c) < n_owned_c:
                touching.add(int(c))
    disp_vals = np.ones(n_owned_c, dtype=np.int32)
    if touching:
        disp_vals[np.array(sorted(touching), dtype=np.int64)] = 2
    display_tags = dolfinx.mesh.meshtags(
        mesh, tdim, np.arange(n_owned_c, dtype=np.int32), disp_vals
    )
    write_setup_figure(
        mesh,
        display_tags,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={
            1: "air box (faces: facet tag 101)",
            2: "tets touching the drawn sheet (facet tag 201)",
        },
        translucent_tags=(1,),
        slice_normal=(0.0, 1.0, 0.0),
        title=("mesh:13 -- a sheet drawn inside a box becomes facet tag 201; "
               "tets touching it shown (GEO-37)"),
        comm=comm,
    )

    # ---- (c3) values and prints ---------------------------------------------
    total = time.perf_counter() - started
    if rank == 0:
        print(f"\n[mesh] cells={n_cells} cells_ghost={n_ghost_cells} "
              f"N_gmsh={n_gmsh} (records: 10150 cells, 192 triangles)")
        print(f"[mesh] seconds: gmsh {t1 - t0:.2f}  import {t2 - t1:.2f}  "
              f"assemble {t4 - t3:.2f}")
        print(f"[G] ghost cells (reduced) = {ghost_total} (>0 required at >1 rank)")
        print(f"[T] reduced facet-tag set = {tag_set}")
        print(f"[T] owned tag-{TAG_RECT} count = {count_201} vs N_gmsh = {n_gmsh}")
        print(f"[T] per-rank owned / local / one-sided tag-{TAG_RECT}: "
              f"{per_rank_owned} / {per_rank_local} / {per_rank_one_sided}")
        print(f"[T] unfiltered rank-local tag-{TAG_RECT} sum = {unfiltered_201} "
              f"beside owned sum = {count_201} (printed, not asserted)")
        print(f"[A] I1 = {I1!r}  rel {_rel(I1, I1_CF):+.3e}  (CF {I1_CF})")
        print(f"[F] I2 = {I2!r}  rel {_rel(I2, I2_CF):+.3e}  (CF {I2_CF})")
        print(f"[F] I3 = {I3!r}  rel {_rel(I3, I3_CF):+.3e}  (CF {I3_CF})")
        print(f"[I] I4 = {I4!r}  rel {_rel(I4, I4_CF):+.3e}  (CF {I4_CF})")
        print(f"[V] absent-tag validator raised: {raised_text!r}")
        print(f"[c1] dS({TAG_ABSENT}) = {c1!r} (asserted exactly 0.0)")
        print(f"[c2] drop_rank={drop_rank} I1' = {I1_drop!r} rel "
              f"{_rel(I1_drop, I1_CF):+.3e}; one-facet share 1/N_gmsh = "
              f"{1.0 / n_gmsh:.3e}; band {STRIP_BAND:.3e}")
        print(f"[c3] one-sided int u('+') dS({TAG_RECT}) = {I3p!r} "
              f"(window [{I1_CF:.2f}, {3 * I1_CF:.2f}]; differs between widths; "
              "printed, not asserted)")
        print(f"\n[paraview] wrote {xdmf_path}")
        print(f"[time] total {total:.1f} s", flush=True)

    # ---- asserts, identical on every rank ------------------------------------
    for name, err_text in (("I1", e1), ("I2", e2), ("I3", e3), ("I4", e4),
                           ("c1", ec1), ("c2", ed)):
        assert err_text is None, f"{name} assembly raised: {err_text}"
    assert tag_set == sorted([TAG_BOX, TAG_RECT]), f"(T) tag set {tag_set}"
    assert count_201 == n_gmsh, (
        f"(T) owned tag-{TAG_RECT} count {count_201} != gmsh triangles {n_gmsh}"
    )
    for name, val, cf in (("A I1", I1, I1_CF), ("F I2", I2, I2_CF),
                          ("F I3", I3, I3_CF), ("I I4", I4, I4_CF)):
        assert abs(_rel(val, cf)) <= STRIP_BAND, (
            f"({name}) {val!r} misses its closed form {cf} by "
            f"{_rel(val, cf):+.3e} (band {STRIP_BAND:.3e})"
        )
    assert raised_text is not None and "missing required port tags: [999]" in raised_text, (
        f"(V) validator did not raise on the absent tag: {raised_text!r}"
    )
    assert c1 == 0.0, f"(c1) dS({TAG_ABSENT}) = {c1!r}, expected exactly 0.0"
    assert abs(_rel(I1_drop, I1_CF)) > STRIP_BAND, (
        f"(c2) the planted dropped facet stayed inside the band: {I1_drop!r}"
    )

    if rank == 0:
        print("\nAll identities and controls hold.", flush=True)


if __name__ == "__main__":
    main()
