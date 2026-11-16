"""Figures that equal their records (`WF-8` step 1, §9 item 328, chain U22).

Three figure kinds, each drawn from a written record with no transform of the
values it plots, so every plotted array read back from the figure object is
bitwise equal (``np.array_equal``, float64) to the source array it draws:

* :func:`plot_sparameters` — ``|S_ij|`` and ``arg S_ij`` (degrees) against
  frequency from a Touchstone file, the arrays being exactly
  ``skrf.Network(path).f`` / ``.s_mag[:, i, j]`` / ``.s_deg[:, i, j]``.
* :func:`plot_xdmf_slice` — a cell-centred field from a DolfinX XDMF/HDF5
  bundle on a slab ``|c_axis - c0| <= h`` of cell centroids.  The slice is a
  **selection, not an interpolation**: the plotted colours are
  ``values[indices]`` and the plotted points are the centroids of ``indices``.
* :func:`plot_probe_series` — one CSV column against ``t_s``.

Rank 0 reads and draws; every other rank returns ``None``; all ranks pass a
``Barrier``.  Style choices are minimal defaults (the operator's to tune).
Matplotlib's ``Agg`` backend is selected before ``pyplot`` is imported.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
from mpi4py import MPI

__all__ = [
    "FigureResult",
    "read_sparameters",
    "read_xdmf_cell_field",
    "read_probe_series",
    "select_slab",
    "draw_sparameters",
    "draw_slice",
    "draw_series",
    "plot_sparameters",
    "plot_xdmf_slice",
    "plot_probe_series",
    "mismatch_indices",
]

_AXES = {"x": 0, "y": 1, "z": 2}


def _pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


@dataclass
class FigureResult:
    """A drawn figure plus handles to the artists that carry the data."""

    figure: object
    path: Optional[Path] = None
    lines: Dict[object, object] = field(default_factory=dict)
    collection: object = None
    indices: Optional[np.ndarray] = None


def mismatch_indices(plotted, source) -> np.ndarray:
    """Flat indices where ``plotted`` and ``source`` differ (exact comparison).

    A shape mismatch returns ``[-1]`` so it can never pass as "no mismatch".
    """
    a = np.asarray(plotted)
    b = np.asarray(source)
    if a.shape != b.shape:
        return np.array([-1])
    return np.flatnonzero(a.ravel() != b.ravel())


# --------------------------------------------------------------------------- readers


def read_sparameters(path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(f, s_mag, s_deg)`` exactly as ``skrf.Network(path)`` gives them."""
    import skrf

    nw = skrf.Network(str(path))
    return nw.f, nw.s_mag, nw.s_deg


def _xdmf_items(xdmf) -> Tuple[Path, Dict[str, str]]:
    text = Path(xdmf).read_text()
    items = {}
    h5 = None
    for m in re.finditer(r">\s*([^<>:\s]+\.h5):(/[^<\s]+)\s*<", text):
        h5 = Path(xdmf).parent / m.group(1)
        items[m.group(2)] = m.group(2)
    if h5 is None:
        raise ValueError(f"no HDF5 data item in {xdmf}")
    return h5, items


def read_xdmf_cell_field(xdmf, name) -> Tuple[np.ndarray, np.ndarray]:
    """``(centroids, values)`` for a cell-centred field ``name`` in a DolfinX bundle.

    ``values`` is ``/Function/<name>/0`` raveled (float64, no copy-transform);
    ``centroids`` is the vertex mean of each cell of ``/Mesh/mesh/topology``.
    """
    import h5py

    h5, _ = _xdmf_items(xdmf)
    with h5py.File(h5, "r") as f:
        geometry = np.asarray(f["/Mesh/mesh/geometry"][...], dtype=np.float64)
        topology = np.asarray(f["/Mesh/mesh/topology"][...], dtype=np.int64)
        values = np.asarray(f[f"/Function/{name}/0"][...]).ravel()
    if values.shape[0] != topology.shape[0]:
        raise ValueError(f"{name!r} is not cell-centred ({values.shape[0]} values, {topology.shape[0]} cells)")
    centroids = geometry[topology].mean(axis=1)
    return centroids, values


def read_probe_series(csv, column) -> Tuple[np.ndarray, np.ndarray]:
    """``(t_s, column)`` parsed from a probe-series CSV with ``np.genfromtxt``."""
    data = np.genfromtxt(str(csv), delimiter=",", names=True, dtype=np.float64)
    return np.asarray(data["t_s"]), np.asarray(data[column])


def select_slab(centroids, axis="z", center=None, half_width=None) -> np.ndarray:
    """Indices of cells whose centroid lies in ``|c_axis - center| <= half_width``.

    Defaults: ``center`` is the mid-plane of the centroid range on ``axis``;
    ``half_width`` is 5 % of that range.
    """
    c = np.asarray(centroids)[:, _AXES[axis]]
    lo, hi = float(c.min()), float(c.max())
    if center is None:
        center = 0.5 * (lo + hi)
    if half_width is None:
        half_width = 0.05 * (hi - lo)
    return np.flatnonzero(np.abs(c - center) <= half_width)


# --------------------------------------------------------------------------- drawers


def draw_sparameters(f, s_mag, s_deg) -> FigureResult:
    plt = _pyplot()
    n = s_mag.shape[1]
    fig, (ax_m, ax_p) = plt.subplots(2, 1, sharex=True, figsize=(7, 7))
    lines = {}
    for i in range(n):
        for j in range(n):
            label = f"S{i + 1}{j + 1}"
            (lines[(i, j, "mag")],) = ax_m.plot(f, s_mag[:, i, j], marker="o", label=label)
            (lines[(i, j, "deg")],) = ax_p.plot(f, s_deg[:, i, j], marker="o", label=label)
    ax_m.set_ylabel("|S| (linear)")
    ax_p.set_ylabel("arg S (deg)")
    ax_p.set_xlabel("frequency (Hz)")
    ax_m.legend(fontsize=6, ncol=4)
    return FigureResult(figure=fig, lines=lines)


def draw_slice(centroids, values, indices, axis="z", label="") -> FigureResult:
    plt = _pyplot()
    a, b = [k for k in range(3) if k != _AXES[axis]]
    fig, ax = plt.subplots(figsize=(6, 5))
    coll = ax.scatter(centroids[indices, a], centroids[indices, b], c=values[indices], s=8)
    fig.colorbar(coll, ax=ax, label=label)
    names = "xyz"
    ax.set_xlabel(names[a])
    ax.set_ylabel(names[b])
    ax.set_aspect("equal")
    return FigureResult(figure=fig, collection=coll, indices=np.asarray(indices))


def draw_series(t, y, label="") -> FigureResult:
    plt = _pyplot()
    fig, ax = plt.subplots(figsize=(6, 4))
    (line,) = ax.plot(t, y, label=label)
    ax.set_xlabel("t (s)")
    ax.set_ylabel(label)
    return FigureResult(figure=fig, lines={label: line})


# --------------------------------------------------------------------------- public


def _save(result: FigureResult, out) -> FigureResult:
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.figure.savefig(out)
    result.path = out
    try:
        from fem_em_solver.io.paraview_utils import adopt_host_ownership

        adopt_host_ownership(out.parent, comm=MPI.COMM_SELF)
    except Exception:  # ownership is a papercut fix, never a failure
        pass
    return result


def _rank0(build, comm):
    comm = comm or MPI.COMM_WORLD
    result = build() if comm.rank == 0 else None
    comm.Barrier()
    return result


def plot_sparameters(path, out, comm=None) -> Optional[FigureResult]:
    """Draw ``|S_ij|`` / ``arg S_ij`` against frequency from a Touchstone file."""
    return _rank0(lambda: _save(draw_sparameters(*read_sparameters(path)), out), comm)


def plot_xdmf_slice(xdmf, name, out, axis="z", center=None, half_width=None, comm=None) -> Optional[FigureResult]:
    """Draw cell-centred ``name`` on a centroid slab normal to ``axis``."""

    def build():
        centroids, values = read_xdmf_cell_field(xdmf, name)
        idx = select_slab(centroids, axis=axis, center=center, half_width=half_width)
        return _save(draw_slice(centroids, values, idx, axis=axis, label=name), out)

    return _rank0(build, comm)


def plot_probe_series(csv, column, out, comm=None) -> Optional[FigureResult]:
    """Draw CSV ``column`` against ``t_s``."""
    return _rank0(lambda: _save(draw_series(*read_probe_series(csv, column), label=column), out), comm)
