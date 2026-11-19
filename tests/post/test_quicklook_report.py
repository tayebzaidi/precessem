"""Tests for quick-look phantom metrics reporting helpers (chunk E3)."""

from __future__ import annotations

from pathlib import Path

from fem_em_solver.post.quicklook import (
    build_phantom_quicklook_report,
    format_phantom_quicklook_report,
    write_phantom_quicklook_report,
)
from tests.complex_mode import complex_only


def _sample_metrics() -> dict:
    return {
        "sampling": {
            "requested_cells": 20,
            "sampling_cells": 16,
            "valid_sample_cells": 16,
            "boundary_adjacent_cells_dropped": 4,
            "invalid_samples_dropped_e": 0,
            "invalid_samples_dropped_b": 0,
        },
        "E_magnitude": {
            "count": 16,
            "min": 1.0,
            "max": 3.0,
            "mean": 2.0,
        },
        "B_magnitude": {
            "count": 16,
            "min": 0.2,
            "max": 0.8,
            "mean": 0.4,
        },
        "consistency": {
            "e_to_b_mean_ratio": 5.0,
            "e_to_b_max_ratio": 6.5,
            "mean_balance_rel_diff": 0.12,
            "warnings": [],
        },
    }


def test_build_phantom_quicklook_report_ok_status_without_warning_flags():
    report = build_phantom_quicklook_report(_sample_metrics())

    assert report["status"] == "OK"
    assert report["warning_flags"] == []
    assert report["consistency_warnings"] == []
    assert report["summary"]["e_mean"] == 2.0
    assert report["summary"]["b_mean"] == 0.4


def test_build_phantom_quicklook_report_warn_status_when_consistency_warnings_present():
    metrics = _sample_metrics()
    metrics["consistency"]["warnings"] = ["Large E/B imbalance detected"]

    report = build_phantom_quicklook_report(metrics)

    assert report["status"] == "WARN"
    assert report["warning_flags"] == []
    assert report["consistency_warnings"] == ["Large E/B imbalance detected"]


def test_quicklook_report_passes_an_unreported_balance_through(tmp_path: Path):
    """`mean_balance_rel_diff` is ``None`` without a reference speed: OK, no flag, n/a."""
    metrics = _sample_metrics()
    metrics["consistency"]["mean_balance_rel_diff"] = None

    report = build_phantom_quicklook_report(metrics)

    assert report["status"] == "OK"
    assert report["warning_flags"] == []
    assert report["summary"]["mean_balance_rel_diff"] is None
    text = format_phantom_quicklook_report(report)
    assert "mean-balance relative diff: n/a" in text

    outputs = write_phantom_quicklook_report(report, output_dir=tmp_path, basename="no_speed")
    json_text = Path(outputs["json"]).read_text(encoding="utf-8")
    assert '"mean_balance_rel_diff": null' in json_text
    assert "NaN" not in json_text


def test_format_and_write_phantom_quicklook_report_produce_expected_artifacts(tmp_path: Path):
    report = build_phantom_quicklook_report(_sample_metrics())
    text = format_phantom_quicklook_report(report)

    assert "quick-look phantom metrics:" in text
    assert "status: OK" in text
    assert "warning flags: none" in text

    outputs = write_phantom_quicklook_report(
        report,
        output_dir=tmp_path,
        basename="mri_coil_phantom",
        write_json=True,
        write_markdown=True,
    )

    json_path = Path(outputs["json"])
    md_path = Path(outputs["markdown"])

    assert json_path.exists()
    assert md_path.exists()

    md_text = md_path.read_text(encoding="utf-8")
    assert "quick-look phantom metrics:" in md_text
    assert "status: OK" in md_text


# ---------------------------------------------------------------------------
# `WF-3` (§9 item 99): the report generated from a field with a closed form.
# `TH-6`'s lossy plane wave E = ẑ e^{−jkx}: on the cell centroids x_c,
# |E| = e^{−αx_c} (α = −Im k) and |B| = (|k|/ω)|E|, so every number the report
# prints is known.  Parsed from the *written markdown*, so a min/max or E/B
# swap in the formatter is visible.  Never reads the mean-balance line.
# ---------------------------------------------------------------------------


def _parse_quicklook_markdown(text: str) -> dict:
    """The eight numbers of the three value lines, keyed like the JSON summary."""
    out = {}
    for line in text.splitlines():
        label, _, rest = line.strip().partition(": ")
        values = [float(tok) for tok in rest.split(" / ")] if " / " in rest else []
        if label == "|E| min/max/mean":
            out["e_min"], out["e_max"], out["e_mean"] = values
        elif label == "|B| min/max/mean":
            out["b_min"], out["b_max"], out["b_mean"] = values
        elif label == "shape ratios, non-physical (mean/max)":
            out["e_to_b_mean_ratio"], out["e_to_b_max_ratio"] = values
        elif label == "status":
            out["status"] = rest
    return out


@complex_only
def test_wf3_quicklook_report_prints_the_plane_wave_closed_form():
    import copy
    import json
    import tempfile

    import numpy as np
    from dolfinx import mesh as dmesh
    from mpi4py import MPI

    from fem_em_solver.post.faraday import magnetic_flux_density_from_e
    from fem_em_solver.post.phantom_fields import compute_phantom_eb_metrics_and_export
    from fem_em_solver.utils.constants import EPSILON_0, MU_0
    from tests.validation.test_field_consistency_metrics import POST2_BAND
    from tests.validation.test_lossy_plane_wave import (
        EPSILON_R,
        FREQUENCY_HZ,
        SIGMA,
        _solve_plane_wave,
        _wavenumber,
    )

    comm = MPI.COMM_WORLD
    n = 24
    _rel, _a, _b, ncells, msh, fields = _solve_plane_wave(n, SIGMA, return_fields=True)
    omega = 2.0 * np.pi * FREQUENCY_HZ
    k = _wavenumber(SIGMA)
    alpha = -k.imag
    b_fn = magnetic_flux_density_from_e(fields.e_complex, omega)

    tdim = msh.topology.dim
    cmap = msh.topology.index_map(tdim)
    n_local = cmap.size_local + cmap.num_ghosts
    tags = dmesh.meshtags(
        msh, tdim, np.arange(n_local, dtype=np.int32), np.ones(n_local, dtype=np.int32)
    )

    root_dir = tempfile.mkdtemp(prefix="wf3_quicklook_") if comm.rank == 0 else None
    output_dir = Path(comm.bcast(root_dir, root=0))
    metrics = compute_phantom_eb_metrics_and_export(
        fields.e_complex,
        b_fn,
        tags,
        phantom_tag=1,
        output_dir=output_dir,
        basename="wf3",
        comm=comm,
    )
    report = build_phantom_quicklook_report(metrics)

    # Rank 0 writes and re-reads both artifacts; every rank asserts on the broadcast.
    parsed = loaded = None
    if comm.rank == 0:
        outputs = write_phantom_quicklook_report(report, output_dir=output_dir, basename="wf3")
        md_text = Path(outputs["markdown"]).read_text(encoding="utf-8")
        parsed = _parse_quicklook_markdown(md_text)
        loaded = json.loads(Path(outputs["json"]).read_text(encoding="utf-8"))
    parsed = comm.bcast(parsed, root=0)
    loaded = comm.bcast(loaded, root=0)

    # Closed form on the owned cell centroids (numpy; extremes MIN/MAX-allreduced).
    owned = np.arange(cmap.size_local, dtype=np.int32)
    cx = (
        np.mean(msh.geometry.x[msh.geometry.dofmap[owned], 0], axis=1)
        if owned.size
        else np.array([])
    )
    x_min = comm.allreduce(float(cx.min()) if cx.size else np.inf, op=MPI.MIN)
    x_max = comm.allreduce(float(cx.max()) if cx.size else -np.inf, op=MPI.MAX)
    count = comm.allreduce(int(cx.size), op=MPI.SUM)
    mean_e = comm.allreduce(float(np.exp(-alpha * cx).sum()), op=MPI.SUM) / count
    b_over_e = abs(k) / omega
    ratio_cf = omega / abs(k)
    closed = {
        "e_min": np.exp(-alpha * x_max),
        "e_max": np.exp(-alpha * x_min),
        "e_mean": mean_e,
        "b_min": b_over_e * np.exp(-alpha * x_max),
        "b_max": b_over_e * np.exp(-alpha * x_min),
        "b_mean": b_over_e * mean_e,
        "e_to_b_mean_ratio": ratio_cf,
        "e_to_b_max_ratio": ratio_cf,
    }
    miss = {key: parsed[key] / closed[key] - 1.0 for key in closed}
    ratio_blind = 1.0 / np.sqrt(MU_0 * EPSILON_0) / np.sqrt(EPSILON_R)
    n1 = parsed["e_to_b_mean_ratio"] / ratio_blind - 1.0
    n2 = parsed["e_min"] / closed["e_max"] - 1.0
    n3 = parsed["b_mean"] / closed["e_mean"] - 1.0
    summary = loaded["summary"]
    r1 = {key: abs(parsed[key] / summary[key] - 1.0) for key in closed}

    if comm.rank == 0:
        print(
            f"\nWF-3 n={n}: cells={ncells} centroids={count} alpha={alpha:.6f} "
            f"x_c in [{x_min:.6f}, {x_max:.6f}] omega/|k|={ratio_cf:.6e}"
        )
        for key in closed:
            print(
                f"  {key}: markdown {parsed[key]:.6e}  closed form {closed[key]:.6e}  "
                f"miss {miss[key] * 100:+.4f} %  |md/json - 1| {r1[key]:.2e}"
            )
        print(
            f"  status line (unmodified): {parsed['status']}; warning_flags={report['warning_flags']}; "
            f"consistency_warnings={report['consistency_warnings']}"
        )
        print(
            f"  controls: n1 mean ratio vs c/sqrt(78) {n1 * 100:+.4f} % (prior -21.1963 %); "
            f"n2 |E| min vs max's closed form {n2 * 100:+.2f} % (predicted ~72 %); "
            f"n3 |B| mean vs |E|'s {n3 * 100:+.4f} % (predicted ~100 %)"
        )

    # (A) every parsed number on the closed form, in POST-2's band.
    for key, value in miss.items():
        assert abs(value) < POST2_BAND, (key, value)
    # (R1) the markdown is the JSON at .6e.
    for key, value in r1.items():
        assert value <= 5e-7, (key, value)
    # (R2) the JSON summary is the metrics dict, exactly; sampling is 6·24³, nothing dropped.
    assert summary["e_min"] == metrics["E_magnitude"]["min"]
    assert summary["e_max"] == metrics["E_magnitude"]["max"]
    assert summary["e_mean"] == metrics["E_magnitude"]["mean"]
    assert summary["b_min"] == metrics["B_magnitude"]["min"]
    assert summary["b_max"] == metrics["B_magnitude"]["max"]
    assert summary["b_mean"] == metrics["B_magnitude"]["mean"]
    assert summary["e_to_b_mean_ratio"] == metrics["consistency"]["e_to_b_mean_ratio"]
    assert summary["e_to_b_max_ratio"] == metrics["consistency"]["e_to_b_max_ratio"]
    sampling = loaded["sampling"]
    assert (
        sampling["requested_cells"],
        sampling["sampling_cells"],
        sampling["valid_sample_cells"],
        sampling["boundary_adjacent_cells_dropped"],
        sampling["invalid_samples_dropped_e"],
        sampling["invalid_samples_dropped_b"],
    ) == (6 * n**3, 6 * n**3, 6 * n**3, 0, 0, 0)
    assert count == 6 * n**3
    # (R3) no shape flag; the consistency warnings pass through; emptied, the status is OK.
    assert loaded["warning_flags"] == []
    assert loaded["consistency_warnings"] == list(metrics["consistency"]["warnings"])
    quiet = copy.deepcopy(metrics)
    quiet["consistency"]["warnings"] = []
    quiet_report = build_phantom_quicklook_report(quiet)
    assert quiet_report["status"] == "OK"
    assert "status: OK" in format_phantom_quicklook_report(quiet_report)
    # Negative controls.
    assert abs(abs(n1) - 0.212) <= 0.001, n1
    assert abs(n2) > POST2_BAND, n2
    assert abs(n3) > POST2_BAND, n3
