"""`PORT-7` — the Touchstone writer, read by an independent reader (scikit-rf).

``tests/io/test_touchstone_export.py`` round-trips a 2-port through our writer
*and our own reader*; both used the same flattening order, so a layout error
cancels and N >= 3 was never exercised. Touchstone 1.x lists a 2-port as
S11 S21 S12 S22 and an N >= 3 port matrix row by row (S11 S12 ... S1N on the
first line(s), then S21 ...). This module judges our writer and reader against
scikit-rf 2.1.0 (pinned, `OPS-60`):

* **Premise (asserted first):** spec-layout 2- and 3-port literals, one matrix
  row per line, read by ``skrf.Network`` equal the known matrix to 1e-15. If
  red, skrf cannot judge and nothing below means anything.
* **Anchor (asserted):** writer -> ``skrf.Network`` and
  ``Network.write_touchstone`` -> ``load_touchstone``, both max|dS| <= 1e-12,
  on (a) the stored `PORT-15` 4x4 records (10 and 64 MHz, one stack), (b) the
  stored 3x3 C-terminated record, (c) seeded non-reciprocal matrices
  N = 2, 3, 4 with |S_ij| <= 0.5; metadata ``nports``, ``f`` to 1e-3 Hz,
  ``z0`` = 50 ohm.
* **Negative control:** (c) as read vs its own transpose -- *asserted* the
  premise max|S - S^T| >= 0.1 on (c), then a miss >= 1e-2 (ceiling 1).
  *Predicted, printed:* (a)'s max|S - S^T| is at round-off -- the reciprocal
  field records cannot see a layout error; only (c) can.

Scope: RI, Hz, one real z0, Touchstone 1.x.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
import skrf
from mpi4py import MPI

from fem_em_solver.ports import PortDefinition
from fem_em_solver.ports.sparameters import (
    SParameterSweepResult,
    summarize_sparameter_sanity,
)
from fem_em_solver.ports.touchstone import export_touchstone, load_touchstone
from tests.validation.test_port_circuit_layer_field import (
    S_10MHZ_EPS0_RECORD,
    S_64MHZ_EPS0_RECORD,
    S_64MHZ_TERMINATED_C_RECORD,
)

PREMISE_TOL = 1.0e-15
ANCHOR_TOL = 1.0e-12
FREQ_ATOL_HZ = 1.0e-3
Z0_OHM = 50.0
TRANSPOSE_PREMISE_FLOOR = 0.1
TRANSPOSE_MISS_FLOOR = 1.0e-2
TRANSPOSE_MISS_CEILING = 1.0
SEED = 20260929


def _rank_dir(tmp_path: Path) -> Path:
    d = tmp_path / f"rank{MPI.COMM_WORLD.rank}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _seeded_nonreciprocal(n_ports: int, n_freq: int) -> np.ndarray:
    rng = np.random.default_rng(SEED + n_ports)
    mag = 0.5 * rng.uniform(0.05, 1.0, size=(n_freq, n_ports, n_ports))
    phase = rng.uniform(-np.pi, np.pi, size=(n_freq, n_ports, n_ports))
    return mag * np.exp(1j * phase)


def _cases() -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """name -> (frequencies_hz, S stack of shape (n_freq, N, N))."""
    cases = {
        "a_records_4x4": (
            np.array([10.0e6, 64.0e6]),
            np.stack([S_10MHZ_EPS0_RECORD, S_64MHZ_EPS0_RECORD]),
        ),
        "b_record_3x3_C": (np.array([64.0e6]), S_64MHZ_TERMINATED_C_RECORD[None, :, :]),
    }
    for n in (2, 3, 4):
        cases[f"c_nonrecip_{n}port"] = (
            np.array([63.87e6, 127.74e6]),
            _seeded_nonreciprocal(n, 2),
        )
    return cases


CASES = _cases()


def _write_ours(freqs: np.ndarray, s_stack: np.ndarray, out_dir: Path) -> Path:
    n = s_stack.shape[1]
    port_ids = tuple(f"P{i + 1}" for i in range(n))
    ports = [
        PortDefinition(
            port_id=pid,
            positive_tag=10 * (i + 1) + 1,
            negative_tag=10 * (i + 1) + 2,
            orientation="cw",
            z0_ohm=Z0_OHM,
        )
        for i, pid in enumerate(port_ids)
    ]
    sweeps = [
        SParameterSweepResult(
            frequency_hz=float(f),
            port_ids=port_ids,
            s_matrix=np.asarray(s, dtype=np.complex128),
            excitation_results={},
            sanity_report=summarize_sparameter_sanity(np.asarray(s, dtype=np.complex128)),
            is_placeholder=False,
        )
        for f, s in zip(freqs, s_stack)
    ]
    return export_touchstone(sweeps, ports, output_dir=out_dir, write_csv_companion=False)


def _max_abs(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))


# --- premise ----------------------------------------------------------------

_S2_LITERAL = np.array(
    [[0.11 + 0.12j, 0.21 + 0.22j], [0.31 + 0.32j, 0.41 + 0.42j]], dtype=np.complex128
)
_S3_LITERAL = np.array(
    [
        [0.11 + 0.12j, 0.13 + 0.14j, 0.15 + 0.16j],
        [0.21 + 0.22j, 0.23 + 0.24j, 0.25 + 0.26j],
        [0.31 + 0.32j, 0.33 + 0.34j, 0.35 + 0.36j],
    ],
    dtype=np.complex128,
)


def _spec_literal(s: np.ndarray, freq_hz: float) -> str:
    """Touchstone 1.x by the spec: 2-port S11 S21 S12 S22 on one line;
    N >= 3 one matrix row per line (N <= 4 here, so no continuation)."""
    n = s.shape[0]
    lines = ["! spec-layout literal", f"# Hz S RI R {Z0_OHM:g}"]
    if n == 2:
        vals = [s[0, 0], s[1, 0], s[0, 1], s[1, 1]]
        lines.append(
            f"{float(freq_hz)!r} " + " ".join(f"{float(v.real)!r} {float(v.imag)!r}" for v in vals)
        )
    else:
        for i in range(n):
            row = " ".join(f"{float(v.real)!r} {float(v.imag)!r}" for v in s[i, :])
            lines.append((f"{float(freq_hz)!r} " if i == 0 else "") + row)
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize("s_literal", [_S2_LITERAL, _S3_LITERAL], ids=["2port", "3port"])
def test_premise_skrf_reads_spec_layout_literals(tmp_path, s_literal):
    n = s_literal.shape[0]
    path = _rank_dir(tmp_path) / f"literal.s{n}p"
    path.write_text(_spec_literal(s_literal, 64.0e6), encoding="utf-8")
    ntw = skrf.Network(str(path))
    err = _max_abs(ntw.s[0], s_literal)
    print(f"\n[PORT-7 premise] N={n} skrf max|dS|={err:.3e} (tol {PREMISE_TOL:g})")
    assert ntw.nports == n
    assert err <= PREMISE_TOL


# --- anchor -----------------------------------------------------------------


@pytest.mark.parametrize("case", list(CASES), ids=list(CASES))
def test_our_writer_read_by_skrf(tmp_path, case):
    freqs, s_stack = CASES[case]
    n = s_stack.shape[1]
    path = _write_ours(freqs, s_stack, _rank_dir(tmp_path))
    ntw = skrf.Network(str(path))
    err = _max_abs(ntw.s, s_stack)
    err_t = _max_abs(ntw.s, np.swapaxes(s_stack, 1, 2))
    print(
        f"\n[PORT-7 writer->skrf] {case}: N={n} max|dS|={err:.3e} "
        f"max|S_read - S^T|={err_t:.3e} nports={ntw.nports} "
        f"f={ntw.f.tolist()} z0={np.unique(ntw.z0).tolist()}"
    )
    assert ntw.nports == n
    np.testing.assert_allclose(ntw.f, freqs, rtol=0.0, atol=FREQ_ATOL_HZ)
    np.testing.assert_allclose(ntw.z0, Z0_OHM, rtol=0.0, atol=0.0)
    assert err <= ANCHOR_TOL


@pytest.mark.parametrize("case", list(CASES), ids=list(CASES))
def test_skrf_writer_read_by_load_touchstone(tmp_path, case):
    freqs, s_stack = CASES[case]
    n = s_stack.shape[1]
    out_dir = _rank_dir(tmp_path)
    frequency = skrf.Frequency.from_f(freqs, unit="Hz")
    ntw = skrf.Network(frequency=frequency, s=s_stack, z0=Z0_OHM, name="skrf_out")
    ntw.write_touchstone(filename="skrf_out", dir=str(out_dir), form="ri")
    written = sorted(out_dir.glob(f"skrf_out.s{n}p"))
    assert written, f"skrf wrote nothing matching skrf_out.s{n}p in {out_dir}"
    f_read, s_read, port_ids = load_touchstone(written[0])
    err = _max_abs(s_read, s_stack)
    print(f"\n[PORT-7 skrf->load_touchstone] {case}: N={n} max|dS|={err:.3e}")
    assert len(port_ids) == n
    np.testing.assert_allclose(f_read, freqs, rtol=0.0, atol=FREQ_ATOL_HZ)
    assert err <= ANCHOR_TOL


# --- negative control -------------------------------------------------------


def test_transpose_negative_control(tmp_path):
    for case in ("a_records_4x4", "b_record_3x3_C"):
        _, s_stack = CASES[case]
        asym = _max_abs(s_stack, np.swapaxes(s_stack, 1, 2))
        print(
            f"\n[PORT-7 control, predicted/printed] {case}: max|S - S^T|={asym:.3e} "
            "(predicted at round-off: the reciprocal records cannot see layout)"
        )
    for n in (2, 3, 4):
        case = f"c_nonrecip_{n}port"
        freqs, s_stack = CASES[case]
        s_t = np.swapaxes(s_stack, 1, 2)
        asym = _max_abs(s_stack, s_t)
        path = _write_ours(freqs, s_stack, _rank_dir(tmp_path) / case)
        ntw = skrf.Network(str(path))
        miss = _max_abs(ntw.s, s_t)
        print(
            f"[PORT-7 control, asserted] {case}: max|S - S^T|={asym:.3e} "
            f"(floor {TRANSPOSE_PREMISE_FLOOR}); read-vs-S^T miss={miss:.3e} "
            f"(floor {TRANSPOSE_MISS_FLOOR}, ceiling {TRANSPOSE_MISS_CEILING})"
        )
        assert asym >= TRANSPOSE_PREMISE_FLOOR
        assert TRANSPOSE_MISS_FLOOR <= miss <= TRANSPOSE_MISS_CEILING
