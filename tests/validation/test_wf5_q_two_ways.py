"""`WF-5` step 1 (§9 item 134, chain T14): Q two ways at f₀ on `ANS-7`'s Route A.

The reactance theorem (Foster's) says a one-port's susceptance slope at a
parallel resonance is twice the stored energy over |V|², so the loaded Q of
the tuned F-small birdcage's input admittance,

    Q_slope = ω₀ (dB/dω) / (2G)          (``parallel_rows`` on Route A),

must equal the field's own

    Q_field = ω₀ (W_m + W_e + W_sheets + W_C1 + W_surf,F) / P_loss

from one Route A solve at f₀ (per column, Cu and PEC):

* W_m = ∫|curl E|²/(4ω²μ₀), W_e = (ε₀/4)∫ε_r′|E|²;
* W_sheets (P2–P4, the C-terminated gaps) = −Re a_sheet(E,E)/(4ω²μ₀) through
  the package's ``lumped_port_bilinear_term`` on the **κ-corrected specs as
  solved** — the form the matrix carried, never the terminal voltage;
* W_C1 = ¼C₁|V₁|² (C₁ is the circuit-layer capacitor across the port);
* W_surf,F = ½ · ¼(Im Z_s/ω)∫|n×E|²/|Z_s|² (Cu; Foster's ½ for X_s ∝ √ω);
* P_loss = ½∫σ|E|² (phantom) + ``_surface_loss_w``; P1's 50 Ω sheet
  dissipation is outside the network and is excluded.

Band (``test_band_candidate``): the secant's own error against a cubic
through fresh Route A solves at 63.8 … 64.1 MHz, ``B_Q = ceil_1sf(max(3·e_sec,
1e-6))`` over both columns, frozen from W1 into ``B_Q_FROZEN``.

Imports from the `ans:7` script by importlib — never ``_assemble`` /
``_write_outputs``.  Tests only.  Degree 1, complex mode.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
import ufl
from dolfinx import fem
from mpi4py import MPI

from fem_em_solver.utils.constants import EPSILON_0, MU_0
from tests.complex_mode import complex_only

_REPO_ROOT = Path(__file__).resolve().parents[2]
ANS7_DIR = _REPO_ROOT / "examples" / "ansys_benchmarks" / "ans7_tuned_copper_birdcage_driven_sweep"
ANS7_SCRIPT = ANS7_DIR / "07_tuned_copper_birdcage_driven_sweep.py"
COLUMNS = ("cu", "pec")

# The cubic's support (§9 item 134 "Band").
CUBIC_MHZ = (63.8, 63.9, 64.0, 64.1)
SECANT_MHZ = (63.9, 64.0)
# Anchor: fresh f0 / Q on the secant reproduce the tracked records (item 79's
# re-run reproduced every metrics.json leaf to <= 4.2e-9).
RECORD_RTOL = 1.0e-8
BAND_FLOOR = 1.0e-6
# *Predicted, printed only* (§9 item 134).
PREDICTED_E_SEC = {"cu": 5.0e-4, "pec": 7.2e-4}
PREDICTED_SHEET_SHARE = (0.20, 0.30)
PREDICTED_W_SURF_SHARE = 4.4e-4
# Frozen from W1 (``20261001T141329Z_WF-5-step1-W1.log:960–962``): e_sec Cu
# 5.8473e-4, PEC 8.2352e-4 ⇒ ceil_1sf(3 × 8.2352e-4) = 3e-3.  Never loosened.
B_Q_FROZEN = 3.0e-3


def _load(path: Path, name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def ceil_1sf(x: float) -> float:
    """Round ``x > 0`` up to one significant figure."""
    e = math.floor(math.log10(x))
    m = x / 10.0**e
    m_up = math.ceil(m - 1e-12)
    return float(m_up * 10.0**e)


def _reduce(comm, form) -> complex:
    return complex(comm.allreduce(fem.assemble_scalar(fem.form(form)), op=MPI.SUM))


@pytest.fixture(scope="module")
def ans7():
    return _load(ANS7_SCRIPT, "ans7_driven_sweep_wf5")


@pytest.fixture(scope="module")
def fixture_build(ans7):
    """`_stage_a:690–700`'s build, κ from the tracked metrics.json."""
    from tests.validation import test_th14_birdcage_copper as th14

    comm = MPI.COMM_WORLD
    tracked = json.loads((ANS7_DIR / "metrics.json").read_text())
    kappa = float(tracked["kappa"])
    msh, cell_tags, facet_tags, _diag, _t_mesh = th14._sheets_build(True, as_hole=True)
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    msh.topology.create_entity_permutations()
    outer_tags, _census = th14._outer_box_tags(msh, facet_tags, comm)
    tags_k, port_defs_k, specs_k = th14._build_ports(
        msh, cell_tags, comm, width_correction_kappa=kappa
    )
    geo_k = (msh, cell_tags, facet_tags, outer_tags, tags_k, port_defs_k)
    if comm.rank == 0:
        print(
            f"[wf5] hole mesh {int(msh.topology.index_map(tdim).size_global)} cells; "
            f"kappa (tracked) = {kappa:.9e}",
            flush=True,
        )
    return {"th14": th14, "geo": geo_k, "specs": specs_k, "tracked": tracked, "comm": comm}


def _y_in(ans7, s11, f_hz):
    return 1.0 / ans7.add_c1(ans7._z_from_s11(s11), 2.0 * np.pi * f_hz)


@pytest.fixture(scope="module")
def sweep4(ans7, fixture_build):
    """Fresh Route A at the cubic's four points, both columns."""
    fb = fixture_build
    out = {}
    for col in COLUMNS:
        y = []
        for f_mhz in CUBIC_MHZ:
            s11, _ = ans7._route_a_solve(fb["th14"], fb["geo"], fb["specs"], f_mhz * 1e6, col)
            y.append(_y_in(ans7, s11, f_mhz * 1e6))
        out[col] = np.array(y)
    return out


def _secant_and_cubic(ans7, y):
    f = np.array(CUBIC_MHZ) * 1e6
    i0, i1 = CUBIC_MHZ.index(SECANT_MHZ[0]), CUBIC_MHZ.index(SECANT_MHZ[1])
    z_sec = 1.0 / y[[i0, i1]]
    rows = ans7.parallel_rows(f[[i0, i1]], z_sec)
    x = (f - 64.0e6) / 1.0e6
    pb = np.polyfit(x, y.imag, 3)
    pg = np.polyfit(x, y.real, 3)
    roots = [r.real for r in np.roots(pb) if abs(r.imag) < 1e-12 and -0.1 <= r.real <= 0.0]
    assert len(roots) == 1, f"cubic B: expected one root in [63.9, 64.0] MHz, got {roots}"
    x0 = roots[0]
    f0_c = 64.0e6 + x0 * 1.0e6
    dbdf = np.polyval(np.polyder(pb), x0) / 1.0e6
    q_cubic = f0_c * dbdf / (2.0 * np.polyval(pg, x0))
    return rows, float(f0_c), float(q_cubic)


@complex_only
def test_band_candidate(ans7, fixture_build, sweep4):
    """Anchors: fresh secant f0 / Q reproduce the records; B_Q from e_sec."""
    comm = fixture_build["comm"]
    tracked = fixture_build["tracked"]["rows_route_a"]
    rank0 = comm.rank == 0
    e_sec = {}
    fails = []
    for col in COLUMNS:
        rows, f0_c, q_cubic = _secant_and_cubic(ans7, sweep4[col])
        df0 = abs(rows["f0_hz"] - tracked[col]["f0_hz"]) / tracked[col]["f0_hz"]
        dq = abs(rows["q_loaded"] - tracked[col]["q_loaded"]) / tracked[col]["q_loaded"]
        e_sec[col] = abs(rows["q_loaded"] - q_cubic) / q_cubic
        if rank0:
            print(
                f"[wf5] {col}: secant f0 {rows['f0_hz']:.6f} Hz (record "
                f"{tracked[col]['f0_hz']:.6f}, rel {df0:.3e}); Q_secant {rows['q_loaded']:.10f} "
                f"(record {tracked[col]['q_loaded']:.10f}, rel {dq:.3e}; ASSERTED <= "
                f"{RECORD_RTOL:.0e}); cubic f0 {f0_c:.6f} Hz, Q_cubic {q_cubic:.10f}; "
                f"e_sec {e_sec[col]:.4e} (PREDICTED {PREDICTED_E_SEC[col]:.1e})",
                flush=True,
            )
        if not (df0 <= RECORD_RTOL and dq <= RECORD_RTOL):
            fails.append((col, df0, dq))
    b_q = ceil_1sf(max(3.0 * max(e_sec.values()), BAND_FLOOR))
    if rank0:
        print(
            f"[wf5] B_Q candidate = ceil_1sf(max(3*{max(e_sec.values()):.4e}, "
            f"{BAND_FLOOR:.0e})) = {b_q:.1e}; frozen B_Q_FROZEN = {B_Q_FROZEN}",
            flush=True,
        )
    assert not fails, f"fresh secant rows miss the records: {fails}"
    if B_Q_FROZEN is not None:
        assert b_q <= B_Q_FROZEN, f"B_Q candidate {b_q} exceeds the frozen {B_Q_FROZEN}"


def _solve_with_result(ans7, fb, f_hz, col):
    """``_route_a_solve`` with ``return_fields=True``; the package's result
    (V₁, I₁, terminal voltages) captured by wrapping the solver entry point the
    function imports, so the solve is byte-for-byte the script's."""
    import fem_em_solver.ports.lumped as lumped

    captured = {}
    original = lumped.run_lumped_sheet_port_case

    def wrapper(*args, **kwargs):
        out = original(*args, **kwargs)
        captured["result"] = out[0] if isinstance(out, tuple) else out
        return out

    lumped.run_lumped_sheet_port_case = wrapper
    try:
        s11, fields = ans7._route_a_solve(
            fb["th14"], fb["geo"], fb["specs"], f_hz, col, return_fields=True
        )
    finally:
        lumped.run_lumped_sheet_port_case = original
    return s11, fields, captured["result"]


def _budget(ans7, fb, f_hz, col):
    """The energy / loss budget of one Route A solve at f_hz, all reduced."""
    from fem_em_solver.core.cavity import surface_resistance_ohm
    from fem_em_solver.ports.lumped import lumped_port_bilinear_term
    from tests.validation.test_birdcage_power_identity import (
        _sheet_field_dissipation_w,
        _source_power_w,
    )

    th14 = fb["th14"]
    comm = fb["comm"]
    msh, cell_tags, facet_tags, _outer, tags_f, _pd = fb["geo"]
    omega = 2.0 * np.pi * f_hz
    s11, fields, result = _solve_with_result(ans7, fb, f_hz, col)
    e = fields.e_complex
    z_c = ans7.z_cap(omega)
    term_idx = ans7.TUNING_TERMINATED_INDICES
    specs = [
        dataclasses.replace(s, port_impedance_ohm=z_c) if i in term_idx else s
        for i, s in enumerate(fb["specs"])
    ]
    dx = ufl.Measure("dx", domain=msh)
    w_m = _reduce(comm, ufl.inner(ufl.curl(e), ufl.curl(e)) * dx).real / (4.0 * omega**2 * MU_0)
    eps_r = ufl.real(fields.epsilon_r_field)
    w_e = 0.25 * EPSILON_0 * _reduce(comm, eps_r * ufl.inner(e, e) * dx).real
    w_e_air_only = 0.25 * EPSILON_0 * _reduce(comm, ufl.inner(e, e) * dx).real
    p_phantom = 0.5 * _reduce(comm, ufl.real(fields.sigma_field) * ufl.inner(e, e) * dx).real
    w_sheet = {}
    w_sheet_blind = {}
    for i in term_idx:
        sp = specs[i]
        form = lumped_port_bilinear_term(
            msh, tags_f, sp.sheet(driven=False), e, e, omega_rad_per_s=omega
        )
        w_sheet[sp.port_id] = -_reduce(comm, form).real / (4.0 * omega**2 * MU_0)
        v_t = result.responses[sp.port_id].voltage_v
        w_sheet_blind[sp.port_id] = 0.25 * ans7.C_F * abs(v_t) ** 2
    drive = result.responses[th14.DRIVEN]
    v1, i1 = complex(drive.voltage_v), complex(drive.current_a)
    w_c1 = 0.25 * ans7.C_F * abs(v1) ** 2
    if col == "cu":
        z_s = (1.0 + 1.0j) * surface_resistance_ohm(omega, th14.COPPER_SIGMA)
        p_surf = th14._surface_loss_w(msh, facet_tags, z_s, e, comm)
        ds_c = ufl.Measure(
            "ds",
            domain=msh,
            subdomain_data=facet_tags,
            subdomain_id=th14.BIRDCAGE_CONDUCTOR_SURFACE_TAG,
        )
        n = ufl.FacetNormal(msh)
        tan2 = _reduce(comm, ufl.inner(ufl.cross(n, e), ufl.cross(n, e)) * ds_c).real
        w_surf_naive = 0.25 * (z_s.imag / omega) * tan2 / abs(z_s) ** 2
    else:
        p_surf, w_surf_naive = 0.0, 0.0
    w_surf = 0.5 * w_surf_naive
    driven_sheet = specs[0].sheet(driven=True)
    assert specs[0].port_id == th14.DRIVEN
    sweep = {"mesh": msh, "facet_tags": tags_f}
    p_src = _source_power_w(sweep, driven_sheet, e, omega)
    p_sheet1 = _sheet_field_dissipation_w(sweep, driven_sheet, e, omega)
    return {
        "omega": omega,
        "s11": s11,
        "w_m": w_m,
        "w_e": w_e,
        "w_e_air_only": w_e_air_only,
        "w_sheet": w_sheet,
        "w_sheet_blind": w_sheet_blind,
        "w_c1": w_c1,
        "w_surf": w_surf,
        "w_surf_naive": w_surf_naive,
        "p_phantom": p_phantom,
        "p_surf": p_surf,
        "p_src": p_src,
        "p_sheet1": p_sheet1,
        "v1": v1,
        "i1": i1,
    }


@complex_only
def test_q_two_ways(ans7, fixture_build, sweep4):
    """Q_field (energy / loss integrals) against Q_slope (Im Y_in slope), Cu and PEC."""
    assert B_Q_FROZEN is not None, "B_Q not frozen yet (run W1 first)"
    from tests.validation.test_birdcage_power_identity import DISCRETE_IDENTITY_RTOL

    fb = fixture_build
    comm = fb["comm"]
    rank0 = comm.rank == 0
    fails = []
    for col in COLUMNS:
        rows, _f0c, _qc = _secant_and_cubic(ans7, sweep4[col])
        f0, q_slope = rows["f0_hz"], rows["q_loaded"]
        b = _budget(ans7, fb, f0, col)
        omega = b["omega"]
        p_loss = b["p_phantom"] + b["p_surf"]
        w_sheets = sum(b["w_sheet"].values())
        w_tot = b["w_m"] + b["w_e"] + w_sheets + b["w_c1"] + b["w_surf"]
        q_field = omega * w_tot / p_loss
        # ε printed before the gate is read.
        eps = 0.5 * (b["v1"] * np.conj(b["i1"])).real / p_loss - 1.0
        id_rel = abs(p_loss - (b["p_src"] - b["p_sheet1"])) / p_loss
        miss = abs(q_field - q_slope) / q_slope
        q_noshe = omega * (w_tot - w_sheets) / p_loss
        miss_noshe = abs(q_noshe - q_slope) / q_slope
        q_blind = omega * (w_tot - w_sheets + sum(b["w_sheet_blind"].values())) / p_loss
        q_naive = omega * (w_tot - b["w_surf"] + b["w_surf_naive"]) / p_loss
        if rank0:
            print(
                f"[wf5] {col}: f0 {f0:.6f} Hz; epsilon = Re(V1 I1*)/2/P_loss - 1 = {eps:+.6e} "
                f"(PRINTED, read before the gate)",
                flush=True,
            )
            print(
                f"[wf5] {col}: W_m {b['w_m']:.9e} J, W_e {b['w_e']:.9e} J (eps_r=1 everywhere "
                f"{b['w_e_air_only']:.9e}), W_sheets {w_sheets:.9e} J "
                f"({ {k: round(v, 15) for k, v in b['w_sheet'].items()} }), W_C1 {b['w_c1']:.9e} J, "
                f"W_surf,F {b['w_surf']:.9e} J (naive {b['w_surf_naive']:.9e}); W_tot {w_tot:.9e} J",
                flush=True,
            )
            print(
                f"[wf5] {col}: shares of W_tot: W_m {b['w_m'] / w_tot:.6f}, W_e "
                f"{b['w_e'] / w_tot:.6f}, W_sheets {w_sheets / w_tot:.6f} (control PREDICTED "
                f"{PREDICTED_SHEET_SHARE}), W_C1 {b['w_c1'] / w_tot:.6f}, W_surf "
                f"{b['w_surf'] / w_tot:.4e} (PREDICTED {PREDICTED_W_SURF_SHARE:.1e}); "
                f"W_m/W_e,tot {b['w_m'] / (b['w_e'] + w_sheets + b['w_c1']):.6f}",
                flush=True,
            )
            print(
                f"[wf5] {col}: P_phantom {b['p_phantom']:.9e} W, P_surf {b['p_surf']:.9e} W, "
                f"P_loss {p_loss:.9e} W; p_src {b['p_src']:.9e} W, p_sheet(P1) "
                f"{b['p_sheet1']:.9e} W; identity rel {id_rel:.3e} (ASSERTED <= "
                f"{DISCRETE_IDENTITY_RTOL:.0e}; PREDICTED ~1e-12)",
                flush=True,
            )
            print(
                f"[wf5] {col}: 1/Q split phantom {b['p_phantom'] / (omega * w_tot):.6e}, "
                f"surface {b['p_surf'] / (omega * w_tot):.6e}",
                flush=True,
            )
            print(
                f"[wf5] {col}: Q_field {q_field:.9f}, Q_slope {q_slope:.9f}; "
                f"|Q_field - Q_slope|/Q_slope = {miss:.4e} (ASSERTED <= B_Q {B_Q_FROZEN:.0e})",
                flush=True,
            )
            print(
                f"[wf5] {col}: control W_sheets omitted Q {q_noshe:.6f}, miss {miss_noshe:.4e} "
                f"(ASSERTED > B_Q); kappa-blind 1/4 C|V_term|^2 Q {q_blind:.6f}, miss "
                f"{abs(q_blind - q_slope) / q_slope:.4e} (printed); naive W_surf Q {q_naive:.6f}, "
                f"miss {abs(q_naive - q_slope) / q_slope:.4e} (printed)",
                flush=True,
            )
        if not id_rel <= DISCRETE_IDENTITY_RTOL:
            fails.append((col, "identity", id_rel))
        if not miss <= B_Q_FROZEN:
            fails.append((col, "Q two ways", miss))
        if not miss_noshe > B_Q_FROZEN:
            fails.append((col, "control inside band", miss_noshe))
    assert not fails, f"WF-5 step 1 gates: {fails}"
