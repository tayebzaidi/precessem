"""Example (`EX-145`): a plane wave in a dispersive (Debye) dielectric.

`MAT-3` step 1's fixture: `TH-6`'s lossy plane wave (a 0.1 m box, the exact wave
imposed on every face, interior decay/phase slopes fitted) solved in a
single-pole Debye dielectric with ionic conductivity,

    eps_hat(f) = eps_inf + d_eps / (1 + j 2 pi f tau) - j sigma / (2 pi f eps0),

evaluated by ``ColeColeMaterial`` and frozen at the solve frequency by
``as_homogeneous(f)``.  The decay alpha and phase beta must match
k0 sqrt(eps_hat(f)) at each Larmor frequency, and eps_hat frozen at the wrong
frequency, or dispersion dropped altogether, must miss by more than the band.

**It asserts, and it does not re-implement** (the `ANS-1` rule).  The
frequencies, the band, the material and the helpers come by name from
``tests/validation/test_mat3_debye_plane_wave.py`` (``_solve_plane_wave`` from
``tests/validation/test_lossy_plane_wave.py``, ``ColeColeMaterial`` from the
package).  Asserted as the gate asserts them: (i) 127.74 MHz at 24^3, rel L2,
alpha and beta inside ``BAND``; (ii) 63.87 MHz at 24^3, alpha and beta inside
``BAND``; (c-a) eps_hat(63.87 MHz) solved at 127.74 MHz misses alpha by more
than ``BAND``; (c-b) dispersion dropped misses alpha by more than ``BAND``.
The L2 rate and every record are printed beside their readings and never
asserted for equality.

Run it::

    ./run_examples.sh -e th:20 -n 2 -t 600

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_20_debye_plane_wave_combined.xdmf`` (CG1 ``E_real``, ``E_imag``,
``E_magnitude`` of the 24^3 solve at 127.74 MHz; a qualitative picture).
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

from fem_em_solver.io.paraview_utils import (  # noqa: E402
    adopt_host_ownership,
    write_xdmf_with_tags,
)
from fem_em_solver.materials import ColeColeMaterial  # noqa: E402
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_lossy_plane_wave import _solve_plane_wave  # noqa: E402
from tests.validation.test_mat3_debye_plane_wave import (  # noqa: E402
    BAND,
    DEBYE,
    DELTA_EPS,
    EPS_INF,
    F_HIGH,
    F_LOW,
    SIGMA_IONIC,
    TAU,
    _solve,
    _truth,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_20_debye_plane_wave"

#: Gate records -- printed beside this run, a difference reported and never asserted.
#: W1 = 20260929T183252Z_MAT-3-step1.log.gz
RECORD = {
    "hi": (61.3993, 84.9196, 0.60348, 12.4703, 24.4047),   # W1:157 eps', eps'', sigma_eff, alpha, beta
    "lo": (71.3078, 154.3256, 0.54836, 9.4035, 14.7038),   # W1:158
    "ratio_closed": 1.3261,                                 # W1:159
    "ratio_measured": 1.3259,                               # W1:163
    "cells_12": 10368,                                      # W1:160
    "cells_24": 82944,                                      # W1:161
    "rel_24": 3.29632e-02,                                  # W1:161
    "rate": 0.9998,                                         # W1:161
    "ea_hi": 0.0262, "eb_hi": 0.0554,                       # percent, W1:161
    "ea_lo": 0.0419, "eb_lo": 0.0475,                       # percent, W1:162
    "miss_a_wf": 11.9990,                                   # percent, W1:164
    "miss_a_nd": 21.9949,                                   # percent, W1:165
}


def _say(comm, msg):
    if comm.rank == 0:
        print(msg, flush=True)


def _paraview_fields(msh, fields):
    """CG1 ``E_real`` / ``E_imag`` / ``E_magnitude`` from the solved phasor."""
    v_cg = fem.functionspace(msh, ("Lagrange", 1, (3,)))
    e_cg = fem.Function(v_cg, name="E_phasor")
    e_cg.interpolate(fields.e_complex)
    e_cg.x.scatter_forward()
    e_re = fem.Function(v_cg, name="E_real")
    e_re.x.array[:] = np.real(e_cg.x.array)
    e_im = fem.Function(v_cg, name="E_imag")
    e_im.x.array[:] = np.imag(e_cg.x.array)
    s_cg = fem.functionspace(msh, ("Lagrange", 1))
    e_mag = fem.Function(s_cg, name="E_magnitude")
    comp = np.abs(e_cg.x.array.reshape(-1, 3))
    e_mag.x.array[:] = np.sqrt(np.sum(comp * comp, axis=1))
    for f in (e_re, e_im, e_mag):
        f.x.scatter_forward()
    return e_re, e_im, e_mag


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
    _say(comm, "EX-145 -- a plane wave in a Debye dielectric: decay = k0 sqrt(eps_hat(f))")
    _say(comm, "=" * 78)
    _say(comm, f"[material] Debye eps_inf={EPS_INF:g}, d_eps={DELTA_EPS:g}, tau={TAU:.6e} s, "
               f"sigma_ionic={SIGMA_IONIC:g} S/m (imported from the MAT-3 step 1 gate module)")

    # ---- closed form at both Larmor frequencies ------------------------------
    alpha_hi, beta_hi, eps_hi = _truth(F_HIGH)
    alpha_lo, beta_lo, eps_lo = _truth(F_LOW)
    sig_hi, _ = DEBYE.sigma_epsilon_at(F_HIGH)
    sig_lo, _ = DEBYE.sigma_epsilon_at(F_LOW)
    for tag, f, eps, sig, al, be in (
        ("hi", F_HIGH, eps_hi, sig_hi, alpha_hi, beta_hi),
        ("lo", F_LOW, eps_lo, sig_lo, alpha_lo, beta_lo),
    ):
        r = RECORD[tag]
        _say(comm, f"[closed form] {f / 1e6:.2f} MHz: eps'={eps.real:.4f} eps''={-eps.imag:.4f} "
                   f"sigma_eff={sig:.5f} alpha={al:.4f} beta={be:.4f}")
        _say(comm, f"[record]      {f / 1e6:.2f} MHz: {r[0]:.4f} / {r[1]:.4f} / {r[2]:.5f} / "
                   f"{r[3]:.4f} / {r[4]:.4f}")
    _say(comm, f"[closed form] alpha ratio 128/64 = {alpha_hi / alpha_lo:.4f} "
               f"(record {RECORD['ratio_closed']})")

    # ---- the gate's solves ---------------------------------------------------
    t = time.perf_counter()
    rel_c, a_c, b_c, n_c = _solve(12, F_HIGH)
    _say(comm, f"[solve] 12^3 @127.74 ({n_c} cells, record {RECORD['cells_12']}): rel L2 {rel_c:.5e}, "
               f"alpha {a_c:.5f}, beta {b_c:.5f} ({time.perf_counter() - t:.1f} s)")

    t = time.perf_counter()
    hm_hi = DEBYE.as_homogeneous(F_HIGH)
    rel_f, a_f, b_f, n_f, msh, fields = _solve_plane_wave(
        24, hm_hi.sigma, return_fields=True, frequency_hz=F_HIGH, epsilon_r=hm_hi.epsilon_r
    )
    fine_s = time.perf_counter() - t

    write_setup_figure(
        msh,
        None,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={0: "box (homogeneous Debye dielectric, no cell-tag regions)"},
        slice_normal=(0.0, 1.0, 0.0),
        title="th:20 -- plane wave in a Debye dielectric: 0.1 m cube, 24^3, 127.74 MHz",
        comm=comm,
    )

    rate = float(np.log(rel_c / rel_f) / np.log(2.0))
    ea_f = abs(a_f - alpha_hi) / alpha_hi
    eb_f = abs(b_f - beta_hi) / beta_hi
    _say(comm, f"[solve] 24^3 @127.74 ({n_f} cells, record {RECORD['cells_24']}): rel L2 {rel_f:.5e} "
               f"(record {RECORD['rel_24']:.5e}), alpha {a_f:.5f}, beta {b_f:.5f} ({fine_s:.1f} s)")
    _say(comm, f"[rate]  L2 rate {rate:.4f} (record {RECORD['rate']}; printed, not asserted)")

    t = time.perf_counter()
    rel_l, a_l, b_l, n_l = _solve(24, F_LOW)
    ea_l = abs(a_l - alpha_lo) / alpha_lo
    eb_l = abs(b_l - beta_lo) / beta_lo
    _say(comm, f"[solve] 24^3 @63.87 ({n_l} cells): rel L2 {rel_l:.5e}, alpha {a_l:.5f}, "
               f"beta {b_l:.5f} ({time.perf_counter() - t:.1f} s)")
    _say(comm, f"[ratio] measured alpha 128/64 = {a_f / a_l:.4f} (record {RECORD['ratio_measured']})")

    t = time.perf_counter()
    _, a_wf, b_wf, _ = _solve(16, F_HIGH, material_f=F_LOW)
    miss_a_wf = abs(a_wf - alpha_hi) / alpha_hi
    miss_b_wf = abs(b_wf - beta_hi) / beta_hi
    _say(comm, f"[control a] eps_hat(63.87 MHz) solved at 127.74 MHz, 16^3: alpha {a_wf:.5f}, "
               f"miss {miss_a_wf:.4%} (record {RECORD['miss_a_wf']:.4f}%); beta miss {miss_b_wf:.4%} "
               f"({time.perf_counter() - t:.1f} s)")

    t = time.perf_counter()
    static = ColeColeMaterial.debye(EPS_INF + DELTA_EPS, sigma_ionic=SIGMA_IONIC)
    _, a_nd, b_nd, _ = _solve(16, F_HIGH, hm=static.as_homogeneous(F_HIGH))
    miss_a_nd = abs(a_nd - alpha_hi) / alpha_hi
    miss_b_nd = abs(b_nd - beta_hi) / beta_hi
    _say(comm, f"[control b] dispersion dropped (eps_r = eps_s), 16^3: alpha {a_nd:.5f}, "
               f"miss {miss_a_nd:.4%} (record {RECORD['miss_a_nd']:.4f}%); beta miss {miss_b_nd:.4%} "
               f"({time.perf_counter() - t:.1f} s)")

    # ---- gates ---------------------------------------------------------------
    _say(comm, f"\n(i)    127.74 MHz 24^3: rel L2 {rel_f:.4e}, alpha {ea_f:.4%} "
               f"(record {RECORD['ea_hi']}%), beta {eb_f:.4%} (record {RECORD['eb_hi']}%); "
               f"each < {BAND:.0%}: ASSERTED")
    _say(comm, f"(ii)   63.87 MHz 24^3: alpha {ea_l:.4%} (record {RECORD['ea_lo']}%), "
               f"beta {eb_l:.4%} (record {RECORD['eb_lo']}%); each < {BAND:.0%}: ASSERTED")
    _say(comm, f"(c-a)  wrong-frequency alpha miss {miss_a_wf:.4%} > {BAND:.0%} "
               f"({miss_a_wf / BAND:.1f}x the band): ASSERTED")
    _say(comm, f"(c-b)  dispersion-dropped alpha miss {miss_a_nd:.4%} > {BAND:.0%} "
               f"({miss_a_nd / BAND:.1f}x the band): ASSERTED")

    # ---- ParaView ------------------------------------------------------------
    e_re, e_im, e_mag = _paraview_fields(msh, fields)
    OUTPUT_DIR.mkdir(exist_ok=True)
    path, _ = write_xdmf_with_tags(
        OUTPUT_DIR / f"{BASENAME}_combined",
        msh,
        None,
        {"E_real": e_re, "E_imag": e_im, "E_magnitude": e_mag},
        comm=comm,
    )
    adopt_host_ownership(OUTPUT_DIR, comm=comm)
    _say(comm, f"[output] {path}")

    assert rel_f < BAND, f"24^3 rel L2 {rel_f:.4e}"
    assert ea_f < BAND, f"alpha @127.74 off by {ea_f:.3%}"
    assert eb_f < BAND, f"beta @127.74 off by {eb_f:.3%}"
    assert ea_l < BAND, f"alpha @63.87 off by {ea_l:.3%}"
    assert eb_l < BAND, f"beta @63.87 off by {eb_l:.3%}"
    assert miss_a_wf > BAND, f"wrong-frequency control inside the band: {miss_a_wf:.3%}"
    assert miss_a_nd > BAND, f"dispersion-dropped control inside the band: {miss_a_nd:.3%}"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
