"""Example (`EX-157`): a plane wave in a broadened two-pole (Cole-Cole) dielectric.

`MAT-3` step 1b's fixture: `TH-6`'s lossy plane wave (a 0.1 m box, the exact wave
imposed on every face, interior decay/phase slopes fitted) solved in a two-pole,
broadened dielectric with ionic conductivity,

    eps_hat(f) = eps_inf + sum_n d_eps_n / (1 + (j 2 pi f tau_n)^(1 - alpha_n))
                 - j sigma / (2 pi f eps0),

a gamma-like pole near 100 MHz (alpha = 0.30) and a large beta-like pole near
1 MHz (alpha = 0.35).  ``th:20`` (`EX-145`) shows one Debye pole (alpha = 0);
this is the shape Gabriel tissues take.  The decay alpha and phase beta must
match k0 sqrt(eps_hat(f)) at each Larmor frequency, and a Debye-reduced twin
(every alpha_n = 0) or eps_hat frozen at the wrong frequency must miss.

**It asserts, and it does not re-implement** (the `ANS-1` rule).  The
frequencies, band, medium and helpers come by name from
``tests/validation/test_mat3_cole_cole_plane_wave.py`` (``_solve_plane_wave``
from ``tests/validation/test_lossy_plane_wave.py``).  Asserted as the gate
asserts them: (Ev) the evaluator against the independent ``cmath`` sum
(the gate's literal 1e-12, ``test_mat3_cole_cole_plane_wave.py:109``);
(i) 127.74 MHz, 12^3 -> 24^3, rel L2 falling, rate > 0.9, rel L2 / alpha / beta
inside ``BAND``; (ii) 63.87 MHz at 24^3, alpha and beta inside ``BAND``;
(a') the Debye-reduced medium and (b') eps_hat(63.87 MHz) solved at
127.74 MHz each miss alpha by more than ``BAND``.

Run it::

    ./run_examples.sh -e th:23 -n 2 -t 420

Output: ``examples/time_harmonic/paraview_output/
time_harmonic_23_cole_cole_plane_wave_combined.xdmf`` (CG1 ``E_real``,
``E_imag``, ``E_magnitude`` of the 24^3 solve at 127.74 MHz; qualitative).
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
from fem_em_solver.post.setup_figure import write_setup_figure  # noqa: E402

from tests.validation.test_lossy_plane_wave import _solve_plane_wave  # noqa: E402
from tests.validation.test_mat3_cole_cole_plane_wave import (  # noqa: E402
    BAND,
    COLE_COLE,
    DEBYE_REDUCED,
    EPS_INF,
    F_HIGH,
    F_LOW,
    POLES,
    SIGMA_IONIC,
    _eps_cmath,
    _solve,
    _truth,
    _wrong_f_truth,
)

OUTPUT_DIR = Path(__file__).resolve().parent / "paraview_output"
FIGURE_DIR = Path(__file__).resolve().parent / "figures"
BASENAME = "time_harmonic_23_cole_cole_plane_wave"

#: The gate's literal evaluator bound, tests/validation/test_mat3_cole_cole_plane_wave.py:109
#: (it has no module name to import; cited, not chosen here).
EV_BOUND = 1e-12

#: The gate's record log, cited for the reader; nothing is compared for equality.
RECORD_NOTE = "W1b = 20261026T094151Z_MAT-3-step1b-W1b.log"


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
    _say(comm, "EX-157 -- a plane wave in a two-pole Cole-Cole dielectric: decay = k0 sqrt(eps_hat(f))")
    _say(comm, "=" * 78)
    _say(comm, f"[material] eps_inf={EPS_INF:g}, sigma_ionic={SIGMA_IONIC:g} S/m "
               "(imported from the MAT-3 step 1b gate module)")
    for i, (d, tau, a) in enumerate(POLES, 1):
        _say(comm, f"[material] pole {i}: d_eps={d:g}, tau={tau:.6e} s "
                   f"(f_c={1.0 / (2.0 * np.pi * tau) / 1e6:.3g} MHz), alpha={a:g}")
    _say(comm, f"[record] gate log: {RECORD_NOTE}")

    # ---- (Ev) + closed form at both Larmor frequencies -----------------------
    ev_rel = {}
    truth = {}
    for f in (F_HIGH, F_LOW):
        ev = COLE_COLE.complex_permittivity(f)
        ref = _eps_cmath(f)
        ev_rel[f] = abs(ev - ref) / abs(ref)
        truth[f] = _truth(f, COLE_COLE)
        _truth(f, DEBYE_REDUCED)
    alpha_hi, beta_hi, _ = truth[F_HIGH]
    alpha_lo, beta_lo, _ = truth[F_LOW]
    a_dr_cf = _truth(F_HIGH, DEBYE_REDUCED)[0]

    _say(comm, "\n[closed form] eps_hat(f) of the Cole-Cole medium and of its Debye-reduced twin "
               "(every alpha_n = 0):")
    for f in (F_HIGH, F_LOW):
        eps_cc = COLE_COLE.complex_permittivity(f)
        eps_db = DEBYE_REDUCED.complex_permittivity(f)
        a, b, _ = truth[f]
        _say(comm, f"  {f / 1e6:7.2f} MHz  Cole-Cole {eps_cc.real:9.3f} {eps_cc.imag:+10.3f}j   "
                   f"Debye-reduced {eps_db.real:9.3f} {eps_db.imag:+10.3f}j   "
                   f"alpha={a:.4f} beta={b:.4f}   |ev-cmath|/|cmath|={ev_rel[f]:.3e}")
    a_wf_cf, _ = _wrong_f_truth()
    _say(comm, f"[closed form] ceilings at 127.74 MHz: Debye-reduced alpha miss "
               f"{abs(a_dr_cf - alpha_hi) / alpha_hi:.4%}, "
               f"wrong-f alpha miss {abs(a_wf_cf - alpha_hi) / alpha_hi:.4%}")

    # ---- the gate's solves ---------------------------------------------------
    t = time.perf_counter()
    rel_c, a_c, b_c, n_c = _solve(12, F_HIGH, COLE_COLE)
    _say(comm, f"\n[solve] 12^3 @127.74 ({n_c} cells): rel L2 {rel_c:.5e}, "
               f"alpha {a_c:.5f}, beta {b_c:.5f} ({time.perf_counter() - t:.1f} s)")

    t = time.perf_counter()
    hm_hi = COLE_COLE.as_homogeneous(F_HIGH)
    rel_f, a_f, b_f, n_f, msh, fields = _solve_plane_wave(
        24, hm_hi.sigma, return_fields=True, frequency_hz=F_HIGH, epsilon_r=hm_hi.epsilon_r
    )
    fine_s = time.perf_counter() - t

    write_setup_figure(
        msh,
        None,
        FIGURE_DIR / f"{BASENAME}_setup.png",
        region_names={0: "box (homogeneous two-pole Cole-Cole dielectric, no cell-tag regions)"},
        slice_normal=(0.0, 1.0, 0.0),
        title="th:23 -- plane wave in a two-pole Cole-Cole dielectric: 0.1 m cube, 24^3, 127.74 MHz",
        comm=comm,
    )

    rate = float(np.log(rel_c / rel_f) / np.log(2.0))
    ea_f = abs(a_f - alpha_hi) / alpha_hi
    eb_f = abs(b_f - beta_hi) / beta_hi
    _say(comm, f"[solve] 24^3 @127.74 ({n_f} cells): rel L2 {rel_f:.5e}, alpha {a_f:.5f}, "
               f"beta {b_f:.5f} ({fine_s:.1f} s)")
    _say(comm, f"[rate]  L2 rate {rate:.4f}")

    t = time.perf_counter()
    rel_l, a_l, b_l, n_l = _solve(24, F_LOW, COLE_COLE)
    ea_l = abs(a_l - alpha_lo) / alpha_lo
    eb_l = abs(b_l - beta_lo) / beta_lo
    _say(comm, f"[solve] 24^3 @63.87 ({n_l} cells): rel L2 {rel_l:.5e}, alpha {a_l:.5f}, "
               f"beta {b_l:.5f} ({time.perf_counter() - t:.1f} s)")

    t = time.perf_counter()
    _, a_dr, b_dr, _ = _solve(16, F_HIGH, DEBYE_REDUCED)
    miss_a_dr = abs(a_dr - alpha_hi) / alpha_hi
    miss_b_dr = abs(b_dr - beta_hi) / beta_hi
    _say(comm, f"[control a'] Debye-reduced, 16^3: alpha {a_dr:.5f}, alpha miss {miss_a_dr:.4%}, "
               f"beta miss {miss_b_dr:.4%} ({time.perf_counter() - t:.1f} s)")

    t = time.perf_counter()
    _, a_wf, b_wf, _ = _solve(16, F_HIGH, COLE_COLE, material_f=F_LOW)
    miss_a_wf = abs(a_wf - alpha_hi) / alpha_hi
    miss_b_wf = abs(b_wf - beta_hi) / beta_hi
    _say(comm, f"[control b'] eps_hat(63.87 MHz) solved at 127.74 MHz, 16^3: alpha {a_wf:.5f}, "
               f"alpha miss {miss_a_wf:.4%}, beta miss {miss_b_wf:.4%} ({time.perf_counter() - t:.1f} s)")

    # ---- gates ---------------------------------------------------------------
    _say(comm, f"\n(Ev)   evaluator vs cmath: {max(ev_rel.values()):.3e} <= {EV_BOUND:g} "
               "(gate's literal bound, :109): ASSERTED")
    _say(comm, f"(i)    127.74 MHz 24^3: rel L2 {rel_f:.4e} < 12^3's {rel_c:.4e}, rate {rate:.4f} > 0.9, "
               f"alpha {ea_f:.4%}, beta {eb_f:.4%}; each < {BAND:.0%}: ASSERTED")
    _say(comm, f"(ii)   63.87 MHz 24^3: alpha {ea_l:.4%}, beta {eb_l:.4%}; each < {BAND:.0%}: ASSERTED")
    _say(comm, f"(a')   Debye-reduced alpha miss {miss_a_dr:.4%} > {BAND:.0%} "
               f"({miss_a_dr / BAND:.1f}x the band): ASSERTED")
    _say(comm, f"(b')   wrong-frequency alpha miss {miss_a_wf:.4%} > {BAND:.0%} "
               f"({miss_a_wf / BAND:.1f}x the band): ASSERTED")

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

    for f, r in ev_rel.items():
        assert r <= EV_BOUND, f"evaluator vs cmath at {f:.4e} Hz: {r:.3e}"
    assert rel_f < rel_c, f"rel L2 did not fall: {rel_c:.4e} -> {rel_f:.4e}"
    assert rate > 0.9, f"L2 rate {rate:.3f}"
    assert rel_f < BAND, f"24^3 rel L2 {rel_f:.4e}"
    assert ea_f < BAND, f"alpha @127.74 off by {ea_f:.3%}"
    assert eb_f < BAND, f"beta @127.74 off by {eb_f:.3%}"
    assert ea_l < BAND, f"alpha @63.87 off by {ea_l:.3%}"
    assert eb_l < BAND, f"beta @63.87 off by {eb_l:.3%}"
    assert miss_a_dr > BAND, f"Debye-reduced control inside the band: {miss_a_dr:.3%}"
    assert miss_a_wf > BAND, f"wrong-frequency control inside the band: {miss_a_wf:.3%}"

    _say(comm, f"\n[timing] total {time.perf_counter() - started:.1f} s at -n {comm.size}")
    _say(comm, "\nAll imported gates hold.")


if __name__ == "__main__":
    main()
