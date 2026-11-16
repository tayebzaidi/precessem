"""Diagnostic only (item 824 park): why the precondition |F| <= 1e-12 of scale is red.

Imports the module as run in W1 (unedited) and prints, per case, the worst
root's |F| / scale, the float floor (min |F| over ±8 ulps around the root,
over the same scale), and |F| over the envelope sqrt(Z0^2 + Z1^2)(b)."""

import numpy as np
from scipy import special

import tests.validation.test_imp8_composite_cylinder_closed_form as m8

for name, kr, rc, hole in (("Ti", m8.ROD_K, m8.ROD_RHO_C, False), ("gel", m8.K_GEL, m8.RHO_C, False),
                           ("hole", m8.ROD_K, m8.ROD_RHO_C, True), ("ctrl", m8.ROD_K, m8.RHO_C, False)):
    md = m8._modes(float(kr), float(rc), m8.N_MODES, hole)
    rel = np.abs(md["fval"]) / md["fscale"]
    i = int(np.argmax(rel))
    sig = md["sig"][i]
    nb = sig + np.spacing(sig) * np.arange(-8, 9)
    f, sc = m8._F(nb, float(kr), float(rc), hole)
    A, B, _, mu2 = m8._coeffs(md["sig"], float(kr), float(rc), hole)
    xb = mu2 * m8.B_WALL
    env = np.sqrt((A * special.j0(xb) + B * special.y0(xb)) ** 2 + (A * special.j1(xb) + B * special.y1(xb)) ** 2)
    renv = np.abs(md["fval"]) / env
    print(f"[diag] {name}: worst root n = {i + 1}, sqrt(s) = {sig:.7g}, mu2 b = {xb[i]:.6g}, "
          f"|F|/scale = {rel[i]:.3e}, float floor (min over +-8 ulp)/scale = {np.min(np.abs(f)) / sc[8]:.3e}, "
          f"sign change in +-8 ulp: {bool(np.any(np.sign(f[:-1]) * np.sign(f[1:]) < 0))}, "
          f"max |F|/envelope = {renv.max():.3e} (n = {int(np.argmax(renv)) + 1}), "
          f"|A|, |B| at worst = {abs(A[i]):.3e}, {abs(B[i]):.3e}, rel bracket = {(md['hi'][i] - md['lo'][i]) / md['hi'][i]:.3e}",
          flush=True)
