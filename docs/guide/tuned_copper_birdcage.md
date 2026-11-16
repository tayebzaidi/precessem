# The tuned copper birdcage

## What you will compute

The 10 g mass-averaged SAR at one point in the gel phantom of the human-scale
16-leg copper birdcage, tuned at 64 MHz and driven in quadrature from two
feeds. Every one of the 32 leg gaps carries a lumped capacitor with a 0.2 Ω
series resistance, modelled as a sheet whose impedance is
`Z = r + 1/(jωC)`; the resistance is where part of the drive power is lost
outside the phantom and the copper.

The copper is not meshed: the coil is cut out of the air and its surface is a
Leontovich (surface-impedance) wall. The two feeds sit 90° apart and are
driven with weights `(1, -1j)` — the co-rotating (`ccw`) sense. A two-feed
drive keeps exactly one mirror of the coil, the plane bisecting the two
feeds, so the readout point sits on a feed azimuth rather than on that plane.

## Before you start

- The complex DolfinX build: `source /usr/local/bin/dolfinx-complex-mode`.
- 8 MPI ranks, ≈ 3–4 min.
- `degree=1` is set explicitly. The `Study` default is degree 2, which this
  coil cannot afford at this width.

## The study

```python tuned-copper-birdcage
import numpy as np

from fem_em_solver import (
    FHumanBirdcageGeometry,
    GelledSalinePhantomMaterial,
    LeontovichWall,
    Study,
    TwoPortQuadratureDrive,
)

geometry = FHumanBirdcageGeometry(as_hole=True)  # copper coil cut out of the air
phantom = GelledSalinePhantomMaterial(0.5, 78.0, 64e6)  # sigma S/m, eps_r, Hz
c_in = 118.2437 * 1e-12  # the in-model Larmor capacitance on copper, F
drive = TwoPortQuadratureDrive(c_in, esr_ohm=0.2)  # 0.2 ohm ESR on every gap

study = Study(geometry, phantom, drive, frequency_hz=64e6,
              conductor=LeontovichWall(), degree=1)
az = np.radians(11.25)  # the P17 feed azimuth
c_star = (0.045 * np.cos(az), 0.045 * np.sin(az), 0.0)  # r = 45 mm, mid-plane
sar = study.sar(mass_g=10.0, centres=[c_star], drives=("ccw",))
print(f"SAR10g(c*) = {sar['ccw'][0]['averaged_sar_w_per_kg']:.4e} W/kg")
```

## Run it

Save the block as a Python script and run it inside the container, in complex
mode, on 8 ranks (`<script>` is the file you saved):

```bash
source /usr/local/bin/dolfinx-complex-mode
PYTHONPATH=/workspace/src mpiexec -n 8 python3 <script>
```

## What you should see

```text tuned-copper-birdcage-output
SAR10g(c*) = 4.8529e-03 W/kg
```

## What this does not show

No tuning-accuracy, `f₀`, Q or match claim; no compliance or absolute SAR
claim (the figure is per the package's unit drive). c\* is the argmax of a
fixed lattice of ball centres, not a global SAR maximum. Degree 1 only.

## How this page is tested

The block above is extracted from this page verbatim and executed at 8 ranks
by
[`test_wf17_tuned_copper_tutorial.py`](../../tests/package/test_wf17_tuned_copper_tutorial.py);
its SAR must reproduce the record of
[`test_wf16_study_tuned_f_human_b1_sar.py`](../../tests/package/test_wf16_study_tuned_f_human_b1_sar.py).
