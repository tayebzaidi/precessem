# Four-port outputs from Python: 10 g SAR, |B₁⁺| and a frequency sweep

## What you will compute

The F-small four-port birdcage (exact port rectangles, κ off, degree 1,
64 MHz) around its saline phantom, driven one port at a time at 1 V, built as
a `Study` in Python rather than from a study file. Three blocks read three
outputs from it: `Study.sar` gives the 10 g mass-averaged SAR in a ball at one
centre per port, placed on a ring at each solved port sheet's azimuth;
`Study.b1_plus` and `Study.b1_plus_at` give `|B₁⁺|` at listed points; and
`Study.solve_sweep` solves the S-matrix at 10, 64 and 128 MHz on one mesh.
The same outputs from a study file are on
[Four-port EM outputs from YAML](four_port_outputs_from_yaml.md) and
[A frequency sweep from YAML](sweep_from_yaml.md).

## Before you start

You need the complex DolfinX build (`source /usr/local/bin/dolfinx-complex-mode`)
and **2 MPI ranks**. The SAR block takes ≈ 2 min; the other two ≈ 1 min each.
Each block is a complete script.

## Run it

The SAR block solves the S-matrix first, so it can read the solved sheets'
azimuths from the result and put one centre per port at radius 15 mm on its
port's azimuth. `Study.sar` then solves the four port drives once more,
keeping their fields, and returns one SAR record per drive per centre.

```python sar-python
import numpy as np

from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial

study = Study(
    FourPortBirdcageGeometry.f_small(),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    LumpedPortDrive(),
    frequency_hz=64.0e6,
    degree=1,
)
result = study.solve()
centres = [
    (0.015 * np.cos(np.radians(s["azimuth_deg"])),
     0.015 * np.sin(np.radians(s["azimuth_deg"])), 0.0)
    for s in result.sheets
]
sar = study.sar(mass_g=10.0, centres=centres)
print(sorted(sar), [len(per_centre) for per_centre in sar.values()])
```

`Study.b1_plus` projects each drive's magnetic flux density onto a continuous
field once; `Study.b1_plus_at` evaluates `|B₁⁺| = |B_x + jB_y|/2` from it at
the points you give, on every rank.

```python b1-python
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial

study = Study(
    FourPortBirdcageGeometry.f_small(),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    LumpedPortDrive(),
    frequency_hz=64.0e6,
    degree=1,
)
maps = study.b1_plus()
b1 = study.b1_plus_at([(0.0, 0.0, 0.0)])
print(sorted(maps), {port: bool(valid.all()) for port, (_, valid) in b1.items()})
```

`Study.solve_sweep` meshes once and solves the 4×4 S-matrix at each frequency;
the phantom keeps its σ and ε_r at every frequency.

```python sweep-python
from fem_em_solver import FourPortBirdcageGeometry, LumpedPortDrive, Study
from fem_em_solver.materials import GelledSalinePhantomMaterial

study = Study(
    FourPortBirdcageGeometry.f_small(),
    GelledSalinePhantomMaterial(0.5, 78.0, 64.0e6),
    LumpedPortDrive(),
    frequency_hz=64.0e6,
    degree=1,
)
sweep = study.solve_sweep([10.0e6, 64.0e6, 128.0e6])
print([float(f) for f in sweep.frequencies_hz], sweep.s_matrices.shape)
```

## What you should see

The C4 symmetry of the birdcage, read from each Study. The 10 g SAR at port
k's own centre under port k's drive is compared cyclically between
neighbours; `|B₁⁺|` under ports 2, 4 and 3 on a point set rotated by +Δ, −Δ
and 2Δ (Δ the port spacing) is compared with port 1 on the unrotated set; and
each S-matrix, the 64 MHz one and the three swept ones, is compared with the
project's degree-1 record at that frequency. The test, not the blocks,
computes these comparisons; the pairing reads the solved sheets' azimuths,
not the port names. The bands they are held to:

```text four-port-python-bands
S vs the P1 records, per frequency: 1e-09
10 g C4 cyclic pairs: 5e-02
|B1+| C4 covariance: 5e-02
```

## What this does not show

These are per-1 V port-drive records, not absolute or compliance SAR or
`|B₁⁺|`; there is no spatial-peak search and no homogeneity figure. The sweep
holds the phantom's 64 MHz σ and ε_r at every frequency. Degree 1 only; the
coil is the solid σ = 800 S/m conductor, not copper (a Study with a
`conductor` is refused by `solve_sweep`); no Touchstone file is written here
(the YAML sweep page writes one); no F-human; no external-solver comparison.

## Worked examples

- [mri:14, the same SAR and `|B₁⁺|` outputs as a worked example](../../examples/mri/14_study_sar_and_b1.md)
- [mri:09, a study from YAML](../../examples/mri/09_study_from_yaml.md)

## How this page is tested

[tests/package/test_wf17_four_port_python_tutorial.py](../../tests/package/test_wf17_four_port_python_tutorial.py)
extracts the three blocks, execs each at 2 ranks in the complex build, and
recomputes the comparisons above from each block's Study at their bands, the
way the gates
[test_wf16_study_sar.py](../../tests/package/test_wf16_study_sar.py),
[test_wf16_study_b1.py](../../tests/package/test_wf16_study_b1.py) and
[test_wf16_study_sweep.py](../../tests/package/test_wf16_study_sweep.py)
do, and the YAML page's test
[test_wf17_four_port_outputs_tutorial.py](../../tests/package/test_wf17_four_port_outputs_tutorial.py)
does for the study files.
