# Your first study: mean SAR in a gelled-saline sphere

## What you will compute

The mean specific absorption rate (SAR) in a 1 cm-radius gelled-saline sphere
(σ = 0.57 S/m, εᵣ = 78) sitting in an imposed, uniform 1 V/m field at 64 MHz
(the 1.5 T Larmor frequency). The package composes the geometry, the material,
the imposed field, a degree-2 time-harmonic solve and the SAR integral behind
one object, [`Study`](../../src/fem_em_solver/study.py). The answer has a
closed form, so you can check it.

## Before you start

- The `fem-em-solver` Docker service must be up
  (`docker compose -f docker/docker-compose.yml up -d`), and you work inside it
  with the repository at `/workspace`.
- The solve is complex-valued: switch the container shell to the complex
  DolfinX build first with `source /usr/local/bin/dolfinx-complex-mode`.
  Real mode raises.
- Cost: ≈ 5 min wall time and ≈ 20 GiB of RAM summed over 2 MPI ranks.

## The first study

```python first-study
from fem_em_solver import (
    GelledSalinePhantomMaterial,
    QuasiStaticSphereField,
    SphereInBoxGeometry,
    Study,
)

R = 0.01  # sphere radius, m

geometry = SphereInBoxGeometry(
    sphere_radius=R,
    box_half_width=0.02,
    resolution_sphere=R / 10,
    resolution_far=R / 5,
    sphere_tag=1,
)
material = GelledSalinePhantomMaterial(0.57, 78.0, 64e6)  # sigma S/m, eps_r, Hz
field = QuasiStaticSphereField(e0=1.0)  # 1 V/m, along z

study = Study(geometry, material, field, frequency_hz=64e6)
sar = study.mean_sar()
print(f"mean SAR = {sar['mean_sar_w_per_kg']:.10e} W/kg")
```

## Run it

Save the block as a Python script and run it inside the container, in complex
mode, on 2 ranks (`<script>` is the file you saved):

```bash
source /usr/local/bin/dolfinx-complex-mode
PYTHONPATH=/workspace/src mpiexec -n 2 python3 <script>
```

## What you should see

```text first-study-output
mean SAR = 8.0116577245e-08 W/kg
```

The last digit or two can differ between runs (≈ 1e-11 relative: the parallel
direct solver is not bitwise reproducible across processes); the test below
holds the value to 1e-9 relative. The closed form for a lossy sphere in a
uniform quasi-static field is 8.008354e-08 W/kg, 0.04 % away.

## What this does not show

No coil, no ports, no h-convergence study: the field is imposed on the box
wall, not produced by an antenna, and the 0.04 % is one mesh at degree 2.

## How this page is tested

The block above is extracted from this page verbatim and executed through the
logging harness by
[`test_wf17_first_study_tutorial.py`](../../tests/package/test_wf17_first_study_tutorial.py);
it must reproduce the record of
[`test_wf15_study_lossy_sphere.py`](../../tests/package/test_wf15_study_lossy_sphere.py)
to 1e-9 relative and stay inside the closed-form band.
