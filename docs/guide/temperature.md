# Steady temperature of a SAR-heated gel sphere

## What you will compute

The steady temperature of the 1 cm-radius gelled-saline sphere of
[the first study](first_study.md), heated by its own absorbed power
q = G · ½σ|E|² and cooled through its surface to 20 °C air (conductivity
k = 0.6 W/(m·K), convection h = 100 W/(m²·K)). G scales the heating as if the
imposed field were 5 kV/m instead of 1 V/m. One object,
[`Study`](../../src/fem_em_solver/study.py), solves the field once and hands it
to the conduction solve; the result has a closed form.

## Before you start

- The `fem-em-solver` Docker service must be up
  (`docker compose -f docker/docker-compose.yml up -d`), and you work inside it
  with the repository at `/workspace`.
- The field solve is complex-valued: switch the container shell to the complex
  DolfinX build first with `source /usr/local/bin/dolfinx-complex-mode`.
- Cost: ≈ 30 s wall time on 2 MPI ranks.
- **Set `degree=1`.** `Study` defaults to degree 2; the recorded temperature
  below is the degree-1 one.

## The study

```python temperature
from fem_em_solver import (
    GelledSalinePhantomMaterial,
    QuasiStaticSphereField,
    RobinCooling,
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

study = Study(geometry, material, field, frequency_hz=64e6, degree=1)
cooling = RobinCooling(0.6, 100.0, 20.0, source_scale=2.5e7)  # k, h, T_inf, G
thermal = study.temperature(cooling)  # solves the field once, then conduction
cells = thermal.submesh.topology.index_map(3).size_global
print(f"sphere cells = {cells}")
```

`thermal.temperature` is the temperature field (°C) on the sphere, and
`thermal.heat_source` the heating q (W/m³).

## Run it

Save the block as a Python script and run it inside the container, in complex
mode, on 2 ranks (`<script>` is the file you saved):

```bash
source /usr/local/bin/dolfinx-complex-mode
PYTHONPATH=/workspace/src mpiexec -n 2 python3 <script>
```

## What you should see

```text temperature-output
sphere cells = 19496
```

The test then compares the temperature rise with the closed form for a
uniformly heated, Robin-cooled sphere: the interior relative L2 error is
3.2464e-02 (band 0.11), and the heating integrates to the absorbed power and
balances the convected power to ≤ 1e-10 and ≤ 1e-9 relative.

## What this does not show

No transient heating, no perfusion, no tissue and no coil: the field is imposed
on the box wall, and this is one mesh at degree 1. Heating over time from the
same `Study` is on [the transient page](sphere_heating_transient.md).

## How this page is tested

The block above is extracted from this page verbatim and executed through the
logging harness by
[`test_wf17_more_tutorials.py`](../../tests/package/test_wf17_more_tutorials.py);
its error against the closed form must equal the record of
[`test_wf15_study_temperature.py`](../../tests/package/test_wf15_study_temperature.py)
to four digits and stay inside the closed-form band.
