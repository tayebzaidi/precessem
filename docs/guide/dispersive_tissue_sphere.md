# A dispersive medium: mean SAR in a Cole–Cole sphere

## What you will compute

The mean specific absorption rate (SAR) in the same 1 cm-radius sphere as
[Your first study](first_study.md), in the same imposed, uniform 1 V/m field
at 64 MHz, but filled with a **dispersive** medium: a two-pole Cole–Cole
permittivity with an ionic conductivity,

ε̂(ω) = ε∞ + Σₖ Δεₖ / (1 + (jωτₖ)^(1 − αₖ)) + σᵢ / (jωε₀).

The medium is a [`ColeColeMaterial`](../../src/fem_em_solver/materials/dispersion.py)
handed to [`Study`](../../src/fem_em_solver/study.py) in place of a fixed
(σ, εᵣ) pair. The numbers below are an illustrative two-pole medium, not a
tissue.

## Before you start

- The `fem-em-solver` Docker service must be up
  (`docker compose -f docker/docker-compose.yml up -d`), and you work inside it
  with the repository at `/workspace`.
- The solve is complex-valued: switch the container shell to the complex
  DolfinX build first with `source /usr/local/bin/dolfinx-complex-mode`.
  Real mode raises.
- Cost: ≈ 5 min wall time on 2 MPI ranks, the same mesh and degree-2 solve as
  the first study.

## The dispersive study

```python dispersive-sphere
import math

from fem_em_solver import QuasiStaticSphereField, SphereInBoxGeometry, Study
from fem_em_solver.materials import ColeColeMaterial

R = 0.01  # sphere radius, m

geometry = SphereInBoxGeometry(
    sphere_radius=R,
    box_half_width=0.02,
    resolution_sphere=R / 10,
    resolution_far=R / 5,
    sphere_tag=1,
)
medium = ColeColeMaterial(
    epsilon_inf=4.0,
    poles=(  # (delta_eps, tau s, alpha)
        (50.0, 1.0 / (2.0 * math.pi * 100e6), 0.30),
        (4000.0, 1.0 / (2.0 * math.pi * 1e6), 0.35),
    ),
    sigma_ionic=0.2,  # S/m
)
field = QuasiStaticSphereField(e0=1.0)  # 1 V/m, along z

study = Study(geometry, medium, field, frequency_hz=64e6)
print(
    f"resolved at 64 MHz: sigma = {study.material.sigma:.6g} S/m, "
    f"epsilon_r = {study.material.epsilon_r:.6g}"
)
sar = study.mean_sar()
print(f"mean SAR = {sar['mean_sar_w_per_kg']:.6e} W/kg")
```

## Run it

Save the block as a Python script and run it inside the container, in complex
mode, on 2 ranks (`<script>` is the file you saved):

```bash
source /usr/local/bin/dolfinx-complex-mode
PYTHONPATH=/workspace/src mpiexec -n 2 python3 <script>
```

## What you should see

```text dispersive-sphere-output
resolved at 64 MHz: sigma = 1.0098 S/m, epsilon_r = 181.279
mean SAR = 3.986757e-08 W/kg
```

The closed form for a lossy sphere in a uniform quasi-static field, evaluated
at the resolved pair, is 3.985081e-08 W/kg, 0.04 % away.

## Why the medium is resolved at the solve frequency

A time-harmonic solve runs at one frequency, so the `Study` evaluates ε̂ once,
at `frequency_hz`, and solves with the real pair σ_eff = −ωε₀ Im ε̂ and
εᵣ = Re ε̂ — here 1.0098 S/m and 181.279, far from σᵢ = 0.2 S/m and
ε∞ = 4 because both poles still contribute at 64 MHz. Two shortcuts
give a different answer: dropping the Cole–Cole broadening (αₖ = 0, a Debye
medium) or dropping the poles (σᵢ and ε∞ alone) puts the closed-form SAR
more than 60 % away from the Cole–Cole one, and the test below holds the page's SAR
away from the Debye one. A sweep over frequency resolves the medium afresh at
each frequency; this page solves at one.

## What this does not show

No tissue table: the Gabriel parameters are a separate, operator-held step.
No coil or port study: a `ColeColeMaterial` on a port geometry is refused
with a `TypeError` before any mesh. No YAML: this medium is Python-only here.

## How this page is tested

The block above is extracted from this page verbatim and executed through the
logging harness by
[`test_wf17_dispersive_sphere_tutorial.py`](../../tests/package/test_wf17_dispersive_sphere_tutorial.py);
it must reproduce the record of
[`test_wf15_study_cole_cole.py`](../../tests/package/test_wf15_study_cole_cole.py)
at every digit that record carries, stay inside the closed-form band, and
miss the Debye-reduced closed form.
