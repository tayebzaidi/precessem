# The SAR-heated sphere over time, from one `Study`

## What you will compute

The temperature rise of the 1 cm-radius gelled-saline sphere of
[the steady temperature page](temperature.md) **over time**, heated by its own
absorbed power q = G · ½σ|E|² from T(0) = T∞ = 0 and cooled through its surface
(k = 0.6 W/(m·K), h = 100 W/(m²·K), G = 2.5e7). One
[`Study`](../../src/fem_em_solver/study.py) solves the field once; two
`temperature(…, transient=TransientSpec(…))` calls then march the conduction
equation on the sphere without a second field solve:

- **A:** Crank–Nicolson (θ = 0.5, two backward-Euler half-steps first), 40 steps
  to t* = 0.2 R²ρc/k ≈ 139 s (Fourier number 0.2);
- **S:** one backward-Euler step of dt₁ = 1e-5 R²ρc/k ≈ 7 ms — so short that the
  sphere has not yet conducted anything away, and the stored heat must equal
  the deposited heat (calorimetry).

Both results have a closed form: the Robin-cooled sphere's eigen-series at t*,
and ρc·mean(θ) = dt₁·mean(q) for the single short step.

## Before you start

- The `fem-em-solver` Docker service must be up
  (`docker compose -f docker/docker-compose.yml up -d`), and you work inside it
  with the repository at `/workspace`.
- The field solve is complex-valued: switch the container shell to the complex
  DolfinX build first with `source /usr/local/bin/dolfinx-complex-mode`.
- Cost: ≈ 30 s wall time on 2 MPI ranks; "At the `Study`'s default degree"
  costs ≈ 5 min and ≈ 21 GiB summed over 2 ranks.
- **Set `degree=1`** (except in "At the `Study`'s default degree"). `Study`
  defaults to degree 2; the records below are the degree-1 ones.
- T∞ = 0 here, so the temperature *is* the temperature rise θ.

## The study

```python sphere-transient
from fem_em_solver import (
    GelledSalinePhantomMaterial,
    QuasiStaticSphereField,
    RobinCooling,
    SphereInBoxGeometry,
    Study,
    TransientSpec,
)

R = 0.01  # sphere radius, m
k = 0.6  # W/(m K)
rho_c = 1000.0 * 4180.0  # density x specific heat, J/(m^3 K)
t_r = R**2 * rho_c / k  # conduction time scale, s
t_star = 0.2 * t_r  # Fourier number 0.2
n_steps = 40
dt_1 = 1.0e-5 * t_r  # the calorimetric step

geometry = SphereInBoxGeometry(
    sphere_radius=R,
    box_half_width=0.02,
    resolution_sphere=R / 10,
    resolution_far=R / 5,
    sphere_tag=1,
)
material = GelledSalinePhantomMaterial(0.57, 78.0, 64e6)  # sigma S/m, eps_r, Hz
field = QuasiStaticSphereField(e0=1.0)  # 1 V/m, along z

study = Study(geometry, material, field, frequency_hz=64e6, density_kg_m3=1000.0, degree=1)
result = study.solve()  # the one field solve
cooling = RobinCooling(k, 100.0, 0.0, source_scale=2.5e7)  # k, h, T_inf, G

heating = study.temperature(cooling, transient=TransientSpec(t_star / n_steps, n_steps, rho_c))
first_step = study.temperature(
    cooling, transient=TransientSpec(dt_1, 1, rho_c, theta=1.0, startup_half_steps=0)
)
cells = heating.mesh.topology.index_map(3).size_global
print(f"sphere cells = {cells}")
```

`heating.temperature` is θ(t*) (K) on the sphere, `heating.q` the heating
(W/m³), and `heating.max_balance_residual` the worst per-step energy-balance
residual; `first_step` carries the same fields for the single short step.

## Run it

Save the block as a Python script and run it inside the container, in complex
mode, on 2 ranks (`<script>` is the file you saved):

```bash
source /usr/local/bin/dolfinx-complex-mode
PYTHONPATH=/workspace/src mpiexec -n 2 python3 <script>
```

## What you should see

```text sphere-transient-output
sphere cells = 19496
(A_t) 3.4708e-02
(S2) 4.9994e-05
```

The block prints the first line. The other two are what the test computes from
the block's results and must reproduce to four digits: **(A_t)** is the
interior relative L2 error of θ(t*) against the closed-form series (band 0.11),
and **(S2)** is |ρc·mean(θ) / (dt₁·mean(q)) − 1| after the single short step
(band 1e-04). The test also asserts the per-step energy balance of both runs to
≤ 1e-9 relative, and that the field was solved once.

## The same heating, pulsed

A `DutyCycle` switches the source on for the first half of each period
P = t\*/5 (8 steps) and off for the second. The block repeats the study above
and runs the continuous heating beside the pulsed one, both from the one field
solve:

```python sphere-transient-duty
from fem_em_solver import (
    DutyCycle,
    GelledSalinePhantomMaterial,
    QuasiStaticSphereField,
    RobinCooling,
    SphereInBoxGeometry,
    Study,
    TransientSpec,
)

R = 0.01  # sphere radius, m
k = 0.6  # W/(m K)
rho_c = 1000.0 * 4180.0  # density x specific heat, J/(m^3 K)
t_r = R**2 * rho_c / k  # conduction time scale, s
t_star = 0.2 * t_r  # Fourier number 0.2
n_steps = 40

geometry = SphereInBoxGeometry(
    sphere_radius=R,
    box_half_width=0.02,
    resolution_sphere=R / 10,
    resolution_far=R / 5,
    sphere_tag=1,
)
material = GelledSalinePhantomMaterial(0.57, 78.0, 64e6)  # sigma S/m, eps_r, Hz
field = QuasiStaticSphereField(e0=1.0)  # 1 V/m, along z

study = Study(geometry, material, field, frequency_hz=64e6, density_kg_m3=1000.0, degree=1)
result = study.solve()  # the one field solve
cooling = RobinCooling(k, 100.0, 0.0, source_scale=2.5e7)  # k, h, T_inf, G

dt = t_star / n_steps
heating = study.temperature(cooling, transient=TransientSpec(dt, n_steps, rho_c))
pulsed = study.temperature(cooling, transient=TransientSpec(dt, n_steps, rho_c), duty=DutyCycle(8 * dt, 0.5))
cells = heating.mesh.topology.index_map(3).size_global
print(f"sphere cells = {cells}")
```

`pulsed.temperature` is θ(t\*) under the pulsed source, and
`pulsed.source_multipliers` the on/off factor the march applied at each step.
The switches fall on step boundaries, which is why the period is written as
`8 * dt`: a period that is not a whole number of steps would switch the source
part-way through a step.

```text sphere-transient-duty-output
sphere cells = 19496
(A_t) 3.4708e-02
```

The block prints the first line; the second is the continuous run's (A_t),
the record of "What you should see". The test checks the pulsed run against the
superposition of ten switched copies of the step-on run (Duhamel, `MAT-17`,
band 1e-12), and checks that it is *not* the continuous closed form.

## At the `Study`'s default degree

Leave `degree=` out and the `Study` solves the field with second-order
elements, its default; the heat equation stays first-order. The block is the
pulsed block above with that one keyword removed, and it prints the degree the
study used:

```python sphere-transient-degree2
from fem_em_solver import (
    DutyCycle,
    GelledSalinePhantomMaterial,
    QuasiStaticSphereField,
    RobinCooling,
    SphereInBoxGeometry,
    Study,
    TransientSpec,
)

R = 0.01  # sphere radius, m
k = 0.6  # W/(m K)
rho_c = 1000.0 * 4180.0  # density x specific heat, J/(m^3 K)
t_r = R**2 * rho_c / k  # conduction time scale, s
t_star = 0.2 * t_r  # Fourier number 0.2
n_steps = 40

geometry = SphereInBoxGeometry(
    sphere_radius=R,
    box_half_width=0.02,
    resolution_sphere=R / 10,
    resolution_far=R / 5,
    sphere_tag=1,
)
material = GelledSalinePhantomMaterial(0.57, 78.0, 64e6)  # sigma S/m, eps_r, Hz
field = QuasiStaticSphereField(e0=1.0)  # 1 V/m, along z

study = Study(geometry, material, field, frequency_hz=64e6, density_kg_m3=1000.0)
result = study.solve()  # the one field solve
cooling = RobinCooling(k, 100.0, 0.0, source_scale=2.5e7)  # k, h, T_inf, G

dt = t_star / n_steps
heating = study.temperature(cooling, transient=TransientSpec(dt, n_steps, rho_c))
pulsed = study.temperature(cooling, transient=TransientSpec(dt, n_steps, rho_c), duty=DutyCycle(8 * dt, 0.5))
cells = heating.mesh.topology.index_map(3).size_global
print(f"sphere cells = {cells}")
print(f"study degree = {study.degree}")
```

```text sphere-transient-degree2-output
sphere cells = 19496
study degree = 2
```

The block prints both lines. The test checks the continuous run against the
same closed-form series and band (0.11) as "What you should see", the pulsed
run against the Duhamel superposition (band 1e-12), and that the pulsed run is
*not* the continuous closed form. This section costs ≈ 5 min wall time and
≈ 21 GiB of memory summed over 2 ranks.

## The same two runs from YAML files

The two runs are also shipped as study files:
[`sphere_mat12_transient_a.yaml`](../../tests/package/wf1/sphere_mat12_transient_a.yaml)
(run A) and
[`sphere_mat12_transient_s.yaml`](../../tests/package/wf1/sphere_mat12_transient_s.yaml)
(run S). Each carries the same sphere, field and cooling as the block above,
with `study.degree: 1` and T∞ = 0, and a `thermal.transient` mapping that
names the end time `t_end_s` and the step count `n_steps` (the step is
`t_end_s / n_steps`) beside ρc, θ and the startup half-steps. Each file is one
field solve and one transient, so the two runs below solve the field twice.

```python sphere-transient-yaml
from fem_em_solver import run_study_file

FILES = ["sphere_mat12_transient_a.yaml", "sphere_mat12_transient_s.yaml"]
runs = {name: run_study_file(name, output_dir=f"femem_out/{name[:-5]}") for name in FILES}
for name, run in runs.items():
    print(name, "cells", run.summary["cells"])
```

`summary["cells"]` counts the whole meshed box; `run.thermal` is the transient
result on the sphere, with the same fields as `heating` and `first_step`
above. The test computes (A_t) from the A run and (S2) from the S run and
requires the records of "What you should see" to four digits.

A duty cycle has no key on this sphere: this file adds `thermal.duty` to the
A file, and `femem check` exits with status 2 and names `thermal.duty`:

```bash sphere-transient-yaml-check
femem check malformed_sphere_transient_duty.yaml
```

## What this does not show

Degree 2 from Python only, with no degree-2 record on the page (its test asserts the bands), one duty schedule only (P = t\*/5, D = 0.5, from Python; no
periodic steady state), no perfusion, no probes, no tissue and
no coil: the field is imposed on the box wall, and this is one mesh at
degree 1. `Study.probe` is refused on this sphere with `transient`.

## How this page is tested

The block above is extracted from this page verbatim and executed through the
logging harness by
[`test_wf17_sphere_transient_tutorial.py`](../../tests/package/test_wf17_sphere_transient_tutorial.py);
its errors against the closed forms must equal the records of
[`test_wf15_study_sphere_transient.py`](../../tests/package/test_wf15_study_sphere_transient.py)
to four digits and stay inside the closed-form bands.
[`test_wf17_sphere_transient_yaml_tutorial.py`](../../tests/package/test_wf17_sphere_transient_yaml_tutorial.py)
runs the YAML block the same way, on the two shipped files, against the same
records and bands.
[`test_wf17_sphere_transient_duty_tutorial.py`](../../tests/package/test_wf17_sphere_transient_duty_tutorial.py)
runs the pulsed block and checks it against the Duhamel superposition.
[`test_wf17_sphere_transient_degree2_tutorial.py`](../../tests/package/test_wf17_sphere_transient_degree2_tutorial.py)
runs the default-degree block and checks the degree, the Duhamel superposition
and the closed-form band.
