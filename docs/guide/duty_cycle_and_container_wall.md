# Duty-cycled heating and a container wall on the gel ball

## What you will compute

The transient temperature rise of a gel ball cooled through its surface by a
heat-transfer coefficient (a Robin condition) and heated by a uniform
volumetric source that a duty cycle switches on and off. Time stepping is
Crank–Nicolson with a Rannacher start-up (two backward-Euler half steps).
Then the same ball, with the source left on, inside a container wall whose
conductivity and heat capacity differ from the gel's. The wall carries no
source, and the cooling acts on its outer surface.

## Before you start

You need the **real** DolfinX build: do not source `dolfinx-complex-mode`,
because a thermal-only study has no complex solve. Run with **2 MPI ranks**.
Each block below takes ≈ 1 min.

## A duty-cycled source

The block builds a ball of radius 0.1 m and gives it a uniform source of
2000 W/m³, gel conductivity 0.6 W/(m·K), a surface coefficient of 10 W/(m²·K)
to 20 °C air, and the gel's volumetric heat capacity. The source is on for
the first half of every 1200 s period, and the run takes 94 steps of 150 s.

```python duty-python
from fem_em_solver import (
    DutyCycle,
    RobinCooling,
    Study,
    ThermalBallGeometry,
    TransientSpec,
    UniformHeatSource,
)

study = Study(ThermalBallGeometry(0.1, 0.005), None, UniformHeatSource(2000.0))
result = study.temperature(
    RobinCooling(0.6, 10.0, 20.0),
    transient=TransientSpec(150.0, 94, 4.18e6),
    duty=DutyCycle(1200.0, 0.5),
)
```

## A container wall

The block puts the same gel ball inside a wall that reaches out to a radius
of 0.11 m, with conductivity 0.19 W/(m·K) and volumetric heat capacity
1.7e6 J/(m³·K). The run takes 40 Crank–Nicolson steps. Their total is one
fifth of the gel's diffusion time a²ρc/k.

```python wall-python
from fem_em_solver import (
    ContainerWall,
    RobinCooling,
    Study,
    ThermalBallGeometry,
    TransientSpec,
    UniformHeatSource,
)

wall = ContainerWall(0.11 - 0.1, 0.19, 1.7e6)
study = Study(ThermalBallGeometry(0.1, 0.005, wall=wall), None, UniformHeatSource(2000.0))
result = study.temperature(
    RobinCooling(0.6, 10.0, 20.0),
    transient=TransientSpec(348.3333333333334, 40, 4.18e6),
)
```

## What you should see

Each block binds `study` and `result`:

- `result.temperature` is the temperature at the final time, a P1 field on
  `result.mesh`.
- `result.max_balance_residual` is the largest discrete energy-balance
  residual over the steps.
- `result.source_multipliers` lists the schedule's on/off multiplier for each
  step and start-up half step. It is empty when no duty cycle is given.

The page states no temperature, error or residual. The test checks the
values against the recorded closed-form comparisons.

## What this does not show

- A coil or another electromagnetic source; see [temperature.md](temperature.md)
  for the imposed-field route.
- Perfusion, or any tissue.
- The same studies from a file:
  [thermal_ball_duty_a.yaml](../../tests/package/wf1/thermal_ball_duty_a.yaml)
  and [thermal_ball_wall.yaml](../../tests/package/wf1/thermal_ball_wall.yaml),
  gated by
  [test_wf1_yaml_thermal_ball.py](../../tests/package/test_wf1_yaml_thermal_ball.py).
- A CEM43 thermal dose.
- Any external-solver comparison.

## Worked examples

- [mri:25, a transient study with a duty cycle and a wall](../../examples/mri/25_study_transient_duty_and_wall.md)
- [materials:11, a duty-cycled gel ball](../../examples/materials/11_duty_cycled_gel_ball.md)
- [materials:12, a gel ball in a container wall](../../examples/materials/12_gel_ball_in_container_wall.md)

## How this page is tested

[tests/package/test_wf17_duty_wall_tutorial.py](../../tests/package/test_wf17_duty_wall_tutorial.py)
extracts both blocks and runs them at 2 ranks in the real build. It checks
that each block's inputs are the validation tests' inputs, that it builds the
same mesh, and that its temperature reproduces the recorded comparison with
the closed form, the way
[test_wf15_study_transient.py](../../tests/package/test_wf15_study_transient.py)
does.
